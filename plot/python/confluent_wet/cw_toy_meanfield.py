#!/usr/bin/env python3
"""Non-spatial toy model of the tau_m-controlled bistability: N cells, synthetic pressure.

Each cell carries an Ornstein-Uhlenbeck pressure p, a memory m and a phenotype chi
(hill convention, chi = 1 active); the cells interact ONLY through the population activity
x = <chi>, which sets the pressure statistics. Time is in the reference tau_c (674.33 steps).

    dp       = -(p - mu(a))/tau_p(a) dt + sigma(a) sqrt(2/tau_p(a)) dW     a = z0 + (1-z0) x
    tau_m dm = (s(p) - m) dt
    m_eff    = mbar + (m - mbar)/sqrt(c_D)        mbar = population mean of m
    tau_chi dchi = (h(m_eff) + beta g(chi)(1 - h(m_eff)) - chi) dt

INPUTS, all from the OPEN-LOOP ladder (six activity rungs, L = 500, no closed-loop data):
mu(a), sigma(a) (Gaussian pressure moments), tau_p(a) (Lagrangian pressure correlation time
along tracers), and c_D. c_D is the single constant that stands for the damping of the
memory width by memory diffusion, which a non-spatial model does not have: it is fitted
once, c_D = v_toy/std_m^2 averaged in log over the rungs, where v_toy is this model's own
memory variance at the rung's activity and tau_m (= 10 tau_c) and std_m is the measured
spatial std. The m_eff map applies it exactly -- same mean, same correlation time, variance
divided by c_D -- which is the constant-c_D closure assumption itself, not a fit to its
consequences.

Integration: the OU step is exact, m and chi use the exponential (exact for an input held
over dt) update, dt = 0.02 tau_c << min tau_p.

Modes
  calibrate                  fit c_D on the rungs
  map                        N = infinity: F(x; tau_m) = <h(m_eff)> at fixed x, fixed points,
                             saddle-nodes; plus the Gaussian/Beta closure for comparison
  cells                      finite-N collective runs mirroring the spatial scans
  figure                     summary figure against the spatial scans
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.special import betainc, ndtr
from scipy.stats import beta as beta_dist

TC = 674.3290333006435
CODE = Path(__file__).resolve().parents[4]
CALIB = CODE/'results/20260904_cw_s3/calib_s3.json'
LADDER = CODE/'results/20260908_gaussian_beta_phase/calibration_inputs.json'
Z0, PMEM, MC = 0.3, 0.016838, 0.2287
TCHI0 = 202.3/TC                      # the scans' tau_chi = 0.3 tau_c
LAWS = {'hill': {'h-n': 16., 'wp': .002}, 'step': {'h-n': 0., 'wp': 0.}}


class Inputs:
    """Open-loop pressure statistics against activity a in [0.3, 1]."""
    def __init__(self, calib=CALIB, ladder=LADDER):
        c = json.loads(Path(calib).read_text())
        lad = json.loads(Path(ladder).read_text())
        p = sorted(lad['pressure'], key=lambda d: d['a'])
        self.a = np.array([d['a'] for d in p])
        rung = lambda a: c['rungs'][min(c['rungs'], key=lambda k: abs(float(k)-a))]
        self.mu_k = np.array([d['mean'] for d in p])
        self.sig_k = np.array([d['std'] for d in p])
        self.tau_k = np.array([rung(a)['tau_c']/TC for a in self.a])
        assert np.allclose(self.sig_k, [rung(a)['sigma_P'] for a in self.a])
        ladder_rows = sorted(lad['theory']['ladder'], key=lambda d: d['a'])
        self.std_m_k = np.array([d['std_m'] for d in ladder_rows])
        self.f_k = np.array([d['f'] for d in ladder_rows])
        self.g_open = float(ladder_rows[0]['g'])       # open-loop tau_m in reference tau_c
        assert all(abs(d['g']-self.g_open) < 1e-9 for d in ladder_rows)
        self.mu = PchipInterpolator(self.a, self.mu_k)
        self.sig = PchipInterpolator(self.a, self.sig_k)
        self.tau = PchipInterpolator(self.a, self.tau_k)

    def at(self, a):
        a = np.clip(a, self.a[0], self.a[-1])
        return self.mu(a), self.sig(a), self.tau(a)


class Law:
    def __init__(self, hn, wp, mc=MC, pmem=PMEM, beta=0., gn=2., gchic=.5):
        self.hn, self.wp, self.mc, self.pmem = hn, wp, mc, pmem
        self.beta, self.gn, self.gchic = beta, gn, gchic
        self.gnorm = 1 + gchic**gn

    def s(self, p):
        if self.wp == 0:
            return (p > self.pmem).astype(float)
        return 1/(1 + np.exp(-(p - self.pmem)/self.wp))

    def h(self, m):
        if self.hn == 0:
            return (m >= self.mc).astype(float)
        with np.errstate(divide='ignore', over='ignore'):
            return np.where(m > 0, 1/(1 + (self.mc/np.maximum(m, 1e-300))**self.hn), 0.)

    def g(self, x):
        with np.errstate(divide='ignore', over='ignore'):
            return np.where(x > 0, self.gnorm/(1 + (self.gchic/np.maximum(x, 1e-300))**self.gn), 0.)

    def target(self, m, chi):
        H = self.h(m)
        return H if self.beta == 0 else H + self.beta*self.g(chi)*(1 - H)

    def s_moments(self, mu, sig):
        """E[s], Var[s] for Gaussian p. Analytic for the step; for the logistic a fine
        quadrature -- Gauss-Hermite nodes cannot resolve a switch of width wp << sigma."""
        mu, sig = np.asarray(mu, float), np.asarray(sig, float)
        if self.wp == 0:
            f = ndtr((mu - self.pmem)/sig)
            return f, f*(1 - f)
        z = np.linspace(-9, 9, 6001)
        w = np.exp(-z*z/2); w /= w.sum()
        sv = self.s(mu[..., None] + sig[..., None]*z)
        e1 = (sv*w).sum(-1)
        return e1, np.maximum((sv**2*w).sum(-1) - e1**2, 0.)


def ou_step(p, mu, sig, tau, dt, rng):
    e = np.exp(-dt/tau)
    return mu + (p - mu)*e + sig*np.sqrt(1 - e*e)*rng.standard_normal(p.shape)


# ---- c_D --------------------------------------------------------------------------------

def calibrate(inp, cells=20000, dt=.02, seed=1):
    """Toy memory variance at each open-loop rung (step s, as the ladder runs) vs measured."""
    rng = np.random.default_rng(seed)
    law = Law(0., 0.)
    a = inp.a[:, None]
    mu, sig, tau = (np.broadcast_to(v, (len(inp.a), cells)) for v in inp.at(a))
    p = mu + sig*rng.standard_normal(mu.shape)
    m = np.broadcast_to(law.s_moments(inp.mu_k, inp.sig_k)[0][:, None], mu.shape).copy()
    tm, burn, window = inp.g_open, 5*inp.g_open + 20, 20*inp.g_open
    em, acc, n = np.exp(-dt/tm), np.zeros(len(inp.a)), 0
    for k in range(int((burn + window)/dt)):
        p = ou_step(p, mu, sig, tau, dt, rng)
        s = law.s(p)
        m = s + (m - s)*em
        if k*dt >= burn and k % 25 == 0:
            acc += m.var(axis=1); n += 1
    v_toy = acc/n
    c = v_toy/inp.std_m_k**2
    f_gauss, var_s = law.s_moments(inp.mu_k, inp.sig_k)
    return {'a': inp.a.tolist(), 'tau_m_open_over_tc': tm,
            'tau_p_over_tc': inp.tau_k.tolist(), 'sigma_P': inp.sig_k.tolist(),
            'f_measured': inp.f_k.tolist(), 'f_gaussian': f_gauss.tolist(),
            'std_m_measured': inp.std_m_k.tolist(), 'std_m_toy_cD1': np.sqrt(v_toy).tolist(),
            'onepole_std': np.sqrt(var_s/(1 + tm/inp.tau_k)).tolist(),
            'c_D_per_rung': c.tolist(), 'c_D': float(np.exp(np.mean(np.log(c)))),
            'c_D_log_spread': float(np.std(np.log(c)))}


# ---- N = infinity ---------------------------------------------------------------------------

def fmap(inp, law, tms, cD, xs, cells=3000, dt=.02, seed=2):
    """F(x; tau_m) = <h(m_eff)> at fixed x (no feedback), cells independent within a row."""
    rng = np.random.default_rng(seed)
    out = np.zeros((len(tms), len(xs)))
    a = (Z0 + (1 - Z0)*np.asarray(xs))[:, None]
    mu, sig, tau = (np.broadcast_to(v, (len(xs), cells)) for v in inp.at(a))
    f0 = law.s_moments(*(v[:, 0] for v in (mu, sig)))[0][:, None]
    shrink = 1/np.sqrt(cD)
    for i, tm in enumerate(tms):
        p = mu + sig*rng.standard_normal(mu.shape)
        m = np.broadcast_to(f0, mu.shape).copy()
        burn, window = 5*tm + 5*inp.tau_k.max(), max(60., 6*tm)
        em, acc, n = np.exp(-dt/tm), np.zeros(len(xs)), 0
        for k in range(int((burn + window)/dt)):
            p = ou_step(p, mu, sig, tau, dt, rng)
            s = law.s(p)
            m = s + (m - s)*em
            if k*dt >= burn and k % 25 == 0:
                mbar = m.mean(axis=1, keepdims=True)
                acc += law.h(mbar + (m - mbar)*shrink).mean(axis=1); n += 1
        out[i] = acc/n
        print(json.dumps({'tau_m': tm, 'F_range': [float(out[i].min()), float(out[i].max())]}),
              flush=True)
    return out


def closure_map(inp, law, tms, cD, xs):
    """Gaussian pressure + Beta memory: mean E[s], variance Var[s]/(c_D (1 + tau_m/tau_p))."""
    a = Z0 + (1 - Z0)*np.asarray(xs)
    mu, sig, tau = inp.at(a)
    f, var_s = law.s_moments(mu, sig)
    out = np.zeros((len(tms), len(xs)))
    grid = np.linspace(0, 1, 4001)[1:-1]
    for i, tm in enumerate(tms):
        v = np.minimum(var_s/(cD*(1 + tm/tau)), f*(1 - f)*(1 - 1e-9))
        k = f*(1 - f)/v - 1
        al, be = f*k, (1 - f)*k
        if law.hn == 0:
            out[i] = 1 - betainc(al, be, law.mc)
        else:
            pdf = beta_dist.pdf(grid[None, :], al[:, None], be[:, None])
            out[i] = np.trapz(pdf*law.h(grid)[None, :], grid, axis=1)/np.trapz(pdf, grid, axis=1)
    return out


def fixed_points(xs, F):
    """Roots of F(x) - x on the grid; stable where the crossing slope of F - x is negative."""
    d = F - np.asarray(xs)
    roots = []
    for j in range(len(xs) - 1):
        if d[j] == 0 or d[j]*d[j+1] < 0:
            t = d[j]/(d[j] - d[j+1]) if d[j] != d[j+1] else 0.
            roots.append({'x': float(xs[j] + t*(xs[j+1] - xs[j])), 'stable': bool(d[j] > d[j+1])})
    if d[-1] > 0:
        roots.append({'x': 1., 'stable': True})
    if d[0] < 0:
        roots.insert(0, {'x': 0., 'stable': True})
    return roots


# ---- finite N --------------------------------------------------------------------------------

def cells_run(inp, law, rows, N, cD, dt=.02, seed=3, every=.5):
    """All rows integrated together: each row is (tau_m, tau_chi, start) with its own N cells.

    Preparation max(100, 10 tau_m) with the chi reaction frozen (as the spatial scans), then
    2000 tau_c; x(t) is recorded every `every`. Starts: active chi=1, m=.3301; passive chi=0,
    m=.0559 (the scans' uniform starts)."""
    rng = np.random.default_rng(seed)
    R = len(rows)
    tm = np.array([r['tau_m'] for r in rows])[:, None]
    tx = np.array([r['tau_chi'] for r in rows])[:, None]
    prep = np.maximum(100., 10*tm[:, 0])
    chi = np.array([1. if r['start'] == 'active' else 0. for r in rows])[:, None]*np.ones((R, N))
    m = np.array([.3301 if r['start'] == 'active' else .0559 for r in rows])[:, None]*np.ones((R, N))
    x = chi.mean(axis=1)
    mu, sig, tau = inp.at(Z0 + (1 - Z0)*x)
    p = mu[:, None] + sig[:, None]*rng.standard_normal((R, N))
    em, ec, shrink = np.exp(-dt/tm), np.exp(-dt/tx), 1/np.sqrt(cD)
    total = prep.max() + 2000.
    stride = int(round(every/dt))
    tt, xt, mt = [], [], []
    for k in range(int(round(total/dt)) + 1):
        t = k*dt
        released = (t - (prep.max() - prep)) >= prep            # per-row release time
        # rows start at different times so that every row ends at its own prep + 2000
        alive = t >= prep.max() - prep
        x = chi.mean(axis=1)
        mu, sig, tau = inp.at(Z0 + (1 - Z0)*x)
        p = ou_step(p, mu[:, None], sig[:, None], tau[:, None], dt, rng)
        s = law.s(p)
        m = np.where(alive[:, None], s + (m - s)*em, m)
        mbar = m.mean(axis=1, keepdims=True)
        T = law.target(mbar + (m - mbar)*shrink, chi)
        chi = np.where((alive & released)[:, None], T + (chi - T)*ec, chi)
        if k % stride == 0:
            tt.append(t); xt.append(chi.mean(axis=1).copy()); mt.append(m.mean(axis=1).copy())
    tt, xt, mt = np.array(tt), np.array(xt).T, np.array(mt).T
    res = []
    for i, r in enumerate(rows):
        t_rel = tt - prep.max()                                 # time after release, all rows
        tail = t_rel >= 2000. - max(500., 20*r['tau_m'])
        half = t_rel >= 2000. - max(500., 20*r['tau_m'])/2
        res.append(dict(r, N=N, tail_mean=float(xt[i, tail].mean()),
                        # Var_t <chi> over the tail at the scans' 0.5 tau_c video cadence
                        tail_var=float(xt[i, tail].var()),
                        tail_m=float(mt[i, tail].mean()),
                        half_drift=float(xt[i, half].mean() - xt[i, tail & ~half].mean()),
                        series_t=t_rel[::4].tolist(), series_x=xt[i, ::4].round(5).tolist()))
    return res


# ---- figure ------------------------------------------------------------------------------------

def spatial_rows(path, legacy=False):
    """Tail means of a spatial scan as (tau_m, start, active fraction), hill convention."""
    if not Path(path).exists():
        return []
    s = json.loads(Path(path).read_text())
    out = []
    for r in s['cases']:
        if r.get('replicate', 1) != 1 or r['initialization'] not in ('0', '1'):
            continue
        if legacy:   # legacy chi = 1 passive; its all-0 start is the active start
            out.append((r['tm_over_tc'], 'active' if r['initialization'] == '0' else 'passive',
                        1 - r['tail_mean']))
        else:
            out.append((r['tm_over_tc'], 'active' if r['initialization'] == '1' else 'passive',
                        r['tail_mean']))
    return out


def figure(out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'pdf.fonttype': 42,
                         'axes.spines.top': False, 'axes.spines.right': False})
    R = CODE/'results'
    spatial = {'hill': [('spatial, τχ = 0.3 τc', spatial_rows(R/'20261013_cw_hill_tm_scan/summary/hill_tm_summary.json'), '#30363B'),
                        ('spatial, τχ = τm', spatial_rows(R/'20261014_cw_hill_tm_scan_tchi1tm/summary/hill_tm_summary.json'), '#8C6BB1')],
               'step': [('spatial legacy step', spatial_rows(R/'20260911_cw_mc02287_four_init/summary/four_init_summary.json', legacy=True), '#30363B')]}
    cal = json.loads((out/'calibration.json').read_text())
    fig, ax = plt.subplots(2, 3, figsize=(16, 9.2))
    for j, law in enumerate(('hill', 'step')):
        a = ax[0, j]
        path = out/f'map_{law}.json'
        if not path.exists():
            continue
        mp = json.loads(path.read_text())
        for tm, roots in zip(mp['tau_m'], mp['roots_closure']):
            for q in roots:
                a.plot(tm, q['x'], 'o', ms=3, color='#9AA3AB' if q['stable'] else 'none',
                       mec='#9AA3AB', zorder=1)
        for tm, roots in zip(mp['tau_m'], mp['roots']):
            for q in roots:
                a.plot(tm, q['x'], 'o', ms=5.5, color='#E07A5F' if q['stable'] else 'white',
                       mec='#E07A5F', mew=1.4, zorder=3)
        a.plot([], [], 'o', color='#E07A5F', label='toy N = ∞, stable')
        a.plot([], [], 'o', color='white', mec='#E07A5F', mew=1.4, label='toy N = ∞, unstable')
        a.plot([], [], 'o', ms=3, color='#9AA3AB', label='Beta closure (open: unstable)')
        for lab, rows, col in spatial[law]:
            for st, mk in (('active', '^'), ('passive', 'v')):
                rr = sorted((r for r in rows if r[1] == st))
                if rr:
                    a.plot([r[0] for r in rr], [r[2] for r in rr], mk, ms=7, mfc='none', mec=col,
                           mew=1.3, label=f'{lab}' if st == 'active' else None, zorder=4)
        if law == 'hill':
            for N, mk, col in ((1024, 's', '#3D7EA6'), (4096, 'D', '#21918c')):
                p = out/f'cells_hill_N{N}_tchi0.3tc.json'
                if p.exists():
                    runs = json.loads(p.read_text())['runs']
                    a.plot([r['tau_m'] for r in runs], [r['tail_mean'] for r in runs], mk, ms=4.5,
                           color=col, alpha=.8, label=f'toy N = {N}, τχ = 0.3 τc', zorder=5)
        a.set(xlabel='τm / τc', ylabel='late active fraction x = ⟨χ⟩', ylim=(-.03, 1.03),
              title=('hill law (n₁ = 16, wp = 0.002)' if law == 'hill' else 'step law (legacy)')
              + f'; c_D = {mp["c_D"]:.2f}, β = 0')
        a.legend(frameon=False, fontsize=7.5, loc='center right')
    # F(x) - x at a few memory times, toy against closure
    a = ax[0, 2]
    if (out/'map_hill.json').exists():
        mp = json.loads((out/'map_hill.json').read_text())
        cmap = plt.get_cmap('viridis')
        pick = [t for t in (2., 5., 7., 8., 9., 10., 12., 16.) if t in mp['tau_m']]
        for k, tm in enumerate(pick):
            i = mp['tau_m'].index(tm); col = cmap(k/max(len(pick)-1, 1))
            a.plot(mp['x'], np.array(mp['F'][i]) - np.array(mp['x']), color=col, lw=1.6, label=f'{tm:g}')
            a.plot(mp['x_closure'], np.array(mp['F_closure'][i]) - np.array(mp['x_closure']),
                   color=col, lw=.9, ls='--')
        a.axhline(0, color='#30363B', lw=.7)
        a.set(xlabel='x', ylabel='F(x) − x', title='hill: toy (solid) vs Beta closure (dashed)')
        a.legend(title='τm/τc', frameon=False, fontsize=7.5, ncol=2)
    # finite-N traces and the tau_chi comparison
    for j, tchi in enumerate(('0.3tc', 'tm')):
        a = ax[1, j]
        p = out/f'cells_hill_N4096_tchi{tchi}.json'
        if not p.exists():
            continue
        runs = json.loads(p.read_text())['runs']
        grid = sorted({r['tau_m'] for r in runs})
        cmap = plt.get_cmap('viridis')
        for r in runs:
            col = cmap(grid.index(r['tau_m'])/max(len(grid)-1, 1))
            a.plot(r['series_t'], r['series_x'], color=col, lw=.8, ls='-' if r['start'] == 'active' else '--')
        a.set(xlabel='t / τc after release', ylabel='x(t)', ylim=(-.03, 1.03),
              title=f'toy N = 4096, τχ = {"0.3 τc" if tchi == "0.3tc" else "τm"} (dashed: passive start)')
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, len(grid)-1))
        cb = fig.colorbar(sm, ax=a, ticks=range(len(grid)), pad=.01)
        cb.ax.set_yticklabels([f'{g:.3g}' for g in grid], fontsize=7); cb.set_label('τm/τc')
    a = ax[1, 2]
    for tchi, col in (('0.3tc', '#E07A5F'), ('tm', '#8C6BB1')):
        for N, mk in ((1024, 's'), (4096, 'D')):
            p = out/f'cells_hill_N{N}_tchi{tchi}.json'
            if p.exists():
                runs = json.loads(p.read_text())['runs']
                a.plot([r['tau_m'] for r in runs], [r['tail_mean'] for r in runs], mk, ms=4.5,
                       color=col, mfc=col if N == 4096 else 'none',
                       label=f'N = {N}, τχ = {"0.3 τc" if tchi == "0.3tc" else "τm"}')
    a.set(xlabel='τm / τc', ylabel='late x', ylim=(-.03, 1.03), title='toy: effect of τχ (β = 0)')
    a.legend(frameon=False, fontsize=7.5)
    fig.suptitle('Non-spatial toy: OU pressure (open-loop μ, σ_P, τ_p), one-pole memory, '
                 f'constant c_D = {cal["c_D"]:.2f} (fitted on six open-loop rungs), m_c = {MC}', fontsize=11)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(out/f'toy_meanfield.{ext}', dpi=180, facecolor='white')
    plt.close(fig)
    print(json.dumps({'figure': str(out/'toy_meanfield.png')}))


# ---- driver ------------------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('mode', choices=('calibrate', 'map', 'cells', 'figure'))
    ap.add_argument('out', type=Path)
    ap.add_argument('--law', choices=tuple(LAWS), default='hill')
    ap.add_argument('--N', type=int, default=1024)
    ap.add_argument('--tchi', choices=('0.3tc', 'tm'), default='0.3tc')
    ap.add_argument('--grid', choices=('scan', 'legacy'), default='scan',
                    help='cells: the hill-scan grid (12 points) or the 20-point legacy grid')
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    if a.mode == 'figure':
        figure(a.out)
        return
    inp = Inputs()
    if a.mode == 'calibrate':
        res = calibrate(inp)
        (a.out/'calibration.json').write_text(json.dumps(res, indent=1))
        print(json.dumps({k: res[k] for k in ('c_D', 'c_D_per_rung', 'c_D_log_spread')}))
        return
    cD = json.loads((a.out/'calibration.json').read_text())['c_D']
    law = Law(LAWS[a.law]['h-n'], LAWS[a.law]['wp'])
    if a.mode == 'map':
        tms = np.r_[0.3, np.arange(1., 6., 1.), np.arange(6., 14.01, .5), 15., 16., 18., 20.]
        xs = np.linspace(0, 1, 41)
        F = fmap(inp, law, tms, cD, xs)
        Fc = closure_map(inp, law, tms, cD, np.linspace(0, 1, 401))
        res = {'law': a.law, **LAWS[a.law], 'mc': MC, 'pmem': PMEM, 'c_D': cD,
               'tau_m': tms.tolist(), 'x': xs.tolist(), 'F': F.tolist(),
               'x_closure': np.linspace(0, 1, 401).tolist(), 'F_closure': Fc.tolist(),
               'roots': [fixed_points(xs, f) for f in F],
               'roots_closure': [fixed_points(np.linspace(0, 1, 401), f) for f in Fc]}
        (a.out/f'map_{a.law}.json').write_text(json.dumps(res))
        for tm, r, rc in zip(tms, res['roots'], res['roots_closure']):
            print(f'tau_m={tm:5.2f} toy {[(round(q["x"],3), q["stable"]) for q in r]}  '
                  f'closure {[(round(q["x"],3), q["stable"]) for q in rc]}')
        return
    if a.mode == 'cells':
        grid = ([0.3, 2.163157895, 4.957894737, 6., 7., 7.752631579, 8.684210526, 9.6,
                 10.5, 12.410526316, 14., 18.] if a.grid == 'scan' else
                [round(0.3 + k*(18 - 0.3)/19, 9) for k in range(20)])
        rows = [{'tau_m': tm, 'tau_chi': TCHI0 if a.tchi == '0.3tc' else tm, 'start': st}
                for tm in grid for st in ('active', 'passive')]
        res = cells_run(inp, law, rows, a.N, cD)
        name = f'cells_{a.law}_N{a.N}_tchi{a.tchi}' + ('_legacygrid' if a.grid == 'legacy' else '') + '.json'
        (a.out/name).write_text(json.dumps({'law': a.law, 'N': a.N, 'tchi': a.tchi, 'c_D': cD,
                                            'runs': res}))
        for r in res:
            print(f'tau_m={r["tau_m"]:6.3f} {r["start"]:7s} x={r["tail_mean"]:.3f} '
                  f'm={r["tail_m"]:.3f} drift={r["half_drift"]:+.3f}')
        return


if __name__ == '__main__':
    main()
