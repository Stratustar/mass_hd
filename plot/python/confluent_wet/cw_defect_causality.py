#!/usr/bin/env python3
"""Does the active | passive interface follow the +1/2 defects, or do the defects follow the interface?

  cw_defect_causality.py case <run.npz> <out_dir> [--tau-chi-tc X]   per-case statistics -> <out_dir>/<case>.json
  cw_defect_causality.py plot <out_dir>                              figures from every <case>.json in out_dir

Input: a dense npz (cw_frames_to_npz.py; chi, m, ux_mat, uy_mat, QQxx, QQyx every 0.2 tau_c). The interface
is chi smoothed (Gaussian sigma 1, periodic) at 1/2, s = signed distance (+ inside active), as in
cw_defect_interface_relative.py. Units grid and tau_c. The two arrows have very different delays:
defect -> interface runs through P -> m (tau_m) -> chi (tau_chi); interface -> defect through
zeta_eff(chi) -> active stress, which is instantaneous.

A  interface speed against distance d to the nearest +1/2 (and -1/2) core, interface pixels |s| <= 1:
   G = [s_{i+2} - s_i]/dt at a fixed pixel (+ = active grows), and the instantaneous budget of the
   smoothed field, G = (d_t chi_s)/|grad chi_s| = G_adv + G_dif + G_sw with
       d_t chi = -u.grad chi + Dchi lap chi + (h(m) - chi)/tau_chi      (beta = 0: chi* = h(m))
   with the measured (chi_s(i+1) - chi_s(i-1))/(2 dt)/|grad chi_s| ('dchi') as the closure check;
   split again by the class of the nearest +1/2 (HA head active / tail passive, HP the reverse,
   rays s = 3..7 as before) for d <= BAND. Also the enrichment of +1/2 cores against s.
B  wake: for +1/2 cores every START_STRIDE frames, chi and m at FIXED lab points of the head frame
   (head at 90 deg, r = 1..12, 10 deg) at lag 0 ... 20 tau_c later. Excess over a null of random
   points and orientations with the same starting value (binned), in three regions: core r <= 3,
   head (r 3..7 within 30 deg of the head), tail (same behind). Classes all / near (|s| <= BAND) / HA / HP.
C  lagged field cross-correlation C(tau) = <drho(x, t) dX(x, t + tau)> / (std std), rho = +1/2 (-1/2)
   core density (Gaussian sigma 3), X = |grad chi_s| (interface density), chi, m; d = deviation from
   the frame mean. tau > 0: X after the defects (interface follows); tau < 0: X before (defects follow).
"""
import argparse
import json
from pathlib import Path
import re
import sys

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cw_defect_interface_relative import BAND, S_RAY, TC, at, grad, signed_distance  # noqa: E402

DCHI, MC, HN = 0.00166667, 0.2287, 16
LAGS = (0, 1, 2, 5, 10, 25, 50, 75, 100)          # frames of 0.2 tau_c: 0 ... 20 tau_c
START_STRIDE = 5
RW = np.arange(1.0, 12.0 + 1e-9, 1.0)
TW = np.radians(np.arange(5.0, 360.0, 10.0))
DBINS = np.array([0, 1, 2, 3, 4, 6, 8, 10, 13, 16, 20, 25, 32, 64.])
SBINS = np.arange(-20.0, 20.0 + 1e-9, 1.0)
NULL, CC_MAX, SIG_RHO = 1500, 100, 3.0
RG, TG = np.meshgrid(RW, np.degrees(TW), indexing='ij')
REGIONS = {'core': RG <= 3, 'head': (RG >= 3) & (RG <= 7) & (np.abs(TG - 90) <= 30),
           'tail': (RG >= 3) & (RG <= 7) & (np.abs(TG - 270) <= 30)}
VBINS = {'chi': np.linspace(0, 1, 21), 'm': np.linspace(0, .5, 26)}


def lap(f):
    return np.roll(f, 1, 0) + np.roll(f, -1, 0) + np.roll(f, 1, 1) + np.roll(f, -1, 1) - 4*f


def cores(qxx, qxy):
    th = .5*np.arctan2(qxy, qxx)
    d = lambda a, b: (b - a + np.pi/2) % np.pi - np.pi/2
    t10 = np.roll(th, -1, 0); t11 = np.roll(t10, -1, 1); t01 = np.roll(th, -1, 1)
    w = (d(th, t10) + d(t10, t11) + d(t11, t01) + d(t01, th))/(2*np.pi)
    return np.argwhere(w > .25) + .5, np.argwhere(w < -.25) + .5


def heads(qxx, qxy, plus):
    a, b = grad(qxx); c, e = grad(qxy)
    px, py = at(a + e, plus), at(c - b, plus)
    return -np.c_[px, py]/np.maximum(np.hypot(px, py), 1e-300)[:, None]


def classes(chi, p, h):
    hh = at(chi, p[:, None, :] + S_RAY[None, :, None]*h[:, None, :]).mean(1)
    tt = at(chi, p[:, None, :] - S_RAY[None, :, None]*h[:, None, :]).mean(1)
    return np.where((hh > .5) & (tt <= .5), 1, np.where((hh <= .5) & (tt > .5), 2, 0))   # 1 HA, 2 HP


def density(p, L):
    r = np.zeros((L, L), np.float32)
    if len(p):
        ij = np.floor(p).astype(int) % L
        np.add.at(r, (ij[:, 0], ij[:, 1]), 1)
    return gaussian_filter(r, SIG_RHO, mode='wrap')


def polar_points(p, h):
    e = np.c_[h[:, 1], -h[:, 0]]
    c, s = np.cos(TW), np.sin(TW)
    return (p[:, None, None, :] + RW[None, :, None, None]*(c[None, None, :, None]*e[:, None, None, :]
                                                          + s[None, None, :, None]*h[:, None, None, :]))


def hill(m):
    return np.where(m > 0, 1/(1 + (MC/np.maximum(m, 1e-12))**HN), 0.)


def case(npz, out, tau_chi_tc):
    name = Path(npz).stem
    if tau_chi_tc is None:
        tau_chi_tc = float(re.search(r'_tc([0-9p]+)_', name).group(1).replace('p', '.'))
    tau_chi = 202.3 if tau_chi_tc == .3 else tau_chi_tc*TC          # steps, as the cards
    z = np.load(npz)
    f32 = lambda k: z[k].astype(np.float32)
    chi, m, ux, uy, qxx, qxy = f32('chi'), f32('m'), f32('ux_mat'), f32('uy_mat'), f32('QQxx'), f32('QQyx')
    T, L = chi.shape[0], chi.shape[1]
    dt = float(z['t'][1] - z['t'][0])/TC
    D = [cores(qxx[i].astype(float), qxy[i].astype(float)) for i in range(T)]
    H = [heads(qxx[i].astype(float), qxy[i].astype(float), D[i][0]) for i in range(T)]
    S = np.stack([signed_distance(chi[i].astype(float)) for i in range(T)])
    chis = np.stack([gaussian_filter(chi[i], 1, mode='wrap') for i in range(T)])
    gn = np.stack([np.hypot(*grad(chis[i])) for i in range(T)])
    textured = bool(np.isfinite(S).all())
    print(name, 'frames', T, 'dt', round(dt, 4), 'textured', textured, 'mean +1/2', np.mean([len(d[0]) for d in D]), flush=True)
    res = {'case': name, 'dt_tc': dt, 'tau_chi_tc': tau_chi_tc, 'textured': textured, 'chi_mean': float(chi.mean()),
           'chi_var_x': float(chi.var(axis=(1, 2)).mean()), 'n_plus_frame': float(np.mean([len(d[0]) for d in D]))}

    # ---- A: interface speed vs distance to the nearest core, and its budget ----
    nb = len(DBINS) - 1
    keys = ('G', 'dchi', 'adv', 'dif', 'sw')
    acc = {sgn: {k: np.zeros(nb) for k in keys + ('n',)} for sgn in ('plus', 'minus')}
    cls_acc = {c: {k: 0. for k in keys + ('n',)} for c in ('HA', 'HP', 'other')}
    glob = {k: 0. for k in keys + ('n',)}
    hs_def, hs_all = np.zeros(len(SBINS) - 1), np.zeros(len(SBINS) - 1)
    for i in range(1, T - 2) if textured else ():
        c = chi[i].astype(float)
        gx, gy = grad(c)
        terms = {'adv': -(ux[i]*gx + uy[i]*gy)*TC, 'dif': DCHI*lap(c)*TC, 'sw': (hill(m[i].astype(float)) - c)/tau_chi*TC}
        mask = (np.abs(S[i]) <= 1) & (gn[i] > .02)
        pix = np.argwhere(mask).astype(float)
        vals = {'G': ((S[i + 2] - S[i])/(2*dt))[mask], 'dchi': ((chis[i + 1] - chis[i - 1])/(2*dt)/np.maximum(gn[i], 1e-12))[mask]}
        for k, f in terms.items():
            vals[k] = (gaussian_filter(f, 1, mode='wrap')/np.maximum(gn[i], 1e-12))[mask]
        for k in keys:
            glob[k] += vals[k].sum()
        glob['n'] += mask.sum()
        for sgn, p in (('plus', D[i][0]), ('minus', D[i][1])):
            if not len(p):
                continue
            d, j = cKDTree(p, boxsize=L).query(pix % L)
            b = np.clip(np.digitize(d, DBINS) - 1, 0, nb - 1)
            for k in keys:
                acc[sgn][k] += np.bincount(b, vals[k], nb)
            acc[sgn]['n'] += np.bincount(b, minlength=nb)
            if sgn == 'plus':
                lab = classes(c, p, H[i])[j]
                near = d <= BAND
                for cl, sel in (('HA', near & (lab == 1)), ('HP', near & (lab == 2)), ('other', near & (lab == 0))):
                    for k in keys:
                        cls_acc[cl][k] += vals[k][sel].sum()
                    cls_acc[cl]['n'] += sel.sum()
        if len(D[i][0]):
            hs_def += np.histogram(at(S[i], D[i][0]), SBINS)[0]
        hs_all += np.histogram(S[i], SBINS)[0]
    if textured:
        res['interface'] = {
            'global': {k: float(glob[k]/max(glob['n'], 1)) for k in keys} | {'n': int(glob['n'])},
            'd_edges': DBINS.tolist(),
            **{sgn: {k: (acc[sgn][k]/np.maximum(acc[sgn]['n'], 1)).tolist() for k in keys} | {'n': acc[sgn]['n'].astype(int).tolist()}
               for sgn in acc},
            'by_class': {cl: {k: float(a[k]/max(a['n'], 1)) for k in keys} | {'n': int(a['n'])} for cl, a in cls_acc.items()},
            's_edges': SBINS.tolist(),
            'enrichment': (hs_def/max(hs_def.sum(), 1)/np.maximum(hs_all/max(hs_all.sum(), 1), 1e-12)).tolist(),
            'n_def_s': hs_def.astype(int).tolist()}

    # ---- B: wake at fixed lab points, excess over a matched null ----
    rng = np.random.default_rng(11)
    fields = {'chi': chi, 'm': m}
    wsum = {cl: {f: np.zeros((len(LAGS), len(RW), len(TW))) for f in fields} for cl in ('all', 'near', 'HA', 'HP')}
    wn = {cl: np.zeros(len(LAGS), int) for cl in wsum}
    dreg = {(f, r, l): [] for f in fields for r in REGIONS for l in range(len(LAGS))}   # (start value, change, class code)
    nreg = {(f, r, l): [] for f in fields for r in REGIONS for l in range(len(LAGS))}
    for i in range(0, T, START_STRIDE):
        p, h = D[i][0], H[i]
        ang = rng.uniform(0, 2*np.pi, NULL)
        q, hq = rng.uniform(0, L, (NULL, 2)), np.c_[np.cos(ang), np.sin(ang)]
        lab = classes(chi[i].astype(float), p, h) if len(p) else np.zeros(0, int)
        near = np.abs(at(S[i], p)) <= BAND if (len(p) and textured) else np.zeros(len(p), bool)
        sel = {'all': np.ones(len(p), bool), 'near': near, 'HA': near & (lab == 1), 'HP': near & (lab == 2)}
        code = np.where(sel['HA'], 1, np.where(sel['HP'], 2, np.where(near, 3, 0)))
        xy, xq = polar_points(p, h), polar_points(q, hq)
        for f, F in fields.items():
            v0, q0 = at(F[i], xy), at(F[i], xq)
            for l, lag in enumerate(LAGS):
                if i + lag >= T:
                    continue
                dv, dq = at(F[i + lag], xy) - v0, at(F[i + lag], xq) - q0
                for cl, s_ in sel.items():
                    wsum[cl][f][l] += dv[s_].sum(0)
                    if f == 'chi':
                        wn[cl][l] += int(s_.sum())
                for r, msk in REGIONS.items():
                    dreg[(f, r, l)].append(np.c_[v0[:, msk].mean(1), dv[:, msk].mean(1), code])
                    nreg[(f, r, l)].append(np.c_[q0[:, msk].mean(1), dq[:, msk].mean(1)])
    wake = {'lags_tc': [lag*dt for lag in LAGS], 'r': RW.tolist(), 'theta_deg': np.degrees(TW).tolist(),
            'n': {cl: wn[cl].tolist() for cl in wn},
            'maps': {cl: {f: (wsum[cl][f]/np.maximum(wn[cl], 1)[:, None, None]).tolist() for f in fields} for cl in wsum},
            'excess': {}}
    for (f, r, l), rows in dreg.items():
        if not rows:
            continue
        A, N = np.concatenate(rows), np.concatenate(nreg[(f, r, l)])
        e = VBINS[f]
        bn, ba = np.clip(np.digitize(N[:, 0], e) - 1, 0, len(e) - 2), np.clip(np.digitize(A[:, 0], e) - 1, 0, len(e) - 2)
        cnt = np.bincount(bn, minlength=len(e) - 1)
        mean_null = np.where(cnt > 20, np.bincount(bn, N[:, 1], len(e) - 1)/np.maximum(cnt, 1), np.nan)
        x = A[:, 1] - mean_null[ba]
        for cl, s_ in (('all', np.ones(len(A), bool)), ('near', A[:, 2] > 0), ('HA', A[:, 2] == 1), ('HP', A[:, 2] == 2)):
            ok = s_ & np.isfinite(x)
            d = wake['excess'].setdefault(cl, {}).setdefault(f, {}).setdefault(r, {'mean': [], 'se': [], 'raw': [], 'n': []})
            d['mean'].append(float(x[ok].mean()) if ok.any() else None)
            d['se'].append(float(x[ok].std()/np.sqrt(ok.sum())) if ok.sum() > 1 else None)
            d['raw'].append(float(A[ok, 1].mean()) if ok.any() else None)
            d['n'].append(int(ok.sum()))
    res['wake'] = wake

    # ---- C: lagged field cross-correlations ----
    def flat(F):
        X = F.reshape(T, -1).astype(np.float32)
        X = X - X.mean(1, keepdims=True)
        return X/max(float(X.std()), 1e-30)
    rp = flat(np.stack([density(D[i][0], L) for i in range(T)]))
    rm = flat(np.stack([density(D[i][1], L) for i in range(T)]))
    X = {'iface': flat(gn), 'chi': flat(chi), 'm': flat(m)}
    lags = np.arange(-CC_MAX, CC_MAX + 1)
    def cc(A, B):
        M = (A @ B.T)/A.shape[1]           # M[t, t'] = <A(t) B(t')>
        return [float(np.diagonal(M, k).mean()) for k in lags]
    res['xcorr'] = {'lags_tc': (lags*dt).tolist(),
                    **{f'plus|{k}': cc(rp, v) for k, v in X.items()},
                    **{f'minus|{k}': cc(rm, v) for k, v in X.items()},
                    'plus|plus': cc(rp, rp), 'iface|iface': cc(X['iface'], X['iface']), 'chi|chi': cc(X['chi'], X['chi'])}
    out.mkdir(parents=True, exist_ok=True)
    (out/f'{name}.json').write_text(json.dumps(res))
    if textured:
        I = res['interface']
        print('global G', {k: round(v, 3) for k, v in I['global'].items()})
        print('G vs d(+1/2)', np.round(I['plus']['G'][:8], 3), 'sw', np.round(I['plus']['sw'][:8], 3), 'adv', np.round(I['plus']['adv'][:8], 3))
        print('by class', {c: {k: round(v, 3) for k, v in a.items()} for c, a in I['by_class'].items()})
    ex = res['wake']['excess']['all']['chi']
    print('wake chi excess core', [None if v is None else round(v, 4) for v in ex['core']['mean']])
    print('xcorr plus|iface at -20, -5, 0, 5, 20 tc',
          [round(res['xcorr']['plus|iface'][CC_MAX + k], 4) for k in (-100, -25, 0, 25, 100)])


def plot(out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    J = [json.loads(f.read_text()) for f in sorted(out.glob('hn*.json'))]
    J.sort(key=lambda j: ('open' in j['case'], j['case']))
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9.5, 'xtick.direction': 'in', 'ytick.direction': 'in',
                         'legend.frameon': False, 'figure.facecolor': 'white', 'savefig.facecolor': 'white', 'pdf.fonttype': 42})
    title = lambda j: j['case'].replace('hn16_wp0p002_', '').replace('_rep', ' rep') + f"\n⟨χ⟩ {j['chi_mean']:.2f}"
    n = len(J)
    fig, ax = plt.subplots(4, n, figsize=(3.6*n, 13), squeeze=False)
    for c, j in enumerate(J):
        a = ax[0, c]
        if j['textured']:
            I = j['interface']; e = np.array(I['d_edges']); x = .5*(e[1:] + e[:-1]); x[-1] = 40
            ok = np.array(I['plus']['n']) >= 200
            for k, col, ls in (('G', 'k', '-'), ('dchi', '#6a3d9a', '-'), ('adv', '#2a9d8f', '--'), ('sw', '#c8382f', '--'), ('dif', '0.6', ':')):
                a.plot(x[ok], np.array(I['plus'][k])[ok], ls, color=col, marker='o' if k == 'G' else None, ms=3,
                       label={'G': 'G from s', 'dchi': '∂tχ/|∇χ| measured', 'adv': 'advection', 'sw': 'switching', 'dif': 'diffusion'}[k])
            a.plot(x[ok], np.array(I['minus']['G'])[ok], '-', color='0.55', lw=.8, label='G vs −1/2')
            a.axhline(I['global']['G'], color='k', lw=.6, ls=':')
            a.set_xscale('log')
        a.axhline(0, color='0.7', lw=.6)
        a.set_xlabel('distance to nearest +1/2 (grid)'); a.set_title(title(j), fontsize=9)
        if c == 0:
            a.set_ylabel('interface speed (grid/τc)\n+ = active grows'); a.legend(fontsize=7)
        W = j['wake']; t = np.array(W['lags_tc'])
        for row, f in ((1, 'chi'), (2, 'm')):
            a = ax[row, c]
            for r, col in (('core', 'k'), ('head', '#c8382f'), ('tail', '#1f5fa8')):
                y = np.array(W['excess']['all'][f][r]['mean'], float)
                a.plot(t[:len(y)], y, '-o', ms=3, color=col, label=f'all: {r}')
            if f == 'chi' and j['textured']:
                for cl, col in (('HA', '#c8382f'), ('HP', '#1f5fa8')):
                    y = np.array(W['excess'][cl][f]['core']['mean'], float)
                    a.plot(t[:len(y)], y, '--', color=col, lw=1, label=f'{cl}: core')
            a.axhline(0, color='0.7', lw=.6); a.set_xlabel('lag (τc)')
            if c == 0:
                a.set_ylabel(f'excess Δ{"χ" if f == "chi" else "m"} at fixed points\n(over matched random points)')
                a.legend(fontsize=7)
        a = ax[3, c]; X = j['xcorr']; t = np.array(X['lags_tc'])
        for k, col in (('plus|iface', 'k'), ('plus|chi', '#c8382f'), ('plus|m', '#e3a51a'), ('minus|iface', '0.6')):
            y = np.array(X[k]); a.plot(t, y, color=col, label=k.replace('|', ' × ').replace('plus', 'ρ+').replace('minus', 'ρ−'))
            a.plot(-t, y, color=col, ls=':', lw=.8)
        a.axhline(0, color='0.7', lw=.6); a.axvline(0, color='0.7', lw=.6)
        a.set_xlabel('τ (τc): X after defects  →'); a.set_xlim(t.min(), t.max())
        if c == 0:
            a.set_ylabel('C(τ) = ⟨δρ(t) δX(t+τ)⟩\n(dotted: mirrored, C(−τ))'); a.legend(fontsize=7)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(out/f'causality_summary.{ext}', dpi=180)
    plt.close(fig)

    show = (1, 5, 25, 50, 100)
    fig = plt.figure(figsize=(2.3*len(show) + 1, 2.4*n))
    gs = fig.add_gridspec(n, len(show) + 1, width_ratios=[1]*len(show) + [.06], wspace=.15, hspace=.35)
    for r_, j in enumerate(J):
        W = j['wake']; M = np.array(W['maps']['all']['chi'])
        idx = [LAGS.index(s) for s in show]
        v = max(np.abs(M[idx]).max(), 1e-6)
        et = np.radians(np.arange(0, 365, 10)); er = np.r_[RW - .5, RW[-1] + .5]
        for k, l in enumerate(idx):
            a = fig.add_subplot(gs[r_, k], projection='polar')
            mm = a.pcolormesh(et, er, M[l], cmap='RdBu_r', vmin=-v, vmax=v, shading='flat')
            a.set_theta_zero_location('E'); a.set_xticks([]); a.set_yticks([]); a.set_ylim(0, er[-1])
            a.set_title(f'lag {W["lags_tc"][l]:.1f} τc' + (f'\n{j["case"].replace("hn16_wp0p002_", "")}' if k == 0 else ''), fontsize=8)
        cb = fig.colorbar(mm, cax=fig.add_subplot(gs[r_, -1])); cb.set_label('Δχ', fontsize=8)
    for ext in ('png', 'pdf'):
        fig.savefig(out/f'wake_maps.{ext}', dpi=180, bbox_inches='tight')
    print('figures', out/'causality_summary.png', out/'wake_maps.png', n, 'cases')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('mode', choices=('case', 'plot'))
    ap.add_argument('paths', nargs='+', type=Path)
    ap.add_argument('--tau-chi-tc', type=float)
    a = ap.parse_args()
    if a.mode == 'case':
        case(a.paths[0], a.paths[1], a.tau_chi_tc)
    else:
        plot(a.paths[0])
