#!/usr/bin/env python3
"""Generate and check the local validation of the confluent-wet hill response law.

generate <dir>        write the run.dat cards
check <outroot>       compare the simulation outputs (<outroot>/<name>[_ref|_new]/)

Cases
  legacy_step, legacy_tanh   legacy law, run on the reference (c144403) AND the new binary:
                             every output byte must agree (bit-for-bit legacy).
  hill_mirror                hill law at its step limits (h-n = 0, wp = 0, beta = 0) with the
                             phenotype flipped: chi_hill = 1 - chi_legacy to roundoff.
  hill_mirror_split          the same with Dbio = 0, Dm = Dchi = the old Dbio: bit-identical
                             to hill_mirror.
  ode_*                      uniform, flow-free layer: chi(t), m(t) against a Python replica
                             of the Heun stages (incl. clamps) to machine precision.
  diff_split                 pure diffusion of a stripe with Dm != Dchi against a 1D replica.
  bad_*                      parameter combinations that must be rejected at startup.
"""
import argparse
import json
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[1]
BASE = REPO/'cases/20260911/cw_mc02287_four_init/L256/mc0p2287_tm0p3_chi0_rep1/run.dat'


def base_card():
    p = {}
    for line in BASE.read_text().splitlines():
        if line.strip() and not line.startswith('#'):
            k, v = line.split('=', 1)
            p[k.strip()] = v.strip()
    return p


def turbulent(**kw):
    q = base_card()
    q.update({'LX': 64, 'LY': 64, 'nsteps': 1500, 'ninfo': 500, 'seed': 7351,
              'chi-config': 'binary-noise', 'chi-length': 6., 'chi-seed': 200901,
              'chi-lo': 0, 'm-lo': .3301, 'chi-hi': 1, 'm-hi': .0559,
              'chi-freeze-steps': 0, 'tau-m': 200., 'tau-chi': 60., 'nresponse': 25,
              'nvideo': 100, 'video-stride': 4, 'frame-light': 0})
    q.update(kw)
    return q


def uniform(**kw):
    q = base_card()
    q.update({'LX': 16, 'LY': 16, 'nsteps': 600, 'ninfo': 600, 'seed': 7351, 'noise': 0.,
              'chi-config': 'uniform', 'chi-freeze-steps': 0, 'mem-freeze-steps': 0,
              'nresponse': 1, 'nvideo': 0, 'frame-light': 0,
              'switch-model': 'hill', 'switch-sign': 1, 'tau-m': 40., 'tau-chi': 15.})
    q.update(kw)
    return q


HILL_STEP = {'switch-model': 'hill', 'switch-sign': 1, 'h-n': 0, 'wp': 0, 'beta': 0}

CASES = {
    'legacy_step': turbulent(),
    'legacy_tanh': turbulent(**{'chi-width': .05, 'pmem-width': .005, 'switch-sign': 1,
                                'pmem-pulse-value': 0., 'pmem-pulse-start': 300,
                                'pmem-pulse-steps': 100}),
    'hill_mirror': turbulent(**HILL_STEP, **{'chi-lo': 1, 'chi-hi': 0}),
    'hill_mirror_split': turbulent(**HILL_STEP, **{'chi-lo': 1, 'chi-hi': 0, 'Dbio': 0,
                                                   'Dm': .00166667, 'Dchi': .00166667}),
    # m decays from 0.6 towards s(0) = 1/(1+e^5): chi is switched on, then must stay locked
    'ode_active_irrev': uniform(**{'beta': .95, 'h-n': 4, 'mc': .3, 'g-n': 2, 'g-chic': .5,
                                   'g-lock': 'active', 'chi0': 0., 'm0': .6,
                                   'pmem': .05, 'wp': .01}),
    'ode_active_rev': uniform(**{'beta': .6, 'h-n': 4, 'mc': .3, 'g-n': 2, 'g-chic': .5,
                                 'g-lock': 'active', 'chi0': 0., 'm0': .6,
                                 'pmem': .05, 'wp': .01}),
    # m rises from exactly 0 towards s(0) = 1/(1+e^-1): chi must stay locked passive
    'ode_passive': uniform(**{'beta': .9, 'h-n': 6, 'mc': .3, 'g-n': 2.5, 'g-chic': .4,
                              'g-lock': 'passive', 'chi0': 0., 'm0': 0.,
                              'pmem': -.01, 'wp': .01}),
    # both hard steps: s = 1 (P = 0 > pmem), m crosses mc from below, h jumps 0 -> 1
    'ode_steps': uniform(**{'beta': .5, 'h-n': 0, 'mc': .3, 'g-n': 3, 'g-chic': .6,
                            'g-lock': 'active', 'chi0': .3, 'm0': .2,
                            'pmem': -.01, 'wp': 0}),
    'diff_split': uniform(**{'LX': 64, 'LY': 8, 'nsteps': 300, 'nresponse': 10,
                             'zeta': 0., 'chi-config': 'stripe', 'chi-hi': 1, 'chi-lo': 0,
                             'm-hi': .8, 'm-lo': .2, 'tau-chi': 0., 'mem-freeze-steps': 300,
                             'Dbio': 0, 'Dm': .05, 'Dchi': .01, 'h-n': 4, 'mc': .3}),
    'bad_legacy_beta': turbulent(**{'nsteps': 10, 'beta': .5}),
    'bad_hill_sign': turbulent(**{**HILL_STEP, 'nsteps': 10, 'switch-sign': -1}),
    'bad_beta': uniform(**{'nsteps': 10, 'beta': 1.2}),
    'bad_gn': uniform(**{'nsteps': 10, 'g-n': 1.}),
    'bad_D': turbulent(**{'nsteps': 10, 'Dm': .01}),
}
REFERENCE = ('legacy_step', 'legacy_tanh')   # also run on the reference binary


def generate(out):
    for name, q in CASES.items():
        path = out/name/'run.dat'
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('# Local validation of switch-model=hill; not a physics result.\n' +
                        ''.join(f'{k} = {v}\n' for k, v in q.items()))
    print(json.dumps({'generated': list(CASES), 'reference': list(REFERENCE),
                      'out': str(out)}))


# ---- replicas -----------------------------------------------------------------------

def hill_fns(q):
    mc, hn = float(q['mc']), float(q['h-n'])
    beta, gn, gc = float(q['beta']), float(q['g-n']), float(q['g-chic'])
    gnorm = 1 + gc**gn
    passive = q['g-lock'] == 'passive'

    def h(m):
        if hn == 0:
            return 1. if m >= mc else 0.
        return 1/(1 + (mc/m)**hn) if m > 0 else 0.

    def g(x):
        return gnorm/(1 + (gc/x)**gn) if x > 0 else 0.

    def target(m, c):
        H = h(m)
        return H*(1 - beta*g(1 - c)) if passive else H + beta*g(c)*(1 - H)

    pmem, wp = float(q['pmem']), float(q['wp'])
    s = 1/(1 + np.exp(-(0. - pmem)/wp)) if wp > 0 else float(0. > pmem)
    return target, s


def ode_replica(q):
    target, s = hill_fns(q)
    tm, tc = float(q['tau-m']), float(q['tau-chi'])
    clamp = lambda x: min(1., max(0., x))
    m, c = float(q['m0']), float(q['chi0'])
    out = [(m, c)]
    for _ in range(int(q['nsteps'])):
        rm, rc = (s - m)/tm, (target(m, c) - c)/tc
        mN, cN = m + .5*rm, c + .5*rc
        mt, ct = clamp(m + rm), clamp(c + rc)
        rm2, rc2 = (s - mt)/tm, (target(mt, ct) - ct)/tc
        m, c = clamp(mN + .5*rm2), clamp(cN + .5*rc2)
        out.append((m, c))
    return np.array(out)


def diffusion_replica(q):
    L, n = int(q['LX']), int(q['nsteps'])
    x = np.arange(L)
    res = {}
    for fld, D, hi, lo in (('chi', float(q['Dchi']), float(q['chi-hi']), float(q['chi-lo'])),
                           ('m', float(q['Dm']), float(q['m-hi']), float(q['m-lo']))):
        u = np.where(x < L//2, hi, lo).astype(float)
        hist = [u.copy()]
        lap = lambda v: np.roll(v, 1) + np.roll(v, -1) - 2*v
        for _ in range(n):
            r = D*lap(u)
            uN, ut = u + .5*r, np.clip(u + r, 0, 1)
            u = np.clip(uN + .5*D*lap(ut), 0, 1)
            hist.append(u.copy())
        res[fld] = np.array(hist)
    return res


# ---- checks -------------------------------------------------------------------------

def response(d):
    return np.genfromtxt(d/'response.csv', delimiter=',', names=True)


def frames(d):
    return sorted(p.name for p in d.glob('frame*.json'))


def frame(d, name):
    j = json.loads((d/name).read_text())['data']
    return {k: np.asarray(v['value'], dtype=float) for k, v in j.items()}


def same_bytes(a, b):
    names = sorted(p.name for p in a.iterdir() if p.name != 'parameters.json')
    assert names == sorted(p.name for p in b.iterdir() if p.name != 'parameters.json'), \
        (names, sorted(p.name for p in b.iterdir()))
    return {n: (a/n).read_bytes() == (b/n).read_bytes() for n in names}


def check(root):
    out = {}

    # 1. legacy: reference binary vs new binary, every output byte
    for name in REFERENCE:
        eq = same_bytes(root/f'{name}_ref', root/f'{name}_new')
        assert all(eq.values()), eq
        out[name] = {'identical_files': sorted(eq)}

    # 2. mirror: hill step limit vs legacy step, chi -> 1 - chi
    a, b = response(root/'legacy_step_new'), response(root/'hill_mirror')
    assert np.array_equal(a['step'], b['step'])
    dev = {'chi_mean': np.abs(b['chi_mean'] - (1 - a['chi_mean'])),
           'chi_target_mean': np.abs(b['chi_target_mean'] - (1 - a['chi_target_mean']))}
    for k in ('chi_std', 'm_mean', 'm_std', 'u_rms', 'P_std', 'source_mean'):
        dev[k] = np.abs(b[k] - a[k])
    early = a['step'] <= 500
    mirror = {k: {'max_first_500_steps': float(v[early].max()), 'max_all': float(v.max())}
              for k, v in dev.items()}
    assert all(v['max_first_500_steps'] < 1e-9 for v in mirror.values()), mirror
    fa, fb = frame(root/'legacy_step_new', 'frame1500.json'), frame(root/'hill_mirror', 'frame1500.json')
    mirror['frame1500_max_abs'] = {
        'chi': float(np.abs(fb['chi'] - (1 - fa['chi'])).max()),
        **{k: float(np.abs(fb[k] - fa[k]).max()) for k in ('m', 'QQxx', 'QQyx', 'pressure')}}
    out['hill_mirror'] = mirror

    # 3. Dm = Dchi = old Dbio, Dbio = 0: bit-identical to the Dbio run
    eq = same_bytes(root/'hill_mirror', root/'hill_mirror_split')
    assert all(eq.values()), eq
    out['hill_mirror_split'] = {'identical_files': sorted(eq)}

    # 4. local ODE against the Heun replica
    for name in ('ode_active_irrev', 'ode_active_rev', 'ode_passive', 'ode_steps'):
        r, q = response(root/name), CASES[name]
        ref = ode_replica(q)
        t = r['step'].astype(int)
        em = float(np.abs(r['m_mean'] - ref[t, 0]).max())
        ec = float(np.abs(r['chi_mean'] - ref[t, 1]).max())
        # chi_std/m_std in response.csv are sqrt(<x^2> - <x>^2): ~1e-7 of pure cancellation
        # on an exactly uniform field, so uniformity is checked on the frames instead.
        flat = float(max(r['u_rms'].max(), np.abs(r['P_mean']).max()))
        for fn in frames(root/name):
            f = frame(root/name, fn)
            flat = max(flat, float(np.ptp(f['chi'])), float(np.ptp(f['m'])))
        assert em < 1e-13 and ec < 1e-13 and flat < 1e-12, (name, em, ec, flat)
        out[name] = {'max_err_m': em, 'max_err_chi': ec, 'max_spatial_or_flow': flat,
                     'chi_start_peak_end': [float(ref[0, 1]), float(ref[:, 1].max()),
                                            float(ref[-1, 1])],
                     'm_end': float(ref[-1, 0])}
    # the physics the cases were built to show
    assert out['ode_active_irrev']['chi_start_peak_end'][2] > .85     # locked on
    assert out['ode_active_rev']['chi_start_peak_end'][2] < .05       # released
    assert out['ode_passive']['chi_start_peak_end'][1] < .2           # held passive

    # 5. split diffusion against the 1D replica
    r, q = response(root/'diff_split'), CASES['diff_split']
    ref, t = diffusion_replica(q), response(root/'diff_split')['step'].astype(int)
    errs = {}
    for fld in ('chi', 'm'):
        errs[fld+'_mean'] = float(np.abs(r[fld+'_mean'] - ref[fld][t].mean(1)).max())
        errs[fld+'_std'] = float(np.abs(r[fld+'_std'] - ref[fld][t].std(1)).max())
    assert max(errs.values()) < 1e-12, errs
    last = frame(root/'diff_split', frames(root/'diff_split')[-1])
    errs['frame_chi_max_abs'] = float(np.abs(last['chi'].reshape(64, 8)[:, 0]
                                             - ref['chi'][-1]).max())
    errs['frame_m_max_abs'] = float(np.abs(last['m'].reshape(64, 8)[:, 0] - ref['m'][-1]).max())
    out['diff_split'] = errs

    print(json.dumps(out, indent=1))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=('generate', 'check'))
    ap.add_argument('path', type=Path)
    ap.add_argument('--json', type=Path)
    a = ap.parse_args()
    if a.mode == 'generate':
        generate(a.path)
    else:
        res = check(a.path)
        if a.json:
            a.json.write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
