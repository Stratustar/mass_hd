#!/usr/bin/env python3
"""Per-run +1/2 defect-centred polar sums of chi and P (steady frames), for scan-wide maps.

  cw_defect_polar.py <run_dir> <out_dir> [--min-tc 1000] [--only hn16_wp0p002]

Frames: every full frame*.json at least MIN_TC reference tau_c after release (chi-freeze-steps).
Per frame: +1/2 cores by plaquette winding, head h = -div Q/|div Q|; raw chi and P sampled
bilinearly (periodic) at core + r (cos(theta) e + sin(theta) h), e = h turned by -90 deg, so the
head sits at theta = 90 deg; r = 0.5..12 (0.5), theta in 5 deg bins. Cores are also classed by the
Gaussian(2)-smoothed chi at the core: passive < .2 <= intermediate <= .8 < active.

Output <out_dir>/polar_case.npz: sums S_chi, S_P, S_chi2, S_P2, S_chiP with shape (4, nr, nth) for
classes (all, passive, intermediate, active), counts N (4,), frame list, frame means of chi and P,
+1/2 count per frame, and the run's tau_m, tau_chi, chi0, seed. Runs whose name does not contain
--only are skipped with exit 0 (the analysis array runs over a whole campaign list).
"""
import argparse
import json
from pathlib import Path
import re
import sys

import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates

TC = 674.3290333006435
R = np.arange(0.5, 12.0 + 1e-9, 0.5)
DTH = 5
TH = np.radians(np.arange(DTH/2, 360, DTH))
CLASS_EDGES = (.2, .8)


def grad(f):
    return (np.roll(f, -1, 0) - np.roll(f, 1, 0))/2, (np.roll(f, -1, 1) - np.roll(f, 1, 1))/2


def plus_defects(qxx, qxy):
    th = .5*np.arctan2(qxy, qxx)
    d = lambda a, b: (b - a + np.pi/2) % np.pi - np.pi/2
    t10 = np.roll(th, -1, 0); t11 = np.roll(t10, -1, 1); t01 = np.roll(th, -1, 1)
    w = (d(th, t10) + d(t10, t11) + d(t11, t01) + d(t01, th))/(2*np.pi)
    return np.argwhere(w > .25) + .5


def at(f, xy):
    s = xy.shape[:-1]
    return map_coordinates(f, xy.reshape(-1, 2).T, order=1, mode='wrap').reshape(s)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('run', type=Path)
    ap.add_argument('out', type=Path)
    ap.add_argument('--min-tc', type=float, default=1000.)
    ap.add_argument('--only', default='hn16_wp0p002')
    a = ap.parse_args()
    if a.only not in a.run.name:
        print(f'skip {a.run.name}')
        return
    p = {k: v['value'] for k, v in json.loads((a.run/'parameters.json').read_text())['data'].items()}
    L, prep = int(p['LX']), int(p['chi_freeze_steps'])
    steps = sorted(int(re.search(r'frame(\d+)\.json', f.name).group(1)) for f in a.run.glob('frame*.json'))
    steps = [s for s in steps if (s - prep)/TC >= a.min_tc]
    if not steps:
        sys.exit(f'no frames >= {a.min_tc} tau_c after release in {a.run}')
    e_r = np.stack([np.cos(TH), np.sin(TH)], -1)
    S = {k: np.zeros((4, len(R), len(TH))) for k in ('chi', 'P', 'chi2', 'P2', 'chiP')}
    N = np.zeros(4, int)
    chi_means, P_means, nplus = [], [], []
    for s in steps:
        d = json.loads((a.run/f'frame{s}.json').read_text())['data']
        F = {k: np.asarray(d[k]['value'], float).reshape(L, L) for k in ('QQxx', 'QQyx', 'chi', 'pressure')}
        plus = plus_defects(F['QQxx'], F['QQyx'])
        ga, gb = grad(F['QQxx']); gc, ge = grad(F['QQyx'])
        px, py = at(ga + ge, plus), at(gc - gb, plus)
        h = -np.c_[px, py]/np.maximum(np.hypot(px, py), 1e-300)[:, None]
        e = np.c_[h[:, 1], -h[:, 0]]
        xy = plus[:, None, None, :] + R[None, :, None, None]*(e_r[None, None, :, 0, None]*e[:, None, None, :]
                                                              + e_r[None, None, :, 1, None]*h[:, None, None, :])
        c, P = at(F['chi'], xy), at(F['pressure'], xy)
        cs = at(gaussian_filter(F['chi'], 2, mode='wrap'), plus)
        cls = np.where(cs < CLASS_EDGES[0], 1, np.where(cs > CLASS_EDGES[1], 3, 2))
        for k in range(4):
            sel = np.ones(len(plus), bool) if k == 0 else cls == k
            S['chi'][k] += c[sel].sum(0); S['P'][k] += P[sel].sum(0)
            S['chi2'][k] += (c[sel]**2).sum(0); S['P2'][k] += (P[sel]**2).sum(0)
            S['chiP'][k] += (c[sel]*P[sel]).sum(0); N[k] += int(sel.sum())
        chi_means.append(float(F['chi'].mean())); P_means.append(float(F['pressure'].mean())); nplus.append(len(plus))
    a.out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.out/'polar_case.npz', r=R, theta_deg=np.degrees(TH), classes=np.array(['all', 'passive', 'intermediate', 'active']),
                        N=N, steps=np.array(steps), chi_mean=np.array(chi_means), P_mean=np.array(P_means), n_plus=np.array(nplus),
                        tau_m_over_tc=float(p['tau_m'])/TC, tau_chi_over_tc=float(p['tau_chi'])/TC, chi0=float(p['chi0']),
                        seed=int(p['seed']), case=a.run.name, **{f'S_{k}': v for k, v in S.items()})
    print(json.dumps({'case': a.run.name, 'frames': len(steps), 'N': N.tolist(), 'chi_mean': chi_means}))


if __name__ == '__main__':
    main()
