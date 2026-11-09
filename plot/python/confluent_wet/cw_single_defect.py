#!/usr/bin/env python3
"""Pressure around a single +1/2 (defect-pair start, uniform open-loop activity), several definitions.

  cw_single_defect.py <run_dir> <out.npz> [--t-min 2] [--min-sep 32]

Uses full frames (frame-light = 0). Per frame: +1/2 / -1/2 cores by plaquette winding; a frame is kept when
t >= T_MIN tau_c (flow spun up), there is exactly one +1/2 and the nearest -1/2 is >= MIN_SEP away
(periodic). Fields, each minus its frame mean, sampled bilinearly on a Cartesian head frame (x along e,
y along the head h = -div Q/|div Q|, |x|, |y| <= 12, step 0.5):
  P        = p_LB - sigma_bulk         the code's mechanical pressure (field 'pressure')
  p_LB     = (n - rho)/3               ('pressure_lb')
  sigma_bulk                           the non-LB isotropic stress ('sigma_bulk'); every other stress is traceless
  p_inc    Poisson pressure of the non-LB stress, lap p = d_i d_j sigma_ij (FFT, periodic): the pressure an
           incompressible fluid would need to balance div sigma (sigma = sigmaXX, XY, YX, YY of the frame)
Q is rotated into the head frame and averaged (director = half the angle of the mean). Output: maps,
the kept frame times, the +1/2 track and its head-projected speed.
"""
import argparse
import json
from pathlib import Path
import re

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.spatial import cKDTree

TC = 674.3290333006435
X = np.arange(-12, 12.01, .5)


def at(f, xy):
    s = xy.shape[:-1]
    return map_coordinates(f, xy.reshape(-1, 2).T, order=1, mode='wrap').reshape(s)


def grad(f):
    return (np.roll(f, -1, 0) - np.roll(f, 1, 0))/2, (np.roll(f, -1, 1) - np.roll(f, 1, 1))/2


def cores(qxx, qxy):
    th = .5*np.arctan2(qxy, qxx)
    d = lambda a, b: (b - a + np.pi/2) % np.pi - np.pi/2
    t10 = np.roll(th, -1, 0); t11 = np.roll(t10, -1, 1); t01 = np.roll(th, -1, 1)
    w = (d(th, t10) + d(t10, t11) + d(t11, t01) + d(t01, th))/(2*np.pi)
    return np.argwhere(w > .25) + .5, np.argwhere(w < -.25) + .5


def poisson(sxx, sxy, syx, syy):
    L = sxx.shape[0]; k = 2*np.pi*np.fft.fftfreq(L)
    kx, ky = np.meshgrid(k, k, indexing='ij'); k2 = kx**2 + ky**2; k2[0, 0] = 1
    rhs = kx*kx*np.fft.fft2(sxx) + kx*ky*np.fft.fft2(sxy + syx) + ky*ky*np.fft.fft2(syy)   # -(d_i d_j sigma_ij)^
    p = rhs/k2; p[0, 0] = 0                                                                   # -k^2 p = -(k k : sigma)
    return np.fft.ifft2(p).real


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run', type=Path); ap.add_argument('out', type=Path)
    ap.add_argument('--t-min', type=float, default=2.); ap.add_argument('--min-sep', type=float, default=32.)
    a = ap.parse_args()
    par = {k: v['value'] for k, v in json.loads((a.run/'parameters.json').read_text())['data'].items()}
    L = int(par['LX'])
    steps = sorted(int(re.search(r'frame(\d+)\.json', f.name).group(1)) for f in a.run.glob('frame*.json'))
    XX, YY = np.meshgrid(X, X, indexing='xy')
    S = {k: 0. for k in ('P', 'p_LB', 'sigma_bulk', 'p_inc', 'Qxx', 'Qxy')}
    kept, track, counts = [], [], []
    for s in steps:
        t = s/TC
        d = json.loads((a.run/f'frame{s}.json').read_text())['data']
        F = {k: np.asarray(d[k]['value'], float).reshape(L, L) for k in
             ('QQxx', 'QQyx', 'pressure', 'pressure_lb', 'sigma_bulk', 'sigmaXX', 'sigmaXY', 'sigmaYX', 'sigmaYY')}
        pl, mi = cores(F['QQxx'], F['QQyx']); counts.append((t, len(pl), len(mi)))
        if len(pl) != 1:
            continue
        ga, gb = grad(F['QQxx']); gc, ge = grad(F['QQyx'])
        px, py = at(ga + ge, pl), at(gc - gb, pl)
        h = -np.array([px[0], py[0]])/max(np.hypot(px[0], py[0]), 1e-300); e = np.array([h[1], -h[0]])
        track.append((t, *pl[0], *h))
        sep = cKDTree(mi, boxsize=L).query(pl)[0][0] if len(mi) else np.inf
        if t < a.t_min or sep < a.min_sep:
            continue
        xy = pl[0][None, None, :] + XX[..., None]*e + YY[..., None]*h
        fields = {'P': F['pressure'], 'p_LB': F['pressure_lb'], 'sigma_bulk': F['sigma_bulk'],
                  'p_inc': poisson(F['sigmaXX'], F['sigmaXY'], F['sigmaYX'], F['sigmaYY'])}
        for k, f in fields.items():
            S[k] = S[k] + at(f - f.mean(), xy)
        psi = np.arctan2(e[1], e[0]); qa, qb = at(F['QQxx'], xy), at(F['QQyx'], xy)
        S['Qxx'] = S['Qxx'] + qa*np.cos(2*psi) + qb*np.sin(2*psi)
        S['Qxy'] = S['Qxy'] + -qa*np.sin(2*psi) + qb*np.cos(2*psi)
        kept.append(t)
    n = max(len(kept), 1)
    tr = np.array(track) if track else np.zeros((0, 5))
    v = []
    for i in range(1, len(tr)):
        dpos = (tr[i, 1:3] - tr[i - 1, 1:3] + L/2) % L - L/2
        v.append((tr[i, 0], float(dpos @ tr[i - 1, 3:5])/(tr[i, 0] - tr[i - 1, 0])))
    a.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.out, x=X, kept_t=np.array(kept), counts=np.array(counts), track=tr, speed=np.array(v),
                        zeta_open=float(par.get('zeta_open', np.nan)), **{k: S[k]/n for k in S})
    print(json.dumps({'run': a.run.name, 'frames': len(steps), 'kept': len(kept),
                      'kept_t': [round(kept[0], 2), round(kept[-1], 2)] if kept else None,
                      'counts_last': counts[-1] if counts else None}))


if __name__ == '__main__':
    main()
