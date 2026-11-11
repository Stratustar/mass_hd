#!/usr/bin/env python3
"""Per-run pressure and phenotype windows around +-1/2 defects: mean maps and the pixelwise P-chi correlation.

  cw_defect_corr.py <run_dir> <out_dir> [--min-tc 1000] [--only hn16_wp0p002]

Frames: every full frame*.json at least MIN_TC reference tau_c after release (chi-freeze-steps).
Per frame, per defect: total pressure P and chi sampled bilinearly (periodic) on a 25 x 25 window (step 1 grid,
half width 12) in the defect frame, map row = along the axis (+y), map column = across.
  plus    +1/2 cores (plaquette winding); axis = head h = -div Q/|div Q| (as cw_defect_polar.py).
  minus   -1/2 cores; axis = one of the three arms along which the director is radial,
          psi = arg[(dx Qxx - dy Qxy) + i (dx Qxy + dy Qxx)]/3 + 2 pi k/3. Each core enters three times (k = 0, 1, 2),
          so the maps are threefold symmetric by construction.
  random  as many random points as +1/2 cores, random axis (baseline: what any point of the field gives).
Each window is referred to a zero in three ways (variants):
  frame   minus the whole-frame mean
  win     minus the mean of that window over the disk r <= 12        (the per-defect normalisation)
  z       'win' divided by that window's standard deviation over the disk (windows with zero spread are dropped)
Sums over windows, per frame: S[frame, kind, variant, q, row, col], q = (P, chi, P^2, chi^2, P chi), and counts N.
Also SQ[frame, kind, (xx, xy), row, col]: the sum of Q turned into the window frame, Q'xx + i Q'xy = (Qxx + i Qxy) exp(-2 i alpha),
alpha = the lab angle of the across direction; its average gives the mean director around the defect (angle from the across axis).
From these, per pixel: mean maps <P>, <chi>, and the Pearson correlation over defects
    r(x, y) = <dP dchi> / (sigma_P sigma_chi),  dP_n = P_n(x, y) - <P>(x, y)
Output <out_dir>/corr_case.npz. `accumulate` and `maps` are imported by analysis/defect_pchi_corr.py.
"""
import argparse
import json
from pathlib import Path
import re
import sys

import numpy as np
from scipy.ndimage import map_coordinates

TC = 674.3290333006435
H = 12
U = np.arange(-H, H + 1.)
UU, VV = np.meshgrid(U, U, indexing='xy')          # UU across the axis, VV along it
DISK = np.hypot(UU, VV) <= H
KINDS = ('plus', 'minus', 'random')
VARIANTS = ('frame', 'win', 'z')
NQ = 5


def grad(f):
    return (np.roll(f, -1, 0) - np.roll(f, 1, 0))/2, (np.roll(f, -1, 1) - np.roll(f, 1, 1))/2


def at(f, xy):
    xy = np.asarray(xy, float)
    return map_coordinates(f, xy.reshape(-1, 2).T, order=1, mode='grid-wrap').reshape(xy.shape[:-1])


def cores(qxx, qxy):
    th = .5*np.arctan2(qxy, qxx)
    d = lambda a, b: (b - a + np.pi/2) % np.pi - np.pi/2
    t10 = np.roll(th, -1, 0); t11 = np.roll(t10, -1, 1); t01 = np.roll(th, -1, 1)
    w = (d(th, t10) + d(t10, t11) + d(t11, t01) + d(t01, th))/(2*np.pi)
    return np.argwhere(w > .25) + .5, np.argwhere(w < -.25) + .5


def heads(qxx, qxy, plus):
    ga, gb = grad(qxx); gc, ge = grad(qxy)
    px, py = at(ga + ge, plus), at(gc - gb, plus)
    return -np.c_[px, py]/np.maximum(np.hypot(px, py), 1e-300)[:, None]


def arms(qxx, qxy, minus):
    """the three arm directions of each -1/2, shape (3, n, 2)."""
    ga, gb = grad(qxx); gc, ge = grad(qxy)
    psi = np.arctan2(at(gc + gb, minus), at(ga - ge, minus))/3
    return np.stack([np.c_[np.cos(psi + k*2*np.pi/3), np.sin(psi + k*2*np.pi/3)] for k in range(3)])


def windows(f, pts, h):
    e = np.c_[h[:, 1], -h[:, 0]]
    return at(f, pts[:, None, None, :] + UU[None, :, :, None]*e[:, None, None, :] + VV[None, :, :, None]*h[:, None, None, :])


def qsum(qxx, qxy, pts, h):
    """sum over windows of Q in the window frame (across, along), shape (2, 25, 25)."""
    if not len(pts):
        return np.zeros((2,) + UU.shape)
    q = (windows(qxx, pts, h) + 1j*windows(qxy, pts, h))*np.exp(-2j*np.arctan2(-h[:, 0], h[:, 1]))[:, None, None]
    return np.stack([q.real.sum(0), q.imag.sum(0)])


def director(SQ, N, step=2):
    """mean director from summed Q: positions (u, v), unit components (du, dv) and the order |<Q>| on every step-th pixel."""
    q = (SQ[0] + 1j*SQ[1])/max(float(N), 1.)
    sl = slice(H % step, None, step)
    th = .5*np.angle(q[sl, sl])
    return UU[sl, sl], VV[sl, sl], np.cos(th), np.sin(th), np.abs(q[sl, sl])


def sums(Pw, Cw, Pm, Cm):
    """Pw, Cw: (n, 25, 25) windows; Pm, Cm: frame means. Returns S (3 variants, 5, 25, 25), N (3,)."""
    S = np.zeros((len(VARIANTS), NQ) + UU.shape)
    N = np.zeros(len(VARIANTS))
    if not len(Pw):
        return S, N
    pw = Pw - Pw[:, DISK].mean(1)[:, None, None]
    cw = Cw - Cw[:, DISK].mean(1)[:, None, None]
    sp, sc = pw[:, DISK].std(1), cw[:, DISK].std(1)
    ok = (sp > 1e-14) & (sc > 1e-14)
    for v, (p, c) in enumerate(((Pw - Pm, Cw - Cm), (pw, cw), (pw[ok]/sp[ok, None, None], cw[ok]/sc[ok, None, None]))):
        S[v] = p.sum(0), c.sum(0), (p*p).sum(0), (c*c).sum(0), (p*c).sum(0)
        N[v] = len(p)
    return S, N


def accumulate(qxx, qxy, chi, P, rng):
    """one frame -> S (3 kinds, 3 variants, 5, 25, 25), N (3 kinds, 3 variants), (n_plus, n_minus), SQ (3 kinds, 2, 25, 25)."""
    L = chi.shape[0]
    plus, minus = cores(qxx, qxy)
    Pm, Cm = P.mean(), chi.mean()
    S = np.zeros((len(KINDS), len(VARIANTS), NQ) + UU.shape)
    N = np.zeros((len(KINDS), len(VARIANTS)))
    SQ = np.zeros((len(KINDS), 2) + UU.shape)
    h = heads(qxx, qxy, plus)
    S[0], N[0] = sums(windows(P, plus, h), windows(chi, plus, h), Pm, Cm)
    SQ[0] = qsum(qxx, qxy, plus, h)
    for a in arms(qxx, qxy, minus):
        s, n = sums(windows(P, minus, a), windows(chi, minus, a), Pm, Cm)
        S[1] += s; N[1] += n
        SQ[1] += qsum(qxx, qxy, minus, a)
    pts = rng.uniform(0, L, (len(plus), 2))
    ang = rng.uniform(0, 2*np.pi, len(plus))
    a = np.c_[np.cos(ang), np.sin(ang)]
    S[2], N[2] = sums(windows(P, pts, a), windows(chi, pts, a), Pm, Cm)
    SQ[2] = qsum(qxx, qxy, pts, a)
    return S, N, (len(plus), len(minus)), SQ


def maps(S, N):
    """S (..., 5, 25, 25), N (...) summed over frames -> mean P, mean chi, r, sigma_P, sigma_chi (each (..., 25, 25))."""
    n = np.asarray(N, float)[..., None, None]
    mp, mc = S[..., 0, :, :]/n, S[..., 1, :, :]/n
    vp, vc = S[..., 2, :, :]/n - mp**2, S[..., 3, :, :]/n - mc**2
    with np.errstate(invalid='ignore', divide='ignore'):
        r = (S[..., 4, :, :]/n - mp*mc)/np.sqrt(vp*vc)
    return mp, mc, r, np.sqrt(np.maximum(vp, 0)), np.sqrt(np.maximum(vc, 0))


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
    rng = np.random.default_rng(int(p['seed']))
    S, N, cnt, SQ, chi_mean, P_mean = [], [], [], [], [], []
    for s in steps:
        d = json.loads((a.run/f'frame{s}.json').read_text())['data']
        F = {k: np.asarray(d[k]['value'], float).reshape(L, L) for k in ('QQxx', 'QQyx', 'chi', 'pressure')}
        s_, n_, c_, q_ = accumulate(F['QQxx'], F['QQyx'], F['chi'], F['pressure'], rng)
        S.append(s_); N.append(n_); cnt.append(c_); SQ.append(q_)
        chi_mean.append(float(F['chi'].mean())); P_mean.append(float(F['pressure'].mean()))
    a.out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.out/'corr_case.npz', S=np.array(S), N=np.array(N), SQ=np.array(SQ), u=U, kinds=np.array(KINDS), variants=np.array(VARIANTS),
                        quantities=np.array(['P', 'chi', 'P2', 'chi2', 'Pchi']), steps=np.array(steps), n_plus_minus=np.array(cnt),
                        chi_mean=np.array(chi_mean), P_mean=np.array(P_mean), tau_m_over_tc=float(p['tau_m'])/TC,
                        tau_chi_over_tc=float(p['tau_chi'])/TC, chi0=float(p['chi0']), seed=int(p['seed']), case=a.run.name)
    print(json.dumps({'case': a.run.name, 'frames': len(steps), 'n_plus_minus': np.array(cnt).tolist(), 'chi_mean': chi_mean}))


if __name__ == '__main__':
    main()
