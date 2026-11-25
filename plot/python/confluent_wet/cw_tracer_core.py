#!/usr/bin/env python3
"""Everything the core heterogeneity figure needs from one tracer run, as a small file.

Over the last 500 tau_c of the run (P at the tracer every ntracer steps; env = the 8 x 8 block mean of chi at the
tracer, bilinear in space and linear in time between video frames; own = h(m_tracer) relaxed with tau_chi):
  cell     one cell, the first of a fixed random permutation of the tracer grid (seed 1; the same grid index in
           every run), for the first 150 tau_c of the window: P, env, own, m and the standard deviation of P in a
           sliding window of 5 tau_c centred on t
  sig_env  sigma_P in ten bins of env (NaN where a bin holds fewer than 2000 samples), and the time share of each bin
  hist     histograms of P while env < .3 and while env > .7 (90 bins on [-0.09, 0.09]), with n, mean, sigma
  acf      Lagrangian autocorrelation of P, m (the model's memory at the tracer) and own, the instantaneous mean
           over the tracers removed, every 3 samples up to 40 tau_c and every 15 up to 150 tau_c

  cw_tracer_core.py <case_in> <case_out>      -> core.npz
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import map_coordinates

from cw_mc_tm_scan import TC
from cw_tracer_drive import chi_path, h_of, load

WINDOW, SHOW, W_SLIDE, LAGMAX = 500., 150., 5., 150.     # tau_c
EDGES = np.linspace(0, 1, 11)
PBINS = np.linspace(-.09, .09, 91)
LOW, HIGH = .3, .7


def acf(x, lags):
    d = x - x.mean(1, keepdims=True)
    v = np.mean(d*d)
    if v <= 0:
        return np.zeros(len(lags))
    return np.array([np.mean(d[:len(d) - k]*d[k:]) for k in lags])/v


def reduce_case(root, out):
    p, t, P, mtr, field = load(root)
    tm, tchi, mc, pmem, hn = (float(p[k]) for k in ('tau_m', 'tau_chi', 'mc', 'pmem', 'h_n'))
    freeze, nsteps, stride, LX = int(p['chi_freeze_steps']), int(p['nsteps']), int(p['video_stride']), int(p['LX'])
    dts = float(t[1] - t[0]); dt = dts/TC; N = P.shape[1]
    own = chi_path(h_of(mtr, mc, hn), float(p['chi0']), dts/tchi, int(np.ceil(freeze/dts)))
    tt = (t - freeze)/TC; obs = (nsteps - freeze)/TC
    sel = np.flatnonzero((tt >= obs - WINDOW) & (tt <= obs))
    nb = LX//stride
    video = np.fromfile(root/'video_chi.u8', np.uint8)
    nf = min(len(video)//(nb*nb), len(field)); video = video[:nf*nb*nb].reshape(nf, nb, nb)
    tv = field['t'][:nf]; dv = float(tv[1] - tv[0])
    xy = np.memmap(root/'tracer.f32', np.float32, 'r').reshape(-1, N, 4)
    fr = np.flatnonzero((tv >= t[sel[0]] - 2*dv) & (tv <= t[sel[-1]] + 2*dv) & (tv <= t[-1]))
    env_f = np.empty((len(fr), N))
    for k, f in enumerate(fr):
        pos = np.asarray(xy[int(round(tv[f]/dts)), :, 2:4], float)
        env_f[k] = map_coordinates(video[f]/255., [pos[:, 0]/stride - .5, pos[:, 1]/stride - .5], order=1, mode='grid-wrap')
    u = (t[sel] - tv[fr[0]])/dv; i0 = np.clip(np.floor(u).astype(int), 0, len(fr) - 2); w = (u - i0)[:, None]
    env = (1 - w)*env_f[i0] + w*env_f[i0 + 1]
    Pw, own, m = P[sel], own[sel], mtr[sel]

    pick = int(np.random.default_rng(1).permutation(N)[0]); ns = int(SHOW/dt)
    x = Pw[:, pick]; h = int(round(.5*W_SLIDE/dt))
    c1, c2 = np.concatenate([[0.], np.cumsum(x)]), np.concatenate([[0.], np.cumsum(x*x)])
    lo, hi = np.maximum(np.arange(len(x)) - h, 0), np.minimum(np.arange(len(x)) + h + 1, len(x))
    slide = np.sqrt(np.maximum((c2[hi] - c2[lo])/(hi - lo) - ((c1[hi] - c1[lo])/(hi - lo))**2, 0.))

    b = np.digitize(env, EDGES) - 1
    sig = np.array([Pw[b == j].std() if (b == j).sum() > 2000 else np.nan for j in range(10)])
    frac = np.array([(b == j).mean() for j in range(10)])
    stats = {}
    hist = {}
    for tag, mask in (('low', env < LOW), ('high', env > HIGH)):
        v = Pw[mask]
        hist[tag] = np.histogram(v, PBINS)[0]
        stats[tag] = dict(n=int(v.size), mean=float(v.mean()) if v.size else None, sigma=float(v.std()) if v.size else None,
                          share_above_pc=float((v > pmem).mean()) if v.size else None)
    fine = int(40/dt)                                   # every 3 samples to 40 tau_c, every 15 beyond (slow chi at large tau_chi)
    lags = np.concatenate([np.arange(0, fine, 3), np.arange(fine, int(LAGMAX/dt), 15)])
    ftt = (field['t'] - freeze)/TC; ftail = (ftt >= obs - WINDOW) & (ftt <= obs)
    open_loop = int(p.get('open_loop', 0)); z0 = float(p['zeta0_frac'])
    chi_field = float(field['chi_mean'][ftail].mean())
    meta = dict(case=root.name, open_loop=open_loop, tm_over_tc=tm/TC, tchi_over_tc=tchi/TC, mc=mc, pmem=pmem, chi0=float(p['chi0']), dt=dt,
                activity=float(p['zeta_open'])/float(p['zeta']) if open_loop else z0 + (1 - z0)*chi_field, field_chi=chi_field,
                sigma_all=float(Pw.std()), cell=pick, window_tc=WINDOW, slide_tc=W_SLIDE, low=LOW, high=HIGH, by_env=stats)
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out/'core.npz', meta=json.dumps(meta), t=(tt[sel] - tt[sel[0]])[:ns], P_cell=x[:ns], env_cell=env[:ns, pick], own_cell=own[:ns, pick],
        m_cell=m[:ns, pick], slide_cell=slide[:ns], edges=EDGES, sig_env=sig, frac_env=frac, pbins=PBINS, hist_low=hist['low'], hist_high=hist['high'],
        lag=lags*dt, acf_P=acf(Pw, lags), acf_m=acf(m, lags), acf_own=acf(own, lags))
    print(json.dumps(dict(case=root.name, sigma_all=meta['sigma_all'], sig_env=[None if np.isnan(v) else round(float(v), 5) for v in sig])))


if __name__ == '__main__':
    q = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    q.add_argument('src', type=Path)
    q.add_argument('out', type=Path)
    a = q.parse_args()
    reduce_case(a.src, a.out)
