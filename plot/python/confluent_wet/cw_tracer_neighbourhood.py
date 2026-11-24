#!/usr/bin/env python3
"""How the pressure a cell feels depends on its own phenotype and on the phenotype around it, along the tracers.

A tracer run stores P, m and the position of N material paths (tracer.f32) and, every nvideo steps, the 8 x 8 block
means of chi (video_chi.u8). At each video time in the terminal window max(500, 20 tau_m) tau_c:
  own chi   chi_o: h(m_tracer) relaxed with tau_chi along the path (the model's chi at the tracer, Dchi dropped)
  around    chi_n: the block-mean chi field interpolated (bilinear between block centres) at the tracer,
            i.e. the mean phenotype within ~8 cells
  pressure  dP = P - <P> over the tracers at that time
Reductions:
  regression  dP on d(chi_o) alone, and on (d chi_o, d chi_n) together
  table       <dP>, the share of P > p_c and sigma_P in bins of (chi_o, chi_n)
  local       <dP>, sigma_P and the share of P > p_c in bins of chi_n alone
  lagged      slopes of dP(t + tau) on chi_o(t) alone, and on (chi_o(t), chi_n(t)) together, tau in [-10, 10] tau_c
An open-loop run (activity uniform, chi not fed back) is the control: there the dependence of P on chi is the
imprint P -> m -> chi alone.

  cw_tracer_neighbourhood.py <case_in> <case_out>      -> neighbourhood.json
  cw_tracer_neighbourhood.py <results_root> <out.json> --collect <out_subdir>
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import map_coordinates

from cw_mc_tm_scan import TC
from cw_tracer_drive import chi_path, h_of, load

EDGES = np.array([0., .1, .3, .7, .9, 1.000001])
LAGS = np.arange(-20, 21)                                # video frames, .5 tau_c each


def reduce_case(root, out):
    p, t, P, mtr, field = load(root)
    tm, tchi = float(p['tau_m']), float(p['tau_chi'])
    mc, pmem, hn = float(p['mc']), float(p['pmem']), float(p['h_n'])
    freeze, nsteps, nvideo = int(p['chi_freeze_steps']), int(p['nsteps']), int(p['nvideo'])
    stride, LX = int(p['video_stride']), int(p['LX'])
    nb = LX//stride
    dts = float(t[1] - t[0]); ns, N = P.shape
    chi_o_all = chi_path(h_of(mtr, mc, hn), float(p['chi0']), dts/tchi, int(np.ceil(freeze/dts)))

    video = np.fromfile(root/'video_chi.u8', np.uint8)
    nf = min(len(video)//(nb*nb), len(field))
    video = video[:nf*nb*nb].reshape(nf, nb, nb)
    tv = field['t'][:nf]
    obs = (nsteps - freeze)/TC
    window = max(500., 20*tm/TC)
    tvt = (tv - freeze)/TC
    idx = np.rint(tv/dts).astype(int)                    # tracer sample nearest to each video frame
    ok = idx < ns
    lo, hi = LAGS.min(), LAGS.max()
    frames = np.flatnonzero(ok & (tvt >= obs - window) & (tvt <= obs))
    frames = frames[(frames + lo >= 0) & (frames + hi < nf) & (idx[np.minimum(frames + hi, nf - 1)] < ns)]

    xy = np.memmap(root/'tracer.f32', np.float32, 'r').reshape(-1, N, 4)
    span = np.arange(frames[0] + lo, frames[-1] + hi + 1)                 # every frame any lag touches
    Pv = P[idx[span]]; dPv = Pv - Pv.mean(1, keepdims=True)
    co = chi_o_all[idx[span]]
    cn = np.empty_like(co)
    for k, f in enumerate(span):
        pos = np.asarray(xy[idx[f], :, 2:4], float)
        cn[k] = map_coordinates(video[f]/255., [pos[:, 0]/stride - .5, pos[:, 1]/stride - .5], order=1, mode='grid-wrap')
    base = frames - span[0]
    dP, O, Nn, Praw = dPv[base], co[base], cn[base], Pv[base]
    dO, dN = O - O.mean(1, keepdims=True), Nn - Nn.mean(1, keepdims=True)

    def fit2(y, a, b):
        """Least squares y = ka a + kb b (all centred per time); returns ka, kb, R2."""
        saa, sbb, sab, say, sby, syy = (float(np.mean(u*v)) for u, v in ((a, a), (b, b), (a, b), (a, y), (b, y), (y, y)))
        det = saa*sbb - sab*sab
        if det <= 1e-14*max(saa*sbb, 1e-300) or syy <= 0:
            return None
        ka, kb = (say*sbb - sby*sab)/det, (sby*saa - say*sab)/det
        return ka, kb, (ka*say + kb*sby)/syy

    def fit1(y, a):
        saa, say, syy = float(np.mean(a*a)), float(np.mean(a*y)), float(np.mean(y*y))
        return (say/saa, say*say/(saa*syy)) if saa > 1e-12 and syy > 0 else None

    f1, f2 = fit1(dP, dO), fit2(dP, dO, dN)
    nbins = len(EDGES) - 1
    bo, bn = np.digitize(O, EDGES) - 1, np.digitize(Nn, EDGES) - 1
    tab = {k: [[None]*nbins for _ in range(nbins)] for k in ('count', 'dP', 'share', 'sigma')}
    for i in range(nbins):
        for j in range(nbins):
            sel = (bo == i) & (bn == j)
            c = int(sel.sum())
            tab['count'][i][j] = c
            if c >= 200:
                tab['dP'][i][j] = float(dP[sel].mean()); tab['share'][i][j] = float((Praw[sel] > pmem).mean())
                tab['sigma'][i][j] = float(Praw[sel].std())
    local = []
    for j in range(nbins):
        sel = bn == j
        c = int(sel.sum())
        local.append(dict(count=c, dP=float(dP[sel].mean()), sigma=float(Praw[sel].std()),
                          share=float((Praw[sel] > pmem).mean())) if c >= 200 else dict(count=c))
    lag1, lag2 = [], []
    for L in LAGS:
        y = dPv[base + L]
        a, b = fit1(y, dO), fit2(y, dO, dN)
        lag1.append(a[0] if a else None); lag2.append(list(b[:2]) if b else None)

    row = dict(case=root.name, tm_over_tc=tm/TC, tchi_over_tc=tchi/TC, mc=mc, chi0=float(p['chi0']),
               open_loop=int(p.get('open_loop', 0)), n_frames=int(len(frames)), n_tracer=int(N), window_tc=window,
               field_chi=float(np.mean(field['chi_mean'][frames])), field_chi_std=float(np.mean(field['chi_std'][frames])),
               chi_own=dict(mean=float(O.mean()), std=float(O.std())), chi_around=dict(mean=float(Nn.mean()), std=float(Nn.std())),
               corr_own_around=float(np.corrcoef(dO.ravel(), dN.ravel())[0, 1]) if dO.std() > 0 and dN.std() > 0 else None,
               sigma_P=float(Praw.std()), share_above_pmem=float((Praw > pmem).mean()),
               fit_own=dict(slope=f1[0], r2=f1[1]) if f1 else None,
               fit_own_around=dict(own=f2[0], around=f2[1], r2=f2[2]) if f2 else None,
               edges=EDGES.tolist(), table=tab, by_around=local,
               lags_tc=(LAGS*nvideo/TC).tolist(), lag_own=lag1, lag_own_around=lag2)
    out.mkdir(parents=True, exist_ok=True)
    (out/'neighbourhood.json').write_text(json.dumps(row, indent=1, allow_nan=False) + '\n')
    print(json.dumps({k: row[k] for k in ('case', 'n_frames', 'corr_own_around', 'fit_own', 'fit_own_around')}))


def collect(results_root, out, sub):
    rows = [json.loads(f.read_text()) for f in sorted(results_root.glob(f'**/{sub}/neighbourhood.json'))]
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rows, indent=1) + '\n')
    print(f'{len(rows)} cases -> {out}')


if __name__ == '__main__':
    q = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    q.add_argument('src', type=Path)
    q.add_argument('out', type=Path)
    q.add_argument('--collect', metavar='OUT_SUBDIR')
    a = q.parse_args()
    if a.collect:
        collect(a.src, a.out, a.collect)
    else:
        reduce_case(a.src, a.out)
