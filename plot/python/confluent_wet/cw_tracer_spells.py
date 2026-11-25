#!/usr/bin/env python3
"""The temporal heterogeneity of the pressure along the tracers of one run, as a handful of numbers.

A cell (tracer) alternates between spells in active and in passive surroundings; in the adaptive nematic the
pressure it feels fluctuates more strongly where the surroundings are active, so its history is a process whose
strength is modulated by a slowly varying environment. Over the last 500 tau_c of the run, with
  P        the pressure at the tracer (every ntracer steps)
  env      chi around the tracer: the 8 x 8 block mean of the video stream interpolated at the tracer (bilinear in
           space, linear in time between frames)
  own      the model's chi at the tracer: h(m_tracer) relaxed with tau_chi
this writes
  by_env          sigma_P, <dP> and the share of P > p_c in fixed bins of env, and the same in the bottom and top
                  quartile of env (always defined, also where chi hardly varies)
  local variance  the variance of P in non-overlapping windows of 2 tau_c: its coefficient of variation and its
                  autocorrelation at 4, 6, 10, 20 tau_c and summed from 4 to 30 tau_c (0 for a process whose
                  strength is not modulated on these times)
  persistence     1/e times of P, of env and of own; the integral time of s(P); mean residence of env above and
                  below .5
  one point       sigma, skewness of P, <s>, the memory width without diffusion at the run's tau_m

  cw_tracer_spells.py <case_in> <case_out>             -> spells.json
  cw_tracer_spells.py <results_root> <out.json> --collect <out_subdir>
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import map_coordinates

from cw_mc_tm_scan import TC
from cw_tracer_drive import chi_path, h_of, load, memory_path, s_of

WINDOW = 500.                 # tau_c
W_LOC = 2.                    # tau_c
EDGES = np.array([0., .2, .4, .6, .8, 1.0001])
LAGMAX = 30.                  # tau_c


def acf(x, nlag):
    d = x - x.mean(1, keepdims=True)
    v = np.mean(d*d)
    if v <= 0:
        return np.zeros(nlag)
    return np.array([np.mean(d[:len(d) - k]*d[k:]) for k in range(nlag)])/v


def e_time(a, dt):
    k = np.flatnonzero(a < np.exp(-1))
    return float(k[0]*dt) if len(k) else None


def reduce_case(root, out):
    p, t, P, mtr, field = load(root)
    tm, tchi, mc, pmem, wp, hn = (float(p[k]) for k in ('tau_m', 'tau_chi', 'mc', 'pmem', 'wp', 'h_n'))
    freeze, nsteps, stride, LX = int(p['chi_freeze_steps']), int(p['nsteps']), int(p['video_stride']), int(p['LX'])
    dts = float(t[1] - t[0]); dt = dts/TC; N = P.shape[1]
    mpath = memory_path(P, mtr[0], dts/tm, pmem, wp)
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
    Pw, own, mpath = P[sel], own[sel], mpath[sel]
    dev = Pw - Pw.mean(1, keepdims=True)
    s = s_of(Pw, pmem, wp)

    def cell(mask):
        if mask.sum() < 500:
            return None
        return dict(frac=float(mask.mean()), dP=float(dev[mask].mean()), sigma=float(Pw[mask].std()), share=float((Pw[mask] > pmem).mean()),
                    own=float(own[mask].mean()), env=float(env[mask].mean()))
    b = np.digitize(env, EDGES) - 1
    q1, q3 = np.quantile(env, [.25, .75])
    lo, hi = cell(env <= q1), cell(env >= q3)

    nlag = int(LAGMAX/dt)
    aP, aS, aE, aO = acf(Pw, nlag), acf(s, nlag), acf(env, nlag), acf(own, nlag)
    nw = int(round(W_LOC/dt)); k = len(Pw)//nw
    lv = Pw[:k*nw].reshape(k, nw, N).var(1)
    av = acf(lv, int(LAGMAX/W_LOC) + 1)
    res = {}
    for above, tag in ((True, 'active'), (False, 'passive')):
        dur = []
        for c in range(N):
            z = (env[:, c] > .5) == above
            e = np.flatnonzero(np.diff(z.astype(np.int8)))
            if len(e) > 1:
                d = np.diff(e)[z[e[:-1] + 1]]
                dur.append(d)
        dur = np.concatenate(dur)*dt if dur else np.array([])
        res[tag] = dict(n=int(len(dur)), mean=float(dur.mean()) if len(dur) > 50 else None,
                        share_over_20=float((dur > 20).mean()) if len(dur) > 50 else None, time_share=float(((env > .5) == above).mean()))
    ftt = (field['t'] - freeze)/TC; ftail = (ftt >= obs - WINDOW) & (ftt <= obs)
    lag = lambda x: int(round(x/W_LOC))
    row = dict(
        case=root.name, tm_over_tc=tm/TC, tchi_over_tc=tchi/TC, mc=mc, chi0=float(p['chi0']), open_loop=int(p.get('open_loop', 0)),
        window_tc=WINDOW, sample_tc=dt, field_chi=float(field['chi_mean'][ftail].mean()), field_chi_std=float(field['chi_std'][ftail].mean()),
        sigma_P=float(Pw.std()), skew_P=float(np.mean(dev**3)/dev.std()**3), s_mean=float(s.mean()), m_path_std=float(mpath.std()),
        env_mean=float(env.mean()), env_std=float(env.std()), corr_own_env=float(np.corrcoef(own.ravel(), env.ravel())[0, 1]) if env.std() > 0 and own.std() > 0 else None,
        env_edges=EDGES.tolist(), by_env=[cell(b == j) for j in range(len(EDGES) - 1)], env_low_quartile=lo, env_high_quartile=hi,
        sigma_ratio_quartiles=(hi['sigma']/lo['sigma']) if lo and hi else None,
        local_var_cv=float(lv.std()/lv.mean()), local_sigma_cv=float(np.sqrt(lv).std()/np.sqrt(lv).mean()),
        local_var_acf={f'{x:g}': float(av[lag(x)]) for x in (2, 4, 6, 10, 20)},
        local_var_persistence_tc=float(av[lag(4):].sum()*W_LOC),
        tau_P_1e=e_time(aP, dt), tau_env_1e=e_time(aE, dt), tau_own_1e=e_time(aO, dt), tau_s_int=float(np.sum(aS)*dt), residence=res)
    out.mkdir(parents=True, exist_ok=True)
    (out/'spells.json').write_text(json.dumps(row, indent=1, allow_nan=False) + '\n')
    print(json.dumps({k: row[k] for k in ('case', 'sigma_ratio_quartiles', 'local_var_acf', 'local_var_persistence_tc', 'tau_env_1e', 'tau_s_int')}))


def collect(results_root, out, sub):
    rows = [json.loads(f.read_text()) for f in sorted(results_root.glob(f'**/{sub}/spells.json'))]
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
