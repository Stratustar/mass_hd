#!/usr/bin/env python3
"""The open-loop map of the single cell driven by a tracer run's pressure, for a whole grid of tau_m at once.

In an OPEN-LOOP run the activity is uniform and chi is not fed back, so the pressure along the tracers depends on
the pinned activity alone -- not on the run's own tau_m or tau_chi. Each such run is therefore one entry of a
library of real pressure histories indexed by the activity a, and the single cell (no diffusion closure,
m' = (s(P) - m)/tau_m, <chi> = <h(m)>) can be driven by it at ANY tau_m and mc:

    F(a; tau_m, mc) = < h(m_path) >   over the tracers and the last `window` tau_c of the run

which is the open-loop map whose fixed points x = F(z0 + (1 - z0) x) are the self-consistent states of the single
cell with the real pressure of uniform turbulence. Every memory starts from the tracer's own m at t = 0 and has
the first half of the run to relax (>= 25 tau_m at the largest tau_m).

  cw_tracer_map.py <case_in> <case_out>                 -> tracer_map.json
  cw_tracer_map.py <results_root> <out.json> --collect <out_subdir>
"""
import argparse
import json
from pathlib import Path

import numpy as np

from cw_mc_tm_scan import TC
from cw_tracer_drive import SUB, h_of, load, s_of

TM = (0.3, 1, 2, 3.5, 5, 6, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 13, 14.5, 16, 18)   # tau_c
MCS = (0.2287, 0.22)
WINDOW = 500.                                           # tau_c


def reduce_case(root, out):
    p, t, P, mtr, field = load(root)
    pmem, wp, hn = float(p['pmem']), float(p['wp']), float(p['h_n'])
    freeze, nsteps = int(p['chi_freeze_steps']), int(p['nsteps'])
    dts = float(t[1] - t[0])
    tt = (t - freeze)/TC; obs = (nsteps - freeze)/TC
    tail = (tt >= obs - WINDOW) & (tt <= obs)
    if t[np.argmax(tail)]/TC < 25*max(TM):
        raise ValueError('run too short for the largest tau_m to relax before the window')
    tm = np.array(TM)[:, None]*TC                       # steps
    e = np.exp(-dts/tm/SUB)                             # (ntm, 1)
    w = (np.arange(SUB) + .5)/SUB
    m = np.repeat(mtr[0][None], len(TM), 0)             # (ntm, N)
    acc = dict(n=0, m=np.zeros(len(TM)), m2=np.zeros(len(TM)), s=0.,
               h={mc: np.zeros(len(TM)) for mc in MCS}, above={mc: np.zeros(len(TM)) for mc in MCS})
    for i in range(len(P) - 1):
        dP = P[i + 1] - P[i]
        for x in w:
            s = s_of(P[i] + x*dP, pmem, wp)
            m = s + (m - s)*e
        if tail[i + 1]:
            acc['n'] += 1; acc['m'] += m.mean(1); acc['m2'] += (m*m).mean(1); acc['s'] += s_of(P[i + 1], pmem, wp).mean()
            for mc in MCS:
                acc['h'][mc] += h_of(m, mc, hn).mean(1); acc['above'][mc] += (m > mc).mean(1)
    n = acc['n']
    mm = acc['m']/n
    zeta, z0 = float(p['zeta']), float(p['zeta0_frac'])
    open_loop = int(p.get('open_loop', 0))
    ftt = (field['t'] - freeze)/TC; ftail = (ftt >= obs - WINDOW) & (ftt <= obs)
    a = float(p['zeta_open'])/zeta if open_loop else z0 + (1 - z0)*float(field['chi_mean'][ftail].mean())
    row = dict(case=root.name, open_loop=open_loop, activity=a, x_in=(a - z0)/(1 - z0), run_tm_over_tc=float(p['tau_m'])/TC,
               run_tchi_over_tc=float(p['tau_chi'])/TC, n_tracer=int(P.shape[1]), n_sample_window=int(n), window_tc=WINDOW,
               tm_grid=list(TM), s_mean=float(acc['s']/n), sigma_P=float(P[tail].std()), mean_P=float(P[tail].mean()),
               m_mean=mm.tolist(), m_std=np.sqrt(np.maximum(acc['m2']/n - mm*mm, 0)).tolist(),
               chi={f'{mc:g}': (acc['h'][mc]/n).tolist() for mc in MCS},
               share_m_above_mc={f'{mc:g}': (acc['above'][mc]/n).tolist() for mc in MCS},
               field_chi=float(field['chi_mean'][ftail].mean()), field_m=float(field['m_mean'][ftail].mean()),
               field_sigma_P=float(field['P_std'][ftail].mean()))
    out.mkdir(parents=True, exist_ok=True)
    (out/'tracer_map.json').write_text(json.dumps(row, indent=1, allow_nan=False) + '\n')
    k = list(TM).index(8)
    print(json.dumps(dict(case=row['case'], x_in=row['x_in'], chi_tm8=row['chi']['0.2287'][k], field_chi=row['field_chi'])))


def collect(results_root, out, sub):
    rows = [json.loads(f.read_text()) for f in sorted(results_root.glob(f'**/{sub}/tracer_map.json'))]
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
