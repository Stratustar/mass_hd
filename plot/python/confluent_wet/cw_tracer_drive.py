#!/usr/bin/env python3
"""The single cell driven by the REAL pressure histories of the nematic model (hill law, beta = 0).

A tracer run (ntracer > 0) stores P and m along N material paths for the whole run. Along a path
the model's memory and phenotype obey

    D_t m   = (s(P) - m)/tau_m       + Dm   lap m
    D_t chi = (h(m) - chi)/tau_chi   + Dchi lap chi      (chi reaction off before chi-freeze-steps)

so a cell that carries its own m and chi and feels the path's P(t) is the nematic model minus the
two diffusion terms. Each path is integrated from the run's own start (m = m0, chi = chi0, the chi
reaction released at chi-freeze-steps), i.e. through the same protocol as the fields, at four
levels:

  T1        m_path from s(P(t)) alone; chi reads h(m_path). The single cell with the real pressure
            and NO diffusion closure.
  T2(c)     the same m_path, but chi reads the single cell's shrunk memory
            m_eff = mbar + (m_path - mbar)/sqrt(c), mbar = the instantaneous mean over the paths;
            on a grid of c (so any c_D table can be applied afterwards) and at
            c = c_meas = Var(m_path)/Var(m_tracer), both variances about the instantaneous mean,
            in the terminal window.
  T2m       T2 at c_meas, but shrunk about the MODEL's mean memory on the same paths,
            m_eff = <m_tracer> + (m_path - mbar)/sqrt(c_meas): mean and variance of the memory both
            taken from the model, only the shape and the time structure left to the path. Needed
            because diffusion also lowers the mean of m along a path (the paths gather where m is
            high), which a shrink about mbar cannot give.
  R         chi reads h(m_tracer), the model's own memory at the same material points: only
            Dchi is dropped. The reference that carries the tracers' sampling of the layer.
  field     the area means of video_meta.csv in the same window.

P is linear between samples, s(P) is evaluated at SUB sub-step midpoints and each sub-step is the
exact exponential for a constant source; chi takes one exact exponential step per sample with the
target at the two ends averaged.

  cw_tracer_drive.py <case_in> <case_out>          -> tracer_drive.json, tracer_drive.npz
  cw_tracer_drive.py <results_root> <out.json> --collect <out_subdir>
                                                     gather every <case>/<out_subdir>/tracer_drive.json
"""
import argparse
import io
import json
from pathlib import Path

import numpy as np

from cw_mc_tm_scan import TC, parameters

SUB = 4
C_GRID = np.round(np.geomspace(1., 30., 35), 4)
BINS = np.linspace(0., 1., 201)


def s_of(p, pmem, wp):
    if wp == 0:
        return (p > pmem).astype(float)
    return 1/(1 + np.exp(-np.clip((p - pmem)/wp, -700, 700)))


def h_of(m, mc, hn):
    if hn == 0:
        return (m >= mc).astype(float)
    with np.errstate(divide='ignore', over='ignore'):
        return np.where(m > 0, 1/(1 + (mc/np.maximum(m, 1e-300))**hn), 0.)


def load(root):
    p = parameters(root)
    if p.get('switch_model') != 'hill' or float(p.get('beta', 0)) != 0:
        raise ValueError('needs a switch-model = hill run with beta = 0')
    if int(p.get('mem_freeze_steps', 0)) != 0:
        raise ValueError('mem-freeze-steps is not handled')
    text = '\n'.join(x for x in (root/'tracer_meta.csv').read_text().splitlines() if not x.startswith('#'))
    meta = np.genfromtxt(io.StringIO(text), delimiter=',', names=True)
    n = int(meta['n_tracer'][0])
    raw = np.fromfile(root/'tracer.f32', np.float32)
    ns = min(len(raw)//(4*n), len(meta))
    raw = raw[:ns*4*n].reshape(ns, n, 4)
    t = meta['t'][:ns]
    if np.any(np.diff(t) != t[1] - t[0]):
        raise ValueError('irregular tracer clock')
    text = '\n'.join(x for x in (root/'video_meta.csv').read_text().splitlines() if not x.startswith('#'))
    field = np.genfromtxt(io.StringIO(text), delimiter=',', names=True)
    return p, t, raw[..., 0].astype(float), raw[..., 1].astype(float), field


def memory_path(P, m0, dt_tm, pmem, wp):
    """m along each path from s(P) alone; P (nsample, N), dt_tm = sample interval / tau_m."""
    m = np.empty_like(P)
    m[0] = m0
    e = np.exp(-dt_tm/SUB)
    w = (np.arange(SUB) + .5)/SUB
    cur = m[0].copy()
    for i in range(len(P) - 1):
        dP = P[i + 1] - P[i]
        for x in w:
            s = s_of(P[i] + x*dP, pmem, wp)
            cur = s + (cur - s)*e
        m[i + 1] = cur
    return m


def chi_path(target, chi0, dt_tchi, release):
    """chi along each path for a target history (nsample, ...); the reaction starts at sample `release`."""
    e = np.exp(-dt_tchi)
    chi = np.empty_like(target)
    chi[0] = chi0
    for i in range(len(target) - 1):
        if i < release:
            chi[i + 1] = chi[i]
        else:
            T = .5*(target[i] + target[i + 1])
            chi[i + 1] = T + (chi[i] - T)*e
    return chi


def reduce_case(root, out):
    p, t, P, mtr, field = load(root)
    tm, tchi = float(p['tau_m']), float(p['tau_chi'])
    mc, pmem, wp, hn = float(p['mc']), float(p['pmem']), float(p['wp']), float(p['h_n'])
    freeze, nsteps = int(p['chi_freeze_steps']), int(p['nsteps'])
    dts = float(t[1] - t[0])
    release = int(np.ceil(freeze/dts))
    obs = (nsteps - freeze)/TC
    window = max(500., 20*tm/TC)                       # the reducer's terminal window, tau_c
    tt = (t - freeze)/TC                               # observation time, tau_c
    tail = (tt >= obs - window) & (tt <= obs)
    ftt = (field['t'] - freeze)/TC
    ftail = (ftt >= obs - window) & (ftt <= obs)

    mpath = memory_path(P, mtr[0], dts/tm, pmem, wp)
    dev = lambda a: a - a.mean(1, keepdims=True)
    v_path, v_tr = float(np.mean(dev(mpath[tail])**2)), float(np.mean(dev(mtr[tail])**2))
    c_meas = v_path/v_tr
    mbar = mpath.mean(1, keepdims=True)

    chi0 = float(p['chi0'])
    chi_T1 = chi_path(h_of(mpath, mc, hn), chi0, dts/tchi, release)
    chi_R = chi_path(h_of(mtr, mc, hn), chi0, dts/tchi, release)
    c_all = np.append(C_GRID, c_meas)
    chi_T2 = np.empty((len(c_all), len(t)))            # instantaneous means only
    std_T2 = np.empty(len(c_all))
    for k, c in enumerate(c_all):
        x = chi_path(h_of(mbar + (mpath - mbar)/np.sqrt(c), mc, hn), chi0, dts/tchi, release)
        chi_T2[k] = x.mean(1)
        std_T2[k] = x[tail].std()
    chi_T2m = chi_path(h_of(mtr.mean(1, keepdims=True) + (mpath - mbar)/np.sqrt(c_meas), mc, hn), chi0, dts/tchi, release)

    def moments(a):
        a = a[tail]
        d = a - a.mean()
        return dict(mean=float(a.mean()), std=float(a.std()), skew=float(np.mean(d**3)/a.std()**3))

    s_tr = s_of(P[tail], pmem, wp)
    row = dict(
        case=root.name, tm_over_tc=tm/TC, tchi_over_tc=tchi/TC, mc=mc, pmem=pmem, wp=wp, h_n=hn, chi0=chi0,
        n_tracer=int(P.shape[1]), n_sample=int(len(t)), sample_tc=dts/TC, sub=SUB,
        observation_tc=obs, window_tc=window,
        field=dict(chi=float(field['chi_mean'][ftail].mean()), chi_std=float(field['chi_std'][ftail].mean()),
                   m=float(field['m_mean'][ftail].mean()), m_std=float(field['m_std'][ftail].mean()),
                   P_std=float(field['P_std'][ftail].mean())),
        tracer=dict(P=moments(P), m=moments(mtr), s_mean=float(s_tr.mean()), share_above_pmem=float((P[tail] > pmem).mean()),
                    h_of_m=float(h_of(mtr[tail], mc, hn).mean()), share_m_above_mc=float((mtr[tail] > mc).mean())),
        m_path=dict(**moments(mpath), share_above_mc=float((mpath[tail] > mc).mean()),
                    corr_with_tracer_m=float(np.corrcoef(dev(mpath[tail]).ravel(), dev(mtr[tail]).ravel())[0, 1])),
        var_m_path=v_path, var_m_tracer=v_tr, c_meas=c_meas,
        chi_T1=dict(mean=float(chi_T1[tail].mean()), std=float(chi_T1[tail].std())),
        chi_R=dict(mean=float(chi_R[tail].mean()), std=float(chi_R[tail].std())),
        chi_T2_meas=dict(mean=float(chi_T2[-1][tail].mean()), std=float(std_T2[-1])),
        chi_T2m=dict(mean=float(chi_T2m[tail].mean()), std=float(chi_T2m[tail].std())),
        c_grid=C_GRID.tolist(), chi_T2_grid=[float(x[tail].mean()) for x in chi_T2[:-1]],
        chi_T2_grid_std=std_T2[:-1].tolist())
    out.mkdir(parents=True, exist_ok=True)
    (out/'tracer_drive.json').write_text(json.dumps(row, indent=1, allow_nan=False) + '\n')
    k = max(1, int(round(.5*TC/dts)))                  # series every ~.5 tau_c
    np.savez_compressed(
        out/'tracer_drive.npz', t_tc=tt[::k], chi_T1=chi_T1.mean(1)[::k], chi_R=chi_R.mean(1)[::k],
        chi_T2_meas=chi_T2[-1][::k], chi_T2m=chi_T2m.mean(1)[::k], m_path=mpath.mean(1)[::k], m_tracer=mtr.mean(1)[::k],
        field_t_tc=ftt, field_chi=field['chi_mean'], field_m=field['m_mean'], bins=BINS,
        hist_m_path=np.histogram(mpath[tail], BINS)[0], hist_m_tracer=np.histogram(mtr[tail], BINS)[0],
        hist_chi_T1=np.histogram(chi_T1[tail], BINS)[0], hist_chi_R=np.histogram(chi_R[tail], BINS)[0])
    print(json.dumps({k: row[k] for k in ('case', 'c_meas')} | {
        'chi': {'field': row['field']['chi'], 'R': row['chi_R']['mean'], 'T1': row['chi_T1']['mean'],
                'T2_meas': row['chi_T2_meas']['mean'], 'T2m': row['chi_T2m']['mean']}}))


def collect(results_root, out, sub):
    rows = [json.loads(f.read_text()) for f in sorted(results_root.glob(f'**/{sub}/tracer_drive.json'))]
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
