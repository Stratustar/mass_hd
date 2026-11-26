#!/usr/bin/env python3
"""One run of the half active / half passive start (cases/20261125/cw_halfhalf, cases/20261126/cw_halfhalf_tm_tchi; scripts/gen_cw_halfhalf.py): a small
file with the interface dynamics, and a video.

At t = 0 the active band is LX/4 <= x < 3 LX/4 (chi = 1), the rest is passive (chi = 0); then everything evolves.
d is the signed distance to the nearest initial boundary, d > 0 inside the initially active band; "deep" is |d| > LX/8.

halfhalf.npz
  t                      frame times in tau_c
  prof_{chi,m,P,speed}   averages over y at every frame, as functions of x (speed in lattice units per tau_c)
  chi_mean               the active share of the box
  {chi,m}_deep_{active,passive}, speed_deep_*, sdP_deep_*
                         means (rms for the speed, standard deviation for P) deep inside each initial band
  boundary               number of lattice links across chi = 1/2, in units of the two straight initial boundaries
  snap, snap_t           chi maps at SNAP times (tau_c)
  mc, z0, tau_m, tau_chi of the run (times in tau_c)
halfhalf.mp4             |u|, P - <P>, m, chi in the colours of the confluent-wet dashboard spec; one video frame per saved frame

  cw_halfhalf.py <case_in> <case_out> [--no-video]
"""
import argparse
import json
from pathlib import Path

import numpy as np

TC = 674.3290333006435
SNAP = (0, 5, 10, 20, 40, 80, 150, 300)
KEYS = ('chi', 'm', 'pressure', 'ux_mat', 'uy_mat')


def read_params(root):
    d = json.loads((root/'parameters.json').read_text())
    d = d.get('data', d)
    return {k: (v['value'] if isinstance(v, dict) and 'value' in v else v) for k, v in d.items()}


def frames(root):
    return sorted(root.glob('frame*.json'), key=lambda f: int(f.stem[5:]))


def load(f, LX, LY):
    d = json.loads(f.read_text())['data']
    return {k: np.asarray(d[k]['value'], float).reshape(LX, LY) for k in KEYS}      # axis 0 is x


def boundary_length(chi):
    a = chi > .5
    return (np.count_nonzero(a != np.roll(a, 1, 0)) + np.count_nonzero(a != np.roll(a, 1, 1)))/(2*chi.shape[1])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('case_in', type=Path)
    ap.add_argument('case_out', type=Path)
    ap.add_argument('--no-video', action='store_true')
    a = ap.parse_args()
    p = read_params(a.case_in)
    LX, LY = int(p['LX']), int(p['LY'])
    mc, z0, tm, tchi = float(p['mc']), float(p['zeta0_frac']), float(p['tau_m'])/TC, float(p['tau_chi'])/TC
    fr = frames(a.case_in)
    assert fr, f'no frames in {a.case_in}'
    a.case_out.mkdir(parents=True, exist_ok=True)
    t = np.array([int(f.stem[5:]) for f in fr])/TC
    d = LX/4 - np.abs(np.arange(LX) - LX/2)
    A, Pz = d > LX/8, d < -LX/8
    prof = {k: np.zeros((len(fr), LX), np.float32) for k in ('chi', 'm', 'P', 'speed')}
    names = ('chi_mean', 'chi_deep_active', 'chi_deep_passive', 'm_deep_active', 'm_deep_passive', 'boundary',
             'speed_deep_active', 'speed_deep_passive', 'sdP_deep_active', 'sdP_deep_passive')
    ser = {k: np.zeros(len(fr)) for k in names}
    want = sorted({int(np.argmin(np.abs(t - s))) for s in SNAP if s <= t[-1] + .5})
    snaps = []

    writer = None
    if not a.no_video:
        import imageio_ffmpeg
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.animation import FFMpegWriter
        plt.rcParams.update({'font.family': 'DejaVu Sans', 'figure.facecolor': 'white', 'savefig.facecolor': 'white',
                             'animation.ffmpeg_path': imageio_ffmpeg.get_ffmpeg_exe()})
        # The confluent-wet dashboard spec (dev_log, default rules): u, P, m, chi in this order; magma 0 .. p99.5,
        # RdBu_r +- p99.5, viridis 0 .. 1, coolwarm 0 .. 1; limits sampled from ~10 frames and then fixed; a
        # one-symbol title and a colour bar per panel. The only addition is the clock in the bottom margin.
        pick = [fr[j] for j in np.unique(np.linspace(0, len(fr) - 1, 10).round().astype(int))]
        smp = [load(f, LX, LY) for f in pick]
        umax = float(np.percentile(np.concatenate([np.hypot(q['ux_mat'], q['uy_mat']).ravel() for q in smp])*TC, 99.5))
        pmax = float(np.percentile(np.abs(np.concatenate([(q['pressure'] - q['pressure'].mean()).ravel() for q in smp])), 99.5))
        fig = plt.figure(figsize=(16.4, 4.35), dpi=100)
        blank = np.zeros((LY, LX))
        spec = (('magma', 0, umax, r'$|u|$'), ('RdBu_r', -pmax, pmax, r'$P$'), ('viridis', 0, 1, r'$m$'), ('coolwarm', 0, 1, r'$\chi$'))
        ims = []
        for j, (cm, lo, hi, s_) in enumerate(spec):
            x0 = .008 + j*.248
            x = fig.add_axes([x0, .075, .205, .205*16.4/4.35])
            ims.append(x.imshow(blank, origin='lower', cmap=cm, vmin=lo, vmax=hi, interpolation='nearest'))
            x.set_xticks([]); x.set_yticks([]); x.set_title(s_, fontsize=14, pad=5)
            for sp_ in x.spines.values():
                sp_.set_visible(False)
            cb = fig.colorbar(ims[-1], cax=fig.add_axes([x0 + .209, .075, .006, .205*16.4/4.35]))
            cb.ax.tick_params(labelsize=8.5)
        clock = fig.text(.008, .02, '', fontsize=10.5)
        writer = FFMpegWriter(fps=12, codec='libx264', extra_args=['-pix_fmt', 'yuv420p', '-crf', '20', '-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2'])
        writer.setup(fig, str(a.case_out/'halfhalf.mp4'), dpi=100)

    for i, f in enumerate(fr):
        q = load(f, LX, LY)
        sp = np.hypot(q['ux_mat'], q['uy_mat'])*TC
        prof['chi'][i], prof['m'][i], prof['P'][i], prof['speed'][i] = q['chi'].mean(1), q['m'].mean(1), q['pressure'].mean(1), sp.mean(1)
        ser['chi_mean'][i] = q['chi'].mean()
        ser['chi_deep_active'][i], ser['chi_deep_passive'][i] = q['chi'][A].mean(), q['chi'][Pz].mean()
        ser['m_deep_active'][i], ser['m_deep_passive'][i] = q['m'][A].mean(), q['m'][Pz].mean()
        ser['speed_deep_active'][i], ser['speed_deep_passive'][i] = np.sqrt((sp[A]**2).mean()), np.sqrt((sp[Pz]**2).mean())
        ser['sdP_deep_active'][i], ser['sdP_deep_passive'][i] = q['pressure'][A].std(), q['pressure'][Pz].std()
        ser['boundary'][i] = boundary_length(q['chi'])
        if i in want:
            snaps.append(q['chi'].astype(np.float32))
        if writer is not None:
            P = q['pressure']
            ims[0].set_data(sp.T); ims[1].set_data((P - P.mean()).T); ims[2].set_data(q['m'].T); ims[3].set_data(q['chi'].T)
            clock.set_text(rf't = {t[i]:5.1f} $\tau_c$')
            writer.grab_frame()
    if writer is not None:
        writer.finish()
    np.savez_compressed(a.case_out/'halfhalf.npz', t=t, snap=np.array(snaps), snap_t=t[want], mc=mc, z0=z0, tau_m=tm, tau_chi=tchi,
                        **{f'prof_{k}': v for k, v in prof.items()}, **ser)
    print(f'{a.case_in.name}: {len(fr)} frames to {t[-1]:.1f} tau_c; <chi> {ser["chi_mean"][0]:.3f} -> {ser["chi_mean"][-1]:.3f}; '
          f'boundary {ser["boundary"][0]:.2f} -> {ser["boundary"][-1]:.2f}')


if __name__ == '__main__':
    main()
