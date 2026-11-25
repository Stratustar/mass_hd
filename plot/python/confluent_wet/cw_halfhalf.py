#!/usr/bin/env python3
"""One run of the half active / half passive start (cases/20261125/cw_halfhalf, scripts/gen_cw_halfhalf.py): a small
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
  mc, z0, tau_m          of the run (tau_m in tau_c)
halfhalf.mp4             chi, m, P - <P>, flow speed; one video frame per saved frame

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
    mc, z0, tm = float(p['mc']), float(p['zeta0_frac']), float(p['tau_m'])/TC
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
        fig, ax = plt.subplots(1, 4, figsize=(15.6, 4.5), dpi=100)
        fig.subplots_adjust(left=.006, right=.994, top=.86, bottom=.07, wspace=.03)
        blank = np.zeros((LY, LX))
        ims = [ax[0].imshow(blank, origin='lower', cmap='RdBu_r', vmin=0, vmax=1, interpolation='nearest'),
               ax[1].imshow(blank, origin='lower', cmap='viridis', vmin=0, vmax=.4, interpolation='nearest'),
               ax[2].imshow(blank, origin='lower', cmap='RdBu_r', vmin=-60, vmax=60, interpolation='nearest'),
               ax[3].imshow(blank, origin='lower', cmap='magma', vmin=0, vmax=15, interpolation='nearest')]
        titles = (r'phenotype $\chi$  (red = active, blue = passive)', 'memory m  (0 to 0.4)',
                  r'pressure $P - \langle P\rangle$  ($\pm$60e-3)', r'flow speed  (0 to 15 per $\tau_c$)')
        for x, s in zip(ax, titles):
            x.set_xticks([]); x.set_yticks([]); x.set_title(s, fontsize=11.5, pad=5)
            for sp in x.spines.values():
                sp.set_visible(False)
        fig.text(.008, .955, rf'half active / half passive start;  $\tau_m$ = {tm:g} $\tau_c$,  $m_c$ = {mc:g},  $z_0$ = {z0:g}',
                 fontsize=12.5, va='center')
        clock = fig.text(.008, .025, '', fontsize=11.5)
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
            ims[0].set_data(q['chi'].T); ims[1].set_data(q['m'].T); ims[2].set_data(((P - P.mean())*1e3).T); ims[3].set_data(sp.T)
            clock.set_text(rf't = {t[i]:5.1f} $\tau_c$ = {t[i]/tm:.2f} $\tau_m$      active share $\langle\chi\rangle$ = {ser["chi_mean"][i]:.2f}')
            writer.grab_frame()
    if writer is not None:
        writer.finish()
    np.savez_compressed(a.case_out/'halfhalf.npz', t=t, snap=np.array(snaps), snap_t=t[want], mc=mc, z0=z0, tau_m=tm,
                        **{f'prof_{k}': v for k, v in prof.items()}, **ser)
    print(f'{a.case_in.name}: {len(fr)} frames to {t[-1]:.1f} tau_c; <chi> {ser["chi_mean"][0]:.3f} -> {ser["chi_mean"][-1]:.3f}; '
          f'boundary {ser["boundary"][0]:.2f} -> {ser["boundary"][-1]:.2f}')


if __name__ == '__main__':
    main()
