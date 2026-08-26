#!/usr/bin/env python
"""Field views for the confluent-wet model: velocity, phenotype chi, memory m, pressure P.

One-off code kept out of the plot library on purpose: the library's routines assume the
lyotropic phase field phi (visualization, counts_vs_t) or read the velocity from ff, i.e.
the bare lattice-Boltzmann moment. This model has no phi, and its physical velocity is
ux_mat/uy_mat = u_code + F/(2n), which near a melted defect core differs from the ff
moment by up to ~50%.

Panels, per frame:
  |u| with streamlines   the material velocity, log-free linear scale
  chi                    phenotype, fixed [0,1] so cases are comparable at a glance
  m                      memory, fixed [0,1] for the same reason
  P                      mechanical pressure -1/2 Tr(Pi), symmetric diverging scale

Usage: confluent_wet_fields.py <case_dir> <outfile.png> [frame_index]
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archive.archive import loadarchive
import plot.confluent as pc


def grid(frame, name):
    lx, ly = frame.parameters["LX"], frame.parameters["LY"]
    return np.array(getattr(frame, name)).reshape((lx, ly))


def _finish(ax, cax, label):
    ax.set_xticks([])
    ax.set_yticks([])
    cb = plt.colorbar(cax, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label(label)


def velocity_panel(ax, frame):
    # the library routine owns the material-velocity choice and the transposes
    cax = pc.velocity(frame, engine=ax, colorbar=False)
    _finish(ax, cax, r"$|u|$")
    ax.set_title(rf"velocity   $u_{{\rm rms}}$ = {pc.velocity_rms(frame):.2e}", fontsize=9)


def scalar_panel(ax, frame, name, cmap, label, vmin=0.0, vmax=1.0,
                 narrow=0.1):
    """chi and m on the fixed [0,1] scale, so cases are comparable at a glance.

    When the field spans less than `narrow` of that range the fixed scale hides everything
    -- a genuinely varying field then looks identical to a frozen one -- so fall back to a
    diverging map centred on the mean and say so in the title. Distinguishing "flat because
    it is frozen" from "flat because the colour scale is wrong" is the whole point.
    """
    f = grid(frame, name)
    span = float(f.max() - f.min())
    if span < narrow and span > 0:
        h = span/2
        cax = ax.imshow((f - f.mean()).T, origin="lower", interpolation="nearest",
                        cmap="RdBu_r", vmin=-h, vmax=h)
        _finish(ax, cax, label + r" $-\ \langle$" + label + r"$\rangle$")
        ax.set_title(f"{name}   [{f.min():.4f}, {f.max():.4f}]   "
                     f"(auto scale, span {span:.1e})", fontsize=9)
    else:
        cax = ax.imshow(f.T, origin="lower", interpolation="nearest",
                        cmap=cmap, vmin=vmin, vmax=vmax)
        _finish(ax, cax, label)
        ax.set_title(f"{name}   [{f.min():.4f}, {f.max():.4f}]", fontsize=9)


def pressure_panel(ax, frame):
    P = grid(frame, "pressure")
    # symmetric about the spatial median: in a periodic box <P> is pinned by the initial
    # density, not by the dynamics, so only P - <P> carries information
    P = P - np.median(P)
    v = np.percentile(np.abs(P), 99.5)
    v = v if v > 0 else 1.0
    cax = ax.imshow(P.T, origin="lower", interpolation="nearest",
                    cmap="coolwarm", vmin=-v, vmax=v)
    _finish(ax, cax, r"$P - \mathrm{med}(P)$")
    ax.set_title(rf"pressure   $\sigma_P$ = {P.std():.2e}", fontsize=9)


def figure(case_dir, outfile, index=None):
    ar = loadarchive(case_dir)
    p = ar.parameters
    n = ar._nframes
    i = n - 1 if index is None else index
    fr = ar[i]

    fig, axes = plt.subplots(1, 4, figsize=(17, 4.4))
    velocity_panel(axes[0], fr)
    scalar_panel(axes[1], fr, "chi", "viridis", r"$\chi$")
    scalar_panel(axes[2], fr, "m", "cividis", r"$m$")
    pressure_panel(axes[3], fr)

    xi_N = np.sqrt(p["LL"]/(2*p["CC"]))
    la = np.sqrt(p["LL"]/p["zeta"]) if p["zeta"] > 0 else float("inf")
    fig.suptitle(
        f"{os.path.basename(case_dir)}   t = {p['nstart'] + i*p['ninfo']}   "
        f"$\\zeta$={p['zeta']:g}, CC={p['CC']:g}, LL={p['LL']:g}  "
        f"($\\xi_N$={xi_N:.2f}, $\\ell_a$={la:.1f})   "
        f"$\\tau_\\chi$={p['tau_chi']:g}, $\\tau_m$={p['tau_m']:g}",
        fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    fig.savefig(outfile, dpi=170)
    plt.close(fig)
    print(f"wrote {outfile}")


if __name__ == "__main__":
    figure(sys.argv[1], sys.argv[2],
           int(sys.argv[3]) if len(sys.argv) > 3 else None)
