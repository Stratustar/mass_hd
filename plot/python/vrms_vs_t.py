#!/usr/bin/env python3
"""RMS velocity vs time across the switch/defect arms (bulk_switchdefect_zOchi).

For each (zeta, Ochi) variant we overlay v_rms(t) for arm1/arm2/arm3:
  arm1 = switching ON, fresh init, 50k steps            (t = 0 .. 50000)
  arm2 = switching OFF, structured chi advected         (branches from arm1 frame20000)
  arm3 = switching OFF, uniform phi-weighted-mean chi   (branches from arm1 frame20000)

arm2/arm3 restart their step counter at 0, so they are shifted by +ARM_OFFSET
(=20000) to share arm1's clock. NB arm2 carries the known m-conservation bug
([0,phi] clamp at sharp chi interfaces) and is drawn dashed as a caveat.

v_rms is the unweighted full-lattice RMS of the LB velocity, reconstructed from
the distribution field `ff` exactly as velocity_divergence.py / plot_hd. The
runs are config=full (uniform phi fill), so there is no colony mask to apply.

Usage:
  python vrms_vs_t.py --base BASE --outdir OUTDIR [--arms arm1 arm2 arm3]

BASE is the scratch dir holding arm0/arm1/arm2/arm3/<variant>/ trees, e.g.
  /scratch/helu/mass_hd/cases/06jul/bulk_switchdefect_zOchi
"""

import matplotlib
matplotlib.use("Agg")

import argparse
import glob
import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archive.archive import loadarchive

ARM_OFFSET = {"arm1": 0, "arm2": 20000, "arm3": 20000}
ARM_STYLE = {
    "arm1": dict(color="#1f77b4", ls="-",  label="arm1 (switch ON)"),
    "arm2": dict(color="#d62728", ls="--", label="arm2 (switch OFF, structured; m-bug)"),
    "arm3": dict(color="#2ca02c", ls="-",  label="arm3 (switch OFF, uniform chi)"),
}


def velocity_grids(frame):
    """Full-lattice (vx, vy) as (LX, LY) arrays, from the LB field ff."""
    lx = frame.parameters["LX"]
    ly = frame.parameters["LY"]
    ff = np.array(frame.ff)
    density = np.sum(ff, axis=1)
    density = np.where(density == 0.0, np.nan, density)
    vx = (ff[:, 1] - ff[:, 2] + ff[:, 5] - ff[:, 6] - ff[:, 7] + ff[:, 8]) / density
    vy = (ff[:, 3] - ff[:, 4] + ff[:, 5] - ff[:, 6] + ff[:, 7] - ff[:, 8]) / density
    vx = np.nan_to_num(vx).reshape((lx, ly))
    vy = np.nan_to_num(vy).reshape((lx, ly))
    return vx, vy


def last_present_index(ar, archive_dir):
    """Highest frame index whose json file actually exists (survive blow-ups)."""
    steps = []
    for f in glob.glob(os.path.join(archive_dir, "frame*.json")):
        m = re.search(r"frame(\d+)\.json$", os.path.basename(f))
        if m:
            steps.append(int(m.group(1)))
    nominal = int((ar.nsteps - ar.nstart) / ar.ninfo)
    if not steps:
        return nominal
    return min(nominal, int((max(steps) - ar.nstart) / ar.ninfo))


def vrms_series(archive_dir, offset=0):
    """(steps, vrms) for one archive; steps shifted onto arm1's clock by offset."""
    ar = loadarchive(archive_dir)
    imax = last_present_index(ar, archive_dir)
    steps, vrms = [], []
    for i in range(imax + 1):
        try:
            frame = ar.read_frame(i)
        except (ValueError, FileNotFoundError, OSError) as exc:
            print(f"       {archive_dir} frame {i} unreadable ({exc}); stopping")
            break
        vx, vy = velocity_grids(frame)
        vrms.append(float(np.sqrt(np.mean(vx**2 + vy**2))))
        steps.append(ar.nstart + i * ar.ninfo + offset)
    return np.array(steps), np.array(vrms)


def parse_variant(name):
    """z0p15_O0p03 -> (zeta, Ochi) floats for grid placement/sorting."""
    m = re.match(r"z([0-9p]+)_O([0-9p]+)", name)
    if not m:
        return (float("nan"), float("nan"))
    z = float(m.group(1).replace("p", "."))
    o = float(m.group(2).replace("p", "."))
    return (z, o)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--base", required=True, help="dir holding arm*/<variant>/ trees")
    p.add_argument("--outdir", required=True, help="where to write the figures")
    p.add_argument("--arms", nargs="+", default=["arm1", "arm2", "arm3"])
    args = p.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # variants come from arm1 (the superset); keep only those present in arm1
    a1 = os.path.join(args.base, "arm1")
    variants = sorted(
        [d for d in os.listdir(a1) if os.path.isdir(os.path.join(a1, d))],
        key=parse_variant)
    zetas = sorted(set(parse_variant(v)[0] for v in variants))
    ochis = sorted(set(parse_variant(v)[1] for v in variants))
    print(f"[vrms] {len(variants)} variants; zeta={zetas} Ochi={ochis}")

    # gather all series first (so we can share a y-axis per row/whole grid)
    data = {}  # (variant, arm) -> (steps, vrms)
    for v in variants:
        for arm in args.arms:
            d = os.path.join(args.base, arm, v)
            if not os.path.isdir(d):
                continue
            steps, vr = vrms_series(d, ARM_OFFSET.get(arm, 0))
            data[(v, arm)] = (steps, vr)
            print(f"[vrms] {arm:5s} {v:14s}: {len(steps):3d} frames, "
                  f"vrms final={vr[-1]:.3e} max={vr.max():.3e}")

    # 3x3 grid: rows = zeta, cols = Ochi
    nr, nc = len(zetas), len(ochis)
    fig, axes = plt.subplots(nr, nc, figsize=(4.2 * nc, 3.4 * nr),
                             sharex=True, squeeze=False)
    for v in variants:
        z, o = parse_variant(v)
        ax = axes[zetas.index(z)][ochis.index(o)]
        for arm in args.arms:
            if (v, arm) not in data:
                continue
            steps, vr = data[(v, arm)]
            ax.plot(steps, vr, **ARM_STYLE[arm], lw=1.4)
        ax.axvline(20000, color="grey", ls=":", lw=0.8, alpha=0.7)  # branch time
        ax.set_title(f"zeta={z:g}, Ochi={o:g}", fontsize=10)
        ax.grid(alpha=0.25)
    for i, z in enumerate(zetas):
        axes[i][0].set_ylabel(r"$v_{\rm rms}$")
    for j, o in enumerate(ochis):
        axes[-1][j].set_xlabel("step")
    axes[0][0].legend(fontsize=7, loc="upper left")
    fig.suptitle("RMS velocity vs time by arm  "
                 "(bulk_switchdefect_zOchi; arm2/arm3 aligned +20000)",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    out = os.path.join(args.outdir, "vrms_vs_t_grid.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[vrms] wrote {out}")


if __name__ == "__main__":
    main()
