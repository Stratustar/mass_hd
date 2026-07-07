#!/usr/bin/env python3
"""Final-frame velocity divergence (div v = d vx/dx + d vy/dy) for MASS archives.

Velocity is reconstructed from the LB distribution field `ff` exactly as
plot_hd.velocity_components does; the divergence is a centered finite difference
(np.gradient, 1 lattice unit spacing). For each archive we render a diverging,
0-centered heatmap cropped to the colony (phi>0.5 bbox), plus a combined
comparison figure across archives on a shared symmetric color scale so the
alpha-dependence of the growth-driven expansion is directly readable.

Usage:
  python velocity_divergence.py --outdir OUTDIR ARCHIVE_DIR [ARCHIVE_DIR ...]

Each ARCHIVE_DIR is a MASS output directory (the *.json frame dump). The label
of a case defaults to the directory basename (e.g. alpha0p0001).
"""

import matplotlib
matplotlib.use("Agg")

import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archive.archive import loadarchive

DIV_CMAP = "RdBu_r"        # positive (expansion) = red, negative (compression) = blue
CLIP_PCTILE = 99.5         # symmetric color limit from this percentile of |div v|
CROP_MARGIN = 30           # lattice cells of padding around the phi>0.5 bbox


def available_frame_count(ar):
    return int((ar.nsteps - ar.nstart) / ar.ninfo) + 1


def _grid(frame, name):
    lx = frame.parameters["LX"]
    ly = frame.parameters["LY"]
    return np.array(getattr(frame, name)).reshape((lx, ly))


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


def divergence(frame):
    """div v on the full lattice, axis 0 = x, axis 1 = y (1 cell spacing)."""
    vx, vy = velocity_grids(frame)
    dvx_dx = np.gradient(vx, axis=0)
    dvy_dy = np.gradient(vy, axis=1)
    return dvx_dx + dvy_dy


def colony_bbox(phi, margin=CROP_MARGIN):
    lx, ly = phi.shape
    mask = phi > 0.5
    if not mask.any():
        return (0, lx - 1, 0, ly - 1)
    xs = np.where(mask.any(axis=1))[0]
    ys = np.where(mask.any(axis=0))[0]
    x0 = max(0, int(xs.min()) - margin)
    x1 = min(lx - 1, int(xs.max()) + margin)
    y0 = max(0, int(ys.min()) - margin)
    y1 = min(ly - 1, int(ys.max()) + margin)
    return (x0, x1, y0, y1)


def load_case(archive_dir):
    ar = loadarchive(archive_dir)
    n = available_frame_count(ar)
    frame = ar.read_frame(n - 1)
    step = ar.nstart + (n - 1) * ar.ninfo
    div = divergence(frame)
    phi = _grid(frame, "phi")
    alpha = frame.parameters.get("alpha", float("nan"))
    return {
        "div": div,
        "phi": phi,
        "step": step,
        "alpha": alpha,
        "bbox": colony_bbox(phi),
    }


def sym_limit(arr_list, pctile=CLIP_PCTILE):
    vals = np.concatenate([np.abs(a).ravel() for a in arr_list])
    vmax = np.percentile(vals, pctile)
    if not np.isfinite(vmax) or vmax == 0.0:
        vmax = max((np.abs(a).max() for a in arr_list), default=1.0) or 1.0
    return vmax


def draw_panel(ax, case, vmax, title):
    x0, x1, y0, y1 = case["bbox"]
    div = case["div"][x0:x1 + 1, y0:y1 + 1]
    phi = case["phi"][x0:x1 + 1, y0:y1 + 1]
    im = ax.imshow(div.T, origin="lower", cmap=DIV_CMAP, vmin=-vmax, vmax=vmax,
                   extent=(x0, x1, y0, y1), interpolation="nearest")
    # colony outline for context
    ax.contour(np.linspace(x0, x1, phi.shape[0]),
               np.linspace(y0, y1, phi.shape[1]),
               phi.T, levels=[0.5], colors="k", linewidths=0.6, alpha=0.6)
    ax.set_title(title, fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_aspect("equal")
    return im


def case_title(case, label):
    return (f"{label}\n"
            f"alpha={case['alpha']:g}, step={case['step']}\n"
            f"max|div v|={np.abs(case['div']).max():.2e}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("archives", nargs="+", help="MASS output dirs (one per case)")
    p.add_argument("--outdir", required=True, help="where to write the figures")
    p.add_argument("--labels", nargs="*", default=None,
                   help="labels per archive (default: dir basename)")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    labels = args.labels or [os.path.basename(os.path.normpath(a)) for a in args.archives]
    if len(labels) != len(args.archives):
        p.error("number of --labels must match number of archives")

    cases = []
    for archive_dir, label in zip(args.archives, labels):
        print(f"[divv] loading {label}: {archive_dir}")
        case = load_case(archive_dir)
        cases.append((label, case))
        print(f"       alpha={case['alpha']:g} step={case['step']} "
              f"max|div v|={np.abs(case['div']).max():.3e} "
              f"mean(div v over phi>0.5)={case['div'][case['phi']>0.5].mean():.3e}")

    vmax = sym_limit([c["div"] for _, c in cases])
    print(f"[divv] shared symmetric color limit (p{CLIP_PCTILE}) = {vmax:.3e}")

    # per-case single figures
    for label, case in cases:
        fig, ax = plt.subplots(figsize=(5, 4.6))
        im = draw_panel(ax, case, vmax, case_title(case, label))
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label(r"$\nabla\!\cdot\mathbf{v}$")
        fig.tight_layout()
        out = os.path.join(args.outdir, f"divv_{label}_frame{case['step']}.png")
        fig.savefig(out, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"[divv] wrote {out}")

    # combined comparison (shared scale)
    n = len(cases)
    fig, axes = plt.subplots(1, n, figsize=(4.4 * n, 4.6))
    if n == 1:
        axes = [axes]
    im = None
    for ax, (label, case) in zip(axes, cases):
        im = draw_panel(ax, case, vmax, case_title(case, label))
    cb = fig.colorbar(im, ax=axes, fraction=0.046 / n, pad=0.02)
    cb.set_label(r"$\nabla\!\cdot\mathbf{v}$  (shared scale)")
    fig.suptitle("final-frame velocity divergence vs alpha (all-grow, switching off)",
                 fontsize=12)
    out = os.path.join(args.outdir, "divv_compare.png")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[divv] wrote {out}")


if __name__ == "__main__":
    main()
