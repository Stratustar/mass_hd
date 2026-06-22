"""Center mechanical pressure versus time for one simulation archive.

Center pressure = p(L/2, L/2) = -sigma_bulk at the box-center node each frame
(the colony is initialised centered, so the box center is the colony center).
Writes center_pressure_vs_t.png and center_pressure_vs_t.csv into outdir.

Usage:
    python center_pressure_vs_t.py <inputdir> <outdir> [--patch N] [--dpi D]

<inputdir> is the raw simulation output dir (loadarchive target);
<outdir> is where the figure/csv go. --patch N averages an N x N block
centered on the box center (default 1 = single node) to damp lattice noise.
"""
import matplotlib
matplotlib.use("Agg")

import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from archive.archive import loadarchive
from plot import plot_hd


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("inputdir", help="raw simulation output directory")
    ap.add_argument("outdir", help="directory for the figure / csv")
    ap.add_argument("--patch", type=int, default=1,
                    help="odd block size averaged around the center (default 1)")
    ap.add_argument("--dpi", type=int, default=150)
    return ap.parse_args()


def pressure_field(frame):
    """Mechanical pressure p = -sigma_bulk (or frame.pressure if present)."""
    if hasattr(frame, "pressure"):
        return plot_hd.grid(frame, "pressure")
    return -plot_hd.grid(frame, "sigma_bulk")


def center_pressure(field, patch):
    lx, ly = field.shape
    cx, cy = lx // 2, ly // 2
    if patch <= 1:
        return field[cx, cy]
    h = patch // 2
    block = field[cx - h:cx + h + 1, cy - h:cy + h + 1]
    return float(np.mean(block))


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    ar = loadarchive(args.inputdir)
    n = plot_hd.frame_count(ar)

    ts, ps = [], []
    for i in range(n):
        try:
            frame = ar.read_frame(i)
        except Exception:
            break
        ts.append(plot_hd.frame_time(ar, i))
        ps.append(center_pressure(pressure_field(frame), args.patch))

    ts = np.asarray(ts)
    ps = np.asarray(ps)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(ts, ps, color="#255c99", lw=1.6)
    ax.axhline(0.0, color="0.6", lw=0.8, ls="--")
    ax.set_xlabel("time")
    label = "center pressure" + ("" if args.patch <= 1 else f" ({args.patch}x{args.patch} mean)")
    ax.set_ylabel(r"$-\sigma_{\mathrm{bulk}}$ at center")
    ax.set_title(label)
    fig.tight_layout()

    png = os.path.join(args.outdir, "center_pressure_vs_t.png")
    fig.savefig(png, dpi=args.dpi)
    plt.close(fig)

    csv = os.path.join(args.outdir, "center_pressure_vs_t.csv")
    np.savetxt(csv, np.column_stack([ts, ps]), delimiter=",",
               header="time,center_pressure", comments="")

    print(f"frames={len(ts)}  p_center[last]={ps[-1] if len(ps) else float('nan'):.6g}")
    print(f"wrote {png}")
    print(f"wrote {csv}")


if __name__ == "__main__":
    main()
