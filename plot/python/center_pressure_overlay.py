"""Overlay center-pressure-vs-time curves across an alpha scan.

Reads center_pressure_vs_t.csv from each variant subdirectory of a results
study dir (produced by center_pressure_vs_t.py), parses alpha from the
directory name (encoding 'alpha0p0002' -> 0.0002), and overlays the curves
sorted by alpha. Writes center_pressure_overlay.png into the study dir.

Usage:
    python center_pressure_overlay.py <study_results_dir> [--dpi D]
"""
import matplotlib
matplotlib.use("Agg")

import argparse
import glob
import os
import re

import matplotlib.pyplot as plt
import numpy as np


def alpha_from_name(name):
    """allgrow_dirrandom_alpha0p00015 -> 0.00015 (None if no alpha token)."""
    m = re.search(r"alpha([0-9p]+)", name)
    if not m:
        return None
    try:
        return float(m.group(1).replace("p", "."))
    except ValueError:
        return None


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("studydir", help="results dir holding variant subdirs")
    ap.add_argument("--dpi", type=int, default=150)
    return ap.parse_args()


def main():
    args = parse_args()

    curves = []
    for csv in glob.glob(os.path.join(args.studydir, "*", "center_pressure_vs_t.csv")):
        name = os.path.basename(os.path.dirname(csv))
        alpha = alpha_from_name(name)
        d = np.loadtxt(csv, delimiter=",", skiprows=1)
        if d.ndim == 1:
            d = d.reshape(1, -1)
        curves.append((alpha, name, d[:, 0], d[:, 1]))

    if not curves:
        raise SystemExit(f"no center_pressure_vs_t.csv under {args.studydir}/*/")

    # sort by alpha (None last), large alpha first for legend readability
    curves.sort(key=lambda c: (-1e9 if c[0] is None else -c[0]))

    cmap = plt.cm.viridis
    alphas = [c[0] for c in curves if c[0] is not None]
    amin, amax = (min(alphas), max(alphas)) if alphas else (0.0, 1.0)

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for alpha, name, t, p in curves:
        if alpha is not None and amax > amin:
            color = cmap((alpha - amin) / (amax - amin))
            label = rf"$\alpha={alpha:g}$"
        else:
            color = "0.4"
            label = name
        ax.plot(t, p, color=color, lw=1.7, label=label)

    ax.axhline(0.0, color="0.6", lw=0.8, ls="--")
    ax.set_xlabel("time")
    ax.set_ylabel(r"$-\sigma_{\mathrm{bulk}}$ at center")
    ax.set_title("Center pressure vs time across alpha")
    ax.legend(title="proliferation rate", frameon=False)
    fig.tight_layout()

    out = os.path.join(args.studydir, "center_pressure_overlay.png")
    fig.savefig(out, dpi=args.dpi)
    plt.close(fig)
    print(f"curves={len(curves)}  wrote {out}")


if __name__ == "__main__":
    main()
