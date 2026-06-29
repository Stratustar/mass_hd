#!/usr/bin/env python3
"""Pressure field + bulk-stress contribution breakdown, one figure per variant.

    P = -sigma_bulk

The bulk stress splits into (go-or-grow.cpp UpdateQuantitiesAtNode):

    sigma_bulk = compress + LdG-nematic + [surface_stress] * surface

    compress = 0.5*B*d+^2 - B*d+*phi ,  d+ = max(phi - phi_critical, 0)
    LdG      = 0.5*CC*term^2 - CC*Snem*term*phi ,  term = Snem*phi - |Q|^2
    surface  = f_surf - mu_surf*phi              (only if surface_stress != 0)

compress needs only phi, so it is computed exactly here.  When surface_stress=0
(the current batch) the rest is pure LdG and is obtained by subtraction --
LdG = sigma_bulk - compress -- with no need to recompute term or any laplacian.
If surface_stress != 0 the remainder is LdG+surface lumped; a warning is printed.

Each variant gets one figure: pressure map (coolwarm) on the left, the three
contributions' shares (mean |.| over the colony, phi>phi_th) as a vertical bar
list on the right.  A summary CSV collects the per-variant percentages.

Usage:
    python pressure_decomposition.py <study_dir> <outdir> [--pattern '*'] [--phi-th 0.5]
"""
import matplotlib
matplotlib.use("Agg")

import argparse
import csv
import glob
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)
from archive.archive import loadarchive

COMPONENTS = ["compress", "LdG-nematic", "surface"]
COLORS = {"compress": "#3a7d44", "LdG-nematic": "#c0653a", "surface": "#888888"}


def grid(frame, field):
    lx = frame.parameters["LX"]
    ly = frame.parameters["LY"]
    return np.array(getattr(frame, field)).reshape((lx, ly))


def decompose(frame):
    """Return the three pressure contributions (P = -sigma_bulk) as fields."""
    p = frame.parameters
    B = float(p["B"])
    phic = float(p["phi_critical"])
    surf_on = float(p["surface_stress"])

    phi = grid(frame, "phi")
    sigma_bulk = grid(frame, "sigma_bulk")

    dplus = np.clip(phi - phic, 0.0, None)
    s_compress = 0.5 * B * dplus ** 2 - B * dplus * phi          # compress stress
    remainder = sigma_bulk - s_compress                          # LdG (+surface if on)

    if surf_on != 0.0:
        print("WARNING: surface_stress=%g -> 'LdG-nematic' share also contains "
              "the surface term (not split)." % surf_on)

    fields = {
        "compress": -s_compress,        # pressure contributions: P_i = -sigma_i
        "LdG-nematic": -remainder,
        "surface": np.zeros_like(phi),  # identically 0 when surface_stress=0
    }
    return phi, -sigma_bulk, fields, surf_on


def shares(phi, fields, phi_th):
    """Percentage share of each contribution = mean |P_i| over the colony."""
    mask = phi > phi_th
    if not mask.any():
        return {k: 0.0 for k in COMPONENTS}
    mag = {k: float(np.mean(np.abs(fields[k][mask]))) for k in COMPONENTS}
    tot = sum(mag.values()) or 1.0
    return {k: 100.0 * mag[k] / tot for k in COMPONENTS}


def plot_variant(name, phi, P, fields, sh, phi_th, outpath):
    mask = phi > phi_th
    fig, (axP, axT) = plt.subplots(
        1, 2, figsize=(7.6, 4.0), gridspec_kw=dict(width_ratios=[3, 1.1]))

    vmax = np.percentile(np.abs(P[mask]), 99) if mask.any() else float(np.abs(P).max())
    vmax = vmax or 1.0
    im = axP.imshow(P.T, origin="lower", cmap="coolwarm", vmin=-vmax, vmax=vmax,
                    interpolation="nearest")
    axP.contour(phi.T, levels=[phi_th], colors="k", linewidths=0.6)
    axP.set_title("pressure  P = -sigma_bulk")
    axP.set_xticks([]); axP.set_yticks([])
    cb = fig.colorbar(im, ax=axP, fraction=0.046, pad=0.04)
    cb.set_label("P")

    # three plain-text rows of the share percentages (no extra plot)
    axT.axis("off")
    labels = {"compress": "compress", "LdG-nematic": "LdG", "surface": "surface"}
    ypos = [0.62, 0.50, 0.38]
    for k, yy in zip(COMPONENTS, ypos):
        axT.text(0.0, yy, "%-9s %5.1f%%" % (labels[k] + ":", sh[k]),
                 transform=axT.transAxes, va="center", ha="left",
                 fontsize=12, family="monospace", color=COLORS[k])

    fig.suptitle(name, fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("studydir")
    ap.add_argument("outdir")
    ap.add_argument("--pattern", default="*")
    ap.add_argument("--phi-th", type=float, default=0.5)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    rows = []
    for dpath in sorted(glob.glob(os.path.join(args.studydir, args.pattern))):
        if not os.path.isdir(dpath):
            continue
        name = os.path.basename(dpath)
        try:
            ar = loadarchive(dpath)
            last = int((ar.nsteps - ar.nstart) / ar.ninfo)
            frame = ar.read_frame(last)
            phi, P, fields, surf_on = decompose(frame)
            sh = shares(phi, fields, args.phi_th)
            out = os.path.join(args.outdir, name + "_pressure_decomp.png")
            plot_variant(name, phi, P, fields, sh, args.phi_th, out)
            rows.append(dict(variant=name, surface_stress=surf_on,
                             compress_pct=sh["compress"], ldg_pct=sh["LdG-nematic"],
                             surface_pct=sh["surface"]))
            print("ok  %-40s compress=%5.1f%%  LdG=%5.1f%%  surface=%4.1f%%"
                  % (name, sh["compress"], sh["LdG-nematic"], sh["surface"]))
        except Exception as e:
            print("FAIL", name, repr(e))

    with open(os.path.join(args.outdir, "pressure_decomp_summary.csv"), "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=["variant", "surface_stress",
                                            "compress_pct", "ldg_pct", "surface_pct"])
        wr.writeheader()
        for r in rows:
            wr.writerow(r)
    print("wrote %d figures + summary to %s" % (len(rows), args.outdir))


if __name__ == "__main__":
    main()
