#!/usr/bin/env python3
"""Full-domain chi field per variant, to check for spurious chi outside the colony.

The standard `visualization` masks chi to phi>0.5 and crops to the colony, hiding
whatever chi does in the low-density halo. Here chi is rendered over the WHOLE
domain, uncropped and unmasked, with the phi=0.5 (colony boundary) and phi=0.01
(near-vacuum edge) contours overlaid. Between them is the halo, where chi = m/phi
of two small tails can go spurious.

Recall the code: chi = 0 where phi <= 1e-12 (deep vacuum), else chi = m/phi (with
m clamped to <= phi, so chi in [0,1]). "Spurious" chi therefore shows up as [0,1]
structure in the halo ring just outside the phi=0.5 contour.

Outputs, per variant, `chi_full_frame<n>.png` (last frame by default) and a study
CSV `chi_full_halo_stats.csv` quantifying chi in the halo (1e-12 < phi < 0.5):
max, mean, and the fraction of halo area with |chi-0.5| > 0.1 (a spurious flag).

Usage: python chi_full_field.py <study_dir> <outdir> [--pattern '*'] [--frames last|all|N ...]
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

HALO_LO = 1e-12       # below this phi, chi is forced to 0 (deep vacuum)
COLONY = 0.5          # official colony boundary


def grid(frame, field):
    lx = frame.parameters["LX"]
    ly = frame.parameters["LY"]
    return np.array(getattr(frame, field)).reshape((lx, ly))


NEARVAC = 0.05        # below this phi there is essentially no material: chi should be ~0

def halo_stats(chi, phi):
    halo = (phi > HALO_LO) & (phi < COLONY)
    nearvac = (phi > HALO_LO) & (phi < NEARVAC)
    if not halo.any():
        return dict(halo_px=0, nearvac_px=0, chi_max=0.0, chi_mean=0.0, spurious_frac=0.0)
    c = chi[halo]
    # spurious = near-vacuum (phi<0.05, no material) cells that still carry chi>0.1
    sp = float(np.mean(chi[nearvac] > 0.1)) if nearvac.any() else 0.0
    return dict(
        halo_px=int(halo.sum()),
        nearvac_px=int(nearvac.sum()),
        chi_max=float(np.max(c)),
        chi_mean=float(np.mean(c)),
        spurious_frac=sp,
    )


def plot_chi_full(name, chi, phi, stats, outpath):
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    im = ax.imshow(chi.T, origin="lower", cmap="bwr", vmin=0.0, vmax=1.0,
                   interpolation="nearest")
    # colony boundary (solid black) and near-vacuum halo edge (dashed grey)
    ax.contour(phi.T, levels=[COLONY], colors="k", linewidths=0.8, origin="lower")
    ax.contour(phi.T, levels=[0.01], colors="0.4", linewidths=0.5,
               linestyles="dashed", origin="lower")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("chi (full domain)\nhalo max=%.2f mean=%.2f | near-vac(phi<0.05) chi>0.1: %.1f%%"
                 % (stats["chi_max"], stats["chi_mean"], 100 * stats["spurious_frac"]),
                 fontsize=8)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04); cb.set_label("chi")
    fig.suptitle(name, fontsize=8)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def resolve_frames(ar, spec):
    last = int((ar.nsteps - ar.nstart) / ar.ninfo)
    out = []
    for s in spec:
        if s == "last":
            out.append(last)
        elif s == "all":
            out.extend(range(0, last + 1))
        else:
            out.append(min(int(s), last))
    return sorted(set(out))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("studydir")
    ap.add_argument("outdir")
    ap.add_argument("--pattern", default="*")
    ap.add_argument("--frames", nargs="+", default=["last"])
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    rows = []
    for dpath in sorted(glob.glob(os.path.join(args.studydir, args.pattern))):
        if not os.path.isdir(dpath):
            continue
        name = os.path.basename(dpath)
        try:
            ar = loadarchive(dpath)
            for fi in resolve_frames(ar, args.frames):
                frame = ar.read_frame(fi)
                chi, phi = grid(frame, "chi"), grid(frame, "phi")
                st = halo_stats(chi, phi)
                vout = os.path.join(args.outdir, name)
                os.makedirs(vout, exist_ok=True)
                plot_chi_full(name, chi, phi, st,
                              os.path.join(vout, "chi_full_frame%d.png" % fi))
                if fi == resolve_frames(ar, ["last"])[0]:
                    rows.append(dict(variant=name, frame=fi, **st))
            print("ok  %-44s halo max=%.2f mean=%.2f spurious=%.1f%%"
                  % (name, rows[-1]["chi_max"], rows[-1]["chi_mean"],
                     100 * rows[-1]["spurious_frac"]))
        except Exception as e:
            print("FAIL", name, repr(e))

    csv_path = os.path.join(args.outdir, "chi_full_halo_stats.csv")
    with open(csv_path, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=["variant", "frame", "halo_px", "nearvac_px",
                                            "chi_max", "chi_mean", "spurious_frac"])
        wr.writeheader()
        for r in rows:
            wr.writerow(r)
    # headline: worst spurious cases
    worst = sorted(rows, key=lambda r: r["spurious_frac"], reverse=True)[:5]
    print("\n=== most spurious halo chi (top 5) ===")
    for r in worst:
        print("  %-44s spurious=%.1f%%  max=%.2f" % (r["variant"], 100 * r["spurious_frac"], r["chi_max"]))
    print("wrote %d figures + %s" % (len(rows), csv_path))


if __name__ == "__main__":
    main()
