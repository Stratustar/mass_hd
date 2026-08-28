"""Show the L_chi/d collapse instead of asserting it.

At a matched ratio tau_m/t_eddy ~ 2 the three activities have raw phenotype correlation
lengths of 37.6 / 25.1 / 14.8 -- a factor 2.5 apart, which is why the raw fields look
different. Divided by each run's own defect spacing they are 0.73 / 0.78 / 0.85. So the claim
is that the fields are statistically the SAME PICTURE at different magnifications, and the
honest way to show that is to crop each one to the same number of defect spacings and blow
them up to a common size. Top row: as simulated. Bottom row: same fields, N*d on a side.

Usage: cw_montage.py <out.png> <ndefect_spacings> <case> [<case> ...]
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cw_common as cw

STUDY = "20260828/cw_loop"


def main():
    out, nd = sys.argv[1], float(sys.argv[2])
    cases = sys.argv[3:]
    scan = {}
    import json
    with open(cw.results_root(*STUDY.split("/"), "cw_scan.json")) as f:
        for r in json.load(f):
            scan[r["case"]] = r

    fig, axes = plt.subplots(2, len(cases), figsize=(4.2 * len(cases), 8.8), dpi=125)
    for j, c in enumerate(cases):
        root = cw.case_root(*STUDY.split("/"), c)
        oa = cw.loadarchive(root)
        fr = cw.load_frame(oa, cw.frame_count(root) - 1)
        chi = fr["chi"]
        r = scan[c]
        d, L = r["d"], chi.shape[0]

        ax = axes[0, j]
        ax.imshow(chi.T, origin="lower", cmap="coolwarm", vmin=0, vmax=1,
                  interpolation="nearest")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"A = {r['A']:.0f}    " + r"$\tau_m/t_{\rm eddy}$ = "
                     f"{r['tau_m_set']/r['t_eddy']:.2f}\n"
                     f"d = {d:.1f}   " + r"$L_\chi$ = " + f"{r['L_chi']:.1f}   "
                     + r"$L_\chi/d$ = " + f"{r['L_chi']/d:.2f}", fontsize=10)
        # the box that the bottom row shows
        w = nd * d
        ax.add_patch(plt.Rectangle(((L - w) / 2, (L - w) / 2), w, w,
                                   fill=False, ec="k", lw=1.6, ls="--"))

        ax = axes[1, j]
        i0 = int((L - w) / 2)
        i1 = int(i0 + w)
        ax.imshow(chi[i0:i1, i0:i1].T, origin="lower", cmap="coolwarm", vmin=0, vmax=1,
                  interpolation="bilinear")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{nd:.0f}d x {nd:.0f}d  =  {w:.0f} x {w:.0f} cells", fontsize=10)

    axes[0, 0].set_ylabel("as simulated\n(400 x 400 lattice)", fontsize=11)
    axes[1, 0].set_ylabel(f"cropped to {nd:.0f} defect spacings\nand rescaled to a common size",
                          fontsize=11)
    fig.suptitle(r"The phenotype field $\chi$ at a matched $\tau_m/t_{\rm eddy}$:  "
                 "different at fixed lattice scale, the same at fixed $d$", fontsize=13)
    fig.tight_layout()
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
