#!/usr/bin/env python3
"""The campaign-level figures of the STEP-SWITCH tau_m scan.

    cw_step_fig.py <b1_results> <b2_results> --calib calib.json --out DIR

Two figures, and between them they are the whole result:

  fates.png   a 2-panel map (B1 | B2). Rows are the four initial conditions, the x axis is
              tau_m/tau_c on a log scale, and each cell is coloured by the record-window
              tail <chi>. Cells whose run is still drifting are hatched -- an UNDECIDED run
              must not be read as a state, and hatching says so without hiding the number.
              tau_x is drawn as a vertical line on both panels: left of it the four rows
              should agree, right of it a and b should separate and c has to pick a side.

  curves.png  the same data as curves, <chi>(tau_m) per start, B1 solid and B2 dashed. The
              map answers "which phase"; this answers "how sharply", and it is where a
              boundary that moved between the two tau_chi rules would be visible.

WHY THE TAIL AND NOT THE RECORD MEAN: a run sliding from one phase to the other has a
record-window mean in the middle, and plotting that invents a third state out of a
transient. Every number here is fate.chi_tail with fate.drifting carried alongside it.
"""
import argparse
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402

STARTS = [("a", "chi = 0\n(active)"), ("b", "chi = 1\n(passive)"),
          ("c1", "noise\nseed 1"), ("c2", "noise\nseed 2")]


def load(d):
    """{(tau_m/tau_c, start): part} from a results tree of tm<g>_<start>/part.json."""
    out = {}
    for p in sorted(glob.glob(os.path.join(d, "*", "part.json"))):
        case = os.path.basename(os.path.dirname(p))
        if not case.startswith("tm") or "_" not in case:
            continue
        gtag, start = case[2:].split("_", 1)
        with open(p) as fh:
            j = json.load(fh)
        out[(float(gtag.replace("p", ".")), start)] = j
    return out


def cell(j):
    """(<chi> tail, drifting) or (nan, False) when the run has no usable fate."""
    f = j.get("fate")
    if not f:
        s = j.get("settled")
        return (float(s["chi_mean_tail"]), False) if s else (float("nan"), False)
    return float(f["chi_tail"]), bool(f["drifting"])


def panel(ax, data, gs, title, tau_x):
    grid = np.full((len(STARTS), len(gs)), np.nan)
    drift = np.zeros_like(grid, dtype=bool)
    for i, (key, _) in enumerate(STARTS):
        for k, g in enumerate(gs):
            if (g, key) in data:
                grid[i, k], drift[i, k] = cell(data[(g, key)])

    x = np.arange(len(gs))
    im = ax.pcolormesh(np.arange(len(gs) + 1) - 0.5, np.arange(len(STARTS) + 1) - 0.5,
                       grid, cmap="RdYlBu", vmin=0.0, vmax=1.0, shading="flat")
    for i in range(len(STARTS)):
        for k in range(len(gs)):
            if drift[i, k]:
                ax.add_patch(plt.Rectangle((k - 0.5, i - 0.5), 1, 1, fill=False,
                                           hatch="///", edgecolor="0.15", linewidth=0.0))
            if np.isfinite(grid[i, k]):
                ax.text(k, i, f"{grid[i, k]:.2f}", ha="center", va="center", fontsize=6,
                        color="0.1" if 0.25 < grid[i, k] < 0.75 else "w")
    if np.isfinite(tau_x) and gs[0] < tau_x < gs[-1]:
        xi = float(np.interp(np.log(tau_x), np.log(gs), x))
        ax.axvline(xi, color="k", lw=1.6, ls=":")
        # axis-fraction y so the label sits ABOVE the panel: the data axis is inverted and
        # the bottom of it is where the x tick labels are, which is what it collided with.
        ax.text(xi, 1.02, r"$\tau_x$", transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{g:g}" for g in gs], fontsize=7, rotation=90)
    ax.set_yticks(range(len(STARTS)))
    ax.set_yticklabels([lab for _, lab in STARTS], fontsize=7)
    ax.set_xlabel(r"$\tau_m/\tau_c$")
    ax.set_title(title, fontsize=9)
    ax.invert_yaxis()
    return im


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("b1")
    ap.add_argument("b2")
    ap.add_argument("--calib", required=True)
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    with open(a.calib) as fh:
        cal = json.load(fh)
    tau_x = float(cal.get("tau_x", float("nan")))

    d1, d2 = load(a.b1), load(a.b2)
    if not d1 and not d2:
        raise RuntimeError("no part.json found under either results tree")
    g1 = sorted(set(g for g, _ in d1))
    g2 = sorted(set(g for g, _ in d2))

    # ---------------------------------------------------------------- the map
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.4), constrained_layout=True)
    im = None
    for ax, data, gs, ttl in ((axes[0], d1, g1, r"B1:  $\tau_\chi = 0.3\,\tau_c$"),
                              (axes[1], d2, g2, r"B2:  $\tau_\chi = \tau_m$")):
        if gs:
            im = panel(ax, data, gs, ttl, tau_x)
        else:
            ax.set_title(ttl + "  (no runs)", fontsize=9); ax.axis("off")
    if im is not None:
        cb = fig.colorbar(im, ax=axes, fraction=0.03, pad=0.01)
        cb.set_label(r"record-window tail $\langle\chi\rangle$   (0 = active, 1 = passive)",
                     fontsize=8)
        cb.ax.tick_params(labelsize=7)
    fig.suptitle(f"step switch, floor 0.3 zeta -- fate against memory time "
                 f"(hatched = still drifting; tau_x = {tau_x:.2f} tau_c)", fontsize=9)
    fig.savefig(os.path.join(a.out, "fates.png"), dpi=200)
    plt.close(fig)

    # ------------------------------------------------------------- the curves
    fig, ax = plt.subplots(figsize=(6.4, 4.0), constrained_layout=True)
    colours = {"a": "#c0392b", "b": "#2471a3", "c1": "#7d3c98", "c2": "#b7950b"}
    for data, gs, ls, tag in ((d1, g1, "-", "B1"), (d2, g2, "--", "B2")):
        for key, lab in STARTS:
            xs = [g for g in gs if (g, key) in data]
            if not xs:
                continue
            ys = [cell(data[(g, key)])[0] for g in xs]
            dr = [cell(data[(g, key)])[1] for g in xs]
            ax.plot(xs, ys, ls, color=colours[key], lw=1.3, marker="o", ms=3.2,
                    label=f"{tag} {key}")
            xd = [x for x, d in zip(xs, dr) if d]
            yd = [y for y, d in zip(ys, dr) if d]
            if xd:
                ax.plot(xd, yd, "o", ms=7, mfc="none", mec="0.2", mew=1.0)
    if np.isfinite(tau_x):
        ax.axvline(tau_x, color="k", ls=":", lw=1.6)
        ax.text(tau_x, 1.02, r"$\tau_x$", transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=9)
    ax.set_xscale("log")
    ax.set_xlabel(r"$\tau_m/\tau_c$")
    ax.set_ylabel(r"tail $\langle\chi\rangle$")
    ax.set_ylim(-0.05, 1.05)
    ax.grid(alpha=0.25)
    # centre, not an edge: once the runs separate the middle of the panel is the one band
    # with no data in it, while every edge has either the active or the passive branch on it
    ax.legend(fontsize=6.5, ncol=2, loc="center", framealpha=0.92)
    ax.set_title("open circles = still drifting (undecided)", fontsize=8)
    fig.savefig(os.path.join(a.out, "curves.png"), dpi=200)
    plt.close(fig)
    print(f"wrote {a.out}/fates.png and {a.out}/curves.png "
          f"({len(d1)} B1 runs, {len(d2)} B2 runs)")


if __name__ == "__main__":
    main()
