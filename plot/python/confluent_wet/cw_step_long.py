#!/usr/bin/env python3
"""Settled or slowly drifting? The 500000-step runs, read as time series.

    cw_step_long.py <long_results>... --calib calib.json --out DIR

Runs unattended behind the array, so the verdict goes into JSON as well as the figure.

WHY TIME SERIES AND NOT A TAIL MEAN. Every other analysis in this campaign reduces a run to
its record-window tail, which is the right summary when the run has settled. These two tau_m
are exactly the cases where that is in doubt: at 3.35 tau_c the three starts landed 0.108
apart in what should be one state, and at 8.32 they agreed at 0.04 one grid point below
where the campaign's operating point collapsed. A tail mean cannot tell a settled state from
a slow slide, so this plots the whole trajectory and then asks whether the last fifth still
moves.

THE TEST. The record window is split into five equal blocks and their means compared. A
settled run has a flat sequence; a drifting one has a monotone one. The block-to-block range
is reported beside the total drift, because a run can wander without going anywhere and that
is a different statement from one that is still on its way somewhere.

The dotted line marks 100000 steps, where the fine scan stopped -- the point of the figure
is what happens to the right of it.
"""
import argparse
import glob
import json
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cw_common as cw                   # noqa: E402
import cw_stream                         # noqa: E402

ARMS = [("a", "chi = 0", "#2166ac"), ("b", "chi = 1", "#b2182b"),
        ("c1", "binary noise", "#7d3c98")]
NBLOCK = 5
SHORT = 100000          # where the fine scan stopped


def series(root, par):
    """(t, <chi>, <m>) from the video stream's per-frame CSV."""
    st = cw_stream.Stream(root, par)
    t = np.asarray(st.steps, float)[:st.n]
    chi = np.asarray(st.meta["chi_mean"], float)[:st.n]
    m = np.asarray(st.meta["m_mean"], float)[:st.n]
    return t, chi, m


def blocks(t, y, t0):
    """Means of NBLOCK equal blocks of the record window, oldest first."""
    k = t >= t0
    tt, yy = t[k], y[k]
    if len(tt) < NBLOCK * 2:
        return []
    edges = np.linspace(tt[0], tt[-1], NBLOCK + 1)
    out = []
    for i in range(NBLOCK):
        sel = (tt >= edges[i]) & (tt <= edges[i + 1])
        out.append(float(yy[sel].mean()) if sel.any() else float("nan"))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("long", nargs="+",
                    help="one or more results trees; they are pooled and each tau_m gets "
                         "its own column")
    ap.add_argument("--calib", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--record-frac", type=float, default=0.2)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    with open(a.calib) as fh:
        cal = json.load(fh)

    runs = {}
    for p in sorted(sum((glob.glob(os.path.join(d, "*", "part.json")) for d in a.long), [])):
        d = os.path.dirname(p)
        c = os.path.basename(d)
        m = re.match(r"tm(.+)_([abc]\d?)$", c)
        if not m:
            continue
        with open(p) as fh:
            part = json.load(fh)
        src = part.get("inputdir") or d
        try:
            par = cw.read_params(src)
            t, chi, mm = series(src, par)
        except Exception as exc:
            print(f"  {c}: no stream ({type(exc).__name__}: {exc})")
            continue
        runs[(float(m.group(1).replace("p", ".")), m.group(2))] = (t, chi, mm, part)

    if not runs:
        raise RuntimeError(f"no usable runs under {' '.join(a.long)}")
    gs = sorted(set(g for g, _ in runs))
    tau_c = cal["tau_c"]

    fig, ax = plt.subplots(2, len(gs), figsize=(6.0 * len(gs), 7.2), squeeze=False,
                           constrained_layout=True)
    summary = {}
    for j, g in enumerate(gs):
        tm = None
        for arm, lab, col in ARMS:
            if (g, arm) not in runs:
                continue
            t, chi, mm, part = runs[(g, arm)]
            tm = part["params"]["tau_m"]
            t0 = a.record_frac * part["params"]["nsteps"]
            ax[0, j].plot(t / tm, chi, "-", color=col, lw=1.0, label=lab)
            ax[1, j].plot(t / tm, mm, "-", color=col, lw=1.0, label=lab)
            bl = blocks(t, chi, t0)
            summary[f"tm{g:g}_{arm}"] = {
                "tau_m": tm, "nsteps": part["params"]["nsteps"],
                "blocks": bl,
                "block_range": (max(bl) - min(bl)) if bl else None,
                "net_drift": (bl[-1] - bl[0]) if bl else None,
                "tail": bl[-1] if bl else None,
                "settled": bool(bl and abs(bl[-1] - bl[0]) < 0.02
                                and (max(bl) - min(bl)) < 0.05),
            }
        for i in (0, 1):
            ax[i, j].axvline(SHORT / tm, color="k", ls=":", lw=1.3)
            ax[i, j].grid(alpha=0.25)
            ax[i, j].set_xlabel(r"$t/\tau_m$")
        ax[0, j].axhline(cal["mc"], color="0.6", ls="--", lw=0.8)
        ax[1, j].axhline(cal["mc"], color="0.6", ls="--", lw=0.8)
        ax[0, j].set_ylim(-0.05, 1.05)
        ax[0, j].set_title(rf"$\tau_m$ = {g:g} $\tau_c$ = {tm:.0f} steps", fontsize=10)
        ax[0, j].legend(fontsize=8)
    ax[0, 0].set_ylabel(r"$\langle\chi\rangle$")
    ax[1, 0].set_ylabel(r"$\langle m\rangle$  (dashed: $m_c$)")
    fig.suptitle("500000 steps at the symmetric point -- dotted line is where the fine "
                 "scan stopped", fontsize=10)
    fig.savefig(os.path.join(a.out, "long.png"), dpi=200)
    plt.close(fig)

    with open(os.path.join(a.out, "long.json"), "w") as fh:
        json.dump({"mc": cal["mc"], "tau_c": tau_c, "record_frac": a.record_frac,
                   "runs": summary}, fh, indent=1)

    print(f"{len(runs)} runs")
    print(f"{'case':>12} {'blocks over the record window':>44} {'range':>7} {'drift':>7} "
          f"{'settled':>8}")
    for k in sorted(summary):
        v = summary[k]
        bl = " ".join(f"{b:.3f}" for b in v["blocks"]) if v["blocks"] else "-"
        print(f"{k:>12} {bl:>44} {v['block_range'] or float('nan'):7.3f} "
              f"{v['net_drift'] or float('nan'):+7.3f} {str(v['settled']):>8}")
    print(f"\nwrote {a.out}/long.png and long.json")


if __name__ == "__main__":
    main()
