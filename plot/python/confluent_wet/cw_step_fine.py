#!/usr/bin/env python3
"""The shape of the transition: a fine tau_m scan at the symmetric operating point.

    cw_step_fine.py <fine_results_dir> --calib calib.json --out DIR

Runs unattended behind the simulation array, so it says what it found in a JSON summary as
well as in the figure -- nobody is watching the log.

WHAT IT MEASURES. Four panels against log tau_m:

  1. tail <chi> for the three starts. Below the coexistence threshold they must agree; the
     question is whether the single state slides continuously or jumps.
  2. the stability margin |<m> - mc| / std(m), computed in each run's OWN std(m). This is
     the quantity the symmetric operating point was chosen to equalise, and it is what
     decides whether a phase survives: the campaign's active phase sat at 2.55 and decayed,
     the symmetric one sits at 5.27 and does not.
  3. std(m) against the one-pole closed form sqrt(f(1-f)) sqrt(tau_c/(tau_c+tau_m)). The
     campaign measured the ratio falling from 0.72 to 0.23 across 1-30 tau_c because Dbio
     smooths m over sqrt(2 Dbio tau_m); this fills in that curve at 20 points.
  4. <m> itself against mc, which is panel 2 without the normalisation -- it separates "the
     memory moved" from "the memory got quieter".

SEPARATION is reported as max over starts minus min over starts at each tau_m, and the
threshold as the first grid point where it exceeds 0.5 and stays there.
"""
import argparse
import glob
import json
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402

STARTS = [("a", "χ ≡ 0", "#2166ac"), ("b", "χ ≡ 1", "#b2182b"),
          ("c1", "二值噪声", "#7d3c98")]
FATE_MARK = {"active": "o", "passive": "s", "mixed": "^", "undecided": "D"}


def load(d):
    out = {}
    for p in sorted(glob.glob(os.path.join(d, "*", "part.json"))):
        c = os.path.basename(os.path.dirname(p))
        m = re.match(r"tm(.+)_([abc]\d?)$", c)
        if not m:
            continue
        with open(p) as fh:
            out[(float(m.group(1).replace("p", ".")), m.group(2))] = json.load(fh)
    return out


def row(d):
    f = d.get("fate") or {}
    s = d.get("settled") or {}
    u = d["params"]
    fl = d["flow"]
    m = s.get("m_mean_tail", float("nan"))
    sd = fl["std_m"]
    dep = f.get("departure_steps", float("nan"))
    lag = d.get("lagrangian") or {}
    f_here = lag.get("f_lag", float("nan"))
    closed = (np.sqrt(f_here * (1 - f_here)) * np.sqrt(1.0 / (1.0 + u["tau_m"] / 225.56))
              if f_here == f_here else float("nan"))
    return dict(chi=f.get("chi_tail", s.get("chi_mean_tail", float("nan"))),
                label=f.get("label", "?"), drift=bool(f.get("drifting", False)),
                m=m, sd=sd, mc=u["mc"], margin=abs(m - u["mc"]) / sd if sd > 0 else float("nan"),
                std_chi=fl["std_chi"], life=dep / u["tau_m"] if dep == dep else float("inf"),
                f_lag=f_here, closed=closed, tau_m=u["tau_m"], N_def=fl["N_def"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fine")
    ap.add_argument("--calib", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    with open(a.calib) as fh:
        cal = json.load(fh)

    data = load(a.fine)
    if not data:
        raise RuntimeError(f"no part.json under {a.fine}")
    gs = sorted(set(g for g, _ in data))
    R = {(g, k): row(d) for (g, k), d in data.items()}
    mc = next(iter(R.values()))["mc"]

    # ---- separation and where it opens
    sep, thr = [], None
    for g in gs:
        vals = [R[(g, k)]["chi"] for k, _, _ in STARTS if (g, k) in R]
        sep.append(max(vals) - min(vals) if len(vals) > 1 else float("nan"))
    for i, g in enumerate(gs):
        if sep[i] > 0.5 and all(s > 0.5 for s in sep[i:] if s == s):
            thr = g
            break

    fig, ax = plt.subplots(2, 2, figsize=(11.5, 7.6), constrained_layout=True)

    # 1. the fates
    for k, lab, col in STARTS:
        x = [g for g in gs if (g, k) in R]
        y = [R[(g, k)]["chi"] for g in x]
        ax[0, 0].plot(x, y, "-", color=col, lw=1.4, label=lab, zorder=2)
        for g, v in zip(x, y):
            r = R[(g, k)]
            ax[0, 0].plot([g], [v], FATE_MARK.get(r["label"], "o"), color=col, ms=4.5,
                          mfc="none" if r["drift"] else col, zorder=3)
    if thr:
        ax[0, 0].axvline(thr, color="k", ls=":", lw=1.4)
        ax[0, 0].text(thr, 1.02, f"分离 {thr:g}", transform=ax[0, 0].get_xaxis_transform(),
                      ha="center", va="bottom", fontsize=8)
    ax[0, 0].set_ylabel(r"tail $\langle\chi\rangle$")
    ax[0, 0].set_ylim(-0.05, 1.05)
    ax[0, 0].legend(fontsize=8, loc="center left")
    ax[0, 0].set_title(f"命运（空心 = 仍在漂移）　mc = {mc:g}", fontsize=9)

    # 2. the margin -- the quantity the symmetric point equalises
    for k, lab, col in STARTS:
        x = [g for g in gs if (g, k) in R]
        ax[0, 1].plot(x, [R[(g, k)]["margin"] for g in x], "-o", color=col, lw=1.4, ms=3.5,
                      label=lab)
    ax[0, 1].axhline(1.0, color="0.5", ls="--", lw=1)
    ax[0, 1].set_ylabel(r"$|\langle m\rangle - m_c| / \mathrm{std}(m)$")
    ax[0, 1].set_title("稳定性余量（该相自己的 σ_m）", fontsize=9)
    ax[0, 1].legend(fontsize=8)

    # 3. std(m) against the one-pole closed form
    for k, lab, col in STARTS:
        x = [g for g in gs if (g, k) in R]
        ax[1, 0].plot(x, [R[(g, k)]["sd"] for g in x], "-o", color=col, lw=1.4, ms=3.5,
                      label=f"{lab} 实测")
    x0 = [g for g in gs if (g, "a") in R]
    cl = [R[(g, "a")]["closed"] for g in x0]
    if any(v == v for v in cl):
        ax[1, 0].plot(x0, cl, "--", color="0.45", lw=1.2, label="单极点闭式 (χ≡0 的 f)")
    ax[1, 0].set_yscale("log")
    ax[1, 0].set_ylabel(r"$\mathrm{std}(m)$")
    ax[1, 0].set_xlabel(r"$\tau_m/\tau_c$")
    ax[1, 0].set_title("记忆涨落：实测 vs 闭式（差值即 Dbio 的平滑）", fontsize=9)
    ax[1, 0].legend(fontsize=7)

    # 4. <m> against the threshold
    for k, lab, col in STARTS:
        x = [g for g in gs if (g, k) in R]
        ax[1, 1].plot(x, [R[(g, k)]["m"] for g in x], "-o", color=col, lw=1.4, ms=3.5,
                      label=lab)
    ax[1, 1].axhline(mc, color="k", ls=":", lw=1.4)
    ax[1, 1].text(gs[0], mc, " mc", va="bottom", fontsize=8)
    ax[1, 1].axhline(cal["f_top"], color="0.6", ls="--", lw=1)
    ax[1, 1].axhline(cal["f_floor"], color="0.6", ls="--", lw=1)
    ax[1, 1].set_ylabel(r"tail $\langle m\rangle$")
    ax[1, 1].set_xlabel(r"$\tau_m/\tau_c$")
    ax[1, 1].set_title("记忆本身（虚线 = f_floor, f_top）", fontsize=9)
    ax[1, 1].legend(fontsize=8)

    for axx in ax.ravel():
        axx.set_xscale("log")
        axx.grid(alpha=0.25)
    fig.suptitle("阶跃开关 · 对称工作点 mc = 0.21 · τ_m 细扫描 0.3–10 τ_c", fontsize=10)
    fig.savefig(os.path.join(a.out, "fine.png"), dpi=200)
    plt.close(fig)

    summary = {
        "mc": mc, "pmem": next(iter(data.values()))["params"]["pmem"],
        "tau_c": cal["tau_c"], "n_runs": len(data), "tau_m_grid": gs,
        "separation": [float(s) for s in sep],
        "separation_threshold": thr,
        "per_run": {f"tm{g:g}_{k}": {kk: (None if isinstance(vv, float) and vv != vv else vv)
                                     for kk, vv in R[(g, k)].items()}
                    for (g, k) in sorted(R)},
    }
    with open(os.path.join(a.out, "fine.json"), "w") as fh:
        json.dump(summary, fh, indent=1)

    print(f"{len(data)} runs, {len(gs)} tau_m points")
    print(f"{'tau_m/tau_c':>11} " + " ".join(f"{lab:>12}" for _, lab, _ in STARTS)
          + f" {'separation':>11}")
    for g, s in zip(gs, sep):
        cells = []
        for k, _, _ in STARTS:
            r = R.get((g, k))
            cells.append(f'{r["chi"]:.3f}{"*" if r["drift"] else " "}{r["label"][:4]:>7}'
                         if r else f'{"-":>12}')
        print(f"{g:11g} " + " ".join(f"{c:>12}" for c in cells) + f" {s:11.3f}")
    print(f"\nseparation opens at tau_m/tau_c = {thr if thr else 'nowhere in this range'}")
    print(f"wrote {a.out}/fine.png and fine.json")


if __name__ == "__main__":
    main()
