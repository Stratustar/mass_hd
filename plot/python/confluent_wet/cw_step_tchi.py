#!/usr/bin/env python3
"""Does the phenotype clock move the bifurcation? Three tau_chi rules on one axis.

    cw_step_tchi.py <tchi_results> <fine_results> <sym_results> --calib c.json --out DIR

Runs unattended behind the simulation array, so the finding goes into a JSON summary as well
as the figure.

THE QUESTION. The campaign predicted the fate depends on sigma_m alone and therefore not on
tau_chi. It does depend on it: at mc = 0.2610 the tau_chi = tau_m group held the active phase
at 22 tau_c where tau_chi = 0.3 tau_c did not. The proposed mechanism is that a slow chi
low-passes its response to m's fluctuations and damps the feedback loop. This puts three
rules -- 0.3 tau_c and 1.0 tau_c fixed, 0.5 tau_m proportional -- on the same five tau_m at
the symmetric operating point, which is the cleanest test available: if the fate is a
function of sigma_m, the three collapse.

The 0.3 tau_c rule comes from the existing fine and sym trees except at tau_m = 3, which the
tchi tree supplies, because the fine grid has 2.79 and 3.35 and an interpolated point is not
worth having when the comparison IS the purpose.
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

TAU_M = [0.3, 3, 10, 15, 30]
# label, colour, and how tau_chi is set
RULES = [("0.3 tau_c", "#2166ac"), ("1.0 tau_c", "#b26a00"), ("0.5 tau_m", "#7d3c98")]
ARMS = [("a", "chi = 0", "-", "o"), ("b", "chi = 1", "--", "s"),
        ("c1", "noise", ":", "^")]
TOL = 0.06          # how close a tau_m has to be to count as one of the five


def scan(d, pat):
    out = []
    for p in sorted(glob.glob(os.path.join(d, "*", "part.json"))):
        c = os.path.basename(os.path.dirname(p))
        m = re.match(pat, c)
        if not m:
            continue
        with open(p) as fh:
            out.append((c, json.load(fh)))
    return out


def rule_of(u, tau_c):
    """Which of the three rules a run belongs to, from its own numbers."""
    gchi = u["tau_chi"] / tau_c
    frac = u["tau_chi"] / u["tau_m"]
    if abs(frac - 0.5) < 0.02:
        return "0.5 tau_m"
    if abs(gchi - 1.0) < 0.05:
        return "1.0 tau_c"
    if abs(gchi - 0.3) < 0.05:
        return "0.3 tau_c"
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tchi")
    ap.add_argument("fine")
    ap.add_argument("sym")
    ap.add_argument("--calib", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    with open(a.calib) as fh:
        cal = json.load(fh)
    tau_c = cal["tau_c"]

    runs = (scan(a.tchi, r"(?:t1|th|t0p3)_tm(.+)_([abc]\d?)$")
            + scan(a.fine, r"tm(.+)_([abc]\d?)$")
            + scan(a.sym, r"s1_tm(.+)_([ab])$"))

    # R[(rule, tau_m, arm)] -> the numbers. A run only enters if its tau_m is one of the
    # five and its tau_chi matches a rule; 0.5 tau_m at tau_m = 0.3 tau_c is checked first
    # so it is not mistaken for a fixed rule.
    R = {}
    for c, d in runs:
        u, fl = d["params"], d["flow"]
        s, f = d.get("settled") or {}, d.get("fate") or {}
        g = u["tau_m"] / tau_c
        hit = [t for t in TAU_M if abs(g - t) / t < TOL]
        r = rule_of(u, tau_c)
        arm = re.search(r"_([abc]\d?)$", c).group(1)
        if not hit or not r:
            continue
        dep = f.get("departure_steps", float("nan"))
        sd = fl["std_m"]
        R[(r, hit[0], arm)] = dict(
            case=c, chi=f.get("chi_tail", s.get("chi_mean_tail", float("nan"))),
            std_chi=fl["std_chi"], label=f.get("label", "?"),
            drift=bool(f.get("drifting", False)),
            life=dep / u["tau_m"] if dep == dep else float("inf"),
            margin=abs(s.get("m_mean_tail", float("nan")) - u["mc"]) / sd if sd > 0
            else float("nan"),
            tau_chi=u["tau_chi"], tau_m=u["tau_m"])

    def sep(r, g):
        v = [R[(r, g, k)]["chi"] for k, _, _, _ in ARMS if (r, g, k) in R]
        return max(v) - min(v) if len(v) > 1 else float("nan")

    fig, ax = plt.subplots(2, 2, figsize=(11.5, 7.6), constrained_layout=True)

    # 1. <chi> per arm, per rule
    for rl, col in RULES:
        for arm, alab, ls, mk in ARMS:
            x = [g for g in TAU_M if (rl, g, arm) in R]
            if not x:
                continue
            ax[0, 0].plot(x, [R[(rl, g, arm)]["chi"] for g in x], ls, color=col, lw=1.3,
                          marker=mk, ms=4, label=f"{rl}, {alab}")
    ax[0, 0].set_ylabel(r"tail $\langle\chi\rangle$")
    ax[0, 0].set_ylim(-0.05, 1.05)
    ax[0, 0].legend(fontsize=6, ncol=3, loc="center left")
    ax[0, 0].set_title("fate, per initial condition and rule", fontsize=9)

    # 2. separation -- the one number that answers the question
    for rl, col in RULES:
        x = [g for g in TAU_M if any((rl, g, k) in R for k, _, _, _ in ARMS)]
        if not x:
            continue
        ax[0, 1].plot(x, [sep(rl, g) for g in x], "-o", color=col, lw=1.6, ms=5, label=rl)
    ax[0, 1].axhline(0.5, color="0.5", ls="--", lw=1)
    ax[0, 1].set_ylabel("separation  max - min over starts")
    ax[0, 1].set_ylim(-0.05, 1.05)
    ax[0, 1].legend(fontsize=8, title=r"$\tau_\chi$", title_fontsize=8)
    ax[0, 1].set_title("if the fate is a function of sigma_m alone, these collapse",
                       fontsize=9)

    # 3. the active arm's lifetime
    for rl, col in RULES:
        x = [g for g in TAU_M if (rl, g, "a") in R]
        if not x:
            continue
        y = [R[(rl, g, "a")]["life"] for g in x]
        fin = [(g, v) for g, v in zip(x, y) if np.isfinite(v)]
        inf = [g for g, v in zip(x, y) if not np.isfinite(v)]
        if fin:
            ax[1, 0].plot([g for g, _ in fin], [v for _, v in fin], "-o", color=col,
                          lw=1.6, ms=5, label=rl)
        for g in inf:
            ax[1, 0].plot([g], [ax[1, 0].get_ylim()[1]], "^", color=col, ms=9,
                          mfc="none", mew=1.6)
    ax[1, 0].set_yscale("log")
    ax[1, 0].set_ylabel(r"active-arm lifetime  [$\tau_m$]")
    ax[1, 0].set_xlabel(r"$\tau_m/\tau_c$")
    ax[1, 0].legend(fontsize=8, title=r"$\tau_\chi$", title_fontsize=8)
    ax[1, 0].set_title("open triangle at the top = never departed", fontsize=9)

    # 4. std(chi): the texture, not the mean
    for rl, col in RULES:
        for arm, alab, ls, mk in (ARMS[0], ARMS[1]):
            x = [g for g in TAU_M if (rl, g, arm) in R]
            if not x:
                continue
            ax[1, 1].plot(x, [R[(rl, g, arm)]["std_chi"] for g in x], ls, color=col,
                          lw=1.3, marker=mk, ms=4, label=f"{rl}, {alab}")
    ax[1, 1].set_ylabel(r"$\mathrm{std}(\chi)$")
    ax[1, 1].set_xlabel(r"$\tau_m/\tau_c$")
    ax[1, 1].legend(fontsize=6, ncol=2)
    ax[1, 1].set_title("phenotype texture", fontsize=9)

    for axx in ax.ravel():
        axx.set_xscale("log")
        axx.grid(alpha=0.25)
    fig.suptitle(r"phenotype clock at the symmetric point ($m_c$ = "
                 f"{cal['mc']:.4g}"
                 r"): three $\tau_\chi$ rules on five $\tau_m$", fontsize=10)
    fig.savefig(os.path.join(a.out, "tchi.png"), dpi=200)
    plt.close(fig)

    summary = {"mc": cal["mc"], "tau_c": tau_c, "tau_m_grid": TAU_M,
               "rules": [r for r, _ in RULES],
               "separation": {r: {str(g): sep(r, g) for g in TAU_M} for r, _ in RULES},
               "runs": {f"{r}|{g}|{k}": {kk: vv for kk, vv in v.items()}
                        for (r, g, k), v in sorted(R.items(), key=lambda x: str(x[0]))}}
    with open(os.path.join(a.out, "tchi.json"), "w") as fh:
        json.dump(summary, fh, indent=1, default=lambda o: None)

    print(f"{len(R)} runs matched into the grid")
    print(f"{'tau_m/tau_c':>11} " + " ".join(f"{r:>26}" for r, _ in RULES))
    for g in TAU_M:
        cells = []
        for rl, _ in RULES:
            bits = []
            for arm, _l, _s, _m in ARMS:
                v = R.get((rl, g, arm))
                bits.append(f"{v['chi']:.3f}" if v else "  -  ")
            cells.append(" ".join(bits) + f"  d={sep(rl, g):.2f}"
                         if any((rl, g, k) in R for k, _, _, _ in ARMS) else " " * 26)
        print(f"{g:11g} " + " ".join(f"{c:>26}" for c in cells))
    print("\n(each cell: a  b  c1, then the separation)")
    print(f"wrote {a.out}/tchi.png and tchi.json")


if __name__ == "__main__":
    main()
