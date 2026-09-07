#!/usr/bin/env python3
"""Plot the normalised structure factors of cw_s3_sk: rows = starts, columns = tau_m, curves = chi / m / P.

    cw_s3_sk_plot.py <dir with sk_*.json> <out dir>

Two figures per pair of starts: S_mode/Var on log-log (the shape) and the k-shell power / Var
on lin-log (where the variance lives). Dotted verticals mark 2 pi / L_P and 2 pi / xi_N."""
import glob, json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = sys.argv[1]                      # dir with sk_*.json
OUT = sys.argv[2]
G = [(0.3, "0p3"), (3.426, "3p43"), (9.679, "9p68"), (15.932, "15p93"), (30.0, "30")]
STARTS = [("chi0", "chi = 0 start"), ("chi1", "chi = 1 start"),
          ("leftright", "left/right start"), ("patches", "patches start")]
COL = {"chi": "C0", "m": "C3", "P": "0.25"}
LAB = {"chi": r"$\chi$", "m": r"$m$", "P": r"$P$"}

def load(tag, st):
    p = os.path.join(D, f"sk_tm{tag}_{st}.json")
    return json.load(open(p)) if os.path.exists(p) else None

def one_figure(starts, fname, key, ylab, ylog):
    fig, axs = plt.subplots(len(starts), len(G), figsize=(3.6 * len(G), 3.0 * len(starts)),
                            squeeze=False, sharex=True, sharey=True)
    for i, (st, slab) in enumerate(starts):
        for j, (g, tag) in enumerate(G):
            ax = axs[i, j]
            d = load(tag, st)
            if d is None:
                ax.text(0.5, 0.5, "missing", transform=ax.transAxes, ha="center"); continue
            k = np.asarray(d["k"])
            for f in ("chi", "m", "P"):
                v = d["fields"][f]
                if v.get("saturated"):
                    ax.text(0.03, 0.05 + 0.08 * list("chi m P".split()).index(f),
                            f"{LAB[f]} saturated", transform=ax.transAxes, fontsize=8, color=COL[f])
                    continue
                S = np.asarray(v[key])
                ax.plot(k[1:], S[1:], color=COL[f], lw=1.4,
                        label=f"{LAB[f]}  $\\ell_{{1/e}}$={v['r_1e']:.1f}")
            ax.set_xscale("log")
            if ylog: ax.set_yscale("log")
            ax.grid(alpha=0.3, which="both")
            ax.legend(fontsize=7, loc="lower left" if ylog else "upper right")
            if i == 0: ax.set_title(rf"$\tau_m/\tau_c = {g:g}$", fontsize=10)
            if j == 0: ax.set_ylabel(f"{slab}\n{ylab}", fontsize=9)
            if i == len(starts) - 1: ax.set_xlabel(r"$k$  [1/lattice unit]")
    # reference scales
    for ax in axs.ravel():
        for kk, lab in ((2 * np.pi / 3.5, r"$L_P$"), (2 * np.pi / 2.0, r"$\xi_N$")):
            ax.axvline(kk, color="0.7", lw=0.7, ls=":")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, fname), dpi=140)
    plt.close(fig)
    print("wrote", fname)

one_figure(STARTS[:2], "sk_mode_uniform.png", "S_mode_norm", r"$S(k)/\mathrm{Var}$", True)
one_figure(STARTS[:2], "sk_shell_uniform.png", "S_shell_norm", r"$k$-shell power / Var", False)
one_figure(STARTS[2:], "sk_mode_mixed.png", "S_mode_norm", r"$S(k)/\mathrm{Var}$", True)
one_figure(STARTS[2:], "sk_shell_mixed.png", "S_shell_norm", r"$k$-shell power / Var", False)

# the numbers
print(f"\n{'tau_m':>6} {'start':>10} | " + " | ".join(f"{f:>26}" for f in ("chi", "m", "P")))
print(f"{'':>6} {'':>10} | " + " | ".join(f"{'var':>8} {'l_1e':>6} {'l_peak':>6}" for _ in range(3)))
for g, tag in G:
    for st, _ in STARTS:
        d = load(tag, st)
        if d is None: continue
        cells = []
        for f in ("chi", "m", "P"):
            v = d["fields"][f]
            if v.get("saturated"): cells.append(f"{'saturated':>22}")
            else: cells.append(f"{v['var']:8.2e} {v['r_1e']:6.2f} {2*np.pi/v['k_peak']:6.1f}")
        print(f"{g:6.3f} {st:>10} | " + " | ".join(cells))
