#!/usr/bin/env python3
"""The one clean figure: measured <chi> of the two uniform starts against the Beta-closure theory.

    cw_s3_fig_clean.py <theory.json> <series root with */series.json> <out.png>

No grid, no title; three legend entries. Theory = variant C (Beta law, m_open): its stable
branches solid, the unstable branch dashed, all in one colour and under one label.
"""
import glob, json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

th = json.load(open(sys.argv[1]))
root, out = sys.argv[2], sys.argv[3]

# measured: closing-100-tau_c mean and in-window std, from the series files
meas = {"chi0": [], "chi1": []}
for p in glob.glob(os.path.join(root, "*", "series.json")):
    d = json.load(open(p))
    if d["start"] in meas:
        meas[d["start"]].append((d["tau_m_over_tau_c"], d["chi_mean_window"], d["chi_std_window"]))
for k in meas:
    meas[k].sort()

# theory: variant C
name = next(k for k in th["branches"] if k.startswith("C"))
br = th["branches"][name]
g = np.asarray(br["g"])
lower = np.array([min(s) if s else np.nan for s in br["stable"]])
upper = np.array([max(s) if len(s) >= 2 and max(s) > 0.5 and min(s) < 0.5 else np.nan for s in br["stable"]])
unst = np.array([u[0] if u else np.nan for u in br["unstable"]])

fig, ax = plt.subplots(figsize=(6.4, 4.4))
ax.plot(g, lower, "-", color="k", lw=1.8, label="theory")
ax.plot(g, upper, "-", color="k", lw=1.8)
ax.plot(g, unst, "--", color="k", lw=1.0, alpha=0.6)
# chi=1 drawn hollow and on top: below the threshold the two starts coincide to 0.008 and
# filled markers of one would simply hide the other
for k, col, mk, mfc, ms in (("chi0", "C3", "o", "C3", 5), ("chi1", "C0", "s", "none", 7)):
    a = np.asarray(meas[k])
    ax.errorbar(a[:, 0], a[:, 1], yerr=a[:, 2], fmt=mk, color=col, mfc=mfc, mew=1.4, ms=ms,
                capsize=2, lw=1, label="chi=0 measured" if k == "chi0" else "chi=1 measured")
ax.set_xscale("log")
ax.set_xlabel(r"$\tau_m/\tau_c$", fontsize=12)
ax.set_ylabel(r"$\langle\chi\rangle$", fontsize=12)
ax.set_ylim(-0.03, 1.03)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.legend(frameon=False, fontsize=10, loc="center left")
fig.tight_layout()
fig.savefig(out, dpi=200)
fig.savefig(os.path.splitext(out)[0] + ".pdf")
print("wrote", out, "and .pdf")
