"""Clean summary figures for the go-grow interface / director study.

Reads director_metrics.csv (variant, alpha, zeta, Ochi, LL, xi_align, ell_d, A,
lam) and produces a small set of publication-style figures, each showing only
the clearest relationship, with proper math labels (alpha, zeta, lambda, ell).

Usage: python summary_figures.py <director_metrics.csv> <outdir>
"""
import matplotlib
matplotlib.use("Agg")

import argparse
import csv
import os

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 12, "axes.grid": True, "grid.alpha": 0.3,
                     "axes.titlesize": 12, "legend.fontsize": 10})
BLUE, RED, ORANGE, GREEN, PURPLE = "#1f77b4", "#d62728", "#ff7f0e", "#2ca02c", "#9467bd"


def load(path):
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            d = {}
            for k, v in r.items():
                if k == "variant":
                    d[k] = v
                else:
                    try:
                        d[k] = float(v)
                    except (TypeError, ValueError):
                        d[k] = np.nan
            rows.append(d)
    return rows


def line(rows, fixed, xkey, ykey, xfun=None):
    pts = []
    for r in rows:
        if any(not np.isfinite(r.get(k, np.nan)) or abs(r[k] - v) > 1e-12 for k, v in fixed.items()):
            continue
        y = r.get(ykey)
        if y is None or not np.isfinite(y):
            continue
        x = xfun(r) if xfun else r.get(xkey)
        if x is None or not np.isfinite(x):
            continue
        pts.append((x, y))
    return sorted(pts)


def plot_line(ax, pts, **kw):
    if pts:
        ax.plot([p[0] for p in pts], [p[1] for p in pts], **kw)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("csv")
    ap.add_argument("outdir")
    ap.add_argument("--dpi", type=int, default=160)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rows = load(args.csv)
    sqf = lambda r: np.sqrt(r["LL"] / r["zeta"])
    refA = dict(alpha=1e-4, zeta=0.05, Ochi=0.01, xi_align=0.3)
    refB = dict(alpha=1e-4, zeta=0.1, Ochi=0.005, xi_align=0.3)

    # ---- Fig 1: wavelength vs elasticity ----
    fig, ax = plt.subplots(figsize=(6, 4.6))
    plot_line(ax, line(rows, refA, "LL", "lam"), marker="o", color=BLUE, lw=1.8, label=r"$\zeta=0.05$")
    plot_line(ax, line(rows, refB, "LL", "lam"), marker="s", color=RED, lw=1.8, label=r"$\zeta=0.1$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"nematic elasticity  $L$"); ax.set_ylabel(r"interface wavelength  $\lambda$")
    ax.set_title(r"Wavelength grows with elasticity  ($\alpha=10^{-4}$)")
    ax.legend(title="activity")
    fig.tight_layout(); fig.savefig(os.path.join(args.outdir, "fig1_wavelength_vs_elasticity.png"), dpi=args.dpi); plt.close(fig)

    # ---- Fig 2: no active-length collapse ----
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    plot_line(ax, line(rows, dict(alpha=1e-4, LL=0.002, Ochi=0.01, xi_align=0.3), "zeta", "lam", sqf),
              marker="o", color=BLUE, lw=1.8, label=r"vary $\zeta$  (fixed $L$)")
    plot_line(ax, line(rows, refA, "LL", "lam", sqf), marker="s", ls="--", color=RED, lw=1.8, label=r"vary $L$  ($\zeta=0.05$)")
    plot_line(ax, line(rows, refB, "LL", "lam", sqf), marker="s", ls="--", color=ORANGE, lw=1.8, label=r"vary $L$  ($\zeta=0.1$)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"active length  $\sqrt{L/\zeta}$"); ax.set_ylabel(r"interface wavelength  $\lambda$")
    ax.set_title(r"$\lambda$ set by $L$, not $\sqrt{L/\zeta}$ (no collapse)")
    ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(args.outdir, "fig2_wavelength_no_active_collapse.png"), dpi=args.dpi); plt.close(fig)

    # ---- Fig 3: director coherence length (vs L, vs zeta) ----
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.4))
    ax = axs[0]
    pA = line(rows, refA, "LL", "ell_d"); pB = line(rows, refB, "LL", "ell_d")
    plot_line(ax, pA, marker="o", color=BLUE, lw=1.8, label=r"$\zeta=0.05$")
    plot_line(ax, pB, marker="s", color=RED, lw=1.8, label=r"$\zeta=0.1$")
    if pA:
        L0, e0 = pA[0]; LL = np.array([p[0] for p in pA])
        ax.plot(LL, e0 * np.sqrt(LL / L0), "k--", alpha=0.5, label=r"$\sqrt{L}$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"nematic elasticity  $L$"); ax.set_ylabel(r"director coherence length  $\ell$")
    ax.set_title(r"$\ell \sim \sqrt{L}$  ($\alpha=10^{-4}$)"); ax.legend()
    ax = axs[1]
    for al, c in [(2e-4, GREEN), (5e-5, PURPLE)]:
        plot_line(ax, line(rows, dict(alpha=al, LL=0.002, Ochi=0.01, xi_align=0.3), "zeta", "ell_d"),
                  marker="o", color=c, lw=1.8, label=r"$\alpha=%g$" % al)
    ax.set_xscale("log")
    ax.set_xlabel(r"activity  $\zeta$"); ax.set_ylabel(r"director coherence length  $\ell$")
    ax.set_title(r"activity shortens $\ell$"); ax.legend(title="proliferation")
    fig.tight_layout(); fig.savefig(os.path.join(args.outdir, "fig3_coherence_length.png"), dpi=args.dpi); plt.close(fig)

    # ---- Fig 4: anchoring vs activity ----
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    for al, c in [(2e-4, GREEN), (1e-4, BLUE), (5e-5, PURPLE)]:
        plot_line(ax, line(rows, dict(alpha=al, LL=0.002, Ochi=0.01, xi_align=0.3), "zeta", "A"),
                  marker="o", color=c, lw=1.8, label=r"$\alpha=%g$" % al)
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_xscale("log"); ax.set_ylim(-0.2, 1.05)
    ax.set_xlabel(r"activity  $\zeta$"); ax.set_ylabel(r"anchoring index  $\langle\cos 2\theta\rangle$")
    ax.set_title("Tangential anchoring destroyed by activity")
    ax.text(0.03, 0.06, "+1 = parallel,  0 = none", transform=ax.transAxes, fontsize=10,
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="0.6"))
    ax.legend(title="proliferation")
    fig.tight_layout(); fig.savefig(os.path.join(args.outdir, "fig4_anchoring_vs_activity.png"), dpi=args.dpi); plt.close(fig)

    print("wrote 4 summary figures to", args.outdir)


if __name__ == "__main__":
    main()
