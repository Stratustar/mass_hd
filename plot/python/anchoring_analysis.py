"""Definitive analysis of go-grow interface anchoring A = <cos 2(theta_dir - tangent)>.

Reads director_metrics.csv. Quantifies how A correlates with the control
parameters (alpha, zeta, Ochi, L, xi), with the director coherence length ell_d,
and with the derived ratio alpha/zeta; runs a standardized regression; classifies
each case as anchored / intermediate / broken; and makes clean figures:
  fig_A_regime        : A on the alpha x zeta plane (mean over Ochi)
  fig_A_vs_zeta       : A vs activity zeta, one line per alpha
  fig_A_vs_ratio      : A vs alpha/zeta (collapse test) + A vs ell_d

Usage: python anchoring_analysis.py <director_metrics.csv> <outdir>
"""
import matplotlib
matplotlib.use("Agg")

import argparse
import csv
import os

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

plt.rcParams.update({"font.size": 12, "axes.grid": True, "grid.alpha": 0.3})


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


def col(rows, k):
    return np.array([(r.get(k) if r.get(k) is not None else np.nan) for r in rows], float)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("csv")
    ap.add_argument("outdir")
    ap.add_argument("--dpi", type=int, default=160)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rows = load(args.csv)
    for r in rows:
        z = r.get("zeta")
        r["ratio"] = (r.get("alpha") / z) if (z and np.isfinite(z) and z != 0) else np.nan

    base = lambda r: (np.isfinite(r.get("LL", np.nan)) and abs(r["LL"] - 0.002) < 1e-9
                      and np.isfinite(r.get("xi_align", np.nan)) and abs(r["xi_align"] - 0.3) < 1e-9)
    MAIN = [r for r in rows if base(r)]

    quantities = ["alpha", "zeta", "Ochi", "LL", "xi_align", "ell_d", "lam", "ratio"]
    lines = []
    for scope, rs in [("all (n=%d)" % len(rows), rows), ("main grid (n=%d)" % len(MAIN), MAIN)]:
        A = col(rs, "A")
        lines.append("Spearman rho of A vs quantity  [%s]:" % scope)
        for q in quantities:
            x = col(rs, q); ok = np.isfinite(A) & np.isfinite(x)
            if ok.sum() >= 5 and len(set(x[ok])) > 1:
                rho, pp = spearmanr(x[ok], A[ok])
                lines.append("    %-9s rho=%+.2f (p=%.2g)" % (q, rho, pp))
        lines.append("")

    # standardized regression A ~ log10(alpha, zeta, Ochi) on main grid
    A = col(MAIN, "A")
    P = {p: col(MAIN, p) for p in ("alpha", "zeta", "Ochi")}
    ok = np.isfinite(A)
    for v in P.values():
        ok &= np.isfinite(v) & (v > 0)
    if ok.sum() >= 6:
        X = np.column_stack([np.log10(P[p][ok]) for p in P])
        X = (X - X.mean(0)) / X.std(0)
        yy = (A[ok] - A[ok].mean()) / A[ok].std()
        coef, *_ = np.linalg.lstsq(np.column_stack([X, np.ones(ok.sum())]), yy, rcond=None)
        lines.append("Standardized regression A ~ log10 params (relative importance):")
        lines.append("    " + ", ".join("%s=%+.2f" % (p, c) for p, c in zip(P, coef[:3])))
        lines.append("")

    # regime classification
    anc = [r for r in rows if np.isfinite(r.get("A", np.nan)) and r["A"] > 0.5]
    brk = [r for r in rows if np.isfinite(r.get("A", np.nan)) and r["A"] < 0.2]
    def rng(rs, k):
        v = sorted({r.get(k) for r in rs if np.isfinite(r.get(k, np.nan))})
        return "{%s}" % ",".join("%g" % x for x in v)
    lines.append("Anchored (A>0.5): n=%d   zeta in %s, alpha in %s" % (len(anc), rng(anc, "zeta"), rng(anc, "alpha")))
    lines.append("Broken   (A<0.2): n=%d   zeta in %s, alpha in %s" % (len(brk), rng(brk, "zeta"), rng(brk, "alpha")))
    summary = "\n".join(lines)
    print(summary)
    with open(os.path.join(args.outdir, "anchoring_summary.txt"), "w") as fh:
        fh.write(summary + "\n")

    alphas = sorted({r["alpha"] for r in MAIN if np.isfinite(r.get("alpha", np.nan))})
    zetas = sorted({r["zeta"] for r in MAIN if np.isfinite(r.get("zeta", np.nan))})

    # Fig 1: regime heatmap A on alpha x zeta (mean over Ochi)
    M = np.full((len(zetas), len(alphas)), np.nan)
    for iz, z in enumerate(zetas):
        for ia, a in enumerate(alphas):
            vals = [r["A"] for r in MAIN if r.get("zeta") == z and r.get("alpha") == a
                    and np.isfinite(r.get("A", np.nan))]
            if vals:
                M[iz, ia] = np.mean(vals)
    fig, ax = plt.subplots(figsize=(6, 4.4))
    im = ax.imshow(M, origin="lower", aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1)
    for iz in range(len(zetas)):
        for ia in range(len(alphas)):
            if np.isfinite(M[iz, ia]):
                ax.text(ia, iz, "%.2f" % M[iz, ia], ha="center", va="center",
                        color="white" if abs(M[iz, ia]) > 0.5 else "black")
    ax.set_xticks(range(len(alphas))); ax.set_xticklabels(["%g" % a for a in alphas])
    ax.set_yticks(range(len(zetas))); ax.set_yticklabels(["%g" % z for z in zetas])
    ax.set_xlabel(r"proliferation  $\alpha$"); ax.set_ylabel(r"activity  $\zeta$")
    ax.set_title(r"anchoring $\langle\cos 2\theta\rangle$  (mean over $O_\chi$)")
    fig.colorbar(im, ax=ax, label=r"$\langle\cos 2\theta\rangle$  (+1 parallel)")
    fig.tight_layout(); fig.savefig(os.path.join(args.outdir, "fig_A_regime.png"), dpi=args.dpi); plt.close(fig)

    # Fig 2: A vs zeta, one line per alpha (Ochi=0.01)
    fig, ax = plt.subplots(figsize=(6, 4.4))
    cmap = plt.cm.viridis
    for i, a in enumerate(alphas):
        pts = sorted((r["zeta"], r["A"]) for r in MAIN if r.get("alpha") == a
                     and np.isfinite(r.get("Ochi", np.nan)) and abs(r["Ochi"] - 0.01) < 1e-12
                     and np.isfinite(r.get("A", np.nan)))
        if pts:
            ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", color=cmap(i / max(1, len(alphas) - 1)),
                    lw=1.7, ms=6, label=r"$\alpha=%g$" % a)
    ax.axhline(0, color="0.6", lw=0.8); ax.set_xscale("log"); ax.set_ylim(-0.2, 1.05)
    ax.set_xlabel(r"activity  $\zeta$"); ax.set_ylabel(r"anchoring  $\langle\cos 2\theta\rangle$")
    ax.set_title(r"Anchoring destroyed by activity ($O_\chi=0.01$)")
    ax.legend(title="proliferation", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(args.outdir, "fig_A_vs_zeta.png"), dpi=args.dpi); plt.close(fig)

    # Fig 3: collapse test A vs alpha/zeta, and A vs ell_d
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.4))
    A_all = col(MAIN, "A"); rat = col(MAIN, "ratio"); ell = col(MAIN, "ell_d")
    ok = np.isfinite(A_all) & np.isfinite(rat) & (rat > 0)
    sc = axs[0].scatter(rat[ok], A_all[ok], c=np.log10(col(MAIN, "zeta")[ok]), cmap="viridis", s=45)
    if ok.sum() >= 5:
        rho, pp = spearmanr(rat[ok], A_all[ok])
        axs[0].set_title(r"$A$ vs $\alpha/\zeta$  (Spearman=%.2f)" % rho)
    axs[0].axhline(0, color="0.6", lw=0.8); axs[0].set_xscale("log")
    axs[0].set_xlabel(r"$\alpha/\zeta$  (growth / activity)"); axs[0].set_ylabel(r"anchoring  $\langle\cos 2\theta\rangle$")
    fig.colorbar(sc, ax=axs[0], label=r"$\log_{10}\zeta$")
    ok2 = np.isfinite(A_all) & np.isfinite(ell)
    sc2 = axs[1].scatter(ell[ok2], A_all[ok2], c=np.log10(col(MAIN, "zeta")[ok2]), cmap="viridis", s=45)
    if ok2.sum() >= 5:
        rho2, _ = spearmanr(ell[ok2], A_all[ok2])
        axs[1].set_title(r"$A$ vs director coherence $\ell$  (Spearman=%.2f)" % rho2)
    axs[1].axhline(0, color="0.6", lw=0.8)
    axs[1].set_xlabel(r"director coherence  $\ell$"); axs[1].set_ylabel(r"anchoring  $\langle\cos 2\theta\rangle$")
    fig.colorbar(sc2, ax=axs[1], label=r"$\log_{10}\zeta$")
    fig.tight_layout(); fig.savefig(os.path.join(args.outdir, "fig_A_vs_ratio.png"), dpi=args.dpi); plt.close(fig)

    print("\nwrote anchoring analysis to", args.outdir)


if __name__ == "__main__":
    main()
