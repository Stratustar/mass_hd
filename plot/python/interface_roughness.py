"""Go-grow interface roughness analysis for the boundary scan.

For each variant (final frame): extract the chi=0.5 contour inside the colony
(phi>0.5), take the longest contour as the main go-grow interface, and
characterise its roughness by LENGTH SCALES:
  w   = RMS radial deviation from the mean radius (amplitude, lattice units)
  lam = characteristic lateral wavelength from the radial fluctuation spectrum,
        lam = 2*pi*<r> / <n>,  <n> = spectral-weighted mean angular mode (n>=2)
  xi  = lateral correlation length (1/e decay of the angular autocorrelation)
plus dimensionless shape metrics Q = L^2/(4 pi A) - 1 and solidity = A / A_hull,
and the dominant mode n_dom (finger count). Relative versions divide by <r>.

Writes interface_roughness.csv, correlation_summary.txt, and figures
(Spearman heatmap, scatter vs params, alpha x zeta heatmaps per Ochi).

Usage:
  python interface_roughness.py <scratch_study_dir> <outdir> [--pattern GLOB]
"""
import matplotlib
matplotlib.use("Agg")

import argparse
import csv
import glob
import os
import re
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from scipy.spatial import ConvexHull

SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)
from archive.archive import loadarchive

TOKEN = re.compile(r"^([A-Za-z]+)([-0-9pe.]+)$")
KEYMAP = {"a": "alpha", "z": "zeta", "O": "Ochi", "LL": "LL", "xi": "xi"}
COLS = ["variant", "alpha", "zeta", "Ochi", "R_eff", "w", "w_rel", "lam",
        "lam_rel", "xi", "xi_rel", "n_dom", "n_mean", "Q", "solidity",
        "n_contours", "n_pts"]


def parse_params(name):
    p = {}
    for tok in name.split("_"):
        m = TOKEN.match(tok)
        if m:
            try:
                p[KEYMAP.get(m.group(1), m.group(1))] = float(m.group(2).replace("p", "."))
            except ValueError:
                pass
    return p


def grid(frame, field):
    lx = frame.parameters["LX"]
    ly = frame.parameters["LY"]
    return np.array(getattr(frame, field)).reshape((lx, ly))


def polygon_area(x, y):
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def extract_interface(chi, phi):
    """Metrics for the main chi=0.5 contour inside phi>0.5, or None."""
    chim = np.ma.masked_where(phi <= 0.5, chi)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cs = ax.contour(chim.T, levels=[0.5])
    segs = list(cs.allsegs[0]) if len(cs.allsegs) else []
    plt.close(fig)
    if not segs:
        return None

    def perim(s):
        d = np.diff(s, axis=0)
        return float(np.sqrt((d ** 2).sum(axis=1)).sum())

    n_contours = len(segs)
    main = max(segs, key=perim)
    if len(main) < 8:
        return None
    x, y = main[:, 0], main[:, 1]
    L = perim(main)
    A = polygon_area(x, y)
    if A <= 0:
        return None
    R_eff = np.sqrt(A / np.pi)
    Q = L * L / (4 * np.pi * A) - 1.0
    try:
        solidity = A / ConvexHull(main).volume   # 2D ConvexHull.volume == area
    except Exception:
        solidity = np.nan

    # radial fluctuation about the contour centroid
    cx, cy = x.mean(), y.mean()
    rr = np.sqrt((x - cx) ** 2 + (y - cy) ** 2)
    th = np.arctan2(y - cy, x - cx)
    rmean = float(rr.mean())
    w = float(np.sqrt(np.mean((rr - rmean) ** 2)))

    # resample r(theta) for the spectrum
    order = np.argsort(th)
    N = 512
    tg = np.linspace(-np.pi, np.pi, N, endpoint=False)
    rg = np.interp(tg, th[order], rr[order], period=2 * np.pi)
    h = rg - rg.mean()
    P = np.abs(np.fft.rfft(h)) ** 2
    n = np.arange(len(P))
    nmax = max(3, int(rmean / 2))                 # pixel-resolution cutoff
    sel = (n >= 2) & (n <= nmax)
    if P[sel].sum() > 0:
        n_mean = float(np.sum(n[sel] * P[sel]) / np.sum(P[sel]))
        n_dom = int(n[sel][np.argmax(P[sel])])
        lam = 2 * np.pi * rmean / n_mean
    else:
        n_mean = n_dom = lam = np.nan

    # angular autocorrelation -> lateral correlation length (arclength)
    ac = np.correlate(h, h, mode="full")[len(h) - 1:]
    ac = ac / ac[0] if ac[0] != 0 else ac
    below = np.where(ac < 1.0 / np.e)[0]
    xi = float(below[0] * (2 * np.pi / N) * rmean) if len(below) else np.nan

    return dict(R_eff=R_eff, w=w, w_rel=w / rmean,
                lam=lam, lam_rel=(lam / rmean if np.isfinite(lam) else np.nan),
                xi=xi, xi_rel=(xi / rmean if np.isfinite(xi) else np.nan),
                n_dom=n_dom, n_mean=n_mean, Q=Q, solidity=solidity,
                n_contours=n_contours, n_pts=len(main))


def standardized_regression(y, params):
    ok = np.isfinite(y)
    for pv in params.values():
        ok &= np.isfinite(pv)
    if ok.sum() < 6:
        return None
    X = np.column_stack([np.log10(params[k][ok]) for k in params])
    X = (X - X.mean(0)) / X.std(0)
    yy = (y[ok] - y[ok].mean()) / y[ok].std()
    coef, *_ = np.linalg.lstsq(np.column_stack([X, np.ones(ok.sum())]), yy, rcond=None)
    return dict(zip(list(params) + ["const"], coef))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("studydir")
    ap.add_argument("outdir")
    ap.add_argument("--pattern", default="*LL0p002_xi0p3")
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    rows = []
    for d in sorted(glob.glob(os.path.join(args.studydir, args.pattern))):
        if not os.path.isdir(d):
            continue
        name = os.path.basename(d)
        m = None
        try:
            ar = loadarchive(d)
            last = int((ar.nsteps - ar.nstart) / ar.ninfo)
            frame = ar.read_frame(last)
            m = extract_interface(grid(frame, "chi"), grid(frame, "phi"))
        except Exception as e:
            print("FAIL", name, repr(e))
        row = dict(variant=name, **{k: parse_params(name).get(k) for k in ("alpha", "zeta", "Ochi")})
        if m:
            row.update(m)
        rows.append(row)
        print("ok " if m else "no-interface ", name)

    with open(os.path.join(args.outdir, "interface_roughness.csv"), "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=COLS)
        wr.writeheader()
        for r in rows:
            wr.writerow({k: r.get(k) for k in COLS})

    def col(k):
        return np.array([(r.get(k) if r.get(k) is not None else np.nan) for r in rows], float)

    params = {"alpha": col("alpha"), "zeta": col("zeta"), "Ochi": col("Ochi")}
    metric_keys = ["w", "w_rel", "lam", "lam_rel", "xi", "n_dom", "Q", "solidity"]
    metrics = {k: col(k) for k in metric_keys}

    # Spearman correlations
    rho_mat = np.full((len(metric_keys), 3), np.nan)
    lines = ["Spearman rho (p) of metric vs parameter (n=%d cases):" % len(rows)]
    for i, mk in enumerate(metric_keys):
        cells = []
        for j, pk in enumerate(params):
            mv, pv = metrics[mk], params[pk]
            ok = np.isfinite(mv) & np.isfinite(pv)
            if ok.sum() >= 5:
                rho, pp = spearmanr(pv[ok], mv[ok])
                rho_mat[i, j] = rho
                cells.append("%s %+.2f(p=%.2g)" % (pk, rho, pp))
            else:
                cells.append("%s n/a" % pk)
        lines.append("  %-9s " % mk + "  ".join(cells))
    lines.append("")
    lines.append("Standardized regression (metric ~ log10 alpha,zeta,Ochi; coeffs = relative importance):")
    for yk in ["lam", "w", "lam_rel", "w_rel", "n_dom"]:
        c = standardized_regression(metrics[yk], params)
        if c:
            lines.append("  %-9s " % yk + ", ".join("%s=%+.2f" % (k, c[k]) for k in ("alpha", "zeta", "Ochi")))
    summary = "\n".join(lines)
    print("\n" + summary)
    with open(os.path.join(args.outdir, "correlation_summary.txt"), "w") as fh:
        fh.write(summary + "\n")

    from matplotlib.colors import LogNorm

    # 1) Spearman correlation heatmap
    fig, ax = plt.subplots(figsize=(5.4, 5.8))
    im = ax.imshow(rho_mat, cmap="coolwarm", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(3)); ax.set_xticklabels(list(params), fontsize=12)
    ax.set_yticks(range(len(metric_keys))); ax.set_yticklabels(metric_keys, fontsize=11)
    for i in range(len(metric_keys)):
        for j in range(3):
            if np.isfinite(rho_mat[i, j]):
                ax.text(j, i, "%.2f" % rho_mat[i, j], ha="center", va="center",
                        fontsize=11, color="white" if abs(rho_mat[i, j]) > 0.55 else "black")
    ax.set_title("Spearman rho: roughness vs parameter")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout(); fig.savefig(os.path.join(args.outdir, "corr_heatmap.png"), dpi=args.dpi); plt.close(fig)

    # 2) marginal trends: metric vs each param, mean +/- std over the other two axes
    trend = [("lam", "wavelength lam"), ("w_rel", "rel. amplitude w/<r>"),
             ("Q", "isoperimetric Q-1"), ("solidity", "solidity A/A_hull")]
    fig, axs = plt.subplots(len(trend), 3, figsize=(11.5, 2.5 * len(trend)), squeeze=False)
    for ri, (mk, ml) in enumerate(trend):
        y = metrics[mk]
        for cj, pk in enumerate(params):
            ax = axs[ri][cj]
            pv = params[pk]
            xs, mu, sd = [], [], []
            for v in sorted({u for u in pv if np.isfinite(u)}):
                sub = y[(pv == v) & np.isfinite(y)]
                if len(sub):
                    xs.append(v); mu.append(sub.mean()); sd.append(sub.std())
            ax.errorbar(xs, mu, yerr=sd, marker="o", ms=5, capsize=3, lw=1.6, color="#255c99")
            ax.set_xscale("log"); ax.grid(True, alpha=0.3)
            if ri == 0:
                ax.set_title(pk, fontsize=11)
            if ri == len(trend) - 1:
                ax.set_xlabel(pk, fontsize=10)
            if cj == 0:
                ax.set_ylabel(ml, fontsize=10)
    fig.suptitle("Roughness vs parameter (mean +/- std over the other two axes)", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(os.path.join(args.outdir, "trends.png"), dpi=args.dpi); plt.close(fig)

    # 3) alpha x zeta heatmaps per Ochi: shared color scale + cell annotations
    fin = lambda k: sorted({v for v in col(k) if np.isfinite(v)})
    alphas, zetas, ochis = fin("alpha"), fin("zeta"), fin("Ochi")

    def heatmaps(yk, label, lognorm):
        vals = col(yk); vals = vals[np.isfinite(vals)]
        if not len(vals):
            return
        vmin, vmax = np.percentile(vals, 5), np.percentile(vals, 95)
        norm = LogNorm(max(vmin, 1e-9), vmax) if lognorm else None
        fig, axs = plt.subplots(1, len(ochis), figsize=(3.5 * len(ochis), 3.6), squeeze=False)
        im = None
        for k, oc in enumerate(ochis):
            M = np.full((len(zetas), len(alphas)), np.nan)
            for r in rows:
                yv = r.get(yk)
                if r.get("Ochi") == oc and yv is not None and np.isfinite(yv):
                    M[zetas.index(r["zeta"]), alphas.index(r["alpha"])] = yv
            ax = axs[0][k]
            im = ax.imshow(M, origin="lower", aspect="auto", cmap="viridis", norm=norm,
                           vmin=None if lognorm else vmin, vmax=None if lognorm else vmax)
            for iy in range(len(zetas)):
                for ix in range(len(alphas)):
                    if np.isfinite(M[iy, ix]):
                        txt = ("%.0f" % M[iy, ix]) if abs(M[iy, ix]) >= 10 else ("%.2g" % M[iy, ix])
                        ax.text(ix, iy, txt, ha="center", va="center", fontsize=7, color="black",
                                bbox=dict(facecolor="white", alpha=0.55, pad=0.4, edgecolor="none"))
            ax.set_xticks(range(len(alphas))); ax.set_xticklabels(["%g" % a for a in alphas], rotation=45, fontsize=7)
            ax.set_yticks(range(len(zetas))); ax.set_yticklabels(["%g" % z for z in zetas], fontsize=8)
            ax.set_title("Ochi=%g" % oc, fontsize=10); ax.set_xlabel("alpha", fontsize=9)
            if k == 0:
                ax.set_ylabel("zeta", fontsize=9)
        fig.suptitle(label, fontsize=12)
        fig.subplots_adjust(left=0.06, right=0.9, top=0.85, bottom=0.18, wspace=0.35)
        fig.colorbar(im, ax=list(axs[0]), fraction=0.025, pad=0.02)
        fig.savefig(os.path.join(args.outdir, "heatmap_%s.png" % yk), dpi=args.dpi); plt.close(fig)

    heatmaps("lam", "characteristic wavelength lam (lattice units; log color, shared scale)", True)
    heatmaps("w_rel", "relative amplitude w/<r> (shared scale)", False)

    n_ok = sum(1 for r in rows if r.get("lam") is not None and np.isfinite(r.get("lam", np.nan)))
    print("\nprocessed %d cases (%d with a measurable interface). wrote to %s" % (len(rows), n_ok, args.outdir))


if __name__ == "__main__":
    main()
