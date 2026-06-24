"""Go-grow interface roughness analysis for the boundary scan.

For each variant (final frame): extract the chi=0.5 contour inside the colony
(phi>0.5), take the longest contour as the main go-grow interface, and
characterise it from the radial fluctuation spectrum via the DOMINANT (peak)
angular mode n* = argmax P(n>=2):
  lam   = characteristic wavelength = 2*pi*<r> / n*
  w_dom = amplitude of the dominant mode = 2*|F_n*| / N
plus dimensionless shape metrics Q = L^2/(4 pi A) - 1 and solidity = A / A_hull,
and n_dom = n*. Relative versions (lam_rel, w_dom_rel) divide by <r>.

Writes interface_roughness.csv, correlation_summary.txt, and figures (Spearman
heatmap, representative-slice trends, alpha x zeta heatmaps, LL/xi study, and
the active-length sqrt(LL/zeta) scaling).

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
COLS = ["variant", "alpha", "zeta", "Ochi", "LL", "xi_align", "R_eff",
        "w_dom", "w_dom_rel", "lam", "lam_rel", "n_dom", "Q", "solidity",
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

    # resample r(theta) and take the angular power spectrum
    order = np.argsort(th)
    N = 512
    tg = np.linspace(-np.pi, np.pi, N, endpoint=False)
    rg = np.interp(tg, th[order], rr[order], period=2 * np.pi)
    h = rg - rg.mean()
    F = np.fft.rfft(h)
    P = np.abs(F) ** 2
    n = np.arange(len(P))
    nmax = max(3, int(rmean / 2))                 # pixel-resolution cutoff
    sel = (n >= 2) & (n <= nmax)
    if P[sel].sum() > 0:
        n_dom = int(n[sel][np.argmax(P[sel])])
        lam = 2 * np.pi * rmean / n_dom            # wavelength from the dominant (peak) mode
        w_dom = float(2.0 * np.abs(F[n_dom]) / N)  # amplitude of the dominant mode
    else:
        n_dom = lam = w_dom = np.nan

    return dict(R_eff=R_eff,
                w_dom=w_dom, w_dom_rel=(w_dom / rmean if np.isfinite(w_dom) else np.nan),
                lam=lam, lam_rel=(lam / rmean if np.isfinite(lam) else np.nan),
                n_dom=n_dom, Q=Q, solidity=solidity,
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
    ap.add_argument("--pattern", default="*")
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
        pp = parse_params(name)
        row = dict(variant=name, alpha=pp.get("alpha"), zeta=pp.get("zeta"),
                   Ochi=pp.get("Ochi"), LL=pp.get("LL"), xi_align=pp.get("xi"))
        if m:
            row.update(m)
        rows.append(row)
        print("ok " if m else "no-interface ", name)

    with open(os.path.join(args.outdir, "interface_roughness.csv"), "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=COLS)
        wr.writeheader()
        for r in rows:
            wr.writerow({k: r.get(k) for k in COLS})

    def col(k, rl):
        return np.array([(r.get(k) if r.get(k) is not None else np.nan) for r in rl], float)

    def is_main(r):
        return (r.get("LL") is not None and abs(r["LL"] - 0.002) < 1e-9
                and r.get("xi_align") is not None and abs(r["xi_align"] - 0.3) < 1e-9)
    MAIN = [r for r in rows if is_main(r)]

    params = {"alpha": col("alpha", MAIN), "zeta": col("zeta", MAIN), "Ochi": col("Ochi", MAIN)}
    metric_keys = ["w_dom", "w_dom_rel", "lam", "lam_rel", "n_dom", "Q", "solidity"]
    metrics = {k: col(k, MAIN) for k in metric_keys}

    # Spearman correlations
    rho_mat = np.full((len(metric_keys), 3), np.nan)
    lines = ["Spearman rho (p) of metric vs parameter (main grid, n=%d cases):" % len(MAIN)]
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
    for yk in ["lam", "w_dom", "lam_rel", "w_dom_rel", "n_dom"]:
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

    # 2) representative slices: metric vs each param, 4 lines = 2x2 extremes of the other two
    def lohi(arr):
        u = sorted({v for v in arr if np.isfinite(v)})
        return [u[0], u[-1]] if len(u) >= 2 else u
    pvals = {p: lohi(params[p]) for p in params}
    colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd"]
    for mk, ml in [("lam", "wavelength lam (n*)"), ("w_dom", "dominant amplitude w_dom"),
                   ("Q", "isoperimetric Q-1"), ("solidity", "solidity A/A_hull")]:
        y = metrics[mk]
        fig, axs = plt.subplots(1, 3, figsize=(13, 3.6), squeeze=False)
        for cj, xp in enumerate(params):
            ax = axs[0][cj]
            others = [p for p in params if p != xp]
            xvals = sorted({v for v in params[xp] if np.isfinite(v)})
            ci = 0
            for a_ in pvals[others[0]]:
                for b_ in pvals[others[1]]:
                    xs, ys = [], []
                    for xv in xvals:
                        sel = ((np.abs(params[xp] - xv) < 1e-12)
                               & (np.abs(params[others[0]] - a_) < 1e-12)
                               & (np.abs(params[others[1]] - b_) < 1e-12) & np.isfinite(y))
                        idx = np.where(sel)[0]
                        if len(idx):
                            xs.append(xv); ys.append(float(np.mean(y[idx])))
                    if xs:
                        ax.plot(xs, ys, "o-", color=colors[ci], lw=1.6, ms=5,
                                label="%s=%g, %s=%g" % (others[0], a_, others[1], b_))
                    ci += 1
            ax.set_xscale("log"); ax.set_xlabel(xp, fontsize=10); ax.grid(True, alpha=0.3)
            if cj == 0:
                ax.set_ylabel(ml, fontsize=10)
            ax.legend(fontsize=7, framealpha=0.6)
        fig.suptitle("%s vs parameters (4 representative slices = extremes of the other two)" % ml, fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        fig.savefig(os.path.join(args.outdir, "trend_%s.png" % mk), dpi=args.dpi); plt.close(fig)

    # 3) alpha x zeta heatmaps per Ochi: shared color scale + cell annotations
    fin = lambda k: sorted({v for v in col(k, MAIN) if np.isfinite(v)})
    alphas, zetas, ochis = fin("alpha"), fin("zeta"), fin("Ochi")

    def heatmaps(yk, label, lognorm):
        vals = col(yk, MAIN); vals = vals[np.isfinite(vals)]
        if not len(vals):
            return
        vmin, vmax = np.percentile(vals, 5), np.percentile(vals, 95)
        norm = LogNorm(max(vmin, 1e-9), vmax) if lognorm else None
        fig, axs = plt.subplots(1, len(ochis), figsize=(3.5 * len(ochis), 3.6), squeeze=False)
        im = None
        for k, oc in enumerate(ochis):
            M = np.full((len(zetas), len(alphas)), np.nan)
            for r in MAIN:
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
    heatmaps("w_dom_rel", "relative dominant amplitude w_dom/<r> (shared scale)", False)

    # 4) LL / xi study at the reference points (groups with multiple LL,xi)
    from collections import defaultdict
    groups = defaultdict(list)
    for r in rows:
        groups[(r.get("alpha"), r.get("zeta"), r.get("Ochi"))].append(r)
    refpts = [(key, rs) for key, rs in sorted(groups.items())
              if len({(rr.get("LL"), rr.get("xi_align")) for rr in rs}) > 1]
    if refpts:
        for mk, ml in [("lam", "wavelength lam (n*)"), ("w_dom", "dominant amplitude w_dom"),
                       ("Q", "isoperimetric Q-1")]:
            fig, axs = plt.subplots(1, len(refpts), figsize=(5.2 * len(refpts), 4), squeeze=False)
            for ci, (key, rs) in enumerate(refpts):
                ax = axs[0][ci]
                for xv in sorted({rr.get("xi_align") for rr in rs if rr.get("xi_align") is not None}):
                    pts = sorted((rr.get("LL"), rr.get(mk)) for rr in rs
                                 if rr.get("xi_align") == xv and rr.get("LL") is not None
                                 and rr.get(mk) is not None and np.isfinite(rr.get(mk, np.nan)))
                    if pts:
                        ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", lw=1.7, ms=6, label="xi=%g" % xv)
                if mk == "lam":
                    LLs = sorted({rr.get("LL") for rr in rs if rr.get("LL")})
                    b = [rr.get(mk) for rr in rs if rr.get("LL") == LLs[0] and rr.get("xi_align") == 0.3
                         and rr.get(mk) is not None and np.isfinite(rr.get(mk, np.nan))]
                    if b and len(LLs) > 1:
                        ax.plot(LLs, [b[0] * np.sqrt(L / LLs[0]) for L in LLs], "k--", alpha=0.5, label="~sqrt(LL)")
                a_, z_, o_ = key
                ax.set_xscale("log"); ax.set_xlabel("LL"); ax.grid(True, alpha=0.3)
                if ci == 0:
                    ax.set_ylabel(ml)
                ax.set_title("alpha=%g, zeta=%g, Ochi=%g" % (a_, z_, o_), fontsize=9)
                ax.legend(fontsize=8)
            fig.suptitle("%s vs LL (lines: xi) at the reference points" % ml, fontsize=12)
            fig.tight_layout(rect=(0, 0, 1, 0.94))
            fig.savefig(os.path.join(args.outdir, "llxi_%s.png" % mk), dpi=args.dpi); plt.close(fig)
    else:
        print("(no LL/xi reference groups found)")

    # 5) characteristic length vs active length sqrt(LL/zeta): controlled lines.
    #    zeta-scan (vary zeta at fixed LL) vs LL-scan (vary LL at fixed zeta);
    #    if both collapse onto one curve, lam depends only on sqrt(LL/zeta).
    A = 1e-4

    def active_line(fixed, mk):
        pts = []
        for r in rows:
            if any(r.get(k) is None or abs(r.get(k) - v) > 1e-12 for k, v in fixed.items()):
                continue
            yv = r.get(mk)
            if r.get("LL") is None or r.get("zeta") is None or yv is None or not np.isfinite(yv):
                continue
            pts.append((np.sqrt(r["LL"] / r["zeta"]), yv))
        return sorted(pts)

    zeta_lines = [dict(alpha=A, LL=0.002, xi_align=0.3, Ochi=o) for o in (0.01, 0.005)]
    LL_lines = [dict(alpha=A, zeta=0.05, Ochi=0.01, xi_align=0.3),
                dict(alpha=A, zeta=0.1, Ochi=0.005, xi_align=0.3)]
    summ2 = []
    for mk, ml in [("lam", "wavelength lam (n*)"), ("w_dom", "dominant amplitude w_dom")]:
        ALA = [r for r in rows if r.get("alpha") is not None and abs(r["alpha"] - A) < 1e-12
               and r.get("LL") and r.get("zeta")]
        Lact = np.array([np.sqrt(r["LL"] / r["zeta"]) for r in ALA])
        yall = np.array([(r.get(mk) if r.get(mk) is not None else np.nan) for r in ALA], float)
        ok = np.isfinite(yall) & np.isfinite(Lact) & (yall > 0) & (Lact > 0)
        if ok.sum() >= 4:
            slope, _ = np.polyfit(np.log10(Lact[ok]), np.log10(yall[ok]), 1)
            rho, pp = spearmanr(Lact[ok], yall[ok])
            summ2.append("%s ~ sqrt(LL/zeta) (alpha=1e-4, n=%d): log-log slope=%.2f, Spearman=%.2f (p=%.2g)"
                         % (mk, int(ok.sum()), slope, rho, pp))
        fig, ax = plt.subplots(figsize=(6.6, 5))
        allp = []
        for fx in zeta_lines:
            pts = active_line(fx, mk); allp += pts
            if len(pts) >= 2:
                ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", color="#1f77b4",
                        lw=1.7, ms=6, label="zeta-scan (LL=0.002, Ochi=%g)" % fx["Ochi"])
        for fx in LL_lines:
            pts = active_line(fx, mk); allp += pts
            if len(pts) >= 2:
                ax.plot([p[0] for p in pts], [p[1] for p in pts], "s--", color="#d62728",
                        lw=1.7, ms=6, label="LL-scan (zeta=%g, Ochi=%g)" % (fx["zeta"], fx["Ochi"]))
        if allp:
            allp.sort(); x0, y0 = allp[0]; xr = np.sort(np.array([p[0] for p in allp]))
            ax.plot(xr, y0 * xr / x0, "k:", alpha=0.6, label="slope 1")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("active length  sqrt(LL/zeta)"); ax.set_ylabel(ml)
        ax.set_title("%s vs active length (alpha=1e-4): zeta-scan vs LL-scan" % ml, fontsize=10)
        ax.legend(fontsize=8); ax.grid(True, which="both", alpha=0.3)
        fig.tight_layout(); fig.savefig(os.path.join(args.outdir, "active_%s.png" % mk), dpi=args.dpi); plt.close(fig)
    if summ2:
        block = "\nActive-length scaling:\n  " + "\n  ".join(summ2)
        print(block)
        with open(os.path.join(args.outdir, "correlation_summary.txt"), "a") as fh:
            fh.write(block + "\n")

    n_ok = sum(1 for r in rows if r.get("lam") is not None and np.isfinite(r.get("lam", np.nan)))
    print("\nprocessed %d cases (%d with a measurable interface). wrote to %s" % (len(rows), n_ok, args.outdir))


if __name__ == "__main__":
    main()
