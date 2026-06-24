"""Director-field analysis at the go-grow interface.

For each variant (final frame):
  ell_d = director coherence length from an exponential fit to the nematic
          correlation C(r) = <cos 2(theta(x)-theta(x'))> inside the colony
          (masked FFT autocorrelation, radial average, fit C=exp(-r/ell_d)).
  A     = interface anchoring order parameter <cos 2*Delta> along the chi=0.5
          contour, Delta = director angle - local interface tangent.
          A=+1 parallel (tangential/planar); A=-1 perpendicular (homeotropic).
  lam   = interface wavelength (dominant mode), for the ell_d-vs-lam comparison.

Outputs CSV, Spearman heatmap, controlled-line trends, LL/xi controlled lines,
an ell_d-vs-lam scatter, and a director-vs-interface demo.

Usage: python director_analysis.py <scratch_study_dir> <outdir> [--pattern G] [--demo VARIANT]
"""
import matplotlib
matplotlib.use("Agg")

import argparse
import csv
import glob
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)
from archive.archive import loadarchive
from interface_roughness import grid, parse_params, extract_interface

COLS = ["variant", "alpha", "zeta", "Ochi", "LL", "xi_align", "ell_d", "A", "lam", "n_dom"]


def director_field(frame):
    qxx = grid(frame, "QQxx")
    qyx = grid(frame, "QQyx")
    S = np.sqrt(qxx ** 2 + qyx ** 2)
    theta = 0.5 * np.arctan2(qyx, qxx)
    return theta, S, qxx, qyx


def coherence_length(frame, phi, maxr=80):
    """ell_d from exp fit to the nematic correlation C(r) inside phi>0.5."""
    _, S, qxx, qyx = director_field(frame)
    mask = phi > 0.5
    Sd = np.where(S > 1e-9, S, np.nan)
    ux = np.nan_to_num(np.where(mask, qxx / Sd, 0.0))
    uy = np.nan_to_num(np.where(mask, qyx / Sd, 0.0))
    m = mask.astype(float)
    ny, nx = m.shape
    sh = (2 * ny, 2 * nx)

    def acf(a):
        A = np.fft.rfft2(a, s=sh)
        return np.fft.irfft2(A * np.conj(A), s=sh)

    num = acf(ux) + acf(uy)
    den = acf(m)
    block = slice(0, maxr)
    num = num[block, block]
    den = den[block, block]
    with np.errstate(invalid="ignore", divide="ignore"):
        C = num / den
    dy, dx = np.indices(C.shape)
    r = np.sqrt(dx ** 2 + dy ** 2)
    rb = r.astype(int)
    good = den > (den.max() * 0.02)
    Cr = np.array([C[(rb == k) & good].mean() if np.any((rb == k) & good) else np.nan
                   for k in range(maxr)])
    if not np.isfinite(Cr[0]) or Cr[0] == 0:
        return np.nan, np.arange(maxr), Cr
    Cr = Cr / Cr[0]
    rs = np.arange(maxr)
    sel = np.isfinite(Cr) & (Cr > 0.05) & (rs > 0)
    if sel.sum() < 3:
        return np.nan, rs, Cr
    b = np.polyfit(rs[sel], np.log(Cr[sel]), 1)[0]
    ell = (-1.0 / b) if b < 0 else np.nan
    return ell, rs, Cr


def main_contour(chi, phi):
    chim = np.ma.masked_where(phi <= 0.5, chi)
    fig = plt.figure(); ax = fig.add_subplot(111)
    segs = list(ax.contour(chim.T, levels=[0.5]).allsegs[0])
    plt.close(fig)
    if not segs:
        return None
    main = max(segs, key=lambda s: np.sqrt((np.diff(s, axis=0) ** 2).sum(axis=1)).sum())
    return main if len(main) >= 8 else None


def anchoring(frame, chi, phi):
    """A = <cos 2(theta_n - tangent)> along chi=0.5; +1 parallel, -1 perpendicular."""
    main = main_contour(chi, phi)
    if main is None:
        return np.nan, None, None
    x, y = main[:, 0], main[:, 1]
    theta, _, _, _ = director_field(frame)
    ttan = np.arctan2(np.gradient(y), np.gradient(x))
    ix = np.clip(np.round(x).astype(int), 0, theta.shape[0] - 1)
    iy = np.clip(np.round(y).astype(int), 0, theta.shape[1] - 1)
    delta = theta[ix, iy] - ttan
    A = float(np.mean(np.cos(2 * delta)))
    return A, delta, main


def col(rows, k):
    return np.array([(r.get(k) if r.get(k) is not None else np.nan) for r in rows], float)


def lohi(arr):
    u = sorted({v for v in arr if np.isfinite(v)})
    return [u[0], u[-1]] if len(u) >= 2 else u


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("studydir")
    ap.add_argument("outdir")
    ap.add_argument("--pattern", default="*")
    ap.add_argument("--demo", default="a0p00005_z0p01_O0p02_LL0p002_xi0p3")
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    rows = []
    for d in sorted(glob.glob(os.path.join(args.studydir, args.pattern))):
        if not os.path.isdir(d):
            continue
        name = os.path.basename(d)
        ell = A = lam = ndom = np.nan
        try:
            ar = loadarchive(d)
            fr = ar.read_frame(int((ar.nsteps - ar.nstart) / ar.ninfo))
            chi, phi = grid(fr, "chi"), grid(fr, "phi")
            ell, _, _ = coherence_length(fr, phi)
            A, _, _ = anchoring(fr, chi, phi)
            m = extract_interface(chi, phi)
            if m:
                lam, ndom = m["lam"], m["n_dom"]
        except Exception as e:
            print("FAIL", name, repr(e))
        pp = parse_params(name)
        rows.append(dict(variant=name, alpha=pp.get("alpha"), zeta=pp.get("zeta"),
                         Ochi=pp.get("Ochi"), LL=pp.get("LL"), xi_align=pp.get("xi"),
                         ell_d=ell, A=A, lam=lam, n_dom=ndom))
        print("ok" if np.isfinite(ell) else "no-ell", name)

    with open(os.path.join(args.outdir, "director_metrics.csv"), "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=COLS); wr.writeheader()
        for r in rows:
            wr.writerow({k: r.get(k) for k in COLS})

    base = lambda r: (r.get("LL") is not None and abs(r["LL"] - 0.002) < 1e-9
                      and r.get("xi_align") is not None and abs(r["xi_align"] - 0.3) < 1e-9)
    MAIN = [r for r in rows if base(r)]
    params = ("alpha", "zeta", "Ochi")

    # Spearman heatmap (ell_d, A vs params), main grid
    mets = ["ell_d", "A", "lam"]
    rho = np.full((len(mets), 3), np.nan)
    lines = ["Spearman rho (main grid, n=%d):" % len(MAIN)]
    for i, mk in enumerate(mets):
        mv = col(MAIN, mk)
        cells = []
        for j, pk in enumerate(params):
            pv = col(MAIN, pk); ok = np.isfinite(mv) & np.isfinite(pv)
            if ok.sum() >= 5:
                rr, pp = spearmanr(pv[ok], mv[ok]); rho[i, j] = rr
                cells.append("%s %+.2f(p=%.2g)" % (pk, rr, pp))
        lines.append("  %-6s " % mk + "  ".join(cells))
    summary = "\n".join(lines)
    print("\n" + summary)
    with open(os.path.join(args.outdir, "director_summary.txt"), "w") as fh:
        fh.write(summary + "\n")

    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    im = ax.imshow(rho, cmap="coolwarm", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(3)); ax.set_xticklabels(params)
    ax.set_yticks(range(len(mets))); ax.set_yticklabels(mets)
    for i in range(len(mets)):
        for j in range(3):
            if np.isfinite(rho[i, j]):
                ax.text(j, i, "%.2f" % rho[i, j], ha="center", va="center",
                        color="white" if abs(rho[i, j]) > 0.55 else "black")
    ax.set_title("Spearman: director metrics vs param")
    fig.colorbar(im, ax=ax, fraction=0.046); fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "corr_heatmap.png"), dpi=args.dpi); plt.close(fig)

    # ell_d vs lam (test lam ~ ell_d)
    e = col(rows, "ell_d"); L = col(rows, "lam"); ok = np.isfinite(e) & np.isfinite(L)
    fig, ax = plt.subplots(figsize=(5.5, 5))
    sc = ax.scatter(e[ok], L[ok], c=np.log10(col(rows, "LL")[ok]), cmap="viridis", s=45)
    lim = [min(e[ok].min(), L[ok].min()), max(e[ok].max(), L[ok].max())]
    ax.plot(lim, lim, "k--", alpha=0.5, label="lam = ell_d")
    if ok.sum() >= 4:
        rr, pp = spearmanr(e[ok], L[ok])
        ax.set_title("lam vs director coherence ell_d  (Spearman=%.2f)" % rr, fontsize=10)
    ax.set_xlabel("director coherence ell_d"); ax.set_ylabel("interface wavelength lam")
    ax.legend(fontsize=9); fig.colorbar(sc, ax=ax, label="log10 LL"); fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "ell_vs_lam.png"), dpi=args.dpi); plt.close(fig)

    # controlled-slice trends for ell_d and A (main grid)
    P = {p: col(MAIN, p) for p in params}
    pv2 = {p: lohi(P[p]) for p in params}
    colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd"]
    for mk, ml in [("ell_d", "director coherence ell_d"), ("A", "anchoring A = <cos2(dir-tan)>")]:
        y = col(MAIN, mk)
        fig, axs = plt.subplots(1, 3, figsize=(13, 3.6), squeeze=False)
        for cj, xp in enumerate(params):
            ax = axs[0][cj]; others = [p for p in params if p != xp]
            xs = sorted({v for v in P[xp] if np.isfinite(v)})
            ci = 0
            for a_ in pv2[others[0]]:
                for b_ in pv2[others[1]]:
                    xx, yy = [], []
                    for xv in xs:
                        sel = ((np.abs(P[xp] - xv) < 1e-12) & (np.abs(P[others[0]] - a_) < 1e-12)
                               & (np.abs(P[others[1]] - b_) < 1e-12) & np.isfinite(y))
                        idx = np.where(sel)[0]
                        if len(idx):
                            xx.append(xv); yy.append(float(np.mean(y[idx])))
                    if xx:
                        ax.plot(xx, yy, "o-", color=colors[ci], lw=1.6, ms=5,
                                label="%s=%g,%s=%g" % (others[0], a_, others[1], b_))
                    ci += 1
            ax.set_xscale("log"); ax.set_xlabel(xp); ax.grid(True, alpha=0.3)
            if cj == 0:
                ax.set_ylabel(ml)
            if mk == "A":
                ax.axhline(0, color="0.6", lw=0.8); ax.set_ylim(-1.05, 1.05)
            ax.legend(fontsize=7, framealpha=0.6)
        fig.suptitle("%s vs parameters (representative slices)" % ml, fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        fig.savefig(os.path.join(args.outdir, "trend_%s.png" % mk), dpi=args.dpi); plt.close(fig)

    # LL/xi controlled lines for ell_d and A
    def lines_at(fixed, mk):
        pts = []
        for r in rows:
            if any(r.get(k) is None or abs(r.get(k) - v) > 1e-12 for k, v in fixed.items()):
                continue
            yv = r.get(mk)
            if r.get("LL") is None or yv is None or not np.isfinite(yv):
                continue
            pts.append((r["LL"], yv))
        return sorted(pts)

    refs = [dict(alpha=1e-4, zeta=0.05, Ochi=0.01), dict(alpha=1e-4, zeta=0.1, Ochi=0.005)]
    for mk, ml in [("ell_d", "director coherence ell_d"), ("A", "anchoring A")]:
        fig, axs = plt.subplots(1, len(refs), figsize=(5.2 * len(refs), 4), squeeze=False)
        for ci, rp in enumerate(refs):
            ax = axs[0][ci]
            for xv in (0.3, 1.2):
                fx = dict(rp, xi_align=xv)
                pts = lines_at(fx, mk)
                if len(pts) >= 2:
                    ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", lw=1.7, ms=6, label="xi=%g" % xv)
            ax.set_xscale("log"); ax.set_xlabel("LL"); ax.grid(True, alpha=0.3)
            if mk == "ell_d":
                LLs = np.array(sorted({r["LL"] for r in rows if r.get("LL")}))
                base_pts = lines_at(dict(rp, xi_align=0.3), mk)
                if base_pts and len(LLs) > 1:
                    y0 = base_pts[0][1]; x0 = base_pts[0][0]
                    ax.plot(LLs, y0 * np.sqrt(LLs / x0), "k--", alpha=0.5, label="~sqrt(LL)")
            if mk == "A":
                ax.axhline(0, color="0.6", lw=0.8); ax.set_ylim(-1.05, 1.05)
            if ci == 0:
                ax.set_ylabel(ml)
            ax.set_title("zeta=%g, Ochi=%g" % (rp["zeta"], rp["Ochi"]), fontsize=9); ax.legend(fontsize=8)
        fig.suptitle("%s vs LL (lines: xi)" % ml, fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        fig.savefig(os.path.join(args.outdir, "llxi_%s.png" % mk), dpi=args.dpi); plt.close(fig)

    # demo: director field + chi=0.5 interface + anchoring coloring
    dpath = os.path.join(args.studydir, args.demo)
    if os.path.isdir(dpath):
        ar = loadarchive(dpath)
        fr = ar.read_frame(int((ar.nsteps - ar.nstart) / ar.ninfo))
        chi, phi = grid(fr, "chi"), grid(fr, "phi")
        theta, S, qxx, qyx = director_field(fr)
        A, delta, main = anchoring(fr, chi, phi)
        fig, ax = plt.subplots(figsize=(6.5, 6))
        disp = np.ma.masked_where(phi <= 0.5, chi)
        cmap = plt.get_cmap("bwr").copy(); cmap.set_bad("0.85")
        ax.imshow(disp.T, origin="lower", cmap=cmap, vmin=0, vmax=1, alpha=0.4, interpolation="nearest")
        st = 6
        xs = np.arange(0, theta.shape[0], st); ys = np.arange(0, theta.shape[1], st)
        XX, YY = np.meshgrid(xs, ys, indexing="ij")
        inside = phi[XX, YY] > 0.5
        u = np.cos(theta[XX, YY]); v = np.sin(theta[XX, YY])
        ax.quiver(XX[inside], YY[inside], u[inside], v[inside], pivot="middle",
                  headwidth=0, headlength=0, headaxislength=0, scale=35, width=0.003, color="k", alpha=0.7)
        if main is not None:
            ax.plot(main[:, 0], main[:, 1], "g-", lw=2, label="chi=0.5 interface")
        cy = phi.shape[1]
        ax.set_aspect("equal"); ax.set_title("director + go-grow interface  (A=%.2f, +1=parallel/-1=perp)" % A)
        ax.legend(loc="upper right", fontsize=9)
        fig.tight_layout(); fig.savefig(os.path.join(args.outdir, "demo_director.png"), dpi=args.dpi); plt.close(fig)

    print("\nprocessed %d cases. wrote to %s" % (len(rows), args.outdir))


if __name__ == "__main__":
    main()
