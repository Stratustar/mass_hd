"""Self-affine (KPZ-type) roughness analysis of the go-grow interface.

For each variant (final frame): extract the chi=0.5 contour inside the colony,
form the radial fluctuation h(theta)=r(theta)-<r>, and its angular structure
factor S(n)=|h_n|^2. Then:
  W           = interface width = RMS(h)  (and W_rel = W/<r>)
  alpha_rough = roughness exponent from the PSD power law S(n) ~ n^-(1+2 alpha)
                over n in [2, ncut]:  slope = d logS/d logn,  alpha = -(slope+1)/2
  R2          = goodness of the power-law fit (how self-affine)
  peak_ratio  = max_n S(n) / S_fit(n)  (excess of the strongest mode over the
                power-law background); has_peak if peak_ratio > THRESH
  lam_peak    = 2*pi*<r>/n_peak  (reported only where has_peak)

Outputs CSV, Spearman heatmap, controlled-slice trends for alpha_rough and W,
and a spectra demo (S(n) with fits) for representative cases.

Usage: python interface_kpz.py <scratch_study_dir> <outdir> [--pattern G]
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
from interface_roughness import grid, parse_params

PEAK_THRESH = 4.0
COLS = ["variant", "alpha", "zeta", "Ochi", "LL", "xi_align", "R_eff", "W", "W_rel",
        "alpha_rough", "R2", "peak_ratio", "n_peak", "lam_peak", "has_peak", "n_pts"]


def spectrum(chi, phi, N=512):
    """Return (rmean, S(n), h, ncut) for the main chi=0.5 contour, or None."""
    chim = np.ma.masked_where(phi <= 0.5, chi)
    fig = plt.figure(); ax = fig.add_subplot(111)
    segs = list(ax.contour(chim.T, levels=[0.5]).allsegs[0])
    plt.close(fig)
    if not segs:
        return None
    main = max(segs, key=lambda s: np.sqrt((np.diff(s, axis=0) ** 2).sum(axis=1)).sum())
    if len(main) < 12:
        return None
    x, y = main[:, 0], main[:, 1]
    cx, cy = x.mean(), y.mean()
    rr = np.sqrt((x - cx) ** 2 + (y - cy) ** 2)
    th = np.arctan2(y - cy, x - cx)
    tg = np.linspace(-np.pi, np.pi, N, endpoint=False)
    order = np.argsort(th)
    rg = np.interp(tg, th[order], rr[order], period=2 * np.pi)
    rmean = float(rg.mean())
    h = rg - rmean
    S = np.abs(np.fft.rfft(h)) ** 2
    ncut = int(max(8, min(rmean / 3.0, len(S) - 1)))
    return rmean, S, h, ncut, len(main)


def averaged_spectrum(ar, N=512, K=12):
    """Average S(n), <r>, W over the last K frames to cut single-frame PSD noise."""
    last = int((ar.nsteps - ar.nstart) / ar.ninfo)
    Ssum = None; nS = 0; rmeans = []; Ws = []; npts = 0
    for fi in range(max(0, last - K + 1), last + 1):
        try:
            fr = ar.read_frame(fi)
            sp = spectrum(grid(fr, "chi"), grid(fr, "phi"), N=N)
        except Exception:
            sp = None
        if sp is None:
            continue
        rmean, S, h, ncut, np_ = sp
        Ssum = S if Ssum is None else Ssum + S
        nS += 1; rmeans.append(rmean); Ws.append(np.sqrt(np.mean(h ** 2))); npts = np_
    if nS == 0:
        return None
    rm = float(np.mean(rmeans))
    ncut = int(max(8, min(rm / 3.0, len(Ssum) - 1)))
    return rm, Ssum / nS, float(np.mean(Ws)), ncut, npts, nS


def logbin(n, S, nmin=2, nmax=None, nbins=14):
    nmax = nmax or n.max()
    edges = np.logspace(np.log10(nmin), np.log10(nmax), nbins + 1)
    nb, Sb = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (n >= a) & (n < b) & (S > 0)
        if sel.sum() > 0:
            nb.append(np.exp(np.mean(np.log(n[sel])))); Sb.append(np.exp(np.mean(np.log(S[sel]))))
    return np.array(nb), np.array(Sb)


def kpz_metrics(ar):
    sp = averaged_spectrum(ar)
    if sp is None:
        return None
    rmean, S, W, ncut, npts, nframes = sp
    n = np.arange(len(S))
    nb, Sb = logbin(n, S, 2, ncut)
    if len(nb) < 4:
        return None
    slope, intc = np.polyfit(np.log(nb), np.log(Sb), 1)
    fit = np.exp(intc) * nb ** slope
    ss_res = np.sum((np.log(Sb) - np.log(fit)) ** 2)
    ss_tot = np.sum((np.log(Sb) - np.log(Sb).mean()) ** 2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    alpha_rough = -(slope + 1) / 2.0
    sel = (n >= 2) & (n <= ncut) & (S > 0)
    ratio = S[sel] / (np.exp(intc) * n[sel] ** slope)
    j = int(np.argmax(ratio))
    peak_ratio = float(ratio[j]); n_peak = int(n[sel][j])
    has_peak = peak_ratio > PEAK_THRESH
    lam_peak = (2 * np.pi * rmean / n_peak) if has_peak else np.nan
    return dict(R_eff=rmean, W=W, W_rel=W / rmean, alpha_rough=alpha_rough, R2=R2,
                peak_ratio=peak_ratio, n_peak=n_peak, lam_peak=lam_peak,
                has_peak=int(has_peak), n_pts=npts)


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
            m = kpz_metrics(ar)
        except Exception as e:
            print("FAIL", name, repr(e))
        pp = parse_params(name)
        row = dict(variant=name, alpha=pp.get("alpha"), zeta=pp.get("zeta"),
                   Ochi=pp.get("Ochi"), LL=pp.get("LL"), xi_align=pp.get("xi"))
        if m:
            row.update(m)
        rows.append(row)
        print("ok" if m else "fail", name)

    with open(os.path.join(args.outdir, "interface_kpz.csv"), "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=COLS); wr.writeheader()
        for r in rows:
            wr.writerow({k: r.get(k) for k in COLS})

    base = lambda r: (r.get("LL") is not None and abs(r["LL"] - 0.002) < 1e-9
                      and r.get("xi_align") is not None and abs(r["xi_align"] - 0.3) < 1e-9)
    MAIN = [r for r in rows if base(r)]
    params = ("alpha", "zeta", "Ochi")
    npeak = sum(1 for r in rows if r.get("has_peak") == 1)
    aok = col(rows, "alpha_rough"); aok = aok[np.isfinite(aok)]
    r2 = col(rows, "R2"); r2 = r2[np.isfinite(r2)]
    lines = ["KPZ roughness: n=%d cases, %d with a clear spectral peak (peak_ratio>%.0f)."
             % (len(rows), npeak, PEAK_THRESH),
             "alpha_rough: mean=%.2f median=%.2f range[%.2f,%.2f]; power-law fit R2 mean=%.2f"
             % (aok.mean(), np.median(aok), aok.min(), aok.max(), r2.mean())]

    # Spearman heatmap (main grid)
    mets = ["alpha_rough", "W_rel", "peak_ratio", "R2"]
    rho = np.full((len(mets), 3), np.nan)
    lines.append("Spearman rho (main grid, n=%d):" % len(MAIN))
    for i, mk in enumerate(mets):
        mv = col(MAIN, mk); cells = []
        for j, pk in enumerate(params):
            pv = col(MAIN, pk); ok = np.isfinite(mv) & np.isfinite(pv)
            if ok.sum() >= 5:
                rr, pp = spearmanr(pv[ok], mv[ok]); rho[i, j] = rr
                cells.append("%s %+.2f(p=%.2g)" % (pk, rr, pp))
        lines.append("  %-11s " % mk + "  ".join(cells))
    summary = "\n".join(lines)
    print("\n" + summary)
    with open(os.path.join(args.outdir, "kpz_summary.txt"), "w") as fh:
        fh.write(summary + "\n")

    fig, ax = plt.subplots(figsize=(4.8, 3.4))
    im = ax.imshow(rho, cmap="coolwarm", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(3)); ax.set_xticklabels([r"$\alpha$", r"$\zeta$", r"$O_\chi$"])
    ax.set_yticks(range(len(mets))); ax.set_yticklabels([r"$\alpha_{rough}$", r"$W/\langle r\rangle$", "peak ratio", r"$R^2$"])
    for i in range(len(mets)):
        for j in range(3):
            if np.isfinite(rho[i, j]):
                ax.text(j, i, "%.2f" % rho[i, j], ha="center", va="center",
                        color="white" if abs(rho[i, j]) > 0.55 else "black")
    ax.set_title("Spearman: KPZ metrics vs parameter")
    fig.colorbar(im, ax=ax, fraction=0.046); fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "corr_heatmap.png"), dpi=args.dpi); plt.close(fig)

    # controlled-slice trends for alpha_rough and W_rel
    P = {p: col(MAIN, p) for p in params}; pv2 = {p: lohi(P[p]) for p in params}
    colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd"]
    labels = {"alpha_rough": r"roughness exponent $\alpha_{rough}$", "W_rel": r"interface width $W/\langle r\rangle$"}
    pl = {"alpha": r"$\alpha$ (proliferation)", "zeta": r"$\zeta$ (activity)", "Ochi": r"$O_\chi$"}
    for mk in ["alpha_rough", "W_rel"]:
        y = col(MAIN, mk)
        fig, axs = plt.subplots(1, 3, figsize=(13, 3.6), squeeze=False)
        for cj, xp in enumerate(params):
            ax = axs[0][cj]; others = [p for p in params if p != xp]
            xs = sorted({v for v in P[xp] if np.isfinite(v)}); ci = 0
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
            ax.set_xscale("log"); ax.set_xlabel(pl[xp]); ax.grid(True, alpha=0.3)
            if cj == 0:
                ax.set_ylabel(labels[mk])
            ax.legend(fontsize=7, framealpha=0.6)
        fig.suptitle(labels[mk] + " vs parameters (representative slices)", fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        fig.savefig(os.path.join(args.outdir, "trend_%s.png" % mk), dpi=args.dpi); plt.close(fig)

    # spectra demo: S(n) with power-law fits for representative cases (peaked vs power-law)
    demo = ["a0p00005_z0p01_O0p02_LL0p002_xi0p3", "a0p0002_z0p01_O0p01_LL0p002_xi0p3",
            "a0p0001_z0p1_O0p01_LL0p002_xi0p3", "a0p0001_z0p05_O0p01_LL0p032_xi0p3"]
    fig, ax = plt.subplots(figsize=(7, 5))
    cols = ["#1f77b4", "#2ca02c", "#d62728", "#9467bd"]
    for k, name in enumerate(demo):
        dpath = os.path.join(args.studydir, name)
        if not os.path.isdir(dpath):
            continue
        try:
            ar = loadarchive(dpath)
            sp = averaged_spectrum(ar)
            mm = kpz_metrics(ar)
        except Exception:
            continue
        if sp is None or mm is None:
            continue
        rmean, S, W, ncut, _, _ = sp
        n = np.arange(len(S))
        nb, Sb = logbin(n, S, 2, ncut)
        lab = name.split("_LL")[0].replace("a0p", r"$\alpha$=0.").replace("_z0p", r", $\zeta$=0.").replace("_O0p", r", $O$=0.")
        ax.loglog(nb, Sb, "o", ms=5, color=cols[k], alpha=0.85,
                  label=lab + (r"  ($\alpha_r$=%.2f%s)" % (mm["alpha_rough"], ", PEAK" if mm["has_peak"] else "")))
        sl, ic = np.polyfit(np.log(nb), np.log(Sb), 1)
        ax.loglog(nb, np.exp(ic) * nb ** sl, "-", color=cols[k], lw=1.4)
        if mm["has_peak"] and np.isfinite(mm["n_peak"]):
            ax.axvline(mm["n_peak"], color=cols[k], ls=":", lw=1, alpha=0.6)
    ax.set_xlabel(r"angular mode  $n$"); ax.set_ylabel(r"structure factor  $S(n)$")
    ax.set_title(r"Interface structure factor: power law (slope $\to\alpha_r$) + peaks")
    ax.legend(fontsize=7, loc="lower left")
    fig.tight_layout(); fig.savefig(os.path.join(args.outdir, "spectra_demo.png"), dpi=args.dpi); plt.close(fig)

    print("\nprocessed %d cases. wrote to %s" % (len(rows), args.outdir))


if __name__ == "__main__":
    main()
