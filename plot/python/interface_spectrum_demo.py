"""Pedagogical figure: how the radial Fourier spectrum of a go-grow interface
yields the two characteristic length scales.

Panels:
  (a) the interface (chi field) with the chi=0.5 contour, centroid, mean-radius
      circle, and radial spokes illustrating r(theta);
  (b) the loop "unrolled" -- radius r(theta) vs angle theta, with <r> and the
      +/- w amplitude band;
  (c) the angular power spectrum P(n), with the spectral-weighted mode <n>
      marked (-> wavelength lam = 2*pi*<r>/<n>) and the total power (-> amplitude
      w = sqrt(sum P), by Parseval).

Usage:
  python interface_spectrum_demo.py <case_dir> <outdir> [--frame N] [--dpi D]
"""
import matplotlib
matplotlib.use("Agg")

import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)
from archive.archive import loadarchive


def grid(frame, field):
    lx = frame.parameters["LX"]; ly = frame.parameters["LY"]
    return np.array(getattr(frame, field)).reshape((lx, ly))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("casedir")
    ap.add_argument("outdir")
    ap.add_argument("--frame", type=int, default=None)
    ap.add_argument("--dpi", type=int, default=160)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    ar = loadarchive(args.casedir)
    idx = args.frame if args.frame is not None else int((ar.nsteps - ar.nstart) / ar.ninfo)
    frame = ar.read_frame(idx)
    chi = grid(frame, "chi"); phi = grid(frame, "phi")

    # main chi=0.5 contour inside the colony
    chim = np.ma.masked_where(phi <= 0.5, chi)
    fig0 = plt.figure(); ax0 = fig0.add_subplot(111)
    segs = list(ax0.contour(chim.T, levels=[0.5]).allsegs[0])
    plt.close(fig0)
    main = max(segs, key=lambda s: np.sqrt((np.diff(s, axis=0) ** 2).sum(axis=1)).sum())
    x, y = main[:, 0], main[:, 1]
    cx, cy = x.mean(), y.mean()
    rr = np.sqrt((x - cx) ** 2 + (y - cy) ** 2)
    th = np.arctan2(y - cy, x - cx)

    # resample r(theta) and take the spectrum
    N = 512
    tg = np.linspace(-np.pi, np.pi, N, endpoint=False)
    order = np.argsort(th)
    rg = np.interp(tg, th[order], rr[order], period=2 * np.pi)
    rmean = rg.mean()
    h = rg - rmean
    w = float(np.sqrt(np.mean(h ** 2)))
    P = np.abs(np.fft.rfft(h)) ** 2
    n = np.arange(len(P))
    nmax = max(6, int(rmean / 2))
    sel = (n >= 2) & (n <= nmax)
    nmean = float(np.sum(n[sel] * P[sel]) / np.sum(P[sel]))
    ndom = int(n[sel][np.argmax(P[sel])])
    lam = 2 * np.pi * rmean / nmean

    cyc = plt.get_cmap("twilight")
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.8))

    # (a) interface
    ax = axs[0]
    disp = np.ma.masked_where(phi <= 0.5, chi)
    cmap = plt.get_cmap("bwr").copy(); cmap.set_bad("0.85")
    ax.imshow(disp.T, origin="lower", cmap=cmap, vmin=0, vmax=1, alpha=0.55, interpolation="nearest")
    ax.scatter(x, y, c=(th + np.pi) / (2 * np.pi), cmap=cyc, s=6, zorder=3)
    tt = np.linspace(0, 2 * np.pi, 240)
    ax.plot(cx + rmean * np.cos(tt), cy + rmean * np.sin(tt), "k--", lw=1.2, label="mean radius <r>")
    ax.plot(cx, cy, "k+", ms=12, mew=2)
    for a in np.linspace(-np.pi, np.pi, 12, endpoint=False):
        ri = np.interp(a, tg, rg)
        ax.plot([cx, cx + ri * np.cos(a)], [cy, cy + ri * np.sin(a)], color=cyc((a + np.pi) / (2 * np.pi)), lw=1.0)
    m = 1.4 * rmean
    ax.set_xlim(cx - m, cx + m); ax.set_ylim(cy - m, cy + m)
    ax.set_aspect("equal"); ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_title("(a) go-grow interface  r(theta) about centroid")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.7)

    # (b) unrolled
    ax = axs[1]
    ax.scatter(tg, rg, c=(tg + np.pi) / (2 * np.pi), cmap=cyc, s=10)
    ax.axhline(rmean, ls="--", color="k", lw=1.2, label="<r> = %.1f" % rmean)
    ax.fill_between(tg, rmean - w, rmean + w, color="gray", alpha=0.2)
    ax.annotate("", xy=(2.6, rmean + w), xytext=(2.6, rmean - w),
                arrowprops=dict(arrowstyle="<->", color="k"))
    ax.text(2.75, rmean, "2w", va="center", fontsize=10)
    ax.set_xlabel("angle theta (loop unrolled)"); ax.set_ylabel("radius r(theta)")
    ax.set_title("(b) unroll the loop: radius vs angle")
    ax.text(0.02, 0.04, "amplitude  w = RMS(r - <r>) = %.1f" % w, transform=ax.transAxes,
            fontsize=10, bbox=dict(facecolor="white", alpha=0.7, edgecolor="0.6"))
    ax.legend(loc="upper right", fontsize=9)

    # (c) spectrum
    ax = axs[2]
    nn = n[1:nmax + 1]
    PP = P[1:nmax + 1]
    bars = ax.bar(nn, PP, color="#1f77b4")
    bars[0].set_color("0.7")   # n=1 = centroid offset, not roughness
    ax.axvline(nmean, color="#d62728", lw=2, label="<n> = %.1f" % nmean)
    ax.axvline(ndom, color="#2ca02c", lw=1.2, ls=":", label="n_dom = %d" % ndom)
    ax.set_xlabel("angular mode  n"); ax.set_ylabel("power  P(n)")
    ax.set_title("(c) Fourier spectrum of r(theta)")
    ax.legend(loc="upper right", fontsize=9)
    ax.text(0.97, 0.62,
            "wavelength\nlam = 2*pi*<r>/<n>\n   = %.0f" % lam,
            transform=ax.transAxes, ha="right", fontsize=10,
            bbox=dict(facecolor="#fff3e0", alpha=0.95, edgecolor="#d62728"))
    ax.text(0.97, 0.40,
            "amplitude\nw = sqrt(sum P(n>=2))\n   = %.1f" % w,
            transform=ax.transAxes, ha="right", fontsize=10,
            bbox=dict(facecolor="#e8f5e9", alpha=0.95, edgecolor="#2ca02c"))

    fig.suptitle("Radial Fourier analysis of the go-grow interface  —  %s  (frame %d)"
                 % (os.path.basename(os.path.normpath(args.casedir)), idx), fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = os.path.join(args.outdir, "spectrum_demo.png")
    fig.savefig(out, dpi=args.dpi); plt.close(fig)
    print("rmean=%.1f w=%.2f <n>=%.2f n_dom=%d lam=%.1f" % (rmean, w, nmean, ndom, lam))
    print("wrote", out)


if __name__ == "__main__":
    main()
