#!/usr/bin/env python3
"""Presentation figures: one panel per file, no parameter text, nothing but the plot.

The per-case sheet exists to let one person diagnose one run, so it carries panel numbering,
a parameter suptitle and six panels at ~500 px each. None of that survives being put on a
slide. These are the same curves replotted large and bare.

Everything comes from cw_scan.json -- the correlation arrays and the full C(r,tau) map are
already in it -- so this never touches the 554 GB of frames and runs in seconds per case.

Emits, per case:
    <case>_crosscorr.png    spatial cross-correlations C_Pm, C_m,chi, C_P,chi
    <case>_lagcorr.png      lagged cross-correlations, k=0 and Eulerian, with tau markers
    <case>_st.png           the chi spatio-temporal autocorrelation C(r, tau)
    <case>_autocorr.png     spatial autocorrelations (added because a cross-correlation is
                            hard to read without the scales its own fields set)

Usage: cw_figs.py <cw_scan.json> <outdir> [--cases ...] [--dpi 200] [--label]
`--label` puts a single short parameter line back on, for when a figure travels alone.
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

FS = 15          # axis label size
TS = 13          # tick size


def bare(ax, xlab, ylab):
    ax.set_xlabel(xlab, fontsize=FS)
    ax.set_ylabel(ylab, fontsize=FS)
    ax.tick_params(labelsize=TS)
    ax.grid(alpha=.25)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def short(r):
    return (r"$\zeta$ = %g   $\tau_m/\tau_{\rm mot}$ = %.2f   "
            r"$\tau_\chi/\tau_{\rm mot}$ = %.2f   s = %+d   $D_{\rm bio}$ = %g"
            % (2 * r["A"] * r["CC"], r["tau_m_set"] / r["t_eddy"],
               r["tau_chi_set"] / r["t_eddy"], r["switch_sign"], r["Dbio"]))


def maybe_label(fig, r, label):
    if label:
        fig.text(0.5, 0.005, short(r), ha="center", fontsize=11, color="0.35")


def fig_cross(r, out, dpi, label):
    rr = np.asarray(r["r"], float)
    fig, ax = plt.subplots(figsize=(7.2, 5.4), dpi=dpi)
    for k, lab, c in (("Pm", r"$C_{Pm}$", "C0"),
                      ("mc", r"$C_{m\chi}$", "C1"),
                      ("Pc", r"$C_{P\chi}$", "C2")):
        y = np.asarray(r["spatial"][k]["C"], float)
        ax.plot(rr[:len(y)], y, c=c, lw=2.4, label=lab)
    ax.axhline(0, c="0.4", lw=.9)
    ax.set_xlim(0, min(rr[-1], 120))
    bare(ax, r"$r$", r"$C(r)$")
    ax.legend(fontsize=FS, frameon=False)
    maybe_label(fig, r, label)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)


def fig_auto(r, out, dpi, label):
    rr = np.asarray(r["r"], float)
    fig, ax = plt.subplots(figsize=(7.2, 5.4), dpi=dpi)
    for k, lab, c in (("uu", r"$C_{uu}$", "C0"), ("PP", r"$C_{PP}$", "C3"),
                      ("mm", r"$C_{mm}$", "C1"), ("cc", r"$C_{\chi\chi}$", "C2")):
        y = np.asarray(r["spatial"][k]["C"], float)
        ax.plot(rr[:len(y)], y, c=c, lw=2.4, label=lab)
    ax.axhline(0, c="0.4", lw=.9)
    ax.set_xlim(0, min(rr[-1], 120))
    bare(ax, r"$r$", r"$C(r)$")
    ax.legend(fontsize=FS, frameon=False)
    maybe_label(fig, r, label)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)


def fig_lag(r, out, dpi, label, eulerian=False):
    """Only the k=0 curves by default.

    An unexplained dotted line is worse than an absent one, and the Eulerian estimator needs a
    paragraph to justify (it also decorrelates on the advection time, so it under-reads the
    lag whenever tau exceeds t_eddy).  On a slide that paragraph is not available.
    """
    tp = r["temporal"]
    lags = np.asarray(tp["lags"], float)
    fig, ax = plt.subplots(figsize=(7.6, 5.4), dpi=dpi)
    for k, lab, c in (("Pm", r"$P \to m$", "C0"), ("mc", r"$m \to \chi$", "C1")):
        if "k0" in tp[k]:
            y = np.asarray(tp[k]["k0"], float)
            ax.plot(lags[:len(y)], y, c=c, lw=2.6, label=lab)
        if eulerian:
            y = np.asarray(tp[k]["euler"], float)
            ax.plot(lags[:len(y)], y, c=c, lw=1.4, ls=":", alpha=.8,
                    label=lab + " (Eulerian)")
    seen = set()
    for tau, c, lab in ((r["tau_m_set"], "C0", r"$\tau_m$"),
                        (r["tau_chi_set"], "C1", r"$\tau_\chi$")):
        if tau and np.isfinite(tau) and tau < lags[-1] and round(tau) not in seen:
            seen.add(round(tau))
            ax.axvline(tau, c=c, ls="--", lw=1.6, alpha=.75)
            ax.text(tau, 1.01, lab, color=c, fontsize=FS, ha="center", va="bottom",
                    transform=ax.get_xaxis_transform())
    ax.axhline(0, c="0.4", lw=.9)
    ax.set_xlim(0, lags[-1])
    bare(ax, r"lag $\tau$", r"$C(\tau)$")
    ax.legend(fontsize=FS, frameon=False)
    maybe_label(fig, r, label)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)


def fig_st(r, out, dpi, label):
    st = np.asarray(r["chi_st"], float)
    rr = np.asarray(r["r"], float)
    lags = np.asarray(r["temporal"]["lags"], float)
    nr = min(len(rr), st.shape[1])
    nt = min(len(lags), st.shape[0])
    # Crop to where the correlation still is. Left whole, this map is ~80% empty: the lag axis
    # runs to 40000 while the 1/e contour closes by 4000, and the eye reads the white, not the
    # ridge. Cropping to a few decay times is what makes the two scales legible at once.
    c0 = st[:nt, 0]
    below = np.where(c0 <= np.exp(-2))[0]
    tmax = lags[min(nt - 1, int(below[0] * 2.0))] if len(below) else lags[nt - 1]
    r0 = st[0, :nr]
    belowr = np.where(r0 <= np.exp(-2))[0]
    rmax = rr[min(nr - 1, int(belowr[0] * 2.0))] if len(belowr) else rr[nr - 1]
    fig, ax = plt.subplots(figsize=(7.0, 5.6), dpi=dpi)
    im = ax.imshow(st[:nt, :nr], origin="lower", aspect="auto", cmap="RdBu_r",
                   vmin=-1, vmax=1, extent=[0, rr[nr - 1], lags[0], lags[nt - 1]])
    cs = ax.contour(np.linspace(0, rr[nr - 1], nr), np.linspace(lags[0], lags[nt - 1], nt),
                    st[:nt, :nr], levels=[np.exp(-1)], colors="k", linewidths=1.4)
    ax.clabel(cs, fmt={np.exp(-1): "1/e"}, fontsize=11)
    ax.set_xlim(0, rmax)
    ax.set_ylim(0, tmax)
    bare(ax, r"$r$", r"lag $\tau$")
    ax.grid(False)
    cb = fig.colorbar(im, ax=ax, fraction=.046, pad=.02)
    cb.set_label(r"$C_{\chi\chi}(r,\tau)$", fontsize=FS)
    cb.ax.tick_params(labelsize=TS)
    maybe_label(fig, r, label)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("scan")
    ap.add_argument("outdir")
    ap.add_argument("--cases", nargs="*", default=None)
    ap.add_argument("--dpi", type=int, default=200)
    ap.add_argument("--label", action="store_true",
                    help="add one short parameter line at the bottom")
    ap.add_argument("--eulerian", action="store_true",
                    help="also draw the Eulerian lagged curves (dotted) on the lag figure")
    a = ap.parse_args()

    with open(a.scan) as f:
        rows = json.load(f)
    rows = rows if isinstance(rows, list) else [rows]
    if a.cases:
        rows = [r for r in rows if r["case"] in a.cases]
    os.makedirs(a.outdir, exist_ok=True)
    n = 0
    for r in rows:
        c = r["case"]
        for tag, fn in (("crosscorr", fig_cross), ("autocorr", fig_auto),
                        ("lagcorr", fig_lag), ("st", fig_st)):
            try:
                kw = {"eulerian": a.eulerian} if tag == "lagcorr" else {}
                fn(r, os.path.join(a.outdir, f"{c}_{tag}.png"), a.dpi, a.label, **kw)
                n += 1
            except Exception as exc:
                print(f"  {c} {tag}: {type(exc).__name__}: {exc}", flush=True)
        print(f"  {c}", flush=True)
    print(f"wrote {n} figures to {a.outdir}")


if __name__ == "__main__":
    main()
