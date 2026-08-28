"""Merge the per-case cw_scan results and draw the campaign figures.

PER CASE (--per-case): one six-panel sheet -- the spatial correlations with their lengths, the
lagged P->m and m->chi correlations against the tau values the runcard SET, the chi
spatio-temporal map C(r,tau), the loop's time series, and its detrended spectrum.

ACROSS THE CAMPAIGN:
  loop_lag     measured lag vs set tau, for both links.  This is the first thing to look at:
               if the measured P->m lag does not track tau_m, nothing downstream means what
               it is supposed to.  The k=0 estimate is the honest one -- the Eulerian one is
               also capped by the advection time and is plotted beside it to show that.
  collapse     the phenotype domain size L_chi against tau_m/t_eddy and tau_chi/t_eddy.  If the
               physics is set by the RATIOS, the three activities fall on one curve here; if
               it is set by absolute times, they do not.  That is the campaign's question.
  feedback     A_eff/A, how far the closed loop pushes the activity away from its nominal
               value, against the same ratios and split by switch sign.  The dbio batch saw
               four different Dbio converge onto one zeta_eff at the instability threshold;
               with three activities this either reproduces as self-organised criticality
               (all three landing on the same A_eff) or it does not.

Usage:
  cw_summary.py --merge <percase_dir> --out cw_scan.json
  cw_summary.py --json cw_scan.json --figdir <dir> [--per-case]
"""
import argparse
import glob
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

ACOL = {1: "C0", 3: "C1", 9: "C3"}
SMARK = {1: "o", -1: "s"}
DFILL = {0.0: "none", 0.03: "full"}


def merge(percase_dir, out):
    rows = []
    for p in sorted(glob.glob(os.path.join(percase_dir, "*.json"))):
        with open(p) as f:
            d = json.load(f)
        rows.extend(d if isinstance(d, list) else [d])
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as f:
        json.dump(rows, f)
    print(f"merged {len(rows)} cases -> {out}")
    return rows


def style(r):
    return dict(color=ACOL.get(int(r.get("A", 0)), "k"),
                marker=SMARK.get(int(r.get("switch_sign", 1)), "^"),
                fillstyle=DFILL.get(float(r.get("Dbio", 0)), "full"),
                ls="none", ms=7, mew=1.4)


def legend_key(ax):
    h = [plt.Line2D([], [], color=c, marker="o", ls="none", label=f"A = {a}")
         for a, c in ACOL.items()]
    h += [plt.Line2D([], [], color="k", marker=m, ls="none",
                     label=f"switch-sign {s:+d}") for s, m in SMARK.items()]
    h += [plt.Line2D([], [], color="k", marker="o", ls="none", fillstyle=f,
                     label=f"Dbio = {d:g}") for d, f in DFILL.items()]
    ax.legend(handles=h, fontsize=8, loc="best", framealpha=.9)


def fig_loop_lag(rows, figdir):
    """Measured lag against the tau the runcard set, for both links of the loop.

    A measured lag of 0 does not mean "instantaneous", it means "shorter than one output
    interval", so those points are drawn at that case's own dt and the resolution floor is
    shaded -- a point sitting in the shaded band carries no information about tau, and
    pretending otherwise by plotting it at zero on a log axis is what collapsed this figure
    the first time it was drawn.
    """
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5), dpi=120)
    for ax, (mk, sk, tk, ttl) in zip(axes, [
            ("lag_Pm_k0", "Pm", "tau_m_set", r"memory:  $P \to m$   (expect $\tau_m$)"),
            ("lag_mchi_k0", "mc", "tau_chi_set", r"phenotype:  $m \to \chi$   (expect $\tau_\chi$)")]):
        xs, ys, dts = [], [], []
        for r in rows:
            x, dt = r.get(tk), r.get("dt") or 1.0
            if x is None or not np.isfinite(x):
                continue
            dts.append(dt)
            for key, kw in ((mk, style(r)),
                            (None, dict(marker="x", color=style(r)["color"], ls="none",
                                        ms=5, alpha=.45))):
                y = r.get(key) if key else \
                    r.get("temporal", {}).get(sk, {}).get("euler_peak_lag")
                if y is None or not np.isfinite(y):
                    continue
                yp = max(y, dt)                      # 0 means "below one output interval"
                ax.plot([x], [yp], **kw)
                if key:
                    xs.append(x); ys.append(yp)
        if xs:
            lo = 0.5 * min([v for v in xs + ys if v > 0] or [1.0])
            hi = 2.0 * max(xs + ys)
            ax.plot([lo, hi], [lo, hi], "k--", lw=1)
            ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
            if dts:
                ax.axhspan(lo, max(dts), color="0.85", zorder=0)
                ax.text(lo * 1.15, max(dts) * 0.8, "below one output interval",
                        fontsize=7, va="top", color="0.35")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel(f"{tk}  (steps, set in the runcard)")
        ax.set_ylabel("measured peak lag  (steps)")
        ax.set_title(ttl, fontsize=11)
        ax.grid(alpha=.3, which="both")
    axes[0].text(.02, .97, "filled/open = k=0 (box average)\nfaint x = Eulerian, point-wise\n"
                           "dashed = measured equals set",
                 transform=axes[0].transAxes, va="top", fontsize=8)
    legend_key(axes[1])
    fig.tight_layout()
    out = os.path.join(figdir, "loop_lag.png")
    fig.savefig(out); plt.close(fig)
    print("wrote", out)


def fig_collapse(rows, figdir):
    """Does the phenotype domain size depend on the RATIOS or on absolute times?"""
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), dpi=120)
    panels = [("tau_m_set", "t_eddy", r"$\tau_m/t_{\rm eddy}$", "L_chi", r"$L_\chi$"),
              ("tau_chi_set", "t_eddy", r"$\tau_\chi/t_{\rm eddy}$", "L_chi", r"$L_\chi$"),
              ("tau_m_set", "t_eddy", r"$\tau_m/t_{\rm eddy}$", "L_m", r"$L_m$"),
              ("tau_chi_set", "t_eddy", r"$\tau_\chi/t_{\rm eddy}$", "chi_bar",
               r"$\langle\chi\rangle$")]
    for ax, (nk, dk, xl, yk, yl) in zip(axes.ravel(), panels):
        for r in rows:
            den = r.get(dk)
            x, y = r.get(nk), r.get(yk)
            if not all(v is not None and np.isfinite(v) for v in (den, x, y)) or den <= 0:
                continue
            ax.plot([x / den], [y], **style(r))
        ax.set_xscale("log"); ax.set_xlabel(xl); ax.set_ylabel(yl); ax.grid(alpha=.3)
    axes[0, 0].set_title("phenotype domain size vs the memory ratio", fontsize=10)
    axes[0, 1].set_title("phenotype domain size vs the phenotype ratio", fontsize=10)
    legend_key(axes[0, 0])
    fig.suptitle("Does the structure collapse on the timescale ratios?", fontsize=12)
    fig.tight_layout()
    out = os.path.join(figdir, "collapse.png")
    fig.savefig(out); plt.close(fig)
    print("wrote", out)


def fig_feedback(rows, figdir):
    """How far the closed loop moves the activity, and whether it converges on one value."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.6), dpi=120)
    for ax, (xk, xl) in zip(axes[:2], [("tau_m_set", r"$\tau_m/t_{\rm eddy}$"),
                                       ("tau_chi_set", r"$\tau_\chi/t_{\rm eddy}$")]):
        for r in rows:
            te, x, ae, A = r.get("t_eddy"), r.get(xk), r.get("A_eff"), r.get("A")
            if not all(v is not None and np.isfinite(v) for v in (te, x, ae, A)) or not A:
                continue
            ax.plot([x / te], [ae / A], **style(r))
        ax.axhline(1.0, ls=":", c="k", lw=.9)
        ax.set_xscale("log"); ax.set_xlabel(xl)
        ax.set_ylabel(r"$A_{\rm eff}/A$"); ax.grid(alpha=.3)
    ax = axes[2]
    for r in rows:
        A, ae = r.get("A"), r.get("A_eff")
        if not all(v is not None and np.isfinite(v) for v in (A, ae)):
            continue
        ax.plot([A], [ae], **style(r))
    lim = [0.5, 20]
    ax.plot(lim, lim, "k--", lw=1)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim)
    ax.set_xlabel("A (nominal)"); ax.set_ylabel(r"$A_{\rm eff}$ (closed loop)")
    ax.set_title("a horizontal band here = self-organised criticality", fontsize=9)
    ax.grid(alpha=.3, which="both")
    legend_key(axes[0])
    fig.suptitle("Closed-loop feedback on the activity", fontsize=12)
    fig.tight_layout()
    out = os.path.join(figdir, "feedback.png")
    fig.savefig(out); plt.close(fig)
    print("wrote", out)


def fig_case(r, figdir):
    """The six-panel per-case sheet."""
    fig = plt.figure(figsize=(15, 8.5), dpi=115)
    gs = fig.add_gridspec(2, 3, hspace=.32, wspace=.28)
    rr = np.asarray(r["r"], float)
    sp = r["spatial"]

    ax = fig.add_subplot(gs[0, 0])
    for k, lab in (("uu", r"$C_{uu}$"), ("PP", r"$C_{PP}$"),
                   ("mm", r"$C_{mm}$"), ("cc", r"$C_{\chi\chi}$")):
        c = np.asarray(sp[k]["C"], float)
        ax.plot(rr[:len(c)], c, label=f"{lab}  $\\ell$={sp[k]['length']:.0f}")
    ax.axhline(0, c="k", lw=.6); ax.set_xlim(0, min(rr[-1], 150))
    ax.set_xlabel("r"); ax.set_ylabel("C(r)"); ax.legend(fontsize=8); ax.grid(alpha=.3)
    ax.set_title("1. spatial autocorrelations", fontsize=10)

    ax = fig.add_subplot(gs[0, 1])
    for k, lab in (("Pm", r"$C_{Pm}$"), ("mc", r"$C_{m\chi}$"), ("Pc", r"$C_{P\chi}$")):
        c = np.asarray(sp[k]["C"], float)
        ax.plot(rr[:len(c)], c, label=lab)
    ax.axhline(0, c="k", lw=.6); ax.set_xlim(0, min(rr[-1], 150))
    ax.set_xlabel("r"); ax.set_ylabel("C(r)"); ax.legend(fontsize=8); ax.grid(alpha=.3)
    ax.set_title("2-3. spatial cross-correlations", fontsize=10)

    ax = fig.add_subplot(gs[0, 2])
    tp = r["temporal"]; lags = np.asarray(tp["lags"], float)
    for k, lab, c in (("Pm", r"$P\to m$", "C0"), ("mc", r"$m\to\chi$", "C1")):
        if "k0" in tp[k]:
            y = np.asarray(tp[k]["k0"], float)
            ax.plot(lags[:len(y)], y, c=c, label=f"{lab} (k=0)")
        y = np.asarray(tp[k]["euler"], float)
        ax.plot(lags[:len(y)], y, c=c, ls=":", alpha=.6, label=f"{lab} (Eulerian)")
    for tau, c, lab in ((r["tau_m_set"], "C0", r"$\tau_m$"),
                        (r["tau_chi_set"], "C1", r"$\tau_\chi$")):
        if tau and np.isfinite(tau) and tau < lags[-1]:
            ax.axvline(tau, c=c, ls="--", lw=1)
            ax.text(tau, .95, lab, color=c, fontsize=8, ha="center", va="top",
                    transform=ax.get_xaxis_transform())
    ax.axhline(0, c="k", lw=.6)
    ax.set_xlabel(r"lag $\tau$ (steps)"); ax.set_ylabel(r"$C(\tau)$")
    ax.legend(fontsize=7); ax.grid(alpha=.3)
    ax.set_title("2-3. lagged cross-correlations", fontsize=10)

    ax = fig.add_subplot(gs[1, 0])
    st = np.asarray(r["chi_st"], float)
    ex = [0, rr[min(len(rr), st.shape[1]) - 1], lags[0], lags[min(len(lags), st.shape[0]) - 1]]
    im = ax.imshow(st, origin="lower", aspect="auto", extent=ex, cmap="RdBu_r",
                   vmin=-1, vmax=1)
    fig.colorbar(im, ax=ax, fraction=.046)
    ax.set_xlim(0, min(rr[-1], 150))
    ax.set_xlabel("r"); ax.set_ylabel(r"lag $\tau$")
    ax.set_title(r"4. $C_{\chi\chi}(r,\tau)$", fontsize=10)

    ax = fig.add_subplot(gs[1, 1])
    S = r["series"]; t = np.asarray(S["t"], float)
    ax.plot(t, S["chibar"], label=r"$\langle\chi\rangle$")
    ax.plot(t, S["mbar"], label=r"$\langle m\rangle$")
    ax2 = ax.twinx(); ax2.plot(t, S["urms"], c="C2", lw=.9, label=r"$u_{\rm rms}$")
    ax2.set_ylabel(r"$u_{\rm rms}$", color="C2")
    ax.axhline(.5, ls=":", c="k", lw=.8)
    ax.set_xlabel("t"); ax.set_ylabel("loop state"); ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=.3); ax.set_title("loop time series", fontsize=10)

    ax = fig.add_subplot(gs[1, 2])
    for k in ("chibar", "mbar", "urms"):
        s = r.get("spectra", {}).get(k)
        if s:
            ax.loglog(s["f"], s["P"], lw=.9, label=f"{k}  peak T={s['period_peak']:.0f}")
    ax.set_xlabel("frequency (1/step)"); ax.set_ylabel("power")
    ax.legend(fontsize=7); ax.grid(alpha=.3, which="both")
    ax.set_title("detrended spectra (a limit cycle shows here)", fontsize=10)

    fig.suptitle(f"{r['case']}    A={r.get('A')}  s={r.get('switch_sign'):+d}  "
                 f"Dbio={r.get('Dbio'):g}  tau_m={r['tau_m_set']:.0f}  "
                 f"tau_chi={r['tau_chi_set']:.0f}  t_eddy={r['t_eddy']:.0f}  "
                 f"A_eff={r['A_eff']:.2f}", fontsize=12)
    out = os.path.join(figdir, f"case_{r['case']}.png")
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--merge", default=None)
    ap.add_argument("--json", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--figdir", default=None)
    ap.add_argument("--per-case", action="store_true")
    a = ap.parse_args()

    if a.merge:
        rows = merge(a.merge, a.out or os.path.join(a.merge, "..", "cw_scan.json"))
    else:
        with open(a.json) as f:
            rows = json.load(f)
    if not a.figdir:
        return
    os.makedirs(a.figdir, exist_ok=True)
    fig_loop_lag(rows, a.figdir)
    fig_collapse(rows, a.figdir)
    fig_feedback(rows, a.figdir)
    if a.per_case:
        for r in rows:
            try:
                print("wrote", fig_case(r, a.figdir))
            except Exception as exc:
                print(f"  case figure {r.get('case')} FAILED: {type(exc).__name__}: {exc}")


if __name__ == "__main__":
    main()
