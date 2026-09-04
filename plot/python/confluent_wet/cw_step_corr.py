#!/usr/bin/env python3
"""Lagged correlations of P, m and chi, in three stages of the run.

    cw_step_corr.py <case_output_dir>... --out DIR [--tau-c 225.6] [--stride 4]

WHAT IS CORRELATED. The model is a cascade, P -> m -> chi, so the informative objects are
the CROSS-correlations between successive links and across the whole chain:

    C_XY(lag) = < dX(t) dY(t+lag) >  / sqrt(var dX var dY)

with d denoting the departure from that frame's spatial mean, averaged over pixels and over
t. A positive lag means Y follows X. The peak position is the transport delay of one link
and should sit near tau_m for P -> m and near tau_chi for m -> chi.

chi's temporal AUTO-correlation is computed alongside, because the memory time it defines is
what the fate labels are ultimately arguing about, and its 1/e crossing is reported.

THREE STAGES, which is the point of doing it this way. The run is cut into three equal
windows -- early, middle, late -- and every correlation is computed separately in each. A
transient that has not finished shows up as a correlation that moves between the windows,
which no single whole-run correlation can reveal. The windows are a third of the record each
so the maximum lag can still reach ~1/4 of a window, i.e. tens of tau_c, which is what a
tau_m of 22 tau_c needs.

The data is the video stream (0.2 tau_c per frame, dequantised), not the full-resolution
frames: 26 frames 17.7 tau_c apart cannot resolve a lag at all. Pixels are subsampled by
--stride, since neighbouring stream pixels are already block averages and are not
independent.
"""
import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cw_common as cw                   # noqa: E402
import cw_stream                         # noqa: E402

PAIRS = [("P", "m"), ("m", "chi"), ("P", "chi")]
STAGES = ["early", "middle", "late"]
MAXLAG_FRAC = 0.25          # of the window


def load_window(st, i0, i1, stride):
    """{field: (nt, npix)} dequantised and spatial-mean removed, for one window."""
    out = {}
    for name in ("P", "m", "chi"):
        a = st.dequant(name, slice(i0, i1))[:, ::stride, ::stride]
        a = a.reshape(a.shape[0], -1).astype(np.float64)
        out[name] = a - a.mean(axis=1, keepdims=True)
    return out


def xcorr(X, Y, maxlag):
    """C(lag) for lag = -maxlag..maxlag; positive lag means Y follows X."""
    n = X.shape[0]
    sx = np.sqrt((X * X).mean())
    sy = np.sqrt((Y * Y).mean())
    if sx <= 0 or sy <= 0:
        return np.full(2 * maxlag + 1, np.nan)
    out = np.empty(2 * maxlag + 1)
    for k, lag in enumerate(range(-maxlag, maxlag + 1)):
        if lag >= 0:
            out[k] = (X[:n - lag] * Y[lag:]).mean()
        else:
            out[k] = (X[-lag:] * Y[:n + lag]).mean()
    return out / (sx * sy)


def autocorr(X, maxlag):
    n = X.shape[0]
    v0 = (X * X).mean()
    if v0 <= 0:
        return np.full(maxlag + 1, np.nan)
    return np.array([(X[:n - l] * X[l:]).mean() / v0 for l in range(maxlag + 1)])


def cross_1e(lag, c):
    k = int(np.argmax(c < 1.0 / np.e))
    if k == 0:
        return float("nan")
    x0, x1, y0, y1 = lag[k - 1], lag[k], c[k - 1], c[k]
    return float(x0 + (y0 - 1 / np.e) * (x1 - x0) / (y0 - y1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tau-c", type=float, default=225.56045493417156)
    ap.add_argument("--stride", type=int, default=4)
    ap.add_argument("--tag", nargs="*", default=None)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    summary = {}
    fig, ax = plt.subplots(len(a.cases), 4, figsize=(17, 3.4 * len(a.cases)),
                           squeeze=False, constrained_layout=True)

    for r, c in enumerate(a.cases):
        tag = a.tag[r] if a.tag and r < len(a.tag) else os.path.basename(c.rstrip("/"))
        par = cw.read_params(c)
        st = cw_stream.Stream(c, par)
        tau_m, tau_chi = float(par["tau_m"]), float(par["tau_chi"])
        dt = float(st.steps[1] - st.steps[0]) if st.n > 1 else float(par["nvideo"])
        w = st.n // 3
        maxlag = max(4, int(MAXLAG_FRAC * w))
        lag_steps = np.arange(-maxlag, maxlag + 1) * dt / a.tau_c
        print(f"{tag}: {st.n} frames, dt = {dt:.0f} steps, window {w}, maxlag "
              f"{maxlag} frames = {maxlag*dt/a.tau_c:.1f} tau_c", flush=True)

        summary[tag] = {"tau_m": tau_m, "tau_chi": tau_chi, "dt": dt,
                        "n_frames": st.n, "window": w, "maxlag": maxlag,
                        "lag_tau_c": lag_steps.tolist(), "stages": {}}

        for s, name in enumerate(STAGES):
            i0, i1 = s * w, (s + 1) * w
            D = load_window(st, i0, i1, a.stride)
            rec = {"frames": [int(st.steps[i0]), int(st.steps[i1 - 1])]}
            for (u, v) in PAIRS:
                cc = xcorr(D[u], D[v], maxlag)
                # A field with no spatial variance -- chi in a saturated active phase is
                # EXACTLY 0 everywhere -- makes every correlation involving it undefined.
                # That is the physics, not a failure: report it as such instead of letting
                # nanargmax raise on an all-NaN slice.
                if not np.any(np.isfinite(cc)):
                    rec[f"{u}_{v}"] = {"c": None, "peak": None, "peak_lag_tau_c": None,
                                       "peak_lag_steps": None,
                                       "undefined": f"{u} or {v} has zero spatial variance"}
                    continue
                k = int(np.nanargmax(np.abs(cc)))
                rec[f"{u}_{v}"] = {"c": cc.tolist(), "peak": float(cc[k]),
                                   "peak_lag_tau_c": float(lag_steps[k]),
                                   "peak_lag_steps": float(lag_steps[k] * a.tau_c)}
                ax[r, PAIRS.index((u, v))].plot(lag_steps, cc, lw=1.2, label=name)
            ac = autocorr(D["chi"], maxlag)
            lg = np.arange(maxlag + 1) * dt / a.tau_c
            if not np.any(np.isfinite(ac)):
                rec["chi_auto"] = {"c": None, "tau_1e_tau_c": None, "tau_1e_steps": None,
                                   "undefined": "chi has zero spatial variance"}
                t1e = float("nan")
            else:
                t1e = cross_1e(lg, ac)
                rec["chi_auto"] = {"c": ac.tolist(), "tau_1e_tau_c": t1e,
                                   "tau_1e_steps": t1e * a.tau_c if t1e == t1e else None}
                ax[r, 3].plot(lg, ac, lw=1.2, label=f"{name}  1/e = {t1e:.2f}"
                              if t1e == t1e else name)
            summary[tag]["stages"][name] = rec
            def _fmt(u, v):
                e = rec[f"{u}_{v}"]
                return (f"{u}->{v} undefined" if e["peak"] is None else
                        f"{u}->{v} peak {e['peak']:+.3f} @ {e['peak_lag_tau_c']:+.2f}")
            print(f"    {name:>6}: " + "  ".join(_fmt(u, v) for u, v in PAIRS)
                  + (f"   chi 1/e = {t1e:.2f} tau_c" if t1e == t1e
                     else "   chi 1/e undefined (chi has no spatial variance)"), flush=True)

        for j, ttl in enumerate([f"{u} -> {v}" for u, v in PAIRS] + ["chi autocorr"]):
            ax[r, j].axhline(0, color="0.6", lw=0.8)
            if j < 3:
                ax[r, j].axvline(0, color="0.6", lw=0.8)
                for tv, cl in ((tau_m / a.tau_c, "#b26a00"),
                               (tau_chi / a.tau_c, "#7d3c98")):
                    if abs(tv) < lag_steps[-1]:
                        ax[r, j].axvline(tv, color=cl, ls=":", lw=1.2)
            ax[r, j].grid(alpha=0.25)
            ax[r, j].set_xlabel(r"lag  [$\tau_c$]")
            ax[r, j].legend(fontsize=7)
            ax[r, j].set_title(f"{tag}   {ttl}", fontsize=9)
        ax[r, 0].set_ylabel("correlation")

    fig.suptitle(r"lagged correlations in three stages   (dotted: $\tau_m$ orange, "
                 r"$\tau_\chi$ purple)", fontsize=11)
    fig.savefig(os.path.join(a.out, "corr.png"), dpi=170)
    plt.close(fig)
    with open(os.path.join(a.out, "corr.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(f"\nwrote {a.out}/corr.png and corr.json")


if __name__ == "__main__":
    main()
