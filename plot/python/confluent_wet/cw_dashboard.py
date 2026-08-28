"""A multi-panel dashboard video for one confluent-wet closed-loop run.

Four fields animated together -- chi, m, P, |u| -- over live time series of the loop's own
state, with a cursor marking the frame being shown. The point is that the loop is a CHAIN:
P drives m drives chi drives the activity drives the flow drives P. A single-field video
shows one link; watching all four at once is what makes the lag visible, because the same
structure appears in P, then in m, then in chi.

TWO THINGS IT DOES NOT DO, deliberately:
  * it does not re-scan the archive for the time series -- those come from the per-case
    cw_scan.json, which already holds them.  A frame is 42 MB of JSON at L = 400 and the
    scan is the expensive part, so a second pass would double the cost of the whole thing
    for data already on disk.  Without that file it falls back to a cheap pre-pass.
  * it does not autoscale per frame.  Colour limits are fixed across the video from a
    sample, because a field that grows by two decades under per-frame autoscaling animates
    as a constant (learned 20260827).

Usage:
  cw_dashboard.py <inputdir> <outdir> [--scan cw_scan.json] [--px 380] [--fps 12]
                  [--bitrate 3000k]
"""
import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cw_common as cw

PANELS = [("chi",      "coolwarm", r"$\chi$   phenotype",  (0, 1)),
          ("m",        "viridis",  r"$m$   memory",        (0, 1)),
          ("P",        "RdBu_r",   r"$P$   pressure",      None),
          ("speed",    "inferno",  r"$|u|$",               None)]


def field(fr, name):
    return np.sqrt(fr["ux"]**2 + fr["uy"]**2) if name == "speed" else fr[name]


def limits(oa, nf, nsample=8):
    """Fixed colour limits from a spread of frames; symmetric for the signed pressure."""
    out = {}
    idx = np.unique(np.linspace(0, nf - 1, nsample).astype(int))
    acc = {n: [] for n, _, _, cl in PANELS if cl is None}
    for i in idx:
        fr = cw.load_frame(oa, int(i))
        for n in acc:
            acc[n].append(np.abs(field(fr, n)).ravel())
    for n, _, _, cl in PANELS:
        if cl is not None:
            out[n] = cl
        else:
            v = float(np.percentile(np.concatenate(acc[n]), 99.5)) or 1.0
            out[n] = (-v, v) if n == "P" else (0.0, v)
    return out


def series_for(root, scan_path, oa, nf):
    """The loop's time series, from the scan JSON if it exists, else a cheap pre-pass."""
    case = os.path.basename(root.rstrip("/"))
    if scan_path and os.path.exists(scan_path):
        with open(scan_path) as f:
            rows = json.load(f)
        rows = rows if isinstance(rows, list) else [rows]
        for r in rows:
            if r.get("case") == case and "series" in r:
                S = dict(r["series"])
                S["_meta"] = r
                return S
    print("  no scan JSON for this case -- pre-scanning for the time series", flush=True)
    S = {k: [] for k in ("t", "chibar", "mbar", "urms", "nd", "Pmed")}
    for i in range(nf):
        fr = cw.load_frame(oa, i)
        S["t"].append(cw.ph.frame_time(oa, i))
        S["chibar"].append(float(fr["chi"].mean()))
        S["mbar"].append(float(fr["m"].mean()))
        S["urms"].append(float(np.sqrt(np.mean(fr["ux"]**2 + fr["uy"]**2))))
        S["nd"].append(float(cw.n_defects(fr["qxx"], fr["qyx"])))
        S["Pmed"].append(float(np.median(fr["P"])))
    S["_meta"] = {}
    return S


def build(root, outdir, scan_path, px, fps, bitrate):
    import imageio_ffmpeg
    nf = cw.frame_count(root)
    oa = cw.loadarchive(root)
    par = oa.parameters
    L = int(par["LX"])
    zeta = float(par["zeta"])
    case = os.path.basename(root.rstrip("/"))
    print(f"dashboard {case}: {nf} frames, L = {L}", flush=True)

    S = series_for(root, scan_path, oa, nf)
    meta = S.pop("_meta", {})
    T = np.asarray(S["t"], float)
    clim = limits(oa, nf)

    dpi = 100.0
    fig = plt.figure(figsize=((4 * px * 1.14) / dpi, (px + 300) / dpi), dpi=dpi)
    gs = fig.add_gridspec(2, 4, height_ratios=[px, 250], hspace=0.28, wspace=0.22,
                          left=0.045, right=0.99, top=0.90, bottom=0.09)

    ims = {}
    for j, (n, cmap, lab, _) in enumerate(PANELS):
        ax = fig.add_subplot(gs[0, j])
        ims[n] = ax.imshow(np.zeros((L, L)), origin="lower", interpolation="nearest",
                           cmap=cmap, vmin=clim[n][0], vmax=clim[n][1])
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(lab, fontsize=12)
        cb = fig.colorbar(ims[n], ax=ax, fraction=0.046, pad=0.02)
        cb.ax.tick_params(labelsize=7)

    # --- time series, drawn once; only the cursor moves
    traces = [
        (r"loop state", [("chibar", r"$\langle\chi\rangle$", "C0"),
                         ("mbar", r"$\langle m\rangle$", "C1")], 0.5),
        (r"$\zeta_{\rm eff}=\zeta(1-\langle\chi\rangle)$",
         [("_zeff", None, "C3")], None),
        (r"$u_{\rm rms}$", [("urms", None, "C2")], None),
        (r"$N_{\rm defect}$", [("nd", None, "C4")], None),
    ]
    S["_zeff"] = [zeta * (1 - c) for c in S["chibar"]]
    cursors = []
    for j, (ylab, keys, hline) in enumerate(traces):
        ax = fig.add_subplot(gs[1, j])
        for k, lab, c in keys:
            ax.plot(T, S[k], c=c, lw=1.1, label=lab)
        if hline is not None:
            ax.axhline(hline, ls=":", c="k", lw=.8)
        if any(lab for _, lab, _ in keys):
            ax.legend(fontsize=9, loc="best")
        ax.set_ylabel(ylab, fontsize=10)
        ax.set_xlabel("t", fontsize=9)
        ax.tick_params(labelsize=8)
        ax.grid(alpha=.3)
        ax.set_xlim(T[0], T[-1])
        cursors.append(ax.axvline(T[0], c="k", lw=1.3, alpha=.75))

    tm = float(par.get("tau_m", np.nan))
    tc = float(par.get("tau_chi", np.nan))
    sub = (f"{case}    A = {meta.get('A', float('nan')):.0f}    "
           f"switch-sign {int(par.get('switch_sign', 1)):+d}    Dbio = {par.get('Dbio', 0):g}    "
           f"$\\tau_m$ = {tm:.0f}    $\\tau_\\chi$ = {tc:.0f}")
    if meta.get("t_eddy"):
        sub += (f"    $\\tau_m/t_{{\\rm eddy}}$ = {tm/meta['t_eddy']:.2f}"
                f"    $L_\\chi/d$ = {meta['L_chi']/meta['d']:.2f}")
    title = fig.suptitle("", fontsize=13)

    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, "dashboard.mp4")
    writer = None
    try:
        for i in range(nf):
            try:
                fr = cw.load_frame(oa, i)
            except Exception:
                print(f"  stopping at unreadable frame {i}", flush=True)
                break
            for n, _, _, _ in PANELS:
                ims[n].set_data(field(fr, n).T)
            t = cw.ph.frame_time(oa, i)
            for cur in cursors:
                cur.set_xdata([t, t])
            title.set_text(sub + f"\nt = {t:.0f}")
            fig.canvas.draw()
            arr = np.asarray(fig.canvas.buffer_rgba())[..., :3]
            h, w = arr.shape[:2]
            arr = np.ascontiguousarray(arr[:h - h % 2, :w - w % 2])
            if writer is None:
                h, w = arr.shape[:2]
                enc = dict(bitrate=bitrate) if bitrate else dict(quality=7)
                writer = imageio_ffmpeg.write_frames(
                    out, (w, h), fps=fps, codec="libx264", macro_block_size=2,
                    pix_fmt_out="yuv420p", **enc)
                writer.send(None)
            writer.send(arr.tobytes())
            if i % 40 == 0:
                print(f"  frame {i}/{nf}", flush=True)
    finally:
        if writer is not None:
            writer.close()
        plt.close(fig)
    mb = os.path.getsize(out) / 1e6 if os.path.exists(out) else 0.0
    print(f"  wrote {out}  ({mb:.1f} MB, {nf} frames @ {fps} fps)", flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputdir")
    ap.add_argument("outdir")
    ap.add_argument("--scan", default=None,
                    help="cw_scan.json holding this case's series (avoids a second pass)")
    ap.add_argument("--px", type=int, default=380)
    ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--bitrate", default="3000k")
    a = ap.parse_args()
    build(a.inputdir, a.outdir, a.scan, a.px, a.fps, a.bitrate)


if __name__ == "__main__":
    main()
