#!/usr/bin/env python3
"""Stage B of the s = 1/3 campaign: the per-run series, the three-panel video, the two figures.

Three stages, because they have different costs and different lifetimes. The positional
contract is the project's usual `inputdir outdir`, so the per-run stages drop straight into
scripts_cluster/submit_analysis_array.sh:

    cw_s3_b.py <sim_dir> <out_dir> --calib C                one run  -> series.txt + chi_P_m.mp4
    cw_s3_b.py <sim_dir> <out_dir> --calib C --stage series one run  -> series.txt only
    cw_s3_b.py <results_root> <out> --calib C --stage figs  all runs -> the two figures

WHY THE SERIES IS NOT COMPUTED FROM THE VIDEO PIXELS. video_meta.csv carries <chi>, std(chi),
<m>, u_rms and std(P) written by the C++ in double precision from the FULL-RESOLUTION field,
on the same clock as the video frames. The uint8 stream is a picture of the field; the CSV is
the measurement of it. Var_x(chi) is std(chi)^2 -- the SPATIAL variance over the box at one
instant, which is what separates a mixed state (large) from a saturated one (~0).

THE TWO FIGURES read only the series files, so they are cheap to re-make and they never
depend on a frame being present. Each point is the mean over the closing WINDOW_TAU_C of the
run and the error bar is the standard deviation IN THAT WINDOW -- a spread in time, not a
standard error: it says how much the state still moves, which is the thing that distinguishes
a settled phase from one that is drifting or oscillating.

READING THEM. Where the chi0 and chi1 curves coincide the outcome does not depend on the
initial condition -- either a mixed state or a single high phase. Where they separate is
tau_freeze. The mixed-to-single-high transition is CONTINUOUS along either curve (chibar_hi
falls slowly with tau_m) and is not a crossing; tau_x, where the low phase first holds a
domain, is not in these two curves at all and has to be read off the leftright and patches
videos.
"""
import argparse
import glob
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cw_common as cw                                        # noqa: E402
import cw_stream                                              # noqa: E402

WINDOW_TAU_C = 100.0          # the averaging window at the end of each run
P_DISPLAY_SIGMA = 2.0         # the video's P window, in campaign sigma_P
STARTS = ("chi0", "chi1", "leftright", "patches")
FIG_STARTS = ("chi0", "chi1")  # the two the figures carry


# --------------------------------------------------------------------- series

def series_of(sim_dir, tau_c=None):
    """(t/tau_c, <chi>, Var_x(chi)) plus the run's identity, from video_meta.csv."""
    par = cw.read_params(sim_dir)
    st = cw_stream.Stream(sim_dir, par)
    tc = float(tau_c) if tau_c else float("nan")
    t = st.steps / tc
    chi = np.asarray(st.meta["chi_mean"][:st.n], float)
    var = np.asarray(st.meta["chi_std"][:st.n], float) ** 2
    return {
        "case": os.path.basename(sim_dir.rstrip("/")),
        "tau_c": tc,
        "tau_m": float(par["tau_m"]), "tau_chi": float(par["tau_chi"]),
        "tau_m_over_tau_c": float(par["tau_m"]) / tc,
        "chi_config": str(par.get("chi_config", "")),
        "chi0": float(par.get("chi0", np.nan)),
        "seed": int(par.get("seed", -1)), "chi_seed": int(par.get("chi_seed", 0)),
        "nsteps": int(par["nsteps"]), "nvideo": int(par["nvideo"]),
        "frames": int(st.n),
        "t": t, "chi_mean": chi, "var_x_chi": var,
        "m_mean": np.asarray(st.meta["m_mean"][:st.n], float),
    }


def start_label(s):
    """The runcard's (chi-config, chi0) back to the campaign's name for the start."""
    cc, c0 = s["chi_config"], s["chi0"]
    if cc == "stripe":
        return "leftright"
    if cc == "binary-noise":
        return "patches"
    return "chi1" if c0 > 0.5 else "chi0"


def write_series(s, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "series.txt")
    with open(path, "w") as fh:
        fh.write(f"# {s['case']}  start={start_label(s)}  tau_m/tau_c={s['tau_m_over_tau_c']:.3f}"
                 f"  tau_c={s['tau_c']:.1f} steps  seed={s['seed']}/{s['chi_seed']}\n")
        fh.write("# t/tau_c   <chi>        Var_x(chi)\n")
        for t, c, v in zip(s["t"], s["chi_mean"], s["var_x_chi"]):
            fh.write(f"{t:10.4f} {c:12.6f} {v:14.8f}\n")
    meta = {k: v for k, v in s.items() if not isinstance(v, np.ndarray)}
    meta["start"] = start_label(s)
    w = s["t"] >= s["t"][-1] - WINDOW_TAU_C
    meta["window_tau_c"] = WINDOW_TAU_C
    meta["window_frames"] = int(w.sum())
    meta["chi_mean_window"] = float(s["chi_mean"][w].mean())
    meta["chi_std_window"] = float(s["chi_mean"][w].std())
    meta["var_mean_window"] = float(s["var_x_chi"][w].mean())
    meta["var_std_window"] = float(s["var_x_chi"][w].std())
    # Drift across the window: the honest flag for "this has not settled". Compared with the
    # window's own scatter, because a slow drift and a large fluctuation look the same in a
    # mean +/- std and are not the same statement.
    half = int(w.sum()) // 2
    cw_ = s["chi_mean"][w]
    meta["chi_drift_window"] = float(cw_[half:].mean() - cw_[:half].mean())
    with open(os.path.join(out_dir, "series.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    return path, meta


# ---------------------------------------------------------------------- video

def render_three(sim_dir, out_path, tau_c, sigma_P, every=4, fps=25, bitrate=None):
    """chi | P | m side by side, one caption, fixed scales.

    `every` subsamples the 0.05 tau_c series to the 0.2 tau_c the campaign asks a video to
    animate: the stream is stored fine enough to MEASURE with and is four times finer than
    anything worth watching.
    """
    import imageio_ffmpeg
    from matplotlib import colormaps
    from PIL import Image, ImageDraw, ImageFont

    par = cw.read_params(sim_dir)
    st = cw_stream.Stream(sim_dir, par)
    panels = [("chi", "viridis", 0.0, 1.0),
              ("P", "RdBu_r", -P_DISPLAY_SIGMA * sigma_P, P_DISPLAY_SIGMA * sigma_P),
              ("m", "magma", 0.0, 1.0)]
    luts = {n: (colormaps[c](np.linspace(0, 1, 256))[:, :3] * 255).astype(np.uint8)
            for n, c, _, _ in panels}
    raws = {n: st.raw(n) for n, _, _, _ in panels}
    gap = 4
    w = 3 * st.nx + 2 * gap
    tm_over_tc = float(par["tau_m"]) / tau_c
    start = start_label({"chi_config": str(par.get("chi_config", "")),
                         "chi0": float(par.get("chi0", 0.0))})
    head = f"tau_m/tau_c = {tm_over_tc:.3f}   {start}"
    # TWO caption rows: the run's identity and clock on top, the three colour windows under
    # the panels they belong to. One row would put the first panel's label on top of the
    # main caption at any realistic width. Fonts are fitted, never truncated -- PIL cuts the
    # right-hand end, which is where <chi> is.
    widest = f"{head}   t/tau_c = {st.steps[-1]/tau_c:7.2f}   <chi> = 0.000"
    font = cw_stream._fit_font(widest, w, start=15)
    lab_font = cw_stream._fit_font("P [-0.0673, 0.0673]", st.nx, start=13)
    bar_h = 44
    bar_h += (st.ny + bar_h) % 2

    idx = range(0, st.n, max(1, int(every)))
    writer = None
    try:
        for i in idx:
            cols = []
            for name, _, dlo, dhi in panels:
                lo, hi = st.limits(name)
                phys = lo + raws[name][i].astype(np.float32) * ((hi - lo) / 255.0)
                x = (phys - dlo) / (dhi - dlo)
                q = np.clip(np.rint(x * 255.0), 0, 255).astype(np.uint8)
                cols.append(luts[name][np.flipud(q.T)])
            sep = np.zeros((st.ny, gap, 3), dtype=np.uint8)
            img = np.concatenate([cols[0], sep, cols[1], sep, cols[2]], axis=1)
            bar = Image.new("RGB", (w, bar_h), (16, 16, 16))
            d = ImageDraw.Draw(bar)
            d.text((6, bar_h // 4),
                   f"{head}   t/tau_c = {st.steps[i]/tau_c:7.2f}   <chi> = "
                   f"{st.meta['chi_mean'][i]:.3f}", fill=(235, 235, 235), font=font,
                   anchor="lm")
            for j, (name, _, dlo, dhi) in enumerate(panels):
                d.text((j * (st.nx + gap) + st.nx // 2, 3 * bar_h // 4),
                       f"{name} [{dlo:.3g}, {dhi:.3g}]", fill=(150, 150, 150),
                       font=lab_font, anchor="mm")
            frame = np.concatenate([img, np.asarray(bar, dtype=np.uint8)], axis=0)
            if writer is None:
                enc = dict(bitrate=bitrate) if bitrate else dict(quality=7)
                writer = imageio_ffmpeg.write_frames(
                    out_path, (frame.shape[1], frame.shape[0]), fps=fps,
                    codec="libx264", macro_block_size=2, pix_fmt_out="yuv420p", **enc)
                writer.send(None)
            writer.send(frame.tobytes())
    finally:
        if writer is not None:
            writer.close()
    return os.path.getsize(out_path)


# --------------------------------------------------------------------- figures

def make_figs(root, out_dir, tau_c):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    runs = []
    for p in sorted(glob.glob(os.path.join(root, "*", "series.json"))):
        with open(p) as fh:
            runs.append(json.load(fh))
    if not runs:
        raise RuntimeError(f"no series.json under {root} -- run the `series` stage first")

    by = {s: sorted([r for r in runs if r["start"] == s],
                    key=lambda r: r["tau_m_over_tau_c"]) for s in STARTS}
    style = {"chi0": ("C3", "o", "chi0  (all active)"),
             "chi1": ("C0", "s", "chi1  (all passive)")}

    out = {}
    for key, ylab, fname, title in (
            ("chi_mean", r"$\langle\chi\rangle$", "fig1_chi_mean.png",
             "domain-averaged phenotype"),
            ("var", r"Var$_x(\chi)$", "fig2_var_chi.png",
             "spatial variance of the phenotype")):
        fig, ax = plt.subplots(figsize=(7.2, 4.6))
        for s in FIG_STARTS:
            rs = by[s]
            if not rs:
                continue
            x = [r["tau_m_over_tau_c"] for r in rs]
            y = [r["chi_mean_window" if key == "chi_mean" else "var_mean_window"] for r in rs]
            e = [r["chi_std_window" if key == "chi_mean" else "var_std_window"] for r in rs]
            c, m, lab = style[s]
            ax.errorbar(x, y, yerr=e, color=c, marker=m, ms=5, lw=1.4, capsize=3,
                        label=lab)
            # a run still moving across its own averaging window is not a settled point
            for r in rs:
                if abs(r.get("chi_drift_window", 0.0)) > 0.05:
                    ax.plot([r["tau_m_over_tau_c"]],
                            [r["chi_mean_window" if key == "chi_mean" else "var_mean_window"]],
                            marker="x", ms=11, color="k", mew=1.6, ls="none")
        ax.set_xlabel(r"$\tau_m / \tau_c$")
        ax.set_ylabel(ylab)
        ax.set_title(f"{title} -- mean over the closing {WINDOW_TAU_C:g} "
                     r"$\tau_c$, bars = std in that window", fontsize=10)
        ax.grid(alpha=0.3)
        ax.legend()
        if key == "chi_mean":
            ax.set_ylim(-0.05, 1.05)
        fig.tight_layout()
        path = os.path.join(out_dir, fname)
        fig.savefig(path, dpi=150)
        plt.close(fig)
        out[key] = path
        print(f"wrote {path}")

    # the numbers behind the figures, so a reading can be checked without re-plotting
    tab = os.path.join(out_dir, "figs.json")
    with open(tab, "w") as fh:
        json.dump({"window_tau_c": WINDOW_TAU_C, "tau_c": tau_c,
                   "runs": {s: [{k: r[k] for k in
                                 ("case", "tau_m_over_tau_c", "chi_mean_window",
                                  "chi_std_window", "var_mean_window", "var_std_window",
                                  "chi_drift_window", "seed")}
                                for r in by[s]] for s in STARTS}}, fh, indent=1)
    print(f"wrote {tab}")

    # the separation, which is what the figures are read for
    a = {r["tau_m_over_tau_c"]: r for r in by["chi0"]}
    b = {r["tau_m_over_tau_c"]: r for r in by["chi1"]}
    print(f"\n{'tau_m/tau_c':>11} {'chi0':>9} {'chi1':>9} {'separation':>11} "
          f"{'Var chi0':>9} {'Var chi1':>9}")
    for g in sorted(set(a) & set(b)):
        d = abs(b[g]["chi_mean_window"] - a[g]["chi_mean_window"])
        print(f"{g:11.3f} {a[g]['chi_mean_window']:9.4f} {b[g]['chi_mean_window']:9.4f} "
              f"{d:11.4f} {a[g]['var_mean_window']:9.5f} {b[g]['var_mean_window']:9.5f}")
    return out


# ------------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path", help="sim dir (series/video/both) or results root (figs)")
    ap.add_argument("out", nargs="?", default=None)
    ap.add_argument("--stage", choices=["series", "video", "both", "figs"], default="both",
                    help="default `both`: the series and the video for one run")
    ap.add_argument("--calib", default=None, help="calib_s3.json; supplies tau_c and sigma_P")
    ap.add_argument("--tau-c", type=float, default=None)
    ap.add_argument("--sigma-p", type=float, default=None)
    ap.add_argument("--every", type=int, default=4, help="video frame stride (0.05 -> 0.2 tau_c)")
    ap.add_argument("--fps", type=int, default=25)
    a = ap.parse_args()

    tau_c, sigma_P = a.tau_c, a.sigma_p
    if a.calib:
        with open(a.calib) as fh:
            c = json.load(fh)
        tau_c = tau_c or float(c["tau_c"])
        sigma_P = sigma_P or float(c["sigma_P"])
    if tau_c is None:
        raise SystemExit("need --calib or --tau-c: every axis of this stage is in tau_c")

    if a.stage in ("series", "both"):
        out = a.out or a.path
        s = series_of(a.path, tau_c)
        path, meta = write_series(s, out)
        print(f"{meta['case']}: start {meta['start']}, tau_m/tau_c {meta['tau_m_over_tau_c']:.3f}, "
              f"{meta['frames']} frames")
        print(f"  window {meta['window_frames']} frames: <chi> = {meta['chi_mean_window']:.4f} "
              f"+/- {meta['chi_std_window']:.4f}, Var_x = {meta['var_mean_window']:.5f} "
              f"+/- {meta['var_std_window']:.5f}, drift {meta['chi_drift_window']:+.4f}")
        print(f"  wrote {path}", flush=True)

    if a.stage in ("video", "both"):
        out = a.out or a.path
        os.makedirs(out, exist_ok=True)
        if sigma_P is None:
            raise SystemExit("the video needs sigma_P (--calib or --sigma-p): its P window "
                             "is the campaign's +/- 2 sigma_P, never the run's own range")
        vp = os.path.join(out, "chi_P_m.mp4")
        nb = render_three(a.path, vp, tau_c, sigma_P, every=a.every, fps=a.fps)
        print(f"wrote {vp} ({nb/1e6:.1f} MB)", flush=True)

    if a.stage == "figs":
        out = a.out or a.path
        os.makedirs(out, exist_ok=True)
        make_figs(a.path, out, tau_c)


if __name__ == "__main__":
    main()
