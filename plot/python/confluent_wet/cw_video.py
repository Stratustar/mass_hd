"""In-job plotting for a confluent-wet closed-loop run: a chi video, stills, and time series.

Same `inputdir outdir [...]` contract as plot_hd.py, so submit_case.sh can select it with
PLOT_SCRIPT.

WHY THIS EXISTS RATHER THAN confluent_wet_plot.py.  That driver builds a NEW matplotlib
figure, axes and colorbar for every frame; at L = 800 x 140 frames it spent 1 h 32 min against
31 min of simulation -- 75% of the wall time was plotting, and at 72 production runs that is
the dominant cost of the whole campaign.  Here one figure, one image artist and one colorbar
are built once and only the pixel data and the title change per frame.

RESOLUTION AND FILE SIZE.  The image is rendered at an explicit pixel size (default 900 px of
field, so a 400-cell lattice is upsampled 2.25x and defect cores are actually visible) with
`interpolation="nearest"` so the upsampling never invents structure.  Size is controlled by
the x264 quality and the frame count, both explicit; a 181-frame 900 px chi video lands in the
few-MB range.

The whole archive is the steady state by construction: the production runcards set nstart past
the transient, so nothing here needs to trim a warm-up.
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cw_common as cw

# `key` is the name the field carries in a cw_common.load_frame dict, which is the ANALYSIS
# name and not the name in the archive -- the pressure is "P" there, not "pressure".
FIELDS = {
    "chi":      dict(key="chi", cmap="coolwarm", label=r"$\chi$  (phenotype)", clim=(0, 1)),
    "m":        dict(key="m",   cmap="viridis",  label=r"$m$  (memory)",       clim=(0, 1)),
    "speed":    dict(key=None,  cmap="inferno",  label=r"$|u|$",               clim=None),
    "pressure": dict(key="P",   cmap="RdBu_r",   label=r"$P$",                 clim=None),
    "q2":       dict(key="q2",  cmap="magma",    label=r"$q^2$ (order)",       clim=(0, 1)),
}


def field_of(fr, name):
    if name == "speed":
        return np.sqrt(fr["ux"]**2 + fr["uy"]**2)
    return fr[FIELDS[name]["key"]]


def sample_limits(oa, name, indices, pct=99.5):
    """Colour limits fixed ACROSS the video from a sample of frames.

    Per-frame autoscaling makes a field that grows by two decades look constant; that lesson
    is already paid for (20260827).  Symmetric limits for the signed pressure.
    """
    spec = FIELDS[name]
    if spec["clim"] is not None:
        return spec["clim"]
    vals = []
    for i in indices[:: max(1, len(indices) // 12)]:
        vals.append(np.abs(field_of(cw.load_frame(oa, i), name)).ravel())
    v = float(np.percentile(np.concatenate(vals), pct))
    return (-v, v) if name == "pressure" else (0.0, v)


class Canvas:
    """One figure, one image, one colorbar -- reused for every frame."""

    def __init__(self, L, name, clim, px):
        spec = FIELDS[name]
        dpi = 100.0
        fig_w = (px * 1.22) / dpi
        fig_h = px / dpi
        self.fig = plt.figure(figsize=(fig_w, fig_h), dpi=dpi)
        ax = self.fig.add_axes((0.01, 0.01, 0.80, 0.92))
        self.im = ax.imshow(np.zeros((L, L)), origin="lower", interpolation="nearest",
                            cmap=spec["cmap"], vmin=clim[0], vmax=clim[1])
        ax.set_xticks([]); ax.set_yticks([])
        self.title = ax.set_title("", fontsize=11)
        cb = self.fig.colorbar(self.im, cax=self.fig.add_axes((0.835, 0.01, 0.028, 0.92)))
        cb.set_label(spec["label"])

    def draw(self, data, title):
        self.im.set_data(data.T)
        self.title.set_text(title)
        self.fig.canvas.draw()
        arr = np.asarray(self.fig.canvas.buffer_rgba())[..., :3]
        h, w = arr.shape[:2]                      # h264 wants even dimensions
        return np.ascontiguousarray(arr[:h - h % 2, :w - w % 2])

    def close(self):
        plt.close(self.fig)


def still_panel(fr, t, outdir, px):
    """A four-panel chi / m / P / |u| snapshot."""
    names = ["chi", "m", "pressure", "speed"]
    fig, axes = plt.subplots(2, 2, figsize=(px / 90.0, px / 90.0), dpi=110)
    for ax, nm in zip(axes.ravel(), names):
        d = field_of(fr, nm)
        spec = FIELDS[nm]
        cl = spec["clim"]
        if cl is None:
            v = float(np.percentile(np.abs(d), 99.5)) or 1.0
            cl = (-v, v) if nm == "pressure" else (0.0, v)
        im = ax.imshow(d.T, origin="lower", interpolation="nearest",
                       cmap=spec["cmap"], vmin=cl[0], vmax=cl[1])
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(spec["label"], fontsize=10)
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.suptitle(f"t = {t:.0f}", fontsize=12)
    fig.tight_layout()
    out = os.path.join(outdir, f"fields_{int(t)}.png")
    fig.savefig(out, dpi=110); plt.close(fig)
    print(f"  wrote {out}", flush=True)


def scalars_figure(T, S, zeta, root, outdir):
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), dpi=110)
    a = axes.ravel()
    a[0].plot(T, S["chibar"], label=r"$\langle\chi\rangle$")
    a[0].plot(T, S["mbar"], label=r"$\langle m\rangle$")
    a[0].axhline(0.5, ls=":", c="k", lw=.8); a[0].legend(); a[0].set_ylabel("loop state")
    a[1].plot(T, [zeta * (1 - c) for c in S["chibar"]], c="C3")
    a[1].set_ylabel(r"$\zeta_{\rm eff}=\zeta(1-\langle\chi\rangle)$")
    a[2].plot(T, S["urms"], c="C2"); a[2].set_ylabel(r"$u_{\rm rms}$")
    a[3].plot(T, S["nd"], c="C4"); a[3].set_ylabel(r"$N_{\rm defect}$")
    for ax in a:
        ax.set_xlabel("t"); ax.grid(alpha=.3)
    fig.suptitle(os.path.basename(root.rstrip("/")), fontsize=12)
    fig.tight_layout()
    out = os.path.join(outdir, "scalars_vs_t.png")
    fig.savefig(out, dpi=110); plt.close(fig)
    print(f"  wrote {out}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputdir")
    ap.add_argument("outdir")
    ap.add_argument("--videos", nargs="*", default=["chi"],
                    help="fields to animate; chi is the deliverable, each extra one costs "
                         "only rendering now, not another pass over the archive")
    ap.add_argument("--px", type=int, default=900)
    ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--quality", type=int, default=7)
    ap.add_argument("--bitrate", default="2000k",
                    help="explicit x264 bitrate; empty string falls back to --quality")
    a = ap.parse_args()

    import imageio_ffmpeg
    os.makedirs(a.outdir, exist_ok=True)
    root = a.inputdir
    nf = cw.frame_count(root)
    oa = cw.loadarchive(root)
    zeta = float(oa.parameters["zeta"])
    L = int(oa.parameters["LX"])
    print(f"archive {root}: {nf} frames, L = {L}", flush=True)

    # ONE pre-pass, and only for the auto-scaled fields: chi and m have fixed limits, so a
    # chi-only run (the production default) touches the archive exactly once.
    clims = {n: (FIELDS[n]["clim"] if FIELDS[n]["clim"] is not None
                 else sample_limits(oa, n, list(range(nf)))) for n in a.videos}
    canvases = {n: Canvas(L, n, clims[n], a.px) for n in a.videos}
    writers = {n: None for n in a.videos}
    still_at = sorted({0, nf // 2, nf - 1})

    T, S = [], {k: [] for k in ("chibar", "mbar", "urms", "nd")}
    try:
        for i in range(nf):
            try:
                fr = cw.load_frame(oa, i)          # THE one read of this frame
            except Exception:
                print(f"  stopping at unreadable frame {i}", flush=True)
                break
            t = cw.ph.frame_time(oa, i)

            for n in a.videos:
                arr = canvases[n].draw(field_of(fr, n), f"{n}    t = {t:.0f}")
                if writers[n] is None:
                    h, w = arr.shape[:2]
                    # bitrate, not quality, is what bounds the file: at quality=7 a 201-frame
                    # 900px turbulent |u| came out at 18.6 MB, against 4.3 MB at 2000k.
                    enc = dict(bitrate=a.bitrate) if a.bitrate else dict(quality=a.quality)
                    writers[n] = imageio_ffmpeg.write_frames(
                        os.path.join(a.outdir, f"{n}.mp4"), (w, h), fps=a.fps,
                        codec="libx264", macro_block_size=2, pix_fmt_out="yuv420p", **enc)
                    writers[n].send(None)
                writers[n].send(arr.tobytes())

            T.append(t)
            S["chibar"].append(float(fr["chi"].mean()))
            S["mbar"].append(float(fr["m"].mean()))
            S["urms"].append(float(np.sqrt(np.mean(fr["ux"]**2 + fr["uy"]**2))))
            S["nd"].append(float(cw.n_defects(fr["qxx"], fr["qyx"])))

            if i in still_at:
                still_panel(fr, t, a.outdir, a.px)
            if i % 40 == 0:
                print(f"  frame {i}/{nf}", flush=True)
    finally:
        for n, wr in writers.items():
            if wr is not None:
                wr.close()
                mb = os.path.getsize(os.path.join(a.outdir, f"{n}.mp4")) / 1e6
                print(f"  wrote {a.outdir}/{n}.mp4  ({mb:.1f} MB, {len(T)} frames @ "
                      f"{a.fps} fps, clim={clims[n][0]:.3g}..{clims[n][1]:.3g})", flush=True)
        for c in canvases.values():
            c.close()

    if T:
        scalars_figure(T, S, zeta, root, a.outdir)


if __name__ == "__main__":
    main()
