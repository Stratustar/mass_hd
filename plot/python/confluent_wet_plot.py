#!/usr/bin/env python
"""Videos and figures for a confluent-wet archive.

CLI-compatible with plot/python/plot_hd.py (`inputdir outdir [...]`) so the cluster job
script can call it in place of the go-or-grow plotter, which reads `phi` and `chi` and
cannot run on confluent-wet frames at all.

Produces, per archive:
  <field>_<first>-<last>_step<n>.mp4   velocity |u|, chi, m, pressure
  <field>_final.png                    the last frame of each, at full dpi
  fields_<t>.png                       the four-panel view at a few times
  scalars_vs_t.png                     u_rms, <chi>, <m>, sigma_P against time

TWO CHOICES WORTH KNOWING ABOUT.

Colour limits are FIXED across each video, computed once from a sample of frames. Per-frame
autoscaling produces a video in which nothing can be compared to anything -- a field that
grows by two decades looks constant, and a field that collapses looks unchanged.

The videos are plain imshow, with no streamlines. Streamlines are integrated curves whose
seeding shifts slightly from frame to frame, so they flicker distractingly in an animation
and trible the render time on an 800x800 lattice; motion is what a video conveys anyway.
The still frames keep them, where they earn their place.
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archive.archive import loadarchive
import plot.confluent as pc

MP4_FPS = 10
PNG_DPI = 200
SAMPLE_FRAMES = 12          # frames used to fix the colour limits


def grid(frame, name):
    lx, ly = frame.parameters["LX"], frame.parameters["LY"]
    return np.array(getattr(frame, name)).reshape((lx, ly))


# ---------------------------------------------------------------- field accessors

def f_speed(frame):
    ux, uy = pc.velocity_field(frame)
    return np.hypot(ux, uy)


def f_chi(frame):
    return grid(frame, "chi")


def f_m(frame):
    return grid(frame, "m")


def f_pressure(frame):
    # <P> is pinned by the initial density in a periodic box and carries no information;
    # only P - <P> does, so every frame is shown relative to its own spatial median.
    P = grid(frame, "pressure")
    return P - np.median(P)


FIELDS = {
    "velocity": dict(fn=f_speed,    cmap="magma",    label=r"$|u|$",                 mode="positive"),
    "chi":      dict(fn=f_chi,      cmap="viridis",  label=r"$\chi$",                mode="unit"),
    "m":        dict(fn=f_m,        cmap="cividis",  label=r"$m$",                   mode="unit"),
    "pressure": dict(fn=f_pressure, cmap="coolwarm", label=r"$P-\mathrm{med}(P)$",   mode="symmetric"),
}


def colour_limits(ar, spec, indices):
    """Fixed (vmin, vmax) for a whole video, from a sample of frames."""
    if spec["mode"] == "unit":
        return 0.0, 1.0
    vals = []
    for i in indices:
        vals.append(spec["fn"](ar[i]).ravel())
    v = np.concatenate(vals)
    if spec["mode"] == "symmetric":
        hi = float(np.percentile(np.abs(v), 99.5)) or 1.0
        return -hi, hi
    return 0.0, float(np.percentile(v, 99.5)) or 1.0


# ---------------------------------------------------------------- rendering

def render(ar, index, name, spec, vmin, vmax, figsize=6.0, dpi=110, streamlines=False):
    frame = ar[index]
    fig = plt.figure(figsize=(figsize*1.18, figsize), dpi=dpi)
    ax = fig.add_axes((0.02, 0.02, 0.80, 0.92))
    if streamlines and name == "velocity":
        cax = pc.velocity(frame, engine=ax, colorbar=False, vmin=vmin, vmax=vmax)
    else:
        cax = ax.imshow(spec["fn"](frame).T, origin="lower", interpolation="nearest",
                        cmap=spec["cmap"], vmin=vmin, vmax=vmax)
    ax.set_xticks([]); ax.set_yticks([])
    p = frame.parameters
    ax.set_title(f"{name}   t = {p['nstart'] + index*p['ninfo']}", fontsize=11)
    cb = fig.colorbar(cax, cax=fig.add_axes((0.845, 0.02, 0.03, 0.92)))
    cb.set_label(spec["label"])
    return fig


def fig_to_rgb(fig):
    fig.canvas.draw()
    arr = np.asarray(fig.canvas.buffer_rgba())[..., :3]
    # h264 wants even dimensions
    h, w = arr.shape[:2]
    return np.ascontiguousarray(arr[:h - h % 2, :w - w % 2])


def write_video(ar, name, spec, outdir, indices, vmin, vmax, fps):
    import imageio_ffmpeg
    out = os.path.join(outdir, f"{name}_{indices[0]}-{indices[-1]}_step"
                               f"{indices[1]-indices[0] if len(indices) > 1 else 1}.mp4")
    writer = None
    try:
        for i in indices:
            fig = render(ar, i, name, spec, vmin, vmax)
            arr = fig_to_rgb(fig)
            plt.close(fig)
            if writer is None:
                h, w = arr.shape[:2]
                writer = imageio_ffmpeg.write_frames(
                    out, (w, h), fps=fps, codec="libx264",
                    quality=8, macro_block_size=2, pix_fmt_out="yuv420p")
                writer.send(None)
            writer.send(arr.tobytes())
            if i % (10*max(1, indices[1]-indices[0] if len(indices) > 1 else 1)) == 0:
                print(f"  {name}: frame {i}/{indices[-1]}", flush=True)
    finally:
        if writer is not None:
            writer.close()
    return out


def four_panel(ar, index, outdir, limits):
    frame = ar[index]
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.4))
    for ax, name in zip(axes, ("velocity", "chi", "m", "pressure")):
        spec = FIELDS[name]
        vmin, vmax = limits[name]
        if name == "velocity":
            cax = pc.velocity(frame, engine=ax, colorbar=False, vmin=vmin, vmax=vmax)
            ax.set_title(rf"velocity   $u_{{\rm rms}}$ = {pc.velocity_rms(frame):.2e}", fontsize=9)
        else:
            f = spec["fn"](frame)
            cax = ax.imshow(f.T, origin="lower", interpolation="nearest",
                            cmap=spec["cmap"], vmin=vmin, vmax=vmax)
            ax.set_title(f"{name}   [{f.min():.4f}, {f.max():.4f}]", fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
        cb = plt.colorbar(cax, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label(spec["label"])
    p = frame.parameters
    xi_N = np.sqrt(p["LL"]/(2*p["CC"]))
    la = np.sqrt(p["LL"]/p["zeta"]) if p["zeta"] > 0 else float("inf")
    fig.suptitle(f"t = {p['nstart'] + index*p['ninfo']}   "
                 f"$\\zeta$={p['zeta']:g}, CC={p['CC']:g}, LL={p['LL']:g}, $\\xi$={p['xi']:g}, "
                 f"$\\tau$={p['tau']:g}  ($\\xi_N$={xi_N:.2f}, $\\ell_a$={la:.1f})   "
                 f"$\\tau_\\chi$={p['tau_chi']:g}, $\\tau_m$={p['tau_m']:g}, Dbio={p['Dbio']:g}",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out = os.path.join(outdir, f"fields_{p['nstart'] + index*p['ninfo']}.png")
    fig.savefig(out, dpi=PNG_DPI); plt.close(fig)
    return out


def scalars_vs_t(ar, outdir):
    n = ar._nframes
    p = ar.parameters
    t, ur, ch, mm, sp = [], [], [], [], []
    for i in range(n):
        fr = ar[i]
        t.append(p["nstart"] + i*p["ninfo"])
        ur.append(pc.velocity_rms(fr))
        ch.append(float(grid(fr, "chi").mean()))
        mm.append(float(grid(fr, "m").mean()))
        sp.append(float(grid(fr, "pressure").std()))
    fig, axes = plt.subplots(1, 3, figsize=(14, 3.8))
    axes[0].semilogy(t, ur); axes[0].set_ylabel(r"$u_{\rm rms}$")
    axes[1].plot(t, ch, label=r"$\langle\chi\rangle$")
    axes[1].plot(t, mm, label=r"$\langle m\rangle$"); axes[1].legend(); axes[1].set_ylim(0, 1)
    axes[2].semilogy(t, sp); axes[2].set_ylabel(r"$\sigma_P$")
    for a in axes: a.set_xlabel("t"); a.grid(alpha=.3)
    fig.tight_layout()
    out = os.path.join(outdir, "scalars_vs_t.png")
    fig.savefig(out, dpi=PNG_DPI); plt.close(fig)
    print(f"  final: u_rms={ur[-1]:.4e}  <chi>={ch[-1]:.4f}  <m>={mm[-1]:.4f}  sigma_P={sp[-1]:.4e}",
          flush=True)
    return out


def main():
    ap = argparse.ArgumentParser(description="Videos and figures for a confluent-wet archive.")
    ap.add_argument("inputdir")
    ap.add_argument("outdir")
    ap.add_argument("--fps", type=int, default=MP4_FPS)
    ap.add_argument("--frame-step", type=int, default=1)
    ap.add_argument("--no-video", action="store_true")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    ar = loadarchive(args.inputdir)
    n = ar._nframes
    print(f"archive {args.inputdir}: {n} frames", flush=True)

    indices = list(range(0, n, max(1, args.frame_step)))
    if indices[-1] != n-1:
        indices.append(n-1)
    sample = sorted(set(np.linspace(0, n-1, min(SAMPLE_FRAMES, n)).astype(int).tolist()))

    limits = {}
    for name, spec in FIELDS.items():
        limits[name] = colour_limits(ar, spec, sample)
        print(f"  {name}: clim = {limits[name]}", flush=True)

    for name, spec in FIELDS.items():
        if not args.no_video:
            print(f"video: {name}", flush=True)
            print("  wrote", write_video(ar, name, spec, args.outdir, indices,
                                         *limits[name], args.fps), flush=True)
        fig = render(ar, n-1, name, spec, *limits[name], figsize=8.0, dpi=PNG_DPI,
                     streamlines=True)
        out = os.path.join(args.outdir, f"{name}_final.png")
        fig.savefig(out, dpi=PNG_DPI); plt.close(fig)
        print("  wrote", out, flush=True)

    for i in sorted(set([0, n//4, n//2, 3*n//4, n-1])):
        print("  wrote", four_panel(ar, i, args.outdir, limits), flush=True)
    print("  wrote", scalars_vs_t(ar, args.outdir), flush=True)


if __name__ == "__main__":
    main()
