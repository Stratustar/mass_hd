#!/usr/bin/env python3
"""THE dashboard video standard for confluent-wet.

Four videos per parameter set, always these four and always in this order:

    u        speed |u| in magma with the flow topology drawn on top as white streamlines
    p        pressure, diverging scale centred on zero
    m        mechanical memory
    chi      phenotype

Each panel carries a one-symbol title and a colour bar. Nothing else: no time stamp, no
parameter line, no case name. The dashboard's buttons say which run this is, and a caption
says what the parameters are; burning either into the pixels means the video cannot be
reused anywhere else.

STREAMLINE SEEDING IS FIXED, and that is the whole reason this exists rather than a flag on
the old driver. matplotlib seeds streamplot from its own density grid, and those seeds move
with the flow, so consecutive frames get different curves and the animation strobes -- which
is why the earlier driver dropped streamlines from the videos and kept them only in the
stills (20260827). Passing an explicit, frame-independent `start_points` grid pins every
curve's origin, so a vortex that persists is drawn by the same streamline in every frame and
the eye can follow it.

Colour limits are fixed across each video from a sample of frames. Per-frame autoscaling
animates a field that grows by two decades as a constant.

Usage:
  cw_dash.py <inputdir> <outdir> [--px 460] [--fps 12] [--bitrate 2200k]
             [--seed-step 24] [--density 1.4] [--no-streamlines]
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

# the standard order; changing it changes the standard
PANELS = [
    ("u",        r"$|u|$",  "magma",   None),
    ("pressure", r"$P$",    "RdBu_r",  None),
    ("m",        r"$m$",    "viridis", (0.0, 1.0)),
    ("chi",      r"$\chi$", "coolwarm", (0.0, 1.0)),
]


def panel_data(fr, name):
    if name == "u":
        return np.hypot(fr["ux"], fr["uy"])
    if name == "pressure":
        return fr["P"]
    return fr[name]


def limits(oa, nf, nsample=10):
    """Fixed colour limits from a spread of frames; symmetric for the signed pressure."""
    out, acc = {}, {}
    idx = np.unique(np.linspace(0, nf - 1, nsample).astype(int))
    for n, _, _, cl in PANELS:
        if cl is None:
            acc[n] = []
    for i in idx:
        fr = cw.load_frame(oa, int(i))
        for n in acc:
            acc[n].append(np.abs(panel_data(fr, n)).ravel())
    for n, _, _, cl in PANELS:
        if cl is not None:
            out[n] = cl
        else:
            v = float(np.percentile(np.concatenate(acc[n]), 99.5)) or 1.0
            out[n] = (-v, v) if n == "pressure" else (0.0, v)
    return out


def seed_grid(L, step):
    """A frame-independent seed lattice, inset so no curve starts on the boundary."""
    a = np.arange(step / 2.0, L - step / 2.0 + 1e-9, step)
    X, Y = np.meshgrid(a, a, indexing="ij")
    return np.column_stack([X.ravel(), Y.ravel()])


class Panel:
    """One figure, one image, one colorbar, reused for every frame."""

    def __init__(self, L, name, label, cmap, clim, px, streams, seeds, density):
        dpi = 100.0
        self.name, self.L = name, L
        self.streams, self.seeds, self.density = streams, seeds, density
        self.fig = plt.figure(figsize=(px * 1.20 / dpi, (px + 34) / dpi), dpi=dpi)
        self.ax = self.fig.add_axes((0.012, 0.012, 0.80, 0.90))
        self.im = self.ax.imshow(np.zeros((L, L)), origin="lower", interpolation="nearest",
                                 cmap=cmap, vmin=clim[0], vmax=clim[1])
        self.ax.set_xticks([]); self.ax.set_yticks([])
        self.ax.set_xlim(0, L - 1); self.ax.set_ylim(0, L - 1)
        self.ax.set_title(label, fontsize=17, pad=7)
        cax = self.fig.add_axes((0.828, 0.012, 0.038, 0.90))
        cb = self.fig.colorbar(self.im, cax=cax)
        cb.ax.tick_params(labelsize=8)

    def _clear_streams(self):
        """Remove last frame's streamlines FROM THE AXES, not from the returned container.

        streamplot's return value is a trap for animation. Its `.lines` really is the
        LineCollection that was added to the axes, but `.arrows` is a PatchCollection built
        from the arrow patches and never added to anything -- the patches themselves went on
        one by one via axes.add_patch(). Calling .remove() on that container therefore clears
        nothing, and every frame's arrowheads stay on the figure: after 181 frames the panel
        is 46000 accumulated white arrows with the current frame's streamlines buried under
        them. Clearing the axes' own collections and patches is what actually works.
        """
        for coll in list(self.ax.collections):
            coll.remove()
        for patch in list(self.ax.patches):
            patch.remove()

    def draw(self, fr):
        self.im.set_data(panel_data(fr, self.name).T)
        if self.streams:
            self._clear_streams()
            self._assert_cleared()
            ux, uy = fr["ux"], fr["uy"]
            if np.hypot(ux, uy).max() > 0:
                self.ax.streamplot(np.arange(self.L), np.arange(self.L), ux.T, uy.T,
                                   start_points=self.seeds, density=self.density,
                                   linewidth=0.7, arrowsize=0.6, color="white")
        self.fig.canvas.draw()
        arr = np.asarray(self.fig.canvas.buffer_rgba())[..., :3]
        h, w = arr.shape[:2]
        return np.ascontiguousarray(arr[:h - h % 2, :w - w % 2])

    def _assert_cleared(self):
        """After clearing, the axes must hold NOTHING but the image.

        The obvious invariant -- 'the artist count is the same every frame' -- is wrong here,
        and saying so cost a run: streamplot returns a different number of curves each frame
        because the number that integrate successfully depends on the flow, so the count
        legitimately moved 77 -> 87 and a constant-count check fires on healthy frames.

        What must be true is the weaker, actually-meaningful statement: once the previous
        frame has been cleared, nothing of it is left. imshow's AxesImage lives in
        self.ax.images and is reused, so collections and patches should both be empty.
        """
        left = len(self.ax.collections) + len(self.ax.patches)
        if left:
            raise RuntimeError(
                f"panel {self.name!r}: {left} artist(s) survived the clear; the previous "
                f"frame is still on the axes and will accumulate over the whole video")

    def close(self):
        plt.close(self.fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputdir")
    ap.add_argument("outdir")
    ap.add_argument("--px", type=int, default=460)
    ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--bitrate", default="2200k")
    ap.add_argument("--seed-step", type=float, default=24.0,
                    help="spacing of the fixed streamline seed lattice, in lattice units")
    ap.add_argument("--density", type=float, default=1.4)
    ap.add_argument("--no-streamlines", action="store_true")
    a = ap.parse_args()

    import imageio_ffmpeg
    nf = cw.frame_count(a.inputdir)
    oa = cw.loadarchive(a.inputdir)
    L = int(oa.parameters["LX"])
    os.makedirs(a.outdir, exist_ok=True)
    print(f"{os.path.basename(a.inputdir.rstrip('/'))}: {nf} frames, L = {L}", flush=True)

    clim = limits(oa, nf)
    seeds = seed_grid(L, a.seed_step)
    print(f"  streamline seeds: {len(seeds)} fixed points every {a.seed_step:g} cells",
          flush=True)

    panels, writers = {}, {}
    for n, lab, cmap, _ in PANELS:
        panels[n] = Panel(L, n, lab, cmap, clim[n], a.px,
                          streams=(n == "u" and not a.no_streamlines),
                          seeds=seeds, density=a.density)
        writers[n] = None
    try:
        for i in range(nf):
            try:
                fr = cw.load_frame(oa, i)
            except Exception:
                print(f"  stopping at unreadable frame {i}", flush=True)
                break
            for n, _, _, _ in PANELS:
                arr = panels[n].draw(fr)
                if writers[n] is None:
                    h, w = arr.shape[:2]
                    enc = dict(bitrate=a.bitrate) if a.bitrate else dict(quality=7)
                    writers[n] = imageio_ffmpeg.write_frames(
                        os.path.join(a.outdir, f"{n}.mp4"), (w, h), fps=a.fps,
                        codec="libx264", macro_block_size=2, pix_fmt_out="yuv420p", **enc)
                    writers[n].send(None)
                writers[n].send(arr.tobytes())
            if i % 20 == 0:
                print(f"  frame {i}/{nf}", flush=True)
    finally:
        for n, wr in writers.items():
            if wr is not None:
                wr.close()
                p = os.path.join(a.outdir, f"{n}.mp4")
                print(f"  wrote {p}  ({os.path.getsize(p)/1e6:.1f} MB, "
                      f"clim={clim[n][0]:.3g}..{clim[n][1]:.3g})", flush=True)
        for p in panels.values():
            p.close()


if __name__ == "__main__":
    main()
