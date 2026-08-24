"""Extra field views for the confluent-memory project, as PNG frames and mp4.

plot_hd.py renders visualization / pressure / director. This adds the two views the
phenotype-blob question needs and that the library does not provide:

  chi_velocity  the phenotype field as a DARK two-tone background with the velocity
                drawn on top as bright arrows. The background is deliberately low
                luminance so the arrows stay readable; a normal viridis/coolwarm map
                is too bright for a quiver overlay.
  memory        the mechanical memory m on its own.

It also records the flocking diagnostics per frame, since "are the arrows aligned over
a long range" is the same question as "is there flocking":

  |<v>| / v_rms    polar order parameter: ~1/sqrt(N) for uncorrelated flow, -> 1 for a
                   uniformly drifting (flocked) state
  l_vv             velocity correlation length (first zero crossing of C_vv)

Usage:
  cm_fields.py <study> [--cases z3_tc1_tm30 ...] [--stride 2] [--skip-video] [--out stats.json]
"""
import argparse, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from cm_common import SCRATCH, frame_count, iter_frames, Radial, parse_variant, ph

MP4_FPS = 10
QUIVER_STRIDE = 12          # one arrow per 12 lattice sites; denser is unreadable at 400^2

# Two dark, well-separated hues for the phenotype: deep teal = go (chi -> 0),
# deep magenta = grow (chi -> 1). Both are kept at low luminance on purpose.
CHI_CMAP = LinearSegmentedColormap.from_list(
    "chi_dark", ["#0d3b4a", "#14606f", "#2a2a3a", "#5c1f4e", "#7d1f52"])
M_CMAP = "magma"


def _fig(L):
    fig = plt.figure(figsize=(6.4, 6.4), dpi=110)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, L); ax.set_ylim(0, L); ax.set_axis_off()
    return fig, ax


def _to_rgb(fig):
    fig.canvas.draw()
    arr = np.asarray(fig.canvas.buffer_rgba())[..., :3]
    return np.ascontiguousarray(arr, dtype=np.uint8)


def render_chi_velocity(fr, L, vscale=None):
    chi = ph.grid(fr, "chi"); ux = ph.grid(fr, "ux"); uy = ph.grid(fr, "uy")
    fig, ax = _fig(L)
    ax.imshow(chi.T, origin="lower", extent=[0, L, 0, L], cmap=CHI_CMAP,
              vmin=0., vmax=1., interpolation="nearest")
    s = QUIVER_STRIDE
    xx, yy = np.meshgrid(np.arange(0, L, s) + s / 2., np.arange(0, L, s) + s / 2.,
                         indexing="ij")
    u = ux[::s, ::s]; v = uy[::s, ::s]
    # one common scale per case so arrow length is comparable across frames
    sc = vscale if vscale else max(np.sqrt(u**2 + v**2).max(), 1e-12)
    ax.quiver(xx, yy, u, v, color="#ffe9a8", angles="xy", scale_units="xy",
              scale=sc / (0.9 * s), width=0.0028, headwidth=3.4, headlength=4.2,
              headaxislength=3.6, alpha=.95)
    return fig


def render_memory(fr, L):
    m = ph.grid(fr, "m")
    fig, ax = _fig(L)
    ax.imshow(m.T, origin="lower", extent=[0, L, 0, L], cmap=M_CMAP,
              vmin=0., vmax=1., interpolation="nearest")
    return fig


RENDER = {"chi_velocity": render_chi_velocity, "memory": render_memory}


def flocking(ux, uy):
    """(|<v>|/v_rms, v_rms). The ratio is the polar order parameter."""
    vrms = float(np.sqrt(np.mean(ux**2 + uy**2)))
    vbar = float(np.hypot(ux.mean(), uy.mean()))
    return (vbar / vrms if vrms > 0 else np.nan), vrms


def process(root, outdir, L, stride, skip_video):
    import imageio_ffmpeg
    os.makedirs(outdir, exist_ok=True)
    radial = Radial(L)
    frames = [(i, t, fr) for i, t, fr in iter_frames(root)]
    if not frames:
        return None
    # common arrow scale from the median frame, so arrows are comparable across the run
    mid = frames[len(frames) // 2][2]
    umid, vmid = ph.grid(mid, "ux"), ph.grid(mid, "uy")
    vscale = max(float(np.percentile(np.hypot(umid, vmid), 99)), 1e-12)

    stats = {"t": [], "polar": [], "vrms": [], "l_vv": []}
    writers, sel = {}, frames[::stride]
    still_at = {0, len(frames) // 2, len(frames) - 1}
    try:
        for n, (i, t, fr) in enumerate(frames):
            ux, uy = ph.grid(fr, "ux"), ph.grid(fr, "uy")
            p, vr = flocking(ux, uy)
            c = radial.mean(radial.corr([ux, uy]))
            c = c / c[0] if c[0] > 0 else c * np.nan
            stats["t"].append(t); stats["polar"].append(p); stats["vrms"].append(vr)
            stats["l_vv"].append(radial.length(c))
            want_still = i in still_at
            want_video = (not skip_video) and (frames[::stride] and (i, t, fr) in sel)
            if not (want_still or want_video):
                continue
            for name, fn in RENDER.items():
                fig = fn(fr, L, vscale) if name == "chi_velocity" else fn(fr, L)
                if want_still:
                    fig.savefig(os.path.join(outdir, f"{name}_frame{i}.png"), dpi=110)
                if want_video:
                    arr = _to_rgb(fig)
                    if name not in writers:
                        h, w = arr.shape[:2]
                        out = os.path.join(outdir, f"{name}_0-{len(frames)-1}_step{stride}.mp4")
                        writers[name] = imageio_ffmpeg.write_frames(
                            out, (w, h), fps=MP4_FPS, codec="libx264", quality=8,
                            macro_block_size=2, pix_fmt_out="yuv420p")
                        writers[name].send(None)
                    writers[name].send(arr.tobytes())
                plt.close(fig)
            if n % 20 == 0:
                print(f"  frame {i}/{len(frames)-1}", flush=True)
    finally:
        for w in writers.values():
            w.close()
    return stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("study")
    ap.add_argument("--cases", nargs="*", default=None)
    ap.add_argument("--stride", type=int, default=2)
    ap.add_argument("--L", type=int, default=400)
    ap.add_argument("--skip-video", action="store_true")
    ap.add_argument("--out", default=None)
    ap.add_argument("--steady", type=float, default=None)
    args = ap.parse_args()

    base = os.path.join(SCRATCH, "cases", args.study)
    results_base = os.path.join(SCRATCH, "results", "cases", args.study)
    cases = args.cases or sorted(os.listdir(base))
    out = []
    for case in cases:
        root = os.path.join(base, case)
        if not os.path.isdir(root) or frame_count(root) < 20:
            print(f"skip {case}", flush=True)
            continue
        print(f"=== {case} ===", flush=True)
        st = process(root, os.path.join(results_base, case), args.L,
                     args.stride, args.skip_video)
        if st is None:
            continue
        steady = args.steady if args.steady is not None else st["t"][-1] / 2.
        sel = [j for j, t in enumerate(st["t"]) if t >= steady]
        rec = dict(study=args.study, case=case, ts=st, steady=steady, **parse_variant(case))
        for k in ("polar", "vrms", "l_vv"):
            rec["mean_" + k] = float(np.nanmean([st[k][j] for j in sel]))
        out.append(rec)
        print(f"  polar={rec['mean_polar']:.4f}  l_vv={rec['mean_l_vv']:.2f}", flush=True)
    if args.out:
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        json.dump(out, open(args.out, "w"))
        print(f"wrote {len(out)} cases -> {args.out}")


if __name__ == "__main__":
    main()
