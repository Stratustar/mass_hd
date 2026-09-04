#!/usr/bin/env python3
"""Six fields side by side, from the full-resolution frames.

    cw_step_fields.py <case_output_dir>... --out DIR [--fps 5] [--scales scales.json]

The video STREAM carries only |u|, P, m and chi. Two of the six fields asked for here --
p_LB and sigma_bulk -- are not in it, and are not in the frames either: these runs were
written with frame-light = 1, which drops sigma_bulk and pressure_lb along with the LB
populations. Both are nonetheless EXACTLY recoverable from what is stored, because they are
algebraic in Q:

    sigma_bulk = 1/2 CC (1 - q^2)^2 ,   q^2 = QQxx^2 + QQyx^2
    p_LB       = P + sigma_bulk

the second because the model defines P = p_LB - sigma_bulk. So nothing has to be re-run;
the six fields are reconstructed from QQxx, QQyx, chi, m and pressure.

WHAT THE FRAME CADENCE COSTS. ninfo = 4000 gives 26 frames over 100000 steps -- 17.7 tau_c
apart. These are comparison clips, not motion: consecutive frames are independent
snapshots, and at 5 fps the clip is ~5 s. The high-cadence stream (0.2 tau_c) exists for
chi, P and m if smooth motion is wanted, but only there, and mixing the two would put the
six panels on different clocks.

COLOUR SCALES ARE SHARED ACROSS EVERY CASE PASSED, computed in one pass before rendering,
so the four tau_m can be compared by eye. Q is greyscale as |Q| = sqrt(q^2), which is 1 in
the ordered state and 0 in a defect core.
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cw_common as cw                    # noqa: E402
from cw_stream import _fit_font, _text_bar  # noqa: E402

# name -> (colormap, label, symmetric-about-zero)
FIELDS = [("chi", "coolwarm", "chi", False),
          ("P", "RdBu_r", "P", True),
          ("p_LB", "RdBu_r", "p_LB", True),
          ("sigma_B", "magma", "sigma_B", False),
          ("Qabs", "gray", "|Q|", False),
          ("m", "viridis", "m", False)]


def frames_of(root):
    """[(t, {field: 2-D array})] for one case, with the two derived fields added."""
    from archive.archive import loadarchive
    oa = loadarchive(root)
    par = cw.read_params(root)
    CC = float(par["CC"])
    ninfo = int(par["ninfo"])
    out = []
    for i in range(oa._nframes):
        fr = cw.load_frame(oa, i)
        q2 = fr["q2"]
        sig = 0.5 * CC * (1.0 - q2) ** 2
        out.append((i * ninfo, {
            "chi": fr["chi"], "m": fr["m"], "P": fr["P"],
            "sigma_B": sig,
            "p_LB": fr["P"] + sig,          # P = p_LB - sigma_bulk, by definition
            "Qabs": np.sqrt(np.clip(q2, 0, None)),
        }))
    return out


def limits(all_frames, pct=99.0):
    """One scale per field over every case, so the clips are comparable."""
    lim = {}
    for name, _cm, _lab, sym in FIELDS:
        vals = np.concatenate([f[1][name][::4, ::4].ravel() for f in all_frames])
        if name in ("chi", "m", "Qabs"):
            lim[name] = (0.0, 1.0)
        elif sym:
            h = float(np.percentile(np.abs(vals), pct))
            lim[name] = (-h, h)
        else:
            lim[name] = (0.0, float(np.percentile(vals, pct)))
    return lim


def render(frames, name, cmap, label, lo, hi, outpath, tau_c, fps):
    import imageio_ffmpeg
    from matplotlib import colormaps
    lut = (colormaps[cmap](np.linspace(0, 1, 256))[:, :3] * 255).astype(np.uint8)
    h, w = frames[0][1][name].shape[1], frames[0][1][name].shape[0]
    widest = f"{label}   t/tau_c = {frames[-1][0]/tau_c:7.2f}"
    font = _fit_font(widest, w)
    bar_h = 22 + (h + 22) % 2
    writer = None
    try:
        for t, f in frames:
            x = (f[name] - lo) / (hi - lo if hi > lo else 1.0)
            idx = np.clip(np.rint(x * 255.0), 0, 255).astype(np.uint8)
            rgb = lut[np.flipud(idx.T)]
            cap = _text_bar(w, bar_h, f"{label}   t/tau_c = {t/tau_c:7.2f}", font)
            frame = np.concatenate([rgb, cap], axis=0)
            if writer is None:
                writer = imageio_ffmpeg.write_frames(
                    outpath, (frame.shape[1], frame.shape[0]), fps=fps,
                    codec="libx264", macro_block_size=2, pix_fmt_out="yuv420p", quality=7)
                writer.send(None)
            writer.send(frame.tobytes())
    finally:
        if writer is not None:
            writer.close()
    return os.path.getsize(outpath)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", nargs="+", help="case OUTPUT dirs (the ones holding frame*.json)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tau-c", type=float, default=225.56045493417156)
    ap.add_argument("--fps", type=int, default=5)
    ap.add_argument("--tag", nargs="*", default=None,
                    help="clip name per case; defaults to the case directory name")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    loaded = []
    for c in a.cases:
        try:
            fr = frames_of(c)
        except Exception as exc:
            print(f"  {c}: SKIPPED ({type(exc).__name__}: {exc})", flush=True)
            continue
        loaded.append((c, fr))
        print(f"  {os.path.basename(c.rstrip('/'))}: {len(fr)} frames", flush=True)
    if not loaded:
        raise RuntimeError("no readable cases")

    lim = limits([f for _c, fr in loaded for f in fr])
    print("\nshared colour scales:")
    for name, _cm, lab, _s in FIELDS:
        print(f"  {lab:>8}: [{lim[name][0]:+.5g}, {lim[name][1]:+.5g}]")
    with open(os.path.join(a.out, "scales.json"), "w") as fh:
        json.dump({k: list(v) for k, v in lim.items()}, fh, indent=1)

    print()
    for i, (c, fr) in enumerate(loaded):
        tag = (a.tag[i] if a.tag and i < len(a.tag)
               else os.path.basename(c.rstrip("/")))
        for name, cmap, lab, _s in FIELDS:
            p = os.path.join(a.out, f"fld_{tag}_{name}.mp4")
            nb = render(fr, name, cmap, lab, lim[name][0], lim[name][1], p, a.tau_c, a.fps)
            print(f"  {tag}/{name}: {nb/1e6:.2f} MB", flush=True)
    print(f"\nwrote {len(loaded)*len(FIELDS)} clips to {a.out}")


if __name__ == "__main__":
    main()
