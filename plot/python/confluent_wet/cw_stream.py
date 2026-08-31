"""The confluent-wet VIDEO STREAM: reader, de-quantiser and renderer.

The C++ model writes two output channels on two different clocks. `frame*.json` carries the
full-resolution fields at the ninfo cadence (the campaign uses 5 tau_c, ~40 frames per run,
~40 MB each); `video_{u,P,m,chi}.u8` carry four block-averaged uint8 fields at the nvideo
cadence (0.2 tau_c, ~10^3 frames per run, 160 kB each). This module owns the second one.

THE STREAM IS NOT ONLY A PICTURE. It is also the only time axis fine enough for the
temporal statistics: a lag of 3(tau_m + tau_chi) at tau_m = 0.3 tau_c is 2.7 tau_c, which
the 5 tau_c frame cadence cannot resolve at all -- it is shorter than one sample. So the
autocorrelation time of chi and the lagged chi/P cross-correlation are measured here, and
only the spatial statistics (correlation lengths, block variances, defect counts) come from
the full frames. Two consequences to keep in mind and to check rather than assume:

  * QUANTISATION. 256 levels across the stored range. For chi and m that is 1/255 = 0.004,
    far below any structure of interest. For P the stored range is deliberately wider than
    the displayed one (typically 6 sigma_P stored, 3 sigma_P displayed), so the resolution
    across the interesting range is ~128 levels, i.e. 0.05 sigma_P.
  * CLIPPING. Values outside the stored range saturate at byte 0 or 255. The C++ writes the
    per-frame clipped fraction into video_meta.csv precisely so this is visible as a number;
    `Stream.clip_fraction()` surfaces it and cw_part turns it into a warning.

The per-frame SCALARS in video_meta.csv (mean and std of chi and m, u_rms, mean and std of
P) are computed in the C++ in double precision on the FULL-resolution field. They are
therefore exact, and every time series that only needs a domain average -- <chi>(t) above
all -- must be taken from there rather than from the quantised, block-averaged pixels.
"""
import os

import numpy as np

FIELDS = ("u", "P", "m", "chi")

# Display windows, as a multiple of the campaign scale, and the colour map. The stored range
# is wider (see the module docstring); these are what the colour bar spans.
DISPLAY = {
    "u":   ("magma",    "|u|"),
    "P":   ("RdBu_r",   "P"),
    "m":   ("viridis",  "m"),
    "chi": ("coolwarm", "chi"),
}


class Stream:
    """One run's video stream, memory-mapped.

    Frames are never all resident: at 10^3 frames x 400 x 400 each field is 160 MB, and four
    of them plus a float copy would be 2.5 GB. np.memmap keeps that on disk and lets the
    statistics walk it in slices.
    """

    def __init__(self, root, params):
        self.root = root
        self.p = params
        self.stride = int(params.get("video_stride", 2))
        self.nvideo = int(params.get("nvideo", 0))
        if self.nvideo <= 0:
            raise RuntimeError(f"{root}: nvideo = {self.nvideo}, this run has no video "
                               f"stream (it was written by pre-2026-09 code, or nvideo was "
                               f"left at 0 in the runcard)")
        self.nx = int(params["LX"]) // self.stride
        self.ny = int(params["LY"]) // self.stride
        self.p_scale = float(params.get("video_p_scale", 0.0))
        self.u_scale = float(params.get("video_u_scale", 0.0))

        self.meta = self._read_meta()
        # A killed job leaves a partial last frame and a flushed meta line that describes
        # it; trust the byte count, which is the thing that can be short.
        npix = self.nx * self.ny
        counts = [os.path.getsize(self._path(f)) // npix for f in FIELDS]
        self.n = int(min(min(counts), len(self.meta["t"])))
        if self.n < 2:
            raise RuntimeError(f"{root}: video stream has {self.n} complete frames")
        self.steps = np.asarray(self.meta["t"][:self.n], float)

    # ------------------------------------------------------------------ files
    def _path(self, name):
        return os.path.join(self.root, f"video_{name}.u8")

    def _read_meta(self):
        path = os.path.join(self.root, "video_meta.csv")
        cols, rows = None, []
        with open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                if cols is None:
                    cols = line.strip().split(",")
                    continue
                parts = line.strip().split(",")
                if len(parts) == len(cols):        # skip a truncated final line
                    rows.append([float(x) for x in parts])
        if not rows:
            raise RuntimeError(f"{path}: no data rows")
        arr = np.asarray(rows, float)
        return {c: arr[:, i] for i, c in enumerate(cols)}

    def raw(self, name):
        """The uint8 frames as a (n, nx, ny) memmap, x-major to match plot.grid()."""
        return np.memmap(self._path(name), dtype=np.uint8, mode="r",
                         shape=(self.n, self.nx, self.ny))

    # ------------------------------------------------------------ conversion
    def limits(self, name):
        """(lo, hi) of the STORED range in physical units."""
        if name == "u":
            return 0.0, self.u_scale
        if name == "P":
            return -self.p_scale, self.p_scale
        return 0.0, 1.0

    def dequant(self, name, sl=slice(None)):
        """Physical values of a slice of frames, float32."""
        lo, hi = self.limits(name)
        return lo + self.raw(name)[sl].astype(np.float32) * ((hi - lo) / 255.0)

    def clip_fraction(self, i0=0):
        """Worst per-frame clipped area fraction of P and |u| from step index i0 on."""
        return {"P": float(np.max(self.meta["P_clip"][i0:self.n])),
                "u": float(np.max(self.meta["u_clip"][i0:self.n]))}

    def window(self, t_start):
        """Frame indices with step >= t_start."""
        return np.nonzero(self.steps >= t_start)[0]


# --------------------------------------------------------------------- render

def _fit_font(text, width, start=13):
    """Largest default font at which `text` still fits in `width` pixels.

    A caption that overflows is silently truncated by PIL, and the part that gets cut is the
    right-hand end -- which is where <chi> is. Shrinking instead of truncating keeps the
    number visible on a small lattice (the smoke tests run at 50 px) and changes nothing at
    the production 400 px.
    """
    from PIL import ImageDraw, Image, ImageFont
    d = ImageDraw.Draw(Image.new("RGB", (8, 8)))
    for size in range(start, 5, -1):
        try:
            f = ImageFont.load_default(size=size)
        except TypeError:                              # Pillow < 10.1: one fixed size
            return ImageFont.load_default()
        if d.textlength(text, font=f) <= width - 10:
            return f
    return ImageFont.load_default(size=6)


def _text_bar(width, height, line, font):
    from PIL import Image, ImageDraw
    img = Image.new("RGB", (width, height), (16, 16, 16))
    d = ImageDraw.Draw(img)
    d.text((6, height // 2), line, fill=(235, 235, 235), font=font, anchor="lm")
    return np.asarray(img, dtype=np.uint8)


def render(stream, name, outpath, disp_lo, disp_hi, tau_c, fps=25, bitrate=None,
           label=None, progress=None):
    """Encode one field to mp4 with a fixed colour scale and a per-frame caption.

    disp_lo/disp_hi are in PHYSICAL units and are normally the campaign-wide convention
    (chi, m in [0,1]; P in +/- 3 sigma_P(zeta); |u| in [0, 3 u_rms(zeta)]), NOT this run's
    own range. A per-run colour scale animates a field that grows by decades as though it
    were constant and makes two runs incomparable by eye, which is the one thing the videos
    exist to allow.

    The caption carries t/tau_c(zeta) and <chi>, both read from video_meta.csv, i.e. exact
    domain averages of the full-resolution field rather than of these pixels.
    """
    import imageio_ffmpeg
    from matplotlib import colormaps

    cmap, sym = DISPLAY[name]
    lut = (colormaps[cmap](np.linspace(0, 1, 256))[:, :3] * 255).astype(np.uint8)

    lo, hi = stream.limits(name)
    raw = stream.raw(name)
    chi_mean = stream.meta["chi_mean"]

    # size the font once, on the widest caption the run can produce
    widest = (f"{label or sym}   t/tau_c = {stream.steps[-1]/tau_c:7.2f}   <chi> = 0.000")
    font = _fit_font(widest, stream.nx)

    bar_h = 22
    h, w = stream.ny, stream.nx                    # image rows = y, columns = x
    bar_h += (h + bar_h) % 2                       # keep the frame height even for x264

    writer = None
    try:
        for i in range(stream.n):
            # stored byte -> physical -> position in the DISPLAY window -> colour index
            phys = lo + raw[i].astype(np.float32) * ((hi - lo) / 255.0)
            x = (phys - disp_lo) / (disp_hi - disp_lo)
            idx = np.clip(np.rint(x * 255.0), 0, 255).astype(np.uint8)
            rgb = lut[np.flipud(idx.T)]            # (ny, nx, 3), y increasing upwards
            cap = _text_bar(w, bar_h,
                            f"{label or sym}   t/tau_c = {stream.steps[i]/tau_c:7.2f}"
                            f"   <chi> = {chi_mean[i]:.3f}", font)
            frame = np.concatenate([rgb, cap], axis=0)
            if writer is None:
                enc = dict(bitrate=bitrate) if bitrate else dict(quality=7)
                writer = imageio_ffmpeg.write_frames(
                    outpath, (frame.shape[1], frame.shape[0]), fps=fps,
                    codec="libx264", macro_block_size=2, pix_fmt_out="yuv420p", **enc)
                writer.send(None)
            writer.send(frame.tobytes())
            if progress and i % 200 == 0:
                print(f"    {name}: frame {i}/{stream.n}", flush=True)
    finally:
        if writer is not None:
            writer.close()
    return os.path.getsize(outpath)
