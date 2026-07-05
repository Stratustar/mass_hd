import matplotlib
matplotlib.use("Agg")

import argparse
import os
import sys

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from PIL import Image

DIRECTOR_STRIDE = 3
VELOCITY_STRIDE = 5
DIRECTOR_WIDTH = 0.0015
VISUALIZATION_CMAP = plt.get_cmap("bwr")  # go (chi=0)=blue, grow (chi=1)=red, boundary (chi=0.5)=white
VISUALIZATION_BACKGROUND = "#7e7e7e"
GIF_DPI = 300  # default plotting spec: video frames rendered at 300 dpi
MP4_FPS = 10   # default plotting spec: animations written as mp4 (h264) at this fps
# Black-and-white cyclic map for the nematic director ANGLE. theta is pi-periodic, so a
# monotonic grayscale would put a false black/white seam at 0/180 deg (same orientation);
# instead we go black -> white -> black across 0..180 deg so the wrap is seamless.
DIRECTOR_ANGLE_CMAP = mcolors.LinearSegmentedColormap.from_list(
    "cyclic_gray", ["black", "white", "black"])

BASE_DIR = os.path.abspath(
    os.path.join(
        os.path.dirname(__file__),
        "..",
        "mass_hd",
        "plot",
        "python",
    )
)
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)

from archive.archive import loadarchive


def require_frame_fields(frame, field_name, required):
    missing = [name for name in required if not hasattr(frame, name)]
    if missing:
        available = sorted(key for key in frame.__dict__.keys() if not key.startswith("_"))
        raise AttributeError(
            f"Frame needs fields {missing} to compute '{field_name}'. "
            f"Available fields: {available}"
        )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot m, chi, phi, and pressure fields from a MASS archive."
    )
    parser.add_argument("inputdir", help="Archive directory, for example ./output")
    parser.add_argument("outdir", help="Directory for generated figures")
    parser.add_argument(
        "--frames",
        nargs="+",
        type=int,
        default=None,
        help="Frame indices to export as PNG. Default: first, middle and last frame.",
    )
    parser.add_argument(
        "--fields",
        nargs="+",
        default=None,
        help=(
            "Field names to render. Default: all standard fields. "
            "Use '--fields visualization' for only the chi-phi composite."
        ),
    )
    parser.add_argument(
        "--skip-png",
        action="store_true",
        help="Skip PNG export and only render GIFs.",
    )
    parser.add_argument(
        "--gif-frames",
        type=int,
        default=None,
        help="Max frames per GIF. Default: auto (all available).",
    )
    parser.add_argument(
        "--gif-frame-start",
        type=int,
        default=0,
        help="First frame index to render into each GIF. Default: 0",
    )
    parser.add_argument(
        "--gif-frame-step",
        type=int,
        default=None,
        help="Frame index step for GIF. Default: auto (~40 frames over the run).",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=600,
        help="PNG output resolution. Default: 600",
    )
    parser.add_argument(
        "--skip-counts",
        action="store_true",
        help="Skip the N0/N1 time-series plot.",
    )
    return parser.parse_args()


def field_as_grid(frame, field_name):
    if field_name == "phi_one_minus_chi":
        if not hasattr(frame, "chi") or not hasattr(frame, "phi"):
            available = sorted(
                key for key in frame.__dict__.keys() if not key.startswith("_")
            )
            raise AttributeError(
                "Frame needs fields 'chi' and 'phi' to compute "
                f"'phi_one_minus_chi'. Available fields: {available}"
            )
        lx = frame.parameters["LX"]
        ly = frame.parameters["LY"]
        chi = np.array(frame.chi).reshape((lx, ly))
        phi = np.array(frame.phi).reshape((lx, ly))
        return phi * (1.0 - chi)

    if field_name == "chi_phi":
        if not hasattr(frame, "chi") or not hasattr(frame, "phi"):
            available = sorted(
                key for key in frame.__dict__.keys() if not key.startswith("_")
            )
            raise AttributeError(
                "Frame needs fields 'chi' and 'phi' to compute 'chi_phi'. "
                f"Available fields: {available}"
            )
        lx = frame.parameters["LX"]
        ly = frame.parameters["LY"]
        chi = np.array(frame.chi).reshape((lx, ly))
        phi = np.array(frame.phi).reshape((lx, ly))
        return chi * phi

    if field_name == "pressure" and not hasattr(frame, "pressure"):
        if hasattr(frame, "sigma_bulk"):
            lx = frame.parameters["LX"]
            ly = frame.parameters["LY"]
            return -np.array(frame.sigma_bulk).reshape((lx, ly))

    if not hasattr(frame, field_name):
        available = sorted(
            key for key in frame.__dict__.keys() if not key.startswith("_")
        )
        raise AttributeError(
            f"Frame has no field '{field_name}'. Available fields: {available}"
        )

    lx = frame.parameters["LX"]
    ly = frame.parameters["LY"]
    return np.array(getattr(frame, field_name)).reshape((lx, ly))


def srgb_to_linear(rgb):
    rgb = np.asarray(rgb, dtype=float)
    return np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)


def linear_to_srgb(rgb):
    rgb = np.asarray(rgb, dtype=float)
    rgb = np.clip(rgb, 0.0, 1.0)
    return np.where(rgb <= 0.0031308, 12.92 * rgb, 1.055 * rgb ** (1.0 / 2.4) - 0.055)


def xyz_to_lab(xyz):
    white = np.array([0.95047, 1.0, 1.08883])
    scaled = xyz / white
    delta = 6.0 / 29.0
    f = np.where(
        scaled > delta**3,
        np.cbrt(scaled),
        scaled / (3.0 * delta**2) + 4.0 / 29.0,
    )
    return np.stack(
        (
            116.0 * f[..., 1] - 16.0,
            500.0 * (f[..., 0] - f[..., 1]),
            200.0 * (f[..., 1] - f[..., 2]),
        ),
        axis=-1,
    )


def lab_to_xyz(lab):
    white = np.array([0.95047, 1.0, 1.08883])
    fy = (lab[..., 0] + 16.0) / 116.0
    fx = fy + lab[..., 1] / 500.0
    fz = fy - lab[..., 2] / 200.0
    f = np.stack((fx, fy, fz), axis=-1)
    delta = 6.0 / 29.0
    xyz_scaled = np.where(f > delta, f**3, 3.0 * delta**2 * (f - 4.0 / 29.0))
    return xyz_scaled * white


def srgb_to_lab(rgb):
    linear = srgb_to_linear(rgb)
    xyz = linear @ np.array(
        [
            [0.4124564, 0.3575761, 0.1804375],
            [0.2126729, 0.7151522, 0.0721750],
            [0.0193339, 0.1191920, 0.9503041],
        ]
    ).T
    return xyz_to_lab(xyz)


def lab_to_srgb(lab):
    xyz = lab_to_xyz(lab)
    linear = xyz @ np.array(
        [
            [3.2404542, -1.5371385, -0.4985314],
            [-0.9692660, 1.8760108, 0.0415560],
            [0.0556434, -0.2040259, 1.0572252],
        ]
    ).T
    return linear_to_srgb(linear)


def visualization_rgba(frame):
    if not hasattr(frame, "chi") or not hasattr(frame, "phi"):
        available = sorted(key for key in frame.__dict__.keys() if not key.startswith("_"))
        raise AttributeError(
            "Frame needs fields 'chi' and 'phi' to compute 'visualization'. "
            f"Available fields: {available}"
        )

    chi = np.clip(field_as_grid(frame, "chi"), 0.0, 1.0)
    phi = np.clip(field_as_grid(frame, "phi"), 0.0, 1.0)

    rgb = VISUALIZATION_CMAP(chi)[..., :3]

    return np.dstack((rgb, phi))


def visualization_counts(frame):
    chi = field_as_grid(frame, "chi")
    phi = field_as_grid(frame, "phi")
    n0 = np.sum(phi * (1.0 - chi))
    n1 = np.sum(phi * chi)
    return n0, n1


def annotate_visualization_counts(ax, frame):
    n0, n1 = visualization_counts(frame)
    ax.text(
        1.05,
        0.5,
        f"N0= {n0:.6g}\nN1= {n1:.6g}",
        transform=ax.transAxes,
        ha="left",
        va="center",
        fontsize=12,
        clip_on=False,
    )


def director_components(frame, stride=DIRECTOR_STRIDE):
    if not hasattr(frame, "QQxx") or not hasattr(frame, "QQyx"):
        available = sorted(
            key for key in frame.__dict__.keys() if not key.startswith("_")
        )
        raise AttributeError(
            "Frame needs fields 'QQxx' and 'QQyx' to compute 'director'. "
            f"Available fields: {available}"
        )

    lx = frame.parameters["LX"]
    ly = frame.parameters["LY"]
    qxx = np.array(frame.QQxx).reshape((lx, ly))
    qyx = np.array(frame.QQyx).reshape((lx, ly))
    theta = 0.5 * np.arctan2(qyx, qxx)
    strength = np.sqrt(qxx * qxx + qyx * qyx)

    xs = np.arange(0, lx, stride)
    ys = np.arange(0, ly, stride)
    xx, yy = np.meshgrid(xs, ys, indexing="ij")
    order = np.clip(strength[xx, yy], 0.0, 1.0)
    uu = np.cos(theta[xx, yy]) * order
    vv = np.sin(theta[xx, yy]) * order
    return xx, yy, uu, vv, strength


def velocity_components(frame, stride=VELOCITY_STRIDE):
    if not hasattr(frame, "ff"):
        available = sorted(
            key for key in frame.__dict__.keys() if not key.startswith("_")
        )
        raise AttributeError(
            "Frame needs field 'ff' to compute velocity. "
            f"Available fields: {available}"
        )

    lx = frame.parameters["LX"]
    ly = frame.parameters["LY"]
    ff = np.array(frame.ff)
    density = np.sum(ff, axis=1)
    density = np.where(density == 0.0, np.nan, density)
    vx = (ff[:, 1] - ff[:, 2] + ff[:, 5] - ff[:, 6] - ff[:, 7] + ff[:, 8]) / density
    vy = (ff[:, 3] - ff[:, 4] + ff[:, 5] - ff[:, 6] + ff[:, 7] - ff[:, 8]) / density
    vx = np.nan_to_num(vx).reshape((lx, ly))
    vy = np.nan_to_num(vy).reshape((lx, ly))
    speed = np.sqrt(vx * vx + vy * vy)
    max_speed = np.max(speed)

    xs = np.arange(0, lx, stride)
    ys = np.arange(0, ly, stride)
    xx, yy = np.meshgrid(xs, ys, indexing="ij")
    if max_speed > 0.0:
        uu = vx[xx, yy] / max_speed
        vv = vy[xx, yy] / max_speed
    else:
        uu = np.zeros_like(xx, dtype=float)
        vv = np.zeros_like(yy, dtype=float)
    return xx, yy, uu, vv, speed


def draw_boundary_contours(ax, frame, lw_scale=1.0):
    """Overlay phi=0.5 (black dashed) and chi=0.5 within phi>0.5 (black solid).

    lw_scale scales both line widths (1.0 for visualization, 0.5 for the
    pressure/velocity/director field panels).
    """
    if not (hasattr(frame, "phi") and hasattr(frame, "chi")):
        return
    phi = field_as_grid(frame, "phi")
    chi = field_as_grid(frame, "chi")
    ax.contour(phi.T, levels=[0.5], colors=["black"], linestyles="dashed",
               linewidths=0.8 * lw_scale, origin="lower")
    chi_interior = np.ma.masked_where(phi <= 0.5, chi)
    ax.contour(chi_interior.T, levels=[0.5], colors=["black"],
               linewidths=1.0 * lw_scale, origin="lower")


def draw_defects(ax, frame, size_scale=1.0):
    """Overlay +/-1/2 nematic defect markers: circle = +1/2, triangle = -1/2.

    Defects are detected on the full-resolution director (theta = 1/2 atan2(QQyx,
    QQxx)) inside the colony (phi>0.5 mask) by the winding method in
    plot/python/plot/defects.py. Markers are drawn in lattice/plot coordinates so
    they sit on the imshow(field.T)+quiver panels. Silently skips if the Q fields
    or scipy are unavailable.
    """
    if not (hasattr(frame, "QQxx") and hasattr(frame, "QQyx")):
        return
    try:
        from plot.defects import defect_markers
    except Exception as exc:  # scipy/plot package missing
        print(f"defect markers unavailable ({exc}); skipping.", flush=True)
        return
    qxx = field_as_grid(frame, "QQxx")
    qyx = field_as_grid(frame, "QQyx")
    mask = field_as_grid(frame, "phi") > 0.5 if hasattr(frame, "phi") else None
    plus, minus = defect_markers(qxx, qyx, mask=mask)
    if plus:
        px, py = zip(*plus)
        ax.scatter(px, py, marker="o", c="#00e5ff", edgecolors="black",
                   s=30 * size_scale, linewidths=0.5, zorder=6)
    if minus:
        mx, my = zip(*minus)
        ax.scatter(mx, my, marker="^", c="#ff00ff", edgecolors="black",
                   s=42 * size_scale, linewidths=0.5, zorder=6)


def director_angle_deg_masked(frame):
    """Nematic director angle in degrees [0, 180), masked to the colony (phi>0.5).

    theta = 1/2 atan2(QQyx, QQxx) is pi-periodic; returned in degrees on [0,180) so a
    cyclic colormap (DIRECTOR_ANGLE_CMAP) with vmin=0/vmax=180 colours orientation
    continuously (0 deg == 180 deg). Outside the colony the director is vacuum noise,
    so those nodes are masked (shown as the neutral background), matching the field
    convention of colouring orientation only where there is material.
    """
    qxx = field_as_grid(frame, "QQxx")
    qyx = field_as_grid(frame, "QQyx")
    deg = np.degrees(0.5 * np.arctan2(qyx, qxx) % np.pi)  # [0, 180)
    if hasattr(frame, "phi"):
        deg = np.ma.masked_where(field_as_grid(frame, "phi") <= 0.5, deg)
    return deg


def draw_field(ax, frame, field_name, cmap, clim=None, show_defects=True):
    if field_name == "visualization":
        rgba = visualization_rgba(frame)
        ax.set_facecolor(VISUALIZATION_BACKGROUND)
        im = ax.imshow(
            np.transpose(rgba, (1, 0, 2)),
            origin="lower",
            interpolation="nearest",
        )
        draw_boundary_contours(ax, frame, 1.0)
        if show_defects:
            draw_defects(ax, frame)
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        annotate_visualization_counts(ax, frame)
        im.skip_colorbar = True
        return im

    if field_name in ("director", "phi_director"):
        xx, yy, uu, vv, strength = director_components(frame)
        # Background: director ORIENTATION via a cyclic colormap (masked to the colony),
        # overlaid with the black director glyphs and the +/-1/2 defect markers.
        ax.set_facecolor(VISUALIZATION_BACKGROUND)
        im = ax.imshow(
            director_angle_deg_masked(frame).T,
            origin="lower",
            interpolation="nearest",
            cmap=DIRECTOR_ANGLE_CMAP,
            vmin=0.0,
            vmax=180.0,
        )
        ax.quiver(
            xx,
            yy,
            uu,
            vv,
            pivot="middle",
            headwidth=0,
            headlength=0,
            headaxislength=0,
            angles="xy",
            scale_units="xy",
            scale=1.0 / DIRECTOR_STRIDE,
            width=DIRECTOR_WIDTH,
            color="black",
        )
        draw_boundary_contours(ax, frame, 0.5)
        if show_defects:
            draw_defects(ax, frame)
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        im.colorbar_label = "director angle (deg)"
        return im

    if field_name == "phi_velocity":
        xx, yy, uu, vv, speed = velocity_components(frame)
        # Background: director orientation (same cyclic map as the director panel),
        # overlaid with the red velocity arrows and the +/-1/2 defect markers.
        ax.set_facecolor(VISUALIZATION_BACKGROUND)
        im = ax.imshow(
            director_angle_deg_masked(frame).T,
            origin="lower",
            interpolation="nearest",
            cmap=DIRECTOR_ANGLE_CMAP,
            vmin=0.0,
            vmax=180.0,
        )
        ax.quiver(
            xx,
            yy,
            uu,
            vv,
            pivot="middle",
            scale=28,
            width=0.003,
            color="red",
        )
        draw_boundary_contours(ax, frame, 0.5)
        if show_defects:
            draw_defects(ax, frame)
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        im.colorbar_label = "director angle (deg)"
        im.max_speed = np.max(speed)
        return im

    field = field_as_grid(frame, field_name)

    im = ax.imshow(
        field.T,
        origin="lower",
        interpolation="nearest",
        cmap=cmap,
        vmin=None if clim is None else clim[0],
        vmax=None if clim is None else clim[1],
    )
    if field_name == "pressure":
        draw_boundary_contours(ax, frame, 0.5)
        if show_defects:
            draw_defects(ax, frame)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    return im


def available_frame_count(ar):
    return int((ar.nsteps - ar.nstart) / ar.ninfo) + 1


def compute_crop_box(ar, margin=30):
    """View window: bbox of phi>0.5 at the LAST frame, expanded by `margin` per
    side (clamped to the lattice). Returns (x0, x1, y0, y1) in plot coords, or None.
    Applied to every field panel so all frames share one colony-tight view."""
    n = available_frame_count(ar)
    if n < 1:
        return None
    frame = ar.read_frame(n - 1)
    if not hasattr(frame, "phi"):
        return None
    phi = field_as_grid(frame, "phi")
    lx, ly = phi.shape
    mask = phi > 0.5
    if not mask.any():
        return None
    xs = np.where(mask.any(axis=1))[0]
    ys = np.where(mask.any(axis=0))[0]
    x0 = max(0, int(xs.min()) - margin)
    x1 = min(lx - 1, int(xs.max()) + margin)
    y0 = max(0, int(ys.min()) - margin)
    y1 = min(ly - 1, int(ys.max()) + margin)
    return (x0, x1, y0, y1)


def save_field_png(ar, frame_index, field_name, outdir, cmap, clim, dpi, crop=None,
                   show_defects=True, name_suffix=""):
    frame = ar.read_frame(frame_index)

    fig, ax = plt.subplots(figsize=(6, 5))
    im = draw_field(ax, frame, field_name, cmap=cmap, clim=clim, show_defects=show_defects)
    if crop is not None:
        ax.set_xlim(crop[0], crop[1])
        ax.set_ylim(crop[2], crop[3])
    ax.set_title(f"{field_name}, frame {frame_index}")
    if not getattr(im, "skip_colorbar", False):
        fig.colorbar(im, ax=ax, label=getattr(im, "colorbar_label", field_name))
    fig.tight_layout()
    if field_name == "visualization":
        fig.subplots_adjust(right=0.75)

    outfile = os.path.join(outdir, f"{field_name}{name_suffix}_frame{frame_index}.png")
    fig.savefig(
        outfile,
        dpi=dpi,
        transparent=getattr(im, "transparent_output", False),
        bbox_inches="tight" if field_name == "visualization" else None,
    )
    plt.close(fig)
    return outfile


def render_field_image(ar, frame_index, field_name, cmap, clim, crop=None, show_defects=True):
    frame = ar.read_frame(frame_index)

    fig, ax = plt.subplots(figsize=(6, 5), dpi=GIF_DPI)
    im = draw_field(ax, frame, field_name, cmap=cmap, clim=clim, show_defects=show_defects)
    if crop is not None:
        ax.set_xlim(crop[0], crop[1])
        ax.set_ylim(crop[2], crop[3])
    ax.set_title(f"{field_name}, frame {frame_index}")
    if not getattr(im, "skip_colorbar", False):
        fig.colorbar(im, ax=ax, label=getattr(im, "colorbar_label", field_name))
    fig.tight_layout()
    if field_name == "visualization":
        fig.subplots_adjust(right=0.75)

    fig.canvas.draw()
    width, height = fig.canvas.get_width_height()
    buf = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8).reshape(
        height, width, 4
    )
    image = Image.fromarray(buf[..., :3].copy())
    plt.close(fig)
    return image


def save_field_mp4(ar, field_name, outdir, cmap, clim, nframes, frame_start, frame_step, crop=None,
                   show_defects=True, name_suffix=""):
    """Write the animation directly as an h264 mp4 (default plotting spec).

    mp4 replaces the old gif path: it is far smaller than a 300-dpi gif, needs no
    downstream transcode (the dashboard plays mp4 with an adjustable-speed <video>),
    and streams frame-by-frame so we never hold all frames in memory. Frame
    selection is identical to the previous gif routine.
    """
    import imageio_ffmpeg  # ffmpeg frame-writer (same dep the dashboard already uses)

    if nframes < 1:
        raise RuntimeError(f"Video needs at least one {field_name} frame.")
    if frame_step < 1:
        raise ValueError("--gif-frame-step must be at least 1.")

    n_available = available_frame_count(ar)
    frame_indices = list(range(frame_start, n_available, frame_step))[:nframes]
    final_frame = n_available - 1
    if frame_indices and frame_indices[-1] != final_frame and final_frame >= frame_start:
        frame_indices.append(final_frame)
    if not frame_indices:
        raise RuntimeError(f"No {field_name} frames selected for video.")

    outfile = os.path.join(
        outdir,
        f"{field_name}{name_suffix}_{frame_indices[0]}-{frame_indices[-1]}_step{frame_step}.mp4",
    )
    writer = None
    try:
        for frame_index in frame_indices:
            print(f"Rendering {field_name} video frame {frame_index}", flush=True)
            img = render_field_image(ar, frame_index, field_name, cmap, clim, crop, show_defects=show_defects)
            arr = np.asarray(img)
            if arr.ndim == 3 and arr.shape[2] == 4:
                arr = arr[..., :3]
            arr = np.ascontiguousarray(arr, dtype=np.uint8)
            if writer is None:
                h, w = arr.shape[:2]
                writer = imageio_ffmpeg.write_frames(
                    outfile, (w, h), fps=MP4_FPS, codec="libx264",
                    quality=8, macro_block_size=2, pix_fmt_out="yuv420p",
                )
                writer.send(None)  # seed the generator
            writer.send(arr.tobytes())
    finally:
        if writer is not None:
            writer.close()
    return outfile


def save_counts_timeseries(ar, outdir, dpi):
    n_available = available_frame_count(ar)
    frame_indices = list(range(n_available))
    n0_values = []
    n1_values = []

    for frame_index in frame_indices:
        print(f"Computing N0/N1 for frame {frame_index}", flush=True)
        frame = ar.read_frame(frame_index)
        n0, n1 = visualization_counts(frame)
        n0_values.append(n0)
        n1_values.append(n1)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(frame_indices, n0_values, color="tab:red", label="N0")
    ax.plot(frame_indices, n1_values, color="tab:blue", label="N1")
    ax.set_xlabel("t")
    ax.set_ylabel("N")
    ax.set_title("N0 and N1 vs t")
    ax.legend()
    fig.tight_layout()

    outfile = os.path.join(outdir, "N0_N1_vs_t.png")
    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)
    return outfile


def defect_variants(style):
    """(show_defects, name_suffix) variants to render for a field.

    'never' -> a single variant with no defect markers (velocity).
    'toggle' -> the default (with-defect, unsuffixed) variant plus a '_nodef'
    variant, so the dashboard can switch defects on/off.
    """
    if style.get("defects", "toggle") == "never":
        return [(False, "")]
    return [(True, ""), (False, "_nodef")]


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print(f"Loading archive from: {args.inputdir}", flush=True)
    ar = loadarchive(args.inputdir)
    n_available = available_frame_count(ar)
    print(f"Archive loaded. Available frame count: {n_available}", flush=True)

    crop = compute_crop_box(ar, margin=30)

    # defects: "toggle" -> render a with-defect (default-named) and a nodef variant
    # so the dashboard can switch them; "never" -> single variant with no defect
    # markers at all (velocity).
    fields = {
        "visualization": {"cmap": None, "clim": None, "defects": "toggle"},
        "pressure": {"cmap": "coolwarm", "clim": None, "defects": "toggle"},
        "phi_velocity": {"cmap": None, "clim": None, "defects": "never"},
        "phi_director": {"cmap": None, "clim": None, "defects": "toggle"},
    }

    if args.fields is not None:
        unknown = [field for field in args.fields if field not in fields]
        if unknown:
            known = ", ".join(fields.keys())
            raise SystemExit(f"Unknown plot field(s): {unknown}. Known fields: {known}")
        fields = {field: fields[field] for field in args.fields}

    if args.frames is not None:
        frames_to_plot = args.frames
    elif n_available > 0:
        mid = int(round(0.5 * (n_available - 1)))
        frames_to_plot = sorted({0, mid, n_available - 1})
    else:
        frames_to_plot = []

    if args.skip_png:
        print("Skipping PNG export.", flush=True)
    else:
        print("Starting PNG export...", flush=True)
        for frame_index in frames_to_plot:
            if frame_index >= n_available:
                print(
                    f"Skipping frame {frame_index}; archive only has {n_available} frames.",
                    flush=True,
                )
                continue

            for field_name, style in fields.items():
                for show_def, suffix in defect_variants(style):
                    print(f"Rendering {field_name}{suffix} PNG frame {frame_index}", flush=True)
                    outfile = save_field_png(
                        ar,
                        frame_index,
                        field_name,
                        args.outdir,
                        style["cmap"],
                        style["clim"],
                        args.dpi,
                        crop,
                        show_defects=show_def,
                        name_suffix=suffix,
                    )
                    print(f"Saved: {outfile}", flush=True)

    gif_step = args.gif_frame_step if args.gif_frame_step is not None else max(1, int(round(n_available / 40.0)))
    gif_count = args.gif_frames if args.gif_frames is not None else n_available
    print("Preparing videos...", flush=True)
    for field_name, style in fields.items():
        for show_def, suffix in defect_variants(style):
            outfile = save_field_mp4(
                ar,
                field_name,
                args.outdir,
                style["cmap"],
                style["clim"],
                gif_count,
                args.gif_frame_start,
                gif_step,
                crop,
                show_defects=show_def,
                name_suffix=suffix,
            )
            print(f"Saved video: {outfile}", flush=True)

    if args.skip_counts:
        print("Skipping N0/N1 time-series plot.", flush=True)
    else:
        print("Preparing N0/N1 time-series plot...", flush=True)
        outfile = save_counts_timeseries(ar, args.outdir, args.dpi)
        print(f"Saved N0/N1 time-series plot: {outfile}", flush=True)

    print("All done.", flush=True)


if __name__ == "__main__":
    main()
