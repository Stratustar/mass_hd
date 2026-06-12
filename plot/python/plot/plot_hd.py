################################################################################
#
# Plotting routines for the go-or-grow / hydrodynamic (hd) models
#
# Follows the same convention as plot/lyotropic.py: every routine takes a
# `frame` (a single output frame) and an `engine` (defaults to pyplot) and
# draws directly onto it. Time-series routines instead take the output
# archive `oa`. Nothing here loads archives or writes files; that stays with
# the caller and with the CLI scripts plot/python/plot_hd.py and
# proliferation_summary.py.
#
################################################################################

import numpy as np
import matplotlib.pyplot as plt

################################################################################
# useful functions

DIRECTOR_STRIDE = 5
DIRECTOR_WIDTH = 0.0015
VISUALIZATION_CMAP = plt.get_cmap("viridis")
VISUALIZATION_BACKGROUND = "#7e7e7e"
BOUNDARY_COLOR = "#d9468f"


def grid(frame, name):
    """Reshape a flat frame field to the (LX, LY) lattice."""
    lx = frame.parameters["LX"]
    ly = frame.parameters["LY"]
    return np.array(getattr(frame, name)).reshape((lx, ly))


def director_components(frame, stride=DIRECTOR_STRIDE):
    """Return (xx, yy, uu, vv, strength) for the nematic director quiver.

    theta = 1/2 atan2(QQyx, QQxx); the glyph length is scaled by the local
    nematic order strength = sqrt(QQxx^2 + QQyx^2).
    """
    lx = frame.parameters["LX"]
    ly = frame.parameters["LY"]
    qxx = grid(frame, "QQxx")
    qyx = grid(frame, "QQyx")
    theta = 0.5 * np.arctan2(qyx, qxx)
    strength = np.sqrt(qxx * qxx + qyx * qyx)

    xs = np.arange(0, lx, stride)
    ys = np.arange(0, ly, stride)
    xx, yy = np.meshgrid(xs, ys, indexing="ij")
    order = np.clip(strength[xx, yy], 0.0, 1.0)
    uu = np.cos(theta[xx, yy]) * order
    vv = np.sin(theta[xx, yy]) * order
    return xx, yy, uu, vv, strength


def visualization_rgba(frame):
    """RGBA composite: hue from phenotype (viridis on 1-chi), alpha from phi."""
    chi = np.clip(grid(frame, "chi"), 0.0, 1.0)
    phi = np.clip(grid(frame, "phi"), 0.0, 1.0)
    rgb = VISUALIZATION_CMAP(1.0 - chi)[..., :3]
    return np.dstack((rgb, phi))


def counts(frame):
    """Return (N0, N1): go-type Sum[phi(1-chi)] and grow-type Sum[phi*chi]."""
    chi = grid(frame, "chi")
    phi = grid(frame, "phi")
    n0 = float(np.sum(phi * (1.0 - chi)))
    n1 = float(np.sum(phi * chi))
    return n0, n1


def frame_count(oa):
    """Number of readable frames in the archive (matches the CLI scripts)."""
    return int((oa.nsteps - oa.nstart) / oa.ninfo) + 1


def frame_time(oa, frame_index):
    """Simulation time of a frame index."""
    step = oa.nstart + frame_index * oa.ninfo
    dt = float(getattr(oa, "time_step", 1.0))
    return step * dt


################################################################################
# actual plotting routines (single frame)

def director(frame, engine=plt, background="order", stride=DIRECTOR_STRIDE):
    """Plot the nematic director as black headless glyphs over a background.

    background = "order" colours by nematic order strength (viridis);
    background = "phi"   colours by the density field phi (summer).
    """
    xx, yy, uu, vv, strength = director_components(frame, stride=stride)
    if background == "phi":
        field, cmap, label = grid(frame, "phi"), "summer", "phi"
    else:
        field, cmap, label = strength, "viridis", "nematic order"

    cax = engine.imshow(field.T, interpolation="nearest", origin="lower", cmap=cmap)
    engine.quiver(
        xx, yy, uu, vv,
        pivot="middle",
        headwidth=0, headlength=0, headaxislength=0,
        angles="xy", scale_units="xy", scale=1.0 / stride,
        width=DIRECTOR_WIDTH, color="black",
    )
    cbar = plt.colorbar(cax)
    cbar.set_label(label)
    return cax


def pressure(frame, engine=plt):
    """Plot the mechanical pressure field p = -sigma_bulk (coolwarm).

    Whether the phase-field surface stress is included is set by the model's
    surface_stress option, not by plotting: surface_stress=0 yields a
    surface-free mechanical pressure, surface_stress=1 includes it.
    """
    field = grid(frame, "pressure") if hasattr(frame, "pressure") else -grid(frame, "sigma_bulk")
    cax = engine.imshow(field.T, interpolation="nearest", origin="lower", cmap="coolwarm")
    cbar = plt.colorbar(cax)
    cbar.set_label("pressure")
    return cax


def visualization(frame, engine=plt, annotate=True):
    """Plot the chi/phi composite with the phi=0.5 colony boundary.

    Hue encodes phenotype (go=bright, grow=dark via viridis on 1-chi) and the
    opacity encodes density phi over a neutral grey background. With annotate
    the go/grow masses N0/N1 are printed next to the panel.
    """
    rgba = visualization_rgba(frame)
    phi = grid(frame, "phi")

    ax = engine.gca() if engine is plt else engine
    ax.set_facecolor(VISUALIZATION_BACKGROUND)
    cax = ax.imshow(
        np.transpose(rgba, (1, 0, 2)), origin="lower", interpolation="nearest"
    )
    ax.contour(phi.T, levels=[0.5], colors=[BOUNDARY_COLOR], linewidths=0.8, origin="lower")
    ax.set_aspect("equal")
    if annotate:
        n0, n1 = counts(frame)
        ax.text(
            1.03, 0.5, f"N0= {n0:.6g}\nN1= {n1:.6g}",
            transform=ax.transAxes, ha="left", va="center", fontsize=11, clip_on=False,
        )
    return cax


################################################################################
# time-series routines (whole archive)

def counts_vs_t(oa, engine=plt, use_time=True):
    """Plot the two biomass curves N0 (go) and N1 (grow) versus time.

    N0 = Sum[phi(1-chi)], N1 = Sum[phi*chi], summed over the lattice each frame.
    use_time=True uses simulation time on the x axis, otherwise the frame index.
    """
    n = frame_count(oa)
    xs, n0s, n1s = [], [], []
    for i in range(n):
        try:
            frame = oa.read_frame(i)
        except Exception:
            break
        n0, n1 = counts(frame)
        xs.append(frame_time(oa, i) if use_time else i)
        n0s.append(n0)
        n1s.append(n1)

    engine.plot(xs, n0s, color="tab:red", label="N0 (go)")
    engine.plot(xs, n1s, color="tab:blue", label="N1 (grow)")
    engine.legend()
    return xs, n0s, n1s


################################################################################
# animation helper

def save_gif(
    oa,
    plotfn,
    outfile,
    frame_start=0,
    frame_step=1,
    frame_stop=None,
    figsize=(6, 5),
    dpi=150,
    duration=200,
    title=None,
):
    """Render a single-frame routine across frames into an animated GIF.

    plotfn(frame) must draw onto the current pyplot figure (the engine=plt
    convention used by director/pressure/visualization). Use functools.partial
    to fix extra arguments, e.g. partial(director, background="phi").
    title, if given, is a callable frame_index -> str for the panel title.
    Requires Pillow and a non-interactive backend (matplotlib.use("Agg")).
    """
    from PIL import Image

    n = frame_count(oa)
    if frame_stop is not None:
        n = min(n, frame_stop)
    indices = list(range(frame_start, n, max(1, frame_step)))
    if not indices:
        raise RuntimeError("No frames selected for GIF.")

    images = []
    for i in indices:
        frame = oa.read_frame(i)
        fig = plt.figure(figsize=figsize)
        plotfn(frame)
        if title is not None:
            plt.title(title(i))
        fig.tight_layout()
        fig.canvas.draw()
        w, h = fig.canvas.get_width_height()
        buf = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8).reshape(h, w, 4)
        images.append(Image.fromarray(buf[..., :3].copy()))
        plt.close(fig)

    images[0].save(
        outfile, save_all=True, append_images=images[1:], duration=duration, loop=0
    )
    return outfile
