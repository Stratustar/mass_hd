import numpy as np
import matplotlib.pyplot as plt

from ..plot_hd import grid

################################################################################
# defaults

VELOCITY_CMAP = "magma"
STREAM_DENSITY = 1.4
STREAM_LINEWIDTH = 0.5
STREAM_ARROWSIZE = 0.5
STREAM_COLOR = "white"


def velocity_field(frame):
    """Return (ux, uy): the MATERIAL velocity, ux_mat/uy_mat.

    Not ux/uy. A confluent-wet frame carries both, and the difference matters. The flow
    comes from a forced lattice Boltzmann, which carries the force trapezoidally in time,
    so the bare moment of the distribution is not the fluid velocity:

        u_phys = (Sum_v f_v c_v + F/2)/n = u_code + F/(2n)

    ux/uy hold the bare moment u_code; ux_mat/uy_mat hold the corrected u_phys, which is
    what advects Q, chi and m inside the model. The gap is largest exactly where the
    interesting physics is -- measured at a melted defect core it reaches 53% of u_rms
    (~4% in developed turbulence), because div(sigma) there is dominated by the Frank
    stress -- so plotting ux/uy would draw a fictitious drift around every core.
    """
    return grid(frame, "ux_mat"), grid(frame, "uy_mat")


def velocity_rms(frame):
    """Root-mean-square material speed over the lattice."""
    ux, uy = velocity_field(frame)
    return float(np.sqrt(np.mean(ux*ux + uy*uy)))


def velocity(frame, engine=plt, cmap=VELOCITY_CMAP, streamlines=True,
             density=STREAM_DENSITY, color=STREAM_COLOR, vmin=None, vmax=None,
             colorbar=True):
    """Plot the speed |u| as a colour map with the flow topology drawn on top as streamlines.

    Streamlines rather than a quiver: in active turbulence the speed spans more than an
    order of magnitude across the box, so arrows scaled by magnitude are dominated by the
    few strongest jets and the vortices are invisible. Separating the two -- topology in the
    white streamlines, magnitude in the background colour -- keeps both readable.

    Note the transposes. Fields are stored flat with k = y + LY*x, so grid() returns
    arr[x][y]; imshow and streamplot both take the FIRST axis as rows, i.e. as y. Passing
    the arrays untransposed silently plots the transpose of the flow, which for an isotropic
    turbulent field looks perfectly plausible.

    vmin/vmax pin the colour scale, e.g. to compare frames of an animation or two cases
    side by side; by default each frame is auto-scaled.
    """
    ux, uy = velocity_field(frame)
    mag = np.hypot(ux, uy)

    cax = engine.imshow(mag.T, interpolation="nearest", origin="lower",
                        cmap=cmap, vmin=vmin, vmax=vmax)

    if streamlines and mag.max() > 0:
        lx, ly = mag.shape
        engine.streamplot(np.arange(lx), np.arange(ly), ux.T, uy.T,
                          density=density, linewidth=STREAM_LINEWIDTH,
                          arrowsize=STREAM_ARROWSIZE, color=color)

    if colorbar:
        cbar = plt.colorbar(cax)
        cbar.set_label(r"$|u|$")
    return cax
