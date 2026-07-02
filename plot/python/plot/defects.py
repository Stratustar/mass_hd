################################################################################
#
# Topological defect detection for the 2D nematic director field
#
# Vectorised numpy port of the winding-number method in
#   plot/matlab/plot/lyotropic/{wang,calcs,calcnoofdefects}.m
#
# The nematic director is headless: the angle theta is defined mod pi.  The
# charge is the winding of the director around the 2x2 plaquette between four
# neighbouring nodes (the same winding principle as wang.m / calcs.m, but
# plaquette-centred instead of node-centred).  This places each disclination on
# a single plaquette, so the charge field is exact -- no 4-node smearing and no
# /4 correction (cf. calcnoofdefects.m).  The charge q (shape (H-1, W-1)) is the
# loop sum / 2pi:
#     q = +/- 1/2  at half-integer disclinations,
#     q = +/- 1    at integer defects.
# Summing q over a region equals the net winding around its boundary (discrete
# Stokes), so net_charge() with the colony (phi>0.5) mask answers questions like
# "is there a +1 defect enclosed by the go-grow loop?".
#
################################################################################

import numpy as np

TWO_PI = 2.0 * np.pi


def director_angle(qxx, qyx):
    """Nematic director angle theta in (-pi/2, pi/2] from Q-tensor components.

    Matches director_field() in director_analysis.py: theta = 0.5 atan2(Qyx, Qxx).
    """
    return 0.5 * np.arctan2(qyx, qxx)


def _wang(ax, ay, bx, by):
    """Signed nematic angle from director a to director b (vectorised wang.m).

    Works on director *vectors*, not the global angle field, so it is immune to
    the atan2 branch cut and treats the +/- n ambiguity symmetrically: b is
    flipped when the two directors are more than 90 deg apart, then the signed
    angle is atan2(cross, dot) in (-pi/2, pi/2).
    """
    dot = ax * bx + ay * by
    flip = dot < 0.0
    bx = np.where(flip, -bx, bx)
    by = np.where(flip, -by, by)
    dot = ax * bx + ay * by
    cross = ax * by - ay * bx
    return np.arctan2(cross, dot)


def _plaquette_valid(mask):
    """True for plaquettes whose four corner nodes all lie inside mask."""
    m = np.asarray(mask, bool)
    return m[:-1, :-1] & m[1:, :-1] & m[1:, 1:] & m[:-1, 1:]


def defect_charge(theta, mask=None):
    """Plaquette topological charge field q, shape (H-1, W-1).

    q[i, j] is the director winding (/2pi) around the 2x2 plaquette whose corners
    are nodes (i,j), (i+1,j), (i+1,j+1), (i,j+1).  With a mask, plaquettes touching
    vacuum (director undefined) are set to 0.
    """
    th = np.asarray(theta, float)
    nx, ny = np.cos(th), np.sin(th)
    # four plaquette corners, traversed CCW in physical (x=col, y=row) axes:
    # (i,j) -> (i,j+1) -> (i+1,j+1) -> (i+1,j) so that a +1/2 texture gives +1/2.
    ax, ay = nx[:-1, :-1], ny[:-1, :-1]
    bx, by = nx[:-1, 1:], ny[:-1, 1:]
    cx, cy = nx[1:, 1:], ny[1:, 1:]
    dx, dy = nx[1:, :-1], ny[1:, :-1]
    inc = (_wang(ax, ay, bx, by) + _wang(bx, by, cx, cy)
           + _wang(cx, cy, dx, dy) + _wang(dx, dy, ax, ay))
    q = inc / TWO_PI
    if mask is not None:
        q = np.where(_plaquette_valid(mask), q, 0.0)
    return q


def net_charge(theta, mask=None, q=None):
    """Net enclosed topological charge = sum of the plaquette charge field."""
    if q is None:
        q = defect_charge(theta, mask)
    return float(np.sum(q))


def find_defects(theta, mask=None, min_charge=0.3, q=None):
    """Locate defects as connected clusters of |q| > min_charge.

    Returns a list of dicts {row, col, charge}: the charge-weighted centroid of
    each cluster and its summed charge (~ +/-0.5 or +/-1). Needs scipy.ndimage.
    """
    if q is None:
        q = defect_charge(theta, mask)
    try:
        from scipy import ndimage
    except ImportError as e:  # pragma: no cover
        raise ImportError("find_defects requires scipy.ndimage") from e

    defects = []
    for sgn in (+1.0, -1.0):
        sel = (sgn * q) > min_charge
        labels, n = ndimage.label(sel)
        for i in range(1, n + 1):
            cl = labels == i
            rows, cols = np.nonzero(cl)
            w = np.abs(q[cl])
            wsum = w.sum()
            defects.append({
                # +0.5: plaquette charge sits between the four corner nodes
                "row": float((rows * w).sum() / wsum) + 0.5,
                "col": float((cols * w).sum() / wsum) + 0.5,
                "charge": float(q[cl].sum()),
            })
    return defects


def count_defects(theta, mask=None, min_charge=0.3, q=None):
    """Return (n_plus, n_minus): number of +1/2 and -1/2 disclinations."""
    defects = find_defects(theta, mask, min_charge, q)
    n_plus = sum(1 for d in defects if d["charge"] > 0)
    n_minus = sum(1 for d in defects if d["charge"] < 0)
    return n_plus, n_minus


def defect_markers(qxx, qyx, mask=None, min_charge=0.35):
    """Detect defects from Q-tensor grids and return plot-coordinate markers.

    qxx, qyx are the Q-tensor component grids in the plotting index convention
    grid()==[x, y] (axis 0 = x, axis 1 = y), as produced by plot_hd.grid /
    field_as_grid.  mask (same shape) restricts detection to real material
    (e.g. phi > 0.5) so the random exterior director does not spawn spurious
    defects.  Returns (plus, minus): two lists of (x, y) positions in lattice /
    plot data coordinates -- i.e. they overlay directly on
    imshow(field.T, origin="lower") and quiver(x_index, y_index).

    The winding routines here assume the image convention theta[row=y, col=x], so
    this helper transposes internally: theta[x,y] -> theta.T[y,x], detects, then
    maps the returned (row=y, col=x) back to (x, y).  Positive charge (+1/2, and
    the rare +1) goes to `plus`, negative to `minus`.
    """
    theta = director_angle(np.asarray(qxx, float), np.asarray(qyx, float))  # [x, y]
    theta_img = theta.T                                                     # [row=y, col=x]
    mask_img = None if mask is None else np.asarray(mask, bool).T
    ds = find_defects(theta_img, mask=mask_img, min_charge=min_charge)
    plus = [(d["col"], d["row"]) for d in ds if d["charge"] > 0]   # (x, y)
    minus = [(d["col"], d["row"]) for d in ds if d["charge"] < 0]  # (x, y)
    return plus, minus
