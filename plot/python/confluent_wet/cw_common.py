"""Shared helpers for the confluent-wet analyses.

Three things differ from the dry confluent-memory line and each of them is a correctness
issue, not a naming detail:

  * THE VELOCITY IS ux_mat/uy_mat, NOT ux/uy.  A forced lattice Boltzmann carries the force
    trapezoidally, so the bare LB moment is not the material velocity; the two differ by 53%
    of u_rms at a melted defect core (measured 20260826).  Anything that advects, or any
    velocity gradient, must use the _mat pair.
  * THE PRESSURE IS ITS OWN FIELD `pressure`.  The dry line reads -sigma_bulk; here the flow
    is incompressible and P = -1/2 Tr(Pi) is written directly.  `pressure_lb` is the
    gauge-dependent diagnostic and must NOT be used, nor reconstructed from `n` (the JSON
    writer emits ~6 significant digits and n ~ rho = 40 quantises delta_n badly).
  * THERE IS NO phi.  Nothing here may restrict a statistic to the interior of a colony;
    the layer is confluent by construction and every node is material.

The generic lattice math (radial binning, strain rate, defect counting) is imported from the
dry line rather than duplicated -- it is model-agnostic, and one implementation is one place
for a bug to live.
"""
import os
import sys

import numpy as np

_PLOT_PYTHON = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for _p in (_PLOT_PYTHON, os.path.join(_PLOT_PYTHON, "confluent_memory")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from cm_common import Radial, strain_rate, n_defects, frame_count, ph   # noqa: E402
from archive.archive import loadarchive                                  # noqa: E402

SCRATCH = os.environ.get("MASS_SCRATCH", "/scratch/helu/mass_hd")

# name in the analysis -> name in the frame
FIELDS = {"ux": "ux_mat", "uy": "uy_mat", "P": "pressure",
          "chi": "chi", "m": "m", "qxx": "QQxx", "qyx": "QQyx"}


def case_root(*parts):
    return os.path.join(SCRATCH, "cases", *parts)


def results_root(*parts):
    return os.path.join(SCRATCH, "results", "cases", *parts)


def read_params(root):
    """The run's own parameters.json, so the analysis never re-derives what the run recorded."""
    import json
    with open(os.path.join(root, "parameters.json")) as f:
        return json.load(f)


def load_frame(oa, i):
    """One frame as a dict of 2-D arrays, under the analysis names of FIELDS.

    Fails loudly on a frame that lacks ux_mat: that is a pre-confluent-wet archive, and
    silently falling back to `ux` would reintroduce the 53% advection error.
    """
    fr = oa.read_frame(i)
    out = {}
    for name, key in FIELDS.items():
        try:
            out[name] = ph.grid(fr, key)
        except Exception as exc:
            raise RuntimeError(f"frame {i} has no field '{key}' -- not a confluent-wet "
                               f"archive, or an incomplete write ({exc})")
    out["q2"] = out["qxx"]**2 + out["qyx"]**2
    return out


def window_indices(nframes, nlast, stride):
    """The last `nlast` frames sampled with `stride`, oldest first.

    This is the agreed convention for the correlation MAPS: independent snapshots matter more
    than a fine time axis there.  The temporal correlations use a contiguous window instead
    (see `contiguous_indices`), because a stride throws away exactly the lag resolution they
    are trying to measure.
    """
    idx = list(range(nframes - 1, -1, -stride))[:nlast][::-1]
    return idx


def contiguous_indices(nframes, nlast_span):
    """The last `nlast_span` frames at stride 1 -- the lag axis for temporal correlations."""
    return list(range(max(0, nframes - nlast_span), nframes))


def steady_frames(root, nlast=50, stride=3, span=None):
    """Load the analysis window once and serve both conventions from it.

    Returns (times, frames, map_positions) where `frames` is the contiguous stride-1 list
    used for the lag axis and `map_positions` indexes into it the stride-`stride` subset of
    `nlast` frames used for the correlation maps.  Loading is the expensive step (frames are
    tens of MB of JSON), so it happens exactly once.
    """
    nf = frame_count(root)
    if nf < 8:
        raise RuntimeError(f"{root}: only {nf} frames, refusing to call this a steady window")
    span = span if span is not None else min(nf, max(nlast * stride, 60))
    idx = contiguous_indices(nf, span)
    oa = loadarchive(root)
    times, frames = [], []
    for i in idx:
        try:
            frames.append(load_frame(oa, i))
        except Exception:
            break                      # stop at the first unreadable frame, keep what we have
        times.append(ph.frame_time(oa, i))
    if len(frames) < 8:
        raise RuntimeError(f"{root}: only {len(frames)} readable frames in the window")
    # the map subset, taken from the END of what actually loaded
    pos = list(range(len(frames) - 1, -1, -stride))[:nlast][::-1]
    return np.asarray(times, float), frames, pos


# --------------------------------------------------------------------------- correlations

def spatial_corr(rad, frames, positions, fa, fb=None):
    """Isotropic C_ab(r), averaged over the snapshot subset and normalised to C(0) = 1.

    Each snapshot is normalised BEFORE averaging, so a frame that happens to have a larger
    variance does not dominate the shape; the amplitude is reported separately.
    """
    acc = np.zeros(len(rad.bins))
    amp = []
    for p in positions:
        fr = frames[p]
        a = [fr[fa]] if isinstance(fa, str) else [fr[k] for k in fa]
        b = a if fb is None else ([fr[fb]] if isinstance(fb, str) else [fr[k] for k in fb])
        c = rad.mean(rad.corr(a, b))
        if not np.isfinite(c[0]) or c[0] == 0:
            continue
        amp.append(float(c[0]))
        acc += c / abs(c[0])
    return acc / max(len(amp), 1), float(np.mean(amp)) if amp else float("nan")


def lagged_corr(frames, fa, fb, maxlag):
    """Eulerian point-wise lagged cross-correlation C_ab(tau) = <a(x,t) b(x,t+tau)>.

    Normalised by sqrt(var_a var_b) so it is a correlation coefficient.  A peak at tau > 0
    means b LAGS a, which is the whole point: the loop predicts m to lag P by tau_m and chi
    to lag m by tau_chi.

    CAVEAT, and it is why `k0_lagged_corr` exists too: this is measured at a fixed lattice
    point, so it also decorrelates on the advection time t_eddy.  When tau_m >> t_eddy the
    advection wins and the lag is not visible here.
    """
    n = len(frames)
    A = np.stack([frames[i][fa] for i in range(n)])
    B = np.stack([frames[i][fb] for i in range(n)])
    A = A - A.mean(axis=(1, 2), keepdims=True)
    B = B - B.mean(axis=(1, 2), keepdims=True)
    out = []
    for lag in range(0, min(maxlag, n - 2) + 1):
        a, b = A[:n - lag], B[lag:]
        num = float(np.mean(a * b))
        den = float(np.sqrt(np.mean(a**2) * np.mean(b**2)))
        out.append(num / den if den > 0 else np.nan)
    return np.asarray(out)


def k0_lagged_corr(series_a, series_b, maxlag):
    """Lagged correlation of two domain-averaged time series (the k = 0 mode).

    Advection cannot decorrelate a box average, so this measures the LOOP's lag cleanly even
    when tau_m >> t_eddy, which is exactly where the Eulerian point-wise version fails.
    """
    a = np.asarray(series_a, float); b = np.asarray(series_b, float)
    a = a - a.mean(); b = b - b.mean()
    n = len(a)
    out = []
    for lag in range(0, min(maxlag, n - 2) + 1):
        x, y = a[:n - lag], b[lag:]
        den = float(np.sqrt(np.mean(x**2) * np.mean(y**2)))
        out.append(float(np.mean(x * y)) / den if den > 0 else np.nan)
    return np.asarray(out)


def spatiotemporal_corr(rad, frames, field, maxlag):
    """C(r, tau) for one field: the domain size on the r axis, its lifetime on the tau axis.

    A ridge that stays at r = 0 means domains sit still and simply fade; a ridge that walks
    out in r means they are advected or propagate.
    """
    n = len(frames)
    F = np.stack([frames[i][field] for i in range(n)])
    F = F - F.mean(axis=(1, 2), keepdims=True)
    L = F.shape[-1]
    out = []
    norm = None
    for lag in range(0, min(maxlag, n - 2) + 1):
        acc = np.zeros((L, L))
        a, b = F[:n - lag], F[lag:]
        for i in range(len(a)):
            A = np.fft.fft2(a[i]); B = np.fft.fft2(b[i])
            acc += np.real(np.fft.ifft2(A * np.conj(B)))
        c = rad.mean(acc / (len(a) * L * L))
        if norm is None:
            norm = abs(c[0]) if np.isfinite(c[0]) and c[0] != 0 else 1.0
        out.append(c / norm)
    return np.asarray(out)


def detrended_spectrum(t, y):
    """Power spectrum of a time series with the linear trend removed.

    A limit cycle in the loop shows up as a peak well inside the resolvable band; a slow
    drift would otherwise pile all the power into the lowest bins and hide it.
    """
    y = np.asarray(y, float)
    if len(y) < 8:
        return np.array([]), np.array([])
    p = np.polyfit(t, y, 1)
    r = y - np.polyval(p, t)
    dt = float(np.mean(np.diff(t)))
    f = np.fft.rfftfreq(len(r), d=dt)
    P = np.abs(np.fft.rfft(r * np.hanning(len(r))))**2
    return f[1:], P[1:]
