"""Shared helpers for the confluent-layer-with-memory analyses.

Paths are derived from this file's own location, so the same scripts run unchanged on
the Mac and on the cluster checkout -- no hard-coded /home/helu anywhere.
"""
import os, sys
import numpy as np

PLOT_PYTHON = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PLOT_PYTHON not in sys.path:
    sys.path.insert(0, PLOT_PYTHON)

from archive.archive import loadarchive          # noqa: E402
import plot.plot_hd as ph                        # noqa: E402
import plot.defects as pd                        # noqa: E402

SCRATCH = os.environ.get("MASS_SCRATCH", "/scratch/helu/mass_hd")


def case_root(study, case):
    return os.path.join(SCRATCH, "cases", study, case)


def results_root(*parts):
    return os.path.join(SCRATCH, "results", "cases", *parts)


def frame_count(root):
    return len([f for f in os.listdir(root)
                if f.startswith("frame") and f.endswith(".json")])


def iter_frames(root):
    """Yield (index, time, frame) for every readable frame, stopping at the first gap."""
    oa = loadarchive(root)
    for i in range(frame_count(root)):
        try:
            fr = oa.read_frame(i)
        except Exception:
            break
        yield i, ph.frame_time(oa, i), fr


def deriv(f, axis):
    """Centred difference on a periodic lattice."""
    return .5 * (np.roll(f, -1, axis) - np.roll(f, 1, axis))


def strain_rate(ux, uy):
    """sqrt(2 E_ij E_ij) with E the symmetric part of grad(u); its mean sets tau_motion."""
    dxux, dyux = deriv(ux, 0), deriv(ux, 1)
    dxuy, dyuy = deriv(uy, 0), deriv(uy, 1)
    exy = .5 * (dxuy + dyux)
    return np.sqrt(2 * (dxux**2 + dyuy**2 + 2 * exy**2))


def n_defects(qxx, qyx):
    npl, nmi = pd.count_defects(pd.director_angle(qxx, qyx))
    return npl + nmi


class Radial:
    """Isotropic radial binning of periodic correlation maps on an L x L lattice."""

    def __init__(self, L):
        self.L = L
        yy, xx = np.meshgrid(np.arange(L), np.arange(L), indexing="ij")
        rr = np.sqrt(np.minimum(xx, L - xx)**2 + np.minimum(yy, L - yy)**2)
        self.bins = np.arange(0, L // 2 + 1)
        self.idx = np.digitize(rr.ravel(), self.bins) - 1
        self.count = np.bincount(self.idx, minlength=len(self.bins))[:len(self.bins)]

    def mean(self, field2d):
        s = np.bincount(self.idx, weights=field2d.ravel(),
                        minlength=len(self.bins))[:len(self.bins)]
        return np.where(self.count > 0, s / np.maximum(self.count, 1), np.nan)

    def corr(self, fs, gs=None):
        """Correlation map of one or several component fields, means removed."""
        acc = np.zeros((self.L, self.L))
        gs = fs if gs is None else gs
        for f, g in zip(fs, gs):
            F = np.fft.fft2(f - f.mean())
            G = np.fft.fft2(g - g.mean())
            acc += np.real(np.fft.ifft2(F * np.conj(G)))
        return acc / (self.L * self.L)

    def length(self, c):
        """Correlation length: first zero crossing, else the 1/e point, else nan."""
        r = self.bins
        for i in range(1, len(c)):
            if np.isfinite(c[i]) and c[i] <= 0:
                a, b = c[i - 1], c[i]
                return float(r[i - 1] + (r[i] - r[i - 1]) * a / (a - b)) if a != b else float(r[i])
        for i in range(1, len(c)):
            if np.isfinite(c[i]) and c[i] <= np.exp(-1):
                return float(r[i])
        return float("nan")


def parse_variant(case):
    """z<zeta>_tc<ratio>_tm<ratio>, or z<zeta>_tcoff for the frozen-chi control."""
    p = case.split("_")
    out = {"zeta": int(p[0][1:]), "tc_ratio": None, "tm_ratio": None}
    if len(p) >= 3 and p[1] != "tcoff":
        out["tc_ratio"] = float(p[1][2:].replace("p", "."))
        out["tm_ratio"] = float(p[2][2:].replace("p", "."))
    return out
