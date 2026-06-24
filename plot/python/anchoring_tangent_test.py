"""Diagnose why anchoring A is low: compare interface-tangent estimators.

For a few cases, recompute A = <cos 2(theta_dir - tangent)> on the chi=0.5 contour
with three tangent estimators:
  raw      : np.gradient on the raw (pixel-jagged) contour vertices (current method)
  smooth-w : finite difference over +/- w contour points (smoothed tangent)
  gradchi  : tangent perpendicular to the chi-field gradient (normal = grad chi)

Usage: python anchoring_tangent_test.py <scratch_study_dir> [variant ...]
"""
import sys
import os
import numpy as np

SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)
from archive.archive import loadarchive
from director_analysis import grid, director_field, main_contour


def A_from(theta, main, tangent):
    x, y = main[:, 0], main[:, 1]
    ix = np.clip(np.round(x).astype(int), 0, theta.shape[0] - 1)
    iy = np.clip(np.round(y).astype(int), 0, theta.shape[1] - 1)
    return float(np.mean(np.cos(2 * (theta[ix, iy] - tangent))))


def main():
    study = sys.argv[1]
    variants = sys.argv[2:] or [
        "a0p0002_z0p01_O0p01_LL0p002_xi0p3",   # expected high A
        "a0p0001_z0p05_O0p01_LL0p002_xi0p3",   # mid
        "a0p0001_z0p1_O0p01_LL0p002_xi0p3",    # expected low A (high activity)
        "a0p00005_z0p01_O0p02_LL0p002_xi0p3",  # small colony
    ]
    print("%-42s %6s %7s %7s %7s %8s" % ("variant", "Araw", "Asm3", "Asm7", "Agchi", "vtx_d"))
    for v in variants:
        d = os.path.join(study, v)
        if not os.path.isdir(d):
            print(v, "MISSING"); continue
        ar = loadarchive(d)
        fr = ar.read_frame(int((ar.nsteps - ar.nstart) / ar.ninfo))
        chi, phi = grid(fr, "chi"), grid(fr, "phi")
        theta, _, _, _ = director_field(fr)
        main_c = main_contour(chi, phi)
        if main_c is None:
            print(v, "no-contour"); continue
        x, y = main_c[:, 0], main_c[:, 1]
        n = len(x)
        vtx = float(np.median(np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2)))

        t_raw = np.arctan2(np.gradient(y), np.gradient(x))

        def smooth(w):
            ip = np.roll(np.arange(n), -w); im = np.roll(np.arange(n), w)
            return np.arctan2(y[ip] - y[im], x[ip] - x[im])

        gx, gy = np.gradient(chi)
        ix = np.clip(np.round(x).astype(int), 0, chi.shape[0] - 1)
        iy = np.clip(np.round(y).astype(int), 0, chi.shape[1] - 1)
        t_gchi = np.arctan2(gy[ix, iy], gx[ix, iy]) + np.pi / 2

        print("%-42s %6.2f %7.2f %7.2f %7.2f %8.2f" % (
            v, A_from(theta, main_c, t_raw), A_from(theta, main_c, smooth(3)),
            A_from(theta, main_c, smooth(7)), A_from(theta, main_c, t_gchi), vtx))


if __name__ == "__main__":
    main()
