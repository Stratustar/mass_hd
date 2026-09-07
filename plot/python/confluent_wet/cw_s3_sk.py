#!/usr/bin/env python3
"""Normalised structure factors of chi, m and P over the closing window of a run.

    cw_s3_sk.py --root SIM_ROOT --cases CASE [CASE ...] --out DIR [--window 100] [--tau-c 674.3]

For every case, every full frame inside the closing `window` tau_c is read at full
resolution, the fluctuation field delta f = f - <f> is Fourier transformed, |F_k|^2 is
averaged over the frames, and the result is binned in rings of integer wavenumber
n = |k| L / 2 pi. Two normalisations are written, because they answer different questions:

    S_mode(k) / Var   the per-mode spectrum, i.e. the isotropic power spectrum divided by
                      the total variance -- the SHAPE of the correlations, comparable
                      between fields whose variances differ by orders of magnitude
    S_shell(k) / Var  the power per ring, summing to 1 -- where the variance LIVES in k;
                      its peak is the dominant scale and its mean the inverse correlation
                      length

plus, per field, the 1/e length of the radial correlation function obtained from the same
spectrum (so chi, m and P get the same estimator; cw_part reports L_P and L_chi but not
L_m). A field whose variance is below 1e-8 -- chi in a saturated phase -- is marked
saturated and skipped: there is no structure to factor.

Parseval fixes the normalisation: with F = fft2(delta f) unnormalised, sum_k |F_k|^2 / N^2
= Var, N = L^2, so P2 = |F|^2 / N^2 sums to the variance over all modes.
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cw_common as cw                                       # noqa: E402

FIELDS = ("chi", "m", "P")
VAR_FLOOR = 1e-8


def ring_index(L):
    n = np.fft.fftfreq(L) * L
    nx, ny = np.meshgrid(n, n, indexing="ij")
    return np.rint(np.sqrt(nx ** 2 + ny ** 2)).astype(int)


def radial_corr(P2, L, nr=None):
    """Radially averaged correlation function C(r)/C(0) from the 2-D power spectrum."""
    C = np.real(np.fft.ifft2(P2)) * L * L          # inverse of P2 = |F|^2/N^2 gives <df df>/N... times N
    C = C / C[0, 0]
    x = np.fft.fftfreq(L) * L
    X, Y = np.meshgrid(x, x, indexing="ij")
    r = np.rint(np.sqrt(X ** 2 + Y ** 2)).astype(int)
    nr = nr or L // 2
    prof = np.array([C[r == i].mean() for i in range(nr)])
    return prof


def length_1e(prof):
    below = np.nonzero(prof < np.exp(-1.0))[0]
    if not len(below):
        return float("nan")
    i = below[0]
    if i == 0:
        return 0.0
    # linear interpolation between i-1 and i
    y0, y1 = prof[i - 1], prof[i]
    return float((i - 1) + (y0 - np.exp(-1.0)) / (y0 - y1))


def one_case(root, window, tau_c):
    from archive.archive import loadarchive
    par = cw.read_params(root)
    L = int(par["LX"])
    nsteps = int(par["nsteps"])
    t_start = nsteps - window * tau_c
    oa = loadarchive(root)
    nf = cw.frame_count(root)
    idx = [i for i in range(nf) if cw.ph.frame_time(oa, i) >= t_start]
    ring = ring_index(L)
    nring = L // 2 + 1
    counts = np.bincount(ring.ravel(), minlength=nring)[:nring]
    acc = {f: np.zeros((L, L)) for f in FIELDS}
    var = {f: [] for f in FIELDS}
    mean = {f: [] for f in FIELDS}
    times = []
    for i in idx:
        try:
            fr = cw.load_frame(oa, i)
        except Exception:
            break
        times.append(float(cw.ph.frame_time(oa, i)))
        for f in FIELDS:
            a = np.asarray(fr[f], float)
            d = a - a.mean()
            F = np.fft.fft2(d)
            acc[f] += (np.abs(F) ** 2) / (L * L) ** 2
            var[f].append(float(d.var()))
            mean[f].append(float(a.mean()))
    if len(times) < 2:
        raise RuntimeError(f"{root}: only {len(times)} frames in the closing {window} tau_c")
    n = len(times)
    k = 2 * np.pi * np.arange(nring) / L
    out = {"case": os.path.basename(root.rstrip("/")),
           "L": L, "tau_m": float(par["tau_m"]), "tau_m_over_tau_c": float(par["tau_m"]) / tau_c,
           "chi_config": str(par.get("chi_config", "")), "chi0": float(par.get("chi0", np.nan)),
           "frames_used": n, "times": times, "window_tau_c": window,
           "k": k.tolist(), "ring_counts": counts.tolist(), "fields": {}}
    for f in FIELDS:
        P2 = acc[f] / n
        v = float(np.mean(var[f]))
        if v < VAR_FLOOR:
            out["fields"][f] = {"var": v, "mean": float(np.mean(mean[f])), "saturated": True}
            continue
        shell = np.bincount(ring.ravel(), weights=P2.ravel(), minlength=nring)[:nring]
        mode = shell / np.maximum(counts, 1)
        total = shell.sum()                       # == v up to the frame average
        shell_n = shell / total
        mode_n = mode / total
        prof = radial_corr(P2, L)
        kk = k[1:]
        out["fields"][f] = {
            "var": v, "mean": float(np.mean(mean[f])), "saturated": False,
            "S_mode_norm": mode_n.tolist(), "S_shell_norm": shell_n.tolist(),
            "k_mean": float((kk * shell_n[1:]).sum() / shell_n[1:].sum()),
            "k_peak": float(kk[int(np.argmax(shell_n[1:]))]),
            "corr_profile": prof.tolist(),
            "r_1e": length_1e(prof),
        }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--cases", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--window", type=float, default=100.0)
    ap.add_argument("--tau-c", type=float, default=674.329)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    for c in a.cases:
        try:
            res = one_case(os.path.join(a.root, c), a.window, a.tau_c)
        except Exception as exc:
            print(f"{c}: FAILED {type(exc).__name__}: {exc}", flush=True)
            continue
        with open(os.path.join(a.out, f"sk_{c}.json"), "w") as fh:
            json.dump(res, fh)
        bits = []
        for f, d in res["fields"].items():
            if d.get("saturated"):
                bits.append(f"{f}: saturated")
            else:
                bits.append(f"{f}: var {d['var']:.2e} r1e {d['r_1e']:.2f} kpeak {d['k_peak']:.3f}")
        line = "  ".join(bits)
        print(f"{c}: {res['frames_used']} frames  {line}", flush=True)


if __name__ == "__main__":
    main()
