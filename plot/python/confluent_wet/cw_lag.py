#!/usr/bin/env python3
"""Wave-1 calibration of the Lagrangian-memory campaign, from ONE open-loop run.

    cw_lag.py <case_output_dir> [--out part_lag.json] [--figs DIR]

It produces the four numbers a closed-loop runcard cannot guess, plus the checks that say
whether they can be trusted:

    sigma_P              -> pmem-width = 0.5 sigma_P
    <g(P)>               -> mc            (NOT 0.5: <g> is what makes chi* = 1/2 a fixed
                                           point, and g is a tanh against a skewed P)
    std(g(P))            -> chi-width = 0.5 std(g)
    tau_L                -> the unit of the whole tau_m ladder

WHY tau_L AND NOT 1/lambda. The memory obeys D_t m = (g(P) - m)/tau_m along a material
path, so the clock tau_m must be compared against is the decorrelation time of P FOLLOWING
A CELL. 1/lambda, which cw_calib names t_eddy, is an rms strain rate -- an instantaneous
local quantity, not a correlation time at all, and the campaign has been quoting rungs in
it. The ratio tau_L/(1/lambda) is reported so every historical rung can be converted.

The Lagrangian correlation subtracts the INSTANTANEOUS ensemble mean of P over tracers at
each sample, not a global constant: a slow drift in <P> would otherwise appear as a long
correlation tail and inflate tau_L, which is the one number this run exists to measure.
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from confluent_wet import cw_common as cw   # noqa: E402

WIDTH_COEFFS = (0.20, 0.34, 0.50, 0.75, 1.00)


def read_tracers(root):
    """(dt, arr) with arr[sample, tracer, {P,m,x,y}]; the count comes from the CSV header."""
    meta = os.path.join(root, "tracer_meta.csv")
    dat = os.path.join(root, "tracer.f32")
    if not (os.path.exists(meta) and os.path.exists(dat)):
        return None, None
    npart, stride = None, None
    for line in open(meta):
        if "tracers x 4 values" in line:
            npart = int(line.split(":")[1].split("tracers")[0].split(",")[-1])
        if "sampled every" in line:
            stride = int(line.split("sampled every")[1].split("steps")[0])
        if not line.startswith("#"):
            break
    a = np.fromfile(dat, dtype=np.float32)
    a = a[:(a.size // (4 * npart)) * 4 * npart].reshape(-1, npart, 4)
    return stride, a


def lagrangian_corr(P, dt, maxlag_frac=0.25):
    """C(tau) = <dP_i(t) dP_i(t+tau)>_{i,t} / <dP^2>, tracer- and time-averaged."""
    P = P.astype(np.float64)
    P = P - P.mean(axis=1, keepdims=True)      # instantaneous ensemble mean, see docstring
    n = P.shape[0]
    nl = max(2, int(maxlag_frac * n))
    v0 = float((P * P).mean())
    C = np.array([float((P[:n - l] * P[l:]).mean()) / v0 for l in range(nl + 1)])
    return np.arange(nl + 1) * dt, C


def cross_1e(lag, C):
    k = int(np.argmax(C < 1.0 / np.e))
    if k == 0:
        return float("nan")
    x0, x1, y0, y1 = lag[k - 1], lag[k], C[k - 1], C[k]
    return float(x0 + (y0 - 1 / np.e) * (x1 - x0) / (y0 - y1))


def integral_time(lag, C):
    """Integral of C up to its first zero crossing -- reported beside the 1/e time because
    the two agree only for an exponential, and disagreeing is itself information."""
    k = int(np.argmax(C < 0.0))
    if k == 0:
        k = len(C)
    return float(np.trapezoid(C[:k], lag[:k])) if k > 1 else float("nan")


def g_moments(P, pmem, width):
    g = 0.5 * (1.0 + np.tanh((P - pmem) / width))
    return float(g.mean()), float(g.std())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--out", default=None)
    ap.add_argument("--figs", default=None)
    ap.add_argument("--nlast", type=int, default=30)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--pmem", type=float, default=0.0)
    a = ap.parse_args()

    par = cw.read_params(a.root)
    L = int(par.get("LX", 0))
    zeta_eff = float(par.get("zeta_open", par.get("zeta", np.nan)))
    tau_m_set = float(par.get("tau_m", np.nan))

    # ---- field statistics, from the steady frames ---------------------------
    times, frames, pos = cw.steady_frames(a.root, nlast=a.nlast, stride=a.stride)
    rad = cw.Radial(L)
    P = np.concatenate([frames[p]["P"].ravel() for p in pos])
    sigma_P = float(P.std())
    lam = float(np.mean([cw.strain_rate(frames[p]["ux"], frames[p]["uy"]).mean()
                         for p in pos]))
    urms = float(np.mean([np.sqrt(np.mean(frames[p]["ux"]**2 + frames[p]["uy"]**2))
                          for p in pos]))
    nd = float(np.mean([cw.n_defects(frames[p]["qxx"], frames[p]["qyx"]) for p in pos]))
    melt = float(np.mean([(frames[p]["q2"] < 0.5).mean() for p in pos]))
    Cpp, _ = cw.spatial_corr(rad, frames, pos, "P")
    Cuu, _ = cw.spatial_corr(rad, frames, pos, ("ux", "uy"))
    L_P, L_u = cw.length_1e(rad, Cpp), cw.length_1e(rad, Cuu)
    std_m_field = float(np.mean([frames[p]["m"].std() for p in pos]))

    # ---- the g(P) moments that set mc and chi-width -------------------------
    gtab = {}
    for c in WIDTH_COEFFS:
        mg, sg = g_moments(P, a.pmem, c * sigma_P)
        gtab[f"{c:.2f}"] = {"pmem_width": c * sigma_P, "g_mean": mg, "g_std": sg}
    adopted = gtab["0.50"]

    # ---- the Lagrangian clock ----------------------------------------------
    dt, tr = read_tracers(a.root)
    lagr = {}
    if tr is not None:
        i0 = tr.shape[0] // 4                      # discard the spin-up quarter
        lag, C = lagrangian_corr(tr[i0:, :, 0], dt)
        tau_L = cross_1e(lag, C)
        # MSD, a free check that the advection is doing what it should
        L_box = float(L)
        xy = tr[i0:, :, 2:4].astype(np.float64)
        d = np.diff(xy, axis=0)
        d -= L_box * np.round(d / L_box)
        disp = np.cumsum(d, axis=0)
        msd = (disp**2).sum(axis=2).mean(axis=1)
        lagr = {
            "dt": dt, "n_tracer": int(tr.shape[1]), "n_sample": int(tr.shape[0]),
            "tau_L": tau_L, "tau_L_integral": integral_time(lag, C),
            "tau_L_over_inv_lambda": tau_L * lam,
            "msd_at_tau_L": float(np.interp(tau_L, np.arange(1, len(msd) + 1) * dt, msd)),
            "std_P_tracer": float(tr[i0:, :, 0].std()),
            "std_m_tracer": float(tr[i0:, :, 1].std()),
            "C_lag": lag.tolist(), "C": C.tolist(),
        }

    out = {
        "case": os.path.basename(os.path.normpath(a.root)),
        "L": L, "zeta_eff": zeta_eff, "open_loop": int(par.get("open_loop", 0)),
        "tau_m_set": tau_m_set,
        "field": {
            "sigma_P": sigma_P, "P_mean": float(P.mean()),
            "P_median": float(np.median(P)), "P_median_over_sigma": float(np.median(P)/sigma_P),
            "lambda": lam, "inv_lambda": 1.0 / lam if lam > 0 else float("nan"),
            "inv_lambda_times_zeta_eff": zeta_eff / lam if lam > 0 else float("nan"),
            "u_rms": urms, "Ma": urms / (1.0 / np.sqrt(3.0)),
            "N_def": nd, "melted_frac": melt, "L_P": L_P, "L_u": L_u,
            "sigma_P_over_zeta_eff": sigma_P / zeta_eff if zeta_eff else float("nan"),
            "std_m_field": std_m_field,
        },
        "g_scan": gtab,
        "adopted": {
            "pmem": a.pmem,
            "pmem_width": adopted["pmem_width"],
            "mc": adopted["g_mean"],
            "chi_width": 0.5 * adopted["g_std"],
        },
        "lagrangian": lagr,
    }

    # std(m)/std(g) is the DESIGN FORMULA under test: an exponential filter of a process
    # with correlation time tau_L gives sqrt(tau_L/(tau_L+tau_m)). The whole tau_m ladder's
    # predicted chi amplitudes rest on it, and this is the first chance to check it.
    if lagr and np.isfinite(lagr["tau_L"]):
        pred = np.sqrt(lagr["tau_L"] / (lagr["tau_L"] + tau_m_set))
        out["design_check"] = {
            "std_m_over_std_g_measured": std_m_field / adopted["g_std"],
            "std_m_over_std_g_predicted": float(pred),
        }

    txt = json.dumps(out, indent=1)
    if a.out:
        open(a.out, "w").write(txt)
    f, g = out["field"], out["adopted"]
    print(f"{out['case']}  zeta_eff={out['zeta_eff']:g}  L={L}")
    print(f"  sigma_P      = {f['sigma_P']:.5f}   (/zeta_eff = {f['sigma_P_over_zeta_eff']:.3f})")
    print(f"  <P>/sigma_P  = {f['P_mean']/f['sigma_P']:+.4f}   median/sigma_P = "
          f"{f['P_median_over_sigma']:+.4f}")
    print(f"  1/lambda     = {f['inv_lambda']:.1f} steps   (x zeta_eff = "
          f"{f['inv_lambda_times_zeta_eff']:.2f}, the A1 law says 80)")
    print(f"  L_P = {f['L_P']:.2f}   L_u = {f['L_u']:.2f}   N_def = {f['N_def']:.0f}   "
          f"melt = {f['melted_frac']*100:.1f}%   Ma = {f['Ma']:.4f}")
    if lagr:
        print(f"  tau_L (1/e)  = {lagr['tau_L']:.1f} steps   integral = "
              f"{lagr['tau_L_integral']:.1f}")
        print(f"  tau_L / (1/lambda) = {lagr['tau_L_over_inv_lambda']:.4f}   <-- the "
              f"conversion every historical rung needs")
    print(f"  ADOPTED: pmem_width = {g['pmem_width']:.5f}   mc = {g['mc']:.4f}   "
          f"chi_width = {g['chi_width']:.4f}")
    if "design_check" in out:
        d = out["design_check"]
        print(f"  design check std(m)/std(g): measured {d['std_m_over_std_g_measured']:.3f} "
              f"vs predicted {d['std_m_over_std_g_predicted']:.3f}")


if __name__ == "__main__":
    main()
