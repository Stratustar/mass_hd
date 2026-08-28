"""Open-loop calibration of a confluent-wet activity: the numbers a closed-loop runcard needs.

Delivers per case
    pmem       = median(P)      the loop threshold.  With chi0 = mc = 0.5 this makes the
                                uniform state an EXACT fixed point, so the loop dynamics comes
                                only from the spatial fluctuations of P.  Getting it wrong is
                                not a small error: the 20260826 batch used pmem = 0 against a
                                true median of +1.3e-2 and the memory was driven by that
                                offset rather than by the physics.
    pmem-width = IQR(P)/2       the dry line's convention (cm_regime), not sigma_P: the tanh
                                should saturate on the tails, not be linear across them.
    t_eddy     = 1/lambda       the material clock, from the measured rms strain rate.
    tau_motion = d/u_rms        the same clock the other way; they agree to ~15% because the
                                velocity structure scale IS the defect spacing.

and reports the saturation check that makes those numbers meaningful, plus L_u and L_P.

Usage:  cw_calib.py <study> [--cases A1 A3 A9] [--out results.json] [--nlast 50] [--stride 3]
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cw_common as cw

CS = 1.0 / np.sqrt(3.0)          # lattice speed of sound


def coarse_series(root, nwant=25):
    """A cheap pass over the WHOLE run, to prove the steady window is actually steady."""
    from archive.archive import loadarchive
    nf = cw.frame_count(root)
    step = max(1, nf // nwant)
    oa = loadarchive(root)
    out = {k: [] for k in ("t", "urms", "nd", "Pmed", "Psig", "melt")}
    for i in range(0, nf, step):
        try:
            fr = cw.load_frame(oa, i)
        except Exception:
            break
        out["t"].append(cw.ph.frame_time(oa, i))
        out["urms"].append(float(np.sqrt(np.mean(fr["ux"]**2 + fr["uy"]**2))))
        out["nd"].append(cw.n_defects(fr["qxx"], fr["qyx"]))
        out["Pmed"].append(float(np.median(fr["P"])))
        out["Psig"].append(float(fr["P"].std()))
        out["melt"].append(float((fr["q2"] < 0.5).mean()))
    return {k: np.asarray(v, float) for k, v in out.items()}


def saturated(series, frac=0.5, tol=0.15):
    """Has u_rms stopped moving?  Compare the last quarter to the preceding quarter.

    Returned as a number, not a verdict: a run that is still drifting must be visible in the
    JSON, because every calibration number downstream assumes saturation.
    """
    u = series["urms"]
    n = len(u)
    if n < 8:
        return {"drift": float("nan"), "ok": False}
    a, b = u[n // 2:3 * n // 4].mean(), u[3 * n // 4:].mean()
    drift = abs(b - a) / max(abs(a), 1e-30)
    return {"drift": float(drift), "ok": bool(drift < tol)}


def analyse(root, nlast, stride):
    par = cw.read_params(root)
    L = int(par.get("LX", 0))
    zeta = float(par.get("zeta", np.nan))
    CC = float(par.get("CC", np.nan))
    chi0 = float(par.get("chi0", 0.5))
    zeta_eff = zeta * (1.0 - chi0)

    ser = coarse_series(root)
    sat = saturated(ser)

    times, frames, pos = cw.steady_frames(root, nlast=nlast, stride=stride)
    rad = cw.Radial(L)

    P = np.concatenate([frames[p]["P"].ravel() for p in pos])
    q1, q3 = np.percentile(P, [25, 75])
    urms = float(np.mean([np.sqrt(np.mean(frames[p]["ux"]**2 + frames[p]["uy"]**2))
                          for p in pos]))
    lam = float(np.mean([cw.strain_rate(frames[p]["ux"], frames[p]["uy"]).mean()
                         for p in pos]))
    nd = float(np.mean([cw.n_defects(frames[p]["qxx"], frames[p]["qyx"]) for p in pos]))
    melt = float(np.mean([(frames[p]["q2"] < 0.5).mean() for p in pos]))
    q2min = float(np.min([frames[p]["q2"].min() for p in pos]))

    Cuu, _ = cw.spatial_corr(rad, frames, pos, ("ux", "uy"))
    Cpp, _ = cw.spatial_corr(rad, frames, pos, "P")
    d = L / np.sqrt(nd) if nd > 0 else float("nan")

    return {
        "case": os.path.basename(root),
        "L": L, "zeta": zeta, "CC": CC, "chi0": chi0,
        "zeta_eff": zeta_eff, "A": zeta_eff / CC if CC else float("nan"),
        "nsteps": float(par.get("nsteps", np.nan)),
        "frames_total": cw.frame_count(root),
        "frames_used": len(pos), "stride": stride,
        "window_t": [float(times[pos[0]]), float(times[pos[-1]])],
        # --- what the runcards need
        "pmem": float(np.median(P)),
        "pmem_width_iqr2": float(0.5 * (q3 - q1)),
        "sigma_P": float(P.std()),
        "t_eddy": 1.0 / lam if lam > 0 else float("nan"),
        "tau_motion": d / urms if urms > 0 else float("nan"),
        # --- context
        "lambda": lam, "u_rms": urms, "Ma": urms / CS,
        "N_def": nd, "d": d,
        "melted_frac": melt, "q2_min": q2min,
        "L_u": cw.length_1e(rad, Cuu), "L_P": cw.length_1e(rad, Cpp),
        "L_u_zero": rad.length(Cuu), "L_P_zero": rad.length(Cpp),
        "xi_N": float(np.sqrt(par.get("LL", np.nan) / (2 * CC))) if CC else float("nan"),
        "saturation": sat,
        "series": {k: v.tolist() for k, v in ser.items()},
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("study")
    ap.add_argument("--cases", nargs="*", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--nlast", type=int, default=50)
    ap.add_argument("--stride", type=int, default=3)
    a = ap.parse_args()

    base = cw.case_root(*a.study.split("/"))
    cases = a.cases or sorted(d for d in os.listdir(base)
                              if os.path.isdir(os.path.join(base, d)))
    res = []
    for c in cases:
        root = os.path.join(base, c)
        try:
            r = analyse(root, a.nlast, a.stride)
        except Exception as exc:
            print(f"  {c}: FAILED -- {type(exc).__name__}: {exc}", flush=True)
            continue
        res.append(r)
        s = "sat" if r["saturation"]["ok"] else f"DRIFT {r['saturation']['drift']:.1%}"
        print(f"  {c}: A={r['A']:.2f}  pmem={r['pmem']:+.4e}  IQR/2={r['pmem_width_iqr2']:.4e}"
              f"  sigma_P={r['sigma_P']:.4e}  t_eddy={r['t_eddy']:.0f}"
              f"  tau_motion={r['tau_motion']:.0f}  N={r['N_def']:.0f} d={r['d']:.1f}"
              f"  L_u={r['L_u']:.1f}({r['L_u_zero']:.0f}) L_P={r['L_P']:.1f}({r['L_P_zero']:.0f})"
              f"  Ma={r['Ma']:.3f}  [{s}]", flush=True)

    out = a.out or cw.results_root(*a.study.split("/"), "cw_calib.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as f:
        json.dump(res, f, indent=1)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
