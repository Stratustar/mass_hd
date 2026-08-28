"""The four correlation analyses of a closed-loop confluent-wet scan, per case.

  1. VELOCITY AUTOCORRELATION      C_uu(r) -> L_u, and C_uu(tau) -> the Eulerian decorrelation
                                   time.  L_u is the yardstick everything else is read against:
                                   in this model the velocity structure scale IS the defect
                                   spacing, so L_u is the size of one "defect encounter".
  2. PRESSURE - MEMORY             C_Pm(r), and the lagged C_Pm(tau).  The loop says
                                   D_t m = (g(P) - m)/tau_m, so m must lag P by ~tau_m.  This
                                   is the first closing link and the one that can be checked
                                   against a number the runcard set.
  3. MEMORY - PHENOTYPE            C_m,chi(r) and lagged C_m,chi(tau).  Same test for the
                                   second link, D_t chi = (chi*(m) - chi)/tau_chi.
  4. PHENOTYPE SPATIO-TEMPORAL     C_chichi(r, tau).  The domain size on the r axis and the
                                   domain lifetime on the tau axis, in one map.  A ridge that
                                   stays at r = 0 means domains sit and fade; a ridge that
                                   walks out means they are advected or propagate.

Every lagged correlation is computed twice, and the pair is the point:
  * Eulerian, point-wise -- also decorrelates on the advection time t_eddy, so it under-reads
    the lag whenever tau_m or tau_chi exceeds t_eddy;
  * k = 0, on the box-averaged time series -- advection cannot decorrelate a box average, so
    this sees the loop's own lag at any tau.  Where the two disagree, the k=0 one is the loop.

Correlation MAPS use the last 50 frames at stride 3; the lag axis uses the contiguous
stride-1 window.  See cw_common.steady_frames.

Usage:
  cw_scan.py <study> [--cases ...] [--out results.json] [--figdir DIR]
             [--nlast 50] [--stride 3] [--maxlag-frac 0.4]
"""
import argparse
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cw_common as cw

VARIANT = re.compile(r"^A(?P<A>[\dp]+)_s(?P<s>[pm])_D(?P<D>[\dp]+)"
                     r"_tm(?P<tm>\d+)_tc(?P<tc>[\dp]+)$")


def parse_variant(name):
    """A<activity>_s<p|m>_D<Dbio>_tm<tau_m>_tc<tau_chi/tau_m>, e.g. A3_sp_D0p03_tm945_tc3."""
    mt = VARIANT.match(name)
    if not mt:
        return {"A": np.nan, "switch_sign": 0, "Dbio": np.nan,
                "tau_m": np.nan, "tc_ratio": np.nan}
    num = lambda s: float(s.replace("p", "."))
    tm = num(mt["tm"])
    return {"A": num(mt["A"]), "switch_sign": 1 if mt["s"] == "p" else -1,
            "Dbio": num(mt["D"]), "tau_m": tm, "tc_ratio": num(mt["tc"]),
            "tau_chi": tm * num(mt["tc"])}


def first_peak(c, lags):
    """Lag of the maximum of |C|, with its sign -- the measured delay between two fields."""
    if len(c) == 0 or not np.any(np.isfinite(c)):
        return float("nan"), float("nan")
    i = int(np.nanargmax(np.abs(c)))
    return float(lags[i]), float(c[i])


def analyse(root, nlast, stride, maxlag_frac):
    par = cw.read_params(root)
    L = int(par["LX"])
    CC, zeta, chi0 = float(par["CC"]), float(par["zeta"]), float(par["chi0"])
    tau_m = float(par["tau-m"]) if "tau-m" in par else float(par.get("tau_m", np.nan))
    tau_chi = float(par["tau-chi"]) if "tau-chi" in par else float(par.get("tau_chi", np.nan))

    times, frames, pos = cw.steady_frames(root, nlast=nlast, stride=stride)
    rad = cw.Radial(L)
    dt = float(np.mean(np.diff(times))) if len(times) > 1 else np.nan
    maxlag = max(4, int(maxlag_frac * len(frames)))
    lags = np.arange(0, min(maxlag, len(frames) - 2) + 1) * dt

    out = {"case": os.path.basename(root), "L": L, "CC": CC, "zeta": zeta,
           "chi0": chi0, "tau_m_set": tau_m, "tau_chi_set": tau_chi,
           "Dbio": float(par.get("Dbio", np.nan)),
           "switch_sign": int(par.get("switch-sign", par.get("switch_sign", 1))),
           "pmem": float(par.get("pmem", np.nan)),
           "pmem_width": float(par.get("pmem-width", par.get("pmem_width", np.nan))),
           "frames_total": cw.frame_count(root), "frames_lag": len(frames),
           "frames_map": len(pos), "stride": stride, "dt": dt,
           "window_t": [float(times[0]), float(times[-1])]}
    out.update(parse_variant(out["case"]))

    # ---- scalar time series (the k=0 mode, and the loop's own state)
    S = {"t": times.tolist()}
    for k, fn in (("chibar", lambda f: f["chi"].mean()),
                  ("chistd", lambda f: f["chi"].std()),
                  ("mbar", lambda f: f["m"].mean()),
                  ("mstd", lambda f: f["m"].std()),
                  ("Pbar", lambda f: f["P"].mean()),
                  ("Pmed", lambda f: np.median(f["P"])),
                  ("Pstd", lambda f: f["P"].std()),
                  ("urms", lambda f: np.sqrt(np.mean(f["ux"]**2 + f["uy"]**2))),
                  ("melt", lambda f: (f["q2"] < 0.5).mean())):
        S[k] = [float(fn(f)) for f in frames]
    S["nd"] = [float(cw.n_defects(f["qxx"], f["qyx"])) for f in frames]
    S["zeta_eff"] = [zeta * (1.0 - c) for c in S["chibar"]]
    out["series"] = S
    out["zeta_eff_mean"] = float(np.mean(S["zeta_eff"]))
    out["A_eff"] = out["zeta_eff_mean"] / CC
    out["chi_bar"] = float(np.mean(S["chibar"]))
    out["u_rms"] = float(np.mean(S["urms"]))
    lam = float(np.mean([cw.strain_rate(f["ux"], f["uy"]).mean() for f in frames]))
    out["lambda"] = lam
    out["t_eddy"] = 1.0 / lam if lam > 0 else float("nan")
    nd = float(np.mean(S["nd"]))
    out["N_def"] = nd
    out["d"] = L / np.sqrt(nd) if nd > 0 else float("nan")

    # ---- 1..3 spatial correlations (maps: last `nlast` frames at stride `stride`)
    sp = {}
    for key, fa, fb in (("uu", ("ux", "uy"), None), ("PP", "P", None),
                        ("mm", "m", None), ("cc", "chi", None), ("QQ", ("qxx", "qyx"), None),
                        ("Pm", "P", "m"), ("mc", "m", "chi"), ("Pc", "P", "chi")):
        c, amp = cw.spatial_corr(rad, frames, pos, fa, fb)
        sp[key] = {"C": c.tolist(), "amp": amp, "length": rad.length(c),
                   "C0": float(c[0]) if len(c) else float("nan")}
    out["spatial"] = sp
    out["r"] = rad.bins.tolist()
    out["L_u"] = sp["uu"]["length"]
    out["L_P"] = sp["PP"]["length"]
    out["L_chi"] = sp["cc"]["length"]
    out["L_m"] = sp["mm"]["length"]

    # ---- lagged correlations, Eulerian and k=0
    tp = {"lags": lags.tolist()}
    for key, fa, fb in (("uu", "ux", "ux"), ("Pm", "P", "m"), ("mc", "m", "chi"),
                        ("Pc", "P", "chi"), ("cc", "chi", "chi")):
        ce = cw.lagged_corr(frames, fa, fb, len(lags) - 1)
        tp[key] = {"euler": ce.tolist()}
        # an autocorrelation peaks at lag 0 by construction, so a "peak lag" is meaningless
        lag, val = (float("nan"), float("nan")) if key in ("uu", "cc") else first_peak(ce, lags)
        tp[key]["euler_peak_lag"] = lag
        tp[key]["euler_peak"] = val
    for key, a, b in (("Pm", "Pbar", "mbar"), ("mc", "mbar", "chibar"),
                      ("Pc", "Pbar", "chibar")):
        ck = cw.k0_lagged_corr(S[a], S[b], len(lags) - 1)
        tp[key]["k0"] = ck.tolist()
        lag, val = first_peak(ck, lags)
        tp[key]["k0_peak_lag"] = lag
        tp[key]["k0_peak"] = val
    # decorrelation times: first 1/e crossing of the autocorrelations
    for key, series in (("uu", tp["uu"]["euler"]), ("cc", tp["cc"]["euler"])):
        c = np.asarray(series, float)
        idx = np.where(c <= np.exp(-1))[0]
        tp[key]["tau_decorr"] = float(lags[idx[0]]) if len(idx) else float("nan")
    out["temporal"] = tp
    out["lag_Pm_k0"] = tp["Pm"]["k0_peak_lag"]
    out["lag_mchi_k0"] = tp["mc"]["k0_peak_lag"]
    out["tau_u_decorr"] = tp["uu"]["tau_decorr"]
    out["tau_chi_decorr"] = tp["cc"]["tau_decorr"]

    # ---- 4. chi spatio-temporal autocorrelation
    st = cw.spatiotemporal_corr(rad, frames, "chi", len(lags) - 1)
    out["chi_st"] = st.tolist()

    # ---- an oscillatory branch would show here
    out["spectra"] = {}
    for k in ("chibar", "mbar", "urms", "Pmed"):
        f, P = cw.detrended_spectrum(times, S[k])
        if len(f):
            out["spectra"][k] = {"period_peak": float(1.0 / f[int(np.argmax(P))]),
                                 "f": f.tolist(), "P": P.tolist()}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("study")
    ap.add_argument("--cases", nargs="*", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--nlast", type=int, default=50)
    ap.add_argument("--stride", type=int, default=3)
    ap.add_argument("--maxlag-frac", type=float, default=0.4)
    a = ap.parse_args()

    base = cw.case_root(*a.study.split("/"))
    cases = a.cases or sorted(d for d in os.listdir(base)
                              if os.path.isdir(os.path.join(base, d)))
    res = []
    for c in cases:
        try:
            r = analyse(os.path.join(base, c), a.nlast, a.stride, a.maxlag_frac)
        except Exception as exc:
            print(f"  {c}: FAILED -- {type(exc).__name__}: {exc}", flush=True)
            continue
        res.append(r)
        print(f"  {c}: A_eff={r['A_eff']:.2f} <chi>={r['chi_bar']:.3f} "
              f"t_eddy={r['t_eddy']:.0f} L_u={r['L_u']:.1f} L_chi={r['L_chi']:.1f} "
              f"lag(P->m)={r['lag_Pm_k0']:.0f}/{r['tau_m_set']:.0f} "
              f"lag(m->chi)={r['lag_mchi_k0']:.0f}/{r['tau_chi_set']:.0f}", flush=True)

    out = a.out or cw.results_root(*a.study.split("/"), "cw_scan.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as f:
        json.dump(res, f)
    print(f"wrote {out}  ({len(res)}/{len(cases)} cases)")


if __name__ == "__main__":
    main()
