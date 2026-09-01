#!/usr/bin/env python3
"""Stage A3 of the STEP-SWITCH campaign: read A1 and A2b, decide every closed-loop number.

    cw_step_pick.py <a1_results_dir> <a2b_results_dir> [--out calib.json] [--figs DIR]

WHAT IS DIFFERENT FROM cw_ladder.py. That script serves the tanh campaign, where chi*(m) is
a smoothed step of width w_chi and the design question is how to pick w_chi. Here w_chi is
identically zero, and the smoothing that decides everything comes from sigma_m instead. So:

  * The fixed point is  chibar = Phi( (mc - f(a(chibar))) / sigma_m ),  the space-average of
    a POINTWISE step over a normally distributed m -- verified against a closed-loop run to
    0.03 in chibar. It is not the mean-field chi*(<m>), which would be a step and could
    never sit between 0 and 1.

  * Bistability is a CUSP condition. Differentiating the fixed point at fixed mc,

        sigma_m = psi(z(chibar)) * (1 - r) * |df/da|  ==  K(chibar) ,   z = Phi^-1(chibar)

    so three roots exist iff sigma_m < K_max = max_chibar K(chibar). The tangency defines
    both the optimal threshold and the crossing:

        mc  = f(chibar*) + K_max z(chibar*)          chibar* = argmax K
        tau_x  such that  sigma_m(tau_x) = K_max

  * sigma_m(tau_m) is MEASURED, from A2b, not predicted. The one-pole closed form
    sqrt(f(1-f)) sqrt(tau_c/(tau_c+tau)) is exact at Dbio = 0 (0.953 in a local test) and
    wrong by 3-4x at the campaign's Dbio, because cell motility smooths m over
    sqrt(2 Dbio tau_m). A2b measures the ratio at four tau_m and two activities; the f
    dependence is divided out through sqrt(f(1-f)), leaving a pure h(tau_m) to interpolate.

  * pmem is picked by CONTRAST, f(zeta) - f(0.3 zeta), which is what the campaign asked for
    and what sets the width of the reachable window. With a hard g there is no pmem-width to
    scan beside it, so the scan is one-dimensional.
"""
import argparse
import glob
import json
import math
import os

import numpy as np
from scipy.interpolate import PchipInterpolator

PMEM_COEFFS = (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
R_FLOOR = 0.3


def Phi(z):
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def psi(z):
    return math.exp(-0.5 * z * z) / math.sqrt(2.0 * math.pi)


def zq(p, lo=-8.0, hi=8.0):
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if Phi(mid) < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def load_parts(d):
    out = []
    for p in sorted(glob.glob(os.path.join(d, "*", "part.json"))):
        with open(p) as fh:
            out.append(json.load(fh))
    if not out:
        raise RuntimeError(f"no part.json under {d}")
    return out


def f_hard(part, pmem):
    """f = fraction of TIME-and-SPACE above pmem, from the percentile table.

    With a HARD g this is exactly <g(P)>, so the two definitions that the tanh campaign had
    to keep apart -- and got wrong once, at the cost of a whole wave -- coincide here. Read
    from the pooled record-window percentile table, which is full-resolution and unquantised,
    never from the video stream whose stored range differs per run.
    """
    lv = np.asarray(part["flow"]["P_pctl_levels"], float)
    va = np.asarray(part["flow"]["P_pctl_values"], float)
    scale = 100.0 if lv.max() > 1.5 else 1.0
    return float(1.0 - np.interp(pmem, va, lv) / scale)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a1")
    ap.add_argument("a2b")
    ap.add_argument("--out", default="calib.json")
    a = ap.parse_args()

    # ---------------------------------------------------------------- A1: the ladder
    parts = load_parts(a.a1)
    parts = [p for p in parts if p["params"]["open_loop"]]
    parts.sort(key=lambda p: p["params"]["zeta_open"])
    zeta = parts[-1]["params"]["zeta"]
    ratios = np.array([p["params"]["zeta_open"] / zeta for p in parts])
    top = parts[int(np.argmax(ratios))]

    def lag_tau(p):
        L = p.get("lagrangian")
        return L["tau_c_lag"] if L and np.isfinite(L.get("tau_c_lag", np.nan)) else np.nan

    tau_c = lag_tau(top)
    if not np.isfinite(tau_c):
        tau_c = top["flow"]["tau_c"]
        print("  WARNING: no tracer tau_c at the top rung; falling back to the Eulerian one")
    sigma_P = top["flow"]["sigma_P"]
    u_rms = top["flow"]["u_rms"]

    print(f"A1 ladder: zeta = {zeta:g}, {len(parts)} rungs, "
          f"tau_c(zeta) = {tau_c:.0f} steps (Lagrangian), sigma_P = {sigma_P:.5f}")
    print(f"{'r':>6} {'zeta_eff':>9} {'tau_c':>8} {'tc/tc0':>7} {'sigma_P':>9} {'u_rms':>10} "
          f"{'N_def':>7} {'melt':>6} {'L_P':>6} {'std(m)':>8}")
    rungs = {}
    for p, r in zip(parts, ratios):
        fl = p["flow"]
        tl = lag_tau(p)
        tl = tl if np.isfinite(tl) else fl["tau_c"]
        print(f"{r:6.2f} {p['params']['zeta_open']:9.4g} {tl:8.0f} {tl/tau_c:7.2f} "
              f"{fl['sigma_P']:9.5f} {fl['u_rms']:10.4e} {fl['N_def']:7.0f} "
              f"{fl['melted_frac']:6.1%} {fl['L_P']:6.2f} {fl['std_m']:8.4f}")
        rungs[f"{r:g}"] = {"tau_c": tl, "sigma_P": fl["sigma_P"], "u_rms": fl["u_rms"],
                           "N_def": fl["N_def"], "L_P": fl["L_P"],
                           "melted_frac": fl["melted_frac"]}

    # ---- checkpoint 1: the floor must still be turbulent, and not too slow
    flo = min(parts, key=lambda p: abs(p["params"]["zeta_open"] / zeta - R_FLOOR))
    tl_floor = lag_tau(flo)
    tl_floor = tl_floor if np.isfinite(tl_floor) else flo["flow"]["tau_c"]
    cp1 = {"rung": flo["params"]["zeta_open"] / zeta, "N_def": flo["flow"]["N_def"],
           "tau_c_ratio": tl_floor / tau_c,
           "pass": bool(flo["flow"]["N_def"] > 20 and tl_floor / tau_c < 5.0)}
    print(f"\nCHECKPOINT 1 (floor r = {cp1['rung']:.2f}): N_def = {cp1['N_def']:.0f}, "
          f"tau_c ratio = {cp1['tau_c_ratio']:.2f} -> "
          f"{'PASS' if cp1['pass'] else 'FAIL -- raise the floor to 0.4 zeta'}")

    # ---------------------------------------------------------------- A2: the f table
    print(f"\nA2  f(pmem) per rung, pmem in units of sigma_P(zeta) = {sigma_P:.5f}")
    print(f"{'pmem':>6} " + " ".join(f"{r:>7.2f}" for r in ratios) + "  contrast")
    table = []
    for c in PMEM_COEFFS:
        fs = [f_hard(p, c * sigma_P) for p in parts]
        ft = float(np.interp(1.0, ratios, fs))
        ff = float(np.interp(R_FLOOR, ratios, fs))
        table.append({"coeff": c, "pmem": c * sigma_P, "f": fs,
                      "f_top": ft, "f_floor": ff, "contrast": ft - ff})
        print(f"{c:6.2f} " + " ".join(f"{v:7.4f}" for v in fs) + f"  {ft - ff:+.4f}")

    # ------------------------------------------------- A2b: sigma_m(tau_m), measured
    a2b = load_parts(a.a2b)
    sm = []
    for p in a2b + parts:                       # A1 supplies the tau_m = 10 tau_c points
        u = p["params"]
        if not u["open_loop"]:
            continue
        r = u["zeta_open"] / zeta
        g = u["tau_m"] / tau_c
        f_here = f_hard(p, u["pmem"])
        norm = math.sqrt(max(f_here * (1 - f_here), 1e-12))
        sm.append({"r": r, "g": g, "std_m": p["flow"]["std_m"], "f": f_here,
                   "h": p["flow"]["std_m"] / norm,
                   "closed": norm * math.sqrt(1.0 / (1.0 + g))})
    sm.sort(key=lambda d: (d["g"], d["r"]))
    print(f"\nA2b sigma_m against tau_m.  h = std(m)/sqrt(f(1-f)) divides out the f dependence,"
          f"\n    leaving the pure filter-plus-diffusion response the two rungs must share.")
    print(f"{'tau_m/tau_c':>11} {'r':>6} {'f':>7} {'std(m)':>8} {'closed':>8} {'ratio':>6} {'h':>7}")
    for d in sm:
        print(f"{d['g']:11.1f} {d['r']:6.2f} {d['f']:7.4f} {d['std_m']:8.5f} "
              f"{d['closed']:8.5f} {d['std_m']/d['closed']:6.3f} {d['h']:7.4f}")

    gs = sorted(set(round(d["g"], 3) for d in sm))
    hs = [float(np.mean([d["h"] for d in sm if abs(d["g"] - g) < 1e-3])) for g in gs]
    spread = [float(np.std([d["h"] for d in sm if abs(d["g"] - g) < 1e-3])) for g in gs]
    print(f"    h(tau_m/tau_c): " +
          "  ".join(f"{g:g}->{h:.4f}+-{s:.4f}" for g, h, s in zip(gs, hs, spread)))
    lg, lh = np.log(gs), np.log(hs)
    slope = float(np.polyfit(lg, lh, 1)[0])
    print(f"    h ~ (tau_m/tau_c)^{slope:+.3f}   (a pure one-pole filter would give -0.5 "
          f"asymptotically; anything steeper is the Dbio smoothing)")
    h_of = PchipInterpolator(lg, lh, extrapolate=True)

    def sigma_m(g, f_b):
        """sigma_m at tau_m = g tau_c and threshold-fraction f_b."""
        return math.sqrt(max(f_b * (1 - f_b), 1e-12)) * math.exp(float(h_of(math.log(g))))

    # ---------------------------------------------------- A3: the cusp, mc and tau_x
    print(f"\nA3  the cusp for r = {R_FLOOR:g}")
    print(f"{'pmem':>6} {'f_top':>7} {'f_flr':>7} {'contrast':>9} {'K_max':>7} {'chi*':>6} "
          f"{'mc':>7} {'tau_x':>7}")
    cands = []
    for row in table:
        order = np.argsort(ratios)
        fi = PchipInterpolator(ratios[order], np.asarray(row["f"])[order], extrapolate=True)
        Kmax, cst, fst = -1.0, float("nan"), float("nan")
        for i in range(1, 1000):
            cb = i / 1000.0
            aa = R_FLOOR + (1 - R_FLOOR) * (1 - cb)
            dfda = float(fi.derivative()(aa))
            K = psi(zq(cb)) * (1 - R_FLOOR) * abs(dfda)
            if K > Kmax:
                Kmax, cst, fst = K, cb, float(fi(aa))
        mc = fst + Kmax * zq(cst)
        # tau_x: the tau_m at which the measured sigma_m falls through K_max
        gl, gr = 0.05, 200.0
        if sigma_m(gl, fst) < Kmax:
            tau_x = float("nan")               # already bistable at the shortest memory
        elif sigma_m(gr, fst) > Kmax:
            tau_x = float("inf")               # never realised in reach
        else:
            for _ in range(80):
                gm = math.sqrt(gl * gr)
                if sigma_m(gm, fst) > Kmax:
                    gl = gm
                else:
                    gr = gm
            tau_x = math.sqrt(gl * gr)
        cands.append(dict(row, K_max=Kmax, chi_star=cst, f_star=fst, mc=mc, tau_x=tau_x))
        print(f"{row['coeff']:6.2f} {row['f_top']:7.4f} {row['f_floor']:7.4f} "
              f"{row['contrast']:9.4f} {Kmax:7.4f} {cst:6.3f} {mc:7.4f} {tau_x:7.2f}")

    pick = max(cands, key=lambda d: d["contrast"])
    # A MAXIMUM ON THE EDGE IS NOT A MAXIMUM. The rule is "take the pmem with the largest
    # f(zeta) - f(0.3 zeta)", and if that lands on the first or last coefficient scanned then
    # the scan, not the physics, chose it -- the true optimum is outside and every downstream
    # number (mc, the window, tau_x) inherits the truncation. Loud, because this runs
    # unattended between two waves of a 100-run campaign.
    edge = pick["coeff"] in (PMEM_COEFFS[0], PMEM_COEFFS[-1])
    print(f"\n    PICKED pmem = {pick['coeff']:g} sigma_P = {pick['pmem']:.5f}  (max contrast)")
    if edge:
        print(f"      *** WARNING: the contrast maximum sits ON THE EDGE of the scanned "
              f"range {PMEM_COEFFS[0]:g}..{PMEM_COEFFS[-1]:g} sigma_P. The optimum is "
              f"outside it and mc, the window and tau_x are all truncated by the scan. "
              f"Widen PMEM_COEFFS and re-run before generating B. ***")
    print(f"      mc = {pick['mc']:.4f}   window (f_floor, f_top) = "
          f"({pick['f_floor']:.4f}, {pick['f_top']:.4f})")
    print(f"      K_max = {pick['K_max']:.4f} at chibar* = {pick['chi_star']:.3f}")
    print(f"      tau_x = {pick['tau_x']:.2f} tau_c = {pick['tau_x']*tau_c:.0f} steps")

    grid = [0.3, 0.45, 0.68, 1, 1.5, 2.2, 3.3, 4.7, 6.8, 10, 15, 22, 30]
    below = [g for g in grid if g < pick["tau_x"]]
    print(f"      the B grid straddles it: {len(below)} points below, "
          f"{len(grid)-len(below)} above; nearest neighbours "
          f"{below[-1] if below else '-'} and "
          f"{[g for g in grid if g >= pick['tau_x']][0] if len(below) < len(grid) else '-'}")

    out = {"tau_c": tau_c, "sigma_P": sigma_P, "u_rms": u_rms,
           "zeta": zeta, "zeta0_frac": R_FLOOR,
           "pmem": pick["pmem"], "pmem_coeff": pick["coeff"], "mc": pick["mc"],
           "f_top": pick["f_top"], "f_floor": pick["f_floor"],
           "K_max": pick["K_max"], "chi_star": pick["chi_star"], "tau_x": pick["tau_x"],
           "rungs": rungs, "checkpoint1": cp1,
           "f_table": table, "sigma_m": {"points": sm, "g": gs, "h": hs,
                                         "slope": slope},
           "pmem_on_scan_edge": bool(edge),
           "note": "step campaign: chi-width = 0 and pmem-width = 0; the smoothing is "
                   "sigma_m, and tau_x is where it crosses K_max"}
    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
