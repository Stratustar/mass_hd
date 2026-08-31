#!/usr/bin/env python3
"""Stages A2 and A3 of the memory campaign: the f table, the loop gain, and the parameters.

A1 measured, open loop, what the layer does at six prescribed activities. This script turns
those six part.json files into the closed-loop runcard, and it is where every number the
production groups need is decided. Nothing downstream re-derives any of it; calib.json is
the single source.

WHAT f IS. f(pmem; zeta_eff) is the fraction of time a lattice point spends above the
pressure threshold pmem -- equivalently, since sites are statistically equivalent, the
fraction of (site, time) samples above it, which is 1 - CDF_P(pmem). It is read from the
percentile table in part.json, i.e. from the FULL-resolution unquantised pressure field, not
from the video stream: the stream's stored range is chosen per run, so a threshold fixed in
absolute units across the ladder would sit at a different byte resolution on every rung.

f IS ALSO THE FIXED POINT OF THE MEMORY. At steady state D_t m = 0 gives m = <g(P)>, and
with a sharp enough g that is exactly the fraction of time above pmem. So the closed loop
is, at mean-field level,

    chi_bar = chi*( f( zeta_eff(chi_bar) ) ) ,   zeta_eff(chi) = zeta[z0 + (1-z0)(1-chi)]

a scalar fixed-point problem on [0,1] whose roots are the candidate phases. This ignores
var(m) -- the spatial spread of the memory -- and that is a real omission, not a technical
one: a wide var(m) smears chi*(m) and can wash out a bistability the mean-field equation
predicts. The A4 pair of runs is what tests it, and the check exists because the mean-field
prediction has failed before.

THE LOOP GAIN. Differentiating the fixed-point map,

    G = max |dF/dchi| = max|df/dchi| / (2 w_chi)     (at m = mc, where sech^2 = 1)

with df/dchi = (df/dzeta_eff) * zeta (1 - z0). G > 1 is the condition for the map to have
three roots at all; the campaign asks for G >= 1.5 and, failing that, sets w_chi so that
G = 2.5. Note which way that cuts: a SMALLER w_chi is a SHARPER switch and a LARGER gain.

Usage:
  cw_ladder.py <study> [--out calib.json] [--g-target 2.5] [--g-min 1.5]
               [--pmem-coeffs 0 0.2 0.3 0.4 0.5 0.6]
"""
import argparse
import glob
import json
import os
import sys

import numpy as np
from scipy.interpolate import PchipInterpolator

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cw_common as cw

TAU_M_GRID = (0.3, 1.0, 3.0, 10.0, 30.0)     # in units of tau_c(zeta)
TAU_CHI_RATIO = 0.3                          # tau_chi = 0.3 tau_c(zeta), fixed throughout
TAU_CHI_MIN_STEPS = 20.0                     # and never fewer than 20 integration steps
BATCHELOR_CELLS = 3.0                        # Dbio set so sqrt(Dbio tau_c) = 3 lattice cells


def load_parts(study_dir):
    """Every part.json under a study, sorted by the activity actually applied."""
    out = []
    for p in sorted(glob.glob(os.path.join(study_dir, "*", "part.json"))):
        with open(p) as fh:
            d = json.load(fh)
        u = d["params"]
        if not u["open_loop"]:
            print(f"  skipping {d['case']}: not an open-loop run", flush=True)
            continue
        out.append(d)
    out.sort(key=lambda d: d["params"]["zeta_open"])
    if len(out) < 3:
        raise RuntimeError(f"{study_dir}: found {len(out)} open-loop part.json files, "
                           f"need at least 3 to interpolate f(zeta_eff)")
    return out


def f_of(part, pmem, pmem_width):
    """f = <g(P)>, the ACTUAL fixed point of the memory -- not the time fraction above pmem.

    THIS IS THE DEFINITION THAT MATTERS, and getting it wrong cost the 20260831 A4 wave.
    The model integrates tau_m D_t m = g(P) - m with the SMOOTH

        g(P) = 1/2 [1 + tanh((P - pmem)/pmem_width)]

    so the steady state of m is <g(P)>, not 1 - CDF_P(pmem). The two agree only when
    pmem_width is small against the width of the pressure distribution. The first pass used
    the hard-step fraction with pmem_width = IQR/4 = 0.34 sigma_P(zeta), which at the
    ACTIVITY FLOOR is 1.0 times that rung's own sigma_P, with the threshold 1.5 sigma out in
    its tail. There the soft step hands partial credit to the whole bulk of the distribution
    and <g> came out 2.9x the hard-step fraction -- 0.126 against a predicted 0.044. The
    error is negligible at full activity (0.332 vs 0.327, threshold near the middle of a
    wider distribution) and enormous at the floor, i.e. it is concentrated exactly on the
    phase whose existence the whole campaign is about. Measured consequence: the chi == 1 arm
    never held, m climbed straight past mc from t = 0, and both ends converged onto the
    active phase.

    Computed from the percentile table, which is the pooled record-window pressure field at
    0.1% resolution -- exact and unquantised, unlike the video stream.
    """
    vals = np.asarray(part["flow"]["P_pctl_values"], float)
    return float(np.mean(0.5 * (1.0 + np.tanh((vals - pmem) / pmem_width))))


def solve_loop(zeta, z0, f_interp, mc, w_chi, sign=-1, n=2001):
    """Roots of chi = chi*(f(zeta_eff(chi))) on [0,1], and the loop gain."""
    chi = np.linspace(0.0, 1.0, n)
    zeff = zeta * (z0 + (1.0 - z0) * (1.0 - chi))
    f = f_interp(zeff)
    F = 0.5 * (1.0 + sign * np.tanh((f - mc) / w_chi))
    g = F - chi
    roots = []
    for i in range(n - 1):
        if g[i] == 0.0:
            roots.append(float(chi[i]))
        elif g[i] * g[i + 1] < 0:
            t = g[i] / (g[i] - g[i + 1])
            roots.append(float(chi[i] + t * (chi[i + 1] - chi[i])))
    # de-duplicate roots that the grid resolved twice
    ded = []
    for r in roots:
        if not ded or abs(r - ded[-1]) > 2.0 / n:
            ded.append(r)
    return ded, chi, f, F


def mc_window(zeta, z0, f_interp, wchi, sign, lo=0.0, hi=0.6, step=0.0025):
    """Range of mc over which the map has three roots -- the window mc* has to live in.

    THIS IS WHAT SETS THE G2 SCAN, and it is not a detail. mc* -- the threshold at which
    neither phase invades -- exists ONLY inside this window: outside it there is a single
    stable state and a front always runs to completion, so a scan point outside measures
    nothing. Measured on the 20260831 ladder the window is 0.060 wide at chi-width = 0.06,
    0.025 at 0.09 and 0.0025 at 0.12, against the campaign's assumed scan half-width of
    0.06 -- i.e. the original grid would have put four of its five points outside.
    """
    win = [mc for mc in np.arange(lo, hi, step)
           if len(solve_loop(zeta, z0, f_interp, mc, wchi, sign)[0]) >= 3]
    if not win:
        return None
    return float(min(win)), float(max(win))


def max_df_dchi(zeta, z0, f_interp, n=2001):
    """max |df/dchi| over the reachable activity range, and where it occurs."""
    chi = np.linspace(0.0, 1.0, n)
    zeff = zeta * (z0 + (1.0 - z0) * (1.0 - chi))
    f = f_interp(zeff)
    d = np.abs(np.gradient(f, chi))
    i = int(np.argmax(d))
    return float(d[i]), float(chi[i]), float(f[i])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("study", help="results dir holding the A1 <case>/part.json tree")
    ap.add_argument("--out", default=None)
    ap.add_argument("--g-min", type=float, default=1.5)
    ap.add_argument("--g-target", type=float, default=2.5)
    ap.add_argument("--pmem-coeffs", type=float, nargs="*",
                    default=[0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    ap.add_argument("--pmem-width-coeffs", type=float, nargs="*",
                    default=[0.05, 0.10, 0.15, 0.20],
                    help="g(P) smoothing widths to scan, in units of sigma_P(zeta)")
    ap.add_argument("--std-m", type=float, default=0.012,
                    help="measured spatial std of m in a closed-loop run; sets the scale a "
                         "phase basin has to beat. Default is the 20260831 A4 value")
    ap.add_argument("--chi-width", type=float, default=None,
                    help="override the ranked pick; the campaign uses this to trade some "
                         "robustness for a wider bistable mc window, which is what the G2 "
                         "scan needs")
    ap.add_argument("--mc-scan-frac", type=float, default=0.67,
                    help="G2 scans this fraction of the bistable mc half-window")
    ap.add_argument("--min-separation", type=float, default=0.4,
                    help="checkpoint 2's threshold; candidates below it are not considered")
    ap.add_argument("--switch-sign", type=int, default=-1)
    a = ap.parse_args()

    parts = load_parts(a.study)
    zeta = parts[-1]["params"]["zeta"]
    z0 = parts[-1]["params"]["zeta0_frac"]
    zeff = np.array([p["params"]["zeta_open"] for p in parts])

    # the reference rung: the top of the ladder, zeta_eff = zeta. tau_c(zeta) is the
    # SMALLEST tau_c on the ladder and is the campaign's time unit.
    top = min(parts, key=lambda p: abs(p["params"]["zeta_open"] - zeta))
    tau_c = top["flow"]["tau_c"]
    sigma_P = top["flow"]["sigma_P"]
    u_rms = top["flow"]["u_rms"]

    print(f"ladder: zeta = {zeta:g}, z0 = {z0:g}, {len(parts)} rungs")
    print(f"{'zeta_eff':>10} {'r':>6} {'tau_c':>8} {'u_rms':>10} {'N_def':>7} "
          f"{'sigma_P':>10} {'Ma':>6} {'melt':>6}")
    for p in parts:
        fl, u = p["flow"], p["params"]
        print(f"{u['zeta_open']:10.4g} {u['zeta_open']/zeta:6.2f} {fl['tau_c']:8.0f} "
              f"{fl['u_rms']:10.4e} {fl['N_def']:7.0f} {fl['sigma_P']:10.4e} "
              f"{fl['Ma']:6.3f} {fl['melted_frac']:6.1%}")

    # ---- checkpoint 1: is the activity FLOOR still turbulent?
    floor = min(parts, key=lambda p: abs(p["params"]["zeta_open"] - z0 * zeta))
    cp1 = {"rung": floor["params"]["zeta_open"], "N_def": floor["flow"]["N_def"],
           "melted": floor["flow"]["melted_frac"],
           "pass": bool(floor["flow"]["N_def"] > 20)}
    print(f"\ncheckpoint 1 (turbulence at the floor zeta_eff = {cp1['rung']:.4g}): "
          f"N_def = {cp1['N_def']:.0f} -> {'PASS' if cp1['pass'] else 'FAIL, raise z0'}")

    # ---- A2: the f table. Thresholds are absolute pressures fixed by sigma_P(zeta), and
    # ---- the SMOOTHING WIDTH is scanned beside them, because f depends on both.
    wid_coeffs = a.pmem_width_coeffs
    table = []
    for c in a.pmem_coeffs:
        for wc in wid_coeffs:
            pmem, pw = c * sigma_P, wc * sigma_P
            fs = [f_of(p, pmem, pw) for p in parts]
            f_top = float(np.interp(zeta, zeff, fs))
            f_flr = float(np.interp(z0 * zeta, zeff, fs))
            table.append({"coeff": c, "width_coeff": wc, "pmem": pmem, "pmem_width": pw,
                          "f": fs, "f_top": f_top, "f_floor": f_flr,
                          "contrast": f_top - f_flr})
    print(f"\nA2  f = <g(P)> per rung; pmem and pmem-width in units of "
          f"sigma_P(zeta) = {sigma_P:.4e}")
    print(f"{'pmem':>6} {'width':>6} " + " ".join(f"{r/zeta:>7.2f}" for r in zeff) + "  contrast")
    for row in table:
        print(f"{row['coeff']:6.2f} {row['width_coeff']:6.2f} "
              + " ".join(f"{v:7.3f}" for v in row["f"])
              + f"  {row['contrast']:+.3f}")

    # ---- A3: pick the operating point by ROBUSTNESS, not by contrast alone.
    #
    # The campaign's rule is "take the pmem with the largest f(zeta) - f(0.3 zeta)", and
    # contrast is reported below because it is the right first instinct: no contrast, no
    # loop. But contrast does not decide whether the two phases SURVIVE, and checkpoint 2
    # tests survival. Two things beyond contrast decide it:
    #
    #   * the BASIN of the weaker phase, min(chi_mid - chi_lo, chi_hi - chi_mid). A pair of
    #     roots with an unstable root hugging one of them is bistable on paper and not in a
    #     simulation.
    #   * the size of a MEMORY FLUCTUATION in the same units. A fluctuation std(m) moves the
    #     phenotype target by std(m)/(2 w_chi), so sharpening the switch widens the basin and
    #     amplifies the fluctuation in exact proportion. That is why "make the switch harder"
    #     is not a free way to buy bistability, and why ranking by basin alone drives w_chi to
    #     zero and lies.
    #
    # So rank by  R = (weaker basin) / (std(m) / (2 w_chi)),  the number of memory-fluctuation
    # widths that fit inside the basin the weaker phase has to hold. R has an interior
    # maximum in w_chi. std(m) is MEASURED (--std-m, from a closed-loop run); the default is
    # the 20260831 A4 value, and it is the one input here that a fresh operating point can
    # change, so it is re-measured with every A4 wave.
    cands = []
    for row in table:
        fs = np.asarray(row["f"], float)
        fi = PchipInterpolator(zeff[np.argsort(zeff)], fs[np.argsort(zeff)], extrapolate=True)
        mc = 0.5 * (row["f_top"] + row["f_floor"])
        slope_r, _, _ = max_df_dchi(zeta, z0, fi)
        for wchi in np.arange(0.02, 0.3001, 0.002):
            rt, _, _, _ = solve_loop(zeta, z0, fi, mc, wchi, a.switch_sign)
            if len(rt) < 3:
                continue
            lo, mid, hi = rt[0], rt[len(rt) // 2], rt[-1]
            sep = hi - lo
            if sep < a.min_separation:
                continue
            basin = min(mid - lo, hi - mid)
            R = basin / (a.std_m / (2 * wchi))
            cands.append({"coeff": row["coeff"], "width_coeff": row["width_coeff"],
                          "pmem": row["pmem"], "pmem_width": row["pmem_width"],
                          "mc_0": mc, "chi_width": float(wchi), "G": slope_r / (2 * wchi),
                          "max_df_dchi": slope_r, "roots": rt, "separation": sep,
                          "basin": basin, "robustness": R, "f": row["f"],
                          "f_top": row["f_top"], "f_floor": row["f_floor"],
                          "contrast": row["contrast"]})
    if not cands:
        raise RuntimeError("no (pmem, pmem-width, chi-width) in the scan gives three roots "
                           f"with separation >= {a.min_separation}; the loop is not bistable "
                           "anywhere on this ladder")
    cands.sort(key=lambda d: -d["robustness"])
    print(f"\nA3  top operating points, ranked by R = basin / (chi* shift per std(m) = "
          f"{a.std_m:g})")
    print(f"{'pmem':>6} {'width':>6} {'w_chi':>6} {'mc_0':>6} {'f_flr':>6} {'f_top':>6} "
          f"{'chi_lo':>7} {'chi_mid':>8} {'chi_hi':>7} {'sep':>6} {'basin':>6} {'G':>5} {'R':>5}")
    seen, shown = set(), []
    for d in cands:
        k = (d["coeff"], d["width_coeff"])
        if k in seen:
            continue
        seen.add(k); shown.append(d)
        r = d["roots"]
        print(f"{d['coeff']:6.2f} {d['width_coeff']:6.2f} {d['chi_width']:6.3f} {d['mc_0']:6.3f} "
              f"{d['f_floor']:6.3f} {d['f_top']:6.3f} {r[0]:7.4f} {r[len(r)//2]:8.4f} "
              f"{r[-1]:7.4f} {d['separation']:6.3f} {d['basin']:6.3f} {d['G']:5.2f} "
              f"{d['robustness']:5.1f}")
        if len(shown) >= 8:
            break

    pick = cands[0]
    best = {"coeff": pick["coeff"], "pmem": pick["pmem"], "f_top": pick["f_top"],
            "f_floor": pick["f_floor"], "contrast": pick["contrast"], "f": pick["f"]}
    pmem, pmem_width = pick["pmem"], pick["pmem_width"]
    mc0, w_chi = pick["mc_0"], pick["chi_width"]
    slope = pick["max_df_dchi"]
    G = pick["G"]
    note = (f"chosen by robustness R = {pick['robustness']:.1f}; contrast {pick['contrast']:+.3f} "
            f"(the largest contrast in the scan is "
            f"{max(t['contrast'] for t in table):+.3f})")
    fs = np.asarray(pick["f"], float)
    order = np.argsort(zeff)
    f_interp = PchipInterpolator(zeff[order], fs[order], extrapolate=True)
    roots = pick["roots"]
    print(f"\n    PICKED  pmem = {pmem:.4e} ({pick['coeff']:g} sigma_P), "
          f"pmem-width = {pmem_width:.4e} ({pick['width_coeff']:g} sigma_P)")
    print(f"            mc_0 = {mc0:.4f}   chi-width = {w_chi:.4f}   G = {G:.2f}   "
          f"R = {pick['robustness']:.1f}")
    print(f"            roots chi = {['%.4f' % r for r in roots]}  "
          f"separation {pick['separation']:.3f}  weaker basin {pick['basin']:.3f}")
    print(f"            {note}")

    # the chi-width scan AT THE PICKED THRESHOLD, so the A4 wave can bracket it
    g_scan = []
    for wchi in (0.04, 0.05, 0.06, 0.07, 0.088, 0.10, 0.12, 0.15, 0.20):
        rt, _, _, _ = solve_loop(zeta, z0, f_interp, mc0, wchi, a.switch_sign)
        b = (min(rt[len(rt)//2] - rt[0], rt[-1] - rt[len(rt)//2])
             if len(rt) >= 3 else 0.0)
        g_scan.append({"G": slope / (2 * wchi), "chi_width": wchi, "n_roots": len(rt),
                       "roots": rt, "separation": (rt[-1]-rt[0]) if len(rt) >= 3 else 0.0,
                       "basin": b, "robustness": b / (a.std_m / (2 * wchi))})
    for e in g_scan:
        win = mc_window(zeta, z0, f_interp, e["chi_width"], a.switch_sign)
        e["mc_window"] = list(win) if win else None
        e["mc_window_width"] = (win[1] - win[0]) if win else 0.0
    print("\n    chi-width scan at the picked threshold:")
    print(f"    {'w_chi':>7} {'G':>5} {'roots':>6} {'chi_lo':>8} {'chi_mid':>8} {'chi_hi':>8} "
          f"{'sep':>6} {'basin':>6} {'R':>5} {'mc window':>10}")
    for e in g_scan:
        r = e["roots"]
        cols = (f"{r[0]:8.4f} {r[len(r)//2]:8.4f} {r[-1]:8.4f} {e['separation']:6.3f} "
                f"{e['basin']:6.3f} {e['robustness']:5.1f}") if e["n_roots"] >= 3 else \
               (f"{r[0]:8.4f} {'--':>8} {'--':>8} {'--':>6} {'--':>6} {'--':>5}" if r else "")
        print(f"    {e['chi_width']:7.3f} {e['G']:5.2f} {e['n_roots']:6d} {cols} "
              f"{e['mc_window_width']:10.4f}")

    # ---- the chi-width the campaign actually runs at, and the G2 scan grid it implies
    if a.chi_width is not None:
        w_chi = a.chi_width
        roots, _, _, _ = solve_loop(zeta, z0, f_interp, mc0, w_chi, a.switch_sign)
        G = slope / (2 * w_chi)
        note += (f"; chi-width OVERRIDDEN to {w_chi:g} on the command line")
    win = mc_window(zeta, z0, f_interp, w_chi, a.switch_sign)
    if win is None:
        raise RuntimeError(f"chi-width = {w_chi} has no bistable mc window at all")
    half = 0.5 * (win[1] - win[0])
    # Scan the middle two thirds of the window: mean-field edges are where bistability
    # DISAPPEARS, and var(m) shrinks the real window inward from them, so a scan point
    # placed on a mean-field edge is a point that measures a single-phase run.
    span = a.mc_scan_frac * half
    mc_offsets = [-span, -span / 2, 0.0, span / 2, span]
    print(f"\n    chi-width in use = {w_chi:g}; bistable mc window "
          f"[{win[0]:.4f}, {win[1]:.4f}] (width {2*half:.4f}, mc_0 sits "
          f"{(mc0-win[0])/(2*half):.0%} of the way across)")
    print(f"    G2 mc scan = mc_0 + " + ", ".join(f"{o:+.4f}" for o in mc_offsets))

    if len(roots) >= 3:
        chi_lo, chi_mid, chi_hi = roots[0], roots[len(roots) // 2], roots[-1]
    else:
        chi_lo = chi_mid = chi_hi = roots[0] if roots else float("nan")

    def f_at_chi(c):
        return float(f_interp(zeta * (z0 + (1 - z0) * (1 - c))))

    # G1 needs a DIFFERENT threshold from the one A3 just chose. Its question is the
    # scaling of chi's statistics with tau_m when the phenotype is a smoothed copy of the
    # pressure, and for that the uniform state has to be an EXACT fixed point -- otherwise
    # the loop is driven by a misplaced threshold rather than by the fluctuations, which is
    # what happened to the 20260826 batch (pmem = 0 against a true median of +1.3e-2). With
    # chi == 0.5 and mc = 0.5 the fixed point needs f = 0.5, i.e. pmem = median(P) at the
    # activity a half-switched layer actually runs at, zeta_eff(0.5), NOT at the top rung.
    zeff_g1 = zeta * (z0 + (1 - z0) * 0.5)
    med = np.array([np.interp(50.0, p["flow"]["P_pctl_levels"], p["flow"]["P_pctl_values"])
                    for p in parts])
    pmem_g1 = float(PchipInterpolator(zeff[order], med[order], extrapolate=True)(zeff_g1))

    tau_chi = max(TAU_CHI_RATIO * tau_c, TAU_CHI_MIN_STEPS)
    Dbio = BATCHELOR_CELLS**2 / tau_c

    calib = {
        "study": a.study,
        "zeta": zeta, "zeta0_frac": z0, "switch_sign": a.switch_sign,
        "tau_c": tau_c, "sigma_P": sigma_P, "u_rms": u_rms,
        "N_def_top": top["flow"]["N_def"], "d_top": top["flow"]["d"],
        "L_P_top": top["flow"]["L_P"], "Ma_top": top["flow"]["Ma"],
        "checkpoint1": cp1,
        "ladder": [{"zeta_eff": p["params"]["zeta_open"], "ratio": p["params"]["zeta_open"]/zeta,
                    "tau_c": p["flow"]["tau_c"], "u_rms": p["flow"]["u_rms"],
                    "sigma_P": p["flow"]["sigma_P"], "N_def": p["flow"]["N_def"],
                    "Ma": p["flow"]["Ma"], "melted": p["flow"]["melted_frac"],
                    "L_P": p["flow"]["L_P"], "case": p["case"]} for p in parts],
        "f_table": table,
        "pmem": pmem, "pmem_coeff": pick["coeff"],
        "pmem_width": pmem_width, "pmem_width_coeff": pick["width_coeff"],
        "robustness": pick["robustness"], "basin": pick["basin"],
        "separation": pick["separation"], "std_m_assumed": a.std_m,
        "mc_0": mc0, "chi_width": w_chi, "chi_width_note": note,
        "mc_window": list(win), "mc_scan_offsets": mc_offsets,
        "candidates": shown,
        "max_df_dchi": slope, "G": G,
        "g_scan": g_scan,
        "roots": roots,
        "chi_lo": chi_lo, "chi_mid": chi_mid, "chi_hi": chi_hi,
        "f_lo": f_at_chi(chi_lo), "f_mid": f_at_chi(chi_mid), "f_hi": f_at_chi(chi_hi),
        "f_at_full_zeta": best["f_top"], "f_at_floor": best["f_floor"],
        "pmem_g1": pmem_g1, "zeta_eff_g1": zeff_g1,
        "tau_chi": tau_chi, "tau_m_grid": [r * tau_c for r in TAU_M_GRID],
        "tau_m_ratios": list(TAU_M_GRID),
        "Dbio": Dbio,
        "video": {"nvideo": max(1, int(round(0.2 * tau_c))),
                  "ninfo": max(1, int(round(5.0 * tau_c))),
                  "p_scale": 6.0 * sigma_P, "u_scale": 6.0 * u_rms},
    }
    out = a.out or os.path.join(a.study, "calib.json")
    with open(out, "w") as fh:
        json.dump(calib, fh, indent=1)
    print(f"\ntau_c(zeta) = {tau_c:.0f} steps   tau_chi = {tau_chi:.0f}   "
          f"Dbio = {Dbio:.4g}\ntau_m grid = "
          + ", ".join(f"{r:g}tc={r*tau_c:.0f}" for r in TAU_M_GRID))
    print(f"video: nvideo = {calib['video']['nvideo']}, ninfo = {calib['video']['ninfo']}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
