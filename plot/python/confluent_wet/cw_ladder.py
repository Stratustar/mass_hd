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


def f_of(part, pmem):
    """f = 1 - CDF_P(pmem) from the run's own percentile table."""
    lv = np.asarray(part["flow"]["P_pctl_levels"], float)
    vals = np.asarray(part["flow"]["P_pctl_values"], float)
    # vals is non-decreasing in lv; np.interp needs the x axis increasing, so invert
    below = float(np.interp(pmem, vals, lv, left=0.0, right=100.0))
    return 1.0 - below / 100.0


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
                    default=[0.0, 0.2, 0.3, 0.4, 0.5, 0.6])
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

    # ---- A2: the f table, on absolute thresholds fixed by sigma_P(zeta)
    table = []
    for c in a.pmem_coeffs:
        pmem = c * sigma_P
        fs = [f_of(p, pmem) for p in parts]
        f_top = float(np.interp(zeta, zeff, fs))
        f_flr = float(np.interp(z0 * zeta, zeff, fs))
        table.append({"coeff": c, "pmem": pmem, "f": fs,
                      "f_top": f_top, "f_floor": f_flr, "contrast": f_top - f_flr})
    print(f"\nA2  f(pmem) per rung, thresholds in units of sigma_P(zeta) = {sigma_P:.4e}")
    print(f"{'c':>5} {'pmem':>11} " + " ".join(f"{r/zeta:>7.2f}" for r in zeff) + "  contrast")
    for row in table:
        print(f"{row['coeff']:5.1f} {row['pmem']:11.4e} "
              + " ".join(f"{v:7.3f}" for v in row["f"])
              + f"  {row['contrast']:+.3f}")

    # ---- A3: pick the threshold with the strongest activity contrast
    best = max(table, key=lambda r: r["contrast"])
    pmem = best["pmem"]
    fs = np.asarray(best["f"], float)
    order = np.argsort(zeff)
    f_interp = PchipInterpolator(zeff[order], fs[order], extrapolate=True)
    mc0 = 0.5 * (best["f_top"] + best["f_floor"])

    slope, chi_at, f_at = max_df_dchi(zeta, z0, f_interp)
    w_chi_0 = float(parts[-1]["params"]["chi_width"])
    G0 = slope / (2 * w_chi_0)
    if G0 < a.g_min:
        w_chi = slope / (2 * a.g_target)
        note = (f"G = {G0:.2f} < {a.g_min} at the runcard chi-width {w_chi_0:g}; "
                f"chi-width reset to |max df/dchi|/(2 G_target) = {w_chi:.4g} for G = {a.g_target}")
    else:
        w_chi, note = w_chi_0, f"G = {G0:.2f} >= {a.g_min}, chi-width kept at {w_chi_0:g}"
    G = slope / (2 * w_chi)

    roots, chig, fg, Fg = solve_loop(zeta, z0, f_interp, mc0, w_chi, a.switch_sign)
    print(f"\nA3  pmem = {pmem:.4e} ({best['coeff']:g} sigma_P), contrast "
          f"{best['contrast']:+.3f}\n    mc_0 = {mc0:.4f}   (f_floor {best['f_floor']:.3f} .. "
          f"f_top {best['f_top']:.3f})")
    print(f"    max |df/dchi| = {slope:.4g} at chi = {chi_at:.3f}\n    {note}\n"
          f"    G = {G:.2f}")
    print(f"    roots chi = {['%.4f' % r for r in roots]}"
          + ("   <- BISTABLE" if len(roots) >= 3 else "   <- single fixed point"))

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
        "pmem": pmem, "pmem_coeff": best["coeff"],
        # IQR/4, the convention carried over from cw_regime / cw_scan: the tanh in g(P) then
        # spans about half an interquartile range -- sharp enough that "above threshold"
        # really selects a fraction of the area, soft enough that the memory source still
        # has a gradient to follow instead of seeing a step.
        "pmem_width": float(0.25 * (np.interp(75.0, top["flow"]["P_pctl_levels"],
                                              top["flow"]["P_pctl_values"])
                                    - np.interp(25.0, top["flow"]["P_pctl_levels"],
                                                top["flow"]["P_pctl_values"]))),
        "mc_0": mc0, "chi_width": w_chi, "chi_width_note": note,
        "max_df_dchi": slope, "G": G,
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
