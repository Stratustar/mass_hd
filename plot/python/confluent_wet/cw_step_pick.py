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


ACTIVE_ROOT, PASSIVE_ROOT = 0.2, 0.8      # what "the two phases coexist" has to mean


def coexists(g, f_of, sigma_m, mc, r, n=2001):
    """Does tau_m = g tau_c admit BOTH an active and a passive stable phase?

    THE FIXED POINT WITH sigma_m CARRIED. chibar = Phi((mc - f(a))/sigma_m(g, a)) with
    a = r + (1-r)(1-chibar). The cusp condition sigma_m < K_max that picks mc treats
    sigma_m as a constant, and it is NOT one: sigma_m grows as chi rises, because the
    activity falls, L_P grows, and the diffusion removes proportionally less. That is a
    DESTABILISING feedback the cusp analysis cannot see, and it is why the cusp's tau_x
    came out well below the measured one.

    AND THREE ROOTS IS NOT ENOUGH. Carrying sigma_m(a) produces triplets that sit entirely
    inside the passive branch (0.73/0.83/0.97 at 4.7 tau_c) -- real structure, perhaps, but
    not the coexistence the campaign is about. The test is that a stable root exists below
    ACTIVE_ROOT and another above PASSIVE_ROOT.

    chibar = 1 is an exact fixed point whenever f(r) < mc, since Phi of a large positive
    number is 1; it is counted as the passive root when it is stable, which a sign change
    cannot detect at the boundary.
    """
    x = np.linspace(0.0, 1.0, n)
    F = np.array([norm_cdf((mc - f_of(r + (1 - r) * (1 - v))) /
                           sigma_m(g, r + (1 - r) * (1 - v))) for v in x])
    d = F - x
    stable = []
    for i in range(n - 1):
        if d[i] * d[i + 1] < 0:
            t = d[i] / (d[i] - d[i + 1])
            root = float(x[i] + t * (x[i + 1] - x[i]))
            if d[i] > 0 and d[i + 1] < 0:          # F crosses from above: stable
                stable.append(root)
    # BOTH ENDS need the boundary treatment, not just chibar = 1. As sigma_m falls the
    # active root collapses onto chibar = 0, and there d[0] = F(0) - 0 is zero or a
    # rounding of it -- no sign change to find, so a pure crossing search silently loses
    # the very root the campaign is about. (That is what made tau_x = 100 tau_c report
    # "not coexisting" while 30 did.) A root inside the first cell is the same statement.
    if F[0] <= 1.0 / (n - 1):
        stable.append(float(F[0]))
    if d[-1] >= -1e-9 and F[-1] > 1 - 1e-6:        # chibar = 1 is a stable boundary root
        stable.append(1.0)
    return (any(v < ACTIVE_ROOT for v in stable) and
            any(v > PASSIVE_ROOT for v in stable)), stable


def norm_cdf(z):
    return Phi(float(z))


def _tau_x_flat(sigma_m, Kmax, f_b):
    """tau_x from the rung-AVERAGED h -- reported only so the correction is visible."""
    gl, gr = 0.05, 200.0
    if sigma_m(gl, f_b) < Kmax:
        return float("nan")
    if sigma_m(gr, f_b) > Kmax:
        return float("inf")
    for _ in range(80):
        gm = math.sqrt(gl * gr)
        if sigma_m(gm, f_b) > Kmax:
            gl = gm
        else:
            gr = gm
    return math.sqrt(gl * gr)


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
    print(f"    h(tau_m/tau_c), averaged over rungs: " +
          "  ".join(f"{g:g}->{h:.4f}+-{s:.4f}" for g, h, s in zip(gs, hs, spread)))
    lg, lh = np.log(gs), np.log(hs)
    slope = float(np.polyfit(lg, lh, 1)[0])
    print(f"    h ~ (tau_m/tau_c)^{slope:+.3f}   (a pure one-pole filter would give -0.5 "
          f"asymptotically; anything steeper is the Dbio smoothing)")

    # h STILL DEPENDS ON ACTIVITY after sqrt(f(1-f)) is divided out, and the residual is not
    # noise. What the diffusion removes is set by l_D/L_P, and L_P runs from 3.7 at full
    # activity to 8.5 at the floor -- so the same tau_m smooths a low-activity field
    # proportionally less. Measured: h(tau_m = tau_c) is 0.700 at a = 0.45 against 0.511 at
    # a = 1.0, and the gap WIDENS with tau_m (0.096 vs 0.041 at 30 tau_c) exactly as
    # l_D ~ sqrt(tau_m) says it should.
    #
    # This matters because the cusp does not sit at either rung: chibar* ~ 0.7 puts it at
    # a* = r + (1-r)(1-chibar*) ~ 0.51. Averaging the rungs would evaluate sigma_m at the
    # wrong activity and pull tau_x low. Interpolate in BOTH variables instead: log h is
    # linear in log a to 1% over the five rungs measured at tau_m = 10 tau_c.
    by_g = {}
    for g in gs:
        pts = sorted(((d["r"], d["h"]) for d in sm if abs(d["g"] - g) < 1e-3))
        by_g[g] = (np.log([r for r, _ in pts]), np.log([h for _, h in pts]))
    a_exp = {g: (float(np.polyfit(la, lh_, 1)[0]) if len(la) > 1 else 0.0)
             for g, (la, lh_) in by_g.items()}
    print(f"    residual activity dependence  h ~ a^p:  " +
          "  ".join(f"{g:g}->{p:+.3f}" for g, p in a_exp.items()))

    def h_at(g, a):
        """log h interpolated in log a at each measured tau_m, then in log tau_m."""
        per = []
        for gg in gs:
            la, lh_ = by_g[gg]
            per.append(float(np.interp(math.log(a), la, lh_)) if len(la) > 1
                       else float(lh_[0]))
        return math.exp(float(PchipInterpolator(lg, per, extrapolate=True)(math.log(g))))

    def sigma_m(g, f_b, a=None):
        """sigma_m at tau_m = g tau_c, threshold-fraction f_b, activity ratio a."""
        h = h_at(g, a) if a is not None else math.exp(float(
            PchipInterpolator(lg, lh, extrapolate=True)(math.log(g))))
        return math.sqrt(max(f_b * (1 - f_b), 1e-12)) * h

    # ---------------------------------------------------- A3: the cusp, mc and tau_x
    print(f"\nA3  the cusp for r = {R_FLOOR:g}")
    print(f"{'pmem':>6} {'f_top':>7} {'f_flr':>7} {'contrast':>9} {'K_max':>7} {'chi*':>6} "
          f"{'a*':>6} {'mc':>7} {'tau_x':>7}")
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
        # tau_x: the tau_m at which sigma_m falls through K_max, evaluated AT THE CUSP --
        # its own activity a_star as well as its own f. Both matter and they pull the same
        # way: the cusp sits at lower activity than the top rung, where L_P is larger and
        # the diffusion has removed less.
        a_star = R_FLOOR + (1 - R_FLOOR) * (1 - cst)
        gl, gr = 0.05, 200.0
        if sigma_m(gl, fst, a_star) < Kmax:
            tau_x = float("nan")               # already bistable at the shortest memory
        elif sigma_m(gr, fst, a_star) > Kmax:
            tau_x = float("inf")               # never realised in reach
        else:
            for _ in range(80):
                gm = math.sqrt(gl * gr)
                if sigma_m(gm, fst, a_star) > Kmax:
                    gl = gm
                else:
                    gr = gm
            tau_x = math.sqrt(gl * gr)
        cands.append(dict(row, K_max=Kmax, chi_star=cst, f_star=fst, mc=mc, tau_x=tau_x,
                          a_star=a_star,
                          tau_x_rung_avg=_tau_x_flat(sigma_m, Kmax, fst)))
        print(f"{row['coeff']:6.2f} {row['f_top']:7.4f} {row['f_floor']:7.4f} "
              f"{row['contrast']:9.4f} {Kmax:7.4f} {cst:6.3f} {a_star:6.3f} {mc:7.4f} "
              f"{tau_x:7.2f}")

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
    print(f"      tau_x = {pick['tau_x']:.2f} tau_c = {pick['tau_x']*tau_c:.0f} steps "
          f"(at the cusp activity a* = {pick['a_star']:.3f}; averaging the rungs instead "
          f"would have said {pick['tau_x_rung_avg']:.2f})")

    # ---- the same question with sigma_m(a) carried, which is the honest version
    order = np.argsort(ratios)
    fi_pick = PchipInterpolator(ratios[order], np.asarray(pick["f"])[order],
                                extrapolate=True)

    def f_of_a(a):
        return float(fi_pick(a))

    def sig(g, a):
        return sigma_m(g, f_of_a(a), a)

    grid = [0.3, 0.45, 0.68, 1, 1.5, 2.2, 3.3, 4.7, 6.8, 10, 15, 22, 30, 47, 68, 100]
    print("\n    the fixed point with sigma_m(tau_m, a) CARRIED, not frozen at the cusp:")
    print(f"    {'tau_m/tau_c':>11} {'sigma_m(a=1)':>13} {'coexists':>9}   stable roots")
    tau_x_co = None
    for g in grid:
        ok, st = coexists(g, f_of_a, sig, pick["mc"], R_FLOOR)
        if ok and tau_x_co is None:
            tau_x_co = g
        print(f"    {g:11g} {sig(g, 1.0):13.4f} {str(ok):>9}   "
              + "  ".join(f"{v:.4f}" for v in st))
    if tau_x_co:
        lo, hi = 1.0, float(tau_x_co)
        for _ in range(40):
            mid = math.sqrt(lo * hi)
            if coexists(mid, f_of_a, sig, pick["mc"], R_FLOOR)[0]:
                hi = mid
            else:
                lo = mid
        tau_x_co = math.sqrt(lo * hi)
        print(f"\n    COEXISTENCE THRESHOLD = {tau_x_co:.1f} tau_c, against the cusp's "
              f"{pick['tau_x']:.2f}. The cusp treats sigma_m as constant; carrying its "
              f"activity dependence -- a destabilising feedback -- moves the threshold up.")
    else:
        print("\n    NO COEXISTENCE anywhere on the grid with sigma_m carried")

    ref = tau_x_co if tau_x_co else pick["tau_x"]
    below = [g for g in grid if g < ref]
    print(f"      the B grid straddles it: {len(below)} points below, "
          f"{len([g for g in grid if g >= ref and g <= 30])} above (within the scanned "
          f"0.3-30); nearest neighbours {below[-1] if below else '-'} and "
          f"{[g for g in grid if g >= ref][0] if len(below) < len(grid) else '-'}")

    out = {"tau_c": tau_c, "sigma_P": sigma_P, "u_rms": u_rms,
           "zeta": zeta, "zeta0_frac": R_FLOOR,
           "pmem": pick["pmem"], "pmem_coeff": pick["coeff"], "mc": pick["mc"],
           "f_top": pick["f_top"], "f_floor": pick["f_floor"],
           "K_max": pick["K_max"], "chi_star": pick["chi_star"], "tau_x": pick["tau_x"],
           "a_star": pick["a_star"], "tau_x_rung_avg": pick["tau_x_rung_avg"],
           "tau_x_coexistence": tau_x_co,
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
