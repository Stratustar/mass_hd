#!/usr/bin/env python3
"""The S-curve: the mean-field fixed point of the step switch, against the 80 measured runs.

    cw_s3_theory.py <bparts_dir> <aparts_dir> --calib calib_s3.json --old calib_step.json
                    --out DIR

THE THEORY, in the order the assumptions enter.

1. Along a material path the memory is a one-pole filter of the telegraph signal g(P(t)):
       D_t m = (g - m)/tau_m ,   g = Theta(P - pmem) in {0, 1}
   with duty cycle f = <g> and correlation time tau_c. Its stationary mean is f and, for an
   exponentially correlated drive, its variance is f(1-f) tau_c/(tau_c + tau_m). Cell
   motility (Dbio) smooths m over l_D = sqrt(2 Dbio tau_m) and takes the variance further
   down; the campaign measures the ratio h = sigma_m/sqrt(f(1-f)) rather than predicting it.

2. The phenotype follows the memory on tau_chi = 0.3 tau_c << tau_m, so chi ~ chi*(m) =
   Theta(mc - m) pointwise, and the space average is a FRACTION:
       chibar = P(m < mc)
   -- exactly, no closure yet: chi relaxes linearly to chi* and diffuses conservatively, so
   <chi> = <chi*> in the steady state. What needs a closure is the DISTRIBUTION of m. Taking
   it normal with mean <m> and spread sigma_m gives the campaign's form
       chibar = Phi((mc - <m>)/sigma_m).

3. The loop closes through the activity. chi = 1 is the passive phenotype, zeta_eff =
   zeta [r + (1-r)(1-chi)], so a homogeneous layer at chibar runs at the mean activity
       abar = r + (1-r)(1-chibar)
   and the S-curve is the map
       F(chibar) = Phi( (mc - <m>(abar)) / sigma_m(tau_m, abar) ) ,
   whose crossings with the diagonal are the fixed points: stable where F crosses from above,
   unstable where it crosses from below. Three crossings = bistable.

4. Bistability is a slope condition. Holding sigma_m fixed,
       F'(chibar) = phi(z) (1-r) |d<m>/da| / sigma_m ,   z = Phi^-1(chibar)
   so three roots need sigma_m < K(chibar) = phi(z)(1-r)|d<m>/da| somewhere. sigma_m falls
   with tau_m; tau_x is where it drops through K_max. Carrying sigma_m's own activity
   dependence (it GROWS as chi rises, because L_P grows and the diffusion removes less) is a
   destabilising feedback the slope condition cannot see, and it raises the threshold.

WHAT IS MEASURED HERE, and what each measurement is allowed to assume.

  * THE CLOSURE ITSELF, no model of anything: chibar_measured against Phi((mc - <m>)/std(m))
    with <m> and std(m) read from the SAME frames of the SAME run, for all 80 closed-loop
    runs and the 6 open-loop rungs. Its residual is the error of the normal assumption in
    step 2 alone.

  * <m> AGAINST f. In the continuum the two are equal in the steady state of a periodic
    incompressible box: average the transport equation and every term but the source drops.
    Measured, <m> sits ~8% below f in the open loop and in the closed loop alike. The
    identity that survives discretisation is
       <m> = <g> + tau_m <m div u>
    and the centred-difference advection of an intensive scalar is not exactly conservative
    on the lattice, so <m div u> is small per step and tau_m multiplies it. The deficit is
    reported per run against tau_m sigma_m, which is the scaling that argument predicts.

  * sigma_m(tau_m, a), data-driven: log h is PCHIP in log tau_m through the a = 1 anchors
    (the L = 800 A2b points at 1 and 3 tau_c -- transferred, the rescaling leaves l_D and
    L_P alone and the transfer checks to 2% at 10 and 30 tau_c -- plus this campaign's
    a = 1 rung and its chi0 branch wherever it sits at chibar < 0.02), with h ~ a^beta(tau_m)
    and beta from the six-rung ladders.

  * THE PREDICTED BRANCHES. Roots of F on a fine tau_m grid, three variants:
       A  normal law, <m>(a) = f(a)        the campaign's A3 form: the theoretical mean of m
       B  normal law, <m>(a) = m_open(a)   the mean the open loop actually produces
       C  Beta law,   <m>(a) = m_open(a)   the stationary law of a filtered telegraph
                                           signal, which is Bernoulli at short tau_m
                                           (chibar -> 1 - f) and normal at long
    and, at tau_m = 10 tau_c only, the DIRECT map F(chibar) = chi_open(abar) read straight
    off the six open-loop rungs -- no distributional assumption, no sigma_m model -- which
    is the cleanest single statement of where the bifurcation sits.

WHAT CAME OUT (2026-09-05, the 80-run scan). The normal closure fails at short tau_m for a
structural reason: m is bounded and nearly two-valued there, and Phi cannot reach the
1 - f limit. The Beta law, with the SAME two moments, halves the closure residual and
removes its bias, and puts the mixed branch within 0.08 of the measurement everywhere.
All three variants and the direct map put the saddle-node at 7.7-10 tau_c; the scan
separates at 12.0 +/- 0.8. The gap is a LIFETIME, not a missing root: chi1 at 8-11 tau_c
plateaus at 0.85-0.9 for tens of tau_c and then collapses, at 47 / 69 / 121 / > 300 tau_c
for 8.1 / 9.7 / 11.2 / 12.8 -- the passive phase exists where mean field says, and is
eroded by the active one until the erosion outlives the run.
"""
import argparse
import glob
import json
import math
import os
import sys

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.stats import norm, beta as beta_dist


def beta_cdf(x, mean, std):
    """P(m < x) for a Beta distribution with the given mean and spread.

    THE ONE-POLE FILTER OF A TWO-STATE MARKOV PROCESS HAS A BETA STATIONARY LAW. With
    switching rates k_on, k_off and memory time tau_m the density of m is
    Beta(k_on tau_m, k_off tau_m): mean f = k_on/(k_on + k_off), variance f(1-f)/(1 + s),
    s = (k_on + k_off) tau_m. It interpolates the two limits a normal law cannot: at
    tau_m << tau_c it is the Bernoulli(f) of g itself, so P(m < mc) -> 1 - f; at
    tau_m >> tau_c it is normal. Matching the measured mean and spread fixes s = f(1-f)/
    sigma^2 - 1, so this closure has NO extra parameter over the normal one."""
    v = std * std
    s = mean * (1 - mean) / v - 1.0
    if s <= 0:                       # spread at the Bernoulli limit or beyond: two-point law
        return 1.0 - mean
    return float(beta_dist.cdf(x, mean * s, (1 - mean) * s))

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

R_FLOOR = 0.3
WINDOW_TAU_C = 100.0
ACTIVE_ROOT, PASSIVE_ROOT = 0.2, 0.8
G_GRID = np.exp(np.linspace(math.log(0.25), math.log(40.0), 240))


# ------------------------------------------------------------------- reading

def load(d):
    out = []
    for p in sorted(glob.glob(os.path.join(d, "*", "part.json"))):
        with open(p) as fh:
            q = json.load(fh)
        q["_case"] = os.path.basename(os.path.dirname(p))
        sj = os.path.join(os.path.dirname(p), "series.json")
        q["_series"] = json.load(open(sj)) if os.path.exists(sj) else None
        out.append(q)
    if not out:
        raise RuntimeError(f"no part.json under {d}")
    return out


def start_of(u):
    cc = u["chi_config"]
    if cc == "stripe":
        return "leftright"
    if cc == "binary-noise":
        return "patches"
    return "chi1" if u["chi0"] > 0.5 else "chi0"


def window_moments(q, tau_c):
    """<m>, std(m), <chi>, std(chi) over the frames in the closing WINDOW_tau_c."""
    t = np.asarray(q["flow"]["times"], float)
    pf = q["flow"]["per_frame"]
    keep = t >= t[-1] - WINDOW_TAU_C * tau_c
    if keep.sum() < 3:
        keep = np.ones_like(t, bool)
    m = np.asarray(pf["m_mean"], float)[keep]
    s = np.asarray(pf["std_m"], float)[keep]
    c = np.asarray(pf["chi_mean"], float)[keep]
    sc = np.asarray(pf["std_chi"], float)[keep]
    return {"m": float(m.mean()), "std_m": float(s.mean()), "chi": float(c.mean()),
            "std_chi": float(sc.mean()), "frames": int(keep.sum()),
            "chi_frames": c, "drift": float(c[len(c)//2:].mean() - c[:len(c)//2].mean())}


# --------------------------------------------------------------- the S-curve

def roots_of(F, n=4001):
    """(stable, unstable) fixed points of F on [0,1], with the boundary cases.

    A crossing from above the diagonal to below is stable. chibar = 0 and 1 are fixed
    points whenever F pins there (F(0) ~ 0 means Phi of a large negative number: the active
    phase saturated), and a pinned end is stable when F stays on its side of the diagonal
    in the first cell."""
    x = np.linspace(0.0, 1.0, n)
    Fx = np.array([F(v) for v in x])
    d = Fx - x
    stable, unstable = [], []
    for i in range(n - 1):
        if d[i] == 0.0:
            continue
        if d[i] * d[i + 1] < 0:
            t = d[i] / (d[i] - d[i + 1])
            r = float(x[i] + t * (x[i + 1] - x[i]))
            (stable if d[i] > 0 else unstable).append(r)
    eps = 1.0 / (n - 1)
    if Fx[0] <= eps and d[1] < 0:
        stable.insert(0, float(Fx[0]))
    if Fx[-1] >= 1 - eps and d[-2] > 0:
        stable.append(float(Fx[-1]))
    return stable, unstable


def coexist(stable):
    return (any(v < ACTIVE_ROOT for v in stable) and any(v > PASSIVE_ROOT for v in stable))


def first_g(G, flags):
    """Geometric midpoint between the last False and the first True."""
    flags = list(flags)
    if True not in flags:
        return float("nan")
    i = flags.index(True)
    return float(G[i]) if i == 0 else float(math.sqrt(G[i - 1] * G[i]))


# ----------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bparts")
    ap.add_argument("aparts")
    ap.add_argument("--calib", required=True)
    ap.add_argument("--old", required=True, help="calib_step.json of the L = 800 campaign")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    cal = json.load(open(a.calib))
    old = json.load(open(a.old))
    tau_c, mc, pmem = float(cal["tau_c"]), float(cal["mc"]), float(cal["pmem"])
    r = R_FLOOR

    # ============================================================ stage A: the ladder
    A = load(a.aparts)
    lad = []
    for q in A:
        u = q["params"]
        aa = u["zeta_open"] / u["zeta"]
        lv = np.asarray(q["flow"]["P_pctl_levels"], float)
        va = np.asarray(q["flow"]["P_pctl_values"], float)
        f_here = float(1.0 - np.interp(pmem, va, lv) / 100.0)
        fl = q["flow"]
        lad.append({"a": aa, "f": f_here, "m": fl["m_mean"], "std_m": fl["std_m"],
                    "chi": fl["chi_mean"], "L_P": fl["L_P"],
                    "h": fl["std_m"] / math.sqrt(f_here * (1 - f_here)),
                    "g": u["tau_m"] / tau_c,
                    "f_lag": (q.get("lagrangian") or {}).get("f_lag", float("nan"))})
    lad.sort(key=lambda d: d["a"])
    la = np.array([d["a"] for d in lad])
    g_A = float(np.mean([d["g"] for d in lad]))
    f_of = PchipInterpolator(la, [d["f"] for d in lad], extrapolate=True)
    m_of = PchipInterpolator(la, [d["m"] for d in lad], extrapolate=True)
    chi_open = PchipInterpolator(la, [d["chi"] for d in lad], extrapolate=True)
    beta_A = float(np.polyfit(np.log(la), np.log([d["h"] for d in lad]), 1)[0])

    print(f"STAGE A ladder at tau_m = {g_A:.2f} tau_c, mc = {mc:.4f}, pmem = {pmem:.6f}")
    print(f"{'a':>5} {'f':>7} {'<m>':>7} {'<m>/f':>6} {'std_m':>7} {'h':>7} {'chi_open':>9} "
          f"{'Phi':>7} {'L_P':>6}")
    for d in lad:
        print(f"{d['a']:5.2f} {d['f']:7.4f} {d['m']:7.4f} {d['m']/d['f']:6.3f} {d['std_m']:7.4f} "
              f"{d['h']:7.4f} {d['chi']:9.4f} {norm.cdf((mc - d['m'])/d['std_m']):7.4f} "
              f"{d['L_P']:6.2f}")
    print(f"  h ~ a^{beta_A:+.3f} at this tau_m;  <m>/f = "
          f"{np.mean([d['m']/d['f'] for d in lad if d['a'] > 0.4]):.3f} averaged over a > 0.4")

    # ============================================================ stage B: the runs
    B = load(a.bparts)
    runs = []
    for q in B:
        u = q["params"]
        w = window_moments(q, tau_c)
        g = u["tau_m"] / tau_c
        ab = r + (1 - r) * (1 - w["chi"])
        runs.append({
            "case": q["_case"], "start": start_of(u), "g": g,
            "chi": w["chi"], "std_chi": w["std_chi"], "m": w["m"], "std_m": w["std_m"],
            "drift": w["drift"], "abar": ab,
            "f_lag": (q.get("lagrangian") or {}).get("f_lag", float("nan")),
            "h": w["std_m"] / math.sqrt(max(w["m"] * (1 - w["m"]), 1e-12)),
            "phi": float(norm.cdf((mc - w["m"]) / w["std_m"])) if w["std_m"] > 0 else float("nan"),
            "beta": beta_cdf(mc, w["m"], w["std_m"]) if w["std_m"] > 0 else float("nan"),
            "f_open": float(f_of(ab)), "m_open": float(m_of(ab)),
        })
    runs.sort(key=lambda d: (d["g"], d["start"]))

    # ---- closure test, all 86 runs
    res_B = np.array([d["chi"] - d["phi"] for d in runs if np.isfinite(d["phi"])])
    res_A = np.array([d["chi"] - norm.cdf((mc - d["m"]) / d["std_m"]) for d in lad])
    print(f"\nCLOSURE chibar = Phi((mc - <m>)/std(m)) with MEASURED moments:")
    print(f"  80 closed-loop runs: rms residual {np.sqrt(np.mean(res_B**2)):.4f}, "
          f"max |res| {np.abs(res_B).max():.4f};  6 open-loop rungs: rms "
          f"{np.sqrt(np.mean(res_A**2)):.4f}, max {np.abs(res_A).max():.4f}")
    mixed = [d for d in runs if 0.05 < d["chi"] < 0.95]
    print(f"  on the {len(mixed)} runs with 0.05 < chibar < 0.95: rms "
          f"{np.sqrt(np.mean([(d['chi']-d['phi'])**2 for d in mixed])):.4f}, "
          f"mean signed {np.mean([d['chi']-d['phi'] for d in mixed]):+.4f}  "
          f"(positive = more nodes below mc than a normal m allows)")
    resb = np.array([d["chi"] - d["beta"] for d in runs if np.isfinite(d["beta"])])
    resbA = np.array([d["chi"] - beta_cdf(mc, d["m"], d["std_m"]) for d in lad])
    print(f"  THE SAME WITH A BETA LAW (same two moments, no extra parameter): 80 runs rms "
          f"{np.sqrt(np.mean(resb**2)):.4f}, max {np.abs(resb).max():.4f};  mixed rms "
          f"{np.sqrt(np.mean([(d['chi']-d['beta'])**2 for d in mixed])):.4f}, mean signed "
          f"{np.mean([d['chi']-d['beta'] for d in mixed]):+.4f};  6 rungs rms {np.sqrt(np.mean(resbA**2)):.4f}")
    print(f"{'g':>7} {'start':>6} {'chibar':>7} {'<m>':>7} {'std_m':>7} {'Phi':>7} {'Beta':>7}")
    for d in runs:
        if d["start"] == "chi0" or (d["start"] == "chi1" and d["chi"] > 0.5):
            print(f"{d['g']:7.3f} {d['start']:>6} {d['chi']:7.4f} {d['m']:7.4f} {d['std_m']:7.4f} "
                  f"{d['phi']:7.4f} {d['beta']:7.4f}")

    # ---- <m> against f, and the deficit against tau_m sigma_m
    print(f"\n<m> AGAINST f.  identity: <m> = <g> + tau_m <m div u>;  deficit = f_lag - <m>")
    print(f"{'g':>7} {'start':>9} {'chibar':>7} {'<m>':>7} {'f_lag':>7} {'deficit':>8} "
          f"{'std_m':>7} {'tau_m*std_m':>11} {'deficit/(tau_m std_m)':>22}")
    defic = []
    for d in runs:
        if not np.isfinite(d["f_lag"]) or d["chi"] > 0.05 and d["chi"] < 0.95:
            continue                       # saturated states only: f_lag is one phase's f
        tm = d["g"] * tau_c
        dd = d["f_lag"] - d["m"]
        defic.append((d["g"], d["start"], dd, tm * d["std_m"], d["chi"]))
        if d["start"] in ("chi0", "chi1"):
            print(f"{d['g']:7.3f} {d['start']:>9} {d['chi']:7.4f} {d['m']:7.4f} {d['f_lag']:7.4f} "
                  f"{dd:8.4f} {d['std_m']:7.4f} {tm*d['std_m']:11.1f} {dd/(tm*d['std_m']):22.3e}")
    dv = np.array([(x[2], x[3]) for x in defic if x[4] < 0.05])
    if len(dv) > 3:
        k = float(np.sum(dv[:, 0] * dv[:, 1]) / np.sum(dv[:, 1] ** 2))
        sc = float(np.std(dv[:, 0] - k * dv[:, 1]))
        print(f"  active phase: deficit = {k:.3e} x (tau_m std_m), scatter {sc:.4f}  "
              f"-> <m div u> ~ {k:.1e} per unit std_m, per step")

    # ============================================================ sigma_m(tau_m, a)
    # a = 1 anchors in (g, h): old A2b at 1 and 3, new ladder, chi0 branch where saturated.
    old_pts = old.get("sigma_m", {}).get("points", [])
    anchors = [(p["g"], p["h"], "L800 A2b") for p in old_pts
               if abs(p["r"] - 1.0) < 1e-6 and p["g"] < 5]
    anchors.append((g_A, [d for d in lad if abs(d["a"] - 1) < 1e-6][0]["h"], "s3 ladder"))
    for d in runs:
        if d["start"] == "chi0" and d["chi"] < 0.02 and d["g"] > 8:
            anchors.append((d["g"], d["h"], "s3 chi0 branch"))
    anchors.sort()
    ga = np.array([x[0] for x in anchors]); ha = np.array([x[1] for x in anchors])
    # merge near-duplicate g (the ladder and the 9.68/11.24 branch points) by averaging log h
    ug, uh = [], []
    for gg, hh in zip(ga, ha):
        if ug and abs(math.log(gg / ug[-1])) < 0.05:
            uh[-1] = 0.5 * (uh[-1] + math.log(hh))
        else:
            ug.append(gg); uh.append(math.log(hh))
    h1 = PchipInterpolator(np.log(ug), uh, extrapolate=True)
    # beta(g): old A2b gives h at a = 1 and 0.45 at four tau_m; the new ladder gives it at 10
    bg, bb = [], []
    for gg in sorted(set(round(p["g"], 2) for p in old_pts)):
        pts = {round(p["r"], 2): p["h"] for p in old_pts if abs(p["g"] - gg) < 0.05}
        if 1.0 in pts and 0.45 in pts:
            bg.append(gg); bb.append(math.log(pts[0.45] / pts[1.0]) / math.log(0.45))
    bg.append(g_A); bb.append(beta_A)
    o = np.argsort(bg)
    beta_of = PchipInterpolator(np.log(np.array(bg)[o]), np.array(bb)[o], extrapolate=True)

    def h_model(g, aa):
        return math.exp(float(h1(math.log(g)))) * aa ** float(beta_of(math.log(g)))

    def sigma_model(g, aa, mean_of):
        mm = float(mean_of(aa))
        return h_model(g, aa) * math.sqrt(max(mm * (1 - mm), 1e-12))

    print(f"\nsigma_m MODEL  h(tau_m, a) = h1(tau_m) a^beta(tau_m)")
    print("  a = 1 anchors: " + "  ".join(f"{g:.2f}->{h:.4f}[{s}]" for g, h, s in anchors))
    print("  beta(tau_m):   " + "  ".join(f"{g:.2f}->{b:+.3f}" for g, b in zip(np.array(bg)[o], np.array(bb)[o])))
    # how the model does on the closed-loop states it did NOT see: the mixed ones and chi1
    print(f"  check on closed-loop runs (chi0/chi1 only):  measured std_m / model at (g, abar)")
    print(f"{'g':>7} {'start':>6} {'chibar':>7} {'abar':>6} {'std_m':>7} {'model':>7} {'ratio':>6}")
    for d in runs:
        if d["start"] not in ("chi0", "chi1"):
            continue
        sm = sigma_model(d["g"], d["abar"], m_of)
        d["sigma_model"] = sm
        if d["start"] == "chi0" or d["chi"] > 0.5:
            print(f"{d['g']:7.3f} {d['start']:>6} {d['chi']:7.4f} {d['abar']:6.3f} {d['std_m']:7.4f} "
                  f"{sm:7.4f} {d['std_m']/sm:6.3f}")

    # ============================================================ the S-curve roots
    variants = {"A: normal, <m> = f(a)": (f_of, "normal"),
                "B: normal, <m> = m_open(a)": (m_of, "normal"),
                "C: Beta, <m> = m_open(a)": (m_of, "beta")}
    branches = {}
    for name, (mean_of, law) in variants.items():
        st_all, un_all, co = [], [], []
        for g in G_GRID:
            def F(cb, g=g, mean_of=mean_of, law=law):
                ab = r + (1 - r) * (1 - cb)
                s = sigma_model(g, ab, mean_of)
                mm = float(mean_of(ab))
                return (float(norm.cdf((mc - mm) / s)) if law == "normal"
                        else beta_cdf(mc, mm, s))
            st, un = roots_of(F)
            st_all.append(st); un_all.append(un); co.append(coexist(st))
        branches[name] = {"stable": st_all, "unstable": un_all, "coexist": co,
                          "tau_freeze": first_g(G_GRID, co)}
        print(f"\nVARIANT {name}:  coexistence from tau_m = {branches[name]['tau_freeze']:.2f} tau_c")
        for gq in (0.3, 1.0, 3.0, 5.0, 8.0, 10.0, 11.24, 12.8, 15.0, 20.0, 30.0):
            i = int(np.argmin(np.abs(np.log(G_GRID / gq))))
            print(f"  g = {G_GRID[i]:6.2f}  stable {', '.join(f'{v:.4f}' for v in st_all[i]):<28} "
                  f"unstable {', '.join(f'{v:.4f}' for v in un_all[i])}")

    # the direct map at tau_m = 10: F = chi_open(abar), nothing assumed
    def F_direct(cb):
        return float(np.clip(chi_open(r + (1 - r) * (1 - cb)), 0.0, 1.0))
    st_d, un_d = roots_of(F_direct)
    x = np.linspace(0, 1, 2001)
    Fd = np.array([F_direct(v) for v in x])
    slope_max = float(np.max(np.gradient(Fd, x)))
    gap = float(np.min(Fd[x > 0.5] - x[x > 0.5]))
    print(f"\nDIRECT MAP at tau_m = {g_A:.2f} tau_c, F(chibar) = chi_open(abar):")
    print(f"  stable roots {st_d}, unstable {un_d};  max slope {slope_max:.2f};  "
          f"min (F - chibar) on chibar > 0.5 = {gap:+.4f}  (a negative gap that small is a "
          f"passive branch about to be born)")

    # ============================================================ measured summary
    meas = {s: sorted([d for d in runs if d["start"] == s], key=lambda d: d["g"])
            for s in ("chi0", "chi1", "leftright", "patches")}
    gm = np.array([d["g"] for d in meas["chi0"]])
    sep = np.array([abs(b["chi"] - a_["chi"]) for a_, b in zip(meas["chi0"], meas["chi1"])])
    tf_meas = first_g(gm, sep > 0.5)
    print(f"\nMEASURED tau_freeze = {tf_meas:.2f} tau_c (separation > 0.5 between chi0 and chi1)")

    # ============================================================ figures
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # ---- fig T1: the closure
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.4))
    ax = axs[0]
    cols = {"chi0": "C3", "chi1": "C0", "leftright": "C2", "patches": "C1"}
    for s, ds in meas.items():
        ax.scatter([d["phi"] for d in ds], [d["chi"] for d in ds], s=22, color=cols[s],
                   label=s, alpha=0.85)
    ax.scatter([norm.cdf((mc - d["m"]) / d["std_m"]) for d in lad], [d["chi"] for d in lad],
               s=60, marker="*", color="k", label=f"open loop, {g_A:.0f} tau_c")
    ax.scatter([d["beta"] for d in runs], [d["chi"] for d in runs], s=40, facecolors="none",
               edgecolors="0.4", lw=0.8, label="same moments, Beta law")
    ax.plot([0, 1], [0, 1], "k--", lw=0.8)
    ax.set_xlabel(r"$P(m<m_c)$ from the run's own two moments (filled: normal, hollow: Beta)")
    ax.set_ylabel(r"measured $\langle\chi\rangle$")
    ax.set_title("the normal closure, per run", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axs[1]
    for s in ("chi0", "chi1"):
        ds = meas[s]
        ax.plot([d["g"] for d in ds], [d["chi"] - d["phi"] for d in ds], "o-", ms=4,
                color=cols[s], label=f"{s}, normal")
        ax.plot([d["g"] for d in ds], [d["chi"] - d["beta"] for d in ds], "o--", ms=4,
                mfc="none", color=cols[s], label=f"{s}, Beta")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xscale("log")
    ax.set_xlabel(r"$\tau_m/\tau_c$")
    ax.set_ylabel(r"$\langle\chi\rangle - \Phi$")
    ax.set_title("closure residual along the scan", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axs[2]
    for s_, col, lab in (("chi0", "C3", "active branch"), ("chi1", "C0", "passive branch")):
        ds = [d for d in meas[s_] if np.isfinite(d["f_lag"]) and (d["chi"] < 0.05 or d["chi"] > 0.95)]
        ax.plot([d["g"] for d in ds], [(d["f_lag"] - d["m"]) / d["f_lag"] for d in ds], "o",
                ms=4, color=col, label=lab)
    ax.plot([d["g"] for d in lad], [(d["f"] - d["m"]) / d["f"] for d in lad], "*", ms=9,
            color="k", ls="none", label="open loop, 6 rungs")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\tau_m/\tau_c$")
    ax.set_ylabel(r"$(f - \langle m\rangle)\,/\,f$")
    ax.set_title(r"the memory-mean deficit, relative to the duty cycle", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "T1_closure.png"), dpi=140)
    plt.close(fig)

    # ---- fig T2: the S-curves at chosen tau_m, with the measured fixed points
    show = [3.426, 9.679, 12.805, 20.621]
    fig, axs = plt.subplots(1, len(show), figsize=(4.0 * len(show), 4.2))
    x = np.linspace(0, 1, 801)
    for ax, gq in zip(axs, show):
        for (name, (mean_of, law)), col, ls in zip(variants.items(), ("C3", "C1", "C4"), ("-", "--", "-")):
            Fx = []
            for v in x:
                ab = r + (1-r)*(1-v); mm = float(mean_of(ab)); sg = sigma_model(gq, ab, mean_of)
                Fx.append(norm.cdf((mc - mm)/sg) if law == "normal" else beta_cdf(mc, mm, sg))
            ax.plot(x, Fx, ls, color=col, lw=1.6, label=name)
        if abs(gq - g_A) < 1.0:
            ax.plot(x, [F_direct(v) for v in x], "-", color="k", lw=1.2,
                    label="direct: chi_open(abar)")
        ax.plot([0, 1], [0, 1], color="0.4", lw=0.8)
        for s in ("chi0", "chi1"):
            d = min(meas[s], key=lambda d: abs(d["g"] - gq))
            ax.plot([d["chi"]], [d["chi"]], "s" if s == "chi1" else "o", color=cols[s], ms=9,
                    mfc="none", mew=2, label=f"measured {s}")
        ax.set_title(rf"$\tau_m/\tau_c = {gq:g}$", fontsize=10)
        ax.set_xlabel(r"$\bar\chi$"); ax.set_ylabel(r"$F(\bar\chi)$")
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.grid(alpha=0.3)
    axs[0].legend(fontsize=7)
    fig.suptitle("the S-curve F and the diagonal: crossings from above are the phases", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "T2_scurves.png"), dpi=140)
    plt.close(fig)

    # ---- fig T3: fig 1 with the predicted branches
    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    for s in ("chi0", "chi1"):
        ds = meas[s]
        ax.errorbar([d["g"] for d in ds], [d["chi"] for d in ds],
                    yerr=[d["_e"] if "_e" in d else 0 for d in ds],
                    fmt="s-" if s == "chi1" else "o-", color=cols[s], ms=5, lw=1.3,
                    label=f"measured {s}")
    for s in ("leftright", "patches"):
        ds = meas[s]
        ax.plot([d["g"] for d in ds], [d["chi"] for d in ds], "^" if s == "patches" else "v",
                color=cols[s], ms=5, ls="none", alpha=0.7, label=f"measured {s}")
    for (name, br), col in zip(branches.items(), ("C3", "C1", "C4")):
        # Roots as POINTS. Above the saddle-node there are two stable branches and one
        # unstable; joining them in g-order draws a zig-zag between 0 and 1, so each root is
        # its own marker and the branches read as dotted curves.
        for kind, mk, ms in (("stable", "o", 2.6), ("unstable", "x", 2.2)):
            gx = [g for g, vs in zip(G_GRID, br[kind]) for _ in vs]
            vy = [v for vs in br[kind] for v in vs]
            if gx:
                ax.plot(gx, vy, mk, color=col, ms=ms, alpha=0.75, ls="none",
                        label=f"theory {name} ({kind})")
        if np.isfinite(br["tau_freeze"]):
            ax.axvline(br["tau_freeze"], color=col, lw=1.0, ls="--", alpha=0.7)
    ax.axvline(tf_meas, color="k", lw=1.0, ls="--", alpha=0.8,
               label=rf"measured $\tau_{{freeze}}$ = {tf_meas:.1f}")
    ax.plot([g_A] * len(st_d), st_d, "*", color="k", ms=13, label="direct map roots (10 tau_c)")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\tau_m/\tau_c$"); ax.set_ylabel(r"$\langle\chi\rangle$")
    ax.set_ylim(-0.03, 1.03)
    ax.set_title("the S-curve prediction against the scan  (dots = stable roots, crosses = unstable)",
                 fontsize=10)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7, loc="center left", ncol=1)
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "T3_prediction.png"), dpi=150)
    plt.close(fig)

    # ---- fig T4: how binary is chi; sigma_m along the branches
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axs[0]
    for s, ds in meas.items():
        ax.plot([d["g"] for d in ds],
                [d["std_chi"]**2 / max(d["chi"] * (1 - d["chi"]), 1e-9) for d in ds],
                "o", ms=4, color=cols[s], label=s)
    ax.set_xscale("log"); ax.set_ylim(0, 1.05)
    ax.set_xlabel(r"$\tau_m/\tau_c$"); ax.set_ylabel(r"Var$_x(\chi)\,/\,\bar\chi(1-\bar\chi)$")
    ax.set_title("how binary the phenotype field is (1 = pure 0/1)", fontsize=10)
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    ax = axs[1]
    for s in ("chi0", "chi1"):
        ds = meas[s]
        ax.plot([d["g"] for d in ds], [d["std_m"] for d in ds], "o-", ms=4, color=cols[s],
                label=f"std(m) measured, {s}")
        ax.plot([d["g"] for d in ds], [d.get("sigma_model", np.nan) for d in ds], "--",
                color=cols[s], alpha=0.7, label=f"model at (tau_m, abar), {s}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$\tau_m/\tau_c$"); ax.set_ylabel(r"$\sigma_m$")
    ax.set_title(r"$\sigma_m$ along the two branches", fontsize=10)
    ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "T4_binary_sigma.png"), dpi=140)
    plt.close(fig)

    # ---- fig T5: the front speed of the leftright start, and where it would vanish
    # Two flat fronts in a periodic box: d<chi>/dt = -2 v / L, v > 0 when the ACTIVE phase
    # advances. Fitted where <chi> runs from 0.45 to 0.25, i.e. after the interface has
    # roughened and before the passive domain breaks up. This is the observable tau_x is
    # made of: the low phase "holds a domain" only where v changes sign.
    fronts = []
    for d in meas["leftright"]:
        sp = os.path.join(a.bparts, d["case"], "series.txt")
        if not os.path.exists(sp):
            continue
        ser = np.loadtxt(sp)
        t, c = ser[:, 0], ser[:, 1]
        if not (c < 0.45).any():
            fronts.append({"g": d["g"], "v": float("nan"), "chi_end": float(c[-1])})
            continue
        i0 = int(np.argmax(c < 0.45))
        i1 = int(np.argmax(c < 0.25)) if (c < 0.25).any() else len(c) - 1
        if i1 - i0 < 20:
            i1 = min(len(c) - 1, i0 + 200)
        A = np.vstack([t[i0:i1], np.ones(i1 - i0)]).T
        sl, ic = np.linalg.lstsq(A, c[i0:i1], rcond=None)[0]
        pred = A @ [sl, ic]
        r2 = 1 - np.sum((c[i0:i1] - pred) ** 2) / max(np.sum((c[i0:i1] - c[i0:i1].mean()) ** 2), 1e-12)
        L = 500.0
        fronts.append({"g": d["g"], "v": float(-L / 2 * sl), "r2": float(r2),
                       "t_window": [float(t[i0]), float(t[i1])], "chi_end": float(c[-1])})
    fr = [f for f in fronts if np.isfinite(f.get("v", np.nan)) and f.get("r2", 0) > 0.9 and f["g"] > 6]
    g_maxwell = float("nan")
    if len(fr) >= 4:
        gx = np.array([f["g"] for f in fr]); vy = np.array([f["v"] for f in fr])
        sel = gx >= 11
        pl = np.polyfit(gx[sel], vy[sel], 1)
        g_maxwell = float(-pl[1] / pl[0]) if pl[0] < 0 else float("nan")
        print(f"\nFRONT SPEED (leftright), v > 0 = active advances:")
        print("  " + "  ".join(f"{f['g']:.1f}->{f['v']:.3f}" for f in fr))
        print(f"  linear fit over tau_m >= 11: v = {pl[1]:.3f} {pl[0]:+.4f} g  -> v = 0 at "
              f"tau_m = {g_maxwell:.1f} tau_c  (the Maxwell point, extrapolated)")
        fig, ax = plt.subplots(figsize=(6.4, 4.2))
        ax.plot(gx, vy, "v", color="C2", ms=7, label="measured, leftright")
        xx = np.linspace(8, max(40, g_maxwell * 1.05 if np.isfinite(g_maxwell) else 40), 50)
        ax.plot(xx, np.polyval(pl, xx), "--", color="0.4", lw=1,
                label=f"linear fit, zero at {g_maxwell:.1f}")
        ax.axhline(0, color="k", lw=0.8)
        ax.axvline(tf_meas, color="k", ls=":", lw=1, label=rf"$\tau_{{freeze}}$ = {tf_meas:.1f}")
        for (name, br), col in zip(branches.items(), ("C3", "C1", "C4")):
            if np.isfinite(br["tau_freeze"]):
                ax.axvline(br["tau_freeze"], color=col, ls="--", lw=0.8, alpha=0.7,
                           label=f"saddle-node, {name.split(':')[0]}")
        ax.set_xlabel(r"$\tau_m/\tau_c$"); ax.set_ylabel(r"front speed $v$  [lattice units / $\tau_c$]")
        ax.set_title("the active phase invades the passive one, slower with memory", fontsize=10)
        ax.set_ylim(bottom=min(-0.2, float(np.nanmin(vy)) - 0.1))
        ax.grid(alpha=0.3); ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(a.out, "T5_fronts.png"), dpi=140)
        plt.close(fig)

    out = {
        "fronts": fronts, "tau_maxwell_extrapolated": g_maxwell,
        "tau_c": tau_c, "mc": mc, "pmem": pmem, "r": r,
        "ladder": lad, "beta_A": beta_A,
        "closure": {"rms_B": float(np.sqrt(np.mean(res_B**2))), "rms_A": float(np.sqrt(np.mean(res_A**2)))},
        "sigma_m_model": {"anchors": anchors, "beta_g": list(map(float, np.array(bg)[o])),
                          "beta": list(map(float, np.array(bb)[o]))},
        "tau_freeze_measured": tf_meas,
        "tau_freeze_theory": {k: v["tau_freeze"] for k, v in branches.items()},
        "direct_map_10": {"stable": st_d, "unstable": un_d, "slope_max": slope_max, "gap": gap},
        "runs": [{k: v for k, v in d.items() if not k.startswith("_")} for d in runs],
    }
    with open(os.path.join(a.out, "theory.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(f"\nwrote {a.out}/T1_closure.png T2_scurves.png T3_prediction.png T4_binary_sigma.png T5_fronts.png theory.json")


if __name__ == "__main__":
    main()
