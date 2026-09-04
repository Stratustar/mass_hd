#!/usr/bin/env python3
"""Runcard generator for the s = 1/3 RESCALED step-switch campaign, L = 500.

    cw_s3_gen.py a [--out cases/20260904]

WHY A RESCALING. The 20260901 campaign ran at L = 800 with Ma = 0.032 at full activity, and
the box-scale acoustic mode (period L/c_s = 1386 steps) sat at Omega_1 tau_c = 2 pi c_s
tau_c / L = 1.02 -- INSIDE the band the memory integrates over, so a sound wave bouncing
round the periodic box was part of the pressure the memory recorded. Dividing every stress
by 3 divides the velocity by 3 (Stokes balance at a fixed viscosity) and multiplies every
time by 3 (tau_c = l/u), while the lattice sound speed stays 1/sqrt(3): Ma falls to 0.011,
and with L = 500 Omega_1 tau_c rises to 4.9. The physics is unchanged up to two things the
lattice Boltzmann does not scale away -- compressibility, which is the point, and inertia:
tau and hence nu = (tau - 1/2)/3 are held, so Re = u l/nu falls 3x, from ~0.8 to ~0.3.
Stage A exists to check that both are small: the five numbers must come back within 5% of
the old ladder once the old ladder is put through the same factors.

THE RULE. Every coefficient in the runcard carries a CLASS, and the class fixes its factor
under stress -> s * stress at fixed lengths (the table is CLASS below; the runcard header
of every case repeats it):

    stress        x s     CC, LL (energy density x length^2), zeta, zeta-open, pmem,
                          video-p-scale
    velocity      x s     video-u-scale
    diffusivity   x s     Dbio (length^2 / time, and time goes as 1/s)
    time          x 1/s   tau-m, tau-chi, and every step count: nsteps, ninfo, nvideo,
                          ntracer, mem-freeze-steps, chi-freeze-steps
    rate          x 1     Gamma (1/(stress x time)); tau, the BGK time, whose fixed nu is
                          what makes u scale as s; friction (0 anyway)
    density       x 1     rho -- the LB density; delta_n = 3 delta_P shrinks with the stress,
                          so the rho = 40 floor of 2026-08-30 is only safer here
    dimensionless x 1     xi, zeta0-frac, mc, m0, chi0, noise, initial-order, switch-sign,
                          chi-lo/hi, m-lo/hi
    geometry      x 1     LX, LY, video-stride, tracer-count, chi-length, chi-block, defect-sep
    step          --      chi-width = 0, pmem-width = 0 (hard steps, no scale)

Invariants worth checking by eye: xi_N = sqrt(LL/2CC) = 2.0, A = zeta_eff/CC = 32 and
l_a = sqrt(LL/zeta_eff) are all unchanged, and sigma_P/zeta_eff should come back at 0.522.

NOT PRESENT in confluent-wet, so nothing to classify: GammaP and the crowding modulus B
(no phase field), and a separate diffusivity for m -- there is ONE Dbio, applied to chi and
m alike, because a cell carries both.

pmem IS NOT SIMPLY DIVIDED BY 3. It is 0.5 sigma_P(a = 1) by construction and stage A
re-measures sigma_P; the value written here is the old one / 3 as a placeholder, which in
open loop only shapes the recorded m and never the flow. mc keeps 0.2100 -- the symmetric
point of cw_step_sym / fine / tchi / long -- a threshold on m, dimensionless.

THE OLD LADDER (L = 800, 20260901 cw_step_a1, calib_step.json) is carried below as
OLD_RUNGS. It sizes the cadence of each rung (a run should not store 4x more video at the
floor to show the same motion) and it is what the check script compares against.
"""
import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cw_step_gen import (tag, write_case, cadence, tracer_interval,       # noqa: E402
                         RECORD_FRAC, STORE_SIGMA, TRACER_COUNT, VIDEO_TAU_C)

S = 1.0 / 3.0                 # the stress factor
L_NEW = 500
CS = 1.0 / math.sqrt(3.0)     # lattice speed of sound

# ------------------------------------------------ the L = 800 working point, with its class
# (value at L = 800, class). Verbatim from cw_step_gen.FROZEN plus the campaign's A1 choices.
OLD = {
    "model":           ("confluent-wet", "flag"),
    "LX":              (800,     "geometry"),
    "LY":              (800,     "geometry"),
    "bc":              (0,       "flag"),
    "nstart":          (0,       "flag"),
    "Gamma":           (0.05,    "rate"),
    "xi":              (0.4,     "dimensionless"),
    "CC":              (0.00625, "stress"),
    "LL":              (0.05,    "stress"),
    "tau":             (2.0,     "rate"),
    "rho":             (40,      "density"),
    "friction":        (0.0,     "rate"),
    "zeta":            (0.2,     "stress"),
    "zeta0-frac":      (0.3,     "dimensionless"),
    "Dbio":            (0.005,   "diffusivity"),
    "angle":           (0,       "dimensionless"),
    "noise":           (0.05,    "dimensionless"),
    "initial-order":   (1,       "dimensionless"),
    "director-config": ("uniform", "flag"),
    "defect-sep":      (0,       "geometry"),
    "video-stride":    (2,       "geometry"),
    "frame-light":     (1,       "flag"),
    "switch-sign":     (-1,      "dimensionless"),
    "chi-width":       (0.0,     "step"),
    "pmem-width":      (0.0,     "step"),
    "mc":              (0.2100,  "dimensionless"),
    "chi0":            (0.5,     "dimensionless"),
    "tracer-count":    (TRACER_COUNT, "geometry"),
}
FACTOR = {"stress": S, "velocity": S, "diffusivity": S, "time": 1.0 / S,
          "rate": 1.0, "density": 1.0, "dimensionless": 1.0, "geometry": 1.0,
          "step": 1.0, "flag": 1.0}
# Keys that are derived per run rather than carried: their class is still recorded so the
# header can list every coefficient in the file.
DERIVED_CLASS = {
    "nsteps": "time", "ninfo": "time", "nvideo": "time", "ntracer": "time",
    "seed": "flag", "open-loop": "flag", "zeta-open": "stress", "chi-config": "flag",
    "tau-chi": "time", "tau-m": "time", "m0": "dimensionless", "pmem": "stress",
    "video-p-scale": "stress", "video-u-scale": "velocity",
}

# The old ladder, MEASURED (calib_step.json, L = 800): a -> tau_c along the tracers (steps),
# sigma_P, u_rms, L_P, f at pmem = 0.5 sigma_P(1) = 0.052177.
OLD_RUNGS = {
    0.30: (899.197, 0.0366630, 0.0142784, 8.4788, 0.056990),
    0.45: (559.017, 0.0528502, 0.0164419, 6.6139, 0.159131),
    0.60: (393.271, 0.0692654, 0.0175779, 5.6413, 0.234395),
    0.80: (291.823, 0.0877946, 0.0186671, 4.4900, 0.291847),
    1.00: (225.560, 0.1043533, 0.0195399, 3.6746, 0.327204),
}
OLD_TAU_C, OLD_SIGMA_P, OLD_U_RMS = OLD_RUNGS[1.0][:3]
OLD_F_TOP = OLD_RUNGS[1.0][4]

# The new campaign clock and scales, PREDICTED (stage A measures them)
TAU_C = OLD_TAU_C / S                   # 676.7 steps
SIGMA_P = OLD_SIGMA_P * S               # 0.034784
U_RMS = OLD_U_RMS * S                   # 0.006513
PMEM = 0.5 * SIGMA_P                    # 0.017392, placeholder until A re-measures sigma_P

A_RATIOS = [1.0, 0.8, 0.65, 0.6, 0.5, 0.3]
SEEDS = {1.0: 3001, 0.8: 3002, 0.65: 3003, 0.6: 3004, 0.5: 3005, 0.3: 3006}
NSTEPS = 100000                         # 150 tau_c in the campaign clock, every rung
ACOUSTIC_COEFF = 0.05                   # the (1,0) mode of P must stay under this x pmem


def old_at(a):
    """(tau_c, sigma_P, u_rms, L_P, f) of the OLD ladder at activity ratio a.

    Power laws between rungs for the three flow scales (they are power laws in a to a few
    per cent), linear for L_P and f."""
    ks = np.array(sorted(OLD_RUNGS))
    cols = np.array([OLD_RUNGS[k] for k in ks])
    la = math.log(a)
    out = []
    for j in range(3):
        out.append(float(math.exp(np.interp(la, np.log(ks), np.log(cols[:, j])))))
    for j in (3, 4):
        out.append(float(np.interp(a, ks, cols[:, j])))
    return tuple(out)


def scaled_frozen():
    v = {}
    for k, (val, cls) in OLD.items():
        f = FACTOR[cls]
        if isinstance(val, float) and f != 1.0:
            v[k] = round(val * f, 8)
        else:
            v[k] = val
    v["LX"] = v["LY"] = L_NEW
    return v


def class_table():
    rows = []
    for k, (val, cls) in OLD.items():
        rows.append((k, cls, f"x{FACTOR[cls]:g}" if FACTOR[cls] != 1 else "x1", val))
    for k, cls in DERIVED_CLASS.items():
        rows.append((k, cls, f"x{FACTOR[cls]:g}" if FACTOR[cls] != 1 else "x1", "derived"))
    w = max(len(r[0]) for r in rows)
    return "\n".join(f"  {k:<{w}}  {cls:<14} {fac:<9} {val}" for k, cls, fac, val in rows)


HDR_A = """20260904 cw_s3_a -- stage A of the s = 1/3 RESCALED campaign: the activity ladder at
L = 500, open loop, six rungs.

WHY. At L = 800 the memory campaign ran at Ma = 0.032 with the box acoustic mode at
Omega_1 tau_c = 1.02, inside the band m integrates over. Every stress is divided by 3, every
time multiplied by 3, the box cut to 500: Ma -> 0.011, Omega_1 tau_c -> 4.9. This ladder
re-measures the five numbers the closed loop needs and checks that the rescaling is a
similarity transform to within 5% -- the two things the LB does not scale away are its
compressibility (the point) and inertia (Re falls 3x, from ~0.8 to ~0.3, because tau is held).

THE FIVE NUMBERS, per rung: sigma_P, tau_c (1/e time of P along a tracer), L_P, u_rms, and
f(a) = fraction of time above pmem once pmem = 0.5 sigma_P(a = 1) is re-taken. PASS if
sigma_P/zeta, L_P and f agree with the old ladder to 5% and tau_c (steps) is 3x the old.
Then the measured tau_c(1) converts every step count of stage B.

THE ACOUSTIC CHECK is read off the video stream, which is why nvideo is the SAME 135 steps
on every rung here (the old ladder let each rung run its own cadence): the box mode has
period L/c_s = 866 steps, so 135 gives 6.4 samples per period, and the (1,0) Fourier
coefficient of the 250x250 block-averaged P is exact for that mode to cos(pi/L). The
criterion is 2|P_hat(1,0)| < {ac:.2f} pmem = {acthr:.2e} throughout; turbulence alone puts
roughly sigma_P L_P sqrt(2 pi)/L ~ {turb:.1e} into that coefficient, so a spectral LINE at
Omega_1 = 2 pi c_s/L, not the raw maximum, is what says "sound".

EVERY COEFFICIENT IN THIS FILE, ITS CLASS AND ITS FACTOR (value column = L = 800 value):
{table}

THIS CASE: a = zeta_eff/zeta = {a:g}, zeta_eff = {ze:.6g}, A = zeta_eff/CC = {A:.1f}, seed {seed}.
Old ladder at this a: tau_c = {tc_old:.0f} steps, sigma_P = {sp_old:.5f}, L_P = {lp_old:.2f},
f(0.5 sigma_P(1)) = {f_old:.4f}.  Expected here: tau_c = {tc:.0f} steps ({ntc:.0f} of its own
tau_c in the record window), sigma_P = {sp:.5f}, Omega_1 tau_c(a) = {om:.2f}."""


def gen_a(out):
    made = []
    base_v = scaled_frozen()
    for a in A_RATIOS:
        tc_old, sp_old, ur_old, lp_old, f_old = old_at(a)
        tc, sp, ur = tc_old / S, sp_old * S, ur_old * S
        ze = base_v["zeta"] * a
        ninfo, nvideo = cadence(tc, TAU_C, NSTEPS, video_ref=True)
        v = dict(base_v)
        v.update({
            "nsteps": NSTEPS, "ninfo": ninfo, "seed": SEEDS[a],
            "open-loop": 1, "zeta-open": round(ze, 8),
            "chi-config": "uniform",
            "tau-chi": round(0.3 * TAU_C, 1),
            "tau-m": round(10.0 * TAU_C, 1),
            "m0": round(OLD_F_TOP, 4),
            "pmem": round(PMEM, 6),
            "ntracer": tracer_interval(tc), "tracer-count": TRACER_COUNT,
            "nvideo": nvideo,
            "video-p-scale": round(STORE_SIGMA * sp, 5),
            "video-u-scale": round(STORE_SIGMA * ur, 5),
        })
        hdr = HDR_A.format(
            ac=ACOUSTIC_COEFF, acthr=ACOUSTIC_COEFF * PMEM,
            turb=2 * SIGMA_P * lp_old * math.sqrt(2 * math.pi) / L_NEW,
            table=class_table(),
            a=a, ze=ze, A=ze / base_v["CC"], seed=SEEDS[a],
            tc_old=tc_old, sp_old=sp_old, lp_old=lp_old, f_old=f_old,
            tc=tc, ntc=RECORD_FRAC * NSTEPS / tc, sp=sp,
            om=2 * math.pi * CS * tc / L_NEW)
        made.append(write_case(os.path.join(out, "cw_s3_a", f"a{tag(a)}"), hdr, v))
    return made


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["a"])
    ap.add_argument("--out", default="cases/20260904")
    args = ap.parse_args()
    print(f"s = {S:.4f}: tau_c {OLD_TAU_C:.1f} -> {TAU_C:.1f} steps, sigma_P {OLD_SIGMA_P:.5f} -> "
          f"{SIGMA_P:.5f}, u_rms {OLD_U_RMS:.5f} -> {U_RMS:.5f}, pmem -> {PMEM:.6f} (placeholder)")
    print(f"Omega_1 tau_c = {2 * math.pi * CS * TAU_C / L_NEW:.2f},  Ma -> "
          f"{U_RMS / CS:.4f},  0.05 pmem = {ACOUSTIC_COEFF * PMEM:.2e}")
    print("class table:\n" + class_table())
    made = gen_a(args.out)
    for p in made:
        print(f"  {p}")
    print(f"{len(made)} runcards")


if __name__ == "__main__":
    main()
