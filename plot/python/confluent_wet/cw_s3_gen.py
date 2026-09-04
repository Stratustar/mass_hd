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
from cw_step_gen import (tag, write_case, cadence, tracer_interval, ORDER,  # noqa: E402
                         saturation_note, RECORD_FRAC, STORE_SIGMA, TRACER_COUNT,
                         VIDEO_TAU_C)

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


# =========================================================================== B

# The runcard key order, plus init-frame. ORDER itself is left alone so every earlier
# campaign's runcards still round-trip byte for byte.
ORDER_B = list(ORDER)
ORDER_B.insert(ORDER_B.index("chi-freeze-steps"), "init-frame")

# The measured stage-A calibration, filled by Calib3 from calib_s3.json.
TAU_M_OVER_TAU_C = [0.300, 1.863, 3.426, 4.989, 6.553, 8.116, 9.679, 11.242, 12.805,
                    14.368, 15.932, 17.495, 19.058, 20.621, 22.184, 23.747, 25.311,
                    26.874, 28.437, 30.000]
T_OVER_TAU_C = 300.0          # every run, one number
SERIES_TAU_C = 0.05           # the video/series cadence, in tau_c
FRAMES_TARGET = 26            # full frames per run
MEM_FREEZE_TAU_C_B = 1.0      # insurance only: init-frame means the flow is already developed
CHI_LENGTH_OVER_L_P = 2.0
FRAME_ROOT = "/scratch/helu/mass_hd/cases/20260904/cw_s3_a"
# Which open-loop snapshot each start is laid on. chi0 and chi1 each get the flow of the
# phase they ARE; the two mixed starts get the ladder's middle rung, which is the closest
# thing to the flow a half-and-half layer generates -- zeta_eff = 0.65 zeta against the
# 0.65 zeta a 50/50 mixture produces at z0 = 0.3 (r + (1-r)(1-0.5) = 0.65 exactly).
SNAPSHOT = {"chi0": "a1", "chi1": "a0p3", "leftright": "a0p65", "patches": "a0p65"}
SEED0, CHI_SEED0 = 4000, 7000


class Calib3:
    """The MEASURED stage-A calibration; without a file, stage A's own predictions."""

    def __init__(self, path=None):
        self.measured = path is not None
        if path:
            import json
            with open(path) as fh:
                d = json.load(fh)
            self.tau_c = float(d["tau_c"]); self.sigma_P = float(d["sigma_P"])
            self.u_rms = float(d["u_rms"]); self.pmem = float(d["pmem"])
            self.mc = float(d["mc"]); self.f_top = float(d["f_top"])
            self.f_floor = float(d["f_floor"]); self.L_P = float(d["L_P"])
        else:
            self.tau_c, self.sigma_P, self.u_rms = TAU_C, SIGMA_P, U_RMS
            self.pmem = PMEM; self.mc = OLD["mc"][0]
            self.f_top, self.f_floor = OLD_F_TOP, OLD_RUNGS[0.3][4]
            self.L_P = OLD_RUNGS[1.0][3] * 1.0


HDR_B = """20260904 cw_s3_b -- the tau_m scan of the s = 1/3 campaign: 20 memory times x 4 starts.

THE QUESTION. Four initial conditions on one tau_m axis. Where the chi == 0 and chi == 1
curves COINCIDE the outcome is independent of how the layer started -- a mixed state, or a
single high phase; where they SEPARATE is tau_freeze. The mixed starts (leftright, patches)
are not summarised by a number: they are watched, because what they show is whether the low
phase can hold a domain at all, which is tau_x.

THE FOUR STARTS. m is set to MATCH chi everywhere, because a phase is a (chi, m) pair and
a half prepared at the wrong m relaxes towards the other phase on tau_m whether or not the
physics is bistable:
  chi0        chi == 0 everywhere, m == f(zeta)     = {f_top:.4f}   the ACTIVE phase
  chi1        chi == 1 everywhere, m == f(0.3 zeta) = {f_floor:.4f}   the PASSIVE phase
  leftright   x < L/2 is chi = 0, the rest chi = 1; two flat interfaces under PBC
  patches     correlated Gaussian noise at {chilen:.1f} = 2 L_P, split at its OWN MEDIAN into
              exactly half and half -- a fixed threshold would let the sample mean scatter
              the initial area fraction, which is the very quantity the loop amplifies
Both uniform starts are EXACT fixed points of the hard step chi* = Theta(mc - m), with
mc = {mc:.4f} between f_floor and f_top, so "did it move" is a clean question.

A NOTE ON THE hi/lo SUFFIXES. The model applies (chi-hi, m-hi) to the left half of a stripe
and to the above-median half of the binary noise. This campaign puts the ACTIVE phase there,
so chi-hi = 0 and m-hi = f_top: the suffix names the model's slot, not the value of chi.

EVERY RUN STARTS FROM A DEVELOPED FLOW. init-frame loads Q, u and P from the last frame of
the matching open-loop rung of stage A ({snap}), and the LB populations are rebuilt as
f_eq(u_code, n) with n = rho + 3(P + sigma_bulk) and u_code = u_mat - div(sigma)/2n. Only
the non-equilibrium part of ff is lost and it regenerates on tau = 2 steps; measured on a
restart test, the flow tracks the continued run to 0.3% over 300 steps with no jump at
t = 0. This replaces the 10 tau_c memory freeze the L = 800 campaign needed: the memory now
integrates the pressure of the turbulent state from the first step, which is what it is
meant to do. A {freeze:.0f}-step freeze is kept purely as insurance.

THE CLOCKS, all from the MEASURED stage-A calibration (calib_s3.json):
  tau_c   = {tau_c:.1f} steps      the material clock, tracers at full activity
  sigma_P = {sigma_P:.5f}        pmem = 0.5 sigma_P = {pmem:.6f}
  T       = {T:g} tau_c = {nsteps} steps, the same for all 80 runs
  tau_chi = 0.3 tau_c = {tau_chi:.0f} steps
  series  = every {nvideo} steps = {series:.3f} tau_c, the video stream's own clock; <chi> and
            std(chi) are written there in double precision from the full-resolution field

THIS CASE: tau_m = {g:g} tau_c = {tau_m:.0f} steps, start {start}, on the {snap} flow.
seed {seed}, chi-seed {chi_seed}.
{sat}"""


def gen_b(out, cal):
    made = []
    nsteps = int(round(T_OVER_TAU_C * cal.tau_c / 100.0)) * 100
    nvideo = max(1, int(round(SERIES_TAU_C * cal.tau_c)))
    ninfo = int(round(nsteps / FRAMES_TARGET / 500.0)) * 500
    tau_chi = 0.3 * cal.tau_c
    freeze = int(round(MEM_FREEZE_TAU_C_B * cal.tau_c))
    chilen = CHI_LENGTH_OVER_L_P * cal.L_P
    base_v = scaled_frozen()
    base_v.update({
        "nsteps": nsteps, "ninfo": ninfo, "nvideo": nvideo,
        "mc": round(cal.mc, 4), "pmem": round(cal.pmem, 6),
        "tau-chi": round(tau_chi, 1),
        "mem-freeze-steps": freeze,
        "ntracer": tracer_interval(cal.tau_c), "tracer-count": TRACER_COUNT,
        "video-p-scale": round(STORE_SIGMA * cal.sigma_P, 5),
        "video-u-scale": round(STORE_SIGMA * cal.u_rms, 5),
    })
    # (chi, m) of the two phases, and the model slots they go into
    starts = {
        "chi0":      {"chi-config": "uniform", "chi0": 0.0, "m0": round(cal.f_top, 4)},
        "chi1":      {"chi-config": "uniform", "chi0": 1.0, "m0": round(cal.f_floor, 4)},
        "leftright": {"chi-config": "stripe", "chi0": 0.5, "m0": 0.5,
                      "chi-hi": 0.0, "m-hi": round(cal.f_top, 4),
                      "chi-lo": 1.0, "m-lo": round(cal.f_floor, 4)},
        "patches":   {"chi-config": "binary-noise", "chi0": 0.5, "m0": 0.5,
                      "chi-noise": 1.0, "chi-length": round(chilen, 2),
                      "chi-hi": 0.0, "m-hi": round(cal.f_top, 4),
                      "chi-lo": 1.0, "m-lo": round(cal.f_floor, 4)},
    }
    for i, g in enumerate(TAU_M_OVER_TAU_C):
        tau_m = g * cal.tau_c
        for j, (name, over) in enumerate(starts.items()):
            seed = SEED0 + 4 * i + j
            chi_seed = CHI_SEED0 + 4 * i + j
            v = dict(base_v)
            v.update(over)
            v.update({
                "seed": seed, "chi-seed": chi_seed, "tau-m": round(tau_m, 1),
                "init-frame": f"{FRAME_ROOT}/{SNAPSHOT[name]}/frame{NSTEPS}.json",
            })
            hdr = HDR_B.format(
                f_top=cal.f_top, f_floor=cal.f_floor, mc=cal.mc, chilen=chilen,
                snap=SNAPSHOT[name], freeze=freeze, tau_c=cal.tau_c, sigma_P=cal.sigma_P,
                pmem=cal.pmem, T=T_OVER_TAU_C, nsteps=nsteps, tau_chi=tau_chi,
                nvideo=nvideo, series=nvideo / cal.tau_c,
                g=g, tau_m=tau_m, start=name, seed=seed, chi_seed=chi_seed,
                sat=saturation_note(nsteps, cal.tau_c, tau_m, tau_chi))
            made.append(write_case(
                os.path.join(out, "cw_s3_b", f"tm{tag(g)}_{name}"), hdr, v, order=ORDER_B))
    return made


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["a", "b"])
    ap.add_argument("--calib", default=None,
                    help="measured calib_s3.json from stage A; stage B refuses "
                         "to run without it")
    ap.add_argument("--out", default="cases/20260904")
    args = ap.parse_args()
    print(f"s = {S:.4f}: tau_c {OLD_TAU_C:.1f} -> {TAU_C:.1f} steps, sigma_P {OLD_SIGMA_P:.5f} -> "
          f"{SIGMA_P:.5f}, u_rms {OLD_U_RMS:.5f} -> {U_RMS:.5f}, pmem -> {PMEM:.6f} (placeholder)")
    print(f"Omega_1 tau_c = {2 * math.pi * CS * TAU_C / L_NEW:.2f},  Ma -> "
          f"{U_RMS / CS:.4f},  0.05 pmem = {ACOUSTIC_COEFF * PMEM:.2e}")
    print("class table:\n" + class_table())
    if args.stage == "b":
        if not args.calib:
            raise SystemExit("stage B needs --calib: sizing 80 runs on stage A's PREDICTED "
                             "scales would put every step count and every m0 on numbers the "
                             "ladder was run to replace")
        cal = Calib3(args.calib)
        print(f"measured: tau_c = {cal.tau_c:.1f}, sigma_P = {cal.sigma_P:.5f}, "
              f"pmem = {cal.pmem:.6f}, f_top = {cal.f_top:.4f}, f_floor = {cal.f_floor:.4f}, "
              f"L_P = {cal.L_P:.2f}, mc = {cal.mc:.4f}")
        made = gen_b(args.out, cal)
    else:
        made = gen_a(args.out)
    for p in made:
        print(f"  {p}")
    print(f"{len(made)} runcards")


if __name__ == "__main__":
    main()
