#!/usr/bin/env python3
"""Runcard generator for the 2026-09 STEP-SWITCH memory campaign.

    cw_step_gen.py a1|a2b|a4 [--out cases/20260901] [--calib calib.json]

WHAT MAKES THIS CAMPAIGN DIFFERENT from the 20260831 one cw_gen.py serves: both switches
are HARD STEPS, g(P) = Theta(P - pmem) and chi*(m) = Theta(mc - m), so chi-width is gone as
a parameter. What replaces it is not a smaller number but a MECHANISM -- the spatial spread
of m smooths the pointwise step when it is averaged over space, and the space-averaged
chi* is just the CDF of m, whose central slope corresponds to an effective width
sqrt(2 pi)/2 * sigma_m. The switch therefore sharpens by itself as tau_m grows and sigma_m
falls, and the campaign's single knob tau_m sets both the memory time AND the sharpness.

THE THRESHOLD tau_x. Bistability is realised when sigma_m drops below the cusp value

    K_max = max_chibar  psi(chibar) * (1 - r) * |df/da| ,   a(chibar) = r + (1-r)(1-chibar)

with psi the standard normal density at the chibar quantile. sigma_m falls with tau_m, so
this defines a crossing tau_x -- the whole point of the tau_m scan. The closed form

    sigma_m(tau) = sqrt(f_b (1 - f_b)) sqrt(tau_c / (tau_c + tau))

(the variance of a one-pole filter driven by exponentially correlated noise) is verified
to 0.953 at Dbio = 0 and NOT usable at the campaign's Dbio: cell motility smooths m over
l_D = sqrt(2 Dbio tau_m), which passes L_P at tau_m ~ 6 tau_c and keeps growing, taking
sigma_m down to 0.36 and 0.24 of the closed form at tau_m/tau_c = 10 and 30. Stage A2b
measures sigma_m(tau_m) instead of predicting it, and A3 reads tau_x off that measurement.

THE PRE-CALIBRATION BELOW is an EXTRAPOLATION of the wave-1 open-loop ladder
(cases/20260831/cw_lag_base, L = 1000, zeta_eff = 0.05 / 0.1 / 0.2), fitted as

    tau_c(a) = 222.3 a^-1.146      sigma_P(a) = 0.10343 a^0.877
    u_rms(a) = 0.01865 a^0.200     a = zeta_eff / zeta

It sizes stage A, and nothing else: A1 measures these on THIS box (L = 800, not 1000) and
A3 re-derives every closed-loop number from the measurement. Passing --calib replaces the
extrapolation with the measured calib.json, which is what the A4 and B stages must use.
"""
import argparse
import json
import math
import os

# ----------------------------------------------------------------- the frozen working point
FROZEN = {
    "model": "confluent-wet", "LX": 800, "LY": 800, "bc": 0, "nstart": 0,
    "Gamma": 0.05, "xi": 0.4, "CC": 0.00625, "LL": 0.05, "tau": 2.0, "rho": 40,
    "friction": 0.0, "zeta": 0.2, "zeta0-frac": 0.3, "Dbio": 0.005,
    "angle": 0, "noise": 0.05, "initial-order": 1, "director-config": "uniform",
    "defect-sep": 0, "video-stride": 2, "frame-light": 1, "switch-sign": -1,
    # BOTH SWITCHES HARD. Not "very small": 0 selects the step branch in the model.
    "chi-width": 0.0, "pmem-width": 0.0,
}

# ---- wave-1 extrapolation; superseded by calib.json from A1 ----
PRE = {
    "tau_c": 222.3, "tau_c_exp": -1.1464,
    "sigma_P": 0.10343, "sigma_P_exp": 0.8768,
    "u_rms": 0.01865, "u_rms_exp": 0.2004,
    "pmem_coeff": 0.5,        # pmem = 0.5 sigma_P(zeta): the contrast maximum on the A3 scan
    "mc": 0.236,              # the cusp mc for r = 0.3
    "f_top": 0.3085,          # f(zeta)   -> m0 of the chi == 0 arm
    "f_floor": 0.0754,        # f(0.3 z)  -> m0 of the chi == 1 arm
}
R_FLOOR = 0.3                 # the one activity floor this campaign runs
NSTEPS = 100000               # stage A default; A4 and the long B points override it
RECORD_FRAC = 0.8             # the record window is the closing 80%
SATURATION_K = 15.0           # record window >= 15 (tau_c + tau_m + tau_chi)
VIDEO_TAU_C = 0.2             # video frame spacing, in tau_c(zeta)
INFO_TAU_C = 18.0             # full-frame spacing, in tau_c(zeta) -- ~25 frames per run
MEM_FREEZE_TAU_C = 10.0       # hold the memory source off this long while the flow spins up
TRACER_COUNT = 2000
STORE_SIGMA = 6.0             # video stores +/- 6 sigma_P, the renderer displays +/- 3

ORDER = ["model", "nsteps", "nstart", "ninfo", "LX", "LY", "bc", "seed", "",
         "Gamma", "xi", "CC", "LL", "tau", "rho", "friction", "zeta", "zeta0-frac",
         "open-loop", "zeta-open", "",
         "angle", "noise", "initial-order", "director-config", "defect-sep", "",
         "Dbio", "chi-config", "chi0", "chi-noise", "chi-length", "chi-seed",
         "chi-lo", "chi-hi", "m-lo", "m-hi", "chi-block",
         "tau-chi", "chi-width", "switch-sign", "mc", "m0", "chi-freeze-steps",
         "mem-freeze-steps", "",
         "tau-m", "pmem", "pmem-width", "",
         "ntracer", "tracer-count", "",
         "nvideo", "video-stride", "video-p-scale", "video-u-scale", "frame-light"]


def tag(x, nd=2):
    """0.45 -> 0p45, 1.0 -> 1, -0.06 -> m0p06; the project's directory-name convention."""
    s = f"{x:.{nd}f}".rstrip("0").rstrip(".")
    if s in ("", "-0"):
        s = "0"
    return s.replace("-", "m").replace(".", "p")


def write_case(path, header, vals):
    os.makedirs(path, exist_ok=True)
    unknown = set(vals) - set(ORDER)
    if unknown:
        raise RuntimeError(f"keys not in ORDER (they would be silently dropped): {unknown}")
    lines = [f"# {h}" for h in header.strip("\n").split("\n")]
    for key in ORDER:
        if key == "":
            lines.append("")
        elif key in vals:
            v = vals[key]
            v = f"{v:.6e}" if isinstance(v, float) and abs(v) < 1e-3 and v != 0 else v
            lines.append(f"{key:<17}= {v}")
    with open(os.path.join(path, "run.dat"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return os.path.join(path, "run.dat")


class Calib:
    """The campaign's scales, from the wave-1 extrapolation or from a measured calib.json."""

    def __init__(self, path=None):
        self.measured = path is not None
        if path:
            with open(path) as fh:
                d = json.load(fh)
            self.tau_c = float(d["tau_c"])
            self.sigma_P = float(d["sigma_P"])
            self.u_rms = float(d["u_rms"])
            self.pmem = float(d["pmem"])
            self.mc = float(d["mc"])
            self.f_top = float(d["f_top"])
            self.f_floor = float(d["f_floor"])
            self.tau_x = float(d.get("tau_x", float("nan")))
            self.tau_x_co = d.get("tau_x_coexistence")
            self.pmem_on_edge = bool(d.get("pmem_on_scan_edge", False))
            self.f_table = d.get("f_table") or []
            self.sigma_P_meas = float(d["sigma_P"])
            self._rung = d.get("rungs")            # measured per-activity table, if present
            # L_P at full activity: sets the mixed start's correlation length
            top = (self._rung or {}).get("1", (self._rung or {}).get("1.0", {}))
            self.L_P = float(top.get("L_P", float("nan")))
        else:
            self.tau_c = PRE["tau_c"]
            self.sigma_P = PRE["sigma_P"]
            self.u_rms = PRE["u_rms"]
            self.pmem = PRE["pmem_coeff"] * PRE["sigma_P"]
            self.mc = PRE["mc"]
            self.f_top = PRE["f_top"]
            self.f_floor = PRE["f_floor"]
            self.tau_x = float("nan")
            self.tau_x_co = None
            self.pmem_on_edge = False
            self.f_table = []
            self.sigma_P_meas = PRE["sigma_P"]
            self.L_P = float("nan")
            self._rung = None

    def at(self, a):
        """(tau_c, sigma_P, u_rms) at activity ratio a, measured where available."""
        if self._rung and str(a) in self._rung:
            r = self._rung[str(a)]
            return r["tau_c"], r["sigma_P"], r["u_rms"]
        return (self.tau_c * a**PRE["tau_c_exp"],
                self.sigma_P * a**PRE["sigma_P_exp"],
                self.u_rms * a**PRE["u_rms_exp"])


MIN_RECORD_FRAMES = 20            # full frames inside the record window, per run


def cadence(tau_c_local, tau_c_ref, nsteps, video_ref):
    """ninfo and nvideo.

    The full-frame clock follows the run's OWN tau_c -- those frames are averaged, not
    watched, so what matters is that consecutive ones are independent, and a rung whose
    turnover is 4x slower would otherwise contribute 4x fewer independent snapshots. But it
    is capped so the record window still holds MIN_RECORD_FRAMES: at the floor rung
    18 tau_c_local is 16000 steps, which left 5 frames to compute std_R, L_P and L_chi
    from. 4 tau_c apart is still independent; 5 frames is not a statistic.

    The video clock is the campaign's tau_c for CLOSED-LOOP runs, so every production run
    animates at one rate and can be compared by eye. Open-loop ladder rungs use their own,
    because nothing compares them frame to frame and a 4x slower rung would otherwise store
    4x more bytes to show the same motion."""
    ninfo = min(int(round(INFO_TAU_C * tau_c_local / 500)) * 500,
                int(nsteps * RECORD_FRAC / MIN_RECORD_FRAMES / 500) * 500)
    tc_vid = tau_c_ref if video_ref else tau_c_local
    return (max(500, ninfo), max(5, int(round(VIDEO_TAU_C * tc_vid / 5)) * 5))


def tracer_interval(tau_c_local):
    """~26 samples per correlation time, rounded to a multiple of 5."""
    return max(5, int(round(tau_c_local / 26.0 / 5.0)) * 5)


def base(cal, a_for_scales, nsteps=NSTEPS, video_ref=True, **over):
    """The frozen block plus the scales at activity ratio `a_for_scales`."""
    tc_local, sp_local, ur_local = cal.at(a_for_scales)
    ninfo, nvideo = cadence(tc_local, cal.tau_c, nsteps, video_ref)
    v = dict(FROZEN)
    v.update({
        "nsteps": nsteps, "ninfo": ninfo,
        "pmem": round(cal.pmem, 6), "mc": round(cal.mc, 4),
        "nvideo": nvideo,
        "video-p-scale": round(STORE_SIGMA * sp_local, 4),
        "video-u-scale": round(STORE_SIGMA * ur_local, 4),
        "ntracer": tracer_interval(tc_local), "tracer-count": TRACER_COUNT,
    })
    v.update(over)
    return v


def saturation_note(nsteps, tau_c, tau_m, tau_chi):
    need = SATURATION_K * (tau_c + tau_m + tau_chi)
    rec = RECORD_FRAC * nsteps
    ok = rec >= need
    return (f"SATURATION: record window {rec:.0f} steps against "
            f"{SATURATION_K:.0f}(tau_c + tau_m + tau_chi) = {need:.0f} -- "
            f"{'satisfied' if ok else f'SHORT by {need/rec:.1f}x'}.")


# ------------------------------------------------------------------------ A1

A1_RATIOS = [1.0, 0.8, 0.6, 0.45, 0.3]

HDR_A1 = """20260901 cw_step_a1 -- stage A1: the ACTIVITY LADDER, open loop.

Five OPEN-LOOP runs. m and chi evolve as they will in production, but the activity is
pinned at zeta-open and chi is not fed back, so this measures the layer's RESPONSE to a
prescribed activity rather than a fixed point. It produces every number the closed loop
needs and cannot guess:

  tau_c(zeta)              the 1/e decorrelation time of P ALONG A TRACER at the top rung,
                           and from here on THE time unit of the campaign. A material clock,
                           not an Eulerian one: the memory integrates g(P) following a cell.
  sigma_P(zeta), u_rms(zeta)   the fixed colour scales of every production video
  f(pmem; zeta_eff)        with a HARD g, simply the fraction of time a lattice point spends
                           above pmem. Stage A3 solves chibar = Phi((mc - f(a(chibar)))/sigma_m)
                           on it -- the step version of the fixed-point equation, in which
                           the smoothing comes from sigma_m and not from a switch width.

CHECKPOINT 1 is the r = 0.30 rung, the activity FLOOR the passive phenotype will sit at. It
must still be turbulent -- a stable non-zero defect count -- and tau_c(0.3 zeta)/tau_c(zeta)
must stay under 5. The wave-1 extrapolation predicts 3.98. If either fails the floor rises
to 0.4 zeta and the ladder is re-run: a passive phase that has stopped stirring generates no
pressure, can never be re-excited, and is an absorbing state rather than a phase.

THE SECOND FLOOR WAS CANCELLED. This campaign was designed with two floors, 0.3 and 0.6
zeta, to test the 1/(1-r)^2 scaling of tau_x. The r = 0.6 half is dropped: its tau_x lands
near 100 tau_c, outside any tau_m the run length can reach, so the whole group would have
returned "mixed" at every point. The ladder still spans a = 1.0 down to 0.3 because A3 needs
df/da across the REACHABLE activity range, which is exactly [r, 1] = [0.3, 1].

THIS CASE: r = zeta_eff/zeta = {r:g}, zeta_eff = {ze:g}, A = zeta_eff/CC = {A:.1f}.
Predicted tau_c = {tc:.0f} steps ({ntc:.0f} tau_c of record), sigma_P = {sp:.4f}.

VIDEO SCALES are PER RUNG here and predicted, because the ladder spans a factor 3 in
sigma_P and one scale would waste most of the byte at the bottom. From A4 on they are the
single measured pair, so that production runs stay comparable by eye."""


def gen_a1(out, cal):
    made = []
    for r in A1_RATIOS:
        ze = FROZEN["zeta"] * r
        tc, sp, ur = cal.at(r)
        v = base(cal, r, video_ref=False,
                 seed=1001, **{"open-loop": 1, "zeta-open": round(ze, 6),
                               "chi-config": "uniform", "chi0": 0.5,
                               "tau-chi": round(0.3 * cal.tau_c, 1),
                               "m0": round(cal.f_top, 4),
                               "tau-m": round(10.0 * cal.tau_c, 1)})
        hdr = HDR_A1.format(r=r, ze=ze, A=ze / FROZEN["CC"], tc=tc, sp=sp,
                            ntc=RECORD_FRAC * NSTEPS / tc)
        made.append(write_case(os.path.join(out, "cw_step_a1", f"r{tag(r)}"), hdr, v))
    return made


# ----------------------------------------------------------------------- A2b

A2B_TAU_M = [1.0, 3.0, 30.0]      # tau_m = 10 tau_c is already in A1, at both these rungs
A2B_RATIOS = [1.0, 0.45]

HDR_A2B = """20260901 cw_step_a2b -- stage A2b: sigma_m(tau_m), MEASURED.

The realisation threshold tau_x is defined by sigma_m(tau_x) = K_max, so the campaign's
central prediction is only as good as sigma_m(tau_m). The closed form

    sigma_m(tau) = sqrt(f(1-f)) sqrt(tau_c/(tau_c+tau))

is the variance of a one-pole filter driven by exponentially correlated noise, and it is
RIGHT -- at Dbio = 0 a local test measured 0.953 of it. At the campaign's Dbio = 0.005 it is
not: measured/closed came out 0.72 / 0.54 / 0.36 / 0.24 at tau_m/tau_c = 1 / 3 / 10 / 30,
because cell motility smooths m over l_D = sqrt(2 Dbio tau_m), which passes L_P at
tau_m ~ 6 tau_c and keeps growing. The measured exponent is -0.73 against the filter's -0.5.

Dbio is a PHYSICAL parameter here, not a numerical smoother, so the answer is to measure
sigma_m rather than to remove the diffusion. Open loop pins f, which is what makes the
tau_m dependence readable at all: in a closed-loop run f moves with chi and the two effects
are inseparable.

TWO RUNGS, because sigma_m depends on f through sqrt(f(1-f)) and the cusp sits near
a = 0.51, between them. a = 1.0 and 0.45 bracket it; A1 supplies the tau_m = 10 tau_c point
at both, so this stage runs only the other three.

THIS CASE: r = {r:g}, zeta_eff = {ze:g}, tau_m = {g:g} tau_c = {tm:.0f} steps.
{sat}"""


def gen_a2b(out, cal):
    made = []
    for r in A2B_RATIOS:
        ze = FROZEN["zeta"] * r
        for g in A2B_TAU_M:
            tm = g * cal.tau_c
            v = base(cal, r, video_ref=False,
                     seed=1001, **{"open-loop": 1, "zeta-open": round(ze, 6),
                                   "chi-config": "uniform", "chi0": 0.5,
                                   "tau-chi": round(0.3 * cal.tau_c, 1),
                                   "m0": round(cal.f_top, 4), "tau-m": round(tm, 1)})
            tc = cal.at(r)[0]
            hdr = HDR_A2B.format(r=r, ze=ze, g=g, tm=tm,
                                 sat=saturation_note(NSTEPS, tc, tm, 0.3 * cal.tau_c))
            made.append(write_case(
                os.path.join(out, "cw_step_a2b", f"r{tag(r)}_tm{tag(g)}"), hdr, v))
    return made


# ------------------------------------------------------------------------ A4

HDR_A4 = """20260901 cw_step_a4 -- stage A4: does the hard switch separate at all?

Two CLOSED-LOOP runs at tau_m = 30 tau_c, deep on the bistable side of the predicted
tau_x ~ 2.5 tau_c, started from the two ends:

  chi == 0, m == f(zeta)     the fully ACTIVE phase, at its own fixed point
  chi == 1, m == f(0.3 zeta) the PASSIVE phase, at its own

Both are EXACT fixed points of the hard step: f(zeta) = {ftop:.4f} > mc = {mc:.4f} gives
chi* = Theta(mc - m) = 0, and f(0.3 zeta) = {ffl:.4f} < mc gives chi* = 1. So the arms hold
unless the dynamics moves them, and the gate is whether sigma_m is small enough for the two
basins to exist.

REQUIRED: the two arms end more than 0.5 apart in <chi>. If they converge, tau_x is larger
than the closed form and the measured sigma_m say, and the B grid has to be re-centred
before 100 runs are spent on it.

WATCH THE APPROACH, NOT ONLY THE END. From rest P = 0, so g = 0 and m decays on tau_m
before there is any turbulence to feed it; crossing mc from f(zeta) takes
ln(f_top/mc) tau_m = {cross:.2f} tau_m = {crossteps:.0f} steps against a flow that reaches
u_rms within ~2 tau_c. The margin is real but it is not large, and it is SMALLER than in the
tanh campaign, where mc sat further below f_top. A local L = 200 test never crossed. The
`departure_steps` field in part.json is what checks it per run.

THIS CASE: {arm}.
{sat}"""


def gen_a4(out, cal):
    made = []
    tm = 30.0 * cal.tau_c
    tchi = 0.3 * cal.tau_c
    nsteps = int(math.ceil(SATURATION_K * (cal.tau_c + tm + tchi) / RECORD_FRAC / 10000)) * 10000
    arms = [("chi0", 0.0, cal.f_top, "chi == 0, m == f(zeta): the ACTIVE arm"),
            ("chi1", 1.0, cal.f_floor, "chi == 1, m == f(0.3 zeta): the PASSIVE arm")]
    for name, chi0, m0, what in arms:
        v = base(cal, 1.0, nsteps=nsteps,
                 seed=1001,
                 **{"chi-config": "uniform", "chi0": chi0, "m0": round(m0, 4),
                    "tau-chi": round(tchi, 1), "tau-m": round(tm, 1)})
        hdr = HDR_A4.format(ftop=cal.f_top, mc=cal.mc, ffl=cal.f_floor, arm=what,
                            cross=math.log(cal.f_top / cal.mc),
                            crossteps=math.log(cal.f_top / cal.mc) * tm,
                            sat=saturation_note(nsteps, cal.tau_c, tm, tchi))
        made.append(write_case(os.path.join(out, "cw_step_a4", name), hdr, v))
    return made


# ------------------------------------------------------------------- B1 / B2

TAU_M_GRID = [0.3, 0.45, 0.68, 1, 1.5, 2.2, 3.3, 4.7, 6.8, 10, 15, 22, 30]
CHI_LENGTH_OVER_L_P = 2.0        # correlation length of the mixed start


def mem_freeze(cal):
    """How long to hold the memory source off, in steps.

    THE FIRST B WAVE WAS INVALIDATED BY NOT HAVING THIS. Every run starts from rest, so
    P = 0, so g = 0, and m decays as exp(-t/tau_m) until there is turbulence to feed it. The
    active start sits at m0 = f(zeta) and crosses mc after tau_m ln(f/mc) = 0.226 tau_m,
    while the flow needs 2-4 tau_c to reach u_rms -- the same number once tau_m ~ 10 tau_c.
    Measured in the cancelled wave: the chi == 0 arm left its start at 0.239 tau_m
    (tau_m = 10 tau_c) and 0.253 tau_m (15 tau_c) against the 0.226 pure decay predicts, and
    all four starts collapsed onto the passive phase at every tau_m >= 10.

    10 tau_c(zeta) is the compromise. It is 2.5-5x the spin-up at full activity and still
    2.5 tau_c at the FLOOR, where the passive arm lives and the clock runs 4x slower. Longer
    would be safer for the flow and worse for the mixed start, whose binary pattern diffuses
    while it waits: over this window Dbio takes its correlation length from 2 L_P = 7.3 to
    sqrt(7.3^2 + 2 Dbio t) = 9.9, a 35% coarsening that is physical (it is cell motility)
    but which there is no reason to inflate. It also sits far inside the 20% warm-up that
    the record window discards.
    """
    return int(round(MEM_FREEZE_TAU_C * cal.tau_c / 10.0)) * 10

HDR_B = """20260901 cw_step_{grp} -- the tau_m scan, {tchi_rule}.

THE QUESTION. Three initial conditions, one tau_m axis: where does each of them end up?
Below tau_x the memory is too noisy for two basins to exist and all three must converge on
one state; above it the two ends should hold apart and the mixed start has to pick a side.
tau_x is MEASURED, not assumed -- {taux_note}

THE THREE STARTS. All at the same `seed`, so the initial director field is bit-identical
across the whole group and the only thing that changes is (chi, m):
  a   chi == 0, m == f(zeta)      = {f_top:.4f}   the ACTIVE phase at its own fixed point
  b   chi == 1, m == f(0.3 zeta)  = {f_floor:.4f}   the PASSIVE phase at its own
  c   binary-noise, two chi-seeds, half and half by construction at correlation length
      {chilen:.1f} = 2 L_P, each half carrying the m of the phase it belongs to
Both uniform starts are EXACT fixed points of the hard step, which is what makes "did it
move" a clean question rather than a matter of degree.

TWO tau_chi RULES, run as two groups. B1 holds tau_chi at 0.3 tau_c so the phenotype
follows the memory almost instantly and tau_m is the only slow clock. B2 sets tau_chi =
tau_m, which is the case where the phenotype lags as much as the memory does. The campaign's
prediction is that the FATE does not depend on which -- tau_x is set by sigma_m alone, and
sigma_m does not know about tau_chi. A boundary that moves between B1 and B2 would refute
that, and the drift flag is what separates a real shift from slow relaxation.

THE MEMORY SOURCE IS HELD OFF for the first {freeze:.0f} steps = 10 tau_c. Without it the
run measures its own spin-up rather than its physics: from rest P = 0, so g = 0, and m
decays as exp(-t/tau_m) until there is turbulence to feed it, crossing mc after
tau_m ln(f_top/mc) = {crossfrac:.3f} tau_m. The first wave of this group was cancelled for
exactly this -- the active arm left its start at 0.239 tau_m against the 0.226 pure decay
predicts, and every start at tau_m >= 10 tau_c collapsed onto the passive phase. Transport
and diffusion of m stay on throughout, and a frozen m holds chi* fixed, so chi does not
need a freeze of its own.

THIS CASE: tau_m = {g:g} tau_c = {tm:.0f} steps, tau_chi = {tchi:.0f} steps, start {start}.
{sat}"""


def gen_b(out, cal, group):
    """group = 'b1' (tau_chi = 0.3 tau_c) or 'b2' (tau_chi = tau_m)."""
    if not cal.measured:
        raise RuntimeError("the B stages must be generated from a MEASURED calib.json "
                           "(cw_step_pick.py), not from the wave-1 extrapolation: mc and "
                           "the two starts' m are all read off the measured f table.")
    if cal.pmem_on_edge:
        raise RuntimeError("calib.json reports the contrast maximum ON THE EDGE of the "
                           "scanned pmem range: the threshold was chosen by the scan "
                           "bounds, not by the f table, and mc, the window and tau_x are "
                           "all truncated with it. Widen PMEM_COEFFS in cw_step_pick.py "
                           "and re-run A3 before spending 100 runs on this.")
    grid = TAU_M_GRID if group == "b1" else TAU_M_GRID[1:]
    L_P = cal.L_P
    chilen = CHI_LENGTH_OVER_L_P * L_P
    taux_note = (f"stage A2b put it at {cal.tau_x:.2f} tau_c."
                 if math.isfinite(cal.tau_x) else
                 "stage A2b could not bracket it inside the scanned range.")
    starts = [
        ("a", {"chi-config": "uniform", "chi0": 0.0, "m0": round(cal.f_top, 4)},
         "chi == 0 (active)"),
        ("b", {"chi-config": "uniform", "chi0": 1.0, "m0": round(cal.f_floor, 4)},
         "chi == 1 (passive)"),
    ]
    for k in (1, 2):
        starts.append((f"c{k}",
                       {"chi-config": "binary-noise", "chi-length": round(chilen, 2),
                        "chi-seed": k, "chi-lo": 0.0, "m-lo": round(cal.f_top, 4),
                        "chi-hi": 1.0, "m-hi": round(cal.f_floor, 4),
                        "chi0": 0.5, "m0": 0.5},
                       f"binary noise, chi-seed {k}"))
    made = []
    for g in grid:
        tm = g * cal.tau_c
        tchi = tm if group == "b2" else 0.3 * cal.tau_c
        nsteps = max(NSTEPS, int(math.ceil(
            SATURATION_K * (cal.tau_c + tm + tchi) / RECORD_FRAC / 10000)) * 10000)
        for name, over, what in starts:
            v = base(cal, 1.0, nsteps=nsteps, seed=1001,
                     **{"tau-m": round(tm, 1), "tau-chi": round(tchi, 1),
                        "mem-freeze-steps": mem_freeze(cal)}, **over)
            hdr = HDR_B.format(
                grp=group, tchi_rule=("tau_chi = 0.3 tau_c" if group == "b1"
                                      else "tau_chi = tau_m"),
                taux_note=taux_note, f_top=cal.f_top, f_floor=cal.f_floor,
                chilen=chilen, g=g, tm=tm, tchi=tchi, start=what,
                freeze=mem_freeze(cal),
                crossfrac=math.log(cal.f_top / cal.mc),
                sat=saturation_note(nsteps, cal.tau_c, tm, tchi))
            made.append(write_case(
                os.path.join(out, f"cw_step_{group}", f"tm{tag(g)}_{name}"), hdr, v))
    return made


# ------------------------------------------------------------------- SYM

# (tag, pmem in units of sigma_P, mc). Both mc values equalise the two phases' normalised
# margins (m_phase - mc)/sigma_m_phase at tau_m = 30 tau_c; s1 keeps the campaign's
# threshold and moves mc alone, s2 raises the threshold too.
SYM_VARIANTS = [("s1", 0.5, 0.2100), ("s2", 0.6, 0.1585)]
SYM_TAU_M = [10, 15, 22, 30]

HDR_SYM = """20260901 cw_step_sym -- can the two phases be made EQUALLY stable?

THE PROBLEM. At the campaign's operating point the two coexisting phases are nowhere near
equally stable. Measured at tau_m = 30 tau_c:

    active   <m> = 0.3097, sigma_m = 0.0191  ->  (m - mc)/sigma =  2.55
    passive  <m> = 0.0554, sigma_m = 0.0295  ->  (mc - m)/sigma =  6.97

so 0.54% of nodes in the active phase sit below the threshold at any instant against 1.6e-12
in the passive one, and it is always the active phase that decays. Two causes stack. mc was
chosen by the CUSP, which maximises K_max -- "three roots exist easily" -- and not symmetry,
leaving it 0.066 below f_top and 0.204 above f_floor. And the two phases have genuinely
different sigma_m, 0.0191 against 0.0295, because the passive phase runs at the floor where
L_P is 8.5 rather than 3.7 and the same Dbio removes proportionally less.

THE FIX, and it is nearly free. Equalising (m_act - mc)/sigma_act = (mc - m_pass)/sigma_pass
on the MEASURED numbers gives mc = 0.21 at the campaign's pmem, which takes both margins to
~5.2 sigma. The cost in reachability is essentially nil: the coexistence threshold moves
11.9 -> 11.7 tau_c, because tau_x has a minimum near mc = 0.24 and rises on both sides, so
the symmetric point sits beside it rather than up a slope. Raising pmem as well improves
BOTH -- at 0.6 sigma_P the symmetric mc = 0.1585 gives 7.1 sigma and tau_x = 9.9 -- because
f_floor falls and sqrt(f(1-f)) collapses with it. Not pushed past 0.6: by 0.7 sigma_P the
floor's f is 0.007, m piles up against 0 and the normal closure every one of these numbers
rests on stops being true.

THIS GROUP. Four tau_m spanning the measured boundary, both pure starts, two candidate
operating points. BOTH ARMS ARE RUN, because symmetry cuts both ways: the passive margin
drops 6.97 -> 5.2 sigma, and a change that rescues the active phase by destabilising the
passive one has not made the two comparable, it has swapped which one decays.

WHAT WOULD SETTLE IT: the active arm holding at tau_m = 15 tau_c, where the campaign's
threshold gives it a lifetime of 2.0 tau_m, with the passive arm still never departing.

THIS CASE: {tag} -- pmem = {pc:g} sigma_P = {pmem:.6f}, mc = {mc:.4f}
  (f_floor, f_top) = ({ff:.4f}, {ft:.4f}), margins {za:.2f} / {zp:.2f} sigma at 30 tau_c
  tau_m = {g:g} tau_c = {tm:.0f} steps, start {start}
{sat}"""


def gen_sym(out, cal):
    if not cal.measured:
        raise RuntimeError("cw_step_sym must be generated from a MEASURED calib.json")
    made = []
    for vtag, pc, mc in SYM_VARIANTS:
        row = [t for t in cal.f_table if abs(t["coeff"] - pc) < 1e-9]
        if not row:
            raise RuntimeError(f"calib.json has no f-table row at pmem = {pc} sigma_P")
        row = row[0]
        ff, ft, pmem = row["f_floor"], row["f_top"], row["pmem"]
        # margins reported in the header, on the same model the choice was made with
        sa = math.sqrt(ft * (1 - ft)) * 0.0407 / math.sqrt(0.3272 * (1 - 0.3272))
        sp = math.sqrt(ff * (1 - ff)) * 0.1279 / math.sqrt(0.0570 * (1 - 0.0570))
        za, zp = (ft - mc) / sa, (mc - ff) / sp
        for g in SYM_TAU_M:
            tm = g * cal.tau_c
            tchi = 0.3 * cal.tau_c
            nsteps = max(NSTEPS, int(math.ceil(
                SATURATION_K * (cal.tau_c + tm + tchi) / RECORD_FRAC / 10000)) * 10000)
            for name, over, what in (
                    ("a", {"chi0": 0.0, "m0": round(ft, 4)}, "chi == 0 (active)"),
                    ("b", {"chi0": 1.0, "m0": round(ff, 4)}, "chi == 1 (passive)")):
                v = base(cal, 1.0, nsteps=nsteps, seed=1001,
                         **{"tau-m": round(tm, 1), "tau-chi": round(tchi, 1),
                            "mem-freeze-steps": mem_freeze(cal),
                            "chi-config": "uniform",
                            "pmem": round(pmem, 6), "mc": round(mc, 4)}, **over)
                hdr = HDR_SYM.format(tag=vtag, pc=pc, pmem=pmem, mc=mc, ff=ff, ft=ft,
                                     za=za, zp=zp, g=g, tm=tm, start=what,
                                     sat=saturation_note(nsteps, cal.tau_c, tm, tchi))
                made.append(write_case(
                    os.path.join(out, "cw_step_sym", f"{vtag}_tm{tag(g)}_{name}"),
                    hdr, v))
    return made


# ------------------------------------------------------------------- FINE

# 20 points, log-uniform from 0.3 to 10 tau_c: the region where the campaign's B1 grid has
# only 10 points and where, at the SYMMETRIC operating point, the whole transition happens.
FINE_TAU_M = [0.30, 0.36, 0.44, 0.53, 0.63, 0.76, 0.92, 1.11, 1.33, 1.60,
              1.93, 2.32, 2.79, 3.35, 4.03, 4.85, 5.83, 7.01, 8.32, 10.0]
FINE_MC = 0.2100          # s1: the symmetric threshold
FINE_PMEM_COEFF = 0.5     # unchanged from the campaign

HDR_FINE = """20260901 cw_step_fine -- HOW the transition happens, at the symmetric point.

WHERE THIS SITS. The campaign scanned 0.3-30 tau_c at mc = 0.2610, where the two phases are
2.7x apart in stability and everything below 22 tau_c collapses to passive. The 16-run
symmetry test then moved mc to 0.2100, equalising the phases (5.27 against 5.26 sigma at
30 tau_c) and dropping the coexistence threshold to 10-15 tau_c -- but it only sampled four
tau_m. Between 0.3 and 10 the symmetric point runs from a mixed state to a nearly pure
ACTIVE one, the opposite direction to the campaign's, and nothing has looked at the shape of
that crossover.

THE GRID is 20 points log-uniform over 0.3-10 tau_c, a factor 1.2 apart, against the
campaign's 10 points over the same range at a factor 1.7. Three starts: both pure ones,
which bracket whatever is reachable, and one binary-noise start, because at the campaign's
mc the mixed starts NEVER chose the active phase and whether symmetry changes that is a
question about the size of the active basin rather than its existence.

WHAT THE SHAPE SHOULD SETTLE. Below the coexistence threshold there is one state and all
three starts must agree on it; the interesting content is then whether <chi> moves
CONTINUOUSLY as tau_m grows -- a single fixed point sliding -- or jumps. Above it the pure
starts separate. The campaign's own data says the single state slides (0.65 -> 0.99 over
0.3-15 tau_c at mc = 0.2610), so a jump here would mean the symmetric point is qualitatively
different and not merely shifted.

THIS CASE: tau_m = {g:g} tau_c = {tm:.0f} steps, tau_chi = {tchi:.0f} steps, start {start}.
mc = {mc:.4f}, pmem = {pmem:.6f} ({pc:g} sigma_P), memory frozen for the first {freeze:.0f}
steps while the flow spins up.
{sat}"""


def gen_fine(out, cal):
    if not cal.measured:
        raise RuntimeError("cw_step_fine must be generated from a MEASURED calib.json")
    row = [t for t in cal.f_table if abs(t["coeff"] - FINE_PMEM_COEFF) < 1e-9]
    if not row:
        raise RuntimeError(f"no f-table row at pmem = {FINE_PMEM_COEFF} sigma_P")
    row = row[0]
    ff, ft, pmem = row["f_floor"], row["f_top"], row["pmem"]
    chilen = CHI_LENGTH_OVER_L_P * cal.L_P
    starts = [
        ("a", {"chi-config": "uniform", "chi0": 0.0, "m0": round(ft, 4)}, "chi == 0 (active)"),
        ("b", {"chi-config": "uniform", "chi0": 1.0, "m0": round(ff, 4)}, "chi == 1 (passive)"),
        ("c1", {"chi-config": "binary-noise", "chi-length": round(chilen, 2),
                "chi-seed": 1, "chi-lo": 0.0, "m-lo": round(ft, 4),
                "chi-hi": 1.0, "m-hi": round(ff, 4), "chi0": 0.5, "m0": 0.5},
         "binary noise, chi-seed 1"),
    ]
    made = []
    for g in FINE_TAU_M:
        tm = g * cal.tau_c
        tchi = 0.3 * cal.tau_c
        nsteps = max(NSTEPS, int(math.ceil(
            SATURATION_K * (cal.tau_c + tm + tchi) / RECORD_FRAC / 10000)) * 10000)
        for name, over, what in starts:
            v = base(cal, 1.0, nsteps=nsteps, seed=1001,
                     **{"tau-m": round(tm, 1), "tau-chi": round(tchi, 1),
                        "mem-freeze-steps": mem_freeze(cal),
                        "pmem": round(pmem, 6), "mc": FINE_MC}, **over)
            hdr = HDR_FINE.format(g=g, tm=tm, tchi=tchi, start=what, mc=FINE_MC,
                                  pmem=pmem, pc=FINE_PMEM_COEFF, freeze=mem_freeze(cal),
                                  sat=saturation_note(nsteps, cal.tau_c, tm, tchi))
            made.append(write_case(
                os.path.join(out, "cw_step_fine", f"tm{tag(g)}_{name}"), hdr, v))
    return made


# ------------------------------------------------------------------- TCHI

TCHI_TAU_M = [0.3, 3, 10, 15, 30]
# (tag, rule, label). "fixed" is in units of tau_c, "prop" a multiple of tau_m.
TCHI_RULES = [("t1", "fixed", 1.0), ("th", "prop", 0.5), ("t0p3", "fixed", 0.3)]
TCHI_ONLY = {"t0p3": [3]}     # 0.3 tau_c already exists everywhere except tau_m = 3

HDR_TCHI = """20260902 cw_step_tchi -- does the phenotype clock move the bifurcation?

WHAT IS ALREADY KNOWN. The campaign predicted the fate would not depend on tau_chi at all,
since tau_x is set by sigma_m and sigma_m does not know about tau_chi. It does depend on it:
at mc = 0.2610 the tau_chi = tau_m group (B2) held the active phase at 22 tau_c where the
tau_chi = 0.3 tau_c group (B1) did not, and B2's active lifetime was longer at EVERY tau_m.
The proposed mechanism is that a slow chi low-passes its response to m's fluctuations and so
damps the positive feedback m -> chi -> activity -> f -> m.

WHAT THIS GROUP ADDS. Two more phenotype clocks on the SYMMETRIC operating point
(mc = 0.2100), where the whole 0.3-30 axis has so far been measured at one tau_chi only:

    tau_chi = 1.0 tau_c    a fixed clock, 3.3x slower than the existing 0.3 tau_c
    tau_chi = 0.5 tau_m    a proportional clock, half the memory time

Together with the existing 0.3 tau_c data that is three rules spanning both kinds -- fixed
and proportional -- at five tau_m. If the fate really is a function of sigma_m alone, all
three collapse; if the low-pass picture is right, the slower the phenotype the more robust
the active phase, and the proportional rule should separate earliest at large tau_m where it
is slowest of all.

The 0.3 tau_c rule is re-run at tau_m = 3 only: the fine scan has 2.79 and 3.35 but not 3,
and an interpolated point is not worth having when the comparison is the whole purpose.

THIS CASE: {rule} -- tau_chi = {tchi:.1f} steps = {gchi:.2f} tau_c = {frac:.2f} tau_m,
tau_m = {g:g} tau_c = {tm:.0f} steps, start {start}.
{sat}"""


def gen_tchi(out, cal):
    if not cal.measured:
        raise RuntimeError("cw_step_tchi must be generated from a MEASURED calib.json")
    row = [t for t in cal.f_table if abs(t["coeff"] - FINE_PMEM_COEFF) < 1e-9][0]
    ff, ft, pmem = row["f_floor"], row["f_top"], row["pmem"]
    chilen = CHI_LENGTH_OVER_L_P * cal.L_P
    starts = [
        ("a", {"chi-config": "uniform", "chi0": 0.0, "m0": round(ft, 4)}, "chi == 0"),
        ("b", {"chi-config": "uniform", "chi0": 1.0, "m0": round(ff, 4)}, "chi == 1"),
        ("c1", {"chi-config": "binary-noise", "chi-length": round(chilen, 2),
                "chi-seed": 1, "chi-lo": 0.0, "m-lo": round(ft, 4),
                "chi-hi": 1.0, "m-hi": round(ff, 4), "chi0": 0.5, "m0": 0.5},
         "binary noise"),
    ]
    made = []
    for rtag, kind, coef in TCHI_RULES:
        for g in TCHI_ONLY.get(rtag, TCHI_TAU_M):
            tm = g * cal.tau_c
            tchi = coef * cal.tau_c if kind == "fixed" else coef * tm
            # tau_chi must resolve on the integration step as well as on the physics
            if tchi < 20.0:
                raise RuntimeError(f"{rtag} at tau_m = {g}: tau_chi = {tchi:.1f} steps is "
                                   f"under the 20-step floor")
            nsteps = max(NSTEPS, int(math.ceil(
                SATURATION_K * (cal.tau_c + tm + tchi) / RECORD_FRAC / 10000)) * 10000)
            for name, over, what in starts:
                v = base(cal, 1.0, nsteps=nsteps, seed=1001,
                         **{"tau-m": round(tm, 1), "tau-chi": round(tchi, 1),
                            "mem-freeze-steps": mem_freeze(cal),
                            "pmem": round(pmem, 6), "mc": FINE_MC}, **over)
                hdr = HDR_TCHI.format(rule=rtag, tchi=tchi, gchi=tchi / cal.tau_c,
                                      frac=tchi / tm, g=g, tm=tm, start=what,
                                      sat=saturation_note(nsteps, cal.tau_c, tm, tchi))
                made.append(write_case(
                    os.path.join(out, "cw_step_tchi", f"{rtag}_tm{tag(g)}_{name}"), hdr, v))
    return made


# ------------------------------------------------------------------- LONG

LONG_TAU_M = [3.35, 8.32]
LONG_NSTEPS = 500000

HDR_LONG = """20260902 cw_step_long -- are these two states settled, or drifting slowly?

The fine scan gave every point 100000 steps, and at these two tau_m that leaves the answer
open in different ways.

  tau_m = 3.35 tau_c   the three starts landed at 0.505 / 0.397 / 0.466 -- a 0.108 spread in
                       what is supposed to be a single state, with points flagged drifting.
                       Either the mixed state has genuine seed-to-seed scatter, or 100000
                       steps is not long enough for it to converge.
  tau_m = 8.32 tau_c   the three agreed to 0.003 at 0.04, just inside the active phase and
                       one grid point below where the campaign's mc collapsed. If the active
                       phase here is metastable rather than stable, 5x the run length is
                       where that shows.

Everything is unchanged from the fine scan except the run length: same mc = 0.2100, same
pmem, same tau_chi = 0.3 tau_c, same three starts, same seed. 500000 steps puts the record
window at 529 tau_m for the shorter memory and 213 for the longer, against 106 and 43
before.

OUTPUT CADENCE IS RESCALED, and it has to be. At the fine scan's nvideo = 45 a 500000-step
run writes 11111 video frames, ~7 GB per run before transcoding; at ntracer = 10 the tracer
file reaches 1.6 GB. The video clock goes to 1 tau_c (2222 frames, the same count the
existing clips have, so they stay comparable to watch) and the tracer interval to 25 steps,
which is still ~9 samples per correlation time -- enough to resolve the Lagrangian 1/e.
ninfo = 20000 keeps 25 full frames, above cw_part's 20-frame floor.

THIS CASE: tau_m = {g:g} tau_c = {tm:.0f} steps, start {start}, {nst} steps.
{sat}"""


def gen_long(out, cal):
    if not cal.measured:
        raise RuntimeError("cw_step_long must be generated from a MEASURED calib.json")
    row = [t for t in cal.f_table if abs(t["coeff"] - FINE_PMEM_COEFF) < 1e-9][0]
    ff, ft, pmem = row["f_floor"], row["f_top"], row["pmem"]
    chilen = CHI_LENGTH_OVER_L_P * cal.L_P
    starts = [
        ("a", {"chi-config": "uniform", "chi0": 0.0, "m0": round(ft, 4)}, "chi == 0"),
        ("b", {"chi-config": "uniform", "chi0": 1.0, "m0": round(ff, 4)}, "chi == 1"),
        ("c1", {"chi-config": "binary-noise", "chi-length": round(chilen, 2),
                "chi-seed": 1, "chi-lo": 0.0, "m-lo": round(ft, 4),
                "chi-hi": 1.0, "m-hi": round(ff, 4), "chi0": 0.5, "m0": 0.5},
         "binary noise"),
    ]
    made = []
    for g in LONG_TAU_M:
        tm = g * cal.tau_c
        tchi = 0.3 * cal.tau_c
        for name, over, what in starts:
            v = base(cal, 1.0, nsteps=LONG_NSTEPS, seed=1001,
                     **{"tau-m": round(tm, 1), "tau-chi": round(tchi, 1),
                        "mem-freeze-steps": mem_freeze(cal),
                        "pmem": round(pmem, 6), "mc": FINE_MC,
                        # rescaled for a 5x longer run -- see the header
                        "ninfo": 20000,
                        "nvideo": int(round(cal.tau_c / 5.0)) * 5,
                        "ntracer": 25}, **over)
            hdr = HDR_LONG.format(g=g, tm=tm, start=what, nst=LONG_NSTEPS,
                                  sat=saturation_note(LONG_NSTEPS, cal.tau_c, tm, tchi))
            made.append(write_case(
                os.path.join(out, "cw_step_long", f"tm{tag(g)}_{name}"), hdr, v))
    return made


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage",
                    choices=["a1", "a2b", "a4", "b1", "b2", "sym", "fine", "tchi", "long"])
    ap.add_argument("--out", default="cases/20260901")
    ap.add_argument("--calib", default=None,
                    help="measured calib.json from stage A3; without it the wave-1 "
                         "extrapolation is used, which is only good enough for stage A")
    a = ap.parse_args()

    cal = Calib(a.calib)
    print(f"scales: tau_c = {cal.tau_c:.1f} steps, sigma_P = {cal.sigma_P:.5f}, "
          f"u_rms = {cal.u_rms:.5f}, pmem = {cal.pmem:.5f}, mc = {cal.mc:.4f}  "
          f"[{'MEASURED' if cal.measured else 'wave-1 extrapolation'}]")
    if a.stage in ("a4",) and not cal.measured:
        print("  NOTE: a closed-loop stage on the extrapolation. A4 is a yes/no test and "
              "tolerates it; the B stages must not be generated this way.")
    if a.stage == "long":
        made = gen_long(a.out, cal)
    elif a.stage == "tchi":
        made = gen_tchi(a.out, cal)
    elif a.stage == "fine":
        made = gen_fine(a.out, cal)
    elif a.stage == "sym":
        made = gen_sym(a.out, cal)
    elif a.stage in ("b1", "b2"):
        made = gen_b(a.out, cal, a.stage)
    else:
        made = {"a1": gen_a1, "a2b": gen_a2b, "a4": gen_a4}[a.stage](a.out, cal)
    for p in made:
        print(f"  {p}")
    print(f"{len(made)} runcards")


if __name__ == "__main__":
    main()
