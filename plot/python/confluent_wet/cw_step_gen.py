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
TRACER_COUNT = 2000
STORE_SIGMA = 6.0             # video stores +/- 6 sigma_P, the renderer displays +/- 3

ORDER = ["model", "nsteps", "nstart", "ninfo", "LX", "LY", "bc", "seed", "",
         "Gamma", "xi", "CC", "LL", "tau", "rho", "friction", "zeta", "zeta0-frac",
         "open-loop", "zeta-open", "",
         "angle", "noise", "initial-order", "director-config", "defect-sep", "",
         "Dbio", "chi-config", "chi0", "chi-noise", "chi-length", "chi-seed",
         "chi-lo", "chi-hi", "m-lo", "m-hi", "chi-block",
         "tau-chi", "chi-width", "switch-sign", "mc", "m0", "chi-freeze-steps", "",
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

THIS CASE: tau_m = {g:g} tau_c = {tm:.0f} steps, tau_chi = {tchi:.0f} steps, start {start}.
{sat}"""


def gen_b(out, cal, group):
    """group = 'b1' (tau_chi = 0.3 tau_c) or 'b2' (tau_chi = tau_m)."""
    if not cal.measured:
        raise RuntimeError("the B stages must be generated from a MEASURED calib.json "
                           "(cw_step_pick.py), not from the wave-1 extrapolation: mc and "
                           "the two starts' m are all read off the measured f table.")
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
                     **{"tau-m": round(tm, 1), "tau-chi": round(tchi, 1)}, **over)
            hdr = HDR_B.format(
                grp=group, tchi_rule=("tau_chi = 0.3 tau_c" if group == "b1"
                                      else "tau_chi = tau_m"),
                taux_note=taux_note, f_top=cal.f_top, f_floor=cal.f_floor,
                chilen=chilen, g=g, tm=tm, tchi=tchi, start=what,
                sat=saturation_note(nsteps, cal.tau_c, tm, tchi))
            made.append(write_case(
                os.path.join(out, f"cw_step_{group}", f"tm{tag(g)}_{name}"), hdr, v))
    return made


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["a1", "a2b", "a4", "b1", "b2"])
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
    if a.stage in ("b1", "b2"):
        made = gen_b(a.out, cal, a.stage)
    else:
        made = {"a1": gen_a1, "a2b": gen_a2b, "a4": gen_a4}[a.stage](a.out, cal)
    for p in made:
        print(f"  {p}")
    print(f"{len(made)} runcards")


if __name__ == "__main__":
    main()
