#!/usr/bin/env python3
"""Generate the runcards of the 2026-09 confluent-wet memory campaign.

Every run.dat in the campaign is written by this script, from one of two sources: the
hand-set activity ladder (stage A1, which measures the campaign's own units) or calib.json
(stage A3's output, which fixes every parameter downstream). Nothing is typed by hand, and
nothing downstream re-derives a number that calib.json already holds -- that is the point.
A parameter that appears in two runcards with two values is the failure mode this replaces.

    cw_gen.py ladder  <outdir> --zeta 0.15 [--tau-c-pred 540] ...   stage A1  (6 runs)
    cw_gen.py a4      <outdir> --calib calib.json                   stage A4  (2 runs)
    cw_gen.py g2      <outdir> --calib calib.json                   group G2  (45 runs)
    cw_gen.py rest    <outdir> --calib calib.json --mc-star MC      G1,G3,G4,G5 (90 runs)

THE UNITS. tau_c(zeta) = d/u_rms at the top of the ladder -- the time the flow takes to
carry material one defect spacing, i.e. the interval between defect encounters along a
material trajectory. Every time in the campaign is a multiple of it: tau_chi = 0.3 tau_c,
tau_m in {0.3, 1, 3, 10, 30} tau_c, the video frame interval 0.2 tau_c, the analysis frame
interval 5 tau_c. A1 cannot know it yet and uses a predicted value for its own cadence; from
A4 on, the measured one is used, and the difference is recorded rather than hidden.
"""
import argparse
import json
import os

# ---- the frozen baseline. Adopted 2026-08-31 (dev_log): Ma inside budget for the first
# ---- time, xi_N = 2.0 the measured sweet spot, rho = 40 the only density our LB survives.
FROZEN = dict(model="confluent-wet", LX=800, LY=800, bc=0, nstart=0, nsteps=100000,
              Gamma=0.05, xi=0.4, CC=0.00625, LL=0.05, tau=2.0, rho=40, friction=0.0,
              angle=0, noise=0.05, initial_order=1, director_config="uniform",
              defect_sep=0, video_stride=2, frame_light=1, zeta0_frac=0.3)

LADDER = (1.0, 0.8, 0.6, 0.45, 0.3, 0.15)
TAU_M_RATIOS = (0.3, 1.0, 3.0, 10.0, 30.0)
SEEDS_SINGLE = (1001, 1002)          # single-phase initial conditions
SEEDS_STRIPE = (1001, 1002)
SEEDS_MIXED = (1001, 1002, 1003)     # the mixed start is the one whose result is a fraction

# measured anchor for the A1 predictions: cw_lowk, L = 500, zeta_eff = 0.05 (dev_log 20260830)
ANCHOR = dict(zeta_eff=0.05, u_rms=0.025 / (3 ** 0.5), sigma_P_over_zeta_eff=0.625)


def tag(x, nd=2):
    """0.15 -> 0p15, -0.06 -> m0p06; the project's directory-name convention."""
    s = f"{x:.{nd}f}".rstrip("0").rstrip(".")
    if s in ("", "-0"):
        s = "0"
    return s.replace("-", "m").replace(".", "p")


ORDER = ["model", "nsteps", "nstart", "ninfo", "LX", "LY", "bc", "seed", "",
         "Gamma", "xi", "CC", "LL", "tau", "rho", "friction", "zeta", "zeta0-frac",
         "open-loop", "zeta-open", "",
         "angle", "noise", "initial-order", "director-config", "defect-sep", "",
         "Dbio", "chi-config", "chi0", "chi-noise", "chi-length",
         "chi-lo", "chi-hi", "m-lo", "m-hi", "chi-block",
         "tau-chi", "chi-width", "switch-sign", "mc", "m0", "chi-freeze-steps", "",
         "tau-m", "pmem", "pmem-width", "",
         "nvideo", "video-stride", "video-p-scale", "video-u-scale", "frame-light"]


def write_case(path, header, vals):
    os.makedirs(path, exist_ok=True)
    lines = [f"# {h}" for h in header.strip("\n").split("\n")]
    for key in ORDER:
        if key == "":
            lines.append("")
        elif key in vals:
            v = vals[key]
            v = f"{v:.6e}" if isinstance(v, float) and (abs(v) < 1e-3 and v != 0) else v
            lines.append(f"{key:<17}= {v}")
    unknown = set(vals) - set(ORDER)
    if unknown:
        raise RuntimeError(f"keys not in ORDER (they would be silently dropped): {unknown}")
    with open(os.path.join(path, "run.dat"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


def base(cal=None, **over):
    """The frozen block plus whatever calib.json fixes; runcard spelling, with dashes."""
    v = {"model": FROZEN["model"], "nsteps": FROZEN["nsteps"], "nstart": FROZEN["nstart"],
         "LX": FROZEN["LX"], "LY": FROZEN["LY"], "bc": FROZEN["bc"],
         "Gamma": FROZEN["Gamma"], "xi": FROZEN["xi"], "CC": FROZEN["CC"],
         "LL": FROZEN["LL"], "tau": FROZEN["tau"], "rho": FROZEN["rho"],
         "friction": FROZEN["friction"], "zeta0-frac": FROZEN["zeta0_frac"],
         "angle": FROZEN["angle"], "noise": FROZEN["noise"],
         "initial-order": FROZEN["initial_order"],
         "director-config": FROZEN["director_config"], "defect-sep": FROZEN["defect_sep"],
         "video-stride": FROZEN["video_stride"], "frame-light": FROZEN["frame_light"]}
    if cal:
        v.update({"zeta": cal["zeta"], "zeta0-frac": cal["zeta0_frac"],
                  "Dbio": round(cal["Dbio"], 6),
                  "tau-chi": round(cal["tau_chi"], 1),
                  "chi-width": round(cal["chi_width"], 5),
                  "switch-sign": cal["switch_sign"],
                  "pmem": cal["pmem"], "pmem-width": cal["pmem_width"],
                  "ninfo": cal["video"]["ninfo"], "nvideo": cal["video"]["nvideo"],
                  "video-p-scale": cal["video"]["p_scale"],
                  "video-u-scale": cal["video"]["u_scale"]})
    v.update(over)
    return v


# ------------------------------------------------------------------ stage A1

HDR_A1 = """20260831 cw_mem_a1 -- stage A1 of the memory campaign: the ACTIVITY LADDER.

Six OPEN-LOOP runs. m and chi evolve exactly as they will in production, but the activity
is pinned at zeta-open and chi is not fed back, so this measures the RESPONSE of the layer
to a prescribed activity rather than a fixed point. What comes out of it:

  tau_c(zeta) = d/u_rms at the top rung -- the smallest tau_c on the ladder, and from here
                on THE time unit of the whole campaign
  sigma_P(zeta), u_rms(zeta)   the fixed colour scales of every video in the campaign
  f(pmem; zeta_eff)            the fraction of time a lattice point spends above a pressure
                               threshold; stage A3 solves chi = chi*(f(zeta_eff(chi))) on it

CHECKPOINT 1 is the r = 0.30 rung, which is the activity FLOOR zeta_0 = 0.3 zeta that the
passive phenotype will sit at. It must still be turbulent -- a stable non-zero defect count.
If it is not, the floor has to be raised to 0.4 zeta and the ladder re-run, because a
passive phase that has stopped stirring generates no pressure, can never be re-excited, and
is an absorbing state rather than a phase.

THE ACTIVITY WAS RAISED for this campaign, from the zeta = 0.1 of the 20260830 working point
to zeta = {zeta:g}. The reason is run length measured in the layer's own clock: at zeta = 0.1,
tau_c was 1541 steps, so a 100000-step run is only 65 tau_c and the tau_m = 30 tau_c corner
gets 1.7 memory times of recording. tau_c falls roughly as 1/zeta (u_rms ~ sqrt(zeta) and the
defect spacing falls as A^-0.45), so 1.5x the activity buys ~2.8x the tau_c count. The price
is Mach: Ma ~ sqrt(zeta_eff), so the top rung is PREDICTED at Ma = {ma:.3f} against the 0.045
weakly-compressible budget -- inside it, but not by much, and the prediction is an
extrapolation. This ladder measures it. If the top rung comes out over budget, the campaign
drops one rung: the r = 0.8 run is already a complete calibration of zeta_eff = {z08:g}, so
the nominal zeta can be lowered WITHOUT re-running anything.

THIS CASE: r = zeta_eff/zeta = {r:g}, zeta_eff = {ze:g}, A = zeta_eff/CC = {A:.1f}.

CADENCE. ninfo and nvideo are set from the PREDICTED tau_c = {tcp:.0f} steps (5 tau_c and
0.2 tau_c). From stage A4 on they are set from the measured value; A1's own cadence being
~20% off costs nothing, because A1 measures rather than compares.

VIDEO SCALES are per rung here and predicted (6 sigma_P, 6 u_rms with sigma_P = 0.625
zeta_eff and u_rms extrapolated from the 20260830 anchor). They are the STORED clipping
limits, not the colour bar: the renderer displays the inner +/- 3 sigma_P using the MEASURED
sigma_P. Storing wide is what keeps the tails -- and hence f near the threshold -- intact."""


def gen_ladder(out, zeta, tau_c_pred):
    u_anchor = ANCHOR["u_rms"]
    for r in LADDER:
        ze = zeta * r
        u_pred = u_anchor * (ze / ANCHOR["zeta_eff"]) ** 0.5
        sP_pred = ANCHOR["sigma_P_over_zeta_eff"] * ze
        v = base(zeta=zeta, ninfo=int(round(5 * tau_c_pred)),
                 nvideo=int(round(0.2 * tau_c_pred)), seed=1001,
                 **{"open-loop": 1, "zeta-open": ze,
                    "Dbio": round(9.0 / tau_c_pred, 6),
                    "chi-config": "uniform", "chi0": 0.5, "m0": 0.5,
                    "chi-noise": 0, "chi-length": 0,
                    "tau-chi": round(0.3 * tau_c_pred, 1), "chi-width": 0.15,
                    "switch-sign": -1, "mc": 0.5,
                    "tau-m": round(3.0 * tau_c_pred, 1),
                    "pmem": 0.0, "pmem-width": round(sP_pred / 4, 6),
                    "video-p-scale": round(6 * sP_pred, 6),
                    "video-u-scale": round(6 * u_pred, 6)})
        hdr = HDR_A1.format(zeta=zeta, r=r, ze=ze, A=ze / FROZEN["CC"], tcp=tau_c_pred,
                            ma=u_anchor * (zeta / ANCHOR["zeta_eff"]) ** 0.5 * (3 ** 0.5),
                            z08=zeta * 0.8)
        write_case(os.path.join(out, f"r{tag(r)}"), hdr, v)
    print(f"A1: {len(LADDER)} runcards under {out}")


# ------------------------------------------------------------------ stage A4

HDR_A4 = """20260831 cw_mem_a4 -- stage A4: the BISTABILITY QUICK CHECK.

Two closed-loop runs at tau_m = 30 tau_c with the stage-A3 parameters, started from the two
ends: chi == 0 (fully active) and chi == 1 (at the activity floor). CHECKPOINT 2 is that
their record-window <chi> differ by more than 0.4. If they do not, the mean-field roots are
not realised and A3 is redone with the loop gain raised to G = 3.5.

WHY THIS CHECK EXISTS. The A3 fixed-point equation chi = chi*(f(zeta_eff(chi))) is a
MEAN-FIELD statement: it uses only <m>, and it is blind to var(m). A wide spread of m across
the box smears chi*(m) and can wash out a bistability the scalar equation predicts. Nothing
in A3 can see that; only running both ends can.

EACH END STARTS AT ITS OWN SELF-CONSISTENT MEMORY. m0 = f(zeta_eff(chi_init)): the chi == 0
run starts at m0 = {m_hi_act:.4f}, the memory a fully active layer actually accumulates, and
the chi == 1 run at m0 = {m_floor:.4f}. Starting both at m = 0.5 would have both memories
drifting towards their own equilibrium for ~tau_m = 30 tau_c, i.e. for a third of the record
window, and a separation that merely reflects that transient would be indistinguishable from
genuine bistability.

A GAIN SCAN, NOT A SINGLE PAIR. The campaign's recipe is "run the two ends; if <chi>
differs by less than 0.4, go back to A3 and raise G to 3.5". On the 20260831 ladder the
mean-field map already answers the first round: at the A3 chi-width the map has ONE root, so
the two ends must converge and checkpoint 2 cannot pass. Running that round anyway to watch
it fail, then running the remedy, is two serial waves for one bit of information. Both are
run here instead, in one wave, plus one stiffer setting as the fallback -- which also turns
the checkpoint into a measurement: at which loop gain does the separation actually appear,
and does it appear where the mean-field says it should?

WHY THE PROXY MISLED. G is built from the MAXIMUM of |df/dchi| over the reachable activity
range, but bistability is decided by the slope where the fixed point SITS. Here the maximum
is at chi = 0.88 (the layer near its activity floor) and the fixed point is at chi = 0.17, so
G = 1.78 cleared the G >= 1.5 test with a map that has a single root.

THIS CASE: chi == {chi0:g}, m0 = {m0:.4f}, G = {G:g} (chi-width = {w:.4f}), mean-field roots
{rootnote}."""


def gen_a4(out, cal, g_list=None, wchi_list=None):
    """One pair of runs per switch width. The width, not the gain, is the knob now.

    G was the campaign's handle on the loop, and it is the wrong one: it is built from the
    maximum of |df/dchi| anywhere in the reachable range, so it says nothing about the slope
    where the fixed points actually sit. The 20260831 point that A3 now picks has G = 1.51 --
    below the campaign's own G >= 1.5 line by a hair -- and three well separated roots, while
    the first pass had G = 1.78 and a single root. chi-width is what is actually varied, so
    it is what these runs are indexed by.
    """
    byw = {round(e["chi_width"], 4): e for e in cal.get("g_scan", [])}
    scan = {round(e["G"], 3): e for e in cal.get("g_scan", [])}
    if wchi_list:
        pairs = [(cal["max_df_dchi"] / (2 * w), w) for w in wchi_list]
    elif g_list:
        pairs = [(G, (scan[round(G, 3)]["chi_width"] if round(G, 3) in scan
                      else cal["max_df_dchi"] / (2 * G))) for G in g_list]
    else:
        pairs = [(cal["G"], cal["chi_width"])]
    n = 0
    for G, w in pairs:
        e = byw.get(round(w, 4)) or scan.get(round(G, 3))
        if e and e["n_roots"] >= 3:
            r = e["roots"]
            note = (f"chi_lo = {r[0]:.4f}, chi_mid = {r[len(r)//2]:.4f}, "
                    f"chi_hi = {r[-1]:.4f} (separation {e['separation']:.3f}, weaker basin "
                    f"{e.get('basin', float('nan')):.3f}, R = "
                    f"{e.get('robustness', float('nan')):.1f})")
        elif e:
            note = f"a SINGLE root at chi = {e['roots'][0]:.4f} -- the two ends must converge"
        else:
            note = "not tabulated"
        for chi0, m0 in ((0.0, cal["f_at_full_zeta"]), (1.0, cal["f_at_floor"])):
            v = base(cal, seed=SEEDS_SINGLE[0],
                     **{"chi-config": "uniform", "chi0": chi0, "chi-noise": 0,
                        "chi-length": 0, "m0": round(m0, 6),
                        "mc": round(cal["mc_0"], 6), "chi-width": round(w, 5),
                        "tau-m": round(30.0 * cal["tau_c"], 1)})
            write_case(os.path.join(out, f"w{tag(w,3)}_chi{int(chi0)}"),
                       HDR_A4.format(chi0=chi0, m0=m0, G=G, w=w, rootnote=note,
                                     m_hi_act=cal["f_at_full_zeta"],
                                     m_floor=cal["f_at_floor"]), v)
            n += 1
    print(f"A4: {n} runcards under {out}")


# ------------------------------------------------------------------- group G2

HDR_G2 = """20260831 cw_mem_g2 -- group G2: LOCATING THE COEXISTENCE POINT mc*.

The memory threshold mc is scanned across mc_0 = {mc0:.4f}, the midpoint of the reachable
memory range (f at the floor {f_flr:.3f} .. f at full activity {f_top:.3f}). Every run starts
from the RANDOM MIXED state -- 20-cell blocks, exactly half at (chi_hi, m_hi) and half at
(chi_lo, m_lo), the two mean-field roots each prepared with its own self-consistent memory --
so neither phase is favoured by the initial condition.

mc* is read off the curve: the value of mc at which the record-window <chi> crosses
(chi_hi + chi_lo)/2 = {mid:.4f}. Away from mc* one phase invades and the box ends up in it;
at mc* the two phases have equal free energy in the loose sense that neither front moves,
and the run stays two-phase. That run IS the coexistence state, and mc* is what G3, G4 and
G5 are all run at.

EXACTLY HALF, NOT HALF ON AVERAGE. The block assignment is a shuffle, not an independent
coin per block: with 1600 blocks a coin leaves the initial area fraction scattered by
+/- 1.25%, and the quantity being measured is precisely an area fraction.

THREE SEEDS per point, because the answer here is a fraction and the run-to-run spread of a
fraction is the error bar on mc*.

THIS CASE: mc = mc_0 {dmc:+.2f} = {mc:.4f}, tau_m = {tmr:g} tau_c = {tm:.0f} steps,
seed {seed}."""


def gen_g2(out, cal):
    mid = 0.5 * (cal["chi_hi"] + cal["chi_lo"])
    n = 0
    for dmc in (-0.06, -0.03, 0.0, 0.03, 0.06):
        for tmr in (3.0, 10.0, 30.0):
            for seed in SEEDS_MIXED:
                mc = cal["mc_0"] + dmc
                v = base(cal, seed=seed,
                         **{"chi-config": "blocks", "chi-block": 20,
                            "chi-lo": round(cal["chi_lo"], 6),
                            "chi-hi": round(cal["chi_hi"], 6),
                            "m-lo": round(cal["f_lo"], 6), "m-hi": round(cal["f_hi"], 6),
                            "chi0": 0.5, "m0": 0.5,
                            "mc": round(mc, 6), "tau-m": round(tmr * cal["tau_c"], 1)})
                name = f"mc{tag(dmc)}_tm{tag(tmr,1)}_s{seed}"
                write_case(os.path.join(out, name),
                           HDR_G2.format(mc0=cal["mc_0"], f_flr=cal["f_at_floor"],
                                         f_top=cal["f_at_full_zeta"], mid=mid, dmc=dmc,
                                         mc=mc, tmr=tmr, tm=tmr * cal["tau_c"], seed=seed),
                           v)
                n += 1
    print(f"G2: {n} runcards under {out}")


# ------------------------------------------------------------- groups G1..G5

HDR_G1 = """20260831 cw_mem_g1 -- group G1: THE PHOTOGRAPH TRACK.

The reference arm. pmem is the MEDIAN of P at zeta_eff(chi = 0.5) = {zg1:.4g}, and mc = 0.5,
which together make the uniform chi == 0.5 state an EXACT fixed point of the loop: f = 0.5
gives m = 0.5 = mc gives chi* = 0.5. The loop therefore has no bias of its own, and
everything chi does comes from the spatial fluctuations of P. This is the arm in which the
memory is a linear time-average of the pressure and the phenotype is a blurred copy of it --
a photograph, not a decision -- and the measurement is how blurred: how std(chi), L_chi and
the autocorrelation time of chi scale with tau_m.

That the threshold differs from the one A3 chose is deliberate and is the whole distinction
between this arm and G3. A3's pmem maximises the contrast f(zeta) - f(0.3 zeta), i.e. it is
chosen to make the loop as nonlinear as possible; this one is chosen to make it as linear as
possible.

THIS CASE: tau_m = {tmr:g} tau_c = {tm:.0f} steps, seed {seed}."""

HDR_G3 = """20260831 cw_mem_g3 -- group G3: THE MAIN FIGURE.

At the coexistence threshold mc* = {mcs:.4f} located by G2, the tau_m grid is run from both
ends: chi == 0 (fully active) and chi == 1 (at the activity floor). The question is whether
a high-activity phase APPEARS and SUSTAINS itself as tau_m grows -- i.e. whether the two
starts stop agreeing, and at what tau_m the disagreement sets in.

The story the campaign is testing lives on this axis. At small tau_m the memory tracks the
instantaneous pressure, chi is a copy of P, and both starts wash into the same state. As
tau_m grows past the flow's own clock the memory averages over many defect encounters, chi
stops being a copy and becomes a slow, coarse variable of its own, and the two starts can
end up in different places -- which is what a phase is.

Each end starts at its own self-consistent memory, m0 = f(zeta_eff(chi_init)); see A4.

THIS CASE: chi == {chi0:g}, m0 = {m0:.4f}, tau_m = {tmr:g} tau_c = {tm:.0f} steps,
seed {seed}."""

HDR_G4 = """20260831 cw_mem_g4 -- group G4: THE STRIPE (FRONT) GROUP.

Half the box at (chi_hi = {chi_hi:.4f}, m = f_hi = {m_hi:.4f}), half at
(chi_lo = {chi_lo:.4f}, m = f_lo = {m_lo:.4f}), split along x -- so two flat interfaces under
the periodic boundary, which check each other. The front direction and speed are read from
the slope of <chi>(t) and from the video; mc is scanned around mc* because the front speed
is what changes sign there, and a speed that crosses zero is a far sharper locator of the
coexistence point than an area fraction is.

THE PHENOTYPE IS FROZEN FOR THE FIRST {frz:.0f} STEPS (5 tau_c). The initial condition
prescribes chi and m but the flow starts from rest with a nearly uniform director, so
without the freeze the front would begin moving while the turbulence that is supposed to
drive it does not yet exist, and the first 5 tau_c of every <chi>(t) slope would be an
artifact of the spin-up. Only the switching SOURCE is frozen; advection and diffusion keep
running, so the interface relaxes to a profile the dynamics can actually produce.

THIS CASE: mc = mc* {dmc:+.2f} = {mc:.4f}, tau_m = {tmr:g} tau_c = {tm:.0f} steps,
seed {seed}."""

HDR_G5 = """20260831 cw_mem_g5{arm} -- group G5: CONTROLS.

Each arm removes one alternative explanation of the G3 result.

  (a) switch-sign = +1. The SAME grid as G3 with the feedback sign flipped: high memory now
      drives chi UP, towards the floor, which is negative feedback. If the two-phase state
      survives a sign flip it is not the feedback that makes it, and the G3 result would be
      an artifact of the initial condition or of the noise.
  (b) tau_chi at 0.1, 1 and 3 tau_c with tau_m fixed at 3 tau_c. The campaign holds
      tau_chi = 0.3 tau_c everywhere, so nothing in G3 distinguishes "tau_m matters" from
      "the ratio tau_m/tau_chi matters". This arm moves the other clock.
  (c) Dbio halved at tau_m = 10 tau_c. The phenotype's structure could be set by its
      diffusive cutoff rather than by the memory; if L_chi is unchanged at half the
      diffusivity it is not the diffusion that sets it.
  (d) chi-width doubled at tau_m = 30 tau_c, from both ends. Halves the loop gain G. If the
      two ends still separate at half the gain, the separation is not sitting on the
      knife-edge of the gain criterion.

THIS CASE: {this}."""


def gen_rest(out, cal, mc_star):
    tc = cal["tau_c"]
    n = 0

    # ---- G1: the photograph track
    for tmr in TAU_M_RATIOS:
        for seed in SEEDS_SINGLE:
            v = base(cal, seed=seed,
                     **{"chi-config": "uniform", "chi0": 0.5, "chi-noise": 0,
                        "chi-length": 0, "m0": 0.5, "mc": 0.5,
                        "pmem": cal["pmem_g1"], "tau-m": round(tmr * tc, 1)})
            write_case(os.path.join(out, "cw_mem_g1", f"tm{tag(tmr,1)}_s{seed}"),
                       HDR_G1.format(zg1=cal["zeta_eff_g1"], tmr=tmr, tm=tmr * tc,
                                     seed=seed), v)
            n += 1

    ends = ((0.0, cal["f_at_full_zeta"]), (1.0, cal["f_at_floor"]))

    def single(mc, tmr, chi0, m0, seed, **extra):
        return base(cal, seed=seed,
                    **{"chi-config": "uniform", "chi0": chi0, "chi-noise": 0,
                       "chi-length": 0, "m0": round(m0, 6), "mc": round(mc, 6),
                       "tau-m": round(tmr * tc, 1), **extra})

    # ---- G3: the main figure
    for tmr in TAU_M_RATIOS:
        for chi0, m0 in ends:
            for seed in SEEDS_SINGLE:
                write_case(os.path.join(out, "cw_mem_g3",
                                        f"tm{tag(tmr,1)}_chi{int(chi0)}_s{seed}"),
                           HDR_G3.format(mcs=mc_star, chi0=chi0, m0=m0, tmr=tmr,
                                         tm=tmr * tc, seed=seed),
                           single(mc_star, tmr, chi0, m0, seed))
                n += 1

    # ---- G4: the stripe group
    frz = int(round(5 * tc))
    for dmc in (-0.04, -0.02, 0.0, 0.02, 0.04):
        for tmr in (10.0, 30.0):
            for seed in SEEDS_STRIPE:
                mc = mc_star + dmc
                v = base(cal, seed=seed,
                         **{"chi-config": "stripe",
                            "chi-lo": round(cal["chi_lo"], 6),
                            "chi-hi": round(cal["chi_hi"], 6),
                            "m-lo": round(cal["f_lo"], 6), "m-hi": round(cal["f_hi"], 6),
                            "chi0": 0.5, "m0": 0.5, "mc": round(mc, 6),
                            "tau-m": round(tmr * tc, 1), "chi-freeze-steps": frz})
                write_case(os.path.join(out, "cw_mem_g4",
                                        f"mc{tag(dmc)}_tm{tag(tmr,1)}_s{seed}"),
                           HDR_G4.format(chi_hi=cal["chi_hi"], chi_lo=cal["chi_lo"],
                                         m_hi=cal["f_hi"], m_lo=cal["f_lo"], frz=frz,
                                         dmc=dmc, mc=mc, tmr=tmr, tm=tmr * tc, seed=seed), v)
                n += 1

    # ---- G5a: the feedback sign flipped
    for tmr in TAU_M_RATIOS:
        for chi0, m0 in ends:
            for seed in SEEDS_SINGLE:
                v = single(mc_star, tmr, chi0, m0, seed, **{"switch-sign": 1})
                write_case(os.path.join(out, "cw_mem_g5a",
                                        f"tm{tag(tmr,1)}_chi{int(chi0)}_s{seed}"),
                           HDR_G5.format(arm="a", this=f"switch-sign = +1, chi == {chi0:g}, "
                                         f"tau_m = {tmr:g} tau_c, seed {seed}"), v)
                n += 1

    # ---- G5b: the other clock
    for tcr in (0.1, 1.0, 3.0):
        for chi0, m0 in ends:
            for seed in SEEDS_SINGLE:
                v = single(mc_star, 3.0, chi0, m0, seed,
                           **{"tau-chi": round(max(tcr * tc, 20.0), 1)})
                write_case(os.path.join(out, "cw_mem_g5b",
                                        f"tc{tag(tcr,1)}_chi{int(chi0)}_s{seed}"),
                           HDR_G5.format(arm="b", this=f"tau_chi = {tcr:g} tau_c = "
                                         f"{max(tcr*tc, 20.0):.0f} steps at tau_m = 3 tau_c, "
                                         f"chi == {chi0:g}, seed {seed}"), v)
                n += 1

    # ---- G5c: half the biomass diffusivity
    for chi0, m0 in ends:
        for seed in SEEDS_SINGLE:
            v = single(mc_star, 10.0, chi0, m0, seed,
                       **{"Dbio": round(cal["Dbio"] / 2, 6)})
            write_case(os.path.join(out, "cw_mem_g5c", f"chi{int(chi0)}_s{seed}"),
                       HDR_G5.format(arm="c", this=f"Dbio = {cal['Dbio']/2:.4g} (half) at "
                                     f"tau_m = 10 tau_c, chi == {chi0:g}, seed {seed}"), v)
            n += 1

    # ---- G5d: half the loop gain
    for chi0, m0 in ends:
        for seed in SEEDS_SINGLE:
            v = single(mc_star, 30.0, chi0, m0, seed,
                       **{"chi-width": round(2 * cal["chi_width"], 5)})
            write_case(os.path.join(out, "cw_mem_g5d", f"chi{int(chi0)}_s{seed}"),
                       HDR_G5.format(arm="d", this=f"chi-width = {2*cal['chi_width']:.4g} "
                                     f"(doubled, G halved) at tau_m = 30 tau_c, "
                                     f"chi == {chi0:g}, seed {seed}"), v)
            n += 1

    print(f"G1/G3/G4/G5: {n} runcards under {out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["ladder", "a4", "g2", "rest"])
    ap.add_argument("outdir")
    ap.add_argument("--zeta", type=float, default=0.15)
    ap.add_argument("--tau-c-pred", type=float, default=540.0)
    ap.add_argument("--calib", default=None)
    ap.add_argument("--mc-star", type=float, default=None,
                    help="the coexistence threshold located by G2 (stage 'rest')")
    ap.add_argument("--g-list", type=float, nargs="*", default=None,
                    help="loop gains to run in stage a4; two runs (chi==0, chi==1) each")
    ap.add_argument("--wchi-list", type=float, nargs="*", default=None,
                    help="chi-widths to run in stage a4; overrides --g-list")
    a = ap.parse_args()
    cal = json.load(open(a.calib)) if a.calib else None
    if a.stage == "ladder":
        gen_ladder(a.outdir, a.zeta, a.tau_c_pred)
    elif a.stage == "a4":
        gen_a4(a.outdir, cal, a.g_list, a.wchi_list)
    elif a.stage == "g2":
        gen_g2(a.outdir, cal)
    else:
        if a.mc_star is None:
            raise SystemExit("stage 'rest' needs --mc-star from the G2 scan")
        gen_rest(a.outdir, cal, a.mc_star)
