#!/usr/bin/env python3
"""Generate the confluent-wet closed-loop production grid.

THE DESIGN, and why each axis is what it is.

ACTIVITY  A = zeta_eff/CC in {1, 3, 9}.  A is the only control parameter for defect density
    (A = energy density of the activity / cost of melting the nematic order), and this ladder
    is geometric in the ONE material clock the model has:
        t_eddy ~ 85/zeta_eff  ->  2833 / 944 / 315 steps, a factor 3 apart each.
    That is the whole reason for three activities rather than two.  With two regimes at A = 1
    and 10 the accessible windows of tau/t_eddy were [0.06, 2.5] and [0.62, 26] -- nearly
    disjoint, so the two families could not be put on a common axis.  With a factor-3 ladder
    ONE absolute tau ladder slides across all three activities and the windows overlap.

TAU LADDER  tau_m in {t_eddy(A=9), t_eddy(A=3), t_eddy(A=1)}, i.e. the measured clocks
    themselves, used as absolute step counts for every activity.  Each activity then sees
    ratios tau_m/t_eddy of {0.11, 0.33, 1}, {0.33, 1, 3} and {1, 3, 9} respectively -- the
    union covers 0.11 to 9 and ratio 1 is realised in all three, which is what makes the
    "does the physics collapse on the ratio?" question answerable.
    The floor is real: tau_m must exceed the order-magnitude relaxation 1/(4*Gamma*CC) = 167
    steps, below which the memory records |Q| transients and sigma_bulk -- part of P -- is a
    function of |Q|.  The smallest rung is above it.

PHENOTYPE LAG  tau_chi/tau_m in {1, 3}: how much further the phenotype lags behind the memory.

SWITCH SIGN  both.  s = +1 is the negative-feedback branch (the dry line's self-limiting one,
    where the wet dbio batch self-organised onto the instability threshold); s = -1 is the
    escape branch.  They are different physics, not a sign convention.

DBIO  {0, 0.03}.  0.03 is not arbitrary: l_B = sqrt(Dbio*t_eddy) = 9.2 / 5.3 / 3.1 against a
    defect spacing d = 49 / 30 / 18, so l_B/d = 0.19 / 0.18 / 0.17 -- the SAME fraction of the
    defect spacing at all three activities, while staying above 3 lattice units (resolved) and
    far below the 0.2 explicit-diffusion stability limit.  Dbio = 0 is the control.

    3 activities x 2 signs x 2 Dbio x 3 tau_m x 2 ratios = 72 runs.

RUN LENGTH AND OUTPUT CADENCE, per case, so that no run is longer or fatter than its own
physics requires:
    t_grow  = 5*eta/zeta_eff = 100/zeta_eff        (eta = rho*(tau-1/2)/3 = 20)
    T_trans = max(13*t_grow, 5*tau_chi)            hydrodynamic saturation AND loop equilibration
    T_fast  = min(tau_m, t_eddy)                   the fastest thing the lag axis must resolve
    T_slow  = max(tau_chi, t_eddy)                 the slowest thing the window must span
    ninfo   = max(T_fast/5, T_slow/36)             5 samples per fast clock, if affordable
    nstart  = T_trans                              <- OUTPUT STARTS AFTER THE TRANSIENT
    nsteps  = nstart + 180*ninfo                   180 steady frames, spanning >= 5*T_slow

nstart is what makes the grid affordable at all: A = 1 needs a 43000-step hydrodynamic
transient, which at the cadence the lag axis wants would be ~700 frames of pure warm-up.
Writing nothing until t >= nstart leaves 181 frames per run (~7.6 GB at L = 400) instead of
~870 (~36 GB).

L = 400, not 800.  The correlation lengths to be measured are d = 49 at worst and L_P ~ 58 at
L = 800, so L = 400 leaves box/L_P ~ 7 and 66 defects at the sparsest activity -- and it is
the box the dry line's own regime scan used for the same reason.  A handful of L = 800
confirmations belong afterwards, driven by whichever corner turns out to matter.

Usage:
  generate_cases.py [--calib results.json] [--outdir .] [--dry-run]
The calibration JSON is 20260828/cw_pmem's cw_calib.json: it supplies pmem = median(P),
pmem-width = IQR(P)/2 and t_eddy per activity.  Without it the predicted values below are
used and the runcards say so.
"""
import argparse
import json
import math
import os

# ----------------------------------------------------------------- frozen baseline
FROZEN = dict(Gamma=0.05, xi=0.4, CC=0.03, LL=0.24, tau=2.0, rho=40, friction=0.0,
              chi0=0.5, mc=0.5, m0=0.5, chi_width=0.15, noise=0.05, seed=1001, L=400)
ETA = FROZEN["rho"] * (FROZEN["tau"] - 0.5) / 3.0            # = 20
TAU_M_FLOOR = 1.0 / (4 * FROZEN["Gamma"] * FROZEN["CC"])      # = 167 steps

ACTIVITIES = [1, 3, 9]
SIGNS = [(1, "sp"), (-1, "sm")]
DBIOS = [(0.0, "D0"), (0.03, "D0p03")]
TC_RATIOS = [1, 3]
NSTEADY = 180

# fallbacks if no calibration JSON is supplied (t_eddy ~ 85/zeta_eff, sigma_P ~ 0.76*zeta_eff)
PREDICTED = {A: dict(t_eddy=85.0 / (A * FROZEN["CC"]), pmem=0.0,
                     pmem_width=0.674 * 0.76 * A * FROZEN["CC"] / 2)
             for A in ACTIVITIES}


def nice(x):
    """Round a cadence to a readable number without moving it more than a few percent."""
    if x <= 0:
        return 1
    mag = 10 ** math.floor(math.log10(x))
    for step in (mag, mag / 2, mag / 5, mag / 10):
        if x / step >= 3:
            return int(max(1, round(x / step) * step))
    return int(max(1, round(x)))


def load_calib(path):
    if not path or not os.path.exists(path):
        return None
    with open(path) as f:
        rows = json.load(f)
    out = {}
    for r in rows:
        A = int(round(r["A"]))
        out[A] = dict(t_eddy=r["t_eddy"], pmem=r["pmem"],
                      pmem_width=r["pmem_width_iqr2"], d=r.get("d"),
                      u_rms=r.get("u_rms"), sigma_P=r.get("sigma_P"),
                      tau_motion=r.get("tau_motion"), N_def=r.get("N_def"),
                      saturated=r.get("saturation", {}).get("ok"))
    return out


def plan(cal):
    """The full grid as a list of dicts, with every derived number computed once."""
    t_eddy = {A: cal[A]["t_eddy"] for A in ACTIVITIES}
    ladder = [nice(t_eddy[A]) for A in sorted(ACTIVITIES, reverse=True)]   # A9, A3, A1
    rows = []
    for A in ACTIVITIES:
        ze = A * FROZEN["CC"]
        t_grow = 5 * ETA / ze
        for tau_m in ladder:
            for ratio in TC_RATIOS:
                tau_chi = tau_m * ratio
                T_trans = max(13 * t_grow, 5 * tau_chi)
                T_fast = min(tau_m, t_eddy[A])
                T_slow = max(tau_chi, t_eddy[A])
                ninfo = nice(max(T_fast / 5.0, T_slow / 36.0))
                nstart = int(math.ceil(T_trans / ninfo) * ninfo)
                nsteps = nstart + NSTEADY * ninfo
                for s, stag in SIGNS:
                    for D, dtag in DBIOS:
                        rows.append(dict(
                            A=A, zeta_eff=ze, zeta=ze / (1 - FROZEN["chi0"]),
                            switch_sign=s, Dbio=D, tau_m=tau_m, ratio=ratio,
                            tau_chi=tau_chi, t_eddy=t_eddy[A], t_grow=t_grow,
                            ninfo=ninfo, nstart=nstart, nsteps=nsteps,
                            frames=NSTEADY + 1,
                            window=NSTEADY * ninfo,
                            pmem=cal[A]["pmem"], pmem_width=cal[A]["pmem_width"],
                            name=f"A{A}_{stag}_{dtag}_tm{tau_m}_tc{ratio}"))
    return rows, ladder


HEADER = """# 20260828 cw_loop: the confluent-wet closed-loop (tau_m, tau_chi) production grid.
#
# THIS CASE
#   A = {A}  (zeta_eff = {zeta_eff:g}, zeta = {zeta:g} at chi0 = 0.5)   t_eddy = {t_eddy:.0f} steps
#   tau_m   = {tau_m}  =  {rm:.2f} x t_eddy        switch-sign = {switch_sign:+d}   Dbio = {Dbio:g}
#   tau_chi = {tau_chi}  =  {ratio} x tau_m  =  {rc:.2f} x t_eddy
#
# WHY THESE NUMBERS.  A in {{1,3,9}} is geometric in the model's single material clock
# (t_eddy ~ 85/zeta_eff -> 2833/944/315), so ONE absolute tau ladder slides across all three
# activities and their ratio windows overlap -- which two regimes at A = 1 and 10 could not do
# (windows [0.06,2.5] and [0.62,26], nearly disjoint).  The ladder IS the three measured
# t_eddy values, so this case sits at ratio {rm:.2f} for the memory and {rc:.2f} for the phenotype.
# tau_m = {tau_m} clears the order-magnitude relaxation floor 1/(4*Gamma*CC) = {floor:.0f} steps, below
# which the memory would record |Q| transients (and sigma_bulk, part of P, is a function of |Q|).
# Dbio = {Dbio:g}: l_B = sqrt(Dbio*t_eddy) = {lB} against a defect spacing d ~ {dsp:.0f}, i.e. the same
# fraction of the defect spacing at every activity, resolved (>3 cells) and far below the 0.2
# explicit-diffusion limit.
#
# THE LOOP, CLOSED
#   D_t m   = (g(P) - m)/tau_m          g(P)    = .5*(1 + tanh((P - pmem)/pmem-width))
#   D_t chi = (chi*(m) - chi)/tau_chi   chi*(m) = .5*(1 + s*tanh((m - mc)/chi-width))
#   zeta_eff = zeta*(1 - chi)
# pmem = {pmem:.6g} is the MEASURED median(P) of the open-loop run at this activity
# (20260828/cw_pmem/A{A}) and pmem-width = {pmem_width:.6g} its measured IQR/2.  With chi0 = mc = 0.5
# this makes the uniform state an EXACT fixed point, so the loop is driven only by the spatial
# fluctuations of P and not by a misplaced threshold -- the error that drove the 20260826 batch.
#
# LENGTH AND CADENCE, from this case's own clocks
#   t_grow  = 5*eta/zeta_eff = {t_grow:.0f}      T_trans = max(13*t_grow, 5*tau_chi) = {nstart}
#   ninfo   = max(min(tau_m,t_eddy)/5, max(tau_chi,t_eddy)/36) = {ninfo}
#   nstart  = {nstart}: NOTHING IS WRITTEN before the transient ends.  Without it A = 1 would
#             spend ~700 frames of a ~870-frame archive warming up.
#   nsteps  = nstart + {nsteady}*ninfo = {nsteps}, i.e. {frames} steady frames spanning {window} steps
#             = {wspan:.1f} x max(tau_chi, t_eddy).
#
# FROZEN PARAMETER SET, unchanged since 20260827 cw_regime -- xi = 0.4 (the flow-aligning
# branch xi >= 1 is numerically unreachable in this scheme), tau = 2 (tau = 1 blows the Mach
# budget before A = 4; tau = 2 also puts Re in the Stokes range a cell layer lives in),
# CC = 0.03 and LL = 0.24 giving xi_N = 2.0 (measured sweet spot: core is one cell at 1.0,
# cores merge into melted walls at 3.0), Gamma = 0.05, rho = 40, friction = 0.
model           = confluent-wet
nsteps          = {nsteps}
nstart          = {nstart}
ninfo           = {ninfo}
LX              = {L}
LY              = {L}
bc              = 0
seed            = {seed}

Gamma           = {Gamma}
xi              = {xi}
CC              = {CC}
LL              = {LL}
tau             = {tau}
rho             = {rho}
friction        = {friction}
zeta            = {zeta:g}

angle           = 0
noise           = {noise}
initial-order   = 1
director-config = uniform
defect-sep      = 0

Dbio            = {Dbio:g}
chi-config      = uniform
chi0            = {chi0}
chi-noise       = 0
chi-length      = 0
tau-chi         = {tau_chi}
chi-width       = {chi_width}
switch-sign     = {switch_sign}
mc              = {mc}
m0              = {m0}
tau-m           = {tau_m}
pmem            = {pmem:.6g}
pmem-width      = {pmem_width:.6g}
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calib", default=None)
    ap.add_argument("--outdir", default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    cal = load_calib(a.calib)
    if cal is None or any(A not in cal for A in ACTIVITIES):
        print("!! no calibration JSON -- using PREDICTED t_eddy/pmem; runcards will be wrong")
        cal = PREDICTED
    rows, ladder = plan(cal)

    print(f"tau_m ladder (= measured t_eddy of A9, A3, A1): {ladder}")
    print(f"{'case':30s} {'ninfo':>6s} {'nstart':>8s} {'nsteps':>8s} {'win/Tslow':>9s} "
          f"{'GB':>6s} {'sim_min':>8s}")
    tot_gb = tot_min = 0.0
    seen = set()
    for r in rows:
        # 168 MB/frame at L=800 (measured) -> scales with L^2
        gb = r["frames"] * 0.168 * (FROZEN["L"] / 800.0) ** 2
        # 31 min for 70000 steps at L=800 on 32 cores (measured) -> scales with L^2 * nsteps
        mn = 31.0 * (FROZEN["L"] / 800.0) ** 2 * (r["nsteps"] / 70000.0)
        tot_gb += gb; tot_min += mn
        key = (r["A"], r["tau_m"], r["ratio"])
        if key not in seen:                      # one line per distinct (A, tau) cell
            seen.add(key)
            print(f"{r['name']:30s} {r['ninfo']:6d} {r['nstart']:8d} {r['nsteps']:8d} "
                  f"{r['window']/max(r['tau_chi'], r['t_eddy']):9.1f} {gb:6.1f} {mn:8.1f}")
    print(f"\n{len(rows)} runs   total {tot_gb:.0f} GB   {tot_min/60:.1f} core-hours-equivalent "
          f"(32 cores each, wall time = this / concurrency)")

    if a.dry_run:
        return
    for r in rows:
        d = os.path.join(a.outdir, r["name"])
        os.makedirs(d, exist_ok=True)
        lB = "/".join(f"{math.sqrt(r['Dbio'] * cal[A]['t_eddy']):.1f}" for A in ACTIVITIES) \
            if r["Dbio"] > 0 else "0 (control)"
        with open(os.path.join(d, "run.dat"), "w") as f:
            f.write(HEADER.format(
                rm=r["tau_m"] / r["t_eddy"], rc=r["tau_chi"] / r["t_eddy"],
                floor=TAU_M_FLOOR, lB=lB, dsp=49 * r["A"] ** -0.45,
                nsteady=NSTEADY, wspan=r["window"] / max(r["tau_chi"], r["t_eddy"]),
                **FROZEN, **r))
    print(f"wrote {len(rows)} runcards under {a.outdir}")


if __name__ == "__main__":
    main()
