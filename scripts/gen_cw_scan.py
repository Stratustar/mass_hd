#!/usr/bin/env python3
"""Emit the runcards for the 20260830 tau_m/tau_chi scan at fixed zeta.

Everything except tau_m, tau_chi, pmem and switch-sign is frozen. Those four are set here
from numbers MEASURED by the open-loop pre-run (20260830/cw_pre), which is why this is a
script and not 248 hand-written files: tau_motion and the pressure percentiles are not known
until that run exists, and hand-copying them 248 times is how a scan ends up with two
different values of a "fixed" constant.

Grid:
  group 0   tau_m/tau_motion = tau_chi/tau_motion = 0.2, one baseline point
  group 1   both ratios drawn from RATIOS, keeping them within a factor 3 of each other
  x switch-sign in {+1, -1}
  x pmem at the 75th and 90th percentile of P -- so only the top 25% / top 10% of the area is
    above threshold and able to accumulate memory
  pmem-width = IQR(P)/4 throughout: the transition then spans about half an IQR, sharp enough
    to mean "only that fraction accumulates" without becoming a step that has no gradient for
    the memory source to follow.

Usage:
  gen_cw_scan.py --tau-motion T --p75 X --p90 Y --iqr Z [--outdir cases/20260830/cw_scan]
"""
import argparse
import os

RATIOS = [0.5, 2, 4, 6, 8, 10, 15, 20, 25, 30, 40]
MAX_SPREAD = 3.0
GROUP0 = (0.2, 0.2)

# frozen since 20260827 cw_regime; see that runcard for the measurement behind each value
FROZEN = dict(Gamma=0.05, xi=0.4, CC=0.03, LL=0.24, tau=2.0, rho=40, friction=0.0,
              zeta=0.54, chi0=0.5, mc=0.5, m0=0.5, chi_width=0.15, Dbio=0.03)
L, NSTEPS, NINFO, SEED = 800, 120000, 1000, 1001


def tag(x):
    return f"{x:g}".replace(".", "p")


def pairs():
    out = [(GROUP0, "g0")]
    for a in RATIOS:
        for b in RATIOS:
            if 1.0 / MAX_SPREAD - 1e-9 <= a / b <= MAX_SPREAD + 1e-9:
                out.append(((a, b), "g1"))
    return out


HEADER = """# 20260830 cw_scan: the tau_m / tau_chi scan at fixed activity.
#
# zeta = 0.54 (A = zeta_eff/CC = 9) and every other parameter is frozen; the only things that
# vary across the 248 runs are tau_m, tau_chi, pmem and switch-sign.
#
# THE TWO AXES ARE RATIOS TO tau_motion = l_d / u_rms = {tmot:.0f} steps, the time the flow
# takes to carry material across one defect spacing, measured open-loop at this activity by
# 20260830/cw_pre. Along a material trajectory the pressure history is a sequence of defect
# encounters and tau_motion is the interval between them, so
#     tau_m   / tau_motion = how many encounters the memory averages over   = {rm:g}
#     tau_chi / tau_motion = how far the phenotype lags behind the flow     = {rc:g}
# The two are kept within a factor {spread:g} of each other; a phenotype relaxing decades away
# from its own memory is a different problem from the one this scan is about.
#
# PMEM IS A PERCENTILE OF THE MEASURED PRESSURE DISTRIBUTION, not a guess: at pmem = P{pct}
# only the top {top}% of the area is above threshold and able to accumulate memory. Measured
# open-loop: P75 = {p75:.5e}, P90 = {p90:.5e}, IQR = {iqr:.5e}.
# pmem-width = IQR/4 = {pw:.5e}: the tanh transition then spans about half an IQR -- sharp
# enough to mean "only that fraction accumulates", not so sharp that g(P) becomes a step and
# leaves the memory source with no gradient to follow.
#
# THIS CASE: {grp}, switch-sign = {s:+d}, pmem = P{pct} ({top}% of area above threshold),
#            tau_m = {taum:.0f}, tau_chi = {tauc:.0f} steps.
#
# FROZEN (unchanged from 20260827 cw_regime / 20260828 cw_loop; each pinned by measurement):
#   xi = 0.4    flow-tumbling; xi >= 1 is numerically unreachable in this scheme
#   tau = 2.0, rho = 40, Gamma = 0.05, friction = 0
#   CC = 0.03, LL = 0.24  ->  xi_N = 2.0, eta = 20
#   Dbio = 0.03  l_B = sqrt(Dbio*tau_motion) is a fixed fraction ~0.18 of the defect spacing
#                at every activity, and stays above 3 lattice cells
#   chi0 = mc = 0.5   with pmem at a percentile the uniform state is NOT a fixed point here,
#                     unlike the 20260828 scan where pmem = median(P) made it one -- that is
#                     deliberate: the memory should be driven by the top tail of P.
#
# ninfo = {ninfo} for EVERY case, so all videos share one sampling interval (dashboard rule,
# 2026-08-28); at L = {L} that is {nframes} frames per run.
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tau-motion", type=float, required=True)
    ap.add_argument("--p75", type=float, required=True)
    ap.add_argument("--p90", type=float, required=True)
    ap.add_argument("--iqr", type=float, required=True)
    ap.add_argument("--outdir", default="cases/20260830/cw_scan")
    a = ap.parse_args()

    pw = a.iqr / 4.0
    n = 0
    for (rm, rc), grp in pairs():
        taum, tauc = rm * a.tau_motion, rc * a.tau_motion
        for s, stag in ((1, "sp"), (-1, "sm")):
            for pct, top, pmem in ((75, 25, a.p75), (90, 10, a.p90)):
                name = f"{stag}_p{pct}_tm{tag(rm)}_tc{tag(rc)}"
                d = os.path.join(a.outdir, name)
                os.makedirs(d, exist_ok=True)
                head = HEADER.format(tmot=a.tau_motion, rm=rm, rc=rc, spread=MAX_SPREAD,
                                     pct=pct, top=top, p75=a.p75, p90=a.p90, iqr=a.iqr,
                                     pw=pw, grp=grp, s=s, taum=taum, tauc=tauc,
                                     ninfo=NINFO, L=L, nframes=NSTEPS // NINFO)
                f = FROZEN
                with open(os.path.join(d, "run.dat"), "w") as fh:
                    fh.write(head + f"""model           = confluent-wet
nsteps          = {NSTEPS}
nstart          = 0
ninfo           = {NINFO}
LX              = {L}
LY              = {L}
bc              = 0
seed            = {SEED}

Gamma           = {f['Gamma']}
xi              = {f['xi']}
CC              = {f['CC']}
LL              = {f['LL']}
tau             = {f['tau']}
rho             = {f['rho']}
friction        = {f['friction']}
zeta            = {f['zeta']}

angle           = 0
noise           = 0.05
initial-order   = 1
director-config = uniform
defect-sep      = 0

Dbio            = {f['Dbio']}
chi-config      = uniform
chi0            = {f['chi0']}
chi-noise       = 0
chi-length      = 0
tau-chi         = {tauc:.1f}
chi-width       = {f['chi_width']}
switch-sign     = {s}
mc              = {f['mc']}
m0              = {f['m0']}
tau-m           = {taum:.1f}
pmem            = {pmem:.6e}
pmem-width      = {pw:.6e}
""")
                n += 1
    print(f"wrote {n} runcards under {a.outdir}")
    print(f"  tau_motion = {a.tau_motion:.1f}   pmem = {a.p75:.4e} (P75) / {a.p90:.4e} (P90)"
          f"   pmem-width = {pw:.4e}")
    print(f"  {len(pairs())} (ratio pairs incl. group 0) x 2 signs x 2 pmem = {n}")


if __name__ == "__main__":
    main()
