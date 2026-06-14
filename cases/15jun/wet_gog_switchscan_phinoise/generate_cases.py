#!/usr/bin/env python3
"""Generate the 15jun WET switching scan with correlated-noise init (48 cases).

Wet (full LB) follow-up to the local wet switching tests. Baseline = the local
wet case (alpha=0.0001, chi correlated-noise, phi-noise=0.35), scaled up to a
400x400 domain so the colony does not fill the box.

Axes:
  director : random (noise=1) | aligned (noise=0, angle=0)
  Achi     : 0 | 0.005 | 0.01 | 0.02   (0 = reversible/smooth chi; larger binarizes)
  Ochi     : 0.01 | 0.02 | 0.05
  pswitch  : 0.05 | 0.1                 (alpha=0.0001 center P only reaches ~0.11)

Fixed: model=go-or-grow (WET), alpha=0.0001, chi correlated-noise (chi0=0.5,
       chi-noise=0.04, chi-length=5), phi-noise=0.35, phi-length=5, Dchi=0,
       Snem=1, zeta=0.05, B=1, phi_c=1.0, 400x400 radius 20, 40000 steps.

2 x 4 x 3 x 2 = 48 cases. One run.dat per leaf. Re-run to regenerate.
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))

TEMPLATE = """\
# 15jun WET switching scan, correlated-noise init. model = go-or-grow (full LB).
# 400x400 radius 20, 40000 steps, alpha=0.0001, chi correlated-noise, phi-noise=0.35.
# Axes: director x Achi x Ochi x pswitch.

model    = go-or-grow
config   = circle
nsteps   = 40000
nstart   = 0000
ninfo    = 100
LX       = 400
LY       = 400
bc       = 0
seed     = 1001
GammaQ   = 0.3
GammaP   = 0.05
zeta     = 0.05
zetaI    = 0
friction = 10
tauNem   = 1.0
tauIso   = 1.0
AA       = 0.03
CC       = 0.1
LL       = 0.002
KK       = 0.06
rho      = 40
xi       = 0.3
B        = 1.0
Snem     = 1
level    = 50
radius   = 20
conc     = 1
noise    = {noise}
initial-order = 1
angle    = 0
relax-steps = 500
relax-dt    = 1.0
relax-phi   = 1
relax-Q     = 1
chi-config = noise
chi0       = 0.5
chi-noise  = 0.04
chi-length = 5
phi-noise  = 0.35
phi-length = 5
Dchi       = 0
Achi      = {Achi}
Ochi      = {Ochi}
pswitch   = {pswitch}
growTogether = 1
beta         = 0
death-rate   = 0
death-time   = 1
death-radius = 0
division-rate   = 1e-1
division-time   = 2
alpha           = 0.0001
phi-critical    = 1.0
division-phi-critical-factor = 0
division-radius = 3
surface-stress = 0
"""

DIRECTOR = [("dirrandom", "1"), ("diraligned", "0")]
ACHI = [("0", "0"), ("0p005", "0.005"), ("0p01", "0.01"), ("0p02", "0.02")]
OCHI = [("0p01", "0.01"), ("0p02", "0.02"), ("0p05", "0.05")]
PSWITCH = [("0p05", "0.05"), ("0p1", "0.1")]


def main():
    n = 0
    for dir_label, noise in DIRECTOR:
        for achi_label, achi in ACHI:
            for ochi_label, ochi in OCHI:
                for ps_label, ps in PSWITCH:
                    variant = f"{dir_label}_Achi{achi_label}_Ochi{ochi_label}_ps{ps_label}"
                    leaf = os.path.join(HERE, variant)
                    os.makedirs(leaf, exist_ok=True)
                    text = TEMPLATE.format(noise=noise, Achi=achi, Ochi=ochi, pswitch=ps)
                    with open(os.path.join(leaf, "run.dat"), "w") as f:
                        f.write(text)
                    n += 1
    print(f"wrote {n} run.dat files under {HERE}")


if __name__ == "__main__":
    main()
