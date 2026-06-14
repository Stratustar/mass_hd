#!/usr/bin/env python3
"""Generate the 15jun switching-parameter scan on the phi-noise init (72 cases).

Follow-up to 14jun/dry_gog_initscan_phinoise: same correlated-noise-phi init and
domain, now scanning the switching/transport knobs with phenotype diffusion on.

Axes:
  chi      : correlated-noise (chi0=0.5, chi-noise=0.04, chi-length=5) | all-grow (chi0=1)
  director : random (noise=1) | aligned (noise=0)
  alpha    : 0.0002 | 0.0001
  Achi     : 0.1 | 0.2 | 0.3          (stronger bistable barrier, the "B" axis)
  Ochi     : 0.005 | 0.01 | 0.02      (pressure-coupling sensitivity)

Fixed: Dchi=0.0005 (phenotype diffusion on, vs 0 in 14jun); phi-noise=0.35,
       phi-length=5; pswitch=0.2; 400x400, radius 20, 20000 steps.

2 x 2 x 2 x 3 x 3 = 72 cases. One run.dat per leaf. Re-run to regenerate.
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))

TEMPLATE = """\
# 15jun switching-parameter scan on the phi-noise init (B+C follow-up).
# Same init/domain as 14jun/dry_gog_initscan_phinoise; Dchi=0.0005 (diffusion on).
# Axes: chi x director x alpha x Achi x Ochi. Fixed pswitch=0.2, phi-noise=0.35.

model    = dry-go-or-grow
config   = circle
nsteps   = 20000
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
chi0       = {chi0}
chi-noise  = {chi_noise}
chi-length = {chi_length}
phi-noise  = 0.35
phi-length = 5
Dchi       = 0.0005
Achi      = {Achi}
Ochi      = {Ochi}
pswitch   = 0.2
growTogether = 1
beta         = 0
death-rate   = 0
death-time   = 1
death-radius = 0
division-rate   = 1e-1
division-time   = 2
alpha           = {alpha}
phi-critical    = 1.0
division-phi-critical-factor = 0
division-radius = 3
surface-stress = 0
"""

# (label, chi0, chi_noise, chi_length)
CHI = [
    ("chinoise", "0.5", "0.04", "5"),
    ("allgrow",  "1",   "0",    "0"),
]
# (label, noise)
DIRECTOR = [
    ("dirrandom",  "1"),
    ("diraligned", "0"),
]
# (label-suffix, value)
ALPHA = [
    ("0p0002", "0.0002"),
    ("0p0001", "0.0001"),
]
ACHI = [
    ("0p1", "0.1"),
    ("0p2", "0.2"),
    ("0p3", "0.3"),
]
OCHI = [
    ("0p005", "0.005"),
    ("0p01",  "0.01"),
    ("0p02",  "0.02"),
]


def main():
    n = 0
    for chi_label, chi0, chi_noise, chi_length in CHI:
        for dir_label, noise in DIRECTOR:
            for alpha_label, alpha in ALPHA:
                for achi_label, achi in ACHI:
                    for ochi_label, ochi in OCHI:
                        variant = (f"{chi_label}_{dir_label}_Achi{achi_label}"
                                   f"_Ochi{ochi_label}_alpha{alpha_label}")
                        leaf = os.path.join(HERE, variant)
                        os.makedirs(leaf, exist_ok=True)
                        text = TEMPLATE.format(
                            noise=noise, chi0=chi0, chi_noise=chi_noise,
                            chi_length=chi_length, Achi=achi, Ochi=ochi, alpha=alpha,
                        )
                        with open(os.path.join(leaf, "run.dat"), "w") as f:
                            f.write(text)
                        n += 1
    print(f"wrote {n} run.dat files under {HERE}")


if __name__ == "__main__":
    main()
