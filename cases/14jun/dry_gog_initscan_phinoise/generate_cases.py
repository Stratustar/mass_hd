#!/usr/bin/env python3
"""Generate the 14jun correlated-noise-phi init scan (24 cases).

Axes:
  chi      : correlated-noise (chi0=0.5, chi-noise=0.04, chi-length=5) | all-grow (chi0=1, chi-noise=0)
  director : random (noise=1) | aligned (noise=0)
  Achi     : 0 | 0.01 | 0.05          (bistable-switching barrier, the "B" axis)
  alpha    : 0.0002 | 0.0001          (baseline growth and half)

Fixed: phi correlated-noise std=0.35, phi-length=5; Ochi=0.02; pswitch=0.2;
       400x400, radius 20, 20000 steps. Everything else from the 12jun baseline.

One run.dat per leaf <variant>/run.dat. Re-run to regenerate.
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))

TEMPLATE = """\
# 14jun correlated-noise-phi init scan (B+C): symmetry-breaking initial condition.
# Baseline = cases/12jun/dry_gog_switching_big (Snem=1, zeta=0.05), 400x400, radius 20.
# Axes: chi x director x Achi x alpha. Fixed Ochi=0.02, pswitch=0.2, phi-noise=0.35.

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
Dchi       = 0
Achi      = {Achi}
Ochi      = 0.02
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
ACHI = [
    ("0",    "0"),
    ("0p01", "0.01"),
    ("0p05", "0.05"),
]
ALPHA = [
    ("0p0002", "0.0002"),
    ("0p0001", "0.0001"),
]


def main():
    n = 0
    for chi_label, chi0, chi_noise, chi_length in CHI:
        for dir_label, noise in DIRECTOR:
            for achi_label, achi in ACHI:
                for alpha_label, alpha in ALPHA:
                    variant = f"{chi_label}_{dir_label}_Achi{achi_label}_alpha{alpha_label}"
                    leaf = os.path.join(HERE, variant)
                    os.makedirs(leaf, exist_ok=True)
                    text = TEMPLATE.format(
                        noise=noise, chi0=chi0, chi_noise=chi_noise,
                        chi_length=chi_length, Achi=achi, alpha=alpha,
                    )
                    with open(os.path.join(leaf, "run.dat"), "w") as f:
                        f.write(text)
                    n += 1
    print(f"wrote {n} run.dat files under {HERE}")


if __name__ == "__main__":
    main()
