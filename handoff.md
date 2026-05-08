# Handoff Summary

## Working Mode

- Discuss and confirm the physics/code plan before implementation.
- Only modify code; do not compile or run unless explicitly requested.
- User validates compilation and simulation.
- Main working branch: `pheno`.
- Related branch: `pheno_nemiso`, where Q order preference is coupled to `phi`.

## Model Base

- New model: `go-or-grow`.
- Based on `LyotropicWithDivision`.
- Keeps original `Q`, `phi`, `u` dynamics, LB hydrodynamics, and division mask logic.
- `phi` is cell/nematic density.
- `Q` is a 2D traceless nematic tensor:

```text
Q = [[Qxx, Qyx],
     [Qyx, -Qxx]]

Q_ab = 2S(n_a n_b - delta_ab/2)
1/2 Tr(Q^2) = Qxx^2 + Qyx^2
```

## Implemented Changes

- Added `src/models/go-or-grow.hpp` and `src/models/go-or-grow.cpp`.
- Registered `model = go-or-grow`.
- `circle` initialization is now a solid disk.
- Initial nematic director is controlled by:

```ini
angle
noise
```

- `angle = 0`: horizontal/x direction.
- `angle = 90`: vertical/y direction.

## Phenotype Variables

- Do not evolve `chi` independently.
- Evolve grow-type density:

```text
m = phi * chi
```

- `chi` is diagnostic/derived:

```text
chi = m/phi  if phi > eps
chi = 0      if phi <= eps
```

- Enforce consistency every step:

```text
0 <= m <= phi
```

## Phenotype Initialization

New parameters:

```ini
chi-config = noise/front
chi0
chi-noise
chi-length
```

- `chi-config = noise`: correlated Gaussian noise initialization.
- `chi-config = front`: left side `chi ~= 1`, right side `chi ~= 0`, with tanh transition.
- In `front` mode, `chi-length` controls transition width.

## Free Energy On `pheno`

```text
F = integral [f_phi + f_Q] dA
```

```text
f_phi =
  1/2 AA phi^2(1-phi)^2
+ 1/2 KK |grad phi|^2
+ 1/2 B (phi - phi_c)^2 Theta(phi - phi_c)
```

```text
f_Q =
  CC/2 (Snem - 1/2 Tr(Q^2))^2
+ 1/2 LL |grad Q|^2
```

Relevant parameters:

```ini
B
Snem
phi-critical
```

## Free Energy On `pheno_nemiso`

The Q bulk term was changed to prefer nematic order only where `phi` is high:

```text
f_Q =
  CC/2 (Snem*phi - 1/2 Tr(Q^2))^2
+ 1/2 LL |grad Q|^2
```

Corresponding chemical-potential contribution:

```text
mu_Q = CC*Snem*(Snem*phi - q2)
q2 = Qxx^2 + Qyx^2
```

## Active Stress

Current active stress:

```text
sigma_active = - zeta * phi * (1 - chi) * Q
```

- `zeta > 0` is extensile.
- `chi = 0` go-like region: high activity.
- `chi = 1` grow-like region: activity suppressed.

## Proliferation

Original division mask mechanism is retained:

```text
division_rate * phi * (1 - phi/phi_critical)
```

selects division centers. Each center creates a patch of radius:

```ini
division-radius
```

lasting:

```ini
division-time
```

Current `phi` source:

```text
R_phi = chi * phi * (alpha*division_mask - beta*death_mask)
```

Current `m` source:

```text
R_m = chi * R_phi
```

## Current m Dynamics

```text
partial_t m =
  - div(m u)
+ GammaP div(chi grad mu)
+ Dchi div(phi grad chi)
+ chi R_phi
```

New phenotype diffusion parameter:

```ini
Dchi
```

- Default: `Dchi = 0`.
- If omitted, old behavior is preserved.

