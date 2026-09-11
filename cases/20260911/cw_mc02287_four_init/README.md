# mc=0.2287: four-initialization memory scan

L=256, 20 linear memory times from 0.3 to 18 reference tau_c, four starts, two paired seeds: 160 runs.

Starts at t=0: all 0, all 1, left 0/right 1 stripe, correlated Gaussian noise around 0.5. Noise std=0.1, smoothing length=7.02 lattice units. Noise m starts uniform at 0.2503; stripe memory matches its two initial activity branches.

The phase-scan preparation protocol is retained. Only the chi reaction is frozen; advection/diffusion can mix stripe/noise before feedback release. The analysis records the actual release mean and spatial standard deviation.

After preparation, every run observes 2000 tau_c (1,348,658 steps). Preparation lasts 100–180 tau_c; total 1,416,091–1,470,037 steps. Tail: 500 tau_c. All physical coefficients and output cadence are inherited from the previous L256 phase scan.

Per-run analysis: cw_four_init_scan.py; aggregate with --summary --manifest. Each initialization remains separately labeled; noise/stripe are never folded into chi0.
