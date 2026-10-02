# tau_chi scan at mc=0.2287, tau_m=0.3 tau_c

14 runs: tau_chi/tau_c = 0.3, 0.6, 1.2, 2.5, 5, 10, 20; uniform chi=0 (m0=.3301) and chi=1 (m0=.0559); one seed (190901/200901).

Copied from the tm0p3 rep1 uniform cards of cases/20260915/cw_mc02287_four_init_fullres; only tau-chi and nsteps change. L256, mc=.2287, pmem=.016838, r=.3, tau_m=202.2987 steps retained.

Preparation 100 tau_c (chi reaction off), observation 1000 tau_c, 741763 steps. Native 256x256 uint8 u/P/m/chi every 337 steps (~0.5 tau_c), 2202 frames/run.

Per-run analysis and summary: plot/python/confluent_wet/cw_tchi_scan.py.
