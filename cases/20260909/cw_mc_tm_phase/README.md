# Long-time empirical mc / tau_m phase map

L=256: 13 memory times × 14 thresholds × 2 chi starts × 2 seeds = 728 runs.
L=384: 5 memory times × 3 thresholds × 2 chi starts × 2 seeds = 60 controls.

See manifest.json for exact grids, fixed coefficients, source hashes, timing, seeds and diagnostic tolerances. Each run prepares its own matching flow and memory with uniform chi frozen, then releases feedback. Preparation time is excluded from observation and terminal statistics. There is no simulation code change.

Per-run analysis: plot/python/confluent_wet/cw_mc_tm_scan.py via PLOT_SCRIPT. Aggregate using --summary and --manifest. The summary preserves missing, incomplete, drifting, switching and seed-sensitive cases instead of forcing a phase label.
