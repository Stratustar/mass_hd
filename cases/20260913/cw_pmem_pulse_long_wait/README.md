# Fixed-age pmem pulses

Complete 160 previously held pulses: 128 near-critical all0/all1 pulses plus 32 all1 pulses at tau_m/tau_c=12/20. Four seeds and all four pulse durations retained.

Near-critical feedback waiting is 2000 tau_c after the original frozen preparation; recovery is at least 500 tau_c after the longest pulse. This is a prescribed age, not a guarantee of stationarity. No drift or velocity gate controls submission.

Run 192 new simulations (32 controls +160 pulses), reuse eight existing 12/20 all1 controls with identical inputs. simulation_cases.txt excludes the reused controls. prepare_cw_pulse_reuse.py validates and reduces those raw archives without copying large videos.

Summary uses cw_pmem_pulse_analysis.py --allow-drift. All paired curves remain available. A decaying paired response can be fitted despite baseline drift, retaining that diagnostic. A persistent late difference receives no forced zero-offset exponential fit.
