# Finite-size scan of the 50% pulse response

L128 (12 seeds) and L512 (3 seeds) at tau_m/tau_c = 6, 6.5, 7, 7.25, 7.5, 7.7, 7.9, 8.1, 8.3, 8.5, 9, 10, both uniform initializations; one pulse and one zero-duration control per seed. L256 is the completed cw_poster_pulse50 campaign (listed in the manifest, not re-run).

All physics and protocol inherited from cw_poster_pulse50: mc=.2287, pc=.016838, pulse pc=.008419 for 3 tau_c, r=.3, tau_chi=202.3, 2000 tau_c feedback wait, 1000 tau_c recovery. Only box size and seeds change.

Metrics (plot/python/confluent_wet/cw_pulse_fss.py): half-recovery time T_half of <chi> after the pulse, and the stationary fluctuation sigma^2 = Var_t <chi> of the controls (S = L^2 sigma^2), plus autocorrelation time, Binder cumulant and escape counts.
