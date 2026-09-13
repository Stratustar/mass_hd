# Temporary pmem=0 pulses at mc=.2287

Eight memory times; binary patches at 1/2.5, all0/all1 at other times; four independent seeds. 224 pulse runs plus 56 controls. Pulse durations .3/1/3/10 reference tau_c.

Baseline pmem=.016838, pulse pmem=0. No density reset or direct flow forcing. Controls and all pulse durations repeat the same seed and history.

Feedback waiting: 100 tau_c away from the transition; 500 tau_c at 7.7/7.9/8.1/8.3. Preceded by the original max(100 tau_c,10 tau_m) reaction-frozen preparation. Post-pulse observation: at least 200/500 tau_c respectively. Short near-critical preparation is a first control check, not a stationarity guarantee.

response.csv samples full-grid moments and source occupancies every 13 steps and at events. Video stores all 256x256 sites from 20 tau_c before onset; dense through onset+110 tau_c, then sparse. Actual timestamps must determine movie speed. Late exponential relaxation is fitted only to resolved return-to-branch responses.

Submit controls with cw_pmem_pulse_analysis.py as the per-case reducer. Then run scripts_cluster/launch_cw_pmem_pulses.py as an afterany analysis job. It requires all four seeds of each tau_m/initialization group to pass the final 50/200 tau_c baseline and control-tail drift checks before submitting the 16 paired pulses. Failed groups remain explicitly held in control_gate/gate.json for locally generated longer preparation cases; individual successful seeds are not selected in isolation.
