# Full-resolution videos: five four-initialization groups

40 runs: tau_m/tau_c = 0.3, 6.821052632, 7.752631579, 8.684210526, 18; all0, all1, stripe01 and Gaussian noise; original paired seeds 190901/190902.

Each input is copied from the original 20260911 campaign with only video-stride changed from 8 to 1. L256, mc=0.2287, pmem=0.016838, tau_chi=202.3 steps and all other dynamics are retained. No pressure pulse.

Preparation 100–180 tau_c; feedback observation 2000 tau_c; total 1,416,091–1,470,037 steps. The original t=0 mixed initializations continue to advect/diffuse during preparation.

Four native256x256 uint8 streams every337steps (~0.499756tau_c), 169400 total stored timestamps,44.407GB raw streams. Exact full-domain scalar statistics retain their existing cadence. Original sparse JSON snapshots and fixed stored value ranges are retained.

Per-run analysis: cw_four_init_scan.py. Summary: same script with --summary --manifest. New output directories keep original results intact.
