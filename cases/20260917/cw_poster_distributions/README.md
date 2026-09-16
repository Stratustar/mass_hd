# Matched pressure and memory distributions

96 runs: 4 fixed activities (.3,.5,.7,1) x8 memory times (.3,1,3,6,8,10,20,30) x3 seeds.

Core physics exactly inherited from mc=.2287, L256 figure. pc=.016838, tau_chi=202.3, r=.3, tau_c=674.3290333006435 simulation steps. Activity is prescribed (open loop).

Memory compared with Beta and physical pressure with Gaussian having the same first and second moments. No pressure recentering, clipped tails, MLE or iid-site inference. The same simulations supply both fields. All seed outcomes remain visible.

Warmup max(200tau_c,20tau_m); measurement >=600tau_c; frames every~5tau_c. Seven native fields retained by frame-light; no video. Analysis: cw_poster_distribution.py.
