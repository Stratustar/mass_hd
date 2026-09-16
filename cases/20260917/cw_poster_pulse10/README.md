# Dense 10% compression-sensing pulse scan

300 runs:25 memory times x2 uniform initial histories x3 seeds x(control+pulse).

g=tau_m/tau_c: 1, 2, 3, 4, 5, 6, 6.5, 7, 7.25, 7.5, 7.7, 7.9, 8.1, 8.3, 8.5, 8.75, 9, 9.5, 10, 11, 12, 14, 16, 18, 20.

Core physics exactly inherited from the mc=.2287 L256 figure. Baseline pc=.016838; pulse pc=.0151542 for3tau_c; tau_chi=202.3, r=.3. This changes memory input only.

Every point uses the same2000tau_c feedback wait and1000tau_c post-pulse observation after max(100tau_c,10tau_m) frozen preparation. These are finite ages, not imposed stationarity. Exact global response sampled every13steps and at events; sparse native fields every~100tau_c; no large video streams.

Three seeds are averaged as signed pulse-control response before fitting. Each initialization is separate. All drift/persistent/missing-fit states retained; no forced zero/infinite relaxation times. Analysis: cw_poster_pulse.py.
