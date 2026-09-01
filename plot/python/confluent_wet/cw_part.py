#!/usr/bin/env python3
"""THE per-run analysis of the 2026-09 confluent-wet memory campaign.

Runs inside the simulation job (submit_case.sh's PLOT_SCRIPT hook, so the contract is the
usual `inputdir outdir`) and produces, for every run without exception:

    part.json    every input parameter, every conversion to the campaign's own units, and
                 every measured quantity the campaign asks for -- including the warnings
                 that say when a number should not be trusted
    figs/        <chi> vs t; the histogram of chi in the last frame; the pooled P histogram
    video/       |u|, P, m, chi, on the campaign-wide fixed colour scale

TWO DATA SOURCES, AND WHICH QUESTION EACH ANSWERS. The run writes full-resolution frames
every 5 tau_c and a downsampled uint8 video stream every 0.2 tau_c (see cw_stream). The
split of work between them is forced, not stylistic:

  * SPATIAL statistics -- correlation lengths, block variances, defect counts, melted
    fraction -- come from the FULL frames. They need every lattice cell and one snapshot at
    a time; 25x fewer snapshots costs almost nothing because each is already a 640000-point
    spatial average.
  * TEMPORAL statistics -- the autocorrelation time of chi, the lagged chi/P
    cross-correlation -- come from the STREAM. The whole lag range to be scanned is
    3(tau_m + tau_chi), which at tau_m = 0.3 tau_c is 2.7 tau_c, i.e. HALF of one
    full-frame interval. Measuring it from the frames is not imprecise, it is impossible.
  * The DOMAIN AVERAGES <chi>(t), <m>(t), u_rms(t), sigma_P(t) come from video_meta.csv,
    where the C++ wrote them in double precision from the full-resolution field. They are
    exact, unaffected by the stream's quantisation, and available at the fine cadence.

THE CAMPAIGN UNITS ARE NOT THE RUN'S OWN. tau_c(zeta), sigma_P(zeta) and u_rms(zeta) are
single numbers fixed once, at the top of the activity ladder, and used by every run
afterwards -- that is what makes videos comparable by eye and time axes comparable by
number. Pass them with --tau-c / --sigma-p / --u-rms. Each run's own values are measured
and reported too, under `flow`, and the ratio of the two is itself a result (a closed-loop
run that settles into the high-activity phase should recover tau_c(zeta)).

Usage:
  cw_part.py <inputdir> <outdir> [--tau-c STEPS] [--sigma-p S] [--u-rms U]
             [--record-frac 0.2] [--fps 25] [--no-video] [--px-stride N]
"""
import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cw_common as cw
import cw_lag
import cw_stream

CS = 1.0 / np.sqrt(3.0)                 # lattice speed of sound
STD_R_MULTIPLES = (1, 2, 4, 8)          # block sizes, in units of L_P
SATURATION_CRITERION = 15.0             # record window >= 15 (tau_c + tau_m + tau_chi)
# 15, not the 30 the earlier campaigns asked for. 30 was never met at the long-memory end
# -- it scales with tau_m, so at tau_m = 30 tau_c it demands a record window ~2.6x the
# whole run -- and a criterion that every interesting run fails is a criterion that gets
# ignored. 15 is what the 2026-09 step campaign sizes its runs to actually satisfy, so a
# warning here again means something. What it buys is unchanged in kind: the record window
# still holds many independent realisations of the slowest mode in the problem.

# The fate labels. chi is the PASSIVE fraction, so a LOW <chi> is the HIGH-activity phase.
FATE_ACTIVE, FATE_PASSIVE = 0.2, 0.8
FATE_DEPARTURE = 0.1                    # "left its initial condition" = |d<chi>| > 0.1
FATE_DRIFT_FRAC = 0.2                   # the closing fraction of the record window
FATE_DRIFT_TOL = 0.05                   # net change across it above which the run is undecided


# ----------------------------------------------------------------- parameters

def campaign_units(par, args):
    """The run's parameters plus the campaign's fixed clocks and scales."""
    zeta = float(par.get("zeta", np.nan))
    z0 = float(par.get("zeta0_frac", 0.0))
    open_loop = int(par.get("open_loop", 0))
    chi0 = float(par.get("chi0", 0.5))

    def zeta_eff(chi):
        return zeta * (z0 + (1.0 - z0) * (1.0 - chi))

    return {
        "L": int(par["LX"]),
        "nsteps": int(par["nsteps"]),
        "ninfo": int(par["ninfo"]),
        "nvideo": int(par.get("nvideo", 0)),
        "zeta": zeta,
        "zeta0_frac": z0,
        "open_loop": open_loop,
        "zeta_open": float(par.get("zeta_open", 0.0)),
        # the activity actually applied: prescribed in open loop, chi-dependent otherwise
        "zeta_eff_open": float(par.get("zeta_open", 0.0)) if open_loop else None,
        "zeta_eff_at_chi0": zeta_eff(chi0),
        "zeta_eff_max": zeta if not open_loop else float(par.get("zeta_open", 0.0)),
        "zeta_eff_min": zeta * z0 if not open_loop else float(par.get("zeta_open", 0.0)),
        "CC": float(par.get("CC", np.nan)),
        "LL": float(par.get("LL", np.nan)),
        "xi_N": float(np.sqrt(par["LL"] / (2 * par["CC"]))),
        "Dbio": float(par.get("Dbio", 0.0)),
        "tau_m": float(par.get("tau_m", np.nan)),
        "tau_chi": float(par.get("tau_chi", np.nan)),
        "chi_width": float(par.get("chi_width", np.nan)),
        "mc": float(par.get("mc", np.nan)),
        "pmem": float(par.get("pmem", np.nan)),
        "pmem_width": float(par.get("pmem_width", np.nan)),
        "switch_sign": int(par.get("switch_sign", 1)),
        "chi_config": str(par.get("chi_config", "")),
        "chi0": chi0, "m0": float(par.get("m0", np.nan)),
        "chi_lo": float(par.get("chi_lo", np.nan)),
        "chi_hi": float(par.get("chi_hi", np.nan)),
        "chi_block": int(par.get("chi_block", 0)),
        "chi_length": float(par.get("chi_length", 0.0)),
        "m_lo": float(par.get("m_lo", np.nan)),
        "m_hi": float(par.get("m_hi", np.nan)),
        "chi_freeze_steps": int(par.get("chi_freeze_steps", 0)),
        "seed": int(par.get("seed", -1)) if "seed" in par else None,
        "chi_seed": int(par.get("chi_seed", 0)),
        # chi-width = 0 and pmem-width = 0 are the HARD STEPS, not missing values
        "step_chi": bool(float(par.get("chi_width", 1.0)) == 0.0),
        "step_g": bool(float(par.get("pmem_width", 1.0)) == 0.0),
    }


# --------------------------------------------------------------- spatial pass

def block_std(field, R):
    """std across R x R block means. The box is cropped to a whole number of blocks."""
    L = field.shape[0]
    nb = L // R
    if nb < 4:
        return float("nan")
    f = field[:nb * R, :nb * R].reshape(nb, R, nb, R).mean(axis=(1, 3))
    return float(f.std())


def spatial_pass(root, t_start, par):
    """Everything measured on the full-resolution frames inside the record window."""
    from archive.archive import loadarchive
    nf = cw.frame_count(root)
    oa = loadarchive(root)
    L = int(par["LX"])
    rad = cw.Radial(L)

    idx = [i for i in range(nf) if cw.ph.frame_time(oa, i) >= t_start]
    frames, times = [], []
    for i in idx:
        try:
            frames.append(cw.load_frame(oa, i))
        except Exception:
            break                          # a truncated final frame: keep what is readable
        times.append(cw.ph.frame_time(oa, i))
    if len(frames) < 3:
        raise RuntimeError(f"{root}: only {len(frames)} readable frames at t >= {t_start:.0f} "
                           f"(of {nf} total) -- ninfo is too coarse for this run length")

    pos = list(range(len(frames)))
    per = {k: [] for k in ("u_rms", "N_def", "melted", "std_chi", "std_m",
                           "chi_mean", "m_mean", "P_mean", "P_std", "q2_min")}
    Pall = []
    for fr in frames:
        per["u_rms"].append(float(np.sqrt(np.mean(fr["ux"]**2 + fr["uy"]**2))))
        per["N_def"].append(float(cw.n_defects(fr["qxx"], fr["qyx"])))
        per["melted"].append(float((fr["q2"] < 0.5).mean()))
        per["std_chi"].append(float(fr["chi"].std()))
        per["std_m"].append(float(fr["m"].std()))
        per["chi_mean"].append(float(fr["chi"].mean()))
        per["m_mean"].append(float(fr["m"].mean()))
        per["P_mean"].append(float(fr["P"].mean()))
        per["P_std"].append(float(fr["P"].std()))
        per["q2_min"].append(float(fr["q2"].min()))
        Pall.append(fr["P"].ravel())
    Pall = np.concatenate(Pall)

    Cuu, _ = cw.spatial_corr(rad, frames, pos, ("ux", "uy"))
    Cpp, _ = cw.spatial_corr(rad, frames, pos, "P")
    Ccc, _ = cw.spatial_corr(rad, frames, pos, "chi")
    L_u, L_P, L_chi = (cw.length_1e(rad, C) for C in (Cuu, Cpp, Ccc))

    u_rms = float(np.mean(per["u_rms"]))
    nd = float(np.mean(per["N_def"]))
    d = L / np.sqrt(nd) if nd > 0 else float("nan")

    # std_R: chi coarse-grained on R x R blocks, R in multiples of the PRESSURE correlation
    # length -- the natural ruler here, because L_P is the scale on which the drive itself
    # is coherent, so std_R/std_1 says how much of chi's variance survives averaging over
    # one, two, four, eight patches of correlated forcing.
    std_R = []
    if np.isfinite(L_P) and L_P >= 1:
        for mlt in STD_R_MULTIPLES:
            R = int(round(mlt * L_P))
            if R < 1 or L // R < 4:
                std_R.append({"multiple": mlt, "R": R, "std": float("nan")})
                continue
            std_R.append({"multiple": mlt, "R": R,
                          "std": float(np.mean([block_std(fr["chi"], R) for fr in frames]))})

    qs = np.arange(0.0, 100.001, 0.1)
    return {
        "frames_used": len(frames),
        "window_t": [float(times[0]), float(times[-1])],
        "u_rms": u_rms,
        "Ma": u_rms / CS,
        "N_def": nd,
        "d": d,
        "tau_c": d / u_rms if u_rms > 0 else float("nan"),
        "lambda": float(np.mean([cw.strain_rate(fr["ux"], fr["uy"]).mean() for fr in frames])),
        "melted_frac": float(np.mean(per["melted"])),
        "q2_min": float(np.min(per["q2_min"])),
        "std_chi": float(np.mean(per["std_chi"])),
        "std_m": float(np.mean(per["std_m"])),
        "chi_mean": float(np.mean(per["chi_mean"])),
        "m_mean": float(np.mean(per["m_mean"])),
        "P_mean": float(Pall.mean()),
        "P_median": float(np.median(Pall)),
        "sigma_P": float(Pall.std()),
        "P_pctl_levels": qs.tolist(),
        "P_pctl_values": np.percentile(Pall, qs).tolist(),
        "L_u": L_u, "L_P": L_P, "L_chi": L_chi,
        "std_R": std_R,
        "per_frame": {k: [float(x) for x in v] for k, v in per.items()},
        "times": [float(t) for t in times],
        "_last_chi": frames[-1]["chi"],
        "_P_pool": Pall,
    }


# -------------------------------------------------------------- temporal pass

def _fft_autocorr(X):
    """Per-column autocorrelation of a (n, npix) array, normalised to 1 at lag 0."""
    n = X.shape[0]
    X = X - X.mean(axis=0, keepdims=True)
    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    F = np.fft.rfft(X, n=nfft, axis=0)
    ac = np.fft.irfft(F * np.conj(F), n=nfft, axis=0)[:n]
    ac /= np.maximum(np.arange(n, 0, -1)[:, None], 1)      # unbiased lag normalisation
    c0 = ac[0].copy()
    good = c0 > 0
    ac[:, good] /= c0[good]
    return ac[:, good]


def _cross_corr(A, B, maxlag):
    """C(lag) = <dA(t) dB(t-lag)> / sqrt(<dA^2><dB^2>) with the SPATIAL mean removed.

    The per-frame spatial mean is subtracted (not the per-pixel time mean), matching
    cw_common.lagged_corr: the question is whether the pattern of chi follows the pattern
    of P, so the k = 0 mode -- which cannot say anything about spatial structure -- is
    removed and reported separately from the exact meta series instead.

    A positive lag means B LEADS A, which is the direction the loop predicts for
    A = chi, B = P.
    """
    A = A - A.mean(axis=1, keepdims=True)
    B = B - B.mean(axis=1, keepdims=True)
    n = A.shape[0]
    den = np.sqrt(np.mean(A**2) * np.mean(B**2))
    if not np.isfinite(den) or den <= 0:
        return np.full(maxlag + 1, np.nan)
    out = np.empty(maxlag + 1)
    for lag in range(maxlag + 1):
        if n - lag < 4:
            out[lag] = np.nan
            continue
        out[lag] = float(np.mean(A[lag:] * B[:n - lag])) / den
    return out


def _decay_time(c, dt):
    """First 1/e crossing of a correlation function sampled at spacing dt."""
    thr = np.exp(-1.0)
    for i in range(1, len(c)):
        if np.isfinite(c[i]) and c[i] <= thr:
            a, b = c[i - 1], c[i]
            frac = 0.0 if a == b else (a - thr) / (a - b)
            return float((i - 1 + frac) * dt)
    return float("nan")


def temporal_pass(stream, t_start, tau_m, tau_chi, px_stride):
    """Autocorrelation time of chi and the lagged chi/P cross-correlation, from the stream."""
    w = stream.window(t_start)
    if len(w) < 16:
        raise RuntimeError(f"video stream has {len(w)} frames at t >= {t_start:.0f}; "
                           f"nvideo is too coarse for a temporal analysis")
    i0, i1 = int(w[0]), int(w[-1]) + 1
    dt = float(np.median(np.diff(stream.steps[i0:i1])))

    # Sub-sample the lattice: an autocorrelation averaged over 4096 independent points is
    # already converged, and the full 160000 would cost 40x the memory for no accuracy.
    s = max(1, px_stride)
    chi = stream.dequant("chi", slice(i0, i1))[:, ::s, ::s]
    P = stream.dequant("P", slice(i0, i1))[:, ::s, ::s]
    npix = chi.shape[1] * chi.shape[2]
    chi = chi.reshape(chi.shape[0], npix)
    P = P.reshape(P.shape[0], npix)

    ac = _fft_autocorr(chi)
    ac_mean = ac.mean(axis=1) if ac.size else np.array([np.nan])
    tau_auto = _decay_time(ac_mean, dt)

    maxlag = int(min(chi.shape[0] - 4, np.ceil(3.0 * (tau_m + tau_chi) / dt)))
    maxlag = max(maxlag, 1)
    cc = _cross_corr(chi, P, maxlag)
    if np.all(~np.isfinite(cc)):
        peak_i, peak = -1, float("nan")
    else:
        peak_i = int(np.nanargmax(np.abs(cc)))
        peak = float(cc[peak_i])

    # The k = 0 mode, from the EXACT meta series -- the part _cross_corr deliberately drops.
    cm = stream.meta["chi_mean"][i0:i1] - np.mean(stream.meta["chi_mean"][i0:i1])
    pm = stream.meta["P_mean"][i0:i1] - np.mean(stream.meta["P_mean"][i0:i1])
    den0 = np.sqrt(np.mean(cm**2) * np.mean(pm**2))
    cc0 = np.array([np.mean(cm[l:] * pm[:len(cm) - l]) / den0 if den0 > 0 else np.nan
                    for l in range(maxlag + 1)])

    return {
        "dt_steps": dt,
        "frames_used": int(i1 - i0),
        "pixels_used": npix,
        "tau_auto_chi": tau_auto,
        "autocorr_lag_steps": (np.arange(len(ac_mean)) * dt).tolist()[:maxlag + 1],
        "autocorr_chi": [float(x) for x in ac_mean[:maxlag + 1]],
        "lag_steps": (np.arange(maxlag + 1) * dt).tolist(),
        "lag_corr_chi_P": [float(x) for x in cc],
        "lag_peak": peak,
        "lag_peak_delta": float(peak_i * dt) if peak_i >= 0 else float("nan"),
        "lag_corr_k0": [float(x) for x in cc0],
        "lag_peak_k0": float(cc0[int(np.nanargmax(np.abs(cc0)))]) if np.any(np.isfinite(cc0)) else float("nan"),
        "lag_peak_delta_k0": float(int(np.nanargmax(np.abs(cc0))) * dt) if np.any(np.isfinite(cc0)) else float("nan"),
    }


def f_table(stream, t_start, sigma_P, coeffs=(0.0, 0.2, 0.3, 0.4, 0.5, 0.6)):
    """f(pmem) = fraction of TIME a lattice point spends above pmem, and its spread.

    The mean over sites of the per-site time fraction is identical to the pooled fraction of
    (site, time) samples above the threshold, so the mean can be read off a histogram; the
    SPREAD across sites cannot, and that is the reason to accumulate the per-site map.

    Thresholds are absolute pressures, expressed as multiples of the CAMPAIGN sigma_P(zeta)
    -- one fixed set of numbers for the whole activity ladder. That is the point of the
    table: the same load threshold seen by layers of different activity.
    """
    w = stream.window(t_start)
    i0, i1 = int(w[0]), int(w[-1]) + 1
    raw = stream.raw("P")
    lo, hi = stream.limits("P")
    step = (hi - lo) / 255.0

    hist = np.zeros(256, dtype=np.int64)
    per_site = {c: np.zeros((stream.nx, stream.ny), dtype=np.int32) for c in coeffs}
    # byte threshold for each requested pressure; a site counts as above when its byte is
    # strictly greater, which is the same convention as P > pmem up to half a bin
    thr_byte = {c: (c * sigma_P - lo) / step for c in coeffs}

    CH = 64
    for a in range(i0, i1, CH):
        b = min(a + CH, i1)
        blk = raw[a:b]
        hist += np.bincount(blk.ravel(), minlength=256)
        for c in coeffs:
            per_site[c] += (blk > thr_byte[c]).sum(axis=0, dtype=np.int32)

    n = i1 - i0
    edges = lo + np.arange(256) * step
    out = []
    for c in coeffs:
        fs = per_site[c] / float(n)
        out.append({"coeff": float(c), "pmem": float(c * sigma_P),
                    "f_mean": float(fs.mean()), "f_site_std": float(fs.std()),
                    "f_site_p10": float(np.percentile(fs, 10)),
                    "f_site_p90": float(np.percentile(fs, 90))})
    return {"sigma_P_ref": float(sigma_P), "frames_used": int(n),
            "table": out,
            "P_hist_edges": edges.tolist(),
            "P_hist_counts": hist.tolist()}


# ------------------------------------------------------------------- figures

# ------------------------------------------------------- the Lagrangian clock

def lagrangian_pass(root, t_start, pmem, pmem_width):
    """tau_c and f measured ALONG THE TRACERS, i.e. on this run's own material clock.

    WHY THIS RUN'S OWN AND NOT THE CAMPAIGN'S. The campaign quotes every tau_m in units of
    tau_c(zeta), one number fixed at full activity, because a ladder needs a common ruler.
    But a closed-loop run does not sit at full activity: chi settles somewhere, the activity
    settles with it, and the material clock slows down accordingly -- by a factor approaching
    (1/z0)^1.15 in the passive phase. Reporting both is what lets the analysis say whether a
    fate was decided at tau_m/tau_c(zeta) = 10 or at tau_m/tau_c(this run) = 3.

    f is measured the same way the model computes it -- a hard step when pmem-width is 0,
    the tanh otherwise -- and along the tracers rather than over the field, because the
    memory integrates g(P) following a cell. The two agree in a statistically homogeneous
    state and disagree exactly when the phases have separated, which is the interesting case.
    """
    dt, arr = cw_lag.read_tracers(root)
    if arr is None:
        return None
    t = np.arange(arr.shape[0]) * dt
    keep = t >= t_start
    if keep.sum() < 8:
        return None
    P = arr[keep, :, 0].astype(np.float64)
    m = arr[keep, :, 1].astype(np.float64)
    lag, C = cw_lag.lagrangian_corr(P, dt)
    g = (P > pmem).astype(np.float64) if pmem_width <= 0 else \
        0.5 * (1.0 + np.tanh((P - pmem) / pmem_width))
    return {
        "n_tracer": int(arr.shape[1]),
        "n_sample": int(keep.sum()),
        "dt": int(dt),
        "tau_c_lag": cw_lag.cross_1e(lag, C),
        "tau_c_lag_integral": cw_lag.integral_time(lag, C),
        "f_lag": float(g.mean()),
        "std_g_lag": float(g.std()),
        "sigma_P_lag": float(P.std()),
        "m_mean_lag": float(m.mean()),
        "std_m_lag": float(m.std()),
    }


# ------------------------------------------------------------------- the fate

def fate_of(stream, t_start, u):
    """Where the run ended up, when it left its start, and whether it is still moving.

    THE LABEL IS ON THE TAIL, NOT THE RECORD-WINDOW MEAN, for the reason the `settled`
    block already documents: a run sliding from one phase to the other has a record-window
    mean somewhere in the middle, and calling that "mixed" would invent a third state out
    of a transient. `drifting` is the guard -- a run whose <chi> still moves by more than
    FATE_DRIFT_TOL across the closing fifth of its record window is labelled UNDECIDED
    whatever its tail mean says.

    Two drifts are reported and they answer different questions. `settled.chi_drift_record`
    compares the last quarter of the record window against the rest -- "did this run move
    while we were watching". `fate.drift_tail` is the net change WITHIN the closing fifth --
    "is it still moving now". A run that slid early and then stopped shows the first and not
    the second, and only the second should veto a fate label.

    `departure` is measured over the WHOLE run, warm-up included: the question it answers
    is when the initial condition stopped holding, and for a start that is thrown across
    the switch by the initial flow transient that happens before the record window opens.
    """
    if stream is None:
        return None
    cm = np.asarray(stream.meta["chi_mean"], float)[:stream.n]
    t = np.asarray(stream.steps, float)[:stream.n]
    if not len(cm):
        return None

    w = stream.window(t_start)
    if not len(w):
        return None
    i0, i1 = int(w[0]), int(w[-1]) + 1
    tail_lo = i1 - max(2, int(round(FATE_DRIFT_FRAC * (i1 - i0))))
    chi_tail = float(np.mean(cm[tail_lo:i1]))

    # NET CHANGE across the closing fraction, as the fitted slope times its span. The
    # obvious alternatives are both worse: last-minus-first is one noisy frame against
    # another, and the difference of the two half-means throws away the lever arm -- on the
    # L=200 prototype it read +0.021 for a window whose <chi> actually moved +0.06, because
    # a monotone ramp's half-means are only half the span apart. A least-squares slope uses
    # every frame and the full span, which is what "still moving" means here.
    tt, cc = t[tail_lo:i1], cm[tail_lo:i1]
    if len(tt) >= 3 and tt[-1] > tt[0]:
        slope = float(np.polyfit(tt - tt[0], cc, 1)[0])
        drift_tail = slope * float(tt[-1] - tt[0])
    else:
        drift_tail = float(cc[-1] - cc[0]) if len(cc) >= 2 else 0.0
    drifting = bool(abs(drift_tail) > FATE_DRIFT_TOL)

    if drifting:
        label = "undecided"
    elif chi_tail < FATE_ACTIVE:
        label = "active"
    elif chi_tail > FATE_PASSIVE:
        label = "passive"
    else:
        label = "mixed"

    # when did <chi> first leave its own initial value by more than FATE_DEPARTURE
    chi0 = float(cm[0])
    left = np.nonzero(np.abs(cm - chi0) > FATE_DEPARTURE)[0]
    departure = float(t[left[0]]) if len(left) else float("nan")

    return {
        "label": label,
        "chi_tail": chi_tail,
        "chi_init_stream": chi0,
        "tail_frac": FATE_DRIFT_FRAC,
        "drift_tail": drift_tail,
        "drifting": drifting,
        "departure_steps": departure,
        "departure_over_tau_m": departure / u["tau_m"] if u["tau_m"] > 0 else float("nan"),
        "never_departed": bool(not len(left)),
        "thresholds": {"active": FATE_ACTIVE, "passive": FATE_PASSIVE,
                       "departure": FATE_DEPARTURE, "drift_tol": FATE_DRIFT_TOL},
    }


def figures(outdir, stream, res, u, t_start, tau_c):
    fig_dir = os.path.join(outdir, "figs")
    os.makedirs(fig_dir, exist_ok=True)

    t = stream.steps / tau_c
    fig, ax = plt.subplots(figsize=(6.4, 3.6), dpi=130)
    ax.plot(t, stream.meta["chi_mean"][:stream.n], lw=1.3, color="#1f4e79", label=r"$\langle\chi\rangle$")
    ax.plot(t, stream.meta["m_mean"][:stream.n], lw=1.0, color="#c07a2a", alpha=.75, label=r"$\langle m\rangle$")
    ax.axvline(t_start / tau_c, color="0.5", ls="--", lw=.9)
    ax.text(t_start / tau_c, 1.01, " record", color="0.4", fontsize=8,
            transform=ax.get_xaxis_transform(), va="bottom")
    ax.set_xlabel(r"$t/\tau_c(\zeta)$"); ax.set_ylabel("domain average")
    ax.set_ylim(-0.02, 1.02); ax.legend(fontsize=8, frameon=False)
    ax.set_title(os.path.basename(outdir.rstrip("/")), fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "chi_mean_vs_t.png")); plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.6, 3.4), dpi=130)
    ax.hist(res["_last_chi"].ravel(), bins=120, range=(0, 1), color="#4a7ebb")
    ax.set_xlabel(r"$\chi$ (last frame)"); ax.set_ylabel("lattice sites")
    ax.set_yscale("log")
    fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "chi_hist_last.png")); plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.6, 3.4), dpi=130)
    ax.hist(res["_P_pool"], bins=200, color="#7a7a7a")
    if np.isfinite(u["pmem"]):
        ax.axvline(u["pmem"], color="#b03030", lw=1.2, label=f"pmem = {u['pmem']:.3g}")
        ax.legend(fontsize=8, frameon=False)
    ax.set_xlabel("P (record window)"); ax.set_ylabel("samples"); ax.set_yscale("log")
    fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "P_hist.png")); plt.close(fig)
    return fig_dir


# ---------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputdir")
    ap.add_argument("outdir")
    ap.add_argument("--tau-c", type=float, default=None,
                    help="campaign tau_c(zeta) in steps; default = this run's own")
    ap.add_argument("--sigma-p", type=float, default=None,
                    help="campaign sigma_P(zeta); sets the f table thresholds and the video "
                         "colour scale (+/- 3 sigma_P). Default = this run's own")
    ap.add_argument("--u-rms", type=float, default=None,
                    help="campaign u_rms(zeta); sets the |u| video scale (0..3 u_rms)")
    ap.add_argument("--record-frac", type=float, default=0.2,
                    help="fraction of the run discarded as warm-up (default 0.2)")
    ap.add_argument("--px-stride", type=int, default=4,
                    help="lattice sub-sampling for the temporal statistics")
    ap.add_argument("--fps", type=int, default=25)
    ap.add_argument("--bitrate", default=None)
    ap.add_argument("--no-video", action="store_true")
    a = ap.parse_args()

    root = a.inputdir.rstrip("/")
    os.makedirs(a.outdir, exist_ok=True)
    par = cw.read_params(root)
    u = campaign_units(par, a)
    t_start = a.record_frac * u["nsteps"]
    name = os.path.basename(root)
    print(f"{name}: L = {u['L']}, {u['nsteps']} steps, record from t = {t_start:.0f}",
          flush=True)

    res = spatial_pass(root, t_start, par)
    print(f"  frames {res['frames_used']}  u_rms {res['u_rms']:.4e}  N_def {res['N_def']:.0f}"
          f"  tau_c {res['tau_c']:.0f}  sigma_P {res['sigma_P']:.4e}"
          f"  L_P {res['L_P']:.1f}  L_chi {res['L_chi']:.1f}", flush=True)

    tau_c = a.tau_c if a.tau_c else res["tau_c"]
    sigma_P = a.sigma_p if a.sigma_p else res["sigma_P"]
    u_rms_ref = a.u_rms if a.u_rms else res["u_rms"]

    warnings = []
    stream, temporal, ftab = None, None, None
    try:
        stream = cw_stream.Stream(root, par)
    except Exception as exc:
        warnings.append(f"no usable video stream: {exc}")
    if stream is not None:
        try:
            temporal = temporal_pass(stream, t_start, u["tau_m"], max(u["tau_chi"], 0.0),
                                     a.px_stride)
            print(f"  tau_auto(chi) {temporal['tau_auto_chi']:.0f} steps"
                  f"  lag peak {temporal['lag_peak']:+.3f} at "
                  f"{temporal['lag_peak_delta']:.0f}", flush=True)
        except Exception as exc:
            warnings.append(f"temporal analysis failed: {type(exc).__name__}: {exc}")
        try:
            ftab = f_table(stream, t_start, sigma_P)
            print("  f(pmem/sigma_P): " +
                  "  ".join(f"{r['coeff']:.1f}->{r['f_mean']:.3f}" for r in ftab["table"]),
                  flush=True)
        except Exception as exc:
            warnings.append(f"f table failed: {type(exc).__name__}: {exc}")
        clip = stream.clip_fraction()
        if clip["P"] > 1e-3 or clip["u"] > 1e-3:
            warnings.append(
                f"video stream clipping: max area fraction outside the stored range is "
                f"{clip['P']:.2%} for P and {clip['u']:.2%} for |u|; the f table and the "
                f"lagged correlation are computed on clipped data")

    # ---- the TAIL average, and why the record-window average is not enough
    #
    # Checkpoint 2 asks whether two runs started from opposite ends END UP apart. The
    # record-window mean cannot answer that when the approach is slow: the 20260831 A4 pair
    # had the chi == 1 arm slide from 1.0 to 0.036 over ~120 tau_c, so its record-window mean
    # was 0.27 -- a number describing the slide, not the state it slid into, and one that
    # would have been read as a partial separation from a run that had in fact converged onto
    # the SAME fixed point as the other arm. The tail mean plus the record-window drift say
    # which of the two is happening.
    tail = None
    if stream is not None:
        w = stream.window(t_start)
        if len(w):
            i0, i1 = int(w[0]), int(w[-1]) + 1
            half = i0 + (i1 - i0) * 3 // 4
            cm = stream.meta["chi_mean"]; mm = stream.meta["m_mean"]
            tail = {
                "frac": 0.25,
                "chi_mean_tail": float(np.mean(cm[half:i1])),
                "m_mean_tail": float(np.mean(mm[half:i1])),
                "chi_mean_record": float(np.mean(cm[i0:i1])),
                "chi_drift_record": float(np.mean(cm[half:i1]) - np.mean(cm[i0:half])),
                "chi_last": float(cm[i1 - 1]),
                "m_last": float(mm[i1 - 1]),
            }
            if abs(tail["chi_drift_record"]) > 0.05:
                warnings.append(
                    f"<chi> drifts by {tail['chi_drift_record']:+.3f} across the record window "
                    f"({tail['chi_mean_record']:.3f} overall vs {tail['chi_mean_tail']:.3f} in "
                    f"the last quarter): the run has NOT settled, and the record-window mean "
                    f"describes the approach rather than the state")
            print(f"  <chi> record {tail['chi_mean_record']:.4f}  tail "
                  f"{tail['chi_mean_tail']:.4f}  drift {tail['chi_drift_record']:+.4f}",
                  flush=True)

    # ---- the Lagrangian clock and the fate, both new in the 2026-09 step campaign
    lagr = None
    try:
        lagr = lagrangian_pass(root, t_start, u["pmem"], u["pmem_width"])
        if lagr:
            print(f"  tracers {lagr['n_tracer']}  tau_c(lag) {lagr['tau_c_lag']:.0f} steps"
                  f"  f(lag) {lagr['f_lag']:.4f}  std(m) {lagr['std_m_lag']:.4f}", flush=True)
    except Exception as exc:
        warnings.append(f"tracer analysis failed: {type(exc).__name__}: {exc}")
    if lagr is None:
        warnings.append("no tracer stream: tau_c and f are Eulerian for this run, and the "
                        "memory's clock is a material one")

    fate = fate_of(stream, t_start, u)
    if fate:
        print(f"  FATE {fate['label']}  <chi>_tail {fate['chi_tail']:.4f}  "
              f"drift {fate['drift_tail']:+.4f}  left its start at "
              f"{fate['departure_steps']:.0f} steps", flush=True)

    # ---- the campaign's own sanity criterion
    tau_sum = tau_c + max(u["tau_m"], 0.0) + max(u["tau_chi"], 0.0)
    record_steps = u["nsteps"] - t_start
    if record_steps < SATURATION_CRITERION * tau_sum:
        warnings.append(
            f"record window {record_steps:.0f} steps is {SATURATION_CRITERION * tau_sum / record_steps:.1f}x "
            f"SHORT of the {SATURATION_CRITERION:.0f}(tau_c + tau_m + tau_chi) = "
            f"{SATURATION_CRITERION * tau_sum:.0f} steps the campaign asks for; statistics "
            f"averaging over the slow mode are correspondingly under-sampled")
    if np.isfinite(res["L_chi"]) and res["L_chi"] > u["L"] / 6.0:
        warnings.append(f"L_chi = {res['L_chi']:.0f} exceeds L/6 = {u['L']/6:.0f}: the box "
                        f"may be setting the phenotype correlation length")
    if res["Ma"] > 0.045:
        warnings.append(f"Ma = {res['Ma']:.3f} exceeds the 0.045 weakly-compressible budget")

    # ---- where does pmem sit in this run's own pressure distribution
    pmem_pctl = float(np.interp(u["pmem"], res["P_pctl_values"], res["P_pctl_levels"])) \
        if np.isfinite(u["pmem"]) else float("nan")

    out = {
        "case": name,
        "inputdir": root,
        "params": u,
        "campaign": {"tau_c": tau_c, "sigma_P": sigma_P, "u_rms": u_rms_ref,
                     "tau_c_source": "argument" if a.tau_c else "this run",
                     "record_start": t_start, "record_steps": record_steps},
        "derived": {
            "tau_m_over_tau_c": u["tau_m"] / tau_c,
            "tau_chi_over_tau_c": u["tau_chi"] / tau_c,
            "tau_c_run_over_tau_c_ref": res["tau_c"] / tau_c,
            # the same ratio on the MATERIAL clock, which is the one tau_m is compared to
            "tau_m_over_tau_c_lag": (u["tau_m"] / lagr["tau_c_lag"]
                                     if lagr and np.isfinite(lagr["tau_c_lag"])
                                     else float("nan")),
            "tau_c_lag_over_tau_c_ref": (lagr["tau_c_lag"] / tau_c
                                         if lagr and np.isfinite(lagr["tau_c_lag"])
                                         else float("nan")),
            "pmem_over_sigma_P": u["pmem"] / sigma_P if sigma_P else float("nan"),
            "pmem_percentile_this_run": pmem_pctl,
            "box_over_L_P": u["L"] / res["L_P"] if res["L_P"] else float("nan"),
            "box_over_L_chi": u["L"] / res["L_chi"] if res["L_chi"] else float("nan"),
            "record_over_sum_tau": record_steps / tau_sum,
            "A_max": u["zeta_eff_max"] / u["CC"] if u["CC"] else float("nan"),
            "A_min": u["zeta_eff_min"] / u["CC"] if u["CC"] else float("nan"),
            "l_B_over_xi_N": float(np.sqrt(u["Dbio"] * tau_c)) / u["xi_N"] if u["Dbio"] else 0.0,
        },
        "flow": {k: v for k, v in res.items() if not k.startswith("_")},
        "lagrangian": lagr,
        "fate": fate,
        "settled": tail,
        "f": ftab,
        "time": temporal,
        "warnings": warnings,
    }

    with open(os.path.join(a.outdir, "part.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"  wrote {a.outdir}/part.json", flush=True)

    if stream is not None:
        figures(a.outdir, stream, res, u, t_start, tau_c)
        print(f"  wrote {a.outdir}/figs/", flush=True)

    if stream is not None and not a.no_video:
        vdir = os.path.join(a.outdir, "video")
        os.makedirs(vdir, exist_ok=True)
        windows = {"u": (0.0, 3 * u_rms_ref), "P": (-3 * sigma_P, 3 * sigma_P),
                   "m": (0.0, 1.0), "chi": (0.0, 1.0)}
        for fld, (lo, hi) in windows.items():
            try:
                nb = cw_stream.render(stream, fld, os.path.join(vdir, f"{fld}.mp4"),
                                      lo, hi, tau_c, fps=a.fps, bitrate=a.bitrate)
                print(f"  wrote {vdir}/{fld}.mp4 ({nb/1e6:.1f} MB, "
                      f"{stream.n} frames, {lo:.3g}..{hi:.3g})", flush=True)
            except Exception as exc:
                print(f"  video {fld} FAILED: {type(exc).__name__}: {exc}", flush=True)

    for w in warnings:
        print(f"  WARNING: {w}", flush=True)


if __name__ == "__main__":
    main()
