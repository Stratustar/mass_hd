#!/usr/bin/env python3
"""Stage A check of the s = 1/3 RESCALED campaign: five numbers, one acoustic test, a verdict.

    cw_s3_check.py <results_dir> [--sim SIM_DIR] [--old calib_step.json] [--out calib_s3.json]
                   [--figs DIR] [--s 0.3333]

<results_dir> holds one <rung>/part.json per open-loop run (written in-job by cw_part.py);
SIM_DIR holds the matching <rung>/video_*.u8 streams and defaults to <results_dir> with
/results/cases/ replaced by /cases/, i.e. the scratch layout submit_case.sh produces.

WHAT IT DECIDES. The rescaling divides every stress by s and multiplies every time by 1/s,
and it is a similarity transform of the incompressible Stokes problem. It is NOT one of the
lattice Boltzmann, which keeps its sound speed and its viscosity: Ma falls (the point of the
exercise) and Re falls with it (a side effect). So the ladder is re-measured and compared
with the OLD ladder put through the same factors:

    sigma_P / zeta_eff       invariant              |ratio - 1| < 5%
    L_P                      invariant              |ratio - 1| < 5%
    f(pmem = 0.5 sigma_P(1)) invariant              |ratio - 1| < 5%
    tau_c (steps)            x 1/s                  |ratio - 1| < 5%
    u_rms                    x s                    reported, not gated

The old rungs are a = 0.3, 0.45, 0.6, 0.8, 1; the new ones add 0.65 and 0.5, which are
compared against a log-log interpolation of the old ladder and flagged as such.

THE ACOUSTIC TEST reads the (1,0) and (0,1) Fourier coefficients of P off the video stream,
frame by frame. The stream is the 2x2 block mean of P quantised to uint8 on +/- 6 sigma_P:
block-averaging attenuates the box mode by cos(pi/L) ~ 1 and the quantisation noise reaches
the coefficient at 0.05 sigma_P/sqrt(12)/sqrt(N_pix) ~ 5e-5 sigma_P, both negligible against
the 0.05 pmem = 0.025 sigma_P criterion. Two readings are reported, because they answer
different questions:

  * the raw amplitude 2|P_hat(1,0)|(t) -- the criterion as posed: it must stay below
    0.05 pmem for the whole run. Turbulence alone puts ~ sigma_P L_P sqrt(2 pi)/L into this
    coefficient (a homogeneous field with correlation length L_P has a flat spectrum at
    k -> 0), which at L = 500 is the same order as the threshold, so a raw exceedance does
    not by itself say "sound".
  * the LINE at Omega_1 = 2 pi c_s / L in the temporal spectrum of P_hat(1,0)(t). A box
    mode is a coherent oscillation at exactly that frequency; turbulence is broadband on
    1/tau_c, and Omega_1 tau_c ~ 4.9 puts the line well outside the broad peak. The
    line's rms amplitude and its ratio to the neighbouring background are what decide.

Everything measured lands in calib_s3.json, in the same shape as calib_step.json so the B
stage generator can read it: tau_c, sigma_P, u_rms, pmem, f_top, f_floor, rungs, f_table.
"""
import argparse
import glob
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cw_s3_gen as gen                                   # noqa: E402
from cw_step_pick import f_hard                            # noqa: E402

CS = 1.0 / math.sqrt(3.0)
TOL = 0.05
PMEM_COEFFS = (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
LINE_BAND = (0.85, 1.15)          # the acoustic line, in units of Omega_1
BACK_BANDS = ((0.5, 0.85), (1.15, 1.5))


# ------------------------------------------------------------------ the parts

def load_parts(d):
    out = []
    for p in sorted(glob.glob(os.path.join(d, "*", "part.json"))):
        with open(p) as fh:
            q = json.load(fh)
        q["_dir"] = os.path.dirname(p)
        out.append(q)
    if not out:
        raise RuntimeError(f"no part.json under {d}")
    return out


def lag_tau(p):
    L = p.get("lagrangian")
    return L["tau_c_lag"] if L and np.isfinite(L.get("tau_c_lag", np.nan)) else np.nan


# --------------------------------------------------------------- the acoustics

def mode_coefficients(frames_u8, lo, hi, stride, L):
    """(c10, c01) per frame: c_k = (1/N) sum P e^{-i k.x}, with x at the block centres.

    A wave P = A cos(2 pi x/L + phase) gives |c10| = A/2, so the physical amplitude is
    2|c|. Works on the raw uint8 memmap in chunks; the mean over the transverse axis is
    taken first so the complex sum is over nx numbers per frame, not nx*ny."""
    n, nx, ny = frames_u8.shape
    scale = (hi - lo) / 255.0
    xc = (np.arange(nx) * stride + 0.5 * (stride - 1)) / L
    yc = (np.arange(ny) * stride + 0.5 * (stride - 1)) / L
    wx = np.exp(-2j * np.pi * xc)
    wy = np.exp(-2j * np.pi * yc)
    c10 = np.empty(n, dtype=np.complex128)
    c01 = np.empty(n, dtype=np.complex128)
    CH = 64
    for a in range(0, n, CH):
        b = min(a + CH, n)
        blk = frames_u8[a:b].astype(np.float64) * scale + lo
        c10[a:b] = blk.mean(axis=2) @ wx / nx
        c01[a:b] = blk.mean(axis=1) @ wy / ny
    return c10, c01


def line_power(c, dt, L, i0):
    """rms amplitude in the acoustic line and the line/background ratio, record window only.

    c is the complex coefficient sampled every dt steps; the spectrum is two-sided because a
    travelling wave has one sign of frequency and a standing one both. Parseval: with
    C_k = (1/n) sum c e^{-2 pi i k t/n}, sum_k |C_k|^2 = <|c|^2>, so the rms of the part of
    c(t) living in a band is sqrt(sum_band |C_k|^2), and the wave amplitude is twice that."""
    x = c[i0:] - c[i0:].mean()
    n = len(x)
    if n < 16:
        return None
    C = np.fft.fft(x) / n
    f = np.fft.fftfreq(n, d=dt)                 # cycles per step, signed
    f1 = CS / L                                  # the box mode
    af = np.abs(f) / f1
    line = (af >= LINE_BAND[0]) & (af <= LINE_BAND[1])
    back = np.zeros_like(line)
    for lo, hi in BACK_BANDS:
        back |= (af >= lo) & (af <= hi)
    P = np.abs(C) ** 2
    nyq = 0.5 / dt
    out = {
        "f1_per_step": f1, "nyquist_per_step": nyq,
        "resolved": bool(nyq > LINE_BAND[1] * f1),
        "n_line_bins": int(line.sum()), "n_back_bins": int(back.sum()),
        "amp_line_rms": 2.0 * float(np.sqrt(P[line].sum())) if line.any() else float("nan"),
        # background at the SAME bandwidth as the line, so the two are comparable
        "amp_back_rms": (2.0 * float(np.sqrt(np.median(P[back]) * line.sum()))
                         if back.any() and line.any() else float("nan")),
        "amp_total_rms": 2.0 * float(np.sqrt(P.sum())),
    }
    out["line_over_back"] = (out["amp_line_rms"] / out["amp_back_rms"]
                             if out["amp_back_rms"] > 0 else float("nan"))
    # PEAK amplitude, the quantity the criterion is stated in. A travelling wave keeps |c|
    # constant, so its peak equals the rms; a standing wave's coefficient oscillates and its
    # peak is sqrt(2) x rms. The box mode is some mixture, so the standing-wave value is the
    # conservative (upper) reading and is the one gated. Verified on a synthetic standing
    # wave: A = 1.20e-3 in, 1.20e-3 out.
    out["amp_line_peak"] = math.sqrt(2.0) * out["amp_line_rms"]
    return out


def acoustic_pass(sim_root, par, t_start, pmem, coeff=gen.ACOUSTIC_COEFF):
    """Everything the acoustic test reports for one run; None without a stream."""
    import cw_stream
    try:
        st = cw_stream.Stream(sim_root, par)
    except Exception as exc:
        return {"error": f"{type(exc).__name__}: {exc}"}
    L = int(par["LX"])
    lo, hi = st.limits("P")
    c10, c01 = mode_coefficients(st.raw("P"), lo, hi, st.stride, L)
    A10, A01 = 2 * np.abs(c10), 2 * np.abs(c01)
    thr = coeff * pmem
    rec = st.steps >= t_start
    i0 = int(np.argmax(rec)) if rec.any() else len(st.steps)
    dt = float(np.median(np.diff(st.steps))) if st.n > 1 else float("nan")
    out = {
        "threshold": thr, "threshold_coeff": coeff, "pmem_used": pmem,
        "nvideo": dt, "frames": int(st.n), "record_from_frame": i0,
        "samples_per_period": (L / CS) / dt if dt else float("nan"),
        "max_full": float(max(A10.max(), A01.max())),
        "max_record": float(max(A10[i0:].max(), A01[i0:].max())) if i0 < st.n else float("nan"),
        "rms_record_10": float(np.sqrt(np.mean(A10[i0:] ** 2))) if i0 < st.n else float("nan"),
        "rms_record_01": float(np.sqrt(np.mean(A01[i0:] ** 2))) if i0 < st.n else float("nan"),
        "frac_above_full": float(np.mean(np.maximum(A10, A01) > thr)),
        "line_10": line_power(c10, dt, L, i0),
        "line_01": line_power(c01, dt, L, i0),
        "_t": st.steps, "_A10": A10, "_A01": A01, "_c10": c10, "_c01": c01, "_i0": i0, "_dt": dt,
    }
    out["pass_raw"] = bool(out["max_full"] < thr)
    lines = [v for v in (out["line_10"], out["line_01"]) if v]
    out["amp_line_max"] = max(v["amp_line_peak"] for v in lines) if lines else float("nan")
    out["line_over_back_max"] = max(v["line_over_back"] for v in lines) if lines else float("nan")
    out["pass_line"] = bool(lines and out["amp_line_max"] < thr)
    return out


# ------------------------------------------------------------------- figures

def fig_ladder(rows, s, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    a = np.array([r["a"] for r in rows])
    panels = [("sigma_P / zeta_eff", "sigma_P_over_zeta", "sigma_P_over_zeta_old", 1.0),
              ("tau_c  [steps]", "tau_c", "tau_c_old", 1.0 / s),
              ("L_P", "L_P", "L_P_old", 1.0),
              ("u_rms", "u_rms", "u_rms_old", s),
              ("f(pmem = 0.5 sigma_P(1))", "f", "f_old", 1.0)]
    fig, axs = plt.subplots(1, 5, figsize=(17, 3.6))
    for ax, (lab, kn, ko, fac) in zip(axs, panels):
        new = np.array([r[kn] for r in rows])
        old = np.array([r[ko] for r in rows]) * fac
        ax.plot(a, old, "o--", color="0.6", label=f"old x {fac:g}")
        ax.plot(a, new, "s-", color="C3", label="L = 500, s = 1/3")
        ax.set_xlabel("a = zeta_eff / zeta")
        ax.set_title(lab, fontsize=10)
        ax.grid(alpha=0.3)
    axs[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def fig_acoustic(rows, tau_c, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = [r for r in rows if r.get("_ac") and "_A10" in r["_ac"]]
    if not rows:
        return
    fig, axs = plt.subplots(len(rows), 2, figsize=(13, 2.4 * len(rows)), squeeze=False)
    for (axt, axf), r in zip(axs, rows):
        ac = r["_ac"]
        t = ac["_t"] / tau_c
        axt.plot(t, ac["_A10"], lw=0.6, label="2|P(1,0)|")
        axt.plot(t, ac["_A01"], lw=0.6, label="2|P(0,1)|", alpha=0.8)
        axt.axhline(ac["threshold"], color="k", ls="--", lw=0.8, label="0.05 pmem")
        axt.axvline(ac["_t"][ac["_i0"]] / tau_c if ac["_i0"] < len(t) else t[-1],
                    color="0.5", lw=0.6)
        axt.set_ylabel(f"a = {r['a']:g}")
        axt.set_xlabel("t / tau_c")
        axt.grid(alpha=0.3)
        if r is rows[0]:
            axt.legend(fontsize=7, ncol=3)
        x = ac["_c10"][ac["_i0"]:]
        n = len(x)
        if n >= 16:
            C = np.abs(np.fft.fft(x - x.mean()) / n)
            f = np.fft.fftfreq(n, d=ac["_dt"])
            f1 = CS / int(r["L"])
            o = np.argsort(f)
            axf.semilogy(f[o] / f1, 2 * C[o], lw=0.6)
            for sgn in (-1, 1):
                axf.axvline(sgn, color="C3", lw=0.8, ls="--")
            axf.set_xlabel("omega / Omega_1   (Omega_1 = 2 pi c_s / L)")
            axf.set_xlim(-2.5, 2.5)
            axf.grid(alpha=0.3)
            lp = ac["line_10"]
            if lp:
                axf.set_title(f"line peak {lp['amp_line_peak']:.2e}, line/back {lp['line_over_back']:.1f}",
                              fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results")
    ap.add_argument("--sim", default=None)
    ap.add_argument("--old", default=None,
                    help="calib_step.json of the L = 800 ladder; default = the table "
                         "carried in cw_s3_gen.OLD_RUNGS")
    ap.add_argument("--out", default=None, help="default <results>/calib_s3.json")
    ap.add_argument("--figs", default=None, help="default <results>")
    ap.add_argument("--s", type=float, default=gen.S, help="stress factor (self-test: 1)")
    ap.add_argument("--record-frac", type=float, default=0.2)
    ap.add_argument("--no-acoustic", action="store_true")
    a = ap.parse_args()

    s = a.s
    res_dir = a.results.rstrip("/")
    sim_dir = a.sim or res_dir.replace("/results/cases/", "/cases/")
    out_path = a.out or os.path.join(res_dir, "calib_s3.json")
    fig_dir = a.figs or res_dir
    os.makedirs(fig_dir, exist_ok=True)

    if a.old:
        with open(a.old) as fh:
            old = json.load(fh)
        rungs = {float(k): v for k, v in old["rungs"].items()}
        f05 = next(t for t in old["f_table"] if abs(t["coeff"] - 0.5) < 1e-9)
        ks = sorted(rungs)
        gen.OLD_RUNGS = {k: (rungs[k]["tau_c"], rungs[k]["sigma_P"], rungs[k]["u_rms"],
                             rungs[k]["L_P"], f05["f"][ks.index(k)]) for k in ks}
    old_as = sorted(gen.OLD_RUNGS)

    parts = [p for p in load_parts(res_dir) if p["params"]["open_loop"]]
    parts.sort(key=lambda p: -p["params"]["zeta_open"])
    zeta = parts[0]["params"]["zeta"]
    top = parts[0]
    a_top = top["params"]["zeta_open"] / zeta
    if abs(a_top - 1.0) > 1e-6:
        print(f"WARNING: top rung is a = {a_top:.3f}, not 1; pmem is taken from it anyway")

    tau_c = lag_tau(top)
    if not np.isfinite(tau_c):
        tau_c = top["flow"]["tau_c"]
        print("WARNING: no tracer tau_c at the top rung; using the Eulerian one")
    sigma_P = top["flow"]["sigma_P"]
    u_rms = top["flow"]["u_rms"]
    pmem = 0.5 * sigma_P
    L = int(top["params"]["L"])
    omega1_tau_c = 2 * math.pi * CS * tau_c / L

    print(f"s = {s:.4f}, L = {L}, zeta = {zeta:g}: tau_c(1) = {tau_c:.1f} steps "
          f"(old {gen.OLD_TAU_C:.1f} x {1/s:.3g} = {gen.OLD_TAU_C/s:.1f}), "
          f"sigma_P(1) = {sigma_P:.5f}, pmem = 0.5 sigma_P = {pmem:.6f}")
    print(f"Omega_1 tau_c = 2 pi c_s tau_c / L = {omega1_tau_c:.2f}     "
          f"Ma(1) = {top['flow']['Ma']:.4f}")

    # ---------------------------------------------------------- the ladder
    rows = []
    for p in parts:
        ap_ = p["params"]["zeta_open"] / zeta
        tc_old, sp_old, ur_old, lp_old, f_old = gen.old_at(ap_)
        fl = p["flow"]
        tl = lag_tau(p)
        tl = tl if np.isfinite(tl) else fl["tau_c"]
        r = {
            "a": ap_, "case": p["case"], "L": L, "seed": p["params"].get("seed"),
            "zeta_eff": p["params"]["zeta_open"],
            "sigma_P": fl["sigma_P"], "sigma_P_over_zeta": fl["sigma_P"] / p["params"]["zeta_open"],
            "tau_c": tl, "tau_c_eulerian": fl["tau_c"], "L_P": fl["L_P"], "u_rms": fl["u_rms"],
            "Ma": fl["Ma"], "N_def": fl["N_def"], "melted_frac": fl["melted_frac"],
            "std_m": fl["std_m"], "f": f_hard(p, pmem),
            "f_lag": (p.get("lagrangian") or {}).get("f_lag"),
            "old_interpolated": bool(min(abs(ap_ - k) for k in old_as) > 1e-6),
            "tau_c_old": tc_old, "sigma_P_old": sp_old,
            "sigma_P_over_zeta_old": sp_old / (gen.OLD["zeta"][0] * ap_),
            "u_rms_old": ur_old, "L_P_old": lp_old, "f_old": f_old,
            "warnings": p.get("warnings", []),
        }
        r["ratio"] = {
            "sigma_P_over_zeta": r["sigma_P_over_zeta"] / r["sigma_P_over_zeta_old"],
            "tau_c": r["tau_c"] / (tc_old / s),
            "L_P": r["L_P"] / lp_old,
            "u_rms": r["u_rms"] / (ur_old * s),
            "f": r["f"] / f_old if f_old > 0 else float("nan"),
        }
        r["pass"] = {k: bool(abs(v - 1) < TOL) for k, v in r["ratio"].items()
                     if k in ("sigma_P_over_zeta", "tau_c", "L_P", "f")}
        rows.append(r)

    print(f"\n{'a':>5} {'sP/z':>7} {'ratio':>6} {'tau_c':>7} {'ratio':>6} {'L_P':>6} {'ratio':>6} "
          f"{'u_rms':>9} {'ratio':>6} {'f':>7} {'ratio':>6} {'Ma':>6} {'N_def':>6} {'melt':>6} {'std_m':>7}")
    for r in rows:
        q = r["ratio"]
        flag = "*" if r["old_interpolated"] else " "
        print(f"{r['a']:5.2f}{flag}{r['sigma_P_over_zeta']:7.4f} {q['sigma_P_over_zeta']:6.3f} "
              f"{r['tau_c']:7.0f} {q['tau_c']:6.3f} {r['L_P']:6.2f} {q['L_P']:6.3f} "
              f"{r['u_rms']:9.3e} {q['u_rms']:6.3f} {r['f']:7.4f} {q['f']:6.3f} "
              f"{r['Ma']:6.4f} {r['N_def']:6.0f} {r['melted_frac']:6.1%} {r['std_m']:7.4f}")
    print("  ratio = new / (old x factor); * = old value interpolated between rungs")
    fails = [(r["a"], k) for r in rows for k, ok in r["pass"].items() if not ok]
    ladder_pass = not fails
    print(f"\nLADDER {'PASS' if ladder_pass else 'FAIL'} at {TOL:.0%}"
          + ("" if ladder_pass else ": " + ", ".join(f"a={a_:g} {k}" for a_, k in fails)))

    # ------------------------------------------------------- the f table (for B)
    f_table = []
    for c in PMEM_COEFFS:
        fs = [f_hard(p, c * sigma_P) for p in parts]
        ratios = [r["a"] for r in rows]
        ft = float(np.interp(1.0, ratios[::-1], fs[::-1]))
        ff = float(np.interp(gen.OLD["zeta0-frac"][0], ratios[::-1], fs[::-1]))
        f_table.append({"coeff": c, "pmem": c * sigma_P, "f": fs, "a": ratios,
                        "f_top": ft, "f_floor": ff, "contrast": ft - ff})
    f_top = next(t for t in f_table if abs(t["coeff"] - 0.5) < 1e-9)["f_top"]
    f_floor = next(t for t in f_table if abs(t["coeff"] - 0.5) < 1e-9)["f_floor"]
    print(f"f_top = f(1) = {f_top:.4f} (old {gen.OLD_F_TOP:.4f}), f_floor = f(0.3) = {f_floor:.4f} "
          f"(old {gen.OLD_RUNGS[0.3][4]:.4f}), duty cycle f(a): "
          + "  ".join(f"{r['a']:g}->{r['f']:.3f}" for r in rows))

    # ------------------------------------------------------------ acoustics
    ac_pass = None
    if not a.no_acoustic:
        import cw_common as cw
        print(f"\nACOUSTIC  threshold 0.05 pmem = {gen.ACOUSTIC_COEFF * pmem:.3e}"
              f"   (turbulent floor estimate sigma_P L_P sqrt(2pi)/L x 2 at a = 1: "
              f"{2 * sigma_P * top['flow']['L_P'] * math.sqrt(2 * math.pi) / L:.1e})")
        print(f"{'a':>5} {'frames':>6} {'smp/per':>7} {'max full':>9} {'max rec':>9} {'rms10':>9} "
              f"{'rms01':>9} {'line10pk':>9} {'l/b':>5} {'line01pk':>9} {'l/b':>5} raw line")
        for r in rows:
            root = os.path.join(sim_dir, r["case"])
            try:
                par = cw.read_params(root)
            except Exception as exc:
                r["acoustic"] = {"error": f"no parameters at {root}: {exc}"}
                print(f"{r['a']:5.2f}  no stream: {exc}")
                continue
            t_start = a.record_frac * int(par["nsteps"])
            ac = acoustic_pass(root, par, t_start, pmem)
            if "error" in ac:
                r["acoustic"] = ac
                print(f"{r['a']:5.2f}  {ac['error']}")
                continue
            r["_ac"] = ac
            r["acoustic"] = {k: v for k, v in ac.items() if not k.startswith("_")}
            l10, l01 = ac["line_10"] or {}, ac["line_01"] or {}
            print(f"{r['a']:5.2f} {ac['frames']:6d} {ac['samples_per_period']:7.1f} "
                  f"{ac['max_full']:9.2e} {ac['max_record']:9.2e} {ac['rms_record_10']:9.2e} "
                  f"{ac['rms_record_01']:9.2e} {l10.get('amp_line_peak', np.nan):9.2e} "
                  f"{l10.get('line_over_back', np.nan):5.1f} {l01.get('amp_line_peak', np.nan):9.2e} "
                  f"{l01.get('line_over_back', np.nan):5.1f} "
                  f"{'ok ' if ac['pass_raw'] else 'X  '} {'ok' if ac['pass_line'] else 'X'}")
        acs = [r["acoustic"] for r in rows if "pass_raw" in r.get("acoustic", {})]
        if acs:
            ac_pass = {"raw": all(x["pass_raw"] for x in acs),
                       "line": all(x["pass_line"] for x in acs)}
            print(f"ACOUSTIC raw-maximum criterion: {'PASS' if ac_pass['raw'] else 'FAIL'};  "
                  f"line at Omega_1: {'PASS' if ac_pass['line'] else 'FAIL'}")
            print("  (a raw failure with a passing line is turbulence in the (1,0) mode, not sound)")
        fig_acoustic(rows, tau_c, os.path.join(fig_dir, "s3_acoustic.png"))

    fig_ladder(rows, s, os.path.join(fig_dir, "s3_ladder.png"))

    # ---------------------------------------------------------------- write
    calib = {
        "s": s, "L": L,
        "tau_c": tau_c, "sigma_P": sigma_P, "u_rms": u_rms,
        "zeta": zeta, "zeta0_frac": gen.OLD["zeta0-frac"][0],
        "pmem": pmem, "pmem_coeff": 0.5,
        "mc": gen.OLD["mc"][0],
        "f_top": f_top, "f_floor": f_floor,
        "L_P": top["flow"]["L_P"],
        "omega1_tau_c": omega1_tau_c,
        "Ma_top": top["flow"]["Ma"],
        "rungs": {f"{r['a']:g}": {k: r[k] for k in ("tau_c", "sigma_P", "u_rms", "N_def", "L_P",
                                                   "melted_frac", "f", "std_m", "Ma", "seed")}
                  for r in rows},
        "f_table": f_table,
        "ladder_check": {"tolerance": TOL, "pass": ladder_pass,
                         "rows": [{k: r[k] for k in ("a", "ratio", "pass", "old_interpolated")}
                                  for r in rows]},
        "acoustic_check": {"threshold": gen.ACOUSTIC_COEFF * pmem, "pass": ac_pass,
                           "rows": {f"{r['a']:g}": r.get("acoustic") for r in rows}},
        "old": {"tau_c": gen.OLD_TAU_C, "sigma_P": gen.OLD_SIGMA_P, "u_rms": gen.OLD_U_RMS,
                "f_top": gen.OLD_F_TOP, "rungs": {f"{k:g}": v for k, v in gen.OLD_RUNGS.items()}},
        "note": "s = 1/3 rescaled step campaign at L = 500; every step count of stage B "
                "converts with tau_c above; pmem = 0.5 sigma_P(1) re-measured here",
    }
    with open(out_path, "w") as fh:
        json.dump(calib, fh, indent=1, default=float)
    print(f"\nwrote {out_path}, {fig_dir}/s3_ladder.png"
          + ("" if a.no_acoustic else f", {fig_dir}/s3_acoustic.png"))


if __name__ == "__main__":
    main()
