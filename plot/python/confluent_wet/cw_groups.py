#!/usr/bin/env python3
"""Cross-run summaries of the production groups: G1, G3, G4, G5.

One script, four questions, because they share every convention -- the same part.json
fields, the same tail-vs-record distinction, the same calib.json for the campaign units.
Splitting them into four scripts would be four places for those conventions to drift.

    cw_groups.py <prod-results-dir> --calib calib.json [--out summary.json] [--figs DIR]

WHICH AVERAGE. Every statement about "what state a run ended in" uses `settled.chi_mean_tail`
(the last quarter), never `flow.chi_mean` (the whole record window). The first A4 wave is why:
an arm that slid from 1.0 to 0.036 over 120 tau_c had a record-window mean of 0.27, which
reads like a partial separation and is in fact the average of the slide. `settled` also
carries the drift, so a run that has not converged is visible rather than silently averaged.

G4 IS THE ONE THAT MEASURES A RATE, not a state. Its observable is the SLOPE of <chi>(t) --
the front speed -- and a slope is well defined precisely when the run has not settled, which
is the regime a front near mc* lives in. So G4 reads the video-stream series rather than the
tail mean, and mc* is where the slope crosses zero.
"""
import argparse
import json
import os
import glob

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load(root):
    """Every part.json under <root>/<group>/<case>/, keyed by (group, case)."""
    out = {}
    for f in sorted(glob.glob(os.path.join(root, "*", "*", "part.json"))):
        grp = os.path.basename(os.path.dirname(os.path.dirname(f)))
        out.setdefault(grp, []).append(json.load(open(f)))
    return out


def tail(j):
    s = j.get("settled") or {}
    return s.get("chi_mean_tail", j["flow"]["chi_mean"])


def drift(j):
    return (j.get("settled") or {}).get("chi_drift_record", float("nan"))


def by(rows, keyfn):
    d = {}
    for r in rows:
        d.setdefault(keyfn(r), []).append(r)
    return dict(sorted(d.items()))


# --------------------------------------------------------------------------- G3

def g3(rows, cal, figs):
    """Does the two-phase separation appear, and does it survive, as tau_m grows?"""
    tc = cal["tau_c"]
    print("\nG3 -- the main figure: separation of the two starts vs tau_m")
    print(f"{'tau_m/tau_c':>11} {'chi==0 arm':>22} {'chi==1 arm':>22} {'separation':>11} "
          f"{'drift0':>8} {'drift1':>8}")
    xs, seps, lo_s, hi_s = [], [], [], []
    for tm, g in by(rows, lambda r: round(r["params"]["tau_m"] / tc, 2)).items():
        a = [r for r in g if r["params"]["chi0"] == 0.0]
        b = [r for r in g if r["params"]["chi0"] == 1.0]
        if not a or not b:
            continue
        ta, tb = np.array([tail(r) for r in a]), np.array([tail(r) for r in b])
        xs.append(tm); seps.append(abs(tb.mean() - ta.mean()))
        lo_s.append((ta.mean(), ta.std())); hi_s.append((tb.mean(), tb.std()))
        print(f"{tm:11.1f} {ta.mean():12.4f} +/- {ta.std():5.4f} "
              f"{tb.mean():12.4f} +/- {tb.std():5.4f} {abs(tb.mean()-ta.mean()):11.4f} "
              f"{np.mean([drift(r) for r in a]):+8.4f} {np.mean([drift(r) for r in b]):+8.4f}")
    if figs and xs:
        fig, ax = plt.subplots(figsize=(6.0, 4.0), dpi=140)
        ax.errorbar(xs, [v[0] for v in lo_s], yerr=[v[1] for v in lo_s], marker="o",
                    color="#b03030", label=r"start $\chi\equiv0$ (active)")
        ax.errorbar(xs, [v[0] for v in hi_s], yerr=[v[1] for v in hi_s], marker="s",
                    color="#3060b0", label=r"start $\chi\equiv1$ (passive)")
        ax.axhline(0.5 * (cal["chi_lo"] + cal["chi_hi"]), color="0.6", ls=":", lw=1)
        ax.set_xscale("log"); ax.set_xlabel(r"$\tau_m/\tau_c(\zeta)$")
        ax.set_ylabel(r"$\langle\chi\rangle$  (last quarter)")
        ax.set_ylim(-0.03, 1.03); ax.legend(fontsize=8, frameon=False)
        ax.set_title("G3: do the two starts stay apart?", fontsize=10)
        fig.tight_layout(); fig.savefig(os.path.join(figs, "g3_separation.png")); plt.close(fig)
    return {"tau_m_over_tau_c": xs, "separation": seps,
            "arm_chi0": lo_s, "arm_chi1": hi_s}


# --------------------------------------------------------------------------- G4

def front_speed(j, cal):
    """d<chi>/dt over the record window, in units of 1/tau_c -- the front speed.

    Read from the video-stream series in part.json's `settled` block: the last-quarter and
    first-three-quarter means of the record window give a clean secant slope without needing
    the raw series. A front that has run to completion has slope ~0 with <chi> at a phase
    value; a front still moving has a slope whose SIGN says which phase is winning.
    """
    s = j.get("settled") or {}
    if "chi_drift_record" not in s:
        return float("nan")
    rec_steps = j["campaign"]["record_steps"]
    # The drift is (mean of the last quarter) - (mean of the first three quarters). Those
    # two segments have their centres of mass at 7/8 and 3/8 of the record window, so the
    # secant they define spans HALF the window, not five eighths.
    return s["chi_drift_record"] / (0.5 * rec_steps / cal["tau_c"])


def g4(rows, cal, figs):
    """The front speed changes sign at mc*: a far sharper locator than an area fraction."""
    tc = cal["tau_c"]
    print("\nG4 -- front speed vs mc; the zero crossing IS mc*")
    out = {}
    for tm, g in by(rows, lambda r: round(r["params"]["tau_m"] / tc, 2)).items():
        print(f"  tau_m = {tm:g} tau_c")
        print(f"    {'mc':>8} {'d<chi>/dt':>12} {'+/-':>8} {'<chi> tail':>11} {'std_chi':>8}")
        mcs, vs, es = [], [], []
        for mc, gg in by(g, lambda r: round(r["params"]["mc"], 6)).items():
            v = np.array([front_speed(r, cal) for r in gg])
            t = np.array([tail(r) for r in gg])
            mcs.append(mc); vs.append(float(v.mean())); es.append(float(v.std()))
            print(f"    {mc:8.4f} {v.mean():+12.5f} {v.std():8.5f} {t.mean():11.4f} "
                  f"{np.mean([r['flow']['std_chi'] for r in gg]):8.4f}")
        star = None
        for i in range(len(mcs) - 1):
            if vs[i] * vs[i + 1] < 0:
                star = mcs[i] + (0 - vs[i]) * (mcs[i + 1] - mcs[i]) / (vs[i + 1] - vs[i])
                break
        print(f"    -> mc* = " + (f"{star:.5f}" if star else "not bracketed by the sign change"))
        out[tm] = {"mc": mcs, "speed": vs, "speed_sd": es, "mc_star": star}
        if figs:
            fig, ax = plt.subplots(figsize=(5.6, 3.8), dpi=140)
            ax.errorbar(mcs, vs, yerr=es, marker="o", color="#1f4e79")
            ax.axhline(0, color="0.5", lw=1)
            if star: ax.axvline(star, color="#b03030", ls="--", lw=1,
                                label=f"$m_c^*$ = {star:.4f}")
            ax.set_xlabel(r"$m_c$"); ax.set_ylabel(r"$d\langle\chi\rangle/dt$  [$1/\tau_c$]")
            ax.set_title(rf"G4 front speed, $\tau_m$ = {tm:g} $\tau_c$", fontsize=10)
            ax.legend(fontsize=8, frameon=False)
            fig.tight_layout()
            fig.savefig(os.path.join(figs, f"g4_front_tm{tm:g}.png")); plt.close(fig)
    return out


# --------------------------------------------------------------------------- G1

def g1(rows, cal, figs):
    """The photograph track: how the phenotype's statistics scale with tau_m."""
    tc = cal["tau_c"]
    print("\nG1 -- the phenotype as a blurred copy of the pressure: scaling with tau_m")
    print(f"{'tau_m/tau_c':>11} {'std(chi)':>10} {'std(m)':>9} {'L_chi':>8} {'L_P':>7} "
          f"{'tau_auto/tau_c':>15} {'lag peak':>9} {'at Delta/tau_c':>15}")
    xs, sd, lc, ta = [], [], [], []
    for tm, g in by(rows, lambda r: round(r["params"]["tau_m"] / tc, 2)).items():
        s = np.mean([r["flow"]["std_chi"] for r in g])
        l = np.mean([r["flow"]["L_chi"] for r in g])
        t = np.mean([(r["time"] or {}).get("tau_auto_chi", np.nan) for r in g]) / tc
        pk = np.mean([(r["time"] or {}).get("lag_peak", np.nan) for r in g])
        dl = np.mean([(r["time"] or {}).get("lag_peak_delta", np.nan) for r in g]) / tc
        xs.append(tm); sd.append(s); lc.append(l); ta.append(t)
        print(f"{tm:11.1f} {s:10.4f} {np.mean([r['flow']['std_m'] for r in g]):9.4f} "
              f"{l:8.2f} {np.mean([r['flow']['L_P'] for r in g]):7.2f} {t:15.2f} "
              f"{pk:9.3f} {dl:15.2f}")
    for name, ys in (("std_chi", sd), ("L_chi", lc), ("tau_auto", ta)):
        good = [(x, y) for x, y in zip(xs, ys) if np.isfinite(y) and y > 0 and x > 0]
        if len(good) >= 3:
            p = np.polyfit(np.log([g[0] for g in good]), np.log([g[1] for g in good]), 1)
            print(f"    {name} ~ (tau_m/tau_c)^{p[0]:+.3f}")
    if figs and xs:
        fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), dpi=140)
        for ax, ys, lab in zip(axes, (sd, lc, ta),
                               (r"std$(\chi)$", r"$L_\chi$ [cells]", r"$\tau_{auto}/\tau_c$")):
            ax.loglog(xs, ys, "o-", color="#1f4e79")
            ax.set_xlabel(r"$\tau_m/\tau_c$"); ax.set_ylabel(lab)
        fig.suptitle("G1: the phenotype's statistics vs the memory time", fontsize=10)
        fig.tight_layout(); fig.savefig(os.path.join(figs, "g1_scaling.png")); plt.close(fig)
    return {"tau_m_over_tau_c": xs, "std_chi": sd, "L_chi": lc, "tau_auto_over_tau_c": ta}


# --------------------------------------------------------------------------- G5

def g5(groups, cal, g3res):
    """Each control removes one alternative explanation of the G3 result."""
    tc = cal["tau_c"]
    print("\nG5 -- controls")
    out = {}
    for arm, rows in sorted(groups.items()):
        if not arm.startswith("cw_mem_g5"):
            continue
        print(f"  {arm}:")
        for key, g in by(rows, lambda r: (round(r["params"]["tau_m"] / tc, 2),
                                          round(r["params"]["tau_chi"] / tc, 3),
                                          round(r["params"]["chi_width"], 4),
                                          round(r["params"]["Dbio"], 6),
                                          r["params"]["switch_sign"])).items():
            tm, tchi, w, db, sgn = key
            a = [r for r in g if r["params"]["chi0"] == 0.0]
            b = [r for r in g if r["params"]["chi0"] == 1.0]
            if not a or not b:
                continue
            ta, tb = np.mean([tail(r) for r in a]), np.mean([tail(r) for r in b])
            base = dict(zip(g3res["tau_m_over_tau_c"], g3res["separation"])).get(tm)
            cmp = f"  (G3 at the same tau_m: {base:.3f})" if base is not None else ""
            print(f"    tau_m={tm:5.1f} tau_chi={tchi:5.2f} w_chi={w:.4f} Dbio={db:.5f} "
                  f"s={sgn:+d}:  {ta:.4f} vs {tb:.4f}  separation {abs(tb-ta):.4f}{cmp}")
            out.setdefault(arm, []).append({"tau_m": tm, "tau_chi": tchi, "chi_width": w,
                                            "Dbio": db, "switch_sign": sgn,
                                            "chi0_arm": ta, "chi1_arm": tb,
                                            "separation": abs(tb - ta),
                                            "g3_separation": base})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--calib", required=True)
    ap.add_argument("--out", default=None)
    ap.add_argument("--figs", default=None)
    a = ap.parse_args()
    cal = json.load(open(a.calib))
    groups = load(a.root)
    if a.figs:
        os.makedirs(a.figs, exist_ok=True)
    print(f"{a.root}: " + ", ".join(f"{k} {len(v)}" for k, v in sorted(groups.items())))

    res = {}
    if "cw_mem_g3" in groups:
        res["g3"] = g3(groups["cw_mem_g3"], cal, a.figs)
    if "cw_mem_g4" in groups:
        res["g4"] = {str(k): v for k, v in g4(groups["cw_mem_g4"], cal, a.figs).items()}
    if "cw_mem_g1" in groups:
        res["g1"] = g1(groups["cw_mem_g1"], cal, a.figs)
    if any(k.startswith("cw_mem_g5") for k in groups):
        res["g5"] = g5(groups, cal, res.get("g3", {"tau_m_over_tau_c": [], "separation": []}))
    if a.out:
        with open(a.out, "w") as fh:
            json.dump(res, fh, indent=1, default=float)
        print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
