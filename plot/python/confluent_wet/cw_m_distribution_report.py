#!/usr/bin/env python3
"""Render diagnostics and a factual report from cw_m_distribution.py outputs.

Usage: conda run --no-capture-output -n env1 python \
    plot/python/confluent_wet/cw_m_distribution_report.py INPUT_ROOT OUTPUT_DIR

Reads distribution.json recursively; does not read simulation frames or refit data.
Histogram densities are derived by interpolating the stored empirical CDF onto
fixed bins. Confidence intervals are copied from the analysis, not recalculated.
"""

import argparse
import csv
import json
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import numpy as np
from scipy.special import betainc, betaincc, ndtr
from scipy.stats import beta as beta_dist
from scipy.stats import norm


COLORS = {"data": "#222222", "gaussian": "#0072B2", "beta": "#D55E00"}
HIST_EDGES = np.unique(np.r_[0., np.geomspace(1e-5, .02, 16),
                             np.linspace(0., 1., 201),
                             1. - np.geomspace(1e-5, .02, 16), 1.])
REPRESENTATIVES = ("tm0p3_chi0", "tm3p43_chi0", "tm9p68_chi0",
                   "tm12_chi0", "tm20_chi0", "tm20_chi1")
METHODS = {
    "precision": "Read serialized floating-point m fields (approximately six significant digits), not 8-bit video. Reported exact endpoint counts refer to stored values and need not imply physical point masses.",
    "target": "Per-case equal-area, equal-saved-frame marginal distribution of m.",
    "window": "Regular saved frames inside the final 100 tau_c of each simulation; the present campaign uses nine frames at approximately t/tau_c = 201.7-296.6, separated by 11.86 tau_c. An additional non-regular end frame is excluded. Actual frame times are retained in each source JSON and diagnostic page.",
    "fits": "Gaussian and Beta use the same empirical mean and variance. Beta support is fixed to [0,1].",
    "grid_ks": "Supremum CDF error on the saved grid, with the 0-left and 1-left limits included; not an exact continuous-grid KS statistic.",
    "difference": "Delta D = D_Gaussian - D_Beta; positive values favour Beta.",
    "intervals": "Pointwise 95% spatial column-block bootstrap intervals, with all selected frames retained within each block and both models refitted in every replicate. These intervals condition on the saved window; residual dependence across blocks is not fully resolved.",
    "block_sizes": "25, 50 and 100 lattice cells; the overview uses 50. Block-size sensitivity is retained in the tables.",
    "crossfit": "Five 100-cell spatial strips; 25-cell periodic training buffer. Held-out binned log score is Beta minus Gaussian. The bootstrap resamples 50-cell column blocks and refits the training moments. Gaussian probabilities are not renormalised to [0,1].",
    "histogram": "Grid-derived histogram: the saved empirical CDF is linearly interpolated onto fixed bins. No smoothing is applied. Model density lines show bin-averaged probabilities on the same bins. Endpoint atoms are included in the first/last bin.",
    "drift_flags": "Drift-flagged window / stationarity not established: preset diagnostics compare early/late window mean, standard deviation, chi and CDF. These flags are not a formal rejection of stationarity; absence of flags does not establish stationarity.",
    "interpretation": "Intervals compare the approximation errors of two moment closures. Better agreement does not establish an exact Beta or Gaussian law. No iid KS p-values or absolute goodness-of-fit claims are reported; intervals are pointwise and are not adjusted for multiple comparisons.",
}


def set_style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.labelsize": 10, "axes.titlesize": 10,
        "xtick.labelsize": 9, "ytick.labelsize": 9,
        "legend.fontsize": 9, "axes.spines.top": False,
        "axes.spines.right": False, "axes.linewidth": .7,
        "lines.linewidth": 1.4, "savefig.dpi": 200,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def load_cases(root):
    cases = []
    for path in sorted(root.rglob("distribution.json")):
        with path.open() as stream:
            item = json.load(stream)
        required = ("case", "tau_m", "chi0", "mean", "std", "alpha", "beta",
                    "edges", "cdf", "gaussian_cdf", "beta_cdf", "bootstrap",
                    "crossfit_logscore", "per_frame", "ks_gaussian", "ks_beta")
        missing = [name for name in required if name not in item]
        if missing:
            raise ValueError(f"{path}: missing {missing}")
        edges = np.asarray(item["edges"], float)
        empirical = np.asarray(item["cdf"], float)
        if (len(edges) != len(empirical) or np.any(np.diff(edges) <= 0)
                or np.any(np.diff(empirical) < -1e-10)):
            raise ValueError(f"{path}: invalid CDF grid")
        item["source_json"] = str(path.resolve())
        cases.append(item)
    if not cases:
        raise ValueError(f"No distribution.json found below {root}")
    return sorted(cases, key=lambda d: (d["tau_m"], d["chi0"],
                                       d.get("study", ""), d["case"]))


def number(value, digits=4):
    return f"{float(value):.{digits}g}" if value is not None and np.isfinite(value) else "NA"


def ci_text(result, digits=4):
    point = number(result["estimate"], digits)
    low, high = (number(x, digits) for x in result["ci95"])
    return f"{point} [{low}, {high}]"


def interval_class(result):
    lo, hi = result["ci95"]
    return "Beta" if lo > 0 else ("Gaussian" if hi < 0 else "unresolved")


def flagged(d):
    return bool(d.get("nonstationary_flags", []))


def model_cdf(d, x):
    return (ndtr((np.asarray(x) - d["mean"]) / d["std"]),
            betainc(d["alpha"], d["beta"], x))


def bin_probabilities(cdf, survival):
    # Use survival differences in the upper tail to avoid catastrophic cancellation.
    return np.maximum(np.where(cdf[1:] < .5, np.diff(cdf), -np.diff(survival)), 0.)


def densities(d):
    empirical = np.interp(HIST_EDGES, d["edges"], d["cdf"])
    probability = np.diff(empirical)
    probability[0] += empirical[0]  # include m == 0, which belongs to the first bin
    widths = np.diff(HIST_EDGES)
    z = (HIST_EDGES - d["mean"]) / d["std"]
    gauss = bin_probabilities(ndtr(z), ndtr(-z))
    beta = bin_probabilities(betainc(d["alpha"], d["beta"], HIST_EDGES),
                             betaincc(d["alpha"], d["beta"], HIST_EDGES))
    return np.maximum(probability, 0.) / widths, gauss / widths, beta / widths


def view_limits(d):
    """Show the central distributions and mc; account explicitly for cropped mass."""
    empirical = np.asarray(d["cdf"])
    edges = np.asarray(d["edges"])
    low_data = edges[max(0, np.searchsorted(empirical, .0005) - 1)]
    high_data = edges[min(len(edges) - 1, np.searchsorted(empirical, .9995))]
    low_model = min(norm.ppf(.0005, d["mean"], d["std"]),
                    beta_dist.ppf(.0005, d["alpha"], d["beta"]))
    high_model = max(norm.ppf(.9995, d["mean"], d["std"]),
                     beta_dist.ppf(.9995, d["alpha"], d["beta"]))
    mc = d.get("mc", .21)
    lo = max(0., min(low_data, low_model, mc) - .015)
    hi = min(1., max(high_data, high_model, mc) + .015)
    if hi - lo > .8:
        return 0., 1.
    return float(lo), float(hi)


def outside_view(d, lo, hi):
    # Display interval includes its endpoints. At lo == 0 the zero atom is visible.
    empirical = (0. if lo == 0 else np.interp(lo, d["edges"], d["cdf"]))
    empirical += 1. - np.interp(hi, d["edges"], d["cdf"])
    g = ndtr((lo - d["mean"]) / d["std"]) + ndtr((d["mean"] - hi) / d["std"])
    b = betainc(d["alpha"], d["beta"], lo) + betaincc(d["alpha"], d["beta"], hi)
    return empirical, g, b


def plot_density(ax, d, annotate=True):
    observed, gaussian, beta = densities(d)
    centres = (HIST_EDGES[:-1] + HIST_EDGES[1:]) / 2
    lo, hi = view_limits(d)
    ax.stairs(observed, HIST_EDGES, color=COLORS["data"], lw=1.2,
              label="Empirical (grid-derived)")
    ax.plot(centres, gaussian, color=COLORS["gaussian"], label="Gaussian (same moments)")
    ax.plot(centres, beta, color=COLORS["beta"], label="Beta (same moments)")
    ax.axvline(d.get("mc", .21), color=".4", ls="--", lw=1)
    ax.set(xlim=(lo, hi), xlabel="$m$", ylabel="Binned density (log scale)", yscale="log")
    visible = (centres >= lo) & (centres <= hi)
    positive = observed[visible & (observed > 0)]
    ymax = max(np.max(y[visible]) for y in (observed, gaussian, beta))
    ymin = max(1e-5, min(.01, .2 * positive.min())) if len(positive) else 1e-4
    ax.set_ylim(ymin, max(ymin * 100, ymax * 1.4))
    if annotate:
        omitted = outside_view(d, lo, hi)
        ax.text(.99, .025, "Outside x view (data/G/B): "
                + "/".join(f"{100*x:.2g}%" for x in omitted),
                transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": .8, "pad": 1})


def model_title(d):
    return (rf"$\tau_m/\tau_c={d['tau_m']:.3g},\ \chi_0={d['chi0']:.1f}$"
            + f"  ({d['case']})")


def representative_figure(cases, out):
    fig, axes = plt.subplots(2, 3, figsize=(14, 8.9))
    found = []
    for ax, name, letter in zip(axes.ravel(), REPRESENTATIVES, "abcdef"):
        matches = [d for d in cases if d["case"] == name]
        if not matches:
            ax.text(.5, .5, f"{name}\nNo input result", transform=ax.transAxes,
                    ha="center", va="center")
            ax.set_axis_off()
            continue
        d = matches[0]
        found.append(name)
        plot_density(ax, d)
        flag_text = " | drift-flagged" if flagged(d) else ""
        ax.set_title(f"({letter}) " + model_title(d) + flag_text + "\n"
                     + f"Grid D (G/B): {d['ks_gaussian']:.4f}/{d['ks_beta']:.4f}\n"
                     + "Pr(m < mc), data/G/B: "
                     + "/".join(f"{d[k]:.4f}" for k in
                                  ("cdf_mc", "cdf_mc_gaussian", "cdf_mc_beta")),
                     fontsize=9, pad=7)
    handles = [Line2D([], [], color=COLORS[key], label=label) for key, label in
               (("data", "Empirical: grid-derived histogram"),
                ("gaussian", "Moment-matched Gaussian"),
                ("beta", "Moment-matched Beta"))]
    handles.append(Line2D([], [], color=".4", ls="--", label="$m_c$"))
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, .945),
               ncol=4, frameon=False, fontsize=9)
    fig.suptitle("Memory distributions: matched moments, measured shape", fontsize=13, y=.985)
    fig.text(.5, .018,
             "Fixed-bin densities from CDF interpolation; no smoothing. Model lines are bin averages. "
             "Log y axes show tails; cropped x mass is disclosed in each panel.\n"
             "Drift flags are preset window diagnostics, not a formal stationarity test.",
             ha="center", va="bottom", fontsize=8)
    fig.tight_layout(rect=(0, .075, 1, .92), h_pad=2.6, w_pad=2)
    for ext in ("png", "pdf"):
        fig.savefig(out / f"representative_distributions.{ext}")
    plt.close(fig)
    return found


def representative_cdfs(cases, out):
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.1), sharey=True)
    for ax, name, letter in zip(axes.ravel(), REPRESENTATIVES, "abcdef"):
        matches = [d for d in cases if d["case"] == name]
        if not matches:
            ax.text(.5, .5, f"{name}\nNo input result", transform=ax.transAxes,
                    ha="center", va="center")
            ax.set_axis_off()
            continue
        d = matches[0]
        x = np.asarray(d["edges"])
        ax.plot(x, d["cdf"], color=COLORS["data"])
        ax.plot(x, d["gaussian_cdf"], color=COLORS["gaussian"])
        ax.plot(x, d["beta_cdf"], color=COLORS["beta"])
        mc = d.get("mc", .21)
        ax.axvline(mc, color=".4", ls="--", lw=1)
        for key, colour in (("cdf_mc", "data"), ("cdf_mc_gaussian", "gaussian"),
                            ("cdf_mc_beta", "beta")):
            ax.plot(mc, d[key], "o", ms=3.5, color=COLORS[colour])
        ax.set(xlim=view_limits(d), ylim=(-.02, 1.02), xlabel="$m$", ylabel="CDF")
        status = " | drift-flagged" if flagged(d) else ""
        ax.set_title(f"({letter}) " + model_title(d) + status + "\n"
                     + f"Grid D (G/B): {d['ks_gaussian']:.4f}/{d['ks_beta']:.4f}",
                     fontsize=9, pad=7)
        if name == "tm0p3_chi0":
            inset = ax.inset_axes([.48, .10, .47, .44])
            for key, colour in (("cdf", "data"), ("gaussian_cdf", "gaussian"),
                                ("beta_cdf", "beta")):
                inset.plot(x[1:], np.asarray(d[key])[1:], color=COLORS[colour], lw=1)
            inset.set(xscale="log", xlim=(1e-8, 1e-2), ylim=(0, .65))
            inset.set_title("Near-zero tail", fontsize=8)
            inset.set_xticks([1e-8, 1e-5, 1e-2])
            inset.tick_params(labelsize=7)
    handles = [Line2D([], [], color=COLORS[key], label=label) for key, label in
               (("data", "Empirical CDF"), ("gaussian", "Moment-matched Gaussian"),
                ("beta", "Moment-matched Beta"))]
    handles.append(Line2D([], [], color=".4", ls="--", label="$m_c$"))
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, .95),
               ncol=4, frameon=False, fontsize=9)
    fig.suptitle("Memory CDFs: direct comparison of the closure error", fontsize=13, y=.99)
    fig.text(.5, .018,
             "Dots mark Pr(m < mc) and both closure predictions. Grid D is computed over the full saved CDF grid, "
             "including endpoint limits.\n"
             "Displayed x ranges follow each distribution and include mc; full-support CDFs are in all_cases.pdf. "
             "No absolute goodness-of-fit claim.",
             ha="center", va="bottom", fontsize=8)
    fig.tight_layout(rect=(0, .075, 1, .92), h_pad=2.5, w_pad=2)
    for ext in ("png", "pdf"):
        fig.savefig(out / f"representative_cdfs.{ext}")
    plt.close(fig)


def all_case_pages(cases, out):
    with PdfPages(out / "all_cases.pdf") as pdf:
        for d in cases:
            fig, axes = plt.subplots(2, 3, figsize=(12, 8.4))
            plot_density(axes[0, 0], d)
            axes[0, 0].set_title("Grid-derived histogram and closures")
            axes[0, 0].legend(frameon=False, fontsize=7, loc="upper right")
            x = np.asarray(d["edges"])
            h = np.asarray(d["cdf"])
            g = np.asarray(d["gaussian_cdf"])
            b = np.asarray(d["beta_cdf"])
            ax = axes[0, 1]
            ax.plot(x, h, color=COLORS["data"], label="Empirical CDF")
            ax.plot(x, g, color=COLORS["gaussian"], label="Gaussian")
            ax.plot(x, b, color=COLORS["beta"], label="Beta")
            ax.axvline(d.get("mc", .21), color=".4", ls="--", lw=1)
            ax.set(xlim=(0, 1), ylim=(-.02, 1.02), xlabel="$m$", ylabel="CDF",
                   title="Full-support CDF; tails retained")
            ax.legend(frameon=False, fontsize=8)
            ax = axes[0, 2]
            ax.plot(x, h - g, color=COLORS["gaussian"], label="ECDF - Gaussian")
            ax.plot(x, h - b, color=COLORS["beta"], label="ECDF - Beta")
            ax.axhline(0, color=".65", lw=.8)
            ax.axvline(d.get("mc", .21), color=".4", ls="--", lw=1)
            ax.set(xlim=(0, 1), xlabel="$m$", ylabel="CDF residual",
                   title=f"Grid D (G/B): {d['ks_gaussian']:.4f}/{d['ks_beta']:.4f}")
            ax.legend(frameon=False, fontsize=8)
            per = d["per_frame"]
            times = [p["t"] for p in per]
            for ax, key, label in zip(axes[1], ("mean", "std", "chi"),
                                      (r"$\langle m\rangle_x$", r"$\sigma_{m,x}$",
                                       r"$\langle\chi\rangle_x$")):
                ax.plot(times, [p[key] for p in per], "o-", color=COLORS["data"], ms=4)
                ax.set(xlabel=r"$t/\tau_c$", ylabel=label, title="Saved-frame diagnostic")
                if key == "chi":
                    ax.set_ylim(-.03, 1.03)
            status = ("Drift-flagged window; stationarity not established: "
                      + "; ".join(d.get("nonstationary_flags", []))) if flagged(d) else (
                      "No preset drift flag; stationarity is not established by this diagnostic.")
            stats = (f"Window: {d['window'][0]:.1f}-{d['window'][1]:.1f} tau_c; "
                     f"{d['frames']} saved frames.  "
                     f"Pr(m < mc), data/G/B: {d['cdf_mc']:.5f}/"
                     f"{d['cdf_mc_gaussian']:.5f}/{d['cdf_mc_beta']:.5f}.  "
                     f"Atoms at 0/1: {d.get('edge_zero', 0):.2g}/{d.get('edge_one', 0):.2g}.")
            ci = ("Delta grid D, 95% CI (50-cell blocks): " + ci_text(d["bootstrap"]["50"])
                  + "; held-out log-score difference: " + ci_text(d["crossfit_logscore"]) + ".")
            fig.suptitle(model_title(d) + "\n" + d.get("study", ""), fontsize=12, y=.985)
            fig.text(.035, .1, stats + "\n" + ci + "\n" + textwrap.fill(status, 155),
                     fontsize=8, va="top")
            fig.text(.035, .015,
                     "Histogram: fixed bins from CDF interpolation; no smoothing. "
                     "Bootstrap intervals condition on the saved window; cross-block dependence remains approximate.",
                     fontsize=7.5)
            fig.tight_layout(rect=(0, .125, 1, .92), h_pad=2, w_pad=2)
            pdf.savefig(fig)
            plt.close(fig)


def overview_figure(cases, out):
    height = max(5.5, .29 * len(cases) + 2.5)
    fig, axes = plt.subplots(1, 2, figsize=(12, height), sharey=True)
    for ax, metric, title, xlabel in (
        (axes[0], "ks", "CDF approximation error", r"$\Delta D=D_G-D_B$ (grid KS)"),
        (axes[1], "score", "Held-out predictive score", "Mean log-score: Beta - Gaussian [nats]"),
    ):
        for i, d in enumerate(cases):
            result = d["bootstrap"]["50"] if metric == "ks" else d["crossfit_logscore"]
            lo, hi = result["ci95"]
            value = result["estimate"]
            ax.plot([lo, hi], [i, i], color=".35", lw=1.2, zorder=2)
            ax.plot(value, i, "o", ms=5, color=".2", mfc="white" if flagged(d) else ".2",
                    markeredgewidth=1.1, zorder=3)
            if i and abs(d["tau_m"] - cases[i - 1]["tau_m"]) > .01:
                ax.axhline(i - .5, color=".84", lw=.7)
        ax.axvline(0, color=".5", ls="--", lw=1)
        ax.set(title=title + "\nPointwise 95% bootstrap CI; positive favours Beta",
               xlabel=xlabel, ylim=(len(cases) - .3, -.7))
        ax.tick_params(axis="y", length=0)
        ax.grid(axis="x", alpha=.12)
        if metric == "score":
            ax.set_xscale("symlog", linthresh=.01)
            ax.set_xlabel("Mean log-score: Beta - Gaussian [nats]\nSymlog axis; linear within +/-0.01")
    axes[0].set_yticks(range(len(cases)),
                      [f"{d['tau_m']:.3g}  |  {d['chi0']:.1f}" for d in cases])
    axes[0].set_ylabel(r"Case: $\tau_m/\tau_c$  |  $\chi_0$")
    handles = [Line2D([], [], marker="o", color=".2", ls="none", mfc=".2",
                      label="No preset drift flag"),
               Line2D([], [], marker="o", color=".2", ls="none", mfc="white",
                      label="Drift-flagged window; stationarity not established")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .047),
               ncol=2, frameon=False, fontsize=9)
    fig.suptitle(f"Gaussian versus Beta memory closure: {len(cases)} cases", fontsize=13, y=.99)
    fig.text(.5, .016,
             "50-cell spatial column blocks retain all saved frames. CIs condition on the saved window; "
             "cross-block dependence is not fully resolved.\n"
             "No iid KS p-values or absolute goodness-of-fit claims. Intervals are pointwise, without multiplicity adjustment.",
             ha="center", va="bottom", fontsize=8)
    fig.tight_layout(rect=(.015, .085, 1, .965), w_pad=3.5)
    for ext in ("png", "pdf"):
        fig.savefig(out / f"overview.{ext}")
    plt.close(fig)


def flat_record(d):
    result = {key: d.get(key) for key in
              ("study", "case", "tau_m", "chi0", "frames", "sites_per_frame", "mean", "std",
               "alpha", "beta", "skewness", "excess_kurtosis", "cdf_mc", "cdf_mc_gaussian",
               "cdf_mc_beta", "ks_gaussian", "ks_beta", "gaussian_outside_support",
               "edge_zero", "edge_one", "temporal_ks", "drift_mean_over_sd", "drift_chi",
               "m_corr_1e", "m_corr_50", "bootstrap_reps", "source_json")}
    result.update(window_start=d["window"][0], window_end=d["window"][1],
                  cdf_mc_error_gaussian=d["cdf_mc_gaussian"] - d["cdf_mc"],
                  cdf_mc_error_beta=d["cdf_mc_beta"] - d["cdf_mc"],
                  drift_flagged=flagged(d),
                  drift_diagnostics="; ".join(d.get("nonstationary_flags", [])))
    for side in ("25", "50", "100"):
        block = d["bootstrap"][side]
        result[f"delta_grid_ks_{side}"] = block["estimate"]
        result[f"delta_grid_ks_ci_low_{side}"], result[f"delta_grid_ks_ci_high_{side}"] = block["ci95"]
        result[f"blocks_{side}"] = block["blocks"]
        threshold = block.get("mc_abs_error_difference")
        if threshold is not None:
            result[f"delta_mc_abs_error_{side}"] = threshold["estimate"]
            result[f"delta_mc_abs_error_ci_low_{side}"], result[f"delta_mc_abs_error_ci_high_{side}"] = threshold["ci95"]
    score = d["crossfit_logscore"]
    result["heldout_logscore_delta"] = score["estimate"]
    result["heldout_logscore_ci_low"], result["heldout_logscore_ci_high"] = score["ci95"]
    result["grid_ks_comparison_ci50"] = interval_class(d["bootstrap"]["50"])
    result["logscore_comparison_ci50"] = interval_class(score)
    result["grid_ks_ci_class_consistent_25_50_100"] = (
        len({interval_class(d["bootstrap"][k]) for k in ("25", "50", "100")}) == 1)
    return result


def clean_json(value):
    if isinstance(value, dict):
        return {k: clean_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean_json(x) for x in value]
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return None
    return value


def write_tables(cases, out, representatives):
    rows = [flat_record(d) for d in cases]
    counts = {
        "cases": len(cases), "drift_flagged_windows": sum(flagged(d) for d in cases),
        "frames_min": min(d["frames"] for d in cases),
        "frames_max": max(d["frames"] for d in cases),
        "window_start_range": [min(d["window"][0] for d in cases),
                               max(d["window"][0] for d in cases)],
        "window_end_range": [min(d["window"][1] for d in cases),
                             max(d["window"][1] for d in cases)],
        "median_frame_spacing_tau_c": float(np.median(
            [gap for d in cases for gap in np.diff([p["t"] for p in d["per_frame"]])])),
        "grid_ks_ci50": {name: sum(interval_class(d["bootstrap"]["50"]) == name for d in cases)
                         for name in ("Beta", "Gaussian", "unresolved")},
        "logscore_ci50": {name: sum(interval_class(d["crossfit_logscore"]) == name for d in cases)
                          for name in ("Beta", "Gaussian", "unresolved")},
        "block_size_ci_class_changes": sum(not r["grid_ks_ci_class_consistent_25_50_100"] for r in rows),
    }
    summary = clean_json({"methods": METHODS, "counts": counts,
                          "representative_cases_found": representatives,
                          "representative_cases_missing": [n for n in REPRESENTATIVES if n not in representatives],
                          "cases": rows})
    (out / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    with (out / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(clean_json(rows))
    lines = ["# Memory distribution comparison", "",
             "This report compares moment-matched Gaussian and Beta closures case by case. "
             "It does not claim either is the exact data-generating distribution.", "",
             "## Methods and limits", ""]
    for name, description in METHODS.items():
        lines.append(f"- **{name.replace('_', ' ').capitalize()}:** {description}")
    lines.extend(["", "## Observed comparison counts", "",
                  f"- Cases: {counts['cases']}; saved frames per case: "
                  f"{counts['frames_min']}-{counts['frames_max']}.",
                  f"- Actual saved-window start range: {counts['window_start_range']}; "
                  f"end range: {counts['window_end_range']} tau_c. Median frame spacing: "
                  f"{counts['median_frame_spacing_tau_c']:.4g} tau_c.",
                  f"- Drift-flagged windows: {counts['drift_flagged_windows']}. "
                  "Stationarity is not established; these are preset diagnostic flags, not formal tests.",
                  "- Pointwise interval classifications count cases descriptively; they are not independent "
                  "replicate experiments and do not provide a multiplicity-controlled campaign conclusion.",
                  f"- Grid KS CI50: {counts['grid_ks_ci50']}.",
                  f"- Held-out log-score CI50: {counts['logscore_ci50']}.",
                  f"- Cases whose grid KS interval classification changes across 25/50/100-cell blocks: "
                  f"{counts['block_size_ci_class_changes']}.", "", "## Case facts", "",
                  "Positive Delta D or log-score difference favours Beta; negative favours Gaussian. "
                  "Brackets give pointwise 95% CIs from 50-cell blocks. 'Flag' denotes a drift-flagged window.", "",
                  "| Case | tau_m/tau_c | chi0 | mean m | std m | Grid D G/B | Delta D [95% CI] | Log-score B-G [95% CI] | Window |",
                  "|---|---:|---:|---:|---:|---:|---|---|---|"])
    for d in cases:
        lines.append(f"| {d['case']} | {d['tau_m']:.3g} | {d['chi0']:.1f} | {d['mean']:.5g} "
                     f"| {d['std']:.4g} | {d['ks_gaussian']:.4g}/{d['ks_beta']:.4g} "
                     f"| {ci_text(d['bootstrap']['50'])} | {ci_text(d['crossfit_logscore'])} "
                     f"| {'Flag' if flagged(d) else 'No preset flag'} |")
    lines.extend(["", "## Drift diagnostics", ""])
    drift_rows = [d for d in cases if flagged(d)]
    lines.extend([f"- **{d['case']}:** {'; '.join(d['nonstationary_flags'])}." for d in drift_rows]
                 or ["No case crossed a preset drift diagnostic threshold."])
    lines.extend(["", "## Outputs", "",
                  "- [Overview](overview.png): both effect sizes and 95% intervals for every case.",
                  "- [Representative distributions](representative_distributions.png): six selected cases; "
                  "all density axes are logarithmic and omitted x-range mass is printed.",
                  "- [Representative CDFs](representative_cdfs.png): direct comparison of the empirical CDF "
                  "with both closures, with threshold probabilities marked.",
                  "- [All-case diagnostics](all_cases.pdf): one page per case, including full-support CDFs, "
                  "residuals and saved-frame time series.",
                  "- [Machine-readable table](summary.csv) and [structured summary](summary.json): "
                  "all three block sizes, threshold errors and source paths. No p-values are exported.", ""])
    if summary["representative_cases_missing"]:
        lines.append("Missing requested representative cases: " + ", ".join(summary["representative_cases_missing"]) + ".")
    (out / "report.md").write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputdir", type=Path)
    parser.add_argument("outdir", type=Path)
    args = parser.parse_args()
    set_style()
    cases = load_cases(args.inputdir)
    args.outdir.mkdir(parents=True, exist_ok=True)
    representatives = representative_figure(cases, args.outdir)
    representative_cdfs(cases, args.outdir)
    overview_figure(cases, args.outdir)
    all_case_pages(cases, args.outdir)
    write_tables(cases, args.outdir, representatives)
    print(f"Wrote {len(cases)} case pages, summary tables, report and figures to {args.outdir}", flush=True)


if __name__ == "__main__":
    main()
