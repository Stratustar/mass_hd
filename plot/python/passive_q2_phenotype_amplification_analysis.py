import matplotlib
matplotlib.use("Agg")

import argparse
import csv
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager


def configure_plot_style():
    arial_font = os.environ.get("MASS_HD_ARIAL_FONT")
    if arial_font and os.path.exists(arial_font):
        font_manager.fontManager.addfont(arial_font)
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial"],
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "axes.titlepad": 6,
            "axes.grid": False,
        }
    )


configure_plot_style()


SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from archive.archive import loadarchive


PHI_THRESHOLD = 0.5
PHI_VISIBLE_THRESHOLD = 0.05
BOUNDARY_COLOR = "#d9468f"
BACKGROUND_COLOR = "#7e7e7e"
COLORS = {
    (0, "weak"): "#255c99",
    (0, "strong"): "#b45309",
    (1, "weak"): "#047857",
    (1, "strong"): "#7c3aed",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Q2 passive proliferation analysis: phenotype amplification."
    )
    parser.add_argument(
        "input_root",
        help="Directory containing B1/B2/B3/B4 simulation archive directories.",
    )
    parser.add_argument("outdir", help="Directory for CSV and PNG outputs.")
    parser.add_argument(
        "--dpi",
        type=int,
        default=220,
        help="PNG resolution. Default: 220.",
    )
    parser.add_argument(
        "--frame-stride",
        type=int,
        default=1,
        help="Analyze every Nth frame. Default: 1.",
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=50,
        help="Number of chi histogram bins. Default: 50.",
    )
    return parser.parse_args()


def archive_marker(files):
    return "parameters.json" in files or "parameters.json.zip" in files


def discover_archives(input_root):
    input_root = os.path.abspath(input_root)
    archives = []
    for root, _, files in os.walk(input_root):
        if archive_marker(files):
            rel = os.path.relpath(root, input_root)
            label = os.path.basename(root) if rel == "." else rel
            archives.append((label, root))
    return sorted(archives, key=lambda item: item[0])


def available_frame_count(ar):
    return int((ar.nsteps - ar.nstart) / ar.ninfo) + 1


def frame_indices(ar, stride):
    n_available = available_frame_count(ar)
    stride = max(1, int(stride))
    indices = list(range(0, n_available, stride))
    if indices[-1] != n_available - 1:
        indices.append(n_available - 1)
    return indices


def simulation_time(ar, frame_index):
    step = ar.nstart + frame_index * ar.ninfo
    return step * float(getattr(ar, "time_step", 1.0))


def param(ar, name, default=np.nan):
    return ar.parameters.get(name, default)


def classify_alpha(alpha):
    if not np.isfinite(alpha):
        return "unknown"
    return "weak" if alpha < 0.001 else "strong"


def display_label(case):
    grow = int(case["growTogether"])
    alpha_label = classify_alpha(case["alpha"])
    return rf"gT={grow}, {alpha_label}, $\alpha={case['alpha']:g}$"


def grid(frame, name):
    lx = frame.parameters["LX"]
    ly = frame.parameters["LY"]
    return np.asarray(getattr(frame, name), dtype=float).reshape((lx, ly))


def optional_grid(frame, name, fallback):
    if hasattr(frame, name):
        return grid(frame, name)
    return np.array(fallback, dtype=float, copy=True)


def field_set(frame):
    phi = grid(frame, "phi")
    m = optional_grid(frame, "m", np.zeros_like(phi))
    chi_raw = optional_grid(frame, "chi", np.full_like(phi, np.nan))
    chi = np.divide(m, phi, out=np.zeros_like(phi), where=phi > PHI_VISIBLE_THRESHOLD)
    if np.any(np.isfinite(chi_raw)):
        chi = chi_raw
    chi = np.clip(chi, 0.0, 1.0)
    return {"phi": phi, "m": m, "chi": chi}


def weighted_mean(values, weights):
    weight_sum = float(np.sum(weights))
    if weight_sum <= 0.0:
        return np.nan
    return float(np.sum(values * weights) / weight_sum)


def weighted_var(values, weights, mean):
    weight_sum = float(np.sum(weights))
    if weight_sum <= 0.0 or not np.isfinite(mean):
        return np.nan
    return float(np.sum(weights * (values - mean) ** 2) / weight_sum)


def pearson_corr(x, y):
    if len(x) < 2:
        return np.nan
    x0 = x - np.mean(x)
    y0 = y - np.mean(y)
    denom = float(np.sqrt(np.sum(x0 * x0) * np.sum(y0 * y0)))
    if denom <= 0.0:
        return np.nan
    return float(np.sum(x0 * y0) / denom)


def chi_bin_phi_curve(phi, chi, material, nbins):
    bins = np.linspace(0.0, 1.0, nbins + 1)
    centers = 0.5 * (bins[:-1] + bins[1:])
    which = np.digitize(chi[material], bins) - 1
    which = np.clip(which, 0, nbins - 1)
    sums = np.bincount(which, weights=phi[material], minlength=nbins)
    counts = np.bincount(which, minlength=nbins)
    means = np.divide(
        sums,
        counts,
        out=np.full(nbins, np.nan),
        where=counts > 0,
    )
    return centers, means


def frame_metrics(ar, label, frame_index, first_mass):
    fields = field_set(ar.read_frame(frame_index))
    phi = fields["phi"]
    chi = fields["chi"]
    material = phi > PHI_THRESHOLD
    weights = np.where(material, np.maximum(phi, 0.0), 0.0)

    chi_bar = weighted_mean(chi, weights)
    var_chi = weighted_var(chi, weights, chi_bar)
    h_chi = weighted_mean(chi * (1.0 - chi), weights)
    corr_phi_chi = pearson_corr(phi[material], chi[material]) if np.any(material) else np.nan
    high = material & (chi >= 0.75)
    low = material & (chi <= 0.25)

    return {
        "case": label,
        "growTogether": int(param(ar, "growTogether", -1)),
        "alpha": float(param(ar, "alpha", np.nan)),
        "frame": frame_index,
        "step": int(ar.nstart + frame_index * ar.ninfo),
        "time": simulation_time(ar, frame_index),
        "M_phi": float(np.sum(phi)),
        "M_phi_norm": float(np.sum(phi) / first_mass) if first_mass > 0.0 else np.nan,
        "M_m": float(np.sum(fields["m"])),
        "chi_bar": chi_bar,
        "Var_chi_phi": var_chi,
        "H_chi": h_chi,
        "C_phi_chi": corr_phi_chi,
        "high_chi_area_fraction": (
            float(np.sum(high) / np.sum(material)) if np.any(material) else np.nan
        ),
        "low_chi_area_fraction": (
            float(np.sum(low) / np.sum(material)) if np.any(material) else np.nan
        ),
        "mean_phi_high_chi": float(np.mean(phi[high])) if np.any(high) else np.nan,
        "mean_phi_low_chi": float(np.mean(phi[low])) if np.any(low) else np.nan,
        "max_phi": float(np.max(phi)),
    }


def analyze_case(label, path, stride):
    ar = loadarchive(path)
    first_fields = field_set(ar.read_frame(0))
    first_mass = float(np.sum(first_fields["phi"]))
    rows = [frame_metrics(ar, label, idx, first_mass) for idx in frame_indices(ar, stride)]
    return {
        "label": label,
        "path": path,
        "ar": ar,
        "rows": rows,
        "alpha": float(param(ar, "alpha", np.nan)),
        "growTogether": int(param(ar, "growTogether", -1)),
    }


def values(rows, key):
    return np.array([row[key] for row in rows], dtype=float)


def case_color(case):
    key = (int(case["growTogether"]), classify_alpha(case["alpha"]))
    return COLORS.get(key, "#374151")


def case_linestyle(case):
    return "-" if int(case["growTogether"]) == 0 else "--"


def style_axis(ax, ylabel, xlabel=None):
    ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(xlabel)


def plot_core(cases, outfile, dpi):
    panels = [
        ("chi_bar", r"$\bar{\chi}=\langle\chi\rangle_{\phi}$"),
        ("H_chi", r"$H_{\chi}=\langle\chi(1-\chi)\rangle_{\phi}$"),
        ("C_phi_chi", r"$C_{\phi,\chi}$"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.8), sharex=True)
    for ax, (key, ylabel) in zip(axes, panels):
        for case in cases:
            rows = case["rows"]
            ax.plot(
                values(rows, "M_phi_norm"),
                values(rows, key),
                color=case_color(case),
                linestyle=case_linestyle(case),
                linewidth=2.0,
                label=display_label(case),
            )
        style_axis(ax, ylabel, r"$M_{\phi}/M_{\phi}(0)$")
    axes[0].legend(frameon=False, fontsize=8, loc="best")
    fig.suptitle("Q2 core result: phenotype amplification vs density-phenotype coupling", fontsize=17, y=1.03)
    fig.tight_layout()
    fig.savefig(outfile, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def histogram_data(ar, frame_index, nbins):
    fields = field_set(ar.read_frame(frame_index))
    phi = fields["phi"]
    chi = fields["chi"]
    material = phi > PHI_THRESHOLD
    weights = np.maximum(phi[material], 0.0)
    hist, edges = np.histogram(
        chi[material],
        bins=np.linspace(0.0, 1.0, nbins + 1),
        weights=weights,
        density=False,
    )
    total = float(np.sum(hist))
    if total > 0.0:
        hist = hist / total
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, hist


def plot_chi_distributions(cases, outfile, dpi, nbins):
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.2), sharex=True, sharey=True)
    axes = axes.flat
    for ax, case in zip(axes, cases):
        ar = case["ar"]
        n_available = available_frame_count(ar)
        for frame_index, label, linestyle in [
            (0, "initial", "-"),
            (n_available // 2, "mid", "--"),
            (n_available - 1, "final", ":"),
        ]:
            centers, hist = histogram_data(ar, frame_index, nbins)
            ax.plot(centers, hist, linestyle=linestyle, linewidth=2.0, label=label)
        ax.set_title(display_label(case), fontsize=12)
        ax.set_xlabel(r"$\chi$")
        ax.set_ylabel(r"$P_{\phi}(\chi)$")
        ax.legend(frameon=False, fontsize=8, loc="best")
    fig.suptitle(r"Q2 phenotype distribution inside $\phi>0.5$", fontsize=17, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)


def draw_boundary(ax, phi, linewidth=0.9):
    if np.nanmin(phi) <= PHI_THRESHOLD <= np.nanmax(phi):
        ax.contour(
            phi.T,
            levels=[PHI_THRESHOLD],
            colors=[BOUNDARY_COLOR],
            linewidths=linewidth,
            origin="lower",
        )


def chi_phi_rgba(fields):
    phi_alpha = np.clip(
        (fields["phi"] - PHI_VISIBLE_THRESHOLD) / (1.0 - PHI_VISIBLE_THRESHOLD),
        0.0,
        1.0,
    )
    rgb = plt.get_cmap("viridis")(np.clip(fields["chi"], 0.0, 1.0))[..., :3]
    return np.dstack((rgb, phi_alpha))


def hide_axes(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal")


def plot_composite_snapshots(cases, outfile, dpi):
    fig, axes = plt.subplots(2, len(cases), figsize=(3.1 * len(cases), 6.4))
    for col, case in enumerate(cases):
        ar = case["ar"]
        n_available = available_frame_count(ar)
        for row, frame_index in enumerate([0, n_available - 1]):
            fields = field_set(ar.read_frame(frame_index))
            phi = fields["phi"]
            ax = axes[row, col]
            ax.set_facecolor(BACKGROUND_COLOR)
            ax.imshow(
                np.transpose(chi_phi_rgba(fields), (1, 0, 2)),
                origin="lower",
                interpolation="nearest",
            )
            draw_boundary(ax, phi)
            hide_axes(ax)
            if row == 0:
                ax.set_title(display_label(case), fontsize=11)
            if col == 0:
                ax.set_ylabel("initial" if row == 0 else "final", fontsize=13)
    fig.suptitle(r"Q2 $\chi$-$\phi$ composites", fontsize=17, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)


def plot_density_chi_coupling(cases, outfile, dpi, nbins):
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.2), sharex=True, sharey=True)
    for ax, case in zip(axes.flat, cases):
        ar = case["ar"]
        n_available = available_frame_count(ar)
        for frame_index, label, linestyle in [
            (0, "initial", "-"),
            (n_available - 1, "final", "--"),
        ]:
            fields = field_set(ar.read_frame(frame_index))
            phi = fields["phi"]
            chi = fields["chi"]
            material = phi > PHI_THRESHOLD
            centers, mean_phi = chi_bin_phi_curve(phi, chi, material, nbins)
            ax.plot(centers, mean_phi, linestyle=linestyle, linewidth=2.0, label=label)
        ax.set_title(display_label(case), fontsize=12)
        ax.set_xlabel(r"$\chi$ bin")
        ax.set_ylabel(r"$\langle\phi\rangle$ in bin")
        ax.legend(frameon=False, fontsize=8, loc="best")
    fig.suptitle("Q2 density-phenotype coupling", fontsize=17, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)


def write_csv(rows, outfile):
    if not rows:
        return
    fieldnames = sorted({key for row in rows for key in row.keys()})
    with open(outfile, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    archive_paths = discover_archives(args.input_root)
    if not archive_paths:
        raise SystemExit(f"No simulation archives found under {args.input_root}")

    cases = [analyze_case(label, path, args.frame_stride) for label, path in archive_paths]
    cases.sort(key=lambda case: (case["growTogether"], case["alpha"], case["label"]))

    all_rows = []
    for case in cases:
        all_rows.extend(case["rows"])

    write_csv(all_rows, os.path.join(args.outdir, "q2_metrics_timeseries.csv"))
    plot_core(cases, os.path.join(args.outdir, "q2_core_phenotype_amplification.png"), args.dpi)
    plot_chi_distributions(cases, os.path.join(args.outdir, "q2_chi_distribution.png"), args.dpi, args.bins)
    plot_composite_snapshots(cases, os.path.join(args.outdir, "q2_chi_phi_snapshots.png"), args.dpi)
    plot_density_chi_coupling(cases, os.path.join(args.outdir, "q2_density_chi_coupling.png"), args.dpi, args.bins)

    print(f"Loaded {len(cases)} Q2 archives from {args.input_root}")
    for case in cases:
        rows = case["rows"]
        print(
            f"{case['label']}: growTogether={case['growTogether']}, "
            f"alpha={case['alpha']:g}, frames={len(rows)}, "
            f"final chi_bar={rows[-1]['chi_bar']:.4g}, "
            f"final H_chi={rows[-1]['H_chi']:.4g}, "
            f"final C_phi_chi={rows[-1]['C_phi_chi']:.4g}"
        )
    print(f"Saved Q2 analysis outputs in: {args.outdir}")


if __name__ == "__main__":
    main()
