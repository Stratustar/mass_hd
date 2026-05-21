import matplotlib
matplotlib.use("Agg")

import argparse
import csv
import os
import sys
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager


SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from archive.archive import loadarchive


Q1_COLORS = {
    "weak": "#255c99",
    "medium": "#047857",
    "strong": "#b45309",
}
Q2_COLORS = {
    "B1": "#255c99",
    "B2": "#b45309",
    "B3": "#047857",
    "B4": "#7c3aed",
}
Q2_STYLES = {
    "B1": "-",
    "B2": "-",
    "B3": "--",
    "B4": "--",
}
PHI_THRESHOLD = 0.5
PHI_VISIBLE_THRESHOLD = 0.05


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
            "axes.grid": False,
            "axes.titlepad": 6,
        }
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Redraw publication-style Q1/Q2 core figures from saved summaries."
    )
    parser.add_argument("--q1-csv", required=True)
    parser.add_argument("--q1-outdir", required=True)
    parser.add_argument("--q2-csv", required=True)
    parser.add_argument("--q2-archives", required=True)
    parser.add_argument("--q2-outdir", required=True)
    parser.add_argument("--dpi", type=int, default=220)
    return parser.parse_args()


def read_rows(path):
    rows = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            converted = {}
            for key, value in row.items():
                try:
                    converted[key] = float(value)
                except ValueError:
                    converted[key] = value
            rows.append(converted)
    return rows


def order_q1_rows(rows):
    return sorted(rows, key=lambda row: row["Mnorm"])


def group_by(rows, key):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row[key]].append(row)
    return grouped


def style_axis(ax, ylabel, xlabel=None):
    ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(xlabel)


def plot_q1_core(csv_path, outdir, dpi):
    rows = read_rows(csv_path)
    grouped = group_by(rows, "case")
    order = ["weak", "medium", "strong"]
    labels = {
        "weak": r"weak  $\alpha=0.00025$",
        "medium": r"medium  $\alpha=0.001$",
        "strong": r"strong  $\alpha=0.004$",
    }

    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.2))
    for case in order:
        case_rows = order_q1_rows(grouped[case])
        x = np.array([row["Mnorm"] for row in case_rows], dtype=float)
        p = np.array([row["p_rms"] for row in case_rows], dtype=float)
        max_phi = np.array([row["max_phi"] for row in case_rows], dtype=float)
        axes[0].plot(
            x,
            p,
            color=Q1_COLORS[case],
            marker="o",
            linewidth=2.8,
            markersize=7.0,
            label=labels[case],
        )
        axes[1].plot(
            x,
            max_phi,
            color=Q1_COLORS[case],
            marker="o",
            linewidth=2.8,
            markersize=7.0,
            label=labels[case],
        )

    for ax in axes:
        ax.axvline(1.2, color="#9ca3af", linestyle="--", linewidth=1.2)
        ax.spines[["top", "right"]].set_visible(True)

    axes[0].annotate(
        "matched\ncomparison",
        xy=(1.2, grouped["strong"][2]["p_rms"]),
        xytext=(1.34, 0.0066),
        arrowprops={"arrowstyle": "->", "linewidth": 1.2, "color": "black"},
        fontsize=11,
    )
    strong_final = max(grouped["strong"], key=lambda row: row["Mnorm"])
    axes[1].annotate(
        "final strong-alpha\nbuildup",
        xy=(strong_final["Mnorm"], strong_final["max_phi"]),
        xytext=(1.76, 1.34),
        arrowprops={"arrowstyle": "->", "linewidth": 1.2, "color": "black"},
        fontsize=11,
    )

    axes[0].set_title("Weak growth stays low-pressure", fontsize=17)
    axes[1].set_title("Strong growth accumulates density", fontsize=17)
    style_axis(axes[0], r"pressure heterogeneity  $p_{rms}$", r"growth fraction  $M_\phi/M_\phi(0)$")
    style_axis(axes[1], r"density buildup  $\max(\phi)$", r"growth fraction  $M_\phi/M_\phi(0)$")
    axes[0].legend(frameon=False, fontsize=11, loc="upper left")
    fig.suptitle(
        "Q1 core result: weak alpha is quasi-static; strong alpha builds pressure and density",
        fontsize=18,
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "q1_core_pressure_density.png"), dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def q2_case_order(case):
    return {"B1": 0, "B2": 1, "B3": 2, "B4": 3}.get(case, 99)


def plot_q2_core(csv_path, outdir, dpi):
    rows = read_rows(csv_path)
    grouped = group_by(rows, "case")
    panels = [
        ("chi_bar", r"grow fraction  $\bar{\chi}$"),
        ("H_chi", r"mixing  $H_\chi=\langle\chi(1-\chi)\rangle_\phi$"),
        ("C_phi_chi", r"density-phenotype corr.  $C_{\phi\chi}$"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 3.85))
    for case in sorted(grouped, key=q2_case_order):
        case_rows = sorted(grouped[case], key=lambda row: row["Mnorm"])
        x = np.array([row["Mnorm"] for row in case_rows], dtype=float)
        label = str(case_rows[0]["label"])
        for ax, (key, _) in zip(axes, panels):
            y = np.array([row[key] for row in case_rows], dtype=float)
            ax.plot(
                x,
                y,
                color=Q2_COLORS[case],
                linestyle=Q2_STYLES[case],
                marker="o",
                linewidth=2.3,
                markersize=5.0,
                label=label,
            )

    for ax, (_, ylabel) in zip(axes, panels):
        style_axis(ax, ylabel, r"$M_\phi/M_\phi(0)$")
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")
    fig.suptitle("Q2 core result: growth rule controls phenotype amplification", fontsize=18, y=1.04)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "q2_core_phenotype_amplification.png"), dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def archive_marker(files):
    return "parameters.json" in files or "parameters.json.zip" in files


def discover_archives(input_root):
    archives = []
    for root, _, files in os.walk(input_root):
        if archive_marker(files):
            label = os.path.basename(root)
            archives.append((label, root))
    return sorted(archives, key=lambda item: q2_case_order(item[0].split("_", 1)[0]))


def available_frame_count(ar):
    return int((ar.nsteps - ar.nstart) / ar.ninfo) + 1


def field_grid(frame, name):
    lx = frame.parameters["LX"]
    ly = frame.parameters["LY"]
    return np.asarray(getattr(frame, name), dtype=float).reshape((lx, ly))


def field_set(frame):
    phi = field_grid(frame, "phi")
    if hasattr(frame, "chi"):
        chi = field_grid(frame, "chi")
    else:
        m = field_grid(frame, "m")
        chi = np.divide(m, phi, out=np.zeros_like(phi), where=phi > PHI_VISIBLE_THRESHOLD)
    return phi, np.clip(chi, 0.0, 1.0)


def histogram_data(ar, frame_index, nbins):
    phi, chi = field_set(ar.read_frame(frame_index))
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


def q2_display_label(label):
    case = label.split("_", 1)[0]
    return {
        "B1": "gT=0 weak",
        "B2": "gT=0 strong",
        "B3": "gT=1 weak",
        "B4": "gT=1 strong",
    }.get(case, label)


def plot_q2_distribution(input_root, outdir, dpi, nbins=50):
    archives = discover_archives(input_root)
    fig, axes = plt.subplots(2, 2, figsize=(10.2, 7.35), sharex=True, sharey=True)
    for ax, (label, path) in zip(axes.flat, archives):
        ar = loadarchive(path)
        final_frame = available_frame_count(ar) - 1
        for frame_index, line_label, linestyle, color in [
            (0, "initial", "-", "#1f77b4"),
            (final_frame, "final", "--", "#ff7f0e"),
        ]:
            centers, hist = histogram_data(ar, frame_index, nbins)
            ax.plot(centers, hist, linestyle=linestyle, color=color, linewidth=2.2, label=line_label)
        ax.set_title(q2_display_label(label), fontsize=13)
        ax.set_xlabel(r"$\chi$")
        ax.set_ylabel(r"$P_\phi(\chi)$")
        ax.legend(frameon=False, fontsize=9, loc="best")
    fig.suptitle("Q2 phenotype distribution: initial vs final", fontsize=18, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(os.path.join(outdir, "q2_chi_distribution_initial_final.png"), dpi=dpi)
    plt.close(fig)


def main():
    configure_plot_style()
    args = parse_args()
    os.makedirs(args.q1_outdir, exist_ok=True)
    os.makedirs(args.q2_outdir, exist_ok=True)
    plot_q1_core(args.q1_csv, args.q1_outdir, args.dpi)
    plot_q2_core(args.q2_csv, args.q2_outdir, args.dpi)
    plot_q2_distribution(args.q2_archives, args.q2_outdir, args.dpi)


if __name__ == "__main__":
    main()
