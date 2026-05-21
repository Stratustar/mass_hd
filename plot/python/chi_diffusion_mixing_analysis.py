import matplotlib
matplotlib.use("Agg")

import argparse
import csv
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager


SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from archive.archive import loadarchive


COLORS = {
    0.0: "#255c99",
    0.005: "#047857",
    0.02: "#b45309",
}
PHI_THRESHOLD = 0.5


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
            "mathtext.cal": "Arial",
            "mathtext.sf": "Arial",
            "mathtext.tt": "Arial",
            "axes.grid": False,
            "axes.titlepad": 6,
        }
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyze chi diffusion-only mixing in the validation sweep."
    )
    parser.add_argument("input_root", help="Directory containing D*_diffusion_* archives.")
    parser.add_argument("outdir", help="Directory for CSV and figure outputs.")
    parser.add_argument("--frame-stride", type=int, default=1)
    parser.add_argument("--phi-threshold", type=float, default=PHI_THRESHOLD)
    parser.add_argument("--dpi", type=int, default=220)
    return parser.parse_args()


def archive_marker(files):
    return "parameters.json" in files or "parameters.json.zip" in files


def discover_diffusion_archives(input_root):
    archives = []
    for root, _, files in os.walk(input_root):
        if not archive_marker(files):
            continue
        label = os.path.basename(root)
        if not label.startswith("D"):
            continue
        ar = loadarchive(root)
        if abs(float(getattr(ar, "zeta", np.nan))) > 1e-12:
            continue
        if abs(float(getattr(ar, "alpha", np.nan))) > 1e-12:
            continue
        if float(getattr(ar, "Dchi", np.nan)) < -1e-12:
            continue
        archives.append((label, root, float(getattr(ar, "Dchi", np.nan))))
    return sorted(archives, key=lambda item: item[2])


def available_frame_count(ar):
    return int((ar.nsteps - ar.nstart) / ar.ninfo) + 1


def simulation_time(ar, frame_index):
    step = ar.nstart + frame_index * ar.ninfo
    return step * float(getattr(ar, "time_step", 1.0))


def weighted_mean(values, weights):
    total = float(np.sum(weights))
    if total <= 0.0:
        return np.nan
    return float(np.sum(values * weights) / total)


def analyze_archive(label, path, phi_threshold, frame_stride):
    ar = loadarchive(path)
    frame_indices = list(range(0, available_frame_count(ar), frame_stride))
    if frame_indices[-1] != available_frame_count(ar) - 1:
        frame_indices.append(available_frame_count(ar) - 1)

    first_mass = None
    rows = []
    for frame_index in frame_indices:
        frame = ar.read_frame(frame_index)
        phi = np.asarray(frame.phi, dtype=float)
        chi = np.asarray(frame.chi, dtype=float)
        material = phi > phi_threshold
        weights = np.where(material, np.maximum(phi, 0.0), 0.0)
        mass = float(np.sum(phi))
        if first_mass is None:
            first_mass = mass

        h_chi = weighted_mean(chi * (1.0 - chi), weights)
        chi_bar = weighted_mean(chi, weights)
        rows.append(
            {
                "case": label,
                "Dchi": float(getattr(ar, "Dchi", np.nan)),
                "frame": frame_index,
                "step": int(ar.nstart + frame_index * ar.ninfo),
                "time": simulation_time(ar, frame_index),
                "M_phi": mass,
                "M_phi_norm": mass / first_mass if first_mass and first_mass > 0.0 else np.nan,
                "H_chi": h_chi,
                "chi_bar": chi_bar,
            }
        )
    return rows


def write_csv(rows, path):
    if not rows:
        return
    keys = list(rows[0].keys())
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def fmt_dchi(value):
    if abs(value) < 1e-14:
        return "0"
    return f"{value:g}"


def plot_mixing(rows, outdir, dpi):
    grouped = {}
    for row in rows:
        grouped.setdefault(row["case"], []).append(row)

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.0))
    for case, case_rows in sorted(grouped.items(), key=lambda item: item[1][0]["Dchi"]):
        case_rows = sorted(case_rows, key=lambda row: row["frame"])
        dchi = float(case_rows[0]["Dchi"])
        color = COLORS.get(dchi, None)
        label = f"Dchi={fmt_dchi(dchi)}"
        mnorm = np.array([row["M_phi_norm"] for row in case_rows], dtype=float)
        time = np.array([row["time"] for row in case_rows], dtype=float)
        mixing = np.array([row["H_chi"] for row in case_rows], dtype=float)
        axes[0].plot(
            mnorm,
            mixing,
            color=color,
            linewidth=2.4,
            marker="o",
            markevery=max(1, len(mixing) // 7),
            markersize=4.5,
            label=label,
        )
        axes[1].plot(
            time,
            mixing,
            color=color,
            linewidth=2.4,
            marker="o",
            markevery=max(1, len(mixing) // 7),
            markersize=4.5,
            label=label,
        )

    all_mnorm = np.array([row["M_phi_norm"] for row in rows], dtype=float)
    if np.nanmax(all_mnorm) - np.nanmin(all_mnorm) < 1e-3:
        center = float(np.nanmean(all_mnorm))
        axes[0].set_xlim(center - 0.004, center + 0.004)

    axes[0].set_title("Mixing at fixed mass", fontsize=16)
    axes[0].set_xlabel("M/M0")
    axes[0].set_ylabel("mixing  Hchi = <chi(1-chi)>_phi")
    axes[0].legend(frameon=False, fontsize=10, loc="upper left")

    axes[1].set_title("Diffusion increases mixing", fontsize=16)
    axes[1].set_xlabel(r"time")
    axes[1].set_ylabel("mixing  Hchi")
    axes[1].legend(frameon=False, fontsize=10, loc="upper left")

    fig.suptitle("Diffusion-only chi dynamics: no mass growth, increasing phenotype mixing", fontsize=18, y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "chi_diffusion_mixing_vs_mass.png"), dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    configure_plot_style()
    os.makedirs(args.outdir, exist_ok=True)

    archives = discover_diffusion_archives(args.input_root)
    if not archives:
        raise SystemExit(f"No diffusion-only archives found under {args.input_root}")

    all_rows = []
    for label, path, _ in archives:
        print(f"Analyzing {label}...")
        all_rows.extend(analyze_archive(label, path, args.phi_threshold, args.frame_stride))

    write_csv(all_rows, os.path.join(args.outdir, "chi_diffusion_mixing_metrics.csv"))
    plot_mixing(all_rows, args.outdir, args.dpi)

    for case in sorted({row["case"] for row in all_rows}):
        case_rows = [row for row in all_rows if row["case"] == case]
        case_rows.sort(key=lambda row: row["frame"])
        m_values = np.array([row["M_phi_norm"] for row in case_rows], dtype=float)
        print(
            f"{case}: Dchi={case_rows[0]['Dchi']:g}, "
            f"M/M0 range={np.nanmin(m_values):.8f}-{np.nanmax(m_values):.8f}, "
            f"H_chi {case_rows[0]['H_chi']:.5f}->{case_rows[-1]['H_chi']:.5f}"
        )
    print(f"Saved diffusion mixing outputs in: {args.outdir}")


if __name__ == "__main__":
    main()
