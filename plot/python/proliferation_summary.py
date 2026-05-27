import matplotlib
matplotlib.use("Agg")

import argparse
import csv
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update(
    {
        "mathtext.fontset": "stix",
        "font.family": "STIXGeneral",
    }
)


SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from archive.archive import loadarchive


PHI_THRESHOLD = 0.5
COMPOSITE_BACKGROUND = "#7e7e7e"
BOUNDARY_COLOR = "#d9468f"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize passive proliferation test output."
    )
    parser.add_argument("inputdir", help="Simulation archive directory")
    parser.add_argument("outdir", help="Directory for plots and CSV output")
    parser.add_argument(
        "--dpi",
        type=int,
        default=250,
        help="Figure resolution. Default: 250",
    )
    parser.add_argument(
        "--snapshot-count",
        type=int,
        default=5,
        help="Number of composite snapshots. Default: 5",
    )
    return parser.parse_args()


def available_frame_count(ar):
    return int((ar.nsteps - ar.nstart) / ar.ninfo) + 1


def simulation_time(ar, frame_index):
    step = ar.nstart + frame_index * ar.ninfo
    dt = float(getattr(ar, "time_step", 1.0))
    return step * dt


def grid(frame, name):
    lx = frame.parameters["LX"]
    ly = frame.parameters["LY"]
    return np.asarray(getattr(frame, name), dtype=float).reshape((lx, ly))


def velocity_fields(frame):
    lx = frame.parameters["LX"]
    ly = frame.parameters["LY"]
    ff = np.asarray(frame.ff, dtype=float)
    density = np.sum(ff, axis=1)
    density = np.where(density == 0.0, np.nan, density)
    vx = (ff[:, 1] - ff[:, 2] + ff[:, 5] - ff[:, 6] - ff[:, 7] + ff[:, 8]) / density
    vy = (ff[:, 3] - ff[:, 4] + ff[:, 5] - ff[:, 6] + ff[:, 7] - ff[:, 8]) / density
    vx = np.nan_to_num(vx).reshape((lx, ly))
    vy = np.nan_to_num(vy).reshape((lx, ly))
    return vx, vy


def vorticity(vx, vy):
    dvdx = 0.5 * (np.roll(vy, -1, axis=0) - np.roll(vy, 1, axis=0))
    dudy = 0.5 * (np.roll(vx, -1, axis=1) - np.roll(vx, 1, axis=1))
    return dvdx - dudy


def perimeter_4(mask):
    padded = np.pad(mask.astype(bool), 1, mode="constant", constant_values=False)
    center = padded[1:-1, 1:-1]
    neighbors = (
        padded[2:, 1:-1],
        padded[:-2, 1:-1],
        padded[1:-1, 2:],
        padded[1:-1, :-2],
    )
    exposed = sum(center & ~neighbor for neighbor in neighbors)
    return float(np.sum(exposed))


def boundary_roughness(mask):
    points = np.argwhere(mask)
    if len(points) == 0:
        return np.nan

    padded = np.pad(mask.astype(bool), 1, mode="constant", constant_values=False)
    center = padded[1:-1, 1:-1]
    interior = (
        center
        & padded[2:, 1:-1]
        & padded[:-2, 1:-1]
        & padded[1:-1, 2:]
        & padded[1:-1, :-2]
    )
    boundary = center & ~interior
    boundary_points = np.argwhere(boundary)
    if len(boundary_points) == 0:
        return np.nan

    weights = mask.astype(float)
    yy, xx = np.indices(mask.shape)
    total = np.sum(weights)
    if total <= 0:
        center_x = 0.5 * (mask.shape[0] - 1)
        center_y = 0.5 * (mask.shape[1] - 1)
    else:
        center_x = np.sum(xx * weights) / total
        center_y = np.sum(yy * weights) / total

    bx = boundary_points[:, 1]
    by = boundary_points[:, 0]
    radii = np.sqrt((bx - center_x) ** 2 + (by - center_y) ** 2)
    return float(np.std(radii))


def frame_metrics(ar, frame_index):
    frame = ar.read_frame(frame_index)
    phi = grid(frame, "phi")
    chi = grid(frame, "chi") if hasattr(frame, "chi") else np.zeros_like(phi)
    m = grid(frame, "m") if hasattr(frame, "m") else phi * chi
    vx, vy = velocity_fields(frame)
    speed2 = vx * vx + vy * vy
    omega = vorticity(vx, vy)

    mass_phi = float(np.sum(phi))
    mass_m = float(np.sum(m))
    phi_weight = np.maximum(phi, 0.0)
    phi_weight_sum = float(np.sum(phi_weight))
    mask = phi > PHI_THRESHOLD
    area = float(np.sum(mask))
    perimeter = perimeter_4(mask)
    circularity = (
        float(4.0 * np.pi * area / (perimeter * perimeter))
        if perimeter > 0.0
        else np.nan
    )

    return {
        "frame": frame_index,
        "step": int(ar.nstart + frame_index * ar.ninfo),
        "time": simulation_time(ar, frame_index),
        "M_phi": mass_phi,
        "M_m": mass_m,
        "chi_bar": mass_m / mass_phi if mass_phi > 0.0 else np.nan,
        "area_phi_05": area,
        "perimeter_phi_05": perimeter,
        "R_eff": np.sqrt(area / np.pi) if area > 0.0 else np.nan,
        "boundary_roughness": boundary_roughness(mask),
        "circularity": circularity,
        "max_phi": float(np.max(phi)),
        "mean_phi_material": float(np.mean(phi[mask])) if np.any(mask) else np.nan,
        "u_rms": (
            float(np.sqrt(np.sum(phi_weight * speed2) / phi_weight_sum))
            if phi_weight_sum > 0.0
            else np.nan
        ),
        "max_speed": float(np.sqrt(np.max(speed2))),
        "enstrophy": (
            float(0.5 * np.sum(phi_weight * omega * omega) / phi_weight_sum)
            if phi_weight_sum > 0.0
            else np.nan
        ),
    }


def add_effective_growth(rows):
    times = np.array([row["time"] for row in rows], dtype=float)
    masses = np.array([row["M_phi"] for row in rows], dtype=float)
    if len(rows) < 2 or np.any(masses <= 0.0):
        for row in rows:
            row["g_eff_mass"] = np.nan
            row["tau_grow_mass"] = np.nan
        return rows

    growth = np.gradient(np.log(masses), times)
    for row, g in zip(rows, growth):
        row["g_eff_mass"] = float(g)
        row["tau_grow_mass"] = float(1.0 / g) if g > 0.0 else np.nan
    return rows


def write_csv(rows, outfile):
    fieldnames = list(rows[0].keys())
    with open(outfile, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def plot_summary(rows, outfile, dpi):
    t = np.array([row["time"] for row in rows], dtype=float)

    panels = [
        ("M_phi", r"$M_{\phi}$"),
        ("area_phi_05", r"$A_{\phi>0.5}$"),
        ("R_eff", r"$R_{\mathrm{eff}}$"),
        ("u_rms", r"$\langle u^2\rangle_{\phi}^{1/2}$"),
        ("enstrophy", r"$\frac{1}{2}\langle\omega^2\rangle_{\phi}$"),
        ("max_phi", r"$\max(\phi)$"),
        ("chi_bar", r"$\bar{\chi}$"),
        ("g_eff_mass", r"$d\log(M_{\phi})/dt$"),
    ]

    fig, axes = plt.subplots(4, 2, figsize=(9.0, 10.5), sharex=True)
    for ax, (key, ylabel) in zip(axes.flat, panels):
        values = np.array([row[key] for row in rows], dtype=float)
        ax.plot(t, values, color="#255c99", linewidth=1.8)
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.25, linewidth=0.6)

    axes[-1, 0].set_xlabel(r"$t$")
    axes[-1, 1].set_xlabel(r"$t$")
    fig.suptitle(
        r"Passive proliferation baseline, $\zeta=0$, $\chi_0=1$",
        fontsize=18,
        y=0.995,
    )
    fig.tight_layout()
    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)


def composite_rgba(phi, chi, gamma=0.8):
    chi = np.clip(chi, 0.0, 1.0)
    phi_alpha = np.clip((phi - 0.05) / (1.0 - 0.05), 0.0, 1.0) ** gamma
    rgb = plt.get_cmap("viridis")(chi)[..., :3]
    return np.dstack((rgb, phi_alpha))


def snapshot_indices(n_available, snapshot_count):
    if snapshot_count <= 1:
        return [0]
    return sorted(set(np.linspace(0, n_available - 1, snapshot_count, dtype=int).tolist()))


def plot_composite_snapshots(ar, outfile, dpi, snapshot_count):
    n_available = available_frame_count(ar)
    indices = snapshot_indices(n_available, snapshot_count)
    fig, axes = plt.subplots(1, len(indices), figsize=(3.1 * len(indices), 3.25))
    if len(indices) == 1:
        axes = [axes]

    for ax, frame_index in zip(axes, indices):
        frame = ar.read_frame(frame_index)
        phi = grid(frame, "phi")
        chi = grid(frame, "chi") if hasattr(frame, "chi") else np.zeros_like(phi)
        ax.set_facecolor(COMPOSITE_BACKGROUND)
        ax.imshow(
            np.transpose(composite_rgba(phi, chi), (1, 0, 2)),
            origin="lower",
            interpolation="nearest",
        )
        ax.contour(
            phi.T,
            levels=[PHI_THRESHOLD],
            colors=[BOUNDARY_COLOR],
            linewidths=0.8,
            origin="lower",
        )
        ax.set_title(rf"$t={simulation_time(ar, frame_index):.0f}$", fontsize=14)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")

    fig.suptitle(r"$\chi$-$\phi$ composite snapshots", fontsize=18, y=0.98)
    fig.tight_layout()
    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    ar = loadarchive(args.inputdir)
    n_available = available_frame_count(ar)
    rows = [frame_metrics(ar, frame_index) for frame_index in range(n_available)]
    rows = add_effective_growth(rows)

    csv_path = os.path.join(args.outdir, "proliferation_timeseries.csv")
    summary_path = os.path.join(args.outdir, "proliferation_summary.png")
    snapshots_path = os.path.join(args.outdir, "chi_phi_snapshots.png")

    write_csv(rows, csv_path)
    plot_summary(rows, summary_path, args.dpi)
    plot_composite_snapshots(ar, snapshots_path, args.dpi, args.snapshot_count)

    print(f"Loaded {n_available} frames from {args.inputdir}")
    print(f"Saved CSV: {csv_path}")
    print(f"Saved summary figure: {summary_path}")
    print(f"Saved composite snapshots: {snapshots_path}")


if __name__ == "__main__":
    main()
