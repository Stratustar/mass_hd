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
            "mathtext.cal": "Arial",
            "mathtext.sf": "Arial",
            "mathtext.tt": "Arial",
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
EDGE_LOW = 0.4
EDGE_HIGH = 0.6
BULK_THRESHOLD = 0.8
BOUNDARY_COLOR = "#d9468f"
VISUALIZATION_BACKGROUND = "#7e7e7e"
LINE_COLORS = ["#255c99", "#b45309", "#047857", "#7c3aed", "#be123c", "#374151"]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Q1 passive proliferation analysis: compare alpha cases at matched "
            "growth fractions."
        )
    )
    parser.add_argument(
        "input_root",
        help="Directory containing one or more simulation archive directories.",
    )
    parser.add_argument("outdir", help="Directory for CSV and PNG outputs.")
    parser.add_argument(
        "--growth-targets",
        type=float,
        nargs="+",
        default=[1.0, 1.1, 1.2],
        help="M_phi/M_phi(0) targets used for matched comparisons.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=220,
        help="PNG resolution. Default: 220.",
    )
    parser.add_argument(
        "--smooth-window",
        type=int,
        default=5,
        help="Odd moving-average window for g_eff. Default: 5.",
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


def simulation_time(ar, frame_index):
    step = ar.nstart + frame_index * ar.ninfo
    return step * float(getattr(ar, "time_step", 1.0))


def param(ar, name, default=np.nan):
    return ar.parameters.get(name, default)


def format_alpha(alpha):
    if not np.isfinite(alpha):
        return "nan"
    if alpha == 0:
        return "0"
    return f"{alpha:.3g}"


def grid(frame, name):
    lx = frame.parameters["LX"]
    ly = frame.parameters["LY"]
    return np.asarray(getattr(frame, name), dtype=float).reshape((lx, ly))


def optional_grid(frame, name, fallback):
    if hasattr(frame, name):
        return grid(frame, name)
    return np.array(fallback, dtype=float, copy=True)


def velocity_fields(frame):
    lx = frame.parameters["LX"]
    ly = frame.parameters["LY"]
    ff = np.asarray(frame.ff, dtype=float)
    density = np.sum(ff, axis=1)
    density = np.where(density == 0.0, np.nan, density)
    vx = (ff[:, 1] - ff[:, 2] + ff[:, 5] - ff[:, 6] - ff[:, 7] + ff[:, 8]) / density
    vy = (ff[:, 3] - ff[:, 4] + ff[:, 5] - ff[:, 6] + ff[:, 7] - ff[:, 8]) / density
    return np.nan_to_num(vx).reshape((lx, ly)), np.nan_to_num(vy).reshape((lx, ly))


def vorticity(vx, vy):
    dvdx = 0.5 * (np.roll(vy, -1, axis=0) - np.roll(vy, 1, axis=0))
    dudy = 0.5 * (np.roll(vx, -1, axis=1) - np.roll(vx, 1, axis=1))
    return dvdx - dudy


def pressure_proxy(frame, phi):
    if hasattr(frame, "pressure"):
        return grid(frame, "pressure"), "pressure"
    if hasattr(frame, "sigma_bulk"):
        return -grid(frame, "sigma_bulk"), "-sigma_bulk"
    return None, "unavailable"


def field_set(frame):
    phi = grid(frame, "phi")
    m = optional_grid(frame, "m", np.zeros_like(phi))
    chi_raw = optional_grid(frame, "chi", np.full_like(phi, np.nan))
    chi_from_m = np.divide(m, phi, out=np.zeros_like(phi), where=phi > PHI_VISIBLE_THRESHOLD)
    chi = chi_raw if np.any(np.isfinite(chi_raw)) else chi_from_m
    chi = np.clip(chi, 0.0, 1.0)

    qxx = optional_grid(frame, "QQxx", np.zeros_like(phi))
    qyx = optional_grid(frame, "QQyx", np.zeros_like(phi))
    qmag = np.sqrt(qxx * qxx + qyx * qyx)
    theta = 0.5 * np.arctan2(qyx, qxx)
    pressure, pressure_source = pressure_proxy(frame, phi)
    vx, vy = velocity_fields(frame)
    speed = np.sqrt(vx * vx + vy * vy)
    omega = vorticity(vx, vy)

    return {
        "phi": phi,
        "m": m,
        "chi": chi,
        "qxx": qxx,
        "qyx": qyx,
        "qmag": qmag,
        "theta": theta,
        "pressure": pressure,
        "pressure_source": pressure_source,
        "vx": vx,
        "vy": vy,
        "speed": speed,
        "omega": omega,
    }


def center_of_mass(field):
    weights = np.maximum(field, 0.0)
    total = float(np.sum(weights))
    lx, ly = field.shape
    if total <= 0.0:
        return 0.5 * (lx - 1), 0.5 * (ly - 1)
    x, y = np.meshgrid(np.arange(lx), np.arange(ly), indexing="ij")
    return float(np.sum(x * weights) / total), float(np.sum(y * weights) / total)


def perimeter_4(mask):
    padded = np.pad(mask.astype(bool), 1, mode="constant", constant_values=False)
    center = padded[1:-1, 1:-1]
    exposed = (
        (center & ~padded[2:, 1:-1]).astype(int)
        + (center & ~padded[:-2, 1:-1]).astype(int)
        + (center & ~padded[1:-1, 2:]).astype(int)
        + (center & ~padded[1:-1, :-2]).astype(int)
    )
    return float(np.sum(exposed))


def boundary_roughness(mask, center):
    padded = np.pad(mask.astype(bool), 1, mode="constant", constant_values=False)
    active = padded[1:-1, 1:-1]
    interior = (
        active
        & padded[2:, 1:-1]
        & padded[:-2, 1:-1]
        & padded[1:-1, 2:]
        & padded[1:-1, :-2]
    )
    boundary = active & ~interior
    points = np.argwhere(boundary)
    if len(points) == 0:
        return np.nan
    radii = np.sqrt((points[:, 0] - center[0]) ** 2 + (points[:, 1] - center[1]) ** 2)
    return float(np.std(radii))


def weighted_mean(field, weights):
    total = float(np.sum(weights))
    if total <= 0.0:
        return np.nan
    return float(np.sum(field * weights) / total)


def moving_average(values, window):
    values = np.asarray(values, dtype=float)
    if window <= 1 or len(values) == 0:
        return values
    window = int(window)
    if window % 2 == 0:
        window += 1
    if len(values) < window:
        return values
    kernel = np.ones(window, dtype=float) / window
    padded = np.pad(values, window // 2, mode="edge")
    return np.convolve(padded, kernel, mode="valid")


def finite_gradient(values, times):
    values = np.asarray(values, dtype=float)
    times = np.asarray(times, dtype=float)
    if len(values) < 2:
        return np.full_like(values, np.nan)
    return np.gradient(values, times)


def frame_metrics(ar, label, frame_index):
    frame = ar.read_frame(frame_index)
    fields = field_set(frame)
    phi = fields["phi"]
    qxx = fields["qxx"]
    qyx = fields["qyx"]
    qmag = fields["qmag"]
    pressure = fields["pressure"]
    speed = fields["speed"]
    omega = fields["omega"]

    material = phi > PHI_THRESHOLD
    edge = (phi > EDGE_LOW) & (phi < EDGE_HIGH)
    bulk = phi > BULK_THRESHOLD
    material_weight = np.where(material, np.maximum(phi, 0.0), 0.0)
    material_weight_sum = float(np.sum(material_weight))
    area = float(np.sum(material))
    perimeter = perimeter_4(material)
    radius = np.sqrt(area / np.pi) if area > 0.0 else np.nan
    center = center_of_mass(np.where(material, phi, 0.0))
    roughness = boundary_roughness(material, center)

    if pressure is not None and np.any(material):
        pressure_material = pressure[material]
        pressure_prime = pressure_material - np.mean(pressure_material)
        p_mean = float(np.mean(pressure_material))
        p_rms = float(np.sqrt(np.mean(pressure_prime * pressure_prime)))
        p95 = float(np.percentile(pressure_material, 95))
        p5 = float(np.percentile(pressure_material, 5))
    else:
        p_mean = p_rms = p95 = p5 = np.nan

    qxx_mean = weighted_mean(qxx, material_weight)
    qyx_mean = weighted_mean(qyx, material_weight)
    s_global = (
        float(np.sqrt(qxx_mean * qxx_mean + qyx_mean * qyx_mean))
        if np.isfinite(qxx_mean) and np.isfinite(qyx_mean)
        else np.nan
    )

    return {
        "case": label,
        "alpha": float(param(ar, "alpha", np.nan)),
        "frame": frame_index,
        "step": int(ar.nstart + frame_index * ar.ninfo),
        "time": simulation_time(ar, frame_index),
        "M_phi": float(np.sum(phi)),
        "M_m": float(np.sum(fields["m"])),
        "chi_bar": weighted_mean(fields["chi"], material_weight),
        "area_phi05": area,
        "perimeter_phi05": perimeter,
        "R_eff": radius,
        "roughness": roughness,
        "roughness_norm": roughness / radius if radius > 0.0 else np.nan,
        "circularity": (
            float(4.0 * np.pi * area / (perimeter * perimeter))
            if perimeter > 0.0
            else np.nan
        ),
        "mean_phi_material": float(np.mean(phi[material])) if np.any(material) else np.nan,
        "mean_phi_bulk": float(np.mean(phi[bulk])) if np.any(bulk) else np.nan,
        "sigma_phi": float(np.std(phi[material])) if np.any(material) else np.nan,
        "max_phi": float(np.max(phi)),
        "pressure_source": fields["pressure_source"],
        "p_mean": p_mean,
        "p_rms": p_rms,
        "p95": p95,
        "p5": p5,
        "u_rms": (
            float(np.sqrt(np.sum(material_weight * speed * speed) / material_weight_sum))
            if material_weight_sum > 0.0
            else np.nan
        ),
        "u_max": float(np.max(speed[material])) if np.any(material) else np.nan,
        "enstrophy": (
            float(0.5 * np.sum(omega[material] * omega[material]))
            if np.any(material)
            else np.nan
        ),
        "enstrophy_density": (
            float(0.5 * np.mean(omega[material] * omega[material]))
            if np.any(material)
            else np.nan
        ),
        "S_mean": weighted_mean(qmag, material_weight),
        "S_global": s_global,
        "S_edge": float(np.mean(qmag[edge])) if np.any(edge) else np.nan,
    }


def analyze_case(label, path, smooth_window):
    ar = loadarchive(path)
    rows = [frame_metrics(ar, label, i) for i in range(available_frame_count(ar))]

    m0 = rows[0]["M_phi"] if rows and rows[0]["M_phi"] != 0.0 else np.nan
    r0 = rows[0]["R_eff"] if rows and rows[0]["R_eff"] != 0.0 else np.nan
    times = np.array([row["time"] for row in rows], dtype=float)
    masses = np.array([row["M_phi"] for row in rows], dtype=float)
    radii = np.array([row["R_eff"] for row in rows], dtype=float)

    g_eff = (
        moving_average(finite_gradient(np.log(masses), times), smooth_window)
        if len(rows) > 1 and np.all(masses > 0.0)
        else np.full(len(rows), np.nan)
    )
    dRdt = moving_average(finite_gradient(radii, times), smooth_window)

    for row, g, radial_speed in zip(rows, g_eff, dRdt):
        row["M_phi_norm"] = row["M_phi"] / m0 if np.isfinite(m0) else np.nan
        row["R_eff_norm"] = row["R_eff"] / r0 if np.isfinite(r0) else np.nan
        row["g_eff"] = float(g)
        row["tau_grow"] = float(1.0 / g) if g > 0.0 else np.nan
        row["dRdt"] = float(radial_speed)

    alpha = float(param(ar, "alpha", np.nan))
    return {"label": label, "path": path, "ar": ar, "alpha": alpha, "rows": rows}


def values(rows, key):
    return np.array([row[key] for row in rows], dtype=float)


def style_axis(ax, ylabel, xlabel=None):
    ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(xlabel)


def case_display_name(case):
    return f"{case['label']}\n$\\alpha={format_alpha(case['alpha'])}$"


def plot_timeseries(cases, outfile, dpi):
    panels = [
        ("M_phi_norm", r"$M_{\phi}/M_{\phi}(0)$"),
        ("g_eff", r"$d\log M_{\phi}/dt$"),
        ("p_rms", r"$p_{\mathrm{rms}}$"),
        ("sigma_phi", r"$\sigma_{\phi}$"),
        ("u_rms", r"$u_{\mathrm{rms}}$"),
        ("enstrophy_density", r"$\frac{1}{2}\langle\omega^2\rangle$"),
        ("roughness_norm", r"$w_R/R_{\mathrm{eff}}$"),
        ("S_global", r"$S_{\mathrm{global}}$"),
    ]
    fig, axes = plt.subplots(4, 2, figsize=(11.0, 12.0), sharex=True)
    for ax, (key, ylabel) in zip(axes.flat, panels):
        for idx, case in enumerate(cases):
            rows = case["rows"]
            ax.plot(
                values(rows, "time"),
                values(rows, key),
                color=LINE_COLORS[idx % len(LINE_COLORS)],
                linewidth=1.8,
                label=case_display_name(case),
            )
        style_axis(ax, ylabel)
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$t$")
    axes[0, 0].legend(frameon=False, fontsize=8, loc="best")
    fig.suptitle("Q1 passive proliferation: time response", fontsize=20, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)


def plot_growth_response(cases, outfile, dpi):
    panels = [
        ("p_rms", r"$p_{\mathrm{rms}}$"),
        ("p95", r"$p_{95}$"),
        ("sigma_phi", r"$\sigma_{\phi}$"),
        ("u_rms", r"$u_{\mathrm{rms}}$"),
        ("enstrophy_density", r"$\frac{1}{2}\langle\omega^2\rangle$"),
        ("roughness_norm", r"$w_R/R_{\mathrm{eff}}$"),
        ("circularity", r"circularity"),
        ("S_global", r"$S_{\mathrm{global}}$"),
        ("S_edge", r"$S_{\mathrm{edge}}$"),
    ]
    fig, axes = plt.subplots(3, 3, figsize=(12.0, 10.2), sharex=True)
    for ax, (key, ylabel) in zip(axes.flat, panels):
        for idx, case in enumerate(cases):
            rows = case["rows"]
            ax.plot(
                values(rows, "M_phi_norm"),
                values(rows, key),
                color=LINE_COLORS[idx % len(LINE_COLORS)],
                linewidth=1.8,
                label=case_display_name(case),
            )
        style_axis(ax, ylabel)
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$M_{\phi}/M_{\phi}(0)$")
    axes[0, 0].legend(frameon=False, fontsize=8, loc="best")
    fig.suptitle("Q1 passive proliferation: response at matched growth", fontsize=20, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)


def nearest_target_row(case, target):
    rows = case["rows"]
    growth = values(rows, "M_phi_norm")
    valid = np.isfinite(growth)
    if not np.any(valid) or target > np.nanmax(growth[valid]) + 1e-9:
        return None
    valid_indices = np.flatnonzero(valid)
    chosen = valid_indices[int(np.argmin(np.abs(growth[valid] - target)))]
    return rows[int(chosen)]


def build_matched_summary(cases, targets):
    summary = []
    for target in targets:
        for case in cases:
            row = nearest_target_row(case, target)
            if row is None:
                summary.append(
                    {
                        "target_M_phi_norm": target,
                        "case": case["label"],
                        "alpha": case["alpha"],
                        "reached": 0,
                    }
                )
                continue

            out = {"target_M_phi_norm": target, "reached": 1}
            for key in [
                "case",
                "alpha",
                "frame",
                "step",
                "time",
                "M_phi_norm",
                "R_eff_norm",
                "p_rms",
                "p95",
                "p5",
                "sigma_phi",
                "u_rms",
                "enstrophy",
                "enstrophy_density",
                "roughness_norm",
                "circularity",
                "S_mean",
                "S_global",
                "S_edge",
            ]:
                out[key] = row.get(key, np.nan)
            summary.append(out)
    return summary


def plot_matched_metrics(summary, targets, outfile, dpi):
    panels = [
        ("p_rms", r"$p_{\mathrm{rms}}$"),
        ("sigma_phi", r"$\sigma_{\phi}$"),
        ("u_rms", r"$u_{\mathrm{rms}}$"),
        ("enstrophy_density", r"$\frac{1}{2}\langle\omega^2\rangle$"),
        ("roughness_norm", r"$w_R/R_{\mathrm{eff}}$"),
        ("circularity", r"circularity"),
        ("S_global", r"$S_{\mathrm{global}}$"),
        ("S_edge", r"$S_{\mathrm{edge}}$"),
    ]
    alphas = sorted({row["alpha"] for row in summary if row.get("reached") and np.isfinite(row["alpha"])})
    if not alphas:
        return

    alpha_to_x = {alpha: idx for idx, alpha in enumerate(alphas)}
    fig, axes = plt.subplots(4, 2, figsize=(10.5, 11.4), sharex=True)
    markers = ["o", "s", "^", "D", "v"]
    for ax, (key, ylabel) in zip(axes.flat, panels):
        for tidx, target in enumerate(targets):
            rows = [
                row
                for row in summary
                if row.get("reached") and np.isclose(row["target_M_phi_norm"], target)
            ]
            x = [alpha_to_x[row["alpha"]] for row in rows]
            y = [row.get(key, np.nan) for row in rows]
            ax.plot(
                x,
                y,
                marker=markers[tidx % len(markers)],
                linewidth=1.5,
                color=LINE_COLORS[tidx % len(LINE_COLORS)],
                label=rf"$M/M_0={target:g}$",
            )
        style_axis(ax, ylabel)
    for ax in axes[-1, :]:
        ax.set_xticks(range(len(alphas)))
        ax.set_xticklabels([format_alpha(alpha) for alpha in alphas])
        ax.set_xlabel(r"$\alpha$")
    axes[0, 0].legend(frameon=False, fontsize=8, loc="best")
    fig.suptitle("Q1 passive proliferation: matched-growth metrics", fontsize=20, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)


def draw_boundary(ax, phi, linewidth=0.85):
    if np.nanmin(phi) <= PHI_THRESHOLD <= np.nanmax(phi):
        ax.contour(
            phi.T,
            levels=[PHI_THRESHOLD],
            colors=[BOUNDARY_COLOR],
            linewidths=linewidth,
            origin="lower",
        )


def hide_axes(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal")


def plot_directors(ax, fields):
    phi = fields["phi"]
    theta = fields["theta"]
    lx, ly = phi.shape
    stride = max(5, int(max(lx, ly) / 28))
    x, y = np.meshgrid(np.arange(lx), np.arange(ly), indexing="ij")
    sample = np.zeros_like(phi, dtype=bool)
    sample[::stride, ::stride] = True
    sample &= phi > PHI_THRESHOLD
    if not np.any(sample):
        return
    ax.quiver(
        x[sample],
        y[sample],
        np.cos(theta)[sample],
        np.sin(theta)[sample],
        pivot="middle",
        headwidth=0,
        headlength=0,
        headaxislength=0,
        angles="xy",
        scale_units="xy",
        scale=1.0 / stride,
        width=0.0017,
        color="black",
        alpha=0.85,
    )


def chi_phi_rgba(fields):
    chi = np.clip(fields["chi"], 0.0, 1.0)
    phi_alpha = np.clip(
        (fields["phi"] - PHI_VISIBLE_THRESHOLD) / (1.0 - PHI_VISIBLE_THRESHOLD),
        0.0,
        1.0,
    )
    rgb = plt.get_cmap("viridis")(chi)[..., :3]
    return np.dstack((rgb, phi_alpha))


def pressure_limit(snapshot_specs):
    data = []
    for case, row in snapshot_specs:
        fields = field_set(case["ar"].read_frame(int(row["frame"])))
        pressure = fields["pressure"]
        if pressure is None:
            continue
        mask = fields["phi"] > PHI_THRESHOLD
        if np.any(mask):
            data.append(pressure[mask])
    if not data:
        return None
    merged = np.concatenate(data)
    p1, p99 = np.nanpercentile(merged, [1, 99])
    return max(abs(p1), abs(p99), 1e-12)


def plot_target_snapshots(cases, target, outfile, dpi):
    snapshot_specs = []
    for case in cases:
        row = nearest_target_row(case, target)
        if row is not None:
            snapshot_specs.append((case, row))
    if not snapshot_specs:
        return

    ncols = len(snapshot_specs)
    p_lim = pressure_limit(snapshot_specs)
    fig = plt.figure(figsize=(3.15 * ncols + 1.0, 8.8))
    grid = fig.add_gridspec(
        3,
        ncols + 1,
        width_ratios=[1.0] * ncols + [0.055],
        wspace=0.08,
        hspace=0.08,
    )
    axes = np.array([[fig.add_subplot(grid[row, col]) for col in range(ncols)] for row in range(3)])
    caxes = [fig.add_subplot(grid[row, ncols]) for row in range(3)]
    row_images = [None, None, None]

    for col, (case, row) in enumerate(snapshot_specs):
        fields = field_set(case["ar"].read_frame(int(row["frame"])))
        phi = fields["phi"]
        pressure = fields["pressure"]

        axes[0, col].set_facecolor(VISUALIZATION_BACKGROUND)
        row_images[0] = axes[0, col].imshow(
            np.transpose(chi_phi_rgba(fields), (1, 0, 2)),
            origin="lower",
            interpolation="nearest",
        )
        draw_boundary(axes[0, col], phi)

        row_images[1] = axes[1, col].imshow(
            phi.T,
            origin="lower",
            interpolation="nearest",
            cmap="summer",
        )
        draw_boundary(axes[1, col], phi)
        plot_directors(axes[1, col], fields)

        if pressure is not None and p_lim is not None:
            row_images[2] = axes[2, col].imshow(
                pressure.T,
                origin="lower",
                interpolation="nearest",
                cmap="coolwarm",
                vmin=-p_lim,
                vmax=p_lim,
            )
            draw_boundary(axes[2, col], phi)
        else:
            axes[2, col].text(
                0.5,
                0.5,
                "pressure unavailable",
                ha="center",
                va="center",
                transform=axes[2, col].transAxes,
            )

        axes[0, col].set_title(
            "\n".join(
                [
                    case["label"],
                    rf"$\alpha={format_alpha(case['alpha'])}$",
                    rf"$t={row['time']:.0f}$, $M/M_0={row['M_phi_norm']:.3f}$",
                ]
            ),
            fontsize=12,
        )
        for row_idx in range(3):
            hide_axes(axes[row_idx, col])

    axes[0, 0].set_ylabel(r"$\chi$-$\phi$", fontsize=14)
    axes[1, 0].set_ylabel(r"$\phi$ + director", fontsize=14)
    axes[2, 0].set_ylabel(r"$p$ proxy", fontsize=14)
    caxes[0].axis("off")
    plt.colorbar(row_images[1], cax=caxes[1], label=r"$\phi$")
    if row_images[2] is not None:
        plt.colorbar(row_images[2], cax=caxes[2], label=r"$p$")
    else:
        caxes[2].axis("off")
    fig.suptitle(rf"Q1 matched snapshots, $M_\phi/M_\phi(0)\approx {target:g}$", fontsize=20, y=0.995)
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

    cases = [analyze_case(label, path, args.smooth_window) for label, path in archive_paths]
    cases.sort(key=lambda case: (case["alpha"], case["label"]))

    all_rows = []
    for case in cases:
        all_rows.extend(case["rows"])
    matched_summary = build_matched_summary(cases, args.growth_targets)

    write_csv(all_rows, os.path.join(args.outdir, "q1_metrics_timeseries.csv"))
    write_csv(matched_summary, os.path.join(args.outdir, "q1_matched_growth_summary.csv"))

    plot_timeseries(cases, os.path.join(args.outdir, "q1_timeseries_by_time.png"), args.dpi)
    plot_growth_response(cases, os.path.join(args.outdir, "q1_response_by_growth_fraction.png"), args.dpi)
    plot_matched_metrics(
        matched_summary,
        args.growth_targets,
        os.path.join(args.outdir, "q1_matched_metrics.png"),
        args.dpi,
    )
    for target in args.growth_targets:
        target_name = str(target).replace(".", "p")
        plot_target_snapshots(
            cases,
            target,
            os.path.join(args.outdir, f"q1_snapshots_M{target_name}.png"),
            args.dpi,
        )

    print(f"Loaded {len(cases)} Q1 archives from {args.input_root}")
    for case in cases:
        rows = case["rows"]
        print(
            f"{case['label']}: alpha={case['alpha']:g}, "
            f"frames={len(rows)}, final M/M0={rows[-1]['M_phi_norm']:.4g}"
        )
    print(f"Saved Q1 analysis outputs in: {args.outdir}")


if __name__ == "__main__":
    main()
