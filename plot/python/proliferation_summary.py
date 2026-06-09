import matplotlib
matplotlib.use("Agg")

import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update(
    {
        "mathtext.fontset": "stix",
        "font.family": "STIXGeneral",
        "axes.titlepad": 6,
    }
)


SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from archive.archive import loadarchive


PHI_THRESHOLD = 0.5
PHI_VISIBLE_THRESHOLD = 0.05
BOUNDARY_COLOR = "#d9468f"
EDGE_COLOR = "#f59e0b"
LINE_COLOR = "#255c99"
SECONDARY_LINE_COLOR = "#b45309"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Visualization-only diagnostics for passive proliferation output."
    )
    parser.add_argument("inputdir", help="Simulation archive directory")
    parser.add_argument("outdir", help="Directory for PNG output")
    parser.add_argument(
        "--dpi",
        type=int,
        default=220,
        help="Figure resolution. Default: 220",
    )
    parser.add_argument(
        "--snapshot-count",
        type=int,
        default=5,
        help="Number of snapshot columns. Default: 5",
    )
    parser.add_argument(
        "--radial-bin-width",
        type=float,
        default=1.0,
        help="Radial bin width in lattice units. Default: 1",
    )
    return parser.parse_args()


def available_frame_count(ar):
    return int((ar.nsteps - ar.nstart) / ar.ninfo) + 1


def simulation_time(ar, frame_index):
    step = ar.nstart + frame_index * ar.ninfo
    dt = float(getattr(ar, "time_step", 1.0))
    return step * dt


def param(ar, name, default=np.nan):
    return ar.parameters.get(name, default)


def title_from_params(ar):
    return (
        rf"$\alpha={param(ar, 'alpha', np.nan):g}$, "
        rf"$\phi_c={param(ar, 'phi_critical', np.nan):g}$, "
        rf"$\zeta={param(ar, 'zeta', np.nan):g}$, "
        rf"relax={param(ar, 'relax_steps', 0)}"
    )


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


def center_of_mass(field):
    weights = np.maximum(field, 0.0)
    total = float(np.sum(weights))
    lx, ly = field.shape
    if total <= 0.0:
        return 0.5 * (lx - 1), 0.5 * (ly - 1)
    x, y = np.meshgrid(np.arange(lx), np.arange(ly), indexing="ij")
    return float(np.sum(x * weights) / total), float(np.sum(y * weights) / total)


def radial_distances(shape, center):
    lx, ly = shape
    x, y = np.meshgrid(np.arange(lx), np.arange(ly), indexing="ij")
    return np.sqrt((x - center[0]) ** 2 + (y - center[1]) ** 2)


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
    boundary_points = np.argwhere(boundary)
    if len(boundary_points) == 0:
        return np.nan

    bx = boundary_points[:, 0]
    by = boundary_points[:, 1]
    radii = np.sqrt((bx - center[0]) ** 2 + (by - center[1]) ** 2)
    return float(np.std(radii))


def material_chi(phi, m, chi):
    from_m = np.divide(m, phi, out=np.zeros_like(phi), where=phi > PHI_VISIBLE_THRESHOLD)
    if np.any(np.isfinite(chi)):
        return np.clip(chi, 0.0, 1.0)
    return np.clip(from_m, 0.0, 1.0)


def finite_gradient(values, times):
    values = np.asarray(values, dtype=float)
    times = np.asarray(times, dtype=float)
    if len(values) < 2 or np.any(~np.isfinite(values)) or np.any(~np.isfinite(times)):
        return np.full_like(values, np.nan)
    return np.gradient(values, times)


def radial_profile(field, radial_bin, nbins, weights=None):
    valid = np.isfinite(field)
    if weights is not None:
        valid &= np.isfinite(weights) & (weights > 0.0)
        weighted_values = np.bincount(
            radial_bin[valid], weights=field[valid] * weights[valid], minlength=nbins
        )
        weight_sum = np.bincount(radial_bin[valid], weights=weights[valid], minlength=nbins)
    else:
        weighted_values = np.bincount(radial_bin[valid], weights=field[valid], minlength=nbins)
        weight_sum = np.bincount(radial_bin[valid], minlength=nbins)

    return np.divide(
        weighted_values,
        weight_sum,
        out=np.full(nbins, np.nan),
        where=weight_sum > 0.0,
    )


def snapshot_indices(n_available, snapshot_count):
    if snapshot_count <= 1:
        return [0]
    return sorted(set(np.linspace(0, n_available - 1, snapshot_count, dtype=int).tolist()))


def field_set(frame):
    phi = grid(frame, "phi")
    m = optional_grid(frame, "m", np.zeros_like(phi))
    chi_raw = optional_grid(frame, "chi", np.full_like(phi, np.nan))
    chi = material_chi(phi, m, chi_raw)
    qqxx = optional_grid(frame, "QQxx", np.zeros_like(phi))
    qqyx = optional_grid(frame, "QQyx", np.zeros_like(phi))
    qmag = np.sqrt(qqxx * qqxx + qqyx * qqyx)
    theta = 0.5 * np.arctan2(qqyx, qqxx)
    division_mask = optional_grid(frame, "division_mask", np.zeros_like(phi))
    vx, vy = velocity_fields(frame)
    speed = np.sqrt(vx * vx + vy * vy)
    omega = vorticity(vx, vy)
    return {
        "phi": phi,
        "m": m,
        "chi": chi,
        "qqxx": qqxx,
        "qqyx": qqyx,
        "qmag": qmag,
        "theta": theta,
        "division_mask": division_mask,
        "vx": vx,
        "vy": vy,
        "speed": speed,
        "omega": omega,
    }


def frame_metrics(ar, frame_index, center, radial_bin, nbins):
    frame = ar.read_frame(frame_index)
    fields = field_set(frame)
    phi = fields["phi"]
    m = fields["m"]
    chi = fields["chi"]
    qmag = fields["qmag"]
    speed = fields["speed"]
    omega = fields["omega"]
    division_mask = fields["division_mask"] > 0.5

    phi_c = param(ar, "phi_critical", np.nan)
    mask = phi > PHI_THRESHOLD
    edge_mask = mask & (
        ~(np.roll(mask, 1, axis=0) & np.roll(mask, -1, axis=0))
        | ~(np.roll(mask, 1, axis=1) & np.roll(mask, -1, axis=1))
    )
    area = float(np.sum(mask))
    perimeter = perimeter_4(mask)
    mass_phi = float(np.sum(phi))
    mass_m = float(np.sum(m))
    positive_phi = np.maximum(phi, 0.0)
    material_weight = np.where(phi > PHI_VISIBLE_THRESHOLD, positive_phi, 0.0)
    positive_phi_sum = float(np.sum(positive_phi))

    metrics = {
        "frame": frame_index,
        "step": int(ar.nstart + frame_index * ar.ninfo),
        "time": simulation_time(ar, frame_index),
        "M_phi": mass_phi,
        "M_m": mass_m,
        "chi_bar": mass_m / mass_phi if mass_phi > 0.0 else np.nan,
        "area_phi_05": area,
        "perimeter_phi_05": perimeter,
        "R_eff": np.sqrt(area / np.pi) if area > 0.0 else np.nan,
        "boundary_roughness": boundary_roughness(mask, center),
        "circularity": (
            float(4.0 * np.pi * area / (perimeter * perimeter))
            if perimeter > 0.0
            else np.nan
        ),
        "mean_phi_material": float(np.mean(phi[mask])) if np.any(mask) else np.nan,
        "max_phi": float(np.max(phi)),
        "overcritical_fraction": (
            float(np.sum((phi > phi_c) & mask) / area)
            if np.isfinite(phi_c) and area > 0.0
            else np.nan
        ),
        "Q_mean": (
            float(np.sum(positive_phi * qmag) / positive_phi_sum)
            if positive_phi_sum > 0.0
            else np.nan
        ),
        "Q_edge": float(np.mean(qmag[edge_mask])) if np.any(edge_mask) else np.nan,
        "u_rms": (
            float(np.sqrt(np.sum(positive_phi * speed * speed) / positive_phi_sum))
            if positive_phi_sum > 0.0
            else np.nan
        ),
        "max_speed": float(np.max(speed)),
        "enstrophy": (
            float(0.5 * np.sum(positive_phi * omega * omega) / positive_phi_sum)
            if positive_phi_sum > 0.0
            else np.nan
        ),
        "division_mask_area": float(np.sum(division_mask)),
        "division_mask_phi_mass": float(np.sum(phi[division_mask])),
    }

    profiles = {
        "phi": radial_profile(phi, radial_bin, nbins),
        "chi": radial_profile(chi, radial_bin, nbins, weights=material_weight),
        "qmag": radial_profile(qmag, radial_bin, nbins, weights=material_weight),
        "speed": radial_profile(speed, radial_bin, nbins, weights=material_weight),
    }
    return metrics, profiles


def analyze_archive(ar, radial_bin_width):
    n_available = available_frame_count(ar)
    first = ar.read_frame(0)
    first_phi = grid(first, "phi")
    center = center_of_mass(first_phi)
    r = radial_distances(first_phi.shape, center)
    nbins = int(np.ceil(np.nanmax(r) / radial_bin_width)) + 1
    radial_bin = np.minimum((r / radial_bin_width).astype(int), nbins - 1)
    radial_r = (np.arange(nbins) + 0.5) * radial_bin_width

    rows = []
    radial = {key: [] for key in ("phi", "chi", "qmag", "speed")}
    for frame_index in range(n_available):
        metrics, profiles = frame_metrics(ar, frame_index, center, radial_bin, nbins)
        rows.append(metrics)
        for key, values in profiles.items():
            radial[key].append(values)

    for key in radial:
        radial[key] = np.vstack(radial[key])

    t = np.array([row["time"] for row in rows], dtype=float)
    masses = np.array([row["M_phi"] for row in rows], dtype=float)
    radii = np.array([row["R_eff"] for row in rows], dtype=float)
    growth = (
        finite_gradient(np.log(masses), t)
        if len(masses) > 1 and np.all(masses > 0.0)
        else np.full_like(masses, np.nan)
    )
    radius_speed = finite_gradient(radii, t)
    for row, g_eff, drdt in zip(rows, growth, radius_speed):
        row["g_eff_mass"] = float(g_eff)
        row["dRdt"] = float(drdt)

    return {
        "rows": rows,
        "radial": radial,
        "radial_r": radial_r,
        "center": center,
        "n_available": n_available,
    }


def values(rows, key):
    return np.array([row[key] for row in rows], dtype=float)


def style_timeseries(ax, ylabel):
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.25, linewidth=0.6)


def plot_dashboard(ar, analysis, outfile, dpi):
    rows = analysis["rows"]
    t = values(rows, "time")
    m_phi = values(rows, "M_phi")
    m0 = m_phi[0] if len(m_phi) and m_phi[0] != 0.0 else 1.0
    phi_c = param(ar, "phi_critical", np.nan)

    panels = [
        ("M_phi_norm", r"$M_{\phi}/M_{\phi}(0)$"),
        ("log_M_phi_norm", r"$\log[M_{\phi}/M_{\phi}(0)]$"),
        ("g_eff_mass", r"$d\log M_{\phi}/dt$"),
        ("R_eff", r"$R_{\mathrm{eff}}$"),
        ("dRdt", r"$dR_{\mathrm{eff}}/dt$"),
        ("area_phi_05", r"$A_{\phi>0.5}$"),
        ("boundary_roughness", r"boundary roughness"),
        ("circularity", r"circularity"),
        ("phi_density", r"$\phi$ density"),
        ("overcritical_fraction", r"$\phi>\phi_c$ area fraction"),
        ("Q_mean", r"$\langle |Q| \rangle_{\phi}$"),
        ("Q_edge", r"$\langle |Q| \rangle_{\mathrm{edge}}$"),
        ("chi_bar", r"$\bar{\chi}=M_m/M_{\phi}$"),
        ("speed", r"speed"),
        ("enstrophy", r"$\frac{1}{2}\langle\omega^2\rangle_{\phi}$"),
    ]

    fig, axes = plt.subplots(5, 3, figsize=(12.0, 13.4), sharex=True)
    for ax, (key, ylabel) in zip(axes.flat, panels):
        if key == "M_phi_norm":
            ax.plot(t, m_phi / m0, color=LINE_COLOR, linewidth=1.9)
        elif key == "log_M_phi_norm":
            ax.plot(t, np.log(m_phi / m0), color=LINE_COLOR, linewidth=1.9)
        elif key == "phi_density":
            ax.plot(t, values(rows, "mean_phi_material"), color=LINE_COLOR, linewidth=1.9, label="mean")
            ax.plot(t, values(rows, "max_phi"), color=SECONDARY_LINE_COLOR, linewidth=1.5, label="max")
            if np.isfinite(phi_c):
                ax.axhline(phi_c, color="#374151", linestyle="--", linewidth=1.0, label=r"$\phi_c$")
            ax.legend(frameon=False, fontsize=8, loc="best")
        elif key == "speed":
            ax.plot(t, values(rows, "u_rms"), color=LINE_COLOR, linewidth=1.9, label=r"$u_{\mathrm{rms}}$")
            ax.plot(t, values(rows, "max_speed"), color=SECONDARY_LINE_COLOR, linewidth=1.5, label=r"$u_{\max}$")
            ax.legend(frameon=False, fontsize=8, loc="best")
        else:
            ax.plot(t, values(rows, key), color=LINE_COLOR, linewidth=1.9)
        style_timeseries(ax, ylabel)

    for ax in axes[-1, :]:
        ax.set_xlabel(r"$t$")
    fig.suptitle(f"Passive proliferation diagnostics: {title_from_params(ar)}", fontsize=20, y=0.996)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)


def image_extent(field):
    lx, ly = field.shape
    return (-0.5, lx - 0.5, -0.5, ly - 0.5)


def draw_boundary(ax, phi, color=BOUNDARY_COLOR, linewidth=0.9):
    if np.nanmin(phi) <= PHI_THRESHOLD <= np.nanmax(phi):
        ax.contour(
            phi.T,
            levels=[PHI_THRESHOLD],
            colors=[color],
            linewidths=linewidth,
            origin="lower",
        )


def hide_ticks(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal")


def masked_material(field, phi):
    cmap_field = np.ma.masked_where(phi <= PHI_VISIBLE_THRESHOLD, field)
    return cmap_field


def add_director_glyphs(ax, fields, stride=None, color="#111827", scale=18):
    phi = fields["phi"]
    theta = fields["theta"]
    lx, ly = phi.shape
    if stride is None:
        stride = max(4, int(max(lx, ly) / 18))
    xs, ys = np.meshgrid(np.arange(lx), np.arange(ly), indexing="ij")
    sample = np.zeros_like(phi, dtype=bool)
    sample[::stride, ::stride] = True
    sample &= phi > PHI_THRESHOLD
    if not np.any(sample):
        return
    ax.quiver(
        xs[sample],
        ys[sample],
        np.cos(theta)[sample],
        np.sin(theta)[sample],
        angles="xy",
        pivot="middle",
        color=color,
        scale=scale,
        width=0.0035,
        headlength=0,
        headaxislength=0,
        headwidth=0,
        alpha=0.8,
    )


def row_colorbar(fig, image, cax, label):
    cbar = fig.colorbar(image, cax=cax)
    cbar.set_label(label)


def plot_field_snapshots(ar, analysis, outfile, dpi, snapshot_count):
    n_available = analysis["n_available"]
    indices = snapshot_indices(n_available, snapshot_count)
    ncols = len(indices)
    fig = plt.figure(figsize=(3.15 * ncols + 1.0, 11.0))
    grid_spec = fig.add_gridspec(
        4,
        ncols + 1,
        width_ratios=[1.0] * ncols + [0.055],
        wspace=0.08,
        hspace=0.08,
    )
    axes = np.array(
        [[fig.add_subplot(grid_spec[row, col]) for col in range(ncols)] for row in range(4)]
    )
    colorbar_axes = [fig.add_subplot(grid_spec[row, ncols]) for row in range(4)]

    row_images = [None, None, None, None]
    for col, frame_index in enumerate(indices):
        frame = ar.read_frame(frame_index)
        fields = field_set(frame)
        phi = fields["phi"]
        chi = fields["chi"]
        qmag = fields["qmag"]
        speed = fields["speed"]
        t = simulation_time(ar, frame_index)

        im = axes[0, col].imshow(phi.T, origin="lower", cmap="magma", interpolation="nearest")
        row_images[0] = im
        draw_boundary(axes[0, col], phi)
        axes[0, col].set_title(rf"$t={t:.0f}$", fontsize=14)

        chi_cmap = plt.get_cmap("viridis").copy()
        chi_cmap.set_bad("#7e7e7e")
        im = axes[1, col].imshow(
            masked_material(chi, phi).T,
            origin="lower",
            cmap=chi_cmap,
            vmin=0.0,
            vmax=1.0,
            interpolation="nearest",
        )
        row_images[1] = im
        draw_boundary(axes[1, col], phi)

        im = axes[2, col].imshow(
            masked_material(qmag, phi).T,
            origin="lower",
            cmap="cividis",
            vmin=0.0,
            interpolation="nearest",
        )
        row_images[2] = im
        draw_boundary(axes[2, col], phi)
        add_director_glyphs(axes[2, col], fields)

        im = axes[3, col].imshow(
            masked_material(speed, phi).T,
            origin="lower",
            cmap="inferno",
            vmin=0.0,
            interpolation="nearest",
        )
        row_images[3] = im
        draw_boundary(axes[3, col], phi)

        for row in range(4):
            hide_ticks(axes[row, col])

    row_labels = [r"$\phi$", r"$\chi$", r"$|Q|$ + director", r"$|u|$"]
    for row, label in enumerate(row_labels):
        axes[row, 0].set_ylabel(label, fontsize=14)
        row_colorbar(fig, row_images[row], colorbar_axes[row], label)

    fig.suptitle(f"Field snapshots: {title_from_params(ar)}", fontsize=20, y=0.995)
    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)


def plot_radial_profiles(ar, analysis, outfile, dpi, snapshot_count):
    indices = snapshot_indices(analysis["n_available"], snapshot_count)
    r = analysis["radial_r"]
    radial = analysis["radial"]
    cmap = plt.get_cmap("plasma")
    colors = cmap(np.linspace(0.08, 0.92, len(indices)))
    panels = [
        ("phi", r"$\langle\phi\rangle_{\theta}$"),
        ("chi", r"$\langle\chi\rangle_{\phi,\theta}$"),
        ("qmag", r"$\langle |Q| \rangle_{\phi,\theta}$"),
        ("speed", r"$\langle |u| \rangle_{\phi,\theta}$"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(10.2, 7.6), sharex=True)
    for ax, (key, ylabel) in zip(axes.flat, panels):
        for color, frame_index in zip(colors, indices):
            ax.plot(
                r,
                radial[key][frame_index],
                color=color,
                linewidth=1.8,
                label=rf"$t={simulation_time(ar, frame_index):.0f}$",
            )
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.25, linewidth=0.6)
    axes[1, 0].set_xlabel(r"$r$")
    axes[1, 1].set_xlabel(r"$r$")
    axes[0, 1].legend(frameon=False, fontsize=8, loc="best")
    fig.suptitle(f"Radial profiles: {title_from_params(ar)}", fontsize=20, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)


def plot_radial_kymographs(ar, analysis, outfile, dpi):
    rows = analysis["rows"]
    t = values(rows, "time")
    r = analysis["radial_r"]
    radial = analysis["radial"]
    panels = [
        ("phi", r"$\phi(r,t)$", "magma"),
        ("chi", r"$\chi(r,t)$", "viridis"),
        ("qmag", r"$|Q|(r,t)$", "cividis"),
        ("speed", r"$|u|(r,t)$", "inferno"),
    ]
    extent = (float(r[0]), float(r[-1]), float(t[0]), float(t[-1]))

    fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.8), sharex=True, sharey=True)
    for ax, (key, label, cmap) in zip(axes.flat, panels):
        data = np.ma.masked_invalid(radial[key])
        im = ax.imshow(
            data,
            origin="lower",
            aspect="auto",
            extent=extent,
            cmap=cmap,
            interpolation="nearest",
        )
        ax.set_title(label, fontsize=14)
        ax.set_xlabel(r"$r$")
        ax.set_ylabel(r"$t$")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)

    fig.suptitle(f"Radial kymographs: {title_from_params(ar)}", fontsize=20, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)


def plot_mask_diagnostics(ar, analysis, outfile, dpi, snapshot_count):
    rows = analysis["rows"]
    t = values(rows, "time")
    indices = snapshot_indices(analysis["n_available"], snapshot_count)
    ncols = max(2, len(indices))

    fig = plt.figure(figsize=(3.1 * ncols + 0.8, 6.5))
    grid_spec = fig.add_gridspec(2, ncols)
    ax_area = fig.add_subplot(grid_spec[0, : ncols // 2])
    ax_mass = fig.add_subplot(grid_spec[0, ncols // 2 :])

    ax_area.plot(t, values(rows, "division_mask_area"), color=LINE_COLOR, linewidth=1.9)
    style_timeseries(ax_area, "mask area")
    ax_area.set_xlabel(r"$t$")
    ax_mass.plot(t, values(rows, "division_mask_phi_mass"), color=SECONDARY_LINE_COLOR, linewidth=1.9)
    style_timeseries(ax_mass, r"$\sum_{\mathrm{mask}}\phi$")
    ax_mass.set_xlabel(r"$t$")

    for col, frame_index in enumerate(indices):
        ax = fig.add_subplot(grid_spec[1, col])
        fields = field_set(ar.read_frame(frame_index))
        phi = fields["phi"]
        mask = np.ma.masked_where(fields["division_mask"] <= 0.5, fields["division_mask"])
        ax.imshow(phi.T, origin="lower", cmap="Greys", interpolation="nearest")
        ax.imshow(mask.T, origin="lower", cmap="autumn", vmin=0.0, vmax=1.0, alpha=0.65)
        draw_boundary(ax, phi)
        ax.set_title(rf"$t={simulation_time(ar, frame_index):.0f}$", fontsize=13)
        hide_ticks(ax)

    fig.suptitle(f"Division-mask diagnostics: {title_from_params(ar)}", fontsize=20, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)


def plot_final_state(ar, analysis, outfile, dpi):
    frame_index = analysis["n_available"] - 1
    fields = field_set(ar.read_frame(frame_index))
    phi = fields["phi"]
    chi = fields["chi"]
    qmag = fields["qmag"]
    speed = fields["speed"]
    omega = fields["omega"]
    phi_c = param(ar, "phi_critical", np.nan)
    overcritical = (phi > phi_c).astype(float) if np.isfinite(phi_c) else np.zeros_like(phi)

    panels = [
        ("phi", phi, "magma", r"$\phi$"),
        ("chi", masked_material(chi, phi), "viridis", r"$\chi$"),
        ("qmag", masked_material(qmag, phi), "cividis", r"$|Q|$"),
        ("speed", masked_material(speed, phi), "inferno", r"$|u|$"),
        ("omega", masked_material(omega, phi), "coolwarm", r"$\omega$"),
        ("overcritical", np.ma.masked_where(phi <= PHI_VISIBLE_THRESHOLD, overcritical), "Reds", r"$\phi>\phi_c$"),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(11.5, 7.4))
    for ax, (key, data, cmap, label) in zip(axes.flat, panels):
        if key in ("chi", "overcritical"):
            vmin, vmax = 0.0, 1.0
        elif key == "omega":
            vmax = np.nanmax(np.abs(omega))
            vmin = -vmax
        else:
            vmin, vmax = None, None
        cmap_obj = plt.get_cmap(cmap).copy()
        cmap_obj.set_bad("#7e7e7e")
        im = ax.imshow(
            data.T,
            origin="lower",
            cmap=cmap_obj,
            vmin=vmin,
            vmax=vmax,
            interpolation="nearest",
        )
        draw_boundary(ax, phi)
        if key == "qmag":
            add_director_glyphs(ax, fields)
        ax.set_title(label, fontsize=14)
        hide_ticks(ax)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)

    fig.suptitle(
        rf"Final state, $t={simulation_time(ar, frame_index):.0f}$: {title_from_params(ar)}",
        fontsize=20,
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    ar = loadarchive(args.inputdir)
    analysis = analyze_archive(ar, args.radial_bin_width)

    outputs = {
        "dashboard": os.path.join(args.outdir, "proliferation_dashboard.png"),
        "snapshots": os.path.join(args.outdir, "field_snapshots.png"),
        "radial_profiles": os.path.join(args.outdir, "radial_profiles.png"),
        "radial_kymographs": os.path.join(args.outdir, "radial_kymographs.png"),
        "mask_diagnostics": os.path.join(args.outdir, "mask_diagnostics.png"),
        "final_state": os.path.join(args.outdir, "final_state_diagnostics.png"),
    }

    plot_dashboard(ar, analysis, outputs["dashboard"], args.dpi)
    plot_field_snapshots(ar, analysis, outputs["snapshots"], args.dpi, args.snapshot_count)
    plot_radial_profiles(ar, analysis, outputs["radial_profiles"], args.dpi, args.snapshot_count)
    plot_radial_kymographs(ar, analysis, outputs["radial_kymographs"], args.dpi)
    plot_mask_diagnostics(ar, analysis, outputs["mask_diagnostics"], args.dpi, args.snapshot_count)
    plot_final_state(ar, analysis, outputs["final_state"], args.dpi)

    print(f"Loaded {analysis['n_available']} frames from {args.inputdir}")
    for label, path in outputs.items():
        print(f"Saved {label}: {path}")


if __name__ == "__main__":
    main()
