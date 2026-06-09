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
VISUALIZATION_BACKGROUND = "#7e7e7e"
LINE_COLOR = "#255c99"
SECONDARY_LINE_COLOR = "#b45309"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Block 1 analysis: proliferation timescale vs relaxation timescale."
    )
    parser.add_argument("inputdir", help="Simulation archive directory")
    parser.add_argument("outdir", help="Directory for PNG output")
    parser.add_argument(
        "--tau-relax",
        type=float,
        default=None,
        help="Relaxation time used to compute P_number for formal growth runs.",
    )
    parser.add_argument("--dpi", type=int, default=220, help="PNG resolution")
    parser.add_argument(
        "--snapshot-count",
        type=int,
        default=5,
        help="Number of snapshot columns.",
    )
    parser.add_argument(
        "--smooth-window",
        type=int,
        default=7,
        help="Odd moving-average window for g_eff. Default: 7.",
    )
    return parser.parse_args()


def available_frame_count(ar):
    return int((ar.nsteps - ar.nstart) / ar.ninfo) + 1


def simulation_time(ar, frame_index):
    return (ar.nstart + frame_index * ar.ninfo) * float(getattr(ar, "time_step", 1.0))


def param(ar, name, default=np.nan):
    return ar.parameters.get(name, default)


def title_from_params(ar):
    return (
        rf"$\alpha={param(ar, 'alpha', np.nan):g}$, "
        rf"$\phi_c={param(ar, 'phi_critical', np.nan):g}$"
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
    return np.nan_to_num(vx).reshape((lx, ly)), np.nan_to_num(vy).reshape((lx, ly))


def vorticity(vx, vy):
    dvdx = 0.5 * (np.roll(vy, -1, axis=0) - np.roll(vy, 1, axis=0))
    dudy = 0.5 * (np.roll(vx, -1, axis=1) - np.roll(vx, 1, axis=1))
    return dvdx - dudy


def pressure_proxy(frame, phi):
    if hasattr(frame, "pressure"):
        return grid(frame, "pressure")
    if hasattr(frame, "sigma_bulk"):
        return -grid(frame, "sigma_bulk")
    phi_c = param_from_frame(frame, "phi_critical", 1.0)
    B = param_from_frame(frame, "B", 0.0)
    return B * np.maximum(phi - phi_c, 0.0)


def param_from_frame(frame, name, default=np.nan):
    return frame.parameters.get(name, default)


def field_set(frame):
    phi = grid(frame, "phi")
    m = optional_grid(frame, "m", np.zeros_like(phi))
    chi_raw = optional_grid(frame, "chi", np.full_like(phi, np.nan))
    chi = np.divide(m, phi, out=np.zeros_like(phi), where=phi > PHI_VISIBLE_THRESHOLD)
    if np.any(np.isfinite(chi_raw)):
        chi = chi_raw
    chi = np.clip(chi, 0.0, 1.0)
    qxx = optional_grid(frame, "QQxx", np.zeros_like(phi))
    qyx = optional_grid(frame, "QQyx", np.zeros_like(phi))
    qmag = np.sqrt(qxx * qxx + qyx * qyx)
    theta = 0.5 * np.arctan2(qyx, qxx)
    pressure = pressure_proxy(frame, phi)
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


def radial_geometry(shape, center, bin_width=1.0):
    lx, ly = shape
    x, y = np.meshgrid(np.arange(lx), np.arange(ly), indexing="ij")
    dx = x - center[0]
    dy = y - center[1]
    r = np.sqrt(dx * dx + dy * dy)
    rhat_x = np.divide(dx, r, out=np.zeros_like(r), where=r > 0.0)
    rhat_y = np.divide(dy, r, out=np.zeros_like(r), where=r > 0.0)
    nbins = int(np.ceil(np.nanmax(r) / bin_width)) + 1
    radial_bin = np.minimum((r / bin_width).astype(int), nbins - 1)
    radial_r = (np.arange(nbins) + 0.5) * bin_width
    return radial_bin, radial_r, rhat_x, rhat_y, nbins


def radial_profile(field, radial_bin, nbins, weights=None):
    valid = np.isfinite(field)
    if weights is not None:
        valid &= np.isfinite(weights) & (weights > 0.0)
        numerator = np.bincount(
            radial_bin[valid], weights=field[valid] * weights[valid], minlength=nbins
        )
        denominator = np.bincount(radial_bin[valid], weights=weights[valid], minlength=nbins)
    else:
        numerator = np.bincount(radial_bin[valid], weights=field[valid], minlength=nbins)
        denominator = np.bincount(radial_bin[valid], minlength=nbins)
    return np.divide(numerator, denominator, out=np.full(nbins, np.nan), where=denominator > 0.0)


def perimeter_4(mask):
    padded = np.pad(mask.astype(bool), 1, mode="constant", constant_values=False)
    center = padded[1:-1, 1:-1]
    exposed = (
        (center & ~padded[2:, 1:-1])
        + (center & ~padded[:-2, 1:-1])
        + (center & ~padded[1:-1, 2:])
        + (center & ~padded[1:-1, :-2])
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


def moving_average(values, window):
    values = np.asarray(values, dtype=float)
    if window <= 1:
        return values
    window = int(window)
    if window % 2 == 0:
        window += 1
    if len(values) < window:
        return values
    kernel = np.ones(window, dtype=float) / window
    padded = np.pad(values, window // 2, mode="edge")
    return np.convolve(padded, kernel, mode="valid")


def analyze_archive(ar, smooth_window):
    n_available = available_frame_count(ar)
    first_fields = field_set(ar.read_frame(0))
    center = center_of_mass(first_fields["phi"])
    radial_bin, radial_r, rhat_x, rhat_y, nbins = radial_geometry(first_fields["phi"].shape, center)
    rows = []
    radial = {key: [] for key in ("phi", "pressure", "ur", "qmag")}

    for frame_index in range(n_available):
        frame = ar.read_frame(frame_index)
        fields = field_set(frame)
        phi = fields["phi"]
        pressure = fields["pressure"]
        qmag = fields["qmag"]
        mask = phi > PHI_THRESHOLD
        material_weight = np.where(phi > PHI_VISIBLE_THRESHOLD, np.maximum(phi, 0.0), 0.0)
        weight_sum = float(np.sum(material_weight))
        pressure_material = pressure[mask]
        pressure_prime = pressure_material - np.mean(pressure_material) if np.any(mask) else np.array([])
        speed2 = fields["speed"] * fields["speed"]
        omega = fields["omega"]
        area = float(np.sum(mask))
        perimeter = perimeter_4(mask)
        qxx_global = float(np.sum(material_weight * fields["qxx"]) / weight_sum) if weight_sum > 0 else np.nan
        qyx_global = float(np.sum(material_weight * fields["qyx"]) / weight_sum) if weight_sum > 0 else np.nan
        edge_mask = (phi > 0.4) & (phi < 0.6)

        rows.append(
            {
                "frame": frame_index,
                "step": int(ar.nstart + frame_index * ar.ninfo),
                "time": simulation_time(ar, frame_index),
                "M_phi": float(np.sum(phi)),
                "R_eff": np.sqrt(area / np.pi) if area > 0.0 else np.nan,
                "dRdt": np.nan,
                "area_phi05": area,
                "perimeter_phi05": perimeter,
                "roughness_norm": (
                    boundary_roughness(mask, center) / np.sqrt(area / np.pi)
                    if area > 0.0
                    else np.nan
                ),
                "circularity": (
                    float(4.0 * np.pi * area / (perimeter * perimeter))
                    if perimeter > 0.0
                    else np.nan
                ),
                "p_rms": float(np.sqrt(np.mean(pressure_prime * pressure_prime)))
                if len(pressure_prime)
                else np.nan,
                "p95": float(np.percentile(pressure_material, 95)) if np.any(mask) else np.nan,
                "p5": float(np.percentile(pressure_material, 5)) if np.any(mask) else np.nan,
                "sigma_phi": float(np.std(phi[mask])) if np.any(mask) else np.nan,
                "u_rms": (
                    float(np.sqrt(np.sum(material_weight * speed2) / weight_sum))
                    if weight_sum > 0
                    else np.nan
                ),
                "enstrophy": (
                    float(0.5 * np.sum(material_weight * omega * omega) / weight_sum)
                    if weight_sum > 0
                    else np.nan
                ),
                "S_mean": (
                    float(np.sum(material_weight * qmag) / weight_sum)
                    if weight_sum > 0
                    else np.nan
                ),
                "S_global": float(np.sqrt(qxx_global * qxx_global + qyx_global * qyx_global)),
                "S_edge": float(np.mean(qmag[edge_mask])) if np.any(edge_mask) else np.nan,
            }
        )

        ur = fields["vx"] * rhat_x + fields["vy"] * rhat_y
        radial["phi"].append(radial_profile(phi, radial_bin, nbins))
        radial["pressure"].append(radial_profile(pressure, radial_bin, nbins, weights=material_weight))
        radial["ur"].append(radial_profile(ur, radial_bin, nbins, weights=material_weight))
        radial["qmag"].append(radial_profile(qmag, radial_bin, nbins, weights=material_weight))

    for key in radial:
        radial[key] = np.vstack(radial[key])

    t = values(rows, "time")
    mass = values(rows, "M_phi")
    radii = values(rows, "R_eff")
    if len(rows) > 1 and np.all(mass > 0):
        g_eff = np.gradient(np.log(mass), t)
        g_eff = moving_average(g_eff, smooth_window)
        dRdt = moving_average(np.gradient(radii, t), smooth_window)
    else:
        g_eff = np.full(len(rows), np.nan)
        dRdt = np.full(len(rows), np.nan)
    for row, g, v in zip(rows, g_eff, dRdt):
        row["g_eff_emp"] = float(g)
        row["tau_grow_emp"] = float(1.0 / g) if g > 0 else np.nan
        row["dRdt"] = float(v)

    return {"rows": rows, "radial": radial, "radial_r": radial_r, "n_available": n_available}


def values(rows, key):
    return np.array([row[key] for row in rows], dtype=float)


def snapshot_indices(n_available, snapshot_count):
    if snapshot_count <= 1:
        return [0]
    return sorted(set(np.linspace(0, n_available - 1, snapshot_count, dtype=int).tolist()))


def draw_boundary(ax, phi):
    if np.nanmin(phi) <= PHI_THRESHOLD <= np.nanmax(phi):
        ax.contour(phi.T, levels=[PHI_THRESHOLD], colors=[BOUNDARY_COLOR], linewidths=0.8, origin="lower")


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
    )


def chi_phi_rgba(fields):
    chi = np.clip(fields["chi"], 0.0, 1.0)
    phi_alpha = np.clip(fields["phi"], 0.0, 1.0)
    rgb = plt.get_cmap("viridis")(1.0 - chi)[..., :3]
    return np.dstack((rgb, phi_alpha))


def robust_pressure_limit(ar, indices):
    values_p = []
    for frame_index in indices:
        fields = field_set(ar.read_frame(frame_index))
        pressure = fields["pressure"]
        mask = fields["phi"] > PHI_THRESHOLD
        if np.any(mask):
            values_p.append(pressure[mask])
    if not values_p:
        return 1.0
    merged = np.concatenate(values_p)
    p1, p99 = np.percentile(merged, [1, 99])
    return max(abs(p1), abs(p99), 1e-12)


def plot_snapshots(ar, analysis, outfile, dpi, snapshot_count):
    indices = snapshot_indices(analysis["n_available"], snapshot_count)
    ncols = len(indices)
    p_lim = robust_pressure_limit(ar, indices)
    fig = plt.figure(figsize=(3.15 * ncols + 1.0, 8.8))
    grid = fig.add_gridspec(3, ncols + 1, width_ratios=[1.0] * ncols + [0.055], wspace=0.08, hspace=0.08)
    axes = np.array([[fig.add_subplot(grid[row, col]) for col in range(ncols)] for row in range(3)])
    caxes = [fig.add_subplot(grid[row, ncols]) for row in range(3)]
    row_images = [None, None, None]

    for col, frame_index in enumerate(indices):
        fields = field_set(ar.read_frame(frame_index))
        phi = fields["phi"]
        pressure = fields["pressure"]
        axes[0, col].set_facecolor(VISUALIZATION_BACKGROUND)
        row_images[0] = axes[0, col].imshow(np.transpose(chi_phi_rgba(fields), (1, 0, 2)), origin="lower", interpolation="nearest")
        draw_boundary(axes[0, col], phi)
        axes[0, col].set_title(rf"$t={simulation_time(ar, frame_index):.0f}$", fontsize=14)

        row_images[1] = axes[1, col].imshow(phi.T, origin="lower", interpolation="nearest", cmap="summer")
        draw_boundary(axes[1, col], phi)
        plot_directors(axes[1, col], fields)

        row_images[2] = axes[2, col].imshow(
            pressure.T,
            origin="lower",
            interpolation="nearest",
            cmap="coolwarm",
            vmin=-p_lim,
            vmax=p_lim,
        )
        draw_boundary(axes[2, col], phi)

        for row in range(3):
            hide_axes(axes[row, col])

    axes[0, 0].set_ylabel(r"$\chi$-$\phi$", fontsize=14)
    axes[1, 0].set_ylabel(r"$\phi$ + director", fontsize=14)
    axes[2, 0].set_ylabel(r"$p=-\sigma_{\mathrm{bulk}}$", fontsize=14)
    caxes[0].axis("off")
    plt.colorbar(row_images[1], cax=caxes[1], label=r"$\phi$")
    plt.colorbar(row_images[2], cax=caxes[2], label=r"$p$")
    fig.suptitle(f"Block 1 snapshots: {title_from_params(ar)}", fontsize=20, y=0.995)
    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)


def style_axis(ax, ylabel):
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.25, linewidth=0.6)


def plot_timescale(ar, analysis, outfile, dpi, tau_relax):
    rows = analysis["rows"]
    t = values(rows, "time")
    g = values(rows, "g_eff_emp")
    tau_g = values(rows, "tau_grow_emp")
    fig, axes = plt.subplots(3 if tau_relax else 2, 1, figsize=(7.4, 7.4 if tau_relax else 5.2), sharex=True)
    axes = np.atleast_1d(axes)
    axes[0].plot(t, g, color=LINE_COLOR, linewidth=1.9)
    style_axis(axes[0], r"$g_{\mathrm{eff,emp}}$")
    axes[1].plot(t, tau_g, color=SECONDARY_LINE_COLOR, linewidth=1.9)
    style_axis(axes[1], r"$\tau_{\mathrm{grow,emp}}$")
    if tau_relax:
        axes[2].plot(t, g * tau_relax, color="#7c3aed", linewidth=1.9)
        axes[2].axhline(0.1, color="#374151", linestyle="--", linewidth=1.0)
        axes[2].axhline(1.0, color="#374151", linestyle=":", linewidth=1.0)
        style_axis(axes[2], r"$P=g_{\mathrm{eff}}\tau_{\mathrm{relax}}$")
    axes[-1].set_xlabel(r"$t$")
    fig.suptitle(f"Timescale response: {title_from_params(ar)}", fontsize=18, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)


def plot_response(rows, panels, outfile, title, dpi):
    t = values(rows, "time")
    fig, axes = plt.subplots(2, 2, figsize=(9.4, 6.8), sharex=True)
    for ax, spec in zip(axes.flat, panels):
        if len(spec) == 2:
            key, ylabel = spec
            ax.plot(t, values(rows, key), color=LINE_COLOR, linewidth=1.9)
        else:
            keys, ylabel, labels = spec
            colors = [LINE_COLOR, SECONDARY_LINE_COLOR, "#7c3aed"]
            for key, label, color in zip(keys, labels, colors):
                ax.plot(t, values(rows, key), color=color, linewidth=1.7, label=label)
            ax.legend(frameon=False, fontsize=8, loc="best")
        style_axis(ax, ylabel)
    axes[1, 0].set_xlabel(r"$t$")
    axes[1, 1].set_xlabel(r"$t$")
    fig.suptitle(title, fontsize=18, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)


def plot_radial(analysis, outfile, title, dpi):
    rows = analysis["rows"]
    t = values(rows, "time")
    r = analysis["radial_r"]
    extent = (float(r[0]), float(r[-1]), float(t[0]), float(t[-1]))
    panels = [
        ("phi", r"$\phi(r,t)$", "magma"),
        ("pressure", r"$p(r,t)$", "coolwarm"),
        ("ur", r"$u_r(r,t)$", "coolwarm"),
        ("qmag", r"$|Q|(r,t)$", "cividis"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.6), sharex=True, sharey=True)
    for ax, (key, label, cmap) in zip(axes.flat, panels):
        data = np.ma.masked_invalid(analysis["radial"][key])
        kwargs = {}
        if key in ("pressure", "ur"):
            lim = np.nanpercentile(np.abs(data.compressed()), 99) if data.count() else 1.0
            kwargs.update({"vmin": -lim, "vmax": lim})
        im = ax.imshow(data, origin="lower", aspect="auto", extent=extent, cmap=cmap, interpolation="nearest", **kwargs)
        ax.set_title(label, fontsize=14)
        ax.set_xlabel(r"$r$")
        ax.set_ylabel(r"$t$")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    fig.suptitle(title, fontsize=18, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    ar = loadarchive(args.inputdir)
    analysis = analyze_archive(ar, args.smooth_window)
    rows = analysis["rows"]

    plot_snapshots(ar, analysis, os.path.join(args.outdir, "snapshots_block1.png"), args.dpi, args.snapshot_count)

    tau_relax = args.tau_relax

    plot_timescale(ar, analysis, os.path.join(args.outdir, "timescale_response.png"), args.dpi, tau_relax)
    plot_response(
        rows,
        [
            (("p95", "p5"), r"pressure percentiles", [r"$p_{95}$", r"$p_5$"]),
            ("p_rms", r"$p_{\mathrm{rms}}$"),
            ("sigma_phi", r"$\sigma_{\phi}$"),
            (("u_rms", "enstrophy"), r"flow response", [r"$u_{\mathrm{rms}}$", r"enstrophy"]),
        ],
        os.path.join(args.outdir, "hydro_pressure_response.png"),
        f"Hydro-pressure response: {title_from_params(ar)}",
        args.dpi,
    )
    plot_response(
        rows,
        [
            ("R_eff", r"$R_{\mathrm{eff}}$"),
            ("dRdt", r"$dR_{\mathrm{eff}}/dt$"),
            ("roughness_norm", r"$w_R/R_{\mathrm{eff}}$"),
            ("circularity", r"circularity"),
        ],
        os.path.join(args.outdir, "morphology_response.png"),
        f"Morphology response: {title_from_params(ar)}",
        args.dpi,
    )
    plot_response(
        rows,
        [
            ("S_mean", r"$S_{\mathrm{mean}}$"),
            ("S_global", r"$S_{\mathrm{global}}$"),
            ("S_edge", r"$S_{\mathrm{edge}}$"),
            (("S_mean", "S_global", "S_edge"), r"nematic order", [r"$S_{\mathrm{mean}}$", r"$S_{\mathrm{global}}$", r"$S_{\mathrm{edge}}$"]),
        ],
        os.path.join(args.outdir, "nematic_response.png"),
        f"Nematic response: {title_from_params(ar)}",
        args.dpi,
    )
    plot_radial(
        analysis,
        os.path.join(args.outdir, "radial_response.png"),
        f"Radial response: {title_from_params(ar)}",
        args.dpi,
    )

    print(f"Loaded {analysis['n_available']} frames from {args.inputdir}")
    if tau_relax is not None:
        print(f"tau_relax: {tau_relax:.8g}")
    print(f"Saved PNG outputs in: {args.outdir}")


if __name__ == "__main__":
    main()
