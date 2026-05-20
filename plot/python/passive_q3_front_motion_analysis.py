import matplotlib
matplotlib.use("Agg")

import argparse
import csv
import json
import os
import sys
import zipfile

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
CHI_FRONT_THRESHOLD = 0.5
BOUNDARY_COLOR = "#d9468f"
BACKGROUND_COLOR = "#7e7e7e"
LINE_COLORS = {
    "C1a": "#255c99",
    "C1b": "#7c3aed",
    "C2a": "#047857",
    "C2b": "#b45309",
    "C3": "#be123c",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Q3 passive proliferation analysis: phenotype front motion."
    )
    parser.add_argument(
        "input_root",
        help="Directory containing Q3 simulation archive directories.",
    )
    parser.add_argument("outdir", help="Directory for CSV and PNG outputs.")
    parser.add_argument(
        "--strip-half-width",
        type=float,
        default=30.0,
        help="Analyze central strip |y-LY/2| <= this width. Default: 30.",
    )
    parser.add_argument(
        "--frame-stride",
        type=int,
        default=1,
        help="Analyze every Nth frame, always including the last frame. Default: 1.",
    )
    parser.add_argument(
        "--xi-bins",
        type=int,
        default=160,
        help="Material-coordinate bins for kymographs. Default: 160.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=220,
        help="PNG resolution. Default: 220.",
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


def case_id(label):
    return os.path.basename(label).split("_", 1)[0]


def format_float(value):
    if not np.isfinite(value):
        return "nan"
    return f"{value:.3g}"


def display_label(case):
    grow = int(case["growTogether"])
    dchi = float(case["Dchi"])
    alpha = float(case["alpha"])
    prefix = case_id(case["label"])
    return rf"{prefix}: gT={grow}, $D_\chi={format_float(dchi)}$, $\alpha={format_float(alpha)}$"


def case_color(case):
    return LINE_COLORS.get(case_id(case["label"]), "#374151")


def case_linestyle(case):
    if int(case["growTogether"]) == 1:
        return "--"
    if float(case["Dchi"]) > 0.0:
        return ":"
    return "-"


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
    chi_from_m = np.divide(m, phi, out=np.zeros_like(phi), where=phi > PHI_VISIBLE_THRESHOLD)
    chi = chi_raw if np.any(np.isfinite(chi_raw)) else chi_from_m
    chi = np.clip(chi, 0.0, 1.0)
    return {"phi": phi, "m": m, "chi": chi}


def read_text_from_archive(path, stem):
    fname = stem + ".json"
    direct = os.path.join(path, fname)
    if os.path.isfile(direct):
        with open(direct) as handle:
            return handle.read()

    zipped = direct + ".zip"
    if os.path.isfile(zipped):
        with zipfile.ZipFile(zipped) as handle:
            return handle.read(fname).decode()

    if path.endswith(".zip") and os.path.isfile(path):
        with zipfile.ZipFile(path) as handle:
            return handle.read(fname).decode()

    raise FileNotFoundError(f"Cannot find {fname} in {path}")


def field_object_from_text(text, name):
    pattern = f'"{name}"'
    key_start = text.find(pattern)
    if key_start < 0:
        return None

    colon = text.find(":", key_start + len(pattern))
    object_start = text.find("{", colon)
    if colon < 0 or object_start < 0:
        raise ValueError(f"Malformed field object for {name}")

    depth = 0
    in_string = False
    escape = False
    for idx in range(object_start, len(text)):
        char = text[idx]
        if in_string:
            if escape:
                escape = False
            elif char == "\\":
                escape = True
            elif char == '"':
                in_string = False
            continue

        if char == '"':
            in_string = True
        elif char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return json.loads(text[object_start : idx + 1])

    raise ValueError(f"Unterminated field object for {name}")


def array_field_from_text(text, name, shape, fallback=None):
    entry = field_object_from_text(text, name)
    if entry is None:
        if fallback is None:
            raise KeyError(f"Missing field {name}")
        return np.array(fallback, dtype=float, copy=True)
    return np.asarray(entry["value"], dtype=float).reshape(shape)


def read_frame_fields(ar, path, frame_index):
    step = int(ar.nstart + frame_index * ar.ninfo)
    text = read_text_from_archive(path, f"frame{step}")
    shape = (int(param(ar, "LX")), int(param(ar, "LY")))
    phi = array_field_from_text(text, "phi", shape)
    m = array_field_from_text(text, "m", shape, fallback=np.zeros(shape, dtype=float))
    chi_entry = field_object_from_text(text, "chi")
    if chi_entry is None:
        chi = np.divide(m, phi, out=np.zeros_like(phi), where=phi > PHI_VISIBLE_THRESHOLD)
    else:
        chi = np.asarray(chi_entry["value"], dtype=float).reshape(shape)
    return {"phi": phi, "m": m, "chi": np.clip(chi, 0.0, 1.0)}


def strip_mask(ly, half_width):
    y = np.arange(ly, dtype=float)
    center = 0.5 * (ly - 1)
    return np.abs(y - center) <= float(half_width)


def strip_profiles(fields, half_width):
    phi = fields["phi"]
    chi = fields["chi"]
    _, ly = phi.shape
    mask = strip_mask(ly, half_width)
    phi_strip = phi[:, mask]
    chi_strip = chi[:, mask]
    mass_x = np.sum(phi_strip, axis=1)
    phi_x = np.mean(phi_strip, axis=1)
    chi_phi_x = np.divide(
        np.sum(phi_strip * chi_strip, axis=1),
        mass_x,
        out=np.full(phi.shape[0], np.nan),
        where=mass_x > 1.0e-12,
    )
    chi_phi_x = np.clip(chi_phi_x, 0.0, 1.0)
    return {
        "phi_x": phi_x,
        "chi_phi_x": chi_phi_x,
        "strip_mass": float(np.sum(phi_strip)),
        "strip_chi_mass": float(np.sum(phi_strip * chi_strip)),
    }


def interpolate_crossing(x0, y0, x1, y1, target):
    if not np.isfinite(y0) or not np.isfinite(y1) or abs(y1 - y0) < 1.0e-14:
        return 0.5 * (x0 + x1)
    return x0 + (target - y0) * (x1 - x0) / (y1 - y0)


def material_edges(phi_x, threshold=PHI_THRESHOLD):
    material = np.asarray(phi_x) >= threshold
    if not np.any(material):
        return np.nan, np.nan

    x = np.arange(len(phi_x), dtype=float)
    first = int(np.flatnonzero(material)[0])
    last = int(np.flatnonzero(material)[-1])

    if first > 0:
        left = interpolate_crossing(
            x[first - 1], phi_x[first - 1], x[first], phi_x[first], threshold
        )
    else:
        left = x[first]

    if last < len(phi_x) - 1:
        right = interpolate_crossing(
            x[last], phi_x[last], x[last + 1], phi_x[last + 1], threshold
        )
    else:
        right = x[last]

    return float(left), float(right)


def front_crossing(chi_x, left, right, threshold=CHI_FRONT_THRESHOLD):
    if not np.isfinite(left) or not np.isfinite(right) or right <= left:
        return np.nan

    x = np.arange(len(chi_x), dtype=float)
    i0 = max(0, int(np.floor(left)))
    i1 = min(len(chi_x) - 1, int(np.ceil(right)))
    section = np.asarray(chi_x[i0 : i1 + 1], dtype=float)
    xs = x[i0 : i1 + 1]
    valid = np.isfinite(section)
    if np.count_nonzero(valid) == 0:
        return np.nan

    if np.nanmin(section[valid]) > threshold:
        return float(right)
    if np.nanmax(section[valid]) < threshold:
        return float(left)

    descending = []
    any_crossing = []
    for idx in range(len(section) - 1):
        y0 = section[idx]
        y1 = section[idx + 1]
        if not np.isfinite(y0) or not np.isfinite(y1):
            continue
        crosses = (y0 - threshold) * (y1 - threshold) <= 0.0 and y0 != y1
        if not crosses:
            continue
        location = interpolate_crossing(xs[idx], y0, xs[idx + 1], y1, threshold)
        any_crossing.append(location)
        if y0 >= threshold and y1 <= threshold:
            descending.append(location)

    candidates = descending if descending else any_crossing
    if candidates:
        center = 0.5 * (left + right)
        return float(candidates[int(np.argmin(np.abs(np.asarray(candidates) - center)))])

    closest = xs[valid][int(np.argmin(np.abs(section[valid] - threshold)))]
    return float(np.clip(closest, left, right))


def resample_material_profile(chi_x, left, right, xi_grid):
    profile = np.full_like(xi_grid, np.nan, dtype=float)
    if not np.isfinite(left) or not np.isfinite(right) or right <= left:
        return profile
    x = np.arange(len(chi_x), dtype=float)
    valid = np.isfinite(chi_x)
    if np.count_nonzero(valid) < 2:
        return profile
    material_x = left + xi_grid * (right - left)
    profile[:] = np.interp(material_x, x[valid], chi_x[valid], left=np.nan, right=np.nan)
    return profile


def finite_gradient(values, times):
    values = np.asarray(values, dtype=float)
    times = np.asarray(times, dtype=float)
    out = np.full_like(values, np.nan)
    valid = np.isfinite(values) & np.isfinite(times)
    if np.count_nonzero(valid) >= 2:
        out[valid] = np.gradient(values[valid], times[valid])
    return out


def frame_metrics(ar, label, frame_index, half_width, first_mass, fields):
    profiles = strip_profiles(fields, half_width)
    left, right = material_edges(profiles["phi_x"])
    x_chi = front_crossing(profiles["chi_phi_x"], left, right)
    width = right - left if np.isfinite(left) and np.isfinite(right) else np.nan
    xi_chi = (x_chi - left) / width if np.isfinite(x_chi) and width > 0.0 else np.nan
    delta_right = right - x_chi if np.isfinite(x_chi) and np.isfinite(right) else np.nan
    mass = float(np.sum(fields["phi"]))

    return {
        "case": label,
        "growTogether": int(param(ar, "growTogether", -1)),
        "Dchi": float(param(ar, "Dchi", np.nan)),
        "alpha": float(param(ar, "alpha", np.nan)),
        "frame": frame_index,
        "step": int(ar.nstart + frame_index * ar.ninfo),
        "time": simulation_time(ar, frame_index),
        "M_phi": mass,
        "M_phi_norm": mass / first_mass if first_mass > 0.0 else np.nan,
        "strip_mass": profiles["strip_mass"],
        "strip_chi_fraction": (
            profiles["strip_chi_mass"] / profiles["strip_mass"]
            if profiles["strip_mass"] > 0.0
            else np.nan
        ),
        "x_phi_left": left,
        "x_phi_right": right,
        "material_width": width,
        "x_chi": x_chi,
        "xi_chi": xi_chi,
        "Delta_right": delta_right,
        "Delta_norm": delta_right / width if width > 0.0 else np.nan,
        "max_phi_x": float(np.nanmax(profiles["phi_x"])),
    }, profiles


def analyze_case(label, path, frame_stride, half_width, xi_grid):
    ar = loadarchive(path)
    first_fields = read_frame_fields(ar, path, 0)
    first_mass = float(np.sum(first_fields["phi"]))
    rows = []
    lab_kymo = []
    material_kymo = []

    for frame_index in frame_indices(ar, frame_stride):
        fields = first_fields if frame_index == 0 else read_frame_fields(ar, path, frame_index)
        row, profiles = frame_metrics(ar, label, frame_index, half_width, first_mass, fields)
        rows.append(row)
        lab_profile = np.array(profiles["chi_phi_x"], dtype=float, copy=True)
        lab_profile[profiles["phi_x"] < PHI_VISIBLE_THRESHOLD] = np.nan
        lab_kymo.append(lab_profile)
        material_kymo.append(
            resample_material_profile(
                profiles["chi_phi_x"],
                row["x_phi_left"],
                row["x_phi_right"],
                xi_grid,
            )
        )

    if rows:
        x0 = rows[0]["x_chi"]
        xi0 = rows[0]["xi_chi"]
        width0 = rows[0]["material_width"]
        times = np.array([row["time"] for row in rows], dtype=float)
        x_chi = np.array([row["x_chi"] for row in rows], dtype=float)
        x_left = np.array([row["x_phi_left"] for row in rows], dtype=float)
        x_right = np.array([row["x_phi_right"] for row in rows], dtype=float)
        xi_chi = np.array([row["xi_chi"] for row in rows], dtype=float)
        delta = np.array([row["Delta_right"] for row in rows], dtype=float)
        for row, v_chi, v_left, v_right, v_xi, v_delta in zip(
            rows,
            finite_gradient(x_chi, times),
            finite_gradient(x_left, times),
            finite_gradient(x_right, times),
            finite_gradient(xi_chi, times),
            finite_gradient(delta, times),
        ):
            row["x_chi_shift"] = row["x_chi"] - x0 if np.isfinite(x0) else np.nan
            row["xi_chi_shift"] = row["xi_chi"] - xi0 if np.isfinite(xi0) else np.nan
            row["material_width_norm"] = (
                row["material_width"] / width0 if np.isfinite(width0) and width0 > 0.0 else np.nan
            )
            row["v_chi"] = float(v_chi)
            row["v_phi_left"] = float(v_left)
            row["v_phi_right"] = float(v_right)
            row["v_xi_chi"] = float(v_xi)
            row["v_delta_right"] = float(v_delta)
            row["v_chi_minus_v_right"] = (
                float(v_chi - v_right) if np.isfinite(v_chi) and np.isfinite(v_right) else np.nan
            )

    return {
        "label": label,
        "path": path,
        "ar": ar,
        "rows": rows,
        "lab_kymo": np.asarray(lab_kymo, dtype=float),
        "material_kymo": np.asarray(material_kymo, dtype=float),
        "alpha": float(param(ar, "alpha", np.nan)),
        "Dchi": float(param(ar, "Dchi", np.nan)),
        "growTogether": int(param(ar, "growTogether", -1)),
    }


def values(rows, key):
    return np.array([row.get(key, np.nan) for row in rows], dtype=float)


def style_axis(ax, ylabel, xlabel=None):
    ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(xlabel)
    ax.grid(True, alpha=0.25, linewidth=0.6)


def plot_core(cases, outfile, dpi):
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.9))
    for case in cases:
        rows = case["rows"]
        color = case_color(case)
        style = case_linestyle(case)
        label = display_label(case)
        axes[0].plot(
            values(rows, "time"),
            values(rows, "x_chi_shift"),
            color=color,
            linestyle=style,
            linewidth=2.0,
            label=label,
        )
        axes[1].plot(
            values(rows, "M_phi_norm"),
            values(rows, "xi_chi"),
            color=color,
            linestyle=style,
            linewidth=2.0,
            label=label,
        )
        axes[2].plot(
            values(rows, "M_phi_norm"),
            values(rows, "Delta_norm"),
            color=color,
            linestyle=style,
            linewidth=2.0,
            label=label,
        )

    style_axis(axes[0], r"$x_\chi(t)-x_\chi(0)$", r"$t$")
    style_axis(axes[1], r"$\xi_\chi=(x_\chi-x_L)/(x_R-x_L)$", r"$M_\phi/M_\phi(0)$")
    style_axis(axes[2], r"$(x_R-x_\chi)/(x_R-x_L)$", r"$M_\phi/M_\phi(0)$")
    axes[1].axhline(0.5, color="black", linewidth=0.8, alpha=0.35)
    axes[2].axhline(0.5, color="black", linewidth=0.8, alpha=0.35)
    axes[0].legend(frameon=False, fontsize=7.5, loc="best")
    fig.suptitle("Q3 core result: front motion in lab and material coordinates", fontsize=17, y=1.03)
    fig.tight_layout()
    fig.savefig(outfile, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def chi_cmap():
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad(BACKGROUND_COLOR)
    return cmap


def imshow_kymo(ax, data, extent, title, cmap):
    image = ax.imshow(
        data,
        origin="lower",
        aspect="auto",
        extent=extent,
        cmap=cmap,
        vmin=0.0,
        vmax=1.0,
        interpolation="nearest",
    )
    ax.set_title(title, fontsize=10)
    ax.set_ylabel(r"$t$")
    return image


def plot_lab_kymograph(cases, outfile, dpi):
    fig, axes = plt.subplots(len(cases), 1, figsize=(10.8, 2.0 * len(cases)), sharex=True)
    if len(cases) == 1:
        axes = [axes]
    cmap = chi_cmap()
    image = None
    for ax, case in zip(axes, cases):
        rows = case["rows"]
        times = values(rows, "time")
        lx = case["lab_kymo"].shape[1]
        extent = [0.0, float(lx - 1), float(times[0]), float(times[-1])]
        image = imshow_kymo(
            ax,
            case["lab_kymo"],
            extent,
            display_label(case),
            cmap,
        )
        ax.plot(values(rows, "x_phi_left"), times, color=BOUNDARY_COLOR, linewidth=0.9)
        ax.plot(values(rows, "x_phi_right"), times, color=BOUNDARY_COLOR, linewidth=0.9)
        ax.plot(values(rows, "x_chi"), times, color="white", linewidth=1.5)
    axes[-1].set_xlabel(r"lab position $x$")
    fig.colorbar(image, ax=axes, fraction=0.018, pad=0.015, label=r"$\langle\chi\rangle_{\phi,y}$")
    fig.suptitle("Q3 lab-frame kymograph: central strip phenotype front", fontsize=17, y=0.995)
    fig.savefig(outfile, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_material_kymograph(cases, outfile, dpi):
    fig, axes = plt.subplots(len(cases), 1, figsize=(8.4, 2.0 * len(cases)), sharex=True)
    if len(cases) == 1:
        axes = [axes]
    cmap = chi_cmap()
    image = None
    for ax, case in zip(axes, cases):
        rows = case["rows"]
        times = values(rows, "time")
        extent = [0.0, 1.0, float(times[0]), float(times[-1])]
        image = imshow_kymo(
            ax,
            case["material_kymo"],
            extent,
            display_label(case),
            cmap,
        )
        ax.plot(values(rows, "xi_chi"), times, color="white", linewidth=1.5)
        ax.axvline(0.5, color="black", linewidth=0.8, alpha=0.35)
        ax.set_xlim(0.0, 1.0)
    axes[-1].set_xlabel(r"material coordinate $\xi=(x-x_L)/(x_R-x_L)$")
    fig.colorbar(image, ax=axes, fraction=0.022, pad=0.015, label=r"$\langle\chi\rangle_{\phi,y}$")
    fig.suptitle("Q3 material-frame kymograph: advection vs invasion", fontsize=17, y=0.995)
    fig.savefig(outfile, dpi=dpi, bbox_inches="tight")
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


def snapshot_indices(ar):
    n_available = available_frame_count(ar)
    return [0, n_available // 2, n_available - 1]


def nearest_row(rows, frame_index):
    frames = values(rows, "frame")
    idx = int(np.argmin(np.abs(frames - frame_index)))
    return rows[idx]


def plot_snapshots(cases, outfile, dpi, half_width):
    fig, axes = plt.subplots(3, len(cases), figsize=(3.0 * len(cases), 8.6))
    if len(cases) == 1:
        axes = axes[:, np.newaxis]

    for col, case in enumerate(cases):
        ar = case["ar"]
        for row_idx, frame_index in enumerate(snapshot_indices(ar)):
            fields = read_frame_fields(ar, case["path"], frame_index)
            row = nearest_row(case["rows"], frame_index)
            ax = axes[row_idx, col]
            ax.set_facecolor(BACKGROUND_COLOR)
            ax.imshow(
                np.transpose(chi_phi_rgba(fields), (1, 0, 2)),
                origin="lower",
                interpolation="nearest",
            )
            draw_boundary(ax, fields["phi"])
            ly = fields["phi"].shape[1]
            y_center = 0.5 * (ly - 1)
            ax.axhline(y_center - half_width, color="white", linewidth=0.6, alpha=0.7)
            ax.axhline(y_center + half_width, color="white", linewidth=0.6, alpha=0.7)
            if np.isfinite(row["x_chi"]):
                ax.axvline(row["x_chi"], color="white", linewidth=1.0)
            hide_axes(ax)
            if row_idx == 0:
                ax.set_title(display_label(case), fontsize=10)
            if col == 0:
                label = ["initial", "middle", "final"][row_idx]
                ax.set_ylabel(label, fontsize=12)

    fig.suptitle(r"Q3 $\chi$-$\phi$ composites with analyzed central strip", fontsize=17, y=0.995)
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
    xi_grid = np.linspace(0.0, 1.0, args.xi_bins)

    archive_paths = discover_archives(args.input_root)
    if not archive_paths:
        raise SystemExit(f"No simulation archives found under {args.input_root}")

    cases = []
    for label, path in archive_paths:
        print(f"Analyzing {label}...", flush=True)
        cases.append(analyze_case(label, path, args.frame_stride, args.strip_half_width, xi_grid))
    cases.sort(key=lambda case: (case_id(case["label"]), case["label"]))

    all_rows = []
    for case in cases:
        all_rows.extend(case["rows"])

    write_csv(all_rows, os.path.join(args.outdir, "q3_front_metrics.csv"))
    plot_core(cases, os.path.join(args.outdir, "q3_core_front_motion.png"), args.dpi)
    plot_lab_kymograph(cases, os.path.join(args.outdir, "q3_lab_kymograph.png"), args.dpi)
    plot_material_kymograph(cases, os.path.join(args.outdir, "q3_material_kymograph.png"), args.dpi)
    plot_snapshots(
        cases,
        os.path.join(args.outdir, "q3_chi_phi_snapshots.png"),
        args.dpi,
        args.strip_half_width,
    )

    print(f"Loaded {len(cases)} Q3 archives from {args.input_root}")
    for case in cases:
        rows = case["rows"]
        first = rows[0]
        final = rows[-1]
        print(
            f"{case['label']}: growTogether={case['growTogether']}, "
            f"Dchi={case['Dchi']:g}, alpha={case['alpha']:g}, frames={len(rows)}, "
            f"x_chi shift={final['x_chi_shift']:.4g}, "
            f"xi_chi {first['xi_chi']:.4g}->{final['xi_chi']:.4g}, "
            f"Delta_norm {first['Delta_norm']:.4g}->{final['Delta_norm']:.4g}"
        )
    print(f"Saved Q3 analysis outputs in: {args.outdir}")


if __name__ == "__main__":
    main()
