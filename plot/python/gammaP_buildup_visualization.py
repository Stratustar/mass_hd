#!/usr/bin/env python3
"""Create phi and pressure field visualizations for the GammaP buildup sweep."""

from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np


OUTPUT_ROOT = Path("output/08jun/gammaP_buildup_sweep_100")
VIS_ROOT = Path("analysis_outputs/08jun/gammaP_buildup_sweep_100/visualizations")
SEED = "seed1001"
FINAL_FRAME = 10000
ANIMATION_GAMMAS = ["GammaP0p001", "GammaP0p02", "GammaP0p1"]
ANIMATION_STRIDE = 500


def parse_gamma(label: str) -> float:
    return float(label.replace("GammaP", "").replace("p", "."))


def frame_number(path: Path) -> int:
    match = re.search(r"frame(\d+)\.json$", path.name)
    if not match:
        raise ValueError(f"Cannot parse frame number from {path}")
    return int(match.group(1))


def load_fields(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open() as handle:
        data = json.load(handle)["data"]
    phi = np.asarray(data["phi"]["value"], dtype=float).reshape(100, 100)
    pressure = -np.asarray(data["sigma_bulk"]["value"], dtype=float).reshape(100, 100)
    return phi, pressure


def gamma_dirs() -> list[Path]:
    return sorted(
        [path for path in OUTPUT_ROOT.glob("GammaP*") if path.is_dir()],
        key=lambda path: parse_gamma(path.name),
    )


def pressure_material(phi: np.ndarray, pressure: np.ndarray) -> np.ndarray:
    return np.where(phi > 0.5, pressure, np.nan)


def field_limits(run_dirs: list[Path]) -> tuple[float, float]:
    phi_max = 0.0
    pressure_max = 0.0
    for gamma_dir in run_dirs:
        frame = gamma_dir / SEED / f"frame{FINAL_FRAME}.json"
        phi, pressure = load_fields(frame)
        phi_max = max(phi_max, float(np.nanmax(phi)))
        pressure_max = max(
            pressure_max, float(np.nanmax(pressure_material(phi, pressure)))
        )
    return max(1.05, phi_max), max(0.01, pressure_max)


def make_final_snapshot_grid() -> None:
    run_dirs = gamma_dirs()
    phi_vmax, pressure_vmax = field_limits(run_dirs)

    pressure_cmap = plt.colormaps["magma"].copy()
    pressure_cmap.set_bad("#f7f7f7")

    fig, axes = plt.subplots(
        len(run_dirs),
        2,
        figsize=(5.6, 2.15 * len(run_dirs)),
        constrained_layout=True,
    )

    phi_image = None
    pressure_image = None
    for row, gamma_dir in enumerate(run_dirs):
        phi, pressure = load_fields(gamma_dir / SEED / f"frame{FINAL_FRAME}.json")
        pmat = pressure_material(phi, pressure)

        ax_phi = axes[row, 0]
        ax_pressure = axes[row, 1]

        phi_image = ax_phi.imshow(phi, origin="lower", vmin=0, vmax=phi_vmax, cmap="viridis")
        pressure_image = ax_pressure.imshow(
            pmat, origin="lower", vmin=0, vmax=pressure_vmax, cmap=pressure_cmap
        )

        ax_phi.set_ylabel(rf"$\Gamma_P={parse_gamma(gamma_dir.name):g}$")
        if row == 0:
            ax_phi.set_title(r"$\phi$", pad=8)
            ax_pressure.set_title(r"$P=-\sigma_{\rm bulk}$", pad=8)
        for axis in (ax_phi, ax_pressure):
            axis.set_xticks([])
            axis.set_yticks([])

    cbar_phi = fig.colorbar(phi_image, ax=axes[:, 0], fraction=0.035, pad=0.02)
    cbar_pressure = fig.colorbar(
        pressure_image, ax=axes[:, 1], fraction=0.035, pad=0.02
    )
    cbar_phi.set_label(r"$\phi$")
    cbar_pressure.set_label(r"$P=-\sigma_{\rm bulk}$")
    fig.savefig(
        VIS_ROOT / "gammaP_phi_pressure_final_seed1001.png",
        dpi=220,
        bbox_inches="tight",
    )
    plt.close(fig)


def selected_frame_numbers() -> list[int]:
    return list(range(0, FINAL_FRAME + 1, ANIMATION_STRIDE))


def selected_limits() -> tuple[float, float]:
    phi_max = 0.0
    pressure_max = 0.0
    for gamma in ANIMATION_GAMMAS:
        run_dir = OUTPUT_ROOT / gamma / SEED
        for number in selected_frame_numbers():
            phi, pressure = load_fields(run_dir / f"frame{number}.json")
            phi_max = max(phi_max, float(np.nanmax(phi)))
            pressure_max = max(
                pressure_max, float(np.nanmax(pressure_material(phi, pressure)))
            )
    return max(1.05, phi_max), max(0.01, pressure_max)


def make_selected_animation() -> None:
    frame_numbers = selected_frame_numbers()
    phi_vmax, pressure_vmax = selected_limits()

    pressure_cmap = plt.colormaps["magma"].copy()
    pressure_cmap.set_bad("#f7f7f7")

    fig, axes = plt.subplots(
        2,
        len(ANIMATION_GAMMAS),
        figsize=(9.0, 5.6),
        constrained_layout=True,
    )

    phi_images = []
    pressure_images = []
    for col, gamma in enumerate(ANIMATION_GAMMAS):
        phi, pressure = load_fields(OUTPUT_ROOT / gamma / SEED / "frame0.json")
        pmat = pressure_material(phi, pressure)

        ax_phi = axes[0, col]
        ax_pressure = axes[1, col]
        phi_image = ax_phi.imshow(phi, origin="lower", vmin=0, vmax=phi_vmax, cmap="viridis")
        pressure_image = ax_pressure.imshow(
            pmat, origin="lower", vmin=0, vmax=pressure_vmax, cmap=pressure_cmap
        )
        ax_phi.set_title(rf"$\Gamma_P={parse_gamma(gamma):g}$")
        if col == 0:
            ax_phi.set_ylabel(r"$\phi$")
            ax_pressure.set_ylabel(r"$P=-\sigma_{\rm bulk}$")
        for axis in (ax_phi, ax_pressure):
            axis.set_xticks([])
            axis.set_yticks([])
        phi_images.append(phi_image)
        pressure_images.append(pressure_image)

    fig.colorbar(phi_images[0], ax=axes[0, :], fraction=0.035, pad=0.02)
    fig.colorbar(pressure_images[0], ax=axes[1, :], fraction=0.035, pad=0.02)
    title = fig.suptitle("")

    def update(frame_index: int):
        number = frame_numbers[frame_index]
        title.set_text(f"{SEED}, t={number}")
        artists = [title]
        for col, gamma in enumerate(ANIMATION_GAMMAS):
            phi, pressure = load_fields(OUTPUT_ROOT / gamma / SEED / f"frame{number}.json")
            phi_images[col].set_data(phi)
            pressure_images[col].set_data(pressure_material(phi, pressure))
            artists.extend([phi_images[col], pressure_images[col]])
        return artists

    ani = animation.FuncAnimation(
        fig, update, frames=len(frame_numbers), interval=220, blit=False
    )
    ani.save(
        VIS_ROOT / "gammaP_phi_pressure_selected_timeseries_seed1001.gif",
        writer=animation.PillowWriter(fps=5),
        dpi=140,
    )
    plt.close(fig)


def main() -> None:
    VIS_ROOT.mkdir(parents=True, exist_ok=True)
    make_final_snapshot_grid()
    make_selected_animation()


if __name__ == "__main__":
    main()
