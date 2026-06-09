#!/usr/bin/env python3
"""Analyze friction dependence of phi and bulk-pressure buildup."""

from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


OUTPUT_ROOT = Path("output/08jun/friction_buildup_sweep_100")
ANALYSIS_ROOT = Path("analysis_outputs/08jun/friction_buildup_sweep_100")
PHI_CRITICAL = 1.2


@dataclass
class RunMetrics:
    friction: float
    seed: int
    phi_p99_buildup: float
    phi_max_peak: float
    phi_critical_excess_peak: float
    pressure_p99_buildup: float
    pressure_max_peak: float
    final_phi_p99: float
    final_phi_max: float
    final_pressure_p99: float
    final_pressure_max: float


def parse_friction(label: str) -> float:
    return float(label.replace("friction", "").replace("p", "."))


def frame_number(path: Path) -> int:
    match = re.search(r"frame(\d+)\.json$", path.name)
    if not match:
        raise ValueError(f"Cannot parse frame number from {path}")
    return int(match.group(1))


def field(data: dict, name: str) -> np.ndarray:
    return np.asarray(data[name]["value"], dtype=float)


def load_frame(path: Path) -> dict:
    with path.open() as handle:
        return json.load(handle)["data"]


def analyze_run(run_dir: Path) -> RunMetrics:
    friction = parse_friction(run_dir.parent.name)
    seed = int(run_dir.name.replace("seed", ""))
    frames = sorted(run_dir.glob("frame*.json"), key=frame_number)
    if not frames:
        raise FileNotFoundError(f"No frames found in {run_dir}")

    phi0_p99 = None
    pressure0_p99 = None
    phi_p99_peak = -np.inf
    phi_max_peak = -np.inf
    pressure_p99_peak = -np.inf
    pressure_max_peak = -np.inf
    phi_critical_excess_peak = 0.0
    final_phi_p99 = np.nan
    final_phi_max = np.nan
    final_pressure_p99 = np.nan
    final_pressure_max = np.nan

    for path in frames:
        data = load_frame(path)
        phi = field(data, "phi")
        pressure = -field(data, "sigma_bulk")
        material = phi > 0.5
        phi_material = phi[material] if np.any(material) else phi
        pressure_material = pressure[material] if np.any(material) else pressure

        phi_p99 = float(np.quantile(phi_material, 0.99))
        phi_max = float(np.max(phi_material))
        pressure_p99 = float(np.quantile(pressure_material, 0.99))
        pressure_max = float(np.max(pressure_material))

        if phi0_p99 is None:
            phi0_p99 = phi_p99
            pressure0_p99 = pressure_p99

        phi_p99_peak = max(phi_p99_peak, phi_p99)
        phi_max_peak = max(phi_max_peak, phi_max)
        pressure_p99_peak = max(pressure_p99_peak, pressure_p99)
        pressure_max_peak = max(pressure_max_peak, pressure_max)
        phi_critical_excess_peak = max(
            phi_critical_excess_peak, max(0.0, phi_max - PHI_CRITICAL)
        )

        if frame_number(path) == 10000:
            final_phi_p99 = phi_p99
            final_phi_max = phi_max
            final_pressure_p99 = pressure_p99
            final_pressure_max = pressure_max

    return RunMetrics(
        friction=friction,
        seed=seed,
        phi_p99_buildup=phi_p99_peak - float(phi0_p99),
        phi_max_peak=phi_max_peak,
        phi_critical_excess_peak=phi_critical_excess_peak,
        pressure_p99_buildup=pressure_p99_peak - float(pressure0_p99),
        pressure_max_peak=pressure_max_peak,
        final_phi_p99=final_phi_p99,
        final_phi_max=final_phi_max,
        final_pressure_p99=final_pressure_p99,
        final_pressure_max=final_pressure_max,
    )


def write_csv(rows: list[RunMetrics], path: Path) -> None:
    fields = list(RunMetrics.__dataclass_fields__)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: getattr(row, name) for name in fields})


def grouped(rows: list[RunMetrics], metric: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    frictions = np.asarray(sorted({row.friction for row in rows}), dtype=float)
    means = []
    sems = []
    for friction in frictions:
        vals = np.asarray([getattr(row, metric) for row in rows if row.friction == friction])
        means.append(float(vals.mean()))
        sems.append(float(vals.std(ddof=1) / np.sqrt(len(vals))) if len(vals) > 1 else 0.0)
    return frictions, np.asarray(means), np.asarray(sems)


def write_summary(rows: list[RunMetrics], path: Path) -> None:
    metrics = [
        "phi_p99_buildup",
        "phi_max_peak",
        "phi_critical_excess_peak",
        "pressure_p99_buildup",
        "pressure_max_peak",
    ]
    frictions = sorted({row.friction for row in rows})
    with path.open("w") as handle:
        handle.write("Friction sweep buildup summary\n")
        handle.write("Metrics use material pixels phi > 0.5; pressure = -sigma_bulk.\n\n")
        handle.write("Empirical log-linear fits over group means:\n")
        handle.write("  y = a + b log10(friction / 10)\n")
        for metric in ("phi_p99_buildup", "pressure_p99_buildup"):
            x, y, _ = grouped(rows, metric)
            scaled_log = np.log10(x / 10.0)
            b, a = np.polyfit(scaled_log, y, 1)
            yhat = a + b * scaled_log
            ss_res = float(np.sum((y - yhat) ** 2))
            ss_tot = float(np.sum((y - np.mean(y)) ** 2))
            r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
            handle.write(f"  {metric}: a={a:.6g}, b={b:.6g}, R2={r2:.4f}\n")
        handle.write("\n")
        for friction in frictions:
            handle.write(f"friction = {friction:g}\n")
            sub = [row for row in rows if row.friction == friction]
            for metric in metrics:
                vals = np.asarray([getattr(row, metric) for row in sub], dtype=float)
                handle.write(
                    f"  {metric}: mean={vals.mean():.6g}, "
                    f"sd={vals.std(ddof=1):.6g}, n={len(vals)}\n"
                )
            handle.write("\n")


def plot_summary(rows: list[RunMetrics], path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.4), constrained_layout=True)

    x, y, e = grouped(rows, "phi_p99_buildup")
    axes[0].errorbar(x, y, yerr=e, marker="o", linewidth=1.6, capsize=3)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("friction")
    axes[0].set_ylabel(r"Peak $\Delta \phi_{99}$")
    axes[0].set_title("Density buildup")

    x, y, e = grouped(rows, "pressure_p99_buildup")
    axes[1].errorbar(x, y, yerr=e, marker="o", linewidth=1.6, capsize=3, color="#b44a3a")
    axes[1].set_xscale("log")
    axes[1].set_xlabel("friction")
    axes[1].set_ylabel(r"Peak $\Delta P_{99}$")
    axes[1].set_title(r"Compressive bulk pressure")

    for axis in axes:
        axis.grid(True, color="0.88", linewidth=0.8)
        axis.tick_params(direction="out")

    fig.savefig(path, dpi=220)
    plt.close(fig)


def main() -> None:
    ANALYSIS_ROOT.mkdir(parents=True, exist_ok=True)
    rows = [
        analyze_run(run_dir)
        for run_dir in sorted(OUTPUT_ROOT.glob("friction*/seed*"))
        if (run_dir / "frame10000.json").exists()
    ]
    if not rows:
        raise RuntimeError(f"No completed runs found under {OUTPUT_ROOT}")

    write_csv(rows, ANALYSIS_ROOT / "friction_buildup_metrics.csv")
    write_summary(rows, ANALYSIS_ROOT / "friction_buildup_summary.txt")
    plot_summary(rows, ANALYSIS_ROOT / "friction_buildup_summary.png")


if __name__ == "__main__":
    main()
