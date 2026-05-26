#!/usr/bin/env python3
"""
Plot Wu fast-start Meep results from the current directory.

Expected files:

  reference_incident.txt
  theta_080_phi_000.txt
  theta_080_phi_020.txt
  ...

Uncomment additional phi values in PHI_VALUES as those runs finish.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt


THETA = 80
PHI_RE = re.compile(r"theta_(?P<theta>\d+)_phi_(?P<phi>\d+)\.txt$")

PHI_VALUES = [
    0,
    20,
    40,
    60,
    #80,
    # 90,
    # 100,
    # 120,
    # 140,
    # 160,
    # 180,
]


def parse_flux_file(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []

    with path.open(newline="") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line.startswith("flux"):
                continue

            _, values_text = line.split(":", 1)
            values = []
            for value in next(csv.reader([values_text])):
                value = value.strip()
                if not value:
                    continue
                try:
                    values.append(float(value))
                except ValueError:
                    continue

            if len(values) < 3:
                continue

            frequency, reflected, transmitted = values[:3]
            rows.append(
                {
                    "frequency": frequency,
                    "wavelength_um": 1.0 / frequency,
                    "reflected": reflected,
                    "transmitted": transmitted,
                }
            )

    if not rows:
        raise ValueError(f"No Meep flux rows found in {path}")

    return rows


def discover_phi_files(input_dir: Path) -> list[tuple[int, Path]]:
    discovered: list[tuple[int, Path]] = []

    for path in sorted(input_dir.glob(f"theta_{THETA:03d}_phi_*.txt")):
        match = PHI_RE.match(path.name)
        if match:
            discovered.append((int(match.group("phi")), path))

    return discovered


def finite_range(values: list[float]) -> str:
    finite = [value for value in values if math.isfinite(value)]
    if not finite:
        return "no finite values"
    return f"{min(finite):.6g} to {max(finite):.6g}"


def normalized_series(
    reference: list[dict[str, float]],
    case: list[dict[str, float]],
    key: str,
    sign: float = 1.0,
) -> tuple[list[float], list[float]]:
    if len(reference) != len(case):
        raise ValueError(
            f"Reference has {len(reference)} rows but case has {len(case)} rows"
        )

    wavelengths = [row["wavelength_um"] for row in case]
    values: list[float] = []

    for ref_row, case_row in zip(reference, case):
        inc_tran = ref_row["transmitted"]
        if inc_tran == 0:
            values.append(float("nan"))
        else:
            values.append(sign * case_row[key] / abs(inc_tran))

    return wavelengths, values


def raw_series(
    case: list[dict[str, float]], key: str
) -> tuple[list[float], list[float]]:
    return (
        [row["wavelength_um"] for row in case],
        [row[key] for row in case],
    )


def plot(input_dir: Path, output: Path, mode: str, use_list: bool) -> None:
    reference_path = input_dir / "reference_incident.txt"
    reference = parse_flux_file(reference_path) if reference_path.exists() else None

    if mode == "normalized" and reference is None:
        print("Warning: reference_incident.txt not found.")
        print("Falling back to raw transmitted/reflected flux.")
        mode = "raw"

    loaded: list[tuple[int, list[dict[str, float]]]] = []
    missing: list[Path] = []

    if use_list:
        phi_files = [
            (phi, input_dir / f"theta_{THETA:03d}_phi_{phi:03d}.txt")
            for phi in PHI_VALUES
        ]
    else:
        phi_files = discover_phi_files(input_dir)

    for phi, path in phi_files:
        if path.exists():
            loaded.append((phi, parse_flux_file(path)))
        else:
            missing.append(path)

    if not loaded:
        found = sorted(input_dir.glob("theta_*_phi_*.txt"))
        found_text = "\n".join(f"  {path}" for path in found[:20])
        raise RuntimeError(
            f"No listed phi files were found in {input_dir}.\n"
            f"First theta/phi files found:\n{found_text or '  none'}"
        )

    print(f"Plot mode: {mode}")
    if reference is not None:
        ref_tran = [row["transmitted"] for row in reference]
        print(f"reference_incident.txt: {len(reference)} rows, tran {finite_range(ref_tran)}")

    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

    for phi, case in loaded:
        label = f"phi={phi} deg"
        raw_tran = [row["transmitted"] for row in case]
        raw_refl = [row["reflected"] for row in case]
        print(
            f"theta_{THETA:03d}_phi_{phi:03d}.txt: "
            f"{len(case)} rows, tran {finite_range(raw_tran)}, refl {finite_range(raw_refl)}"
        )

        if mode == "raw":
            wavelengths, transmitted = raw_series(case, "transmitted")
            axes[0].plot(wavelengths, transmitted, label=label)

            wavelengths, reflected = raw_series(case, "reflected")
            axes[1].plot(wavelengths, reflected, label=label)
        else:
            wavelengths, extinction = normalized_series(
                reference, case, "transmitted", sign=-1.0
            )
            axes[0].plot(wavelengths, extinction, label=label)

            wavelengths, backscatter = normalized_series(reference, case, "reflected")
            axes[1].plot(wavelengths, backscatter, label=label)

    if mode == "raw":
        axes[0].set_title(f"Theta {THETA} deg, raw available phi results")
        axes[0].set_ylabel("raw transmitted flux")
        axes[1].set_ylabel("raw reflected flux")
    else:
        axes[0].set_title(f"Theta {THETA} deg, available phi results")
        axes[0].set_ylabel("-tran / |inc_tran|")
        axes[1].set_ylabel("refl / |inc_tran|")

    axes[0].grid(True, alpha=0.3)
    axes[0].legend(ncols=2, fontsize=9)
    axes[1].set_xlabel("Wavelength (um)")
    axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(output, dpi=200)
    print(f"Wrote {output}")

    if missing:
        print("Listed but not found yet:")
        for path in missing:
            print(f"  {path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--output",
        default="wu_fast_partial_results.png",
        type=Path,
        help="Output image path",
    )
    parser.add_argument(
        "--input-dir",
        default=Path("."),
        type=Path,
        help="Directory containing reference_incident.txt and theta/phi files",
    )
    parser.add_argument(
        "--mode",
        choices=["raw", "normalized"],
        default="raw",
        help="Plot raw fluxes, or normalized proxies using reference_incident.txt",
    )
    parser.add_argument(
        "--use-list",
        action="store_true",
        help="Use PHI_VALUES instead of auto-discovering theta/phi files",
    )
    args = parser.parse_args()
    plot(args.input_dir, args.output, args.mode, args.use_list)


if __name__ == "__main__":
    main()
