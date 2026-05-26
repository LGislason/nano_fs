#!/usr/bin/env python3
"""
Plot available Wu fast-start Meep results.

Copy or run this script from the directory containing:

  results/reference_incident.txt
  results/theta_080_phi_000.txt
  results/theta_080_phi_020.txt
  ...

Uncomment additional phi values in PHI_VALUES as those runs finish.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


DEFAULT_RESULTS_DIR = Path(".")
THETA = 80

# Completed so far, based on your current run status.
# Uncomment the remaining values as their files appear in results/.
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
    """Return rows from Meep display-fluxes output.

    Meep writes CSV-like lines beginning with "flux1:" or similar. For this
    control file the columns are expected to be:

      frequency, reflected_flux, transmitted_flux
    """
    rows: list[dict[str, float]] = []

    with path.open(newline="") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line.startswith("flux"):
                continue

            _, values_text = line.split(":", 1)
            values = [
                float(value.strip())
                for value in next(csv.reader([values_text]))
                if value.strip()
            ]

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


def load_case(results_dir: Path, phi: int) -> list[dict[str, float]]:
    path =  f"theta_{THETA:03d}_phi_{phi:03d}.txt"
    if not path.exists():
        raise FileNotFoundError(path)
    return parse_flux_file(path)


def extinction_proxy(
    reference: list[dict[str, float]], case: list[dict[str, float]]
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
            # See README: -tran / inc_tran is an extinction-like proxy.
            values.append(-case_row["transmitted"] / abs(inc_tran))

    return wavelengths, values


def backscatter_proxy(
    reference: list[dict[str, float]], case: list[dict[str, float]]
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
            # See README: refl / inc_tran is a backscatter proxy.
            values.append(case_row["reflected"] / abs(inc_tran))

    return wavelengths, values


def plot_available(results_dir: Path, output: Path) -> None:
    reference_file = results_dir / "reference_incident.txt"
    reference = None
    location_name = "current directory" if results_dir == Path(".") else str(results_dir)

    if reference_file.exists():
        reference = parse_flux_file(reference_file)
    else:
        print(f"Warning: missing {reference_file}")
        print("Plotting raw case fluxes instead of normalized proxies.")

    loaded: list[tuple[int, list[dict[str, float]]]] = []
    missing: list[Path] = []

    for phi in PHI_VALUES:
        path = results_dir / f"theta_{THETA:03d}_phi_{phi:03d}.txt"
        if path.exists():
            loaded.append((phi, parse_flux_file(path)))
        else:
            missing.append(path)

    if not loaded:
        found = sorted(results_dir.glob("*.txt"))
        found_text = "\n".join(f"  {path}" for path in found[:20])
        raise RuntimeError(
            f"No listed phi files were found in {location_name}.\n"
            f"First .txt files found there:\n{found_text or '  none'}"
        )

    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

    for phi, case in loaded:
        if reference is None:
            wavelengths = [row["wavelength_um"] for row in case]
            axes[0].plot(
                wavelengths,
                [row["transmitted"] for row in case],
                label=f"phi={phi} deg",
            )
            axes[1].plot(
                wavelengths,
                [row["reflected"] for row in case],
                label=f"phi={phi} deg",
            )
        else:
            wavelengths, extinction = extinction_proxy(reference, case)
            axes[0].plot(wavelengths, extinction, label=f"phi={phi} deg")

            wavelengths, backscatter = backscatter_proxy(reference, case)
            axes[1].plot(wavelengths, backscatter, label=f"phi={phi} deg")

    if reference is None:
        axes[0].set_ylabel("raw transmitted flux")
        axes[1].set_ylabel("raw reflected flux")
        axes[0].set_title(f"Theta {THETA} deg, raw available phi results")
    else:
        axes[0].set_ylabel("-tran / |inc_tran|")
        axes[1].set_ylabel("refl / |inc_tran|")
        axes[0].set_title(f"Theta {THETA} deg, available phi results")

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
        "--results-dir",
        default=DEFAULT_RESULTS_DIR,
        type=Path,
        help="Directory containing reference_incident.txt and theta/phi result files",
    )
    args = parser.parse_args()

    plot_available(args.results_dir, args.output)


if __name__ == "__main__":
    main()
