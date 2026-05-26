#!/usr/bin/env python3
"""Create a vertically stacked spectrum plot from Wu fast-start result files."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt


PHI_RE = re.compile(r"theta_(?P<theta>\d+)_phi_(?P<phi>\d+)\.txt$")


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
                try:
                    values.append(float(value.strip()))
                except ValueError:
                    pass

            if len(values) >= 3:
                frequency, reflected, transmitted = values[:3]
                rows.append(
                    {
                        "wavelength_um": 1.0 / frequency,
                        "reflected": reflected,
                        "transmitted": transmitted,
                    }
                )

    if not rows:
        raise ValueError(f"No flux rows found in {path}")

    return rows


def discover_cases(input_dir: Path) -> tuple[str, list[tuple[int, Path]]]:
    raw_cases: list[tuple[int, int, Path]] = []

    for path in sorted(input_dir.glob("theta_*_phi_*.txt")):
        match = PHI_RE.match(path.name)
        if match:
            raw_cases.append(
                (int(match.group("theta")), int(match.group("phi")), path)
            )

    if not raw_cases:
        raise RuntimeError(f"No theta/phi result files found in {input_dir}")

    thetas = {theta for theta, _, _ in raw_cases}
    phis = {phi for _, phi, _ in raw_cases}

    if len(thetas) > 1 and len(phis) == 1:
        return "theta", sorted((theta, path) for theta, _, path in raw_cases)
    if len(phis) > 1 and len(thetas) == 1:
        return "phi", sorted((phi, path) for _, phi, path in raw_cases)

    return "case", [(index, path) for index, (_, _, path) in enumerate(raw_cases)]


def normalized_proxy(
    reference: list[dict[str, float]],
    case: list[dict[str, float]],
    quantity: str,
) -> tuple[list[float], list[float]]:
    if len(reference) != len(case):
        raise ValueError(
            f"Reference has {len(reference)} rows but case has {len(case)} rows"
        )

    wavelengths = [row["wavelength_um"] for row in case]
    values: list[float] = []

    for ref_row, case_row in zip(reference, case):
        inc = abs(ref_row["transmitted"])
        if quantity == "extinction":
            values.append(-case_row["transmitted"] / inc)
        elif quantity == "backscatter":
            values.append(case_row["reflected"] / inc)
        else:
            raise ValueError(f"Unknown quantity: {quantity}")

    return wavelengths, values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("."))
    parser.add_argument(
        "--quantity",
        choices=["extinction", "backscatter"],
        default="backscatter",
    )
    parser.add_argument(
        "--offset",
        type=float,
        default=None,
        help="Vertical offset between traces. Defaults to 1.15 times max curve height.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("stacked_spectra.png"),
    )
    args = parser.parse_args()

    reference = parse_flux_file(args.input_dir / "reference_incident.txt")
    label_kind, cases = discover_cases(args.input_dir)

    curves: list[tuple[int, list[float], list[float]]] = []
    max_height = 0.0

    for label_value, path in cases:
        wavelengths, values = normalized_proxy(
            reference, parse_flux_file(path), args.quantity
        )
        peak = max(values)
        max_height = max(max_height, peak)
        curves.append((label_value, wavelengths, values))

    offset = args.offset if args.offset is not None else max_height * 1.15
    if offset <= 0:
        offset = 1.0

    fig_height = max(5.5, 0.75 * len(curves) + 1.6)
    fig, ax = plt.subplots(figsize=(7.2, fig_height))

    for index, (label_value, wavelengths, values) in enumerate(curves):
        y_offset = index * offset
        shifted = [value + y_offset for value in values]
        ax.plot(wavelengths, shifted, color="black", linewidth=1.4)
        label = f"{label_value} deg" if label_kind in {"theta", "phi"} else str(label_value)
        ax.text(
            max(wavelengths) + 0.01,
            y_offset,
            label,
            va="center",
            fontsize=9,
        )

    ax.set_xlabel("Wavelength (um)")
    ylabel = {
        "extinction": "Extinction proxy (offset)",
        "backscatter": "Backscatter proxy (offset)",
    }[args.quantity]
    ax.set_ylabel(ylabel)
    sweep_label = {"theta": "theta", "phi": "phi", "case": "case"}[label_kind]
    ax.set_title(f"Stacked {args.quantity} proxy spectra by {sweep_label}")
    ax.set_yticks([])
    ax.grid(True, axis="x", alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(min(curves[0][1]), max(curves[0][1]) + 0.08)

    fig.tight_layout()
    fig.savefig(args.output, dpi=250)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
