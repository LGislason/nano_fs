#!/usr/bin/env python3
"""Analyze theta-fixed phi sweeps for Wu-style modal angle dependence."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path


CASE_RE = re.compile(r"theta_(?P<theta>\d+)_phi_(?P<phi>\d+)\.txt$")


def parse_flux_file(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []

    with path.open(newline="") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line.startswith("flux") or ":" not in line:
                continue

            values: list[float] = []
            for value in next(csv.reader([line.split(":", 1)[1]])):
                try:
                    values.append(float(value.strip()))
                except ValueError:
                    continue

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


def discover_cases(input_dir: Path, theta: int) -> list[tuple[int, Path]]:
    cases: list[tuple[int, Path]] = []
    for path in sorted(input_dir.glob(f"theta_{theta:03d}_phi_*.txt")):
        match = CASE_RE.match(path.name)
        if match:
            cases.append((int(match.group("phi")), path))
    if not cases:
        raise RuntimeError(f"No theta_{theta:03d}_phi_*.txt files found in {input_dir}")
    return cases


def nearest_value(
    wavelengths: list[float], values: list[float], target_um: float
) -> tuple[float, float]:
    index = min(range(len(wavelengths)), key=lambda i: abs(wavelengths[i] - target_um))
    return wavelengths[index], values[index]


def peak_in_window(
    wavelengths: list[float],
    values: list[float],
    center_um: float,
    half_width_um: float,
) -> tuple[float, float]:
    candidates = [
        i for i, wavelength in enumerate(wavelengths)
        if abs(wavelength - center_um) <= half_width_um
    ]
    if not candidates:
        return nearest_value(wavelengths, values, center_um)
    index = max(candidates, key=lambda i: values[i])
    return wavelengths[index], values[index]


def solve_3x3(a: list[list[float]], b: list[float]) -> list[float]:
    matrix = [row[:] + [rhs] for row, rhs in zip(a, b)]
    n = 3
    for col in range(n):
        pivot = max(range(col, n), key=lambda row: abs(matrix[row][col]))
        matrix[col], matrix[pivot] = matrix[pivot], matrix[col]
        if abs(matrix[col][col]) < 1e-30:
            raise ValueError("Singular fit matrix")
        scale = matrix[col][col]
        matrix[col] = [value / scale for value in matrix[col]]
        for row in range(n):
            if row == col:
                continue
            factor = matrix[row][col]
            matrix[row] = [
                current - factor * pivot_value
                for current, pivot_value in zip(matrix[row], matrix[col])
            ]
    return [matrix[row][n] for row in range(n)]


def fit_cos2_family(phis: list[int], values: list[float]) -> tuple[float, float, float]:
    """Fit y = c + a cos(2 phi) + b sin(2 phi).

    This is the general form of an offset cosine-squared response. The returned
    R2 is a compact quality check for Wu-style angle dependence.
    """

    normal = [[0.0 for _ in range(3)] for _ in range(3)]
    rhs = [0.0 for _ in range(3)]

    for phi, value in zip(phis, values):
        angle = math.radians(2.0 * phi)
        row = [1.0, math.cos(angle), math.sin(angle)]
        for i in range(3):
            rhs[i] += row[i] * value
            for j in range(3):
                normal[i][j] += row[i] * row[j]

    c, a, b = solve_3x3(normal, rhs)
    predictions = []
    for phi in phis:
        angle = math.radians(2.0 * phi)
        predictions.append(c + a * math.cos(angle) + b * math.sin(angle))

    mean = sum(values) / len(values)
    ss_res = sum((value - prediction) ** 2 for value, prediction in zip(values, predictions))
    ss_tot = sum((value - mean) ** 2 for value in values)
    r2 = 1.0 - ss_res / ss_tot if ss_tot else float("nan")
    amplitude = math.hypot(a, b)
    return c, amplitude, r2


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--theta", type=int, default=80)
    parser.add_argument(
        "--quantity",
        choices=["backscatter", "extinction"],
        default="backscatter",
    )
    parser.add_argument(
        "--modes",
        type=float,
        nargs="+",
        default=[0.672, 0.853],
        help="Mode probe wavelengths in microns.",
    )
    parser.add_argument(
        "--window",
        type=float,
        default=0.025,
        help="Half-width around each mode wavelength for local peak extraction.",
    )
    args = parser.parse_args()

    reference = parse_flux_file(args.input_dir / "reference_incident.txt")
    rows: list[dict[str, float]] = []

    for phi, path in discover_cases(args.input_dir, args.theta):
        wavelengths, values = normalized_proxy(
            reference, parse_flux_file(path), args.quantity
        )
        global_peak_index = max(range(len(values)), key=lambda index: values[index])
        row = {
            "phi": float(phi),
            "peak_wavelength_um": wavelengths[global_peak_index],
            "peak_value": values[global_peak_index],
        }
        for mode in args.modes:
            mode_wavelength, mode_value = peak_in_window(
                wavelengths, values, mode, args.window
            )
            row[f"mode_{mode:.3f}_wavelength_um"] = mode_wavelength
            row[f"mode_{mode:.3f}_value"] = mode_value
        rows.append(row)

    print(f"Input: {args.input_dir}")
    print(f"Theta: {args.theta} deg")
    print(f"Quantity: {args.quantity}")
    print()

    headers = ["phi", "peak_wavelength_um", "peak_value"]
    for mode in args.modes:
        headers.extend([f"mode_{mode:.3f}_wavelength_um", f"mode_{mode:.3f}_value"])
    print(",".join(headers))
    for row in rows:
        print(",".join(f"{row[header]:.6g}" for header in headers))

    print()
    phis = [int(row["phi"]) for row in rows]
    for mode in args.modes:
        key = f"mode_{mode:.3f}_value"
        values = [row[key] for row in rows]
        offset, amplitude, r2 = fit_cos2_family(phis, values)
        print(
            f"mode {mode:.3f} um: min={min(values):.6g}, max={max(values):.6g}, "
            f"offset={offset:.6g}, cos2_amplitude={amplitude:.6g}, R2={r2:.3f}"
        )


if __name__ == "__main__":
    main()
