#!/usr/bin/env python3
"""Quick sanity checks for Wu fast-start Meep text outputs."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def parse_flux_file(path: Path) -> list[tuple[float, float, float]]:
    rows: list[tuple[float, float, float]] = []

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
                rows.append((values[0], values[1], values[2]))

    return rows


def finite_range(values: list[float]) -> str:
    finite = [value for value in values if math.isfinite(value)]
    if not finite:
        return "no finite values"
    return f"{min(finite):.6g} to {max(finite):.6g}"


def has_bad_scale(values: list[float]) -> bool:
    finite = [abs(value) for value in values if math.isfinite(value)]
    if not finite:
        return True
    vmax = max(finite)
    return vmax < 1e-40 or vmax > 1e40


def check_file(path: Path, is_reference: bool) -> bool:
    rows = parse_flux_file(path)
    if not rows:
        print(f"BAD {path}: no flux rows found")
        return False

    reflected = [row[1] for row in rows]
    transmitted = [row[2] for row in rows]
    print(
        f"{path}: {len(rows)} rows, "
        f"refl {finite_range(reflected)}, tran {finite_range(transmitted)}"
    )

    ok = True
    if has_bad_scale(transmitted):
        print(f"BAD {path}: transmitted flux scale is suspicious")
        ok = False

    if not is_reference and has_bad_scale(reflected):
        print(f"BAD {path}: reflected flux scale is suspicious")
        ok = False

    return ok


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "files",
        nargs="*",
        type=Path,
        help="Files to check. Defaults to reference_incident.txt and theta_080_phi_*.txt",
    )
    args = parser.parse_args()

    files = args.files
    if not files:
        files = [Path("reference_incident.txt")]
        files.extend(sorted(Path(".").glob("theta_080_phi_*.txt")))

    ok = True
    for path in files:
        ok = check_file(path, path.name == "reference_incident.txt") and ok

    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
