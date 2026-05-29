#!/usr/bin/env python3
"""Plot old-Meep HDF5 field slices from nanorod_wu_field_snapshot.ctl."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SLICES = [
    ("xy_rod", "XY at rod midplane", "x (um)", "y (um)"),
    ("xy_above", "XY above rods", "x (um)", "y (um)"),
    ("xz", "XZ at y=0", "x (um)", "z (um)"),
]


def first_dataset(group: h5py.Group) -> np.ndarray:
    datasets: list[np.ndarray] = []

    def visitor(_name: str, obj: h5py.Dataset) -> None:
        if isinstance(obj, h5py.Dataset):
            datasets.append(np.asarray(obj))

    group.visititems(visitor)
    if not datasets:
        raise ValueError("No datasets found")

    # Prefer multidimensional field arrays over scalar metadata.
    datasets.sort(key=lambda data: data.ndim, reverse=True)
    return datasets[0]


def read_h5_array(path: Path) -> np.ndarray:
    with h5py.File(path, "r") as handle:
        return np.squeeze(first_dataset(handle))


def find_h5(input_dir: Path, stem: str) -> Path:
    matches = sorted(input_dir.glob(f"{stem}*.h5"))
    if not matches:
        raise FileNotFoundError(f"No HDF5 file matching {stem}*.h5 in {input_dir}")
    return matches[0]


def load_slice(input_dir: Path, prefix: str) -> tuple[np.ndarray, np.ndarray]:
    eps = read_h5_array(find_h5(input_dir, f"{prefix}_eps"))
    ex = read_h5_array(find_h5(input_dir, f"{prefix}_ex"))
    ey = read_h5_array(find_h5(input_dir, f"{prefix}_ey"))
    ez = read_h5_array(find_h5(input_dir, f"{prefix}_ez"))

    # Old Meep output snapshots are real time-domain fields. If complex arrays
    # are encountered, abs handles them correctly.
    e2 = np.abs(ex) ** 2 + np.abs(ey) ** 2 + np.abs(ez) ** 2
    return e2, eps


def edge_median(field: np.ndarray) -> float:
    edges = np.concatenate([field[0, :], field[-1, :], field[:, 0], field[:, -1]])
    median = float(np.median(edges))
    if median > 0:
        return median
    positive = field[field > 0]
    return float(np.median(positive)) if positive.size else 1.0


def normalize(field: np.ndarray, contrast: str, percentile: float) -> np.ndarray:
    if contrast == "relative":
        field = np.abs(field / edge_median(field) - 1.0)
    vmax = float(np.percentile(field, percentile))
    if vmax <= 0:
        vmax = float(field.max()) or 1.0
    return np.clip(field / vmax, 0, 1)


def add_contours(ax: plt.Axes, eps: np.ndarray) -> None:
    eps_plot = np.rot90(np.real(eps))
    levels = [level for level in (1.5, 3.0, 5.0) if eps_plot.min() < level < eps_plot.max()]
    if levels:
        ax.contour(eps_plot, levels=levels, colors="white", linewidths=0.7)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("-o", "--output", type=Path, default=Path("field_snapshot.png"))
    parser.add_argument(
        "--contrast",
        choices=["total", "relative"],
        default="relative",
        help="Plot total |E|^2 or relative contrast against panel edge median.",
    )
    parser.add_argument("--vmax-percentile", type=float, default=99.5)
    args = parser.parse_args()

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
    im = None
    for ax, (prefix, title, xlabel, ylabel) in zip(axes, SLICES):
        field, eps = load_slice(args.input_dir, prefix)
        field = normalize(field, args.contrast, args.vmax_percentile)
        im = ax.imshow(np.rot90(field), origin="lower", cmap="inferno", vmin=0, vmax=1)
        add_contours(ax, eps)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

    fig.colorbar(im, ax=axes, shrink=0.88, label=f"{args.contrast} field, p{args.vmax_percentile:g}=1")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=250)
    plt.close(fig)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
