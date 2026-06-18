#!/usr/bin/env python3
"""Plot old-Meep HDF5 field slices from nanorod_wu_field_snapshot.ctl."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


SLICES = [
    ("xy_rod", "Rod midplane (z = 0)", "x (nm)", "y (nm)"),
    ("xy_above", "30 nm above rods", "x (nm)", "y (nm)"),
    ("xz", "Vertical cut (y = 0)", "x (nm)", "z (nm)"),
]

# These match the output volumes in nanorod_wu_field_snapshot.ctl:
# sx/sy/sz minus the two PML layers and a small margin, centered at zero.
FIELD_SX = 0.52
FIELD_SY = 0.52
FIELD_SZ = 0.92
ROD_LENGTH = 0.069
ROD_WIDTH = 0.024
NM = 1000.0  # um -> nm for axis display

# Publication-style defaults (kept local so the script stays self-contained).
PLOT_STYLE = {
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "axes.linewidth": 0.8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
}


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
        matches = sorted(input_dir.glob(f"*{stem}*.h5"))
    if not matches:
        raise FileNotFoundError(f"No HDF5 file matching {stem}*.h5 in {input_dir}")
    return matches[0]


def read_h5_optional(input_dir: Path, stem: str) -> np.ndarray | None:
    """Read a slice if a matching *stem*.h5 exists, else return None."""
    matches = sorted(input_dir.glob(f"*{stem}*.h5"))
    return read_h5_array(matches[0]) if matches else None


def load_slice(input_dir: Path, prefix: str) -> tuple[np.ndarray, np.ndarray]:
    eps = read_h5_array(find_h5(input_dir, f"{prefix}_eps"))
    ex = read_h5_array(find_h5(input_dir, f"{prefix}_ex"))
    ey = read_h5_array(find_h5(input_dir, f"{prefix}_ey"))
    ez = read_h5_array(find_h5(input_dir, f"{prefix}_ez"))

    # If quarter-period companion snapshots ("_q") exist, combine them to recover
    # the phase-independent steady-state amplitude:  |E_amp|^2 = E(t)^2 +
    # E(t+T/4)^2  (since cos^2 + sin^2 = 1, per component).  This is the correct
    # quantity for an |E|^2/|E0|^2 enhancement map.  Otherwise fall back to a
    # single instantaneous snapshot, which is phase-dependent and only
    # qualitatively correct.
    ex_q = read_h5_optional(input_dir, f"{prefix}_ex_q")
    ey_q = read_h5_optional(input_dir, f"{prefix}_ey_q")
    ez_q = read_h5_optional(input_dir, f"{prefix}_ez_q")

    # Old Meep output snapshots are real time-domain fields. If complex arrays
    # are encountered, abs handles them correctly.
    if ex_q is not None and ey_q is not None and ez_q is not None:
        e2 = (
            np.abs(ex) ** 2 + np.abs(ex_q) ** 2
            + np.abs(ey) ** 2 + np.abs(ey_q) ** 2
            + np.abs(ez) ** 2 + np.abs(ez_q) ** 2
        )
        print(f"[{prefix}] using quarter-period amplitude (true |E|^2)")
    else:
        e2 = np.abs(ex) ** 2 + np.abs(ey) ** 2 + np.abs(ez) ** 2
        print(f"[{prefix}] single instantaneous snapshot (phase-dependent)")
    return e2, eps


def run_params(input_dir: Path) -> dict[str, float]:
    params = {"theta": 80.0, "phi": 0.0, "gap": 0.006}
    log_path = input_dir / "run.log"
    if not log_path.exists():
        return params
    cli_pattern = re.compile(r"command-line param: (theta|phi|gap)=([0-9.]+)")
    tip_gap_pattern = re.compile(r"rod-tip-gap=([0-9.eE+-]+)")
    for line in log_path.read_text(errors="replace").splitlines():
        match = cli_pattern.search(line)
        if match:
            params[match.group(1)] = float(match.group(2))
        match = tip_gap_pattern.search(line)
        if match:
            params["rod_tip_gap"] = float(match.group(1))
    return params


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


def image_data(field: np.ndarray) -> np.ndarray:
    return np.real(field).T


def image_extent(prefix: str) -> tuple[float, float, float, float]:
    if prefix.startswith("xy_"):
        return (-FIELD_SX / 2 * NM, FIELD_SX / 2 * NM, -FIELD_SY / 2 * NM, FIELD_SY / 2 * NM)
    if prefix == "xz":
        return (-FIELD_SX / 2 * NM, FIELD_SX / 2 * NM, -FIELD_SZ / 2 * NM, FIELD_SZ / 2 * NM)
    raise ValueError(f"Unknown slice prefix: {prefix}")


def rotate_xy(x: float, y: float, angle: float) -> tuple[float, float]:
    return x * np.cos(angle) - y * np.sin(angle), x * np.sin(angle) + y * np.cos(angle)


def capsule_outline(center: np.ndarray, angle: float) -> tuple[np.ndarray, np.ndarray]:
    radius = ROD_WIDTH / 2
    body_length = ROD_LENGTH - 2 * radius
    u = np.array([np.cos(angle), np.sin(angle)])
    front = center + 0.5 * body_length * u
    back = center - 0.5 * body_length * u
    front_angles = np.linspace(angle - np.pi / 2, angle + np.pi / 2, 48)
    back_angles = np.linspace(angle + np.pi / 2, angle + 3 * np.pi / 2, 48)
    front_arc = front[:, None] + radius * np.vstack((np.cos(front_angles), np.sin(front_angles)))
    back_arc = back[:, None] + radius * np.vstack((np.cos(back_angles), np.sin(back_angles)))
    points = np.hstack((front_arc, back_arc, front_arc[:, :1]))
    return points[0], points[1]


def rod_geometry(params: dict[str, float]) -> tuple[list[tuple[np.ndarray, float]], tuple[np.ndarray, np.ndarray]]:
    theta = np.deg2rad(params["theta"])
    phi = np.deg2rad(params["phi"])
    radius = ROD_WIDTH / 2
    tip_gap = params.get("rod_tip_gap", params["gap"] + 2 * radius * (1 - np.sin(theta / 2)))

    rod1_ang0 = np.pi
    rod2_ang0 = np.pi - theta
    gap_axis = 0.5 * (rod1_ang0 + rod2_ang0) - np.pi / 2
    gap_axis_vec = np.array([np.cos(gap_axis), np.sin(gap_axis)])
    rod1_tip = -0.5 * tip_gap * gap_axis_vec
    rod2_tip = 0.5 * tip_gap * gap_axis_vec
    rod1_center0 = rod1_tip + 0.5 * ROD_LENGTH * np.array([np.cos(rod1_ang0), np.sin(rod1_ang0)])
    rod2_center0 = rod2_tip + 0.5 * ROD_LENGTH * np.array([np.cos(rod2_ang0), np.sin(rod2_ang0)])
    shift = np.array([0.0, 0.0])

    rods = []
    for center0, angle0 in ((rod1_center0, rod1_ang0), (rod2_center0, rod2_ang0)):
        cx, cy = rotate_xy(*(center0 + shift), phi)
        rods.append((np.array([cx, cy]), angle0 + phi))

    tip1 = np.array(rotate_xy(*(rod1_tip + shift), phi))
    tip2 = np.array(rotate_xy(*(rod2_tip + shift), phi))
    return rods, (tip1, tip2)


def add_analytic_rods(ax: plt.Axes, params: dict[str, float], linewidth: float = 0.9) -> None:
    rods, (tip1, tip2) = rod_geometry(params)
    for center, angle in rods:
        x, y = capsule_outline(center, angle)
        ax.plot(x * NM, y * NM, color="white", linewidth=linewidth, solid_capstyle="round")
    ax.plot(
        [tip1[0] * NM, tip2[0] * NM],
        [tip1[1] * NM, tip2[1] * NM],
        color="white",
        linewidth=max(linewidth, 1.0),
        solid_capstyle="round",
    )


def add_contours(ax: plt.Axes, eps: np.ndarray, prefix: str, linewidth: float = 0.45) -> None:
    eps_plot = image_data(eps)
    levels = [level for level in (1.5, 3.0, 5.0) if eps_plot.min() < level < eps_plot.max()]
    if levels:
        xmin, xmax, ymin, ymax = image_extent(prefix)
        x = np.linspace(xmin, xmax, eps_plot.shape[1])
        y = np.linspace(ymin, ymax, eps_plot.shape[0])
        ax.contour(
            x,
            y,
            eps_plot,
            levels=levels,
            colors="white",
            linewidths=linewidth,
        )


def plot_slice(
    ax: plt.Axes,
    prefix: str,
    title: str,
    xlabel: str,
    ylabel: str,
    field: np.ndarray,
    eps: np.ndarray,
    params: dict[str, float],
) -> plt.AxesImage:
    extent = image_extent(prefix)
    im = ax.imshow(
        image_data(field),
        origin="lower",
        extent=extent,
        cmap="inferno",
        vmin=0,
        vmax=1,
        interpolation="nearest",
    )
    if prefix == "xy_rod":
        add_analytic_rods(ax, params)
    else:
        add_contours(ax, eps, prefix)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_aspect("equal")
    return im


def add_gap_zoom(
    ax: plt.Axes,
    field: np.ndarray,
    eps: np.ndarray,
    params: dict[str, float],
    zoom_half_width: float,
) -> plt.AxesImage:
    extent = image_extent("xy_rod")
    im = ax.imshow(
        image_data(field),
        origin="lower",
        extent=extent,
        cmap="inferno",
        vmin=0,
        vmax=1,
        interpolation="nearest",
    )
    add_analytic_rods(ax, params, linewidth=1.0)
    ax.set_xlim(-zoom_half_width * NM, zoom_half_width * NM)
    ax.set_ylim(-zoom_half_width * NM, zoom_half_width * NM)
    ax.set_title("Gap zoom")
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.set_aspect("equal")
    return im


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
    parser.add_argument(
        "--zoom-half-width",
        type=float,
        default=0.065,
        help="Half-width in um for the rod-midplane gap zoom panel.",
    )
    args = parser.parse_args()

    plt.rcParams.update(PLOT_STYLE)
    params = run_params(args.input_dir)

    # Layout: top row keeps the rod midplane and ITS gap zoom side by side
    # (the zoomed-out view next to the zoom); bottom row holds the other cuts.
    #   [ rod midplane ][ gap zoom ]
    #   [ 30nm above   ][  y=0 cut ]
    fig, axes = plt.subplots(2, 2, figsize=(10.6, 9.4), constrained_layout=True)
    ax_full, ax_zoom = axes[0, 0], axes[0, 1]
    ax_above, ax_xz = axes[1, 0], axes[1, 1]
    panel_axis = {"xy_rod": ax_full, "xy_above": ax_above, "xz": ax_xz}

    slices: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    im = None
    for prefix, title, xlabel, ylabel in SLICES:
        field, eps = load_slice(args.input_dir, prefix)
        field = normalize(field, args.contrast, args.vmax_percentile)
        slices[prefix] = (field, eps)
        im = plot_slice(panel_axis[prefix], prefix, title, xlabel, ylabel, field, eps, params)

    xy_field, xy_eps = slices["xy_rod"]
    im = add_gap_zoom(ax_zoom, xy_field, xy_eps, params, args.zoom_half_width)
    # Mark the zoomed region on the full rod-midplane panel (units are nm).
    zb = args.zoom_half_width * NM
    ax_full.add_patch(Rectangle((-zb, -zb), 2 * zb, 2 * zb, fill=False,
                                edgecolor="white", linewidth=0.8, linestyle="--"))

    fig.suptitle(args.input_dir.name, fontsize=12)
    label = ("relative near-field contrast" if args.contrast == "relative"
             else r"$|E|^2$ (a.u.)") + f"  (p{args.vmax_percentile:g} = 1)"
    fig.colorbar(im, ax=axes, shrink=0.85, pad=0.02, label=label)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output)
    plt.close(fig)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
