#!/usr/bin/env python3
"""Plot a quick Meep dielectric map of the nanorod dimer geometry."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def import_meep():
    try:
        import meep as mp
    except ImportError as exc:
        raise SystemExit(
            "This script requires the Python Meep package. Run it in the same "
            "environment where `python -c 'import meep'` works."
        ) from exc
    return mp


def rotate_xy(x: float, y: float, angle: float) -> tuple[float, float]:
    return (
        x * math.cos(angle) - y * math.sin(angle),
        x * math.sin(angle) + y * math.cos(angle),
    )


def nanorod_centers_and_axes(
    theta_deg: float,
    phi_deg: float,
    gap: float,
    rod_length: float,
    rod_width: float,
) -> list[tuple[float, float, float, float]]:
    theta = math.radians(theta_deg)
    phi = math.radians(phi_deg)
    rod_radius = rod_width / 2.0
    rod_tip_gap = gap + 2.0 * rod_radius * (1.0 - math.sin(theta / 2.0))

    # theta is the internal opening angle between the two rod directions
    # pointing away from the gap.
    rod1_out_angle = math.pi
    rod2_out_angle = math.pi - theta
    gap_axis_angle = 0.5 * (rod1_out_angle + rod2_out_angle) - math.pi / 2.0
    gap_axis_x = math.cos(gap_axis_angle)
    gap_axis_y = math.sin(gap_axis_angle)

    rod1_tip_x = -0.5 * rod_tip_gap * gap_axis_x
    rod1_tip_y = -0.5 * rod_tip_gap * gap_axis_y
    rod2_tip_x = 0.5 * rod_tip_gap * gap_axis_x
    rod2_tip_y = 0.5 * rod_tip_gap * gap_axis_y

    rod1_bx = rod1_tip_x + 0.5 * rod_length * math.cos(rod1_out_angle)
    rod1_by = rod1_tip_y + 0.5 * rod_length * math.sin(rod1_out_angle)
    rod2_bx = rod2_tip_x + 0.5 * rod_length * math.cos(rod2_out_angle)
    rod2_by = rod2_tip_y + 0.5 * rod_length * math.sin(rod2_out_angle)

    dimer_shift_x = 0.0
    dimer_shift_y = 0.0

    rod1_cx0 = rod1_bx + dimer_shift_x
    rod1_cy0 = rod1_by + dimer_shift_y
    rod2_cx0 = rod2_bx + dimer_shift_x
    rod2_cy0 = rod2_by + dimer_shift_y

    rod1_cx, rod1_cy = rotate_xy(rod1_cx0, rod1_cy0, phi)
    rod2_cx, rod2_cy = rotate_xy(rod2_cx0, rod2_cy0, phi)

    rod1_angle = phi + rod1_out_angle
    rod2_angle = phi + rod2_out_angle

    return [
        (rod1_cx, rod1_cy, math.cos(rod1_angle), math.sin(rod1_angle)),
        (rod2_cx, rod2_cy, math.cos(rod2_angle), math.sin(rod2_angle)),
    ]


def build_geometry(args, mp):
    rod_radius = args.rod_width / 2.0
    rod_body_length = args.rod_length - 2.0 * rod_radius
    rod_material = mp.Medium(epsilon=args.rod_epsilon)

    geometry = []
    for cx, cy, ux, uy in nanorod_centers_and_axes(
        args.theta, args.phi, args.gap, args.rod_length, args.rod_width
    ):
        axis = mp.Vector3(ux, uy, 0)
        center = mp.Vector3(cx, cy, 0)
        geometry.append(
            mp.Cylinder(
                radius=rod_radius,
                height=rod_body_length,
                center=center,
                axis=axis,
                material=rod_material,
            )
        )
        geometry.append(
            mp.Sphere(
                radius=rod_radius,
                center=mp.Vector3(
                    cx + 0.5 * rod_body_length * ux,
                    cy + 0.5 * rod_body_length * uy,
                    0,
                ),
                material=rod_material,
            )
        )
        geometry.append(
            mp.Sphere(
                radius=rod_radius,
                center=mp.Vector3(
                    cx - 0.5 * rod_body_length * ux,
                    cy - 0.5 * rod_body_length * uy,
                    0,
                ),
                material=rod_material,
            )
        )

    return geometry


def plot_geometry(args) -> None:
    mp = import_meep()

    cell_size = mp.Vector3(args.sx, args.sy, args.sz)
    sim = mp.Simulation(
        cell_size=cell_size,
        boundary_layers=[mp.PML(args.dpml)],
        geometry=build_geometry(args, mp),
        default_material=mp.Medium(index=args.nbg),
        resolution=args.res,
    )

    sim.init_sim()
    eps = sim.get_array(
        center=mp.Vector3(0, 0, 0),
        size=mp.Vector3(args.sx, args.sy, 0),
        component=mp.Dielectric,
    )

    eps = np.rot90(eps)
    extent = [
        -args.sx / 2.0,
        args.sx / 2.0,
        -args.sy / 2.0,
        args.sy / 2.0,
    ]

    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    image = ax.imshow(eps, extent=extent, origin="lower", cmap="viridis")
    fig.colorbar(image, ax=ax, label="dielectric epsilon")

    pml_box = plt.Rectangle(
        (-args.sx / 2.0 + args.dpml, -args.sy / 2.0 + args.dpml),
        args.sx - 2.0 * args.dpml,
        args.sy - 2.0 * args.dpml,
        fill=False,
        linestyle="--",
        linewidth=1.0,
        color="white",
        alpha=0.9,
        label="non-PML region",
    )
    ax.add_patch(pml_box)

    ax.set_aspect("equal")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title(
        f"Nanorod dimer geometry: theta={args.theta:g} deg, "
        f"phi={args.phi:g} deg, gap={args.gap * 1000:g} nm"
    )
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(args.output, dpi=250)
    print(f"Wrote {args.output}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--theta", type=float, default=80.0)
    parser.add_argument("--phi", type=float, default=0.0)
    parser.add_argument("--gap", type=float, default=0.010, help="Gap in um")
    parser.add_argument("--res", type=int, default=400)
    parser.add_argument("--nbg", type=float, default=1.28)
    parser.add_argument("--sx", type=float, default=0.80)
    parser.add_argument("--sy", type=float, default=0.80)
    parser.add_argument("--sz", type=float, default=1.20)
    parser.add_argument("--dpml", type=float, default=0.12)
    parser.add_argument("--rod-length", type=float, default=0.069)
    parser.add_argument("--rod-width", type=float, default=0.024)
    parser.add_argument(
        "--rod-epsilon",
        type=float,
        default=9.0,
        help="Simple display material epsilon for the rods.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("nanodimer_geometry_theta080_phi000.png"),
    )
    args = parser.parse_args()
    plot_geometry(args)


if __name__ == "__main__":
    main()
