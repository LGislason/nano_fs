#!/usr/bin/env python3
"""Plot a simple top-view nanorod dimer schematic with the internal theta angle."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Rectangle
from matplotlib.transforms import Affine2D


def rod_centers_axes(theta_deg: float, phi_deg: float, gap: float, rod_length: float):
    theta = math.radians(theta_deg)
    phi = math.radians(phi_deg)

    rod1_out = math.pi
    rod2_out = math.pi - theta

    tips = [(-gap / 2.0, 0.0), (gap / 2.0, 0.0)]
    angles = [rod1_out, rod2_out]

    centers = []
    axes = []
    for tip, angle in zip(tips, angles):
        cx = tip[0] + 0.5 * rod_length * math.cos(angle)
        cy = tip[1] + 0.5 * rod_length * math.sin(angle)
        centers.append((cx, cy))
        axes.append(angle)

    dimer_shift_x = -0.5 * (centers[0][0] + centers[1][0])
    dimer_shift_y = -0.5 * (centers[0][1] + centers[1][1])

    rotated = []
    for (cx, cy), angle in zip(centers, axes):
        x = cx + dimer_shift_x
        y = cy + dimer_shift_y
        xr = x * math.cos(phi) - y * math.sin(phi)
        yr = x * math.sin(phi) + y * math.cos(phi)
        rotated.append((xr, yr, angle + phi))

    return rotated


def add_rod(ax, cx, cy, angle, length, width, color):
    patch = Rectangle(
        (-length / 2.0, -width / 2.0),
        length,
        width,
        facecolor=color,
        edgecolor="black",
        linewidth=1.0,
        alpha=0.85,
    )
    transform = Affine2D().rotate(angle).translate(cx, cy) + ax.transData
    patch.set_transform(transform)
    ax.add_patch(patch)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--theta", type=float, default=80.0)
    parser.add_argument("--phi", type=float, default=0.0)
    parser.add_argument("--gap", type=float, default=0.010)
    parser.add_argument("--rod-length", type=float, default=0.069)
    parser.add_argument("--rod-width", type=float, default=0.024)
    parser.add_argument("-o", "--output", type=Path, default=Path("nanodimer_topview_theta080.png"))
    args = parser.parse_args()

    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    rods = rod_centers_axes(args.theta, args.phi, args.gap, args.rod_length)

    for cx, cy, angle in rods:
        add_rod(ax, cx, cy, angle, args.rod_length, args.rod_width, "#d8a020")

    arc_radius = 0.055
    arc = Arc(
        (0, 0),
        2 * arc_radius,
        2 * arc_radius,
        angle=math.degrees(math.radians(args.phi)),
        theta1=180 - args.theta,
        theta2=180,
        color="tab:red",
        linewidth=1.8,
    )
    ax.add_patch(arc)

    ax.scatter([-args.gap / 2.0, args.gap / 2.0], [0, 0], s=16, color="black", zorder=3)
    ax.text(0, arc_radius * 1.15, f"theta = {args.theta:g} deg", ha="center", color="tab:red")

    lim = 0.095
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title("Nanodimer top view")
    fig.tight_layout()
    fig.savefig(args.output, dpi=250)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
