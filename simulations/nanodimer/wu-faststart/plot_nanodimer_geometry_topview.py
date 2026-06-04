#!/usr/bin/env python3
"""Plot a simple top-view nanorod dimer schematic with the internal theta angle."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Arc


def rod_tip_gap(theta_deg: float, gap: float, rod_radius: float) -> float:
    theta = math.radians(theta_deg)
    return gap + 2.0 * rod_radius * (1.0 - math.sin(theta / 2.0))


def rod_centers_axes(theta_deg: float, phi_deg: float, gap: float, rod_length: float, rod_width: float):
    theta = math.radians(theta_deg)
    phi = math.radians(phi_deg)
    tip_gap = rod_tip_gap(theta_deg, gap, rod_width / 2.0)

    rod1_out = math.pi
    rod2_out = math.pi - theta
    gap_axis = 0.5 * (rod1_out + rod2_out) - math.pi / 2.0
    gap_axis_x = math.cos(gap_axis)
    gap_axis_y = math.sin(gap_axis)

    tips = [
        (-0.5 * tip_gap * gap_axis_x, -0.5 * tip_gap * gap_axis_y),
        (0.5 * tip_gap * gap_axis_x, 0.5 * tip_gap * gap_axis_y),
    ]
    angles = [rod1_out, rod2_out]

    centers = []
    axes = []
    for tip, angle in zip(tips, angles):
        cx = tip[0] + 0.5 * rod_length * math.cos(angle)
        cy = tip[1] + 0.5 * rod_length * math.sin(angle)
        centers.append((cx, cy))
        axes.append(angle)

    dimer_shift_x = 0.0
    dimer_shift_y = 0.0

    rotated = []
    rotated_tips = []
    for point in tips:
        x = point[0] + dimer_shift_x
        y = point[1] + dimer_shift_y
        rotated_tips.append(
            (
                x * math.cos(phi) - y * math.sin(phi),
                x * math.sin(phi) + y * math.cos(phi),
            )
        )
    for (cx, cy), angle in zip(centers, axes):
        x = cx + dimer_shift_x
        y = cy + dimer_shift_y
        xr = x * math.cos(phi) - y * math.sin(phi)
        yr = x * math.sin(phi) + y * math.cos(phi)
        rotated.append((xr, yr, angle + phi))

    return rotated, rotated_tips, tip_gap


def capsule_outline(cx: float, cy: float, angle: float, length: float, width: float):
    radius = width / 2.0
    body_length = length - 2.0 * radius
    ux = math.cos(angle)
    uy = math.sin(angle)
    front = np.array([cx + 0.5 * body_length * ux, cy + 0.5 * body_length * uy])
    back = np.array([cx - 0.5 * body_length * ux, cy - 0.5 * body_length * uy])
    front_angles = np.linspace(angle - math.pi / 2.0, angle + math.pi / 2.0, 80)
    back_angles = np.linspace(angle + math.pi / 2.0, angle + 3.0 * math.pi / 2.0, 80)
    front_arc = front[:, None] + radius * np.vstack((np.cos(front_angles), np.sin(front_angles)))
    back_arc = back[:, None] + radius * np.vstack((np.cos(back_angles), np.sin(back_angles)))
    points = np.hstack((front_arc, back_arc, front_arc[:, :1]))
    return points[0], points[1]


def add_rod(ax, cx, cy, angle, length, width, color):
    x, y = capsule_outline(cx, cy, angle, length, width)
    ax.fill(x, y, facecolor=color, edgecolor="black", linewidth=1.2, alpha=0.9)


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
    rods, tips, tip_gap = rod_centers_axes(args.theta, args.phi, args.gap, args.rod_length, args.rod_width)

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

    tip_x = [tip[0] for tip in tips]
    tip_y = [tip[1] for tip in tips]
    ax.scatter(tip_x, tip_y, s=16, color="black", zorder=3, label="tip reference")
    ax.plot(tip_x, tip_y, color="black", linewidth=0.8, alpha=0.7)
    ax.text(0, arc_radius * 1.15, f"theta = {args.theta:g} deg", ha="center", color="tab:red")
    ax.text(
        0,
        -0.085,
        f"surface gap = {args.gap * 1000:g} nm; tip spacing = {tip_gap * 1000:.2f} nm",
        ha="center",
        fontsize=9,
    )

    lim = 0.095
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title(f"Nanodimer top view: theta={args.theta:g} deg, phi={args.phi:g} deg")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(args.output, dpi=250)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
