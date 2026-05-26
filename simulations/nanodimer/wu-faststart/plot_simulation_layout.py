#!/usr/bin/env python3
"""Plot the current source/flux/PML layout used by nanorod_wu_fast.ctl."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sx", type=float, default=0.80)
    parser.add_argument("--sz", type=float, default=1.20)
    parser.add_argument("--dpml", type=float, default=0.12)
    parser.add_argument("--rod-length", type=float, default=0.069)
    parser.add_argument("--rod-width", type=float, default=0.024)
    parser.add_argument("-o", "--output", type=Path, default=Path("simulation_layout_xz.png"))
    args = parser.parse_args()

    halfx = args.sx / 2.0
    halfz = args.sz / 2.0

    src_z = halfz - args.dpml - 0.10
    refl_z = halfz - args.dpml - 0.18
    tran_z = -halfz + args.dpml + 0.18

    fig, ax = plt.subplots(figsize=(6.2, 4.6))

    # Full cell.
    ax.add_patch(
        Rectangle(
            (-halfx, -halfz),
            args.sx,
            args.sz,
            fill=False,
            linewidth=1.6,
            color="black",
            label="cell boundary",
        )
    )

    # PML regions in z, shown because this cross-section is about plane spacing.
    ax.add_patch(
        Rectangle(
            (-halfx, halfz - args.dpml),
            args.sx,
            args.dpml,
            facecolor="0.85",
            edgecolor="none",
            label="PML",
        )
    )
    ax.add_patch(
        Rectangle(
            (-halfx, -halfz),
            args.sx,
            args.dpml,
            facecolor="0.85",
            edgecolor="none",
        )
    )

    # Approximate rod envelope at z=0.
    ax.add_patch(
        Rectangle(
            (-args.rod_length / 2.0, -args.rod_width / 2.0),
            args.rod_length,
            args.rod_width,
            facecolor="#d8a020",
            edgecolor="#6f4e00",
            linewidth=1.2,
            label="nanorod envelope",
        )
    )

    planes = [
        ("source plane, z=+0.38 um", src_z, "#1f77b4", "-"),
        ("reflection flux plane, z=+0.30 um", refl_z, "#d62728", "--"),
        ("transmission flux plane, z=-0.30 um", tran_z, "#2ca02c", "--"),
    ]

    for label, z, color, linestyle in planes:
        ax.axhline(z, color=color, linestyle=linestyle, linewidth=1.7, label=label)

    ax.axhline(0, color="0.25", linewidth=0.8, label="dimer center")

    ax.set_xlim(-halfx - 0.02, halfx + 0.02)
    ax.set_ylim(-halfz - 0.04, halfz + 0.04)
    ax.set_aspect("equal")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("z (um)")
    ax.set_title("Simulation layout, x-z cross-section")
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.14),
        fontsize=8,
        frameon=True,
        ncols=2,
    )
    fig.tight_layout()
    fig.savefig(args.output, dpi=250)
    print(f"Wrote {args.output}")
    print(f"source z = {src_z:.3f} um, distance from dimer = {abs(src_z):.3f} um")
    print(f"reflection flux z = {refl_z:.3f} um, distance from dimer = {abs(refl_z):.3f} um")
    print(f"transmission flux z = {tran_z:.3f} um, distance from dimer = {abs(tran_z):.3f} um")


if __name__ == "__main__":
    main()
