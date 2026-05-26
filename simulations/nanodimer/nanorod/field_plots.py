#!/usr/bin/env python3
"""Generate local Meep field-intensity plots for the nanorod validation model."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import meep as mp
import numpy as np


BASE_DIR = Path(__file__).resolve().parent
DEFAULT_OUTDIR = BASE_DIR / "field_plots"


def gold_medium() -> mp.Medium:
    omega_meep = 1.88365e15

    au_epsilon_inf = 4.8929
    au_omega_d = 1.2944e16 / omega_meep
    au_omega_1 = 1.3617e15 / omega_meep
    au_omega_2 = 4.1636e15 / omega_meep
    au_omega_3 = 5.0753e15 / omega_meep

    au_gamma_d = 1.0003e09 / omega_meep
    au_gamma_1 = (4.7356e14 * 2) / omega_meep
    au_gamma_2 = (4.4931e14 * 2) / omega_meep
    au_gamma_3 = (5.8469e14 * 2) / omega_meep

    return mp.Medium(
        epsilon=au_epsilon_inf,
        E_susceptibilities=[
            mp.DrudeSusceptibility(
                frequency=1e-20,
                gamma=au_gamma_d,
                sigma=(1e20 * au_omega_d) ** 2,
            ),
            mp.LorentzianSusceptibility(
                frequency=au_omega_1,
                gamma=au_gamma_1,
                sigma=4.7282,
            ),
            mp.LorentzianSusceptibility(
                frequency=au_omega_2,
                gamma=au_gamma_2,
                sigma=0.72996,
            ),
            mp.LorentzianSusceptibility(
                frequency=au_omega_3,
                gamma=au_gamma_3,
                sigma=1.5103,
            ),
        ],
    )


def capsule_rod(
    center: tuple[float, float, float],
    axis: tuple[float, float, float],
    rod_length: float,
    rod_radius: float,
    material: mp.Medium,
) -> list[mp.GeometricObject]:
    cx, cy, cz = center
    ux, uy, uz = axis
    body_length = rod_length - 2 * rod_radius

    return [
        mp.Cylinder(
            radius=rod_radius,
            height=body_length,
            center=mp.Vector3(cx, cy, cz),
            axis=mp.Vector3(ux, uy, uz),
            material=material,
        ),
        mp.Sphere(
            radius=rod_radius,
            center=mp.Vector3(
                cx + 0.5 * body_length * ux,
                cy + 0.5 * body_length * uy,
                cz + 0.5 * body_length * uz,
            ),
            material=material,
        ),
        mp.Sphere(
            radius=rod_radius,
            center=mp.Vector3(
                cx - 0.5 * body_length * ux,
                cy - 0.5 * body_length * uy,
                cz - 0.5 * body_length * uz,
            ),
            material=material,
        ),
    ]


def nanorod_geometry(geom: str, gap: float, rod_angle_deg: float) -> tuple[list[mp.GeometricObject], float]:
    sx = sy = 0.80
    substrate_thickness = 1.20
    sio2 = mp.Medium(epsilon=2.1)
    gold = gold_medium()

    substrate = mp.Block(
        size=mp.Vector3(sx, sy, substrate_thickness),
        center=mp.Vector3(0, 0, -0.5 * substrate_thickness),
        material=sio2,
    )

    rod_length = 0.069
    rod_width = 0.024
    rod_radius = 0.5 * rod_width
    rod_z = rod_radius

    if geom == "single":
        return [substrate, *capsule_rod((0, 0, rod_z), (1, 0, 0), rod_length, rod_radius, gold)], rod_z

    angle = np.deg2rad(rod_angle_deg)
    rod2_ux = float(np.cos(angle))
    rod2_uy = float(np.sin(angle))

    rod1_cx_raw = -0.5 * gap - 0.5 * rod_length
    rod1_cy_raw = 0.0
    rod2_cx_raw = 0.5 * gap + 0.5 * rod_length * rod2_ux
    rod2_cy_raw = 0.5 * rod_length * rod2_uy

    shift_x = -0.5 * (rod1_cx_raw + rod2_cx_raw)
    shift_y = -0.5 * (rod1_cy_raw + rod2_cy_raw)

    rod1 = capsule_rod(
        (rod1_cx_raw + shift_x, rod1_cy_raw + shift_y, rod_z),
        (1, 0, 0),
        rod_length,
        rod_radius,
        gold,
    )
    rod2 = capsule_rod(
        (rod2_cx_raw + shift_x, rod2_cy_raw + shift_y, rod_z),
        (rod2_ux, rod2_uy, 0),
        rod_length,
        rod_radius,
        gold,
    )
    return [substrate, *rod1, *rod2], rod_z


def get_field_intensity(sim: mp.Simulation, dft: object) -> np.ndarray:
    ex = sim.get_dft_array(dft, mp.Ex, 0)
    ey = sim.get_dft_array(dft, mp.Ey, 0)
    ez = sim.get_dft_array(dft, mp.Ez, 0)
    return np.abs(ex) ** 2 + np.abs(ey) ** 2 + np.abs(ez) ** 2


def add_contours(ax: plt.Axes, eps: np.ndarray, extent: list[float]) -> None:
    eps_plot = np.rot90(np.real(eps))
    levels = [level for level in (1.5, 3.0) if eps_plot.min() < level < eps_plot.max()]
    if levels:
        ax.contour(eps_plot, levels=levels, colors="white", linewidths=0.7, extent=extent)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geom", choices=["single", "dimer"], default="dimer")
    parser.add_argument("--pol", choices=["Ex", "Ey"], default="Ex")
    parser.add_argument("--wavelength", type=float, default=0.760, help="Wavelength in microns.")
    parser.add_argument("--resolution", type=int, default=80, help="Pixels per micron.")
    parser.add_argument("--run-until", type=float, default=200)
    parser.add_argument("--gap", type=float, default=0.010, help="Dimer gap in microns.")
    parser.add_argument("--angle", type=float, default=80.0, help="Dimer rod angle in degrees.")
    parser.add_argument(
        "--normalize",
        choices=["panel", "global"],
        default="panel",
        help="Normalize each panel separately or all panels together.",
    )
    parser.add_argument(
        "--vmax-percentile",
        type=float,
        default=99.5,
        help="Percentile used as the color scale maximum before clipping.",
    )
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    sx = sy = 0.80
    sz = 2.40
    dpml = 0.25
    halfz = 0.5 * sz
    frequency = 1.0 / args.wavelength
    src_comp = mp.Ex if args.pol == "Ex" else mp.Ey

    geometry, rod_z = nanorod_geometry(args.geom, args.gap, args.angle)
    source_z = halfz - dpml - 0.12
    field_sx = sx - 2 * dpml - 0.02
    field_sy = sy - 2 * dpml - 0.02

    sim = mp.Simulation(
        cell_size=mp.Vector3(sx, sy, sz),
        boundary_layers=[
            mp.PML(dpml, direction=mp.X),
            mp.PML(dpml, direction=mp.Y),
            mp.PML(dpml, direction=mp.Z),
        ],
        geometry=geometry,
        sources=[
            mp.Source(
                src=mp.ContinuousSource(frequency=frequency),
                component=src_comp,
                center=mp.Vector3(0, 0, source_z),
                size=mp.Vector3(field_sx, field_sy, 0),
            )
        ],
        resolution=args.resolution,
        ensure_periodicity=False,
    )

    z_metal = rod_z
    z_above = rod_z + 0.030

    dft_xy_metal = sim.add_dft_fields(
        [mp.Ex, mp.Ey, mp.Ez],
        frequency,
        0,
        1,
        center=mp.Vector3(0, 0, z_metal),
        size=mp.Vector3(field_sx, field_sy, 0),
    )
    dft_xy_above = sim.add_dft_fields(
        [mp.Ex, mp.Ey, mp.Ez],
        frequency,
        0,
        1,
        center=mp.Vector3(0, 0, z_above),
        size=mp.Vector3(field_sx, field_sy, 0),
    )
    dft_xz = sim.add_dft_fields(
        [mp.Ex, mp.Ey, mp.Ez],
        frequency,
        0,
        1,
        center=mp.Vector3(0, 0, 0),
        size=mp.Vector3(field_sx, 0, sz - 2 * dpml),
    )

    sim.run(until=args.run_until)

    e2_xy_metal = get_field_intensity(sim, dft_xy_metal)
    e2_xy_above = get_field_intensity(sim, dft_xy_above)
    e2_xz = get_field_intensity(sim, dft_xz)

    eps_xy_metal = sim.get_array(
        center=mp.Vector3(0, 0, z_metal),
        size=mp.Vector3(field_sx, field_sy, 0),
        component=mp.Dielectric,
    )
    eps_xy_above = sim.get_array(
        center=mp.Vector3(0, 0, z_above),
        size=mp.Vector3(field_sx, field_sy, 0),
        component=mp.Dielectric,
    )
    eps_xz = sim.get_array(
        center=mp.Vector3(0, 0, 0),
        size=mp.Vector3(field_sx, 0, sz - 2 * dpml),
        component=mp.Dielectric,
    )

    fields = [e2_xy_metal, e2_xy_above, e2_xz]
    if max(field.max() for field in fields) <= 0:
        raise RuntimeError("Field intensity is zero everywhere.")

    if args.normalize == "global":
        vmax = np.percentile(np.concatenate([field.ravel() for field in fields]), args.vmax_percentile)
        fields = [np.clip(field / vmax, 0, 1) for field in fields]
        colorbar_label = f"Global-normalized |E|^2, p{args.vmax_percentile:g}=1"
    else:
        normalized_fields = []
        for field in fields:
            vmax = np.percentile(field, args.vmax_percentile)
            normalized_fields.append(np.clip(field / vmax, 0, 1))
        fields = normalized_fields
        colorbar_label = f"Panel-normalized |E|^2, p{args.vmax_percentile:g}=1"

    e2_xy_metal, e2_xy_above, e2_xz = fields

    xy_extent = [-0.5 * field_sx, 0.5 * field_sx, -0.5 * field_sy, 0.5 * field_sy]
    xz_extent = [-0.5 * field_sx, 0.5 * field_sx, -0.5 * (sz - 2 * dpml), 0.5 * (sz - 2 * dpml)]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
    images = [
        (axes[0], e2_xy_metal, eps_xy_metal, xy_extent, f"XY at rod midplane z={z_metal:.3f} um"),
        (axes[1], e2_xy_above, eps_xy_above, xy_extent, f"XY above rods z={z_above:.3f} um"),
        (axes[2], e2_xz, eps_xz, xz_extent, "XZ slice at y=0"),
    ]

    im = None
    for ax, field, eps, extent, title in images:
        im = ax.imshow(
            np.rot90(field),
            extent=extent,
            origin="lower",
            interpolation="nearest",
            aspect="auto" if "XZ" in title else "equal",
            cmap="inferno",
            vmin=0,
            vmax=1,
        )
        add_contours(ax, eps, extent)
        ax.set_title(title)
        ax.set_xlabel("x (um)")
        ax.set_ylabel("z (um)" if "XZ" in title else "y (um)")

    fig.colorbar(im, ax=axes, shrink=0.88, label=colorbar_label)
    fig.suptitle(
        f"Nanorod {args.geom} field intensity: {args.pol}, "
        f"lambda={args.wavelength:.3f} um, res={args.resolution}"
    )

    args.outdir.mkdir(parents=True, exist_ok=True)
    outfile = args.outdir / (
        f"field_{args.geom}_{args.pol}_lambda{args.wavelength:.3f}um_"
        f"gap{args.gap:.3f}um_res{args.resolution}_{args.normalize}.png"
    )
    fig.savefig(outfile, dpi=250)
    plt.close(fig)
    print(f"Wrote {outfile}")


if __name__ == "__main__":
    main()
