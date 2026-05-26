import math
import os
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent

MPLCONFIGDIR = BASE_DIR / ".mplconfig"
MPLCONFIGDIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR.resolve()))

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np


OUTDIR = BASE_DIR / "geometry_plots"
OUTDIR.mkdir(exist_ok=True)


# Cell and geometry parameters from nanorod_validation.ctl
sx = 0.80
sy = 0.80
sz = 2.40

substrate_epsilon = 2.1
substrate_thickness = 1.20
substrate_top_z = 0.0
substrate_bottom_z = -substrate_thickness

rod_length = 0.069
rod_width = 0.024
rod_radius = rod_width / 2
rod_body_length = rod_length - 2 * rod_radius
rod_z = rod_radius

gap = 0.010
rod_angle_deg = 80.0
rod_angle_rad = math.radians(rod_angle_deg)

rod1_axis = np.array([1.0, 0.0])
rod2_axis = np.array([math.cos(rod_angle_rad), math.sin(rod_angle_rad)])

rod1_cx_raw = -(gap / 2) - (rod_length / 2)
rod1_cy_raw = 0.0
rod2_cx_raw = (gap / 2) + 0.5 * rod_length * rod2_axis[0]
rod2_cy_raw = 0.5 * rod_length * rod2_axis[1]

dimer_shift_x = -0.5 * (rod1_cx_raw + rod2_cx_raw)
dimer_shift_y = -0.5 * (rod1_cy_raw + rod2_cy_raw)

rod1_center = np.array([rod1_cx_raw + dimer_shift_x, rod1_cy_raw + dimer_shift_y, rod_z])
rod2_center = np.array([rod2_cx_raw + dimer_shift_x, rod2_cy_raw + dimer_shift_y, rod_z])
single_center = np.array([0.0, 0.0, rod_z])


def capsule_outline_xy(center, axis_xy, body_length, radius, n_arc=80):
    axis_xy = np.asarray(axis_xy, dtype=float)
    axis_xy = axis_xy / np.linalg.norm(axis_xy)
    normal_xy = np.array([-axis_xy[1], axis_xy[0]])

    p1 = center[:2] - 0.5 * body_length * axis_xy
    p2 = center[:2] + 0.5 * body_length * axis_xy

    top_edge = np.vstack([p1 + radius * normal_xy, p2 + radius * normal_xy])
    bottom_edge = np.vstack([p2 - radius * normal_xy, p1 - radius * normal_xy])

    theta_right = np.linspace(0.0, math.pi, n_arc)
    theta_left = np.linspace(math.pi, 2 * math.pi, n_arc)

    right_arc = p2 + radius * (
        np.outer(np.cos(theta_right), normal_xy) + np.outer(np.sin(theta_right), axis_xy)
    )
    left_arc = p1 + radius * (
        np.outer(np.cos(theta_left), normal_xy) + np.outer(np.sin(theta_left), axis_xy)
    )

    return np.vstack([top_edge, right_arc[1:], bottom_edge, left_arc[1:]])


def add_box(ax, x0, x1, y0, y1, z0, z1, color, alpha):
    verts = [
        [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0)],
        [(x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)],
        [(x0, y0, z0), (x1, y0, z0), (x1, y0, z1), (x0, y0, z1)],
        [(x0, y1, z0), (x1, y1, z0), (x1, y1, z1), (x0, y1, z1)],
        [(x0, y0, z0), (x0, y1, z0), (x0, y1, z1), (x0, y0, z1)],
        [(x1, y0, z0), (x1, y1, z0), (x1, y1, z1), (x1, y0, z1)],
    ]
    ax.add_collection3d(
        Poly3DCollection(verts, facecolors=color, edgecolors=color, linewidths=0.4, alpha=alpha)
    )


def add_cylinder(ax, center, axis, length, radius, color, n_theta=36):
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    center = np.asarray(center, dtype=float)

    helper = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(axis, helper)) > 0.95:
        helper = np.array([0.0, 1.0, 0.0])

    u = np.cross(axis, helper)
    u = u / np.linalg.norm(u)
    v = np.cross(axis, u)

    t = np.linspace(0.0, 2 * math.pi, n_theta)
    offsets = radius * (np.outer(np.cos(t), u) + np.outer(np.sin(t), v))

    p0 = center - 0.5 * length * axis
    p1 = center + 0.5 * length * axis
    ring0 = p0 + offsets
    ring1 = p1 + offsets

    faces = []
    for i in range(n_theta - 1):
        faces.append([ring0[i], ring0[i + 1], ring1[i + 1], ring1[i]])
    faces.append([ring0[-1], ring0[0], ring1[0], ring1[-1]])

    ax.add_collection3d(
        Poly3DCollection(faces, facecolors=color, edgecolors="none", alpha=0.9)
    )


def add_sphere(ax, center, radius, color, n_theta=24, n_phi=12):
    center = np.asarray(center, dtype=float)
    theta = np.linspace(0, 2 * math.pi, n_theta)
    phi = np.linspace(0, math.pi, n_phi)
    theta, phi = np.meshgrid(theta, phi)

    x = center[0] + radius * np.cos(theta) * np.sin(phi)
    y = center[1] + radius * np.sin(theta) * np.sin(phi)
    z = center[2] + radius * np.cos(phi)

    ax.plot_surface(x, y, z, color=color, linewidth=0, antialiased=True, alpha=0.9)


def draw_capsule_3d(ax, center, axis_xy, body_length, radius, color):
    axis = np.array([axis_xy[0], axis_xy[1], 0.0], dtype=float)
    axis = axis / np.linalg.norm(axis)

    add_cylinder(ax, center, axis, body_length, radius, color)
    end1 = center - 0.5 * body_length * axis
    end2 = center + 0.5 * body_length * axis
    add_sphere(ax, end1, radius, color)
    add_sphere(ax, end2, radius, color)


def plot_top_view(ax, rods, title):
    substrate = Rectangle(
        (-sx / 2, -sy / 2),
        sx,
        sy,
        facecolor="#d9ecff",
        edgecolor="#4a6fa5",
        linewidth=1.2,
        alpha=0.5,
    )
    ax.add_patch(substrate)

    for center, axis_xy in rods:
        xy = capsule_outline_xy(center, axis_xy, rod_body_length, rod_radius)
        ax.fill(xy[:, 0], xy[:, 1], color="#d4a017", alpha=0.9, ec="#8a6510", lw=1.0)

    ax.set_title(title)
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")
    ax.set_aspect("equal")
    ax.set_xlim(-0.12, 0.12)
    ax.set_ylim(-0.12, 0.12)
    ax.grid(True, alpha=0.25)


def plot_side_view(ax, rods, title):
    ax.add_patch(
        Rectangle(
            (-sx / 2, substrate_bottom_z),
            sx,
            substrate_thickness,
            facecolor="#d9ecff",
            edgecolor="#4a6fa5",
            linewidth=1.2,
            alpha=0.6,
        )
    )

    for center, axis_xy in rods:
        x_min = center[0] - 0.5 * rod_length * abs(axis_xy[0]) - rod_radius * abs(axis_xy[1])
        x_max = center[0] + 0.5 * rod_length * abs(axis_xy[0]) + rod_radius * abs(axis_xy[1])
        z_min = center[2] - rod_radius
        z_max = center[2] + rod_radius
        ax.add_patch(
            Rectangle(
                (x_min, z_min),
                x_max - x_min,
                z_max - z_min,
                facecolor="#d4a017",
                edgecolor="#8a6510",
                linewidth=1.0,
                alpha=0.9,
            )
        )

    ax.axhline(substrate_top_z, color="#4a6fa5", linewidth=1.0)
    ax.set_title(title)
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("z (µm)")
    ax.set_xlim(-0.12, 0.12)
    ax.set_ylim(-1.28, 0.08)
    ax.grid(True, alpha=0.25)


def plot_3d_geometry(rods, title, filename):
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")

    add_box(ax, -sx / 2, sx / 2, -sy / 2, sy / 2, substrate_bottom_z, substrate_top_z, "#9cc4f5", 0.18)

    for center, axis_xy in rods:
        draw_capsule_3d(ax, center, axis_xy, rod_body_length, rod_radius, "#d4a017")

    ax.set_title(title)
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")
    ax.set_zlabel("z (µm)")
    ax.set_xlim(-0.12, 0.12)
    ax.set_ylim(-0.12, 0.12)
    ax.set_zlim(-0.08, 0.08)
    ax.view_init(elev=24, azim=-58)
    ax.set_box_aspect((1, 1, 0.7))
    plt.tight_layout()
    plt.savefig(OUTDIR / filename, dpi=300)
    plt.close(fig)


def plot_case(rods, stem, label):
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
    plot_top_view(axes[0], rods, f"{label}: top view")
    plot_side_view(axes[1], rods, f"{label}: side view")
    plt.tight_layout()
    plt.savefig(OUTDIR / f"{stem}_views.png", dpi=300)
    plt.close(fig)

    plot_3d_geometry(rods, f"{label}: 3D geometry", f"{stem}_3d.png")


def main():
    single_rods = [(single_center, rod1_axis)]
    dimer_rods = [(rod1_center, rod1_axis), (rod2_center, rod2_axis)]

    plot_case(single_rods, "single_nanorod", "Single nanorod")
    plot_case(dimer_rods, "nanorod_dimer", "Nanorod dimer")

    print("Saved geometry plots to:", OUTDIR.resolve())


if __name__ == "__main__":
    main()
