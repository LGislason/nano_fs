import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# USER SETTINGS
# ============================================================
resolution = 100
period = 0.80
sz = 5.0
dpml = 1.5

bar_length = 0.40       # um
bar_width = 0.12        # um
metal_thickness = 0.06  # um

gold_z_center = -1.25

substrate_thickness = 1.0
substrate_center_z = -1.8

wavelength = .76  # um
frequency = 1.0 / wavelength

polarization = "Ey"
run_until = 500

# Set desired edge-to-edge gap here
gap = 0.1  # um = 200 nm

# y-position of both bars
bar_y = 0.0

# z planes for slices
z_metal = gold_z_center
z_above = gold_z_center + metal_thickness/2 + .03 # 1 um above top of metal

# ============================================================
# MATERIALS
# ============================================================
gold = mp.metal
sio2 = mp.Medium(epsilon=2.1)

# ============================================================
# RECTANGLE POSITIONS
# ============================================================
# Rectangle 1: horizontal bar (long in x)
# Rectangle 2: vertical bar (long in y)

half_x_horiz = bar_length / 2.0
half_x_vert  = bar_width / 2.0

# center-to-center separation needed for exact edge gap
center_sep = half_x_horiz + half_x_vert + gap

blk1_x = -center_sep / 2.0
blk2_x =  center_sep / 2.0

blk1_y = bar_y
blk2_y = bar_y

print(f"Center-to-center separation = {center_sep:.3f} um")
print(f"Target edge-to-edge gap     = {gap:.3f} um")

# ============================================================
# GEOMETRY
# ============================================================
geometry = [
    mp.Block(
        size=mp.Vector3(period, period, substrate_thickness),
        center=mp.Vector3(0, 0, substrate_center_z),
        material=sio2
    ),

    # Rectangle 1: horizontal (x-oriented)
    mp.Block(
        size=mp.Vector3(bar_length, bar_width, metal_thickness),
        center=mp.Vector3(blk1_x, blk1_y, gold_z_center),
        material=gold
    ),

    # Rectangle 2: vertical (y-oriented)
    mp.Block(
        size=mp.Vector3(bar_width, bar_length, metal_thickness),
        center=mp.Vector3(blk2_x, blk2_y, gold_z_center),
        material=gold
    )
]

# ============================================================
# CELL / SOURCE
# ============================================================
cell = mp.Vector3(period, period, sz)
pml_layers = [mp.PML(dpml, direction=mp.Z)]

halfz = sz / 2.0
src_z = halfz - dpml - 0.2

if polarization == "Ex":
    src_component = mp.Ex
elif polarization == "Ey":
    src_component = mp.Ey
else:
    raise ValueError("polarization must be 'Ex' or 'Ey'")

sources = [
    mp.Source(
        src=mp.ContinuousSource(frequency=frequency),
        component=src_component,
        center=mp.Vector3(0, 0, src_z),
        size=mp.Vector3(period, period, 0)
    )
]

sim = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    boundary_layers=pml_layers,
    sources=sources,
    resolution=resolution,
    k_point=mp.Vector3()
)

# ============================================================
# DFT FIELD MONITORS
# ============================================================
dft_xy_metal = sim.add_dft_fields(
    [mp.Ex, mp.Ey, mp.Ez],
    frequency, 0, 1,
    center=mp.Vector3(0, 0, z_metal),
    size=mp.Vector3(period, period, 0)
)

dft_xy_above = sim.add_dft_fields(
    [mp.Ex, mp.Ey, mp.Ez],
    frequency, 0, 1,
    center=mp.Vector3(0, 0, z_above),
    size=mp.Vector3(period, period, 0)
)

dft_xz = sim.add_dft_fields(
    [mp.Ex, mp.Ey, mp.Ez],
    frequency, 0, 1,
    center=mp.Vector3(0, bar_y, 0),
    size=mp.Vector3(period, 0, sz)
)

# ============================================================
# RUN
# ============================================================
sim.run(until=run_until)

# ============================================================
# EXTRACT DFT FIELDS
# ============================================================
Ex_xy_m = sim.get_dft_array(dft_xy_metal, mp.Ex, 0)
Ey_xy_m = sim.get_dft_array(dft_xy_metal, mp.Ey, 0)
Ez_xy_m = sim.get_dft_array(dft_xy_metal, mp.Ez, 0)

Ex_xy_a = sim.get_dft_array(dft_xy_above, mp.Ex, 0)
Ey_xy_a = sim.get_dft_array(dft_xy_above, mp.Ey, 0)
Ez_xy_a = sim.get_dft_array(dft_xy_above, mp.Ez, 0)

Ex_xz = sim.get_dft_array(dft_xz, mp.Ex, 0)
Ey_xz = sim.get_dft_array(dft_xz, mp.Ey, 0)
Ez_xz = sim.get_dft_array(dft_xz, mp.Ez, 0)

# ============================================================
# TOTAL FIELD INTENSITY
# ============================================================
E2_xy_metal = np.abs(Ex_xy_m)**2 + np.abs(Ey_xy_m)**2 + np.abs(Ez_xy_m)**2
E2_xy_above = np.abs(Ex_xy_a)**2 + np.abs(Ey_xy_a)**2 + np.abs(Ez_xy_a)**2
E2_xz = np.abs(Ex_xz)**2 + np.abs(Ey_xz)**2 + np.abs(Ez_xz)**2

eps_xy_metal = sim.get_array(
    center=mp.Vector3(0, 0, z_metal),
    size=mp.Vector3(period, period, 0),
    component=mp.Dielectric
)

eps_xy_above = sim.get_array(
    center=mp.Vector3(0, 0, z_above),
    size=mp.Vector3(period, period, 0),
    component=mp.Dielectric
)

eps_xz = sim.get_array(
    center=mp.Vector3(0, bar_y, 0),
    size=mp.Vector3(period, 0, sz),
    component=mp.Dielectric
)

# ============================================================
# NORMALIZE
# ============================================================
global_max = max(
    np.max(E2_xy_metal),
    np.max(E2_xy_above),
    np.max(E2_xz)
)

if global_max == 0:
    raise RuntimeError("Field intensity is zero everywhere.")

E2_xy_metal /= global_max
E2_xy_above /= global_max
E2_xz /= global_max

# ============================================================
# AXES EXTENTS
# ============================================================
xy_extent = [-period/2, period/2, -period/2, period/2]
xz_extent = [-period/2, period/2, -sz/2, sz/2]

# ============================================================
# PLOTTING
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

im0 = axes[0].imshow(
    np.rot90(E2_xy_metal),
    extent=xy_extent,
    origin="lower",
    interpolation="nearest",
    aspect="equal"
)
axes[0].contour(
    np.rot90(eps_xy_metal),
    levels=[1.5],
    colors="k",
    linewidths=1.0,
    extent=xy_extent
)
axes[0].set_title(f"XY at metal plane\n{polarization}, λ = {wavelength:.3f} µm")
axes[0].set_xlabel("x (µm)")
axes[0].set_ylabel("y (µm)")

im1 = axes[1].imshow(
    np.rot90(E2_xy_above),
    extent=xy_extent,
    origin="lower",
    interpolation="nearest",
    aspect="equal"
)
axes[1].contour(
    np.rot90(eps_xy_above),
    levels=[1.5],
    colors="k",
    linewidths=1.0,
    extent=xy_extent
)
axes[1].set_title(f"XY above metal\n{polarization}, λ = {wavelength:.3f} µm")
axes[1].set_xlabel("x (µm)")
axes[1].set_ylabel("y (µm)")

im2 = axes[2].imshow(
    np.rot90(E2_xz),
    extent=xz_extent,
    origin="lower",
    interpolation="nearest",
    aspect="auto"
)
axes[2].contour(
    np.rot90(eps_xz),
    levels=[1.5],
    colors="k",
    linewidths=1.0,
    extent=xz_extent
)
axes[2].set_title(f"XZ slice at y = {bar_y:.3f} µm\n{polarization}, λ = {wavelength:.3f} µm")
axes[2].set_xlabel("x (µm)")
axes[2].set_ylabel("z (µm)")

cbar = fig.colorbar(im2, ax=axes.ravel().tolist(), shrink=0.92)
cbar.set_label("Normalized |E|^2")

outfile = f"fields_gap{gap}_{polarization}_{wavelength:.3f}um.png"
plt.savefig(outfile, dpi=300)
plt.show()

print(f"Saved figure to {outfile}")