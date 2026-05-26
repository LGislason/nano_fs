import meep as mp
import matplotlib.pyplot as plt
import numpy as np

# =====================================================
# Geometry Verification Plot (Consistent with Main Sim)
# Units: microns (µm)
# =====================================================

resolution = 60

# --- Unit cell (µm) ---
period = 0.80
cell = mp.Vector3(period, period, 0)

# --- Bar dimensions (µm) ---
bar_length = 0.40
bar_width  = 0.12

cx = 0.12
cy = -0.1

geometry = [
    mp.Block(
        size=mp.Vector3(bar_length, bar_width, mp.inf),
        center=mp.Vector3(cx - bar_length/2,
                          cy + bar_length/4,
                          0),
        material=mp.Medium(epsilon=4)
    ),
    mp.Block(
        size=mp.Vector3(bar_width, bar_length, mp.inf),
        center=mp.Vector3(cx + bar_length/2,
                          cy + bar_length/4,
                          0),
        material=mp.Medium(epsilon=4)
    )
]

sim = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    resolution=resolution
)

sim.init_sim()

eps = sim.get_array(
    center=mp.Vector3(),
    size=cell,
    component=mp.Dielectric
)

# --- Plot with physical units ---
plt.figure(figsize=(6,5))

extent = [-period/2, period/2,
          -period/2, period/2]

plt.imshow(eps.T,
           origin="lower",
           extent=extent,
           cmap="viridis")

plt.xlabel("x (µm)")
plt.ylabel("y (µm)")
plt.colorbar(label="Relative Permittivity (ε)")
plt.title("Asymmetric '- |' Unit Cell (800 nm Period)")

plt.savefig("asymmetric_bar_geometry.png", dpi=300)
plt.show()

