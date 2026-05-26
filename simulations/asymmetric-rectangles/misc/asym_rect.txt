import meep as mp
import meep.materials as materials
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Stable Gold Asymmetric "- |" Unit Cell
# Normalized Reflectance (Ex & Ey)
# ============================================================

# ----------------------------
# Geometry parameters (µm)
# ----------------------------
period = 0.80
bar_length = 0.40
bar_width  = 0.12
thickness  = 0.06

cx = 0.12
cy = -0.1

# ----------------------------
# Simulation controls (metal-safe)
# ----------------------------
resolution = 60     # start conservative
dpml = 1.0
sz = 6.0            # more vertical separation for metal
run_time = 600

cell = mp.Vector3(period, period, sz)
boundary_layers = [mp.PML(dpml, direction=mp.Z)]

# ----------------------------
# Spectral band (inside gold validity)
# Gold model is reliable ~0.4–1.24 µm
# ----------------------------
wvl_min = 0.45
wvl_max = 1.20
nfreq = 120

fmin = 1 / wvl_max
fmax = 1 / wvl_min
fcen = 0.5 * (fmin + fmax)
df   = (fmax - fmin)

# ----------------------------
# Geometry definition
# ----------------------------
def make_geometry():
    return [
        mp.Block(
            size=mp.Vector3(bar_length, bar_width, thickness),
            center=mp.Vector3(cx - bar_length/2,
                              cy + bar_length/4,
                              0),
            material=materials.Au
        ),
        mp.Block(
            size=mp.Vector3(bar_width, bar_length, thickness),
            center=mp.Vector3(cx + bar_length/2,
                              cy + bar_length/4,
                              0),
            material=materials.Au
        ),
    ]

# ----------------------------
# Reflectance routine
# ----------------------------
def run_reflectance(pol_component, tag):

    src_z  = sz/2 - dpml - 0.2
    refl_z = sz/2 - dpml - 0.8

    sources = [
        mp.Source(
            mp.GaussianSource(frequency=fcen, fwidth=df),
            component=pol_component,
            center=mp.Vector3(0, 0, src_z),
            size=mp.Vector3(period, period, 0),
        )
    ]

    refl_region = mp.FluxRegion(
        center=mp.Vector3(0, 0, refl_z),
        size=mp.Vector3(period, period, 0),
    )

    # ---- Empty cell (incident flux) ----
    sim_empty = mp.Simulation(
        cell_size=cell,
        geometry=[],
        boundary_layers=boundary_layers,
        sources=sources,
        resolution=resolution,
    )

    refl_empty = sim_empty.add_flux(fcen, df, nfreq, refl_region)
    sim_empty.run(until=run_time)

    inc_flux = np.array(mp.get_fluxes(refl_empty))
    freqs = np.array(mp.get_flux_freqs(refl_empty))
    wvls = 1 / freqs

# Sort so wavelength increases left → right
    sort_idx = np.argsort(wvls)
    wvls = wvls[sort_idx]
    inc_flux = inc_flux[sort_idx]
    print(f"[{tag}] wavelength range (µm): {wvls.min():.4f} to {wvls.max():.4f}")


    print(f"[{tag}] wavelength range (µm): {wvls.min():.4f} to {wvls.max():.4f}")

    empty_flux_data = sim_empty.get_flux_data(refl_empty)

    # ---- Structure run ----
    sim = mp.Simulation(
        cell_size=cell,
        geometry=make_geometry(),
        boundary_layers=boundary_layers,
        sources=sources,
        resolution=resolution,
    )

    refl = sim.add_flux(fcen, df, nfreq, refl_region)
    sim.load_minus_flux_data(refl, empty_flux_data)

    sim.run(until=run_time)

    refl_flux = np.array(mp.get_fluxes(refl))
    R = np.abs(refl_flux) / np.maximum(np.abs(inc_flux), 1e-30)

    order = np.argsort(wvls)
    return wvls[order], R[order]


# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":

    w_ex, R_ex = run_reflectance(mp.Ex, "Ex")
    w_ey, R_ey = run_reflectance(mp.Ey, "Ey")

np.savetxt("reflectance_ex.txt",
           np.column_stack((w_ex, R_ex)),
           header="wavelength(um)  R_ex")

np.savetxt("reflectance_ey.txt",
           np.column_stack((w_ey, R_ey)),
           header="wavelength(um)  R_ey")

plt.figure(figsize=(7,5))

plt.plot(w_ex, R_ex, label="Ex incidence", linewidth=2)
plt.plot(w_ey, R_ey, label="Ey incidence", linewidth=2)

plt.xlabel("Wavelength (µm)")
plt.ylabel("Reflectance")
plt.xlim(min(w_ex), max(w_ex))
plt.ylim(0, 1.1)

plt.title("Gold Asymmetric '- |' Unit Cell (800 nm Period)")
plt.legend()
plt.tight_layout()

plt.savefig("gold_asymmetric_reflectance.png", dpi=300)
plt.show()

print("Saved reflectance_gold_Ex_Ey.png")

