import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# ===============================
# Units: microns (um)
# ===============================

# ---- Metasurface unit cell ----
period = 0.8            # 800 nm
bar_length = 0.40       # 400 nm
bar_width  = 0.12       # 120 nm
thickness  = 0.06       # 60 nm

# Motif placement (break symmetry)
offset_x = 0.10
offset_y = -0.05

# ---- Materials ----
eps_struct = 4.0        # placeholder dielectric; swap later for realistic n(λ)

# ---- Simulation controls (laptop-friendly defaults) ----
resolution = 50         # start 40-60 on laptop
dpml = 1.0              # PML thickness (um)
sz = 4.0                # total z height (um). Keep ~3-5 um for speed
run_time = 400          # time units; increase if spectra look noisy

# ---- Spectrum settings ----
wvl_min = 0.8           # um
wvl_max = 1.6           # um
nfreq = 80              # 50-120 is fine on laptop

fmin = 1 / wvl_max   # lowest frequency corresponds to longest wavelength
fmax = 1 / wvl_min   # highest frequency corresponds to shortest wavelength
fcen = 0.5 * (fmin + fmax)
df   = (fmax - fmin)

# ===============================
# Geometry: "- |" motif
# ===============================
def make_geometry():

    cx = 0.15
    cy = -0.1

    return [
        # Left horizontal bar
        mp.Block(
            size=mp.Vector3(bar_length, bar_width, thickness),
            center=mp.Vector3(cx - bar_length/2,
                              cy + bar_length/4,
                              0),
            material=gold_3term
        ),

        # Right vertical bar
        mp.Block(
            size=mp.Vector3(bar_width, bar_length, thickness),
            center=mp.Vector3(cx + bar_length/2,
                              cy + bar_length/4,
                              0),
            material=gold_3term
        ),
    ]


# ===============================
# Core routine: normalized reflectance for a given polarization
# ===============================
def run_reflectance(pol_component, tag):
    cell = mp.Vector3(period, period, sz)
    boundary_layers = [mp.PML(dpml, direction=mp.Z)]

    # Source plane near top (launch downward)
    src_z = + (sz/2 - dpml - 0.2)

    # Reflection monitor plane above structure (between source and structure)
    refl_z = + (sz/2 - dpml - 0.6)

    sources = [
        mp.Source(
            mp.GaussianSource(frequency=fcen, fwidth=df),
            component=pol_component,
            center=mp.Vector3(0, 0, src_z),
            size=mp.Vector3(period, period, 0)  # uniform over unit cell
        )
    ]

    # ---------- 1) Empty-cell run (incident flux) ----------
    sim_empty = mp.Simulation(
        cell_size=cell,
        geometry=[],  # empty
        boundary_layers=boundary_layers,
        sources=sources,
        resolution=resolution,
    )

    refl_region = mp.FluxRegion(center=mp.Vector3(0, 0, refl_z),
                                size=mp.Vector3(period, period, 0))
    refl_empty = sim_empty.add_flux(fcen, df, nfreq, refl_region)

    sim_empty.run(until=run_time)

    inc_flux = np.array(mp.get_fluxes(refl_empty))
    freqs = np.array(mp.get_flux_freqs(refl_empty))
    print("Wavelength range (um):",
       (1/freqs).min(),
       (1/freqs).max())
 

    # Save the "incident" flux data (useful for debugging)
    np.savetxt(f"{tag}_incident_flux.txt",
               np.column_stack([1/freqs, inc_flux]),
               header="wavelength_um incident_flux")

    # Store flux data for subtraction
    empty_flux_data = sim_empty.get_flux_data(refl_empty)

    # ---------- 2) Structure run (reflected flux, subtract incident) ----------
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
    freqs2 = np.array(mp.get_flux_freqs(refl))

    # Sanity: freqs should match
    if not np.allclose(freqs, freqs2):
        raise RuntimeError("Frequency grids do not match between empty and structure runs.")

    # In MEEP, after load_minus_flux_data, this flux is the scattered (reflected) contribution.
    # Reflectance: R = reflected / incident
    # Use absolute value to be robust to sign conventions.
    R = np.abs(refl_flux) / np.maximum(np.abs(inc_flux), 1e-30)
    wvls = 1 / freqs

    # Sort by wavelength (optional, but nice)
    order = np.argsort(wvls)
    wvls = wvls[order]
    R = R[order]

    np.savetxt(f"{tag}_reflectance.txt",
               np.column_stack([wvls, R]),
               header="wavelength_um reflectance")

    return wvls, R

# ===============================
# Run Ex and Ey, plot together
# ===============================
if __name__ == "__main__":
    # Keep these two runs separate so you can see progress
    w_ex, R_ex = run_reflectance(mp.Ex, "Ex")
    w_ey, R_ey = run_reflectance(mp.Ey, "Ey")

    plt.figure(figsize=(6,4))
    plt.plot(w_ex, R_ex, label="Ex incidence")
    plt.plot(w_ey, R_ey, label="Ey incidence")
    plt.xlabel("Wavelength (μm)")
    plt.ylabel("Reflectance R")
    plt.title("Normalized Reflectance: Asymmetric - | Unit Cell")
    plt.gca().invert_xaxis()
    plt.legend()
    plt.tight_layout()
    plt.savefig("reflectance_Ex_Ey.png", dpi=300)
    plt.show()

    print("Saved: Ex_reflectance.txt, Ey_reflectance.txt, reflectance_Ex_Ey.png")

