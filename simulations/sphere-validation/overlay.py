import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# User inputs
# ------------------------------------------------------------
MEEP_FILE = "sphere_gold.txt"   # change to your actual sphere output file

sphere_radius_um = 0.050        # 50 nm radius. change if needed
n_medium = 1.45                  # air = 1.0, SiO2-ish = 1.45
normalize_for_shape = True

# ------------------------------------------------------------
# Load Meep output robustly
# Works with Meep text files that contain non-numeric lines
# and flux lines like:
# flux1:, frequency, flux
# ------------------------------------------------------------
def load_meep_flux(filename):
    """
    Reads only Meep flux output lines like:

        flux1:, frequency, flux

    Ignores all other Meep log text.
    Keeps only frequency and the first flux column.
    """
    rows = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            # Only use actual flux output lines
            if not line.lower().startswith("flux"):
                continue

            # Remove label like "flux1:"
            if ":" in line:
                line = line.split(":", 1)[1]

            parts = line.replace(",", " ").split()

            vals = []
            for p in parts:
                try:
                    vals.append(float(p))
                except ValueError:
                    pass

            # Need frequency + one flux value
            if len(vals) >= 2:
                rows.append([vals[0], vals[1]])

    if len(rows) == 0:
        raise ValueError(f"No usable flux rows found in {filename}")

    data = np.array(rows, dtype=float)

    freq = data[:, 0]
    flux = data[:, 1]
    wavelength = 1.0 / freq

    order = np.argsort(wavelength)

    return wavelength[order], flux[order]

# ------------------------------------------------------------
# Simple Drude-like gold estimate for Mie overlay
# This is NOT your exact Meep gold model.
# It is only for checking approximate resonance location.
# ------------------------------------------------------------
def gold_eps_meep_3term(wavelength_um):
    """
    Gold dielectric function matching the Meep/Scheme model:

    epsilon_inf + Drude + 3 Lorentz terms

    Wavelength is in microns.
    Meep frequency unit assumes 1 um length scale, so:
        omega_meep_dimensionless = 1 / wavelength_um
    """

    # ----------------------------
    # Same constants as ctl file
    # ----------------------------
    omega_meep_si = 1.88365e15

    au_epsilon_inf = 4.8929

    au_omega_d = 1.2944e16 / omega_meep_si
    au_omega_1 = 1.3617e15 / omega_meep_si
    au_omega_2 = 4.1636e15 / omega_meep_si
    au_omega_3 = 5.0753e15 / omega_meep_si

    au_gamma_d = 1.0003e09 / omega_meep_si
    au_gamma_1 = (4.7356e14 * 2.0) / omega_meep_si
    au_gamma_2 = (4.4931e14 * 2.0) / omega_meep_si
    au_gamma_3 = (5.8469e14 * 2.0) / omega_meep_si

    au_sigma_1 = 4.7282
    au_sigma_2 = 0.72996
    au_sigma_3 = 1.5103

    # Meep dimensionless angular frequency
    omega = 1.0 / wavelength_um

    # ----------------------------
    # Meep susceptibility form:
    # eps = eps_inf + sigma*w0^2 / (w0^2 - w^2 - i*gamma*w)
    # ----------------------------

    # Drude term:
    # In ctl:
    # omega = 1e-20
    # sigma = (1e20 * au_omega_d)^2
    # Therefore sigma * omega0^2 = au_omega_d^2
    drude = (au_omega_d**2) / (
        0.0 - omega**2 - 1j * au_gamma_d * omega
    )

    lorentz_1 = (au_sigma_1 * au_omega_1**2) / (
        au_omega_1**2 - omega**2 - 1j * au_gamma_1 * omega
    )

    lorentz_2 = (au_sigma_2 * au_omega_2**2) / (
        au_omega_2**2 - omega**2 - 1j * au_gamma_2 * omega
    )

    lorentz_3 = (au_sigma_3 * au_omega_3**2) / (
        au_omega_3**2 - omega**2 - 1j * au_gamma_3 * omega
    )

    eps = au_epsilon_inf + drude + lorentz_1 + lorentz_2 + lorentz_3

    return eps

def mie_scattering_curve(wavelength_um, radius_um, n_medium=1.0):
    try:
        import miepython
    except ImportError:
        print("miepython not installed. Run: pip install miepython")
        return None

    qsca = []

    for wl in wavelength_um:
        eps = gold_eps_meep_3term(wl)

        # Convert dielectric function to refractive index
        n_particle = np.sqrt(eps)

        # miepython convention: absorbing material should be n - i*k
        if np.imag(n_particle) > 0:
            n_particle = np.conj(n_particle)

        m = n_particle / n_medium
        x = 2 * np.pi * n_medium * radius_um / wl

        qext, qscat, qback, g = miepython.efficiencies_mx(m, x)
        qsca.append(qscat)

    return np.array(qsca)


# ------------------------------------------------------------
# Load Meep curve
# ------------------------------------------------------------
wvl_meep, refl_meep = load_meep_flux(MEEP_FILE)

# Your reflectance sign convention:
# If the curve is negative, flip it.
if np.nanmean(refl_meep) < 0:
    refl_meep = -refl_meep

# Optional cutoff to avoid weird long-wavelength behavior
mask = (wvl_meep >= 0.40) & (wvl_meep <= 0.80)
wvl_meep = wvl_meep[mask]
refl_meep = refl_meep[mask]

# ------------------------------------------------------------
# Mie curve on same wavelength grid
# ------------------------------------------------------------
mie_qsca = mie_scattering_curve(wvl_meep,.03)

if normalize_for_shape:
    refl_plot = refl_meep / np.nanmax(np.abs(refl_meep))
    mie_plot = mie_qsca / np.nanmax(np.abs(mie_qsca))
    ylabel = "Normalized response"
else:
    refl_plot = refl_meep
    mie_plot = mie_qsca
    ylabel = "Response"

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
plt.figure(figsize=(9, 5))
plt.plot(wvl_meep, refl_plot, label="Meep sphere reflectance")
plt.plot(wvl_meep, mie_plot, label="Mie scattering, normalized")

plt.xlabel("Wavelength (µm)")
plt.ylabel(ylabel)
plt.title("Gold Sphere Validation: Meep vs Mie")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("sphere_meep_vs_mie.png", dpi=300)
plt.show()

print("Saved: sphere_meep_vs_mie.png")
print(f"Meep peak wavelength: {wvl_meep[np.argmax(refl_plot)]:.4f} µm")
print(f"Mie peak wavelength:  {wvl_meep[np.argmax(mie_plot)]:.4f} µm")