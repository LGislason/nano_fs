import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


# ============================================================
# User settings
# ============================================================

INCIDENT_FILE = "sphere_inc.txt"
SPHERE_FILE = "sphere_gold.txt"

OUTDIR = Path("sphere_plots")
OUTDIR.mkdir(exist_ok=True)

# Cut off wavelengths where incident flux is too small
INCIDENT_CUTOFF_FRACTION = 0.001

# Mie overlay settings
SPHERE_RADIUS_UM = 0.03
N_MEDIUM = 1  # air = 1.0, SiO2-ish = 1.45, water = 1.33
NORMALIZE_FOR_SHAPE = True


# ============================================================
# Load one-column Meep reflectance flux
# ============================================================

def load_reflectance_flux(filename):
    """
    Reads Meep output lines like:
        flux1:, frequency, flux

    Skips all non-flux log lines.
    """
    rows = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if "flux" not in line.lower():
                continue

            if ":" in line:
                line = line.split(":", 1)[1]

            parts = line.replace(",", " ").split()

            vals = []
            for p in parts:
                try:
                    vals.append(float(p))
                except ValueError:
                    pass

            if len(vals) >= 2:
                rows.append(vals[:2])

    if len(rows) == 0:
        raise ValueError(f"No reflectance flux rows found in {filename}")

    data = np.array(rows, dtype=float)

    freq = data[:, 0]
    flux = data[:, 1]
    wavelength = 1.0 / freq

    order = np.argsort(wavelength)

    return wavelength[order], flux[order]


# ============================================================
# Gold dielectric function for Mie overlay
# ============================================================

def gold_eps_meep_3term(wavelength_um):
    """
    Gold dielectric function matching the Meep/Scheme model:
    epsilon_inf + Drude + 3 Lorentz terms
    """

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

    omega = 1.0 / wavelength_um

    drude = (au_omega_d**2) / (0.0 - omega**2 - 1j * au_gamma_d * omega)

    lorentz_1 = (au_sigma_1 * au_omega_1**2) / (
        au_omega_1**2 - omega**2 - 1j * au_gamma_1 * omega
    )

    lorentz_2 = (au_sigma_2 * au_omega_2**2) / (
        au_omega_2**2 - omega**2 - 1j * au_gamma_2 * omega
    )

    lorentz_3 = (au_sigma_3 * au_omega_3**2) / (
        au_omega_3**2 - omega**2 - 1j * au_gamma_3 * omega
    )

    return au_epsilon_inf + drude + lorentz_1 + lorentz_2 + lorentz_3


def mie_scattering_curve(wavelength_um, radius_um, n_medium=1.0):
    try:
        import miepython
    except ImportError:
        print("miepython not installed. Run: pip install miepython")
        return None

    qsca = []

    for wl in wavelength_um:
        eps = gold_eps_meep_3term(wl)
        n_particle = np.sqrt(eps)

        # miepython expects absorbing materials with negative imaginary part.
        if np.imag(n_particle) > 0:
            n_particle = np.conj(n_particle)

        m = n_particle / n_medium
        x = 2 * np.pi * n_medium * radius_um / wl

        qext, qscat, qback, g = miepython.efficiencies_mx(m, x)
        qsca.append(qscat)

    return np.array(qsca)


# ============================================================
# Load files
# ============================================================

wl_inc, inc_flux = load_reflectance_flux(INCIDENT_FILE)
wl_sphere, sphere_flux = load_reflectance_flux(SPHERE_FILE)

if len(wl_inc) != len(wl_sphere):
    raise ValueError("Incident and sphere files have different numbers of points.")

if not np.allclose(wl_inc, wl_sphere, rtol=1e-9, atol=1e-12):
    raise ValueError("Incident and sphere wavelength grids do not match.")

wavelength = wl_inc


# ============================================================
# Calculate reflectance
# ============================================================

cutoff = INCIDENT_CUTOFF_FRACTION * np.max(np.abs(inc_flux))
valid = np.abs(inc_flux) > cutoff

reflectance = sphere_flux / inc_flux

wavelength = wavelength[valid]
reflectance = reflectance[valid]


# ============================================================
# Build Mie curve
# ============================================================

mie_qsca = mie_scattering_curve(wavelength, SPHERE_RADIUS_UM, N_MEDIUM)

if NORMALIZE_FOR_SHAPE:
    reflectance_plot = reflectance / np.nanmax(np.abs(reflectance))
    ylabel = "Normalized response"
else:
    reflectance_plot = reflectance
    ylabel = "Reflectance"

if mie_qsca is not None:
    if NORMALIZE_FOR_SHAPE:
        mie_plot = mie_qsca / np.nanmax(np.abs(mie_qsca))
    else:
        mie_plot = mie_qsca


# ============================================================
# Print summary
# ============================================================

print("Gold sphere reflectance with Mie overlay")
print(f"Incident file: {INCIDENT_FILE}")
print(f"Sphere file:   {SPHERE_FILE}")
print(f"Points kept: {np.sum(valid)} / {len(valid)}")
print(f"Incident cutoff: {cutoff:.4e}")
print(f"Wavelength range: {wavelength.min():.4f} to {wavelength.max():.4f} µm")
print(f"Reflectance range: {reflectance.min():.4e} to {reflectance.max():.4e}")
if mie_qsca is not None:
    print(f"Reflectance peak wavelength: {wavelength[np.argmax(reflectance_plot)]:.4f} µm")
    print(f"Mie peak wavelength:         {wavelength[np.argmax(mie_plot)]:.4f} µm")


# ============================================================
# Plot
# ============================================================

plt.figure(figsize=(10, 6))
plt.plot(wavelength, reflectance_plot, linewidth=2, label="Meep sphere reflectance")

if mie_qsca is not None:
    plt.plot(wavelength, mie_plot, linewidth=2, label="Mie scattering")

plt.xlabel("Wavelength (µm)")
plt.ylabel(ylabel)
plt.title("Gold Sphere Reflectance Validation")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

outfile = OUTDIR / "sphere_reflectance_with_mie.png"
plt.savefig(outfile, dpi=300)
plt.show()

print(f"Saved: {outfile}")
