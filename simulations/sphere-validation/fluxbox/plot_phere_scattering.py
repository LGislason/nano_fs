import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# User settings
# ============================================================

INC_FILE = "sphere_inc.txt"
SCAT_FILE = "sphere_scat.txt"

OUTDIR = Path("sphere_plots")
OUTDIR.mkdir(exist_ok=True)

# Cut off weak source regions
INCIDENT_CUTOFF_FRACTION = 0.01

# Set this to your actual sphere radius from the ctl
SPHERE_RADIUS_UM = 0.030  # 30 nm radius

# Optional Mie overlay
DO_MIE_OVERLAY = True


# ============================================================
# Load Meep flux-box output
# ============================================================

def load_flux_box(filename):
    """
    Reads Meep flux output lines like:

    flux1:, freq, xplus, xminus, yplus, yminus, zplus, zminus

    Skips all non-flux log lines.
    """
    rows = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if not line.lower().startswith("flux"):
                continue

            # Remove "flux1:" label
            if ":" in line:
                line = line.split(":", 1)[1]

            parts = line.replace(",", " ").split()

            vals = []
            for p in parts:
                try:
                    vals.append(float(p))
                except ValueError:
                    pass

            # Need frequency + six flux columns
            if len(vals) >= 7:
                rows.append(vals[:7])

    if len(rows) == 0:
        raise ValueError(f"No usable flux-box rows found in {filename}")

    data = np.array(rows, dtype=float)

    freq = data[:, 0]
    flux = data[:, 1:7]
    wavelength = 1.0 / freq

    order = np.argsort(wavelength)

    return wavelength[order], flux[order, :]


# ============================================================
# Optional Mie curve
# ============================================================

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

def mie_qsca(wavelength_um, radius_um, n_medium=1.0):
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


# ============================================================
# Main calculation
# ============================================================

wavelength_inc, inc_flux = load_flux_box(INC_FILE)
wavelength_scat, scat_flux = load_flux_box(SCAT_FILE)

if len(wavelength_inc) != len(wavelength_scat):
    raise ValueError("Incident and scattered files have different numbers of points.")

if not np.allclose(wavelength_inc, wavelength_scat, rtol=1e-9, atol=1e-12):
    raise ValueError("Incident and scattered wavelength grids do not match.")

wavelength = wavelength_scat

# Flux columns:
# 0 = xplus
# 1 = xminus
# 2 = yplus
# 3 = yminus
# 4 = zplus
# 5 = zminus

xplus = scat_flux[:, 0]
xminus = scat_flux[:, 1]
yplus = scat_flux[:, 2]
yminus = scat_flux[:, 3]
zplus = scat_flux[:, 4]
zminus = scat_flux[:, 5]

# Outward flux box sum.
# Plus-side monitors are outward positive.
# Minus-side monitors must be subtracted because outward normal is negative.
scattered_power = xplus - xminus + yplus - yminus + zplus - zminus

# Use incident flux magnitude to cut off weak source regions.
# The incident run has nonzero flux mainly through the z faces.
incident_strength = np.max(np.abs(inc_flux), axis=1)
cutoff = INCIDENT_CUTOFF_FRACTION * np.max(incident_strength)
valid = incident_strength > cutoff

wavelength = wavelength[valid]
scattered_power = scattered_power[valid]

# Make positive if sign is globally flipped
if np.nanmean(scattered_power) < 0:
    scattered_power = -scattered_power

# Normalize for shape comparison
scattered_norm = scattered_power / np.nanmax(np.abs(scattered_power))

print("Gold sphere flux-box scattering")
print(f"Incident file: {INC_FILE}")
print(f"Scattered file: {SCAT_FILE}")
print(f"Points kept: {np.sum(valid)} / {len(valid)}")
print(f"Wavelength range: {wavelength.min():.4f} to {wavelength.max():.4f} µm")
print(f"Scattered power range: {scattered_power.min():.4e} to {scattered_power.max():.4e}")
print(f"Peak wavelength: {wavelength[np.argmax(scattered_norm)]:.4f} µm")


# ============================================================
# Plot Meep scattering
# ============================================================

plt.figure(figsize=(10, 6))
plt.plot(wavelength, scattered_norm, linewidth=2, label="Meep flux-box scattering")

plt.xlabel("Wavelength (µm)")
plt.ylabel("Normalized scattered power")
plt.title("Gold Sphere Validation: Meep Flux-Box Scattering")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(OUTDIR / "sphere_flux_box_scattering.png", dpi=300)
plt.show()


# ============================================================
# Plot Meep vs Mie, normalized shape
# ============================================================

if DO_MIE_OVERLAY:
    qsca = mie_qsca(wavelength, SPHERE_RADIUS_UM)

    if qsca is not None:
        qsca_norm = qsca / np.nanmax(np.abs(qsca))

        print(f"Mie peak wavelength: {wavelength[np.argmax(qsca_norm)]:.4f} µm")

        plt.figure(figsize=(10, 6))
        plt.plot(wavelength, scattered_norm, linewidth=2, label="Meep flux-box scattering")
        plt.plot(wavelength, qsca_norm, linewidth=2, label="Mie scattering, normalized")

        plt.xlabel("Wavelength (µm)")
        plt.ylabel("Normalized response")
        plt.title("Gold Sphere Validation: Meep vs Mie")
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(OUTDIR / "sphere_meep_vs_mie_flux_box.png", dpi=300)
        plt.show()