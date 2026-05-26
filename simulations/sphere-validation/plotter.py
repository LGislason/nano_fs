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

            # Remove "flux1:" or similar
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

R = sphere_flux / inc_flux

wavelength = wavelength[valid]
R = R[valid]

# ============================================================
# Print summary
# ============================================================

print("Gold sphere reflectance")
print(f"Incident file: {INCIDENT_FILE}")
print(f"Sphere file:   {SPHERE_FILE}")
print(f"Points kept: {np.sum(valid)} / {len(valid)}")
print(f"Incident cutoff: {cutoff:.4e}")
print(f"Wavelength range: {wavelength.min():.4f} to {wavelength.max():.4f} µm")
print(f"Reflectance range: {R.min():.4e} to {R.max():.4e}")

# ============================================================
# Plot
# ============================================================

plt.figure(figsize=(10, 6))
plt.plot(wavelength, R, linewidth=2)



plt.xlabel("Wavelength (µm)")
plt.ylabel("Reflectance")
plt.title("Gold Sphere Reflectance Validation")
plt.grid(True, alpha=0.3)
plt.tight_layout()

plt.savefig(OUTDIR / "sphere_reflectance.png", dpi=300)
plt.show()