import numpy as np
import matplotlib.pyplot as plt


def load_flux(filename, tag):
    freq, flux = [], []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(tag):
                parts = line.split(",", maxsplit=2)
                freq.append(float(parts[1]))
                flux.append(float(parts[2]))

    return np.array(freq), np.array(flux)


def reflectance(inc_file, refl_file):
    f_inc, inc = load_flux(inc_file, "flux1:,")
    if len(f_inc) == 0:
        f_inc, inc = load_flux(inc_file, "flux2:,")

    f_refl, refl = load_flux(refl_file, "flux2:,")
    if len(f_refl) == 0:
        f_refl, refl = load_flux(refl_file, "flux1:,")

    if len(f_inc) != len(f_refl):
        raise ValueError(
            f"Length mismatch:\n"
            f"{inc_file}: {len(f_inc)} points\n"
            f"{refl_file}: {len(f_refl)} points"
        )

    if not np.allclose(f_inc, f_refl):
        raise ValueError(f"Frequency mismatch: {inc_file}, {refl_file}")

    wavelength = 1.0 / f_inc
    R = refl / np.abs(inc)
    
    return wavelength, R


# ----------------------------
# Load completed runs
# ----------------------------
wl_h_ex, R_h_ex = reflectance("inc_ex_h.txt", "refl_ex_h.txt")
wl_h_ey, R_h_ey = reflectance("inc_ey_h.txt", "refl_ey_h.txt")

wl_v_ex, R_v_ex = reflectance("inc_ex_v.txt", "refl_ex_v.txt")
wl_v_ey, R_v_ey = reflectance("inc_ey_v.txt", "refl_ey_v.txt")


# ----------------------------
# Horizontal only
# ----------------------------
plt.figure(figsize=(7, 4.5))
plt.plot(wl_h_ex, R_h_ex, label="Horizontal bar, Ex")
plt.plot(wl_h_ey, R_h_ey, label="Horizontal bar, Ey")
plt.xlabel("Wavelength (µm)")
plt.ylabel("Reflectance")
plt.title("Centered Horizontal Bar Only")
plt.legend()
plt.tight_layout()
plt.savefig("centered_horizontal_reflectance.png", dpi=300)
plt.show()


# ----------------------------
# Vertical only
# ----------------------------
plt.figure(figsize=(7, 4.5))
plt.plot(wl_v_ex, R_v_ex, label="Vertical bar, Ex")
plt.plot(wl_v_ey, R_v_ey, label="Vertical bar, Ey")
plt.xlabel("Wavelength (µm)")
plt.ylabel("Reflectance")
plt.title("Centered Vertical Bar Only")
plt.legend()
plt.tight_layout()
plt.savefig("centered_vertical_reflectance.png", dpi=300)
plt.show()


# ----------------------------
# All single-bar curves
# ----------------------------
plt.figure(figsize=(8, 5))
plt.plot(wl_h_ex, R_h_ex, label="Horizontal, Ex")
plt.plot(wl_h_ey, R_h_ey, label="Horizontal, Ey")
plt.plot(wl_v_ex, R_v_ex, label="Vertical, Ex")
plt.plot(wl_v_ey, R_v_ey, label="Vertical, Ey")
plt.xlabel("Wavelength (µm)")
plt.ylabel("Reflectance")
plt.title("Centered Single-Bar Comparison")
plt.legend()
plt.tight_layout()
plt.savefig("centered_single_all.png", dpi=300)
plt.show()


# ----------------------------
# Polarization-aligned comparison
# ----------------------------
plt.figure(figsize=(7, 4.5))
plt.plot(wl_h_ex, R_h_ex, label="Horizontal bar, Ex")
plt.plot(wl_v_ey, R_v_ey, label="Vertical bar, Ey")
plt.xlabel("Wavelength (µm)")
plt.ylabel("Reflectance")
plt.title("Polarization-Aligned Single Bars")
plt.legend()
plt.tight_layout()
plt.savefig("centered_aligned_comparison.png", dpi=300)
plt.show()
