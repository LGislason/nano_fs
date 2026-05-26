import numpy as np
import matplotlib.pyplot as plt


def load_flux(filename, tag="flux1:,"):
    freq = []
    flux = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(tag):
                parts = line.split(",", maxsplit=2)
                freq.append(float(parts[1]))
                flux.append(float(parts[2]))

    return np.array(freq), np.array(flux)


def compute_reflectance(inc_file, refl_file):
    freq_inc, inc_flux = load_flux(inc_file, "flux1:,")
    freq_refl, refl_flux = load_flux(refl_file, "flux1:,")

    print(f"{inc_file}: {len(freq_inc)} points")
    print(f"{refl_file}: {len(freq_refl)} points")

    if len(freq_inc) == 0 or len(freq_refl) == 0:
        raise ValueError("No flux data found in one of the log files")

    if len(freq_inc) != len(freq_refl):
        raise ValueError("Length mismatch between incident and reflected runs")

    if not np.allclose(freq_inc, freq_refl, rtol=1e-10, atol=1e-12):
        raise ValueError("Frequency arrays do not match")

    wavelength = 1.0 / freq_inc

    inc_flux = (inc_flux)
    refl_flux = (refl_flux)

    R = -refl_flux / inc_flux

    return wavelength, R


wl_ex, R_ex = compute_reflectance("inc_ex.log", "refl_ex.log")
wl_ey, R_ey = compute_reflectance("inc_ey.log", "refl_ey.log")

i_peak = np.nanargmax(R_ex)   # or R_ey
peak_wavelength = wl_ex[i_peak]
peak_value = R_ex[i_peak]

print("Peak wavelength =", peak_wavelength, "um")
print("Peak reflectance =", peak_value)
print("Ex max R =", np.max(R_ex))
print("Ey max R =", np.max(R_ey))

plt.figure(figsize=(7, 4.5))
plt.plot(wl_ex, R_ex, label="Ex Polarization")
plt.plot(wl_ey, R_ey, label="Ey Polarization")
plt.xlabel("Wavelength (µm)")
plt.ylabel("Reflectance")
plt.title("Asymmetric Rectangle Reflectance")
plt.legend()
plt.tight_layout()
plt.savefig("100res_Sweep2.png",dpi=300)
plt.show()