import numpy as np
import matplotlib.pyplot as plt

def load_flux(filename):
    freqs = []
    flux = []

    with open(filename, "r") as f:
        for line in f:
            if line.startswith("flux1"):
                parts = line.strip().split(",")
                freqs.append(float(parts[1]))
                flux.append(float(parts[2]))

    freqs = np.array(freqs)
    flux = np.array(flux)

    wavelengths = 1 / freqs

    # sort wavelength increasing
    idx = np.argsort(wavelengths)
    return wavelengths[idx], flux[idx]


w_ex, F_ex = load_flux("ex_output.txt")
w_ey, F_ey = load_flux("ey_output.txt")

plt.figure(figsize=(7,5))
plt.plot(w_ex, -F_ex, label="Ex", linewidth=2)
plt.plot(w_ey, -F_ey, label="Ey", linewidth=2)

plt.xlabel("Wavelength (um)")
plt.ylabel("Reflected Flux (arb.)")
plt.title("Raw Reflected Flux Asymmetric Gold Bars")
plt.legend()
plt.tight_layout()
plt.savefig("excitonspec.png", dpi=300)
plt.show()
