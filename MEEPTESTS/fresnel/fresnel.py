import numpy as np
import matplotlib.pyplot as plt

# Load data using numpy.loadtxt (equivalent to dlmread)
f0 = np.loadtxt("flux0.txt", delimiter=",")
f = np.loadtxt("flux.txt", delimiter=",")

# Frequency is in the first column, flux in the second
# wvls = 1 / frequency
wvls = 1.0 / f0[:, 0]

# Calculate MEEP reflectance
# Note: Ensure the indices match your data structure. 
# Usually, it's (reflected_flux / incident_flux)
R_meep = f0[:,1]
#R_meep = -f[:, 1] / f0[:, 1]

# Define the Sellmeier equation for Quartz
def eps_quartz(l):
    return (1 + (0.6961663 * l**2) / (l**2 - 0.0684043**2) + 
            (0.4079426 * l**2) / (l**2 - 0.1162414**2) + 
            (0.8974794 * l**2) / (l**2 - 9.896161**2))

# Analytic Fresnel reflectance at normal incidence
def R_fresnel(l):
    n = np.sqrt(eps_quartz(l))
    return np.abs((1 - n) / (1 + n))**2

# Plotting
plt.figure(figsize=(8, 6))
plt.plot(wvls, R_meep, 'bo-', label="meep")
plt.plot(wvls, R_fresnel(wvls), 'rs-', label="analytic")

plt.xlabel("wavelength (μm)")
plt.ylabel("reflectance")
plt.legend()
plt.grid(True)
plt.show()
