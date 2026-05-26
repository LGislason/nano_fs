import numpy as np
import matplotlib.pyplot as plt

f0 = np.loadtxt("flux0.txt", delimiter=",")
f = np.loadtxt("flux.txt", delimiter=",")

wvls = 1.0 / f0[:, 0]

# Reflected / incident
R_meep = -f[:, 1] / f0[:, 1]

# Constant-epsilon quartz/SiO2 model used in the .ctl
n = 1.45
R_analytic = ((1 - n) / (1 + n))**2

plt.figure(figsize=(8, 6))
plt.plot(wvls, R_meep, "bo-", label="Meep")
plt.axhline(R_analytic, color="r", linestyle="--",
            label=f"Fresnel analytic eps = {R_analytic:.4f}")

plt.xlabel("wavelength (μm)")
plt.ylabel("reflectance")
plt.legend()
plt.grid(True)
plt.show()
x=np.average(R_meep)
err=abs(x-R_analytic)/R_analytic
print(err)