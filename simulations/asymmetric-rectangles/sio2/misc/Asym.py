import numpy as np
import matplotlib.pyplot as plt

freq=[]
inc_flux=[]
freq_x = []
refl_flux_x = []
freq_y = []
refl_flux_y = []

# Inc run
with open("inc_ey.log","r") as l:
    for lineinc in l:
        lineinc=lineinc.strip()
        if lineinc.startswith("flux1:,"):
            partinc=lineinc.split(",",maxsplit=2)
            freq.append(float(partinc[1]))
            inc_flux.append(float(partinc[2]))

freq = np.array(freq)
inc_flux = np.array(inc_flux)
inc_flux = -inc_flux
wavelength_inc = 1/freq
print(inc_flux)
# Ex run
with open("refl_ey.log","r") as f:
    for line in f:
        line=line.strip()
        if line.startswith("flux1:,"):
            parts=line.split(",",maxsplit=2)
            freq_x.append(float(parts[1]))
            refl_flux_x.append(float(parts[2]))

freq_x = np.array(freq_x)
refl_flux_x = np.array(refl_flux_x[:120]) # skip first 121 points to match inc_flux size
R_x = -refl_flux_x/inc_flux
wavelength = 1/freq_x
R_x[R_x < 1e-5] = np.nan

#Ey run
with open("Ey_1400s.log","r") as f2:
    for lines in f2:
        lines=lines.strip()
        if lines.startswith("flux2:,"):
            parts=lines.split(",",maxsplit=2)
            freq_y.append(float(parts[1]))
            refl_flux_y.append(float(parts[2]))

freq_y = np.array(freq_y)
refl_flux_y = np.array(refl_flux_y) # skip first 121 points to match inc_flux size
R_y = -refl_flux_y/inc_flux
wavelength_y = 1/freq_y
R_y[R_y < 1e-5] = np.nan


plt.plot(wavelength, R_x,c='g', label="Ex polarization")
#plt.plot(wavelength_y, R_y,c='b', label="Ey polarization")
plt.xlabel("wavelength (µm)")
plt.ylabel("reflectance")
plt.title("Reflectance Spectrum")
plt.legend()

plt.savefig("new.png",dpi=300)
plt.show()