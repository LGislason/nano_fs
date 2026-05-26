import numpy as np
import matplotlib.pyplot as plt
import glob

# Find all reference CSVs
ref_files = sorted(glob.glob("ref_*.csv"))

plt.figure()

for ref_file in ref_files:
    pitch = ref_file.split("_")[1].split(".")[0]
    obj_file = f"obj_{pitch}.csv"

    ref = np.loadtxt(ref_file, delimiter=",")
    obj = np.loadtxt(obj_file, delimiter=",")

    wl_ref = ref[:, 0]
    flux_ref = ref[:, 1]

    wl_obj = obj[:, 0]
    flux_obj = obj[:, 1]

    # Sort for safety
    iref = np.argsort(wl_ref)
    iobj = np.argsort(wl_obj)

    wl_ref = wl_ref[iref]
    flux_ref = flux_ref[iref]

    wl_obj = wl_obj[iobj]
    flux_obj = flux_obj[iobj]

    # Interpolate object flux onto reference grid
    flux_obj_interp = np.interp(wl_ref, wl_obj, flux_obj)

    # Reflectivity
    R = -flux_obj_interp / flux_ref

    plt.plot(wl_ref, R, label=f"pitch = {pitch} nm")

plt.xlabel("Wavelength (nm)")
plt.ylabel("Reflectivity")
plt.title("Reflectivity vs Wavelength (200 nm PPMA)")
plt.legend()
plt.grid(True)
plt.tight_layout()

plt.savefig("reflectivity_all_pitches.png", dpi=300)
plt.show()

