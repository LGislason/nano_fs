"""
Diameter Sweep

Plots reflectivity vs wavelength for all swept diameters,
and also shows individual plots per diameter.

Leif Gislason
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import os

from utils import load_two_col_csv


# -------------------------
# Collect files
# -------------------------

ref_files = sorted(glob.glob("ref_*.csv"))

if len(ref_files) == 0:
    raise RuntimeError("No reference files found.")


# Combined figure
plt.figure()
plt.title("Radius Sweep Reflectivity (All)")
plt.xlabel("Wavelength")
plt.ylabel("Reflectivity")
plt.grid(True)


# -------------------------
# Main loop
# -------------------------

for ref_file in ref_files:
    diameter = ref_file.split("_")[1].split(".")[0]
    obj_file = f"obj_{diameter}.csv"

    if not os.path.exists(obj_file):
        print(f"Skipping missing {obj_file}")
        continue

    wl_ref, flux_ref = load_two_col_csv(ref_file)
    wl_obj, flux_obj = load_two_col_csv(obj_file)

    # Interpolate object onto reference wavelength grid
    flux_obj_interp = np.interp(wl_ref, wl_obj, flux_obj)

    reflectivity = flux_obj_interp / flux_ref

    # ---- Individual plot ----
    plt.figure()
    plt.plot(wl_ref, -reflectivity)
    plt.xlabel("Wavelength")
    plt.ylabel("Reflectivity")
    plt.title(f"Radius {diameter} nm")
    plt.grid(True)

    plt.savefig(f"reflectivity_{diameter}nm.png", dpi=200)
    plt.close()


    # ---- Add to combined plot ----
    plt.figure(1)
    plt.plot(wl_ref, -reflectivity, label=f"{diameter} nm")


# -------------------------
# Finish combined plot
# -------------------------

plt.figure(1)
plt.legend()
plt.tight_layout()
plt.savefig("reflectivity_all_diameters.png", dpi=200)
plt.show()


