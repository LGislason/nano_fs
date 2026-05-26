import numpy as np
import matplotlib.pyplot as plt
import glob
import re

def get_pitch_from_filename(fname: str) -> int:
    m = re.search(r"ref_(\d+)\.csv$", fname)
    if not m:
        raise ValueError(f"Could not parse pitch from filename: {fname}")
    return int(m.group(1))

def load_two_col_csv(path: str):
    data = np.loadtxt(path, delimiter=",")
    wl = data[:, 0]
    y = data[:, 1]
    # sort by wavelength
    idx = np.argsort(wl)
    return wl[idx], y[idx]

# Find all pitches from ref_*.csv
ref_files = sorted(glob.glob("ref_*.csv"))
if not ref_files:
    raise SystemExit("No files matching ref_*.csv found in this directory.")

pitches = []
res_wls = []

for ref_file in ref_files:
    pitch = get_pitch_from_filename(ref_file)
    obj_file = f"obj_{pitch}.csv"

    wl_ref, flux_ref = load_two_col_csv(ref_file)
    wl_obj, flux_obj = load_two_col_csv(obj_file)

    # Interpolate object flux onto reference wavelength grid
    flux_obj_interp = np.interp(wl_ref, wl_obj, flux_obj)

    # Reflectivity
    R = -flux_obj_interp / flux_ref

    # Optional: trim edge noise by restricting to central wavelength window
    # Uncomment and adjust if needed
    # mask = (wl_ref > 500) & (wl_ref < 1000)
    # wl_use = wl_ref[mask]
    # R_use = R[mask]
    # Otherwise use all:
    wl_use = wl_ref
    R_use = R

    # Resonance definition:
    # Peak reflectivity wavelength
    res_idx = np.argmax(R_use)

    # If your resonance is a dip, use this instead:
    # res_idx = np.argmin(R_use)

    res_wl = wl_use[res_idx]

    pitches.append(pitch)
    res_wls.append(res_wl)

# Sort by pitch for nice plotting
pitches = np.array(pitches)
res_wls = np.array(res_wls)
order = np.argsort(pitches)
pitches = pitches[order]
res_wls = res_wls[order]

# Save results
out = np.column_stack([pitches, res_wls])
np.savetxt("pitch_vs_resonance.csv", out, delimiter=",", header="pitch_nm,resonance_wavelength_nm", comments="")
print("Wrote pitch_vs_resonance.csv")

# Plot
plt.figure()
plt.plot(pitches, res_wls, marker="o")
plt.xlabel("Pitch (nm)")
plt.ylabel("Resonance wavelength (nm)")
plt.title("Pitch vs Resonance (100 nm pitch)")
plt.grid(True)
plt.tight_layout()
plt.savefig("pitch_vs_resonance.png", dpi=300)
plt.show()

