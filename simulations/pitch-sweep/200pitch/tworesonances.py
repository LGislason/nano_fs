import numpy as np
import matplotlib.pyplot as plt
import glob
import re

# --------- USER SETTINGS ----------
# Choose a window where the resonances you care about live.
# Based on your plot, these are reasonable starting values.
WL_MIN = 600
WL_MAX = 950

# Minimum wavelength separation between the two picked dips (nm)
MIN_SEP_NM = 25

# Set True to lightly smooth reflectivity to avoid picking tiny ripples
SMOOTH = True
SMOOTH_WINDOW = 9  # must be odd
# ----------------------------------

def get_pitch_from_filename(fname: str) -> int:
    m = re.search(r"ref_(\d+)\.csv$", fname)
    if not m:
        raise ValueError(f"Could not parse pitch from filename: {fname}")
    return int(m.group(1))

def load_two_col_csv(path: str):
    data = np.loadtxt(path, delimiter=",")
    wl = data[:, 0]
    y = data[:, 1]
    idx = np.argsort(wl)
    return wl[idx], y[idx]

def moving_average(y, w):
    if w <= 1:
        return y
    if w % 2 == 0:
        raise ValueError("SMOOTH_WINDOW must be odd.")
    kernel = np.ones(w) / w
    return np.convolve(y, kernel, mode="same")

def find_local_minima(wl, y):
    # local minima: y[i] < y[i-1] and y[i] < y[i+1]
    # returns indices
    return np.where((y[1:-1] < y[:-2]) & (y[1:-1] < y[2:]))[0] + 1

ref_files = sorted(glob.glob("ref_*.csv"))
if not ref_files:
    raise SystemExit("No files matching ref_*.csv found.")

results = []

for ref_file in ref_files:
    pitch = get_pitch_from_filename(ref_file)
    obj_file = f"obj_{pitch}.csv"

    wl_ref, flux_ref = load_two_col_csv(ref_file)
    wl_obj, flux_obj = load_two_col_csv(obj_file)

    # Interpolate object flux onto reference wavelength grid
    flux_obj_interp = np.interp(wl_ref, wl_obj, flux_obj)

    # Reflectivity
    R = -flux_obj_interp / flux_ref

    # Window to avoid edges and focus on the resonances of interest
    mask = (wl_ref >= WL_MIN) & (wl_ref <= WL_MAX)
    wl = wl_ref[mask]
    Rw = R[mask]

    if SMOOTH:
        Rw_use = moving_average(Rw, SMOOTH_WINDOW)
    else:
        Rw_use = Rw

    # We want dips, so work with minima of R
    mins = find_local_minima(wl, Rw_use)

    if len(mins) < 2:
        # Fall back: just take two lowest points if local minima are scarce
        cand = np.argsort(Rw_use)[:2]
        mins = cand

    # Sort candidate minima by depth (lowest reflectivity first)
    mins_sorted = mins[np.argsort(Rw_use[mins])]

    # Pick first dip
    i1 = mins_sorted[0]
    wl1 = wl[i1]
    depth1 = Rw_use[i1]

    # Pick second dip that is sufficiently separated
    wl2 = np.nan
    depth2 = np.nan
    for idx in mins_sorted[1:]:
        if abs(wl[idx] - wl1) >= MIN_SEP_NM:
            wl2 = wl[idx]
            depth2 = Rw_use[idx]
            break

    results.append([pitch, wl1, depth1, wl2, depth2])

# Sort by pitch
results = np.array(results, dtype=float)
order = np.argsort(results[:, 0])
results = results[order]

# Save CSV
header = "pitch_nm,res1_wavelength_nm,res1_reflectivity,res2_wavelength_nm,res2_reflectivity"
np.savetxt("pitch_vs_resonance_two_branches.csv", results, delimiter=",", header=header, comments="")
print("Wrote pitch_vs_resonance_two_branches.csv")

# Plot
p = results[:, 0]
wl1 = results[:, 1]
wl2 = results[:, 3]

plt.figure()
plt.plot(p, wl1, marker="o", label="Dip 1 (global minimum)")
plt.plot(p, wl2, marker="o", label="Dip 2 (local minimum)")
plt.xlabel("Pitch (nm)")
plt.ylabel("Resonance Wavelength (nm)")
plt.title("Pitch vs Resonance (200 nm PPMA)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("pitch_vs_resonance_two_branches.png", dpi=300)
plt.show()

