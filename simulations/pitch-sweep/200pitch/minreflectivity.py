import numpy as np
import matplotlib.pyplot as plt
import glob
import re

WL_MIN = 700
WL_MAX = 950

def pitch_from(fname):
    return int(re.search(r"ref_(\d+)", fname).group(1))

results = []

for ref in sorted(glob.glob("ref_*.csv")):
    pitch = pitch_from(ref)
    obj = f"obj_{pitch}.csv"

    wl_r, fr = np.loadtxt(ref, delimiter=",").T
    wl_o, fo = np.loadtxt(obj, delimiter=",").T

    fo_i = np.interp(wl_r, wl_o, fo)
    R = -fo_i / fr

    mask = (wl_r > WL_MIN) & (wl_r < WL_MAX)
    Rw = R[mask]

    results.append([pitch, np.percentile(Rw, 1)])

results = np.array(results)
results = results[results[:,0].argsort()]

np.savetxt("pitch_vs_min_reflectivity.csv", results,
           delimiter=",",
           header="pitch_nm,min_reflectivity",
           comments="")

plt.plot(results[:,0], results[:,1], marker="o")
plt.xlabel("Pitch (nm)")
plt.ylabel("Minimum reflectivity")
plt.title("Pitch vs minimum reflectivity (200 nm PMMA)")
plt.grid(True)
plt.tight_layout()
plt.savefig("pitch_vs_min_reflectivity.png", dpi=300)
plt.show()

