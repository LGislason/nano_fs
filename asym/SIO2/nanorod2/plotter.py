import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


# ============================================================
# USER FILES
# Edit these names to match your actual files
# ============================================================

CASES = {
    "Single nanorod, Ex": {
        "incident": "inc_ex_single.txt",
        "rod": "rod_ex_single.txt",
    },
    "Single nanorod, Ey": {
        "incident": "inc_ey_single.txt",
        "rod": "rod_ey_single.txt",
    },
    "Nanorod dimer, Ex": {
        "incident": "inc_ex_dimer.txt",
        "rod": "rod_ex_dimer.txt",
    },
    "Nanorod dimer, Ey": {
        "incident": "inc_ey_dimer.txt",
        "rod": "rod_ey_dimer.txt",
    },
}


OUTDIR = Path("nanorod_plots")
OUTDIR.mkdir(exist_ok=True)


# ============================================================
# FILE LOADER
# ============================================================

def load_flux_file(filename):
    """
    Reads a Meep flux text file.

    Expected numeric data after removing possible Meep text/log lines:

        frequency   reflection_flux   transmission_flux

    Also supports Meep-style lines like:

        flux1:, frequency, reflection_flux, transmission_flux

    It skips non-numeric lines like:
        Using MPI...
        meep: ...
        elapsed run time...
    """

    rows = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            # If line is like "flux1:, 1.0, -0.003, -0.003"
            # remove the "flux1:" part.
            if ":" in line:
                line = line.split(":", 1)[1]

            # Allow comma-separated or whitespace-separated data
            parts = line.replace(",", " ").split()

            if len(parts) < 3:
                continue

            try:
                nums = [float(x) for x in parts[:3]]
                rows.append(nums)
            except ValueError:
                continue

    if len(rows) == 0:
        raise ValueError(f"No numeric flux rows found in {filename}")

    data = np.array(rows, dtype=float)

    freq = data[:, 0]
    refl_flux = data[:, 1]
    trans_flux = data[:, 2]
    wavelength = 1.0 / freq

    # Put wavelength in increasing order for plotting
    order = np.argsort(wavelength)

    return {
        "freq": freq[order],
        "wavelength": wavelength[order],
        "refl_flux": refl_flux[order],
        "trans_flux": trans_flux[order],
    }


# ============================================================
# CALCULATIONS
# ============================================================

def calculate_raw_rt(inc, rod):
    """
    Calculates R/T/A only where the incident normalization flux is strong enough.

    Column convention:
        refl_flux = reflection monitor
        trans_flux = transmission monitor

    Sign convention from your sign-check:
        R = -rod_refl / inc_refl
        T =  rod_trans / inc_trans
    """

    wavelength = rod["wavelength"]

    inc_refl = inc["refl_flux"]
    inc_trans = inc["trans_flux"]

    rod_refl = rod["refl_flux"]
    rod_trans = rod["trans_flux"]

    # --------------------------------------------------------
    # Cutoff weak incident-flux regions
    # --------------------------------------------------------
    # Keep only points where incident flux is at least this fraction
    # of the maximum incident flux.
    cutoff_fraction = 5e-5

    refl_cutoff = cutoff_fraction * np.nanmax(np.abs(inc_refl))
    trans_cutoff = cutoff_fraction * np.nanmax(np.abs(inc_trans))

    valid = np.abs(inc_refl) > refl_cutoff
    valid &= np.abs(inc_trans) > trans_cutoff

    wavelength = wavelength[valid]
    inc_refl = inc_refl[valid]
    inc_trans = inc_trans[valid]
    rod_refl = rod_refl[valid]
    rod_trans = rod_trans[valid]

    # --------------------------------------------------------
    # Normalize
    # --------------------------------------------------------
    R_plus = rod_refl / inc_refl
    R_minus = -rod_refl / inc_refl

    T = rod_trans / inc_trans

    A_using_R_plus = 1.0 - R_plus - T
    A_using_R_minus = 1.0 - R_minus - T

    return {
        "wavelength": wavelength,
        "R_plus": R_plus,
        "R_minus": R_minus,
        "T": T,
        "A_using_R_plus": A_using_R_plus,
        "A_using_R_minus": A_using_R_minus,
        "inc_refl": inc_refl,
        "inc_trans": inc_trans,
        "rod_refl": rod_refl,
        "rod_trans": rod_trans,
        "refl_cutoff": refl_cutoff,
        "trans_cutoff": trans_cutoff,
        "points_kept": len(wavelength),
        "points_total": len(valid),
    }


# ============================================================
# PLOTTING FUNCTIONS
# ============================================================

def plot_raw_flux(case_name, inc, rod):
    wl = inc["wavelength"]

    plt.figure(figsize=(10, 6))
    plt.plot(wl, inc["refl_flux"], label="Incident reflection monitor")
    plt.plot(wl, inc["trans_flux"], label="Incident transmission monitor")
    plt.plot(wl, rod["refl_flux"], label="Rod reflection monitor")
    plt.plot(wl, rod["trans_flux"], label="Rod transmission monitor")

    plt.axhline(0, linewidth=1)
    plt.xlabel("Wavelength (µm)")
    plt.ylabel("Raw Meep flux")
    plt.title(f"{case_name}: raw monitor flux")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    safe_name = case_name.replace(",", "").replace(" ", "_").lower()
    plt.savefig(OUTDIR / f"{safe_name}_raw_flux.png", dpi=300)
    plt.close()


def plot_reflection_sign_check(case_name, result):
    wl = result["wavelength"]

    plt.figure(figsize=(10, 6))
    plt.plot(wl, result["R_plus"], label="R = rod_refl / inc_refl")
    plt.plot(wl, result["R_minus"], label="R = -rod_refl / inc_refl")

    plt.axhline(0, linewidth=1)
    plt.xlabel("Wavelength (µm)")
    plt.ylabel("Reflectance candidate")
    plt.title(f"{case_name}: reflectance sign check")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    safe_name = case_name.replace(",", "").replace(" ", "_").lower()
    plt.savefig(OUTDIR / f"{safe_name}_reflectance_sign_check.png", dpi=300)
    plt.close()


def plot_rta_case(case_name, result):
    wl = result["wavelength"]

    # Based on sign check, this is the physical candidate
    R = result["R_minus"]
    T = result["T"]
    A = result["A_using_R_minus"]

    plt.figure(figsize=(10, 6))
    plt.plot(wl, R, label="Reflectance")
    plt.plot(wl, T, label="Transmission")
    plt.plot(wl, A, label="Absorption / loss candidate")

    plt.axhline(0, linewidth=1)
    plt.xlabel("Wavelength (µm)")
    plt.ylabel("Fraction")
    plt.title(case_name)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    safe_name = case_name.replace(",", "").replace(" ", "_").lower()
    plt.savefig(OUTDIR / f"{safe_name}_rta.png", dpi=300)
    plt.close()


def plot_comparison(results, quantity, ylabel, title, filename):
    plt.figure(figsize=(10, 6))

    for case_name, result in results.items():
        plt.plot(result["wavelength"], result[quantity], label=case_name)

    plt.axhline(0, linewidth=1)
    plt.xlabel("Wavelength (µm)")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / filename, dpi=300)
    plt.close()


# ============================================================
# MAIN
# ============================================================

all_results = {}

for case_name, files in CASES.items():
    print("=" * 70)
    print(case_name)
    print("=" * 70)

    inc_file = files["incident"]
    rod_file = files["rod"]

    inc = load_flux_file(inc_file)
    rod = load_flux_file(rod_file)

    result = calculate_raw_rt(inc, rod)
    all_results[case_name] = result

    print(f"Incident file: {inc_file}")
    print(f"Rod file:      {rod_file}")
    print(f"Wavelength range: {result['wavelength'].min():.4f} to {result['wavelength'].max():.4f} µm")
    print()
    print("Raw flux ranges:")
    print(f"  inc reflection:   {result['inc_refl'].min():.6e} to {result['inc_refl'].max():.6e}")
    print(f"  inc transmission: {result['inc_trans'].min():.6e} to {result['inc_trans'].max():.6e}")
    print(f"  rod reflection:   {result['rod_refl'].min():.6e} to {result['rod_refl'].max():.6e}")
    print(f"  rod transmission: {result['rod_trans'].min():.6e} to {result['rod_trans'].max():.6e}")
    print()
    print("Candidate ranges:")
    print(f"  R_plus:           {result['R_plus'].min():.6e} to {result['R_plus'].max():.6e}")
    print(f"  R_minus:          {result['R_minus'].min():.6e} to {result['R_minus'].max():.6e}")
    print(f"  T:                {result['T'].min():.6e} to {result['T'].max():.6e}")
    print(f"  A using R_minus:  {result['A_using_R_minus'].min():.6e} to {result['A_using_R_minus'].max():.6e}")
    print()

    plot_raw_flux(case_name, inc, rod)
    plot_reflection_sign_check(case_name, result)
    plot_rta_case(case_name, result)


# ============================================================
# COMPARISON PLOTS
# ============================================================

plot_comparison(
    all_results,
    quantity="R_minus",
    ylabel="Reflectance",
    title="Nanorod validation: Reflectance",
    filename="comparison_reflectance.png",
)

plot_comparison(
    all_results,
    quantity="T",
    ylabel="Transmission",
    title="Nanorod validation: Transmission",
    filename="comparison_transmission.png",
)

plot_comparison(
    all_results,
    quantity="A_using_R_minus",
    ylabel="Absorption / loss candidate",
    title="Nanorod validation: Absorption / loss candidate",
    filename="comparison_absorption_loss_candidate.png",
)

print("Saved plots to:", OUTDIR.resolve())