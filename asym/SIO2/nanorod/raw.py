import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def read_meep_flux_txt(filename):
    """
    Reads raw Meep text output and extracts lines beginning with 'flux1:'.

    Expected line format:
    flux1:, freq, flux_col_1, flux_col_2, ...

    Returns:
        freqs: frequency array
        fluxes: 2D array with shape (n_points, n_flux_columns)
    """
    freqs = []
    flux_rows = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            if not line.startswith("flux1:"):
                continue

            parts = line.replace("flux1:", "").replace(",", " ").split()
            nums = [float(p) for p in parts]

            freq = nums[0]
            flux_vals = nums[1:]

            freqs.append(freq)
            flux_rows.append(flux_vals)

    freqs = np.array(freqs)
    fluxes = np.array(flux_rows)

    return freqs, fluxes


def plot_raw_file(filename, outdir="raw_plots"):
    filename = Path(filename)
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)

    freqs, fluxes = read_meep_flux_txt(filename)

    if len(freqs) == 0:
        print(f"No flux1 lines found in {filename}")
        return

    wavelengths = 1.0 / freqs

    print(f"\nFile: {filename}")
    print(f"Number of flux points: {len(freqs)}")
    print(f"Number of flux columns: {fluxes.shape[1]}")
    print(f"Frequency range: {freqs.min():.6g} to {freqs.max():.6g}")
    print(f"Wavelength range: {wavelengths.min():.6g} to {wavelengths.max():.6g} µm")

    for i in range(fluxes.shape[1]):
        col = fluxes[:, i]
        print(f"Column {i + 1} range: {col.min():.6e} to {col.max():.6e}")

    # Sort by wavelength so plots go left-to-right naturally
    sort_idx = np.argsort(wavelengths)
    wavelengths = wavelengths[sort_idx]
    fluxes = fluxes[sort_idx, :]

    plt.figure(figsize=(9, 5))

    for i in range(fluxes.shape[1]):
        plt.plot(wavelengths, fluxes[:, i], label=f"Raw flux column {i + 1}")

    plt.xlabel("Wavelength (µm)")
    plt.ylabel("Raw Meep flux")
    plt.title(f"Raw Meep flux: {filename.name}")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    save_path = outdir / f"{filename.stem}_raw_flux.png"
    plt.savefig(save_path, dpi=300)
    plt.show()

    print(f"Saved: {save_path}")


if __name__ == "__main__":
    files = [
        "inc_ex_dimer.txt",
        "rod_ey_single.txt",
        # add more files here:
        # "rod_ex_single.txt",
        # "rod_ex_dimer.txt",
        # "rod_ey_dimer.txt",
    ]

    for file in files:
        if Path(file).exists():
            plot_raw_file(file)
        else:
            print(f"Missing file: {file}")