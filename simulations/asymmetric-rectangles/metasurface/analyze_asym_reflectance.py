#!/usr/bin/env python3
"""Reflectance spectra for the asymmetric nanobar metasurface sweep.

Reads the output of run_asym_sweep.sh and computes, per shape and polarization,
the reflectance R(λ) = reflected flux / incident flux (vacuum reference), plus
transmittance T and absorptance A = 1 − R − T. Produces two comparison figures:
Ex vs Ey for the asymmetric cell, and asymmetric vs single-bar controls.

Usage:
    python analyze_asym_reflectance.py <results_dir>
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SHAPE_NAME = {0: "horizontal bar", 1: "vertical bar", 2: "asymmetric pair"}
POL_NAME = {0: "Ex", 1: "Ey"}
CBLUE, CGOLD = "#1C7293", "#C08028"


def parse(path: Path):
    wl, refl, tran = [], [], []
    for raw in path.read_text().splitlines():
        s = raw.strip()
        if not s.startswith("flux") or ":" not in s:
            continue
        v = []
        for x in next(csv.reader([s.split(":", 1)[1]])):
            try:
                v.append(float(x.strip()))
            except ValueError:
                pass
        if len(v) >= 3:
            wl.append(1.0 / v[0]); refl.append(v[1]); tran.append(v[2])
    if not wl:
        raise ValueError(f"no flux rows in {path}")
    o = np.argsort(wl)
    return np.array(wl)[o], np.array(refl)[o], np.array(tran)[o]


def rta(results: Path, shape: int, pol: int):
    """Return wavelength, R, T, A for a (shape, pol) case."""
    wl, refl, tran = parse(results / f"shape{shape}_pol{pol}.txt")
    _, iref, itran = parse(results / f"reference_incident_pol_{pol}.txt")
    inc = np.abs(itran)                      # incident power (downward)
    inc[inc == 0] = np.nan
    R = refl / inc
    T = -tran / inc
    A = 1.0 - R - T
    return wl, R, T, A


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("results_dir", type=Path)
    ap.add_argument("--band", type=float, nargs=2, default=[0.55, 1.15],
                    help="reliable wavelength window (um) to display; source power "
                         "is low outside this, so R can get noisy near the edges")
    args = ap.parse_args()
    d = args.results_dir
    lo, hi = args.band

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "savefig.dpi": 300, "savefig.bbox": "tight"})

    # ---- Figure 1: Ex vs Ey for the asymmetric cell ----
    fig, ax = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
    for pol, c in ((0, CBLUE), (1, CGOLD)):
        try:
            wl, R, _, _ = rta(d, 2, pol)
        except FileNotFoundError:
            continue
        m = (wl >= lo) & (wl <= hi)
        ax.plot(wl[m], R[m], color=c, lw=2, label=f"{POL_NAME[pol]}")
    ax.set_xlabel("wavelength (µm)"); ax.set_ylabel("reflectance R")
    ax.set_title("Asymmetric metasurface — polarization dependence")
    ax.grid(True, color="#E2E8F0"); ax.legend(frameon=False)
    fig.savefig(d / "reflectance_asym_ExEy.png"); plt.close(fig)

    # ---- Figure 2: asymmetric vs single-bar controls (per polarization) ----
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True, sharey=True)
    for pol, ax in zip((0, 1), axes):
        for shape, c, ls in ((0, "#6D2E46", "--"), (1, "#2C5F2D", "--"), (2, CBLUE, "-")):
            try:
                wl, R, _, _ = rta(d, shape, pol)
            except FileNotFoundError:
                continue
            m = (wl >= lo) & (wl <= hi)
            ax.plot(wl[m], R[m], color=c, ls=ls, lw=2 if shape == 2 else 1.4,
                    label=SHAPE_NAME[shape])
        ax.set_xlabel("wavelength (µm)"); ax.set_title(f"{POL_NAME[pol]} excitation")
        ax.grid(True, color="#E2E8F0"); ax.legend(frameon=False, fontsize=9)
    axes[0].set_ylabel("reflectance R")
    fig.suptitle("Asymmetric cell vs single-bar controls")
    fig.savefig(d / "reflectance_asym_vs_controls.png"); plt.close(fig)

    # ---- text summary ----
    print(f"Input: {d}   display band: {lo}-{hi} µm")
    print("shape,pol,R_peak,lambda_peak_um,R_min,lambda_min_um")
    for shape in (0, 1, 2):
        for pol in (0, 1):
            try:
                wl, R, _, _ = rta(d, shape, pol)
            except FileNotFoundError:
                continue
            m = (wl >= lo) & (wl <= hi)
            wlm, Rm = wl[m], R[m]
            i, j = int(np.argmax(Rm)), int(np.argmin(Rm))
            print(f"{shape},{pol},{Rm[i]:.3f},{wlm[i]:.3f},{Rm[j]:.3f},{wlm[j]:.3f}")
    print(f"\nWrote {d/'reflectance_asym_ExEy.png'}")
    print(f"Wrote {d/'reflectance_asym_vs_controls.png'}")


if __name__ == "__main__":
    main()
