#!/usr/bin/env python3
"""Plot the asymmetric nanobar metasurface sweep.

Reads the output of run_asym_sweep.sh. Each shape{S}_pol{P}.txt file has a
labelled reflected block and a transmitted block:

    --- BEGIN REFLECTED FLUX ---
    flux1:, freq, value
    ...
    --- END REFLECTED FLUX ---
    --- BEGIN TRANSMITTED FLUX ---
    ...

The reflected flux is DIFFERENTIAL (metasurface minus bare substrate), so it is
the metal's polarization-dependent contribution rather than an absolute R. Since
the same source is used for every run, the raw differential flux is directly
comparable across shapes and polarizations. Two figures are produced: Ex vs Ey
for the asymmetric cell, and asymmetric vs single-bar controls.

Usage:
    python analyze_asym_reflectance.py <results_dir>
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SHAPE_NAME = {0: "horizontal bar", 1: "vertical bar", 2: "asymmetric pair"}
POL_NAME = {0: "Ex", 1: "Ey"}
CBLUE, CGOLD = "#1C7293", "#C08028"


def parse_blocks(path: Path):
    """Return {'reflected': (wl, flux), 'transmitted': (wl, flux)} from a run file."""
    text = path.read_text()
    out = {}
    for key, tag in (("reflected", "REFLECTED"), ("transmitted", "TRANSMITTED")):
        m = re.search(rf"BEGIN {tag} FLUX(.*?)END {tag} FLUX", text, re.S)
        if not m:
            continue
        wl, fx = [], []
        for line in m.group(1).splitlines():
            s = line.strip()
            if not s.startswith("flux") or ":" not in s:
                continue
            v = []
            for x in next(csv.reader([s.split(":", 1)[1]])):
                try:
                    v.append(float(x.strip()))
                except ValueError:
                    pass
            if len(v) >= 2:
                wl.append(1.0 / v[0]); fx.append(v[1])
        if wl:
            o = np.argsort(wl)
            out[key] = (np.array(wl)[o], np.array(fx)[o])
    return out


def load(results: Path, shape: int, pol: int, which: str):
    return parse_blocks(results / f"shape{shape}_pol{pol}.txt")[which]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("results_dir", type=Path)
    ap.add_argument("--which", choices=["reflected", "transmitted"], default="reflected")
    ap.add_argument("--band", type=float, nargs=2, default=[0.55, 1.15],
                    help="reliable wavelength window (µm); source power is low near the edges")
    args = ap.parse_args()
    d = args.results_dir
    lo, hi = args.band
    ylab = ("differential reflected flux (metal − substrate)"
            if args.which == "reflected" else "transmitted flux")

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "savefig.dpi": 300, "savefig.bbox": "tight"})

    # ---- Figure 1: Ex vs Ey for the asymmetric cell ----
    fig, ax = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
    for pol, c in ((0, CBLUE), (1, CGOLD)):
        try:
            wl, f = load(d, 2, pol, args.which)
        except (FileNotFoundError, KeyError):
            continue
        m = (wl >= lo) & (wl <= hi)
        ax.plot(wl[m], f[m], color=c, lw=2, label=POL_NAME[pol])
    ax.axhline(0, color="0.7", lw=0.8)
    ax.set_xlabel("wavelength (µm)"); ax.set_ylabel(ylab)
    ax.set_title("Asymmetric metasurface — polarization dependence")
    ax.grid(True, color="#E2E8F0"); ax.legend(frameon=False)
    fig.savefig(d / f"{args.which}_asym_ExEy.png"); plt.close(fig)

    # ---- Figure 2: asymmetric vs single-bar controls ----
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True, sharey=True)
    for pol, ax in zip((0, 1), axes):
        for shape, c, ls in ((0, "#6D2E46", "--"), (1, "#2C5F2D", "--"), (2, CBLUE, "-")):
            try:
                wl, f = load(d, shape, pol, args.which)
            except (FileNotFoundError, KeyError):
                continue
            m = (wl >= lo) & (wl <= hi)
            ax.plot(wl[m], f[m], color=c, ls=ls, lw=2 if shape == 2 else 1.4,
                    label=SHAPE_NAME[shape])
        ax.axhline(0, color="0.7", lw=0.8)
        ax.set_xlabel("wavelength (µm)"); ax.set_title(f"{POL_NAME[pol]} excitation")
        ax.grid(True, color="#E2E8F0"); ax.legend(frameon=False, fontsize=9)
    axes[0].set_ylabel(ylab)
    fig.suptitle("Asymmetric cell vs single-bar controls")
    fig.savefig(d / f"{args.which}_asym_vs_controls.png"); plt.close(fig)

    # ---- text summary ----
    print(f"Input: {d}   quantity: {args.which}   band: {lo}-{hi} µm")
    print("shape,pol,peak,lambda_peak_um")
    for shape in (0, 1, 2):
        for pol in (0, 1):
            try:
                wl, f = load(d, shape, pol, args.which)
            except (FileNotFoundError, KeyError):
                continue
            m = (wl >= lo) & (wl <= hi)
            wlm, fm = wl[m], f[m]
            i = int(np.argmax(np.abs(fm)))
            print(f"{shape},{pol},{fm[i]:.4g},{wlm[i]:.3f}")
    print(f"\nWrote {d / (args.which + '_asym_ExEy.png')}")
    print(f"Wrote {d / (args.which + '_asym_vs_controls.png')}")


if __name__ == "__main__":
    main()
