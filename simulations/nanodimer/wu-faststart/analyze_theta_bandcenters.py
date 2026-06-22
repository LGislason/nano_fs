#!/usr/bin/env python3
"""Track the two coupled-mode band centers vs structure angle theta.

Reads a two-polarization theta sweep (run_theta_polsum_background.sh): for each
theta it sums the Ex (pol 0) and Ey (pol 1) proxy spectra to get an
orientation-independent response in which BOTH modes are visible, finds each band
center, and plots band-center-vs-theta. Wu et al. report that the resonance
wavelengths are independent of the structure angle, so a near-flat band-center vs
theta is the clean validation statement (as opposed to the phi=0 / Ex-only sweep,
where the dominant peak appears to slide because changing theta also reorients the
rods relative to the fixed field).

No scipy dependency (parabolic sub-grid refinement of the windowed peak).
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

NAVY = "#12243B"
GOLD = "#C08028"
TEAL = "#1C7293"


def parse_flux_file(path: Path):
    """Return (wavelength_um, reflected, transmitted) arrays from a meep flux file."""
    wl, refl, tran = [], [], []
    with path.open(newline="") as fh:
        for raw in fh:
            line = raw.strip()
            if not line.startswith("flux") or ":" not in line:
                continue
            vals = []
            for v in next(csv.reader([line.split(":", 1)[1]])):
                try:
                    vals.append(float(v.strip()))
                except ValueError:
                    pass
            if len(vals) >= 3:
                f, r, t = vals[:3]
                wl.append(1.0 / f); refl.append(r); tran.append(t)
    if not wl:
        raise ValueError(f"No flux rows in {path}")
    order = np.argsort(wl)
    return np.array(wl)[order], np.array(refl)[order], np.array(tran)[order]


def proxy(ref_path: Path, case_path: Path, quantity: str):
    wl, refl, tran = parse_flux_file(case_path)
    _, _, ref_tran = parse_flux_file(ref_path)
    inc = np.abs(ref_tran)
    inc[inc == 0] = np.nan
    if quantity == "backscatter":
        return wl, refl / inc
    return wl, -tran / inc  # extinction


def band_center(wl, y, lo, hi):
    """Sub-grid peak wavelength within [lo, hi] by 3-point parabolic fit."""
    sel = np.where((wl >= lo) & (wl <= hi))[0]
    if sel.size == 0:
        return float("nan"), float("nan")
    i = sel[np.argmax(y[sel])]
    if 0 < i < len(wl) - 1:
        y0, y1, y2 = y[i - 1], y[i], y[i + 1]
        denom = (y0 - 2 * y1 + y2)
        delta = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
        delta = max(-1.0, min(1.0, delta))
        dx = (wl[i + 1] - wl[i - 1]) / 2.0
        return float(wl[i] + delta * dx), float(y1)
    return float(wl[i]), float(y[i])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--quantity", choices=["backscatter", "extinction"], default="backscatter")
    ap.add_argument("--short-window", type=float, nargs=2, default=[0.64, 0.73],
                    help="wavelength window (um) for the short band")
    ap.add_argument("--main-window", type=float, nargs=2, default=[0.73, 0.85],
                    help="wavelength window (um) for the main band")
    ap.add_argument("-o", "--output", type=Path, default=None)
    args = ap.parse_args()

    d = args.input_dir
    thetas = sorted(int(m.group(1)) for m in
                    (re.match(r"theta_(\d+)_pol_0\.txt$", p.name) for p in d.glob("theta_*_pol_0.txt")) if m)
    if not thetas:
        raise SystemExit(f"No theta_*_pol_0.txt files in {d}")

    ref0 = d / "reference_incident_pol_0.txt"
    ref1 = d / "reference_incident_pol_1.txt"

    spectra = {}   # theta -> (wl, summed proxy)
    table = []
    for th in thetas:
        wl0, p0 = proxy(ref0, d / f"theta_{th:03d}_pol_0.txt", args.quantity)
        wl1, p1 = proxy(ref1, d / f"theta_{th:03d}_pol_1.txt", args.quantity)
        ysum = p0 + p1  # orientation-independent (Ex + Ey)
        spectra[th] = (wl0, ysum)
        sc, sv = band_center(wl0, ysum, *args.short_window)
        mc, mv = band_center(wl0, ysum, *args.main_window)
        table.append((th, sc, sv, mc, mv))

    # ---- report ----
    print(f"Input: {d}  quantity: {args.quantity}  (Ex+Ey polarization sum)")
    print("theta_deg, short_center_um, short_val, main_center_um, main_val")
    for th, sc, sv, mc, mv in table:
        print(f"{th:6d}, {sc:.4f}, {sv:.4g}, {mc:.4f}, {mv:.4g}")
    scs = [r[1] for r in table]; mcs = [r[3] for r in table]
    def spread(v):
        v = [x for x in v if np.isfinite(x)]
        return (max(v) - min(v)) if v else float("nan")
    print(f"\nshort-band center spread over theta: {spread(scs)*1000:.1f} nm")
    print(f"main-band  center spread over theta: {spread(mcs)*1000:.1f} nm")
    print("(small spread => resonance wavelengths ~independent of structure angle, per Wu)")

    # ---- figure ----
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.titlesize": 12, "savefig.dpi": 300, "savefig.bbox": "tight"})
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)

    # left: stacked polarization-summed spectra
    off = 0.0
    step = 1.1 * max(np.nanmax(y) for _, y in spectra.values())
    for th in thetas:
        wl, y = spectra[th]
        axL.plot(wl, y + off, color=NAVY, lw=1.6)
        axL.text(wl[-1], off, f"  {th} deg", va="center", fontsize=9, color=NAVY)
        off += step
    axL.set_xlabel("wavelength (um)"); axL.set_ylabel(f"{args.quantity} proxy (Ex+Ey, offset)")
    axL.set_title("Polarization-summed spectra by theta")
    axL.set_yticks([])

    # right: band centers vs theta
    axR.plot(thetas, [c * 1000 for c in scs], "o-", color=TEAL, label="short band")
    axR.plot(thetas, [c * 1000 for c in mcs], "s-", color=GOLD, label="main band")
    axR.set_xlabel("structure angle theta (deg)")
    axR.set_ylabel("band center (nm)")
    axR.set_title("Band center vs theta (flat = Wu-consistent)")
    axR.grid(True, color="#E2E8F0")
    axR.legend(frameon=False)
    fig.suptitle("Structure-angle dependence of the coupled-mode wavelengths", fontsize=12)

    out = args.output or d / "band_centers_vs_theta.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
