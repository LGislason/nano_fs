#!/usr/bin/env python3
"""Plot per-mode strength vs phi for the Wu-style nanodimer phi sweep.

Reads the table written by analyze_phi_sweep.py (the lorentzian_mode_fit_*.txt
file), plots each mode's amplitude vs phi, and fits the offset cos^2 form
    y = c + a cos(2 phi) + b sin(2 phi).

--mirror exploits the symmetric dimer's exact symmetry  I(phi) = I(theta - phi)
(mod 180): every computed angle has an equivalent partner at (theta - phi),
so the 0-180 deg curve can be filled WITHOUT re-running those redundant angles.
The mirrored points are drawn open / dashed and labelled as symmetry partners.
The cos^2 fit uses only the independent (computed) angles so redundant points do
not inflate R^2.

Usage:
    python plot_phi_angle_response.py lorentzian_mode_fit_backscatter.txt --mirror
    python plot_phi_angle_response.py <results_dir>/lorentzian_mode_fit_backscatter.txt
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

NAVY = "#12243B"
COLORS = ["#1C7293", "#C08028", "#6D2E46", "#2C5F2D"]


def parse_fit_table(path: Path):
    """Return (theta_deg, quantity, phis, {mode_label: values}). Reads the CSV
    block produced by analyze_phi_sweep.py."""
    text = path.read_text().splitlines()
    theta = 80.0
    quantity = "backscatter"
    header = None
    rows = []
    for line in text:
        m = re.search(r"Theta:\s*([\d.]+)", line)
        if m:
            theta = float(m.group(1))
        m = re.search(r"Quantity:\s*(\w+)", line)
        if m:
            quantity = m.group(1)
        if line.startswith("phi,"):
            header = [h.strip() for h in line.split(",")]
            continue
        if header and line.strip() and line[0].isdigit():
            rows.append([x.strip() for x in line.split(",")])
        elif header and not line.strip():
            if rows:
                break
    if not header or not rows:
        raise SystemExit(f"No fit table found in {path}")

    cols = {h: i for i, h in enumerate(header)}
    phis = [float(r[cols["phi"]]) for r in rows]
    mode_cols = [h for h in header if h.endswith("_value") and h.startswith("mode_")]
    modes = {}
    for mc in mode_cols:
        label = mc[len("mode_"):-len("_value")]  # e.g. "0.690"
        modes[label] = [float(r[cols[mc]]) for r in rows]
    return theta, quantity, phis, modes


def cos2_fit(phis, values):
    """Least-squares y = c + a cos2phi + b sin2phi (no scipy)."""
    A = []
    for p in phis:
        ang = math.radians(2 * p)
        A.append([1.0, math.cos(ang), math.sin(ang)])
    A = np.array(A); y = np.array(values)
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ coef
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot else float("nan")
    return coef, r2


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fit_file", type=Path, help="lorentzian_mode_fit_*.txt from analyze_phi_sweep.py")
    ap.add_argument("--mirror", action="store_true",
                    help="fill 0-180 deg using the dimer symmetry I(phi)=I(theta-phi)")
    ap.add_argument("--theta", type=float, default=None, help="override theta (deg) for the mirror")
    ap.add_argument("-o", "--output", type=Path, default=None)
    args = ap.parse_args()

    theta, quantity, phis, modes = parse_fit_table(args.fit_file)
    if args.theta is not None:
        theta = args.theta

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12,
                         "savefig.dpi": 300, "savefig.bbox": "tight"})
    fig, ax = plt.subplots(figsize=(8.2, 5.0), constrained_layout=True)

    fine = np.linspace(0, 180, 361)
    for i, (label, vals) in enumerate(modes.items()):
        c = COLORS[i % len(COLORS)]
        ax.plot(phis, vals, "o", color=c, ms=9, label=f"{float(label)*1000:.0f} nm (computed)")
        if args.mirror:
            mphi = [(theta - p) % 180 for p in phis]
            ax.plot(mphi, vals, "o", color=c, ms=8, mfc="none", mew=1.5,
                    label=f"{float(label)*1000:.0f} nm (symmetry partner)")
        # cos^2 fit on the INDEPENDENT computed angles only
        coef, r2 = cos2_fit(phis, vals)
        ang = np.radians(2 * fine)
        ax.plot(fine, coef[0] + coef[1] * np.cos(ang) + coef[2] * np.sin(ang),
                "-", color=c, lw=1.8, alpha=0.8,
                label=f"{float(label)*1000:.0f} nm cos²φ fit (R²={r2:.3f})")

    if args.mirror:
        ax.axvline(theta / 2, color="0.6", ls=":", lw=1)
        ax.text(theta / 2, ax.get_ylim()[1], f" mirror axis φ={theta/2:.0f}°",
                va="top", ha="left", fontsize=9, color="0.4")

    ax.set_xlabel("in-plane orientation φ (deg)")
    ax.set_ylabel(f"{quantity} proxy, per-mode amplitude")
    ax.set_xlim(-5, 185); ax.set_xticks(range(0, 181, 30))
    ax.grid(True, color="#E2E8F0")
    ax.set_title(f"Angle-resolved mode strength (θ = {theta:.0f}°)"
                 + ("  —  symmetry-mirrored to 0–180°" if args.mirror else ""))
    ax.legend(frameon=False, fontsize=8, ncol=len(modes))

    out = args.output or args.fit_file.with_name("phi_angle_response.png")
    fig.savefig(out)
    plt.close(fig)
    print(f"theta={theta:g}  modes={list(modes)}  mirror={args.mirror}")
    for label, vals in modes.items():
        _, r2 = cos2_fit(phis, vals)
        print(f"  mode {label}: cos2phi R^2 = {r2:.3f} (from {len(phis)} independent angles)")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
