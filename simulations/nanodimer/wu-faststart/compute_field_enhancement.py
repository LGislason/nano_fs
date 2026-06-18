#!/usr/bin/env python3
"""Absolute field-intensity enhancement |E|^2 / |E0|^2 from nanorod field snapshots.

Implements the procedure in field_methodology_notes.txt, with the fixes that
note calls for:

  * Caveat 3 (staircasing): metal and the 1-pixel ring around it are masked out
    using the eps slice, so spurious single-pixel spikes on the metal boundary
    cannot inflate the reported peak. Headline numbers use percentiles + a gap
    region, never the raw maximum.
  * Caveats 1 & 2 (source-dependent / contaminated reference): |E0|^2 can be
    taken from a clean no-dimer reference run (--reference-dir, produced with
    dimer=0 in the CTL) instead of the panel border.
  * Reproducibility note: every case is exported with the identical procedure,
    and --shared-vmax lets two modes (e.g. phi=40 vs phi=120) share one absolute
    color scale for an honest comparison.

Requires the quarter-period "_q" companion snapshots (current CTL). Without them
the true amplitude cannot be formed and the script refuses to guess.

Usage
-----
    # single case, border reference
    python compute_field_enhancement.py field_sym_..._phi120_wvl0792

    # clean reference from a no-dimer run (dimer=0):
    python compute_field_enhancement.py field_sym_..._phi120_wvl0792 \
        --reference-dir ref_nodimer_wvl0792

    # two modes on a shared absolute scale:
    python compute_field_enhancement.py CASE_A --save-vmax /tmp/vmax.txt
    python compute_field_enhancement.py CASE_B --shared-vmax $(cat /tmp/vmax.txt)
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import re

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

FIELD_SX = 0.52  # um, lateral half-extent * 2 of the xy output volumes
ZOOM_NM = 65.0   # half-width of the gap-zoom panel, nm

# Publication-style defaults (kept local so the script stays self-contained).
PLOT_STYLE = {
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "axes.linewidth": 0.8,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
}


def parse_case(name: str) -> str:
    """Build a human-readable parameter subtitle from a results-dir name."""
    bits = []
    for label, pat, fmt in (
        ("theta", r"theta(\d+)", lambda v: f"θ = {int(v)}°"),
        ("phi", r"phi(\d+)", lambda v: f"φ = {int(v)}°"),
        ("gap", r"gap(\d+)nm", lambda v: f"gap = {int(v)} nm"),
        ("wvl", r"wvl(\d+)", lambda v: f"λ = {int(v)} nm"),
    ):
        m = re.search(pat, name)
        if m:
            bits.append(fmt(m.group(1)))
    return ",  ".join(bits)


# --------------------------------------------------------------------------- #
# HDF5 reading (same conventions as plot_field_snapshot_h5.py)
# --------------------------------------------------------------------------- #
def _first_dataset(group: h5py.Group) -> np.ndarray:
    found: list[np.ndarray] = []
    group.visititems(lambda _n, o: found.append(np.asarray(o)) if isinstance(o, h5py.Dataset) else None)
    if not found:
        raise ValueError("No datasets found")
    found.sort(key=lambda d: d.ndim, reverse=True)
    return found[0]


def _read(path: Path) -> np.ndarray:
    with h5py.File(path, "r") as h:
        return np.squeeze(_first_dataset(h))


def _find(input_dir: Path, stem: str) -> Path | None:
    m = sorted(input_dir.glob(f"*{stem}*.h5"))
    return m[0] if m else None


def amplitude_intensity(input_dir: Path, prefix: str) -> tuple[np.ndarray, np.ndarray]:
    """Return (|E|^2, eps) for one slice, using the quarter-period _q pairs."""
    eps_p = _find(input_dir, f"{prefix}_eps")
    eps = _read(eps_p) if eps_p else None
    e2 = np.zeros(())
    for comp in ("ex", "ey", "ez"):
        main = _find(input_dir, f"{prefix}_{comp}")
        quad = _find(input_dir, f"{prefix}_{comp}_q")
        if main is None or quad is None:
            raise FileNotFoundError(
                f"Missing {prefix}_{comp}(_q) in {input_dir}. This script needs the "
                "quarter-period '_q' snapshots from the current CTL."
            )
        a = np.abs(_read(main)) ** 2
        b = np.abs(_read(quad)) ** 2
        e2 = a + b if e2.ndim == 0 else e2 + a + b
    return e2, eps


# --------------------------------------------------------------------------- #
# Masking and reference
# --------------------------------------------------------------------------- #
def dilate(mask: np.ndarray) -> np.ndarray:
    out = mask.copy()
    out[1:, :] |= mask[:-1, :]
    out[:-1, :] |= mask[1:, :]
    out[:, 1:] |= mask[:, :-1]
    out[:, :-1] |= mask[:, 1:]
    return out


def metal_mask(eps: np.ndarray | None, threshold: float, grow: int) -> np.ndarray:
    """True where a pixel is metal or within `grow` pixels of the metal surface."""
    if eps is None:
        return np.zeros((1, 1), dtype=bool)
    m = eps > threshold
    for _ in range(max(grow, 0)):
        m = dilate(m)
    return m


def incident_reference(
    field: np.ndarray, border: int, reference_dir: Path | None, prefix: str
) -> tuple[float, str]:
    if reference_dir is not None:
        ref2, _ = amplitude_intensity(reference_dir, prefix)
        return float(np.median(ref2)), f"no-dimer run ({reference_dir.name})"
    edges = np.concatenate(
        [field[:border, :].ravel(), field[-border:, :].ravel(),
         field[:, :border].ravel(), field[:, -border:].ravel()]
    )
    return float(np.median(edges)), f"panel border ({border}px median)"


# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_dir", type=Path)
    ap.add_argument("--slice", default="xy_rod", choices=["xy_rod", "xy_above", "xz"])
    ap.add_argument("--reference-dir", type=Path, default=None,
                    help="A dimer=0 run to use for a clean |E0|^2 (resolves caveats 1 & 2).")
    ap.add_argument("--border", type=int, default=8, help="Border width for the fallback reference.")
    ap.add_argument("--metal-threshold", type=float, default=2.5,
                    help="eps above this is treated as metal (bg=nbg^2~1.64, Au eps_inf~4.89).")
    ap.add_argument("--grow", type=int, default=1, help="Pixels to grow the metal mask (kills boundary spikes).")
    ap.add_argument("--gap-halfwidth", type=float, default=0.02, help="Half-width (um) of the gap stats box.")
    ap.add_argument("--shared-vmax", type=float, default=None, help="Fix the color-scale max for cross-case comparison.")
    ap.add_argument("--save-vmax", type=Path, default=None, help="Write the chosen vmax here (feed to --shared-vmax).")
    ap.add_argument("-o", "--output", type=Path, default=None)
    args = ap.parse_args()

    field, eps = amplitude_intensity(args.input_dir, args.slice)
    e0_sq, ref_desc = incident_reference(field, args.border, args.reference_dir, args.slice)
    if e0_sq <= 0:
        raise SystemExit("Non-positive |E0|^2; cannot normalize.")
    enh = field / e0_sq

    masked = metal_mask(eps, args.metal_threshold, args.grow)
    if masked.shape != enh.shape:
        masked = np.zeros_like(enh, dtype=bool)
    valid = enh[~masked]

    # gap-region box (dielectric pixels only)
    ny, nx = enh.shape
    xs = np.linspace(-FIELD_SX / 2, FIELD_SX / 2, nx)
    ys = np.linspace(-FIELD_SX / 2, FIELD_SX / 2, ny)
    gx = np.abs(xs) <= args.gap_halfwidth
    gy = np.abs(ys) <= args.gap_halfwidth
    gapbox = np.outer(gy, gx) & ~masked
    gapvals = enh[gapbox]

    p999 = float(np.percentile(valid, 99.9)) if valid.size else float("nan")
    summary = {
        "case": args.input_dir.name,
        "slice": args.slice,
        "reference": ref_desc,
        "E0_squared": e0_sq,
        "raw_max_all": float(np.max(enh)),
        "raw_max_excluding_metal": float(np.max(valid)) if valid.size else float("nan"),
        "p99.9_excluding_metal": p999,
        "p99.5_excluding_metal": float(np.percentile(valid, 99.5)) if valid.size else float("nan"),
        "gap_box_p99.9": float(np.percentile(gapvals, 99.9)) if gapvals.size else float("nan"),
        "gap_box_median": float(np.median(gapvals)) if gapvals.size else float("nan"),
        "n_pixels_over_1e4_on_metal_ring": int(np.sum((enh > 1e4) & masked)),
        "metal_pixels_masked": int(np.sum(masked)),
    }

    print("=== field enhancement |E|^2/|E0|^2 ===")
    for k, v in summary.items():
        print(f"  {k:32s}: {v:.4g}" if isinstance(v, float) else f"  {k:32s}: {v}")
    print(f"  headline (robust): ~{p999:.2g}x intensity enhancement (99.9 pct, metal masked)")

    vmax = args.shared_vmax if args.shared_vmax else p999
    vmin = max(vmax / 1e4, 1e-2)
    if args.save_vmax:
        args.save_vmax.write_text(f"{vmax:.6g}\n")

    # ---- publication-style figure: full slice | gap zoom, side by side ----
    plt.rcParams.update(PLOT_STYLE)
    out = args.output or args.input_dir / f"enhancement_{args.slice}.png"
    half = FIELD_SX / 2 * 1000.0  # nm
    extent = (-half, half, -half, half)
    enh_img = np.where(masked, np.nan, enh).T  # hide metal, transpose to image orientation
    eps_img = eps.T if eps is not None else None

    cmap = plt.cm.inferno.copy()
    cmap.set_bad("0.55")  # neutral grey for masked metal
    norm = mcolors.LogNorm(vmin=vmin, vmax=vmax)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.9), constrained_layout=True)
    im = None
    for ax, zoom, title in ((axes[0], None, "Full slice"), (axes[1], ZOOM_NM, "Gap zoom")):
        im = ax.imshow(enh_img, origin="lower", extent=extent, cmap=cmap, norm=norm,
                       interpolation="nearest", rasterized=True)
        if eps_img is not None and eps_img.max() > args.metal_threshold:  # outline the metal
            xs = np.linspace(-half, half, eps_img.shape[1])
            ys = np.linspace(-half, half, eps_img.shape[0])
            ax.contour(xs, ys, eps_img, levels=[args.metal_threshold],
                       colors="white", linewidths=0.9, alpha=0.9)
        ax.set_aspect("equal")
        ax.set_xlabel("x (nm)")
        ax.set_ylabel("y (nm)")
        ax.set_title(title)
        if zoom:
            ax.set_xlim(-zoom, zoom)
            ax.set_ylim(-zoom, zoom)
        else:
            ax.add_patch(Rectangle((-ZOOM_NM, -ZOOM_NM), 2 * ZOOM_NM, 2 * ZOOM_NM,
                                   fill=False, edgecolor="white", linewidth=0.8, linestyle="--"))
            ax.text(0.04, 0.96, f"peak ≈ {p999:.0f}×", transform=ax.transAxes,
                    va="top", ha="left", color="white", fontsize=10,
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="black", alpha=0.55, edgecolor="none"))

    cbar = fig.colorbar(im, ax=axes, shrink=0.92, pad=0.02,
                        label=r"intensity enhancement  $|E|^2 / |E_0|^2$")
    cbar.ax.tick_params(labelsize=9)
    subtitle = parse_case(args.input_dir.name)
    fig.suptitle(f"Near-field intensity enhancement\n{subtitle}" if subtitle
                 else "Near-field intensity enhancement", fontsize=12)
    fig.savefig(out)
    plt.close(fig)

    np.savetxt(args.input_dir / f"enhancement_{args.slice}.txt", enh,
               header=f"|E|^2/|E0|^2  {args.slice}  E0^2={e0_sq:.6g}  ref={ref_desc}")
    (args.input_dir / f"enhancement_{args.slice}_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"Wrote {out}")
    print(f"Wrote {args.input_dir / f'enhancement_{args.slice}.txt'} and *_summary.json")


if __name__ == "__main__":
    main()
