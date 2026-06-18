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
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

FIELD_SX = 0.52  # um, lateral half-extent * 2 of the xy output volumes


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

    out = args.output or args.input_dir / f"enhancement_{args.slice}.png"
    enh_plot = np.where(masked, np.nan, enh).T  # hide metal, transpose to image orientation
    extent = (-FIELD_SX / 2, FIELD_SX / 2, -FIELD_SX / 2, FIELD_SX / 2)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)
    cmap = plt.cm.inferno.copy()
    cmap.set_bad("dimgray")
    norm = mcolors.LogNorm(vmin=vmin, vmax=vmax)
    for ax, zoom in ((axes[0], None), (axes[1], 0.065)):
        im = ax.imshow(enh_plot, origin="lower", extent=extent, cmap=cmap, norm=norm, interpolation="nearest")
        ax.set_aspect("equal")
        ax.set_xlabel("x (um)"); ax.set_ylabel("y (um)")
        if zoom:
            ax.set_xlim(-zoom, zoom); ax.set_ylim(-zoom, zoom); ax.set_title("Gap zoom")
        else:
            ax.set_title(f"{args.slice}  |E|^2/|E0|^2  (peak~{p999:.2g}x)")
    fig.colorbar(im, ax=axes, shrink=0.85, label="|E|^2 / |E0|^2  (metal masked, log scale)")
    fig.suptitle(args.input_dir.name)
    fig.savefig(out, dpi=220)
    plt.close(fig)

    np.savetxt(args.input_dir / f"enhancement_{args.slice}.txt", enh,
               header=f"|E|^2/|E0|^2  {args.slice}  E0^2={e0_sq:.6g}  ref={ref_desc}")
    (args.input_dir / f"enhancement_{args.slice}_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"Wrote {out}")
    print(f"Wrote {args.input_dir / f'enhancement_{args.slice}.txt'} and *_summary.json")


if __name__ == "__main__":
    main()
