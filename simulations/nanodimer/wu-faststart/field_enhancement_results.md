# Field Enhancement Analysis — results paragraph

Draft text for the "Field Enhancement Analysis" section of the weekly report.
Numbers are the reference-normalized, metal-masked values from
`compute_field_enhancement.py` (reference run `ref_wvl0672`, |E0|^2 = 0.079).

---

To quantify the near-field response, the quarter-period field snapshots were
converted to an absolute intensity-enhancement map, |E|^2 / |E0|^2, following the
procedure recorded in the field methodology note. The incident reference |E0|^2
was taken from a separate no-metal run using the identical source and cell with
the dimer removed; because the background medium is non-dispersive, this single
reference (|E0|^2 = 0.079 in source-normalized units) applies at every
wavelength. Metal pixels are masked using the dielectric slice and isolated
single-pixel staircasing spikes are removed by a neighbour comparison; because
the 6 nm gap is only ~2 pixels wide at this resolution, the single hottest pixels
sit on the gold tip surface and are resolution-limited, so enhancements are
reported as a robust percentile of the field and as a typical (median) gap value
rather than as the raw maximum.

At the main near-infrared band (lambda ~ 0.79 um, phi = 120 deg — the bonding
orientation), the dimer junction shows an intensity enhancement of order 10^3: a
robust 99.9th-percentile of about 840x across the slice and a typical (median)
gap enhancement near 170x. The single hottest pixels reach ~10^4 but lie on the
gold tip surface and are resolution-limited, so they are not quoted as the
enhancement. The short-wavelength band (lambda ~ 0.67 um, phi = 40 deg) is
weaker, with a robust peak near 330x and a median gap value near 50x. Because
both maps are normalized to the same incident reference, the bonding mode being
roughly three times stronger in the junction is a direct, like-for-like
comparison and the expected signature of the bonding mode, in which charge
accumulates across the gap. Away from the dimer the enhancement relaxes to unity,
which confirms the reference normalization.

Because second-harmonic generation scales as the square of the local intensity, a
junction intensity enhancement of order 10^3 implies a local SHG enhancement of
order 10^6, which is the quantitative motivation for placing a 2D emitter such as
hBN in the gap. These values should be read as lower bounds: at the present
resolution (~3.3 nm voxels) the 6 nm gap is spanned by only about two cells, so
the field is spatially under-sampled, and the reference geometry's ~1 nm gap
would couple more strongly still.

---

**Figure caption (for the enhancement figure).** Near-field intensity
enhancement |E|^2 / |E0|^2 at the rod midplane for the bonding mode
(theta = 80 deg, gap = 6 nm, phi = 120 deg, lambda = 0.79 um), normalized to a
no-metal incident reference and plotted on a logarithmic scale with the gold
masked (grey). Left: full slice, with the dashed box marking the zoom region.
Right: gap zoom showing the junction hotspot. The far-field level relaxes to
~1, confirming the normalization.

---

## Quick-reference numbers (not for the report body)

| metric (|E|^2/|E0|^2) | phi = 40 deg (short, 0.67 um) | phi = 120 deg (bonding, 0.79 um) |
| --- | --- | --- |
| robust slice peak (99.9 pct) | 331x | 838x |
| typical gap value (median) | 51x | 173x |
| on-surface max (resolution-limited; do NOT quote) | ~0.9x10^3 | ~1.9x10^4 |

Reference: |E0|^2 = 0.0793 (no-dimer run `ref_wvl0672`, same for both cases).
Values from `--grow 0` (gap kept visible) with default despike. The on-surface
maximum sits on the under-resolved gold tip and is not grid-robust; headline the
99.9-percentile and median instead. The bonding/short ratio is ~2.5x on the
robust peak and ~3.4x on the median.
