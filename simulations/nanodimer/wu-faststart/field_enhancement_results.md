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
wavelength. Metal pixels and the one-pixel ring around them were masked using the
dielectric slice so that grid-staircasing spikes on the gold boundary cannot
inflate the reported values, and peak enhancements are quoted as the 99.9th
percentile of the masked field and within a +/- 20 nm gap box rather than as the
raw maximum.

At the main near-infrared band (lambda ~ 0.79 um, phi = 120 deg — the bonding
orientation), the dimer junction shows an intensity enhancement of order 10^3,
with a robust gap-region peak near 4.5 x 10^3 and a median gap enhancement near
90x. The short-wavelength band (lambda ~ 0.67 um, phi = 40 deg) is markedly
weaker in the junction, with a gap-region peak near 2.9 x 10^2 and a median near
50x. Because both maps are normalized to the same incident reference, the
roughly fifteen-fold stronger junction field of the phi = 120 deg mode is a
direct, like-for-like comparison and is the expected signature of the bonding
mode, in which charge accumulates across the gap. Away from the dimer the
enhancement relaxes to unity, which confirms the reference normalization.

Because second-harmonic generation scales as the square of the local intensity,
a junction intensity enhancement of order 4.5 x 10^3 implies a local SHG
enhancement of order 10^7, which is the quantitative motivation for placing a 2D
emitter such as hBN in the gap. These values should be read as lower bounds: at
the present resolution (~3.3 nm voxels) the 6 nm gap is spanned by only about two
cells, so the true peak is spatially under-sampled, and the reference geometry's
~1 nm gap would couple more strongly still.

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
| robust slice peak (99.9 pct) | 181x | 446x |
| gap-box peak (99.9 pct) | 291x | 4540x |
| gap-box median | 49x | 91x |
| spurious >1e4 metal-ring px (masked) | 0 | 2 (raw max 19393x before masking) |

Reference: |E0|^2 = 0.0793 (no-dimer run `ref_wvl0672`, same for both cases).
