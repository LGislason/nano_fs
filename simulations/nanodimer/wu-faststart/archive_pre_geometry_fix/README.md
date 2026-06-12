# Archived pre-geometry-fix runs

These nanodimer Wu phi/theta-sweep runs were all produced **before** the
symmetric-geometry bug was fixed (commit `2314da1`, "Fix symmetric nanodimer
geometry", 2026-06-04) and are therefore superseded. They are kept here on disk
for reference and remain recoverable from git history.

The current, post-fix golden run lives one level up at:

    ../results_sym_gap0006nm_res200_nfreq121_theta080_phi_sweep/

## Archived runs (all dated 2026-05-19 .. 2026-05-28, pre-fix)

- `results_smoke_res80_nfreq31_theta80_phi000_090_180/` — early smoke test
- `results052026/` — phi sweep, 05/20
- `results_res400_nfreq121_thetasweep_phi00/` — theta sweep, res 400
- `results_theta_sweep_res200_nfreq121_phi0/` — theta sweep, res 200
- `results_res200_nfreq121_theta80_phi30/` — single phi=30 point
- `results_phi_res300_gap6/` — refined phi sweep, res 300 (see its ANALYSIS.md)

All use the old geometry and should not be used for publication or comparison.

## Archived field snapshots (also pre-fix, dated 2026-06-03)

- `field_phi040_wvl0672/` — gap near-field snapshot, old geometry
- `field_phi120_wvl0950/` — gap near-field snapshot, old geometry (0.950 um is
  not one of the bands the corrected runs actually produce)

The corrected, post-fix field snapshot lives one level up at
`../field_sym_gap0006nm_theta080_phi120_wvl0792/`.
