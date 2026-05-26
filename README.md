# SDSMT Simulation Workspace

This repository contains Meep/FDTD simulation scripts, control files, and
analysis utilities for plasmonic and nanophotonic structures.

## Layout

- `asym/`: asymmetric rectangle and related SiO2 simulation studies.
- `nanodimer/`: nanorod and nanorod-dimer validation, plotting, and Wu-inspired
  angle-sweep workflows.
- `diamsweep/`: diameter sweep scripts, reference/object CSV outputs, and
  reflectivity plots.
- `Pitchtest/`: pitch sweep scripts, CSV outputs, and resonance plots.
- `MEEPTESTS/`: small Meep/Fresnel validation tests.
- `*.ctl`: top-level Meep control files kept from earlier experiments.

## Git Policy

The repository is intended to track source code, control files, documentation,
small tabular outputs, and selected plots. Large/generated simulation outputs
such as HDF5 field/flux files, run logs, and Python caches are ignored by
default.

If a generated result is needed for publication or reproducibility, add a short
note describing how it was produced and consider storing the large raw data
outside Git.

## Existing Notes

The most complete workflow note currently lives at:

- `nanodimer/nanorod_wu_faststart/WU_TO_DESIGN_FLOW.md`
