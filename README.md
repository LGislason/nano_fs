# SDSMT Simulation Workspace

This repository contains Meep/FDTD simulation scripts, control files, and
analysis utilities for plasmonic and nanophotonic structures.

## Layout

- `simulations/asymmetric-rectangles/`: asymmetric rectangle and related SiO2
  simulation studies.
- `simulations/nanodimer/`: nanorod and nanorod-dimer validation, plotting, and
  Wu-inspired angle-sweep workflows.
- `simulations/diameter-sweep/`: diameter sweep scripts, reference/object CSV
  outputs, and reflectivity plots.
- `simulations/pitch-sweep/`: pitch sweep scripts, CSV outputs, and resonance
  plots.
- `simulations/fresnel-tests/`: small Meep/Fresnel validation tests.
- `simulations/sphere-validation/`: gold sphere validation and Mie comparison
  scripts/results.
- `examples/meep-control-files/`: standalone Meep control files kept from
  earlier experiments.
- `docs/`: workflow notes and project status documents.

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

- `docs/wu-validation-to-design-flow.md`
