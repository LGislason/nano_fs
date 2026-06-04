# Wu Fast Start

This directory is a low-cost Meep starting point for qualitatively reproducing
the angle-dependent nanorod-dimer trends in:

Jian Wu et al., "Angle-Resolved Plasmonic Properties of Single Gold Nanorod
Dimers", Nano-Micro Letters 6, 372-380 (2014)

It is intentionally not a strict paper-faithful reproduction yet. It trades
accuracy for runtime so you can do initial sweeps on a remote machine without
waiting days per attempt.

## What This Model Preserves

- Capsule-like Au nanorods
- Rod size close to the paper: 69 x 24 nm
- Dimer structure angle `theta` as a parameter
- In-plane orientation angle `phi` as a parameter
- Broadband response around the reported bonding/antibonding resonances

## What This Model Simplifies

- Uses a homogeneous background medium instead of a substrate
- Uses normal-incidence excitation instead of the paper's oblique dark-field setup
- Uses a larger default surface gap than the paper because 1 nm is not practical for a
  cheap 3D FDTD sweep
- Produces fast directional flux proxies, not a rigorous scattering cross section

## Why It Is Much Cheaper

Compared with the old run:

- smaller cell
- lower default resolution
- narrower wavelength window
- only one reference run for all `theta/phi` cases
- no substrate block

The old setup had about 24 million Yee cells. The default fast-start setup is
about 1.7 million, so each run should be substantially cheaper.

## Files

- `nanorod_wu_fast.ctl`: parameterized Meep control file
- `run_wu_fast_background.sh`: launch script for quick sweeps
- `nanorod_wu_field_snapshot.ctl`: single-wavelength field snapshot control file
- `run_field_snapshots_background.sh`: cluster launcher for high-resolution field snapshots

## Default Interpretation

The structure run subtracts the homogeneous-medium reference on both flux
planes:

- `refl` is backward-scattered flux proxy
- `tran` is forward differential flux proxy

For quick qualitative analysis you can treat:

- `-tran / inc_tran` as an extinction-like proxy
- `refl / inc_tran` as a backscatter proxy

but these are still not paper-grade dark-field scattering cross sections.

## Default Sweep

The run script defaults to:

- `theta = 80 deg`
- `phi = 0, 20, 40, 60, 80, 90, 100, 120, 140, 160, 180 deg`
- `gap = 6 nm` true metal-to-metal surface gap
- `resolution = 200 px/um`
- `background index = 1.28`

This is intended as the corrected-geometry orientation sweep at fixed
`theta=80 deg`.

## Recommended Refinement Path

1. Confirm the angle dependence with the default fast sweep.
2. Increase `resolution` to `250` or `300`.
3. Reduce `gap` from `10 nm` toward `6 nm`, then lower if runtime allows.
4. Only after the trends look stable, consider reintroducing a substrate or a
   more faithful illumination model.

## Usage

From this directory on your external machine:

```bash
chmod +x run_wu_fast_background.sh
./run_wu_fast_background.sh
```

Useful environment overrides:

```bash
NP=8 RES=250 GAP=0.006 THETA_LIST="80" PHI_LIST="0 30 60 90 120 150 180" ./run_wu_fast_background.sh
```

The script writes results into a geometry-tagged directory such as
`results_sym_gap0006nm_res200_nfreq121_theta080_phi_sweep`, unless `RESULTS_DIR`
is set explicitly.
