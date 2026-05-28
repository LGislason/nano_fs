# Wu Validation to Design Flow

This workflow keeps the Wu et al. comparison useful without turning the project
into a full experimental reproduction.

## 1. Baseline Geometry Check

Use `theta` for the nanorod dimer structure angle and `phi` for the in-plane
orientation of the whole dimer relative to the fixed source polarization.

For the Wu-style validation, keep:

```bash
THETA_LIST="80"
```

Do not use the theta sweep as the main paper comparison. The paper's primary
angle-resolved result is closer to a `phi` or polarization-angle sweep.

## 2. Fast Phi Sweep

Run a cheap qualitative pass first:

```bash
NP=4 RES=200 GAP=0.010 NFREQ=121 \
THETA_LIST="80" \
PHI_LIST="0 20 40 60 80 90 100 120 140 160 180" \
RESULTS_DIR="results_wu_phi_res200_gap10" \
./run_wu_fast_background.sh
```

Check the run:

```bash
python3 check_wu_results.py results_wu_phi_res200_gap10/reference_incident.txt \
  results_wu_phi_res200_gap10/theta_080_phi_*.txt
```

Plot the spectra:

```bash
python3 plot_stacked_spectra.py \
  --input-dir results_wu_phi_res200_gap10 \
  --quantity backscatter \
  -o results_wu_phi_res200_gap10/stacked_backscatter_phi_sweep.png

python3 plot_stacked_spectra.py \
  --input-dir results_wu_phi_res200_gap10 \
  --quantity extinction \
  -o results_wu_phi_res200_gap10/stacked_extinction_phi_sweep.png
```

Extract mode intensities and cosine-squared fit quality:

```bash
python3 analyze_phi_sweep.py \
  --input-dir results_wu_phi_res200_gap10 \
  --theta 80 \
  --quantity backscatter \
  --modes 0.672 0.950 \
  --window 0.04
```

For the refined `RES=300`, `GAP=6 nm` run, the useful long-wavelength angular
feature appears closer to `0.950 um` than `0.853 um`. Treat `0.853 um` as a
poor probe for this simplified model unless a later geometry or illumination
change moves the second feature back into that window.

## 3. Pass Criteria

Treat the model as good enough for a Wu-inspired baseline if:

- `phi = 0` and `phi = 180` are similar.
- The two modal intensities vary strongly with `phi`.
- One mode strengthens while the other weakens over part of the sweep.
- The cosine-squared-family fit has a useful trend, roughly `R2 > 0.7`.
- If a nominal literature mode gives a poor fit, inspect the stacked spectra
  and move the probe to the simulated peak before changing flux monitors.

If the fit is poor but the spectra clearly move with angle, the model is still
useful for design exploration, but it is not a strong Wu reproduction.

## 4. One Refinement Step

If the fast pass looks promising, run one refinement:

```bash
NP=8 RES=300 GAP=0.006 NFREQ=161 \
THETA_LIST="80" \
PHI_LIST="0 20 40 60 80 90 100 120 140 160 180" \
RESULTS_DIR="results_wu_phi_res300_gap6" \
./run_wu_fast_background.sh
```

Then repeat the plotting and `analyze_phi_sweep.py` commands on the new results.

## 5. Move to Your Own Simulations

After one refined validation, stop chasing exact Wu agreement unless the paper
match is the actual project goal. For design work, choose one objective:

- maximize near-field intensity in the gap
- maximize polarization contrast between two `phi` angles
- tune a resonance to a target wavelength
- increase directional scattering or backscatter proxy
- compare rods against asymmetric rectangles

Then sweep the design variables that directly control that objective: gap,
rod length, rod width, `theta`, and orientation/polarization angle.
