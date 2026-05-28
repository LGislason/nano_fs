# Phi Sweep Analysis: RES=300, GAP=6 nm

This run uses `theta = 80 deg`, `gap = 0.006 um`, `resolution = 300`, and
`NFREQ = 161`.

## Result

The short-wavelength mode is centered around `0.67-0.70 um` and shows a strong
cosine-squared-family angular trend.

| quantity | probe | min | max | R2 |
| --- | ---: | ---: | ---: | ---: |
| backscatter | 0.672 um | 0.000464759 | 0.00909607 | 0.973 |
| extinction | 0.672 um | 0.0004655 | 0.00908595 | 0.974 |

The originally probed `0.853 um` wavelength is not a good second-mode probe for
this refined run.

| quantity | probe | min | max | R2 |
| --- | ---: | ---: | ---: | ---: |
| backscatter | 0.853 um | 0.00027492 | 0.000507679 | 0.053 |
| extinction | 0.853 um | 0.000276806 | 0.000511561 | 0.044 |

The meaningful long-wavelength feature is closer to `0.950 um`. It is weaker
than the short-wavelength mode, but it has a much cleaner angular trend.

| quantity | probe | min | max | R2 |
| --- | ---: | ---: | ---: | ---: |
| backscatter | 0.950 um | 0.000168846 | 0.00102738 | 0.886 |
| extinction | 0.950 um | 0.00016615 | 0.00103271 | 0.887 |

## Interpretation

Track `0.672 um` and `0.950 um` for this model. Do not treat the weak
`0.853 um` response as evidence that the flux planes are wrong; the existing
planes capture a strong angular response at `0.672 um` and a usable secondary
response near `0.950 um`.

Before changing monitor placement, inspect the stacked spectra and update the
probe wavelengths to the simulated peaks.

## Reproduction

From `simulations/nanodimer/wu-faststart`:

```bash
conda run -n meep-env python analyze_phi_sweep.py \
  --input-dir results_phi_res300_gap6 \
  --theta 80 \
  --quantity backscatter \
  --modes 0.672 0.950 \
  --window 0.04

conda run -n meep-env python analyze_phi_sweep.py \
  --input-dir results_phi_res300_gap6 \
  --theta 80 \
  --quantity extinction \
  --modes 0.672 0.950 \
  --window 0.04
```
