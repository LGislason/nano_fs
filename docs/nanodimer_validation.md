# Nanorod-Dimer FDTD Validation

**Goal.** Validate the FDTD simulation workflow against the gold nanorod dimer of
Wu et al. (Nano-Micro Lett. **6**, 372–380, 2014) before applying it to the
asymmetric-rectangle metasurface. This is a **qualitative** validation at a 6 nm
gap under normal incidence — strong enough to trust the workflow, not a 1:1
reproduction of Wu's 1 nm dark-field experiment.

**Bottom line.** Three independent checks pass: (1) the orientation (φ) sweep
reproduces the cos²φ two-mode polarization dependence Wu reports; (2) the
near-field maps give a physically sensible bonding-mode junction enhancement of
order 10³ on a clean reference; (3) the structure-angle (θ) sweep shows the
expected effective-length redshift. The φ result is converged (res 200 → 300).

---

## 1. Methods

- **Solver:** Meep 1.1.1 (Scheme/libctl), MPI on a 32-core workstation.
- **Excitation:** broadband Gaussian plane wave at normal incidence; separate
  Ex / Ey polarizations (`pol` switch). Reflection and transmission flux monitors
  are placed away from the near field; runs stop on `stop-when-fields-decayed`
  rather than a fixed clock, so spectra are not truncation-limited.
- **Material:** gold as Drude + 3-Lorentz, independently validated against Mie
  theory for a 30 nm sphere (resonance position and long-wavelength decay agree).
- **Geometry:** two capsule-shaped Au rods, 69 × 24 nm, in a homogeneous
  background (n = 1.28, no substrate) for cheap sweeps. θ is the internal opening
  angle between the rods; φ is the in-plane rotation of the whole dimer relative
  to the fixed Ex source.
- **Corrected gap:** the symmetric construction sets the input gap equal to the
  true closest metal-to-metal surface separation. Verified analytically — for
  θ = 80°, gap = 6 nm the rod-tip reference spacing is 14.573 nm and the true
  surface gap is 6.000 nm, holding across θ = 30–90°.
- **Resolution:** 200 px/µm for spectra, 300 px/µm for field snapshots (a 6 nm
  gap is ~1.8 px at res 300).
- **Proxies:** "backscatter" = reflected flux / |incident|, "extinction" =
  −transmitted / |incident|. These are directional flux proxies, not dark-field
  scattering cross sections.

---

## 2. Result 1 — Orientation (φ) sweep reproduces cos²φ

Holding the dimer fixed (θ = 80°, gap = 6 nm) and rotating it through φ against a
fixed Ex field, two bands appear — a short band near 0.69 µm and a main NIR band
near 0.79 µm — whose strengths trade off with orientation. Extracting each band
amplitude (two-Lorentzian fit) and fitting the offset cos²φ form
`y = c + a·cos2φ + b·sin2φ` gives anti-correlated curves: the short band peaks
near φ = 20°/60° and the main band near φ = 120°/140°. This is the signature of
two coupled dimer modes driven by orthogonal projections of the incident field,
consistent with the cos²φ polarization dependence Wu reports.

**Symmetry.** The symmetric dimer obeys an exact relation, `I(φ) = I(θ − φ)`
(mod 180°), confirmed in the data (φ = 60° ≡ 20°, 140° ≡ 120°, 160° ≡ 100°,
80°/180° ≡ 0°). The independent angle set is therefore {0, 20, 40, 90, 100, 120};
the 0–180° curve is filled by reflecting these to their symmetry partners
(`plot_phi_angle_response.py --mirror`) without redundant runs.

**Figures:** `results_sym_gap0006nm_res300_nfreq121_phi_sweep/phi_angle_response.png`
(angle response), `.../stacked_backscatter_phi_sweep.png` (stacked spectra).

---

## 3. Result 2 — Near-field enhancement (bonding mode)

Field snapshots were converted to absolute intensity-enhancement maps,
|E|² / |E₀|², with three safeguards: the true amplitude from a quarter-period
snapshot pair (phase-independent), a clean no-metal reference run for the incident
|E₀|² (= 0.079, identical for all wavelengths in the non-dispersive background),
and masking of the metal plus removal of isolated boundary staircasing spikes so
artifacts do not set the reported peak.

| metric (|E|²/|E₀|²) | φ = 40° (short, 0.67 µm) | φ = 120° (bonding, 0.79 µm) |
| --- | --- | --- |
| robust peak (99.9 pct) | 331× | **838×** |
| typical gap value (median) | 51× | **173×** |

The bonding mode concentrates an intensity enhancement of **order 10³** in the
junction (~840× robust peak, ~170× median) — about **3× the short band** on the
same reference, the expected bonding signature. Far from the dimer the
enhancement relaxes to 1, confirming the normalization. Because second-harmonic
generation scales as intensity², this implies a local SHG enhancement of order
10⁶, motivating an emitter (e.g. hBN) in the gap.

**Caveat:** at 6 nm the gap is only ~2 pixels wide, so the single hottest pixels
(~10⁴) sit on the under-resolved gold tip surface and are resolution-limited —
they are reported as a floor, not the headline.

**Figure:** `.../enhancement_field_sym_gap0006nm_theta080_phi120_wvl0790.png`.

---

## 4. Result 3 — Structure-angle (θ) sweep

Sweeping the opening angle θ = 30–150° with a two-polarization (Ex + Ey) sum —
which removes orientation and shows both modes at every angle — the dominant
coupled mode **redshifts ~0.69 → 0.82 µm (~130 nm)** and strengthens ~2.5× as the
dimer opens toward collinear. This is an effective-antenna-length effect.

It does **not** reproduce Wu's finding that the resonance wavelengths are
independent of the structure angle — and that is expected: Wu's statement holds
in a tight (~1 nm) gap, small-angle regime where the gap coupling pins the mode,
whereas this sweep is at a 6 nm gap over a wide angle range, where the overall
geometry controls the wavelength. Reaching the 1 nm regime is compute-limited.

**Figure:**
`results_sym_gap0006nm_res200_nfreq121_theta_polsum/band_centers_vs_theta.png`.

---

## 5. Convergence

The φ sweep was repeated at res 300 (independent angle set). Comparing the
two-Lorentzian fits, res 200 → 300:

- cos²φ fit quality improves: main band R² 0.97 → **0.996**, short band
  0.77 → **0.87**.
- per-angle Lorentzian fits improve: minimum R² 0.53 → **0.99** (res 200 had one
  poor fit at φ = 40°; res 300 fits every angle at R² > 0.99).
- band centers shift only ~7–12 nm (~1–2 %), i.e. ~1–2 spectral bins.

So the φ validation is converged: res 200 was already qualitatively correct and
res 300 tightens it. Quote ~1–2 % wavelength uncertainty. (θ-sweep convergence at
res 250/300 is in progress.)

---

## 6. Caveats / scope

1. **Qualitative, not quantitative.** 6 nm gap vs Wu's ~1 nm; normal-incidence
   backscatter proxy vs Wu's oblique dark-field scattering cross section. Absolute
   peak positions (~0.69 / 0.79 µm here vs Wu's 0.672 / 0.853 µm) differ from the
   gap, the absent substrate, and the n = 1.28 background.
2. **Resolution-limited gap.** The 6 nm junction is ~2 px at res 300, so on-tip
   field values are not grid-robust; headline the percentile/median, not the max.
3. **Compute bound.** The 1 nm gap and full θ sweeps need a modern Meep build or
   HPC; on the current 2009-era stack they are not reachable in practical runtime.

---

## 7. Figure manifest

| # | Topic | File |
| - | ----- | ---- |
| 1 | Simulation layout | `simulation_layout_xz.png` |
| 2 | Dimer geometry (θ=80, φ=120) | `nanodimer_topview_theta080_phi120_gap6nm_corrected.png` |
| 3 | φ stacked spectra | `results_sym_gap0006nm_res300_nfreq121_phi_sweep/stacked_backscatter_phi_sweep.png` |
| 4 | φ angle response (cos²φ, mirrored) | `results_sym_gap0006nm_res300_nfreq121_phi_sweep/phi_angle_response.png` |
| 5 | Bonding-mode enhancement map | `enhancement_field_sym_gap0006nm_theta080_phi120_wvl0790.png` |
| 6 | θ band-center vs angle | `results_sym_gap0006nm_res200_nfreq121_theta_polsum/band_centers_vs_theta.png` |

All paths are under `simulations/nanodimer/wu-faststart/`.
