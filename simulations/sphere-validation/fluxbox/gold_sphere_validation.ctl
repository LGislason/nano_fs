; ============================================================
; Gold nanosphere Mie validation
; Meep 1.1.1 compatible
;
; mode = 0 incident/reference run: no sphere
; mode = 1 structure run: gold sphere with incident field subtracted
;
; Goal:
;   Compare total scattered power from a gold sphere to Mie scattering.
;   This uses a closed flux box around the sphere, not a reflectance plane.
; ============================================================

(if (not (defined? 'mode)) (define mode 0))

; ----------------------------
; Gold model: Drude + 3 Lorentz
; ----------------------------
(define omega_meep 1.88365e15)

(define au_epsilon_inf 4.8929)

(define au_omega_d (/ 1.2944e+16 omega_meep))
(define au_omega_1 (/ 1.3617e+15 omega_meep))
(define au_omega_2 (/ 4.1636e+15 omega_meep))
(define au_omega_3 (/ 5.0753e+15 omega_meep))

(define au_gamma_d (/ 1.0003e+09 omega_meep))
(define au_gamma_1 (/ (* 4.7356e+14 2) omega_meep))
(define au_gamma_2 (/ (* 4.4931e+14 2) omega_meep))
(define au_gamma_3 (/ (* 5.8469e+14 2) omega_meep))

(define au_sigma_1 4.7282)
(define au_sigma_2 0.72996)
(define au_sigma_3 1.5103)

(define gold-3term
  (make dielectric
        (epsilon au_epsilon_inf)
        (E-polarizations
         (make polarizability
               (omega 1e-20)
               (gamma au_gamma_d)
               (sigma (* (* 1e20 au_omega_d)
                         (* 1e20 au_omega_d))))
         (make polarizability
               (omega au_omega_1)
               (gamma au_gamma_1)
               (sigma au_sigma_1))
         (make polarizability
               (omega au_omega_2)
               (gamma au_gamma_2)
               (sigma au_sigma_2))
         (make polarizability
               (omega au_omega_3)
               (gamma au_gamma_3)
               (sigma au_sigma_3)))))

; ----------------------------
; Cell
; ----------------------------
(define sx 1.4)
(define sy 1.4)
(define sz 1.4)

(define dpml 0.30)
(define resolution 120)

(set! geometry-lattice
      (make lattice (size (vector3 sx sy sz))))

(set! pml-layers
      (list (make pml (thickness dpml))))

; ----------------------------
; Geometry
; ----------------------------
(define sphere-radius 0.030) ; 30 nm

(define sphere-geometry
  (list
   (make sphere
         (radius sphere-radius)
         (center (vector3 0 0 0))
         (material gold-3term))))

(if (= mode 0)
    (set! geometry '())
    (set! geometry sphere-geometry))

; ----------------------------
; Spectrum
; ----------------------------
(define wvl-min 0.40)
(define wvl-max 0.80)
(define nfreq 200)

(define fmin (/ 1 wvl-max))
(define fmax (/ 1 wvl-min))
(define fcen (* 0.5 (+ fmin fmax)))
(define df (- fmax fmin))

; ----------------------------
; Source
; Propagation along +z
; Ex-polarized plane wave
; ----------------------------
(define halfz (/ sz 2))
(define src-z (- (+ (- halfz) dpml) 0.05))

(set! sources
      (list
       (make source
             (src (make gaussian-src
                        (frequency fcen)
                        (fwidth df)))
             (component Ex)
             (center (vector3 0 0 src-z))
             (size (vector3 (- sx (* 2 dpml))
                            (- sy (* 2 dpml))
                            0)))))

; ----------------------------
; Scattering flux box
; Box must fully surround sphere, but stay away from PML
; ----------------------------
(define box-half 0.25)
(define box-size (* 2 box-half))

(define xplus-region
  (make flux-region
        (center (vector3 box-half 0 0))
        (size (vector3 0 box-size box-size))
        (direction X)))

(define xminus-region
  (make flux-region
        (center (vector3 (- box-half) 0 0))
        (size (vector3 0 box-size box-size))
        (direction X)))

(define yplus-region
  (make flux-region
        (center (vector3 0 box-half 0))
        (size (vector3 box-size 0 box-size))
        (direction Y)))

(define yminus-region
  (make flux-region
        (center (vector3 0 (- box-half) 0))
        (size (vector3 box-size 0 box-size))
        (direction Y)))

(define zplus-region
  (make flux-region
        (center (vector3 0 0 box-half))
        (size (vector3 box-size box-size 0))
        (direction Z)))

(define zminus-region
  (make flux-region
        (center (vector3 0 0 (- box-half)))
        (size (vector3 box-size box-size 0))
        (direction Z)))

(define xplus (add-flux fcen df nfreq xplus-region))
(define xminus (add-flux fcen df nfreq xminus-region))
(define yplus (add-flux fcen df nfreq yplus-region))
(define yminus (add-flux fcen df nfreq yminus-region))
(define zplus (add-flux fcen df nfreq zplus-region))
(define zminus (add-flux fcen df nfreq zminus-region))

; ----------------------------
; Run
; ----------------------------
(if (= mode 0)
    (begin
      (run-sources+ 800)

      (save-flux "sphere_inc_xplus" xplus)
      (save-flux "sphere_inc_xminus" xminus)
      (save-flux "sphere_inc_yplus" yplus)
      (save-flux "sphere_inc_yminus" yminus)
      (save-flux "sphere_inc_zplus" zplus)
      (save-flux "sphere_inc_zminus" zminus)

      (display-fluxes xplus xminus yplus yminus zplus zminus))

    (begin
      (load-minus-flux "sphere_inc_xplus" xplus)
      (load-minus-flux "sphere_inc_xminus" xminus)
      (load-minus-flux "sphere_inc_yplus" yplus)
      (load-minus-flux "sphere_inc_yminus" yminus)
      (load-minus-flux "sphere_inc_zplus" zplus)
      (load-minus-flux "sphere_inc_zminus" zminus)

      (run-sources+ 800)

      (display-fluxes xplus xminus yplus yminus zplus zminus)))
