; ============================================================
; Gold nanorod / dimer validation simulation
; Updated legacy-style Meep control file
;
; mode = 0 incident run: substrate only
; mode = 1 structure run: substrate + Au rod(s)
;
; pol = 0 Ex
; pol = 1 Ey
;
; geom = 0 single rod
; geom = 1 rod dimer
;
; Example:
; meep-openmpi mode=0 pol=0 geom=0 nanorod_validation.ctl > inc_ex_single.txt
; meep-openmpi mode=1 pol=0 geom=0 nanorod_validation.ctl > rod_ex_single.txt
; ============================================================

; ----------------------------
; Defaults
; ----------------------------
(if (not (defined? 'mode)) (define mode 0))
(if (not (defined? 'pol))  (define pol 0))
(if (not (defined? 'geom)) (define geom 0))

; ----------------------------
; Gold model: Drude + 3 Lorentz terms
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
; Materials
; ----------------------------
(define substrate-epsilon 2.1)
(define SiO2 (make dielectric (epsilon substrate-epsilon)))

; ----------------------------
; Utility helpers
; ----------------------------
(define (deg->rad deg)
  (* deg (/ pi 180.0)))

; Capsule-like nanorod: cylinder + hemispherical endcaps.
(define (make-capsule-rod cx cy cz ux uy uz)
  (list
   (make cylinder
         (center (vector3 cx cy cz))
         (axis (vector3 ux uy uz))
         (height rod-body-length)
         (radius rod-radius)
         (material gold-3term))
   (make sphere
         (radius rod-radius)
         (center (vector3 (+ cx (* 0.5 rod-body-length ux))
                          (+ cy (* 0.5 rod-body-length uy))
                          (+ cz (* 0.5 rod-body-length uz))))
         (material gold-3term))
   (make sphere
         (radius rod-radius)
         (center (vector3 (- cx (* 0.5 rod-body-length ux))
                          (- cy (* 0.5 rod-body-length uy))
                          (- cz (* 0.5 rod-body-length uz))))
         (material gold-3term))))

; ----------------------------
; Cell
; ----------------------------
; Lateral periodicity is disabled so this behaves like an isolated
; particle above a finite-width piece of substrate rather than an array.
(define sx 0.80)
(define sy 0.80)
(define sz 2.40)
(define dpml 0.25)
(define resolution 250)

(set! geometry-lattice
      (make lattice (size (vector3 sx sy sz))))

(set! pml-layers
      (list (make pml (thickness dpml) (direction X))
            (make pml (thickness dpml) (direction Y))
            (make pml (thickness dpml) (direction Z))))

(set! ensure-periodicity false)

; ----------------------------
; Nanorod dimensions
; ----------------------------
(define rod-length 0.069) ; 69 nm end-to-end
(define rod-width  0.024) ; 24 nm diameter
(define rod-height 0.024) ; 24 nm diameter

(define rod-radius (/ rod-width 2))
(define rod-body-length (- rod-length (* 2 rod-radius)))

; Dimer geometry
; The second rod is rotated 80 degrees in the x-y plane.
(define gap 0.010)
(define rod-angle-deg 80.0)
(define rod-angle-rad (deg->rad rod-angle-deg))

; substrate top at z = 0
; rod sits on top of substrate
(define rod-z rod-radius)

; axis vectors
(define rod1-ux 1.0)
(define rod1-uy 0.0)
(define rod2-ux (cos rod-angle-rad))
(define rod2-uy (sin rod-angle-rad))

; Place the closest tips gap apart, then shift so the dimer is centered.
(define rod1-cx-raw (- 0 (/ gap 2) (/ rod-length 2)))
(define rod1-cy-raw 0.0)

(define rod2-cx-raw (+ (/ gap 2) (* 0.5 rod-length rod2-ux)))
(define rod2-cy-raw (* 0.5 rod-length rod2-uy))

(define dimer-shift-x (* -0.5 (+ rod1-cx-raw rod2-cx-raw)))
(define dimer-shift-y (* -0.5 (+ rod1-cy-raw rod2-cy-raw)))

(define rod1-cx (+ rod1-cx-raw dimer-shift-x))
(define rod1-cy (+ rod1-cy-raw dimer-shift-y))
(define rod2-cx (+ rod2-cx-raw dimer-shift-x))
(define rod2-cy (+ rod2-cy-raw dimer-shift-y))

; ----------------------------
; Substrate block
; ----------------------------
; This still approximates an infinite substrate with a finite lateral
; chunk, but it avoids the stronger error of periodic nanorod copies.
(define substrate-thickness 1.20)
(define substrate
  (make block
        (size (vector3 sx sy substrate-thickness))
        (center (vector3 0 0 (- (/ substrate-thickness 2))))
        (material SiO2)))

; ----------------------------
; Gold rods
; ----------------------------
(define single-rod
  (make-capsule-rod 0 0 rod-z 1 0 0))

(define rod-1
  (make-capsule-rod rod1-cx rod1-cy rod-z rod1-ux rod1-uy 0))

(define rod-2
  (make-capsule-rod rod2-cx rod2-cy rod-z rod2-ux rod2-uy 0))

(define rod-geometry
  (if (= geom 0)
      (append (list substrate) single-rod)
      (append (list substrate) rod-1 rod-2)))

; ----------------------------
; Spectrum
; ----------------------------
(define wvl-min 0.45)
(define wvl-max 1.20)
(define nfreq 250)

(define fmin (/ 1 wvl-max))
(define fmax (/ 1 wvl-min))
(define fcen (* 0.5 (+ fmin fmax)))
(define df (- fmax fmin))

; ----------------------------
; Source and monitors
; ----------------------------
(define src-comp
  (if (= pol 0) Ex Ey))

(define halfz (/ sz 2))
(define monitor-sx (- sx (* 2 dpml) 0.02))
(define monitor-sy (- sy (* 2 dpml) 0.02))

(define src-z  (- halfz dpml 0.12))
(define refl-z (- halfz dpml 0.24))
(define tran-z (+ (- halfz) dpml 0.18))

(define flux-tag-refl
  (if (= pol 0)
      (if (= geom 0) "inc_refl_ex_single" "inc_refl_ex_dimer")
      (if (= geom 0) "inc_refl_ey_single" "inc_refl_ey_dimer")))

(define flux-tag-tran
  (if (= pol 0)
      (if (= geom 0) "inc_tran_ex_single" "inc_tran_ex_dimer")
      (if (= geom 0) "inc_tran_ey_single" "inc_tran_ey_dimer")))

(set! sources
      (list
       (make source
             (src (make gaussian-src (frequency fcen) (fwidth df)))
             (component src-comp)
             (center (vector3 0 0 src-z))
             (size (vector3 monitor-sx monitor-sy 0)))))

(define refl-region
  (make flux-region
        (center (vector3 0 0 refl-z))
        (size (vector3 monitor-sx monitor-sy 0))))

(define tran-region
  (make flux-region
        (center (vector3 0 0 tran-z))
        (size (vector3 monitor-sx monitor-sy 0))))

(define refl 0)
(define tran 0)

; ----------------------------
; Run
; ----------------------------
; Use field decay rather than a fixed stop time to reduce truncation error.
(define decay-point (vector3 0 0 (+ rod-z rod-radius 0.02)))

(if (= mode 0)
    (begin
      (set! geometry (list substrate))
      (set! refl (add-flux fcen df nfreq refl-region))
      (set! tran (add-flux fcen df nfreq tran-region))
      (run-sources+
       (stop-when-fields-decayed 50 src-comp decay-point 1e-6))
      (save-flux flux-tag-refl refl)
      (save-flux flux-tag-tran tran)
      (display-fluxes refl tran))
    (begin
      (set! geometry rod-geometry)
      (set! refl (add-flux fcen df nfreq refl-region))
      (set! tran (add-flux fcen df nfreq tran-region))
      (load-minus-flux flux-tag-refl refl)
      (run-sources+
       (stop-when-fields-decayed 50 src-comp decay-point 1e-6))
      (display-fluxes refl tran)))
