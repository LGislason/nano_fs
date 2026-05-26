; ============================================================
; Low-cost fast-start proxy for Wu et al. (2014)
; "Angle-Resolved Plasmonic Properties of Single Gold Nanorod Dimers"
;
; This is designed for rapid qualitative sweeps, not paper-grade
; dark-field scattering reproduction.
;
; mode = 0  homogeneous-medium reference
; mode = 1  dimer structure
;
; theta = angle between rod long axes in degrees
; phi   = in-plane rotation of the whole dimer in degrees
;
; A fixed Ex plane wave is used. Sweeping phi changes the angle between
; the excitation polarization and the dimer bisector, which is a cheap
; proxy for the paper's angle-resolved measurements.
; ============================================================

(if (not (defined? 'mode))      (define mode 0))
(if (not (defined? 'theta))     (define theta 80.0))
(if (not (defined? 'phi))       (define phi 0.0))
(if (not (defined? 'res))       (define res 200))
(if (not (defined? 'gap))       (define gap 0.010))
(if (not (defined? 'nbg))       (define nbg 1.28))
(if (not (defined? 'nfreq))     (define nfreq 121))
(if (not (defined? 'wvl_min))   (define wvl_min 0.55))
(if (not (defined? 'wvl_max))   (define wvl_max 0.95))
(if (not (defined? 'decay_time)) (define decay_time 30))
(if (not (defined? 'decay_tol))  (define decay_tol 1e-4))

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
; Helpers
; ----------------------------
(define (deg->rad deg)
  (* deg (/ pi 180.0)))

(define (rot-x x y ang)
  (- (* x (cos ang)) (* y (sin ang))))

(define (rot-y x y ang)
  (+ (* x (sin ang)) (* y (cos ang))))

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
(define sx 0.80)
(define sy 0.80)
(define sz 1.20)
(define dpml 0.12)

(set! geometry-lattice
      (make lattice (size (vector3 sx sy sz))))

(set! resolution res)

(set! pml-layers
      (list (make pml (thickness dpml) (direction X))
            (make pml (thickness dpml) (direction Y))
            (make pml (thickness dpml) (direction Z))))

(set! ensure-periodicity false)
(set! default-material (make dielectric (epsilon (* nbg nbg))))

; ----------------------------
; Rod geometry
; ----------------------------
(define rod-length 0.069)
(define rod-width 0.024)
(define rod-radius (/ rod-width 2))
(define rod-body-length (- rod-length (* 2 rod-radius)))
(define rod-z 0.0)

(define theta-rad (deg->rad theta))
(define phi-rad (deg->rad phi))

; Base dimer: theta is the internal opening angle between rods, measured
; between the two rod directions pointing away from the gap.
(define rod1-out-ang pi)
(define rod2-out-ang (- pi theta-rad))

(define rod1-tip-x (- (/ gap 2)))
(define rod1-tip-y 0.0)
(define rod2-tip-x (/ gap 2))
(define rod2-tip-y 0.0)

(define rod1-bx (+ rod1-tip-x (* 0.5 rod-length (cos rod1-out-ang))))
(define rod1-by (+ rod1-tip-y (* 0.5 rod-length (sin rod1-out-ang))))
(define rod2-bx (+ rod2-tip-x (* 0.5 rod-length (cos rod2-out-ang))))
(define rod2-by (+ rod2-tip-y (* 0.5 rod-length (sin rod2-out-ang))))

(define dimer-shift-x (* -0.5 (+ rod1-bx rod2-bx)))
(define dimer-shift-y (* -0.5 (+ rod1-by rod2-by)))

(define rod1-cx0 (+ rod1-bx dimer-shift-x))
(define rod1-cy0 (+ rod1-by dimer-shift-y))
(define rod2-cx0 (+ rod2-bx dimer-shift-x))
(define rod2-cy0 (+ rod2-by dimer-shift-y))

; Rotate the whole structure by phi.
(define rod1-cx (rot-x rod1-cx0 rod1-cy0 phi-rad))
(define rod1-cy (rot-y rod1-cx0 rod1-cy0 phi-rad))
(define rod2-cx (rot-x rod2-cx0 rod2-cy0 phi-rad))
(define rod2-cy (rot-y rod2-cx0 rod2-cy0 phi-rad))

(define rod1-ang (+ phi-rad rod1-out-ang))
(define rod2-ang (+ phi-rad rod2-out-ang))

(define rod1-ux (cos rod1-ang))
(define rod1-uy (sin rod1-ang))
(define rod2-ux (cos rod2-ang))
(define rod2-uy (sin rod2-ang))

(define dimer-geometry
  (append
   (make-capsule-rod rod1-cx rod1-cy rod-z rod1-ux rod1-uy 0)
   (make-capsule-rod rod2-cx rod2-cy rod-z rod2-ux rod2-uy 0)))

; ----------------------------
; Spectrum
; ----------------------------
(define fmin (/ 1 wvl_max))
(define fmax (/ 1 wvl_min))
(define fcen (* 0.5 (+ fmin fmax)))
(define df (- fmax fmin))

; ----------------------------
; Source and monitors
; ----------------------------
(define halfz (/ sz 2))
(define monitor-sx (- sx (* 2 dpml) 0.04))
(define monitor-sy (- sy (* 2 dpml) 0.04))

(define src-z  (- halfz dpml 0.10))
(define refl-z (- halfz dpml 0.18))
(define tran-z (+ (- halfz) dpml 0.18))

(set! sources
      (list
       (make source
             (src (make gaussian-src (frequency fcen) (fwidth df)))
             (component Ex)
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

(define decay-point (vector3 0 0 0))
(define flux-tag-refl "wu_fast_refl_ref")
(define flux-tag-tran "wu_fast_tran_ref")

; ----------------------------
; Run
; ----------------------------
(if (= mode 0)
    (begin
      (set! geometry (list))
      (set! refl (add-flux fcen df nfreq refl-region))
      (set! tran (add-flux fcen df nfreq tran-region))
      (run-sources+
       (stop-when-fields-decayed decay_time Ex decay-point decay_tol))
      (save-flux flux-tag-refl refl)
      (save-flux flux-tag-tran tran)
      (display-fluxes refl tran))
    (begin
      (set! geometry dimer-geometry)
      (set! refl (add-flux fcen df nfreq refl-region))
      (set! tran (add-flux fcen df nfreq tran-region))
      (load-minus-flux flux-tag-refl refl)
      (load-minus-flux flux-tag-tran tran)
      (run-sources+
       (stop-when-fields-decayed decay_time Ex decay-point decay_tol))
      (display-fluxes refl tran)))
