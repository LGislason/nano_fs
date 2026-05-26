; ============================================================
; Gold nanosphere validation
; Meep 1.1.1 compatible
;
; mode = 0 incident run: no sphere
; mode = 1 structure run: gold sphere
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
(define sx 1.2)
(define sy 1.2)
(define sz 1.2)
(define dpml 0.25)
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
(define wvl-max 0.90)
(define nfreq 200)

(define fmin (/ 1 wvl-max))
(define fmax (/ 1 wvl-min))
(define fcen (* 0.5 (+ fmin fmax)))
(define df (- fmax fmin))

; ----------------------------
; Source and monitor
; Propagation along z, Ex polarization
; ----------------------------
(define halfz (/ sz 2))
(define src-z (- halfz dpml 0.10))
(define tran-z (+ (- halfz) dpml 0.10))

(set! sources
      (list
       (make source
             (src (make gaussian-src (frequency fcen) (fwidth df)))
             (component Ex)
             (center (vector3 0 0 src-z))
             (size (vector3 (- sx (* 2 dpml)) (- sy (* 2 dpml)) 0)))))

(define tran-region
  (make flux-region
        (center (vector3 0 0 tran-z))
        (size (vector3 (- sx (* 2 dpml)) (- sy (* 2 dpml)) 0))))

(define tran (add-flux fcen df nfreq tran-region))

(if (= mode 0)
    (begin
      (run-sources+ 600)
      (save-flux "sphere_inc_tran" tran)
      (display-fluxes tran))
    (begin
      (load-minus-flux "sphere_inc_tran" tran)
      (run-sources+ 600)
      (display-fluxes tran)))
