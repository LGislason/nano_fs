; ============================================================
; Asymmetric gold metasurface reflectance
; Meep 1.1.1 compatible
; ============================================================

; ----------------------------
; Gold model
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

(define SiO2 (make dielectric (epsilon 2.1)))

; ----------------------------
; Geometry parameters
; ----------------------------

(define period 0.80)
(define bar-length 0.40)
(define bar-width 0.12)
(define metal-z 0.06)

(define cx 0.12)
(define cy -0.10)

; ----------------------------
; Simulation setup
; ----------------------------

(set! geometry-lattice
      (make lattice (size (vector3 period period 4.0))))

(define dpml 1.5)
(define resolution 80)
(define run-time 1200)

(set! pml-layers (list (make pml (thickness dpml) (direction Z))))
(set! ensure-periodicity true)

; ----------------------------
; Spectrum
; ----------------------------

(define wvl-min 0.45)
(define wvl-max 1.20)
(define nfreq 120)

(define fmin (/ 1 wvl-max))
(define fmax (/ 1 wvl-min))
(define fcen (* 0.5 (+ fmin fmax)))
(define df (- fmax fmin))

; ----------------------------
; Geometry placement
; ----------------------------

(define half-length (/ bar-length 2))
(define quarter-length (/ bar-length 4))

(define blk1-x (- cx half-length))
(define blk1-y (+ cy quarter-length))

(define blk2-x (+ cx half-length))
(define blk2-y blk1-y)

(define gold-z-center -1.25)

; ----------------------------
; Geometry
; ----------------------------

(define gold-geometry
  (list
   (make block
         (size (vector3 period period infinity))
         (center (vector3 0 0 -1.8))
         (material SiO2))

   (make block
         (size (vector3 bar-length bar-width metal-z))
         (center (vector3 blk1-x blk1-y gold-z-center))
         (material gold-3term))

   (make block
         (size (vector3 bar-width bar-length metal-z))
         (center (vector3 blk2-x blk2-y gold-z-center))
         (material gold-3term))))
; ============================================================
; Reflectance calculation
; ============================================================
(define refl 0)
(define refl2 0)
(define halfz (/ 4.0 2))

(let*
 (
  (src-z (- halfz dpml 0.2))
  (refl-z (- halfz dpml 0.3))

  (src
   (list
    (make source
          (src (make gaussian-src (frequency fcen) (fwidth df)))
          (component Ey)
          (center (vector3 0 0 src-z))
          (size (vector3 period period 0)))))

  (refl-region
   (make flux-region
         (center (vector3 0 0 refl-z))
         (size (vector3 period period 0))))
 )

 ; EMPTY RUN
 (set! geometry '())
 (set! sources src)

 (set! refl (add-flux fcen df nfreq refl-region))

 (run-until run-time)

 (save-flux "incident" refl)

 ; STRUCTURE RUN
 (set! geometry gold-geometry)
 (set! sources src)

 (set! refl2 (add-flux fcen df nfreq refl-region))

 (load-minus-flux "incident" refl2)

 (run-sources+ 300)
 ; WRITE DATA
 ;(display-fluxes refl)
 (display-fluxes refl2)
)

