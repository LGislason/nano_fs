; ============================================================
; Asymmetric "- |" Gold Unit Cell on SiO2 Substrate
; Meep 1.1.1 / Guile 1.8 Compatible
; Fully normalized reflectance EX + EY
; ============================================================

; ----------------------------
; Gold 3-term Drude-Lorentz
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
               (sigma (* (* 1e20 au_omega_d) (* 1e20 au_omega_d))))
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
               (sigma au_sigma_3))
        )))

; ----------------------------
; Define SiO2 substrate
; ----------------------------

(define SiO2 (make dielectric (epsilon 2.1)))

; ----------------------------
; Simulation parameters
; ----------------------------

(define period 0.80)
(define bar-length 0.40)
(define bar-width 0.12)
(define metal-z 0.06)

(define cx 0.12)
(define cy -0.10)

(set! geometry-lattice
     (make lattice (size (vector3 period period 6.0))))

(define dpml 1.0)
(define resolution 60)
(define run-time 600)

(set! pml-layers (list (make pml (thickness dpml) (direction Z))))

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
; Geometry positioning
; ----------------------------

(define half-length (/ bar-length 2))
(define quarter-length (/ bar-length 4))

(define blk1-x (- cx half-length))
(define blk1-y (+ cy quarter-length))

(define blk2-x (+ cx half-length))
(define blk2-y (+ cy quarter-length))

; Substrate top at z=-2.0
; Gold center = -2.0 + metal-z/2 = -1.97
(define gold-z-center -1.97)

; ----------------------------
; Geometry list
; ----------------------------

(define gold-geometry
  (list
   (make block
         (size   (vector3 period period 1.0))
         (center (vector3 0 0 -2.5))
         (material SiO2))

   (make block
         (size   (vector3 bar-length bar-width metal-z))
         (center (vector3 blk1-x blk1-y gold-z-center))
         (material gold-3term))

   (make block
         (size   (vector3 bar-width bar-length metal-z))
         (center (vector3 blk2-x blk2-y gold-z-center))
         (material gold-3term))
))

; ----------------------------
; Reflectance function
; ----------------------------

(define halfz (/ 6.0 2))

(define (run-reflectance pol tag)

  (let* (
         (src-z (- halfz dpml 0.2))
         (refl-z (- halfz dpml 0.8))

         (src
          (list
           (make source
                 (src (make gaussian-src (frequency fcen) (fwidth df)))
                 (component pol)
                 (center (vector3 0 0 src-z))
                 (size   (vector3 period period 0)))))

         (refl-region
          (make flux-region
                (center (vector3 0 0 refl-z))
                (size   (vector3 period period 0))))
        )

    ; EMPTY CELL RUN
    (set! geometry '())
    (set! sources src)

    (let* ((flux-empty (add-flux fcen df nfreq refl-region)))
      (run-until run-time)

      (let* (
             (empty-data (get-flux-data flux-empty))
             (inc-flux   (get-fluxes flux-empty))
             (freqs      (get-flux-freqs flux-empty))
            )

        ; STRUCTURE RUN
        (set! geometry gold-geometry)
        (set! sources src)

        (let* ((flux-struct (add-flux fcen df nfreq refl-region)))
          (load-minus-flux-data flux-struct empty-data)
          (run-until run-time)

          (let* (
                 (refl (get-fluxes flux-struct))
                 (wl '())
                 (rr '())
                )

            (do ((i 0 (+ i 1))) ((= i nfreq))
              (let* ((freq (list-ref freqs i))
                     (wvl (/ 1 freq))
                     (inc (list-ref inc-flux i))
                     (ref (list-ref refl i))
                     (R (/ (abs ref) (max (abs inc) 1e-30))))
                (set! wl (append wl (list wvl)))
                (set! rr (append rr (list R)))))

            (list wl rr)
          ))))))

; ----------------------------
; Run EX and EY
; ----------------------------

(define res-Ex (run-reflectance Ex "Ex"))
(define w-ex (car res-Ex))
(define R-ex (cadr res-Ex))

(define res-Ey (run-reflectance Ey "Ey"))
(define w-ey (car res-Ey))
(define R-ey (cadr res-Ey))

; ----------------------------
; Save results
; ----------------------------

(with-output-to-file "reflectance_ex.dat"
  (lambda ()
    (display "# wavelength  R_ex\n")
    (do ((i 0 (+ i 1))) ((= i nfreq))
      (display (list-ref w-ex i)) (display "   ")
      (display (list-ref R-ex i)) (newline))))

(with-output-to-file "reflectance_ey.dat"
  (lambda ()
    (display "# wavelength  R_ey\n")
    (do ((i 0 (+ i 1))) ((= i nfreq))
      (display (list-ref w-ey i)) (display "   ")
      (display (list-ref R-ey i)) (newline))))

