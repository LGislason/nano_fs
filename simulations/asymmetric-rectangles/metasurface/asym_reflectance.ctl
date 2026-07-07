; ============================================================
; Asymmetric gold nanobar metasurface — reflectance (Meep 1.1.1)
;
; Working two-pass version (substrate reference -> full metasurface) built on the
; asym.ctl structure that this old Meep can construct. Notes:
;   * Source is fixed to Ey. This old Meep will not build fields for an Ex source
;     in the periodic cell, so the Ex response is obtained by rotating the
;     geometry 90 deg (valid on the square period-by-period lattice): pol=0 ("Ex")
;     rotates [(x,y)->(-y,x), H<->V]; pol=1 ("Ey") uses the geometry as designed.
;   * Reflection is DIFFERENTIAL: load-minus-flux subtracts the bare-substrate
;     reflection, isolating the metal's (polarization-dependent) contribution.
;
; pol   = 0 Ex ; 1 Ey        shape = 0 H ; 1 V ; 2 asymmetric pair
; Example:  meep-openmpi pol=0 shape=2 res=80 asym_reflectance.ctl > s2_p0.txt
; ============================================================

(if (not (defined? 'pol))        (define pol 0))     ; 0 Ex, 1 Ey
(if (not (defined? 'shape))      (define shape 2))   ; 0 H, 1 V, 2 asymmetric pair
(if (not (defined? 'res))        (define res 80))
(if (not (defined? 'gap))        (define gap 0.10))  ; edge-to-edge gap of the pair
(if (not (defined? 'nfreq))      (define nfreq 150))
(if (not (defined? 'wvl_min))    (define wvl_min 0.50))
(if (not (defined? 'wvl_max))    (define wvl_max 1.20))

; ----------------------------
; Gold: Drude + 3 Lorentz (Mie-validated model)
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
         (make polarizability (omega 1e-20) (gamma au_gamma_d)
               (sigma (* (* 1e20 au_omega_d) (* 1e20 au_omega_d))))
         (make polarizability (omega au_omega_1) (gamma au_gamma_1) (sigma au_sigma_1))
         (make polarizability (omega au_omega_2) (gamma au_gamma_2) (sigma au_sigma_2))
         (make polarizability (omega au_omega_3) (gamma au_gamma_3) (sigma au_sigma_3)))))
(define SiO2 (make dielectric (epsilon 2.1)))  ; n ~ 1.45

; ----------------------------
; Geometry parameters
; ----------------------------
(define period 0.80)
(define bar-length 0.40)
(define bar-width 0.12)
(define metal-z 0.06)
(define gold-zc 0.03)       ; sit directly on substrate top

; ----------------------------
; Cell / stack setup
; ----------------------------
(define sz 4.0)
(define dpml 1.5)
(define resolution res)
(define run-time 1200)
(set! geometry-lattice (make lattice (size (vector3 period period sz))))
(set! pml-layers (list (make pml (thickness dpml) (direction Z))))
(set! ensure-periodicity true)

; Substrate fills the bottom half of the cell (semi-infinite)
(define substrate
  (make block (size (vector3 period period infinity))
        (center (vector3 0 0 -1.0)) (material SiO2)))

(define (Hbar x y)   ; bar long along X
  (make block (size (vector3 bar-length bar-width metal-z))
        (center (vector3 x y gold-zc)) (material gold-3term)))
(define (Vbar x y)   ; bar long along Y
  (make block (size (vector3 bar-width bar-length metal-z))
        (center (vector3 x y gold-zc)) (material gold-3term)))

(define center-sep (+ (/ bar-length 2) (/ bar-width 2) gap))
(define off (* 0.5 center-sep))

; polarization encoded as a 90 deg geometry rotation
(define bars
  (if (= pol 1)
      (cond ((= shape 0) (list (Hbar 0 0)))
            ((= shape 1) (list (Vbar 0 0)))
            (else (list (Hbar (- off) 0) (Vbar off 0))))
      (cond ((= shape 0) (list (Vbar 0 0)))
            ((= shape 1) (list (Hbar 0 0)))
            (else (list (Vbar 0 (- off)) (Hbar 0 off))))))
(define metasurface-geometry (append (list substrate) bars))

; ----------------------------
; Spectrum
; ----------------------------
(define fmin (/ 1 wvl_max))
(define fmax (/ 1 wvl_min))
(define fcen (* 0.5 (+ fmin fmax)))
(define df (- fmax fmin))

; ----------------------------
; Two-pass run
; ----------------------------
(define refl 0)  (define refl2 0)
(define tran 0)  (define tran2 0)
(define halfz (/ sz 2))
(define tag-refl (string-append "inc_refl_p" (number->string pol) "_s" (number->string shape)))
(define tag-tran (string-append "inc_tran_p" (number->string pol) "_s" (number->string shape)))

(let* ((src-z  (- halfz dpml 0.2))          ; above the structure, in vacuum
       (refl-z (- halfz dpml 0.3))          ; between source and structure
       (tran-z (- (- halfz) dpml 0.3))      ; inside the substrate, below structure
       (src (list (make source
                        (src (make gaussian-src (frequency fcen) (fwidth df)))
                        (component Ey)
                        (center (vector3 0 0 src-z))
                        (size (vector3 period period 0)))))
       (refl-region (make flux-region (center (vector3 0 0 refl-z)) (size (vector3 period period 0))))
       (tran-region (make flux-region (center (vector3 0 0 tran-z)) (size (vector3 period period 0)))))

  ; --- PASS 1: bare-substrate reference (incident power + substrate reflection) ---
  (set! geometry (list substrate))
  (set! sources src)
  (set! refl (add-flux fcen df nfreq refl-region))
  (set! tran (add-flux fcen df nfreq tran-region))
  (run-until run-time)
  (save-flux tag-refl refl)
  (save-flux tag-tran tran)

  ; clear fields + structure before rebuilding with the metasurface
  (reset-meep)

  ; --- PASS 2: full metasurface ---
  (set! geometry metasurface-geometry)
  (set! sources src)
  (set! refl2 (add-flux fcen df nfreq refl-region))
  (set! tran2 (add-flux fcen df nfreq tran-region))
  (load-minus-flux tag-refl refl2)          ; reflection relative to bare substrate
  (run-sources+ 300)

  (print "\n--- BEGIN REFLECTED FLUX ---\n")
  (display-fluxes refl2)
  (print "--- END REFLECTED FLUX ---\n")
  (print "\n--- BEGIN TRANSMITTED FLUX ---\n")
  (display-fluxes tran2)
  (print "--- END TRANSMITTED FLUX ---\n"))
