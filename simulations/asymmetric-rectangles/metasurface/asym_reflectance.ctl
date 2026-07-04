; ============================================================
; Asymmetric gold nanobar metasurface — reflectance (Meep 1.1.1 Scheme)
;
; Corrected version of the old centered_single.ctl. Fixes:
;   * The metal + substrate are now in the NON-PML interior. (The old file put
;     gold-z-center = -1.25 with sz=5, dpml=1.5, so the interior was only
;     [-1,+1] and the gold sat ~1 um deep inside the bottom PML — it was being
;     absorbed instead of scattering.)
;   * Field-decay stopping instead of a fixed run time, so spectra are not
;     truncation-limited (removes the short-wavelength ringing).
;   * shape = 2 builds the full asymmetric (H + V) unit cell; 0 / 1 are the
;     single-bar controls.
;
; Coordinates: substrate top at z = 0; gold bars sit ON it (z in [0, metal_z]);
; source and reflection monitor are in the vacuum above; a transmission monitor
; sits in the substrate; PML caps the top and bottom.
;
; mode  = 0 empty (vacuum) incident reference ; 1 structure
; pol   = 0 Ex ; 1 Ey
; shape = 0 horizontal bar ; 1 vertical bar ; 2 asymmetric H+V pair
;
; Example:
;   meep-openmpi mode=0 pol=0 asym_reflectance.ctl > inc_ex.txt
;   meep-openmpi mode=1 pol=0 shape=2 asym_reflectance.ctl > asym_ex.txt
; ============================================================

(if (not (defined? 'mode))       (define mode 0))
(if (not (defined? 'pol))        (define pol 0))     ; 0 Ex, 1 Ey
(if (not (defined? 'shape))      (define shape 2))   ; 0 H, 1 V, 2 asymmetric pair
(if (not (defined? 'res))        (define res 80))
(if (not (defined? 'gap))        (define gap 0.10))  ; edge-to-edge gap of the pair
(if (not (defined? 'nfreq))      (define nfreq 150))
(if (not (defined? 'wvl_min))    (define wvl_min 0.50))
(if (not (defined? 'wvl_max))    (define wvl_max 1.20))
(if (not (defined? 'decay_time)) (define decay_time 30))
(if (not (defined? 'decay_tol))  (define decay_tol 1e-4))

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
(define gold-zc (/ metal-z 2))       ; bars sit ON the substrate (bottom at z=0)

; ----------------------------
; Cell / stack  (substrate top at z = 0)
; ----------------------------
(define sz 4.0)
(define dpml 1.0)
(set! dimensions 3)              ; force 3D; old Meep can misinfer periodic cells
(set! geometry-lattice (make lattice (size (vector3 period period sz))))
(set! pml-layers (list (make pml (thickness dpml) (direction Z))))
(set! ensure-periodicity true)
(set! k-point (vector3 0 0 0))   ; normal incidence; required for the periodic (x,y) cell

(define substrate
  (make block (size (vector3 period period 2.0))
        (center (vector3 0 0 -1.0)) (material SiO2)))  ; z in [-2, 0]

(define (hbar cx cy)
  (make block (size (vector3 bar-length bar-width metal-z))
        (center (vector3 cx cy gold-zc)) (material gold-3term)))
(define (vbar cx cy)
  (make block (size (vector3 bar-width bar-length metal-z))
        (center (vector3 cx cy gold-zc)) (material gold-3term)))

; asymmetric pair: horizontal bar and vertical bar, edge-to-edge gap = gap
(define center-sep (+ (/ bar-length 2) (/ bar-width 2) gap))
(define bars
  (cond ((= shape 0) (list (hbar 0 0)))
        ((= shape 1) (list (vbar 0 0)))
        (else (list (hbar (* -0.5 center-sep) 0)
                    (vbar (*  0.5 center-sep) 0)))))
(define structure (append (list substrate) bars))

; ----------------------------
; Spectrum / source
; ----------------------------
(define fmin (/ 1 wvl_max))
(define fmax (/ 1 wvl_min))
(define fcen (* 0.5 (+ fmin fmax)))
(define df (- fmax fmin))
(define src-comp (if (= pol 0) Ex Ey))

(set! resolution res)

(define halfz (/ sz 2))
(define src-z  (- halfz dpml 0.30))       ; 0.70  (vacuum, above the bars)
(define refl-z (- halfz dpml 0.50))       ; 0.50  (between source and structure)
(define tran-z (+ (- halfz) dpml 0.40))   ; -0.60 (inside the substrate)

(set! sources
      (list (make source
                  (src (make gaussian-src (frequency fcen) (fwidth df)))
                  (component src-comp)
                  (center (vector3 0 0 src-z))
                  (size (vector3 period period 0)))))

(define refl-region (make flux-region (center (vector3 0 0 refl-z)) (size (vector3 period period 0))))
(define tran-region (make flux-region (center (vector3 0 0 tran-z)) (size (vector3 period period 0))))
(define refl 0)
(define tran 0)
(define decay-pt (vector3 0 0 refl-z))
(define tag-r (string-append "asym_refl_ref_p" (if (= pol 0) "0" "1")))
(define tag-t (string-append "asym_tran_ref_p" (if (= pol 0) "0" "1")))

(print "asym-metasurface: mode=" mode " shape=" shape " pol=" pol
       " gap=" gap " res=" res " wvl=[" wvl_min "," wvl_max "]\n")

; ----------------------------
; Run
; ----------------------------
(if (= mode 0)
    (begin
      (set! geometry (list))                 ; vacuum reference = incident power
      (set! refl (add-flux fcen df nfreq refl-region))
      (set! tran (add-flux fcen df nfreq tran-region))
      (run-sources+
       (stop-when-fields-decayed decay_time src-comp decay-pt decay_tol))
      (save-flux tag-r refl)
      (save-flux tag-t tran)
      (display-fluxes refl tran))
    (begin
      (set! geometry structure)
      (set! refl (add-flux fcen df nfreq refl-region))
      (set! tran (add-flux fcen df nfreq tran-region))
      (load-minus-flux tag-r refl)           ; subtract incident -> reflected only
      (run-sources+
       (stop-when-fields-decayed decay_time src-comp decay-pt decay_tol))
      (display-fluxes refl tran)))
