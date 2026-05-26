; Asymmetric gold metasurface reflectance
; Meep 1.1.1 compatible
;
; mode = 0  -> empty incident run
; mode = 1  -> structure run
;
; pol  = 0  -> Ex polarization
; pol  = 1  -> Ey polarization
; shape = 0 horizontal only
; shape = 1 vertical only

; Example:
;   meep-openmpi mode=0 pol=0 asym_param.ctl > inc_ex.log
;   meep-openmpi mode=1 pol=0 asym_param.ctl > refl_ex.log
;   meep-openmpi mode=0 pol=1 asym_param.ctl > inc_ey.log
;   meep-openmpi mode=1 pol=1 asym_param.ctl > refl_ey.log
; ============================================================

; ----------------------------
; Default run parameters
; ----------------------------
(if (not (defined? 'mode)) (define mode 0))
(if (not (defined? 'pol))  (define pol 1))
(if (not (defined? 'shape)) (define shape 0))

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

(define gap 0.10) ; 100 nm edge-to-edge gap

(define bar-y 0.0)

; ----------------------------
; Rectangle positions
; ----------------------------
(define half-x-horiz (/ bar-length 2))
(define half-x-vert  (/ bar-width 2))

; center-to-center separation for exact edge gap
(define center-sep (+ half-x-horiz half-x-vert gap))

(define blk1-x (/ (- center-sep) 2))
(define blk2-x (/ center-sep 2))

(define blk1-y bar-y)
(define blk2-y bar-y)

(define gold-z-center -1.25)


; ----------------------------
; Simulation setup
; ----------------------------
(define sz 5.0)
(define dpml 1.5)
(define resolution 100)

(set! geometry-lattice
      (make lattice (size (vector3 period period sz))))

(set! pml-layers
      (list (make pml (thickness dpml) (direction Z))))

(set! ensure-periodicity true)

; ----------------------------
; Spectrum
; ----------------------------
(define wvl-min 0.50)
(define wvl-max 1.20)
(define nfreq 120)

(define fmin (/ 1 wvl-max))
(define fmax (/ 1 wvl-min))
(define fcen (* 0.5 (+ fmin fmax)))
(define df (- fmax fmin))



; ----------------------------
; Geometry
; ----------------------------

(define substrate-block
  (make block
        (size (vector3 period period 3.0))
        (center (vector3 0 0 -1.75))
        (material SiO2)))

(define horizontal-bar
  (make block
        (size (vector3 bar-length bar-width metal-z))
        (center (vector3 0 0 gold-z-center))
        (material gold-3term)))

(define vertical-bar
  (make block
        (size (vector3 bar-width bar-length metal-z))
        (center (vector3 0 0 gold-z-center))
        (material gold-3term)))

(define gold-geometry
  (if (= shape 0)
      (list substrate-block horizontal-bar)
      (list substrate-block vertical-bar)))

; ----------------------------
; Polarization selection
; ----------------------------
(define src-comp
  (if (= pol 0) Ex Ey))

(define flux-tag
  (if (= pol 0) "incident_ex" "incident_ey"))

; ----------------------------
; Reflectance setup
; ----------------------------
(define halfz (/ sz 2))
(define src-z (- halfz dpml 0.2))
(define refl-z (- halfz dpml 0.35))

(set! sources
      (list
       (make source
             (src (make gaussian-src (frequency fcen) (fwidth df)))
             (component src-comp)
             (center (vector3 0 0 src-z))
             (size (vector3 period period 0)))))

(define refl-region
  (make flux-region
        (center (vector3 0 0 refl-z))
        (size (vector3 period period 0))))

(define refl 0)

; ----------------------------
; Run selection
; ----------------------------
(if (= mode 0)
    (begin
      ; empty incident run
      (set! geometry '())
      (set! refl (add-flux fcen df nfreq refl-region))
      (run-sources+ 600)
      (save-flux flux-tag refl)
      (display-fluxes refl))
    (begin
      ; structure run
      (set! geometry gold-geometry)
      (set! refl (add-flux fcen df nfreq refl-region))
      (load-minus-flux flux-tag refl)
      (run-sources+ 1200)
      (display-fluxes refl)))
