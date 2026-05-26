;Asymmetric "-|" array reflectance
; Built to match nano-pillar.ctl style exactly for Meep 1.1.1
; Two-pass workflow with objects? false/true
; ============================================================

; ---------------- switches ----------------
(define-param objects? false)
(define-param pol 0) ; 0 -> Ex, 1 -> Ey

; ---------------- units ----------------
; nano-pillar.ctl uses cell=1 nm units
(define-param cell 1)              ; 1 nm


; ---------------- geometry params (nm) ----------------
(define-param pitch 800)           ; period (nm)
(define-param sy    2000)          ; propagation span (nm) like nano-pillar
(define-param tAu   60)            ; gold thickness along Y (nm)
(define-param L     400)           ; long bar length (nm)
(define-param W     120)           ; bar width (nm)

; offsets in the periodic plane (X-Z) in nm
(define-param ox 120)              ; x offset
(define-param oz -100)             ; z offset

; ---------------- resolution/PML ----------------
(define-param res 1)
(define-param dpml 12)

; ---------------- spectrum (match nano-pillar approach) ----------------
; choose a center wavelength and a bandwidth via gaussian_width + flux_width
(define-param wavelen 800)         ; center wavelength in nm (within gold model range)
(define wavelength-c (/ wavelen cell))

(define-param gaussian_width 700)  ; nm, controls time pulse width
(define gaussian_width-c (/ gaussian_width cell))

(define-param flux_width 1000)     ; nm span you want sampled (roughly)
(define flux_width-c (/ flux_width cell))

(define freq (/ 1 wavelength-c))
(define df   (/ 1 gaussian_width-c))
(define dflux (/ 1 flux_width-c))
(define-param nfreq 500)

; ---------------- derived cell sizes in Meep units ----------------
(define sx pitch)
(define sz pitch)

(define sx-c (/ sx cell))
(define sy-c (/ sy cell))
(define sz-c (/ sz cell))

(define ypml-c dpml)

(set! dimensions 3)
(set! resolution res)

; nano-pillar lattice style: numeric size for lattice
(set! geometry-lattice (make lattice (size sx-c sy-c sz-c)))

; propagate along +Y, PML only in Y (same as nano-pillar)
(set! pml-layers
  (list (make pml (thickness ypml-c) (direction Y))))

(set! ensure-periodicity true)

; ---------------- gold-3term model (COPIED from nano-pillar.ctl) ----------------
(define um 1e-6)
(define c0 299792458.0)
(define omega_meep (* 2 pi (/ c0 um)))

(define au_epsilon_inf 4.8929e+00)

(define au_omega_d (/ 1.2944e+16 omega_meep))
(define au_omega_1 (/ 1.3617e+15 omega_meep))
(define au_omega_2 (/ 4.1636e+15 omega_meep))
(define au_omega_3 (/ 5.0753e+15 omega_meep))

(define au_gamma_d (/ 1.0003e+09 omega_meep))
(define au_gamma_1 (/ (* 4.7356e+14 2) omega_meep))
(define au_gamma_2 (/ (* 4.4931e+14 2) omega_meep))
(define au_gamma_3 (/ (* 5.8469e+14 2) omega_meep))

(define au_sigma_1 4.7282e+00)
(define au_sigma_2 7.2996e-01)
(define au_sigma_3 1.5103e+00)

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
        (sigma au_sigma_3)))))

; ---------------- source + monitor (match nano-pillar conventions) ----------------
(define src-comp (if (= pol 0) Ex Ey))

; source near bottom Y boundary
(define src-y (+ 1 (* -0.5 sy-c)))

; reflection plane above source
(define refl-y (+ 5 (* -0.5 sy-c)))

(set! sources
  (list
    (make source
      (src (make gaussian-src (frequency freq) (fwidth df)))
      (component src-comp)
      (center (vector3 0 src-y 0))
      (size   (vector3 sx-c 0 sz-c)))))

(define flux1
  (add-flux freq dflux nfreq
    (make flux-region
      (size   (vector3 sx-c 0 sz-c))
      (center (vector3 0 refl-y 0)))))

; output directory fixed so BOTH passes share the same prefix
(use-output-directory "meep-out-asym")

(define tag (if (= pol 0) "Ex" "Ey"))
(define flux-tag (string-append "ASYM_" tag))

; ---------------- geometry: only when objects? true ----------------
; Place "-|" in X-Z plane, thickness along Y = tAu
(define tAu-c (/ tAu cell))
(define L-c   (/ L   cell))
(define W-c   (/ W   cell))
(define ox-c  (/ ox  cell))
(define oz-c  (/ oz  cell))

(define halfL (/ L-c 2))
(define quarterL (/ L-c 4))

; centers (precomputed, no math inside vector3 beyond symbols)
(define x1 (- ox-c halfL))
(define z1 (+ oz-c quarterL))
(define x2 (+ ox-c halfL))
(define z2 z1)

(if objects?
  (set! geometry
    (list
      ; horizontal bar (long in X, narrow in Z)
      (make block
        (center (vector3 x1 0 z1))
        (size   (vector3 L-c tAu-c W-c))
        (material gold-3term))
      ; vertical bar (long in Z, narrow in X)
      (make block
        (center (vector3 x2 0 z2))
        (size   (vector3 W-c tAu-c L-c))
        (material gold-3term)))))

; ---------------- nano-pillar reflectance workflow ----------------
; subtract incident fields at the reflection monitor
(if objects? (load-minus-flux flux-tag flux1))

(run-sources+
  (stop-when-fields-decayed 50 src-comp (vector3 0 0 0) 1e-6))

(if (not objects?) (save-flux flux-tag flux1))

; write out spectrum each run for easy normalization
(define P (get-fluxes flux1))
(define freqs (get-flux-freqs flux1))
(define outname
  (if objects?
      (string-append "reflected_" tag ".txt")
      (string-append "incident_" tag ".txt")))

(with-output-to-file outname
  (lambda ()
    (display "lambda_nm  P\n")
    (let loop ((i 0))
      (if (< i (length freqs))
        (let* ((f (list-ref freqs i))
               (lam (/ 1 f))          ; because cell=1 nm -> lambda in nm
               (pi (list-ref P i)))
          (display lam) (display "  ") (display pi) (display "\n")
          (loop (+ i 1)))
        'done))))

(display-fluxes flux1)
(print "Wrote " outname)
