; ============================================================
; Normalized Reflectance – Asymmetric Gold Bars on SiO2
; Meep 1.1.1 compatible
; ============================================================

; ----------------------------
; Lattice
; ----------------------------
(set! geometry-lattice (make lattice (size 0.80 0.80 7.06)))
(set! resolution 100)

(set! pml-layers (list (make pml (thickness 1.5) (direction Z))))

; ============================================================
; Materials
; ============================================================

(define sio2 (make dielectric (index 1.45)))

(define omega_meep (* 2 pi 2.99792458e14))

(define au_epsilon_inf 4.8929)

(define au_omega_d (/ 1.2944e16 omega_meep))
(define au_omega_1 (/ 1.3617e15 omega_meep))
(define au_omega_2 (/ 4.1636e15 omega_meep))
(define au_omega_3 (/ 5.0753e15 omega_meep))

(define au_gamma_d (/ 1.0003e9 omega_meep))
(define au_gamma_1 (/ (* 2 4.7356e14) omega_meep))
(define au_gamma_2 (/ (* 2 4.4931e14) omega_meep))
(define au_gamma_3 (/ (* 2 5.8469e14) omega_meep))

(define au_sigma_1 4.7282)
(define au_sigma_2 0.72996)
(define au_sigma_3 1.5103)

(define gold
  (make dielectric
    (epsilon au_epsilon_inf)
    (E-polarizations
      (make polarizability
        (omega 1e-20)
        (gamma au_gamma_d)
        (sigma (* (* 1e20 au_omega_d)
                  (* 1e20 au_omega_d))))
      (make polarizability (omega au_omega_1) (gamma au_gamma_1) (sigma au_sigma_1))
      (make polarizability (omega au_omega_2) (gamma au_gamma_2) (sigma au_sigma_2))
      (make polarizability (omega au_omega_3) (gamma au_gamma_3) (sigma au_sigma_3)))))

; ============================================================
; Broadband Source
; ============================================================

(define fmin (/ 1 1.20))
(define fmax (/ 1 0.45))
(define fcen (/ (+ fmin fmax) 2))
(define df (- fmax fmin))
(define nfreq 120)

(define pol Ex)  ; CHANGE TO Ey FOR EY RUN

(set! sources
  (list
    (make source
      (src (make gaussian-src (frequency fcen) (fwidth df)))
      (component pol)
      (center 0 0 2.23)
      (size 0.80 0.80 0))))

(define refl-region
  (make flux-region
    (center 0 0 1.53)
    (size 0.80 0.80 0)))

; ============================================================
; 1) EMPTY RUN (incident flux)
; ============================================================

(set! geometry
  (list
    (make block (size infinity infinity 2.0)
                (center 0 0 -1.0)
                (material sio2))))

(define refl-empty (add-flux fcen df nfreq refl-region))

(run-until 800)

(define inc-flux (get-fluxes refl-empty))
(define empty-flux-data (get-flux-data refl-empty))

; ============================================================
; 2) STRUCTURE RUN
; ============================================================

(reset-meep)

(set! geometry
  (list
    (make block (size infinity infinity 2.0)
                (center 0 0 -1.0)
                (material sio2))

    (make block (size 0.40 0.12 0.06)
                (center -0.08 0.00 0.03)
                (material gold))

    (make block (size 0.12 0.40 0.06)
                (center 0.32 0.00 0.03)
                (material gold))))

(set! sources
  (list
    (make source
      (src (make gaussian-src (frequency fcen) (fwidth df)))
      (component pol)
      (center 0 0 2.23)
      (size 0.80 0.80 0))))

(define refl-struct (add-flux fcen df nfreq refl-region))

(load-minus-flux-data refl-struct empty-flux-data)

(run-until 1200)

(define refl-flux (get-fluxes refl-struct))
(define freqs (get-flux-freqs refl-struct))

; ============================================================
; Output R = -refl / inc
; ============================================================

(define out (open-output-file "R_output.txt"))

(define (loop i)
  (if (< i nfreq)
      (begin
        (define f (list-ref freqs i))
        (define wvl (/ 1 f))
        (define inc (list-ref inc-flux i))
        (define rfl (list-ref refl-flux i))
        (define R (/ (- rfl) inc))
        (fprintf out "~g ~g ~g\n" f wvl R)
        (loop (+ i 1)))))

(loop 0)
(close-output-port out)

(print "Normalized reflectance written to R_output.txt\n")
