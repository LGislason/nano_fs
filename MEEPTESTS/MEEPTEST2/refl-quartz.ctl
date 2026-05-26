(set-param! resolution 200) ; pixels/um

(define-param sz 10)
(set! geometry-lattice (make lattice (size no-size no-size sz)))
(set! dimensions 1)

; wavelength range: 0.4 to 0.8 um
(define lambda-min 0.4)
(define lambda-max 0.8)
(define fmax (/ 1 lambda-min))
(define fmin (/ 1 lambda-max))
(define fcen (* 0.5 (+ fmax fmin)))
(define df (- fmax fmin))

(define dpml 1.0)
(set! pml-layers (list (make pml (thickness dpml))))

(set! k-point (vector3 0 0 0))

; source placed away from PML
(set! sources
      (list
       (make source
             (src (make gaussian-src (frequency fcen) (fwidth df)))
             (component Ex)
             (center 0 0 -3.5))))

; Run first with true, then second with false
(define-param empty? true)

; air/quartz interface at z = 0
; quartz modeled as constant epsilon = n^2 = 1.45^2 = 2.1025
(if (not empty?)
    (set! geometry
          (list
           (make block
                 (size infinity infinity (* 0.5 sz))
                 (center 0 0 (* 0.25 sz))
                 (material (make dielectric (epsilon 2.1025)))))))

(define nfreq 50)

; reflection monitor between source and interface
(define refl
  (add-flux fcen df nfreq
            (make flux-region (center 0 0 -2.5))))

; subtract incident fields during structure run
(if (not empty?) (load-minus-flux "refl-flux" refl))

(run-sources+
 (stop-when-fields-decayed 50 Ex (vector3 0 0 -3.5) 1e-9))

; save incident flux during empty run
(if empty? (save-flux "refl-flux" refl))

(display-fluxes refl)
