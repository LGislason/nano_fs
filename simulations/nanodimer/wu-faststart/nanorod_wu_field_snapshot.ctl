; ============================================================
; Nanorod dimer field snapshot for old Scheme Meep clusters
;
; This uses the same simplified homogeneous-background geometry as
; nanorod_wu_fast.ctl, but runs a single continuous-wave wavelength and writes
; 2D HDF5 field slices for local plotting.
;
; Example:
; mpirun -np 4 meep-openmpi \
;   res=300 gap=0.006 theta=80 phi=40 wvl=0.672 run_time=300 \
;   nanorod_wu_field_snapshot.ctl > run.log
; ============================================================

; ----------------------------
; Defaults
; ----------------------------
(if (not (defined? 'theta))    (define theta 80.0))
(if (not (defined? 'phi))      (define phi 40.0))
(if (not (defined? 'res))      (define res 300))
(if (not (defined? 'gap))      (define gap 0.006))
(if (not (defined? 'nbg))      (define nbg 1.28))
(if (not (defined? 'wvl))      (define wvl 0.672))
(if (not (defined? 'run_time)) (define run_time 300))
(if (not (defined? 'pol))      (define pol 0)) ; 0 Ex, 1 Ey

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
; Dimer geometry
; ----------------------------
(define rod-length 0.069)
(define rod-width 0.024)
(define rod-radius (/ rod-width 2))
(define rod-body-length (- rod-length (* 2 rod-radius)))
(define rod-z 0.0)

(define theta-rad (deg->rad theta))
(define phi-rad (deg->rad phi))

; Base dimer: theta is the internal opening angle between rods.  The
; user-facing gap is the closest metal-to-metal surface separation.  The
; rounded capsule end centers must therefore be farther apart by 2 radii.
(define rod1-out-ang pi)
(define rod2-out-ang (- pi theta-rad))
(define gap-axis-ang (- (* 0.5 (+ rod1-out-ang rod2-out-ang)) (/ pi 2)))
(define gap-axis-x (cos gap-axis-ang))
(define gap-axis-y (sin gap-axis-ang))
(define rod-tip-gap (+ gap (* 2 rod-radius (- 1 (sin (/ theta-rad 2))))))

(define rod1-tip-x (* -0.5 rod-tip-gap gap-axis-x))
(define rod1-tip-y (* -0.5 rod-tip-gap gap-axis-y))
(define rod2-tip-x (* 0.5 rod-tip-gap gap-axis-x))
(define rod2-tip-y (* 0.5 rod-tip-gap gap-axis-y))

(define rod1-bx (+ rod1-tip-x (* 0.5 rod-length (cos rod1-out-ang))))
(define rod1-by (+ rod1-tip-y (* 0.5 rod-length (sin rod1-out-ang))))
(define rod2-bx (+ rod2-tip-x (* 0.5 rod-length (cos rod2-out-ang))))
(define rod2-by (+ rod2-tip-y (* 0.5 rod-length (sin rod2-out-ang))))

; Keep the gap reference centered at the cell origin so centered field slices
; pass through the dimer junction.
(define dimer-shift-x 0.0)
(define dimer-shift-y 0.0)

(define rod1-cx0 (+ rod1-bx dimer-shift-x))
(define rod1-cy0 (+ rod1-by dimer-shift-y))
(define rod2-cx0 (+ rod2-bx dimer-shift-x))
(define rod2-cy0 (+ rod2-by dimer-shift-y))

; Rotate the whole dimer by phi relative to the fixed source polarization.
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

(set! geometry
      (append
       (make-capsule-rod rod1-cx rod1-cy rod-z rod1-ux rod1-uy 0)
       (make-capsule-rod rod2-cx rod2-cy rod-z rod2-ux rod2-uy 0)))

; ----------------------------
; Source
; ----------------------------
(define freq (/ 1 wvl))
(define src-comp (if (= pol 0) Ex Ey))

(define halfz (/ sz 2))
(define field-sx (- sx (* 2 dpml) 0.04))
(define field-sy (- sy (* 2 dpml) 0.04))
(define field-sz (- sz (* 2 dpml) 0.04))
(define src-z (- halfz dpml 0.10))

(set! sources
      (list
       (make source
             (src (make continuous-src (frequency freq)))
             (component src-comp)
             (center (vector3 0 0 src-z))
             (size (vector3 field-sx field-sy 0)))))

; ----------------------------
; Field-slice volumes
; ----------------------------
(define xy-rod-volume
  (volume (center (vector3 0 0 rod-z))
          (size (vector3 field-sx field-sy 0))))

(define xy-above-volume
  (volume (center (vector3 0 0 (+ rod-z 0.030)))
          (size (vector3 field-sx field-sy 0))))

(define xz-volume
  (volume (center (vector3 0 0 0))
          (size (vector3 field-sx 0 field-sz))))

; ----------------------------
; Run and output
; ----------------------------
(print "field-snapshot: wvl=" wvl " phi=" phi " theta=" theta
       " gap=" gap " rod-tip-gap=" rod-tip-gap
       " res=" res " nbg=" nbg "\n")

(run-until run_time
  (at-end
   (to-appended "xy_rod_eps" (in-volume xy-rod-volume output-epsilon)))
  (at-end
   (to-appended "xy_rod_ex" (in-volume xy-rod-volume output-efield-x)))
  (at-end
   (to-appended "xy_rod_ey" (in-volume xy-rod-volume output-efield-y)))
  (at-end
   (to-appended "xy_rod_ez" (in-volume xy-rod-volume output-efield-z)))

  (at-end
   (to-appended "xy_above_eps" (in-volume xy-above-volume output-epsilon)))
  (at-end
   (to-appended "xy_above_ex" (in-volume xy-above-volume output-efield-x)))
  (at-end
   (to-appended "xy_above_ey" (in-volume xy-above-volume output-efield-y)))
  (at-end
   (to-appended "xy_above_ez" (in-volume xy-above-volume output-efield-z)))

  (at-end
   (to-appended "xz_eps" (in-volume xz-volume output-epsilon)))
  (at-end
   (to-appended "xz_ex" (in-volume xz-volume output-efield-x)))
  (at-end
   (to-appended "xz_ey" (in-volume xz-volume output-efield-y)))
  (at-end
   (to-appended "xz_ez" (in-volume xz-volume output-efield-z))))
