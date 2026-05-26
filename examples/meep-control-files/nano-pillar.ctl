
; %----------------------------------------------------------------------------%
; |                                                                            |
; |      1. DEFINE BASIC UNITS AND CONSTANTS                                   |
; |                                                                            |
; %----------------------------------------------------------------------------%
; %----------------------------------------------------------------------------%
; |                                                                            |
; | -----DEFINE MEEP UNITS-----                                                |
; |                                                                            |
; | In meep, units are defined such that the speed of light, c = 1.  Thus,     | 
; | velocities are defined relative to the speed of light in units of c.       |
; | From the relation 'velocity = distance / time' we can define the basic     |
; | units of length and time within the simulation.  Since c=1, we see that:   |
; |   c = a / t => a = t in meep units, i.e. one unit of time is equivalent    |
; |   to one unit of distance.                                                 |
; | To convert this to real world units, we choose our basic unit of           |
; | length, a.  Now, our basic unit of time is related to our chosen unit of   |
; | length thusly, t = a / c.                                                  |
; | Similarly, the basic unit of frequency in meep is just the inverse of      |
; | our basic unit of time:  omega is specified in units of 2*pi*(c / a).      |
; | Thus, it is not necessary to multiply frequencies by 2*pi.  Note also      |
; | that w_vacuum / w_meep equal to f_vacuum / f_meep.  Thus our               |
; | dimensionless frequency is the same with regard to either frequency        |
; | or angular frequency.                                                      |
; %----------------------------------------------------------------------------%
;             ;define base cell dimension in nanometers
(define-param cell                    1)

(define       a                       (* cell 1e-9)) ; 1nm
(define       c                       299762458.0)   ; speed of light (m/s)

;             ;basic Meep unit of time (sec)
;             ;i.e. time for light to move one "unit" a in a vacuum
(define       T_meep                  (/ a c))

;             ;basic Meep unit of frequency (rad / sec)
(define       omega_meep              (* (* 2 pi) (/ c a)))




; %----------------------------------------------------------------------------%
; |                                                                            |
; |      2. DEFINE WAVELENGTH AND FREQUENCY                                    |
; |                                                                            |
; %----------------------------------------------------------------------------%
;(define-param wavelen                 750) ; Wavelength in nanometers
;(define-param wavelen                 667) ; Wavelength in nanometers
(define-param wavelen                 706) ; Wavelength in nanometers
(define       wavelength-c            (/ wavelen  cell)) ;Wavelength in cells

(define-param gaussian_width          500) ; Full gaussian width
(define       gaussian_width-c        (/ gaussian_width cell))

(define       period                  wavelength-c)

; If we want a width of 500nm (500nm-1000nm)
;(define-param flux_width              1000) ; just a number that gives the
;                                            ; correct frequency width
; If we want a width of 700nm (500nm-1200nm)
(define-param flux_width              857 ) ; just a number that gives the
                                            ; correct frequency width
(define       flux_width-c            (/ flux_width cell))

;             ;frequency of our wavelength in Meep units of inverse time
(define       freq                    (/ 1 wavelength-c))
(define       df                      (/ 1 gaussian_width-c))
(define       dflux                   (/ 1 flux_width-c))
;(define-param nfreq                   100)
(define-param nfreq                   500)



; %----------------------------------------------------------------------------%
; |                                                                            |
; |      3. DEFINE SIMULATION OUTPUT CONTROLS                                  |
; |                                                                            |
; %----------------------------------------------------------------------------%
;(define-param largevol?               false)
(define-param objects?                true)
;(define-param output?                 true)
;(define-param output-to-screen?       true)

;(define-param number-of-run-periods   12)
;(define-param output-after-period     11)
;(define-param points-per-period       50)


; %----------------------------------------------------------------------------%
; |                                                                            |
; |      4. DEFINE PARAMETERS FOR COMPUTATIONAL CELL                           |
; |                                                                            |
; %----------------------------------------------------------------------------%
; |  NOTE:                                                                     |
; |   Make sure that SX, SY, and SZ are evenly divisible by CELL               |
; |   DPML is given in number of grid cells.                                   |
; |   ALL other parameters have dimensions of nanometers.                      |
; %----------------------------------------------------------------------------%
;(define-param pitch                   618)
;(define-param pitch                   615)
(define-param pitch                   700)
(define-param sx                      pitch) ;Make sure evenly divisible by cell
;(define-param sy                      900)
(define-param sy                      2000)
(define-param sz                      pitch)

(define-param res                     1)
(define-param dpml                    12) ;try a fixed num (/ wavelen 8))
(define-param xpml-c                  0)
(define-param ypml-c                  dpml)
(define-param zpml-c                  0)
(define flux-tag
  (string-append "refl-flux-pitch-" (number->string pitch)))

(define       sx-c                    (/ sx                cell))
(define       sy-c                    (/ sy                cell))
(define       sz-c                    (/ sz                cell))
(define       slx-c                   (+ sx-c (* 2 xpml-c)))
(define       sly-c                   (+ sy-c (* 2 ypml-c)))
(define       slz-c                   (+ sz-c (* 2 zpml-c)))



; %----------------------------------------------------------------------------%
; |                                                                            |
; |      5. DEFINE PARAMETERS FOR GEOMETRIC MODEL                              |
; |                                                                            |
; %----------------------------------------------------------------------------%
; |                                                                            |
; |  NOTE2:                                                                    |
; |    The origin for the simulation coordinate system is located at           |
; |    the center of the computational cell.                                   |
; |                                                                            |
; %----------------------------------------------------------------------------%


;-------------------------------------------------------------------------------
;  A. First, define the parameters in units of nanometers.
;-------------------------------------------------------------------------------
(define-param offset_x                 0)
(define-param offset_y                 0)
(define-param gold_depth               100)
(define-param pmma_depth               100)
(define-param cylinder_length          70)
(define-param cylinder_radius          157.5)

;-------------------------------------------------------------------------------
;  B. Next, express the parameters in units of cells.
;-------------------------------------------------------------------------------
(define       offset_x-c               (/ offset_x         cell))
(define       offset_y-c               (/ offset_y         cell))
(define       gold_depth-c             (/ gold_depth       cell))
(define       pmma_depth-c             (/ pmma_depth       cell))
(define       cylinder_length-c        (/ cylinder_length  cell))
(define       cylinder_radius-c        (/ cylinder_radius  cell))

;-------------------------------------------------------------------------------
;  C. Finally, calculate quantities needed for the geometric models.
;-------------------------------------------------------------------------------
(define       pmma_block_size_x-c      sx-c)
(define       pmma_block_size_y-c      (+ pmma_depth-c gold_depth-c))
(define       pmma_block_size_z-c      sz-c)
(define       pmma_block_center_y-c    (+ (*  0.5 sy-c)
                                          (* -0.5 pmma_block_size_y-c)))

(define       gold_block_size_x-c      sx-c)
(define       gold_block_size_y-c      gold_depth-c)
(define       gold_block_size_z-c      sz-c)
(define       gold_block_center_y-c    (+ (*  0.5 sy-c)
                                          (* -0.5 gold_block_size_y-c)))

(define       cylinder_center_y-c      (+ (*  0.5 sy-c)
                                          (* -1.0 gold_block_size_y-c)
                                          (* -0.5 cylinder_length-c)))



; %----------------------------------------------------------------------------%
; |                                                                            |
; |      6. DEFINE PARAMETERS FOR MATERIAL MODELS                              |
; |                                                                            |
; %----------------------------------------------------------------------------%
(define-param epsilon-real            1)
(define-param epsilon-imag            0)


;-----Material using conductivity-----
; See http://ab-initio.mit.edu/wiki/index.php/Materials_in_Meep
;-------------------------------------
;(define material-with-conductivity
;    (make medium
;        (epsilon            epsilon-real)
;        (D-conductivity     (/ (* 2 pi freq epsilon-imag) epsilon-real)) ))

;-----Silver-1term-----
;  Model good from 1.14eV (1088nm) to 3.0eV (413nm)
;----------------------
;ag_epsilon_inf=3.9845e+00
;ag_omega_d    =1.3985e+16
;ag_gamma_d    =1.4348e+13
;ag_d_epsilon_1=2.2048e-01
;ag_omega_1    =3.4315e+15
;ag_gamma_1    =1.1473e+15

(define       ag_epsilon_inf          3.9845e+00)

(define       ag_omega_d              (/ 1.3985e+16 omega_meep))
(define       ag_omega_1              (/ 3.4315e+15 omega_meep))

(define       ag_gamma_d              (/ 1.4348e+13 omega_meep))
(define       ag_gamma_1              (/ (* 1.1473e+15 2) omega_meep))

(define       ag_sigma_1              2.2048e-01)

(define silver-1term
    (make dielectric
        (epsilon ag_epsilon_inf)
        (E-polarizations
            ;Drude term
            (make polarizability
                (omega 1e-20)
                (gamma ag_gamma_d)
                (sigma (* 1e20 ag_omega_d 1e20 ag_omega_d ))
            )
            ;Lorentz 1st term
            (make polarizability
                (omega ag_omega_1)
                (gamma ag_gamma_1)
                (sigma ag_sigma_1)
            )
        )
    )
)


;-----Gold-3term-----
;  Model good from 1.0eV (1240nm) to 3.0eV (413nm)
;--------------------
;au_epsilon_inf=4.8929e+00
;au_omega_d    =1.2944e+16
;au_gamma_d    =1.0003e+09
;au_d_epsilon_1=4.7282e+00
;au_omega_1    =1.3617e+15
;au_gamma_1    =4.7356e+14
;au_d_epsilon_2=7.2996e-01
;au_omega_2    =4.1636e+15
;au_gamma_2    =4.4931e+14
;au_d_epsilon_3=1.5103e+00
;au_omega_3    =5.0753e+15
;au_gamma_3    =5.8469e+14
(define       au_epsilon_inf          4.8929e+00)

; plasma freq in Meep units
(define       au_omega_d              (/ 1.2944e+16 omega_meep))
(define       au_omega_1              (/ 1.3617e+15 omega_meep))
(define       au_omega_2              (/ 4.1636e+15 omega_meep))
(define       au_omega_3              (/ 5.0753e+15 omega_meep))

; collision freq in Meep units
(define       au_gamma_d              (/ 1.0003e+09 omega_meep))
(define       au_gamma_1              (/ (* 4.7356e+14 2) omega_meep))
(define       au_gamma_2              (/ (* 4.4931e+14 2) omega_meep))
(define       au_gamma_3              (/ (* 5.8469e+14 2) omega_meep))

(define       au_sigma_1              4.7282e+00)
(define       au_sigma_2              7.2996e-01)
(define       au_sigma_3              1.5103e+00)

(define gold-3term
    (make dielectric
        (epsilon au_epsilon_inf)
        (E-polarizations
            ;Drude term
            (make polarizability
                (omega 1e-20)
                (gamma au_gamma_d)
                (sigma (* (* 1e20 au_omega_d) (* 1e20 au_omega_d) ))
            )
            ;Lorentz 1st term
            (make polarizability
                (omega au_omega_1)
                (gamma au_gamma_1)
                (sigma au_sigma_1)
            )
            ;Lorentz 2nd term
            (make polarizability
                (omega au_omega_2)
                (gamma au_gamma_2)
                (sigma au_sigma_2)
            )
            ;Lorentz 3rd term
            (make polarizability
                (omega au_omega_3)
                (gamma au_gamma_3)
                (sigma au_sigma_3)
            )
        )
    )
)

;Generic dielectric
(define       generic                 (make dielectric (epsilon 25)));
(define       glass                   (make dielectric (epsilon 2.25)));
;(define       pmma                    (make dielectric (epsilon 1.48)));
;(define       pmma                    (make dielectric (epsilon 2.20)));
(define       pmma                    (make dielectric (index 1.48236)));
;-------------------------------------------------------------------------------



; %----------------------------------------------------------------------------%
; |                                                                            |
; |      7. DEFINE PARAMETERS FOR PLANE-WAVE EXCITATION                        |
; |                                                                            |
; %----------------------------------------------------------------------------%
; | Note:  Periodicity is defined by the k-vector of the excitation, not by a  |
; | periodicity on the y-dimension.  This produces a excitation field which    |
; | appears to be truly infinite.  See the experiment-geom3d ctl files         |
; %----------------------------------------------------------------------------%
(define-param theta                   0)
(define-param theta-x                 theta)
(define-param theta-y                 theta)

;; k-vector in x direction (left to right), y-polarized
;(define       kx                      (* freq (cos (deg->rad theta-x))))
;(define       ky                      (* freq (sin (deg->rad theta-x))))
;(define       kz                      0)

; k-vector in -y direction (bottom to top), x-polarized
(define       kx                      (* freq (sin (deg->rad theta-y))))
(define       ky                      (* freq (cos (deg->rad theta-y))))
(define       kz                      0)

(define       k                       (vector3 kx ky kz))
(define       (my-amp-func p)         (exp (* 0+2i pi (vector3-dot k p))))




; %----------------------------------------------------------------------------%
; |                                                                            |
; |      8. DEFINE BOUNDARY CONDITIONS                                         |
; |                                                                            |
; %----------------------------------------------------------------------------%
(set!         k-point                 k)
(set!         ensure-periodicity      true)
;(set!         ensure-periodicity      false)




; %----------------------------------------------------------------------------%
; |                                                                            |
; |      9. DEFINE COMPUTATIONAL CELL                                          |
; |                                                                            |
; %----------------------------------------------------------------------------%
(set! geometry-lattice (make lattice (size slx-c sly-c slz-c)))

;(set! default-material pmma)


; %----------------------------------------------------------------------------%
; |                                                                            |
; |      10. DEFINE GEOMETRIC MODEL                                            |
; |                                                                            |
; %----------------------------------------------------------------------------%

(if objects?
    (set! geometry
        (list
;-------------------------------------------------------------------------------
; Materials-
;   gold-3term     --Model good from 1.00eV (1240nm) to 3.00eV (413nm)
;   silver-2term   --Model good from 1.14eV (1088nm) to 4.12eV (301nm)
;   silver-1term   --Model good from 1.14eV (1088nm) to 3.00eV (413nm)
;   glass          --epsilon = 2.25
;-------------------------------------------------------------------------------
            (make block
                (center (vector3 0 pmma_block_center_y-c 0))
                (size   (vector3 pmma_block_size_x-c 
                                 pmma_block_size_y-c
                                 pmma_block_size_z-c))
                (material pmma)
            )
            (make block
                (center (vector3 0 gold_block_center_y-c 0))
                (size   (vector3 gold_block_size_x-c 
                                 gold_block_size_y-c
                                 gold_block_size_z-c))
                (material gold-3term)
            )
            (make cylinder
                (center (vector3 0 cylinder_center_y-c 0))
                (axis   (vector3 0 1 0))
                (height cylinder_length-c)
                (radius cylinder_radius-c)
                (material gold-3term)
            )
        )
    )
)


; %----------------------------------------------------------------------------%
; |                                                                            |
; |      11. DEFINE PML LAYERS                                                 |
; |                                                                            |
; %----------------------------------------------------------------------------%
; y-direction propagation
(set! pml-layers
    (list
        (make pml
            (thickness ypml-c)
            (direction Y)
        )
    )
)



; %----------------------------------------------------------------------------%
; |                                                                            |
; |      12. DEFINE EXCITATION SOURCES                                         |
; |                                                                            |
; %----------------------------------------------------------------------------%
;Top of cell (y-axis propagation)
(set! sources
    (list
        (make source
            (src (make gaussian-src 
;                     (cutoff 6.0)
                     (frequency freq)
                     (fwidth (/ 1 gaussian_width-c))))
                     ;(wavelength wavelength-c)
                     ;(width gaussian_width-c)))
            (component Hz)
            (center (vector3 0 (+ 1 (* -0.5 sy-c)) 0))
            (size   (vector3 sx-c 0 sz-c))
            (amp-func my-amp-func)
        )
    )
)


; %----------------------------------------------------------------------------%
; |                                                                            |
; |      13. DEFINE RESOLUTION                                                 |
; |                                                                            |
; %----------------------------------------------------------------------------%
(set! resolution res)


; %----------------------------------------------------------------------------%
; |                                                                            |
; |      14. DEFINE COMPUTATIONAL SYMMETRIES                                   |
; |                                                                            |
; %----------------------------------------------------------------------------%
;(set! symmetries
;    (list
;        (make mirror-sym
;            (direction Z)
;            (phase 1)
;        )
;    )
;)




;(define refl ;reflected flux, located near the source
(define flux1 ;reflected flux, located near the source
    (add-flux freq dflux nfreq
        (make flux-region
            (size   (vector3 sx-c 0 sz-c))
            ;(center (vector3 0 (+ 2 (* -0.5 sy-c)) 0))
            (center (vector3 0 (+ 5 (* -0.5 sy-c)) 0))
        )
    )
)

;(define flux2
;    (add-flux freq dflux nfreq
;        (make flux-region
;            (size   (vector3 sx-c 0 sz-c))
;            (center (vector3 0 (* -0.25 sy-c) 0))
;        )
;    )
;)

;(define flux3
;    (add-flux freq dflux nfreq
;        (make flux-region
;            (size   (vector3 sx-c 0 sz-c))
;            (center (vector3 0 0 0))
;        )
;    )
;)

;(define flux4
;    (add-flux freq dflux nfreq
;        (make flux-region
;            (size   (vector3 sx-c 0 sz-c))
;            (center (vector3 0 (*  0.25 sy-c) 0))
;        )
;    )
;)

;(define trans ;transmitted flux, located below the substrate
;(define flux5 ;transmitted flux, located below the substrate
;    (add-flux freq dflux nfreq
;        (make flux-region
;            (size   (vector3 sx-c 0 sz-c))
;            ;(center (vector3 0 (+ -2 (* 0.5 sy-c)) 0))
;            (center (vector3 0 (+ -5 (* 0.5 sy-c)) 0))
;        )
;    )
;)



(set! eps-averaging?            true)
(set! output-single-precision?  true)
(set! progress-interval         30)
(use-output-directory
  (string-append "meep-out-pitch-" (number->string pitch)))

(if objects? (load-minus-flux flux-tag flux1))
(run-sources+
;    (stop-when-fields-decayed 50 Ex (vector3 0 (+ 2 (* -0.5 sy-c)) 0) 1e-6))
    (stop-when-fields-decayed 50 Ex (vector3 0 0 0) 1e-6))
(if (not objects?) (save-flux flux-tag flux1))

;(display-fluxes flux1 flux2 flux3 flux4 flux5)
(display-fluxes flux1)





