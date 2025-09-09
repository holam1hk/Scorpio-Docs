Welcome to Scorpio's Documentation!
=================================================
**Scorpio** is a grid-based finite-volume code developed by the *Star Formation Group (SFG)* at the Chinese University of Hong Kong (CUHK). It is designed for astrophysical simulations, including star formation and planet formation, supporting Hydrodynamics (HD), Magnetohydrodynamics (MHD), and Ambipolar Diffusion MHD (ADMHD).

The code solves fluid conservation equations using Godunov-type methods, employing a second-order accurate Total Variation Diminishing (TVD) scheme. Key physics modules include:

*   **Flux Computation:** Utilizes a combination of the Harten-Lax-van Leer-Contact (HLLC) and Harten-Lax-van Leer-Discontinuities (HLLD) approximate Riemann solvers. A Piecewise Linear Method (PLM) for cell interpolation and various limiters are used to ensure the TVD property.
*   **Magnetic Fields:** The divergence-free constraint of the magnetic field is maintained through the Constrained Transport (CT) method.
*   **Self-Gravity:** Self-gravitating fluids are modeled using a Fast Fourier Transform (FFT)-based Poisson solver.
*   **Turbulence:** Turbulence is driven with a power spectrum using FFT.
*   **Ambipolar Diffusion:** The two-fluid model for ion-neutral drift is implemented with the second-order accurate TR-BDF2 (Trapezoidal Rule with second-order Backward Differentiation Formula) scheme, following the methodology of Tilley et al. (2012).


More details can be found in `our paper on Scorpio <https://www.phy.cuhk.edu.hk/sfg/publications/Scorpio_RASTI.pdf>`_.

Contents
--------

.. note::

	This is still under development.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   getting_started
   hydro
   problem_file
   code_structure 
   selfgravity
   turbulence
   AD
   time_integration
   riemann
   limiter
   parallelization






