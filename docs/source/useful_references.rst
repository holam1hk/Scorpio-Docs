.. _ch:useful_references:

*****************
Useful References
*****************

MHD and astrophysical simulation codes
======================================

- `RAMSES <https://bitbucket.org/rteyssie/ramses/>`_
  - Teyssier, R. (2002). Cosmological hydrodynamics with adaptive mesh refinement: a new high resolution code called RAMSES. Astronomy & Astrophysics, 385, 337–364. DOI: https://doi.org/10.1051/0004-6361:20011817

- `PLUTO <http://plutocode.ph.unito.it/>`_
  - Mignone, A., Bodo, G., Massaglia, S., Matsakos, T., Tesileanu, O., Zanni, C., & Ferrari, A. (2007). PLUTO: A Numerical Code for Computational Astrophysics. The Astrophysical Journal Supplement Series, 170(1), 228–242. DOI: https://doi.org/10.1086/513316

- `ENZO <https://enzo-project.org/>`_
  - Bryan, G. L., et al. (2014). ENZO: An Adaptive Mesh Refinement Code for Astrophysics. The Astrophysical Journal Supplement Series, 211(2), 19. DOI: https://doi.org/10.1088/0067-0049/211/2/19

- `FLASH <https://flash.rochester.edu/site/flashcode/>`_
  - Fryxell, B., et al. (2000). FLASH: An Adaptive Mesh Hydrodynamics Code for Modeling Astrophysical Thermonuclear Flashes. The Astrophysical Journal Supplement Series, 131(1), 273–334. DOI: https://doi.org/10.1086/317361

- `Athena <https://github.com/PrincetonUniversity/athena-public-version>`_
  - Stone, J. M., Gardiner, T. A., Teuben, P., Hawley, J. F., & Simon, J. B. (2008). Athena: A New Code for Astrophysical MHD. The Astrophysical Journal Supplement Series, 178(1), 137–177. DOI: https://doi.org/10.1086/588755

- `Athena++ <https://github.com/PrincetonUniversity/athena>`_
  - Stone, J. M., Tomida, K., White, C. J., & Felker, K. G. (2020). The Athena++ Adaptive-Framework: Design and Magnetohydrodynamic Solvers. The Astrophysical Journal Supplement Series, 249(1), 4. DOI: https://doi.org/10.3847/1538-4365/ab929b

- `AthenaK <https://github.com/IAS-Astrophysics/athenak>`_
  - Stone, J. M., Mullen, P. D., Fielding, D., Grete, P., Guo, M., Kempski, P., Most, E. R., White, C. J., & Wong, G. N. (2024). AthenaK: A Performance-Portable Version of the Athena++ AMR Framework. arXiv e-prints. DOI: https://doi.org/10.48550/arXiv.2409.16053 (See `ADS abstract <https://ui.adsabs.harvard.edu/abs/2024arXiv240916053S/abstract>`_)

- `ZEUS <https://www2.astro.psu.edu/xray/astro511/Zeus/>`_
  - Stone, J. M., & Norman, M. L. (1992). ZEUS-2D: A radiation magnetohydrodynamics code for astrophysical flows in two space dimensions. I. The hydrodynamic algorithms and tests. The Astrophysical Journal Supplement Series, 80, 753–790. DOI: https://doi.org/10.1086/191680

- `GIZMO <http://www.tapir.caltech.edu/~phopkins/Site/GIZMO.html>`_
  - Hopkins, P. F. (2015). A new class of accurate, mesh-free hydrodynamic simulation methods. Monthly Notices of the Royal Astronomical Society, 450(1), 53–110. DOI: https://doi.org/10.1093/mnras/stv195

- `PHANTOM <https://phantomsph.readthedocs.io>`_
  - Price, D. J., et al. (2018). Phantom: A Smoothed Particle Hydrodynamics and Magnetohydrodynamics Code for Astrophysics. Publications of the Astronomical Society of Australia, 35, e031. DOI: https://doi.org/10.1017/pasa.2018.25

- `AREPO <https://arepo-code.org>`_
  - Springel, V. (2010). E pur si muove: Galilean-invariant cosmological hydrodynamical simulations on a moving mesh. Monthly Notices of the Royal Astronomical Society, 401(2), 791–851. DOI: https://doi.org/10.1111/j.1365-2966.2009.15715.x


Methods used in Scorpio
=======================

Scorpio itself

- Cheng, H. L., Wang, H.-H., Zeng, W.-G., Cao, Z., Luk, S. S., Tsang, M. H., & Li, H.-b. (2025). Scorpio: a two-fluid code for ambipolar diffusion simulations informed by recent observational constraints. RAS Techniques and Instruments, 4, 1–17. DOI: https://doi.org/10.1093/rasti/rzaf054 (`PDF <https://www.phy.cuhk.edu.hk/sfg/publications/Scorpio_RASTI.pdf>`_)

Riemann solvers and reconstruction

- Toro, E. F., Spruce, M., & Speares, W. (1994). Restoration of the contact surface in the HLL-Riemann solver. Shock Waves, 4, 25–34. — HLLC
- Miyoshi, T., & Kusano, K. (2005). A multi-state HLL approximate Riemann solver for ideal magnetohydrodynamics. Journal of Computational Physics, 208, 315–344. — HLLD, and the degenerate-state prescription used by the hardened kernel
- Mignone, A. (2007). A simple and accurate Riemann solver for isothermal MHD. Journal of Computational Physics, 225, 1427–1441. — isothermal HLLD
- Colella, P., & Woodward, P. R. (1984). The piecewise parabolic method (PPM) for gas-dynamical simulations. Journal of Computational Physics, 54, 174–201. — PPM
- Colella, P., & Sekora, M. D. (2008). A limiter for PPM that preserves accuracy at smooth extrema. Journal of Computational Physics, 227, 7069–7076. — tested, not adopted
- Zhang, X., & Shu, C.-W. (2010). On positivity-preserving high order discontinuous Galerkin schemes for compressible Euler equations on rectangular meshes. Journal of Computational Physics, 229, 8918–8934. — the θ positivity guard
- Skinner, M. A., & Ostriker, E. C. (2010). The Athena astrophysical MHD code in cylindrical geometry. The Astrophysical Journal Supplement Series, 188, 290–311. — limiters in cylindrical coordinates

Constrained transport and positivity

- Balsara, D. S., & Spicer, D. S. (1999). A staggered mesh algorithm using high order Godunov fluxes to ensure solenoidal magnetic fields in magnetohydrodynamic simulations. Journal of Computational Physics, 149, 270–292. — centred CT EMF
- Gardiner, T. A., & Stone, J. M. (2005). An unsplit Godunov method for ideal MHD via constrained transport. Journal of Computational Physics, 205, 509–539. — upwinded corner EMF (default)
- Gardiner, T. A., & Stone, J. M. (2008). An unsplit Godunov method for ideal MHD via constrained transport in three dimensions. Journal of Computational Physics, 227, 4123–4141.
- Londrillo, P., & Del Zanna, L. (2004). On the divergence-free condition in Godunov-type schemes for ideal magnetohydrodynamics: the upwind constrained transport method. Journal of Computational Physics, 195, 17–48.
- Mignone, A., & Del Zanna, L. (2021). Systematic construction of upwind constrained transport schemes for MHD. Journal of Computational Physics, 424, 109748.
- Balsara, D. S. (2001). Divergence-free adaptive mesh refinement for magnetohydrodynamics. Journal of Computational Physics, 174, 614–648. — divergence-preserving prolongation on the AMR mesh
- Stone, J. M., Tomida, K., White, C. J., & Felker, K. G. (2020). The Athena++ adaptive mesh refinement framework: design and magnetohydrodynamic solvers. The Astrophysical Journal Supplement Series, 249, 4. — first-order flux correction (FOFC)

Self-gravity

- James, R. A. (1977). The solution of Poisson's equation for isolated source distributions. Journal of Computational Physics, 25, 71–93. — screening-charge isolated boundary of the multigrid solver
- Ricker, P. M. (2008). A direct multigrid Poisson solver for oct-tree adaptive meshes. The Astrophysical Journal Supplement Series, 176, 293–300.
- Moon, S., Kim, W.-T., & Ostriker, E. C. (2019). A fast Poisson solver of second-order accuracy for isolated systems in three-dimensional Cartesian and cylindrical coordinates. The Astrophysical Journal Supplement Series, 241, 24.
- Tomida, K., & Stone, J. M. (2023). The Athena++ adaptive mesh refinement framework: multigrid solvers for self-gravity. The Astrophysical Journal Supplement Series, 266, 7.
- Chandrasekhar, S. (1969). Ellipsoidal Figures of Equilibrium. Yale University Press. — Maclaurin-spheroid gate
- Truelove, J. K., et al. (1997). The Jeans condition: a new constraint on spatial resolution in simulations of isothermal self-gravitational hydrodynamics. The Astrophysical Journal, 489, L179–L183. — Truelove stop / Jeans refinement

Adaptive mesh refinement

- Löhner, R. (1987). An adaptive finite element scheme for transient problems in CFD. Computer Methods in Applied Mechanics and Engineering, 61, 323–338. — refinement indicator
- Fryxell, B., et al. (2000). FLASH (see above). — block-structured octree design
- Berger, M. J., & Colella, P. (1989). Local adaptive mesh refinement for shock hydrodynamics. Journal of Computational Physics, 82, 64–84. — refluxing

Ion–neutral coupling

- Tilley, D. A., Balsara, D. S., & Meyer, C. (2012). A numerical scheme and benchmark tests for non-isothermal two-fluid ambipolar diffusion. New Astronomy, 17, 368–376. — TR-BDF2 drag, C-shock benchmark
- Pareschi, L., & Russo, G. (2005). Implicit-explicit Runge–Kutta schemes and applications to hyperbolic systems with relaxation. Journal of Scientific Computing, 25, 129–155. — IMEX-SSP2(2,2,2), (3,2,2)
- Krapp, L., Garrido-Deutelmoser, J., Benítez-Llambay, P., & Kratter, K. M. (2024). A fast second-order solver for stiff multifluid dust and gas hydrodynamics. The Astrophysical Journal Supplement Series (arXiv:2310.04435). — the monotone IMEX(4,3,2) scheme (AthenaK "imex2+")
- Kulsrud, R., & Pearce, W. P. (1969). The effect of wave-particle interactions on the propagation of cosmic rays. The Astrophysical Journal, 156, 445. — two-fluid Alfvén-wave dispersion relation (damping gate)

Turbulence driving

- Federrath, C., Roman-Duval, J., Klessen, R. S., Schmidt, W., & Mac Low, M.-M. (2010). Comparing the statistics of interstellar turbulence in simulations and observations: solenoidal versus compressive turbulence forcing. Astronomy & Astrophysics, 512, A81. — the projection operator (ζ)
- Otto, F., Ji, W., & Li, H.-b. (2017). Velocity anisotropy in self-gravitating molecular clouds. I. Simulation. The Astrophysical Journal, 836, 95. — the ``expo`` driving spectrum



.. BibTeX
.. ======

.. .. code-block:: bibtex

..    @article{Teyssier2002RAMSES,
..      author  = {Teyssier, R.},
..      title   = {Cosmological hydrodynamics with adaptive mesh refinement: a new high resolution code called RAMSES},
..      journal = {Astronomy & Astrophysics},
..      volume  = {385},
..      pages   = {337-364},
..      year    = {2002},
..      doi     = {10.1051/0004-6361:20011817}
..    }

..    @article{Mignone2007PLUTO,
..      author  = {Mignone, A. and Bodo, G. and Massaglia, S. and Matsakos, T. and Tesileanu, O. and Zanni, C. and Ferrari, A.},
..      title   = {PLUTO: A Numerical Code for Computational Astrophysics},
..      journal = {The Astrophysical Journal Supplement Series},
..      volume  = {170},
..      number  = {1},
..      pages   = {228-242},
..      year    = {2007},
..      doi     = {10.1086/513316}
..    }

..    @article{Bryan2014ENZO,
..      author  = {Bryan, G. L. and others},
..      title   = {ENZO: An Adaptive Mesh Refinement Code for Astrophysics},
..      journal = {The Astrophysical Journal Supplement Series},
..      volume  = {211},
..      number  = {2},
..      pages   = {19},
..      year    = {2014},
..      doi     = {10.1088/0067-0049/211/2/19}
..    }

..    @article{Fryxell2000FLASH,
..      author  = {Fryxell, B. and others},
..      title   = {FLASH: An Adaptive Mesh Hydrodynamics Code for Modeling Astrophysical Thermonuclear Flashes},
..      journal = {The Astrophysical Journal Supplement Series},
..      volume  = {131},
..      number  = {1},
..      pages   = {273-334},
..      year    = {2000},
..      doi     = {10.1086/317361}
..    }

..    @article{Stone2020AthenaPP,
..      author  = {Stone, J. M. and Tomida, K. and White, C. J. and Felker, K. G.},
..      title   = {The Athena++ Adaptive-Framework: Design and Magnetohydrodynamic Solvers},
..      journal = {The Astrophysical Journal Supplement Series},
..      volume  = {249},
..      number  = {1},
..      pages   = {4},
..      year    = {2020},
..      doi     = {10.3847/1538-4365/ab929b}
..    }

..    @article{Stone1992ZEUS2D,
..      author  = {Stone, J. M. and Norman, M. L.},
..      title   = {ZEUS-2D: A radiation magnetohydrodynamics code for astrophysical flows in two space dimensions. I. The hydrodynamic algorithms and tests},
..      journal = {The Astrophysical Journal Supplement Series},
..      volume  = {80},
..      pages   = {753-790},
..      year    = {1992},
..      doi     = {10.1086/191680}
..    }

..    @article{Hopkins2015GIZMO,
..      author  = {Hopkins, P. F.},
..      title   = {A new class of accurate, mesh-free hydrodynamic simulation methods},
..      journal = {Monthly Notices of the Royal Astronomical Society},
..      volume  = {450},
..      number  = {1},
..      pages   = {53-110},
..      year    = {2015},
..      doi     = {10.1093/mnras/stv195}
..    }

..    @article{Price2018PHANTOM,
..      author  = {Price, D. J. and Wurster, J. and Tricco, T. S. and Nixon, C. and Toupin, S. and Pettitt, A. R. and Chan, C. and Mentiplay, D. and Laibe, G. and Throssell, K. and Nealon, R. and Liptai, D. and Bate, M. R. and Pinte, C. and Lesur, G. and Forgan, D. and Ballabio, G. and Hutchison, M. and Commer{ 1}on, B. and Young, A. K. and Dobbs, C. L. and Worpel, H. and Hirsh, K. and Lodato, G. and Dipierro, G. and Ragusa, E. and Alexander, R. D. and Cuello, N. and Humphries, R. and Cleeves, L. I. and Gonzalez, J.-F. and Casassus, S. and Zurlo, A. and Price, M. A. and P 0E9rez, S.},
..      title   = {Phantom: A Smoothed Particle Hydrodynamics and Magnetohydrodynamics Code for Astrophysics},
..      journal = {Publications of the Astronomical Society of Australia},
..      volume  = {35},
..      pages   = {e031},
..      year    = {2018},
..      doi     = {10.1017/pasa.2018.25}
..    }

..    @article{Springel2010AREPO,
..      author  = {Springel, V.},
..      title   = {E pur si muove: Galilean-invariant cosmological hydrodynamical simulations on a moving mesh},
..      journal = {Monthly Notices of the Royal Astronomical Society},
..      volume  = {401},
..      number  = {2},
..      pages   = {791-851},
..      year    = {2010},
..      doi     = {10.1111/j.1365-2966.2009.15715.x}
..    }


