.. _ch:AD:

*******************
Ambipolar Diffusion
*******************

Introduction
============

Ambipolar diffusion is modelled with two fluids: neutrals (hydrodynamic) and
ions (magnetohydrodynamic), coupled by the collisional drag
:math:`\alpha\rho_n\rho_i(\boldsymbol{u}_i - \boldsymbol{u}_n)` and the
associated frictional heating and thermal exchange (:ref:`ch:hydro`,
two-fluid equations). In the code the two fluids are two ``grid`` objects,
``gn`` and ``gi``, with their own ``gridID`` s (a *pair* of case numbers,
e.g. 609/616), evolved together by one of the two-fluid drivers. The
coupling coefficient is ``alpha_ad`` in code units
(:math:`{\rm pc^2\,km\,s^{-1}}\,M_\odot^{-1}`; :math:`3.7\times10^{13}`
cm\ :sup:`3` g\ :sup:`-1` s\ :sup:`-1` ≈ :math:`7.7\times10^{4}` code units)
and ``mu_ad`` the molecular weight of each species (2.3 / 29 for a molecular
cloud), both set with ``setADparams`` on both grids after ``enableAD``.

Two ways to integrate the drag
==============================

The drag is stiff whenever the coupling time
:math:`\tau = 1/(\alpha\rho_i)` is short compared with the time step. Two
families of integrators are available, chosen at run time with
``SCORPIO_AD_SCHEME`` (:ref:`ch:problem_file`):

.. list-table::
   :header-rows: 1
   :widths: 24 22 54

   * - value
     - driver
     - scheme
   * - ``split`` (default)
     - ``rk2AD_3D_HSHSMD`` (``rk2.f03``); ``rk2AD_2D_HSHSMD`` in 2D
     - operator splitting: in every RK2 stage the hyperbolic sweeps of both fluids, then the
       drag update ``evolveAD3D_MD`` (``evolveAmbipolarDiffusion.f03``) with the TR-BDF2 scheme of
       Tilley, Balsara & Meyer (2012) — the method of the paper
   * - ``imex``
     - ``rk2AD_3D_IMEX`` (``evolveAD_imex.f03``)
     - IMEX-SSP2(2,2,2) of Pareschi & Russo (2005), :math:`\gamma = 1 - 1/\sqrt{2}`: the drag is
       solved implicitly *inside* each RK stage (DIRK), the fluxes stay explicit
   * - ``imex322`` / ``athenak``
     - ``rk2AD_3D_IMEX322``
     - IMEX-SSP2(3,2,2) (Pareschi & Russo 2005, Table III) — the scheme of AthenaK's ion–neutral
       module; L-stable and stiffly accurate
   * - ``imexpp`` / ``imex2+`` / ``krapp``
     - ``rk2AD_3D_IMEXPP``
     - IMEX(4,3,2) of Krapp et al. (2024) with the monotone :math:`\gamma_+ = 1 + 1/\sqrt{2}`
       root (AthenaK "imex2+"); no ringing in stiff relaxation; CFL ≤ 0.35 (its first explicit
       stage advances :math:`1.707\,\Delta t`)

The explicit part of the IMEX schemes is exactly Scorpio's Heun/SSP-RK2, so
the hyperbolic machinery (PLM/PPM, HLLD hardening, CT, dual energy, FOFC) is
identical in all four drivers; only the placement of the drag differs. What
the choice changes is the ion–neutral *drift* :math:`\boldsymbol{u}_i -
\boldsymbol{u}_n`, i.e. the residual of the stiff balance and the science
observable: the split scheme damps it spuriously at :math:`\Delta t/\tau
\gtrsim 6`, the IMEX schemes do not. The IMEX drivers exist in 3D; the 2D
driver refuses them rather than silently running the split scheme, and an
unrecognised value aborts.

Call path
=========

.. code-block:: text

   ADMHD3D time loop                                          testSuiteMPI.f03
   ├─ gn%griddt(); gi%griddt();  dt = min → gn%dt = gi%dt = dt
   ├─ SCORPIO_AD_SCHEME → rk2AD_3D_HSHSMD | rk2AD_3D_IMEX | rk2AD_3D_IMEX322 | rk2AD_3D_IMEXPP
   └─ ghost exchange, boundaries and writeGrid on both grids

   rk2AD_3D_HSHSMD(gn, qn, qn1, qn2, gi, qi, qi1, qi2)          rk2.f03
   ├─ rieSolvern => solverAdi3D | solverIso3D ;  rieSolveri => solverAdiMHD3D | solverIsoMHD3D
   ├─ do (failsafe loop)
   │   stage 1: rieSolvern(gn, qn, …, dd=1..3);  rieSolveri(gi, qi, …, dd=1..3)
   │            evolveAD3D_MD(gn, gi, qn, qi, qn1, qi1)       drag, frictional heating, thermal exchange (TR-BDF2)
   │            post-drag dual-energy re-sync of the ion energy
   │            exchangeBdryMPI / setBoundary on qn1, qi1
   │   stage 2: the same from (qn1, qi1) into (qn2, qi2)
   │   negative-state check on both fluids → FOFC / Δt halving / abort  (as in rk2_3D)
   │  enddo
   └─ Heun average on both grids, ghost fill, t += dt

   rk2AD_3D_IMEXPP(…)                                          evolveAD_imex.f03
   ├─ per stage: explicit sweeps (same rieSolver pointers)
   │            adImexDragStage(gn, gi, …, γΔt)   per-cell exact 2×2 solve of u = u_in + γΔt R(u)
   │                                             momentum exchange + frictional heating + thermal exchange
   │            ad_imex_resync(gi, …)             dual-energy re-sync of the ion energy
   └─ closed-form Butcher assembly of u^{n+1}; total pair momentum and energy conserved to round-off

``evolveAD1D/2D/3D`` are the older per-dimension drag routines;
``evolveAD3D_hdt`` a half-step variant; ``calcDT3D_AD`` the driving for the
two-fluid setup (:ref:`ch:turbulence`).

Test cases
==========

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - ``gridID`` pair
     - case
   * - 22/23, 122/123, 52/53
     - C-shock in 1D, 2D, 3D (the Tilley, Balsara & Meyer 2012 setup; 52/53 can start from the semi-analytic steady profile and includes the AthenaK comparison hooks ``SCORPIO_CSK_*``)
   * - 26/27
     - Wardle instability
   * - 28/29
     - 2D core collapse with AD
   * - 54/55
     - two-fluid spherical blast (the original low-β negative-pressure case; 0 negatives with dual energy)
   * - 500/501, 600/601
     - Kelvin–Helmholtz and Richtmyer–Meshkov instabilities with AD
   * - 609/616
     - 3D driven ion–neutral turbulence (``SCORPIO_TURB_*``, ``SCORPIO_COUPLE_T``)
   * - 620/621
     - Alfvén-wave damping against the Kulsrud & Pearce (1969) dispersion relation (``SCORPIO_ADW_*``)

The C-shock and Alfvén-damping harnesses are the reference tests for the
coupling schemes; the measured properties of each scheme are summarised in
:ref:`ch:methods` and, in full, in ``CHANGES.md`` §8 and §10.
