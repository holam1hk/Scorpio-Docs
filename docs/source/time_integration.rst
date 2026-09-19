.. _ch:time_integration:

****************
Time Integration
****************

Introduction
============

One time step of a uniform-grid run is two calls made by the problem driver,
``griddt`` (choose :math:`\Delta t`) and ``evolveGridRK2`` (advance the
state), plus an optional ``writeGrid``. This page follows those calls all the
way down to the flux kernels, so that every routine on the path is named with
the file it lives in. All names are the ones in the current source.

The full call path of a run
===========================

.. code-block:: text

   program main                                              main.f03
   └─ runProblemFromNamelist                                 problemRegistry.f03
      └─ dispatchProblemByGridID(gridID)                     problemRegistry.f03  (select case on gridID)
         └─ <problem driver>, e.g. cloud_20pc3_3DMHD         hinnyCloud.f03 / testSuiteMPI.f03 / …
            ├─ grid setup: setGridID, enableDrivingTurbulence, enableSelfgravitySolver,
            │  setTopologyMPI, setMesh, setVariable, setMPIWindows, setEoS, …, setTime,
            │  initSelfgravitySolver, initVariable, exchangeBdryMPI, setBoundary, writeGrid
            │                                                 (order and meaning: Problem File page)
            └─ time loop:  do while (g%t < g%tend)
                 ├─ g%griddt()                                gridModule.f03 → gridModule_dt.f03
                 ├─ [DT_mode = 1] g%calcDrivingTurbulence_MD  gridModule.f03 → gridModule_dtcalc.f03
                 ├─ g%evolveGridRK2()                         gridModule.f03 → gridModule_rk2.f03
                 └─ if (g%writeFlag) g%writeGrid()            gridModule.f03 → gridModule_io3d.f03

Step 1: the time step — ``griddt``
==================================

``griddt`` (``gridModule.f03``) dispatches on ``ndim`` to ``dt1D`` /
``dt2D`` / ``dt3D`` (``gridModule_dt.f03``). ``dt3D`` does, in order:

1. Loop over the interior cells and take the minimum of
   :math:`{\rm CFL}\,\min(\Delta x,\Delta y,\Delta z)/w` with the fastest
   signal speed :math:`w = |\boldsymbol{u}| + c`:

   .. list-table::
      :header-rows: 1
      :widths: 30 70

      * - ``eosType`` / ``solverType``
        - :math:`c`
      * - isothermal HD (1, 1–2)
        - :math:`c_s`
      * - isothermal MHD (1, 4–5)
        - :math:`c_f = \sqrt{\tfrac{1}{2}\left(c_s^2 + b^2 + \sqrt{(c_s^2+b^2)^2 - 4 c_s^2 b_{\min}^2/\rho}\right)}`, :math:`b^2 = |\boldsymbol{B}|^2/\rho`
      * - adiabatic HD (2, 2–3)
        - :math:`\sqrt{\Gamma P/\rho}`
      * - adiabatic MHD (2, 4–5)
        - :math:`c_f = \sqrt{\left(\Gamma P + b^2\rho + \sqrt{(\Gamma P + b^2\rho)^2 - 4\Gamma P\, b_{\min}^2}\right)/2\rho}`

   (:math:`b_{\min}` is the smallest of the three field components — a
   conservative upper bound on the fast speed.)
2. With self-gravity on, also the acceleration limit
   :math:`\Delta t \le 0.2\left(-|\boldsymbol{u}|/g + \sqrt{|\boldsymbol{u}|^2/g^2 + 2\min(\Delta x)/g}\right)`
   with :math:`g = |\boldsymbol{g}|` from the last gravity solve.
3. ``MPI_ALLREDUCE(…, MPI_MIN)`` over all ranks — every rank uses the same
   :math:`\Delta t`.
4. Clamp to the output clock: if :math:`t + \Delta t` would pass the next
   output time ``toutput``, shorten the step to land on it exactly, set
   ``writeFlag``, advance ``toutput`` by ``dtout`` and increment ``fnum``;
   likewise for ``tend``.

The result is stored in ``this%dt``; nothing else changes.

Step 2: the update — ``evolveGridRK2`` → ``rk2_3D``
====================================================

``evolveGridRK2`` dispatches on ``ndim`` to ``rk2_1D`` / ``rk2_2D`` /
``rk2_3D`` (``gridModule_rk2.f03``). The scheme is Heun's method (second-order
Runge–Kutta, TVD): two stages, then the average. Within a stage the flux
divergence is accumulated sweep by sweep on a *frozen* state (method of
lines), so a stage is one consistent evaluation of the right-hand side, not a
sequence of one-dimensional updates.

.. code-block:: text

   rk2_3D(this, q, q1, q2)                                        gridModule_rk2.f03
   │
   │  choose the sweep routine from (eosType, solverType):
   │     eosType=1: solverType 1,2 → solverIso3D      4,5 → solverIsoMHD3D
   │     eosType=2: solverType 2,3 → solverAdi3D      4,5 → solverAdiMHD3D     (riemannSolverModule.f03)
   │  reset the FOFC flags (fofc_flag3d = 0, fofc_pass = .false.), ndthalve = 0
   │
   ├─ do  ── the failsafe loop; a clean step runs it exactly once ───────────────────────────
   │   │
   │   │  STAGE 1  (q → q1)
   │   ├─ rieSolver(this, q, q,  q1, dd=1)      x-sweep: q1 = q  − Δt ∂F_x/∂x
   │   ├─ rieSolver(this, q, q1, q1, dd=2)      y-sweep: q1 = q1 − Δt ∂F_y/∂y   (fluxes from q)
   │   ├─ rieSolver(this, q, q1, q1, dd=3)      z-sweep: q1 = q1 − Δt ∂F_z/∂z, then the CT update of B
   │   ├─ [enable_sg] this%calcSelfgravity(q)   → gphi, sgfx, sgfy, sgfz from the stage state
   │   │              q1(mom) += ρ g Δt ;  q1(ene) += (ρu·g) Δt   (energy only for eosType=2)
   │   ├─ this%exchangeBdryMPI(q1, winq1)       MPI ghost exchange (gridModule_exchange.f03)
   │   ├─ this%setBoundary(q1)                  physical boundaries (gridModule_boundary.f03 → bdry<Case>)
   │   │
   │   │  STAGE 2  (q1 → q2)
   │   ├─ rieSolver(this, q1, q1, q2, dd=1..3)  same three sweeps from q1
   │   ├─ [enable_sg] this%calcSelfgravity(q1) → sources on q2
   │   ├─ scan q2 for min pressure (adiabatic)  → this%neg_pressure
   │   ├─ this%exchangeBdryMPI(q2, winq2);  this%setBoundary(q2)
   │   ├─ MPI_ALLREDUCE(changeSolver, neg_pressure, MPI_LOR)   → global flags
   │   │
   │   │  DECISION
   │   ├─ clean (no negative state anywhere)            → exit
   │   ├─ SCORPIO_LEGACY_FAILSAFE=1                     → global HLLD→HLL switch / global Δt halving, redo
   │   ├─ FOFC: pass ≤ 4 and solverType = 5             → flag cells with ρ≤0 or P≤0 in q2
   │   │        (ghost-inclusive), fofc_pass = .true., redo the step: faces touching a
   │   │        flagged cell use first-order HLL inside the solver (see below)
   │   ├─ FOFC exhausted, SCORPIO_DT_HALVE on, < 10 halvings → this%dt = dt/2, redo from stage 1
   │   └─ otherwise                                    → FATAL stop (never integrates P<0)
   │  enddo
   │
   ├─ q(interior) = ½ (q + q2)                         Heun average
   ├─ this%exchangeBdryMPI(q, winq);  this%setBoundary(q)
   └─ this%t = this%t + this%dt

Notes on the stages:

- The three sweeps of a stage all reconstruct from the *same* input state
  (first argument) and accumulate into the output state (third argument);
  the middle argument is the running accumulator. That is why the
  constrained-transport update can be assembled once per stage at the last
  sweep from the EMFs cached by all three sweeps (:ref:`ch:methods`).
- Gravity is evaluated once per stage from the stage's input state; the
  source terms are added to the stage output *after* the sweeps and before
  the ghost exchange. ``calcSelfgravity`` dispatches on ``sgSolverType`` and
  ``ndim`` to ``calcSG3D`` / ``calcSG3Dperiodic`` (FFT) or ``calcSG3D_MG``
  (multigrid), see :ref:`ch:selfgravity`.
- The turbulence kick is *not* part of the step: the driver applies it before
  ``evolveGridRK2`` on the ``DT_mode = 1`` schedule (:ref:`ch:turbulence`).
- ``rk2_2D`` has the same structure with two sweeps; ``rk2_1D`` has one sweep
  and no failsafe ladder (1D has no FOFC/dual energy).

Inside a sweep — ``solverAdiMHD3D``
===================================

Each ``rieSolver`` call is one directional sweep. Taking the 3D adiabatic MHD
solver (``riemannSolverModule.f03``) as the reference — the isothermal and
hydro variants follow the same skeleton with fewer variables:

.. code-block:: text

   solverAdiMHD3D(this, q, q1, q2, dd)                             riemannSolverModule.f03
   ├─ select the flux kernel:   solverType 4 → fluxPtr => fluxHLLAdiMHD1D
   │                            solverType 5 → fluxPtr => fluxHLLDAdiMHD1D
   ├─ select the limiter:       limiterType 0/1/2/3 → slope => zslop / vslop / fslop / minmod   (limiterModule.f03)
   ├─ sweep frame (cx,cy,cz) and the index map f1..f11 / signs coef(:) for dd = 1,2,3
   │     the sweep always works on "normal" and "transverse" components: for dd=2 the
   │     roles of x and y are exchanged (slots 2↔3, 5↔6, 9↔10), for dd=3 x and z (2↔4, 5↔7, 9↔11),
   │     with sign flips on the rotated transverse components to keep the frame right-handed
   ├─ primitive variables W = (ρ, u_n, u_t1, u_t2, B_n, B_t1, B_t2, P) per cell, P from E
   │     (or from the entropy variable in low-β cells — dual energy)
   ├─ slopes:  SL(i,j,k,m) = slope(W_{m-1}, W_m, W_{m+1}) for every variable  [PLM]
   │     [SCORPIO_PPM=1] CW84 face values + parabola monotonisation → SL, SLR (left/right excursions)
   ├─ positivity guard: scale a cell's slopes by θ∈[0,1] so both face states keep ρ>0, P>0
   ├─ for every face: ql = W_i, qr = W_{i+1};  UL = ql + ½ SL, UR = qr − ½ SLR
   │     fofc_pass and (cell i or i+1 flagged) → zero slopes + fluxHLLAdiMHD1D (first order)
   │     otherwise                             → fluxPtr(ql, qr, slopeL, slopeR, gamma, flux, nvar)
   │     entropy (dual energy): flux = (σ/ρ)_upwind × mass flux
   │     [AMR] capture the block-edge fluxes (amr_flo3/amr_fhi3) when amr_capture
   ├─ update:  q2(i,·) = q1(i,·) − Δt/Δx_n (F_{i+1/2} − F_{i-1/2})   for the hydro variables
   ├─ CT:      cache this sweep's face EMFs (E = −F_B components) and upwind weights;
   │           at dd = 3 assemble the corner EMFs (upwinded by default) and update the six
   │           face-B slots with the discrete curl                           (:ref:`ch:methods`)
   └─ dual energy: recover P from σ where e_int < thr·E and re-sync E, else re-sync σ from E

The flux kernels all share one interface,
``flux<name>(ql, qr, slopeL, slopeR, gam, flux, nvar)``: they build the
reconstructed states ``UL = ql + 0.5*slopeL`` and ``UR = qr − 0.5*slopeR``
and return the interface flux. ``fluxHLLDAdiMHD1D`` calls
``fluxHLLAdiMHD1D`` itself for the degenerate fan (and, with
``SCORPIO_HLLD_STARCHECK=1``, for inadmissible star states). The full list of
solvers and kernels is in :ref:`ch:riemann`.

Two-fluid (ion–neutral) drivers
===============================

For the ambipolar-diffusion cases the driver owns two grids, ``gn``
(neutrals, HD) and ``gi`` (ions, MHD), and calls one of the two-fluid
integrators instead of ``evolveGridRK2``. The 3D driver ``ADMHD3D``
(``testSuiteMPI.f03``) picks it from ``SCORPIO_AD_SCHEME``:

.. code-block:: text

   ADMHD3D time loop                                         testSuiteMPI.f03
   ├─ gn%griddt(); gi%griddt(); dt = min of both, set on both grids
   ├─ select case (SCORPIO_AD_SCHEME)
   │    split (default) → rk2AD_3D_HSHSMD(gn,…,gi,…)      rk2.f03
   │    imex            → rk2AD_3D_IMEX                    evolveAD_imex.f03
   │    imex322/athenak → rk2AD_3D_IMEX322                 evolveAD_imex.f03
   │    imexpp/imex2+   → rk2AD_3D_IMEXPP                  evolveAD_imex.f03
   └─ exchangeBdryMPI / setBoundary on both grids, writeGrid on both

   rk2AD_3D_HSHSMD (operator-split TR-BDF2)                  rk2.f03
   ├─ do (failsafe loop, as in rk2_3D)
   │   stage 1: rieSolvern(gn, qn, …, dd=1..3)   neutral sweeps (solverAdi3D / solverIso3D)
   │            rieSolveri(gi, qi, …, dd=1..3)   ion sweeps     (solverAdiMHD3D / solverIsoMHD3D)
   │            evolveAD3D_MD(gn, gi, qn, qi, qn1, qi1)   drag + heating, TR-BDF2   (evolveAmbipolarDiffusion.f03)
   │            post-drag dual-energy re-sync of the ion energy
   │            exchangeBdryMPI / setBoundary on qn1 and qi1
   │   stage 2: the same from (qn1, qi1) into (qn2, qi2)
   │   negative-state check on both fluids → FOFC / Δt halving / abort
   │  enddo
   └─ Heun average on both grids, ghost fill, t += dt

   rk2AD_3D_IMEXPP (Krapp et al. 2024 IMEX(4,3,2)) — same sweeps, but the drag is
   solved implicitly inside each stage by adImexDragStage(…) and the entropy is
   re-synchronised by ad_imex_resync(…); IMEX and IMEX322 differ only in the
   Butcher tableau (:ref:`ch:AD`).

AMR runs
========

On the adaptive mesh the driver calls ``amrDt`` and ``amrStep``
(``amrModule.f03``) instead of ``griddt`` / ``evolveGridRK2``. ``amrStep``
dispatches to ``amrStep3D``, which applies the same Heun scheme to every leaf
block with the same sweep routines, plus the inter-level operations:

.. code-block:: text

   amrDt(m)         global CFL minimum over all leaf blocks and ranks, output clamp (as dt3D)
   amrStep(m) → amrStep3D(m)
   ├─ stage 1: amrGhostFill(q)  → amrExchangeLevel3D (same-level copies, coarse→fine prolongation),
   │                              amrApplyBCLevel (physical boundaries)
   │           per local block:  rieSolver(blk%g, q, q, q1, dd=1..3) with amr_capture = .true.
   │                              (captureCopyOut3, emfAccum3 store the block-edge fluxes / EMFs)
   │           amrReflux3 (coarse–fine flux and EMF matching), restriction fine→coarse
   │           [gravOn] amrSolveGravity3 → amrGatherDensity → gio%calcSelfgravity → amrScatterForce;
   │                    amrApplyGravSource3
   │           amrGhostFill(q1)
   ├─ stage 2: the same into q2
   └─ Heun average per block, restriction, amrGhostFill(q)

   driver:  every regridEvery steps → amrRegrid3 (Löhner/Jeans flags, refine/derefine, move data,
            rebuild the exchange/reflux/IO plans);  on writeFlag → amrOutput (base-grid projection
            + checkpoint)

Details of the AMR machinery are in :ref:`ch:methods`.

Restart
=======

``setTime(fstart = N)`` (``gridModule.f03``) is the restart entry: it reads
``g<gridID>_<N>.h5`` with ``read1d`` / ``read3d`` (``gridModule_io*.f03``),
calling ``setMesh``, ``setVariable`` and ``setMPIWindows`` itself, restores
``t`` and the solver settings, fills the fields, exchanges the ghost zones and
re-seeds the dual-energy entropy. After it the driver continues with
``initSelfgravitySolver`` and the time loop exactly as for a fresh start; see
the Problem File page for the driver-side flag and the turbulence-driving
bookkeeping.
