.. _ch:code_structure:

**************
Code Structure
**************

.. warning:: Column-major order is used in Fortran: in ``q(i,j,k,n)`` the
   first index varies fastest in memory, and the arrays come out of ``h5py``
   transposed, as ``(n, k, j, i)``.

Source files
============

All sources are in ``src/`` and are built by the top-level ``Makefile``
(``make``; the executable is ``./Scorpio``). One module, ``gridModule``,
holds the grid object; its implementation is split into *submodules* (one
file per topic), which is why most files are called ``gridModule_*.f03``.

.. list-table::
   :header-rows: 1
   :widths: 30 12 58

   * - file
     - unit
     - contents
   * - ``main.f03``
     - program
     - ``MPI_INIT``, ``setMPI(np)``, ``runProblemFromNamelist``, dual-energy fire-rate summary, ``MPI_FINALIZE``
   * - ``problemRegistry.f03``
     - module
     - reads ``problem.nml``; ``dispatchProblemByGridID`` (the ``select case`` on ``gridID``) and ``dispatchProblemByName``; the ``applyProblemInit*`` / ``applyProblemBoundary*`` hooks
   * - ``scorpio_kinds.f03``
     - module
     - kind parameters
   * - ``limiterModule.f03``
     - module
     - slope limiters ``zslop``, ``vslop`` (van Leer), ``fslop`` (MC), ``minmod``
   * - ``gridModule.f03``
     - module
     - the ``grid`` derived type and its type-bound procedures (interfaces of all ``module subroutine`` s), run-option globals (dual energy, FOFC, PPM, HLLD, EMF), ``setVariable``, ``setTime`` (restart reader), ``writeGrid``, ``setTopologyMPI``, gravity enable/init, driving enable, ``TrueloveCondition``, ``setRunOption``
   * - ``gridModule_setters.f03``
     - submodule
     - one-line setters: ``setGridID``, ``setEoS``, ``setSoundSpeed``, ``setAdiGamma``, ``setCFL``, ``setSolverType``, ``setSlopeLimiter``, ``setBoundaryType``, ``enableAD``, ``setADparams``
   * - ``gridModule_coordinates.f03``
     - submodule
     - ``setCoordinates``: cell faces/centres/widths from ``leftBdry``, ``rightBdry``, ``nMesh``, ``nbuf`` (Cartesian and the two cylindrical options)
   * - ``gridModule_windows.f03``
     - submodule
     - ``initMPIWindows1D/2D/3D``: one-sided MPI windows on ``q``, ``q1``, ``q2``; ``initSGWindows2D/3D``
   * - ``gridModule_exchange.f03``
     - submodule
     - ``exchgBdryMPI1D/2D/3D``: ghost-zone exchange with ``MPI_WIN_FENCE`` + ``MPI_GET``/``MPI_PUT``
   * - ``gridModule_boundary.f03``
     - submodule
     - ``setBdry1D/2D/3D``: physical boundaries — registry hook first, then ``select case (gridID)`` → the case's ``bdry<Case>`` routine
   * - ``gridModule_init.f03``
     - submodule
     - ``init1d/2d/3d``: initial condition — registry hook first, then ``select case (gridID)`` → ``init<Case>``; seeds the dual-energy entropy
   * - ``gridModule_dt.f03``
     - submodule
     - ``dt1D/2D/3D``: CFL time step, gravity limit, global minimum, output clamp
   * - ``gridModule_rk2.f03``
     - submodule
     - ``rk2_1D/2D/3D``: the single-fluid Heun step with gravity sources and the failsafe ladder (:ref:`ch:time_integration`)
   * - ``gridModule_source.f03``
     - submodule
     - ``source2D``: per-case external source hook called by ``rk2_2D`` (e.g. the Rayleigh–Taylor gravity, cases 44/45)
   * - ``gridModule_io1d/2d/3d.f03``
     - submodules
     - parallel HDF5 writers/readers ``output1d/2d/3d``, ``read1d/2d/3d`` (hyperslab per rank)
   * - ``gridModule_pvts.f03``
     - submodule
     - VTK/PVTS writers and readers (``writeGrid_HD_vtk``, ``writeGrid_MHD_vtk``, ``readGrid_*_vtk``)
   * - ``gridModule_fft_plans.f03``
     - submodule
     - FFTW-MPI plans ``sgPlan2D/3D`` (gravity) and ``DTPlan2D/3D`` (driving), created in ``setVariable``
   * - ``gridModule_sgcalc.f03``
     - submodule
     - FFT self-gravity: ``calcSG2D/3D`` (isolated), ``calcSG2D/3Dperiodic``, the hydro↔FFT-slab remaps (``remap_*``)
   * - ``calcSG_MG.f03``
     - submodule
     - multigrid self-gravity ``calcSG2D_MG`` / ``calcSG3D_MG``: V-cycles, James and multipole boundaries, periodic gauge
   * - ``exchgBdryMPI_sgMG.f03``
     - module
     - ghost exchange for the multigrid levels
   * - ``gridModule_dtcalc.f03``
     - submodule
     - turbulence driving: ``calcDT2D/3D_MD`` (spectrum, projection, remap, momentum/energy normalisation), driving spectra
   * - ``riemannSolverModule.f03``
     - module
     - the directional sweep solvers ``solver{Iso,IsoMHD,Adi,AdiMHD,Poly}{1D,2D,3D}`` with reconstruction (PLM/PPM), positivity guard, FOFC faces, CT EMF assembly, dual energy; the flux kernels ``flux*1D``
   * - ``rk2.f03``
     - external procedures
     - two-fluid drivers: ``rk2AD_1D/2D/3D``, ``rk2AD_3D_HSHSMD`` / ``rk2AD_2D_HSHSMD`` (operator-split TR-BDF2, the production ones), the older ``rk2ADsg_*`` variants, and ``rk2MHD_3D`` (a single-fluid variant used by two Richtmyer–Meshkov cases)
   * - ``evolveAmbipolarDiffusion.f03``
     - external procedures
     - the ion–neutral drag/heating update ``evolveAD1D/2D/3D``, ``evolveAD3D_MD`` (Tilley et al. 2012)
   * - ``evolveAD_imex.f03``
     - external procedures
     - IMEX two-fluid drivers ``rk2AD_3D_IMEX``, ``rk2AD_3D_IMEX322``, ``rk2AD_3D_IMEXPP`` with ``adImexDragStage``, ``ad_imex_resync``
   * - ``amrModule.f03``
     - module
     - block-structured AMR: tree, ghost fill, prolongation/restriction, reflux, EMF matching, regridding, base-grid and FAC gravity, projected output, checkpoints
   * - ``testSuiteMPI.f03``
     - module
     - the standard test cases: for every case a driver, an ``init<Case>`` and a ``bdry<Case>`` routine
   * - ``ShiboTestSuite.f03``
     - module
     - additional cases (``load_struct03``, ``MHD3DTurbDriven``, ``clumpRerun``)
   * - ``testSuiteAMR.f03``
     - module
     - the AMR test cases (700–727) and the ``&amr_*`` namelist readers
   * - ``spectrumCompensation.f03``
     - module
     - the spectrum-compensated driven-turbulence case (401)
   * - ``hinnyCloud.f03``
     - module
     - the 20 pc cloud: uniform driver ``cloud_20pc3_3DMHD``, its ``initcloud_``/``bdrycloud_`` routines, ``&cloud_nml``, the AMR driver ``cloud_20pc3_3DMHD_AMR``

Module dependency order (this is the order the Makefile compiles the
module files in): ``scorpio_kinds`` → ``limiterModule`` → ``gridModule`` →
``riemannSolverModule`` → ``exchangeBdryMPI_sgMG`` → ``testSuiteMPI`` →
``amrModule`` → ``testSuiteAMR`` → ``problemRegistry``; the case modules
(``ShiboTestSuite``, ``HinnyTestSuite``, ``SpectrumCompensationSuite``) sit
between ``gridModule`` and ``problemRegistry``, and the submodules after
``gridModule``.

Program flow
============

.. code-block:: text

   main.f03
     MPI_INIT → setMPI(np) [sets nprocs, myid] → setTestOnOff(.true.)
     runProblemFromNamelist                             problemRegistry.f03
       read &problem_config from problem.nml
       dispatchProblemByGridID(gridID)  or  dispatchProblemByName(problem_name)
         → one problem driver, e.g. cloud_20pc3_3DMHD (hinnyCloud.f03)
             grid setup through the grid API            (Problem File page, "Grid setup API")
             time loop: griddt → [driving] → evolveGridRK2 → writeGrid   (Time Integration page)
     dual-energy fire-rate summary → MPI_FINALIZE

Two dispatch tables besides the registry decide what a case does: the
initial condition (``init3d`` in ``gridModule_init.f03``) and the physical
boundary (``setBdry3D`` in ``gridModule_boundary.f03``) both do
``select case (gridID)`` and call the case's ``init<Case>`` /
``bdry<Case>`` routines — after first offering the ``applyProblemInit3D`` /
``applyProblemBoundary3D`` hooks of the registry (used by the ``problem_name``
cases). So a new case needs an entry in three places: ``problemRegistry.f03``
(driver), ``gridModule_init.f03`` (initial condition) and
``gridModule_boundary.f03`` (boundary), plus the ``use`` lines that import
the routines.

The grid object
===============

Everything about one mesh lives in one ``type(grid)`` variable (``g1`` in the
drivers; ``gn`` and ``gi`` for the two fluids). The components you will meet
when reading the code:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - component
     - meaning
   * - ``nMesh(ndim)``, ``nMesh_global(ndim)``, ``meshStart(ndim)``
     - local cells per direction, global cells, offset of this rank's block in the global mesh
   * - ``nbuf``, ``ndim``, ``nvar``, ``variable(8)``, ``coordType``, ``gridID``
     - ghost width, dimensions, number of stored variables, the variable flags, coordinate system, case ID
   * - ``leftBdry``, ``rightBdry`` (+ ``_global``)
     - box edges of the local block and of the whole domain
   * - ``xl(d)%coords``, ``xr(d)%coords``, ``xc(d)%coords``, ``dx(d)%coords``
     - left face, right face, centre and width of every cell along direction ``d``, indexed ``1-nbuf : nMesh(d)+nbuf``
   * - ``q``, ``q1``, ``q2``
     - the state and the two RK2 stage states. Stored as flat 1D arrays and passed to the sweeps as 4D
       ``(1-nbuf:nx+nbuf, 1-nbuf:ny+nbuf, 1-nbuf:nz+nbuf, nvar)`` by sequence association; slot order
       in :ref:`ch:hydro`
   * - ``winq``, ``winq1``, ``winq2``, ``databuf1/2``
     - MPI one-sided windows on the three states and the exchange buffers
   * - ``t``, ``dt``, ``tend``, ``dtout``, ``toutput``, ``fnum``, ``fstart``, ``writeFlag``
     - the clock and the output bookkeeping
   * - ``CFL``, ``snd``, ``adiGamma``, ``eosType``, ``solverType``, ``limiterType``, ``boundaryType``
     - the numerical settings (also stored in every snapshot)
   * - ``enable_sg``, ``sgSolverType``, ``sgBdryType``, ``GravConst``, ``gphi``, ``sgfx/y/z``
     - self-gravity switch, solver (0 FFT / 1 MG), boundary (0 isolated / 1 periodic), :math:`G`, potential and acceleration
   * - ``sg*Kernel``, ``sg*Cmplx``, ``sgPlan*``, ``sgDensityBuffer``
     - FFT gravity work arrays and plans
   * - ``enable_DT``, ``DT_mode``, ``drivingWN_DT``, ``Energy_DT``, ``zeta_DT``, ``netmom*_DT``, ``DTenergyfaction``
     - turbulence-driving switch and parameters
   * - ``enable_ad``, ``mu_ad``, ``alpha_ad``
     - two-fluid coupling parameters (set on both grids)
   * - ``vu_mpi``, ``mpiCoord``, ``dims_mpi``, ``left/right/up/down/top/bottom_mpi``
     - Cartesian communicator, this rank's coordinates and its six neighbours
   * - ``changeSolver``, ``neg_pressure``
     - failsafe flags raised by the sweeps and consumed by ``rk2_*``

Where to look for what
======================

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - I want to change …
     - look in
   * - the initial condition of a case
     - ``init<Case>`` in the case's module (cloud: ``initcloud_20pc3_3DMHD`` in ``hinnyCloud.f03``)
   * - the physical boundary of a case
     - ``bdry<Case>`` (called from ``setBdry3D``)
   * - the time step
     - ``dt3D`` in ``gridModule_dt.f03``
   * - the update / sources / failsafes
     - ``rk2_3D`` in ``gridModule_rk2.f03``
   * - reconstruction, Riemann flux, CT
     - ``solverAdiMHD3D`` / ``solverIsoMHD3D`` and the ``flux*`` kernels in ``riemannSolverModule.f03``
   * - a limiter
     - ``limiterModule.f03``
   * - gravity
     - ``gridModule_sgcalc.f03`` (FFT), ``calcSG_MG.f03`` (multigrid), ``amrModule.f03`` (AMR coupling, FAC)
   * - turbulence driving
     - ``gridModule_dtcalc.f03``
   * - what is written to the snapshot
     - ``writeGrid`` in ``gridModule.f03`` and ``output3d`` in ``gridModule_io3d.f03``
   * - the restart reader
     - ``setTime`` in ``gridModule.f03``
   * - MPI decomposition and ghost exchange
     - ``setTopologyMPI`` (``gridModule.f03``), ``gridModule_windows.f03``, ``gridModule_exchange.f03``
   * - the run-option switches (``SCORPIO_*``)
     - top of ``gridModule.f03`` (declarations and documentation) and ``setVariable`` (where they are read)
   * - AMR
     - ``amrModule.f03``, design notes in ``docs/AMR_DESIGN.md`` and ``docs/FAC_GRAVITY_DESIGN.md``

Conventions in the source
=========================

- Comments tagged ``[OPUS]`` / ``[FABLE]`` mark the 2026 additions (dual
  energy, FOFC, PPM, CT, IMEX, AMR) and usually carry the date and the reason;
  ``CHANGES.md`` in the repository root is the narrative index of them.
- Arrays with ghost cells are always declared ``1-nbuf : nMesh+nbuf``;
  interior loops run ``1 : nMesh``.
- Real constants are written ``1.d0``; every routine has ``implicit none``.
- Rank 0 prints; use ``if (myid .eq. 0)`` for new diagnostics.
