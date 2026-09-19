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
the routines. A complete worked example is in :ref:`sec:new_case`.

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
     - ``amrModule.f03`` (:ref:`ch:methods` for how it works; the design notes ``docs/AMR_DESIGN.md`` and ``docs/FAC_GRAVITY_DESIGN.md`` live in the code repository)

.. _sec:new_case:

Adding a new case
=================

A case is three routines plus three registrations. The example adds a 3D
isothermal MHD case with ``gridID = 802`` called ``myCloud``; put the
routines into an existing case module (here ``HinnyTestSuite`` in
``src/hinnyCloud.f03``) so that no Makefile change is needed.

**1. The driver** — builds the grid, chooses the physics, runs the loop.
This is the complete minimal version (fresh start only; the cloud driver
shows the optional pieces: gravity, driving, restart, Truelove stop):

.. code-block:: fortran

   subroutine myCloud(gridID)
       integer, intent(in) :: gridID
       type(grid) :: g1
       integer :: ndim, nbuf, coordType, variable(8), nMesh(3), dims(3), ierr
       double precision :: leftBdry(3), rightBdry(3)
       logical :: periods(3), reorder

       ndim = 3;  nbuf = 2;  coordType = 1                  ! 3D, two ghost cells, Cartesian
       variable = 1                                         ! den, momx/y/z, bxl/byl/bzl, ene -> MHD
       nMesh = (/32, 32, 64/)
       leftBdry  = (/-5.d0, -5.d0, -10.d0/)                 ! pc
       rightBdry = (/ 5.d0,  5.d0,  10.d0/)
       dims = 0;  call MPI_DIMS_CREATE(nprocs, ndim, dims, ierr)
       periods = .true.;  reorder = .true.

       call g1%setGridID(gridID = gridID)
       call g1%setTopologyMPI(ndim, dims, periods, reorder)
       call g1%setMesh(nMesh, leftBdry, rightBdry, nbuf, coordType, gridID)
       call g1%setVariable(variable)
       call g1%setMPIWindows()
       call g1%setEoS(eosType = 1);  call g1%setSoundSpeed(snd = 0.3d0)   ! isothermal, 0.3 km/s
       call g1%setCFL(CFL = 0.4d0)
       call g1%setSlopeLimiter(limiterType = 3)                          ! minmod
       call g1%setSolverType(solverType = 5)                             ! HLLD
       call g1%setBoundaryType(boundaryType = 3)                         ! periodic
       call g1%setTime(fstart = 0, tend = 1.d0, dtout = 0.1d0)
       call g1%initVariable()                                            ! -> initmyCloud (via init3d)
       call g1%exchangeBdryMPI(g1%q, g1%winq)
       call g1%setBoundary(g1%q)                                         ! -> bdrymyCloud (via setBdry3D)
       call g1%writeGrid()                                               ! g0802_0000.h5
       g1%writeFlag = .false.

       do while (g1%t .lt. g1%tend)
           call g1%griddt()
           call g1%evolveGridRK2()
           if (myid .eq. 0) print *, "myCloud: t =", g1%t, " dt =", g1%dt
           if (g1%writeFlag) then
               call g1%writeGrid()
               g1%writeFlag = .false.
           end if
       end do
   end subroutine myCloud

**2. The initial condition** — fills the interior cells of ``q`` on this
rank; the coordinates ``this%xc(d)%coords`` are already the global positions
of the local cells:

.. code-block:: fortran

   subroutine initmyCloud(this, q)
       class(grid) :: this
       double precision, dimension(1-this%nbuf:this%nMesh(1)+this%nbuf, 1-this%nbuf:this%nMesh(2)+this%nbuf, &
                                   1-this%nbuf:this%nMesh(3)+this%nbuf, this%nvar) :: q
       integer :: i, j, k
       double precision :: x, y, z, rho, b0

       b0 = 35.d0                                            ! code units (x 2.9 for microgauss)
       do k = 1, this%nMesh(3)
           do j = 1, this%nMesh(2)
               do i = 1, this%nMesh(1)
                   x = this%xc(1)%coords(i);  y = this%xc(2)%coords(j);  z = this%xc(3)%coords(k)
                   rho = 1.d0 + 10.d0*exp(-(x*x + y*y + z*z))
                   q(i,j,k,1)   = rho                       ! density        [Msun/pc^3]
                   q(i,j,k,2:4) = 0.d0                      ! momentum = rho * velocity
                   q(i,j,k,5) = 0.d0;  q(i,j,k,9)  = 0.d0   ! Bx on the left face / right face
                   q(i,j,k,6) = 0.d0;  q(i,j,k,10) = 0.d0   ! By
                   q(i,j,k,7) = b0;    q(i,j,k,11) = b0     ! Bz: uniform field along z
                   q(i,j,k,8)   = 0.d0                      ! energy (isothermal: not used)
               end do
           end do
       end do
   end subroutine initmyCloud

Set both face slots of every field component (:ref:`sec:conserved`); the
ghost cells are filled afterwards by ``exchangeBdryMPI`` / ``setBoundary``.

**3. The boundary routine** — ``setBdry3D`` calls it after the MPI
exchange on every ghost fill. It is entirely the case's responsibility:
``boundaryType`` is only a label the routine may test, and ``left_mpi < 0``
(and ``right_mpi``, ``up_mpi``, ``down_mpi``, ``top_mpi``, ``bottom_mpi``)
tells whether this rank sits at a physical edge. With a periodic box every
face has an MPI neighbour and the routine does nothing:

.. code-block:: fortran

   subroutine bdrymyCloud(this, q)
       class(grid) :: this
       double precision, dimension(1-this%nbuf:this%nMesh(1)+this%nbuf, 1-this%nbuf:this%nMesh(2)+this%nbuf, &
                                   1-this%nbuf:this%nMesh(3)+this%nbuf, this%nvar) :: q
       integer :: i, nx

       if (this%boundaryType .ne. 1) return                  ! periodic: the MPI exchange did it
       nx = this%nMesh(1)
       if (this%left_mpi .lt. 0) then                        ! zero gradient at the physical x edges
           do i = 1, this%nbuf
               q(1-i, :, :, :) = q(1, :, :, :)
           end do
       end if
       if (this%right_mpi .lt. 0) then
           do i = 1, this%nbuf
               q(nx+i, :, :, :) = q(nx, :, :, :)
           end do
       end if
   end subroutine bdrymyCloud

**4. Register the three routines** — one line in each of three files:

.. code-block:: fortran

   ! src/problemRegistry.f03 — the driver
   use HinnyTestSuite, only: cloud_20pc3_3DMHD, cloud_20pc3_3DMHD_AMR, myCloud
   ...
       case(802)
         call myCloud(gridID=802)

   ! src/gridModule_init.f03, subroutine init3d — the initial condition
   use HinnyTestSuite, only: initcloud_20pc3_3DMHD, initmyCloud
   ...
       case(802)
         call initmyCloud(this,q)

   ! src/gridModule_boundary.f03, subroutine setBdry3D — the boundary
   use HinnyTestSuite, only: bdrycloud_20pc3_3DMHD, bdrymyCloud
   ...
       case(802)
         call bdrymyCloud(this,q)

Then ``make``, put ``gridID = 802`` in ``problem.nml`` and run. Two-fluid
cases register a *pair* of IDs (``case(26, 27)`` calling one driver with
``gridIDn`` and ``gridIDi``) and provide ``init``/``bdry`` routines for each
grid.

If the routines go into a **new file**, that file is a new module: add it
to ``MOD_SRCS`` in the ``Makefile`` and give it the dependency lines the
existing case modules have (``problemRegistry.o``, ``gridModule_init.o`` and
``gridModule_boundary.o`` must depend on it, and it on ``gridModule.o``).

Checklist: a ``gridID`` nobody else uses; ``variable(5:7)`` all 1 or all 0,
``variable(8) = 1`` for any MHD run; ``nMesh`` divisible by the rank counts
you will use; ``periods`` matching ``boundaryType``; a ``.h5`` file appears
after the first step.

.. _sec:source_terms:

Adding a source term
====================

Where a source goes depends on how fast it acts and whether it must be part
of the Runge–Kutta stages:

.. list-table::
   :header-rows: 1
   :widths: 26 34 40

   * - kind of source
     - where
     - pattern in the code
   * - body force / heating that is smooth on the time step (gravity, an external potential, a cooling function)
     - inside ``rk2_3D`` (``gridModule_rk2.f03``), once per stage, after the three sweeps and before ``exchangeBdryMPI``; evaluated from the stage's *input* state, added to the stage's *output* state; interior cells only
     - self-gravity: ``q1(mom) += q(rho)*g*dt``, ``q1(ene) += (q(mom)·g)*dt`` in stage 1, the same from ``q1`` into ``q2`` in stage 2
   * - the same, but per case and in 2D
     - ``source2D`` (``gridModule_source.f03``): ``rk2_2D`` calls ``source2D(this, q, q1)`` and ``source2D(this, q1, q2)`` after the sweeps of each stage; it dispatches on ``gridID``
     - the Rayleigh–Taylor gravity (cases 44/45). There is no ``source3D`` yet — adding one means a ``module subroutine source3D`` in the same file, its interface in ``gridModule.f03``, and the two calls in ``rk2_3D``
   * - stiff relaxation (a rate faster than the time step, e.g. the ion–neutral drag)
     - either operator-split after the sweeps of each stage (``evolveAD3D_MD``, called from ``rk2AD_3D_HSHSMD``) or implicitly inside the stage (``adImexDragStage`` in the IMEX drivers)
     - solve the stiff term per cell exactly or implicitly; the IMEX placement avoids the splitting error when the relaxation time is much shorter than ``dt``
   * - something that changes the state between steps (a turbulence kick, a re-seeding)
     - in the driver's time loop, before ``evolveGridRK2``
     - the ``DT_mode = 1`` kick: modify ``g%q``, then ``exchangeBdryMPI`` + ``setBoundary``

Rules that apply to all of them:

- Work on the interior ``1:nMesh`` and refresh the ghost zones afterwards;
  never write into the ghost cells directly.
- Add the source to both stages, from the right input state (``q`` → ``q1``,
  then ``q1`` → ``q2``); the Heun average then makes it second order.
- Momentum sources need the matching energy source in adiabatic runs
  (``rho u · a`` for a force); isothermal runs carry no usable energy.
- In adiabatic MHD runs the dual-energy entropy ``σ = P/ρ^(γ-1)`` (the last
  slot of ``q``) is re-synchronised inside the sweeps. A source that changes
  the internal energy *after* the sweeps must update σ as well, otherwise
  low-β cells — which take their pressure from σ — will not see the heating.
  The two-fluid drivers do exactly this after the drag step ("post-drag
  re-sync").
- If the source has its own time scale, add its limit to ``dt3D``, as the
  gravity acceleration limit is.
- Keep it deterministic and rank-independent: no random numbers without a
  reproducible seed, no dependence on the decomposition.

.. _sec:globals:

Run-option globals in ``gridModule``
====================================

Besides the ``grid`` type, ``gridModule.f03`` declares module-level variables
that act on every grid. Most are the run options a user sets through
``SCORPIO_*`` environment variables (or ``setRunOption`` from a driver);
they are read **once**, in ``setVariable`` (and again in ``amrConfigure3D``
for AMR runs), which is why phase-1 calls must precede it.

.. list-table::
   :header-rows: 1
   :widths: 22 12 22 44

   * - variable
     - default
     - set by
     - read in / effect
   * - ``nprocs``, ``myid``
     - —
     - ``setMPI`` (``main.f03``)
     - everywhere; ``if (myid .eq. 0)`` guards prints
   * - ``dual_energy_on``
     - ``.true.``
     - ``SCORPIO_DUAL_ENERGY``
     - the entropy branch in ``solverAdiMHD2D/3D``, the two-fluid drivers, ``amrModule``
   * - ``dual_energy_thr``
     - 1e-3
     - ``SCORPIO_DE_THR``
     - same: fire when ``e_int < thr * E``
   * - ``dual_energy_nfire``, ``dual_energy_ntot``
     - 0
     - counters
     - summed and printed by ``main.f03`` at the end of the run
   * - ``de_print``
     - ``.true.``
     - ``SCORPIO_DE_PRINT``
     - the rank-local ``[dual-energy]`` line when cells are rescued
   * - ``health_monitor``, ``de_rescued``
     - ``.false.``, 0
     - ``SCORPIO_HEALTH``
     - ``opus_report_health`` (one line per eventful step, one ``MPI_ALLREDUCE``)
   * - ``use_legacy_failsafe``
     - ``.false.``
     - ``SCORPIO_LEGACY_FAILSAFE``
     - the decision block of ``rk2_2D/3D``: global HLL switch / global dt halving instead of FOFC
   * - ``dt_halve_on``
     - ``.true.``
     - ``SCORPIO_DT_HALVE``
     - the dt-halving rung of the failsafe ladder (``rk2_*``, two-fluid and IMEX drivers)
   * - ``disable_failsafe``
     - ``.false.``
     - set in code by the raw-diagnostic case (``SCORPIO_RAW``)
     - ``rk2_*`` accept every step as is (no FOFC, no retry)
   * - ``fofc_pass``, ``fofc_npass``, ``fofc_flag`` (2D), ``fofc_flag3d`` (3D)
     - ``.false.``, 0, unallocated
     - the RK drivers
     - read by the sweep solvers: faces touching a flagged cell use first-order HLL during a redo
   * - ``use_ppm``, ``ppm_warned``
     - ``.false.``
     - ``SCORPIO_PPM``
     - reconstruction branch of the 2D/3D MHD sweeps (needs ``nbuf = 3``; one warning otherwise)
   * - ``use_upwind_emf``
     - ``.true.``
     - ``SCORPIO_UPWIND_EMF``
     - CT corner-EMF assembly in the MHD sweeps (``0`` = centred average)
   * - ``hlld_eps``
     - 1e-8
     - ``SCORPIO_HLLD_EPS``
     - degeneracy threshold inside ``fluxHLLDAdiMHD1D`` / ``fluxHLLDIsoMHD1D``
   * - ``hlld_star_check``, ``hlld_fallback``, ``hlld_degen``
     - ``.false.``, 0, 0
     - ``SCORPIO_HLLD_STARCHECK``; counters
     - per-face admissibility check with HLL fallback; counts reported by the health monitor
   * - ``ad_imex``, ``ad_imex322``, ``ad_imexpp``
     - ``.false.``
     - ``SCORPIO_AD_SCHEME``
     - which two-fluid driver ``ADMHD3D`` calls (split / IMEX variants)
   * - ``amr_capture``, ``amr_flo``/``amr_fhi``, ``amr_flo3``/``amr_fhi3``, ``amr_e*``
     - ``.false.``, unallocated
     - ``amrStep`` / ``amrStep3D``
     - when true the sweeps copy their block-edge fluxes and edge EMFs into these arrays for refluxing and EMF matching; inert on uniform grids
   * - ``SG_SOLVER_FFT/MG``, ``SG_BDRY_ISOLATED/PERIODIC``
     - 0/1, 0/1
     - parameters
     - the values of ``sgSolverType`` and ``sgBdryType``

Two things follow from this layout. First, an option changed in the
environment after ``setVariable`` has run has no effect — set it before
launching, or from the top of the driver with ``setRunOption``. Second, the
flags are shared by every grid in the process, so in a two-fluid run both
fluids necessarily use the same reconstruction, EMF and dual-energy settings.

Conventions in the source
=========================

- Comments tagged ``[OPUS]`` / ``[FABLE]`` mark the 2026 additions (dual
  energy, FOFC, PPM, CT, IMEX, AMR) and usually carry the date and the reason;
  ``CHANGES.md`` in the repository root is the narrative index of them.
- Arrays with ghost cells are always declared ``1-nbuf : nMesh+nbuf``;
  interior loops run ``1 : nMesh``.
- Real constants are written ``1.d0``; every routine has ``implicit none``.
- Rank 0 prints; use ``if (myid .eq. 0)`` for new diagnostics.
