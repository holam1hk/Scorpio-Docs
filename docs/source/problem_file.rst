.. _ch:problem_file:

************
Problem File
************

A Scorpio run is configured on three levels:

1. **The problem driver** — a Fortran subroutine in ``src/`` (for example
   ``cloud_20pc3_3DMHD`` in ``src/hinnyCloud.f03``, or the cases in
   ``src/testSuiteMPI.f03``). It sets the mesh, the physics switches and the
   initial condition through the ``grid`` API described below. Changing it
   requires ``make``.
2. **problem.nml** — a Fortran namelist file in the directory where you run.
   It selects *which* driver runs (``&problem_config``) and, for drivers that
   read one, a per-case group of overrides (for example ``&cloud_nml``,
   ``&amr_config``). No recompilation.
3. **Run options** — ``SCORPIO_*`` environment variables (or the in-code
   equivalent ``setRunOption``) that switch numerical options such as the
   reconstruction, dual energy, the failsafe ladder, or the ion–neutral
   coupling scheme. No recompilation.

All names below are the ones in the current source; the tables list the value
each setting takes in the 20 pc cloud driver as "cloud default" where that is
useful.


Selecting the case: ``problem.nml``
===================================

The executable reads ``problem.nml`` from the current working directory and
writes all output there. The minimum is:

.. code-block:: fortran

   &problem_config
     gridID = 800
   /

.. list-table::
   :header-rows: 1
   :widths: 28 22 50

   * - entry
     - values
     - notes
   * - ``gridID`` (alias ``grid_id``)
     - integer
     - Case number; dispatched by ``src/problemRegistry.f03``. Also used in the
       output file names ``g<gridID>_<fnum>.h5``. An unknown value prints the
       list of supported IDs and aborts. Full map: ``docs/LegacyGridID.md``.
   * - ``problem_name``
     - ``'sod'``, ``'load_struct03'``, ``'load_struct03_withAD'``
     - Named cases; used only when ``gridID`` is absent or ≤ 0.
   * - ``gridIDn``, ``gridIDi`` (aliases ``grid_id_n``, ``grid_id_i``)
     - integer
     - Neutral / ion grid IDs for the two-fluid named case
       (``load_struct03_withAD``, defaults 1691 / 1692).

Frequently used ``gridID`` values:

.. list-table::
   :header-rows: 1
   :widths: 16 84

   * - ``gridID``
     - case
   * - 800
     - 20 pc\ :sup:`3` magnetized cloud, uniform grid (``cloud_20pc3_3DMHD``); reads ``&cloud_nml``
   * - 801
     - the same cloud on the block-structured AMR mesh (``cloud_20pc3_3DMHD_AMR``); reads ``&cloud_nml`` and the ``&amr_*`` groups
   * - 1–61
     - the standard 1D/2D/3D HD and MHD test suite (shock tubes, Orszag–Tang, blast waves, field loop, KH/RT, linear waves, …)
   * - 22–29, 52–57, 122–123, 500–501, 600–601, 609/616, 620/621
     - two-fluid ion–neutral (ambipolar diffusion) cases; each is a *pair* of grid IDs (neutrals, ions)
   * - 100, 399–401, 1688
     - driven-turbulence cases
   * - 21, 35, 49, 357–361, 557–558
     - self-gravity tests (FFT and multigrid, isolated and periodic, analytic gates)
   * - 700–727
     - AMR test cases (``src/testSuiteAMR.f03``)

.. note::

   Every group in ``problem.nml`` must end with ``/``. A misspelled entry name
   makes Fortran reject the *whole* group silently in the per-case readers, so
   check the start-up echo (for the cloud case the ``HinnyTestSuite: nMesh =``
   line) to confirm your values arrived. If the file was edited on Windows,
   strip the carriage returns (``sed -i -e 's/\r$//' -e '$a\' problem.nml``);
   the registry prints a hint when it detects this.


Numerical setup
===============

These integers are set in the driver through ``setMesh`` (``coordType``),
``setBoundaryType``, ``setEoS``, ``setSolverType`` and ``setSlopeLimiter``.
For a normal user the coordinate system, boundary condition and equation of
state must be chosen; the Riemann solver and limiter can be left at the cloud
defaults (HLLD, minmod).

.. _sec:coordtype:

Coordinate system
-----------------

.. list-table::
   :header-rows: 1
   :widths: 40 20 40

   * - system
     - value
     - notes
   * - Cartesian
     - 1
     - default; required for AMR, PPM and dual energy
   * - cylindrical, logarithmic :math:`r`
     - 2
     - 2D polar cases (``HDBlastWavePolar2D`` …)
   * - cylindrical, uniform :math:`r`
     - 3
     -

Boundary condition
------------------

.. list-table::
   :header-rows: 1
   :widths: 40 20 40

   * - condition
     - value
     - notes
   * - user defined
     - 0
     - the case's ``bdry<Case>`` routine fills the ghost zones
   * - zero gradient (outflow)
     - 1
     -
   * - reflective
     - 2
     - the cloud case's boundary routine adds a no-inflow clamp here
   * - periodic
     - 3
     - requires ``periods(:) = .true.`` in the MPI topology; required by turbulence driving

Equation of state
-----------------

.. list-table::
   :header-rows: 1
   :widths: 40 20 40

   * - EOS
     - value
     - notes
   * - isothermal, :math:`P=c_s^2\rho`
     - 1
     - ``setSoundSpeed(snd=...)``; cloud default
   * - adiabatic, :math:`P=(\Gamma-1)e_{\rm int}`
     - 2
     - ``setAdiGamma(gam=...)``; dual energy available in 2D/3D Cartesian MHD

Riemann solver
--------------

.. list-table::
   :header-rows: 1
   :widths: 40 20 40

   * - solver
     - value
     - notes
   * - exact (HD)
     - 1
     -
   * - HLL (HD)
     - 2
     -
   * - HLLC (HD)
     - 3
     -
   * - HLL (MHD)
     - 4
     -
   * - HLLD (MHD)
     - 5
     - Miyoshi & Kusano (2005); cloud default. Both isothermal and adiabatic
       variants exist and are chosen from ``eosType``.

Slope limiter and reconstruction
--------------------------------

.. list-table::
   :header-rows: 1
   :widths: 40 20 40

   * - limiter
     - value
     - notes
   * - zero (first order)
     - 0
     -
   * - van Leer
     - 1
     -
   * - fslop
     - 2
     -
   * - minmod
     - 3
     - cloud default

Reconstruction is piecewise-linear (PLM) by default. Piecewise-parabolic
reconstruction (PPM, Colella & Woodward 1984) is selected at run time with
``SCORPIO_PPM=1`` and needs ``nbuf = 3`` ghost cells (with ``nbuf = 2`` the
code falls back to PLM with a one-time warning). The constrained-transport
corner EMF is the Gardiner & Stone (2005) upwinded form by default
(``SCORPIO_UPWIND_EMF=0`` restores the centred average).

.. note::
   Not every combination exists: 1D problems have no MHD CT update, the polar
   coordinates carry no PPM / dual-energy path, and AMR is Cartesian only.


Basic parameters
================

Set in the driver before the grid is built. For a normal user: the number of
dimensions, the mesh, the box, the end time and the sound speed (or
:math:`\Gamma`); the rest can stay at the defaults.

.. list-table::
   :header-rows: 1
   :widths: 24 22 54

   * - variable
     - values
     - notes
   * - ``gridID``
     - integer
     - Problem ID (``setGridID``); selects the initial-condition and boundary
       routines in ``gridModule_init`` / ``gridModule_boundary`` and names the output files.
   * - ``ndim``
     - 1 / 2 / 3
     - Number of dimensions.
   * - ``nbuf``
     - 2 (default), 3 for PPM
     - Ghost cells on every side; ``q`` is indexed ``1-nbuf : nMesh+nbuf``.
   * - ``coordType``
     - 1 / 2 / 3
     - see :ref:`sec:coordtype`.
   * - ``variable(1:8)``
     - 0 / 1
     - Which conserved variables exist (:ref:`sec:conserved`).
       MHD needs 5–7 all 1; ``variable(8)`` must be 1 for any MHD run, isothermal included.
   * - ``nMesh(1:ndim)``
     - integer
     - Global cells per direction. Each direction is split among the MPI ranks
       of that direction, so keep it divisible (powers of two are safest).
   * - ``leftBdry(1:ndim)``, ``rightBdry(1:ndim)``
     - double precision, pc
     - Box edges (cloud default −5..5, −5..5, −10..10).
   * - ``sndspd``
     - double precision, km/s
     - Isothermal sound speed (cloud default 0.3 ≈ 25 K).
   * - ``gam`` / ``adiGamma``
     - 5/3 (default)
     - Adiabatic index; unused when ``eosType = 1``.
   * - ``CFL``
     - 0 < CFL < 1, 0.4 default
     - Courant number for ``griddt``. Use ≤ 0.35 with ``SCORPIO_AD_SCHEME=imexpp``.
   * - ``time_end`` (``tend``)
     - double precision, code time
     - End of the run (cloud: ``tend_tff * tff``).
   * - ``dt_out`` (``dtout``)
     - double precision, code time
     - Snapshot interval (cloud: ``dtout_tff * tff``).
   * - ``tff``
     - 1.5353 (cloud)
     - Free-fall time used only to scale ``tend`` and ``dtout`` in the cloud driver;
       a fixed number, not recomputed from the density.

Run-time state you will see in the log (not user settings): ``t`` current
time, ``dt`` current step (``griddt``, CFL-limited and clipped to hit the next
output time exactly), ``nstep``, ``fnum`` next file number, ``writeFlag``.

.. _sec:grid_api:

Grid setup API (call order)
---------------------------

The drivers all follow the same sequence. The two-phase self-gravity setup
and the restart branch are the parts that are easy to get wrong.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - call
     - what it does / when
   * - ``setRunOption("SCORPIO_...", "value")``
     - Optional; first thing in the driver. Sets a run option from inside the
       code (overrides the shell environment; empty string = no-op).
   * - ``g%setGridID(gridID)``
     - Problem ID.
   * - ``g%enableDrivingTurbulence(DT_mode)``
     - If driving: sets ``enable_DT``, ``DT_mode``, initialises FFTW-MPI. Before ``setMesh``.
   * - ``g%enableSelfgravitySolver(sgSolver, ndim, dims, nprocs)``
       ``g%setSgBdryType(sgBdryType)``, ``g%setGravConst(GravConst)``
     - If gravity: **phase 1**, must precede ``setTopologyMPI``/``setMesh``.
       ``sgSolver`` 0 = FFT, 1 = multigrid.
   * - ``g%setTopologyMPI(ndim, dims, periods, reorder)``
     - Cartesian communicator; ``dims`` from ``MPI_DIMS_CREATE`` or by hand.
   * - ``g%setMesh(nMesh, leftBdry, rightBdry, nbuf, coordType, gridID)``
     - Allocates the local mesh and coordinates. **Not** on a restart (``setTime`` does it).
   * - ``g%setVariable(variable)``
     - Allocates ``q`` (reads ``SCORPIO_DUAL_ENERGY``, ``SCORPIO_PPM`` and the other run options here). Not on a restart.
   * - ``g%setMPIWindows()``
     - One-sided MPI windows for ghost exchange. Not on a restart.
   * - ``g%setEoS(eosType)``, ``g%setSoundSpeed(snd)`` or ``g%setAdiGamma(gam)``,
       ``g%setCFL(CFL)``, ``g%setSlopeLimiter(limiterType)``, ``g%setSolverType(solverType)``,
       ``g%setBoundaryType(boundaryType)``
     - Plain settings; on a restart they are read back from the snapshot instead.
   * - ``g%setTime(fstart, tend, dtout)``
     - ``fstart = 0``: fresh start (``t = 0``, ``fnum = 0``).
       ``fstart = N``: **restart** — opens ``g<gridID>_<N>.h5``, calls
       ``setMesh``/``setVariable``/``setMPIWindows`` itself, restores ``t``,
       the settings above and all fields, exchanges ghost zones.
       Only ``tend`` and ``dtout`` are taken from the call.
   * - ``g%initSelfgravitySolver()``
     - If gravity: **phase 2**, after the mesh exists on either path (FFT
       kernels and windows; multigrid windows). Returns at once if gravity is off.
   * - ``g%initVariable()``
     - Calls the case's initial-condition routine (fresh start only).
   * - ``g%exchangeBdryMPI(g%q, g%winq)``, ``g%setBoundary(g%q)``
     - Fill ghost zones (MPI, then physical boundaries). Repeat after anything that changes interior cells (e.g. a turbulence kick).
   * - ``g%calcSelfgravity(g%q)``
     - Potential/force for the first snapshot; inside the time loop RK2 calls it itself.
   * - ``g%writeGrid()``, ``g%writeGrid_HD_vtk()``
     - HDF5 snapshot (and optional VTK).
   * - loop: ``g%griddt()``, ``g%evolveGridRK2()``, ``g%writeGrid()`` when ``g%writeFlag``
     - Time step, second-order RK update (with gravity, dual energy, FOFC failsafes inside), output.
   * - ``g%TrueloveCondition(g%q, stop_job)``
     - Optional isothermal collapse stop: ``stop_job = .true.`` when the Jeans length at the global maximum density is below 5 cells.
   * - ``g%enableAD(.true.)``, ``g%setADparams(mu_ad, alpha_ad)``
     - Two-fluid cases only, on both the neutral and the ion grid.


File output and restart
=======================

.. list-table::
   :header-rows: 1
   :widths: 24 22 54

   * - variable
     - values
     - notes
   * - ``dt_out``
     - double precision, code time
     - Snapshot interval; the time step is clipped so that snapshots land exactly on multiples of it.
   * - ``toutput``
     - (state)
     - Time of the next snapshot.
   * - ``fnum``
     - (state)
     - Number of the next snapshot; file name ``g<gridID>_<fnum>.h5`` (4 digits each, e.g. ``g0800_0042.h5``).
   * - ``write_vtk``
     - ``.true.`` / ``.false.``
     - Also write VTK/PVTS files (``writeGrid_HD_vtk`` / ``writeGrid_MHD_vtk``) for ParaView.
   * - ``Restart`` (cloud driver)
     - ``.true.`` / ``.false.``
     - The **flag** that selects a restart. ``.false.``: build the grid, apply the initial condition, write ``#0000``.
       ``.true.``: ``setTime(fstart=file_start)`` restores everything from the snapshot; no initial condition and no repeat of the t = 0 turbulence kick
       (``DT_mode = 0``: the run continues undriven; ``DT_mode = 1``: the kick schedule resumes from the ``.turb`` side file, see turbulence driving).
       Namelist equivalent ``cloud_restart``.
   * - ``file_start``
     - integer ≥ 0
     - Snapshot number to continue from; only meaningful with ``Restart = .true.``
       (namelist ``cloud_fstart``). The driver aborts with a message if the flag and the number disagree or the file is missing.
   * - AMR checkpoints
     - ``c<caseID>_<fnum>.h5``
     - Written next to every AMR output unless ``SCORPIO_AMR_CHECKPOINT=0``; restart with ``&amr_restart amr_fstart = N``.

A restarted run reproduces the uninterrupted run bit for bit (verified for
the cloud case with the same and with a different number of MPI ranks); the
rank count may change between the original run and the restart because the
reader builds each rank's hyperslab from the current decomposition.

Datasets in a snapshot
----------------------

Every ``g<gridID>_<fnum>.h5`` is self-describing:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - dataset
     - content
   * - ``nMesh``, ``nbuf``, ``ndim``, ``nvar``, ``coordType``, ``gridID``, ``variable``
     - Mesh size, ghost width, dimensions, number of variables, coordinate system, case ID, the ``variable(1:8)`` flags.
   * - ``leftBdry``, ``rightBdry``
     - Box edges [pc].
   * - ``t``, ``CFL``, ``snd``, ``adiGamma``, ``eosType``, ``solverType``, ``limiterType``, ``boundaryType``
     - Time [code] and the numerical settings — this is what a restart reads back.
   * - ``xl1``, ``xr1``, ``xc1``, ``dx1`` (and ``…2``, ``…3``)
     - Left face, right face, centre and width of every cell along x (y, z), **including** the ``nbuf`` ghost cells at each end.
   * - ``den``, ``momx``, ``momy``, ``momz``
     - Density and momentum density (:ref:`ch:hydro` for units); ghost cells included; 3D arrays appear as ``(z, y, x)`` in ``h5py``.
   * - ``bxl``, ``byl``, ``bzl``, ``bxr``, ``byr``, ``bzr``
     - Face-centred magnetic field (left/right faces); cell value = average of the pair.
   * - ``ene``
     - Total energy density (adiabatic runs; carried but unused for isothermal).
   * - ``gphi``, ``sgfx``, ``sgfy``, ``sgfz``
     - Gravitational potential and acceleration when self-gravity is on (written from the last RK2 stage that evaluated gravity, so consistent with the state to :math:`O(\Delta t)`).


MPI setup
=========

.. note::
   If the periodic boundary condition ``boundaryType = 3`` is used, the MPI
   periodic flags ``periods(1)/(2)/(3)`` must be ``.true.``. The remaining
   parameters can be left as default.

.. list-table::
   :header-rows: 1
   :widths: 24 22 54

   * - variable
     - values
     - notes
   * - ``periods(1:ndim)``
     - ``.true.`` / ``.false.``
     - Periodicity of the Cartesian communicator per direction.
   * - ``dims(1:ndim)``
     - integer, 0 = automatic
     - Ranks per direction. ``dims = 0`` then ``MPI_DIMS_CREATE(nprocs, ndim, dims, ierr)`` lets MPI choose.
       The FFT self-gravity and turbulence solvers no longer require a slab layout (``dims = (1,1,nprocs)``):
       they remap to FFTW's slab decomposition internally with ``MPI_ALLTOALLV``.
   * - ``reorder``
     - ``.true.`` / ``.false.``
     - MPI may renumber ranks in the Cartesian communicator.
   * - ``myid``, ``nprocs``
     - (module variables)
     - Rank and size, set by ``setMPI(np)`` in ``main.f03``; use ``if (myid .eq. 0)`` for printing.
   * - ``ierr``, ``errcode``
     - integer
     - MPI status arguments.

Each direction of ``nMesh`` must be divisible by the number of ranks assigned
to it; with ``MPI_DIMS_CREATE`` and rank counts 1, 2, 4, 8, 16 … this is
automatic for meshes that are multiples of 16.


Physics modules
===============

Turbulence driving
------------------

Driving injects a divergence-free / compressive velocity perturbation with a
prescribed power spectrum (Federrath et al. 2010 projection,
:math:`P_{ij} = \zeta\,\delta_{ij} + (1-2\zeta)\,k_i k_j/k^2`) and rescales it so
the injected kinetic energy equals the requested value.

.. list-table::
   :header-rows: 1
   :widths: 26 20 54

   * - variable
     - values
     - notes
   * - ``DriveTurbulence``
     - ``.true.`` / ``.false.``
     - Master switch in the driver (``enableDrivingTurbulence``); cloud default on; namelist ``cloud_driving``.
   * - ``DT_mode``
     - 0 / 1
     - 0: one kick at :math:`t=0`. 1: kick at :math:`t=0` and then every ``dt_turb`` until ``n_turb`` kicks are done.
   * - ``E_turb_tot``
     - double precision, :math:`M_\odot` km\ :sup:`2` s\ :sup:`-2`
     - Total kinetic energy injected over all kicks (cloud default 4.5 = :math:`9.0\times10^{43}` erg).
   * - ``n_turb``
     - integer
     - Number of kicks; each injects ``E_turb = E_turb_tot / n_turb`` (= ``Energy_DT``).
   * - ``turn_over_time``, ``dt_turb``
     - double precision, code time
     - ``dt_turb = turn_over_time / n_turb`` is the interval between kicks in mode 1.
   * - ``t_count_turb``, ``t_accum_turb``
     - integer / double precision
     - Kicks done so far and time since the last kick. Not part of the HDF5 snapshot: a ``DT_mode = 1`` run writes them to ``g<gridID>_<fnum>.turb`` next to every snapshot,
       and a restart reads that file so the next kick fires on the same step as in an uninterrupted run (if the file is missing they are estimated from ``t`` with a warning).
       The kick phases are drawn from the OS-seeded random generator, so the realisation after a restart differs from the uninterrupted run.
   * - ``zeta`` → ``zeta_DT``
     - 0 ≤ ζ ≤ 1
     - 1 = purely solenoidal, 0 = purely compressive, 0.5 = natural mix.
   * - ``DT_scale`` → ``drivingWN_DT``
     - double precision
     - Driving wavenumber :math:`k_d` in units of :math:`2\pi/L`; 2 = eddies of half the box.
   * - ``Energy_DT``, ``DTenergyfaction``
     - double precision
     - Energy per kick, and the fraction of it applied by the next call
       (``setDTenergyfaction(t_accum_turb, dt_turb)`` sets the fraction to ``t_accum_turb/dt_turb``).
   * - ``netmomx_DT``, ``netmomy_DT``, ``netmomz_DT``
     - double precision
     - Net momentum removed from the perturbation (normally 0).
   * - ``SCORPIO_DRIVING_SPECTRUM`` / ``OPT_DRIVING_SPECTRUM``
     - ``expo`` (default), ``kolmogorov``, ``burgers``
     - Shape of the driving spectrum per mode :math:`k`:
       ``expo``: :math:`k^6 e^{-8k/k_d}` (Otto et al. 2017);
       ``kolmogorov``: :math:`k^{-11/3}` for :math:`k\ge2` (shell :math:`E(k)\propto k^{-5/3}`);
       ``burgers``: :math:`k^{-4}` for :math:`k\ge2` (shell :math:`E(k)\propto k^{-2}`). The cloud driver sets ``burgers`` in code.

.. note::
   Only the periodic boundary condition is supported for turbulence driving;
   the MPI periodic flags ``periods(1)/(2)/(3)`` must be ``.true.``. After a
   kick, call ``exchangeBdryMPI`` and ``setBoundary`` again — the kick fills
   interior cells only.

Self-gravity
------------

.. list-table::
   :header-rows: 1
   :widths: 26 20 54

   * - variable
     - values
     - notes
   * - ``SelfGravity``
     - ``.true.`` / ``.false.``
     - Master switch in the driver; gates ``enableSelfgravitySolver`` (phase 1), ``initSelfgravitySolver`` (phase 2) and ``calcSelfgravity``.
   * - ``sgSolverType`` (``sg_solver``)
     - 0 / 1
     - 0 = FFTW convolution (Green-function kernels for isolated, spectral Poisson for periodic); 1 = multigrid. Argument of ``enableSelfgravitySolver``; namelist ``cloud_sgsolver``. AMR is multigrid only.
   * - ``sgBdryType`` (``sg_bdry``)
     - 0 / 1
     - 0 = isolated (open) boundaries, 1 = periodic Poisson; ``setSgBdryType``; namelist ``cloud_sgbdry``.
   * - ``GravConst``
     - 4.3011e-3
     - :math:`G` in code units, ``setGravConst``; do not change.
   * - ``gphi``
     - output
     - Potential :math:`\Phi` [km\ :sup:`2` s\ :sup:`-2`]; periodic solutions are mean-zero.
   * - ``sgfx``, ``sgfy``, ``sgfz``
     - output
     - :math:`\boldsymbol{g}=-\nabla\Phi`; the momentum source is :math:`\rho\boldsymbol{g}` and the energy source :math:`\rho\boldsymbol{u}\cdot\boldsymbol{g}`.
   * - ``SCORPIO_MG_ISOBC``
     - unset / ``multipole``
     - Isolated multigrid boundary: James (1977) screening-charge boundary by default; ``multipole`` selects the multipole Dirichlet boundary.
   * - ``do_truelove`` (``cloud_truelove``)
     - ``.true.`` / ``.false.``
     - Cloud driver: stop (after a final snapshot) when :math:`\lambda_J < 5\,\Delta x` at the densest cell.

.. note::
   Periodic self-gravity ``sgBdryType = 1`` needs ``periods(1)/(2)/(3) = .true.``.
   Phase 1 (``enableSelfgravitySolver``) must be called *before*
   ``setTopologyMPI``/``setMesh`` and phase 2 (``initSelfgravitySolver``)
   *after* the mesh exists — on a restart that means after ``setTime``.

Ambipolar diffusion (two-fluid)
-------------------------------

.. list-table::
   :header-rows: 1
   :widths: 26 20 54

   * - variable
     - values
     - notes
   * - ``enable_ad``
     - ``.true.`` / ``.false.``
     - ``enableAD`` on both grids.
   * - ``mu_ad``
     - double precision, amu
     - Mean molecular weight of the species (2.3 neutrals, 29 ions in the cloud setups); ``setADparams``.
   * - ``alpha_ad``
     - double precision, pc\ :sup:`2` km s\ :sup:`-1` :math:`M_\odot^{-1}`
     - Collision coefficient :math:`\alpha` in code units; :math:`3.7\times10^{13}` cm\ :sup:`3` g\ :sup:`-1` s\ :sup:`-1` = :math:`7.7\times10^{4}` code units (:ref:`ch:hydro`).
   * - ``gridIDn``, ``gridIDi``
     - integer pair
     - Neutral and ion grid IDs of a two-fluid case (e.g. 609/616, 620/621).
   * - ``SCORPIO_AD_SCHEME``
     - ``split`` (default), ``imex``, ``imex322`` | ``athenak``, ``imexpp`` | ``imex2+`` | ``krapp``
     - Drag integration: operator-split TR-BDF2 (default); IMEX-SSP2(2,2,2); IMEX-SSP2(3,2,2) as in AthenaK; IMEX(4,3,2) of Krapp et al. (2024, CFL ≤ 0.35). An unrecognised value aborts.

Dual energy and robustness options
----------------------------------

These are read once, in ``setVariable``, from the environment (or
``setRunOption``). Defaults are the production settings. What each
method does is described in :ref:`ch:methods`.

.. list-table::
   :header-rows: 1
   :widths: 30 18 52

   * - variable
     - default
     - effect
   * - ``SCORPIO_DUAL_ENERGY``
     - on
     - ``0`` disables the entropy-based pressure recovery in low-:math:`\beta` cells (pure total-energy update).
   * - ``SCORPIO_DE_THR``
     - 1e-3
     - Dual energy fires when :math:`e_{\rm int} < {\rm thr}\times E`.
   * - ``SCORPIO_DE_PRINT``
     - on
     - ``0`` silences the immediate per-step report of cells rescued from negative pressure.
   * - ``SCORPIO_PPM``
     - off
     - ``1`` = piecewise-parabolic reconstruction (needs ``nbuf = 3``).
   * - ``SCORPIO_UPWIND_EMF``
     - on
     - ``0`` = centred corner EMF instead of the Gardiner & Stone (2005) upwinded CT EMF.
   * - ``SCORPIO_LEGACY_FAILSAFE``
     - off
     - ``1`` = whole-domain HLLD→HLL switch / global dt halving instead of the localized first-order flux correction (FOFC).
   * - ``SCORPIO_DT_HALVE``
     - on
     - ``0`` removes the last rung of the failsafe ladder (dual energy → FOFC → dt halving → abort): after FOFC fails the run aborts immediately.
   * - ``SCORPIO_HLLD_EPS``
     - 1e-8
     - Relative degeneracy threshold inside the HLLD star states (0 = legacy exact-equality tests).
   * - ``SCORPIO_HLLD_STARCHECK``
     - off
     - ``1`` = per-face admissibility check of the HLLD fan with HLL fallback on that face.
   * - ``SCORPIO_HEALTH``
     - off
     - ``1`` = per-step ``[health]`` line (rescued cells, FOFC faces, HLLD fallbacks); costs one reduction per step.
   * - ``SCORPIO_DRIVING_SPECTRUM``
     - ``expo``
     - See turbulence driving.
   * - ``SCORPIO_AD_SCHEME``
     - ``split``
     - See ambipolar diffusion.
   * - ``SCORPIO_MG_ISOBC``
     - James
     - See self-gravity.

The cloud driver fixes three of these in code through ``setRunOption`` at
its top (``OPT_DRIVING_SPECTRUM = 'burgers'``, ``OPT_DUAL_ENERGY``,
``OPT_PPM`` together with ``nbuf_opt``), so the shell environment cannot
change them by accident; an empty string keeps the environment value.


Case-specific namelist groups
=============================

The 20 pc cloud (``gridID`` 800 / 801): ``&cloud_nml``
-------------------------------------------------------

Every entry is optional; a missing entry keeps the value in the driver.

.. list-table::
   :header-rows: 1
   :widths: 24 22 54

   * - entry
     - default (driver)
     - notes
   * - ``cloud_nmesh``
     - 32, 32, 64
     - ``nMesh(1:3)``; 801: the base grid.
   * - ``cloud_tend_tff``
     - 10.0
     - End time in units of ``tff`` (= 1.5353 code time).
   * - ``cloud_dtout_tff``
     - 0.1
     - Snapshot interval in units of ``tff``.
   * - ``cloud_truelove``
     - ``.true.``
     - Truelove stop (800 only; 801 uses ``amr_jeans_n``).
   * - ``cloud_sgsolver``
     - 0
     - 0 FFT / 1 multigrid (800 only; 801 is multigrid).
   * - ``cloud_sgbdry``
     - 0
     - 0 isolated / 1 periodic gravity.
   * - ``cloud_driving``
     - ``.true.``
     - Turbulence kick (800 only; not implemented on the AMR path).
   * - ``cloud_restart``
     - ``.false.``
     - Continue from a snapshot (800 only).
   * - ``cloud_fstart``
     - 0
     - Snapshot number for ``cloud_restart = .true.``.

.. code-block:: fortran

   &problem_config
     gridID = 800
   /
   &cloud_nml
     cloud_nmesh = 64, 64, 128
     cloud_tend_tff = 3.0
     cloud_dtout_tff = 0.05
     cloud_driving = .false.
   /

AMR (``gridID`` 700–727, 801)
-----------------------------

The block-structured AMR driver (``src/amrModule.f03``, design in
``docs/AMR_DESIGN.md``) reads these groups. Refinement uses the Löhner (1987)
second-derivative indicator on density (plus total energy and :math:`B^2` in
3D), optionally combined with a Jeans-length criterion.

.. list-table::
   :header-rows: 1
   :widths: 18 24 18 40

   * - group
     - entry
     - default
     - notes
   * - ``&amr_config``
     - ``amr_bsx``, ``amr_bsy``
     - case-specific (16)
     - Block size in cells along x, y.
   * -
     - ``amr_max_level``
     - case-specific
     - 0 = base grid only; each level halves :math:`\Delta x`.
   * -
     - ``amr_regrid_every``
     - 0
     - 0 = static tree; N = re-evaluate the refinement every N steps.
   * -
     - ``amr_refine_thr``, ``amr_deref_thr``
     - 0.8, 0.2
     - Löhner indicator thresholds for refining / derefining a block.
   * - ``&amr_config_z``
     - ``amr_bsz``
     - 16
     - Block size along z (3D).
   * - ``&amr_tuning``
     - ``amr_jeans_n``
     - 0 (off)
     - Refine wherever :math:`\lambda_J < ` ``amr_jeans_n`` :math:`\Delta x` (Truelove et al. 1997).
   * -
     - ``amr_deref_hyst``
     - 2
     - Consecutive derefine votes required before a block is coarsened.
   * -
     - ``amr_flag_buffer``
     - 0
     - Dilate refine flags by N block shells.
   * -
     - ``amr_fac``, ``amr_fac_iters``
     - 0, 1
     - Enable the FAC multi-level gravity solve (``docs/FAC_GRAVITY_DESIGN.md``) and its number of composite passes.
   * - ``&amr_restart``
     - ``amr_fstart``
     - 0
     - Restart from checkpoint ``c<caseID>_<N>.h5``.
   * - ``&amr_tend_nml``
     - ``amr_tend``
     - case-specific
     - Override the end time.
   * - ``&amr_blob_nml``, ``&amr_refbox_nml``
     - ``amr_blob_c``, ``amr_refbox``
     - test-specific
     - Parameters of individual AMR test cases.

AMR-specific environment variables: ``SCORPIO_AMR_CHECKPOINT=0`` (no
checkpoint files), ``SCORPIO_AMR_CFL`` (override the CFL number),
``SCORPIO_FAC_SOLVER=sor`` (pointwise SOR instead of multigrid inside FAC),
``SCORPIO_AMR_AUDIT=1`` (regrid audit output).

Test-suite knobs
----------------

Individual test drivers read further ``SCORPIO_*`` variables so that
parameter scans can be run without editing code. They only affect the case
named:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - variable(s)
     - case
   * - ``SCORPIO_SG_SOLVER``, ``SCORPIO_SPHERE_RADIUS``, ``SCORPIO_SPHERE_OFFSET``
     - self-gravity tests (357–361): solver choice, sphere geometry
   * - ``SCORPIO_WAVE_N``, ``SCORPIO_WAVE_LIMITER``, ``SCORPIO_WAVE_RHO``, ``SCORPIO_WAVE_B0``, ``SCORPIO_WAVE_GAM``, ``SCORPIO_WAVE_P``, ``SCORPIO_WAVE_AMP``, ``SCORPIO_SOLVER``
     - circularly polarised Alfvén wave tests (18, 58, 60): resolution, limiter, background state, amplitude, solver
   * - ``SCORPIO_TURB_N``, ``SCORPIO_TURB_TEND``, ``SCORPIO_TURB_DTOUT``, ``SCORPIO_COUPLE_T``
     - two-fluid driven turbulence (609/616): resolution, end time, output cadence, ion coupling time
   * - ``SCORPIO_ADW_ALPHA``, ``SCORPIO_ADW_N``, ``SCORPIO_ADW_TEND``
     - two-fluid Alfvén damping (620/621)
   * - ``SCORPIO_CSK_NY``, ``SCORPIO_CSK_NZ``, ``SCORPIO_CSK_SOLVER``, ``SCORPIO_CSK_CFL``, ``SCORPIO_CSK_TEND``, ``SCORPIO_CSK_SMOOTH``, ``SCORPIO_CSK_ATHENAK``, ``SCORPIO_CSK_ATHDIR``
     - 3D C-shock (52/53), including the AthenaK comparison hooks
   * - ``SCORPIO_SPEC_NMESH``, ``SCORPIO_SPEC_TEND``, ``SCORPIO_SPEC_TPHASE2``
     - spectrum-compensated driven turbulence (401)
   * - ``SCORPIO_AMR_NMESH``, ``SCORPIO_AMR_TESTG``, ``SCORPIO_AMR_TESTBG``
     - AMR self-gravity tests (720/721)
   * - ``SCORPIO_NBUF``, ``SCORPIO_TEST_LIMITER``, ``SCORPIO_RAW``
     - 3D isothermal MHD shock tube (48), field-loop reconstruction study (61), raw (no failsafe) diagnostic mode of the 3D ion blast (59)
