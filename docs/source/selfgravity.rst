.. _ch:selfgravity:

************
Self-gravity
************

Introduction
============

The potential is obtained from the Poisson equation
:math:`\nabla^2\Phi = 4\pi G\rho`, and the acceleration
:math:`\boldsymbol{g} = -\nabla\Phi` enters the momentum equation as
:math:`\rho\boldsymbol{g}` and the energy equation as
:math:`\rho\boldsymbol{u}\cdot\boldsymbol{g}`. Two solvers exist behind one
interface — an FFT convolution/spectral solver and a multigrid solver — and
both support isolated (open) and periodic boundaries. The solver is chosen
with ``sgSolverType`` (0 = FFT, 1 = multigrid) and the boundary with
``sgBdryType`` (0 = isolated, 1 = periodic); :math:`G` is ``GravConst``
(4.3011e-3 in code units, :ref:`ch:hydro`).

Where gravity enters the time step
==================================

.. code-block:: text

   rk2_3D  (gridModule_rk2.f03)             — once per RK2 stage
   ├─ sweeps of the stage
   ├─ if (this%enable_sg)  call this%calcSelfgravity(<stage input state>)
   │      q_out(momx) += ρ sgfx Δt,  … momy, momz
   │      q_out(ene)  += (momx sgfx + momy sgfy + momz sgfz) Δt      (adiabatic only)
   └─ ghost exchange, boundaries

   dt3D  (gridModule_dt.f03)                — the acceleration also limits Δt:
          Δt ≤ 0.2 ( −|u|/g + sqrt(|u|²/g² + 2 min(Δx)/g) )

So gravity is recomputed from the stage's input state at every stage
(twice per step); the snapshot fields ``gphi``, ``sgfx/y/z`` are whatever the
last such call left behind.

Setup: the two-phase API
========================

.. code-block:: text

   driver
   ├─ g%enableSelfgravitySolver(sgSolver, ndim, dims, nprocs)     phase 1 — BEFORE setTopologyMPI/setMesh
   │     sgSolver = 0 → enableSelfgravity        (enable_sg, fftw_mpi_init, FFT-slab layout)
   │     sgSolver = 1 → enableSelfgravity_MG     (enable_sg, sgSolverType = MG)
   ├─ g%setSgBdryType(0 | 1),  g%setGravConst(G)
   ├─ g%setTopologyMPI … g%setMesh … g%setVariable                 (setVariable creates the FFT plans: sgPlan3D)
   ├─ g%initSelfgravitySolver()                                    phase 2 — AFTER the mesh exists
   │     FFT: setSelfgravityKernel → sgkernel3d (Green kernels, transformed once)
   │          setSGMPIWindows      → initSGWindows3D
   │     MG:  setSGMPIWindows_MG
   └─ g%calcSelfgravity(g%q)                                       first force, for snapshot #0

On a restart ``setTime(fstart>0)`` does the ``setMesh``/``setVariable``
part, and phase 2 comes after it. With ``SelfGravity = .false.`` none of this
is called and ``initSelfgravitySolver`` returns at once if it is.

FFT solver
==========

.. code-block:: text

   calcSelfgravity(q)                                   gridModule.f03 — dispatch on sgSolverType, sgBdryType, ndim
   └─ calcSG3D(q, …)             isolated             gridModule_sgcalc.f03
      ├─ remap_density_hydro_to_fft3d     hydro decomposition → FFTW slab layout, zero-padded box (MPI_ALLTOALLV)
      ├─ fftw_mpi_execute_dft (forward)   ρ̂(k)
      ├─ multiply by the three force kernels and the potential kernel (sgfxKernel …, sgPhiKernel)
      ├─ fftw_mpi_execute_dft (inverse) ×4
      ├─ remap_force_fft_to_hydro3d       g_x, g_y, g_z back to the hydro layout
      └─ remap_scalar_fft_to_hydro3d      Φ back to the hydro layout (gphi)
   └─ calcSG3Dperiodic(q, …)     periodic             gridModule_sgcalc.f03
      same remaps; spectral solve Φ̂ = −4πG ρ̂ / k², mean-zero, forces from ik Φ̂

The isolated solver is the zero-padded Green-function convolution (exact
free-space boundary): the kernels are the FFTs of the :math:`-1/r`
potential kernel (zero self-cell contribution) and of the force kernels,
computed once in ``sgkernel3d``. Because FFTW-MPI transforms a slab decomposition, the hydro
decomposition (any ``dims``) is remapped to slabs and back with packed
rectangular ``MPI_ALLTOALLV`` exchanges; the old requirement
``dims = (1, 1, nprocs)`` no longer exists.

Multigrid solver
================

.. code-block:: text

   calcSelfgravity(q) → calcSG3D_MG(q, gphi, sgfx, sgfy, sgfz)        calcSG_MG.f03
   ├─ build the level hierarchy: halve the local mesh while it stays even and ≥ 4 cells; the
   │  level count is the minimum over ranks and every level is collective; rhs = 4πG ρ
   ├─ periodic:  subtract the global mean density; project_level_mean keeps every level mean-free
   ├─ isolated:  boundary values for the ghost layer of the finest level
   │     james (default): pass 1 solves with zero Dirichlet ghosts; compute_james_boundary reads the
   │                      screening charge off the boundary layer of that solution, sums its
   │                      potential onto the ghost positions (psi_sum); pass 2 solves with those values
   │     multipole (SCORPIO_MG_ISOBC=multipole): compute_multipole_moments → monopole, dipole,
   │                      quadrupole Dirichlet values
   ├─ run_cycles → v_cycle (recursive): smooth (weighted Jacobi, ω = 0.8, 3 pre / 3 post sweeps,
   │              60 on the coarsest) → compute_residual → restrict_residual → coarse solve →
   │              prolong_add; apply_mg_boundary and the level ghost exchange (exchgBdryMPI_sgMG)
   │              at every smoothing sweep; stop at relative residual 1e-6 (max 200 cycles)
   └─ forces by centred differences of Φ; periodic Φ is written in the zero-mean gauge

References: James (1977) for the screening-charge boundary (cf. Ricker
2008; Moon, Kim & Ostriker 2019); the periodic sinusoid gate follows Tomida
& Stone (2023). The 2D solver (``calcSG2D_MG``) uses a complex multipole
expansion to order 4 for the isolated boundary.

On the AMR mesh
===============

AMR runs use the multigrid solver only. The base-level covering grid
``gio`` receives the block densities, solves, and scatters the forces:

.. code-block:: text

   amrStep3D → amrSolveGravity3(m, sel)                            amrModule.f03
   ├─ amrGatherDensity3      restrict every leaf block's density onto gio
   ├─ m%gio%calcSelfgravity  → calcSG3D_MG on the base grid
   ├─ amrScatterForce3       forces (and Φ) back to the blocks
   └─ [amr_fac = 1]  facRefineForces3   nested level solves with Dirichlet data prolonged from below
                     (block-local multigrid, red-black block colouring; SCORPIO_FAC_SOLVER=sor for SOR);
                     amr_fac_iters > 1: facTauCorrect3 → base re-solve → scatter → refine again
   amrApplyGravSource3(m, sel, tgt)   momentum / energy sources per block

Design and gate results are in ``docs/FAC_GRAVITY_DESIGN.md``; the summary
is in :ref:`ch:methods`.

Tests
=====

.. list-table::
   :header-rows: 1
   :widths: 16 84

   * - ``gridID``
     - case
   * - 21, 35
     - 2D / 3D FFT self-gravity tests (periodic Poisson path)
   * - 49
     - 3D isothermal MHD with FFT gravity
   * - 357, 358
     - 2D / 3D multigrid, isolated (multipole / James boundary), uniform sphere
   * - 557, 558
     - 2D / 3D multigrid, periodic
   * - 359
     - Maclaurin spheroid against the analytic potential (Chandrasekhar 1969; Ricker 2008)
   * - 360
     - pressureless free-fall collapse against the cycloid solution
   * - 361
     - periodic sinusoidal density against its exact potential
   * - 720, 721, 722, 724–727
     - AMR gravity: 2D/3D uniform sphere, Jeans collapse, Jeans wave, two spheres, Evrard collapse, periodic sinusoid

``validation/gravity_analytic.sh`` runs the analytic gates and reports the
force and potential errors; ``validation/gravity_deep.sh`` runs the deep
collapse ladder.
