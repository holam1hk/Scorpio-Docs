.. _ch:parallelization:

***************
Parallelization
***************

.. note::
   If you are not familiar with MPI, the `MPI tutorial <https://mpitutorial.com/>`_
   is a good start. Nothing on this page is needed to *run* the code; it
   explains what happens under the hood.

Introduction
============

Scorpio is a pure-MPI code (flat MPI, no OpenMP). It has been built and run
with MPICH and OpenMPI; the parallel HDF5 and FFTW-MPI libraries it links
must be built against the same MPI. The domain decomposition is a Cartesian
block decomposition in all directions; each rank owns one block of the
uniform mesh (or a set of AMR blocks) with ``nbuf`` ghost cells on every side.

Domain decomposition
====================

.. code-block:: text

   driver:  dims = 0;  MPI_DIMS_CREATE(nprocs, ndim, dims)    (or dims set by hand)
   g%setTopologyMPI(ndim, dims, periods, reorder)                gridModule.f03
   ├─ MPI_CART_CREATE   → the Cartesian communicator g%vu_mpi
   ├─ MPI_CART_COORDS   → g%mpiCoord (this rank's position)
   └─ MPI_CART_SHIFT    → g%left_mpi/right_mpi (x), up_mpi/down_mpi (y), top_mpi/bottom_mpi (z)
                          (MPI_PROC_NULL at a non-periodic physical boundary)
   g%setMesh(nMesh, …)   splits nMesh(d) over dims(d): nx_list/ny_list/nz_list hold every rank's
                          local size, meshStart the offset of this rank's block in the global mesh

Each direction of ``nMesh`` must be divisible by the number of ranks along
it. ``periods(d) = .true.`` makes the communicator periodic in that
direction, which is what the periodic hydro boundary (``boundaryType = 3``),
turbulence driving and periodic gravity require.

Ghost-zone exchange
===================

The three states ``q``, ``q1``, ``q2`` are exposed as one-sided MPI windows,
and the exchange is done with fences and remote gets/puts:

.. code-block:: text

   g%setMPIWindows() → initMPIWindows3D(this, q, q1, q2, databuf1, databuf2)     gridModule_windows.f03
       MPI_WIN_CREATE on each state (winq, winq1, winq2) and on the two transfer buffers

   g%exchangeBdryMPI(q, winq) → exchgBdryMPI3D                                    gridModule_exchange.f03
       for each direction:  MPI_WIN_FENCE → MPI_GET / MPI_PUT of the nbuf-wide strips
                            from/to the neighbour rank → MPI_WIN_FENCE
       (non-contiguous strips are packed through databuf1/databuf2 first)

   g%setBoundary(q) → setBdry3D → bdry<Case>(this, q)                              gridModule_boundary.f03
       physical boundaries on the faces whose neighbour is MPI_PROC_NULL

The sequence *exchange, then physical boundary* appears after every stage
of the RK2 step, after a turbulence kick, and after the initial condition.
Ghost widths are generic (``nbuf = 2`` or 3); an old hard-coded ``nbuf = 2``
in the RMA strips was fixed when PPM was added.

Global reductions
=================

- the time step: ``MPI_ALLREDUCE(…, MPI_MIN)`` in ``dt3D``;
- the failsafe flags ``changeSolver`` and ``neg_pressure``:
  ``MPI_ALLREDUCE(…, MPI_LOR)`` in ``rk2_3D`` — so all ranks redo a step
  together;
- the Truelove criterion and the periodic-gravity mean density:
  ``MPI_ALLREDUCE`` of the maximum density / the density sum;
- the optional health monitor (``SCORPIO_HEALTH=1``): one small
  ``MPI_ALLREDUCE`` per step; the default dual-energy print is rank-local
  and needs no communication.

FFT solvers: slab remap
=======================

FFTW-MPI transforms a **slab** decomposition (the last dimension split over
ranks). The hydro decomposition is not required to be a slab any more: the
self-gravity and turbulence-driving solvers remap between the two layouts
with packed rectangular ``MPI_ALLTOALLV`` exchanges, using only small
ownership metadata (no rank replicates a full global array):

.. code-block:: text

   calcSG3D / calcSG3Dperiodic          gridModule_sgcalc.f03
     remap_density_hydro_to_fft3d  → FFTW (sgPlanDen, sgPlanFx/Fy/Fz) → remap_force_fft_to_hydro3d,
                                                                       remap_scalar_fft_to_hydro3d
   calcDT3D_MD                          gridModule_dtcalc.f03
     create_phase_space_profile_3D → FFTW (DTV*_Plan_C2R) → remap_dt_vector_fft_to_hydro3d

The plans are created once in ``setVariable`` (``sgPlan3D``, ``DTPlan3D`` in
``gridModule_fft_plans.f03``); the FFT-side local sizes come from
``fftw_mpi_local_size``. The multigrid gravity solver does not use FFTs; its
levels follow the hydro decomposition and exchange their ghost layers with
``exchangeBdryMPI_sgMG``.

Parallel I/O
============

Snapshots are single HDF5 files written collectively: ``output3d``
(``gridModule_io3d.f03``) opens the file with the MPI-IO driver
(``h5pset_fapl_mpio_f`` on ``vu_mpi``), and every rank writes its own block
as a hyperslab (``h5sselect_hyperslab_f`` with ``meshStart`` as offset,
``H5FD_MPIO_COLLECTIVE_F``). ``read3d`` mirrors this, which is why a restart
may use a different number of ranks. Scalars and coordinate axes are
written by ``output1d``. The VTK/PVTS writers (``gridModule_pvts.f03``)
write one piece per rank plus a ``.pvts`` index.

AMR
===

On the adaptive mesh the tree metadata (levels, positions, neighbours,
owner ranks) is replicated on every rank, so neighbours are found without
communication; the block data are distributed in Morton order and
re-partitioned at each regrid (``regridMoveData3``). Ghost fill
(``amrExchangeLevel3D``), refluxing (``amrReflux3``) and EMF matching use
non-blocking point-to-point messages (``MPI_ISEND`` / ``MPI_IRECV`` /
``MPI_WAITALL``) following exchange plans built once per mesh
(``buildXchgPlan3D``, ``buildRefluxPlan3``, ``buildIoPlan3D``). Output
gathers the leaf blocks onto the base-level grid (``amrScatterToIo3D``) and
then uses the same collective HDF5 writer; checkpoints are written and read
with per-block hyperslabs (``amrWriteCheckpoint3``, ``amrReadCheckpoint3``).
The AMR gate battery checks ``np = 1`` against ``np = 4`` bitwise.

Reproducibility
===============

Results have been verified to be independent of the rank count to the last
bit for the cases in the validation gates and for the cloud case (shared
faces are computed from identical stencil data on both sides; the FFT remap
and the reductions are deterministic). The one exception is
the random turbulence kick, whose Fortran-intrinsic random numbers are
seeded per process from the OS (:ref:`ch:turbulence`).

Running
=======

.. code-block:: bash

   mpirun -n 8 ./Scorpio            # in the directory that holds problem.nml
   MPIRUN=srun ./validation/validate.sh   # the validation scripts honour MPIRUN

Use 1, 2, 4, 8, 16, … ranks and mesh sizes that are multiples of 16; on a
Windows/WSL mount export ``HDF5_USE_FILE_LOCKING=FALSE``. Hybrid MPI+OpenMP
is not implemented; on nodes with many cores flat MPI is the intended mode.
