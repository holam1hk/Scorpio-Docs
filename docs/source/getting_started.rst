.. _ch:getting_started:

***************
Getting Started
***************

Requirements
============

Scorpio needs three libraries, all built against the **same** MPI:

- an MPI library (MPICH or OpenMPI) with a Fortran compiler (gfortran);
- **parallel HDF5** with the Fortran interface — the compiler wrapper
  ``h5pfc`` it installs is what the ``Makefile`` calls, and it supplies the MPI
  and HDF5 flags automatically;
- **FFTW3 with MPI support** (``libfftw3_mpi``).

Python 3 with ``numpy`` and ``h5py`` (and ``matplotlib``) is needed for the
validation scripts and for looking at the output.

On Ubuntu / WSL with MPICH:

.. code-block:: bash

   sudo apt install build-essential gfortran make mpich libhdf5-mpich-dev hdf5-tools \
                    libfftw3-dev python3-numpy python3-h5py python3-matplotlib

Ubuntu's ``libfftw3-mpi-dev`` is built against OpenMPI, so with MPICH build
FFTW-MPI once from source (about five minutes):

.. code-block:: bash

   cd /tmp && wget http://www.fftw.org/fftw-3.3.10.tar.gz && tar xzf fftw-3.3.10.tar.gz && cd fftw-3.3.10
   ./configure --enable-shared --enable-threads --enable-mpi MPICC=mpicc
   make -j4 && sudo make install && sudo ldconfig

which installs it under ``/usr/local`` (hence ``FFTW_PREFIX=/usr/local``
below). With OpenMPI throughout (``openmpi-bin``, ``libhdf5-openmpi-dev``,
``libfftw3-mpi-dev``) no source build is needed.

On a cluster load the corresponding modules (``module load <mpi>
<hdf5-parallel> <fftw>``); if the parallel-HDF5 wrapper has another name
than ``h5pfc``, pass it with ``make FC=<wrapper>``.

Build
=====

From the repository root:

.. code-block:: bash

   make print-config                     # shows the compiler, MPI, HDF5 and FFTW that will be used
   make FFTW_PREFIX=/usr/local -j4       # builds ./Scorpio (drop FFTW_PREFIX if FFTW is in a standard path)
   make clean && make FFTW_PREFIX=/usr/local -j4     # clean rebuild

The build writes objects and module files under ``build/``; the only
product is the executable ``./Scorpio``. Run ``make`` again after any change
in ``src/`` — only the changed files and their dependants are recompiled.

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - compile error
     - fix
   * - ``Cannot open included file 'fftw3-mpi.f03'``
     - add ``FFTW_PREFIX=/path/to/fftw``
   * - ``h5pfc: command not found``
     - install the parallel-HDF5 development package / load its module, or ``make FC=<wrapper>``
   * - ``undefined reference to fftw_mpi_…``
     - FFTW-MPI was built against a different MPI than the compiler wrapper uses

Quick checks that the build works:

.. code-block:: bash

   make test              # FFTW-MPI smoke test (4 ranks)
   make gridid-smoke      # runs two tiny built-in cases (gridID 1 and 2) in a temporary directory

Run
===

The executable reads ``problem.nml`` from the directory it is started in and
writes all output there, so use one directory per experiment:

.. code-block:: bash

   mkdir -p ~/runs/test1 && cd ~/runs/test1
   printf '&problem_config\n  gridID = 800\n/\n' > problem.nml
   mpirun -n 4 /path/to/scorpio_modern/Scorpio 2>&1 | tee run.log

``-n 4`` is the number of MPI ranks; use 1, 2, 4, 8, 16, … and mesh sizes
that are multiples of 16. ``gridID`` selects the case (800 is the 20 pc
cloud; the map of all cases is in ``docs/LegacyGridID.md`` and the most used
ones are listed on the Problem File page). Long runs:
``nohup mpirun -n 8 /path/to/Scorpio > run.log 2>&1 &``.

.. note::
   On WSL, when the run directory is on a Windows drive (``/mnt/c/…``),
   ``export HDF5_USE_FILE_LOCKING=FALSE`` first, otherwise the parallel HDF5
   writes fail. If ``problem.nml`` was edited on Windows, strip the carriage
   returns (``sed -i -e 's/\r$//' -e '$a\' problem.nml``).

Output
======

Snapshots are HDF5 files ``g<gridID>_<nnnn>.h5`` (``g0800_0000.h5`` is the
initial condition, then one file per ``dtout``), self-describing and
readable with ``h5py``; the dataset list and a minimal reader are on the
Problem File and Hydrodynamics pages. Optional VTK/PVTS output for ParaView
is switched on per case with ``write_vtk``. AMR runs additionally write
checkpoints ``c<caseID>_<nnnn>.h5``.

A run ends when ``t`` reaches ``tend``; the cloud case can also stop itself
when the collapse becomes unresolved (the Truelove criterion — a final
snapshot is written first) or when its wall-clock limit is reached. All of
these can be continued with the restart options (Problem File page).

Validation
==========

.. code-block:: bash

   ./validation/validate.sh                 # quick gate (~5 min): builds, then Tier 0–2 tests + low-β robustness
   ./validation/validate.sh --full          # + the heavy Tier 3 cases (~15–20 min)
   ./validation/validate.sh --update-refs   # once on a new machine: reference values are machine-local
   ./validation/amr_gates.sh --np "1 8"     # AMR battery
   ./validation/gravity_analytic.sh         # analytic self-gravity gates

Machine-independent invariants (no NaN, :math:`\max|\nabla\cdot\boldsymbol{B}|<10^{-9}`,
mass drift :math:`<10^{-11}`, positivity, convergence order) are enforced on
every run; ``validation/README.md`` lists the tests.

Where to go next
================

- ``GETTING_STARTED.md`` in the repository: a step-by-step guide to running
  and modifying the 20 pc cloud case without knowing the solver internals.
- :ref:`ch:problem_file` — every setting, namelist group and ``SCORPIO_*``
  option.
- :ref:`ch:hydro` — equations, code units with cgs/SI conversions, the
  variable layout.
- :ref:`ch:methods` and :ref:`ch:time_integration` — what the code does in
  a step and which routine does it.

Platforms used by the group
===========================

- tianhexy
- cluster2
- hbli-s2
- hbli-s1
- scorpio
- gpu2/mike
- paulsr
- stor2
- nas3

Create a new Python environment (optional)
==========================================

#. Initialize conda

   .. code-block:: bash

      conda init

#. Activate the base environment

   .. code-block:: bash

      conda activate

#. Create a new environment from scratch with Python 3, or by cloning base

   .. code-block:: bash

      conda create -n my_new_env python=3 numpy h5py matplotlib
      conda create -n my_new_env --clone base

#. Activate it and install what you need

   .. code-block:: bash

      conda activate my_new_env
      conda install <package_name>

#. Update conda and all packages

   .. code-block:: bash

      conda update --all

For more information, visit
https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html
