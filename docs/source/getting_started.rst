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

With OpenMPI throughout (``openmpi-bin``, ``libhdf5-openmpi-dev``,
``libfftw3-mpi-dev``) no source build is needed.

On a cluster load the corresponding modules (``module load <mpi>
<hdf5-parallel> <fftw>``); if the parallel-HDF5 wrapper has another name
than ``h5pfc``, pass it with ``make FC=<wrapper>``.

Build
=====

From the repository root:

.. code-block:: bash

   make print-config      # shows the compiler, MPI, HDF5 and FFTW that will be used
   make                   # builds ./Scorpio
   make clean && make     # clean rebuild

The ``Makefile`` finds the libraries through the ``h5pfc`` wrapper and the
system paths; nothing has to be passed on the command line on a machine that
is set up. The build writes objects and module files under ``build/``; the
only product is the executable ``./Scorpio``. Run ``make`` again after any
change in ``src/`` — only the changed files and their dependants are
recompiled.

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - compile error
     - fix
   * - ``Cannot open included file 'fftw3-mpi.f03'``
     - FFTW is not in a standard path: ``make FFTW_PREFIX=/path/to/fftw`` (e.g. ``/usr/local`` for a source build)
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

The case to run is chosen in ``problem.nml``, a small text file in the
directory you start the code from (all output is written there too):

.. code-block:: fortran

   &problem_config
     gridID = 800
   /

Then run with MPI:

.. code-block:: bash

   mpirun -n 4 ./Scorpio

``-n 4`` is the number of MPI ranks (CPU cores); use 1, 2, 4, 8, 16, … and
mesh sizes that are multiples of 16. ``gridID`` selects the case — 800 is
the 20 pc cloud; every case is listed in :ref:`sec:all_cases`, and the
per-case settings (``&cloud_nml`` etc.) are on the Problem File page. For a
long run that should survive logging out:
``nohup mpirun -n 8 ./Scorpio > run.log 2>&1 &``.

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

``./validation/validate.sh`` rebuilds the code and runs the standard test
battery in about five minutes; run it before committing a change. The
batteries, what they check and how to read their reports are on the
:ref:`ch:validation` page.

Where to go next
================

- :ref:`ch:quickstart` — run the 20 pc cloud, look at the output, change
  its settings and its physics, restart it.
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
