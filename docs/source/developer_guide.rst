.. _ch:developer:

***************
Developer Guide
***************

Introduction
============

This page is for **setting up a new machine** and for **working on the
code**. On the group's machines Scorpio and its libraries are already
installed — if you only want to run simulations, go to
:ref:`ch:getting_started` and :ref:`ch:quickstart` instead.

Toolchain
=========

Scorpio needs three libraries, and they must all be built against the
**same** MPI:

- an MPI library (MPICH or OpenMPI) with a Fortran compiler (gfortran);
- **parallel HDF5** with the Fortran interface. Its compiler wrapper
  ``h5pfc`` is what the ``Makefile`` calls, and it supplies the MPI and HDF5
  flags automatically;
- **FFTW3 with MPI support** (``libfftw3_mpi``), used by the FFT self-gravity
  solver and by the turbulence driving.

Python 3 with ``numpy`` and ``h5py`` (and ``matplotlib``) is needed for the
validation scripts and for looking at the output.

Ubuntu / WSL with MPICH
-----------------------

.. code-block:: bash

   sudo apt install build-essential gfortran make mpich libhdf5-mpich-dev hdf5-tools \
                    libfftw3-dev python3-numpy python3-h5py python3-matplotlib

Ubuntu's ``libfftw3-mpi-dev`` package is built against OpenMPI, so with MPICH
build FFTW-MPI once from source (about five minutes):

.. code-block:: bash

   cd /tmp && wget http://www.fftw.org/fftw-3.3.10.tar.gz && tar xzf fftw-3.3.10.tar.gz && cd fftw-3.3.10
   ./configure --enable-shared --enable-threads --enable-mpi MPICC=mpicc
   make -j4 && sudo make install && sudo ldconfig

That installs it under ``/usr/local``, which the build then needs to be told
about once: ``make FFTW_PREFIX=/usr/local``. Mixing the two MPIs is the
classic trap — if a duplicate ``libfftw3_mpi.so.3`` from the OpenMPI package
is still on the system, remove it so ``ldconfig`` resolves to the MPICH
build.

With OpenMPI throughout (``openmpi-bin``, ``libhdf5-openmpi-dev``,
``libfftw3-mpi-dev``) no source build is needed.

On a cluster
------------

.. code-block:: bash

   module load <mpi> <hdf5-parallel> <fftw>      # names are site-specific; ask the admin
   make FFTW_PREFIX=$FFTW_ROOT                   # or whatever variable the FFTW module sets

If the parallel-HDF5 wrapper is not called ``h5pfc`` on the site, pass its
name: ``make FC=<wrapper>``.

Check what the build will use
-----------------------------

.. code-block:: bash

   make print-config

prints the compiler, the include and library flags, and the detected MPI,
HDF5 and FFTW versions. ``mpirun``, ``h5pfc`` and ``gfortran`` must all be
found.

Build options
=============

Plain ``make`` is enough on a machine that is set up; these are the knobs for
one that is not, and for development builds.

.. list-table::
   :header-rows: 1
   :widths: 24 22 54

   * - variable
     - default
     - meaning
   * - ``FC``
     - ``h5pfc``
     - the compiler wrapper (it brings MPI and parallel HDF5); e.g. ``make FC=mpifort``
   * - ``FFTW_PREFIX``
     - empty (system paths)
     - where FFTW lives, e.g. ``/usr/local`` for a source build
   * - ``HDF5_PREFIX``
     - empty
     - only needed when ``FC`` is not an HDF5 wrapper
   * - ``FFLAGS``
     - ``-O3 -g -fimplicit-none -Wall -Wextra -ffree-line-length-none``
     - compiler flags; use ``-O0 -g -fcheck=all -fbacktrace`` for a debug build
   * - ``BUILD_DIR``
     - ``build``
     - where objects and ``.mod`` files go
   * - ``EXEC``
     - ``Scorpio``
     - name of the executable

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - target
     - what it does
   * - ``make`` / ``make all``
     - build ``./Scorpio`` (serial is safe; ``-j4`` works because the module dependencies are declared)
   * - ``make clean``
     - remove ``build/``, the ``.mod`` files and the executable
   * - ``make print-config``
     - show the toolchain and the flags (``COLOR=0`` for plain output)
   * - ``make smoke`` / ``make smoke-run`` / ``make test``
     - build and run the FFTW-MPI smoke test on 4 ranks
   * - ``make gridid-smoke``
     - run two tiny built-in cases (``gridID`` 1 and 2) in a temporary directory.
       ``GRIDID_SMOKE_IDS="1 2 5"``, ``GRIDID_SMOKE_NP=4`` and
       ``KEEP_GRIDID_SMOKE_DIR=1`` change the set, the rank count and whether the
       outputs are kept

Compile errors
--------------

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - message
     - cause / fix
   * - ``Cannot open included file 'fftw3-mpi.f03'``
     - FFTW is not in a standard path: ``make FFTW_PREFIX=/path/to/fftw``
   * - ``h5pfc: command not found``
     - parallel HDF5 not installed or its module not loaded; or pass ``make FC=<wrapper>``
   * - ``undefined reference to fftw_mpi_…``
     - FFTW-MPI was built against a different MPI than the compiler wrapper uses
   * - ``Error: Symbol '…' at (1) has no IMPLICIT type``
     - a typo in a new edit: every routine has ``implicit none``
   * - ``make: Nothing to be done for 'all'``
     - nothing changed since the last build — not an error

The build writes objects and module files under ``build/``; the only product
is the executable ``./Scorpio``. After any change in ``src/``, run ``make``
again — only the changed files and their dependants are recompiled. A
change in ``gridModule.f03`` rebuilds almost everything, a change in a case
module only that module and the two dispatch submodules.

Working on the code
===================

- :ref:`ch:code_structure` — every source file and what it holds, the
  program flow, the ``grid`` object, and where to look for what. It also
  carries the two recipes you are most likely to need,
  :ref:`sec:new_case` and :ref:`sec:source_terms`, and the table of
  :ref:`module-level run options <sec:globals>`.
- :ref:`ch:time_integration` — which routine calls which, from ``main``
  down to the flux kernels.
- :ref:`ch:methods` — what the newer machinery (CT, PPM, dual energy, the
  positivity ladder, multigrid gravity, AMR, IMEX coupling) does and how it
  was validated.

Conventions worth keeping: arrays with ghost cells are declared
``1-nbuf : nMesh+nbuf`` and interior loops run ``1 : nMesh``; real constants
are written ``1.d0``; every routine has ``implicit none``; only rank 0
prints (``if (myid .eq. 0)``); additions since 2026 are tagged ``[OPUS]`` /
``[FABLE]`` in the source, with ``CHANGES.md`` in the code repository as
their index.

Before committing
-----------------

.. code-block:: bash

   ./validation/validate.sh          # ~5 min: build + the standard battery and the invariants

:ref:`ch:validation` describes the batteries, the machine-local reference
values (``--update-refs`` once per machine) and the pre-commit hook. For a
change that touches the AMR or gravity paths, also run
``./validation/amr_gates.sh`` and ``./validation/gravity_analytic.sh``.

This documentation
==================

The site is Sphinx with the readthedocs theme; the sources are the ``.rst``
files in ``docs/source/`` of the documentation repository, and readthedocs
rebuilds the site from its ``main`` branch on every push
(``.readthedocs.yaml`` → ``docs/source/conf.py``,
``docs/requirements.txt``). To preview a change locally:

.. code-block:: bash

   pip install -r docs/requirements.txt
   python -m sphinx -b html docs/source _build/html    # then open _build/html/index.html

Page-specific styling lives in ``docs/source/_static/custom.css``.
