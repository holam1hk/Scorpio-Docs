.. _ch:getting_started:

***************
Getting Started
***************

Scorpio is a grid-based finite-volume MPI code for hydrodynamics,
magnetohydrodynamics and two-fluid ambipolar diffusion. This page is what you
need on a machine where it is already installed: build it, run a case, find
the output.

.. note::
   The libraries Scorpio needs (MPI, parallel HDF5, FFTW-MPI) are already
   installed on the group's machines, so there is nothing to install. If you
   are setting up a new machine, or building with a different compiler or
   library paths, see :ref:`ch:developer`.

You will also want Python 3 with ``numpy`` and ``h5py`` (and ``matplotlib``)
to look at the output; see the appendix below if you prefer your own
environment.

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

If ``make`` fails with a missing library or a missing ``h5pfc``, the
machine is not set up as expected — the build variables and the common
compile errors are in :ref:`ch:developer`.

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
- :ref:`ch:developer` — setting up a new machine, build options, and
  working on the code.

Machines where Scorpio is set up
================================

The group's machines, with the libraries already installed:

- tianhexy
- cluster2
- hbli-s2
- hbli-s1
- scorpio
- gpu2/mike
- paulsr
- stor2
- nas3

Appendix: your own Python environment (optional)
================================================

Only needed if the system Python does not have ``numpy``, ``h5py`` and
``matplotlib``, or you want them isolated:

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
