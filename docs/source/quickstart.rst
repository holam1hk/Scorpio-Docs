.. _ch:quickstart:

**********************************
Quickstart: your first cloud run
**********************************

This page takes you from a compiled ``Scorpio`` to a finished simulation of
the 20 pc magnetized cloud (``gridID = 800``), its output, and the first
changes you will want to make — without needing the MHD details. It
assumes the build from :ref:`ch:getting_started` is done. Everything you
touch is in two places:

.. list-table::
   :header-rows: 1
   :widths: 55 30 15

   * - what
     - where
     - recompile?
   * - run settings: resolution, run length, gravity solver, driving on/off, restart
     - ``problem.nml`` in the run directory
     - no
   * - the physics of the cloud: density profile, field strength, sound speed, box
     - ``src/hinnyCloud.f03``
     - yes (``make``)

1. A two-minute test run
========================

Scorpio reads ``problem.nml`` from, and writes everything to, the directory
it is started in. Keep the runs out of the source tree: make a directory for
each one next to the code directory, in the folder above it.

.. code-block:: bash

   cd ..                        # the folder that holds the code directory
   mkdir cloud_test && cd cloud_test

Write ``problem.nml`` there with any editor:

.. code-block:: fortran

   &problem_config
     gridID = 800
   /
   &cloud_nml
     cloud_nmesh = 32, 32, 64        ! small mesh for a quick test
     cloud_tend_tff = 0.05           ! run length in free-fall times (production: 10)
     cloud_dtout_tff = 0.01          ! snapshot interval in free-fall times (production: 0.1)
     cloud_truelove = .false.        ! don't stop on the resolution criterion in a test
   /

and run:

.. code-block:: bash

   mpirun -n 4 ../scorpio_modern/Scorpio 2>&1 | tee run.log

``-n 4`` is the number of MPI ranks — use 1, 2, 4, 8, 16, … with mesh sizes
that are multiples of 16. ``../scorpio_modern/Scorpio`` is the executable
``make`` built; adjust the path if your code directory has another name. For
a long run that should survive logging out:
``nohup mpirun -n 8 ../scorpio_modern/Scorpio > run.log 2>&1 &``.

.. note::
   On WSL with the run directory on a Windows drive (``/mnt/c/…``), first
   ``export HDF5_USE_FILE_LOCKING=FALSE``. If ``problem.nml`` was edited on
   Windows and the code reports "invalid namelist content", strip the
   carriage returns: ``sed -i -e 's/\r$//' -e '$a\' problem.nml``.

What the log shows
------------------

At start-up, check that your settings arrived:

.. code-block:: text

    HinnyTestSuite: turbulence driving =  T
    HinnyTestSuite: nMesh =          32          32          64  tend =   7.6765000000000000E-002  dtout =   1.5353000000000000E-002
    HinnyTestSuite: gravity solver =FFT
    calcDT: driving spectrum = burgers (per-mode km^-4)
    HinnyTestSuite: turbulence kick #           1  of            1  at t =   0.0000000000000000
    gridModule.f03: writing data to g0800_0000.h5

then one line per step, ``HinnyTestSuite.f03: gridID= 800 t= … dt= …``, and
``writing data to g0800_0001.h5`` at every snapshot. A misspelled entry in
``&cloud_nml`` makes Fortran silently ignore the *whole* group — the
``nMesh =`` line is where you notice.

How a run ends
--------------

1. ``t`` reaches ``tend`` — normal.
2. ``Truelove criteria violated`` — the densest cell can no longer resolve
   its own collapse (Jeans length below 5 cells). A final snapshot is written,
   then the run aborts; this is the intended end of a collapse run, not a
   crash. At 32×32×64 the initial Jeans length is only ≈ 5.1 cells, so it
   fires almost immediately unless the resolution is raised (64×64×128
   gives ≈ 10) or ``cloud_truelove = .false.``.
3. ``Times up!`` — the built-in wall-clock limit (1 day 5 h; ``dd``/``hh``
   in the driver) was reached.
4. Ctrl-C / ``kill`` — the snapshots already written are fine.

After 2–4 you can continue from the last snapshot (section 4).

2. Look at the output
=====================

``g0800_0000.h5`` is the initial condition, ``g0800_0001.h5`` the state
after ``dt_out``, and so on. Each file is self-describing HDF5: ``t``,
``nMesh``, ``nbuf`` (ghost cells padded on *every* side of every array),
the coordinates ``xc1/xc2/xc3`` and widths ``dx1/dx2/dx3``, and the fields
``den``, ``momx/momy/momz``, the face-centred field ``bxl/bxr/byl/byr/bzl/bzr``,
``ene``, and with gravity ``gphi``, ``sgfx/sgfy/sgfz`` (units in
:ref:`ch:hydro`). Arrays come out of ``h5py`` in ``(z, y, x)`` order.

.. code-block:: python

   import h5py, numpy as np, matplotlib.pyplot as plt

   f = h5py.File("g0800_0005.h5", "r")
   nx, ny, nz = [int(v) for v in f["nMesh"][()]]
   nb = int(np.asarray(f["nbuf"][()]).reshape(-1)[0])
   t  = float(np.asarray(f["t"][()]).reshape(-1)[0])

   def cube(name):                              # (z,y,x) on disk -> (x,y,z), ghost cells removed
       return np.transpose(f[name][()])[nb:nx+nb, nb:ny+nb, nb:nz+nb]
   def axis(name, n):
       return np.asarray(f[name][()]).reshape(-1)[nb:n+nb]

   rho = cube("den")                            # Msun / pc^3
   vx  = cube("momx") / rho                     # km / s
   bz  = 0.5 * (cube("bzl") + cube("bzr"))      # code units; x 2.9 for microgauss
   x, z, dy = axis("xc1", nx), axis("xc3", nz), axis("dx2", ny)

   colden = (rho * dy[None, :, None]).sum(axis=1)          # integrate along y -> (x, z), Msun/pc^2
   plt.imshow(np.log10(colden).T, origin="lower", aspect="equal",
              extent=[x[0], x[-1], z[0], z[-1]])
   plt.colorbar(label="log10 column density [Msun/pc^2]")
   plt.xlabel("x [pc]"); plt.ylabel("z [pc]"); plt.title(f"t = {t:.3f} code units = {0.978*t:.3f} Myr")
   plt.show()

   dx, dz = axis("dx1", nx), axis("dx3", nz)
   vol = dx[:, None, None] * dy[None, :, None] * dz[None, None, :]
   print("max density:", rho.max(), " total mass [Msun]:", (rho * vol).sum())

While a run is going, ``grep "t=" run.log | tail`` shows where it is.

3. Change settings without recompiling
======================================

Every entry of ``&cloud_nml`` is optional; a missing one keeps the value in
the driver.

.. list-table::
   :header-rows: 1
   :widths: 22 18 60

   * - entry
     - default
     - meaning
   * - ``cloud_nmesh``
     - 32, 32, 64
     - cells in x, y, z (keep multiples of 16)
   * - ``cloud_tend_tff``
     - 10.0
     - run length in units of ``tff`` (a fixed 1.5353 code time units ≈ 1.5 Myr)
   * - ``cloud_dtout_tff``
     - 0.1
     - snapshot interval in units of ``tff``
   * - ``cloud_truelove``
     - ``.true.``
     - stop when the collapse becomes unresolved (see above)
   * - ``cloud_sgsolver``
     - 0
     - gravity solver: 0 = FFT, 1 = multigrid
   * - ``cloud_sgbdry``
     - 0
     - gravity boundary: 0 = isolated cloud, 1 = periodic
   * - ``cloud_driving``
     - ``.true.``
     - the initial turbulence kick on/off
   * - ``cloud_restart``
     - ``.false.``
     - continue from a snapshot (section 4)
   * - ``cloud_fstart``
     - 0
     - the snapshot number to continue from

A production-like run at twice the default resolution, without driving:

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

4. Restart a run
================

In the directory that holds the snapshot (say ``g0800_0042.h5``):

.. code-block:: fortran

   &cloud_nml
     cloud_restart = .true.
     cloud_fstart  = 42               ! continue from g0800_0042.h5
     cloud_tend_tff = 20.0            ! optional: a later end time than the original run
   /

and launch as before. The code prints ``RESTART from g0800_0042.h5``, takes
the mesh, the fields and the solver settings from that file, and continues
with ``g0800_0043.h5``, ``0044``, …

- The flag decides; ``cloud_fstart`` is only the number. The run stops with
  a clear message if the two disagree or the snapshot is missing.
- Only ``cloud_tend_tff`` / ``cloud_dtout_tff`` are taken from the new
  ``problem.nml``; mesh, sound speed, solver, etc. come from the snapshot.
- The number of MPI ranks may differ. A restarted run reproduces the
  uninterrupted run bit for bit.
- The initial turbulence kick is part of the snapshot and is never applied
  again. With repeated kicks (``DT_mode = 1`` in the driver) keep the small
  ``g0800_NNNN.turb`` file next to the snapshot — it carries the kick
  schedule.

5. Change the physics in the code
=================================

The driver ``cloud_20pc3_3DMHD`` (in ``src/hinnyCloud.f03``) sets the run;
the routine ``initcloud_20pc3_3DMHD`` below it defines the cloud. Search
for the variable name, change the value, run ``make``, run in a fresh
directory. Two things to know first:

- **Some numbers exist twice.** The AMR version of the case (``gridID =
  801``) has its own copies of the cloud numbers in ``amrHinnyIC3`` and of
  the sound speed in ``cloud_20pc3_3DMHD_AMR``. Change both, or the two
  cases simulate different clouds.
- **Units.** Lengths in pc, masses in M☉, velocities in km/s, so time is
  0.978 Myr, density M☉/pc³ (≈ 15 H₂ cm⁻³), and one unit of magnetic field
  is ≈ 2.9 μG (the code absorbs the 4π: ``B_code = B_Gauss/√4π``). The
  temperature follows from the sound speed: ≈ 25 K × (``sndspd``/0.3 km/s)².
  Full conversions in :ref:`ch:hydro`.

Run setup (in ``cloud_20pc3_3DMHD``)
-------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 22 48

   * - to change …
     - variable
     - notes
   * - resolution
     - ``nMesh(1..3)``
     - or ``cloud_nmesh``; multiples of 16
   * - box
     - ``leftBdry``, ``rightBdry``
     - −5..5, −5..5, −10..10 pc; the cloud shape assumes this box (``L`` and the ``5d0`` in ``rho_max``)
   * - run length, snapshot spacing
     - ``tend_tff``, ``dtout_tff``
     - in units of ``tff``; or the namelist
   * - free-fall time used for the two above
     - ``tff``
     - a fixed number (1.5353), not recomputed from the density: :math:`t_{\rm ff} = \sqrt{3\pi/(32 G \rho)}`
   * - sound speed / temperature
     - ``sndspd``
     - 0.3 km/s ≈ 25 K; also ``snd=0.3d0`` in the AMR driver
   * - gravity on/off, solver, boundary
     - ``SelfGravity``, ``sg_solver``, ``sg_bdry``
     - or ``cloud_sgsolver`` / ``cloud_sgbdry``
   * - turbulence kick on/off, strength, scale
     - ``DriveTurbulence``, ``E_turb_tot``, ``DT_scale``
     - energy in M☉ km² s⁻² (4.5 ≈ 9×10⁴³ erg); ``DT_scale = 2`` = eddies of half the box
   * - one kick or repeated kicks
     - ``DT_mode``, ``n_turb``, ``turn_over_time``
     - 0: one kick at t = 0; 1: ``n_turb`` kicks, one every ``turn_over_time/n_turb``
   * - solenoidal / compressive mix
     - ``zeta``
     - 1 solenoidal, 0 compressive, 0.5 natural
   * - driving spectrum
     - ``OPT_DRIVING_SPECTRUM``
     - ``'burgers'``, ``'kolmogorov'``, ``'expo'``
   * - wall-clock limit
     - ``dd``, ``hh``
     - days + hours; set it below the queue limit
   * - restart
     - ``Restart``, ``file_start``
     - or the namelist (section 4)
   * - VTK output for ParaView
     - ``write_vtk``
     -
   * - CFL number
     - ``CFL``
     - 0.4; lower is safer and slower
   * - solver, limiter, EOS, boundary type
     - ``solverType``, ``limiterType``, ``eosType``, ``boundaryType``
     - 5 (HLLD), 3 (minmod), 1 (isothermal), 3 (periodic): leave unless you know why

The cloud (in ``initcloud_20pc3_3DMHD``)
-----------------------------------------

The density is a cylinder along z with a flat core:
``rho = rho_c / (1 + (r/r_flat)**2)`` with :math:`r=\sqrt{x^2+y^2}`,
multiplied by ``exp(-alpha*(|z| - L/2))`` beyond ``|z| > L/2``, and floored
at ``rho_max`` (the value at the box edge, r = 5 pc).

.. list-table::
   :header-rows: 1
   :widths: 30 22 48

   * - to change …
     - variable
     - notes
   * - central density
     - ``rho_c``
     - 25.8866 M☉/pc³ ≈ 375 H₂ cm⁻³; also in ``amrHinnyIC3``
   * - core radius, cloud length, end taper
     - ``r_flat``, ``L``, ``alpha``
     - 2.5 pc, 15 pc, 1.5 pc⁻¹; also in ``amrHinnyIC3``
   * - density floor
     - the ``5d0`` in ``rho_max``
     - the box half-width; change it with the box
   * - field strength
     - ``b0``
     - 35.0842 ≈ 100 μG; also in ``amrHinnyIC3``
   * - field direction
     - ``b_dir``
     - 0 = along z, 1 = along x, other = no field (``amrHinnyIC3`` is z only)
   * - initial bulk velocity
     - ``u0``, ``v0``, ``w0``
     - km/s; if non-zero, change the momentum lines from ``rho_c * u0`` to ``q(i, j, k, 1) * u0`` — momentum is the *local* density times velocity
   * - profile exponent
     - ``p``
     - only affects the floor; the profile line has the exponent 1 written out — edit ``rho = …`` to change the power law

Each field component is set twice (``q(…,5)`` and ``q(…,9)`` for Bx, 6/10
for By, 7/11 for Bz): the two faces of the cell. Keep them equal for a
uniform field.

Recipes
-------

*Add solid-body rotation about z:* declare ``omega``, set it (km s⁻¹ pc⁻¹),
and replace the three momentum lines after ``q(i, j, k, 1) = DMAX1(rho, rho_max)`` with

.. code-block:: fortran

   q(i, j, k, 2) = q(i, j, k, 1) * (-omega * yc(j))
   q(i, j, k, 3) = q(i, j, k, 1) * ( omega * xc(i))
   q(i, j, k, 4) = 0.d0

*A different density profile:* replace the ``rho = (rho_c / …)`` line with
your formula in ``r``, ``xc(i)``, ``yc(j)``, ``zc(k)``; keep the
``DMAX1(rho, rho_max)`` floor or set ``rho_max`` to what you want.

*Keep the original and make your own case:* copy the three routines under
new names and register a new ``gridID`` — the worked example is in
:ref:`sec:new_case`.

Fortran reminders: real constants are ``2.5d0`` / ``1.d-3`` (the ``d``
makes them double precision), logicals ``.true.``/``.false.``, comments
start with ``!``. ``git diff src/hinnyCloud.f03`` shows what you changed.

6. Troubleshooting
==================

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - symptom
     - cause / fix
   * - ``main.f03: no runnable problem.nml entry``
     - no ``problem.nml`` in the current directory, or no ``gridID`` in it
   * - ``problemRegistry: invalid namelist content``
     - Windows line endings, or a group without its closing ``/``
   * - my ``&cloud_nml`` values are ignored
     - a misspelled entry name fails the whole group silently; compare with the table in section 3 and check the ``nMesh =`` start-up line
   * - ``Truelove criteria violated`` after a few steps
     - expected at 32×32×64; raise the resolution or set ``cloud_truelove = .false.`` for tests
   * - ``Times up!``
     - the 29-hour wall-clock limit; restart from the last snapshot
   * - HDF5 write errors on WSL
     - ``export HDF5_USE_FILE_LOCKING=FALSE`` on ``/mnt/…``
   * - the run got much slower after a mesh change
     - cost ∝ cells × steps and the step shrinks with the cell size: 2× resolution ≈ 16× the time
   * - ``dt`` shrinks towards zero
     - something is collapsing or a velocity is blowing up — look at the last snapshot; check the units of what you changed (``2.5`` instead of ``2.5d0``, a velocity in the wrong unit)
   * - I changed the code but nothing happened
     - run ``make`` after editing, and run the freshly built ``Scorpio`` (check its timestamp)

Cheat sheet
===========

.. code-block:: bash

   # build, in the code directory (after every edit in src/)
   make

   # a run directory next to the code directory, with its own problem.nml
   cd .. && mkdir exp1 && cd exp1
   export HDF5_USE_FILE_LOCKING=FALSE          # WSL on /mnt only
   mpirun -n 4 ../scorpio_modern/Scorpio 2>&1 | tee run.log

   # watch
   grep "t=" run.log | tail -3
   ls g0800_*.h5 | wc -l

   # continue a stopped run: add cloud_restart/cloud_fstart to problem.nml, same directory
   mpirun -n 4 ../scorpio_modern/Scorpio 2>&1 | tee -a run.log
