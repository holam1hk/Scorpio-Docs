.. _ch:quickstart_spectrum:

********************************************************
Quickstart: driven turbulence with spectrum compensation
********************************************************

``gridID = 401`` is a driven hydrodynamic turbulence box that can, in a
second phase, push the velocity spectrum towards a prescribed
:math:`E(k) \propto k^{-2}` shell slope. The driver is
``SpectrumCompensation3DHD`` in ``src/spectrumCompensation.f03``; the
compensation itself is ``apply_spectrum_to_BGV_3D`` in
``src/gridModule_dtcalc.f03``. This page is the equivalent of
:ref:`ch:quickstart` for that case.

What the case is
================

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - setting
     - value
   * - physics
     - hydrodynamics, no magnetic field, no self-gravity (``GravConst`` is set in the driver but gravity is never enabled)
   * - equation of state
     - adiabatic with :math:`\Gamma = 1.00001` — near-isothermal
   * - solver
     - HLLC (``solverType = 3``), minmod, ``CFL = 0.3``, ``nbuf = 2``
   * - box
     - :math:`[-0.5, 0.5]^3` (1 pc across), periodic
   * - mesh
     - 256³ (``SCORPIO_SPEC_NMESH``)
   * - initial state
     - uniform :math:`\rho_0 = 17.6827` :math:`M_\odot`/pc³ (≈ 300 particles cm⁻³) at rest, :math:`c_s = 0.19196` km/s, :math:`P_0 = c_s^2\rho_0/\Gamma`
   * - driving
     - repeated kicks (``DT_mode = 1``), ``n_turb = 500`` of them, one every
       ``dt_turb`` = (crossing time)/500 ≈ 0.0104; solenoidal (:math:`\zeta = 1`),
       driving wavenumber 2, total energy :math:`\tfrac{1}{2}(10\,c_s)^2` — an rms Mach number of about 10
   * - run length
     - ``tend = 5.0`` code time ≈ one crossing time (5.21), 50 snapshots (``dtout = tend/50``)
   * - phase 2
     - **dormant by default**: ``t_phase2 = 5.5`` is later than ``tend``, so the default run is plain Otto (2017) driving

Times are in code units, 1 unit = 0.978 Myr (:ref:`ch:hydro`).

Run it
======

Make a run directory next to the code directory, put ``problem.nml`` in it
and start the code there:

.. code-block:: bash

   cd .. && mkdir turb_test && cd turb_test

.. code-block:: fortran

   &problem_config
     gridID = 401
   /

.. code-block:: bash

   mpirun -n 8 ../scorpio_modern/Scorpio 2>&1 | tee run.log

The full 256³ run is a production job (a few hundred core-hours). For a
first look, shrink it with the three environment variables this case reads:

.. code-block:: bash

   SCORPIO_SPEC_NMESH=64 SCORPIO_SPEC_TEND=0.5 mpirun -n 8 ../scorpio_modern/Scorpio 2>&1 | tee run.log

.. list-table::
   :header-rows: 1
   :widths: 30 14 56

   * - variable
     - default
     - meaning
   * - ``SCORPIO_SPEC_NMESH``
     - 256
     - cells per direction (the box stays 1 pc; it also sets ``k_end`` = mesh/2)
   * - ``SCORPIO_SPEC_TEND``
     - 5.0
     - end time; the snapshot interval follows as ``tend/50``
   * - ``SCORPIO_SPEC_TPHASE2``
     - 5.5
     - when phase 2 starts. Larger than ``tend`` keeps the compensation dormant;
       set it below ``tend`` to exercise it

``SCORPIO_DRIVING_SPECTRUM`` (``expo`` | ``kolmogorov`` | ``burgers``) chooses
the shape of the *driving* spectrum, and the usual run options apply
(:ref:`sec:run_options`). The driver calls ``setRunOption`` for
``SCORPIO_DRIVING_SPECTRUM``, ``SCORPIO_DUAL_ENERGY`` and ``SCORPIO_PPM``
with empty strings, which keeps whatever the environment says; put a value
there to fix a choice in code (with PPM also set ``nbuf = 3``).

What the log shows
==================

.. code-block:: text

    calcDT: driving spectrum = expo (default, km^6 exp(-8 km/k_d))
    gridModule.f03: writing data to g0401_0000.h5
    Applying Otto 2017 driving at t=   1.0419019746883021E-002
    testSuite.f03: gridID =          401  t =    1.1E-002  dt =    6.9E-004  t_count_turb =            1
    ...
    Spectrum compensation test: elapse time =    43.7  secs

one ``Applying Otto 2017 driving`` line per kick (500 of them over a
crossing time), one line per step, and a timing line at every snapshot.
With phase 2 active you also get, per compensation:

.. code-block:: text

    Fitted spectrum slope=  -1.78 ; applying prescribed slope=  -2.00
    Phase 2: Applying scheduled spectrum compensation at t=   2.0040000000000002
    Snapshot spectrum-compensation age =   1.0E-003 at snapshot t=   2.1

The fitted slope is a diagnostic only — the target always uses the
prescribed :math:`-2`.

Phase 2: the spectrum compensation
==================================

After ``t_phase2`` the driver stops the clock exactly on a schedule of
compensation times (every ``dt_turb_phase2`` = 5×10⁻³, offset so that each
snapshot is written a fixed lag ``snapshot_comp_lag`` = 10⁻³ after the last
compensation) and calls ``applySpectrumToBGV`` there. One compensation does:

1. gather the velocity field into the FFTW slab layout and transform it;
2. bin :math:`|\hat{\boldsymbol{v}}|^2` into unit-width shells with
   half-integer centres (interior :math:`k_x` planes carry the Hermitian
   weight 2);
3. normalise a :math:`k^{-2}` target on the *uncorrected* resolved shells and
   form a correction factor per shell; a straight-line fit of the measured
   spectrum is printed but not used as the target;
4. multiply every mode with :math:`k_{\rm start} < k \le k_{\rm end}` by the
   square root of the interpolated correction factor, raised to
   ``relaxation_fraction`` and clipped to ``max_amplitude_factor``;
5. transform back, scatter to the hydro decomposition, and add the change in
   kinetic energy to the total energy so the **internal energy — and therefore
   the pressure — is unchanged**.

The knobs are in-code, at the top of the driver:

.. list-table::
   :header-rows: 1
   :widths: 32 18 50

   * - variable
     - default
     - meaning
   * - ``ApplySpectrumCompensation``
     - ``.true.``
     - master switch. ``.false.`` turns 401 into a plain driven-turbulence run: no schedule, no dt clipping, no final-step compensation
   * - ``ApplyFinalSpectrumCompensation``
     - ``.false.``
     - one extra compensation on the last step (only with the master switch on)
   * - ``dt_turb_phase2``
     - 5e-3
     - interval between compensations
   * - ``snapshot_comp_lag``
     - 1e-3
     - how long after a compensation the snapshot is written
   * - ``spectrum_k_start``
     - 10
     - modes below this wavenumber are left alone (the driving range)
   * - ``spectrum_k_end``
     - mesh/2
     - largest complete isotropic shell
   * - ``spectrum_relaxation``
     - 1.0
     - exponent on the correction factor: 1 = apply it fully, <1 = approach the target gradually. Must be in (0, 1]
   * - ``spectrum_factor_max``
     - unlimited
     - clip on the per-mode amplitude factor; must be ≥ 1

To see the compensation act, run with ``SCORPIO_SPEC_TPHASE2`` below the end
time, e.g.

.. code-block:: bash

   SCORPIO_SPEC_NMESH=64 SCORPIO_SPEC_TEND=1.0 SCORPIO_SPEC_TPHASE2=0.5 \
       mpirun -n 8 ../scorpio_modern/Scorpio 2>&1 | tee run.log

Changing the physics
====================

Two routines in ``src/spectrumCompensation.f03``, and one trap:

- ``SpectrumCompensation3DHD`` — the driver: mesh, box, EOS, solver, CFL,
  end time, the driving schedule and all the phase-2 knobs above.
- ``initSpectrumCompensation3DHD`` — the initial condition: ``rho0``, the
  sound speed, the initial velocity (zero) and the pressure.
- **The trap:** the driving parameters that actually take effect
  (``drivingWN_DT``, ``Energy_DT``, ``zeta_DT``) are set in the *initial-condition*
  routine and overwrite what the driver computed. ``n_turb`` appears in both
  places — change it in both, or the energy per kick and the number of kicks
  will disagree.

``bdrySpectrumCompensation3DHD`` exists but does nothing: the box is
periodic, so the MPI exchange fills every ghost zone.

After an edit, ``make`` and run in a fresh directory. Adding your own case
instead of editing this one: :ref:`sec:new_case`.

Output and the spectrum
=======================

Snapshots are ``g0401_0000.h5`` … with the hydrodynamic variable set:
``den``, ``momx``, ``momy``, ``momz`` and ``ene`` (in a hydro run the energy
sits in slot 5, not 8; :ref:`sec:conserved`). To measure the shell spectrum
the compensation targets — :math:`E(k) = \sum_{|k|\in{\rm shell}} \tfrac{1}{2}|\hat{\boldsymbol{v}}|^2`:

.. code-block:: python

   import h5py, numpy as np, matplotlib.pyplot as plt

   f = h5py.File("g0401_0020.h5", "r")
   n  = int(f["nMesh"][()][0])
   nb = int(np.asarray(f["nbuf"][()]).reshape(-1)[0])
   cube = lambda name: np.transpose(f[name][()])[nb:n+nb, nb:n+nb, nb:n+nb]

   rho = cube("den")
   v   = [cube(m) / rho for m in ("momx", "momy", "momz")]

   kf  = np.fft.fftfreq(n, d=1.0/n)                      # integer wavenumbers
   kx, ky, kz = np.meshgrid(kf, kf, kf, indexing="ij")
   kmag = np.sqrt(kx**2 + ky**2 + kz**2)

   power = sum(np.abs(np.fft.fftn(c) / n**3)**2 for c in v) * 0.5
   kbin  = np.arange(0.5, n//2 + 1)                      # unit-width shells
   idx   = np.digitize(kmag.ravel(), np.arange(1, n//2 + 1))
   Ek    = np.bincount(idx, weights=power.ravel(), minlength=len(kbin))[:len(kbin)]

   plt.loglog(kbin, Ek, label="measured")
   plt.loglog(kbin, Ek[9] * (kbin/kbin[9])**-2, "--", label=r"$k^{-2}$")
   plt.xlabel("k"); plt.ylabel("E(k)"); plt.legend(); plt.show()

With phase 2 active the measured curve should approach the dashed line above
``spectrum_k_start`` within a few compensations; below it the driving range
is untouched by construction.

Cost
====

The default 256³ run with 500 driving kicks over a crossing time is a
cluster job; 64³ finishes in minutes on a workstation and is enough to see
the machinery work. Cost scales as cells × steps, so doubling the mesh costs
roughly 16×. The driving and the compensation both remap between the hydro
decomposition and the FFTW slab layout internally, so any ``mpirun -n``
works — the old requirement of a slab decomposition is gone
(:ref:`ch:parallelization`).
