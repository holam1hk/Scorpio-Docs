.. _ch:turbulence:

******************
Turbulence Driving
******************

Introduction
============

Turbulence is driven by adding a random, divergence-free and/or compressive
velocity perturbation with a prescribed power spectrum to the whole box —
a *kick*. The perturbation is built in Fourier space with FFTW-MPI,
transformed to real space, shifted to zero net momentum, and rescaled so
that it injects exactly the requested kinetic energy. Kicks are applied
either once at :math:`t = 0` (``DT_mode = 0``) or repeatedly on a schedule
(``DT_mode = 1``). The driving needs a periodic box (``boundaryType = 3``,
``periods = .true.``).

Call path
=========

.. code-block:: text

   driver (e.g. cloud_20pc3_3DMHD, hinnyCloud.f03)
   ├─ g%drivingWN_DT = DT_scale;  g%Energy_DT = E_turb;  g%zeta_DT = zeta;  g%netmom*_DT = 0
   ├─ g%enableDrivingTurbulence(DT_mode)           gridModule.f03 — enable_DT, DT_mode, fftw_mpi_init
   ├─ … setMesh, setVariable (creates the driving FFT plans: DTPlan3D, gridModule_fft_plans.f03) …
   ├─ initial kick (fresh start only):
   │     g%setDTenergyfaction(dt_turb, dt_turb)      DTenergyfaction = 1
   │     g%calcDrivingTurbulence_MD(g%q)             the kick
   │     g%exchangeBdryMPI(g%q, g%winq); g%setBoundary(g%q)
   └─ time loop, DT_mode = 1 only, before evolveGridRK2:
         if (t_accum_turb ≤ dt_turb .and. t_accum_turb + dt > dt_turb .and. t_count_turb < n_turb)
            g%setDTenergyfaction(t_accum_turb, dt_turb); g%calcDrivingTurbulence_MD(g%q)
            exchangeBdryMPI; setBoundary; t_accum_turb = 0; t_count_turb += 1
         t_accum_turb += dt

   calcDrivingTurbulence_MD(q)                         gridModule.f03 — dispatch on ndim
   └─ calcDT3D_MD(this, …, q)                          gridModule_dtcalc.f03
      ├─ create_phase_space_profile_3D                 the perturbation in k-space (below)
      ├─ fftw_mpi_execute_dft_c2r ×3                   δv_x, δv_y, δv_z on the FFTW slab layout
      ├─ remap_dt_vector_fft_to_hydro3d                back to the hydro decomposition (MPI_ALLTOALLV)
      ├─ momentum_shift_3D                             subtract the mass-weighted mean → zero net momentum
      ├─ energy_shift_3D                               scale so the injected kinetic energy = Energy_DT × DTenergyfaction
      └─ q(mom) += ρ δv                                 (density and B untouched)

The perturbation in Fourier space
=================================

``create_phase_space_profile_3D`` fills every mode :math:`\boldsymbol{k}`
(wavenumbers in units of :math:`2\pi/L`, :math:`k = |\boldsymbol{k}|`) with

.. math::

   \hat{\boldsymbol{v}}(\boldsymbol{k}) = A_0(k)\,\mathsf{P}(\boldsymbol{k})\,\boldsymbol{a}\;e^{2\pi i\,\phi},
   \qquad
   \mathsf{P}_{ij} = \zeta\,\delta_{ij} + (1-2\zeta)\,\frac{k_i k_j}{k^2},

where :math:`\boldsymbol{a}` are three standard normal random numbers
(``randn``, Box–Muller), :math:`\phi` a uniform random phase (``rand1``),
:math:`\mathsf{P}` the Federrath et al. (2010) projection operator —
:math:`\zeta = 1` purely solenoidal, :math:`\zeta = 0` purely compressive,
:math:`\zeta = 0.5` the natural mix — and

.. math::

   A_0(k) = \sigma\sqrt{\frac{E_d(k)}{(3\zeta-2)\zeta+1}}

with the per-mode driving spectrum :math:`E_d` = ``driving_spectrum_3D(k, k_d)``
selected by ``SCORPIO_DRIVING_SPECTRUM`` (:math:`k_d` = ``drivingWN_DT``);
:math:`\sigma` is an overall normalisation that the energy rescaling below
makes irrelevant:

.. list-table::
   :header-rows: 1
   :widths: 22 40 38

   * - value
     - per-mode power :math:`E_d(k)`
     - notes
   * - ``expo`` (default)
     - :math:`k^6 e^{-8k/k_d}`
     - peaked near :math:`k_d`; Otto et al. (2017)
   * - ``kolmogorov``
     - :math:`k^{-11/3}` for :math:`k \ge 2`, else 0
     - shell spectrum :math:`E(k) \propto k^{-5/3}`
   * - ``burgers``
     - :math:`k^{-4}` for :math:`k \ge 2`, else 0
     - shell spectrum :math:`E(k) \propto k^{-2}`; the cloud driver selects this in code

The random numbers come from Fortran's intrinsic generator with its default
(OS-dependent) seed, so two runs produce different realisations of the
same statistics, and a restarted run's later kicks differ from an
uninterrupted run's (the schedule and the injected energy are identical).

Normalisation and schedule
==========================

- ``momentum_shift_3D`` subtracts the mass-weighted mean of the
  perturbation (minus any requested ``netmom*_DT``) so the kick adds no net
  momentum.
- ``energy_shift_3D`` scales the perturbation so that the *change* in kinetic
  energy, including the cross term with the existing velocity field, equals
  ``Energy_DT × DTenergyfaction`` (Lester's thesis, eqs. 3.11–3.13).
  ``Energy_DT`` is the energy per kick, ``E_turb_tot / n_turb`` in the cloud
  driver, in :math:`M_\odot\,{\rm km^2\,s^{-2}}`; ``setDTenergyfaction(a, b)``
  sets the fraction to :math:`a/b` (1 for the initial kick).
- ``DT_mode = 0``: one kick at :math:`t = 0`. ``DT_mode = 1``: the first kick
  at :math:`t = 0`, then a kick on the step before ``t_accum_turb`` would
  cross ``dt_turb = turn_over_time / n_turb``, until ``n_turb`` kicks are
  done. The counters ``t_count_turb`` / ``t_accum_turb`` are the driver's;
  the cloud driver stores them in ``g<gridID>_<fnum>.turb`` next to every
  snapshot so that a restart resumes the schedule on the exact step
  (:ref:`ch:problem_file`).
- After every kick the ghost zones must be refreshed
  (``exchangeBdryMPI`` + ``setBoundary``): the kick fills interior cells only.

Variants
========

``calcDT3D`` / ``calcDT2D`` are the older versions without the
momentum/energy shift; ``calcDT3D_AD`` drives the two-fluid case (neutral
grid, optionally both); ``apply_spectrum_to_BGV_3D`` (case 401,
``spectrumCompensation.f03``) re-shapes the *existing* velocity spectrum
instead of adding a perturbation. Cases 100, 399, 400, 401 and 1688 exercise
the driving on its own; the cloud case (800) uses it as the initial
turbulence.
