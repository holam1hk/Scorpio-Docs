.. _ch:validation:

**********
Validation
**********

Introduction
============

Scorpio is validated by scripted test batteries under ``validation/`` in the
code repository. Every gate prints one line, ``PASS name …`` or
``FAIL name …``, with the measured numbers and the criterion; a script's exit
code is the number of failed gates, so the batteries can be wired into a
commit hook or a job that runs on the server. The batteries never look at
pictures — the figures on the results dashboard illustrate gates that are
decided numerically.

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - script
     - time
     - what it covers
   * - ``validation/validate.sh``
     - ~5 min (``--full`` ~20 min)
     - the uniform-grid code: shock tubes, MHD standards, convergence order, low-β robustness
   * - ``validation/amr_gates.sh``
     - 30–60 min on 8 cores (``--full`` 1–3 h)
     - the AMR machinery: parity with the uniform code, conservation, div B, np-invariance, restart, multi-level
   * - ``validation/gravity_analytic.sh``
     - minutes (``--deep`` longer)
     - self-gravity against closed-form solutions
   * - ``validation/gravity_deep.sh``
     - hours (server)
     - deep AMR gravity runs of the cloud case
   * - ``validation/server_pipeline.sh``
     - 30–60 min (``--full`` 3–5 h)
     - build check + ``amr_gates`` + the analytic tier, in one report for a server run

The standard battery: ``validate.sh``
=====================================

.. code-block:: bash

   ./validation/validate.sh                 # quick gate: build + Tier 0/1/2 + low-β robustness
   ./validation/validate.sh --full          # + Tier 3: Orszag–Tang, oblique-B, low-β Alfvén, two-fluid blast, 3D blast
   ./validation/validate.sh --np 8          # MPI ranks (default 4); MPIRUN=srun is honoured
   ./validation/validate.sh --update-refs   # (re)measure the machine-local reference values

.. list-table::
   :header-rows: 1
   :widths: 10 34 12 44

   * - tier
     - test
     - ``gridID``
     - validates
   * - 0/2
     - Brio–Wu 2D
     - 12
     - shock capturing vs reference
   * - 0/2
     - circularly polarised Alfvén wave
     - 18
     - dispersion / dissipation (exact nonlinear solution)
   * - 0/2
     - rotor
     - 16
     - torsional Alfvén fronts, symmetry
   * - 0/2
     - MHD blast
     - 14
     - strong-shock robustness, positivity
   * - 0/2
     - field loop
     - 15
     - CT / div B, field dissipation
   * - 1
     - CP Alfvén convergence (N = 16 / 32 / 64)
     - 60
     - formal order of accuracy in smooth flow
   * - H
     - Sod shock tube 2D
     - 5
     - exact-Riemann anchor, HLL hydro path
   * - H
     - Einfeldt rarefaction
     - 8
     - near-vacuum positivity (the hydro analogue of low β)
   * - H
     - Shu–Osher
     - 9
     - shock–entropy interaction, reconstruction fidelity
   * - H
     - Woodward–Colella two-blast
     - 6
     - extreme-shock robustness, wall heating
   * - H*
     - double Mach reflection
     - 7
     - canonical strong-shock 2D structure
   * - H*
     - HD blast (Cartesian)
     - 37
     - Sedov-like symmetry, energy conservation
   * - 3
     - low-β ion blast
     - 56
     - dual-energy positivity at β ~ 4×10⁻⁴
   * - 3*
     - Orszag–Tang (t = 0.5)
     - 13
     - integrator regression standard
   * - 3*
     - oblique-B low-β blast
     - 57
     - CT / energy consistency off-axis
   * - 3*
     - low-β CP Alfvén wave
     - 58
     - Alfvén-wave retention where dual energy fires
   * - 3*
     - two-fluid AD blast
     - 54
     - ion–neutral coupling + AD source
   * - 3*
     - 3D low-β blast
     - 59
     - 3D dual energy + FOFC

(* = ``--full`` only.)

Two kinds of checks are made. **Reference comparisons** hold each test to
values measured on a known-good commit; because floating-point results differ
between compilers and machines at round-off, and nonlinear tests amplify that
to the 4th–6th digit, the references are *machine-local*: run
``--update-refs`` once on a new machine and commit the resulting
``validation/references.dat``. **Invariants** are machine-independent and
enforced on every run:

.. list-table::
   :header-rows: 1
   :widths: 60 40

   * - invariant
     - threshold
   * - no NaN in the output
     - 0 cells
   * - :math:`\max|\nabla\cdot\boldsymbol{B}|`
     - :math:`< 10^{-9}`
   * - mass conservation
     - drift :math:`< 10^{-11}`
   * - positivity in the final state
     - 0 cells with :math:`P<0` or :math:`\rho\le0`
   * - legacy failsafe engagements (global HLL switch, dt halving)
     - 0
   * - entropy-branch fire rate on high-β tests
     - exactly 0 (bitwise-identity guard for dual energy)
   * - CP-Alfvén :math:`L_1` convergence order, N = 16 → 64
     - > 1.4 (the second-order scheme gives ≈ 2)

To make it a commit gate, copy ``validation/pre-commit.sample`` to
``.git/hooks/pre-commit``; ``git commit --no-verify`` bypasses it in an
emergency.

The AMR battery: ``amr_gates.sh`` and ``server_pipeline.sh``
=============================================================

.. code-block:: bash

   ./validation/amr_gates.sh                 # quick set
   ./validation/amr_gates.sh --full          # + long 3D cases (718/719) and the deep three-level ladder
   ./validation/amr_gates.sh --np "1 4 8"    # rank counts; the first and last are compared bitwise
   ./validation/server_pipeline.sh --np "1 8" [--full]   # build check + amr_gates + analytic tier, one report

The gates run the AMR cases 700–727 (:ref:`sec:all_cases`) and check:

- single-level **parity** with the uniform-grid driver (bitwise);
- **conservation** of mass and energy through refluxing and regridding
  (drift :math:`\lesssim 10^{-11}`);
- :math:`\max|\nabla\cdot\boldsymbol{B}|` at round-off on every level with the
  CT/EMF matching;
- **np-invariance**: the first and last rank counts of ``--np`` give
  identical results, bitwise (for zero-field backgrounds in 2D the tolerance
  is :math:`10^{-30}`, denormal dust);
- **restart**: continuing from a checkpoint reproduces the uninterrupted run;
- **multi-level**: static and dynamic refinement to two (MHD) and three
  (hydro) levels, far-field invariance.

The deliverable is ``<work>/amr_gate_report.txt`` (or ``PIPELINE_REPORT.txt``),
a few kB; ``logs.tar.gz`` is produced only when something failed. The
snapshots stay on the machine that ran the battery.

Gravity: ``gravity_analytic.sh`` and ``gravity_deep.sh``
=========================================================

``gravity_analytic.sh`` holds the self-gravity solvers to closed-form
answers rather than to reference runs:

.. list-table::
   :header-rows: 1
   :widths: 8 24 12 56

   * - gate
     - test
     - case
     - against
   * - G1
     - uniform sphere, centred and offset by 0.25, FFT and multigrid (James boundary), 64³
     - 358
     - Newton's field; the offset sphere is the sharp isolated-boundary test (a centred sphere cannot probe the boundary)
   * - G2
     - FAC composite sphere and a sphere cut by the level-1 boundary (AMR)
     - 721
     - block-level force error vs the analytic sphere
   * - G3
     - smooth Poisson convergence 32 → 64
     - 727
     - exact sinusoidal potential; full second order
   * - G4
     - Jeans linear growth (FAC forces, periodic gravity)
     - 724
     - exact growth rate σ = 2.50663 (Jeans 1902)
   * - G5
     - two off-centre spheres (AMR)
     - 725
     - exact superposed field; isolated boundary beyond the monopole

``--deep`` extends G3 to 128 → 256, G5 to 128³ and G2 to a 64 → 128
convergence order. ``gravity_deep.sh`` is the server-scale battery of the
cloud case on AMR (np-invariance, the L0/L1/L2 ladder at 64×64×128, the
isolated James-boundary twin, smooth Poisson at 128 → 256); its
``GRAVITY_DEEP_REPORT.txt`` is the only file that needs to come back.

Reading a report
================

.. code-block:: text

   PASS G1 sphere_mg_offset fmax=4.2059e-02 fl2=1.6904e-03 phimax=2.0289e-03 (max 5.5e-2/2.2e-3/3.0e-3; certified <=4.21e-2/1.69e-3/2.03e-3)
   PASS G2 fac_sphere max=2.0077E-02 l2=1.2309E-03 (max 2.5e-2/1.6e-3; certified 2.01e-2/1.23e-3)
   PASS G3 poisson_order_32_64 order=2.01 (min 1.8; certified 2.01 - full 2nd order through the interface)

Each line is one gate (from ``GRAVITY_ANALYTIC_REPORT.txt``): the measured
numbers, then in parentheses the threshold the script applied and the value
certified on the reference machine; a ``SUMMARY`` line closes the report. A FAIL in a
*reference* gate on a new machine usually means the references have not been
measured there yet (``--update-refs``); a FAIL in an *invariant* is a real
regression. The methods and the numbers behind the certified thresholds are
summarised in :ref:`ch:methods`.

Results dashboard
=================

The report files are collected into a history and rendered as one
self-contained web page with the figures of the calibration runs
(``validation/dashboard/collect.py`` → ``docs/validation/runs/*.json``,
``validation/dashboard/build.py`` → ``docs/validation/index.html``; only the
Python standard library is needed):

.. code-block:: bash

   python3 validation/dashboard/collect.py validation/pipeline_work/PIPELINE_REPORT.txt \
                                           validation/gravity_analytic_work/GRAVITY_ANALYTIC_REPORT.txt
   python3 validation/dashboard/build.py

The page lists every gate's latest verdict, the per-run detail and the
per-gate history, and shows the attached figures (Brio–Wu, Sod, Alfvén
wave, Sedov, uniform and two-sphere gravity, Jeans dispersion, collapse,
momentum conservation). It will be linked from here once it is published
alongside this manual.
