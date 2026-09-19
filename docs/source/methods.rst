.. _ch:methods:

*****************
Numerical methods
*****************

Introduction
============

This page describes the numerical machinery that was added to Scorpio after
the method paper (Cheng et al. 2025, RASTI 4, 1): the constrained-transport update and its
upwinded EMF, the piecewise-parabolic reconstruction, the dual-energy
formulation with its positivity ladder, the multigrid self-gravity solver, the
block-structured adaptive mesh refinement, and the IMEX ion–neutral coupling
schemes. For every method it gives what it does, when it acts, how it is
switched on or off, and what it was validated against. The switches themselves
are listed in :ref:`ch:problem_file`; the underlying equations and code units
are in :ref:`ch:hydro`.

All of these additions are *additive*: with the default settings every run
that does not use the new feature reproduces the previous results bit for
bit, which is checked by the validation gate (``validation/validate.sh``).


Time integration and reconstruction
===================================

Scorpio integrates the conservation laws with a second-order Heun (RK2)
scheme. Within each stage the flux divergence is evaluated sweep by sweep
(method of lines: every directional sweep reads the same state ``q``), so a
stage is *not* a sequence of one-dimensional updates and the constrained
transport can be assembled once per stage from all sweeps (below). The time
step is the global CFL minimum over all cells and ranks (``griddt``), clipped
so that snapshots land exactly on the output times.

Interface states are reconstructed from the primitive variables in a frame
rotated to the sweep direction:

- **PLM** (default): piecewise-linear with one of the slope limiters
  ``limiterType`` 0–3 (zero, van Leer, MC — called ``fslop`` in the code —,
  minmod). On the 2D field-loop test the limiter is the dominant lever on
  magnetic-energy retention: minmod 0.66, van Leer 0.84, MC 0.88 after eight
  crossings, all with :math:`\max|\nabla\cdot\boldsymbol{B}|\sim10^{-16}`.
- **PPM** (``SCORPIO_PPM=1``): the piecewise-parabolic method of Colella &
  Woodward (1984) — MC-limited differences, fourth-order face values (their
  eq. 1.6), parabola monotonisation (eq. 1.10) — expressed as effective left/
  right slopes so the unchanged Riemann interface consumes the PPM states. It
  needs **three ghost cells** (``nbuf = 3``): the face interpolation and the
  neighbour parabola read :math:`i-2\ldots i+3`, and with only two ghost cells
  the two ranks sharing a face would build different parabolae, breaking
  conservation and :math:`\nabla\cdot\boldsymbol{B}=0`. With ``nbuf = 2`` the
  code prints one warning and falls back to PLM. Field-loop retention rises
  from 0.884 to 0.931 (2D) and from 0.941 to 0.962 (3D); ``np = 1`` and
  ``np = 4`` agree to every digit. Available in the Cartesian adiabatic and
  isothermal MHD solvers (2D/3D), single- and two-fluid; not in 1D or polar
  coordinates.

A positivity-preserving reconstruction guard (Zhang & Shu type :math:`\theta`
limiting) drops a cell to first order if its reconstructed states would have
negative density or pressure. The Riemann solvers are HLL/HLLC (hydro) and
HLL/HLLD (MHD); the HLLD kernel hardening is described with the failsafes
below.


Constrained transport (CT)
==========================

The magnetic field is stored on cell faces — the slots ``bxl``/``bxr`` etc.
of the state vector (:ref:`sec:conserved`) — and advanced with the discrete
curl of an edge-centred electromotive force,
:math:`\partial_t \boldsymbol{B}_{\rm face} = -\nabla\times\boldsymbol{E}_{\rm edge}`.
Because the discrete divergence of a discrete curl vanishes identically,
:math:`\nabla\cdot\boldsymbol{B}` stays at round-off for the whole run
(:math:`10^{-16}` to :math:`10^{-12}` in every test) — there is no cleaning
step and nothing to tune.

The edge EMF is assembled from the face-centred EMFs that the Riemann solver
returns in each sweep. Two assemblies exist:

- **Upwinded corner EMF** (default): the contact-mode upwinding of Gardiner &
  Stone (2005, 2008), the scheme used by Athena/Athena++. Each sweep caches its
  face EMFs and donor-cell weights (:math:`w = 1, 0, \tfrac{1}{2}` by the sign
  of the face mass flux); the corner EMF is completed when the last sweep has
  run (:math:`E_z` at the second sweep in 2D; :math:`E_x, E_y, E_z` at the
  third sweep in 3D) and the full staggered curl is applied once with the
  correct per-direction spacings, so non-cubic cells are handled correctly.
  Default since 2026-07-09 in the 2D and 3D adiabatic and isothermal solvers.
- **Centred average** (``SCORPIO_UPWIND_EMF=0``): the plain arithmetic
  average of the four face EMFs around an edge (Balsara & Spicer 1999). More
  diffusive; kept to reproduce older runs.

What the A/B tests showed: on smooth advection (field loop) the upwinded EMF
changes the retained magnetic energy by only +0.3 %, but at shocks it keeps
sharper field structures and therefore *increases* the number of cells whose
thermal energy becomes negative (2D MHD blast: 30 → 248 events). It is an
accuracy option, not a robustness one — the robustness comes from the dual
energy and the flux correction below, which are always on. The upwinded EMF
also removed a transverse instability of the split-sweep centred scheme that
appeared on multi-level AMR meshes at CFL 0.8 (see AMR).

On AMR meshes the fine-level edge EMFs are restricted onto the coincident
coarse edges at every coarse–fine interface, so the coarse-side face field is
updated with the same EMF as the fine side, and newly refined blocks receive
their face field by the divergence-preserving prolongation of Balsara (2001).


Dual energy
===========

**The problem.** In strongly magnetized (low-:math:`\beta`) flow the thermal
energy is a small difference of large numbers,
:math:`P = (\Gamma-1)\left(E - \tfrac{1}{2}\rho u^2 - \tfrac{1}{2}B^2\right)`,
and round-off in the conservative update makes it negative: thousands of
events per run in the two-fluid blast, sub-Alfvénic turbulence and the
low-:math:`\beta` ion blast. Controlled experiments showed that this is a
cancellation problem, not a flaw of the Riemann solver, the CT scheme, the
reconstruction or the ambipolar source: a single-fluid replica of the ion state
reproduces it exactly.

**The method.** In 2D/3D Cartesian adiabatic MHD an entropy density
:math:`\sigma = P/\rho^{\Gamma-1}` is carried as one extra state variable (the
last slot, ``nvar``; :ref:`sec:conserved`) and advected as a passive scalar.
After each update, in cells where the thermal energy is a small fraction of the
total (:math:`e_{\rm int} < {\rm thr}\times E`, ``SCORPIO_DE_THR`` = 10\ :sup:`-3`),
the pressure is taken from :math:`\sigma` instead of from the energy difference,
and the total energy is re-synchronised to it. Elsewhere :math:`\sigma` is reset
from the conservative pressure, so the scheme stays conservative wherever the
energy equation is well conditioned. Mass, momentum and :math:`\boldsymbol{B}`
are never touched: it is a switch between two representations of the same
state, **not a floor**. The entropy is seeded generically in the
initial-condition routines and re-seeded on restart, and in the two-fluid
drivers it is re-synchronised after the ion–neutral drag step.

**Behaviour.** At high :math:`\beta` the branch never fires, so every standard
test is bitwise unchanged; where it fires it removes the negatives completely
(two-fluid AD blast 6240 → 0 events, low-:math:`\beta` blast 4296 → 0, 3D blast
13864 → 0) with mass exact and :math:`\max|\nabla\cdot\boldsymbol{B}|\sim10^{-12}`.

**Diagnostics.** Each rank prints ``[dual-energy] t=… rank=R cells_rescued(p<0)=N``
the moment it heals cells (rank-local, no communication; ``SCORPIO_DE_PRINT=0``
silences it), and at the end of the run rank 0 prints the cumulative fire rate
(entropy-branch updates / all updates). ``SCORPIO_DUAL_ENERGY=0`` disables the
recovery for A/B comparisons.


Positivity failsafe ladder
==========================

When a cell still ends a step with :math:`P \le 0` or :math:`\rho \le 0` the
step is *not* accepted. The response is a ladder; each rung is tried before the
next:

1. **Dual energy** (above) — for adiabatic MHD.
2. **Localized first-order flux correction (FOFC)** — the troubled cells are
   flagged (over the ghost-inclusive range, so the flags agree across MPI
   boundaries), and the step is redone with the faces touching flagged cells
   computed from piecewise-constant states and the positivity-safe HLL flux;
   every other face keeps HLLD with its high-order reconstruction. The CT EMFs
   are rebuilt from the corrected fluxes, so the correction is conservative and
   :math:`\nabla\cdot\boldsymbol{B}`-preserving. At most 4–5 passes (first order
   spreads by one cell per pass). Follows Stone et al. (2020, Athena++) and
   PLUTO.
3. **Time-step halving** — the step is redone at :math:`\Delta t/2`, up to ten
   times; the CFL step is restored on the next step, so healthy runs pay nothing.
   ``SCORPIO_DT_HALVE=0`` removes this rung.
4. **Abort** — the run stops with a message rather than integrating an
   unphysical state.

``SCORPIO_LEGACY_FAILSAFE=1`` restores the previous global response (switch the
whole domain from HLLD to HLL, or halve the global time step), which is robust
but over-diffusive and made strong-shock runs crawl.

Two further hardenings act inside the HLLD kernel itself, before any of the
above is needed: the degeneracy switches of the star and double-star states
(Miyoshi & Kusano 2005, eqs. 44–47) use a relative threshold
(``SCORPIO_HLLD_EPS``, default 10\ :sup:`-8`, as in Athena) instead of exact
floating-point equality — at low :math:`\beta` the exact test essentially never
fired and the code divided by a tiny cancellation residual; and an optional
per-face admissibility check of the HLLD fan (``SCORPIO_HLLD_STARCHECK=1``) can
fall back to HLL on that face alone.

``SCORPIO_HEALTH=1`` prints one line per step (only when something happened):
cells rescued, FOFC passes, HLLD fallback and degeneracy counts, globally
summed.


Self-gravity: FFT and multigrid solvers
=======================================

The potential satisfies :math:`\nabla^2\Phi = 4\pi G\rho`, and the force
:math:`\boldsymbol{g} = -\nabla\Phi` enters the momentum and energy equations
as source terms at every RK2 stage (``gridModule_rk2.f03`` calls
``calcSelfgravity`` on the stage state). Two solvers share one interface
(``enableSelfgravitySolver`` / ``initSelfgravitySolver`` / ``calcSelfgravity``,
:ref:`sec:grid_api`) and one boundary flag ``sgBdryType``:

.. list-table::
   :header-rows: 1
   :widths: 18 41 41

   * -
     - FFT (``sgSolverType = 0``)
     - multigrid (``sgSolverType = 1``)
   * - isolated (``sgBdryType = 0``)
     - zero-padded Green-function convolution (exact free-space boundary);
       kernels for the three force components and the potential
       (:math:`-1/r`, zero self-cell)
     - 3D: **James (1977) screening-charge boundary** (default): a first solve
       with zero Dirichlet ghosts, the induced surface charge is read off its
       boundary layer, its potential is summed onto the ghost positions and a
       second solve uses those values (cf. Ricker 2008; Moon, Kim & Ostriker
       2019). ``SCORPIO_MG_ISOBC=multipole`` selects the older multipole
       Dirichlet boundary (monopole + dipole + quadrupole). 2D: complex
       multipole expansion to order 4.
   * - periodic (``sgBdryType = 1``)
     - spectral Poisson solve; mean-zero potential
     - the mean density is subtracted, residuals are projected to zero mean and
       the potential is written in the zero-mean gauge
   * - parallel layout
     - FFTW-MPI slab decomposition; the hydro decomposition is remapped to it
       (and the force back) with packed ``MPI_ALLTOALLV`` exchanges, so the
       hydro grid no longer has to use ``dims = (1,1,nprocs)``
     - the multigrid levels follow the hydro decomposition; ghost exchange per
       level (``exchgBdryMPI_sgMG``)
   * - algorithm
     - one forward transform of the (zero-padded) density and one inverse
       transform per force component and for the potential
     - V-cycles (up to 200) with, in 3D, a weighted-Jacobi smoother (:math:`\omega = 0.8`,
       3 pre- and 3 post-smoothing sweeps, 60 on the coarsest level) to a
       relative residual of 10\ :sup:`-6`

Both solvers write ``gphi`` and ``sgfx``, ``sgfy``, ``sgfz`` to the snapshots.
On the uniform grid the FFT solver is the default and the multigrid is
selectable (``cloud_sgsolver`` for the cloud case); on AMR the multigrid is the
only solver (policy: all AMR gravity gates were calibrated against it, and it
scales better at large size, Tomida & Stone 2023).

Validation uses analytic gates rather than reference runs: the uniform sphere
and the Maclaurin spheroid (case 359, Chandrasekhar 1969 / Ricker 2008), the
pressureless free-fall collapse against the cycloid solution (360), a periodic
sinusoidal density against its exact potential (361), and the 2D/3D periodic and
isolated multigrid cases (357, 358, 557, 558). Because the isolated FFT solve
uses the exact free-space boundary, comparing it with the isolated multigrid
solve on the same problem isolates the boundary-truncation part of the
multigrid error (the two-sphere gate).


Adaptive mesh refinement (AMR)
==============================

Mesh
----

Scorpio's AMR is a **block-structured octree** in the FLASH/Paramesh style:
the domain is tiled with blocks of ``amr_bsx × amr_bsy × amr_bsz`` cells
(default 16; the base mesh must be a multiple of the block size and the block
size a multiple of :math:`2^{\rm max\_level}`). Refining a block replaces it
by 4 (2D) or 8 (3D) children with the same number of cells and half the
spacing; the refinement ratio is 2 and neighbouring leaves may differ by at
most one level (2:1 proper nesting, enforced by flag dilation). The tree
metadata (level, position, parent/children, owner rank, neighbours) is
replicated on every rank, so neighbours are found without communication; the
block *data* are distributed in Morton order and re-partitioned at every
regrid. Every block is an ordinary ``grid`` object, so the uniform-grid solver
kernels run on it unchanged.

Time stepping and coupling between levels
-----------------------------------------

All levels advance with **one global time step** (the CFL minimum over all
blocks; no subcycling), with the same Heun RK2 as the uniform code, so after
every stage all levels are at the same time and the coarse–fine ghost fill
needs only spatial interpolation:

- ghost cells at same-level faces are copies of the neighbour's interior;
  at coarse–fine faces the coarse data are prolonged with conservative
  limited slopes, and for the face-centred field with the divergence-preserving
  scheme of Balsara (2001);
- after each stage the fine solution is **restricted** onto its parent
  (volume average; area average for the face field);
- **refluxing**: at every coarse–fine face the coarse cell next to the
  interface is corrected with the difference between the area-weighted sum of
  the fine fluxes and its own flux, so mass, momentum and energy are conserved
  exactly across levels;
- **EMF matching**: the coarse-side face field at interfaces is updated with
  the restricted fine edge EMF, which keeps :math:`\nabla\cdot\boldsymbol{B}`
  at round-off on every level.

Refinement criteria
-------------------

Blocks are flagged with the Löhner (1987) second-derivative indicator (the
FLASH form) on density — and on total energy and :math:`B^2` in 3D — and
refined when the indicator exceeds ``amr_refine_thr`` (0.8), derefined when all
siblings are below ``amr_deref_thr`` (0.2) for ``amr_deref_hyst`` (2)
consecutive regrids. Optionally the **Jeans criterion** of Truelove et al.
(1997) refines wherever :math:`\lambda_J <` ``amr_jeans_n`` :math:`\Delta x`
(off by default). Regridding runs every ``amr_regrid_every`` steps; the initial
condition is built by refining ``amr_max_level`` times and re-evaluating the
analytic initial condition on each new block, so the :math:`t=0` fine data are
exact rather than interpolated. ``amr_flag_buffer`` dilates the refined region
by whole block shells. ``SCORPIO_AMR_AUDIT=1`` prints every flag decision and
the conservation totals before and after each regrid.

Gravity on the AMR mesh
-----------------------

Self-gravity is solved with the multigrid solver on the base-level covering
grid (``gio``): the block densities are restricted onto it, the potential is
solved there, and the forces are scattered back to the blocks. This reproduces
the uniform run exactly for ``amr_max_level = 0`` and is the certified default.
With ``amr_fac = 1`` the **FAC** (fast adaptive composite) extension solves the
refined levels as well: nested level solves with Dirichlet data prolonged from
the level below (stage 1, block-local multigrid within a red-black block
colouring — a multiplicative Schwarz iteration; ``SCORPIO_FAC_SOLVER=sor``
selects the older pointwise SOR), and the two-way :math:`\tau`-correction
feedback of the fine solution onto the coarse equation (stage 2,
``amr_fac_iters`` composite passes). On the 3D uniform-sphere gate this cuts the
force error on the refined region by about 40 % relative to the base-grid
solve. Known limits: refined blocks touching a *physical* boundary get a
zero-gradient potential boundary there (keep refined regions interior for
isolated problems), and levels that fill most of the domain converge slowly
(the block-coupling limit of the Schwarz iteration; a coarse-space correction
is the planned fix), which is why deep, domain-filling refinement of the cloud
case is expensive.

Output, checkpoints, restart
----------------------------

Each output writes the solution **projected onto the base grid**
(``g<caseID>_<fnum>.h5``, readable exactly like a uniform snapshot) plus a
**checkpoint** with every leaf block (``c<caseID>_<fnum>.h5``;
``SCORPIO_AMR_CHECKPOINT=0`` disables). A run is restarted from checkpoint
:math:`N` with ``&amr_restart amr_fstart = N``; the rank count may change. Per
step the AMR driver prints the conservation report and, for the cloud case,
the ``HINNY801`` line with the maximum density, the Jeans resolution
:math:`\min(\lambda_J/\Delta x)` and the block count per level.

Validation and limits
---------------------

The gate battery (``validation/amr_gates.sh``, cases 700–727) checks
single-level parity with the uniform code (bitwise), conservation
(:math:`\lesssim 10^{-11}` drift), :math:`\max|\nabla\cdot\boldsymbol{B}|` at
round-off, far-field invariance, and ``np = 1`` versus ``np = 4`` bitwise
equality through dynamic regridding. 3D hydro is certified to three levels of
refinement (case 716, conservation :math:`3.5\times10^{-12}`); 3D MHD to two
levels — the multi-level MHD runs exposed a transverse instability of the
centred CT scheme at CFL 0.8, which is what motivated making the upwinded EMF
the default (with it, the static and dynamic two-level MHD cases conserve to
:math:`10^{-13}`–:math:`10^{-12}` with :math:`\nabla\cdot\boldsymbol{B}` exactly
zero).

Current limits of the AMR path: Cartesian coordinates only; no turbulence
driving and no two-fluid ambipolar diffusion on AMR; PPM is wired but not
validated there; the global
time step makes a three-level run cost roughly 4–8 times the wall time of its
base grid (subcycling is the designed mitigation).


Ion–neutral coupling integrators
================================

The two-fluid drag term is stiff when the coupling time
:math:`\tau = 1/(\alpha\rho_i)` is short compared with the time step. Four
integrators are available through ``SCORPIO_AD_SCHEME``:

.. list-table::
   :header-rows: 1
   :widths: 22 30 48

   * - value
     - scheme
     - properties
   * - ``split`` (default)
     - operator-split TR-BDF2 drag after each RK stage (the method of the
       paper; Tilley et al. 2012)
     - second order; at :math:`\Delta t/\tau \gtrsim 6` the splitting error
       damps the ion–neutral drift spuriously (:math:`\sim 1.1` e-foldings per
       period in the Alfvén-damping test)
   * - ``imex``
     - IMEX-SSP2(2,2,2), drag solved implicitly inside the RK stages
     - removes the splitting error
   * - ``imex322`` / ``athenak``
     - IMEX-SSP2(3,2,2) of Pareschi & Russo (2005, Table III), the scheme of
       AthenaK's ion–neutral module
     - L-stable and stiffly accurate (the stiff limit lands on the drag
       equilibrium); rings in stiff relaxation
   * - ``imexpp`` / ``imex2+`` / ``krapp``
     - IMEX(4,3,2) of Krapp et al. (2024) with the monotone
       :math:`\gamma_+ = 1 + 1/\sqrt{2}` root (AthenaK "imex2+")
     - same accuracy without the ringing, best ion-momentum conservation;
       needs CFL :math:`\le 0.35` (its first explicit stage advances
       :math:`1.707\,\Delta t`; a runtime guard refuses larger values)

The IMEX schemes exist for the 3D two-fluid driver; the 2D driver refuses them
rather than silently running the split scheme, and an unrecognised value
aborts. The validated production combination for driven ion–neutral turbulence
is HLLD + PPM + ``imexpp`` at CFL 0.3. The C-shock and Alfvén-damping harnesses
(cases 52/53 and 620/621, with an AthenaK comparison mode) are the reference
tests; see :ref:`ch:AD`.


Restart
=======

Uniform-grid runs restart from a snapshot through ``setTime(fstart = N)``,
which rebuilds the mesh, the variables and the solver settings from
``g<gridID>_<N>.h5`` and restores :math:`t` and all fields (the dual-energy
entropy is re-seeded from the restored state). The cloud driver exposes this as
the ``Restart`` flag / ``cloud_restart`` (:ref:`ch:problem_file`); a restarted
run reproduces the uninterrupted run bit for bit, also with a different number
of ranks. The t = 0 turbulence kick is part of the snapshot and is not
repeated; for periodic driving (``DT_mode = 1``) the kick counters travel in a
side file next to each snapshot so the schedule resumes on the exact step,
while the random phases of later kicks are new. AMR runs restart from their
checkpoints as described above.


Validation gate
===============

``validation/validate.sh`` rebuilds the code and runs the Tier 0–2 battery
(Brio–Wu, CP Alfvén wave with its convergence order, rotor, MHD blast, field
loop, Sod, Einfeldt rarefaction, Shu–Osher, Woodward–Colella, …) in about five
minutes; ``--full`` adds Orszag–Tang, the oblique and low-:math:`\beta` cases,
the two-fluid blast and the 3D blast. Machine-independent invariants are
enforced on every run — no NaN, :math:`\max|\nabla\cdot\boldsymbol{B}| < 10^{-9}`,
mass drift :math:`< 10^{-11}`, no negative pressure or density in the final
state, no legacy failsafe engagements, zero entropy-branch fires on
high-:math:`\beta` tests, Alfvén-wave convergence order :math:`> 1.4` — while
the reference values themselves are machine-local (``--update-refs`` once on a
new machine). ``validation/amr_gates.sh`` and ``validation/gravity_analytic.sh``
are the corresponding AMR and gravity batteries.
