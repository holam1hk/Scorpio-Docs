.. _ch:riemann:

***************
Riemann Solvers
***************

Introduction
============

``riemannSolverModule.f03`` contains two layers:

- the **sweep solvers** ``solver<EOS><Physics><Dim>`` — one call is one
  directional sweep over the whole local block: rotate into the sweep frame,
  reconstruct the interface states, call the flux kernel on every face,
  update the output state, and (for MHD) assemble the constrained-transport
  update. They are called by ``rk2_1D/2D/3D`` through the procedure pointer
  ``rieSolver`` (:ref:`ch:time_integration`);
- the **flux kernels** ``flux<Solver><EOS><Physics>1D`` — the approximate
  Riemann solver proper, acting on one face.

Which routine runs
==================

``rk2_*`` picks the sweep solver from the equation of state and the
solver type; the sweep solver then picks the flux kernel from the solver
type and the slope function from the limiter type:

.. list-table::
   :header-rows: 1
   :widths: 16 16 26 42

   * - ``eosType``
     - ``solverType``
     - sweep solver (1D / 2D / 3D)
     - flux kernel
   * - 1 isothermal
     - 1
     - ``solverIso1D/2D/3D``
     - ``fluxExactIsoHD1D`` (exact isothermal Riemann solver; ``…2D``/``…3D`` variants in the 2D/3D sweeps)
   * - 1 isothermal
     - 2
     - ``solverIso1D/2D/3D``
     - ``fluxHLLIsoHD1D`` (``…2D``/``…3D``)
   * - 1 isothermal
     - 4
     - ``solverIsoMHD1D/2D/3D``
     - ``fluxHLLIsoMHD1D``
   * - 1 isothermal
     - 5
     - ``solverIsoMHD1D/2D/3D``
     - ``fluxHLLDIsoMHD1D`` (isothermal HLLD)
   * - 2 adiabatic
     - 2
     - ``solverAdi1D/2D/3D``
     - ``fluxHLLAdiHD1D`` (``…2D``/``…3D``)
   * - 2 adiabatic
     - 3
     - ``solverAdi1D/2D/3D``
     - ``fluxHLLCAdiHD1D`` (``…2D``/``…3D``)
   * - 2 adiabatic
     - 4
     - ``solverAdiMHD1D/2D/3D``
     - ``fluxHLLAdiMHD1D``
   * - 2 adiabatic
     - 5
     - ``solverAdiMHD1D/2D/3D``
     - ``fluxHLLDAdiMHD1D`` (Miyoshi & Kusano 2005), falling back to ``fluxHLLAdiMHD1D`` on a degenerate fan
   * - polytropic
     - —
     - ``solverPoly1D/2D``
     - ``fluxHLLPolyHD1D/2D`` (``polyGamma``, ``polyK``; 1D/2D only)

The cloud case runs isothermal HLLD: ``solverIsoMHD3D`` + ``fluxHLLDIsoMHD1D``.

The kernel interface
====================

All kernels have the same signature,

.. code-block:: fortran

   subroutine flux<name>(ql, qr, slopeL, slopeR, gam, flux, nvar)

and begin by building the two interface states from the cell values and the
limited slopes handed to them by the sweep:

.. code-block:: fortran

   UL = ql + 0.5d0*slopeL     ! state on the left of the face  (right face of cell i)
   UR = qr - 0.5d0*slopeR     ! state on the right of the face (left face of cell i+1)

so the reconstruction (PLM or PPM) is entirely the sweep's business — a
kernel only ever sees two states. With PPM the sweep passes different slopes
for the two faces of a cell (``SL`` and ``SLR``), which is how the parabola
is expressed without changing the kernels.

The sweep frame
===============

A sweep solver always works on "normal" and "transverse" components, so the
2D/3D solvers map the physical variable slots onto a sweep frame before the
reconstruction and back after the flux is computed. The index maps for the
11-variable MHD layout (``den, momx, momy, momz, bxl, byl, bzl, ene, bxr,
byr, bzr``) are:

.. list-table::
   :header-rows: 1
   :widths: 25 25 25 25

   * - sweep-frame slot
     - ``dd = 1`` (x-sweep)
     - ``dd = 2`` (y-sweep)
     - ``dd = 3`` (z-sweep)
   * - 1 density
     - 1
     - 1
     - 1
   * - 2 normal momentum
     - 2
     - 3
     - 4
   * - 3, 4 transverse momenta
     - 3, 4
     - 2, 4
     - 3, 2
   * - 5 normal B (left face)
     - 5
     - 6
     - 7
   * - 6, 7 transverse B (left faces)
     - 6, 7
     - 5, 7
     - 6, 5
   * - 8 energy
     - 8
     - 8
     - 8
   * - 9 normal B (right face)
     - 9
     - 10
     - 11
   * - 10, 11 transverse B (right faces)
     - 10, 11
     - 9, 11
     - 10, 9

The exchanged components carry a sign (the ``coef``/``g`` factors in the
code) so that the rotated frame stays right-handed; the fluxes are rotated
back with the same factors before they are applied. ``(cx, cy, cz)`` is the
unit step along the sweep, used for all neighbour look-ups.

What a sweep does, step by step
===============================

Taking ``solverAdiMHD3D`` as the reference (the other sweeps are subsets):

1. Select ``fluxPtr`` and ``slope`` as in the table above; allocate the slope
   arrays ``SL`` (and ``SLR`` for PPM).
2. Convert every cell (ghosts included) to primitive variables in the sweep
   frame; for adiabatic runs the pressure comes from the energy, or from the
   dual-energy entropy where the energy difference is unreliable.
3. Slopes: ``SL(i,j,k,m) = slope(W(i-1), W(i), W(i+1))`` with the limiter
   function (:ref:`ch:limiter`); with ``SCORPIO_PPM=1`` the CW84 face values
   and monotonised parabola replace them (needs ``nbuf = 3``). The normal
   field component keeps zero slope (it is continuous across the face by
   construction).
4. Positivity guard: if the reconstructed face states of a cell would have
   :math:`\rho \le 0` or :math:`P \le 0`, its slopes are scaled by the largest
   :math:`\theta \in [0,1]` that keeps them positive (Zhang & Shu type).
5. For every face, call the flux kernel on ``(ql, qr, slopeL, slopeR)``.
   During a first-order flux-correction redo (``fofc_pass``), faces touching
   a flagged cell use zero slopes and the HLL kernel instead. For the
   dual-energy variable the flux is the upwinded :math:`\sigma/\rho` times the
   mass flux. On AMR meshes the physical block-edge fluxes are captured for
   the reflux step.
6. Update the hydro variables of the output state,
   ``q2 = q1 − dt/dx_normal · (F(i+½) − F(i−½))`` (in the cylindrical
   coordinates with the metric factors and ``source2D``).
7. Constrained transport: cache the face EMFs (the transverse components of
   the field flux) and the upwind weights of this sweep; on the last sweep
   assemble the corner EMFs — Gardiner & Stone (2005) upwinded by default,
   arithmetic average with ``SCORPIO_UPWIND_EMF=0`` — and update the six
   face-field slots with the discrete curl. This keeps
   :math:`\nabla\cdot\boldsymbol{B}` at round-off (:ref:`ch:methods`).
8. Dual energy: recover :math:`P` from :math:`\sigma` in cells with
   :math:`e_{\rm int} < {\rm thr}\,E` and re-synchronise :math:`E`; elsewhere
   re-synchronise :math:`\sigma` from :math:`E`.

HLLD robustness
===============

Two hardenings live inside ``fluxHLLDAdiMHD1D`` / ``fluxHLLDIsoMHD1D`` and
are on by default:

- the degeneracy tests of the star and double-star states use the relative
  threshold ``SCORPIO_HLLD_EPS`` (10\ :sup:`-8`, as in Athena) instead of
  exact floating-point equality — at low :math:`\beta` the exact tests never
  fired and the kernel divided by a cancellation residual;
- with ``SCORPIO_HLLD_STARCHECK=1`` the kernel additionally checks that the
  contact speed lies strictly inside :math:`(S_L, S_R)` and the star total
  pressure is positive, and otherwise returns the HLL flux for that face (off
  by default: the check also triggers on strong but benign expansion fans).

References: Toro, Spruce & Speares (1994) for HLLC; Miyoshi & Kusano (2005)
for HLLD; Gardiner & Stone (2005, 2008) for the CT EMF; Mignone (2007) for
the isothermal HLLD variant.
