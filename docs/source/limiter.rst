.. _ch:limiter:

****************************
Slope Limiter, PLM and PPM
****************************

Introduction
============

The interface states handed to the Riemann kernels are reconstructed from
the primitive variables of three neighbouring cells. The default
reconstruction is piecewise linear (PLM) with a slope limiter from
``limiterModule.f03``; the piecewise-parabolic method (PPM) can be selected at
run time. Both live inside the sweep solvers of ``riemannSolverModule.f03``
(:ref:`ch:riemann`); the limiter module only supplies the slope functions.

Call path
=========

.. code-block:: text

   solver<…>3D(this, q, q1, q2, dd)                       riemannSolverModule.f03
   ├─ select case (this%limiterType)
   │     0 → slope => zslop      1 → slope => vslop
   │     2 → slope => fslop      3 → slope => minmod     (limiterModule.f03)
   ├─ per cell and variable:  SL(i,j,k,m) = slope(W(i-1), W(i), W(i+1))       [PLM]
   │   [use_ppm] Dp = fslop(...) (MC-limited differences)
   │             AFc = ½(W_i + W_{i+1}) − (Dp_{i+1} − Dp_i)/6   (CW84 eq. 1.6, 4th-order face value)
   │             (aL, aR) monotonised per cell (CW84 eq. 1.10)
   │             SL  = 2 (aR − a),  SLR = 2 (a − aL)             (effective slopes for the two faces)
   ├─ θ positivity guard on SL / SLR
   └─ flux kernel receives  UL = ql + ½ SL(i),  UR = qr − ½ SLR(i+1)

The limiters
============

With :math:`a = W_{i-1}`, :math:`b = W_i`, :math:`c = W_{i+1}` (the code's
``aa, bb, cc``), each function returns the limited slope over one cell:

.. list-table::
   :header-rows: 1
   :widths: 14 16 46 24

   * - ``limiterType``
     - function
     - formula
     - character
   * - 0
     - ``zslop``
     - :math:`0`
     - first order (piecewise constant)
   * - 1
     - ``vslop``
     - :math:`\dfrac{[{\rm sgn}(b-a) + {\rm sgn}(c-b)]\,|b-a|\,|c-b|}{|b-a| + |c-b| + 10^{-7}}` — van Leer's harmonic mean
     - smooth, moderately diffusive
   * - 2
     - ``fslop``
     - :math:`{\rm sgn}(c-a)\,\min\!\left(\tfrac{1}{2}|c-a|,\ 2|b-a|,\ 2|c-b|\right)` if :math:`(b-a)(c-b)>0`, else 0 — monotonized central (MC)
     - sharpest of the three TVD choices
   * - 3
     - ``minmod``
     - the smaller-magnitude of :math:`b-a` and :math:`c-b` if they agree in sign, else 0
     - most diffusive, and the safest at strong shocks

The choice matters more than one might expect: on the 2D field-loop test the
retained magnetic energy after eight crossings is 0.66 (minmod), 0.84 (van
Leer), 0.88 (MC) — and 0.93 with PPM. ``minmod`` is the safest at strong
shocks; for advection-dominated problems MC (``limiterType = 2``) is the
better PLM choice.

``zslop3`` and ``minmod3`` (interface ``limiter3``) are the variants for the
cylindrical coordinates (``solverIsoMHD2D``, ``solverIso2D``): they take the
radii of the three cells as well and apply the minmod limiter in the
logarithmic/uniform-:math:`r` metric of Skinner & Ostriker (2010).

PPM
===

``SCORPIO_PPM=1`` switches the 2D/3D Cartesian adiabatic and isothermal MHD
sweeps to the piecewise-parabolic method of Colella & Woodward (1984):
MC-limited differences, fourth-order face values, and the CW84
monotonisation of each cell's parabola. The parabola is passed to the
unchanged Riemann kernels as two *effective slopes* (``SL`` for the left
face, ``SLR`` for the right face), so the kernels reproduce the PPM face
states exactly. The normal field component keeps zero slope; the positivity
guard and the first-order flux correction apply to PPM unchanged.

PPM needs **three ghost cells** (``nbuf = 3``): the face value at
:math:`i+\tfrac{1}{2}` uses :math:`W_{i-1}\ldots W_{i+2}` and the neighbour's
parabola one more, so with two ghost cells the two ranks sharing a face would
build different parabolae — different fluxes on one face — and break
conservation and :math:`\nabla\cdot\boldsymbol{B}=0`. With ``nbuf = 2`` the
code prints one warning and runs PLM. The Colella & Sekora (2008)
smooth-extremum limiter was tried and not kept: no measurable accuracy gain
here (the linear-wave error is set by the RK2 dispersion) and its discontinuous
switch amplified round-off into an ``np``-dependence.

Not available: PPM in 1D, in the cylindrical coordinates, and (wired but not
validated) on the AMR mesh.
