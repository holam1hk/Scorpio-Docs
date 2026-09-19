.. _ch:hydro:

*************
Hydrodynamics
*************

Introduction
============

Governing equations, the code units, and the layout of the variable array
used in the code are presented here. Everything on this page is derived from
the current source (``src/gridModule.f03``, ``src/riemannSolverModule.f03``,
``src/hinnyCloud.f03``); where an older version of this page gave a different
number, the difference is pointed out explicitly.

-------------------------

Equations
=========

We begin with defining the conserved variables for MHD equations,
:math:`\boldsymbol{U} = (\rho, \rho \boldsymbol{u}, E, \boldsymbol{B})`,
where :math:`\rho` is the mass density, :math:`\boldsymbol{u}` is the velocity
vector, :math:`E` is the total energy density, and :math:`\boldsymbol{B}` is
the magnetic field vector.

Single-fluid (M)HD
------------------

.. math::

   \begin{align}
   \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \boldsymbol{u})&= 0 , \\
   \frac{\partial (\rho \boldsymbol{u})}{\partial t} + \nabla \cdot (\rho \boldsymbol{u} \boldsymbol{u} - \boldsymbol{B}\boldsymbol{B})
   + \nabla \left(P + \tfrac{1}{2}|\boldsymbol{B}|^2\right) &= - \rho \nabla \Phi, \\
   \frac{\partial E}{\partial t} + \nabla \cdot \left[ \left(E + P + \tfrac{1}{2}|\boldsymbol{B}|^2\right) \boldsymbol{u}
   - \boldsymbol{B}(\boldsymbol{B}\cdot\boldsymbol{u}) \right] &= - \rho \boldsymbol{u} \cdot \nabla \Phi, \\
   \frac{\partial \boldsymbol{B}}{\partial t} - \nabla \times (\boldsymbol{u} \times \boldsymbol{B}) &= 0, \qquad
   \nabla \cdot \boldsymbol{B} = 0 ,
   \end{align}

with :math:`E = \frac{P}{\Gamma-1} + \frac{1}{2}\rho|\boldsymbol{u}|^2 + \frac{1}{2}|\boldsymbol{B}|^2`.
Setting :math:`\boldsymbol{B}=0` gives the hydrodynamic (HD) system. The
gravitational potential :math:`\Phi` is either absent or obtained from the
self-gravity solver (:ref:`ch:selfgravity`), :math:`\nabla^2 \Phi = 4\pi G \rho`.

.. note::

   The magnetic pressure is written :math:`\tfrac{1}{2}|\boldsymbol{B}|^2`,
   i.e. the factor :math:`4\pi` of Gaussian units is absorbed into the code's
   magnetic field: :math:`B_{\rm code} = B_{\rm Gauss}/\sqrt{4\pi}`. This
   matters when converting field strengths, see the :ref:`derived units <tbl-derived-units>`.
   The face-centred storage of :math:`\boldsymbol{B}` and the constrained-transport
   update are described in :ref:`ch:methods`.

Equation of state (closure)
---------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - ``eosType``
     - pressure
     - notes
   * - 1 isothermal
     - :math:`P = c_s^2\,\rho`
     - ``sndspd`` :math:`= c_s` in km/s sets the temperature,
       :math:`T = \mu m_{\rm H} c_s^2 / k_B \simeq 282\,{\rm K}\,(\mu/2.33)\,(c_s/{\rm km\,s^{-1}})^2`.
       ``ene`` is carried but not used for the pressure.
   * - 2 adiabatic
     - :math:`P = (\Gamma-1)\left(E - \tfrac{1}{2}\rho|\boldsymbol{u}|^2 - \tfrac{1}{2}|\boldsymbol{B}|^2\right)`
     - ``adiGamma`` :math:`= \Gamma` (default 5/3). In 2D/3D Cartesian
       adiabatic MHD a *dual-energy* entropy variable
       :math:`\sigma = P/\rho^{\Gamma-1}` is advected as well and used to
       recover :math:`P` cancellation-free in low-:math:`\beta` cells
       (``SCORPIO_DUAL_ENERGY``, :ref:`ch:problem_file`).

Two-fluid ion-neutral (M)HD
---------------------------

For the two-fluid HD-MHD equations (neutrals :math:`n`, ions :math:`i`), we have:

.. math::

   \begin{align}
   \frac{\partial \rho_n}{\partial t} + \nabla \cdot (\rho_n \boldsymbol{u}_n)&= 0 , \\
   \frac{\partial (\rho_n \boldsymbol{u}_n)}{\partial t} + \nabla \cdot (\rho_n \boldsymbol{u}_n \boldsymbol{u}_n) + \nabla P_n +\rho_n \nabla \Phi
   &= - \alpha \rho_n \rho_i (\boldsymbol{u}_n - \boldsymbol{u}_i), \\
   \frac{\partial E_n}{\partial t} + \nabla \cdot (\boldsymbol{u}_n E_n + P_n \boldsymbol{u}_n)
   &= - \rho_n \boldsymbol{u}_n \cdot \nabla \Phi + \frac{3 \alpha}{\mu_n + \mu_i} (\mu_i \rho_n P_i - \mu_n \rho_i P_n) + \frac{\alpha \mu_i}{\mu_n + \mu_i} \rho_n \rho_i (\boldsymbol{u}_n - \boldsymbol{u}_i)^2, \\
   \frac{\partial \rho_i}{\partial t} + \nabla \cdot (\rho_i \boldsymbol{u}_i) &= 0 , \\
   \frac{\partial (\rho_i \boldsymbol{u}_i)}{\partial t} + \nabla \cdot ( \rho_i \boldsymbol{u}_i \boldsymbol{u}_i - \boldsymbol{B} \boldsymbol{B} )
   + \nabla \left(P_i + \tfrac{1}{2} |\boldsymbol{B}|^2 \right) + \rho_i \nabla \Phi
   &= - \alpha \rho_n \rho_i (\boldsymbol{u}_i - \boldsymbol{u}_n), \\
   \frac{\partial E_i}{\partial t} + \nabla \cdot \left[\left(E_i + P_i +\tfrac{1}{2} |\boldsymbol{B}|^2\right) \boldsymbol{u}_i - \boldsymbol{B}(\boldsymbol{B}\cdot\boldsymbol{u}_i)\right]
   &= - \rho_i \boldsymbol{u}_i \cdot \nabla \Phi + \frac{3 \alpha}{\mu_n + \mu_i} (\mu_n \rho_i P_n - \mu_i \rho_n P_i) + \frac{\alpha \mu_i}{\mu_n + \mu_i} \rho_n \rho_i (\boldsymbol{u}_n - \boldsymbol{u}_i)^2,\\
   \frac{\partial \boldsymbol{B}}{\partial t} - \nabla \times (\boldsymbol{u}_i \times \boldsymbol{B}) &= 0,\\
   \nabla \cdot \boldsymbol{B} &= 0.
   \end{align}

In the above formulas :math:`\rho_n, \mu_n, P_n, \boldsymbol{u}_n, \Gamma_n`,
and :math:`E_n=\frac{P_n}{\Gamma_n-1}+\frac{1}{2}\rho_n|\boldsymbol{u}_n|^2`
respectively denote the mass density, molecular weight, thermal pressure,
velocity vector, specific heat ratio, and total energy density of neutrals,
while :math:`\rho_i, \mu_i, P_i, \boldsymbol{u}_i, \Gamma_i`, and
:math:`E_i=\frac{P_i}{\Gamma_i-1}+\frac{1}{2}\rho_i|\boldsymbol{u}_i|^2 + \frac{1}{2}|\boldsymbol{B}|^2`
are the corresponding variables of the ions, which carry the magnetic field.
The collision coefficient is
:math:`\alpha_0=\frac{1.9\times10^{-19}}{m_n+m_i}\ {\rm cm^3\,s^{-1}}`
(where :math:`m_n=\mu_n m_H`, :math:`m_i=\mu_i m_H`). For a molecular cloud
with :math:`\mu_n=2.3` and :math:`\mu_i=29` this gives
:math:`\alpha = 3.7\times10^{13}\ {\rm cm^3\,s^{-1}\,g^{-1}}`, i.e.
:math:`\alpha \simeq 7.7\times 10^{4}` in code units (:ref:`derived units <tbl-derived-units>`).

In the code the two fluids are two ``grid`` objects (``gn``, ``gi``) that are
evolved together; the drag term is integrated either operator-split (TR-BDF2)
or inside the RK stages (IMEX), selected with ``SCORPIO_AD_SCHEME``
(:ref:`ch:problem_file`, :ref:`ch:AD`).

Again, the finite-volume method always solves for the conserved variables
(density, momenta, total energy, and the left/right interface magnetic
fields). Primitive variables are recovered from them (:ref:`sec:primitive`).
The conserved variables are stored in the array ``q(i,j,k,nvar)``, where
``i,j,k`` are the cell indices and the last index selects the variable
(:ref:`sec:conserved`).


.. _sec:units:

Code units and conversions
==========================

The code is dimensionless; the units below are the *molecular-cloud* set that
all built-in cases and the constant ``GravConst`` assume. Three base units are
chosen and everything else follows from them:

.. rubric:: Base code units

.. list-table::
   :header-rows: 1
   :widths: 14 22 22 22 20

   * - dimension
     - code unit
     - cgs
     - SI
     - useful value
   * - length :math:`[L]`
     - 1 pc
     - :math:`3.0857\times10^{18}` cm
     - :math:`3.0857\times10^{16}` m
     - :math:`2.063\times10^{5}` AU
   * - mass :math:`[M]`
     - 1 :math:`M_\odot`
     - :math:`1.989\times10^{33}` g
     - :math:`1.989\times10^{30}` kg
     -
   * - velocity :math:`[V]`
     - 1 km/s
     - :math:`10^{5}` cm/s
     - :math:`10^{3}` m/s
     -
   * - time :math:`[T]=[L]/[V]`
     - pc / (km/s)
     - :math:`3.0857\times10^{13}` s
     - :math:`3.0857\times10^{13}` s
     - :math:`9.778\times10^{5}` yr = 0.978 Myr

So ``t = 1.0`` in the code is 0.978 Myr, and 1 Myr is ``t = 1.0227``.

.. rubric:: Physical constants used for the conversions

.. list-table::
   :header-rows: 1
   :widths: 22 18 30 30

   * - constant
     - symbol
     - cgs
     - in code units
   * - gravitational constant (``GravConst`` in the code)
     - :math:`G`
     - :math:`6.674\times10^{-8}` cm\ :sup:`3` g\ :sup:`-1` s\ :sup:`-2`
     - **4.3011e-3** pc km\ :sup:`2` s\ :sup:`-2` :math:`M_\odot^{-1}`
   * - parsec
     - pc
     - :math:`3.0857\times10^{18}` cm
     - 1
   * - solar mass
     - :math:`M_\odot`
     - :math:`1.989\times10^{33}` g
     - 1
   * - year
     - yr
     - :math:`3.156\times10^{7}` s
     - :math:`1.0227\times10^{-6}`
   * - hydrogen mass
     - :math:`m_{\rm H}`
     - :math:`1.6735\times10^{-24}` g
     - :math:`8.41\times10^{-58}`
   * - Boltzmann constant
     - :math:`k_B`
     - :math:`1.3807\times10^{-16}` erg/K
     - (use the :math:`T`–:math:`c_s` relation below)

.. note::

   An older version of this page (and the linked gist) used
   :math:`M_\odot = 1.9\times10^{33}` g, which gives
   :math:`[\rho] = 6.465\times10^{-23}` g cm\ :sup:`-3` and
   :math:`[B] = 8.0405\times10^{-7}` (the value quoted in a comment in
   ``hinnyCloud.f03``). The value of ``GravConst`` hard-coded in the code,
   4.3011e-3, corresponds to :math:`M_\odot = 1.989\times10^{33}` g, so that is
   the mass unit used consistently in the tables here. The difference is 4.5 %
   in every quantity that contains the mass unit.

.. _tbl-derived-units:

.. rubric:: Derived units of the common variables

.. list-table::
   :header-rows: 1
   :widths: 20 20 22 22 16

   * - quantity (code variable)
     - code unit
     - cgs
     - SI
     - notes
   * - density ``den``, :math:`\rho`
     - :math:`M_\odot` pc\ :sup:`-3`
     - :math:`6.770\times10^{-23}` g cm\ :sup:`-3`
     - :math:`6.770\times10^{-20}` kg m\ :sup:`-3`
     - 17.4 particles cm\ :sup:`-3` for :math:`\mu=2.33`; 14.5 H\ :sub:`2` cm\ :sup:`-3`; :math:`n_{\rm H}=28.7` cm\ :sup:`-3`
   * - column density :math:`\Sigma=\int\rho\,dl`
     - :math:`M_\odot` pc\ :sup:`-2`
     - :math:`2.089\times10^{-4}` g cm\ :sup:`-2`
     - :math:`2.089\times10^{-3}` kg m\ :sup:`-2`
     - :math:`N_{\rm H_2} = 4.46\times10^{19}` cm\ :sup:`-2`
   * - velocity, :math:`\boldsymbol{u}` = ``mom``/``den``
     - km s\ :sup:`-1`
     - :math:`10^{5}` cm s\ :sup:`-1`
     - :math:`10^{3}` m s\ :sup:`-1`
     -
   * - momentum density ``momx`` …, :math:`\rho\boldsymbol{u}`
     - :math:`M_\odot` pc\ :sup:`-3` km s\ :sup:`-1`
     - :math:`6.770\times10^{-18}` g cm\ :sup:`-2` s\ :sup:`-1`
     - :math:`6.770\times10^{-17}` kg m\ :sup:`-2` s\ :sup:`-1`
     -
   * - pressure, energy density ``ene``, :math:`P,\ E`
     - :math:`M_\odot` pc\ :sup:`-3` km\ :sup:`2` s\ :sup:`-2`
     - :math:`6.770\times10^{-13}` erg cm\ :sup:`-3` (dyn cm\ :sup:`-2`)
     - :math:`6.770\times10^{-14}` Pa
     - :math:`P/k_B = 4.90\times10^{3}` K cm\ :sup:`-3`
   * - energy (integrated), e.g. ``E_turb_tot``
     - :math:`M_\odot` km\ :sup:`2` s\ :sup:`-2`
     - :math:`1.989\times10^{43}` erg
     - :math:`1.989\times10^{36}` J
     -
   * - specific energy, potential ``gphi``, :math:`\Phi`
     - km\ :sup:`2` s\ :sup:`-2`
     - :math:`10^{10}` erg g\ :sup:`-1`
     - :math:`10^{6}` J kg\ :sup:`-1`
     -
   * - acceleration ``sgfx`` …, :math:`\boldsymbol{g}=-\nabla\Phi`
     - km s\ :sup:`-1` / [T]
     - :math:`3.241\times10^{-9}` cm s\ :sup:`-2`
     - :math:`3.241\times10^{-11}` m s\ :sup:`-2`
     -
   * - magnetic field ``bxl`` …, :math:`\boldsymbol{B}`
     - :math:`\sqrt{M_\odot\,{\rm pc^{-3}}}` km s\ :sup:`-1`
     - :math:`8.228\times10^{-7}` :math:`\sqrt{\rm g\,cm^{-1}}` s\ :sup:`-1` (Heaviside–Lorentz)
       = **2.917 μG** (Gaussian, :math:`B_{\rm Gauss}=\sqrt{4\pi}\,B_{\rm code}`)
     - :math:`2.917\times10^{-10}` T
     - e.g. ``b0 = 35.08`` in the cloud case is 102 μG
   * - temperature (isothermal, via ``sndspd``)
     - :math:`T = \mu m_{\rm H} c_s^2/k_B`
     - 282 K × :math:`(\mu/2.33)(c_s/{\rm km\,s^{-1}})^2`
     - same
     - :math:`c_s = 0.2` → 11 K; :math:`0.3` → 25 K
   * - gravitational constant ``GravConst``
     - pc km\ :sup:`2` s\ :sup:`-2` :math:`M_\odot^{-1}`
     - :math:`6.674\times10^{-8}` cgs
     - :math:`6.674\times10^{-11}` SI
     - 4.3011e-3 in code units
   * - ion–neutral coupling ``alpha_ad``, :math:`\alpha`
     - pc\ :sup:`2` km s\ :sup:`-1` :math:`M_\odot^{-1}`
     - :math:`4.787\times10^{8}` cm\ :sup:`3` g\ :sup:`-1` s\ :sup:`-1`
     - :math:`4.787\times10^{5}` m\ :sup:`3` kg\ :sup:`-1` s\ :sup:`-1`
     - :math:`3.7\times10^{13}` cgs = :math:`7.73\times10^{4}` code
   * - time ``t``, ``dt``, ``tend``, ``dtout``
     - pc / (km/s)
     - :math:`3.0857\times10^{13}` s
     - same
     - 0.978 Myr

Handy formulas in code units (:math:`G` = ``GravConst``, :math:`\rho` in
:math:`M_\odot` pc\ :sup:`-3`, :math:`c_s` in km/s):

.. math::

   \begin{align}
   t_{\rm ff} &= \sqrt{\frac{3\pi}{32 G \rho}} = \frac{8.275}{\sqrt{\rho}}\ [T] = \frac{8.09\ {\rm Myr}}{\sqrt{\rho}}, \\
   \lambda_J &= c_s\sqrt{\frac{\pi}{G \rho}} = 27.0\,\frac{c_s}{\sqrt{\rho}}\ {\rm pc}, \\
   v_A &= \frac{|\boldsymbol{B}|}{\sqrt{\rho}}\ {\rm km\,s^{-1}}, \qquad
   \beta = \frac{2P}{|\boldsymbol{B}|^2} = \frac{2 c_s^2 \rho}{|\boldsymbol{B}|^2}\ ({\rm isothermal}).
   \end{align}

The Truelove stop criterion used by the cloud driver
(``TrueloveCondition``) compares :math:`\lambda_J` at the global maximum
density with the cell size and stops the run when :math:`\lambda_J < 5\,\Delta x`;
the AMR driver instead refines where :math:`\lambda_J < ` ``amr_jeans_n``
:math:`\Delta x`.


.. _sec:conserved:

Conserved variables
===================

Which variables are evolved is chosen with the flag array ``variable(1:8)``
(1 = evolve, 0 = absent) in the problem driver. The storage order in
``q(i,j,k,:)`` — and in the HDF5 output — is fixed:

.. list-table::
   :header-rows: 1
   :widths: 14 14 14 12 22 24

   * - ``variable()``
     - dataset
     - quantity
     - slot (MHD)
     - slot (HD, ``variable(5:7)=0``)
     - units
   * - ``variable(1)``
     - ``den``
     - :math:`\rho`
     - 1
     - 1
     - :math:`M_\odot` pc\ :sup:`-3`
   * - ``variable(2)``
     - ``momx``
     - :math:`\rho u_x`
     - 2
     - 2
     - :math:`M_\odot` pc\ :sup:`-3` km s\ :sup:`-1`
   * - ``variable(3)``
     - ``momy``
     - :math:`\rho u_y`
     - 3
     - 3
     -
   * - ``variable(4)``
     - ``momz``
     - :math:`\rho u_z`
     - 4
     - 4
     -
   * - ``variable(5)``
     - ``bxl``
     - :math:`B_x` on the left (:math:`-x`) cell face
     - 5
     - —
     - :math:`\sqrt{M_\odot\,{\rm pc^{-3}}}` km s\ :sup:`-1`
   * - ``variable(6)``
     - ``byl``
     - :math:`B_y` on the left (:math:`-y`) face
     - 6
     - —
     -
   * - ``variable(7)``
     - ``bzl``
     - :math:`B_z` on the left (:math:`-z`) face
     - 7
     - —
     -
   * - ``variable(8)``
     - ``ene``
     - :math:`E` (total energy density)
     - 8
     - 5
     - :math:`M_\odot` pc\ :sup:`-3` km\ :sup:`2` s\ :sup:`-2`
   * - (with ``variable(5)``)
     - ``bxr``
     - :math:`B_x` on the right (:math:`+x`) face
     - 9
     - —
     -
   * - (with ``variable(6)``)
     - ``byr``
     - :math:`B_y` on the right face
     - 10
     - —
     -
   * - (with ``variable(7)``)
     - ``bzr``
     - :math:`B_z` on the right face
     - 11
     - —
     -
   * - (automatic)
     - not written
     - :math:`\sigma = P/\rho^{\Gamma-1}` dual-energy entropy
     - 12
     - —
     - only 2D/3D Cartesian with ``variable(8)=variable(5)=1``

The number of stored variables is

.. code-block:: fortran

   nvar = den + vx + vy + vz + 2*(bx + by + bz) + ene   ! + 1 for the dual-energy slot

so a full 3D MHD run has ``nvar = 11`` (12 with dual energy), an HD run with
energy has ``nvar = 5``, and an isothermal HD run without ``ene`` has
``nvar = 4``. In the two-fluid case the neutral grid is HD (``ene`` in slot 5)
and the ion grid is MHD (``ene`` in slot 8), which is what the older
"8 (ions) / 5 (neutrals)" note meant.

- Each magnetic-field component is stored twice, on the two cell faces normal
  to it, because the constrained-transport (CT) update keeps
  :math:`\nabla\cdot\boldsymbol{B}=0` to machine precision on the faces. The
  cell-centred field is the average, e.g. :math:`B_x = \tfrac{1}{2}(` ``bxl`` + ``bxr`` :math:`)`.
  When you write an initial condition, set both faces (``q(i,j,k,5)`` and
  ``q(i,j,k,9)`` for :math:`B_x`, etc.); for a uniform field they are equal.
- For isothermal MHD ``variable(8)`` must still be 1 (the solver layout expects
  the slot); the value is not used for the pressure.
- The indices ``i,j,k`` run from ``1-nbuf`` to ``nMesh+nbuf``: ``nbuf`` ghost
  cells (default 2, 3 with PPM reconstruction) pad every side. The HDF5 output
  writes the arrays *including* the ghost cells, so strip ``nbuf`` cells from
  each end when you analyse them.

.. _sec:primitive:

Primitive variables
===================

.. math::

   \begin{align}
   \boldsymbol{u} &= \frac{(\texttt{momx},\texttt{momy},\texttt{momz})}{\texttt{den}}, \qquad
   \boldsymbol{B} = \tfrac{1}{2}\left(\boldsymbol{B}_{\rm left} + \boldsymbol{B}_{\rm right}\right), \\
   P &= c_s^2\,\rho \quad (\texttt{eosType}=1), \qquad
   P = (\Gamma-1)\left(\texttt{ene} - \tfrac{1}{2}\rho|\boldsymbol{u}|^2 - \tfrac{1}{2}|\boldsymbol{B}|^2\right) \quad (\texttt{eosType}=2), \\
   T &= \frac{\mu m_{\rm H}}{k_B}\frac{P}{\rho} = 282\ {\rm K}\,\frac{\mu}{2.33}\,\frac{P}{\rho}\Big|_{\rm code}.
   \end{align}

A minimal reader for a 3D snapshot (arrays come out of ``h5py`` in ``(z, y, x)``
order):

.. code-block:: python

   import h5py, numpy as np
   f  = h5py.File("g0800_0010.h5", "r")
   nx, ny, nz = [int(v) for v in f["nMesh"][()]]
   nb = int(np.asarray(f["nbuf"][()]).reshape(-1)[0])
   cube = lambda name: np.transpose(f[name][()])[nb:nx+nb, nb:ny+nb, nb:nz+nb]   # -> (x,y,z), no ghosts
   rho = cube("den")                                  # Msun / pc^3
   vx  = cube("momx") / rho                           # km / s
   bz  = 0.5 * (cube("bzl") + cube("bzr"))            # code units; x 2.917 for microgauss
   rho_cgs = rho * 6.770e-23                          # g / cm^3

.. note:: `Gist for more useful units and constants <https://gist.github.com/S-Yuan137/33a1489bfc5d697e0748b76e0228fdf8>`_
   (it uses :math:`M_\odot = 1.9\times10^{33}` g; see the note above).
