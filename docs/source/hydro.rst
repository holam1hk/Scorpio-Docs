.. _ch:hydro:

*************
Hydrodynamics
*************

Introduction
============

Governing equations and variables used in the code are presented here.

-------------------------
   
Equations 
==================
We begin with defining the conserved variables for MHD equations,
:math:`\boldsymbol{U} = (\rho, \rho \boldsymbol{u}, \rho E, \boldsymbol{B}),`
where :math:`\rho` is the mass density, :math:`\boldsymbol{u}` is the velocity vector, :math:`E` is the total energy density, and :math:`\boldsymbol{B}` is the magnetic field vector.
For single-fluid MHD equations, we have:

.. math::

   \begin{align}
   \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \boldsymbol{u})&= 0 , \\
   \frac{\partial (\rho \boldsymbol{u})}{\partial t} + \nabla \cdot (\rho \boldsymbol{u} \boldsymbol{u}) + \nabla P &= - \rho \nabla \Phi, \\
   \frac{\partial (E)}{\partial t} + \nabla \cdot (\boldsymbol{u} E + P \boldsymbol{u}) &= - \rho \boldsymbol{u} \nabla \Phi. \label{eq:compressible-equations}
   \end{align}

For two-fluid HD-MHD equations, we have:

.. math::

   \begin{align}
   \frac{\partial \rho_n}{\partial t} + \nabla \cdot (\rho_n \boldsymbol{u}_n)&= 0 , \\
   \frac{\partial (\rho_n \boldsymbol{u}_n)}{\partial t} + \nabla \cdot (\rho_n \boldsymbol{u}_n \boldsymbol{u}_n) + \nabla P_n +\rho_n \nabla \Phi 
   &= - \alpha \rho_n \rho_i (\boldsymbol{v}_n - \boldsymbol{v}_i), \\
   \frac{\partial E_n}{\partial t} + \nabla \cdot (\boldsymbol{u}_n E_n + P_n \boldsymbol{u}_n) 
   &= - \rho_n \boldsymbol{u}_n \nabla \Phi + \frac{3 \alpha}{\mu_n + \mu_i} (\mu_i \rho_n P_i - \mu_n \rho_i P_n) + \frac{\alpha \mu_i}{\mu_n + \mu_i} \rho_n \rho_i (\boldsymbol{v}_n - \boldsymbol{v}_i)^2, \\
   \frac{\partial \rho_i}{\partial t} + \nabla \cdot (\rho_i \boldsymbol{u}_i) &= 0 , \\
   \frac{\partial (\rho_i \boldsymbol{u}_i)}{\partial t} + \nabla \cdot ( (\rho_i \boldsymbol{u}_i \boldsymbol{u}_i) - \boldsymbol{B} \boldsymbol{B} )
   + \nabla (P_i + \frac{1}{2} |\boldsymbol{B}|^2 ) + \rho_i \nabla \Phi 
   &= - \alpha \rho_n \rho_i (\boldsymbol{v}_i - \boldsymbol{v}_n), \\
   \frac{\partial E_i}{\partial t} + \nabla \cdot [(E_i + P_i +\frac{1}{2} |\boldsymbol{B}|^2) \boldsymbol{u}_i ] 
   &= - \rho_i \boldsymbol{u}_i \nabla \Phi + \frac{3 \alpha}{\mu_n + \mu_i} (\mu_n \rho_i P_n - \mu_i \rho_n P_i) + \frac{\alpha \mu_i}{\mu_n + \mu_i} \rho_n \rho_i (\boldsymbol{v}_n - \boldsymbol{v}_i)^2,\\
   \frac{\partial \boldsymbol{B}}{\partial t} + \nabla \times (v_i \times \boldsymbol{B}) &= 0,\\
   \nabla \boldsymbol{B} &= 0. \\  
   \label{eq:admhd-equations}
   \end{align}

In the above formulas :math:`\rho_n, \mu_n, P_n, \boldsymbol{v_n}, \Gamma_n`, and :math:`E_n=\frac{P_n}{\Gamma_n-1}+\frac{1}{2}\rho_n|\boldsymbol{u_n}|^2` respectively denote the mass density, molecular weight, thermal pressure,
velocity vector, specific heat ratio, and total energy density of neutrals. While :math:`\rho_i, \mu_i, P_i, \boldsymbol{v_i}, \Gamma_i`, and :math:`E_i=\frac{P_i}{\Gamma_i-1}+\frac{1}{2}\rho_i|\boldsymbol{v_i}|^2 + \frac{1}{2}|\boldsymbol{B}|^2` are the corresponding physical variables of ions, where :math:`\boldsymbol{B}` is the magnetic field.  
:math:`\alpha=\alpha_0 max\left(1,\frac{|\boldsymbol{v_i}-\boldsymbol{v_n}|}{\boldsymbol{v_{\alpha}}}\right)` is the collision coefficient, in which :math:`\alpha_0=\frac{1.9\times10^{-19}}{m_n+m_i} cm^3 s^{-1}` (where :math:`m_n=\mu_n m_H, m_i=\mu_i m_H` with :math:`m_H` being the mass of hydrogen atom). In molecular cloud, we assume :math:`\mu_n=2.3` and :math:`\mu_i= 29` and the
corresponding collision coefficient :math:`\alpha = 3.7\times10^{13} cm^3 s^{-1} g^{-1}`. 

Again, FVM always solves for conserved variables (density, momenta, total energy, and left-and-right interface magnetic fields). Primitive variables can be obtained easily.
The conserved variables are stored in the variable array :math:`q(i,j,k,nvar)`, where :math:`i,j,k` are the cell indices and :math:`nvar` is the index of the variable.

   
   
Useful constants and units
=============================

.. _table:Useful constants:
.. table:: Useful constants:
    
   +-----------------------+-----------------------+-----------+-------------------------------+
   | **variable**          | **quantity**          | **value** |  **units**                    |
   +=======================+=======================+===========+===============================+
   | ``GravConst``         | :math:`G`             | 4.3011d-3 | km^2*pc/(s^2*Msun)            |
   +-----------------------+-----------------------+-----------+-------------------------------+
   | ``year``              | :math:`year`          | 3.156e7   | s                             |
   +-----------------------+-----------------------+-----------+-------------------------------+
   | ``pc``                | :math:`pc`            | 3.086e18  | cm                            |
   +-----------------------+-----------------------+-----------+-------------------------------+
   | ``Msun``              | :math:`M_{sun}`       | 1.9e33    | g                             |
   +-----------------------+-----------------------+-----------+-------------------------------+
   |                       |                       |           |                               |
   +-----------------------+-----------------------+-----------+-------------------------------+

The code units are taken for molecular cloud simulations. 

.. _table:Code units:
.. table:: code units:

   +-----------------------+-------------------------------------------+-------------------+-------------------------------+
   | **Dimension**         | **Code Unit**                             | **Useful Value**  | **CGS Equivalent**            |
   +=======================+===========================================+===================+===============================+
   | ``[T]``               | pc * s / km                               | 9.78e5 yrs        | 3.086e13 ms                   |
   +-----------------------+-------------------------------------------+-------------------+-------------------------------+
   | ``[L]``               | pc                                        | 3.086e13 km       | 3.086e18 cm                   |
   +-----------------------+-------------------------------------------+-------------------+-------------------------------+
   | ``[M]``               | Solar mass :math:`M_{sun}`                |                   | 1.9e33 g                      |
   +-----------------------+-------------------------------------------+-------------------+-------------------------------+

.. _table:Variables:
.. table:: Common Variables:

   +-----------------------+-------------------------------------------+-------------------+--------------------------------+
   | **Variables**         | **Code Unit**                             | **Useful Value**  | **CGS Equivalent**             |
   +=======================+===========================================+===================+================================+
   | ``rho``               | :math:`M_{sun}` / pc:math:`^3`            |                   | 6.46496e-23 g cm :math:`^-3`   |
   +-----------------------+-------------------------------------------+-------------------+--------------------------------+
   | ``v``                 | km/s                                      |                   | 10 :math:`^5` cm/s             |
   +-----------------------+-------------------------------------------+-------------------+--------------------------------+
   | ``B``                 | :math:`B_{code}=B_{phy}/\sqrt{4\pi}`      |                   | 8.0405e-7 :math:`\sqrt{g/cm}`/s|
   +-----------------------+-------------------------------------------+-------------------+--------------------------------+
   | ``alpha``             | pc:math:`^2` * km / (s * M)               |                   | 5.012314e8 cm :math:`^3`/g/s   |
   +-----------------------+-------------------------------------------+-------------------+--------------------------------+

..  +-----------------------+-------------------------------------------+-------------------+-------------------------------+ 
   | ``tff``               | :math:`\sqrt{3\pi/(32G\rho)}`             |                   | free-fall time                |

Conserved variables
===================

.. _table:Conserved variables:
.. table:: Conserved variables:
   
   +-----------------------+-----------------------+-----------------------+-------------------------------+
   | **variable()**        | **quantity**          | **nvar**              |  **units**                    |
   +=======================+=======================+=======================+===============================+
   | ``density``           | :math:`\rho`          | 1                     | [M_sun / pc^3]                |
   +-----------------------+-----------------------+-----------------------+-------------------------------+
   | ``momx``              | :math:`\rho v_x`      | 2                     | [M_sun / pc^3 * km/s]         |
   +-----------------------+-----------------------+-----------------------+-------------------------------+
   | ``momy``              | :math:`\rho v_y`      | 3                     | [M_sun / pc^3 * km/s]         |
   +-----------------------+-----------------------+-----------------------+-------------------------------+
   | ``momz``              | :math:`\rho v_z`      | 4                     | [M_sun / pc^3 * km/s]         |
   +-----------------------+-----------------------+-----------------------+-------------------------------+
   | ``bx``                | :math:`B_x`           | 5                     | [sqrt(M_sun / pc^3) * km/s]   |
   +-----------------------+-----------------------+-----------------------+-------------------------------+
   | ``by``                | :math:`B_y`           | 6                     | [sqrt(M_sun / pc^3) * km/s]   |
   +-----------------------+-----------------------+-----------------------+-------------------------------+
   | ``bz``                | :math:`B_z`           | 7                     | [sqrt(M_sun / pc^3) * km/s]   |
   +-----------------------+-----------------------+-----------------------+-------------------------------+
   | ``ene``               | :math:`E`             | 8 (ions) 5 (neutrals) | [M_sun / pc^3 * km^2/s^2]     |
   +-----------------------+-----------------------+-----------------------+-------------------------------+

-  ``ene`` is only required for adiabatic EOS, but not for isothermal EOS. 
-  The indices of variable q(i,j,k,nvar) counts from 1-nbuf to nMesh+nbuf. The last index refers to different variables in :numref:`table:Conserved variables`. 

Numerical setup
=================

.. _table:Coordinate Type:
.. table:: coordType:
   
   +---------------------------+-----------+
   | **Coordinate System**     | **Values**|
   +===========================+===========+
   | ``Cartesian``             | 1         |                               
   +---------------------------+-----------+
   | ``Cylindrical log r``     | 2         |                              
   +---------------------------+-----------+
   | ``Cylindrical uniform r`` | 3         |                              
   +---------------------------+-----------+

.. _table:Boundary Type:
.. table:: boundaryType:
   
   +---------------------------+-----------+
   | **Boundary Condition**    | **Values**|
   +---------------------------+-----------+
   | ``user defined``          | 0         | 
   +===========================+===========+
   | ``zero gradient``         | 1         |                               
   +---------------------------+-----------+
   | ``reflective``            | 2         |                              
   +---------------------------+-----------+
   | ``periodic``              | 3         |                              
   +---------------------------+-----------+

.. _table:solverType:
.. table:: solverType:
   
   +---------------------------+-----------+
   | **Riemann Solver**        | **Values**|
   +===========================+===========+
   | ``exactHD``               | 1         |                               
   +---------------------------+-----------+
   | ``HLL`` for HD            | 2         |                              
   +---------------------------+-----------+
   | ``HLLC``                  | 3         |                              
   +---------------------------+-----------+
   | ``HLL`` for MHD           | 4         |                              
   +---------------------------+-----------+
   | ``HLLD``                  | 5         |                              
   +---------------------------+-----------+
   
.. _table:limiterType:
.. table:: limiterType:
   
   +---------------------------+-----------+
   | **Slope Limiter**         | **Values**|
   +===========================+===========+
   | ``zero``                  | 0         |                               
   +---------------------------+-----------+
   | ``van Leer``              | 1         |                              
   +---------------------------+-----------+
   | ``fslop``                 | 2         |                              
   +---------------------------+-----------+  
   | ``minmod``                | 3         |                              
   +---------------------------+-----------+
   
.. _table:eosType:
.. table:: eosType:
   
   +---------------------------+-----------+
   | **Equation of State**     | **Values**|
   +===========================+===========+
   | ``isothermal``            | 1         |                               
   +---------------------------+-----------+
   | ``adiabatic``             | 2         |                              
   +---------------------------+-----------+

.. note::
   Not all EOS and solvers are supported.


.. _table:Problem setup:
.. table:: Problem setup:
   
   +-----------------------------+----------------------+-------------------------------------------------------------+
   | **Variables**               | **Values**           | **notes**                                                   |
   +=============================+======================+=============================================================+
   | ``gridID``                  | integer              | Problem ID                                                  |
   +-----------------------------+----------------------+-------------------------------------------------------------+
   | ``ndim``                    | 1 / 2 / 3            | Number of dimensions                                        |
   +-----------------------------+----------------------+-------------------------------------------------------------+
   | ``nbuf``                    | 2 (default)          | Number of buffer zones                                      |
   +-----------------------------+----------------------+-------------------------------------------------------------+
   | ``variable``                | 0 (default)          | Variables in :numref:`table:variables`, set 1 to initialize |
   +-----------------------------+----------------------+-------------------------------------------------------------+
   | ``nMesh(1)/(2)/(3)``        | integer              | Number of meshes in each dimension                          |
   +-----------------------------+----------------------+-------------------------------------------------------------+
   | ``leftBdry(1)/(2)/(3)``     | double precision     | Boundary position in each dimension                         |
   +-----------------------------+----------------------+-------------------------------------------------------------+
   | ``rightBdry(1)/(2)/(3)``    | double precision     | Boundary position in each dimension                         |
   +-----------------------------+----------------------+-------------------------------------------------------------+
   | ``sndspd``                  | double precision     | Soundspeed [km/s]                                           |
   +-----------------------------+----------------------+-------------------------------------------------------------+
   | ``CFL``                     | 0 < CFL < 1          | CFL condition                                               |
   +-----------------------------+----------------------+-------------------------------------------------------------+
   | ``gam`` or ``adiGamma``     | 5.d0/3.d0 (default)  | Useless if isothermal                                       |
   +-----------------------------+----------------------+-------------------------------------------------------------+
   | ``nstep``                   | integer              | Number of time steps                                        |
   +-----------------------------+----------------------+-------------------------------------------------------------+
   | ``t``                       | double precision     | Current time in [code unit]                                 |
   +-----------------------------+----------------------+-------------------------------------------------------------+
   | ``dt``                      | double precision     | Time step in [code unit]                                    |
   +-----------------------------+----------------------+-------------------------------------------------------------+
   | ``time_end``                | double precision     | Simulation end time in [code unit]                          |
   +-----------------------------+----------------------+-------------------------------------------------------------+


File output
=================

.. _table:File output:
.. table:: File output:

   +-----------------------------+----------------------+-------------------------------------------------------------+
   | **Variables**               | **Values**           | **notes**                                                   |
   +=============================+======================+=============================================================+
   | ``file_start``              | 0                    | start file number                                           |
   +-----------------------------+----------------------+-------------------------------------------------------------+
   | ``dt_out``                  | double precision     | output interval in code unit                                |
   +-----------------------------+----------------------+-------------------------------------------------------------+
   | ``toutput``                 | double precision     | output time in code unit                                    |
   +-----------------------------+----------------------+-------------------------------------------------------------+
   | ``fnum``                    | integer              | outputfile number                                           |
   +-----------------------------+----------------------+-------------------------------------------------------------+
   | ``write_vtk``               | .true./.false.       | vtk format output                                           |
   +-----------------------------+----------------------+-------------------------------------------------------------+


MPI setup
=================

.. _table:MPI setup:
.. table:: MPI setup:

   +-----------------------------+----------------------+-------------------------------------------------------------+
   | **Variables**               | **Values**           | **notes**                                                   |
   +=============================+======================+=============================================================+
   | ``periods(1)/(2)/(3)``      | .true./.false.       | MPI periodic boundary conditions                            |
   +-----------------------------+----------------------+-------------------------------------------------------------+
   | ``dims(1)/(2)/(3)``         | integer              | Number of processes (CPU cores) in each dimension           |
   +-----------------------------+----------------------+-------------------------------------------------------------+
   | ``reorder``                 | .true./.false.       | MPI reorder flag                                            |
   +-----------------------------+----------------------+-------------------------------------------------------------+
   | ``myid``                    | integer              | MPI rank                                                    |
   +-----------------------------+----------------------+-------------------------------------------------------------+
   | ``nprocs``                  | integer              | Total number of MPI processes                               |
   +-----------------------------+----------------------+-------------------------------------------------------------+
   | ``MPI_COMM_WORLD``          | integer              | MPI communicator                                            |
   +-----------------------------+----------------------+-------------------------------------------------------------+
   | ``ierr``                    | integer              | MPI error code                                              |
   +-----------------------------+----------------------+-------------------------------------------------------------+


Physics Modules   
=================
.. _table:TurbulenceDriving:
.. table:: TurbulenceDriving:

   +---------------------------+------------------+------------------------------------+
   | **Variables**             | **Values**       | **notes**                          |
   +===========================+==================+====================================+
   | ``DriveTurbulence``       | .true./.false.   | Drive turbulence or not            |
   +---------------------------+------------------+------------------------------------+
   | ``DT_mode``               | 0                | Drive at begining                  |            
   +---------------------------+------------------+------------------------------------+
   | ``DT_mode``               | 1                | Drive periodically                 |             
   +---------------------------+------------------+------------------------------------+
   | ``dt_turb``               | double precision | Energy injection time interval     |        
   +---------------------------+------------------+------------------------------------+
   | ``t_count_turb``          | double precision | Energy injection count time        |
   +---------------------------+------------------+------------------------------------+
   | ``t_accum_turb``          | double precision | Energy injection accumulation time |                                 
   +---------------------------+------------------+------------------------------------+
   | ``zeta``                  | 0<=zeta<=1       | Soledoidal or compressive driving  |       
   +---------------------------+------------------+------------------------------------+
   | ``DT_scale``              | double precision | Driving scale                      |           
   +---------------------------+------------------+------------------------------------+
   | ``n_turb``                | integer          | Number of driving                  |       
   +---------------------------+------------------+------------------------------------+
   | ``drivingWN_DT``          | double precision | Driving wavenumber                 |
   +---------------------------+------------------+------------------------------------+
   | ``Energy_DT``             | double precision | Energy injected in each driving    |
   +---------------------------+------------------+------------------------------------+
   | ``netmomx_DT``            | double precision | Net momentum in x direction        |
   +---------------------------+------------------+------------------------------------+
   | ``netmomy_DT``            | double precision | Net momentum in y direction        |
   +---------------------------+------------------+------------------------------------+
   | ``netmomz_DT``            | double precision | Net momentum in z direction        |
   +---------------------------+------------------+------------------------------------+


   
.. _table:SelfGravity:
.. table:: SelfGravity:

   +---------------------------+------------------------------------+------------------------------------+
   | **Variables**             | **Values**                         | **notes**                          |
   +===========================+====================================+====================================+
   | ``SelfGravity``           | .true./.false.                     |                                    |
   +---------------------------+------------------------------------+------------------------------------+
   | ``sgBdryType``            | 0                                  | isolated                           |            
   +---------------------------+------------------------------------+------------------------------------+
   |                           | 1                                  | periodic                           |             
   +---------------------------+------------------------------------+------------------------------------+
   | ``sgfx`` ``sgfy`` ``sgfz``| :math:`g=-\nabla \phi`             | gravitational acceleration         |             
   +---------------------------+------------------------------------+------------------------------------+
   |                           | :math:`\rho v \cdot \boldsymbol{g}`| gravitational energy               |             
   +---------------------------+------------------------------------+------------------------------------+


.. note::
   If the periodic boundary condition for self-gravity ``sgBdryType = 1`` is used, the MPI periodic boundary condition ``periods(1)/(2)/(3)`` must be ``.true.``.

.. _table:Ambipolar Diffusion:
.. table:: Ambipolar Diffusion:

   +---------------------------+------------------------------------+------------------------------------+
   | **Variables**             | **Values**                         | **notes**                          |
   +===========================+====================================+====================================+
   | ``enable_ad``             | .true./.false.                     |                                    |
   +---------------------------+------------------------------------+------------------------------------+
   | ``mu_ad``                 | double precision                   | molecular mass in [amu]            |
   +---------------------------+------------------------------------+------------------------------------+
   | ``alpha_ad``              | double precision                   | coupling coefficient in [code unit]|
   +---------------------------+------------------------------------+------------------------------------+


