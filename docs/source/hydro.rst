.. _ch:hydro:

*************
Hydrodynamics
*************

Introduction
============

Governing equations and corresponding variables used in the code are presented here.

-------------------------
   
Conservation Forms 
==================

We begin with the compressible equations of pure hydrodynamics without viscosity for the conserved state vector,
:math:`\boldsymbol{U} = (\rho, \rho \boldsymbol{u}, \rho E):`

.. math::

   \begin{align}
   \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \boldsymbol{u})&= 0 , \\
   \frac{\partial (\rho \boldsymbol{u})}{\partial t} + \nabla \cdot (\rho \boldsymbol{u} \boldsymbol{u}) + \nabla P +\rho \nabla \Phi &=0, \\
   \frac{\partial (\rho E)}{\partial t} + \nabla \cdot (\rho \boldsymbol{u} E + P \boldsymbol{u}) &= - \rho \boldsymbol{u} \nabla \Phi. \label{eq:compressible-equations}
   \end{align}

Here :math:`\rho, \boldsymbol{u}, T, p`
velocity, temperature, pressure
respectively, and :math:`E = e + \boldsymbol{u} \cdot \boldsymbol{u} / 2` is the total
energy with :math:`e` representing the internal energy. 

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
velocity vector, specific heat ratio, and total energy density of neutrals. While :math:`\rho_i, \mu_i, P_i, \boldsymbol{v_i}, \Gamma_i`, and :math:`E_i=\frac{P_i}{\Gamma_i-1}+\frac{1}{2}\rho_i|\boldsymbol{v_i}|^2 + \frac{1}{2}|\boldsymbol{B}|^2` are the corresponding physical variables of ions, where :math:`|\boldsymbol{B}|` is the magnetic field vector.  
:math:`\alpha=\alpha_0 max\left(1,\frac{|\boldsymbol{v_i}-\boldsymbol{v_n}|}{\boldsymbol{v_{\alpha}}}\right)` is the collision coefficient, in which :math:`\alpha_0=\frac{1.9\times10^{-19}}{m_n+m_i} cm^3 s^{-1}` (where :math:`m_n=\mu_n m_H, m_i=\mu_i m_H` with :math:`m_H` being the mass of hydrogen atom in the unit of g)
, and :math:`\boldsymbol{v_{\alpha}} = 19.0 km s^{-`}` is the value of the drift velocity at which the Langevin approximation breaks down. Here we assume :math:`\mu_n=2.3` and :math:`\mu_i= 29` and the
corresponding collision coefficient is :math:`3.7\times10^{13} cm^3 s^{-1} g^{-1}`. 

The finite volume method solve for conserved variables (density, momenta, total energy, and cell-centered magnetic fields). Primitive variables can be obtained easily (density, velocities, pressure, and cell-centered magnetic fields).

   
   
Hydrodynamics Data Structures
=============================


.. _table:constants:
.. table:: constants:
    
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
   
.. _table:code units:   
.. table:: code units:

   +-----------------------+------------------------------+-----------+-------------------------------+
   | **code units**        | **quantity**                 | **value** |  **units**                    |
   +=======================+==============================+===========+===============================+
   | ``[T]``               | 9.78e5 yrs                   | 3.086e13s | pc*s/km                       |
   +-----------------------+------------------------------+-----------+-------------------------------+
   | ``[L]``               |                              |           | pc                            |
   +-----------------------+------------------------------+-----------+-------------------------------+
   | ``[M]``               |                              |           | Msun                          |
   +-----------------------+------------------------------+-----------+-------------------------------+
   | ``tff``               | sqrt(3*pi/(32*GravConst*rho))|           | free-fall time                |
   +-----------------------+------------------------------+-----------+-------------------------------+
   | ``B``                 | :math:`B_{phy}/\sqrt{4\pi}'  |           |                               |
   +-----------------------+------------------------------+-----------+-------------------------------+
   | ``B_{phy}``           | 8.0405e-7 [G,sqrt(g/cm)/s]   |           | 1.0 [sqrt(Msun)*km/(s*pc^1.5)]|
   +-----------------------+------------------------------+-----------+-------------------------------+

.. _table:variables:
.. table:: variables:
    
   +-----------------------+-----------------------+-----------+-------------------------------+
   | **variable()**        | **quantity**          | **int**   |  **units**                    |
   +=======================+=======================+===========+===============================+
   | ``density``           | :math:`\rho`          | 1         | [M_sun / pc^3]                |
   +-----------------------+-----------------------+-----------+-------------------------------+
   | ``momx``              | :math:`\rho v_x`      | 2         | [M_sun / pc^3 * km/s]         |
   +-----------------------+-----------------------+-----------+-------------------------------+
   | ``momy``              | :math:`\rho v_y`      | 3         | [M_sun / pc^3 * km/s]         |
   +-----------------------+-----------------------+-----------+-------------------------------+
   | ``momz``              | :math:`\rho v_z`      | 4         | [M_sun / pc^3 * km/s]         |
   +-----------------------+-----------------------+-----------+-------------------------------+
   | ``bx``                | :math:`B_x`           | 5         |                               |
   +-----------------------+-----------------------+-----------+-------------------------------+
   | ``by``                | :math:`B_y`           | 6         |                               |
   +-----------------------+-----------------------+-----------+-------------------------------+
   | ``bz``                | :math:`B_z`           | 7         |                               |
   +-----------------------+-----------------------+-----------+-------------------------------+
   | ``ene``               | :math:`E`             | 8(i) 5(n) |                               |
   +-----------------------+-----------------------+-----------+-------------------------------+

.. _table:coordType:
.. table:: coordType:
   
   +---------------------------+-----------+
   | **coordType**             | **int**   |
   +===========================+===========+
   | ``Cartesian``             | 1         |                               
   +---------------------------+-----------+
   | ``Cylindrical log r``     | 2         |                              
   +---------------------------+-----------+
   | ``Cylindrical uniform r`` | 3         |                              
   +---------------------------+-----------+

.. _table:boundaryType:
.. table:: boundaryType:
   
   +---------------------------+-----------+
   | **boundaryType**          | **int**   |
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
   | **solverType**            | **int**   |
   +===========================+===========+
   | ``exactHD``               | 1         |                               
   +---------------------------+-----------+
   | ``HLLHD``                 | 2         |                              
   +---------------------------+-----------+
   | ``HLLCHD``                | 3         |                              
   +---------------------------+-----------+
   | ``HLLMHD``                | 4         |                              
   +---------------------------+-----------+
   | ``HLLDMHD``               | 5         |                              
   +---------------------------+-----------+
   
.. _table:limiterType:
.. table:: limiterType:
   
   +---------------------------+-----------+
   | **limiterType**           | **int**   |
   +===========================+===========+
   | ``zero``                  | 0         |                               
   +---------------------------+-----------+
   | ``van Leer``              | 1         |                              
   +---------------------------+-----------+
   | ``fslop``                 | 2         |                              
   +---------------------------+-----------+  
   | ``minmod``                | 2         |                              
   +---------------------------+-----------+
   
.. _table:eosType:
.. table:: eosType:
   
   +---------------------------+-----------+
   | **eosType**               | **int**   |
   +===========================+===========+
   | ``isothermal``            | 1         |                               
   +---------------------------+-----------+
   | ``adiabatic``             | 2         |                              
   +---------------------------+-----------+
 

.. _table:setting:
.. table:: setting:
   
   +---------------------------+----------------+----------------------+
   | **setting**               | **int**        | **notes**            |
   +===========================+================+======================+
   | ``ndim``                  | 1/2/3          |                      |         
   +---------------------------+----------------+----------------------+
   | ``nbuf``                  | (2)            | buffer zone          |                   
   +---------------------------+----------------+----------------------+
   | ``nstep``                 | 0              | for MPI              |               
   +---------------------------+----------------+----------------------+
   | ``variable``              | 0              | default              |                
   +---------------------------+----------------+----------------------+
   | ``nMesh(1)/(2)/(3)``      | int            |                      |        
   +---------------------------+----------------+----------------------+
   | ``leftBdry(1)/(2)/(3)``   |                |                      |        
   +---------------------------+----------------+----------------------+
   | ``rightBdry(1)/(2)/(3)``  |                |                      |        
   +---------------------------+----------------+----------------------+
   | ``periods(1)/(2)/(3)``    | .true./.false. |                      |        
   +---------------------------+----------------+----------------------+
   | ``reorder``               | .true./.false. | ?????                |             
   +---------------------------+----------------+----------------------+
   | ``sndspd``                | int            | soundspeed   [km/s]  |                           
   +---------------------------+----------------+----------------------+
   | ``CFL``                   |                | CFL condition        |                     
   +---------------------------+----------------+----------------------+
   | ``gam``                   | 5.d0 / 3.d0    | useless if isothermal|                             
   +---------------------------+----------------+----------------------+
   | ``file_start``            |                | start time [Myrs]    |                         
   +---------------------------+----------------+----------------------+
   | ``time_end``              |                | end time   [Myrs]    |                         
   +---------------------------+----------------+----------------------+
   | ``dt_out``                |                | output time step     |                        
   +---------------------------+----------------+----------------------+
   | ``write_vtk``             | .true./.false. | vtk format output    |                        
   +---------------------------+----------------+----------------------+
   
.. _table:TurbulenceDriving:
.. table:: TurbulenceDriving:

   +---------------------------+---------------+----------------------------------+
   | **setting**               | **int**       | **notes**                        |
   +===========================+===============+==================================+
   | ``DriveTurbulence``       | .true./.false.|                                  |
   +---------------------------+---------------+----------------------------------+
   | ``DT_mode``               | 0             | Drive at begining                |            
   +---------------------------+---------------+----------------------------------+
   | ``DT_mode``               | 1             | Drive periodically               |             
   +---------------------------+---------------+----------------------------------+
   | ``dt_turb``               |               | turnover time                    |        
   +---------------------------+---------------+----------------------------------+
   | ``t_count_turb``          |               | ?????                            |
   +---------------------------+---------------+----------------------------------+
   | ``t_accum_turb``          |               | turnover time/??                 |         
   +---------------------------+---------------+----------------------------------+
   | ``E_turb_tot``            |               | Total injected turbulence energy |                          
   +---------------------------+---------------+----------------------------------+
   | ``E_turb``                |               | energy injected in each driving  |                          
   +---------------------------+---------------+----------------------------------+
   | ``zeta``                  | 0-1           | # of driving                     |       
   +---------------------------+---------------+----------------------------------+
   | ``DT_scale``              |               | k0 driving scale                 |           
   +---------------------------+---------------+----------------------------------+
   | ``n_turb``                |               | # of driving                     |       
   +---------------------------+---------------+----------------------------------+
   
.. _table:SelfGravity:
.. table:: SelfGravity:

   +---------------------------+---------------+----------------------------------+
   | **setting**               | **int**       | **notes**                        |
   +===========================+===============+==================================+
   | ``SelfGravity``           | .true./.false.|                                  |
   +---------------------------+---------------+----------------------------------+
   | ``sgBdryType``            | 0             | isolated[default]                |            
   +---------------------------+---------------+----------------------------------+
   |                           | 1             | periodic                         |             
   +---------------------------+---------------+----------------------------------+
   
problem setup
p = NkT  !!! thermal energy per volume, N number density
u = p/(gamma-1), thermal energy density, gamma adiabatic index Cp/Cv
u = rho*kT/((gamma-1)*\mu*mH)  (energy per volume)
K = Vu = M/rho*u = M*k*T/((gamma-1)*\mu*mH) !!! M total mass, V total volume
U = -3/5*GravConst*M^2/(Jeans length) potential Energy of a uniform density cloud
Virial Equilibrium => 2K = abs(U)
Jeans length: sqrt(5.0*c_s^2/(2.0*pi*GravConst*rho))
Jeans mass = 4.0*pi/3.0*rho*(Jean's length)^3
The corresponding volume density rho=68.4765 [Msun/pc^3]
Jeans length: Jlen = 0.5654 [pc]
Jeans mass: Jmass = 51.8439 [Msun]
critical mass-to-flux ratio: M/(\Phi_B)=sqrt(5/2)/(3*pi*sqrt(GravConst))=2.5581
The factor sqrt(4*pi) is absorbed in the definition of Bcode
alpha_code = pc^2km/sM = 3.669e13 cgs
   
Within the routines that implement the hydrodynamics, there are
several main data structures that hold the state.
Units and initial conditions used in :numref:`table:constants`.

-  interface variables: these are the time-centered interface states
   returned by the Riemann solver. They are used to discretize some
   non-conservative terms in the equations. These arrays are generally
   called ``q1``, ``q2``, and ``q3`` for the x, y, and z
   interfaces respectively. 


