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
:math:`\boldsymbol{U} = (\rho, \rho \boldsymbol{u}, \rho E, \boldsymbol{B}):`
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
   | **code units**        | **quantity**                 |  **cgs**  |  **cgs**                      |
   +=======================+==============================+===========+===============================+
   | ``[T]``               | pc*s/km = 9.78e5 yrs         |           | 3.086e13 ms                   |
   +-----------------------+------------------------------+-----------+-------------------------------+
   | ``[L]``               | pc                           |           | 3.086e18 cm                   |
   +-----------------------+------------------------------+-----------+-------------------------------+
   | ``[M]``               | Solar mass Msun              |           | 1.9e33 g                      |
   +-----------------------+------------------------------+-----------+-------------------------------+
   | ``rho``               | [M_sun / pc^3]               |           | 6.46496e-23 g/cm^-3           |
   +-----------------------+------------------------------+-----------+-------------------------------+
   | ``v``                 | km/s                         |           | 10^5        cm/s              |
   +-----------------------+------------------------------+-----------+-------------------------------+
   | ``tff``               | sqrt(3*pi/(32*GravConst*rho))|           | free-fall time                |
   +-----------------------+------------------------------+-----------+-------------------------------+
   | ``B``                 | :math:`B_{phy}/\sqrt{4\pi}'  |           |                               |
   +-----------------------+------------------------------+-----------+-------------------------------+
   | ``B_{phy}``           | 8.0405e-7 [G,sqrt(g/cm)/s]   |           | 1.0 [sqrt(Msun)*km/(s*pc^1.5)]|
   +-----------------------+------------------------------+-----------+-------------------------------+
   | ``alpha``             | 5.012314e8 cm^3/g/s          |           | 1.0 [pc^2km/sM]               |
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

-  ``ene`` is only required for adiabatic EOS. 
-  The indices of variable q(i,j,k,nvar) counts from 1-nbuf to nMesh+nbuf. The last index refers to different int in :numref:`table:variables`. 

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
   | ``minmod``                | 3         |                              
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

.. _table:Variables:
.. table:: Variables:

   +---------------------------+-------------------------------------+----------------------------------+
   | **setting**               | **int**                             | **notes**                        |
   +===========================+=====================================+==================================+
   | ``sgfx`` ``sgfy`` ``sgfz``| :math:`g=-\nabla \phi`              | gravitational acceleration       |
   +---------------------------+-------------------------------------+----------------------------------+
   |                           | :math:`\rho v \cdot \boldsymbol{g}` | gravitational energy             |                       
   +---------------------------+-------------------------------------+----------------------------------+

problem setup
p = NkT  !!! thermal energy per volume, N number density

u = p/(gamma-1), thermal energy density, gamma adiabatic index Cp/Cv

u = rho*kT/((gamma-1)*\mu*mH)  (energy per volume)

Units and initial conditions used in :numref:`table:constants`.
