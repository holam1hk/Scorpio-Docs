.. _ch:AD:

*************
Ambipolar Diffusion
*************

Introduction
============
Two-fluid source term
====================

.. math::

   \begin{align}
   \frac{\partial \rho_n}{\partial t} + \nabla \cdot (\rho_n \boldsymbol{u}_n)&= 0 , \\
   \frac{\partial (\rho_n \boldsymbol{u}_n)}{\partial t} + \nabla \cdot (\rho_n \boldsymbol{u}_n \boldsymbol{u}_n) + \nabla P_n +\rho_n \nabla \Phi 
   &= - \alpha \rho_n \rho_i (\boldsymbol{v}_n - \boldsymbol{v}_i), \\
   \frac{\partial E_n}{\partial t} + \nabla \cdot (\boldsymbol{u}_n E_n + P_n \boldsymbol{u}_n) 
   &= - \alpha \rho_n \rho_i (\boldsymbol{v}_n - \boldsymbol{v}_i) \cdot \boldsymbol{u}_n.
   \end{align}  
