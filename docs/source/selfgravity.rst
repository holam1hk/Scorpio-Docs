.. _ch:selfgravity:

************
Self-gravity
************

Introduction
============

evolveGridRK2()
Call calcSelfgravity after riemann

subroutine calcSelfgravity in gridModule
call calcSG3D/calcSG3Dperiodic


``calcSG.f90``
``sgPlan.f90``
The gravitational potential :math:`\Phi` is defined by the Poisson equation:

DCMPLX(X [,Y]) returns a double complex number where X is converted to the real component. If Y is present it is converted to the imaginary component. If Y is not present then the imaginary component is set to 0.0. If X is complex then Y must not be present. 

The Kernel-based method37-39 is adopted to solve the Poisson equation
37. Yen, Chien-Chang, et al. "Self-gravitational force calculation of infinitesimally thin gaseous disks." Journal of Computational Physics 231.24 (2012): 8246-8261.
38. Wang, Hsiang-Hsu, Ronald E. Taam, and David CC Yen. "Self-Gravitational Force Calculation of Infinitesimally Thin Gaseous Disks on Nested Grids." The Astrophysical Journal Supplement Series 224.2 (2016): 16.
39. Wang, Hsiang-Hsu, et al. "Self-gravitational Force Calculation of High-order Accuracy for Infinitesimally Thin Gaseous Disks." The Astrophysical Journal Supplement Series 242.2 (2019): 17.
