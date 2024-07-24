.. _ch:turbulence:

******************
Turbulence Driving
******************

Introduction
============
``DTPlan.f90``
``calcDT.f90``

dt_turb

   call g1%calcDrivingTurbulence_MD(g1%q) --drive turb
   DT_mode=0 called once. 
   DT_mode=1 periodic driving

``driving_spectrum_3D`` includes ::

   driving_spectrum_3D = (km**6) * exp(-8.0 * km / drivingWN) !! power spectrum P(k)=E(k)/k^2
   driving_spectrum_3D = km ** (-2.0-2.0) ! burgers spectrum P = k^-2 / k^2
   driving_spectrum_3D = km ** (-5.0/3.0-2.0) ! Kolmogorov spectrum P = k^(-5/3) / k^2
