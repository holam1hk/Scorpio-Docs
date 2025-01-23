.. _ch:selfgravity:

************
Self-gravity
************

Introduction
============

``rk2.f90`` includes ::  

   rk2ADsg_3D(nthis, qn, qn1, qn2, ithis, qi, qi1, qi2)
   call this%calcSelfgravity(totalden)  !! this is called after Riemann solvers
   subroutine calcSelfgravity(this, q) in gridModule

   Choosing 2D/3D and calling subroutine calcSG3D in ``calcSG.f90``
   if (this%sgBdryType .eq. 0) then
                call calcSG3D(this, q, this%sgDensityBuffer, this%sgDenCmplx, this%sgfxKernel, this%sgfyKernel, &
                              this%sgfzKernel, this%sgfxCmplx, this%sgfyCmplx, this%sgfzCmplx, this%sgfx, this%sgfy, this%sgfz)
   elseif (this%sgBdryType .eq. 1) then
                call calcSG3Dperiodic(this, q, this%sgDenCmplx, this%sgfxCmplx, this%sgfyCmplx, this%sgfzCmplx, &
                                      this%sgfx, this%sgfy, this%sgfz)
   endif

The gravity contributes to momentum and energy equations ::

   q(i, j, k, 2) = q(i, j, k, 2) + den*sgfx*dt
   q(i, j, k, 5) = q(i, j, k, 5) + (momx*sgfx + momy*sgfy + momz*sgfz)*dt


``calcSG.f90`` includes ::  
   subroutine calcSG3Dperiodic(this, q, sgDenCmplx, sgfxCmplx, sgfyCmplx, sgfzCmplx, sgfx, sgfy, sgfz)
   subroutine calcSG3D(this, q, sgDensityBuffer, sgDenCmplx, sgfxKernel, sgfyKernel, sgfzKernel, sgfxCmplx, &
                    sgfyCmplx, sgfzCmplx, sgfx, sgfy, sgfz)
   

subroutine calcSG3Dperiodic(this, q, sgDenCmplx, sgfxCmplx, sgfyCmplx, sgfzCmplx, sgfx, sgfy, sgfz)

``sgPlan.f90`` includes ::
``sgKernel.f90`` includes ::
``initSGWindows3D.f90`` includes ::

The gravitational potential :math:`\Phi` is defined by the Poisson equation:



DCMPLX(X [,Y]) returns a double complex number where X is converted to the real component. If Y is present it is converted to the imaginary component. If Y is not present then the imaginary component is set to 0.0. If X is complex then Y must not be present. 

