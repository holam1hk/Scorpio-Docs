.. _ch:time_integration:

****************
Time Integration
****************



Introduction
============
``gridModule.f90`` calling ::   

    subroutine evolveGridRK2(this)
    subroutine evolveGridRK2sg(this)
    call rk2_1D(this,this%q,this%q1,this%q2)  from rk2.f03
    call rk2_2D(this,this%q,this%q1,this%q2)
    call rk2_3D(this,this%q,this%q1,this%q2)

End

``rk2.f03`` calling ::  

   call rieSolver(this,q,q,q1,dd=1)
   call rieSolver(this,q,q1,q1,dd=2)
   call rieSolver(this,q,q1,q1,dd=3)
   rieSolver is a pointer calling subroutine riemannSolver3D(this,q,q1,q2,dd) defined in ``riemannSolverModule.f03``
   then reference to solverAdiMHD3D e.g. rieSolver=>solverAdiMHD3D

End

``subroutine solverAdiMHD3D(this,q,q1,q2,dd)`` ::
     !!!!! dd=1 ==>x, dd=2 ==>y, dd=3 ==>z !!!!!
   procedure(limiter),pointer::slope=>null()    !! in ``riemannSolverModule.f03``
   case(0)
     slope=>zslop
   case(1)
     slope=>vslop
   case(2)
     slope=>fslop
   case(3)
     slope=>minmod
   procedure(fluxSolver),pointer::fluxPtr=>null()   !! in ``riemannSolverModule.f03``
   case (4)
     fluxPtr=>fluxHLLAdiMHD1D
   case (5)
     fluxPtr=>fluxHLLDAdiMHD1D

    Primitive variables like pressure, instead of energy, is used to compute the slope 
    SL(i,j,k,nvar)=slope(rhoL,rhoM,rhoR) !! slope
    ql(nvar) !!left interface state
    qr(nvar) !!right interface state

    Then calculate flux
    calling fluxSolver 
    which is subroutine fluxHLLDAdiMHD1D(ql,qr,slopeL,slopeR,gam,flux,nvar)

    UL=ql+0.5d0*slopeL   !!reconstructed left state
    UR=qr-0.5d0*slopeR   !!reconstructed right state

End
