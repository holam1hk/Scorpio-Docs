.. _ch:code_structure:

**************
Code Structure
**************
Getting Started
===============
Procedures ::

    make  //Compiling the code
    mpiexec -np 8 ./Scorpio > log &
    make clean
	
- Compiler
- Makefile, several flags too
- libraries

Different platforms
- cluster2
- gpu2/mike
- scorpio
- tianhe

conda init

conda activate

in (base) now

# To create an env named `my_new_env` by clone `base`
conda create -n my_new_env --clone base

# To create a new env named `my_new_env` from scratch
conda create -n my_new_env python=3

conda activate my_new_env

https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html

conda install <package_name>

Updating conda

Update conda update --all

Ramses (Teyssier 2002), 
PLUTO (Mignone et al. 2007), 
ENZO (Wang & Abel 2009), and 
FLASH (Fryxell et al. 2000).
Athena ++
Zeus
GIZMO 

SPH?


Input file
``TestSuite.f90``


.. warning:: Don't use too many cpu cores

and compute the fluxes of these quantities through the interfaces of the zones (this is a finite-volume approach). 
Parallelization is achieved by domain decomposition. We divide our domain into many smaller boxes, and distributed these across processors. 

Introduction
============
The file ``main.f90`` includes::

    use testSuiteMPI
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)
    call setMPI(np) !nprocs=np=number of processes this subroutine is in gridmodule.f30
    call setTestOnOff(.true.)
    call MPI_FINALIZE(ierr)
	
End

Stardard case setting ``testSuiteMPI.f90``::

	integer :: gridID
    type(grid) :: g1
    integer :: ndim, nbuf, coordType, variable(8)
    integer :: nMesh(3), dims(3)
    double precision :: leftBdry(3), rightBdry(3)
    logical :: periods(3), reorder
    double precision :: t0, t1
    integer :: nstep, ierr, N
    gridID !! case ID
	nstep = 0
    variable = 0


    dims = (/0, 0, 0/)    ????
    call MPI_DIMS_CREATE(nprocs, ndim, dims, ierr)  
    periods(1) = .true.  ????
    periods(2) = .true.
    periods(3) = .true.
    reorder = .true.
	call g1%setTopologyMPI(ndim, dims, periods, reorder)  ???
	call g1%setGridID(gridID = gridID)  ???
	call g1%setTime(fstart = 0, tend = 0.02d0, dtout = 0.01d0)  !! time interval for data output
    call g1%setMesh(nMesh, leftBdry, rightBdry, nbuf, coordType, gridID)   !!calling ``setCoordinates.f90``
    call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene !!calling ``sgPlan.f90``  !!!!!!!!!!! claim memory for variables !!!!!!!!!
    call g1%setMPIWindows()
    call g1%setEoS(eosType = 2) !! isothermal, 2 adiabatic
    call g1%setadiGamma(gam = 5.d0 / 3.d0) !! ratio of heat capacity
    call g1%setCFL(CFL = 0.4d0) !! courant number
    call g1%setSlopeLimiter(limiterType = 3)  !! 0=>zero,1=>van Leer, 2=>fslop, 3=>minmod
    call g1%setSolverType(solverType = 5) !! 1=>exactHD,2=HLLHD, 3=HLLC, 4=AdiHLLMHD, 5=AdiHLLDMHD ?????
    call g1%setBoundaryType(boundaryType = 3)  !! 1=>zero gradient, 2=>reflective, 3=>periodic
    call g1%initVariable()  !!  choose init1/2/3d.f90 by ndim
    call g1%exchangeBdryMPI(g1%q, g1%winq)  !! calling exchgBdryMPI.f03
    call g1%setBoundary(g1%q)  !!  calling ``setBdry3D.f03``
    call g1%writeGrid()  !! calling output3d in ``inout.f90``
	call g1%griddt()  !! calling ``dt3D.f03``
    call g1%evolveGridRK2() !! calling rk2_3D.f03
	Initial
	Boundary

    call g1%enableDrivingTurbulence(DT_mode=1) !!! initialize fftw mpi
	isRestart=0 !!!Default unless setRestart is called
	
End	
``setCoordinates.f90`` includes ::
	
	#remarks: fortran can take negative indices. always define q[1-nbuf:nMesh+nbuf]
	dx=(rightBdry(i)-leftBdry(i))/dble(nMesh(i))
    do j=1-nbuf, nMesh(i)+nbuf  !! divide the grids from left-nbuf to right+nbuf
    dx(i)=dx
    xl(i)=leftBdry(i)+dble(j-1)*dx  !! leftmost cell left interface are nbuf away from the left bounday
    xr(i)=leftBdry(i)+dble(j  )*dx  !! cell right interface is dx away from left interface
    xc(i)=0.5d0*(xl(i)+xr(i))  !! cell center = average of left and right interface
	
End
	
	
``sgPlan.f90``	includes ::

	study more about 'fftw3-mpi.f03'
	


``init3D.f90`` includes ::

	init3d
	init3d_for_FFTW
	
End

``exchgBdryMPI.f03`` includes ::

    subroutine initMPIWindows3D(this,q,q1,q2,databuf1,databuf2)
    call MPI_SIZEOF(q(1,1,1,1),sizedouble,ierr)  ??????????
    datasize=(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)*nvar*sizedouble
    call MPI_WIN_CREATE(q ,datasize,sizedouble,MPI_INFO_NULL,MPI_COMM_WORLD,this%winq,ierr)
	
End

``setBdry3D.f03`` includes ::

    call back in testsuite
	
End

``inout.f90`` includes ::

    read3d

End	

``dt3D.f03`` includes ::

    dt_temp=1.d10
	dt_pressure=1.d10
    EOS=1
	solverType = 1,2
	vtot  !! total v
	wavespd=vtot+snd
    dt_temp=dmin1(dt_temp,dmin1(dmin1( dx(1)(i), dx(2)(j)), dx(3)(k))/wavespd*CFL)
	solverType = 4,5
	vtot=dsqrt(vsq)
    bsq=(bxc**2+byc**2+bzc**2)/rho
    bmin=dmin1(dmin1(dabs(bxc),dabs(byc)),dabs(bzc))  !! min b
    cfast=dsqrt(0.5d0*(snd**2+bsq+dsqrt((snd**2+bsq)**2-4.d0*snd**2*bmin**2/rho)))  !! ????
    wavespd=vtot+cfast
    dt_temp=dmin1(dt_temp,dmin1(dmin1(dx(1),dx(2)),dx(3))/wavespd*CFL)
	
	EOS=2
	solverType = 2,3
	pressure=(gam-1.d0)*(ene-0.5d0*rho*(vx**2+vy**2+vz**2))
    wavespd=vtot+dsqrt(gam*pressure/rho)
    dt_temp=dmin1(dt_temp,dmin1(dmin1(dx(1),dx(2)),dx(3))/wavespd*CFL)
	
	solver=4,5
	pressure=(gam-1.d0)*(ene-0.5d0*rho*vsq-0.5d0*bsq)
    bmin=dmin1(dmin1(dabs(bxc),dabs(byc)),dabs(bzc))
    cfast=dsqrt((gam*pressure+bsq+dsqrt((gam*pressure+bsq)**2.d0-4.d0*gam*pressure*bmin**2.d0))/(2.d0*rho))  !! !! ????
    wavespd=vtot+cfast
    dt_temp=dmin1(dt_temp,dmin1(dmin1(dx(1),dx(2)),dx(3))/wavespd*CFL)
	
	SG !! avoid large self gravity
    sgftot=dsqrt(sgfx**2+sgfy**2+sgfz**2)
    dt_temp=dmin1(dt_temp,0.2d0*(-vtot/sgftot+dsqrt(vtot**2/sgftot**2+2.d0*dmin1(dmin1(dx(1),dx(2)),dx(3))/sgftot)))  !! ????

	call MPI_ALLREDUCE(dt_temp,global_dt,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD,ierr)
	
	if (dt > toutput-t) then  !! check if dt is larger than 
        global_dt = toutput-t
        toutput=toutput+dtout
        fnum=fnum+1
    elseif(dt > tend-t) then
        global_dt = tend-t
        fnum=fnum+1
    endif
    dt=global_dt
	
End

``rk2.f90`` includes ::

    subroutine rk2ADsg_3D(nthis,qn,qn1,qn2,ithis,qi,qi1,qi2)
    use gridModule
    use riemannSolverModule
    use mpi
	
	solverAdiMHD3D  !! includes ``riemannSolverModule.f90``
	calcSelfgravity !!!!!! apply gravity !!!!!!! has ``calcSG.f90``???
    evolveAD3D  !! ``evolveAmbipolarDiffusion.f90``

    call nthis%exchangeBdryMPI(nthis%q1,nthis%winq1)
    call nthis%setBoundary(nthis%q1)
    call ithis%exchangeBdryMPI(ithis%q1,ithis%winq1)
    call ithis%setBoundary(ithis%q1)
   
    again for rk2 step 2 
    call MPI_ALLREDUCE(nthis%changeSolver,global_changeSolvern,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(ithis%changeSolver,global_changeSolveri,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)

End
 
``calcSG.f90``& ``sgKernel.f90`` & ``initSGWindows3D.f90`` includes ::

    i dont know
 
End
 
``riemannSolverModule.f90`` includes ::  
   
    !!!!!!!!!!!ask for more about this
    !!!!!!!!!The Harten-Lax-van Leer-Contact (HLLC) Riemann solver neutral and Harten-Lax-van Leer-Discontinuities (HLLD) Riemann solver ion
    !!!!!!30. Toro, E. F., Spruce, M., & Speares, W. (1994). Restoration of the contact surface in the HLL-Riemann solver. Shock waves, 4(1), 25-34.
    31. Miyoshi, T., & Kusano, K. (2005). A multi-state HLL approximate Riemann solver for ideal magnetohydrodynamics. Journal of Computational Physics, 208(1), 315-344.
    how about !!!!!!!! Gardiner & Stone, JCP, 2005, 205, 509?
   
End   
      
``evolveAmbipolarDiffusion.f90`` includes ::  

    evolveAD3D  !!!! D. A. Tilly, D. S. Balsara, C. Meyer, 2012, New Astronomy, 17, 368 !!!!

End  
   
``limiterModule.f90`` includes ::

    Minimod limiter 3  !! Bryan, Greg L., et al. Enzo: An adaptive mesh refinement code for astrophysics. The Astrophysical Journal Supplement Series, 2014, 211.2: 19.?
    !! Skinner & Ostriker, 2010, ApJS, 188, 290 ??????????????????

End

To keep B-field divergence-free ( ), the constrained-transport algorithm33-35 is adopted

33. Balsara, D. S., & Spicer, D. S. (1999). A staggered mesh algorithm using high order Godunov fluxes to ensure solenoidal magnetic fields in magnetohydrodynamic simulations. Journal of Computational Physics, 149(2), 270-292.
34. Gardiner, T. A., & Stone, J. M. (2005). An unsplit Godunov method for ideal MHD via constrained transport. Journal of Computational Physics, 205(2), 509-539.
35. Gardiner, T. A., & Stone, J. M. (2008). An unsplit Godunov method for ideal MHD via constrained transport in three dimensions. Journal of Computational Physics, 227(8), 4123-4141.

The TR-BDF2 (trapezoidal rule and backward-difference formula of order two) scheme36 is employed to overcome the numerical stiffness of the ion-neutral collision source term.
36. Tilley, D. A., Balsara, D. S., & Meyer, C. (2012). A numerical scheme and benchmark tests for non-isothermal two-fluid ambipolar diffusion. New Astronomy, 17(3), 368-376.


The above procedures are coupled into the second-order Runge-Kutta (RK2) time-integration method40
40. Mignone, A., et al. "PLUTO: a numerical code for computational astrophysics." The Astrophysical Journal Supplement Series 170.1 (2007): 228.

alpha_ad = drag cofficient
mu_ad(ion) = 29.0
mu_ad(neutral) = 2.3	
	
