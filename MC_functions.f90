! Suite of functions to be used by the LightMatterMC "simulationLoss.f90"
! subroutine for the simulation of a single atom inside an optical tweezer
! in the presence of a light field
! Daniel S. Grun, Innsbruck 2023

module MC_functions
  use physical_parameters
  implicit none
  
contains

  real(8) function trapPotential(r,P,w0,alpha,lambd)
    implicit none
    real(8) :: P,w0,alpha,lambd,zR,x,y,z,U0
    real(8), dimension(3) :: r

    x=r(1); y=r(2); z=r(3)
    U0 = P*alpha / (pi*e0*c*w0**2)
    zR = pi * w0**2 / lambd

    trapPotential = -U0/(1+(z/zR)**2)*exp(-2/(w0**2)*(x**2+y**2)/(1+(z/zR)**2))
  end function trapPotential

  
  function trapPotDerivs(r,U,w0,lambd) result(derivs)
    implicit none
    real(8) :: U,w0,lambd,ax,ay,az,x,y,z,z_term,zR
    real(8), dimension(3) :: r
    real(8), dimension(3) :: derivs

    x=r(1); y=r(2); z=r(3)
    zR = pi*w0**2/lambd
    z_term = sqrt(1+(z/zR)**2)

    ax = 4/m *x/(w0**2*z_term**2)*U
    ay = 4/m *y/(w0**2*z_term**2)*U
    az = 2/m *z/(zR**4*w0**2*z_term**4)*(zR**2*(w0**2-2*(x**2+y**2))+w0**2*z**2)*U
    
    derivs(1)=ax; derivs(2)=ay; derivs(3)=az

  end function trapPotDerivs


  real(8) function kineticEnergy(v)
    implicit none
    real(8), dimension(3) :: v
    kineticEnergy = m/2 * (v(1)**2+v(2)**2+v(3)**2)
  end function kineticEnergy


  real(8) function checkLost(rvVector, P, w0)
    implicit none
    real(8), dimension(6) :: rvVector
    real(8) :: P, w0, potEnergy, kinEnergy, totalEnergy
    !real(8), external :: trapPotential, kineticEnergy
    real(8), dimension(3) :: r,v

    r = rvVector(1:3); v = rvVector(4:6)

    potEnergy = trapPotential(r,P,w0,alpha_GS,lambd_trap)
    kinEnergy = kineticEnergy(v)
    totalEnergy = potEnergy + kinEnergy
    
    if (totalEnergy >= 0) then
       checkLost = 1.0
    else if (totalEnergy < 0) then
       checkLost = 0.0
    endif
  end function checkLost
  
    
  function trapEvol(rvVector, P, w0, alpha)
    implicit none
    real(8), dimension(6) :: rvVector, trapEvol
    real(8), dimension(3) :: r,v, accels
    !real(8), external :: trapPotential
    !real(8), dimension(3), external :: trapPotDerivs
    real(8) :: P, w0, ax, ay, az, U, alpha
    integer :: i

    r = rvVector(1:3)
    v = rvVector(4:6)

    U = trapPotential(r,P,w0, alpha,lambd_trap)
    accels = trapPotDerivs(r,U,w0,lambd_trap)

    ax=accels(1); ay=accels(2); az=accels(3)
    
    do i=1,3
       trapEvol(i) = v(i)
    enddo

    trapEvol(4)=ax; trapEvol(5)=ay; trapEvol(6)=az-g

  end function trapEvol


  function odeRK4_solver(y0, dt, P, w0, alpha)
    implicit none
    real(8) :: dt, P, w0, alpha
    !real(8), dimension(6), external :: trapEvol
    real(8), dimension(6) :: y0, dtArray, k1, k2, k3, k4
    real(8), dimension(6) :: odeRK4_solver

    dtArray = [dt,dt,dt,dt,dt,dt]
    
    k1 = dtArray * trapEvol(y0, P, w0, alpha)
    k2 = dtArray * trapEvol(y0+0.5*k1, P, w0, alpha)
    k3 = dtArray * trapEvol(y0+0.5*k2, P, w0, alpha)
    k4 = dtArray * trapEvol(y0+k3, P, w0, alpha)
    
    odeRK4_solver = y0 + (k1 + 2*k2 + 2*k3 + k4)/6.0

  end function odeRK4_solver



  real(8) function getAcStarkShift(rvVector, P, w0, alpha_E)
    implicit none
    real(8), dimension(3) :: r
    real(8), dimension(6) :: rvVector
    !real(8), external :: trapPotential
    real(8) :: P, w0, alpha_E, U_g, U_e

    r = rvVector(1:3)
    U_g = trapPotential(r,P,w0,alpha_GS,lambd_trap)
    U_e = trapPotential(r,P,w0,alpha_E,lambd_trap)
  
    getAcStarkShift = -(U_e-U_g)/hbar
  end function getAcStarkShift



  function recoilVel(absProj, lambd, case)
    implicit none
    real(8), dimension(3) :: absProj, recoilDir, spontDir
    real(8) :: lambd, k, absProjLength, spontDirLength, random
    real(8), dimension(3) :: random3D
    character(len=3) :: case
    integer :: i
    real(8), dimension(3) :: recoilVel
    
    absProjLength = 0.0

    call random_number(random3D)
    
    spontDir = 1.0-2*random3D
    spontDirLength = 0.0
    
    do i=1,3
       absProjLength = absProjLength + absProj(i)**2
       spontDirLength = spontDirLength + spontDir(i)**2
    enddo
    
    spontDirLength = sqrt(spontDirLength)
    absProjLength = sqrt(absProjLength)
    
    do i=1,3
       absProj(i) = absProj(i)/absProjLength
       spontDir(i) = spontDir(i)/spontDirLength
    enddo

    if (case=="abs") then
       recoilDir = absProj
       !print*, "absorption"
    else if (case=="emi") then
       recoilDir = spontDir
       !print*, "emission"
    endif

    k = 2*pi/lambd

    recoilVel = hbar*k/m * recoilDir

    !print*, recoilDir

  end function recoilVel
       
end module MC_functions

    
  
  
