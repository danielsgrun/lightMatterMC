module photon_functions
  use physical_parameters
  implicit none

contains

  real(8) function Rscatt(Gamma, s, Delta)
    implicit none
    real(8) :: Gamma, s, Delta
    Rscatt = Gamma/2. * (s/(1+s+4*(Delta/Gamma)**2))
  end function Rscatt


  real(8) function getDopplerShift(rvVector, absProj, lambd)
    implicit none
    real(8), dimension(3) :: absProj, v
    real(8), dimension(6) :: rvVector
    real(8) :: absProjLength, k, lambd, DopplerShift
    integer :: i
    
    v = rvVector(4:6)
    absProjLength = 0.0
    do i=1,3
       absProjLength = absProjLength + absProj(i)**2
    enddo

    absProjLength = sqrt(absProjLength)
    do i=1,3
       absProj(i) = absProj(i)/absProjLength
    enddo

    k = 2*pi/lambd
    DopplerShift = dot_product(absProj, v)

    getDopplerShift = -k*DopplerShift
  end function getDopplerShift

end module photon_functions
  
