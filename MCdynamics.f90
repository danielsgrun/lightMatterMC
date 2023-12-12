!  Monte-Carlo simulation of dynamics of a single atom inside the Tweezer during
! imaging with near-resonant light.
! Includes photon absorption and spontaneous emission in the presence of
! a light field
! Daniel S. Grun, Innsbruck 2023

subroutine create_simulation_loss(firstInitialCond, P,T,w0,titf,npts,s0,nBeams,lambd,Gammas,absProj,delta,alpha,solution)

  use physical_parameters
  use MC_functions
  use photon_functions

  implicit none
  
  integer :: i, index, count, npts
  real(8), dimension(3) :: recoil
  integer, intent(in) :: nBeams
  real(8), intent(in) ::  P,T,w0
  real(8), dimension(3), intent(in) :: titf
  real(8), dimension(nBeams, 3), intent(in) :: absProj
  real(8), dimension(6), intent(in) :: firstInitialCond
  real(8), dimension(nBeams), intent(in) :: s0, lambd, Gammas, delta, alpha
  real(8), dimension(3,npts), intent(out) :: solution
  real(8), dimension(6) :: sols, initialCond
  real(8) :: scattProb, DopplerShift, acStarkShift, deltaTotal, ti, tf, dt
  real(8), dimension(nBeams) :: auxRatios, auxMasks
  
  real(8) :: auxNum, auxRatio, waitTime, passedTime, phScatt, lost, lostTime, currentTime
  !real(8), external :: getDopplerShift, getAcStartShift, Rscatt
  !integer, external :: checkLost
  !real(8), dimension(6), external :: odeRK4_solver
  !real(8), dimension(3), external :: recoilVel

  call random_seed
  
  ti=titf(1); tf=titf(2); dt=titf(3)

  count = 1
  phScatt = 0.0
  currentTime = 0.0

  initialCond = firstInitialCond

  solution(:,count) = initialCond(1:3)
  
  do while (currentTime < tf) 
     sols = odeRK4_solver(initialCond, dt, P, w0, alpha_GS)
     auxMasks = 0.0
     auxRatios = 0.0

     do i=1,nBeams
        DopplerShift = getDopplerShift(sols, absProj(i,:), lambd(i))
        acStarkShift = getAcStarkShift(sols, P, w0, alpha(i))
        deltaTotal = 2*pi*delta(i) + DopplerShift + acStarkShift

        scattProb = Rscatt(Gammas(i), s0(i), deltaTotal)*dt
        
        call random_number(auxNum)

        auxRatios(i) = scattProb/auxNum

        !print*, auxNum, scattProb, scattProb/auxNum

     enddo

     auxRatio = maxval(auxRatios)
     index = findloc(auxRatios, auxRatio, 1)

     initialCond = sols
     
     currentTime = currentTime + dt

     lost = checkLost(initialCond, P, w0)

     count = count + 1

     solution(:,count) = sols(1:3)
     
     if (auxRatio >= 1) then
        recoil = recoilVel(absProj(index,:), lambd(index), "abs")
        sols(4:6) = sols(4:6) + recoil
        initialCond = sols
        phScatt = phScatt + 1

        call random_number(auxNum)

        waitTime = -log(1-auxNum)/Gammas(index)
        !print*, "Time in excited state:", waitTime*1e6, "us"
        
        passedTime = 0.0
        do while (passedTime < waitTime)
           count = count + 1
           sols = odeRK4_solver(initialCond, dt, P, w0, alpha(index))

           currentTime = currentTime + dt
           passedTime = passedTime + dt

           initialCond = sols
           solution(:,count) = sols(1:3)
        enddo
        recoil = recoilVel(absProj(index,:), lambd(index), "emi")
        !print*, recoil
        sols(4:6) = sols(4:6) + recoil
        initialCond = sols
     endif
     
  end do

end subroutine create_simulation_loss
