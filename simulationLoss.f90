! Monte-Carlo simulation of survival of a single atom inside the Tweezer during
! imaging with near-resonant light.
! Includes photon absorption and spontaneous emission in the presence of
! a light field
! Daniel S. Grun, Innsbruck 2023

subroutine create_simulation_loss(firstInitialCond, P,T,w0,titf,s0,nBeams,lambd,Gammas,absProj,delta,alpha,solution)

  use physical_parameters
  use MC_functions
  use photon_functions

  implicit none
  
  integer :: i, index
  real(8), dimension(3) :: recoil
  integer, intent(in) :: nBeams
  real(8), intent(in) ::  P,T,w0
  real(8), dimension(3), intent(in) :: titf
  real(8), dimension(nBeams, 3), intent(in) :: absProj
  real(8), dimension(6), intent(in) :: firstInitialCond
  real(8), dimension(nBeams), intent(in) :: s0, lambd, Gammas, delta, alpha
  real(8), dimension(3), intent(out) :: solution
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
  
  phScatt = 0.0
  currentTime = 0.0

  initialCond = firstInitialCond
  
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
           sols = odeRK4_solver(initialCond, dt, P, w0, alpha(index))

           currentTime = currentTime + dt
           passedTime = passedTime + dt

           initialCond = sols
        enddo
        recoil = recoilVel(absProj(index,:), lambd(index), "emi")
        !print*, recoil
        sols(4:6) = sols(4:6) + recoil
        initialCond = sols
     endif

     if (lost == 1.0) then
        lostTime = currentTime
        currentTime = tf+1
        solution(1)=1-lost; solution(2)=lostTime; solution(3)=phScatt
     else
        if (currentTime >= tf) then
           solution(1)=1-lost; solution(2)=lostTime; solution(3)=phScatt
        endif
     endif

  end do

end subroutine create_simulation_loss
        
     
   



    



    

    
  
  
    
  
