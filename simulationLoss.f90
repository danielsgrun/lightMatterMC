! Monte-Carlo simulation of N atoms inside an optical dipode trap along 
! with illumination with near-resonant light.
! Includes photon absorption and spontaneous emission in the presence of
! a light field along with light-induced 2-body interactions between atoms.

subroutine create_simulation_loss(firstInitialCond,P,T,w0,titf,C3vals,nAtoms,s0,nBeams,lambd,Gammas,absProj,&
     delta,alpha,spontCase,forcedSisyphusInput,modFreqInput,solution)

  use physical_parameters
  use MC_functions
  
  implicit none
  
  integer :: i, stopIndex, indexC3, j, scatteredPhotons, forcedSisyphus
  integer :: index_beam, indexAtom_excited, indexAtom_int, indexAtom_solo
  integer, parameter :: arrayLength = 1E4
  real(8), parameter :: arraySize = 1E4, C3null = 0
  real(8), dimension(arrayLength, 3) :: spontEmissionArray
  real(8), dimension(3) :: recoil, externalSpontDir, auxNum3D, lost
  real(8), dimension(2), intent(in) :: C3vals
  real(8), dimension(2) :: scattProbC3
  integer, intent(in) :: nBeams, nAtoms
  integer, intent(in), optional :: forcedSisyphusInput
  real(8), intent(in) ::  P,T,w0
  real(8), intent(in), optional :: modFreqInput
  real(8), dimension(3), intent(in) :: titf
  real(8), dimension(nBeams, 3), intent(in) :: absProj
  real(8), dimension(nAtoms,6), intent(in) :: firstInitialCond
  real(8), dimension(nBeams), intent(in) :: s0, lambd, Gammas, delta, alpha
  real(8), dimension(nAtoms+1), intent(out) :: solution
  real(8), dimension(nAtoms, 6) :: sols, initialCond ! 6-dimensional for matching the correct phase-space
  real(8) :: scattProb, DopplerShift, acStarkShift, deltaTotal, ti, tf, dt
  real(8), dimension(nBeams) :: auxRatios, auxMasks
  real(8), dimension(nAtoms) :: maxRatio_atoms, distances
  integer, dimension(nAtoms) :: index_atoms
  character(len=1), intent(in) :: spontCase
  real(8) :: auxNum, auxRatio, waitTime, passedTime, phScatt, minSep
  real(8) :: lostTime, currentTime, sum3D, maxProbC3, deltaC3, modFreq
  
  call random_seed  

  ! Reading from the dipole emission pattern file ! 
  spontEmissionArray = readSpontEmiFile(spontCase)


  ! Defining simulation time-step and initial and final simulation times !
  ti=titf(1); tf=titf(2); dt=titf(3)

  
  ! Initializing variables for scattered photons, simulation time and excited-state time !
  phScatt = 0.0
  currentTime = 0.0
  passedTime = 0.0


  ! initializing the initial-condition variable !
  initialCond = firstInitialCond


  ! Probably needs cleaning?? !
  scatteredPhotons = 0

  !lost = 0.


  ! Checking whether "forcedSisyphus" and "modFreq" were parsed !
  if (.not. present(forcedSisyphusInput)) then
     forcedSisyphus = 0 ! set to zero, so there's no modulation if forcedSisyphus is not parsed as 1
  else
     forcedSisyphus = forcedSisyphusInput
  end if

  if (.not. present(modFreqInput)) then
     forcedSisyphus = 0 ! set to zero, so there's no modulation if no frequency is provided
     modFreq = 1 ! doesn't matter, there will be no modulation with forcedSisyphus = 0
  else
     forcedSisyphus = forcedSisyphusInput
     modFreq = modFreqInput
  end if

  
  do while (currentTime < tf)

     do j=1,nAtoms

        lost(j) = checkLost(initialCond(j,:), P, w0, C3null)
        
        sols(j,:) = odeRK4_solver(initialCond(j,:),initialCond(j,:),dt,P,w0,alpha_GS,C3null)
        auxMasks = 0.0
        auxRatios = 0.0

     enddo

     !print*, sols

     ! print*, lost

     initialCond = sols

     do j=1,nAtoms

        do i=1,nBeams
        
           DopplerShift = getDopplerShift(sols(j,:), absProj(i,:), lambd(i))
           acStarkShift = getAcStarkShift(sols(j,:), P, w0, alpha(i))
           deltaTotal = 2*pi*delta(i) + DopplerShift + acStarkShift

           scattProb = Rscatt(Gammas(i), s0(i), deltaTotal)*dt
        
           call random_number(auxNum)

           auxRatios(i) = scattProb/auxNum

           !print*, auxNum, scattProb, scattProb/auxNum

        enddo

        maxRatio_atoms(j) = maxval(auxRatios) ! what is the biggest auxiliary ratio, and...
        index_atoms(j) = findloc(auxRatios, maxRatio_atoms(j), 1) ! ... for which beam it occurs, for each atom.

     enddo

     auxRatio = maxval(maxRatio_atoms) ! largest ratio within all atoms
     indexAtom_excited = findloc(maxRatio_atoms, auxRatio, 1) ! which atom it corresps. to
     index_beam = index_atoms(indexAtom_excited) ! which beam it corresponds to 

     initialCond = sols
     
     currentTime = currentTime + dt

     do j=1,nAtoms
        lost(j) = checkLost(initialCond(j,:), P, w0, C3null)
     enddo

     !lost = checkLost(initialCond, P, w0, C3null)
     
     if (auxRatio >= 1) then

        !print*, "Problem here?"

        scatteredPhotons = scatteredPhotons + 1 ! Decide on which of them
        
        recoil = recoilVel(absProj(index_beam,:), lambd(index_beam), "abs", absProj(index_beam,:)) ! parsing absProj(index,:)
        sols(indexAtom_excited,4:6) = sols(indexAtom_excited,4:6) + recoil          ! on last argument only for
                                                                                    ! only for type consistenty
        initialCond = sols                                                         
        phScatt = phScatt + 1 ! Decide on which of them 

        !        currentTime = currentTime + dt
        !        passedTime = dt
        
        call random_number(auxNum)

        waitTime = -log(1-auxNum)/Gammas(index_beam)
        !print*, "Time in excited state:", waitTime*1e6, "us"

        distances = 1

        do j=1,nAtoms
           if (j==indexAtom_excited) then
              distances(j) = 1E6 ! if comparing atom to itself, set very large distance artificially
           else
              distances(j) = getR12(sols(j,:), sols(indexAtom_excited,:))
           endif
        enddo

        minSep = minval(distances)
        indexAtom_int = findloc(distances,minSep,1)

        if (lost(indexAtom_int) == 1) then
           distances(indexAtom_int) = 1E5
           minSep = minval(distances)
           indexAtom_int = findloc(distances, minSep, 1)
        endif
        

        do j=1,nAtoms
           if ((j .ne. indexAtom_excited) .and. (j .ne. indexAtom_int)) then
              indexAtom_solo = j
           endif
        enddo
 
        do j=1,2
           DopplerShift = getDopplerShift(sols(indexAtom_excited,:), absProj(index_beam,:), lambd(index_beam))
           acStarkShift = getAcStarkShift(sols(indexAtom_excited,:), P, w0, alpha(index_beam))
           deltaC3 = getDetuningFromC3(C3vals(j), sols(indexAtom_excited,:), sols(indexAtom_int,:))
           deltaTotal = 2*pi*delta(index_beam) + DopplerShift + acStarkShift + deltaC3

           !print*, deltaC3

           call random_number(auxNum)
           
           scattProbC3(j) = Rscatt(Gammas(index_beam), s0(index_beam), deltaTotal)*dt / auxNum
        enddo
        
        maxProbC3 = maxval(scattProbC3)
        indexC3 = findloc(scattProbC3,maxProbC3,1)

        !print*, indexC3, scattProbC3

        !print*, "atom indices", indexAtom_solo, indexAtom_int, indexAtom_excited
        !print*, "beam index", index_beam

        passedTime = 0.0
        
        do while (passedTime < waitTime)
           sols(indexAtom_solo,:) = odeRK4_solver(initialCond(indexAtom_solo,:), &
                initialCond(indexAtom_solo,:), dt, P, w0, alpha(index_beam), C3null)
           
           sols(indexAtom_excited,:) = odeRK4_solver(initialCond(indexAtom_excited,:), &
                initialCond(indexAtom_int,:), dt, P, w0, alpha(index_beam), C3vals(indexC3))

           sols(indexAtom_int,:) = odeRK4_solver(initialCond(indexAtom_int,:), &
                initialCond(indexAtom_excited,:), dt, P, w0, alpha(index_beam), C3vals(indexC3))
           
           currentTime = currentTime + dt
           passedTime = passedTime + dt

           initialCond = sols

           !print*, indexAtom_excited, indexAtom_int, w0
        enddo

        call random_number(auxNum)

        externalSpontDir = spontEmissionArray(floor(arraySize*auxNum),:)
        
        recoil = recoilVel(absProj(index_beam,:), lambd(index_beam), "emi", externalSpontDir)

        sols(indexAtom_excited, 4:6) = sols(indexAtom_excited, 4:6) + recoil

        currentTime = currentTime + dt
        initialCond = sols
     endif

     
     if (sum(lost) == nAtoms) then
        lostTime = currentTime
        currentTime = tf+1
        do j=1,nAtoms
           solution(j) = 1-lost(j)
        enddo
        solution(nAtoms+1) = scatteredPhotons
     else
        if (currentTime >= tf) then
           do j=1,nAtoms
              solution(j) = 1-lost(j)
           enddo
           solution(nAtoms+1) = scatteredPhotons
        endif
     endif

  end do

end subroutine create_simulation_loss
        
     
   



    



    

    
  
  
    
  
