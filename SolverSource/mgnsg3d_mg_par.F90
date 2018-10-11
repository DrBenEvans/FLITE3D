!**********************************************
!*mgns3d_mg.f90                               *
!*                                            *
!*Solver file at multigrid abstraction level  *
!*                                            *
!*           mgns3d v 2.3.1                   *
!*                                            *
!*                                            *
!*16.08.99-                                   *
!**********************************************
!*Description:                                *
!* This file contains the main driver sub-    *
!* routines for the multigrid solver, i.e.    *
!* it creates and sets up the grids and does  *
!* the multigrid solution procedure.          *
!**********************************************
!
! for moving boundary...
!
module MultigridSolver
 use GridSolver
 use BoundarySolver
 use RestrictionOperator
 use ProlongationOperator
 use InputVariables

 type(GridSolverData),pointer :: grp(:) ! contains the single grid solvers
 type(InputVariablesData) :: ivd        ! contains the input variables 

 real(KIND=4) :: userandsysTime(2)
 real :: startTime,stopTime,firstTime

 type(GridParallelizationData),pointer :: baseParallelData

 integer :: numberOfGrids
 real :: averageTimePerCycle
 integer :: time0(8),time1(8),time2(8),time3(8),time4(8)
 character*120 :: computationFileName,inputFileName,resultFileName
 character*120 :: startupFileName,startupFileName2,residualFileName 
 character*120 :: liftFileName
 character*200 :: fullResultFileName,fullResidualFileName
 character*120 :: problemName
 character*30 :: timeBuffer
 integer :: COMPINFILE,RESOUTFILE,INPINFILE,LIFTFILE,nameLength,allocateStatus
 integer :: STARTFILE,STARTFILE2,RESIDUALOUTFILE,RUNFILE,residualNameLength,resultNameLength
 integer :: timestepNameLength,liftNameLength
 integer :: inputNameLength,computationNameLength,startupNameLength,startupNameLength2
 integer :: directoryNameLength,fullResultNameLength,fullResidualNameLength
 integer :: problemNameLength,outputNumber

 real :: initialResidual,currentResidual,stopAtResidualReduction
 integer :: physicalTimestepNumber,numberOfPhysicalTimesteps,inputTimestepNumber
 real :: physicalTimestep

 character :: itExt*5,itForm*4,itExt2*5,itForm2*4
 integer :: itLen,itLen2

 character :: eFormat*4,fileExtension*8
 integer :: extLen,solverCommandNameLen
 integer :: IDbuffs(1024) 
 logical :: isMaster,startFromFile,startFromFile2

 character(len=120) :: FileName

contains

!----------------------------------------------------------------------
 subroutine executePar(mast)
 ! main procedure
 IMPLICIT NONE

 logical :: mast

  isMaster = mast

  call communicatePar()
  call multiGridCycles()
 end subroutine executePar
!----------------------------------------------------------------------
 subroutine executeSeq()
 ! main procedure
 IMPLICIT NONE
 
  isMaster = .true.
  
  print *,' in executeSeq'
  call communicateSeq()
  call multiGridCycles()
 end subroutine executeSeq
!----------------------------------------------------------------------
 subroutine communicatePar()
 IMPLICIT NONE

  character :: tmp*200
  integer :: shift


 if(isMaster) then  
  write(*,*) ""
  write(*,'(A)') "*************************************************"
  write(*,'(A)') "   WELCOME TO MGNS3DPAR - 3D MULTIGRID N-S SOLVER   "
  write(*,'(A)') "*************************************************"
  write(*,*) ""


  write(*,'(A)',advance="no") "Enter filename (no extension): "
  OPEN(1,FILE='run.inp') !MODIF
!!  read(*,"(a)") FileName
  read(1,"(a)") FileName !MODIF
  CLOSE(1) !MODIF

  inputFileName = trim(FileName)//".inp"
  inputNameLength = nameLen(inputFileName)
  write(*,*) 'input file:', inputFileName(1:inputNameLength)

  if(inputNameLength>0) then
   INPINFILE = 21
   print *,' Input file name is: ', inputFileName
   open(INPINFILE,&
       file=inputFileName(1:inputNameLength),form='formatted',status='old')
   write(*,*) 'File '//inputFileName(1:inputNameLength) // ' opened ...'
  end if
  call readInputVariables(ivd,INPINFILE)
  if(inputNameLength>0) then
   if(ivd%restartNumber>0) then
    close(INPINFILE)
   endif
  end if

! spawn

  solverCommandNameLen = nameLen(ivd%solverCommand)

  baseParallelData%numberOfProcesses = ivd%numberOfProcesses

  baseParallelData%isMaster = .true.

  allocate(baseParallelData%processorIDs(baseParallelData%numberOfProcesses),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in executePar"
  allocate(baseParallelData%receivedMessageFromProcess(baseParallelData%numberOfProcesses),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in executePar"

  if(baseParallelData%numberOfProcesses>1) then 
   call master(baseParallelData%numberOfProcesses,baseParallelData%currentTID,&
    baseParallelData%currentDomain,baseParallelData%processorIDs,& 
    ivd%solverCommand,solverCommandNameLen)
  else 
   baseParallelData%currentDomain = 1
  end if

  if (baseParallelData%currentDomain.le.9) then
    eFormat='(i1)'
    extLen = 1
  else if (baseParallelData%currentDomain.le.99) then
    eFormat='(i2)'
    extLen = 2
  else if (baseParallelData%currentDomain.le.999) then
    eFormat ='(i3)'
    extLen = 3
  endif

  write(fileExtension,eFormat) baseParallelData%currentDomain

  if(ivd%flowType.ne.0) then
    inputTimestepNumber=ivd%restartNumber
    if (inputTimestepNumber.le.9) then
     itForm='(i1)'
     itLen = 1
    else if (inputTimestepNumber.le.99) then
     itForm='(i2)'
     itLen = 2
    else if (inputTimestepNumber.le.999) then
     itForm ='(i3)'
     itLen = 3
    endif
    write(itExt,itForm) inputTimestepNumber 
                    
    inputTimestepNumber=mod(ivd%restartNumber,ivd%numberOfStepsPerCycle)
    if(inputTimestepNumber==0)inputTimestepNumber=ivd%numberOfStepsPerCycle
    if (inputTimestepNumber.le.9) then
     itForm2='(i1)'
     itLen2 = 1
    else if (inputTimestepNumber.le.99) then
     itForm2='(i2)'
     itLen2 = 2
      else if (inputTimestepNumber.le.999) then
     itForm2 ='(i3)'
     itLen2 = 3
    endif
    write(itExt2,itForm2) inputTimestepNumber                                  ! itExt2 : time step number in cycle 

    if(ivd%flowType.eq.2) then
      computationFileName = trim(FileName)//'_'//itExt2(1:itLen2)// ".sol"
      computationNameLength = nameLen(computationFileName)
      COMPINFILE = 20
      write(*,*) 'Trying to open '//computationFileName(1:computationNameLength)//'_'//fileExtension(1:extLen) // ' ...'
      open(COMPINFILE,&
         file=computationFileName(1:computationNameLength)//'_'//fileExtension(1:extLen),form='unformatted',status='old')
      write(*,*) 'File '//computationFileName(1:computationNameLength)//'_'//fileExtension(1:extLen) // ' opened ...'
    else
      computationFileName = trim(FileName)// ".sol"
      computationNameLength = nameLen(computationFileName)
      COMPINFILE = 20
      write(*,*) 'Trying to open '//computationFileName(1:computationNameLength)//'_'//fileExtension(1:extLen) // ' ...'
      open(COMPINFILE,&
         file=computationFileName(1:computationNameLength)//'_'//fileExtension(1:extLen),form='unformatted',status='old')
      write(*,*) 'File '//computationFileName(1:computationNameLength)//'_'//fileExtension(1:extLen) // ' opened ...'
    end if

    STARTFILE = 0
    if(ivd%restartNumber.eq.0)then                     
     startupFileName = ""
    else        
     startupFileName = trim(FileName) // ".rst"         
    endif            
    
    startupNameLength = nameLen(startupFileName)
    if(startupNameLength>0) then
     STARTFILE = 22
     open(STARTFILE,&
         file=startupFileName(1:startupNameLength)//'_'//fileExtension(1:extLen),form='unformatted',status='old')
     write(*,*) 'File '//startupFileName(1:startupNameLength)//'_'//fileExtension(1:extLen) // ' opened ...'
    end if 

    RESOUTFILE=26
    resultFileName = trim(FileName)//'_'//itExt2(1:itLen2)// ".res"
    resultNameLength = nameLen(resultFileName)
  
    RESIDUALOUTFILE=0
    residualFileName = trim(FileName)//'_'//itExt2(1:itLen2)// ".rsd"

    residualNameLength = nameLen(residualFileName)
    if(residualNameLength>0) then
     RESIDUALOUTFILE = 25
    end if

  else

    computationFileName = trim(FileName)// ".sol"
    computationNameLength = nameLen(computationFileName)
    COMPINFILE = 20
    write(*,*) 'Trying to open '//computationFileName(1:computationNameLength)//'_'//fileExtension(1:extLen) // ' ...'
    open(COMPINFILE,&
     file=computationFileName(1:computationNameLength)//'_'//fileExtension(1:extLen),form='unformatted',status='old')
    write(*,*) 'File '//computationFileName(1:computationNameLength)//'_'//fileExtension(1:extLen) // ' opened ...'

    STARTFILE = 0
    if(ivd%restartNumber.eq.0)then
     startupFileName = ""
    else
     startupFileName = trim(FileName) // ".rst"
    endif

    startupNameLength = nameLen(startupFileName)
    if(startupNameLength>0) then
     STARTFILE = 22
     open(STARTFILE,&
         file=startupFileName(1:startupNameLength)//'_'//fileExtension(1:extLen),form='unformatted',status='old')
     write(*,*) 'File '//startupFileName(1:startupNameLength)//'_'//fileExtension(1:extLen) // ' opened ...'
    end if

    RESOUTFILE=26
    resultFileName = trim(FileName) // ".res"
    resultNameLength = nameLen(resultFileName)

    RESIDUALOUTFILE=0
    residualFileName = trim(FileName)//".rsd"

    residualNameLength = nameLen(residualFileName)
    if(residualNameLength>0) then
      RESIDUALOUTFILE = 25
    end if

  end if

  directoryNameLength = nameLen(ivd%dataDirectory)


  call mp_init_buffer(baseParallelData%processorIDs,-1,baseParallelData%sendFlag)
  position = 1

  call mp_pak_intg(baseParallelData%processorIDs,-1,baseParallelData%integerType,inputNameLength,baseParallelData%sendFlag)
  call mp_pak_string(baseParallelData%processorIDs,-1,inputFileName,inputNameLength,baseParallelData%sendFlag)
  call mp_pak_intg(baseParallelData%processorIDs,-1,baseParallelData%integerType,computationNameLength,baseParallelData%sendFlag)
  call mp_pak_string(baseParallelData%processorIDs,-1,computationFileName,computationNameLength,baseParallelData%sendFlag)
  call mp_pak_intg(baseParallelData%processorIDs,-1,baseParallelData%integerType,startupNameLength,baseParallelData%sendFlag)
  call mp_pak_string(baseParallelData%processorIDs,-1,startupFileName,startupNameLength,baseParallelData%sendFlag)
  call mp_pak_intg(baseParallelData%processorIDs,-1,baseParallelData%integerType,resultNameLength,baseParallelData%sendFlag)
  call mp_pak_string(baseParallelData%processorIDs,-1,resultFileName,resultNameLength,baseParallelData%sendFlag)
  call mp_pak_intg(baseParallelData%processorIDs,-1,baseParallelData%integerType,residualNameLength,baseParallelData%sendFlag)
  call mp_pak_string(baseParallelData%processorIDs,-1,residualFileName,residualNameLength,baseParallelData%sendFlag)
  call mp_pak_intg(baseParallelData%processorIDs,-1,baseParallelData%integerType,directoryNameLength,baseParallelData%sendFlag)
  call mp_pak_string(baseParallelData%processorIDs,-1,ivd%dataDirectory,directoryNameLength,baseParallelData%sendFlag)

  call mp_snd_others(baseParallelData%numberOfProcesses,baseParallelData%currentDomain,&
          baseParallelData%processorIDs,2,baseParallelData%sendFlag)
  call mp_wait_comms( baseParallelData%sendFlag )

  call readComputationData(COMPINFILE,0)
  close(COMPINFILE)

  write(*,*) "YYIOP: ",physicalTimestep

  call setUpInitialField2(ivd,grp(1),STARTFILE)
  if(STARTFILE>0) then
   close(STARTFILE)
  end if

 else 
  !slaves read from same file

  write(*,*) "Hello - I'm a happy slave"


  call slave_get_tids (baseParallelData%numberOfProcesses,&
                       baseParallelData%currentDomain,IDbuffs,baseParallelData%currentTID)


  baseParallelData%isMaster = .false.

  allocate(baseParallelData%processorIDs(baseParallelData%numberOfProcesses),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in main"
  allocate(baseParallelData%receivedMessageFromProcess(baseParallelData%numberOfProcesses),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in main"

  baseParallelData%processorIDs(1:baseParallelData%numberOfProcesses) = IDbuffs(1:baseParallelData%numberOfProcesses)


!  receive info from master

  call init_recv_others(baseParallelData%numberOfProcesses,baseParallelData%currentDomain,&
                        baseParallelData%receivedMessageFromProcess)

  computationFileName = ' '
  startupFileName = ' '
  inputFileName = ' '
 

  call mp_recv( baseParallelData%processorIDs, 1, 2, baseParallelData%sendFlag )

  call mp_wait_comms( baseParallelData%sendFlag )

   call mp_upak_intg(baseParallelData%processorIDs,1,baseParallelData%integerType,inputNameLength,baseParallelData%sendFlag)
   call mp_upak_string(baseParallelData%processorIDs,1,inputFileName,inputNameLength,baseParallelData%sendFlag)
   call mp_upak_intg(baseParallelData%processorIDs,1,baseParallelData%integerType,computationNameLength,baseParallelData%sendFlag)
   call mp_upak_string(baseParallelData%processorIDs,1,computationFileName,computationNameLength,baseParallelData%sendFlag)
   call mp_upak_intg(baseParallelData%processorIDs,1,baseParallelData%integerType,startupNameLength,baseParallelData%sendFlag)
   call mp_upak_string(baseParallelData%processorIDs,1,startupFileName,startupNameLength,baseParallelData%sendFlag)
   call mp_upak_intg(baseParallelData%processorIDs,1,baseParallelData%integerType,resultNameLength,baseParallelData%sendFlag)
   call mp_upak_string(baseParallelData%processorIDs,1,resultFileName,resultNameLength,baseParallelData%sendFlag)
   call mp_upak_intg(baseParallelData%processorIDs,1,baseParallelData%integerType,residualNameLength,baseParallelData%sendFlag)
   call mp_upak_string(baseParallelData%processorIDs,1,residualFileName,residualNameLength,baseParallelData%sendFlag)
   call mp_upak_intg(baseParallelData%processorIDs,1,baseParallelData%integerType,directoryNameLength,baseParallelData%sendFlag)
   call mp_upak_string(baseParallelData%processorIDs,1,ivd%dataDirectory,directoryNameLength,baseParallelData%sendFlag)


  if(inputNameLength>0) then
   INPINFILE = 21
   tmp(1:directoryNameLength) = ivd%dataDirectory(1:directoryNameLength)
   shift = directoryNameLength+1
   tmp(shift:shift+inputNameLength) = inputFileName(1:inputNameLength)
   shift = shift+inputNameLength-1
   open(INPINFILE,&
       file=tmp(1:shift),form='formatted',status='old')
   write(*,*) 'File '//tmp(1:shift) // ' opened ...'
  end if
  call readInputVariables(ivd,INPINFILE)
  if(inputNameLength>0) then
   close(INPINFILE)
  end if



  if (baseParallelData%currentDomain.le.9) then
    eFormat='(i1)'
    extLen = 1
  else if (baseParallelData%currentDomain.le.99) then
    eFormat='(i2)'
    extLen = 2
  else if (baseParallelData%currentDomain.le.999) then
    eFormat ='(i3)'
    extLen = 3
  endif


  write(fileExtension,eFormat) baseParallelData%currentDomain


  tmp(1:directoryNameLength) = ivd%dataDirectory(1:directoryNameLength)
  shift = directoryNameLength+1
  tmp(shift:shift+computationNameLength-1) = computationFileName(1:computationNameLength)
  shift = shift+computationNameLength
  tmp(shift:shift) = '_'
  shift = shift+1
  tmp(shift:shift+extLen-1) = fileExtension(1:extLen)
  shift = shift+extLen-1
  COMPINFILE = 20
  open(COMPINFILE,&
     file=tmp(1:shift),&
     form='unformatted',status='old')
  write(*,*) 'File '//tmp(1:shift) // ' opened ...'

  call readComputationData(COMPINFILE,0)
  close(COMPINFILE)

  STARTFILE = 0
  write(*,*) "SA: ",startupNameLength,startupFileName,fileExtension(1:extLen)
  if(startupNameLength>0) then
   STARTFILE = 22
   tmp(1:directoryNameLength) = ivd%dataDirectory(1:directoryNameLength)
   shift = directoryNameLength+1
   tmp(shift:shift+startupNameLength-1) = startupFileName(1:startupNameLength)
   shift = shift+startupNameLength
   tmp(shift:shift) = '_'
   shift = shift+1
   tmp(shift:shift+extLen-1) = fileExtension(1:extLen)
   shift = shift+extLen-1

   write(*,*) 'opening '//tmp(1:shift) // ' ...'
   open(STARTFILE,&
       file=tmp(1:shift),&
       form='unformatted',status='old')
   write(*,*) 'File '//tmp(1:shift) // ' opened ...'
  end if

  call setUpInitialField2(ivd,grp(1),STARTFILE)
  if(startupNameLength>0) then
   close(STARTFILE)
  end if

  RESOUTFILE=26

  RESIDUALOUTFILE=0
  if(residualNameLength>0) then
   RESIDUALOUTFILE = 25
  end if
 end if 

  
  fullResultFileName(1:directoryNameLength) = ivd%dataDirectory(1:directoryNameLength)
  shift = directoryNameLength+1
  fullResultFileName(shift:shift+resultNameLength-1) = resultFileName(1:resultNameLength)
  shift = shift+resultNameLength
  fullResultFileName(shift:shift) = '_'
  shift = shift+1
  fullResultFileName(shift:shift+extLen-1) = fileExtension(1:extLen)
  shift = shift+extLen 
  fullResultNameLength = shift-1
  
  fullResidualFileName(1:directoryNameLength) = ivd%dataDirectory(1:directoryNameLength)
  shift = directoryNameLength+1
  fullResidualFileName(shift:shift+residualNameLength-1) = residualFileName(1:residualNameLength) 
  shift = shift+residualNameLength
  fullResidualNameLength = shift-1


 end subroutine communicatePar 
!----------------------------------------------------------------------
 subroutine multiGridCycles()
! handles entire multigrid process  
 IMPLICIT NONE
 
 integer :: i,itime,ii
 character :: tmp*200
 integer :: shift

 real(KIND=4), external :: etime

 real :: lift,drag,latF,frictionDrag,latM,longM,longMl,longMd,rollM
 real :: lift2 ,fpre,fshe, totalArea
!START-GUST
 integer :: ispan
 real :: localspan
 character :: eFormat1*4,fileExtension1*8
 integer :: extLen1
!END-GUST
 
 initialResidual = ivd%initialResidual
 if(ivd%flowType.eq.0) then

   ivd%restartNumber = 1
   ivd%numberOfphysicalTimesteps = 1

 end if

 do physicalTimestepNumber = ivd%restartNumber,ivd%numberOfphysicalTimesteps
!START-TURB
! initialResidual = 10000.0
!END-TURB

  if(isMaster)write(*,*)''
  write(*,*)">>>>>>>>>----> number of timestep = ",physicalTimestepNumber,grp(1)%pdp%currentDomain

  if(physicalTimestepNumber>ivd%restartNumber)then

    !update new file name and read new .sol ....

    inputTimestepNumber=physicalTimestepNumber
    if (inputTimestepNumber.le.9) then
     itForm='(i1)'
     itLen = 1
    else if (inputTimestepNumber.le.99) then
     itForm='(i2)'
     itLen = 2
    else if (inputTimestepNumber.le.999) then
     itForm ='(i3)'
     itLen = 3
    endif
    write(itExt,itForm) inputTimestepNumber

    inputTimestepNumber=mod(physicalTimestepNumber,ivd%numberOfStepsPerCycle)
    if(inputTimestepNumber==0)inputTimestepNumber=ivd%numberOfStepsPerCycle
    if (inputTimestepNumber.le.9) then
     itForm2='(i1)'
     itLen2 = 1
    else if (inputTimestepNumber.le.99) then
     itForm2='(i2)'
     itLen2 = 2
    else if (inputTimestepNumber.le.999) then
     itForm2 ='(i3)'
     itLen2 = 3
    endif
    write(itExt2,itForm2) inputTimestepNumber

    if(ivd%flowType.eq.2) then
      tmp(1:directoryNameLength) = ivd%dataDirectory(1:directoryNameLength)
      shift = directoryNameLength
      tmp(shift+1:shift+inputNameLength-4) = inputFileName(1:inputNameLength-4)
      shift = shift+inputNameLength-4
      tmp(shift+1:shift+1) = '_'
      shift = shift+1
      tmp(shift+1:shift+itLen2) = itExt2(1:itLen2)
      shift = shift+itLen2
      tmp(shift+1:shift+4) = '.sol'
      shift = shift+4
      tmp(shift+1:shift+1) = '_'
      shift = shift+1
      tmp(shift+1:shift+extLen) = fileExtension(1:extLen)
      shift = shift+extLen

      COMPINFILE = 20
      open(COMPINFILE,&
       file=tmp(1:shift),&
       form='unformatted',status='old')
      write(*,*) 'new data File '//tmp(1:shift) // ' opened ...'

      call readComputationData(COMPINFILE,1)
      close(COMPINFILE)
    end if

    fullResultFileName(1:directoryNameLength) = ivd%dataDirectory(1:directoryNameLength)
    shift = directoryNameLength
    fullResultFileName(shift+1:shift+inputNameLength-4) = inputFileName(1:inputNameLength-4)
    shift = shift+inputNameLength-4
    fullResultFileName(shift+1:shift+1) = '_'
    shift = shift+1
    fullResultFileName(shift+1:shift+itLen2) = itExt2(1:itLen2)
    shift = shift+itLen2
    fullResultFileName(shift+1:shift+4) = '.res'
    shift = shift+4
    fullResultFileName(shift+1:shift+1) = '_'
    shift = shift+1
    fullResultFileName(shift+1:shift+extLen) = fileExtension(1:extLen)
    shift = shift+extLen
    fullResultNameLength = shift

    fullResidualFileName(1:directoryNameLength) = ivd%dataDirectory(1:directoryNameLength)
    shift = directoryNameLength
    fullResidualFileName(shift+1:shift+inputNameLength-4) = inputFileName(1:inputNameLength-4)
    shift = shift+inputNameLength-4
    fullResidualFileName(shift+1:shift+1) = '_'
    shift = shift+1
    fullResidualFileName(shift+1:shift+itLen2) = itExt2(1:itLen2)
    shift = shift+itLen2
    fullResidualFileName(shift+1:shift+4) = '.rsd'
    shift = shift+4
    fullResidualNameLength = shift

  end if

  averageTimePerCycle = 0.0
  i = 0

 ! check what kind of iteration scheme to use

  if(ivd%ReynoldsNumber>0.0) then 
    do i=1,numberOfGrids
     grp(i)%laminarViscosity = 1./ivd%ReynoldsNumber
     if(i>1) then 
       call mapFineToCoarse(grp(i-1)%u,grp(i),grp(i)%u)
       call mapFineToCoarse(grp(i-1)%uold,grp(i),grp(i)%uold)
     end if
    end do
  else
    do i=1,numberOfGrids
     grp(i)%laminarViscosity = 0.0 
     if(i>1) then 
       call mapFineToCoarse(grp(i-1)%u,grp(i),grp(i)%u)
       call mapFineToCoarse(grp(i-1)%uold,grp(i),grp(i)%uold)
     end if
   end do
  end if
 
  currentResidual = 10000.0 
  stopAtResidualReduction = -1.0 
  if(ivd%numberOfMGIterations<0) then 
    stopAtResidualReduction = abs(dble(ivd%numberOfMGIterations))
  end if

! call calculateLiftAndDrag(grp(1),ivd,lift,drag)

  write(*,*) "Starting timestepping..."
  call mp_barrier()

  i = 1
  call date_and_time(timeBuffer(1:8),timeBuffer(9:18),timeBuffer(18:23),time0)
  call wtime(firstTime)

  if(stopAtResidualReduction<0) then 
    do i=1,ivd%numberOfMGIterations 
     call wtime(startTime)
     call oneCycle(i,physicalTimestepNumber)!START-GUST
     call wtime(stopTime)

    ! write residuals

     if(ivd%multigridScheme>0) then
!START-TURB
       call writeResiduals(ivd,grp(1),i)
     else
      call writeResiduals(ivd,grp(numberOfGrids),i)
!END-TURB   
     end if

   ! write to file if needed 

     if(mod(i,ivd%writeToFileInterval) == 0) then 
       call writeResultsToFile()
       if (ivd%explicit) then
	 if(isMaster) then
           write(*,*) "Current time: ", ivd%currentTime, "- Time step: ", ivd%physicalTimestep 
	 end if
!START-GUST
         if(ivd%nspan.ne.0) then
           Do ispan = 1,ivd%nspan
            localspan = (ivd%tipspan-ivd%rootspan)/ivd%nspan*ispan+ivd%rootspan
            call calculateLiftAndDrag2(grp(1),ivd,lift,drag,latF,frictionDrag,latM,longM,ispan,localspan)
           end do
         else
           if(ivd%forceCalcType==1)then
              call calculateLiftAndDrag(grp(1),ivd,lift,drag,latF,frictionDrag,latM,longM)
           elseif(ivd%forceCalcType==2)then
              call calculateLiftAndDrag(grp(1),ivd,lift,drag,latF,frictionDrag,latM,longM)
              call calculateLiftAndDrag2(grp(1),ivd,lift,drag,latF,frictionDrag,latM,longM,0,localspan)
           elseif(ivd%forceCalcType==3)then
              call calculateLiftAndDrag(grp(1),ivd,lift,drag,latF,frictionDrag,latM,longM)
              call calculateLiftAndDrag3(grp(1),ivd,lift,drag,latF,latM,longMl,longMd,rollM)
           endif
         endif

	 if(isMaster) then
           if(ivd%nspan.ne.0) then
             Do ispan=1,ivd%nspan
              if (ispan.le.9) then
               eFormat1 ='(i1)'
               extLen1 = 1
              else if (ispan.le.99) then
               eFormat1 ='(i2)'
               extLen1 = 2
              endif
              write(fileExtension1,eFormat1) ispan
 
              lift = grp(1)%localAerParam(ispan,1)
              drag = grp(1)%localAerParam(ispan,2)
              latF = grp(1)%localAerParam(ispan,3)
              frictionDrag = grp(1)%localAerParam(ispan,6)
              latM = grp(1)%localAerParam(ispan,7)
              longM = grp(1)%localAerParam(ispan,8)
              LIFTFILE = 28
	      if (i>1) then
                open(LIFTFILE,&
                file='liftdrag_'//fileExtension1(1:extLen1)//'.res',form='formatted',status='unknown',position='append')
              else
                open(LIFTFILE,&
                file='liftdrag_'//fileExtension1(1:extLen1)//'.res',form='formatted',status='unknown')
              end if
            
              write(LIFTFILE,'(7E14.5)') ivd%currentTime,lift,drag,latF,frictionDrag,latM,longM
              close(LIFTFILE)
             end do
           else
	     LIFTFILE = 28
	     if (i>1) then
	       open(LIFTFILE,&
	       file='liftdrag.res',form='formatted',status='unknown',position='append')
             else
	       open(LIFTFILE,&
	       file='liftdrag.res',form='formatted',status='unknown')
             end if
             if(ivd%forceCalcType==3)then
               write(LIFTFILE,'(9E14.5)') ivd%currentTime,lift,drag,latF,frictionDrag,latM,longMl,longMd,rollM
             else
	       write(LIFTFILE,'(7E14.5)') ivd%currentTime,lift,drag,latF,frictionDrag,latM,longM
             end if
	     close(LIFTFILE)
           end if
         end if
       end if
     end if 
     
     if (ivd%explicit.and.ivd%flowType.ne.0) then				  
!Updating current time
       ivd%currentTime = ivd%currentTime+ivd%physicalTimestep
!Stop at stopTime
       if (ivd%currentTime.ge.ivd%stopTime) then 
         call writeResultsToFile()
!START-GUST
         if(ivd%nspan.ne.0) then
           Do ispan = 1,ivd%nspan
            localspan = (ivd%tipspan-ivd%rootspan)/ivd%nspan*ispan+ivd%rootspan
            call calculateLiftAndDrag2(grp(1),ivd,lift,drag,latF,frictionDrag,latM,longM,ispan,localspan)
           end do
         else
           if(ivd%forceCalcType==1)then
              call calculateLiftAndDrag(grp(1),ivd,lift,drag,latF,frictionDrag,latM,longM)
           elseif(ivd%forceCalcType==2)then
              call calculateLiftAndDrag(grp(1),ivd,lift,drag,latF,frictionDrag,latM,longM)
              call calculateLiftAndDrag2(grp(1),ivd,lift,drag,latF,frictionDrag,latM,longM,0,localspan)
           elseif(ivd%forceCalcType==3)then
              call calculateLiftAndDrag(grp(1),ivd,lift,drag,latF,frictionDrag,latM,longM)
              call calculateLiftAndDrag3(grp(1),ivd,lift,drag,latF,latM,longMl,longMd,rollM)
           endif
         endif

         if(isMaster) then
           write(*,*) "Current time: ", ivd%currentTime, "- Time step: ", ivd%physicalTimestep 
         end if
         if(isMaster) then
           if(ivd%nspan.ne.0) then
             Do ispan=1,ivd%nspan
              if (ispan.le.9) then
               eFormat1 ='(i1)'
               extLen1 = 1
              else if (ispan.le.99) then
               eFormat1 ='(i2)'
               extLen1 = 2
              endif
              write(fileExtension1,eFormat1) ispan

              lift = grp(1)%localAerParam(ispan,1)
              drag = grp(1)%localAerParam(ispan,2)
              latF = grp(1)%localAerParam(ispan,3)
              frictionDrag = grp(1)%localAerParam(ispan,6)
              latM = grp(1)%localAerParam(ispan,7)
              longM = grp(1)%localAerParam(ispan,8)

              LIFTFILE = 28
              open(LIFTFILE,&
              file='liftdrag_'//fileExtension1(1:extLen1)//'.res',form='formatted',status='unknown',position='append')
              write(LIFTFILE,'(7E14.5)') ivd%currentTime,lift,drag,latF,frictionDrag,latM,longM
              close(LIFTFILE)
             end do
           else
             LIFTFILE = 28
             open(LIFTFILE,&
             file='liftdrag.res',form='formatted',status='unknown',position='append')
             if(ivd%forceCalcType==3)then
               write(LIFTFILE,'(9E14.5)') ivd%currentTime,lift,drag,latF,frictionDrag,latM,longMl,longMd,rollM
             else
               write(LIFTFILE,'(7E14.5)') ivd%currentTime,lift,drag,latF,frictionDrag,latM,longM
             end if
             close(LIFTFILE)
           end if
         end if
!Making .unk file
         if(isMaster) then
           write(*,*) ''
           write(*,*) 'Merging results to .unk file ...'
!START-TURB
           call makeplot(itExt,itLen,ivd)
!END-TURB
         end if 
         return
       end if
     end if
     averageTimePerCycle = averageTimePerCycle + stopTime - startTime
    end do
  else
    do while(log(initialResidual)-log(currentResidual)&
                <log(10.0)*stopAtResidualReduction)
     call wtime(startTime)
     call oneCycle(i,physicalTimestepNumber)!START-GUST
     call wtime(stopTime)

   ! write residuals

     if(ivd%multigridScheme>0) then
!START-TURB
       call writeResiduals(ivd,grp(1),i)
     else
       call writeResiduals(ivd,grp(numberOfGrids),i)
!END-TURB
     end if

   ! write to file if needed 

     if(mod(i,ivd%writeToFileInterval) == 0) then 
       call writeResultsToFile()
     end if 

     i = i+1
     if(i>ivd%maximumnumberofcycles) then
       print *,'EXITING JWJ'
       exit
     end if
    end do
  end if

 ! always write results at end of calculations

  call writeResultsToFile()
	
 ! write lift file

  if(ivd%nspan.ne.0) then
    Do ispan = 1,ivd%nspan
     localspan = (ivd%tipspan-ivd%rootspan)/ivd%nspan*ispan+ivd%rootspan
     call calculateLiftAndDrag2(grp(1),ivd,lift,drag,latF,frictionDrag,latM,longM,ispan,localspan)
    end do
  else
    if(ivd%forceCalcType==1)then
       call calculateLiftAndDrag(grp(1),ivd,lift,drag,latF,frictionDrag,latM,longM)
    elseif(ivd%forceCalcType==2)then
       call calculateLiftAndDrag(grp(1),ivd,lift,drag,latF,frictionDrag,latM,longM)
       call calculateLiftAndDrag2(grp(1),ivd,lift,drag,latF,frictionDrag,latM,longM,0,localspan)
    elseif(ivd%forceCalcType==3)then
       call calculateLiftAndDrag(grp(1),ivd,lift,drag,latF,frictionDrag,latM,longM)
       call calculateLiftAndDrag3(grp(1),ivd,lift,drag,latF,latM,longMl,longMd,rollM)
    endif
  end if
  if(isMaster) then
    if(ivd%nspan.ne.0) then
      Do ispan=1,ivd%nspan
       if (ispan.le.9) then
        eFormat1 ='(i1)'
        extLen1 = 1
       else if (ispan.le.99) then
        eFormat1 ='(i2)'
        extLen1 = 2
       endif
       write(fileExtension1,eFormat1) ispan

       lift = grp(1)%localAerParam(ispan,1)
       drag = grp(1)%localAerParam(ispan,2)
       latF = grp(1)%localAerParam(ispan,3)
       frictionDrag = grp(1)%localAerParam(ispan,6)
       latM = grp(1)%localAerParam(ispan,7)
       longM = grp(1)%localAerParam(ispan,8)
       LIFTFILE = 28
       if(physicalTimestepNumber>1) then
         open(LIFTFILE,&
         file='liftdrag_'//fileExtension1(1:extLen1)//'.res',form='formatted',status='unknown',position='append')
       else
         open(LIFTFILE,&
         file='liftdrag_'//fileExtension1(1:extLen1)//'.res',form='formatted',status='unknown')
       end if

       write(LIFTFILE,'(7E14.5)') physicalTimestepNumber*ivd%physicalTimeStep,lift,drag,latF,frictionDrag,latM,longM
       close(LIFTFILE)
      end do
    else
      LIFTFILE = 28
      if(physicalTimestepNumber>1) then
        open(LIFTFILE,&
        file='liftdrag.res',form='formatted',status='unknown',position='append')
      else
        open(LIFTFILE,&
        file='liftdrag.res',form='formatted',status='unknown')
      end if
      if(ivd%forceCalcType==3)then
        write(LIFTFILE,'(9E14.5)') physicalTimestepNumber*ivd%physicalTimeStep,lift,drag,latF,frictionDrag,latM,longMl,longMd,rollM
      else
        write(LIFTFILE,'(7E14.5)') physicalTimestepNumber*ivd%physicalTimeStep,lift,drag,latF,frictionDrag,latM,longM
      end if
      close(LIFTFILE)
    end if
  end if

  grp(1)%uold2=grp(1)%uold
  grp(1)%uold=grp(1)%u

  call mp_barrier()
  if(isMaster) then
    write(*,*) ''
    write(*,*) 'Merging results to .unk file ...'
!START-TURB
!   call makeplot(itExt,itLen,ivd)
!END-TURB
  end if
  
  if(isMaster) write(*,*) "Average time per cycle: ",averageTimePerCycle/dble(i-1)," seconds"

 end do

 end subroutine multiGridCycles
!-----------------------------------------------------------------------
 subroutine oneCycle(cycleNo,physicalTimestepNumber)!START-GUST
! one multigrid cycle
 IMPLICIT NONE

 integer :: cycleNo,ii,k

 integer :: i,j,numberOfSweeps,finestSweepLevel,localNumberOfGrids
 logical :: calculateTurbulence

 type(LinkedInteger),pointer :: currentInteger
 real,pointer :: diffield(:,:)
 real :: maxDiff
 integer :: maxDiffInd
 integer :: physicalTimestepNumber !START-GUST

 if(cycleNo.le.ivd%numberOfTurbulenceSteps) then 
  calculateTurbulence = .true.
 else
  calculateTurbulence = .false.
 end if


 if(ivd%multigridScheme.le.0) then 

 ! restriction part

  do i=1,numberOfGrids-1
   grp(i)%uprev = grp(i)%u
   call mapFineToCoarse(grp(i)%u,grp(i+1),grp(i+1)%u)
  end do

  grp(numberOfGrids)%sourceTerm = 0.0
  call doTimeIterations(grp(numberOfGrids),ivd,ivd%numberOfRelaxationSteps,.true.,calculateTurbulence,cycleNo,physicalTimestepNumber)!START-GUST
 
 
 ! prolongation part
 
  do i=numberOfGrids,2,-1
   call mapCoarseToFine(grp(i)%u,grp(i),grp(i-1),grp(i-1)%u,ivd)
  end do
 else if(ivd%multigridScheme==1) then 
 ! (V-cycle)

 ! restriction part
  do i=1,numberOfGrids-1
  ! call grid solver
   call doTimeIterations(grp(i),ivd,ivd%numberOfRelaxationSteps,i==1,calculateTurbulence,cycleNo,physicalTimestepNumber)!START-GUST
   call restrictGrid(grp(i),grp(i+1),ivd,i,calculateTurbulence)
   call mapFineToCoarse(grp(i)%u,grp(i+1),grp(i+1)%u)
  end do

  call doTimeIterations(grp(numberOfGrids),ivd,ivd%numberOfRelaxationSteps,numberOfGrids==1,calculateTurbulence,cycleNo,physicalTimestepNumber)!START-GUST
 
 ! prolongation part
 
  do i=numberOfGrids,2,-1
   call prolongateGrid(grp(i),grp(i-1),ivd)
   call doTimeIterations(grp(i-1),ivd,ivd%numberOfRelaxationSteps,(i-1)==1,calculateTurbulence,cycleNo,physicalTimestepNumber)!START-GUST
  end do

 else if(ivd%multigridScheme==2) then 
! W-cycle

  if(numberOfGrids<=2) then 
   numberOfSweeps = 1
  else
   numberOfSweeps = 2*numberOfGrids - 4
  end if

  finestSweepLevel = 1 
  do i=1,numberOfSweeps

   ! restriction part

   do j=finestSweepLevel,numberOfGrids-1
   ! call grid solver
    call doTimeIterations(grp(j),ivd,ivd%numberOfRelaxationSteps,j==finestSweepLevel,calculateTurbulence,cycleNo,physicalTimestepNumber)!START-GUST
    call restrictGrid(grp(j),grp(j+1),ivd,j,calculateTurbulence)
   end do
   call doTimeIterations(grp(numberOfGrids),ivd,ivd%numberOfRelaxationSteps,numberOfGrids==1,calculateTurbulence,cycleNo,physicalTimestepNumber)!START-GUST

   ! set start and end points for sweep

   if(finestSweepLevel==1) then 
    if(numberOfGrids>2) then
     finestSweepLevel=numberOfGrids-1
    else
     finestSweepLevel=numberOfGrids 
    end if
   else if(i<=numberOfSweeps/2) then
    finestSweepLevel = finestSweepLevel - 1 
   else
    finestSweepLevel = finestSweepLevel + 1 
   end if 

   ! prolongation part

   do j=numberOfGrids,finestSweepLevel+1,-1
    call prolongateGrid(grp(j),grp(j-1),ivd)
    call doTimeIterations(grp(j-1),ivd,ivd%numberOfRelaxationSteps,(j-1)==1,calculateTurbulence,cycleNo,physicalTimestepNumber)!START-GUST
   end do

  end do
   ! prolongate back to finest grid again 

   do j=numberOfGrids,2,-1
    call prolongateGrid(grp(j),grp(j-1),ivd)
    call doTimeIterations(grp(j-1),ivd,ivd%numberOfRelaxationSteps,(j-1)==1,calculateTurbulence,cycleNo,physicalTimestepNumber)!START-GUST
   end do

 else if(ivd%multigridScheme==3) then  
 ! growing V-cycle

  if(numberOfGrids==1) then 
   numberOfSweeps = 1
  else
   numberOfSweeps = numberOfGrids - 1 
  end if

  finestSweepLevel = 1 
  do i=1,numberOfSweeps

   ! restriction part

   do j=finestSweepLevel,numberOfGrids-1
   ! call grid solver
    call doTimeIterations(grp(j),ivd,ivd%numberOfRelaxationSteps,j==finestSweepLevel,calculateTurbulence,cycleNo,physicalTimestepNumber)!START-GUST
    call restrictGrid(grp(j),grp(j+1),ivd,j,calculateTurbulence)
   end do

   call doTimeIterations(grp(numberOfGrids),ivd,ivd%numberOfRelaxationSteps,numberOfGrids==1,calculateTurbulence,cycleNo,physicalTimestepNumber)!START-GUST

   ! set start and end points for sweep

   if(finestSweepLevel==1) then 
    if(numberOfGrids>=2) then
     finestSweepLevel=numberOfGrids-1
    else
     finestSweepLevel=numberOfGrids 
    end if
   else
    finestSweepLevel = finestSweepLevel - 1 
   end if 

   ! prolongation part

   do j=numberOfGrids,finestSweepLevel+1,-1
    call prolongateGrid(grp(j),grp(j-1),ivd)
 ! call doTimeIterations(grp(j-1),ivd,ivd%numberOfRelaxationSteps,(j-1)==1,calculateTurbulence,cycleNo,physicalTimestepNumber)!START-GUST
   end do

  end do
 
 else if(ivd%multigridScheme==4) then 
 ! decaying V-cycle

  if(numberOfGrids==1) then 
   numberOfSweeps = 1
  else
   numberOfSweeps = numberOfGrids - 1 
  end if

  finestSweepLevel = 1 
  do i=1,numberOfSweeps

   ! restriction part

   do j=finestSweepLevel,numberOfGrids-1
   ! call grid solver
    call doTimeIterations(grp(j),ivd,ivd%numberOfRelaxationSteps,j==finestSweepLevel,calculateTurbulence,cycleNo,physicalTimestepNumber)!START-GUST
    call restrictGrid(grp(j),grp(j+1),ivd,j,calculateTurbulence)
   end do
   call doTimeIterations(grp(numberOfGrids),ivd,ivd%numberOfRelaxationSteps,numberOfGrids==1,calculateTurbulence,cycleNo,physicalTimestepNumber)!START-GUST

   ! set start and end points for sweep

   finestSweepLevel = finestSweepLevel + 1 

   ! prolongation part

   do j=numberOfGrids,finestSweepLevel+1,-1
    call prolongateGrid(grp(j),grp(j-1),ivd)
    call doTimeIterations(grp(j-1),ivd,ivd%numberOfRelaxationSteps,(j-1)==1,calculateTurbulence,cycleNo,physicalTimestepNumber)!START-GUST
   end do
  end do

  ! prolongate back to finest grid again 

   do j=numberOfGrids,2,-1
    call prolongateGrid(grp(j),grp(j-1),ivd)
    call doTimeIterations(grp(j-1),ivd,ivd%numberOfRelaxationSteps,(j-1)==1,calculateTurbulence,cycleNo,physicalTimestepNumber)!START-GUST
   end do

 else if(ivd%multigridScheme==5) then 
! weighted W-cycle

  if(numberOfGrids<=2) then 
   numberOfSweeps = 1
  else
   numberOfSweeps = 2*numberOfGrids - 4
  end if

  finestSweepLevel = 1 
  do i=1,numberOfSweeps

   ! restriction part

   do j=finestSweepLevel,numberOfGrids-1
   ! call grid solver
    call doTimeIterations(grp(j),ivd,j,j==finestSweepLevel,calculateTurbulence,cycleNo,physicalTimestepNumber)!START-GUST
    call restrictGrid(grp(j),grp(j+1),ivd,j,calculateTurbulence)
   end do
   call doTimeIterations(grp(numberOfGrids),ivd,numberOfGrids,numberOfGrids==1,calculateTurbulence,cycleNo,physicalTimestepNumber)!START-GUST

   ! set start and end points for sweep

   if(finestSweepLevel==1) then 
    if(numberOfGrids>2) then
     finestSweepLevel=numberOfGrids-1
    else
     finestSweepLevel=numberOfGrids 
    end if
   else if(i<=numberOfSweeps/2) then
    finestSweepLevel = finestSweepLevel - 1 
   else
    finestSweepLevel = finestSweepLevel + 1 
   end if 

   ! prolongation part

   do j=numberOfGrids,finestSweepLevel+1,-1
    call prolongateGrid(grp(j),grp(j-1),ivd)
    call doTimeIterations(grp(j-1),ivd,ivd%numberOfRelaxationSteps,(j-1)==1,calculateTurbulence,cycleNo,physicalTimestepNumber)!START-GUST
   end do

  end do
   ! prolongate back to finest grid again 

   do j=numberOfGrids,2,-1
    call prolongateGrid(grp(j),grp(j-1),ivd)
    call doTimeIterations(grp(j-1),ivd,ivd%numberOfRelaxationSteps,(j-1)==1,calculateTurbulence,cycleNo,physicalTimestepNumber)!START-GUST
   end do
 else if(ivd%multigridScheme==7) then  
 ! growing V-cycle, every other grid
  localNumberOfGrids = numberOfGrids - 0.5*(numberOfGrids-1)
  if(numberOfGrids==1) then 
   numberOfSweeps = 1
  else
   numberOfSweeps = numberOfGrids - 1 
  end if

  finestSweepLevel = 1 
  do i=1,numberOfSweeps

   ! restriction part

   do j=finestSweepLevel,numberOfGrids-1
   ! call grid solver
    call doTimeIterations(grp(j),ivd,ivd%numberOfRelaxationSteps,j==finestSweepLevel,calculateTurbulence,cycleNo,physicalTimestepNumber)!START-GUST
    call restrictGrid(grp(j),grp(j+1),ivd,j,calculateTurbulence)
   end do
   call doTimeIterations(grp(numberOfGrids),ivd,ivd%numberOfRelaxationSteps,numberOfGrids==1,calculateTurbulence,cycleNo,physicalTimestepNumber)!START-GUST

   ! set start and end points for sweep

   if(finestSweepLevel==1) then 
    if(numberOfGrids>=2) then
     finestSweepLevel=numberOfGrids-1
    else
     finestSweepLevel=numberOfGrids 
    end if
   else
    finestSweepLevel = finestSweepLevel - 1 
   end if 

   ! prolongation part

   do j=numberOfGrids,finestSweepLevel+1,-1
    call prolongateGrid(grp(j),grp(j-1),ivd)
    call doTimeIterations(grp(j-1),ivd,ivd%numberOfRelaxationSteps,(j-1)==1,calculateTurbulence,cycleNo,physicalTimestepNumber)!START-GUST
   end do

  end do
 
 else if(ivd%multigridScheme==8) then 
 else
  write(*,*) "I: ",ivd%multigridScheme
  STOP "ERROR: Specified multigrid scheme not implemented"
 end if
 end subroutine oneCycle
!-----------------------------------------------------------------------
 subroutine readComputationData(INFILE,ia)
! reads data from file 
  IMPLICIT NONE

  integer :: INFILE,ia,ii

  integer :: i,j,k,allocateStatus,ibuff,numberOfGridDomains
  type(LinkedInteger),pointer :: currentInteger
  type(LinkedReal),pointer :: currentReal
  
 ! get number of grids to use
 
  read(INFILE) numberOfGrids
  
  if(ivd%numberOfGridsToUse>0) then
   numberOfGrids = min(numberOfGrids,ivd%numberOfGridsToUse)
  end if
  write(*,*) "Number of grids: ",numberOfGrids,ivd%numberOfGridsToUse
  
  read(INFILE) i,numberOfGridDomains
  

  if(ia.eq.0) then
    allocate(grp(numberOfGrids),stat=allocateStatus)
    if (allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
    grp(1)%pdp=baseParallelData 
  end if


  if(ia.eq.0) then
    allocate(grp(1)%pdp%processorIDs(grp(1)%pdp%numberOfProcesses),stat=allocateStatus)
    if (allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
    grp(1)%pdp%processorIDs=baseParallelData%processorIDs
    allocate(grp(1)%pdp%receivedMessageFromProcess(grp(1)%pdp%numberOfProcesses),stat=allocateStatus)
    if (allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"

    grp(1)%gridNumber = 1
    grp(1)%physicalTimeStep = ivd%physicalTimeStep

    if(numberOfGridDomains.ne.grp(1)%pdp%numberOfProcesses) then
     write(*,*) "ERROR: Number of processes specified differs from computation files",&
       grp(1)%pdp%numberOfProcesses,numberOfGridDomains
     stop
    end if
  end if

  write(*,*) "reading in grid number 1", grp(1)%gridNumber

  call readGridComputationData(grp(1),grp(1),grp(1),ivd,INFILE,ia)

! if(ia.eq.0) call makeIORegister(grp(1),ia)

  do i=2,numberOfGrids
   if(ia.eq.0) then
     grp(i)%pdp=baseParallelData
     allocate(grp(i)%pdp%processorIDs(grp(i)%pdp%numberOfProcesses),stat=allocateStatus)
     if (allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
     grp(i)%pdp%processorIDs=baseParallelData%processorIDs
     allocate(grp(i)%pdp%receivedMessageFromProcess(grp(i)%pdp%numberOfProcesses),stat=allocateStatus)
     if (allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"

     grp(i)%gridNumber = i
     grp(i)%physicalTimeStep = ivd%physicalTimeStep
   end if

   write(*,'(A,I2)') "reading in grid number ",i
   call readGridComputationData(grp(1),grp(i-1),grp(i),ivd,INFILE,ia)
!  if(ia.eq.0) call makeIORegister(grp(i),ia)
  end do

 end subroutine readComputationData
!----------------------------------------------------------------------
!START-TURB
 subroutine writeResiduals(ivd,grp,iterationNumber)
 IMPLICIT NONE

 type(InputVariablesData) :: ivd
 type(GridSolverData) :: grp
 integer :: iterationNumber

 integer :: i,j,ii,maxresInd(7)
 real :: res(7),maxres(7),cres
 real :: sec,lift,drag,ftmp,latF,frictionDrag,latM,longM,longMl,longMd,rollM
 real :: maxresCoor(7,3)
 integer :: maxresDomain(7),indtmp
 real :: maxtmp,coortmp(3)
 real :: currentRes
 integer :: nres

!START-GUST
integer :: ispan
real :: localspan
!END-GUST

if(ivd%turbulenceModel==0)nres=5 
if(ivd%turbulenceModel==1)nres=6 
if(ivd%turbulenceModel==2)nres=7
if(ivd%turbulenceModel==3)nres=7
!END-TURB

 res = 0.0
!START-TURB
 do i=1,nres
!END-TURB
   res(i) = sum((grp%u(:,i)-grp%uprev(:,i))&
            *(grp%u(:,i)-grp%uprev(:,i)))
 end do

#ifdef PARALLEL


 call mp_init_buffer(grp%pdp%processorIDs,-1,grp%pdp%sendFlag)
 position = 1
!START-TURB
 do i = 1,nres
!END-TURB
  call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,res(i),grp%pdp%sendFlag)
 end do

 call mp_snd_others(grp%pdp%numberOfProcesses,grp%pdp%currentDomain,grp%pdp%processorIDs,25,grp%pdp%sendFlag)

 call init_recv_others(grp%pdp%numberOfProcesses,grp%pdp%currentDomain,grp%pdp%receivedMessageFromProcess)

  call mp_recv_others(grp%pdp%numberOfProcesses,grp%pdp%currentDomain,grp%pdp%processorIDs,25,grp%pdp%sendFlag)

  call mp_wait_comms( grp%pdp%sendFlag )

  do i=1,grp%pdp%numberOfProcesses
    if( i.ne.grp%pdp%currentDomain ) then
!START-TURB
    do ii = 1,nres
!END-TURB
     call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
     res(ii)  = res(ii) + ftmp
    end do
  endif
 end do

#endif PARALLEL

!START-TURB
 do i=1,nres
!END-TURB
  res(i) = sqrt(res(i))
 end do

 currentResidual = res(1)
 if(initialResidual==10000.0) then 
  initialResidual = res(1) 
  write(*,*) ' Initial Residual is : ', initialResidual
 end if

 
 if(grp%gridNumber==1) then 
 !START-GUST
   if(ivd%nspan.ne.0) then
     ispan = ivd%nspan
     localspan = (ivd%tipspan-ivd%rootspan)/ivd%nspan*ispan+ivd%rootspan
!    call calculateLiftAndDrag2(grp,ivd,lift,drag,latF,frictionDrag,latM,longM,ispan,localspan)
     lift = grp%localAerParam(ispan,1)
     drag = grp%localAerParam(ispan,2)
     latF = grp%localAerParam(ispan,3)
     frictionDrag = grp%localAerParam(ispan,6)
     latM = grp%localAerParam(ispan,7)
     longM = grp%localAerParam(ispan,8)
   else
     if(ivd%forceCalcType==1)then
       call calculateLiftAndDrag(grp,ivd,lift,drag,latF,frictionDrag,latM,longM) 
     elseif(ivd%forceCalcType==2)then
       call calculateLiftAndDrag(grp,ivd,lift,drag,latF,frictionDrag,latM,longM) 
       call calculateLiftAndDrag2(grp,ivd,lift,drag,latF,frictionDrag,latM,longM,0,localspan)
     elseif(ivd%forceCalcType==3)then
       call calculateLiftAndDrag(grp,ivd,lift,drag,latF,frictionDrag,latM,longM) 
       call calculateLiftAndDrag3(grp,ivd,lift,drag,latF,latM,longMl,longMd,rollM)
     endif
   end if
 else
   lift = 0.0
   drag = 0.0
   latF = 0.0
   frictionDrag = 0.0
   latM = 0.0
   longM = 0.0
   longMl = 0.0
   longMd = 0.0
   rollM = 0.0
 end if

! find maximum residual 

 maxres = 0.
 maxresInd = 1
 do i=1,grp%numberOfNodes
!START-TURB
  do j=1,nres
!END-TURB
   cres = abs(grp%u(i,j)-grp%uprev(i,j))
   if(cres>maxres(j)) then 
    maxres(j) = cres
    maxresInd(j) = i
    maxResCoor(j,:) = grp%coordinates(maxresInd(j),:)
   end if
  end do
 end do

#ifdef PARALLEL

 call mp_init_buffer(grp%pdp%processorIDs,-1,grp%pdp%sendFlag)
 position = 1
!START-TURB
 do i = 1,nres
!END-TURB
  call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,maxres(i),grp%pdp%sendFlag)
  call mp_pak_intg(grp%pdp%processorIDs,-1,grp%pdp%integerType,maxresInd(i),grp%pdp%sendFlag)
  call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,grp%coordinates(maxresInd(i),1),grp%pdp%sendFlag)
  call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,grp%coordinates(maxresInd(i),2),grp%pdp%sendFlag)
  call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,grp%coordinates(maxresInd(i),3),grp%pdp%sendFlag)
 end do

 call mp_snd_others(grp%pdp%numberOfProcesses,grp%pdp%currentDomain,grp%pdp%processorIDs,107,grp%pdp%sendFlag)

 call init_recv_others(grp%pdp%numberOfProcesses,grp%pdp%currentDomain,grp%pdp%receivedMessageFromProcess)

  call mp_recv_others(grp%pdp%numberOfProcesses,grp%pdp%currentDomain,grp%pdp%processorIDs,107,grp%pdp%sendFlag)

  call mp_wait_comms( grp%pdp%sendFlag )

  do i=1,grp%pdp%numberOfProcesses
    if( i.ne.grp%pdp%currentDomain ) then
!START-TURB
    do ii = 1,nres
!END-TURB
     call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,maxtmp,grp%pdp%sendFlag)
     call mp_upak_intg(grp%pdp%processorIDs,i,grp%pdp%integerType,indtmp,grp%pdp%sendFlag)
     call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,coortmp(1),grp%pdp%sendFlag)
     call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,coortmp(2),grp%pdp%sendFlag)
     call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,coortmp(3),grp%pdp%sendFlag)
     if(maxtmp>maxRes(ii)) then
      maxres(ii) = maxtmp
      maxresInd(ii) = indtmp
      maxResCoor(ii,:) = coortmp
      maxResDomain(ii) = i
     end if
    end do
  endif
 end do

#endif PARALLEL


! write to screen
! VT: added a check if there is NAN
 if(isMaster) then 
 write(*,501) iterationNumber,res(1)*grp%residualScalingFactor," (",maxres(1),") ",&
                               res(2)*grp%residualScalingFactor," (",maxres(2),") ",&
                               res(3)*grp%residualScalingFactor," (",maxres(3),") ",&
                               res(4)*grp%residualScalingFactor," (",maxres(4),") ",&
                               res(5)*grp%residualScalingFactor," (",maxres(5),") ",&
                               lift,drag
 
! write(*,*) "mres: ",maxresInd,maxResDomain
! write(*,*) "mrescoor1: ",maxResCoor(1,:)
! write(*,*) "mrescoor2: ",maxResCoor(2,:)

!This is the check for NAN
!START-TURB
 do i=1,nres
 	if(ISNAN(res(i))) then
		write(*,*) "---> NAN in residual, exiting now ...."
		write(*,*) "---> stop in sub writeResiduals"
		STOP 
	end if
 end do
!END-TURB


 if(res(1)/initialResidual>1.0e-15) then 
  currentRes = log(res(1)/initialResidual)
 else
  currentRes = -100.0
 end if

! write to file
  if(RESIDUALOUTFILE>0) then 
   if(iterationNumber==1) then 
    open(RESIDUALOUTFILE,&
       file=fullResidualFileName(1:fullResidualNameLength),form='formatted',status='unknown')
   else
    open(RESIDUALOUTFILE,&
       file=fullResidualFileName(1:fullResidualNameLength),form='formatted',status='unknown',&
                                                            position='append')
   end if
   if(ivd%useTimeResidual) then 
     sec = stopTime - firstTime
!START-TURB
    if(ivd%forceCalcType==3)then
      write(RESIDUALOUTFILE,503) sec,currentRes/log(10.0),lift,drag,latF,frictionDrag,latM,longMl,longMd,rollM
    else
      write(RESIDUALOUTFILE,503) sec,currentRes/log(10.0),lift,drag,latF,frictionDrag,latM,longM,(log(res(i))/log(10.0),i=2,nres)
    endif
   else
    if(ivd%forceCalcType==3)then
      write(RESIDUALOUTFILE,504) iterationNumber,currentRes/log(10.0),lift,drag,latF,frictionDrag,latM,longMl,longMd,rollM
    else
      write(RESIDUALOUTFILE,504) iterationNumber,currentRes/log(10.0),lift,drag,latF,frictionDrag,latM,longM,(log(res(i))/log(10.0),i=2,nres)
    endif
   end if
!END-TURB
   close(RESIDUALOUTFILE)
  end if
 end if

501 format(I7,E14.5,A,E8.2,A,E14.5,A,E8.2,A,E14.5,A,E8.2,A,E14.5,A,E8.2,A,E14.5,A,E8.2,A,E12.4,E12.4)
502 format(I7,E14.5,E14.5,E14.5) 
!START-TURB
503 format(E14.5,E14.5,E14.5,12(f13.6,1x)) 
504 format(I7,E14.5,E14.5,E14.5,12(f13.6,1x))
!END-TURB
 end subroutine writeResiduals
!-----------------------------------------------------------------------
 subroutine writeResultsToFile()
 ! writes the solution fields of the finest grid to file
 IMPLICIT NONE

 integer :: i 
 real :: turbulenceCoefficient

 open(RESOUTFILE,&
     file=fullResultFileName(1:fullResultNameLength),form='unformatted',status='replace')
 write(*,*) 'File '//fullResultFileName(1:fullResultNameLength)// ' opened ...'

 call writeGridResults(grp(1),ivd,RESOUTFILE)
 close(RESOUTFILE)

 end subroutine writeResultsToFile
!----------------------------------------------------------------------
 subroutine communicateSeq()
 IMPLICIT NONE

  write(*,*) ""
  write(*,'(A)') "*************************************************"
  write(*,'(A)') "   WELCOME TO MGNS3DSEQ - 3D MULTIGRID N-S SOLVER   "
  write(*,'(A)') "*************************************************"
  write(*,*) ""

  write(*,'(A)',advance="no") "Enter control filename: "
  read(*,'(A)') inputFileName
  inputNameLength = nameLen(inputFileName)
  if(inputNameLength>0) then
   INPINFILE = 21
   open(INPINFILE,&
       file=inputFileName(1:inputNameLength),form='formatted',status='old')
   write(*,*) 'File '//inputFileName(1:inputNameLength) // ' opened ...'
  end if
  call readInputVariables(ivd,INPINFILE)
  if(inputNameLength>0) then
   close(INPINFILE)
  end if

  baseParallelData%isMaster = .true.

  write(*,'(A)',advance="no") "Enter computation filename: "
  read(*,'(A)') computationFileName
  computationNameLength = nameLen(computationFileName)
  COMPINFILE = 20
  open(COMPINFILE,&
     file=computationFileName(1:computationNameLength),form='unformatted',status='old')
  write(*,*) 'File '//computationFileName(1:computationNameLength)// ' opened ...'
  call readComputationData(COMPINFILE,0)
  close(COMPINFILE)

  STARTFILE = 0
  write(*,'(A)',advance="no") "Enter startup filename: "
  read(*,'(A)') startupFileName
  startupNameLength = nameLen(startupFileName)
  if(startupNameLength>0) then
   STARTFILE = 22
   open(STARTFILE,&
       file=startupFileName(1:startupNameLength),form='unformatted',status='old')
   write(*,*) 'File '//startupFileName(1:startupNameLength)// ' opened ...'
  end if
  call setUpInitialField(ivd,grp(1),STARTFILE)
  if(startupNameLength>0) then
   close(STARTFILE)
  end if

  RESOUTFILE=26
  write(*,'(A)',advance="no") "Enter result filename: "
  read(*,'(A)') resultFileName
  resultNameLength = nameLen(resultFileName)

  RESIDUALOUTFILE=0
  write(*,'(A)',advance="no") "Enter residual filename: "
  read(*,'(A)') residualFileName
  residualNameLength = nameLen(residualFileName)
  if(residualNameLength>0) then
   RESIDUALOUTFILE = 25
  end if
 end subroutine communicateSeq
!----------------------------------------------------------------------
!START-TURB
 subroutine makeplot(itExt,itLen,ivd)

 implicit none

 type(InputVariablesData) :: ivd

 character :: tmp*200,fm*4,char_ipr*4
 integer :: numberOfNodes,i,j,k,xlen,numberOfDomains
 integer :: numberOfLocalNodes,numberOfToNodes,switchLength
 real,pointer :: u(:,:),unk(:,:)
 integer,pointer :: transArray(:,:)
 integer :: problemNameLength,physicalTimestepNumber,numberOfPhysicalTimesteps

 character :: itExt*5
 integer :: itLen,numberOfSwitches,buffi,shift

 integer :: nvar

 if(ivd%turbulenceModel==0) nvar=5
 if(ivd%turbulenceModel==1) nvar=6
 if(ivd%turbulenceModel==2 .or. ivd%turbulenceModel==3) nvar=7 

! first read uold

  write(*,*) "opening ",'plotreg.reg'
  open(15,file='plotreg.reg',form='unformatted',status='old')

  read(15) numberOfNodes,numberOfDomains

  write(*,*) "number of nodes and domains: ",numberOfNodes,numberOfDomains

  allocate(u(numberOfNodes,nvar))

  write(*,*) "starting merging..."

  ! read translation array

 do i=1,numberOfDomains
  read(15) numberOfLocalNodes

  write(*,*) "number of local nodes for domain: ",numberOfLocalNodes

  allocate(transArray(numberOfLocalNodes,2))

  transArray = 0
 

  do j=1,numberOfLocalNodes
   read(15) transArray(j,:)
  end do

  if (i.le.9) then
   fm='(i1)'
   xlen = 1
  else if (i.le.99) then
   fm='(i2)'
   xlen = 2
  else if (i.le.999) then
   fm='(i3)'
   xlen = 3
  endif

  write(char_ipr,fm) i

  shift=fullResultNameLength-1
  tmp(1:shift)=fullResultFileName(1:shift)
  tmp(shift+1:shift+xlen)=char_ipr
  shift=shift+xlen

  write(*,*) "opening ",tmp(1:shift)
  open(16,file=tmp(1:shift),form='unformatted',status='old')

  read(16) numberOfToNodes

  allocate(unk(numberOfToNodes,nvar))
  read(16) ((unk(j,k),j=1,numberOfToNodes),k=1,nvar)
  close(16)

  do j=1,numberOfLocalNodes
   u(transArray(j,2),:) = unk(transArray(j,1),:)
  end do

  deallocate(unk)
 end do

 close(15)

 deallocate(transArray)


  shift = fullResidualNameLength-4
  tmp(1:shift)=fullResidualFileName(1:shift)

  tmp(shift+1:shift+4)='.unk'
  shift=shift+4

 write(*,*) "opening ",tmp(1:shift)
 open(17,file=tmp(1:shift),form='unformatted',status='unknown')

 write(17) numberOfNodes

 write(17) ((u(i,j),i=1,numberOfNodes),j=1,nvar)

 close(17)

 deallocate(u)

 write(*,*) "Merging done"
 end subroutine makeplot
!END-TURB
!----------------------------------------------------------------------
end module MultigridSolver






