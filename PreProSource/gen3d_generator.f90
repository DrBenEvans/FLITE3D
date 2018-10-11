!**********************************************
!*Aggl3d_Agglomerator.f90                     *
!*                                            *
!*Main program module                         *
!*                                            * 
!*   Aggl3d v 2.3.1                           * 
!*                                            *
!**********************************************
!*Description:                                *
!* File controls the datastructure            *
!* creation and agglomeration for the pre-    *
!* processor.                                 *
!**********************************************

module Agglomerator3d
! main program module, handles everything

! Include the definitions of other modules. 
! These could be changed if other types are 
! wanted
 use BoundaryRegister ! handles boundary faces
 use TetrahedralElements ! handles elements
 use PrismElements
 use PyramidElements
 use HexahedralElements
 use CoordinateRegister ! provides coordinates for nodes in grid
 use Grid ! handles a grid   

 IMPLICIT NONE ! to disallow implicit variable declarations

 integer,parameter:: NSD = 3 ! number of space dimensions

 type(CoordinateRegisterData),pointer :: crp
 type(CoordinateRegisterData),pointer :: crpPrev
 type(CoordinateRegisterData),pointer :: crpPrev2
 type(TetrahedralElementData),pointer :: tep
 type(PrismElementData),pointer :: prp 
 type(PyramidElementData),pointer :: pyp 
 type(HexahedralElementData),pointer :: hep
 type(GridData),pointer :: grp(:)
 type(BoundaryRegisterData), pointer :: brp

 character :: problemName*80


 integer :: numberOfNodes,numberOfElements,numberOfSurfaceSegments,numberOfLineSegments
 integer :: numberOfBoundaryNodes,numberOfBoundarySides,numberOfBoundaryFaces
 integer :: numberOfTriangularBFs,numberOfRectangularBFs,numberOfSurfaceNormals
 integer :: numberOfBoundarySurfaceBlocks,numberOfBoundarySegmentBlocks                         
 integer :: numberOfHexEl,numberOfPriEl,numberOfPyrEl,numberOfTetEl
 integer :: numberOfGrids ! number of grids to produce (1: just fine grid)
 integer :: numberOfSides

 integer :: numberOfDomains ! for parallelization

 integer :: numberOfCycles,timestepNumber,stepInCycle,numberOfStepsPerCycle
 integer :: inputNumber,inputNumber2

 integer :: problemNameLength 
 integer :: InvVis
 real :: directionalityParameter,minimumAspectRatio,groundAngle
 real :: normalSmoothingFactor
 real :: wallDistanceLimit
 logical :: doVisualization,readPart,mergeLonesomeNodes,isHybrid,isRollGround
 integer :: visualizationMode,numberOfDirectionalAggl,normalSmoothingIterations
 contains

!--------------------------------------------------------------------------
  subroutine main()
  ! runs the show, calls other functions

  call communicate() ! communicates with user - calls procedures to get and
                     ! set up data from file
  call readFineGrid() ! reads fine grid from file and prepares data structure
  call setUpInitialGrid() ! prepares fine grid data structure for agglomeration
  call agglomerate() ! conducts agglomeration
  write(*,*) "Ha en god dag!" ! "Have a nice day!" in Norwegian
  end subroutine main
!--------------------------------------------------------------------------
  subroutine communicate()

  IMPLICIT NONE
 
  integer :: j,allocateStatus,RUNFILE

  write(*,*) "" 
  write(*,'(A)') "************************************************"
  write(*,'(A)') "    WELCOME TO AGGL3D - 3D GRID AGGLOMERATOR    "  
  write(*,'(A)') "************************************************"
  write(*,*) "" 

  ! read runfile
  problemNameLength = -1
  do while(problemNameLength.le.0)
   write(*,'(A)',advance="no") "Enter problem name: "
   read(*,'(A)') problemName
   problemNameLength = nameLen(problemName)
  end do
  write(*,'(A)',advance="no") " 1: Inviscid , 2: Viscous:   "
  read(*,*) InvVis
  write(*,'(A)',advance="no") "Number of grids to generate: "
  read(*,*) numberOfGrids 
  write(*,'(A)',advance="no") "Is mesh hybrid? (T/F): "
  read(*,*) isHybrid
  write(*,'(A)',advance="no") "Number of parallel domains: "
  read(*,*) numberOfDomains

  write(*,'(A)',advance="no") "Is there a rolling ground surface in the mesh? (T/F)"
  read(*,*) isRollGround

  if(isRollGround)then
    write(*,'(A)',advance="no") "Input the ground roll angle (in degrees):"
    read(*,*) groundAngle
  else
    groundAngle = 0.0
  endif

  write(*,'(A)',advance="no") " Input starting step and Number of step per cycles: "
  read(*,*) stepInCycle,numberOfStepsPerCycle
  
  inputNumber = stepInCycle

  write(*,*) '>>>>>>>----------->  Timestep No. : ',inputNumber


  directionalityParameter = 0.0

  normalSmoothingIterations = 0
  normalSmoothingFactor = 0.0

  doVisualization = .false.


  allocate(grp(numberOfGrids),stat=allocateStatus)
  if(allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary register"

  allocate(grp(1)%gid,stat=allocateStatus)
  if(allocateStatus /= 0) STOP "ERROR: Not enough memory in communicate"

  do j=2,numberOfGrids
   nullify(grp(j)%gid)
  end do


      write(*,*) 'checking the input data'    			!cw
      write(*,*) 'problemName',problemName		            !cw
      write(*,*) 'numberOfDomains',numberOfDomains                  !cw
      write(*,*) 'numberOfGrids',numberOfGrids		            !cw
      write(*,*) 'directionalityParameter',directionalityParameter  !cw
      write(*,*) 'isHybrid',isHybrid                                !cw

  end subroutine communicate  
!--------------------------------------------------------------------------
  subroutine agglomerate()

! Only need to store the current and previous
! grids, the rest can be written to file and deleted 
! from memory. The two first alone will then decide 
! the size limitations for the agglomerator

  integer :: i,visNameLength,allocateStatus
  integer :: VISOUTFILE,COMPOUTFILE  ! visualization and calculation output
                                     ! files 
  logical :: visualizeCurrentGrid
  character :: num*3

! set directionality parameter
  do i=1,numberOfGrids
   if(i-1<numberOfDirectionalAggl) then 
    grp(i)%directionalityParameter = directionalityParameter
   else
    grp(i)%directionalityParameter = 0.0 
   end if 
   grp(i)%minimumAspectRatio = minimumAspectRatio
  end do

  mergeLonesomeNodes = .true.

  write(*,*) "Writing grid data #",1," to file...",inputNumber
  grp(1)%numberOfGridDomains = numberOfDomains
  call writeParallelComputationFiles(.true.,numberOfGrids,1,problemName,grp(1),grp(1),crp,crpPrev,crpPrev2,inputNumber,numberOfStepsPerCycle)

  VISOUTFILE = -1 
  do i=1,numberOfGrids-1
   write(*,*) ""
   write(*,*) "Grid #",i+1,":" 
   if(doVisualization) then 
    write(*,'(A)',advance="no") "Visualize grid (T/F): "
    read(*,*) visualizeCurrentGrid
   else
    visualizeCurrentGrid = .false.
   end if
   if(visualizeCurrentGrid) then
    VISOUTFILE = 27
    write(num,*) i
    open(VISOUTFILE,file=problemName(1:problemNameLength)//num(1:i/10)//'.vis',form='formatted',status='unknown')
   end if
   write(*,*) "Agglomerating grid #",i
   allocate(grp(i+1)%brp,stat=allocateStatus)
   if(allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary register"
   call setUpNewBoundaryData(grp(i)%brp,grp(i+1)%brp)
   call doGridAgglomeration(grp(i),grp(i+1),crp,doVisualization,normalSmoothingFactor,normalSmoothingIterations)

   grp(i)%numberOfGridDomains = numberOfDomains
   grp(i+1)%numberOfGridDomains = numberOfDomains

   call makeDistanceFieldsForCoarseMesh(grp(i),grp(i+1))

   call writeParallelComputationFiles(.false.,numberOfGrids,i+1,problemName,grp(i),grp(i+1),crp,crpPrev,crpPrev2,inputNumber,numberOfStepsPerCycle)
  end do

  end subroutine agglomerate
!--------------------------------------------------------------------------
  subroutine setUpInitialGrid()
! This procedure sets up the initial grid for use. 
! The initial grid is given from the grid generator
! as a set of elements (possibly hybrid), 
! this is first transformed into a side-based structure.

  character :: itExt*5,itForm*4
  integer :: itLen,NEWPLTFILE
  
  call setBoundaryFirst(grp(1),crp,crpPrev,crpPrev2,hep,pyp,prp,tep,grp(1)%brp) !MOVE

  if(stepInCycle.eq.numberOfStepsPerCycle.and.numberOfStepsPerCycle.eq.1) then

   NEWPLTFILE = 51
   open(NEWPLTFILE,file='base.plt',form='unformatted',status='unknown')
   write(*,*) 'File base.plt'//' opened for rewrite...' 
  
  else

   write(*,*) "Writing renumbered plt file",stepInCycle,inputNumber

  if (inputNumber.le.9) then
    itForm='(i1)'
    itLen = 1
  else if (inputNumber.le.99) then
    itForm='(i2)'
    itLen = 2
  else if (inputNumber.le.999) then
    itForm ='(i3)'
    itLen = 3
  endif
  write(itExt,itForm) inputNumber

   NEWPLTFILE = 51
  !open(NEWPLTFILE,file=problemName(1:problemNameLength)//'_'//itExt(1:itLen)//'.plt',form='unformatted',status='old')
  !write(*,*) 'File '//problemName(1:problemNameLength)//'_'//itExt(1:itLen)//'.plt'//' opened for rewrite...' 
   open(NEWPLTFILE,file='base_'//itExt(1:itLen)//'.plt',form='unformatted',status='unknown')
   write(*,*) 'File base_'//itExt(1:itLen)//'.plt'//' opened for rewrite...' 
  
  end if

   call makeNewPltFile(grp(1),crp,hep,pyp,prp,tep,grp(1)%brp,NEWPLTFILE)

   close(NEWPLTFILE)


!  call makeBoundaryNormals(grp(1)%brp,crp)
!  call createBoundaryNormalRegister(grp(1)%brp)

  if(InvVis.eq.1) call findTrailingEdges(grp(1)%brp,crp)
  call setBoundIndForViscousMesh(grp(1)%brp)
  call makeNodeNormals(grp(1)%brp,crp)

  call makeCalculationData(crp,crpPrev,crpPrev2,hep,prp,pyp,tep,grp(1),grp(1)%brp,doVisualization) !MOVE

  call findBoundaryNodeArea(grp(1)%brp)
  call findNodeDistanceFromWall(grp(1),crp)
  end subroutine setUpInitialGrid 
!-------------------------------------------------------------------------- 
  subroutine readFineGrid()
! communicates with user and reads fine grid

  IMPLICIT NONE

  integer, parameter :: SURFINFILE = 25,VOLINFILE = 26,SEGINFILE = 27,WALLFILE = 28 
  integer :: i,j,ip,allocateStatus,dummy,numberOfBoundaryFacesSurf
  real :: dummyr
  character :: text*10

  integer :: DCOORFILE,DCOORFILE2,itLen,itLen2
  logical :: coorMov,readWALLFile
  character :: itExt*5,itForm*4,itExt2*5,itForm2*4
  

  open(SEGINFILE,file=problemName(1:problemNameLength)//'.bco',form='formatted',status='old')
  write(*,*) 'File '//problemName(1:problemNameLength)//'.bco'//' opened...' 
  write(*,*) "" 

  inquire(file=problemName(1:problemNameLength)//'.wid',exist=readWALLFile)
  if(readWALLFile) then
   open(WALLFILE,file=problemName(1:problemNameLength)//'.wid',form='formatted',status='old')
   write(*,*) 'File '//problemName(1:problemNameLength)//'.wid'//' opened...'
   write(*,*) "reading wall input data..."  
   call readGridInput(grp(1)%gid,WALLFILE)
   close(WALLFILE)
  else
   grp(1)%gid%numberOfTripLines = 0
  end if

  if(stepInCycle.eq.numberOfStepsPerCycle.and.numberOfStepsPerCycle.eq.1) then

   open(VOLINFILE,file=problemName(1:problemNameLength)//'.plt',form='unformatted',status='old')
   write(*,*) 'File '//problemName(1:problemNameLength)//'.plt'//' opened...' 

   DCOORFILE = -1
   DCOORFILE2= -1
  else

   if (inputNumber.le.9) then
     itForm='(i1)'
     itLen = 1
   else if (inputNumber.le.99) then
     itForm='(i2)'
     itLen = 2
   else if (inputNumber.le.999) then
     itForm ='(i3)'
     itLen = 3
   endif
   write(itExt,itForm) inputNumber

   write(*,*) "input number: ",inputNumber,itExt
   open(VOLINFILE,file=problemName(1:problemNameLength)//'_'//itExt(1:itLen)//'.plt',form='unformatted',status='old')
   write(*,*) 'File '//problemName(1:problemNameLength)//'_'//itExt(1:itLen)//'.plt'//' opened...' 

!MOVE

   inquire(file=problemName(1:problemNameLength)//'_'//itExt(1:itLen)//'.dcr',exist=coorMov)
   write(*,*) "OPENING: ",problemName(1:problemNameLength)//'_'//itExt(1:itLen)//'.dcr',coorMov
   if(coorMov) then
    DCOORFILE = 39
    write(*,*) "found dcr1 file"
    open(DCOORFILE,file=problemName(1:problemNameLength)//'_'//itExt(1:itLen)//'.dcr',form='unformatted',status='old')
   else
    DCOORFILE = -1
   end if


   if(inputNumber==1)then
     inputNumber2=numberOfStepsPerCycle
   else
     inputNumber2=inputNumber-1
   end if 

   if(inputNumber2.le.9) then
     itForm2 = '(i1)'
     itLen2 = 1
   else if (inputNumber2.le.99) then
     itForm2 = '(i2)'
     itLen2 = 2
   else if (inputNumber2.le.999) then
     itForm2 ='(i3)'
     itLen2 = 3
   endif
 
   write(itExt2,itForm2) inputNumber2
 
 
   inquire(file=problemName(1:problemNameLength)//'_'//itExt2(1:itLen2)//'.dcr',exist=coorMov)
   write(*,*) "OPENING: ",problemName(1:problemNameLength)//'_'//itExt2(1:itLen2)//'.dcr',coorMov
   if(coorMov) then
    DCOORFILE2 = 49
    write(*,*) "found dcr2 file"
    open(DCOORFILE2,file=problemName(1:problemNameLength)//'_'//itExt2(1:itLen2)//'.dcr',form='unformatted',status='old')
   else
    DCOORFILE2 = -1
   end if
 
!MOVE

  end if

! Start reading files

 if(isHybrid) then 
  read(VOLINFILE) numberOfElements, numberOfNodes, numberOfBoundaryFaces,&
                  numberOfHexEl,numberOfPriEl,numberOfPyrEl,numberOfTetEl,&
                  numberOfRectangularBFs,numberOfTriangularBFs,numberOfSides 
 else
  read(VOLINFILE) numberOfElements, numberOfNodes, numberOfBoundaryFaces
  numberOfTetEl = numberOfElements
  numberOfTriangularBFs = numberOfBoundaryFaces
  numberOfRectangularBFs = 0
  numberOfHexEl = 0
  numberOfPriEl = 0
  numberOfPyrEl = 0
 end if



  read(SEGINFILE,*) text
  read(SEGINFILE,*) numberOfSurfaceSegments,numberOfLineSegments 

  write(*,*) "" 
  write(*,'(''Elements         :'',i8)') numberOfElements  
  write(*,'(''Nodes            :'',i8)') numberOfNodes 
  write(*,'(''Boundary faces   :'',i8)') numberOfBoundaryFaces 
  write(*,'(''Hexahedra        :'',i8)') numberOfHexEl
  write(*,'(''Prisms           :'',i8)') numberOfPriEl
  write(*,'(''Pyramids         :'',i8)') numberOfPyrEl
  write(*,'(''Tetrahedra       :'',i8)') numberOfTetEl
  write(*,'(''Triangular faces :'',i8)') numberOfTriangularBFs 
  write(*,'(''Rectangular faces:'',i8)') numberOfRectangularBFs 
  write(*,'(''Surface segments :'',i8)') numberOfSurfaceSegments 
  write(*,'(''Line segments    :'',i8)') numberOfLineSegments 
  write(*,'(''Number of edges  :'',i8)') numberOfSides
  write(*,*) "" 
  
! initiates data structure for input

  allocate(grp(1)%brp,stat=allocateStatus)
  if(allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary register"

  do i=1,numberOfGrids
   grp(i)%gridNumber = i   
  end do
  grp(1)%numberOfNodes = numberOfNodes
  
  call constructBoundaryRegister(numberOfSurfaceSegments,numberOfLineSegments,&
                                 numberOfTriangularBFs,numberOfRectangularBFs,&
                                 numberOfNodes,&
                                 numberOfSurfaceNormals,grp(1)%brp)  

  allocate(hep,stat=allocateStatus)
  if(allocateStatus /= 0) STOP "ERROR: Not enough memory to create element register"
  allocate(prp,stat=allocateStatus)
  if(allocateStatus /= 0) STOP "ERROR: Not enough memory to create element register"
  allocate(pyp,stat=allocateStatus)
  if(allocateStatus /= 0) STOP "ERROR: Not enough memory to create element register"
  allocate(tep,stat=allocateStatus)
  if(allocateStatus /= 0) STOP "ERROR: Not enough memory to create element register"

  call constructHexElArray(numberOfHexEl,hep) ! calls the constructor of the current element 
                                              ! module 
  call constructPriElArray(numberOfPriEl,prp)  
  call constructPyrElArray(numberOfPyrEl,pyp)  
  call constructTetElArray(numberOfTetEl,tep)  

  allocate(crp,stat=allocateStatus)
  if(allocateStatus /= 0) STOP "ERROR: Not enough memory to create coordinate register"
  allocate(crpPrev,stat=allocateStatus)                            !MOVE
  if(allocateStatus /= 0) STOP "ERROR: Not enough memory to create coordinate register"
  allocate(crpPrev2,stat=allocateStatus)                           !MOVE
  if(allocateStatus /= 0) STOP "ERROR: Not enough memory to create coordinate register"

  call constructCoordinateRegister(numberOfNodes,NSD,crp) ! calls the constructor of the point
                                                          ! register
  call constructCoordinateRegister(numberOfNodes,NSD,crpPrev)      !MOVE
  call constructCoordinateRegister(numberOfNodes,NSD,crpPrev2)     !MOVE

! read element data 
  write(*,*) "Reading element data..." 
  if(isHybrid) then 
   call readHexElData(VOLINFILE,hep) ! calls element module to read it's own data 
   call readPriElData(VOLINFILE,prp) 
   call readPyrElData(VOLINFILE,pyp) 
   call readTetElData(VOLINFILE,tep) 
  else
   call readTetElData(VOLINFILE,tep)
  end if
! read coordinate data 
  write(*,*) "Reading coordinate data..."
  call readCoordinateData(VOLINFILE,crp) ! calls point register module to read data 


!MOVE

  if(DCOORFILE>0) then
   write(*,*) "Reading first .dcr file"
   call readCoordinateData(DCOORFILE,crpPrev)
   call setUpPreviousCoordinates(crpPrev,crp,.true.)
   close(DCOORFILE)
  else
   write(*,*) "Setting dcr1 zero "
   call setUpPreviousCoordinates(crpPrev,crp,.false.)
  end if

  if(DCOORFILE2>0) then
   write(*,*) "Reading second .dcr file"
   call readCoordinateData(DCOORFILE2,crpPrev2)
   call setUpPreviousCoordinates(crpPrev2,crpPrev,.true.)
   close(DCOORFILE2)
  else
   write(*,*) "Setting dcr2 zero"
   call setUpPreviousCoordinates(crpPrev2,crpPrev,.false.)
  end if

!MOVE


! read boundary face data 
  write(*,*) "Reading boundary data..."
  call readBoundData(numberOfBoundaryFacesSurf,SEGINFILE,VOLINFILE,grp(1)%brp,isHybrid,groundAngle) 
  close(SEGINFILE)
  close(VOLINFILE) 

  write(*,*) ""

  end subroutine readFineGrid
!--------------------------------------------------------------------------
  subroutine cleanUpElements()
  ! deallocates element data

  IMPLICIT NONE
 
  call cleanUpHexElData(hep)
  call cleanUpPriElData(prp)
  call cleanUpPyrElData(pyp)
  call cleanUpTetElData(tep)  

  end subroutine cleanUpElements
!--------------------------------------------------------------------------
  integer function nameLen(fn)
  IMPLICIT NONE
 
  character*80 :: fn
  integer :: i

  do i = 80,1,-1
   nameLen = i
   if(fn(i:i)/=' ') GOTO 77 ! EXIT
   end do 
   nameLen = 0
77  end function nameLen
!--------------------------------------------------------------------------
 subroutine membreak()
 IMPLICIT NONE
  
 character*20 :: tex
 
 write(*,*) "Type enter to continue"
 read(*,*) tex
 
 end subroutine membreak
!--------------------------------------------------------------------------

end module Agglomerator3d

