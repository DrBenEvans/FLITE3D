!**************************************
!*mgnsg3d_def_par.F90                 *
!*                                    *
!*Module definition file for program  *
!*                                    *
!*         mgnsg3d v 3.0              *
!*                                    *
!*16.08.99-15.12.2001                 *
!**************************************
!*Description:                        *
!* This file includes the single grid *
!* solver together with prolongation  *
!* and restriction operators needed   *
!* for the multigrid algorithm.       *
!**************************************

module Toolbox
  use Comms
  ! various handy datastructures

  type LinkedInteger
     integer :: int  
     type(linkedInteger),pointer :: next
  end type linkedInteger

  type LinkedIntegerP
     type(LinkedInteger),pointer :: first
  end type LinkedIntegerP

  type LinkedIntegerPArray
     type(linkedIntegerPArray),pointer :: lipa(:)
  end type linkedIntegerPArray

  type LinkedReal
     real :: re 
     type(LinkedReal),pointer :: next 
  end type LinkedReal

  type LinkedRealP
     type(LinkedReal),pointer :: first 
  end type LinkedRealP

  type LinkedRealPArray
     type(LinkedRealP),pointer :: lrpa(:)
  end type LinkedRealPArray

  type RealArrayPointers
     real,pointer :: ar(:,:)
  end type realArrayPointers

  type IntegerArrayPointers
     integer,pointer :: ar(:,:)
  end type IntegerArrayPointers

end module Toolbox

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

module InputVariables
  ! Holds the input variables specified by the user, such as the
  ! Reynolds, Mach, CFL numbers etc. 

  type InputVariablesData
     integer :: numberOfProcesses ! for parallel 
     character :: dataDirectory*120
     character :: solverCommand*200

     integer :: numberOfMGIterations
     integer :: numberOfSGIterations
     integer :: numberOfCFLIncrements
     integer :: writeToFileInterval
     real :: ReynoldsNumber,MachNumber,PrandtlNumber,CFLNumber
     real :: turbulentPrandtlNumber,turbulentCFLNumber
     real :: inflowTemperature
     real :: gamma ! ratio of specific heats
     real :: alpha ! angle of attack
     real :: beta  ! yaw angle
     real :: wallTemperature ! for isothermal wall
     real :: HartensCorrectionFactor ! used in roe matrix construction for
     ! dissipation terms
     logical :: useMatrixDissipation ! decides if matrix or scalar smoothing
     ! is to be applied
     real :: secondOrderDissipationFactor,fourthOrderDissipationFactor
     real :: coarseGridDissipationFactor ! for multigrid coarse grids
     integer :: numberOfRelaxationSteps ! number of iterations for each grid
     integer :: numberOfRSSteps ! number of residual smoothing steps
     real :: inflowField(6,2) ! inflow variables, 2d: 1-6 where 6 is pressure

     !VT: Patch initialization
     logical :: patchInitialization;
     integer :: patchType;
     real :: patchBox(3,2); ! lowerleft and upperight point
     real :: patchSphere(4); ! center and radius of the sphere 

     integer :: HLLCFluxOrder ! HLLC order
     integer :: HLLCICLS
     integer :: Automaticlipping
     integer :: GradientMethod
     real :: HLLCClipping 
     
     logical :: explicit

     real :: residualSmoothingFactor
     real :: prolongationSmoothingFactor
     integer :: numberOfPSSteps ! number of prolongation smoothing steps
     integer :: prolongationMappingScheme
     integer :: viscosityScheme
     integer :: numberOfTurbulenceSteps
     real :: turbulenceSmoothingFactor
     real :: maxIncrementFactor ! maximum percentage of unknown to be added
     real :: prolongationRelaxation

     integer :: engineFlowType   !bje jan11
     real :: enginesFrontMassFlow(2)
     real :: enginesRearParameters(7,2)
     real :: engineBCRelaxation

     real :: datum(3)   !datum for moment calculations (bje jan11)

     real :: nwheels    !number of wheels (max 4)  (bje jan11)
     real :: wheelCentre(4,3)  !coordinates of wheel centres for viscous torqe calcs
     real :: wheelDimensions(4,2)   !wheels' radius and width
     real :: groundplane       !z-coord of groundplane

     real :: tripFactor  ! introduces a possibility to reduce or increase the trip influence
     real :: maxTurbulenceValue  ! maximum allowable turbulence value 

!START-TURB
!New variables for turbulence models (Araya 2009) 
     real :: firstoffwallpoint !Distance of the first off-wall point to the wall 
     real :: kinfoverUinf2 !freestream of kinetic energy
     real :: maxTurbulenceValueK,maxTurbulenceValueW,minTurbulenceValueW,maxTurbulenceValueMU !setting limits for k-omega
     real :: winf !setting the freestream value of omega
     real :: nuwall !value of nu at wall for omega
     logical :: SAS
!END-TURB
!START-GUST
     real :: XPER(2),YPER(2),ZPER(2) !Defining the points where the perturbations are going to be prescribed
     logical :: GUSTMODEL
     real :: UDS
     integer :: typegust !Type of gust
     integer :: TTS !duration of the gust in timesteps
     integer :: nspan !number of spanwise stations for local aer. parameters
     integer :: spanCor !Coordinates of the number of spanwise stations for local aer. parameters
     real :: tipspan,rootspan !spanwise coordinates of the root and tip 
!END-GUST

     logical :: wallsAreIsentropic  ! if false, isothermal wall is used
     logical :: movingwall   !if true, moving walls are present (bje jan11)

     real :: turbProlongationRelaxation  ! prolongation relaxation for turbulence

     logical :: addTimeStep
     integer :: flowType
     integer ::  timeScheme
     integer :: numberOfStepsPerCycle
     integer :: cycleNumber
     integer :: numberofphysicaltimesteps
     integer :: maximumnumberofcycles 
     integer :: restartNumber
     real :: physicalTimeStep
     real :: stopTime
     real :: currentTime

     ! Added fluxType to define the flux type used: Vinh-Tan
     integer :: fluxType; !< 1: Central differencing; 2: LLF flux; 3: Roe flux

     real :: referenceArea  ! for lift and drag calculations
     integer :: forceCalcType    !bje jan11

     integer :: boundaryTerm ! 1: FE, 2: FV

     real :: multigridBCRelaxation  ! boundary condition relaxation on multigrid correction

     real :: minimumDensity  ! minimum allowable density

     logical :: useDissipationWeighting ! scale dissipation with inverse edge length

     integer :: numberOfTurbulenceGrids ! number of grids to use for turbulence equation

     integer :: numberOfEORelaxationSteps ! number of engine outlet relaxation steps

     integer :: numberOfTriggerSteps ! number of steps to perform triggering of turbulence
     real :: turbulenceTriggerValue 
     real :: triggerRadius ! radius from trip line where triggering is to occur

     real :: turbK   ! factor for turbulent timesteps 
     integer :: multigridScheme
     integer :: turbulenceModel
     integer :: dissipationScheme
     integer :: coarseGridDissipationScheme
     integer :: numberOfGridsToUse
     logical :: useTimeResidual ! if true, plots residual vs. time instead of iterations
     real    :: initialResidual
  end type InputVariablesData

contains

  !-----------------------------------------------------------------------
  subroutine readInputVariables(ivd,INFILE)
    ! as the name says
    IMPLICIT NONE

    type(InputVariablesData) :: ivd
    integer :: INFILE

    integer :: i

    integer :: allocateStatus
    real :: PI,alphaRad,betaRad

    namelist /InputVariables/ ivd

    ! set default values 

    ivd%numberOfMGIterations = 100 
    ivd%numberOfSGIterations = 0 
    ivd%numberOfCFLIncrements = 0
    ivd%writeToFileInterval = 25
    ivd%gamma = 1.4
    ivd%MachNumber = 0.5
    ivd%alpha = 0.0
    ivd%beta = 0.0
    ivd%CFLNumber = 1.0 
    ivd%ReynoldsNumber = 1.0E06 
    ivd%viscosityScheme = 1
!START-TURB
!The employed parameters are computed at the end of this subroutine
    ivd%inflowField(1,1) = 1.0 ! density 
    ivd%inflowField(2,1) = 1.0 ! x-momentum
    ivd%inflowField(3,1) = 0.0 ! y-momentum
    ivd%inflowField(4,1) = 0.0 ! z-momentum
    ivd%inflowField(5,1) = 1.0 ! total energy
    ivd%inflowField(6,1) = ivd%inflowField(1,1)*(ivd%MachNumber**2)/ivd%gamma ! pressure 
    ivd%inflowField(:,2) = ivd%inflowField(:,1) !Specified at the input file
!END-TURB
    ivd%inflowTemperature = 528.0
    ivd%wallsAreIsentropic = .false.
    ivd%movingwall = .false.     ! (bje jan11)
    ivd%wallTemperature = 528.0

    ivd%patchInitialization = .false.
    ivd%patchType = 1;
    ivd%patchSphere(1) = 0.0 ! x-center
    ivd%patchSphere(2) = 0.0 ! y-center
    ivd%patchSphere(3) = 0.0 ! z-center
    ivd%patchSphere(4) = 1.0 ! radius
    ivd%patchBox(1,1) = 0.0
    ivd%patchBox(2,1) = 0.0
    ivd%patchBox(3,1) = 0.0
    ivd%patchBox(1,2) = 0.0
    ivd%patchBox(2,2) = 0.0
    ivd%patchBox(3,2) = 0.0

    ivd%engineFlowType = 1  !bje jan11
    ivd%enginesFrontMassFlow = 0.0
    ivd%enginesRearParameters = 0.0
    ivd%numberOfEORelaxationSteps = 0
    ivd%engineBCRelaxation = 1.0
    ivd%dissipationScheme = 2
    ivd%coarseGridDissipationScheme = 2
    ivd%secondOrderDissipationFactor = 0.4
    ivd%fourthOrderDissipationFactor = 0.2
    ivd%coarseGridDissipationFactor = 0.75
    ivd%HartensCorrectionFactor = 0.0
    ivd%useMatrixDissipation = .false.
    ivd%residualSmoothingFactor = 0.0 
    ivd%multigridScheme = 3
    ivd%numberOfRelaxationSteps = 1 
    ivd%numberOfGridsToUse = 1 
    ivd%numberOfRSSteps = 0 
    ivd%numberOfPSSteps = 0 
    ivd%multigridBCRelaxation = 1.0
    ivd%turbulenceModel = 0
    ivd%numberOfTurbulenceGrids = 1 
    ivd%turbulentCFLNumber = 1.0
    ivd%tripFactor = 1.0
    ivd%triggerRadius = 9999999.9
    ivd%turbulenceTriggerValue = 0.0
    ivd%referenceArea = 1.0
    ivd%forceCalcType = 1        !bje jan11
    ivd%useTimeResidual = .false.

    ivd%flowType = 0
    ivd%addTimeStep = .true.
    ivd%timeScheme=2
    ivd%physicaltimestep=0.1
    ivd%numberofphysicaltimesteps=1
    ivd%maximumnumberofcycles=10000
    ivd%restartNumber = 1
    ivd%initialResidual = 10000

    ! other constants
    ! flux type
    ivd%fluxType = 1 
    ivd%HLLCFluxOrder = 2
    ivd%HLLCClipping = 1.0
    ivd%Automaticlipping = 0 ! Automatic clpping if # 0
    ivd%GradientMethod = 1   ! 1 for Galerkin method and LS else
    ivd%HLLCICLS = 1
    
    ivd%explicit = .false.
    
    ivd%stopTime = 1.0
    ivd%currentTime = 0.0

    ivd%turbProlongationRelaxation = 1.0 
    ivd%prolongationSmoothingFactor = 0.0
    ivd%prolongationmappingScheme = 1
    ivd%prolongationRelaxation = 1.0
    ivd%PrandtlNumber = 0.72 
    ivd%turbulenceSmoothingFactor = 0.002
    ivd%turbulentPrandtlNumber = 0.9
    ivd%numberOfTurbulenceSteps = 9999999
    ivd%useDissipationWeighting = .false.
    ivd%maxTurbulenceValue = 3000
    ivd%numberOfTriggerSteps = 0 
    ivd%maxIncrementFactor = 0.25
    ivd%minimumDensity = 0.1
    ivd%gamma = 1.4
    ivd%turbK = 0.1
    ivd%boundaryTerm = 1

!START-TURB
  ivd%maxTurbulenceValueK = 1.5 
  ivd%maxTurbulenceValueW = 1.0E7
  ivd%maxTurbulenceValueMU = 3000
  ivd%minTurbulenceValueW = 5.0  
  ivd%firstoffwallpoint = 1.0E-4 !for wall-condition of omega 
  ivd%kinfoverUinf2 = 1.0E-6 !freestream value of k
  ivd%winf = 5.0 !freestream value of omega
  ivd%nuwall = 1.0 !value of nu laminar at wall
  ivd%SAS = .false.
!END-TURB
!START-GUST
  ivd%XPER(1) = 0.0
  ivd%XPER(2) = 0.0
  ivd%YPER(1) = 0.0
  ivd%YPER(2) = 0.0
  ivd%ZPER(1) = 0.0
  ivd%ZPER(2) = 0.0
  ivd%GUSTMODEL = .false.
  ivd%UDS = 0.0
  ivd%typegust = 0
  ivd%TTS = 0
  ivd%nspan = 0   ! 8
  ivd%spanCor = 2   ! 8
  ivd%rootspan = 0.0
  ivd%tipspan = 15.0
!END-GUST

  ivd%groundplane = -9999999.9  ! (bje jan11)
  ivd%datum(1) = 0.0
  ivd%datum(2) = 0.0
  ivd%datum(3) = 0.0

  ivd%nwheels = 0.0              ! (bje jan11)
  ivd%wheelCentre(1,1) = 0.0
  ivd%wheelCentre(1,2) = 0.0
  ivd%wheelCentre(1,3) = 0.0
  ivd%wheelCentre(2,1) = 0.0
  ivd%wheelCentre(2,2) = 0.0
  ivd%wheelCentre(2,3) = 0.0
  ivd%wheelCentre(3,1) = 0.0
  ivd%wheelCentre(3,2) = 0.0
  ivd%wheelCentre(3,3) = 0.0
  ivd%wheelCentre(4,1) = 0.0
  ivd%wheelCentre(4,2) = 0.0
  ivd%wheelCentre(4,3) = 0.0
  ivd%wheelDimensions(1,1) = 0.0
  ivd%wheelDimensions(1,2) = 0.0
  ivd%wheelDimensions(2,1) = 0.0
  ivd%wheelDimensions(2,2) = 0.0
  ivd%wheelDimensions(3,1) = 0.0
  ivd%wheelDimensions(3,2) = 0.0
  ivd%wheelDimensions(4,1) = 0.0
  ivd%wheelDimensions(4,2) = 0.0


    ! read input variables file

    if(INFILE>0) then  
       read(INFILE,InputVariables) 
    end if

    if(ivd%flowType.eq.0) then
       ivd%maximumnumberofcycles=1000000
       ivd%addTimeStep = .false.
    end if
    
    if (ivd%explicit) then
       ivd%numberofphysicaltimesteps=1
       ivd%addTimeStep = .false.   
    end if
    
    if (ivd%fluxType.eq.4) then
    !No dissipation in HLLC flux computation
       ivd%secondOrderDissipationFactor = 0.0
       ivd%fourthOrderDissipationFactor = 0.0
    end if

    PI = 4.0*atan(1.0)
    if (.not.ivd%patchInitialization) then
     ivd%wallTemperature = (1.0/((ivd%gamma-1.0)*(ivd%MachNumber**2)))*(ivd%wallTemperature/ivd%inflowTemperature) 
!Scaling the variables
     alphaRad = ivd%alpha*PI/180.
     betaRad = ivd%beta*PI/180. 
     ivd%inflowField(1,1) = 1.0
     ivd%inflowField(2,1) = cos(alphaRad)*cos(betaRad) 
     ivd%inflowField(3,1) = sin(betaRad)
     ivd%inflowField(4,1) = sin(alphaRad) 
     ivd%inflowField(5,1) = (1./((ivd%gamma-1.)*ivd%gamma*ivd%MachNumber**2))+0.5
     ivd%inflowField(6,1) = 1./(ivd%gamma*ivd%MachNumber**2)
    end if  
      
    write(*,*) "Reading inputVariables"

    write(*,'(6(1e15.7))') (ivd%inflowField(i,1),i=1,6)
    write(*,'(6(1e15.7))') (ivd%inflowField(i,2),i=1,6)

  end subroutine readInputVariables
  !-----------------------------------------------------------------------
end module InputVariables

!-----------------------------------------------------------------------

module ProlongationOperator
  use Toolbox
  ! responsible for the coarse-to-fine mappings

  type ProlongationOperatorData
     integer :: numberOfCoarseNodes,numberOfFineNodes
     integer :: numberOfInternalCoarseNodes,numberOfInternalFineNodes
     integer,pointer :: prolongationArray(:) 

     integer :: numberOfComNodes
     integer,pointer :: comProlongationArray(:)
  end type ProlongationOperatorData

contains

  !-----------------------------------------------------------------------
  subroutine readProlongationOperatorData(pod,INFILE,ia,idp)
    IMPLICIT NONE

    type(ProlongationOperatorData) :: pod
    integer :: INFILE,ia,idp

    type(LinkedInteger),pointer :: currentInteger
    integer :: allocateStatus,ibuff,i,j,numberOfSubNodes

    ! prolongation operator (coarse to fine)

    read(INFILE) pod%numberOfInternalFineNodes
    read(INFILE) pod%numberOfInternalCoarseNodes

    if(ia.eq.0) then
       allocate(pod%prolongationArray(pod%numberOfFineNodes),stat=allocateStatus)
       if(allocateStatus/=0)&
            STOP "ERROR: Not enough memory in readProlongationOperatorData "
    end if
    pod%prolongationArray = 0
    do i=1,pod%numberOfInternalFineNodes
       read(INFILE) ibuff,pod%prolongationArray(ibuff)
       if(ibuff>pod%numberOfFineNodes) then
          write(*,*) "SSF: ",i,ibuff,pod%numberOfFineNodes
          pause
       end if
    end do

    read(INFILE) pod%numberOfComNodes

    if(ia.eq.0) then
       allocate(pod%comProlongationArray(pod%numberOfComNodes),stat=allocateStatus)
       if(allocateStatus/=0)&
            STOP "ERROR: Not enough memory in readProlongationOperatorData "
    end if

    do i=1,pod%numberOfComNodes
       read(INFILE)  pod%comProlongationArray(i)
    end do


  end subroutine readProlongationOperatorData
  !-----------------------------------------------------------------------
end module ProlongationOperator

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

module RestrictionOperator
  use Toolbox
  ! responsible for fine-to-coarse mappings

  type RestrictionOperatorData
     integer :: numberOfCoarseNodes,numberOfFineNodes
     integer :: numberOfInternalCoarseNodes,numberOfInternalFineNodes
     real,pointer :: restrictionArray(:)

     integer :: numberOfComNodes
     real,pointer :: comRestrictionArray(:)
  end type RestrictionOperatorData

contains

  !-----------------------------------------------------------------------
  subroutine readRestrictionOperatorData(rod,INFILE,ia,idp)
    IMPLICIT NONE
    type(RestrictionOperatorData) :: rod 
    integer :: INFILE,ia,idp

    integer :: allocateStatus,i,j,numberOfSubNodes,ibuff
    real :: rbuff 

    ! restriction operator (fine to coarse)
    read(INFILE) rod%numberOfInternalFineNodes
    read(INFILE) rod%numberOfInternalCoarseNodes

    if(ia.eq.0) then
       allocate(rod%restrictionArray(rod%numberOfFineNodes),stat=allocateStatus)
       if(allocateStatus/=0)&
            STOP "ERROR: Not enough memory in readRestrictionOperatorData "
    end if
    rod%restrictionArray = 0.0
    do i=1,rod%numberOfInternalFineNodes
       read(INFILE) ibuff,rod%restrictionArray(ibuff)
    end do

    read(INFILE) rod%numberOfComNodes

    if(ia.eq.0) then
       allocate(rod%comRestrictionArray(rod%numberOfComNodes),stat=allocateStatus)
       if(allocateStatus/=0)&
            STOP "ERROR: Not enough memory in readProlongationOperatorData "
    end if

    do i=1,rod%numberOfComNodes
       read(INFILE)  rod%comRestrictionArray(i)
    end do


  end subroutine readRestrictionOperatorData
  !-----------------------------------------------------------------------
end module RestrictionOperator

!-----------------------------------------------------------------------

module GridParallelization
  ! contains the datastructure for parallelization

  type GridParallelizationData
     integer :: currentDomain
     integer :: currentTID
     integer :: numberOfProcesses
     integer :: noInternalCommunicationNodes
     integer :: noBoundaryCommunicationNodes
     integer :: sendFlag
     logical :: isMaster

     integer,pointer :: processorIDs(:)

     integer,pointer :: numberOfSendNodes(:)
     integer,pointer :: numberOfReceiveNodes(:)
     integer,pointer :: numberOfBoundarySendNodes(:)
     integer,pointer :: numberOfBoundaryReceiveNodes(:)

     integer,pointer :: sendNodeRegister(:,:)
     integer,pointer :: receiveNodeRegister(:,:)
     integer,pointer :: boundarySendNodeRegister(:,:)
     integer,pointer :: boundaryReceiveNodeRegister(:,:)

     integer :: numberOfComNodes,numberOfComSNodes,numberOfComRNodes
     integer,pointer :: numberOfComReceiveNodes(:)
     integer,pointer :: numberOfComSendNodes(:)
     integer,pointer :: comReceiveRegister(:,:)
     integer,pointer :: comSendRegister(:,:)
     real,pointer :: comBuffer(:) 

     real,pointer :: buffer(:)

     integer,pointer :: receivedMessageFromProcess(:)

     integer :: integerType,realType,maxComm
     integer :: numberOfComSides,numberOfBoundaryComSides

     character*120 :: dataDirectory
     integer :: dDL
  end type GridParallelizationData

contains

  !-----------------------------------------------------------------------
  subroutine readGridParallelizationData(INFILE,pdp,ia)
    IMPLICIT NONE
    ! reads the parallel communication arrays from file
    integer :: INFILE,ia
    type(GridParallelizationData) :: pdp


    integer :: allocateStatus,i,k,maxSend,maxRecieve,maxComm

    if(pdp%numberOfProcesses.ne.1) then

       if(ia.eq.0) then

          allocate(pdp%numberOfSendNodes(pdp%numberOfProcesses),stat=allocateStatus)
          if(allocateStatus/=0) STOP "ERROR: Not enough memory in readGridParallelizationData"

          allocate(pdp%numberOfReceiveNodes(pdp%numberOfProcesses),stat=allocateStatus)
          if(allocateStatus/=0) STOP "ERROR: Not enough memory in readGridParallelizationData"

       endif

       read(INFILE) pdp%numberOfComSides,pdp%numberOfBoundaryComSides

       write(*,*) "Number of communication sides: ",pdp%numberOfComSides
       write(*,*) "Number of boundary communication sides: ",pdp%numberOfBoundaryComSides

       maxSend = 0
       do i=1,pdp%numberOfProcesses
          read(INFILE) pdp%numberOfSendNodes(i)
          if(pdp%numberOfSendNodes(i)>maxSend) maxSend =  pdp%numberOfSendNodes(i)
       end do


       if(ia.eq.0) then
          allocate(pdp%sendNodeRegister(pdp%numberOfProcesses,maxSend),stat=allocateStatus)
          if(allocateStatus/=0) STOP "ERROR: Not enough memory in readGridParallelizationData"
       end if
       pdp%sendNodeRegister = 0

       do i=1,pdp%numberOfProcesses
          do k=1,pdp%numberOfSendNodes(i)
             read(INFILE) pdp%sendNodeRegister(i,k)
          end do
       end do

       maxRecieve = 0
       do i=1,pdp%numberOfProcesses
          read(INFILE) pdp%numberOfReceiveNodes(i)
          if(pdp%numberOfReceiveNodes(i)>maxRecieve) maxRecieve = pdp%numberOfReceiveNodes(i)
       end do


       if(ia.eq.0) then
          allocate(pdp%receiveNodeRegister(pdp%numberOfProcesses,maxRecieve),stat=allocateStatus)
          if(allocateStatus/=0) STOP "ERROR: Not enough memory in readGridParallelizationData"
       end if

       pdp%receiveNodeRegister = 0

       do i=1,pdp%numberOfProcesses
          do k=1,pdp%numberOfReceiveNodes(i)
             read(INFILE) pdp%receiveNodeRegister(i,k)
          end do
       end do

       maxComm = max(maxSend,maxRecieve)
       pdp%maxComm = maxComm  !KS

       if(ia.eq.0) then
          allocate(pdp%buffer(pdp%maxComm),stat=allocateStatus)
          if(allocateStatus/=0) STOP "ERROR: Not enough memory in readGridParallelizationData"
       end if

    else

       read(INFILE) pdp%numberOfComSides,pdp%numberOfBoundaryComSides

       if(ia.eq.0) then
          allocate(pdp%numberOfSendNodes(pdp%numberOfProcesses),stat=allocateStatus)
          if(allocateStatus/=0) STOP "ERROR: Not enough memory in readGridParallelizationData"
       end if

       read(INFILE) pdp%numberOfSendNodes(1)

       if(ia.eq.0) then
          allocate(pdp%numberOfReceiveNodes(pdp%numberOfProcesses),stat=allocateStatus)
          if(allocateStatus/=0) STOP "ERROR: Not enough memory in readGridParallelizationData"
       end if

       read(INFILE) pdp%numberOfReceiveNodes(1) 

       if(ia.eq.0) then
          nullify(pdp%sendNodeRegister)
          nullify(pdp%receiveNodeRegister)
       end if
    end if


  end subroutine readGridParallelizationData
  !-----------------------------------------------------------------------
  subroutine readComParallelizationData(INFILE,pdp,ia)
    ! reads the parallelization datastructure used for the intergrid mapings (global agglomeration)
    IMPLICIT NONE

    integer :: INFILE,ia
    type(GridParallelizationData) :: pdp

    integer :: allocateStatus,i,j

    read(INFILE) pdp%numberOfComSNodes

    if(ia.eq.0) then
       allocate(pdp%numberOfComSendNodes(pdp%numberOfProcesses),stat=allocateStatus)
       if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComParallelizationData"
       allocate(pdp%comSendRegister(pdp%numberOfProcesses,pdp%numberOfComSNodes),stat=allocateStatus)
       if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComParallelizationData"
    end if

    pdp%comSendRegister = 0

    do i=1,pdp%numberOfProcesses
       read(INFILE) pdp%numberOfComSendNodes(i)
       do j=1,pdp%numberOfComSendNodes(i)
          read(INFILE) pdp%comSendRegister(i,j)
       end do
    end do

    read(INFILE) pdp%numberOfComRNodes

    if(ia.eq.0) then
       allocate(pdp%numberOfComReceiveNodes(pdp%numberOfProcesses),stat=allocateStatus)
       if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComParallelizationData"
       allocate(pdp%comReceiveRegister(pdp%numberOfProcesses,pdp%numberOfComRNodes),stat=allocateStatus)
       if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComParallelizationData"
    end if

    pdp%comReceiveRegister = 0

    do i=1,pdp%numberOfProcesses
       read(INFILE) pdp%numberOfComReceiveNodes(i)
       do j=1,pdp%numberOfComReceiveNodes(i)
          read(INFILE) pdp%comReceiveRegister(i,j)
       end do
    end do

    pdp%numberOfComNodes = max(pdp%numberOfComSNodes,pdp%numberOfComRNodes) 

    if(ia.eq.0) then
       allocate(pdp%comBuffer(pdp%numberOfComNodes),stat=allocateStatus)
       if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComParallelizationData"
    end if


  end subroutine readComParallelizationData
  !-----------------------------------------------------------------------
end module GridParallelization


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

module BoundarySolver
  ! contains BC resters etc
  use InputVariables 

  type BoundarySolverData
     integer,pointer :: sideIndexArray(:,:)
     integer,pointer :: nodeIndicatorRegister(:)
     integer,pointer :: nodeIndicatorArray(:)
     integer,pointer :: trailingEdgeIndexes(:)
     integer,pointer :: engineInletSideIndexes(:,:)
     real,pointer :: sideWeightsArray(:,:)
     real,pointer :: nodeNormalArray(:,:)
     real,pointer :: sideWeightNorms(:)
     real,pointer :: engineInletSideCoefficients(:,:)
     real,pointer :: normtan(:,:)

     real,pointer :: GCLWeightsArray(:,:)

     real,pointer :: sideLengthArray(:)

     integer,pointer :: IORegister(:) ! for internal outflow 
     integer :: numberOfIONodes

     real :: engineInletAreas(2)
     real :: engineInletNormals(2,3)

     integer :: numberOfDoubleNodes,numberOfTripleNodes

     integer :: numberOfBoundarySides,numberOfBoundaryNodes
     integer :: numberOfTrailingEdges
     integer :: numberOfEngineInletSides
  end type BoundarySolverData

contains

  !-----------------------------------------------------------------------
  subroutine readBoundarySolverData(ivd,bsd,INFILE,isInitial,ia)
    IMPLICIT NONE
    type(BoundarySolverData) :: bsd
    type(InputVariablesData) :: ivd
    integer :: INFILE,ia
    logical :: isInitial

    integer :: allocateStatus,i,j,k,ib,nodeNumber,dummy,numberOfNormals
    integer :: test
    real :: r1(3),r2(3),norm,GCLWeightsPrev(2)


    read(INFILE) bsd%numberOfBoundarySides
    write(*,*) ' number of boundary sides=', bsd%numberOfBoundarySides

    if(ia.eq.0) then
       allocate(bsd%sideIndexArray(bsd%numberOfBoundarySides,3),stat=allocateStatus)
       if(allocateStatus/=0) STOP "ERROR: Not enough memory in readBoundarySolverData"
       allocate(bsd%sideWeightsArray(bsd%numberOfBoundarySides,3),stat=allocateStatus)
       if(allocateStatus/=0) STOP "ERROR: Not enough memory in readBoundarySolverData"
       allocate(bsd%GCLWeightsArray(bsd%numberOfBoundarySides,2),stat=allocateStatus)
       if(allocateStatus/=0) STOP "ERROR: Not enough memory in readBoundarySolverData"
    end if

    do i=1,bsd%numberOfBoundarySides
       read(INFILE) bsd%sideIndexArray(i,1:3),bsd%sideWeightsArray(i,1:3),bsd%GCLWeightsArray(i,1:2),GCLWeightsPrev(1:2)
       if(ivd%timeScheme==2.and.isInitial) then
          bsd%GCLWeightsArray(i,1:2) = 1.5*bsd%GCLWeightsArray(i,1:2)-0.5*GCLWeightsPrev(1:2)
       end if
    end do

    read(INFILE) bsd%numberOfBoundaryNodes,numberOfNormals
    write(*,*) ' nymber of boundary nodes=', bsd%numberOfBoundarynodes,numberOfNormals

    if(ia.eq.0) then
       allocate(bsd%nodeIndicatorRegister(-20:0),stat=allocateStatus)
       if(allocateStatus/=0) STOP "ERROR: Not enough memory in readBoundarySolverData"
    end if
    do i=-20,0
       read(INFILE) bsd%nodeIndicatorRegister(i)
    end do

    if(ia.eq.0) then
       allocate(bsd%nodeIndicatorArray(bsd%numberOfBoundaryNodes),stat=allocateStatus)
       if(allocateStatus/=0) STOP "ERROR: Not enough memory in readBoundarySolverData"
!START-TURB
       allocate(bsd%normtan(numberOfNormals,7),stat=allocateStatus)
!END-TURB
       if(allocateStatus/=0) STOP "ERROR: Not enough memory in readBoundarySolverData"
    end if


    do i=1,numberOfNormals
     if(ivd%movingwall==.true.)then
       read(INFILE) bsd%nodeIndicatorArray(i),(bsd%normtan(i,ib),ib=1,6)
     else
       read(INFILE) bsd%nodeIndicatorArray(i)
     endif
       if(bsd%nodeIndicatorArray(i)>bsd%numberOfBoundaryNodes) then 
          write(*,*) "axc: ",i,bsd%nodeIndicatorArray(i),bsd%numberOfBoundaryNodes
       end if
    end do
    read(INFILE) bsd%numberOfEngineInletSides
    write(*,*) 'Engine:', bsd%numberOfEngineInletSides
    if(ia.eq.0) then
       allocate(bsd%engineInletSideIndexes(bsd%numberOfEngineInletSides,3),stat=allocateStatus)
       if(allocateStatus/=0) STOP "ERROR: Not enough memory in readBoundarySolverData"
       allocate(bsd%engineInletSideCoefficients(bsd%numberOfEngineInletSides,3),stat=allocateStatus)
       if(allocateStatus/=0) STOP "ERROR: Not enough memory in readBoundarySolverData"
    end if

    do i=1,bsd%numberOfEngineInletSides
       read(INFILE) bsd%engineInletSideIndexes(i,1:3),bsd%engineInletSideCoefficients(i,1:3)
    end do
  end subroutine readBoundarySolverData
  !-----------------------------------------------------------------------
end module BoundarySolver

!-----------------------------------------------------------------------

module GridSolver
  use Toolbox
  use RestrictionOperator
  use ProlongationOperator
  use BoundarySolver
  use InputVariables 
  use GridParallelization

  ! basically a one-grid solver with the added functionality of 
  ! intergrid mappings. 

  type GridSolverData
     real,pointer :: sideWeightsArray(:,:)
     integer,pointer :: sideIndexArray(:,:)       
     real,pointer :: sideLengthArray(:)
     real,pointer :: nodeHelpArray(:,:) ! node-based help array 
     real,pointer :: uprev(:,:) ! unknowns at previous time step
     real,pointer :: u(:,:)     ! unknowns at current time step
     real,pointer :: uold(:,:) ! uknowns at previous time step
     real,pointer :: uold2(:,:) ! uknowns at previous previous time step
     real,pointer :: dissipation(:,:)
     real,pointer :: laminarViscosity(:)
     real,pointer :: vorticity(:)
     real,pointer :: divergence(:)
     real,pointer :: coordinates(:,:) ! the fine mesh coordinates
     real,pointer :: nodeConnectivityArray(:) ! number of sides per node 
     real,pointer :: rhs(:,:) ! right hand side in calculation
     real,pointer :: sourceTerm(:,:) ! used in FAS multigrid
     real,pointer :: p(:) ! pressure field
     real,pointer :: nodeVolume(:) ! the nodal volume (lumped mass matrix) 
     real,pointer :: nodeVolumeP(:) ! at previous time step 
     real,pointer :: nodeVolumeP2(:) ! st previous previous time step
     real,pointer :: localTimeSteps(:) ! time steps used
     real,pointer :: localTurbulentTimeSteps(:) ! local time steps for turbulent equation
     real,pointer :: wallDistance(:) ! used in turbulence modeling
     real,pointer :: wallArea(:) ! used in calculating lift and drag
     real,pointer :: wallStress(:,:) ! used in lift and drag calculations

     real,pointer :: coordinateMovement(:,:)

     real,pointer :: GCLWeightsArray(:) ! ALE flux for geometric conservation

     integer,pointer :: wallDistanceBoundaryNodeArray(:) ! for turbulence modeling
     integer,pointer :: seedNodeRegister(:) 

     real,pointer :: mapComBuffer(:,:) ! buffer for global agglomeration strategy

!START-TURB
     real,pointer :: turbulenceDiffusionTermK(:),turbulenceDiffusionTermW(:)
     real,pointer :: SK(:),SW(:)
     real,pointer :: dkdw(:),F1(:),F2(:),TVR(:),imsr(:), msrt(:)
     real,pointer :: D1U(:,:),D2U(:,:)
     real,pointer  :: sas(:,:) !auxiliary array for SAS variables
!END-TURB
     INTEGER,pointer :: IPER(:)!START-GUST gust indexes
     INTEGER      :: GUSTPAR(3)!START-GUST gust coordinates
     real :: localAerParam(10,8)! START-GUST auxiliary array for local aerodynamic parameters


     real,pointer :: coefficientSideLengthArray(:)

     real,pointer :: turbulenceDiffusionTerm(:)

     integer,pointer :: tripNodeFieldIndexes(:,:)
     real,pointer :: tripNodeFieldDistances(:)
     integer,pointer :: tripLineIndexes(:)
     real,pointer :: tripWallLength(:)
     integer :: numberOfTripNodes
     integer :: numberOfTripFieldNodes

     type(RestrictionOperatorData) :: rod
     type(ProlongationOperatorData) :: pod 
     type(BoundarySolverData) :: brp
     type(GridParallelizationData) :: pdp

     integer :: numberOfSides,numberOfNodes
     real :: residualScalingFactor
     real,pointer :: RKCoefficients(:) ! for timestepping
     integer :: gridNumber,numberOfSeparationPoints,sizeOfSeparationField
     integer :: cycleNumber 
     real :: physicalTimeStep
     real :: physicalTime

     real :: startTime,stopTime

  end type GridSolverData

  integer :: position,bufferSize          
  parameter (bufferSize=2000000)            
  character :: bufferzz(bufferSize)          
  character :: bufferoo(bufferSize)          

contains

  !-----------------------------------------------------------------------
  subroutine doTimeIterations(grp,ivd,relaxationSteps,calculateDissipation,calculateTurbulence,iterationNumber,physicalTimestepNumber)!START-GUST
    IMPLICIT NONE

    type(GridSolverData) :: grp
    type(InputVariablesData) :: ivd
    integer :: relaxationSteps
    logical :: calculateDissipation,calculateTurbulence
    integer :: iterationNumber

    integer :: i ,ii
    integer :: physicalTimestepNumber !START-GUST

    ! time loop
    
    grp%cycleNumber = iterationNumber 
  
    do i=1,relaxationSteps
       grp%uprev(:,1) = grp%u(:,1)
       grp%uprev(:,2) = grp%u(:,2)
       grp%uprev(:,3) = grp%u(:,3)
       grp%uprev(:,4) = grp%u(:,4)
       grp%uprev(:,5) = grp%u(:,5)
!START-TURB
  if(ivd%turbulenceModel==1) then
  grp%uprev(:,6) = grp%u(:,6)
 end if
 if(ivd%turbulenceModel==2 .or. ivd%turbulenceModel==3) then
  grp%uprev(:,6) = grp%u(:,6)
  grp%uprev(:,7) = grp%u(:,7) 
 end if
!END-TURB
       call singleIteration(grp,ivd,calculateDissipation,calculateTurbulence,iterationNumber,physicalTimestepNumber)!START-GUST
    end do

  end subroutine doTimeIterations
  !-----------------------------------------------------------------------
  subroutine singleIteration(grp,ivd,calculateDissipation,calculateTurbulence,iterationNumber,physicalTimestepNumber)!START-GUST
    ! main driving subroutine for solver at a given grid level
    IMPLICIT NONE

    type(GridSolverData) :: grp
    type(InputVariablesData) :: ivd
    logical :: calculateDissipation,calculateTurbulence
    integer :: iterationNumber
!START-TURB
    real :: mv,mva(7)
    integer :: i,j,k,ind,mvi,mvia(7),ist,ien,ib,ip
!END-TURB

    integer :: ii

    integer :: maxind
    real :: maxincr,maxTurb
    integer :: physicalTimestepNumber !START-GUST

    !write(grp%pdp%currentDomain+500,'(i7,6E15.7)') (k,grp%u(k,:),k=1,grp%numberOfNodes)
    if(calculateDissipation) then 
       ! make pressure field 
       call makePressureField(grp,ivd,grp%u)
       grp%p = abs(grp%p)
       ! make local timesteps
       call makeTimeSteps(grp,ivd)

       ! make artificial dissipation
       if(ivd%dissipationScheme==1.and.grp%gridNumber==1) then 
          call makeArtificialDissipation1(grp,ivd)
       else if(ivd%dissipationScheme==3) then 
          call makeArtificialDissipation3(grp,ivd)
       else
          call makeArtificialDissipation2(grp,ivd)
       end if

       if(ivd%ReynoldsNumber>0.0) then  ! means viscous flow
          if(ivd%viscosityScheme==2) then
             call makeViscosityTerm2(grp,ivd,grp%dissipation,.true.)
          else
             call makeViscosityTerm(grp,ivd,grp%dissipation,.true.)
          end if
          if(ivd%turbulenceModel==1.and.calculateTurbulence) then
             call makeTurbulentTimeSteps(grp,ivd)
             call makeSADiffusion(grp,ivd)
          end if
!START-TURB
          if(ivd%turbulenceModel==2.and.calculateTurbulence) then !MOD
             call makeTurbulentTimeStepsKW(grp,ivd)
             call makeKWDiffusion(grp,ivd)    
          end if
          if(ivd%turbulenceModel==3.and.calculateTurbulence) then !SST
             call makeTurbulentTimeStepsKW(grp,ivd)
             if(ivd%SAS)call makeSASd2u(grp,ivd)
             call makeSSTdkdw(grp,ivd)
             call makeSSTDiffusion(grp,ivd)
             call makeSSTCrossDiff(grp,ivd) !Cross-diffusion    
          end if
!END-TURB
       end if
    end if

    do i=1,3  ! Runge-Kutta iterations 
     ! make pressure field 
     call makePressureField(grp,ivd,grp%u)

       ! make RHS
     grp%rhs = 0.0

     call makeRHS(grp,ivd)

     call addTimeDifferences(grp,ivd)

!START-TURB    
     if(ivd%turbulenceModel==1.and.calculateTurbulence) then
       if(iterationNumber<ivd%numberOfTriggerSteps.and.grp%gridNumber==1) then 
         call triggerTurbulenceField(grp,ivd)
       end if
     end if
     if(ivd%turbulenceModel==1.and.calculateTurbulence) then
       call makeSARHS(grp,ivd)
       call makeSASourceTerm(grp,ivd)
       if(grp%gridNumber==1) then 
         call makeSATripTerm(grp,ivd)
       end if
     end if
     if(ivd%turbulenceModel==2.and.calculateTurbulence) then !MOD
       call makeKWRHS(grp,ivd)
       call makeKWSourceTerm(grp,ivd)
     end if
     if(ivd%turbulenceModel==3.and.calculateTurbulence) then 
       call makeSSTRHS(grp,ivd)
       call makeSSTSourceTerm(grp,ivd)
     end if
!END-TURB

 ! multiply by timestep and invert mass matrix 

     grp%rhs(:,1) = grp%rhs(:,1) + grp%dissipation(:,1)
     grp%rhs(:,2) = grp%rhs(:,2) + grp%dissipation(:,2)
     grp%rhs(:,3) = grp%rhs(:,3) + grp%dissipation(:,3)
     grp%rhs(:,4) = grp%rhs(:,4) + grp%dissipation(:,4)
     grp%rhs(:,5) = grp%rhs(:,5) + grp%dissipation(:,5)
     if(ivd%turbulenceModel==1.and.calculateTurbulence) then
       grp%rhs(:,6) = grp%rhs(:,6) + grp%turbulenceDiffusionTerm(:)
     end if
!START-TURB
     if(ivd%turbulenceModel==2.and.calculateTurbulence) then !MOD
       grp%rhs(:,6) = grp%rhs(:,6) + grp%turbulenceDiffusionTermK(:)
       grp%rhs(:,7) = grp%rhs(:,7) + grp%turbulenceDiffusionTermW(:)
     end if
     if(ivd%turbulenceModel==3.and.calculateTurbulence) then !SST
       grp%rhs(:,6) = grp%rhs(:,6) + grp%turbulenceDiffusionTermK(:)
       grp%rhs(:,7) = grp%rhs(:,7) + grp%turbulenceDiffusionTermW(:)
     end if
!END-TURB
       
     call solveSystem(grp,ivd,i,grp%rhs)
       
    ! fix boundary conditions on increment field
     call setBCsOnIncrementField(grp,ivd,grp%rhs)
    ! smooth residual if wanted
     if(ivd%residualSmoothingFactor>0.0) then 
       do j=1,ivd%numberOfRSSteps
        call smoothResidual2(grp,ivd,grp%rhs,ivd%residualSmoothingFactor)
       end do
       call setBCsOnIncrementField(grp,ivd,grp%rhs)
     end if
 
    ! nullify increment at outer boundary for coarse grid 
     if(grp%gridNumber>1) then
       call nullifyOuterBoundary(grp,grp%rhs)
     end if

    ! update unknown field
     call updateSolutionVector(grp%u,grp%uprev,grp%rhs,grp,ivd,i)
    ! fix boundary conditions on solution field
     grp%GUSTPAR(1) = iterationNumber !START-GUST
     grp%GUSTPAR(2) = i !START-GUST
     call setBCsOnSolutionField(grp,ivd,grp%u,physicalTimestepNumber)!START-GUST
    ! add multigrid source term
     if(grp%gridNumber>1) then 
      ! multigrid source term
       grp%u(:,1) = grp%u(:,1) - grp%RKCoefficients(i)*grp%sourceTerm(:,1)
       grp%u(:,2) = grp%u(:,2) - grp%RKCoefficients(i)*grp%sourceTerm(:,2)
       grp%u(:,3) = grp%u(:,3) - grp%RKCoefficients(i)*grp%sourceTerm(:,3)
       grp%u(:,4) = grp%u(:,4) - grp%RKCoefficients(i)*grp%sourceTerm(:,4)
       grp%u(:,5) = grp%u(:,5) - grp%RKCoefficients(i)*grp%sourceTerm(:,5)
!START-TURB
       if(ivd%turbulenceModel==1.and.calculateTurbulence) then 
         grp%u(:,6) = grp%u(:,6) - grp%RKCoefficients(i)*grp%sourceTerm(:,6)
       end if
       if(ivd%turbulenceModel==2.and.calculateTurbulence) then 
         grp%u(:,6) = grp%u(:,6) - grp%RKCoefficients(i)*grp%sourceTerm(:,6)
         grp%u(:,7) = grp%u(:,7) - grp%RKCoefficients(i)*grp%sourceTerm(:,7)
       end if
       if(ivd%turbulenceModel==3.and.calculateTurbulence) then 
         grp%u(:,6) = grp%u(:,6) - grp%RKCoefficients(i)*grp%sourceTerm(:,6)
         grp%u(:,7) = grp%u(:,7) - grp%RKCoefficients(i)*grp%sourceTerm(:,7)
       end if
!END-TURB
     end if
       
    ! make sure that density is above a user--defined value
     call limitDensity(grp,grp%u(:,1),ivd)
       
    ! make sure that the turbulent viscosity is within user--defined bounds
     if (calculateTurbulence) then
!START-TURB
       if(ivd%turbulenceModel==1) then
         call limitTurbulence(grp,grp%u(:,6),ivd)
       end if
     end if
!END-TURB

       ! fix boundary conditions at outer boundary
     call setOuterBoundaryConditions(grp,ivd,grp%u,grp%p)
    end do

  end subroutine singleIteration
  !-----------------------------------------------------------------------
  
  subroutine addTimeDifferences(grp,ivd)
    ! add time differences to fine grid residuals

    IMPLICIT NONE

    type(GridSolverData) :: grp
    type(InputVariablesData) :: ivd

    integer :: i

!START-TURB
 if(ivd%addTimeStep) then
  if(grp%gridNumber==1.and.ivd%timeScheme==2)then
   do i=1,grp%numberOfNodes
    if(ivd%turbulenceModel==0) then
      grp%rhs(i,1:5) = grp%rhs(i,1:5) + &
      (1.5*grp%nodeVolume(i)*grp%u(i,1:5) - 2.0*grp%nodeVolumeP(i)*grp%uold(i,1:5) &
      + 0.5*grp%nodeVolumeP2(i)*grp%uold2(i,1:5))/grp%physicalTimestep
    end if
    if(ivd%turbulenceModel==1) then
      grp%rhs(i,1:6) = grp%rhs(i,1:6) + &
      (1.5*grp%nodeVolume(i)*grp%u(i,1:6) - 2.0*grp%nodeVolumeP(i)*grp%uold(i,1:6) &
      + 0.5*grp%nodeVolumeP2(i)*grp%uold2(i,1:6))/grp%physicalTimestep
    end if
    if(ivd%turbulenceModel==2) then
      grp%rhs(i,1:7) = grp%rhs(i,1:7) + &
      (1.5*grp%nodeVolume(i)*grp%u(i,1:7) - 2.0*grp%nodeVolumeP(i)*grp%uold(i,1:7) &
      + 0.5*grp%nodeVolumeP2(i)*grp%uold2(i,1:7))/grp%physicalTimestep
    end if
    if(ivd%turbulenceModel==3) then
      grp%rhs(i,1:7) = grp%rhs(i,1:7) + &
      (1.5*grp%nodeVolume(i)*grp%u(i,1:7) - 2.0*grp%nodeVolumeP(i)*grp%uold(i,1:7) &
      + 0.5*grp%nodeVolumeP2(i)*grp%uold2(i,1:7))/grp%physicalTimestep
    end if
   end do
  else
   do i=1,grp%numberOfNodes
    if(ivd%turbulenceModel==0) then
      grp%rhs(i,1:5) = grp%rhs(i,1:5) + grp%nodeVolume(i)*               &
                      (grp%u(i,1:5) - grp%uold(i,1:5))/grp%physicalTimestep
    end if
    if(ivd%turbulenceModel==1) then
      grp%rhs(i,1:6) = grp%rhs(i,1:6) + grp%nodeVolume(i)*               &
                      (grp%u(i,1:6) - grp%uold(i,1:6))/grp%physicalTimestep
    end if
    if(ivd%turbulenceModel==2) then
      grp%rhs(i,1:7) = grp%rhs(i,1:7) + grp%nodeVolume(i)*               &
                      (grp%u(i,1:7) - grp%uold(i,1:7))/grp%physicalTimestep
    end if
    if(ivd%turbulenceModel==3) then
      grp%rhs(i,1:7) = grp%rhs(i,1:7) + grp%nodeVolume(i)*               &
                      (grp%u(i,1:7) - grp%uold(i,1:7))/grp%physicalTimestep
    end if
   end do
  end if
 end if
!END-TURB
  end subroutine addTimeDifferences
  !-----------------------------------------------------------------------
  subroutine limitDensity(grp,u,ivd)
    IMPLICIT NONE
    ! ensures that the density is above a certain value

    type(GridSolverData) :: grp
    real :: u(:)
    type(InputVariablesData) :: ivd

    integer :: i

    do i=1,grp%numberOfNodes
       u(i) = max(u(i),ivd%minimumDensity)
    end do

  end subroutine limitDensity
  !-----------------------------------------------------------------------
  subroutine limitTurbulence(grp,u,ivd)
    IMPLICIT NONE
    ! keeps the density above zero and below a specified value

    type(GridSolverData) :: grp
    real :: u(:)
    type(InputVariablesData) :: ivd

    integer :: i
    integer :: maxTurbInd
    real :: maxTurb

    maxTurb = 0.0
    maxTurbInd = 0

    do i=1,grp%numberOfNodes
       u(i) = max(u(i),0.0) 
       if(grp%u(i,6)>ivd%maxTurbulenceValue) then 
          if(grp%u(i,6)>maxTurb) then
             maxTurb = grp%u(i,6)
             maxTurbInd = i
          end if
          grp%u(i,6) = ivd%maxTurbulenceValue
       end if
    end do

    if(maxTurbInd>0.and.grp%gridNumber==1) then
       write(*,*) "Warning: large turbulence: ",maxTurb,maxTurbInd,grp%coordinates(maxTurbInd,:)
    end if


  end subroutine limitTurbulence
  !-----------------------------------------------------------------------
  subroutine updateSolutionVector(u,uprev,delu,grp,ivd,RKIteration)
    ! updates solution increment, relaxes if update is to large
    IMPLICIT NONE 

    real :: u(:,:),uprev(:,:),delu(:,:)
    type(GridSolverData) :: grp
    type(InputVariablesData) :: ivd
    integer :: RKIteration

    integer :: i,maxNode
    real :: increment,sign,maxIncr
    logical :: hasRelaxed
    real :: dtJ

    if(RKIteration > 0) then
       dtJ = ivd%CFLNumber*grp%RKCoefficients(RKIteration)
    else
       dtJ = ivd%CFLNumber
    end if

    delu(:,1) = dtJ*delu(:,1)
    delu(:,2) = dtJ*delu(:,2)
    delu(:,3) = dtJ*delu(:,3)
    delu(:,4) = dtJ*delu(:,4)
    delu(:,5) = dtJ*delu(:,5)


!START-TURB
    if(ivd%turbulenceModel > 0) then
       if(RKIteration > 0) then
          dtJ = ivd%turbulentCFLNumber*grp%RKCoefficients(RKIteration)
       else
          dtJ = ivd%turbulentCFLNumber
       end if

       if(ivd%turbulenceModel==1) then 
          delu(:,6) = dtJ*delu(:,6)
       end if
       if(ivd%turbulenceModel==2) then 
          delu(:,6) = dtJ*delu(:,6)
          delu(:,7) = dtJ*delu(:,7)
       end if
       if(ivd%turbulenceModel==3) then 
          delu(:,6) = dtJ*delu(:,6)
          delu(:,7) = dtJ*delu(:,7)
       end if
!END-TURB
    end if

    maxIncr = -1.0
    hasRelaxed = .false.
    do i=1,grp%numberOfNodes
       ! relax density if necessary
       increment = delu(i,1)
       if(abs(increment)>maxIncr) then 
          if(abs(increment+1)<0.1) then 
!            write(*,*) "Kaal: ",i,increment
          end if
          maxIncr = increment
          maxNode = i
       end if
       if(abs(increment/uprev(i,1))<ivd%maxIncrementFactor) then 
          u(i,1) = uprev(i,1)  - increment
       else 
          sign = increment/abs(increment)
          u(i,1) = uprev(i,1)  - sign*ivd%maxIncrementFactor*uprev(i,1)
          hasRelaxed = .true.
       end if
       ! velocities are not relaxed
       u(i,2) = uprev(i,2)  - delu(i,2)
       u(i,3) = uprev(i,3)  - delu(i,3)
       u(i,4) = uprev(i,4)  - delu(i,4)

       ! total energy not relaxed
       u(i,5) = uprev(i,5)  - delu(i,5)

 !START-TURB
      if(ivd%turbulenceModel==1) then 
       u(i,6) = uprev(i,6) - delu(i,6)
       u(i,6) = max(u(i,6),0.0)
      end if

      if(ivd%turbulenceModel==2) then 
       u(i,6) = uprev(i,6) - delu(i,6)
       u(i,7) = uprev(i,7) - delu(i,7)
       u(i,6) = max(u(i,6),0.0)!zero value if negative
       u(i,6) = min(u(i,6),ivd%maxTurbulenceValueK)
! Maximum condition for omega: stress limiter (Wilcox)
! By using the mean strain rate tensor
       u(i,7) = max(u(i,7),7./8.*sqrt(2*grp%msrt(i)*100./9.))
       u(i,7) = max(u(i,7),ivd%minTurbulenceValueW)
       u(i,7) = min(u(i,7),ivd%maxTurbulenceValueW)
      end if

      if(ivd%turbulenceModel==3) then 
       u(i,6) = uprev(i,6) - delu(i,6)
       u(i,7) = uprev(i,7) - delu(i,7)
       u(i,6) = max(u(i,6),0.0)!zero value if negative
       u(i,6) = min(u(i,6),ivd%maxTurbulenceValueK)
       u(i,7) = max(u(i,7),ivd%minTurbulenceValueW)
       u(i,7) = min(u(i,7),ivd%maxTurbulenceValueW)
      end if
!END-TURB

    end do

    do i=1,grp%numberOfNodes
       if(grp%u(i,1)<ivd%minimumDensity) then 
          if(grp%pdp%currentDomain==1) write(*,*) "Small density: ",i,grp%u(i,1),grp%gridNumber
          grp%u(i,1) = ivd%minimumDensity
       end if
    end do

  end subroutine updateSolutionVector
  !-----------------------------------------------------------------------
  subroutine solveSystem(grp,ivd,RKIteration,rhs)
    ! inverts mass matrix to get unknown from RHS
    IMPLICIT NONE

    type(GridSolverData) :: grp
    type(InputVariablesData) :: ivd
    integer :: RKIteration
    real :: rhs(:,:)

    integer :: i
    real :: dt,dtJ,timeDiffCoeff
    ! multiply with time increments and divide by lumped mass

    if(RKIteration > 0) then 
       dtJ = ivd%CFLNumber*grp%RKCoefficients(RKIteration)
    else
       dtJ = ivd%CFLNumber
    end if

    timeDiffCoeff = 0.0
    if(ivd%flowType.ne.0) timeDiffCoeff = 1.5

    do i=1,grp%numberOfNodes
       if (ivd%flowType.ne.0) then 
 if (ivd%explicit) then
  dt = ivd%physicalTimeStep
 else
  dt = grp%localTimeSteps(i)
  dt = dt/(1.0+timeDiffCoeff*dt/grp%physicalTimeStep)
 end if
       else
 dt = grp%localTimeSteps(i)
       end if
       rhs(i,1:5) = dt*rhs(i,1:5)/grp%nodeVolume(i)  
    end do

    if(ivd%turbulenceModel==1) then
       call solveTurbulenceSystem(grp,ivd,RKIteration,rhs(:,6))
    end if
!START-TURB
    if(ivd%turbulenceModel==2) then 
     do i=1,grp%numberOfNodes
       if (ivd%flowType.ne.0) then 
 if (ivd%explicit) then
  dt = ivd%physicalTimeStep
 else
  dt = grp%localTurbulentTimeSteps(i)
 end if
       else
 dt = grp%localTurbulentTimeSteps(i)
       end if
      rhs(i,6) = (dt/(1.0+0.5*grp%SK(i)*dt))*rhs(i,6)/grp%nodeVolume(i)/grp%u(i,1)
      rhs(i,7) = (dt/(1.0+0.5*grp%SW(i)*dt))*rhs(i,7)/grp%nodeVolume(i)/grp%u(i,1)
     end do
    end if
   if(ivd%turbulenceModel==3) then
    do i=1,grp%numberOfNodes
       if (ivd%flowType.ne.0) then 
 if (ivd%explicit) then
  dt = ivd%physicalTimeStep
 else
  dt = grp%localTurbulentTimeSteps(i)
 end if
       else
 dt = grp%localTurbulentTimeSteps(i)
       end if
     rhs(i,6) = (dt/(1.0+0.5*grp%SK(i)*dt))*rhs(i,6)/grp%nodeVolume(i)/grp%u(i,1)
     rhs(i,7) = (dt/(1.0+0.5*grp%SW(i)*dt))*rhs(i,7)/grp%nodeVolume(i)/grp%u(i,1)
    end do
   end if
!END-TURB

  end subroutine solveSystem
  !-----------------------------------------------------------------------
  subroutine solveTurbulenceSystem(grp,ivd,RKIteration,rhs)
    ! inverts mass matrix to get unknown from RHS
    IMPLICIT NONE

    type(GridSolverData) :: grp
    type(InputVariablesData) :: ivd
    integer :: RKIteration
    real :: rhs(:)

    integer :: i
    real :: dt,dtJ
!START-TURB
    real :: timeDiffCoeff
!END-TURB

    ! multiply with time increments and divide by lumped mass

    if(RKIteration > 0) then
       dtJ = ivd%turbulentCFLNumber*grp%RKCoefficients(RKIteration)
    else
       dtJ = ivd%turbulentCFLNumber
    end if

!START-TURB
    timeDiffCoeff = 0.0
    if(ivd%flowType.ne.0) timeDiffCoeff = 1.5

    do i=1,grp%numberOfNodes
       if (ivd%flowType.ne.0) then 
 if (ivd%explicit) then
  dt = ivd%physicalTimeStep
 else
  dt = grp%localTurbulentTimeSteps(i)
 end if
       else
 dt = grp%localTurbulentTimeSteps(i)
       end if
       rhs(i) = dt*rhs(i)/grp%nodeVolume(i)
    end do
!END-TURB

  end subroutine solveTurbulenceSystem
  !-----------------------------------------------------------------------

  subroutine makeRHS(grp,ivd)
    ! makes right hand side by utilizing a side-based structure
    ! with coefficients created in preprocessor

    IMPLICIT NONE

    include 'mpif.h'

    type(GridSolverData) :: grp
    type(InputVariablesData) :: ivd

    integer :: i,j,ind1,ind2,faceIndicator,ii

    real :: pr1,pr2,rinv1,rinv2,ru1,rv1,rw1,ru2,rv2,rw2,u1,v1,u2,v2,w1,w2,rh1,rh2
    real :: f1,f2,f3,f4,f5,wx,wy,wz,rkx,rky,rkz,rnorm,unorm1,unorm2,sc

    real :: rnx,rny,rnz,rho1,normalVelocity1
    real :: rho2,normalVelocity2,vml,aml
    real :: ro,ro1,uo,vo,eo,wo,po,unr,rhor,uxr,vyr,wzr,epsr,presr,hr
    real :: f11,f21,f31,f41,f51,f12,f22,f32,f42,f52,eps1,eps1l,eps2l
    real :: fr11,fr21,fr31,fr41,fr51,fr12,fr22,fr32,fr42,fr52
    real :: fl11,fl21,fl31,fl41,fl51,fl12,fl22,fl32,fl42,fl52
    real :: di,d1,ui,vi,wi,hi,ci2,ci,af,ucp,rh1l,rh2l,h1l,h2l,h1,h2
    real :: du1,du2,du3,du4,du5,rlam1,rlam2,rlam3,rlam,gam1,epslm
    real :: s1,s2,al1x,al2x,cc1,cc2,f1l,f2l,f3l,f4l,f5l
    real :: pr1l,rho1l,rinv1l,ru1l,rv1l,rw1l,u1l,v1l,w1l,normalVelocity1l
    real :: pr2l,rho2l,rinv2l,ru2l,rv2l,rw2l,u2l,v2l,w2l,normalVelocity2l
    integer :: ind
    real :: wgcl

    real :: r1,E1,r2,E2,pu1,pv1,pw1,pu2,pv2,pw2
    real :: c1,c2,rlam11, rlam12, rlam21, rlam22


    integer :: kk,istart,ifinish,autoClip,GradientCompute

    !For HLLC flux!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real ::roe1,roe2,re1,re2,rroe,utild,vtild,wtild,qtild,htild,ctild,c,sl,sr,sm
    real ::fact,fl1,fl2,fl3,fl4,fl5,fr1,fr2,fr3,fr4,fr5,plstar,prstar
    real ::factl,factr,roelstar,roerstar,m1lstar,m2lstar,m3lstar,elstar
    real ::m1rstar,m2rstar,m3rstar,erstar
    real ::flstar1,flstar2,flstar3,flstar4,flstar5,frstar1,frstar2,frstar3,frstar4,frstar5,rnorm1
    real ::dx,dy,dz,ulstar,vlstar,wlstar,urstar,vrstar,wrstar,pstar
    real :: beta0,hsize,hessgrad,gmin,gmax
    real :: dr1,dv1,dw1,de1,dr2,dv2,dw2,de2,her1,her2,her3,herder,gradl,gradr,hlength,clip
    real :: dr11,du11,dv11,dw11,de11,dr22,du22,dv22,dw22,de22,drl,dul,dvl,dwl,del,drr,dur,dvr,dwr,der

    integer :: hllcflux
    
    real :: a,b,a1,bmin,bmax,dq,invmax,invmin,dpv1,dpv11
  integer :: nsig
  real :: dpvl(5)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    gam1 = ivd%gamma-1.0
    epslm = 0.
    eps1 = 1./(max(epslm,1.e-5))

    ! initialize
    !Compute the gradient of primitive variables
    !if the Hllc second order solver is selected 
    !L. Remaki
    if (ivd%fluxType.eq.4) then
       hllcflux = ivd%HLLCFluxOrder
       clip = ivd%HLLCClipping
       autoClip = ivd%Automaticlipping
       GradientCompute = ivd%GradientMethod
       if(ivd%ReynoldsNumber<1) then
           nsig =4  
       else
           nsig =2
       end if 
       if( hllcflux .gt. 1) then
          !Computing gradient of physical variables
          if(GradientCompute .eq. 1) then
            call computeGradient(grp,ivd,grp%u)  ! Using Galerkin method
          else
            call LSComputeGradient(grp,ivd,grp%u)  ! Using Leats-Square method
          end if
   !Using automatic local clipping (LED-TVD Scheme)    
          if( autoClip .eq. 1)  call LocalMinMax(grp,ivd,grp%u)
       end if
    end if


#ifdef PARALLEL

    do kk=1,2
       if(kk==1) then
          istart = 1
          ifinish = grp%pdp%numberOfComSides
       else
          istart = grp%pdp%numberOfComSides + 1
          ifinish = grp%numberOfSides
       end if


       do i=istart,ifinish

#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL

  ind1 = grp%sideIndexArray(i,1)
  ind2 = grp%sideIndexArray(i,2)
  
  pr1 = grp%p(ind1)             ! pressure
  r1 = grp%u(ind1,1)
  rinv1 = 1./grp%u(ind1,1)      ! rho inverse
  ru1 = grp%u(ind1,2)           ! x-momentum
  rv1 = grp%u(ind1,3)           ! y-momentum
  rw1 = grp%u(ind1,4)           ! z-momentum
  u1 = rinv1*ru1                ! x-velocity
  v1 = rinv1*rv1                ! y-velocity
  w1 = rinv1*rw1                ! z-velocity
  E1 = grp%u(ind1,5)
  rh1 = grp%u(ind1,5) + pr1     ! enthalpy
   
  pr2 = grp%p(ind2)
  r2 = grp%u(ind2,1)
  rinv2 = 1./grp%u(ind2,1)
  ru2 = grp%u(ind2,2)
  rv2 = grp%u(ind2,3)
  rw2 = grp%u(ind2,4)
  u2 = rinv2*ru2
  v2 = rinv2*rv2
  w2 = rinv2*rw2
  E2 = grp%u(ind2,5)
  rh2 = grp%u(ind2,5) + pr2 

! side weights from preprocessor 
  wx = grp%sideWeightsArray(i,1)
  wy = grp%sideWeightsArray(i,2)  
  wz = grp%sideWeightsArray(i,3)  

  wgcl = grp%GCLWeightsArray(i)/grp%physicalTimestep

  if(ivd%fluxType==1) then
     f1 = wx*(ru1+ru2)+wy*(rv1+rv2)+wz*(rw1+rw2) - (r1+r2)*wgcl
     f2 = wx*(ru1*u1+pr1+ru2*u2+pr2) + wy*(ru1*v1+ru2*v2)+wz*(ru1*w1+ru2*w2)-(ru1+ru2)*wgcl
     f3 = wx*(rv1*u1+rv2*u2) + wy*(rv1*v1+pr1+rv2*v2+pr2)+wz*(rv1*w1+rv2*w2)-(rv1+rv2)*wgcl
     f4 = wx*(rw1*u1+rw2*u2) + wy*(rw1*v1+rw2*v2)+wz*(rw1*w1+pr1+rw2*w2+pr2)-(rw1+rw2)*wgcl
     f5 = wx*(rh1*u1+rh2*u2) + wy*(rh1*v1+rh2*v2) + wz*(rh1*w1+rh2*w2) - (E1+E2)*wgcl
  else if(ivd%fluxType==2) then !LLF flux
     rnorm=sqrt(wx*wx+wy*wy+wz*wz)
     if(rnorm>1.0e-10) then
        rnx = wx/rnorm
        rny = wy/rnorm
        rnz = wz/rnorm

        ! Compute Fn1
        h1  = rh1*rinv1
        c1  = sqrt(ivd%gamma*pr1*rinv1)
        normalVelocity1 = u1*rnx+v1*rny+w1*rnz

        f11 = rho1*normalVelocity1
        f21 = pr1*rnx+ru1*normalVelocity1
        f31 = pr1*rny+rv1*normalVelocity1
        f41 = pr1*rnz+rw1*normalVelocity1
        f51 = rh1*normalVelocity1

        ! Compute Fn2
        h2  = rh2*rinv2
        c2  = sqrt(ivd%gamma*pr2*rinv2)
        normalVelocity2 = u2*rnx+v2*rny+w2*rnz

        f12 = rho2*normalVelocity2
        f22 = pr2*rnx+ru2*normalVelocity2
        f32 = pr2*rny+rv2*normalVelocity2
        f42 = pr2*rnz+rw2*normalVelocity2
        f52 = rh2*normalVelocity2

        du1   = r1-r2
        du2   = ru1-ru2
        du3   = rv1-rv2
        du4   = rw1-rw2
        du5   = E1-E2

        rlam11 = abs(normalVelocity1+c1)
        rlam12 = abs(normalVelocity2+c2)
        rlam1 = rlam11
        if(rlam1.lt.rlam12) rlam1=rlam12

        rlam21 = abs(normalVelocity1-c1)
        rlam22 = abs(normalVelocity2-c2)
        rlam2 = rlam21
        if(rlam2.lt.rlam22) rlam2=rlam22

        rlam = rlam1
        if(rlam.lt.rlam2) rlam = rlam2

        f1   = rnorm*((f11+f12) - rlam*du1) - (r1+r2)*wgcl
        f2   = rnorm*((f21+f22) - rlam*du2) - (ru1+ru2)*wgcl
        f3   = rnorm*((f31+f32) - rlam*du3) - (rv1+rv2)*wgcl
        f4   = rnorm*((f41+f42) - rlam*du4) - (rw1+rw2)*wgcl
        f5   = rnorm*((f51+f52) - rlam*du5) - (E1+E2)*wgcl

     end if ! if(rnorm>1.0e-10)
  else if (ivd%fluxType.eq.3) then ! Roe flux
     rnorm=sqrt(wx*wx+wy*wy+wz*wz)
     if(rnorm>1.0e-10) then
        rnx = wx/rnorm
        rny = wy/rnorm
        rnz = wz/rnorm

        ! Compute Fn1
        u1  = ru1*rinv1
        v1  = rv1*rinv1
        w1  = rw1*rinv1
        normalVelocity1 = u1*rnx+v1*rny+w1*rnz

        f11 = rho1*normalVelocity1
        f21 = pr1*rnx+ru1*normalVelocity1
        f31 = pr1*rny+rv1*normalVelocity1
        f41 = pr1*rnz+rw1*normalVelocity1
        f51 = rh1*normalVelocity1

        ! Compute Fn2

        normalVelocity2 = u2*rnx+v2*rny+w2*rnz

        f12 = rho2*normalVelocity2
        f22 = pr2*rnx+ru2*normalVelocity2
        f32 = pr2*rny+rv2*normalVelocity2
        f42 = pr2*rnz+rw2*normalVelocity2
        f52 = rh2*normalVelocity2

        ! Compute chractersitic
        di    = sqrt(r1/r2)
        d1    = 1.0/(di+1.0)
        ui    = (di*u1+u2)*d1
        vi    = (di*v1+v2)*d1
        wi    = (di*w1+w2)*d1
        hi    = (di*h1+h2)*d1
        ci2   = gam1*(hi-0.5*(ui*ui+vi*vi+wi*wi))
        ci2   = max(ci2,1.0e-5)
        ci    = sqrt(ci2)
        af    = 0.5*(ui*ui+vi*vi+wi*wi)
        ucp   = ui*rnx+vi*rny+wi*rnz

        du1   = r1-r2
        du2   = ru1-ru2
        du3   = rv1-rv2
        du4   = rw1-rw2
        du5   = E1-E2

        rlam1 = abs(ucp+ci)
        rlam2 = abs(ucp-ci)
        rlam3 = abs(ucp)

        if(rlam1.lt.epslm) rlam1 = 0.5*(rlam1*rlam1*eps1+epslm)
        if(rlam2.lt.epslm) rlam2 = 0.5*(rlam2*rlam2*eps1+epslm)
        if(rlam3.lt.epslm) rlam3 = 0.5*(rlam3*rlam3*eps1+epslm)

        s1    = 0.5*(rlam1+rlam2)
        s2    = 0.5*(rlam1-rlam2)
        al1x  = gam1*(af*du1-ui*du2-vi*du3-wi*du4+du5)
        al2x  = -ucp*du1+du2*rnx+du3*rny+du4*rnz
        cc1   = ((s1-rlam3)*al1x/ci2)+(s2*al2x/ci)
        cc2   = (s2*al1x/ci)+(s1-rlam3)*al2x
        f1   = rnorm*((f11+f12)-(rlam3*du1+cc1           )) - (r1+r2)*wgcl
        f2   = rnorm*((f21+f22)-(rlam3*du2+cc1*ui+cc2*rnx)) - (ru1+ru2)*wgcl
        f3   = rnorm*((f31+f32)-(rlam3*du3+cc1*vi+cc2*rny)) - (rv1+rv2)*wgcl
        f4   = rnorm*((f41+f42)-(rlam3*du4+cc1*wi+cc2*rnz)) - (rw1+rw2)*wgcl
        f5   = rnorm*((f51+f52)-(rlam3*du5+cc1*hi+cc2*ucp)) - (E1+E2)*wgcl
     end if ! if(rnorm>1.0e-10)
     !HLLC flux
  else if (ivd%fluxType.eq.4) then
     !CASE (4)
     ! HLLC Riemann Solver to Approximate invisid Fluxes
     ! L.Remaki
     f1 = wx*(ru1+ru2)+wy*(rv1+rv2)+wz*(rw1+rw2)
     f2 = wx*(ru1*u1+pr1+ru2*u2+pr2) + wy*(ru1*v1+ru2*v2)+wz*(ru1*w1+ru2*w2)
     f3 = wx*(rv1*u1+rv2*u2) + wy*(rv1*v1+pr1+rv2*v2+pr2)+wz*(rv1*w1+rv2*w2)
     f4 = wx*(rw1*u1+rw2*u2) + wy*(rw1*v1+rw2*v2)+wz*(rw1*w1+pr1+rw2*w2+pr2)
     f5 = wx*(rh1*u1+rh2*u2) + wy*(rh1*v1+rh2*v2) + wz*(rh1*w1+rh2*w2)
     roe1 = grp%u(ind1,1)
     roe2 =  grp%u(ind2,1)
     re1 = grp%u(ind1,5)
     re2 = grp%u(ind2,5)
        
     wx = 2.0*wx
     wy = 2.0*wy
     wz = 2.0*wz
     rnorm=sqrt(wx*wx+wy*wy+wz*wz)
     rnorm1 =rnorm
     if(rnorm > 0.0) rnorm = 1.0/rnorm
     rnx = wx*rnorm
     rny = wy*rnorm
     rnz = wz*rnorm

     ! Add second order terms
     if( hllcflux .gt. 1) then

        normalvelocity1 = u1*rnx+v1*rny+w1*rnz
        normalvelocity2 = u2*rnx+v2*rny+w2*rnz

        !Compute the acoustic waves sl, sr and the contact wave sm

        beta0 = sqrt((ivd%gamma - 1)/(2*ivd%gamma))
        rroe = sqrt(roe2/roe1)
        utild = (u1 + u2*rroe)/(1 + rroe)
        vtild = (v1 + v2*rroe)/(1 + rroe)
        wtild = (w1 + w2*rroe)/(1 + rroe)
        qtild = utild*rnx+vtild*rny+wtild*rnz
        htild = (rh1*rinv1 + rh2*rinv2*rroe)/(1 + rroe)
        ctild = (ivd%gamma - 1)*(htild - 0.5*(utild*utild + vtild*vtild + wtild*wtild))
        ctild = sqrt(abs(ctild)) 
        c1 = sqrt((ivd%gamma - 1)*abs(rh1*rinv1- 0.5*(u1*u1 + v1*v1 + w1*w1)))
        c2 = sqrt((ivd%gamma - 1)*abs(rh2*rinv2- 0.5*(u2*u2 + v2*v2 + w2*w2)))
        sl = amin1(normalvelocity1 - beta0*c1, qtild - ctild)
        sr = amax1(normalvelocity2 + beta0*c2, qtild + ctild)

        dx = grp%coordinates(ind2,1) - grp%coordinates(ind1,1)
        dy = grp%coordinates(ind2,2) - grp%coordinates(ind1,2)
        dz = grp%coordinates(ind2,3) - grp%coordinates(ind1,3)
        hsize = 0.5   
        hlength = 0.5*clip*sqrt(dx*dx +dy*dy +dz*dz)
        
        !limiting derivatives
        do j=1,5
         dpv1=hsize*(grp%nodeHelpArray(ind1,j)*dx + grp%nodeHelpArray(ind1,j+5)*dy + grp%nodeHelpArray(ind1,j+10)*dz)
         dpv11=hsize*(grp%nodeHelpArray(ind2,j)*dx + grp%nodeHelpArray(ind2,j+5)*dy + grp%nodeHelpArray(ind2,j+10)*dz)

     !Automatic clipping
         if (autoClip .eq. 1) then 
          a=dpv1
          bmin=grp%nodeHelpArray(ind1,j+15)
          bmax=grp%nodeHelpArray(ind1,j+20)   
          c = a
          a1 = 2.0*a
          if(a1 .gt. 0.0) then
           invmax =0.0
           if(bmax .gt. 0.0) invmax =1.0/bmax
           c= 0.5*bmax*((a1*invmax))/((1.+(a1*invmax)**nsig)**(1.0/nsig))
     
          else if(a1 .lt. 0.0) then
           invmin =0.0
           if(bmin .lt. 0.0) invmin =1.0/bmin
           c= 0.5*bmin*((a1*invmin))/((1.+(a1*invmin)**nsig)**(1.0/nsig))
          end if
         
          dpv1 = c    
      
          bmin=grp%nodeHelpArray(ind2,j+15)
          bmax=grp%nodeHelpArray(ind2,j+20)
          a=-dpv11   
          c = a
          a1 = 2.0*a
          if(a1 .gt. 0.0) then
           invmax =0.0
           if(bmax .gt. 0.0) invmax =1.0/bmax
           c= 0.5*bmax*((a1*invmax))/((1.+(a1*invmax)**nsig)**(1.0/nsig))
     
          else if(a1 .lt. 0.0) then
           invmin =0.0
           if(bmin .lt. 0.0) invmin =1.0/bmin
           c= 0.5*bmin*((a1*invmin))/((1.+(a1*invmin)**nsig)**(1.0/nsig))
          end if
         
          dpv11 = -c   
     
         end if !end if (autoClip)
         
         !Using Minmod function to ensure conservation
         a = dpv1
         b = dpv11
         if( (a*b .le. 0.0)  ) then
          c = 0.0
         else
          if (abs(a) .le. abs(b)) c = a
          if (abs(a) .gt. abs(b)) c = b 
          if((abs(c) .gt. hlength) .and. (hlength .gt. 0.0)) then
           c= hlength*c/abs(c)
          end if
         end if
         dpvl(j) = c
     
        end do
        

        dr1 = dpvl(1)
        du1 = dpvl(2)
        dv1 = dpvl(3)
        dw1 = dpvl(4)
        de1 = dpvl(5)
        dr2 = dpvl(1)
        du2 = dpvl(2)
        dv2 = dpvl(3)
        dw2 = dpvl(4)
        de2 = dpvl(5)
        
        !Corrected variables (add hight order terms)
        if( (roe1+ dr1) .gt. 0.1 ) roe1 = roe1+ dr1
        u1   = u1  + du1
        v1   = v1  + dv1
        w1   = w1  + dw1 
        if( (re1 + de1) .gt. 0.0001 ) re1  = re1 + de1
        pr1  = (ivd%gamma - 1)*(re1 - 0.5*(u1*u1 + v1*v1 + w1*w1)*roe1)
        pr1  = max(pr1, 0.1)
        rh1  = re1 + pr1
        rinv1= 1/roe1
        ru1  = roe1*u1
        rv1  = roe1*v1
        rw1  = roe1*w1

        if( (roe2 - dr2) .gt. 0.1 ) roe2 = roe2- dr2
        u2   = u2 - du2
        v2   = v2 - dv2
        w2   = w2 - dw2
        if( (re2 - de2) .gt. 0.0001 ) re2  = re2 - de2
        pr2  = (ivd%gamma - 1)*(re2 - 0.5*(u2*u2 + v2*v2 + w2*w2)*roe2)
        pr2  = max(pr2, 0.1)
        rh2  = re2 + pr2
        rinv2= 1/roe2
        ru2  = roe2*u2
        rv2  = roe2*v2
        rw2  = roe2*w2    

     end if !if (hllcflux.gt.1)  

     normalvelocity1 = u1*rnx+v1*rny+w1*rnz
     normalvelocity2 = u2*rnx+v2*rny+w2*rnz

     !Compute the acoustic waves sl, sr and the contact wave sm

     beta0 = sqrt((ivd%gamma - 1)/(2*ivd%gamma))
     rroe = sqrt(roe2/roe1)
     utild = (u1 + u2*rroe)/(1 + rroe)
     vtild = (v1 + v2*rroe)/(1 + rroe)
     wtild = (w1 + w2*rroe)/(1 + rroe)
     qtild = utild*rnx+vtild*rny+wtild*rnz
     htild = (rh1*rinv1 + rh2*rinv2*rroe)/(1 + rroe)
     ctild = (ivd%gamma - 1)*(htild - 0.5*(utild*utild + vtild*vtild + wtild*wtild))
     ctild = sqrt(abs(ctild)) 
     c1 = sqrt((ivd%gamma - 1)*abs(rh1*rinv1- 0.5*(u1*u1 + v1*v1 + w1*w1)))
     c2 = sqrt((ivd%gamma - 1)*abs(rh2*rinv2- 0.5*(u2*u2 + v2*v2 + w2*w2)))
     sl = amin1(normalvelocity1 - beta0*c1, qtild - ctild)
     sr = amax1(normalvelocity2 + beta0*c2, qtild + ctild)
     sl = amin1(sl,normalvelocity2 - beta0*c2)
     sr = amax1(sr,normalvelocity1 + beta0*c1)
     sm = roe2*normalvelocity2*(sr - normalvelocity2) - roe1*normalvelocity1*(sl - normalvelocity1)
     sm = sm + pr1 - pr2
     fact = roe2*(sr - normalvelocity2) - roe1*(sl - normalvelocity1)
     if(fact .ne. 0.0 ) fact = 1.0/fact
     sm = fact*sm
  
 
     ! Right and Left Fluxes
     fl1 = rnx*(ru1)+rny*(rv1)+rnz*(rw1)       
     fl2 = rnx*(ru1*u1+pr1) + rny*(ru1*v1) + rnz*(ru1*w1) 
     fl3 = rnx*(rv1*u1) + rny*(rv1*v1+pr1) + rnz*(rv1*w1)
     fl4 = rnx*(rw1*u1) + rny*(rw1*v1) + rnz*(rw1*w1+pr1)
     fl5 = rnx*(rh1*u1) + rny*(rh1*v1) + rnz*(rh1*w1)
     fr1 = rnx*(ru2)+rny*(rv2)+rnz*(rw2)   
     fr2 = rnx*(ru2*u2+pr2) + rny*(ru2*v2) + rnz*(ru2*w2) 
     fr3 = rnx*(rv2*u2) + rny*(rv2*v2+pr2) + rnz*(rv2*w2)
     fr4 = rnx*(rw2*u2) + rny*(rw2*v2) + rnz*(rw2*w2+pr2)
     fr5 = rnx*(rh2*u2) + rny*(rh2*v2) + rnz*(rh2*w2)

     !HLLC Riemann problem solution
     plstar = roe1*(normalvelocity1 - sl)*(normalvelocity1 - sm) + pr1
     prstar = roe2*(normalvelocity2 - sr)*(normalvelocity2 - sm) + pr2 
     factl = (sl - sm)
     if(factl .ne. 0.0) then 
     factl = 1.0/factl
     roelstar = roe1*(sl - normalvelocity1)*factl
     m1lstar =( (sl - normalvelocity1)*roe1*u1 + (plstar - pr1)*rnx )*factl
     m2lstar =( (sl - normalvelocity1)*roe1*v1 + (plstar - pr1)*rny )*factl
     m3lstar =( (sl - normalvelocity1)*roe1*w1 + (plstar - pr1)*rnz )*factl
     elstar = ( (sl - normalvelocity1)*re1 - pr1*normalvelocity1 + plstar*sm )*factl
     else
     roelstar = roe1
     m1lstar  = ru1
     m2lstar  = rv1
     m3lstar  = rw1
     elstar   = re1
     end  if
     factr = (sr - sm)
     if(factr .ne. 0.0) then
     factr = 1.0/factr   
     roerstar = roe2*(sr - normalvelocity2)*factr
     m1rstar =( (sr - normalvelocity2)*roe2*u2 + (prstar - pr2)*rnx )*factr
     m2rstar =( (sr - normalvelocity2)*roe2*v2 + (prstar - pr2)*rny )*factr
     m3rstar =( (sr - normalvelocity2)*roe2*w2 + (prstar - pr2)*rnz )*factr
     erstar = ( (sr - normalvelocity2)*re2 - pr2*normalvelocity2 + prstar*sm )*factr
     else
     roerstar = roe2
     m1rstar = ru2
     m2rstar = rv2
     m3rstar = rw2
     erstar  = re2 
     end if 

     flstar1 = fl1 + sl*(roelstar - roe1)
     flstar2 = fl2 + sl*(m1lstar - ru1)
     flstar3 = fl3 + sl*(m2lstar - rv1)
     flstar4 = fl4 + sl*(m3lstar - rw1)
     flstar5 = fl5 + sl*(elstar - re1)

     frstar1 = fr1 + sr*(roerstar - roe2)
     frstar2 = fr2 + sr*(m1rstar - ru2)
     frstar3 = fr3 + sr*(m2rstar - rv2)
     frstar4 = fr4 + sr*(m3rstar - rw2)
     frstar5 = fr5 + sr*(erstar - re2)

     !HLLC fluxes at t = 0
     if( sl .gt. 0.0) then
        f1 = fl1*rnorm1
        f2 = fl2*rnorm1
        f3 = fl3*rnorm1
        f4 = fl4*rnorm1
        f5 = fl5*rnorm1
     else if((sl .le. 0.0) .and. (sm .gt. 0.0)) then
        f1 = flstar1*rnorm1
        f2 = flstar2*rnorm1
        f3 = flstar3*rnorm1
        f4 = flstar4*rnorm1
        f5 = flstar5*rnorm1
  r1  = roelstar
  ru1 = m1lstar
  rv1 = m2lstar
  rw1 = m3lstar
  E1  = elstar 
     else if((sm .le. 0.0) .and. (sr .ge. 0.0)) then
        f1 = frstar1*rnorm1
        f2 = frstar2*rnorm1
        f3 = frstar3*rnorm1
        f4 = frstar4*rnorm1
        f5 = frstar5*rnorm1
  r1  = roerstar
  ru1 = m1rstar
  rv1 = m2rstar
  rw1 = m3rstar
  E1  = erstar 
     else if( sr .lt. 0.0) then
        f1 = fr1*rnorm1
        f2 = fr2*rnorm1
        f3 = fr3*rnorm1
        f4 = fr4*rnorm1
        f5 = fr5*rnorm1
  r1  = r2 
  ru1 = ru2
  rv1 = rv2
  rw1 = rw2
  E1  = E2 
     end if
     !UPDATE FLUX FOR ALE GCL
    ! f1 = f1 - (r1+r2)*wgcl
    ! f2 = f2-(ru1+ru2)*wgcl
    ! f3 = f3-(rv1+rv2)*wgcl
    ! f4 = f4-(rw1+rw2)*wgcl
    ! f5 = f5 - (E1+E2)*wgcl
      ! Use HLLC solution
      f1 = f1 - r1*wgcl
      f2 = f2 - ru1*wgcl
      f3 = f3 - rv1*wgcl
      f4 = f4 - rw1*wgcl
      f5 = f5 - E1*wgcl
  end if !end fluxtype options
  !END SELECT
 
    
  ! assemble right hand side
  grp%rhs(ind1,1) = grp%rhs(ind1,1) + f1
  grp%rhs(ind1,2) = grp%rhs(ind1,2) + f2 
  grp%rhs(ind1,3) = grp%rhs(ind1,3) + f3  
  grp%rhs(ind1,4) = grp%rhs(ind1,4) + f4 
  grp%rhs(ind1,5) = grp%rhs(ind1,5) + f5 
  grp%rhs(ind2,1) = grp%rhs(ind2,1) - f1 
  grp%rhs(ind2,2) = grp%rhs(ind2,2) - f2
  grp%rhs(ind2,3) = grp%rhs(ind2,3) - f3 
  grp%rhs(ind2,4) = grp%rhs(ind2,4) - f4 
  grp%rhs(ind2,5) = grp%rhs(ind2,5) - f5 
 
end do 

#ifdef PARALLEL
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfBoundaryComSides
 else
  istart = grp%pdp%numberOfBoundaryComSides + 1
  ifinish = grp%brp%numberOfBoundarySides
 end if
 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL

 ind1  = grp%brp%sideIndexArray(i,1)
 ind2  = grp%brp%sideIndexArray(i,2)
 faceIndicator = grp%brp%sideIndexArray(i,3)

! add GCL terms

 r1 = grp%u(ind1,1)
 rinv1 = 1./grp%u(ind1,1)      ! rho inverse
 ru1 = grp%u(ind1,2)           ! x-momentum
 rv1 = grp%u(ind1,3)           ! y-momentum
 rw1 = grp%u(ind1,4)           ! z-momentum
 u1 = rinv1*ru1
 v1 = rinv1*rv1
 w1 = rinv1*rw1
 E1 = grp%u(ind1,5)

 r2 = grp%u(ind2,1)
 rinv2 = 1./grp%u(ind2,1)
 ru2 = grp%u(ind2,2)
 rv2 = grp%u(ind2,3)
 rw2 = grp%u(ind2,4)
 u2 = rinv2*ru2
 v2 = rinv2*rv2
 w2 = rinv2*rw2
 E2 = grp%u(ind2,5)

 if(ivd%boundaryTerm==1) then
  wgcl = grp%brp%GCLWeightsArray(i,1)/grp%physicalTimestep
  
  grp%rhs(ind1,1) = grp%rhs(ind1,1) - (3.*r1+r2)*wgcl
  grp%rhs(ind1,2) = grp%rhs(ind1,2) - (3.*ru1+ru2)*wgcl
  grp%rhs(ind1,3) = grp%rhs(ind1,3) - (3.*rv1+rv2)*wgcl
  grp%rhs(ind1,4) = grp%rhs(ind1,4) - (3.*rw1+rw2)*wgcl
  grp%rhs(ind1,5) = grp%rhs(ind1,5) - (3.*E1+E2)*wgcl
  
  wgcl = grp%brp%GCLWeightsArray(i,2)/grp%physicalTimestep
  
  grp%rhs(ind2,1) = grp%rhs(ind2,1) - (r1+3.*r2)*wgcl
  grp%rhs(ind2,2) = grp%rhs(ind2,2) - (ru1+3.*ru2)*wgcl
  grp%rhs(ind2,3) = grp%rhs(ind2,3) - (rv1+3.*rv2)*wgcl
  grp%rhs(ind2,4) = grp%rhs(ind2,4) - (rw1+3.*rw2)*wgcl
  grp%rhs(ind2,5) = grp%rhs(ind2,5) - (E1+3.*E2)*wgcl
 else
  wgcl = 4.0*grp%brp%GCLWeightsArray(i,1)/grp%physicalTimestep
  
  grp%rhs(ind1,1) = grp%rhs(ind1,1) - r1*wgcl
  grp%rhs(ind1,2) = grp%rhs(ind1,2) - ru1*wgcl
  grp%rhs(ind1,3) = grp%rhs(ind1,3) - rv1*wgcl
  grp%rhs(ind1,4) = grp%rhs(ind1,4) - rw1*wgcl
  grp%rhs(ind1,5) = grp%rhs(ind1,5) - E1*wgcl

  wgcl = 4.0*grp%brp%GCLWeightsArray(i,2)/grp%physicalTimestep

  grp%rhs(ind2,1) = grp%rhs(ind2,1) - r2*wgcl
  grp%rhs(ind2,2) = grp%rhs(ind2,2) - ru2*wgcl
  grp%rhs(ind2,3) = grp%rhs(ind2,3) - rv2*wgcl
  grp%rhs(ind2,4) = grp%rhs(ind2,4) - rw2*wgcl
  grp%rhs(ind2,5) = grp%rhs(ind2,5) - E2*wgcl

 end if

 ! Boundary conditions
 if(faceIndicator.le.2.or.faceIndicator.eq.9) then
  ! symmetry or wall

  wx  = grp%brp%sideWeightsArray(i,1)
  wy  = grp%brp%sideWeightsArray(i,2)
  wz  = grp%brp%sideWeightsArray(i,3)
  rnorm=sqrt(wx*wx+wy*wy+wz*wz)
  if(rnorm>1.0e-10) then
   rnx = wx/rnorm
   rny = wy/rnorm
   rnz = wz/rnorm

   pr1 = grp%p(ind1) ! pressure
   rho1 = grp%u(ind1,1)
   rinv1  = 1./rho1
   ru1 = grp%u(ind1,2)
   rv1 = grp%u(ind1,3)
   rw1 = grp%u(ind1,4)
   rh1 = grp%u(ind1,5) + pr1

   u1  = ru1*rinv1
   v1  = rv1*rinv1
   w1  = rw1*rinv1
   normalVelocity1 = u1*rnx+v1*rny+w1*rnz

   f11 = rho1*normalVelocity1
   f21 = pr1*rnx+rho1*u1*normalVelocity1
   f31 = pr1*rny+rho1*v1*normalVelocity1
   f41 = pr1*rnz+rho1*w1*normalVelocity1
   f51 = rh1*normalVelocity1

   pr2 = grp%p(ind2) ! pressure
   rho2 = grp%u(ind2,1)
   rinv2  = 1./rho2
   ru2 = grp%u(ind2,2)
   rv2 = grp%u(ind2,3)
   rw2 = grp%u(ind2,4)
   rh2 = grp%u(ind2,5) + pr2

   u2  = ru2*rinv2
   v2  = rv2*rinv2
   w2  = rw2*rinv2
   normalVelocity2 = u2*rnx+v2*rny+w2*rnz

   f12 = rho2*normalVelocity2
   f22 = pr2*rnx+rho2*u2*normalVelocity2
   f32 = pr2*rny+rho2*v2*normalVelocity2
   f42 = pr2*rnz+rho2*w2*normalVelocity2
   f52 = rh2*normalVelocity2

   if(ivd%boundaryTerm==1) then 
    grp%rhs(ind1,1) = grp%rhs(ind1,1)+rnorm*(3.*f11+f12)
    grp%rhs(ind1,2) = grp%rhs(ind1,2)+rnorm*(3.*f21+f22)
    grp%rhs(ind1,3) = grp%rhs(ind1,3)+rnorm*(3.*f31+f32)
    grp%rhs(ind1,4) = grp%rhs(ind1,4)+rnorm*(3.*f41+f42)
    grp%rhs(ind1,5) = grp%rhs(ind1,5)+rnorm*(3.*f51+f52)

    grp%rhs(ind2,1) = grp%rhs(ind2,1)+rnorm*(f11+3.*f12)
    grp%rhs(ind2,2) = grp%rhs(ind2,2)+rnorm*(f21+3.*f22)
    grp%rhs(ind2,3) = grp%rhs(ind2,3)+rnorm*(f31+3.*f32)
    grp%rhs(ind2,4) = grp%rhs(ind2,4)+rnorm*(f41+3.*f42)
    grp%rhs(ind2,5) = grp%rhs(ind2,5)+rnorm*(f51+3.*f52)
   else
    grp%rhs(ind1,1) = grp%rhs(ind1,1)+rnorm*4.*f11
    grp%rhs(ind1,2) = grp%rhs(ind1,2)+rnorm*4.*f21
    grp%rhs(ind1,3) = grp%rhs(ind1,3)+rnorm*4.*f31
    grp%rhs(ind1,4) = grp%rhs(ind1,4)+rnorm*4.*f41
    grp%rhs(ind1,5) = grp%rhs(ind1,5)+rnorm*4.*f51

    grp%rhs(ind2,1) = grp%rhs(ind2,1)+rnorm*4.*f12
    grp%rhs(ind2,2) = grp%rhs(ind2,2)+rnorm*4.*f22
    grp%rhs(ind2,3) = grp%rhs(ind2,3)+rnorm*4.*f32
    grp%rhs(ind2,4) = grp%rhs(ind2,4)+rnorm*4.*f42
    grp%rhs(ind2,5) = grp%rhs(ind2,5)+rnorm*4.*f52
   end if
  end if
 else if(faceIndicator.le.4) then ! faceIndicator >2
  wx  = grp%brp%sideWeightsArray(i,1)
  wy  = grp%brp%sideWeightsArray(i,2)
  wz  = grp%brp%sideWeightsArray(i,3)
  rnorm=sqrt(wx*wx+wy*wy+wz*wz)
  if(rnorm>1.0e-10) then
   rnx = wx/rnorm
   rny = wy/rnorm
   rnz = wz/rnorm

   ro    = ivd%inflowField(1,faceIndicator-2)
   ro1   = 1./ro
   uo    = ivd%inflowField(2,faceIndicator-2)
   vo    = ivd%inflowField(3,faceIndicator-2)
   wo    = ivd%inflowField(4,faceIndicator-2)
   eo    = ivd%inflowField(5,faceIndicator-2)
   po    = ivd%inflowField(6,faceIndicator-2)

   pr1l = grp%p(ind1) ! pressure
   rho1l = grp%u(ind1,1)
   rinv1l  = 1./rho1l
   ru1l = grp%u(ind1,2)
   rv1l = grp%u(ind1,3)
   rw1l = grp%u(ind1,4)
   rh1l = grp%u(ind1,5) + pr1l
   u1l  = ru1l*rinv1l
   v1l  = rv1l*rinv1l
   w1l  = rw1l*rinv1l
   eps1l = grp%u(ind1,5)
   h1l  = rh1l*rinv1l
   normalVelocity1l = u1l*rnx+v1l*rny+w1l*rnz

   fl11 = rho1l*normalVelocity1l
   fl21 = pr1l*rnx+rho1l*u1l*normalVelocity1l
   fl31 = pr1l*rny+rho1l*v1l*normalVelocity1l
   fl41 = pr1l*rnz+rho1l*w1l*normalVelocity1l
   fl51 = rh1l*normalVelocity1l

   rhor  = ro
   uxr   = uo
   vyr   = vo
   wzr   = wo
   epsr  = eo*ro
   presr = po
   hr    = eo + po*ro1

   unr   = uxr*rnx+vyr*rny+wzr*rnz
   fr11  = rhor*unr
   fr21  = rnx*presr+rhor*uxr*unr
   fr31  = rny*presr+rhor*vyr*unr
   fr41  = rnz*presr+rhor*wzr*unr
   fr51  = (epsr+presr)*unr

   di    = sqrt(rhor/rho1l)
   d1    = 1.0/(di+1.0)
   ui    = (di*uxr+u1l)*d1
   vi    = (di*vyr+v1l)*d1
   wi    = (di*wzr+w1l)*d1
   hi    = (di*hr+h1l)*d1
   ci2   = gam1*(hi-0.5*(ui*ui+vi*vi+wi*wi))
   ci2   = max(ci2,1.0e-5)
   ci    = sqrt(ci2)
   af    = 0.5*(ui*ui+vi*vi+wi*wi)
   ucp   = ui*rnx+vi*rny+wi*rnz

   du1   = rhor-rho1l
   du2   = rhor*uxr-rho1l*u1l
   du3   = rhor*vyr-rho1l*v1l
   du4   = rhor*wzr-rho1l*w1l
   du5   = epsr-eps1l

   rlam1 = abs(ucp+ci)
   rlam2 = abs(ucp-ci)
   rlam3 = abs(ucp)
   if(rlam1.lt.epslm) rlam1 = 0.5*(rlam1*rlam1*eps1+epslm)
   if(rlam2.lt.epslm) rlam2 = 0.5*(rlam2*rlam2*eps1+epslm)
   if(rlam3.lt.epslm) rlam3 = 0.5*(rlam3*rlam3*eps1+epslm)

   s1    = 0.5*(rlam1+rlam2)
   s2    = 0.5*(rlam1-rlam2)
   al1x  = gam1*(af*du1-ui*du2-vi*du3-wi*du4+du5)
   al2x  = -ucp*du1+du2*rnx+du3*rny+du4*rnz
   cc1   = ((s1-rlam3)*al1x/ci2)+(s2*al2x/ci)
   cc2   = (s2*al1x/ci)+(s1-rlam3)*al2x
   f11   = 0.5*(fr11+fl11)-0.5*(rlam3*du1+cc1           )
   f21   = 0.5*(fr21+fl21)-0.5*(rlam3*du2+cc1*ui+cc2*rnx)
   f31   = 0.5*(fr31+fl31)-0.5*(rlam3*du3+cc1*vi+cc2*rny)
   f41   = 0.5*(fr41+fl41)-0.5*(rlam3*du4+cc1*wi+cc2*rnz)
   f51   = 0.5*(fr51+fl51)-0.5*(rlam3*du5+cc1*hi+cc2*ucp)


   pr2l = grp%p(ind2) ! pressure
   rho2l = grp%u(ind2,1)
   rinv2l  = 1./rho2l
   ru2l = grp%u(ind2,2)
   rv2l = grp%u(ind2,3)
   rw2l = grp%u(ind2,4)
   rh2l = grp%u(ind2,5) + pr2l
   u2l = ru2l*rinv2l
   v2l  = rv2l*rinv2l
   w2l  = rw2l*rinv2l
   h2l  = rh2l*rinv2l
   eps2l = grp%u(ind2,5)
   normalVelocity2l = u2l*rnx+v2l*rny+w2l*rnz

   fl12 = rho2l*normalVelocity2l
   fl22 = pr2l*rnx+rho2l*u2l*normalVelocity2l
   fl32 = pr2l*rny+rho2l*v2l*normalVelocity2l
   fl42 = pr2l*rnz+rho2l*w2l*normalVelocity2l
   fl52 = rh2l*normalVelocity2l

   rhor  = ro
   uxr   = uo
   vyr   = vo
   wzr   = wo
   epsr  = eo*ro
   presr = po
   hr    = eo + po*ro1

   unr   = uxr*rnx+vyr*rny+wzr*rnz
   fr12  = rhor*unr
   fr22  = rnx*presr+rhor*uxr*unr
   fr32  = rny*presr+rhor*vyr*unr
   fr42  = rnz*presr+rhor*wzr*unr
   fr52  = (epsr+presr)*unr

   di    = sqrt(rhor/rho2l)
   d1    = 1.0/(di+1.0)
   ui    = (di*uxr+u2l)*d1
   vi    = (di*vyr+v2l)*d1
   wi    = (di*wzr+w2l)*d1
   hi    = (di*hr+h2l)*d1
   ci2   = gam1*(hi-0.5*(ui*ui+vi*vi+wi*wi))
   ci2   = max(ci2,1.0e-5)
   ci    = sqrt(ci2)
   af    = 0.5*(ui*ui+vi*vi+wi*wi)
   ucp   = ui*rnx+vi*rny+wi*rnz

   du1   = rhor-rho2l
   du2   = rhor*uxr-rho2l*u2l
   du3   = rhor*vyr-rho2l*v2l
   du4   = rhor*wzr-rho2l*w2l
   du5   = epsr-eps2l

   rlam1 = abs(ucp+ci)
   rlam2 = abs(ucp-ci)
   rlam3 = abs(ucp)
   if(rlam1.lt.epslm) rlam1 = 0.5*(rlam1*rlam1*eps1+epslm)
   if(rlam2.lt.epslm) rlam2 = 0.5*(rlam2*rlam2*eps1+epslm)
   if(rlam3.lt.epslm) rlam3 = 0.5*(rlam3*rlam3*eps1+epslm)

   s1    = 0.5*(rlam1+rlam2)
   s2    = 0.5*(rlam1-rlam2)
   al1x  = gam1*(af*du1-ui*du2-vi*du3-wi*du4+du5)
   al2x  = -ucp*du1+du2*rnx+du3*rny+du4*rnz
   cc1   = ((s1-rlam3)*al1x/ci2)+(s2*al2x/ci)
   cc2   = (s2*al1x/ci)+(s1-rlam3)*al2x
   f12   = 0.5*(fr12+fl12)-0.5*(rlam3*du1+cc1           )
   f22   = 0.5*(fr22+fl22)-0.5*(rlam3*du2+cc1*ui+cc2*rnx)
   f32   = 0.5*(fr32+fl32)-0.5*(rlam3*du3+cc1*vi+cc2*rny)
   f42   = 0.5*(fr42+fl42)-0.5*(rlam3*du4+cc1*wi+cc2*rnz)
   f52   = 0.5*(fr52+fl52)-0.5*(rlam3*du5+cc1*hi+cc2*ucp)

   if(ivd%boundaryTerm==1) then 
    grp%rhs(ind1,1) = grp%rhs(ind1,1)+rnorm*(3.*f11+f12)
    grp%rhs(ind1,2) = grp%rhs(ind1,2)+rnorm*(3.*f21+f22)
    grp%rhs(ind1,3) = grp%rhs(ind1,3)+rnorm*(3.*f31+f32)
    grp%rhs(ind1,4) = grp%rhs(ind1,4)+rnorm*(3.*f41+f42)
    grp%rhs(ind1,5) = grp%rhs(ind1,5)+rnorm*(3.*f51+f52)

    grp%rhs(ind2,1) = grp%rhs(ind2,1)+rnorm*(f11+3.*f12)
    grp%rhs(ind2,2) = grp%rhs(ind2,2)+rnorm*(f21+3.*f22) 
    grp%rhs(ind2,3) = grp%rhs(ind2,3)+rnorm*(f31+3.*f32)
    grp%rhs(ind2,4) = grp%rhs(ind2,4)+rnorm*(f41+3.*f42)
    grp%rhs(ind2,5) = grp%rhs(ind2,5)+rnorm*(f51+3.*f52)
   else
    grp%rhs(ind1,1) = grp%rhs(ind1,1)+rnorm*4.*f11
    grp%rhs(ind1,2) = grp%rhs(ind1,2)+rnorm*4.*f21
    grp%rhs(ind1,3) = grp%rhs(ind1,3)+rnorm*4.*f31
    grp%rhs(ind1,4) = grp%rhs(ind1,4)+rnorm*4.*f41
    grp%rhs(ind1,5) = grp%rhs(ind1,5)+rnorm*4.*f51

    grp%rhs(ind2,1) = grp%rhs(ind2,1)+rnorm*4.*f12
    grp%rhs(ind2,2) = grp%rhs(ind2,2)+rnorm*4.*f22
    grp%rhs(ind2,3) = grp%rhs(ind2,3)+rnorm*4.*f32
    grp%rhs(ind2,4) = grp%rhs(ind2,4)+rnorm*4.*f42
    grp%rhs(ind2,5) = grp%rhs(ind2,5)+rnorm*4.*f52
   end if
  end if
 else
  pr1 = grp%p(ind1)             ! pressure
  r1 = grp%u(ind1,1)
  rinv1 = 1./grp%u(ind1,1)      ! rho inverse
  ru1 = grp%u(ind1,2)           ! x-momentum
  rv1 = grp%u(ind1,3)           ! y-momentum
  rw1 = grp%u(ind1,4)           ! z-momentum
  u1 = rinv1*ru1
  v1 = rinv1*rv1
  w1 = rinv1*rw1
  E1 = grp%u(ind1,5)

  pr2 = grp%p(ind2)
  r2 = grp%u(ind2,1)
  rinv2 = 1./grp%u(ind2,1)
  ru2 = grp%u(ind2,2)
  rv2 = grp%u(ind2,3)
  rw2 = grp%u(ind2,4)
  u2 = rinv2*ru2
  v2 = rinv2*rv2
  w2 = rinv2*rw2
  E2 = grp%u(ind2,5)

  pu1 = pr1*rinv1*ru1
  pv1 = pr1*rinv1*rv1
  pw1 = pr1*rinv1*rw1
  pu2 = pr2*rinv2*ru2
  pv2 = pr2*rinv2*rv2
  pw2 = pr2*rinv2*rw2


  wx  = grp%brp%sideWeightsArray(i,1)
  wy  = grp%brp%sideWeightsArray(i,2)
  wz  = grp%brp%sideWeightsArray(i,3)

  f11 = wx*r1*u1+wy*r1*v1+wz*r1*w1
  f21 = wx*(ru1*u1+pr1) + wy*ru1*v1 + wz*ru1*w1
  f31 = wx*rv1*u1 + wy*(rv1*v1+pr1) + wz*rv1*w1
  f41 = wx*rw1*u1 + wy*rw1*v1 + wz*(rw1*w1+pr1)
  f51 = wx*(E1*u1+pu1) + wy*(E1*v1+pv1) + wz*(E1*w1+pw1)

  f12 = wx*r2*u2+wy*r2*v2+wz*r2*w2
  f22 = wx*(ru2*u2+pr2) + wy*ru2*v2 + wz*ru2*w2
  f32 = wx*rv2*u2 + wy*(rv2*v2+pr2) + wz*rv2*w2
  f42 = wx*rw2*u2 + wy*rw2*v2 + wz*(rw2*w2+pr2)
  f52 = wx*(E2*u2+pu2) + wy*(E2*v2+pv2) + wz*(E2*w2+pw2)

  f1 = f11 + f12
  f2 = f21 + f22
  f3 = f31 + f32
  f4 = f41 + f42
  f5 = f51 + f52

! assemble right hand side
  grp%rhs(ind1,1) = grp%rhs(ind1,1) + f1+2.0*f11
  grp%rhs(ind1,2) = grp%rhs(ind1,2) + f2+2.0*f21
  grp%rhs(ind1,3) = grp%rhs(ind1,3) + f3+2.0*f31
  grp%rhs(ind1,4) = grp%rhs(ind1,4) + f4+2.0*f41
  grp%rhs(ind1,5) = grp%rhs(ind1,5) + f5+2.0*f51
  grp%rhs(ind2,1) = grp%rhs(ind2,1) + f1+2.0*f12
  grp%rhs(ind2,2) = grp%rhs(ind2,2) + f2+2.0*f22
  grp%rhs(ind2,3) = grp%rhs(ind2,3) + f3+2.0*f32
  grp%rhs(ind2,4) = grp%rhs(ind2,4) + f4+2.0*f42
  grp%rhs(ind2,5) = grp%rhs(ind2,5) + f5+2.0*f52
 end if
end do

#ifdef PARALLEL

! send domain boundary residuals

    if (kk.eq.1) then
      do i=1,grp%pdp%numberOfProcesses
        if (i.ne.grp%pdp%currentDomain.and.grp%pdp%numberOfSendNodes(i)>0) then
          call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
          position = 1
          do j=1,5
            call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                               grp%pdp%realType,grp%rhs(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
          end do
          call mp_sendv(grp%pdp%processorIDs,i,27,grp%pdp%sendFlag)
        end if
      end do
    end if
  end do

  call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,27,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )


  do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
        do j=1,5
          call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                                 grp%pdp%realType,grp%rhs(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
        end do
    end if
  end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )


  do i=1,grp%pdp%numberOfProcesses
    if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i)>0)) then
      call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
      position = 1
      do j=1,5
        call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                          grp%pdp%realType,grp%rhs(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
      end do
      call mp_sendv(grp%pdp%processorIDs,i,28,grp%pdp%sendFlag)
    end if
  end do

  call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,28,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )


  do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
        do j = 1,5
          call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                             grp%pdp%realType,grp%rhs(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
        end do
    end if
  end do

#endif PARALLEL

end subroutine makeRHS 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!HLLC FLUX ROUTINES
subroutine LSComputeGradient(grp,ivd,u)
!===================================================
! Compute the gradient of the primitive variables
! to be used in HLLC hight order solver. 
! this is a Least Square method L. Remaki
!===================================================
 
 IMPLICIT NONE
 include 'mpif.h'
 
 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd
 real,pointer :: u(:,:),v(:),grad(:,:)
 integer :: i,j,i1,i2,kk,istart,ifinish
 real :: dx, dy,dz,dr2,du2,dv2,dw2,de2,det,norm,norm1,norm2
 real :: a(5)

integer ::ip1,ip2,faceIndicator,ib,ip,ist,ien,ind
real :: r1,r2,u1,u2,v1,v2,e1,e2,wx,wy,wz,drx,dux,dvx,dex,dry,duy,dvy,dey,oneOverReynoldsNumber
real :: dwx,dwy,dwz,dwdx,dwdy,dwdz,dTdz,tz,t61,t71,t81,t91,t62,t72,t82,t92,drz,duz,dvz,dez
real :: w1,w2,dudz,dvdz,u31,u32,dt6,dt7,dt8,dt9,dut3,f5,dut31,dut32
real :: mu,k,T,T0,f2,f3,f4,rx,ry,vt,t11,t21,t31,t41,t51,t12,t22,t32
real :: t42,t52,dudx,dvdx,dTdx,dudy,dvdy,dTdy,dt1,dt2,dt3,dt4,dt5,oneOverPrandtlNumber
real :: twoThirds,fourThirds,gammaOverGammaMinusOne,u11,u12,u21,u22,dut1,dut2
real :: wallTangent(3),wallNormal(3),tempConv,dut11,dut12,dut21,dut22,inflowTemp
real :: adTemp1,adTemp2,turbDiffusionCoefficient
real :: viscosityCoefficient,diffusionCoefficient,oneOverTurbPrandtlNumber
real :: uu1,uu2,duux,duuy,duuz,locu,dt,increment,sign, eps,fact,sum
double precision :: turbViscosityCoefficient,Ksi

 integer :: mind, allocateStatus
 real :: mres

allocate(v(grp%numberOfNodes),stat=allocateStatus)
    if (allocateStatus/=0) STOP "ERROR: Not enough memory in computehessianhat"
    
allocate(grad(grp%numberOfNodes,15),stat=allocateStatus)
    if (allocateStatus/=0) STOP "ERROR: Not enough memory in computehessianhat" 
    
    
grp%nodeHelpArray(1:grp%numberOfNodes,1:15) = 0.0
do i=1,grp%numberOfNodes
v(i) = u(i,1)
end do
 call computegradhat(grp,ivd,v)
do i=1,grp%numberOfNodes
grad(i,1) = grp%nodeHelpArray(i,10)
grad(i,6) = grp%nodeHelpArray(i,11)
grad(i,11) = grp%nodeHelpArray(i,12)
end do

do i=1,grp%numberOfNodes
v(i) = u(i,2)/u(i,1)
end do
 call computegradhat(grp,ivd,v)
do i=1,grp%numberOfNodes
grad(i,2) = grp%nodeHelpArray(i,10)
grad(i,7) = grp%nodeHelpArray(i,11)
grad(i,12) = grp%nodeHelpArray(i,12)
end do

do i=1,grp%numberOfNodes
v(i) = u(i,3)/u(i,1)
end do
 call computegradhat(grp,ivd,v)
do i=1,grp%numberOfNodes
grad(i,3) = grp%nodeHelpArray(i,10)
grad(i,8) = grp%nodeHelpArray(i,11)
grad(i,13) = grp%nodeHelpArray(i,12)
end do

do i=1,grp%numberOfNodes
v(i) = u(i,4)/u(i,1)
end do
 call computegradhat(grp,ivd,v)
 do i=1,grp%numberOfNodes
grad(i,4) = grp%nodeHelpArray(i,10)
grad(i,9) = grp%nodeHelpArray(i,11)
grad(i,14) = grp%nodeHelpArray(i,12)
end do
 
 
do i=1,grp%numberOfNodes
v(i) = u(i,5)
end do
 call computegradhat(grp,ivd,v) 
 do i=1,grp%numberOfNodes
grad(i,5) = grp%nodeHelpArray(i,10)
grad(i,10) = grp%nodeHelpArray(i,11)
grad(i,15) = grp%nodeHelpArray(i,12)
end do
do i=1,grp%numberOfNodes
 grp%nodeHelpArray(i,1) = grad(i,1) 
 grp%nodeHelpArray(i,2) = grad(i,2) 
 grp%nodeHelpArray(i,3) = grad(i,3) 
 grp%nodeHelpArray(i,4) = grad(i,4) 
 grp%nodeHelpArray(i,5) = grad(i,5) 
 
 grp%nodeHelpArray(i,6) = grad(i,6) 
 grp%nodeHelpArray(i,7) = grad(i,7) 
 grp%nodeHelpArray(i,8) = grad(i,8) 
 grp%nodeHelpArray(i,9) = grad(i,9) 
 grp%nodeHelpArray(i,10) = grad(i,10) 
 
 grp%nodeHelpArray(i,11) = grad(i,11) 
 grp%nodeHelpArray(i,12) = grad(i,12) 
 grp%nodeHelpArray(i,13) = grad(i,13) 
 grp%nodeHelpArray(i,14) = grad(i,14) 
 grp%nodeHelpArray(i,15) = grad(i,15) 

 end do
 
 deallocate(v,stat=allocateStatus)
    if (allocateStatus/=0) STOP "ERROR: Could not dellocate computehessianhat"
    
 deallocate(grad,stat=allocateStatus)
    if (allocateStatus/=0) STOP "ERROR: Could not dellocate computehessianhat"
 end subroutine LSComputeGradient
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Compute a LS reconstructed gradient   L. Remaki
!-----------------------------------------------------------------------
 subroutine computegradhat(grp,ivd,u)
 
  IMPLICIT NONE
 include 'mpif.h'

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd
 real,pointer :: u(:)

 integer :: i,j,i1,i2,kk,istart,ifinish, ismoothgrad, ismoo
 real :: dx, dy,dz,dr2,du2,dv2,dw2,de2,det,dist,dist1,dist2,wx,wy,wz,rnorm,gradl,gradr,sum
 real :: a(3,3), c(3,3), b(3)
 
 grp%nodeHelpArray(1:grp%numberOfNodes,1:9) = 0.0
 
#ifdef PARALLEL

 do kk=1,2
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfComSides
 else
  istart = grp%pdp%numberOfComSides + 1
  ifinish = grp%numberOfSides
 end if


 do i=istart,ifinish
 
#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL


 i1 = grp%sideIndexArray(i,1)
 i2 = grp%sideIndexArray(i,2)
 
 dx = grp%coordinates(i2,1) - grp%coordinates(i1,1)
 dy = grp%coordinates(i2,2) - grp%coordinates(i1,2)
 dz = grp%coordinates(i2,3) - grp%coordinates(i1,3)
 
 dist1 = 1./((dx*dx +dy*dy +dz*dz))
 dist2 = dist1
 dr2 = (u(i2) - u(i1))

 grp%nodeHelpArray(i1,1) =  grp%nodeHelpArray(i1,1) + dx*dx*dist1
 grp%nodeHelpArray(i1,2) =  grp%nodeHelpArray(i1,2) + dx*dy*dist1
 grp%nodeHelpArray(i1,3) =  grp%nodeHelpArray(i1,3) + dx*dz*dist1
 grp%nodeHelpArray(i1,4) =  grp%nodeHelpArray(i1,4) + dy*dy*dist1
 grp%nodeHelpArray(i1,5) =  grp%nodeHelpArray(i1,5) + dy*dz*dist1
 grp%nodeHelpArray(i1,6) =  grp%nodeHelpArray(i1,6) + dz*dz*dist1
 
 grp%nodeHelpArray(i2,1) =  grp%nodeHelpArray(i2,1) + dx*dx*dist2
 grp%nodeHelpArray(i2,2) =  grp%nodeHelpArray(i2,2) + dx*dy*dist2
 grp%nodeHelpArray(i2,3) =  grp%nodeHelpArray(i2,3) + dx*dz*dist2
 grp%nodeHelpArray(i2,4) =  grp%nodeHelpArray(i2,4) + dy*dy*dist2
 grp%nodeHelpArray(i2,5) =  grp%nodeHelpArray(i2,5) + dy*dz*dist2
 grp%nodeHelpArray(i2,6) =  grp%nodeHelpArray(i2,6) + dz*dz*dist2
 
 grp%nodeHelpArray(i1,7) =  grp%nodeHelpArray(i1,7) + dr2*dx*dist1
 grp%nodeHelpArray(i1,8) =  grp%nodeHelpArray(i1,8) + dr2*dy*dist1
 grp%nodeHelpArray(i1,9) =  grp%nodeHelpArray(i1,9) + dr2*dz*dist1
 
 grp%nodeHelpArray(i2,7) =  grp%nodeHelpArray(i2,7) + dr2*dx*dist2
 grp%nodeHelpArray(i2,8) =  grp%nodeHelpArray(i2,8) + dr2*dy*dist2
 grp%nodeHelpArray(i2,9) =  grp%nodeHelpArray(i2,9) + dr2*dz*dist2
  
end do

#ifdef PARALLEL

!    write(get_unit(),*)  'Starting first ArtDis comms',kk
 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses
   if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
    do j=1,9
     call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
    call mp_sendv(grp%pdp%processorIDs,i,108,grp%pdp%sendFlag)
   end if
  end do
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,108,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1239 continue

 do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
!  grp%pdp%sendFlag = 0
!  call mp_check_mesg_arrived(grp%pdp%processorIDs,i,159,grp%pdp%sendFlag,bufferzz,bufferSize)
!  if (grp%pdp%sendFlag.gt.0) then
!   grp%pdp%receivedMessageFromProcess(i) = 1
!   position = 1
    do j=1,9
     call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                            grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
!  end if
  end if
 end do

!do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1239
!end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

!  write(get_unit(),*)  'Starting second ArtDis comms'
 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position = 1
   do j=1,9
    call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                      grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
   end do
   call mp_sendv(grp%pdp%processorIDs,i,109,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,109,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1339 continue

 do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
!  grp%pdp%sendFlag = 0
!  call mp_check_mesg_arrived(grp%pdp%processorIDs,i,60,grp%pdp%sendFlag,bufferzz,bufferSize)
!  if (grp%pdp%sendFlag.gt.0) then
!   grp%pdp%receivedMessageFromProcess(i) = 1
!   position = 1
    do j = 1,9
     call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
!  end if
  end if
 end do

!do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1339
!end do
#endif PARALLEL

  do i=1,grp%numberOfNodes
   a(1,1) = grp%nodeHelpArray(i,1)
   a(1,2) = grp%nodeHelpArray(i,2)
   a(1,3) = grp%nodeHelpArray(i,3)
   a(2,1) = grp%nodeHelpArray(i,2)
   a(2,2) = grp%nodeHelpArray(i,4)
   a(2,3) = grp%nodeHelpArray(i,5)
   a(3,1) = grp%nodeHelpArray(i,3)
   a(3,2) = grp%nodeHelpArray(i,5)
   a(3,3) = grp%nodeHelpArray(i,6)
   c(1,1) = a(2,2)*a(3,3) - a(2,3)*a(3,2)
   c(2,1) = a(2,3)*a(3,1) - a(2,1)*a(3,3)
   c(3,1) = a(2,1)*a(3,2) - a(2,2)*a(3,1)
   det  = a(1,1)*c(1,1) + a(1,2)*c(2,1) + a(1,3)*c(3,1)
   
   !write(*,*) "...........det              ", det
   if(abs(det) .gt. 0.0) det = 1.0/det
   
   c(1,1) =  det*c(1,1)
   c(2,1) =  det*c(2,1)
   c(3,1) =  det*c(3,1)
   c(1,2) =  det*( a(1,3) * a(3,2) - a(1,2) * a(3,3) )
   c(2,2) =  det*( a(1,1) * a(3,3) - a(1,3) * a(3,1) )
   c(3,2) =  det*( a(1,2) * a(3,1) - a(1,1) * a(3,2) )
   c(1,3) =  det*( a(1,2) * a(2,3) - a(1,3) * a(2,2) )
   c(2,3) =  det*( a(1,3) * a(2,1) - a(2,3) * a(1,1) )
   c(3,3) =  det*( a(1,1) * a(2,2) - a(1,2) * a(2,1) )
    
   b(1) = grp%nodeHelpArray(i,7)
   b(2) = grp%nodeHelpArray(i,8)
   b(3) = grp%nodeHelpArray(i,9)   
   grp%nodeHelpArray(i,10) = c(1,1)*b(1) + c(1,2)*b(2) + c(1,3)*b(3)
   grp%nodeHelpArray(i,11) = c(2,1)*b(1) + c(2,2)*b(2) + c(2,3)*b(3)
   grp%nodeHelpArray(i,12) = c(3,1)*b(1) + c(3,2)*b(2) + c(3,3)*b(3)
  end do

 end subroutine computegradhat
!===================================================

subroutine computeGradient(grp,ivd,u)
!===================================================
! Compute the gradient of the primitive variables
! to be used in HLLC hight order solver. 
! this is a Galerkin method . L. Remaki
!===================================================
 
 IMPLICIT NONE
 include 'mpif.h'
 
 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd
 real,pointer :: u(:,:)
 integer :: i,j,i1,i2,kk,istart,ifinish

real :: r1,r2,u1,u2,v1,v2,w1,w2,e1,e2,wx,wy,wz
real :: drx,dux,dvx,dwx,dex,dry,duy,dvy,dwy,dey,drz,duz,dvz,dwz,dez,vt

grp%nodeHelpArray(1:grp%numberOfNodes,1:15) = 0.0

#ifdef PARALLEL

do kk=1,2
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfComSides
 else
  istart = grp%pdp%numberOfComSides + 1
  ifinish = grp%numberOfSides
 end if


 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL
  

! indexes of nodes in side
 i1 = grp%sideIndexArray(i,1)
 i2 = grp%sideIndexArray(i,2)

! side weights from preprocessor 
 wx = grp%sideWeightsArray(i,1)  
 wy = grp%sideWeightsArray(i,2)  
 wz = grp%sideWeightsArray(i,3)  

  r1   = u(i1,1) 
  u1   = u(i1,2)/r1 
  v1   = u(i1,3)/r1 
  w1   = u(i1,4)/r1 
  e1   = u(i1,5)
  r2   = u(i2,1)
  u2   = u(i2,2)/r2
  v2   = u(i2,3)/r2
  w2   = u(i2,4)/r2
  e2   = u(i2,5)

 drx   = (r1+r2)*wx
 dux   = (u1+u2)*wx
 dvx   = (v1+v2)*wx
 dwx   = (w1+w2)*wx
 dex   = (e1+e2)*wx
 dry   = (r1+r2)*wy
 duy   = (u1+u2)*wy
 dvy   = (v1+v2)*wy
 dwy   = (w1+w2)*wy
 dey   = (e1+e2)*wy
 drz   = (r1+r2)*wz
 duz   = (u1+u2)*wz
 dvz   = (v1+v2)*wz
 dwz   = (w1+w2)*wz
 dez   = (e1+e2)*wz

 grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + drx
 grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + dux
 grp%nodeHelpArray(i1,3) = grp%nodeHelpArray(i1,3) + dvx
 grp%nodeHelpArray(i1,4) = grp%nodeHelpArray(i1,4) + dwx
 grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) + dex
 
 grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) + dry
 grp%nodeHelpArray(i1,7) = grp%nodeHelpArray(i1,7) + duy
 grp%nodeHelpArray(i1,8) = grp%nodeHelpArray(i1,8) + dvy
 grp%nodeHelpArray(i1,9) = grp%nodeHelpArray(i1,9) + dwy
 grp%nodeHelpArray(i1,10) = grp%nodeHelpArray(i1,10) + dey
 
 grp%nodeHelpArray(i1,11) = grp%nodeHelpArray(i1,11) + drz
 grp%nodeHelpArray(i1,12) = grp%nodeHelpArray(i1,12) + duz
 grp%nodeHelpArray(i1,13) = grp%nodeHelpArray(i1,13) + dvz
 grp%nodeHelpArray(i1,14) = grp%nodeHelpArray(i1,14) + dwz
 grp%nodeHelpArray(i1,15) = grp%nodeHelpArray(i1,15) + dez

 grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) - drx
 grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) - dux
 grp%nodeHelpArray(i2,3) = grp%nodeHelpArray(i2,3) - dvx
 grp%nodeHelpArray(i2,4) = grp%nodeHelpArray(i2,4) - dwx
 grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) - dex

 grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) - dry
 grp%nodeHelpArray(i2,7) = grp%nodeHelpArray(i2,7) - duy
 grp%nodeHelpArray(i2,8) = grp%nodeHelpArray(i2,8) - dvy
 grp%nodeHelpArray(i2,9) = grp%nodeHelpArray(i2,9) - dwy
 grp%nodeHelpArray(i2,10) = grp%nodeHelpArray(i2,10) - dey

 grp%nodeHelpArray(i2,11) = grp%nodeHelpArray(i2,11) - drz
 grp%nodeHelpArray(i2,12) = grp%nodeHelpArray(i2,12) - duz
 grp%nodeHelpArray(i2,13) = grp%nodeHelpArray(i2,13) - dvz
 grp%nodeHelpArray(i2,14) = grp%nodeHelpArray(i2,14) - dwz
 grp%nodeHelpArray(i2,15) = grp%nodeHelpArray(i2,15) - dez

end do


#ifdef PARALLEL
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfBoundaryComSides
 else
  istart = grp%pdp%numberOfBoundaryComSides + 1
  ifinish = grp%brp%numberOfBoundarySides
 end if
 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL

 i1 = grp%brp%sideIndexArray(i,1)
 i2 = grp%brp%sideIndexArray(i,2)
 wx   = grp%brp%sideWeightsArray(i,1)
 wy   = grp%brp%sideWeightsArray(i,2) 
 wz   = grp%brp%sideWeightsArray(i,3) 
  r1   = u(i1,1) 
  u1   = u(i1,2)/r1 
  v1   = u(i1,3)/r1 
  w1   = u(i1,4)/r1 
  e1   = u(i1,5)
  r2   = u(i2,1)
  u2   = u(i2,2)/r2
  v2   = u(i2,3)/r2
  w2   = u(i2,4)/r2
  e2   = u(i2,5)
 drx   = r1+r2
 dux   = u1+u2
 dvx   = v1+v2
 dwx   = w1+w2
 dex   = e1+e2
 dry   = r1+r2
 duy   = u1+u2
 dvy   = v1+v2
 dwy   = w1+w2
 dey   = e1+e2
 drz   = r1+r2
 duz   = u1+u2
 dvz   = v1+v2
 dwz   = w1+w2
 dez   = e1+e2

! adding boundary face contributions        

 grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + (2.*r1+drx)*wx
 grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + (2.*u1+dux)*wx
 grp%nodeHelpArray(i1,3) = grp%nodeHelpArray(i1,3) + (2.*v1+dvx)*wx
 grp%nodeHelpArray(i1,4) = grp%nodeHelpArray(i1,4) + (2.*w1+dwx)*wx
 grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) + (2.*e1+dex)*wx
 
 grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) + (2.*r1+dry)*wy
 grp%nodeHelpArray(i1,7) = grp%nodeHelpArray(i1,7) + (2.*u1+duy)*wy
 grp%nodeHelpArray(i1,8) = grp%nodeHelpArray(i1,8) + (2.*v1+dvy)*wy
 grp%nodeHelpArray(i1,9) = grp%nodeHelpArray(i1,9) + (2.*w1+dwy)*wy
 grp%nodeHelpArray(i1,10) = grp%nodeHelpArray(i1,10) + (2.*e1+dey)*wy
 
 grp%nodeHelpArray(i1,11) = grp%nodeHelpArray(i1,11) + (2.*r1+drz)*wz
 grp%nodeHelpArray(i1,12) = grp%nodeHelpArray(i1,12) + (2.*u1+duz)*wz
 grp%nodeHelpArray(i1,13) = grp%nodeHelpArray(i1,13) + (2.*v1+dvz)*wz
 grp%nodeHelpArray(i1,14) = grp%nodeHelpArray(i1,14) + (2.*w1+dwz)*wz
 grp%nodeHelpArray(i1,15) = grp%nodeHelpArray(i1,15) + (2.*e1+dez)*wz
        
 
 grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + (2.*r2+drx)*wx
 grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) + (2.*u2+dux)*wx
 grp%nodeHelpArray(i2,3) = grp%nodeHelpArray(i2,3) + (2.*v2+dvx)*wx
 grp%nodeHelpArray(i2,4) = grp%nodeHelpArray(i2,4) + (2.*w2+dwx)*wx
 grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) + (2.*e2+dex)*wx

 grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) + (2.*r2+dry)*wy
 grp%nodeHelpArray(i2,7) = grp%nodeHelpArray(i2,7) + (2.*u2+duy)*wy
 grp%nodeHelpArray(i2,8) = grp%nodeHelpArray(i2,8) + (2.*v2+dvy)*wy
 grp%nodeHelpArray(i2,9) = grp%nodeHelpArray(i2,9) + (2.*w2+dwy)*wy
 grp%nodeHelpArray(i2,10) = grp%nodeHelpArray(i2,10) + (2.*e2+dey)*wy

 grp%nodeHelpArray(i2,11) = grp%nodeHelpArray(i2,11) + (2.*r2+drz)*wz
 grp%nodeHelpArray(i2,12) = grp%nodeHelpArray(i2,12) + (2.*u2+duz)*wz
 grp%nodeHelpArray(i2,13) = grp%nodeHelpArray(i2,13) + (2.*v2+dvz)*wz
 grp%nodeHelpArray(i2,14) = grp%nodeHelpArray(i2,14) + (2.*w2+dwz)*wz
 grp%nodeHelpArray(i2,15) = grp%nodeHelpArray(i2,15) + (2.*e2+dez)*wz
 
end do

#ifdef PARALLEL

!    write(get_unit(),*)  'Starting first ArtDis comms',kk
 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses
   if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
    do j=1,15
     call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
    call mp_sendv(grp%pdp%processorIDs,i,59,grp%pdp%sendFlag)
   end if
  end do
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,59,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1239 continue

 do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
!  grp%pdp%sendFlag = 0
!  call mp_check_mesg_arrived(grp%pdp%processorIDs,i,59,grp%pdp%sendFlag,bufferzz,bufferSize)
!  if (grp%pdp%sendFlag.gt.0) then
!   grp%pdp%receivedMessageFromProcess(i) = 1
!   position = 1
    do j=1,15
     call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                            grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
!  end if
  end if
 end do

!do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1239
!end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

!  write(get_unit(),*)  'Starting second ArtDis comms'
 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position = 1
   do j=1,15
    call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                      grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
   end do
   call mp_sendv(grp%pdp%processorIDs,i,60,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,60,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1339 continue

 do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
!  grp%pdp%sendFlag = 0
!  call mp_check_mesg_arrived(grp%pdp%processorIDs,i,60,grp%pdp%sendFlag,bufferzz,bufferSize)
!  if (grp%pdp%sendFlag.gt.0) then
!   grp%pdp%receivedMessageFromProcess(i) = 1
!   position = 1
    do j = 1,15
     call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
!  end if
  end if
 end do

!do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1339
!end do
#endif PARALLEL

! divide by volume to remove mass matrix on LHS
do i=1,grp%numberOfNodes
 vt = 1./grp%nodeVolume(i)
 grp%nodeHelpArray(i,1)   = grp%nodeHelpArray(i,1)*vt
 grp%nodeHelpArray(i,2)   = grp%nodeHelpArray(i,2)*vt
 grp%nodeHelpArray(i,3)   = grp%nodeHelpArray(i,3)*vt
 grp%nodeHelpArray(i,4)   = grp%nodeHelpArray(i,4)*vt
 grp%nodeHelpArray(i,5)   = grp%nodeHelpArray(i,5)*vt
 grp%nodeHelpArray(i,6)   = grp%nodeHelpArray(i,6)*vt
 grp%nodeHelpArray(i,7)   = grp%nodeHelpArray(i,7)*vt
 grp%nodeHelpArray(i,8)   = grp%nodeHelpArray(i,8)*vt
 grp%nodeHelpArray(i,9)   = grp%nodeHelpArray(i,9)*vt
 grp%nodeHelpArray(i,10)  = grp%nodeHelpArray(i,10)*vt
 grp%nodeHelpArray(i,11)  = grp%nodeHelpArray(i,11)*vt
 grp%nodeHelpArray(i,12)  = grp%nodeHelpArray(i,12)*vt
 grp%nodeHelpArray(i,13)  = grp%nodeHelpArray(i,13)*vt
 grp%nodeHelpArray(i,14)  = grp%nodeHelpArray(i,14)*vt
 grp%nodeHelpArray(i,15)  = grp%nodeHelpArray(i,15)*vt
 end do

 end subroutine computeGradient
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine LocalMinMax(grp,ivd,u)
!===================================================
! Compute the local Max and Min of varaibles to use 
! for the authomatic local clipping. Remaki
!===================================================
 
 IMPLICIT NONE
 include 'mpif.h'
 
 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd
 real,pointer :: u(:,:)
 integer :: i,j,i0,i1,i2,kk,istart,ifinish,ind1,ind2

real :: r1,r2,u1,u2,v1,v2,w1,w2,e1,e2
real :: dr1,du1,dv1,dw1,de1,dr2,du2,dv2,dw2,de2

 grp%nodeHelpArray(1:grp%numberOfNodes,16:25) = 1.0e15
 
#ifdef PARALLEL

do kk=1,2
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfComSides
 else
  istart = grp%pdp%numberOfComSides + 1
  ifinish = grp%numberOfSides
 end if


 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL
  

! indexes of nodes in side
 ind1 = grp%sideIndexArray(i,1)
 ind2 = grp%sideIndexArray(i,2)

! side weights from preprocessor 
  
  r1   = grp%u(ind1,1) 
  u1   = grp%u(ind1,2)/r1
  v1   = grp%u(ind1,3)/r1
  w1   = grp%u(ind1,4)/r1
  e1   = grp%u(ind1,5)
 
  r2   = grp%u(ind2,1)
  u2   = grp%u(ind2,2)/r2
  v2   = grp%u(ind2,3)/r2
  w2   = grp%u(ind2,4)/r2
  e2   = grp%u(ind2,5)

  dr1 = amin1(r1,r2)
  du1 = amin1(u1,u2)
  dv1 = amin1(v1,v2)
  dw1 = amin1(w1,w2)
  de1 = amin1(e1,e2)
    
  dr2 = amin1(-r1,-r2)
  du2 = amin1(-u1,-u2)
  dv2 = amin1(-v1,-v2)
  dw2 = amin1(-w1,-w2)
  de2 = amin1(-e1,-e2)
  
 grp%nodeHelpArray(ind1,16) = amin1(grp%nodeHelpArray(ind1,16),dr1)
 grp%nodeHelpArray(ind1,17) = amin1(grp%nodeHelpArray(ind1,17),du1)
 grp%nodeHelpArray(ind1,18) = amin1(grp%nodeHelpArray(ind1,18),dv1)
 grp%nodeHelpArray(ind1,19) = amin1(grp%nodeHelpArray(ind1,19),dw1)
 grp%nodeHelpArray(ind1,20) = amin1(grp%nodeHelpArray(ind1,20),de1)

 grp%nodeHelpArray(ind2,16) = amin1(grp%nodeHelpArray(ind2,16),dr1)
 grp%nodeHelpArray(ind2,17) = amin1(grp%nodeHelpArray(ind2,17),du1)
 grp%nodeHelpArray(ind2,18) = amin1(grp%nodeHelpArray(ind2,18),dv1)
 grp%nodeHelpArray(ind2,19) = amin1(grp%nodeHelpArray(ind2,19),dw1)
 grp%nodeHelpArray(ind2,20) = amin1(grp%nodeHelpArray(ind2,20),de1)
 
 grp%nodeHelpArray(ind1,21) = amin1(grp%nodeHelpArray(ind1,21),dr2)
 grp%nodeHelpArray(ind1,22) = amin1(grp%nodeHelpArray(ind1,22),du2)
 grp%nodeHelpArray(ind1,23) = amin1(grp%nodeHelpArray(ind1,23),dv2)
 grp%nodeHelpArray(ind1,24) = amin1(grp%nodeHelpArray(ind1,24),dw2)
 grp%nodeHelpArray(ind1,25) = amin1(grp%nodeHelpArray(ind1,25),de2)

 grp%nodeHelpArray(ind2,21) = amin1(grp%nodeHelpArray(ind2,21),dr2)
 grp%nodeHelpArray(ind2,22) = amin1(grp%nodeHelpArray(ind2,22),du2)
 grp%nodeHelpArray(ind2,23) = amin1(grp%nodeHelpArray(ind2,23),dv2)
 grp%nodeHelpArray(ind2,24) = amin1(grp%nodeHelpArray(ind2,24),dw2)
 grp%nodeHelpArray(ind2,25) = amin1(grp%nodeHelpArray(ind2,25),de2)

 
end do
#ifdef PARALLEL

!    write(get_unit(),*)  'Starting first ArtDis comms',kk
 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses
   if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
    do j=16,25
     call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
    call mp_sendv(grp%pdp%processorIDs,i,96,grp%pdp%sendFlag)
   end if
  end do  
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,96,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1239 continue

 do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
!  grp%pdp%sendFlag = 0
!  call mp_check_mesg_arrived(grp%pdp%processorIDs,i,59,grp%pdp%sendFlag,bufferzz,bufferSize)
!  if (grp%pdp%sendFlag.gt.0) then
!   grp%pdp%receivedMessageFromProcess(i) = 1
!   position = 1
    do j=16,25
     call mp_upakv_min_bdry_jj(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                            grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
!  end if
  end if
 end do

!do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1239
!end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

!  write(get_unit(),*)  'Starting second ArtDis comms'
 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position = 1
   do j=16,25
    call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                      grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
   end do
   call mp_sendv(grp%pdp%processorIDs,i,97,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,97,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1339 continue

 do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
!  grp%pdp%sendFlag = 0
!  call mp_check_mesg_arrived(grp%pdp%processorIDs,i,60,grp%pdp%sendFlag,bufferzz,bufferSize)
!  if (grp%pdp%sendFlag.gt.0) then
!   grp%pdp%receivedMessageFromProcess(i) = 1
!   position = 1
    do j = 16,25
     call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
!  end if
  end if
 end do

!do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1339
!end do
#endif PARALLEL

do i=1,grp%numberOfNodes
 

 grp%nodeHelpArray(i,16) =  grp%nodeHelpArray(i,16) -  grp%u(i,1)
 grp%nodeHelpArray(i,17) =  grp%nodeHelpArray(i,17) -  (grp%u(i,2)/grp%u(i,1))
 grp%nodeHelpArray(i,18) =  grp%nodeHelpArray(i,18) -  (grp%u(i,3)/grp%u(i,1))
 grp%nodeHelpArray(i,19) =  grp%nodeHelpArray(i,19) -  (grp%u(i,4)/grp%u(i,1))
 grp%nodeHelpArray(i,20) =  grp%nodeHelpArray(i,20) -  grp%u(i,5)

 grp%nodeHelpArray(i,21) =  -grp%nodeHelpArray(i,21) - grp%u(i,1)
 grp%nodeHelpArray(i,22) =  -grp%nodeHelpArray(i,22) - (grp%u(i,2)/grp%u(i,1))
 grp%nodeHelpArray(i,23) =  -grp%nodeHelpArray(i,23) - (grp%u(i,3)/grp%u(i,1))
 grp%nodeHelpArray(i,24) =  -grp%nodeHelpArray(i,24) - (grp%u(i,4)/grp%u(i,1))
 grp%nodeHelpArray(i,25) =  -grp%nodeHelpArray(i,25) - grp%u(i,5)

 end do

 end subroutine LocalMinMax
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
real function minmod(a,b,hlength)
 ! A minmod fucntion used as a limiter in HLLC solver 
 ! to drop the order to first order in a discontinuities vicinity
 ! L.Remaki 
 IMPLICIT NONE
 real :: a,b,c,hlength
   if( (a*b .le. 0.0) .or. (hlength .lt. 0.0) ) then
    c = 0.0
   else
    if (abs(a) .le. abs(b)) c = a
    if (abs(a) .gt. abs(b)) c = b 
    !if((abs(c) .gt. hlength) .and. (hlength .gt. 0.0)) then
    !c= hlength*c/abs(c)
    !end if
   end if
 minmod = c
 end function minmod
!-----------------------------------------------------------------------
real function minmod1(a,b,hlength,hessgrad)
 ! A minmod fucntion used as a limiter in HLLC solver 
 ! to drop the order to first order in a discontinuities vicinity
 ! L.Remaki 
 IMPLICIT NONE
 real :: a,b,c,hlength,hessgrad
    if( (a*b .le. 0.0)  ) then
    c = 0.0
   else
   if (abs(a) .le. abs(b)) c = a
   if (abs(a) .gt. abs(b)) c = b 
   if((abs(c) .gt. hlength) .and. (hlength .gt. 0.0)) then
   c= hlength*c/abs(c)
   end if
  ! if((abs(c) .gt. hessgrad)) then
  ! c= hessgrad*c/abs(c)
  ! end if
   end if
 minmod1 = c
 end function minmod1
!-----------------------------------------------------------------------
!----------------------------------------------------------------
real function ledim(a,bmin,bmax)
 ! A minmod fucntion used as a limiter in HLLC solver 
 ! to drop the order to first order in a discontinuities vicinity
 ! L.Remaki 
 IMPLICIT NONE
 real :: a,a1,bmin,bmax,dq,c,invmax,invmin
 integer :: nsig
   nsig =2   
    c = a
    a1 = 2.0*a
    if( (a1 .gt. 0.0) .and.(a1 .gt. bmax)) then
     invmax =0.0
       if(bmax .gt. 0.0) invmax =1.0/bmax
       c= 0.5*bmax*((a1*invmax))/((1.+(a1*invmax)**nsig)**(1.0/nsig))
    end if
    
    if( (a1 .lt. 0.0) .and.(a1 .lt. bmin) ) then
     invmin =0.0
       if(bmin .lt. 0.0) invmin =1.0/bmin
       c= 0.5*bmin*((a1*invmin))/((1.+(a1*invmin)**nsig)**(1.0/nsig))
    end if
    
 ledim = c
 end function ledim

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
subroutine setUpInitialField(ivd,grp,INFILE)!START-GUST
! initializes the unknown fields, either by reading from file or inserting 
! freestream values
IMPLICIT NONE

type(InputVariablesData) :: ivd
type(GridSolverData) :: grp
integer :: INFILE
!START-TURB
double precision :: xu(7)
!END-TURB
integer :: i,j,number

integer :: ind

real dist
INTEGER :: physicalTimestepNumber!START-GUST
!START-TURB
real :: betaomega
betaomega = 3.0/40.0
!END-TURB

if(INFILE>0) then 
! start up from given field
 write(*,*) "starting up from given field"

  read(INFILE) number 
  if(number.ne.grp%numberOfNodes) then 
   write(*,*) number,grp%numberOfNodes
   STOP "ERROR: Startup file is not compatible"
  end if

!START-TURB
 if(ivd%turbulenceModel==0) then 
  read(INFILE) (grp%u(i,1),i=1,grp%numberOfNodes),&
                 ((grp%u(i,j),i=1,grp%numberOfNodes),j=2,5)
  do i=1,grp%numberOfNodes
   grp%u(i,2:5) = grp%u(i,1)*grp%u(i,2:5)
  end do
 end if

 if(ivd%turbulenceModel==1) then 
  read(INFILE) (grp%u(i,1),i=1,grp%numberOfNodes),&
                 ((grp%u(i,j),i=1,grp%numberOfNodes),j=2,6)
  do i=1,grp%numberOfNodes
   grp%u(i,2:5) = grp%u(i,1)*grp%u(i,2:5)
  end do
 end if
 if(ivd%turbulenceModel==2) then 
  read(INFILE) (grp%u(i,1),i=1,grp%numberOfNodes),&
                 ((grp%u(i,j),i=1,grp%numberOfNodes),j=2,7)
  do i=1,grp%numberOfNodes
   grp%u(i,2:5) = grp%u(i,1)*grp%u(i,2:5)
  end do
 end if

if(ivd%turbulenceModel==3) then 
  read(INFILE) (grp%u(i,1),i=1,grp%numberOfNodes),&
                 ((grp%u(i,j),i=1,grp%numberOfNodes),j=2,7)
  do i=1,grp%numberOfNodes
   grp%u(i,2:5) = grp%u(i,1)*grp%u(i,2:5)
  end do
 end if
!END-TURB
 
 grp%uold = grp%u
 
 grp%uold2 = grp%u

!Start up from freestream condition on 2nd array
!while using the first as inflow boundary condition
else
! start up from freestream
 
 do i=1,grp%numberOfNodes
  grp%u(i,1) = ivd%inflowField(1,1)
  grp%u(i,2) = ivd%inflowField(2,1)
  grp%u(i,3) = ivd%inflowField(3,1)
  grp%u(i,4) = ivd%inflowField(4,1)
  grp%u(i,5) = ivd%inflowField(5,1)
  grp%uold(i,1) = ivd%inflowField(1,1)
  grp%uold(i,2) = ivd%inflowField(2,1)
  grp%uold(i,3) = ivd%inflowField(3,1)
  grp%uold(i,4) = ivd%inflowField(4,1)
  grp%uold(i,5) = ivd%inflowField(5,1)
  grp%uold2(i,1) = ivd%inflowField(1,1)
  grp%uold2(i,2) = ivd%inflowField(2,1)
  grp%uold2(i,3) = ivd%inflowField(3,1)
  grp%uold2(i,4) = ivd%inflowField(4,1)
  grp%uold2(i,5) = ivd%inflowField(5,1)
  enddo
  
  if (ivd%patchInitialization) then
  SELECT CASE (ivd%patchType)
  CASE (1)
   write(*,*) "initialize box patch"
   do i=1,grp%numberOfNodes
    dist = -MIN(MIN(MIN(grp%coordinates(i,1)-ivd%patchBox(1,1),-grp%coordinates(i,1)+ivd%patchBox(1,2)), &
        MIN(grp%coordinates(i,2)-ivd%patchBox(2,1),-grp%coordinates(i,2)+ivd%patchBox(2,2))), &
        MIN(grp%coordinates(i,3)-ivd%patchBox(3,1),-grp%coordinates(i,3)+ivd%patchBox(3,2)))
    if (dist .le. 0.0) then 
     grp%u(i,1) = ivd%inflowField(1,2)
     grp%u(i,2) = ivd%inflowField(2,2)
     grp%u(i,3) = ivd%inflowField(3,2)
     grp%u(i,4) = ivd%inflowField(4,2)
     grp%u(i,5) = ivd%inflowField(5,2)
     grp%uold(i,1) = ivd%inflowField(1,2)
     grp%uold(i,2) = ivd%inflowField(2,2)
     grp%uold(i,3) = ivd%inflowField(3,2)
     grp%uold(i,4) = ivd%inflowField(4,2)
     grp%uold(i,5) = ivd%inflowField(5,2)
     grp%uold2(i,1) = ivd%inflowField(1,2)
     grp%uold2(i,2) = ivd%inflowField(2,2)
     grp%uold2(i,3) = ivd%inflowField(3,2)
     grp%uold2(i,4) = ivd%inflowField(4,2)
     grp%uold2(i,5) = ivd%inflowField(5,2)
    endif
   enddo
  CASE (2)
   write(*,*) "initialize spherical patch"
   do i=1,grp%numberOfNodes
    !Spherical pressure
    dist = sqrt((grp%coordinates(i,1)-ivd%patchSphere(1))**2 &
     +(grp%coordinates(i,2)-ivd%patchSphere(2))**2 &
     +(grp%coordinates(i,3)-ivd%patchSphere(3))**2)
    if (dist .le. ivd%patchSphere(4)) then
     grp%u(i,1) = ivd%inflowField(1,2)
     grp%u(i,2) = ivd%inflowField(2,2)
     grp%u(i,3) = ivd%inflowField(3,2)
     grp%u(i,4) = ivd%inflowField(4,2)
     grp%u(i,5) = ivd%inflowField(5,2)
     grp%uold(i,1) = ivd%inflowField(1,2)
     grp%uold(i,2) = ivd%inflowField(2,2)
     grp%uold(i,3) = ivd%inflowField(3,2)
     grp%uold(i,4) = ivd%inflowField(4,2)
     grp%uold(i,5) = ivd%inflowField(5,2)
     grp%uold2(i,1) = ivd%inflowField(1,2)
     grp%uold2(i,2) = ivd%inflowField(2,2)
     grp%uold2(i,3) = ivd%inflowField(3,2)
     grp%uold2(i,4) = ivd%inflowField(4,2)
     grp%uold2(i,5) = ivd%inflowField(5,2)
    endif
   enddo
  CASE DEFAULT
   write(*,*) "invalid patch type", ivd%patchType
  END SELECT
  endif

 
 if(ivd%ReynoldsNumber>0.0) then 
  grp%laminarViscosity = 1./ivd%ReynoldsNumber 
 else
  grp%laminarViscosity = 0.0
 end if

!START-TURB
 if(ivd%turbulenceModel==1) then
  grp%u(:,6) = 0.1
 end if 
if(ivd%turbulenceModel==2) then !MOD
  do i=1,grp%numberOfNodes
   grp%u(i,7) = 6.0*ivd%nuwall/(betaomega*(grp%wallDistance(i)**2.0)*ivd%ReynoldsNumber)
   grp%u(i,7) = max(grp%u(i,7),ivd%winf)
  end do
  grp%u(:,6) = ivd%kinfoverUinf2
 end if 
 if(ivd%turbulenceModel==3) then 
  do i=1,grp%numberOfNodes
   grp%u(i,7) = 60.0*ivd%nuwall/(betaomega*(grp%wallDistance(i)**2.0)*ivd%ReynoldsNumber)
   grp%u(i,7) = max(grp%u(i,7),ivd%winf)
  end do
  grp%u(:,6) = ivd%kinfoverUinf2
 end if 
!END-TURB 
end if

call makePressureField(grp,ivd,grp%u)
grp%cycleNumber = 1
physicalTimestepNumber = 0
call setBCsOnSolutionField(grp,ivd,grp%u,physicalTimestepNumber)!START-GUST

end subroutine setUpInitialField
!-------------------------------------------------------------------------

subroutine setUpInitialField2(ivd,grp,INFILE)!START-GUST
! parallel version 
IMPLICIT NONE
include 'mpif.h'

type(InputVariablesData) :: ivd
type(GridSolverData) :: grp
integer :: INFILE
!START-TURB
double precision :: xu(7)
!END-TURB
integer :: i,j,number,locnumb
integer :: allocateStatus

integer,pointer :: reg(:)
integer :: startupType

integer :: ind

real :: dist
INTEGER :: physicalTimestepNumber!START-GUST
!START-TURB
integer :: icom
real :: betaomega
betaomega = 3.0/40.0
!END-TURB

if(INFILE>0) then
! start up from given field
  write(*,*) "starting up from given field"

  if(grp%pdp%numberOfProcesses==1) then
    read(INFILE) number
    if(number.ne.grp%numberOfNodes) then
      write(*,*) number,grp%numberOfNodes
      STOP "ERROR: Startup file is not compatible"
    end if
  !START-TURB
   if(ivd%turbulenceModel==0) then 
    read(INFILE) (grp%u(i,1),i=1,grp%numberOfNodes),&
                 ((grp%u(i,j),i=1,grp%numberOfNodes),j=2,5)
   end if
   if(ivd%turbulenceModel==1) then 
    read(INFILE) (grp%u(i,1),i=1,grp%numberOfNodes),&
                 ((grp%u(i,j),i=1,grp%numberOfNodes),j=2,6)
   end if
   if(ivd%turbulenceModel==2) then 
    read(INFILE) (grp%u(i,1),i=1,grp%numberOfNodes),&
                 ((grp%u(i,j),i=1,grp%numberOfNodes),j=2,7)
   end if
   if(ivd%turbulenceModel==3) then 
    read(INFILE) (grp%u(i,1),i=1,grp%numberOfNodes),&
                 ((grp%u(i,j),i=1,grp%numberOfNodes),j=2,7)
   end if
!END-TURB
  else
  
    read(INFILE) locnumb,number

    allocate(reg(locnumb),stat=allocateStatus)
    if (allocateStatus/=0) STOP "ERROR: Not enough memory in setUpInitialField"

    read(INFILE) (reg(i),i=1,locnumb)

!START-TURB
    if(ivd%turbulenceModel==0) then !MOD
    icom=5
!    read(INFILE) ((grp%u(reg(i),j),i=1,locnumb),j=1,5)
    end if
    if(ivd%turbulenceModel==1) then !MOD
    icom=6
!    read(INFILE) ((grp%u(reg(i),j),i=1,locnumb),j=1,6)
    end if
    if(ivd%turbulenceModel==2) then !MOD
    icom=7
!    read(INFILE) ((grp%u(reg(i),j),i=1,locnumb),j=1,7)
    end if
    if(ivd%turbulenceModel==3) then 
    icom=7
!    read(INFILE) ((grp%u(reg(i),j),i=1,locnumb),j=1,7)
    end if

    read(INFILE) ((grp%u(reg(i),j),i=1,locnumb),j=1,icom)
!   do i=1,locnumb                       
!    read(INFILE) (grp%u(reg(i),j),j=1,icom)
!   end do

!END-TURB       

    write(*,*) "finished reading..."
 

    deallocate(reg,stat=allocateStatus)
    if (allocateStatus/=0) STOP "ERROR: Could not dellocate in setUpInitialField"
 
 ! communicate
 
#ifdef PARALLEL
 
    do i=1,grp%pdp%numberOfProcesses
      if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i)>0)) then
        call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
        position = 1 
!START-TURB
        do j=1,icom
!END-TURB
          call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                            grp%pdp%realType,grp%u(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
        end do
        call mp_sendv(grp%pdp%processorIDs,i,163,grp%pdp%sendFlag)
      end if
    end do
 
    call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,&
                        grp%pdp%receivedMessageFromProcess)
 
    do i  = 1, grp%pdp%numberOfProcesses
      if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
        call mp_recv( grp%pdp%processorIDs,i,163,grp%pdp%sendFlag )
      endif
    enddo

    call mp_wait_comms( grp%pdp%sendFlag )

    do i=1,grp%pdp%numberOfProcesses
      if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
!START-TURB
          do j = 1,icom
!END-TURB
            call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                               grp%pdp%realType,grp%u(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
          end do
      end if
    end do
 
    call MPI_BARRIER( MPI_COMM_WORLD, i )
#endif PARALLEL
 
    write(*,*) "finished communicating" 
  end if
 
!START-TURB
if(ivd%turbulenceModel==0) then 
  do i=1,grp%numberOfNodes
    grp%u(i,2:5) = grp%u(i,1)*grp%u(i,2:5)
  end do
end if
if(ivd%turbulenceModel==1) then 
  do i=1,grp%numberOfNodes
    grp%u(i,2:6) = grp%u(i,1)*grp%u(i,2:6)
  end do
end if
if(ivd%turbulenceModel==2 .or. ivd%turbulenceModel==3) then 
  do i=1,grp%numberOfNodes
    grp%u(i,2:5) = grp%u(i,1)*grp%u(i,2:5)
  end do
end if
!END-TURB

  grp%uold=grp%u
  grp%uold2=grp%u
 
else
 ! start up from freestream
 write(*,*) "starting up from freestream"
 write(*,'(5(1e15.7))') (ivd%inflowField(i,1),i=1,5)
 write(*,'(5(1e15.7))') (ivd%inflowField(i,2),i=1,5)
 do i=1,grp%numberOfNodes
  grp%u(i,1) = ivd%inflowField(1,1)
  grp%u(i,2) = ivd%inflowField(2,1)
  grp%u(i,3) = ivd%inflowField(3,1)
  grp%u(i,4) = ivd%inflowField(4,1)
  grp%u(i,5) = ivd%inflowField(5,1)
  grp%uold(i,1) = ivd%inflowField(1,1)
  grp%uold(i,2) = ivd%inflowField(2,1)
  grp%uold(i,3) = ivd%inflowField(3,1)
  grp%uold(i,4) = ivd%inflowField(4,1)
  grp%uold(i,5) = ivd%inflowField(5,1)
  grp%uold2(i,1) = ivd%inflowField(1,1)
  grp%uold2(i,2) = ivd%inflowField(2,1)
  grp%uold2(i,3) = ivd%inflowField(3,1)
  grp%uold2(i,4) = ivd%inflowField(4,1)
  grp%uold2(i,5) = ivd%inflowField(5,1)
        enddo
  
  if (ivd%patchInitialization) then
  SELECT CASE (ivd%patchType)
  CASE (1)
   write(*,*) "initialize box patch"
   do i=1,grp%numberOfNodes
    dist = -MIN(MIN(MIN(grp%coordinates(i,1)-ivd%patchBox(1,1),-grp%coordinates(i,1)+ivd%patchBox(1,2)), &
        MIN(grp%coordinates(i,2)-ivd%patchBox(2,1),-grp%coordinates(i,2)+ivd%patchBox(2,2))), &
        MIN(grp%coordinates(i,3)-ivd%patchBox(3,1),-grp%coordinates(i,3)+ivd%patchBox(3,2)))
    if (dist .le. 0.0) then 
     grp%u(i,1) = ivd%inflowField(1,2)
     grp%u(i,2) = ivd%inflowField(2,2)
     grp%u(i,3) = ivd%inflowField(3,2)
     grp%u(i,4) = ivd%inflowField(4,2)
     grp%u(i,5) = ivd%inflowField(5,2)
     grp%uold(i,1) = ivd%inflowField(1,2)
     grp%uold(i,2) = ivd%inflowField(2,2)
     grp%uold(i,3) = ivd%inflowField(3,2)
     grp%uold(i,4) = ivd%inflowField(4,2)
     grp%uold(i,5) = ivd%inflowField(5,2)
     grp%uold2(i,1) = ivd%inflowField(1,2)
     grp%uold2(i,2) = ivd%inflowField(2,2)
     grp%uold2(i,3) = ivd%inflowField(3,2)
     grp%uold2(i,4) = ivd%inflowField(4,2)
     grp%uold2(i,5) = ivd%inflowField(5,2)
    endif
   enddo
  CASE (2)
   write(*,*) "initialize spherical patch"
   do i=1,grp%numberOfNodes
    !Spherical pressure
    dist = sqrt((grp%coordinates(i,1)-ivd%patchSphere(1))**2 &
     +(grp%coordinates(i,2)-ivd%patchSphere(2))**2 &
     +(grp%coordinates(i,3)-ivd%patchSphere(3))**2)
    if (dist .le. ivd%patchSphere(4)) then
     grp%u(i,1) = ivd%inflowField(1,2)
     grp%u(i,2) = ivd%inflowField(2,2)
     grp%u(i,3) = ivd%inflowField(3,2)
     grp%u(i,4) = ivd%inflowField(4,2)
     grp%u(i,5) = ivd%inflowField(5,2)
     grp%uold(i,1) = ivd%inflowField(1,2)
     grp%uold(i,2) = ivd%inflowField(2,2)
     grp%uold(i,3) = ivd%inflowField(3,2)
     grp%uold(i,4) = ivd%inflowField(4,2)
     grp%uold(i,5) = ivd%inflowField(5,2)
     grp%uold2(i,1) = ivd%inflowField(1,2)
     grp%uold2(i,2) = ivd%inflowField(2,2)
     grp%uold2(i,3) = ivd%inflowField(3,2)
     grp%uold2(i,4) = ivd%inflowField(4,2)
     grp%uold2(i,5) = ivd%inflowField(5,2)
    endif
   enddo
  CASE DEFAULT
   write(*,*) "invalid patch type", ivd%patchType
  END SELECT
  endif
  
  if(ivd%ReynoldsNumber>0.0) then
    grp%laminarViscosity = 1./ivd%ReynoldsNumber
  else
    grp%laminarViscosity = 0.0
  end if

!START-TURB
  if(ivd%turbulenceModel==1) then
    grp%u(:,6) = 0.1
  end if
  if(ivd%turbulenceModel==2) then !MOD
  do i=1,grp%numberOfNodes
   grp%u(i,7) = 6.0*ivd%nuwall/(betaomega*(grp%wallDistance(i)**2.0)*ivd%ReynoldsNumber)
   grp%u(i,7) = max(grp%u(i,7),ivd%winf)
  end do
   grp%u(:,6) = ivd%kinfoverUinf2
  end if
  if(ivd%turbulenceModel==3) then
  do i=1,grp%numberOfNodes
   grp%u(i,7) = 60.0*ivd%nuwall/(betaomega*(grp%wallDistance(i)**2.0)*ivd%ReynoldsNumber)
   grp%u(i,7) = max(grp%u(i,7),ivd%winf)
  end do
   grp%u(:,6) = ivd%kinfoverUinf2
  end if
!END-TURB

  grp%uold=grp%u
  grp%uold2=grp%u

  call makePressureField(grp,ivd,grp%u)

 if(ivd%turbulenceModel==1) then
  call triggerTurbulenceField(grp,ivd)
 end if
end if

call makePressureField(grp,ivd,grp%u)

grp%cycleNumber = 1
physicalTimestepNumber = 0
call setBCsOnSolutionField(grp,ivd,grp%u,physicalTimestepNumber)!START-GUST

end subroutine setUpInitialField2

!-------------------------------------------------------------------------
subroutine triggerTurbulenceField(grp,ivd)
IMPLICIT NONE
! set turbulence field to a specified value around the trigger lines 

type(GridSolverData) :: grp
type(InputVariablesData) :: ivd

 integer :: i,j,ind

! write(*,*) "triggering"

 do i=1,grp%numberOfTripFieldNodes
  ind = grp%tripNodeFieldIndexes(i,1)
  if(grp%tripNodeFieldDistances(i)<ivd%triggerRadius) then 
   grp%u(ind,6) = ivd%turbulenceTriggerValue
  end if
 end do
end subroutine triggerTurbulenceField
!-------------------------------------------------------------------------
subroutine setOuterBoundaryConditions(grp,ivd,u,p)
! sets BC's at outer boundaries
IMPLICIT NONE

type(GridSolverData) :: grp
type(InputVariablesData) :: ivd
real :: u(:,:),p(:)
 
integer :: ist,ien,ib,ip,i1,i2,j,k,ri
real :: gammaMinusOne,inpRho,inpU1,inpU2,inpU3,inpE,inpPressure,normal(3),normalVelocity
real :: rho,U1,U2,U3,E,ww,s,t,xs,xc,x1,x2,ys,yc,y1,y2,sx,sy,rx,ry
real :: outRho,outU1,outU2,outU3,outE,outPressure
logical :: intersectionIsValid,cantFind,doPlot,solutionFound
real :: A,B,C,radicand,v1,v2,t1,t2

!nodeIndicatorRegistor
!-20,-19: wall
!-18,-17: sym
!-16,-15: outflow 3
!-14,-13: outflow 4
!-12,-11: engin in 1
!-10,-9: engin in 2
!-8,-7: engin out 1
!-6,-5: engin out 2
!-4,-3: inviscid wall (slip)

gammaMinusOne = ivd%gamma-1.0
ist = grp%brp%nodeIndicatorRegister(-16)
ien = grp%brp%nodeIndicatorRegister(-15)

!INFLOW CONDITIONS
 inpRho = ivd%inflowField(1,1)
 inpU1 = ivd%inflowField(2,1)
 inpU2 = ivd%inflowField(3,1)
 inpU3 = ivd%inflowField(4,1)
 inpE = ivd%inflowField(5,1)
 inpPressure = ivd%inflowField(6,1)
   
do ib =ist,ien
 ip = grp%brp%nodeIndicatorArray(ib)
 normal = grp%brp%nodeNormalArray(ip,:)
 
 rho = u(ip,1) 
 U1 = u(ip,2) 
 U2 = u(ip,3) 
 U3 = u(ip,4) 
 normalVelocity = inpU1*normal(1) + inpU2*normal(2) + inpU3*normal(3)
 !normalVelocity = U1*normal(1) + U2*normal(2) + U3*normal(3)
 if(normalVelocity > 0.00) then 
 ! inflow
  if(ivd%MachNumber>1.0) then 
   ! supersonic inflow
   u(ip,1) = inpRho 
   u(ip,2) = inpU1*inpRho
   u(ip,3) = inpU2*inpRho
   u(ip,4) = inpU3*inpRho
   u(ip,5) = inpE*inpRho
!START-TURB
   if(ivd%turbulenceModel==1) then
    u(ip,6) = 0.1
   end if
   if(ivd%turbulenceModel==2) then
    u(ip,6) = ivd%kinfoverUinf2
    u(ip,7) = ivd%winf
   end if
   if(ivd%turbulenceModel==3) then 
    u(ip,6) = ivd%kinfoverUinf2
    u(ip,7) = ivd%winf
   end if
!END-TURB
  else
  ! subsonic inflow  - assume inviscid at farfield
   u(ip,1) = inpRho
   u(ip,2) = inpU1*inpRho
   u(ip,3) = inpU2*inpRho
   u(ip,4) = inpU3*inpRho
   !VT: FIXED here for computing total energy
   u(ip,5) = p(ip)/gammaMinusOne&
             +0.5*(inpU1**2+inpU2**2+inpU3**2)*inpRho
!START-TURB
   if(ivd%turbulenceModel==1) then
    u(ip,6) = 0.1
   end if
   if(ivd%turbulenceModel==2) then
    u(ip,6) = ivd%kinfoverUinf2
    u(ip,7) = ivd%winf
   end if
   if(ivd%turbulenceModel==3) then 
    u(ip,6) = ivd%kinfoverUinf2
    u(ip,7) = ivd%winf
   end if
!END-TURB
  end if
 else
 ! outflow 
  if(ivd%MachNumber<1.0) then 
  ! subsonic outflow
   u(ip,5) = inpPressure/gammaMinusOne + 0.5*(U1**2+U2**2+U3**2)/rho
!START-TURB
   if(ivd%turbulenceModel==1) then
    u(ip,6) = 0.1
   end if
   if(ivd%turbulenceModel==2) then
    u(ip,6) = ivd%kinfoverUinf2
    u(ip,7) = ivd%winf
   end if
   if(ivd%turbulenceModel==3) then
    u(ip,6) = ivd%kinfoverUinf2
    u(ip,7) = ivd%winf
   end if
!END-TURB
  end if

 ! for supersonic outflow, do nothing (assuming inviscid at farfield) 
!START-TURB 
  if(ivd%turbulenceModel==1) then 
   u(ip,6) = 0.1
  end if
  if(ivd%turbulenceModel==2) then
    u(ip,6) = ivd%kinfoverUinf2
    u(ip,7) = ivd%winf
  end if
  if(ivd%turbulenceModel==3) then 
    u(ip,6) = ivd%kinfoverUinf2
    u(ip,7) = ivd%winf
  end if
!END-TURB
 end if
end do

ist = grp%brp%nodeIndicatorRegister(-14)
ien = grp%brp%nodeIndicatorRegister(-13)

!INFLOW CONDITIONS
 inpRho = ivd%inflowField(1,2)
 inpU1 = ivd%inflowField(2,2)
 inpU2 = ivd%inflowField(3,2)
 inpU3 = ivd%inflowField(4,2)
 inpE = ivd%inflowField(5,2)
 inpPressure = ivd%inflowField(6,2)
   
do ib =ist,ien
 ip = grp%brp%nodeIndicatorArray(ib)
 normal = grp%brp%nodeNormalArray(ip,:)
 
 rho = u(ip,1) 
 U1 = u(ip,2) 
 U2 = u(ip,3) 
 U3 = u(ip,4) 
 normalVelocity = inpU1*normal(1) + inpU2*normal(2) + inpU3*normal(3)
 !normalVelocity = U1*normal(1) + U2*normal(2) + U3*normal(3)
 if(normalVelocity > 0.00) then 
 ! inflow
  if(ivd%MachNumber>1.0) then 
   ! supersonic inflow
   u(ip,1) = inpRho 
   u(ip,2) = inpU1*inpRho
   u(ip,3) = inpU2*inpRho
   u(ip,4) = inpU3*inpRho
   u(ip,5) = inpE*inpRho
!START-TURB 
  if(ivd%turbulenceModel==1) then 
   u(ip,6) = 0.1
  end if
  if(ivd%turbulenceModel==2) then
    u(ip,6) = ivd%kinfoverUinf2
    u(ip,7) = ivd%winf
  end if
  if(ivd%turbulenceModel==3) then 
    u(ip,6) = ivd%kinfoverUinf2
    u(ip,7) = ivd%winf
  end if
!END-TURB
  else
  ! subsonic inflow  - assume inviscid at farfield
   u(ip,1) = inpRho
   u(ip,2) = inpU1*inpRho
   u(ip,3) = inpU2*inpRho
   u(ip,4) = inpU3*inpRho
   !VT: FIXED here for computing total energy
   u(ip,5) = p(ip)/gammaMinusOne&
             +0.5*(inpU1**2+inpU2**2+inpU3**2)*inpRho
!START-TURB 
  if(ivd%turbulenceModel==1) then 
   u(ip,6) = 0.1
  end if
  if(ivd%turbulenceModel==2) then
    u(ip,6) = ivd%kinfoverUinf2
    u(ip,7) = ivd%winf
  end if
  if(ivd%turbulenceModel==3) then 
    u(ip,6) = ivd%kinfoverUinf2
    u(ip,7) = ivd%winf
  end if
!END-TURB
  end if
 else
 ! outflow 
  if(ivd%MachNumber<1.0) then 
  ! subsonic outflow
   u(ip,5) = inpPressure/gammaMinusOne + 0.5*(U1**2+U2**2+U3**2)/rho
!START-TURB 
  if(ivd%turbulenceModel==1) then 
   u(ip,6) = 0.1
  end if
  if(ivd%turbulenceModel==2) then
    u(ip,6) = ivd%kinfoverUinf2
    u(ip,7) = ivd%winf
  end if
  if(ivd%turbulenceModel==3) then 
    u(ip,6) = ivd%kinfoverUinf2
    u(ip,7) = ivd%winf
  end if
!END-TURB
  end if

 ! for supersonic outflow, do nothing (assuming inviscid at farfield) 
!START-TURB 
  if(ivd%turbulenceModel==1) then 
   u(ip,6) = 0.1
  end if
  if(ivd%turbulenceModel==2) then
    u(ip,6) = ivd%kinfoverUinf2
    u(ip,7) = ivd%winf
  end if
  if(ivd%turbulenceModel==3) then 
    u(ip,6) = ivd%kinfoverUinf2
    u(ip,7) = ivd%winf
  end if
!END-TURB
 end if
end do

end subroutine setOuterBoundaryConditions
!-------------------------------------------------------------------------

subroutine nullifyOuterBoundary(grp,u)
! set variable to zero on boundary (used for multigrid)
IMPLICIT NONE

type(GridSolverData) :: grp
type(InputVariablesData) :: ivd
real :: u(:,:)

integer :: ist,ien,ib,ip

ist = grp%brp%nodeIndicatorRegister(-16)
ien = grp%brp%nodeIndicatorRegister(-15)

do ib =ist,ien
 ip = grp%brp%nodeIndicatorArray(ib)

!START-TURB
if(ivd%turbulenceModel==1) then
 u(ip,1:6) = 0.0
end if
if(ivd%turbulenceModel==2) then 
 u(ip,1:6) = 0.0 !omega is not set to zero
end if
if(ivd%turbulenceModel==3) then 
 u(ip,1:6) = 0.0 !omega is not set to zero
end if
!END-TURB

end do


ist = grp%brp%nodeIndicatorRegister(-12)
ien = grp%brp%nodeIndicatorRegister(-11)

do ib =ist,ien
 ip = grp%brp%nodeIndicatorArray(ib)

!START-TURB
if(ivd%turbulenceModel==1) then
 u(ip,1:6) = 0.0
end if
if(ivd%turbulenceModel==2) then
 u(ip,1:6) = 0.0 !omega is not set to zero
end if
if(ivd%turbulenceModel==3) then 
 u(ip,1:6) = 0.0 !omega is not set to zero
end if
!END-TURB

end do

ist = grp%brp%nodeIndicatorRegister(-10)
ien = grp%brp%nodeIndicatorRegister(-9)

do ib =ist,ien
 ip = grp%brp%nodeIndicatorArray(ib)

!START-TURB
if(ivd%turbulenceModel==1) then
 u(ip,1:6) = 0.0
end if
if(ivd%turbulenceModel==2) then 
 u(ip,1:6) = 0.0 !omega is not set to zero
end if
if(ivd%turbulenceModel==3) then 
 u(ip,1:6) = 0.0 !omega is not set to zero
end if
!END-TURB

end do

ist = grp%brp%nodeIndicatorRegister(-8)
ien = grp%brp%nodeIndicatorRegister(-7)

do ib =ist,ien
 ip = grp%brp%nodeIndicatorArray(ib)
!START-TURB
if(ivd%turbulenceModel==1) then
 u(ip,1:6) = 0.0
end if
if(ivd%turbulenceModel==2) then  
 u(ip,1:6) = 0.0 !omega is not set to zero
end if
if(ivd%turbulenceModel==3) then 
 u(ip,1:6) = 0.0 !omega is not set to zero
end if
!END-TURB
end do

ist = grp%brp%nodeIndicatorRegister(-6)
ien = grp%brp%nodeIndicatorRegister(-5)

do ib =ist,ien
 ip = grp%brp%nodeIndicatorArray(ib)

!START-TURB
if(ivd%turbulenceModel==1) then
 u(ip,1:6) = 0.0
end if
if(ivd%turbulenceModel==2) then  
 u(ip,1:6) = 0.0 !omega is not set to zero
end if
if(ivd%turbulenceModel==3) then 
 u(ip,1:6) = 0.0 !omega is not set to zero
end if
!END-TURB
end do

end subroutine nullifyOuterBoundary
!-------------------------------------------------------------------------
subroutine setBCsOnSolutionField(grp,ivd,u,physicalTimestepNumber)
! sets boundary conditions on the solution field
IMPLICIT NONE
include 'mpif.h'

type(GridSolverData) :: grp
type(InputVariablesData) :: ivd
real :: u(:,:)

integer :: ist,ien,ib,ip,nbf,ind1,ind2,i,inletNumber,j
real :: normal(3),tangent(3),normalVelocity,tangentialVelocityVector(3),ww,vectorNorm
real :: tangentialVelocity,press1,press2,press,temp1,temp2,temp,ptfunc(2),area(2)
real :: wallVelocity(3),normalWallVelocityComp
real :: engineMassFlows(2),wx(3),u1(3),u2(3),engineMassFactor(2)
real :: intE,xs(3),ftmp,massflowfunction(2)
real :: factV(3),relaxationParam,factor1,factor2,storeMassFactor(2)
INTEGER :: physicalTimestepNumber,IPRINT,nper!START-GUST
REAL :: xmax,PI,Ug
!START-TURB
real :: betaomega
betaomega = 3.0/40.0
!END-TURB

! Boundary Condition relaxation

  if(ivd%numberOfCFLIncrements>0) then
    relaxationParam = dble(grp%cycleNumber-1)/dble(ivd%numberOfCFLIncrements)
    relaxationParam = min(relaxationParam,1.0)
  else
    relaxationParam = 1.0
  end if

! wall 

  ist = grp%brp%nodeIndicatorRegister(-20)
  ien = grp%brp%nodeIndicatorRegister(-19)
  
  if(ivd%ReynoldsNumber>0) then  ! this means viscous flow
    do ib =ist,ien
      ip = grp%brp%nodeIndicatorArray(ib)
      if(ivd%movingwall==.true.)then              ! bje jan11
        tangent=grp%brp%normtan(ib,4:6)
        u(ip,2:4)=(1.-relaxationParam)*u(ip,2:4)+relaxationParam*tangent*u(ip,1)+u(ip,1)*grp%coordinateMovement(ip,1:3)/grp%physicalTimestep
      else
        u(ip,2:4) = u(ip,1)*grp%coordinateMovement(ip,1:3)/grp%physicalTimestep
      endif
      if(.not.ivd%wallsAreIsentropic) then 
        ! set temperature 
        u(ip,5) = u(ip,1)*ivd%inflowTemperature/ivd%gamma 
      end if

!START-TURB
   if(ivd%turbulenceModel==1) then
      u(ip,6) = 0.0
   end if
   if(ivd%turbulenceModel==2) then 
      u(ip,6) = 0.0
      u(ip,7) = 6.0*ivd%nuwall/(betaomega*(ivd%firstoffwallpoint**2.0)*ivd%ReynoldsNumber)
   end if
  if(ivd%turbulenceModel==3) then 
      u(ip,6) = 0.0
      u(ip,7) = 60.0*ivd%nuwall/(betaomega*(ivd%firstoffwallpoint**2.0)*ivd%ReynoldsNumber)
   end if
!END-TURB

    end do
  else
    do ib =ist,ien
      ip = grp%brp%nodeIndicatorArray(ib)
      normal = grp%brp%nodeNormalArray(ip,:)
      if(grp%gridNumber==1.and.ivd%addTimeStep) then 
       wallVelocity(1) = grp%coordinateMovement(ip,1)/grp%physicalTimestep
       wallVelocity(2) = grp%coordinateMovement(ip,2)/grp%physicalTimestep
       wallVelocity(3) = grp%coordinateMovement(ip,3)/grp%physicalTimestep
      else
       wallVelocity = 0.0
      end if
      normalWallVelocityComp = sum(wallVelocity*normal)
      normalVelocity = sum(u(ip,2:4)*normal(:))
      tangentialVelocityVector(:) = u(ip,2:4)-normalVelocity*normal(:) 
      u(ip,2:4) = tangentialVelocityVector(:)+grp%u(ip,1)*normalWallVelocityComp*normal(:)
!START-TURB
    if(ivd%turbulenceModel > 0) then
      u(ip,6) = 0.0
    end if
!END-TURB
    end do
  end if

! Inviscid Wall

  ist = grp%brp%nodeIndicatorRegister(-4)
  ien = grp%brp%nodeIndicatorRegister(-3)
  
  do ib =ist,ien
    ip = grp%brp%nodeIndicatorArray(ib)
    normal = grp%brp%nodeNormalArray(ip,:)
    normalVelocity = sum(u(ip,2:4)*normal(:))
    tangentialVelocityVector(:) = u(ip,2:4)-relaxationParam*normalVelocity*normal(:)
    u(ip,2:4) = tangentialVelocityVector(:) 
    if(.not.ivd%wallsAreIsentropic) then
      ! set temperature 
      u(ip,5) = u(ip,1)*ivd%inflowTemperature/ivd%gamma + 0.5 *sum(u(ip,2:4)*u(ip,2:4))/u(ip,1)              
    end if
!START-TURB
    if(ivd%turbulenceModel > 0) then
      u(ip,6) = 0.0
    end if
!END-TURB
  end do

! symmetry

  ist = grp%brp%nodeIndicatorRegister(-18)
  ien = grp%brp%nodeIndicatorRegister(-17)

  do ib =ist,ien
    ip = grp%brp%nodeIndicatorArray(ib)
    normal = grp%brp%nodeNormalArray(ip,:)
    normalVelocity = sum(u(ip,2:4)*normal(:))
    tangentialVelocityVector(:) = u(ip,2:4)-normalVelocity*normal(:) 
    u(ip,2:4) = tangentialVelocityVector(:) 
  end do

! engine outlet relaxation

  if(ivd%numberOfEORelaxationSteps>0) then 
    relaxationParam = dble(grp%cycleNumber-1)/dble(ivd%numberOfEORelaxationSteps)
    relaxationParam = min(relaxationParam,1.0)
  else
    relaxationParam = 1.0
  end if

! engine outflow - 1

  ist = grp%brp%nodeIndicatorRegister(-8)
  ien = grp%brp%nodeIndicatorRegister(-7)

  do ib =ist,ien
    ip = grp%brp%nodeIndicatorArray(ib)
    u(ip,1:5) = (1.0-relaxationParam)*ivd%inflowField(1:5,1)+relaxationParam*ivd%enginesRearParameters(1:5,1)
    grp%p(ip) = (1.0-relaxationParam)*ivd%inflowField(6,1)+relaxationParam*ivd%enginesRearParameters(6,1)
!START-TURB
  if(ivd%turbulenceModel==1) then
    u(ip,6) = relaxationParam*ivd%enginesRearParameters(7,1)
  end if
  if(ivd%turbulenceModel==2 .or. ivd%turbulenceModel==3) then
 !   u(ip,6) = relaxationParam*ivd%kinfoverUinf2
 !   u(ip,7) = relaxationParam*ivd%winf
    u(ip,6) = relaxationParam*ivd%enginesRearParameters(7,2) !BJE (JAN 11)
  end if
!END-TURB
  end do

! engine outflow - 2

  ist = grp%brp%nodeIndicatorRegister(-6)
  ien = grp%brp%nodeIndicatorRegister(-5)

  do ib =ist,ien
    ip = grp%brp%nodeIndicatorArray(ib)
    u(ip,1:5) = (1.0-relaxationParam)*ivd%inflowField(1:5,1)+relaxationParam*ivd%enginesRearParameters(1:5,2)
    grp%p(ip) = (1.0-relaxationParam)*ivd%inflowField(6,1)+relaxationParam*ivd%enginesRearParameters(6,2)
!START-TURB
  if(ivd%turbulenceModel==1) then
    u(ip,6) = relaxationParam*ivd%enginesRearParameters(7,2)
  end if
  if(ivd%turbulenceModel==2 .or. ivd%turbulenceModel==3) then
 !   u(ip,6) = relaxationParam*ivd%kinfoverUinf2
 !   u(ip,7) = relaxationParam*ivd%winf
    u(ip,6) = relaxationParam*ivd%enginesRearParameters(7,2) !BJE (jan 11)
  end if
!END-TURB
  end do

! engine inlet - original: imposing an absolute mass flow rate (Kg/s)

  if(ivd%engineFlowType==1)then
   engineMassFlows = 0.
   if(grp%gridNumber==1) then
     do i=1,grp%brp%numberOfEngineInletSides
       ind1 = grp%brp%engineInletSideIndexes(i,1)
       ind2 = grp%brp%engineInletSideIndexes(i,2)
       inletNumber = grp%brp%engineInletSideIndexes(i,3)
       wx = grp%brp%engineInletSideCoefficients(i,:)
       u1 = u(ind1,2:4)
       u2 = u(ind2,2:4)
       engineMassFlows(inletNumber) = engineMassFlows(inletNumber) + sum(u1*wx) + sum(u2*wx)
      end do
   else
     do i=1,grp%brp%numberOfEngineInletSides
       ind1 = grp%brp%engineInletSideIndexes(i,1)
       ind2 = grp%brp%engineInletSideIndexes(i,2)
       inletNumber = grp%brp%engineInletSideIndexes(i,3)
       wx = grp%brp%engineInletSideCoefficients(i,:)
       engineMassFlows(inletNumber) = engineMassFlows(inletNumber) + 2.0*sqrt(sum(wx*wx))
     end do
   end if

#ifdef PARALLEL

   call mp_init_buffer(grp%pdp%processorIDs,-1,grp%pdp%sendFlag)
   position = 1
   call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,engineMassFlows(1),grp%pdp%sendFlag)
   call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,engineMassFlows(2),grp%pdp%sendFlag)

   call mp_snd_others(grp%pdp%numberOfProcesses,grp%pdp%currentDomain,grp%pdp%processorIDs,367,grp%pdp%sendFlag)

   call init_recv_others(grp%pdp%numberOfProcesses,grp%pdp%currentDomain,grp%pdp%receivedMessageFromProcess)

   call mp_recv_others(grp%pdp%numberOfProcesses,grp%pdp%currentDomain,grp%pdp%processorIDs,367,grp%pdp%sendFlag)

   call mp_wait_comms( grp%pdp%sendFlag )


   do i=1,grp%pdp%numberOfProcesses
     if( i.ne.grp%pdp%currentDomain ) then
         call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
         engineMassFlows(1) = engineMassFlows(1) + ftmp
         call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
         engineMassFlows(2) = engineMassFlows(2) + ftmp
     end if
   end do


   call MPI_BARRIER( MPI_COMM_WORLD, i )
  
#endif PARALLEL

   if(engineMassFlows(1)>0.0) then
       engineMassFactor(1) = ivd%enginesFrontMassFlow(1)/engineMassFlows(1)
   end if
   if(engineMassFlows(2)>0.0) then
       engineMassFactor(2) = ivd%enginesFrontMassFlow(2)/engineMassFlows(2)
   end if


   if(grp%pdp%currentDomain==1) then 
     if(grp%gridNumber==1) then
       if(ivd%enginesFrontMassFlow(1)>0.0) then
         write(*,*) "Engine mass flow 1 (prescribed): ",engineMassFlows(1),ivd%enginesFrontMassFlow(1)
       end if
 
       if(ivd%enginesFrontMassFlow(2)>0.0) then
         write(*,*) "Engine mass flow 2 (prescribed): ",engineMassFlows(2),ivd%enginesFrontMassFlow(2)
       end if
     end if
   end if

! engine inlet - 1

   if(grp%gridNumber==1) then
     ist = grp%brp%nodeIndicatorRegister(-12)
     ien = grp%brp%nodeIndicatorRegister(-11)

     do ib =ist,ien
       ip = grp%brp%nodeIndicatorArray(ib)
       factV = grp%brp%engineInletNormals(1,:)*ivd%enginesFrontMassFlow(1)/grp%brp%engineInletAreas(1)
       intE = grp%u(ip,5) - 0.5*sum(grp%u(ip,2:4)*grp%u(ip,2:4))/grp%u(ip,1)
       u(ip,2:4) = u(ip,2:4)*engineMassFactor(1)
       u(ip,5) = intE + 0.5*sum(u(ip,2:4)*u(ip,2:4))/u(ip,1)
!START-TURB
     if(ivd%turbulenceModel > 0) then
       u(ip,6) = 0.0
     end if
!END-TURB
     end do

! engine inlet - 2

     ist = grp%brp%nodeIndicatorRegister(-10)
     ien = grp%brp%nodeIndicatorRegister(-9)

     do ib =ist,ien
       ip = grp%brp%nodeIndicatorArray(ib)
       factV = grp%brp%engineInletNormals(2,:)*ivd%enginesFrontMassFlow(2)/grp%brp%engineInletAreas(2)
       intE = grp%u(ip,5) - 0.5*sum(grp%u(ip,2:4)*grp%u(ip,2:4))/grp%u(ip,1)
       u(ip,2:4) = u(ip,2:4)*engineMassFactor(2)
       u(ip,5) = intE + 0.5*sum(u(ip,2:4)*u(ip,2:4))/u(ip,1)
!START-TURB
     if(ivd%turbulenceModel > 0) then
       u(ip,6) = 0.0
     end if
!END-TURB
     end do
   end if

! make sure there are no negative velocities in engine inlet

   if(grp%gridNumber==1) then
     ist = grp%brp%nodeIndicatorRegister(-12)
     ien = grp%brp%nodeIndicatorRegister(-11)

     do ib =ist,ien
       ip = grp%brp%nodeIndicatorArray(ib)
       normal = grp%brp%nodeNormalArray(ip,:)
       normalVelocity = sum(-normal*u(ip,2:4))/u(ip,1)
       if(normalVelocity<0.00) then
         if(grp%pdp%currentDomain==1) &
           write(*,*) "WARNING: Normal velocity for engine inlet node ",ip," is positive ",normalVelocity,normal
         u(ip,2) = u(ip,2)+u(ip,1)*(normalVelocity-0.01)*normal(1)
         u(ip,3) = u(ip,3)+u(ip,1)*(normalVelocity-0.01)*normal(2)
         u(ip,4) = u(ip,4)+u(ip,1)*(normalVelocity-0.01)*normal(3)
       end if
     end do

! engine inlet - 2

     ist = grp%brp%nodeIndicatorRegister(-10)
     ien = grp%brp%nodeIndicatorRegister(-9)

     do ib =ist,ien
       ip = grp%brp%nodeIndicatorArray(ib)
       normal = grp%brp%nodeNormalArray(ip,:)
       normalVelocity = sum(-normal*u(ip,2:4))/u(ip,1)
       if(normalVelocity<0.00) then
         write(*,*) "WARNING: Normal velocity for engine inlet node ",ip," is positive"
         u(ip,2) = u(ip,2)+u(ip,1)*(normalVelocity-0.01)*normal(1)
         u(ip,3) = u(ip,3)+u(ip,1)*(normalVelocity-0.01)*normal(2)
         u(ip,4) = u(ip,4)+u(ip,1)*(normalVelocity-0.01)*normal(3)
       end if
     end do
   end if
! engine inlet - enforce the given mass flow function (mass flow * sqrt(temp) over total pressure)

  elseif(ivd%engineFlowType==2)then

! mass flow function iteration
   storeMassFactor(1) = 99999999999.
   storeMassFactor(2) = 99999999999.
 ! do j=1,100
    
   engineMassFlows = 0.
   area = 0.
   ptfunc = 0.
   if(grp%gridNumber==1) then
     do i=1,grp%brp%numberOfEngineInletSides
       ind1 = grp%brp%engineInletSideIndexes(i,1)
       ind2 = grp%brp%engineInletSideIndexes(i,2)
       inletNumber = grp%brp%engineInletSideIndexes(i,3)
       wx = grp%brp%engineInletSideCoefficients(i,:)
       u1 = u(ind1,2:4)
       u2 = u(ind2,2:4)
       engineMassFlows(inletNumber) = engineMassFlows(inletNumber) + sum(u1*wx) + sum(u2*wx)
       press1 = grp%p(ind1) + 0.5*sum(u(ind1,2:4)*u(ind1,2:4))/u(ind1,1)
       press2 = grp%p(ind2) + 0.5*sum(u(ind2,2:4)*u(ind2,2:4))/u(ind2,1)
       press = 0.5*(press1 + press2)
       temp1 = (ivd%gamma/u(ind1,1))*(u(ind1,5)-(0.5*sum(u(ind1,2:4)*u(ind1,2:4)))/u(ind1,1))
       temp2 = (ivd%gamma/u(ind2,1))*(u(ind2,5)-(0.5*sum(u(ind2,2:4)*u(ind2,2:4)))/u(ind2,1))
       temp = 0.5*(temp1+temp2)
       ptfunc(inletNumber) = ptfunc(inletNumber) + (sqrt(temp)/press)*2.0*sqrt(sum(wx*wx))
       area(inletNumber) = area(inletNumber) + 2.0*sqrt(sum(wx*wx))
     end do
    else
     do i=1,grp%brp%numberOfEngineInletSides
       ind1 = grp%brp%engineInletSideIndexes(i,1)
       ind2 = grp%brp%engineInletSideIndexes(i,2)
       inletNumber = grp%brp%engineInletSideIndexes(i,3)
       wx = grp%brp%engineInletSideCoefficients(i,:)
       engineMassFlows(inletNumber) = engineMassFlows(inletNumber) + 2.0*sqrt(sum(wx*wx))
     end do
   end if
   
#ifdef PARALLEL
  
!  write(get_unit(),*)  'Starting engine comms'
   call mp_init_buffer(grp%pdp%processorIDs,-1,grp%pdp%sendFlag)
   position = 1
   call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,engineMassFlows(1),grp%pdp%sendFlag)
   call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,engineMassFlows(2),grp%pdp%sendFlag)
   call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,ptfunc(1),grp%pdp%sendFlag)
   call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,ptfunc(2),grp%pdp%sendFlag)
   call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,area(1),grp%pdp%sendFlag)
   call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,area(2),grp%pdp%sendFlag)
    
   call mp_snd_others(grp%pdp%numberOfProcesses,grp%pdp%currentDomain,grp%pdp%processorIDs,367,grp%pdp%sendFlag)
    
   call init_recv_others(grp%pdp%numberOfProcesses,grp%pdp%currentDomain,grp%pdp%receivedMessageFromProcess)
    
   call mp_recv_others(grp%pdp%numberOfProcesses,grp%pdp%currentDomain,grp%pdp%processorIDs,367,grp%pdp%sendFlag)
 
   call mp_wait_comms( grp%pdp%sendFlag )
    
!2971 continue

   do i=1,grp%pdp%numberOfProcesses
!   if(grp%pdp%receivedMessageFromProcess(i).eq.0) then
     if( i.ne.grp%pdp%currentDomain ) then
!     grp%pdp%sendFlag = 0
!     call mp_check_mesg_arrived(grp%pdp%processorIDs,i,367,grp%pdp%sendFlag,bufferzz,bufferSize)
  
!     if(grp%pdp%sendFlag.gt.0) then
!       grp%pdp%receivedMessageFromProcess(i) = 1
!       position = 1
         call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
         engineMassFlows(1) = engineMassFlows(1) + ftmp
         call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
         engineMassFlows(2) = engineMassFlows(2) + ftmp
         call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
         ptfunc(1) = ptfunc(1) + ftmp
         call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
         ptfunc(2) = ptfunc(2) + ftmp
         call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
         area(1) = area(1) + ftmp
         call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
         area(2) = area(2) + ftmp
!     end if
     end if
   end do

   call MPI_BARRIER( MPI_COMM_WORLD, i )
! write(get_unit(),*)  'Finished engine comms'

#endif PARALLEL

   if(engineMassFlows(1)>0.0) then
!   if(grp%gridNumber==1) then
       massflowfunction(1) = engineMassFlows(1)*(ptfunc(1)/area(1))
       engineMassFactor(1) = ivd%enginesFrontMassFlow(1)/massflowfunction(1)
! if exceeding the max possible mass flow factor for this duct cross section print warning and leave the loop
       if(engineMassFactor(1).gt.storeMassFactor(1))then
         print*,'!!!WARNING!!! You are trying to exceed the maximum achievable mass flow function for engine 1'
         goto 9999
       endif
        storeMassFactor(1) = engineMassFactor(1)
!     else
!       engineMassFactor(1) = 1.0
!     end if
   end if
   if(engineMassFlows(2)>0.0) then
!     if(grp%gridNumber==1) then
         massflowfunction(2) = engineMassFlows(2)*(ptfunc(2)/area(2))
         engineMassFactor(2) = ivd%enginesFrontMassFlow(2)/massflowfunction(2)
! if exceeding the max possible mass flow factor for this duct cross section print warning and leave the loop
     if(engineMassFactor(2).gt.storeMassFactor(2))then
       print*,'!!!WARNING!!! You are trying to exceed the maximum achievable mass flow function for engine 2'
       goto 9999
     endif
        storeMassFactor(2) = engineMassFactor(2)
!    else
!       engineMassFactor(2) = 1.0
!    end if
   end if


   if(grp%pdp%currentDomain==1) then
     if(grp%gridNumber==1) then
       if(ivd%enginesFrontMassFlow(1)>0.0) then
         write(*,*) 'ENGINES INFLOW ITERATION NUMBER',j
         write(*,*) 'engine mass factor =',engineMassFactor(1)
         write(*,*) "Engine mass flow function 1 (prescribed): ",massflowfunction(1),ivd%enginesFrontMassFlow(1)
         write(*,*) "Engine 1 inlet area :", area(1)
       end if

       if(ivd%enginesFrontMassFlow(2)>0.0) then
         write(*,*) "Engine mass flow function 2 (prescribed): ",massflowfunction(2),ivd%enginesFrontMassFlow(2)
         write(*,*) "Engine 2 inlet area:", area(2)
       end if
     end if
   end if

! engine inlet - 1

   if(grp%gridNumber==1) then
     ist = grp%brp%nodeIndicatorRegister(-12)
     ien = grp%brp%nodeIndicatorRegister(-11)

     do ib =ist,ien
       ip = grp%brp%nodeIndicatorArray(ib)
       factV = grp%brp%engineInletNormals(1,:)*ivd%enginesFrontMassFlow(1)/grp%brp%engineInletAreas(1)
       intE = grp%u(ip,5) - 0.5*sum(grp%u(ip,2:4)*grp%u(ip,2:4))/grp%u(ip,1)
       u(ip,2:4) = u(ip,2:4)*engineMassFactor(1)
       u(ip,5) = intE + 0.5*sum(u(ip,2:4)*u(ip,2:4))/u(ip,1)
!START-SST
       if(ivd%turbulenceModel > 0) then
         u(ip,6) = 0.0
       end if
!END-SST
     end do

! engine inlet - 2

     ist = grp%brp%nodeIndicatorRegister(-10)
     ien = grp%brp%nodeIndicatorRegister(-9)

     do ib =ist,ien
       ip = grp%brp%nodeIndicatorArray(ib)
       factV = grp%brp%engineInletNormals(2,:)*ivd%enginesFrontMassFlow(2)/grp%brp%engineInletAreas(2)
       intE = grp%u(ip,5) - 0.5*sum(grp%u(ip,2:4)*grp%u(ip,2:4))/grp%u(ip,1)
       u(ip,2:4) = u(ip,2:4)*engineMassFactor(2)
       u(ip,5) = intE + 0.5*sum(u(ip,2:4)*u(ip,2:4))/u(ip,1)
!START-SST
       if(ivd%turbulenceModel > 0) then
         u(ip,6) = 0.0
       end if
     end do


! recompute the pressure field

     call makePressureField(grp,ivd,grp%u)

     if(ivd%enginesFrontMassFlow(1)>0.0)then
       factor1 = abs(massflowfunction(1) - ivd%enginesFrontMassFlow(1))/ivd%enginesFrontMassFlow(1)
      else
       factor1 = 0.0
     endif
     if(ivd%enginesFrontMassFlow(2)>0.0)then
       factor2 = abs(massflowfunction(2) - ivd%enginesFrontMassFlow(2))/ivd%enginesFrontMassFlow(2)
      else
       factor2 = 0.0
     endif

     if((factor1.lt.0.01).and.(factor2.lt.0.01))goto 9999

!  end do
   end if
9999 continue

! make sure there are no negative velocities in engine inlet

   if(grp%gridNumber==1) then
     ist = grp%brp%nodeIndicatorRegister(-12)
     ien = grp%brp%nodeIndicatorRegister(-11)

     do ib =ist,ien
       ip = grp%brp%nodeIndicatorArray(ib)
       normal = grp%brp%nodeNormalArray(ip,:)
       normalVelocity = sum(-normal*u(ip,2:4))/u(ip,1)
       if(normalVelocity<0.00) then
         if(grp%pdp%currentDomain==1) &
           write(*,*) "WARNING: Normal velocity for engine inlet node ",ip," is positive ",normalVelocity,normal
         u(ip,2) = u(ip,2)+u(ip,1)*(normalVelocity-0.01)*normal(1)
         u(ip,3) = u(ip,3)+u(ip,1)*(normalVelocity-0.01)*normal(2)
         u(ip,4) = u(ip,4)+u(ip,1)*(normalVelocity-0.01)*normal(3)
      end if
     end do

! engine inlet - 2

     ist = grp%brp%nodeIndicatorRegister(-10)
     ien = grp%brp%nodeIndicatorRegister(-9)

     do ib =ist,ien
       ip = grp%brp%nodeIndicatorArray(ib)
       normal = grp%brp%nodeNormalArray(ip,:)
       normalVelocity = sum(-normal*u(ip,2:4))/u(ip,1)
       if(normalVelocity<0.00) then
         write(*,*) "WARNING: Normal velocity for engine inlet node ",ip," is positive"
         u(ip,2) = u(ip,2)+u(ip,1)*(normalVelocity-0.01)*normal(1)
         u(ip,3) = u(ip,3)+u(ip,1)*(normalVelocity-0.01)*normal(2)
         u(ip,4) = u(ip,4)+u(ip,1)*(normalVelocity-0.01)*normal(3)
       end if
     end do
   end if
  end if
!START-GUST
! Gust generator models !4/2/2011
! Discrete turbulence
  PI = 4.0*atan(1.0)

  IPRINT=1
!!PRINT*,"physicalTimestepNumber",physicalTimestepNumber,grp%GUSTPAR(1),grp%GUSTPAR(2)
  if(ivd%GUSTMODEL)then
    if(physicalTimestepNumber==ivd%restartNumber) then
      if(grp%GUSTPAR(1)==1) then
        if(grp%GUSTPAR(2)==1) then
          nper=0
!!   PRINT*,"COORD",ivd%XPER(1),ivd%XPER(2),ivd%YPER(1),ivd%YPER(2),ivd%ZPER(1),ivd%ZPER(2)
          Do i= 1,grp%numberOfNodes
!          xmax=ivd%XPER(2) !!+grp%coordinates(i,2)*0.577 !considering an inclination of 30grad
           if (grp%coordinates(i,1).GE.ivd%XPER(1) .AND. grp%coordinates(i,1).LE.ivd%XPER(2)) then
             if (grp%coordinates(i,2).GE.ivd%YPER(1) .AND. grp%coordinates(i,2).LE.ivd%YPER(2)) then
               if (grp%coordinates(i,3).GE.ivd%ZPER(1) .AND. grp%coordinates(i,3).LE.ivd%ZPER(2)) then
                 nper=nper+1
                 grp%IPER(nper)=i
               end if
             end if
           end if
          end do
          grp%GUSTPAR(3)=nper
        end if
      end if
    end if
  
    if(physicalTimestepNumber==ivd%restartNumber .AND. grp%GUSTPAR(1)==1 .AND. grp%GUSTPAR(2)==1)PRINT*,"NPER",grp%GUSTPAR(3)," ",grp%numberOfNodes," ",grp%pdp%currentDomain," ",physicalTimestepNumber

    Ug = 0.0
!Prescribing the vertical gust velocity
    If(ivd%typegust==1) then
      Ug = ivd%UDS
      Do i= 1,grp%GUSTPAR(3)
       u(grp%IPER(i),4) = (ivd%inflowField(4,1)+Ug)*u(grp%IPER(i),1)
      end do
    end if
    if(ivd%typegust==2) then
      if(physicalTimestepNumber.le.ivd%TTS) then
        Ug = ivd%UDS/2.*(1.-cos(2.*PI*physicalTimestepNumber/ivd%TTS))
        Do i= 1,grp%GUSTPAR(3)
         u(grp%IPER(i),4) = (ivd%inflowField(4,1)+Ug)*u(grp%IPER(i),1)
        end do
      else
        Ug = 0.0
      end if
    end if

    if(grp%pdp%currentDomain==1.and.grp%GUSTPAR(1)==1.and.grp%GUSTPAR(2)==1)print*,"Ug=",Ug,"step=",physicalTimestepNumber
    if(grp%pdp%currentDomain==1.and.grp%GUSTPAR(1)==1.and.grp%GUSTPAR(2)==1)print*,"TTS=",ivd%TTS,"Uz(1and2)",u(grp%IPER(1),4),u(grp%IPER(2),4)
!END-GUST

  end if

end subroutine setBCsOnSolutionField 
!-------------------------------------------------------------------------
subroutine setBCsOnIncrementField(grp,ivd,rhs)
! sets BC's on the solution increment
IMPLICIT NONE

type(GridSolverData) :: grp
type(InputVariablesData) :: ivd
real :: rhs(:,:)

integer :: i,ist,ien,ib,ip
real :: normal(3),normalIncrement,tangentialIncrementVector(3),ww,vectorNorm,Tw
real :: tangent(3),tangentialIncrement

!
! wall
!

ist = grp%brp%nodeIndicatorRegister(-20)
ien = grp%brp%nodeIndicatorRegister(-19)

if(ivd%ReynoldsNumber>0) then 
 do ib =ist,ien
  ip = grp%brp%nodeIndicatorArray(ib)
  rhs(ip,2:4) = 0.0 
  if(.not.ivd%wallsAreIsentropic) then
   rhs(ip,5) = 0.0
  end if

!START-TURB
  if(ivd%turbulenceModel==1) then  
   rhs(ip,6) = 0.0
  end if 
  if(ivd%turbulenceModel==2) then 
   rhs(ip,6) = 0.0
   rhs(ip,7) = 0.0
  end if
if(ivd%turbulenceModel==3) then 
   rhs(ip,6) = 0.0
   rhs(ip,7) = 0.0
  end if
!END-TURB 
   
 end do
else
 do ib =ist,ien
  ip = grp%brp%nodeIndicatorArray(ib)
  normal = grp%brp%nodeNormalArray(ip,:)
  normalIncrement = sum(rhs(ip,2:4)*normal(:))
  tangentialIncrementVector(:) = rhs(ip,2:4)-normalIncrement*normal(:) 
  rhs(ip,2:4) = tangentialIncrementVector(:) 
!START-TURB
  if(ivd%turbulenceModel==1) then  
   rhs(ip,6) = 0.0
  end if 
  if(ivd%turbulenceModel==2) then 
   rhs(ip,6) = 0.0
   rhs(ip,7) = 0.0
  end if
if(ivd%turbulenceModel==3) then 
   rhs(ip,6) = 0.0
   rhs(ip,7) = 0.0
  end if
!END-TURB
 end do
end if

! Inviscid Wall

ist = grp%brp%nodeIndicatorRegister(-4)
ien = grp%brp%nodeIndicatorRegister(-3)

do ib =ist,ien
 ip = grp%brp%nodeIndicatorArray(ib)
 normal = grp%brp%nodeNormalArray(ip,:)
 normalIncrement = sum(rhs(ip,2:4)*normal(:))
 tangentialIncrementVector(:) = rhs(ip,2:4)-normalIncrement*normal(:) 
 rhs(ip,2:4) = tangentialIncrementVector(:) 
 !START-TURB
  if(ivd%turbulenceModel==1) then  
   rhs(ip,6) = 0.0
  end if 
  if(ivd%turbulenceModel==2) then
   rhs(ip,6) = 0.0
   rhs(ip,7) = 0.0
  end if
if(ivd%turbulenceModel==3) then 
   rhs(ip,6) = 0.0
   rhs(ip,7) = 0.0
  end if
!END-TURB
end do

! symmetry

ist = grp%brp%nodeIndicatorRegister(-18)
ien = grp%brp%nodeIndicatorRegister(-17)

do ib =ist,ien
 ip = grp%brp%nodeIndicatorArray(ib)

 normal = grp%brp%nodeNormalArray(ip,:)


 normalIncrement = sum(rhs(ip,2:4)*normal(:))
 tangentialIncrementVector(:) = rhs(ip,2:4)-normalIncrement*normal(:)
 rhs(ip,2:4) = tangentialIncrementVector(:)
end do


! relax engine inlets

! engine inflow - 1

ist = grp%brp%nodeIndicatorRegister(-12)
ien = grp%brp%nodeIndicatorRegister(-11)

do ib =ist,ien
 ip = grp%brp%nodeIndicatorArray(ib)
 rhs(ip,:) = ivd%engineBCRelaxation*rhs(ip,:)
end do

! engine inflow - 2

ist = grp%brp%nodeIndicatorRegister(-10)
ien = grp%brp%nodeIndicatorRegister(-9)

do ib =ist,ien
 ip = grp%brp%nodeIndicatorArray(ib)
 rhs(ip,:) = ivd%engineBCRelaxation*rhs(ip,:)
end do

! engine outflow - 1 

ist = grp%brp%nodeIndicatorRegister(-8)
ien = grp%brp%nodeIndicatorRegister(-7)

do ib =ist,ien
 ip = grp%brp%nodeIndicatorArray(ib)
!START-TURB
 rhs(ip,:) = 0.0
!END-TURB
end do

! engine outflow - 2 

ist = grp%brp%nodeIndicatorRegister(-6)
ien = grp%brp%nodeIndicatorRegister(-5)

do ib =ist,ien
 ip = grp%brp%nodeIndicatorArray(ib)
!START-TURB
 rhs(ip,:) = 0.0 
!END-TURB
end do

if(grp%brp%numberOfIONodes>0) then  
 ist = grp%brp%nodeIndicatorRegister(-4)
 ien = grp%brp%nodeIndicatorRegister(-3)
 do ib =ist,ien
  ip = grp%brp%nodeIndicatorArray(ib)
  rhs(ip,:) = 0.0
 end do
end if
end subroutine setBCsOnIncrementField
!-------------------------------------------------------------------------
subroutine setBCsOnCoarseIncrementField(grp,ivd,delu)
! nullifies coarse mesh increments at inner boundary

IMPLICIT NONE

type(GridSolverData) :: grp
type(InputVariablesData) :: ivd
real :: delu(:,:)

integer :: ist,ien,ib,ip

ist = grp%brp%nodeIndicatorRegister(-12)
ien = grp%brp%nodeIndicatorRegister(-11)

do ib =ist,ien
 ip = grp%brp%nodeIndicatorArray(ib)
 delu(ip,:) = 0.0
end do 

ist = grp%brp%nodeIndicatorRegister(-10)
ien = grp%brp%nodeIndicatorRegister(-9)

do ib =ist,ien
 ip = grp%brp%nodeIndicatorArray(ib)
 delu(ip,:) = 0.0
end do 
end subroutine setBCsOnCoarseIncrementField
!-------------------------------------------------------------------------
subroutine setTemperatureForAdiabaticWall(grp,dTx,dTy,dTz)
! sets zero gradient of the temperature on the walls
IMPLICIT NONE

type(GridSolverData) :: grp
real :: dTx(:),dTy(:),dTz(:)
integer :: ib,ist,ien,ip
real :: normalTemperatureGradient,vectorNorm,tx,ty,normal(3),ww
real :: tangTempGradientVector(3)
real :: dT(3)

ist = grp%brp%nodeIndicatorRegister(-20)
ien = grp%brp%nodeIndicatorRegister(-19)

do ib =ist,ien
 ip = grp%brp%nodeIndicatorArray(ib)
 normal = grp%brp%nodeNormalArray(ip,:)
 dT(1) = dTx(ip)
 dT(2) = dTy(ip)
 dT(3) = dTz(ip)
 normalTemperatureGradient = sum(dT*normal)
 tangTempGradientVector(:) = dT-normalTemperatureGradient*normal(:)
 dTx(ip) = tangTempGradientVector(1)
 dTy(ip) = tangTempGradientVector(2)
 dTz(ip) = tangTempGradientVector(3)
end do
end subroutine setTemperatureForAdiabaticWall
!-------------------------------------------------------------------------
subroutine makeViscosityTerm(grp,ivd,rhs,recalculateTurbulence)
! makes the viscosity fluxes using the wide stencil 
IMPLICIT NONE
include 'mpif.h'

type(GridSolverData) :: grp
type(InputVariablesData) :: ivd
real :: rhs(:,:)
logical :: recalculateTurbulence

integer :: i1,i2,i,j,ip1,ip2,faceIndicator,ib,ip,ist,ien,ind
real :: r1,r2,u1,u2,v1,v2,t1,t2,wx,wy,wz,dux,dvx,dtx,duy,dvy,dty,oneOverReynoldsNumber
real :: dwx,dwy,dwz,dwdx,dwdy,dwdz,dTdz,tz,t61,t71,t81,t91,t62,t72,t82,t92,duz,dvz,dtz
real :: w1,w2,dudz,dvdz,u31,u32,dt6,dt7,dt8,dt9,dut3,f5,dut31,dut32
real :: mu,k,T,T0,f2,f3,f4,rx,ry,vt,t11,t21,t31,t41,t51,t12,t22,t32
real :: t42,t52,dudx,dvdx,dTdx,dudy,dvdy,dTdy,dt1,dt2,dt3,dt4,dt5,oneOverPrandtlNumber
real :: twoThirds,fourThirds,gammaOverGammaMinusOne,u11,u12,u21,u22,dut1,dut2
real :: wallTangent(3),wallNormal(3),tempConv,dut11,dut12,dut21,dut22,inflowTemp
real :: adTemp1,adTemp2,turbDiffusionCoefficient
real :: viscosityCoefficient,diffusionCoefficient,oneOverTurbPrandtlNumber
real :: uu1,uu2,duux,duuy,duuz,locu,dt,increment,sign
!START-TURB
real :: onehalf,fournineth
integer :: noutput
real :: nut,a1
!END-TURB
double precision :: turbViscosityCoefficient,Ksi

 integer :: kk,istart,ifinish

 integer :: mind
 real :: mres

  grp%nodeHelpArray(1:grp%numberOfNodes,1:12) = 0.0

  gammaOverGammaMinusOne = ivd%gamma/(ivd%gamma-1.0)
  tempConv = (ivd%gamma-1.0)*ivd%MachNumber**2
  inflowTemp = ivd%inflowTemperature
!START-TURB
a1 = 0.31 
!END-TURB

#ifdef PARALLEL

  do kk=1,2
    if(kk==1) then
      istart = 1
      ifinish = grp%pdp%numberOfComSides
    else
      istart = grp%pdp%numberOfComSides + 1
      ifinish = grp%numberOfSides
    end if


    do i=istart,ifinish

#else PARALLEL

    do i=1,grp%numberOfSides

#endif PARALLEL

! indexes of nodes in side
      i1 = grp%sideIndexArray(i,1)
      i2 = grp%sideIndexArray(i,2)

! side weights from preprocessor 
      wx = grp%sideWeightsArray(i,1)  
      wy = grp%sideWeightsArray(i,2)  
      wz = grp%sideWeightsArray(i,3)  

      r1   = grp%u(i1,1)                                  ! density
      u1   = grp%u(i1,2)/r1                               ! x - velocity
      v1   = grp%u(i1,3)/r1                               ! y - velocity
      w1   = grp%u(i1,4)/r1                               ! z - velocity  
      t1 = gammaOverGammaMinusOne*grp%p(i1)/r1
      r2   = grp%u(i2,1)
      u2   = grp%u(i2,2)/r2
      v2   = grp%u(i2,3)/r2
      w2   = grp%u(i2,4)/r2
      t2 = gammaOverGammaMinusOne*grp%p(i2)/r2
      dux   = (u1+u2)*wx
      dvx   = (v1+v2)*wx
      dwx   = (w1+w2)*wx
      dtx   = (t1+t2)*wx
      duy   = (u1+u2)*wy
      dvy   = (v1+v2)*wy
      dwy   = (w1+w2)*wy
      dty   = (t1+t2)*wy
      duz   = (u1+u2)*wz
      dvz   = (v1+v2)*wz
      dwz   = (w1+w2)*wz
      dtz   = (t1+t2)*wz

      if(i1.le.grp%brp%numberOfBoundaryNodes.or.i2.le.grp%brp%numberOfBoundaryNodes) then 
        uu1 = sqrt(u1*u1 + v1*v1 + w1*w1)
        uu2 = sqrt(u2*u2 + v2*v2 + w2*w2)

        duux = (uu1+uu2)*wx
        duuy = (uu1+uu2)*wy
        duuz = (uu1+uu2)*wz
  
        grp%nodeHelpArray(i1,15) = grp%nodeHelpArray(i1,15) + duux
        grp%nodeHelpArray(i1,16) = grp%nodeHelpArray(i1,16) + duuy
        grp%nodeHelpArray(i1,17) = grp%nodeHelpArray(i1,17) + duuz

        grp%nodeHelpArray(i2,15) = grp%nodeHelpArray(i2,15) - duux
        grp%nodeHelpArray(i2,16) = grp%nodeHelpArray(i2,16) - duuy
        grp%nodeHelpArray(i2,17) = grp%nodeHelpArray(i2,17) - duuz
      end if 

      grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + dux
      grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + dvx
      grp%nodeHelpArray(i1,3) = grp%nodeHelpArray(i1,3) + dwx
      grp%nodeHelpArray(i1,4) = grp%nodeHelpArray(i1,4) + dtx
      grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) + duy
      grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) + dvy
      grp%nodeHelpArray(i1,7) = grp%nodeHelpArray(i1,7) + dwy
      grp%nodeHelpArray(i1,8) = grp%nodeHelpArray(i1,8) + dty
      grp%nodeHelpArray(i1,9) = grp%nodeHelpArray(i1,9) + duz
      grp%nodeHelpArray(i1,10) = grp%nodeHelpArray(i1,10) + dvz
      grp%nodeHelpArray(i1,11) = grp%nodeHelpArray(i1,11) + dwz
      grp%nodeHelpArray(i1,12) = grp%nodeHelpArray(i1,12) + dtz

      grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) - dux
      grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) - dvx
      grp%nodeHelpArray(i2,3) = grp%nodeHelpArray(i2,3) - dwx
      grp%nodeHelpArray(i2,4) = grp%nodeHelpArray(i2,4) - dtx
      grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) - duy
      grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) - dvy
      grp%nodeHelpArray(i2,7) = grp%nodeHelpArray(i2,7) - dwy
      grp%nodeHelpArray(i2,8) = grp%nodeHelpArray(i2,8) - dty
      grp%nodeHelpArray(i2,9) = grp%nodeHelpArray(i2,9) - duz
      grp%nodeHelpArray(i2,10) = grp%nodeHelpArray(i2,10) - dvz
      grp%nodeHelpArray(i2,11) = grp%nodeHelpArray(i2,11) - dwz
      grp%nodeHelpArray(i2,12) = grp%nodeHelpArray(i2,12) - dtz
    end do

#ifdef PARALLEL
    if(kk==1) then
      istart = 1
      ifinish = grp%pdp%numberOfBoundaryComSides
    else
      istart = grp%pdp%numberOfBoundaryComSides + 1
      ifinish = grp%brp%numberOfBoundarySides
    end if
    do i=istart,ifinish

#else PARALLEL

    do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL

      i1 = grp%brp%sideIndexArray(i,1)
      i2 = grp%brp%sideIndexArray(i,2)
      wx   = grp%brp%sideWeightsArray(i,1)
      wy   = grp%brp%sideWeightsArray(i,2) 
      wz   = grp%brp%sideWeightsArray(i,3) 

      r1   = grp%u(i1,1)
      u1   = grp%u(i1,2)/r1
      v1   = grp%u(i1,3)/r1
      w1   = grp%u(i1,4)/r1
      t1 = gammaOverGammaMinusOne*grp%p(i1)/r1
      r2   = grp%u(i2,1)
      u2   = grp%u(i2,2)/r2
      v2   = grp%u(i2,3)/r2
      w2   = grp%u(i2,4)/r2
      t2 = gammaOverGammaMinusOne*grp%p(i2)/r2
      dux   = u1+u2
      dvx   = v1+v2
      dwx   = w1+w2
      dtx   = t1+t2
      duy   = u1+u2
      dvy   = v1+v2
      dwy   = w1+w2
      dty   = t1+t2
      duz   = u1+u2
      dvz   = v1+v2
      dwz   = w1+w2
      dtz   = t1+t2

! adding boundary face contributions        

      grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + (2.*u1+dux)*wx
      grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + (2.*v1+dvx)*wx
      grp%nodeHelpArray(i1,3) = grp%nodeHelpArray(i1,3) + (2.*w1+dwx)*wx
      grp%nodeHelpArray(i1,4) = grp%nodeHelpArray(i1,4) + (2.*t1+dtx)*wx
      grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) + (2.*u1+duy)*wy
      grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) + (2.*v1+dvy)*wy
      grp%nodeHelpArray(i1,7) = grp%nodeHelpArray(i1,7) + (2.*w1+dwy)*wy
      grp%nodeHelpArray(i1,8) = grp%nodeHelpArray(i1,8) + (2.*t1+dty)*wy
      grp%nodeHelpArray(i1,9) = grp%nodeHelpArray(i1,9) + (2.*u1+duz)*wz
      grp%nodeHelpArray(i1,10) = grp%nodeHelpArray(i1,10) + (2.*v1+dvz)*wz
      grp%nodeHelpArray(i1,11) = grp%nodeHelpArray(i1,11) + (2.*w1+dwz)*wz
      grp%nodeHelpArray(i1,12) = grp%nodeHelpArray(i1,12) + (2.*t1+dtz)*wz
        
      grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + (2.*u2+dux)*wx
      grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) + (2.*v2+dvx)*wx
      grp%nodeHelpArray(i2,3) = grp%nodeHelpArray(i2,3) + (2.*w2+dwx)*wx
      grp%nodeHelpArray(i2,4) = grp%nodeHelpArray(i2,4) + (2.*t2+dtx)*wx
      grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) + (2.*u2+duy)*wy
      grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) + (2.*v2+dvy)*wy
      grp%nodeHelpArray(i2,7) = grp%nodeHelpArray(i2,7) + (2.*w2+dwy)*wy
      grp%nodeHelpArray(i2,8) = grp%nodeHelpArray(i2,8) + (2.*t2+dty)*wy
      grp%nodeHelpArray(i2,9) = grp%nodeHelpArray(i2,9) + (2.*u2+duz)*wz
      grp%nodeHelpArray(i2,10) = grp%nodeHelpArray(i2,10) + (2.*v2+dvz)*wz
      grp%nodeHelpArray(i2,11) = grp%nodeHelpArray(i2,11) + (2.*w2+dwz)*wz
      grp%nodeHelpArray(i2,12) = grp%nodeHelpArray(i2,12) + (2.*t2+dtz)*wz
    end do

#ifdef PARALLEL

    if (kk.eq.1) then
      do i=1,grp%pdp%numberOfProcesses
        if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
          call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
          position = 1
          do j=1,12
            call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                               grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
          end do
          call mp_sendv(grp%pdp%processorIDs,i,29,grp%pdp%sendFlag)
        end if
      end do
    end if
  end do

  call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,29,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

  do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
        do j=1,12
          call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                                 grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
        end do
    end if
  end do

  do i=1,grp%pdp%numberOfProcesses
    if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
      call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
      position = 1
      do j=1,12
        call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                          grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
      end do
      call mp_sendv(grp%pdp%processorIDs,i,30,grp%pdp%sendFlag)
    end if
  end do

  call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,30,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1333 continue

  do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
        do j = 1,12
          call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                             grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
        end do
    end if
  end do

#endif PARALLEL

  if(grp%gridNumber==1) then
    do i=1,grp%brp%numberOfBoundaryNodes
      vt = 1./grp%nodeVolume(i)
      wallNormal = grp%brp%nodeNormalArray(i,:)
      grp%wallStress(i,4) = vt*sum(grp%nodeHelpArray(i,15:17)*wallNormal(:))
    end do
  end if

!START-TURB
noutput=0
onehalf = 1./2.
fournineth = 4./9.
!END-TURB

  twoThirds = 2./3.
  fourThirds = 4./3.
  T0 = 1./((ivd%gamma-1)*ivd%MachNumber**2) 
  oneOverReynoldsNumber = 1./ivd%ReynoldsNumber
  oneOverPrandtlNumber = 1./ivd%PrandtlNumber
  oneOverTurbPrandtlNumber = 1./ivd%turbulentPrandtlNumber
! divide by volume to remove mass matrix on LHS
  do i=1,grp%numberOfNodes
    vt = 1./grp%nodeVolume(i)
    dudx = grp%nodeHelpArray(i,1)*vt
    dvdx = grp%nodeHelpArray(i,2)*vt
    dwdx = grp%nodeHelpArray(i,3)*vt
    dTdx = grp%nodeHelpArray(i,4)*vt
    dudy = grp%nodeHelpArray(i,5)*vt
    dvdy = grp%nodeHelpArray(i,6)*vt
    dwdy = grp%nodeHelpArray(i,7)*vt
    dTdy = grp%nodeHelpArray(i,8)*vt
    dudz = grp%nodeHelpArray(i,9)*vt
    dvdz = grp%nodeHelpArray(i,10)*vt
    dwdz = grp%nodeHelpArray(i,11)*vt
    dTdz = grp%nodeHelpArray(i,12)*vt

!START-TURB
!Saving the first derivative of velocity
  if(ivd%turbulenceModel==3 .and. ivd%SAS) then
    grp%D1U(1,i) = dudx
    grp%D1U(2,i) = dudy
    grp%D1U(3,i) = dudz
    grp%D1U(4,i) = dvdx
    grp%D1U(5,i) = dvdy
    grp%D1U(6,i) = dvdz
    grp%D1U(7,i) = dwdx
    grp%D1U(8,i) = dwdy
    grp%D1U(9,i) = dwdz
  end if
!END-TURB

! create the viscosity matrix  
    T = gammaOverGammaMinusOne*grp%p(i)/grp%u(i,1)
    if(T.le.0.001) then 
 !     write(*,*) "T negative ",T," at node ",i,grp%u(i,1),grp%p(i)
      T = 0.001
    end if
 ! Sunderlands law of viscosity 
    mu = oneOverReynoldsNumber*(((tempConv*T)**1.5)*((inflowTemp+198.6)/(tempConv*inflowTemp*T+198.6))) 
!   mu = oneOverReynoldsNumber
    k = oneOverPrandtlNumber*mu
    grp%laminarViscosity(i) = mu 
    grp%vorticity(i) = sqrt((dudy-dvdx)**2 + (dvdz-dwdy)**2 + (dwdx-dudz)**2)  ! vorticity array
    grp%divergence(i) = dudx+dvdy+dwdz

!START-TURB
 if(ivd%turbulenceModel==3 .and. ivd%SAS) then
    grp%sas(1,i)= (dudy-dvdx) + (dvdz-dwdy) + (dwdx-dudz) !instantaneous vorticity
 end if
    grp%msrt(i) = onehalf*(fournineth*(dudx**2+dvdy**2+dwdz**2)+0.5*((dudy+dvdx)**2+(dudz+dwdx)**2+(dvdz+dwdy)**2)) 
    grp%imsr(i) = sqrt((2*dudx**2+(dudy+dvdx)*dudy+(dudz+dwdx)*dudz+(dvdx+dudy)*dvdx+2*dvdy**2+(dvdz+dwdy)*dvdz+(dwdx+dudz)*dwdx&
                  +(dwdy+dvdz)*dwdy+2*dwdz**2)-twoThirds*(dudx+dvdy+dwdz)**2) 
!END-TURB

 ! add turbulence effects
    if(ivd%turbulenceModel==1) then ! Spalart-Allmaras

      turbViscosityCoefficient = oneOverReynoldsNumber*grp%u(i,1)*grp%u(i,6)

      Ksi = turbViscosityCoefficient/mu
      turbViscosityCoefficient =&
              turbViscosityCoefficient*(Ksi**3)/(Ksi**3 + 357.9)

      mu = mu + turbViscosityCoefficient
      k = k + turbViscosityCoefficient*oneOverTurbPrandtlNumber 
    end if

!START-TURB
     if(ivd%turbulenceModel==2) then ! k-omega model
       turbViscosityCoefficient = grp%u(i,1)*grp%u(i,6)/grp%u(i,7)
!     Limiting turbulence   
       turbViscosityCoefficient=max(0.0,turbViscosityCoefficient)
       if(abs(turbViscosityCoefficient*ivd%ReynoldsNumber) .gt. ivd%maxTurbulenceValueMU/a1)then
       if (noutput==0.and.grp%gridNumber==1)write(*,*)'WARNING, high mut',turbViscosityCoefficient*ivd%ReynoldsNumber,grp%coordinates(i,:)
       noutput=1
       turbViscosityCoefficient=ivd%maxTurbulenceValueMU*oneOverReynoldsNumber
      end if     
      mu = mu + turbViscosityCoefficient
      k = k + turbViscosityCoefficient*oneOverTurbPrandtlNumber 
     end if

     if(ivd%turbulenceModel==3) then !SST
 !    Eddy viscosity limiter
       nut=a1*grp%u(i,6)/(max(a1*grp%u(i,7),grp%imsr(i)*grp%F2(i)))
       turbViscosityCoefficient = grp%u(i,1)*nut
       grp%TVR(i)= turbViscosityCoefficient/grp%laminarViscosity(i)
 !    Limiting turbulence
       turbViscosityCoefficient=max(0.0,turbViscosityCoefficient)
       if(abs(turbViscosityCoefficient*ivd%ReynoldsNumber) .gt. ivd%maxTurbulenceValueMU/a1)then
       if (noutput==0.and.grp%gridNumber==1)write(*,*)'WARNING, high mut',turbViscosityCoefficient*ivd%ReynoldsNumber,grp%coordinates(i,:)
       noutput=1
       turbViscosityCoefficient=ivd%maxTurbulenceValueMU*oneOverReynoldsNumber
       end if
       mu = mu + turbViscosityCoefficient
       k = k + turbViscosityCoefficient*oneOverTurbPrandtlNumber 
     end if
!END-TURB

! stress tensor
    grp%nodeHelpArray(i,1) = mu*(fourThirds*dudx-twoThirds*(dvdy+dwdz)) 
    grp%nodeHelpArray(i,2) = mu*(dudy+dvdx) 
    grp%nodeHelpArray(i,3) = mu*(dudz+dwdx) 
    grp%nodeHelpArray(i,4) = mu*(fourThirds*dvdy-twoThirds*(dudx+dwdz)) 
    grp%nodeHelpArray(i,5) = mu*(dvdz+dwdy) 
    grp%nodeHelpArray(i,6) = mu*(fourThirds*dwdz-twoThirds*(dudx+dvdy)) 
    grp%nodeHelpArray(i,7) = k*dTdx 
    grp%nodeHelpArray(i,8) = k*dTdy 
    grp%nodeHelpArray(i,9) = k*dTdz
  end do

! for adiabatic wall, set temperature gradient equal to zero
  if(grp%gridNumber==1) then
    call setTemperatureForAdiabaticWall(grp,grp%nodeHelpArray(1:grp%numberOfNodes,7),grp%nodeHelpArray(1:grp%numberOfNodes,8),&
                                        grp%nodeHelpArray(1:grp%numberOfNodes,9))
  end if



! make stress tensor along boundary

 
  if(grp%gridNumber==1) then 
    do ib=1,grp%brp%numberOfBoundaryNodes
      ip = grp%brp%nodeIndicatorArray(ib)
      wallNormal = grp%brp%nodeNormalArray(ip,:)

      grp%wallStress(ip,1) = (grp%nodeHelpArray(ip,1)*wallNormal(1)&
                           +  grp%nodeHelpArray(ip,2)*wallNormal(2)&
                           +  grp%nodeHelpArray(ip,3)*wallNormal(3))
      grp%wallStress(ip,2) = (grp%nodeHelpArray(ip,2)*wallNormal(1)&
                           +  grp%nodeHelpArray(ip,4)*wallNormal(2)&
                           +  grp%nodeHelpArray(ip,5)*wallNormal(3))
      grp%wallStress(ip,3) = (grp%nodeHelpArray(ip,3)*wallNormal(1)&
                           +  grp%nodeHelpArray(ip,5)*wallNormal(2)&
                           +  grp%nodeHelpArray(ip,6)*wallNormal(3))
      grp%wallStress(ip,4) = 0.0
    end do
  end if

  grp%nodeHelpArray(1:grp%numberOfNodes,10:13) = 0.0


#ifdef PARALLEL

  do kk=1,2
    if(kk==1) then
      istart = 1
      ifinish = grp%pdp%numberOfComSides
    else
      istart = grp%pdp%numberOfComSides + 1
      ifinish = grp%numberOfSides
    end if
 
 
    do i=istart,ifinish

#else PARALLEL

    do i=1,grp%numberOfSides

#endif PARALLEL

      i1 = grp%sideIndexArray(i,1)
      i2 = grp%sideIndexArray(i,2)

! side weights from preprocessor 
      wx = grp%sideWeightsArray(i,1)
      wy = grp%sideWeightsArray(i,2)  
      wz = grp%sideWeightsArray(i,3)  

      t11 = grp%nodeHelpArray(i1,1)
      t21 = grp%nodeHelpArray(i1,2)
      t31 = grp%nodeHelpArray(i1,3)
      t41 = grp%nodeHelpArray(i1,4)
      t51 = grp%nodeHelpArray(i1,5)
      t61 = grp%nodeHelpArray(i1,6)
      t71 = grp%nodeHelpArray(i1,7)
      t81 = grp%nodeHelpArray(i1,8)
      t91 = grp%nodeHelpArray(i1,9)
      u11 = grp%u(i1,2)/grp%u(i1,1)
      u21 = grp%u(i1,3)/grp%u(i1,1) 
      u31 = grp%u(i1,4)/grp%u(i1,1) 

      t12 = grp%nodeHelpArray(i2,1)
      t22 = grp%nodeHelpArray(i2,2)
      t32 = grp%nodeHelpArray(i2,3)
      t42 = grp%nodeHelpArray(i2,4)
      t52 = grp%nodeHelpArray(i2,5)
      t62 = grp%nodeHelpArray(i2,6)
      t72 = grp%nodeHelpArray(i2,7)
      t82 = grp%nodeHelpArray(i2,8)
      t92 = grp%nodeHelpArray(i2,9)
      u12 = grp%u(i2,2)/grp%u(i2,1)
      u22 = grp%u(i2,3)/grp%u(i2,1)
      u32 = grp%u(i2,4)/grp%u(i2,1)

      dt1 = t11 + t12 
      dt2 = t21 + t22 
      dt3 = t31 + t32 
      dt4 = t41 + t42
      dt5 = t51 + t52
      dt6 = t61 + t62
      dt7 = t71 + t72
      dt8 = t81 + t82
      dt9 = t91 + t92
      dut1 = u11*t11 + u12*t12 + u21*t21 + u22*t22 + u31*t31 + u32*t32 ! dissipation 
      dut2 = u11*t21 + u12*t22 + u21*t41 + u22*t42 + u31*t51 + u32*t52 
      dut3 = u11*t31 + u12*t32 + u21*t51 + u22*t52 + u31*t61 + u32*t62 


      f2 = wx*dt1 + wy*dt2 + wz*dt3
      f3 = wx*dt2 + wy*dt4 + wz*dt5
      f4 = wx*dt3 + wy*dt5 + wz*dt6
      f5 = wx*(dt7+dut1) + wy*(dt8+dut2) + wz*(dt9+dut3)

      grp%nodeHelpArray(i1,10) = grp%nodeHelpArray(i1,10) + f2 
      grp%nodeHelpArray(i1,11) = grp%nodeHelpArray(i1,11) + f3  
      grp%nodeHelpArray(i1,12) = grp%nodeHelpArray(i1,12) + f4 
      grp%nodeHelpArray(i1,13) = grp%nodeHelpArray(i1,13) + f5

      grp%nodeHelpArray(i2,10) = grp%nodeHelpArray(i2,10) - f2 
      grp%nodeHelpArray(i2,11) = grp%nodeHelpArray(i2,11) - f3  
      grp%nodeHelpArray(i2,12) = grp%nodeHelpArray(i2,12) - f4 
      grp%nodeHelpArray(i2,13) = grp%nodeHelpArray(i2,13) - f5

    end do

#ifdef PARALLEL
    if(kk==1) then
      istart = 1
      ifinish = grp%pdp%numberOfBoundaryComSides
    else
      istart = grp%pdp%numberOfBoundaryComSides + 1
      ifinish = grp%brp%numberOfBoundarySides
    end if
    do i=istart,ifinish

#else PARALLEL

    do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL

      faceIndicator = grp%brp%sideIndexArray(i,3)
      if(faceIndicator.ne.2.and.faceIndicator.ne.9) then 
        i1 = grp%brp%sideIndexArray(i,1)
        i2 = grp%brp%sideIndexArray(i,2)
        wx   = grp%brp%sideWeightsArray(i,1)
        wy   = grp%brp%sideWeightsArray(i,2) 
        wz   = grp%brp%sideWeightsArray(i,3) 
 
        t11 = grp%nodeHelpArray(i1,1)
        t21 = grp%nodeHelpArray(i1,2)
        t31 = grp%nodeHelpArray(i1,3)
        t41 = grp%nodeHelpArray(i1,4)
        t51 = grp%nodeHelpArray(i1,5)
        t61 = grp%nodeHelpArray(i1,6)
        t71 = grp%nodeHelpArray(i1,7)
        t81 = grp%nodeHelpArray(i1,8)
        t91 = grp%nodeHelpArray(i1,9)
        u11 = grp%u(i1,2)/grp%u(i1,1)
        u21 = grp%u(i1,3)/grp%u(i1,1)
        u31 = grp%u(i1,4)/grp%u(i1,1)

        t12 = grp%nodeHelpArray(i2,1)
        t22 = grp%nodeHelpArray(i2,2)
        t32 = grp%nodeHelpArray(i2,3)
        t42 = grp%nodeHelpArray(i2,4)
        t52 = grp%nodeHelpArray(i2,5)
        t62 = grp%nodeHelpArray(i2,6)
        t72 = grp%nodeHelpArray(i2,7)
        t82 = grp%nodeHelpArray(i2,8)
        t92 = grp%nodeHelpArray(i2,9)
        u12 = grp%u(i2,2)/grp%u(i2,1)
        u22 = grp%u(i2,3)/grp%u(i2,1)
        u32 = grp%u(i2,4)/grp%u(i2,1)

        dt1 = t11 + t12
        dt2 = t21 + t22
        dt3 = t31 + t32
        dt4 = t41 + t42
        dt5 = t51 + t52
        dt6 = t61 + t62
        dt7 = t71 + t72
        dt8 = t81 + t82
        dt9 = t91 + t92

        dut11 = u11*t11 + u21*t21 + u31*t31
        dut12 = u12*t12 + u22*t22 + u32*t32
        dut21 = u11*t21 + u21*t41 + u31*t51 
        dut22 = u12*t22 + u22*t42 + u32*t52 
        dut31 = u11*t31 + u21*t51 + u31*t61 
        dut32 = u12*t32 + u22*t52 + u32*t62 

        dut1 = dut11 + dut12
        dut2 = dut21 + dut22
        dut3 = dut31 + dut32

        f2 = wx*dt1 + wy*dt2 + wz*dt3
        f3 = wx*dt2 + wy*dt4 + wz*dt5
        f4 = wx*dt3 + wy*dt5 + wz*dt6
        f5 = wx*(dt7+dut1) + wy*(dt8+dut2) + wz*(dt9+dut3)

        grp%nodeHelpArray(i1,10) = grp%nodeHelpArray(i1,10) + f2 + 2.0*(wx*t11+wy*t21+wz*t31) 
        grp%nodeHelpArray(i1,11) = grp%nodeHelpArray(i1,11) + f3 + 2.0*(wx*t21+wy*t41+wz*t51) 
        grp%nodeHelpArray(i1,12) = grp%nodeHelpArray(i1,12) + f4 + 2.0*(wx*t31+wy*t51+wz*t61) 
        grp%nodeHelpArray(i1,13) = grp%nodeHelpArray(i1,13) + f5 + 2.0*(wx*(t71+dut11) + wy*(t81+dut21)+wz*(t91+dut31))

        grp%nodeHelpArray(i2,10) = grp%nodeHelpArray(i2,10) + f2 + 2.0*(wx*t12+wy*t22+wz*t32) 
        grp%nodeHelpArray(i2,11) = grp%nodeHelpArray(i2,11) + f3 + 2.0*(wx*t22+wy*t42+wz*t52) 
        grp%nodeHelpArray(i2,12) = grp%nodeHelpArray(i2,12) + f4 + 2.0*(wx*t32+wy*t52+wz*t62) 
        grp%nodeHelpArray(i2,13) = grp%nodeHelpArray(i2,13) + f5 + 2.0*(wx*(t72+dut12) + wy*(t82+dut22)+wz*(t92+dut32)) 
      end if
    end do

#ifdef PARALLEL

    if (kk.eq.1) then
      do i=1,grp%pdp%numberOfProcesses
        if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
          call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
          position = 1
          do j=10,13
            call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                               grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
          end do
          call mp_sendv(grp%pdp%processorIDs,i,31,grp%pdp%sendFlag)
        end if
      end do
    end if
  end do

  call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,31,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )


  do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
        do j=10,13
          call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                                 grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
        end do
    end if
  end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

  do i=1,grp%pdp%numberOfProcesses
    if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
      call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
      position = 1
      do j=10,13
        call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                          grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
      end do
      call mp_sendv(grp%pdp%processorIDs,i,32,grp%pdp%sendFlag)
    end if
  end do

  call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,32,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

  do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
        do j =10,13
          call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                             grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
        end do
    end if
  end do

#endif PARALLEL

! add viscosity to dissipation vector

  do i=1,grp%numberOfNodes
    rhs(i,2:5) = rhs(i,2:5) - grp%nodeHelpArray(i,10:13)
  end do
end subroutine makeViscosityTerm

!-------------------------------------------------------------------------
subroutine makeViscosityTerm2(grp,ivd,rhs,recalculateTurbulence)
! makes the viscosity fluxes using the compact stencil
IMPLICIT NONE
include 'mpif.h'

type(GridSolverData) :: grp
type(InputVariablesData) :: ivd
real :: rhs(:,:)
logical :: recalculateTurbulence

integer :: i1,i2,i,j,ip1,ip2,faceIndicator,ib,ip,ist,ien,ind
real :: r1,r2,u1,u2,v1,v2,t1,t2,wx,wy,wz,dux,dvx,dtx,duy,dvy,dty,oneOverReynoldsNumber
real :: dwx,dwy,dwz,dwdx,dwdy,dwdz,dTdz,tz,t61,t71,t81,t91,t62,t72,t82,t92,duz,dvz,dtz
real :: w1,w2,dudz,dvdz,u31,u32,dt6,dt7,dt8,dt9,dut3,f5,dut31,dut32
real :: mu,k,T,T0,f2,f3,f4,rx,ry,vt,t11,t21,t31,t41,t51,t12,t22,t32
real :: t42,t52,dudx,dvdx,dTdx,dudy,dvdy,dTdy,dt1,dt2,dt3,dt4,dt5,oneOverPrandtlNumber
real :: twoThirds,fourThirds,gammaOverGammaMinusOne,u11,u12,u21,u22,dut1,dut2
real :: wallTangent(3),wallNormal(3),tempConv,dut11,dut12,dut21,dut22,inflowTemp
real :: adTemp1,adTemp2,turbDiffusionCoefficient
real :: viscosityCoefficient,diffusionCoefficient,oneOverTurbPrandtlNumber
real :: uu1,uu2,duux,duuy,duuz,locu,dt,increment,sign

!START-TURB
real :: onehalf,fournineth
integer :: noutput
real :: nut,a1
!END-TURB

double precision :: turbViscosityCoefficient,Ksi

real :: c1,c2,c3,c4,c5,c6,c7,c8,c9,mmu1,mmu2,k1,k2,rhoinv1,rhoinv2
real :: r(3),dudx1,dudy1,dudz1,dudx2,dudy2,dudz2
real :: dvdx1,dvdy1,dvdz1,dvdx2,dvdy2,dvdz2,dwdx1,dwdy1,dwdz1,dwdx2,dwdy2,dwdz2
real :: dTdx1,dTdy1,dTdz1,dTdx2,dTdy2,dTdz2,dudr,dvdr,dwdr,dTdr,inverseSideLength
real :: du,dv,dw

 integer :: kk,istart,ifinish

  grp%nodeHelpArray(1:grp%numberOfNodes,1:14) = 0.0

  gammaOverGammaMinusOne = ivd%gamma/(ivd%gamma-1.0)
  tempConv = (ivd%gamma-1.0)*ivd%MachNumber**2
  inflowTemp = ivd%inflowTemperature

!START-TURB
  a1 = 0.31
!END-TURB
 
#ifdef PARALLEL

  do kk=1,2
    if(kk==1) then
      istart = 1
      ifinish = grp%pdp%numberOfComSides
    else
      istart = grp%pdp%numberOfComSides + 1
      ifinish = grp%numberOfSides
    end if


    do i=istart,ifinish

#else PARALLEL

    do i=1,grp%numberOfSides

#endif PARALLEL


! indexes of nodes in side
      i1 = grp%sideIndexArray(i,1)
      i2 = grp%sideIndexArray(i,2)

! side weights from preprocessor 
      wx = grp%sideWeightsArray(i,1)  
      wy = grp%sideWeightsArray(i,2)  
      wz = grp%sideWeightsArray(i,3)  

      r1   = grp%u(i1,1)                                  ! density
      u1   = grp%u(i1,2)/r1                               ! x - velocity
      v1   = grp%u(i1,3)/r1                               ! y - velocity
      w1   = grp%u(i1,4)/r1                               ! z - velocity  
      t1 = gammaOverGammaMinusOne*grp%p(i1)/r1
      r2   = grp%u(i2,1)
      u2   = grp%u(i2,2)/r2
      v2   = grp%u(i2,3)/r2
      w2   = grp%u(i2,4)/r2
      t2 = gammaOverGammaMinusOne*grp%p(i2)/r2
      dux   = (u1+u2)*wx
      dvx   = (v1+v2)*wx
      dwx   = (w1+w2)*wx
      dtx   = (t1+t2)*wx
      duy   = (u1+u2)*wy
      dvy   = (v1+v2)*wy
      dwy   = (w1+w2)*wy
      dty   = (t1+t2)*wy
      duz   = (u1+u2)*wz
      dvz   = (v1+v2)*wz
      dwz   = (w1+w2)*wz
      dtz   = (t1+t2)*wz

      if(i1.le.grp%brp%numberOfBoundaryNodes.or.i2.le.grp%brp%numberOfBoundaryNodes) then 
        uu1 = sqrt(u1*u1 + v1*v1 + w1*w1)
        uu2 = sqrt(u2*u2 + v2*v2 + w2*w2)

        duux = (uu1+uu2)*wx
        duuy = (uu1+uu2)*wy
        duuz = (uu1+uu2)*wz
  
        grp%nodeHelpArray(i1,15) = grp%nodeHelpArray(i1,15) + duux
        grp%nodeHelpArray(i1,16) = grp%nodeHelpArray(i1,16) + duuy
        grp%nodeHelpArray(i1,17) = grp%nodeHelpArray(i1,17) + duuz

        grp%nodeHelpArray(i2,15) = grp%nodeHelpArray(i2,15) - duux
        grp%nodeHelpArray(i2,16) = grp%nodeHelpArray(i2,16) - duuy
        grp%nodeHelpArray(i2,17) = grp%nodeHelpArray(i2,17) - duuz
      end if 

      grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + dux
      grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + dvx
      grp%nodeHelpArray(i1,3) = grp%nodeHelpArray(i1,3) + dwx
      grp%nodeHelpArray(i1,4) = grp%nodeHelpArray(i1,4) + dtx
      grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) + duy
      grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) + dvy
      grp%nodeHelpArray(i1,7) = grp%nodeHelpArray(i1,7) + dwy
      grp%nodeHelpArray(i1,8) = grp%nodeHelpArray(i1,8) + dty
      grp%nodeHelpArray(i1,9) = grp%nodeHelpArray(i1,9) + duz
      grp%nodeHelpArray(i1,10) = grp%nodeHelpArray(i1,10) + dvz
      grp%nodeHelpArray(i1,11) = grp%nodeHelpArray(i1,11) + dwz
      grp%nodeHelpArray(i1,12) = grp%nodeHelpArray(i1,12) + dtz

      grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) - dux
      grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) - dvx
      grp%nodeHelpArray(i2,3) = grp%nodeHelpArray(i2,3) - dwx
      grp%nodeHelpArray(i2,4) = grp%nodeHelpArray(i2,4) - dtx
      grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) - duy
      grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) - dvy
      grp%nodeHelpArray(i2,7) = grp%nodeHelpArray(i2,7) - dwy
      grp%nodeHelpArray(i2,8) = grp%nodeHelpArray(i2,8) - dty
      grp%nodeHelpArray(i2,9) = grp%nodeHelpArray(i2,9) - duz
      grp%nodeHelpArray(i2,10) = grp%nodeHelpArray(i2,10) - dvz
      grp%nodeHelpArray(i2,11) = grp%nodeHelpArray(i2,11) - dwz
      grp%nodeHelpArray(i2,12) = grp%nodeHelpArray(i2,12) - dtz
    end do

#ifdef PARALLEL
    if(kk==1) then
      istart = 1
      ifinish = grp%pdp%numberOfBoundaryComSides
    else
      istart = grp%pdp%numberOfBoundaryComSides + 1
      ifinish = grp%brp%numberOfBoundarySides
    end if
    do i=istart,ifinish

#else PARALLEL

    do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL

      i1 = grp%brp%sideIndexArray(i,1)
      i2 = grp%brp%sideIndexArray(i,2)
      wx   = grp%brp%sideWeightsArray(i,1)
      wy   = grp%brp%sideWeightsArray(i,2) 
      wz   = grp%brp%sideWeightsArray(i,3) 

      r1   = grp%u(i1,1)
      u1   = grp%u(i1,2)/r1
      v1   = grp%u(i1,3)/r1
      w1   = grp%u(i1,4)/r1
      t1 = gammaOverGammaMinusOne*grp%p(i1)/r1
      r2   = grp%u(i2,1)
      u2   = grp%u(i2,2)/r2
      v2   = grp%u(i2,3)/r2
      w2   = grp%u(i2,4)/r2
      t2 = gammaOverGammaMinusOne*grp%p(i2)/r2
      dux   = u1+u2
      dvx   = v1+v2
      dwx   = w1+w2
      dtx   = t1+t2
      duy   = u1+u2
      dvy   = v1+v2
      dwy   = w1+w2
      dty   = t1+t2
      duz   = u1+u2
      dvz   = v1+v2
      dwz   = w1+w2
      dtz   = t1+t2

! adding boundary face contributions        

      if(ivd%boundaryTerm==1) then 
        grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + (2.*u1+dux)*wx
        grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + (2.*v1+dvx)*wx
        grp%nodeHelpArray(i1,3) = grp%nodeHelpArray(i1,3) + (2.*w1+dwx)*wx
        grp%nodeHelpArray(i1,4) = grp%nodeHelpArray(i1,4) + (2.*t1+dtx)*wx
        grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) + (2.*u1+duy)*wy
        grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) + (2.*v1+dvy)*wy
        grp%nodeHelpArray(i1,7) = grp%nodeHelpArray(i1,7) + (2.*w1+dwy)*wy
        grp%nodeHelpArray(i1,8) = grp%nodeHelpArray(i1,8) + (2.*t1+dty)*wy
        grp%nodeHelpArray(i1,9) = grp%nodeHelpArray(i1,9) + (2.*u1+duz)*wz
        grp%nodeHelpArray(i1,10) = grp%nodeHelpArray(i1,10) + (2.*v1+dvz)*wz
        grp%nodeHelpArray(i1,11) = grp%nodeHelpArray(i1,11) + (2.*w1+dwz)*wz
        grp%nodeHelpArray(i1,12) = grp%nodeHelpArray(i1,12) + (2.*t1+dtz)*wz
        
        grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + (2.*u2+dux)*wx
        grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) + (2.*v2+dvx)*wx
        grp%nodeHelpArray(i2,3) = grp%nodeHelpArray(i2,3) + (2.*w2+dwx)*wx
        grp%nodeHelpArray(i2,4) = grp%nodeHelpArray(i2,4) + (2.*t2+dtx)*wx
        grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) + (2.*u2+duy)*wy
        grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) + (2.*v2+dvy)*wy
        grp%nodeHelpArray(i2,7) = grp%nodeHelpArray(i2,7) + (2.*w2+dwy)*wy
        grp%nodeHelpArray(i2,8) = grp%nodeHelpArray(i2,8) + (2.*t2+dty)*wy
        grp%nodeHelpArray(i2,9) = grp%nodeHelpArray(i2,9) + (2.*u2+duz)*wz
        grp%nodeHelpArray(i2,10) = grp%nodeHelpArray(i2,10) + (2.*v2+dvz)*wz
        grp%nodeHelpArray(i2,11) = grp%nodeHelpArray(i2,11) + (2.*w2+dwz)*wz
        grp%nodeHelpArray(i2,12) = grp%nodeHelpArray(i2,12) + (2.*t2+dtz)*wz
      else
        grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + 4.*u1*wx
        grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + 4.*v1*wx
        grp%nodeHelpArray(i1,3) = grp%nodeHelpArray(i1,3) + 4.*w1*wx
        grp%nodeHelpArray(i1,4) = grp%nodeHelpArray(i1,4) + 4.*t1*wx
        grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) + 4.*u1*wy
        grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) + 4.*v1*wy
        grp%nodeHelpArray(i1,7) = grp%nodeHelpArray(i1,7) + 4.*w1*wy
        grp%nodeHelpArray(i1,8) = grp%nodeHelpArray(i1,8) + 4.*t1*wy
        grp%nodeHelpArray(i1,9) = grp%nodeHelpArray(i1,9) + 4.*u1*wz
        grp%nodeHelpArray(i1,10) = grp%nodeHelpArray(i1,10) + 4.*v1*wz
        grp%nodeHelpArray(i1,11) = grp%nodeHelpArray(i1,11) + 4.*w1*wz
        grp%nodeHelpArray(i1,12) = grp%nodeHelpArray(i1,12) + 4.*t1*wz
  
        grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + 4.*u2*wx
        grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) + 4.*v2*wx
        grp%nodeHelpArray(i2,3) = grp%nodeHelpArray(i2,3) + 4.*w2*wx
        grp%nodeHelpArray(i2,4) = grp%nodeHelpArray(i2,4) + 4.*t2*wx
        grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) + 4.*u2*wy
        grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) + 4.*v2*wy
        grp%nodeHelpArray(i2,7) = grp%nodeHelpArray(i2,7) + 4.*w2*wy
        grp%nodeHelpArray(i2,8) = grp%nodeHelpArray(i2,8) + 4.*t2*wy
        grp%nodeHelpArray(i2,9) = grp%nodeHelpArray(i2,9) + 4.*u2*wz
        grp%nodeHelpArray(i2,10) = grp%nodeHelpArray(i2,10) + 4.*v2*wz
        grp%nodeHelpArray(i2,11) = grp%nodeHelpArray(i2,11) + 4.*w2*wz
        grp%nodeHelpArray(i2,12) = grp%nodeHelpArray(i2,12) + 4.*t2*wz
      end if
    end do

#ifdef PARALLEL

    if (kk.eq.1) then
      do i=1,grp%pdp%numberOfProcesses
        if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
          call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
          position = 1
          do j=1,12
            call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                               grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
          end do
          call mp_sendv(grp%pdp%processorIDs,i,29,grp%pdp%sendFlag)
        end if
      end do
    end if
  end do

  call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,29,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )


  do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
        do j=1,12
          call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                                 grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
        end do
    end if
  end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

  do i=1,grp%pdp%numberOfProcesses
    if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
      call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
      position = 1
      do j=1,12
        call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                          grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
      end do
      call mp_sendv(grp%pdp%processorIDs,i,30,grp%pdp%sendFlag)
    end if
  end do

  call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,30,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

  do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
        do j = 1,12
          call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                             grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
        end do
    end if
  end do

#endif PARALLEL

!START-TURB
noutput=0
onehalf = 1./2.
fournineth = 4./9.
!END-TURB

  twoThirds = 2./3.
  fourThirds = 4./3.
  T0 = 1./((ivd%gamma-1)*ivd%MachNumber**2) 
  oneOverReynoldsNumber = 1./ivd%ReynoldsNumber
  oneOverPrandtlNumber = 1./ivd%PrandtlNumber
  oneOverTurbPrandtlNumber = 1./ivd%turbulentPrandtlNumber
  ! divide by volume to remove mass matrix on LHS
  do i=1,grp%numberOfNodes
    vt = 1./grp%nodeVolume(i)

  
    grp%nodeHelpArray(i,1:12) = grp%nodeHelpArray(i,1:12)*vt

    dudx = grp%nodeHelpArray(i,1)
    dvdx = grp%nodeHelpArray(i,2)
    dwdx = grp%nodeHelpArray(i,3)
    dTdx = grp%nodeHelpArray(i,4)
    dudy = grp%nodeHelpArray(i,5)
    dvdy = grp%nodeHelpArray(i,6)
    dwdy = grp%nodeHelpArray(i,7)
    dTdy = grp%nodeHelpArray(i,8)
    dudz = grp%nodeHelpArray(i,9)
    dvdz = grp%nodeHelpArray(i,10)
    dwdz = grp%nodeHelpArray(i,11)
    dTdz = grp%nodeHelpArray(i,12)

!START-TURB
!Saving the first derivative of velocity
  if(ivd%turbulenceModel==3 .and. ivd%SAS) then
    grp%D1U(1,i) = dudx
    grp%D1U(2,i) = dudy
    grp%D1U(3,i) = dudz
    grp%D1U(4,i) = dvdx
    grp%D1U(5,i) = dvdy
    grp%D1U(6,i) = dvdz
    grp%D1U(7,i) = dwdx
    grp%D1U(8,i) = dwdy
    grp%D1U(9,i) = dwdz
  end if
!END-TURB

! create the viscosity matrix  
    T = gammaOverGammaMinusOne*grp%p(i)/grp%u(i,1)
    if(T.le.0.001) then 
      write(*,*) "T negative ",T," at node ",i,grp%u(i,1),grp%p(i)
      T = 0.001
    end if


 ! Sunderlands law of viscosity 
    mu = oneOverReynoldsNumber*(((tempConv*T)**1.5)*((inflowTemp+198.6)/(tempConv*inflowTemp*T+198.6))) 
!   mu = oneOverReynoldsNumber
    k = oneOverPrandtlNumber*mu
    grp%laminarViscosity(i) = mu 
    grp%vorticity(i) = sqrt((dudy-dvdx)**2 + (dvdz-dwdy)**2 + (dwdx-dudz)**2)  ! vorticity array
    grp%divergence(i) = dudx+dvdy+dwdz

!START-TURB
  if(ivd%turbulenceModel==3 .and. ivd%SAS) then
    grp%sas(1,i)= (dudy-dvdx) + (dvdz-dwdy) + (dwdx-dudz) !instantaneous vorticity
  end if
    grp%msrt(i) = onehalf*(fournineth*(dudx**2+dvdy**2+dwdz**2)+0.5*((dudy+dvdx)**2+(dudz+dwdx)**2+(dvdz+dwdy)**2)) 
    grp%imsr(i) = sqrt((2*dudx**2+(dudy+dvdx)*dudy+(dudz+dwdx)*dudz+(dvdx+dudy)*dvdx+2*dvdy**2+(dvdz+dwdy)*dvdz+(dwdx+dudz)*dwdx&
                  +(dwdy+dvdz)*dwdy+2*dwdz**2)-twoThirds*(dudx+dvdy+dwdz)**2) 
!END-TURB

 ! add turbulence effects
    if(ivd%turbulenceModel==1) then ! Spalart-Allmaras

      turbViscosityCoefficient = oneOverReynoldsNumber*grp%u(i,1)*grp%u(i,6)

      Ksi = turbViscosityCoefficient/mu
      turbViscosityCoefficient =&
          turbViscosityCoefficient*(Ksi**3)/(Ksi**3 + 357.9)

      mu = mu + turbViscosityCoefficient
      k = k + turbViscosityCoefficient*oneOverTurbPrandtlNumber 
    end if

!START-TURB
     if(ivd%turbulenceModel==2) then ! k-omega model
       turbViscosityCoefficient = grp%u(i,1)*grp%u(i,6)/grp%u(i,7)
!     Limiting turbulence   
       turbViscosityCoefficient=max(0.0,turbViscosityCoefficient)
       if(abs(turbViscosityCoefficient*ivd%ReynoldsNumber) .gt. ivd%maxTurbulenceValueMU/a1)then
       if (noutput==0.and.grp%gridNumber==1)write(*,*)'WARNING, high mut',turbViscosityCoefficient*ivd%ReynoldsNumber,grp%coordinates(i,:)
       noutput=1
       turbViscosityCoefficient=ivd%maxTurbulenceValueMU*oneOverReynoldsNumber
      end if     
      mu = mu + turbViscosityCoefficient
      k = k + turbViscosityCoefficient*oneOverTurbPrandtlNumber 
     end if

     if(ivd%turbulenceModel==3) then !SST
 !    Eddy viscosity limiter
       nut=a1*grp%u(i,6)/(max(a1*grp%u(i,7),grp%imsr(i)*grp%F2(i)))
       turbViscosityCoefficient = grp%u(i,1)*nut
       grp%TVR(i)= turbViscosityCoefficient/grp%laminarViscosity(i)
 !    Limiting turbulence
       turbViscosityCoefficient=max(0.0,turbViscosityCoefficient)
       if(abs(turbViscosityCoefficient*ivd%ReynoldsNumber) .gt. ivd%maxTurbulenceValueMU/a1)then
       if (noutput==0.and.grp%gridNumber==1)write(*,*)'WARNING, high mut',turbViscosityCoefficient*ivd%ReynoldsNumber,grp%coordinates(i,:)
       noutput=1
       turbViscosityCoefficient=ivd%maxTurbulenceValueMU*oneOverReynoldsNumber
       end if
       mu = mu + turbViscosityCoefficient
       k = k + turbViscosityCoefficient*oneOverTurbPrandtlNumber 
     end if
!END-TURB


    grp%nodeHelpArray(i,1:3) = mu*grp%nodeHelpArray(i,1:3)
    grp%nodeHelpArray(i,5:7) = mu*grp%nodeHelpArray(i,5:7)
    grp%nodeHelpArray(i,9:11) = mu*grp%nodeHelpArray(i,9:11)
    grp%nodeHelpArray(i,4) = k*grp%nodeHelpArray(i,4)
    grp%nodeHelpArray(i,8) = k*grp%nodeHelpArray(i,8)
    grp%nodeHelpArray(i,12) = k*grp%nodeHelpArray(i,12)


 
    grp%nodeHelpArray(i,13) = mu
    grp%nodeHelpArray(i,14) = k
  end do

! for adiabatic wall, set temperature gradient equal to zero
  if(grp%gridNumber==1) then
    call setTemperatureForAdiabaticWall(grp,grp%nodeHelpArray(1:grp%numberOfNodes,4),grp%nodeHelpArray(1:grp%numberOfNodes,8),&
                                        grp%nodeHelpArray(1:grp%numberOfNodes,12))
  end if


! make stress tensor along boundary

 
  if(grp%gridNumber==1) then 
    ist = grp%brp%nodeIndicatorRegister(-20)
    ien = grp%brp%nodeIndicatorRegister(-19)
    do ib =ist,ien
      ip = grp%brp%nodeIndicatorArray(ib)
      wallNormal = grp%brp%nodeNormalArray(ip,:)

      dudx = grp%nodeHelpArray(ip,1)
      dvdx = grp%nodeHelpArray(ip,2)
      dwdx = grp%nodeHelpArray(ip,3)
      dTdx = grp%nodeHelpArray(ip,4)
      dudy = grp%nodeHelpArray(ip,5)
      dvdy = grp%nodeHelpArray(ip,6)
      dwdy = grp%nodeHelpArray(ip,7)
      dTdy = grp%nodeHelpArray(ip,8)
      dudz = grp%nodeHelpArray(ip,9)
      dvdz = grp%nodeHelpArray(ip,10)
      dwdz = grp%nodeHelpArray(ip,11)
      dTdz = grp%nodeHelpArray(ip,12)

      c1 = fourThirds*dudx-twoThirds*(dvdy+dwdz)
      c2 = dudy+dvdx
      c3 = dudz+dwdx
      c4 = fourThirds*dvdy-twoThirds*(dudx+dwdz)
      c5 = dvdz+dwdy
      c6 = fourThirds*dwdz-twoThirds*(dudx+dvdy)
      c7 = dTdx
      c8 = dTdy
      c9 = dTdz

      grp%wallStress(ip,1) = (c1*wallNormal(1)&
                           +  c2*wallNormal(2)&
                           +  c3*wallNormal(3))
      grp%wallStress(ip,2) = (c2*wallNormal(1)&
                           +   c4*wallNormal(2)&
                           +   c5*wallNormal(3))
      grp%wallStress(ip,3) = (c3*wallNormal(1)&
                           +   c5*wallNormal(2)&
                           +   c6*wallNormal(3))
      grp%wallStress(ip,4) = 0.0
    end do
  end if


  if(ivd%wallsAreIsentropic) then
    ! this is OK since the tangential component of the temperature gradient
    ! will be calculated using the finite difference term 
    do i=1,grp%brp%numberOfBoundaryNodes
      grp%nodeHelpArray(i,4) = 0.0
      grp%nodeHelpArray(i,8) = 0.0
      grp%nodeHelpArray(i,12) = 0.0
    end do
  end if

  grp%nodeHelpArray(:,15:18) = 0.0

#ifdef PARALLEL

do kk=1,2
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfComSides
 else
  istart = grp%pdp%numberOfComSides + 1
  ifinish = grp%numberOfSides
 end if


 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL

 i1 = grp%sideIndexArray(i,1)
 i2 = grp%sideIndexArray(i,2)
 inverseSideLength = 1.0/grp%sideLengthArray(i)

! side weights from preprocessor 
 wx = grp%sideWeightsArray(i,1)
 wy = grp%sideWeightsArray(i,2)  
 wz = grp%sideWeightsArray(i,3)  

 mmu1 = grp%nodeHelpArray(i1,13)
 mmu2 = grp%nodeHelpArray(i2,13)

 k1 = grp%nodeHelpArray(i1,14)
 k2 = grp%nodeHelpArray(i2,14)


 rhoinv1 = 1.0/grp%u(i1,1)
 u1 = grp%u(i1,2)*rhoinv1
 v1 = grp%u(i1,3)*rhoinv1
 w1 = grp%u(i1,4)*rhoinv1
 T1 = gammaOverGammaMinusOne*grp%p(i1)*rhoinv1

 rhoinv2 = 1.0/grp%u(i2,1)
 u2 = grp%u(i2,2)*rhoinv2
 v2 = grp%u(i2,3)*rhoinv2
 w2 = grp%u(i2,4)*rhoinv2
 T2 = gammaOverGammaMinusOne*grp%p(i2)*rhoinv2

 r(1) = grp%coordinates(i2,1)-grp%coordinates(i1,1)
 r(2) = grp%coordinates(i2,2)-grp%coordinates(i1,2)
 r(3) = grp%coordinates(i2,3)-grp%coordinates(i1,3)
 r = r*inverseSideLength

 dudx1 = grp%nodeHelpArray(i1,1)
 dudy1 = grp%nodeHelpArray(i1,5)
 dudz1 = grp%nodeHelpArray(i1,9)
 dudx2 = grp%nodeHelpArray(i2,1)
 dudy2 = grp%nodeHelpArray(i2,5)
 dudz2 = grp%nodeHelpArray(i2,9)
 dvdx1 = grp%nodeHelpArray(i1,2)
 dvdy1 = grp%nodeHelpArray(i1,6)
 dvdz1 = grp%nodeHelpArray(i1,10)
 dvdx2 = grp%nodeHelpArray(i2,2)
 dvdy2 = grp%nodeHelpArray(i2,6)
 dvdz2 = grp%nodeHelpArray(i2,10)
 dwdx1 = grp%nodeHelpArray(i1,3)
 dwdy1 = grp%nodeHelpArray(i1,7)
 dwdz1 = grp%nodeHelpArray(i1,11)
 dwdx2 = grp%nodeHelpArray(i2,3)
 dwdy2 = grp%nodeHelpArray(i2,7)
 dwdz2 = grp%nodeHelpArray(i2,11)
 dTdx1 = grp%nodeHelpArray(i1,4)
 dTdy1 = grp%nodeHelpArray(i1,8)
 dTdz1 = grp%nodeHelpArray(i1,12)
 dTdx2 = grp%nodeHelpArray(i2,4)
 dTdy2 = grp%nodeHelpArray(i2,8)
 dTdz2 = grp%nodeHelpArray(i2,12)

 dudx = dudx1 + dudx2
 dudy = dudy1 + dudy2
 dudz = dudz1 + dudz2
 dvdx = dvdx1 + dvdx2
 dvdy = dvdy1 + dvdy2
 dvdz = dvdz1 + dvdz2
 dwdx = dwdx1 + dwdx2
 dwdy = dwdy1 + dwdy2
 dwdz = dwdz1 + dwdz2
 dTdx = dTdx1 + dTdx2
 dTdy = dTdy1 + dTdy2
 dTdz = dTdz1 + dTdz2

! remove derivative component along side

 dudr = dudx*r(1)+dudy*r(2)+dudz*r(3)
 dvdr = dvdx*r(1)+dvdy*r(2)+dvdz*r(3)
 dwdr = dwdx*r(1)+dwdy*r(2)+dwdz*r(3)
 dTdr = dTdx*r(1)+dTdy*r(2)+dTdz*r(3)

 dudx = dudx - dudr*r(1)
 dudy = dudy - dudr*r(2)
 dudz = dudz - dudr*r(3)
 dvdx = dvdx - dvdr*r(1)
 dvdy = dvdy - dvdr*r(2)
 dvdz = dvdz - dvdr*r(3)
 dwdx = dwdx - dwdr*r(1)
 dwdy = dwdy - dwdr*r(2)
 dwdz = dwdz - dwdr*r(3)
 dTdx = dTdx - dTdr*r(1)
 dTdy = dTdy - dTdr*r(2)
 dTdz = dTdz - dTdr*r(3)


! add finite difference along edge

 du = (mmu1+mmu2)*(u2-u1)*inverseSideLength  ! factor of 2 since wx,wy,wz have factor of 0.5
 dv = (mmu1+mmu2)*(v2-v1)*inverseSideLength
 dw = (mmu1+mmu2)*(w2-w1)*inverseSideLength
 dT = (k1+k2)*(T2-T1)*inverseSideLength


 dudx = dudx + du*r(1)
 dudy = dudy + du*r(2)
 dudz = dudz + du*r(3)
 dvdx = dvdx + dv*r(1)
 dvdy = dvdy + dv*r(2)
 dvdz = dvdz + dv*r(3)
 dwdx = dwdx + dw*r(1)
 dwdy = dwdy + dw*r(2)
 dwdz = dwdz + dw*r(3)
 dTdx = dTdx + dT*r(1)
 dTdy = dTdy + dT*r(2)
 dTdz = dTdz + dT*r(3)

! stress tensor

 dt1 = fourThirds*dudx-twoThirds*(dvdy+dwdz)
 dt2 = dudy+dvdx
 dt3 = dudz+dwdx
 dt4 = fourThirds*dvdy-twoThirds*(dudx+dwdz)
 dt5 = dvdz+dwdy
 dt6 = fourThirds*dwdz-twoThirds*(dudx+dvdy)
 dt7 = dTdx
 dt8 = dTdy
 dt9 = dTdz

 dut1 = 0.5*(u1+u2)*dt1+0.5*(v1+v2)*dt2+0.5*(w1+w2)*dt3
 dut2 = 0.5*(u1+u2)*dt2+0.5*(v1+v2)*dt4+0.5*(w1+w2)*dt5
 dut3 = 0.5*(u1+u2)*dt3+0.5*(v1+v2)*dt5+0.5*(w1+w2)*dt6

 f2 = wx*dt1 + wy*dt2 + wz*dt3
 f3 = wx*dt2 + wy*dt4 + wz*dt5
 f4 = wx*dt3 + wy*dt5 + wz*dt6
 f5 = wx*(dt7+dut1) + wy*(dt8+dut2) + wz*(dt9+dut3)

 grp%nodeHelpArray(i1,15) = grp%nodeHelpArray(i1,15) + f2 
 grp%nodeHelpArray(i1,16) = grp%nodeHelpArray(i1,16) + f3  
 grp%nodeHelpArray(i1,17) = grp%nodeHelpArray(i1,17) + f4 
 grp%nodeHelpArray(i1,18) = grp%nodeHelpArray(i1,18) + f5

 grp%nodeHelpArray(i2,15) = grp%nodeHelpArray(i2,15) - f2 
 grp%nodeHelpArray(i2,16) = grp%nodeHelpArray(i2,16) - f3  
 grp%nodeHelpArray(i2,17) = grp%nodeHelpArray(i2,17) - f4 
 grp%nodeHelpArray(i2,18) = grp%nodeHelpArray(i2,18) - f5
end do


#ifdef PARALLEL
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfBoundaryComSides
 else
  istart = grp%pdp%numberOfBoundaryComSides + 1
  ifinish = grp%brp%numberOfBoundarySides
 end if
 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL

 faceIndicator = grp%brp%sideIndexArray(i,3)
 if(faceIndicator.ne.2.and.faceIndicator.ne.9) then 
 i1 = grp%brp%sideIndexArray(i,1)
 i2 = grp%brp%sideIndexArray(i,2)

 wx = grp%brp%sideWeightsArray(i,1)
 wy = grp%brp%sideWeightsArray(i,2)
 wz = grp%brp%sideWeightsArray(i,3)

 mmu1 = grp%nodeHelpArray(i1,13)
 mmu2 = grp%nodeHelpArray(i2,13)

 k1 = grp%nodeHelpArray(i1,14)
 k2 = grp%nodeHelpArray(i2,14)


 rhoinv1 = 1.0/grp%u(i1,1)
 u1 = grp%u(i1,2)*rhoinv1
 v1 = grp%u(i1,3)*rhoinv1
 w1 = grp%u(i1,4)*rhoinv1
 T1 = gammaOverGammaMinusOne*grp%p(i1)*rhoinv1

 rhoinv2 = 1.0/grp%u(i2,1)
 u2 = grp%u(i2,2)*rhoinv2
 v2 = grp%u(i2,3)*rhoinv2
 w2 = grp%u(i2,4)*rhoinv2
 T2 = gammaOverGammaMinusOne*grp%p(i2)*rhoinv2

 dudx1 = grp%nodeHelpArray(i1,1)
 dudy1 = grp%nodeHelpArray(i1,5)
 dudz1 = grp%nodeHelpArray(i1,9)
 dudx2 = grp%nodeHelpArray(i2,1)
 dudy2 = grp%nodeHelpArray(i2,5)
 dudz2 = grp%nodeHelpArray(i2,9)
 dvdx1 = grp%nodeHelpArray(i1,2)
 dvdy1 = grp%nodeHelpArray(i1,6)
 dvdz1 = grp%nodeHelpArray(i1,10)
 dvdx2 = grp%nodeHelpArray(i2,2)
 dvdy2 = grp%nodeHelpArray(i2,6)
 dvdz2 = grp%nodeHelpArray(i2,10)
 dwdx1 = grp%nodeHelpArray(i1,3)
 dwdy1 = grp%nodeHelpArray(i1,7)
 dwdz1 = grp%nodeHelpArray(i1,11)
 dwdx2 = grp%nodeHelpArray(i2,3)
 dwdy2 = grp%nodeHelpArray(i2,7)
 dwdz2 = grp%nodeHelpArray(i2,11)
 dTdx1 = grp%nodeHelpArray(i1,4)
 dTdy1 = grp%nodeHelpArray(i1,8)
 dTdz1 = grp%nodeHelpArray(i1,12)
 dTdx2 = grp%nodeHelpArray(i2,4)
 dTdy2 = grp%nodeHelpArray(i2,8)
 dTdz2 = grp%nodeHelpArray(i2,12)

 dudx = dudx1 + dudx2
 dudy = dudy1 + dudy2
 dudz = dudz1 + dudz2
 dvdx = dvdx1 + dvdx2
 dvdy = dvdy1 + dvdy2
 dvdz = dvdz1 + dvdz2
 dwdx = dwdx1 + dwdx2
 dwdy = dwdy1 + dwdy2
 dwdz = dwdz1 + dwdz2
 dTdx = dTdx1 + dTdx2
 dTdy = dTdy1 + dTdy2
 dTdz = dTdz1 + dTdz2

! stress tensor

 t11 = fourThirds*dudx1-twoThirds*(dvdy1+dwdz1)
 t21 = dudy1+dvdx1
 t31 = dudz1+dwdx1
 t41 = fourThirds*dvdy1-twoThirds*(dudx1+dwdz1)
 t51 = dvdz1+dwdy1
 t61 = fourThirds*dwdz1-twoThirds*(dudx1+dvdy1)
 t71 = dTdx1
 t81 = dTdy1
 t91 = dTdz1

 t12 = fourThirds*dudx2-twoThirds*(dvdy2+dwdz2)
 t22 = dudy2+dvdx2
 t32 = dudz2+dwdx2
 t42 = fourThirds*dvdy2-twoThirds*(dudx2+dwdz2)
 t52 = dvdz2+dwdy2
 t62 = fourThirds*dwdz2-twoThirds*(dudx2+dvdy2)
 t72 = dTdx2
 t82 = dTdy2
 t92 = dTdz2

 dt1 = fourThirds*dudx-twoThirds*(dvdy+dwdz)
 dt2 = dudy+dvdx
 dt3 = dudz+dwdx
 dt4 = fourThirds*dvdy-twoThirds*(dudx+dwdz)
 dt5 = dvdz+dwdy
 dt6 = fourThirds*dwdz-twoThirds*(dudx+dvdy)
 dt7 = dTdx
 dt8 = dTdy
 dt9 = dTdz

! isentropic boundary

 if(ivd%wallsAreIsentropic) then 
  dt7 = 0.0
  dt8 = 0.0
  dt9 = 0.0
  t71 = 0.0
  t81 = 0.0
  t91 = 0.0 
  t72 = 0.0
  t82 = 0.0
  t92 = 0.0 
 end if

 dut11 = u1*t11 + v1*t21 + w1*t31
 dut21 = u1*t21 + v1*t41 + w1*t51
 dut31 = u1*t31 + v1*t51 + w1*t61

 dut12 = u2*t12 + v2*t22 + w2*t32
 dut22 = u2*t22 + v2*t42 + w2*t52
 dut32 = u2*t32 + v2*t52 + w2*t62

 dut1 = 0.5*(u1+u2)*dt1+0.5*(v1+v2)*dt2+0.5*(w1+w2)*dt3
 dut2 = 0.5*(u1+u2)*dt2+0.5*(v1+v2)*dt4+0.5*(w1+w2)*dt5
 dut3 = 0.5*(u1+u2)*dt3+0.5*(v1+v2)*dt5+0.5*(w1+w2)*dt6

 grp%nodeHelpArray(i1,15) = grp%nodeHelpArray(i1,15) + 4.0*(wx*t11+wy*t21+wz*t31)
 grp%nodeHelpArray(i1,16) = grp%nodeHelpArray(i1,16) + 4.0*(wx*t21+wy*t41+wz*t51)
 grp%nodeHelpArray(i1,17) = grp%nodeHelpArray(i1,17) + 4.0*(wx*t31+wy*t51+wz*t61)
 grp%nodeHelpArray(i1,18) = grp%nodeHelpArray(i1,18) + 4.0*(wx*(t71+dut11) + wy*(t81+dut21)+wz*(t91+dut31))

 grp%nodeHelpArray(i2,15) = grp%nodeHelpArray(i2,15) + 4.0*(wx*t12+wy*t22+wz*t32)
 grp%nodeHelpArray(i2,16) = grp%nodeHelpArray(i2,16) + 4.0*(wx*t22+wy*t42+wz*t52)
 grp%nodeHelpArray(i2,17) = grp%nodeHelpArray(i2,17) + 4.0*(wx*t32+wy*t52+wz*t62)
 grp%nodeHelpArray(i2,18) = grp%nodeHelpArray(i2,18) + 4.0*(wx*(t72+dut12) + wy*(t82+dut22)+wz*(t92+dut32))
 end if
end do


#ifdef PARALLEL

 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses
   if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
    do j=15,18
     call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
    call mp_sendv(grp%pdp%processorIDs,i,31,grp%pdp%sendFlag)
   end if
  end do
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,31,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
    do j=15,18
     call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                            grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
  end if
 end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position = 1
   do j=15,18
    call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                      grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
   end do
   call mp_sendv(grp%pdp%processorIDs,i,32,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,32,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )


 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
    do j=15,18
     call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
  end if
 end do

#endif PARALLEL

! add viscosity to dissipation vector

 do i=1,grp%numberOfNodes
  rhs(i,2:5) = rhs(i,2:5) - grp%nodeHelpArray(i,15:18)
 end do
end subroutine makeViscosityTerm2
!-------------------------------------------------------------------------
subroutine makeSARHS(grp,ivd) 
IMPLICIT NONE 
include 'mpif.h'
! convective term for the Spalart-Allmaras turbulence model

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

 real :: turbVisc1,turbVisc2,ru1,rv1,rw1,u1,v1,w1,rinv1,rinv2
 real :: ru2,rv2,rw2,u2,v2,w2,wx,wy,wz,fconv,ax,ay,az

 integer :: i,ind1,ind2,ind

 integer :: kk,istart,ifinish

#ifdef PARALLEL

do kk=1,2
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfComSides
 else
  istart = grp%pdp%numberOfComSides + 1
  ifinish = grp%numberOfSides
 end if


 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL

  ind1 = grp%sideIndexArray(i,1)
  ind2 = grp%sideIndexArray(i,2)

  rinv1 = 1./grp%u(ind1,1)                 ! rho inverse
  turbVisc1 = grp%u(ind1,6)
  ru1 = grp%u(ind1,2)                      ! x-momentum
  rv1 = grp%u(ind1,3)                      ! y-momentum
  rw1 = grp%u(ind1,4)                      ! z-momentum
  u1 = rinv1*ru1                           ! x-velocity
  v1 = rinv1*rv1                           ! y-velocity
  w1 = rinv1*rw1                           ! z-velocity
 
  rinv2 = 1./grp%u(ind2,1)
  turbVisc2 = grp%u(ind2,6)
  ru2 = grp%u(ind2,2)
  rv2 = grp%u(ind2,3)
  rw2 = grp%u(ind2,4)
  u2 = rinv2*ru2
  v2 = rinv2*rv2
  w2 = rinv2*rw2

! side weights from preprocessor

  wx = grp%sideWeightsArray(i,1)
  wy = grp%sideWeightsArray(i,2)
  wz = grp%sideWeightsArray(i,3)

  ! convection terms

  fconv = wx*(u1*turbVisc1 + u2*turbVisc2) + wy*(v1*turbVisc1 + v2*turbVisc2) + wz*(w1*turbvisc1 + w2*turbvisc2)
  if(abs(turbvisc1-turbvisc2)>1.0e-20) then
   ax = (u2*turbVisc2 - u1*turbVisc1)/(turbvisc1-turbvisc2)
   ay = (v2*turbVisc2 - v1*turbVisc1)/(turbvisc1-turbvisc2)
   az = (w2*turbVisc2 - w1*turbVisc1)/(turbvisc1-turbvisc2)

   fconv = fconv + abs(wx*ax+wy*ay+wz*az)*(turbvisc1-turbvisc2)
  end if

  grp%rhs(ind1,6) = grp%rhs(ind1,6) + fconv
  grp%rhs(ind2,6) = grp%rhs(ind2,6) - fconv
 end do

#ifdef PARALLEL

!    write(get_unit(),*)  'Starting first SAHRS comms',kk
 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses
   if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
    call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%rhs(:,6),grp%pdp%buffer,grp%pdp%sendFlag)
    call mp_sendv(grp%pdp%processorIDs,i,33,grp%pdp%sendFlag)
   end if
  end do
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,33,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
    call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                           grp%pdp%realType,grp%rhs(:,6),grp%pdp%buffer,grp%pdp%sendFlag)
  end if
 end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position = 1
   call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                     grp%pdp%realType,grp%rhs(:,6),grp%pdp%buffer,grp%pdp%sendFlag)
   call mp_sendv(grp%pdp%processorIDs,i,34,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,34,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
    call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%rhs(:,6),grp%pdp%buffer,grp%pdp%sendFlag)
  end if
 end do

#endif PARALLEL

end subroutine makeSARHS
!-----------------------------------------------------------------------
!START-TURB
subroutine makeKWRHS(grp,ivd) !MOD 
IMPLICIT NONE 
include 'mpif.h'
! convective terms for the k-omega turbulence model equation

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

 real :: kw1,kw2,ru1,rv1,rw1,u1,v1,w1,rinv1,rinv2
 real :: ru2,rv2,rw2,u2,v2,w2,wx,wy,wz,fconv,ax,ay,az,fconv1,fconv2

 integer :: i,ind1,ind2,ind,ij

 integer :: kk,istart,ifinish,IEST,IBOUND
IEST=0 !0 includes artificial diss
IBOUND=1 !0 includes face contributions


!Loop for k(ij=6) and omega(ij=7)
do ij=6,7 

#ifdef PARALLEL

do kk=1,2
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfComSides
 else
  istart = grp%pdp%numberOfComSides + 1
  ifinish = grp%numberOfSides
 end if


 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL
 
  ind1 = grp%sideIndexArray(i,1)
  ind2 = grp%sideIndexArray(i,2)

  kw1 = grp%u(ind1,ij)
  ru1 = grp%u(ind1,2)                      ! x-momentum
  rv1 = grp%u(ind1,3)                      ! y-momentum
  rw1 = grp%u(ind1,4)                      ! z-momentum
 
  kw2 = grp%u(ind2,ij)
  ru2 = grp%u(ind2,2)
  rv2 = grp%u(ind2,3)
  rw2 = grp%u(ind2,4)

! side weights from preprocessor

  wx = grp%sideWeightsArray(i,1)
  wy = grp%sideWeightsArray(i,2)
  wz = grp%sideWeightsArray(i,3)

  fconv = wx*(ru1*kw1 + ru2*kw2) + wy*(rv1*kw1 + rv2*kw2) + wz*(rw1*kw1 + rw2*kw2)
 
IF(IEST.EQ.0)THEN
  if(abs(kw1-kw2)>1.0e-20) then
   ax = (ru2*kw2 - ru1*kw1)/(kw1-kw2)
   ay = (rv2*kw2 - rv1*kw1)/(kw1-kw2)
   az = (rw2*kw2 - rw1*kw1)/(kw1-kw2)

   fconv = fconv + abs(wx*ax+wy*ay+wz*az)*(kw1-kw2)
   
  end if
END IF
  grp%rhs(ind1,ij) = grp%rhs(ind1,ij) + fconv
  grp%rhs(ind2,ij) = grp%rhs(ind2,ij) - fconv
 

 end do 

IF(IBOUND.EQ.0)THEN
!adding boundary face contribution

#ifdef PARALLEL
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfBoundaryComSides
 else
  istart = grp%pdp%numberOfBoundaryComSides + 1
  ifinish = grp%brp%numberOfBoundarySides
 end if
 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL


  ind1 = grp%sideIndexArray(i,1)
  ind2 = grp%sideIndexArray(i,2)

  kw1 = grp%u(ind1,ij)
  ru1 = grp%u(ind1,2)                      ! x-momentum
  rv1 = grp%u(ind1,3)                      ! y-momentum
  rw1 = grp%u(ind1,4)                      ! z-momentum

  kw2 = grp%u(ind2,ij)
  ru2 = grp%u(ind2,2)
  rv2 = grp%u(ind2,3)
  rw2 = grp%u(ind2,4)

  wx = grp%sideWeightsArray(i,1)
  wy = grp%sideWeightsArray(i,2)
  wz = grp%sideWeightsArray(i,3)

  if(ivd%boundaryTerm==1) then
   fconv1 = wx*(3.*ru1*kw1+ru2*kw2)+wy*(3.*rv1*kw1+rv2*kw2)+wz*(3.*rw1*kw1+rw2*kw2)
   fconv2 = wx*(ru1*kw1+3.*ru2*kw2)+wy*(rv1*kw1+3.*rv2*kw2)+wz*(rw1*kw1+3.*rw2*kw2)
  else 
   fconv1 = 4.0*(wx*ru1*kw1+wy*rv1*kw1+wz*rw1*kw1)
   fconv2 = 4.0*(wx*ru2*kw2+wy*rv2*kw2+wz*rw2*kw2)
  end if

    grp%rhs(ind1,ij) = grp%rhs(ind1,ij) + fconv1
    grp%rhs(ind2,ij) = grp%rhs(ind2,ij) + fconv2

IF(IEST.EQ.0)THEN
  if(abs(kw1-kw2)>1.0e-20) then
   ax = (ru2*kw2 - ru1*kw1)/(kw1-kw2)
   ay = (rv2*kw2 - rv1*kw1)/(kw1-kw2)
   az = (rw2*kw2 - rw1*kw1)/(kw1-kw2)
   fconv = abs(wx*ax+wy*ay+wz*az)*(kw1-kw2)

   grp%rhs(ind1,ij) = grp%rhs(ind1,ij) + fconv
   grp%rhs(ind2,ij) = grp%rhs(ind2,ij) + fconv
  end if
END IF

!! end if 
 end do

END IF !ADDING BOUNDARY CONTRIBUTION

#ifdef PARALLEL

!    write(get_unit(),*)  'Starting first SAHRS comms',kk
 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses
   if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
    call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%rhs(:,ij),grp%pdp%buffer,grp%pdp%sendFlag)
    call mp_sendv(grp%pdp%processorIDs,i,33,grp%pdp%sendFlag)
   end if
  end do
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,33,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1235 continue

 do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
!  grp%pdp%sendFlag = 0
!  call mp_check_mesg_arrived(grp%pdp%processorIDs,i,33,grp%pdp%sendFlag,bufferzz,bufferSize)
!  if (grp%pdp%sendFlag.gt.0) then
!   grp%pdp%receivedMessageFromProcess(i) = 1
!   position = 1
    call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                           grp%pdp%realType,grp%rhs(:,ij),grp%pdp%buffer,grp%pdp%sendFlag)
!  end if
  end if
 end do

!do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1235
!end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

!  write(get_unit(),*)  'Starting second SAHRS comms'

 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position = 1
   call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                     grp%pdp%realType,grp%rhs(:,ij),grp%pdp%buffer,grp%pdp%sendFlag)
   call mp_sendv(grp%pdp%processorIDs,i,34,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,34,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1335 continue

 do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
!  grp%pdp%sendFlag = 0
!  call mp_check_mesg_arrived(grp%pdp%processorIDs,i,34,grp%pdp%sendFlag,bufferzz,bufferSize)
!  if (grp%pdp%sendFlag.gt.0) then
!   grp%pdp%receivedMessageFromProcess(i) = 1
!   position = 1
    call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%rhs(:,ij),grp%pdp%buffer,grp%pdp%sendFlag)
!  end if
  end if
 end do

!do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1335
!end do
#endif PARALLEL

end do!end loop for k and omega

end subroutine makeKWRHS
!END-TURB
!-----------------------------------------------------------------------
!START-TURB
subroutine makeSSTRHS(grp,ivd) !SST
IMPLICIT NONE 
include 'mpif.h'
! convective terms for the SST turbulence model equation

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

 real :: kw1,kw2,ru1,rv1,rw1,u1,v1,w1,rinv1,rinv2
 real :: ru2,rv2,rw2,u2,v2,w2,wx,wy,wz,fconv,ax,ay,az,fconv1,fconv2

 integer :: i,ind1,ind2,ind,ij

 integer :: kk,istart,ifinish,IEST,IBOUND
IEST=0 !0 includes artificial diss
IBOUND=1 !0 includes face contributions


!Loop for k(ij=6) and omega(ij=7)
do ij=6,7 

#ifdef PARALLEL

do kk=1,2
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfComSides
 else
  istart = grp%pdp%numberOfComSides + 1
  ifinish = grp%numberOfSides
 end if


 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL
 
  ind1 = grp%sideIndexArray(i,1)
  ind2 = grp%sideIndexArray(i,2)

  kw1 = grp%u(ind1,ij)
  ru1 = grp%u(ind1,2)                      ! x-momentum
  rv1 = grp%u(ind1,3)                      ! y-momentum
  rw1 = grp%u(ind1,4)                      ! z-momentum
 
  kw2 = grp%u(ind2,ij)
  ru2 = grp%u(ind2,2)
  rv2 = grp%u(ind2,3)
  rw2 = grp%u(ind2,4)

! side weights from preprocessor

  wx = grp%sideWeightsArray(i,1)
  wy = grp%sideWeightsArray(i,2)
  wz = grp%sideWeightsArray(i,3)

  fconv = wx*(ru1*kw1 + ru2*kw2) + wy*(rv1*kw1 + rv2*kw2) + wz*(rw1*kw1 + rw2*kw2)
 

!artificial dissipation term!
IF(IEST.EQ.0)THEN
  if(abs(kw1-kw2)>1.0e-20) then
   ax = (ru2*kw2 - ru1*kw1)/(kw1-kw2)
   ay = (rv2*kw2 - rv1*kw1)/(kw1-kw2)
   az = (rw2*kw2 - rw1*kw1)/(kw1-kw2)

   fconv = fconv + abs(wx*ax+wy*ay+wz*az)*(kw1-kw2)
   
  end if
END IF
  grp%rhs(ind1,ij) = grp%rhs(ind1,ij) + fconv
  grp%rhs(ind2,ij) = grp%rhs(ind2,ij) - fconv
 

 end do 

IF(IBOUND.EQ.0)THEN
!adding boundary face contribution

#ifdef PARALLEL
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfBoundaryComSides
 else
  istart = grp%pdp%numberOfBoundaryComSides + 1
  ifinish = grp%brp%numberOfBoundarySides
 end if
 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL


  ind1 = grp%sideIndexArray(i,1)
  ind2 = grp%sideIndexArray(i,2)

  kw1 = grp%u(ind1,ij)
  ru1 = grp%u(ind1,2)                      ! x-momentum
  rv1 = grp%u(ind1,3)                      ! y-momentum
  rw1 = grp%u(ind1,4)                      ! z-momentum

  kw2 = grp%u(ind2,ij)
  ru2 = grp%u(ind2,2)
  rv2 = grp%u(ind2,3)
  rw2 = grp%u(ind2,4)

  wx = grp%sideWeightsArray(i,1)
  wy = grp%sideWeightsArray(i,2)
  wz = grp%sideWeightsArray(i,3)

  if(ivd%boundaryTerm==1) then
   fconv1 = wx*(3.*ru1*kw1+ru2*kw2)+wy*(3.*rv1*kw1+rv2*kw2)+wz*(3.*rw1*kw1+rw2*kw2)
   fconv2 = wx*(ru1*kw1+3.*ru2*kw2)+wy*(rv1*kw1+3.*rv2*kw2)+wz*(rw1*kw1+3.*rw2*kw2)
  else 
   fconv1 = 4.0*(wx*ru1*kw1+wy*rv1*kw1+wz*rw1*kw1)
   fconv2 = 4.0*(wx*ru2*kw2+wy*rv2*kw2+wz*rw2*kw2)
  end if

    grp%rhs(ind1,ij) = grp%rhs(ind1,ij) + fconv1
    grp%rhs(ind2,ij) = grp%rhs(ind2,ij) + fconv2

!artificial dissipation term!
IF(IEST.EQ.0)THEN
  if(abs(kw1-kw2)>1.0e-20) then
   ax = (ru2*kw2 - ru1*kw1)/(kw1-kw2)
   ay = (rv2*kw2 - rv1*kw1)/(kw1-kw2)
   az = (rw2*kw2 - rw1*kw1)/(kw1-kw2)
   fconv = abs(wx*ax+wy*ay+wz*az)*(kw1-kw2)

   grp%rhs(ind1,ij) = grp%rhs(ind1,ij) + fconv
   grp%rhs(ind2,ij) = grp%rhs(ind2,ij) + fconv
  end if
END IF

!! end if 
 end do

END IF !ADDING BOUNDARY CONTRIBUTION

#ifdef PARALLEL

!    write(get_unit(),*)  'Starting first SAHRS comms',kk
 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses
   if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
    call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%rhs(:,ij),grp%pdp%buffer,grp%pdp%sendFlag)
    call mp_sendv(grp%pdp%processorIDs,i,33,grp%pdp%sendFlag)
   end if
  end do
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,33,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1235 continue

 do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
!  grp%pdp%sendFlag = 0
!  call mp_check_mesg_arrived(grp%pdp%processorIDs,i,33,grp%pdp%sendFlag,bufferzz,bufferSize)
!  if (grp%pdp%sendFlag.gt.0) then
!   grp%pdp%receivedMessageFromProcess(i) = 1
!   position = 1
    call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                           grp%pdp%realType,grp%rhs(:,ij),grp%pdp%buffer,grp%pdp%sendFlag)
!  end if
  end if
 end do

!do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1235
!end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

!  write(get_unit(),*)  'Starting second SAHRS comms'

 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position = 1
   call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                     grp%pdp%realType,grp%rhs(:,ij),grp%pdp%buffer,grp%pdp%sendFlag)
   call mp_sendv(grp%pdp%processorIDs,i,34,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,34,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1335 continue

 do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
!  grp%pdp%sendFlag = 0
!  call mp_check_mesg_arrived(grp%pdp%processorIDs,i,34,grp%pdp%sendFlag,bufferzz,bufferSize)
!  if (grp%pdp%sendFlag.gt.0) then
!   grp%pdp%receivedMessageFromProcess(i) = 1
!   position = 1
    call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%rhs(:,ij),grp%pdp%buffer,grp%pdp%sendFlag)
!  end if
  end if
 end do

!do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1335
!end do
#endif PARALLEL


end do!end loop for k and omega

end subroutine makeSSTRHS
!END-TURB
!-----------------------------------------------------------------------
subroutine makeSADiffusion(grp,ivd)
IMPLICIT NONE
include 'mpif.h'
! diffusion term for the Spalart-Allmaras turbulence model

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

 integer :: i,j,ind1,ind2,ind
 real :: ru1,ru2,rv1,rv2,rw1,rw2,fact,vt,visc1,visc2,fdiff1,fdiff2,fdiff3
 real :: turbVisc1,turbVisc2,tvx1,tvx2,tvy1,tvy2,tvz1,tvz2,dnudx,dnudy,dnudz
 real :: sigma,cb2,wx,wy,wz,dnudx1,dnudx2,dnudy1,dnudy2,dnudz1,dnudz2
 real :: vecProd,oneOverReynoldsNumber

 integer :: kk,istart,ifinish


 sigma = 2./3.
 cb2 = 0.622

 oneOverReynoldsNumber = 1./ivd%ReynoldsNumber

! diffusion terms

 grp%turbulenceDiffusionTerm = 0.0
 grp%nodeHelpArray(1:grp%numberOfNodes,9:14) = 0.0

#ifdef PARALLEL

do kk=1,2
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfComSides
 else
  istart = grp%pdp%numberOfComSides + 1
  ifinish = grp%numberOfSides
 end if


 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL

  ind1 = grp%sideIndexArray(i,1)
  ind2 = grp%sideIndexArray(i,2)

  turbVisc1 = grp%u(ind1,6)
  turbVisc2 = grp%u(ind2,6)

! side weights from preprocessor
  wx = grp%sideWeightsArray(i,1)
  wy = grp%sideWeightsArray(i,2)
  wz = grp%sideWeightsArray(i,3)


  fdiff1 = wx*(turbVisc1+turbVisc2)
  fdiff2 = wy*(turbVisc1+turbVisc2)
  fdiff3 = wz*(turbVisc1+turbVisc2)

  grp%nodeHelpArray(ind1,9) = grp%nodeHelpArray(ind1,9) + fdiff1
  grp%nodeHelpArray(ind2,9) = grp%nodeHelpArray(ind2,9) - fdiff1
  grp%nodeHelpArray(ind1,10) = grp%nodeHelpArray(ind1,10) + fdiff2
  grp%nodeHelpArray(ind2,10) = grp%nodeHelpArray(ind2,10) - fdiff2
  grp%nodeHelpArray(ind1,11) = grp%nodeHelpArray(ind1,11) + fdiff3
  grp%nodeHelpArray(ind2,11) = grp%nodeHelpArray(ind2,11) - fdiff3
 end do
 ! boundary faces

#ifdef PARALLEL
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfBoundaryComSides
 else
  istart = grp%pdp%numberOfBoundaryComSides + 1
  ifinish = grp%brp%numberOfBoundarySides
 end if
 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL

  ind1 = grp%brp%sideIndexArray(i,1)
  ind2 = grp%brp%sideIndexArray(i,2)
  wx   = grp%brp%sideWeightsArray(i,1)
  wy   = grp%brp%sideWeightsArray(i,2)
  wz   = grp%brp%sideWeightsArray(i,3)

  turbVisc1 = grp%u(ind1,6)
  turbVisc2 = grp%u(ind2,6)

 ! add boundary face contributions

  if(ivd%boundaryTerm==1) then 
   tvx1 = (3.*turbVisc1 + turbVisc2)*wx
   tvx2 = (turbVisc1 + 3.*turbVisc2)*wx
   tvy1 = (3.*turbVisc1 + turbVisc2)*wy
   tvy2 = (turbVisc1 + 3.*turbVisc2)*wy
   tvz1 = (3.*turbVisc1 + turbVisc2)*wz
   tvz2 = (turbVisc1 + 3.*turbVisc2)*wz
 
   grp%nodeHelpArray(ind1,9) = grp%nodeHelpArray(ind1,9) + tvx1
   grp%nodeHelpArray(ind2,9) = grp%nodeHelpArray(ind2,9) + tvx2
   grp%nodeHelpArray(ind1,10) = grp%nodeHelpArray(ind1,10) + tvy1
   grp%nodeHelpArray(ind2,10) = grp%nodeHelpArray(ind2,10) + tvy2
   grp%nodeHelpArray(ind1,11) = grp%nodeHelpArray(ind1,11) + tvz1
   grp%nodeHelpArray(ind2,11) = grp%nodeHelpArray(ind2,11) + tvz2
  else
   grp%nodeHelpArray(ind1,9) = grp%nodeHelpArray(ind1,9) + 4.*turbVisc1*wx
   grp%nodeHelpArray(ind2,9) = grp%nodeHelpArray(ind2,9) + 4.*turbVisc2*wx
   grp%nodeHelpArray(ind1,10) = grp%nodeHelpArray(ind1,10) + 4.*turbVisc1*wy
   grp%nodeHelpArray(ind2,10) = grp%nodeHelpArray(ind2,10) + 4.*turbvisc2*wy
   grp%nodeHelpArray(ind1,11) = grp%nodeHelpArray(ind1,11) + 4.*turbvisc1*wz
   grp%nodeHelpArray(ind2,11) = grp%nodeHelpArray(ind2,11) + 4.*turbvisc2*wz
  end if
 end do
#ifdef PARALLEL

 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses

   if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
    do j=9,11
     call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
    call mp_sendv(grp%pdp%processorIDs,i,35,grp%pdp%sendFlag)
   end if
  end do
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,35,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
    do j=9,11
     call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                            grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
  end if
 end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position = 1
   do j=9,11
    call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                      grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
   end do
   call mp_sendv(grp%pdp%processorIDs,i,36,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,36,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
    do j=9,11
     call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
  end if
 end do

#endif PARALLEL

 grp%nodeHelpArray(:,1) = 0.0

 fact = oneOverReynoldsNumber*cb2/sigma
 do i=1,grp%numberOfNodes
  vt = 1./grp%nodeVolume(i)
  dnudx = grp%nodeHelpArray(i,9)*vt
  dnudy = grp%nodeHelpArray(i,10)*vt
  dnudz = grp%nodeHelpArray(i,11)*vt
  grp%nodeHelpArray(i,12) = dnudx
  grp%nodeHelpArray(i,13) = dnudy
  grp%nodeHelpArray(i,14) = dnudz

  vecProd = grp%nodeHelpArray(i,9)*dnudx + grp%nodeHelpArray(i,10)*dnudy + grp%nodeHelpArray(i,11)*dnudz

  ! add vector producs
  grp%nodeHelpArray(i,1) = grp%nodeHelpArray(i,1) - fact*vecProd
 end do

 fact = oneOverReynoldsNumber/sigma

#ifdef PARALLEL

do kk=1,2
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfComSides
 else
  istart = grp%pdp%numberOfComSides + 1
  ifinish = grp%numberOfSides
 end if


 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL

  ind1 = grp%sideIndexArray(i,1)
  ind2 = grp%sideIndexArray(i,2)

  wx = grp%sideWeightsArray(i,1)
  wy = grp%sideWeightsArray(i,2)
  wz = grp%sideWeightsArray(i,3)

  dnudx1 = grp%nodeHelpArray(ind1,12)
  dnudy1 = grp%nodeHelpArray(ind1,13)
  dnudz1 = grp%nodeHelpArray(ind1,14)
  dnudx2 = grp%nodeHelpArray(ind2,12)
  dnudy2 = grp%nodeHelpArray(ind2,13)
  dnudz2 = grp%nodeHelpArray(ind2,14)

  visc1 = ivd%ReynoldsNumber*grp%laminarViscosity(ind1)/grp%u(ind1,1)+grp%u(ind1,6)
  visc2 = ivd%ReynoldsNumber*grp%laminarViscosity(ind2)/grp%u(ind2,1)+grp%u(ind2,6)

  fdiff1 = visc1*(wx*dnudx1+wy*dnudy1+wz*dnudz1)
  fdiff2 = visc2*(wx*dnudx2+wy*dnudy2+wz*dnudz2)

  ! add second derivatives
  grp%turbulenceDiffusionTerm(ind1) = grp%turbulenceDiffusionTerm(ind1) - fact*(fdiff1+fdiff2)
  grp%turbulenceDiffusionTerm(ind2) = grp%turbulenceDiffusionTerm(ind2) + fact*(fdiff1+fdiff2)
 end do

 ! boundary faces

#ifdef PARALLEL
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfBoundaryComSides
 else
  istart = grp%pdp%numberOfBoundaryComSides + 1
  ifinish = grp%brp%numberOfBoundarySides
 end if
 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL

  ind1 = grp%brp%sideIndexArray(i,1)
  ind2 = grp%brp%sideIndexArray(i,2)
  wx   = grp%brp%sideWeightsArray(i,1)
  wy   = grp%brp%sideWeightsArray(i,2)
  wz   = grp%brp%sideWeightsArray(i,3)

  dnudx1 = grp%nodeHelpArray(ind1,12)
  dnudy1 = grp%nodeHelpArray(ind1,13)
  dnudz1 = grp%nodeHelpArray(ind1,14)
  dnudx2 = grp%nodeHelpArray(ind2,12)
  dnudy2 = grp%nodeHelpArray(ind2,13)
  dnudz2 = grp%nodeHelpArray(ind2,14)

  visc1 = ivd%ReynoldsNumber*grp%laminarViscosity(ind1)/grp%u(ind1,1)+grp%u(ind1,6)
  visc2 = ivd%ReynoldsNumber*grp%laminarViscosity(ind2)/grp%u(ind2,1)+grp%u(ind2,6)

  fdiff1 = visc1*(wx*dnudx1+wy*dnudy1+wz*dnudz1)
  fdiff2 = visc2*(wx*dnudx2+wy*dnudy2+wz*dnudz2)

  ! adding boundary face contributions

  if(ivd%boundaryTerm==1) then 
   grp%turbulenceDiffusionTerm(ind1) = grp%turbulenceDiffusionTerm(ind1) - fact*(3.*fdiff1+fdiff2)
   grp%turbulenceDiffusionTerm(ind2) = grp%turbulenceDiffusionTerm(ind2) - fact*(fdiff1+3.*fdiff2)
  else
   grp%turbulenceDiffusionTerm(ind1) = grp%turbulenceDiffusionTerm(ind1) - fact*4.*fdiff1
   grp%turbulenceDiffusionTerm(ind2) = grp%turbulenceDiffusionTerm(ind2) - fact*4.*fdiff2
  end if
 end do
#ifdef PARALLEL

!    write(get_unit(),*)  'Starting first SA 2 comms',kk
 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses
   if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
    call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%turbulenceDiffusionTerm(:),grp%pdp%buffer,grp%pdp%sendFlag)
    call mp_sendv(grp%pdp%processorIDs,i,37,grp%pdp%sendFlag)
   end if
  end do
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,37,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
    call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                           grp%pdp%realType,grp%turbulenceDiffusionTerm(:),grp%pdp%buffer,grp%pdp%sendFlag)
  end if
 end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position =  1
   call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                     grp%pdp%realType,grp%turbulenceDiffusionTerm(:),grp%pdp%buffer,grp%pdp%sendFlag)
   call mp_sendv(grp%pdp%processorIDs,i,38,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,38,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
    call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%turbulenceDiffusionTerm(:),grp%pdp%buffer,grp%pdp%sendFlag)
  end if
 end do

#endif PARALLEL

 grp%turbulenceDiffusionTerm = grp%turbulenceDiffusionTerm + grp%nodeHelpArray(1:grp%numberOfNodes,1)

 end subroutine makeSADiffusion
!-----------------------------------------------------------------------
!START-TURB
subroutine makeKWDiffusion(grp,ivd)!MOD
IMPLICIT NONE
include 'mpif.h'
! diffusion terms for the k-omega turbulence model equation

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

 integer :: i,j,ind1,ind2,ind
 real :: ru1,ru2,rv1,rv2,rw1,rw2,fact,vt,visc1,visc2,fdiff1,fdiff2,fdiff3
 real :: turbVisc1,turbVisc2,tvx1,tvx2,tvy1,tvy2,tvz1,tvz2,dnudx,dnudy,dnudz
 real :: sigma,wx,wy,wz,dnudx1,dnudx2,dnudy1,dnudy2,dnudz1,dnudz2
 real :: oneOverReynoldsNumber,mut1,mut2

 integer :: kk,istart,ifinish,ij

 sigma = 2.
 oneOverReynoldsNumber = 1./ivd%ReynoldsNumber
 fact = oneOverReynoldsNumber/sigma

! Computation of derivatives for k and omega


!Loop for k(ij=6) and omega(ij=7)
do ij=6,7
 grp%nodeHelpArray(1:grp%numberOfNodes,1) = 0.0
 grp%nodeHelpArray(1:grp%numberOfNodes,9:14) = 0.0

#ifdef PARALLEL

do kk=1,2
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfComSides
 else
  istart = grp%pdp%numberOfComSides + 1
  ifinish = grp%numberOfSides
 end if


 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL

  ind1 = grp%sideIndexArray(i,1)
  ind2 = grp%sideIndexArray(i,2)

  turbVisc1 = grp%u(ind1,ij)
  turbVisc2 = grp%u(ind2,ij)

! side weights from preprocessor
  wx = grp%sideWeightsArray(i,1)
  wy = grp%sideWeightsArray(i,2)
  wz = grp%sideWeightsArray(i,3)


  fdiff1 = wx*(turbVisc1+turbVisc2)
  fdiff2 = wy*(turbVisc1+turbVisc2)
  fdiff3 = wz*(turbVisc1+turbVisc2)

  grp%nodeHelpArray(ind1,9) = grp%nodeHelpArray(ind1,9) + fdiff1
  grp%nodeHelpArray(ind2,9) = grp%nodeHelpArray(ind2,9) - fdiff1
  grp%nodeHelpArray(ind1,10) = grp%nodeHelpArray(ind1,10) + fdiff2
  grp%nodeHelpArray(ind2,10) = grp%nodeHelpArray(ind2,10) - fdiff2
  grp%nodeHelpArray(ind1,11) = grp%nodeHelpArray(ind1,11) + fdiff3
  grp%nodeHelpArray(ind2,11) = grp%nodeHelpArray(ind2,11) - fdiff3
 end do

! boundary faces


#ifdef PARALLEL
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfBoundaryComSides
 else
  istart = grp%pdp%numberOfBoundaryComSides + 1
  ifinish = grp%brp%numberOfBoundarySides
 end if
 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL

  ind1 = grp%brp%sideIndexArray(i,1)
  ind2 = grp%brp%sideIndexArray(i,2)
  wx   = grp%brp%sideWeightsArray(i,1)
  wy   = grp%brp%sideWeightsArray(i,2)
  wz   = grp%brp%sideWeightsArray(i,3)

  turbVisc1 = grp%u(ind1,ij)
  turbVisc2 = grp%u(ind2,ij)

 ! add boundary face contributions

  if(ivd%boundaryTerm==1) then 
   tvx1 = (3.*turbVisc1 + turbVisc2)*wx
   tvx2 = (turbVisc1 + 3.*turbVisc2)*wx
   tvy1 = (3.*turbVisc1 + turbVisc2)*wy
   tvy2 = (turbVisc1 + 3.*turbVisc2)*wy
   tvz1 = (3.*turbVisc1 + turbVisc2)*wz
   tvz2 = (turbVisc1 + 3.*turbVisc2)*wz
 
   grp%nodeHelpArray(ind1,9) = grp%nodeHelpArray(ind1,9) + tvx1
   grp%nodeHelpArray(ind2,9) = grp%nodeHelpArray(ind2,9) + tvx2
   grp%nodeHelpArray(ind1,10) = grp%nodeHelpArray(ind1,10) + tvy1
   grp%nodeHelpArray(ind2,10) = grp%nodeHelpArray(ind2,10) + tvy2
   grp%nodeHelpArray(ind1,11) = grp%nodeHelpArray(ind1,11) + tvz1
   grp%nodeHelpArray(ind2,11) = grp%nodeHelpArray(ind2,11) + tvz2
  else
   grp%nodeHelpArray(ind1,9) = grp%nodeHelpArray(ind1,9) + 4.*turbVisc1*wx
   grp%nodeHelpArray(ind2,9) = grp%nodeHelpArray(ind2,9) + 4.*turbVisc2*wx
   grp%nodeHelpArray(ind1,10) = grp%nodeHelpArray(ind1,10) + 4.*turbVisc1*wy
   grp%nodeHelpArray(ind2,10) = grp%nodeHelpArray(ind2,10) + 4.*turbvisc2*wy
   grp%nodeHelpArray(ind1,11) = grp%nodeHelpArray(ind1,11) + 4.*turbvisc1*wz
   grp%nodeHelpArray(ind2,11) = grp%nodeHelpArray(ind2,11) + 4.*turbvisc2*wz
  end if
 end do
#ifdef PARALLEL

!    write(get_unit(),*)  'Starting first KW comms',kk
 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses
   if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
    do j=9,11
     call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
    call mp_sendv(grp%pdp%processorIDs,i,35,grp%pdp%sendFlag)
   end if
  end do
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,35,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1236 continue

 do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
!  grp%pdp%sendFlag = 0
!  call mp_check_mesg_arrived(grp%pdp%processorIDs,i,35,grp%pdp%sendFlag,bufferzz,bufferSize)
!  if (grp%pdp%sendFlag.gt.0) then
!   grp%pdp%receivedMessageFromProcess(i) = 1
!   position = 1
    do j=9,11
     call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                            grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
!  end if 
  end if
 end do

!do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1236
!end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

!  write(get_unit(),*)  'Starting second KW comms'
 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position = 1
   do j=9,11
    call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                      grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
   end do
   call mp_sendv(grp%pdp%processorIDs,i,36,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,36,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1336 continue

 do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
!  grp%pdp%sendFlag = 0
!  call mp_check_mesg_arrived(grp%pdp%processorIDs,i,36,grp%pdp%sendFlag,bufferzz,bufferSize)
!  if (grp%pdp%sendFlag.gt.0) then
!   grp%pdp%receivedMessageFromProcess(i) = 1
!   position = 1
    do j=9,11
     call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
!  end if
  end if
 end do

!do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1336
!end do
#endif PARALLEL


 do i=1,grp%numberOfNodes
  vt = 1./grp%nodeVolume(i)
  dnudx = grp%nodeHelpArray(i,9)*vt
  dnudy = grp%nodeHelpArray(i,10)*vt
  dnudz = grp%nodeHelpArray(i,11)*vt
  grp%nodeHelpArray(i,12) = dnudx
  grp%nodeHelpArray(i,13) = dnudy
  grp%nodeHelpArray(i,14) = dnudz
 end do


#ifdef PARALLEL

do kk=1,2
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfComSides
 else
  istart = grp%pdp%numberOfComSides + 1
  ifinish = grp%numberOfSides
 end if


 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL

  ind1 = grp%sideIndexArray(i,1)
  ind2 = grp%sideIndexArray(i,2)

  wx = grp%sideWeightsArray(i,1)
  wy = grp%sideWeightsArray(i,2)
  wz = grp%sideWeightsArray(i,3)

  dnudx1 = grp%nodeHelpArray(ind1,12)
  dnudy1 = grp%nodeHelpArray(ind1,13)
  dnudz1 = grp%nodeHelpArray(ind1,14)
  dnudx2 = grp%nodeHelpArray(ind2,12)
  dnudy2 = grp%nodeHelpArray(ind2,13)
  dnudz2 = grp%nodeHelpArray(ind2,14)

  mut1=ivd%ReynoldsNumber*grp%u(ind1,1)*grp%u(ind1,6)/grp%u(ind1,7)
  mut2=ivd%ReynoldsNumber*grp%u(ind2,1)*grp%u(ind2,6)/grp%u(ind2,7)

  mut1=max(0.0,mut1)
  mut1=min(mut1,ivd%maxTurbulenceValueMU)
  mut2=max(0.0,mut2)
  mut2=min(mut2,ivd%maxTurbulenceValueMU)

  visc1 = ivd%ReynoldsNumber*sigma*grp%laminarViscosity(ind1)+mut1
  visc2 = ivd%ReynoldsNumber*sigma*grp%laminarViscosity(ind2)+mut2

  fdiff1 = visc1*(wx*dnudx1+wy*dnudy1+wz*dnudz1)
  fdiff2 = visc2*(wx*dnudx2+wy*dnudy2+wz*dnudz2)

  ! add second derivatives
  grp%nodeHelpArray(ind1,1) = grp%nodeHelpArray(ind1,1) + fact*(fdiff1+fdiff2)
  grp%nodeHelpArray(ind2,1) = grp%nodeHelpArray(ind2,1) - fact*(fdiff1+fdiff2)
 end do


 ! boundary faces

#ifdef PARALLEL
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfBoundaryComSides
 else
  istart = grp%pdp%numberOfBoundaryComSides + 1
  ifinish = grp%brp%numberOfBoundarySides
 end if
 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL

  ind1 = grp%brp%sideIndexArray(i,1)
  ind2 = grp%brp%sideIndexArray(i,2)
  wx   = grp%brp%sideWeightsArray(i,1)
  wy   = grp%brp%sideWeightsArray(i,2)
  wz   = grp%brp%sideWeightsArray(i,3)

  dnudx1 = grp%nodeHelpArray(ind1,12)
  dnudy1 = grp%nodeHelpArray(ind1,13)
  dnudz1 = grp%nodeHelpArray(ind1,14)
  dnudx2 = grp%nodeHelpArray(ind2,12)
  dnudy2 = grp%nodeHelpArray(ind2,13)
  dnudz2 = grp%nodeHelpArray(ind2,14)

  mut1=ivd%ReynoldsNumber*grp%u(ind1,1)*grp%u(ind1,6)/grp%u(ind1,7)
  mut2=ivd%ReynoldsNumber*grp%u(ind2,1)*grp%u(ind2,6)/grp%u(ind2,7)

  mut1=max(0.0,mut1)
  mut1=min(mut1,ivd%maxTurbulenceValueMU)
  mut2=max(0.0,mut2)
  mut2=min(mut2,ivd%maxTurbulenceValueMU)

  visc1 = ivd%ReynoldsNumber*sigma*grp%laminarViscosity(ind1)+mut1
  visc2 = ivd%ReynoldsNumber*sigma*grp%laminarViscosity(ind2)+mut2

  fdiff1 = visc1*(wx*dnudx1+wy*dnudy1+wz*dnudz1)
  fdiff2 = visc2*(wx*dnudx2+wy*dnudy2+wz*dnudz2)

  ! adding boundary face contributions

  if(ivd%boundaryTerm==1) then 
   grp%nodeHelpArray(ind1,1) = grp%nodeHelpArray(ind1,1) + fact*(3.*fdiff1+fdiff2)
   grp%nodeHelpArray(ind2,1) = grp%nodeHelpArray(ind2,1) + fact*(fdiff1+3.*fdiff2)
  else
   grp%nodeHelpArray(ind1,1) = grp%nodeHelpArray(ind1,1) + fact*4.*fdiff1
   grp%nodeHelpArray(ind2,1) = grp%nodeHelpArray(ind2,1) + fact*4.*fdiff2
  end if
 end do
#ifdef PARALLEL

!    write(get_unit(),*)  'Starting first KW 2 comms',kk
 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses
   if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
    call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
    call mp_sendv(grp%pdp%processorIDs,i,37,grp%pdp%sendFlag)
   end if
  end do
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,37,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1237 continue

 do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
!  grp%pdp%sendFlag = 0
!  call mp_check_mesg_arrived(grp%pdp%processorIDs,i,37,grp%pdp%sendFlag,bufferzz,bufferSize)
!  if (grp%pdp%sendFlag.gt.0) then
!   grp%pdp%receivedMessageFromProcess(i) = 1
!   position = 1
    call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                           grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
!  end if
  end if
 end do

!do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1237
!end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

!  write(get_unit(),*)  'Starting second KW 2 comms'
 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position =  1
   call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                     grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
   call mp_sendv(grp%pdp%processorIDs,i,38,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,38,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1337 continue

 do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
!  grp%pdp%sendFlag = 0
!  call mp_check_mesg_arrived(grp%pdp%processorIDs,i,38,grp%pdp%sendFlag,bufferzz,bufferSize)
!  if (grp%pdp%sendFlag.gt.0) then
!   grp%pdp%receivedMessageFromProcess(i) = 1
!   position = 1
    call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
!  end if
  end if
 end do

!do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1337
!end do
#endif PARALLEL

 if (ij.eq.6) then
  grp%turbulenceDiffusionTermK = -grp%nodeHelpArray(1:grp%numberOfNodes,1)
 else
  grp%turbulenceDiffusionTermW = -grp%nodeHelpArray(1:grp%numberOfNodes,1)
 end if
 end do !end loop for k and omega

 end subroutine makeKWDiffusion
!END-TURB
!-----------------------------------------------------------------------
!START-TURB
subroutine makeSSTDiffusion(grp,ivd)!SST
IMPLICIT NONE
include 'mpif.h'
! diffusion terms for the SST turbulence model equation

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

 integer :: i,j,ind1,ind2,ind
 real :: ru1,ru2,rv1,rv2,rw1,rw2,fact,vt,visc1,visc2,fdiff1,fdiff2,fdiff3
 real :: turbVisc1,turbVisc2,tvx1,tvx2,tvy1,tvy2,tvz1,tvz2,dnudx,dnudy,dnudz
 real :: wx,wy,wz,dnudx1,dnudx2,dnudy1,dnudy2,dnudz1,dnudz2
 real :: oneOverReynoldsNumber,mut(2),a1,nut
 real :: sigma1,sigma2,sigmak1,sigmak2,sigmaw1,sigmaw2

 integer :: kk,istart,ifinish,ij,ijk,in

 oneOverReynoldsNumber = 1./ivd%ReynoldsNumber
 fact = oneOverReynoldsNumber
 a1 = 0.31

! Computation of derivatives for k and omega

!Loop for k(ij=6) and omega(ij=7)
do ij=6,7

 grp%nodeHelpArray(1:grp%numberOfNodes,1) = 0.0
 grp%nodeHelpArray(1:grp%numberOfNodes,9:14) = 0.0

#ifdef PARALLEL

do kk=1,2
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfComSides
 else
  istart = grp%pdp%numberOfComSides + 1
  ifinish = grp%numberOfSides
 end if


 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL

  ind1 = grp%sideIndexArray(i,1)
  ind2 = grp%sideIndexArray(i,2)

  turbVisc1 = grp%u(ind1,ij)
  turbVisc2 = grp%u(ind2,ij)

! side weights from preprocessor
  wx = grp%sideWeightsArray(i,1)
  wy = grp%sideWeightsArray(i,2)
  wz = grp%sideWeightsArray(i,3)


  fdiff1 = wx*(turbVisc1+turbVisc2)
  fdiff2 = wy*(turbVisc1+turbVisc2)
  fdiff3 = wz*(turbVisc1+turbVisc2)

  grp%nodeHelpArray(ind1,9) = grp%nodeHelpArray(ind1,9) + fdiff1
  grp%nodeHelpArray(ind2,9) = grp%nodeHelpArray(ind2,9) - fdiff1
  grp%nodeHelpArray(ind1,10) = grp%nodeHelpArray(ind1,10) + fdiff2
  grp%nodeHelpArray(ind2,10) = grp%nodeHelpArray(ind2,10) - fdiff2
  grp%nodeHelpArray(ind1,11) = grp%nodeHelpArray(ind1,11) + fdiff3
  grp%nodeHelpArray(ind2,11) = grp%nodeHelpArray(ind2,11) - fdiff3
 end do

! boundary faces


#ifdef PARALLEL
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfBoundaryComSides
 else
  istart = grp%pdp%numberOfBoundaryComSides + 1
  ifinish = grp%brp%numberOfBoundarySides
 end if
 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL

  ind1 = grp%brp%sideIndexArray(i,1)
  ind2 = grp%brp%sideIndexArray(i,2)
  wx   = grp%brp%sideWeightsArray(i,1)
  wy   = grp%brp%sideWeightsArray(i,2)
  wz   = grp%brp%sideWeightsArray(i,3)

  turbVisc1 = grp%u(ind1,ij)
  turbVisc2 = grp%u(ind2,ij)

 ! add boundary face contributions

  if(ivd%boundaryTerm==1) then 
   tvx1 = (3.*turbVisc1 + turbVisc2)*wx
   tvx2 = (turbVisc1 + 3.*turbVisc2)*wx
   tvy1 = (3.*turbVisc1 + turbVisc2)*wy
   tvy2 = (turbVisc1 + 3.*turbVisc2)*wy
   tvz1 = (3.*turbVisc1 + turbVisc2)*wz
   tvz2 = (turbVisc1 + 3.*turbVisc2)*wz
 
   grp%nodeHelpArray(ind1,9) = grp%nodeHelpArray(ind1,9) + tvx1
   grp%nodeHelpArray(ind2,9) = grp%nodeHelpArray(ind2,9) + tvx2
   grp%nodeHelpArray(ind1,10) = grp%nodeHelpArray(ind1,10) + tvy1
   grp%nodeHelpArray(ind2,10) = grp%nodeHelpArray(ind2,10) + tvy2
   grp%nodeHelpArray(ind1,11) = grp%nodeHelpArray(ind1,11) + tvz1
   grp%nodeHelpArray(ind2,11) = grp%nodeHelpArray(ind2,11) + tvz2
  else
   grp%nodeHelpArray(ind1,9) = grp%nodeHelpArray(ind1,9) + 4.*turbVisc1*wx
   grp%nodeHelpArray(ind2,9) = grp%nodeHelpArray(ind2,9) + 4.*turbVisc2*wx
   grp%nodeHelpArray(ind1,10) = grp%nodeHelpArray(ind1,10) + 4.*turbVisc1*wy
   grp%nodeHelpArray(ind2,10) = grp%nodeHelpArray(ind2,10) + 4.*turbvisc2*wy
   grp%nodeHelpArray(ind1,11) = grp%nodeHelpArray(ind1,11) + 4.*turbvisc1*wz
   grp%nodeHelpArray(ind2,11) = grp%nodeHelpArray(ind2,11) + 4.*turbvisc2*wz
  end if
 end do
#ifdef PARALLEL

!    write(get_unit(),*)  'Starting first KW comms',kk
 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses
   if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
    do j=9,11
     call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
    call mp_sendv(grp%pdp%processorIDs,i,35,grp%pdp%sendFlag)
   end if
  end do
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,35,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1236 continue

 do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
!  grp%pdp%sendFlag = 0
!  call mp_check_mesg_arrived(grp%pdp%processorIDs,i,35,grp%pdp%sendFlag,bufferzz,bufferSize)
!  if (grp%pdp%sendFlag.gt.0) then
!   grp%pdp%receivedMessageFromProcess(i) = 1
!   position = 1
    do j=9,11
     call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                            grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
!  end if 
  end if
 end do

!do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1236
!end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

!  write(get_unit(),*)  'Starting second KW comms'
 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position = 1
   do j=9,11
    call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                      grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
   end do
   call mp_sendv(grp%pdp%processorIDs,i,36,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,36,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1336 continue

 do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
!  grp%pdp%sendFlag = 0
!  call mp_check_mesg_arrived(grp%pdp%processorIDs,i,36,grp%pdp%sendFlag,bufferzz,bufferSize)
!  if (grp%pdp%sendFlag.gt.0) then
!   grp%pdp%receivedMessageFromProcess(i) = 1
!   position = 1
    do j=9,11
     call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
!  end if
  end if
 end do

!do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1336
!end do
#endif PARALLEL


 do i=1,grp%numberOfNodes
  vt = 1./grp%nodeVolume(i)
  dnudx = grp%nodeHelpArray(i,9)*vt
  dnudy = grp%nodeHelpArray(i,10)*vt
  dnudz = grp%nodeHelpArray(i,11)*vt
  grp%nodeHelpArray(i,12) = dnudx
  grp%nodeHelpArray(i,13) = dnudy
  grp%nodeHelpArray(i,14) = dnudz
 end do


#ifdef PARALLEL

do kk=1,2
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfComSides
 else
  istart = grp%pdp%numberOfComSides + 1
  ifinish = grp%numberOfSides
 end if


 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL

  ind1 = grp%sideIndexArray(i,1)
  ind2 = grp%sideIndexArray(i,2)

  wx = grp%sideWeightsArray(i,1)
  wy = grp%sideWeightsArray(i,2)
  wz = grp%sideWeightsArray(i,3)

  dnudx1 = grp%nodeHelpArray(ind1,12)
  dnudy1 = grp%nodeHelpArray(ind1,13)
  dnudz1 = grp%nodeHelpArray(ind1,14)
  dnudx2 = grp%nodeHelpArray(ind2,12)
  dnudy2 = grp%nodeHelpArray(ind2,13)
  dnudz2 = grp%nodeHelpArray(ind2,14)

 ! Eddy viscosity limiter for edge points
 do ijk=1,2
  if(ijk .eq. 1) then
   in=ind1
  else
   in=ind2
  end if
   nut=a1*grp%u(in,6)/(max(a1*grp%u(in,7),grp%imsr(in)*grp%F2(in)))
   mut(ijk)=ivd%ReynoldsNumber*grp%u(in,1)*nut
   mut(ijk)=max(0.0,mut(ijk))
   mut(ijk)=min(mut(ijk),ivd%maxTurbulenceValueMU)
 end do

  if(ij==6)then !k
   sigma1 = grp%F1(ind1)*sigmak1+(1.-grp%F1(ind1))*sigmak2
   sigma2 = grp%F1(ind2)*sigmak1+(1.-grp%F1(ind2))*sigmak2
   else !omega
   sigma1 = grp%F1(ind1)*sigmaw1+(1.-grp%F1(ind1))*sigmaw2
   sigma2 = grp%F1(ind2)*sigmaw1+(1.-grp%F1(ind2))*sigmaw2
  end if

  visc1 = ivd%ReynoldsNumber*grp%laminarViscosity(ind1)+sigma1*mut(1)
  visc2 = ivd%ReynoldsNumber*grp%laminarViscosity(ind2)+sigma2*mut(2)


  fdiff1 = visc1*(wx*dnudx1+wy*dnudy1+wz*dnudz1)
  fdiff2 = visc2*(wx*dnudx2+wy*dnudy2+wz*dnudz2)

  ! add second derivatives
  grp%nodeHelpArray(ind1,1) = grp%nodeHelpArray(ind1,1) + fact*(fdiff1+fdiff2)
  grp%nodeHelpArray(ind2,1) = grp%nodeHelpArray(ind2,1) - fact*(fdiff1+fdiff2)
 end do


 ! boundary faces

#ifdef PARALLEL
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfBoundaryComSides
 else
  istart = grp%pdp%numberOfBoundaryComSides + 1
  ifinish = grp%brp%numberOfBoundarySides
 end if
 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL

  ind1 = grp%brp%sideIndexArray(i,1)
  ind2 = grp%brp%sideIndexArray(i,2)
  wx   = grp%brp%sideWeightsArray(i,1)
  wy   = grp%brp%sideWeightsArray(i,2)
  wz   = grp%brp%sideWeightsArray(i,3)

  dnudx1 = grp%nodeHelpArray(ind1,12)
  dnudy1 = grp%nodeHelpArray(ind1,13)
  dnudz1 = grp%nodeHelpArray(ind1,14)
  dnudx2 = grp%nodeHelpArray(ind2,12)
  dnudy2 = grp%nodeHelpArray(ind2,13)
  dnudz2 = grp%nodeHelpArray(ind2,14)


 ! Eddy viscosity limiter for edge points
 do ijk=1,2
  if(ijk .eq. 1) then
   in=ind1
  else
   in=ind2
  end if
   nut=a1*grp%u(in,6)/(max(a1*grp%u(in,7),grp%imsr(in)*grp%F2(in)))
   mut(ijk)=ivd%ReynoldsNumber*grp%u(in,1)*nut
   mut(ijk)=max(0.0,mut(ijk))
   mut(ijk)=min(mut(ijk),ivd%maxTurbulenceValueMU)
 end do

  if(ij==6)then !k
   sigma1 = grp%F1(ind1)*sigmak1+(1.-grp%F1(ind1))*sigmak2
   sigma2 = grp%F1(ind2)*sigmak1+(1.-grp%F1(ind2))*sigmak2
   else !omega
   sigma1 = grp%F1(ind1)*sigmaw1+(1.-grp%F1(ind1))*sigmaw2
   sigma2 = grp%F1(ind2)*sigmaw1+(1.-grp%F1(ind2))*sigmaw2
  end if

  visc1 = ivd%ReynoldsNumber*grp%laminarViscosity(ind1)+sigma1*mut(1)
  visc2 = ivd%ReynoldsNumber*grp%laminarViscosity(ind2)+sigma2*mut(2)

  fdiff1 = visc1*(wx*dnudx1+wy*dnudy1+wz*dnudz1)
  fdiff2 = visc2*(wx*dnudx2+wy*dnudy2+wz*dnudz2)

  ! adding boundary face contributions

  if(ivd%boundaryTerm==1) then 
   grp%nodeHelpArray(ind1,1) = grp%nodeHelpArray(ind1,1) + fact*(3.*fdiff1+fdiff2)
   grp%nodeHelpArray(ind2,1) = grp%nodeHelpArray(ind2,1) + fact*(fdiff1+3.*fdiff2)
  else
   grp%nodeHelpArray(ind1,1) = grp%nodeHelpArray(ind1,1) + fact*4.*fdiff1
   grp%nodeHelpArray(ind2,1) = grp%nodeHelpArray(ind2,1) + fact*4.*fdiff2
  end if
 end do
#ifdef PARALLEL

!    write(get_unit(),*)  'Starting first KW 2 comms',kk
 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses
   if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
    call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
    call mp_sendv(grp%pdp%processorIDs,i,37,grp%pdp%sendFlag)
   end if
  end do
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,37,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1237 continue

 do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
!  grp%pdp%sendFlag = 0
!  call mp_check_mesg_arrived(grp%pdp%processorIDs,i,37,grp%pdp%sendFlag,bufferzz,bufferSize)
!  if (grp%pdp%sendFlag.gt.0) then
!   grp%pdp%receivedMessageFromProcess(i) = 1
!   position = 1
    call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                           grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
!  end if
  end if
 end do

!do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1237
!end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

!  write(get_unit(),*)  'Starting second KW 2 comms'
 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position =  1
   call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                     grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
   call mp_sendv(grp%pdp%processorIDs,i,38,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,38,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1337 continue

 do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
!  grp%pdp%sendFlag = 0
!  call mp_check_mesg_arrived(grp%pdp%processorIDs,i,38,grp%pdp%sendFlag,bufferzz,bufferSize)
!  if (grp%pdp%sendFlag.gt.0) then
!   grp%pdp%receivedMessageFromProcess(i) = 1
!   position = 1
    call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
!  end if
  end if
 end do

!do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1337
!end do
#endif PARALLEL

 if (ij.eq.6) then !k
  grp%turbulenceDiffusionTermK = - grp%nodeHelpArray(1:grp%numberOfNodes,1)
 else !w
  grp%turbulenceDiffusionTermW = grp%turbulenceDiffusionTermW - grp%nodeHelpArray(1:grp%numberOfNodes,1)
 end if
 end do !end loop for k and omega

 end subroutine makeSSTDiffusion
!END-TURB
!-----------------------------------------------------------------------
!START-TURB
subroutine makeSSTdkdw(grp,ivd)!SST
 IMPLICIT NONE
 include 'mpif.h'
! diffusion terms for the k-omega turbulence model equation

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

 integer :: i,j,ind1,ind2,ind
 real :: ru1,ru2,rv1,rv2,rw1,rw2,vt,visc1,visc2,fdiff1,fdiff2,fdiff3
 real :: turbVisc1,turbVisc2,tvx1,tvx2,tvy1,tvy2,tvz1,tvz2,dnudx,dnudy,dnudz
 real :: wx,wy,wz,sigmaw2,CDkw,arg1,arg2,arga,argb,argc
 real :: sita2,kappa,delta,Cs,L,Lvk,betast,term2,modlap,beta1,beta2,QSAS
 real :: sigmatheta,sigmaw1,gama1,gama2,beta,gama,FSAS,C

 integer :: kk,istart,ifinish,ij

 sita2 = 3.51 !Version 2007 !1.755(Version 2005)
 FSAS = 1.0 !Version 2007 !1.25(Version 2005)
 C = 2.0 !Version 2007
 kappa = 0.41
 sigmatheta = 2./3.
 Cs = 0.11 !Smagorinsky constant
 betast = 0.09
 beta1=0.075
 beta2=0.0828
 sigmaw1=0.5
 sigmaw2=0.856
 gama1=beta1/betast-sigmaw1*kappa**2/sqrt(betast)
 gama2=beta2/betast-sigmaw2*kappa**2/sqrt(betast)

 grp%nodeHelpArray(1:grp%numberOfNodes,12:17) = 0.0

! Computation of derivatives for k and omega
! Loop for k(ij=6) and omega(ij=7)
do ij=6,7
grp%nodeHelpArray(1:grp%numberOfNodes,9:11) = 0.0 

#ifdef PARALLEL

do kk=1,2
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfComSides
 else
  istart = grp%pdp%numberOfComSides + 1
  ifinish = grp%numberOfSides
 end if


 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL

  ind1 = grp%sideIndexArray(i,1)
  ind2 = grp%sideIndexArray(i,2)

  turbVisc1 = grp%u(ind1,ij)
  turbVisc2 = grp%u(ind2,ij)

! side weights from preprocessor
  wx = grp%sideWeightsArray(i,1)
  wy = grp%sideWeightsArray(i,2)
  wz = grp%sideWeightsArray(i,3)


  fdiff1 = wx*(turbVisc1+turbVisc2)
  fdiff2 = wy*(turbVisc1+turbVisc2)
  fdiff3 = wz*(turbVisc1+turbVisc2)

  grp%nodeHelpArray(ind1,9) = grp%nodeHelpArray(ind1,9) + fdiff1
  grp%nodeHelpArray(ind2,9) = grp%nodeHelpArray(ind2,9) - fdiff1
  grp%nodeHelpArray(ind1,10) = grp%nodeHelpArray(ind1,10) + fdiff2
  grp%nodeHelpArray(ind2,10) = grp%nodeHelpArray(ind2,10) - fdiff2
  grp%nodeHelpArray(ind1,11) = grp%nodeHelpArray(ind1,11) + fdiff3
  grp%nodeHelpArray(ind2,11) = grp%nodeHelpArray(ind2,11) - fdiff3
 end do

! boundary faces

#ifdef PARALLEL
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfBoundaryComSides
 else
  istart = grp%pdp%numberOfBoundaryComSides + 1
  ifinish = grp%brp%numberOfBoundarySides
 end if
 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL

  ind1 = grp%brp%sideIndexArray(i,1)
  ind2 = grp%brp%sideIndexArray(i,2)
  wx   = grp%brp%sideWeightsArray(i,1)
  wy   = grp%brp%sideWeightsArray(i,2)
  wz   = grp%brp%sideWeightsArray(i,3)

  turbVisc1 = grp%u(ind1,ij)
  turbVisc2 = grp%u(ind2,ij)

 ! add boundary face contributions

  if(ivd%boundaryTerm==1) then 
   tvx1 = (3.*turbVisc1 + turbVisc2)*wx
   tvx2 = (turbVisc1 + 3.*turbVisc2)*wx
   tvy1 = (3.*turbVisc1 + turbVisc2)*wy
   tvy2 = (turbVisc1 + 3.*turbVisc2)*wy
   tvz1 = (3.*turbVisc1 + turbVisc2)*wz
   tvz2 = (turbVisc1 + 3.*turbVisc2)*wz
 
   grp%nodeHelpArray(ind1,9) = grp%nodeHelpArray(ind1,9) + tvx1
   grp%nodeHelpArray(ind2,9) = grp%nodeHelpArray(ind2,9) + tvx2
   grp%nodeHelpArray(ind1,10) = grp%nodeHelpArray(ind1,10) + tvy1
   grp%nodeHelpArray(ind2,10) = grp%nodeHelpArray(ind2,10) + tvy2
   grp%nodeHelpArray(ind1,11) = grp%nodeHelpArray(ind1,11) + tvz1
   grp%nodeHelpArray(ind2,11) = grp%nodeHelpArray(ind2,11) + tvz2
  else
   grp%nodeHelpArray(ind1,9) = grp%nodeHelpArray(ind1,9) + 4.*turbVisc1*wx
   grp%nodeHelpArray(ind2,9) = grp%nodeHelpArray(ind2,9) + 4.*turbVisc2*wx
   grp%nodeHelpArray(ind1,10) = grp%nodeHelpArray(ind1,10) + 4.*turbVisc1*wy
   grp%nodeHelpArray(ind2,10) = grp%nodeHelpArray(ind2,10) + 4.*turbvisc2*wy
   grp%nodeHelpArray(ind1,11) = grp%nodeHelpArray(ind1,11) + 4.*turbvisc1*wz
   grp%nodeHelpArray(ind2,11) = grp%nodeHelpArray(ind2,11) + 4.*turbvisc2*wz
  end if
 end do
#ifdef PARALLEL

!    write(get_unit(),*)  'Starting first dkdxj or dw/dxj comms',kk
 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses
   if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
    do j=9,11
     call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
    call mp_sendv(grp%pdp%processorIDs,i,35,grp%pdp%sendFlag)
   end if
  end do
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,35,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1236 continue

 do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
!  grp%pdp%sendFlag = 0
!  call mp_check_mesg_arrived(grp%pdp%processorIDs,i,35,grp%pdp%sendFlag,bufferzz,bufferSize)
!  if (grp%pdp%sendFlag.gt.0) then
!   grp%pdp%receivedMessageFromProcess(i) = 1
!   position = 1
    do j=9,11
     call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                            grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
!  end if 
  end if
 end do

!do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1236
!end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

!  write(get_unit(),*)  'Starting second KW comms'
 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position = 1
   do j=9,11
    call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                      grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
   end do
   call mp_sendv(grp%pdp%processorIDs,i,36,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,36,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1336 continue

 do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
!  grp%pdp%sendFlag = 0
!  call mp_check_mesg_arrived(grp%pdp%processorIDs,i,36,grp%pdp%sendFlag,bufferzz,bufferSize)
!  if (grp%pdp%sendFlag.gt.0) then
!   grp%pdp%receivedMessageFromProcess(i) = 1
!   position = 1
    do j=9,11
     call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
!  end if
  end if
 end do

!do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1336
!end do
#endif PARALLEL

 do i=1,grp%numberOfNodes
  vt = 1./grp%nodeVolume(i)
  if (ij.eq.6) then !derivatives of k
   grp%nodeHelpArray(i,12) = grp%nodeHelpArray(i,9)*vt
   grp%nodeHelpArray(i,13) = grp%nodeHelpArray(i,10)*vt
   grp%nodeHelpArray(i,14) = grp%nodeHelpArray(i,11)*vt
  else !derivatives of w
   grp%nodeHelpArray(i,15) = grp%nodeHelpArray(i,9)*vt
   grp%nodeHelpArray(i,16) = grp%nodeHelpArray(i,10)*vt
   grp%nodeHelpArray(i,17) = grp%nodeHelpArray(i,11)*vt
  end if
 end do

end do !end loop for k and omega

 do i=1,grp%numberOfNodes
!Computing the product dk/dxj*dw/dxj at each point 
  grp%dkdw(i)=grp%nodeHelpArray(i,12)*grp%nodeHelpArray(i,15)+grp%nodeHelpArray(i,13)*grp%nodeHelpArray(i,16)+grp%nodeHelpArray(i,14)*grp%nodeHelpArray(i,17)   

arga=0.0
argb=0.0
argc=0.0

!Calculating the functions F1 and F2
 CDkw=max(2*grp%u(i,1)*sigmaw2*grp%dkdw(i)/grp%u(i,7),1E-10)
!Calculating the argument arg1
if(grp%wallDistance(i) .eq. 0.0 .or. grp%u(i,7) .eq. 0.0) then
 grp%F1(i)= 1.0
 grp%F2(i)= 1.0
else
 arga=sqrt(abs(grp%u(i,6)))/(0.09*grp%u(i,7)*grp%wallDistance(i))
 argb=500.*grp%laminarViscosity(i)/grp%u(i,1)/grp%wallDistance(i)**2/grp%u(i,7)
 argc=4.*grp%u(i,1)*sigmaw2*grp%u(i,6)/CDkw/grp%wallDistance(i)**2
 arg1=min(max(arga,argb),argc)
 grp%F1(i)=tanh(arg1**4)
 arg2=max(2.*arga,argb)
 grp%F2(i)=tanh(arg2**2)
end if 

!Computing the Qsas term in the SST-SAS model
 if(ivd%SAS) then
  if(grp%u(i,7) .ne. 0.0) then
   L = sqrt(abs(grp%u(i,6)))/(betast**(0.25)*grp%u(i,7)) !length scale
  else
   L = 0.0
  end if
 delta = (grp%nodeVolume(i))**(1./3.)
 term2=C*2.*grp%u(i,6)/sigmatheta*max((grp%nodeHelpArray(i,15)**2.+grp%nodeHelpArray(i,16)**2.+grp%nodeHelpArray(i,17)**2.)/grp%u(i,7)**2.,(grp%nodeHelpArray(i,12)**2.+grp%nodeHelpArray(i,13)**2.+grp%nodeHelpArray(i,14)**2.)/grp%u(i,6)**2.)
 modlap = sqrt(grp%D2U(1,i)**2.+ grp%D2U(2,i)**2.+ grp%D2U(3,i)**2.)
 beta=grp%F1(i)*beta1+(1.-grp%F1(i))*beta2 
 gama=grp%F1(i)*gama1+(1.-grp%F1(i))*gama2
 Lvk = max(kappa*grp%imsr(i)/modlap,Cs*sqrt(kappa*sita2/(beta/betast-gama))*delta)
 QSAS = grp%u(i,1)*FSAS*max(sita2*kappa*grp%imsr(i)**2.*(L/Lvk)**2.-term2,0.0) !Version 2007
  vt = grp%nodeVolume(i)
!For output purposes
 !! grp%sas(1,i) = vt*QSAS (now, inst. vorticity is saved in this array)
  grp%sas(2,i) = L
  grp%sas(3,i) = Lvk
  grp%sas(4,i) = grp%vorticity(i)**2.-grp%imsr(i)**2. 
!end output purposes
  grp%turbulenceDiffusionTermW(i) = - vt*QSAS
 else
  grp%turbulenceDiffusionTermW(i) = 0.0
 end if

 end do

 end subroutine makeSSTdkdw
!END-TURB
!-----------------------------------------------------------------------
!START-TURB
subroutine makeSASd2u(grp,ivd)
! makes the viscosity fluxes using the wide stencil 
IMPLICIT NONE
include 'mpif.h'

type(GridSolverData) :: grp
type(InputVariablesData) :: ivd

integer :: i1,i2,i,j,ip1,ip2,faceIndicator,ib,ip,ist,ien,ind
real :: r1,r2,u1,u2,v1,v2,t1,t2,wx,wy,wz,dux,dvx,dtx,duy,dvy,dty,oneOverReynoldsNumber
real :: dwx,dwy,dwz,dwdx,dwdy,dwdz,dTdz,tz,t61,t71,t81,t91,t62,t72,t82,t92,duz,dvz,dtz
real :: w1,w2,dudz,dvdz,u31,u32,dt6,dt7,dt8,dt9,dut3,f5,dut31,dut32
real :: mu,k,T,T0,f2,f3,f4,rx,ry,vt,t11,t21,t31,t41,t51,t12,t22,t32
real :: t42,t52,dudx,dvdx,dTdx,dudy,dvdy,dTdy,dt1,dt2,dt3,dt4,dt5,oneOverPrandtlNumber
real :: twoThirds,fourThirds,gammaOverGammaMinusOne,u11,u12,u21,u22,dut1,dut2
real :: wallTangent(3),wallNormal(3),tempConv,dut11,dut12,dut21,dut22,inflowTemp
real :: adTemp1,adTemp2,turbDiffusionCoefficient
real :: viscosityCoefficient,diffusionCoefficient,oneOverTurbPrandtlNumber
real :: uu1,uu2,duux,duuy,duuz,locu,dt,increment,sign
real :: onehalf,fournineth
integer :: noutput
real :: nut,a1
double precision :: turbViscosityCoefficient,Ksi

 integer :: kk,istart,ifinish

 integer :: mind
 real :: mres

  grp%nodeHelpArray(1:grp%numberOfNodes,1:9) = 0.0


#ifdef PARALLEL

  do kk=1,2
    if(kk==1) then
      istart = 1
      ifinish = grp%pdp%numberOfComSides
    else
      istart = grp%pdp%numberOfComSides + 1
      ifinish = grp%numberOfSides
    end if


    do i=istart,ifinish

#else PARALLEL

    do i=1,grp%numberOfSides

#endif PARALLEL

! indexes of nodes in side
      i1 = grp%sideIndexArray(i,1)
      i2 = grp%sideIndexArray(i,2)

! side weights from preprocessor 
      wx = grp%sideWeightsArray(i,1)  
      wy = grp%sideWeightsArray(i,2)  
      wz = grp%sideWeightsArray(i,3)  

      dux   = (grp%D1U(1,i1)+grp%D1U(1,i2))*wx
      duy   = (grp%D1U(2,i1)+grp%D1U(2,i2))*wy
      duz   = (grp%D1U(3,i1)+grp%D1U(3,i2))*wz
      dvx   = (grp%D1U(4,i1)+grp%D1U(4,i2))*wx
      dvy   = (grp%D1U(5,i1)+grp%D1U(5,i2))*wy
      dvz   = (grp%D1U(6,i1)+grp%D1U(6,i2))*wz
      dwx   = (grp%D1U(7,i1)+grp%D1U(7,i2))*wx
      dwy   = (grp%D1U(8,i1)+grp%D1U(8,i2))*wy
      dwz   = (grp%D1U(9,i1)+grp%D1U(9,i2))*wz

      grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + dux
      grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + duy
      grp%nodeHelpArray(i1,3) = grp%nodeHelpArray(i1,3) + duz
      grp%nodeHelpArray(i1,4) = grp%nodeHelpArray(i1,4) + dvx
      grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) + dvy
      grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) + dvz
      grp%nodeHelpArray(i1,7) = grp%nodeHelpArray(i1,7) + dwx
      grp%nodeHelpArray(i1,8) = grp%nodeHelpArray(i1,8) + dwy
      grp%nodeHelpArray(i1,9) = grp%nodeHelpArray(i1,9) + dwz
  
      grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) - dux
      grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) - duy
      grp%nodeHelpArray(i2,3) = grp%nodeHelpArray(i2,3) - duz
      grp%nodeHelpArray(i2,4) = grp%nodeHelpArray(i2,4) - dvx
      grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) - dvy
      grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) - dvz  
      grp%nodeHelpArray(i2,7) = grp%nodeHelpArray(i2,7) - dwx
      grp%nodeHelpArray(i2,8) = grp%nodeHelpArray(i2,8) - dwy
      grp%nodeHelpArray(i2,9) = grp%nodeHelpArray(i2,9) - dwz

    end do

#ifdef PARALLEL
    if(kk==1) then
      istart = 1
      ifinish = grp%pdp%numberOfBoundaryComSides
    else
      istart = grp%pdp%numberOfBoundaryComSides + 1
      ifinish = grp%brp%numberOfBoundarySides
    end if
    do i=istart,ifinish

#else PARALLEL

    do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL

      i1 = grp%brp%sideIndexArray(i,1)
      i2 = grp%brp%sideIndexArray(i,2)
      wx   = grp%brp%sideWeightsArray(i,1)
      wy   = grp%brp%sideWeightsArray(i,2) 
      wz   = grp%brp%sideWeightsArray(i,3) 

      dux   = (grp%D1U(1,i1)+grp%D1U(1,i2))
      duy   = (grp%D1U(2,i1)+grp%D1U(2,i2))
      duz   = (grp%D1U(3,i1)+grp%D1U(3,i2))
      dvx   = (grp%D1U(4,i1)+grp%D1U(4,i2))
      dvy   = (grp%D1U(5,i1)+grp%D1U(5,i2))
      dvz   = (grp%D1U(6,i1)+grp%D1U(6,i2))
      dwx   = (grp%D1U(7,i1)+grp%D1U(7,i2))
      dwy   = (grp%D1U(8,i1)+grp%D1U(8,i2))
      dwz   = (grp%D1U(9,i1)+grp%D1U(9,i2))
     

! adding boundary face contributions        

      grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + (2.*grp%D1U(1,i1)+dux)*wx
      grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + (2.*grp%D1U(2,i1)+duy)*wy
      grp%nodeHelpArray(i1,3) = grp%nodeHelpArray(i1,3) + (2.*grp%D1U(3,i1)+duz)*wz
      grp%nodeHelpArray(i1,4) = grp%nodeHelpArray(i1,4) + (2.*grp%D1U(4,i1)+dvx)*wx
      grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) + (2.*grp%D1U(5,i1)+dvy)*wy
      grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) + (2.*grp%D1U(6,i1)+dvz)*wz
      grp%nodeHelpArray(i1,7) = grp%nodeHelpArray(i1,7) + (2.*grp%D1U(7,i1)+dwx)*wx
      grp%nodeHelpArray(i1,8) = grp%nodeHelpArray(i1,8) + (2.*grp%D1U(8,i1)+dwy)*wy
      grp%nodeHelpArray(i1,9) = grp%nodeHelpArray(i1,9) + (2.*grp%D1U(9,i1)+dwz)*wz
           
      grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + (2.*grp%D1U(1,i2)+dux)*wx
      grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) + (2.*grp%D1U(2,i2)+duy)*wy
      grp%nodeHelpArray(i2,3) = grp%nodeHelpArray(i2,3) + (2.*grp%D1U(3,i2)+duz)*wz
      grp%nodeHelpArray(i2,4) = grp%nodeHelpArray(i2,4) + (2.*grp%D1U(4,i2)+dvx)*wx
      grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) + (2.*grp%D1U(5,i2)+dvy)*wy
      grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) + (2.*grp%D1U(6,i2)+dvz)*wz
      grp%nodeHelpArray(i2,7) = grp%nodeHelpArray(i2,7) + (2.*grp%D1U(7,i2)+dwx)*wx
      grp%nodeHelpArray(i2,8) = grp%nodeHelpArray(i2,8) + (2.*grp%D1U(8,i2)+dwy)*wy
      grp%nodeHelpArray(i2,9) = grp%nodeHelpArray(i2,9) + (2.*grp%D1U(9,i2)+dwz)*wz

    end do

#ifdef PARALLEL

    if (kk.eq.1) then
      do i=1,grp%pdp%numberOfProcesses
        if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
          call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
          position = 1
          do j=1,9
            call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                               grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
          end do
          call mp_sendv(grp%pdp%processorIDs,i,29,grp%pdp%sendFlag)
        end if
      end do
    end if
  end do

  call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,29,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

  do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
        do j=1,9
          call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                                 grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
        end do
    end if
  end do

  do i=1,grp%pdp%numberOfProcesses
    if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
      call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
      position = 1
      do j=1,9
        call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                          grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
      end do
      call mp_sendv(grp%pdp%processorIDs,i,30,grp%pdp%sendFlag)
    end if
  end do

  call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,30,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1333 continue

  do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
        do j = 1,9
          call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                             grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
        end do
    end if
  end do

#endif PARALLEL

do i=1,grp%numberOfNodes
  vt = 1./grp%nodeVolume(i)
   grp%D2U(1,i) = (grp%nodeHelpArray(i,1)+grp%nodeHelpArray(i,2)+grp%nodeHelpArray(i,3))*vt
   grp%D2U(2,i) = (grp%nodeHelpArray(i,4)+grp%nodeHelpArray(i,5)+grp%nodeHelpArray(i,6))*vt
   grp%D2U(3,i) = (grp%nodeHelpArray(i,7)+grp%nodeHelpArray(i,8)+grp%nodeHelpArray(i,9))*vt
end do

end subroutine makeSASd2u
!END-TURB
!-----------------------------------------------------------------------
!START-TURB
 subroutine makeSSTCrossDiff(grp,ivd)!SST
 IMPLICIT NONE
 include 'mpif.h'
! cross-diffusion term for the SST model equation

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

 real :: sigmaw2,vt
 integer :: i

 sigmaw2 = 0.856
 do i=1,grp%numberOfNodes
  vt = grp%nodeVolume(i)
  grp%turbulenceDiffusionTermW(i) = grp%turbulenceDiffusionTermW(i) - vt*2.*(1.-grp%F1(i))*sigmaw2*grp%u(i,1)/grp%u(i,7)*grp%dkdw(i)
 end do

 end subroutine makeSSTCrossDiff
!END-TURB
!-----------------------------------------------------------------------
 subroutine makeSASourceTerm(grp,ivd)
 ! source term for the Spalart-Allmaras turbulence model
 IMPLICIT NONE


 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

 integer :: i,j,ind
 double precision :: vt,mucurl,Ksi,ftt,fv1,fv2,fv3,Scurl,d,dt,KsiInThird,r
 double precision :: dSquared,muOverd,fw,deltau,deltauSquared,tripVorticity,tripDeltax
 double precision :: ft1,ft2,gt,cb1,cb2,kappa,kappaSquared,sigma,cw1,cw2,cw3,cv1,cv2
 double precision :: ct1,ct2,ct3,ct4,cv1Tripled,oneOverReynoldsNumber
 double precision :: radicand,g,cw3InSixth,oneOverSix,cFact

 ind = 26063

 oneOverReynoldsNumber = 1./ivd%ReynoldsNumber
 cb1 = 0.1355
 cb2 = 0.622
 kappa = 0.41
 kappaSquared = kappa*kappa
 sigma = 2./3.
 cw1 = cb1/kappaSquared + (1.0+cb2)/sigma
 cw2 = 0.3
 cw3 = 2.0
 cv1 = 7.1
 cv2 = 5.0
 ct1 = 1.0
 ct2 = 2.0
 ct3 = 1.1 ! old value
 ct4 = 2.0 ! old value

 cv1Tripled = cv1**3
 cw3InSixth =cw3**6
 oneOverSix = 1./6.
 cFact = exp(oneOverSix*log(1.0+cw3InSixth))

 do i=1,grp%numberOfNodes
  vt = grp%nodeVolume(i)
! mass matrices are lumped for source terms
  mucurl = grp%u(i,6)
  mucurl = max(mucurl,0.0) 
  if(grp%laminarViscosity(i)<0.1/ivd%ReynoldsNumber) then 
   grp%laminarViscosity(i) = 0.1*ivd%ReynoldsNumber
  end if
  Ksi = grp%u(i,1)*oneOverReynoldsNumber*mucurl/grp%laminarViscosity(i)

  Ksi = max(Ksi,0.001)
  ft2 = ct3*exp(-ct4*Ksi*Ksi)
  KsiInThird = Ksi**3
  d = grp%wallDistance(i)
  d = max(d,1.0e-10)

  muOverd = mucurl/d

  dSquared = d*d
  fv1 = KsiInThird/(KsiInThird+cv1Tripled)
  fv2 = 1.0-(Ksi/(1.0+Ksi*fv1)) ! old formulation
  fv3 = 1.0                        ! old formulation 

  Scurl = fv3*grp%vorticity(i) + oneOverReynoldsNumber*mucurl*fv2/(kappaSquared*dSquared)

  if(abs(Scurl)>1.0e-10) then
   r = oneOverReynoldsNumber*mucurl/(Scurl*kappaSquared*dSquared)
   if(r<0.0) r = 0.0
   if(r>100.0) r = 100.0
  else
   r = 100.0
  end if

  g = r + cw2*(r**6-r)
  radicand = (1.0+cw3InSixth)/(g**6+cw3InSixth)
  if(radicand>1.0e-10) then
   fw = g*(exp(oneOverSix*log(radicand)))
  else
   fw = cFact
  end if
! add source terms

  grp%rhs(i,6) = grp%rhs(i,6) - vt*cb1*(1.0-ft2)*Scurl*mucurl
  grp%rhs(i,6) = grp%rhs(i,6) + oneOverReynoldsNumber*vt*(cw1*fw-cb1*ft2/kappaSquared)*muOverd*muOverd
  grp%rhs(i,6) = grp%rhs(i,6) - vt*mucurl*grp%divergence(i)
 end do

 end subroutine makeSASourceTerm
!-----------------------------------------------------------------------
!START-TURB
subroutine makeKWSourceTerm(grp,ivd)!MOD
 ! source term for k in the k-omega turbulence model equation
 ! the viscosity fluxes are computed by using the wide stencil 
IMPLICIT NONE
include 'mpif.h'

type(GridSolverData) :: grp
type(InputVariablesData) :: ivd

integer :: i1,i2,i,j,ip1,ip2,faceIndicator,ib,ip,ist,ien,ind
real :: r1,r2,u1,u2,v1,v2,wx,wy,wz,oneOverReynoldsNumber
real :: dudx,dvdx,dwdx,dudy,dvdy,dwdy,dudz,dvdz,dwdz
real :: dux,dvx,dwx,duy,dvy,dwy,duz,dvz,dwz
real :: dudx1,dvdx1,dwdx1,dudy1,dvdy1,dwdy1,dudz1,dvdz1,dwdz1
real :: dudx2,dvdx2,dwdx2,dudy2,dvdy2,dwdy2,dudz2,dvdz2,dwdz2
real :: w1,w2,dt1,dt2,dt3,dt4,dt5,dt6,dut1,dut2,dut3
real :: mu,vt,t11,t21,t31,t41,t51,t61,t12,t22,t32,t42,t52,t62
real :: twoThirds,fourThirds
real :: prod,betast,beta,alfa,betast0,beta0,Mto,H
real :: dt,prodk,dissk,prodw,dissw
double precision :: turbViscosityCoefficient

 integer :: kk,istart,ifinish,ii,IAUX,IAREA

  grp%nodeHelpArray(1:grp%numberOfNodes,1:16) = 0.0

!Constant definition
 betast0 = 9./100. !Wilcox 2006 and 1988 (incompressible value)
! beta0 = 0.0708    !Wilcox 2006 (incompressible value)
 beta0 = 3./40.    !Wilcox 1988 (incompressible value)
! alfa = 13./25. !Wilcox 2006
 alfa = 5./9. !Wilcox 1988
 IAREA = 0



#ifdef PARALLEL

  do kk=1,2
    if(kk==1) then
      istart = 1
      ifinish = grp%pdp%numberOfComSides
    else
      istart = grp%pdp%numberOfComSides + 1
      ifinish = grp%numberOfSides
    end if


    do i=istart,ifinish

#else PARALLEL

    do i=1,grp%numberOfSides

#endif PARALLEL


! indexes of nodes in side
      i1 = grp%sideIndexArray(i,1)
      i2 = grp%sideIndexArray(i,2)

! side weights from preprocessor 
      wx = grp%sideWeightsArray(i,1)  
      wy = grp%sideWeightsArray(i,2)  
      wz = grp%sideWeightsArray(i,3)  

      r1   = grp%u(i1,1)                                  ! density
      u1   = grp%u(i1,2)/r1                               ! x - velocity
      v1   = grp%u(i1,3)/r1                               ! y - velocity
      w1   = grp%u(i1,4)/r1                               ! z - velocity  
  
      r2   = grp%u(i2,1)
      u2   = grp%u(i2,2)/r2
      v2   = grp%u(i2,3)/r2
      w2   = grp%u(i2,4)/r2

      dux   = (u1+u2)*wx
      dvx   = (v1+v2)*wx
      dwx   = (w1+w2)*wx
      duy   = (u1+u2)*wy
      dvy   = (v1+v2)*wy
      dwy   = (w1+w2)*wy
      duz   = (u1+u2)*wz
      dvz   = (v1+v2)*wz
      dwz   = (w1+w2)*wz


      grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + dux
      grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + dvx
      grp%nodeHelpArray(i1,3) = grp%nodeHelpArray(i1,3) + dwx
      grp%nodeHelpArray(i1,4) = grp%nodeHelpArray(i1,4) + duy
      grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) + dvy
      grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) + dwy
      grp%nodeHelpArray(i1,7) = grp%nodeHelpArray(i1,7) + duz
      grp%nodeHelpArray(i1,8) = grp%nodeHelpArray(i1,8) + dvz
      grp%nodeHelpArray(i1,9) = grp%nodeHelpArray(i1,9) + dwz

      grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) - dux
      grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) - dvx
      grp%nodeHelpArray(i2,3) = grp%nodeHelpArray(i2,3) - dwx
      grp%nodeHelpArray(i2,4) = grp%nodeHelpArray(i2,4) - duy
      grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) - dvy
      grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) - dwy
      grp%nodeHelpArray(i2,7) = grp%nodeHelpArray(i2,7) - duz
      grp%nodeHelpArray(i2,8) = grp%nodeHelpArray(i2,8) - dvz
      grp%nodeHelpArray(i2,9) = grp%nodeHelpArray(i2,9) - dwz
    end do

#ifdef PARALLEL
    if(kk==1) then
      istart = 1
      ifinish = grp%pdp%numberOfBoundaryComSides
    else
      istart = grp%pdp%numberOfBoundaryComSides + 1
      ifinish = grp%brp%numberOfBoundarySides
    end if
    do i=istart,ifinish

#else PARALLEL

    do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL

      i1 = grp%brp%sideIndexArray(i,1)
      i2 = grp%brp%sideIndexArray(i,2)
      wx   = grp%brp%sideWeightsArray(i,1)
      wy   = grp%brp%sideWeightsArray(i,2) 
      wz   = grp%brp%sideWeightsArray(i,3) 

      r1   = grp%u(i1,1)
      u1   = grp%u(i1,2)/r1
      v1   = grp%u(i1,3)/r1
      w1   = grp%u(i1,4)/r1
      r2   = grp%u(i2,1)
      u2   = grp%u(i2,2)/r2
      v2   = grp%u(i2,3)/r2
      w2   = grp%u(i2,4)/r2
      dux   = u1+u2
      dvx   = v1+v2
      dwx   = w1+w2
      duy   = u1+u2
      dvy   = v1+v2
      dwy   = w1+w2
      duz   = u1+u2
      dvz   = v1+v2
      dwz   = w1+w2

! adding boundary face contributions        

      grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + (2.*u1+dux)*wx
      grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + (2.*v1+dvx)*wx
      grp%nodeHelpArray(i1,3) = grp%nodeHelpArray(i1,3) + (2.*w1+dwx)*wx
      grp%nodeHelpArray(i1,4) = grp%nodeHelpArray(i1,4) + (2.*u1+duy)*wy
      grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) + (2.*v1+dvy)*wy
      grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) + (2.*w1+dwy)*wy
      grp%nodeHelpArray(i1,7) = grp%nodeHelpArray(i1,7) + (2.*u1+duz)*wz
      grp%nodeHelpArray(i1,8) = grp%nodeHelpArray(i1,8) + (2.*v1+dvz)*wz
      grp%nodeHelpArray(i1,9) = grp%nodeHelpArray(i1,9) + (2.*w1+dwz)*wz
        
      grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + (2.*u2+dux)*wx
      grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) + (2.*v2+dvx)*wx
      grp%nodeHelpArray(i2,3) = grp%nodeHelpArray(i2,3) + (2.*w2+dwx)*wx
      grp%nodeHelpArray(i2,4) = grp%nodeHelpArray(i2,4) + (2.*u2+duy)*wy
      grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) + (2.*v2+dvy)*wy
      grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) + (2.*w2+dwy)*wy
      grp%nodeHelpArray(i2,7) = grp%nodeHelpArray(i2,7) + (2.*u2+duz)*wz
      grp%nodeHelpArray(i2,8) = grp%nodeHelpArray(i2,8) + (2.*v2+dvz)*wz
      grp%nodeHelpArray(i2,9) = grp%nodeHelpArray(i2,9) + (2.*w2+dwz)*wz 
    end do

#ifdef PARALLEL

!    write(get_unit(),*)  'Starting first viscosity comms',kk
    if (kk.eq.1) then
      do i=1,grp%pdp%numberOfProcesses
        if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
          call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
          position = 1
          do j=1,9
            call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                               grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
          end do
          call mp_sendv(grp%pdp%processorIDs,i,29,grp%pdp%sendFlag)
        end if
      end do
    end if
  end do

  call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,29,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1233 continue

  do i=1,grp%pdp%numberOfProcesses
!   if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
!     grp%pdp%sendFlag = 0
!     call mp_check_mesg_arrived(grp%pdp%processorIDs,i,29,grp%pdp%sendFlag,bufferzz,bufferSize)
!     if (grp%pdp%sendFlag.gt.0) then
!       grp%pdp%receivedMessageFromProcess(i) = 1
!       position = 1
        do j=1,9
          call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                                 grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
        end do
!     end if
    end if
  end do

! do i=1,grp%pdp%numberOfProcesses
!   if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1233
! end do

!  write(get_unit(),*)  'Starting second viscosity comms'
  do i=1,grp%pdp%numberOfProcesses
    if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
      call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
      position = 1
      do j=1,9
        call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                          grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
      end do
      call mp_sendv(grp%pdp%processorIDs,i,30,grp%pdp%sendFlag)
    end if
  end do

  call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,30,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1333 continue

  do i=1,grp%pdp%numberOfProcesses
!   if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
!     grp%pdp%sendFlag = 0
!     call mp_check_mesg_arrived(grp%pdp%processorIDs,i,30,grp%pdp%sendFlag,bufferzz,bufferSize)
!     if (grp%pdp%sendFlag.gt.0) then
!       grp%pdp%receivedMessageFromProcess(i) = 1
!       position = 1
        do j = 1,9
          call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                             grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
        end do
!     end if
    end if
  end do

! do i=1,grp%pdp%numberOfProcesses
!   if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1333
! end do
#endif PARALLEL


  twoThirds = 2./3.
  fourThirds = 4./3.
  oneOverReynoldsNumber = 1./ivd%ReynoldsNumber

! divide by volume to remove mass matrix on LHS
IAUX=0
  do i=1,grp%numberOfNodes
    vt = 1./grp%nodeVolume(i)
    
    grp%nodeHelpArray(i,1) = grp%nodeHelpArray(i,1)*vt
    grp%nodeHelpArray(i,2) = grp%nodeHelpArray(i,2)*vt
    grp%nodeHelpArray(i,3) = grp%nodeHelpArray(i,3)*vt
    grp%nodeHelpArray(i,4) = grp%nodeHelpArray(i,4)*vt
    grp%nodeHelpArray(i,5) = grp%nodeHelpArray(i,5)*vt
    grp%nodeHelpArray(i,6) = grp%nodeHelpArray(i,6)*vt
    grp%nodeHelpArray(i,7) = grp%nodeHelpArray(i,7)*vt
    grp%nodeHelpArray(i,8) = grp%nodeHelpArray(i,8)*vt
    grp%nodeHelpArray(i,9) = grp%nodeHelpArray(i,9)*vt


    dudx = grp%nodeHelpArray(i,1)
    dvdx = grp%nodeHelpArray(i,2)
    dwdx = grp%nodeHelpArray(i,3)
    dudy = grp%nodeHelpArray(i,4)
    dvdy = grp%nodeHelpArray(i,5)
    dwdy = grp%nodeHelpArray(i,6)
    dudz = grp%nodeHelpArray(i,7)
    dvdz = grp%nodeHelpArray(i,8)
    dwdz = grp%nodeHelpArray(i,9)

 ! add turbulence effects
      turbViscosityCoefficient = grp%u(i,1)*grp%u(i,6)/grp%u(i,7)

!!    Limiting turbulence
      IF(turbViscosityCoefficient .lt. 0.0)THEN
!!       IF(IAUX==0)PRINT*,'NEGATIVE mut',turbViscosityCoefficient,grp%u(i,1),grp%u(i,6),grp%u(i,7)
       IAUX=1
      END IF 
    turbViscosityCoefficient=max(0.0,turbViscosityCoefficient)
    turbViscosityCoefficient=min(turbViscosityCoefficient,ivd%maxTurbulenceValueMU*oneOverReynoldsNumber)
    mu = turbViscosityCoefficient

! stress tensor
    grp%nodeHelpArray(i,10) = mu*(fourThirds*dudx-twoThirds*(dvdy+dwdz))-twoThirds*grp%u(i,1)*grp%u(i,6) 
    grp%nodeHelpArray(i,11) = mu*(dudy+dvdx) 
    grp%nodeHelpArray(i,12) = mu*(dudz+dwdx) 
    grp%nodeHelpArray(i,13) = mu*(fourThirds*dvdy-twoThirds*(dudx+dwdz))-twoThirds*grp%u(i,1)*grp%u(i,6) 
    grp%nodeHelpArray(i,14) = mu*(dvdz+dwdy) 
    grp%nodeHelpArray(i,15) = mu*(fourThirds*dwdz-twoThirds*(dudx+dvdy))-twoThirds*grp%u(i,1)*grp%u(i,6) 
  end do


IF (IAREA==1)THEN
#ifdef PARALLEL

  do kk=1,2
    if(kk==1) then
      istart = 1
      ifinish = grp%pdp%numberOfComSides
    else
      istart = grp%pdp%numberOfComSides + 1
      ifinish = grp%numberOfSides
    end if
 
 
    do i=istart,ifinish

#else PARALLEL

    do i=1,grp%numberOfSides

#endif PARALLEL

      i1 = grp%sideIndexArray(i,1)
      i2 = grp%sideIndexArray(i,2)

! side weights from preprocessor 
      wx = grp%sideWeightsArray(i,1)
      wy = grp%sideWeightsArray(i,2)  
      wz = grp%sideWeightsArray(i,3) 

      t11 = grp%nodeHelpArray(i1,10)
      t21 = grp%nodeHelpArray(i1,11)
      t31 = grp%nodeHelpArray(i1,12)
      t41 = grp%nodeHelpArray(i1,13)
      t51 = grp%nodeHelpArray(i1,14)
      t61 = grp%nodeHelpArray(i1,15)
      dudx1 = grp%nodeHelpArray(i1,1)
      dvdx1 = grp%nodeHelpArray(i1,2)
      dwdx1 = grp%nodeHelpArray(i1,3)
      dudy1 = grp%nodeHelpArray(i1,4)
      dvdy1 = grp%nodeHelpArray(i1,5)
      dwdy1 = grp%nodeHelpArray(i1,6)
      dudz1 = grp%nodeHelpArray(i1,7)
      dvdz1 = grp%nodeHelpArray(i1,8)
      dwdz1 = grp%nodeHelpArray(i1,9)

      t12 = grp%nodeHelpArray(i2,10)
      t22 = grp%nodeHelpArray(i2,11)
      t32 = grp%nodeHelpArray(i2,12)
      t42 = grp%nodeHelpArray(i2,13)
      t52 = grp%nodeHelpArray(i2,14)
      t62 = grp%nodeHelpArray(i2,15)
      dudx2 = grp%nodeHelpArray(i2,1)
      dvdx2 = grp%nodeHelpArray(i2,2)
      dwdx2 = grp%nodeHelpArray(i2,3)
      dudy2 = grp%nodeHelpArray(i2,4)
      dvdy2 = grp%nodeHelpArray(i2,5)
      dwdy2 = grp%nodeHelpArray(i2,6)
      dudz2 = grp%nodeHelpArray(i2,7)
      dvdz2 = grp%nodeHelpArray(i2,8)
      dwdz2 = grp%nodeHelpArray(i2,9)

      dut1 = dudx1*t11 + dudx2*t12 + dvdx1*t21 + dvdx2*t22 + dwdx1*t31 + dwdx2*t32 
      dut2 = dudy1*t21 + dudy2*t22 + dvdy1*t41 + dvdy2*t42 + dwdy1*t51 + dwdy2*t52 
      dut3 = dudz1*t31 + dudz2*t32 + dvdz1*t51 + dvdz2*t52 + dwdz1*t61 + dwdz2*t62 
      prod = wx*(dut1) + wy*(dut2) + wz*(dut3)
      grp%nodeHelpArray(i1,16) = grp%nodeHelpArray(i1,16) + prod 
      grp%nodeHelpArray(i2,16) = grp%nodeHelpArray(i2,16) - prod

    end do

!Contribution from the boundaries

#ifdef PARALLEL
    if(kk==1) then
      istart = 1
      ifinish = grp%pdp%numberOfBoundaryComSides
    else
      istart = grp%pdp%numberOfBoundaryComSides + 1
      ifinish = grp%brp%numberOfBoundarySides
    end if
    do i=istart,ifinish

#else PARALLEL

    do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL
 
!!     faceIndicator = grp%brp%sideIndexArray(i,3)
!!     if(faceIndicator.ne.2.and.faceIndicator.ne.9) then !symmetry or wall??
      i1 = grp%brp%sideIndexArray(i,1)
      i2 = grp%brp%sideIndexArray(i,2)

      wx   = grp%brp%sideWeightsArray(i,1)
      wy   = grp%brp%sideWeightsArray(i,2) 
      wz   = grp%brp%sideWeightsArray(i,3) 

      t11 = grp%nodeHelpArray(i1,10)
      t21 = grp%nodeHelpArray(i1,11)
      t31 = grp%nodeHelpArray(i1,12)
      t41 = grp%nodeHelpArray(i1,13)
      t51 = grp%nodeHelpArray(i1,14)
      t61 = grp%nodeHelpArray(i1,15)
      dudx1 = grp%nodeHelpArray(i1,1)
      dvdx1 = grp%nodeHelpArray(i1,2)
      dwdx1 = grp%nodeHelpArray(i1,3)
      dudy1 = grp%nodeHelpArray(i1,4)
      dvdy1 = grp%nodeHelpArray(i1,5)
      dwdy1 = grp%nodeHelpArray(i1,6)
      dudz1 = grp%nodeHelpArray(i1,7)
      dvdz1 = grp%nodeHelpArray(i1,8)
      dwdz1 = grp%nodeHelpArray(i1,9)

      t12 = grp%nodeHelpArray(i2,10)
      t22 = grp%nodeHelpArray(i2,11)
      t32 = grp%nodeHelpArray(i2,12)
      t42 = grp%nodeHelpArray(i2,13)
      t52 = grp%nodeHelpArray(i2,14)
      t62 = grp%nodeHelpArray(i2,15)
      dudx2 = grp%nodeHelpArray(i2,1)
      dvdx2 = grp%nodeHelpArray(i2,2)
      dwdx2 = grp%nodeHelpArray(i2,3)
      dudy2 = grp%nodeHelpArray(i2,4)
      dvdy2 = grp%nodeHelpArray(i2,5)
      dwdy2 = grp%nodeHelpArray(i2,6)
      dudz2 = grp%nodeHelpArray(i2,7)
      dvdz2 = grp%nodeHelpArray(i2,8)
      dwdz2 = grp%nodeHelpArray(i2,9)

      dut1 = dudx1*t11 + dudx2*t12 + dvdx1*t21 + dvdx2*t22 + dwdx1*t31 + dwdx2*t32 
      dut2 = dudy1*t21 + dudy2*t22 + dvdy1*t41 + dvdy2*t42 + dwdy1*t51 + dwdy2*t52 
      dut3 = dudz1*t31 + dudz2*t32 + dvdz1*t51 + dvdz2*t52 + dwdz1*t61 + dwdz2*t62 

      prod = wx*(dut1) + wy*(dut2) + wz*(dut3)

      grp%nodeHelpArray(i1,16) = grp%nodeHelpArray(i1,16)+prod+wx*(2*dudx1*t11+dudx2*t12+2*dvdx1*t21+dvdx2*t22+2*dwdx1*t31+dwdx2*t32)+wy*(2*dudy1*t21+dudy2*t22+2*dvdy1*t41+dvdy2*t42+2*dwdy1*t51+dwdy2*t52)+wz*(2*dudz1*t31+dudz2*t32+2*dvdz1*t51+dvdz2*t52+2*dwdz1*t61+dwdz2*t62)
      grp%nodeHelpArray(i2,16) = grp%nodeHelpArray(i2,16)+prod+wx*(dudx1*t11+2*dudx2*t12+dvdx1*t21+2*dvdx2*t22+dwdx1*t31+2*dwdx2*t32)+wy*(dudy1*t21+2*dudy2*t22+dvdy1*t41+2*dvdy2*t42+dwdy1*t51+2*dwdy2*t52)+wz*(dudz1*t31+2*dudz2*t32+dvdz1*t51+2*dvdz2*t52+dwdz1*t61+2*dwdz2*t62) 
!!     end if
   end do

#ifdef PARALLEL

!    write(get_unit(),*)  'Starting first viscosity2 comms',kk
    if (kk.eq.1) then
      do i=1,grp%pdp%numberOfProcesses
        if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
          call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
          position = 1
          
            call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                               grp%pdp%realType,grp%nodeHelpArray(:,16),grp%pdp%buffer,grp%pdp%sendFlag)
          
          call mp_sendv(grp%pdp%processorIDs,i,31,grp%pdp%sendFlag)
        end if
      end do
    end if
  end do

  call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,31,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1234 continue

  do i=1,grp%pdp%numberOfProcesses
!   if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
!     grp%pdp%sendFlag = 0
!     call mp_check_mesg_arrived(grp%pdp%processorIDs,i,31,grp%pdp%sendFlag,bufferzz,bufferSize)
!     if (grp%pdp%sendFlag.gt.0) then
!       grp%pdp%receivedMessageFromProcess(i) = 1
!       position = 1
        
          call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                                 grp%pdp%realType,grp%nodeHelpArray(:,16),grp%pdp%buffer,grp%pdp%sendFlag)
        
!     end if
    end if
  end do

! do i=1,grp%pdp%numberOfProcesses
!   if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1234
! end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

!  write(get_unit(),*)  'Starting second viscosity2 comms'

  do i=1,grp%pdp%numberOfProcesses
    if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
      call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
      position = 1
      
        call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                          grp%pdp%realType,grp%nodeHelpArray(:,16),grp%pdp%buffer,grp%pdp%sendFlag)
      
      call mp_sendv(grp%pdp%processorIDs,i,32,grp%pdp%sendFlag)
    end if
  end do

  call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,32,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1334 continue

  do i=1,grp%pdp%numberOfProcesses
!   if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
!     grp%pdp%sendFlag = 0
!     call mp_check_mesg_arrived(grp%pdp%processorIDs,i,32,grp%pdp%sendFlag,bufferzz,bufferSize)
!     if (grp%pdp%sendFlag.gt.0) then
!       grp%pdp%receivedMessageFromProcess(i) = 1
!       position = 1
        
          call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                             grp%pdp%realType,grp%nodeHelpArray(:,16),grp%pdp%buffer,grp%pdp%sendFlag)
        
!     end if
    end if
  end do

! do i=1,grp%pdp%numberOfProcesses
!   if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1334
! end do
#endif PARALLEL

  do i=1,grp%numberOfNodes 
   vt = grp%nodeVolume(i)     
!   Assemble source term ( = production - dissipation)
! Compressibility effects
! Sarkar's model
   betast=betast0*(1.+2.*grp%u(i,6)/ivd%MachNumber**2)
   beta=beta0*(1.-betast0/beta0*2.*grp%u(i,6)/ivd%MachNumber**2)
! Wilcox's model
! Mto=0.25
! H=0.0
! if(sqrt(2.*abs(grp%u(i,6))/ivd%MachNumber**2) .gt. Mto)H=1.0
! betast=betast0*(1.+1.5*(2.*grp%u(i,6)/ivd%MachNumber**2-Mto**2)*H)
! beta=beta0*(1.-betast0/beta0*1.5*(2.*grp%u(i,6)/ivd%MachNumber**2-Mto**2)*H)

!   for k  
    prodk=grp%nodeHelpArray(i,16)
    dissk=vt*betast*grp%u(i,1)*grp%u(i,6)*grp%u(i,7)
    grp%rhs(i,6) = grp%rhs(i,6)-prodk+dissk

!   for omega
    if(grp%u(i,6) .ne. 0.0)then !at the wall, k is zero but omega is prescribed
    prodw=alfa*grp%u(i,7)/grp%u(i,6)*grp%nodeHelpArray(i,16)
    dissw=vt*beta*grp%u(i,1)*grp%u(i,7)*grp%u(i,7)
    grp%rhs(i,7) = grp%rhs(i,7)-prodw+dissw
     
    end if 

  end do !end for numberOfNodes

ELSE !assuming that the product of vel. derivatives and stress tensor is constant in the volume
 do i=1,grp%numberOfNodes
   vt = grp%nodeVolume(i)
      t11 = grp%nodeHelpArray(i,10)
      t21 = grp%nodeHelpArray(i,11)
      t31 = grp%nodeHelpArray(i,12)
      t41 = grp%nodeHelpArray(i,13)
      t51 = grp%nodeHelpArray(i,14)
      t61 = grp%nodeHelpArray(i,15)
      dudx1 = grp%nodeHelpArray(i,1)
      dvdx1 = grp%nodeHelpArray(i,2)
      dwdx1 = grp%nodeHelpArray(i,3)
      dudy1 = grp%nodeHelpArray(i,4)
      dvdy1 = grp%nodeHelpArray(i,5)
      dwdy1 = grp%nodeHelpArray(i,6)
      dudz1 = grp%nodeHelpArray(i,7)
      dvdz1 = grp%nodeHelpArray(i,8)
      dwdz1 = grp%nodeHelpArray(i,9)

  
!   Multiplying the stress tensor by the velocity gradient
   grp%nodeHelpArray(i,16)=(t11*dudx1+t21*dvdx1+t31*dwdx1)+(t41*dvdy1+t21*dudy1+t51*dwdy1)+(t61*dwdz1+t31*dudz1+t51*dvdz1)
      
!   Assemble source term ( = production - dissipation)  
! Compressibility effects
! Sarkar's model
   betast=betast0*(1.+2.*grp%u(i,6)/ivd%MachNumber**2)
   beta=beta0*(1.-betast0/beta0*2.*grp%u(i,6)/ivd%MachNumber**2)
! Wilcox's model
! Mto=0.25
! H=0.0
! if(sqrt(2.*abs(grp%u(i,6))/ivd%MachNumber**2) .gt. Mto)H=1.0
! betast=betast0*(1.+1.5*(2.*grp%u(i,6)/ivd%MachNumber**2-Mto**2)*H)
! beta=beta0*(1.-betast0/beta0*1.5*(2.*grp%u(i,6)/ivd%MachNumber**2-Mto**2)*H)


!   for k
    prodk=vt*grp%nodeHelpArray(i,16)
    dissk=vt*betast*grp%u(i,1)*grp%u(i,6)*grp%u(i,7)
    grp%rhs(i,6) = grp%rhs(i,6)-prodk+dissk

!   for omega
    if(grp%u(i,6) .ne. 0.0)then !at the wall, k is zero but omega is prescribed
     prodw=vt*alfa*grp%u(i,7)/grp%u(i,6)*grp%nodeHelpArray(i,16)
     dissw=vt*beta*grp%u(i,1)*grp%u(i,7)*grp%u(i,7)
     grp%rhs(i,7) = grp%rhs(i,7)-prodw+dissw   
    end if 

!Source correction for local time stepping
    if(grp%u(i,6) .ne. 0.0)grp%SK(i) = abs(grp%divergence(i))+grp%vorticity(i)+abs(prodk-dissk)/(grp%u(i,1)*grp%u(i,6))
    if(grp%u(i,7) .ne. 0.0)grp%SW(i) = abs(grp%divergence(i))+grp%vorticity(i)+abs(prodw-dissw)/(grp%u(i,1)*grp%u(i,7))

 end do !end for numberOfNodes 

END IF !for IAREA

 end subroutine makeKWSourceTerm
!END-TURB
!-----------------------------------------------------------------------
!START-TURB
subroutine makeSSTSourceTerm(grp,ivd)!SST
 ! source term for k in the k-omega turbulence model equation
 ! the viscosity fluxes are computed by using the wide stencil 
IMPLICIT NONE
include 'mpif.h'

type(GridSolverData) :: grp
type(InputVariablesData) :: ivd

integer :: i1,i2,i,j,ip1,ip2,faceIndicator,ib,ip,ist,ien,ind
real :: r1,r2,u1,u2,v1,v2,wx,wy,wz,oneOverReynoldsNumber
real :: dudx,dvdx,dwdx,dudy,dvdy,dwdy,dudz,dvdz,dwdz
real :: dux,dvx,dwx,duy,dvy,dwy,duz,dvz,dwz
real :: dudx1,dvdx1,dwdx1,dudy1,dvdy1,dwdy1,dudz1,dvdz1,dwdz1
real :: dudx2,dvdx2,dwdx2,dudy2,dvdy2,dwdy2,dudz2,dvdz2,dwdz2
real :: w1,w2,dt1,dt2,dt3,dt4,dt5,dt6,dut1,dut2,dut3
real :: mu,vt,t11,t21,t31,t41,t51,t61,t12,t22,t32,t42,t52,t62
real :: twoThirds,fourThirds
real :: prod,Mto,H,betast,beta,beta1,beta2,gama,gama1,gama2
real :: sigmaw1,sigmaw2
real :: dt,prodk,dissk,prodw,dissw,nut,a1,klog
double precision :: turbViscosityCoefficient

 integer :: kk,istart,ifinish,ii,IAUX,IAREA

  grp%nodeHelpArray(1:grp%numberOfNodes,1:16) = 0.0

!Constant definition
betast=0.09
beta1=0.075
beta2=0.0828
sigmaw1=0.5
sigmaw2=0.856
klog=0.41
gama1=beta1/betast-sigmaw1*klog**2/sqrt(betast)
gama2=beta2/betast-sigmaw2*klog**2/sqrt(betast)
a1=0.31
IAREA = 0

#ifdef PARALLEL

  do kk=1,2
    if(kk==1) then
      istart = 1
      ifinish = grp%pdp%numberOfComSides
    else
      istart = grp%pdp%numberOfComSides + 1
      ifinish = grp%numberOfSides
    end if


    do i=istart,ifinish

#else PARALLEL

    do i=1,grp%numberOfSides

#endif PARALLEL


! indexes of nodes in side
      i1 = grp%sideIndexArray(i,1)
      i2 = grp%sideIndexArray(i,2)

! side weights from preprocessor 
      wx = grp%sideWeightsArray(i,1)  
      wy = grp%sideWeightsArray(i,2)  
      wz = grp%sideWeightsArray(i,3)  

      r1   = grp%u(i1,1)                                  ! density
      u1   = grp%u(i1,2)/r1                               ! x - velocity
      v1   = grp%u(i1,3)/r1                               ! y - velocity
      w1   = grp%u(i1,4)/r1                               ! z - velocity  
  
      r2   = grp%u(i2,1)
      u2   = grp%u(i2,2)/r2
      v2   = grp%u(i2,3)/r2
      w2   = grp%u(i2,4)/r2

      dux   = (u1+u2)*wx
      dvx   = (v1+v2)*wx
      dwx   = (w1+w2)*wx
      duy   = (u1+u2)*wy
      dvy   = (v1+v2)*wy
      dwy   = (w1+w2)*wy
      duz   = (u1+u2)*wz
      dvz   = (v1+v2)*wz
      dwz   = (w1+w2)*wz


      grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + dux
      grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + dvx
      grp%nodeHelpArray(i1,3) = grp%nodeHelpArray(i1,3) + dwx
      grp%nodeHelpArray(i1,4) = grp%nodeHelpArray(i1,4) + duy
      grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) + dvy
      grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) + dwy
      grp%nodeHelpArray(i1,7) = grp%nodeHelpArray(i1,7) + duz
      grp%nodeHelpArray(i1,8) = grp%nodeHelpArray(i1,8) + dvz
      grp%nodeHelpArray(i1,9) = grp%nodeHelpArray(i1,9) + dwz

      grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) - dux
      grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) - dvx
      grp%nodeHelpArray(i2,3) = grp%nodeHelpArray(i2,3) - dwx
      grp%nodeHelpArray(i2,4) = grp%nodeHelpArray(i2,4) - duy
      grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) - dvy
      grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) - dwy
      grp%nodeHelpArray(i2,7) = grp%nodeHelpArray(i2,7) - duz
      grp%nodeHelpArray(i2,8) = grp%nodeHelpArray(i2,8) - dvz
      grp%nodeHelpArray(i2,9) = grp%nodeHelpArray(i2,9) - dwz
    end do

#ifdef PARALLEL
    if(kk==1) then
      istart = 1
      ifinish = grp%pdp%numberOfBoundaryComSides
    else
      istart = grp%pdp%numberOfBoundaryComSides + 1
      ifinish = grp%brp%numberOfBoundarySides
    end if
    do i=istart,ifinish

#else PARALLEL

    do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL

      i1 = grp%brp%sideIndexArray(i,1)
      i2 = grp%brp%sideIndexArray(i,2)
      wx   = grp%brp%sideWeightsArray(i,1)
      wy   = grp%brp%sideWeightsArray(i,2) 
      wz   = grp%brp%sideWeightsArray(i,3) 

      r1   = grp%u(i1,1)
      u1   = grp%u(i1,2)/r1
      v1   = grp%u(i1,3)/r1
      w1   = grp%u(i1,4)/r1
      r2   = grp%u(i2,1)
      u2   = grp%u(i2,2)/r2
      v2   = grp%u(i2,3)/r2
      w2   = grp%u(i2,4)/r2
      dux   = u1+u2
      dvx   = v1+v2
      dwx   = w1+w2
      duy   = u1+u2
      dvy   = v1+v2
      dwy   = w1+w2
      duz   = u1+u2
      dvz   = v1+v2
      dwz   = w1+w2

! adding boundary face contributions        

      grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + (2.*u1+dux)*wx
      grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + (2.*v1+dvx)*wx
      grp%nodeHelpArray(i1,3) = grp%nodeHelpArray(i1,3) + (2.*w1+dwx)*wx
      grp%nodeHelpArray(i1,4) = grp%nodeHelpArray(i1,4) + (2.*u1+duy)*wy
      grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) + (2.*v1+dvy)*wy
      grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) + (2.*w1+dwy)*wy
      grp%nodeHelpArray(i1,7) = grp%nodeHelpArray(i1,7) + (2.*u1+duz)*wz
      grp%nodeHelpArray(i1,8) = grp%nodeHelpArray(i1,8) + (2.*v1+dvz)*wz
      grp%nodeHelpArray(i1,9) = grp%nodeHelpArray(i1,9) + (2.*w1+dwz)*wz
        
      grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + (2.*u2+dux)*wx
      grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) + (2.*v2+dvx)*wx
      grp%nodeHelpArray(i2,3) = grp%nodeHelpArray(i2,3) + (2.*w2+dwx)*wx
      grp%nodeHelpArray(i2,4) = grp%nodeHelpArray(i2,4) + (2.*u2+duy)*wy
      grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) + (2.*v2+dvy)*wy
      grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) + (2.*w2+dwy)*wy
      grp%nodeHelpArray(i2,7) = grp%nodeHelpArray(i2,7) + (2.*u2+duz)*wz
      grp%nodeHelpArray(i2,8) = grp%nodeHelpArray(i2,8) + (2.*v2+dvz)*wz
      grp%nodeHelpArray(i2,9) = grp%nodeHelpArray(i2,9) + (2.*w2+dwz)*wz 
    end do

#ifdef PARALLEL

!    write(get_unit(),*)  'Starting first viscosity comms',kk
    if (kk.eq.1) then
      do i=1,grp%pdp%numberOfProcesses
        if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
          call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
          position = 1
          do j=1,9
            call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                               grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
          end do
          call mp_sendv(grp%pdp%processorIDs,i,29,grp%pdp%sendFlag)
        end if
      end do
    end if
  end do

  call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,29,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1233 continue

  do i=1,grp%pdp%numberOfProcesses
!   if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
!     grp%pdp%sendFlag = 0
!     call mp_check_mesg_arrived(grp%pdp%processorIDs,i,29,grp%pdp%sendFlag,bufferzz,bufferSize)
!     if (grp%pdp%sendFlag.gt.0) then
!       grp%pdp%receivedMessageFromProcess(i) = 1
!       position = 1
        do j=1,9
          call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                                 grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
        end do
!     end if
    end if
  end do

! do i=1,grp%pdp%numberOfProcesses
!   if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1233
! end do

!  write(get_unit(),*)  'Starting second viscosity comms'
  do i=1,grp%pdp%numberOfProcesses
    if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
      call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
      position = 1
      do j=1,9
        call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                          grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
      end do
      call mp_sendv(grp%pdp%processorIDs,i,30,grp%pdp%sendFlag)
    end if
  end do

  call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,30,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1333 continue

  do i=1,grp%pdp%numberOfProcesses
!   if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
!     grp%pdp%sendFlag = 0
!     call mp_check_mesg_arrived(grp%pdp%processorIDs,i,30,grp%pdp%sendFlag,bufferzz,bufferSize)
!     if (grp%pdp%sendFlag.gt.0) then
!       grp%pdp%receivedMessageFromProcess(i) = 1
!       position = 1
        do j = 1,9
          call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                             grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
        end do
!     end if
    end if
  end do

! do i=1,grp%pdp%numberOfProcesses
!   if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1333
! end do
#endif PARALLEL


  twoThirds = 2./3.
  fourThirds = 4./3.
  oneOverReynoldsNumber = 1./ivd%ReynoldsNumber

! divide by volume to remove mass matrix on LHS
IAUX=0
  do i=1,grp%numberOfNodes
    vt = 1./grp%nodeVolume(i)
    
    grp%nodeHelpArray(i,1) = grp%nodeHelpArray(i,1)*vt
    grp%nodeHelpArray(i,2) = grp%nodeHelpArray(i,2)*vt
    grp%nodeHelpArray(i,3) = grp%nodeHelpArray(i,3)*vt
    grp%nodeHelpArray(i,4) = grp%nodeHelpArray(i,4)*vt
    grp%nodeHelpArray(i,5) = grp%nodeHelpArray(i,5)*vt
    grp%nodeHelpArray(i,6) = grp%nodeHelpArray(i,6)*vt
    grp%nodeHelpArray(i,7) = grp%nodeHelpArray(i,7)*vt
    grp%nodeHelpArray(i,8) = grp%nodeHelpArray(i,8)*vt
    grp%nodeHelpArray(i,9) = grp%nodeHelpArray(i,9)*vt


    dudx = grp%nodeHelpArray(i,1)
    dvdx = grp%nodeHelpArray(i,2)
    dwdx = grp%nodeHelpArray(i,3)
    dudy = grp%nodeHelpArray(i,4)
    dvdy = grp%nodeHelpArray(i,5)
    dwdy = grp%nodeHelpArray(i,6)
    dudz = grp%nodeHelpArray(i,7)
    dvdz = grp%nodeHelpArray(i,8)
    dwdz = grp%nodeHelpArray(i,9)

 ! add turbulence effects
!!   nut=a1*grp%u(i,6)/(max(a1*grp%u(i,7),grp%vorticity(i)*grp%F2(i))) !Eddy viscosity limiter
   nut=a1*grp%u(i,6)/(max(a1*grp%u(i,7),grp%imsr(i)*grp%F2(i)))
   turbViscosityCoefficient = grp%u(i,1)*nut

!!    Limiting turbulence
      IF(turbViscosityCoefficient .lt. 0.0)THEN
!!       IF(IAUX==0)PRINT*,'NEGATIVE mut',turbViscosityCoefficient,grp%u(i,1),grp%u(i,6),grp%u(i,7)
       IAUX=1
      END IF 
    turbViscosityCoefficient=max(0.0,turbViscosityCoefficient)
    turbViscosityCoefficient=min(turbViscosityCoefficient,ivd%maxTurbulenceValueMU*oneOverReynoldsNumber)
    mu = turbViscosityCoefficient

! stress tensor
    grp%nodeHelpArray(i,10) = mu*(fourThirds*dudx-twoThirds*(dvdy+dwdz))-twoThirds*grp%u(i,1)*grp%u(i,6) 
    grp%nodeHelpArray(i,11) = mu*(dudy+dvdx) 
    grp%nodeHelpArray(i,12) = mu*(dudz+dwdx) 
    grp%nodeHelpArray(i,13) = mu*(fourThirds*dvdy-twoThirds*(dudx+dwdz))-twoThirds*grp%u(i,1)*grp%u(i,6) 
    grp%nodeHelpArray(i,14) = mu*(dvdz+dwdy) 
    grp%nodeHelpArray(i,15) = mu*(fourThirds*dwdz-twoThirds*(dudx+dvdy))-twoThirds*grp%u(i,1)*grp%u(i,6) 
  end do


IF (IAREA==1)THEN
#ifdef PARALLEL

  do kk=1,2
    if(kk==1) then
      istart = 1
      ifinish = grp%pdp%numberOfComSides
    else
      istart = grp%pdp%numberOfComSides + 1
      ifinish = grp%numberOfSides
    end if
 
 
    do i=istart,ifinish

#else PARALLEL

    do i=1,grp%numberOfSides

#endif PARALLEL

      i1 = grp%sideIndexArray(i,1)
      i2 = grp%sideIndexArray(i,2)

! side weights from preprocessor 
      wx = grp%sideWeightsArray(i,1)
      wy = grp%sideWeightsArray(i,2)  
      wz = grp%sideWeightsArray(i,3) 

      t11 = grp%nodeHelpArray(i1,10)
      t21 = grp%nodeHelpArray(i1,11)
      t31 = grp%nodeHelpArray(i1,12)
      t41 = grp%nodeHelpArray(i1,13)
      t51 = grp%nodeHelpArray(i1,14)
      t61 = grp%nodeHelpArray(i1,15)
      dudx1 = grp%nodeHelpArray(i1,1)
      dvdx1 = grp%nodeHelpArray(i1,2)
      dwdx1 = grp%nodeHelpArray(i1,3)
      dudy1 = grp%nodeHelpArray(i1,4)
      dvdy1 = grp%nodeHelpArray(i1,5)
      dwdy1 = grp%nodeHelpArray(i1,6)
      dudz1 = grp%nodeHelpArray(i1,7)
      dvdz1 = grp%nodeHelpArray(i1,8)
      dwdz1 = grp%nodeHelpArray(i1,9)

      t12 = grp%nodeHelpArray(i2,10)
      t22 = grp%nodeHelpArray(i2,11)
      t32 = grp%nodeHelpArray(i2,12)
      t42 = grp%nodeHelpArray(i2,13)
      t52 = grp%nodeHelpArray(i2,14)
      t62 = grp%nodeHelpArray(i2,15)
      dudx2 = grp%nodeHelpArray(i2,1)
      dvdx2 = grp%nodeHelpArray(i2,2)
      dwdx2 = grp%nodeHelpArray(i2,3)
      dudy2 = grp%nodeHelpArray(i2,4)
      dvdy2 = grp%nodeHelpArray(i2,5)
      dwdy2 = grp%nodeHelpArray(i2,6)
      dudz2 = grp%nodeHelpArray(i2,7)
      dvdz2 = grp%nodeHelpArray(i2,8)
      dwdz2 = grp%nodeHelpArray(i2,9)

      dut1 = dudx1*t11 + dudx2*t12 + dvdx1*t21 + dvdx2*t22 + dwdx1*t31 + dwdx2*t32 
      dut2 = dudy1*t21 + dudy2*t22 + dvdy1*t41 + dvdy2*t42 + dwdy1*t51 + dwdy2*t52 
      dut3 = dudz1*t31 + dudz2*t32 + dvdz1*t51 + dvdz2*t52 + dwdz1*t61 + dwdz2*t62 
      prod = wx*(dut1) + wy*(dut2) + wz*(dut3)
      grp%nodeHelpArray(i1,16) = grp%nodeHelpArray(i1,16) + prod 
      grp%nodeHelpArray(i2,16) = grp%nodeHelpArray(i2,16) - prod

    end do

!Contribution from the boundaries

#ifdef PARALLEL
    if(kk==1) then
      istart = 1
      ifinish = grp%pdp%numberOfBoundaryComSides
    else
      istart = grp%pdp%numberOfBoundaryComSides + 1
      ifinish = grp%brp%numberOfBoundarySides
    end if
    do i=istart,ifinish

#else PARALLEL

    do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL
 
!!     faceIndicator = grp%brp%sideIndexArray(i,3)
!!     if(faceIndicator.ne.2.and.faceIndicator.ne.9) then !symmetry or wall??
      i1 = grp%brp%sideIndexArray(i,1)
      i2 = grp%brp%sideIndexArray(i,2)

      wx   = grp%brp%sideWeightsArray(i,1)
      wy   = grp%brp%sideWeightsArray(i,2) 
      wz   = grp%brp%sideWeightsArray(i,3) 

      t11 = grp%nodeHelpArray(i1,10)
      t21 = grp%nodeHelpArray(i1,11)
      t31 = grp%nodeHelpArray(i1,12)
      t41 = grp%nodeHelpArray(i1,13)
      t51 = grp%nodeHelpArray(i1,14)
      t61 = grp%nodeHelpArray(i1,15)
      dudx1 = grp%nodeHelpArray(i1,1)
      dvdx1 = grp%nodeHelpArray(i1,2)
      dwdx1 = grp%nodeHelpArray(i1,3)
      dudy1 = grp%nodeHelpArray(i1,4)
      dvdy1 = grp%nodeHelpArray(i1,5)
      dwdy1 = grp%nodeHelpArray(i1,6)
      dudz1 = grp%nodeHelpArray(i1,7)
      dvdz1 = grp%nodeHelpArray(i1,8)
      dwdz1 = grp%nodeHelpArray(i1,9)

      t12 = grp%nodeHelpArray(i2,10)
      t22 = grp%nodeHelpArray(i2,11)
      t32 = grp%nodeHelpArray(i2,12)
      t42 = grp%nodeHelpArray(i2,13)
      t52 = grp%nodeHelpArray(i2,14)
      t62 = grp%nodeHelpArray(i2,15)
      dudx2 = grp%nodeHelpArray(i2,1)
      dvdx2 = grp%nodeHelpArray(i2,2)
      dwdx2 = grp%nodeHelpArray(i2,3)
      dudy2 = grp%nodeHelpArray(i2,4)
      dvdy2 = grp%nodeHelpArray(i2,5)
      dwdy2 = grp%nodeHelpArray(i2,6)
      dudz2 = grp%nodeHelpArray(i2,7)
      dvdz2 = grp%nodeHelpArray(i2,8)
      dwdz2 = grp%nodeHelpArray(i2,9)

      dut1 = dudx1*t11 + dudx2*t12 + dvdx1*t21 + dvdx2*t22 + dwdx1*t31 + dwdx2*t32 
      dut2 = dudy1*t21 + dudy2*t22 + dvdy1*t41 + dvdy2*t42 + dwdy1*t51 + dwdy2*t52 
      dut3 = dudz1*t31 + dudz2*t32 + dvdz1*t51 + dvdz2*t52 + dwdz1*t61 + dwdz2*t62 

      prod = wx*(dut1) + wy*(dut2) + wz*(dut3)

      grp%nodeHelpArray(i1,16) = grp%nodeHelpArray(i1,16)+prod+wx*(2*dudx1*t11+dudx2*t12+2*dvdx1*t21+dvdx2*t22+2*dwdx1*t31+dwdx2*t32)+wy*(2*dudy1*t21+dudy2*t22+2*dvdy1*t41+dvdy2*t42+2*dwdy1*t51+dwdy2*t52)+wz*(2*dudz1*t31+dudz2*t32+2*dvdz1*t51+dvdz2*t52+2*dwdz1*t61+dwdz2*t62)
      grp%nodeHelpArray(i2,16) = grp%nodeHelpArray(i2,16)+prod+wx*(dudx1*t11+2*dudx2*t12+dvdx1*t21+2*dvdx2*t22+dwdx1*t31+2*dwdx2*t32)+wy*(dudy1*t21+2*dudy2*t22+dvdy1*t41+2*dvdy2*t42+dwdy1*t51+2*dwdy2*t52)+wz*(dudz1*t31+2*dudz2*t32+dvdz1*t51+2*dvdz2*t52+dwdz1*t61+2*dwdz2*t62) 
!!     end if
   end do

#ifdef PARALLEL

!    write(get_unit(),*)  'Starting first viscosity2 comms',kk
    if (kk.eq.1) then
      do i=1,grp%pdp%numberOfProcesses
        if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
          call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
          position = 1
          
            call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                               grp%pdp%realType,grp%nodeHelpArray(:,16),grp%pdp%buffer,grp%pdp%sendFlag)
          
          call mp_sendv(grp%pdp%processorIDs,i,31,grp%pdp%sendFlag)
        end if
      end do
    end if
  end do

  call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,31,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1234 continue

  do i=1,grp%pdp%numberOfProcesses
!   if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
!     grp%pdp%sendFlag = 0
!     call mp_check_mesg_arrived(grp%pdp%processorIDs,i,31,grp%pdp%sendFlag,bufferzz,bufferSize)
!     if (grp%pdp%sendFlag.gt.0) then
!       grp%pdp%receivedMessageFromProcess(i) = 1
!       position = 1
        
          call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                                 grp%pdp%realType,grp%nodeHelpArray(:,16),grp%pdp%buffer,grp%pdp%sendFlag)
        
!     end if
    end if
  end do

! do i=1,grp%pdp%numberOfProcesses
!   if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1234
! end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

!  write(get_unit(),*)  'Starting second viscosity2 comms'

  do i=1,grp%pdp%numberOfProcesses
    if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
      call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
      position = 1
      
        call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                          grp%pdp%realType,grp%nodeHelpArray(:,16),grp%pdp%buffer,grp%pdp%sendFlag)
      
      call mp_sendv(grp%pdp%processorIDs,i,32,grp%pdp%sendFlag)
    end if
  end do

  call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,32,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!1334 continue

  do i=1,grp%pdp%numberOfProcesses
!   if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
!     grp%pdp%sendFlag = 0
!     call mp_check_mesg_arrived(grp%pdp%processorIDs,i,32,grp%pdp%sendFlag,bufferzz,bufferSize)
!     if (grp%pdp%sendFlag.gt.0) then
!       grp%pdp%receivedMessageFromProcess(i) = 1
!       position = 1
        
          call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                             grp%pdp%realType,grp%nodeHelpArray(:,16),grp%pdp%buffer,grp%pdp%sendFlag)
        
!     end if
    end if
  end do

! do i=1,grp%pdp%numberOfProcesses
!   if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 1334
! end do
#endif PARALLEL

  do i=1,grp%numberOfNodes 
   vt = grp%nodeVolume(i) 
!  for k    
    prodk = grp%nodeHelpArray(i,16)
    dissk=vt*betast*grp%u(i,1)*grp%u(i,6)*grp%u(i,7) 
    prodk = min(prodk,10*dissk) !Production limiter
    grp%rhs(i,6) = grp%rhs(i,6)-prodk+dissk

!   for omega
     gama=grp%F1(i)*gama1+(1.-grp%F1(i))*gama2
     nut=a1*grp%u(i,6)/(max(a1*grp%u(i,7),grp%imsr(i)*grp%F2(i)))
     if (nut .ne. 0.0) then
      prodw = gama/nut*prodk
     else
      prodw=0.0
     end if
     beta=grp%F1(i)*beta1+(1.-grp%F1(i))*beta2
     dissw=vt*beta*grp%u(i,1)*grp%u(i,7)*grp%u(i,7)
     grp%rhs(i,7) = grp%rhs(i,7)-prodw+dissw

  end do !end for numberOfNodes

ELSE !assuming that the product of vel. derivatives and stress tensor is constant in the volume

 do i=1,grp%numberOfNodes
   vt = grp%nodeVolume(i)
      t11 = grp%nodeHelpArray(i,10)
      t21 = grp%nodeHelpArray(i,11)
      t31 = grp%nodeHelpArray(i,12)
      t41 = grp%nodeHelpArray(i,13)
      t51 = grp%nodeHelpArray(i,14)
      t61 = grp%nodeHelpArray(i,15)
      dudx1 = grp%nodeHelpArray(i,1)
      dvdx1 = grp%nodeHelpArray(i,2)
      dwdx1 = grp%nodeHelpArray(i,3)
      dudy1 = grp%nodeHelpArray(i,4)
      dvdy1 = grp%nodeHelpArray(i,5)
      dwdy1 = grp%nodeHelpArray(i,6)
      dudz1 = grp%nodeHelpArray(i,7)
      dvdz1 = grp%nodeHelpArray(i,8)
      dwdz1 = grp%nodeHelpArray(i,9)

  
!   Multiplying the stress tensor by the velocity gradient
   grp%nodeHelpArray(i,16)=(t11*dudx1+t21*dvdx1+t31*dwdx1)+(t41*dvdy1+t21*dudy1+t51*dwdy1)+(t61*dwdz1+t31*dudz1+t51*dvdz1)
      
!   for k
    prodk = vt*grp%nodeHelpArray(i,16)
    dissk = vt*betast*grp%u(i,1)*grp%u(i,6)*grp%u(i,7)
    prodk = min(prodk,10*dissk) !Production limiter
    grp%rhs(i,6) = grp%rhs(i,6)-prodk+dissk

!   for omega 
    nut=a1*grp%u(i,6)/(max(a1*grp%u(i,7),grp%imsr(i)*grp%F2(i)))
    if(nut .ne. 0.0) then
     gama=grp%F1(i)*gama1+(1.-grp%F1(i))*gama2
     prodw = gama/nut*prodk
    else
     prodw=0.0
    end if
    beta=grp%F1(i)*beta1+(1.-grp%F1(i))*beta2
    dissw=vt*beta*grp%u(i,1)*grp%u(i,7)*grp%u(i,7)    

   grp%rhs(i,7) = grp%rhs(i,7)-prodw+dissw

!Source correction for local time stepping
    if(grp%u(i,6) .ne. 0.0)grp%SK(i) = abs(grp%divergence(i))+grp%vorticity(i)+abs(prodk-dissk)/(grp%u(i,1)*grp%u(i,6))
    if(grp%u(i,7) .ne. 0.0)grp%SW(i) = abs(grp%divergence(i))+grp%vorticity(i)+abs(prodw-dissw)/(grp%u(i,1)*grp%u(i,7))   

 end do !end for numberOfNodes

END IF !for IAREA

 end subroutine makeSSTSourceTerm
!END-TURB
!-----------------------------------------------------------------------
 subroutine makeSATripTerm(grp,ivd)
! trip term for the Spalart-Allmaras turbulence model
 IMPLICIT NONE

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

 integer :: i,j,ind,wallInd,wallIndNum
 real :: vt,tripVorticity,tripDeltax,deltau,deltausquared,d,dt,gt,ft1,ct1,ct2,ReynoldsNumber

 ReynoldsNumber = ivd%ReynoldsNumber

 ct1 = 1.0
 ct2 = 2.0

! add trip sources

  do i=1,grp%numberOfTripFieldNodes
   ind = grp%tripNodeFieldIndexes(i,1)
   wallIndNum = grp%tripNodeFieldIndexes(i,2)
   wallInd = grp%tripLineIndexes(wallIndNum)
   vt = grp%nodeVolume(ind)
   deltau = sqrt(sum(grp%u(ind,2:4)*grp%u(ind,2:4)))/grp%u(ind,1)
   deltausquared = deltau*deltau
   tripVorticity = grp%vorticity(wallInd)
   tripDeltax = grp%tripWallLength(wallIndNum)

   d = grp%wallDistance(ind)
   dt = grp%tripNodeFieldDistances(i)
   gt = min(0.1,deltau/(tripVorticity*tripDeltax))
   if(deltausquared>1.0e-20) then
    ft1 = ct1*gt*exp(-ct2*((tripVorticity**2)/deltausquared)*(d**2+(gt**2)*(dt**2)))
   else
    if(d>1.0e-20) then
     ft1 = 0.0
    else
     ft1 = grp%u(ind,1)*ct1*gt
    end if
   end if

   grp%rhs(ind,6) = grp%rhs(ind,6) - ivd%tripFactor*ft1*ReynoldsNumber*vt*deltausquared 
  end do

 end subroutine makeSATripTerm
!-----------------------------------------------------------------------
subroutine makeArtificialDissipation1(grp,ivd)
 ! makes artificial dissipation using a linearly preserving scheme
 IMPLICIT NONE
include 'mpif.h'

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

! a bunch of help variables
 integer :: i,j,i1,i2,ist,ien,ip,ib
 real :: r1,r2,ru1,ru2,u1,u2,h1,h2,p1,p2,wex,wey,wez,rx,ux,vx,wx,tx,ry,uy,vy,wy,ty
 real :: rz,uz,vz,wz,tz,w1,w2,dz,dw4,dw2,w,dw,secondOrderDis,fourthOrderDis
 real :: dr,du,dv,dh,dl,vt,epslm,eps1,gam1,rnx,rny,al,vmo1,rhol,rho1
 real :: uxl,vyl,rhor,uxr,vyr,epsl,presl,hl,epsr,presr,hr,di,d1,rnz,wzl,wzr,wi
 real :: ui,vi,hi,rlam1,rlam2,rlam3,s1,s2,al1x,al2x,cc14,cc24,cc12,cc22
 real :: v1,t1,dr2,du2,dv2,de2,dcf,def,dr4,du4,dv4,de4,d2,d4,d41,d42
 real :: v2,t2,de,dx,dy,df,dp,sp,d0,ci,ci2,af,ucp,u,v,c1,c2,weightNorm,dist
 integer :: ind1,ind2,vind

 real :: maxl(3)
 integer :: maxi(3)

 integer :: kk,istart,ifinish

 ! assumes that rhs is initiated

 ! Linearly preserving artificial dissipation by Crumpton et al
 ! Tenth int. conf. on num. meth. for lam. and turb. flow 

 grp%dissipation = 0.0 
 grp%nodeHelpArray = 0.0


#ifdef PARALLEL

do kk=1,2
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfComSides
 else
  istart = grp%pdp%numberOfComSides + 1
  ifinish = grp%numberOfSides
 end if


 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL


! indexes of nodes in side
  i1 = grp%sideIndexArray(i,1)
  i2 = grp%sideIndexArray(i,2)

! side weights from preprocessor 
  wex = grp%sideWeightsArray(i,1)
  wey = grp%sideWeightsArray(i,2)  
  wez = grp%sideWeightsArray(i,3)  


  r1   = grp%u(i1,1)
  u1   = grp%u(i1,2)
  v1   = grp%u(i1,3)
  w1   = grp%u(i1,4)
  t1   = grp%u(i1,5)
  r2   = grp%u(i2,1)
  u2   = grp%u(i2,2)
  v2   = grp%u(i2,3)
  w2   = grp%u(i2,4)
  t2   = grp%u(i2,5)
  rx   = (r1+r2)*wex
  ux   = (u1+u2)*wex
  vx   = (v1+v2)*wex
  wx   = (w1+w2)*wex
  tx   = (t1+t2)*wex
  ry   = (r1+r2)*wey
  uy   = (u1+u2)*wey
  vy   = (v1+v2)*wey
  wy   = (w1+w2)*wey
  ty   = (t1+t2)*wey
  rz   = (r1+r2)*wez
  uz   = (u1+u2)*wez
  vz   = (v1+v2)*wez
  wz   = (w1+w2)*wez
  tz   = (t1+t2)*wez

! create dU_i/dx_j
  grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + rx
  grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + ux
  grp%nodeHelpArray(i1,3) = grp%nodeHelpArray(i1,3) + vx
  grp%nodeHelpArray(i1,4) = grp%nodeHelpArray(i1,4) + wx
  grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) + tx 

  grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) + ry
  grp%nodeHelpArray(i1,7) = grp%nodeHelpArray(i1,7) + uy
  grp%nodeHelpArray(i1,8) = grp%nodeHelpArray(i1,8) + vy
  grp%nodeHelpArray(i1,9) = grp%nodeHelpArray(i1,9) + wy
  grp%nodeHelpArray(i1,10) = grp%nodeHelpArray(i1,10) + ty

  grp%nodeHelpArray(i1,11) = grp%nodeHelpArray(i1,11) + rz
  grp%nodeHelpArray(i1,12) = grp%nodeHelpArray(i1,12) + uz
  grp%nodeHelpArray(i1,13) = grp%nodeHelpArray(i1,13) + vz
  grp%nodeHelpArray(i1,14) = grp%nodeHelpArray(i1,14) + wz
  grp%nodeHelpArray(i1,15) = grp%nodeHelpArray(i1,15) + tz

  grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) - rx
  grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) - ux
  grp%nodeHelpArray(i2,3) = grp%nodeHelpArray(i2,3) - vx
  grp%nodeHelpArray(i2,4) = grp%nodeHelpArray(i2,4) - wx
  grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) - tx 

  grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) - ry
  grp%nodeHelpArray(i2,7) = grp%nodeHelpArray(i2,7) - uy
  grp%nodeHelpArray(i2,8) = grp%nodeHelpArray(i2,8) - vy
  grp%nodeHelpArray(i2,9) = grp%nodeHelpArray(i2,9) - wy
  grp%nodeHelpArray(i2,10) = grp%nodeHelpArray(i2,10) - ty

  grp%nodeHelpArray(i2,11) = grp%nodeHelpArray(i2,11) - rz
  grp%nodeHelpArray(i2,12) = grp%nodeHelpArray(i2,12) - uz
  grp%nodeHelpArray(i2,13) = grp%nodeHelpArray(i2,13) - vz
  grp%nodeHelpArray(i2,14) = grp%nodeHelpArray(i2,14) - wz
  grp%nodeHelpArray(i2,15) = grp%nodeHelpArray(i2,15) - tz

! make L(U) and L_j(x)
      
  dr = grp%u(i1,1) - grp%u(i2,1)
  du = grp%u(i1,2) - grp%u(i2,2)
  dv = grp%u(i1,3) - grp%u(i2,3)
  dw = grp%u(i1,4) - grp%u(i2,4)
  de = grp%u(i1,5) - grp%u(i2,5)

  if(ivd%useDissipationWeighting) then 
   dl = grp%sideLengthArray(i) 
   if(dl.gt.0) then
    df = 1./dl    
   else
    write(*,*) 'ERROR: Grid points coincide'
    stop
   end if
  else
   df = 1.0
  end if
  
! set help variables for L_j(U)
  grp%nodeHelpArray(i1,16) = grp%nodeHelpArray(i1,16) - dr*df
  grp%nodeHelpArray(i1,17) = grp%nodeHelpArray(i1,17) - du*df
  grp%nodeHelpArray(i1,18) = grp%nodeHelpArray(i1,18) - dv*df
  grp%nodeHelpArray(i1,19) = grp%nodeHelpArray(i1,19) - dw*df
  grp%nodeHelpArray(i1,20) = grp%nodeHelpArray(i1,20) - de*df
  grp%nodeHelpArray(i2,16) = grp%nodeHelpArray(i2,16) + dr*df
  grp%nodeHelpArray(i2,17) = grp%nodeHelpArray(i2,17) + du*df
  grp%nodeHelpArray(i2,18) = grp%nodeHelpArray(i2,18) + dv*df
  grp%nodeHelpArray(i2,19) = grp%nodeHelpArray(i2,19) + dw*df
  grp%nodeHelpArray(i2,20) = grp%nodeHelpArray(i2,20) + de*df

  grp%nodeHelpArray(i1,21) = grp%nodeHelpArray(i1,21) + df
  grp%nodeHelpArray(i2,21) = grp%nodeHelpArray(i2,21) + df 


  if(grp%gridNumber==1) then 
   dx = grp%coordinates(i1,1) - grp%coordinates(i2,1)
   dy = grp%coordinates(i1,2) - grp%coordinates(i2,2)
   dz = grp%coordinates(i1,3) - grp%coordinates(i2,3)
  else
   ind1 = grp%seedNodeRegister(i1)
   ind2 = grp%seedNodeRegister(i2)

   dx = grp%coordinates(ind1,1) - grp%coordinates(ind2,1)
   dy = grp%coordinates(ind1,2) - grp%coordinates(ind2,2)
   dz = grp%coordinates(ind1,3) - grp%coordinates(ind2,3)
  end if

! for linearly preserving factor, L_j(x)
  grp%nodeHelpArray(i1,22) = grp%nodeHelpArray(i1,22) - dx*df
  grp%nodeHelpArray(i1,23) = grp%nodeHelpArray(i1,23) - dy*df
  grp%nodeHelpArray(i1,24) = grp%nodeHelpArray(i1,24) - dz*df
  grp%nodeHelpArray(i2,22) = grp%nodeHelpArray(i2,22) + dx*df
  grp%nodeHelpArray(i2,23) = grp%nodeHelpArray(i2,23) + dy*df     
  grp%nodeHelpArray(i2,24) = grp%nodeHelpArray(i2,24) + dz*df     
 end do

#ifdef PARALLEL
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfBoundaryComSides
 else
  istart = grp%pdp%numberOfBoundaryComSides + 1
  ifinish = grp%brp%numberOfBoundarySides
 end if
 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL

 i1 = grp%brp%sideIndexArray(i,1)
 i2 = grp%brp%sideIndexArray(i,2)
 wex   = grp%brp%sideWeightsArray(i,1)
 wey   = grp%brp%sideWeightsArray(i,2) 
 wez   = grp%brp%sideWeightsArray(i,3) 
 r1   = grp%u(i1,1)
 u1   = grp%u(i1,2)
 v1   = grp%u(i1,3)
 w1   = grp%u(i1,4)
 t1   = grp%u(i1,5)
 r2   = grp%u(i2,1)
 u2   = grp%u(i2,2)
 v2   = grp%u(i2,3)
 w2   = grp%u(i2,4)
 t2   = grp%u(i2,5)

! adding boundary face contributions        

 if(ivd%boundaryTerm==1) then 
  grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + (3.*r1+r2)*wex
  grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + (3.*u1+u2)*wex
  grp%nodeHelpArray(i1,3) = grp%nodeHelpArray(i1,3) + (3.*v1+v2)*wex
  grp%nodeHelpArray(i1,4) = grp%nodeHelpArray(i1,4) + (3.*w1+w2)*wex
  grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) + (3.*t1+t2)*wex
  grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) + (3.*r1+r2)*wey
  grp%nodeHelpArray(i1,7) = grp%nodeHelpArray(i1,7) + (3.*u1+u2)*wey
  grp%nodeHelpArray(i1,8) = grp%nodeHelpArray(i1,8) + (3.*v1+v2)*wey
  grp%nodeHelpArray(i1,9) = grp%nodeHelpArray(i1,9) + (3.*w1+w2)*wey
  grp%nodeHelpArray(i1,10) = grp%nodeHelpArray(i1,10) + (3.*t1+t2)*wey
  grp%nodeHelpArray(i1,11) = grp%nodeHelpArray(i1,11) + (3.*r1+r2)*wez
  grp%nodeHelpArray(i1,12) = grp%nodeHelpArray(i1,12) + (3.*u1+u2)*wez
  grp%nodeHelpArray(i1,13) = grp%nodeHelpArray(i1,13) + (3.*v1+v2)*wez
  grp%nodeHelpArray(i1,14) = grp%nodeHelpArray(i1,14) + (3.*w1+w2)*wez
  grp%nodeHelpArray(i1,15) = grp%nodeHelpArray(i1,15) + (3.*t1+t2)*wez
        
  grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + (r1+3.*r2)*wex
  grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) + (u1+3.*u2)*wex
  grp%nodeHelpArray(i2,3) = grp%nodeHelpArray(i2,3) + (v1+3.*v2)*wex
  grp%nodeHelpArray(i2,4) = grp%nodeHelpArray(i2,4) + (w1+3.*w2)*wex
  grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) + (t1+3.*t2)*wex
  grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) + (r1+3.*r2)*wey
  grp%nodeHelpArray(i2,7) = grp%nodeHelpArray(i2,7) + (u1+3.*u2)*wey
  grp%nodeHelpArray(i2,8) = grp%nodeHelpArray(i2,8) + (v1+3.*v2)*wey
  grp%nodeHelpArray(i2,9) = grp%nodeHelpArray(i2,9) + (w1+3.*w2)*wey
  grp%nodeHelpArray(i2,10) = grp%nodeHelpArray(i2,10) + (t1+3.*t2)*wey
  grp%nodeHelpArray(i2,11) = grp%nodeHelpArray(i2,11) + (r1+3.*r2)*wez
  grp%nodeHelpArray(i2,12) = grp%nodeHelpArray(i2,12) + (u1+3.*u2)*wez
  grp%nodeHelpArray(i2,13) = grp%nodeHelpArray(i2,13) + (v1+3.*v2)*wez
  grp%nodeHelpArray(i2,14) = grp%nodeHelpArray(i2,14) + (w1+3.*w2)*wez
  grp%nodeHelpArray(i2,15) = grp%nodeHelpArray(i2,15) + (t1+3.*t2)*wez
 else
  grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + 4.*r1*wex
  grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + 4.*u1*wex
  grp%nodeHelpArray(i1,3) = grp%nodeHelpArray(i1,3) + 4.*v1*wex
  grp%nodeHelpArray(i1,4) = grp%nodeHelpArray(i1,4) + 4.*w1*wex
  grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) + 4.*t1*wex
  grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) + 4.*r1*wey
  grp%nodeHelpArray(i1,7) = grp%nodeHelpArray(i1,7) + 4.*u1*wey
  grp%nodeHelpArray(i1,8) = grp%nodeHelpArray(i1,8) + 4.*v1*wey
  grp%nodeHelpArray(i1,9) = grp%nodeHelpArray(i1,9) + 4.*w1*wey
  grp%nodeHelpArray(i1,10) = grp%nodeHelpArray(i1,10) + 4.*t1*wey
  grp%nodeHelpArray(i1,11) = grp%nodeHelpArray(i1,11) + 4.*r1*wez
  grp%nodeHelpArray(i1,12) = grp%nodeHelpArray(i1,12) + 4.*u1*wez
  grp%nodeHelpArray(i1,13) = grp%nodeHelpArray(i1,13) + 4.*v1*wez
  grp%nodeHelpArray(i1,14) = grp%nodeHelpArray(i1,14) + 4.*w1*wez
  grp%nodeHelpArray(i1,15) = grp%nodeHelpArray(i1,15) + 4.*t1*wez
  
  grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + 4.*wex
  grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) + 4.*wex
  grp%nodeHelpArray(i2,3) = grp%nodeHelpArray(i2,3) + 4.*wex
  grp%nodeHelpArray(i2,4) = grp%nodeHelpArray(i2,4) + 4.*wex
  grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) + 4.*wex
  grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) + 4.*wey
  grp%nodeHelpArray(i2,7) = grp%nodeHelpArray(i2,7) + 4.*wey
  grp%nodeHelpArray(i2,8) = grp%nodeHelpArray(i2,8) + 4.*wey
  grp%nodeHelpArray(i2,9) = grp%nodeHelpArray(i2,9) + 4.*wey
  grp%nodeHelpArray(i2,10) = grp%nodeHelpArray(i2,10) + 4.*wey
  grp%nodeHelpArray(i2,11) = grp%nodeHelpArray(i2,11) + 4.*wez
  grp%nodeHelpArray(i2,12) = grp%nodeHelpArray(i2,12) + 4.*wez
  grp%nodeHelpArray(i2,13) = grp%nodeHelpArray(i2,13) + 4.*wez
  grp%nodeHelpArray(i2,14) = grp%nodeHelpArray(i2,14) + 4.*wez
  grp%nodeHelpArray(i2,15) = grp%nodeHelpArray(i2,15) + 4.*wez
 end if
end do
#ifdef PARALLEL

 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses
   if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
    do j=1,24
     call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
    call mp_sendv(grp%pdp%processorIDs,i,41,grp%pdp%sendFlag)
   end if
  end do
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,41,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
    do j=1,24
     call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                            grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
  end if
 end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

!  write(get_unit(),*)  'Starting second ArtDis comms'
 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position = 1
   do j=1,24
    call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                      grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
   end do
   call mp_sendv(grp%pdp%processorIDs,i,42,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,42,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
    do j = 1,24
     call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
  end if
 end do

#endif PARALLEL

! scale if dissipation weighting is used

 if(ivd%useDissipationWeighting) then
  do i=1,grp%numberOfNodes
   grp%nodeHelpArray(i,21) = 3.0*grp%nodeHelpArray(i,21)
  end do
 end if

! divide by volume to remove mass matrix on LHS
do i=1,grp%numberOfNodes
 vt = 1./grp%nodeVolume(i)
 grp%nodeHelpArray(i,1) = grp%nodeHelpArray(i,1)*vt
 grp%nodeHelpArray(i,2) = grp%nodeHelpArray(i,2)*vt
 grp%nodeHelpArray(i,3) = grp%nodeHelpArray(i,3)*vt
 grp%nodeHelpArray(i,4) = grp%nodeHelpArray(i,4)*vt
 grp%nodeHelpArray(i,5) = grp%nodeHelpArray(i,5)*vt
 grp%nodeHelpArray(i,6) = grp%nodeHelpArray(i,6)*vt
 grp%nodeHelpArray(i,7) = grp%nodeHelpArray(i,7)*vt
 grp%nodeHelpArray(i,8) = grp%nodeHelpArray(i,8)*vt
 grp%nodeHelpArray(i,9) = grp%nodeHelpArray(i,9)*vt
 grp%nodeHelpArray(i,10) = grp%nodeHelpArray(i,10)*vt
 grp%nodeHelpArray(i,11) = grp%nodeHelpArray(i,11)*vt
 grp%nodeHelpArray(i,12) = grp%nodeHelpArray(i,12)*vt
 grp%nodeHelpArray(i,13) = grp%nodeHelpArray(i,13)*vt
 grp%nodeHelpArray(i,14) = grp%nodeHelpArray(i,14)*vt
 grp%nodeHelpArray(i,15) = grp%nodeHelpArray(i,15)*vt
end do

! make L^hat*grp%nodeHelpArray(21,x) 
!               = (L(U)- dU_i/dx_j L_j(x))*grp%nodeHelpArray(21,x)

 grp%nodeHelpArray(1:grp%brp%numberOfBoundaryNodes,22:24) = 0.0
 do i=1,grp%numberOfNodes                        
  grp%nodeHelpArray(i,16) = grp%nodeHelpArray(i,16)&
    +grp%nodeHelpArray(i,1)*grp%nodeHelpArray(i,22)&
    +grp%nodeHelpArray(i,6)*grp%nodeHelpArray(i,23)&       
    +grp%nodeHelpArray(i,11)*grp%nodeHelpArray(i,24)       
  grp%nodeHelpArray(i,17) = grp%nodeHelpArray(i,17)&
    +grp%nodeHelpArray(i,2)*grp%nodeHelpArray(i,22)&
    +grp%nodeHelpArray(i,7)*grp%nodeHelpArray(i,23)& 
    +grp%nodeHelpArray(i,12)*grp%nodeHelpArray(i,24) 
  grp%nodeHelpArray(i,18) = grp%nodeHelpArray(i,18)& 
    +grp%nodeHelpArray(i,3)*grp%nodeHelpArray(i,22)&
    +grp%nodeHelpArray(i,8)*grp%nodeHelpArray(i,23)&
    +grp%nodeHelpArray(i,13)*grp%nodeHelpArray(i,24)
  grp%nodeHelpArray(i,19) = grp%nodeHelpArray(i,19)&  
    +grp%nodeHelpArray(i,4)*grp%nodeHelpArray(i,22)&
    +grp%nodeHelpArray(i,9)*grp%nodeHelpArray(i,23)&           
    +grp%nodeHelpArray(i,14)*grp%nodeHelpArray(i,24)           
  grp%nodeHelpArray(i,20) = grp%nodeHelpArray(i,20)&  
    +grp%nodeHelpArray(i,5)*grp%nodeHelpArray(i,22)&
    +grp%nodeHelpArray(i,10)*grp%nodeHelpArray(i,23)&           
    +grp%nodeHelpArray(i,15)*grp%nodeHelpArray(i,24)           
 end do

! remove fourth order viscosity at farfield boundary

ist = grp%brp%nodeIndicatorRegister(-16)
ien = grp%brp%nodeIndicatorRegister(-15)

do ib =ist,ien
 ip = grp%brp%nodeIndicatorArray(ib)
 grp%nodeHelpArray(ip,9:12) = 0.0
end do

! make pressure switch

! no need for gradient matrix anymore
do i = 1,grp%numberOfNodes
 grp%nodeHelpArray(i,1) = 0.0
 grp%nodeHelpArray(i,2) = 0.0
end do

#ifdef PARALLEL

do kk=1,2
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfComSides
 else
  istart = grp%pdp%numberOfComSides + 1
  ifinish = grp%numberOfSides
 end if


 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL


 i1 = grp%sideIndexArray(i,1)
 i2 = grp%sideIndexArray(i,2)
 dp = grp%p(i1) - grp%p(i2)
 sp = grp%p(i1) + grp%p(i2)
! set help variables
 grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) - dp
 grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) + sp
 grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + dp
 grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) + sp
end do
#ifdef PARALLEL

 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses
   if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
    do j=1,2
     call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
    call mp_sendv(grp%pdp%processorIDs,i,43,grp%pdp%sendFlag)
   end if
  end do
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,43,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
    do j=1,2
     call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                            grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
  end if
 end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position = 1
   do j=1,2
    call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                      grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
   end do
   call mp_sendv(grp%pdp%processorIDs,i,44,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,44,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
    do j = 1,2
     call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
  end if
 end do

#endif PARALLEL

! pressure switch, see Hirsh vol2 p.280
do i=1,grp%numberOfNodes
 grp%nodeHelpArray(i,9) = 12.*abs(grp%nodeHelpArray(i,1)/grp%nodeHelpArray(i,2))
end do

! assemble smoother

secondOrderDis = ivd%secondOrderDissipationFactor
if(grp%gridNumber==1) then
 fourthOrderDis = ivd%fourthOrderDissipationFactor
else
 fourthOrderDis = ivd%coarseGridDissipationFactor
end if


#ifdef PARALLEL

do kk=1,2
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfComSides
 else
  istart = grp%pdp%numberOfComSides + 1
  ifinish = grp%numberOfSides
 end if


 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL

 i1 = grp%sideIndexArray(i,1)
 i2 = grp%sideIndexArray(i,2)
! making L^hat differences
 dr4 = grp%nodeHelpArray(i1,16) - grp%nodeHelpArray(i2,16)
 du4 = grp%nodeHelpArray(i1,17) - grp%nodeHelpArray(i2,17)
 dv4 = grp%nodeHelpArray(i1,18) - grp%nodeHelpArray(i2,18)
 dw4 = grp%nodeHelpArray(i1,19) - grp%nodeHelpArray(i2,19)
 de4 = grp%nodeHelpArray(i1,20) - grp%nodeHelpArray(i2,20) 

! making second derivative differences
 dr2 = grp%u(i1,1) - grp%u(i2,1)
 du2 = grp%u(i1,2) - grp%u(i2,2)
 dv2 = grp%u(i1,3) - grp%u(i2,3)
 dw2 = grp%u(i1,4) - grp%u(i2,4)
 de2 = grp%u(i1,5) - grp%u(i2,5)        

 if(ivd%useMatrixDissipation) then 
 ! scalar dissipation coefficient
  d0  = amin1(grp%nodeVolume(i1)/grp%localTimeSteps(i1),&
              grp%nodeVolume(i2)/grp%localTimeSteps(i2))/&
     (grp%nodeConnectivityArray(i1)+grp%nodeConnectivityArray(i2))
  d2  = d0*secondOrderDis&
           *amax1(grp%nodeHelpArray(i1,9),grp%nodeHelpArray(i2,9))
  d4  = d0*fourthOrderDis
  d4  = dim(d4,d2)
  d41  = d4/grp%nodeHelpArray(i1,21)
  d42  = d4/grp%nodeHelpArray(i2,21)                    
 else
! *** One-dimensional Riemann solver along Sij

  eps1 = 1./(max(ivd%HartensCorrectionFactor,1.e-5))
  gam1 = ivd%gamma-1.

  rnx   = grp%sideWeightsArray(i,1)
  rny   = grp%sideWeightsArray(i,2)
  rnz   = grp%sideWeightsArray(i,3)

  al    = sqrt(rnx*rnx+rny*rny+rnz*rnz)
  vmo1  = 1./al
  rnx   = rnx*vmo1
  rny   = rny*vmo1
  rnz   = rnz*vmo1

  rhol  = grp%u(i1,1)
  rho1  = 1./rhol
  uxl   = grp%u(i1,2)*rho1
  vyl   = grp%u(i1,3)*rho1
  wzl   = grp%u(i1,4)*rho1
  epsl  = grp%u(i1,5)
  presl = grp%p(i1)
  hl    = (epsl + presl)*rho1

  rhor  = grp%u(i2,1)
  rho1  = 1./rhor
  uxr   = grp%u(i2,2)*rho1
  vyr   = grp%u(i2,3)*rho1
  wzr   = grp%u(i2,4)*rho1
  epsr  = grp%u(i2,5)
  presr = grp%p(i2)
  hr    = (epsr + presr)*rho1

  di    = sqrt(rhor/rhol)
  d1    = 1.0/(di+1.0)
  ui    = (di*uxr+uxl)*d1
  vi    = (di*vyr+vyl)*d1
  wi    = (di*wzr+wzl)*d1
  hi    = (di*hr+hl)*d1
  ci2   = gam1*(hi-0.5*(ui*ui+vi*vi+wi*wi))
  ci    = sqrt(ci2)
  af    = 0.5*(ui*ui+vi*vi+wi*wi)
  ucp   = ui*rnx+vi*rny+wi*rnz

! *** eigenvalues

  rlam1 = abs(ucp+ci)
  rlam2 = abs(ucp-ci)
  rlam3 = abs(ucp)

! *** Harten's correction

  if(rlam1.lt.ivd%HartensCorrectionFactor) rlam1 = 0.5*(rlam1*rlam1*eps1+ivd%HartensCorrectionFactor)
  if(rlam2.lt.ivd%HartensCorrectionFactor) rlam2 = 0.5*(rlam2*rlam2*eps1+ivd%HartensCorrectionFactor)
  if(rlam3.lt.ivd%HartensCorrectionFactor) rlam3 = 0.5*(rlam3*rlam3*eps1+ivd%HartensCorrectionFactor)

! *** dissipation terms "a la Turkel"

  s1    = 0.5*(rlam1+rlam2)
  s2    = 0.5*(rlam1-rlam2)

  al1x  = gam1*(af*dr4-ui*du4-vi*dv4-wi*dw4+de4)
  al2x  = -ucp*dr4+du4*rnx+dv4*rny+dw4*rnz
  cc14   = ((s1-rlam3)*al1x/ci2)+(s2*al2x/ci)
  cc24   = (s2*al1x/ci)+(s1-rlam3)*al2x

  al1x  = gam1*(af*dr2-ui*du2-vi*dv2-wi*dw2+de2)
  al2x  = -ucp*dr2+du2*rnx+dv2*rny+dw2*rnz
  cc12   = ((s1-rlam3)*al1x/ci2)+(s2*al2x/ci)
  cc22   = (s2*al1x/ci)+(s1-rlam3)*al2x

! *** [A]*(Uj-Ui)

  dcf  =  al
  dr4   = dcf*(rlam3*dr4+cc14           )
  du4   = dcf*(rlam3*du4+cc14*ui+cc24*rnx)
  dv4   = dcf*(rlam3*dv4+cc14*vi+cc24*rny)
  dw4   = dcf*(rlam3*dw4+cc14*wi+cc24*rnz)
  de4   = dcf*(rlam3*de4+cc14*hi+cc24*ucp)

  dr2   = dcf*(rlam3*dr2+cc12           )
  du2   = dcf*(rlam3*du2+cc12*ui+cc22*rnx)
  dv2   = dcf*(rlam3*dv2+cc12*vi+cc22*rny)
  dw2   = dcf*(rlam3*dw2+cc12*wi+cc22*rnz)
  de2   = dcf*(rlam3*de2+cc12*hi+cc22*ucp)

! set dissipation coefficients

  d2  = secondOrderDis&
           *amax1(grp%nodeHelpArray(i1,9),grp%nodeHelpArray(i2,9))
  d4  = fourthOrderDis
  d4  = dim(d4,d2)
  d41  = d4
  d42  = d4
 end if

! update dissipation matrix 
 grp%dissipation(i1,1) = grp%dissipation(i1,1) - dr4*d41+dr2*d2
 grp%dissipation(i1,2) = grp%dissipation(i1,2) - du4*d41+du2*d2
 grp%dissipation(i1,3) = grp%dissipation(i1,3) - dv4*d41+dv2*d2
 grp%dissipation(i1,4) = grp%dissipation(i1,4) - dw4*d41+dw2*d2
 grp%dissipation(i1,5) = grp%dissipation(i1,5) - de4*d41+de2*d2
 grp%dissipation(i2,1) = grp%dissipation(i2,1) + dr4*d42-dr2*d2
 grp%dissipation(i2,2) = grp%dissipation(i2,2) + du4*d42-du2*d2
 grp%dissipation(i2,3) = grp%dissipation(i2,3) + dv4*d42-dv2*d2
 grp%dissipation(i2,4) = grp%dissipation(i2,4) + dw4*d42-dw2*d2    
 grp%dissipation(i2,5) = grp%dissipation(i2,5) + de4*d42-de2*d2    
end do

#ifdef PARALLEL

 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses
   if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
    do j=1,5
     call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%dissipation(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
    call mp_sendv(grp%pdp%processorIDs,i,45,grp%pdp%sendFlag)
   end if
  end do
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,45,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
    do j=1,5
     call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                            grp%pdp%realType,grp%dissipation(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
  end if
 end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position = 1
   do j=1,5
    call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                      grp%pdp%realType,grp%dissipation(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
   end do
   call mp_sendv(grp%pdp%processorIDs,i,46,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,46,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
    do j = 1,5
     call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%dissipation(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
  end if
 end do

#endif PARALLEL


end subroutine makeArtificialDissipation1
!--------------------------------------------------------------------------
subroutine makeArtificialDissipation2(grp,ivd)
! makes artificial dissipation of the JST type
 IMPLICIT NONE
include 'mpif.h'

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

! a bunch of help variables
 integer :: i,j,i1,i2,ist,ien,ib,coarseNodeNumber,ip,ind
 real :: r1,r2,ru1,ru2,u1,u2,h1,h2,p1,p2,wx,wy,rx,ux,vx,tx,ry,uy,vy,ty 
 real :: dr,du,dv,dh,dl,vt,epslm,eps1,gam1,rnx,rny,rnz,al,vmo1,rhol,rho1
 real :: uxl,vyl,wzl,rhor,uxr,vyr,wzr,epsl,presl,hl,epsr,presr,hr,di,d1,dw,dz
 real :: ui,vi,wi,hi,rlam1,rlam2,rlam3,s1,s2,al1x,al2x,cc14,cc24,cc12,cc22
 real :: v1,t1,dr2,du2,dv2,de2,dcf,def,dr4,du4,dv4,de4,d2,d4,d41,d42
 real :: v2,t2,de,dx,dy,df,dp,sp,d0,ci,ci2,af,ucp,dw4,dw2,fact,secondOrderDis
 real :: fourthOrderDis

 integer :: kk,istart,ifinish

 ! Jameson d2-d4 dissipation 

 grp%dissipation = 0.0 
 grp%nodeHelpArray = 0.0

 if (ivd%fluxType.eq.4.and.grp%gridNumber.eq.1) return
 
 if(ivd%fourthOrderDissipationFactor.eq.0.and.ivd%SecondOrderDissipationFactor.eq.0) return
 
#ifdef PARALLEL

do kk=1,2
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfComSides
 else
  istart = grp%pdp%numberOfComSides + 1
  ifinish = grp%numberOfSides
 end if


 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL


! indexes of nodes in side
  i1 = grp%sideIndexArray(i,1)
  i2 = grp%sideIndexArray(i,2)

! make L(U) 
      
  dr = grp%u(i1,1) - grp%u(i2,1)
  du = grp%u(i1,2) - grp%u(i2,2)
  dv = grp%u(i1,3) - grp%u(i2,3)
  dw = grp%u(i1,4) - grp%u(i2,4)
  de = grp%u(i1,5) - grp%u(i2,5) + grp%p(i1) - grp%p(i2)

  if(ivd%useDissipationWeighting) then
   dl = grp%sideLengthArray(i)
   if(dl.gt.0) then
    df = 1./dl
   else
    write(*,*) 'ERROR: Grid points coincide'
    stop
   end if
  else
   df = 1.0
  end if

  grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) - dr*df
  grp%nodeHelpArray(i1,2) = grp%nodeHelpArray(i1,2) - du*df
  grp%nodeHelpArray(i1,3) = grp%nodeHelpArray(i1,3) - dv*df
  grp%nodeHelpArray(i1,4) = grp%nodeHelpArray(i1,4) - dw*df
  grp%nodeHelpArray(i1,5) = grp%nodeHelpArray(i1,5) - de*df
  grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + dr*df
  grp%nodeHelpArray(i2,2) = grp%nodeHelpArray(i2,2) + du*df
  grp%nodeHelpArray(i2,3) = grp%nodeHelpArray(i2,3) + dv*df
  grp%nodeHelpArray(i2,4) = grp%nodeHelpArray(i2,4) + dw*df
  grp%nodeHelpArray(i2,5) = grp%nodeHelpArray(i2,5) + de*df

  grp%nodeHelpArray(i1,6) = grp%nodeHelpArray(i1,6) + df
  grp%nodeHelpArray(i2,6) = grp%nodeHelpArray(i2,6) + df 
end do
#ifdef PARALLEL

 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses
   if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
    do j=1,6
     call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
    call mp_sendv(grp%pdp%processorIDs,i,47,grp%pdp%sendFlag)
   end if
  end do
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,47,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
    do j=1,6
     call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                            grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
  end if
 end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position = 1
   do j=1,6
    call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                      grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
   end do
   call mp_sendv(grp%pdp%processorIDs,i,48,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,48,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
    do j = 1,6
     call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
  end if
 end do

#endif PARALLEL


#ifdef PARALLEL

do kk=1,2
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfComSides
 else
  istart = grp%pdp%numberOfComSides + 1
  ifinish = grp%numberOfSides
 end if


 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL


 i1 = grp%sideIndexArray(i,1)
 i2 = grp%sideIndexArray(i,2)
 dp = grp%p(i1) - grp%p(i2)
 sp = grp%p(i1) + grp%p(i2)
! set help variables
 grp%nodeHelpArray(i1,7) = grp%nodeHelpArray(i1,7) - dp
 grp%nodeHelpArray(i1,8) = grp%nodeHelpArray(i1,8) + sp
 grp%nodeHelpArray(i2,7) = grp%nodeHelpArray(i2,7) + dp
 grp%nodeHelpArray(i2,8) = grp%nodeHelpArray(i2,8) + sp
end do
#ifdef PARALLEL

 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses
   if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
    do j=7,8
     call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
    call mp_sendv(grp%pdp%processorIDs,i,49,grp%pdp%sendFlag)
   end if
  end do
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,49,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
    do j=7,8
     call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                            grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
  end if
 end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position = 1
   do j=7,8
    call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                      grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
   end do
   call mp_sendv(grp%pdp%processorIDs,i,50,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,50,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
    do j = 7,8
     call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%nodeHelpArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
  end if
 end do

#endif PARALLEL

! pressure switch, see Hirsch vol2 p.280
do i=1,grp%numberOfNodes
 fact = 1./grp%nodeConnectivityArray(i)
 grp%nodeHelpArray(i,9) = 12.*abs(grp%nodeHelpArray(i,7)/grp%nodeHelpArray(i,8))
 grp%nodeHelpArray(i,1:5) = fact*grp%nodeHelpArray(i,1:5)
end do

! assemble smoother

secondOrderDis = ivd%secondOrderDissipationFactor
if(grp%gridNumber==1) then 
 fourthOrderDis = ivd%fourthOrderDissipationFactor
else
 fourthOrderDis = ivd%coarseGridDissipationFactor 
end if

#ifdef PARALLEL

do kk=1,2
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfComSides
 else
  istart = grp%pdp%numberOfComSides + 1
  ifinish = grp%numberOfSides
 end if


 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL


 i1 = grp%sideIndexArray(i,1)
 i2 = grp%sideIndexArray(i,2)

! making L differences
 dr4 = grp%nodeHelpArray(i1,1) - grp%nodeHelpArray(i2,1)
 du4 = grp%nodeHelpArray(i1,2) - grp%nodeHelpArray(i2,2)
 dv4 = grp%nodeHelpArray(i1,3) - grp%nodeHelpArray(i2,3)
 dw4 = grp%nodeHelpArray(i1,4) - grp%nodeHelpArray(i2,4)
 de4 = grp%nodeHelpArray(i1,5) - grp%nodeHelpArray(i2,5) 

! making second derivative differences
 dr2 = grp%u(i1,1) - grp%u(i2,1)
 du2 = grp%u(i1,2) - grp%u(i2,2)
 dv2 = grp%u(i1,3) - grp%u(i2,3)
 dw2 = grp%u(i1,4) - grp%u(i2,4)
 de2 = grp%u(i1,5) - grp%u(i2,5) !+ grp%p(i1) - grp%p(i2)       

  ! scalar dissipation coefficient
 if(.not.ivd%useMatrixDissipation) then 
  d0  = amin1(grp%nodeVolume(i1)/grp%localTimeSteps(i1),&
              grp%nodeVolume(i2)/grp%localTimeSteps(i2))/&
     (grp%nodeConnectivityArray(i1)+grp%nodeConnectivityArray(i2)) 
  d2  = d0*secondOrderDis&
           *amax1(grp%nodeHelpArray(i1,9),grp%nodeHelpArray(i2,9))
  d4  = d0*fourthOrderDis
  d4  = dim(d4,d2)
  d41  = d4
  d42  = d4                    


 else
! *** One-dimensional Riemann solver along Sij

  eps1 = 1./(max(ivd%HartensCorrectionFactor,1.e-5))
  gam1 = ivd%gamma-1.

  rnx   = grp%sideWeightsArray(i,1)
  rny   = grp%sideWeightsArray(i,2)
  rnz   = grp%sideWeightsArray(i,3)

  al    = sqrt(rnx*rnx+rny*rny+rnz*rnz)
  vmo1  = 1./al
  rnx   = rnx*vmo1
  rny   = rny*vmo1
  rnz   = rnz*vmo1

  rhol  = grp%u(i1,1)
  rho1  = 1./rhol
  uxl   = grp%u(i1,2)*rho1
  vyl   = grp%u(i1,3)*rho1
  wzl   = grp%u(i1,4)*rho1
  epsl  = grp%u(i1,5)
  presl = grp%p(i1)
  hl    = (epsl + presl)*rho1

  rhor  = grp%u(i2,1)
  rho1  = 1./rhor
  uxr   = grp%u(i2,2)*rho1
  vyr   = grp%u(i2,3)*rho1
  wzr   = grp%u(i2,4)*rho1
  epsr  = grp%u(i2,5)
  presr = grp%p(i2)
  hr    = (epsr + presr)*rho1

  di    = sqrt(rhor/rhol)
  d1    = 1.0/(di+1.0)
  ui    = (di*uxr+uxl)*d1
  vi    = (di*vyr+vyl)*d1
  wi    = (di*wzr+wzl)*d1
  hi    = (di*hr+hl)*d1
  ci2   = gam1*(hi-0.5*(ui*ui+vi*vi+wi*wi))
  ci    = sqrt(ci2)
  af    = 0.5*(ui*ui+vi*vi+wi*wi)
  ucp   = ui*rnx+vi*rny+wi*rnz

! *** eigenvalues

  rlam1 = abs(ucp+ci)
  rlam2 = abs(ucp-ci)
  rlam3 = abs(ucp)

! *** Harten's correction

  if(rlam1.lt.ivd%HartensCorrectionFactor) rlam1 = 0.5*(rlam1*rlam1*eps1+ivd%HartensCorrectionFactor)
  if(rlam2.lt.ivd%HartensCorrectionFactor) rlam2 = 0.5*(rlam2*rlam2*eps1+ivd%HartensCorrectionFactor)
  if(rlam3.lt.ivd%HartensCorrectionFactor) rlam3 = 0.5*(rlam3*rlam3*eps1+ivd%HartensCorrectionFactor)

! *** dissipation terms "a la Turkel"

  s1    = 0.5*(rlam1+rlam2)
  s2    = 0.5*(rlam1-rlam2)

  al1x  = gam1*(af*dr4-ui*du4-vi*dv4-wi*dw4+de4)
  al2x  = -ucp*dr4+du4*rnx+dv4*rny+dw4*rnz
  cc14   = ((s1-rlam3)*al1x/ci2)+(s2*al2x/ci)
  cc24   = (s2*al1x/ci)+(s1-rlam3)*al2x

  al1x  = gam1*(af*dr2-ui*du2-vi*dv2-wi*dw2+de2)
  al2x  = -ucp*dr2+du2*rnx+dv2*rny+dw2*rnz
  cc12   = ((s1-rlam3)*al1x/ci2)+(s2*al2x/ci)
  cc22   = (s2*al1x/ci)+(s1-rlam3)*al2x

! *** [A]*(Uj-Ui)

  dcf  =  al
  dr4   = dcf*(rlam3*dr4+cc14           )
  du4   = dcf*(rlam3*du4+cc14*ui+cc24*rnx)
  dv4   = dcf*(rlam3*dv4+cc14*vi+cc24*rny)
  dw4   = dcf*(rlam3*dw4+cc14*wi+cc24*rnz)
  de4   = dcf*(rlam3*de4+cc14*hi+cc24*ucp)

  dr2   = dcf*(rlam3*dr2+cc12           )
  du2   = dcf*(rlam3*du2+cc12*ui+cc22*rnx)
  dv2   = dcf*(rlam3*dv2+cc12*vi+cc22*rny)
  dw2   = dcf*(rlam3*dw2+cc12*wi+cc22*rnz)
  de2   = dcf*(rlam3*de2+cc12*hi+cc22*ucp)

! set dissipation coefficients

  d2  = secondOrderDis&
           *amax1(grp%nodeHelpArray(i1,9),grp%nodeHelpArray(i2,9))
  d4  = fourthOrderDis
  d4  = dim(d4,d2)
  d41  = d4
  d42  = d4
 end if

! update dissipation matrix 
 grp%dissipation(i1,1) = grp%dissipation(i1,1) - dr4*d41+dr2*d2
 grp%dissipation(i1,2) = grp%dissipation(i1,2) - du4*d41+du2*d2
 grp%dissipation(i1,3) = grp%dissipation(i1,3) - dv4*d41+dv2*d2
 grp%dissipation(i1,4) = grp%dissipation(i1,4) - dw4*d41+dw2*d2
 grp%dissipation(i1,5) = grp%dissipation(i1,5) - de4*d41+de2*d2
 grp%dissipation(i2,1) = grp%dissipation(i2,1) + dr4*d42-dr2*d2
 grp%dissipation(i2,2) = grp%dissipation(i2,2) + du4*d42-du2*d2
 grp%dissipation(i2,3) = grp%dissipation(i2,3) + dv4*d42-dv2*d2
 grp%dissipation(i2,4) = grp%dissipation(i2,4) + dw4*d42-dw2*d2
 grp%dissipation(i2,5) = grp%dissipation(i2,5) + de4*d42-de2*d2  
end do
#ifdef PARALLEL

 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses
   if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
    do j=1,5
     call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%dissipation(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
    call mp_sendv(grp%pdp%processorIDs,i,51,grp%pdp%sendFlag)
   end if
  end do
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,51,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
    do j=1,5
     call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                            grp%pdp%realType,grp%dissipation(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
  end if
 end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position = 1
   do j=1,5
    call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                      grp%pdp%realType,grp%dissipation(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
   end do
   call mp_sendv(grp%pdp%processorIDs,i,52,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,52,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
    do j = 1,5
     call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%dissipation(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
  end if
 end do

#endif PARALLEL

end subroutine makeArtificialDissipation2
!---------------------------------------------------------------------------
subroutine makeArtificialDissipation3(grp,ivd)
 ! makes first order artificial dissipation
 IMPLICIT NONE
include 'mpif.h'

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

! a bunch of help variables
 integer :: i,j,i1,i2
 real :: r1,r2,ru1,ru2,u1,u2,h1,h2,p1,p2,wx,wy,rx,ux,vx,tx,ry,uy,vy,ty 
 real :: dr,du,dv,dh,dl,vt,epslm,eps1,gam1,rnx,rny,al,vmo1,rhol,rho1
 real :: uxl,vyl,rhor,uxr,vyr,epsl,presl,hl,epsr,presr,hr,di,d1
 real :: ui,vi,hi,rlam1,rlam2,rlam3,s1,s2,al1x,al2x,cc14,cc24,cc12,cc22
 real :: v1,t1,dr2,du2,dv2,de2,dcf,def,dr4,du4,dv4,de4,d2,d4,d41,d42
 real :: v2,t2,de,dx,dy,df,dp,sp,d0,ci,ci2,af,ucp,dw2

 integer :: kk,istart,ifinish


 ! d2 dissipation only 
grp%dissipation = 0.0 

! assemble smoother

#ifdef PARALLEL

do kk=1,2
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfComSides
 else
  istart = grp%pdp%numberOfComSides + 1
  ifinish = grp%numberOfSides
 end if


 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL


 i1 = grp%sideIndexArray(i,1)
 i2 = grp%sideIndexArray(i,2)

! making second derivative differences
 dr2 = grp%u(i1,1) - grp%u(i2,1)
 du2 = grp%u(i1,2) - grp%u(i2,2)
 dv2 = grp%u(i1,3) - grp%u(i2,3)
 dw2 = grp%u(i1,4) - grp%u(i2,4)
 de2 = grp%u(i1,5) - grp%u(i2,5)        
 ! scalar dissipation coefficient
 d0  = amin1(grp%nodeVolume(i1)/grp%localTimeSteps(i1),&
             grp%nodeVolume(i2)/grp%localTimeSteps(i2))/&
    (grp%nodeConnectivityArray(i1)+grp%nodeConnectivityArray(i2))
 d2  = d0*ivd%coarseGridDissipationFactor
 ! update dissipation matrix 
 grp%dissipation(i1,1) = grp%dissipation(i1,1)+dr2*d2
 grp%dissipation(i1,2) = grp%dissipation(i1,2)+du2*d2
 grp%dissipation(i1,3) = grp%dissipation(i1,3)+dv2*d2
 grp%dissipation(i1,4) = grp%dissipation(i1,4)+dw2*d2
 grp%dissipation(i1,5) = grp%dissipation(i1,5)+de2*d2
 grp%dissipation(i2,1) = grp%dissipation(i2,1)-dr2*d2
 grp%dissipation(i2,2) = grp%dissipation(i2,2)-du2*d2
 grp%dissipation(i2,3) = grp%dissipation(i2,3)-dv2*d2
 grp%dissipation(i2,4) = grp%dissipation(i2,4)-dw2*d2
 grp%dissipation(i2,5) = grp%dissipation(i2,5)-de2*d2    
end do
#ifdef PARALLEL

 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses
   if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
    do j=1,5
     call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%dissipation(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
    call mp_sendv(grp%pdp%processorIDs,i,53,grp%pdp%sendFlag)
   end if
  end do
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,53,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
    do j=1,5
     call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                            grp%pdp%realType,grp%dissipation(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
  end if
 end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position = 1
   do j=1,5
    call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                      grp%pdp%realType,grp%dissipation(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
   end do
   call mp_sendv(grp%pdp%processorIDs,i,54,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,54,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
    do j = 1,5
     call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%dissipation(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
  end if
 end do

#endif PARALLEL


end subroutine makeArtificialDissipation3
!-----------------------------------------------------------------------
 subroutine restrictGrid(finegrp,coarsegrp,ivd,gridnum,calculateTurbulence)
 ! Does a Full Approximation Storage (FAS) multigrid fine to coarse 
 ! transformation
 IMPLICIT NONE

 type(GridSolverData) :: finegrp,coarsegrp
 type(InputVariablesData) :: ivd
 integer :: gridnum
 logical :: calculateTurbulence

!START-TURB
 integer :: i,j,ind,ind2,mvia(7),ib,ip,ist,ien,indc
 real :: fact,mva(7) !,maxval(5)
!END-TURB

 real :: maxVal
 integer :: maxInd

 ! create L_H I_h^H u_h
 call mapFineToCoarse(finegrp%u,coarsegrp,coarsegrp%u)

 call mapScalarFineToCoarse(finegrp%laminarViscosity,coarsegrp,coarsegrp%laminarViscosity)

 call makePressureField(coarsegrp,ivd,coarsegrp%u)
 call setBCsOnSolutionField(coarsegrp,ivd,coarsegrp%u,0)
 call setOuterBoundaryConditions(coarsegrp,ivd,coarsegrp%u,coarsegrp%p)
 coarsegrp%uprev(:,1) = coarsegrp%u(:,1)
 coarsegrp%uprev(:,2) = coarsegrp%u(:,2)
 coarsegrp%uprev(:,3) = coarsegrp%u(:,3)
 coarsegrp%uprev(:,4) = coarsegrp%u(:,4)
 coarsegrp%uprev(:,5) = coarsegrp%u(:,5)
!START-TURB
 if(ivd%turbulenceModel==1) then 
 coarsegrp%uprev(:,6) = coarsegrp%u(:,6)
end if
if(ivd%turbulenceModel==2) then !MOD
 coarsegrp%uprev(:,6) = coarsegrp%u(:,6)
 coarsegrp%uprev(:,7) = coarsegrp%u(:,7)
end if
if(ivd%turbulenceModel==3) then 
 coarsegrp%uprev(:,6) = coarsegrp%u(:,6)
 coarsegrp%uprev(:,7) = coarsegrp%u(:,7)
end if
!END-TURB
 call makePressureField(coarsegrp,ivd,coarsegrp%u)
 coarsegrp%p = abs(coarsegrp%p)
 call makeTimeSteps(coarsegrp,ivd)
 if(ivd%coarseGridDissipationScheme==1) then 
  call makeArtificialDissipation1(coarsegrp,ivd)
 else if(ivd%coarseGridDissipationScheme==2) then
  call makeArtificialDissipation2(coarsegrp,ivd)
 else if(ivd%coarseGridDissipationScheme==3) then
  call makeArtificialDissipation3(coarsegrp,ivd)
 end if


 if(ivd%ReynoldsNumber>0.0) then 
  if(ivd%viscosityScheme==2) then
   call makeViscosityTerm2(coarsegrp,ivd,coarsegrp%dissipation,.true.)
  else
   call makeViscosityTerm(coarsegrp,ivd,coarsegrp%dissipation,.true.)
  end if
 end if


 coarsegrp%rhs = 0.0
 call makeRHS(coarsegrp,ivd)
 call addTimeDifferences(coarsegrp,ivd)
 if(ivd%turbulenceModel==1.and.calculateTurbulence) then 
  call makeTurbulentTimeSteps(coarsegrp,ivd)
  call makeSADiffusion(coarsegrp,ivd)
  call makeSARHS(coarsegrp,ivd)
  call makeSASourceTerm(coarsegrp,ivd)
 end if
!START-TURB
 if(ivd%turbulenceModel==2.and.calculateTurbulence) then 
  call makeTurbulentTimeStepsKW(coarsegrp,ivd)
  call makeKWDiffusion(coarsegrp,ivd)
  call makeKWRHS(coarsegrp,ivd)
  call makeKWSourceTerm(coarsegrp,ivd)
 end if
 if(ivd%turbulenceModel==3.and.calculateTurbulence) then 
  call makeTurbulentTimeStepsKW(coarsegrp,ivd)
  call makeSSTdkdw(coarsegrp,ivd)
  call makeSSTDiffusion(coarsegrp,ivd)
  call makeSSTCrossDiff(coarsegrp,ivd)!Cross-diffusion
  call makeSSTRHS(coarsegrp,ivd)
  call makeSSTSourceTerm(coarsegrp,ivd)
 end if
!END-TURB

 coarsegrp%rhs(:,1) = coarsegrp%rhs(:,1) + coarsegrp%dissipation(:,1)
 coarsegrp%rhs(:,2) = coarsegrp%rhs(:,2) + coarsegrp%dissipation(:,2)
 coarsegrp%rhs(:,3) = coarsegrp%rhs(:,3) + coarsegrp%dissipation(:,3)
 coarsegrp%rhs(:,4) = coarsegrp%rhs(:,4) + coarsegrp%dissipation(:,4)
 coarsegrp%rhs(:,5) = coarsegrp%rhs(:,5) + coarsegrp%dissipation(:,5)
 !START-TURB
 if(ivd%turbulenceModel==1) then
 coarsegrp%rhs(:,6) = coarsegrp%rhs(:,6) + coarsegrp%turbulenceDiffusionTerm(:)
end if
if(ivd%turbulenceModel==2) then !MOD
 coarsegrp%rhs(:,6) = coarsegrp%rhs(:,6) + coarsegrp%turbulenceDiffusionTermK(:)
 coarsegrp%rhs(:,7) = coarsegrp%rhs(:,7) + coarsegrp%turbulenceDiffusionTermW(:)
end if
if(ivd%turbulenceModel==3) then 
 coarsegrp%rhs(:,6) = coarsegrp%rhs(:,6) + coarsegrp%turbulenceDiffusionTermK(:)
 coarsegrp%rhs(:,7) = coarsegrp%rhs(:,7) + coarsegrp%turbulenceDiffusionTermW(:)
end if
!END-TURB

 call solveSystem(coarsegrp,ivd,0,coarsegrp%rhs)
 call setBCsOnIncrementField(coarsegrp,ivd,coarsegrp%rhs)
 call nullifyOuterBoundary(coarsegrp,coarsegrp%rhs)

 call updateSolutionVector(coarsegrp%u,coarsegrp%uprev,coarsegrp%rhs,coarsegrp,ivd,0)

 call setBCsOnSolutionField(coarsegrp,ivd,coarsegrp%u,0) 

 call setOuterBoundaryConditions(coarsegrp,ivd,coarsegrp%u,coarsegrp%p)

 coarsegrp%sourceTerm = coarsegrp%u - coarsegrp%uprev 

 coarsegrp%u = coarsegrp%uprev

 ! create I_h^H L_h u_h
 call makePressureField(finegrp,ivd,finegrp%u)
 finegrp%p = abs(finegrp%p)

 if (ivd%dissipationScheme==1.and.finegrp%gridNumber==1) then 
  call makeArtificialDissipation1(finegrp,ivd)
 else 
  call makeArtificialDissipation2(finegrp,ivd)
 end if

 if(ivd%ReynoldsNumber>0.0) then 
  if(ivd%viscosityScheme==2) then
   call makeViscosityTerm2(finegrp,ivd,finegrp%dissipation,.true.)
  else
   call makeViscosityTerm(finegrp,ivd,finegrp%dissipation,.true.)
  end if
 end if

 finegrp%rhs = 0.0
 call makeRHS(finegrp,ivd)
 call addTimeDifferences(finegrp,ivd)
 if(ivd%turbulenceModel==1.and.calculateTurbulence) then 
  call makeSADiffusion(finegrp,ivd)
  call makeSARHS(finegrp,ivd)
  call makeSASourceTerm(finegrp,ivd)
  if(finegrp%gridNumber==1) then 
   call makeSATripTerm(finegrp,ivd)
  end if
 end if
!START-TURB
 if(ivd%turbulenceModel==2) then !MOD
  call makeKWDiffusion(finegrp,ivd)
  call makeKWRHS(finegrp,ivd)
  call makeKWSourceTerm(finegrp,ivd)
 end if
 if(ivd%turbulenceModel==3) then 
  call makeSSTdkdw(finegrp,ivd)
  call makeSSTDiffusion(finegrp,ivd)
  call makeSSTCrossDiff(finegrp,ivd)!Cross-diffusion
  call makeSSTRHS(finegrp,ivd)
  call makeSSTSourceTerm(finegrp,ivd)
 end if
!END-TURB

 finegrp%rhs(:,1) = finegrp%rhs(:,1) + finegrp%dissipation(:,1)
 finegrp%rhs(:,2) = finegrp%rhs(:,2) + finegrp%dissipation(:,2)
 finegrp%rhs(:,3) = finegrp%rhs(:,3) + finegrp%dissipation(:,3)
 finegrp%rhs(:,4) = finegrp%rhs(:,4) + finegrp%dissipation(:,4)
 finegrp%rhs(:,5) = finegrp%rhs(:,5) + finegrp%dissipation(:,5)
!START-TURB
 if(ivd%turbulenceModel==1) then 
 finegrp%rhs(:,6) = finegrp%rhs(:,6) + finegrp%turbulenceDiffusionTerm(:)
end if
if(ivd%turbulenceModel==2) then 
 finegrp%rhs(:,6) = finegrp%rhs(:,6) + finegrp%turbulenceDiffusionTermK(:)
 finegrp%rhs(:,7) = finegrp%rhs(:,7) + finegrp%turbulenceDiffusionTermW(:)
end if
if(ivd%turbulenceModel==3) then 
 finegrp%rhs(:,6) = finegrp%rhs(:,6) + finegrp%turbulenceDiffusionTermK(:)
 finegrp%rhs(:,7) = finegrp%rhs(:,7) + finegrp%turbulenceDiffusionTermW(:)
end if
!END-TURB 

 call solveSystem(finegrp,ivd,0,finegrp%rhs)
 if(associated(finegrp%sourceTerm)) then 
  call setBCsOnIncrementField(finegrp,ivd,finegrp%rhs) ! CHANGED 01.09.99 
  call nullifyOuterBoundary(finegrp,finegrp%rhs)
 else
  call setBCsOnIncrementField(finegrp,ivd,finegrp%rhs) 
 end if

 call updateSolutionVector(finegrp%rhs,finegrp%uprev,finegrp%rhs,finegrp,ivd,0)

 call setBCsOnSolutionField(finegrp,ivd,finegrp%rhs,0)
 if(finegrp%gridNumber>1) then 
  finegrp%rhs(:,1) = finegrp%rhs(:,1) - finegrp%sourceTerm(:,1)
  finegrp%rhs(:,2) = finegrp%rhs(:,2) - finegrp%sourceTerm(:,2)
  finegrp%rhs(:,3) = finegrp%rhs(:,3) - finegrp%sourceTerm(:,3)
  finegrp%rhs(:,4) = finegrp%rhs(:,4) - finegrp%sourceTerm(:,4)
  finegrp%rhs(:,5) = finegrp%rhs(:,5) - finegrp%sourceTerm(:,5)
!START-TURB
 if(ivd%turbulenceModel==1) then
  finegrp%rhs(:,6) = finegrp%rhs(:,6) - finegrp%sourceTerm(:,6)
 end if
 if(ivd%turbulenceModel==2) then !MOD
  finegrp%rhs(:,6) = finegrp%rhs(:,6) - finegrp%sourceTerm(:,6)
  finegrp%rhs(:,7) = finegrp%rhs(:,7) - finegrp%sourceTerm(:,7)
 end if
 if(ivd%turbulenceModel==3) then 
  finegrp%rhs(:,6) = finegrp%rhs(:,6) - finegrp%sourceTerm(:,6)
  finegrp%rhs(:,7) = finegrp%rhs(:,7) - finegrp%sourceTerm(:,7)
 end if
!END-TURB
 end if

!START-TURB
 if(ivd%turbulenceModel==1) then 
  call limitTurbulence(finegrp,finegrp%rhs(:,6),ivd)
 end if
!END-TURB

 call setOuterBoundaryConditions(finegrp,ivd,finegrp%rhs,finegrp%p)
 finegrp%rhs(:,1) = finegrp%rhs(:,1) - finegrp%uprev(:,1)
 finegrp%rhs(:,2) = finegrp%rhs(:,2) - finegrp%uprev(:,2)
 finegrp%rhs(:,3) = finegrp%rhs(:,3) - finegrp%uprev(:,3)
 finegrp%rhs(:,4) = finegrp%rhs(:,4) - finegrp%uprev(:,4)
 finegrp%rhs(:,5) = finegrp%rhs(:,5) - finegrp%uprev(:,5)
 !START-TURB
 if(ivd%turbulenceModel==1) then
 finegrp%rhs(:,6) = finegrp%rhs(:,6) - finegrp%uprev(:,6)
 end if
 if(ivd%turbulenceModel==2) then 
 finegrp%rhs(:,6) = finegrp%rhs(:,6) - finegrp%uprev(:,6)
 finegrp%rhs(:,7) = finegrp%rhs(:,7) - finegrp%uprev(:,7)
 end if
 if(ivd%turbulenceModel==3) then 
 finegrp%rhs(:,6) = finegrp%rhs(:,6) - finegrp%uprev(:,6)
 finegrp%rhs(:,7) = finegrp%rhs(:,7) - finegrp%uprev(:,7)
 end if
!END-TURB

 call mapFineToCoarse(finegrp%rhs,coarsegrp,coarsegrp%rhs)

 call setBCsOnIncrementField(coarsegrp,ivd,coarsegrp%rhs)  !KKK


 coarsegrp%sourceTerm(:,1) = coarsegrp%sourceTerm(:,1)-coarsegrp%rhs(:,1)
 coarsegrp%sourceTerm(:,2) = coarsegrp%sourceTerm(:,2)-coarsegrp%rhs(:,2)
 coarsegrp%sourceTerm(:,3) = coarsegrp%sourceTerm(:,3)-coarsegrp%rhs(:,3)
 coarsegrp%sourceTerm(:,4) = coarsegrp%sourceTerm(:,4)-coarsegrp%rhs(:,4)
 coarsegrp%sourceTerm(:,5) = coarsegrp%sourceTerm(:,5)-coarsegrp%rhs(:,5)
!START-TURB
if(ivd%turbulenceModel==1) then
 coarsegrp%sourceTerm(:,6) = coarsegrp%sourceTerm(:,6)-coarsegrp%rhs(:,6)
end if
if(ivd%turbulenceModel==2) then
 coarsegrp%sourceTerm(:,6) = coarsegrp%sourceTerm(:,6)-coarsegrp%rhs(:,6)
 coarsegrp%sourceTerm(:,7) = coarsegrp%sourceTerm(:,7)-coarsegrp%rhs(:,7)
end if
if(ivd%turbulenceModel==3) then 
 coarsegrp%sourceTerm(:,6) = coarsegrp%sourceTerm(:,6)-coarsegrp%rhs(:,6)
 coarsegrp%sourceTerm(:,7) = coarsegrp%sourceTerm(:,7)-coarsegrp%rhs(:,7)
end if
!END-TURB

 ! As a smart move, initialize coarse grid vector using second part of source term
 coarsegrp%u(:,1) = coarsegrp%u(:,1) + coarsegrp%rhs(:,1) 
 coarsegrp%u(:,2) = coarsegrp%u(:,2) + coarsegrp%rhs(:,2) 
 coarsegrp%u(:,3) = coarsegrp%u(:,3) + coarsegrp%rhs(:,3) 
 coarsegrp%u(:,4) = coarsegrp%u(:,4) + coarsegrp%rhs(:,4) 
 coarsegrp%u(:,5) = coarsegrp%u(:,5) + coarsegrp%rhs(:,5) 
!START-TURB
if(ivd%turbulenceModel==1) then
 coarsegrp%u(:,6) = coarsegrp%u(:,6) + coarsegrp%rhs(:,6) 
end if
if(ivd%turbulenceModel==2) then 
 coarsegrp%u(:,6) = coarsegrp%u(:,6) + coarsegrp%rhs(:,6) 
 coarsegrp%u(:,7) = coarsegrp%u(:,7) + coarsegrp%rhs(:,7)
end if
if(ivd%turbulenceModel==3) then 
 coarsegrp%u(:,6) = coarsegrp%u(:,6) + coarsegrp%rhs(:,6) 
 coarsegrp%u(:,7) = coarsegrp%u(:,7) + coarsegrp%rhs(:,7)
end if
!END-TURB

 end subroutine restrictGrid
!-----------------------------------------------------------------------
 subroutine prolongateGrid(coarsegrp,finegrp,ivd) 
 ! Does a Full Approximation Storage (FAS) multigrid coarse to fine
 ! transformation
 IMPLICIT NONE

 type(GridSolverData) :: coarsegrp,finegrp
 type(InputVariablesData) :: ivd
!START-TURB
 integer :: i,j,ind,ind2,maxind(7),numbComp,k
!END-TURB

  real :: maxRes

 ! first get the fine grid approximation mapped to the coarse grid 
 call mapFineToCoarse(finegrp%u,coarsegrp,coarsegrp%uprev)

 call setBCsOnSolutionField(coarsegrp,ivd,coarsegrp%uprev,0)
 call setOuterBoundaryConditions(coarsegrp,ivd,coarsegrp%uprev,coarsegrp%p)

 ! now create the coarse grid correction v_H = u_H-I_h^H u_h
 coarsegrp%rhs(:,1) = coarsegrp%u(:,1)  - coarsegrp%uprev(:,1)
 coarsegrp%rhs(:,2) = coarsegrp%u(:,2)  - coarsegrp%uprev(:,2)
 coarsegrp%rhs(:,3) = coarsegrp%u(:,3)  - coarsegrp%uprev(:,3)
 coarsegrp%rhs(:,4) = coarsegrp%u(:,4)  - coarsegrp%uprev(:,4)
 coarsegrp%rhs(:,5) = coarsegrp%u(:,5)  - coarsegrp%uprev(:,5)
!START-TURB
if(ivd%turbulenceModel==1) then
 coarsegrp%rhs(:,6) = coarsegrp%u(:,6)  - coarsegrp%uprev(:,6)
end if
if(ivd%turbulenceModel==2) then 
 coarsegrp%rhs(:,6) = coarsegrp%u(:,6)  - coarsegrp%uprev(:,6)
 coarsegrp%rhs(:,7) = coarsegrp%u(:,7)  - coarsegrp%uprev(:,7)
end if
if(ivd%turbulenceModel==3) then 
 coarsegrp%rhs(:,6) = coarsegrp%u(:,6)  - coarsegrp%uprev(:,6)
 coarsegrp%rhs(:,7) = coarsegrp%u(:,7)  - coarsegrp%uprev(:,7)
end if
!END-TURB

 call nullifyOuterBoundary(coarsegrp,coarsegrp%rhs)

 call setBCsOnIncrementField(coarsegrp,ivd,coarsegrp%rhs)  !KKK

 call mapCoarseToFine(coarsegrp%rhs,coarsegrp,finegrp,finegrp%rhs,ivd)

 if(ivd%prolongationSmoothingFactor>0.0) then 
  do i=1,ivd%numberOfPSSteps
   call smoothResidual2(finegrp,ivd,finegrp%rhs,ivd%prolongationSmoothingFactor)
  end do
 end if

 call nullifyOuterBoundary(finegrp,finegrp%rhs)

 finegrp%rhs(:,1) = ivd%prolongationRelaxation*finegrp%rhs(:,1)
 finegrp%rhs(:,2) = ivd%prolongationRelaxation*finegrp%rhs(:,2)
 finegrp%rhs(:,3) = ivd%prolongationRelaxation*finegrp%rhs(:,3)
 finegrp%rhs(:,4) = ivd%prolongationRelaxation*finegrp%rhs(:,4)
 finegrp%rhs(:,5) = ivd%prolongationRelaxation*finegrp%rhs(:,5)
!START-TURB
 if(ivd%turbulenceModel==1) then
 finegrp%rhs(:,6) = ivd%turbProlongationRelaxation*finegrp%rhs(:,6)
 end if
if(ivd%turbulenceModel==2) then 
 finegrp%rhs(:,6) = ivd%turbProlongationRelaxation*finegrp%rhs(:,6)
 finegrp%rhs(:,7) = ivd%turbProlongationRelaxation*finegrp%rhs(:,7)
end if
if(ivd%turbulenceModel==3) then 
 finegrp%rhs(:,6) = ivd%turbProlongationRelaxation*finegrp%rhs(:,6)
 finegrp%rhs(:,7) = ivd%turbProlongationRelaxation*finegrp%rhs(:,7)
end if
!END-TURB

 finegrp%rhs(1:finegrp%brp%numberOfBoundaryNodes,1) =&
      ivd%multigridBCRelaxation*finegrp%rhs(1:finegrp%brp%numberOfBoundaryNodes,1)
 finegrp%rhs(1:finegrp%brp%numberOfBoundaryNodes,2) =&
      ivd%multigridBCRelaxation*finegrp%rhs(1:finegrp%brp%numberOfBoundaryNodes,2)
 finegrp%rhs(1:finegrp%brp%numberOfBoundaryNodes,3) =&
      ivd%multigridBCRelaxation*finegrp%rhs(1:finegrp%brp%numberOfBoundaryNodes,3)
 finegrp%rhs(1:finegrp%brp%numberOfBoundaryNodes,4) =&
      ivd%multigridBCRelaxation*finegrp%rhs(1:finegrp%brp%numberOfBoundaryNodes,4)
 finegrp%rhs(1:finegrp%brp%numberOfBoundaryNodes,5) =&
      ivd%multigridBCRelaxation*finegrp%rhs(1:finegrp%brp%numberOfBoundaryNodes,5)

!START-TURB
if(ivd%turbulenceModel==1) then
 if(finegrp%gridNumber==1) then
  do i=1,finegrp%numberOfTripFieldNodes
   ind = finegrp%tripNodeFieldIndexes(i,1)
   if(finegrp%tripNodeFieldDistances(i)<ivd%triggerRadius) then
    finegrp%rhs(ind,6) = 0.0
   end if
  end do
 end if
end if
!END-TURB

 if(finegrp%gridNumber<ivd%numberOfTurbulenceGrids) then 
!START-TURB
  if(ivd%turbulenceModel==1) then
   numbComp = 6
  end if
  if(ivd%turbulenceModel==2 .or. ivd%turbulenceModel==3) then
   numbComp = 7
  end if
 else
   numbComp = 5
 end if
!END-TURB

 do i=1,finegrp%numberOfNodes
  do j=1,numbComp
   finegrp%u(i,j) = finegrp%u(i,j)+finegrp%rhs(i,j)
  end do 
 end do


 call setBCsOnSolutionField(finegrp,ivd,finegrp%u,0)
 call setOuterBoundaryConditions(finegrp,ivd,finegrp%u,finegrp%p)


 do i=1,finegrp%numberOfNodes
!START-TURB
  if(ivd%turbulenceModel==1) then
   finegrp%u(i,6) = max(finegrp%u(i,6),0.0)
  end if
  if(ivd%turbulenceModel==2) then
   finegrp%u(i,6) = max(finegrp%u(i,6),0.0) 
   finegrp%u(i,7) = max(finegrp%u(i,7),ivd%minTurbulenceValueW)
  end if
  if(ivd%turbulenceModel==3) then 
   finegrp%u(i,6) = max(finegrp%u(i,6),0.0) 
   finegrp%u(i,7) = max(finegrp%u(i,7),ivd%minTurbulenceValueW)
  end if
!END-TURB
  finegrp%u(i,1) = max(finegrp%u(i,1),ivd%minimumDensity)
 end do


 end subroutine prolongateGrid 
!-----------------------------------------------------------------------
 subroutine mapFineToCoarse(fineu,coarsegrp,u)
 ! maps fine grid unknown fields to coarse grid
 IMPLICIT NONE 
include 'mpif.h'

 real :: fineu(:,:)
 type(GridSolverData) :: coarsegrp
 real :: u(:,:) 
 
 integer :: i,j,k,currentInteger,numComp
 real :: currentReal

 integer :: ssf

 integer :: ind,indc

 ! initialize coarse grid
 u = 0.0

 do i=1,coarsegrp%rod%numberOfFineNodes
  if(coarsegrp%pod%prolongationArray(i)>0) then 
   u(coarsegrp%pod%prolongationArray(i),:) = u(coarsegrp%pod%prolongationArray(i),:) +& 
                   coarsegrp%rod%restrictionArray(i)*fineu(i,:)
  end if
 end do

! com nodes

 do i=1,coarsegrp%rod%numberOfComNodes
  coarsegrp%mapComBuffer(i,:) = coarsegrp%rod%comRestrictionArray(i)*fineu(coarsegrp%pod%comProlongationArray(i),:)
 end do

#ifdef PARALLEL

 numComp = size(u,2)

 do i=1,coarsegrp%pdp%numberOfProcesses
  if (i.ne.coarsegrp%pdp%currentDomain.and.(coarsegrp%pdp%numberOfComSendNodes(i)>0)) then
   call mp_init_buffer(coarsegrp%pdp%processorIDs,i,coarsegrp%pdp%sendFlag)
   position = 1
   do j=1,numComp
    call mp_pakv_bdry(coarsegrp%pdp%processorIDs,i,coarsegrp%pdp%numberOfComSendNodes(i),coarsegrp%pdp%comSendRegister(i,:),&
                  coarsegrp%pdp%realType,coarsegrp%mapComBuffer(:,j),coarsegrp%pdp%comBuffer,coarsegrp%pdp%sendFlag)
   end do
   call mp_sendv(coarsegrp%pdp%processorIDs,i,778,coarsegrp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(coarsegrp%pdp%numberOfProcesses,coarsegrp%pdp%numberOfComReceiveNodes,&
                     coarsegrp%pdp%receivedMessageFromProcess)

  do i  = 1, coarsegrp%pdp%numberOfProcesses
    if(i.ne.coarsegrp%pdp%currentDomain.and.(coarsegrp%pdp%numberOfComReceiveNodes(i) > 0 ) ) then
      call mp_recv( coarsegrp%pdp%processorIDs,i,778,coarsegrp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( coarsegrp%pdp%sendFlag )

 do i=1,coarsegrp%pdp%numberOfProcesses
    if(i.ne.coarsegrp%pdp%currentDomain.and.(coarsegrp%pdp%numberOfComReceiveNodes(i) > 0 ) ) then
    do j = 1,numComp
     call mp_upakv_add_bdry(coarsegrp%pdp%processorIDs,i,coarsegrp%pdp%numberOfComReceiveNodes(i),coarsegrp%pdp%comReceiveRegister(i,:),&
                        coarsegrp%pdp%realType,u(:,j),coarsegrp%pdp%comBuffer,coarsegrp%pdp%sendFlag)
    end do

  end if
 end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

!  write(get_unit(),*)  'Starting second F2C comms'
 do i=1,coarsegrp%pdp%numberOfProcesses
  if (i.ne.coarsegrp%pdp%currentDomain.and.(coarsegrp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(coarsegrp%pdp%processorIDs,i,coarsegrp%pdp%sendFlag)
   position = 1
   do j=1,numComp
    call mp_pakv_bdry(coarsegrp%pdp%processorIDs,i,coarsegrp%pdp%numberOfReceiveNodes(i),coarsegrp%pdp%receiveNodeRegister(i,:),&
                      coarsegrp%pdp%realType,u(:,j),coarsegrp%pdp%buffer,coarsegrp%pdp%sendFlag)
   end do
   call mp_sendv(coarsegrp%pdp%processorIDs,i,78,coarsegrp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(coarsegrp%pdp%numberOfProcesses,coarsegrp%pdp%numberOfSendNodes,&
                     coarsegrp%pdp%receivedMessageFromProcess)

  do i  = 1, coarsegrp%pdp%numberOfProcesses
    if(i.ne.coarsegrp%pdp%currentDomain.and.(coarsegrp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( coarsegrp%pdp%processorIDs,i,78,coarsegrp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( coarsegrp%pdp%sendFlag )

 do i=1,coarsegrp%pdp%numberOfProcesses
    if(i.ne.coarsegrp%pdp%currentDomain.and.(coarsegrp%pdp%numberOfSendNodes(i) > 0 ) ) then
    do j = 1,numComp
     call mp_upakv_nabr(coarsegrp%pdp%processorIDs,i,coarsegrp%pdp%numberOfSendNodes(i),coarsegrp%pdp%sendNodeRegister(i,:),&
                        coarsegrp%pdp%realType,u(:,j),coarsegrp%pdp%buffer,coarsegrp%pdp%sendFlag)
    end do
  end if
 end do

#endif PARALLEL


 end subroutine mapFineToCoarse
!-----------------------------------------------------------------------
 subroutine mapCoarseToFine(ucoarse,coarsegrp,finegrp,u,ivd)
 ! maps coarse grid fields to fine grid
 ! here simple injection is used
 IMPLICIT NONE

 type(GridSolverData) :: coarsegrp,finegrp
 real :: ucoarse(:,:),u(:,:)
 type(InputVariablesData) :: ivd

 integer :: i,j,currentInteger,i1,i2,numComp
 real :: x0(3),x1(3),x2(3),x3(3),xp(3),coeff0
 real :: u0(5),u1(5),u2(5),u3(5),oneOverDet,val,um(5)

 if(ivd%prolongationMappingScheme==1) then 
  ! injectional prolongation
  do i=1,coarsegrp%pod%numberOfFineNodes
   if(coarsegrp%pod%prolongationArray(i)>0) then 
    u(i,:) = ucoarse(coarsegrp%pod%prolongationArray(i),:)
   end if
  end do
 end if


#ifdef PARALLEL

 numComp = size(u,2)

 do i=1,coarsegrp%pdp%numberOfProcesses
  if (i.ne.coarsegrp%pdp%currentDomain.and.(coarsegrp%pdp%numberOfComReceiveNodes(i)>0)) then
   call mp_init_buffer(coarsegrp%pdp%processorIDs,i,coarsegrp%pdp%sendFlag)
   position = 1
   do j=1,numComp
    call mp_pakv_bdry(coarsegrp%pdp%processorIDs,i,coarsegrp%pdp%numberOfComReceiveNodes(i),coarsegrp%pdp%comReceiveRegister(i,:),&
                      coarsegrp%pdp%realType,ucoarse(:,j),coarsegrp%pdp%comBuffer,coarsegrp%pdp%sendFlag)
   end do
   call mp_sendv(coarsegrp%pdp%processorIDs,i,779,coarsegrp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(coarsegrp%pdp%numberOfProcesses,coarsegrp%pdp%numberOfComSendNodes,&
                     coarsegrp%pdp%receivedMessageFromProcess)

  do i  = 1, coarsegrp%pdp%numberOfProcesses
    if(i.ne.coarsegrp%pdp%currentDomain.and.(coarsegrp%pdp%numberOfComSendNodes(i) > 0 ) ) then
      call mp_recv( coarsegrp%pdp%processorIDs,i,779,coarsegrp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( coarsegrp%pdp%sendFlag )

 do i=1,coarsegrp%pdp%numberOfProcesses
    if(i.ne.coarsegrp%pdp%currentDomain.and.(coarsegrp%pdp%numberOfComSendNodes(i) > 0 ) ) then
    do j = 1,numComp
     call mp_upakv_nabr(coarsegrp%pdp%processorIDs,i,coarsegrp%pdp%numberOfComSendNodes(i),coarsegrp%pdp%comSendRegister(i,:),&
                        coarsegrp%pdp%realType,coarsegrp%mapComBuffer(:,j),coarsegrp%pdp%comBuffer,coarsegrp%pdp%sendFlag)
    end do
  end if
 end do

#endif PARALLEL

! com nodes

 do i=1,coarsegrp%rod%numberOfComNodes
  u(coarsegrp%pod%comProlongationArray(i),:) = coarsegrp%mapComBuffer(i,:)
 end do

#ifdef PARALLEL
 do i=1,finegrp%pdp%numberOfProcesses
  if (i.ne.finegrp%pdp%currentDomain.and.(finegrp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(finegrp%pdp%processorIDs,i,finegrp%pdp%sendFlag)
   position = 1
   do j=1,numComp
    call mp_pakv_bdry(finegrp%pdp%processorIDs,i,finegrp%pdp%numberOfReceiveNodes(i),finegrp%pdp%receiveNodeRegister(i,:),&
                      finegrp%pdp%realType,u(:,j),finegrp%pdp%buffer,finegrp%pdp%sendFlag)
   end do
   call mp_sendv(finegrp%pdp%processorIDs,i,79,finegrp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(finegrp%pdp%numberOfProcesses,finegrp%pdp%numberOfSendNodes,&
                     finegrp%pdp%receivedMessageFromProcess)

  do i  = 1, finegrp%pdp%numberOfProcesses
    if(i.ne.finegrp%pdp%currentDomain.and.(finegrp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( finegrp%pdp%processorIDs,i,79,finegrp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( finegrp%pdp%sendFlag )

 do i=1,finegrp%pdp%numberOfProcesses
    if(i.ne.finegrp%pdp%currentDomain.and.(finegrp%pdp%numberOfSendNodes(i) > 0 ) ) then
    do j = 1,numComp
     call mp_upakv_nabr(finegrp%pdp%processorIDs,i,finegrp%pdp%numberOfSendNodes(i),finegrp%pdp%sendNodeRegister(i,:),&
                        finegrp%pdp%realType,u(:,j),finegrp%pdp%buffer,finegrp%pdp%sendFlag)
    end do
  end if
 end do

#endif PARALLEL

 end subroutine mapCoarseToFine
!-----------------------------------------------------------------------
 subroutine mapScalarFineToCoarse(fineu,coarsegrp,u)
 ! maps fine grid unknown fields to coarse grid
 IMPLICIT NONE
include 'mpif.h'

 real :: fineu(:)
 type(GridSolverData) :: coarsegrp
 real :: u(:)

 integer :: i,j,currentInteger,numComp
 real :: currentReal

 ! initialize coarse grid
 u = 0.0

 do i=1,coarsegrp%rod%numberOfFineNodes
  if(coarsegrp%pod%prolongationArray(i)>0) then
   u(coarsegrp%pod%prolongationArray(i)) = u(coarsegrp%pod%prolongationArray(i)) +&
                   coarsegrp%rod%restrictionArray(i)*fineu(i)
  end if
 end do

! com nodes

 do i=1,coarsegrp%rod%numberOfComNodes
  coarsegrp%mapComBuffer(i,1) = coarsegrp%rod%comRestrictionArray(i)*fineu(coarsegrp%pod%comProlongationArray(i))
 end do

#ifdef PARALLEL

 do i=1,coarsegrp%pdp%numberOfProcesses
  if (i.ne.coarsegrp%pdp%currentDomain.and.(coarsegrp%pdp%numberOfComSendNodes(i)>0)) then
   call mp_init_buffer(coarsegrp%pdp%processorIDs,i,coarsegrp%pdp%sendFlag)
   position = 1
   call mp_pakv_bdry(coarsegrp%pdp%processorIDs,i,coarsegrp%pdp%numberOfComSendNodes(i),coarsegrp%pdp%comSendRegister(i,:),&
                 coarsegrp%pdp%realType,coarsegrp%mapComBuffer(:,1),coarsegrp%pdp%comBuffer,coarsegrp%pdp%sendFlag)
   call mp_sendv(coarsegrp%pdp%processorIDs,i,878,coarsegrp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(coarsegrp%pdp%numberOfProcesses,coarsegrp%pdp%numberOfComReceiveNodes,&
                     coarsegrp%pdp%receivedMessageFromProcess)

  do i  = 1, coarsegrp%pdp%numberOfProcesses
    if(i.ne.coarsegrp%pdp%currentDomain.and.(coarsegrp%pdp%numberOfComReceiveNodes(i) > 0 ) ) then
      call mp_recv( coarsegrp%pdp%processorIDs,i,878,coarsegrp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( coarsegrp%pdp%sendFlag )

 do i=1,coarsegrp%pdp%numberOfProcesses
    if(i.ne.coarsegrp%pdp%currentDomain.and.(coarsegrp%pdp%numberOfComReceiveNodes(i) > 0 ) ) then
    call mp_upakv_add_bdry(coarsegrp%pdp%processorIDs,i,coarsegrp%pdp%numberOfComReceiveNodes(i),coarsegrp%pdp%comReceiveRegister(i,:),&
                       coarsegrp%pdp%realType,u,coarsegrp%pdp%comBuffer,coarsegrp%pdp%sendFlag)
  end if
 end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

!  write(get_unit(),*)  'Starting second SF2C comms'
 do i=1,coarsegrp%pdp%numberOfProcesses
  if (i.ne.coarsegrp%pdp%currentDomain.and.(coarsegrp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(coarsegrp%pdp%processorIDs,i,coarsegrp%pdp%sendFlag)
   position = 1
   call mp_pakv_bdry(coarsegrp%pdp%processorIDs,i,coarsegrp%pdp%numberOfReceiveNodes(i),coarsegrp%pdp%receiveNodeRegister(i,:),&
                     coarsegrp%pdp%realType,u,coarsegrp%pdp%buffer,coarsegrp%pdp%sendFlag)
   call mp_sendv(coarsegrp%pdp%processorIDs,i,88,coarsegrp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(coarsegrp%pdp%numberOfProcesses,coarsegrp%pdp%numberOfSendNodes,&
                     coarsegrp%pdp%receivedMessageFromProcess)

  do i  = 1, coarsegrp%pdp%numberOfProcesses
    if(i.ne.coarsegrp%pdp%currentDomain.and.(coarsegrp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( coarsegrp%pdp%processorIDs,i,88,coarsegrp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( coarsegrp%pdp%sendFlag )

 do i=1,coarsegrp%pdp%numberOfProcesses
    if(i.ne.coarsegrp%pdp%currentDomain.and.(coarsegrp%pdp%numberOfSendNodes(i) > 0 ) ) then
    call mp_upakv_nabr(coarsegrp%pdp%processorIDs,i,coarsegrp%pdp%numberOfSendNodes(i),coarsegrp%pdp%sendNodeRegister(i,:),&
                       coarsegrp%pdp%realType,u,coarsegrp%pdp%buffer,coarsegrp%pdp%sendFlag)
  end if
 end do

#endif PARALLEL

 end subroutine mapScalarFineToCoarse
!-----------------------------------------------------------------------
 subroutine makeTurbulentTimeSteps(grp,ivd)
 ! make local time steps for SA turbulence equation
 IMPLICIT NONE
include 'mpif.h'

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd


 real :: S,dt,wx,wy,wz,al,ro1,ro2,vn1,vn2,vis1,vis2,dv1,dv2
 real :: oneOverSigma,cb2,dvm,vnm,dInvSquared,oneOverKappaSquared
 real :: cb1,cw1,ct3,oneOverReynoldsNumber

 integer :: i,istart,ifinish,kk,i1,i2

 real :: avRatio,minRatio,maxRatio
 integer :: maxRatioInd,minRatioInd

 oneOverSigma = 3.0/2.0
 oneOverKappaSquared = 1.0/(0.41*0.41)
 cb1 = 0.1355
 cb2 = 0.622
 ct3 = 1.1
 cw1 = cb1*oneOverKappaSquared + (1.0+cb2)*oneOverSigma

 oneOverReynoldsNumber = 1.0/ivd%ReynoldsNumber 

 grp%localTurbulentTimeSteps = 0.0
 grp%nodeHelpArray(:,1) = 0.0


! run over sides
#ifdef PARALLEL

do kk=1,2
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfComSides
 else
  istart = grp%pdp%numberOfComSides + 1
  ifinish = grp%numberOfSides
 end if


 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL

  i1  = grp%sideIndexArray(i,1)
  i2  = grp%sideIndexArray(i,2)
  wx  = grp%sideWeightsArray(i,1)
  wy  = grp%sideWeightsArray(i,2)
  wz  = grp%sideWeightsArray(i,3)
  al  = wx*wx + wy*wy + wz*wz
  ro1 = 1./grp%u(i1,1)
  ro2 = 1./grp%u(i2,1)
  vn1 = (wx*grp%u(i1,2) + wy*grp%u(i1,3) + wz*grp%u(i1,4) )*ro1
  vn2 = (wx*grp%u(i2,2) + wy*grp%u(i2,3) + wz*grp%u(i2,4) )*ro2

  vn1 = ivd%turbK*sqrt(al) + abs(vn1)
  vn2 = ivd%turbK*sqrt(al) + abs(vn2)

  grp%localTurbulentTimeSteps(i1) = grp%localTurbulentTimeSteps(i1) + vn1
  grp%localTurbulentTimeSteps(i2) = grp%localTurbulentTimeSteps(i2) + vn2 

  vis1 = grp%laminarViscosity(i1)
  vis2 = grp%laminarViscosity(i2)

  vis1 = oneOverSigma*(vis1+3.0*cb2*grp%u(i1,6)*oneOverReynoldsNumber)
  vis2 = oneOverSigma*(vis2+3.0*cb2*grp%u(i2,6)*oneOverReynoldsNumber)


  dv1 = 2.0*vis1*al
  dv2 = 2.0*vis2*al
  grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + dv2
  grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + dv1
 end do

! run over faces

#ifdef PARALLEL
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfBoundaryComSides
 else
  istart = grp%pdp%numberOfBoundaryComSides + 1
  ifinish = grp%brp%numberOfBoundarySides
 end if
 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL

  i1  = grp%brp%sideIndexArray(i,1)
  i2  = grp%brp%sideIndexArray(i,2)
  wx  = grp%brp%sideWeightsArray(i,1)
  wy  = grp%brp%sideWeightsArray(i,2)
  wz  = grp%brp%sideWeightsArray(i,3)
  al  = wx*wx + wy*wy + wz*wz
  ro1 = 1./grp%u(i1,1)
  ro2 = 1./grp%u(i2,1)
  vn1 = (wx*grp%u(i1,2) + wy*grp%u(i1,3) + wz*grp%u(i1,4))*ro1
  vn2 = (wx*grp%u(i2,2) + wy*grp%u(i2,3) + wz*grp%u(i2,4))*ro2

  vn1 = ivd%turbK*sqrt(al) + abs(vn1)
  vn2 = ivd%turbK*sqrt(al) + abs(vn2)

  vnm = vn1 + vn2

  grp%localTurbulentTimeSteps(i1) = grp%localTurbulentTimeSteps(i1) + vn1 + vnm 
  grp%localTurbulentTimeSteps(i2) = grp%localTurbulentTimeSteps(i2) + vn2 + vnm 

  vis1 = grp%laminarViscosity(i1)
  vis2 = grp%laminarViscosity(i2)

  vis1 = oneOverSigma*(vis1+3.0*cb2*grp%u(i1,6)*oneOverReynoldsNumber)
  vis2 = oneOverSigma*(vis2+3.0*cb2*grp%u(i2,6)*oneOverReynoldsNumber)

  dv1 = 2.0*vis1*al
  dv2 = 2.0*vis2*al
  dvm = dv1 + dv2
  grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + dv1 + dvm
  grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + dv2 + dvm
 end do
#ifdef PARALLEL

 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses
   if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
    call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%localTurbulentTimeSteps(:),grp%pdp%buffer,grp%pdp%sendFlag)
    call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
    call mp_sendv(grp%pdp%processorIDs,i,355,grp%pdp%sendFlag)
   end if
  end do
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,355,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
    call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                           grp%pdp%realType,grp%localTurbulentTimeSteps(:),grp%pdp%buffer,grp%pdp%sendFlag)
    call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                           grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
  end if
 end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position = 1
   call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                     grp%pdp%realType,grp%localTurbulentTimeSteps(:),grp%pdp%buffer,grp%pdp%sendFlag)
   call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                     grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
   call mp_sendv(grp%pdp%processorIDs,i,356,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,356,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )


 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
    call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%localTurbulentTimeSteps(:),grp%pdp%buffer,grp%pdp%sendFlag)
    call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
  end if
 end do

#endif PARALLEL


 do i = 1,grp%numberOfNodes
  grp%localTurbulentTimeSteps(i) = 0.5/&
    (grp%localTurbulentTimeSteps(i)/grp%nodeVolume(i)&
     +grp%nodeHelpArray(i,1)/(grp%nodeVolume(i)**2))
  dt = grp%localTurbulentTimeSteps(i)
  if(grp%wallDistance(i)>1.0e-10) then 
   dInvSquared = 1.0/grp%wallDistance(i)
  else
   dInvSquared = 1.0e10
  end if
  dInvSquared = dInvSquared*dInvSquared
  S = abs(grp%divergence(i))+cb1*(grp%vorticity(i)+&
      (1.0+ct3)*grp%u(i,6)*oneOverKappaSquared*dInvSquared*oneOverReynoldsNumber)+&
      cw1*grp%u(i,6)*dInvSquared*oneOverReynoldsNumber
  grp%localTurbulentTimeSteps(i) = dt/(1.0+0.5*S*dt)
 end do

 end subroutine makeTurbulentTimeSteps
!-----------------------------------------------------------------------
!START-TURB
subroutine makeTurbulentTimeStepsKW(grp,ivd)
 ! make local time steps for KW turbulence equation
 ! Time stepping based on MacCormack(1982) and Soltani (1991)
 IMPLICIT NONE
include 'mpif.h'

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd


 real :: wx,wy,wz,al,ro1,ro2,vn1,vn2,vis1,vis2,dv1,dv2
 real :: dvm,vnm,oneOverReynoldsNumber,oneOverMachNumber

 integer :: i,istart,ifinish,kk,i1,i2

 real :: avRatio,minRatio,maxRatio,gammaPr
 integer :: maxRatioInd,minRatioInd

 oneOverReynoldsNumber = 1.0/ivd%ReynoldsNumber 
 oneOverMachNumber     = 1.0/ivd%MachNumber
 gammaPr = max(2.0,ivd%gamma/ivd%PrandtlNumber)

 grp%localTurbulentTimeSteps = 0.0
 grp%nodeHelpArray(:,1) = 0.0


! run over sides
#ifdef PARALLEL

do kk=1,2
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfComSides
 else
  istart = grp%pdp%numberOfComSides + 1
  ifinish = grp%numberOfSides
 end if


 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL

  i1  = grp%sideIndexArray(i,1)
  i2  = grp%sideIndexArray(i,2)
  wx  = grp%sideWeightsArray(i,1)
  wy  = grp%sideWeightsArray(i,2)
  wz  = grp%sideWeightsArray(i,3)
  al  = wx*wx + wy*wy + wz*wz
  ro1 = 1./grp%u(i1,1)
  ro2 = 1./grp%u(i2,1)
  vn1 = (wx*grp%u(i1,2) + wy*grp%u(i1,3) + wz*grp%u(i1,4) )*ro1
  vn2 = (wx*grp%u(i2,2) + wy*grp%u(i2,3) + wz*grp%u(i2,4) )*ro2

  vn1 = oneOverMachNumber*sqrt(al) + abs(vn1)
  vn2 = oneOverMachNumber*sqrt(al) + abs(vn2)

  grp%localTurbulentTimeSteps(i1) = grp%localTurbulentTimeSteps(i1) + vn1
  grp%localTurbulentTimeSteps(i2) = grp%localTurbulentTimeSteps(i2) + vn2 

  vis1 = grp%laminarViscosity(i1)
  vis2 = grp%laminarViscosity(i2)

  vis1 = max(vis1*300.0,grp%u(i1,1)*grp%u(i1,6)/grp%u(i1,7))
  vis2 = max(vis2*300.0,grp%u(i2,1)*grp%u(i2,6)/grp%u(i2,7))


  dv1 = 2.0*al*vis1*gammaPr/grp%u(i1,1)
  dv2 = 2.0*al*vis2*gammaPr/grp%u(i2,1)
  grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + dv1
  grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + dv2
 end do

! run over faces

#ifdef PARALLEL
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfBoundaryComSides
 else
  istart = grp%pdp%numberOfBoundaryComSides + 1
  ifinish = grp%brp%numberOfBoundarySides
 end if
 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL

  i1  = grp%brp%sideIndexArray(i,1)
  i2  = grp%brp%sideIndexArray(i,2)
  wx  = grp%brp%sideWeightsArray(i,1)
  wy  = grp%brp%sideWeightsArray(i,2)
  wz  = grp%brp%sideWeightsArray(i,3)
  al  = wx*wx + wy*wy + wz*wz
  ro1 = 1./grp%u(i1,1)
  ro2 = 1./grp%u(i2,1)
  vn1 = (wx*grp%u(i1,2) + wy*grp%u(i1,3) + wz*grp%u(i1,4))*ro1
  vn2 = (wx*grp%u(i2,2) + wy*grp%u(i2,3) + wz*grp%u(i2,4))*ro2

  vn1 = oneOverMachNumber*sqrt(al) + abs(vn1)
  vn2 = oneOverMachNumber*sqrt(al) + abs(vn2)

  vnm = vn1 + vn2

  grp%localTurbulentTimeSteps(i1) = grp%localTurbulentTimeSteps(i1) + vn1 + vnm 
  grp%localTurbulentTimeSteps(i2) = grp%localTurbulentTimeSteps(i2) + vn2 + vnm 

  vis1 = grp%laminarViscosity(i1)
  vis2 = grp%laminarViscosity(i2)

  vis1 = max(vis1*300.0,grp%u(i1,1)*grp%u(i1,6)/grp%u(i1,7))
  vis2 = max(vis2*300.0,grp%u(i2,1)*grp%u(i2,6)/grp%u(i2,7))

  dv1 = 2.0*al*vis1*gammaPr/grp%u(i1,1)
  dv2 = 2.0*al*vis2*gammaPr/grp%u(i2,1)
  dvm = dv1 + dv2
  grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + dv1 + dvm
  grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + dv2 + dvm
 end do
#ifdef PARALLEL

!    write(get_unit(),*)  'Starting first TurbTS comms',kk
 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses
   if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
    call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%localTurbulentTimeSteps(:),grp%pdp%buffer,grp%pdp%sendFlag)
    call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
    call mp_sendv(grp%pdp%processorIDs,i,355,grp%pdp%sendFlag)
   end if
  end do
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,355,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!3250 continue

 do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
!  grp%pdp%sendFlag = 0
!  call mp_check_mesg_arrived(grp%pdp%processorIDs,i,355,grp%pdp%sendFlag,bufferzz,bufferSize)
!  if (grp%pdp%sendFlag.gt.0) then
!   grp%pdp%receivedMessageFromProcess(i) = 1
!   position = 1
    call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                           grp%pdp%realType,grp%localTurbulentTimeSteps(:),grp%pdp%buffer,grp%pdp%sendFlag)
    call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                           grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
!  end if
  end if
 end do

!do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 3250
!end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

!  write(get_unit(),*)  'Starting second TurbTS comms'
 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position = 1
   call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                     grp%pdp%realType,grp%localTurbulentTimeSteps(:),grp%pdp%buffer,grp%pdp%sendFlag)
   call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                     grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
   call mp_sendv(grp%pdp%processorIDs,i,356,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,356,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

!3350 continue

 do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
!  grp%pdp%sendFlag = 0
!  call mp_check_mesg_arrived(grp%pdp%processorIDs,i,356,grp%pdp%sendFlag,bufferzz,bufferSize)
!  if (grp%pdp%sendFlag.gt.0) then
!   grp%pdp%receivedMessageFromProcess(i) = 1
!   position = 1
    call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%localTurbulentTimeSteps(:),grp%pdp%buffer,grp%pdp%sendFlag)
    call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
!  end if
  end if
 end do

!do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0) goto 3350
!end do
#endif PARALLEL


 do i = 1,grp%numberOfNodes
  grp%localTurbulentTimeSteps(i) = 0.5*grp%nodeVolume(i)/(grp%localTurbulentTimeSteps(i)&
   +grp%nodeHelpArray(i,1)/grp%nodeVolume(i))
 end do

 end subroutine makeTurbulentTimeStepsKW
!END-TURB
!-----------------------------------------------------------------------
 subroutine makeTimeSteps(grp,ivd)
 ! make local time steps
 IMPLICIT NONE
include 'mpif.h'

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd

 integer :: i,j,i1,i2
 double precision :: c1,c2,c3,c7,wx,wy,wz,al,ro1,ro2,cc1,cc2,vn1,vn2,t1,t2
 double precision :: vis1,vis2,dv1,dv2,vnm,ccm,dvm,xi1,xi2,fv11,fv12 
 real :: am1,am2,visc1,visc2,minval

 integer :: kk,istart,ifinish,ind
!START-TURB
 real :: tke1,tke2,KE1,KE2
!END-TURB

 ! initialize
 grp%localTimeSteps = 0.0
 grp%nodeHelpArray(:,1) = 0.0

 c1 = ivd%gamma
 c2 = c1 - 1.0
 c3 = 0.
 if(ivd%ReynoldsNumber.gt.1.0e-6) c3 = c2**1.5*ivd%MachNumber**3/ivd%ReynoldsNumber
 c7 = ivd%MachNumber**2*c2*ivd%inflowTemperature

! run over sides 
#ifdef PARALLEL

do kk=1,2
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfComSides
 else
  istart = grp%pdp%numberOfComSides + 1
  ifinish = grp%numberOfSides
 end if


 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL

  i1  = grp%sideIndexArray(i,1) 
  i2  = grp%sideIndexArray(i,2) 
  wx  = grp%sideWeightsArray(i,1)
  wy  = grp%sideWeightsArray(i,2)
  wz  = grp%sideWeightsArray(i,3)
  al  = wx*wx + wy*wy + wz*wz
  ro1 = 1./grp%u(i1,1)
  ro2 = 1./grp%u(i2,1)
  cc1 = max(c1*grp%p(i1)*ro1,0.001)
  cc2 = max(c1*grp%p(i2)*ro2,0.001)
  vn1 = (wx*grp%u(i1,2) + wy*grp%u(i1,3) + wz*grp%u(i1,4) )*ro1 
  vn2 = (wx*grp%u(i2,2) + wy*grp%u(i2,3) + wz*grp%u(i2,4) )*ro2 
  grp%localTimeSteps(i1) = grp%localTimeSteps(i1) + abs(vn1) + sqrt(cc1*al)
  grp%localTimeSteps(i2) = grp%localTimeSteps(i2) + abs(vn2) + sqrt(cc2*al)


  if(ivd%turbulenceModel>0) then 
!START-TURB
   if(ivd%turbulenceModel==1) then 
    visc1 = grp%laminarViscosity(i1) + grp%u(i1,6)/ivd%ReynoldsNumber 
    visc2 = grp%laminarViscosity(i2) + grp%u(i2,6)/ivd%ReynoldsNumber 
   end if
   if(ivd%turbulenceModel==2) then 
    visc1 = grp%laminarViscosity(i1) + grp%u(i1,1)*grp%u(i1,6)/grp%u(i1,7) 
    visc2 = grp%laminarViscosity(i2) + grp%u(i2,1)*grp%u(i2,6)/grp%u(i2,7) 
   end if
   if(ivd%turbulenceModel==3) then
    tke1 = grp%u(i1,6)
    tke2 = grp%u(i2,6)
    KE1 = (0.5*(grp%u(i1,2)**2.+grp%u(i1,3)**2.+grp%u(i1,4)**2.))
    KE2 = (0.5*(grp%u(i2,2)**2.+grp%u(i2,3)**2.+grp%u(i2,4)**2.))
    if(grp%u(i1,6).eq.0.0 .and. KE1 .gt. 0.0)tke1=KE1
    if(grp%u(i2,6).eq.0.0 .and. KE2 .gt. 0.0)tke2=KE2
    visc1 = grp%laminarViscosity(i1) + grp%u(i1,1)*tke1/grp%u(i1,7) 
    visc2 = grp%laminarViscosity(i2) + grp%u(i2,1)*tke2/grp%u(i2,7) 
   end if
!END-TURB
  else
   visc1 = grp%laminarViscosity(i1) 
   visc2 = grp%laminarViscosity(i2) 
  end if
  if(ivd%ReynoldsNumber>0.0) then 
   am1 = 1.0/ivd%ReynoldsNumber
   am2 = 1.0/ivd%ReynoldsNumber
   vis1= amax1(am1,visc1)
   vis2= amax1(am2,visc2)
   dv1 = 2.0*vis1*al/grp%u(i1,1)
   dv2 = 2.0*vis2*al/grp%u(i2,1)
   grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + dv2
   grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + dv1
  end if
 end do

! run over faces 

#ifdef PARALLEL
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfBoundaryComSides
 else
  istart = grp%pdp%numberOfBoundaryComSides + 1
  ifinish = grp%brp%numberOfBoundarySides
 end if
 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL

  i1  = grp%brp%sideIndexArray(i,1)
  i2  = grp%brp%sideIndexArray(i,2)
  wx  = grp%brp%sideWeightsArray(i,1)
  wy  = grp%brp%sideWeightsArray(i,2)
  wz  = grp%brp%sideWeightsArray(i,3)
  al  = wx*wx + wy*wy + wz*wz
  ro1 = 1./grp%u(i1,1)
  ro2 = 1./grp%u(i2,1)
  cc1 = max(c1*grp%p(i1)*ro1,0.001)
  cc2 = max(c1*grp%p(i2)*ro2,0.001)
  vn1 = (wx*grp%u(i1,2) + wy*grp%u(i1,3) + wz*grp%u(i1,4) )*ro1 
  vn2 = (wx*grp%u(i2,2) + wy*grp%u(i2,3) + wz*grp%u(i2,4) )*ro2 
  vn1 = abs(vn1)
  vn2 = abs(vn2)
  vnm = vn1 + vn2 
  cc1 = sqrt(cc1*al)
  cc2 = sqrt(cc2*al)
  ccm = cc1 + cc2 


  grp%localTimeSteps(i1) = grp%localTimeSteps(i1) + vn1 + vnm + cc1 + ccm 
  grp%localTimeSteps(i2) = grp%localTimeSteps(i2) + vn2 + vnm + cc2 + ccm 

  if(ivd%turbulenceModel>0) then 
!START-TURB
   if(ivd%turbulenceModel==1) then 
    visc1 = grp%laminarViscosity(i1) + grp%u(i1,6)/ivd%ReynoldsNumber 
    visc2 = grp%laminarViscosity(i2) + grp%u(i2,6)/ivd%ReynoldsNumber 
   end if
   if(ivd%turbulenceModel==2) then 
    visc1 = grp%laminarViscosity(i1) + grp%u(i1,1)*grp%u(i1,6)/grp%u(i1,7) 
    visc2 = grp%laminarViscosity(i2) + grp%u(i2,1)*grp%u(i2,6)/grp%u(i2,7) 
   end if
   if(ivd%turbulenceModel==3) then
    tke1 = grp%u(i1,6)
    tke2 = grp%u(i2,6)
    KE1 = (0.5*(grp%u(i1,2)**2.+grp%u(i1,3)**2.+grp%u(i1,4)**2.))
    KE2 = (0.5*(grp%u(i2,2)**2.+grp%u(i2,3)**2.+grp%u(i2,4)**2.))
    if(grp%u(i1,6).eq.0.0 .and. KE1 .gt. 0.0)tke1=KE1
    if(grp%u(i2,6).eq.0.0 .and. KE2 .gt. 0.0)tke2=KE2
    visc1 = grp%laminarViscosity(i1) + grp%u(i1,1)*tke1/grp%u(i1,7)
    visc2 = grp%laminarViscosity(i2) + grp%u(i2,1)*tke2/grp%u(i2,7) 
   end if
!END-TURB
  else
   visc1 = grp%laminarViscosity(i1) 
   visc2 = grp%laminarViscosity(i2) 
  end if

  if(ivd%ReynoldsNumber>0) then 
   am1 = 1.0/ivd%ReynoldsNumber
   am2 = 1.0/ivd%ReynoldsNumber

   vis1= amax1(am1,visc1)
   vis2= amax1(am2,visc2)

   dv1 = 2.0*vis1*al/grp%u(i1,1)
   dv2 = 2.0*vis2*al/grp%u(i2,1)
   dvm = dv1 + dv2
   grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) + dv1 + dvm
   grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + dv2 + dvm
  end if
 end do
#ifdef PARALLEL

 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses
   if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
    call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%localTimeSteps(:),grp%pdp%buffer,grp%pdp%sendFlag)
    call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
    call mp_sendv(grp%pdp%processorIDs,i,55,grp%pdp%sendFlag)
   end if
  end do
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,55,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
    call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                           grp%pdp%realType,grp%localTimeSteps(:),grp%pdp%buffer,grp%pdp%sendFlag)
    call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                           grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
  end if
 end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

!  write(get_unit(),*)  'Starting second TS comms'
 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position = 1
   call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                     grp%pdp%realType,grp%localTimeSteps(:),grp%pdp%buffer,grp%pdp%sendFlag)
   call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                     grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
   call mp_sendv(grp%pdp%processorIDs,i,56,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,56,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
    call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%localTimeSteps(:),grp%pdp%buffer,grp%pdp%sendFlag)
    call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
  end if
 end do

#endif PARALLEL


 do i = 1,grp%numberOfNodes
  grp%localTimeSteps(i) = 0.5/&
    (grp%localTimeSteps(i)/grp%nodeVolume(i)&
     +grp%nodeHelpArray(i,1)/(grp%nodeVolume(i)**2))
 end do


 minVal = 9999999.9

 do i=1,grp%numberOfNodes
  if(grp%localTimeSteps(i)<minVal) minVal = grp%localTimeSteps(i)
 end do
 
 
 if (ivd%explicit) then
         ivd%physicalTimeStep=minVal
!        write(*,*) "Physical Timestep: ", minVal
 end if

 end subroutine makeTimeSteps
!-----------------------------------------------------------------------
 subroutine makePressureField(grp,ivd,u)
 ! makes pressure field from the ideal gas equation
 IMPLICIT NONE

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd
 real,pointer :: u(:,:)
 integer :: i
 real :: gammaMinusOne,twoThirds


 gammaMinusOne = ivd%gamma-1.0

 do i=1,grp%numberOfNodes
  grp%p(i) = gammaMinusOne*(u(i,5)&
                 -0.5*(u(i,2)**2+u(i,3)**2+u(i,4)**2)/u(i,1)) 
  grp%p(i) = max(grp%p(i),0.1) 
 end do
 end subroutine makePressureField
!-------------------------------------------------------------------------
 subroutine smoothResidual2(grp,ivd,rhs,sfactor)  
 ! smooths residual to increase allowable CFL number

  IMPLICIT NONE
include 'mpif.h'

  type(GridSolverData) :: grp
  type(InputVariablesData) :: ivd
  real :: rhs(:,:)
  real :: sfactor,totVol

  integer :: i,i1,i2,j,k
  real :: dw

  integer :: kk,istart,ifinish

  grp%nodeHelpArray(1:grp%numberOfNodes,1) = 0.0 ! initialize 

! make help variables
 do k=1,5 ! number of PDE's 

#ifdef PARALLEL

  do kk=1,2
   if(kk==1) then
    istart = 1
    ifinish = grp%pdp%numberOfComSides
   else
    istart = grp%pdp%numberOfComSides + 1
    ifinish = grp%numberOfSides
   end if


   do i=istart,ifinish

#else PARALLEL

   do i=1,grp%numberOfSides

#endif PARALLEL

    i1 = grp%sideIndexArray(j,1)
    i2 = grp%sideIndexArray(j,2)
    dw = rhs(i1,k) - rhs(i2,k)
    totVol = grp%nodeVolume(i1)+grp%nodeVolume(i2)
    
    grp%nodeHelpArray(i1,1) = grp%nodeHelpArray(i1,1) - grp%nodeVolume(i2)*dw/totVol
    grp%nodeHelpArray(i2,1) = grp%nodeHelpArray(i2,1) + grp%nodeVolume(i1)*dw/totVol
   end do

#ifdef PARALLEL

   if (kk.eq.1) then
    do i=1,grp%pdp%numberOfProcesses
     if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
      call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
      position = 1
      call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                         grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
      call mp_sendv(grp%pdp%processorIDs,i,57,grp%pdp%sendFlag)
     end if
    end do
   end if
  end do

   call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,57,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

   do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                             grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
    end if
   end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

   do i=1,grp%pdp%numberOfProcesses
    if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
     call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
     position = 1
     call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                       grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
     call mp_sendv(grp%pdp%processorIDs,i,58,grp%pdp%sendFlag)
    end if
   end do

   call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,58,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

   do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
       call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                         grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
    end if
   end do

#endif PARALLEL


   do j=1,grp%numberOfNodes
    rhs(j,i) = rhs(j,i)+sfactor*grp%nodeHelpArray(j,i)
   end do
  end do


 end subroutine smoothResidual2
!-----------------------------------------------------------------------
 subroutine makeIORegister(grp,ia)
 ! for internal outflow boundaries
 IMPLICIT NONE
 
 type(GridSolverData) :: grp

 integer :: i,allocateStatus,ist,ien,ind1,ind2,numberOfIONodes,ia
 real :: directionVector(3),innerProduct,len
 
 ist = grp%brp%nodeIndicatorRegister(-4)
 ien = grp%brp%nodeIndicatorRegister(-3)
 numberOfIONodes = ien-ist
 if(numberOfIONodes>0) numberOfIONodes = numberOfIONodes + 1
 grp%brp%numberOfIONodes = numberOfIONodes
 write(*,*) "Number of IO nodes: ",numberOfIONodes

 grp%nodeHelpArray(:,1) = -99999.9
 grp%nodeHelpArray(:,2) = 0.0
 directionVector(1) = -1.0 
 directionVector(2) = 0.0 
 directionVector(3) = 0.0 
 
 if(numberOfIONodes>0) then 
  if(ia.eq.0) then
  allocate(grp%brp%IORegister(numberOfIONodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: makeIORegister out of memory"
  end if
   
  do i=1,grp%numberOfSides
   ind1 = grp%sideIndexArray(i,1)
   ind2 = grp%sideIndexArray(i,2)
  
   len = sqrt(sum(grp%sideWeightsArray(i,:)*grp%sideWeightsArray(i,:))) 
   if(len>0) then 
    innerProduct = abs(sum(grp%sideWeightsArray(i,:)*directionVector)/len)
    if(innerProduct>grp%nodeHelpArray(ind1,1)) then 
     grp%nodeHelpArray(ind1,1) = innerProduct   
     grp%nodeHelpArray(ind1,2) = ind2
    end if 
    if(innerProduct>grp%nodeHelpArray(ind2,1)) then 
     grp%nodeHelpArray(ind2,1) = innerProduct   
     grp%nodeHelpArray(ind2,2) = ind1   
    end if 
   end if
  end do 

  do i=ist,ien
   if(grp%nodeHelpArray(grp%brp%nodeIndicatorArray(i),2)>0.1) then 
    grp%brp%IORegister(i-ist+1) = grp%nodeHelpArray(grp%brp%nodeIndicatorArray(i),2)
   else
    grp%brp%IORegister(i-ist+1) = grp%brp%nodeIndicatorRegister(i)
   end if
  end do
 else
  nullify(grp%brp%IORegister)
 end if
 end subroutine makeIORegister
!-----------------------------------------------------------------------
 subroutine readGridComputationData(finegrp,hgrp,grp,ivd,INFILE,ia)
! reads data from file, created by preprocessor 
  IMPLICIT NONE
include 'mpif.h'

  type(GridSolverData) :: finegrp,hgrp,grp
  type(InputVariablesData) :: ivd
  integer :: INFILE

  integer :: i,ind1,ind2,j,k,allocateStatus,ip,coarseNodeNumber,ist,ien,ib,dummy,blockSize
  real :: dx,dy,nx,ny,nz,ww,areaIncrement,dummyr
  real :: x(2),x0(3),x1(3),x2(3),x3(3),det1,det2,det3,det4
  real :: maxError,sideError
  real :: xm(3),xp(3)
  integer :: sideNumber,blockNumber,baseNumber
  real :: wx(3),ftmp,fact,dxv(3)
  integer :: inletNumber,ia
  real :: prevGCLWeight
  real :: coordinatesP(3),coordinatesP2(3)
  real :: mcoorm
  integer :: mcoorind

  integer :: kk,istart,ifinish 

   write(*,*) "reading boundary..."
  if(grp%gridNumber==1) then ! finest mesh

   call readBoundarySolverData(ivd,grp%brp,INFILE,.true.,ia)

   ! read sides
   write(*,*) "reading sides..."
   read(INFILE) grp%numberOfSides
   write(*,*) "Number of sides: ",grp%numberOfSides
   if(ia.eq.0) then
   allocate(grp%sideIndexArray(grp%numberOfSides,2),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   allocate(grp%sideWeightsArray(grp%numberOfSides,3),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   allocate(grp%sideLengthArray(grp%numberOfSides),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   allocate(grp%GCLWeightsArray(grp%numberOfSides),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   end if

   do i=1,grp%numberOfSides
    read(INFILE) grp%sideIndexArray(i,1:2),grp%sideWeightsArray(i,1:3),grp%GCLWeightsArray(i),&
                      prevGCLWeight,grp%sideLengthArray(i)
    if(ivd%timeScheme==2) then
     grp%GCLWeightsArray(i) = 1.5*grp%GCLWeightsArray(i)-0.5*prevGCLWeight
    end if
   end do

   nullify(grp%coefficientSideLengthArray)

  ! read node volume and other stuff
   read(INFILE) grp%numberOfNodes
   write(*,*) "Number of nodes: ",grp%numberOfNodes
   if(ia.eq.0) then
   allocate(grp%nodeVolume(grp%numberOfNodes),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   allocate(grp%nodeVolumeP(grp%numberOfNodes),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   allocate(grp%nodeVolumeP2(grp%numberOfNodes),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   allocate(grp%wallDistance(grp%numberOfNodes),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   allocate(grp%wallDistanceBoundaryNodeArray(grp%numberOfNodes),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   allocate(grp%coordinates(grp%numberOfNodes,3),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   allocate(grp%coordinateMovement(grp%numberOfNodes,3),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   end if
!START-GUST
   allocate(grp%IPER(grp%numberOfNodes),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory to allocate gust pertubation indexes"
!END-GUST

   do i=1,grp%numberOfNodes
    read(INFILE) dummy,grp%nodeVolume(dummy),grp%nodeVolumeP(dummy),grp%nodeVolumeP2(dummy),&
                 grp%wallDistance(dummy),&
                 grp%wallDistanceBoundaryNodeArray(dummy),grp%coordinates(dummy,1:3),&
                 coordinatesP(1:3),coordinatesP2(1:3)

   if(ivd%timeScheme==1) then
     grp%coordinateMovement(dummy,1:3) = grp%coordinates(dummy,1:3)-coordinatesP(1:3)
    else
     grp%coordinateMovement(dummy,1:3) = 1.5*grp%coordinates(dummy,1:3)-2.0*coordinatesP(1:3)&
                                        +0.5*coordinatesP2(1:3)
    end if
    if(sum(grp%coordinateMovement(dummy,1:3)*grp%coordinateMovement(dummy,1:3))>mcoorm) then
      mcoorm = sum(grp%coordinateMovement(dummy,1:3)*grp%coordinateMovement(dummy,1:3))
      mcoorind = i
    end if
   end do
   write(*,*) "MAX MOVE: ",mcoorm,mcoorind

   if(ia.eq.0) then
   allocate(grp%brp%sideLengthArray(grp%brp%numberOfBoundarySides),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   end if

   do i=1,grp%brp%numberOfBoundarySides
    xm = grp%coordinates(grp%brp%sideIndexArray(i,1),:)-grp%coordinates(grp%brp%sideIndexArray(i,2),:)
    grp%brp%sideLengthArray(i) = sqrt(sum(xm*xm))
   end do


   ! read wall trip stuff

   write(*,*) "reading wall trip data..."

   read(INFILE) grp%numberOfTripNodes
   write(*,*) 'number of trip nodes=', grp%numberOfTripNodes
   if(ia.eq.0) then
   if(grp%numberOfTripNodes>0) then 
    allocate(grp%tripLineIndexes(grp%numberOfTripNodes),stat=allocateStatus)
    if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
    allocate(grp%tripWallLength(grp%numberOfTripNodes),stat=allocateStatus)
    if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   end if
   end if
   do i=1,grp%numberOfTripNodes
    read(INFILE) grp%tripLineIndexes(i),grp%tripWallLength(i)
   end do

   read(INFILE) grp%numberOfTripFieldNodes
   if(ia.eq.0) then
   if(grp%numberOfTripNodes>0.and.grp%numberOfTripFieldNodes>0) then 
    allocate(grp%tripNodeFieldIndexes(grp%numberOfTripFieldNodes,2),stat=allocateStatus)
    if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
    allocate(grp%tripNodeFieldDistances(grp%numberOfTripFieldNodes),stat=allocateStatus)
    if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
    grp%tripNodeFieldIndexes = 0
    grp%tripNodeFieldDistances = 0.0
   end if
   end if
   do i=1,grp%numberOfTripFieldNodes
    read(INFILE) grp%tripNodeFieldIndexes(i,1:2),grp%tripNodeFieldDistances(i)
   end do 


  if(ia.eq.0) nullify(grp%seedNodeRegister)

   write(*,*) "reading parallelization data..."

   call readGridParallelizationData(INFILE,grp%pdp,ia)

   call makeBoundaryNormals(grp)

   if(ia.eq.0) nullify(grp%mapComBuffer)
   
 
  else ! coarser meshes
   call readBoundarySolverData(ivd,grp%brp,INFILE,.false.,ia)

   ! read sides
   read(INFILE) grp%numberOfSides
   if(grp%pdp%currentDomain==1) write(*,*) "number of sides: ",grp%numberOfSides
   if(ia.eq.0) then
   allocate(grp%sideIndexArray(grp%numberOfSides,2),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   allocate(grp%sideWeightsArray(grp%numberOfSides,3),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   allocate(grp%GCLWeightsArray(grp%numberOfSides),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   allocate(grp%sideLengthArray(grp%numberOfSides),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   allocate(grp%coefficientSideLengthArray(grp%numberOfSides),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   end if
   do i=1,grp%numberOfSides
    read(INFILE) grp%sideIndexArray(i,1:2),grp%sideWeightsArray(i,1:3),grp%GCLWeightsArray(i),dummyr,grp%sideLengthArray(i)
   end do

   do i=1,grp%numberOfSides
    grp%coefficientSideLengthArray(i) = sqrt(sum(grp%sideWeightsArray(i,:)*grp%sideWeightsArray(i,:)))
   end do
   ! read node volume and some other stuff
   read(INFILE) grp%numberOfNodes
   if(ia.eq.0) then
   allocate(grp%nodeVolume(grp%numberOfNodes),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   allocate(grp%wallDistance(grp%numberOfNodes),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   allocate(grp%wallDistanceBoundaryNodeArray(grp%numberOfNodes),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   allocate(grp%seedNodeRegister(grp%numberOfNodes),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   end if
   grp%coordinates=>finegrp%coordinates

   do i=1,grp%numberOfNodes
    read(INFILE) dummy,grp%nodeVolume(dummy),grp%wallDistance(dummy),&
                 grp%wallDistanceBoundaryNodeArray(dummy),grp%seedNodeRegister(i)
   end do


   if(ia.eq.0) then
   nullify(grp%coordinates)
   nullify(grp%nodeVolumeP)
   nullify(grp%nodeVolumeP2)
   nullify(grp%tripNodeFieldIndexes)
   nullify(grp%tripNodeFieldDistances)
   nullify(grp%tripLineIndexes)
   nullify(grp%tripWallLength) 
   end if

   grp%numberOfTripNodes = 0
   grp%numberOfTripFieldNodes = 0 

   call readGridParallelizationData(INFILE,grp%pdp,ia)

   ! read prolongation operator (coarse to fine)

   grp%pod%numberOfFineNodes = finegrp%numberOfNodes 
   grp%pod%numberOfCoarseNodes = grp%numberOfNodes 
   call readProlongationOperatorData(grp%pod,INFILE,ia,grp%pdp%currentDomain)


   ! read restriction operator (fine to coarse)

   grp%rod%numberOfFineNodes = finegrp%numberOfNodes
   grp%rod%numberOfCoarseNodes = grp%numberOfNodes  

   call readRestrictionOperatorData(grp%rod,INFILE,ia,grp%pdp%currentDomain)
   call readComParallelizationData(INFILE,grp%pdp,ia)

   if(ia.eq.0) then
 !START-TURB
  if(ivd%turbulenceModel==0 .or. ivd%turbulenceModel==1) then
   allocate(grp%mapComBuffer(grp%rod%numberOfComNodes,6),stat=allocateStatus)
   if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate unknown field"
  end if
  if(ivd%turbulenceModel==2) then
    allocate(grp%mapComBuffer(grp%rod%numberOfComNodes,7),stat=allocateStatus)
   if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate unknown field"
  end if
  if(ivd%turbulenceModel==3) then
    allocate(grp%mapComBuffer(grp%rod%numberOfComNodes,7),stat=allocateStatus)
   if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate unknown field"
  end if
!END-TURB
   end if

   call makeBoundaryNormals(grp)

   if(ia.eq.0) then
   allocate(grp%coordinateMovement(grp%numberOfNodes,3),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   end if

   grp%coordinateMovement = 0.0
!  call mapFineToCoarse(hgrp%coordinateMovement,grp,grp%coordinateMovement)

   if(ia.eq.0) then
   allocate(grp%coordinates(grp%numberOfNodes,3),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory to allocate coordinate field"
   end if

  end if

! do some allocation

  ! unknowns at this and previous time levels

  if(ia.eq.0) then
!START-TURB
 if(ivd%turbulenceModel==0 .or. ivd%turbulenceModel==1) then
  allocate(grp%u(grp%numberOfNodes,6),stat=allocateStatus)
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate unknown field"
 
  allocate(grp%uprev(grp%numberOfNodes,6),stat=allocateStatus)
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate unknown field"

  allocate(grp%uold(grp%numberOfNodes,6),stat=allocateStatus)
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate unknown field"

  allocate(grp%uold2(grp%numberOfNodes,6),stat=allocateStatus)
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate unknown field"

  allocate(grp%rhs(grp%numberOfNodes,6),stat=allocateStatus)
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate unknown field"
 end if
 if(ivd%turbulenceModel==2 .or. ivd%turbulenceModel==3) then
  allocate(grp%u(grp%numberOfNodes,7),stat=allocateStatus)
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate unknown field"
 
  allocate(grp%uprev(grp%numberOfNodes,7),stat=allocateStatus)
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate unknown field"

  allocate(grp%uold(grp%numberOfNodes,7),stat=allocateStatus)
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate unknown field"

  allocate(grp%uold2(grp%numberOfNodes,7),stat=allocateStatus)
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate unknown field"

  allocate(grp%rhs(grp%numberOfNodes,7),stat=allocateStatus)
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate unknown field"
 end if

  allocate(grp%dissipation(grp%numberOfNodes,5),stat=allocateStatus)
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate dissipation field"
!END-TURB

 ! pressure field
  allocate(grp%p(grp%numberOfNodes),stat=allocateStatus)
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate pressure field"

  allocate(grp%localTimeSteps(grp%numberOfNodes),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
  allocate(grp%localTurbulentTimeSteps(grp%numberOfNodes),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"

  allocate(grp%laminarViscosity(grp%numberOfNodes),stat=allocateStatus)
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate viscosity field"
  grp%laminarViscosity = 0.0
  allocate(grp%vorticity(grp%numberOfNodes),stat=allocateStatus)
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory in readGridData"
  grp%vorticity = 0.0

  allocate(grp%divergence(grp%numberOfNodes),stat=allocateStatus)
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory in readGridData"
  grp%divergence = 0.0

!START-TURB
  allocate(grp%msrt(grp%numberOfNodes),stat=allocateStatus)!mean-strain-rate tensor
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory in readGridData"
  grp%msrt = 0.0
  allocate(grp%imsr(grp%numberOfNodes),stat=allocateStatus)!invariant measure of the strain rate
  if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory in readGridData"
  grp%imsr = 0.0
 if(ivd%turbulenceModel==1) then 
  allocate(grp%turbulenceDiffusionTerm(grp%numberOfNodes),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
 end if
 if(ivd%turbulenceModel==2 .or. ivd%turbulenceModel==3) then !MOD OPT
  allocate(grp%turbulenceDiffusionTermK(grp%numberOfNodes),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
  allocate(grp%turbulenceDiffusionTermW(grp%numberOfNodes),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
  allocate(grp%SK(grp%numberOfNodes),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
  grp%SK(:) = 0.0
  allocate(grp%SW(grp%numberOfNodes),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
  grp%SW(:) = 0.0
 end if
 if(ivd%turbulenceModel==3) then !SST
  allocate(grp%dkdw(grp%numberOfNodes),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
  allocate(grp%F1(grp%numberOfNodes),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
  allocate(grp%F2(grp%numberOfNodes),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
  allocate(grp%TVR(grp%numberOfNodes),stat=allocateStatus)!turbulent viscosity ratio
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
   if (ivd%SAS) then
    allocate(grp%D1U(9,grp%numberOfNodes),stat=allocateStatus)!first derivative of velocity
    if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
    allocate(grp%D2U(3,grp%numberOfNodes),stat=allocateStatus)!second derivative of velocity
    if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
! Allocation of auxiliary array for printing out (not necessary for computation)
    allocate(grp%sas(4,grp%numberOfNodes),stat=allocateStatus)!second derivative of velocity
    if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
    grp%sas(:,:) = 0.0
   end if
 end if 
!END-TURB
 
  if(grp%gridNumber==1) then 
   allocate(grp%wallStress(grp%brp%numberOfBoundaryNodes,4),stat=allocateStatus)
   if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
  else
   nullify(grp%wallStress) 
  end if

  ! if initial grid, allocate help arrays, else point to original 
  if(grp%gridNumber==1) then 
   nullify(grp%sourceTerm)
   !nodeHelpArray: x25 for HLLC
   allocate(grp%nodeHelpArray(grp%numberOfNodes,25),stat=allocateStatus)
   if(allocateStatus/=0)&
      STOP "ERROR: Not enough memory to allocate help field"
  else
   ! allocate FAS source term
!START-TURB
 if(ivd%turbulenceModel==0 .or. ivd%turbulenceModel==1) then
   allocate(grp%sourceTerm(grp%numberOfNodes,6),stat=allocateStatus)
   if(allocateStatus/=0)&
      STOP "ERROR: Not enough memory to allocate source term field"
 end if
 if(ivd%turbulenceModel==2 .or. ivd%turbulenceModel==3) then
   allocate(grp%sourceTerm(grp%numberOfNodes,7),stat=allocateStatus)
   if(allocateStatus/=0)&
      STOP "ERROR: Not enough memory to allocate source term field"
 end if
!END-TURB
   ! to avoid using memory to create help variables for coarse grids
   grp%nodeHelpArray => finegrp%nodeHelpArray
   nullify(grp%tripLineIndexes)
   nullify(grp%tripNodeFieldDistances)
   nullify(grp%tripNodeFieldIndexes)
   nullify(grp%tripwallLength)
  end if

! set Runge-Kutta timesteps
  
  allocate(grp%RKCoefficients(3),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"
  grp%RKCoefficients(1) = 0.6
  grp%RKCoefficients(2) = 0.6 
  grp%RKCoefficients(3) = 1.0 
  


  allocate(grp%wallArea(grp%numberOfNodes),stat=allocateStatus)
   if(allocateStatus/=0)&
     STOP "ERROR: Not enough memory to allocate in readGridData"

 allocate(grp%nodeConnectivityArray(grp%numberOfNodes),stat=allocateStatus)
 if(allocateStatus/=0) STOP "ERROR: Not enough memory in readComputationData"

 end if


  grp%residualScalingFactor = 1.0/float(ivd%numberOfProcesses*grp%numberOfNodes)
  grp%nodeConnectivityArray = 0

  call makeSurfaceArea(grp)

#ifdef PARALLEL

 do kk=1,2
  if(kk==1) then
   istart = 1
   ifinish = grp%pdp%numberOfComSides
  else
   istart = grp%pdp%numberOfComSides + 1
   ifinish = grp%numberOfSides
  end if

  do i=istart,ifinish

#else PARALLEL

  do i=1,grp%numberOfSides

#endif PARALLEL
   ind1 = grp%sideIndexArray(i,1)
   ind2 = grp%sideIndexArray(i,2)
   grp%nodeConnectivityArray(ind1) = grp%nodeConnectivityArray(ind1) + 1.0
   grp%nodeConnectivityArray(ind2) = grp%nodeConnectivityArray(ind2) + 1.0
  end do

#ifdef PARALLEL

! send domain boundary residuals

 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses
   if (i.ne.grp%pdp%currentDomain.and.grp%pdp%numberOfSendNodes(i)>0) then
    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
    call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%nodeConnectivityArray,grp%pdp%buffer,grp%pdp%sendFlag)
    call mp_sendv(grp%pdp%processorIDs,i,139,grp%pdp%sendFlag)
   end if
  end do
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,139,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
    call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                           grp%pdp%realType,grp%nodeConnectivityArray,grp%pdp%buffer,grp%pdp%sendFlag)
  end if
 end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position = 1
   call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                     grp%pdp%realType,grp%nodeConnectivityArray,grp%pdp%buffer,grp%pdp%sendFlag)
   call mp_sendv(grp%pdp%processorIDs,i,129,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,129,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )


 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
    call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%nodeConnectivityArray,grp%pdp%buffer,grp%pdp%sendFlag)
  end if
 end do

#endif PARALLEL

! set engine inlet parameters

 grp%brp%engineInletAreas = 0.0
 grp%brp%engineInletNormals = 0.0
 if(grp%gridNumber==1) then
  do i=1,grp%brp%numberOfEngineInletSides
   ind1 = grp%brp%engineInletSideIndexes(i,1)
   ind2 = grp%brp%engineInletSideIndexes(i,2)
   inletNumber = grp%brp%engineInletSideIndexes(i,3)
   wx = grp%brp%engineInletSideCoefficients(i,:)
   grp%brp%engineInletAreas(inletNumber) = grp%brp%engineInletAreas(inletNumber) + 2.0*sqrt(sum(wx*wx)) 
   grp%brp%engineInletNormals(inletNumber,:) = grp%brp%engineInletNormals(inletNumber,:) + wx 
  end do
 else
  do i=1,grp%brp%numberOfEngineInletSides
   ind1 = grp%brp%engineInletSideIndexes(i,1)
   ind2 = grp%brp%engineInletSideIndexes(i,2)
   inletNumber = grp%brp%engineInletSideIndexes(i,3)
   wx = grp%brp%engineInletSideCoefficients(i,:)
   grp%brp%engineInletAreas(inletNumber) = grp%brp%engineInletAreas(inletNumber) + 2.0*sqrt(sum(wx*wx))
   grp%brp%engineInletNormals(inletNumber,:) = grp%brp%engineInletNormals(inletNumber,:) + wx 
  end do
 end if

#ifdef PARALLEL

  call mp_init_buffer(grp%pdp%processorIDs,-1,grp%pdp%sendFlag)
  position = 1
  call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,grp%brp%engineInletAreas(1),grp%pdp%sendFlag)
  call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,grp%brp%engineInletAreas(2),grp%pdp%sendFlag)
  call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,grp%brp%engineInletNormals(1,1),grp%pdp%sendFlag)
  call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,grp%brp%engineInletNormals(1,2),grp%pdp%sendFlag)
  call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,grp%brp%engineInletNormals(1,3),grp%pdp%sendFlag)
  call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,grp%brp%engineInletNormals(2,1),grp%pdp%sendFlag)
  call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,grp%brp%engineInletNormals(2,2),grp%pdp%sendFlag)
  call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,grp%brp%engineInletNormals(2,3),grp%pdp%sendFlag)

  call mp_snd_others(grp%pdp%numberOfProcesses,grp%pdp%currentDomain,grp%pdp%processorIDs,347,grp%pdp%sendFlag)

  call init_recv_others(grp%pdp%numberOfProcesses,grp%pdp%currentDomain,grp%pdp%receivedMessageFromProcess)

  call mp_recv_others(grp%pdp%numberOfProcesses,grp%pdp%currentDomain,grp%pdp%processorIDs,347,grp%pdp%sendFlag)

  call mp_wait_comms( grp%pdp%sendFlag )

  do i=1,grp%pdp%numberOfProcesses
    if( i.ne.grp%pdp%currentDomain ) then
     call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
     grp%brp%engineInletAreas(1) = grp%brp%engineInletAreas(1) + ftmp
     call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
     grp%brp%engineInletAreas(2) = grp%brp%engineInletAreas(2) + ftmp
     call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
     grp%brp%engineInletNormals(1,1) = grp%brp%engineInletNormals(1,1) + ftmp 
     call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
     grp%brp%engineInletNormals(1,2) = grp%brp%engineInletNormals(1,2) + ftmp 
     call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
     grp%brp%engineInletNormals(1,3) = grp%brp%engineInletNormals(1,3) + ftmp 
     call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
     grp%brp%engineInletNormals(2,1) = grp%brp%engineInletNormals(2,1) + ftmp 
     call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
     grp%brp%engineInletNormals(2,2) = grp%brp%engineInletNormals(2,2) + ftmp 
     call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
     grp%brp%engineInletNormals(2,3) = grp%brp%engineInletNormals(2,3) + ftmp 
   end if
  end do

#endif PARALLEL

  fact = sqrt(sum(grp%brp%engineInletNormals(1,:)*grp%brp%engineInletNormals(1,:)))
  if(fact>1.0e-12) then 
   grp%brp%engineInletNormals(1,:) = grp%brp%engineInletNormals(1,:)/fact
  end if
  fact = sqrt(sum(grp%brp%engineInletNormals(2,:)*grp%brp%engineInletNormals(2,:)))
  if(fact>1.0e-12) then 
   grp%brp%engineInletNormals(2,:) =  grp%brp%engineInletNormals(2,:)/fact
  end if

! Check that coefficients add up

  grp%uprev(:,1:4) = 0.0
#ifdef PARALLEL

do kk=1,2
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfComSides
 else
  istart = grp%pdp%numberOfComSides + 1
  ifinish = grp%numberOfSides
 end if


 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL

   ind1 = grp%sideIndexArray(i,1)
   ind2 = grp%sideIndexArray(i,2)
   grp%uprev(ind1,1:3) = grp%uprev(ind1,1:3) + grp%sideWeightsArray(i,:)
   grp%uprev(ind2,1:3) = grp%uprev(ind2,1:3) - grp%sideWeightsArray(i,:)
   grp%uprev(ind1,4) = max(sqrt(sum(grp%sideWeightsArray(i,:)*grp%sideWeightsArray(i,:))),grp%uprev(ind1,4))
   grp%uprev(ind2,4) = max(sqrt(sum(grp%sideWeightsArray(i,:)*grp%sideWeightsArray(i,:))),grp%uprev(ind2,4))
  end do

#ifdef PARALLEL
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfBoundaryComSides
 else
  istart = grp%pdp%numberOfBoundaryComSides + 1
  ifinish = grp%brp%numberOfBoundarySides
 end if
 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL

   ind1 = grp%brp%sideIndexArray(i,1)
   ind2 = grp%brp%sideIndexArray(i,2)
   grp%uprev(ind1,1:3) = grp%uprev(ind1,1:3) + 2.0*grp%brp%sideWeightsArray(i,:)
   grp%uprev(ind2,1:3) = grp%uprev(ind2,1:3) + 2.0*grp%brp%sideWeightsArray(i,:)
   grp%uprev(ind1,4) = max(sqrt(sum(grp%brp%sideWeightsArray(i,:)*grp%brp%sideWeightsArray(i,:))),grp%uprev(ind1,4))
   grp%uprev(ind2,4) = max(sqrt(sum(grp%brp%sideWeightsArray(i,:)*grp%brp%sideWeightsArray(i,:))),grp%uprev(ind2,4))
  end do

#ifdef PARALLEL

 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses
   if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then

    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
    do j=1,4
     call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%uprev(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
    call mp_sendv(grp%pdp%processorIDs,i,61,grp%pdp%sendFlag)
   end if
  end do
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,61,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
    do j=1,4
     call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                            grp%pdp%realType,grp%uprev(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
  end if
 end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position = 1
   do j=1,4
    call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                      grp%pdp%realType,grp%uprev(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
   end do
   call mp_sendv(grp%pdp%processorIDs,i,62,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,62,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
    do j = 1,4
     call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%uprev(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
    end do
  end if
 end do

#endif PARALLEL

  maxError = 0.0

  do i=1,grp%numberOfNodes
   if(grp%uprev(i,4)>1.0e-15) then
    sideError = sqrt(sum(grp%uprev(i,1:3)*grp%uprev(i,1:3)))/grp%uprev(i,4)
    if(sideError>1.0e-6) then
     write(*,*) "Side error: ",i,sideError,grp%uprev(i,1:3),grp%uprev(i,4)
    end if
    if(sideError>maxError) maxError = sideError
   end if
  end do
 if(grp%pdp%currentDomain==1) write(*,*) "max side error: ",maxError

! check conservation

 write(*,*) "checking conservation..."

 if(grp%gridNumber==1) then
 grp%uprev(:,1) = 0.0

#ifdef PARALLEL

do kk=1,2
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfComSides
 else
  istart = grp%pdp%numberOfComSides + 1
  ifinish = grp%numberOfSides
 end if


 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%numberOfSides

#endif PARALLEL
  ind1 = grp%sideIndexArray(i,1)
  ind2 = grp%sideIndexArray(i,2)

  ww = grp%GCLWeightsArray(i)

  grp%uprev(ind1,1) = grp%uprev(ind1,1) + 2.0*ww
  grp%uprev(ind2,1) = grp%uprev(ind2,1) - 2.0*ww
 end do
#ifdef PARALLEL
 if(kk==1) then
  istart = 1
  ifinish = grp%pdp%numberOfBoundaryComSides
 else
  istart = grp%pdp%numberOfBoundaryComSides + 1
  ifinish = grp%brp%numberOfBoundarySides
 end if
 do i=istart,ifinish

#else PARALLEL

 do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL
  ind1 = grp%brp%sideIndexArray(i,1)
  ind2 = grp%brp%sideIndexArray(i,2)


  grp%uprev(ind1,1) = grp%uprev(ind1,1) + 4.0*grp%brp%GCLWeightsArray(i,1)
  grp%uprev(ind2,1) = grp%uprev(ind2,1) + 4.0*grp%brp%GCLWeightsArray(i,2)
 end do

#ifdef PARALLEL

 if (kk.eq.1) then
  do i=1,grp%pdp%numberOfProcesses
   if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then

    call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
    position = 1
     call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%uprev(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
    call mp_sendv(grp%pdp%processorIDs,i,145,grp%pdp%sendFlag)
   end if
  end do
 end if
end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,145,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
     call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                            grp%pdp%realType,grp%uprev(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
  end if
 end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position = 1
    call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                      grp%pdp%realType,grp%uprev(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
   call mp_sendv(grp%pdp%processorIDs,i,146,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,146,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
     call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%uprev(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
  end if
 end do
#endif PARALLEL
 if(ivd%flowtype.eq.2) then
   do i=1,grp%numberOfNodes
    if(ivd%timeScheme==1) then
      if(abs(grp%uprev(i,1)-(grp%nodeVolume(i)-grp%nodeVolumeP(i)))>1.0e-8) then
        write(*,*) "GCL error: ",i,grp%uprev(i,1),grp%nodeVolume(i)-grp%nodeVolumeP(i)
        write(*,*) "SS: ",grp%nodeVolume(i),grp%nodeVolumeP(i)
!   PAUSE
      end if
    else
      if(abs(grp%uprev(i,1)-(1.5*grp%nodeVolume(i)-2.0*grp%nodeVolumeP(i)+0.5*grp%nodeVolumeP2(i)))>1.0e-8) then
        write(*,*) "GCL error: ",i,grp%uprev(i,1),(1.5*grp%nodeVolume(i)-2.0*grp%nodeVolumeP(i)+0.5*grp%nodeVolumeP2(i))
        write(*,*) "RR: ",grp%nodeVolume(i),grp%nodeVolumeP(i),grp%nodeVolumeP2(i)
!   PAUSE
      end if
    end if
   end do
 end if
 end if

 end subroutine readGridComputationData
!------------------------------------------------------------------------
 subroutine calculateLiftAndDrag(grp,ivd,lift,drag,latF,frictionDRag,latM,longM)
 IMPLICIT NONE

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd
 real :: lift,drag,latF   !lift,drag,lateral force
 real :: latM,longM

 integer :: i,ist,ien,ip,ind1,ind2,j,MPI_IERR
 real :: flowTangent(3),flowNormal1(3),flowNormal2(3),wallTangent(3),wallNormal(3)
 real :: p,tau_t,tau_n,PI,alphaRad,betaRad
 real :: buff(3),normalStressComp,tangentialStress(3),T(2)

 real :: u2,v2,w2,dx(3),dxn,unorm,pos,neg,Minf,gamma
 real :: pressureDrag,pressureLift,frictionDrag
 real :: ftmp,height,xCGdist,yCGdist,zCGdist

 Minf = ivd%MachNumber
 gamma = ivd%gamma
 if(grp%gridNumber>1) STOP "ERROR: Shouldn't be here - 55"

 lift = 0.0
 drag = 0.0
 latF = 0.0 !initialise lateral force
 pressureDrag = 0.0
 pressureLift = 0.0
 frictionDrag = 0.0   !initialise skin friction drag component
 latM = 0.0 !initialise lateral moment 
 longM = 0.0 !initialise longitudinal moment


 PI = 4.0*atan(1.0)
 alphaRad = ivd%alpha*PI/180.
 betaRad = ivd%beta*PI/180.
 pos=0
 neg=0

 flowTangent(1) = cos(alphaRad)*cos(betaRad)
 flowTangent(2) = sin(betaRad)
 flowTangent(3) = sin(alphaRad)
 flowNormal1(1) = - flowTangent(3)   !normal in the lift direction
 flowNormal1(2) = 0.0
 flowNormal1(3) = flowTangent(1)
 flowNormal2(1) = - flowTangent(2)   !normal in the lateral direction
 flowNormal2(2) = flowTangent(1)
 flowNormal2(3) = 0.0

 ist = grp%brp%nodeIndicatorRegister(-20)
 ien = grp%brp%nodeIndicatorRegister(-19)
 do i=ist,ien
  ip = grp%brp%nodeIndicatorArray(i)
  height=grp%coordinates(ip,3)    ! height above ground BJE
  xCGdist = grp%coordinates(ip,1) - ivd%datum(1) !x-dist from CG
  yCGdist = grp%coordinates(ip,2) - ivd%datum(2) !y-dist from CG
  zCGdist = grp%coordinates(ip,3) - ivd%datum(3) !z-dist from CG
  if(height>(ivd%groundplane+0.01))then   !only include if above the ground

   ! pressure contribution
   if(ip>grp%numberOfNodes.or.ip<0) STOP
   lift = lift - grp%wallArea(ip)*(grp%p(ip)-(1./(gamma*Minf*Minf)))*sum(grp%brp%nodeNormalArray(ip,:)*flowNormal1(1:3))
   latF = latF - grp%wallArea(ip)*(grp%p(ip)-(1./(gamma*Minf*Minf)))*sum(grp%brp%nodeNormalArray(ip,:)*flowNormal2(1:3))
   drag = drag - grp%wallArea(ip)*(grp%p(ip)-(1./(gamma*Minf*Minf)))*sum(grp%brp%nodeNormalArray(ip,:)*flowTangent(1:3))

   latM = latM + xCGdist*grp%wallArea(ip)*(grp%p(ip)-(1./(gamma*Minf*Minf)))*sum(grp%brp%nodeNormalArray(ip,:)*flowNormal2(1:3))
   latM = latM - yCGdist*grp%wallArea(ip)*(grp%p(ip)-(1./(gamma*Minf*Minf)))*sum(grp%brp%nodeNormalArray(ip,:)*flowTangent(1:3))
   longM = longM + xCGdist*grp%wallArea(ip)*(grp%p(ip)-(1./(gamma*Minf*Minf)))*sum(grp%brp%nodeNormalArray(ip,:)*flowNormal1(1:3))
   longM = longM - zCGdist*grp%wallArea(ip)*(grp%p(ip)-(1./(gamma*Minf*Minf)))*sum(grp%brp%nodeNormalArray(ip,:)*flowTangent(1:3))

   pressureLift = pressureLift - grp%p(ip)*(grp%wallArea(ip)-(1./(gamma*Minf*Minf)))*sum(grp%brp%nodeNormalArray(ip,:)*flowNormal1(1:3))
   pressureDrag = pressureDrag - grp%p(ip)*(grp%wallArea(ip)-(1./(gamma*Minf*Minf)))*sum(grp%brp%nodeNormalArray(ip,:)*flowTangent(1:3))

  ! stress contribution

   if(ivd%ReynoldsNumber>0.0) then
    lift = lift +  grp%wallArea(ip)*sum(grp%wallStress(ip,1:3)*flowNormal1(:))
    latF = latF +  grp%wallArea(ip)*sum(grp%wallStress(ip,1:3)*flowNormal2(:))
    drag = drag +  grp%wallArea(ip)*sum(grp%wallStress(ip,1:3)*flowTangent(:))
    frictionDrag = frictionDrag +  grp%wallArea(ip)*sum(grp%wallStress(ip,1:3)*flowTangent(:))

    latM = latM - xCGdist*grp%wallArea(ip)*sum(grp%wallStress(ip,1:3)*flowNormal2(1:3))
    latM = latM + yCGdist*grp%wallArea(ip)*sum(grp%wallStress(ip,1:3)*flowTangent(:))
    longM = longM - xCGdist*grp%wallArea(ip)*sum(grp%wallStress(ip,1:3)*flowNormal1(1:3))
    longM = longM + zCGdist*grp%wallArea(ip)*sum(grp%wallStress(ip,1:3)*flowTangent(:))
   end if
  end if
 end do

! convert to lift,drag,LatF coefficients
 lift = 2.0*lift/ivd%referenceArea
 drag = 2.0*drag/ivd%referenceArea
 latF = 2.0*latF/ivd%referenceArea
 pressureLift = 2.0*pressureLift/ivd%referenceArea
 pressureDrag = 2.0*pressureDrag/ivd%referenceArea
 frictionDrag = 2.0*frictionDrag/ivd%referenceArea
 
 ! send to other processes

#ifdef PARALLEL

!  write(get_unit(),*)  'Starting Lift and Drag comms'
 call mp_init_buffer(grp%pdp%processorIDs,-1,grp%pdp%sendFlag)
 position = 1
 call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,lift,grp%pdp%sendFlag)
 call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,drag,grp%pdp%sendFlag)
 call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,latF,grp%pdp%sendFlag)
 call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,pressureLift,grp%pdp%sendFlag)
 call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,pressureDrag,grp%pdp%sendFlag)
 call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,frictionDrag,grp%pdp%sendFlag)
 call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,latM,grp%pdp%sendFlag)
 call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,longM,grp%pdp%sendFlag)
 
 call mp_snd_others(grp%pdp%numberOfProcesses,grp%pdp%currentDomain,grp%pdp%processorIDs,362,grp%pdp%sendFlag)
 
 call init_recv_others(grp%pdp%numberOfProcesses,grp%pdp%currentDomain,grp%pdp%receivedMessageFromProcess)
  
  call mp_recv_others(grp%pdp%numberOfProcesses,grp%pdp%currentDomain,grp%pdp%processorIDs,362,grp%pdp%sendFlag)
  
  call mp_wait_comms( grp%pdp%sendFlag )
!3771 continue
  do i=1,grp%pdp%numberOfProcesses
! if(grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if( i.ne.grp%pdp%currentDomain ) then
!  grp%pdp%sendFlag = 0
!  call mp_check_mesg_arrived(grp%pdp%processorIDs,i,362,grp%pdp%sendFlag,bufferzz,bufferSize)
!  if(grp%pdp%sendFlag.gt.0) then
!   grp%pdp%receivedMessageFromProcess(i) = 1
!   position = 1
    call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
    lift = lift + ftmp
    call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
    drag = drag + ftmp
    call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
    latF = latF + ftmp
    call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
    pressureLift = pressureLift + ftmp
    call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
    pressureDrag = pressureDrag + ftmp
    call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
    frictionDrag = frictionDrag + ftmp
    call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
    latM = latM + ftmp
    call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
    longM = longM + ftmp
!  endif
  endif
 end do
 if( grp%pdp%currentDomain.eq.1) then
 write(*,*) 'LIFT And Drag 1:'
 write(*,*) 'pressureDrag=',pressureDrag,'frictionDrag=',frictionDrag,'drag=',drag
 write(*,*) 'lift=',lift,'latF=',latF,' latM=',latM,' longM=',longM
 end if
 !  CALL MPI_ABORT(MPI_IERR) 
 !  STOP 

!do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0)                                  goto 3771
!end do
#endif PARALLEL
!
!
 end subroutine calculateLiftAndDrag
!------------------------------------------------------------------------
 subroutine calculateLiftAndDrag2(grp,ivd,lift,drag,latF,frictionDrag,latM,longM,ispan,localspan)!START-GUST
 IMPLICIT NONE

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd
 real :: lift,drag,latF   !lift,drag,lateral force
 real :: latM,longM

 integer :: i,ist,ien,ip,ind1,ind2,j,count,iw
 real :: flowTangent(3),flowNormal1(3),flowNormal2(3),wallTangent(3),wallNormal(3)
 real :: p,tau_t,tau_n,PI,alphaRad,betaRad,dist
 real :: buff(3),normalStressComp,tangentialStress(3),T(2)

 real :: u2,v2,w2,dx(3),dxn,unorm,frictiontest
 real :: frictionDrag,pressureDrag,pressureLift,height1,height2
 real :: ftmp,normArea(3),torque(4,3),vec(3)
 real :: p1,p2,area,totalArea,centre(3),radius,width,Minf,gamma
 real :: height,xCGdist,yCGdist,zCGdist,x,y,z
 real :: localspan !START-GUST
 integer :: ispan !START-GUST
 integer :: faceIndicator

 grp%nodeHelpArray(:,1:9) = 0.0
 torque(:,:) = 0.0
 pressureDrag = 0.0
 pressureLift = 0.0
 frictionDrag = 0.0   !initialise skin friction drag component

 totalArea = 0.0
 count = 0

 gamma = ivd%gamma
 Minf = ivd%MachNumber

 do i=1,grp%brp%numberOfBoundarySides
  ind1  = grp%brp%sideIndexArray(i,1)
  ind2  = grp%brp%sideIndexArray(i,2)

   faceIndicator = grp%brp%sideIndexArray(i,3)

!  if(faceIndicator==1.or.faceIndicator.eq.9) then
   if(faceIndicator==1) then
     height1 = grp%coordinates(ind1,3)   ! height above ground  BJE
     height2 = grp%coordinates(ind2,3)
    if((height1>(ivd%groundplane+0.01)).OR.(height2>(ivd%groundplane+0.01)))then   ! BJE modification

   ! pressure contribution

      count = count + 1

      p1 = grp%p(ind1) - (1./(gamma*Minf*Minf))
      p2 = grp%p(ind2) - (1./(gamma*Minf*Minf))
      grp%nodeHelpArray(ind1,1:3) = grp%nodeHelpArray(ind1,1:3) + 4.0*grp%brp%sideWeightsArray(i,:)*p1
      grp%nodeHelpArray(ind2,1:3) = grp%nodeHelpArray(ind2,1:3) + 4.0*grp%brp%sideWeightsArray(i,:)*p2
      grp%nodeHelpArray(ind1,4:6) = grp%nodeHelpArray(ind1,4:6) + 4.0*grp%brp%sideWeightsArray(i,:)*p1
      grp%nodeHelpArray(ind2,4:6) = grp%nodeHelpArray(ind2,4:6) + 4.0*grp%brp%sideWeightsArray(i,:)*p2

   ! viscous contribution
      if(faceIndicator.eq.1) then

        area = sqrt(sum(grp%brp%sideWeightsArray(i,1:3)*grp%brp%sideWeightsArray(i,1:3)))
        totalArea = totalArea + 4.0*area

        if(ivd%ReynoldsNumber>0.0) then
          grp%nodeHelpArray(ind1,1:3) = grp%nodeHelpArray(ind1,1:3) + 4.0*grp%wallStress(ind1,1:3)*area
          grp%nodeHelpArray(ind2,1:3) = grp%nodeHelpArray(ind2,1:3) + 4.0*grp%wallStress(ind2,1:3)*area
          grp%nodeHelpArray(ind1,7:9) = grp%nodeHelpArray(ind1,7:9) + 4.0*grp%wallStress(ind1,1:3)*area
          grp%nodeHelpArray(ind2,7:9) = grp%nodeHelpArray(ind2,7:9) + 4.0*grp%wallStress(ind2,1:3)*area
        end if
      end if

    end if
  end if
 end do


! find local lift and drag contribution

 lift = 0.0
 drag = 0.0
 latF = 0.0 !initialise lateral force
 latM = 0.0 !initialise lateral moment 
 longM = 0.0 !initialise longitudinal moment
 frictiontest = 0.0

 PI = 4.0*atan(1.0)
 alphaRad = ivd%alpha*PI/180.
 betaRad  = ivd%beta*PI/180.

 flowTangent(1) = cos(alphaRad)*cos(betaRad)
 flowTangent(2) = sin(betaRad)
 flowTangent(3) = sin(alphaRad)
 flowNormal1(1) = - flowTangent(3)   !normal in the lift direction
 flowNormal1(2) = 0.0
 flowNormal1(3) = flowTangent(1)
 flowNormal2(1) = - flowTangent(2)   !normal in the lateral direction
 flowNormal2(2) = flowTangent(1)
 flowNormal2(3) = 0.0


!do i=1,grp%brp%numberOfBoundaryNodes
 ist = grp%brp%nodeIndicatorRegister(-20)
 ien = grp%brp%nodeIndicatorRegister(-19)
 do i=ist,ien
  ip = grp%brp%nodeIndicatorArray(i)
!START-GUST
  if(grp%coordinates(ip,ivd%spanCor) .LE. localspan) then
    xCGdist = grp%coordinates(ip,1) - ivd%datum(1) !x-dist from CG
    yCGdist = grp%coordinates(ip,2) - ivd%datum(2) !y-dist from CG
    zCGdist = grp%coordinates(ip,3) - ivd%datum(3) !z-dist from CG
    lift = lift + sum(grp%nodeHelpArray(ip,1:3)*flowNormal1(1:3))
    drag = drag + sum(grp%nodeHelpArray(ip,1:3)*flowTangent(1:3))
    latF = latF + sum(grp%nodeHelpArray(ip,1:3)*flowNormal2(1:3))
    pressureLift = pressureLift + sum(grp%nodeHelpArray(ip,4:6)*flowNormal1(1:3))
    pressureDrag = pressureDrag + sum(grp%nodeHelpArray(ip,4:6)*flowTangent(1:3))
    frictionDrag = frictionDrag + sum(grp%nodeHelpArray(ip,7:9)*flowTangent(1:3))
    latM = latM - xCGdist*sum(grp%nodeHelpArray(ip,1:3)*flowNormal2(1:3))
    latM = latM + yCGdist*sum(grp%nodeHelpArray(ip,1:3)*flowTangent(1:3))
    longM = longM - xCGdist*sum(grp%nodeHelpArray(ip,1:3)*flowNormal1(1:3))
    longM = longM + zCGdist*sum(grp%nodeHelpArray(ip,1:3)*flowTangent(1:3))
  end if
!END-GUST

  if(ivd%nwheels.gt.0)then    !if there are wheels present and calculating the viscous torque
    do iw=1,ivd%nwheels
      centre(:)=ivd%wheelCentre(iw,:)
      radius = ivd%wheelDimensions(iw,1)
      width = ivd%wheelDimensions(iw,2)
      x = grp%coordinates(ip,1)
      y = grp%coordinates(ip,2)
      z = grp%coordinates(ip,3)
      dist = sqrt((x-centre(1))**2+(z-centre(3))**2)
      vec(:) = grp%coordinates(ip,:) - centre(:)
      if((y.le.(centre(2)+(width/2.0))).and.(y.ge.(centre(2)-(width/2.0))).and.(dist.le.radius))then
        torque(iw,1) =  torque(iw,1) + vec(2)*grp%nodeHelpArray(ip,9) - vec(3)*grp%nodeHelpArray(ip,8)
        torque(iw,2) =  torque(iw,2) + vec(3)*grp%nodeHelpArray(ip,7) - vec(1)*grp%nodeHelpArray(ip,9)
        torque(iw,3) =  torque(iw,3) + vec(1)*grp%nodeHelpArray(ip,8) - vec(2)*grp%nodeHelpArray(ip,7)
        frictiontest = frictiontest + sum(grp%nodeHelpArray(ip,7:9)*flowTangent(1:3))
      endif
    enddo
  endif

 end do

 ! communicate

#ifdef PARALLEL

!  write(get_unit(),*)  'Starting Lift and Drag 2 comms'
 call mp_init_buffer(grp%pdp%processorIDs,-1,grp%pdp%sendFlag)
 position = 1
 call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,lift,grp%pdp%sendFlag)
 call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,drag,grp%pdp%sendFlag)
 call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,latF,grp%pdp%sendFlag)
 call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,latM,grp%pdp%sendFlag)
 call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,longM,grp%pdp%sendFlag)
 call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,totalArea,grp%pdp%sendFlag)
 call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,pressureLift,grp%pdp%sendFlag)
 call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,pressureDrag,grp%pdp%sendFlag)
 call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,frictionDrag,grp%pdp%sendFlag)
 do iw = 1,ivd%nwheels
   call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,torque(iw,1),grp%pdp%sendFlag)
   call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,torque(iw,2),grp%pdp%sendFlag)
   call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,torque(iw,3),grp%pdp%sendFlag)
 enddo

 call mp_snd_others(grp%pdp%numberOfProcesses,grp%pdp%currentDomain,grp%pdp%processorIDs,95,grp%pdp%sendFlag)

 call init_recv_others(grp%pdp%numberOfProcesses,grp%pdp%currentDomain,grp%pdp%receivedMessageFromProcess)

  call mp_recv_others(grp%pdp%numberOfProcesses,grp%pdp%currentDomain,grp%pdp%processorIDs,95,grp%pdp%sendFlag)

  call mp_wait_comms( grp%pdp%sendFlag )
!3091 continue

  do i=1,grp%pdp%numberOfProcesses
! if(grp%pdp%receivedMessageFromProcess(i).eq.0) then
    if( i.ne.grp%pdp%currentDomain ) then
!  grp%pdp%sendFlag = 0
!  call mp_check_mesg_arrived(grp%pdp%processorIDs,i,95,grp%pdp%sendFlag,bufferzz,bufferSize)
!  if(grp%pdp%sendFlag.gt.0) then
!   grp%pdp%receivedMessageFromProcess(i) = 1
!   position = 1
    call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
    lift  = lift + ftmp
    call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
    drag  = drag + ftmp
    call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
    latF  = latF + ftmp
    call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
    latM  = latM + ftmp
    call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
    longM  = longM + ftmp
    call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
    totalArea  = totalArea + ftmp
    call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
    pressureLift  = pressureLift + ftmp
    call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
    pressureDrag  = pressureDrag + ftmp
    call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
    frictionDrag  = frictionDrag + ftmp
    do iw=1,ivd%nwheels
      call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
      torque(iw,1)  = torque(iw,1) + ftmp
      call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
      torque(iw,2)  = torque(iw,2) + ftmp
      call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag)
      torque(iw,3)  = torque(iw,3) + ftmp
    enddo
!  endif
  endif
 end do

!do i=1,grp%pdp%numberOfProcesses
! if (grp%pdp%receivedMessageFromProcess(i).eq.0)                                  goto 3091
!end do

#endif PARALLEL


! scale to give lift,drag,LatF coefficients

 lift = 2.0*lift/ivd%referenceArea
 drag = 2.0*drag/ivd%referenceArea
 latF = 2.0*latF/ivd%referenceArea
 latM = 2.0*latM/ivd%referenceArea
 longM = 2.0*longM/ivd%referenceArea
 pressureLift = 2.0*pressureLift/ivd%referenceArea
 pressureDrag = 2.0*pressureDrag/ivd%referenceArea
 frictionDrag = 2.0*frictionDrag/ivd%referenceArea
 torque = 2.0*torque/ivd%referenceArea
!START-GUST
 if(ispan.ne.0) then
   grp%localAerParam(ispan,1) = lift
   grp%localAerParam(ispan,2) = drag
   grp%localAerParam(ispan,3) = latF
   grp%localAerParam(ispan,4) = pressureLift
   grp%localAerParam(ispan,5) = pressureDrag
   grp%localAerParam(ispan,6) = frictionDrag
   grp%localAerParam(ispan,7) = latM !latM
   grp%localAerParam(ispan,8) = longM !longM
 end if
!END-GUST

 if( grp%pdp%currentDomain.eq.1) then
 write(*,*) 'LIFT And Drag 2:'
 write(*,*) 'pressureDrag=',pressureDrag,'frictionDrag=',frictionDrag,'drag=',drag
 write(*,*) 'lift=',lift,'latF=',latF,' latM=',latM,' longM=',longM
 do iw = 1,ivd%nwheels
 write(*,*) 'wheel ',iw
 write(*,*) 'torque =',torque(iw,:)
 enddo
 end if


 end subroutine calculateLiftAndDrag2
!-----------------------------------------------------------------------
subroutine calculateLiftAndDrag3(grp,ivd,lift,drag,latF,latM,longMl,longMd,rollM)
 IMPLICIT NONE
 
 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd
 real :: lift,drag,latF   !lift,drag,lateral force
 real :: latM,longMl,longMd,rollM,gamma,Minf

 integer :: i,ist,ien,ip,ind1,ind2,j,k,count
 real :: flowTangent(3),flowNormal1(3),flowNormal2(3),wallTangent(3),wallNormal(3)
 real :: p,tau_t,tau_n,PI,alphaRad,betaRad
 real :: buff(3),normalStressComp,tangentialStress(3),T(2)
 real :: lift1,lift2

 real :: u2,v2,w2,dx(3),dxn,unorm
 real :: pressureDrag,pressureLift,height1,height2
 real :: ftmp,normArea(3)
 real :: p1,p2,area,totalArea
 real :: height,xCGdist,yCGdist,zCGdist,x,y,z,r,angle
 integer :: faceIndicator
 
 grp%nodeHelpArray(:,1:3) = 0.0
 totalArea = 0.0
 count = 0 
 
 Minf = ivd%MachNumber
 gamma = ivd%gamma
 
 do i=1,grp%brp%numberOfBoundarySides
  ind1  = grp%brp%sideIndexArray(i,1)
  ind2  = grp%brp%sideIndexArray(i,2)
  
   faceIndicator = grp%brp%sideIndexArray(i,3)

   if(faceIndicator==1.or.faceIndicator.eq.9) then
     height1 = grp%coordinates(ind1,3)   ! height above ground  BJE
     height2 = grp%coordinates(ind2,3)   
    if((height1>(ivd%groundplane+0.01)).OR.(height2>(ivd%groundplane+0.01)))then   ! BJE modification
    
   ! pressure contribution

       count = count + 1

       p1 = grp%p(ind1) - (1./(gamma*Minf*Minf))
       p2 = grp%p(ind2) - (1./(gamma*Minf*Minf))
       grp%nodeHelpArray(ind1,1:3) = grp%nodeHelpArray(ind1,1:3) + 4.0*grp%brp%sideWeightsArray(i,:)*p1
       grp%nodeHelpArray(ind2,1:3) = grp%nodeHelpArray(ind2,1:3) + 4.0*grp%brp%sideWeightsArray(i,:)*p2

   ! viscous contribution
      if(ivd%ReynoldsNumber>0.0)then
       if(faceIndicator.eq.1) then

         area = sqrt(sum(grp%brp%sideWeightsArray(i,1:3)*grp%brp%sideWeightsArray(i,1:3)))
         totalArea = totalArea + 4.0*area

         grp%nodeHelpArray(ind1,1:3) = grp%nodeHelpArray(ind1,1:3) + 4.0*grp%wallStress(ind1,1:3)*area
         grp%nodeHelpArray(ind2,1:3) = grp%nodeHelpArray(ind2,1:3) + 4.0*grp%wallStress(ind2,1:3)*area

       end if
      endif

     end if
   end if
  end do


! find local lift and drag contribution

 lift = 0.0
 drag = 0.0
 latF = 0.0 !initialise lateral force
 latM = 0.0 !initialise lateral moment 
 longMl = 0.0 !initialise longitudinal moment
 longMd = 0.0 !initialise longitudinal moment
 rollM = 0.0 !initialise roll moment

 PI = 4.0*atan(1.0)
 alphaRad = ivd%alpha*PI/180.
 betaRad  = ivd%beta*PI/180.
 flowTangent(1) = cos(alphaRad)*cos(betaRad) 
 flowTangent(2) = sin(betaRad) 
 flowTangent(3) = sin(alphaRad) 
 flowNormal1(1) = - flowTangent(3)   !normal in the lift direction 
 flowNormal1(2) = 0.0 
 flowNormal1(3) = flowTangent(1) 
 flowNormal2(1) = - flowTangent(2)   !normal in the lateral direction 
 flowNormal2(2) = flowTangent(1) 
 flowNormal2(3) = 0.0 
 
 
 !do i=1,grp%brp%numberOfBoundaryNodes 
    ist = grp%brp%nodeIndicatorRegister(-20) 
    ien = grp%brp%nodeIndicatorRegister(-19) 
    do i =ist,ien 
  ip = grp%brp%nodeIndicatorArray(i) 
  xCGdist = grp%coordinates(ip,1) - ivd%datum(1) !x-dist from CG 
  yCGdist = grp%coordinates(ip,2) - ivd%datum(2) !y-dist from CG 
  zCGdist = grp%coordinates(ip,3) - ivd%datum(3) !z-dist from CG 
  x = grp%coordinates(ip,1) 
  y = grp%coordinates(ip,2) 
  z = grp%coordinates(ip,3) 
  lift = lift + sum(grp%nodeHelpArray(ip,1:3)*flowNormal1(1:3)) 
  drag = drag + sum(grp%nodeHelpArray(ip,1:3)*flowTangent(1:3)) 
  latF = latF + sum(grp%nodeHelpArray(ip,1:3)*flowNormal2(1:3)) 
  latM = latM - xCGdist*sum(grp%nodeHelpArray(ip,1:3)*flowNormal2(1:3)) 
  latM = latM + yCGdist*sum(grp%nodeHelpArray(ip,1:3)*flowTangent(1:3)) 
  rollM = rollM - yCGdist*sum(grp%nodeHelpArray(ip,1:3)*flowNormal1(1:3)) 
  rollM = rollM + zCGdist*sum(grp%nodeHelpArray(ip,1:3)*flowNormal2(1:3)) 
  longMl = longMl - xCGdist*sum(grp%nodeHelpArray(ip,1:3)*flowNormal1(1:3)) 
  longMd = longMd + zCGdist*sum(grp%nodeHelpArray(ip,1:3)*flowTangent(1:3)) 
! 
 end do 
! communicate 
 
#ifdef PARALLEL 
 
!  write(get_unit(),*)  'Starting Lift and Drag 2 comms' 
 call mp_init_buffer(grp%pdp%processorIDs,-1,grp%pdp%sendFlag) 
 position = 1 
 call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,lift,grp%pdp%sendFlag) 
 call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,drag,grp%pdp%sendFlag) 
 call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,latF,grp%pdp%sendFlag) 
 call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,latM,grp%pdp%sendFlag)
 call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,longMl,grp%pdp%sendFlag) 
 call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,longMd,grp%pdp%sendFlag) 
 call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,totalArea,grp%pdp%sendFlag) 
 call mp_pak_real(grp%pdp%processorIDs,-1,grp%pdp%realType,rollM,grp%pdp%sendFlag)  
 
 
 call mp_snd_others(grp%pdp%numberOfProcesses,grp%pdp%currentDomain,grp%pdp%processorIDs,95,grp%pdp%sendFlag) 
 
 call init_recv_others(grp%pdp%numberOfProcesses,grp%pdp%currentDomain,grp%pdp%receivedMessageFromProcess)
 
  call mp_recv_others(grp%pdp%numberOfProcesses,grp%pdp%currentDomain,grp%pdp%processorIDs,95,grp%pdp%sendFlag) 
 
  call mp_wait_comms( grp%pdp%sendFlag ) 
!3091 continue 
 
  do i=1,grp%pdp%numberOfProcesses 
! if(grp%pdp%receivedMessageFromProcess(i).eq.0) then 
    if( i.ne.grp%pdp%currentDomain ) then 
!  grp%pdp%sendFlag = 0 
!  call mp_check_mesg_arrived(grp%pdp%processorIDs,i,95,grp%pdp%sendFlag,bufferzz,bufferSize) 
!  if(grp%pdp%sendFlag.gt.0) then 
!   grp%pdp%receivedMessageFromProcess(i) = 1 
!   position = 1 
    call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag) 
    lift  = lift + ftmp 
    call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag) 
    drag  = drag + ftmp 
    call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag) 
    latF  = latF + ftmp 
    call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag) 
    latM  = latM + ftmp 
    call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag) 
    longMl  = longMl + ftmp 
    call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag) 
    longMd  = longMd + ftmp 
    call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag) 
    totalArea  = totalArea + ftmp 
    call mp_upak_real(grp%pdp%processorIDs,i,grp%pdp%realType,ftmp,grp%pdp%sendFlag) 
    rollM  = rollM + ftmp 
!  endif
  endif 
 end do 
 
!do i=1,grp%pdp%numberOfProcesses 
! if (grp%pdp%receivedMessageFromProcess(i).eq.0)                                  goto 3091 
!end do 
 
#endif PARALLEL 
 
 
! scale to give lift,drag,latF coefficients 
 
 lift = 2.0*lift/ivd%referenceArea 
 drag = 2.0*drag/ivd%referenceArea 
 latF = 2.0*latF/ivd%referenceArea 
 latM = 2.0*latM/ivd%referenceArea 
 longMl = 2.0*longMl/ivd%referenceArea 
 longMd = 2.0*longMd/ivd%referenceArea 
 rollM = 2.0*rollM/ivd%referenceArea 
! 
 
 end subroutine calculateLiftAndDrag3 
!----------------------------------------------------------------------- 
 subroutine writeGridResults(grp,ivd,OUTFILE)

 IMPLICIT NONE

 type(GridSolverData) :: grp
 type(InputVariablesData) :: ivd
 integer :: OUTFILE

 integer :: i,j,ist,ien,ip,ind1,ind2
 real :: u2,v2,w2,unorm,dx(3),dxn


! do i=1,grp%numberOfNodes
!  if(grp%u(i,1)<ivd%minimumDensity) write(*,*) "small density: ",i,grp%u(i,1),grp%coordinates(i,:)
! end do

  write(OUTFILE) grp%numberOfNodes
!START-TURB
 if(ivd%turbulenceModel==0) then
  write(OUTFILE) (grp%u(i,1),i=1,grp%numberOfNodes),&
                 ((grp%u(i,j)/grp%u(i,1),i=1,grp%numberOfNodes),j=2,5)
end if
if(ivd%turbulenceModel==1) then
  write(OUTFILE) (grp%u(i,1),i=1,grp%numberOfNodes),&
                 ((grp%u(i,j)/grp%u(i,1),i=1,grp%numberOfNodes),j=2,6)
end if

if(ivd%turbulenceModel==2 .or. ivd%turbulenceModel==3) then
  write(OUTFILE)(grp%u(i,1),i=1,grp%numberOfNodes),&
                ((grp%u(i,j)/grp%u(i,1),i=1,grp%numberOfNodes),j=2,5),((grp%u(i,j),i=1,grp%numberOfNodes),j=6,7)
end if
!END-TURB


 if(.false.) then
 if(ivd%ReynoldsNumber>0.0) then

  open(34,file="skinfr.res",form='formatted')
  ist = grp%brp%nodeIndicatorRegister(-20)
  ien = grp%brp%nodeIndicatorRegister(-19)

  grp%nodeHelpArray(:,1:2) = 0.0

  do i=1,grp%numberOfSides
   ind1 = grp%sideIndexArray(i,1)
   ind2 = grp%sideIndexArray(i,2)

   if(ind1.le.grp%brp%numberOfBoundaryNodes.and.ind2>grp%brp%numberOfBoundaryNodes) then
    u2 = grp%u(ind2,2)/grp%u(ind2,1)
    v2 = grp%u(ind2,3)/grp%u(ind2,1)
    w2 = grp%u(ind2,4)/grp%u(ind2,1)

    unorm = sqrt(u2*u2+v2*v2+w2*w2)

    if(u2<0.0) unorm = -unorm

    dx = grp%coordinates(ind1,:) - grp%coordinates(ind2,:)

    dxn = sqrt(sum(dx*dx))

    grp%nodeHelpArray(ind1,1) = grp%nodeHelpArray(ind1,1) + unorm/dxn
    grp%nodeHelpArray(ind1,2) = grp%nodeHelpArray(ind1,2) + 1.0
   end if
  end do
 end if

  do i=ist,ien
   ip = grp%brp%nodeIndicatorArray(i)
   if(grp%nodeHelpArray(ip,2)>1.0e-10) then
    write(34,'(I4,E14.5,E14.5,E14.5,E14.5)') ip,grp%coordinates(ip,:),2.0*grp%laminarViscosity(ip)*grp%nodeHelpArray(ip,1)/&
                                                  grp%nodeHelpArray(ip,2)
   end if
  end do
 end if


 end subroutine writeGridResults
!-----------------------------------------------------------------------

  subroutine makeBoundaryNormals(grp)
! makes normals on boundaries using boundary edge coefficients
  IMPLICIT NONE
  include 'mpif.h'
  type(GridSolverData) :: grp
 
  integer :: i,j,ind1,ind2,allocateStatus
  real :: norm

  integer :: kk,istart,ifinish

  allocate(grp%brp%nodeNormalArray(grp%numberOfNodes,3),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: Not enough memory in readBoundarySolverData"
 
  if(grp%pdp%currentDomain==1) write(*,*) "making boundary normals..."

  grp%brp%nodeNormalArray = 0.0


#ifdef PARALLEL
  do kk=1,2
    if(kk==1) then
      istart = 1
      ifinish = grp%pdp%numberOfBoundaryComSides
    else
      istart = grp%pdp%numberOfBoundaryComSides + 1
      ifinish = grp%brp%numberOfBoundarySides
    end if

    do i=istart,ifinish

#else PARALLEL

    do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL
      ind1  = grp%brp%sideIndexArray(i,1)
      ind2  = grp%brp%sideIndexArray(i,2)
      grp%brp%nodeNormalArray(ind1,:) = grp%brp%nodeNormalArray(ind1,:) - grp%brp%sideWeightsArray(i,:)
      grp%brp%nodeNormalArray(ind2,:) = grp%brp%nodeNormalArray(ind2,:) - grp%brp%sideWeightsArray(i,:)
    end do

#ifdef PARALLEL

    if (kk.eq.1) then
      do i=1,grp%pdp%numberOfProcesses
        if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i)>0)) then
          call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
          position = 1
          do j=1,3
            call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                               grp%pdp%realType,grp%brp%nodeNormalArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
          end do
          call mp_sendv(grp%pdp%processorIDs,i,25,grp%pdp%sendFlag)
        end if
      end do
    end if
  end do

  call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,25,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

  do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
        do j=1,3
          call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                                 grp%pdp%realType,grp%brp%nodeNormalArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
        end do
    end if
  end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

  do i=1,grp%pdp%numberOfProcesses
    if (i.ne.grp%pdp%currentDomain.and. (grp%pdp%numberOfReceiveNodes(i)>0)) then
      call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
      position = 1
      do j=1,3
        call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                          grp%pdp%realType,grp%brp%nodeNormalArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
      end do
      call mp_sendv(grp%pdp%processorIDs,i,26,grp%pdp%sendFlag)
    end if
  end do

  call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,26,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

  do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
        do j = 1,3
          call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                             grp%pdp%realType,grp%brp%nodeNormalArray(:,j),grp%pdp%buffer,grp%pdp%sendFlag)
        end do
    end if
  end do

#endif PARALLEL


  do i=1,grp%brp%numberOfBoundaryNodes
   norm = sqrt(sum(grp%brp%nodeNormalArray(i,:)*grp%brp%nodeNormalArray(i,:)))
   if(norm>1.0e-10) then
    grp%brp%nodeNormalArray(i,:) = grp%brp%nodeNormalArray(i,:)/norm
   else
    write(*,*) "WARNING: Small normal ",i,norm
    grp%brp%nodeNormalArray(i,:)  =  0.0
   end if
  end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )
!  call MPI_ABORT( MPI_COMM_WORLD, i )

 end subroutine makeBoundaryNormals
!------------------------------------------------------------------------
 subroutine makeSurfaceArea(grp)
IMPLICIT NONE
include 'mpif.h'

 type(GridSolverData) :: grp

 integer :: i,ind1,ind2,kk,istart,ifinish
 real :: areaIncrement

 real :: xc(3)

 grp%wallArea = 0.0
 grp%nodeHelpArray(:,1) = 0.0

! make surface area
 if(grp%gridNumber==1) then
#ifdef PARALLEL
  do kk=1,2
   if(kk==1) then
    istart = 1
    ifinish = grp%pdp%numberOfBoundaryComSides
   else
    istart = grp%pdp%numberOfBoundaryComSides + 1
    ifinish = grp%brp%numberOfBoundarySides
   end if
   do i=istart,ifinish

#else PARALLEL

   do i=1,grp%brp%numberOfBoundarySides

#endif PARALLEL
    ind1 = grp%brp%sideIndexArray(i,1)
    ind2 = grp%brp%sideIndexArray(i,2)
    areaIncrement = sqrt(sum(grp%brp%sideWeightsArray(i,:)*grp%brp%sideWeightsArray(i,:)))
    grp%nodeHelpArray(ind1,1) = grp%nodeHelpArray(ind1,1) +  4.0*areaIncrement
    grp%nodeHelpArray(ind2,1) = grp%nodeHelpArray(ind2,1) +  4.0*areaIncrement
    grp%wallArea(ind1) =  grp%wallArea(ind1) +  4.0*areaIncrement
    grp%wallArea(ind2) =  grp%wallArea(ind2) +  4.0*areaIncrement
   end do

! send domain boundary residuals

  if (kk.eq.1) then
   do i=1,grp%pdp%numberOfProcesses
    if (i.ne.grp%pdp%currentDomain.and.grp%pdp%numberOfSendNodes(i)>0) then
     call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
     position = 1
     call mp_pakv_nabor(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                        grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
     call mp_sendv(grp%pdp%processorIDs,i,74,grp%pdp%sendFlag)
    end if
   end do
  end if
 end do

 call init_recv_bdry(grp%pdp%numberOfProcesses,grp%pdp%numberOfReceiveNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,74,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )


 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i) > 0 ) ) then
    call mp_upakv_add_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                           grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
  end if
 end do

  call MPI_BARRIER( MPI_COMM_WORLD, i )

 do i=1,grp%pdp%numberOfProcesses
  if (i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfReceiveNodes(i)>0)) then
   call mp_init_buffer(grp%pdp%processorIDs,i,grp%pdp%sendFlag)
   position = 1
   call mp_pakv_bdry(grp%pdp%processorIDs,i,grp%pdp%numberOfReceiveNodes(i),grp%pdp%receiveNodeRegister(i,:),&
                     grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
   call mp_sendv(grp%pdp%processorIDs,i,75,grp%pdp%sendFlag)
  end if
 end do

 call init_recv_nabr(grp%pdp%numberOfProcesses,grp%pdp%numberOfSendNodes,grp%pdp%receivedMessageFromProcess)

  do i  = 1, grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
      call mp_recv( grp%pdp%processorIDs,i,75,grp%pdp%sendFlag )
    endif
  enddo

  call mp_wait_comms( grp%pdp%sendFlag )

 do i=1,grp%pdp%numberOfProcesses
    if(i.ne.grp%pdp%currentDomain.and.(grp%pdp%numberOfSendNodes(i) > 0 ) ) then
    call mp_upakv_nabr(grp%pdp%processorIDs,i,grp%pdp%numberOfSendNodes(i),grp%pdp%sendNodeRegister(i,:),&
                       grp%pdp%realType,grp%nodeHelpArray(:,1),grp%pdp%buffer,grp%pdp%sendFlag)
  end if
 end do

! place in wallArea array

!do i=1,grp%brp%numberOfBoundaryNodes
! grp%wallArea(i) = grp%nodeHelpArray(i,1)
!end do
end if


 end subroutine makeSurfaceArea
!------------------------------------------------------------------------
 integer function nameLen(fn)
 IMPLICIT NONE

 character*120 :: fn
 integer :: i

 do i=120,1,-1
  nameLen = i
  if(fn(i:i)/=' ') GOTO 77 ! EXIT
 end do
 nameLen = 0
77  end function nameLen
!-----------------------------------------------------------------------
 subroutine startTiming(grp,number)
 ! for timing purposes
 IMPLICIT NONE

  type(GridSolverData) :: grp
  integer :: number

  real(KIND=4), external :: etime
  real(KIND=4) :: userandsysTime(2)
  integer :: message

!  write(*,*) "starting timing number: ",number
!  grp%startTime = etime(userandsysTime)

 end subroutine startTiming
!-----------------------------------------------------------------------
 subroutine stopTiming(grp,number)
 ! for timing purposes
 IMPLICIT NONE
 type(GridSolverData) :: grp
 integer :: number
 real(KIND=4), external :: etime
 real(KIND=4) :: userandsysTime(2)

!  grp%stopTime = etime(userandsysTime)
!  write(*,*) "timing number ",number,": ",grp%stopTime-grp%startTime

 end subroutine stopTiming
!-----------------------------------------------------------------------
end module GridSolver

