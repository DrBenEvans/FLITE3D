!**********************************************
!*Aggl3d_Struct.f90                           *
!*                                            *
!* Side and coordinate modules definition file*
!*                                            *
!*   Aggl3d v 2.3.2                           *
!*                                            *
!*24.06.99-                                   *
!**********************************************
!*Description:                                *
!* File contains coordinate and side          *
!* datastructure                              *
!**********************************************

module CoordinateRegister
! keeps track of coordinate points

 type CoordinateRegisterData ! coordinate register data
  real, pointer :: points(:,:)
  integer :: numberOfPoints,NSD
 end type CoordinateRegisterData
 
 integer, parameter, private :: NSD = 3

 contains

!-------------------------------------------------------------------------
 subroutine cleanUpCoordinateData(crp)
 ! deallocates coordinate data

 IMPLICIT NONE

 type(CoordinateRegisterData) :: crp

 integer :: allocateStatus

  deallocate(crp%points,stat=allocateStatus)
  if(allocateStatus>0) STOP "ERROR: could not deallocate in cleanUpCoordinateData"
  
 end subroutine cleanUpCoordinateData
!-------------------------------------------------------------------------
 subroutine constructCoordinateRegister(np,nsd,crp)
! sets up register

 IMPLICIT NONE

 integer :: np,nsd
 type(CoordinateRegisterData) :: crp

 integer :: allocateStatus

 crp%NSD = nsd
 crp%numberOfPoints = np

 allocate(crp%points(np,nsd),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create point array"

 end subroutine constructCoordinateRegister
!-------------------------------------------------------------------------
 subroutine readCoordinateData(INFILE,crp)
! reads data from file

 IMPLICIT NONE

 integer :: INFILE
 type(CoordinateRegisterData) :: crp

 integer ::  ip,i,dummy
 
 integer NSD,NP

 NSD = crp%NSD
 NP = crp%numberOfPoints

 read(INFILE) ((crp%points(ip,i),ip=1,NP),i=1,NSD)

 
 end subroutine readCoordinateData
!-------------------------------------------------------------------------
 subroutine setUpPreviousCoordinates(crpPrev,crp,hasData)
 IMPLICIT NONE

 type(CoordinateRegisterData) :: crpPrev,crp
 logical :: hasData

 integer :: i


 
 if(hasData) then
  do i=1,crp%numberOfPoints
   crpPrev%points(i,:) = -crpPrev%points(i,:) + crp%points(i,:)
  end do
 else
  do i=1,crp%numberOfPoints
   crpPrev%points(i,:) = crp%points(i,:)
  end do
 end if
 end subroutine setUpPreviousCoordinates
!------------------------------------------------------------------------- 
 subroutine writeCoordinateData(OUTFILE,crp)
! writes coordinate data for computation file

 IMPLICIT NONE

 type(CoordinateRegisterData) :: crp
 integer :: OUTFILE

 integer :: i,j

  write(OUTFILE) ((crp%points(i,j),i=1,crp%numberOfPoints),j=1,NSD)

 end subroutine writeCoordinateData
!-------------------------------------------------------------------------
 function getCoor(i,crp) result (c)
! returns coordinates related to a node index
 IMPLICIT NONE

 integer :: i
 type(CoordinateRegisterData) :: crp

 real :: c(NSD)

 c = crp%points(i,:)

 end function getCoor
!-------------------------------------------------------------------------
 subroutine setCoor(i,c,crp) 
! returns coordinates related to a node index
 IMPLICIT NONE

 integer :: i
 type(CoordinateRegisterData) :: crp
 real :: c(NSD)

 crp%points(i,:) = c
 end subroutine setCoor
!-------------------------------------------------------------------------
end module CoordinateRegister

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

module SideModule

 type SideData
  integer :: sideIndexes(2) ! mapping between side indexes and point indexes
  real :: sideCoefficients(3)
  real :: GCLCoefficient,GCLCoefficientP
  real :: sideLength
  integer :: numberOfSideComponents ! for use to calculate average side length
 end type SideData

 contains

!--------------------------------------------------------------------------
  subroutine sideConstruct(n,isFirst,sd)
! initializes side

  IMPLICIT NONE

  integer :: n
  logical :: isFirst
  type(SideData),pointer :: sd

  integer :: i,allocateStatus

  sd%sideCoefficients = 0.

  end subroutine
!-------------------------------------------------------------------------

end module SideModule

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

module VisualizationModule
! does the visualization bit
 type VisualizationSide
  integer :: sideIndexes(2)
  integer :: neighbouringElements(2)
  integer :: controlVolumes(2)
 end type VisualizationSide

 type VisualizationSidePointer
  type(VisualizationSide),pointer :: vs
  type(VisualizationSidePointer),pointer :: next
 end type VisualizationSidePointer

 type VisualizationSidePointers
  type(VisualizationSidePointer),pointer :: first,last
 end type VisualizationSidePointers

 type VisualizationData
  type(VisualizationSidePointers),pointer :: vsp(:)
  integer :: numberOfVisualizationSides
 end type VisualizationData

end module VisualizationModule

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

module SideRegister
! register for sides
 use SideModule

 type SideDataPointer
  type(SideData),pointer :: sd
 end type

 type SideDataBlockInstance
  type(SideData),pointer :: block(:)
  type(SideDataBlockInstance),pointer :: next
 end type SideDataBlockInstance

 type OrderedSideDataBlockChain
  type(SideDataBlockInstance),pointer :: first
  type(SideDataBlockInstance),pointer :: last
 end type OrderedSideDataBlockChain

 type NodeBlockInstance
  integer,pointer :: block(:,:)
  type(NodeBlockInstance),pointer :: next
 end type NodeBlockInstance

 type OrderedNodeBlockChain
  type(NodeBlockInstance),pointer :: first
  type(NodeBlockInstance),pointer :: last
 end type OrderedNodeBlockChain


 type SideRegisterData
  type(OrderedSideDataBlockChain),pointer :: dataBlockChain
  type(OrderedNodeBlockChain),pointer :: nodeBlockChain
  integer,pointer :: nodeStartRegister(:)
  integer,pointer :: nodeSideArray(:)
  integer,pointer :: nodeSideArrayRegister(:)
  integer,pointer :: parallelizationColorArray(:)
  integer :: dataBlockSize
  integer :: numberOfDataBlocks
  integer :: numberOfSides
  integer :: dataBlockIndex
  integer :: nodeBlockSize
  integer :: numberOfNodeBlocks
  integer :: numberOfNodes
  integer :: nodeBlockIndex
 end type SideRegisterData

 contains
!-------------------------------------------------------------------------
 subroutine setUpSideRegister(srd,numberOfNodes,dataBlockSize,nodeBlockSize)
 IMPLICIT NONE

 type(SideRegisterData) :: srd
 integer :: numberOfNodes,dataBlockSize,nodeBlockSize,allocateStatus

  srd%dataBlockSize = dataBlockSize
  srd%nodeBlockSize = nodeBlockSize
  srd%numberOfDataBlocks = 0
  srd%numberOfNodeBlocks = 0
  srd%numberOfSides = 0
  srd%dataBlockIndex = dataBlockSize
  srd%nodeBlockIndex = nodeBlockSize
  srd%numberOfNodes = numberOfNodes
  allocate(srd%nodeStartRegister(numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory in setUpSideRegister"
  allocate(srd%dataBlockChain,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory in setUpSideRegister"
  allocate(srd%nodeBlockChain,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory in setUpSideRegister"
  nullify(srd%dataBlockChain%first)
  nullify(srd%dataBlockChain%last)
  nullify(srd%nodeBlockChain%first)
  nullify(srd%nodeBlockChain%last)
  srd%nodeStartRegister = 0
 end subroutine setUpSideRegister
!-------------------------------------------------------------------------
 type(SideData) function getSideNumber(sideno,srd)
 IMPLICIT NONE

 integer :: sideno
 type(SideRegisterData) :: srd

 integer :: block,index,i
 type(SideDataBlockInstance),pointer :: currentSideBlock

  block = (sideno-1)/srd%dataBlockSize + 1
  index = sideno - (block-1)*srd%dataBlockSize
  ! shift to correct block
  currentSideBlock => srd%dataBlockChain%first
  do i=2,block
   currentSideBlock => currentSideBlock%next
  end do
  getSideNumber = currentSideBlock%block(index)
 end function getSideNumber
!-------------------------------------------------------------------------
 subroutine addSideToList(sd,srd)
 IMPLICIT NONE

 type(SideData) :: sd
 type(SideRegisterData) :: srd

 type(SideDataBlockInstance),pointer :: currentData,searchData
 type(NodeBlockInstance),pointer :: currentNode
 integer :: allocateStatus,i,j,ind1,ind2

 ! first add side to list

 if(srd%dataBlockIndex == srd%dataBlockSize) then
  ! do some allocation
  allocate(currentData,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory in addSideToList"
  allocate(currentData%block(srd%dataBlockSize),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory in addSideToList"
  ! add new block to end of chain
  nullify(currentData%next)
  if(associated(srd%dataBlockChain%first)) then
   srd%dataBlockChain%last%next=>currentData
   srd%dataBlockChain%last=>currentData
  else
   srd%dataBlockChain%first=>currentData
   srd%dataBlockChain%last=>currentData
  end if
  srd%numberOfDataBlocks = srd%numberOfDataBlocks + 1
  srd%dataBlockIndex = 1
 else
  srd%dataBlockIndex = srd%dataBlockIndex + 1
 end if
 srd%numberOfSides = srd%numberOfSides + 1

 srd%dataBlockChain%last%block(srd%dataBlockIndex) = sd
 ! now update register where the sides connected to each node are
 ! registered (to speed up constroction and agglomeration)

 do i=1,2
  if(srd%nodeBlockIndex == srd%nodeBlockSize) then
   ! do some allocation
   allocate(currentNode,stat=allocateStatus)
   if (allocateStatus /= 0) STOP "ERROR: Not enough memory in addSideToList"
   allocate(currentNode%block(srd%nodeBlockSize,2),stat=allocateStatus)
   if (allocateStatus /= 0) STOP "ERROR: Not enough memory in addSideToList"
   currentNode%block = 0
   nullify(currentNode%next)
   if(associated(srd%nodeBlockChain%first)) then
    srd%nodeBlockChain%last%next=>currentNode
    srd%nodeBlockChain%last=>currentNode
   else
    srd%nodeBlockChain%first=>currentNode
    srd%nodeBlockChain%last=>currentNode
   end if

   srd%numberOfNodeBlocks = srd%numberOfNodeBlocks + 1
   srd%nodeBlockIndex = 1
  else
   srd%nodeBlockIndex = srd%nodeBlockIndex + 1
  end if
  srd%nodeblockChain%last%block(srd%nodeBlockIndex,1) = srd%numberOfSides
  srd%nodeBlockChain%last%block(srd%nodeBlockIndex,2) = srd%nodeStartRegister(sd%sideIndexes(i))
  srd%nodeStartRegister(sd%sideIndexes(i)) = srd%nodeBlockSize*(srd%numberOfNodeBlocks-1) + srd%nodeBlockIndex
 end do

 end subroutine addSideToList
!-------------------------------------------------------------------------
 integer function registerSideData(isd,srd)
 IMPLICIT NONE

 type(SideData) :: isd
 type(SideRegisterData) :: srd

 integer :: ind1,ind2
 real :: sideLength
 integer :: i1,i2,dataSearchIndex,nextSidePosition,block,index,i
 type(SideData) :: sd
 real :: C12(3),GCL,GCLP
 type(SideDataBlockInstance),pointer :: currentSideBlock
 type(NodeBlockInstance),pointer :: currentNodeBlock


 ind1 = isd%sideIndexes(1)
 ind2 = isd%sideIndexes(2)
 sideLength = isd%sideLength

 i1 = min(ind1,ind2)
 i2 = ind1+ind2-i1

 if(i1==ind1) then
  C12 = isd%sideCoefficients
  GCL = isd%GCLCoefficient
  GCLP = isd%GCLCoefficientP 
 else
  C12 = -isd%sideCoefficients
  GCL = -isd%GCLCoefficient
  GCLP = -isd%GCLCoefficientP 
 end if

 dataSearchIndex = 0
 nextSidePosition = srd%nodeStartRegister(i1)
 if(nextSidePosition>0) then
  block = (nextSidePosition-1)/srd%dataBlockSize + 1
  index = nextSidePosition - (block-1)*srd%nodeBlockSize
  ! shift to correct block
  currentNodeBlock => srd%nodeBlockChain%first
  do i=2,block
   currentNodeBlock => currentNodeBlock%next
  end do
  dataSearchIndex = currentNodeBlock%block(index,1)
  nextSidePosition = currentNodeBlock%block(index,2)

  do while(dataSearchIndex>0)
   block = (dataSearchIndex-1)/srd%dataBlockSize + 1
   index = dataSearchIndex - (block-1)*srd%dataBlockSize
   if(index.le.0) STOP "ERROR: Something's wrong in registerSideData - 1"
   ! shift to correct block
   currentSideBlock => srd%dataBlockChain%first
   do i=2,block
    currentSideBlock => currentSideBlock%next
   end do
   if( currentSideBlock%block(index)%sideIndexes(2)==i2) then
    ! side exists
    currentSideBlock%block(index)%sideCoefficients = &
      currentSideBlock%block(index)%sideCoefficients + C12
    currentSideBlock%block(index)%GCLCoefficient = &
      currentSideBlock%block(index)%GCLCoefficient + GCL
    currentSideBlock%block(index)%GCLCoefficientP = &
      currentSideBlock%block(index)%GCLCoefficientP + GCLP
    currentSideBlock%block(index)%numberOfSideComponents =&
     currentSideBlock%block(index)%numberOfSideComponents + 1
    registerSideData = dataSearchIndex  
    exit
   end if

   if(nextSidePosition==0) then
    dataSearchIndex = 0
    exit
   end if
   block = (nextSidePosition-1)/srd%dataBlockSize + 1
   index = nextSidePosition - (block-1)*srd%nodeBlockSize
   ! shift to correct block
   currentNodeBlock => srd%nodeBlockChain%first
   do i=2,block
    currentNodeBlock => currentNodeBlock%next
   end do
   dataSearchIndex = currentNodeBlock%block(index,1)
   nextSidePosition = currentNodeBlock%block(index,2)
  end do
 end if

 if(dataSearchIndex==0) then
  ! must create new side

  sd%sideIndexes(1) = i1
  sd%sideIndexes(2) = i2
  sd%sideCoefficients = C12
  sd%GCLCoefficient = GCL
  sd%GCLCoefficientP = GCLP
  sd%sideLength = sideLength
  sd%numberOfSideComponents = 1
  call addSideToList(sd,srd)
  registerSideData = 0
 end if
 end function registerSideData
!-------------------------------------------------------------------------
 subroutine prepareDataStructure(srd)
 ! put all side numbers connected to nodes in one big array
 ! create a register which points to the starting indices for each node in the above array
 IMPLICIT NONE

  type(SideRegisterData) :: srd

  integer :: allocateStatus,counter,i,j,registerSize,value,ind1,ind2
  type(NodeBlockInstance),pointer :: currentNodeBlock
  type(SideDataBlockInstance),pointer :: currentDataBlock

  ! first remove block register


  do while(associated(srd%nodeBlockChain%first))
   deallocate(srd%nodeBlockChain%first%block,stat=allocateStatus)       
   if (allocateStatus /= 0) STOP "ERROR: Could not dellocate srd%nodeBlockChain%first%block in prepareDataStructure"
   srd%nodeBlockChain%first=>srd%nodeBlockChain%first%next      
  end do

  deallocate(srd%nodeBlockChain,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Could not dellocate srd%nodeBlockChain in prepareDataStructure"
  nullify(srd%nodeBlockChain)
  deallocate(srd%nodeStartRegister,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Could not dellocate srd%nodeStartRegister in prepareDataStructure"
  nullify(srd%nodeStartRegister)

  allocate(srd%nodeSideArrayRegister(srd%numberOfNodes+1),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory in prepareDataStructure"
  allocate(srd%nodeSideArray(2*srd%numberOfSides),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory in prepareDataStructure"

  registerSize = 2*srd%numberOfSides

  ! first count

  counter = 0
  srd%nodeSideArrayRegister = 0
  currentDataBlock => srd%dataBlockChain%first
  do while(counter<srd%numberOfSides)
   do i=1,srd%dataBlockSize
    counter = counter + 1
    srd%nodeSideArrayRegister(currentDataBlock%block(i)%sideIndexes(1)) =&
     srd%nodeSideArrayRegister(currentDataBlock%block(i)%sideIndexes(1)) + 1
    srd%nodeSideArrayRegister(currentDataBlock%block(i)%sideIndexes(2)) =&
     srd%nodeSideArrayRegister(currentDataBlock%block(i)%sideIndexes(2)) + 1
    if(counter.ge.srd%numberOfSides) exit
   end do
   currentDataBlock => currentDataBlock%next
  end do

 ! set start indices

 counter = 1
 do i=1,srd%numberOfNodes
  value = srd%nodeSideArrayRegister(i)
  srd%nodeSideArrayRegister(i) = counter
  counter = counter + value
 end do
 srd%nodeSideArrayRegister(srd%numberOfNodes+1) = counter

 ! transfer data

  srd%nodeSideArray = 0
  counter = 0
  currentDataBlock => srd%dataBlockChain%first
  do while(counter<srd%numberOfSides)
   do i=1,srd%dataBlockSize
    counter = counter + 1
    ind1 = currentDataBlock%block(i)%sideIndexes(1)
    ind2 = currentDataBlock%block(i)%sideIndexes(2)
    do j=srd%nodeSideArrayRegister(ind1),srd%nodeSideArrayRegister(ind1+1)-1
     if(srd%nodeSideArray(j)==0) then
      srd%nodeSideArray(j) = counter
      ind1 = 0
      exit
     end if
    end do
    do j=srd%nodeSideArrayRegister(ind2),srd%nodeSideArrayRegister(ind2+1)-1
     if(srd%nodeSideArray(j)==0) then
      srd%nodeSideArray(j) = counter
      ind2 = 0
      exit
     end if
    end do
    if(ind1>0.or.ind2>0) STOP "ERROR: Something's wrong in prepareDataStructure - 1"
    if(counter.ge.srd%numberOfSides) exit
   end do
   currentDataBlock => currentDataBlock%next
  end do

 ! sideData now contains information on which sides that are connected to a given node

 end subroutine prepareDataStructure
!-------------------------------------------------------------------------
 subroutine meshPartition(numberOfPartitions,srd)
! colors nodes for parallelization
 IMPLICIT NONE

 integer :: numberOfPartitions
 type(SideRegisterData) :: srd


  integer,pointer :: frontToNodeList(:),nodeToFrontList(:),numberOfLocalConnectivities(:)

  integer :: i,j,k,l,currentInd,currentPartition,numberOfSidesInPartition,ind2
  integer :: allocateStatus,sidesPerPartition,currentSide,currentConnectivity,searchConnectivity,shiftStop
  integer :: numberOfNodesInFront,numberOfNodesInPartition,numberOfInternalPartitionSides

  type(SideData) :: sd

  allocate(srd%parallelizationColorArray(srd%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory in meshPartition"
  allocate(frontToNodeList(srd%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory in meshPartition"
  allocate(numberOfLocalConnectivities(srd%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory in meshPartition"
  allocate(nodeToFrontList(srd%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory in meshPartition"

  write(*,*) "Partitioning mesh..."

   sidesPerPartition = (numberOfPartitions+srd%numberOfSides)/numberOfPartitions

   currentPartition = 1
   currentInd = 1
   numberOfInternalPartitionSides = 0
   srd%parallelizationColorArray = 0
   nodeToFrontList = 0
   frontToNodeList = 0
   nodeToFrontList(currentInd) = 1
   frontToNodeList(1) = currentInd
   numberOfSidesInPartition = 0 
   numberOfNodesInPartition = 0
   numberOfLocalConnectivities = 0
   numberOfNodesInFront = 1

   do 
    numberOfNodesInPartition = numberOfNodesInPartition + 1

    ! remove currentInd from front list
    do k=1,numberOfNodesInFront-1
     frontToNodeList(k) = frontToNodeList(k+1) 
     nodeToFrontList(frontToNodeList(k)) = k
    end do
    numberOfNodesInFront = numberOfNodesInFront - 1
    frontToNodeList(numberOfNodesInFront+1) = 0
    srd%parallelizationColorArray(currentInd) = currentPartition
 
    do j=srd%nodeSideArrayRegister(currentInd),srd%nodeSideArrayRegister(currentInd+1)-1
     currentSide = srd%nodeSideArray(j) 
     sd = getSideNumber(currentSide,srd)
     if(sd%sideIndexes(1)==currentInd) then 
      ind2 = sd%sideIndexes(2)
     else if(sd%sideIndexes(2)==currentInd) then
      ind2 = sd%sideIndexes(1)
     else  
      STOP "ERROR: Something's gone wrong in meshPartition - 1" 
     end if
     if(srd%parallelizationColorArray(ind2)==srd%parallelizationColorArray(currentInd)) then 
       numberOfInternalPartitionSides = numberOfInternalPartitionSides + 1 
     else
      if(srd%parallelizationColorArray(ind2)==0) then 
       numberOfSidesInPartition = numberOfSidesInPartition + 1
       numberOfLocalConnectivities(ind2) = numberOfLocalConnectivities(ind2) + 1  
       
       ! place new node in front list ranked according to number of external connectivities
      
       currentConnectivity = srd%nodeSideArrayRegister(ind2+1)-srd%nodeSideArrayRegister(ind2)
       currentConnectivity = currentConnectivity - numberOfLocalConnectivities(ind2) ! remove connectivities in current partition
       if(nodeToFrontList(ind2)==0) then 
        ! node is not already in front
        numberOfNodesInFront = numberOfNodesInFront + 1
        frontToNodeList(numberOfNodesInFront) = 0
        shiftStop = numberOfNodesInFront

        do k=1,numberOfNodesInFront
         if(frontToNodeList(k)>0) then
          searchConnectivity = srd%nodeSideArrayRegister(frontToNodeList(k)+1)-srd%nodeSideArrayRegister(frontToNodeList(k))
          searchConnectivity = searchConnectivity - numberOfLocalConnectivities(frontToNodeList(k)) ! remove connectivities in current partition
         else
          searchConnectivity = currentConnectivity + 1
         end if
         if(currentConnectivity<searchConnectivity) then
          ! shift list down
          do l=numberOfNodesInFront,k+1,-1
           if(l>1) then
            frontToNodeList(l) = frontToNodeList(l-1)
            nodeToFrontList(frontToNodeList(l)) = l
           end if
          end do
          frontToNodeList(k) = ind2
          nodeToFrontList(ind2) = k
          exit
         end if
        end do
       else 
        shiftStop = nodeToFrontList(ind2)  
        do k=1,shiftStop-1
         searchConnectivity = srd%nodeSideArrayRegister(frontToNodeList(k)+1)-srd%nodeSideArrayRegister(frontToNodeList(k))
         searchConnectivity = searchConnectivity - numberOfLocalConnectivities(frontToNodeList(k)) ! remove connectivities in current partition
         if(currentConnectivity<searchConnectivity) then
         ! shift list down
          do l=shiftStop,k+1,-1
           if(l>1) then
            frontToNodeList(l) = frontToNodeList(l-1)
            nodeToFrontList(frontToNodeList(l)) = l
           end if
          end do
          frontToNodeList(k) = ind2
          nodeToFrontList(ind2) = k
          exit
         end if
        end do
       end if
      end if
     end if
    end do 

    if(numberOfSidesInPartition.ge.sidesPerPartition.or.frontToNodeList(1)==0) then 
     write(*,'(A,I3,A,I8,A,I8,A,I8)') " Number of sides (internal) and nodes in partition ",currentPartition,": ",&
             numberOfSidesInPartition," (",numberOfInternalPartitionSides,") ",numberOfNodesInPartition
     currentPartition = currentPartition + 1
     numberOfSidesInPartition = 0 
     numberOfNodesInPartition = 0 
     numberOfInternalPartitionSides = 0
     numberOfLocalConnectivities = 0         
    end if 

    currentInd = frontToNodeList(1) 
    if(currentInd==0) then 
     ! front has stopped
     do i=1,srd%numberOfNodes
      if(srd%parallelizationColorArray(i)==0) then 
       write(*,'(A,I8,A)') "ERROR: Node ",i," is not connected to a region" 
       STOP
      end if
     end do 
     exit 
    end if
   end do

  deallocate(frontToNodeList,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in meshPartition"
  deallocate(nodeToFrontList,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in meshPartition"
  deallocate(numberOfLocalConnectivities,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in meshPartition"

  write(*,*) "Mesh partition completed"
 end subroutine meshPartition
!-------------------------------------------------------------------------
 subroutine doCacheRemap(srd,oldToNewNodeMapping)
 IMPLICIT NONE
  type(SideRegisterData) :: srd
  integer,pointer :: oldToNewNodeMapping(:)

  integer,pointer :: newNumbering(:)
  integer,pointer :: frontToNodeList(:)
  integer,pointer :: numberOfLocalConnectivities(:)
  integer,pointer :: nodeToFrontList(:)

  integer :: numberOfNodesInPartition,currentNodeNumber,numberOfNodesInFront

  type(SideData) :: sd

  integer :: i,j,k,l,allocateStatus,shiftStop,currentInd,ind2,currentSide,currentConnectivity
  integer :: searchConnectivity

  write(*,*) "Doing cache remap..."
 
! count number of nodes in domain 

  numberOfNodesInPartition = 0 
  do i=1,srd%numberOfNodes
   if(oldToNewNodeMapping(i)>0) then 
    numberOfNodesInPartition = numberOfNodesInPartition + 1 
   end if
  end do   


  allocate(newNumbering(srd%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory in meshPartition"
  allocate(frontToNodeList(srd%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory in meshPartition"
  allocate(numberOfLocalConnectivities(srd%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory in meshPartition"
  allocate(nodeToFrontList(srd%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory in meshPartition"

   currentNodeNumber = 0  
   do i=1,srd%numberOfNodes
    if(oldToNewNodeMapping(i)>0) then 
     currentInd = i
    end if
   end do
   newNumbering = 0
   nodeToFrontList = 0
   frontToNodeList = 0
   nodeToFrontList(currentInd) = 1
   frontToNodeList(1) = currentInd
   numberOfLocalConnectivities = 0
   numberOfNodesInFront = 1

   do
    currentNodeNumber = currentNodeNumber + 1
 
    ! remove currentInd from front list
    do k=1,numberOfNodesInFront-1
     frontToNodeList(k) = frontToNodeList(k+1)
     nodeToFrontList(frontToNodeList(k)) = k
    end do
    numberOfNodesInFront = numberOfNodesInFront - 1
    frontToNodeList(numberOfNodesInFront+1) = 0

    newNumbering(currentInd) = currentNodeNumber

    do j=srd%nodeSideArrayRegister(currentInd),srd%nodeSideArrayRegister(currentInd+1)-1
     currentSide = srd%nodeSideArray(j)
     sd = getSideNumber(currentSide,srd)
     if(sd%sideIndexes(1)==currentInd) then
      ind2 = sd%sideIndexes(2)
     else if(sd%sideIndexes(2)==currentInd) then
      ind2 = sd%sideIndexes(1)
     else 
      STOP "ERROR: Something's gone wrong in meshPartition - 1"
     end if
     if(oldToNewNodeMapping(ind2)>0) then 
      if(newNumbering(ind2)==0) then
       numberOfLocalConnectivities(ind2) = numberOfLocalConnectivities(ind2) + 1

       ! place new node in front list ranked according to number of external connectivities
 
       currentConnectivity = srd%nodeSideArrayRegister(ind2+1)-srd%nodeSideArrayRegister(ind2)
       currentConnectivity = currentConnectivity - numberOfLocalConnectivities(ind2) 
       if(nodeToFrontList(ind2)==0) then
        ! node is not already in front
        numberOfNodesInFront = numberOfNodesInFront + 1
        frontToNodeList(numberOfNodesInFront) = 0
        shiftStop = numberOfNodesInFront

        do k=1,numberOfNodesInFront
         if(frontToNodeList(k)>0) then
          searchConnectivity = srd%nodeSideArrayRegister(frontToNodeList(k)+1)-srd%nodeSideArrayRegister(frontToNodeList(k))
          searchConnectivity = searchConnectivity - numberOfLocalConnectivities(frontToNodeList(k)) 
         else
          searchConnectivity = currentConnectivity + 1
         end if
         if(currentConnectivity<searchConnectivity) then
          ! shift list down
          do l=numberOfNodesInFront,k+1,-1
           if(l>1) then
            frontToNodeList(l) = frontToNodeList(l-1)
            nodeToFrontList(frontToNodeList(l)) = l
           end if
          end do
          frontToNodeList(k) = ind2
          nodeToFrontList(ind2) = k
          exit
         end if
        end do
       else
        shiftStop = nodeToFrontList(ind2)
        do k=1,shiftStop-1
         searchConnectivity = srd%nodeSideArrayRegister(frontToNodeList(k)+1)-srd%nodeSideArrayRegister(frontToNodeList(k))
         searchConnectivity = searchConnectivity - numberOfLocalConnectivities(frontToNodeList(k)) 
         if(currentConnectivity<searchConnectivity) then
         ! shift list down
          do l=shiftStop,k+1,-1
           if(l>1) then
            frontToNodeList(l) = frontToNodeList(l-1)
            nodeToFrontList(frontToNodeList(l)) = l
           end if
          end do
          frontToNodeList(k) = ind2
          nodeToFrontList(ind2) = k
          exit
         end if
        end do
       end if
      end if
     end if
    end do

    currentInd = frontToNodeList(1)
    if(currentInd==0) then
     ! front has stopped
     do i=1,srd%numberOfNodes
      if((newNumbering(i)==0.and.oldToNewNodeMapping(i).ne.0)) then
       currentInd = i 
       frontToNodeList(1) = currentInd
       nodeToFrontList(currentInd) = 1
       numberOfNodesInFront = 1
       exit
      end if
     end do
     if(currentInd==0) then
      exit
     end if
    end if
   end do

   do i=1,srd%numberOfNodes
    if((newNumbering(i)==0.and.oldToNewNodeMapping(i).ne.0).or.(newNumbering(i).ne.0.and.oldToNewNodeMapping(i)==0)) then
     write(*,*) "ERROR: Renumbering failed for node ",i,newNumbering(i),oldToNewNodeMapping(i)
     STOP
    end if
    oldToNewNodeMapping(i) = newNumbering(i)
   end do

  write(*,*) "Number of nodes in partition: ",currentNodeNumber,numberOfNodesInPartition

  deallocate(newNumbering,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in meshPartition"
  deallocate(frontToNodeList,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in meshPartition"
  deallocate(nodeToFrontList,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in meshPartition"
  deallocate(numberOfLocalConnectivities,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in meshPartition"

  write(*,*) "Cache remap completed" 
  
 end subroutine doCacheRemap
!-------------------------------------------------------------------------
 subroutine greedyMeshPartition(npart,srd)
! colors nodes for parallelization using greedy type alogorithm
 IMPLICIT NONE

 integer :: npart
 type(SideRegisterData) :: srd

  integer :: allocateStatus,nside,npoin,nsm,ncl,nsc,nei,ip,is,i,iss,il,i1,i2,ist,in
  integer :: isp,sideNumber,ij

  integer,pointer :: lph(:),list(:),lct(:),lbf(:),lpw(:),icone(:),iside(:,:),lcp(:)

  type(SideDataBlockInstance),pointer :: currentSideBlock

! make iside array

  allocate(iside(srd%numberOfSides,2),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory in greedyMeshPartition"

  sideNumber = 0
  currentSideBlock => srd%DataBlockChain%first
  do while(associated(currentSideBlock))
   do i=1,srd%dataBlockSize
    sideNumber = sideNumber + 1
    iside(sideNumber,1) = currentSideBlock%block(i)%sideIndexes(1)
    iside(sideNumber,2) = currentSideBlock%block(i)%sideIndexes(2)
    if(sideNumber==srd%numberOfSides) exit
   end do
   currentSideBlock => currentSideBlock%next
  end do

  allocate(srd%parallelizationColorArray(srd%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory in greedyMeshPartition"

! make help arrays

  allocate(lph(srd%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory in greedyMeshPartition"
  allocate(list(srd%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory in greedyMeshPartition"
  allocate(lct(srd%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory in greedyMeshPartition"
  allocate(lbf(srd%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory in greedyMeshPartition"
  allocate(lpw(srd%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory in greedyMeshPartition"

  allocate(icone(2*srd%numberOfSides),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory in greedyMeshPartition"

  nside = srd%numberOfSides
  npoin = srd%numberOfNodes

  lcp => srd%parallelizationColorArray

   nsm = (nside+npart)/npart
   ncl = 1
   nsc = 0
   nei = 0

   do ip = 1 , npoin
    lph(ip) = 0
    lcp(ip) = 0
    lct(ip) = 0
    list(ip)= 0
    lbf(ip) = 0
   end do

   do is = 1 , nside
    lph(iside(is,1)) = lph(iside(is,1)) + 1
    lph(iside(is,2)) = lph(iside(is,2)) + 1
   end do
   lpw(1)     = 1
   iss        = 1
   do i = 2 , npoin
    if(lph(i).lt.lph(iss)) iss = i
    lpw(i)    = lpw(i-1) + lph(i-1)
   end do

   do i = 1 , npoin
    lph(i) = 0
   end do

   do is = 1 , nside
    i1        = iside(is,1)
    lph(i1)   = lph(i1) + 1
    il        = lpw(i1) + lph(i1) -1
    icone(il) = is
    i2        = iside(is,2)
    lph(i2)   = lph(i2) + 1
    il        = lpw(i2) + lph(i2) -1
    icone(il) = is
   end do

   ist        = iss

10 lcp(ist) = ncl
   do ij = lpw(ist) , lpw(ist)+lph(ist)-1
    is = icone(ij)
    i2 = iside(is,2) + iside(is,1) - ist
    if(lcp(ist).eq.lcp(i2)) then
     nei       = nei + 1
    else
     lct(i2) = lct(i2) + 1
     if(lcp(i2).eq.0) then
! i2 not touched yet

      nsc       = nsc + 1
      if(lbf(i2).ne.0) then
       list(lbf(i2)) = 0
      end if
      isp= ist
12    in = list(isp)
      if(in.eq.0) then
       list(isp) = i2
       lbf( i2 ) = isp
      else if((lph(i2)-2*lct(i2)).lt.(lph(in)-2*lct(in))) then
       list (isp) = i2
       lbf ( i2 ) = isp
       isp        = i2
       i2         = in
       goto 12
      else
       isp        = in
       goto 12
      end if
     end if
    end if
   end do

   ist = list(ist)
   if(nsc.ge.nsm.or.ist.eq.0) then
    print *, 'N. of side for part: ', ncl,nsc,nei
    ncl = ncl + 1
    nsc = 0
    nei = 0
    do i = 1 , npoin
     lct(i) = 0
    end do
   end if


   if(ist.ne.0) then
    goto 10
   else
    do ip = 1 , npoin
     if(lcp(ip).eq.0) print *,ip,' has no region'
     lcp(ip) = lcp(ip) -  1
    end do
   end if

  deallocate(lph,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in greedyMeshPartition"
  deallocate(list,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in greedyMeshPartition"
  deallocate(lct,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in greedyMeshPartition"
  deallocate(lbf,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in greedyMeshPartition"
  deallocate(lpw,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in greedyMeshPartition"

  deallocate(icone,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in greedyMeshPartition"

 end subroutine greedyMeshPartition
!-------------------------------------------------------------------------
 subroutine meshPartition2(numberOfPartitions,srd)
! colors nodes for parallelization using MeTiS
 IMPLICIT NONE
 integer :: numberOfPartitions
 type(SideRegisterData) :: srd


  real*4 ,pointer :: partitionWeighting(:)
  integer :: options(5),numberOfCutEdges,allocateStatus

  integer,pointer :: connectivityRegister(:),connectivityArray(:)
  integer,pointer :: nodesInPartition(:)

  integer :: j,sideNumber,ind1,ind2
  type(SideDataBlockInstance),pointer :: currentSideBlock

  ! consider deallocating srd%nodeSideArray

  ! set up connectivity array

  if(numberOfPartitions>1) then 
   allocate(connectivityRegister(srd%numberOfNodes),stat=allocateStatus)
   if (allocateStatus /= 0) STOP "ERROR: Not enough memory in meshPartition2"
   allocate(connectivityArray(2*srd%numberOfSides),stat=allocateStatus)
   if (allocateStatus /= 0) STOP "ERROR: Not enough memory in meshPartition2"

   connectivityRegister = 0
 
   sideNumber = 0
   currentSideBlock => srd%DataBlockChain%first
   do while(associated(currentSideBlock))
    do j=1,srd%dataBlockSize
     sideNumber = sideNumber+1
     ind1 = currentSideBlock%block(j)%sideIndexes(1)
     ind2 = currentSideBlock%block(j)%sideIndexes(2)
     ! place in array
     connectivityArray(srd%nodeSideArrayRegister(ind1)+connectivityRegister(ind1))=ind2
     connectivityArray(srd%nodeSideArrayRegister(ind2)+connectivityRegister(ind2))=ind1
     connectivityRegister(ind1) = connectivityRegister(ind1)+1
     connectivityRegister(ind2) = connectivityRegister(ind2)+1
     if(sideNumber==srd%numberOfSides) exit
    end do
    currentSideBlock => currentSideBlock%next
   end do
  end if

  allocate(srd%parallelizationColorArray(srd%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory in meshPartition2"

  if(numberOfPartitions>1) then
 
   allocate(partitionWeighting(numberOfPartitions),stat=allocateStatus)
   if (allocateStatus /= 0) STOP "ERROR: Not enough memory in meshPartition2"
   srd%parallelizationColorArray = 0

   partitionWeighting = 1.0/float(numberOfPartitions)
 
   options = 0

   write(*,*) "calling metis..."
 
   call METIS_WPartGraphKway(srd%numberOfNodes,srd%nodeSideArrayRegister,connectivityArray,0,0,0,1,&
                             numberOfPartitions,partitionWeighting,options,numberOfCutEdges,&
                             srd%parallelizationColorArray)

   write(*,*) "metis finished, number of cut sides: ",numberOfCutEdges

   deallocate(connectivityRegister,stat=allocateStatus)
   if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in meshPartition2"
   deallocate(connectivityArray,stat=allocateStatus)
   if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in meshPartition2"
   deallocate(partitionWeighting,stat=allocateStatus)
   if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in meshPartition2"
  else
   srd%parallelizationColorArray = 1
  end if


  allocate(nodesInPartition(numberOfPartitions),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory in meshPartition2"

  nodesInPartition = 0

  do j=1,srd%numberOfNodes
   nodesInPartition(srd%parallelizationColorArray(j)) = nodesInPartition(srd%parallelizationColorArray(j)) + 1
  end do

  do j=1,numberOfPartitions
   write(*,*) "Number of nodes for partition ",j,": ",nodesInPartition(j)
  end do


  deallocate(nodesInPartition,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in meshPartition2"
 end subroutine meshPartition2
!-------------------------------------------------------------------------
end module SideRegister



