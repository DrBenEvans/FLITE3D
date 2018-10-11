!**********************************************
!*Aggl3d_Grid.f90                             *
!*                                            *
!*Grid module definition file                 *
!*                                            * 
!*   Aggl3d v 2.3.2                           * 
!*                                            *
!*24.06.99-                                   *
!**********************************************
!*Description:                                *
!* File contains grid datastructure and       *
!* agglomeration subroutines                  *
!**********************************************

module GridInput
 type GridInputData
  real :: tripRadius
  real :: wallDistanceThickness
  integer :: numberOfTripLines
  real :: tripLineCoordinates(100,6)
 end type GridInputData

contains

 subroutine readGridInput(gid,INFILE)
  IMPLICIT NONE
 
  type(GridInputData) gid
  integer :: INFILE

  namelist /InputVariables/ gid

  gid%tripRadius = 1.0
  gid%wallDistanceThickness = 0.0
  gid%numberOfTripLines = 0
  gid%tripLineCoordinates = 0.0

  read(INFILE,InputVariables)

 end subroutine readGridInput
end module GridInput 

module Grid
! holds variables and procedures for a single grid
 use SideRegister
 use CoordinateRegister
 use TetrahedralElements
 use PyramidElements
 use PrismElements
 use HexahedralElements
 use BoundaryRegister
 use VisualizationModule
 use ToolBox
 use GridInput

 type GridData
  type(BoundaryRegisterData),pointer :: brp ! holds boundary data in grid
  type(SideRegisterData),pointer :: sdp
!  type(VisualizationData),pointer :: vdp

  type(GridInputData),pointer :: gid

  real,pointer :: controlVolumeNodeSize(:) ! contains area of control volumes
                                           ! to be used for intergrid mappings

  real,pointer :: prevControlVolume(:),prevControlVolume2(:)

  real,pointer :: wallDistance(:)
  integer,pointer :: wallDistanceBoundaryNodeArray(:)

  integer,pointer :: nodeSwitchRegister(:)
  integer,pointer :: fullOldToNewNodeMapping(:)
  integer,pointer :: fineToCoarseIndexMapping(:)
  integer,pointer :: seedNodeRegister(:)

  integer,pointer :: numberOfIntNodesInPart(:)

  integer :: numberOfFineComNodes
  integer,pointer :: fineComNodes(:)

  integer :: numberOfTripNodes
  integer :: numberOfTripFieldNodes
  integer,pointer :: tripLineIndexes(:)
  real,pointer :: tripNodeFieldDistances(:)
  integer,pointer :: tripNodeFieldIndexes(:,:)
  real,pointer :: tripWallLength(:)

  integer :: numberOfSides,numberOfNodes
  real :: directionalityParameter,minimumAspectRatio 
  integer :: gridNumber,numberOfGridDomains

  
  logical :: isHybrid
  logical :: hasSwitched
 end type GridData

 contains

!-------------------------------------------------------------------------
 subroutine cleanUpGridsAfterAgglomeration(finegrp,coarsegrp)
 ! deallocate arrays in grid
 
 IMPLICIT NONE

 type(GridData) :: finegrp,coarsegrp

 type(linkedInteger),pointer :: currentInteger
 integer :: allocateStatus,i

! deallocate(finegrp%controlVolumeNodeSize,stat=allocateStatus)
! if (allocateStatus /= 0) STOP "ERROR: cleanUpGridAfterAgglomeration couldn't deallocate"

 allocate(coarsegrp%brp%boundaryFaceNodeMappings(finegrp%brp%numberOfBoundaryNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: cleanUpGridAfterAgglomeration couldn't deallocate"

 do i=1,finegrp%brp%numberOfBoundaryNodes
  coarsegrp%brp%boundaryFaceNodeMappings(i) = coarsegrp%fineToCoarseIndexMapping(i)
 end do

 if(associated(finegrp%fineToCoarseIndexMapping)) then 
  deallocate(finegrp%fineToCoarseIndexMapping,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: cleanUpGridAfterAgglomeration couldn't deallocate"
 end if

 end subroutine cleanUpGridsAfterAgglomeration
!-------------------------------------------------------------------------
 subroutine setBoundaryFirst(grp,crp,crpPrev,crpPrev2,hep,pyp,prp,tep,brp)
 ! subroutine to make sure that the boundary nodes are numbered first

 IMPLICIT NONE

 type(GridData) :: grp
 type(CoordinateRegisterData) :: crp,crpPrev,crpPrev2
 type(HexahedralElementData) :: hep
 type(PyramidElementData) :: pyp
 type(PrismElementData) :: prp
 type(TetrahedralElementData) :: tep
 type(BoundaryRegisterData) :: brp

 integer :: i,jstart,kstart,j,k,ind1,switchNode,l,allocateStatus,searchNode
 real :: coor(3)
 logical,pointer :: nodeOK(:)
 logical :: notOnBoundary,foundOne
 integer,pointer :: swapList(:,:)
 integer :: numberOfSwitches,ind,count1,count2,buff
 logical :: itWentAsPlanned

 type(STEData),pointer :: tetTemp
 type(SPRData),pointer :: priTemp
 type(SPYData),pointer :: pyrTemp
 type(SHEData),pointer :: hexTemp

 integer :: checkind

 write(*,*) "Reordering..."

 allocate(nodeOK(grp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: setBoundaryFirst out of memory"


 nodeOK =.false.

 do i=1,brp%numberOfTriangularBFs
  do j=1,3
   ind = brp%triangularFaces(i)%faceIndexes(j)
   nodeOK(ind) = .true.
  end do
 end do


 do i=1,brp%numberOfRectangularBFs
  do j=1,4
   ind = brp%rectangularFaces(i)%faceIndexes(j)
   nodeOK(ind) = .true.
  end do
 end do

 brp%numberOfBoundaryNodes = 0
 do i=1,grp%numberOfNodes
  if(nodeOK(i)) brp%numberOfBoundaryNodes = brp%numberOfBoundaryNodes + 1
 end do
 brp%numberOfBoundaryNormals = 0
 write(*,*) "Number of boundary nodes: ",brp%numberOfBoundaryNodes

 allocate(swapList(brp%numberOfBoundaryNodes,2),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: setBoundaryFirst out of memory"

 nodeOK = .true.
 nodeOK(1:brp%numberOfBoundaryNodes) = .false.
 count1 = 0
 count2 = 0
 do i=1,brp%numberOfTriangularBFs
  do j=1,3
   ind = brp%triangularFaces(i)%faceIndexes(j)
   if(ind.le.brp%numberOfBoundaryNodes) then
    if(.not.nodeOK(ind)) count2 = count2 + 1
    nodeOK(ind) = .true.
   else
    if(nodeOK(ind)) count1 = count1 + 1
    nodeOK(ind) = .false.
   end if
  end do
 end do

 count1 = 0
 count2 = 0
 do i=1,brp%numberOfRectangularBFs
  do j=1,4
   ind = brp%rectangularFaces(i)%faceIndexes(j)
   if(ind.le.brp%numberOfBoundaryNodes) then
    if(.not.nodeOK(ind)) count2 = count2 + 1
    nodeOK(ind) = .true.
   else
    if(nodeOK(ind)) count1 = count1 + 1
    nodeOK(ind) = .false.
   end if
  end do
 end do

 ! now look for nodes that are larger

 numberOfSwitches = 0

 do i=1,brp%numberOfTriangularBFs
  do j=1,3
   ind = brp%triangularFaces(i)%faceIndexes(j)
   if(.not.nodeOK(ind)) then
    ! find node to switch with
    itWentAsPlanned = .false.
    do k=1,brp%numberOfBoundaryNodes
     if(.not.nodeOK(k)) then
      numberOfSwitches = numberOfSwitches + 1
      swapList(numberOfSwitches,1) = ind
      swapList(numberOfSwitches,2) = k
      nodeOK(k) = .true.
      nodeOK(ind) = .true.
      itWentAsPlanned = .true.
      exit
     end if
    end do
    if(.not.itWentAsPlanned) write(*,*) "Switches: ",numberOfSwitches
    if(.not.itWentAsPlanned) STOP "ERROR: Something's gone wrong in setBoundaryFirst - 1"
   end if
  end do
 end do

do i=1,brp%numberOfRectangularBFs
  do j=1,4
   ind = brp%rectangularFaces(i)%faceIndexes(j)
   if(.not.nodeOK(ind)) then
    ! find node to switch with
    itWentAsPlanned = .false.
    do k=1,brp%numberOfBoundaryNodes
     if(.not.nodeOK(k)) then
      numberOfSwitches = numberOfSwitches + 1
      swapList(numberOfSwitches,1) = ind
      swapList(numberOfSwitches,2) = k
      nodeOK(k) = .true.
      nodeOK(ind) = .true.
      itWentAsPlanned = .true.
      exit
     end if
    end do
    if(.not.itWentAsPlanned) then
     write(*,*) "Switches: ",numberOfSwitches,i,j,k,ind
     do k=1,grp%numberOfNodes
      if(.not.nodeOK(k)) write(*,*) "morna: ",k,nodeOK(k)
     end do
    end if

    if(.not.itWentAsPlanned) STOP "ERROR: Something's gone wrong in setBoundaryFirst - 2"
   end if
  end do
 end do

 deallocate(nodeOK,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: setBoundaryFirst couldn't deallocate"


! switch coordinates

 do i=1,numberOfSwitches
  coor = getCoor(swapList(i,2),crp)
  call setCoor(swapList(i,2),getCoor(swapList(i,1),crp),crp)
  call setCoor(swapList(i,1),coor,crp)
  coor = getCoor(swapList(i,2),crpPrev)
  call setCoor(swapList(i,2),getCoor(swapList(i,1),crpPrev),crpPrev)
  call setCoor(swapList(i,1),coor,crpPrev)
  coor = getCoor(swapList(i,2),crpPrev2)
  call setCoor(swapList(i,2),getCoor(swapList(i,1),crpPrev2),crpPrev2)
  call setCoor(swapList(i,1),coor,crpPrev2)
 end do


! do i=1,grp%numberOfNodes
!  write(557,*) i,getCoor(i,crp)
! end do


 ! to speed things up, make a full node register
 allocate(grp%nodeSwitchRegister(grp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: setBoundaryFirst out of memory"

 do i=1,grp%numberOfNodes
  grp%nodeSwitchRegister(i) = i
 end do

 do i=1,numberOfSwitches
  buff = grp%nodeSwitchRegister(swapList(i,1))
  grp%nodeSwitchRegister(swapList(i,1)) = grp%nodeSwitchRegister(swapList(i,2))
  grp%nodeSwitchRegister(swapList(i,2)) = buff
 end do

! switch element indexes

 do j=1,hep%numberOfElements
  hexTemp => getHexElement(j,hep)
  do k=1,8
   hexTemp%pointIndexes(k) = grp%nodeSwitchRegister(hexTemp%pointIndexes(k))
  end do
 end do

 do j=1,prp%numberOfElements
  priTemp => getPriElement(j,prp)
  do k=1,6
   priTemp%pointIndexes(k) = grp%nodeSwitchRegister(priTemp%pointIndexes(k))
  end do
 end do

 do j=1,pyp%numberOfElements
  pyrTemp => getPyrElement(j,pyp)
  do k=1,5
   pyrTemp%pointIndexes(k) = grp%nodeSwitchRegister(pyrTemp%pointIndexes(k))
  end do
 end do

 do j=1,tep%numberOfElements
  tetTemp => getTetElement(j,tep)
  if(j==546991) write(*,*) "TTR: ",tetTemp%pointIndexes
  do k=1,4
   tetTemp%pointIndexes(k) = grp%nodeSwitchRegister(tetTemp%pointIndexes(k))
  end do
 end do

! finally do switching in boundary register

 do j=1,brp%numberOfTriangularBFs
  do k=1,3
   brp%triangularFaces(j)%faceIndexes(k) = grp%nodeSwitchRegister(brp%triangularFaces(j)%faceIndexes(k))
   if(brp%triangularFaces(j)%faceIndexes(k)>brp%numberOfBoundaryNodes.or.&
      brp%triangularFaces(j)%faceIndexes(k)<1) STOP "ERROR: Something's wrong in setBoundaryFirst - 3"
  end do
 end do
 do j=1,brp%numberOfRectangularBFs
  do k=1,4
   brp%rectangularFaces(j)%faceIndexes(k) = grp%nodeSwitchRegister(brp%rectangularFaces(j)%faceIndexes(k))
   if(brp%rectangularFaces(j)%faceIndexes(k)>brp%numberOfBoundaryNodes.or.&
      brp%rectangularFaces(j)%faceIndexes(k)<1) STOP "ERROR: Something's wrong in setBoundaryFirst - 4"
  end do
 end do

 ! write switch list

 if(numberOfSwitches>0) then
  grp%hasSwitched = .true.
 else
  grp%hasSwitched = .false.
 end if

 write(*,*) "Writing switch file..."
 open(21,file="switch.reg",form='formatted',status='unknown')
 write(21,*) numberOfSwitches
 do i=1,numberOfSwitches
  write(21,*) swapList(i,:)
 end do
 close(21)
 deallocate(swapList,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: setBoundaryFirst couldn't deallocate"
 
! pyrTemp => getPyrElement(2700,pyp)
! write(*,*) "K: ",2700,pyrTemp%pointIndexes

 end subroutine setBoundaryFirst
!-------------------------------------------------------------------------
 subroutine makeNewPltFile(grp,crp,hep,pyp,prp,tep,brp,OUTFILE)
 IMPLICIT NONE

 type(GridData) :: grp
 type(CoordinateRegisterData) :: crp
 type(HexahedralElementData) :: hep
 type(PyramidElementData) :: pyp
 type(PrismElementData) :: prp
 type(TetrahedralElementData) :: tep
 type(BoundaryRegisterData) :: brp

 integer,pointer :: elementRegister(:,:),faceRegister(:,:),segmentRegister(:)
 integer,pointer :: elementContainingFace(:)

 integer :: numberOfTetEl,numberOfTriBFs,p1,p2,p3,p4,p5,p6,p7,p8
 integer :: index,i,j,q1,q2,q3,q4,q5,q6,q7,q8,allocateStatus
 integer :: OUTFILE

 write(*,*) "Writing base.plt file..."

 numberOfTetEl = tep%numberOfElements + 2*pyp%numberOfElements + 3*prp%numberOfElements + 6*hep%numberOfElements
 numberOfTriBFs = brp%numberOfTriangularBFs + 2*brp%numberOfRectangularBFs

 allocate(elementRegister(numberOfTetEl,4),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: makeNewPltFile out of memory"
 allocate(faceRegister(numberOfTriBFs,3),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: makeNewPltFile out of memory"
 allocate(segmentRegister(numberOfTriBFs),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: makeNewPltFile out of memory"
 allocate(elementContainingFace(numberOfTriBFs),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: makeNewPltFile out of memory"

 do i=1,tep%numberOfElements
  elementRegister(i,:) = tep%elements(i)%pointIndexes(:)
 end do

 index = tep%numberOfElements

 do i=1,pyp%numberOfElements
  index = index + 1
  elementRegister(index,1) = pyp%elements(i)%pointIndexes(1)
  elementRegister(index,2) = pyp%elements(i)%pointIndexes(2)
  elementRegister(index,3) = pyp%elements(i)%pointIndexes(3)
  elementRegister(index,4) = pyp%elements(i)%pointIndexes(5)
  index = index + 1
  elementRegister(index,1) = pyp%elements(i)%pointIndexes(1)
  elementRegister(index,2) = pyp%elements(i)%pointIndexes(3)
  elementRegister(index,3) = pyp%elements(i)%pointIndexes(4)
  elementRegister(index,4) = pyp%elements(i)%pointIndexes(5)
 end do

 do i=1,prp%numberOfElements
  p1=prp%elements(i)%pointIndexes(1)
  p2=prp%elements(i)%pointIndexes(2)
  p3=prp%elements(i)%pointIndexes(3)
  p4=prp%elements(i)%pointIndexes(4)
  p5=prp%elements(i)%pointIndexes(5)
  p6=prp%elements(i)%pointIndexes(6)

  q1=grp%nodeSwitchRegister(p1)
  q2=grp%nodeSwitchRegister(p2)
  q3=grp%nodeSwitchRegister(p3)
  q4=grp%nodeSwitchRegister(p4)
  q5=grp%nodeSwitchRegister(p5)
  q6=grp%nodeSwitchRegister(p6)

  index = index + 1

  if(q4<q5) then
   elementRegister(index,1) = prp%elements(i)%pointIndexes(1)
   elementRegister(index,2) = prp%elements(i)%pointIndexes(2)
   elementRegister(index,3) = prp%elements(i)%pointIndexes(3)
   elementRegister(index,4) = prp%elements(i)%pointIndexes(4)
   if(q5<q6) then
    index = index + 1
    elementRegister(index,1) = prp%elements(i)%pointIndexes(2)
    elementRegister(index,2) = prp%elements(i)%pointIndexes(5)
    elementRegister(index,3) = prp%elements(i)%pointIndexes(3)
    elementRegister(index,4) = prp%elements(i)%pointIndexes(4)
    index = index + 1
    elementRegister(index,1) = prp%elements(i)%pointIndexes(3)
    elementRegister(index,2) = prp%elements(i)%pointIndexes(5)
    elementRegister(index,3) = prp%elements(i)%pointIndexes(6)
    elementRegister(index,4) = prp%elements(i)%pointIndexes(4)
   else
    index = index + 1
    elementRegister(index,1) = prp%elements(i)%pointIndexes(4)
    elementRegister(index,2) = prp%elements(i)%pointIndexes(6)
    elementRegister(index,3) = prp%elements(i)%pointIndexes(5)
    elementRegister(index,4) = prp%elements(i)%pointIndexes(2)
    index = index + 1
    elementRegister(index,1) = prp%elements(i)%pointIndexes(2)
    elementRegister(index,2) = prp%elements(i)%pointIndexes(6)
    elementRegister(index,3) = prp%elements(i)%pointIndexes(3)
    elementRegister(index,4) = prp%elements(i)%pointIndexes(4)
   end if
  else
   elementRegister(index,1) = prp%elements(i)%pointIndexes(1)
   elementRegister(index,2) = prp%elements(i)%pointIndexes(2)
   elementRegister(index,3) = prp%elements(i)%pointIndexes(3)
   elementRegister(index,4) = prp%elements(i)%pointIndexes(5)
   if(q4<q6) then
    index = index + 1
    elementRegister(index,1) = prp%elements(i)%pointIndexes(1)
    elementRegister(index,2) = prp%elements(i)%pointIndexes(3)
    elementRegister(index,3) = prp%elements(i)%pointIndexes(4)
    elementRegister(index,4) = prp%elements(i)%pointIndexes(5)
    index = index + 1
    elementRegister(index,1) = prp%elements(i)%pointIndexes(3)
    elementRegister(index,2) = prp%elements(i)%pointIndexes(6)
    elementRegister(index,3) = prp%elements(i)%pointIndexes(4)
    elementRegister(index,4) = prp%elements(i)%pointIndexes(5)
   else
    index = index + 1
    elementRegister(index,1) = prp%elements(i)%pointIndexes(1)
    elementRegister(index,2) = prp%elements(i)%pointIndexes(5)
    elementRegister(index,3) = prp%elements(i)%pointIndexes(6)
    elementRegister(index,4) = prp%elements(i)%pointIndexes(4)
    index = index + 1
    elementRegister(index,1) = prp%elements(i)%pointIndexes(3)
    elementRegister(index,2) = prp%elements(i)%pointIndexes(5)
    elementRegister(index,3) = prp%elements(i)%pointIndexes(6)
    elementRegister(index,4) = prp%elements(i)%pointIndexes(1)
   end if
  end if
 end do

 do i=1,hep%numberOfElements
  p1=hep%elements(i)%pointIndexes(1)
  p2=hep%elements(i)%pointIndexes(2)
  p3=hep%elements(i)%pointIndexes(3)
  p4=hep%elements(i)%pointIndexes(4)
  p5=hep%elements(i)%pointIndexes(5)
  p6=hep%elements(i)%pointIndexes(6)
  p7=hep%elements(i)%pointIndexes(7)
  p8=hep%elements(i)%pointIndexes(8)

  q1=grp%nodeSwitchRegister(p1)
  q2=grp%nodeSwitchRegister(p2)
  q3=grp%nodeSwitchRegister(p3)
  q4=grp%nodeSwitchRegister(p4)
  q5=grp%nodeSwitchRegister(p5)
  q6=grp%nodeSwitchRegister(p6)
  q7=grp%nodeSwitchRegister(p7)
  q8=grp%nodeSwitchRegister(p8)

  index = index + 1

  if(q5<q6) then
   elementRegister(index,1) = hep%elements(i)%pointIndexes(1)
   elementRegister(index,2) = hep%elements(i)%pointIndexes(2)
   elementRegister(index,3) = hep%elements(i)%pointIndexes(3)
   elementRegister(index,4) = hep%elements(i)%pointIndexes(5)
   if(q6<q7) then
    index = index + 1
    elementRegister(index,1) = hep%elements(i)%pointIndexes(2)
    elementRegister(index,2) = hep%elements(i)%pointIndexes(3)
    elementRegister(index,3) = hep%elements(i)%pointIndexes(5)
    elementRegister(index,4) = hep%elements(i)%pointIndexes(6)
    index = index + 1
    elementRegister(index,1) = hep%elements(i)%pointIndexes(6)
    elementRegister(index,2) = hep%elements(i)%pointIndexes(3)
    elementRegister(index,3) = hep%elements(i)%pointIndexes(5)
    elementRegister(index,4) = hep%elements(i)%pointIndexes(7)
   else
    index = index + 1
    elementRegister(index,1) = hep%elements(i)%pointIndexes(5)
    elementRegister(index,2) = hep%elements(i)%pointIndexes(7)
    elementRegister(index,3) = hep%elements(i)%pointIndexes(6)
    elementRegister(index,4) = hep%elements(i)%pointIndexes(2)
    index = index + 1
    elementRegister(index,1) = hep%elements(i)%pointIndexes(2)
    elementRegister(index,2) = hep%elements(i)%pointIndexes(7)
    elementRegister(index,3) = hep%elements(i)%pointIndexes(3)
    elementRegister(index,4) = hep%elements(i)%pointIndexes(5)
   end if
  else
   elementRegister(index,1) = hep%elements(i)%pointIndexes(5)
   elementRegister(index,2) = hep%elements(i)%pointIndexes(7)
   elementRegister(index,3) = hep%elements(i)%pointIndexes(6)
   elementRegister(index,4) = hep%elements(i)%pointIndexes(1)
   if(q6<q7) then
    index = index + 1
    elementRegister(index,1) = hep%elements(i)%pointIndexes(1)
    elementRegister(index,2) = hep%elements(i)%pointIndexes(2)
    elementRegister(index,3) = hep%elements(i)%pointIndexes(3)
    elementRegister(index,4) = hep%elements(i)%pointIndexes(6)
    index = index + 1
    elementRegister(index,1) = hep%elements(i)%pointIndexes(6)
    elementRegister(index,2) = hep%elements(i)%pointIndexes(7)
    elementRegister(index,3) = hep%elements(i)%pointIndexes(3)
    elementRegister(index,4) = hep%elements(i)%pointIndexes(1)
   else
    index = index + 1
    elementRegister(index,1) = hep%elements(i)%pointIndexes(1)
    elementRegister(index,2) = hep%elements(i)%pointIndexes(6)
    elementRegister(index,3) = hep%elements(i)%pointIndexes(2)
    elementRegister(index,4) = hep%elements(i)%pointIndexes(7)
    index = index + 1
    elementRegister(index,1) = hep%elements(i)%pointIndexes(2)
    elementRegister(index,2) = hep%elements(i)%pointIndexes(7)
    elementRegister(index,3) = hep%elements(i)%pointIndexes(3)
    elementRegister(index,4) = hep%elements(i)%pointIndexes(1)
   end if
  end if


  index = index + 1

  if(q5<q8) then
   elementRegister(index,1) = hep%elements(i)%pointIndexes(1)
   elementRegister(index,2) = hep%elements(i)%pointIndexes(3)
   elementRegister(index,3) = hep%elements(i)%pointIndexes(4)
   elementRegister(index,4) = hep%elements(i)%pointIndexes(5)
   if(q7<q8) then
    index = index + 1
    elementRegister(index,1) = hep%elements(i)%pointIndexes(5)
    elementRegister(index,2) = hep%elements(i)%pointIndexes(8)
    elementRegister(index,3) = hep%elements(i)%pointIndexes(7)
    elementRegister(index,4) = hep%elements(i)%pointIndexes(4)
    index = index + 1
    elementRegister(index,1) = hep%elements(i)%pointIndexes(3)
    elementRegister(index,2) = hep%elements(i)%pointIndexes(7)
    elementRegister(index,3) = hep%elements(i)%pointIndexes(4)
    elementRegister(index,4) = hep%elements(i)%pointIndexes(5)
   else
    index = index + 1
    elementRegister(index,1) = hep%elements(i)%pointIndexes(5)
    elementRegister(index,2) = hep%elements(i)%pointIndexes(8)
    elementRegister(index,3) = hep%elements(i)%pointIndexes(7)
    elementRegister(index,4) = hep%elements(i)%pointIndexes(3)
    index = index + 1
    elementRegister(index,1) = hep%elements(i)%pointIndexes(3)
    elementRegister(index,2) = hep%elements(i)%pointIndexes(8)
    elementRegister(index,3) = hep%elements(i)%pointIndexes(4)
    elementRegister(index,4) = hep%elements(i)%pointIndexes(5)
   end if
  else
   elementRegister(index,1) = hep%elements(i)%pointIndexes(5)
   elementRegister(index,2) = hep%elements(i)%pointIndexes(8)
   elementRegister(index,3) = hep%elements(i)%pointIndexes(7)
   elementRegister(index,4) = hep%elements(i)%pointIndexes(1)
   if(q7<q8) then
    index = index + 1
    elementRegister(index,1) = hep%elements(i)%pointIndexes(3)
    elementRegister(index,2) = hep%elements(i)%pointIndexes(4)
    elementRegister(index,3) = hep%elements(i)%pointIndexes(1)
    elementRegister(index,4) = hep%elements(i)%pointIndexes(7)
    index = index + 1
    elementRegister(index,1) = hep%elements(i)%pointIndexes(4)
    elementRegister(index,2) = hep%elements(i)%pointIndexes(7)
    elementRegister(index,3) = hep%elements(i)%pointIndexes(8)
    elementRegister(index,4) = hep%elements(i)%pointIndexes(1)
   else
    index = index + 1
    elementRegister(index,1) = hep%elements(i)%pointIndexes(1)
    elementRegister(index,2) = hep%elements(i)%pointIndexes(3)
    elementRegister(index,3) = hep%elements(i)%pointIndexes(4)
    elementRegister(index,4) = hep%elements(i)%pointIndexes(8)
    index = index + 1
    elementRegister(index,1) = hep%elements(i)%pointIndexes(3)
    elementRegister(index,2) = hep%elements(i)%pointIndexes(7)
    elementRegister(index,3) = hep%elements(i)%pointIndexes(8)
    elementRegister(index,4) = hep%elements(i)%pointIndexes(1)
   end if
  end if
 end do

 do i=1,brp%numberOfTriangularBFs
  faceRegister(i,:) = brp%triangularFaces(i)%faceIndexes
 end do

 index = brp%numberOfTriangularBFs
 do i=1,brp%numberOfRectangularBFs
  index = index + 1
  faceRegister(index,1) = brp%rectangularFaces(i)%faceIndexes(1)
  faceRegister(index,2) = brp%rectangularFaces(i)%faceIndexes(2)
  faceRegister(index,3) = brp%rectangularFaces(i)%faceIndexes(3)
  index = index + 1
  faceRegister(index,1) = brp%rectangularFaces(i)%faceIndexes(1)
  faceRegister(index,2) = brp%rectangularFaces(i)%faceIndexes(3)
  faceRegister(index,3) = brp%rectangularFaces(i)%faceIndexes(4)
 end do

 do i=1,brp%numberOfTriangularBFs
  segmentRegister(i) = brp%triangularFaces(i)%surfaceSegmentNumber
  elementContainingFace(i) = brp%triangularFaces(i)%elementContainingFace
 end do

 index = brp%numberOfTriangularBFs
 do i=1,brp%numberOfRectangularBFs
  index = index + 1
   segmentRegister(index) = brp%rectangularFaces(i)%surfaceSegmentNumber
  index = index + 1
   segmentRegister(index) = brp%rectangularFaces(i)%surfaceSegmentNumber
 end do

 write(OUTFILE) numberOfTetEl,crp%numberOfPoints,numberOfTriBFs
 write(OUTFILE) ((elementRegister(i,j),i=1,numberOfTetEl),j=1,4)
 write(OUTFILE) ((crp%points(i,j),i=1,crp%numberOfPoints),j=1,3)
 write(OUTFILE) ((faceRegister(i,j),i=1,numberOfTriBFs),j=1,3),&
  (elementContainingFace(i),i=1,numberOfTriBFs),(segmentRegister(i),i=1,numberOfTriBFs)

 deallocate(elementRegister,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: makeNewPltFile couldn't deallocate"
 deallocate(faceRegister,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: makeNewPltFile couldn't deallocate"
 deallocate(segmentRegister,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: makeNewPltFile couldn't deallocate"

 write(*,*) "Finished writing new plt file"

 end subroutine makeNewPltFile
!-------------------------------------------------------------------------
 subroutine makeCalculationData(crp,crpPrev,crpPrev2,hep,prp,pyp,tep,grp,brp,doVisualization)

 IMPLICIT NONE

 type(CoordinateRegisterData) :: crp,crpPrev,crpPrev2 ! contains coordinate registers for grid
 type(HexahedralElementData),pointer :: hep ! contains element register for hexahedra
 type(PrismElementData),pointer :: prp ! contains element register for prisms
 type(PyramidElementData),pointer :: pyp ! contains element register for pyramids 
 type(TetrahedralElementData),pointer :: tep ! contains element register for tetrahedra
 type(GridData) ::  grp ! contains grid data  
 type(BoundaryRegisterData) :: brp
 logical :: doVisualization

 integer :: allocateStatus,i,dataSize,nodeSize

 type(SideData) :: sd

  grp%numberOfNodes = crp%numberOfPoints 
  nullify(grp%seedNodeRegister)

  write(*,*) ""
  write(*,*) "Making calculation data..."

! allocate

  allocate(grp%sdp,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: makeCalculationData out of memory" 
 
  dataSize = max(grp%numberOfNodes*2,100)
  nodeSize = dataSize

  call setUpSideRegister(grp%sdp,grp%numberOfNodes,dataSize,nodeSize) 

  allocate(grp%controlVolumeNodeSize(grp%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: makeCalculationData out of memory" 
  allocate(grp%prevControlVolume(grp%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: makeCalculationData out of memory"
  allocate(grp%prevControlVolume2(grp%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: makeCalculationData out of memory"


!  if(doVisualization) then 
!  end if

! initialize

  grp%numberOfSides = 0
  grp%brp%numberOfBoundarySides = 0
  grp%controlVolumeNodeSize = 0
  grp%prevControlVolume = 0.0
  grp%prevControlVolume2 = 0.0

  write(*,*) "Processing hexahedral elements..."
  call makeCalcDataForHex(crp,crpPrev,crpPrev2,hep,grp,brp,doVisualization)
  call cleanUpHexElData(hep)
  write(*,*) "Processing prism elements..."
  call makeCalcDataForPri(crp,crpPrev,crpPrev2,prp,grp,brp,doVisualization)
  call cleanUpPriElData(prp)
  write(*,*) "Processing pyramid elements..."
  call makeCalcDataForPyr(crp,crpPrev,crpPrev2,pyp,grp,brp,doVisualization)
  call cleanUpPyrElData(pyp)
  write(*,*) "Processing tetrahedral elements..."
  call makeCalcDataForTet(crp,crpPrev,crpPrev2,tep,grp,brp,doVisualization)
  call cleanUpTetElData(tep)

  write(*,*) "Processing boundary faces..."
  brp%numberOfEngineInletSideLimit = 0
  call makeCalcDataForBoundary(crp,crpPrev,crpPrev2,brp)
  write(*,*) "Processing engine inlet faces..."
  call makeEngineInflowData(crp,brp)

  if(associated(brp%triangularFaces)) then 
   deallocate(brp%triangularFaces,stat=allocateStatus)
   if (allocateStatus /= 0) STOP "ERROR: makeCalculationData couldn't deallocate" 
  end if
  if(associated(brp%rectangularFaces)) then 
   deallocate(brp%rectangularFaces,stat=allocateStatus)
   if (allocateStatus /= 0) STOP "ERROR: makeCalculationData couldn't deallocate" 
  end if

  call prepareDataStructure(grp%sdp)
  
  grp%numberOfSides = grp%sdp%numberOfSides  
  write(*,*) "Number of sides in fine grid: ",grp%numberOfSides
  write(*,*) "Number of blocks: ",grp%sdp%numberOfDataBlocks," (",grp%sdp%dataBlockSize,")" 
 end subroutine makeCalculationData
!-------------------------------------------------------------------------
 subroutine registerSideGCL(p1,p2,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
 ! creates or updates side

 IMPLICIT NONE

 integer :: p1,p2
 double precision :: r1(3,3),r2(3,3),r3(3,3),r4(3,3),rv(3,3),vr1(3,2),vr2(3,2),nodeVolume(3)
 type(CoordinateRegisterData) :: crp
 type(GridData) :: grp
 logical :: doVisualization

 integer :: newSideCreated
 type(SideData) :: sd
 double precision :: Cxx1(3,3),Cxx2(3,3),Cyy1(3,2),Cyy2(3,2),GCL(2),Cxx(3,3)

 Cxx1(1,:) = 0.5*(r1(2,:)*r2(3,:)-r1(3,:)*r2(2,:)) 
 Cxx2(1,:) = 0.5*(r3(2,:)*r4(3,:)-r3(3,:)*r4(2,:))
 Cxx1(2,:) = 0.5*(r1(3,:)*r2(1,:)-r1(1,:)*r2(3,:)) 
 Cxx2(2,:) = 0.5*(r3(3,:)*r4(1,:)-r3(1,:)*r4(3,:))
 Cxx1(3,:) = 0.5*(r1(1,:)*r2(2,:)-r1(2,:)*r2(1,:)) 
 Cxx2(3,:) = 0.5*(r3(1,:)*r4(2,:)-r3(2,:)*r4(1,:))

 Cyy1(1,1) = 0.25*(r1(2,1)*r2(3,2)-r1(3,1)*r2(2,2)) + 0.25*(r1(2,2)*r2(3,1)-r1(3,2)*r2(2,1))
 Cyy2(1,1) = 0.25*(r3(2,1)*r4(3,2)-r3(3,1)*r4(2,2)) + 0.25*(r3(2,2)*r4(3,1)-r3(3,2)*r4(2,1))
 Cyy1(2,1) = 0.25*(r1(3,1)*r2(1,2)-r1(1,1)*r2(3,2)) + 0.25*(r1(3,2)*r2(1,1)-r1(1,2)*r2(3,1))
 Cyy2(2,1) = 0.25*(r3(3,1)*r4(1,2)-r3(1,1)*r4(3,2)) + 0.25*(r3(3,2)*r4(1,1)-r3(1,2)*r4(3,1))
 Cyy1(3,1) = 0.25*(r1(1,1)*r2(2,2)-r1(2,1)*r2(1,2)) + 0.25*(r1(1,2)*r2(2,1)-r1(2,2)*r2(1,1))
 Cyy2(3,1) = 0.25*(r3(1,1)*r4(2,2)-r3(2,1)*r4(1,2)) + 0.25*(r3(1,2)*r4(2,1)-r3(2,2)*r4(1,1))

 Cyy1(1,2) = 0.25*(r1(2,2)*r2(3,3)-r1(3,2)*r2(2,3)) + 0.25*(r1(2,3)*r2(3,2)-r1(3,3)*r2(2,2))
 Cyy2(1,2) = 0.25*(r3(2,2)*r4(3,3)-r3(3,2)*r4(2,3)) + 0.25*(r3(2,3)*r4(3,2)-r3(3,3)*r4(2,2))
 Cyy1(2,2) = 0.25*(r1(3,2)*r2(1,3)-r1(1,2)*r2(3,3)) + 0.25*(r1(3,3)*r2(1,2)-r1(1,3)*r2(3,2))
 Cyy2(2,2) = 0.25*(r3(3,2)*r4(1,3)-r3(1,2)*r4(3,3)) + 0.25*(r3(3,3)*r4(1,2)-r3(1,3)*r4(3,2))
 Cyy1(3,2) = 0.25*(r1(1,2)*r2(2,3)-r1(2,2)*r2(1,3)) + 0.25*(r1(1,3)*r2(2,2)-r1(2,3)*r2(1,2))
 Cyy2(3,2) = 0.25*(r3(1,2)*r4(2,3)-r3(2,2)*r4(1,3)) + 0.25*(r3(1,3)*r4(2,2)-r3(2,3)*r4(1,2))

 GCL(1) = (sum((Cxx1(:,1)+Cyy1(:,1)+Cxx1(:,2))*vr1(:,1))+sum((Cxx2(:,1)+Cyy2(:,1)+Cxx2(:,2))*vr2(:,1)))/3.
 GCL(2) = (sum((Cxx1(:,2)+Cyy1(:,2)+Cxx1(:,3))*vr1(:,2))+sum((Cxx2(:,2)+Cyy2(:,2)+Cxx2(:,3))*vr2(:,2)))/3.

 Cxx = Cxx1+Cxx2

 nodeVolume(:) = (rv(1,:)*Cxx(1,:)+rv(2,:)*Cxx(2,:)+rv(3,:)*Cxx(3,:))/3.

 grp%controlVolumeNodeSize(p1) = grp%controlVolumeNodeSize(p1)+nodeVolume(1)
 grp%controlVolumeNodeSize(p2) = grp%controlVolumeNodeSize(p2)+nodeVolume(1)
 grp%prevControlVolume(p1) = grp%prevControlVolume(p1)+nodeVolume(2)
 grp%prevControlVolume(p2) = grp%prevControlVolume(p2)+nodeVolume(2)
 grp%prevControlVolume2(p1) = grp%prevControlVolume2(p1)+nodeVolume(3)
 grp%prevControlVolume2(p2) = grp%prevControlVolume2(p2)+nodeVolume(3)

 if(nodeVolume(1).le.0) write(*,*) "WARNING: Control volume for side ",p1,p2," has negative contribution"

 sd%sideCoefficients = 0.5*Cxx(:,1)  ! store half the side coefficient
 sd%sideLength = sqrt(sum((getCoor(p2,crp)-getCoor(p1,crp))*(getCoor(p2,crp)-getCoor(p1,crp))))
 sd%sideIndexes(1) = p1
 sd%sideIndexes(2) = p2
 sd%GCLCoefficient = 0.5*GCL(1)
 sd%GCLCoefficientP = 0.5*GCL(2)

 if(p1==p2) write(*,*) "WARNING: Side has identical indexes in registerSide: ",p1,p2
 newSideCreated = registerSideData(sd,grp%sdp)

 if(newSideCreated==0) then
  if(p1.le.grp%brp%numberOfBoundaryNodes.and.p2.le.grp%brp%numberOfBoundaryNodes) then
   ! count number of boundary sides
   grp%brp%numberOfBoundarySides = grp%brp%numberOfBoundarySides + 1
  end if
 end if

 end subroutine registerSideGCL
!-------------------------------------------------------------------------
 subroutine registerSide(p1,p2,r1,r2,r3,r4,rv,crp,grp,doVisualization,nodeVolume)
 ! creates or updates side

 IMPLICIT NONE

 integer :: p1,p2
 double precision :: r1(3),r2(3),r3(3),r4(3),rv(3),nodeVolume
 type(CoordinateRegisterData) :: crp
 type(GridData) :: grp
 logical :: doVisualization

 integer :: newSideCreated
 type(SideData) :: sd
 double precision :: Cxx(3)

 Cxx(1) = 0.5*(r1(2)*r2(3)-r1(3)*r2(2)) + 0.5*(r3(2)*r4(3)-r3(3)*r4(2))
 Cxx(2) = 0.5*(r1(3)*r2(1)-r1(1)*r2(3)) + 0.5*(r3(3)*r4(1)-r3(1)*r4(3))
 Cxx(3) = 0.5*(r1(1)*r2(2)-r1(2)*r2(1)) + 0.5*(r3(1)*r4(2)-r3(2)*r4(1))

 nodeVolume = (rv(1)*Cxx(1)+rv(2)*Cxx(2)+rv(3)*Cxx(3))/3.

 grp%controlVolumeNodeSize(p1) = grp%controlVolumeNodeSize(p1)+nodeVolume
 grp%controlVolumeNodeSize(p2) = grp%controlVolumeNodeSize(p2)+nodeVolume
 if(nodeVolume.le.0) then 
  write(*,*) "WARNING: Control volume for side ",p1,p2," has negative contribution"
  write(*,*) grp%nodeSwitchRegister(p1),grp%nodeSwitchRegister(p2) 
  pause
 end if
 sd%sideCoefficients = 0.5*Cxx  ! store half the side coefficient
 sd%sideLength = sqrt(sum((getCoor(p2,crp)-getCoor(p1,crp))*(getCoor(p2,crp)-getCoor(p1,crp)))) 
 sd%sideIndexes(1) = p1
 sd%sideIndexes(2) = p2
 if(p1==p2) write(*,*) "WARNING: Side has identical indexes in registerSide: ",p1,p2
 newSideCreated = registerSideData(sd,grp%sdp)

 if(newSideCreated==0) then 
  if(p1.le.grp%brp%numberOfBoundaryNodes.and.p2.le.grp%brp%numberOfBoundaryNodes) then 
   ! count number of boundary sides
   grp%brp%numberOfBoundarySides = grp%brp%numberOfBoundarySides + 1
  end if
 end if

 ! do visualization registration
 end subroutine registerSide
!-------------------------------------------------------------------------
 subroutine makeCalcDataForHex(crp,crpPrev,crpPrev2,hep,grp,brp,doVisualization)
 ! makes side coefficients from hexahedral elements and splits
 ! data structure into side based

 IMPLICIT NONE
 
 type(CoordinateRegisterData) :: crp,crpPrev,crpPrev2 ! contains coordinate register for grid
 type(HexahedralElementData) :: hep ! contains element register for rectangular grid
 type(GridData) ::  grp ! contains grid data  
 type(BoundaryRegisterData) :: brp
 logical :: doVisualization

 integer :: i
 integer :: p1,p2,p3,p4,p5,p6,p7,p8,q1,q2,q3,q4,q5,q6,q7,q8
 double precision :: xmm(3,3),xm1(3,3),xm2(3,3),xm3(3,3),xb1(3,3),xb2(3,3),xb3(3,3)
 double precision :: x1(3,3),x2(3,3),x3(3,3),x4(3,3),x5(3,3),x6(3,3),x7(3,3),x8(3,3) 
 double precision :: r1(3,3),r2(3,3),r3(3,3),r4(3,3),rv(3,3)
 double precision :: vr1(3,2),vr2(3,2)
 double precision :: elementVolume(3),nodeVolume(3)

 real :: alpha

 type(SHEData),pointer :: temp

 do i = 1,hep%numberOfElements ! do for each element
  temp => getHexElement(i,hep) 

  p1=temp%pointIndexes(1)
  p2=temp%pointIndexes(2)
  p3=temp%pointIndexes(3)
  p4=temp%pointIndexes(4)
  p5=temp%pointIndexes(5)
  p6=temp%pointIndexes(6)
  p7=temp%pointIndexes(7)
  p8=temp%pointIndexes(8)


  q1=grp%nodeSwitchRegister(p1)
  q2=grp%nodeSwitchRegister(p2)
  q3=grp%nodeSwitchRegister(p3)
  q4=grp%nodeSwitchRegister(p4)
  q5=grp%nodeSwitchRegister(p5)
  q6=grp%nodeSwitchRegister(p6)
  q7=grp%nodeSwitchRegister(p7)
  q8=grp%nodeSwitchRegister(p8)

  x1(:,1) = getCoor(p1,crp)
  x2(:,1) = getCoor(p2,crp)
  x3(:,1) = getCoor(p3,crp)
  x4(:,1) = getCoor(p4,crp)
  x5(:,1) = getCoor(p5,crp)
  x6(:,1) = getCoor(p6,crp)
  x7(:,1) = getCoor(p7,crp)
  x8(:,1) = getCoor(p8,crp)
  x1(:,2) = getCoor(p1,crpPrev)
  x2(:,2) = getCoor(p2,crpPrev)
  x3(:,2) = getCoor(p3,crpPrev)
  x4(:,2) = getCoor(p4,crpPrev)
  x5(:,2) = getCoor(p5,crpPrev)
  x6(:,2) = getCoor(p6,crpPrev)
  x7(:,2) = getCoor(p7,crpPrev)
  x8(:,2) = getCoor(p8,crpPrev)
  x1(:,3) = getCoor(p1,crpPrev2)
  x2(:,3) = getCoor(p2,crpPrev2)
  x3(:,3) = getCoor(p3,crpPrev2)
  x4(:,3) = getCoor(p4,crpPrev2)
  x5(:,3) = getCoor(p5,crpPrev2)
  x6(:,3) = getCoor(p6,crpPrev2)
  x7(:,3) = getCoor(p7,crpPrev2)
  x8(:,3) = getCoor(p8,crpPrev2)

  xmm = 0.25*(x1+x3+x5+x7)

  ! node 1:

  xm1 = 0.5*(x1+x2)
  xm2 = 0.5*(x1+x4)
  xm3 = 0.5*(x1+x5)
  xb1 = 0.5*(x1+x3)
  if(q5<q6) then
   xb2 = 0.5*(x2+x5)
  else
   xb2 = 0.5*(x1+x6)
  end if
  if(q5<q8) then
   xb3 = 0.5*(x4+x5)
  else
   xb3 = 0.5*(x1+x8)
  end if

  ! side 1-2

  r1 = xmm - xm1
  r2 = xb2 - xm1
  r3 = xb1 - xm1
  r4 = xmm - xm1
  rv = xm1-x1 

  vr1(:,1) = xm1(:,1) - xm1(:,2) 
  vr1(:,2) = xm1(:,2) - xm1(:,3) 
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2) 
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3) 
  vr1(:,1) = vr1(:,1) + xb2(:,1) - xb2(:,2) 
  vr1(:,2) = vr1(:,2) + xb2(:,2) - xb2(:,3) 
  vr1 = vr1/3.0

  vr2(:,1) = xm1(:,1) - xm1(:,2) 
  vr2(:,2) = xm1(:,2) - xm1(:,3) 
  vr2(:,1) = vr2(:,1) + xb1(:,1) - xb1(:,2) 
  vr2(:,2) = vr2(:,2) + xb1(:,2) - xb1(:,3) 
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2) 
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3) 
  vr2 = vr2/3.0

  call registerSideGCL(p1,p2,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = 2.*nodeVolume

! side 1-4

  r1 = xmm - xm2
  r2 = xb1 - xm2
  r3 = xb3 - xm2
  r4 = xmm - xm2
  rv = xm2-x1 

  vr1(:,1) = xm2(:,1) - xm2(:,2)
  vr1(:,2) = xm2(:,2) - xm2(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb1(:,1) - xb1(:,2)
  vr1(:,2) = vr1(:,2) + xb1(:,2) - xb1(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm2(:,1) - xm2(:,2)
  vr2(:,2) = xm2(:,2) - xm2(:,3)
  vr2(:,1) = vr2(:,1) + xb3(:,1) - xb3(:,2)
  vr2(:,2) = vr2(:,2) + xb3(:,2) - xb3(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0

  call registerSideGCL(p1,p4,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

! side 1-5

  r1 = xmm - xm3 
  r2 = xb3 - xm3
  r3 = xb2 - xm3
  r4 = xmm - xm3
  rv = xm3-x1 

  vr1(:,1) = xm3(:,1) - xm3(:,2)
  vr1(:,2) = xm3(:,2) - xm3(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb3(:,1) - xb3(:,2)
  vr1(:,2) = vr1(:,2) + xb3(:,2) - xb3(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm3(:,1) - xm3(:,2)
  vr2(:,2) = xm3(:,2) - xm3(:,3)
  vr2(:,1) = vr2(:,1) + xb2(:,1) - xb2(:,2)
  vr2(:,2) = vr2(:,2) + xb2(:,2) - xb2(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0


  call registerSideGCL(p1,p5,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

 ! node 6:

  xm1 = 0.5*(x6+x5)
  xm2 = 0.5*(x6+x7)
  xm3 = 0.5*(x6+x2)
  xb1 = 0.5*(x5+x7)
  if(q5<q6) then
   xb2 = 0.5*(x2+x5)
  else
   xb2 = 0.5*(x1+x6)
  end if
  if(q6<q7) then
   xb3 = 0.5*(x6+x3)
  else
   xb3 = 0.5*(x2+x7)
  end if

  ! side 6-5 

  r1 = xmm - xm1
  r2 = xb2 - xm1
  r3 = xb1 - xm1
  r4 = xmm - xm1
  rv = xm1-x6

  vr1(:,1) = xm1(:,1) - xm1(:,2)
  vr1(:,2) = xm1(:,2) - xm1(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb2(:,1) - xb2(:,2)
  vr1(:,2) = vr1(:,2) + xb2(:,2) - xb2(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm1(:,1) - xm1(:,2)
  vr2(:,2) = xm1(:,2) - xm1(:,3)
  vr2(:,1) = vr2(:,1) + xb1(:,1) - xb1(:,2)
  vr2(:,2) = vr2(:,2) + xb1(:,2) - xb1(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0

  call registerSideGCL(p6,p5,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

! side 6-7 

  r1 = xmm - xm2
  r2 = xb1 - xm2
  r3 = xb3 - xm2
  r4 = xmm - xm2
  rv = xm2-x6

  vr1(:,1) = xm2(:,1) - xm2(:,2)
  vr1(:,2) = xm2(:,2) - xm2(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb1(:,1) - xb1(:,2)
  vr1(:,2) = vr1(:,2) + xb1(:,2) - xb1(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm2(:,1) - xm2(:,2)
  vr2(:,2) = xm2(:,2) - xm2(:,3)
  vr2(:,1) = vr2(:,1) + xb3(:,1) - xb3(:,2)
  vr2(:,2) = vr2(:,2) + xb3(:,2) - xb3(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0


  call registerSideGCL(p6,p7,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

! side 6-2 

  r1 = xmm - xm3
  r2 = xb3 - xm3
  r3 = xb2 - xm3
  r4 = xmm - xm3
  rv = xm3-x6

  vr1(:,1) = xm3(:,1) - xm3(:,2)
  vr1(:,2) = xm3(:,2) - xm3(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb3(:,1) - xb3(:,2)
  vr1(:,2) = vr1(:,2) + xb3(:,2) - xb3(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm3(:,1) - xm3(:,2)
  vr2(:,2) = xm3(:,2) - xm3(:,3)
  vr2(:,1) = vr2(:,1) + xb2(:,1) - xb2(:,2)
  vr2(:,2) = vr2(:,2) + xb2(:,2) - xb2(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0


  call registerSideGCL(p6,p2,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

  ! node 3:

  xm1 = 0.5*(x3+x4)
  xm2 = 0.5*(x2+x3)
  xm3 = 0.5*(x3+x7)
  xb1 = 0.5*(x1+x3)
  if(q7<q8) then
   xb2 = 0.5*(x7+x4)
  else
   xb2 = 0.5*(x3+x8)
  end if
  if(q6<q7) then
   xb3 = 0.5*(x6+x3)
  else
   xb3 = 0.5*(x2+x7)
  end if

  ! side 3-4 

  r1 = xmm - xm1
  r2 = xb2 - xm1
  r3 = xb1 - xm1
  r4 = xmm - xm1
  rv = xm1-x3

  vr1(:,1) = xm1(:,1) - xm1(:,2)
  vr1(:,2) = xm1(:,2) - xm1(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb2(:,1) - xb2(:,2)
  vr1(:,2) = vr1(:,2) + xb2(:,2) - xb2(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm1(:,1) - xm1(:,2)
  vr2(:,2) = xm1(:,2) - xm1(:,3)
  vr2(:,1) = vr2(:,1) + xb1(:,1) - xb1(:,2)
  vr2(:,2) = vr2(:,2) + xb1(:,2) - xb1(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0


  call registerSideGCL(p3,p4,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

! side 3-2 

  r1 = xmm - xm2
  r2 = xb1 - xm2
  r3 = xb3 - xm2
  r4 = xmm - xm2
  rv = xm2-x3

  vr1(:,1) = xm2(:,1) - xm2(:,2)
  vr1(:,2) = xm2(:,2) - xm2(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb1(:,1) - xb1(:,2)
  vr1(:,2) = vr1(:,2) + xb1(:,2) - xb1(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm2(:,1) - xm2(:,2)
  vr2(:,2) = xm2(:,2) - xm2(:,3)
  vr2(:,1) = vr2(:,1) + xb3(:,1) - xb3(:,2)
  vr2(:,2) = vr2(:,2) + xb3(:,2) - xb3(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0


  call registerSideGCL(p3,p2,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

! side 3-7 

  r1 = xmm - xm3 
  r2 = xb3 - xm3
  r3 = xb2 - xm3
  r4 = xmm - xm3
  rv = xm3-x3

  vr1(:,1) = xm3(:,1) - xm3(:,2)
  vr1(:,2) = xm3(:,2) - xm3(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb3(:,1) - xb3(:,2)
  vr1(:,2) = vr1(:,2) + xb3(:,2) - xb3(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm3(:,1) - xm3(:,2)
  vr2(:,2) = xm3(:,2) - xm3(:,3)
  vr2(:,1) = vr2(:,1) + xb2(:,1) - xb2(:,2)
  vr2(:,2) = vr2(:,2) + xb2(:,2) - xb2(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0

  call registerSideGCL(p3,p7,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

 ! node 8:

  xm1 = 0.5*(x8+x7)
  xm2 = 0.5*(x8+x5)
  xm3 = 0.5*(x8+x4)
  xb1 = 0.5*(x5+x7)
  if(q7<q8) then
   xb2 = 0.5*(x7+x4)
  else
   xb2 = 0.5*(x3+x8)
  end if
  if(q5<q8) then
   xb3 = 0.5*(x5+x4)
  else
   xb3 = 0.5*(x1+x8)
  end if

  ! side 8-7 

  r1 = xmm - xm1
  r2 = xb2 - xm1
  r3 = xb1 - xm1
  r4 = xmm - xm1
  rv = xm1-x8

  vr1(:,1) = xm1(:,1) - xm1(:,2)
  vr1(:,2) = xm1(:,2) - xm1(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb2(:,1) - xb2(:,2)
  vr1(:,2) = vr1(:,2) + xb2(:,2) - xb2(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm1(:,1) - xm1(:,2)
  vr2(:,2) = xm1(:,2) - xm1(:,3)
  vr2(:,1) = vr2(:,1) + xb1(:,1) - xb1(:,2)
  vr2(:,2) = vr2(:,2) + xb1(:,2) - xb1(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0


  call registerSideGCL(p8,p7,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

! side 8-5 

  r1 = xmm - xm2
  r2 = xb1 - xm2
  r3 = xb3 - xm2
  r4 = xmm - xm2
  rv = xm2-x8

  vr1(:,1) = xm2(:,1) - xm2(:,2)
  vr1(:,2) = xm2(:,2) - xm2(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb1(:,1) - xb1(:,2)
  vr1(:,2) = vr1(:,2) + xb1(:,2) - xb1(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm2(:,1) - xm2(:,2)
  vr2(:,2) = xm2(:,2) - xm2(:,3)
  vr2(:,1) = vr2(:,1) + xb3(:,1) - xb3(:,2)
  vr2(:,2) = vr2(:,2) + xb3(:,2) - xb3(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0


  call registerSideGCL(p8,p5,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

! side 8-4 

  r1 = xmm - xm3
  r2 = xb3 - xm3
  r3 = xb2 - xm3
  r4 = xmm - xm3
  rv = xm3-x8

  vr1(:,1) = xm3(:,1) - xm3(:,2)
  vr1(:,2) = xm3(:,2) - xm3(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb3(:,1) - xb3(:,2)
  vr1(:,2) = vr1(:,2) + xb3(:,2) - xb3(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm3(:,1) - xm3(:,2)
  vr2(:,2) = xm3(:,2) - xm3(:,3)
  vr2(:,1) = vr2(:,1) + xb2(:,1) - xb2(:,2)
  vr2(:,2) = vr2(:,2) + xb2(:,2) - xb2(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0


  call registerSideGCL(p8,p4,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

  if(elementVolume(1).le.0) write(*,*) "WARNING: Negative volume for hexahedral element number: ",i," (",elementvolume(1),") "

 end do

 end subroutine makeCalcDataForHex
!-------------------------------------------------------------------------
 subroutine makeCalcDataForPri(crp,crpPrev,crpPrev2,prp,grp,brp,doVisualization)
 ! makes side coefficients from prism elements and splits
 ! data structure into side based

 IMPLICIT NONE
 
 type(CoordinateRegisterData) :: crp,crpPrev,crpPrev2 ! contains coordinate register for grid
 type(PrismElementData) :: prp ! contains element register for rectangular grid
 type(GridData) ::  grp ! contains grid data  
 type(BoundaryRegisterData) :: brp
 logical :: doVisualization

 integer :: i 
 integer :: p1,p2,p3,p4,p5,p6,q1,q2,q3,q4,q5,q6
 double precision :: xmm(3,3),xm1(3,3),xm2(3,3),xm3(3,3),xb1(3,3),xb2(3,3),xb3(3,3)
 double precision :: x1(3,3),x2(3,3),x3(3,3),x4(3,3),x5(3,3),x6(3,3)
 double precision :: r1(3,3),r2(3,3),r3(3,3),r4(3,3),rv(3,3),vr1(3,2),vr2(3,2)
 double precision :: elementVolume(3),nodeVolume(3)

 type(SPRData),pointer :: temp

 do i = 1,prp%numberOfElements ! do for each element
  temp => getPriElement(i,prp) 

  p1=temp%pointIndexes(1)
  p2=temp%pointIndexes(2)
  p3=temp%pointIndexes(3)
  p4=temp%pointIndexes(4)
  p5=temp%pointIndexes(5)
  p6=temp%pointIndexes(6)

  q1=grp%nodeSwitchRegister(p1)
  q2=grp%nodeSwitchRegister(p2)
  q3=grp%nodeSwitchRegister(p3)
  q4=grp%nodeSwitchRegister(p4)
  q5=grp%nodeSwitchRegister(p5)
  q6=grp%nodeSwitchRegister(p6)

  x1(:,1) = getCoor(p1,crp)
  x2(:,1) = getCoor(p2,crp)
  x3(:,1) = getCoor(p3,crp)
  x4(:,1) = getCoor(p4,crp)
  x5(:,1) = getCoor(p5,crp)
  x6(:,1) = getCoor(p6,crp)
  x1(:,2) = getCoor(p1,crpPrev)
  x2(:,2) = getCoor(p2,crpPrev)
  x3(:,2) = getCoor(p3,crpPrev)
  x4(:,2) = getCoor(p4,crpPrev)
  x5(:,2) = getCoor(p5,crpPrev)
  x6(:,2) = getCoor(p6,crpPrev)
  x1(:,3) = getCoor(p1,crpPrev2)
  x2(:,3) = getCoor(p2,crpPrev2)
  x3(:,3) = getCoor(p3,crpPrev2)
  x4(:,3) = getCoor(p4,crpPrev2)
  x5(:,3) = getCoor(p5,crpPrev2)
  x6(:,3) = getCoor(p6,crpPrev2)

  if(q4<q5) then
   xb1 = 0.5*(x4+x2)
  else
   xb1 = 0.5*(x1+x5)
  end if
  if(q4<q6) then
   xb2 = 0.5*(x4+x3)
  else
   xb2 = 0.5*(x1+x6)
  end if
  if(q5<q6) then
   xb3 = 0.5*(x5+x3)
  else
   xb3 = 0.5*(x2+x6)
  end if

  xmm = (xb1+xb2+xb3)/3.

  ! node 1:

  xm1 = 0.5*(x1+x2)
  xm2 = 0.5*(x1+x3)
  xm3 = 0.5*(x1+x4)
  xb1 = (x1+x2+x3)/3.
  if(q4<q5) then
   xb2 = 0.5*(x4+x2)
  else
   xb2 = 0.5*(x1+x5)
  end if
  if(q4<q6) then
   xb3 = 0.5*(x4+x3)
  else
   xb3 = 0.5*(x1+x6)
  end if

  ! side 1-2

  r1 = xmm - xm1
  r2 = xb2 - xm1
  r3 = xb1 - xm1
  r4 = xmm - xm1
  rv = xm1-x1

  vr1(:,1) = xm1(:,1) - xm1(:,2)
  vr1(:,2) = xm1(:,2) - xm1(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb2(:,1) - xb2(:,2)
  vr1(:,2) = vr1(:,2) + xb2(:,2) - xb2(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm1(:,1) - xm1(:,2)
  vr2(:,2) = xm1(:,2) - xm1(:,3)
  vr2(:,1) = vr2(:,1) + xb1(:,1) - xb1(:,2)
  vr2(:,2) = vr2(:,2) + xb1(:,2) - xb1(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0


  call registerSideGCL(p1,p2,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = 2.*nodeVolume

! side 1-3

  r1 = xmm - xm2
  r2 = xb1 - xm2
  r3 = xb3 - xm2
  r4 = xmm - xm2
  rv = xm2-x1

  vr1(:,1) = xm2(:,1) - xm2(:,2)
  vr1(:,2) = xm2(:,2) - xm2(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb1(:,1) - xb1(:,2)
  vr1(:,2) = vr1(:,2) + xb1(:,2) - xb1(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm2(:,1) - xm2(:,2)
  vr2(:,2) = xm2(:,2) - xm2(:,3)
  vr2(:,1) = vr2(:,1) + xb3(:,1) - xb3(:,2)
  vr2(:,2) = vr2(:,2) + xb3(:,2) - xb3(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0


  call registerSideGCL(p1,p3,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

! side 1-4

  r1 = xmm - xm3 
  r2 = xb3 - xm3
  r3 = xb2 - xm3
  r4 = xmm - xm3
  rv = xm3-x1

  vr1(:,1) = xm3(:,1) - xm3(:,2)
  vr1(:,2) = xm3(:,2) - xm3(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb3(:,1) - xb3(:,2)
  vr1(:,2) = vr1(:,2) + xb3(:,2) - xb3(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm3(:,1) - xm3(:,2)
  vr2(:,2) = xm3(:,2) - xm3(:,3)
  vr2(:,1) = vr2(:,1) + xb2(:,1) - xb2(:,2)
  vr2(:,2) = vr2(:,2) + xb2(:,2) - xb2(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0


  call registerSideGCL(p1,p4,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

 ! node 6:

  xm1 = 0.5*(x6+x5)
  xm2 = 0.5*(x6+x4)
  xm3 = 0.5*(x6+x3)
  xb1 = (x4+x5+x6)/3.
  if(q5<q6) then
   xb2 = 0.5*(x5+x3)
  else
   xb2 = 0.5*(x2+x6)
  end if
  if(q4<q6) then
   xb3 = 0.5*(x4+x3)
  else
   xb3 = 0.5*(x1+x6)
  end if

  ! side 6-5 

  r1 = xmm - xm1
  r2 = xb2 - xm1
  r3 = xb1 - xm1
  r4 = xmm - xm1
  rv = xm1-x6

  vr1(:,1) = xm1(:,1) - xm1(:,2)
  vr1(:,2) = xm1(:,2) - xm1(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb2(:,1) - xb2(:,2)
  vr1(:,2) = vr1(:,2) + xb2(:,2) - xb2(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm1(:,1) - xm1(:,2)
  vr2(:,2) = xm1(:,2) - xm1(:,3)
  vr2(:,1) = vr2(:,1) + xb1(:,1) - xb1(:,2)
  vr2(:,2) = vr2(:,2) + xb1(:,2) - xb1(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0



  call registerSideGCL(p6,p5,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

! side 6-4 

  r1 = xmm - xm2
  r2 = xb1 - xm2
  r3 = xb3 - xm2
  r4 = xmm - xm2
  rv = xm2-x6

  vr1(:,1) = xm2(:,1) - xm2(:,2)
  vr1(:,2) = xm2(:,2) - xm2(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb1(:,1) - xb1(:,2)
  vr1(:,2) = vr1(:,2) + xb1(:,2) - xb1(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm2(:,1) - xm2(:,2)
  vr2(:,2) = xm2(:,2) - xm2(:,3)
  vr2(:,1) = vr2(:,1) + xb3(:,1) - xb3(:,2)
  vr2(:,2) = vr2(:,2) + xb3(:,2) - xb3(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0

  call registerSideGCL(p6,p4,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

! side 6-3 

  r1 = xmm - xm3
  r2 = xb3 - xm3
  r3 = xb2 - xm3
  r4 = xmm - xm3
  rv = xm3-x6

  vr1(:,1) = xm3(:,1) - xm3(:,2)
  vr1(:,2) = xm3(:,2) - xm3(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb3(:,1) - xb3(:,2)
  vr1(:,2) = vr1(:,2) + xb3(:,2) - xb3(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm3(:,1) - xm3(:,2)
  vr2(:,2) = xm3(:,2) - xm3(:,3)
  vr2(:,1) = vr2(:,1) + xb2(:,1) - xb2(:,2)
  vr2(:,2) = vr2(:,2) + xb2(:,2) - xb2(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0


  call registerSideGCL(p6,p3,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

  ! node 2:

  xm1 = 0.5*(x2+x5)
  xm2 = 0.5*(x2+x3)
  xm3 = 0.5*(x1+x2)
  xb3 = (x1+x2+x3)/3.
  if(q5<q6) then
   xb1 = 0.5*(x5+x3)
  else
   xb1 = 0.5*(x2+x6)
  end if
  if(q4<q5) then
   xb2 = 0.5*(x4+x2)
  else
   xb2 = 0.5*(x1+x5)
  end if

  ! side 2-3 

  r1 = xmm - xm2
  r2 = xb1 - xm2
  r3 = xb3 - xm2
  r4 = xmm - xm2
  rv = xm2-x2

  vr1(:,1) = xm2(:,1) - xm2(:,2)
  vr1(:,2) = xm2(:,2) - xm2(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb1(:,1) - xb1(:,2)
  vr1(:,2) = vr1(:,2) + xb1(:,2) - xb1(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm2(:,1) - xm2(:,2)
  vr2(:,2) = xm2(:,2) - xm2(:,3)
  vr2(:,1) = vr2(:,1) + xb3(:,1) - xb3(:,2)
  vr2(:,2) = vr2(:,2) + xb3(:,2) - xb3(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0


  call registerSideGCL(p2,p3,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

! side 2-5 

  r1 = xmm - xm1
  r2 = xb2 - xm1
  r3 = xb1 - xm1
  r4 = xmm - xm1
  rv = xm1-x2

  vr1(:,1) = xm1(:,1) - xm1(:,2)
  vr1(:,2) = xm1(:,2) - xm1(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb2(:,1) - xb2(:,2)
  vr1(:,2) = vr1(:,2) + xb2(:,2) - xb2(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm1(:,1) - xm1(:,2)
  vr2(:,2) = xm1(:,2) - xm1(:,3)
  vr2(:,1) = vr2(:,1) + xb1(:,1) - xb1(:,2)
  vr2(:,2) = vr2(:,2) + xb1(:,2) - xb1(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0


  call registerSideGCL(p2,p5,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

 ! node 4:

  xm1 = 0.5*(x4+x5)
  xm2 = 0.5*(x4+x6)
  xm3 = 0.5*(x1+x4)
  xb1 = (x4+x5+x6)/3.
  if(q4<q5) then
   xb2 = 0.5*(x4+x2)
  else
   xb2 = 0.5*(x1+x5)
  end if
  if(q4<q6) then
   xb3 = 0.5*(x4+x3)
  else
   xb3 = 0.5*(x1+x6)
  end if

  ! side 4-5 

  r1 = xmm - xm1
  r2 = xb1 - xm1
  r3 = xb2 - xm1
  r4 = xmm - xm1
  rv = xm1-x4

  vr1(:,1) = xm1(:,1) - xm1(:,2)
  vr1(:,2) = xm1(:,2) - xm1(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb1(:,1) - xb1(:,2)
  vr1(:,2) = vr1(:,2) + xb1(:,2) - xb1(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm1(:,1) - xm1(:,2)
  vr2(:,2) = xm1(:,2) - xm1(:,3)
  vr2(:,1) = vr2(:,1) + xb2(:,1) - xb2(:,2)
  vr2(:,2) = vr2(:,2) + xb2(:,2) - xb2(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0

  call registerSideGCL(p4,p5,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

  if(elementVolume(1).le.0) write(*,*) "WARNING: Negative volume for prism element number: ",i," (",elementvolume,") "
 end do
 end subroutine makeCalcDataForPri
!-------------------------------------------------------------------------
 subroutine makeCalcDataForPyr(crp,crpPrev,crpPrev2,pyp,grp,brp,doVisualization)
 ! makes side coefficients from pyramid elements and splits
 ! data structure into side based

 IMPLICIT NONE
 
 type(CoordinateRegisterData) :: crp,crpPrev,crpPrev2 ! contains coordinate register for grid
 type(PyramidElementData) :: pyp ! contains element register for rectangular grid
 type(GridData) ::  grp ! contains grid data  
 type(BoundaryRegisterData) :: brp
 logical :: doVisualization

 integer :: i
 integer :: p1,p2,p3,p4,p5
 double precision :: xmm(3,3),xm1(3,3),xm2(3,3),xm3(3,3),xm4(3,3),xb1(3,3),xb2(3,3),xb3(3,3),xb4(3,3)
 double precision :: x1(3,3),x2(3,3),x3(3,3),x4(3,3),x5(3,3)
 double precision :: r1(3,3),r2(3,3),r3(3,3),r4(3,3),rv(3,3),vr1(3,2),vr2(3,2)
 double precision :: elementVolume(3),nodeVolume(3)

 type(SPYData),pointer :: temp

 do i = 1,pyp%numberOfElements ! do for each element
  temp => getPyrElement(i,pyp) 

  p1=temp%pointIndexes(1)
  p2=temp%pointIndexes(2)
  p3=temp%pointIndexes(3)
  p4=temp%pointIndexes(4)
  p5=temp%pointIndexes(5)

  x1(:,1) = getCoor(p1,crp)
  x2(:,1) = getCoor(p2,crp)
  x3(:,1) = getCoor(p3,crp)
  x4(:,1) = getCoor(p4,crp)
  x5(:,1) = getCoor(p5,crp)
  x1(:,2) = getCoor(p1,crpPrev)
  x2(:,2) = getCoor(p2,crpPrev)
  x3(:,2) = getCoor(p3,crpPrev)
  x4(:,2) = getCoor(p4,crpPrev)
  x5(:,2) = getCoor(p5,crpPrev)
  x1(:,3) = getCoor(p1,crpPrev2)
  x2(:,3) = getCoor(p2,crpPrev2)
  x3(:,3) = getCoor(p3,crpPrev2)
  x4(:,3) = getCoor(p4,crpPrev2)
  x5(:,3) = getCoor(p5,crpPrev2)



  xmm = (2.0*(x1+x3)+x5)/5.


  ! node 1:

  xm1 = 0.5*(x1+x2)
  xm2 = 0.5*(x1+x4)
  xm3 = 0.5*(x1+x5)
  xb1 = 0.5*(x1+x3)
  xb2 = (x1+x2+x5)/3.
  xb3 = (x1+x4+x5)/3.

  ! side 1-2

  r1 = xmm - xm1
  r2 = xb2 - xm1
  r3 = xb1 - xm1
  r4 = xmm - xm1
  rv = xm1-x1

  vr1(:,1) = xm1(:,1) - xm1(:,2)
  vr1(:,2) = xm1(:,2) - xm1(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb2(:,1) - xb2(:,2)
  vr1(:,2) = vr1(:,2) + xb2(:,2) - xb2(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm1(:,1) - xm1(:,2)
  vr2(:,2) = xm1(:,2) - xm1(:,3)
  vr2(:,1) = vr2(:,1) + xb1(:,1) - xb1(:,2)
  vr2(:,2) = vr2(:,2) + xb1(:,2) - xb1(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0


  call registerSideGCL(p1,p2,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = 2.*nodeVolume

! side 1-4

  r1 = xmm - xm2
  r2 = xb1 - xm2
  r3 = xb3 - xm2
  r4 = xmm - xm2
  rv = xm2-x1

  vr1(:,1) = xm2(:,1) - xm2(:,2)
  vr1(:,2) = xm2(:,2) - xm2(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb1(:,1) - xb1(:,2)
  vr1(:,2) = vr1(:,2) + xb1(:,2) - xb1(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm2(:,1) - xm2(:,2)
  vr2(:,2) = xm2(:,2) - xm2(:,3)
  vr2(:,1) = vr2(:,1) + xb3(:,1) - xb3(:,2)
  vr2(:,2) = vr2(:,2) + xb3(:,2) - xb3(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0

  call registerSideGCL(p1,p4,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

! side 1-5

  r1 = xmm - xm3 
  r2 = xb3 - xm3
  r3 = xb2 - xm3 
  r4 = xmm - xm3  
  rv = xm3-x1

  vr1(:,1) = xm3(:,1) - xm3(:,2)
  vr1(:,2) = xm3(:,2) - xm3(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb3(:,1) - xb3(:,2)
  vr1(:,2) = vr1(:,2) + xb3(:,2) - xb3(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm3(:,1) - xm3(:,2)
  vr2(:,2) = xm3(:,2) - xm3(:,3)
  vr2(:,1) = vr2(:,1) + xb2(:,1) - xb2(:,2)
  vr2(:,2) = vr2(:,2) + xb2(:,2) - xb2(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0

  call registerSideGCL(p1,p5,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

 ! node 3:

  xm1 = 0.5*(x2+x3)
  xm2 = 0.5*(x3+x4)
  xm3 = 0.5*(x3+x5)
  xb1 = 0.5*(x1+x3)
  xb2 = (x2+x3+x5)/3.
  xb3 = (x3+x4+x5)/3.

  ! side 3-2 

  r1 = xmm - xm1
  r2 = xb1 - xm1
  r3 = xb2 - xm1 
  r4 = xmm - xm1 
  rv = xm1-x3

  vr1(:,1) = xm1(:,1) - xm1(:,2)
  vr1(:,2) = xm1(:,2) - xm1(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb1(:,1) - xb1(:,2)
  vr1(:,2) = vr1(:,2) + xb1(:,2) - xb1(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm1(:,1) - xm1(:,2)
  vr2(:,2) = xm1(:,2) - xm1(:,3)
  vr2(:,1) = vr2(:,1) + xb2(:,1) - xb2(:,2)
  vr2(:,2) = vr2(:,2) + xb2(:,2) - xb2(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0


  call registerSideGCL(p3,p2,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

! side 3-4

  r1 = xmm - xm2
  r2 = xb3 - xm2
  r3 = xb1 - xm2
  r4 = xmm - xm2
  rv = xm2-x3

  vr1(:,1) = xm2(:,1) - xm2(:,2)
  vr1(:,2) = xm2(:,2) - xm2(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb3(:,1) - xb3(:,2)
  vr1(:,2) = vr1(:,2) + xb3(:,2) - xb3(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm2(:,1) - xm2(:,2)
  vr2(:,2) = xm2(:,2) - xm2(:,3)
  vr2(:,1) = vr2(:,1) + xb1(:,1) - xb1(:,2)
  vr2(:,2) = vr2(:,2) + xb1(:,2) - xb1(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0


  call registerSideGCL(p3,p4,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

! side 3-5 

  r1 = xmm - xm3 
  r2 = xb2 - xm3 
  r3 = xb3 - xm3
  r4 = xmm - xm3 
  rv = xm3-x3

  vr1(:,1) = xm3(:,1) - xm3(:,2)
  vr1(:,2) = xm3(:,2) - xm3(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb2(:,1) - xb2(:,2)
  vr1(:,2) = vr1(:,2) + xb2(:,2) - xb2(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm3(:,1) - xm3(:,2)
  vr2(:,2) = xm3(:,2) - xm3(:,3)
  vr2(:,1) = vr2(:,1) + xb3(:,1) - xb3(:,2)
  vr2(:,2) = vr2(:,2) + xb3(:,2) - xb3(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0


  call registerSideGCL(p3,p5,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

  ! node 2:

  xm1 = 0.5*(x2+x3)
  xm2 = 0.5*(x1+x2)
  xm3 = 0.5*(x2+x5)
  xb1 = 0.5*(x1+x3)
  xb2 = (x2+x3+x5)/3.
  xb3 = (x1+x2+x5)/3.

  ! side 2-5 

  r1 = xmm - xm3
  r2 = xb3 - xm3
  r3 = xb2 - xm3
  r4 = xmm - xm3
  rv = xm3-x2

  vr1(:,1) = xm3(:,1) - xm3(:,2)
  vr1(:,2) = xm3(:,2) - xm3(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb3(:,1) - xb3(:,2)
  vr1(:,2) = vr1(:,2) + xb3(:,2) - xb3(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm3(:,1) - xm3(:,2)
  vr2(:,2) = xm3(:,2) - xm3(:,3)
  vr2(:,1) = vr2(:,1) + xb2(:,1) - xb2(:,2)
  vr2(:,2) = vr2(:,2) + xb2(:,2) - xb2(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0


  call registerSideGCL(p2,p5,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

 ! node 4:

  xm1 = 0.5*(x1+x4)
  xm2 = 0.5*(x3+x4)
  xm3 = 0.5*(x4+x5)
  xb1 = 0.5*(x1+x3)
  xb2 = (x1+x4+x5)/3.
  xb3 = (x3+x4+x5)/3.

  ! side 4-5 

  r1 = xmm - xm3
  r2 = xb3 - xm3
  r3 = xb2 - xm3
  r4 = xmm - xm3
  rv = xm3-x4

  vr1(:,1) = xm3(:,1) - xm3(:,2)
  vr1(:,2) = xm3(:,2) - xm3(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb3(:,1) - xb3(:,2)
  vr1(:,2) = vr1(:,2) + xb3(:,2) - xb3(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm3(:,1) - xm3(:,2)
  vr2(:,2) = xm3(:,2) - xm3(:,3)
  vr2(:,1) = vr2(:,1) + xb2(:,1) - xb2(:,2)
  vr2(:,2) = vr2(:,2) + xb2(:,2) - xb2(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0

  call registerSideGCL(p4,p5,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

  if(elementVolume(1).le.0) write(*,*) "WARNING: Negative volume for pyramid element number: ",i," (",elementvolume,") ",&
                                    "indexes: ",p1,p2,p3,p4,p5," coordinates: ",x1,x2,x3,x4,x5
 end do
 end subroutine makeCalcDataForPyr
!-------------------------------------------------------------------------
 subroutine makeCalcDataForTet(crp,crpPrev,crpPrev2,tep,grp,brp,doVisualization)
 IMPLICIT NONE

 type(CoordinateRegisterData) :: crp,crpPrev,crpPrev2 ! contains coordinate register for grid
 type(TetrahedralElementData) :: tep ! contains element register for rectangular grid
 type(GridData) ::  grp ! contains grid data  
 type(BoundaryRegisterData) :: brp
 logical :: doVisualization

 integer :: i 
 integer :: p1,p2,p3,p4 
 double precision :: xmm(3,3),xb1(3,3),xb2(3,3),xb3(3,3),xm1(3,3),xm2(3,3),xm3(3,3)  
 double precision :: x1(3,3),x2(3,3),x3(3,3),x4(3,3)
 double precision :: r1(3,3),r2(3,3),r3(3,3),r4(3,3),rv(3,3),vr1(3,2),vr2(3,2)
 double precision :: elementVolume(3),nodeVolume(3)

 type(STEData),pointer :: temp

 do i = 1,tep%numberOfElements ! do for each element
  temp => getTetElement(i,tep) 

  p1=temp%pointIndexes(1)
  p2=temp%pointIndexes(2)
  p3=temp%pointIndexes(3)
  p4=temp%pointIndexes(4)

!  write(*,*) "B: ",p1,p2,p3,p4

  x1(:,1) = getCoor(p1,crp)
  x2(:,1) = getCoor(p2,crp)
  x3(:,1) = getCoor(p3,crp)
  x4(:,1) = getCoor(p4,crp)
  x1(:,2) = getCoor(p1,crpPrev)
  x2(:,2) = getCoor(p2,crpPrev)
  x3(:,2) = getCoor(p3,crpPrev)
  x4(:,2) = getCoor(p4,crpPrev)
  x1(:,3) = getCoor(p1,crpPrev2)
  x2(:,3) = getCoor(p2,crpPrev2)
  x3(:,3) = getCoor(p3,crpPrev2)
  x4(:,3) = getCoor(p4,crpPrev2)


  xmm = 0.25*(x1+x2+x3+x4)

 ! node 1:

  xm1 = 0.5*(x1+x2)
  xm2 = 0.5*(x1+x3)
  xm3 = 0.5*(x1+x4)
  xb1 = (x1+x2+x3)/3.
  xb2 = (x1+x2+x4)/3.
  xb3 = (x1+x3+x4)/3.

  ! side 1-2

  r1 = xmm - xm1
  r2 = xb2 - xm1
  r3 = xb1 - xm1
  r4 = xmm - xm1
  rv = xm1-x1

  vr1(:,1) = xm1(:,1) - xm1(:,2)
  vr1(:,2) = xm1(:,2) - xm1(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb2(:,1) - xb2(:,2)
  vr1(:,2) = vr1(:,2) + xb2(:,2) - xb2(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm1(:,1) - xm1(:,2)
  vr2(:,2) = xm1(:,2) - xm1(:,3)
  vr2(:,1) = vr2(:,1) + xb1(:,1) - xb1(:,2)
  vr2(:,2) = vr2(:,2) + xb1(:,2) - xb1(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0
 
  call registerSideGCL(p1,p2,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = 2.*nodeVolume

! side 1-3

  r1 = xmm - xm2
  r2 = xb1 - xm2
  r3 = xb3 - xm2
  r4 = xmm - xm2
  rv = xm2-x1

  vr1(:,1) = xm2(:,1) - xm2(:,2)
  vr1(:,2) = xm2(:,2) - xm2(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb1(:,1) - xb1(:,2)
  vr1(:,2) = vr1(:,2) + xb1(:,2) - xb1(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm2(:,1) - xm2(:,2)
  vr2(:,2) = xm2(:,2) - xm2(:,3)
  vr2(:,1) = vr2(:,1) + xb3(:,1) - xb3(:,2)
  vr2(:,2) = vr2(:,2) + xb3(:,2) - xb3(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0

  call registerSideGCL(p1,p3,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

! side 1-4

  r1 = xmm - xm3
  r2 = xb3 - xm3
  r3 = xb2 - xm3
  r4 = xmm - xm3
  rv = xm3-x1


  vr1(:,1) = xm3(:,1) - xm3(:,2)
  vr1(:,2) = xm3(:,2) - xm3(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb3(:,1) - xb3(:,2)
  vr1(:,2) = vr1(:,2) + xb3(:,2) - xb3(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm3(:,1) - xm3(:,2)
  vr2(:,2) = xm3(:,2) - xm3(:,3)
  vr2(:,1) = vr2(:,1) + xb2(:,1) - xb2(:,2)
  vr2(:,2) = vr2(:,2) + xb2(:,2) - xb2(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0


  call registerSideGCL(p1,p4,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

 ! node 2:

  xm1 = 0.5*(x2+x3)
  xm2 = 0.5*(x2+x4)
  xm3 = 0.5*(x2+x1)
  xb1 = (x1+x2+x3)/3.
  xb2 = (x2+x3+x4)/3.
  xb3 = (x1+x2+x4)/3.

  ! side 2-3 

  r1 = xmm - xm1
  r2 = xb2 - xm1
  r3 = xb1 - xm1
  r4 = xmm - xm1
  rv = xm1-x2

  vr1(:,1) = xm1(:,1) - xm1(:,2)
  vr1(:,2) = xm1(:,2) - xm1(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb2(:,1) - xb2(:,2)
  vr1(:,2) = vr1(:,2) + xb2(:,2) - xb2(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm1(:,1) - xm1(:,2)
  vr2(:,2) = xm1(:,2) - xm1(:,3)
  vr2(:,1) = vr2(:,1) + xb1(:,1) - xb1(:,2)
  vr2(:,2) = vr2(:,2) + xb1(:,2) - xb1(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0


  call registerSideGCL(p2,p3,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

! side 2-4 

  r1 = xmm - xm2
  r2 = xb3 - xm2
  r3 = xb2 - xm2
  r4 = xmm - xm2
  rv = xm2-x2

  vr1(:,1) = xm2(:,1) - xm2(:,2)
  vr1(:,2) = xm2(:,2) - xm2(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb3(:,1) - xb3(:,2)
  vr1(:,2) = vr1(:,2) + xb3(:,2) - xb3(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm2(:,1) - xm2(:,2)
  vr2(:,2) = xm2(:,2) - xm2(:,3)
  vr2(:,1) = vr2(:,1) + xb2(:,1) - xb2(:,2)
  vr2(:,2) = vr2(:,2) + xb2(:,2) - xb2(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0

  call registerSideGCL(p2,p4,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

! node 3:

  xm1 = 0.5*(x3+x4)
  xm2 = 0.5*(x1+x3)
  xm3 = 0.5*(x2+x3)
  xb1 = (x1+x3+x4)/3.
  xb2 = (x2+x3+x4)/3.
  xb3 = (x1+x2+x3)/3.

  ! side 3-4 

  r1 = xmm - xm1
  r2 = xb2 - xm1
  r3 = xb1 - xm1
  r4 = xmm - xm1
  rv = xm1-x3

  vr1(:,1) = xm1(:,1) - xm1(:,2)
  vr1(:,2) = xm1(:,2) - xm1(:,3)
  vr1(:,1) = vr1(:,1) + xmm(:,1) - xmm(:,2)
  vr1(:,2) = vr1(:,2) + xmm(:,2) - xmm(:,3)
  vr1(:,1) = vr1(:,1) + xb2(:,1) - xb2(:,2)
  vr1(:,2) = vr1(:,2) + xb2(:,2) - xb2(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = xm1(:,1) - xm1(:,2)
  vr2(:,2) = xm1(:,2) - xm1(:,3)
  vr2(:,1) = vr2(:,1) + xb1(:,1) - xb1(:,2)
  vr2(:,2) = vr2(:,2) + xb1(:,2) - xb1(:,3)
  vr2(:,1) = vr2(:,1) + xmm(:,1) - xmm(:,2)
  vr2(:,2) = vr2(:,2) + xmm(:,2) - xmm(:,3)
  vr2 = vr2/3.0


  call registerSideGCL(p3,p4,r1,r2,r3,r4,rv,vr1,vr2,crp,grp,nodeVolume)
  elementVolume = elementVolume + 2.*nodeVolume

  if(elementVolume(1).le.0) write(*,*) "WARNING: Negative volume for tetrahedral element number: ",i," (",elementvolume,") "

 end do
 end subroutine makeCalcDataForTet 
!-------------------------------------------------------------------------
 subroutine doGridAgglomeration(finegrp,coarsegrp,crp,doVisualization,smoothingFactor,smoothingIterations)
 IMPLICIT NONE
  
 type(GridData) :: finegrp,coarsegrp
 type(CoordinateRegisterData) :: crp
 logical :: doVisualization
 real :: smoothingFactor
 integer :: smoothingIterations

 type(BoundaryRegisterData),pointer :: finebrp
 integer :: i

 call mergeNodes(finegrp,coarsegrp) 
 call createSides(finegrp,coarsegrp)   
 finebrp => finegrp%brp
 call cleanUpGridsAfterAgglomeration(finegrp,coarsegrp)
 call prepareDataStructure(coarsegrp%sdp)
 call agglomerateBoundary(finebrp,coarsegrp%brp)
 call cleanUpBoundaryData(finebrp)
 
 write(*,*) "Number of blocks: ",coarsegrp%sdp%numberOfDataBlocks," (",coarsegrp%sdp%dataBlockSize,")" 

 end subroutine doGridAgglomeration
!-------------------------------------------------------------------------
 subroutine mergeNodes(finegrp,coarsegrp)
 IMPLICIT NONE

 type(GridData) :: finegrp,coarsegrp

 integer :: i,j,k,l,m,numberOfBoundarySuperNodes,numberOfSwitches,firstNeighbour
 integer :: allocateStatus,nodeNumber,currentSideNumber,ind1,ind,neighbour,ind2
 integer :: sideNumber,numberOfLonesomeNodes,numberOfSingleNodes,sind1,sind2
 logical :: nodeIsOnBoundary,isSwitched,isLonesome
 integer,pointer :: switchList(:,:),inverseConnectivityList(:)
 integer,pointer :: numberOfSidesConnected(:),availableConnectivity(:),connectivityList(:)
 integer,pointer :: nodeShiftArray(:),numberOfNodesMerged(:)
 logical,pointer :: nodeIsAgglomerated(:),nodeOK(:)
 type(SideDataBlockInstance),pointer :: currentSideBlock
 type(SideData) :: sd,sd2,sds
 real :: averageSideWeight,currentSideWeight
 integer :: numberOfNodeSides,currentSearchSideNumber,searchConnectivity
 logical :: hasPlaced
 integer :: ista,isto,currentConnectivity,currentNode
 integer :: domainNumber,numberOfFineComNodes

 integer,pointer :: maxNodeValue(:,:)

 allocate(nodeIsAgglomerated(finegrp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"
 allocate(coarsegrp%fineToCoarseIndexMapping(finegrp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"
 allocate(availableConnectivity(finegrp%numberOfNodes),stat=allocateStatus) ! contains the connectivity to non-merged nodes
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"
 allocate(connectivityList(finegrp%numberOfNodes),stat=allocateStatus) ! contains the connectivity to non-merged nodes
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"
 allocate(inverseConnectivityList(finegrp%numberOfNodes),stat=allocateStatus) ! contains the connectivity to non-merged nodes
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"

 ! set initial connectivities

 do i=1,finegrp%numberOfNodes
  availableConnectivity(i) = finegrp%sdp%nodeSideArrayRegister(i+1)-finegrp%sdp%nodeSideArrayRegister(i)
  if(availableConnectivity(i).le.0) then
   write(*,*) "ERROR: Something's wrong in mergeNodes - 71"
   STOP
  end if
 end do

 ! set up initial priority list
 ! nodes with largest available connectivities are of highest priority

 connectivityList = -1


 currentConnectivity = -1
 currentNode = 0
 do
  currentConnectivity = currentConnectivity + 1
  do i=1,finegrp%brp%numberOfBoundaryNodes
   if(availableConnectivity(i)==currentConnectivity) then
    currentNode = currentNode + 1
    connectivityList(finegrp%brp%numberOfBoundaryNodes-currentNode+1) = i
   end if
  end do
  if(currentNode==finegrp%brp%numberOfBoundaryNodes) exit
  if(currentConnectivity==500) write(*,*) "WARNING: Current connectivity exceeds 500"
 end do

 currentConnectivity = -1
 currentNode = 0
 do
  currentConnectivity = currentConnectivity + 1
  do i=finegrp%brp%numberOfBoundaryNodes+1,finegrp%numberOfNodes
   if(availableConnectivity(i)==currentConnectivity) then
    currentNode = currentNode + 1
    connectivityList(finegrp%numberOfNodes-currentNode+1) = i
   end if
  end do
  if(currentNode==finegrp%numberOfNodes-finegrp%brp%numberOfBoundaryNodes) exit
  if(currentConnectivity==500) write(*,*) "WARNING: Current connectivity exceeds 500"
 end do

 inverseConnectivityList = 0
 do i=1,finegrp%numberOfNodes
  if(inverseConnectivityList(connectivityList(i))>0) then
   write(*,*) "ERROR: Something's wrong in merge nodes 59"
   STOP
  end if
  inverseConnectivityList(connectivityList(i)) = i
 end do

 nodeIsAgglomerated = .false.
 coarsegrp%fineToCoarseIndexMapping = 0
 ! merge nodes with neighbours
 numberOfBoundarySuperNodes = 0
 nodeNumber = 0
 do i=1,finegrp%numberOfNodes
  ind = connectivityList(i)
  if(.not.nodeIsAgglomerated(ind)) then
   ! find average size of side coefficients connected to node
   numberOfNodeSides = 0
   averageSideWeight = 0.0
   domainNumber = finegrp%sdp%parallelizationColorArray(ind)
   do j=finegrp%sdp%nodeSideArrayRegister(ind),finegrp%sdp%nodeSideArrayRegister(ind+1)-1
    currentSideNumber = finegrp%sdp%nodeSideArray(j)
    sd = getSideNumber(currentSideNumber,finegrp%sdp)
    numberOfNodeSides = numberOfNodeSides + 1
    averageSideWeight = averageSideWeight + sum(sd%sideCoefficients*sd%sideCoefficients)
   end do
   averageSideWeight = sqrt(averageSideWeight)/float(numberOfNodeSides)


   nodeNumber = nodeNumber + 1
   nodeIsAgglomerated(ind) = .true.
   availableConnectivity(ind) = 0
   nodeIsOnBoundary = .false.
   coarsegrp%fineToCoarseIndexMapping(ind) = -nodeNumber  !NOTICE : negative index for seed
   if(ind.le.finegrp%brp%numberOfBoundaryNodes) then
    nodeIsOnBoundary = .true.
    numberOfBoundarySuperNodes = numberOfBoundarySuperNodes + 1
   end if
   ! start merging neighbouring nodes
   do j=finegrp%sdp%nodeSideArrayRegister(ind),finegrp%sdp%nodeSideArrayRegister(ind+1)-1
    currentSideNumber = finegrp%sdp%nodeSideArray(j)
    sd = getSideNumber(currentSideNumber,finegrp%sdp)
    ind1 = sd%sideIndexes(1)
    ind2 = sd%sideIndexes(2)
    if(ind1==ind) then
     neighbour = ind2
    else if(ind2==ind) then
     neighbour = ind1
    else
     STOP "ERROR: Something's wrong in mergeNodes - 1"
    end if
    currentSideWeight = sqrt(sum(sd%sideCoefficients*sd%sideCoefficients))

    if((.not.nodeIsAgglomerated(neighbour)).and.&
       currentSideWeight>coarsegrp%directionalityParameter*averageSideWeight) then 
     ! merge
     nodeIsAgglomerated(neighbour) = .true.
     coarsegrp%fineToCoarseIndexMapping(neighbour) = nodeNumber
     availableConnectivity(neighbour) = 0

     if(.not.nodeIsOnBoundary) then
      if(ind1.le.finegrp%brp%numberOfBoundaryNodes) then
       nodeIsOnBoundary = .true.
       numberOfBoundarySuperNodes = numberOfBoundarySuperNodes + 1
      else if(ind2.le.finegrp%brp%numberOfBoundaryNodes) then
       nodeIsOnBoundary = .true.
       numberOfBoundarySuperNodes = numberOfBoundarySuperNodes + 1
      end if
     end if

     ! now sort out the connectivity list
     do k=finegrp%sdp%nodeSideArrayRegister(neighbour),finegrp%sdp%nodeSideArrayRegister(neighbour+1)-1
      currentSearchSideNumber = finegrp%sdp%nodeSideArray(k)
      sds = getSideNumber(currentSearchSideNumber,finegrp%sdp)
      sind1 = sds%sideIndexes(1)
      sind2 = sds%sideIndexes(2)
      if(.not.nodeIsAgglomerated(sind1)) then
       availableConnectivity(sind1) = availableConnectivity(sind1) - 1
       searchConnectivity = availableConnectivity(sind1)

       if(searchConnectivity<0) then
        write(*,*) "ERROR: Something's worng in merge nodes - 56"
        STOP
       end if

       ! update list

       ista = inverseConnectivityList(sind1)
       if(sind1.le.finegrp%brp%numberOfBoundaryNodes) then
        isto = finegrp%brp%numberOfBoundaryNodes
       else
        isto = finegrp%numberOfNodes
       end if
       do l=ista+1,isto
        if(availableConnectivity(connectivityList(l)).le.searchConnectivity) then
        ! shift
        do m=ista,l-2
         connectivityList(m) = connectivityList(m+1)
         inverseConnectivityList(connectivityList(m)) = m
        end do
        connectivityList(l-1) = sind1
        inverseConnectivityList(sind1) = l-1
        exit
        end if
       end do
      end if
      if(.not.nodeIsAgglomerated(sind2)) then
       availableConnectivity(sind2) = availableConnectivity(sind2) - 1
       searchConnectivity = availableConnectivity(sind2)

       if(searchConnectivity<0) then
        write(*,*) "ERROR: Something's worng in merge nodes - 57"
        STOP
       end if

       ! update list

       ista = inverseConnectivityList(sind2)
       if(sind2.le.finegrp%brp%numberOfBoundaryNodes) then
        isto = finegrp%brp%numberOfBoundaryNodes
       else
        isto = finegrp%numberOfNodes
       end if
       do l=ista+1,isto
        if(availableConnectivity(connectivityList(l))<searchConnectivity) then
        ! shift
        do m=ista,l-2
         connectivityList(m) = connectivityList(m+1)
         inverseConnectivityList(connectivityList(m)) = m
        end do
        connectivityList(l-1) = sind2
        inverseConnectivityList(sind2) = l-1
        exit
        end if
       end do
      end if
     end do
    end if
   end do
  end if
 end do

 do i=1,finegrp%numberOfNodes
  if(.not.nodeIsAgglomerated(i)) then
   write(*,*) "FILLERN: ",i
   PAUSE
  end if
 end do

 deallocate(nodeIsAgglomerated,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes couldn't deallocate"
 deallocate(availableConnectivity,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes couldn't deallocate"
 deallocate(connectivityList,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes couldn't deallocate"
 deallocate(inverseConnectivityList,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes couldn't deallocate"



 ! look for single nodes

 if(.false.) then
 allocate(numberOfNodesMerged(nodeNumber),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"
 numberOfNodesMerged = 0
 do i=1,finegrp%numberOfNodes
  numberOfNodesMerged(coarsegrp%fineToCoarseIndexMapping(i)) = numberOfNodesMerged(coarsegrp%fineToCoarseIndexMapping(i)) + 1
 end do

 deallocate(numberOfNodesMerged,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes couldn't deallocate"
 end if


 ! look for lonesome nodes
 allocate(numberOfSidesConnected(nodeNumber),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"


 ! loop over fine grid sides

 numberOfSidesConnected = 0
 sideNumber = 0
 currentSideBlock => finegrp%sdp%DataBlockChain%first
 do i=1,finegrp%sdp%numberOfDataBlocks
  do j=1,finegrp%sdp%dataBlockSize
   sideNumber = sideNumber + 1
   if(sideNumber>finegrp%sdp%numberOfSides) exit

   sd = currentSideBlock%block(j)
   ind1 = abs(coarsegrp%fineToCoarseIndexMapping(sd%sideIndexes(1)))
   ind2 = abs(coarsegrp%fineToCoarseIndexMapping(sd%sideIndexes(2)))
   if(numberOfSidesConnected(ind1)==0) then
    numberOfSidesConnected(ind1) = -sideNumber
   else if(numberOfSidesConnected(ind1)<0) then
    ! check if side is different from old one
    sd2 = getSideNumber(-numberOfSidesConnected(ind1),finegrp%sdp)
    sind1 = abs(coarsegrp%fineToCoarseIndexMapping(sd2%sideIndexes(1)))
    sind2 = abs(coarsegrp%fineToCoarseIndexMapping(sd2%sideIndexes(2)))
    if(.not.((sind1==ind1.and.sind2==ind2).or.(sind1==ind2.and.sind2==ind1))) then
     numberOfSidesConnected(ind1) = 2
    end if
   end if
   if(numberOfSidesConnected(ind2)==0) then
    numberOfSidesConnected(ind2) = -sideNumber
   else if(numberOfSidesConnected(ind2)<0) then
    sd2 = getSideNumber(-numberOfSidesConnected(ind2),finegrp%sdp)
    sind1 = abs(coarsegrp%fineToCoarseIndexMapping(sd2%sideIndexes(1)))
    sind2 = abs(coarsegrp%fineToCoarseIndexMapping(sd2%sideIndexes(2)))
    if(.not.((sind1==ind1.and.sind2==ind2).or.(sind1==ind2.and.sind2==ind1))) then
     numberOfSidesConnected(ind2) = 2
    end if
   end if
  end do
  currentSideBlock => currentSideBlock%next
 end do

 do i=1,nodeNumber
  if(numberOfSidesConnected(i)==0) then
   write(*,*) "ERROR: Something's wrong in merge nodes - 88: ",i
  end if
 end do

 ! numberOfSidesConnected(ind) now contains the negative side number if node ind is lonesome

 numberOfLonesomeNodes = 0
 do i=1,nodeNumber
  if(numberOfSidesConnected(i)<0) then
   ! lonesome
   numberOfLonesomeNodes = numberOfLonesomeNodes + 1
   sd = getSideNumber(-numberOfSidesConnected(i),finegrp%sdp)
   if(abs(coarsegrp%fineToCoarseIndexMapping(sd%sideIndexes(1)))==i) then
    coarsegrp%fineToCoarseIndexMapping(sd%sideIndexes(1)) = coarsegrp%fineToCoarseIndexMapping(sd%sideIndexes(2))
   else if(abs(coarsegrp%fineToCoarseIndexMapping(sd%sideIndexes(2)))==i) then
    coarsegrp%fineToCoarseIndexMapping(sd%sideIndexes(2)) = coarsegrp%fineToCoarseIndexMapping(sd%sideIndexes(1))
   else
    STOP "ERROR: Something's wrong in merge nodes - 2"
   end if
  else if(numberOfSidesConnected(i)==0) then
   STOP "ERROR: Something's wrong in merge nodes - 3"
  end if
 end do

 ! now shift node numbers down to remove gaps in numbering

 if(numberOfLonesomeNodes>0) then
  allocate(nodeShiftArray(numberOfLonesomeNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"
  j = 0
  do i=1,nodeNumber
   if(numberOfSidesConnected(i)<0) then
    j = j + 1
    nodeShiftArray(j) = i
   end if
  end do

 ! now shift

  do i=1,finegrp%numberOfNodes
   if(abs(coarsegrp%fineToCoarseIndexMapping(i))>nodeShiftArray(1)) then
    do j=1,numberOfLonesomeNodes
     if(abs(coarsegrp%fineToCoarseIndexMapping(i))<nodeShiftArray(j)) exit
    end do
    if(coarsegrp%fineToCoarseIndexMapping(i)>0) then
     coarsegrp%fineToCoarseIndexMapping(i) = coarsegrp%fineToCoarseIndexMapping(i) - j + 1
    else
     coarsegrp%fineToCoarseIndexMapping(i) = coarsegrp%fineToCoarseIndexMapping(i) + j - 1
    end if
   end if
  end do
 end if

 do i=1,numberOfLonesomeNodes
  if(nodeShiftArray(i).le.numberOfBoundarySuperNodes) then
   numberOfBoundarySuperNodes = numberOfBoundarySuperNodes - 1
  else
   exit
  end if
 end do

 deallocate(numberOfSidesConnected,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes couldn't deallocate"

 nodeNumber = nodeNumber - numberOfLonesomeNodes


 coarsegrp%numberOfNodes = nodeNumber
 coarsegrp%brp%numberOfBoundaryNodes = numberOfBoundarySuperNodes
 write(*,*) "Number of nodes in coarse mesh: ",coarsegrp%numberOfNodes
 write(*,*) "Number of boundary nodes in coarse mesh: ",numberOfBoundarySuperNodes
 write(*,*) "Number of lonesome nodes removed: ",numberOfLonesomeNodes

! deallocate stuff that isn't needed anymore

 deallocate(finegrp%sdp%nodeSideArrayRegister,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes couldn't deallocate"
 deallocate(finegrp%sdp%nodeSideArray,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes couldn't deallocate"
 nullify(finegrp%sdp%nodeSideArrayRegister)
 nullify(finegrp%sdp%nodeSideArray)

! reorder to get boundary nodes first

 allocate(switchList(numberOfBoundarySuperNodes,2),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"
 allocate(nodeOK(numberOfBoundarySuperNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"

 switchList = 0
 numberOfSwitches = 0
 nodeOK = .false.
 do i=1,finegrp%brp%numberOfBoundaryNodes
  if(abs(coarsegrp%fineToCoarseIndexMapping(i)) > numberOfBoundarySuperNodes) then
   ! check if node number is already treated
   isSwitched = .false.
   do j=1,numberOfSwitches
    if(switchList(j,1)==coarsegrp%fineToCoarseIndexMapping(i)) then
     isSwitched = .true.
     exit
    end if
   end do
   if(.not.isSwitched) then
    numberOfSwitches = numberOfSwitches + 1
    switchList(numberOfSwitches,1) = coarsegrp%fineToCoarseIndexMapping(i)
   end if
  else
   nodeOK(abs(coarsegrp%fineToCoarseIndexMapping(i))) = .true.
  end if
 end do

 ! find nodes to switch with

 do i=1,numberOfSwitches
  do j=1,numberOfBoundarySuperNodes
   if(.not.nodeOK(j)) then
    switchList(i,2) = j
    nodeOK(j) = .true.
    exit
   end if
  end do
 end do

 ! do switching

 do i=1,numberOfSwitches
  do j=1,finegrp%numberOfNodes
   if(switchList(i,1) == coarsegrp%fineToCoarseIndexMapping(j)) then
    coarsegrp%fineToCoarseIndexMapping(j) = switchList(i,2)
   else if(switchList(i,2) == coarsegrp%fineToCoarseIndexMapping(j)) then
    coarsegrp%fineToCoarseIndexMapping(j) = switchList(i,1)
   end if
  end do
 end do

 deallocate(switchList,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes couldn't deallocate"
 deallocate(nodeOK,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes couldn't deallocate"


! make seed node list

 allocate(coarsegrp%seedNodeRegister(coarsegrp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"

 if(coarsegrp%gridNumber==2) then
  do i=1,finegrp%numberOfNodes
   if(coarsegrp%fineToCoarseIndexMapping(i)<0) then
    coarsegrp%fineToCoarseIndexMapping(i) = -coarsegrp%fineToCoarseIndexMapping(i)
    coarsegrp%seedNodeRegister(coarsegrp%fineToCoarseIndexMapping(i)) = i
   end if
  end do
 else
  do i=1,finegrp%numberOfNodes
   if(coarsegrp%fineToCoarseIndexMapping(i)<0) then
    coarsegrp%fineToCoarseIndexMapping(i) = -coarsegrp%fineToCoarseIndexMapping(i)
    coarsegrp%seedNodeRegister(coarsegrp%fineToCoarseIndexMapping(i)) = finegrp%seedNodeRegister(i)
   end if
  end do
 end if


! set up control volume array

 allocate(coarsegrp%controlVolumeNodeSize(coarsegrp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"

 coarsegrp%controlVolumeNodeSize = 0.0

 do i=1,finegrp%numberOfNodes
   nodeNumber = coarsegrp%fineToCoarseIndexMapping(i)
   coarsegrp%controlVolumeNodeSize(nodeNumber) = coarsegrp%controlVolumeNodeSize(nodeNumber) +&
                                                  finegrp%controlVolumeNodeSize(i)
 end do


 ! make parallelization data for coarse mesh

 allocate(coarsegrp%sdp,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"
 allocate(coarsegrp%sdp%parallelizationColorArray(coarsegrp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"
 coarsegrp%sdp%parallelizationColorArray = 0


! coarse node belongs to domain of which most fine mesh nodes associated with it belong

 allocate(maxNodeValue(coarsegrp%numberOfNodes,2),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"

 maxNodeValue = 0
 write(*,*) "AA: ",finegrp%numberOfGridDomains
 do i=1,finegrp%numberOfGridDomains
  maxNodeValue(:,2) = 0
  do j=1,finegrp%numberOfNodes
   if(finegrp%sdp%parallelizationColorArray(j)==i) then
    maxNodeValue(coarsegrp%fineToCoarseIndexMapping(j),2) = maxNodeValue(coarsegrp%fineToCoarseIndexMapping(j),2) + 1
   end if
  end do

  do j=1,coarsegrp%numberOfNodes
   if(maxNodeValue(j,2)>maxNodeValue(j,1)) then
    maxNodeValue(j,1) = maxNodeValue(j,2)
    coarsegrp%sdp%parallelizationColorArray(j) = i
   end if
  end do
 end do

 deallocate(maxNodeValue,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes couldn't deallocate"

! make register of fine grid nodes that are associated with coarse grid nodes
! in a different parallel domain

! first count

 numberOfFineComNodes = 0

 do i=1,finegrp%numberOfNodes
  ind2 = coarsegrp%fineToCoarseIndexMapping(i) ! coarse node number
  if(finegrp%sdp%parallelizationColorArray(i).ne.coarsegrp%sdp%parallelizationColorArray(ind2)) then
   numberOfFineComNodes = numberOfFineComNodes + 1
  end if
 end do

 coarsegrp%numberOfFineComNodes = numberOfFineComNodes

 if(numberOfFineComNodes>0) then
  allocate(coarsegrp%fineComNodes(numberOfFineComNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"

  numberOfFineComNodes = 0
  do i=1,finegrp%numberOfNodes
   ind2 = coarsegrp%fineToCoarseIndexMapping(i) ! coarse node number
   if(finegrp%sdp%parallelizationColorArray(i).ne.coarsegrp%sdp%parallelizationColorArray(ind2)) then
    numberOfFineComNodes = numberOfFineComNodes + 1
    coarsegrp%fineComNodes(numberOfFineComNodes) = i
   end if
  end do
 else
  nullify(coarsegrp%fineComNodes)
 end if

 end subroutine mergeNodes
!-------------------------------------------------------------------------
 subroutine mergeNodes4(finegrp,coarsegrp)
 IMPLICIT NONE

 ! gives one number to each fine grid node indicating the coarse grid node it belongs to

 type(GridData) :: finegrp,coarsegrp

 integer :: i,j,numberOfBoundarySuperNodes,numberOfSwitches,firstNeighbour
 integer :: allocateStatus,nodeNumber,currentSideNumber,ind1,ind,neighbour,ind2
 integer :: sideNumber,numberOfLonesomeNodes,numberOfSingleNodes,sind1,sind2
 logical :: nodeIsOnBoundary,isSwitched,isLonesome
 integer,pointer :: switchList(:,:)
 integer,pointer :: numberOfSidesConnected(:)
 integer,pointer :: nodeShiftArray(:),numberOfNodesMerged(:)
 logical,pointer :: isBoundaryNode(:),nodeOK(:)
 type(SideDataBlockInstance),pointer :: currentSideBlock
 type(SideData) :: sd,sd2
 real :: averageSideWeight,currentSideWeight,minimumSideWeight,maximumSideWeight
 integer :: numberOfNodeSides,maximumAvailableNeighbour,domainNumber
 logical :: passedDirectionalTest
 character :: gr,prevgr
 integer :: numberOfFinestNodes
 integer,pointer :: prevMapping(:)
 integer :: buffi
 logical,pointer :: hasSeed(:)


 allocate(coarsegrp%fineToCoarseIndexMapping(finegrp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"

 coarsegrp%fineToCoarseIndexMapping = 0

 write(gr,'(I1)') coarsegrp%gridNumber

 open(5,file='aggl_'//gr//'.reg',status="old")
 if(coarsegrp%gridNumber==2) then 
  do i=1,finegrp%numberOfNodes
   read(5,*) coarsegrp%fineToCoarseIndexMapping(i) 
  end do 
 else
  write(prevgr,'(I1)') finegrp%gridNumber
  open(7,file='aggreg_'//prevgr//'.reg',form="unformatted",status="unknown")
  read(7) numberOfFinestNodes
  allocate(prevMapping(numberOfFinestNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"
  read(7) (prevMapping(i),i=1,numberOfFinestNodes)
  close(7)

  do i=1,numberOfFinestNodes
   read(5,*) buffi
   coarsegrp%fineToCoarseIndexMapping(abs(prevMapping(i))) = buffi
  end do

 end if
 close(5)

 do i=1,finegrp%numberOfNodes
  if(coarsegrp%fineToCoarseIndexMapping(i)==0) write(*,*) "pottit: ",i
 end do 

! find number of nodes in coarse mesh

 coarsegrp%numberOfNodes = 0
 do i=1,finegrp%numberOfNodes
  if(abs(coarsegrp%fineToCoarseIndexMapping(i))>coarsegrp%numberOfNodes)&
       coarsegrp%numberOfNodes = abs(coarsegrp%fineToCoarseIndexMapping(i))
 end do

 write(*,*) "Number of nodes in mesh: ",coarsegrp%numberOfNodes



! find number of boundary nodes in coarse mesh

 allocate(isBoundaryNode(coarsegrp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"

 isBoundaryNode = .false.

 do i=1,finegrp%brp%numberOfBoundaryNodes
  isBoundaryNode(abs(coarsegrp%fineToCoarseIndexMapping(i))) = .true.
 end do

 numberOfBoundarySuperNodes = 0
 do i=1,coarsegrp%numberOfNodes
  if(isBoundaryNode(i)) numberOfBoundarySuperNodes = numberOfBoundarySuperNodes + 1
 end do
 coarsegrp%brp%numberOfBoundaryNodes = numberOfBoundarySuperNodes

 write(*,*) "Number of boundary supernodes: ",numberOfBoundarySuperNodes

 deallocate(isBoundaryNode,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes couldn't deallocate"

! reorder to get boundary nodes first

 allocate(switchList(numberOfBoundarySuperNodes,2),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"
 allocate(nodeOK(numberOfBoundarySuperNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"

 switchList = 0
 numberOfSwitches = 0
 nodeOK = .false.
 do i=1,finegrp%brp%numberOfBoundaryNodes
  if(abs(coarsegrp%fineToCoarseIndexMapping(i)) > numberOfBoundarySuperNodes) then
   ! check if node number is already treated
   isSwitched = .false.
   do j=1,numberOfSwitches
    if(switchList(j,1)==abs(coarsegrp%fineToCoarseIndexMapping(i))) then
     isSwitched = .true.
     exit
    end if
   end do
   if(.not.isSwitched) then
    numberOfSwitches = numberOfSwitches + 1
    switchList(numberOfSwitches,1) = abs(coarsegrp%fineToCoarseIndexMapping(i))
   end if
  else
   nodeOK(abs(coarsegrp%fineToCoarseIndexMapping(i))) = .true.
  end if
 end do

 ! find nodes to switch with

 do i=1,numberOfSwitches
  do j=1,numberOfBoundarySuperNodes
   if(.not.nodeOK(j)) then
    switchList(i,2) = j
    nodeOK(j) = .true.
    exit
   end if
  end do
 end do

 ! do switching

 do i=1,numberOfSwitches
  do j=1,finegrp%numberOfNodes
   ind = coarsegrp%fineToCoarseIndexMapping(j)
   if(switchList(i,1) == abs(ind)) then 
    coarsegrp%fineToCoarseIndexMapping(j) = switchList(i,2)
    if(ind<0) coarsegrp%fineToCoarseIndexMapping(j) = -coarsegrp%fineToCoarseIndexMapping(j)
   else if(switchList(i,2) == abs(ind)) then
    coarsegrp%fineToCoarseIndexMapping(j) = switchList(i,1)
    if(ind<0) coarsegrp%fineToCoarseIndexMapping(j) = -coarsegrp%fineToCoarseIndexMapping(j)
   end if
  end do
 end do

 deallocate(switchList,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes couldn't deallocate"
 deallocate(nodeOK,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes couldn't deallocate"

! make seed node list

 allocate(coarsegrp%seedNodeRegister(coarsegrp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"
 coarsegrp%seedNodeRegister = 0

 if(coarsegrp%gridNumber==2) then
  do i=1,finegrp%numberOfNodes
   if(coarsegrp%fineToCoarseIndexMapping(i)<0) then
    coarsegrp%fineToCoarseIndexMapping(i) = -coarsegrp%fineToCoarseIndexMapping(i)
    coarsegrp%seedNodeRegister(coarsegrp%fineToCoarseIndexMapping(i)) = i
   end if
  end do
 else
  do i=1,finegrp%numberOfNodes
   if(coarsegrp%fineToCoarseIndexMapping(i)<0) then
    coarsegrp%fineToCoarseIndexMapping(i) = -coarsegrp%fineToCoarseIndexMapping(i)
    coarsegrp%seedNodeRegister(coarsegrp%fineToCoarseIndexMapping(i)) = finegrp%seedNodeRegister(i)
   end if
  end do
 end if

! set up control volume array

 allocate(coarsegrp%controlVolumeNodeSize(coarsegrp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"

 coarsegrp%controlVolumeNodeSize = 0.0

 do i=1,finegrp%numberOfNodes
   nodeNumber = coarsegrp%fineToCoarseIndexMapping(i)
   coarsegrp%controlVolumeNodeSize(nodeNumber) = coarsegrp%controlVolumeNodeSize(nodeNumber) +&
                                                  finegrp%controlVolumeNodeSize(i)
 end do


 ! make parallelization data for coarse mesh

 allocate(coarsegrp%sdp,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"
 allocate(coarsegrp%sdp%parallelizationColorArray(coarsegrp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"
 coarsegrp%sdp%parallelizationColorArray = 0


 do i=1,finegrp%numberOfNodes
  domainNumber = finegrp%sdp%parallelizationColorArray(i)
  ind2 = coarsegrp%fineToCoarseIndexMapping(i)
  if(coarsegrp%sdp%parallelizationColorArray(ind2)==0) then
   coarsegrp%sdp%parallelizationColorArray(ind2) = domainNumber
  else if(coarsegrp%sdp%parallelizationColorArray(ind2).ne.domainNumber) then
   write(*,*) i,ind2,coarsegrp%sdp%parallelizationColorArray(ind2),domainNumber
   STOP "ERROR: Something's wrong in merge nodes - 31"
  end if
 end do

 if(coarsegrp%gridNumber==2) then 
  open(7,file='aggreg_'//gr//'.reg',form="unformatted",status="unknown")
  write(7) finegrp%numberOfNodes
  write(7) (coarsegrp%fineToCoarseIndexMapping(i),i=1,finegrp%numberOfNodes)
  close(7)
 else

  open(7,file='aggreg_'//gr//'.reg',form="unformatted",status="unknown")
  write(7) numberOfFinestNodes
  write(7) (coarsegrp%fineToCoarseIndexMapping(prevMapping(i)),i=1,numberOfFinestNodes)
  close(7)
  deallocate(prevMapping,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: mergeNodes couldn't deallocate"
 end if

 end subroutine mergeNodes4
!-------------------------------------------------------------------------
 subroutine mergeNodes3(finegrp,coarsegrp)
 IMPLICIT NONE

 ! gives one number to each fine grid node indicating the coarse grid node it belongs to

 type(GridData) :: finegrp,coarsegrp

 integer :: i,j,numberOfBoundarySuperNodes,numberOfSwitches,firstNeighbour
 integer :: allocateStatus,nodeNumber,currentSideNumber,ind1,ind,neighbour,ind2
 integer :: sideNumber,numberOfLonesomeNodes,numberOfSingleNodes,sind1,sind2
 logical :: nodeIsOnBoundary,isSwitched,isLonesome
 integer,pointer :: switchList(:,:)
 integer,pointer :: numberOfSidesConnected(:)
 integer,pointer :: nodeShiftArray(:),numberOfNodesMerged(:)
 logical,pointer :: nodeIsAgglomerated(:),nodeOK(:)
 type(SideDataBlockInstance),pointer :: currentSideBlock
 type(SideData) :: sd,sd2
 real :: averageSideWeight,currentSideWeight,minimumSideWeight,maximumSideWeight
 integer :: numberOfNodeSides,maximumAvailableNeighbour,domainNumber
 logical :: passedDirectionalTest

 allocate(nodeIsAgglomerated(finegrp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"
 allocate(coarsegrp%fineToCoarseIndexMapping(finegrp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"


 nodeIsAgglomerated = .false.
 coarsegrp%fineToCoarseIndexMapping = 0
 ! merge nodes with neighbours
 numberOfBoundarySuperNodes = 0
 nodeNumber = 0

 do i=1,finegrp%numberOfNodes
  if(.not.nodeIsAgglomerated(i)) then
   ! find average size of side coefficients connected to node
   numberOfNodeSides = 0
   averageSideWeight = 0.0
   minimumSideWeight = 99999999.9
   maximumSideWeight = -1.0
   maximumAvailableNeighbour = 0
   domainNumber = finegrp%sdp%parallelizationColorArray(i)
   do j=finegrp%sdp%nodeSideArrayRegister(i),finegrp%sdp%nodeSideArrayRegister(i+1)-1
    currentSideNumber = finegrp%sdp%nodeSideArray(j)
    sd = getSideNumber(currentSideNumber,finegrp%sdp)
    numberOfNodeSides = numberOfNodeSides + 1
    currentSideWeight = sum(sd%sideCoefficients*sd%sideCoefficients)
    if(currentSideWeight>maximumSideWeight) then
     maximumSideWeight = currentSideWeight
     sind1 = sd%sideIndexes(1)
     sind2 = sd%sideIndexes(2)
     if(sind1==i) then
      if(.not.nodeIsAgglomerated(sind2)) then
       maximumAvailableNeighbour = sind2
      else
       ! optional
       maximumAvailableNeighbour = 0
      end if
     else if(sind2==i) then
      if(.not.nodeIsAgglomerated(sind1)) then
       maximumAvailableNeighbour = sind1
      else
       maximumAvailableNeighbour = 0
      end if
     else
      STOP "ERROR: Something's wrong in mergeNodes - 7"
     end if
    end if
    if(currentSideWeight<minimumSideWeight) then
      minimumSideWeight = currentSideWeight
    end if
    averageSideWeight = averageSideWeight + currentSideWeight
   end do
   averageSideWeight = sqrt(averageSideWeight)/float(numberOfNodeSides)

   nodeNumber = nodeNumber + 1
   nodeIsAgglomerated(i) = .true.
   nodeIsOnBoundary = .false.
   coarsegrp%fineToCoarseIndexMapping(i) = -nodeNumber  !NOTICE : negative index for seed

   if(i.le.finegrp%brp%numberOfBoundaryNodes) then
    nodeIsOnBoundary = .true.
    numberOfBoundarySuperNodes = numberOfBoundarySuperNodes + 1
   end if
   do j=finegrp%sdp%nodeSideArrayRegister(i),finegrp%sdp%nodeSideArrayRegister(i+1)-1
    currentSideNumber = finegrp%sdp%nodeSideArray(j)
    sd = getSideNumber(currentSideNumber,finegrp%sdp)
    ind1 = sd%sideIndexes(1)
    ind2 = sd%sideIndexes(2)
    if(ind1==i) then
     neighbour = ind2
    else if(ind2==i) then
     neighbour = ind1
    else
     STOP "ERROR: Something's wrong in mergeNodes - 1"
    end if
    currentSideWeight = sqrt(sum(sd%sideCoefficients*sd%sideCoefficients))

    if(maximumSideWeight>4.0*minimumSideWeight) then
     if(currentSideWeight>coarsegrp%directionalityParameter*minimumSideWeight) then
      passedDirectionalTest = .true.
     else
      passedDirectionalTest = .false.
     end if
    else
     passedDirectionalTest = .true.
    end if

   if(.not.nodeIsAgglomerated(neighbour).and.passedDirectionalTest.and.&
      domainNumber==finegrp%sdp%parallelizationColorArray(neighbour)) then
     ! merge
     nodeIsAgglomerated(neighbour) = .true.
     coarsegrp%fineToCoarseIndexMapping(neighbour) = nodeNumber

     if(.not.nodeIsOnBoundary) then
      if(ind1.le.finegrp%brp%numberOfBoundaryNodes) then
       nodeIsOnBoundary = .true.
       numberOfBoundarySuperNodes = numberOfBoundarySuperNodes + 1
      else if(ind2.le.finegrp%brp%numberOfBoundaryNodes) then
       nodeIsOnBoundary = .true.
       numberOfBoundarySuperNodes = numberOfBoundarySuperNodes + 1
      end if
     end if
    end if
   end do

  ! now merge strongest connected neighbour's neighbours
   if(.false.) then 
   if(maximumAvailableNeighbour>0) then
    if(nodeIsAgglomerated(maximumAvailableNeighbour)) then
     do j=finegrp%sdp%nodeSideArrayRegister(maximumAvailableNeighbour),&
           finegrp%sdp%nodeSideArrayRegister(maximumAvailableNeighbour+1)-1
       currentSideNumber = finegrp%sdp%nodeSideArray(j)
      sd = getSideNumber(currentSideNumber,finegrp%sdp)
      ind1 = sd%sideIndexes(1)
      ind2 = sd%sideIndexes(2)
      if(ind1==maximumAvailableNeighbour) then
       neighbour = ind2
      else if(ind2==maximumAvailableNeighbour) then
       neighbour = ind1
      else
       STOP "ERROR: Something's wrong in mergeNodes - 1"
      end if
!      if(.not.nodeIsAgglomerated(neighbour)) then
      if(.not.nodeIsAgglomerated(neighbour).and.&
      domainNumber==finegrp%sdp%parallelizationColorArray(neighbour)) then
       ! merge
       nodeIsAgglomerated(neighbour) = .true.
       coarsegrp%fineToCoarseIndexMapping(neighbour) = nodeNumber
       if(.not.nodeIsOnBoundary) then
        if(ind1.le.finegrp%brp%numberOfBoundaryNodes) then
         nodeIsOnBoundary = .true.
         numberOfBoundarySuperNodes = numberOfBoundarySuperNodes + 1
        else if(ind2.le.finegrp%brp%numberOfBoundaryNodes) then
         nodeIsOnBoundary = .true.
         numberOfBoundarySuperNodes = numberOfBoundarySuperNodes + 1
        end if
       end if
      end if
     end do
    end if
   end if
   end if
  end if
 end do

 do i=1,finegrp%numberOfNodes
  if(.not.nodeIsAgglomerated(i)) write(*,*) "FILLERN: ",i
 end do


 deallocate(nodeIsAgglomerated,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes couldn't deallocate"

 ! look for single nodes

 if(.false.) then
 allocate(numberOfNodesMerged(nodeNumber),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"
 numberOfNodesMerged = 0
 do i=1,finegrp%numberOfNodes
  numberOfNodesMerged(coarsegrp%fineToCoarseIndexMapping(i)) = numberOfNodesMerged(coarsegrp%fineToCoarseIndexMapping(i)) + 1
 end do

 deallocate(numberOfNodesMerged,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes couldn't deallocate"
 end if


 ! look for lonesome nodes
 allocate(numberOfSidesConnected(nodeNumber),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"


 ! loop over fine grid sides

 numberOfSidesConnected = 0
 sideNumber = 0
 currentSideBlock => finegrp%sdp%DataBlockChain%first
 do i=1,finegrp%sdp%numberOfDataBlocks
  do j=1,finegrp%sdp%dataBlockSize
   sideNumber = sideNumber + 1
   if(sideNumber>finegrp%sdp%numberOfSides) exit

   sd = currentSideBlock%block(j)
   ind1 = abs(coarsegrp%fineToCoarseIndexMapping(sd%sideIndexes(1)))
   ind2 = abs(coarsegrp%fineToCoarseIndexMapping(sd%sideIndexes(2)))
   if(numberOfSidesConnected(ind1)==0) then
    numberOfSidesConnected(ind1) = -sideNumber
   else if(numberOfSidesConnected(ind1)<0) then
    ! check if side is different from old one
    sd2 = getSideNumber(-numberOfSidesConnected(ind1),finegrp%sdp)
    sind1 = abs(coarsegrp%fineToCoarseIndexMapping(sd2%sideIndexes(1)))
    sind2 = abs(coarsegrp%fineToCoarseIndexMapping(sd2%sideIndexes(2)))
    if(.not.((sind1==ind1.and.sind2==ind2).or.(sind1==ind2.and.sind2==ind1))) then
     numberOfSidesConnected(ind1) = 2
    end if
   end if
   if(numberOfSidesConnected(ind2)==0) then
    numberOfSidesConnected(ind2) = -sideNumber
   else if(numberOfSidesConnected(ind2)<0) then
    sd2 = getSideNumber(-numberOfSidesConnected(ind2),finegrp%sdp)
    sind1 = abs(coarsegrp%fineToCoarseIndexMapping(sd2%sideIndexes(1)))
    sind2 = abs(coarsegrp%fineToCoarseIndexMapping(sd2%sideIndexes(2)))
    if(.not.((sind1==ind1.and.sind2==ind2).or.(sind1==ind2.and.sind2==ind1))) then
     numberOfSidesConnected(ind2) = 2
    end if
   end if
  end do
  currentSideBlock => currentSideBlock%next
 end do

 ! numberOfSidesConnected(ind) now contains the negative side number if node ind is lonesome

 numberOfLonesomeNodes = 0
 do i=1,nodeNumber
  if(numberOfSidesConnected(i)<0) then
   ! lonesome
   numberOfLonesomeNodes = numberOfLonesomeNodes + 1
   sd = getSideNumber(-numberOfSidesConnected(i),finegrp%sdp)
   if(abs(coarsegrp%fineToCoarseIndexMapping(sd%sideIndexes(1)))==i) then
    coarsegrp%fineToCoarseIndexMapping(sd%sideIndexes(1)) = coarsegrp%fineToCoarseIndexMapping(sd%sideIndexes(2))
   else if(abs(coarsegrp%fineToCoarseIndexMapping(sd%sideIndexes(2)))==i) then
    coarsegrp%fineToCoarseIndexMapping(sd%sideIndexes(2)) = coarsegrp%fineToCoarseIndexMapping(sd%sideIndexes(1))
   else
    STOP "ERROR: Something's wrong in merge nodes - 2"
   end if
  else if(numberOfSidesConnected(i)==0) then
   STOP "ERROR: Something's wrong in merge nodes - 3"
  end if
 end do

 ! now shift node numbers down to remove gaps in numbering

 if(numberOfLonesomeNodes>0) then
  allocate(nodeShiftArray(numberOfLonesomeNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"
  j = 0
  do i=1,nodeNumber
   if(numberOfSidesConnected(i)<0) then
    j = j + 1
    nodeShiftArray(j) = i
   end if
  end do

 ! now shift

  do i=1,finegrp%numberOfNodes
   if(abs(coarsegrp%fineToCoarseIndexMapping(i))>nodeShiftArray(1)) then
    do j=1,numberOfLonesomeNodes
     if(abs(coarsegrp%fineToCoarseIndexMapping(i))<nodeShiftArray(j)) exit
    end do
    if(coarsegrp%fineToCoarseIndexMapping(i)>0) then
     coarsegrp%fineToCoarseIndexMapping(i) = coarsegrp%fineToCoarseIndexMapping(i) - j + 1
    else
     coarsegrp%fineToCoarseIndexMapping(i) = coarsegrp%fineToCoarseIndexMapping(i) + j - 1
    end if
   end if
  end do
 end if

 do i=1,numberOfLonesomeNodes
  if(nodeShiftArray(i).le.numberOfBoundarySuperNodes) then
   numberOfBoundarySuperNodes = numberOfBoundarySuperNodes - 1
  else
   exit
  end if
 end do

 deallocate(numberOfSidesConnected,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes couldn't deallocate"

 nodeNumber = nodeNumber - numberOfLonesomeNodes


 coarsegrp%numberOfNodes = nodeNumber
 coarsegrp%brp%numberOfBoundaryNodes = numberOfBoundarySuperNodes
 write(*,*) "Number of nodes in coarse mesh: ",coarsegrp%numberOfNodes
 write(*,*) "Number of boundary nodes in coarse mesh: ",numberOfBoundarySuperNodes
 write(*,*) "Number of lonesome nodes removed: ",numberOfLonesomeNodes

! deallocate stuff that isn't needed anymore

 deallocate(finegrp%sdp%nodeSideArrayRegister,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes couldn't deallocate"
 deallocate(finegrp%sdp%nodeSideArray,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes couldn't deallocate"
 nullify(finegrp%sdp%nodeSideArrayRegister)
 nullify(finegrp%sdp%nodeSideArray)

! reorder to get boundary nodes first

 allocate(switchList(numberOfBoundarySuperNodes,2),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"
 allocate(nodeOK(numberOfBoundarySuperNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"

 switchList = 0
 numberOfSwitches = 0
 nodeOK = .false.
 do i=1,finegrp%brp%numberOfBoundaryNodes
  if(abs(coarsegrp%fineToCoarseIndexMapping(i)) > numberOfBoundarySuperNodes) then
   ! check if node number is already treated
   isSwitched = .false.
   do j=1,numberOfSwitches
    if(switchList(j,1)==abs(coarsegrp%fineToCoarseIndexMapping(i))) then
     isSwitched = .true.
     exit
    end if
   end do
   if(.not.isSwitched) then
    numberOfSwitches = numberOfSwitches + 1
    switchList(numberOfSwitches,1) = abs(coarsegrp%fineToCoarseIndexMapping(i))
   end if
  else
   nodeOK(abs(coarsegrp%fineToCoarseIndexMapping(i))) = .true.
  end if
 end do

 ! find nodes to switch with

 do i=1,numberOfSwitches
  do j=1,numberOfBoundarySuperNodes
   if(.not.nodeOK(j)) then
    switchList(i,2) = j
    nodeOK(j) = .true.
    exit
   end if
  end do
 end do

 ! do switching

 do i=1,numberOfSwitches
  do j=1,finegrp%numberOfNodes
   ind = coarsegrp%fineToCoarseIndexMapping(j)
   if(switchList(i,1) == abs(ind)) then
    coarsegrp%fineToCoarseIndexMapping(j) = switchList(i,2)
    if(ind<0) coarsegrp%fineToCoarseIndexMapping(j) = -coarsegrp%fineToCoarseIndexMapping(j)
   else if(switchList(i,2) == abs(ind)) then
    coarsegrp%fineToCoarseIndexMapping(j) = switchList(i,1)
    if(ind<0) coarsegrp%fineToCoarseIndexMapping(j) = -coarsegrp%fineToCoarseIndexMapping(j)
   end if
  end do
 end do

 deallocate(switchList,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes couldn't deallocate"
 deallocate(nodeOK,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes couldn't deallocate"

! make seed node list

 allocate(coarsegrp%seedNodeRegister(coarsegrp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"

 if(coarsegrp%gridNumber==2) then
  do i=1,finegrp%numberOfNodes
   if(coarsegrp%fineToCoarseIndexMapping(i)<0) then
    coarsegrp%fineToCoarseIndexMapping(i) = -coarsegrp%fineToCoarseIndexMapping(i)
    coarsegrp%seedNodeRegister(coarsegrp%fineToCoarseIndexMapping(i)) = i
   end if
  end do
 else
  do i=1,finegrp%numberOfNodes
   if(coarsegrp%fineToCoarseIndexMapping(i)<0) then
    coarsegrp%fineToCoarseIndexMapping(i) = -coarsegrp%fineToCoarseIndexMapping(i)
    coarsegrp%seedNodeRegister(coarsegrp%fineToCoarseIndexMapping(i)) = finegrp%seedNodeRegister(i)
   end if
  end do
 end if

! set up control volume array

 allocate(coarsegrp%controlVolumeNodeSize(coarsegrp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"

 coarsegrp%controlVolumeNodeSize = 0.0

 do i=1,finegrp%numberOfNodes
   nodeNumber = coarsegrp%fineToCoarseIndexMapping(i)
   coarsegrp%controlVolumeNodeSize(nodeNumber) = coarsegrp%controlVolumeNodeSize(nodeNumber) +&
                                                  finegrp%controlVolumeNodeSize(i)
 end do


 ! make parallelization data for coarse mesh

 allocate(coarsegrp%sdp,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"
 allocate(coarsegrp%sdp%parallelizationColorArray(coarsegrp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: mergeNodes out of memory"
 coarsegrp%sdp%parallelizationColorArray = 0


 do i=1,finegrp%numberOfNodes
  domainNumber = finegrp%sdp%parallelizationColorArray(i)
  ind2 = coarsegrp%fineToCoarseIndexMapping(i)
  if(coarsegrp%sdp%parallelizationColorArray(ind2)==0) then
   coarsegrp%sdp%parallelizationColorArray(ind2) = domainNumber
  else if(coarsegrp%sdp%parallelizationColorArray(ind2).ne.domainNumber) then
   write(*,*) i,ind2,coarsegrp%sdp%parallelizationColorArray(ind2),domainNumber
   STOP "ERROR: Something's wrong in merge nodes - 31"
  end if
 end do
 end subroutine mergeNodes3
!-------------------------------------------------------------------------
 subroutine createSides(finegrp,coarsegrp)
 IMPLICIT NONE

 type(GridData) :: finegrp,coarsegrp 

 integer :: i,j,allocateStatus,sideNumber,blockSize,numberOfSidesCreated
 type(SideDataBlockInstance),pointer :: currentSideBlock
 type(SideData) :: sd
 integer :: sideCreated 

 blockSize = max(coarsegrp%numberOfNodes*2,100) 
 call setUpSideRegister(coarsegrp%sdp,coarsegrp%numberOfNodes,blockSize,blockSize)

 ! loop over fine grid sides

  sideNumber = 0
  numberOfSidesCreated = 0
  do i=1,finegrp%sdp%numberOfDataBlocks 
   currentSideBlock => finegrp%sdp%DataBlockChain%first
   do j=1,finegrp%sdp%dataBlockSize
    sideNumber = sideNumber + 1
    if(sideNumber>finegrp%sdp%numberOfSides) exit

    ! start creating supersides
    
    sd = currentSideBlock%block(j)
   
    sd%sideIndexes(1) = coarsegrp%fineToCoarseIndexMapping(sd%sideIndexes(1))   
    sd%sideIndexes(2) = coarsegrp%fineToCoarseIndexMapping(sd%sideIndexes(2))   

    if(sd%sideIndexes(1).ne.sd%sideIndexes(2)) then 
     ! if the side indexes are the same the side is internal and is not
     ! taken into account
     
     ! add side to register
     
     sideCreated = registerSideData(sd,coarsegrp%sdp)
     if(sideCreated==0) then
      numberOfSidesCreated = numberOfSidesCreated + 1
      sideCreated = numberOfSidesCreated
     end if
    end if
   
   end do
   finegrp%sdp%DataBlockChain%first => currentSideBlock%next
   ! remove this side block - not needed anymore 
   deallocate(currentSideBlock%block,stat=allocateStatus)
   if (allocateStatus /= 0) STOP "ERROR: mergeNodes couldn't deallocate" 
  end do

  write(*,*) "Number of sides in coarse mesh: ",coarsegrp%sdp%numberOfSides
  coarsegrp%numberOfSides = coarsegrp%sdp%numberOfSides  

 end subroutine createSides
!-------------------------------------------------------------------------
 subroutine writeComputationFile(isInitial,crp,finegrp,grp,brp,OUTFILE)
! write file for solver input
!use CoordinateRegister
!use BoundaryRegister

 IMPLICIT NONE

 logical :: isInitial 
 type(GridData) :: finegrp,grp
 type(CoordinateRegisterData) :: crp
 type(BoundaryRegisterData) :: brp
 integer :: OUTFILE 
 integer,pointer :: helpArray(:)
 
 integer :: i,j,ind1,allocateStatus,sideNumber,blockNumber,maxv,maxi,nums
 real :: x1(2),x2(2)
 real :: help
 real,pointer :: check(:,:)
 type(SideDataBlockInstance),pointer :: currentSideBlock

 if(isInitial) then
! write face data
  write(OUTFILE) brp%numberOfBoundaryNodes
  call writeBoundaryData(OUTFILE,brp,isInitial) 

! write sides
  write(OUTFILE) grp%numberOfSides,grp%sdp%dataBlockSize
  sideNumber = 0 
  blockNumber = 0
  currentSideBlock => grp%sdp%DataBlockChain%first
  do while(associated(currentSideBlock)) 
   do i=1,grp%sdp%dataBlockSize
    sideNumber = sideNumber + 1
    if(sideNumber==grp%numberOfSides) exit
   end do
  
   blockNumber = blockNumber + 1 
   if(blockNumber*grp%sdp%dataBlockSize>grp%numberOfSides) then 
    sideNumber = grp%numberOfSides - (blockNumber-1)*grp%sdp%dataBlockSize
   else
    sideNumber = grp%sdp%dataBlockSize
   end if 
   write(OUTFILE) ((currentSideBlock%block(i)%sideIndexes(j),i=1,sideNumber),j=1,2)
   write(OUTFILE) ((currentSideBlock%block(i)%sideCoefficients(j),i=1,sideNumber),j=1,3)
   write(OUTFILE) ((currentSideBlock%block(i)%sideLength),i=1,sideNumber)  
   currentSideBlock => currentSideBlock%next
  end do
 
 allocate(check(grp%numberOfNodes,3),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: getNodeVolume out of memory"

 check = 0.0

 write(*,*) "Checking..."
 nums = 0
 sideNumber = 0 
 currentSideBlock => grp%sdp%DataBlockChain%first
 do while(associated(currentSideBlock)) 
  do i=1,grp%sdp%dataBlockSize
   sideNumber = sideNumber + 1
   check(currentSideBlock%block(i)%sideIndexes(1),:) = check(currentSideBlock%block(i)%sideIndexes(1),:) +&
                    2.0*currentSideBlock%block(i)%sideCoefficients(:)
   check(currentSideBlock%block(i)%sideIndexes(2),:) = check(currentSideBlock%block(i)%sideIndexes(2),:) -&
                    2.0*currentSideBlock%block(i)%sideCoefficients(:)
   if(sideNumber==grp%numberOfSides) exit
  end do
  currentSideBlock => currentSideBlock%next
 end do

 do i=1,grp%brp%numberOfBoundarySides
  check(grp%brp%srp(i)%sideIndexes(1),:) = check(grp%brp%srp(i)%sideIndexes(1),:) + 4.0*grp%brp%srp(i)%sideCoefficients(:)
  check(grp%brp%srp(i)%sideIndexes(2),:) = check(grp%brp%srp(i)%sideIndexes(2),:) + 4.0*grp%brp%srp(i)%sideCoefficients(:)
 end do 

 maxv = 0.0
 do i=1,grp%numberOfNodes
  if(sum(check(i,:)*check(i,:))>0.0000001) then 
   write(*,*) "WARNING: Coefficients do not add up ",i,check(i,:)
   if(sum(check(i,:)*check(i,:))>maxv) then 
    maxv = sum(check(i,:)*check(i,:))
    maxi = i
   end if
   write(555,*) i,check(i,:)!,getCoor(i,crp)
  end if
 end do

  help = 0.0
  write(OUTFILE) grp%numberOfNodes   
  write(OUTFILE) (grp%controlVolumeNodeSize(i),i=1,grp%numberOfNodes) 
  do i=1,grp%numberOfNodes
   help = help + grp%controlVolumeNodeSize(i)
  end do
  write(*,*) "total volume: ",help

  call writeCoordinateData(OUTFILE,crp) 
 else

  write(OUTFILE) brp%numberOfBoundaryNodes
  call writeBoundaryData(OUTFILE,brp,isInitial) 

  write(OUTFILE) grp%numberOfSides,grp%sdp%dataBlockSize
  sideNumber = 0
  blockNumber = 0
  currentSideBlock => grp%sdp%DataBlockChain%first
  do while(associated(currentSideBlock)) 

   blockNumber = blockNumber + 1
   if(blockNumber*grp%sdp%dataBlockSize>grp%numberOfSides) then
    sideNumber = grp%numberOfSides - (blockNumber-1)*grp%sdp%dataBlockSize
   else
    sideNumber = grp%sdp%dataBlockSize
   end if
   write(OUTFILE) ((currentSideBlock%block(i)%sideIndexes(j),i=1,sideNumber),j=1,2)
   write(OUTFILE) ((currentSideBlock%block(i)%sideCoefficients(j),i=1,sideNumber),j=1,3)
   write(OUTFILE) ((currentSideBlock%block(i)%sideLength),i=1,sideNumber)

   currentSideBlock => currentSideBlock%next
  end do

  write(OUTFILE) grp%numberOfNodes   
  write(OUTFILE) (grp%controlVolumeNodeSize(i),i=1,grp%numberOfNodes) 

  help = 0.0
  do i=1,grp%numberOfNodes
   help = help + grp%controlVolumeNodeSize(i)
  end do
  write(*,*) "total volume: ",help

! write seed nodes

  write(OUTFILE) (grp%seedNodeRegister(i),i=1,grp%numberOfNodes)

! write mappings 
  call writeMappingToFile(OUTFILE,crp,finegrp,grp,brp,1) 

 end if

 end subroutine writeComputationFile
!-------------------------------------------------------------------------
 subroutine writeMappingToFile(OUTFILE,crp,finegrp,coarsegrp,brp,scheme)
! writes inter-grid transformations to computation file

 IMPLICIT NONE

 integer :: OUTFILE
 type(CoordinateRegisterData) :: crp
 type(GridData) :: finegrp,coarsegrp
 type(BoundaryRegisterData) :: brp
 integer :: scheme ! decides which mapping scheme to use

 integer :: i,j,allocateStatus,numberOfNodes,count,ind1,ind2
  
 write(OUTFILE) (coarsegrp%fineToCoarseIndexMapping(i),i=1,finegrp%numberOfNodes)

! fine to coarse

! use linear interpolation

 write(OUTFILE) (finegrp%controlVolumeNodeSize(i)/&
  coarsegrp%controlVolumeNodeSize(coarsegrp%fineToCoarseIndexMapping(i)),i=1,finegrp%numberOfNodes)

 end subroutine writeMappingToFile
!-------------------------------------------------------------------------
 logical function nodeIsInTetrahedra(xn,xp0,xp1,xp2,xp3)
 IMPLICIT NONE

 real :: xn(3),xp0(3),xp1(3),xp2(3),xp3(3)

 real :: xs(3),x1(3),x2(3),x3(3),buffx(3) 
 real :: r1(3),r2(3),r3(3),prod

 x1 = xp1
 x2 = xp2
 x3 = xp3
 xs = xn 

 ! make sure tetrahedra is numbered according to right hand rule

 r1 = x1-xp0
 r2 = x2-xp0
 r3 = x3-xp0
 prod = r3(1)*(r1(2)*r2(3)-r1(3)*r2(2))+r3(2)*(r1(3)*r2(1)-r1(1)*r2(3))+r3(3)*(r1(1)*r2(2)-r1(2)*r2(1))

 if(prod<0.0) then 
  buffx = x3
  x3 = x2
  x2 = buffx  
 else if(prod==0.0) then 
  nodeIsInTetrahedra = .false. 
  goto 311
 end if 

 r1 = x1-xp0
 r2 = x2-xp0
 r3 = xs-xp0
 prod = r3(1)*(r1(2)*r2(3)-r1(3)*r2(2))+r3(2)*(r1(3)*r2(1)-r1(1)*r2(3))+r3(3)*(r1(1)*r2(2)-r1(2)*r2(1))
 if(prod<0.0) then 
  nodeIsInTetrahedra = .false.
  goto 311
 end if
 
 r1 = x3-xp0
 r2 = x1-xp0
 r3 = xs-xp0
 prod = r3(1)*(r1(2)*r2(3)-r1(3)*r2(2))+r3(2)*(r1(3)*r2(1)-r1(1)*r2(3))+r3(3)*(r1(1)*r2(2)-r1(2)*r2(1))
 if(prod<0.0) then 
  nodeIsInTetrahedra = .false.
  goto 311
 end if
 
 r1 = x2-xp0
 r2 = x3-xp0
 r3 = xs-xp0
 prod = r3(1)*(r1(2)*r2(3)-r1(3)*r2(2))+r3(2)*(r1(3)*r2(1)-r1(1)*r2(3))+r3(3)*(r1(1)*r2(2)-r1(2)*r2(1))
 if(prod<0.0) then 
  nodeIsInTetrahedra = .false.
  goto 311
 end if

 r1 = x3-x1
 r2 = x2-x1 
 r3 = xs-x1
 prod = r3(1)*(r1(2)*r2(3)-r1(3)*r2(2))+r3(2)*(r1(3)*r2(1)-r1(1)*r2(3))+r3(3)*(r1(1)*r2(2)-r1(2)*r2(1))
 if(prod<0.0) then 
  nodeIsInTetrahedra = .false.
  goto 311
 end if

 nodeIsInTetrahedra = .true.

311 continue

 end function nodeIsInTetrahedra
!-------------------------------------------------------------------------
 logical function nodeIsInTriangle(xn,xp0,xp1,xp2)
 IMPLICIT NONE

 real :: xn(3),xp0(3),xp1(3),xp2(3)

 real :: xs(3),x1(3),x2(3),buffx(3)
 real :: r1(3),r2(3),r3(3),prod


 x1 = xp1
 x2 = xp2
 xs = xn

 ! project node to surface and check if it is within

 r1 = x1-xp0
 r2 = x2-xp0
 r3(1) = r1(2)*r2(3)-r1(3)*r2(2)
 r3(2) = r1(3)*r2(1)-r1(1)*r2(3)
 r3(3) = r1(1)*r2(2)-r1(2)*r2(1) 

 prod = sqrt(sum(r3*r3))
 if(prod>1.0e-15) then 
  r3 = r3/prod
 else
  nodeIsInTriangle = .false.
  goto 312
 end if

 xs = xs-xp0
 xs = xs-sum(xs*r3)*r3  

! write(*,*) "S1: ",xs

 prod = xs(1)*(r3(2)*r1(3)-r3(3)*r1(2))+xs(2)*(r3(3)*r1(1)-r3(1)*r1(3))+xs(3)*(r3(1)*r1(2)-r3(2)*r1(1)) 

 if(prod<0.0) then 
  nodeIsInTriangle = .false.
  goto 312
 end if


 prod = xs(1)*(r2(2)*r3(3)-r2(3)*r3(2))+xs(2)*(r2(3)*r3(1)-r2(1)*r3(3))+xs(3)*(r2(1)*r3(2)-r2(2)*r3(1))

 if(prod<0.0) then 
  nodeIsInTriangle = .false.
  goto 312
 end if

 r1 = xp2 -xp1
 xs = xn-xp1
 xs = xs-sum(xs*r3)*r3

 prod = xs(1)*(r3(2)*r1(3)-r3(3)*r1(2))+xs(2)*(r3(3)*r1(1)-r3(1)*r1(3))+xs(3)*(r3(1)*r1(2)-r3(2)*r1(1))

 if(prod<0.0) then 
  nodeIsInTriangle = .false.
  goto 312
 end if

 nodeIsInTriangle = .true.

312 continue 

 end function nodeIsInTriangle
!-------------------------------------------------------------------------
 logical function makeVolumeMappingCoefficients(xp,x0v,x1v,x2v,x3v,coeff)
 IMPLICIT NONE

  real :: x0v(3),x1v(3),x2v(3),x3v(3),xp(3),coeff(3)

  real :: x0,x1,x2,x3,y0,y1,y2,y3,z0,z1,z2,z3
  real :: D,A00,A10,A20,A30,A01,A11,A21,A31,A02,A12,A22,A32
  real :: A03,A13,A23,A33
  integer :: i,j,count


   makeVolumeMappingCoefficients = .true.

  ! full volume mapping

   x0 = x0v(1)
   y0 = x0v(2)
   z0 = x0v(3)
   x1 = x1v(1)
   y1 = x1v(2)
   z1 = x1v(3)
   x2 = x2v(1)
   y2 = x2v(2)
   z2 = x2v(3)
   x3 = x3v(1)
   y3 = x3v(2)
   z3 = x3v(3)
  
   A00 = y1*(z2-z3)-z1*(y2-y3)+(y2*z3-z2*y3)
   A10 = -(y0*(z2-z3)-z0*(y2-y3)+(y2*z3-z2*y3))
   A20 = y0*(z1-z3)-z0*(y1-y3)+(y1*z3-y3*z1)
   A30 = -(y0*(z1-z2)-z0*(y1-y2)+(y1*z2-y2*z1))

   A01 = -(x1*(z2-z3)-z1*(x2-x3)+(x2*z3-x3*z2))
   A11 = x0*(z2-z3)-z0*(x2-x3)+(x2*z3-x3*z2)
   A21 = -(x0*(z1-z3)-z0*(x1-x3)+(x1*z3-x3*z1))
   A31 = x0*(z1-z2)-z0*(x1-x2)+(x1*z2-x2*z1)

   A02 = x1*(y2-y3)-y1*(x2-x3)+(x2*y3-x3*y2)
   A12 = -(x0*(y2-y3)-y0*(x2-x3)+(x2*y3-x3*y2))
   A22 = x0*(y1-y3)-y0*(x1-x3)+(x1*y3-x3*y1)
   A32 = -(x0*(y1-y2)-y0*(x1-x2)+(x1*y2-y1*x2))

   A03 = -(x1*(y2*z3-y3*z2)-y1*(x2*z3-x3*z2)+z1*(x2*y3-x3*y2))
   A13 = x0*(y2*z3-y3*z2)-y0*(x2*z3-x3*z2)+z0*(x2*y3-x3*y2)
   A23 = -(x0*(y1*z3-y3*z1)-y0*(x1*z3-x3*z1)+z0*(x1*y3-x3*y1))
   A33 = x0*(y1*z2-y2*z1)-y0*(x1*z2-x2*z1)+z0*(x1*y2-x2*y1)

   D = x0*A00+y0*A01+z0*A02+A03

   if(abs(D)>0.0) then 
    coeff(1) = (A10*xp(1) + A11*xp(2) + A12*xp(3)+A13)/D
    coeff(2) = (A20*xp(1) + A21*xp(2) + A22*xp(3)+A23)/D
    coeff(3) = (A30*xp(1) + A31*xp(2) + A32*xp(3)+A33)/D
   else 
    makeVolumeMappingCoefficients = .false.
   end if

   ! check whether coefficients are positive
 
  if(coeff(1)<-1.0e-12) then
   makeVolumeMappingCoefficients = .false.
  end if
  if(coeff(2)<-1.0e-12) then
   makeVolumeMappingCoefficients = .false.
  end if
  if(coeff(3)<-1.0e-12) then
   makeVolumeMappingCoefficients = .false.
  end if
  if(1.0-sum(coeff)<-1.0e-12) then
   makeVolumeMappingCoefficients = .false.
  end if
 
    
 end function makeVolumeMappingCoefficients
!-------------------------------------------------------------------------
 logical function makeSurfaceMappingCoefficients(xp,x0,x1,x2,coeff)
 IMPLICIT NONE

  real :: x0(3),x1(3),x2(3),xp(3),coeff(3)

  real :: r1(3),r2(3),normal(3),xs(3),s,t,det,det1,det2,det3
  integer :: usedet,usel

  makeSurfaceMappingCoefficients = .true.

  coeff(3) = 0.0 

  r1 = x1-x0
  r2 = x2-x0 
  normal(1) = r1(2)*r2(3)-r1(3)*r2(2)
  normal(2) = r1(3)*r2(1)-r1(1)*r2(3)
  normal(3) = r1(1)*r2(2)-r1(2)*r2(1)

  det = sqrt(sum(normal*normal))
  if(det>1.0e-15) then 
   normal = normal/det
  else
   makeSurfaceMappingCoefficients = .false.
   goto 575
  end if

 ! map point onto surface

  xs = xp-x0
  xs = xs-sum(xs*normal)*normal

  det1 = r1(1)*r2(2)-r1(2)*r2(1)
  det2 = r1(1)*r2(3)-r1(3)*r2(1)
  det3 = r1(2)*r2(3)-r1(3)*r2(2)

  if(abs(det1)>abs(det2).and.abs(det1)>abs(det3)) then 
   usedet = 1
   det = det1 
  else if(abs(det2)>abs(det1).and.abs(det2)>abs(det3)) then
   usedet = 2
   det = det2 
  else
   usedet = 3
   det = det3
  end if

  if(abs(det)>1.0e-15) then 
   if(usedet==1) then 
    t = (r1(1)*xs(2)-r1(2)*xs(1))/det
   else if(usedet==2) then 
    t = (r1(1)*xs(3)-r1(3)*xs(1))/det
   else
    t = (r1(2)*xs(3)-r1(3)*xs(2))/det 
   end if 
   if(abs(r1(1))>abs(r1(2)).and.abs(r1(1))>abs(r1(3))) then 
    usel = 1
   else if(abs(r1(2))>abs(r1(1)).and.abs(r1(2))>abs(r1(3))) then 
    usel = 2 
   else
    usel = 3 
   end if
   if(abs(r1(usel))>1.0e-15) then 
    s = (xs(usel)-r2(usel)*t)/r1(usel)
   else
    makeSurfaceMappingCoefficients = .false.
    goto 575
   end if
   coeff(1) = s
   coeff(2) = t 
  else
   makeSurfaceMappingCoefficients = .false.
   goto 575
  end if                                   
  ! check whether coefficients are positive

  if(coeff(1)<-1.0e-12) then 
   makeSurfaceMappingCoefficients = .false.
   goto 575
  end if   
  if(coeff(2)<-1.0e-12) then 
   makeSurfaceMappingCoefficients = .false.
   goto 575
  end if   
  if(1.0-sum(coeff)<-1.0e-12) then 
   makeSurfaceMappingCoefficients = .false.
  end if   


 575 continue
!  write(*,*) "P: ",coeff,makeSurfaceMappingCoefficients
  
 end function makeSurfaceMappingCoefficients
!-------------------------------------------------------------------------
 subroutine writeParallelComputationFiles(isInitial,numberOfGrids,currentGrid,problemName,finegrp,grp,crp,crpPrev,crpPrev2,inputNumber,numberOfStepsPerCycle)
 IMPLICIT NONE

 logical :: isInitial 
 type(CoordinateRegisterData) :: crp,crpPrev,crpPrev2
 integer :: numberOfGrids,currentGrid
 character*80 :: problemName
 type(GridData) :: grp
 type(GridData) :: finegrp

 integer,pointer :: oldToNewNodeMapping(:)
 integer,pointer :: ghostNodeFlagArray(:)
 integer,pointer :: numberOfBoundarySendNodes(:)
 integer,pointer :: numberOfBoundaryReceiveNodes(:)
 integer,pointer :: numberOfSendNodes(:)
 integer,pointer :: numberOfReceiveNodes(:)
 integer,pointer :: numberOfDomainSides(:)
 integer,pointer :: numberOfDomainBoundarySides(:)
 integer,pointer :: nodeStartIndex(:)
 integer,pointer :: nodeRegister(:)
 real,pointer :: nodeVecArray(:,:)

 integer,pointer,save :: grid1FullOld2NewMapping(:)
 integer :: inputNumber,itLen,numberOfStepsPerCycle
 character :: itExt*5,itForm*4

 character*3 :: char_ipr
 character*4 :: fm
 integer :: OUTFILE,REGFILE

 integer :: i,j,k,ind1,ind2,ind,indicator,domain1,domain2,sideNumber,allocateStatus,buffi,numberOfComSides
 integer :: localBoundaryNodesRegister(-20:0),localInd1,localInd2,numberOfLocalNodes
 integer :: numberOfTEs,numberOfEIS,currentIndicator,numberOfNodesInDomain,count,numberOfInternalSides
 real :: weights(3),sideLen,GCLWeight,GCLWeightP
 integer :: nameLength,extl,nodesInDomainCount,boundaryNodesInDomainCount,nn,numberOfBoundaryComSides
 integer :: boundaryNormalsInDomainCount,totalNumberOfBoundaryNodes,numberOfBoundaryGhostNodes
 integer :: boundaryNodeCount,maxNumberOfComNodes

 integer :: numberOfLocalTripNodes,numberOfLocalTripFieldNodes

 integer,pointer :: localTripMapping(:)

 logical,pointer :: sideIsInternal(:)

 type(SideDataBlockInstance),pointer :: currentSideBlock

 integer :: numberOfLocalComNodes,numberOfComNodes,numberOfFineNonComNodes
 integer,pointer :: localComMapping(:)
 logical,pointer :: isComNode(:)

 integer :: l,m

! Start of OH Modification
 integer :: numberOfInterfaceSides
 integer,pointer :: CommSidesArray(:,:)
! End of OH Modification

 write(*,*) "writing parallel files..."


 if(currentGrid==1) then
  call meshPartition2(grp%numberOfGridDomains,grp%sdp)
 end if

! flag whether sides are communication sides or not

 allocate(sideIsInternal(grp%sdp%numberOfSides),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeParallelComputationFiles out of memory"
 allocate(CommSidesArray(2,grp%sdp%numberOfSides),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeParallelComputationFiles out of memory"

! Start of OH Modification
  numberOfInterfaceSides = 0
! End of OH Modification
  sideIsInternal = .true.
  sideNumber = 0
  currentSideBlock => grp%sdp%DataBlockChain%first
  do while(associated(currentSideBlock))
   do j=1,grp%sdp%dataBlockSize
    sideNumber = sideNumber+1
    ind1 = currentSideBlock%block(j)%sideIndexes(1)
    ind2 = currentSideBlock%block(j)%sideIndexes(2)
    domain1 = grp%sdp%parallelizationColorArray(ind1)
    domain2 = grp%sdp%parallelizationColorArray(ind2)
    if(domain1.ne.domain2) then
! Start of OH Modification
     numberOfInterfaceSides = numberOfInterfaceSides + 1
     CommSidesArray(1,numberOfInterfaceSides) = ind1
     CommSidesArray(2,numberOfInterfaceSides) = ind2
! End of OH Modification
     sideIsInternal(sideNumber) = .false.
    end if
    if(sideNumber==grp%numberOfSides) exit
   end do
   currentSideBlock => currentSideBlock%next
  end do

  write(*,*) 'number of sides, number of Comm sides'
  write(*,*) sideNumber,numberOfInterfaceSides


 nameLength = len_trim(problemName)

 OUTFILE = 73
 REGFILE = 77

 allocate(grp%fullOldToNewNodeMapping(grp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeParallelComputationFiles out of memory"
 allocate(oldToNewNodeMapping(grp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeParallelComputationFiles out of memory"
 allocate(ghostNodeFlagArray(grp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeParallelComputationFiles out of memory"
 allocate(numberOfBoundarySendNodes(grp%numberOfGridDomains),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeParallelComputationFiles out of memory"
 allocate(numberOfBoundaryReceiveNodes(grp%numberOfGridDomains),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeParallelComputationFiles out of memory"
 allocate(numberOfSendNodes(grp%numberOfGridDomains),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeParallelComputationFiles out of memory"
 allocate(numberOfReceiveNodes(grp%numberOfGridDomains),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeParallelComputationFiles out of memory"
 allocate(numberOfDomainSides(grp%numberOfGridDomains),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeParallelComputationFiles out of memory"
 allocate(numberOfDomainBoundarySides(grp%numberOfGridDomains),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeParallelComputationFiles out of memory"
 allocate(nodeStartIndex(grp%numberOfNodes+1),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeMappingToFile out of memory"
 allocate(nodeRegister(finegrp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeMappingToFile out of memory"
 allocate(grp%numberOfIntNodesInPart(grp%numberOfGridDomains),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: Not enough memory in meshPartition"
 if(currentGrid==1) then
  allocate(localTripMapping(grp%numberOfTripNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: writeMappingToFile out of memory"
  nullify(localComMapping)
  nullify(isComNode)
  allocate(grid1FullOld2NewMapping(grp%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: writeParallelComputationFiles out of memory"
 else
  nullify(localTripMapping)
  allocate(localComMapping(grp%numberOfFineComNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: writeMappingToFile out of memory"
  allocate(isComNode(finegrp%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: writeMappingToFile out of memory"
 end if


 ! set up is com node register
 if(currentGrid>1) then
  isComNode = .false.
  do k=1,grp%numberOfFineComNodes
   isComNode(grp%fineComNodes(k)) = .true.
  end do
 end if

  grp%fullOldToNewNodeMapping = 0

  sideNumber = 0
  numberOfDomainSides = 0
  currentSideBlock => grp%sdp%DataBlockChain%first
  do while(associated(currentSideBlock))
   do i=1,grp%sdp%dataBlockSize
    sideNumber = sideNumber + 1
    ind1 = currentSideBlock%block(i)%sideIndexes(1)
    ind2 = currentSideBlock%block(i)%sideIndexes(2)
   
    domain1 = grp%sdp%parallelizationColorArray(ind1)
    domain2 = grp%sdp%parallelizationColorArray(ind2)

    if(domain1.ne.domain2) then
     ! count number of edges originating from domain
     if(domain1<domain2) then
      numberOfDomainSides(domain1) = numberOfDomainSides(domain1) + 1
     else
      numberOfDomainSides(domain2) = numberOfDomainSides(domain2) + 1
     end if
    else
     numberOfDomainSides(domain1) = numberOfDomainSides(domain1) + 1
    end if
    if(sideNumber==grp%numberOfSides) exit
   end do
   currentSideBlock => currentSideBlock%next
  end do

  ! count number of boundary sides in domains

  numberOfDomainBoundarySides = 0
  do i=1,grp%brp%numberOfBoundarySides
   ind1 = grp%brp%srp(i)%sideIndexes(1)
   ind2 = grp%brp%srp(i)%sideIndexes(2)
   domain1 = grp%sdp%parallelizationColorArray(ind1)
   domain2 = grp%sdp%parallelizationColorArray(ind2)
   if(domain1.le.domain2) then
    numberOfDomainBoundarySides(domain1) = numberOfDomainBoundarySides(domain1) + 1
   else
    numberOfDomainBoundarySides(domain2) = numberOfDomainBoundarySides(domain2) + 1
   end if
  end do

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
  write(itExt,itForm) inputNumber    ! time step number

  if(currentGrid==1) then
   open(REGFILE,file='plotreg.reg',form='unformatted',status='unknown')
   write(REGFILE) grp%numberOfNodes,grp%numberOfGridDomains
  end if

  do i=1,grp%numberOfGridDomains
   ! make renumber register for each grid domain between 1-numberOfNodes

   nodesInDomainCount = 0
   boundaryNodesInDomainCount = 0
   boundaryNormalsInDomainCount = 0
   oldToNewNodeMapping = 0
   do j=1,grp%numberOfNodes
    if(grp%sdp%parallelizationColorArray(j)==i) then
     nodesInDomainCount = nodesInDomainCount + 1
     if(j.le.grp%brp%numberOfBoundaryNodes) then
      boundaryNodesInDomainCount = &
         boundaryNodesInDomainCount + 1
      if(grp%brp%boundaryNodes(j)%indicator>0) then
       boundaryNormalsInDomainCount = &
        boundaryNormalsInDomainCount + 1
      else
       write(*,*) "TT: ",j,grp%brp%boundaryNodes(j)%indicator
      end if
     end if
     oldToNewNodeMapping(j) = nodesInDomainCount
    end if
   end do
   grp%numberOfIntNodesInPart(i) = nodesInDomainCount
   write(*,*) "nodes in domain ",i,": ",nodesInDomainCount

  ! write to file


   if (i.le.9) then
    fm='(i1)'
    extl = 1
   else if (i.le.99) then
    fm='(i2)'
    extl = 2
   else if (i.le.999) then
    fm='(i3)'
    extl = 3
   endif

   write(char_ipr,fm) i           ! domain number

   if(currentGrid==1) then
    if(inputNumber.eq.numberOfStepsPerCycle.and.numberOfStepsPerCycle.eq.1) then
     open(OUTFILE,file=problemName(1:nameLength)//'.sol'//'_'//char_ipr(1:extl),&
               form='unformatted',status='unknown')
    else
     open(OUTFILE,file=problemName(1:nameLength)//'_'//itExt(1:itLen)//'.sol'//'_'//char_ipr(1:extl),&
               form='unformatted',status='unknown')
    end if
    write(OUTFILE) numberOfGrids
    write(OUTFILE) i,grp%numberOfGridDomains
   else
    if(inputNumber.eq.numberOfStepsPerCycle.and.numberOfStepsPerCycle.eq.1) then
     open(OUTFILE,file=problemName(1:nameLength)//'.sol'//'_'//char_ipr(1:extl),&
               form='unformatted',status='old',&
               position='append')
    else
     open(OUTFILE,file=problemName(1:nameLength)//'_'//itExt(1:itLen)//'.sol'//'_'//char_ipr(1:extl),&
               form='unformatted',status='old',&
               position='append')
    end if
   end if

   ! write boundary sides

   ! first count number of boundary ghost nodes

   sideNumber = 0
   ghostNodeFlagArray = 0
   numberOfBoundaryGhostNodes = 0
   currentSideBlock => grp%sdp%DataBlockChain%first
   do while(associated(currentSideBlock))
    do j=1,grp%sdp%dataBlockSize
     sideNumber = sideNumber+1
     ind1 = currentSideBlock%block(j)%sideIndexes(1)
     ind2 = currentSideBlock%block(j)%sideIndexes(2)
     domain1 = grp%sdp%parallelizationColorArray(ind1)
     domain2 = grp%sdp%parallelizationColorArray(ind2)
     if(domain1<domain2) then
      if(domain1==i) then
       if(ind2.le.grp%brp%numberOfBoundaryNodes) then
        if(ghostNodeFlagArray(ind2)==0) then
         ghostNodeFlagArray(ind2) = 1
         numberOfBoundaryGhostNodes = numberOfBoundaryGhostNodes + 1
        end if
       end if
      end if
     else if(domain1>domain2) then
      if(domain2==i) then
       if(ind1.le.grp%brp%numberOfBoundaryNodes) then
        if(ghostNodeFlagArray(ind1)==0) then
         ghostNodeFlagArray(ind1) = 1
         numberOfBoundaryGhostNodes = numberOfBoundaryGhostNodes + 1
        end if
       end if
      end if
     end if
     if(sideNumber==grp%numberOfSides) exit
    end do
    currentSideBlock => currentSideBlock%next
   end do

   totalNumberofBoundaryNodes = boundaryNodesInDomainCount + numberOfBoundaryGhostNodes


   ! shift internal nodes

   do j=1,grp%numberOfNodes
    if(oldToNewNodeMapping(j)>boundaryNodesInDomainCount) then
     oldToNewNodeMapping(j) = oldToNewNodeMapping(j) + numberOfBoundaryGhostNodes
    end if
   end do

   ! now put ghost nodes

   count = 0
   sideNumber = 0
   numberOfLocalNodes = numberOfBoundaryGhostNodes+nodesInDomainCount
   currentSideBlock => grp%sdp%DataBlockChain%first
   do while(associated(currentSideBlock))
    do j=1,grp%sdp%dataBlockSize
     sideNumber = sideNumber+1
     ind1 = currentSideBlock%block(j)%sideIndexes(1)
     ind2 = currentSideBlock%block(j)%sideIndexes(2)
     domain1 = grp%sdp%parallelizationColorArray(ind1)
     domain2 = grp%sdp%parallelizationColorArray(ind2)
     if(domain1<domain2) then
      if(domain1==i) then
       if(ind2.le.grp%brp%numberOfBoundaryNodes) then
        if(oldToNewNodeMapping(ind2)==0) then
         if(grp%brp%boundaryNodes(ind2)%indicator>0) then
          boundaryNormalsInDomainCount = &
          boundaryNormalsInDomainCount + 1
         end if
         count = count + 1
         oldToNewNodeMapping(ind2) = boundaryNodesInDomainCount+count
        end if
       else
        if(oldToNewNodeMapping(ind2)==0) then
         numberOfLocalNodes = numberOfLocalNodes + 1
         oldToNewNodeMapping(ind2) = numberOfLocalNodes
        end if
       end if
      end if
     else if(domain1>domain2) then
      if(domain2==i) then
       if(ind1.le.grp%brp%numberOfBoundaryNodes) then
        if(oldToNewNodeMapping(ind1)==0) then
         if(grp%brp%boundaryNodes(ind1)%indicator>0) then
          boundaryNormalsInDomainCount = &
          boundaryNormalsInDomainCount + 1
         end if
         count = count + 1
         oldToNewNodeMapping(ind1) = boundaryNodesInDomainCount+count
        end if
       else
        if(oldToNewNodeMapping(ind1)==0) then
         numberOfLocalNodes = numberOfLocalNodes + 1
         oldToNewNodeMapping(ind1) = numberOfLocalNodes
        end if
       end if
      end if
     end if
     if(sideNumber==grp%numberOfSides) exit
    end do
    currentSideBlock => currentSideBlock%next
   end do


   ! remap nodes for cache efficiency

   ! write communication boundary sides

   write(OUTFILE) numberOfDomainBoundarySides(i)
   write(*,*) "number of domain boundary sides: ",numberOfDomainBoundarySides(i)
   ghostNodeFlagArray = 0
   numberOfBoundarySendNodes = 0
   numberOfBoundaryReceiveNodes = 0
   numberOfBoundaryComSides = 0
   boundaryNodeCount = 0
   do j=1,grp%brp%numberOfBoundarySides
    ind1 = grp%brp%srp(j)%sideIndexes(1)
    ind2 = grp%brp%srp(j)%sideIndexes(2)
    domain1 = grp%sdp%parallelizationColorArray(ind1)
    domain2 = grp%sdp%parallelizationColorArray(ind2)
    localInd1 = oldToNewNodeMapping(ind1)
    localInd2 = oldToNewNodeMapping(ind2)
    if(domain1<domain2) then
     if(domain1==i) then
      if(ghostNodeFlagArray(ind2)==0) then
       boundaryNodeCount = boundaryNodeCount + 1
       ghostNodeFlagArray(ind2) = 1
       numberOfBoundarySendNodes(domain2) = numberOfBoundarySendNodes(domain2) + 1
      end if
      write(OUTFILE) localInd1,localInd2,grp%brp%srp(j)%indicator,grp%brp%srp(j)%sideCoefficients,&
                                 grp%brp%srp(j)%GCLCoefficients,grp%brp%srp(j)%GCLCoefficientsP
      numberOfBoundaryComSides = numberOfBoundaryComSides + 1
     else if(domain2==i) then
      if(ghostNodeFlagArray(ind2)==0) then
       ghostNodeFlagArray(ind2) = 1
       boundaryNodeCount = boundaryNodeCount + 1
       numberOfBoundaryReceiveNodes(domain1) = numberOfBoundaryReceiveNodes(domain1) + 1
      end if
     end if
    else if(domain1>domain2) then
     if(domain2==i) then
      ! ind1 is ghost
      if(ghostNodeFlagArray(ind1)==0) then
       boundaryNodeCount = boundaryNodeCount + 1
       ghostNodeFlagArray(ind1) = 1
       numberOfBoundarySendNodes(domain1) = numberOfBoundarySendNodes(domain1) + 1
      end if
      write(OUTFILE) localInd1,localInd2,grp%brp%srp(j)%indicator,grp%brp%srp(j)%sideCoefficients,&
                                 grp%brp%srp(j)%GCLCoefficients,grp%brp%srp(j)%GCLCoefficientsP
      numberOfBoundaryComSides = numberOfBoundaryComSides + 1
     else if(domain1==i) then
      if(ghostNodeFlagArray(ind1)==0) then
       ghostNodeFlagArray(ind1) = 1
       boundaryNodeCount = boundaryNodeCount + 1
       numberOfBoundaryReceiveNodes(domain2) = numberOfBoundaryReceiveNodes(domain2) + 1
      end if
     end if
    end if
   end do

  ! write internal boundary sides

   count = 0
   do j=1,grp%brp%numberOfBoundarySides
    ind1 = grp%brp%srp(j)%sideIndexes(1)
    ind2 = grp%brp%srp(j)%sideIndexes(2)

    domain1 = grp%sdp%parallelizationColorArray(ind1)
    domain2 = grp%sdp%parallelizationColorArray(ind2)
    if(domain1==i.and.domain2==i) then
     localInd1 = oldToNewNodeMapping(ind1)
     localInd2 = oldToNewNodeMapping(ind2)
     weights = grp%brp%srp(j)%sideCoefficients
     write(OUTFILE) localInd1,localInd2,grp%brp%srp(j)%indicator,weights,&
          grp%brp%srp(j)%GCLCoefficients,grp%brp%srp(j)%GCLCoefficientsP
     count = count + 1
    end if
   end do

  if(count+numberOfBoundaryComSides.ne.numberOfDomainBoundarySides(i)) then
    write(*,*) count,numberOfBoundaryComSides,numberOfDomainBoundarySides(i)
    STOP "ERROR: Something's gone wrong in writePar - 12"
  end if


   ! make local boundary indicator register

   localBoundaryNodesRegister = 0
   do j=1,grp%brp%numberOfBoundaryNodes
    if(oldToNewNodeMapping(j)>0) then
     indicator = grp%brp%boundaryNodes(j)%indicator
     if(indicator>0) then
      localBoundaryNodesRegister(-21+indicator*2) = localBoundaryNodesRegister(-21+indicator*2) + 1
     end if
    end if
   end do

   localBoundaryNodesRegister(-20) = 1
   do j=-18,-2,2
    localBoundaryNodesRegister(j+1) = localBoundaryNodesRegister(j+1) + localBoundaryNodesRegister(j-1)
    localBoundaryNodesRegister(j) = localBoundaryNodesRegister(j-1) + 1
   end do

   write(OUTFILE) totalNumberOfBoundaryNodes,boundaryNormalsInDomainCount
   write(*,*) "number of boundary nodes, normals: ",totalNumberOfBoundaryNodes,boundaryNormalsInDomainCount

   do j=-20,0
    write(OUTFILE)  localBoundaryNodesRegister(j)
   end do

!   write(OUTFILE) boundaryNormalsInDomainCount
  
   count = 0
   do currentIndicator=1,9
    do j=1,grp%brp%numberOfBoundaryNodes
     if(oldToNewNodeMapping(j)>0) then
      if(grp%brp%boundaryNodes(j)%indicator==currentIndicator) then
       write(OUTFILE) oldToNewNodeMapping(j) ,grp%brp%boundaryNodes(j)%normal
       count = count + 1
       if(oldToNewNodeMapping(j)>totalNumberOfBoundaryNodes) then 
        write(*,*) "rea: ",count,oldToNewNodeMapping(j),j
       end if
      else if(grp%brp%boundaryNodes(j)%indicator>9) then
       write(*,*) grp%brp%boundaryNodes(j)%indicator
       STOP "ERROR: Boundary indicator out of range in writePar"
      end if
     end if
    end do
   end do

   if(count.ne.boundaryNormalsInDomainCount) then
    write(*,*) count,boundaryNormalsInDomainCount
    STOP "ERROR: Something's wrong in writePar - 10"
   end if

   ! first count

   numberOfEIS = 0
   do j=1,grp%brp%numberOfEngineInletSides
    ind1 = grp%brp%engineInletSideIndexes(j,1)
    ind2 = grp%brp%engineInletSideIndexes(j,2)
    domain1 = grp%sdp%parallelizationColorArray(ind1)
    domain2 = grp%sdp%parallelizationColorArray(ind2)
!   if(domain1==i.and.domain2==i)  numberOfEIS = numberOfEIS + 1
    if((domain1.le.domain2.and.domain1.eq.i).or.(domain2.le.domain1.and.domain2.eq.i))   numberOfEIS = numberOfEIS + 1
!   if(oldToNewNodeMapping(ind1)>0) numberOfEIS = numberOfEIS + 1
!   if(oldToNewNodeMapping(ind1)>0.and.oldToNewNodeMapping(ind2)>0) numberOfEIS = numberOfEIS + 1
   end do

   write(OUTFILE) numberOfEIS
   write(*,*) "number of engine inlet sides: ",numberOfEIS
   do j=1,grp%brp%numberOfEngineInletSides
    ind1 = grp%brp%engineInletSideIndexes(j,1)
    ind2 = grp%brp%engineInletSideIndexes(j,2)
    domain1 = grp%sdp%parallelizationColorArray(ind1)
    domain2 = grp%sdp%parallelizationColorArray(ind2)
!   if(domain1==i.and.domain2==i)  then
    if((domain1.le.domain2.and.domain1.eq.i).or.(domain2.le.domain1.and.domain2.eq.i)) then
!   if(oldToNewNodeMapping(ind1)>0) then
!   if(oldToNewNodeMapping(ind1)>0.and.oldToNewNodeMapping(ind2)>0) then
     localInd1 = oldToNewNodeMapping(ind1)
     localInd2 = oldToNewNodeMapping(ind2)
 if(localInd1.eq.0.or.localInd2.eq.0) then
  write(*,'(5I8)') j,localInd1,localInd2, ind1, ind2
 end if
     write(OUTFILE) localInd1,localInd2,grp%brp%engineInletSideIndexes(j,3),grp%brp%engineInletSideCoefficients(j,1:3)
    end if
   end do

 ! write sides
   write(OUTFILE) numberOfDomainSides(i)

   ! first write communication sides

   sideNumber = 0
   numberOfComSides = 0
   currentSideBlock => grp%sdp%DataBlockChain%first
   count = 0
   do while(associated(currentSideBlock))
    do j=1,grp%sdp%dataBlockSize
     sideNumber = sideNumber+1
     ind1 = currentSideBlock%block(j)%sideIndexes(1)
     ind2 = currentSideBlock%block(j)%sideIndexes(2)
     domain1 = grp%sdp%parallelizationColorArray(ind1)
     domain2 = grp%sdp%parallelizationColorArray(ind2)

     if(domain1<domain2) then
      ! nodes are in different domains
      if(domain1==i) then
       ! ind2 is ghost
       localInd1 = oldToNewNodeMapping(ind1)
       localInd2 = oldToNewNodeMapping(ind2)
       if(ghostNodeFlagArray(ind2)==0) then
        ghostNodeFlagArray(ind2) = 1
        if(ind2.le.grp%brp%numberOfBoundaryNodes) then
         boundaryNodeCount = boundaryNodeCount + 1
        end if
       end if
       weights = currentSideBlock%block(j)%sideCoefficients
       GCLweight = currentSideBlock%block(j)%GCLCoefficient
       GCLweightP = currentSideBlock%block(j)%GCLCoefficientP
       sideLen = currentSideBlock%block(j)%sideLength
       numberOfComSides = numberOfComSides + 1
       count = count + 1
       write(OUTFILE) localInd1,localInd2,weights,GCLweight,GCLweightP,sideLen
      else if(domain2==i) then
       if(ghostNodeFlagArray(ind2)==0) then
        ghostNodeFlagArray(ind2) = 1
       end if
      end if
     else if(domain1>domain2) then
      if(domain2==i) then
       ! ind1 is ghost
       localInd1 = oldToNewNodeMapping(ind1)
       localInd2 = oldToNewNodeMapping(ind2)
       if(ghostNodeFlagArray(ind1)==0) then
        ghostNodeFlagArray(ind1) = 1
        if(ind1.le.grp%brp%numberOfBoundaryNodes) then
         boundaryNodeCount = boundaryNodeCount + 1
        end if
       end if
       weights = currentSideBlock%block(j)%sideCoefficients
       GCLweight = currentSideBlock%block(j)%GCLCoefficient
       GCLweightP = currentSideBlock%block(j)%GCLCoefficientP
       sideLen = currentSideBlock%block(j)%sideLength
       numberOfComSides = numberOfComSides + 1
       count = count + 1
       write(OUTFILE) localInd1,localInd2,weights,GCLweight,GCLweightP,sideLen
      else if(domain1==i) then
       if(ghostNodeFlagArray(ind1)==0) then
        ghostNodeFlagArray(ind1) = 1
       end if
      end if
     end if
     if(sideNumber==grp%numberOfSides) exit
    end do
    currentSideBlock => currentSideBlock%next
   end do

   numberOfNodesInDomain = numberOfLocalNodes ! this number includes ghost nodes

   ! now write internal sides

   sideNumber = 0
   numberOfInternalSides = 0
   currentSideBlock => grp%sdp%DataBlockChain%first
   do while(associated(currentSideBlock))
    do j=1,grp%sdp%dataBlockSize
     sideNumber = sideNumber+1
     ind1 = currentSideBlock%block(j)%sideIndexes(1)
     ind2 = currentSideBlock%block(j)%sideIndexes(2)
     domain1 = grp%sdp%parallelizationColorArray(ind1)
     domain2 = grp%sdp%parallelizationColorArray(ind2)
     if(domain1==domain2) then
      if(grp%sdp%parallelizationColorArray(ind1)==i) then
       localInd1 = oldToNewNodeMapping(ind1)
       localInd2 = oldToNewNodeMapping(ind2)
       weights = currentSideBlock%block(j)%sideCoefficients
       sideLen = currentSideBlock%block(j)%sideLength
       numberOfInternalSides = numberOfInternalSides + 1
       count = count + 1
       write(OUTFILE) localInd1,localInd2,weights,currentSideBlock%block(j)%GCLCoefficient,&
                  currentSideBlock%block(j)%GCLCoefficientP,sideLen
      end if
     end if
     if(sideNumber==grp%numberOfSides) exit
    end do
    currentSideBlock => currentSideBlock%next
   end do

   if(count.ne.numberOfDomainSides(i)) then
    write(*,*) count,numberOfDomainSides(i)
    STOP "ERROR: Something's gone wrong in writePar - 21"
   end if



   ! write node volume and distance


   write(OUTFILE) numberOfLocalNodes
   if(currentGrid==1) then 
    do j=1,grp%numberOfNodes
     if(oldToNewNodeMapping(j).ne.0) then
      if(j==17667) write(*,*) "MM: ",getCoor(j,crp),getCoor(j,crpPrev),getCoor(j,crpPrev2)
      write(OUTFILE) oldToNewNodeMapping(j),grp%controlVolumeNodeSize(j),&
        grp%prevControlVolume(j),grp%prevControlVolume2(j),&
        grp%wallDistance(j),oldToNewNodemapping(grp%wallDistanceBoundaryNodeArray(j)),&
        getCoor(j,crp),getCoor(j,crpPrev),getCoor(j,crpPrev2)
     end if
    end do
   else
    do j=1,grp%numberOfNodes
     if(oldToNewNodeMapping(j).ne.0) then
!      write(*,*) "TTX: ",grp%seedNodeRegister(j)
      write(OUTFILE) oldToNewNodeMapping(j),grp%controlVolumeNodeSize(j),&
        grp%wallDistance(j),oldToNewNodeMapping(grp%wallDistanceBoundaryNodeArray(j)),&
        grid1FullOld2NewMapping(grp%seedNodeRegister(j))
     end if
    end do
   end if

 ! write trip line stuff

   if(currentGrid==1) then
    ! first count
    localTripMapping = 0
    numberOfLocalTripNodes = 0
    do j=1,grp%numberOfTripNodes
     if(oldToNewNodeMapping(grp%tripLineIndexes(j)).ne.0) then
      numberOfLocalTripNodes = numberOfLocalTripNodes + 1
      localTripMapping(j) = numberOfLocalTripNodes
     end if
    end do

    write(OUTFILE) numberOfLocalTripNodes
    do j=1,grp%numberOfTripNodes
     if(oldToNewNodeMapping(grp%tripLineIndexes(j)).ne.0) then
      write(OUTFILE) oldToNewNodeMapping(grp%tripLineIndexes(j)),grp%tripWallLength(j)
     end if
    end do

    ! first count
    numberOfLocalTripFieldNodes = 0
    do j=1,grp%numberOfTripFieldNodes
     if(oldToNewNodeMapping(grp%tripNodeFieldIndexes(j,1)).ne.0.and.&
        localTripMapping(grp%tripNodeFieldIndexes(j,2)).ne.0) then
      numberOfLocalTripFieldNodes = numberOfLocalTripFieldNodes + 1
     end if
    end do

    write(OUTFILE) numberOfLocalTripFieldNodes
    do j=1,grp%numberOfTripFieldNodes
     if(oldToNewNodeMapping(grp%tripNodeFieldIndexes(j,1)).ne.0.and.&
        localTripMapping(grp%tripNodeFieldIndexes(j,2)).ne.0) then
       write(OUTFILE) oldToNewNodeMapping(grp%tripNodeFieldIndexes(j,1)),&
                      localTripMapping(grp%tripNodeFieldIndexes(j,2)),&
                      grp%tripNodeFieldDistances(j)
     end if
    end do
   end if

 

 ! write plot register

   if(currentGrid==1) then
    write(REGFILE) nodesInDomainCount

    do j=1,grp%numberOfNodes
     if(grp%sdp%parallelizationColorArray(j)==i) then
      write(REGFILE) oldToNewNodeMapping(j),j
     end if
    end do
   end if


  ! write communication arrays

  ! first count

   numberOfSendNodes = 0
   do j=1,grp%numberOfGridDomains
    ghostNodeFlagArray = 0

! Start of OH Modification
!   sideNumber = 0
!   currentSideBlock => grp%sdp%DataBlockChain%first
!   do while(associated(currentSideBlock))
!    do k=1,grp%sdp%dataBlockSize
     do k=1,numberOfInterfaceSides
!     sideNumber = sideNumber+1
!     if(.not.sideIsInternal(sideNumber)) then
!      ind1 = currentSideBlock%block(k)%sideIndexes(1)
!      ind2 = currentSideBlock%block(k)%sideIndexes(2)
       ind1 = CommSidesArray(1,k)
       ind2 = CommSidesArray(2,k)
! End of OH Modification
       domain1 = grp%sdp%parallelizationColorArray(ind1)
       domain2 = grp%sdp%parallelizationColorArray(ind2)
       if(domain1<domain2) then
        if(domain1==i.and.domain2==j) then
        ! write send array
         if(ghostNodeFlagArray(oldToNewNodeMapping(ind2))==0) then
          numberOfSendNodes(domain2) = numberOfSendNodes(domain2) + 1
          ghostNodeFlagArray(oldToNewNodeMapping(ind2)) = 1
         end if
        end if
       else if(domain1>domain2) then
        if(domain2==i.and.domain1==j) then
         ! write send array
         if(ghostNodeFlagArray(oldToNewNodeMapping(ind1))==0) then
          numberOfSendNodes(domain1) = numberOfSendNodes(domain1) + 1
          ghostNodeFlagArray(oldToNewNodeMapping(ind1)) = 1
         end if
        end if
       end if
!     end if
!     if(sideNumber==grp%numberOfSides) exit
     end do
!    currentSideBlock => currentSideBlock%next
!   end do
   end do


   write(OUTFILE) numberOfComSides,numberOfBoundaryComSides

   do j=1,grp%numberOfGridDomains
    write(OUTFILE) numberOfSendNodes(j)
    write(299,*) j,numberOfSendNodes(j)
   end do


   ! send arrays
   do j=1,grp%numberOfGridDomains
    count = 0
    ghostNodeFlagArray = 0
! Start of OH Modification 
!   sideNumber = 0
!   currentSideBlock => grp%sdp%DataBlockChain%first
!   do while(associated(currentSideBlock))
!    do k=1,grp%sdp%dataBlockSize
     do k=1,numberOfInterfaceSides
!     sideNumber = sideNumber+1
!     if(.not.sideIsInternal(sideNumber)) then
!      ind1 = currentSideBlock%block(k)%sideIndexes(1)
!      ind2 = currentSideBlock%block(k)%sideIndexes(2)
       ind1 = CommSidesArray(1,k)
       ind2 = CommSidesArray(2,k)
! End of OH Modification
       domain1 = grp%sdp%parallelizationColorArray(ind1)
       domain2 = grp%sdp%parallelizationColorArray(ind2)
       if(domain1<domain2) then
        if(domain1==i.and.domain2==j) then
        ! write send array
         if(ghostNodeFlagArray(oldToNewNodeMapping(ind2))==0) then
          count = count + 1
          write(OUTFILE)  oldToNewNodeMapping(ind2)
          ghostNodeFlagArray(oldToNewNodeMapping(ind2)) = 1
         end if
        end if
       else if(domain1>domain2) then
        if(domain2==i.and.domain1==j) then
         ! write send array
         if(ghostNodeFlagArray(oldToNewNodeMapping(ind1))==0) then
          count = count + 1
          write(OUTFILE) oldToNewNodeMapping(ind1)
          ghostNodeFlagArray(oldToNewNodeMapping(ind1)) = 1
         end if
        end if
       end if
!     end if
!     if(sideNumber==grp%numberOfSides) exit
     end do
!    currentSideBlock => currentSideBlock%next
!   end do
    write(299,*)j,count
    if(count.ne.numberOfSendNodes(j)) then
     write(*,*) "XX: ",i
     write(*,*) i,j,count,numberOfSendNodes(j)
     STOP "ERROR: Something's wrong in writePar - 7"
    end if
   end do



  ! receive arrays

  ! first count

   numberOfReceiveNodes = 0
   do j=1,grp%numberOfGridDomains
    ghostNodeFlagArray = 0
! Start of OH Modification 
!   sideNumber = 0
!   currentSideBlock => grp%sdp%DataBlockChain%first
!   do while(associated(currentSideBlock))
!    do k=1,grp%sdp%dataBlockSize
     do k=1,numberOfInterfaceSides
!     sideNumber = sideNumber+1
!     if(.not.sideIsInternal(sideNumber)) then
!      ind1 = currentSideBlock%block(k)%sideIndexes(1)
!      ind2 = currentSideBlock%block(k)%sideIndexes(2)
       ind1 = CommSidesArray(1,k)
       ind2 = CommSidesArray(2,k)
! End of OH Modification
       domain1 = grp%sdp%parallelizationColorArray(ind1)
       domain2 = grp%sdp%parallelizationColorArray(ind2)
       if(domain1<domain2) then
        if(domain1==j.and.domain2==i) then
         if(ghostNodeFlagArray(oldToNewNodeMapping(ind2))==0) then
          numberOfReceiveNodes(domain1) = numberOfReceiveNodes(domain1) + 1
          ghostNodeFlagArray(oldToNewNodeMapping(ind2)) = 1
         end if
        end if
       else if(domain1>domain2) then
        if(domain2==j.and.domain1==i) then
         if(ghostNodeFlagArray(oldToNewNodeMapping(ind1))==0) then
          numberOfReceiveNodes(domain2) = numberOfReceiveNodes(domain2) + 1
          ghostNodeFlagArray(oldToNewNodeMapping(ind1)) = 1
         end if
        end if
       end if
!     end if
!     if(sideNumber==grp%numberOfSides) exit
     end do
!    currentSideBlock => currentSideBlock%next
!   end do
   end do

  do j=1,grp%numberOfGridDomains
   write(OUTFILE) numberOfReceiveNodes(j)
   write(299,*) j,numberOfReceiveNodes(j)
!   write(*,*) "Number of receive nodes for domain: ",j,numberOfReceiveNodes(j)
  end do

   do j=1,grp%numberOfGridDomains
    count = 0
    ghostNodeFlagArray = 0
! Start of OH Modification 
!   sideNumber = 0
!   currentSideBlock => grp%sdp%DataBlockChain%first
!   do while(associated(currentSideBlock))
!    do k=1,grp%sdp%dataBlockSize
     do k=1,numberOfInterfaceSides
!     sideNumber = sideNumber+1
!     if(.not.sideIsInternal(sideNumber)) then
!      ind1 = currentSideBlock%block(k)%sideIndexes(1)
!      ind2 = currentSideBlock%block(k)%sideIndexes(2)
       ind1 = CommSidesArray(1,k)
       ind2 = CommSidesArray(2,k)
! End of OH Modification
       domain1 = grp%sdp%parallelizationColorArray(ind1)
       domain2 = grp%sdp%parallelizationColorArray(ind2)
       if(domain1<domain2) then
        if(domain1==j.and.domain2==i) then
         if(ghostNodeFlagArray(oldToNewNodeMapping(ind2))==0) then

          write(OUTFILE) oldToNewNodeMapping(ind2)
          count = count + 1
          ghostNodeFlagArray(oldToNewNodeMapping(ind2)) = 1
         end if
        end if
       else if(domain1>domain2) then
        if(domain2==j.and.domain1==i) then
         if(ghostNodeFlagArray(oldToNewNodeMapping(ind1))==0) then

          write(OUTFILE) oldToNewNodeMapping(ind1)
          count = count + 1
          ghostNodeFlagArray(oldToNewNodeMapping(ind1)) = 1
         end if
        end if
       end if
!     end if
!     if(sideNumber==grp%numberOfSides) exit
     end do
!    currentSideBlock => currentSideBlock%next
!   end do
    write(299,*) j,count
    if(count.ne.numberOfReceiveNodes(j)) then
     write(*,*) i,j,count,numberOfReceiveNodes(j)
     STOP "ERROR: Something's wrong in writePar - 8"
    end if
   end do

  if(currentGrid>1) then

  ! write mapping data

    numberOfFineNonComNodes = 0
    do j=1,finegrp%numberOfNodes
     if(finegrp%sdp%parallelizationColorArray(j)==i.and.(.not.isComNode(j))) then
      numberOfFineNonComNodes = numberOfFineNonComNodes + 1
     end if
    end do

    numberOfLocalComNodes = 0
    do k=1,grp%numberOfFineComNodes
     if(finegrp%sdp%parallelizationColorArray(grp%fineComNodes(k))==i) then
      numberOfLocalComNodes = numberOfLocalComNodes + 1
      localComMapping(k) = numberOfLocalComNodes
     end if
    end do

    count = 0
    write(OUTFILE) numberOfFineNonComNodes
    write(OUTFILE) grp%numberOfIntNodesInPart(i)


    do j=1,finegrp%numberOfNodes
     if(finegrp%sdp%parallelizationColorArray(j)==i.and.(.not.isComNode(j))) then
      count = count + 1
      write(OUTFILE) finegrp%fullOldToNewNodeMapping(j),oldToNewNodeMapping(grp%fineToCoarseIndexMapping(j))
     end if
    end do

! write mapping for com nodes

   write(OUTFILE) numberOfLocalComNodes

   do k=1,grp%numberOfFineComNodes
    if(finegrp%sdp%parallelizationColorArray(grp%fineComNodes(k))==i) then
     write(OUTFILE) finegrp%fullOldToNewNodeMapping(grp%fineComNodes(k))
    end if
   end do


! fine to coarse

! use linear interpolation

    write(OUTFILE) numberOfFineNonComNodes
    write(OUTFILE) grp%numberOfIntNodesInPart(i)

    count = 0 
    do j=1,finegrp%numberOfNodes
     if(finegrp%sdp%parallelizationColorArray(j)==i.and.grp%sdp%parallelizationColorArray(grp%fineToCoarseIndexMapping(j))==i) then
      write(OUTFILE) finegrp%fullOldToNewNodeMapping(j),(finegrp%controlVolumeNodeSize(j)/&
       grp%controlVolumeNodeSize(grp%fineToCoarseIndexMapping(j)))
      count = count + 1
     end if
    end do

! write mapping for com nodes

   write(OUTFILE) numberOfLocalComNodes

   do k=1,grp%numberOfFineComNodes
    if(finegrp%sdp%parallelizationColorArray(grp%fineComNodes(k))==i) then
     write(OUTFILE) (finegrp%controlVolumeNodeSize(grp%fineComNodes(k))/&
                  grp%controlVolumeNodeSize(grp%fineToCoarseIndexMapping(grp%fineComNodes(k))))

    end if
   end do
  end if

  do j=1,grp%numberOfNodes
   if(grp%sdp%parallelizationColorArray(j)==i) then
    grp%fullOldToNewNodeMapping(j) = oldToNewNodeMapping(j)
   end if
  end do

  if(currentGrid>1) then
! write mapping send arrays

   maxNumberOfComNodes = 0
   do j=1,grp%numberOfGridDomains
    numberOfComNodes = 0
    ! first count
    if(i.ne.j) then
     do k=1,grp%numberOfFineComNodes
      if(grp%sdp%parallelizationColorArray(grp%fineToCoarseIndexMapping(grp%fineComNodes(k)))==j.and.&
         finegrp%sdp%parallelizationColorArray(grp%fineComNodes(k))==i) then
        numberOfComNodes = numberOfComNodes + 1
      end if
     end do
     if(numberOfComNodes>maxNumberOfComNodes) maxNumberOfComNodes = numberOfComNodes
    end if
   end do

   write(OUTFILE) maxNumberOfComNodes


   do j=1,grp%numberOfGridDomains
    numberOfComNodes = 0
    ! first count
    if(i.ne.j) then
     do k=1,grp%numberOfFineComNodes
      if(grp%sdp%parallelizationColorArray(grp%fineToCoarseIndexMapping(grp%fineComNodes(k)))==j.and.&
         finegrp%sdp%parallelizationColorArray(grp%fineComNodes(k))==i) then
        numberOfComNodes = numberOfComNodes + 1
      end if
     end do
    end if
    write(OUTFILE) numberOfComNodes
    ! now write
    if(i.ne.j) then
     do k=1,grp%numberOfFineComNodes
      if(grp%sdp%parallelizationColorArray(grp%fineToCoarseIndexMapping(grp%fineComNodes(k)))==j.and.&
         finegrp%sdp%parallelizationColorArray(grp%fineComNodes(k))==i) then
       write(OUTFILE) localComMapping(k)
      end if
     end do
    end if
   end do


! write mapping receive arrays

   maxNumberOfComNodes = 0
   do j=1,grp%numberOfGridDomains
    numberOfComNodes = 0
    ! first count
    if(i.ne.j) then
     do k=1,grp%numberOfFineComNodes
      if(grp%sdp%parallelizationColorArray(grp%fineToCoarseIndexMapping(grp%fineComNodes(k)))==i.and.&
         finegrp%sdp%parallelizationColorArray(grp%fineComNodes(k))==j) then
        numberOfComNodes = numberOfComNodes + 1
      end if
     end do
     if(numberOfComNodes>maxNumberOfComNodes) maxNumberOfComNodes = numberOfComNodes
    end if
   end do

   write(OUTFILE) maxNumberOfComNodes

   do j=1,grp%numberOfGridDomains
    numberOfComNodes = 0
    ! first count
    if(i.ne.j) then
     do k=1,grp%numberOfFineComNodes
      if(grp%sdp%parallelizationColorArray(grp%fineToCoarseIndexMapping(grp%fineComNodes(k)))==i.and.&
         finegrp%sdp%parallelizationColorArray(grp%fineComNodes(k))==j) then
        numberOfComNodes = numberOfComNodes + 1
      end if
     end do
    end if
    write(OUTFILE) numberOfComNodes
    ! now write
    if(i.ne.j) then
     do k=1,grp%numberOfFineComNodes
      if(grp%sdp%parallelizationColorArray(grp%fineToCoarseIndexMapping(grp%fineComNodes(k)))==i.and.&
         finegrp%sdp%parallelizationColorArray(grp%fineComNodes(k))==j) then
       write(OUTFILE) oldToNewNodeMapping(grp%fineToCoarseIndexMapping(grp%fineComNodes(k)))
      end if
     end do
    end if
   end do
  end if

  close(OUTFILE)
 end do
 close(REGFILE)

  if(currentGrid==1)then
   do j=1,grp%numberOfNodes
    grid1FullOld2NewMapping(j)=grp%fullOldToNewNodeMapping(j)
   end do
  end if


  if(grp%gridNumber>1) then
   deallocate(finegrp%fullOldToNewNodeMapping,stat=allocateStatus)
   if (allocateStatus /= 0) STOP "ERROR: writeParallelComputationFiles couldn't deallocate"
   nullify(finegrp%fullOldToNewNodeMapping)
  end if


 deallocate(nodeStartIndex,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeParallelComputationFiles couldn't deallocate"
 deallocate(nodeRegister,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeParallelComputationFiles couldn't deallocate"
 deallocate(ghostNodeFlagArray,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeParallelComputationFiles couldn't deallocate"
 deallocate(numberOfBoundarySendNodes,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeParallelComputationFiles couldn't deallocate"
 deallocate(numberOfBoundaryReceiveNodes,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeParallelComputationFiles couldn't deallocate"
 deallocate(numberOfSendNodes,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeParallelComputationFiles couldn't deallocate"
 deallocate(numberOfReceiveNodes,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeParallelComputationFiles couldn't deallocate"
 deallocate(numberOfDomainSides,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeParallelComputationFiles couldn't deallocate"
 deallocate(numberOfDomainBoundarySides,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeParallelComputationFiles couldn't deallocate"


 end subroutine writeParallelComputationFiles
!-------------------------------------------------------------------------
 subroutine findNodeDistanceFromWall(grp,crp)
 ! finds the shortest distance from the nodes to a wall boundary

 IMPLICIT NONE

 type(GridData) :: grp
 type(CoordinateRegisterData) :: crp

 integer :: i,j,k,indicator,faceIndex,allocateStatus,nodeOnFace
 integer :: count,maxCount,ip,ind,ind1,ind2,i1,i2,ib,ilayer,maxSearchIndex
 integer :: minimumNeighbourInd,numberOfTripNodes,tripLineNumber,minSNode
 integer :: searchNode,ap,fieldSize,minDistInd,nameLength,currentInd
 real :: minimumDistance,dist1,dist2,dist3,xp(3),x1(3),x2(3),tmin,lengthSum
 real :: numberOfNeighbourNodes,n(3),searchCos,maxSearchCos,searchNormal(3),xv(3)
 real :: minimumNeighbourDistance,x(3),t(3),s,smax,minS,p(3),xprev(3)
 real :: searchDistance,searchVector(3),averageNeighbourDistance,numberOfNeighbours
 real :: normalSum(3),numberOfWallNeighbours,minDist,dist,xlen,currentDistance,cdist
 real :: tripRadiusSquared
 logical :: gridIsViscid,hasChanged,uncertainNormal
 integer :: numberOfRefinements,numberOfTripFieldNodes

 real :: coor(3)

 integer :: sideNumber

 character*80 :: distanceFileName

 integer,pointer :: neighbourStartIndex(:),neighbourRegister(:)

 integer,pointer :: tripFieldIndexes(:)
 real,pointer :: tripFieldDistances(:)
 integer,pointer :: tripLineIndexes(:)
 real,pointer :: tripLineDistances(:)

 logical,pointer :: isWallNode(:)
 real,pointer :: turbData(:,:)

 type(SideDataBlockInstance),pointer :: currentSideBlock 

 write(*,*) "Finding node distance from wall..."

 allocate(grp%wallDistance(grp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: findDistanceFromWall out of memory"
 allocate(grp%wallDistanceBoundaryNodeArray(grp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: findDistanceFromWall out of memory"


  ! first make connectivity register

 allocate(neighbourStartIndex(grp%numberOfNodes+1),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeInterpolationDataToFile out of memory"

 write(*,*) "making connectivity register..."



 neighbourStartIndex = 0


 sideNumber = 0
 currentSideBlock => grp%sdp%DataBlockChain%first
 do while(associated(currentSideBlock))
  do i=1,grp%sdp%dataBlockSize
   sideNumber = sideNumber + 1
   neighbourStartIndex(currentSideBlock%block(i)%sideIndexes(1)) =&
    neighbourStartIndex(currentSideBlock%block(i)%sideIndexes(1)) + 1
   neighbourStartIndex(currentSideBlock%block(i)%sideIndexes(2)) =&
    neighbourStartIndex(currentSideBlock%block(i)%sideIndexes(2)) + 1
   if(sideNumber==grp%numberOfSides) exit
  end do
  currentSideBlock=>currentSideBlock%next
 end do



 count = 1
 maxCount = 0
 do i=1,grp%numberOfNodes
  if(maxCount<neighbourStartIndex(i)) maxCount = neighbourStartIndex(i)
 end do


 do i=1,grp%numberOfNodes+1
  count = count + neighbourStartIndex(i)
  neighbourStartIndex(i) = count - neighbourStartIndex(i)
 end do

 allocate(neighbourRegister(neighbourStartIndex(grp%numberOfNodes+1)),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: writeInterpolationDataToFile out of memory"

 allocate(isWallNode(grp%numberOfNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: findDistanceFromWall out of memory"
 isWallNode = .false.

 do i=1,grp%brp%numberOfBoundaryNodes
  if(grp%brp%boundaryNodes(i)%indicator==1) then 
   isWallNode(i) = .true.
  end if
 end do

 neighbourRegister = 0
 sideNumber = 0
 currentSideBlock => grp%sdp%DataBlockChain%first
 do while(associated(currentSideBlock))
  do i=1,grp%sdp%dataBlockSize
   sideNumber = sideNumber + 1  
   i1 = currentSideBlock%block(i)%sideIndexes(1)
   i2 = currentSideBlock%block(i)%sideIndexes(2)
   ind1 = neighbourStartIndex(i1)
   ind2 = neighbourStartIndex(i1+1)-1
   do j=ind1,ind2
    if(neighbourRegister(j)==0) then
     neighbourRegister(j) = i2
     exit
    end if
   end do
   ind1 = neighbourStartIndex(i2)
   ind2 = neighbourStartIndex(i2+1)-1
   do j=ind1,ind2
    if(neighbourRegister(j)==0) then
     neighbourRegister(j) = i1
     exit
    end if
   end do
   if(sideNumber==grp%numberOfSides) exit
  end do
  currentSideBlock=>currentSideBlock%next
 end do
 
 ! now make lines

 grp%wallDistance = -100.0
 grp%wallDistanceBoundaryNodeArray = -1

! ist = grp%brp%nodeIndicatorRegister(-20)
! ien = grp%brp%nodeIndicatorRegister(-19)

! do ib =ist,ien
 do j=1,grp%brp%numberOfBoundaryNodes
  if(isWallNode(j)) then 


   ip = j

   n = grp%brp%boundaryNodes(ip)%wallNormal(1:3)

  ! do a check to find if the normal differs significantly with it's
  ! neighbours

   ind1 = neighbourStartIndex(ip)
   ind2 = neighbourStartIndex(ip+1)-1
   normalSum = 0.0
   numberOfWallNeighbours = 0.0
   uncertainNormal = .false.
   do i=ind1,ind2
    if(neighbourRegister(i).le.grp%brp%numberOfBoundaryNodes) then
     if(isWallNode(neighbourRegister(i))) then
!      if(sum(grp%brp%nodeNormalArray(neighbourRegister(i),:)*n)<0.95) uncertainNormal = .true.
      if(sum(grp%brp%boundaryNodes(neighbourRegister(i))%wallNormal(1:3)*n)<0.25) uncertainNormal = .true.
!      normalSum = normalSum + grp%brp%nodeNormalArray(neighbourRegister(i),:)
      normalSum = normalSum + grp%brp%boundaryNodes(neighbourRegister(i))%wallNormal(1:3)
      numberOfWallNeighbours = numberOfWallNeighbours + 1.0
     end if
    end if
   end do

   if(.not.uncertainNormal) then
    ilayer = ip
   else
    ilayer = -1
   end if

   if(numberOfWallNeighbours>0.1) then
    normalSum = normalSum/numberOfWallNeighbours
   end if

   do while(ilayer>0)
    ind1 = neighbourStartIndex(ilayer)
    ind2 = neighbourStartIndex(ilayer+1)-1

    searchNormal = 0.0
    searchCos = -1.
    maxSearchCos = -1.
    maxSearchIndex = 0

    do i=ind1,ind2
     ! choose node in the direction closest to the surface normal
     ! if it is within a right angle, if not terminate line
!     searchNormal = grp%coordinates(neighbourRegister(i),:)-grp%coordinates(ilayer,:)
     searchNormal = getCoor(neighbourRegister(i),crp)-getCoor(ilayer,crp)
     searchNormal = searchNormal/sqrt(sum(searchNormal*searchNormal))
     searchCos = sum(searchNormal*n)
     if(searchCos>maxSearchCos) then
      maxSearchCos = searchCos
      maxSearchIndex = neighbourRegister(i)
     end if
    end do
    if(maxSearchCos>0.0) then
     ! find distance
!     xv = grp%coordinates(maxSearchIndex,:)-grp%coordinates(ip,:)
     xv = getCoor(maxSearchIndex,crp) - getCoor(ip,crp)
!     cdist = sqrt(sum(xv*xv))
     cdist = sum(xv*xv)
     if(grp%wallDistance(maxSearchIndex)>cdist.or.grp%wallDistance(maxSearchIndex)<0.0) then
      grp%wallDistance(maxSearchIndex) = cdist
      grp%wallDistanceBoundaryNodeArray(maxSearchIndex) = ip
     end if
     ilayer = maxSearchIndex
    else
     ilayer = -1
    end if
   end do
  end if
 end do
 
 do i=1,grp%brp%numberOfBoundaryNodes
  if(isWallNode(i)) then 
   grp%wallDistance(i) = 0.0
   grp%wallDistanceBoundaryNodeArray(i) = i
  end if
 end do  

 ! now loop over all nodes to find average length of neighbours for nodes
 ! that haven't been touched. This might take several passes

 count = 0
 do i=1,grp%numberOfNodes
  if(grp%wallDistance(i)<0.0) count = count + 1
 end do


 write(*,*) "number of undefined nodes: ",count,grp%numberOfNodes

 write(*,*) "checking remains..."

 hasChanged = .true.
 do while(hasChanged)
  hasChanged = .false.
  count = 0
  do i=1,grp%numberOfNodes
   if(grp%wallDistance(i)<-1.0) then
    count = count + 1
    ind1 = neighbourStartIndex(i)
    ind2 = neighbourStartIndex(i+1)-1
    minimumNeighbourDistance = 9.0E20
    averageNeighbourDistance = 0.0
    numberOfNeighbours = 0.0
    minimumNeighbourInd = 0
    do j=ind1,ind2
     if(grp%wallDistance(neighbourRegister(j)).ge.0.0) then
      currentInd = grp%wallDistanceBoundaryNodeArray(neighbourRegister(j))
      xv = getCoor(i,crp)-getCoor(currentInd,crp)
      currentDistance = sum(xv*xv)
      if(currentDistance<minimumNeighbourDistance) then
       minimumNeighbourDistance = currentDistance
       minimumNeighbourInd = currentInd
      end if
     end if
    end do
    if(minimumNeighbourInd>0) then
!     grp%wallDistance(i) = sqrt(minimumNeighbourDistance)
     grp%wallDistance(i) = minimumNeighbourDistance
     grp%wallDistanceBoundaryNodeArray(i) = minimumNeighbourInd
    end if
    hasChanged = .true.
   end if
  end do
!  write(*,*) "number of undefined nodes: ",count
 end do

 write(*,*) "Doing refinement steps: "

! check to see if distance to any of the neighbouring distance nodes is smaller than the current
! distance

! only check twice if distance to a neighbour has been changed

 numberOfRefinements = 1
! grp%wallDistanceBoundaryNodeArray = -grp%wallDistanceBoundaryNodeArray
 do while(numberOfRefinements>0) 
  numberOfRefinements = 0
  do i=1,grp%numberOfNodes
   ind1 = neighbourStartIndex(i)
   ind2 = neighbourStartIndex(i+1)-1
!   hasChanged = .false.
   do j=ind1,ind2   
    currentInd = grp%wallDistanceBoundaryNodeArray(neighbourRegister(j))
!     if(currentInd<0) then 
!     xv = getCoor(i,crp)-getCoor(-currentInd,crp)
     xv = getCoor(i,crp)-getCoor(currentInd,crp)
     currentDistance = sum(xv*xv)
     if(currentDistance<grp%wallDistance(i)) then
      grp%wallDistance(i) = currentDistance
!      grp%wallDistanceBoundaryNodeArray(i) = -currentInd 
      grp%wallDistanceBoundaryNodeArray(i) = currentInd 
      numberOfRefinements = numberOfRefinements + 1
      hasChanged = .true.
     end if
!    end if
   end do
!   if(.not.hasChanged) then 
!    grp%wallDistanceBoundaryNodeArray(i) = abs(grp%wallDistanceBoundaryNodeArray(i))
!   end if
  end do
  write(*,*) "Number of refinements: ",numberOfRefinements
 end do

! take square roots

 do i=1,grp%numberOfNodes
  grp%wallDistanceBoundaryNodeArray(i) = abs(grp%wallDistanceBoundaryNodeArray(i))
  grp%wallDistance(i) = sqrt(grp%wallDistance(i))
 end do

 ! now find trip fields

 write(*,*) "Finding trip lines..."

 ! trip lines are taken as the leading edges

 numberOfTripNodes = 0

 if(grp%gid%numberOfTripLines>0) then
  allocate(tripLineDistances(grp%brp%numberOfBoundaryNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: findNodeDistanceFromWall out of memory"


  ! tripLineIndexes has value 0 if it is not a trip or else trip line number
  allocate(tripLineIndexes(grp%brp%numberOfBoundaryNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: findNodeDistanceFromWall out of memory"
  tripLineIndexes = 0

  do tripLineNumber=1,grp%gid%numberOfTripLines
   x1 = grp%gid%tripLineCoordinates(tripLineNumber,1:3)
   x2 = grp%gid%tripLineCoordinates(tripLineNumber,4:6)
   t = x2-x1

   sMax = sqrt(sum(t*t))
   t = t/sMax

   ! first find field of distances from current trip line to wall boundary

   minSNode = 0
   minS = 99999999.9
   do ib=1,grp%brp%numberOfBoundaryNodes
    if(isWallNode(ib)) then
 !   ip = grp%brp%nodeIndicatorArray(ib)
     ip = ib
 !    x = grp%coordinates(ip,:)
     x = getCoor(ip,crp)
     s = sum((x-x1)*t)
     if(s.ge.0.0.and.s.le.sMax) then
      p = x-x1-s*t
     else
      if(s<0.0) then
       s = 0.0
       p = x-x1
      else
       s = sMax
       p = x-x2
      end if
     end if
     tripLineDistances(ip) = sum(p*p)

     if(s.le.minS) then
      if(minSNode>0) then
       if(tripLineDistances(ip)<tripLineDistances(minSNode)) then
        minS = s
        minSNode = ip
       end if
      else
       minS = s
       minSNode = ip
      end if
     end if
    end if
   end do

  ! start from node closest to the start of the trip line and work our way through edges that are
  ! closest to line and at an acute angle to it

   ilayer = minSNode

 !  xprev = grp%coordinates(ilayer,:)
   xprev = getCoor(ilayer,crp)
   if(tripLineIndexes(ilayer)==0) then
    tripLineIndexes(ilayer) = tripLineNumber
    numberOfTripNodes = numberOfTripNodes + 1
   end if
   do while(ilayer>0)
    ind1 = neighbourStartIndex(ilayer)
    ind2 = neighbourStartIndex(ilayer+1)-1

    searchDistance = 999999999.9
    searchNode = 0
    searchCos = 0.0

    do i=ind1,ind2
     if(isWallNode(neighbourRegister(i))) then
      ! choose closest node to line in its general direction
 !     x = grp%coordinates(neighbourRegister(i),:)
      x = getCoor(neighbourRegister(i),crp)
      searchVector = x-xprev
      if(sum(searchVector*t)>0) then
       s = sum((x-x1)*t)
       if(s.le.sMax) then
        if(tripLineDistances(neighbourRegister(i))<searchDistance) then
         searchDistance = tripLineDistances(neighbourRegister(i))
         searchNode = neighbourRegister(i)
        end if
       end if
      end if
     end if
    end do
    if(searchNode>0) then
 !    xprev = grp%coordinates(searchNode,:)
     xprev = getCoor(searchNode,crp)
     ilayer = searchNode
     if(tripLineIndexes(searchNode)==0) then
      numberOfTripNodes = numberOfTripNodes + 1
      tripLineIndexes(searchNode) = tripLineNumber
     end if
    else
     ! end of line
     ilayer = -1
    end if
   end do
  end do

  deallocate(isWallNode,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: findNodeDistanceFromWall couldn't deallocate"
  deallocate(tripLineDistances,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: findNodeDistanceFromWall couldn't deallocate"
  allocate(grp%tripLineIndexes(numberOfTripNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: findNodeDistanceFromWall out of memory"

  grp%numberOfTripNodes = numberOfTripNodes
  write(*,*) "Number of trip nodes: ",numberOfTripNodes

  count = 0
  do i=1,grp%brp%numberOfBoundaryNodes
   if(tripLineIndexes(i)>0) then
    count = count + 1
    grp%tripLineIndexes(count) = i
   end if
  end do

  if(count.ne.numberOfTripNodes) then
   write(*,*) "ERROR: Something's wrong in findNodeDistanceFromWall - 1: ",count,numberOfTripNodes
   STOP
  end if


  deallocate(tripLineIndexes,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: findNodeDistanceFromWall couldn't deallocate"

  ! now grow layers from each trip node

  allocate(tripFieldIndexes(grp%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: findNodeDistanceFromWall out of memory"
  allocate(tripFieldDistances(grp%numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: findNodeDistanceFromWall out of memory"

  write(*,*) "Finding trip fields..."

  tripRadiusSquared = grp%gid%tripRadius*grp%gid%tripRadius
  tripFieldIndexes = 0
  tripFieldDistances = 0.0
  numberOfTripFieldNodes = 0
  do i=1,grp%numberOfNodes
   minDist = 99999999.9
   minDistInd = 0
   x1 = getCoor(i,crp)
   do j=1,grp%numberOfTripNodes
    x2 = getCoor(grp%tripLineIndexes(j),crp)
    x = x2-x1
    dist = sum(x*x)
    if(dist<minDist) then
     minDist = dist
     minDistInd = j
    end if
   end do
   if(minDist<tripRadiusSquared) then
    tripFieldIndexes(i) = minDistInd
    tripFieldDistances(i) = minDist
    numberOfTripFieldNodes = numberOfTripFieldNodes + 1
   end if
  end do
 
  write(*,*) "trip fields found"
 
  allocate(grp%tripNodeFieldIndexes(numberOfTripFieldNodes,2),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: findNodeDistanceFromWall out of memory"
  allocate(grp%tripNodeFieldDistances(numberOfTripFieldNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: findNodeDistanceFromWall out of memory"
 
 
  grp%tripNodeFieldIndexes = 0.0
  grp%tripNodeFieldDistances = 0.0
 
  count = 0
  do i=1,grp%numberOfNodes
   if(tripFieldIndexes(i)>0) then
    count = count+1
    grp%tripNodeFieldIndexes(count,1) = i
    grp%tripNodeFieldIndexes(count,2) = tripFieldIndexes(i)
    grp%tripNodeFieldDistances(count) = sqrt(tripFieldDistances(i))
   end if
  end do
 
 
  grp%numberOfTripFieldNodes = numberOfTripFieldNodes
 
  deallocate(tripFieldIndexes,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: findNodeDistanceFromWall couldn't deallocate"
  deallocate(tripFieldDistances,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: findNodeDistanceFromWall couldn't deallocate"
 
 
 
  allocate(grp%tripWallLength(numberOfTripNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: findNodeDistanceFromWall out of memory"
 
  ! find trip node wall length
 
  do i=1,numberOfTripNodes
   ind = grp%tripLineIndexes(i)
   lengthSum = 0.
   numberOfNeighbourNodes = 0.
   ind1 = neighbourStartIndex(ind)
   ind2 = neighbourStartIndex(ind+1)-1
 
   do j=ind1,ind2
    if(neighbourRegister(j).le.grp%brp%numberOfBoundaryNodes) then
     x = getCoor(neighbourRegister(j),crp) - getCoor(ind,crp)
     lengthSum = lengthSum + sqrt(sum(x*x))
     numberOfNeighbourNodes = numberOfNeighbourNodes + 1.
    end if
   end do
 
   if(numberOfNeighbourNodes>0.) then
    grp%tripWallLength(i) = lengthSum/numberOfNeighbourNodes
   else
    grp%tripWallLength(i) = 0.1
   end if
  end do
 
 
  deallocate(neighbourStartIndex,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: findNodeDistanceFromWall couldn't deallocate"
  deallocate(neighbourRegister,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: findNodeDistanceFromWall couldn't deallocate"
 
 
  allocate(turbData(grp%numberOfNodes,3),stat=allocateStatus)
  if(allocateStatus/=0) STOP "ERROR: findNodeDistanceFromWall out of memory"
 
 
  turbData(:,1:3) = 0.0
  turbData(:,2) = grp%gid%tripRadius
  do i=1,grp%numberOfTripNodes
   turbData(grp%tripLineIndexes(i),1) = 1.0
  end do
  do i=1,grp%numberOfTripFieldNodes
   turbData(grp%tripNodeFieldIndexes(i,1),2) = grp%tripNodeFieldDistances(i)
  end do
  do i=1,grp%numberOfTripNodes
   turbData(grp%tripLineIndexes(i),3) = grp%tripWallLength(i)
  end do

  open(32,file='turbulenceData.res',form='unformatted')

  write(32) grp%numberOfNodes
  write(32) (1.0,i=1,grp%numberOfNodes),((turbData(i,j),i=1,grp%numberOfNodes),j=1,3),(grp%wallDistance(i),i=1,grp%numberOfNodes)

  close(32)

  deallocate(turbData,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: findNodeDistanceFromWall couldn't deallocate"
 else
  grp%numberOfTripNodes = 0
  grp%numberOfTripFieldNodes = 0
  nullify(grp%tripLineIndexes)
  nullify(grp%tripNodeFieldIndexes)
 end if
 write(*,*) "node distances found..."


end subroutine findNodeDistanceFromWall
!-------------------------------------------------------------------------
 subroutine makeDistanceFieldsForCoarseMesh(finegrp,coarsegrp)
 IMPLICIT NONE
 type(GridData) :: finegrp,coarsegrp

 integer :: i,allocateStatus

 write(*,*) "Finding node distances for coarse mesh..."

 allocate(coarsegrp%wallDistance(coarsegrp%numberOfNodes),stat=allocateStatus)
 if(allocateStatus/=0) STOP "ERROR: makeDistanceFieldsForCoarseMesh out of memory"
 allocate(coarsegrp%wallDistanceBoundaryNodeArray(coarsegrp%numberOfNodes),stat=allocateStatus)
 if(allocateStatus/=0) STOP "ERROR: makeDistanceFieldsForCoarseMesh out of memory"

 coarsegrp%wallDistanceBoundaryNodeArray = 0
 coarsegrp%wallDistance = 0.0

 do i=1,finegrp%numberOfNodes
  coarsegrp%wallDistance(coarsegrp%fineToCoarseIndexMapping(i)) =&
   coarsegrp%wallDistance(coarsegrp%fineToCoarseIndexMapping(i)) +&
   finegrp%controlVolumeNodeSize(i)*finegrp%wallDistance(i)/&
   coarsegrp%controlVolumeNodeSize(coarsegrp%fineToCoarseIndexMapping(i)) 
  if(coarsegrp%wallDistanceBoundaryNodeArray(coarsegrp%fineToCoarseIndexMapping(i))==0) then 
   coarsegrp%wallDistanceBoundaryNodeArray(coarsegrp%fineToCoarseIndexMapping(i))=&
    coarsegrp%fineToCoarseIndexMapping(finegrp%wallDistanceBoundaryNodeArray(i))
  end if
 end do 

 end subroutine makeDistanceFieldsForCoarseMesh
!-----------------------------------------------------------------------
end module Grid 


