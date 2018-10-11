!**********************************************
!*Aggl3d_Boundary.f90                         *
!*                                            *
!*Boundary module definition file             *
!*                                            * 
!*   Aggl3d v 2.3.2                           * 
!*                                            *
!*                                            *
!*24.06.99-                                   *
!**********************************************
!*Description:                                *
!* File contains boundary datastructure and   *
!* agglomeration subroutines                  *
!**********************************************

module TriangularBoundaryFace3D
! module for a side on boundary

 type TriangularBoundaryFace3DData
  integer :: faceIndexes(3)
!  integer,pointer :: oldFaceIndexes(:)
  integer :: surfaceSegmentNumber
  integer :: elementContainingFace
  real :: surfaceNormals(3,3)
  integer,pointer :: indicator
  real :: faceCoefficients(3)
 end type TriangularBoundaryFace3DData
end module TriangularBoundaryFace3D

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

module RectangularBoundaryFace3D
! module for a side on boundary

 type RectangularBoundaryFace3DData
  integer :: faceIndexes(4)
!  integer,pointer :: oldFaceIndexes(:)
  integer :: surfaceSegmentNumber
  real :: surfaceNormals(4,3)
  integer,pointer :: indicator
  real :: faceCoefficients(3)
 end type RectangularBoundaryFace3DData

end module RectangularBoundaryFace3D

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

module BoundaryRegister
! register for boundaries
 use TriangularBoundaryFace3D
 use RectangularBoundaryFace3D
 use SideRegister
 use CoordinateRegister
 use Toolbox

 integer,parameter,private:: NSD = 3

 type BoundarySideData
  integer :: sideIndexes(2) ! mapping between side indexes and point indexes
  real    :: sideCoefficients(3)
  real    :: GCLCoefficients(2)
  real    :: GCLCoefficientsP(2)
  integer :: indicator
 end type BoundarySideData

 type BoundarySideDataPointer
  type(BoundarySideData),pointer :: sd
 end type boundarySideDataPointer

 type LinkedIntegerPointer
  integer :: index
  type(LinkedIntegerPointer), pointer :: next
 end type

 type LinkedTriangularFace
  type(TriangularBoundaryFace3DData),pointer :: bd
  type(LinkedTriangularFace),pointer :: next
 end type LinkedTriangularFace

 type LinkedRectangularFace
  type(RectangularBoundaryFace3DData),pointer :: bd
  type(LinkedRectangularFace),pointer :: next
 end type LinkedRectangularFace

 type BoundaryNodeData
  real :: normal(6)
  real :: wallNormal(6)
  real :: area
  integer :: indicator
 end type BoundaryNodeData

 type BoundaryNodeDataP
  type(BoundaryNodeData),pointer :: first
 end type BoundaryNodeDataP

 type SurfaceSegmentData
  integer,pointer :: nodeIndexes(:)
  integer :: numberOfNodes,indicator,wallboundindicator
 end type SurfaceSegmentData

 type LineSegmentData
  integer,pointer :: indexes(:)
  integer :: numberOfNodes,indicator
 end type LineSegmentData

 type BoundaryRegisterData
  type(TriangularBoundaryFace3DData), pointer :: triangularFaces(:)
  type(RectangularBoundaryFace3DData), pointer :: rectangularFaces(:)
  type(BoundarySideData),pointer :: srp(:)
  type(LineSegmentData),pointer :: lineSegmentRegister(:)
  type(SurfaceSegmentData),pointer :: surfaceSegmentRegister(:)
  type(BoundaryNodeData),pointer :: boundaryNodes(:)
  integer :: boundaryNodesRegister(-20:0)

  integer :: numberOfBoundaryNodes,numberOfTriangularBFs,numberOfRectangularBFs
  integer :: numberOfSurfaceSegments,numberOfLineSegments,numberOfSurfaceNormals
  integer :: numberOfBoundarySides,numberOfBoundaryNormals,numberOfTrailingEdges
  integer :: maxrot

  integer,pointer :: engineInletSideIndexes(:,:)
  real,pointer :: engineInletSideCoefficients(:,:)

  real,pointer :: rotwall(:,:)
  real :: angle


!  real,pointer :: nodeNormals(:,:)

  integer,pointer :: boundaryFaceNodeMappings(:)

  integer,pointer :: boundarySearchRegister(:) ! to speed up boundary creation
  integer,pointer :: boundarySearchRegisterStarter(:)

  logical,pointer :: isNodeTrailingEdge(:)
  integer,pointer :: trailingEdgeNodeRegister(:)
  integer,pointer :: internalOutflowRegister(:)  ! for internal outflow chains
  integer,pointer :: IORegisterLengths(:)
  integer :: numberOfIORegisters,numberOfIONodes
  integer :: numberOfDoubleNodes,numberOfTripleNodes
  integer :: numberOfBSidesCreated
  integer :: numberOfEngineInletSideLimit
  integer :: numberOfEngineInletSides
 end type BoundaryRegisterData

 contains

!------------------------------------------------------------------------
 subroutine cleanUpBoundaryData(brp)
 ! deallocates data for boundary register

 IMPLICIT NONE

 type(BoundaryRegisterData),pointer :: brp

 integer :: i,allocateStatus

  write(*,*) "WW: ",associated(brp%boundaryNodes)

  deallocate(brp%boundaryNodes,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in cleanUpBoundaryData"

  deallocate(brp%srp,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in cleanUpBoundaryData"

  if(.false.) then 
  if(associated(brp%lineSegmentRegister)) then
   do i=1,brp%numberOfLineSegments
    deallocate(brp%lineSegmentRegister(i)%indexes,stat=allocateStatus)
    if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in cleanUpBoundaryData"
   end do
   deallocate(brp%lineSegmentRegister,stat=allocateStatus)
   if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in cleanUpBoundaryData"
   nullify(brp%lineSegmentRegister)
  end if

  if(associated(brp%surfaceSegmentRegister)) then
   do i=1,brp%numberOfSurfaceSegments
    deallocate(brp%surfaceSegmentRegister(i)%nodeIndexes,stat=allocateStatus)
    if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in cleanUpBoundaryData"
   end do
   deallocate(brp%surfaceSegmentRegister,stat=allocateStatus)
   if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in cleanUpBoundaryData"
   nullify(brp%surfaceSegmentRegister)
  end if
  end if

  if(associated(brp%boundaryFaceNodeMappings)) then
   deallocate(brp%boundaryFaceNodeMappings,stat=allocateStatus)
   if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in cleanUpBoundaryData"
   nullify(brp%boundaryFaceNodeMappings)
  end if

  nullify(brp%trailingEdgeNodeRegister) 

  if(associated(brp%internalOutflowRegister)) then
   nullify(brp%internalOutflowRegister)
   nullify(brp%IORegisterLengths)
  end if
 end subroutine cleanUpBoundaryData
!------------------------------------------------------------------------
 subroutine constructBoundaryRegister(numberOfSurfaces,numberOfLines,nTF,nRF,numberOfNodes,&
                                                       numberOfNormals,brp)
! intializes module

  IMPLICIT NONE

  integer :: numberOfSurfaces,numberOfLines,nTF,nRF,numberOfNodes,numberOfNormals
  integer :: numberOfWallnodes
  type(BoundaryRegisterData) :: brp

  integer :: i,allocateStatus

  if(nTF>0) then
   allocate(brp%triangularFaces(nTF),stat=allocateStatus)
   if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"
  else
   nullify(brp%triangularFaces)
  end if

  if(nRF>0) then
   allocate(brp%rectangularFaces(nRF),stat=allocateStatus)
   if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"
  else
   nullify(brp%rectangularFaces)
  end if

  allocate(brp%lineSegmentRegister(numberOfLines),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"

  allocate(brp%surfaceSegmentRegister(numberOfSurfaces),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"

  allocate(brp%isNodeTrailingEdge(numberOfNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"

   do i=1,numberOfNodes
    brp%isNodeTrailingEdge(i) = .false. 
    end do
   print *,'initialised trailing edges: ',numberOfNodes

  do i=1,nTF
   brp%triangularFaces(i)%faceCoefficients = 0.0
  end do

  do i=1,nRF
   brp%rectangularFaces(i)%faceCoefficients = 0.0
  end do


  brp%numberOfSurfaceSegments = numberOfSurfaces
  brp%numberOfLineSegments = numberOfLines
  brp%numberOfBoundaryNodes = numberOfNodes
  brp%numberOfSurfaceNormals = numberOfNormals
  brp%numberOfTriangularBFs = nTF
  brp%numberOfRectangularBFs = nRF
  brp%numberOfBoundarySides = 0
 end subroutine constructBoundaryRegister
!-------------------------------------------------------------------------
 subroutine setUpNewBoundaryData(brp,newbrp)
 ! sets up boundary data for agglomerated grid
 IMPLICIT NONE

 type(BoundaryRegisterData) :: brp,newbrp

 integer :: i,allocateStatus

 nullify(newbrp%triangularFaces)

 nullify(newbrp%rectangularFaces)

!allocate(newbrp%isNodeTrailingEdge(brp%numberOfBoundaryNodes),stat=allocateStatus)
!if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"

 newbrp%numberOfBoundaryNodes = 0

 newbrp%numberOfTriangularBFs = 0
 newbrp%numberOfRectangularBFs = 0
 newbrp%numberOfIORegisters = 0
 newbrp%numberOfIONodes = 0
 newbrp%numberOfTrailingEdges = 0
 nullify(newbrp%IORegisterLengths)
 nullify(newbrp%internalOutflowRegister)
 nullify(newbrp%lineSegmentRegister)
 nullify(newbrp%surfaceSegmentRegister) 
 end subroutine setUpNewBoundaryData
!-------------------------------------------------------------------------
 subroutine registerBoundarySideGCL(p1,p2,r1,r2,vr1,vr2,segmentNumber,brp)
 ! creates or updates boundary side

 integer :: p1,p2,segmentNumber
 double precision :: r1(3,3),r2(3,3),r3(3,3),vr1(3,2),vr2(3,2)
 type(BoundaryRegisterData) :: brp

 type(LinkedInteger),pointer :: current
 integer :: allocateStatus,currentSide,minp
 double precision :: Cxx(3,3),Cyy(3,2),GCL(2,2)

  Cxx(1,:) = 0.5*(r1(2,:)*r2(3,:)-r1(3,:)*r2(2,:))
  Cxx(2,:) = 0.5*(r1(3,:)*r2(1,:)-r1(1,:)*r2(3,:))
  Cxx(3,:) = 0.5*(r1(1,:)*r2(2,:)-r1(2,:)*r2(1,:))

  Cyy(1,1) = 0.25*(r1(2,1)*r2(3,2)-r1(3,1)*r2(2,2))+0.25*(r1(2,2)*r2(3,1)-r1(3,2)*r2(2,1))
  Cyy(2,1) = 0.25*(r1(3,1)*r2(1,2)-r1(1,1)*r2(3,2))+0.25*(r1(3,2)*r2(1,1)-r1(1,2)*r2(3,1))
  Cyy(3,1) = 0.25*(r1(1,1)*r2(2,2)-r1(2,1)*r2(1,2))+0.25*(r1(1,2)*r2(2,1)-r1(2,2)*r2(1,1))

  Cyy(1,2) = 0.25*(r1(2,2)*r2(3,3)-r1(3,2)*r2(2,3))+0.25*(r1(2,3)*r2(3,2)-r1(3,3)*r2(2,2))
  Cyy(2,2) = 0.25*(r1(3,2)*r2(1,3)-r1(1,2)*r2(3,3))+0.25*(r1(3,3)*r2(1,2)-r1(1,3)*r2(3,2))
  Cyy(3,2) = 0.25*(r1(1,2)*r2(2,3)-r1(2,2)*r2(1,3))+0.25*(r1(1,3)*r2(2,2)-r1(2,3)*r2(1,2))

  GCL(1,1) = sum((Cxx(:,1)+Cyy(:,1)+Cxx(:,2))*vr1(:,1))/3.
  GCL(2,1) = sum((Cxx(:,1)+Cyy(:,1)+Cxx(:,2))*vr2(:,1))/3.
  GCL(1,2) = sum((Cxx(:,2)+Cyy(:,2)+Cxx(:,3))*vr1(:,2))/3.
  GCL(2,2) = sum((Cxx(:,2)+Cyy(:,2)+Cxx(:,3))*vr2(:,2))/3.

 if(brp%surfaceSegmentRegister(segmentNumber)%indicator==5.or.brp%surfaceSegmentRegister(segmentNumber)%indicator==6) then
  brp%numberOfEngineInletSideLimit = brp%numberOfEngineInletSideLimit + 1
 end if

 if(p1.eq.p2) write(*,*) "WARNING: Indexes of boundary face are indentical: ",p1,p2
 currentSide =  findBoundarySide(p1,p2,brp)
 if(currentSide>0) then ! side exists
  brp%srp(currentSide)%sideCoefficients = brp%srp(currentSide)%sideCoefficients + 0.25*Cxx(:,1)
  if(p1<p2) then 
   brp%srp(currentSide)%GCLCoefficients(1) = brp%srp(currentSide)%GCLCoefficients(1) + 0.25*GCL(1,1)
   brp%srp(currentSide)%GCLCoefficientsP(1) = brp%srp(currentSide)%GCLCoefficientsP(1) + 0.25*GCL(1,2)
   brp%srp(currentSide)%GCLCoefficients(2) = brp%srp(currentSide)%GCLCoefficients(2) + 0.25*GCL(2,1)
   brp%srp(currentSide)%GCLCoefficientsP(2) = brp%srp(currentSide)%GCLCoefficientsP(2) + 0.25*GCL(2,2)
  else
   brp%srp(currentSide)%GCLCoefficients(1) = brp%srp(currentSide)%GCLCoefficients(1) + 0.25*GCL(2,1)
   brp%srp(currentSide)%GCLCoefficientsP(1) = brp%srp(currentSide)%GCLCoefficientsP(1) + 0.25*GCL(2,2)
   brp%srp(currentSide)%GCLCoefficients(2) = brp%srp(currentSide)%GCLCoefficients(2) + 0.25*GCL(1,1)
   brp%srp(currentSide)%GCLCoefficientsP(2) = brp%srp(currentSide)%GCLCoefficientsP(2) + 0.25*GCL(1,2)
  end if
 else
!  create new side
  brp%numberOfBSidesCreated =  brp%numberOfBSidesCreated + 1
  currentSide =  brp%numberOfBsidesCreated
! insert data into side

  if(p1<p2) then
   brp%srp(currentSide)%sideIndexes(1) = p1
   brp%srp(currentSide)%sideIndexes(2) = p2
  else
   brp%srp(currentSide)%sideIndexes(1) = p2
   brp%srp(currentSide)%sideIndexes(2) = p1
  end if

  brp%srp(currentSide)%sideCoefficients = 0.25*Cxx(:,1)

  if(p1<p2) then 
   brp%srp(currentSide)%GCLCoefficients(1) = 0.25*GCL(1,1)
   brp%srp(currentSide)%GCLCoefficientsP(1) = 0.25*GCL(1,2)
   brp%srp(currentSide)%GCLCoefficients(2) = 0.25*GCL(2,1)
   brp%srp(currentSide)%GCLCoefficientsP(2) = 0.25*GCL(2,2)
  else
   brp%srp(currentSide)%GCLCoefficients(1) = 0.25*GCL(2,1)
   brp%srp(currentSide)%GCLCoefficientsP(1) = 0.25*GCL(2,2)
   brp%srp(currentSide)%GCLCoefficients(2) = 0.25*GCL(1,1)
   brp%srp(currentSide)%GCLCoefficientsP(2) = 0.25*GCL(1,2)
  end if

  brp%srp(currentSide)%indicator = brp%surfaceSegmentRegister(segmentNumber)%indicator

! update search register

  minp = min(p1,p2)
  brp%boundarySearchRegister(brp%numberOfBsidesCreated) = brp%boundarySearchRegisterStarter(minp)
  brp%boundarySearchRegisterStarter(minp) = brp%numberOfBsidesCreated
 end if
 end subroutine registerBoundarySideGCL
!-------------------------------------------------------------------------
 subroutine registerBoundarySide(p1,p2,Cxx,segmentNumber,brp)
 ! creates or updates boundary side

 integer :: p1,p2,segmentNumber
 double precision :: Cxx(3)
 real :: sign
 type(CoordinateRegisterData) :: crp
 type(BoundaryRegisterData) :: brp
 logical :: doVisualization

 type(LinkedInteger),pointer :: current
 integer :: allocateStatus,currentSide,minp
 real :: C12(3)

 C12 = Cxx


 if(p1>p2) then
  sign = -1.
 else
  sign = 1.
 end if

 if(brp%surfaceSegmentRegister(segmentNumber)%indicator==5.or.brp%surfaceSegmentRegister(segmentNumber)%indicator==6) then
  brp%numberOfEngineInletSideLimit = brp%numberOfEngineInletSideLimit + 1
 end if


 if(p1.eq.p2) write(*,*) "WARNING: Indexes of boundary face are indentical: ",p1,p2
 currentSide =  findBoundarySide(p1,p2,brp)
 if(currentSide>0) then ! side exists
  brp%srp(currentSide)%sideCoefficients = brp%srp(currentSide)%sideCoefficients + 0.25*C12
 else
!  create new side
  brp%numberOfBSidesCreated =  brp%numberOfBSidesCreated + 1
  currentSide =  brp%numberOfBsidesCreated
! insert data into side

  if(sign>0) then
   brp%srp(currentSide)%sideIndexes(1) = p1
   brp%srp(currentSide)%sideIndexes(2) = p2
  else
   brp%srp(currentSide)%sideIndexes(1) = p2
   brp%srp(currentSide)%sideIndexes(2) = p1
  end if

  brp%srp(currentSide)%sideCoefficients = 0.25*C12
  brp%srp(currentSide)%indicator = brp%surfaceSegmentRegister(segmentNumber)%indicator

! update search register

  minp = min(p1,p2)
  brp%boundarySearchRegister(brp%numberOfBsidesCreated) = brp%boundarySearchRegisterStarter(minp)
  brp%boundarySearchRegisterStarter(minp) = brp%numberOfBsidesCreated
 end if
 end subroutine registerBoundarySide
!-------------------------------------------------------------------------
 subroutine registerEngineInletSide(p1,p2,Cxx,inletNumber,brp)
 IMPLICIT NONE

 integer :: p1,p2,inletNumber
 double precision :: Cxx(3)
 type(BoundaryRegisterData) :: brp

 double precision :: C12(3)
 integer :: xp1,xp2,i,sideNumber,sp1,sp2

 C12 = Cxx

 if(p1>p2) then
  xp1 = p2
  xp2 = p1
 else
  xp1 = p1
  xp2 = p2
 end if

 ! find out if side exists
 sideNumber = 1
 do i=1,brp%numberOfEngineInletSideLimit
  sp1 = brp%engineInletSideIndexes(i,1)
  sp2 = brp%engineInletSideIndexes(i,2)
  if(sp1>0) then
   if(xp1==sp1.and.xp2==sp2) then
    exit
   else
    sideNumber = sideNumber + 1
   end if
  else
   sideNumber = -sideNumber
   brp%numberOfEngineInletSides = brp%numberOfEngineInletSides + 1
   exit
  end if
 end do

 if(sideNumber>0) then
  ! side exists
   brp%engineInletSideCoefficients(sideNumber,:) = brp%engineInletSideCoefficients(sideNumber,:) + C12
 else
  ! create new side
  sideNumber = -sideNumber
  if(sideNumber.le.brp%numberOfEngineInletSideLimit) then
   brp%engineInletSideIndexes(sideNumber,1) = xp1
   brp%engineInletSideIndexes(sideNumber,2) = xp2
   brp%engineInletSideIndexes(sideNumber,3) = inletNumber
   brp%engineInletSideCoefficients(sideNumber,:) = C12
  else
   STOP "ERROR: Something's wrong in registerEngineInletSide - 1"
  end if
 end if

 end subroutine registerEngineInletSide
!-------------------------------------------------------------------------
 integer function findBoundarySide(p1,p2,brp)
 IMPLICIT NONE

  integer :: p1,p2
  type(BoundaryRegisterData) :: brp

  integer :: i,i1,i2,current

  if(p1<p2) then
   i1 = p1
   i2 = p2
  else
   i1 = p2
   i2 = p1
  end if

  current = brp%boundarySearchRegisterStarter(i1)
  do while(current>0)
   if(i2==brp%srp(current)%sideIndexes(2)) exit
   current = brp%boundarySearchRegister(current)
  end do
  findBoundarySide = current
 end function findBoundarySide
!-------------------------------------------------------------------------
 subroutine makeCalcDataForBoundary(crp,crpPrev,crpPrev2,brp)
 ! makes boundary coefficients for faces

 IMPLICIT NONE

 type(CoordinateRegisterData) :: crp,crpPrev,crpPrev2 ! contains coordinate register for grid
 type(BoundaryRegisterData) :: brp

 integer :: i,p1,p2,p3,p4,allocateStatus,segmentNumber
 double precision :: x1(3,3),x2(3,3),x3(3,3),x4(3,3),xm(3,3),xb1(3,3),xb2(3,3),xb3(3,3),xb4(3,3)
 double precision :: r1(3,3),r2(3,3),vr1(3,2),vr2(3,2)

 allocate(brp%srp(brp%numberOfBoundarySides),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"
 allocate(brp%boundarySearchRegister(brp%numberOfBoundarySides),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"
 allocate(brp%boundarySearchRegisterStarter(brp%numberOfBoundaryNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"
 allocate(brp%srp(brp%numberOfBoundarySides),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"

 brp%boundarySearchRegister = 0
 brp%boundarySearchRegisterStarter = 0


 write(*,*) "Making calculation data for triangles..."
 ! first triangular boundaries
 do i=1,brp%numberOfTriangularBFs
  p1=brp%triangularFaces(i)%faceIndexes(1)
  p2=brp%triangularFaces(i)%faceIndexes(2)
  p3=brp%triangularFaces(i)%faceIndexes(3)
  segmentNumber = brp%triangularFaces(i)%surfaceSegmentNumber

  x1(:,1) = getCoor(p1,crp)
  x2(:,1) = getCoor(p2,crp)
  x3(:,1) = getCoor(p3,crp)
  x1(:,2) = getCoor(p1,crpPrev)
  x2(:,2) = getCoor(p2,crpPrev)
  x3(:,2) = getCoor(p3,crpPrev)
  x1(:,3) = getCoor(p1,crpPrev2)
  x2(:,3) = getCoor(p2,crpPrev2)
  x3(:,3) = getCoor(p3,crpPrev2)

  xm = (x1+x2+x3)/3.

  xb1 = 0.5*(x1+x2)
  xb2 = 0.5*(x2+x3)
  xb3 = 0.5*(x3+x1)

  ! node 1:

  r1 = xm-x1
  r2 = xb1-x1

  vr1(:,1) = x1(:,1) - x1(:,2)
  vr1(:,2) = x1(:,2) - x1(:,3)
  vr1(:,1) = vr1(:,1) + xm(:,1) - xm(:,2)
  vr1(:,2) = vr1(:,2) + xm(:,2) - xm(:,3)
  vr1(:,1) = vr1(:,1) + xb1(:,1) - xb1(:,2)
  vr1(:,2) = vr1(:,2) + xb1(:,2) - xb1(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = x2(:,1) - x2(:,2)
  vr2(:,2) = x2(:,2) - x2(:,3)
  vr2(:,1) = vr2(:,1) + xm(:,1) - xm(:,2)
  vr2(:,2) = vr2(:,2) + xm(:,2) - xm(:,3)
  vr2(:,1) = vr2(:,1) + xb1(:,1) - xb1(:,2)
  vr2(:,2) = vr2(:,2) + xb1(:,2) - xb1(:,3)
  vr2 = vr2/3.0


  call registerBoundarySideGCL(p1,p2,r1,r2,vr1,vr2,segmentNumber,brp)

  r1 = xb3-x1
  r2 = xm-x1

  vr1(:,1) = x1(:,1) - x1(:,2)
  vr1(:,2) = x1(:,2) - x1(:,3)
  vr1(:,1) = vr1(:,1) + xm(:,1) - xm(:,2)
  vr1(:,2) = vr1(:,2) + xm(:,2) - xm(:,3)
  vr1(:,1) = vr1(:,1) + xb3(:,1) - xb3(:,2)
  vr1(:,2) = vr1(:,2) + xb3(:,2) - xb3(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = x3(:,1) - x3(:,2)
  vr2(:,2) = x3(:,2) - x3(:,3)
  vr2(:,1) = vr2(:,1) + xm(:,1) - xm(:,2)
  vr2(:,2) = vr2(:,2) + xm(:,2) - xm(:,3)
  vr2(:,1) = vr2(:,1) + xb3(:,1) - xb3(:,2)
  vr2(:,2) = vr2(:,2) + xb3(:,2) - xb3(:,3)
  vr2 = vr2/3.0


  call registerBoundarySideGCL(p1,p3,r1,r2,vr1,vr2,segmentNumber,brp)

  ! node 2:

  r1 = xm-x2
  r2 = xb2-x2

  vr1(:,1) = x2(:,1) - x2(:,2)
  vr1(:,2) = x2(:,2) - x2(:,3)
  vr1(:,1) = vr1(:,1) + xm(:,1) - xm(:,2)
  vr1(:,2) = vr1(:,2) + xm(:,2) - xm(:,3)
  vr1(:,1) = vr1(:,1) + xb2(:,1) - xb2(:,2)
  vr1(:,2) = vr1(:,2) + xb2(:,2) - xb2(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = x3(:,1) - x3(:,2)
  vr2(:,2) = x3(:,2) - x3(:,3)
  vr2(:,1) = vr2(:,1) + xm(:,1) - xm(:,2)
  vr2(:,2) = vr2(:,2) + xm(:,2) - xm(:,3)
  vr2(:,1) = vr2(:,1) + xb2(:,1) - xb2(:,2)
  vr2(:,2) = vr2(:,2) + xb2(:,2) - xb2(:,3)
  vr2 = vr2/3.0

  call registerBoundarySideGCL(p2,p3,r1,r2,vr1,vr2,segmentNumber,brp)

 end do


 ! now rectangular boundaries

 write(*,*) "Making calculation data for rectangles..."

 do i=1,brp%numberOfRectangularBFs
  p1=brp%rectangularFaces(i)%faceIndexes(1)
  p2=brp%rectangularFaces(i)%faceIndexes(2)
  p3=brp%rectangularFaces(i)%faceIndexes(3)
  p4=brp%rectangularFaces(i)%faceIndexes(4)

!  if(p1==137.or.p2==137.or.p3==137.or.p4==137) write(*,*) "TT: ",i,p1,p2,p3,p4

  segmentNumber = brp%rectangularFaces(i)%surfaceSegmentNumber

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

  xm = 0.5*(x1+x3)

  xb1 = 0.5*(x1+x2)
  xb2 = 0.5*(x2+x3)
  xb3 = 0.5*(x3+x4)
  xb4 = 0.5*(x1+x4)

  ! node 1:
  r1 = xm-x1
  r2 = xb1-x1

  vr1(:,1) = x1(:,1) - x1(:,2)
  vr1(:,2) = x1(:,2) - x1(:,3)
  vr1(:,1) = vr1(:,1) + xm(:,1) - xm(:,2)
  vr1(:,2) = vr1(:,2) + xm(:,2) - xm(:,3)
  vr1(:,1) = vr1(:,1) + xb1(:,1) - xb1(:,2)
  vr1(:,2) = vr1(:,2) + xb1(:,2) - xb1(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = x2(:,1) - x2(:,2)
  vr2(:,2) = x2(:,2) - x2(:,3)
  vr2(:,1) = vr2(:,1) + xm(:,1) - xm(:,2)
  vr2(:,2) = vr2(:,2) + xm(:,2) - xm(:,3)
  vr2(:,1) = vr2(:,1) + xb1(:,1) - xb1(:,2)
  vr2(:,2) = vr2(:,2) + xb1(:,2) - xb1(:,3)
  vr2 = vr2/3.0

  call registerBoundarySideGCL(p1,p2,r1,r2,vr1,vr2,segmentNumber,brp)

  r1 = xb4-x1
  r2 = xm-x1

  vr1(:,1) = x1(:,1) - x1(:,2)
  vr1(:,2) = x1(:,2) - x1(:,3)
  vr1(:,1) = vr1(:,1) + xm(:,1) - xm(:,2)
  vr1(:,2) = vr1(:,2) + xm(:,2) - xm(:,3)
  vr1(:,1) = vr1(:,1) + xb4(:,1) - xb4(:,2)
  vr1(:,2) = vr1(:,2) + xb4(:,2) - xb4(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = x4(:,1) - x4(:,2)
  vr2(:,2) = x4(:,2) - x4(:,3)
  vr2(:,1) = vr2(:,1) + xm(:,1) - xm(:,2)
  vr2(:,2) = vr2(:,2) + xm(:,2) - xm(:,3)
  vr2(:,1) = vr2(:,1) + xb4(:,1) - xb4(:,2)
  vr2(:,2) = vr2(:,2) + xb4(:,2) - xb4(:,3)
  vr2 = vr2/3.0


  call registerBoundarySideGCL(p1,p4,r1,r2,vr1,vr2,segmentNumber,brp)

  ! node 3:

  r1 = xb2-x3
  r2 = xm-x3

  vr1(:,1) = x3(:,1) - x3(:,2)
  vr1(:,2) = x3(:,2) - x3(:,3)
  vr1(:,1) = vr1(:,1) + xm(:,1) - xm(:,2)
  vr1(:,2) = vr1(:,2) + xm(:,2) - xm(:,3)
  vr1(:,1) = vr1(:,1) + xb2(:,1) - xb2(:,2)
  vr1(:,2) = vr1(:,2) + xb2(:,2) - xb2(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = x2(:,1) - x2(:,2)
  vr2(:,2) = x2(:,2) - x2(:,3)
  vr2(:,1) = vr2(:,1) + xm(:,1) - xm(:,2)
  vr2(:,2) = vr2(:,2) + xm(:,2) - xm(:,3)
  vr2(:,1) = vr2(:,1) + xb2(:,1) - xb2(:,2)
  vr2(:,2) = vr2(:,2) + xb2(:,2) - xb2(:,3)
  vr2 = vr2/3.0


  call registerBoundarySideGCL(p3,p2,r1,r2,vr1,vr2,segmentNumber,brp)

  r1 = xm-x3
  r2 = xb3-x3

  vr1(:,1) = x3(:,1) - x3(:,2)
  vr1(:,2) = x3(:,2) - x3(:,3)
  vr1(:,1) = vr1(:,1) + xm(:,1) - xm(:,2)
  vr1(:,2) = vr1(:,2) + xm(:,2) - xm(:,3)
  vr1(:,1) = vr1(:,1) + xb3(:,1) - xb3(:,2)
  vr1(:,2) = vr1(:,2) + xb3(:,2) - xb3(:,3)
  vr1 = vr1/3.0

  vr2(:,1) = x4(:,1) - x4(:,2)
  vr2(:,2) = x4(:,2) - x4(:,3)
  vr2(:,1) = vr2(:,1) + xm(:,1) - xm(:,2)
  vr2(:,2) = vr2(:,2) + xm(:,2) - xm(:,3)
  vr2(:,1) = vr2(:,1) + xb3(:,1) - xb3(:,2)
  vr2(:,2) = vr2(:,2) + xb3(:,2) - xb3(:,3)
  vr2 = vr2/3.0


  call registerBoundarySideGCL(p3,p4,r1,r2,vr1,vr2,segmentNumber,brp)

 end do

 if(brp%numberOfBSidesCreated<brp%numberOfBoundarySides) write(*,*) "there are internal boundary sides"
 brp%numberOfBoundarySides = brp%numberOfBSidesCreated
 write(*,*) "Number of boundary sides: ",brp%numberOfBoundarySides

 deallocate(brp%boundarySearchRegister,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in makeCalcDataForBoundary"
 deallocate(brp%boundarySearchRegisterStarter,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in makeCalcDataForBoundary"
 end subroutine makeCalcDataForBoundary
!-------------------------------------------------------------------------
 subroutine makeEngineInflowData(crp,brp)

 IMPLICIT NONE

 type(CoordinateRegisterData) :: crp ! contains coordinate register for grid
 type(BoundaryRegisterData) :: brp

 integer :: i,p1,p2,p3,p4,allocateStatus,segmentNumber,inletNumber
 double precision :: x1(3),x2(3),x3(3),x4(3),xm(3),xb1(3),xb2(3),xb3(3),xb4(3)
 double precision :: r1(3),r2(3),C12(3),C13(3),C14(3),C23(3),C32(3),C34(3)

 allocate(brp%engineInletSideIndexes(brp%numberOfEngineInletSideLimit,3),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"
 allocate(brp%engineInletSideCoefficients(brp%numberOfEngineInletSideLimit,3),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"

 brp%numberOfEngineInletSides = 0
 brp%engineInletSideCoefficients = 0.0
 brp%engineInletSideIndexes = 0

 write(*,*) "Making engine inlet registers..."
 ! first triangular boundaries
 do i=1,brp%numberOfTriangularBFs
  segmentNumber = brp%triangularFaces(i)%surfaceSegmentNumber
  if(brp%surfaceSegmentRegister(segmentNumber)%indicator==5.or.brp%surfaceSegmentRegister(segmentNumber)%indicator==6) then
   if(brp%surfaceSegmentRegister(segmentNumber)%indicator==5) then
    inletNumber = 1
   else
    inletNumber = 2
   end if

   p1=brp%triangularFaces(i)%faceIndexes(1)
   p2=brp%triangularFaces(i)%faceIndexes(2)
   p3=brp%triangularFaces(i)%faceIndexes(3)

   x1 = getCoor(p1,crp)
   x2 = getCoor(p2,crp)
   x3 = getCoor(p3,crp)

   xm = (x1+x2+x3)/3.

   xb1 = 0.5*(x1+x2)
   xb2 = 0.5*(x2+x3)
   xb3 = 0.5*(x3+x1)

   ! node 1:

   r1 = xm-x1
   r2 = xb1-x1

   C12(1) = 0.5*(r1(2)*r2(3)-r1(3)*r2(2))
   C12(2) = 0.5*(r1(3)*r2(1)-r1(1)*r2(3))
   C12(3) = 0.5*(r1(1)*r2(2)-r1(2)*r2(1))

   call registerEngineInletSide(p1,p2,C12,inletNumber,brp)

   r1 = xb3-x1
   r2 = xm-x1

   C13(1) = 0.5*(r1(2)*r2(3)-r1(3)*r2(2))
   C13(2) = 0.5*(r1(3)*r2(1)-r1(1)*r2(3))
   C13(3) = 0.5*(r1(1)*r2(2)-r1(2)*r2(1))

   call registerEngineInletSide(p1,p3,C13,inletNumber,brp)

   ! node 2:

   r1 = xm-x2
   r2 = xb2-x2

   C23(1) = 0.5*(r1(2)*r2(3)-r1(3)*r2(2))
   C23(2) = 0.5*(r1(3)*r2(1)-r1(1)*r2(3))
   C23(3) = 0.5*(r1(1)*r2(2)-r1(2)*r2(1))

   call registerEngineInletSide(p2,p3,C23,inletNumber,brp)
  end if
 end do


 do i=1,brp%numberOfRectangularBFs
  segmentNumber = brp%rectangularFaces(i)%surfaceSegmentNumber
  if(brp%surfaceSegmentRegister(segmentNumber)%indicator==5.or.brp%surfaceSegmentRegister(segmentNumber)%indicator==6) then
   if(brp%surfaceSegmentRegister(segmentNumber)%indicator==5) then
    inletNumber = 1
   else
    inletNumber = 2
   end if

   p1=brp%rectangularFaces(i)%faceIndexes(1)
   p2=brp%rectangularFaces(i)%faceIndexes(2)
   p3=brp%rectangularFaces(i)%faceIndexes(3)
   p4=brp%rectangularFaces(i)%faceIndexes(4)


   x1 = getCoor(p1,crp)
   x2 = getCoor(p2,crp)
   x3 = getCoor(p3,crp)
   x4 = getCoor(p4,crp)

   xm = 0.25*(x1+x2+x3+x4)  ! perhaps another definition here
!   xm = 0.5*(x1+x3)


   xb1 = 0.5*(x1+x2)
   xb2 = 0.5*(x2+x3)
   xb3 = 0.5*(x3+x4)
   xb4 = 0.5*(x1+x4)

   ! node 1:
   r1 = xm-x1
   r2 = xb1-x1

   C12(1) = 0.5*(r1(2)*r2(3)-r1(3)*r2(2))
   C12(2) = 0.5*(r1(3)*r2(1)-r1(1)*r2(3))
   C12(3) = 0.5*(r1(1)*r2(2)-r1(2)*r2(1))


   call registerEngineInletSide(p1,p2,C12,inletNumber,brp)

   r1 = xb4-x1
   r2 = xm-x1

   C14(1) = 0.5*(r1(2)*r2(3)-r1(3)*r2(2))
   C14(2) = 0.5*(r1(3)*r2(1)-r1(1)*r2(3))
   C14(3) = 0.5*(r1(1)*r2(2)-r1(2)*r2(1))

   call registerEngineInletSide(p1,p4,C14,inletNumber,brp)

   ! node 3:

   r1 = xb2-x3
   r2 = xm-x3

   C32(1) = 0.5*(r1(2)*r2(3)-r1(3)*r2(2))
   C32(2) = 0.5*(r1(3)*r2(1)-r1(1)*r2(3))
   C32(3) = 0.5*(r1(1)*r2(2)-r1(2)*r2(1))

   call registerEngineInletSide(p3,p2,C32,inletNumber,brp)

   r1 = xm-x3
   r2 = xb3-x3

   C34(1) = 0.5*(r1(2)*r2(3)-r1(3)*r2(2))
   C34(2) = 0.5*(r1(3)*r2(1)-r1(1)*r2(3))
   C34(3) = 0.5*(r1(1)*r2(2)-r1(2)*r2(1))

   call registerEngineInletSide(p3,p4,C34,inletNumber,brp)
  end if
 end do

 write(*,*) "Number of engine inlet sides: ",brp%numberOfEngineInletSides

 end subroutine makeEngineInflowData
!-------------------------------------------------------------------------
 subroutine agglomerateBoundary(finebrp,coarsebrp)
 IMPLICIT NONE
 type(BoundaryRegisterData),pointer :: finebrp,coarsebrp

 integer :: currentSide,ind,minind,maxind,ind1,ind2,allocateStatus,numberOfFarfieldNodes(2)
 integer :: i,j,numberOfBNForIndicator(8),nodeIndex,numberOfSidesCreated,blockSize,shiftNumber
 integer :: numberOfTripleNodes,numberOfDoubleNodes,numberOfWallNodes,numberOfSymmetryNodes
 integer :: numberOfEngineInflowNodes(2),numberOfEngineOutflowNodes(2)
 logical :: isWall,isSymmetry
 type(LinkedIntegerBlockChain),pointer :: libc
 real :: NN,N(3)
 integer :: count,fineNodeIndicator,numberOfIONodes,sideNumber 


 allocate(libc,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"
 blockSize = max(coarsebrp%numberOfBoundaryNodes/2,100)
 call LIBCSetUp(coarsebrp%numberOfBoundaryNodes,blockSize,libc)

 ! agglomerates the boundary faces

 allocate(coarsebrp%isNodeTrailingEdge(coarsebrp%numberOfBoundaryNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"

 ! use boundaryFaceNodeMappings array to agglomerate the sides

 coarsebrp%numberOfBoundarySides = 0

 ! first count number of sides - CHANGE THIS

 do i=1,finebrp%numberOfBoundarySides
  ind1 = coarsebrp%boundaryFaceNodeMappings(finebrp%srp(i)%sideIndexes(1))
  ind2 = coarsebrp%boundaryFaceNodeMappings(finebrp%srp(i)%sideIndexes(2))

  if(ind1<ind2) then
   minind = ind1
   maxind = ind2
  else
   minind = ind2
   maxind = ind1
  end if

  ! find out if side has been accounted for
  if(.not.LIBCFindInstance(maxind,minind,libc)) then
   call LIBCAddInstance(maxind,minind,libc)
   coarsebrp%numberOfBoundarySides = coarsebrp%numberOfBoundarySides + 1
  end if
 end do
 call LIBCDestroy(libc)

 write(*,*) "Number of boundary sides in coarse mesh: ",coarsebrp%numberOfBoundarySides

 ! now do the agglomeration

 allocate(coarsebrp%srp(coarsebrp%numberOfBoundarySides),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"
 allocate(coarsebrp%boundarySearchRegister(coarsebrp%numberOfBoundarySides),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"
 allocate(coarsebrp%boundarySearchRegisterStarter(coarsebrp%numberOfBoundaryNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"

 coarsebrp%boundarySearchRegister = 0
 coarsebrp%boundarySearchRegisterStarter = 0
 numberOfSidesCreated = 0
 do i=1,finebrp%numberOfBoundarySides
  ind1 = coarsebrp%boundaryFaceNodeMappings(finebrp%srp(i)%sideIndexes(1))
  ind2 = coarsebrp%boundaryFaceNodeMappings(finebrp%srp(i)%sideIndexes(2))

  currentSide = findBoundarySide(ind1,ind2,coarsebrp)

  if(currentSide>0) then
  ! side already exists
   coarsebrp%srp(currentSide)%sideCoefficients = coarsebrp%srp(currentSide)%sideCoefficients +&
                                               finebrp%srp(i)%sideCoefficients
   coarsebrp%srp(currentSide)%GCLCoefficients = coarsebrp%srp(currentSide)%GCLCoefficients +&
                                               finebrp%srp(i)%GCLCoefficients
   coarsebrp%srp(currentSide)%GCLCoefficientsP = coarsebrp%srp(currentSide)%GCLCoefficientsP +&
                                               finebrp%srp(i)%GCLCoefficientsP

  else
  ! create new side
   numberOfSidesCreated = numberOfSidesCreated + 1
   currentSide = numberOfSidesCreated

   if(ind1<ind2) then
    coarsebrp%srp(currentSide)%sideIndexes(1) = ind1
    coarsebrp%srp(currentSide)%sideIndexes(2) = ind2
   else
    coarsebrp%srp(currentSide)%sideIndexes(1) = ind2
    coarsebrp%srp(currentSide)%sideIndexes(2) = ind1
   end if

   coarsebrp%srp(currentSide)%sideCoefficients  = finebrp%srp(i)%sideCoefficients
   coarsebrp%srp(currentSide)%GCLCoefficients  = finebrp%srp(i)%GCLCoefficients
   coarsebrp%srp(currentSide)%GCLCoefficientsP  = finebrp%srp(i)%GCLCoefficientsP
   coarsebrp%srp(currentSide)%indicator = finebrp%srp(i)%indicator

   ! set up register

   minind = min(ind1,ind2)
   coarsebrp%boundarySearchRegister(numberOfSidesCreated) = coarsebrp%boundarySearchRegisterStarter(minind)
   coarsebrp%boundarySearchRegisterStarter(minind) = numberOfSidesCreated
  end if
 end do
 deallocate(coarsebrp%boundarySearchRegister,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in agglomerateBoundary"
 deallocate(coarsebrp%boundarySearchRegisterStarter,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in agglomerateBoundary"

 ! set up trailing edges

 coarsebrp%isNodeTrailingEdge = .false.
 do i=1,finebrp%numberOfTrailingEdges
  coarsebrp%isNodeTrailingEdge(coarsebrp%boundaryFaceNodeMappings(finebrp%trailingEdgeNodeRegister(i))) = .true.
 end do

 ! count number of trailing edges

 coarsebrp%numberOfTrailingEdges = 0
 do i=1,coarsebrp%numberOfBoundaryNodes
  if(coarsebrp%isNodeTrailingEdge(i)) coarsebrp%numberOfTrailingEdges = coarsebrp%numberOfTrailingEdges + 1
 end do

 write(*,*) "Number of trailing edges in coarse mesh: ",coarsebrp%numberOfTrailingEdges

 if(coarsebrp%numberOfTrailingEdges>0) then
  allocate(coarsebrp%trailingEdgeNodeRegister(coarsebrp%numberOfTrailingEdges),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create boundary array"
 else
  nullify(coarsebrp%trailingEdgeNodeRegister)
 end if
 j = 0
 do i=1,coarsebrp%numberOfBoundaryNodes
  if(coarsebrp%isNodeTrailingEdge(i)) then
   j = j + 1
   coarsebrp%trailingEdgeNodeRegister(j) = i
  end if
 end do

 ! set up engine inflow boundaries

 allocate(coarsebrp%engineInletSideIndexes(finebrp%numberOfEngineInletSides,3),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: agglomerateBoundary out of memory"
 allocate(coarsebrp%engineInletSideCoefficients(finebrp%numberOfEngineInletSides,3),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: agglomerateBoundary out of memory"

 coarsebrp%numberOfEngineInletSides = 0

 do i=1,finebrp%numberOfEngineInletSides
  ind1 = coarsebrp%boundaryFaceNodeMappings(finebrp%engineInletSideIndexes(i,1))
  ind2 = coarsebrp%boundaryFaceNodeMappings(finebrp%engineInletSideIndexes(i,2))
  minInd = min(ind1,ind2)
  maxInd = max(ind1,ind2)

  ! check whether side already exists

  sideNumber = 0
  do j=1,coarsebrp%numberOfEngineInletSides
   if(coarsebrp%engineInletSideIndexes(j,1)==minInd.and.coarsebrp%engineInletSideIndexes(j,2)==maxInd) then
    sideNumber=j
    exit
   end if
  end do

  if(sideNumber>0) then
   coarsebrp%engineInletSideCoefficients(sideNumber,:) = coarsebrp%engineInletSideCoefficients(sideNumber,:) + &
    finebrp%engineInletSideCoefficients(i,:)
  else
   coarsebrp%numberOfEngineInletSides = coarsebrp%numberOfEngineInletSides + 1
   sideNumber = coarsebrp%numberOfEngineInletSides
   coarsebrp%engineInletSideCoefficients(sideNumber,:) = finebrp%engineInletSideCoefficients(i,:)
   coarsebrp%engineInletSideIndexes(sideNumber,1) = minInd
   coarsebrp%engineInletSideIndexes(sideNumber,2) = maxInd
   coarsebrp%engineInletSideIndexes(sideNumber,3) = finebrp%engineInletSideIndexes(i,3)
  end if
 end do

 write(*,*) "Number of engine inlet sides: ",coarsebrp%numberOfEngineInletSides

 ! set up boundaryNodes register

 allocate(coarsebrp%boundaryNodes(coarsebrp%numberOfBoundaryNodes),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: agglomerateBoundary out of memory"

 do i=1,coarsebrp%numberOfBoundaryNodes
  coarsebrp%boundaryNodes(i)%indicator = 0
  coarsebrp%boundaryNodes(i)%normal(1) = 0.0
  coarsebrp%boundaryNodes(i)%normal(2) = 0.0
  coarsebrp%boundaryNodes(i)%normal(3) = 0.0
 end do

 do i=1,finebrp%numberOfBoundaryNodes
  ind = coarsebrp%boundaryFaceNodeMappings(i)
  coarsebrp%boundaryNodes(ind)%normal(1) = coarsebrp%boundaryNodes(ind)%normal(1) + finebrp%boundaryNodes(i)%normal(1)
  coarsebrp%boundaryNodes(ind)%normal(2) = coarsebrp%boundaryNodes(ind)%normal(2) + finebrp%boundaryNodes(i)%normal(2)
  coarsebrp%boundaryNodes(ind)%normal(3) = coarsebrp%boundaryNodes(ind)%normal(3) + finebrp%boundaryNodes(i)%normal(3)
 end do


 ! normalize normals

 do i=1,coarsebrp%numberOfBoundaryNodes
  N(1:3) = coarsebrp%boundaryNodes(i)%normal(1:3)
  NN = sqrt(sum(N*N)) 
  if(NN>1.0e-20) then
   coarsebrp%boundaryNodes(i)%normal(1) = coarsebrp%boundaryNodes(i)%normal(1)/NN
   coarsebrp%boundaryNodes(i)%normal(2) = coarsebrp%boundaryNodes(i)%normal(2)/NN
   coarsebrp%boundaryNodes(i)%normal(3) = coarsebrp%boundaryNodes(i)%normal(3)/NN
  else
    write(*,*) "Warning: Small normal at node: ",i
   coarsebrp%boundaryNodes(i)%normal(1) = 0.0
   coarsebrp%boundaryNodes(i)%normal(2) = 0.0
   coarsebrp%boundaryNodes(i)%normal(3) = 0.0
  end if
 end do
 deallocate(coarsebrp%isNodeTrailingEdge,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in agglomerateBoundary"

 do i=1,finebrp%numberOfBoundaryNodes
  ind = coarsebrp%boundaryFaceNodeMappings(i)
  fineNodeIndicator = finebrp%boundaryNodes(i)%indicator
  if(fineNodeIndicator==0) then 
!  write(*,*) "A: ",i,fineNodeIndicator
!   pause
  end if
  if(fineNodeIndicator.ne.coarsebrp%boundaryNodes(ind)%indicator) then
   if(coarsebrp%boundaryNodes(ind)%indicator>0) then
    if(coarsebrp%boundaryNodes(ind)%indicator==9) then ! lowest priority to inviscid nodes
     coarsebrp%boundaryNodes(ind)%indicator = fineNodeIndicator
    else if(coarsebrp%boundaryNodes(ind)%indicator==2.and.fineNodeIndicator.ne.9) then ! lowest priority to symmetry nad inviscid nodes
!   if(coarsebrp%boundaryNodes(ind)%indicator==2) then ! lowest priority to symmetry nodes
     coarsebrp%boundaryNodes(ind)%indicator = fineNodeIndicator 
    else if(coarsebrp%boundaryNodes(ind)%indicator==1.and.fineNodeIndicator>2) then
    ! outer and engine nodes are of higher priority than wall
     coarsebrp%boundaryNodes(ind)%indicator = fineNodeIndicator 
    else if(coarsebrp%boundaryNodes(ind)%indicator.le.4.and.fineNodeIndicator>4) then
    ! engine nodes of higher priority than outer
     coarsebrp%boundaryNodes(ind)%indicator = fineNodeIndicator 
    end if
   else
    coarsebrp%boundaryNodes(ind)%indicator = fineNodeIndicator 
   end if
  end if
 end do

 ! set boundaryNodesChainsRegister
 numberOfTripleNodes = 0
 numberOfDoubleNodes = 0
 numberOfWallNodes = 0
 numberOfSymmetryNodes = 0
 numberOfFarfieldNodes = 0
 numberOfEngineInflowNodes = 0
 numberOfEngineOutflowNodes = 0
 numberOfIONodes = 0
 do i=1,coarsebrp%numberOfBoundaryNodes
  ind = coarsebrp%boundaryNodes(i)%indicator
  if(ind>99) then
   numberOfTripleNodes = numberOfTripleNodes + 1
  else if(ind>9) then
   numberOfDoubleNodes = numberOfDoubleNodes + 1
  else if(ind==1) then
   numberOfWallNodes = numberOfWallNodes + 1
  else if(ind==2) then
   numberOfSymmetryNodes = numberOfSymmetryNodes + 1
  else if(ind==3) then
   numberOfFarfieldNodes(1) = numberOfFarfieldNodes(1) + 1
  else if(ind==4) then
   numberOfFarfieldNodes(2) = numberOfFarfieldNodes(2) + 1
  else if(ind==5) then
   numberOfEngineInflowNodes(1) = numberOfEngineInflowNodes(1) + 1
  else if(ind==6) then
   numberOfEngineInflowNodes(2) = numberOfEngineInflowNodes(2) + 1
  else if(ind==7) then
   numberOfEngineOutflowNodes(1) = numberOfEngineOutflowNodes(1) + 1
  else if(ind==8) then
   numberOfEngineOutflowNodes(2) = numberOfEngineOutflowNodes(2) + 1
  else if(ind==9) then 
   numberOfIONodes = numberOfIONodes + 1
  end if
 end do


 write(*,*) "wall nodes: ",numberOfWallNodes
 write(*,*) "symmetry nodes: ",numberOfSymmetryNodes
 write(*,*) "farfield nodes: ",numberOfFarfieldNodes
 write(*,*) "engine inflow nodes: ",numberOfEngineInflowNodes
 write(*,*) "engine outflow nodes: ",numberOfEngineOutflowNodes

 coarsebrp%boundaryNodesRegister(-20:-1) = 0
 coarsebrp%boundaryNodesRegister(-20) = 1
 coarsebrp%boundaryNodesRegister(-19) = numberOfWallNodes

 shiftNumber = numberOfWallNodes + 1

 coarsebrp%boundaryNodesRegister(-18) = shiftNumber
 coarsebrp%boundaryNodesRegister(-17) = shiftNumber + numberOfSymmetryNodes - 1

 shiftNumber = shiftNumber + numberOfSymmetryNodes

 coarsebrp%boundaryNodesRegister(-16) = shiftNumber
 coarsebrp%boundaryNodesRegister(-15) = shiftNumber + numberOfFarfieldNodes(1) - 1

 shiftNumber = shiftNumber + numberOfFarfieldNodes(1)

 coarsebrp%boundaryNodesRegister(-14) = shiftNumber
 coarsebrp%boundaryNodesRegister(-13) = shiftNumber + numberOfFarfieldNodes(2) - 1

 shiftNumber = shiftNumber + numberOfFarfieldNodes(2)

 coarsebrp%boundaryNodesRegister(-12) = shiftNumber
 coarsebrp%boundaryNodesRegister(-11) = shiftNumber + numberOfEngineInflowNodes(1) - 1

 shiftNumber = shiftNumber + numberOfEngineInflowNodes(1)

 coarsebrp%boundaryNodesRegister(-10) = shiftNumber
 coarsebrp%boundaryNodesRegister(-9) = shiftNumber + numberOfEngineInflowNodes(2) - 1

 shiftNumber = shiftNumber + numberOfEngineInflowNodes(2)

 coarsebrp%boundaryNodesRegister(-8) = shiftNumber
 coarsebrp%boundaryNodesRegister(-7) = shiftNumber + numberOfEngineOutflowNodes(1) - 1

 shiftNumber = shiftNumber + numberOfEngineOutflowNodes(1)

 coarsebrp%boundaryNodesRegister(-6) = shiftNumber
 coarsebrp%boundaryNodesRegister(-5) = shiftNumber + numberOfEngineOutflowNodes(2) - 1

 shiftNumber = shiftNumber + numberOfEngineOutflowNodes(2)

 coarsebrp%boundaryNodesRegister(-4) = shiftNumber
 coarsebrp%boundaryNodesRegister(-3) = shiftNumber + numberOfIONodes - 1

 coarsebrp%numberOfBoundaryNormals = 0

 coarsebrp%numberOfTripleNodes = numberOfTripleNodes
 coarsebrp%numberOfDoubleNodes = numberOfDoubleNodes

 end subroutine agglomerateBoundary
!-------------------------------------------------------------------------
 subroutine smoothNormals(brp,factor)
 IMPLICIT NONE

 type(BoundaryRegisterData) :: brp
 real :: factor

 real :: f1,f2,N1(3),N2(3),N(3),NN
 integer :: i,i1,i2,ind1,ind2

 f1 = factor
 f2 = 1.0 - factor
 do i=1,brp%numberOfBoundarySides
  i1 = brp%srp(i)%sideIndexes(1)  
  i2 = brp%srp(i)%sideIndexes(2)  

  ind1 = brp%boundaryNodes(i1)%indicator
  ind2 = brp%boundaryNodes(i2)%indicator
  if(ind1.le.9.and.ind2.le.9.and.ind1>0.and.ind2>0) then 
!  if(ind1.eq.1.and.ind2.eq.1) then 
   N1(1:3) = brp%boundaryNodes(i1)%normal(1:3)
   N2(1:3) = brp%boundaryNodes(i2)%normal(1:3)
   if(ind1.eq.1) then 
    brp%boundaryNodes(i1)%normal(1:3) = N1 + f1*N2
   end if
   if(ind2.eq.1) then 
    brp%boundaryNodes(i2)%normal(1:3) = N2 + f1*N1
   end if
  end if
 end do

 do i=1,brp%numberOfBoundaryNodes
  if(brp%boundaryNodes(i)%indicator==1) then 
   N(1:3) = brp%boundaryNodes(i)%normal(1:3)
   NN = sqrt(sum(N*N)) 
   if(NN>1.0e-10) then 
    brp%boundaryNodes(i)%normal(1:3) = N/NN
   else
     brp%boundaryNodes(i)%normal(1:3) =0.0
   end if
  end if
 end do
  
 end subroutine smoothNormals
!-------------------------------------------------------------------------
 subroutine findBoundaryNodeArea(brp)
 IMPLICIT NONE

 type(BoundaryRegisterData) :: brp

 integer :: i,ind1,ind2
 real :: areaIncrement

 brp%boundaryNodes(:)%area = 0.0

  do i=1,brp%numberOfBoundarySides
   ind1 = brp%srp(i)%sideIndexes(1)
   ind2 = brp%srp(i)%sideIndexes(2)
   areaIncrement = sqrt(sum(brp%srp(i)%sideCoefficients*brp%srp(i)%sideCoefficients))
   brp%boundaryNodes(ind1)%area = brp%boundaryNodes(ind1)%area + areaIncrement
   brp%boundaryNodes(ind2)%area = brp%boundaryNodes(ind2)%area + areaIncrement
  end do

 end subroutine findBoundaryNodeArea
!-------------------------------------------------------------------------
 subroutine readBoundData(nbound,INFILE1,INFILE2,brp,isHybrid,groundAngle)
! Reads boundary data form file
 IMPLICIT NONE

 integer :: nbound,INFILE1,INFILE2
 type(BoundaryRegisterData) :: brp

 integer :: ip,jp,allocateStatus
 integer :: nbd,ind,i,j,dummy,if
 integer :: segmentNumber,nodesInSegment,nodesInUDirection,nodesInVDirection
 real :: dummyCoor,groundAngle
 real :: rotwalldata(10,6)
 character :: buff*10
 logical :: isHybrid

 integer :: i1,i2,i3,dummyi

 ! read surface and line segment indicators

 read(INFILE1,*) buff
 brp%maxrot=0
 do i=1,brp%numberOfSurfaceSegments
  read(INFILE1,*) ind,brp%surfaceSegmentRegister(ind)%indicator,&
            brp%surfaceSegmentRegister(ind)%wallboundindicator
  if(brp%surfaceSegmentRegister(ind)%wallboundindicator>brp%maxrot)then
    brp%maxrot=brp%maxrot+1
    write(*,*) 'Surf =',ind,'maxrot=',brp%maxrot
  endif
  if(brp%surfaceSegmentRegister(ind)%indicator>9)&
    write(*,*) "indicator out of bounds: ",brp%surfaceSegmentRegister(ind)%indicator
  if(brp%surfaceSegmentRegister(ind)%indicator>9)&
    write(*,*) "wall boundary indicator out of bounds: ",brp%surfaceSegmentRegister(ind)%wallboundindicator
 end do
    write(*,*)
    write(*,*)'number of rotating objects = ',brp%maxrot
    write(*,*)
 read(INFILE1,*) buff
 do i=1,brp%numberOfLineSegments
  read(INFILE1,*) ind,brp%lineSegmentRegister(ind)%indicator
 end do
 ! read face indexes and surface segment number

 if(brp%maxrot>0)then
  allocate(brp%rotwall(brp%maxrot,6),stat=allocateStatus)
  if(allocateStatus/=0)STOP'ERROR: Not enough memory to create rotating object array'
 endif
 if(brp%maxrot>0) read(INFILE1,*) buff
 do i=1,brp%maxrot
  read(INFILE1,*) ind,(rotwalldata(i,j),j=1,6)
  write(*,*) 'rotating object',i
  write(*,*) 'rotation coordinates:',rotwalldata(i,1),rotwalldata(i,2),rotwalldata(i,3)
  write(*,*) 'rotation vector:',rotwalldata(i,4),rotwalldata(i,5),rotwalldata(i,6)
 end do
 
 do i=1,brp%maxrot
  do j=1,6
   brp%rotwall(i,j)=rotwalldata(i,j)
  end do
 end do

 brp%angle = groundAngle

 write(*,*) "R1: ",brp%numberOfRectangularBFs

 if(isHybrid.or.brp%numberOfRectangularBFs>0) then 
  read(INFILE2) ((brp%rectangularFaces(if)%faceIndexes(i),if=1,brp%numberOfRectangularBFs),i=1,4),&
                 (brp%rectangularFaces(if)%surfaceSegmentNumber,if=1,brp%numberOfRectangularBFs)
 end if

 write(*,*) "R2: ",brp%numberOfTriangularBFs
 if(brp%numberOfTriangularBFs>0) then 
  read(INFILE2) ((brp%triangularFaces(if)%faceIndexes(i),if=1,brp%numberOfTriangularBFs),i=1,3),&
   (brp%triangularFaces(if)%elementContainingFace,if=1,brp%numberOfTriangularBFs),&
   (brp%triangularFaces(if)%surfaceSegmentNumber,if=1,brp%numberOfTriangularBFs)
 end if 
 PRINT*,'Triangle 934414 surfaceSegmentNumber:',brp%triangularFaces(934414)%surfaceSegmentNumber
 PRINT*,'Triangle 934415 surfaceSegmentNumber:',brp%triangularFaces(934415)%surfaceSegmentNumber

 end subroutine readBoundData
!-------------------------------------------------------------------------
 integer function findNumberOfSurfaceNormals(INFILE,boundaryNodes,numberOfBFs,lineSegments,surfaceSegments)
 IMPLICIT NONE

 integer :: INFILE,boundaryNodes,numberOfBFs,lineSegments,surfaceSegments

 integer :: dummy,i,ip,jp,nodesInSegment
 real :: dummyr

 do i=1,boundaryNodes
  read(INFILE,*) dummy
 end do

 do i=1,numberOfBFs
  read(INFILE,*) dummy
 end do

 ! read line segments
 ip = 0
 do i=1,lineSegments
  read(INFILE,*) dummy,nodesInSegment ! at the beginning of each segment list, segment number is stored
  ip = ip + 1
  read(INFILE,*) (dummy,dummyr,jp=ip+1,ip+nodesInSegment)
  ip = ip + nodesInSegment + 1
 end do

 ! read surface normals and node indexes related to them

 ip = 0
 do i=1,surfaceSegments
  read(INFILE,*) dummy, nodesInSegment
  ip = ip + 1
  findNumberOfSurfaceNormals = findNumberOfSurfaceNormals + nodesInSegment
  read(INFILE,*) (dummy,dummyr,dummyr,jp=ip+1,ip+nodesInSegment)
  ip = ip + nodesInSegment
 end do

 rewind(INFILE)
 read(INFILE,*) dummy

 end function findNumberOfSurfaceNormals
!-------------------------------------------------------------------------
 subroutine findTrailingEdges(brp,crp)
 IMPLICIT NONE

 type(BoundaryRegisterData) :: brp
 type(CoordinateRegisterData) :: crp

 integer :: i,j,i1,i2,i3,i4,allocateStatus
 real :: N(3),NN,r1(3),r2(3)

 real :: alpha

 real,pointer :: boundaryNormals(:,:)

 integer :: numberOfTrailingEdges

 alpha = -0.65  ! trailing edge angle is 130 degrees

print * ,' In find trailing edges'

 allocate(boundaryNormals(brp%numberOfBoundaryNodes,3),stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: Not enough memory in findTrailingEdges"
print * ,' allocated Boundary normals'

print * ,' initialised Boundary normals', brp%numberOfBoundaryNodes
 do i=1,brp%numberOfBoundaryNodes
  boundaryNormals(i,:) = 0.0
  brp%isNodeTrailingEdge(i) = .false.
 end do
print * ,' initialised Boundary normals', brp%numberOfBoundaryNodes

 do j=1,brp%numberOfTriangularBFs
  i1 = brp%triangularFaces(j)%faceIndexes(1)
  i2 = brp%triangularFaces(j)%faceIndexes(2)
  i3 = brp%triangularFaces(j)%faceIndexes(3)

  r1 = getCoor(i2,crp)-getCoor(i1,crp)
  r2 = getCoor(i3,crp)-getCoor(i1,crp)
  N(1) = r1(2)*r2(3)-r1(3)*r2(2)
  N(2) = r1(3)*r2(1)-r1(1)*r2(3)
  N(3) = r1(1)*r2(2)-r1(2)*r2(1)


  NN = sqrt(sum(N*N))
  n = N/NN
  NN = sqrt(sum(boundaryNormals(i1,:)*boundaryNormals(i1,:)))
  if(sum(boundaryNormals(i1,:)*boundaryNormals(i1,:))>0.1) then
   if(sum(boundaryNormals(i1,:)*n)<alpha*NN) then
    brp%isNodeTrailingEdge(i1) = .true.
   end if
  end if
  boundaryNormals(i1,:) = boundaryNormals(i1,:) + n
  NN = sqrt(sum(boundaryNormals(i2,:)*boundaryNormals(i2,:)))
  if(sum(boundaryNormals(i2,:)*boundaryNormals(i2,:))>0.1) then
   if(sum(boundaryNormals(i2,:)*n)<alpha*NN) then
    brp%isNodeTrailingEdge(i2) = .true.
   end if
  end if
  boundaryNormals(i2,:) = boundaryNormals(i2,:) + n
  NN = sqrt(sum(boundaryNormals(i3,:)*boundaryNormals(i3,:)))
  if(sum(boundaryNormals(i3,:)*boundaryNormals(i3,:))>0.1) then
   if(sum(boundaryNormals(i3,:)*n)<alpha*NN) then
    brp%isNodeTrailingEdge(i3) = .true.
   end if
  end if
  boundaryNormals(i3,:) = boundaryNormals(i3,:) + n

 end do

 do j=1,brp%numberOfRectangularBFs
  i1 = brp%rectangularFaces(j)%faceIndexes(1)
  i2 = brp%rectangularFaces(j)%faceIndexes(2)
  i3 = brp%rectangularFaces(j)%faceIndexes(3)
  i4 = brp%rectangularFaces(j)%faceIndexes(4)

  r1 = getCoor(i2,crp)-getCoor(i1,crp)
  r2 = getCoor(i4,crp)-getCoor(i1,crp)
  N(1) = r1(2)*r2(3)-r1(3)*r2(2)
  N(2) = r1(3)*r2(1)-r1(1)*r2(3)
  N(3) = r1(1)*r2(2)-r1(2)*r2(1)

  r1 = getCoor(i3,crp)-getCoor(i2,crp)
  r2 = getCoor(i1,crp)-getCoor(i2,crp)
  N(1) = N(1) + r1(2)*r2(3)-r1(3)*r2(2)
  N(2) = N(2) + r1(3)*r2(1)-r1(1)*r2(3)
  N(3) = N(3) + r1(1)*r2(2)-r1(2)*r2(1)

  r1 = getCoor(i4,crp)-getCoor(i3,crp)
  r2 = getCoor(i2,crp)-getCoor(i3,crp)
  N(1) = N(1) + r1(2)*r2(3)-r1(3)*r2(2)
  N(2) = N(2) + r1(3)*r2(1)-r1(1)*r2(3)
  N(3) = N(3) + r1(1)*r2(2)-r1(2)*r2(1)

  r1 = getCoor(i1,crp)-getCoor(i4,crp)
  r2 = getCoor(i3,crp)-getCoor(i4,crp)
  N(1) = N(1) + r1(2)*r2(3)-r1(3)*r2(2)
  N(2) = N(2) + r1(3)*r2(1)-r1(1)*r2(3)
  N(3) = N(3) + r1(1)*r2(2)-r1(2)*r2(1)

  NN = sqrt(sum(N*N))
  n = N/NN


  if(sum(boundaryNormals(i1,:)*boundaryNormals(i1,:))>0.1) then
   if(sum(boundaryNormals(i1,:)*n)>alpha) then
    brp%isNodeTrailingEdge(i1) = .true.
   end if
  end if
  boundaryNormals(i1,:) = boundaryNormals(i1,:) + n
  if(sum(boundaryNormals(i2,:)*boundaryNormals(i2,:))>0.1) then
   if(sum(boundaryNormals(i2,:)*n)>alpha) then
    brp%isNodeTrailingEdge(i2) = .true.
   end if
  end if
  boundaryNormals(i2,:) = boundaryNormals(i2,:) + n
  if(sum(boundaryNormals(i3,:)*boundaryNormals(i3,:))>0.1) then
   if(sum(boundaryNormals(i3,:)*n)>alpha) then
    brp%isNodeTrailingEdge(i3) = .true.
   end if
  end if
  boundaryNormals(i3,:) = boundaryNormals(i3,:) + n
  if(sum(boundaryNormals(i4,:)*boundaryNormals(i4,:))>0.1) then
   if(sum(boundaryNormals(i4,:)*n)>alpha) then
    brp%isNodeTrailingEdge(i4) = .true.
   end if
  end if
  boundaryNormals(i4,:) = boundaryNormals(i4,:) + n
 end do

 deallocate(boundaryNormals,stat=allocateStatus)
 if (allocateStatus /= 0) STOP "ERROR: Couldn't deallocate in findTrailingEdges"

 numberOfTrailingEdges = 0
 do i=1,brp%numberOfBoundaryNodes
  if(brp%isNodeTrailingEdge(i)) numberOfTrailingEdges = numberOfTrailingEdges + 1
 end do

  write(*,*) "number of trailing edges found: ",numberOfTrailingEdges

 end subroutine findTrailingEdges
!
!-----------------------------------------------------------------------------------
 subroutine setBoundIndForViscousMesh(brp)
 ! sets boundary indicators for a viscous mesh (where the normals aren't calculated) 
 IMPLICIT NONE

 type(BoundaryRegisterData) :: brp

 integer :: i,j,surfaceNumber,surfaceIndicator,ind,shiftNumber,allocateStatus
 integer :: numberOfWallNodes,numberOfSymmetryNodes,numberOfFarfieldNodes(2)
 integer :: numberOfEngineInflowNodes(2),numberOfEngineOutflowNodes(2)
 integer :: numberOfInternalOutflowNodes,numberOfInviscidWallNodes

  allocate(brp%boundaryNodes(brp%numberOfBoundaryNodes),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: setBoundIndForViscousMesh out of memory"
  do i=1,brp%numberOfBoundaryNodes
   brp%boundaryNodes(i)%indicator = 0
  end do

  do i=1,brp%numberOfTriangularBFs
   surfaceNumber = brp%triangularFaces(i)%surfaceSegmentNumber
   surfaceIndicator = brp%surfaceSegmentRegister(surfaceNumber)%indicator
   do j=1,3
    ind = brp%triangularFaces(i)%faceIndexes(j)
    if(surfaceIndicator>6) brp%isNodeTrailingEdge(ind) = .false.
    if(.not.brp%isNodeTrailingEdge(ind)) then
      if(surfaceIndicator.ne.brp%boundaryNodes(ind)%indicator) then
        if(brp%boundaryNodes(ind)%indicator>0) then
          if(brp%boundaryNodes(ind)%indicator==9) then ! lowest priority to Inviscid nodes 
            brp%boundaryNodes(ind)%indicator = surfaceIndicator
          else if(brp%boundaryNodes(ind)%indicator==2.and.surfaceIndicator.ne.9) then ! lowest priority to symmetry nodes 
            brp%boundaryNodes(ind)%indicator = surfaceIndicator
          else if(brp%boundaryNodes(ind)%indicator==1.and.surfaceIndicator>2) then
          ! outer and engine nodes are of higher priority than wall
            brp%boundaryNodes(ind)%indicator = surfaceIndicator
          else if(brp%boundaryNodes(ind)%indicator.le.4.and.surfaceIndicator>4) then
          ! engine nodes of higher priority than outer
           brp%boundaryNodes(ind)%indicator = surfaceIndicator
          end if
        else
          brp%boundaryNodes(ind)%indicator = surfaceIndicator
        end if
      end if
    else
      brp%boundaryNodes(ind)%indicator = 0
    end if
   end do
  end do

  do i=1,brp%numberOfRectangularBFs
   surfaceNumber = brp%rectangularFaces(i)%surfaceSegmentNumber
   surfaceIndicator = brp%surfaceSegmentRegister(surfaceNumber)%indicator
   do j=1,4
    ind = brp%rectangularFaces(i)%faceIndexes(j)
    if(surfaceIndicator>6) brp%isNodeTrailingEdge(ind) = .false.
    if(.not.brp%isNodeTrailingEdge(ind)) then
      if(surfaceIndicator.ne.brp%boundaryNodes(ind)%indicator) then
        if(brp%boundaryNodes(ind)%indicator>0) then
          if(brp%boundaryNodes(ind)%indicator==9) then ! lowest priority to Inviscid nodes 
            brp%boundaryNodes(ind)%indicator = surfaceIndicator
          else if(brp%boundaryNodes(ind)%indicator==2.and.surfaceIndicator.ne.9) then ! lowest priority to symmetry nodes
            brp%boundaryNodes(ind)%indicator = surfaceIndicator
          else if(brp%boundaryNodes(ind)%indicator==1.and.surfaceIndicator>2) then
           ! outer and engine nodes are of higher priority than wall
            brp%boundaryNodes(ind)%indicator = surfaceIndicator
          else if(brp%boundaryNodes(ind)%indicator.le.4.and.surfaceIndicator>4) then
           ! engine nodes of higher priority than outer
            brp%boundaryNodes(ind)%indicator = surfaceIndicator
          end if
        else
          brp%boundaryNodes(ind)%indicator = surfaceIndicator
        end if
      end if
    else
      brp%boundaryNodes(ind)%indicator = 0
    end if
   end do
  end do


 brp%numberOfTrailingEdges = 0
 do i = 1 , brp%numberOfBoundaryNodes
  if(brp%isNodeTrailingEdge(i)) then
    brp%numberOfTrailingEdges = brp%numberOftrailingEdges + 1
  end if
 end do

 write(*,*) ' Number of trailing edges: ',brp%numberOfTrailingEdges

  allocate(brp%trailingEdgeNodeRegister(brp%numberOfTrailingEdges),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: setBoundIndicators out of memory"

  j = 0
  do i=1,brp%numberOfBoundaryNodes
   if(brp%isNodeTrailingEdge(i)) then
    j = j + 1
    brp%trailingEdgeNodeRegister(j) = i
   end if
  end do

  numberOfWallNodes = 0
  numberOfSymmetryNodes = 0
  numberOfFarfieldNodes = 0
  numberOfEngineInflowNodes = 0
  numberOfEngineOutflowNodes = 0 
  numberOfInviscidWallNodes  = 0
  numberOfInternalOutflowNodes = 0 
  do i=1,brp%numberOfBoundaryNodes
   ind = brp%boundaryNodes(i)%indicator
   if(ind==1) then
    numberOfWallNodes = numberOfWallNodes + 1
   else if(ind==2) then
    numberOfSymmetryNodes = numberOfSymmetryNodes + 1
   else if(ind==3) then
    numberOfFarfieldNodes(1) = numberOfFarfieldNodes(1) + 1
   else if(ind==4) then
    numberOfFarfieldNodes(2) = numberOfFarfieldNodes(2) + 1
   else if(ind==5) then 
    numberOfEngineInflowNodes(1) = numberOfEngineInflowNodes(1) + 1
   else if(ind==6) then 
    numberOfEngineInflowNodes(2) = numberOfEngineInflowNodes(2) + 1
   else if(ind==7) then 
    numberOfEngineOutflowNodes(1) = numberOfEngineOutflowNodes(1) + 1
   else if(ind==8) then 
    numberOfEngineOutflowNodes(2) = numberOfEngineOutflowNodes(2) + 1
   else if(ind==9) then 
   numberOfInviscidWallNodes = numberOfInviscidWallNodes + 1
   end if
  end do

  write(*,*) "wall nodes: ",numberOfWallNodes
  write(*,*) "symmetry nodes: ",numberOfSymmetryNodes
  write(*,*) "farfield nodes: ",numberOfFarfieldNodes

  brp%boundaryNodesRegister(-20:-1) = 0 
  brp%boundaryNodesRegister(-20) = 1
  brp%boundaryNodesRegister(-19) = numberOfWallNodes

  shiftNumber = numberOfWallNodes + 1

  brp%boundaryNodesRegister(-18) = shiftNumber
  brp%boundaryNodesRegister(-17) = shiftNumber + numberOfSymmetryNodes - 1

  shiftNumber = shiftNumber + numberOfSymmetryNodes

  brp%boundaryNodesRegister(-16) = shiftNumber
  brp%boundaryNodesRegister(-15) = shiftNumber + numberOfFarfieldNodes(1) - 1

  shiftNumber = shiftNumber + numberOfFarfieldNodes(1)

  brp%boundaryNodesRegister(-14) = shiftNumber
  brp%boundaryNodesRegister(-13) = shiftNumber + numberOfFarfieldNodes(2) - 1

  shiftNumber = shiftNumber + numberOfFarfieldNodes(2)

  brp%boundaryNodesRegister(-12) = shiftNumber
  brp%boundaryNodesRegister(-11) = shiftNumber + numberOfEngineInflowNodes(1) - 1

  shiftNumber = shiftNumber + numberOfEngineInflowNodes(1)

  brp%boundaryNodesRegister(-10) = shiftNumber
  brp%boundaryNodesRegister(-9) = shiftNumber + numberOfEngineInflowNodes(2) - 1

  shiftNumber = shiftNumber + numberOfEngineInflowNodes(2)

  brp%boundaryNodesRegister(-8) = shiftNumber
  brp%boundaryNodesRegister(-7) = shiftNumber + numberOfEngineOutflowNodes(1) - 1

  shiftNumber = shiftNumber + numberOfEngineOutflowNodes(1)

  brp%boundaryNodesRegister(-6) = shiftNumber
  brp%boundaryNodesRegister(-5) = shiftNumber + numberOfEngineOutflowNodes(2) - 1

  shiftNumber = shiftNumber + numberOfEngineOutflowNodes(2)

  brp%boundaryNodesRegister(-4) = shiftNumber
  brp%boundaryNodesRegister(-3) = shiftNumber + numberOfInviscidWallNodes - 1

  brp%numberOfTripleNodes = 0
  brp%numberOfDoubleNodes = 0

 end subroutine setBoundIndForViscousMesh
!-------------------------------------------------------------------------
 subroutine writeBoundaryData(OUTFILE,brp,isInitial)
 IMPLICIT NONE

 integer :: OUTFILE
 type(BoundaryRegisterData) :: brp
 logical :: isInitial
 
 integer :: i,j,currentIndicator,count(9)
 type(BoundaryNodeData),pointer :: currentNode

 write(*,*) "Writing boundary data..."

! write indicators
 do i=-20,0
  write(OUTFILE)  brp%boundaryNodesRegister(i)
 end do

 count = 0
 do currentIndicator=1,9
  do i=1,brp%numberOfBoundaryNodes
   if(brp%boundaryNodes(i)%indicator==currentIndicator) then 
    write(OUTFILE) i 
    count(currentIndicator) = count(currentIndicator) + 1
   else if(brp%boundaryNodes(i)%indicator>9) then 
    write(*,*) "unknown indicator: ",i,brp%boundaryNodes(i)%indicator
   end if
  end do
 end do 


 write(*,*) "Indicators: ",count

 ! write boundary face coefficients

 write(OUTFILE) brp%numberOfBoundarySides
 write(OUTFILE) ((brp%srp(i)%sideIndexes(j),i=1,brp%numberOfBoundarySides),j=1,2)
 write(OUTFILE) (brp%srp(i)%indicator,i=1,brp%numberOfBoundarySides)
 write(OUTFILE) ((brp%srp(i)%sideCoefficients(j),i=1,brp%numberOfBoundarySides),j=1,3)

 write(OUTFILE) brp%numberOfEngineInletSides
 write(OUTFILE) ((brp%engineInletSideIndexes(i,j),i=1,brp%numberOfEngineInletSides),j=1,2)
 write(OUTFILE) ((brp%engineInletSideCoefficients(i,j),i=1,brp%numberOfEngineInletSides),j=1,3)


 end subroutine writeBoundaryData
!-------------------------------------------------------------------------
 subroutine makeNodeNormals(brp,crp)
 IMPLICIT NONE

 type(BoundaryRegisterData) :: brp
 type(CoordinateRegisterData) :: crp

 integer :: i,j,i1,i2,i3,i4,allocateStatus,irot
 real :: flag(brp%numberOfBoundaryNodes)
 real :: mag,SP1,SP2,SP3,SP4,comp1(3),comp2(3),comp3(3),comp4(3)
 real :: N(3),NN,r1(3),r2(3),xrot(3),wrot(3),t1(3),t2(3),t3(3),t4(3)
 real :: rb1(3),rb2(3),rb3(3),rb4(3),wnorm(3),d1(3),d2(3),d3(3),d4(3)
 real :: rc1(3),rc2(3),rc3(3)
 real :: nwall1(6),nwall2(6),nwall3(6),nwall4(6),norm(6)
 real :: test1(6),test2(6),test3(6),test4(6),magtan1,magtan2,magtan3,magtan4
 real, parameter :: pi = 3.1416
! INITIALISE
 do i=1,brp%numberOfBoundaryNodes
  brp%boundaryNodes(i)%normal = 0.0
  brp%boundaryNodes(i)%wallNormal = 0.0
  flag(i)=0
 end do

 do j=1,brp%numberOfTriangularBFs
  i1 = brp%triangularFaces(j)%faceIndexes(1)
  i2 = brp%triangularFaces(j)%faceIndexes(2)
  i3 = brp%triangularFaces(j)%faceIndexes(3)

  r1 = getCoor(i2,crp)-getCoor(i1,crp)
  r2 = getCoor(i3,crp)-getCoor(i1,crp)
  N(1) = r1(2)*r2(3)-r1(3)*r2(2)
  N(2) = r1(3)*r2(1)-r1(1)*r2(3)
  N(3) = r1(1)*r2(2)-r1(2)*r2(1)

  NN = sqrt(sum(N*N))
  n = N/NN
  do i=1,3
   norm(i)=n(i)
   norm(i+3)=0.0
  enddo
! IF ROTATING SURFACE
  if(brp%surfaceSegmentRegister(brp%triangularFaces(j)%surfaceSegmentNumber)%wallboundindicator>0)then
  flag(i1)=flag(i1)+1
  flag(i2)=flag(i2)+1
  flag(i3)=flag(i3)+1
! DETERMINE ROTATION OBJECT NUMBER
   irot=brp%surfaceSegmentRegister(brp%triangularFaces(j)%surfaceSegmentNumber)%wallboundindicator
! GET THE ROTATION COORDINATE AND ROTATION VECTOR
   do i=1,3
    xrot(i)=brp%rotwall(irot,i)
    wrot(i)=brp%rotwall(irot,i+3)
   end do
! GET THE COORDINATES OF THE BOUNDARY NODES
   rb1 = getCoor(i1,crp)
   rb2 = getCoor(i2,crp)
   rb3 = getCoor(i3,crp)
! GET THE VECTOR FROM WHEEL CENTRE TO BOUNDARY
   rc1=rb1-xrot
   rc2=rb2-xrot
   rc3=rb3-xrot
! NORMALISE THE ROTATION VECTOR
   mag=SQRT(wrot(1)**2+wrot(2)**2+wrot(3)**2)
   wnorm=wrot/mag
! CALCULATE THE DISTANCE OF EACH BOUNDARY POINT FROM THE AXIS OF ROTATION
   SP1=wnorm(1)*rc1(1)+wnorm(2)*rc1(2)+wnorm(3)*rc1(3)
   SP2=wnorm(1)*rc2(1)+wnorm(2)*rc2(2)+wnorm(3)*rc2(3)
   SP3=wnorm(1)*rc3(1)+wnorm(2)*rc3(2)+wnorm(3)*rc3(3)
   comp1=SP1*wnorm
   comp2=SP2*wnorm
   comp3=SP3*wnorm
   d1=rc1-comp1
   d2=rc2-comp2
   d3=rc3-comp3
! CALCULATE THE TANGENT VECTORS (cross product of wrot and d)
   t1(1)=wrot(2)*d1(3)-wrot(3)*d1(2)
   t1(2)=wrot(3)*d1(1)-wrot(1)*d1(3)
   t1(3)=wrot(1)*d1(2)-wrot(2)*d1(1)
   t2(1)=wrot(2)*d2(3)-wrot(3)*d2(2)
   t2(2)=wrot(3)*d2(1)-wrot(1)*d2(3)
   t2(3)=wrot(1)*d2(2)-wrot(2)*d2(1)
   t3(1)=wrot(2)*d3(3)-wrot(3)*d3(2)
   t3(2)=wrot(3)*d3(1)-wrot(1)*d3(3)
   t3(3)=wrot(1)*d3(2)-wrot(2)*d3(1)
! CREATE THE MOVING WALL NORMAL + TANGENT VECTORS
   test1=brp%boundaryNodes(i1)%wallnormal
   magtan1=SQRT(test1(4)**2+test1(5)**2+test1(6)**2)
   test2=brp%boundaryNodes(i2)%wallnormal
   magtan2=SQRT(test2(4)**2+test2(5)**2+test2(6)**2)
   test3=brp%boundaryNodes(i3)%wallNormal
   magtan3=SQRT(test3(4)**2+test3(5)**2+test3(6)**2)
   
   do i=1,3
    nwall1(i)=n(i)
    nwall2(i)=n(i)
    nwall3(i)=n(i)
    nwall1(i+3)=t1(i)/xrot(3) !scale tangent vector based on wheel radius
    nwall2(i+3)=t2(i)/xrot(3)
    nwall3(i+3)=t3(i)/xrot(3)
   end do
    print*,'Wheel tangent:',nwall1(4),nwall1(5),nwall1(6)
! IF NO SLIP WALL 
  elseif(brp%surfaceSegmentRegister(brp%triangularFaces(j)%surfaceSegmentNumber)%wallboundindicator==0)then   
  flag(i1)=flag(i1)+1
  flag(i2)=flag(i2)+1
  flag(i3)=flag(i3)+1
   do i=1,3
    nwall1(i)=n(i)
    nwall2(i)=n(i)
    nwall3(i)=n(i)
    nwall1(i+3)=0.0
    nwall2(i+3)=0.0
    nwall3(i+3)=0.0
   end do
! IF ROLLING GROUND
  elseif(brp%surfaceSegmentRegister(brp%triangularFaces(j)%surfaceSegmentNumber)%wallboundindicator==-1)then

  flag(i1)=flag(i1)+1
  flag(i2)=flag(i2)+1
  flag(i3)=flag(i3)+1
   test1=brp%boundaryNodes(i1)%wallnormal
   magtan1=SQRT(test1(4)**2+test1(5)**2+test1(6)**2)
   test2=brp%boundaryNodes(i2)%wallnormal
   magtan2=SQRT(test2(4)**2+test2(5)**2+test2(6)**2)
   test3=brp%boundaryNodes(i3)%wallNormal
   magtan3=SQRT(test3(4)**2+test3(5)**2+test3(6)**2)
   do i=1,3
    nwall1(i)=n(i)
    nwall2(i)=n(i)
    nwall3(i)=n(i)
   end do
   nwall1(4)=cos((brp%angle/360.)*2.0*pi)
   nwall1(5)=sin((brp%angle/360.)*2.0*pi)
   nwall1(6)=0.0
   nwall2(4)=cos((brp%angle/360.)*2.0*pi)
   nwall2(5)=sin((brp%angle/360.)*2.0*pi)
   nwall2(6)=0.0
   nwall3(4)=cos((brp%angle/360.)*2.0*pi)
   nwall3(5)=sin((brp%angle/360.)*2.0*pi)
   nwall3(6)=0.0
  endif

  if(brp%surfaceSegmentRegister(brp%triangularFaces(j)%surfaceSegmentNumber)%indicator==1) then 
   brp%boundaryNodes(i1)%normal = brp%boundaryNodes(i1)%normal + nwall1
   brp%boundaryNodes(i2)%normal = brp%boundaryNodes(i2)%normal + nwall2
   brp%boundaryNodes(i3)%normal = brp%boundaryNodes(i3)%normal + nwall3
   brp%boundaryNodes(i1)%wallNormal = brp%boundaryNodes(i1)%wallNormal + nwall1
   brp%boundaryNodes(i2)%wallNormal = brp%boundaryNodes(i2)%wallNormal + nwall2
   brp%boundaryNodes(i3)%wallNormal = brp%boundaryNodes(i3)%wallNormal + nwall3
  else
   brp%boundaryNodes(i1)%normal = brp%boundaryNodes(i1)%normal + norm
   brp%boundaryNodes(i2)%normal = brp%boundaryNodes(i2)%normal + norm
   brp%boundaryNodes(i3)%normal = brp%boundaryNodes(i3)%normal + norm
  end if
 end do

 do j=1,brp%numberOfRectangularBFs
  i1 = brp%rectangularFaces(j)%faceIndexes(1)
  i2 = brp%rectangularFaces(j)%faceIndexes(2)
  i3 = brp%rectangularFaces(j)%faceIndexes(3)
  i4 = brp%rectangularFaces(j)%faceIndexes(4)

  r1 = getCoor(i2,crp)-getCoor(i1,crp)
  r2 = getCoor(i4,crp)-getCoor(i1,crp)
  N(1) = r1(2)*r2(3)-r1(3)*r2(2)
  N(2) = r1(3)*r2(1)-r1(1)*r2(3)
  N(3) = r1(1)*r2(2)-r1(2)*r2(1)

  r1 = getCoor(i3,crp)-getCoor(i2,crp)
  r2 = getCoor(i1,crp)-getCoor(i2,crp)
  N(1) = N(1) + r1(2)*r2(3)-r1(3)*r2(2)
  N(2) = N(2) + r1(3)*r2(1)-r1(1)*r2(3)
  N(3) = N(3) + r1(1)*r2(2)-r1(2)*r2(1)

  r1 = getCoor(i4,crp)-getCoor(i3,crp)
  r2 = getCoor(i2,crp)-getCoor(i3,crp)
  N(1) = N(1) + r1(2)*r2(3)-r1(3)*r2(2)
  N(2) = N(2) + r1(3)*r2(1)-r1(1)*r2(3)
  N(3) = N(3) + r1(1)*r2(2)-r1(2)*r2(1)

  r1 = getCoor(i1,crp)-getCoor(i4,crp)
  r2 = getCoor(i3,crp)-getCoor(i4,crp)
  N(1) = N(1) + r1(2)*r2(3)-r1(3)*r2(2)
  N(2) = N(2) + r1(3)*r2(1)-r1(1)*r2(3)
  N(3) = N(3) + r1(1)*r2(2)-r1(2)*r2(1)

  NN = sqrt(sum(N*N))

  n = N/NN
  do i=1,3
   norm(i)=n(i)
   norm(i+3)=0.0
  enddo

! IF ROTATING SURFACE
  if(brp%surfaceSegmentRegister(brp%rectangularFaces(j)%surfaceSegmentNumber)%wallboundindicator>0)then
  flag(i1)=flag(i1)+1
  flag(i2)=flag(i2)+1
  flag(i3)=flag(i3)+1
  flag(i4)=flag(i4)+1
! DETERMINE ROTATION OBJECT NUMBER
   irot=brp%surfaceSegmentRegister(brp%rectangularFaces(j)%surfaceSegmentNumber)%wallboundindicator
! GET THE ROTATION COORDINATE AND ROTATION VECTOR
   do i=1,3
    xrot(i)=brp%rotwall(irot,i)
    wrot(i)=brp%rotwall(irot,i+3)
   end do
! GET THE COORDINATES OF THE BOUNDARY NODES
   rb1 = getCoor(i1,crp)
   rb2 = getCoor(i2,crp)
   rb3 = getCoor(i3,crp)
   rb4 = getCoor(i4,crp)
! NORMALISE THE ROTATION VECTOR
   mag=SQRT(wrot(1)**2+wrot(2)**2+wrot(3)**2)
   wnorm=wrot/mag
! CALCULATE THE DISTANCE OF EACH BOUNDARY POINT FROM THE AXIS OF ROTATION
   SP1=wnorm(1)*rb1(1)+wnorm(2)*rb1(2)+wnorm(3)*rb1(3)
   SP2=wnorm(1)*rb2(1)+wnorm(2)*rb2(2)+wnorm(3)*rb2(3)
   SP3=wnorm(1)*rb3(1)+wnorm(2)*rb3(2)+wnorm(3)*rb3(3)
   SP4=wnorm(1)*rb4(1)+wnorm(2)*rb4(2)+wnorm(3)*rb4(3)
   comp1=SP1*wnorm
   comp2=SP2*wnorm
   comp3=SP3*wnorm
   comp4=SP4*wnorm
   d1=rb1-comp1
   d2=rb2-comp2
   d3=rb3-comp3
   d4=rb4-comp4
! CALCULATE THE TANGENT VECTORS (cross product of wrot and d)
   t1(1)=wrot(2)*d1(3)-wrot(3)*d1(2)
   t1(2)=wrot(3)*d1(1)-wrot(1)*d1(3)
   t1(3)=wrot(1)*d1(2)-wrot(2)*d1(1)
   t2(1)=wrot(2)*d2(3)-wrot(3)*d2(2)
   t2(2)=wrot(3)*d2(1)-wrot(1)*d2(3)
   t2(3)=wrot(1)*d2(2)-wrot(2)*d2(1)
   t3(1)=wrot(2)*d3(3)-wrot(3)*d3(2)
   t3(2)=wrot(3)*d3(1)-wrot(1)*d3(3)
   t3(3)=wrot(1)*d3(2)-wrot(2)*d3(1)
   t4(1)=wrot(2)*d4(3)-wrot(3)*d4(2)
   t4(2)=wrot(3)*d4(1)-wrot(1)*d4(3)
   t4(3)=wrot(1)*d4(2)-wrot(2)*d4(1)
! CREATE THE MOVING WALL NORMAL + TANGENT VECTORS
   test1=brp%boundaryNodes(i1)%wallnormal
   magtan1=SQRT(test1(4)**2+test1(5)**2+test1(6)**2)
   test2=brp%boundaryNodes(i2)%wallnormal
   magtan2=SQRT(test2(4)**2+test2(5)**2+test2(6)**2)
   test3=brp%boundaryNodes(i3)%wallNormal
   magtan3=SQRT(test3(4)**2+test3(5)**2+test3(6)**2)
   test4=brp%boundaryNodes(i4)%wallNormal
   magtan4=SQRT(test4(4)**2+test4(5)**2+test4(6)**2)
   
   do i=1,3
    nwall1(i)=n(i)
    nwall2(i)=n(i)
    nwall3(i)=n(i)
    nwall4(i)=n(i)
    nwall1(i+3)=t1(i)/xrot(3)  !scale tangent based on wheel radius
    nwall2(i+3)=t2(i)/xrot(3)
    nwall3(i+3)=t3(i)/xrot(3)
    nwall4(i+3)=t4(i)/xrot(3)
   end do
! IF NO SLIP WALL 
  elseif(brp%surfaceSegmentRegister(brp%rectangularFaces(j)%surfaceSegmentNumber)%wallboundindicator==0)then   
  flag(i1)=flag(i1)+1
  flag(i2)=flag(i2)+1
  flag(i3)=flag(i3)+1
  flag(i4)=flag(i4)+1
   do i=1,3
    nwall1(i)=n(i)
    nwall2(i)=n(i)
    nwall3(i)=n(i)
    nwall4(i)=n(i)
    nwall1(i+3)=0.0
    nwall2(i+3)=0.0
    nwall3(i+3)=0.0
    nwall4(i+3)=0.0
   end do
! IF ROLLING GROUND
  elseif(brp%surfaceSegmentRegister(brp%rectangularFaces(j)%surfaceSegmentNumber)%wallboundindicator==-1)then  
  flag(i1)=flag(i1)+1
  flag(i2)=flag(i2)+1
  flag(i3)=flag(i3)+1
  flag(i4)=flag(i4)+1
   test1=brp%boundaryNodes(i1)%wallnormal
   magtan1=SQRT(test1(4)**2+test1(5)**2+test1(6)**2)
   test2=brp%boundaryNodes(i2)%wallnormal
   magtan2=SQRT(test2(4)**2+test2(5)**2+test2(6)**2)
   test3=brp%boundaryNodes(i3)%wallNormal
   magtan3=SQRT(test3(4)**2+test3(5)**2+test3(6)**2)
   test4=brp%boundaryNodes(i4)%wallNormal
   magtan4=SQRT(test4(4)**2+test4(5)**2+test4(6)**2)
   do i=1,3
    nwall1(i)=n(i)
    nwall2(i)=n(i)
    nwall3(i)=n(i)
    nwall4(i)=n(i)
   end do
   nwall1(4)=1.0
   nwall1(5)=0.0
   nwall1(6)=0.0
   nwall2(4)=1.0
   nwall2(5)=0.0
   nwall2(6)=0.0
   nwall3(4)=1.0
   nwall3(5)=0.0
   nwall3(6)=0.0
   nwall4(4)=1.0
   nwall4(5)=0.0
   nwall4(6)=0.0
  endif

  if(brp%surfaceSegmentRegister(brp%rectangularFaces(j)%surfaceSegmentNumber)%indicator==1) then
   brp%boundaryNodes(i1)%normal = brp%boundaryNodes(i1)%normal + nwall1
   brp%boundaryNodes(i2)%normal = brp%boundaryNodes(i2)%normal + nwall2
   brp%boundaryNodes(i3)%normal = brp%boundaryNodes(i3)%normal + nwall3
   brp%boundaryNodes(i4)%normal = brp%boundaryNodes(i4)%normal + nwall4
   brp%boundaryNodes(i1)%wallNormal = brp%boundaryNodes(i1)%wallNormal + nwall1
   brp%boundaryNodes(i2)%wallNormal = brp%boundaryNodes(i2)%wallNormal + nwall2
   brp%boundaryNodes(i3)%wallNormal = brp%boundaryNodes(i3)%wallNormal + nwall3
   brp%boundaryNodes(i4)%wallNormal = brp%boundaryNodes(i4)%wallNormal + nwall4
  else
   brp%boundaryNodes(i1)%normal = brp%boundaryNodes(i1)%normal + norm
   brp%boundaryNodes(i2)%normal = brp%boundaryNodes(i2)%normal + norm
   brp%boundaryNodes(i3)%normal = brp%boundaryNodes(i3)%normal + norm
   brp%boundaryNodes(i4)%normal = brp%boundaryNodes(i4)%normal + norm
  end if
 end do


 do j=1,brp%numberOfBoundaryNodes
  N(1:3) = brp%boundaryNodes(j)%normal(1:3)
  NN = sqrt(sum(N*N))
  if(NN>1.0e-20) then 
   brp%boundaryNodes(j)%normal(1) = brp%boundaryNodes(j)%normal(1)/NN
   brp%boundaryNodes(j)%normal(2) = brp%boundaryNodes(j)%normal(2)/NN
   brp%boundaryNodes(j)%normal(3) = brp%boundaryNodes(j)%normal(3)/NN
  else
   write(*,*) "WARNING: Small node: ",j,getCoor(j,crp),N
   brp%boundaryNodes(j)%normal = 0.0
  end if 

  N(1:3) = brp%boundaryNodes(j)%wallNormal(1:3)
  NN = sqrt(sum(N*N))
  if(NN>1.0e-20) then
   brp%boundaryNodes(j)%wallNormal(1) = brp%boundaryNodes(j)%wallNormal(1)/NN
   brp%boundaryNodes(j)%wallNormal(2) = brp%boundaryNodes(j)%wallNormal(2)/NN
   brp%boundaryNodes(j)%wallNormal(3) = brp%boundaryNodes(j)%wallNormal(3)/NN
  end if
   if(flag(j).gt.0.5)then
   brp%boundaryNodes(j)%wallNormal(4) = brp%boundaryNodes(j)%wallNormal(4)/flag(j)
   brp%boundaryNodes(j)%wallNormal(5) = brp%boundaryNodes(j)%wallNormal(5)/flag(j)
   brp%boundaryNodes(j)%wallNormal(6) = brp%boundaryNodes(j)%wallNormal(6)/flag(j)
   brp%boundaryNodes(j)%Normal(4) = brp%boundaryNodes(j)%Normal(4)/flag(j)
   brp%boundaryNodes(j)%Normal(5) = brp%boundaryNodes(j)%Normal(5)/flag(j)
   brp%boundaryNodes(j)%Normal(6) = brp%boundaryNodes(j)%Normal(6)/flag(j)
   endif
  
   write(175,*) j,getCoor(j,crp),brp%boundaryNodes(j)%Normal(4:6)
   write(175,*) j,getCoor(j,crp),brp%boundaryNodes(j)%wallNormal(4:6)
   
 end do

 end subroutine makeNodeNormals
!-------------------------------------------------------------------------
end module BoundaryRegister

