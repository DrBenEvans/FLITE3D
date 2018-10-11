!**********************************************
!*Aggl3d_Element.f90                          *
!*                                            *
!* Element module definition file             *
!*                                            * 
!*   gen3d  v 2.3.2                           * 
!*                                            *
!*                                            *
!*24.06.99-                                   *
!**********************************************
!*Description:                                *
!* File contains element datastructure        *
!**********************************************

module HexahedralElements
 integer,parameter,private:: NSD = 3

 type :: SHEData ! data structure for a single element
  integer :: pointIndexes(8) ! indexes to points in element
!  real :: area 
 end type SHEData

 type HexahedralElementData ! data structure for element array
   type(SHEData), pointer :: elements(:) 
   integer ::  numberOfElements 
 end type HexahedralElementData

 contains

!-------------------------------------------------------------------------
 subroutine cleanUpHexElData(hep)
 ! deallocates memory used by element data
 
 IMPLICIT NONE

 type(HexahedralElementData),pointer :: hep
 
 integer :: allocateStatus

  deallocate(hep%elements,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: cleanUpHexElData couldn't deallocate"

  deallocate(hep,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: cleanUpHexElData couldn't deallocate"
  
  nullify(hep)
 end subroutine cleanUpHexElData
!-------------------------------------------------------------------------
 subroutine constructHexElArray(n,hep)
! Sets up a REData type for use

  IMPLICIT NONE 

  integer :: n
  type(HexahedralElementData) :: hep

  integer :: allocateStatus

! allocate memory for element array
 
  allocate(hep%elements(n),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create element array"

  hep%numberOfElements=n  
 
 end subroutine constructHexElArray
!-------------------------------------------------------------------------
 subroutine readHexElData(INFILE,hep)
! Reads element data 
 IMPLICIT NONE

 integer :: INFILE
 type(HexahedralElementData) :: hep

 integer :: i,ip
 integer :: nel 
 logical :: found

  nel = hep%numberOfElements

  read(INFILE) ((hep%elements(ip)%pointIndexes(i),ip=1,nel),i=1,8)

 end subroutine readHexElData
!-------------------------------------------------------------------------
  function getHexElement(i,hep) result(c)
  
  integer :: i
  type(HexahedralElementData) :: hep
  type(SHEData),pointer :: c

  c => hep%elements(i)

 end function getHexElement
!-------------------------------------------------------------------------
end module HexahedralElements

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

module PrismElements
 integer,parameter,private:: NSD = 3

 type :: SPRData ! data structure for a single element    
  integer :: pointIndexes(6) ! indexes to points in element 
!  real :: area 
 end type SPRData

 type PrismElementData ! data structure for element array
   type(SPRData), pointer :: elements(:) 
   integer ::  numberOfElements
 end type PrismElementData

 contains

!-------------------------------------------------------------------------
 subroutine cleanUpPriElData(prp)
 ! deallocates memory used by element data

 IMPLICIT NONE

 type(PrismElementData),pointer :: prp
 
 integer :: allocateStatus

  deallocate(prp%elements,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: cleanUpPriElData couldn't deallocate"

  deallocate(prp,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: cleanUpPriElData couldn't deallocate"

  nullify(prp)
 end subroutine cleanUpPriElData
!-------------------------------------------------------------------------
 subroutine constructPriElArray(n,prp)
! Sets up a REData type for use

  IMPLICIT NONE

  integer :: n
  type(PrismElementData) :: prp

  integer :: allocateStatus

! allocate memory for element array

  allocate(prp%elements(n),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create element array"

  prp%numberOfElements=n

 end subroutine constructPriElArray
!-------------------------------------------------------------------------
 subroutine readPriElData(INFILE,prp)
! Reads element data
 IMPLICIT NONE

 integer :: INFILE
 type(PrismElementData) :: prp

 integer :: i,ip
 integer :: nel
 logical :: found

  nel = prp%numberOfElements

  read(INFILE) ((prp%elements(ip)%pointIndexes(i),ip=1,nel),i=1,6)

 end subroutine readPriElData
!-------------------------------------------------------------------------
  function getPriElement(i,prp) result(c)

  integer :: i
  type(PrismElementData) :: prp
  type(SPRData),pointer :: c

  c => prp%elements(i)

 end function getPriElement
!-------------------------------------------------------------------------
end module PrismElements

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

module PyramidElements
 integer,parameter,private:: NSD = 3

 type :: SPYData ! data structure for a single element
  integer :: pointIndexes(5) ! indexes to points in element
!  real :: area
 end type SPYData

 type PyramidElementData ! data structure for element array
   type(SPYData), pointer :: elements(:)
   integer ::  numberOfElements
 end type PyramidElementData

 contains

!-------------------------------------------------------------------------
 subroutine cleanUpPyrElData(pyp)
 ! deallocates memory used by element data

 IMPLICIT NONE

 type(PyramidElementData),pointer :: pyp

 integer :: allocateStatus

  deallocate(pyp%elements,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: cleanUpPyrElData couldn't deallocate"

  deallocate(pyp,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: cleanUpPyrElData couldn't deallocate"

  nullify(pyp)

 end subroutine cleanUpPyrElData
!-------------------------------------------------------------------------
 subroutine constructPyrElArray(n,pyp)
! Sets up a REData type for use

  IMPLICIT NONE

  integer :: n
  type(PyramidElementData) :: pyp

  integer :: allocateStatus

! allocate memory for element array

  allocate(pyp%elements(n),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create element array"

  pyp%numberOfElements=n

 end subroutine constructPyrElArray
!-------------------------------------------------------------------------
 subroutine readPyrElData(INFILE,pyp)
! Reads element data
 IMPLICIT NONE

 integer :: INFILE
 type(PyramidElementData) :: pyp
 logical :: found

 integer :: i,ip
 integer :: nel

  nel = pyp%numberOfElements

  read(INFILE) ((pyp%elements(ip)%pointIndexes(i),ip=1,nel),i=1,5)

 end subroutine readPyrElData
!-------------------------------------------------------------------------
  function getPyrElement(i,pyp) result(c)

  integer :: i
  type(PyramidElementData) :: pyp
  type(SPYData),pointer :: c

  c => pyp%elements(i)

 end function getPyrElement
!-------------------------------------------------------------------------
end module PyramidElements

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

module TetrahedralElements
 integer,parameter,private:: NSD = 3

 type :: STEData ! data structure for a single tetrhedral element
  integer :: pointIndexes(4) ! indexes to points in element
 end type STEData

 type TetrahedralElementData ! data structure for element array
   type(STEData), pointer :: elements(:)
   integer ::  numberOfElements
 end type TetrahedralElementData

 contains

!-------------------------------------------------------------------------
 subroutine cleanUpTetElData(tep)
 ! deallocates memory used by element data

 IMPLICIT NONE

 type(TetrahedralElementData),pointer :: tep

 integer :: allocateStatus

  deallocate(tep%elements,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: cleanUpTetElData couldn't deallocate"

  deallocate(tep,stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: cleanUpTetElData couldn't deallocate"

  nullify(tep)
 end subroutine cleanUpTetElData
!-------------------------------------------------------------------------
 subroutine constructTetElArray(n,tep)
! Sets up a TEData type for use

  IMPLICIT NONE

  integer :: n
  type(TetrahedralElementData) :: tep

  integer :: allocateStatus

! allocate memory for element array

  allocate(tep%elements(n),stat=allocateStatus)
  if (allocateStatus /= 0) STOP "ERROR: Not enough memory to create element array"

  tep%numberOfElements=n

 end subroutine constructTetElArray
!-------------------------------------------------------------------------
 subroutine readTetElData(INFILE,tep)
! Reads element data
 IMPLICIT NONE

 integer :: INFILE
 type(TetrahedralElementData) :: tep

 integer :: i,ip
 integer :: nel,buff
 logical :: found
 nel = tep%numberOfElements

! do i=1,nel
!  read(INFILE) tep%elements(i)%pointIndexes(1:4)
! end do


 read(INFILE) ((tep%elements(ip)%pointIndexes(i),ip=1,nel),i=1,4)

 end subroutine readTetElData
!-------------------------------------------------------------------------
  function getTetElement(i,tep) result(c)

  integer :: i
  type(TetrahedralElementData) :: tep
  type(STEData),pointer :: c

  c => tep%elements(i)

 end function getTetElement
!-------------------------------------------------------------------------
end module TetrahedralElements


