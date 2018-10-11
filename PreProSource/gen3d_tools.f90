!**********************************************
!*Aggl3d_def.f90                              *
!*                                            *
!*Module definition file for program          *
!*                                            * 
!*   Aggl3d v 2.3.1                           * 
!*                                            *
!*24.06.99-                                   *
!**********************************************
!*Description:                                *
!* This file includes procedures for file     *
!* communication, construction of side based  *
!* datastructure from elements, coloring for  *
!* for vectorization and construction of      *
!* coarser grids from the given original one  *
!* by agglomeration - for use in multigrid    *
!**********************************************
module Toolbox
! various handy datastructures

 type LinkedInteger
  integer :: int  
  type(linkedInteger),pointer :: next
 end type linkedInteger

 type IntegerBlockInstance
  integer,pointer :: block(:,:)
  type(IntegerBlockInstance),pointer :: next 
 end type IntegerBlockInstance

 type LinkedIntegerBlockChain
 ! unordered block chain of integer sets
  type(IntegerBlockInstance),pointer :: first,last
  integer,pointer :: nodeIndexes(:)
  integer :: blockSize
  integer :: numberOfNodes
  integer :: numberOfBlocks
  integer :: blockIndex  
 end type LinkedIntegerBlockChain


 type LinkedDoubleInteger
   integer :: index(2)
   type(LinkedDoubleInteger),pointer :: next
 end type LinkedDoubleInteger

 type LinkedIntegerP
  type(LinkedInteger),pointer :: first
 end type LinkedIntegerP

 type LinkedDoubleIntegerP
  type(LinkedDoubleInteger),pointer :: first,last 
 end type LinkedDoubleIntegerP

 type LinkedIntegerArray
  integer,pointer :: arr(:)   
  type(LinkedIntegerArray),pointer :: next
 end type LinkedIntegerArray

 type LinkedIntegerPArray
  type(linkedIntegerPArray),pointer :: lipa(:)
 end type linkedIntegerPArray

 type LinkedIntegerArrayP
  type(LinkedIntegerArray),pointer :: first
 end type LinkedIntegerArrayP

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

contains
 
!-----------------------------------------------------------------------
 subroutine LIBCSetUp(numberOfNodes,blockSize,libc)
 IMPLICIT NONE

 integer :: numberOfNodes,blockSize
 type(LinkedIntegerBlockChain),pointer :: libc

 integer :: allocateStatus
 
  allocate(libc,stat=allocateStatus)
  if(allocateStatus>0) STOP "ERROR: Not enough memory in LIBCSetUp"

  libc%blockSize = blockSize
  libc%numberOfNodes = numberOfNodes
  libc%numberOfBlocks = 1
  libc%blockIndex = 0
 
  allocate(libc%nodeIndexes(numberOfNodes),stat=allocateStatus)
  if(allocateStatus>0) STOP "ERROR: Not enough memory in LIBCSetUp"
  libc%nodeIndexes = 0
  allocate(libc%first,stat=allocateStatus)
  if(allocateStatus>0) STOP "ERROR: Not enough memory in LIBCSetUp"
  allocate(libc%first%block(blockSize,2),stat=allocateStatus)
  if(allocateStatus>0) STOP "ERROR: Not enough memory in LIBCSetUp"
  libc%last=>libc%first
  libc%first%block = 0
  nullify(libc%last%next)
 end subroutine LIBCSetUp 
!-----------------------------------------------------------------------
 subroutine LIBCDestroy(libc)
 IMPLICIT NONE
 type(LinkedIntegerBlockChain),pointer :: libc

 integer :: allocateStatus
 
 type(IntegerBlockInstance),pointer :: current

  deallocate(libc%nodeIndexes,stat=allocateStatus)
  if(allocateStatus>0) STOP "ERROR: LIBCDestroy couldn't deallocate"

  current => libc%first
  do while(associated(current)) 
   libc%first => current
   current => current%next
   deallocate(libc%first%block,stat=allocateStatus)
   if(allocateStatus>0) STOP "ERROR: LIBCDestroy couldn't deallocate"
   deallocate(libc%first,stat=allocateStatus)
   if(allocateStatus>0) STOP "ERROR: LIBCDestroy couldn't deallocate"
  end do

  deallocate(libc,stat=allocateStatus)
  if(allocateStatus>0) STOP "ERROR: LIBCDestroy couldn't deallocate"
  nullify(libc)
 end subroutine LIBCDestroy
!-----------------------------------------------------------------------
 subroutine LIBCAddInstance(int,node,libc)
 IMPLICIT NONE

 integer :: int,node
 type(LinkedIntegerBlockChain) :: libc

 integer :: allocateStatus

  if(libc%blockIndex==libc%blockSize) then
   ! do some allocation
   allocate(libc%last%next,stat=allocateStatus)
   if(allocateStatus>0) STOP "ERROR: Not enough memory in LIBCAddInstance"
   libc%last=>libc%last%next
   allocate(libc%last%block(libc%blockSize,2),stat=allocateStatus)
   if(allocateStatus>0) STOP "ERROR: Not enough memory in LIBCAddInstance"
   nullify(libc%last%next)
   libc%blockIndex=1
   libc%numberOfBlocks = libc%numberOfBlocks + 1
  else
   libc%blockIndex = libc%blockIndex + 1
  end if

  libc%last%block(libc%blockIndex,1) = int
  libc%last%block(libc%blockIndex,2) = libc%nodeIndexes(node)
  libc%nodeIndexes(node) = libc%blockSize*(libc%numberOfBlocks-1) + libc%blockIndex
 end subroutine LIBCAddInstance
!-----------------------------------------------------------------------
 subroutine getInstance(sb,ind,libc)
 IMPLICIT NONE

 integer :: ind,sb(2)
 type(LinkedIntegerBlockChain) :: libc

 integer :: block,index,i
 type(IntegerBlockInstance),pointer :: current

  block = (ind-1)/libc%blockSize + 1
  index = ind - (block-1)*libc%blockSize

  current => libc%first
  do i=2,block
   current=>current%next
  end do
  sb = current%block(index,:)
 end subroutine getInstance
!-----------------------------------------------------------------------
 logical function LIBCFindInstance(int,node,libc)
 IMPLICIT NONE
 
 integer :: int,node
 type(LinkedIntegerBlockChain) :: libc

 integer :: searchIndex,se(2)
 
 searchIndex = libc%nodeIndexes(node)
 do while(searchIndex>0)
  call getInstance(se,searchIndex,libc)
  if(se(1)==int) exit
  searchIndex = se(2)
 end do 

 LIBCFindInstance = (searchIndex>0)
 end function LIBCFindInstance
!-----------------------------------------------------------------------
end module Toolbox

