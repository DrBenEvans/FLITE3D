!********************
!*mgnsg3d_main.f90  *
!*                  *
!*Main file for     *
!*  mgnsg3d v 2.5par*
!*                  *
!********************

      program main
      use MultigridSolver
      IMPLICIT NONE
      !#include <mpif.h>      

      character*100 solverCommand 
 
      integer, parameter :: realType = MPI_REAL8    
      integer, parameter :: intType = MPI_INTEGER4  
   

#ifdef PARALLEL
 
#endif PARALLEL

      integer :: idbg
      common/debug/idbg

      

      idbg = +1

#ifdef PARALLEL

      allocate(baseParallelData,stat=allocateStatus)
      if(allocateStatus/=0) STOP "ERROR: Not enough memory in main 1" 
      baseParallelData%integerType = intType 
      baseParallelData%realType = realType 

      call mp_init_message (baseParallelData%currentTID,&
                   baseParallelData%numberOfProcesses)

      if  (is_master(baseParallelData%currentTID)) then

       call software_version

       call executePar(.true.) 
      else
       call executePar(.false.) 
      end if
#else PARALLEL

      baseParallelData%numberOfProcesses = 1

      call executeSeq() 
#endif PARALLEL


#ifdef PARALLEL

      call mp_stop

#endif PARALLEL

      print *, ' TERMINATION : O.K.'

      stop
      end
