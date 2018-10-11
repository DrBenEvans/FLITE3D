module Comms
  include 'mpif.h' 
  integer    numprocjj, buffersizejj,SENDREQ,RECVREQ
  parameter  (numprocjj = 256)
  parameter  (buffersizejj = 1000000 )
  parameter  (SENDREQ = 1 )
  parameter  (RECVREQ = 2 )

  integer    myrank

  integer   buffersjj(buffersizejj,2,numprocjj+1)
  integer   bufferposjj(2,numprocjj+1)
  integer   requestsjj(2,numprocjj+1)

  integer statuses(MPI_STATUS_SIZE,2,numprocjj+1)     
  integer saved(2,numprocjj+1)
contains

!----------------------------------------------------------------------
      subroutine software_version

      implicit none

      print 10000

10000 format(1x, //,23x,'==================================', &
                 /,23x,'*                                *', &
                 /,23x,'*  MGNS3D (Version 1.0)          *', &
                 /,23x,'*  Unsteady Compressible Flows   *', &
                 /,23x,'*  ----------------------        *', &
                 /,23x,'*                                *', &
                 /,23x,'*                                *', &
                 /,23x,'*         Flite  CODE            *', &
                 /,23x,'*                                *', &
                 /,23x,'*  AUTHORS: K.S  &  O.H          *', &
                 /,23x,'*                                *', &
                 /,23x,'*  copyright (c) Univ. Swansea   *', &
                 /,23x,'*                                *', &
                 /,23x,'*                                *', &
                 /,23x,'==================================', &
                //)
      return
      end subroutine
!----------------------------------------------------------------------
      logical function is_master( mytid )

      implicit none

      include 'mpif.h'

      integer :: mytid,jj

      is_master = (mytid .eq. 0)

      return
      end function
!----------------------------------------------------------------------
      subroutine master(np,mytid,my_dm,tids,solverCommand,len)

      implicit none

      integer  mytid,np,my_dm,tids(*)

      integer  nslave, ipar_dbg,len
      character*200 solverCommand
      
      integer i

      print *, ' # running processes = ',  np,mytid

      print *,solverCommand(1:len),mytid

!ZZ
      do i=1,np
        tids(i) = i-1
      end do
      my_dm = 1
!ZZ

      return
      end subroutine
!----------------------------------------------------------------------
      subroutine slave_get_tids(np,my_dm,tids,mytid)

      implicit none

      include 'mpif.h'

      integer np,mytid,tids(*) ,my_dm

      integer param_tag
      parameter(param_tag=0)

      integer info, inum, i


!ZZ
      do i=1,np
        tids(i)=i-1
      end do
!ZZ

      my_dm=-1

      do i=2,np
        if( mytid .eq. tids(i) ) my_dm = i
      end do
      if(my_dm.eq.-1)stop 'wrong domain number in slave_get_tids'
         
      return
      end subroutine
!----------------------------------------------------------------------
      subroutine second(time)


      implicit none
      real time

      integer it_sys(3)
      external itime


      call itime(it_sys)


      time = real(it_sys(1))*3600+real(it_sys(2))*60+real(it_sys(3))

      return
      end subroutine
!----------------------------------------------------------------------
      subroutine mp_init_message (mytid,np)

      implicit none



      include 'mpif.h'

      integer mytid

      integer info,np,i,j


      call MPI_INIT ( info )
      call MPI_COMM_RANK(MPI_COMM_WORLD,mytid,info)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,np,info)

      do j = 1,2
        do i = 1,numprocjj+1
          requestsjj(j,i) = MPI_REQUEST_NULL
          bufferposjj(j,i) = 0
        enddo
      enddo

      myrank = mytid

      print *, 'STARTING ...',mytid,'of',np

!#ifndef CRAYT3D
!      call pvmfsetopt(PVMAUTOERR,PVMROUTEDIRECT,info)
!      call mp_test('mp_init_message failed', mytid)
!#endif CRAYT3D

      return
      end subroutine

!----------------------------------------------------------------------
      subroutine mp_stop

      implicit none

      include 'mpif.h'

      integer info


      print *, ' NOW ... EXIT MPI '

      call MPI_FINALIZE(info)

      call mp_test('mp_stop failed', info)

      return
      end subroutine

!----------------------------------------------------------------------
      subroutine mp_init_buffer(tids,dest,info)

      implicit none


      include 'mpif.h'

      integer encoding

      integer info,tids(*),dest,tid,i

!      external  pvmfinitsend, mp_test
!      call pvmfinitsend (encoding,info)
!      call mp_test('mp_init_buffer failed',info)

      if( dest.eq.-1 ) then
        tid = numprocjj + 1
      else
        tid = tids(dest) + 1
      endif
      do i = 1,buffersizejj
        buffersjj(i,SENDREQ,tid) = 0
        buffersjj(i,RECVREQ,tid) = 0
      enddo
      requestsjj(SENDREQ,tid) = MPI_REQUEST_NULL
      bufferposjj(SENDREQ,tid) = 0
      requestsjj(RECVREQ,tid) = MPI_REQUEST_NULL
      bufferposjj(RECVREQ,tid) = 0
      return
      end subroutine
!----------------------------------------------------------------------
      subroutine mp_recv(tids,from,tag,info)

      implicit none


      include 'mpif.h'


      integer from,tag, info,tids(*), a, tid
    
      if( from.eq.-1 )   call abort()
      tid = tids(from)+1
      requestsjj(RECVREQ,tid) = MPI_REQUEST_NULL
      bufferposjj(RECVREQ,tid) = 0
!     write(200+myrank,*) 'MPI_IRECV from ',tids(from),buffersizejj
      call MPI_IRECV( buffersjj(:,RECVREQ,tid), &
                      buffersizejj, MPI_PACKED, tids(from), tag, &
                      MPI_COMM_WORLD, requestsjj(RECVREQ,tid), info )

!       call MPI_IRECV( buffersjj(:,RECVREQ,tid),&
!                       buffersizejj, MPI_PACKED, tids(from), tag, &
!                       MPI_COMM_WORLD, requestsjj(RECVREQ,tid), info )
      call mp_test('mp_recv failed',info)

      return
      end subroutine

!----------------------------------------------------------------------
      subroutine mp_sendv(tids,dest,tag,info)

      implicit none


      include 'mpif.h'

      integer dest, tag, info,tids(*), tid

      integer STATUS(MPI_STATUS_SIZE)     

      tid = tids(dest)+1
      call MPI_ISEND(buffersjj(:,SENDREQ,tid),&
                     bufferposjj(SENDREQ,tid),&
                     MPI_PACKED,tids(dest),tag, &
                     MPI_COMM_WORLD,requestsjj(SENDREQ,tid),info)

      call mp_test('mp_sendv failed',info)

      return
      end subroutine

!----------------------------------------------------------------------
      integer function get_unit()
      implicit none

      get_unit = 200 + myrank
      return
      end function

!----------------------------------------------------------------------

      subroutine mp_wait_comms(info)
      implicit none
!
      include 'mpif.h'

      integer info, i, msgsize, msgsize2


!     write(200+myrank,*)  'mp_wait_for_comms started'
      do i = 1,numprocjj+1
!       write(200+myrank,*)  'Waiting for comms from ',i-1,requestsjj(SENDREQ,i),requestsjj(RECVREQ,i)
!        saved(SENDREQ,i) = requestsjj(SENDREQ,i)
!        saved(RECVREQ,i) = requestsjj(RECVREQ,i)
      enddo
      call MPI_WAITALL(2*(numprocjj+1), requestsjj, statuses, info )
      call mp_test('mp_wait_comms failed',info)

!      do i = 1,numprocjj+1
!        write(200+myrank,*)  'Done comms from ',i-1,statuses(MPI_ERROR,SENDREQ,i),&
!                              statuses(MPI_ERROR,RECVREQ,i)
!      enddo
!      do i = 1,numprocjj+1
!        if( saved(RECVREQ,i).eq.0) then
!          msgsize = -1111
!        else
!          call MPI_GET_COUNT( statuses( 1, RECVREQ, i ), MPI_PACKED, msgsize )
!        endif
!        write(200+myrank,*)  'Comms size from ',i-1,msgsize
!      enddo
!     write(200+myrank,*)  'mp_wait_for_comms finished'
      call mp_test('mp_wait_comms failed',info)

!      call MPI_BARRIER( MPI_COMM_WORLD, info )
!      write(200+myrank,*)  'mp_wait_for_comms barrier finished'

      return
      end subroutine

!----------------------------------------------------------------------
      subroutine mp_pak_real(tids,dest,type,scalar,info)

      implicit none


      include 'mpif.h'

      integer type,info,tids(*),dest,tid

      real   scalar

      if( dest.eq.-1 ) then
        tid = numprocjj + 1
      else
        tid = tids(dest) + 1
      endif

!     write(200+myrank,*)  'mp_pak_real ',scalar,' into pos ',bufferposjj(SENDREQ,tid)

      call MPI_PACK(scalar,1,type,buffersjj(:,SENDREQ,tid),&
                    buffersizejj,bufferposjj(SENDREQ,tid),&
                    MPI_COMM_WORLD,info)

      call mp_test('mp_pak_real failed',info)

      return
      end subroutine
!----------------------------------------------------------------------
      subroutine mp_pak_string(tids,dest,str,strnlen,info )

      implicit none

      include 'mpif.h'

      integer info,strnlen,tids(*),dest,tid
      character*80 str

      if( dest.eq.-1 ) then
        tid = numprocjj + 1
      else
        tid = tids(dest) + 1
      endif

!     write(200+myrank,*)  'mp_pak_string ',str,strnlen,' into pos ',bufferposjj(SENDREQ,tid)

      call MPI_PACK(str,strnlen,MPI_CHARACTER,buffersjj(:,SENDREQ,tid),&
                    buffersizejj,bufferposjj(SENDREQ,tid),&
                    MPI_COMM_WORLD,info)
 
      call mp_test('mp_pak_string failed',info) 

      return 
      end subroutine
!----------------------------------------------------------------------
      subroutine mp_pak_intg(tids,dest,type,scalar,info)

      implicit none


      include 'mpif.h'

      integer type,scalar,info,tids(*),dest,tid

      if( dest.eq.-1 ) then
        tid = numprocjj + 1
      else
        tid = tids(dest) + 1
      endif
       
!     write(200+myrank,*)  'mp_pak_intg ',scalar,' into pos ',bufferposjj(SENDREQ,tid)

      call MPI_PACK(scalar,1,type,buffersjj(:,SENDREQ,tid),&
                    buffersizejj,bufferposjj(SENDREQ,tid),&
                    MPI_COMM_WORLD,info)

      call mp_test('mp_pak_real failed',info)

      return
      end subroutine

!----------------------------------------------------------------------
      subroutine mp_pakv_intg_jj(type,array,nitem,stride,info)

      implicit none

      include 'mpif.h'

      integer type,nitem,stride,info

      integer array(*)

!      external pvmfpack, mp_test

!      call pvmfpack (type,array,nitem,stride,info)

      call mp_test('mp_pak_real failed',info)

      return
      end subroutine

!----------------------------------------------------------------------
      subroutine mp_upakv_intg_jj(type,array,nitem,stride,info)

      implicit none


!#include <fpvm3.h>

      integer type,nitem,stride,info

      integer array(*)

!      external pvmfpack, mp_test


!      call pvmfunpack (type,array,nitem,stride,info)

      call mp_test('mp_upakv_intg failed',info)

      return
      end subroutine

!----------------------------------------------------------------------
      subroutine mp_pakv_real_jj(type,array,nitem,stride,info)

      implicit none


!#include <fpvm3.h>

      integer type,nitem,stride,info

      real    array(*)

!      external pvmfpack, mp_test

!      call pvmfpack (type,array,nitem,stride,info)

      call mp_test('mp_pakv_real failed',info)

      return
      end subroutine
!----------------------------------------------------------------------

      subroutine mp_upakv_real_jj(type,array,nitem,stride,info)

      implicit none


!#include <fpvm3.h>

      integer type,nitem,stride,info

      real    array(*)

!      external pvmfunpack, mp_test

!      call pvmfunpack (type,array,nitem,stride,info)

      call mp_test('mp_upakv_real failed',info)

      return
      end subroutine

!----------------------------------------------------------------------
      subroutine mp_snd_others(np,my_dm,tids,tag,info)

      implicit none


      include 'mpif.h'

      integer np, tag, info, tids(*) ,my_dm, tid
      integer i, j, k

      integer STATUS(MPI_STATUS_SIZE)     

      do i=1,np
         tid = tids(i)+1
         if (i.ne.my_dm) then
            do j=1,buffersizejj
              buffersjj(j,SENDREQ,tid) = buffersjj(j,SENDREQ,numprocjj+1)
            enddo
            bufferposjj(SENDREQ,tid) = bufferposjj(SENDREQ,numprocjj+1)
!           write(200+myrank,*)'inq',i,buffersjj(1,SENDREQ,tid)
!      j = 0
!      call MPI_UNPACK(buffersjj(:,SENDREQ,tid),&
!                      buffersizejj,j,&
!                      k,1,MPI_INTEGER,MPI_COMM_WORLD,info)
!       write(200+myrank,*) 'test = ',k
         endif
      enddo

      do i=1,np
         tid = tids(i)+1
         if (i.ne.my_dm) then
!           write(200+myrank,*) 'MPI_ISEND to ',tids(i),bufferposjj(SENDREQ,tid)
            call MPI_ISEND(buffersjj(:,SENDREQ,tid),&
                           bufferposjj(SENDREQ,tid),&
                           MPI_PACKED,tids(i),tag,MPI_COMM_WORLD,&
                           requestsjj(SENDREQ,tid),info)
			
            call mp_test('mp_snd_others failed',info)

         endif
      enddo

      return
      end subroutine
!----------------------------------------------------------------------
      subroutine mp_recv_others(np,my_dm,tids,tag,info)

      implicit none

      include 'mpif.h'

      integer np, tag, info, tids(*) ,my_dm
      integer i, j, tid

      integer STATUS(MPI_STATUS_SIZE)     

      do i=1,np
         tid = tids(i)+1
         if (i.ne.my_dm) then
!           write(200+myrank,*) 'MPI_IRECV from ',tids(i)
            requestsjj(RECVREQ,tid) = MPI_REQUEST_NULL
            bufferposjj(RECVREQ,tid) = 0
            call MPI_IRECV(buffersjj(:,RECVREQ,tid),&
                           buffersizejj,MPI_PACKED,tids(i),tag, &
                           MPI_COMM_WORLD,requestsjj(RECVREQ,tid),info)
            call mp_test('mp_recv_others failed',info)

         endif
      enddo

      return
      end subroutine

!----------------------------------------------------------------------
      subroutine mp_snd_master(np,my_dm,tids,tag,info,SENDBUF,BUFSIZE)

      implicit none


      include 'mpif.h'

      integer np, tag, info, tids(*) ,my_dm

      integer i

      integer BUFSIZE,request
      character SENDBUF(*)
      integer STATUS(MPI_STATUS_SIZE)     

      request = 1
      if (1.ne.my_dm) then  
!            call MPI_SEND(SENDBUF,BUFSIZE,MPI_PACKED,tids(1),tag,&
!                        MPI_COMM_WORLD,info)
!non-blocking
        call MPI_ISEND(SENDBUF,BUFSIZE,MPI_PACKED,tids(1),tag, &
                        MPI_COMM_WORLD,request,info)
        call MPI_REQUEST_FREE (request,info)

        call mp_test('mp_snd_master failed',info)

      endif

      return
      end subroutine

!----------------------------------------------------------------------
      subroutine mp_upak_real(tids,from,type,scalar,info)

      implicit none


      include 'mpif.h'

      integer type, info,tids(*),from, tid

      real    scalar

      tid = tids(from)+1
      if( from.eq.-1 )   call abort()
      call MPI_UNPACK(buffersjj(:,RECVREQ,tid),&
                      buffersizejj,bufferposjj(RECVREQ,tid),&
                      scalar,1,type, MPI_COMM_WORLD,info)
!     write(200+myrank,*)  'mp_upak_real ',scalar,' from pos ',bufferposjj(RECVREQ,tid)

      call mp_test('mp_upak_real failed',info)

      return
      end subroutine
!----------------------------------------------------------------------  
      subroutine mp_upak_string(tids,from,str,strlen,info)

      implicit none

      include 'mpif.h'

      integer strlen,info,tids(*),from, tid
      character*80 str

      tid = tids(from)+1
      if( from.eq.-1 )   call abort()
      call MPI_UNPACK(buffersjj(:,RECVREQ,tid),&
                      buffersizejj,bufferposjj(RECVREQ,tid),&
                      str,strlen,MPI_CHARACTER,MPI_COMM_WORLD,info)
!     write(200+myrank,*)  'mp_upak_string ',str,strlen,' from pos ',bufferposjj(RECVREQ,tid),tids(from)
       
      call mp_test('mp_upak_string failed',info)

      return 
      end subroutine
 
!----------------------------------------------------------------------
      subroutine mp_upak_intg(tids,from,type,scalar,info)

      implicit none


      include 'mpif.h'

      integer type, info,tids(*),from, tid

      integer scalar
      
      tid = tids(from)+1
      if( from.eq.-1 )   call abort()
!      write(200+myrank,*)'inq',from,buffersjj(1,RECVREQ,tid)
      call MPI_UNPACK(buffersjj(:,RECVREQ,tid),&
                      buffersizejj,bufferposjj(RECVREQ,tid),&
                      scalar,1,type,MPI_COMM_WORLD,info)
!     write(200+myrank,*)  'mp_upak_intg ',scalar,' from pos ',bufferposjj(RECVREQ,tid),tids(from)

      call mp_test('mp_upak_intg failed',info)

      return
      end subroutine

!----------------------------------------------------------------------
      subroutine mp_check_mesg_arrived(tids,from,tag,info,&
                                    RECVBUF,BUFSIZE)

      implicit none


      include 'mpif.h'

      integer from,tag, info, tids(*)

      character RECVBUF(*)
      integer BUFSIZE,request
      integer STATUS(MPI_STATUS_SIZE)     
      request = 1
!      call MPI_RECV(RECVBUF,BUFSIZE,MPI_PACKED,tids(from),tag,&
!                 MPI_COMM_WORLD,STATUS,info)
!non-blocking
      call MPI_IRECV(RECVBUF,BUFSIZE,MPI_PACKED,tids(from),tag, &
                 MPI_COMM_WORLD,request,info)
      call MPI_WAIT(request, STATUS, info)
      call mp_test('mp_check_mesg_arrived failed',info)

      info=1

      return
      end subroutine
!----------------------------------------------------------------------
      subroutine mp_pakv_nabor(tids,dest,nosp,lsp,dtype,data,buf,info)

      implicit none


      include 'mpif.h'

      integer dtype, info,tids(*),dest,tid
      integer nosp,lsp(*)

      real      data(*), buf(*)
      integer i,j,k

      i = nosp

      if( dest.eq.-1 ) then
        tid = numprocjj + 1
      else
        tid = tids(dest) + 1
      endif

      if (i.gt.0) then

         do j=1,i
            k = lsp(j)
            buf(j) = data(k)
         enddo
      
!        write(200+myrank,*) 'mp_pakv_nabor ',i,'items at pos',bufferposjj(SENDREQ,tid)
         call MPI_PACK(buf,i,dtype,buffersjj(:,SENDREQ,tid),&
                       buffersizejj,bufferposjj(SENDREQ,tid),&
                       MPI_COMM_WORLD,info)

      endif

      return
      end subroutine
!----------------------------------------------------------------------
      subroutine mp_pakv_bdry(tids,dest,norp,lrp,dtype,data,buf,info)

      implicit none


      include 'mpif.h'

      integer dtype, info, tids(*), dest,tid
      integer norp ,lrp(*)

      real data(*), buf(*)
      integer i,j,k

      i = norp

      if( dest.eq.-1 ) then
        tid = numprocjj + 1
      else
        tid = tids(dest) + 1
      endif

      if (i.gt.0) then
         do j=1,i
            k = lrp(j) 
            buf(j) = data(k)
         enddo

      call MPI_PACK(buf,i,dtype,buffersjj(:,SENDREQ,tid),&
                    buffersizejj,bufferposjj(SENDREQ,tid),&
                    MPI_COMM_WORLD,info)

      endif

      return
      end subroutine

!----------------------------------------------------------------------
      subroutine mp_upakv_add_bdry(tids,from,norp,lrp,dtype,data,buf,info)

      implicit none


      include 'mpif.h'

      integer dtype, info 
      integer norp, lrp(*), tids(*), from, tid

      real data(*), buf(*)
      integer i,j,k

      tid = tids(from)+1
      if( from.eq.-1 )   call abort()
      i = norp
      if (i.gt.0) then
         call MPI_UNPACK(buffersjj(:,RECVREQ,tid),&
                         buffersizejj,bufferposjj(RECVREQ,tid),&
                         buf,i,dtype,&
                         MPI_COMM_WORLD,info)
         do j=1,i
            k = lrp(j) 
            data(k) = data(k) + buf(j)
         enddo
      endif


      return
      end subroutine

!----------------------------------------------------------------------
      subroutine mp_upakv_min_bdry_jj(tids,from,norp,lrp,dtype,data,buf,info)

      implicit none


      include 'mpif.h'

      integer dtype, info , tids(*), from

      real data(*), buf(*)
      integer i,j,k, tid
      integer norp , lrp(*)

      if( from.eq.-1 )   call abort()
      i = norp

      tid = tids(from)+1
      if (i.gt.0) then

         call MPI_UNPACK(buffersjj(:,RECVREQ,tid),&
                         buffersizejj,bufferposjj(RECVREQ,tid),&
                         buf,i,dtype,MPI_COMM_WORLD,info)

         do j=1,i
            k = lrp(j) 
            data(k) = min(data(k),buf(j))
         enddo
      endif


      return
      end subroutine

!----------------------------------------------------------------------
      subroutine mp_upakv_max_bdry_jj(tids,from,norp,lrp,dtype,data,buf,info)

      implicit none

      include 'mpif.h'

      integer dtype, info , tids(*), from 
      integer norp , lrp(*)

      real      data(*), buf(*)
      integer i,j,k,tid

      i = norp

      tid = tids(from)+1
      if (i.gt.0) then

         call MPI_UNPACK(buffersjj(:,RECVREQ,tid),&
                         buffersizejj,bufferposjj(RECVREQ,tid),&
                         buf,i,dtype,MPI_COMM_WORLD,info)

         do j=1,i
            k = lrp(j) 
            data(k) = max(data(k),buf(j))
         enddo
      endif


      return
      end subroutine

!----------------------------------------------------------------------
      subroutine mp_upakv_nabr(tids,from,nosp,lsp,dtype,data,buf,info)

      implicit none


      include 'mpif.h'

      integer dtype, info, tids(*), from 
      integer nosp , lsp(*)

      real      data(*), buf(*)
      integer i,j,k,tid

      i = nosp

      if( from.eq.-1 )  call abort()
      tid = tids(from)+1
      if (i.gt.0) then

         call MPI_UNPACK(buffersjj(:,RECVREQ,tid),&
                         buffersizejj,bufferposjj(RECVREQ,tid),&
                         buf,i,dtype,MPI_COMM_WORLD,info)

         do j=1,i
            k = lsp(j) 
            data(k) =  buf(j)
         enddo
      endif

      return
      end subroutine
!----------------------------------------------------------------------
        subroutine mp_test(message, ret_code)

        implicit none
        include 'mpif.h'

        character message*(*)
        integer ret_code

        if(ret_code .eq. MPI_SUCCESS) then
                return
        else
          print *, message
          call abort()
        endif

        end subroutine
!----------------------------------------------------------------------

      subroutine init_recv_others(np,my_dm,got_msg)

      implicit none

      integer i, np, got_msg(*) , my_dm

      do i=1,np
         got_msg(i) = 0
      enddo
         got_msg(my_dm) = 1
      return
      end subroutine

!----------------------------------------------------------------------
      subroutine init_recv_bdry(np,norp,got_msg)

      implicit none

      integer i,np,got_msg(*),norp(*)

      do i=1,np
         if (norp(i).eq.0) then
            got_msg(i) = 1
         else
            got_msg(i) = 0
         endif
      enddo
      return
      end subroutine

!----------------------------------------------------------------------
      subroutine init_recv_nabr(np,nosp,got_msg)


      implicit none

      integer i,np,nosp(*),got_msg(*)

      do i=1,np
         if (nosp(i).eq.0) then
            got_msg(i) = 1
         else
            got_msg(i) = 0
         endif
      enddo
      return
      end subroutine
!---------------------------------------------------------------
      subroutine wtime(timenow)
      
      implicit none
      include 'mpif.h'

      double precision :: timenow

      timenow=MPI_WTIME()

      return
      end subroutine
!---------------------------------------------------------------
      subroutine mp_barrier()
     
      implicit none
      include 'mpif.h'

      integer info

      call MPI_BARRIER(MPI_COMM_WORLD,info)

      return
      end subroutine

end module
