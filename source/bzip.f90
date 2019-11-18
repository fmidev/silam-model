!******************************************************************************
!****m* /bzip
! NAME
! module bzip
! NOTES
! This module provides a routine for reading and writing bz2-compressed files.
!******************************************************************************
module bzip

  !$ use OMP_LIB
  use, intrinsic :: ISO_C_BINDING, only: c_int, c_char, c_ptr, c_null_ptr, c_associated, c_size_t, c_loc, c_funptr, c_null_funptr
  use toolbox

  implicit none
  private

  public :: bzReadAll
  public :: bzReadPipe
  public :: bzReadPar

  interface

    ! Fortran interfaces for the C functions 'fopen' and 'fclose', defined in stdio.h

#ifdef VS2012
    type(c_ptr) function pOpen (path, mode) bind (c, NAME='_popen')
#else
    type(c_ptr) function pOpen (path, mode) bind (c, NAME='popen')
#endif
      use, intrinsic :: ISO_C_BINDING
      character(kind=c_char) :: path(*), mode(*)
    end function

#ifdef VS2012
    integer(c_int) function pClose (file) bind (c, NAME='_pclose')
#else
    integer(c_int) function pClose (file) bind (c, NAME='pclose')
#endif
      use, intrinsic :: ISO_C_BINDING
      type(c_ptr), value :: file
    end function


    integer(c_int) function fRead (buf, itemsz, nitems, file) bind (c, NAME='fread')
      use, intrinsic :: ISO_C_BINDING
      integer(c_size_t), value :: itemsz, nitems
      type(c_ptr), value :: buf, file
    end function

    type(c_ptr) function fOpen (path, mode) bind (c, NAME='fopen')
      use, intrinsic :: ISO_C_BINDING
      character(kind=c_char) :: path(*), mode(*)
    end function

    integer(c_int) function fClose (file) bind (c, NAME='fclose')
      use, intrinsic :: ISO_C_BINDING
      type(c_ptr), value :: file
    end function

    integer(c_int) function fEOF (file) bind (c, NAME='feof')
      use, intrinsic :: ISO_C_BINDING
      type(c_ptr), value :: file
    end function

    integer(c_int) function c_ftell (f) bind (c, NAME='ftell')
      use, intrinsic :: ISO_C_BINDING
      type(c_ptr), value :: f
    end function

    subroutine c_memcpy (dst,src,n) bind(c, NAME='memcpy')
      use, intrinsic :: ISO_C_BINDING
      INTEGER(c_intptr_t), intent(in):: dst, src   ! seems to be the standard
!      type(c_ptr), value :: dst, src       ! not allowed by Intel
      integer(kind=c_size_t), value :: n
    end subroutine

    ! Fortran interfaces for the functions from libbz2, defined in bzlib.h
#ifdef WITH_BZIP2


    type(c_ptr) function bzReadOpen (bzerror, f, verbosity, small, unused, nUnused) bind(c, NAME='BZ2_bzReadOpen')
      use, intrinsic :: ISO_C_BINDING
      integer(c_int) :: bzerror
      type(c_ptr), value :: f, unused
      integer(c_int), value :: verbosity, small, nUnused
    end function


    integer(c_int) function bzRead (bzerror, f, buf, len) bind (c, NAME='BZ2_bzRead')
      use, intrinsic :: ISO_C_BINDING
      integer(c_int) :: bzerror
      type(c_ptr), value :: f
      type(c_ptr), value :: buf
      integer(c_int), value :: len
    end function

    subroutine bzReadClose (bzerror, f) bind (c, NAME='BZ2_bzReadClose')
      use, intrinsic :: ISO_C_BINDING
      integer(c_int) :: bzerror
      type(c_ptr), value :: f
    end subroutine

    subroutine bzReadGetUnused (bzerror, f, unused, nUnused) bind (c, NAME='BZ2_bzReadGetUnused')
      use, intrinsic :: ISO_C_BINDING
      integer(c_int) :: bzerror
      type(c_ptr), value :: f
      type(c_ptr) :: unused
      integer(c_int) :: nUnused
    end subroutine

    type(c_ptr) function bzWriteOpen (bzerror, f, blockSize100k, verbosity, workFactor) bind (c, NAME='BZ2_bzWriteOpen')
      use, intrinsic :: ISO_C_BINDING
      integer(c_int) :: bzerror
      type(c_ptr), value :: f
      integer(c_int), value :: blockSize100k, verbosity, workFactor
    end function

    subroutine bzWrite (bzerror, f, buf, len) bind (c, NAME='BZ2_bzWrite')
      use, intrinsic :: ISO_C_BINDING
      integer(c_int) :: bzerror
      type(c_ptr), value :: f
      character(kind=c_char) :: buf(*)
      integer(c_int), value :: len
    end subroutine

    subroutine bzWriteClose (bzerror, f, abandon, nbytes_in, nbytes_out) bind (c, NAME='BZ2_bzWriteClose')
      use, intrinsic :: ISO_C_BINDING
      integer(c_int) :: bzerror, nbytes_in, nbytes_out
      type(c_ptr), value :: f
      integer(c_int), value :: abandon
    end subroutine

    integer (c_int) function bzDecompInit (strm, verbosity, small) bind (c, NAME='BZ2_bzDecompressInit')
      use, intrinsic :: ISO_C_BINDING
      type(c_ptr), value :: strm
      integer(c_int), value :: verbosity, small
    end function

    integer (c_int) function bzDecomp (strm) bind (c, NAME='BZ2_bzDecompress')
      use, intrinsic :: ISO_C_BINDING
      type(c_ptr), value :: strm
    end function

    integer (c_int) function bzDecompEnd (strm) bind (c, NAME='BZ2_bzDecompressEnd')
      use, intrinsic :: ISO_C_BINDING
      type(c_ptr), value :: strm
    end function



#endif



  end interface

!!!!!! from bzlib.h
!!!!!typedef
!!!!!   struct {
!!!!!      char *next_in;
!!!!!      unsigned int avail_in;
!!!!!      unsigned int total_in_lo32;
!!!!!      unsigned int total_in_hi32;
!!!!!
!!!!!      char *next_out;
!!!!!      unsigned int avail_out;
!!!!!      unsigned int total_out_lo32;
!!!!!      unsigned int total_out_hi32;
!!!!!
!!!!!      void *state;
!!!!!
!!!!!      void *(*bzalloc)(void *,int,int);
!!!!!      void (*bzfree)(void *,void *);
!!!!!      void *opaque;
!!!!!   } 
!!!!!   bz_stream;

  INTEGER(C_INT), parameter :: BZ_OK = 0, BZ_STREAM_END = 4


  TYPE, bind(C) ::  BZ_STREAM
    type (C_PTR)  :: next_in = C_NULL_PTR
    INTEGER(C_INT) :: avail_in=0, total_in_lo32=0, total_in_hi32=0
    type (C_PTR)  :: next_out = C_NULL_PTR
    INTEGER(C_INT) :: avail_out=0, total_out_lo32=0, total_out_hi32=0
    type (C_PTR)  :: state = C_NULL_PTR
    type (C_FUNPTR)  :: bzalloc = C_NULL_FUNPTR, bzfree = C_NULL_FUNPTR , opaque=C_NULL_FUNPTR
  END TYPE BZ_STREAM
  
  TYPE (BZ_STREAM), private, parameter :: BZ_STREAM_missing= BZ_STREAM(&
          & C_NULL_PTR, 0,0,0, C_NULL_PTR, 0,0,0, C_NULL_PTR, C_NULL_FUNPTR,C_NULL_FUNPTR, C_NULL_FUNPTR)



contains




  !****************************************************************************
  subroutine bzReadAll (fname, buf, length )
    character(len=*), intent(in) :: fname
    character, dimension(:), allocatable, target, intent(inout) :: buf !! Allocation status is "in"
    integer(kind=8), intent(out) :: length

    integer(c_size_t) :: bufsz
    integer(c_int) :: err, nUnused, cnt
    integer(kind=8) :: size8

    type(c_ptr) :: handle, bzHandle           ! file handles
    character(kind=c_char,len=1000) :: cFname   ! C file name
    logical :: eof = .false.                  ! end of file reached?
    integer, parameter :: inisize=100000
    character, dimension(:), allocatable :: buftmp
    integer, parameter :: BZ_MAX_UNUSED=5000
    type(c_ptr) :: unusedptr
    character(kind=C_CHAR, len=BZ_MAX_UNUSED),target :: unused


    character(kind=c_char,len=2), parameter :: mode = "r"//achar(0)


#ifdef WITH_BZIP2
    cFname=trim(fname)//achar(0)
    handle = fOpen(cFname,mode)

    if (.not. c_associated(handle)) then
      write(*,*) "Error opening file ", fname
      stop
    end if


    if (allocated(buf)) then
        bufsz=size(buf)
    else
        bufsz=inisize
        allocate(buf(inisize))
    endif

    nUnused = 0
    length = 0
    err = 4  

    do while (err == 0 .or. err == 4)
      
        if (err==4) then
          bzHandle = bzReadOpen(err,handle,0,0,c_loc(unused),nUnused)
          if (.not. c_associated(bzHandle) .or. err/=0) then
            write(*,*) "Error opening bz file ",fname
            write(*,*) "Error code: ",err
            stop
         end if
        endif

        cnt = bzRead(err,bzHandle,c_loc(buf(length+1)),min(int(bufsz-length), huge(err)))
        length = length + cnt

!        print *, "err, cnt", err,cnt

        if (err /= 0 .and. err /= 4) then
          write(*,*) "Error reading bz file", fname
          write(*,*) "Error code: ",err
          stop
        endif

        if ( length == bufsz) then
            bufsz = bufsz * 3 / 2
          !  print *, "New size", bufsz
            allocate(buftmp(bufsz))
            buftmp(1:length) = buf(1:length)
            deallocate(buf)
            call move_alloc(buftmp,buf)
          !  print *, "relocation done!"
        endif

        if (err==4) then
          call bzReadGetUnused ( err, bzHandle, unusedptr,nUnused)
          size8 = nUnused
          !
          !
          call c_memcpy(loc(unused), loc(unusedptr), size8) 
          ! Seems to be a GNU extention or alike. DOES NOT WORK with Intel
!          call c_memcpy(c_loc(unused), unusedptr, size8) 

          if (err/=0) then
            write(*,*) "Error bzReadGetUnused ",fname
            write(*,*) "Error code: ",err
            stop
          end if
          call bzReadClose(err, bzHandle)
          if (err/=0) then
            write(*,*) "Error closing bz file ",fname
            write(*,*) "Error code: ",err
            stop
          end if
          if (fEOF(handle) /= 0) exit
          err = 4 !Restore err 

        endif
     enddo
     err = fClose(handle)
     if (err /= 0) then
       write(*,*) "Error closing file ",fname
       write(*,*) "Error code: ",err
       stop
     end if
#else
    write(*,*) "Compiled with no bzip2 support"
    stop
#endif 

  end subroutine


  !****************************************************************************
  subroutine bzReadPipe (fname, buf, length)
    character(len=*), intent(in) :: fname
    character, dimension(:), allocatable, target, intent(inout) :: buf
    integer (kind=8), intent(out) :: length

    integer(kind=c_size_t) :: bufsz, cnt
    integer :: err, itry, ierrno

    type(c_ptr) :: handle           ! file handles
    character(kind=c_char,len=1000) :: cFname   ! C file name
    character(len=700) :: errorstr
    integer, parameter :: inisize=100000
    character, dimension(:), allocatable :: buftmp

    character(kind=c_char,len=2), parameter :: mode = "r"//achar(0)
    integer(kind=c_size_t), parameter :: one = 1


    !cFname="/lustre/apps/silam/bin/lbzip2 -v -d -c "//trim(fname)//achar(0)
    cFname="lbzip2 -d -c "//trim(fname)//achar(0)

    if (allocated(buf)) then
        bufsz=size(buf)
    else 
        print *, "non-allocated buffer"
        bufsz=inisize
        allocate(buf(inisize))
    endif


    handle = c_null_ptr

    !
    ! in Voima lbzip2 sometimes segfaults in arbitratry times 
    ! We have to retry....
   

    do itry=1,10
      write(*,*) "Popen....", itry
      handle = pOpen(cFname,mode)
      write(*,*) "trying handle"
      if (.not. c_associated(handle)) continue
      write(*,*) "handle ok"


      length=0
      do while (fEOF(handle) == 0)
        cnt = fRead(c_loc(buf(length+1)), one, bufsz-length, handle)
        write (*,*) "Read coun", cnt
        length = length + cnt


        if ( length == bufsz) then
            bufsz = bufsz * 3 / 2
          !  print *, "New size", bufsz
            allocate(buftmp(bufsz))
            buftmp(1:length) = buf(1:length)
            deallocate(buf)
            call move_alloc(buftmp,buf)
          !  print *, "relocation done!"
        endif


     enddo
     err = pClose(handle)
     if (err == 0) exit 

       write(*,*) "Error closing file ",cFname
     write(*,*) "pClose returned: ",err
     err= ierrno()
     write(*,*) "After poprn errno", err
     call gerror(errorstr)
     write(*,*) "Error string:", trim(errorstr)
     write(*,*) "ls -l /lustre/apps/silam/bin/lbzip2"
     call system("ls -l /lustre/apps/silam/bin/lbzip2")
     call sleep(1)
     
    enddo

    if (itry>9) then
      write(*,*) "Error opening file ", cFname
       length = -1
      return
     end if

    write(*,*) "Read ok, nbytes", length

  end subroutine bzReadPipe



  !****************************************************************************
  subroutine bzReadPar (fname, buf, length )
    character(len=*), intent(in) :: fname
    character, dimension(:), allocatable, target, intent(inout) :: buf !! Allocation status is "in"
    integer(kind=8), intent(out) :: length !!Number of bytes read


    character, dimension(:), allocatable, target :: bzraw, bufOut
    type (BZ_STREAM), target :: bzstream
    integer (kind=8) :: file_size

    character, dimension(:), allocatable :: buftmp
    integer(kind=8), dimension(:), allocatable :: blockOffsets
    integer ::  totalblocks, nthreads, ithread, max_offsets
    integer, dimension(max_threads) ::  threadblocks, threadblockoff
    integer(kind=8) :: blocksize, chunksz, bufsz, cnt, lTmp, length_loc
    type (c_ptr) :: inloc,outloc
    integer :: iTmp, jTmp, kTmp
    integer :: istatus, iUnit
    integer, parameter :: inisize=100000
    integer, parameter :: max_blocks=10000
    character, dimension(6), parameter :: blockmagick = (/'1','A','Y','&','S','Y'/)
    character, dimension(3), parameter :: filemagick = (/'B','Z','h'/)
    character(len=*), parameter :: sub_name = "bzReadPar";



#ifdef WITH_BZIP2

    if (allocated(buf)) then
        bufsz=size(buf)
    else 
        print *, "non-allocated buffer"
        bufsz=inisize
        allocate(buf(inisize))
    endif

!
!
!   Get the whole file into array


     INQUIRE(FILE=fname, SIZE=file_size)
     allocate(bzraw(file_size))
     iUnit = fu_next_free_unit()
     if(error)return
     open(iUnit, file=fname, status='old', action='read', form="unformatted",  ACCESS='STREAM', iostat= iStatus)
     if(iStatus /= 0)then
       call set_error('Failed to open GRIB file:' + fName, sub_name)
       return
     endif
     read(iUnit) bzraw(1:file_size)
     close (iUnit)

     if (all(bzraw(1:3) == filemagick )) then
        blocksize = (ichar(bzraw(4))-ichar('0'))*100000  !Size of decompressed block
     else
        call set_error("Bad magick in file "//trim(fname), sub_name)
        return
     endif

     allocate(blockOffsets(max_blocks))

     length = 0
     

     !$OMP PARALLEL if (.true.) default(none) shared(blockOffsets,totalblocks, file_size, &
     !$OMP       &threadblockoff, threadblocks, bzraw, buf, buftmp, bufsz, cnt, bufOut, blocksize, length, error) &
     !$OMP & private(nthreads, ithread, iTmp, jTmp, max_offsets, chunksz, bzstream, iStatus, lTmp, length_loc,inloc,outloc)
    
       iThread = 0
       nthreads = 1 
       !$ iThread = OMP_GET_THREAD_NUM()
       !$ nthreads = omp_get_num_threads()
       length_loc = 0

     max_offsets =  max_blocks / nthreads ! max block offsets per thread
     threadblockoff(iThread+1) = (max_offsets * ithread) 
     threadblocks(iThread+1) = 0

     chunksz = (file_size - 1) / nthreads + 1
     do lTmp = 1+iThread*chunksz,min(file_size-5,(iThread+1)*chunksz)
        if (all(bzraw(lTmp:lTmp+5) == blockmagick)) then
          threadblocks(iThread+1) = threadblocks(iThread+1) + 1
          blockOffsets(threadblockoff(iThread+1) + threadblocks(iThread+1)) = lTmp
        endif
     enddo     

     !$OMP BARRIER 
      ! Make joint list of block offsets
     !$OMP MASTER
       allocate(bufOut(nthreads*blocksize))

       totalblocks = threadblocks(1)  !master blocks of masterthread  are in place
       do iTmp=1,nthreads-1  !other thread
         do jTmp = 1,threadblocks(iTmp+1)
            totalblocks = totalblocks + 1
            blockOffsets(totalblocks) = blockOffsets(threadblockoff(iTmp+1) + jTmp )
         enddo
       enddo
       blockOffsets(totalblocks+1) = file_size ! End of last block
     !!  write (*,*) "Total blocks:", totalblocks
     !  write (*,*) "Offsets:", blockOffsets(1:totalblocks+1)


       !!!Allocate needed buffer for output
       lTmp = blocksize*totalblocks
       if (bufsz <  lTmp) then
         call msg(sub_name//": Enlarging buf from, to, kB ", int(bufsz/1024), int(lTmp/1024))
         deallocate(buf)
         allocate(buf(lTmp))
!            allocate(buftmp(lTmp))
!            buftmp(1:bufsz) = buf(1:bufsz)
!            deallocate(buf)
!            call move_alloc(buftmp,buf)
         bufsz =  lTmp
       endif
     !$OMP END MASTER
     !$OMP BARRIER
     !!! Now blockOffsets(1:totalblocks+1)  contains offsets of bzipped blocks 

     bzstream = bz_stream_missing 
     iStatus = bzDecompInit (c_loc(bzstream), 0, 0)

     bzstream%next_in = c_loc(bzraw(1))
     inloc = c_loc(bzraw(1))
     bzstream%avail_in = 4 ! Only header

     outloc=c_loc(bufOut(iThread*blocksize+1))
     bzstream%next_out  = c_loc(bufOut(iThread*blocksize+1))
     bzstream%avail_out = 0
     bzstream%total_out_lo32 = 0
     iStatus = bzDecomp(c_loc(bzstream)) 

!     print *, "After init Thread, inloc, outloc", ithread, diffptr(bzstream%next_in , inloc), diffptr(bzstream%next_out , outloc)


     do iTmp=0, (totalblocks-1) / nthreads  !Workc chunk loop, starting from 0!
          if (error) exit
          threadblocks(ithread+1) = 0 

          jTmp = iTmp*nthreads+ithread + 1  !Number of my block

          if (jTmp <= totalblocks) then
            bzstream%next_in = c_loc(bzraw(blockOffsets(jtmp)))
            !bzstream%avail_in = blockOffsets(jtmp+1) - blockOffsets(jtmp)
            bzstream%avail_in = blockOffsets(jtmp+1) - blockOffsets(jtmp)
            inloc=bzstream%next_in

            bzstream%next_out  = c_loc(bufOut(iThread*blocksize+1))
            outloc = bzstream%next_out
            bzstream%avail_out = blocksize
            bzstream%total_out_lo32 = 0

 !           lTmp = c_loc(bzstream%next_out)
            iStatus = bzDecomp(c_loc(bzstream))
            ! Length of decompressed chunk for this thread
            if (iStatus >= 0) then 
              threadblocks(ithread+1) = bzstream%total_out_lo32
            else
              call set_error("Bad bzDecomp status",sub_name)
            endif
            !print *, "Thread, block status", ithread, jTmp, iStatus,  threadblocks(ithread+1), loc(bzstream%next_out) - lTmp
 !           print *, "Thread, inloc, inloc1, outadvance", ithread, diffptr(inloc, c_loc(bzraw(1))) - blockOffsets(jtmp), &
 !             &  diffptr(bzstream%next_in, c_loc(bzraw(1)))- blockOffsets(jtmp+1) , diffptr(bzstream%next_out,  outloc),  bzstream%avail_out
          endif

          !$OMP BARRIER
          ! threadblocks(1:nthreads) must be coherent
!          if (iThread == 0 ) print *, "Chunk, threadblocks", iTmp, threadblocks(1:nThreads)

          if (jTmp <= totalblocks .and. threadblocks(ithread+1) > 0) then
            lTmp=length_loc + sum(threadblocks(1:ithread))  !!offset of first byte of the chunk
         !   print *, "Thread, block status", ithread, jTmp, iStatus,  threadblocks(ithread+1), lTmp+1,":",lTmp+threadblocks(ithread+1)
            buf(lTmp+1:lTmp+threadblocks(ithread+1)) = bufOut(iThread*blocksize+1: iThread*blocksize+threadblocks(ithread+1))
         !   if (any(buf(lTmp+1:lTmp+threadblocks(ithread+1)) /= &
         !       &bufOut(iThread*blocksize+1:iThread*blocksize+threadblocks(ithread+1)))) then
         !       do kTmp=1,threadblocks(ithread+1)
         !           if (buf(lTmp+kTmp) /= bufOut(iThread*blocksize+kTmp)) call ooops("Found!")
         !       enddo
         !   
            !endif
          endif
          length_loc = length_loc + sum(threadblocks(1:nthreads))
          !threadblocks(1:nthreads) must be inact until now
          !$OMP BARRIER  

     enddo !iTmp 0, totalblocks / nthreads

     iStatus = bzDecompEnd (c_loc(bzstream))
     if (iThread == 0 ) length = length_loc

     !$OMP END PARALLEL
     
     deallocate(blockOffsets)
     deallocate(bufOut)
     deallocate(bzraw)


     if (error) call set_error("Error after decompression of "//trim(fname),sub_name)
#else
    call set_error("Compiled with no bzip2 support",sub_name)
#endif 

  end subroutine

  integer (kind=8) function diffptr(a,b)
    type (c_ptr) :: a,b


    diffptr = (transfer(a, diffptr) - transfer(b, diffptr))
      
  endfunction
end module
