MODULE grib_api_io
  !
  ! Contains the tools needed for handling the GRIB-records and 
  ! GRIB-files in SILAM. Suitable both for reading and writing GRIB files.
  ! 
  ! Author: Mikhail Sofiev e-mail mikhail.sofiev@fmi.fi
  !
  ! The module uses GRIB API of ECMWF
  !
  ! All units: SI
  ! 
  ! All degrees: real numbers in degrees (degrees and parts per
  ! hundred/parts per thousand , NOT in degrees and minutes
  ! 
  ! Language: ANSI standard Fortran 90
  ! 
  ! Modules used:

  !use netcdf_io
  USE grib_api
  USE grads_io
  USE grib_code_table
  use bzip

  IMPLICIT NONE

  private


  private fu_asimof_field          ! READING
  private fu_grib_grid             ! GRIB grid to the SILAM one
  public print_gribheadings
  public grib_tests


  ! Reading stuff from GRIB
  public open_grib_file_i
  public id_list_from_grib_file
  public read_field_from_grib_file
  PUBLIC close_grib_file_i
  public fu_grib_filename

!  PUBLIC read_next_header_from_gribfile ! with NO data decoding
  private get_data_from_grib_index    ! get the data for the field whose header is known
!  public release_grib_buffer

  ! Write stuff to GRIB
  PUBLIC open_gribfile_o                ! file for GRIB output 
  public switch_grib_binary_o           ! New binary for GRIB output
  PUBLIC write_next_field_to_gribfile   ! with coding
  PUBLIC close_gribfile_o


  ! The private functions and subroutines not to be used elsewhere:
  ! Reading
  private gribheadings_to_field_id     ! from GRIB internal to field_id
  PRIVATE get_grib_times
  PRIVATE fu_grib_level
  private fix_strange_quantities 
  private fu_if_replace_real_missing
  private fix_scan_mode
  ! writing
  PRIVATE field_id_to_gribheadings
  PRIVATE set_field_level_to_grib
  PRIVATE set_field_time_to_grib
!  PRIVATE fu_get_10_scale
  private set_field_grid_to_grib


  integer, private, parameter :: pbopen_string_len = 250

#ifdef DEBUG_GRIB  
  logical, private, parameter :: ifDebug = .True.
#else
  logical, private, parameter :: ifDebug = .False.
#endif


  type grib_input_file
    private
    character, dimension(:), allocatable :: grib_raw
    integer(kind=8) :: raw_len
    integer :: nmsgs
    integer(kind=8), dimension(max_2d_fields) :: grib_offsets
    TYPE(silja_field_id), dimension(max_2d_fields) :: id
    type(silja_interval) :: obstime_interval
    logical :: idList_ready = .False.
    integer, dimension(max_2d_fields) :: indGrib ! gribIndex of the field
    real,    dimension(max_2d_fields) :: scale_factor !Scale to apply on reading
    logical, dimension(max_2d_fields) :: ifZeroNegatives ! ZeroNegativesonReading
    real, dimension(:), allocatable   :: vals_my !! Cached values
    real, dimension(:), allocatable   :: vals_others !! Cached values from other members
    integer, dimension(max_2d_fields) :: who_has !smpi_adv_rank of a member who has the values cached
    integer, dimension(max_2d_fields) :: vals_offset !MPI member who has the values cached
    CHARACTER (LEN=fnlen) :: fname
    logical :: defined = .False.
  end type grib_input_file

  INTEGER, public, PARAMETER :: max_nbr_of_iGRIB_files = 10

 type(grib_input_file), dimension(max_nbr_of_iGRIB_files), target :: iGRIB





CONTAINS

!************************************************************************
!
!     Reading of the GRIB file (including private stuff)
!
!************************************************************************


  SUBROUTINE open_grib_file_i(fname, grib_unit, obstime_interval)

    ! Description:
    ! Opens a (compressed) gribfile and returns it's unit (output !!).
    ! Reads its contents to grib_raw array and selects individual messages
    !
    ! ATTENTION. Grib_unit is NOT a FORTRAN unit, but an index in iGRIB array
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    CHARACTER (LEN=*), INTENT(in) :: fname
    type(silja_interval), intent(in) :: obstime_interval

    ! Imported parameters with intent OUT:
    INTEGER, INTENT(out) :: grib_unit

    ! Local declarations:
    INTEGER :: io_status, namelen

    INTEGER :: inq_unit, iRet, iUnit,iStatus
    integer (kind=8) :: file_size, i, j, msg_size
    character, dimension(:), pointer :: grib_raw
    integer(kind=8), dimension(:), pointer :: grib_offsets
    type(grib_input_file), pointer :: gf
    integer(kind=8) :: count0, count1, count2, count3, count4 !!!Timers

    character (len=*), parameter :: sub_name="open_grib_file_i"

    integer :: grib_mpi_rank, grib_mpi_tasks

    grib_mpi_rank  = smpi_adv_cart_rank
    grib_mpi_tasks = smpi_adv_tasks
    !grib_mpi_rank  = 0
    !grib_mpi_tasks = 1

    ! 


    CALL SYSTEM_CLOCK(count0)

    do grib_unit = 1,max_nbr_of_iGRIB_files
      if (.not. iGRIB(grib_unit)%defined) exit
    enddo

    if (grib_unit > max_nbr_of_iGRIB_files) then
      call set_error("max_nbr_of_iGRIB_files exceeded", sub_name)
      call msg("")
      call msg("Already opened iGRIB_files:")
      do grib_unit = 1,max_nbr_of_iGRIB_files
        call msg(iGRIB(grib_unit)%fname)
      enddo
      call msg("")
      grib_unit = int_missing
      return
    endif
    gf => iGRIB(grib_unit)
    gf%fname = fname
    grib_offsets => gf%grib_offsets

    namelen=(len(trim(fname)))
    !get everyone a copy of the (decompressed) grib

    count1=count0
    count2=count0
    if (grib_mpi_rank == 0) then
       if (fname(namelen-3:namelen) == '.bz2') then
         !BZIP file
#ifdef WITH_BZIP2
        call msg('Opening GRIB.bz2 file i with bzReadPar: '//trim(fname))
         call bzReadPar(fname, gf%grib_raw, file_size)
#else         
        call msg('Opening GRIB.bz2 file i with bzReadPipe: '//trim(fname))
         call bzReadPipe(fname, gf%grib_raw, file_size)
#endif         
       else
         ! regular grib file
        call msg('Opening GRIB file i: '//trim(FName))
         INQUIRE(FILE=fname, SIZE=file_size)
         if (allocated(gf%grib_raw)) then
           if (size(gf%grib_raw) < file_size) then
              call msg("reallocating gf%grib_raw to (kB)", int((file_size/10+1)*11/1024))
             deallocate(gf%grib_raw)
             allocate(gf%grib_raw((file_size/10+1)*11)) !Some 10% more, so next file fits
           endif
         else
             call msg ("allocating gf%grib_raw size (kB)", int((file_size/10+1)*11/1024))
             allocate(gf%grib_raw((file_size/10+1)*11))
         endif
         
         iUnit = fu_next_free_unit()
         if(error)return
         open(iUnit, file=fname, status='old', action='read', form="unformatted",  ACCESS='STREAM', iostat= iStatus)
         if(iStatus /= 0)then
           call set_error('Failed to open GRIB file:' + fName, sub_name)
           return
         endif
         do i = 0,100  !!!Or 250G should be enough
            if (i*MAX_INT32 > file_size) exit
            j=min((i+1)*MAX_INT32, file_size)
!            print *, "reading total bytes", i*recsize+1, "to", j 
            read(iUnit) gf%grib_raw(i*MAX_INT32+1:j)
         enddo
         close (iUnit)

       endif
       CALL SYSTEM_CLOCK(count1) !!Grib stream reading complete
       gf%raw_len = file_size
       if (file_size < 1) then
          call set_error('Empty GRIB file:' + fName, sub_name)
          return
       endif

       !
       ! Find messages in the stream
       !

       j = 1 !index of message start
       do i=1,max_2d_fields-2 !Message counter
           do while (j < file_size - 8 ) 
              grib_raw => gf%grib_raw(j:j+8)
              if (all(grib_raw(1:4) == (/'G','R','I','B'/))) then

                    if (grib_raw(8) == char(1)) then !grib1
                            !Here is the trick missing from any grib1 standards i have seen
                            !have found. It is, however, implemented in libcdo
                            !and grib_api. Here is just a lose copy of that

                      msg_size = iand(iachar(grib_raw(5)),127)  !!Highest bit means rescaling
                      msg_size =  msg_size*256 + iachar(grib_raw(6)) 
                      msg_size =  msg_size*256 + iachar(grib_raw(7))
                      if (grib_raw(5) > char(127)) then
                              msg_size = msg_size*120
                      endif
                    elseif (grib_raw(8) == char(2)) then !grib2
                      grib_raw => gf%grib_raw(j:j+16)
                      msg_size = iachar(grib_raw(9))*256 + ichar(grib_raw(10))
                      msg_size = (msg_size*256 + iachar(grib_raw(12)))*256 + ichar(grib_raw(11))
                      msg_size = (msg_size*256 + iachar(grib_raw(13)))*256 + ichar(grib_raw(14))
                      msg_size = (msg_size*256 + iachar(grib_raw(15)))*256 + ichar(grib_raw(16))
                    else
                      Call msg("Oops, grib recognition failed")
                      call set_error("Broken GRIB file?", sub_name)
                      return
                    endif
                    grib_offsets(i) = j-1 !Strat of current message
                    j = j + msg_size
                    grib_offsets(i+1) = j-1 !First offset after current message
                    gf%nmsgs = i
                    exit ! next i
              else
           !       call msg("Seeking", i, j)
                    j = j + 1 !seek
              endif
           enddo

           if (j >= file_size - 8 )  exit
       enddo

       ! Check and return.....
       if (i > max_2d_fields-2) then
          call msg("max_2d_fields",max_2d_fields)
          call msg("chFName="+fname)
          call set_error("Too many mesages in grib file",sub_name)
          return
       endif
       !!! grib_offsets(gf%nmsgs+1) is the next byte after the end of last message
       grib_offsets(gf%nmsgs+2:max_2d_fields)=int_missing 
       grib_offsets(max_2d_fields-1) = gf%nmsgs
       grib_offsets(max_2d_fields) = file_size
     endif !!rank0
     CALL SYSTEM_CLOCK(count2) !!Grib stream reading complete
     
     if (grib_mpi_tasks > 1) then
         call smpi_bcast_int8_aray(grib_offsets, max_2d_fields, smpi_adv_comm, 0) 
         if (grib_mpi_rank /= 0) file_size = grib_offsets(max_2d_fields)
     endif

     CALL SYSTEM_CLOCK(count3) !!Grib stream labels echange


     if (file_size<1) then
       call set_error("zero-sized file",sub_name)
       return
     endif

     if (grib_mpi_rank /= 0) then
         file_size = grib_offsets(max_2d_fields)
         gf%nmsgs =  grib_offsets(max_2d_fields-1)
         if (allocated(gf%grib_raw)) then
           if (size(gf%grib_raw) < file_size) then
             deallocate(gf%grib_raw)
             allocate(gf%grib_raw((file_size/10+1)*11)) !Some 10% more, so next file fits
           endif
          else
              allocate(gf%grib_raw((file_size/10+1)*11))
          endif
     endif
     if (grib_mpi_tasks > 1) &
          call smpi_bcast_char(gf%grib_raw, file_size,  0)

     gf%idList_ready = .false.
     gf%who_has(1:gf%nmsgs) = int_missing
     gf%obstime_interval = obstime_interval !Needed to fix_strange_quantities

     gf%defined = .true.

      CALL SYSTEM_CLOCK(count4)
      call msg("Opened GRIB file "//trim(fu_str(grib_unit))//": "//trim(fname)//" Memory usage (kB), time (ms)", &
           &  fu_system_mem_usage(), 1e-6*(count4-count0))

      call msg("open_grib_file_i timings (ms): read_stream, parse_stream, lables_xcg, stream_xcg", &
            & (/1e-6*(count1-count0), 1e-6*(count2-count1), 1e-6*(count3-count2), 1e-6*(count4-count3)/))
  
  !   if (allocated(gf%grib_raw)) then
  !       call msg("size(gf%grib_raw)",real(size(gf%grib_raw, kind=8),kind=8))
  !   else
  !     call msg ("gf%grib_raw not allocated")
  !   endif


  END SUBROUTINE open_grib_file_i


  subroutine id_list_from_grib_file(grib_unit, idList, iCount)

    implicit none
    integer, intent(in) :: grib_unit
    integer, intent(out) :: iCount
  
    type(silja_field_id),dimension(:), pointer :: idList  !Returns a pointer

    character, dimension(:), pointer :: grib_raw
    integer(kind=8), dimension(:), pointer :: grib_offsets
    type(grib_input_file), pointer :: gf
    integer :: iMsg, io_status

    character (len=*), parameter :: sub_name="id_list_from_grib_file"

    !$ if (omp_in_parallel()) then
    !$  call set_error("Can not be called from parallel", sub_name)
    !$ endif


    idList => null()
    iCount = int_missing

    if (grib_unit < 0 .or. grib_unit > max_nbr_of_iGRIB_files) then
        call set_error("bad grib_unit = "//trim(fu_str(grib_unit)), sub_name)
        return
    endif
    gf => iGRIB(grib_unit)
    if ( .not. gf%defined ) then
        call set_error("Not opened grib_unit = "//trim(fu_str(grib_unit)), sub_name)
        return
    endif




      !Make idList

    if (.not. gf%idList_ready) then
      grib_offsets => gf%grib_offsets
 !       call msg("")
 !       call msg("")
 !       call msg("")
#ifdef VS2012 
      ! No parallel GRIB handling at Windows
      !$OMP PARALLEL if (.false.) &  
#elif DEBUG_GRIB
      ! Parallel debug messages create a mess
      !$OMP PARALLEL if (.false.) &  
#else
      !$OMP PARALLEL if (gf%nmsgs > 10) & 
#endif
      !$OMP &  default(none) shared(gf,error,grib_offsets)&
      !$OMP & private(grib_raw, iMsg, io_status)

      !$OMP DO
      do iMsg = 1,gf%nmsgs
        if (error) cycle

        grib_raw => gf%grib_raw(grib_offsets(iMsg)+1:grib_offsets(iMsg+1))

        call grib_new_from_message(gf%indGrib(iMsg), grib_raw, io_status)

        if (io_status /= GRIB_SUCCESS) then
          !$OMP CRITICAL (id_list_from_grib_file_error)
          call msg("Trouble parsing message", iMsg)
          call msg("grib_file: "//trim(gf%fname))
          call grib_report_error(io_status)
          call set_error("after grib_new_from_message", sub_name)
          !$OMP END CRITICAL  (id_list_from_grib_file_error)

          cycle
        endif

        !if (iMsg==33) call ooops("33")
        call gribheadings_to_field_id (gf%indGrib(iMsg), .true., gf%obstime_interval,&
             &  gf%id(iMsg), gf%scale_factor(iMsg), gf%ifZeroNegatives(iMsg))

!        call msg("Message from GRIB:", iMsg, gf%nmsgs)
!        call report(gf%id(iMsg))

      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      
    endif
    if (error) return

    gf%idList_ready=.True.
    idList => gf%id(1:gf%nmsgs)
    iCount = gf%nmsgs

  end subroutine id_list_from_grib_file

  function fu_grib_filename(grib_unit)

    implicit none

    integer, intent(in) :: grib_unit
    CHARACTER (LEN=fnlen) :: fu_grib_filename
    character (len=*), parameter :: sub_name="fu_grib_filename"
    
    if (grib_unit < 0 .or. grib_unit > max_nbr_of_iGRIB_files) then
        call set_error("bad grib_unit = "//trim(fu_str(grib_unit)), sub_name)
        return
    endif

    if ( .not.  iGRIB(grib_unit)%defined ) then
        call set_error("Not opened grib_unit = "//trim(fu_str(grib_unit)), sub_name)
        return
    endif
    fu_grib_filename = iGRIB(grib_unit)%fname

  end function fu_grib_filename

  !******************************************************************

  subroutine read_field_from_grib_file(grib_unit, field_id, iFldNo, grid_data)
    ! Reads the requested field from grib. Either Id of iFldNo shoudl 

   implicit none
    ! Imported parameters
    integer, intent(in) :: grib_unit
    TYPE(silja_field_id), INTENT(in) :: field_id ! Field id 
    integer, INTENT(in) :: iFldNo ! Hint on field_index
    real, dimension(:), intent(out) :: grid_data  ! must exist

    type(silja_field_id),dimension(:), pointer :: idList  !helper pointer
    integer :: iTmp, fs

    type(grib_input_file), pointer :: gf

    character (len=*), parameter :: sub_name="read_field_from_grib_file"

    if (grib_unit < 0 .or. grib_unit > max_nbr_of_iGRIB_files) then
        call set_error("bad grib_unit = "//trim(fu_str(grib_unit)), sub_name)
        return
    endif

    gf => iGRIB(grib_unit)

    if (gf%idList_ready) then
      idList => gf%id(1:gf%nmsgs)
    else
         call id_list_from_grib_file(grib_unit, idList, iTmp)
         if (error) return
    endif

    if (.not. (idList(iFldNo) == field_id)) then
      call msg("Requested ID:")
      call report(field_id)
      call msg("ID No "//trim(fu_str(iFldNo))//":")
      call report(idList(iFldNo))
      call set_error("ID number does not match requested ID", sub_name)
    endif

    fs = fu_number_of_gridpoints(fu_grid(field_id))
    if (gf%who_has(iFldNo) == int_missing) then
      !good old parsing...
      call get_data_from_grib_index(gf%indGrib(iFldNo), field_id, &
        & gf%scale_factor(iFldNo), gf%ifZeroNegatives(iFldNo), grid_data, .True.)
    else
      call set_error("Field sharing is not implemented yet", sub_name)
!      if (smpi_adv_rank == gf%who_has(iFldNo)) then
!         call smpi_bcast_aray( gf%vals(iFldNo)%vals, fs, smpi_adv_comm, smpi_adv_rank)
!         grid_data(1:fs) = gf%vals(iFldNo)%vals(1:fs)
!      else
!         call smpi_bcast_aray( grid_data, fs, smpi_adv_comm, gf%who_has(iFldNo))
!      endif
    endif

  end subroutine read_field_from_grib_file

  !******************************************************************

  logical function fu_asimof_field(indGrib)
    !
    ! Checks if the given field is an Asimof field. 
    ! ALL criteria must be met:
    ! 1. The field must be small - less than 1000 octets. Note that this option causes problems in Windows
    ! 2. -----
    ! 3. Grid must be 255
    ! 4. Date must be "strange" -  20190101 01:01
    !
    implicit none

    ! Imported parameter
    INTEGER, intent(in) ::  indGrib

    ! Local variables
    integer :: iTmp, io_status
    character(len=20) :: chTmp
    !
    ! Problem: c_size_t in Windows has funny values, plus stack frame check fails. So,
    ! the message size is checked only in Linux. Other checkes are fine.
    integer(kind=4) :: msg_size

    fu_asimof_field = .false.

#ifndef VS2012
    call grib_get_message_size(indGrib, msg_size, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status, &
                                                         & 'call grib_get_message_size',&
                                                         & 'fu_asimof_field')
    if(msg_size > 1000) return                      ! Length of the message

#endif

    ! check for paramId used to be here.  Removed.
    ! It was triggering fu_asimof_field for parameters unknown to grib_api
    ! (might be known to silam). 

    call grib_get(indGrib,'typeOfGrid',iTmp, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'typeOfGrid','fu_asimof_field')
    if(iTmp /= 255) return                      ! Grid definition
    
    
    call grib_get(indGrib,'dataDate',iTmp, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'dataDate','fu_asimof_field')
    if(iTmp /= 20190101) return                     ! 2019-01-01

    call grib_get(indGrib,'dataTime',iTmp, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'dataTime','fu_asimof_field')
    if(iTmp /= 101) return                        ! Hour=1, min=1

    fu_asimof_field = .true.

  end function fu_asimof_field



  ! ***************************************************************

  SUBROUTINE read_next_header_from_gribfile(grib_unit, wdr,&
                                          & fix_oddities, &
                                          & check_asimof, &
                                          & indGrib, &   ! GRIB index of the field found
                                          & id,&         ! ID fo the field found
                                          & scale_factor, & ! factor to SI units
                                          & eof, ifZeroNegatives)         ! whether file is over
    !
    ! Returns the next field from an open grib-file, to which
    ! a unit-number grib_unit points. At this phase the field is separated
    ! into two parts: the identification (silja_field_id) and grid-data
    ! (real array).
    !
    ! Method:
    ! gridrecord_to_field to decode and find administrative data
    ! fix_oddities (module NWPM_ADMINISTRATIONS) 
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    INTEGER, INTENT(in) :: grib_unit ! from which to read
    type(silja_wdr),  INTENT(in)   :: wdr
    LOGICAL, INTENT(in) :: fix_oddities ! if true, then this module tries to get rid of
                                        ! strange field definitions, or "grib-dialects"
    LOGICAL, INTENT(in) :: check_asimof ! if true, then this module tries
                                        ! to check if the field is a header of an Asimof GRIB file
    
    ! Imported parameters with intent OUT:
    TYPE(silja_field_id), INTENT(out) :: id
    integer, intent(out) :: indGrib
    real, intent(out) :: scale_factor
    LOGICAL, INTENT(out) :: eof
    logical,  intent(out) :: ifZeroNegatives

    ! Local declarations:
    LOGICAL :: ifAsimof
    INTEGER :: io_status

    ifAsimof = .true.  ! For a while

    do while(ifAsimof)
      ! 
      ! Read from file
      !
      !if (grib_unit > 0) then
        call grib_new_from_file(grib_unit, indGrib, io_status)
      !else
      !  if (imessage <= nMessages) then
      !    call grib_new_from_message(indGrib, &
      !                  & grib_raw(grib_offsets(iMessage)+1:grib_offsets(iMessage+1)), &
      !                  & io_status)
      !  else
      !    io_status = GRIB_END_OF_FILE
      !  endif
      !endif


      IF (io_status /= GRIB_SUCCESS) THEN
        IF (io_status == GRIB_END_OF_FILE) THEN
        !  nMessages = int_missing
        !  imessage = 0
          eof=.true.
        ELSE
          call msg('Problem with reading GRIB file, io_status:',io_status)
          call set_error('Problem with reading GRIB file', 'read_next_header_from_gribfile')
        END IF
        RETURN
      ELSE
 !       imessage = imessage + 1 !Got this one, switch to the next one!
        eof=.false.
      END IF
      !
      ! Check if this field is an Asimof field
      !
!      if(imessage == 0)then
        ifAsimof = fu_asimof_field(indGrib)

!      else
!        ifAsimof = .false. 
!      endif

      if(ifAsimof)then
        call msg('Seems to be an Asimof field')
        call release_grib_buffer(indGrib)
      endif

    end do  ! while Asimof

    !
    ! Field is reasonable => work it out!
    !
    call gribheadings_to_field_id (indGrib, .true., fu_obstime_interval(wdr), id, scale_factor, ifZeroNegatives)
    !
    ! If an error occurs, print the contents of headings, otherwise
    !    decode the GRIB data into the data array
    !    GRIB data can be in a funny scan mode => someone has to re-arrange 
    !    them to SILAM internal representation
    !
    IF (error) THEN
      if (debug_level > 1) CALL print_gribheadings(indGrib)
      call set_missing(id)
      call unset_error('read_next_header_from_gribfile')
    endif
    
    if(defined(id))then
      if(fu_number_of_gridpoints(fu_grid(id)) <2)then
        call set_error('Strange number of grid points','read_next_header_from_gribfile')
        return
      end if
    END IF

  END SUBROUTINE read_next_header_from_gribfile


  !*********************************************************************
  
  subroutine get_data_from_grib_index(indGribField, id, scale_factor, ifZeroNegatives, grid_data, ifFixOddities)
    !
    ! Obtains the data from the GRIB message whose index is given.
    ! Takes care of the scanning mode in the grib file.
    ! Sometimes, "empty" values are replaced with GRIB-missing ones. Rubbish but
    ! have to keep in mind this possibility
    !
    implicit none
    
    ! Imported parameters
    integer, intent(in) :: indGribField
    type(silja_field_id), intent(in) :: id
    real, intent(in) ::  scale_factor
    real, dimension(:), intent(out), target :: grid_data
    logical, intent(in) :: ifZeroNegatives, ifFixOddities

    ! Local variables
    integer :: nValues, ioStat, i,j, index_1D, nx,ny, scan_mode
    real, dimension(:), pointer :: gribrecord
    real :: grib_real_missing, real_missing_replacement, grib_precision
    !
    ! Chack for stupidity
    !
    call grib_get_size(indGribField,'values',nValues, ioStat)
    if(ioStat /= GRIB_SUCCESS)then
      call process_grib_status(ioStat,'get_size: values','get_data_from_grib_index')
      return
    endif
    if(size(grid_data) < nValues)then
      call msg('Too small data array:', size(grid_data), nValues)
      call set_error('Too small data array','get_data_from_grib_index')
      return
    endif
    !
    ! OK, can get the data. However, have to take into account the scanning mode, which
    ! is described by three upper bits of a byte in section 2. Get them, compile the 
    ! old-fashipned scanning mode and decode the array.
    ! Note the stupidity: iScansNegatively means that rows enumerated descendingly, i.e. 
    !                     the second array index is descreasing.
    !                     Same for jSacnsPositively: columns are enumerated ascendingly, i.e.
    !                     the first index is increasing
    !                     Same for jPointsAreConsequtive: the first index should be scanned first
    ! This is probably because of lat-lon usual flip.
    !
    call grib_get(indGribField,'iScansNegatively',i,ioStat) ! 0 if rows enumerated DEscendingly
    if(ioStat /= GRIB_SUCCESS) then
      call process_grib_status(ioStat,'iScansNegatively','get_data_from_grib_index')
      return
    endif
    
    call grib_get(indGribField,'jScansPositively',j,ioStat)  ! 0 if columns enumerated Ascendingly
    if(ioStat /= GRIB_SUCCESS) then
      call process_grib_status(ioStat,'jScansPositively','get_data_from_grib_index')
      return
    endif
    call grib_get(indGribField,'jPointsAreConsecutive',index_1D,ioStat) ! 0 if consequtive values along row
    if(ioStat /= GRIB_SUCCESS) then
      call process_grib_status(ioStat,'jPointsAreConsecutive','get_data_from_grib_index')
      return
    endif
    !
    ! Re-arrange gribrecord to the grid_data array. Needed
    ! for the cases of some funny scan modes are used. Scan mode means 
    ! the rule how 2D (i,j) array is mapped to 1D gribrecord
    ! Options are:
    ! the first index - column-index grows / drops:  i++/--
    ! the second index - row-index grows / drops:  j++/--
    ! consequtive elements are i and i+1 / j and j+1: row/column-consequtive, respectively
    !
    scan_mode = i * 128 + j * 64 + index_1D * 32
    if(ifFixOddities)call fix_scan_mode(fu_met_src(id), scan_mode)
    !call msg('grib scan mode: ', scan_mode)
    if(scan_mode == 64)then
      gribrecord => grid_data
    else
      gribrecord => fu_work_array(nValues)
      if(error)return
    endif



    call grib_get(indGribField, 'values', gribrecord, ioStat)
    if(ioStat /= GRIB_SUCCESS)then
      call process_grib_status(ioStat,'values','get_data_from_grib_index')
      return
    endif
    !
    ! If there is something funny, take care of it here
    !
    if(ifFixOddities)then
      if(fu_if_replace_real_missing(fu_met_src(id)))then
        real_missing_replacement = fu_real_missing_replacement(fu_quantity(id))
        call grib_get(indGribField,'missingValue',grib_real_missing,ioStat)
        if(ioStat /= GRIB_SUCCESS)then
          call process_grib_status(ioStat,'missingValue','get_data_from_grib_index')
          return
        endif
        do i = 1, nValues
          if(gribrecord(i) == grib_real_missing) gribrecord(i) = real_missing_replacement
        end do
      endif
      
      if (fu_if_fix_packing(fu_met_src(id), fu_quantity(id))) then
        call grib_get(indGribField,'packingType',i,ioStat)
        if(ioStat /= GRIB_SUCCESS)then
          call process_grib_status(ioStat,'packingType','get_data_from_grib_index')
          return
        endif
        if (i == 0) then !Simple packing
          call grib_get(indGribField,'binaryScaleFactor',i,ioStat)
          if(ioStat /= GRIB_SUCCESS)then
            call process_grib_status(ioStat,'packingType','get_data_from_grib_index')
            return
          endif
          grib_precision = 2 * (2.0 ** i) ! Twice the discrete of the GRIB encoding
          where (abs(gribrecord(1:nValues)) < grib_precision) gribrecord(1:nValues) = 0.
        endif
      endif
    endif

    if(scan_mode == 64)then  ! No work arrays used

      if(abs(scale_factor - 1.0) > 1e-6 )then
        grid_data(1:nValues) = grid_data(1:nValues) * scale_factor
      endif
      
    else
      !
      ! Have to re-arrange the data
      !
      
      call grid_dimensions(fu_grid(id),nx,ny)
      if(error) return

      select case(scan_mode) ! bits 1,2,3 of the vaule define. Other bits=0
        case (0)   !----------- 000 column++, row--, row-consesqutive
          do j=1,ny
            do i=1,nx
              index_1D = i + (ny-j)*nx
              grid_data(i+(j-1)*nx) = gribrecord(index_1D) * scale_factor
            end do
          end do

        case (32)  !----------- 001 column++, row--, column-consequtive
          do j=1,ny
            do i=1,nx
              index_1D = ny-j+1 + (i-1)*ny
              grid_data(i+(j-1)*nx) = gribrecord(index_1D) * scale_factor
            end do
          end do

! Normal scanning mode is treated above
!        case (64)  !------------ 010 column++, row++, row-consequtive "NORMAL" FORTRAN array mapping
!          grid_data(1:min(SIZE(grid_data),SIZE(gribrecord))) = &
!                                              & gribrecord(1:min(SIZE(grid_data),SIZE(gribrecord)))

        case (96)  !------------ 011 column++, row++, column-consequtive
          do j=1,ny
            do i=1,nx
              index_1D = j + (i-1)*ny
              grid_data(i+(j-1)*nx) = gribrecord(index_1D) * scale_factor
            end do
          end do

        case (128) !------------ 100 column--, row--, row-consequtive
          do j=1,ny
            do i=1,nx
              index_1D = nx-i+1 + (ny-j)*nx
              grid_data(i+(j-1)*nx) = gribrecord(index_1D) * scale_factor
            end do
          end do

        case (160) !------------ 101 column--, row--, column-consequtive
          do j=1,ny
            do i=1,nx
              index_1D = ny-j+1 + (nx-i)*ny
              grid_data(i+(j-1)*nx) = gribrecord(index_1D) * scale_factor
            end do
          end do

        case (196) !------------ 110 column--, row++, row-consequtive
          do j=1,ny
            do i=1,nx
              index_1D = nx-i+1 + (j-1)*nx
              grid_data(i+(j-1)*nx) = gribrecord(index_1D) * scale_factor
            end do
          end do

        case (224) !------------ 111 column--, row++, column-consequtive
          do j=1,ny
            do i=1,nx
              index_1D = ny-j+1 + (i-1)*ny
              grid_data(i+(j-1)*nx) = gribrecord(index_1D) * scale_factor
            end do
          end do
        case default
          call set_error('Unknown scanning mode','get_data_from_grib_index')
          return
      end select

      call free_work_array(gribrecord)
  
    endif  ! whether the datarearrangement is needed

    ! Force non-negative if needed
    if (ifZeroNegatives) where (grid_data(1:nValues) < 0.) grid_data(1:nValues) = 0
  
  end subroutine get_data_from_grib_index


  !*************************************************************
  
  subroutine release_grib_buffer(indexToRelease)
    !
    ! A feature of the split data acquisition: the buffers are never reused unless
    ! explicitly released. This function does exactly that - releases the buffer
    !
    implicit none
    
    ! Imported parameter
    integer, intent(in) :: indexToRelease
    
    ! Local variables
    integer :: iRet
    
    call grib_release(indexToRelease, iRet)
    IF (iRet /= GRIB_SUCCESS) THEN
      call set_error('Cannot release the grib field','release_grib_buffer')
      CALL unset_error('release_grib_buffer')
      RETURN
    END IF

  end subroutine release_grib_buffer


  !*************************************************************

  subroutine grib_report_error(iErr)
    !Just reports GRIB erro messge
    integer, intent(in) :: iErr
    integer :: stat
    character (len=fnlen) ::  errmsg

    call grib_get_error_string(iErr,errmsg,stat)
    IF (stat == GRIB_SUCCESS) THEN
      call msg("GRIB_ERROR "//trim(fu_str(iErr))//":"//trim(errmsg))
    else
      call msg("Could not get GRIB error for code", stat )
    endif
  end subroutine grib_report_error


  SUBROUTINE close_grib_file_i(grib_unit)

    ! Frees gribfile structure

    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: grib_unit

    !local
    integer :: iErr, iMsg
    type(grib_input_file), pointer :: gf
    character (len=*), parameter :: sub_name="close_grib_file_i"

    if (grib_unit < 0 .or. grib_unit > max_nbr_of_iGRIB_files) then
        call set_error("bad grib_unit = "//trim(fu_str(grib_unit)), sub_name)
        return
    endif
    gf => iGRIB(grib_unit)
    if ( .not. gf%defined ) then
        call set_error("Not opened grib_unit = "//trim(fu_str(grib_unit)), sub_name)
        return
    endif

    do iMsg=1,gf%nmsgs
      if(gf%indGrib(iMsg) == int_missing) cycle
      call grib_release(gf%indGrib(iMsg), iErr)
      IF (iErr /= GRIB_SUCCESS) THEN
        call msg("Trouble closing GRIB file: "//trim(gf%fname))
        call grib_report_error(iErr) 
        call msg("Failed message ", iMsg, gf%nMsgs)
        if (gf%idList_ready) then  
          call msg("ID of failed message:")
          call report(gf%id(iMsg))
        else
          call msg("No ID generated")
        endif
        call msg("gf%indGrib(iMsg)",gf%indGrib(iMsg))
        call set_error('Cannot release the grib field','close_grib_file_i')
        RETURN
      END IF
      gf%indGrib(iMsg) = int_missing

      !free storage
!      if (gf%who_has(iMsg) == smpi_adv_rank) then
!        deallocate(gf%vals(imsg)%vals)
!        gf%who_has(iMsg) = int_missing
!      endif
    enddo

    gf%defined = .False.
    gf%fname = ""
    gf%idList_ready = .false.
    gf%nmsgs=0
    gf%obstime_interval = interval_missing
    
  END SUBROUTINE close_grib_file_i

  ! ***************************************************************
  !
  !       Private functions and subroutines (reading)
  !
  ! ***************************************************************

  subroutine gribheadings_to_field_id (indGrib, fix_oddities, obstime_interval, id, scale_factor, ifZeroNegatives)
    ! 
    ! Creates a complete field-identification from the grib-headings.
    !
    ! GRIB-standard: 
    ! WMO FM-92-IX Ext. GRIB, 1988 edition.
    ! Library: GRIB_API v.1.8
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent OUT
    TYPE(silja_field_id), intent(out) :: id

    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: indGrib
    LOGICAL, INTENT(in) :: fix_oddities ! if true, then some strange nwp-dependent 
                                        ! field-definiots are fixed. If false then the
                                        ! fields are returned just as they are.
    TYPE(silja_interval), intent(in) :: obstime_interval
    real, intent(out) :: scale_factor
    logical,  intent(out) :: ifZeroNegatives

    ! Local declarations:
    INTEGER :: grib_edition, quantity, field_kind, iTMp, jTmp, kTmp, io_status
    TYPE(silja_grid) :: grid
    TYPE(silja_level) :: level
    TYPE(silja_time) :: analysis_time
    TYPE(silja_interval) :: forecast_length, length_of_accumulation
    type(meteo_data_source) :: source_of_data
    character(len=clen) :: chSpeciesString
    character(len=clen) :: cfName, shortName
    real :: fTmp
    type(Taerosol_mode) :: aerosol_mode
!DEBUG only
    integer :: discipline, paramCategory, paramNbr
    character (len=*), parameter :: sub_name = "gribheadings_to_field_id"

    call set_missing(id)

    !
    ! Check the edition number and call for the corresponding processing
    !
    !
    ! Find the parts of identification.
    ! 
    call grib_get(indGrib,'centre',iTmp, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'centre', sub_name)
    if(error)return
    call grib_get(indGrib,'subCentre',kTmp, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'subCentre', sub_name)
    if(error)return
    call grib_get(indGrib,'generatingProcessIdentifier',jTmp, io_status)
!    call grib_get(indGrib,'subCentre',jTmp, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'generatingProcessIdentifier', sub_name)
    if(error)return

    ! This might be needed for fix_strange_quantities
    source_of_data = fu_set_met_src(int_missing, & ! SrcIdx, 
                                  & iTmp, &    ! centre, 
                                  & kTmp, &    ! subCentre
                                  & jTmp, &    ! model, 
                                  & '') !      ! name


       call grib_get(indGrib,'editionNumber', grib_edition, io_status)
       if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'editionNumber', sub_name)

    if (grib_name_table_enabled) then

       !
       !Attempt to make edition-agnostic intepretation
       call grib_get(indGrib,'cfName', cfName, io_status)
       if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'editionNumber', sub_name)

       call grib_get(indGrib,'shortName', shortName, io_status)
       if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'editionNumber', sub_name)
!
!      if (grib_edition == 2) then
!              call grib_get(indGrib,'discipline', discipline, io_status)
!              if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfParameter','process_grib_2_header')
!              if(error)return
!              call grib_get(indGrib,'parameterCategory', paramCategory, io_status)
!              if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfParameter','process_grib_2_header')
!              if(error)return
!              call grib_get(indGrib,'parameterNumber', paramNbr, io_status)
!              if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfParameter','process_grib_2_header')
!              if(error)return
!
!              call msg("shortName:"//trim(shortName))
!              call msg("cfName:"//trim(cfName))
!              if (discipline == 2 .and. paramCategory == 0 .and. paramNbr == 0) call ooops("Landfr")
!       endif



       call get_silam_params_for_grib(cfName, shortName, quantity, chSpeciesString, level, scale_factor)

       if(error) then
         call set_error("After get_silam_params_for_grib", sub_name)
         return
       endif
    else
       scale_factor = 1.0  !!default to feed to fix_strange_quantities
       select case(grib_edition)
         case (1)
           call process_grib_1_header(indGrib, source_of_data, quantity, chSpeciesString, level)
         case (2)
           call process_grib_2_header(indGrib, quantity, chSpeciesString, level)
         case default
           CALL msg('Unknown grib-edition:', grib_edition)
           CALL set_error('Unknown grib-edition', sub_name)
       end select
       if(error) then
         call set_error("After process_grib_X_header", sub_name)
         return
       endif
    endif   

    if ( quantity < 1) return  !!No quantity found  

    !
    ! Processing is over, the main variables are set. Now the SILAM part of the story 
    !
    ! Find the times of grib. Here there is an important point. In GRIB the time
    ! component is defined a bit clumsy, so that SILAM uses another one. This 
    ! very function transforms the GRIB definition to the SILAM one. It may seem 
    ! to be also clumsy, but there is a constraint - many fields are produced
    ! inside SILAM via fu_set_field_id, so that part has to be left untouched.
    !
    CALL get_grib_times(indGrib, &
                      & field_kind, &  ! Already SILAM-definition
                      & analysis_time,&
                      & forecast_length,&
                      & length_of_accumulation)
    if (field_kind == int_missing) return
    IF (error) RETURN 

    !
    ! Find the grid of grib
    !
    grid = fu_grib_grid(indGrib)
    IF (error) RETURN    

    !
    ! Find the level of grib unless it is already set by the GRIB code table above.
    !
    if( .not. defined (level) )then
        level = fu_grib_level(indGrib, source_of_data)
        if( (.not. defined(level)) .or. error) return
    endif
    !
    ! Source of data - reset if something is known. This assignment sets the string name in the data source
    ! Comparison does not check the name of the MDS whereas assignment sets it.
    !
    if(source_of_data == fmi_hirlam_src) then
       source_of_data = fmi_hirlam_src
    elseif(source_of_data == fmi_silam_src) then
      source_of_data = fmi_silam_src
    elseif(source_of_data == silam_internal_src) then
      source_of_data = silam_internal_src
    elseif(source_of_data == centre_model_ecmwf_src) then
      source_of_data = centre_model_ecmwf_src
    endif
    IF (error) RETURN
    !
    ! Depending on the issuing centre and the model, there may be quite a lot
    ! different oddities and specifics.
    !

if(ifDebug)then
 !  if(smpi_global_rank == 0) call grib_dump(indGrib, io_status)
  call msg('*** Before fixing the strange quantities************')
  call report(source_of_data)
  call msg('Level and grid for:' + fu_quantity_string(quantity), quantity)
  call report(level)
  call report(grid)
  call msg('Field_kind',field_kind)
  call msg('Analysis time:'+ fu_str(analysis_time) + &
         & ', accumulation length:' + fu_str(length_of_accumulation) + &
         & ', forecast length:' + fu_str(forecast_length))
endif
    
   

    CALL fix_strange_quantities( obstime_interval, &  !Might need to fix them differently 
                              & fu_centre(source_of_data), & ! centre
                              & fu_subCentre(source_of_data), & ! subCentre
                              & fu_model(source_of_data), & ! model / issuing process id
                              & quantity, level, &
                              & field_kind, &
                              & analysis_time, &
                              & forecast_length, &
                              & length_of_accumulation, &
                              & grid, &
                              & scale_factor, ifZeroNegatives)
if(ifDebug)then
  call msg('*** After fixing the strange quantities************')
  call report(source_of_data)
  call msg('Level and grid for:' + fu_quantity_string(quantity), quantity)
  call report(level)
  call report(grid)
  call msg('Field_kind , scale_factor',field_kind, scale_factor)
  call msg('Analysis time:'+ fu_str(analysis_time) + &
         & ', accumulation length:' + fu_str(length_of_accumulation) + &
         & ', forecast length:' + fu_str(forecast_length))
endif

    if(error)return
    !
    ! Finally set the identification.
    !
    if (field_kind == forecast_flag) length_of_accumulation = zero_interval

    id = fu_set_field_id (source_of_data,&
                            & quantity, &
                            & analysis_time,&
                            & forecast_length, &
                            & grid,&
                            & level,&
                            & length_of_accumulation, &
                            & zero_interval, &
                            & field_kind, &
                            & species = fu_species_from_short_string(chSpeciesString))

if(ifDebug)then
  call msg('============= Final ID: ========================')
  call report(id)
  call msg('')
endif

  END subroutine gribheadings_to_field_id


    !============================================================================

    subroutine process_grib_1_header(indGrib, source_of_data, quantity, chSpeciesString, level)
      !
      ! GRIB 1 processing using GRIB-API
      !
      implicit none

      integer, intent(in) :: indGrib
      type(meteo_data_source), intent(in) :: source_of_data
      integer, intent(out) :: quantity                   ! Returned SILAM quantity
      character(len=*), intent(out) :: chSpeciesString   ! Returned SILAM substance name
      type(silja_level), intent(out) :: level            ! Returned SILAM level

      integer :: iTmp, jTmp, kTmp, tabVersion, io_status
      real :: fTmp
      !
      ! Find the parts of identification.
      ! 
      !
      ! Find the quantity of grib. New function based on grib_code_table module
      ! It also returns the level if it is fully determined by the quantity
      ! Otherwise, level_missing returned
      !
      call grib_get(indGrib, 'table2Version', tabVersion, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'table2Version','process_grib_1_header')
      if(error)return
      call grib_get(indGrib,'indicatorOfParameter', iTmp, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfParameter','process_grib_1_header')
      if(error)return
      call grib_get(indGrib,'indicatorOfTypeOfLevel', jTmp, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfTypeOfLevel','process_grib_1_header')
      if(error)return
      call grib_get(indGrib,'level', fTmp, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'level','process_grib_1_header')
      if(error)return

!      call msg("fu_centre fu_subCentre fu_model", (/fu_centre(source_of_data),fu_subCentre(source_of_data),fu_model(source_of_data)/))
!      call msg("tabVersion, paramID, levType, LevVal", (/tabVersion,iTmp,jTmp,nint(fTmp)/))
      call get_silam_params_for_grib_1(tabVersion, &  ! GRIB 1, table version
                                     & fu_centre(source_of_data), &  ! issuing centre
                                     & fu_subCentre(source_of_data), &  ! issuing sub-centre
                                     & fu_model(source_of_data), &  ! originating model
                                     & iTmp, &  ! GRIB parameter code
                                     & jTmp, &  ! level type
                                     & fTmp, &  ! level 1 value
                                     & quantity, &  ! Returned SILAM quantity
                                     & chSpeciesString, &   ! Returned SILAM species IO string
                                     & level)       ! Returned SILAM level
      IF (error) THEN
        call msg('*************************************************')
        call msg('unknown combination of parameters::')
        call msg ('parameter', iTmp)
        call msg('levtyp, real(level value):', jTmp, fTmp)
        RETURN
      END IF
      if(quantity == int_missing)then
!        call msg('Unknown parameter value (GRIB 1):',iTmp)
!        if(smpi_global_rank == 0) call grib_dump(indGrib, io_status)
      endif

    end subroutine process_grib_1_header


    !============================================================================

    subroutine process_grib_2_header(indGrib, &
                                   & quantity, &  ! Returned SILAM quantity
                                   & chSpeciesString, &   ! Returned SILAM substance name
                                   & level)
      !
      ! GRIB 1 processing using GRIB-API
      !
      implicit none

      integer, intent(in) :: indGrib
      integer, intent(out) :: quantity                   ! Returned SILAM quantity
      character(len=*), intent(out) :: chSpeciesString   ! Returned SILAM substance name
      type(silja_level), intent(out) :: level            ! Returned SILAM level

      ! Local declarations
      integer :: discipline, paramCategory, paramNbr, iTmp, jTmp, kTmp, io_status, tabVersion
      real :: fTmp
      character(len=clen) :: shortNm
      character(len=fnlen) :: nm

      !
      ! Find the quantity of grib. New function based on grib_code_table module
      ! It also returns the level if it is fully determined by the quantity
      ! Otherwise, level_missing returned
      !
      call grib_get(indGrib,'discipline', discipline, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfParameter','process_grib_2_header')
      if(error)return
      call grib_get(indGrib,'parameterCategory', paramCategory, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfParameter','process_grib_2_header')
      if(error)return
      call grib_get(indGrib,'parameterNumber', paramNbr, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfParameter','process_grib_2_header')
      if(error)return

      !!!typeOfLevel and level are generated by GribAPI and might or might not coinside with 
      !!! typeOfFirstFixedSurface and scaledValueOfFirstFixedSurface
      !!! The latter two are parameters of a grib message
      call grib_get(indGrib,'typeOfFirstFixedSurface', jTmp, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'typeOfFirstFixedSurface','process_grib_2_header')
      if(error)return
      call grib_get(indGrib,'scaledValueOfFirstFixedSurface', fTmp, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'scaledValueOfFirstFixedSurface','process_grib_2_header')
      if(error)return
      call grib_get(indGrib,'name', nm, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'name','process_grib_2_header')
      if(error)return
      call grib_get(indGrib,'shortName', shortNm, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'shortName','process_grib_2_header')
      if(error)return
      call get_silam_params_for_grib_2(discipline, paramCategory, paramNbr, &
                                     & jTmp, &  ! level type
                                     & fTmp, &  ! level 1 value
                                     & quantity, &  ! Returned SILAM quantity
                                     & chSpeciesString, &   ! Returned SILAM substance name
                                     & level)       ! Returned SILAM level
      IF (error) THEN
        call msg('*************************************************')
        call msg('unknown combination of parameters:')
        call msg ('parameter', iTmp)
        call msg('levtyp, real(level value):', jTmp, fTmp)
        RETURN
      END IF
#ifdef DEBUG      
      if(quantity == int_missing)then
        call msg('unknown combination of (GRIB-2 discipline, parameterCategory, parameterNumber):', &
                                             & (/discipline, paramCategory, paramNbr/))
        call msg('"'//trim(nm)//'"  aka "'//trim(shortNm)//'"             levtyp, real(level value):', jTmp, fTmp)
!        if(smpi_global_rank == 0) call grib_dump(indGrib, io_status)
!        call msg('')
!        call msg('')
      endif
#endif      

    end subroutine process_grib_2_header




  ! ***************************************************************

  SUBROUTINE fix_strange_quantities(obstime_interval, centre, subCentre, model, quantity, &
                                  & level, &
                                  & field_kind, &
                                  & analysis_time, &
                                  & forecast_length, &
                                  & accumulation_length, &
                                  & grid, &
                                  & scale_factor_inout, ifZeroNegatives)
    !
    ! Fixes the strange parameters of the fields, but does NOT change quantities
    ! Reason: correct establishing of quantities is done via GRIB code table,
    ! so there is nothing to do here.
    ! What is handled here is, for example, wrong accumulation times
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silja_interval), intent(in) :: obstime_interval
    INTEGER, INTENT(in) :: centre, subCentre, model

    ! Imported parameters with intent INOUT:
    INTEGER, INTENT(inout) :: quantity, field_kind
    TYPE(silja_level), INTENT(inout) :: level
    type(silja_time), intent(inout) :: analysis_time
    type(silja_interval), intent(inout) :: forecast_length, accumulation_length
    type(silja_grid), intent(inout) :: grid

    ! Imported parameters with intent OUT
    real, intent(inout) :: scale_factor_inout ! OverrideIf units do not correspond to SILAM one
    logical, intent(out) :: ifZeroNegatives ! If units do not correspond to SILAM ones

    integer :: iTmp
    real :: scale_factor

    scale_factor = scale_factor_inout

    ifZeroNegatives = .false.

    if(.not.defined(level)) then
      !call set_error('undefined level given','fix_strange_quantities')
      level = fu_set_level(top_of_the_atmosphere)
      quantity = int_missing
      return
    end if

!    if (large_scale_accum_rain_flag == quantity) then
!            call msg ("ls_rain")
!    endif
!call msg('',quantity)
!if(quantity == int_missing)then
!call msg('*')
!endif


    SELECT CASE (centre)
      CASE (centre_moscow)
             select case(quantity)
             case(total_cloud_cover_flag)
               scale_factor = 0.01 !% to fraction
               
             case default ! Do nothing

             end select ! model

      CASE (centre_MEPS)

            if (model == model_meps) then
                  if (any(quantity == (/leaf_area_index_flag, &
                    & leaf_area_indexhv_flag, leaf_area_indexlv_flag/))) then !Kill missing lai
                      ifZeroNegatives = .True.
                  endif
            endif

      CASE (centre_hirlam, centre_fmi)

        select case(model)
          case (model_hirlam_eno, model_hirlam_ata, model_silam)

          case (model_harmonie, model_hirlam)

             !
             ! There is a hack in grib_code_table 
             ! Should remind about it, otherwise it will bae forgotten
             if (quantity == large_scale_accum_rain_flag .and. model == model_harmonie) then
                if (field_kind /= accumulated_flag) quantity = int_missing
                call msg_warning("Treating harmonie snow as ls_precip","fix_strange_quantities")
             endif
             if (quantity == convective_accum_rain_flag .and. model == model_harmonie) then
                if (field_kind /= accumulated_flag) quantity = int_missing
                call msg_warning("Treating harmonie rain as cnv_precip","fix_strange_quantities")
             endif
             !
             !  Set the accumulation_length for accumulated quantities
             !
             select case(quantity)
             case(large_scale_accum_rain_flag, &
                & convective_accum_rain_flag, &
                & total_precipitation_acc_flag, &
                & NWP_sensible_heatflux_ac_flag, &
                & NWP_latent_heatflux_ac_flag, &
                & surf_sw_down_radiation_ac_flag, &
                & surf_lw_down_radiation_ac_flag, &
                & surf_sw_net_radiation_ac_flag, &
                & surf_lw_net_radiation_ac_flag)
               
               field_kind = accumulated_flag
 !              call msg("Making it accumulated")
               if(.not.defined(accumulation_length))then
 !                call msg("Setting accumulation length")
                 if(defined(forecast_length))then
                   accumulation_length = forecast_length
                 else
                   call set_error('Undefined forecast len','fix_strange_quantities')
                   return
                 end if
               end if
               if(accumulation_length == zero_interval)then
                iTmp = fu_sec(obstime_interval)
                call msg("Changing accumulation_length", fu_sec(accumulation_length), iTmp)
                  accumulation_length = obstime_interval
               endif
             case (leaf_area_index_flag)
                 ifZeroNegatives = .True.

               
             case default ! Do nothing

             end select ! model

          case default
            call msg('ATTENTION, ATTENTION, ATTENTION')
            call msg('ATTENTION, ATTENTION, ATTENTION')
            call msg_warning('Unknown model','fix_strange_quantities')
            call msg('model ', model)
            call msg('ATTENTION, ATTENTION, ATTENTION')
            call msg('ATTENTION, ATTENTION, ATTENTION')
            return
        end select

        select case (fu_leveltype(level))
            ! the strange levels in hirlam md file with undefined hybrid coefficients have to be handled here
            case(hybrid)
                !if (fu_hybrid_level_number(level) == 1 .or.  fu_hybrid_level_number(level) == 60)then
                !    call msg('hybrid level nr', fu_hybrid_level_number(level))
                !    call msg('hybrid a, b', fu_hybrid_level_coeff_a(level), fu_hybrid_level_coeff_b(level))
                !endif
                if (fu_hybrid_level_number(level) > 1 .and. &
                &   (fu_hybrid_level_coeff_a(level) .eps. 0.) .and. & 
                &   (fu_hybrid_level_coeff_b(level) .eps. 0.)) then
                   
                   call msg_warning('Strange level found','fix_strange_quantities')
                   call msg('Strange level found for quantity:'+fu_quantity_string(quantity), quantity)
                   call report(level)
                   quantity = int_missing
                   level = fu_set_level(top_of_the_atmosphere)
                endif
            
            case(layer_btw_2_hybrid)

                if ((fu_hybrid_level_number(fu_lower_boundary_of_layer(level)) > 1 .and. &
                 &   (fu_hybrid_level_coeff_a(fu_lower_boundary_of_layer(level)) .eps. 0.) .and. & 
                 &   (fu_hybrid_level_coeff_b(fu_lower_boundary_of_layer(level)) .eps. 0.)) .or. &
                 &  (fu_hybrid_level_number(fu_upper_boundary_of_layer(level)) > 1 .and. &
                 &   (fu_hybrid_level_coeff_a(fu_upper_boundary_of_layer(level)) .eps. 0.) .and. & 
                 &   (fu_hybrid_level_coeff_b(fu_upper_boundary_of_layer(level)) .eps. 0.)))then
                   
                   call msg_warning('Strange level found','fix_strange_quantities')
                   call msg('Strange level found for quantity:'+fu_quantity_string(quantity), quantity)
                   call report(level)
                   quantity = int_missing
                   level = fu_set_level(top_of_the_atmosphere)
                endif
            case default
        end select

      case(centre_ecmwf)
        !
        !  Set the accumulation_length for accumulated quantities
        !
        select case(quantity)
          case (large_scale_accum_rain_flag, &
              & convective_accum_rain_flag, &
              & total_precipitation_acc_flag, &
              & NWP_sensible_heatflux_ac_flag, &
              & NWP_latent_heatflux_ac_flag, &
              & surf_sw_down_radiation_ac_flag, &
              & surf_lw_down_radiation_ac_flag, &
              & surf_sw_net_radiation_ac_flag, &
              & surf_lw_net_radiation_ac_flag)

            field_kind = accumulated_flag
            if(.not.defined(accumulation_length))then
              if(defined(forecast_length))then
                accumulation_length = forecast_length
              else
                call set_error('Undefined forecast len','fix_strange_quantities')
                return
              end if
            end if
          case (cloud_water_flag, cloud_ice_flag, cloud_cond_water_flag, specific_humidity_flag)
            ifZeroNegatives = .True. ! Some negatives come due to interpolation
                                     ! from harmonics

          case (u_flag, v_flag)
                 !Ignore surface winds     
             if (fu_leveltype(level) == constant_height ) then
               if(fu_level_height(level) == 0. ) then
                 quantity = int_missing
               endif
             endif

          case default ! Do nothing
        end select ! quantity

        !
        ! Correct scaling for precipitation
        !
        select case(quantity)
          case (large_scale_accum_rain_flag, &
              & convective_accum_rain_flag, &
              & total_precipitation_acc_flag)

            scale_factor = 1000. ! From [m] to [kg m**-2]=[mm]
          case default ! Do nothing
        end select ! quantity

      case (centre_smhi)   ! SMHI
        !
        ! SMHI centre with specific subcentre and model require proper precipitation accumulation to be set
        !
        if(subCentre == 98 .and. model == 1)then
           select case(quantity)
             case (large_scale_accum_rain_flag, convective_accum_rain_flag) 
        !            accumulation_length = one_hour * 3.0   ! it is per-time-step, which is 3hr in the climate model of SMHI

            !Dirty hack for ECHAM crap  (ensclim project)
            accumulation_length = one_hour * 6.0   ! Joana has got 6-hourly ECHAM files wth 3-hour accumulation time
            scale_factor = 2                       ! Double it

            field_kind = accumulated_flag

             case default ! Do nothing


          end select ! quantity


        endif

      case (centre_UERRA)
        if (any(quantity== (/cloud_cond_water_flag, &
              & specific_humidity_flag, &
              & large_scale_accum_rain_flag, &
              & convective_accum_rain_flag, &
              & total_precipitation_acc_flag &
              & /))  ) then
          ifZeroNegatives = .True.
        endif

        
      CASE default

    END SELECT ! centre selection

#ifdef DEBUG_GRIB    
    if (scale_factor_inout /= scale_factor) then
      call msg( "Chnaging scale factor for "//trim(fu_quantity_string(quantity))//" from,to", &
                & scale_factor_inout, scale_factor )
      call msg_warning("Overriding GRIB scale factor",  'fix_strange_quantities')
    endif
#endif
   scale_factor_inout = scale_factor


  END SUBROUTINE fix_strange_quantities


  !****************************************************************
  
  logical function fu_if_replace_real_missing(mds)
    !
    ! Returns true if the meteo data source is known or suspected in replacing meaningful "zeroes"
    ! with grib-real missing
    !
    implicit none
    
    ! Imported parameter
    type(meteo_data_source), intent(in) :: mds
    
    fu_if_replace_real_missing = .false.

    if(fu_centre(mds) == centre_smhi .and. fu_subCentre(mds) == 98 .and. fu_model(mds) == model_hirlam)& 
      & fu_if_replace_real_missing = .true.

    if(fu_centre(mds) == centre_MEPS .and. fu_model(mds) == model_meps)& 
      & fu_if_replace_real_missing = .true.
  end function fu_if_replace_real_missing

  !****************************************************************
  logical function fu_if_fix_packing(mds,quantity)
    !
    ! Returns true if the meteo data source is known or suspected in having
    ! packing, that causes negatives in strictly-positive quantities
    !
    implicit none
    
    ! Imported parameter
    type(meteo_data_source), intent(in) :: mds
    integer, intent(in) :: quantity
    
    fu_if_fix_packing = any(quantity==(/cloud_water_flag,cloud_ice_flag,specific_humidity_flag/))

  end function fu_if_fix_packing

  !****************************************************************
  
  subroutine fix_scan_mode(mds, scMod)
    !
    ! Scan mode for ECHAM data needs fixing
    !
    implicit none
    
    ! Imported parameter
    type(meteo_data_source), intent(in) :: mds
    integer, intent(inout):: scMod
    
    if(fu_centre(mds) == 98 .and. fu_subCentre(mds) == 232 .and. fu_model(mds) == 64)then
      scMod = 64
    endif
    
  end subroutine fix_scan_mode


  ! ***************************************************************

  SUBROUTINE get_grib_times(indGrib,&
                          & field_kind, &
                          & analysis_time,&
                          & forecast_length,&
                          & length_of_accumulation)
    ! 
    ! Translates GRIB time definition (field_kind + reference_time + P1 + P2)
    ! to SILAM one (analysis_time, forecast_length, accumulation_length, 
    ! valid_time_end). Some parts below look clumsy but...
    ! 
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    INTEGER, INTENT(in) :: indGrib

    ! Imported parameters with intent OUT:
    integer, intent(out) :: field_kind
    TYPE(silja_time), INTENT(out) :: analysis_time
    TYPE(silja_interval), INTENT(out) :: forecast_length
    TYPE(silja_interval), INTENT(out) :: length_of_accumulation ! the length of interval from 
                                                                ! valid time backwars that it took to
                                                                ! accumulate this field.
    ! Local declarations:
    INTEGER :: year, month, day, hour, minute, iStepUnit, iEndStep, iStartStep, io_status
    TYPE(silja_interval) :: start_of_accumulation
    character(len=20) :: chTmp
    integer :: iTmp
    !
    ! Basic reference time
    !
    call grib_get(indGrib,'dataDate',iTmp, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'dataDate','get_grib_times')

    year = iTmp / 10000
    IF (year < 1000) THEN
      IF (year > 50) THEN ! oh this is stupid!
        year = year + 1900
      ELSE
        year = year + 2000
      END IF
    END IF
    month = mod (iTmp / 100, 100 ) 
    day = mod (iTmp, 100 ) 


    call grib_get(indGrib,'dataTime',iTmp, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'dataTime','get_grib_times')
    
    hour = iTmp / 100 
    minute = mod (iTmp, 100) 


    length_of_accumulation = interval_missing

    !
    ! Here the GRIB reference_time is supposed to be the analysis time (seems to
    ! correspond to GRIB). Note that period-valid, average and difference products 
    ! use the same SILAM time components and differ only via field_kind
    !
    analysis_time = fu_set_time_utc(year,month,day,hour,minute,0.)

    ! Time range, accumulation and validity all depend on GRIB time range 
    ! indicator, p1 and P2
    !
    call grib_get(indGrib,'stepType',chTmp, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'stepType','get_grib_times')

    iStepUnit = 0  ! Force minutes, otherwise hours forced...
    call grib_set(indGrib,'stepUnits',iStepUnit, io_status)  !!Sic! SET it here!
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'set stepUnits','get_grib_times')

    call grib_get(indGrib,'endStep',iEndStep, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'endStep','get_grib_times')

    call grib_get(indGrib,'startStep',iStartStep, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'startStep','get_grib_times')

    forecast_length = fu_set_interval_min(iEndStep) ! Start of validity

    if (iStartStep==iEndStep) then ! valid for reference + P1
      !! normally  (chTmp == 'instant'), but can be something else 
        field_kind = forecast_flag
    elseif (chTmp == 'avg') then ! Average (or period-valid) field from ref+p1 till ref+p2, valid from r+P1 till ref+P2
        field_kind = averaged_flag
        length_of_accumulation = fu_set_interval_min(iEndStep-iStartStep)
    elseif (chTmp == 'accum') then ! accumulated from P1 to P2, valid for P2:
        field_kind = accumulated_flag
        length_of_accumulation = fu_set_interval_min(iEndStep-iStartStep)
    else
      field_kind = int_missing
      call msg('Non-supported stepType: '//trim(chTmp))
    endif

  END SUBROUTINE get_grib_times


  ! ***************************************************************

  FUNCTION fu_grib_grid(indGrib) RESULT(grid)
    !
    ! Finds the grid of grib-data.
    !
    ! Language: ANSI Fortran 90
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_grid) :: grid
    !
    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: indGrib

    ! Local declarations:
    TYPE(silam_pole) :: pole
    LOGICAL :: corner_in_geographical_latlon, i_grow, j_grow
    REAL :: grid_dist_x_deg, grid_dist_y_deg
    REAL :: corner_lat, corner_lon, end_corner_lon, end_corner_lat, pole_lon, pole_lat, fTmp
    real :: lat0, lat1, lat2, lon0  !!!Lambert parameters
    INTEGER :: scanning_mode, iTypeOfGrid, io_status, is_missing, nx, ny
    character (len=1024) :: strTmp ! 
    character (len=*), parameter :: sub_name="fu_grib_grid"
    !
    grid = grid_missing


    ! All (useful) grids have sizes, FirstGridPointInDegrees and gridType

    call grib_get(indGrib,'Ni',nx, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'Ni', &
                                                         & sub_name)
    call grib_get(indGrib,'Nj',ny, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'nJ', &
                                                         & sub_name)
                                                       
    call grib_get(indGrib,'longitudeOfFirstGridPointInDegrees',corner_lon, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'longitudeOfFirstGridPointInDegrees', &
                                                         & sub_name)

    call grib_get(indGrib,'latitudeOfFirstGridPointInDegrees',corner_lat, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'latitudeOfFirstGridPointInDegrees', &
                                                             & sub_name)
    
    call grib_get(indGrib,'gridType',strTmp, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'gridType',sub_name)

    if (error) then
      call set_error("Failed to get sizes/type for grid",sub_name)
      return
    endif

    select case(strTmp) !! gridType


      case("regular_ll")
      !
      ! 1. Set "normal" lat-lon grid (also called equidistant cylindrical)
      !    or Plate Carree projection grid. It can also be rotated if someone sets the pole
      !
      
        pole = pole_geographical

        corner_in_geographical_latlon = .True.

        ! Calculate dx and dy:
        ! A danger exists: direction increment (dx and dy) can have insufficient
        ! precision. Thus, ECMWF grid with 0.5625 degree resolution does not fit
        ! into the GRIB 3-digit accuracy. For small areas nobody cares but for 
        ! half-global or global ranges the error can be one grid cell: 
        ! 640 for 0.5625 and 639.431616 for 0.563 -- ignore them from GRIB!!!
        !
        call grib_get(indGrib,'longitudeOfLastGridPointInDegrees',end_corner_lon, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'longitudeOfLastGridPointInDegrees',&
                                                             & sub_name)

        call grib_get(indGrib,'latitudeOfLastGridPointInDegrees',end_corner_lat, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'latitudeOfLastGridPointInDegrees',&
                                                             & sub_name)

        ! A simple check: if data order is not "standard", lat start and end can be swapped
        ! Also, 360 degrees rotation may be needed
        if(corner_lon >= 180.0)corner_lon = corner_lon - 360.0 !! ECMWF global GRIB starts form 180
        if(end_corner_lon > 180.0)end_corner_lon = end_corner_lon - 360.0
        if(corner_lon < -180.0)corner_lon = corner_lon + 360.0
        if(end_corner_lon < -180.0)end_corner_lon = end_corner_lon + 360.0
        if(corner_lon >= end_corner_lon) end_corner_lon = end_corner_lon + 360.
        if(corner_lat > end_corner_lat)then  ! If south to north
          fTmp = corner_lat
          corner_lat = end_corner_lat
          end_corner_lat = fTmp
        endif

        is_missing=0
        call grib_is_missing(indGrib,'DiInDegrees', is_missing, io_status)
!        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'DiInDegrees', sub_name)

        IF (is_missing == 1) THEN  
          !
          ! Direction increments are not given, have to compute them explicitly
          !
          grid_dist_x_deg = (end_corner_lon-corner_lon)/REAL(nx-1)
          grid_dist_y_deg = (end_corner_lat-corner_lat)/REAL(ny-1)
        else
          !
          ! Check if compatible with corners
          !
          call grib_get(indGrib,'DiInDegrees',grid_dist_x_deg, io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'DiInDegrees',sub_name)

          call grib_get(indGrib,'DjInDegrees',grid_dist_y_deg, io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'DjInDegrees',sub_name)

          fTmp = end_corner_lon-corner_lon
!          if(fTmp < 0.0)fTmp = fTmp + 360.0
          if(abs(1.-fTmp/REAL(nx-1)/grid_dist_x_deg) < 0.05)then ! Numerical difference allowed in grid size
            grid_dist_x_deg = fTmp/REAL(nx-1)
          else
            !Try swapping the ends in lon direction (highly unlikely to be the correct thing to do ...)    
            call msg_warning('Swapping grid lon corners',sub_name)  
            fTmp = corner_lon
            corner_lon = end_corner_lon
            end_corner_lon = fTmp
            fTmp = end_corner_lon-corner_lon
            if(fTmp < 0.0)fTmp = fTmp + 360.0
            if(abs(1.-fTmp/REAL(nx-1)/grid_dist_x_deg) < 0.05)then
              grid_dist_x_deg = fTmp/REAL(nx-1)
            else
              call msg("Trouble with GRIB message")
              if(smpi_global_rank == 0) call grib_dump(indGrib, io_status)
              call set_error('Incompatible grid features1', sub_name)
              return
            endif
          endif

          if(abs(1.-(end_corner_lat-corner_lat)/REAL(ny-1)/grid_dist_y_deg) < 0.05) then
            grid_dist_y_deg = (end_corner_lat-corner_lat)/REAL(ny-1)
          else
            ! Either something strange in file or grid crossing the pole
              call msg("Trouble with GRIB message")
              if(smpi_global_rank == 0) call grib_dump(indGrib, io_status)
              call set_error('Incompatible grid features2', sub_name)
            return
          endif
        END IF  ! whether the increment is given

        grid = fu_set_lonlat_grid ('latlon_grid_from_grib', &
                                 & corner_lon, corner_lat, &
                                 & corner_in_geographical_latlon, &
                                 & nx, ny, &
                                 & pole, & 
                                 & grid_dist_x_deg, grid_dist_y_deg)

      !----------------------------------------
      !
      ! 2. Set modified (rotated) lat-lon grid.
      !
      case("rotated_ll")

        corner_in_geographical_latlon = .FALSE.

        !
        ! Depending on scanning mode, the SILAM corner can differ from
        ! the GRIB-claimed one
        !


        call grib_get(indGrib,'longitudeOfLastGridPointInDegrees',end_corner_lon, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'longitudeOfLastGridPointInDegrees', &
                                                             & sub_name)

        call grib_get(indGrib,'latitudeOfLastGridPointInDegrees',end_corner_lat, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'latitudeOfLastGridPointInDegrees', &
                                                             & sub_name)

        ! A simple check: if data order is not "standard", the start and end can be swapped
        ! Also, 360 degrees rotation may be needed
        if(corner_lon > 180.0)corner_lon = corner_lon - 360.0
        if(end_corner_lon > 180.0)end_corner_lon = end_corner_lon - 360.0
        if(corner_lon < -180.0)corner_lon = corner_lon + 360.0
        if(end_corner_lon < -180.0)end_corner_lon = end_corner_lon + 360.0
        if(corner_lat > end_corner_lat)then  ! If south to north
          fTmp = corner_lat
          corner_lat = end_corner_lat
          end_corner_lat = fTmp
        endif
        
        call grib_get(indGrib,'longitudeOfSouthernPoleInDegrees',pole_lon, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'longitudeOfSouthernPoleInDegrees', &
                                                             & sub_name)

        call grib_get(indGrib,'latitudeOfSouthernPoleInDegrees',pole_lat, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'latitudeOfSouthernPoleInDegrees', &
                                                             & sub_name)

          ! it is southpole given in grib
        pole = fu_set_pole(south_flag,  pole_lat,  pole_lon )
        IF (error)RETURN

        is_missing=0
        call grib_is_missing(indGrib,'DiInDegrees', is_missing)
        IF (is_missing == 1) THEN  
          !
          ! Direction increments are not given, have to compute them explicitly
          !
          grid_dist_x_deg = abs(corner_lon-end_corner_lon)/REAL(nx-1)
          grid_dist_y_deg = abs(corner_lat-end_corner_lat)/REAL(ny-1)
        else
          !
          ! Check if compatible with corners
          !
          call grib_get(indGrib,'DiInDegrees',grid_dist_x_deg, io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'DiInDegrees',sub_name)

          call grib_get(indGrib,'DjInDegrees',grid_dist_y_deg, io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'DjInDegrees',sub_name)

          fTmp = end_corner_lon-corner_lon
          if(fTmp < 0.0)fTmp = fTmp + 360.0
          if(abs(1.-fTmp/REAL(nx-1)/grid_dist_x_deg) < 0.05)then ! Numerical difference allowed in grid size
            grid_dist_x_deg = fTmp/REAL(nx-1)
          else
            !Try swapping the ends in lon direction (highly unlikely to be the correct thing to do ...)    
            call msg_warning('Swapping grid lon corners',sub_name)  
            fTmp = corner_lon
            corner_lon = end_corner_lon
            end_corner_lon = fTmp
            fTmp = end_corner_lon-corner_lon
            if(fTmp < 0.0)fTmp = fTmp + 360.0
            if(abs(1.-fTmp/REAL(nx-1)/grid_dist_x_deg) < 0.05)then
              grid_dist_x_deg = fTmp/REAL(nx-1)
            else
              call msg("Trouble with GRIB message")
              if(smpi_global_rank == 0) call grib_dump(indGrib, io_status)
              call set_error('Incompatible grid features3', sub_name)
              return
            endif
          endif

          if(abs(1.-(end_corner_lat-corner_lat)/REAL(ny-1)/grid_dist_y_deg) < 0.05) then
            grid_dist_y_deg = (end_corner_lat-corner_lat)/REAL(ny-1)
          else
            ! Either something strange in file or grid crossing the pole
              call msg("Trouble with GRIB message")
              if(smpi_global_rank == 0) call grib_dump(indGrib, io_status)
              call set_error('Incompatible grid features4', sub_name)
            return
          endif
        END IF ! whether the increment is given


      !     CALL report(grid)

      case("lambert") !!Dirty hack:
          !Pretend it is  rotated lon-lat grid with equator at lat0
          ! and 0-th meridian along LoVInDegrees

        call grib_get(indGrib,'LoVInDegrees',lon0, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'LoVInDegrees', &
                                                             & sub_name)
        call grib_get(indGrib,'LaDInDegrees',lat0, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'LaDInDegrees', &
                                                             & sub_name)
        call grib_get(indGrib,'Latin1InDegrees',lat1, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'Latin1InDegrees', &
                                                             & sub_name)
        call grib_get(indGrib,'Latin2InDegrees',lat2, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'Latin2InDegrees', &
                                                             & sub_name)
        if (lat0 /= lat1 .or. lat1/=lat2) then
            call msg("LoVInDegrees, LaDInDegrees, Latin1InDegrees, Latin2InDegrees", &
              & (/lon0,lat0, lat1, lat2/))
            call msg("Only case of all equal latitudes implemented so far")
            call set_error("Wrong (unimplemented) latitudes for Lambert conic grid", sub_name)
        endif          

        ! South pole of rll grid
        if  (lat0 > 0) then 
          pole_lat = lat0 - 90
          pole_lon = lon0
        else
          pole_lat = -90 - lat0
          pole_lon = mod(lon0, 360.) - 180
        endif 
                                                            
        pole = fu_set_pole(south_flag,  pole_lat,  pole_lon )
        IF (error)RETURN


        !Conversion from meters to degrees at equator
        fTmp=radians_to_degrees/earth_radius ! degrees per meter
        call grib_get(indGrib, 'DxInMetres', grid_dist_x_deg, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'DxInMetres',sub_name)
        grid_dist_x_deg = grid_dist_x_deg*fTmp

        call grib_get(indGrib,'DyInMetres',grid_dist_y_deg, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'DyInMetres ',sub_name)
        grid_dist_y_deg = grid_dist_y_deg*fTmp

        corner_in_geographical_latlon = .TRUE.

      !----------------------------------------
      !
      ! Unknown grid-type.
      !
      case default
        call set_error('Unknown grid type: "'//trim(strTmp)//'"',sub_name)
        return
    END SELECT


    grid = fu_set_lonlat_grid ('latlon_grid_from_grib', &
                   & corner_lon, corner_lat, &
                   & corner_in_geographical_latlon, &
                   & nx, ny, &
                   & pole, & 
                   & grid_dist_x_deg, grid_dist_y_deg)

  END FUNCTION fu_grib_grid


  ! ***************************************************************

  FUNCTION fu_grib_level(indGrib, source_of_data) RESULT(level)
    !
    ! Finds the level of grib-data.
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_level) :: level
    !
    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: indGrib
    type(meteo_data_source), intent(in) :: source_of_data

    ! Local variables
    integer :: iLevel, iLevelBottom, nCoef, io_status
    real :: a, b, a2, b2, fLevel, fLevelBottom
    real ::  p0sl, vcflat, sigma
    real, dimension(2*max_levels+2) :: fCoefs
    character(len=clen) :: chLevType, chTmp

    !
    ! All depends on the level type
    !
    call grib_get(indGrib,'typeOfLevel',chLevType, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfTypeOfLevel','fu_grib_level')

    call grib_get(indGrib,'level',fLevel, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'level','fu_grib_level')

    if(chLevType == 'surface')then
      level = ground_level

    elseif(chLevType == 'isobaricInhPa')then
      level = fu_set_pressure_level(fLevel*100.)

    elseif(chLevType == 'isobaricInPa')then
      level = fu_set_pressure_level(fLevel)

    elseif(chLevType == 'isobaricLayer')then
      call grib_get(indGrib,'bottomLevel',fLevelBottom, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'bottomLevel','fu_grib_level')
      level= fu_set_layer_between_two(layer_btw_2_pressure,fLevel,fLevelBottom)
      
    elseif(chLevType == 'meanSea')then
      level = mean_sea_level

    elseif(chLevType == 'heightAboveSea')then
      level = fu_set_constant_altitude_level(fLevel)

    elseif(chLevType == 'heightAboveSeaLayer')then
      call grib_get(indGrib,'bottomLevel',fLevelBottom, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'bottomLevel','fu_grib_level')
      level= fu_set_layer_between_two(layer_btw_2_altitude,fLevel*100,fLevelBottom*100)
      
    elseif(chLevType == 'heightAboveGround')then
      level = fu_set_constant_height_level(fLevel)

    elseif(chLevType == 'heightAboveGroundLayer')then
      call grib_get(indGrib,'bottomLevel',fLevelBottom, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'bottomLevel','fu_grib_level')
      level= fu_set_layer_between_two(layer_btw_2_height,fLevel*100,fLevelBottom*100)
      
    elseif(chLevType == 'sigma')then
      level = fu_set_sigma_level(fLevel / 10000.)

    elseif(chLevType == 'sigmaLayer')then
      call grib_get(indGrib,'bottomLevel',fLevelBottom, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'bottomLevel','fu_grib_level')
      level= fu_set_layer_between_two(layer_btw_2_sigma,fLevel*0.01,fLevelBottom*0.01)
      
    elseif(chLevType == 'hybrid')then
      !
      ! ATTENTION !!
      ! Either 2 parameters: a, b for the given level 
      ! or all a-s and b-s for HALF-LEVELS. So, half-sum is mandatory
      ! ATTENTION !!
      !
      call grib_get(indGrib,'level',chTMP, io_status)
      if (chTMP == 'MISSING') then  !! Asimoff message has hybrid missing level
        level = level_missing
        return
      endif


      call grib_get(indGrib,'level',iLevel, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'level','fu_grib_level')
      if (iLevel < 1) then
        call msg("iLevel from indGrib,", indGrib, iLevel)
        if(smpi_global_rank == 0) call grib_dump(indGrib, io_status)
        if(smpi_global_rank == 0) call backtrace_md()
        level = level_missing
        return
      endif

      call grib_get(indGrib,'numberOfVerticalCoordinateValues', nCoef, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,&
                                                           & 'numberOfVerticalCoordinateValues', &
                                                           & 'fu_grib_level')

      if(error)return
      call grib_get(indGrib,'pv',fCoefs, status=io_status) !, nCoef, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'pv', 'fu_grib_level')

      if(nCoef == 2)then
        a = fCoefs(1)
        b = fCoefs(2)

      elseif(fu_centre(source_of_data) == centre_moscow .and. fu_model(source_of_data) == model_cosmo) then
          !! Dirty hack for cosmo "hybrid-pressure" levels
          !! Unlike people are uesd to the surface pressure in cosmo is "standard" one,
          !! not the actual one, so levels stay fixed in space. 
          !! See "vertical coordinate parameters in the cosmo model grib"
          !! http://www.cosmo-model.org/content/model/modules/coding/grib/gribVerticalCoordinates.pdf
          !! 
          !! Below we pretend that it is ordinary hybrid sigma-pressure levels, though, coded
          !! in some wierd way...
          if (all(fCoefs(1) /= (/1,101/))) then
            call msg("pv values hyb level:",fCoefs(1:10))
            call set_error("Can't handle this verical of cosmo","fu_grib_level")
          endif
          p0sl=fCoefs(3)
          vcflat=fCoefs(6) 
          sigma=fCoefs(6+iLevel)
          if (sigma > vcflat) then
            a=vcflat*p0sl*(1.-sigma)/(1-vcflat)
            b=(sigma-vcflat)/(1-vcflat)
          else
            a = sigma*p0sl
            b = 0
          endif
      elseif(MOD(ncoef, 2) == 0)then ! Even number
        if (source_of_data == smhi_echam_src) then
          if(iLevel*2 > nCoef)then
            call set_error('Insufficient  123123 full-level coefs','fu_grib_level')
          else
            a = fCoefs(iLevel)
            b = fCoefs(iLevel+nCoef/2)
          endif
          
        elseif (ABS(fCoefs(nCoef-1)) >= 1.0e-9 .or. (source_of_data == centre_model_ecmwf_src))THEN
          if(iLevel*2 + 2 > nCoef)then
            call set_error('Insufficient  123123 half-level coefs','fu_grib_level')
          else
            a = 0.5 * (fCoefs(iLevel) + fCoefs(iLevel+1))
            b = 0.5 * (fCoefs(iLevel+nCoef/2) + fCoefs(iLevel+1+nCoef/2))
          endif
        else
          if(iLevel*2 + 2 > nCoef)then
            call set_error('Insufficient 112233 half-lev coefs','fu_grib_level')
          else
            a = 0.5 * (fCoefs(2*iLevel-1) + fCoefs(2*iLevel+1))
            b = 0.5 * (fCoefs(2*iLevel)   + fCoefs(2*iLevel+2))
          endif
        end if
        if (error) then !Report the trouble
          call msg_warning('Error set','fu_grib_level')
          call msg('Number of coefs: ', nCoef)
          call msg('Level number:', iLevel)
          call report(source_of_data)
          call msg("Coeffs array", fCoefs(1:nCoef))
          if(smpi_global_rank == 0) call grib_dump(indGrib, io_status)
          call set_error('The level number is bigger then the number of levels (DUP)','fu_grib_level')
          return

        endif
          
        if(b > 1)then
          call set_error('strange hybrid coefs','fu_grib_level')
          return
        endif
      else
        call set_error('Strange number of vertical parameters of hybrid level','fu_grib_level')
        call msg('Nbr of vertical parameters: ',nCoef)
        return
      end if
      level = fu_set_hybrid_level(iLevel, a, b) ! Level number, a, b

    elseif(chLevType == 'hybridLayer')then                 ! All coefs must be in psec

      call grib_get(indGrib, 'topLevel', iLevel, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status, 'topLevel','fu_grib_level')
      call grib_get(indGrib, 'bottomLevel', iLevelBottom, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'bottomLevel','fu_grib_level')
      call grib_get(indGrib,'numberOfVerticalCoordinateValues',nCoef, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status, & 
                                                           & 'numberOfVerticalCoordinateValues', &
                                                           & 'fu_grib_level')
      call grib_get(indGrib,'pv',fCoefs, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'pv','fu_grib_level')

      if(error)return

      if (source_of_data == hacked_cosmo_src) then 
        !Dirty hack for fake levels 
        if(iLevel*2 > nCoef)then
          call set_error('Insufficient 112233 Cosmo coefs','fu_grib_level')
        else
          a =  fCoefs(2*iLevel-1) 
          b =  fCoefs(2*iLevel)  
        endif
        level = fu_set_hybrid_level(iLevel, a, b) ! Level number, a, b
      else

      if(nCoef == 2)then
        call set_error('Layer between two hybrid requires all coefs','fu_grib_level')
        return
      elseif(nCoef == 4)then
        a = fCoefs(1)
        a2 = fCoefs(2)
        b = fCoefs(3)
        b2 = fCoefs(4)
      elseif(fu_centre(source_of_data) == centre_moscow .and. fu_model(source_of_data) == model_cosmo) then
          !! Dirty hack for cosmo "hybrid-pressure" levels see above
          if (all(fCoefs(1) /= (/1,101/))) then
            call msg("pv values hyb level:",fCoefs(1:10))
            call set_error("Can't handle this verical of cosmo","fu_grib_level")
          endif
          p0sl=fCoefs(3)
          vcflat=fCoefs(6) 
          sigma=0.5*(fCoefs(6+iLevel)+fCoefs(6+iLevel+1))
          if (sigma > vcflat) then
            a=vcflat*p0sl*(1.-sigma)/(1-vcflat)
            b=(sigma-vcflat)/(1-vcflat)
          else
            a = sigma*p0sl
            b = 0
          endif
          level = fu_set_hybrid_level(iLevel, a, b) ! Yet another dirty hack!
                ! (stupid) Silam not handle layers in meteo 
          return 
      elseif(MOD(ncoef, 2) == 0)then ! Even number
        if(iLevel*2 > nCoef)then
          if(smpi_global_rank == 0) call grib_dump(indGrib, io_status)
          call set_error('The level number is bigger then the number of levels2','fu_grib_level')
          return
        end if
        a = fCoefs(iLevel)
        a2 = fCoefs(iLevelBottom)
        b = fCoefs(iLevel+nCoef/2)
        b2 = fCoefs(iLevelBottom+nCoef/2)
      else
        call set_error('Strange number of vertical params for lyr btw 2 hybrids','fu_grib_level')
        call msg('Nbr of vertical parameters:',nCoef)
        return
      end if

      level = fu_set_layer_between_two(layer_btw_2_hybrid, &     ! Level type
                                     & iLevel,&        ! top
                                     & iLevelBottom, & ! bottom
                                     & a, b, a2, b2)
      endif !Hack source_of_data == hacked_cosmo_src

    elseif(chLevType == 'depthBelowLandLayer')then
      call grib_get(indGrib, 'topLevel', fLevel, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status, 'topLevel','fu_grib_level')
      call grib_get(indGrib, 'bottomLevel', fLevelBottom, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'bottomLevel','fu_grib_level')

      level = fu_set_layer_between_two(layer_btw_2_depth,fLevel*0.01,fLevelBottom*0.01)

!      case(191, 192, 193, 194, 195, 196, 197, 198) ! slope for radiation scheme: fake
!        level =  entire_atmosphere_mean_level 

    elseif(chLevType == 'entireAtmosphere')then
      level =  entire_atmosphere_mean_level 

    elseif(chLevType == 'nominalTop')then
      ! Pass quietly
      level = level_missing 

    elseif(chLevType == 'depthBelowSea')then
      ! Pass quietly
      level = level_missing 


    else
      if(chLevType == 'unknown')then !!last chance
        call grib_get(indGrib, 'levelType', chTmp, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'levelType','fu_grib_level')
        if (chTmp == 'sfc') then
          call msg("resetting sfc level to surface")
          level = ground_level
        else
          CALL msg('unknown level type found in grib:' + chLevType)
          level = level_missing
        endif
      else
          CALL msg('unknown level type found in grib:' + chLevType)
          level = level_missing
      endif
    endif  ! level type

  END FUNCTION fu_grib_level



!******************************************************************************
!
!  Writing of the GRIB file (including the private stuff)
!
!******************************************************************************


    !*************************************************************

    SUBROUTINE open_gribfile_o(dir, fname, grib_index, chTemplateStr, ifAppend)
    !
    ! Opens a new gribfile and returns it's unit (output !!). Should the directory
    ! dos not exist - tries to create it.
    !
    ! ATTENTION. inq_unit is NOT a FORTRAN unit, but rather C file handler
    ! needed for PBGRIB and GRIBWRITE
    ! grib_index is an index of the gfile structure containing all
    ! information about the GRIB file
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    CHARACTER (LEN=*), INTENT(in) :: dir, fname
    CHARACTER (LEN=*), INTENT(in), optional :: chTemplateStr
    logical, intent(in), optional :: ifAppend

    ! Imported parameters with intent OUT:
    INTEGER, INTENT(out) :: grib_index

    ! Local declarations:
    LOGICAL :: file_exists, file_open
    INTEGER :: io_status, inq_unit
    character(len=pbopen_string_len) :: pbopen_fname


    !----------------------------------------
    !
    ! 1. Check file. If opened - rewind, if not - open
    !    No check for existence !!

   
    if(len_trim(dir)+ len_trim(fname) == pbopen_string_len)then
      call set_error('File name should be <250 symbols','open_grib_file_o')
      return
    end if

    if(len_trim(dir) > 0)then  ! may be, the directory is empty
      pbopen_fname = trim(fu_connect_strings(dir,dir_slash,fname))
    else
      pbopen_fname = trim(fname)
    endif

    INQUIRE(file=pbopen_fname, exist=file_exists, opened=file_open, number=inq_unit)

    IF (file_open) THEN
      call msg('open_gribfile: Already open, rewinding...')
      close(inq_unit)
    END IF

    if(present(ifAppend))then
      if(ifAppend)then
        CALL grib_open_file(inq_unit, pbopen_fname, 'ab', io_status)
!        CALL grib_open_file_a_md(inq_unit, pbopen_fname)
!        call set_error('Cannot append yet','open_gribfile_o')
!        return
      else
        CALL grib_open_file(inq_unit, pbopen_fname, 'wb', io_status)
!        CALL grib_open_file_o_md(inq_unit, pbopen_fname)
      endif
    else
      CALL grib_open_file(inq_unit, pbopen_fname, 'wb', io_status)
!     CALL grib_open_file_o_md(inq_unit, pbopen_fname)
    endif

    IF (io_status /= 0) THEN
      call process_grib_status(io_status,'Error while opening gribfile:' + pbopen_fname,'open_gribfile_o')
      RETURN
    END IF

    IF (io_status /= 0) THEN
      !
      ! May be, there is no directory ? Or the whole sub-tree is missing ? Let's first
      ! try to create it.
      !
      call msg_warning('Can not open file. Checking the directory tree...','open_gribfile_o')
      call create_directory_tree(pbopen_fname(1:index(pbopen_fname,dir_slash,.true.)-1))
      if(error)then
        return
      else
        CALL grib_open_file(inq_unit, pbopen_fname, 'wb', io_status)
        IF (io_status /= 0) THEN
          call process_grib_status(io_status,'Error while opening gribfile:' + pbopen_fname,'open_gribfile_o')
          RETURN
        END IF
!        CALL grib_open_file_o_md(inq_unit,pbopen_fname)
      endif
    END IF

    if(present(chTemplateStr))then
      grib_index = init_ctl_for_grib(inq_unit, dir, fname, chTemplateStr)
    else
      grib_index = init_ctl_for_grib(inq_unit, dir, fname)
    endif

    END SUBROUTINE open_gribfile_o


  !*************************************************************

  SUBROUTINE switch_grib_binary_o(dir, fname, grib_index)
    !
    ! Closes existing GRIB binary file and opens a new one, keeping the ctl
    ! structure unchanged. Returns new file unit (output !!). Should directory
    ! dos not exist - tries to create it.
    !
    ! ATTENTION. inq_unit is NOT a FORTRAN unit, but rather C file handler
    ! needed for PBGRIB and GRIBWRITE
    ! grib_index is an index of the gfile structure containing all
    ! information about the GRIB file
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    CHARACTER (LEN=*), INTENT(in) :: dir, fname

    ! Imported parameters with intent INOUT:
    INTEGER, INTENT(inout) :: grib_index

    ! Local declarations:
    LOGICAL :: file_exists, file_open
    INTEGER :: io_status, inq_unit
    character(len=pbopen_string_len) :: pbopen_fname
    integer(kind=kindOfInt) :: ifile, status

    !
    ! Stupidity check first...
    !
    IF(io_status /= 0)THEN
      CALL msg_warning('Can not close GRIB file','switch_grib_binary_o')
    END IF
    if(len_trim(dir)+ len_trim(fname) == pbopen_string_len)then
      call set_error('File name should be <250 symbols','switch_grib_binary_o')
      return
    end if

    if(len_trim(dir) > 0)then  ! may be, the directory is empty
      pbopen_fname = trim(fu_connect_strings(dir,dir_slash,fname))
    else
      pbopen_fname = trim(fname)
    endif

    !
    ! First, closes the existing file
    !
    ifile = fu_unit_bin(grib_index)
    
    CALL grib_close_file(ifile,status)
    if(status == GRIB_SUCCESS)then
      io_status = 0
    else
      call process_grib_status(status,'call grib_close_file','switch_grib_binary_o')
      io_status = -1
    endif

    IF(io_status /= 0)THEN
      CALL msg_warning('Can not close GRIB file','switch_grib_binary_o')
    END IF
    if(error)return

    !
    ! Now, check the new file - if open, rewind
    !
    INQUIRE(file=pbopen_fname, exist=file_exists, opened=file_open, number=inq_unit)

    IF (file_open) THEN
      call msg('switch_grib_binary_o: Already open, rewinding...')
      close(inq_unit)
    END IF

    !
    ! Finally, OPEN the file
    !
    CALL grib_open_file(inq_unit, pbopen_fname, 'wb', io_status)
    IF (io_status /= 0) THEN
      call process_grib_status(io_status,'Error while opening gribfile:' + pbopen_fname,'switch_grib_binary_o')
      RETURN
    END IF
!    CALL grib_open_file_o_md(inq_unit,pbopen_fname)

    IF (io_status /= 0) THEN
      !
      ! May be, there is no directory ? Or the whole sub-tree is missing ? Let's first
      ! try to create it.
      !
      call msg_warning('Can not open file. Checking the directory tree...','switch_grib_binary_o')
      call create_directory_tree(pbopen_fname(1:index(pbopen_fname,dir_slash,.true.)-1))
      if(error)then
        return
      else
        CALL grib_open_file(inq_unit, pbopen_fname, 'wb', io_status)
        IF (io_status /= 0) THEN
          call process_grib_status(io_status,'Error while opening gribfile:' + pbopen_fname,'switch_grib_binary_o')
          RETURN
        END IF
!        CALL grib_open_file_o_md(inq_unit,pbopen_fname) !,'wb',io_status)
      endif
    END IF

    call set_grib_binary_unit(grib_index, inq_unit)

  END SUBROUTINE switch_grib_binary_o


  ! ***************************************************************

  SUBROUTINE write_next_field_to_gribfile (iGribType, igf_index, id, grid_data)

    ! Puts the field to an open grib-file, to which
    ! a unit-number grib_unit points. The field must have two parts:
    ! the identification (silja_field_id) and grid-data (real array).
    !
    ! Method:
    ! field_to_gridrecord to encode 
    ! grib_write to put it in
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    INTEGER, INTENT(in) :: iGribType, igf_index ! gfile where to write
    TYPE(silja_field_id), INTENT(in) :: id
    REAL, DIMENSION(:), INTENT(in) :: grid_data

    !Local variables
    INTEGER :: iRet, iTmp, indGribOut, iGribQuantity, iLevelType, iTimeRangeIndicator
    real(r8k), DIMENSION(worksize) :: gribrecord
    real :: fLevel
    
    if(igf_index == int_missing)then
      call set_error('Undefined GRIB index ', 'write_next_field_to_gribfile')
      return
    endif
    !
    ! Create the grib-headings from the field-identification.
    !
    CALL field_id_to_gribheadings (id, grid_data, &
                                    & iGribType, &
                                    & indGribOut, &
                                    & iGribQuantity, iLevelType, fLevel, iTimeRangeIndicator)
    IF(error)THEN
      call msg_warning('Failed GRIB coding but continue','write_next_field_to_gribfile')
      CALL unset_error('write_next_field_to_gribfile')
      RETURN
    endif

    !
    ! NUmber of cells must be given explicitly, otherwise the whole temporary array is packed. 
    ! Also, real*4 must not be sent - it will be copied to real(r8k) array
    !
    call grib_set(indGribOut, 'values', grid_data(1:fu_number_of_gridpoints(fu_grid(id))), iRet)
    if(iRet /= GRIB_SUCCESS)then
      call process_grib_status(iRet,'values','write_next_field_to_gribfile')
      return
    endif
    
    
!    call grib_get(indGribOut, 'getNumberOfValues', iTmp, iRet)
!    call msg('getNumberOfValues', iTmp)
    
    
    !
    !  Fill-in the ctl file structures with this field data
    !
    CALL fill_structures_for_grib(id, &
                                & iGribQuantity, iLevelType, fLevel, iTimeRangeIndicator, &
                                & igf_index)

    ! 
    ! Write to file and release the structure
    !
    CALL grib_write(indGribOut, &             ! message index in grib_api
                  & fu_unit_bin(igf_index), & ! file binary unit
                  & iRet)
    IF (iret /= GRIB_SUCCESS) THEN
      call set_error('Writing to GRIB failed','write_next_field_to_gribfile')
      CALL unset_error('write_next_field_to_gribfile')
      RETURN
    END IF

    call grib_release(indGribOut, iRet)
    IF (iRet /= GRIB_SUCCESS) THEN
      call set_error('Cannot release the grib field','write_next_field_to_gribfile')
      CALL unset_error('write_next_field_to_gribfile')
      RETURN
    END IF

  END SUBROUTINE write_next_field_to_gribfile


  ! ***************************************************************

  SUBROUTINE set_field_grid_to_grib(grid, indGrib) 

    ! Encodes the grid to grib-header.
    ! NOTE. ksec2(12) is already set during PDS section coding
    !  - see levels definition
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters:
    TYPE(silja_grid), INTENT(in) :: grid
    INTEGER, INTENT(in) :: indGrib

    ! Local declarations:
    TYPE(silam_pole) :: pole
    LOGICAL :: corner_in_geographical_latlon
    REAL :: grid_dist_x_deg, grid_dist_y_deg, &
          & southpole_lat_N, southpole_lon_E, &
          & corner_lat, corner_lon, dx_deg, dy_deg
    INTEGER :: nx, ny, io_status

    SELECT CASE (fu_gridtype(grid))
    
      CASE(lonlat) !-------------------- normal or rotated lon-lat grid

        CALL lonlat_grid_parameters(grid,&
                                  & corner_lon, corner_lat, corner_in_geographical_latlon, &
                                  & nx,  ny, &
                                  & southpole_lon_E, southpole_lat_N, & 
                                  & dx_deg, dy_deg)
        if(error)return
        !
        ! Set basic grid parameters
        !
        call grib_set(indGrib,'numberOfPointsAlongAParallel',nx,io_status)
        if(io_status/= GRIB_SUCCESS) call process_grib_status(io_status,'numberOfPointsAlongAParallel','set_field_grid_to_grib')

        call grib_set(indGrib,'numberOfPointsAlongAMeridian',ny, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'numberOfPointsAlongAMeridian','set_field_grid_to_grib')

        call grib_set(indGrib,'longitudeOfFirstGridPointInDegrees',corner_lon, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'longitudeOfFirstGridPointInDegrees','set_field_grid_to_grib')

        call grib_set(indGrib,'latitudeOfFirstGridPointInDegrees',corner_lat, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'latitudeOfFirstGridPointInDegrees','set_field_grid_to_grib')

        call grib_set(indGrib,'longitudeOfLastGridPointInDegrees',corner_lon+(nx-1)*dx_deg, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'longitudeOfLastGridPointInDegrees','set_field_grid_to_grib')

        call grib_set(indGrib,'latitudeOfLastGridPointInDegrees',corner_lat+(ny-1)*dy_deg, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'latitudeOfLastGridPointInDegrees','set_field_grid_to_grib')

        call grib_set(indGrib,'DiInDegrees',dx_deg, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'DiInDegrees','set_field_grid_to_grib')

        call grib_set(indGrib,'DjInDegrees',dy_deg, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'DjInDegrees','set_field_grid_to_grib')

        call grib_set(indGrib,'resolutionAndComponentFlags',128, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'resolutionAndComponentFlags','set_field_grid_to_grib')

        call grib_set(indGrib,'scanningMode',64, io_status)  ! normal FORTRAN mode: i++, j++, i-consequtive
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'scanningMode','set_field_grid_to_grib')

        !
        ! Parameters are set depending on the pole type
        !
        if(fu_pole(grid) == pole_geographical)then
          call grib_set(indGrib,'dataRepresentationType',0, io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'dataRepresentationType','set_field_grid_to_grib')
        else
          call grib_set(indGrib,'dataRepresentationType',10, io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'dataRepresentationType','set_field_grid_to_grib')
          call grib_set(indGrib,'longitudeOfSouthernPoleInDegrees',southpole_lon_E, io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'longitudeOfSouthernPoleInDegrees','set_field_grid_to_grib')
          call grib_set(indGrib,'latitudeOfSouthernPoleInDegrees',southpole_lat_N, io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'latitudeOfSouthernPoleInDegrees','set_field_grid_to_grib')
        endif

        call grib_set(indGrib,'earthIsOblate',0, io_status)  ! Earth is spherical, R=6367.47 km
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'earthIsOblate','set_field_grid_to_grib')

        call grib_set(indGrib,'uvRelativeToGrid',8, io_status)  ! u,v components are resolved to grid directions
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'uvRelativeToGrid','set_field_grid_to_grib')

      CASE(polarstereographic)

        call grib_set(indGrib,'dataRepresentationType',5, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'dataRepresentationType','set_field_grid_to_grib')

        call grib_set(indGrib,'scanningMode',64, io_status)  ! normal FORTRAN mode: i++, j++, i-consequtive
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'scanningMode','set_field_grid_to_grib')

        CALL set_error('Polar stereography is not working yet','set_field_grid_to_grib')

      CASE DEFAULT
        CALL set_error('Unknown grid type','set_field_grid_to_grib')
        RETURN
    END SELECT

  END SUBROUTINE set_field_grid_to_grib


  ! ***************************************************************

  SUBROUTINE field_id_to_gribheadings (id, grid_data, &
                                     & iGribType, indGribOut, &
                                     & iGribQuantity, iLevelType, fLevel, iTimeRangeIndicator)
    !
    ! Creates a complete GRIB-heading set from the field-identification.
    !
    IMPLICIT NONE

    ! Imported parameters:
    TYPE(silja_field_id), INTENT(in) :: id
    REAL, DIMENSION(:), INTENT(in) :: grid_data
    integer, intent(in) :: iGribType
    INTEGER, intent(out) :: indGribOut
    integer, intent(out) :: iGribQuantity, iLevelType, iTimeRangeIndicator
    REAL, intent(out) :: fLevel

    ! Local declarations:
    INTEGER :: io_status, iLevelType_required, iScale10, iNBits
    real :: fLevelValueRequired, fLevelScale

    !
    ! Create the message from a ready-made sample
    !
    call grib_new_from_samples(indGribOut,'GRIB1', io_status) !regular_ll_pl_grib1', io_status)
    if(io_status /= GRIB_SUCCESS)then
      call process_grib_status(io_status,'GRIB1','field_id_to_gribheadings')
      return
    endif
    !
    ! Product identification.
    !
    call grib_set(indGribOut,'editionNumber', iGribType, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'editionNumber','field_id_to_gribheadings')
    if(error)return
    !
    ! Set the basic parts of identification.
    ! 
!    call grib_set(indGribOut,'identificationOfOriginatingGeneratingCentre',centre_fmi, io_status)
!    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'identificationOfOriginatingGeneratingCentre','field_id_to_gribheadings')
    call grib_set(indGribOut,'centre',centre_fmi, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'centre','field_id_to_gribheadings');if(error)return

    call grib_set(indGribOut,'subCentre',int_missing, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'subCentre','make_grib_1_header');if(error)return

    call grib_set(indGribOut,'generatingProcessIdentifier',model_silam, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'generatingProcessIdentifier','field_id_to_gribheadings');if(error)return

    select case(iGribType)
      case (1)
        call make_grib_1_header(indGribOut, id, iGribQuantity)
      case (2)
        call make_grib_2_header(indGribOut, id, iGribQuantity)
      case default
        CALL msg('Unknown grib-edition:', iGribType)
        CALL set_error('Unknown grib-edition','field_id_to_gribheadings')
    end select
    if(error)return

!    call grib_set(indGribOut,'gridType',255, io_status)
!    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'gridType','field_id_to_gribheadings')  ! Unknown grid, see grid definition section
    
    SELECT CASE(fu_quantity(id))
      CASE(particle_counter_flag, &
         & areas_of_risk_flag, &
         & concentration_flag, &
         & mass_in_air_flag, &
         & drydep_flag, &
         & wetdep_flag, &
         & advection_moment_X_flag, advection_moment_Y_flag, advection_moment_Z_flag, &
         & emission_intensity_flag)
        call grib_set(indGribOut,'section1Flags',192, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'section1Flags','field_id_to_gribheadings')   !- both sections 2 and 3(bitmap) present
        if(error)return
        call grib_set(indGribOut,'bitmapPresent',0, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'bitmapPresent','field_id_to_gribheadings')   ! Bitmap follows
        if(error)return
        call grib_set(indGribOut,'missingValue', int_missing, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'missingValue','field_id_to_gribheadings')  ! Missing integer in the binary data array
        if(error)return
        call grib_set(indGribOut,'missingValue',fu_grib_missing_real(fu_quantity(id)), io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'missingValue','field_id_to_gribheadings') ! Missing real 
        if(error)return                                                                             ! in the binary data array.
      CASE DEFAULT
        call grib_set(indGribOut,'section1Flags',128, io_status)         !- only section 2 (GDS) present
    END SELECT

    !
    ! ATTENTION.
    ! Here there is an ambiguity connected with the fact that the GRIB -> SILAM coding
    ! involves level value in some cases. Here, however, we take SILAM quantity and level
    ! parameters from the field ID. Nothing assures that they are compatible with the GRIB
    ! assumptions. In principle, it should be checked. Otherwise the decoding algorithm
    ! will reject this very field due to incompatible level.
    !
!      & iLevelType_required, fLevelValueRequired, fLevelScale) ! univ, OUT
    !
    ! We are not going to rearrange the arrays, so have to overwrite the scanning mode
    !
    call grib_set(indGribOut,'iScansNegatively',1,io_status) ! 0 if rows enumerated DEscendingly
    if(io_status /= GRIB_SUCCESS) then
      call process_grib_status(io_status,'iScansNegatively','field_id_to_gribheadings')
!      return
    endif
    
    call grib_set(indGribOut,'jScansPositively',0,io_status)  ! 0 if columns enumerated Ascendingly
    if(io_status /= GRIB_SUCCESS) then
      call process_grib_status(io_status,'jScansPositively','field_id_to_gribheadings')
      return
    endif
    call grib_set(indGribOut,'jPointsAreConsecutive',0,io_status) ! 0 if consequtive values along row
    if(io_status /= GRIB_SUCCESS) then
      call process_grib_status(io_status,'jPointsAreConsecutive','field_id_to_gribheadings')
      return
    endif

    CALL set_field_level_to_grib(id, indGribOut, iLevelType, fLevel) ! Level indicators
    if(error)return

    CALL set_field_time_to_grib(id, indGribOut, iTimeRangeIndicator)   ! Time indicators
    if(error)return

    CALL set_field_grid_to_grib(fu_grid(id), indGribOut) 
    if(error)return

!    CALL fu_get_10_scale(id, grid_data, iScale10, iNBits) ! decimal scale factor + Nbr of bits per value
!    if(error)return
!    call grib_set(indGribOut,'decimalScaleFactor',iScale10,io_status)
!    if(io_status /= GRIB_SUCCESS) then
!      call process_grib_status(io_status,'decimalScaleFactor','field_id_to_gribheadings')
!      return
!    endif
!    call grib_set(indGribOut,'decimalPrecision',iNBits,io_status)
!    if(io_status /= GRIB_SUCCESS) then
!      call process_grib_status(io_status,'decimalPrecision','field_id_to_gribheadings')
!      return
!    endif

    CONTAINS


    !============================================================================

    subroutine make_grib_1_header(indGrib, id, iGribQuantity)
      !
      ! GRIB 1 processing using GRIB-API
      !
      implicit none

      integer, intent(in) :: indGrib
      type(silja_field_id), intent(in) :: id
      integer, intent(out) :: iGribQuantity

      ! Local declarations
      integer :: iDiscipline, iParamCategory, iParamNbr, iLevelType_required, iTableVersion
      real :: fLevelValueRequired
      !
      ! Set the quantity of grib. New function based on grib_code_table module
      !
      iTableVersion = accept_all
      call get_GRIB_codes(fu_quantity(id), fu_str(fu_species(id)), &  ! SILAM variable  IN
                        & 1, &                                     ! 1 or 2      IN
                        & centre_fmi, accept_all, model_silam, &   ! universal   IN
                        & iTableVersion, &                            ! for GRIB 1  INOUT
                              & iGribQuantity, &                             ! GRIB 1    OUT
                              & iDiscipline, iParamCategory, iParamNbr, &    ! GRIB 2    OUT
                              & iLevelType_required, fLevelValueRequired)    ! univ, OUT
      if(error)return
      !
      ! Find the quantity of grib. New function based on grib_code_table module
      ! It also returns the level if it is fully determined by the quantity
      ! Otherwise, level_missing returned
      !
      call grib_set(indGrib, 'table2Version', iTableVersion, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'table2Version','make_grib_1_header')
      if(error)return
      call grib_set(indGrib,'indicatorOfParameter', iGribQuantity, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfParameter','make_grib_1_header')
      if(error)return
      call grib_set(indGrib,'indicatorOfTypeOfLevel', iLevelType_required, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfTypeOfLevel','make_grib_1_header')
      if(error)return
      call grib_set(indGrib,'level', fLevelValueRequired, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'level','make_grib_1_header')
      if(error)return

    end subroutine make_grib_1_header


    !============================================================================

    subroutine make_grib_2_header(indGrib, id, iGribQuantity)
      !
      ! GRIB 1 processing using GRIB-API
      !
      implicit none

      integer, intent(in) :: indGrib
      type(silja_field_id), intent(in) :: id
      integer, intent(out) :: iGribQuantity

      ! Local declarations
      integer :: iDiscipline, iParamCategory, iParamNbr, iLevelType_required, iTableVersion
      real :: fLevelValueRequired

      !
      ! Set the quantity of grib. New function based on grib_code_table module
      !
      iTableVersion = accept_all
      call get_GRIB_codes(fu_quantity(id), fu_str(fu_species(id)), &  ! SILAM variable  IN
                        & 2, &                                     ! 1 or 2      IN
                        & centre_fmi, accept_all, model_silam, &   ! universal   IN
                        & iTableVersion, &                            ! for GRIB 1  INOUT
                              & iGribQuantity, &                             ! GRIB 1    OUT
                              & iDiscipline, iParamCategory, iParamNbr, &    ! GRIB 2    OUT
                              & iLevelType_required, fLevelValueRequired)    ! univ, OUT
      if(error)return
      !
      ! Find the quantity of grib. New function based on grib_code_table module
      ! It also returns the level if it is fully determined by the quantity
      ! Otherwise, level_missing returned
      !
      call grib_set(indGrib, 'tablesVersion', iTableVersion, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'tablesVersion','make_grib_2_header');if(error)return

!      call grib_get(indGrib, 'localTablesVersion', localTableVersion, io_status)
!      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'localTablesVersion','make_grib_2_header')
!      if(error)return

      call grib_set(indGrib,'discipline', iDiscipline, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfParameter','make_grib_2_header')
      if(error)return
      call grib_set(indGrib,'parameterCategory', iParamCategory, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfParameter','make_grib_2_header')
      if(error)return
      call grib_set(indGrib,'parameterNumber', iParamNbr, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfParameter','make_grib_2_header')
      if(error)return
      call grib_set(indGrib,'typeOfLevel', iLevelType_required, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'typeOfLevel','make_grib_2_header')
      if(error)return
      call grib_set(indGrib,'level', fLevelValueRequired, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'level','make_grib_2_header')
      if(error)return

    end subroutine make_grib_2_header


  END SUBROUTINE field_id_to_gribheadings


  ! ***************************************************************

  SUBROUTINE set_field_level_to_grib(id, indGrib, iLevelType, fLevel)

    ! Extracts the field level parameters and puts them into GRIB arrays
    !
    ! SILAM level flags are the same as WMO ones. So, they can be copied
    ! directly and then necessary numerical values defined (e.g. pressure
    ! for pressure-level-type, etc.)
    !
    IMPLICIT NONE

    ! Imported parameters
    TYPE(silja_field_id), INTENT(in) :: id
    integer, intent(in) :: indGrib
    integer, intent(out) :: iLevelType
    real, intent(out) :: fLevel
    
    ! Local variables
    real, dimension(2) :: arCoef
    integer :: io_status

    ! Checking
    !
    IF (.not.defined(fu_level(id)).or. fu_leveltype(fu_level(id)).eq.no_level) THEN
      CALL set_error('No-level flag detected in field id','set_field_level_to_grib')
      RETURN
    END IF

    ! Copy flag and, if needed, setup the numerical values for the level
    !
    iLevelType = fu_leveltype(fu_level(id))
    if(error)return
    call grib_set(indGrib,'indicatorOfTypeOfLevel',fu_leveltype(fu_level(id)), io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfTypeOfLevel','set_field_level_to_grib')

    
    SELECT CASE (fu_leveltype(fu_level(id)))
      CASE(constant_pressure)
        fLevel = fu_pr_level_pressure(fu_level(id))
        call grib_set(indGrib,'level', fLevel, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'level','set_field_level_to_grib')
        call grib_set(indGrib,'numberOfVerticalCoordinateValues',0, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'numberOfVerticalCoordinateValues','set_field_level_to_grib')

      CASE(constant_altitude)
        fLevel = fu_level_altitude(fu_level(id))
        call grib_set(indGrib,'level',fLevel , io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'level','set_field_level_to_grib')
        call grib_set(indGrib,'numberOfVerticalCoordinateValues',0, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'numberOfVerticalCoordinateValues','set_field_level_to_grib')
        
      CASE(constant_height)
        fLevel = fu_level_height(fu_level(id))
        call grib_set(indGrib,'level',fLevel , io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'level','set_field_level_to_grib')
        call grib_set(indGrib,'numberOfVerticalCoordinateValues',0, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'numberOfVerticalCoordinateValues','set_field_level_to_grib')

      CASE(sigma_level) ! For sigma level it is a sigma coefficient
        fLevel = fu_level_height(fu_level(id))*10000
        call grib_set(indGrib,'level',fLevel , io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'level','set_field_level_to_grib')
        call grib_set(indGrib,'numberOfVerticalCoordinateValues',0, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'numberOfVerticalCoordinateValues','set_field_level_to_grib')

      CASE(hybrid)
        fLevel = fu_hybrid_level_number(fu_level(id))
        call grib_set(indGrib,'level',fLevel , io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'level','set_field_level_to_grib')
        call grib_set(indGrib,'numberOfVerticalCoordinateValues',2, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'numberOfVerticalCoordinateValues','set_field_level_to_grib')
        if(error)return
        arCoef(1) = fu_hybrid_level_coeff_a(fu_level(id))
        arCoef(2) = fu_hybrid_level_coeff_b(fu_level(id))
        call grib_set(indGrib,'pv',arCoef, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'pv','set_field_level_to_grib')

      CASE(surface, top_of_the_atmosphere, entire_atmosphere_single_layer)

     CASE DEFAULT
        CALL set_error('Unknown field level','set_field_level_to_grib')
    END SELECT

  END SUBROUTINE set_field_level_to_grib


  ! ***************************************************************

  SUBROUTINE set_field_time_to_grib(id, indGrib, iTimeRangeIndicator)   ! Time indicators
    !
    ! Extracts time stuff from the field id and puts it into ksec1
    ! NOTE. Very tough problem with time indicators - field_id does not
    ! have enough flags to correctly interpret all possible options.
    ! The best-available guess is used in several cases (see comments below).
    ! 
    ! Code owner: Mikhail Sofiev, FMI

    IMPLICIT NONE

    TYPE(silja_field_id), INTENT(in) :: id
    INTEGER, INTENT(in) :: indGrib
    integer, intent(out) :: iTimeRangeIndicator

    ! Local declarations
    TYPE(silja_interval) UTC_diff
    integer :: yr, mon, day, hr, min, io_status, iDateTime
    REAL :: fTmp
!    type(silam_sp) :: sp
    character(len=20) :: chTmp
    
!    sp%sp => fu_work_string()
!    if(error)return
    !
    !  Setup the reference time - exact correspondence
    !
    chTmp = fu_time_to_date_string(fu_analysis_time(id))
    read(unit=chTmp,fmt='(i8)')iDateTime                    ! yyyymmdd
    call grib_set(indGrib,'dataDate',iDateTime, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'dataDate','set_field_time_to_grib')

    chTmp = fu_time_to_daily_time_string(fu_analysis_time(id))
    read(chTmp, fmt='(i6)')iDateTime                        ! hhmmss
    call grib_set(indGrib,'dataTime',iDateTime, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'dataTime','set_field_time_to_grib')
    IF (error) RETURN

    call grib_set(indGrib,'indicatorOfUnitOfTimeRange','m', io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfUnitOfTimeRange','set_field_time_to_grib')

    !
    ! Setup the P1, P2 and time range indicator
    ! Also a number of products included in the average is set to 1 if averaged_flag
    !
    iTimeRangeIndicator = fu_field_kind(id)
    SELECT CASE(fu_field_kind(id))

    CASE(accumulated_flag, averaged_flag)

      if(fu_days(fu_forecast_length(id)-fu_accumulation_length(id)) > 255 .or. &
       & fu_days(fu_forecast_length(id)) > 255)then
        call grib_set(indGrib,'indicatorOfUnitOfTimeRange','M', io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfUnitOfTimeRange','set_field_time_to_grib')
        call grib_set(indGrib,'periodOfTime', nint(fu_days(fu_forecast_length(id)-fu_accumulation_length(id))/30.5), io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'periodOfTime','set_field_time_to_grib')
        call grib_set(indGrib,'periodOfTimeIntervals', nint(fu_days(fu_forecast_length(id))/30.5), io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'periodOfTimeIntervals','set_field_time_to_grib')
        
      elseif(fu_hours(fu_forecast_length(id)-fu_accumulation_length(id)) > 255 .or. &
           & fu_hours(fu_forecast_length(id)) > 255)then
        call grib_set(indGrib,'indicatorOfUnitOfTimeRange','D', io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfUnitOfTimeRange','set_field_time_to_grib')
        call grib_set(indGrib,'periodOfTime', nint(fu_days(fu_forecast_length(id)-fu_accumulation_length(id))), io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'periodOfTime','set_field_time_to_grib')
        call grib_set(indGrib,'periodOfTimeIntervals', nint(fu_days(fu_forecast_length(id))), io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'periodOfTimeIntervals','set_field_time_to_grib')
        
      elseif(fu_mins(fu_forecast_length(id)-fu_accumulation_length(id)) > 255 .or. &
       & fu_mins(fu_forecast_length(id)) > 255)then
        call grib_set(indGrib,'indicatorOfUnitOfTimeRange','h', io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfUnitOfTimeRange','set_field_time_to_grib')
        call grib_set(indGrib,'periodOfTime', nint(fu_hours(fu_forecast_length(id)-fu_accumulation_length(id))), io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'periodOfTime','set_field_time_to_grib')
        call grib_set(indGrib,'periodOfTimeIntervals', nint(fu_hours(fu_forecast_length(id))), io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'periodOfTimeIntervals','set_field_time_to_grib')
        
      else
        call grib_set(indGrib,'indicatorOfUnitOfTimeRange','m', io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfUnitOfTimeRange','set_field_time_to_grib')
        call grib_set(indGrib,'periodOfTime', nint(fu_mins(fu_forecast_length(id)-fu_accumulation_length(id))), io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'periodOfTime','set_field_time_to_grib')
        call grib_set(indGrib,'periodOfTimeIntervals', nint(fu_mins(fu_forecast_length(id))), io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'periodOfTimeIntervals','set_field_time_to_grib')
      endif
      
      if(fu_field_kind(id) == accumulated_flag)then
        call grib_set(indGrib,'stepType','accum', io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'stepType','set_field_time_to_grib')
        call grib_set(indGrib,'numberIncludedInAverage',0, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'numberIncludedInAverage','set_field_time_to_grib')
        
      elseif(fu_field_kind(id) == averaged_flag)then
        call grib_set(indGrib,'stepType','avg', io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'stepType','set_field_time_to_grib')
        call grib_set(indGrib,'numberIncludedInAverage',1, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'numberIncludedInAverage','set_field_time_to_grib')
        
      else
        call set_error('Unknown type fo accumulated field','set_field_time_to_grib')
        return
      endif
      
      if(fu_validity_length(id) > zero_interval) &
        & call msg_warning('Period-valid field cannot be accumulated in GRIB', 'set_field_time_to_grib')

    CASE(forecast_flag)
      !
      ! Can be a forecast field or period_valid field: GRIB has one more parameter
      !
      if(fu_validity_length(id) == zero_interval)then
        !
        ! Normal forecast field
        !
        if(fu_days(fu_forecast_length(id)) > 255)then
          call grib_set(indGrib,'indicatorOfUnitOfTimeRange','M', io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfUnitOfTimeRange','set_field_time_to_grib')
          call grib_set(indGrib,'periodOfTime', nint(fu_days(fu_forecast_length(id))/30.5), io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'periodOfTime','set_field_time_to_grib')

        elseif(fu_hours(fu_forecast_length(id)) > 255)then
          call grib_set(indGrib,'indicatorOfUnitOfTimeRange','D', io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfUnitOfTimeRange','set_field_time_to_grib')
          call grib_set(indGrib,'periodOfTime', nint(fu_days(fu_forecast_length(id))), io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'periodOfTime','set_field_time_to_grib')

        elseif(fu_mins(fu_forecast_length(id)) > 255)then
          call grib_set(indGrib,'indicatorOfUnitOfTimeRange','h', io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfUnitOfTimeRange','set_field_time_to_grib')
          call grib_set(indGrib,'periodOfTime', nint(fu_hours(fu_forecast_length(id))), io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'periodOfTime','set_field_time_to_grib')

        else
          call grib_set(indGrib,'indicatorOfUnitOfTimeRange','m', io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfUnitOfTimeRange','set_field_time_to_grib')
          call grib_set(indGrib,'periodOfTime', nint(fu_mins(fu_forecast_length(id))), io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'periodOfTime','set_field_time_to_grib')
        endif

        call grib_set(indGrib,'stepType','instant', io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'stepType','set_field_time_to_grib')
        call grib_set(indGrib,'numberIncludedInAverage',0, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'numberIncludedInAverage','set_field_time_to_grib')
        
      else
        !
        ! Period valid field
        !
        if(fu_days(fu_forecast_length(id)) > 255 .or. &
         & fu_days(fu_forecast_length(id) + fu_validity_length(id)) > 255)then
          call grib_set(indGrib,'indicatorOfUnitOfTimeRange','M', io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfUnitOfTimeRange','set_field_time_to_grib')
          call grib_set(indGrib,'periodOfTime', nint(fu_days(fu_forecast_length(id))/30.5), io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'periodOfTime','set_field_time_to_grib')
          call grib_set(indGrib,'periodOfTimeIntervals', nint(fu_days(fu_forecast_length(id) + fu_validity_length(id))/30.5), io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'periodOfTimeIntervals','set_field_time_to_grib')

        elseif(fu_hours(fu_forecast_length(id)) > 255 .or. &
             & fu_hours(fu_forecast_length(id) + fu_validity_length(id)) > 255)then
          call grib_set(indGrib,'indicatorOfUnitOfTimeRange','D', io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfUnitOfTimeRange','set_field_time_to_grib')
          call grib_set(indGrib,'periodOfTime', nint(fu_days(fu_forecast_length(id))), io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'periodOfTime','set_field_time_to_grib')
          call grib_set(indGrib,'periodOfTimeIntervals', nint(fu_days(fu_forecast_length(id) + fu_validity_length(id))), io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'periodOfTimeIntervals','set_field_time_to_grib')

        elseif(fu_mins(fu_forecast_length(id)) > 255 .or. &
             & fu_mins(fu_forecast_length(id) + fu_validity_length(id)) > 255)then
          call grib_set(indGrib,'indicatorOfUnitOfTimeRange','h', io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfUnitOfTimeRange','set_field_time_to_grib')
          call grib_set(indGrib,'periodOfTime', nint(fu_hours(fu_forecast_length(id))), io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'periodOfTime','set_field_time_to_grib')
          call grib_set(indGrib,'periodOfTimeIntervals', nint(fu_hours(fu_forecast_length(id) + fu_validity_length(id))), io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'periodOfTimeIntervals','set_field_time_to_grib')

        else
          call grib_set(indGrib,'indicatorOfUnitOfTimeRange','m', io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfUnitOfTimeRange','set_field_time_to_grib')
          call grib_set(indGrib,'periodOfTime', nint(fu_mins(fu_forecast_length(id))), io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'periodOfTime','set_field_time_to_grib')
          call grib_set(indGrib,'periodOfTimeIntervals', nint(fu_mins(fu_forecast_length(id) + fu_validity_length(id))), io_status)
          if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'periodOfTimeIntervals','set_field_time_to_grib')
        endif
        call grib_set(indGrib,'timeRangeIndicator',2, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'timeRangeIndicator','set_field_time_to_grib')
        call grib_set(indGrib,'numberIncludedInAverage',0, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'numberIncludedInAverage','set_field_time_to_grib')
      endif

    CASE DEFAULT
      CALL set_error('Unknown field kind','set_field_time_to_grib')
    END SELECT

  END SUBROUTINE set_field_time_to_grib


!  ! ***************************************************************
!
!  SUBROUTINE fu_get_10_scale (id, grid_data, scale_10, NbrOfBits)
!    !
!    ! Detects the nesessary conversion 10-factor for GRIB code procedure
!    ! Method:
!    ! All what survives after NINT operation will be coded, so 10-factor 
!    ! is selected so that the whole range can be
!    ! covered with 4-digits precision (up to ~16000, to be exact). 
!    ! And then the number of bits for the maximum value is set 16000=>14
!    !
!    ! 
!    IMPLICIT NONE
!
!    REAL, PARAMETER :: log_10_16000 = 4.2!
!
!   TYPE(silja_field_id), INTENT(in) :: id
!    REAL, DIMENSION (:), INTENT(in) :: grid_data
!
!    INTEGER, INTENT(out) :: scale_10, NbrOfBits
!
!    ! Local variables:
!    REAL :: fMax = -int_missing
!    REAL :: fMin = int_missing
!    INTEGER :: iNbrOfPoints, iTmp
!
!    IF(.not.defined(id).or..not.defined(fu_grid(id)))THEN
!      CALL set_error('Undefined field or grid','fu_get_ao_scale')
!    ENDIF
!
!    iNbrOfPoints = fu_number_of_gridpoints(fu_grid(id))
!
!    DO iTmp = 1, iNbrOfPoints
!      IF(.not. (grid_data(iTmp) .eps. real_missing))THEN
!        fMax = MAX(grid_data(iTmp),fMax)
!        fMin = MIN(grid_data(iTmp),fMin)
!      END IF
!    END DO
!    !
!    ! Check - may be field is completely undefined
!    !
!    IF((fMin .eps. real_missing) .or. (fMax .eps. real_missing) .or. (fMax .eps. fMin))THEN
!        CALL msg_warning('completely undefined field','fu_get_10_scale')
!        scale_10 = 0
!        RETURN
!    ENDIF
!
!    scale_10 = - INT(log_10_16000 - LOG10(fMax-fMin))
!    NbrOfBits = 14
!
!  END SUBROUTINE fu_get_10_scale


    !*************************************************************

  SUBROUTINE close_gribfile_o(grib_index)

    ! Closes a gribfile. grib_index is an index of the gfile structure
    ! with ctl information for this file
    !
    ! Author: Mikhail Sofiev, FMI
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    INTEGER, INTENT(in) :: grib_index
    ! Local declarations:
    INTEGER io_status

    CALL grib_close_file(fu_unit_bin(grib_index), io_status)
    IF(io_status /= 0)THEN
      CALL msg_warning('Can not close GRIB file','close_gribfile_o')
    END IF

    CALL write_ctl_for_grib(grib_index)
    if(error)call unset_error('close_gribfile_o')

    CALL release_index(grib_index)
    
  END SUBROUTINE close_gribfile_o


  !***********************************************************************

  subroutine process_grib_status(io_status, chRequest, chPlace)
    !
    ! Prints the informative messages about the GRIB data. Mainly Used for 
    ! non-error informative mesages
    !
    implicit none

    ! Imported parameters with intent IN
    integer, intent(in) :: io_status
    character(len=*), intent(in) :: chRequest, chPlace

    ! local variables
    integer :: status_local
    character(len=2000) :: error_message

    select case(io_status)

      case (GRIB_SUCCESS)
        call msg('GRIB-API success')
        
      case default
        CALL msg_warning('Non-zero GRIB-API status at:' + chPlace + ',while requesting:' + chRequest)
        if(error)return
        error_message = ""
        call grib_get_error_string (io_status, error_message, status_local)
        if(status_local == GRIB_SUCCESS)then
          call set_error(error_message,'process_grib_status')
        else
          call set_error('Cannot retrieve the error text','process_grib_status')
        endif

    end select

  end subroutine process_grib_status


  ! ***************************************************************
  ! 
  !
  !                 TESTS
  !
  !
  ! ***************************************************************

  SUBROUTINE print_gribheadings(indGrib)    ! data themselves
    !
    ! Print grib-heading to a tmp-file using GRIBEX library routines.
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: indGrib
    
    ! Local variables
    integer :: io_status

    call msg('Writing gribheadings')

    if(smpi_global_rank == 0) then 
         call grib_dump(indGrib, io_status)
         if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'call grib_dump','print_gribheadings')
    endif

  END SUBROUTINE print_gribheadings


  ! ***************************************************************

  SUBROUTINE grib_tests(io_switch)

    ! Language: ANSI Fortran 90
    !
    IMPLICIT NONE
    
    ! Imported parameter
    integer, intent(in) :: io_switch
    
    !
    ! Local declarations:
    INTEGER :: unit
    CHARACTER (LEN=clen) :: fname1, fname2
    TYPE(silja_time) :: time
!    TYPE(silja_stack), TARGET :: stack
!    TYPE(silja_shopping_list) :: list
    TYPE(silja_grid) :: grid_in
    TYPE(silja_field), TARGET :: field
    TYPE(silja_field), POINTER :: fieldpointer
    TYPE(silja_field_id) :: id
    INTEGER :: grib_unit, iTmp, io_status
    LOGICAL :: eof
    CHARACTER (LEN=fnlen) :: fname
    INTEGER :: i, j, indGribOut, iTimeRangeIndicator, nMsg
    REAL, DIMENSION(:), pointer :: grid_data
    !    INTEGER :: GetGribMsgByFn
    CHARACTER (LEN=worksize) :: gribstring
    real :: scale_factor
    logical :: ifZeroNegatives
    type(silja_field_id),dimension(:), pointer :: idList  !Returns a pointer
    

    if(io_switch == 1)then


    fname = '/home/met/silja/program_root/hirlam_test_data/pp/GRIB_86_1_11_100_700_0_9804290000_07.170'

    CALL open_grib_file_i(fname,grib_unit, interval_missing)
    IF (error) RETURN
    call id_list_from_grib_file(grib_unit, idList, nMsg)

    DO i = 1,nMsg

      call report(idList(i))


    END DO

    CALL close_grib_file_i(grib_unit) 

    else
      !
      ! Test for writing
      !

id = fu_set_field_id(met_src_missing,&
                                     & temperature_flag, &
                                     & fu_wallclock(),&
                                     & one_hour, &
                                     & geo_global_grid,&
                                     & surface_level)

grid_data => fu_work_array(fu_number_of_gridpoints(geo_global_grid))
iTmp = 1
do io_status = 1, fu_number_of_gridpoints(geo_global_grid)
  grid_data(io_status) = iTMp
  iTmp = iTmp + 1
end do

    CALL grib_open_file(grib_unit, 'try_grib', 'wb', io_status)
    if(io_status /= GRIB_SUCCESS)then
      call process_grib_status(io_status,'Open_file','field_id_to_gribheadings')
    endif
    !
    ! Create the message from a ready-made sample
    !
    call grib_new_from_samples(indGribOut,'GRIB1', io_status)
    if(io_status /= GRIB_SUCCESS)then
      call process_grib_status(io_status,'GRIB1','field_id_to_gribheadings')
    endif
    !
    ! Product identification.
    !
    call grib_set(indGribOut,'editionNumber', 1, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'editionNumber','field_id_to_gribheadings')
    !
    ! Set the basic parts of identification.
    ! 
!    call grib_set(indGribOut,'identificationOfOriginatingGeneratingCentre',centre_fmi, io_status)
!    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'identificationOfOriginatingGeneratingCentre','field_id_to_gribheadings')
    call grib_set(indGribOut,'centre',centre_fmi, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'centre','field_id_to_gribheadings')

    call grib_set(indGribOut,'subCentre',int_missing, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'subCentre','make_grib_1_header')

    call grib_set(indGribOut,'generatingProcessIdentifier',model_silam, io_status)
    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'generatingProcessIdentifier','field_id_to_gribheadings')

!    call grib_set(indGribOut,'gridType',255, io_status)
!    if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'gridType','field_id_to_gribheadings')  ! Unknown grid, see grid definition section
    
        call grib_set(indGribOut,'section1Flags',192, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'section1Flags','field_id_to_gribheadings')   !- both sections 2 and 3(bitmap) present
        call grib_set(indGribOut,'bitmapPresent',0, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'bitmapPresent','field_id_to_gribheadings')   ! Bitmap follows
        call grib_set(indGribOut,'missingValue', int_missing, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'missingValue','field_id_to_gribheadings')  ! Missing integer in the binary data array
        call grib_set(indGribOut,'missingValue',real_missing, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'missingValue','field_id_to_gribheadings') ! Missing real 

    !
    ! We are not going to rearrange the arrays, so have to overwrite the scanning mode
    !
    call grib_set(indGribOut,'iScansNegatively',1,io_status) ! 0 if rows enumerated DEscendingly
    if(io_status /= GRIB_SUCCESS) then
      call process_grib_status(io_status,'iScansNegatively','field_id_to_gribheadings')
    endif
    
    call grib_set(indGribOut,'jScansPositively',0,io_status)  ! 0 if columns enumerated Ascendingly
    if(io_status /= GRIB_SUCCESS) then
      call process_grib_status(io_status,'jScansPositively','field_id_to_gribheadings')
    endif
    call grib_set(indGribOut,'jPointsAreConsecutive',0,io_status) ! 0 if consequtive values along row
    if(io_status /= GRIB_SUCCESS) then
      call process_grib_status(io_status,'jPointsAreConsecutive','field_id_to_gribheadings')
    endif

        call grib_set(indGribOut,'level', 10.0, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'level','set_field_level_to_grib')

        call grib_set(indGribOut,'numberOfVerticalCoordinateValues',0, io_status)
        if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'numberOfVerticalCoordinateValues','set_field_level_to_grib')

    CALL set_field_time_to_grib(id, indGribOut, iTimeRangeIndicator)   ! Time indicators

!    CALL set_field_grid_to_grib(geo_global_grid, indGribOut) 
!    call grib_set(indGribOut,'decimalScaleFactor',1,io_status)
!    if(io_status /= GRIB_SUCCESS) then
!      call process_grib_status(io_status,'decimalScaleFactor','field_id_to_gribheadings')
!    endif
!    call grib_set(indGribOut,'decimalPrecision',16,io_status)
!    if(io_status /= GRIB_SUCCESS) then
!      call process_grib_status(io_status,'decimalPrecision','field_id_to_gribheadings')
!    endif

    call grib_set(indGribOut,'complexPacking',0,io_status)
    if(io_status /= GRIB_SUCCESS) then
      call process_grib_status(io_status,'complexPacking','field_id_to_gribheadings')
    endif

      call grib_set(indGribOut, 'table2Version', 130, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'table2Version','make_grib_1_header')
      call grib_set(indGribOut,'indicatorOfParameter', 123, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfParameter','make_grib_1_header')
      call grib_set(indGribOut,'indicatorOfTypeOfLevel', 1, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'indicatorOfTypeOfLevel','make_grib_1_header')
      call grib_set(indGribOut,'level', 10.0, io_status)
      if(io_status /= GRIB_SUCCESS) call process_grib_status(io_status,'level','make_grib_1_header')

!    io_status = grib_f_set_real4_array(indGribOut, 'values', grid_data, fu_number_of_gridpoints(fu_grid(id)))
    call grib_set(indGribOut, 'values', grid_data(1:fu_number_of_gridpoints(fu_grid(id))), io_status)
    if(io_status /= GRIB_SUCCESS)then
      call process_grib_status(io_status,'values','write_next_field_to_gribfile')
    endif
    
    
    call grib_get(indGribOut, 'getNumberOfValues', iTmp, io_status)
    call msg('getNumberOfValues', iTmp)
    
    ! 
    ! Write to file and release the structure
    !
    CALL grib_write(indGribOut, &             ! message index in grib_api
                  & grib_unit, & ! file binary unit
                  & io_status)
    IF (io_status /= GRIB_SUCCESS) THEN
      call set_error('Writing to GRIB failed','write_next_field_to_gribfile')
      CALL unset_error('write_next_field_to_gribfile')
    END IF

    call grib_release(indGribOut, io_status)
    IF (io_status /= GRIB_SUCCESS) THEN
      call set_error('Cannot release the grib field','write_next_field_to_gribfile')
      CALL unset_error('write_next_field_to_gribfile')
    END IF

    CALL grib_close_file(grib_unit,io_status)



    endif


  END SUBROUTINE grib_tests



END MODULE grib_api_io
