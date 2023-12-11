module optimisation

!===========================================================================
!
!  Current module contains a simple set for time tracking inside the 
!  SILAM model. 
!  The second part is the process-gathering and handling routines, together
!  with progress meter rotines
!
!  All units: SI
!
!  Language: ANSI FORTRAN-90
!
!  Code owner: Mikhail Sofiev, FMI  mikhail.sofiev@fmi.fi
!
!============================================================================

use silam_times

implicit none

private

public start_count  ! Counters
public stop_count
public report_time


public init_progress_file    ! progress file
public write_progress_file
public fu_compute_progress_time

private report_time_counter
private report_all_times
private stop_count_full
private stop_count_short
private start_count_full
private start_count_short


interface report_time
  module procedure report_time_counter
  module procedure report_all_times
end interface

interface start_count
   module procedure start_count_short
   module procedure start_count_full
end interface

interface stop_count
   module procedure stop_count_short
   module procedure stop_count_full
end interface


!================================================================
!
! The basic type for the counting of the consumed time
!
type Ttime_counter
  private
           ! silja_true => counter is counting (time_start is reasonable)
           ! silja_false => occupied (time_passed is reasonable)
           ! silja_undefined => free unit
  type(silja_logical) :: busy=silja_undefined 

  character(len=clen) :: chNm=''  ! Name of the counter
  type(silja_time) :: time_start=time_missing       ! When the counter is turned on
  type(silja_interval) :: time_passed=zero_interval ! How much time accumulated
end type Ttime_counter

private Ttime_counter

integer, private, parameter :: n_counters = 100

type(Ttime_counter), dimension(n_counters), private, save :: TimeCounters


CONTAINS


!**********************************************************************
  
subroutine start_count_short(counter_name)
  ! A shorthand for calling stop_count(chCounterNm=counter_name),
  ! counter_name can now be a constant.
  implicit none
  character(len=*), intent(in) :: counter_name
  
  character(len=len(counter_name)) :: counter_name_tmp
  counter_name_tmp = counter_name
  call start_count_full(chCounterNm=counter_name_tmp)
end subroutine start_count_short

!************************************************************************************

subroutine start_count_full(indCounter, chCounterNm)
  !
  ! Starts or continues counter pointed by index or by name
  ! Should index meaningfull - returns the name, otherwise
  ! chacks the name and, if meaningfull - returns the index
  ! If both meaningless - creates new counter and returns both
  !
  implicit none

  ! Imported parameters
  integer, intent(inout), optional :: indCounter
  character(len=*), intent(inout), optional :: chCounterNm

  ! Local variables
  integer :: i

  ! Stupidity check - at least one of parameters must be present
  !
  if(.not.(present(indCounter) .or. present(chCounterNm)))then
    call set_error('At least one of index or name must be present','start_count')
    return
  endif

  ! If index is reasonable
  if(present(indCounter))then
    if(indCounter > 0 .and. indCounter < n_counters)then
      if(present(chCounterNm))chCounterNm = TimeCounters(indCounter)%chNm
      TimeCounters(indCounter)%busy = silja_true
      TimeCounters(indCounter)%time_start = fu_wallclock()
      return
    endif
  endif
  !
  ! Index is irrelevant - try name
  !
  if(present(chCounterNm))then
    if(chCounterNm /= '')then  
      !
      ! name is non-empty, try to find
      !
      do i=1,n_counters
        if(TimeCounters(i)%chNm == chCounterNM)then  ! counter found
          if(present(indCounter)) indCounter = i
          TimeCounters(i)%busy = silja_true
          TimeCounters(i)%time_start = fu_wallclock()
          return
        endif
      end do
      !
      ! Counter is not found - create new one with this given name
      !
      do i=1, n_counters
        if(TimeCounters(i)%busy == silja_undefined)then
          if(present(indCounter)) indCounter = i
          TimeCounters(i)%chNm = chCounterNm
          TimeCounters(i)%busy = silja_true
          TimeCounters(i)%time_passed = zero_interval
          TimeCounters(i)%time_start = fu_wallclock()
          return
        endif
      end do
      call msg_warning('Sorry, no free counters','start_count')
      return
    else
      !
      ! Name and index are empty - create the new counter with artificial name
      !
      do i=1, n_counters
        if(TimeCounters(i)%busy == silja_undefined)then
          if(present(indCounter)) indCounter = i
          write(unit=TimeCounters(i)%chNm, fmt='(A,I12)')'Time_Counter_',i
          if(present(chCounterNm)) chCounterNm = trim(TimeCounters(i)%chNm)
          TimeCounters(i)%busy = silja_true
          TimeCounters(i)%time_passed = zero_interval
          TimeCounters(i)%time_start = fu_wallclock()
          return
        endif
      end do
      call msg_warning('Sorry, no free counters','start_count')
      return
    endif  ! If name is empty
  endif  ! if name is present

end subroutine start_count_full

!************************************************************************************

subroutine stop_count_short(counter_name)
  ! A shorthand for calling stop_count(chCounterNm=counter_name),
  ! counter_name can now be a constant.
  implicit none
  character(len=*), intent(in) :: counter_name
  
  character(len=len(counter_name)) :: counter_name_tmp
  
  counter_name_tmp = counter_name
  call stop_count_full(chCounterNm=counter_name_tmp)
end subroutine stop_count_short


!*********************************************************************

subroutine stop_count_full(indCounter, chCounterNm)
  !
  ! Stops the time counting and adds passed time to the interval.
  ! First index is checked, then name. At least one of them should point 
  ! to the right counter, otherwise nothing happens. Correct values
  ! for the counter are always returned if found
  !
  implicit none

  ! Imported parameters
  integer, intent(inout), optional :: indCounter
  character(len=*), intent(inout), optional :: chCounterNm

  ! Local variables
  integer :: i

  ! Stupidity check - at least one of parameters must be present
  !
  if(.not.(present(indCounter) .or. present(chCounterNm)))then
    call set_error('At least one of index or name must be present','start_count')
    return
  endif

  ! If index is reasonable
  !
  if(present(indCounter))then
    if(indCounter > 0 .and. indCounter < n_counters)then
      if(TimeCounters(indCounter)%busy == silja_true)then
        if(present(chCounterNm)) chCounterNm = TimeCounters(indCounter)%chNm
        TimeCounters(indCounter)%busy = silja_false
        TimeCounters(indCounter)%time_passed = TimeCounters(indCounter)%time_passed + &
               & (fu_wallclock() - TimeCounters(indCounter)%time_start)
      else
        call msg_warning(fu_connect_strings('Counter is not counting:',chCounterNm),'stop_count')
        return
      endif
      return
    endif
  endif

  ! Index is irrelevant - search for the name
  !
  if(present(chCounterNm))then
    do i=1, n_counters
      if(TimeCounters(i)%chNm == chCounterNm)then
        if(present(indCounter))indCounter = i
        if(TimeCounters(i)%busy == silja_true)then
          TimeCounters(i)%busy = silja_false
          TimeCounters(i)%time_passed = TimeCounters(i)%time_passed + &
                                      & (fu_wallclock() - TimeCounters(i)%time_start)
          return
        else
          call msg_warning(fu_connect_strings('Counter is not counting:',chCounterNm),'stop_count')
          return
        endif
      endif
    end do  ! Search for the name
    call msg_warning('Sorry, counter is not found','stop_count')
    return
  endif  ! If name is present

end subroutine stop_count_full


!*******************************************************************

subroutine report_all_times(time_passed)
  !
  ! Reports time passed for all active counters. Should the time_passed present
  ! the sum of the time counted will be returned in it
  ! 
  implicit none

  ! Imported parameters
  type(silja_interval), intent(out), optional :: time_passed

  ! Local variables
  integer :: i, iTmp
  type(silja_interval) :: intervTmp
  character(len=clen) :: strTmp=''
  !
  ! Just scan all the counters, check each that it is reasonable
  ! and then report it
  !
  if(present(time_passed)) time_passed= zero_interval
  do i=1, n_counters
    if(.not. TimeCounters(i)%busy == silja_undefined)then
      strTmp = ''
      iTmp = i
      if(present(time_passed))then
        call report_time_counter(iTmp,strTmp,intervTmp)
        time_passed = time_passed + intervTmp
      else
        call report_time_counter(iTmp, strTmp)
      endif
    end if
  end do

end subroutine report_all_times


!*******************************************************************

subroutine report_time_counter(indCounter, chCounterNm, time_passed)
  !
  ! Stops the time counting and adds passed time to the interval.
  ! First index is checked, then name. At least one of them should point 
  ! to the right counter, otherwise nothing happens. Correct values
  ! for the counter are always returned if found
  !
  implicit none

  ! Imported parameters
  integer, intent(inout) :: indCounter
  character(len=*), intent(inout) :: chCounterNm
  type(silja_interval), intent(out), optional :: time_passed

  ! Local variables
  integer :: i

  ! If index is reasonable
  if(indCounter > 0 .and. indCounter < n_counters)then
    if(.not. TimeCounters(indCounter)%busy == silja_undefined)then
      chCounterNm = TimeCounters(indCounter)%chNm
      if(present(time_passed))time_passed = TimeCounters(indCounter)%time_passed
      call msg(fu_connect_strings('======>>>>  Time counter:',TimeCounters(indCounter)%chNm))
      call report(TimeCounters(indCounter)%time_passed)
      return
    else   ! Free unit
      call msg_warning('Sorry, counter is free','report_time')
      return
    endif   ! If counter is busy
  else   ! index is irrelevant, search for the name
    do i=1, n_counters
      if(TimeCounters(i)%chNm == chCounterNm)then
        indCounter = i
        if(.not.TimeCounters(indCounter)%busy == silja_undefined)then
          if(present(time_passed))time_passed = TimeCounters(i)%time_passed
          call msg(fu_connect_strings('======>>>>  Time counter:',TimeCounters(i)%chNm))
          call report(TimeCounters(i)%time_passed)
          return
        else
          call msg_warning('Sorry, counter is not counting','report_time')
          return
        endif
      endif
    end do  ! Search for the name
    call msg_warning('Sorry, counter is not found','report_time')
    return
  endif  ! If index is reasonable

end subroutine report_time_counter



!*********************************************************************************
!*********************************************************************************
!
! Progress meter
!
! Creates a file with progress percentage embedded in name and written inside
!
!*********************************************************************************
!*********************************************************************************

subroutine init_progress_file(chFNm)
  !
  ! Just cleans out the old progress files written during previous runs, if any
  !
  implicit none

  ! imported parameters
  character(len=*), intent(in) :: chFNm

  ! local vaiables
  integer :: iStat, uFile, iFNm, iNbrOfNames
  type(fnlen_str), DIMENSION(:), allocatable :: fnames
  
  if (smpi_global_rank /= 0) return ! Only one has to do it

  if (fu_str_u_case(trim(chFNm)) == 'SMS') then
    call msg("")
    call msg("")
    call msg("")
    call msg("")
    call msg("SMS is depricated!")
    call msg('Please replace "SMS" with "${SMSBIN}/smsmeter progress"')
    call msg_warning("Wrong progress_file parameter!","init_progress_file")
    return
  else if (fu_str_u_case(trim(chFNm(1:3))) == 'CMD') then
    call write_progress_file(chFNm, 0) 
  else
    uFile = fu_next_free_unit()
    !
    ! Find out all the files satisfying the template given in the chFNm. Note: progress has three 
    ! characters varying from 000 to 100
    !
    call string_to_fnames(chFNm + '_???', fnames, iNbrOfNames)
    do iFNm = 1, iNbrOfNames
      open(uFile, file=fnames(iFNm)%s, iostat = iStat)
      close(uFile,status='delete',iostat=iStat)
    enddo
  endif

end subroutine init_progress_file


!******************************************************************************

subroutine write_progress_file(chFNm, iProgr)
  !
  ! Takes the chFNm, adds the progress value and then writes the file with 
  ! this name and with progress inside
  !
  implicit none

  ! imported parameters
  character(len=*), intent(in) :: chFNm
  integer, intent(in) :: iProgr

  ! local vaiables
  integer :: iStat, iStat1, uFile, cnt
  CHARACTER(len=fnlen) :: envstr, mesg
  character(len=4), save :: chProgress=''
  integer, save :: lastprogr = -1

  if (iProgr == lastprogr) return
  lastprogr = iProgr

  if (smpi_global_rank /= 0) return ! Only one has to do it


  !check if we can do anything smarter than writing files
  if (fu_str_u_case(trim(chFNm)) == 'SMS') then
    call msg("SMS failed for progress", iProgr)
    call msg_warning('Please replace "SMS" with "${SMSBIN}/smsmeter progress"')

  elseif (fu_str_u_case(trim(chFNm(1:3))) == 'CMD') then
    envstr=trim(chFNm(4:))//' '//trim(fu_str(iProgr))
    call msg("Executing "//trim(envstr))
    iStat=0
    CALL EXECUTE_COMMAND_LINE(envstr, EXITSTAT=iStat, CMDSTAT=iStat1, CMDMSG=mesg)
    if (iStat1 /= 0) then
        call msg_warning("Progress command failed!", "write_progress_file")
        call msg(mesg) 
    elseif (iStat /= 0) then
       call msg_warning("Progress command returned: "//trim(fu_str(iStat)), "write_progress_file")
    endif
    return
  endif

  ! Not sms or sms failed for some reason....


  uFile = fu_next_free_unit()

  !
  ! Non-empty progress means that the previous file must be deleted
  !
  if(len_trim(chProgress) > 1)then
    open(uFile, file=trim(chFNm)//chProgress, iostat = iStat)
    close(uFile,status='delete',iostat=iStat)
  endif
  !
  ! Carefully fill-in the progress meter
  !
  if(iProgr > 999)then
    chProgress='____'
  else
    write(unit=chProgress,fmt='(A,I3.3)') '_',iProgr
  endif
  !
  ! Write the file
  !
  open(uFile, file=trim(chFNm)//chProgress, status='replace', iostat = iStat)
  if(iStat /= 0)then
    call set_error('Failed to open progress file: '//trim(chFNm)//chProgress, &
        & 'write_progress_file')
    return
  endif
  write(uFile,*) iProgr
  close(uFile)

end subroutine write_progress_file


!*********************************************************************************

integer function fu_compute_progress_time(time_start, time_end, time_now)
  !
  ! Computes the progress by compatring the time_now with start and end times
  !
  implicit none

  ! Imported parameters
  type(silja_time), intent(in) :: time_start, time_end, time_now

  fu_compute_progress_time = int(100.*((time_now-time_start)/(time_end-time_start)))
  if(fu_compute_progress_time < 0 .or. fu_compute_progress_time > 100) then
    call msg_warning('Strange progress index. time_start, time_end, time_now:')
    call msg('Strange progress index:',fu_compute_progress_time)
    call report(time_start)
    call report(time_end)
    call report(time_now)
    fu_compute_progress_time = min(max(fu_compute_progress_time,0),100)
  endif

end function fu_compute_progress_time


end module optimisation
