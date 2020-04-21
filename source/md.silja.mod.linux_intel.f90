MODULE md ! The machine dependent code of silja on Intel platform

  ! Description:
  ! This module contains all the code that might be machine-dependent. 
  ! 
  ! Current code owner: Mikhail Sofiev, FMI, email Mikhail.Sofiev@fmi.fi
  ! Large part of code are taken from other md units made by Mika Salonoja
  ! 
  ! All units: SI
  ! 
  ! Language: ANSI standard Fortran 90
  !
  ! Original code: Mika Salonoja
  ! Current code owner: Mikhail Sofiev, FMI mikhail.sofiev@fmi.fi
  !
  ! Modules used:

  use iso_c_binding, only : c_size_t
  use globals !, only : r4k, r8k
  use silam_mpi
  use revision
  use ifport
  use ifcore

  IMPLICIT NONE
  private
  ! The public functions and subroutines available in this module:

  PUBLIC fu_total_cpu_time_sec ! cpu-time used from the start
  PUBLIC wait
  PUBLIC shell_command ! executes one shell-level command in sh.
  public create_directory_tree
  public fu_current_directory
  public fu_pid
  public display_file_parameters
  public check_dir_slash
  public exit_with_status
  public open_binary_md
  public list_files
  public backtrace_md
  public system_mem_usage

  ! Local declarations:
  CHARACTER (LEN=1), PARAMETER, PUBLIC :: dir_slash = '/'

  integer, public, parameter :: sizeof_sp = 32 ! bits, as returned by sizeof

  public c_size_t
  public RENAME 
  public FTELL
  public revision_str

  CONTAINS 

  subroutine open_binary_md(unit, file, recl, status, action, convert, access, iostat)
    implicit none
    integer, intent(in) :: unit, recl
    character(len=*), intent(in) :: file, status, convert, action, access
    integer, intent(out) :: iostat

    !print *, 'Open binary:', unit, trim(file), recl, status, action, convert, iostat
    if (access == 'direct') then
      open(unit=unit, file=file, recl=recl, access='direct', form='binary', &
         & status=status, convert=convert, iostat=iostat)
    else
      open(unit=unit, file=file, access=access, form='unformatted', &
         & status=status, convert=convert, iostat=iostat)
    end if

  end subroutine open_binary_md

  ! *****************************************************************
  ! ****************************************************************
  ! 
  !
  !        CPU-time calculations
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  REAL FUNCTION fu_total_cpu_time_sec()
    !
    ! Returns the total cpu time [s] used by all the routines since the
    ! initalization of the time counter, or initialises if it zero.
        ! 
        ! Method:
        ! System function SECNDS return the number of seconds since the midnight
        ! minus the argument.
    !
    IMPLICIT NONE
 
  REAL(r4k), save :: start_sec = 0.0

    IF (start_sec > 0.001) THEN ! Positive initial value => already initialized
      fu_total_cpu_time_sec = SECNDS(start_sec)
    ELSE
      start_sec=0.
      start_sec = SECNDS(start_sec)  ! Initialisation
      fu_total_cpu_time_sec = 0.0
    END IF
  
  END FUNCTION fu_total_cpu_time_sec

  ! *****************************************************************

  SUBROUTINE init_time_counter()
    
    ! Description:
    ! Initializes the time counter by void call of fu_total_cpu_time_sec. 
    !
    ! Current code owner: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE
 
    REAL(KIND=r4k) t_int

    t_int = fu_total_cpu_time_sec()

  END SUBROUTINE init_time_counter


  !************************************************************************************

  SUBROUTINE wait(interval_sec)
    ! Description:
    ! Hibernates the program and waits a given interval.
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    REAL, INTENT(in) :: interval_sec ! in seconds

    call sleep(int(interval_sec))

  END SUBROUTINE wait


  ! ***************************************************************
  ! ***************************************************************
  !
  !
  !         Execute shell level commands.
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  SUBROUTINE shell_command(command_string, command_string_2)

    ! Description:
    ! Executs one shell-level command in sh. If an error occurs in
    ! shell, and error is set.
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    CHARACTER (LEN=*), INTENT(in) :: command_string

    ! Imported parameters with intent IN:
    CHARACTER (LEN=*), INTENT(in), OPTIONAL :: command_string_2

    ! Local declarations:
    INTEGER :: istat !, ishell
    CHARACTER (LEN=fnlen) :: cmd, strTmp1, strTmp2

    cmd = ' '
    strTmp1 = ADJUSTL(command_string)

    IF (PRESENT(command_string_2)) THEN
      strTmp2 = ADJUSTL(command_string_2)
      WRITE(unit = cmd, fmt = '(A,A,A)')&
          & TRIM(strTmp1),&
          & ' ', &
          & TRIM(strTmp2)
    ELSE
      cmd = TRIM(strTmp1)
    END IF

    call msg_test(fu_connect_strings('executing:', trim(cmd)))

    istat = SYSTEM(cmd)

    IF (istat < 0) THEN
      CALL set_error(fu_connect_strings('error while executing:',&
          & cmd), 'shell_command')
    END IF

  END SUBROUTINE shell_command

  !************************************************************************************

  subroutine list_files(pattern, file_list_name)
  !
    ! Execute a shell command to list all files matching the given pattern and direct the
    ! output into a given file.
    !
    implicit none
    character(len=*), intent(in) :: pattern, file_list_name
    character(len=*), parameter :: sub_name = 'list_files'

    character(len=1000) :: cmd
  
    if (fu_fails(file_list_name /= '', 'Empty file_list_name', sub_name)) return
    if (fu_fails(pattern /= '', 'Empty pattern', sub_name)) return
  
    cmd = ' '
    WRITE(unit = cmd, fmt = '(4A)') 'ls -1 ', trim(pattern), '>', file_list_name

    CALL shell_command(cmd)

  end subroutine list_files

  !*****************************************************************

  function fu_current_directory() result(dir)
    !
    ! Returns the current directory
    !
    implicit none

    ! Return value
    character(len=fnlen) :: dir
    integer :: stat

    stat = getcwd(dir)
    if (stat /= 0) then
       call set_error('Error with getcwd','fu_current_directory')
      dir = ''
    endif

  end function fu_current_directory

  !*******************************************************************

  subroutine create_directory_tree(chDir)
    !
    ! Creates the directory tree if it does not exist.
    !
    implicit none

    ! Imported parameter
    character(len=*), intent(in) :: chDir

    ! Local variables
    logical :: res_make_dir
    integer :: pos_slash, posTmp, iStart, systat

    ! Nothing to create
    if (len(chDir) .le. 1) return
    
    !
    ! Try simply to create what is requested
    !
    call msg(fu_connect_strings('Directory to create:',chDir,'*********'))
    
    systat = system('mkdir -p '//chDir)
    
    if(systat /= 0)then
       call msg_warning('Could not create directory: '//chDir)
    end if

  end subroutine create_directory_tree

  !************************************************************************************

  subroutine display_file_parameters(chFNm)
    !
    ! Prints file information to stdout and log.
    !
    implicit none

    ! Argument
    character(len=fnlen), intent(in) :: chFNm
    
    ! Local variables
    integer :: unit, filestat
    logical :: is_open
    integer(4), dimension(13) :: stats
    character(len=fnlen) :: fn_full
    character(len=30) :: date

    if(chfnm /= '')then
      do unit = 50, 1000
        inquire(unit=unit, opened=is_open)
        if (.not. is_open) exit
      end do
      if (is_open) then
        call set_error('Too many open files (>1000)', 'display_file_parameters')
          return
        endif

      open(unit, file=chfnm, action='read')
      inquire(unit, name=fn_full)
      filestat = fstat(unit, stats)
      close(unit)
      
      call msg('File name: ' // trim(fn_full))
      call msg('File size, kBytes:', real(stats(8)/1000.0))
      
      call msg('Last modification: ' // ctime(stats(10)))
      call msg('Last access: ' // ctime(stats(9)))
        else
      call msg_warning('Empty file name  given','display_file_parameters')
    endif

  end subroutine display_file_parameters

  !************************************************************************************

  integer function fu_pid()
    implicit none
    !
    ! Return the process ID.
    !
    fu_pid = getpid()
  end function fu_pid

  !*****************************************************************

  subroutine check_dir_slash(chFNm)
    ! Silly windows trick that does absolutely nothing here.
    implicit none
    
    character(len=*), intent(inout) :: chFNm

  end subroutine check_dir_slash

  subroutine exit_with_status(status)
    implicit none
    integer, intent(in), optional :: status
    if (present(status)) then
      call exit(status)
    else
      call exit(0)
    end if
  end subroutine exit_with_status

  subroutine backtrace_md()
    implicit none
      call tracebackqq(user_exit_code=-1)
  end subroutine backtrace_md
  

subroutine system_mem_usage(valueRSS)
        !
        ! Get valueRSS of the current process (in kB)
        ! stolen (with some fixes) from
        ! https://stackoverflow.com/questions/22028571/track-memory-usage-in-fortran-90
        implicit none
        integer, intent(out) :: valueRSS

        character(len=200):: filename=' '
        character(len=80) :: line
        character(len=8)  :: pid_char=' '
        integer :: iunit, stat
       

        valueRSS=-1    ! return negative number if not found

        !--- get proc filename

        write(filename,'(A,I0,A)') '/proc/',getpid(),'/status'
        iunit = 1001 ! something never returned by fu_next_free_unit()

        !--- read system file

        open(unit=iunit, file=filename, action='read', IOSTAT=stat)
        if (stat == 0) then
          do
            read (iunit,'(a)',end=120) line
            if (line(1:6).eq.'VmRSS:') then
               read (line(7:),*) valueRSS
               exit
            endif
          enddo
          120 continue
          close(iunit)
        else
          call msg("Faileled to open:"//trim(filename)//' status:',stat)
        endif

        return
end subroutine system_mem_usage



END MODULE md

