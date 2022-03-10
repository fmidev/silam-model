MODULE md ! The machine dependent code of silja on Intel platform

  ! Description:
  ! This module contains all the code that might be machine
  ! -dependent. Such things are memory-sizes, and calls to performacs
  ! -measuring libraries etc.
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

  use ISO_C_BINDING, only : c_size_t
#ifdef VS2012
  use psapi
#endif
  use globals !, only : r4k, r8k
  use revision
#ifndef __GFORTRAN__
  use ifport
  use ifcore
#endif
  
  IMPLICIT NONE
  private
  ! The public functions and subroutines available in this module:

  PUBLIC fu_total_cpu_time_sec ! cpu-time used from the start
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
  public fu_system_mem_usage
  public fseek_md
  public fu_irand
  public get_hostname

  ! Local declarations:

#ifdef VS2012 
  CHARACTER (LEN=1), PARAMETER, PUBLIC :: dir_slash = '\' ! 
#else
  CHARACTER (LEN=1), PARAMETER, PUBLIC :: dir_slash = '/' ! 
#endif

  integer, public, parameter :: sizeof_sp = 32 ! bits, as returned by sizeof

  public c_size_t
#ifndef __GFORTRAN__
  public RENAME !FROM IFPORT
  public FTELL  !FROM IFPORT
#endif
  public revision_str

#ifdef VS2012 
  !
  ! Define a few variables that are in the NETcdf4 library but not in the netcdf3 one
  ! To be removed when netcdf4 made through to Windows
  !
!  integer, parameter, public :: NF90_USHORT = int_missing, NF90_UINT = int_missing+1, &
!                              & NF90_INT64 = int_missing+2, NF90_UINT64 = int_missing+3, &
!                              & NF90_NETCDF4 = int_missing+4

!  interface
!  subroutine GetCurrentProces (proc) bind (c, NAME='GetCurrentProcess')
!      use, intrinsic :: ISO_C_BINDING
!      integer(c_int), value :: proc
!  end subroutine
!  end interface
#else

  interface
    integer(c_int) function getpagesize () bind (c, NAME='getpagesize')
      use, intrinsic :: ISO_C_BINDING
    end function
  end interface

#endif
  
  
  CONTAINS 

  subroutine open_binary_md(unit, file, recl, status, action, convert, access, iostat)
    implicit none
    integer, intent(in) :: unit, recl
    character(len=*), intent(in) :: file, status, convert, action, access
    integer, intent(out) :: iostat

#ifdef VS2012
!! Non-standard
!!! https://community.intel.com/t5/Intel-Fortran-Compiler/Difference-between-form-unformatted-and-form-binary/td-p/1108571
    character (len=*), parameter :: form='binary'
#else
    character (len=*), parameter :: form='unformatted'
#endif

    !print *, 'Open binary:', unit, trim(file), recl, status, action, convert, iostat
    if (access == 'direct') then
      open(unit=unit, file=file, recl=recl, access='direct', form=form, &
         & status=status, convert=convert, iostat=iostat)
    else
      open(unit=unit, file=file, access=access, form='unformatted', &
         & status=status, convert=convert, iostat=iostat)
    end if

  end subroutine open_binary_md

!*****************************************************************


  subroutine get_hostname(hostname)
    implicit none
    character(len=*), intent(out) :: hostname
    integer :: iSt
    
    iSt = 0
#ifdef __GFORTRAN__
!! Could be some other version
#if    __GNUC__ > 4 || __GNUC__ == 4 &&  __GNUC_MINOR__ > 7   
      call HOSTNM(hostname,STATUS=iSt)
#else
      call msg('(HOSTNM not available in this compiler (?)')
#endif
#else
    iSt = HOSTNAM (hostname)
#endif
    if (iSt /= 0) then
        call set_error("failed to get hostname, status:"//trim(fu_str(iSt)), 'get_hostname')
    endif


  end subroutine get_hostname

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

    call msg_test('executing:' + trim(cmd))

    istat = SYSTEM(cmd)

    IF (istat < 0) THEN
      CALL set_error('error while executing:' + cmd, 'shell_command')
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
    
    character(len=65535) :: cmd
#ifdef VS2012
    character(len=*), parameter :: ls='dir /B /S '
#else    
    character(len=*), parameter :: ls='ls -1 '
#endif

    if (fu_fails(file_list_name /= '', 'Empty file_list_name', sub_name)) return
    if (fu_fails(pattern /= '', 'Empty pattern', sub_name)) return
  
    cmd = ' '
    WRITE(unit = cmd, fmt = '(4A)') ls, trim(pattern), '>', trim(file_list_name)
  
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
    integer :: length, stat

#ifdef VS2012
    dir = FILE$CURDRIVE  ! Current drive
    length = getDriveDirQQ(dir)
    if (length == 0) then
#else
    stat = getcwd(dir)
    if (stat /= 0) then
#endif
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

#ifdef VS2012
    ! Local variables
    logical :: res_make_dir
    integer :: pos_slash, posTmp, iStart, systat
    !
    ! First, try simply to create what is requested
    !
    res_make_dir = makedirqq(chDir)

    ! If failed, let's create one-by-one the whole tree, requiring the last
    ! makedirqq command to be successfull
    !
    if(res_make_dir)then
      return
    else
!    GETLASTERRORQQ to retrieve the error message. Possible errors include: !
!
      !ERR$ACCESS: The directory was not created. The given name is the name of an
      !existing file, directory, or device.  ERR$NOENT:

      iStart = getlasterrorqq()
      if(iStart == ERR$ACCES .or. iStart == ERR$EXIST)then
        call msg_warning(fu_connect_strings('Directory already exists:',chDir), &
                       & 'create_directory_tree')
        return
      elseif(iStart == ERR$NOENT)then
        !
        ! If there is no sub-directories in the requested path - set an error
        !
        if(index(chDir,dir_slash) == 0)then
          call set_error(fu_connect_strings('Failed to create directory:',chDir), &
                       & 'create_directory_tree')
          return
        endif
        !
        ! If there is a leading slash or disk name - skip them
        !
        posTmp=0
        if(index(chDir(1:3),':') > 0 .or. index(chDir(1:1),dir_slash) > 0)then
          pos_slash = index(chDir(1:len_trim(chDir)),dir_slash)+1
        else
          pos_slash = 1
        endif
        !
        ! Add directories one-by-one trying to create them
        !
        do while (pos_slash > 0)
          posTmp = pos_slash + posTmp
          pos_slash = index(chDir(posTmp+1:len(chDir)),dir_slash)
          if(pos_slash > 0) res_make_dir = makedirqq(chDir(1:posTmp + pos_slash-1))
        end do
        !
        ! Check trailing dir slash - if it exists, the pathe is already created.
        ! If not - the last sub-directory is not yet made
        !
        if(index(chDir(len_trim(chDir):len_trim(chDir)),dir_slash) == 0) &
                                          & res_make_dir = makedirqq(chDir)
        if(.not.res_make_dir)then
          call set_error(fu_connect_strings('Failed to create directory:',chDir), &
                       & 'create_directory_tree')
        endif
      elseif(iStart == ERR$NOMEM)then
        call set_error('No memory error happened','create_directory_tree')
      else
        call set_error('Unknown error. Failed to make the tree:' + chDir,'create_directory_tree')
      endif
    end if
#else
    ! Local variables
    logical :: exists
    integer :: pos_slash, posTmp, iStart, EXITSTAT, CMDSTAT
    integer :: itry
    character(len=fnlen) :: CMDMSG
    integer, parameter :: maxtries = 10
    character(len=*), parameter :: sub_name = 'create_directory_tree'

    ! Nothing to create
    if (chDir == '') return

    if (chDir(1:1) == '-') then
       call set_error('Directories, staring with "-" disabled. chDir = '//trim(chDir), sub_name)
       return
    endif

    INQUIRE (file=chDir, EXIST=exists)
    !GNU fortran does not support directory keyword
    ! If the _file_ exists with the name we are doomed anyway..
    if (exists) return !!Calling external programs is risky at voima, so better check...
     
    
    !
    ! Try simply to create what is requested
    !
    do itry=1,maxtries
      call msg('Creating directory: "'+chDir+'"')

      CALL EXECUTE_COMMAND_LINE('/bin/mkdir -p '//trim(chDir), wait=.TRUE.,EXITSTAT= EXITSTAT, &
           & CMDSTAT=CMDSTAT, CMDMSG=CMDMSG )
      
      !!!systat = system('/bin/mkdir -p '//trim(chDir))

      if(EXITSTAT == 0 .and. CMDSTAT == 0) exit

      call msg("EXITSTAT, CMDSTAT", EXITSTAT, CMDSTAT)
      call msg("CMDMSG:" // CMDMSG)
      call msg("len", len(chDir))
      call msg("mkdir -p  failed")

      call msg("Attempt "//trim(fu_str(itry))//" failed. Sleeping 5s")
      call sleep(5)
      if (itry == maxtries) then
        call set_error('Could not create directory: '//trim(chDir), sub_name)
      endif
      
     enddo

#endif
  end subroutine create_directory_tree

  !************************************************************************************

  subroutine display_file_parameters(chFNm)
    !
    ! Prints file information to stdout and log.
    !
    implicit none

    ! Argument
    character(len=*), intent(in) :: chFNm

#ifdef VS2012
    ! Local variables
    integer :: unit, filestat, iTmp
    integer (2) :: iyr, imon, iday, ihr, imin, isec
    integer(8) :: handle
    logical :: is_open
    integer(4), dimension(13) :: stats
    character(len=fnlen) :: fn_full, chTmp
    character(len=30) :: date
    type(FILE$INFOI8) :: info

    if(len_trim(chFNm) > 0)then
      handle = FILE$FIRST
      iTmp = GETFILEINFOQQ(chFNm, info, handle)

      call msg('File name:' + chFNm)
      call msg('File size, kBytes:', real(info%length)/1000.0)
      
      if(info%creation > 0)then
        CALL UNPACKTIMEQQ (info%creation, iyr, imon, iday, ihr, imin, isec)
        WRITE(unit=chTMp, fmt ='(A,1x,I4, A, I2, A, I2, 2x, I2, A, I2, A, I2)') &
             & 'Creation time: ', iyr, '.',imon,'.', iday, ihr,':', imin,':', isec
        call msg(chTmp)
      endif
      
      CALL UNPACKTIMEQQ (info%lastwrite, iyr, imon, iday, ihr, imin, isec)
      WRITE(unit=chTMp, fmt ='(A,1x,I4, A, I2, A, I2, 2x, I2, A, I2, A, I2)') &
           & 'Last write: ', iyr, '.',imon,'.', iday, ihr,':', imin,':', isec
      call msg(chTmp)
      
      if(info%lastaccess > 0)then
        CALL UNPACKTIMEQQ (info%lastaccess, iyr, imon, iday, ihr, imin, isec)
        WRITE(unit=chTMp, fmt ='(A,1x,I4, A, I2, A, I2, 2x, I2, A, I2, A, I2)') &
                         & 'Last access: ', iyr, '.',imon,'.', iday, ihr,':', imin,':', isec
        call msg(chTmp)
      endif
      
    else
      call msg_warning('Empty file name given','display_file_parameters')
    endif
#else

    !!! This should work also on Intel fortran with ifport 
    !!! (probably in Windows as well)
    ! Local variables
    integer(4), dimension(13) :: stats !! Can also be integer8 in intel, but must integer4 in GNU
    integer :: istat
    character(len=30) :: date


      istat =  STAT(chfnm, stats)

      if (istat /= 0) then
        call msg('Failed to get binary status binary by name: "'// trim(chfnm)//'"')
        return
      endif
      
      
      call msg('File name: ' // trim(chfnm))
      call msg('File size, kBytes:', real(stats(8)/1000.0))
      call ctime(int(stats(10), 8), date)
      call msg('Last modification: ' // date)
      call ctime(int(stats(9), 8), date)
      call msg('Last access: ' // date)

#endif
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
    !
    ! Scans the directory searching for backslash (/)
    ! and replacing it forward (\).
    !
    implicit none

    character(len=*), intent(inout) :: chFNm

#ifdef VS2012    
    integer :: iTmp
    
    do iTmp = 1,len_trim(chFNm)
      if(chFNm(iTmp:iTmp) == '/') chFNm(iTmp:iTmp) = dir_slash
    enddo
#endif
    
  end subroutine check_dir_slash

  subroutine exit_with_status(status)
    implicit none
    integer, intent(in), optional :: status
#ifdef VS2012    
    if (present(status)) then
      stop "Status given" !call  exit(status)
    else
      stop ! call exit(0)
    end if
#endif
  end subroutine exit_with_status

  !*****************************************************************

  subroutine backtrace_md()
    implicit none
#ifdef __GFORTRAN__
#if    __GNUC__ > 4 || __GNUC__ == 4 &&  __GNUC_MINOR__ > 7
      call backtrace()
#else
      call msg('(backtrace() not available in this compiler)')
#endif
#else
    call tracebackqq(user_exit_code=-1)
#endif
  end subroutine backtrace_md
  

  !***********************************************************************
  
  integer function fu_system_mem_usage()result(mem_usage) ! valueRSS)
    !
    ! Get valueRSS of the current process (in kB)
    !
    implicit none

    character(len=200):: filename=' '
    character(len=80) :: line
    integer :: iunit, stat
    integer (kind=8) size, rss, share, text, lib, data, dt, my_pagesize
      
#ifdef VS2012    
    character(len=8)  :: pid_char=' '
    integer(HANDLE) :: hProcess
    type (T_PROCESS_MEMORY_COUNTERS_EX) :: pmc
    integer(BOOL) :: ret
    integer :: iTmp

    mem_usage = 0
    ! Print information about the memory usage of the process

    hProcess = NULL
    
!    call GetCurrentProcess(hProcess)
    if (hProcess == NULL) return
    
    iTmp = sizeof(pmc)
    ret = GetProcessMemoryInfoEx (hProcess, pmc, iTmp)

!    print *, pmc
!    stop
    print *, mem_usage
        !!! Something smarter can be invented here
#else
    character(len=*), parameter  :: statm='/proc/self/statm'
    ! Stolen from slurm source code (seems to be the way Slurm counts memory)
    mem_usage = -1    ! return negative number if not found
    iunit = 1001 ! something never returned by fu_next_free_unit()

    open(unit=iunit, file=statm, action='read', IOSTAT=stat)
    if (stat == 0) then
        read (iunit,*) size, rss, share, text, lib, data, dt
        close(iunit)
        my_pagesize = getpagesize()
        mem_usage = int( (rss) * my_pagesize / 1024_8, kind=4) !!!To kB
    else
      call msg("Faileled to open: "//statm//'  status:',stat)
    endif

#endif

end function fu_system_mem_usage

!****************************************************************************  
  
subroutine fseek_md(iUnit, iOff,  from)
  !
  ! A wrapper for FSEEK function/subroutine
  !
  implicit none
  
  ! Imported parameters
  integer, intent(in) :: iUnit, iOff,  from
  
  ! Local variab;es
  integer :: iStat

#ifdef VS2012  
  iStat = FSEEK(iUnit, iOff,  from)
#else
  call FSEEK(iUnit, iOff,  from, iStat)
#endif
  if(iStat == 0)then
    return
  else
    call set_error('FSEEK failed, unit=' + fu_str(iUnit),'fseek_md')
  endif
  
end subroutine fseek_md


!*****************************************************************************
  
integer function fu_irand()
  !
  ! A wrapper for IRAND, which in Intel compiler sits in IFPORT, at least in Windows
  !
  fu_irand = irand()
end function fu_irand


END MODULE md

