PROGRAM silam_main

  ! Description: 
  ! Silam version 5.3.1 model main program
  !
  ! Input information:
  ! One input file name
  ! Options for the input arguments:
  ! - No arguments. File 'control.ini' is read, the output directory name is time-made
  ! - 1 argument. Treated as an ini file name instead of 'control.ini'
  !
  ! Possible output files:
  ! GRIB file, GrADS file, NetCDF file, Trajectory file
  !
  ! ATTENTION. Model is memory-hungry. It requires the stack size of at least 30,000,000 reserved
  !            and 25,000,000 committed.
  !
  ! Author: Mikhail Sofiev, FMI.
  ! email: Mikhail.Sofiev@fmi.fi
  !
  ! Language: ANSI Fortran 90
  !
  ! All units: SI
 
  USE dispersion_models
  use optimisation
  use silam_namelist
  use silam_times
  use revision
  use svndiff
  use silam_partitioning 
  use field_identifications
  use ascii_io

  use, intrinsic :: iso_fortran_env, only: compiler_version, compiler_options
!  use source_terms_wind_blown_dust
  
!  use source_term_fires
use tangent_linear
  
  
  IMPLICIT NONE

  CHARACTER (LEN=fnlen) :: control_fname, hostname, tmpfile
  character(len=clen) :: chTmp, rankTmp 
  real :: fTmp, fTmp2
  integer :: iStatus, i,j, iTmp, iz, it
  logical :: had_error
  
  integer, parameter :: exit_ok = 0, exit_bad = 9
  integer :: ierr

  real, dimension(:), pointer :: x, y
  integer :: nP, nMissing, MyPID
  real :: slope, intercept, r_value, p_value, slope_err, stderr
  !  real, dimension(2) :: fNbrFlux, fVolFlux, fMassMeanDiam
  type(Taerosol_mode), dimension(50) :: aerModes

  
    type(grads_template) :: templ
    type(silja_time) :: valid_time
    type(silam_sp), dimension(:), pointer :: fnames ! inout
    integer :: nbrOfFiles
    type(silja_time) :: anal_time 
    logical :: ifStrict = .false., ifAdd = .false. ,ifAllowZeroFcLen = .false. , ifWait = .false.
    type(silja_interval) :: fc_step = one_hour
    type(silja_interval) :: max_hole_length = zero_interval   ! if some time slot is missing
  
  
#ifdef __GFORTRAN__
    !! Updated with /usr/include/x86_64-linux-gnu/bits/signum-generic.h
     integer, parameter :: SIGHUP = 1
     integer, parameter :: SIGTERM = 15
     integer, parameter :: SIGUSR2 = 31
     integer, parameter :: SIGUSR1 = 30
     call signal(SIGUSR1, die_sigusr1)
     call signal(SIGUSR2,  die_sigusr2)
     call signal(SIGHUP,  die_sighup)
     call signal(SIGTERM,  die_sigterm)
#endif



!   call test_advect_mass_many()  !test_advect_mass()
!  run_log_name = "/dev/null"
!  run_log_funit = fu_next_free_unit()
!  open(run_log_funit, file=run_log_name, iostat = iStatus)
!  if(iStatus /= 0)then  ! The problem is serious, randomised file name does not help
!    call set_error('Failed to open: "'//trim(run_log_name)//'"','grads_2_grib2_main')
!  endif
!   call test_sort()
!   stop

  call smpi_init()

!  call init_random_seed()
!
!  call msg('Allocating...')
!  allocate(arPtr(36000,18000),stat=i)
!  call msg('status:',i)
!  do i=1,36000
!    do j = 1, 18000
!      arPtr(i,j)= fu_random_number_center(90.,10.)
!      if(mod(i,2000) ==0 .and. mod(j,2000) ==0) call msg('val:',arPtr(i,j))
!    enddo
!  enddo
!  
!stop
  
  !
  ! Start the global time counter
  !
  chTmp = 'Overall_run_time'
  call start_count(chCounterNm = chTmp)
  call init_random_seed()

  !----------------------------------------
  !
  !  Open the global error and warning recording file
  !  This very file will be used for error messages, warnings and, later,
  !  for other log functions
  !  However, when the output directory becomes known - it will be transferred to it
  !  and its number will be changed
  !
  MyPID = fu_pid()
  call get_hostname(hostname) !!! Gotcha! PID is not unique among nodes in slurm
  if (smpi_global_rank == 0) then
    write(unit=chTmp,fmt='(I8.8,A)') MyPID, '_'//trim(hostname)
  endif
  call smpi_bcast_string(chTmp, MPI_COMM_WORLD) !!set it from master
  !! Global variable
  write(unit=proc_ID_string,fmt='(A,A1,I3.3)') trim(chTmp), '_', smpi_global_rank

  tmpfile = "run_tmp_"// trim(proc_ID_string) // '.log'
  run_log_funit = fu_next_free_unit()
  open(run_log_funit, file=tmpfile, iostat = iStatus)
  if(iStatus /= 0)then
    call set_error(fu_connect_strings('Failed to open:',chTmp),'silam_main')
    call smpi_finalize()
    stop
  endif

  
  call smpi_global_barrier()

  if (smpi_global_rank == 0) then
    PRINT * ;   PRINT * ;   PRINT * ;
  endif

  CALL msg ('Hello world! This is ' // revision_str // ' speaking. PID = '//trim(fu_str(MyPID)))
  call msg('Running at HostName: '//trim(hostname))
  call msg('compiler_version: '// compiler_version())
  call msg('compiler_options: '// compiler_options())
  call msg('Log file name (on start): '//trim(tmpfile))
  
!  allocate(x(200))
!  open(20,file='d:/data/emis/fakes/test_area_src.grads',form='binary',recl=4*30,access='direct')
!  close(20,status='delete')
!  open(20,file='d:/data/emis/fakes/test_area_src.grads',form='binary',recl=4*30,access='direct')
!
!  iStatus = 1
!  do it = 1, 100
!    do iz = 1, 3
!      do i = 1, 20
!        do j = 1, 30
!          x(j) = (mod(i,5) * 30 + j) * iz * min(max(0.0, sin(it * 0.1) + 0.8), 1.5)
!!          x(j) = (i * 300 + j) * iz * it
!!          x(j) = (i + j) * iz * it
!        end do
!        write(20,rec=iStatus) x(1:30)
!        iStatus = iStatus + 1
!      end do
!    end do
!  end do
!  close (20)
!  stop
  
  
#ifndef SUPRESS_SVN_DIFF
  call print_svn_diff()
#endif

  if(smpi_is_mpi_version())then
     write(unit=rankTmp,fmt='(A25,I3,A6)') 'MPI version running with ', smpi_global_tasks, ' tasks'
  else
     write(unit=rankTmp,fmt='(A14)') 'Serial version'
  endif
  call msg(rankTmp)
  call msg('SILAM executable parameters:')
  !
  ! Display the current executable parameters
  !
  call get_command_argument(0, control_fname)    ! temporary use of the variable
  call display_file_parameters(control_fname)
  call msg('')
  call msg("Current directory:"+fu_current_directory())
  call msg('')
  !
  ! Time marks
  !
!  call msg('Local time now: ' + fu_time_string_original(fu_wallclock()))
  call msg('Local time now: ' + fu_computer_time_string())
  call msg('UTC time now: ' + fu_str(fu_wallclock()))
  call msg('Memory usage, kB',fu_system_mem_usage())
  !
  ! Get the control file name
  !
  SELECT CASE (command_argument_count())
    CASE(0) !--------- No arguments, read all from the default ini file
      call msg('No arguments - use default silam.ini')
      IF(.not.read_ini_file('silam.ini'))THEN
        CALL set_error('Can not read silam.ini','silam_main')
        CALL smpi_abort(exit_bad)
      END IF

    CASE(1) !-------- One argument - ini file name, read from this ini file
      call msg('1 argument - use it as control file name')
      CALL get_command_argument(1, control_fname)

    CASE DEFAULT
      CALL set_error('Too many input arguments','silam_main')
  END SELECT

  call msg('Control file name: '// trim(control_fname))!!Space must be here!
  call msg('')

  CALL start_silam_v5(control_fname, had_error)

  call msg("At the end of the run memory usage (kB)", fu_system_mem_usage())
  
  chTmp = 'Overall_run_time'
  call stop_count(chCounterNm = chTmp)
  call report_time()

  ! Before renaming the file, it must be closed (Intel compiler requirement)
  !
  
  call msg('Local time now: ' // fu_computer_time_string())
  call msg('UTC time now: ' // fu_str(fu_wallclock()))
  
  if ( had_error) then
     call set_error("HAD Errors, not renaming log file",'silam_main')
     call msg("Keeping log file: "//trim(run_log_tmp_name) )
     close(run_log_funit)
  else
    close(run_log_funit)
    if ( RENAME( run_log_tmp_name, run_log_name) /=0) then
      open(run_log_funit, file = run_log_tmp_name, status='old', position='append')
      call msg("TMP    file: "//trim(run_log_tmp_name))
      call msg("target file: "//trim(run_log_name))
      call set_error("Rename log file tmp failed",'silam_main')
      close(run_log_funit)
    else
     print *, "Log file: "//trim(run_log_name) 
     close(run_log_funit)
    endif
  endif

  if (smpi_is_mpi_version()) then
    if (had_error) then
      call sleep(5) !! Let others a chance to finsh
      call smpi_abort(exit_bad)
    else
      call smpi_finalize()
    end if
  else if (had_error) then
    call exit_with_status(exit_bad)
  else
    call exit_with_status(exit_ok)
  end if
  
  CONTAINS

  !----------------------------------------------------------------------------

  LOGICAL FUNCTION read_ini_file(ini_fnm)
    !
    ! Just reads the ini file
    !
    IMPLICIT NONE

    ! Import parameters with intent IN
    CHARACTER(len=*), INTENT(in) :: ini_fnm

    ! Local stuff
    INTEGER :: unit, status
    type(Tsilam_namelist), pointer :: nlIni
    logical :: lExist

      read_ini_file = .false.

      unit = fu_next_free_unit()

      OPEN(unit, file=ini_fnm, action='read', status='old', iostat=status)

      IF(fu_fails(status == 0, 'Ini file does not exist:' + ini_fnm, 'read_ini_file')) return

      nlIni => fu_read_namelist(unit, .true.) ! Just two items in the namelist
      close(unit)

      !
      ! Get the control file name and check that it exists
      !
      control_fname = fu_content(nlIni,'control_file')
      inquire(file=control_fname, exist = lExist)
      if(fu_fails(lExist,'Control file does not exist:' + control_fname, 'read_ini_file')) return

      call destroy_namelist(nlIni)

      read_ini_file = .true.

  END FUNCTION read_ini_file


 subroutine die_gracefully(reason)
    character(len=*), intent(in) :: reason
     !$omp master
     if (smpi_global_tasks > 1) then 
       call msg(trim(reason)//", sleeping before setting error. rank=", smpi_global_rank)
       call sleep(smpi_global_rank*20/smpi_global_tasks) !no more than 10s in totoal
       ! Slurm waits up to 32 seconds on timeout after sighup before killing 
       write (6, '(a,x,i4)') "backtrace from rank",smpi_global_rank
       flush (6) !
       call backtrace_md()
     endif
     flush (run_log_funit) !make sure that last message is there
     call set_error(reason,"warning_sigusr1")
     flush (run_log_funit) !make sure that the error messgae is there
     !$omp end master
 end subroutine die_gracefully

 subroutine die_sigusr1
   call die_gracefully("Got sigusr1")
 end subroutine die_sigusr1

 subroutine die_sigusr2
   call die_gracefully("Got sigusr2")
 end subroutine die_sigusr2
    
 subroutine die_sigterm
   call die_gracefully("Got sigterm")
 end subroutine die_sigterm
    
 subroutine die_sighup
   call die_gracefully("Got sighup")
 end subroutine die_sighup
    
END PROGRAM silam_main

