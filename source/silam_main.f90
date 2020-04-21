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
!  use source_terms_wind_blown_dust
  
!  use source_term_fires
  
  
  IMPLICIT NONE

  CHARACTER (LEN=fnlen) :: control_fname
  character(len=clen) :: chTmp, rankTmp
  real :: fTmp, fTmp2
  integer :: iStatus, i,j, iTmp, iz, it
  logical :: had_error
  
  integer, parameter :: exit_ok = 0, exit_bad = 9
  integer :: ierr

  real, dimension(:), pointer :: x, y
  integer :: nP, nMissing
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
     integer, parameter :: SIGUSR1 = 31 
     integer, parameter :: SIGINT = 15
     call signal(SIGUSR1, warning_sigusr1)
     call signal(SIGINT, warning_sigint)
#endif



!   call test_advect_mass_many()  !test_advect_mass()
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
  run_log_funit = fu_next_free_unit()
  !call random_number(fTmp)
  !if(fTmp < 0.1) fTmp = fTmp + 0.1
  call smpi_allreduce_max_int(fu_pid(), iTmp, MPI_COMM_WORLD)
  write(unit=chTmp,fmt='(A8,I8.8,A1,I3.3,A4)') 'run_tmp_', iTmp, '_', smpi_global_rank, '.log'
  open(run_log_funit, file=chTmp, iostat = iStatus)
  if(iStatus /= 0)then
    call set_error(fu_connect_strings('Failed to open:',chTmp),'silam_main')
    call smpi_finalize()
    stop
  endif

  
!!!  
!!!  
!!!
!!!  anal_time = ref_time_01012000
!!!  valid_time = anal_time + one_hour
!!!  call decode_template_string('d:\model\silam_v5_7\ini\non_existing_%y4.ini', templ)
!!!  ifStrict = .true.
!!!  ifAdd = .false.
!!!  ifAllowZeroFcLen = .true.
!!!  ifWait = .false.
!!!  nullify(fnames)
!!!  do i = 1, 100000
!!!    call msg('Memory usage before the call', fu_system_mem_usage()) 
!!!    call fnm_from_single_template(templ, valid_time, fnames, &
!!!         & ifStrict = .true., &
!!!         & ifadd = .false., &
!!!         & ifWait = .false., &
!!!         & ifAllowZeroFcLen = .true.)
!!!
!!!!    call FNm_from_single_template(templ, &
!!!!                              & valid_time,fnames, &  !nbrOfFiles, anal_time, &
!!!!                              & ifStrict,ifAdd,ifAllowZeroFcLen,ifWait)  !, &
!!!!!                              & fc_step, max_hole_length)
!!!  
!!!  call msg('nUMBER_OF_FILES', nbrOfFiles)
!!!  call msg('Memory usage after the call', fu_system_mem_usage()) 
!!!  enddo
!!!  stop
!!!  
  
  
  
  
  
  
  
!!!  
!!!  
!!!  fTmp = 1e-9
!!!  do i = 1, 50
!!!    aerModes(i) = fu_set_mode(fixed_diameter_flag, fTmp, 1.3*fTmp, 1.2*fTmp)
!!!    print *, fTmp
!!!    
!!!    call fires_flux4mode(aerModes(i), &        ! definition of the spectrum band
!!!!                                 & 6120, &  !private lognorm_3_modes_dust__analytic_flag, &
!!!                                 & 6121, &  !private lognorm_3_modes_dust__numeric_flag, &
!!!                                 & fNbrFlux, fVolFlux, fMassMeanDiam)
!!!    call msg('Mode min,mean,max diam, [um], fluxNbr, fluxVol, meanActualDiam: ', &
!!!           & (/fu_min_d(aerModes(i)), fu_mean_d(aerModes(i)), fu_max_d(aerModes(i)), &
!!!             & fNbrFlux, fVolFlux, fMassMeanDiam, &
!!!             & fNbrFlux / (fu_max_d(aerModes(i)) - fu_min_d(aerModes(i))), &
!!!             & fVolFlux / (fu_max_d(aerModes(i)) - fu_min_d(aerModes(i)))  /))
!!!    fTmp = fTmp * 1.3
!!!  end do
!!!stop  
!!!  
!!!  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
!!!!!    call set_aerosol_mode(aerModes(1), '', &                         ! mode, chNm, 
!!!!!                                & 1e-6, &           ! fp1
!!!!!                                & 1., &  ! fp2
!!!!!                                & 1.5e-6, 2.6e3, &   ! mass_mean_d, dens, 
!!!!!                                & lognormal_flag, 1) 
!!!!!    
!!!!!    call msg('Absolute norm for lognorm nbr:',fu_integrate_number(1e-9, 1e-3, aerModes(1)))
!!!!!    call msg('Absolute norm for lognorm vol:',fu_integrate_volume(1e-9, 1e-3, aerModes(1)))
!!!!    
!!!!    aerModes(1) = fu_set_mode(fixed_diameter_flag, 1e-10, 1e-2, 1e-6)
!!!!    call wind_blown_dust_flux4mode(aerModes(1), &        ! definition of the spectrum band
!!!!!                                 & 6110, &  !private lognorm_4_modes_dust_flag, &
!!!!!                                 & 6111, &  !private Kok spectrumlognorm_4_modes_dust_flag, &
!!!!                                 & 6112, &  !private lognorm_4_modes_dust__numeric_flag, &
!!!!                                 & fNbrFlux, fVolFlux, fMassMeanDiam)
!!!!    call msg('Mode min,mean,max diam, [um], fluxNbr, fluxVol, meanActualDiam: ', &
!!!!           & (/fu_min_d(aerModes(1))*1e6, fu_mean_d(aerModes(1))*1e6, fu_max_d(aerModes(1))*1e6, &
!!!!             & fNbrFlux, fVolFlux, fMassMeanDiam*1e6/))
!!!!  
!!!!  
!!!!  fTmp = 1e-9
!!!!  do i = 1, 50
!!!!    print *, i
!!!!    aerModes(i) = fu_set_mode(fixed_diameter_flag, fTmp, 1.3*fTmp, 1.2*fTmp)
!!!!    
!!!!    call wind_blown_dust_flux4mode(aerModes(i), &        ! definition of the spectrum band
!!!!!                                 & 6110, &  !private lognorm_4_modes_dust_flag, &
!!!!                                 & 6111, &  !private Kok spectrumlognorm_4_modes_dust_flag, &
!!!!!                                 & 6112, &  !private lognorm_4_modes_dust__numeric_flag, &
!!!!                                 & fNbrFlux, fVolFlux, fMassMeanDiam)
!!!!    call msg('Mode min,mean,max diam, [um], fluxNbr, fluxVol, meanActualDiam, unitfluxNbr, unitfluxVol: ', &
!!!!           & (/fu_min_d(aerModes(i))*1e6, fu_mean_d(aerModes(i))*1e6, fu_max_d(aerModes(i))*1e6, &
!!!!             & fNbrFlux, fVolFlux, fMassMeanDiam, fNbrFlux / (fu_max_d(aerModes(i)) - fu_min_d(aerModes(i))), &
!!!!                                                & fVolFlux / (fu_max_d(aerModes(i)) - fu_min_d(aerModes(i)))/))
!!!!    fTmp = fTmp * 1.3
!!!!  end do
!!!!stop  
  
  
  
  
  
!  call system_mem_usage(run_log_funit)
  
  
  
  
  
  
  call smpi_global_barrier()

  if (smpi_global_rank == 0) then
    PRINT * ;   PRINT * ;   PRINT * ;
  endif

  CALL msg ('Hello world! This is ' // revision_str // ' speaking.')

  
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

  call msg('Control file name:' + control_fname)
  call msg('')

  CALL start_silam_v5(control_fname, had_error)
  
  chTmp = 'Overall_run_time'
  call stop_count(chCounterNm = chTmp)
  call report_time()

  ! Before renaming the file, it must be closed (Intel compiler requirement)
  !
  
  call msg('Local time now: ' + fu_computer_time_string())
  call msg('UTC time now: ' + fu_str(fu_wallclock()))
  
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

 subroutine warning_sigusr1
     !$omp master
     if (smpi_global_tasks > 1) then 
       call msg("got sigusr1, sleeping before setting error", smpi_global_rank)
       call sleep(smpi_global_rank) ! do not mix-up mpi members
       write (6, '(a,x,i4)') "backtrace from rank",smpi_global_rank
       flush (6) !

       call backtrace_md()
     endif
     flush (run_log_funit) !make sure that last message is there
     call set_error("silam got a signal (sigusr1)","warning_sigusr1")
     !$omp end master
 end subroutine warning_sigusr1
    
 subroutine warning_sigint
     !$OMP MASTER
     if (smpi_global_tasks > 1) then 
       call msg("Got SIGINT, sleeping before setting error", smpi_global_rank)
       call sleep(smpi_global_rank) ! Do not mix-up mpi members
       write (6, '(A,X,I4)') "Backtrace from rank",smpi_global_rank
       FLUSH (6) !

       call backtrace_md()
     endif
     FLUSH (run_log_funit) !Make sure that last message is there
     call set_error("Silam got a signal (SIGINT)","warning_sigint")
     !$OMP END MASTER
 end subroutine warning_sigint
    
END PROGRAM silam_main

