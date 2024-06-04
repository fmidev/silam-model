PROGRAM test_modules

  ! Description:
  ! Silam version 4.7 model main program
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
  ! Current code owner Mikhail Sofiev, FMI.
  ! email: Mikhail.Sofiev@fmi.fi
  !
  ! Language: ANSI Fortran 90
  !
  ! All units: SI

  USE dispersion_models

  IMPLICIT NONE

  CHARACTER (LEN=fnlen) :: control_fname
  character(len=16) :: chTmp
  real :: fTmp
  integer :: iStatus, iOut, iRec, iTmp

  real, dimension(:,:), pointer :: fGrd
  real, dimension(:), pointer :: fCentre, fWind, garbage, vCellBorder, fLowMass

  !
  ! Start the global time counter
  !
  chTmp = 'Overall_run_time'
  call start_count(chTmp)
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
  call random_number(fTmp)
  if(fTmp < 0.1) fTmp = fTmp + 0.1
  write(unit=chTmp,fmt='(A8,I4,A4)')'run_tmp_', int(fTmp*10000.), '.log'
  open(run_log_funit, file=chTmp, iostat = iStatus)
  if(iStatus /= 0)then
    call set_error(fu_connect_strings('Failed to open:',chTmp),'silam_main')
    stop
  endif

  PRINT * ;   PRINT * ;   PRINT * ; 
  CALL msg ('Hello world! This is test SILAM speaking')
  call msg('SILAM executable parameters:')
  !
  ! Display the current executable parameters
  !
  call fu_getarg(0, control_fname)    ! temporary use of the variable
  call display_file_parameters(control_fname)
  call msg('')
  !
  ! Time marks
  !
  call msg(fu_connect_strings('Local time now: ',fu_computer_time_string() ))
  call msg(fu_connect_strings('UTC time now: ', fu_time_string_UTC(fu_wallclock())))


  !
  ! Get the control file name
  !
!  SELECT CASE (fu_iargc())
!    CASE(0) !--------- No arguments, read all from the default ini file
!      call msg('No arguments - use default silam.ini')
!      IF(.not.read_ini_file('silam.ini'))THEN
!        CALL set_error('Can not read silam.ini','silam_main')
!        STOP
!      END IF
!
!    CASE(1) !-------- One argument - ini file name, read from this ini file
!      call msg('1 argument - use it as control file name')
!      CALL fu_getarg(1, control_fname)
!
!    CASE DEFAULT
!      CALL set_error('Too many input arguments','silam_main')
!  END SELECT
!
!  call msg(fu_connect_strings('Control file name:',control_fname))
!  call msg('')
!
!  CALL start_silam_v5(control_fname)


  !
  ! Advection tests
  !
!  fGrd => fu_work_array_2d()
!  fCentre => fu_work_array()
!  fWind => fu_work_array()
!  vCellBorder => fu_work_array()
!  garbage => fu_work_array()
!  fLowMass => fu_work_array()

  iOut = fu_next_free_unit()
  if(error)stop

  allocate(fGrd(1,0:201), fCentre(0:201), vCellBorder(0:201),fLowMass(1),garbage(1), fWind(200))


  fGrd = 0.0
  fLowMass = 0.001
  garbage = 0.0
  do iTmp = 0, 201
    vCellBorder(iTmp) = real(iTmp)-0.5
    fCentre(iTmp) = iTmp
  end do


  open(iOut,file='output\adv.tst_fixwind_newadv',form='unformatted',recl=200,access='direct')
  iRec = 1
  !
  ! Arbitrary profile
  !
  do iTmp = 10,75
    fGrd(1,iTmp) = sqrt(real(iTmp))
  end do
  !
  ! Constant wind
  !
  fWind = 0.5  ! in relative units
  !
  ! Growing wind
  !
  do iTmp = 0, 201
!    fWind(iTmp) = 0.004*real(iTmp)
!    fWind(iTmp) = 0.8-0.004*real(iTmp)
  end do


  do iTmp = 1, 500
    call advect_Galperin_bulk_1d_abs(200, 1, &
                                   & fGrd, fCentre, &
                                   & fWind, vCellBorder, 1.0, &
                                   & fLowMass, garbage)
!    !
!    ! A semi-Lagrangian scheme of M.Galperin, as in (Galperin, 2000)
!    ! presented at NATO/CCMS and published in IGCE series. 
!    !
!    ! Imported parameters
!    integer, intent(in) :: nCells, nSpecies
!    real, dimension(1:nSpecies, 0:nCells+1), intent(inout) :: vMass
!    real, dimension(0:nCells+1), intent(inout) :: vCentre, vCellBorder
!    real, dimension(:), pointer :: vWind_abs, & ! must be same linear units!
!                                 & garbage, vMinAdvectedMass
    if(mod(iTmp,10) == 0)print *,iTmp

    write(iOut,rec=iRec)fGrd(1:1,1:200)
    iRec = iRec+1
  end do




  chTmp = 'Overall_run_time'
  call stop_count(chCounterNm = chTmp)
  call report_time()


  STOP


END PROGRAM test_modules

