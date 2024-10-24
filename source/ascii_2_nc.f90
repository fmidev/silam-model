!********************************************************************************
!
! This programme converts silamASCII fileld to Silam flavour of NetCDF
!
!********************************************************************************

program ascii_2_nc
  !
  ! The main program will only find out the name of the ini file and call the
  ! actual procedure
  !
  use grads_io
  use chemical_setup
  use grib_api
  use revision
  use silam_times
  use netcdf_io
  use supermarket_of_fields

  implicit none

  ! Local variables
  character(len=fnlen) :: chInFNm, chOutFNm, chChemFNm
  integer :: iStatus, unit, j, iTmp
  type(Tsilam_namelist_group), pointer :: nlGrpIni
  character(len=*), parameter :: sub_name = 'ascii_2_nc'

  !----------------------------------------
  !



  !! No need for a separate log file. hush it here
  run_log_funit = 6 !Default
  run_log_name = "/dev/null"
  run_log_funit = fu_next_free_unit()
  open(run_log_funit, file=run_log_name, iostat = iStatus)
  if(iStatus /= 0)then  ! The
    call set_error('Failed to open: "'//trim(run_log_name)//'"', sub_name)
    stop
  endif


  CALL msg ('Hello world! This is ascii_2_nc converter speaking.')
  call msg(fu_connect_strings('Local time now: ', fu_computer_time_string()))
  call msg(fu_connect_strings('UTC time now: ', fu_str(fu_wallclock())))

  iTmp = command_argument_count()
  SELECT CASE (iTmp)
    CASE(2) !- Try to guess the location  of the chemical database
      CALL get_command_argument(0, chInFNm) !! program name
      iTmp = index(chInFNm,'/', back=.TRUE.)
      chChemFNm = chInFNm(1:iTmp)//"../ini/silam_chemicals_95_OC.dat"
      CALL get_command_argument(1, chInFNm)
      CALL get_command_argument(2, chOutFNm)

    CASE(3) !-
      CALL get_command_argument(1, chChemFNm)
      CALL get_command_argument(2, chInFNm)
      CALL get_command_argument(3, chOutFNm)

    CASE DEFAULT
      do j=1,iTmp
        call get_command_argument(j, chInFNm)
        call msg("argument "//trim(fu_str(j))//": "//trim(chInFNm) )
      enddo

      CALL get_command_argument(0, chInFNm) !! program name
      CALL set_error('Usage: '//trim(chInFNm)//" silam_chemicals.dat  input.asc output.nc", sub_name)
  END SELECT

  if(.not. error) call init_chemical_materials(chChemFNm, "")

  if (.not. error) then
    iTmp = len_trim(chInFNm)
    if (index(chInFNm, ".super_ctl") + 9 == iTmp ) then
      CALL grads_2_nc(chInFNm, chOutFNm)
    else
       CALL convert_2_nc(chInFNm, chOutFNm)
    endif
  endif
  if (error) then
    print *, "Failed!"
    CALL exit_with_status(9)
  else
    print *, "Done!"
    CALL exit_with_status(0)
  end if


  CONTAINS


  !**************************************************************************

  subroutine convert_2_nc(chInFNm, chOutFNm)
    !
    !
    implicit none

    ! Imported parameter
    character(len=*), intent(in) :: chInFNm, chOutFNm


    !!!!
    integer :: iUnit, iStat
    logical :: eof
    real, dimension(:), allocatable :: grid_data
    type(silja_field_id) :: field_id
    type(silam_vertical) :: vertical
    type(TOutputList), dimension(3) :: OutLst
    character(len=*), parameter :: sub_name = 'convert_2_nc'



    call open_ascii_file_i(chInFNm, iUnit)

    allocate(grid_data(worksize), OutLst(1)%ptrItem(1), & 
       & OutLst(2)%ptrItem(0), OutLst(3)%ptrItem(0), stat=iStat)
    if (iStat /= 0) then
      call set_error('allocate failed!', sub_name)
      return
    end if

    !! Use default filling from super_ctl
    call read_next_field_from_ascii_file(iUnit, eof, field_id, grid_data)

    call close_ascii_file_i(iUnit)


    OutLst(1)%ptrItem(1)%quantity = fu_quantity(field_id)
    OutLst(1)%ptrItem(1)%AvType = fu_field_kind(field_id)
    OutLst(1)%ptrItem(1)%AvPeriod = fu_accumulation_length(field_id)
    OutLst(1)%ptrItem(1)%species = fu_species(field_id)  ! missing if no substance defined
    OutLst(1)%ptrItem(1)%chSpecies_string = fu_str(fu_species(field_id))  ! just input string
    OutLst(1)%ptrItem(1)%iSpeciesListType = int_missing  ! can be iNoSubstanceRelation
    OutLst(1)%ptrItem(1)%iVerticalTreatment = int_missing
    OutLst(1)%ptrItem(1)%if3D = .false.

    call set_vertical(fu_level(field_id), vertical)

    iUnit = open_netcdf_file_o(chOutFNm, fu_grid(field_id), vertical, fu_analysis_time(field_id), &
                                    & OutLst, &
                                    & "", .True., 4, .false., real_missing)

    call write_next_field_to_netcdf_file(iUnit, field_id, grid_data)

    call  close_netcdf_file(iUnit)

  end subroutine convert_2_nc


  !**************************************************************************

  subroutine grads_2_nc(chInFNm, chOutFNm)
    !
    !
    implicit none

    ! Imported parameter
    character(len=*), intent(in) :: chInFNm, chOutFNm


    !!!!
    integer :: igFile, iUnit, iStat 
    integer :: nvars, ntimes, nlevels
    integer :: iVar, iTime, iLev

    logical :: eof
    real, dimension(:), allocatable :: grid_data
    type(silja_field_id) :: field_id
    type(silam_vertical), pointer :: vertical
    type(TOutputList), dimension(3) :: OutLst
    character(len=*), parameter :: sub_name = 'grads_2_nc'

    call init_grads_io(1) !! Only one file needed
    call msg("Opening super_ctl file: "//trim(chInFNm))
    igFile = fu_open_gradsfile_i(chInFNm)

    if (error) return

    nvars = fu_n_gvars(igFile)


    allocate(grid_data(worksize), OutLst(1)%ptrItem(nvars), & 
       & OutLst(2)%ptrItem(0), OutLst(3)%ptrItem(0), stat=iStat)
    if (iStat /= 0) then
      call set_error('allocate failed!', sub_name)
      return
    end if


    do iVar = 1,nvars  
        call  get_grads_var_metadata(igFile, iVar, 1, 1, field_id) ! indices: gfile, gvar, glev, gtime; SILAM-id
        OutLst(1)%ptrItem(iVar)%quantity = fu_quantity(field_id)
        OutLst(1)%ptrItem(iVar)%AvType = fu_field_kind(field_id)
        OutLst(1)%ptrItem(iVar)%AvPeriod = fu_accumulation_length(field_id)
        OutLst(1)%ptrItem(iVar)%species = fu_species(field_id)  ! missing if no substance defined
        OutLst(1)%ptrItem(iVar)%chSpecies_string = fu_str(fu_species(field_id))  ! just input string
        OutLst(1)%ptrItem(iVar)%iSpeciesListType = int_missing  ! can be iNoSubstanceRelation
        OutLst(1)%ptrItem(iVar)%iVerticalTreatment = int_missing
        OutLst(1)%ptrItem(iVar)%if3D = (fu_n_gVar_levs(igFile,iVar) > 1 ) 
    enddo

    nTimes = fu_n_gtimes(igFile)

    iUnit = open_netcdf_file_o(chOutFNm, fu_silamGrid_of_grads(igFile), &
                                      &  fu_silamVert_of_grads(igFile), &
                                      &  fu_time_of_grads(igFile, 1), & !! First time -- analysis
                                      & OutLst, &
                                      & "", .True., 4, .false., real_missing)

    do iTime = 1, nTimes                                 
      do iVar = 1,nvars 
          do iLev = 1,fu_n_gVar_levs(igFile,iVar) 
            call  get_grads_var_metadata(igFile, iVar, iLev, iTime, field_id) ! indices: gfile, gvar, glev, gtime; SILAM-id
            call  read_field_from_grads_indices(igFile, iVar, iLev, iTime, grid_data, real_missing)

            call write_next_field_to_netcdf_file(iUnit, field_id, grid_data)
          enddo
      enddo
    enddo

    call  close_netcdf_file(iUnit)

  end subroutine grads_2_nc


end program ascii_2_nc

