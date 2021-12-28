!********************************************************************************
!
! This programme converts GrADS file written with super_ctl 
! to GRIB 2 format using GRIB-API interface
!
! Code owner: Mikhail Sofiev, FMI
!
! Updated for SILAM5.1, GRIB-API 1.9.9 - JV 4/2012.
!
! All units: NOT NESESSARILY SI. However, they stay SI unless opposite is stated
!
!********************************************************************************

program silam_2_grib2
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
  use silam_mpi

  implicit none
  
  ! Local variables
  character(len=fnlen) :: chIniFNm, defultInputFnm, rankTmp,chTmp
  real :: fTmp
  integer :: iStatus, unit, j, iTmp
  type(Tsilam_namelist_group), pointer :: nlGrpIni

  !----------------------------------------
  !
  !  Open the global error and warning recording file
  !  This very file will be used for error messages, warnings and, later,
  !  for other log functions
  !  However, when the output directory becomes known - it will be transferred to it
  !  and its number will be changed
  !


  run_log_funit = 6 !Default
  
 !   call grib_tests(1)
 !   stop

  call  smpi_init()
  

  if(smpi_is_mpi_version())then
    ! Just use what we have to ensure that the string is the same
    call msg("Synchronizing tmp file prefix")
    call smpi_allreduce_max_int(fu_pid(), j, MPI_COMM_WORLD)
    write(unit=run_log_name,fmt='(A,I5.5,A,I2.2,A)') "run_",j,"_",smpi_global_rank,".log"
  else
    j = fu_pid()
    write(unit=run_log_name,fmt='(A,I5.5,A)') "run_",j,".log"
  endif
  run_log_name = "/dev/null"

  run_log_funit = fu_next_free_unit()
  open(run_log_funit, file=run_log_name, iostat = iStatus)
  if(iStatus /= 0)then  ! The problem is serious, randomised file name does not help
    call set_error('Failed to open: "'//trim(run_log_name)//'"','grads_2_grib2_main')
  endif

  if (smpi_global_rank == 0) then
    PRINT * ;   PRINT * ;   PRINT * ;
  endif

  if(smpi_is_mpi_version())then
     write(unit=rankTmp,fmt='(A25,I3,A6,I3)') 'MPI version running with ', smpi_global_tasks, ' tasks, task ', smpi_global_rank
  else
     write(unit=rankTmp,fmt='(A14)') 'Serial version'
  endif
  call msg(rankTmp)

  CALL msg ('Hello world! This is silam_2_GRIB2 converter speaking.')
  call msg(fu_connect_strings('Local time now: ', fu_computer_time_string()))
  call msg(fu_connect_strings('UTC time now: ', fu_str(fu_wallclock())))

  defultInputFnm=''
  iTmp = command_argument_count()
  SELECT CASE (iTmp)
    CASE(0) !--------- No arguments, read all from the default ini file
      call msg('No arguments - use default grads_2_grib2.ini')
      chIniFNm = 'grads_2_grib2.ini'

    CASE(1) !-------- One argument - ini file name, read from this ini file
      call msg('1 argument - use it as ini file name')
      CALL get_command_argument(1, chIniFNm) 

    CASE(2) !-
      call msg('2 arguments - use it as ini file name and default input name')
      CALL get_command_argument(1, chIniFNm) 
      CALL get_command_argument(2, defultInputFnm) 

    CASE DEFAULT
      do j=1,iTmp
        call get_command_argument(j, chIniFNm)
        call msg("argument "//trim(fu_str(j))//": "//trim(chIniFNm) ) 
      enddo
      
      CALL set_error('Too many input arguments','grads_2_grib2_main')

  END SELECT

  if (.not. error) then
      call msg(fu_connect_strings('Ini file name:',chIniFNm))
      call msg('')

      !
      ! Initialize what is needed
      !
      call init_grads_io(10)  ! number of GrADS files
  endif

  if(.not. error) then
    unit = fu_next_free_unit()
    open(unit, file=chIniFNm, action='read', status='old', iostat=iStatus)
    IF(iStatus /= 0) THEN
      CALL set_error(fu_connect_strings('Ini file does not exist:',chIniFNm),'convert_grads_2_grib2')
    END IF
  endif


  if(.not. error) nlGrpIni => fu_read_namelist_group(unit,.true.)        ! get the ini namelist group
  close(unit)
  


  if(.not. error) CALL convert_2_grib2(nlGrpIni,defultInputFnm)

  if (smpi_is_mpi_version() .and. smpi_global_tasks > 1) then
    print *, "Exiting MPI"
    if (error) then
      call smpi_abort(9)
    else
      call smpi_finalize()
    end if
  else if (error) then
    print *, "Failed!"
    call exit_with_status(9)
  else
    print *, "Done!"
    call exit_with_status(0)
  end if


  CONTAINS


  !**************************************************************************

  subroutine convert_2_grib2(nlGrpIni,defultInputFnm)
    !
    ! This subroutine performs actual conversion of GrADS file to GRIB-2
    ! Steps: read the GrADS super-ctl file
    !
    implicit none

    ! Imported parameter
    type(Tsilam_namelist_group), pointer :: nlGrpIni
    character(len=*), intent(in) :: defultInputFnm

    ! Local variables
    INTEGER :: status, iGF, nNamelists, iNl, uGribTempl, uGribOut, iLev, nInLevs, indTime
    integer :: iGribMessageIndex, iGribMessageIndexDefTmpl, iCheckGribMessage, iStartIndex, &
         & iEndIndex, ix, iy, ind_lev_in, num_grads_times
    type(Tsilam_namelist), pointer :: nlIni, nl_common, nl_grib2_defaults, nl_species_map
    logical :: lExist
    character(len=fnlen) :: chInputFNm, chTmp, chMsg, OpenedFileName, chOutTemplate
    type(silam_fformat) ::  fform ! GRIB file, ASCII file, TEST_FIELD, .....
    type(silja_grid) :: gGrid
    type(silja_grid), dimension(:), pointer :: gridPtr
    type(silam_vertical), pointer :: gVert
    type(silam_vertical), dimension(:), pointer :: vertPtr
    type(silja_level) :: output_level, input_level
    type(silja_field_id) :: field_id
    real, dimension(:), pointer :: grads_data
    real :: fScaling
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: itemInLevels
    type(silja_time) :: an_time, an_time_force, now, timeStart, timeEnd
    real :: corner_lon_E, corner_lat_N, southpole_lon_E, southpole_lat_N, dx_deg, dy_deg
    logical :: corner_in_geographical_latlon, if_south_pole
    integer :: number_of_points_x, number_of_points_y, iIn, iOut, ixStart, ixEnd, iyStart, &
            & iyEnd, fsIn, fsOut, nGrids, nVerts, nLevels
    type(grads_template) :: outputGribTemplate
    real(8), allocatable :: output_data(:)
    real(8) :: fSum, f8
    type(Taerosol_mode) :: aerosol_mode
    real :: mode_value, fraction
    !integer :: ixStart_def, ixEnd_def, iyStart_def, iyEnd_def, num_times_def
    !real :: scaling_def
    !character(len=clen) :: quantity_name_def
    character(len=fnlen) :: template_def, template_filename
    character(len=*), parameter :: sub_name = 'convert_2_grib2'
    character(len=clen) :: date_mode
    real ::  findex
    real, dimension(2) :: fract_levs_in
    integer, dimension(2) :: ind_levs_in
    type(silja_time), dimension(:), pointer :: grads_times
    type(silam_species) :: species
    integer :: constituent_code, parameter_code

    vertPtr => null()
    gridPtr => null()
    OpenedFileName = ""

    allocate(output_data(worksize), stat=status)
    if (status /= 0) then
      call set_error('allocate failed!', sub_name)
      return
    end if

    nl_common => fu_namelist(nlGrpIni, 'common')
    if (.not. associated(nl_common)) then
       call set_error('No COMMON namelist', sub_name)
       return
    end if
    call init_chemical_materials(fu_expand_environment(fu_content(nl_common,'chemical_database_fnm')), "")
    call init_netcdf_io(nl_common,nl_common)
    if (error) return

    ! Make template index to iGribMessageIndexDefTmpl
    template_def = fu_content(nl_common, 'grib2_template_file_default')
    iGribMessageIndexDefTmpl = int_missing
    if (template_def /= '') then
       template_filename = fu_process_filepath(fu_content(nl_common,'grib2_template_file_default'))
       if (fu_fails(template_filename /= '', 'Missing  grib2_template_file_default', sub_name)) return
       call grib_open_file(uGribTempl, template_filename, 'r')
       call grib_new_from_file(uGribTempl, iGribMessageIndexDefTmpl)
       call grib_close_file(uGribTempl)
    end if




    ! grid subsetting was overcomplicated and cdo can do regridding much better (both for grib and for .nc)
    if (fu_content_real(nl_common, 'lon_start') /= real_missing) then
      call msg_warning("Looks like grid spec is present the namelist..", sub_name)
      call msg_warning("Grid subsetting no longer supported. Pease remove the request for it..", sub_name)
      return
    endif

    ! just recycle variable
    date_mode = fu_content(nl_common, 'force_an_time')
    if (len_trim(date_mode) > 0) then
      an_time_force = fu_io_string_to_time(date_mode)
      if(error) then
        call set_error("failed to parse force_an_time", sub_name)
        return
      endif
      call msg("Forcing analysis time: "//trim(date_mode))
    else
      an_time_force = time_missing
    endif


    date_mode = fu_content(nl_common, 'date_mode')
    if (fu_fails(date_mode /= '', 'Missing date_mode', sub_name)) return

    nl_grib2_defaults => fu_namelist(nlGrpIni, 'grib2_defaults')
    if (associated(nl_grib2_defaults)) then
      call msg('Global GRIB2 keys and values:')
      call report(nl_grib2_defaults)
    end if

    nl_species_map => fu_namelist(nlGrpIni, 'species_map')

    chOutTemplate = fu_content(nl_common,'grib2_output_file')
    if (smpi_global_tasks > 1) then  ! Either forecast length or all others together
        if(index(chOutTemplate,trim(fu_forecast_length_templ_str()))==0)then
          if((index(chOutTemplate,trim(fu_valid_time_hour_templ_str()))==0).or. &
           & (index(chOutTemplate,trim(fu_valid_time_day_templ_str()))==0).or. &
           & (index(chOutTemplate,trim(fu_valid_time_month_templ_str()))==0).or. &
           & (index(chOutTemplate,trim(fu_valid_time_year_templ_str()))==0))then
            call set_error('Not enough varying items in template for hourly series', &
                         & sub_name)
            return
          endif
        endif
    endif
    call decode_template_string(chOutTemplate, outputGribTemplate)


    nNamelists = fu_nbr_of_namelists(nlGrpIni)      ! get the number of namelists
    call msg('Sets to be processed:',nNamelists)
    if(error .or. nNamelists < 1)return

    grads_data => fu_work_array()
    if(error)return

    !------------------------------------------------------------------------
    !
    ! Start the grand cycle over the namelists
    !
    do iNl = 1, nNamelists
      
      nlIni => fu_namelist(nlGrpIni,iNl)
      if (fu_name(nlIni) == 'COMMON') cycle
      if (fu_name(nlIni) == 'GRIB2_DEFAULTS') cycle
      if (fu_name(nlIni) == 'SPECIES_MAP') cycle
      call msg('Namelist name:' // trim(fu_name(nlIni)))
      if(error .or. .not. associated(nlIni))return
      !
      ! Get the scaling
      !
      if(error)return
      
      !
      ! Open the GrADS data file. Note that the GrADS super-ctl file name is sitting
      ! inside the list file, which contains the last created file. Get it!
      !

      chInputFNm = fu_content(nlIni,"grads_super_ctl_file")
      if (chInputFNm /= '') then
         call msg("Please use input_file = GRADS "//trim(chInputFNm))
         call msg_warning("Old grads file requested", sub_name)
         FForm%iformat = grads_file_flag
      else
         chInputFNm = fu_content(nlIni,'input_file')
         if (fu_fails(chInputFNm /= '', 'Missing input_file', sub_name))return
         fform = fu_input_file_format(chInputFNm)
         if(error)return
         
         chInputFNm=adjustl(chInputFNm)
         chInputFNm=adjustl(chInputFNm(index(chInputFNm,' ')+1:))
         if (chInputFNm == '-') chInputFNm = defultInputFnm
      endif

!      iGF = fu_open_gradsfile_i(chInputFNm)
      if (trim(OpenedFileName) /= trim(chInputFNm) ) then

          !Close if not needed
          if (trim(OpenedFileName) /= '') then
                call msg("Closing input file: "//trim(OpenedFileName))
                if (FForm%iformat == netcdf_file_flag) then
                    call close_netcdf_file(iGF)
               else
                    call close_gradsfile_i(iGF)
               endif
          endif
          call open_input_file(chInputFNm, FForm, iGf, int_missing)
          OpenedFileName = chInputFNm
      
          if(error .or. iGF < 0)return

          ! Get the grid, vertical and time list
          !
          select case(FForm%iformat)
            case(netcdf_file_flag)
               !          call get_content_from_netcdf_file(inputName, inputFormat, q_list, InputContent)
               call get_netcdf_grids(iGf, int_missing, species_missing, gridPtr, nGrids)
               if ( nGrids < 1) then
                 call msg("Wrong number of grids "//trim(fu_str(nGrids))//" in file: "//trim(OpenedFileName))
                 call set_error("At least one grid should be", sub_name)
                 return
               endif
               gGrid = gridPtr(1) !Take first grid

               call get_netcdf_verticals(iGf, int_missing, species_missing, vertPtr, nVerts)
               do ix=1,nVerts
                  gVert => vertPtr(ix)
                  !call report(gVert, .true.)
                  if (defined(gVert)) then  
                    !!! Dirty hack. We try to find first vertical with more than one level if possible
                    ! first valid vertical might be single-level, then no way to extract anything above surface
                    if (fu_NbrOfLevels(gVert) > 1) exit
                  endif
               enddo
               !if ( ix>nVerts ) then
               !   call set_error("All "//trim(fu_str(nVerts))//" verticals in NetCDF file undefined", sub_name )
               !   return
               !endif 

               call timelst_from_netcdf_file(iGf, grads_times, num_grads_times, an_time = an_time)
               if(error) return

            case(grads_file_flag)
               gGrid = fu_silamGrid_of_grads(iGF)
               if(error)return

               gVert => fu_silamVert_of_grads(iGF)
               if(error)return
               
               call get_grads_times(iGF, grads_times, num_grads_times, ifEnvelope=.false.)
               an_time = grads_times(1)
               if (error) return

            case default
              call msg ("only grads and netcdf supported")
              call msg('Unknown input format:',FForm%iformat)
              call set_error('Unknown input format',sub_name)
              return
          end select
          if (defined(an_time_force)) an_time = an_time_force
          call report(gVert,.true.)
          call report(gGrid)

          call lonlat_grid_parameters(gGrid, corner_lon_E, corner_lat_N, corner_in_geographical_latlon, &
                                         & number_of_points_x, number_of_points_y, &
                                         & southpole_lon_E, southpole_lat_N, &
                                         & dx_deg, dy_deg)
          if(error)return
          if (abs(southpole_lat_N +90) > 0.1) then
            call msg_warning("Pole not in geographical pole",sub_name)
            call report(gGrid)
            call set_error("Rotated grids not supported for GRIB output (yet?)", sub_name)
            return
          endif


          fsIn = number_of_points_x * number_of_points_y
          nLevels = fu_NbrOfLevels(gVert)
          
          ! ifHorizInterp = (gGrid == outGrid)


      endif


      ! Open GRIB-2 template file
      !
      if (fu_content(nlIni,'grib2_template_file') == '') then
         !Try to use the default one
         if (iGribMessageIndexDefTmpl > 0) then 
            call grib_clone(iGribMessageIndexDefTmpl, iGribMessageIndex)
         else
            call msg_warning('No grib2_template_file and no default_grib2_template_file', sub_name)
            call set_error('Wrong iGribMessageIndexDefTmpl='+fu_str(iGribMessageIndexDefTmpl), sub_name)
         endif
      else
         template_filename = fu_process_filepath(fu_content(nlIni,'grib2_template_file'))
         if (fu_fails(template_filename /= '', 'Missing grib2_template_file', sub_name)) return
         call grib_open_file(uGribTempl, template_filename, 'r')
         call grib_new_from_file(uGribTempl, iGribMessageIndex)
      end if
      
!      call msg('message extracted',iGribMessageIndex)

      call msg('Variable:' // fu_content(nlIni,'quantity_name')+'_'+fu_content(nlIni,'substance_name'))


      ! Levels. Either output_level and automatic interpolation or
      ! input_level_and_fraction.
      nullify(itemInLevels)
      call get_items(nlIni,'output_level_and_scaling', itemInLevels, nInLevs)
      if(error)return
      if (nInLevs /= 1) then
        call msg("Several output_level_and_scaling given", nInLevs)
        call set_error('Must have only one output_level_and_scaling', sub_name)
        return
      endif

      call set_named_level_with_fract(itemInLevels(1), output_level, fScaling)
      fIndex = fu_project_level_crude(output_level, gVert, clip=.false.)
      call msg('Output level:')
      call report(output_level)
      if (fu_fails(fIndex > 0, 'Failed to project level', sub_name)) return


      !Interp coeffs
      nInLevs = 1 !One level by defalut
      fract_levs_in(1) = 1.
      if (fIndex <= 1.) then !Could be below mid-point of upper layer
        ind_levs_in(1) = nint(fIndex + 1e-6) !should be 1
      elseif(fIndex > nLevels) then !Could be above mid-point of upper layer
        ind_levs_in(1) = nint(fIndex - 1e-5)
      else  !Normal interpolation (Should nearest point be here?)
        nInLevs = 2
        ind_levs_in(1) = floor(fIndex + 1e-6)
        ind_levs_in(2) = ind_levs_in(1) + 1
        fract_levs_in(2) = mod(fIndex, 1.)
        fract_levs_in(1) = 1. - fract_levs_in(2)
      endif
      call msg('Lower interpolated level '//trim(fu_str(ind_levs_in(1)))//' frac:', fract_levs_in(1))
      call report(fu_level(gVert, ind_levs_in(1)))
      if (nInLevs > 1) then
              call msg('Upper interpolated level '//trim(fu_str(ind_levs_in(2)))//' frac:', fract_levs_in(2))
              call report(fu_level(gVert, ind_levs_in(2)))
      endif
        
      !
      ! Determine start time, start and end indices in the grds file
      !
      if(len_trim(fu_content(nlIni,'start_time')) > 0)then
        !
        ! If start_time is defined, get it together with end_time and compute the grads time indices
        !
        timeStart = fu_io_string_to_time(fu_content(nlIni,'start_time'))
        if(error .or. .not.defined(timeStart))then
          call set_error('start_time item is wrong',sub_name)
         return
        endif
        timeEnd = fu_io_string_to_time(fu_content(nlIni,'end_time'))
        if(error .or. .not.defined(timeEnd))then
          call set_error('end_time item is wrong',sub_name)
          return
        endif
        iStartIndex = fu_index( timeStart, grads_times)
        iEndIndex = fu_index(timeEnd, grads_times)

      elseif(len_trim(fu_content(nlIni,'number_of_latest_times')) > 0)then
        !
        ! Possibly, we just need N last moments from the file
        !
        iEndIndex = fu_n_gtimes(num_grads_times)
        iStartIndex = iEndIndex - fu_content_int(nlIni,'number_of_latest_times') + 1
        timeStart = grads_times(iStartIndex)
      else if (fu_content_int(nlIni, 'start_time_index') /= int_missing &
             & .and. fu_content_int(nlIni, 'end_time_index') /= int_missing) then
        !
        ! indices are defined explicitly
        !
        iStartIndex = fu_content_int(nlIni,'start_time_index')
        iEndIndex = fu_content_int(nlIni,'end_time_index')
        if(iStartIndex == int_missing .or. iEndIndex == int_missing)then
          call set_error('iStartIndex == int_missing .or. iEndIndex == int_missing', sub_name)
          return
        endif
        timeStart = grads_times(iStartIndex)
      else
        ! Take all times
        !
        iStartIndex = 1
        iEndIndex = num_grads_times
        timeStart = grads_times(1)
      endif

      ! And check the stuff
      !
      if(iStartIndex < 1 .or. iStartIndex > iEndIndex .or. iEndIndex > num_grads_times)then
        call set_error('Failed to determine the time indices',sub_name)
        call report(nlIni)
        return
      endif
      call msg(fu_connect_strings('TimeStart:',fu_time_to_io_string(timeStart)))

      ! 
      ! Find the grid limits. Global lat/lon limits are in the common namelist, but can be
      ! overriden by the local start indices.
      ixStart = 1
      ixEnd = number_of_points_x
      iyStart = 1
      iyEnd = number_of_points_y
      
      fsOut = fsIn


      if (fu_fails(fsOut <= size(output_data), 'fsOut too small', sub_name)) return
      !
      ! Cycle over time for a single GRIB-2 file
      !
      do indTime = iStartIndex, iEndIndex
        if (mod(indTime,smpi_global_tasks) /= smpi_global_rank) cycle 
        now = grads_times(indTime)

        call msg(fu_connect_strings('Now:',fu_time_to_io_string(now),','), indTime)


        !
        ! Output GRIB file has own time-related template. Get it and open the file
        !
        if (len(trim(fu_content(nlIni,'grib2_output_file'))) > 0) then
          call msg ("grib2_output_file = "//trim(fu_content(nlIni,'grib2_output_file')))
          call set_error("No grib2_output_file is allowed in pitput namelists", &
                                & sub_name)
        endif


        chTmp = fu_FNm(outputGribTemplate, an_time, &  ! ini_time
                     & an_time, &                      ! anal_time
                     & now - an_time) ! forec_len
        call grib_open_file(uGribOut, trim(chTmp), 'a')
        if(error)return

        !
        ! Chemical features may or may not be given
        !
        if(len_trim(fu_content(nlIni, 'substance_name')) > 0)then
!          aerosol_mode = fu_set_mode(nlIni)
!          if (error) return
          call set_species(species, fu_get_material_ptr(fu_content(nlIni, 'substance_name')), &
                                  & fu_set_mode(nlIni))
        else
          species = species_missing
        endif

        !
        ! Since vertical interpolation can take place, have to cycle over vertical levels
        !
        output_data(1:size(output_data)) = 0.

        do iLev = 1, nInLevs
          ind_lev_in = ind_levs_in(ilev)
          fraction = fract_levs_in(ilev)
          input_level = fu_level(gVert, ind_lev_in)
          if (error) return
          field_id = fu_set_field_id(met_src_missing,&
                                   & fu_get_silam_quantity(fu_content(nlIni,'quantity_name')), &
                                   & timeStart, &               ! analysis time
                                   & now - timeStart, &         ! forecast length
                                   & gGrid,&
                                   & input_level,&
                                   & species=species)
          if(error)return

!          call report(field_id)

         select case(FForm%iformat)
           case(netcdf_file_flag)
             call read_field_from_netcdf_file(iGF, field_id, grads_data, real_missing)
           case(grads_file_flag)
             call read_field_from_grads_id(iGF, field_id, grads_data, fill_value_=real_missing)

           case default
              call set_error("Should not be here! Write file select", sub_name)
         endselect
         if(error)return
  
         !call msg('Input sum:',sum(grads_data(1:fsIn)))

          iOut = 1
          do iy = iyStart, iyEnd
            do ix = ixStart, ixEnd
              iIn = ix + (iy - 1) * number_of_points_x
              output_data(iOut) = output_data(iOut) + grads_data(iIn) * fraction
              iOut = iOut + 1
            end do
          end do  ! iy
          !call msg('Output sum after copy:',sum(output_data(1:fsOut)))
!          output_data(1:fs) = output_data(1:fs) + grads_data(:) * fraction

        end do   ! vertical grads levels

        !
        ! Update the scale
        if (.not. (fScaling .eps. 1.0)) output_data(1:fsOut) = output_data(1:fsout) * fScaling
        !call msg('Output sum after scaling:',sum(output_data(1:fsOut)))
        !
        ! If the scanning goes from the north, we have to flip the field
        ! use grads_data as a temporary array
        !
        grads_data(1:fsOut) = output_data(1:fsout)
        do iy = 1, iyEnd - iyStart + 1  !number_of_points_y
          do ix = 1, ixEnd - ixStart + 1  !number_of_points_x
             iIn = ix+(iy-1)*(ixEnd - ixStart + 1)
             iOut = ix+((iyEnd - iyStart + 1) - iy)*(ixEnd - ixStart + 1)
            output_data(iIn) = grads_data(iOut)
          end do
        end do
        
        !-----------------------------------------------------------------------------
        !
        ! Output field is generated. Now update the GRIB message and write it down
        !

        call update_grib2_defaults(iGribMessageIndex, nl_grib2_defaults)
        if (error) return

        ! Time
        !

        if (date_mode == 'forecast') then
          call grib_set_int(iGribMessageIndex,'dataDate', &
              & fu_year(an_time)*10000 + fu_mon(an_time)*100 +fu_day(an_time))

          call grib_set_int(iGribMessageIndex,'dataTime', &  !HHMMSS
              & fu_hour(an_time) * 10000 + fu_min(an_time)*100 + int(fu_sec(an_time)))

          call grib_set_int(iGribMessageIndex, 'startStep',nint(fu_hours(now-an_time)))
          call grib_set_int(iGribMessageIndex, 'significanceOfReferenceTime',1)
        else if (date_mode == 'analysis') then
          call grib_set_int(iGribMessageIndex,'dataDate', &
              & fu_year(now)*10000 + fu_mon(now)*100 +fu_day(now))

          call grib_set_int(iGribMessageIndex, 'startStep',0)
          call grib_set_int(iGribMessageIndex,'hour',  fu_hour(now) )
          call grib_set_int(iGribMessageIndex, 'significanceOfReferenceTime',0)
        else 
          call set_error('Bad date_mode', sub_name)
          return
        end if

        !
        ! Grid parameters
        !
        call grib_set_int(iGribMessageIndex,'numberOfDataPoints', fsOut)
        call grib_set_int(iGribMessageIndex,'numberOfPointsAlongAParallel',(ixEnd - ixStart + 1))
        call grib_set_int(iGribMessageIndex,'numberOfPointsAlongAMeridian',(iyEnd - iyStart + 1))

        f8 = modulo(corner_lon_E + (ixStart - 1) * dx_deg + 360D0, 360D0) 
        call grib_set_real8(iGribMessageIndex,'longitudeOfFirstGridPointInDegrees',f8)

        f8 = modulo(corner_lon_E + (ixEnd - 1) * dx_deg + 360D0, 360D0)
        call grib_set_real8(iGribMessageIndex,'longitudeOfLastGridPointInDegrees',f8)

        !
        ! For first-corner-in-the-south:
        !
!        f8 = corner_lat_N
!        call grib_check(grib_set_real8(iGribMessageIndex,'latitudeOfFirstGridPointInDegrees',f8))
!        f8 = corner_lat_N + (number_of_points_y-1) * dy_deg
!        call grib_check(grib_set_real8(iGribMessageIndex,'latitudeOfLastGridPointInDegrees', f8))

        !
        ! For first-corner-in-the-north:
        !
        f8 = corner_lat_N + (iyEnd - 1) * dy_deg
        call grib_set_real8(iGribMessageIndex,'latitudeOfFirstGridPointInDegrees',f8)
        f8 = corner_lat_N + (iyStart - 1) * dy_deg
        call grib_set_real8(iGribMessageIndex,'latitudeOfLastGridPointInDegrees', f8)


!        f8 = corner_lon_E + (ixEnd - 1) * dx_deg
!        call grib_set_real8(iGribMessageIndex,'longitudeOfLastGridPointInDegrees', f8)
        f8 = 0.001D0*nint(dx_deg*1000,8)
        call grib_set_real8(iGribMessageIndex,'iDirectionIncrementInDegrees',f8)
        f8 = 0.001D0*nint(dy_deg*1000,8)
        call grib_set_real8(iGribMessageIndex,'jDirectionIncrementInDegrees',f8)

        ! Note: The Scanning Mode is set to 0 [00000000] (see Flag Table 3.4) in the example data files
        ! This setting assumes that the ordering of the lat/lon points is as follows:
        !   - first point is the North-West corner point
        !   - last point is the South-East corner point
        !   - scanning is done along lines of latitude from West to East
        !
        ! Scanning mode that starts from the south:
        !
!        call grib_check(grib_set_int(iGribMessageIndex,"iScansNegatively",0))
!        call grib_check(grib_set_int(iGribMessageIndex,"jScansPositively",1))
!        call grib_check(grib_set_int(iGribMessageIndex,"jPointsAreConsecutive",0))
        !
        ! ... or from the north (mind the above flipping of the field):
        !
        call grib_set_int(iGribMessageIndex,"iScansNegatively",0)
        call grib_set_int(iGribMessageIndex,"jScansPositively",0)
        call grib_set_int(iGribMessageIndex,"jPointsAreConsecutive",0)
                
        ! Vertical
        !
        if (fu_leveltype(output_level) /= constant_height) then
          call set_error('Output level type not supported', sub_name)
          return
        end if
        if (fu_level_height(output_level) > 0) then
          ! Note: the level codes for grib2 differ from constants in levels module.
          call grib_set_int(iGribMessageIndex, 'typeOfFirstFixedSurface', 103)
          call grib_set_int(iGribMessageIndex, 'scaleFactorOfFirstFixedSurface', 0)
          call grib_set_int(iGribMessageIndex, 'scaledValueOfFirstFixedSurface', &
                          & int(fu_level_height(output_level)))
        else
          call grib_set_int(iGribMessageIndex, 'typeOfFirstFixedSurface', 1)
        end if

        ! Chemical constituent 
        ! 
        if(defined(species))then
          if (associated(nl_species_map)) then
            constituent_code = fu_content_int(nl_species_map, fu_content(nlIni, 'substance_name'))
          else
            constituent_code = int_missing
          end if
          if (constituent_code == int_missing) constituent_code = fu_content_int(nlIni, 'output_constituent_code')
          if (constituent_code == int_missing) then
            call set_error('Cannot determine constituent_code', sub_name)
            return
          end if
          call grib_set_int(iGribMessageIndex, 'constituentType', constituent_code)
        endif
        !
        ! Output parameter
        !
        parameter_code = fu_content_int(nlIni, 'output_parameter_code')
         ! 0 for mass concentration, 59 for number concentration 
        ! http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-2-0-20.shtm
        if (parameter_code == int_missing) parameter_code = fu_content_int(nl_grib2_defaults, 'parameterNumber')
        if (parameter_code == int_missing) then
           call set_error('Cannot determine output_parameter_code', sub_name)
           return
        end if
        call grib_set_int(iGribMessageIndex, 'parameterNumber', parameter_code)



        !  Update number of values for Section 5 of message
        !
        call grib_set_int(iGribMessageIndex,'numberOfValues',fsOut)

        !
        ! Finally, update the values and write the message
        !
        ! The slowest part
        call grib_set_real8_array(iGribMessageIndex,"values", output_data(1:fsout))

        !  Write modified message to a file
        !
        call grib_write(iGribMessageIndex, uGribOut)

        call grib_close_file(uGribOut)
                
        !
        ! Check the sum
        !
        fSum = sum(output_data(1:fsOut))
        write(unit=chMsg,fmt='(A,1x,4(F15.3,1x),2(I0,1x),2(F16.3,1x),E16.3)')&
             & 'Field:x0,y0,xE,yE,nx,ny,dx,dy,fSum', &
             & corner_lon_E + (ixStart - 1) * dx_deg, &
             & corner_lat_N + (iyEnd - 1) * dy_deg, &
             & corner_lon_E + (ixEnd - 1) * dx_deg, &
             & corner_lat_N + (iyStart - 1) * dy_deg, &
             & ixEnd - ixStart + 1, &
             & iyEnd - iyStart + 1, &
             & dx_deg, &
             & dy_deg, &
             & fSum
        !call msg(chMsg)
        call msg('The output field sum is:', real(fSum))

!        !--------------------------------------------------------------------------------
!        !
!        ! Now let's check immediately what has been just written
!        ! NO SENSE: for multi-message file below stuff just takes the first message
!        !
!        output_data = 0.
!        call grib_check(grib_open_file(uGribOut,trim(chTmp),'r'))
!        call grib_check(grib_new_from_file(uGribOut, iCheckGribMessage))
!        call grib_check(grib_get_real8_array(iCheckGribMessage,"values", output_data, fs))
!        call msg('Data written and re-read. Actual and read sums:',int(fSum+0.5),real(sum(output_data(1:fs))))
!        call grib_check(grib_release(iCheckGribMessage))
!        call grib_check(grib_close_file(uGribOut))
                
      end do   ! time

      !
      ! Release field that is just written
      !
      call grib_release(iGribMessageIndex)

      call msg('Next namelist')

    end do ! cycle over namelists <==> data sets

    !Close input file
    if (trim(OpenedFileName) /= '') then
          call msg("Closing input file: "//trim(OpenedFileName))
          if (FForm%iformat == netcdf_file_flag) then
                 call close_netcdf_file(iGF)
         else
                  call close_gradsfile_i(iGF)
         endif
    endif

    call free_work_array(grads_data)
    
  end subroutine convert_2_grib2

  subroutine update_grib2_defaults(ind_grib, nl_defaults)
    implicit none
    integer, intent(in) :: ind_grib
    type(Tsilam_namelist), pointer :: nl_defaults

    integer :: num_items, ind_item, int_content
    
    if (.not. associated(nl_defaults)) return
    
    do ind_item = 1, fu_nbr_of_items(nl_defaults)
      int_content = fu_content_int(fu_get_item(nl_defaults, ind_item))
      if (error) return
      if (fu_fails(int_content /= int_missing, 'Bad int_content', 'update_grib2_defaults')) return
      !call msg('set:' // trim(fu_name(fu_get_item(nl_defaults, ind_item))), int_content)
      call grib_set_int(ind_grib, fu_name(fu_get_item(nl_defaults, ind_item)), int_content)
      if  (error) return
    end do
        
  end subroutine update_grib2_defaults




  subroutine find_grid_ind(val_to_cover, val_first, step, num_vals, start_or_end, ind_found)
    implicit none
    real, intent(in) :: val_to_cover, val_first, step
    integer, intent(in) :: num_vals
    character(len=*), intent(in) :: start_or_end
    integer, intent(out) :: ind_found

    integer :: ii
    real :: val_next, val_cur

    if (val_to_cover .eps. real_missing) then
      ind_found = int_missing
      return
    end if


    if (start_or_end == 'start') then
      do ii = 1, num_vals
        if (val_to_cover < val_first + ii*step) exit
      end do
      ind_found = ii

    else if (start_or_end == 'end') then
      do ii = 1, num_vals
        if (val_to_cover <= val_first + ii*step) exit
      end do
      ind_found = ii + 1

    else
      call set_error('Bad start_or_end', 'find_grid_ind') 
      return
    end if
    if (ind_found > num_vals) then
      call set_error('Cannot cover requested area', 'find_grid_ind')
      return
    end if

  end subroutine find_grid_ind

end program silam_2_grib2

