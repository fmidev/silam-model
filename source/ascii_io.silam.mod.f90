module ascii_io
  !
  ! This module is responsible for all non-trivial input-output operations
  ! with ASCII files (text files). In particular, it reads various sources,
  ! which are defined here.
  ! 
  ! Author: Mikhail Sofiev,  mikhail.sofive@fmi.fi
  !
  ! All units: SI
  ! 
  ! All degrees: real numbers in degrees (degrees and parts per
  ! hundred/parts per thousand , NOT in degrees and minutes
  ! 
  ! In this type we define:
  ! NORTHERN LATITUDES ARE POSITIVE
  ! EASTERN LONGITUDES ARE POSITIVE
  !
  ! Language: ANSI standard Fortran 90
  ! 
  !use fields

  use stacks 

  implicit none
  
  private

  public read_next_field_from_ascii_file
  public open_ascii_file_i
  public close_ascii_file_i

  CONTAINS


  !**********************************************************************

  subroutine open_ascii_file_i(chFNm, iUnit)
    !
    ! Just opens the file. Even md module is not needed
    !
    implicit none

    ! Imported parametere
    character(len=*), intent(in) :: chFNm
    integer, intent(out) :: iUnit

    ! Local variables
    integer :: iStatus

    iUnit = fu_next_free_unit()

    open(unit=iUnit, file=chFNm, status='old', iostat=iStatus)

    if(iStatus /= 0)then
      call set_error('Failed to open file:' + chFNm,'open_ascii_file_i')
      iUnit = int_missing
    endif

  end subroutine open_ascii_file_i


  !**********************************************************************

  recursive subroutine read_next_field_from_ascii_file(uFile, eof, id, fieldData, fFillValue)
    !
    ! The file is supposed to be similar to GRIB but ASCII rather than binary
    ! A lot of common exists with area source but time and values are somewhat 
    ! different. 
    ! The whole file is just a namelist, which must contain enough information for 
    ! defining the grid, level/layer, time and a set of values.
    ! This routine sends the namelist to the field_id module in order to create the id
    ! and then takes care of the array of data
    !
    implicit none

    ! Imported parameters
    integer, intent(in) :: uFile
    logical, intent(out) :: eof
    type(silja_field_id), intent(out) :: id
    real, dimension(:), intent(out), optional :: fieldData
    real, intent(in), optional :: fFillValue

    ! Local variables
    type(Tsilam_namelist), pointer :: nlField
    character(len=fnlen) :: spContent
    integer :: iLocal, iTmp, nVals, nx, ny, iStatus, ix, iy
    real :: fLon, fLat, fx, fy, fVal, fScaling, fMissingValue
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: ptrItems
    logical :: ifGeoCoords
!    type(silam_grid_position) :: posTmp
    type(silja_grid) :: grid

        !
    ! Find the starting line in the ASCII file
    !
    eof = .false.
    do while (.not.eof)
      call next_line_from_input_file(uFile, spContent, eof)
      if(error .or. index(spContent,'ASCII_FIELD_1') /= 0) exit
    end do
    if(eof.or.error)return

    !
    ! Read the next field to the namelist. If data are not needed, skip them
    !
    if(present(fieldData))then
      nlField => fu_read_namelist(uFile, .false.,'END_ASCII_FIELD_1')
    else 
      nlField => fu_read_namelist(uFile, .false., 'END_ASCII_FIELD_1', chIgnoreItem='val')
    endif

    if(error)return
    if(empty(nlField))then
      eof = .true.
      return
    else
      eof = .false.
    endif

    !-------------------------------------------------------------
	!
	! A trick: it may happen that the field is, in fact, a separate file
	! Read it then
	!
    spContent = fu_content(nlField,'source_file_name')
	if(spContent /= '')then
        iLocal = fu_next_free_unit()
        open(unit=iLocal, file=spContent, status='old', iostat=iTmp)
        if(iTmp == 0)then
          call read_next_field_from_ascii_file(iLocal, eof, id, fieldData) !, fMissingValue)
        else
          call set_error('Field file does not exist:' + spContent, &
                       & 'read_next_field_from_ascii_file')
          call set_missing(id)
        endif
        call destroy_namelist(nlField)
        return
	endif

    !------------------------------------------------------------
    !
    ! Create a field from the namelist
    !
    id = fu_set_field_id(nlField)
    if(.not. defined(id) .or. error)then
      call set_error('Failed to determine the field id','read_next_field_from_ascii_file')
      call destroy_namelist(nlField)
      return
    endif
    call msg(fu_connect_strings('Found:',fu_quantity_string(fu_quantity(id))))

    !------------------------------------------------------------
    !
    ! Check the id and, if it is reasonable, fill the array with the missing values 
    ! and then store the meaningfull values from the namelist
    !
    fMissingValue = fu_content_real(nlField,'missing_value')
    if(error)then
      fMissingValue = real_missing
      call unset_error('read_next_field_from_ascii_file')
    endif

    if(present(fieldData)) then
      !
      ! Sufficiently large to receive the data?
      !
      if(size(fieldData) < fu_number_of_gridpoints(fu_grid(id)))then
        call msg('Input grid bigger then work array size:', fu_number_of_gridpoints(fu_grid(id)), &
                                                          & size(fieldData))
        call set_error('Too large input grid. GrADS or NetCDF formats must be used', &
                     & 'read_next_field_from_ascii_file')
        return
      endif
      !
      ! Get the scaling factor - either via unit or scaling_factor lines
      !
      if(len_trim(fu_content(nlField,'unit')) > 0)then
        fScaling = fu_set_named_value('1.0  ' // fu_content(nlField,'unit'))
        if(len_trim(fu_content(nlFIeld,'scaling_factor')) > 0)then
          call set_error('Both scaling_factor and unit are given','read_next_field_from_ascii_file')
          return
        endif
      elseif(len_trim(fu_content(nlField,'scaling_factor')) > 0)then
        fScaling = fu_content_real(nlField,'scaling_factor')
      else
        fScaling = 1.0
      endif
      if(error .or. (fScaling .eps. real_missing))then
        call set_error('cannot get the scaling factor from the namelist', 'read_next_field_from_ascii_file')
        return
      endif
      !
      ! and the grid
      !
      grid = fu_grid(id)

      if(present(fFillValue))then
        fieldData(:) = fFillValue
      else
        fieldData(:) = fMissingValue
      endif

      nullify(ptrItems)
      call get_items(nlField, 'val', ptrItems, nVals)

      call grid_dimensions(grid, nx, ny)
      if(nx*ny > size(fieldData))then
        call report(id)
        call msg('Too small size of data array:',size(fieldData))
        call set_error('Too small size of the data array','read_next_field_from_ascii_file')
        call destroy_namelist(nlField)
        return
      endif

      ifGeoCoords = fu_content(nlField,'coordinate_of_values') == 'GEOGRAPHICAL'

      do iTmp = 1, nVals
        spContent = fu_content(ptrItems(iTmp))

        if(ifGeoCoords)then
          !
          ! Have to reproject the val line to the source grid
          !
          read(unit=spContent, iostat=iStatus,fmt=*) fLon, fLat, fVal
          if(iStatus /= 0)then
            call set_error(fu_connect_strings('Cannot read val string:',spContent), &
                         & 'read_next_field_from_ascii_file')
            cycle
          endif
          if(fVal .eps. fMissingValue)cycle

          call project_point_to_grid(fLon, fLat, grid, fX, fY)
          if(error)return

!          posTmp = fu_set_pos_in_geo_global_grid(fLon,fLat,100000.,fu_valid_time(id))
!          call project_point_to_grid(geo_global_grid, fu_x(posTmp), fu_y(posTmp), &
!                                   & fu_grid(id), fx, fy)
          !
          ! Check that the reprojected point is indeed inside the source grid
          !
          if(fx < 0.5 .or. fy < 0.5 .or. fx > nx+0.5 .or. fy > ny+0.5)then
            call msg(fu_connect_strings('Geo-coord cell:',spContent, &
                                      & ', grid co-ordinates:'), int(fx+0.5), fy)
            call msg_warning(fu_connect_strings('The above cell of area field:', &
                                              & fu_quantity_short_string(fu_quantity(id)), &
                                              & '-is outside the claimed grid'), &
                           & 'read_next_field_from_ascii_file')
            call report(grid)
            cycle
          endif
        else
          !
          ! No reprojection, the value are given in grid coordinates
          !
          read(unit=spContent, iostat=iLocal,fmt=*) fx, fy, fVal
          if(iLocal /= 0)then
            call set_error(fu_connect_strings('Cannot read val string:',spContent), &
                        & 'read_next_field_from_ascii_file')
            cycle
          endif
          if(fVal .eps. fMissingValue)cycle
          !
          ! Check that the point is inside the source grid
          !
          if(fx < 0.5 .or. fy < 0.5 .or. fx > nx+0.5 .or. fy > ny+0.5)then
            call msg(fu_connect_strings('Geo-coord cell:',spContent, &
                                      & ', grid co-ordinates:'), int(fx+0.5), fy)
            call msg_warning(fu_connect_strings('The above cell of area field:', &
                                              & fu_quantity_short_string(fu_quantity(id)), &
                                              & '-is outside the claimed grid'), &
                           & 'read_next_field_from_ascii_file')
            call report(fu_grid(id))
            cycle
          endif
        endif  ! ifGeoCoords

        fieldData(int(fx+0.5)+(int(fy+0.5)-1)*nx) = fVal * fScaling
          
      end do ! cycle over val strings

      call destroy_items(ptrItems)

    endif ! present fieldData
    !
    ! Destroy the temporary string and the namelist 
    !
    call destroy_namelist(nlField)

  end subroutine read_next_field_from_ascii_file


  !**********************************************************************

  subroutine close_ascii_file_i(iUnit)
    !
    implicit none

    ! Imported parametere
    integer, intent(in) :: iUnit

    ! Local variables
    integer :: iStatus

    close(unit=iUnit, iostat=iStatus)

    if(iStatus /= 0)then
      call msg('Failed to close unit:', iUnit)
      call set_error('Failed to close unit','close_ascii_file_i')
    endif

  end subroutine close_ascii_file_i


end module ascii_io
