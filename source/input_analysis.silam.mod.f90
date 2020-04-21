MODULE input_analysis
!
! Contains necessary structures and routines for the reading-through (via grib_io module
! routines) and analysis of the content of a grib file.
! The task of this module is - to scan the file and advise about the best choice
! for the model system grid, the most appropriate vertical structure, etc.
! Previous approach "take-the-first" does not work if GRIB file contains data in
! different grids, vertical structures, etc. An example of such files is - ECMWF
! operational ones. With such files one has to be careful selecting the system grid
! in order to minimise interpolations and possibly avoid them before the derived
! quantities are made.
!
! IMPORTANT. The structure does not deal with (sometimes) clumsy GRIB parameters.
!            Instead, SILAM-own classification (also sometimes clumsy) is used.
! Be also careful if the GRIB files are different (e.g., from two or more sources)
! or can be different for different time steps. In this case selection made at the
! beginning may be not optimal later.
!
! Language: FORTRAN-90
!
! All units: SI
!
! Author: M.Sofiev
!
!use grib_io
use derived_field_quantities

public read_input_content ! Reads file, fills content structure
public analyse_input_content ! Analyses the GRIB file content and gives recomendations

public fu_meteo_grid    ! Selects the system grid
public fu_meteo_vertical  ! Advises about the vertical structure from meteofields
public fu_shopping_list   ! INFO. Picks a sub-set of shopping-list quantities available
public fu_nbr_of_grid_types ! INFO. Tells the number of different grid types
public fu_nbr_of_verticals ! INFO. Tells the number of different vertical structures
public fu_nbr_of_fields
public set_met_srcs  ! INFO. Drops the list of the meteo sources in the wdr structure
public id_lst

private check_metsrcs_list_size
private compact_ids  ! Eliminates the holes in id list

private count_grids  ! Counts the number of grids and fills the Tgrid_lst list
private exchange_grd_lsts ! Exchanges two Tgrid_lst objects
private check_coverage ! Computes the coverage of all counted grids
private remove_grid  ! Removes grid and corresponding fields from the structure
private fu_gridPtr ! Finds the index of the given grid in the grid list
private check_gridlist_size ! Memory allocation routine

private count_data_sources ! Counts data sources
private exchange_met_src_lsts ! exchanges two TMet_src_lst objects
private remove_met_src 
!private check_MetSrc_list_size ! Memory allocation
private fu_MetSrcPtr ! Finds the index of the data source in the data source list

private count_verticals
private check_vertlist_size ! Memory allocation
private exchange_vert_lsts
private fu_vertPtr

private count_variables ! Counting variables in the content.
private check_varlist_size ! Memory allocation
private fu_compare_vars_quality ! Compares quality parameters of two variables
private remove_variable  ! And all related ids


!---------------------------------------------------------------------
!
! Structures describing the GRIB file content.
!
real, private, parameter :: coverage_threshold = 0.01
!
! Elements of lists of data sources, grids, verticals and variables
!
! grid
!
type Tgrid_lst 
  private 
  type(silja_grid) :: grid ! grid 
  integer :: NFlds_grd  ! Number of fields having one of these grids
  real :: out_cover, cell_size ! the % of the target area covered and the grid cell size
  integer, dimension(max_2d_fields) :: flds ! Indices of the fields with this grid
end type Tgrid_lst

!
! vertical
!
type Tvert_lst
  private
  type(silam_vertical) :: vert
  real :: top=0.
!  integer :: n_layers=0
!  type(silja_level), dimension(max_levels) :: levs
  integer :: gridPtr, srcPtr ! Vertical belongs to source and grid
  integer :: NFlds_vert=0  ! Number of fields with this vertical
  integer, dimension(max_2d_fields) :: flds=0 ! Indices of fields with this vertical
end type Tvert_lst


! data source
!
type Tmet_src_lst
  private
  type(meteo_data_source) :: src ! The data source
  integer :: NFlds_src ! The number of the fields having this source
  real :: best_cell_size ! The smallest available size of the cell
  integer, dimension(max_2d_fields) :: flds ! Indices of the fields from this source
end type Tmet_src_lst

!
! variable. One variable is a combination of quantity, grid, vertical and data source
! This combination must be unique
!
type Tvar_lst
  private
  integer :: quantity
  integer :: NFlds_var
  integer, dimension(max_2d_fields) :: flds
  integer :: gridPtr, MetSrcPtr, vertPtr, vertLevNbr
end type Tvar_lst



type Tinput_content ! Full content of the GRIB file(s) plus results of analysis
  private
  integer :: NbrFields=0, NbrGrids=0, NbrVerts=0, NbrTimes=0, NbrMetSrcs=0, NbrVars=0
  integer :: SystemGridPtr=0, SystemLevPtr=0
  type(Tgrid_lst),dimension(:),pointer :: grids => null() ! List of grids
  type(Tvert_lst),dimension(:),pointer :: verts => null() ! List of vertical co-ordinates
  type(Tmet_src_lst), dimension(:), pointer :: MetSrcs => null()! List of meteodata sources
  type(Tvar_lst), dimension(:), pointer :: vars => null()! List of variables (quantities)
  type(silja_logical) :: defined=silja_false
  logical :: analysed=.false. ! Whether the content was analysed
  type(silja_field_id),dimension(:),pointer :: ids => null() ! List of fields in GRIB
  type(silja_shopping_list) :: shopping_list
end type Tinput_content

type (Tinput_content), parameter, public :: input_content_missing = &
   & Tinput_content( 0,0,0,0,0,0,0,0,&
                     &null(), null(), null(), null(), silja_false, .false.,  &
                     & null(), shopping_list_missing )


CONTAINS



!*********************************************************************

subroutine read_input_content(inputName, wdr, inputFormat, q_list, InputContent)
  !
  ! Reads an input, fills the content structure. 
  ! ATTENTION!! No analysis is done here, just structure is filled.
  ! If the structure is not empty - it is extended, NOT overwritten
  ! Only quantities from the shopping list q_list are considered,
  ! the others are ignored
  !
  ! Code owner Mikhail Sofiev, FMI
  !
  implicit none
  
  ! Imported parameters
  character(len=*), intent(in) :: inputName ! Name of the GRIB file 
  type(silam_fformat), intent(in) :: inputFormat ! GRIB file, ASCII file, TEST_FIELD, .....
  type(silja_wdr), intent(in) :: wdr
  type(silja_shopping_list), intent(in) :: q_list ! List of quantities to consider
  type(Tinput_content), intent(inout) :: InputContent ! GRIB-content structure


  ! Local variables
  type(silja_grid) :: gridTmp
  integer :: j

  !
  ! Depending on the file format, we should read some file or create a test field
  !

  select case(inputFormat%iformat)
    case(grib_file_flag)
      call get_content_from_grib_file(inputName, wdr, q_list, InputContent)
    case(ascii_file_flag)
      call get_content_from_ascii_file(inputName, q_list, InputContent)
    case(netcdf_file_flag)
      call get_content_from_netcdf_file(inputName, inputFormat, q_list, InputContent)
    case(grads_file_flag)
      call get_content_from_grads_file(inputName, q_list, InputContent)
    case(test_field_value_flag)
      call put_test_field_to_content(inputName, q_list, InputContent)
    case default
      call msg('Unknown input format:',inputFormat%iformat)
      call set_error('Unknown input format','read_input_content')
      return
  end select

  if(InputContent%NbrFields < 1)then
    call msg_warning(fu_connect_strings('No usable fields found in:',inputName), &
                   & 'read_input_content')
    return
  end if

  do while(.not.defined(InputContent%ids(InputContent%NbrFields)))
    InputContent%NbrFields = InputContent%NbrFields - 1
    if(InputContent%NbrFields < 1)then
      call msg_warning(fu_connect_strings('No usable fields found in:',inputName), &
                     & 'read_input_content')
      return
    end if
  end do


  ! There can be one special case: a global grid determined within the longitude 
  ! range (0:360), while SILAM requires (-180:180).
  ! This has to be dealt with via shifting the starting point of the meteo_grid.
  ! As a result, each new field will be reprojected (actually, just renumbered)
  ! during the file reading. Not too nice but the only way to handle the situation
  ! without major rebuilding the model.
  ! There is another point where this has to be repeated: when getting fields with
  ! such a grid, the field has to be re-arranged. This is handled within 
  ! grid_data_horizontal_select
  !
  do j=1, InputContent%NbrFields
    gridTmp = fu_grid(InputContent%ids(j))
    if(fu_ifLonGlobal(gridTmp))then
      if(.not.fu_stdSilamGrid(gridTmp))then
!        call msg('repositioning the grid along the longitude')
!        call report(gridTmp)
        call reposition_global_grid_lon(gridTmp)
      endif
    endif
    if(fu_ifLatGlobal(gridTmp))then
      if(.not.fu_stdSilamGrid(gridTmp))then
!        call msg('repositioning the grid along the longitude')
!        call report(IC%grids(IC%SystemGridPtr)%grid)
        call reposition_global_grid_lat(gridTmp)
      endif
    endif
    if(.not.fu_stdSilamGrid(gridTmp))then
      call set_error('The meteo grid is not SILAM-standard and cannot be made to comply', &
                   & 'read_input_content')
      call report(gridTmp)
      return
    endif
    call set_grid(InputContent%ids(j), gridTmp)
  end do  ! cycle over all IDs


end subroutine read_input_content

  !============================================================================

  subroutine get_content_from_grib_file(FName, wdr, q_list, InputContent)
    !
    ! Reads a GRIB file and fills in the InputContent
    !
    implicit none

    character(len=*), intent(in) :: FName ! Name of the GRIB file 
    type(silja_shopping_list), intent(in) :: q_list ! List of quantities to consider
    type(silja_wdr), intent(in) :: wdr
    type(Tinput_content), intent(inout) :: InputContent ! GRIB-content structure

    ! Local variables
    integer :: grib_unit, io_status, nIDs, iID, quantity, iTmp
    logical :: ifAccept
    real :: scale_factor
    logical :: ifZeroNegatives
    TYPE(silja_field_id), dimension(:),pointer :: IDlist

    !
    ! Open the GRIB file and other preparations
    !
    call open_grib_file_i(fname, grib_unit, fu_obstime_interval(wdr))
    if (error) return
    !
    ! Cycle through the gribfile.
    !
    call id_list_from_grib_file(grib_unit,IDlist, nIDs)

    InputContent%analysed = .false.

    do iID = 1,nIDs
      !
      ! Prepare the space for the next field and fill it in
      !


#ifdef DEBUG
call msg('ID from file:', iID)
call report(IDlist(iID))
#endif

      if (.not. defined(IDlist(iID))) cycle

      quantity=fu_quantity(IDlist(iID))
      !
      ! Conditions for the field to be in the GribContent structure are:
      ! - reasonable field
      ! - field quantity is in the shopping list
      ! - field has not been stored already e.g. from another input file
      !
      if ( quantity ==  int_missing) cycle
      if (.not. fu_quantity_in_list(quantity,q_list)) cycle

      !
      ! Try to find a similar ID to avoid duplications
      do iTmp = 1,InputContent%NbrFields
        if(InputContent%ids(iTmp) == IDlist(iID)) exit
      end do
      if (iTmp <= InputContent%NbrFields) cycle !loop above was exited

      InputContent%defined = silja_true  ! Something is there

      InputContent%NbrFields = InputContent%NbrFields + 1

      call chk_content_size(InputContent, InputContent%NbrFields)
      if(error) return

      InputContent%ids(InputContent%NbrFields) = IDlist(iID)
#ifdef DEBUG
      call msg('Added as IC', InputContent%NbrFields)
#endif
      
    END DO ! loop_over_fields

    ! Final arrangements and exit.
    !
    CALL close_grib_file_i(grib_unit)

  end subroutine get_content_from_grib_file


  !============================================================================

  subroutine get_content_from_ascii_file(FName, q_list, InputContent)
    !
    ! Reads a ASCII file and fills in the InputContent
    !
    implicit none

    ! Imported parameters
    !
    character(len=*), intent(in) :: FName ! Name of the GRIB file 
    type(silja_shopping_list), intent(in) :: q_list ! List of quantities to consider
    type(Tinput_content), intent(inout) :: InputContent ! GRIB-content structure

    ! Local variables
    !
    integer :: iUnit, iID
    logical :: ifAccept, eof
    real :: fMissingValue

    !
    ! Open the GRIB file and other preparations
    !
    call open_ascii_file_i(FName, iUnit)
    if (error) return

    !
    ! Cycle through the asciifile.
    !
    eof=.false.
    DO while (.not.eof)
      !
      ! Prepare the space for the next field and fill it in
      !
      InputContent%NbrFields = InputContent%NbrFields + 1

      call chk_content_size(InputContent, InputContent%NbrFields)
      if(error) return
      !
      ! Get the new field
      !
      call read_next_field_from_ascii_file(iUnit, eof, InputContent%ids(InputContent%NbrFields))

      if(error)then
        call unset_error('get_content_from_ascii_file')
        InputContent%NbrFields = InputContent%NbrFields - 1
      else
        !
        ! Conditions for the field to be in the GribContent structure are:
        ! - reasonable field
        ! - field quantity is in the shopping list
        ! - field has not been stored already e.g. from another input file
        !
        ifAccept = fu_quantity(InputContent%ids(InputContent%NbrFields)) /= int_missing .and. &
             & fu_quantity_in_list(fu_quantity(InputContent%ids(InputContent%NbrFields)),q_list)
        if(ifAccept)then
          do iID = 1,InputContent%NbrFields-1
            if(InputContent%ids(iID) == InputContent%ids(InputContent%NbrFields))then
              ifAccept = .false.
              exit
            endif
          end do
        endif
        if(ifAccept)then
          InputContent%defined = silja_true
          InputContent%analysed = .false.
        else
          call set_missing(InputContent%ids(InputContent%NbrFields))
          InputContent%NbrFields = InputContent%NbrFields - 1
        end if
      end if   ! error

    END DO ! loop_over_fields

    ! Final arrangements and exit.
    !
    CALL close_ascii_file_i(iUnit)

  end subroutine get_content_from_ascii_file

  !===========================================================================

  subroutine get_content_from_netcdf_file(FName, fformat, q_list, InputContent)
    !
    ! Gets input content from the NetCDF file header
    !
    implicit none

    character(len=*), intent(in) :: FName
    type(silam_fformat) :: fformat
    type(silja_shopping_list), intent(in) :: q_list ! List of quantities to consider
    type(Tinput_content), intent(inout) :: InputContent ! GRIB-content structure

    ! Local variables
    !
    integer :: iUnit, iID, iVar, nVars
    logical :: ifAccept
    type(silja_field_id),dimension(:), pointer :: idList

    !
    ! Open the file and other preparations
    !
    iUnit = open_netcdf_file_i(FName, fformat)
    if (error) return

    !
    ! Cycle through the file.
    !
    call id_list_from_netcdf_file(iUnit, idList, nVars)
    if (error) return
    call chk_content_size(InputContent, inputContent%nbrFields + nVars)
    if (error) return

!    InputContent%NbrFields = 0
    DO iVar = 1, nVars

        ! Conditions for the field to be in the GribContent structure are:
        ! - reasonable field
        ! - field quantity is in the shopping list
        ! - field has not been stored already e.g. from another input file
     

      if(fu_quantity(idList(iVar)) /= int_missing .and. fu_quantity_in_list(fu_quantity(idList(iVar)),q_list))then
        ifAccept = .true.
        do iID = 1,InputContent%NbrFields-1
          if(InputContent%ids(iID) == idList(iVar))then
            ifAccept = .false.
            exit
          endif
        end do
        if(ifAccept)then
          InputContent%NbrFields = InputContent%NbrFields + 1
          InputContent%ids(InputContent%NbrFields) = idList(iVar)
          InputContent%defined = silja_true
          InputContent%analysed = .false.
        end if
      end if   ! error

    END DO ! loop_over_fields

    ! Final arrangements and exit.
    !
    CALL close_netcdf_file(iUnit)
    deallocate(idList)

  end subroutine get_content_from_netcdf_file


  !============================================================================
  subroutine get_content_from_grads_file(FName, q_list, InputContent)
    !
    ! Gets an input content from the grads file - actually, from ctl/super-ctl
    !
    implicit none

    character(len=*), intent(in) :: FName
    type(silja_shopping_list), intent(in) :: q_list ! List of quantities to consider
    type(Tinput_content), intent(inout) :: InputContent ! GRIB-content structure

    ! Local variables
    !
    integer :: iUnit, iID, iVar, nVars, nLevs, iLev
    logical :: ifAccept
    type(silja_field_id) :: id

    !
    ! Open the file and other preparations
    !
    iUnit = fu_open_gradsfile_i(FName)
    if (error) return
    !
    ! Cycle through the file.
    !
    nVars = fu_n_gvars(iUnit)
    nLevs = fu_n_glevs(iUnit)
   
    iLev = 0
    do iVar = 1, nVars
      iLev = iLev + fu_n_gVar_levs(iUnit, iVar)     !temporary use of variable iLev
    enddo

    call chk_content_size(InputContent, iLev)
    if (error) return

!    InputContent%NbrFields = 0
    DO iVar = 1, nVars
      do iLev = 1, fu_n_gVar_levs(iUnit, iVar)
        call get_grads_var_metadata(iUnit, iVar, iLev, 1, id) 

        ! Conditions for the field to be in the GribContent structure are:
        ! - reasonable field
        ! - field quantity is in the shopping list
        ! - field has not been stored already e.g. from another input file
     
        if(fu_quantity(id) /= int_missing .and. fu_quantity_in_list(fu_quantity(id),q_list))then
          ifAccept = .true.
          do iID = 1,InputContent%NbrFields-1
            if(InputContent%ids(iID) == id)then
              ifAccept = .false.
              exit
            endif
          end do
          if(ifAccept)then
            InputContent%NbrFields = InputContent%NbrFields + 1
            InputContent%ids(InputContent%NbrFields) = id
            InputContent%defined = silja_true
            InputContent%analysed = .false.
          end if
        end if   ! error
        if(.not. fu_multi_level_quantity(fu_quantity(id)))exit
      enddo
    END DO ! loop_over_fields

    ! Final arrangements and exit.
    !
    CALL close_gradsfile_i(iUnit)

  end subroutine get_content_from_grads_file


  !============================================================================

  subroutine put_test_field_to_content(chFieldDescr, q_list, InputContent)
    !
    ! Just creates one more ID in the InputContent. No data arrays are filled.
    ! Attention. The system_grid is not yet known, so we have to ensure that 
    ! this field does not affect anything in terms of grid selection or other
    ! parameters. Its grid will be global and with a very coarse resolution.
    !
    implicit none

    ! Imported parameters
    !
    character(len=*), intent(in) :: chFieldDescr ! Descriptor of the field
    type(silja_shopping_list), intent(in) :: q_list ! List of quantities to consider
    type(Tinput_content), intent(inout) :: InputContent ! GRIB-content structure

    ! Local variables
    !
    integer :: iStatus, quantity
    type(silam_sp) :: spTmp
    type(silja_level) :: level

    spTmp%sp => fu_work_string()
    !
    ! Let's start from quantity
    !
    read(unit=chFieldDescr, iostat=iStatus, fmt=*) spTmp%sp
    if(iStatus /= 0)then
      call msg_warning(fu_connect_strings('Failed to read quantity name from:',chFieldDescr), &
                     & 'put_test_field_to_content')
      return
    endif

    quantity = fu_get_silam_quantity(spTmp%sp)
    !
    ! If the quantity is strange or not interesting 
    if(quantity == int_missing .or. .not. fu_quantity_in_list(quantity,q_list)) return

    !
    ! Now set the level
    !
    if(index(chFieldDescr,'SURFACE_LEVEL') > 0)then
      level = surface_level
    elseif(index(chFieldDescr,'2M_LEVEL') > 0)then
      level = level_2m_above_ground
    elseif(index(chFieldDescr,'10M_LEVEL') > 0)then
      level = level_10m_above_ground
    else
      call msg('Allowed levels: SURFACE_LEVEL, 2M_LEVEL, 10M_LEVEL')
      call msg_warning(fu_connect_strings('Unknown level in:',chFieldDescr), &
                     & 'put_test_field_to_content')
      return
    endif
    !
    ! Prepare the space for the field and set the id
    !
    InputContent%NbrFields = InputContent%NbrFields + 1

    call chk_content_size(InputContent, InputContent%NbrFields)
    if(error) return
    !
    ! Get the new field
    !
    InputContent%ids(InputContent%NbrFields) = fu_set_field_id(silam_internal_src, &
                                                             & quantity, &
                                                             & fu_start_time(q_list), &
                                                             & zero_interval, &
                                                             & coarse_geo_global_grid, &
                                                             & level)
    call free_work_array(spTmp%sp)

  end subroutine put_test_field_to_content


  !============================================================================

  subroutine chk_content_size(InContent, iSize)
    !
    ! Checks the size of the local GribContent and enlarges it if needed
    !
    implicit none

    type(Tinput_content), intent(inout) :: InContent ! GRIB-content structure
    integer, intent(in) :: iSize
    !
    ! Local variables
    !
    integer :: iSize2alloc, iStat

    type(silja_field_id),dimension(:),allocatable :: idTmp

    if(.not. fu_true(inContent%defined) .or. .not.associated(InContent%ids))then
      allocate(InContent%ids(10*(int(real(iSize)/10.)+1)),stat=iStat)
      if(iStat /= 0)then
        call set_error('No more space available','chk_content_size')
        return
      end if
      return
    end if

    if(iSize > size(InContent%ids))then
      !
      ! Create temporary space and store the content there
      !
      allocate(idTmp(size(InContent%ids)),stat=iStat)
      if(iStat /= 0)then
        call set_error('No more space available','chk_content_size')
        return
      end if
      idTmp(:)=InContent%ids

      deallocate(InContent%ids)
      iSize2alloc = 10*(int(real(iSize)/10.)+1)
      allocate(InContent%ids(iSize2alloc),stat=iStat)
      if(iStat /= 0)then
        call set_error('No more space available','chk_content_size')
        return
      end if
      InContent%ids(1:size(idTmp)) = idTmp(1:size(idTmp))
      deallocate(idTmp)

    end if ! iSize > existing content size

  end subroutine chk_content_size



!*********************************************************************

function fu_meteo_grid(GribContent) result(grid)
!
! Returns the advisable system grid. Makes an analysis if necessary
!
  implicit none

  ! Result value
  type(silja_grid)::grid

  ! Imported parameter
  type(Tinput_content), intent(inout) :: GribContent

  grid = grid_missing
  if(.not.(GribContent%defined == silja_true)) then
    call set_error('Sorry, empty content','fu_system_grid')
    return
  end if

  if(.not.GribContent%analysed)then
    call set_error('The content is not yet analysed','fu_system_grid')
    return
  end if

  if(.not.error) grid = GribContent%grids(GribContent%SystemGridPtr)%grid

end function fu_meteo_grid


!*********************************************************************

function fu_meteo_vertical(GribContent) result (vertical)
!
! Returns an advisable vertical structure type.  Makes an analysis if necessary
!
  implicit none

  ! Returned value
  type(silam_vertical) :: vertical

  ! Imported parameter with intent INOUT
  type(Tinput_content), intent(inout) :: GribContent

  call set_missing(vertical, .true.)

  if(.not.(GribContent%defined == silja_true)) then
    call set_error('Sorry, empty content','fu_best_vertical')
    return
  end if

  if(.not.GribContent%analysed)then
    call set_error('The content is not yet analysed','fu_meteo_vertical')
    return
  end if

  if(.not.error) vertical = GribContent%verts(GribContent%SystemLevPtr)%vert

end function fu_meteo_vertical


!*********************************************************************

function fu_shopping_list(GribContent) result(list)
!
! Information only. Returns a shopping_list created during the analysis of 
! the GRIB content structure
!
implicit none

  ! Return value
  type(silja_shopping_list) :: list

  ! Imported parameter
  type(Tinput_content), intent(in) :: GribContent

  ! local variables
  integer :: i
  integer, dimension(2) :: qTmp
  
  qTmp = int_missing

  if(.not.(GribContent%defined == silja_true)) then
    call set_error('Empty content','fu_list_to_shop')
    call set_missing(list)
    return
  end if

  if(.not.GribContent%analysed) then
    call set_error('Content is not yet analysed','fu_shopping_list')
    call set_missing(list)
    return
  end if

  list = GribContent%shopping_list

end function fu_shopping_list


!*********************************************************************

integer function fu_nbr_of_grid_types (GribContent)
!
! Information only. Tells the number of different grid types
!
  implicit none
  
  ! Imported parameter
  type(Tinput_content), intent(in) :: GribContent

  if(.not.(GribContent%defined == silja_true)) then
    call set_error('Sorry, empty content','fu_nbr_of_grid_types')
    fu_nbr_of_grid_types = -1
    return
  end if

  fu_nbr_of_grid_types = GribContent%NbrGrids

end function fu_nbr_of_grid_types 


!*********************************************************************

integer function fu_nbr_of_verticals (GribContent)
!
! Information only. Tells the number of different vertical structures
!
  implicit none
  
  ! Imported parameter
  type(Tinput_content), intent(in) :: GribContent

  if(.not.(GribContent%defined == silja_true)) then
    call set_error('Sorry, empty content','fu_nbr_of_verticals')
    fu_nbr_of_verticals = -1
    return
  end if

  fu_nbr_of_verticals = GribContent%NbrVerts

end function fu_nbr_of_verticals


!*********************************************************************

integer function fu_nbr_of_fields (GribContent)
!
! Information only. Tells the number of fields
!
  implicit none
  
  ! Imported parameter
  type(Tinput_content), intent(in) :: GribContent

  if(.not.(GribContent%defined == silja_true)) then
    call set_error('Sorry, empty content','fu_nbr_of_fields')
    fu_nbr_of_fields = -1
    return
  end if

  fu_nbr_of_fields = GribContent%NbrFields

end function fu_nbr_of_fields 

!**********************************************************************

subroutine id_lst(IC, lst)

  implicit none

  type(Tinput_content), intent(in) ::IC
  type(silja_field_id),dimension(:),pointer :: lst

  if(.not.(IC%defined == silja_true)) then
    call set_error('Sorry, empty content','id_lst')
    return
  end if

  lst => IC%ids

end subroutine id_lst


!**********************************************************************

subroutine set_met_srcs(IC, wdr)
  !
  ! Stores the list of the meteodata sources into the wdr structure
  !
  implicit none

  ! Imported parameters
  !
  type(Tinput_content), intent(inout) :: IC
  type(silja_wdr), intent(inout) :: wdr
  !
  ! Local declarations
  integer :: i
  
  ! Some stupidity check first
  !
  if(.not.(IC%defined == silja_true)) then
    call set_error('Sorry, empty content','set_met_srcs')
    return
  end if
  ! if(.not.defined(wdr))then
  !  call set_error('Undefined wdr given','set_met_srcs')
  !  return
  !end if

  if(.not.IC%analysed)then
    call set_error('The content is not yet analysed','set_met_srcs')
    return
  end if

  ! ATTENTION. A dangerous but necessary trick: if in wdr the field 
  ! ifDisregardMDS == .true., all MDSs from input_content are ignored and wdr
  ! gets only one MDS: met_src_missing.
  ! The reason is: ECMWF changes models from time to time and each time sets 
  ! new MDS (actually, new model number) in the GRIB files. Unless the MDS is ignored,
  ! SILAM is unable to pass through the moment of new MDS appearing first time.
  ! A DANGER: if there were several MDSs from the very beginning with 
  ! overlapping content, it will be impossible to select the proper one
  ! unless their data can be distibguished somehow else - grid, vertical, etc.
  !
  if(fu_ifDisregardMDS(wdr))then
    call add_mds(wdr, met_src_missing)
  else
    do i=1, IC%NbrMetSrcs
      call add_mds(wdr, IC%MetSrcs(i)%src)
    end do
  endif

end subroutine set_met_srcs


!**********************************************************************

subroutine analyse_input_content(IC, out_grid, input_list, wdr)
  !
  ! Analyses the GribContent and gives an advise what to use as a system_grid
  ! vertical level type and meteo source. Checks the meteo data sources and sorts
  ! them so that the first source contains the most critical fields like wind
  ! temperature and pressure.
  !
  ! Idea:
  ! As long as we have a full list of input fields - we can decide the most
  ! convenient types of grid and verticals to be used. This should help to avoid 
  ! or reduce interpolations. Specifically, it will draw the maximum possible area
  ! size (of any type !), which is covered by all fields.
  !
  ! Algorithm in brief:
  ! 1. Priors: 
  !      - all model input quantities
  !      - GRIB content structure with all field id-s stored
  !  FILTERING:
  ! 2. Compile the lists of grids covering the target area. Eliminate all fields
  !    not covering it or covering insufficiently
  ! 3. Sort the grids in order of coverage and try to eliminate the bottom
  !    of the list. Criteria - satisfaction of the model needs. Compact the ids.
  ! 4. Compile the list of data sources together with their grids 
  ! 5. Try to find a list of data source solely satisfying the model needs. If not
  !    any source is rich enough - skip the next point
  ! 6. If one or more satisfying data source is found - take the best-resolution source
  !    and remove the others. Compact ids.
  !  OPTIMISATION / SHOPPING LIST FORMATION.
  ! 7. Count the variables. For each model input quanity from the QHierarchy - find 
  !    the shopping_list of initial quantities and then select the
  !    best-vertical-coverage/3d-resolution variable from the variable list.
  ! 8. Count the grids in the shopping list, sort them with regard to resolution. 
  !    Make the system_grid as the smallest-coverage and best-resolution one.
  ! 9. Count the verticals from the shopping list. Make the meteo_vertical as the
  !    largest-coverage, best-resolution one
  !
  implicit none

  ! Imported parameter
  type(Tinput_content), intent(inout) :: IC
  type(silja_grid), intent(in) :: out_grid
  type(silja_shopping_list), intent(inout) :: input_list
  type(silja_wdr), intent(in) :: wdr

  ! Local variables
  integer :: i, j, k, l, nx, ny, grid_type, iSrc, iBadId, n_cons, &
                & varPtr, n_grids_ini, iQtmp
  real :: best_resolution, fTmp
  logical :: found, OK
  real, dimension(max_quantities) :: gridCellSz  !Should be max_grids
  integer, dimension(4) :: flag_to_search = (/temperature_flag, pressure_flag, &
                                            & geopotential_flag, ground_pressure_flag/)
  logical, dimension(max_quantities) :: dataUsed, QCovered
  integer, dimension(max_quantities):: q_avail, & ! Available quantities from the file
                        &q_avail_st   !To-be the same for static files, Not yet..
  integer, dimension(max_quantities) :: q_shop,  q_shop_static
  type(silja_level) :: levelTmp
  type(silam_vertical) :: vertTmp


  !-----------------------------------------------------------------
  !
  ! 1. If GRIB content is defined - let's hope that all priors are OK
  !
  !-----------------------------------------------------------------
  !
  if(.not.IC%defined == silja_true)then
    call set_error('Undefined content given','analyse_input_content')
    return
  end if

  call set_missing(IC%shopping_list)
  q_shop=int_missing
  q_shop_static=int_missing



  !------------------------------------------------------------------
  !
  ! 2. Count the number of grids. Eliminate those covering less than 0.1
  !    of the target area
  !
  !------------------------------------------------------------------
  !
  call compact_ids(IC) ! Just removes the undefined ids - a precaution
  if(error)return

  if(test_messages)then
    call msg('Existing variables in the file (any grid and vertical)')
    q_avail(:) = int_missing
    do i=1, IC%NbrFields
      k=fu_merge_integer_to_array(fu_quantity(IC%ids(i)), q_avail)
    end do
    do i=1, size(q_avail)
      if(q_avail(i) == int_missing) exit
      call msg(fu_quantity_string(q_avail(i)), q_avail(i))
    end do
  end if


  call count_grids(IC) ! Make a list of grids
  if(error)return

  if(test_messages)then
    call msg('All grids from meteofile:')
    do i=1, IC%NbrGrids
      call report(IC%grids(i)%grid)
    end do
  endif

  call check_coverage(IC, coverage_threshold, out_grid) ! Removes insufficiently covering grids
  if(error)return
  !
  ! Is there any single grid left with above-threshold coverage ?
  !
  found = .false.
  do i=1, IC%NbrGrids
    if(defined(IC%grids(i)%grid)) then
      found = .true.
      exit
    endif
  end do
  if(.not. found)then
    call msg('')
    call msg('Requested output area:')
    call report(out_grid)
    call set_error('No grid covers the requested area','analyse_input_content')
    return
  endif


  !-------------------------------------------------------------------
  !
  ! 3a. Sort the grids in accordance with the coverage of the target area.
!!!!!!!!!!!     If the coverage is the same - in accordance with the number of fields
  !
  !-------------------------------------------------------------------
  !
  found = .true.
  do while (found)
    found = .false.
    do i=1, IC%NbrGrids-1
      if(IC%grids(i)%out_cover < IC%grids(i+1)%out_cover) then
        found  = .true.
        call exchange_grd_lsts(IC%grids(i), IC%grids(i+1))

!      elseif(IC%grids(i)%out_cover .eps. IC%grids(i+1)%out_cover)then
!        if(IC%grids(i)%NFlds_grd < IC%grids(i+1)%NFlds_grd) then
!          found  = .true.
!          call exchange_grd_lsts(IC%grids(i), IC%grids(i+1))
!        end if

      end if  ! coverage(i) <>== coverage(i+1)
    end do  ! cycle through the grids
  end do  ! while sorting 

  !--------------------------------------------------------------------
  !
  ! 3b. Eliminate unnesessary grids with low coverage. We want to keep only 
  !     grids with maximum coverage of the target area plus grids of 
  !     unique variables required for the model.
  !
  !--------------------------------------------------------------------
  !
  OK = .false.  ! All variables are covered by the grids
  do while(.not.OK)
    q_avail(:) = int_missing
    q_avail_st(:) = int_missing
    found = .false.
    !
    ! Count all grids with maximum coverage
    !
    do i=1, IC%NbrGrids
      if(IC%grids(i)%out_cover <= 0.9* IC%grids(1)%out_cover)exit
      do j=1, IC%grids(i)%NFlds_grd
!        print *, fu_quantity(IC%ids(IC%grids(i)%flds(j))), j
        k=fu_merge_integer_to_array(fu_quantity(IC%ids(IC%grids(i)%flds(j))), q_avail)
      end do
    end do
    i=i-1
    !
    ! Add grids until the list of variables is rich enough. However, it may happen that 
    ! some variable does not exist in the file at all or exists only in already rejected 
    ! grids. Then, either set error, or exclude variable, if it is not mandatory for
    ! the output.
    !
    mdl_q: do j=1, size(fu_quantities(input_list))
      if(fu_request(input_list,j) < 1) cycle ! There may be holes

      found = .false.
      do while(.not.found)
call msg('checking:' + fu_quantity_short_string(fu_quantity(input_list,j)))

          iQtmp = fu_quantity(input_list,j)
          call check_input_quantity(iQtmp,& 
                                  & fu_realtime_quantity(iQtmp), &
                                & q_avail, q_avail_st, q_shop, q_shop_static, &
                                  & found, wdr)
call msg('Done:' + fu_quantity_short_string(fu_quantity(input_list,j)))


        if(.not.found)then
          i=i+1       ! Take the next grid
          if(i>IC%NbrGrids)then ! If the list of grids is expired
            if(fu_request(input_list,j) == 2)then ! mandatory
              call set_error('Not all quantities found in data files','analyse_input_content')
              call msg(fu_connect_strings(' *** FAILED *** ', &
                             & fu_quantity_string(fu_quantity(input_list,j))))
              return
            else
              call set_request(input_list,j, 0)  ! Exclude this variable from the list
              call msg_warning(fu_connect_strings('Variable excluded from output:', &
                 &   fu_quantity_string(fu_quantity(input_list,j))),'analyse_input_content')
              exit mdl_q
            endif
          end if
          do l=1, IC%grids(i)%NFlds_grd
            k=fu_merge_integer_to_array(fu_quantity(IC%ids(IC%grids(i)%flds(l))), q_avail)
          end do
        endif ! Variable found in existing grids
      end do  ! while .not.found

    end do mdl_q ! cycle through QHierarchy
    OK = found
  end do  ! While all variables are covered

  !---------------------------------------------------------------
  !
  ! 3c. Remove unused grids, compact the list of ids
  !
  !---------------------------------------------------------------
  !
  if(i < IC%NbrGrids)then ! Not all grids were used - remove the loosers
    do j=i+1,IC%NbrGrids
      call remove_grid(IC, j)
    end do
    call compact_ids(IC)
    call count_grids(IC)
  end if
  if(error)return

  !---------------------------------------------------------------
  !
  ! 4. Compile the list of data sources together with their grids 
  !    Sort them in accordance with the resolution of the "best" grid
  !
  !---------------------------------------------------------------

  
  call count_data_sources(IC)

  !
  !  Sorting the data sources according to the resolution
  !
  found = .true.
  do while (found)
    found = .false.
    do i=1, IC%NbrMetSrcs-1
      if(IC%MetSrcs(i)%best_cell_size > IC%MetSrcs(i+1)%best_cell_size) then
        found  = .true.
        call exchange_met_src_lsts(IC%MetSrcs(i), IC%MetSrcs(i+1))
      end if
    end do
  end do  ! while sorting 

!  !
!  ! Sorting the data sources in accordance with number of fields using each source
!  !
!  found = .true.
!  do while (found)
!    found = .false.
!    do i=1, IC%NbrMetSrcs-1
!      if(IC%MetSrcs(i)%NFlds_src < IC%MetSrcs(i+1)%NFlds_src) then
!        found  = .true.
!        call exchange_met_src_lsts(IC%MetSrcs(i), IC%MetSrcs(i+1))
!      end if
!    end do
!  end do  ! while sorting 



  !-----------------------------------------------------------------
  !  
  ! 5. Try to find a list of data source solely satisfying the model needs.
  !
  !-----------------------------------------------------------------

  if(IC%NbrMetSrcs > 1)then
    dataUsed = .false.

    do i=1,IC%NbrMetSrcs
      QCovered = .false.
      do j=1,IC%MetSrcs(i)%NFlds_src
        k=fu_merge_integer_to_array(fu_quantity(IC%ids(IC%MetSrcs(i)%flds(j))), q_avail)
      end do

      do j = 1, size(fu_quantities(input_list))
        if(fu_request(input_list,j) < 1) then
          QCovered(i) = .true.  ! force it !
        else
          iQtmp = fu_quantity(input_list,j)
          call check_input_quantity(iQtmp,& 
                & fu_realtime_quantity(iQtmp), &
                        & q_avail, q_avail_st, q_shop, q_shop_static, &
                        & found, wdr)
          if(found) then  ! Model quantity is available or can be derived
            QCovered(i) = .true.
          end if
        endif
      end do  ! size(QHierarchy)

      if(all(QCovered(1:size(fu_quantities(input_list)))))then
        dataUsed(i) = .true.
      end if
    end do

    !--------------------------------------------------------------------
    !
    ! 6. If one or more solely sufficient sources exist - remove the others. Also, 
    !    thanks to the sorting order - we take the first accepted data source and
    !    eliminate the others. Resolution will be the best.
    !    After done - compact ids
    !
    !--------------------------------------------------------------------
    !
    if(any(dataUsed))then
      do i=1, size(dataUsed)
        if(dataUsed(i)) then
          !
          ! The source is self-sufficient. Put him to the number 1 and remove the rest
          !
          if( i /= 1) call exchange_met_src_lsts(IC%MetSrcs(i), IC%MetSrcs(1))
          do j=2, IC%NbrMetSrcs
            call remove_met_src(IC, j)
          end do
          exit
        end if
      end do
      call compact_ids(IC)
    end if                  ! Solely-sufficient data source found 

  end if ! if IC%NbrMetSrcs > 1

  !-------------------------------------------------------------------------
  !
  ! 7a. Count variables. A variable contains references to the corresponding
  !     grid, source and vertical. So, they all should be counted first
  !
  !-------------------------------------------------------------------------

  call compact_ids(IC) ! Just precaution
  call count_grids(IC)
  call check_coverage(IC, coverage_threshold, out_grid) ! Compute the coverage
  call count_data_sources(IC) ! Just precaution
  call count_verticals(IC)
  if(error)return

  call count_variables(IC)
  if(error) return

  dataUsed = .false.

  !-------------------------------------------------------------------------
  !
  ! 7b. For each model input quanity from the QHierarchy - find 
  !    the shopping_list of initial quantities and then select the
  !    best-vertical-coverage/3d-resolution variable from the variable list.
  !
  !-------------------------------------------------------------------------
  !
  ! Fill-in the array of what is available now
  !
  q_avail = int_missing
  do i=1, IC%NbrVars
    k=fu_merge_integer_to_array(IC%vars(i)%quantity, q_avail)
  end do
  !
  ! Actual checking starts
  !
  do i=1, size(fu_quantities(input_list))
    if(fu_request(input_list,i) < 1)cycle  ! If the quantity was excluded - just skip it
    iQtmp = fu_quantity(input_list,i)
    call check_input_quantity(iQtmp,& 
                & fu_realtime_quantity(iQtmp), &
                & q_avail, q_avail_st, q_shop, q_shop_static, &
                            & found, wdr)
    if(.not.found)then
      call set_error('Strange, missing model input quantity','analyse_input_content')
      return
    end if
    !
    ! Run through the q_shop, searching for each element the best available from the
    ! list of variables
    !
!    call report(IC%shopping_list)

    do j=1, size(q_shop)
      if(q_shop(j) == int_missing) exit
      !
      ! Find all such quantities in list of vars and compute quality_value for each
      !
      varPtr = 0
      do k = 1, IC%NbrVars
        if(IC%vars(k)%quantity == q_shop(j)) then
          if(varPtr == 0) varPtr = k
          varPtr = fu_compare_vars_quality(IC, k, varPtr)
        end if
      end do ! cycle through vars
      if(varPtr == 0)then
        call set_error('Strange, zero varPtr','analyse_input_content')
        return
      end if

      vertTmp = IC%verts(IC%vars(varPtr)%vertPtr)%vert
      if(fu_nbrOfLevels(vertTmp) > 1) call arrange_levels_in_vertical(vertTmp)
      call add_shopping_variable(IC%shopping_list, &
                               & IC%vars(varPtr)%quantity, species_missing, & ! quaintity, species
                               & IC%grids(IC%vars(varPtr)%gridPtr)%grid, &
!                               & IC%verts(IC%vars(varPtr)%vertPtr)%vert, &
                               & vertTmp, &                                      !IC%verts(IC%vars(varPtr)%vertPtr)%vert, &
                               & IC%vars(varPtr)%vertLevNbr, &
                               & IC%MetSrcs(IC%vars(varPtr)%MetSrcPtr)%src)
      if(error)return


!      call report(IC%shopping_list)


      dataUsed(varPtr) = .true.

    end do ! cycle through the q_shop

  end do ! Cycle through the QHIerarchy

  !--------------------------------------------------------------------------
  !
  ! 7c. Last cleaning. Remove all ids, which were not included into the 
  !     shopping-list. Then - recount everything
  !
  !-------------------------------------------------------------------------
  !
  do i = 1, IC%NbrVars
    if(.not.dataUsed(i)) call remove_variable(IC, i)
  end do

  call compact_ids(IC)
  call count_grids(IC)
  call check_coverage(IC, coverage_threshold, out_grid) ! Compute the coverage
  call count_data_sources(IC) ! Just precaution
  call count_verticals(IC)
  if(error)return

!  if(test_messages)then
!    call msg('All grids before making the system_grid from a common area:')
!    do i=1, IC%NbrGrids
!      call report(IC%grids(i)%grid)
!    end do
!  endif

  !----------------------------------------------------------------------------
  !
  ! 8. Make the meteo_grid as the smallest-coverage and best-resolution one.
  !
  !----------------------------------------------------------------------------
  !
  if(IC%NbrGrids < 1) then ! Error somewhere above. Data do not cover the target area 
    call set_error('Strange, no grids for system_grid selection','analyse_input_content')
    return
  end if

  ! Now, let's select the meteo_grid.
  ! If there is a grid serving the bulk of the fields - let say, 10 times the others, just use it.
  ! In case of no clear preference, find the best resolution.
  !
  if(IC%NbrGrids > 1) then
    !
    ! First, the most used grid
    !
    j = 0
    do i=1,IC%NbrGrids
      if(.not.defined(IC%grids(i)%grid))cycle
      if(IC%grids(i)%NFlds_grd > j)then
        j = IC%grids(i)%NFlds_grd
        IC%SystemGridPtr = i ! Not finally, however...
      end if
    end do  
    !
    ! ... and check that it is by far the most-used one
    !
    found = .true.
    do i=1,IC%NbrGrids
      if(.not.defined(IC%grids(i)%grid) .or. i == IC%SystemGridPtr)cycle
      if(IC%grids(i)%NFlds_grd * 10 > IC%grids(IC%SystemGridPtr)%NFlds_grd)then
        found = .false.
        exit
      end if
    end do  
    
    if(.not. found)then
      !
      ! No clear mostly-used grid. Will select it based on resolution
      !
      best_resolution = 1.0e+35 ! [m**2] Quite poor resolution indeed

      do i=1,IC%NbrGrids
        if(.not.defined(IC%grids(i)%grid))cycle
        gridCellSz(i) =  fu_cell_size(IC%grids(i)%grid)
        if(gridCellSz(i) < best_resolution)then ! In meters**2 !!!
          best_resolution = gridCellSz(i)         
          IC%SystemGridPtr = i ! Not finally, however...
        end if
      end do  
      !
      ! If there are several grids with the best resolution (e.g. Arakawa grids)
      ! - the grid of temperature, pressure or geopotential wins. However, if
      ! the best-reslution does not correspond to one of the above basic fields
      ! the selection of the grid becomes arbitrary - the first-met best-resolution grid
      ! will be taken
      !
      found = .false.
      do i=1,IC%NbrGrids
        if(best_resolution < 0.99*gridCellSz(i)) cycle ! There are better...
        if(i == IC%SystemGridPtr) cycle        ! This grid is already claimed as winner
        found = .true.  ! We got another grid, which has the same resolution
        exit
      end do

      if(found)then
        !
        ! Determine what flag we shall search
        !
        found=.false.
        do j=1,size(flag_to_search)
          do i=1,IC%NbrGrids
            do k=1,IC%grids(i)%NFlds_grd
              if(fu_quantity(IC%ids(IC%grids(i)%flds(k))) == flag_to_search(j))then
                flag_to_search(1) = flag_to_search(j)
                found = .true.
                exit
              end if
            end do  ! fields for grid
          end do ! grids
          if(found)exit
        end do  ! cycle through flag_to_search
        if(.not.found)then
          call set_error('No basic quantities found for selected grids','analyse_input_content')
          return
        end if
        !
        ! Now scan all best-resolution grids looking for the flag_to_search
        !
        found = .false.
        do i=1,IC%NbrGrids
          if(best_resolution < 0.99*gridCellSz(i)) cycle ! There are better...
          if(i == IC%SystemGridPtr) cycle        ! This grid is already claimed as winner

          do j=1,IC%grids(i)%NFlds_grd
            if(fu_quantity(IC%ids(IC%grids(i)%flds(j))) == flag_to_search(1))then
              found = .true.
              exit
            end if
          end do
          if(found)then ! This very grid has the best resolution plus flag_to_search
            IC%SystemGridPtr = i  ! The winning grid is changed
            exit
          end if
        end do

      end if ! Several grids with the same resolution
    endif ! if there is a clearly most-used grid


!! FIXME I could not imagine why this check might be needed
!! It was killingly slow on huge (2k x 2k = 4M) indian WRF grids
    !
    ! OK, the meteo_grid is determined, now - the next question: 
    ! What is the largest area covered by all grids ? Just check-and-cut
    ! the meteo_grid if its original size is bigger than the others
    !
!!!!      do i=1,IC%NbrGrids
!!!!  !      call msg('Grid to cut:')
!!!!  !      call report(IC%grids(i)%grid)
!!!!        call cut_grid_size(IC%grids(IC%SystemGridPtr)%grid, IC%grids(i)%grid, inside_the_grid_area)
!!!!  !      call msg('Grid after cut:')
!!!!  !      call report(IC%grids(i)%grid)
!!!!        if(error)then
!!!!          call msg('METEO_GRID (under cutting):')
!!!!          call report(IC%grids(IC%SystemGridPtr)%grid)
!!!!          call msg('')
!!!!          call msg('Extra variables are in this grid:')
!!!!          call report(IC%grids(i)%grid)
!!!!          call msg('It looks like the file contains non-overlapping grids, which')
!!!!          call msg('are both needed for the calculations. No way to continue')
!!!!          return
!!!!        end if
!!!!      end do

  else 
    IC%SystemGridPtr = 1 ! Just one grid exists

  end if ! If more than one grid exists

  if(.not.fu_stdSilamGrid(IC%grids(IC%SystemGridPtr)%grid))then
    call set_error('The meteo grid is not SILAM-standard and cannot be made to comply', &
                 & 'analyse_input_content')
    call report(IC%grids(IC%SystemGridPtr)%grid)
    return
  endif


  !----------------------------------------------------------------------------
  !
  ! 9. Make the meteo_vertical as the largest-coverage, best-resolution one
  !
  !----------------------------------------------------------------------------

  call msg('All verticals in the meteo file')
  do i=1,IC%NbrVerts
    call report(IC%verts(i)%vert, .true.)
  end do

  if(IC%NbrVerts > 1)then ! Analysis continues only if >1 vertical types found
    
    ! What field can be used as a reference one for the vertical strucutre ?
    ! The first choice is evidently wind as the most sensitive to vertical
    ! interpolations. Find it.
    !
    found = .false.
    do i=1,IC%NbrFields
      if(fu_quantity(IC%ids(i)) == u_flag .or. &
       & fu_quantity(IC%ids(i)) == v_flag) then
        
        do j=1,IC%NbrVerts
          if( fu_leveltype(IC%verts(j)%vert) == fu_leveltype(fu_level(IC%ids(i))) .and. &
             & any(  fu_leveltype(fu_level(IC%ids(i))) == (/hybrid, sigma_level/))    )then
            ! Meteo vertical must be hybrid... otherwise we are in trouble
            found = .true.
            IC%SystemLevPtr = j
            call msg("Selecting vertical")
            call report(IC%verts(j)%vert)
            exit
          end if
        end do
      end if ! If wind field
      if (found) exit
    end do ! Cycle through the fields searching for the wind field
    !
    ! If no wind in the content, it is either an error or a pure-meteo processing without wind.
    ! Meteo vertical then will be the one with max number of levels
    !
    if(.not. found)then
      call msg_warning('No wind fields are found in the input content','analyse_input_content')
      IC%SystemLevPtr = 1
      do j=1,IC%NbrVerts
        if(fu_NbrOfLevels(IC%verts(j)%vert) > fu_NbrOfLevels(IC%verts(IC%SystemLevPtr)%vert)) &
                                                            & IC%SystemLevPtr = j
      end do
    endif   ! not found  wind
    
  else ! Just one vertical type found => nothing to do
    IC%SystemLevPtr = 1

  end if ! If >1 vertical types found

  !
  ! Having the vertical selected, it has to be brought to the SILAM standards
  !
  if(IC%SystemLevPtr > 0)call arrange_levels_in_vertical(IC%verts(IC%SystemLevPtr)%vert)
  
  IC%analysed = .true. ! Everything is OK.


end subroutine analyse_input_content


!**********************************************************************
!**********************************************************************
!**********************************************************************
!
!            PRIVATE stuff
!
!**********************************************************************
!**********************************************************************


!**********************************************************************

subroutine compact_ids(IC)
  !
  ! Removes the missing ids from the list. Also nullifies all lists - grids,
  ! data sources, verticals, variables
  !
  implicit none

  ! Imported variables 
  type(Tinput_content), intent(inout) :: IC

  ! Local variables
  integer :: i, j

  i=1
  do while (i <= IC%NbrFields)

    if (.not.(defined(IC%ids(i)) .and. &
            & defined(fu_grid(IC%ids(i))) .and. & 
            & defined(fu_level(IC%ids(i))) .and. &
            & fu_known_quantity(fu_quantity(IC%ids(i))))) then
      do j = i, IC%NbrFields-1
        IC%ids(j) = IC%ids(j+1)
      end do
      do j = IC%NbrFields, size(IC%ids)
        call set_missing(IC%ids(j))
      end do
      IC%NbrFields = IC%NbrFields - 1
      i=i-1 ! May be, the i+1 field was also empty
    end if
    i=i+1
  end do ! Through fields

  do i=1, IC%NbrGrids
    IC%grids(i)%grid = grid_missing
    IC%grids(i)%NFlds_grd = int_missing
    IC%grids(i)%flds(:) = int_missing
    IC%grids(i)%out_cover = 0.
  end do
  IC%NbrGrids = 0

  do i=1, IC%NbrMetSrcs
    IC%MetSrcs(i)%src = met_src_missing
    IC%MetSrcs(i)%NFlds_src = int_missing
    IC%MetSrcs(i)%flds(:) = int_missing
    IC%MetSrcs(i)%best_cell_size = real_missing
  end do
  IC%NbrMetSrcs = 0

  do i=1, IC%NbrMetSrcs
    call set_missing(IC%verts(i)%vert, .false.)
    IC%verts(i)%NFlds_vert = int_missing
    IC%verts(i)%flds(:) = int_missing
  end do
  IC%NbrVerts = 0

  do i=1, IC%NbrVars
    IC%vars(i)%quantity = int_missing
    IC%vars(i)%NFlds_var = int_missing
    IC%vars(i)%flds(:) = int_missing
  end do
  IC%NbrVars = 0

end subroutine compact_ids


!************************************************************************

subroutine count_grids(IC)
  !
  ! Counts the number of grids in the GRIB content structure. Also
  ! fills-in the % coverage of the output grid
  !
  implicit none

  ! Imported variables 
  type(Tinput_content), intent(inout) :: IC

  ! Local variables
  integer :: i, j
  logical :: found
!  type(silja_grid) :: gridTmp

  do j = 1, IC%NbrFields
    !
    ! Some inconsistency check 
    !
    if(.not.defined(fu_grid(IC%ids(j)))) then
      call set_error('Undefined grid found','count_grids')
      return
    end if
    !
    ! Field is OK - analyse the grid
    !
    found = .false.
    do i=1,IC%NbrGrids    ! Is it already stored ?
      if(IC%grids(i)%grid == fu_grid(IC%ids(j))) then
!      if(IC%grids(i)%grid == gridTmp) then
        found = .true.
        exit
      end if
    end do
    if(found) then ! Store the index of the field with this grid
      IC%grids(i)%NFlds_grd = IC%grids(i)%NFlds_grd + 1
      IC%grids(i)%flds(IC%grids(i)%NFlds_grd) = j
    else
      IC%NbrGrids = IC%NbrGrids + 1  ! Store the new grid

      call check_gridlist_size(IC, IC%NbrGrids)
      if(error)return

!      IC%grids(IC%NbrGrids)%grid = gridTmp
      IC%grids(IC%NbrGrids)%grid = fu_grid(IC%ids(j))
      IC%grids(IC%NbrGrids)%NFlds_grd = 1
      IC%grids(IC%NbrGrids)%flds(1) = j
      IC%grids(IC%NbrGrids)%out_cover = 0.
      !
      ! If the grid is non-standard but global, it should be rotated
      !

    end if
  end do ! Through fields

  ! Clean the rest of the list
  do j= IC%NbrGrids+1, size(IC%grids)
    IC%grids(j)%grid = grid_missing
    IC%grids(j)%NFlds_grd = 0
    IC%grids(j)%flds(:) = int_missing
    IC%grids(j)%out_cover = 0.
  end do

end subroutine count_grids


!**********************************************************************

subroutine check_coverage(IC, threshold, out_grid)
  !
  ! Cheks the fraction of the output grid covered by the grids in IC.
  ! Fields covering less than threshold are removed from the list of ods
  !
  implicit none

  ! Imported variables 
  type(Tinput_content), intent(inout) :: IC
  type(silja_grid), intent(in) :: out_grid
  real, intent(in) :: threshold

  ! Local variables
  integer :: i, j
  real :: cover
  type(silja_grid) :: out_stag

  integer :: nx_gtm, ny_gtm, iy
  integer, dimension(:), pointer :: iWork
  integer, dimension(:,:), pointer :: covermask

  if(defined(out_grid))then
    i=1
    
    out_stag = fu_staggered_grid("out_stag",out_grid, .True., .True.)

    do while (i <= IC%NbrGrids)
!      cover = fu_area_coverage(out_grid, IC%grids(i)%grid)
!     if(cover < threshold)then
!       call remove_grid(IC, i) 
!     else
!       IC%grids(i)%out_cover = cover
!     end if

#ifdef DEBUG
    call grid_dimensions(IC%grids(i)%grid, nx_gtm, ny_gtm)
    iWork => fu_work_int_array(nx_gtm*ny_gtm)
    covermask(1:nx_gtm,1:ny_gtm) =>  iWork(1:nx_gtm*ny_gtm)

    if( .not. fu_if_mesh_intrpolatable(out_stag,  IC%grids(i)%grid, covermask, ifSpk = .True.) )then
       call msg("Grid is NOT covered! Cover mask:")
       do iy = ny_gtm , 1, -1 !Top-bottom, so map looks normal
         call msg("",covermask(:,iy) > 0)
       enddo

#else
    if( .not. fu_if_mesh_intrpolatable(out_stag,  IC%grids(i)%grid))then
#endif

        call remove_grid(IC, i) 
      else
        IC%grids(i)%out_cover = 1.
      end if
#ifdef DEBUG
    call free_work_array(iWork)
#endif
      i=i+1
    end do ! Through fields
  else
    do i = 1, IC%NbrGrids
      IC%grids(i)%out_cover = 1.
    end do
  endif
  
end subroutine check_coverage


!***************************************************************************

subroutine remove_grid (IC, indexGrd)
  ! 
  ! Removes the grid and all fields with this grid from the GRIB content 
  ! structure
  !
  implicit none

  ! Imported variables with intent IN
  integer, intent(in) :: indexGrd

  ! Imported parameters with intent INOUT
  type(Tinput_content), intent(inout) :: IC

  ! Local variable
  integer :: iTmp

  if(.not. IC%defined == silja_true) then
    call set_error('Undefined content given','remove_grid')
    return
  end if

  if(indexGrd < 1 .or. indexGrd > IC%NbrGrids) then
    call msg_warning('Strange index given. Skipped','remove_grid')
    return
  end if

  ! Kill all ids with this grid, but not move ids - otherwise you
  ! loose the connection with grids
  !
  do iTMp = 1, IC%grids(indexGrd)%NFlds_grd
    call set_missing(IC%ids(IC%grids(indexGrd)%flds(iTmp)))
  end do

  IC%grids(indexGrd)%flds(:) = int_missing
  IC%grids(indexGrd)%grid = grid_missing
  IC%grids(indexGrd)%NFlds_grd = 0

end subroutine remove_grid


!**********************************************************************

subroutine exchange_grd_lsts(grd_lst1, grd_lst2) 
  !
  ! Exchanges two Tgrid_lst objects
  !
  implicit none
  !
  ! Imported parameters with intent INOUT
  type(Tgrid_lst), intent(inout) :: grd_lst1, grd_lst2

  ! Local variables
  type(Tgrid_lst) :: glTmp

  glTmp%grid = grd_lst1%grid
  grd_lst1%grid = grd_lst2%grid
  grd_lst2%grid = glTmp%grid

  glTmp%NFlds_grd = grd_lst1%NFlds_grd
  grd_lst1%NFlds_grd = grd_lst2%NFlds_grd
  grd_lst2%NFlds_grd = glTmp%NFlds_grd

  glTmp%out_cover = grd_lst1%out_cover
  grd_lst1%out_cover = grd_lst2%out_cover
  grd_lst2%out_cover = glTmp%out_cover

  glTmp%flds(1:size(glTmp%flds)) = grd_lst1%flds(1:size(glTmp%flds))
  grd_lst1%flds(1:size(glTmp%flds)) = grd_lst2%flds(1:size(glTmp%flds))
  grd_lst2%flds(1:size(glTmp%flds)) = glTmp%flds(1:size(glTmp%flds))

end subroutine exchange_grd_lsts


!*******************************************************************

subroutine check_gridlist_size(IC, iSize)
  !
  ! Checks the size of the local list of grids and enlarges it if needed
  !
  implicit none

  type(Tinput_content), intent(inout) :: IC
  integer, intent(in) :: iSize

  integer :: kret, i
  type(Tgrid_lst),dimension(:), allocatable :: lstTmp

  if(.not.associated(IC%grids))then
    allocate(IC%grids(10*(int(real(iSize)/10.)+1)),stat=kret)
    if(kret /= 0)then
      call set_error('No more space available','analyse_grib_content')
      return
    end if
    return
  end if

  if(iSize > size(IC%grids))then
    !
    ! Create temporary space and store the content there
    !
    allocate(lstTmp(size(IC%grids)),stat=kret)
    if(kret /= 0)then
      call set_error('No more space available','analyse_grib_content')
      return
    end if
    do i=1,size(IC%grids)
      lstTmp(i)%grid=IC%grids(i)%grid
      lstTmp(i)%NFlds_grd = IC%grids(i)%NFlds_grd
      lstTmp(i)%out_cover = IC%grids(i)%out_cover
      lstTmp(i)%cell_size = IC%grids(i)%cell_size
      lstTmp(i)%flds(:)=IC%grids(i)%flds(1:size(IC%grids(i)%flds))
    end do

    !
    ! Resize the grid list and return the content back
    !
    deallocate(IC%grids)
    allocate(IC%grids(10*(int(real(iSize)/10.)+1)),stat=kret)
    if(kret /= 0)then
      call set_error('No more space available','analyse_grib_content')
      return
    end if

    IC%grids(:)%grid = grid_missing
    do i=1,size(lstTmp)
      IC%grids(i)%grid=lstTmp(i)%grid
      IC%grids(i)%NFlds_grd = lstTmp(i)%NFlds_grd
      IC%grids(i)%out_cover = lstTmp(i)%out_cover
      IC%grids(i)%cell_size = lstTmp(i)%cell_size
      IC%grids(i)%flds(1:size(IC%grids(i)%flds)) = lstTmp(i)%flds(:)
    end do

    deallocate(lstTmp)

  end if
end subroutine check_gridlist_size



!*******************************************************************

integer function fu_gridPtr(IC, grid)
  !
  ! Finds the index of the given grid in the list of grids. Returns 
  ! int_missing if not found
  !
  implicit none

  ! Imported parameters with intent IN
  type(Tinput_content), intent(in) :: IC
  type(silja_grid), intent(in) :: grid

  ! Local variables
  integer :: i
  !
  ! Some stupidity check first
  !
  fu_gridPtr = int_missing
  if(IC%NbrGrids < 1)then
    call set_error('No grids in the list','fu_gridPtr')
    return
  end if

  do i=1, IC%NbrGrids
    if(IC%grids(i)%grid == grid)then
      fu_gridPtr = i
      return
    end if
  end do
  call set_error('Grid is not in the list','fu_gridPtr')

end function fu_gridPtr


!*******************************************************************

subroutine count_data_sources(IC)
  !
  ! Make a list of the data sources in the grib file
  !
  implicit none

  ! Imported parameters with intent INOUT
  type(Tinput_content), intent(inout) :: IC

  ! Local variables
  integer :: i, j
  logical :: found

  do j=1,IC%NbrFields
    found = .false.
    do i=1,IC%NbrMetSrcs    ! Is it already stored ?
      if(IC%MetSrcs(i)%src == fu_met_src(IC%ids(j))) then
        found = .true.
        exit
      end if
    end do
    if(found) then ! Store the index of the field with this met_src
      IC%MetSrcs(i)%NFlds_src = IC%MetSrcs(i)%NFlds_src + 1
      IC%MetSrcs(i)%flds(IC%MetSrcs(i)%NFlds_src) = j
      IC%MetSrcs(i)%best_cell_size = &
          & min(IC%MetSrcs(i)%best_cell_size, fu_cell_size(fu_grid(IC%ids(j))))
    else
      IC%NbrMetSrcs = IC%NbrMetSrcs + 1  ! Store the new source

      call check_MetSrcs_list_size(IC, IC%NbrMetSrcs)
      if(error)return

      IC%MetSrcs(IC%NbrMetSrcs)%src = fu_met_src(IC%ids(j))
      IC%MetSrcs(IC%NbrMetSrcs)%NFlds_src = 1
      IC%MetSrcs(IC%NbrMetSrcs)%flds(1) = j
      IC%MetSrcs(IC%NbrMetSrcs)%best_cell_size = fu_cell_size(fu_grid(IC%ids(j)))
    end if
  end do
  if(error)return

end subroutine count_data_sources


!**********************************************************************

subroutine exchange_met_src_lsts(src_lst1, src_lst2) 
  !
  ! Exchanges two TMet_src_lst objects
  !
  implicit none
  !
  ! Imported parameters with intent INOUT
  type(TMet_src_lst), intent(inout) :: src_lst1, src_lst2

  ! Local variables
  type(TMet_src_lst) :: mdsTmp

  mdsTmp%src = src_lst1%src
  src_lst1%src = src_lst2%src
  src_lst2%src = mdsTmp%src

  mdsTmp%NFlds_src = src_lst1%NFlds_src
  src_lst1%NFlds_src = src_lst2%NFlds_src
  src_lst2%NFlds_src = mdsTmp%NFlds_src

  mdsTmp%best_cell_size = src_lst1%best_cell_size
  src_lst1%best_cell_size = src_lst2%best_cell_size
  src_lst2%best_cell_size = mdsTmp%best_cell_size

  mdsTmp%flds(1:size(mdsTmp%flds)) = src_lst1%flds(:)
  src_lst1%flds(1:size(mdsTmp%flds)) = src_lst2%flds(:)
  src_lst2%flds(1:size(mdsTmp%flds)) = mdsTmp%flds(:)

end subroutine exchange_met_src_lsts


!***************************************************************************

subroutine remove_met_src (IC, indexMDS)
  ! 
  ! Removes the grid and all fields with this grid from the GRIB content 
  ! structure
  !
  implicit none

  ! Imported variables with intent IN
  integer, intent(in) :: indexMDS

  ! Imported parameters with intent INOUT
  type(Tinput_content), intent(inout) :: IC

  ! Local variables
  integer :: iTmp

  if(.not. IC%defined == silja_true) then
    call set_error('Undefined content given','remove_met_src')
    return
  end if

  if(indexMDS < 1 .or. indexMDS > IC%NbrMetSrcs) then
    call msg_warning('Strange index given. Skipped','remove_met_src')
    return
  end if

  ! Kill all ids with this data source, but do not move ids - otherwise you
  ! loose the connection with the other sources
  !
  do iTmp = 1, IC%MetSrcs(indexMDS)%NFlds_src
    call set_missing(IC%ids(IC%MetSrcs(indexMDS)%flds(iTmp)))
  end do

  IC%MetSrcs(indexMDS)%flds(:) = int_missing
  IC%MetSrcs(indexMDS)%src = met_src_missing
  IC%MetSrcs(indexMDS)%NFlds_src = int_missing

end subroutine remove_met_src


!******************************************************************************

subroutine check_MetSrcs_list_size(IC, iSize)
  !
  ! Checks the size of the local list of meteodata sources and enlarges it if needed
  !
  implicit none

  ! Imported parameters
  type(Tinput_content), intent(inout) :: IC
  integer, intent(in) :: iSize

  integer :: kret, i
  type(Tmet_src_lst), dimension(:), allocatable :: lstTmp

  if(.not.associated(IC%MetSrcs))then
    allocate(IC%MetSrcs(2*(int(real(iSize)/2.)+1)),stat=kret)
    if(kret /= 0)then
      call set_error('No more space available','analyse_grib_content')
      return
    end if
    return
  end if

  if(iSize > size(IC%MetSrcs))then
    !
    ! Create temporary space and store the content there
    !
    allocate(lstTmp(size(IC%MetSrcs)),stat=kret)
    if(kret /= 0)then
      call set_error('No more space available','analyse_grib_content')
      return
    end if
    do i=1,size(IC%MetSrcs)
      lstTmp(i)%src=IC%MetSrcs(i)%src
      lstTmp(i)%NFlds_src = IC%MetSrcs(i)%NFlds_src
      lstTmp(i)%best_cell_size = IC%MetSrcs(i)%best_cell_size
      lstTmp(i)%flds(:)=IC%MetSrcs(i)%flds(1:size(IC%MetSrcs(i)%flds))
    end do
    !
    ! Resize the main list and drop the content back
    !
    deallocate(IC%MetSrcs)
    allocate(IC%MetSrcs(2*(int(real(iSize)/2.)+1)),stat=kret)
    if(kret /= 0)then
      call set_error('No more space available','analyse_grib_content')
      return
    end if

    IC%MetSrcs(:)%src = met_src_missing
    do i=1,size(lstTmp)
      IC%MetSrcs(i)%src = lstTmp(i)%src
      IC%MetSrcs(i)%NFlds_src = lstTmp(i)%NFlds_src
      IC%MetSrcs(i)%best_cell_size = lstTmp(i)%best_cell_size
      IC%MetSrcs(i)%flds(1:size(IC%MetSrcs(i)%flds)) = lstTmp(i)%flds(:)
    end do

    deallocate (lstTmp)
  end if
end subroutine check_MetSrcs_list_size



!*******************************************************************

integer function fu_MetSrcPtr(IC, mds)
  !
  ! Finds the index of the given meteo data source in the list of sources.
  ! Returns int_missing if not found
  !
  implicit none

  ! Imported parameters with intent IN
  type(Tinput_content), intent(in) :: IC
  type(meteo_data_source), intent(in) :: mds

  ! Local variables
  integer :: i
  !
  ! Some stupidity check first
  !
  if(IC%NbrMetSrcs < 1)then
    call set_error('No data sources in the list','fu_MetSrcPtr')
    return
  end if

  do i=1, IC%NbrMetSrcs
    if(IC%MetSrcs(i)%src == mds)then
      fu_MetSrcPtr = i
      return
    end if
  end do
  call set_error('Data source is not in the list','fu_MetSrcPtr')

end function fu_MetSrcPtr



!********************************************************************

subroutine count_verticals(IC)
  !
  ! Counts the number of vertical structures, taking into account 
  ! a multi-level feature of each structure. Extra dimensions of the
  ! split are - sources and grids. So, the vertical structure belongs
  ! to a specific source and to a specific grid (contrary to grids and
  ! sources, by the way).
  !
  implicit none

  ! Imported parameters
  type(Tinput_content), intent(inout) :: IC

  ! Local variables
  integer :: i,j, k
  logical :: found

  ! Since vertical belongs to the source and may exist only over grid - they
  ! have to be defined before
  !
  if(IC%NbrGrids < 1 .or. IC%NbrMetSrcs < 1)then
    call set_error('Sources or grids are not counted yet','count_verticals')
    return
  end if

  do j=1,IC%NbrFields
    if(.not.defined(fu_level(IC%ids(j)))) then
      call set_error('Undefined level given','count_verticals')
      cycle
    end if

!    call msg('Checking the field:',j)
!    if(j== 167)then
!      call msg('167')
!    endif
    
    found = .false.
    !
    ! Is this vertical already in the list ? Check vertical type, source and grid
    !
    do i=1,IC%NbrVerts
      if(fu_leveltype(IC%verts(i)%vert) == fu_leveltype(fu_level(IC%ids(j)))) then
        if(IC%grids(IC%verts(i)%gridPtr)%grid == fu_grid(IC%ids(j)))then
          if(IC%MetSrcs(IC%verts(i)%srcPtr)%src == fu_met_src(IC%ids(j)))then
            found = .true.
            exit
          end if
        end if
      end if
    end do

    if(found) then 
      !
      ! Store the index of the field with this vertical type, grid and data source
      !
      IC%verts(i)%NFlds_vert = IC%verts(i)%NFlds_vert + 1
      IC%verts(i)%flds(IC%verts(i)%NFlds_vert) = j

      found = .false.
      do k=1,fu_NbrOfLevels(IC%verts(i)%vert)  ! Is the level value stored ?
        if(fu_cmp_levs_eq(fu_level(IC%verts(i)%vert,k), fu_level(IC%ids(j)))) then
          found = .true.
          exit
        end if
      end do
      if(.not.found) then ! Store new level value
        call add_level(IC%verts(i)%vert, fu_level(IC%ids(j)))
      end if
    else
      IC%NbrVerts = IC%NbrVerts + 1  ! Store the new vertical

      call check_vertlist_size(IC, IC%NbrVerts)
      if(error)return

      call set_vertical(fu_level(IC%ids(j)), IC%verts(IC%NbrVerts)%vert)
      if(error)return
      IC%verts(IC%NbrVerts)%gridPtr = fu_gridPtr(IC, fu_grid(IC%ids(j)))
      IC%verts(IC%NbrVerts)%srcPtr = fu_MetSrcPtr(IC, fu_met_src(IC%ids(j)))
      IC%verts(IC%NbrVerts)%NFlds_vert = 1
      IC%verts(IC%NbrVerts)%flds(1) = j
      if(error)then
        call set_error('Failed','count_verticals')
        return
      endif
    end if
  end do   ! j=1,NbrOfFlds
  !
  !  Sorting the vertical structures
  !
  found = .true.
  do while (found)
    found = .false.
    do i=1, IC%NbrVerts-1
      if(IC%verts(i)%NFlds_vert < IC%verts(i+1)%NFlds_vert) then
        found  = .true.
        call exchange_vert_lsts(IC%verts(i), IC%verts(i+1))
      end if
    end do
  end do  ! while sorting 

end subroutine count_verticals


!********************************************************************

subroutine check_vertlist_size(IC, iSize)
  !
  ! Checks the size of the local list of vertical structures and enlarges it if needed
  !
  implicit none

  ! Imported parameters
  type(Tinput_content), intent(inout) :: IC
  integer, intent(in) :: iSize

  integer :: kret, i
  type(TVert_lst),dimension(:),allocatable :: lstTmp

  if(.not.associated(IC%verts))then
    allocate(IC%verts(5*(int(real(iSize)/5.)+1)),stat=kret)
    if(kret /= 0)then
      call set_error('No more space available','analyse_grib_content')
      return
    end if
    return
  end if

  if(iSize > size(IC%verts))then
    !
    ! Create temporary space and store the content there
    !
    allocate(lstTmp(size(IC%verts)),stat=kret)
    if(kret /= 0)then
      call set_error('No more space available','analyse_grib_content')
      return
    end if
    do i=1,size(IC%verts)
     call set_vertical(fu_level(IC%verts(i)%vert,1), lstTmp(i)%vert) ! Create new
     if(error)return
     do kret = 2,fu_NbrOfLevels(IC%verts(i)%vert)
        call add_level(lstTmp(i)%vert,fu_level(IC%verts(i)%vert,kret))
      end do
      lstTmp(i)%top = IC%verts(i)%top
!      lstTmp(i)%n_layers = IC%verts(i)%n_layers
      lstTmp(i)%gridPtr= IC%verts(i)%gridPtr
      lstTmp(i)%srcPtr = IC%verts(i)%srcPtr
      lstTmp(i)%NFlds_vert = IC%verts(i)%NFlds_vert
      lstTmp(i)%flds(:) = IC%verts(i)%flds(:)
    end do
    !
    ! Resize the vertical structure list and drop the content back
    !
    deallocate(IC%verts)
    allocate(IC%verts(10*(int(real(iSize)/10.)+1)),stat=kret)
    if(kret /= 0)then
      call set_error('No more space available','analyse_grib_content')
      return
    end if

    do i=1,size(lstTmp)
      call set_vertical(fu_level(lstTmp(i)%vert,1), IC%verts(i)%vert) ! Create new
      if(error)return
      do kret = 2,fu_NbrOfLevels(lstTmp(i)%vert)
        call add_level(IC%verts(i)%vert,fu_level(lstTmp(i)%vert,kret))
      end do
      IC%verts(i)%top = lstTmp(i)%top
!      IC%verts(i)%n_layers = lstTmp(i)%n_layers
      IC%verts(i)%gridPtr = lstTmp(i)%gridPtr
      IC%verts(i)%srcPtr = lstTmp(i)%srcPtr
      IC%verts(i)%NFlds_vert = lstTmp(i)%NFlds_vert
      IC%verts(i)%flds(:) = lstTmp(i)%flds(:)
    end do

    deallocate(lstTmp)

  end if
end subroutine check_vertlist_size



!***********************************************************************

subroutine exchange_vert_lsts(vert1, vert2)
  !
  ! Exchanges two vertical structures 
  !
  implicit none

  ! Imported parameters
  type(Tvert_lst), intent(inout) :: vert1, vert2

  ! Local variables
  type(Tvert_lst) :: vTmp

  vTmp%vert = vert1%vert
  vert1%vert = vert2%vert
  vert2%vert = vTmp%vert

  vTmp%NFlds_vert = vert1%NFlds_vert
  vert1%NFlds_vert = vert2%NFlds_vert
  vert2%NFlds_vert = vTmp%NFlds_vert

  vTmp%top = vert1%top
  vert1%top = vert2%top
  vert2%top = vTmp%top

!  vTmp%n_layers = vert1%n_layers
!  vert1%n_layers = vert2%n_layers
!  vert2%n_layers = vTmp%n_layers

  vTmp%srcPtr = vert1%srcPtr
  vert1%srcPtr = vert2%srcPtr
  vert2%srcPtr = vTmp%srcPtr

  vTmp%gridPtr = vert1%gridPtr
  vert1%gridPtr = vert2%gridPtr
  vert2%gridPtr = vTmp%gridPtr

!  vTmp%levs(:) = vert1%levs(:)
!  vert1%levs(:) = vert2%levs(:)
!  vert2%levs(:) = vTmp%levs(:)

  vTmp%flds(1:size(vTmp%flds)) = vert1%flds(:)
  vert1%flds(1:size(vTmp%flds)) = vert2%flds(:)
  vert2%flds(1:size(vTmp%flds)) = vTmp%flds(:)

end subroutine exchange_vert_lsts



!*******************************************************************

integer function fu_vertPtr(IC, level, grid, mds)
  !
  ! Finds the index of the vertical for the given vertica level, grid and data 
  ! source . Returns int_missing if not found
  !
  implicit none

  ! Imported parameters with intent IN
  type(Tinput_content), intent(in) :: IC
  type(silja_grid), intent(in) :: grid
  type(silja_level), intent(in) :: level
  type(meteo_data_source), intent(in) :: mds

  ! Local variables
  integer :: i
  !
  ! Some stupidity check first
  !
  fu_vertPtr = int_missing
  if(IC%NbrVerts < 1)then
    call set_error('No verticals in the list','fu_vertPtr')
    return
  end if

  do i=1, IC%NbrVerts
    if(fu_leveltype(IC%verts(i)%vert) == fu_leveltype(level) .and. &
     & IC%grids(IC%verts(i)%gridPtr)%grid == grid .and. &
     & IC%MetSrcs(IC%verts(i)%srcPtr)%src == mds) then
      fu_vertPtr = i
      return
    end if
  end do
  call set_error('Vertical is not in the list','fu_gridPtr')

end function fu_vertPtr


!***********************************************************************

subroutine count_variables(IC)
  !
  ! Makes a list of variables in the grib content structure. 
  ! Variable is a combination of the quantity, data soruce,
  ! grid and vertical structure. In fact, this combination is not
  ! unique because there may be several levels in the vertical,
  ! so some care should be taken in use of this combination as an 
  ! index
  !
  implicit none

  ! Imported variable
  type(Tinput_content), intent(inout) :: IC

  ! Local variables
  integer :: i,j, k, iVar
  logical :: found

  ! Since variable belongs to the source and may exist only over grid with some
  ! concrete vertical structure - they all have to be defined before
  !
  if(IC%NbrGrids < 1 .or. IC%NbrMetSrcs < 1 .or. IC%NbrVerts < 1)then
    call set_error('Sources, grids or verticals are not counted yet','count_variables')
    return
  end if

  do j=1,IC%NbrFields
    if(.not.fu_known_quantity(fu_quantity(IC%ids(j)))) then
      call set_error('Strange quantity given','count_variables')
      cycle
    end if

    found = .false.
    !
    ! Is this variable already in the list ? 
    ! Check quantity, vertical type, source and grid
    !
    do i=1,IC%NbrVars
      if(IC%vars(i)%quantity == fu_quantity(IC%ids(j)))then
        if(fu_leveltype(IC%verts(IC%vars(i)%vertPtr)%vert) == &
                             & fu_leveltype(fu_level(IC%ids(j)))) then
          if(IC%grids(IC%vars(i)%gridPtr)%grid == fu_grid(IC%ids(j)))then
            if(IC%MetSrcs(IC%vars(i)%MetSrcPtr)%src == fu_met_src(IC%ids(j)))then
              found = .true.
              iVar = i
              exit
            end if
          end if
        end if
      end if
    end do

    if(found) then 
      !
      ! Just different level of already stored variable. Another option is:
      ! the same field was loaded twice e.g. from two files. In this case
      ! the level will correspond too
      !
      if(IC%vars(iVar)%vertLevNbr /= int_missing)then ! So far, this var was treated as a single-level one

        if(.not. fu_cmp_levs_eq(fu_level(IC%verts(IC%vars(iVar)%vertPtr)%vert, IC%vars(iVar)%vertLevNbr), &
               & fu_level(IC%ids(j)))) then
          IC%vars(iVar)%NFlds_var = IC%vars(iVar)%NFlds_var + 1
          IC%vars(iVar)%flds(IC%vars(iVar)%NFlds_var) = j
          IC%vars(iVar)%vertLevNbr = int_missing ! More than one level exists
        endif

      else  ! Already known to be multi-level quantity
        IC%vars(iVar)%NFlds_var = IC%vars(iVar)%NFlds_var + 1
        IC%vars(iVar)%flds(IC%vars(iVar)%NFlds_var) = j
      endif

    else
      IC%NbrVars = IC%NbrVars + 1  ! Store the new variable

      call check_varlist_size(IC, IC%NbrVars)
      if(error)return

      IC%vars(IC%NbrVars)%quantity = fu_quantity(IC%ids(j))
      IC%vars(IC%NbrVars)%gridPtr = fu_gridPtr(IC, fu_grid(IC%ids(j)))
      IC%vars(IC%NbrVars)%MetSrcPtr = fu_MetSrcPtr(IC, fu_met_src(IC%ids(j)))
      IC%vars(IC%NbrVars)%vertPtr = fu_vertPtr(IC, fu_level(IC%ids(j)), &
                                                 & fu_grid(IC%ids(j)), &
                                                 & fu_met_src(IC%ids(j)))
      IC%vars(IC%NbrVars)%vertLevNbr = fu_level_nbr(IC%verts(IC%vars(IC%NbrVars)%vertPtr)%vert, &
                                                  & fu_level(IC%ids(j)))
      if(IC%vars(IC%NbrVars)%vertLevNbr < 0)then
        call msg('')
        call report(IC%ids(j))
        call report(IC%verts(IC%vars(IC%NbrVars)%vertPtr)%vert)
        call set_error('Failed to find the level in vertical','count_variables')
      endif
      IC%vars(IC%NbrVars)%NFlds_var = 1
      IC%vars(IC%NbrVars)%flds(1) = j
    end if
  end do    ! Cycle through fields

end subroutine count_variables



!********************************************************************

subroutine check_varlist_size(IC, iSize)
  !
  ! Checks the size of the local list of vertical structures and enlarges it if needed
  !
  implicit none

  ! Imported parameters
  type(Tinput_content), intent(inout) :: IC
  integer, intent(in) :: iSize

  integer :: kret, i
  type(TVar_lst),dimension(:),allocatable :: lstTmp

  if(.not.associated(IC%vars))then
    allocate(IC%vars(5*(int(real(iSize)/5.)+1)),stat=kret)
    if(kret /= 0)then
      call set_error('No more space available','check_varlist_size')
      return
    end if
    return
  end if

  if(iSize > size(IC%vars))then
    !
    ! Create temporary space and store the content there
    !
    allocate(lstTmp(size(IC%vars)),stat=kret)
    if(kret /= 0)then
      call set_error('No more space available','check_varlist_size')
      return
    end if
    do i=1,size(IC%vars)
      lstTmp(i)%quantity = IC%vars(i)%quantity
      lstTmp(i)%gridPtr = IC%vars(i)%gridPtr
      lstTmp(i)%MetSrcPtr = IC%vars(i)%MetSrcPtr
      lstTmp(i)%vertPtr = IC%vars(i)%vertPtr
      lstTmp(i)%NFlds_var = IC%vars(i)%NFlds_var
      lstTmp(i)%vertLevNbr = IC%vars(i)%vertLevNbr
      lstTmp(i)%flds(:) = IC%vars(i)%flds(:)
    end do
    !
    ! Resize the vertical structure list and drop the content back
    !
    deallocate(IC%vars)
    allocate(IC%vars(10*(int(real(iSize)/10.)+1)),stat=kret)
    if(kret /= 0)then
      call set_error('No more space available','check_varlist_size')
      return
    end if

    do i=1,size(lstTmp)
      IC%vars(i)%quantity = lstTmp(i)%quantity
      IC%vars(i)%gridPtr = lstTmp(i)%gridPtr
      IC%vars(i)%MetSrcPtr = lstTmp(i)%MetSrcPtr
      IC%vars(i)%vertPtr = lstTmp(i)%vertPtr
      IC%vars(i)%NFlds_var = lstTmp(i)%NFlds_var
      IC%vars(i)%vertLevNbr = lstTmp(i)%vertLevNbr
      IC%vars(i)%flds(:) = lstTmp(i)%flds(:)
    end do

    deallocate(lstTmp)

  end if
end subroutine check_varlist_size


!***********************************************************************

integer function fu_compare_vars_quality(IC, ind1, ind2)
  !
  ! Compares two variables and selects the best of them, which index
  ! is then returned. The set of criteria is based on the resolution 
  ! of horizontal grid, coverage and resolution of the vertical grid
  ! and the number of fields with similar features
  ! Horizontal coverage does not play any role because the system_grid
  ! will be taken to be the smallest of those passed filtration
  !
  implicit none

  ! Imported parameters with intent IN
  type(Tinput_content), intent(in) :: IC
  integer, intent(in) :: ind1, ind2

  ! Local variables
  real :: grid_resolution_weight = 1., &
        & vertical_coverage_weight = 0.5, &
        & vertical_resolution_weight = 0.7, &
        & number_of_similar_fields_weight = 1.1, &
        & quality1, quality2, factor1, factor2
  type(Tvar_lst), pointer :: var1, var2
  integer :: i

  var1 => IC%vars(ind1);  var2 => IC%vars(ind2)
  quality1=0.;  quality2 = 0.

  !
  ! Horizontal resolution - the most important parameter.
  ! The smaller cell size the better grid is
  !
  factor1 = fu_cell_size(IC%grids(var1%gridPtr)%grid)
  factor2 = fu_cell_size(IC%grids(var2%gridPtr)%grid)

  if(factor1 < 0.8 * factor2) then
    quality1 = quality1 + grid_resolution_weight * 10.

  elseif(factor1 < 0.9 * factor2)then
    quality1 = quality1 + grid_resolution_weight * 2.
    
  elseif(factor1 < 0.999 * factor2)then
    quality1 = quality1 + grid_resolution_weight

  elseif(factor1 .eps. factor2)then

  elseif(factor1 < 1.11 * factor2)then
    quality2 = quality2 + grid_resolution_weight

  elseif(factor1 < 1.25 * factor2)then
    quality2 = quality2 + grid_resolution_weight * 2.

  else
    quality2 = quality2 + grid_resolution_weight * 10.

  end if

  if(error)return

  !
  ! Verticals are not always comparable.
  !
  if(fu_verts_comparable(IC%verts(var1%vertPtr)%vert, IC%verts(var2%vertPtr)%vert))then
    !
    ! Vertical coverage. Computed via absolute ranges between top and bottom
    ! levels. Note - the units can depend on the vertical type but must be the same
    ! for both variables. The larger range the better vertical is
    !
    factor1 = 1.0 / (fu_range(IC%verts(var1%vertPtr)%vert) + 0.0001)
    factor2 = 1.0 / (fu_range(IC%verts(var2%vertPtr)%vert) + 0.0001)
    if(factor1 < 0.5 * factor2) then
      quality1 = quality1 + vertical_coverage_weight * 20.

    elseif(factor1 < 0.7 * factor2)then
      quality1 = quality1 + vertical_coverage_weight * 5.
    
    elseif(factor1 < 0.999 * factor2)then
      quality1 = quality1 + vertical_coverage_weight

    elseif (factor1 .eps. factor2)then

    elseif(factor1 < 1.3 * factor2)then
      quality2 = quality2 + vertical_coverage_weight

    elseif(factor1 < 2. * factor2)then
      quality2 = quality2 + vertical_coverage_weight * 5.

    else
      quality2 = quality2 + vertical_coverage_weight * 20.

    end if
    !
    ! Vertical resolution is determined as a covered vertical range divided 
    ! by the number of vertical levels. This, of course, contradicts to the
    ! coverage criteria, but this contradiction seems to be natural.
    ! The smaller mean layer thickness the better vertical is
    !
    factor1 = 1. / (fu_NbrOfLevels(IC%verts(var1%vertPtr)%vert) * factor1 + 0.0001)
    factor2 = 1. / (fu_NbrOfLevels(IC%verts(var2%vertPtr)%vert) * factor2  + 0.0001)

    if(factor1 < 0.5 * factor2) then
      quality1 = quality1 + vertical_resolution_weight * 20.

    elseif(factor1 < 0.7 * factor2)then
      quality1 = quality1 + vertical_resolution_weight * 5.
    
    elseif(factor1 < 0.999 * factor2)then
      quality1 = quality1 + vertical_resolution_weight

    elseif (factor1 .eps. factor2)then

    elseif(factor1 < 1.3 * factor2)then
      quality2 = quality2 + vertical_resolution_weight

    elseif(factor1 < 2. * factor2)then
      quality2 = quality2 + vertical_resolution_weight * 5.

    else
      quality2 = quality2 + vertical_resolution_weight * 20.

    end if

  end if

  !
  ! Number of similar fields is somewhat strange parameter indirectly controlling
  ! the amount of interpolations (first of all, horizontal) and internal 
  ! inconsistencies between the data from different sources
  !
  factor1=0.; factor2=0.

  do i=1, IC%NbrVars
    if(IC%vars(i)%gridPtr == var1%gridPtr .and. &
     & IC%vars(i)%MetSrcPtr == var1%MetSrcPtr) factor1 = factor1 + 1.
    if(IC%vars(i)%gridPtr == var2%gridPtr .and. &
     & IC%vars(i)%MetSrcPtr == var2%MetSrcPtr) factor2 = factor2 + 1.
  end do

  factor1 = 1./(factor1+0.001)
  factor2 = 1./(factor2+0.001)

  if(factor1 < 0.1 * factor2) then
    quality1 = quality1 + number_of_similar_fields_weight * 100.

  elseif(factor1 < 0.5 * factor2)then
    quality1 = quality1 + number_of_similar_fields_weight * 10.
    
  elseif(factor1 < 0.999 * factor2)then
    quality1 = quality1 + number_of_similar_fields_weight

  elseif (factor1 .eps. factor2)then

  elseif(factor1 < 2. * factor2)then
    quality2 = quality2 + number_of_similar_fields_weight

  elseif(factor1 < 10. * factor2)then
    quality2 = quality2 + number_of_similar_fields_weight * 10.

  else
    quality2 = quality2 + number_of_similar_fields_weight * 100.
  end if

  !
  ! Final comparison
  !
  if(quality1 > quality2)then
    fu_compare_vars_quality = ind1
  else
    fu_compare_vars_quality = ind2
  end if

end function fu_compare_vars_quality


!***************************************************************************

subroutine remove_variable (IC, indexVar)
  ! 
  ! Removes the variable and all related fields from the GRIB content structure
  !
  implicit none

  ! Imported variables with intent IN
  integer, intent(in) :: indexVar

  ! Imported parameters with intent INOUT
  type(Tinput_content), intent(inout) :: IC

  ! local variable
  integer :: iTmp

  if(.not. IC%defined == silja_true) then
    call set_error('Undefined content given','remove_variable')
    return
  end if

  if(indexVar < 1 .or. indexVar > IC%NbrVars) then
    call msg_warning('Strange index given. Skipped','remove_variable')
    return
  end if

  ! Kill all ids pointed by this variable, but not move ids - otherwise you
  ! destroy the other structural links
  !
  do iTmp = 1, IC%vars(indexVar)%NFlds_var
    call set_missing(IC%ids(IC%vars(indexVar)%flds(iTmp)))
  end do
  IC%vars(indexVar)%flds(:) = int_missing
  IC%vars(indexVar)%quantity = int_missing
  IC%vars(indexVar)%gridPtr = int_missing
  IC%vars(indexVar)%MetSrcPtr = int_missing
  IC%vars(indexVar)%vertPtr = int_missing
  IC%vars(indexVar)%NFlds_var = int_missing

end subroutine remove_variable

END MODULE input_analysis
