MODULE fields ! 2-D scalar-data with identification

  ! Description: 
  ! Contains the definition of silja_field -type which contains one
  ! fully definedfield. 
  ! Horizontal grid or vertical level can be of any type defined in
  ! modules grids and levels.
  !
  ! Author: Mika Salonoja, FMI, email Mika.Salonoja@fmi.fi
  !
  ! All u its: SI
  ! 
  ! Language: ANSI standard Fortran 90
  !
  ! 
  ! Modules used:

  USE field_identifications
  USE positions_v2

  IMPLICIT NONE

  ! The public functions and subroutines available in this module:
  PUBLIC set_field
  PUBLIC set_field_empty
  PUBLIC set_new_field
  public make_test_field
  PUBLIC fu_id
  PUBLIC defined
  PUBLIC fu_field_defined
  PUBLIC fu_size
  PUBLIC report
  PUBLIC fu_interpolate_to_position
  PUBLIC fu_met_src
  PUBLIC fu_grid
  PUBLIC fu_valid_time
  PUBLIC fu_analysis_time
  PUBLIC fu_forecast_length
  PUBLIC fu_accumulation_length
  PUBLIC fu_validity_length
  PUBLIC fu_accumulation_start_time
  PUBLIC fu_level
  PUBLIC fu_quantity
  PUBLIC fu_grid_data
  PUBLIC fu_accumulated
  public fu_period_valid
  public fu_species
  public fu_substance_name
  public fu_mode
  public fu_optical_wave_length
  public interpolate_fields_in_time ! Handles different type of field averaging
  public set_valid_time
  public update_field_data
  public fu_field_stats

  ! The private functions and subroutines:
  private set_field_with_data
  private set_field_no_data
  PRIVATE fu_field_id
  PRIVATE fu_met_src_of_field
  PRIVATE fu_grid_of_field
  PRIVATE fu_valid_time_of_field
  PRIVATE fu_forecast_length_of_field
  PRIVATE fu_analysis_time_of_field
  PRIVATE fu_level_of_field
  PRIVATE fu_quantity_of_field
  PRIVATE fu_accumulation_of_field
  private fu_validity_length_of_field
  PRIVATE fu_fieldpointer_defined
  PRIVATE fu_interpolate_scalar_to_pos
  PRIVATE fu_fieldsize
  PRIVATE print_fieldpointer_report
  PRIVATE fu_accum_start_time_field
  PRIVATE fu_accumulated_field
  PRIVATE fu_period_valid_field
  private set_valid_time_of_field
  private fu_species_of_field
  private fu_substance_name_of_field
  private fu_mode_of_field
  private fu_optic_wavelen_of_field

  interface set_field
    module procedure set_field_with_data
    module procedure set_field_no_data
  end interface

  INTERFACE defined
    MODULE PROCEDURE fu_fieldpointer_defined
  END INTERFACE

  INTERFACE report
    MODULE PROCEDURE print_fieldpointer_report
  END INTERFACE

  INTERFACE fu_interpolate_to_position
    MODULE PROCEDURE fu_interpolate_scalar_to_pos
  END INTERFACE

  INTERFACE fu_id 
    MODULE PROCEDURE fu_field_id
  END INTERFACE

  INTERFACE fu_size
    MODULE PROCEDURE fu_fieldsize
  END INTERFACE

  INTERFACE fu_met_src  
    MODULE PROCEDURE fu_met_src_of_field
  END INTERFACE

  INTERFACE fu_grid 
    MODULE PROCEDURE fu_grid_of_field
  END INTERFACE

  INTERFACE fu_valid_time
    MODULE PROCEDURE fu_valid_time_of_field
  END INTERFACE

  INTERFACE fu_analysis_time
    MODULE PROCEDURE fu_analysis_time_of_field
  END INTERFACE

  INTERFACE fu_forecast_length
    MODULE PROCEDURE fu_forecast_length_of_field
  END INTERFACE

  INTERFACE fu_accumulation_length
    MODULE PROCEDURE fu_accumulation_of_field
  END INTERFACE

  INTERFACE fu_validity_length
    module procedure fu_validity_length_of_field
  END INTERFACE

  INTERFACE fu_accumulated
    MODULE PROCEDURE fu_accumulated_field
  END INTERFACE

  INTERFACE fu_period_valid
    MODULE PROCEDURE fu_period_valid_field
  END INTERFACE

  INTERFACE fu_accumulation_start_time
    MODULE PROCEDURE fu_accum_start_time_field
  END INTERFACE

  INTERFACE fu_level
    MODULE PROCEDURE fu_level_of_field
  END INTERFACE

  INTERFACE fu_quantity
    MODULE PROCEDURE fu_quantity_of_field
  END INTERFACE

  interface fu_species
    module procedure fu_species_of_field
  end interface

  interface fu_substance_name
    module procedure fu_substance_name_of_field
  end interface

  interface fu_mode
    module procedure fu_mode_of_field
  end interface

  interface fu_optical_wave_length
    module procedure fu_optic_wavelen_of_field
  endinterface

  interface set_valid_time
    module procedure set_valid_time_of_field
  end interface

  ! Public types with private components defined in this module:
  TYPE silja_field ! normal scalar 2D-field
    PRIVATE
    TYPE(silja_field_id) :: id
    REAL, DIMENSION(:), POINTER :: grid_data ! scalar data.
    TYPE(silja_logical) :: defined
  END TYPE silja_field

  ! Public types with public components defined in this module:
  TYPE silja_fieldpointer
    TYPE(silja_field), POINTER :: fp
  END TYPE silja_fieldpointer

  integer, public, parameter :: create_field = 701, &          ! if absent create; if present error
                              & overwrite_field = 702, &       ! if absent create; if present overwrite
                              & overwrite_field_anytime=703,&  ! aggressive search ignoring time stamp
                              & add_field = 704, &             ! if absent create; if present add. 
                              & create_if_absent = 705, &      ! if absent create; if present warning and do nothing
                              & overwrite_forced = 706         ! if absent error
CONTAINS 

  ! ***************************************************************


  SUBROUTINE set_field_with_data(id, grid_data, field, ifNew)
    
    ! Description:
    ! Sets a value and status for a field. Allocates memory for the
    ! grid-data. If there's alredy something in the field, it is
    ! overwritten. In this case memory area is not deallocated and
    ! then allocated again, if the sizes of the old and new data
    ! happen to be the same.
    !
    ! Negative values are also filtered from some weather quantites
    ! (for example humidity and rain). There sometimes are negative
    ! values in nwp-fields due to numerics.
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_field_id), INTENT(in) :: id
    REAL, DIMENSION(:), intent(in) :: grid_data
    logical, intent (in) :: ifNew
    
    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_field), INTENT(inout) :: field

    ! Local declarations:
    INTEGER :: gridpoints, i, n_out_of_range, n_patched, nx, ny, n_failed

    !--------------------------------------------------------------
    !
    ! 1. Check input parameters id and grid separately (grid can still be vois in the defined id).
    !
    IF (fu_fails(defined(id),'identification given is not ok', 'set_field_with_data'))RETURN
    IF (fu_fails(defined(fu_grid(id)),'the grid of field-id not defined', 'set_field_with_data'))RETURN
 
 
    !--------------------------------------------------------------
    !
    ! 2. See that there's memory for the grid-data. 
    !    -----------------------------------------
  
    ! The number of gridpoints to be stored is found in the grid:
    call grid_dimensions(fu_grid(id), nx, ny)
    gridpoints = nx*ny

!print *, 'set field: 2: ', gridpoints

    IF (gridpoints > SIZE(grid_data)) THEN
      CALL set_error('data-array smaller than number of gridpoints', 'set_field_with_data')
      RETURN
    END IF

!print *, 'set field: 2.1: ', SIZE(field%grid_data)
    if(ifNew) nullify(field%grid_data)
    CALL set_array_size(field%grid_data, gridpoints)
!print *, 'set field: 2.2: ', SIZE(field%grid_data)
 
    IF (error) RETURN

    !    PRINT *, ' Field size set: ', SIZE(field%grid_data)

    !
    field%grid_data(1:gridpoints) = grid_data(1:gridpoints)

    ! Set the stuff filtering the out-of-range values if needed
    call check_quantity_range(fu_quantity(id), field%grid_data, gridpoints, nx, &
                                     & .false., .false., &  ! ifRequireValidity, ifSilent
                                     & n_out_of_range, n_patched, n_failed)
    if(error)return
    
#ifdef DEBUG
    if (n_out_of_range > 0) then
      call msg('Out-of-range values removed:' + fu_quantity_string(fu_quantity(id)), n_out_of_range)
    end if
#endif
    
    field%id = id

    field%defined = fu_set_true()
    
!        CALL report(field%id)


  END SUBROUTINE set_field_with_data
  
  
  ! ***************************************************************


  SUBROUTINE set_field_no_data(id, pField, ifNew, ifResetVals)
    ! 
    ! Sets the field metadata and reserves the space for real data. 
    ! If there's alredy something in the field, memory area is not 
    ! deallocated and then allocated again, if the sizes of the old and new data
    ! happen to be the same.
    !
    IMPLICIT NONE
    !
    ! Imported parameters
    TYPE(silja_field_id), INTENT(in) :: id
    logical, intent (in) :: ifNew
    logical, intent(in), optional :: ifResetVals
    TYPE(silja_field), target :: pField
    !
    ! Check input parameters.
    !
    IF (fu_fails(defined(id),'field id is undefined', 'set_field_no_data'))return
    IF (fu_fails(defined(fu_grid(id)),'the grid of field-id is undefined', 'set_field_no_data'))return
    !
    ! Set the grid_data size, if turns necessary
    !
    if(ifNew)nullify(pField%grid_data)
    CALL set_array_size(pField%grid_data, fu_number_of_gridpoints(fu_grid(id)))
    IF (error) RETURN
    !
    ! Values are NOT reset only if you forbid it
    !
    if(present(ifResetVals))then
      if(ifResetVals) pField%grid_data (1:fu_number_of_gridpoints(fu_grid(id))) = real_missing
    else
      pField%grid_data (1:fu_number_of_gridpoints(fu_grid(id))) = real_missing
    endif
    !
    ! Finally, set the id and raise the defined flag
    !
    pField%id = id
    pField%defined = fu_set_true()

  END SUBROUTINE set_field_no_data


  ! ***************************************************************

  SUBROUTINE set_field_empty(field, deallocate_memory)
    
    ! Description:
    ! Sets field empty. If deallocate_memory, then the actual data
    ! -field is deallocated. If false, then only identifications are
    ! set indicating a missing, undefined field.
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent INOUT or POINTER:
    LOGICAL, INTENT(in) :: deallocate_memory

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_field), INTENT(inout) :: field
    
    ! Local declarations:
    
!    print *,'set_field:',deallocate_memory,ASSOCIATED(field%grid_data),SIZE(field%grid_data)

    IF(deallocate_memory)THEN
      IF (ASSOCIATED(field%grid_data)) THEN
        IF(SIZE(field%grid_data)>0 .and. size(field%grid_data) <= worksize+1)THEN
          CALL free_array_memory(field%grid_data)
        END IF
      END IF
      NULLIFY(field%grid_data) !Even non-associated pointer may be wrong-filled
    END IF

!    print *,'done',deallocate_memory,ASSOCIATED(field%grid_data),SIZE(field%grid_data)

    call set_missing(field%id)
    
    field%defined = fu_set_false()
    
  END SUBROUTINE set_field_empty


  ! ***************************************************************


  SUBROUTINE set_new_field(field)
    
    ! Sets new field. The reason for the specific function is - sometimes new
    ! field can be filled with garbage causing erroneous ASSOCIATED and SIZE
    ! answers. Indeed, it requires the object already allocated or nullified, otherwise 
    ! result is unpredicteable. For new field the data_array is not
    ! allocated and not associated, which has to be explicitly set.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent INOUT:
    TYPE(silja_field), INTENT(inout) :: field
    
    
!    print *,'set_new_field:',ASSOCIATED(field%grid_data),SIZE(field%grid_data)

    NULLIFY(field%grid_data) !Even non-associated pointer may be wrong-filled

!    print *,'set_new_field done',ASSOCIATED(field%grid_data),SIZE(field%grid_data)

    call set_missing(field%id)
    
    field%defined = fu_set_false()
    
  END SUBROUTINE set_new_field


  !***********************************************************************

  subroutine make_test_field(chFName, id, grid_data, valid_time, fc_len, grid)
    !
    ! Sometimes a test field filled with a constant value may be needed.
    ! This sub creates a simple field from its string description.
    ! All what is needed is - quantity and level names, as well as the value
    ! Time has to be given. Later it also can come from describing string
    ! The grid is assumed to be meteorological
    !
    implicit none

    ! Imported parameters
    character(len=*), intent(in) :: chFName
    type(silja_field_id), intent(out) :: id
    real, dimension(:), intent(out) :: grid_data
    type(silja_time), intent(in) :: valid_time
    type(silja_interval), intent(in) :: fc_len
    type(silja_grid), intent(in) :: grid

    ! Local variables
    integer :: iStatus, quantity, fs
    type(silam_species) :: species
    type(silja_level) :: level
    character (len=fnlen) :: spTmp

    !
    ! A bit of preparation
    !
    if(.not. defined(grid))then
      call set_error('Undefined grid','make_test_field')
      return
    endif
    fs = fu_number_of_gridpoints(grid)
    if (size(grid_data) < fs) then
      call msg("size(grid_data)", size(grid_data))
      call report(grid)
      call set_error('Too small storage for grid','make_test_field')
      return
    endif

    !
    ! Let's start from quantity
    !
    read(unit=chFName, iostat=iStatus, fmt=*) spTmp
    if(iStatus /= 0)then
      call set_error('Failed to read quantity name from:' + chFName, 'make_test_field')
      return
    endif

    call decode_id_params_from_io_str(spTmp, &     ! string to decode
                                    & .false., &            ! ifMultiLevel
                                    & quantity, &     ! decoded quantity
                                    & species, &      ! decoded species
                                    & .true.)         ! scream if fail
    if(quantity == int_missing)then
      call set_error('Strange quantity name in:' + chFName, 'make_test_field')
      return
    endif
    !
    ! Now set the level
    !
    if(index(chFName,'SURFACE_LEVEL') > 0)then
      level = surface_level
    elseif(index(chFName,'2M_LEVEL') > 0)then
      level = level_2m_above_ground
    elseif(index(chFName,'10M_LEVEL') > 0)then
      level = level_10m_above_ground
    elseif(index(chFName,'ALL_LEVELS') > 0)then
      level = entire_atmosphere_mean_level
    else
      call msg('Allowed levels: SURFACE_LEVEL, 2M_LEVEL, 10M_LEVEL, ALL_LEVELS')
      call set_error(fu_connect_strings('Unknown level in:',chFName),'make_test_field')
      return
    endif
    !
    ! Finally, the value itself
    !
    read(unit=chFName, iostat=iStatus, fmt=*) spTmp, spTmp, grid_data(1)
    if(iStatus /= 0)then
      call set_error('Failed to read value from:' + chFName, 'make_test_field')
      return
    endif

    grid_data(2:fs) = grid_data(1)
    !
    ! Finally, set the id
    !
    id = fu_set_field_id(fmi_silam_src, &   !silam_internal_src, &
                       & quantity, &
                       & valid_time - fc_len, &  ! analysis
                       & fc_len, & 
                       & grid, &
                       & level,&
                       & species = species)

  end subroutine make_test_field

  ! ***************************************************************

  
  LOGICAL FUNCTION fu_fieldpointer_defined(field)
    !
    ! Description:
    ! Returns a true value, if the field has been filled with data and
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    
    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_field), POINTER :: field
    TYPE(silja_field), POINTER :: fieldTmp

    fieldTmp => field

    IF (ASSOCIATED(field)) THEN    
      fu_fieldpointer_defined = fu_true(field%defined)
    ELSE
      fu_fieldpointer_defined = .false.
    END IF

  END FUNCTION fu_fieldpointer_defined

  
  ! ***************************************************************

  
  LOGICAL FUNCTION fu_field_defined(field)
    !
    ! Returns a true value, if the field has been filled with data
    ! NOTE. It is an addition to the fu_fieldpointer_defined, but with
    ! different argument. It is totally unclear why the difference exists
    ! but it seems to be so
    !
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE
    
    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_field), INTENT(in) :: field

    fu_field_defined = fu_true(field%defined)

  END FUNCTION fu_field_defined
  

  ! ***************************************************************

  
  FUNCTION fu_grid_data(field)
    !
    ! Description:
    ! Returns the level of a field-id.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    REAL, DIMENSION(:), POINTER :: fu_grid_data
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_field), target :: field

    IF (fu_field_defined(field)) THEN  
      fu_grid_data => field%grid_data
    ELSE
      CALL set_error('undefined field given','fu_grid_data')
    END IF

  END FUNCTION fu_grid_data


  !**********************************************************************

  subroutine interpolate_fields_in_time(fld_id_past, data_past, &
                                      & fld_id_future, data_future, &
                                      & weight_past, &
                                      & fld_id_res, data_res, &
                                      & now)
    !
    ! Performs a time interpolation of two fields, taking into account
    ! their types - instant, averaged, accumulated, period-valid, etc.
    !
    implicit none

    ! Imported parameters
    !
    type(silja_field_id), intent(in) :: fld_id_past, fld_id_future
    type(silja_field_id), intent(out) :: fld_id_res
    real, dimension(:), intent(in) :: data_past, data_future
    real, dimension(:),intent(out) ::data_res
    real, intent(in) :: weight_past
    type(silja_time), intent(in) :: now

    ! Local variables
    !
    type(silja_interval) :: accLenPast, accLenFuture, validLenPast, validLenFuture, & 
                          & lenValidity, lenAccumulation
    type(silja_time) :: analysisTime
    real :: fSeconds
    integer :: fs

    ! Some preparations and stupidity check
    !
    call set_missing(fld_id_res)

    if(.not.(defined(fld_id_past) .and.defined(fld_id_future)))then
      call set_error('Undefined field(s) given','interpolate_fields_in_time')
      return
    endif
    if(fu_quantity(fld_id_past) /= fu_quantity(fld_id_future))then
      call set_error('Different field quantities','interpolate_fields_in_time')
      return
    endif
    if(.not. fu_grid(fld_id_past) == fu_grid(fld_id_future))then
      call set_error('Different grids given','interpolate_fields_in_time')
      return
    endif
    if(.not.defined(now))then
      call set_error('Undefined valid time(s)','interpolate_fields_in_time')
      return
    endif
    fs = fu_number_of_gridpoints(fu_grid(fld_id_past))
    if(error)return
    !
    ! One can also call these for single-time fields, which have no time dimension.
    ! Then the IDs are identical. To avoid comparison of anyway identical objects,
    ! let's just check the valid time. If identical, no interpolation anyway.
    !
!    if(fu_valid_time(fld_id_past) == fu_valid_time(fld_id_future))then
    if(fld_id_past == fld_id_future)then
      call msg_warning('Time interpolation is called for single-time field:' + &
                     & fu_quantity_string(fu_quantity(fld_id_past)),'interpolate_fields_in_time')
      data_res(1:fs) = data_past(1:fs)
      fld_id_res = fld_id_past
      return
    endif
    
    analysisTime = fu_analysis_time(fld_id_past)

    !
    ! Validity length and accumulation length need certain care
    ! Fields may have different start of accumulation, they do have
    ! different valid time and may have different length of validity
    !
    if(fu_accumulated_id(fld_id_past) .or. fu_validity_length(fld_id_past) > zero_interval)then
      accLenPast = fu_accumulation_length(fld_id_past)
      accLenFuture = fu_accumulation_length(fld_id_future)
      validLenPast = fu_validity_length(fld_id_past)
      validLenFuture = fu_validity_length(fld_id_future)

      select case (fu_field_kind(fld_id_past))
     
        case(forecast_flag)  ! In fact, never come here
          !
          ! Simple fields
          !
          lenValidity = interval_missing
          lenAccumulation = interval_missing
          data_res(1:fs) = weight_past * data_past + (1.-weight_past)* data_future

        case(accumulated_flag, averaged_flag)
          !
          ! All sorts of accumulated fields
          !
          lenValidity = interval_missing

          if(accLenPast == accLenFuture)then ! Same integration length
            !
            ! Same integration length, different start-end. Keep the length
            !
            lenAccumulation = accLenPast
            data_res(1:fs) = weight_past * data_past + (1.-weight_past)* data_future

          elseif(fu_accumulation_start_time(fld_id_past) == &
               & fu_accumulation_start_time(fld_id_future))then 
            !
            ! Same start of integration, different lengths. Keep the starting point
            !
            lenAccumulation = fu_set_interval_sec(real(nint(fu_sec(accLenPast) * weight_past + &
                            & fu_sec(accLenFuture) * (1.-weight_past))))
            data_res(1:fs) = weight_past * data_past + (1.-weight_past)* data_future

          elseif(now == fu_accumulation_start_time(fld_id_future) .and. &
               & now == fu_valid_time(fld_id_past))then
            !
            ! We are exactly at the point of starting the accumulation for the future field
            ! If should not be used as it is not yet informative. Instead, take the past field
            ! since its validity time is now (e.g., midnight for daily-mean field)
            !
            lenAccumulation = accLenPast
            analysisTime = fu_analysis_time(fld_id_past)
            data_res(1:fs) = data_past(1:fs)
                 
          elseif(fu_accumulation_start_time(fld_id_past) < &
               & fu_accumulation_start_time(fld_id_future))then
            !
            ! New period of accumulation started; totally neglect past and start from 0
            !
            lenAccumulation = now - fu_accumulation_start_time(fld_id_future)
            analysisTime = fu_analysis_time(fld_id_future)
            data_res(1:fs) = data_future * fu_sec(lenAccumulation) / (fu_sec(accLenFuture)+1)

          else
            !
            ! Both length and start/end of accumulation are different
            ! Then we create the field with accumulation starting from the 
            ! accumulation_past and lasting till now. Time after the valid_time_past
            ! is covered via future field
            !
            lenAccumulation = accLenPast + (now - fu_valid_time(fld_id_past))

            if(now < fu_valid_time(fld_id_future) - accLenFuture)then 
              !
              ! Future field is totally in the future, neglect it
              !
              data_res(1:fs) = data_past

            elseif(fu_valid_time(fld_id_future) - accLenFuture < &
                 & fu_valid_time(fld_id_past))then
              fSeconds = fu_sec(now - fu_valid_time(fld_id_past))
              data_res(1:fs) = data_past + data_future * fSeconds / (fu_sec(accLenFuture)+1)
            else
              fSeconds = fu_sec(now - (fu_valid_time(fld_id_future)-accLenFuture))
              data_res(1:fs) = data_past + data_future * fSeconds / (fu_sec(accLenFuture)+1)
            endif
          endif ! comparison of accumulation length and start points in past and future

!        case(period_valid_flag)
!          lenAccumulation = interval_missing
!          lenValidity = validLenPast * weight_past + validLenFuture * (1.-weight_past)

!        case(difference_flag) ! time-difference of two fields
!          call set_error('Time difference of two fields is not supported', &
!                       & 'interpolate_fields_in_time')
!           return

        case default
          call set_error('Unknown field type','interpolate_fields_in_time')
          return
      end select

    else ! accumulated field
            
      lenAccumulation = interval_missing
      data_res(1:fs) = weight_past * data_past + (1.-weight_past)* data_future
    endif ! Accumulated or period valid field

    if(defined(fu_validity_length(fld_id_past)) .and. &
     & defined(fu_validity_length(fld_id_future)))then
      lenValidity = fu_validity_length(fld_id_past) * weight_past + &
                  & fu_validity_length(fld_id_future) * (1.-weight_past)
    else
      lenValidity = interval_missing
    endif

    if(error)return

    !
    ! Just to avoid warning...
    !
    if(lenAccumulation == interval_missing)then
      fld_id_res = fu_set_field_id(fu_met_src(fld_id_past),&
                                 & fu_quantity(fld_id_past), &
                                 & analysisTime, &
                                 & now - analysisTime, & ! forecast_length
                                 & fu_grid(fld_id_past),&
                                 & fu_level(fld_id_past),&
!                                 & lenAccumulation, &
                                 & length_of_validity = lenValidity, &
                                 & field_kind = fu_field_kind(fld_id_past), &
                                 & species = fu_species(fld_id_past), &
                                 & chCocktail = fu_cocktail_name(fld_id_past))
!                                 & chSubstNm = fu_substance_name(fld_id_past))
    else
      fld_id_res = fu_set_field_id(fu_met_src(fld_id_past),&
                                 & fu_quantity(fld_id_past), &
                                 & analysisTime, &
                                 & now - analysisTime, & ! forecast_length
                                 & fu_grid(fld_id_past),&
                                 & fu_level(fld_id_past),&
                                 & lenAccumulation, &
                                 & lenValidity, &
                                 & fu_field_kind(fld_id_past), &
                                 & species = fu_species(fld_id_past), &
                                 & chCocktail = fu_cocktail_name(fld_id_past))
!                                 & fu_substance_name(fld_id_past))
    endif

  end subroutine interpolate_fields_in_time

  ! ***************************************************************

  ! ***************************************************************
  !
  !
  !      Private functions and subroutines.
  !
  !
  ! ***************************************************************
  ! ***************************************************************
  
  FUNCTION fu_grid_of_field(field)

    ! Returns the grid of a field-id.
    !
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_grid) :: fu_grid_of_field
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_field), POINTER :: field

    if(defined(field))then
      fu_grid_of_field = fu_grid(field%id)
    else
      call set_error('undefined field','fu_grid_of_field')
      fu_grid_of_field = grid_missing
    end if
    
  END FUNCTION fu_grid_of_field



  ! ***************************************************************

  
  FUNCTION fu_met_src_of_field(field) result (met_src)

    ! Returns the data source of a field-id.
    !
    IMPLICIT NONE
    !
    type(meteo_data_source) :: met_src
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_field), POINTER :: field

    if(defined(field))then
      met_src = fu_met_src(field%id)
    else
      call set_error('undefined field','fu_met_src_of_field')
      met_src = met_src_missing
    end if
    
  END FUNCTION fu_met_src_of_field



  ! ***************************************************************

  FUNCTION fu_level_of_field(field)
    !
    ! Description:
    ! Returns the level of a field-id.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_level) :: fu_level_of_field
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_field), POINTER :: field
  
    if(defined(field))then
      fu_level_of_field = fu_level(field%id)
    else
      call set_error('undefined field','fu_level_of_field')
      fu_level_of_field = level_missing
    end if

  END FUNCTION fu_level_of_field


  
  ! ***************************************************************

  FUNCTION fu_valid_time_of_field(field)

    ! Description:
    ! Returns the analysis time of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_time) :: fu_valid_time_of_field
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_field), POINTER :: field

    IF (defined(field)) THEN
      fu_valid_time_of_field = fu_valid_time(field%id)
    else
      call set_error('undefined field','fu_valid_time_of_field')
      fu_valid_time_of_field = time_missing
    END IF

  END FUNCTION fu_valid_time_of_field



  ! ***************************************************************

  subroutine set_valid_time_of_field(field, valid_time)
    ! 
    ! Sets the valid time of a field
    !
    IMPLICIT NONE
    !
    ! Imported parameters
    TYPE(silja_time), intent(in) :: valid_time
    TYPE(silja_field), POINTER :: field

    IF (defined(field)) THEN
      call set_valid_time(field%id, valid_time)
    else
      call set_error('undefined field','set_valid_time_of_field')
    END IF

  END subroutine set_valid_time_of_field

   
  
  ! ***************************************************************

  FUNCTION fu_analysis_time_of_field(field)
    
    ! Description:
    ! Returns the analysis time of a field.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_time) :: fu_analysis_time_of_field
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_field), POINTER :: field
    
    IF (defined(field)) THEN
      fu_analysis_time_of_field = fu_analysis_time(field%id)
    else
      call set_error('undefined field','fu_analysis_time_of_field')
      fu_analysis_time_of_field = time_missing
    END IF

  END FUNCTION fu_analysis_time_of_field


   
  ! ***************************************************************

  FUNCTION fu_accum_start_time_field(field)
    
    ! Description:
    ! Returns the start time of accumulation.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_time) :: fu_accum_start_time_field
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_field), POINTER :: field
    
    IF (defined(field)) THEN
      fu_accum_start_time_field = fu_accumulation_start_time(field%id)
    else
      call set_error('undefined field','fu_accum_start_of_field')
      fu_accum_start_time_field = time_missing
    END IF

  END FUNCTION fu_accum_start_time_field


  
  ! ***************************************************************

  FUNCTION fu_forecast_length_of_field(field)
   
    ! The return value of this function:
    TYPE(silja_interval) :: fu_forecast_length_of_field
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_field), POINTER :: field
    
    IF (defined(field)) THEN
      fu_forecast_length_of_field = fu_forecast_length(field%id)
    else
      call set_error('undefined field','fu_forecast_length_of_field')
      fu_forecast_length_of_field= interval_missing
    END IF

  END FUNCTION fu_forecast_length_of_field

   
  ! ***************************************************************

  FUNCTION fu_accumulation_of_field(field)
    
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_interval) :: fu_accumulation_of_field
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_field), POINTER :: field
    
    if(defined(field))then
      fu_accumulation_of_field = fu_accumulation_length(field%id)
    else
      call set_error('undefined field','fu_accumulation_of_field')
      fu_accumulation_of_field = interval_missing
    end if
    
  END FUNCTION fu_accumulation_of_field

   
  ! ***************************************************************

  FUNCTION fu_validity_length_of_field(field)
    
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_interval) :: fu_validity_length_of_field
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_field), POINTER :: field
    
    if(defined(field))then
      fu_validity_length_of_field = fu_validity_length(field%id)
    else
      call set_error('undefined field','fu_validity_length_of_field_id')
      fu_validity_length_of_field = interval_missing
    end if
    
  END FUNCTION fu_validity_length_of_field

 
  ! ***************************************************************

  INTEGER FUNCTION fu_quantity_of_field(field)
    
    ! Description:
    ! Returns the quantity of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_field), POINTER :: field
    
    if(defined(field))then
      fu_quantity_of_field = fu_quantity(field%id)
    else
      call set_error('undefined field','fu_quantity_of_field')
      fu_quantity_of_field = int_missing
    end if
    
  END FUNCTION fu_quantity_of_field


  ! ***************************************************************

  FUNCTION fu_species_of_field (field) result(species)
    !
    ! Returns the substance name of a field-id.
    !
    IMPLICIT NONE
    
    ! Return value
    type(silam_species) :: species

    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_field), POINTER :: field
    
    if(defined(field))then
      species = fu_species(field%id)
    else
      call set_error('undefined field','fu_species_of_field')
      species = species_missing
    end if
    
  END FUNCTION fu_species_of_field


  ! ***************************************************************

  FUNCTION fu_substance_name_of_field (field) result(name)
    
    ! Description:
    ! Returns the substance name of a field-id.
    !
    IMPLICIT NONE
    
    ! Return value
    character(len=clen) :: name

    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_field), POINTER :: field
    
    if(defined(field))then
      name = fu_substance_name(field%id)
    else
      call set_error('undefined field','fu_substance_name_of_field')
      name = ''
    end if
    
  END FUNCTION fu_substance_name_of_field

  ! ***************************************************************

  FUNCTION fu_mode_of_field(field) result(mode)

    IMPLICIT NONE
    TYPE(silja_field), POINTER :: field
    type(Taerosol_mode)::mode

    if(defined(field))then
      mode = fu_mode(field%id)
    else
      call set_error('undefined field','fu_mode_of_field')
      mode = aerosol_mode_missing
    end if

  end FUNCTION fu_mode_of_field

  ! ***************************************************************

  FUNCTION fu_optic_wavelen_of_field(field)

    IMPLICIT NONE
    TYPE(silja_field), POINTER :: field
    real :: fu_optic_wavelen_of_field

    if(defined(field))then
      fu_optic_wavelen_of_field =  fu_optical_wave_length(field%id)
    else
      call set_error('undefined field','fu_optic_wavelen_of_field')
      fu_optic_wavelen_of_field = real_missing
    end if

  end FUNCTION fu_optic_wavelen_of_field

  ! ***************************************************************

  FUNCTION fu_field_id(field) result(id)
    
    ! Description:
    ! Returns the identification section of a field.
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_field_id), POINTER :: id
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_field), target :: field
    
    if(fu_field_defined(field))then
      id => field%id
    else
      call set_error('undefined field','fu_field_id')
      NULLIFY(id)
    end if
    
  END FUNCTION fu_field_id

  

  ! ***************************************************************

  INTEGER FUNCTION fu_fieldsize(field)
    
    ! Description:
    ! Returns the number of real numbers stored in field.
    ! 
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_field), POINTER :: field

    IF (defined(field)) THEN
      fu_fieldsize = SIZE(field%grid_data)
    else
      call set_error('undefined field','fu_fieldsize')
      fu_fieldsize = 0
    END IF
    
  END FUNCTION fu_fieldsize
  


  ! *****************************************************************

  LOGICAL FUNCTION fu_accumulated_field(field)

    ! Description:
    ! Returns a true value, if field is time-accumulated type.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silja_field), POINTER :: field

    IF (defined(field)) THEN
      fu_accumulated_field = fu_accumulated_id(field%id)
    else
      call set_error('undefined field','fu_accumulated_field')
      fu_accumulated_field = .false.
    END IF

  END FUNCTION fu_accumulated_field



  ! *****************************************************************

  LOGICAL FUNCTION fu_period_valid_field(field)

    ! Returns a true value, if field is valid through some time interval.
    !
    IMPLICIT NONE

    ! Imported parameters:
    TYPE(silja_field), POINTER :: field

    IF (defined(field)) THEN
      fu_period_valid_field = fu_validity_length(field%id) > zero_interval
    else
      call set_error('undefined field','fu_period_valid_field')
      fu_period_valid_field = .false.
    END IF

  END FUNCTION fu_period_valid_field


  ! ***************************************************************
  ! ***************************************************************
  !
  !
  !      Handle the values in fields.
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  subroutine update_field_data(newId, newGridData, originalField, iUpdateType)
    !
    ! Takes the new field Id and its data and sets the values in the original
    ! field to the values
    !
    implicit none

    ! Imported parameters
    type(silja_field_id), intent(in):: newId
    type(silja_field), pointer :: originalField
    real, dimension (:), pointer :: newGridData
    integer, intent(in) :: iUpdateType

    ! Local variables
    integer :: iCell,  n_out_of_range, n_patched, nx, ny, n_failed, fs_new
    TYPE(silja_level) :: levTmp
    character(len=*), parameter :: sub_name = 'update_field_data'
    !
    ! Basically the only thing, which we should check is that the sizes of the grids are
    ! the same. Who knows why we would like to exchange data. The quantities and other 
    ! parameters can be pretty different.
    !
    if(.not.defined(newId))then
      call set_error('Undefined new field id',sub_name)
      return
    endif
    if(.not.defined(originalField))then
      call set_error('Undefined original field',sub_name)
      return
    endif
    if(.not.associated(newGridData))then
      call set_error('Undefined new grid data',sub_name)
      return
    endif
    if(.not. fu_grid(originalField%id) == fu_grid(newId))then
      call report(fu_grid(originalField%id))
      call report(fu_grid(newId))
      call set_error('Grids are not equal', sub_name)
      return
    endif
      
    call grid_dimensions(fu_grid(newId), nx, ny)
    fs_new = nx * ny
    select case(iUpdateType)
      case(overwrite_field)
        
        originalField%grid_data(1:fs_new) = newGridData(1:fs_new)
      case(overwrite_field_anytime)
        ! FIXME
        !  The levels recovered from Silam NetCDF (as of now)
        ! come as surface-level ones. Initialisation of single-time fields created
        ! with other levels caused change of level, which is wnough  to break the
        ! consistency of updates and output
        ! Workaround: keep the original level
        ! Proper fix would be to have levels properly stored
        ! and a way to correct legacy files at the netcdf_io level
        levTmp = fu_level(originalField%id)
        if (fu_leveltype(levTmp) /= fu_leveltype(fu_level(newId))) then
          call msg("")
          call msg("originalField%id:")
          call report(originalField%id)
          call msg("newId:")
          call report(newId)
          call msg_warning("Attempt to overwrite with new levwltype, forcing old level", sub_name)
          originalField%id = newId        ! do not miss that !!
          call set_level(originalField%id,levTmp)
        else
          originalField%id = newId        ! do not miss that !!
        endif

        originalField%grid_data(1:fs_new)  = newGridData(1:fs_new)
      case(add_field)
        do iCell = 1, fs_new
          originalField%grid_data(iCell) = originalField%grid_data(iCell) + newGridData(iCell)
        enddo
      case default
        call msg('Strange field update requested: ', iUpdateType)
        call set_error('Strange update requested',sub_name)
        return
    end select

    call check_quantity_range(fu_quantity(newId), originalField%grid_data, fs_new, nx, &
                            & .false., .false., &  ! ifRequireValidity, ifSilent
                            & n_out_of_range, n_patched, n_failed)
    !n_negative = 0
    !IF ( .not. fu_negative_values_possible(fu_quantity(newId))) THEN  
    !  do iCell = 1, fu_number_of_gridpoints(fu_grid(newId))
    !    if(originalField%grid_data(iCell) < 0 .and. .not. (originalField%grid_data(iCell) .eps. real_missing))then
    !      n_negative = n_negative + 1
    !      originalField%grid_data(iCell) = 0.
    !    endif
    !  end do
    !END IF
    if (n_out_of_range > 0) then
      call msg('Out-of-range removed for the following field (update_field_data):' + &
                     & fu_quantity_string(fu_quantity(originalField%id)), n_out_of_range)
    end if
    
  end subroutine update_field_data


  ! ****************************************************************
  
  function fu_field_stats(field) result(stats)
    !
    ! Returns statistics of the field data
    !
    implicit none
    
    ! returned array
    real, dimension(3) :: stats
    
    ! imported parameter
    TYPE(silja_field), POINTER :: field
    
    if(defined(field))then
      stats(1:3) = (/MINVAL(field%grid_data), &
                  & (SUM(field%grid_data)/REAL(SIZE(field%grid_data))),&
                   & MAXVAL(field%grid_data)/)
    else
      stats(:) = real_missing
    endif
      
  end function fu_field_stats


  ! ****************************************************************
  ! 
  !
  !           INTERPOLATION
  !
  !
  ! ***************************************************************
  
  REAL FUNCTION fu_interpolate_scalar_to_pos(field, position, method) result(value)
    
    ! Description:
    ! Interpolates the field to given position.
    ! 
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silam_grid_position), INTENT(in) :: position
    INTEGER, INTENT(in) :: method ! of interpolation
    
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_field), POINTER :: field

    ! Local declarations:
    REAL, DIMENSION(:), POINTER :: grid_data_pointer
    INTEGER :: local_method
    ! ----------------------------------------------
    ! 
    ! 1. Check the input parameters.
    !    --------------------------
    
    IF (.NOT.defined(field)) THEN
      CALL set_error('undefined field given'&
	  & ,'fu_interpolate_scalar_to_pos')
      RETURN
    END IF

    IF (.NOT.defined(position)) THEN
      CALL set_error('undefined position given'&
	  & ,'fu_interpolate_scalar_to_pos')
      RETURN
    END IF
    
    ! -----------------------------------------------
    ! 
    ! 2. Interpolate.
    !    ------------

    grid_data_pointer => field%grid_data
    
    value = fu_2d_interpolation(grid_data_pointer, &
                                    & fu_grid(field), position, method)
    
    IF(error) THEN
      call msg(' % undefined value in fu_int_scalar_to_pos %')
      ! CALL print_fieldpointer_report(field) 
      local_method = nearest_point
      CALL unset_error('% Am unsetting errror and ')      
      call msg(' trying the nearest point...  %', value)
      value = fu_2d_interpolation(grid_data_pointer, &
                                      & fu_grid(field), position, local_method)
      
      IF(error) THEN
        call msg(' % Value still undefined , I give up !!!')
      END IF
    END IF

    
  END FUNCTION fu_interpolate_scalar_to_pos



  

  ! *****************************************************************
  ! ****************************************************************
  ! 
  !
  !        Tests
  !
  !
  ! ***************************************************************
  ! ***************************************************************
  
  SUBROUTINE print_fieldpointer_report(field)
    !
    ! Description:
    ! Prints a report of the field. All the identifications,
    ! status values and max/min values are printed.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_field), POINTER :: field
    !
    ! Local declarations:
    INTEGER :: i, n 
    character(len=worksize_string) :: str
    
    call msg(' %%%%%%%%%%%%%%%% Printing a field report %%%%%%%%%%%%%%%%')
    if(.not. associated(field))then
      call msg ('field not associated')
      RETURN
    END IF

    IF (.NOT.defined(field)) THEN
      call msg('Undefined field')
      RETURN
    END IF

    CALL report(field%id)

    if(.not. associated(field%grid_data))then
      call msg('grid data not associated')
      return
    endif

    call msg (' Field size: ',  SIZE(field%grid_data))

    WRITE (str, *) ' Field min, average and max: ', &
            & MINVAL(field%grid_data),&
            & (SUM(field%grid_data)/REAL(SIZE(field%grid_data))),&
            & MAXVAL(field%grid_data)
    call msg(str)

    n = 0
    DO i = 1, SIZE(field%grid_data)
      IF (field%grid_data(i) .eps. real_missing) n = n+1
    END DO

    if(fu_quantity(field%id) == concentration_flag .or. fu_quantity(field%id) == volume_mixing_ratio_flag)then
      DO i = 1, SIZE(field%grid_data)
        IF (field%grid_data(i) < 0) call msg('Negative value found: ', i, field%grid_data(i))
      END DO
    endif

    call msg(' number of missing values: ', n)

    call msg('%%%%%%%%% end of field report %%%%%%%%%%%%%%%%%%%%')
    

  END SUBROUTINE print_fieldpointer_report

END MODULE fields
