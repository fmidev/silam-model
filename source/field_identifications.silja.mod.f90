MODULE field_identifications

  ! Description: 
  ! Contains the definition of silja_field_id -type which defines
  ! fully one 2-D data on a level in a horizontal grid.
  ! Identification can defined either scalar-data, vector-data (wind)
  ! or material-field-data.
  !
  ! Original code: Mika Salonoja
  ! Author: Mikhail Sofiev, FMI, e-mail mikhail.sofiev@fmi.fi
  !
  ! All units: SI
  ! 
  ! Language: ANSI standard Fortran 90
  !
  ! 
  ! Modules used:
  
!  USE grids_geo
!  USE times 
  USE input_data_rules
  use chemical_setup
!  USE names_of_quantities  
    
  IMPLICIT NONE

  ! The public functions and subroutines available in this module:
  PUBLIC fu_set_field_id
  PUBLIC fu_set_field_id_simple
  PUBLIC fu_set_windfield_id
  PUBLIC fu_set_windfield_id_simple
  public set_missing
  PUBLIC fu_windfield_id_from_uvw_id
  PUBLIC set_quantity
  PUBLIC defined
  PUBLIC report
  public fu_field_id_covers_request
  PUBLIC fu_grid
  PUBLIC fu_u_grid
  PUBLIC fu_v_grid
  PUBLIC fu_w_grid
  PUBLIC fu_grids_of_id
  PUBLIC fu_met_src
  PUBLIC fu_level
  PUBLIC fu_valid_time
  public fu_field_valid
  PUBLIC fu_analysis_time
  PUBLIC fu_forecast_length
  PUBLIC fu_accumulated
  PUBLIC fu_accumulation_length
  PUBLIC fu_validity_length
  PUBLIC fu_accumulation_start_time
  PUBLIC fu_quantity
  public fu_species
  public fu_cocktail_name
  PUBLIC fu_substance_name
  public fu_mode
  public fu_optical_wave_length
  PUBLIC fu_field_kind
  PUBLIC fu_full_id
  public set_grid
  public set_field_kind
  public set_level
  public set_accumulation_length
  PUBLIC set_validity_length
  public set_analysis_time
  public set_valid_time
  public set_met_src
  PUBLIC set_species
  public set_cocktail_name
  public fu_if_internal_silam_field

  ! The private functions and subroutines not to be used elsewhere:
  private fu_set_field_id_from_params
  private fu_set_field_id_from_namelist
  private set_missing_field_id
  PRIVATE fu_compare_field_ids_eq ! interface ==
  PRIVATE fu_analysis_time_of_field_id
  PRIVATE fu_forecast_length_field_id
  PRIVATE fu_valid_time_of_field_id
  PRIVATE fu_level_field_id
  PRIVATE fu_quantity_of_field_id
  PRIVATE fu_chSubstNm_of_field_id
  PRIVATE fu_grid_of_field_id
  PRIVATE fu_met_src_of_field_id
  PRIVATE fu_field_id_defined
  PRIVATE print_field_id_report
  PRIVATE fu_accumulation_field_id
  private fu_validity_length_of_field_id
  PRIVATE fu_accum_start_time_id
  PRIVATE fu_accumulated_id
  private set_acc_len_of_field_id
  private set_valid_time_of_field_id
  private set_analysis_time_of_field_id
  private fu_species_of_field_id
  private fu_cocktail_name_of_field_id
  private fu_mode_of_field_id
  private fu_optical_wave_len_field_id
  private set_quantity_of_field_id
  private set_level_of_field_id
  private set_grid_of_field_id
  private set_met_src_of_field_id
  private fu_u_grid_of_windfield_id
  private fu_v_grid_of_windfield_id
  private set_valid_len_of_field_id
  private fu_w_grid_of_windfield_id
  private set_species_of_field_id
  !private set_species_from_basic_param

  ! Generic names and operator-interfaces of some functions:

  interface fu_set_field_id
    module procedure fu_set_field_id_from_params
    module procedure fu_set_field_id_from_namelist
  end interface

  interface set_missing
    module procedure set_missing_field_id
  end interface

  INTERFACE defined
    MODULE PROCEDURE fu_field_id_defined
  END INTERFACE

  INTERFACE report
    MODULE PROCEDURE print_field_id_report
  END INTERFACE

  INTERFACE operator(==)
    MODULE PROCEDURE fu_compare_field_ids_eq
  END INTERFACE

  INTERFACE fu_met_src  
    MODULE PROCEDURE fu_met_src_of_field_id
  END INTERFACE

  INTERFACE fu_grid 
    MODULE PROCEDURE fu_grid_of_field_id
  END INTERFACE

  INTERFACE fu_u_grid 
    MODULE PROCEDURE fu_u_grid_of_windfield_id
  END INTERFACE

  INTERFACE fu_v_grid 
    MODULE PROCEDURE fu_v_grid_of_windfield_id
  END INTERFACE

  INTERFACE fu_w_grid 
    MODULE PROCEDURE fu_w_grid_of_windfield_id
  END INTERFACE

  INTERFACE fu_analysis_time
    MODULE PROCEDURE fu_analysis_time_of_field_id
  END INTERFACE

  INTERFACE fu_valid_time
    MODULE PROCEDURE fu_valid_time_of_field_id
  END INTERFACE

  INTERFACE fu_forecast_length
    MODULE PROCEDURE fu_forecast_length_field_id
  END INTERFACE

  INTERFACE fu_accumulation_length
    MODULE PROCEDURE fu_accumulation_field_id
  END INTERFACE

  INTERFACE fu_validity_length
    module procedure fu_validity_length_of_field_id
  END INTERFACE

  INTERFACE fu_accumulation_start_time
    MODULE PROCEDURE fu_accum_start_time_id
  END INTERFACE

  INTERFACE fu_accumulated
    MODULE PROCEDURE fu_accumulated_id
  END INTERFACE

  INTERFACE fu_level
    MODULE PROCEDURE fu_level_field_id
  END INTERFACE

  INTERFACE fu_quantity
    MODULE PROCEDURE fu_quantity_of_field_id
  END INTERFACE

  interface fu_species
    module procedure fu_species_of_field_id
  end interface

  interface fu_cocktail_name
    module procedure fu_cocktail_name_of_field_id
  end interface

  INTERFACE fu_substance_name
    MODULE PROCEDURE fu_chSubstNm_of_field_id
  END INTERFACE

  INTERFACE fu_mode
    MODULE PROCEDURE fu_mode_of_field_id
  END INTERFACE

  INTERFACE fu_optical_wave_length
    MODULE PROCEDURE fu_optical_wave_len_field_id
  END INTERFACE

  interface set_grid
    module procedure set_grid_of_field_id
  end interface

  interface set_level
    module procedure set_level_of_field_id
  end interface

  interface set_accumulation_length
    module procedure set_acc_len_of_field_id
  end interface

  interface set_validity_length
    module procedure set_valid_len_of_field_id
  end interface

  interface set_field_kind
    module procedure set_field_kind_of_field_id
  end interface

  interface set_valid_time
    module procedure set_valid_time_of_field_id
  end interface

  interface set_analysis_time
    module procedure set_analysis_time_of_field_id
  end interface

  interface set_met_src
    module procedure set_met_src_of_field_id
  end interface

  interface set_quantity
    module procedure set_quantity_of_field_id
  end interface

  interface set_species
    module procedure set_species_of_field_id
  end interface

  interface set_cocktail_name
    module procedure set_cocktail_name_of_field_id
  end interface

  ! Possible kinds of the field. 
  ! NOTE. Kinds are as close as possible to WMO GRIB types
  ! included in Code Table 5. However, there are differences
  ! connected with strange definitions in that table

     ! instantaneous, valid at valid time analysis_time+forecast_length
  INTEGER, PUBLIC, PARAMETER :: forecast_flag = 3000

     ! Accumulated for length_of_accumulation before valid time
  INTEGER, PUBLIC, PARAMETER :: accumulated_flag = 3002

     ! same time period as above, but averaging
  INTEGER, PUBLIC, PARAMETER :: averaged_flag = 3003
  integer, public, dimension(3), parameter :: legal_field_kinds = &
       & (/forecast_flag,  accumulated_flag, averaged_flag/)


  ! Public types with private components defined in this module:

  TYPE silja_field_id ! identification of a field-data in a field
    PRIVATE 
    type(meteo_data_source) :: met_src !see weather_data_rules

    INTEGER :: quantity = int_missing! whose values are in the field (see globals)
      ! Obs. can also have value cocktail_flag, which means that field
      ! contains cocktails, or wind_flag, which means that field
      ! contains all available wind-components.

    type(silam_species) :: species           ! if one single species
    character(len=substNmLen) :: chCocktail  ! if a mixture of species

    TYPE(silja_time) :: analysis_time, valid_time
    TYPE(silja_interval) :: forecast_length
    INTEGER :: field_kind ! Analysis, averaged, accumulated, etc.

    TYPE(silja_interval) :: length_of_accumulation ! from valid time backwards.
                      ! If accmulation starts from analysis, then this and forecast length
                      ! are the same. Obs. this definition is dfferent from GRIB!

    type(silja_interval) :: length_of_validity ! Starting from valid_time
    TYPE(silja_grid), DIMENSION(3) :: grids ! three needed for wind
    TYPE(silja_level) :: level
    LOGICAL :: full
    TYPE(silja_logical) :: defined = silja_false
  END TYPE silja_field_id


  TYPE(silja_field_id), PARAMETER, public :: field_id_missing = &
                             & silja_field_id (met_src_missing, int_missing, &
                                             & species_missing, &
                                             & '', &
                                             & time_missing, time_missing, interval_missing,&
                                             & int_missing, &
                                             & interval_missing, interval_missing, &
                                             & (/grid_missing, grid_missing, grid_missing/) , &
                                             & level_missing, &
                                             & .false., &
                                             & silja_false)

CONTAINS 

  ! ***************************************************************

  FUNCTION fu_set_field_id_from_params(met_src,&
                                     & quantity, &
                                     & analysis_time,&
                                     & forecast_length, &
                                     & grid,&
                                     & level,&
                                     & length_of_accumulation, & ! optional
                                     & length_of_validity, &     !optional
                                     & field_kind, &             ! optional
                                     & species, &               ! optional, alternative to the previous
                                     & chCocktail) &            ! optional, if a mixture of species
                                     & result(id)
    ! 
    ! Sets the value for the identification sector desribing a field.
    ! The function is heavily modified to cover different types of fields
    ! However, requirements of the backward compatibility forced several
    ! clumsy places, first of all set/corrections/checking of the field_kind
    ! parameter.
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_field_id) :: id

    ! Imported parameters with intent(in):
    type(meteo_data_source), INTENT(in) :: met_src
    integer, intent(in) :: quantity
    TYPE(silja_time), INTENT(in) :: analysis_time
    TYPE(silja_interval), INTENT(in) :: forecast_length
    TYPE(silja_grid), INTENT(in) :: grid
    TYPE(silja_level), INTENT(in) :: level

    ! Optional imported parameters:
    TYPE(silja_interval), INTENT(in), OPTIONAL :: length_of_accumulation 
    ! from valid time backwards.
    ! If accmulation starts from analysis, then it == forecast length
    ! Obs. this definition is dfferent from GRIB.

    TYPE(silja_interval), INTENT(in), OPTIONAL :: length_of_validity 
    ! For period-valid fields. If presents then the field validity 
    ! goes from valid_time over this length_of_validity

    integer, intent(in), optional :: field_kind ! SILAM type as defined above
       ! optional only for the backward compatibility. HIGHLY DESIRABLE.

!    character(len=*), intent(in), optional :: chSubstNm ! For chemical/radioactive
!       ! species, where this determines the substance in the cocktail
!
!    type(Taerosol_mode), intent(in), optional :: aerMode ! aerosol mode size, optical wave length
!    real, intent(in), optional :: fWaveLen   ! aerosol mode size, optical wave length
    type(silam_species), intent(in), optional :: species
    character(len=*), intent(in), optional :: chCocktail ! If a mixture of species
    
    !----------------------------------------
    !
    ! 1. Set general values
    !
    id%met_src = met_src

    id%quantity = quantity 
    !
    ! Chemical material can be absent, can be a single species or a cocktail of species
    !
    if(present(species))then
      id%species = species
    else
      id%species = species_missing
    endif
    if(present(chCocktail))then
      id%chCocktail = chCocktail
    else
      id%chCocktail = ''
    endif

    if(defined(id%species) .and. len_trim(id%chCocktail) > 0) then
      call set_error('Inconsistent arguments: crossed species and cocktail', &
                   & 'fu_set_field_id_from_params')
      return
    end if

    !----------------------------------------
    !
    ! 2. Set times.
    !
    IF (defined(analysis_time)) THEN
      id%analysis_time = analysis_time
    ELSE
      CALL set_error('analysis time not defined', 'fu_set_field_id_from_params')
      RETURN
    END IF

    IF (defined(forecast_length)) THEN 
      id%forecast_length = forecast_length
    ELSE
      CALL set_error('forecast length not defined', 'fu_set_field_id_from_params')
      RETURN
    END IF

    id%valid_time = analysis_time + forecast_length
    IF (error) RETURN

    !----------------------------------------
    !
    ! If the field is non-trivial (period-lasting, averaged, difference)
    ! this has to be stated explicitly via the field_kind
    ! Otherwise we do not know how to treat the length_of_validity,
    ! which starts to play a role. Also a role of the valid_time
    ! may change.
    ! The same is true for the length_of_accumulation. All this stuff has 
    ! to be checked
    !

    if(present(field_kind))then

      id%field_kind = field_kind

      select case (field_kind)
        case (forecast_flag) 
          !
          ! No accumulation and validity periods. Field_kind is already set
          !
          id%field_kind = forecast_flag

          if(present(length_of_accumulation))then
            if(defined(length_of_accumulation))then
              if(.not. (length_of_accumulation == zero_interval))then
                call set_error('Forecast/analysis field must have zero accumulation length', &
                             & 'fu_set_field_id_from_params')
                return
              endif
            endif
          endif
          id%length_of_accumulation = zero_interval
!          id%length_of_validity = zero_interval

        case (accumulated_flag, averaged_flag)
          !
          ! Accumulation must be non-zero
          !
!          if(present(length_of_validity))then
!            if(defined(length_of_validity))then
!              if(.not.(length_of_validity == zero_interval))then
!                call set_error('Accumulated field must have zero validity length', &
!                             & 'fu_set_field_id_from_params')
!                return
!              endif
!            endif
!          endif
!          id%length_of_validity = zero_interval
          IF (PRESENT(length_of_accumulation)) THEN
            IF (defined(length_of_accumulation)) THEN
              id%length_of_accumulation = length_of_accumulation
              IF (error) RETURN
            ELSE
              call set_error('Undefined accumulation of accumulated field', &
                           & 'fu_set_field_id_from_params')
              return
            END IF
          else
            call set_error('Missing accumulation of accumulated field', &
                         & 'fu_set_field_id_from_params')
            return
          endif


        case default 
          call set_error('Non-supported field type','fu_set_field_id_from_params')
          return
      end select

    else  ! present(field_kind)
       id%field_kind = forecast_flag
      
    end if ! present(field_kind)

    !----------------------------------------
    !
    ! Whatever the field-kind is, it can be valid during more than 1 nanosecond
    !
    IF (PRESENT(length_of_validity)) THEN
      IF (defined(length_of_validity)) THEN
        id%length_of_validity = length_of_validity
        IF (error) RETURN
      else
        id%length_of_validity = zero_interval
      endif
    else
      id%length_of_validity = zero_interval
    endif

    !----------------------------------------
    !
    ! Set grid and level.
    !
    IF (defined(grid)) THEN
      id%grids(1) = grid
      id%grids(2) = grid_missing 
      id%grids(3) = grid_missing 
    ELSE
      CALL set_error('grid not defined','fu_set_field_id_from_params')
      RETURN
    END IF

    IF (defined(level)) THEN
      id%level = level
    ELSE
      CALL set_error('LEVEL not defined', 'fu_set_field_id_from_params')
      RETURN
    END IF

    id%full = .true.
    id%defined = fu_set_true()

  END FUNCTION fu_set_field_id_from_params


  !*********************************************************************

  function fu_set_field_id_from_namelist(nlId)result(id)
    !
    ! Sets a full-size field id from the namelist
    !
    implicit none

    ! Return value of the function
    type(silja_field_id) :: id

    ! Imported parameters
    type(Tsilam_namelist), pointer :: nlId

    ! Local variables
    type(silja_level), dimension(:), pointer :: levels
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: pItems
    integer :: nItems, year, month, day, hour, min, iStatus
    real :: sec, dens, min_d, max_d, mean_d
    type(silam_sp) :: spTmp

    call set_missing(id)    ! If something goes wrong
    spTmp%sp => fu_work_string()

    !----------------------------------------
    !
    ! Set general values
    !

    id%met_src = fu_set_met_src(nlId) ! Can be undefined

    id%quantity = fu_get_silam_quantity(fu_content(nlId, 'quantity_name'))
    if(error .or. id%quantity == int_missing)then
      call set_error('Undefined quantity','fu_set_field_id_from_namelist')
      call free_work_array(spTmp%sp)
      return  ! Must be defined
    endif
    !
    ! Chemical material can be absent, can be a single species or a cocktail of species
    !
    spTmp%sp = fu_content(nlId,'substance_name') ! Can be empty
    if(len_trim(spTmp%sp) > 0)then  ! substance given
    
      call set_species(id%species, fu_get_material_ptr(fu_content(nlId,'substance_name')), &
                       & fu_set_mode(nlId))
    
    else
      id%species = species_missing
      id%chCocktail = fu_content(nlId,'cocktail_name') ! Can be empty
      
    endif
    if(error)return

    !
    ! Set grid and level.
    !
    id%grids(1) = fu_set_grid(nlId)
    if(.not. defined(id%grids(1)))then
      call free_work_array(spTmp%sp)
      return  ! Must be defined
    endif

    id%grids(2) = grid_missing 
    id%grids(3) = grid_missing 

    nullify(levels)
    call create_levels_from_namelist_v2(levels, nlId)
    if(error)return
    
    if(size(levels) /= 1)then
      call msg('Strange number of levels in namelist',size(levels))
      call set_error('Strange number of levels in namelist','fu_set_field_id_from_namelist')
      call free_work_array(spTmp%sp)
      return
    endif

    id%level = levels(1)
    deallocate(levels)
    if(.not. defined(id%level))then
      call free_work_array(spTmp%sp)
      return  ! Must be defined
    endif

    !----------------------------------------
    !
    ! Set times.
    !
    nullify(pItems)
    call get_items(nlId, 'analysis_time', pItems, nItems)
    if(nItems /= 1)then
      call msg('Strange number of analysis times',nItems)
      call set_error('Strange number of analysis times','fu_set_field_id_from_namelist')
      return
    endif
    spTmp%sp = fu_str_u_case(fu_content(pItems(1)))

    if(spTmp%sp == 'UNDEFINED_TIME')then
      id%analysis_time = time_missing
      id%forecast_length = interval_missing
      id%valid_time = time_missing
    else
      if(index(spTmp%sp,'UTC') == 0)then
         call msg('No UTC indicator. Assuming UTC time. "utc_difference" is not used anymore')
      endif
      if(error)return

      read(unit=spTmp%sp, iostat=iStatus, fmt=*) year, month, day, hour, min, sec
      if(iStatus /= 0)then
        call set_error(fu_connect_strings('Failed to get analysis time from:',spTmp%sp), &
                     & 'fu_set_field_id_from_namelist')
        return
      endif
      id%analysis_time = fu_set_time_utc(year, month, day, hour, min, sec)
      if(.not.defined(id%analysis_time)) return
      !
      ! Set forecast length
      !
      id%forecast_length = fu_set_named_interval(fu_content(nlId,'forecast_length'))
      if(defined(id%forecast_length))then
        id%valid_time = id%analysis_time + id%forecast_length
      else
        id%valid_time = time_missing
      endif
    endif
    IF (error) RETURN
    !
    ! Validity length may be zero or non-zero for any type of the field
    !
    IF (fu_content(nlId,'length_of_validity') /= '') THEN
      id%length_of_validity = fu_set_named_interval(fu_content(nlId, &
                                                             & 'length_of_validity'))
      IF (.not.defined(id%length_of_validity)) THEN
        call set_error('Undefined length of validity','fu_set_field_id_from_namelist')
        return
      END IF
    else
      id%length_of_validity = zero_interval
    endif

    !----------------------------------------
    !
    ! If the field is non-trivial (period-lasting, averaged, difference)
    ! this has to be stated explicitly via the field_kind
    ! Otherwise we do not know how to treat the length_of_validity,
    ! which starts to play a role. Also a role of the valid_time
    ! may change.
    ! The same is true for the length_of_accumulation. All this stuff has 
    ! to be checked
    !
    if(fu_content(nlId,'field_kind') /= '')then

      select case (fu_content(nlId,'field_kind'))
        case ('forecast_field', 'analysis_field') !!backward compatibility
          !
          ! No accumulation period, ignored even if present
          !
          id%field_kind = forecast_flag
          id%length_of_accumulation = zero_interval

        case ('accumulated_field', 'averaged_field')
          !
          ! Accumulation must exist
          !
          IF (fu_content(nlId,'length_of_accumulation') /= '') THEN
            id%length_of_accumulation = fu_set_named_interval(fu_content(nlId, &
                                                                 & 'length_of_accumulation'))
            IF (.not. defined(id%length_of_accumulation)) THEN
              call set_error('Undefined accumulation of accumulated field', &
                           & 'fu_set_field_id_from_namelist')
              return
            END IF
          else
            call set_error('Missing accumulation of accumulated field', &
                         & 'fu_set_field_id_from_namelist')
            return
          endif

        case default 
          call set_error('Non-supported field type','fu_set_field_id_from_namelist')
          return
      end select

    else  ! present non-empty field_kind
      !
      ! If field_kind is not present, some default is obtained from the forecast length
      !
      id%field_kind = forecast_flag
      
    end if ! present non-empty field_kind

    id%full = .true.
    id%defined = fu_set_true()

    call free_work_array(spTmp%sp)

  end function fu_set_field_id_from_namelist


  ! ***************************************************************

  FUNCTION fu_set_field_id_simple(met_src, quantity, valid_time, level, species, chCockt) result(id)
    !
    ! Sets the value for the identification sector desribing a
    ! required field, and parts of the identification may be undefined.
    !
    ! If level is undefined, it is assumed that the quantity itself
    ! defines its level, and level is not required for fully defining
    ! the field while retrieving it.
    !
    ! If time is undefined, it is assumed that the quantity is time-indepemdent
    ! (for example topography etc.) and thus time
    ! is not required for fully defining the field while retrieving it.
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_field_id) :: id

    ! Imported parameters with intent(in):
    type(meteo_data_source), INTENT(in) :: met_src
    integer, intent(in) :: quantity
    TYPE(silja_time), INTENT(in) :: valid_time
    TYPE(silja_level), INTENT(in) :: level
    type(silam_species),intent(in), optional :: species
    character(len=*), intent(in), optional :: chCockt

    id%met_src = met_src
    id%quantity = quantity 
    if(present(species))then
      id%species = species
    else
      id%species = species_missing
    endif
    if(present(chCockt))then
      id%chCocktail = chCockt
    else
      id%chCocktail = ''
    endif
    id%valid_time = valid_time
    id%length_of_validity = zero_interval
    id%analysis_time = time_missing
    id%forecast_length = interval_missing
    id%grids = grid_missing 
    id%level = level
    id%full = .false.
    id%field_kind = forecast_flag
    id%defined = fu_set_true()

  END FUNCTION fu_set_field_id_simple


  !**************************************************************************

  subroutine set_missing_field_id(id)
    !
    ! Sets the given field id as missing
    !
    implicit none

    type(silja_field_id), intent(out) :: id

    id%met_src = met_src_missing
    id%quantity = int_missing
    id%species = species_missing
    id%chCocktail = ''
    id%analysis_time = time_missing
    id%valid_time = time_missing
    id%forecast_length = interval_missing
    id%field_kind = int_missing
    id%length_of_accumulation = interval_missing
    id%length_of_validity = interval_missing
    id%grids(1) = grid_missing
    id%grids(2) = grid_missing
    id%grids(3) = grid_missing
    id%level = level_missing
    id%full = .false.
    id%defined = silja_false

  end subroutine set_missing_field_id


  ! ***************************************************************

  SUBROUTINE set_quantity_of_field_id(id, quantity, chSubstNm, chCocktail)
  !
  ! Sets a quantity to a requested flag
  !
  ! Code owner Mikhail Sofiev, FMI
  !
  IMPLICIT NONE

  ! Imported parameters, intent IN
  INTEGER, INTENT(in) :: quantity
  character(len=*), intent(in), optional :: chSubstNm, chCocktail

  ! Both in and out parameter
  TYPE(silja_field_id), INTENT(inout) :: id
    
    IF(.not.fu_known_quantity(quantity))THEN
      CALL set_error('Unknown quantity','fu_set_quantity')
      RETURN
    END IF

    id%quantity = quantity
    if(present(chSubstNm))then
      if(present(chCocktail))then
        call set_error('Substance name and cocktail are both present','set_quantity_of_field_id')
        return
      endif
      call set_species(id%species, fu_get_material_ptr(chSubstNm), aerosol_mode_missing)
      id%chCocktail = ''
    else
      id%species = species_missing
      if(present(chCocktail)) id%chCocktail = chCocktail
    endif
  
  END SUBROUTINE set_quantity_of_field_id



  ! ***************************************************************

  FUNCTION fu_set_windfield_id(quantity, &
                             & met_src,&
                             & kind, &
                             & analysis_time,&
                             & forecast_length, &
                             & u_grid,&
                             & v_grid,&
                             & w_grid,&
                             & level,&
                             & length_of_accumulation) result(id)

    ! Description:
    ! Sets the value for the identification sector desribing a field.
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_field_id) :: id

    ! Imported parameters with intent(in):
    type(meteo_data_source), INTENT(in) :: met_src
    integer, intent(in) :: kind, quantity
    TYPE(silja_time), INTENT(in) :: analysis_time
    TYPE(silja_interval), INTENT(in) :: forecast_length
    TYPE(silja_grid), INTENT(in) :: u_grid
    TYPE(silja_grid), INTENT(in) :: v_grid
    TYPE(silja_grid), INTENT(in) :: w_grid
    TYPE(silja_level), INTENT(in) :: level

    ! Optional imported parameters:
    TYPE(silja_interval), INTENT(in), OPTIONAL :: length_of_accumulation

    !----------------------------------------
    !
    ! 1. Set general values
    !    ------------------

    id%met_src = met_src

    if(quantity == u_10m_flag .or. quantity == v_10m_flag) THEN
      id%quantity = wind_10m_flag
    ELSE
      id%quantity = wind_flag
    END IF
    id%species = species_missing
    id%chCocktail = ''

    !----------------------------------------
    !
    ! 2. Set times.
    !    ----------

    id%field_kind = kind

    IF (defined(analysis_time)) THEN
      id%analysis_time = analysis_time
    ELSE
      CALL set_error('analysis time not defined','fu_set_windfield_id')
      RETURN
    END IF

    IF (defined(forecast_length)) THEN 
      id%forecast_length = forecast_length
    ELSE
      CALL set_error('forecast length not defined','fu_set_windfield_id')
      RETURN
    END IF

    IF (PRESENT(length_of_accumulation)) THEN
      IF (defined(length_of_accumulation)) THEN
   id%length_of_accumulation = length_of_accumulation
      ELSE
   CALL set_error('length_of_accumulation not defined',&
       & 'fu_set_windfield_id')
   RETURN
      END IF
    END IF

    id%valid_time = analysis_time + forecast_length


    !----------------------------------------
    !
    ! 3. Set grid and level.
    !    -------------------

    IF (defined(u_grid).and.defined(v_grid)) THEN
      id%grids(1) = u_grid
      id%grids(2) = v_grid
    ELSE
      CALL set_error('uv grids not defined','fu_set_windfield_id')
      RETURN
    END IF

    IF (defined(w_grid)) THEN
      id%grids(3) = w_grid
    ELSE
      id%grids(3) = grid_missing
    END IF

    IF (defined(level)) THEN
      id%level = level
    ELSE
      CALL set_error('LEVEL not defined', 'fu_set_windfield_id')
      RETURN
    END IF


    !----------------------------------------
    !
    ! 4. Everything ok.
    !    --------------

    id%full = .true.
    id%defined = fu_set_true()


  END FUNCTION fu_set_windfield_id


  ! ***************************************************************

  FUNCTION fu_set_windfield_id_simple(quantity, met_src, valid_time, level) result(id)

    ! Description:
    ! Sets the value for the identification sector desribing a
    ! required field, do id may be partially undefined.
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_field_id) :: id

    ! Imported parameters:
    integer, intent(in) :: quantity
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), INTENT(in) :: valid_time
    TYPE(silja_level), INTENT(in) :: level

    ! Local declarations:
    ! 
    !----------------------------------------
    !
    ! 1. Set general values
    !    ------------------

    id%met_src = met_src

    if(quantity == u_10m_flag .or. quantity == v_10m_flag) THEN
      id%quantity = wind_10m_flag
    ELSE 
      id%quantity = wind_flag
    END IF
    id%species = species_missing
    id%chCocktail = ''

    !----------------------------------------
    !
    ! 2. Set times.
    !    ----------

    IF (defined(valid_time)) THEN
      id%valid_time = valid_time
    ELSE
      CALL set_error('valid time not defined',&
     & 'fu_set_windfield_id_simple')
      RETURN
    END IF

    id%analysis_time = time_missing
    id%forecast_length = interval_missing

    !----------------------------------------
    !
    ! 3. Set grid and level.
    !    -------------------

    id%grids = grid_missing
    id%level = level


    !----------------------------------------
    !
    ! 4. Everything ok.
    !    --------------

    id%full = .false.
    id%defined = fu_set_true()

  END FUNCTION fu_set_windfield_id_simple


  ! ***************************************************************

  FUNCTION fu_windfield_id_from_uvw_id(u_id, v_id, w_id) result(id)

    ! Description:
    ! Creates one windfield-id from corresponding u,v and w
    ! identifications. W-id is optional. No checkings are done, so be
    ! careful. Times etc. are defined by u-field. 
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_field_id) :: id
    !
    ! Imported parameters with intent(in):
    TYPE(silja_field_id), INTENT(in) :: u_id, v_id

    ! Optional parameters with intent(in):
    TYPE(silja_field_id), INTENT(in), OPTIONAL :: w_id

    IF (PRESENT(w_id)) THEN
      id = fu_set_windfield_id(u_id%quantity, &
     & u_id%met_src ,&
      & u_id%field_kind, &
     & u_id%analysis_time ,&
     & u_id%forecast_length ,&
     & u_id%grids(1) ,&
     & v_id%grids(1) ,&
     & w_id%grids(1) ,&
     & u_id%level)

    ELSE
      id = fu_set_windfield_id(u_id%quantity, &
     & u_id%met_src ,&
      & u_id%field_kind, &
     & u_id%analysis_time ,&
     & u_id%forecast_length ,&
     & u_id%grids(1) ,&
     & v_id%grids(1) ,&
     & grid_missing ,&
     & u_id%level)
    END IF

  END FUNCTION fu_windfield_id_from_uvw_id

  
  
!!!$  ! ***************************************************************
!!!$  
!!!$  SUBROUTINE set_w_grid_to_id(id, w_grid)
!!!$    
!!!$    ! Description:
!!!$    ! Replaces the grid in id.
!!!$    ! 
!!!$    ! Language: ANSI Fortran 90
!!!$    !
!!!$    ! Author: Mika Salonoja, FMI
!!!$    ! 
!!!$    IMPLICIT NONE
!!!$
!!!$    ! Imported parameters with intent IN:
!!!$    TYPE(silja_grid), INTENT(in) :: w_grid
!!!$
!!!$    ! Imported parameters with intent INout:
!!!$    TYPE(silja_field_id), INTENT(inout) :: id
!!!$
!!!$    id%grids(3) = w_grid
!!!$
!!!$  END SUBROUTINE set_w_grid_to_id
!!!$


  ! ***************************************************************

  FUNCTION fu_grids_of_id(field_id) result(grids) 

    ! Returns the grid of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_grid), DIMENSION(3) :: grids

    ! Imported parameters with intent(in):
    TYPE(silja_field_id), INTENT(in) :: field_id

    grids = field_id%grids

  END FUNCTION fu_grids_of_id


  ! ***************************************************************

  LOGICAL FUNCTION fu_full_id(id)

    ! Returns the grid of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_field_id), INTENT(in) :: id

    IF (defined(id)) THEN
      fu_full_id = id%full
    ELSE
      fu_full_id = .false.
    END IF

  END FUNCTION fu_full_id





  ! ***************************************************************
  ! ***************************************************************
  !
  !
  !      Private functions and subroutines.
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  LOGICAL FUNCTION fu_compare_field_ids_eq(id1, id2) result(eq)

    ! Description:
    ! Returns a true value if identification sections
    ! are  the same. The comparision is more rigorous for full id:s.
    ! For simple id:s both time and level can be undefined.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent(in):
    TYPE(silja_field_id), INTENT(in) :: id1, id2

    ! --------------------------------------------------------
    !
    ! 1. Check id:s.
    !
    IF (.NOT.defined(id1)) THEN
      eq = .NOT.defined(id2)
      RETURN
    END IF

    IF (.NOT.defined(id2)) THEN
      eq = .NOT.defined(id1)
      RETURN
    END IF

    !
    ! 2. Compare full id:s
    !
    IF (id1%full.and.id2%full) THEN
      eq = (  ((id1%met_src == id2%met_src) .or. &
             & (id1%met_src == met_src_missing) .or. &
             & (id2%met_src == met_src_missing)) .and. &
          & (id1%quantity == id2%quantity) .and. &
          & (id1%valid_time == id2%valid_time) .and. &
          & fu_cmp_levs_eq(id1%level, id2%level))

      ! Crap with species vs cocktails: only one can be defined. 
      if(len_trim(id1%chCocktail) > 0) then
        eq = eq .and. (trim(id1%chCocktail) == trim(id2%chCocktail))
      else
        eq = eq .and. (id1%species == id2%species)
      endif
      RETURN
    END IF

    !
    ! 3. Compare in case at least one is a simple id.
    !
    eq = .false.

    if(.not.(id1%met_src == met_src_missing) .and. .not.(id2%met_src == met_src_missing))then
      IF (.not.id1%met_src == id2%met_src) RETURN
    end if

    IF (id1%quantity /= id2%quantity) RETURN

    if(defined(id1%species) .or. defined(id2%species))then
      IF (.not.id1%species == id2%species) RETURN
    end if

    if(id1%chCocktail /= id2%chCocktail)return

    IF (defined(id1%valid_time).and. defined(id2%valid_time)) THEN
      IF (.NOT.(id1%valid_time == id2%valid_time)) RETURN
    END IF

    IF (defined(id1%level).and.defined(id2%level)) THEN
      IF (.NOT.(fu_cmp_levs_eq(id1%level, id2%level))) RETURN
    END IF

    eq = .true.     

  END FUNCTION fu_compare_field_ids_eq


  !******************************************************************************

  LOGICAL FUNCTION fu_field_id_covers_request(idIn, idRequest, ifDisregardGrid) result(accept_field)

    ! Returns a true value, if the given field-idIn "covers" the requested id, i.e.
    ! this field can be provided by someone who requested the idRequested
    !
    ! The horizontal grid is not cheked here since reprojections are possible - only coverage
    !
    IMPLICIT NONE

    ! Imported parameters with intent(in):
    TYPE(silja_field_id), INTENT(in) :: idIn, idRequest
    logical, intent(in) :: ifDisregardGrid

    ! Local declarations:
    INTEGER :: i, quantity_in_field
    INTEGER :: lower_type, upper_type, level_type
    LOGICAL :: time_accept, quantity_accept
    LOGICAL :: lower_accept, upper_accept
    TYPE(silja_level) :: level
    type(silja_grid) :: gridTmp

    !
    ! Stupidity checking first.
    !
    accept_field = .false.

    IF (.NOT.defined(idIn)) THEN
      CALL set_error('field-idIn not defined','fu_field_id_covers_request')
      RETURN
    END IF

    IF (.NOT. (idRequest%defined == silja_true)) THEN
      CALL set_error('Requested id is not defined','fu_field_id_covers_request')
      RETURN
    END IF


    !
    ! Time checking: the idIn must cover the validity range of the request
    !
    if(defined(idRequest%valid_time))then
      if (defined(idIn%valid_time)) then
        if(defined(idRequest%length_of_validity))then
          if(idIn%valid_time > idRequest%valid_time .or. &
           & idIn%valid_time + idIn%length_of_validity < &
           & idRequest%valid_time + idRequest%length_of_validity)return
        else
          if(.not. fu_between_times(idRequest%valid_time, idIn%valid_time, &
                                  & idIn%valid_Time+idIn%length_of_validity, &
                                  & .true.))return
        endif
      else
        ! time specified in request, but missing from the field....
        return
      endif
    endif

    !
    ! Check source.
    !
    IF (.not. idRequest%met_src == met_src_missing) THEN
      IF (.not. idRequest%met_src == idIn%met_src) RETURN
    END IF

    !
    ! 4. Check quantity. Only known quantities are accepted.
    !
    if(fu_known_quantity(idRequest%quantity))then
      if(.not.fu_known_quantity(idIn%quantity))return
      if(idIn%quantity /= idRequest%quantity)return
    endif

    !----------------------------------------
    !
    ! 5. Check grid. Grids should not nesessarily be equal - it is
    !    enough if they Arakawa-correspond to each other and the id's grid
    !    covers the list's one. Arakawa-correspondance means that they are
    !    equal, except for dimensions nx, ny and possible shift up to one 
    !    element in x and/or y dimension.
    !    Another extension: the input grid might be non-standard for SILAM.
    !    In particular, it can be global with (0:360) coverage instead of -180,180.
    !    So far, a clumsy way to handle this is to use temporary standard grid,
    !    analogous to the input one. Later, the grid adjustment will be done
    !    if the field is accepted and stored into the stack.
    !
    if(.not. ifDisregardGrid)then
      if(.not.(idRequest%grids(1) == grid_missing))then
        !
        ! Get the grid and reposition it if needed
        !
        gridTmp = idIn%grids(1)
        if(.not.fu_stdSilamGrid(gridTmp))then
          if(fu_ifLonGlobal(gridTmp))then
            if(.not.fu_stdSilamGrid(gridTmp)) call reposition_global_grid_lon(gridTmp)
          endif
          if(fu_ifLatGlobal(gridTmp))then
            if(.not.fu_stdSilamGrid(gridTmp)) call reposition_global_grid_lat(gridTmp)
          endif
          if(.not.fu_stdSilamGrid(gridTmp))return
        endif  ! input grid is non-standard

        if(.not.(fu_grids_arakawa_correspond(idRequest%grids(1),gridTmp)))then
          if(.not.(fu_if_grid_covered(idRequest%grids(1),gridTmp)))return  ! failed to cover the requested grid
        endif
      end if  ! if grid checking is needed
    end if  ! if grid checking is needed

    !
    ! Check levels.
    !
    if(defined(idRequest%level))then
      IF (.NOT. fu_cmp_levs_eq(idIn%level, idRequest%level)) RETURN
    endif

    ! Check cocktail & species. Only one is allowed to be defined.
    ! 
    if (defined(fu_species(idIn)) .and. fu_cocktail_name(idIn) /= '') then
      call set_error('idIn has both species and cocktail defined', 'fu_field_id_covers_request')
    end if
    if (defined(fu_species(idRequest)) .and. fu_cocktail_name(idRequest) /= '') then
      call set_error('idRequest has both species and cocktail defined', 'fu_field_id_covers_request')
    end if
    if (error) return

    if (defined(fu_species(idRequest))) then
      if (.not. defined(fu_species(idIn))) return
      if (.not.(fu_species(idRequest) == fu_species(idIn))) return
    end if
    if (fu_cocktail_name(idRequest) /= '') then
      if (fu_cocktail_name(idRequest) /= fu_cocktail_name(idIn)) return
    end if

    ! Accept this field-idIn
    !
    accept_field = .true.

  END FUNCTION fu_field_id_covers_request


  ! ***************************************************************

  FUNCTION fu_grid_of_field_id(field_id) result(grid) 

    ! Returns the grid of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_grid) :: grid
    !
    ! Imported parameters with intent(in):
    TYPE(silja_field_id), INTENT(in) :: field_id

    IF(field_id%quantity == wind_flag) THEN
      CALL msg_warning('asked single grid from windfield','fu_grid_of_field_id')
    END IF

    grid = field_id%grids(1)

  END FUNCTION fu_grid_of_field_id



  ! ***************************************************************

  FUNCTION fu_u_grid_of_windfield_id(id)

    ! Returns the u-grid of a field-id when field contains wind.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_grid) :: fu_u_grid_of_windfield_id
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_field_id), INTENT(in) :: id

    fu_u_grid_of_windfield_id = id%grids(1)

  END FUNCTION fu_u_grid_of_windfield_id



  ! ***************************************************************

  FUNCTION fu_v_grid_of_windfield_id(id)

    ! Returns the u-grid of a field-id when field contains wind.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_grid) :: fu_v_grid_of_windfield_id
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_field_id), INTENT(in) :: id

    fu_v_grid_of_windfield_id = id%grids(2)

  END FUNCTION fu_v_grid_of_windfield_id



  ! ***************************************************************

  FUNCTION fu_w_grid_of_windfield_id(id)

    ! Returns the u-grid of a field-id when field contains wind.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_grid) :: fu_w_grid_of_windfield_id
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_field_id), INTENT(in) :: id

    fu_w_grid_of_windfield_id = id%grids(3)

  END FUNCTION fu_w_grid_of_windfield_id



  ! ***************************************************************

  FUNCTION fu_met_src_of_field_id(field_id) result(met_src) 

    ! Returns the data source of a field-id.
    !
    IMPLICIT NONE
    !
    type(meteo_data_source) :: met_src
    ! Imported parameters with intent(in):
    TYPE(silja_field_id), INTENT(in) :: field_id

    met_src = field_id%met_src

  END FUNCTION fu_met_src_of_field_id



  ! ***************************************************************

  FUNCTION fu_level_field_id(field_id) result(level)
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
    TYPE(silja_level) :: level
    !
    ! Imported parameters with intent(in):
    TYPE(silja_field_id), INTENT(in) :: field_id

    level = field_id%level

  END FUNCTION fu_level_field_id



  ! ***************************************************************

  FUNCTION fu_analysis_time_of_field_id(field_id)

    ! Description:
    ! Returns the analysis time of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_time) :: fu_analysis_time_of_field_id
    !
    ! Imported parameters with intent(in):
    TYPE(silja_field_id), INTENT(in) :: field_id

    fu_analysis_time_of_field_id = field_id%analysis_time

  END FUNCTION fu_analysis_time_of_field_id



  ! ***************************************************************

  FUNCTION fu_valid_time_of_field_id(field_id)

    ! Description:
    ! Returns the valid time of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_time) :: fu_valid_time_of_field_id
    !
    ! Imported parameters with intent(in):
    TYPE(silja_field_id), INTENT(in) :: field_id

    fu_valid_time_of_field_id = field_id%valid_time

  END FUNCTION fu_valid_time_of_field_id


  !****************************************************************

  logical function fu_field_valid(field_id, valid_time)
    !
    ! Checks of the field validity period covers the given valid time
    !
    implicit none

    ! Imported parameters with intent(in):
    TYPE(silja_field_id), INTENT(in) :: field_id
    type(silja_time), intent(in) :: valid_time
    
    fu_field_valid = fu_between_times(valid_time, &            ! time to check
                                    & field_id%valid_time, &   ! lim 1
                                    & field_id%valid_time+field_id%length_of_validity, & ! lim 2
                                    & .true.)                  ! if accept boundaries

  end function fu_field_valid
  

  ! ***************************************************************

  FUNCTION fu_accum_start_time_id(field_id)

    ! Description:
    ! Returns the time from which accumulation starts (end time
    ! of accumulation = valid time).
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_time) :: fu_accum_start_time_id

    ! Imported parameters with intent(in):
    TYPE(silja_field_id), INTENT(in) :: field_id

    fu_accum_start_time_id= field_id%valid_time-field_id%length_of_accumulation

  END FUNCTION fu_accum_start_time_id


  ! ***************************************************************

  FUNCTION fu_forecast_length_field_id(field_id)

    ! Description:
    ! Returns the analysis time of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_interval) :: fu_forecast_length_field_id
    !
    ! Imported parameters with intent(in):
    TYPE(silja_field_id), INTENT(in) :: field_id

    fu_forecast_length_field_id = field_id%forecast_length

  END FUNCTION fu_forecast_length_field_id



  ! ***************************************************************

  FUNCTION fu_accumulation_field_id(field_id)

    ! Description:
    ! Returns the of accumulation of a field-id. That is
    ! the interval from valid time backward that is took the
    ! quantity to valid time. If accumulation starts from analysis
    ! (as usually is) then accumulation and forecast length are the same.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_interval) :: fu_accumulation_field_id
    !
    ! Imported parameters with intent(in):
    TYPE(silja_field_id), INTENT(in) :: field_id

    fu_accumulation_field_id = field_id%length_of_accumulation

  END FUNCTION fu_accumulation_field_id


  ! ***************************************************************

  FUNCTION fu_validity_length_of_field_id(field_id)

    ! Description:
    ! Returns the of accumulation of a field-id. That is
    ! the interval from valid time backward that is took the
    ! quantity to valid time. If accumulation starts from analysis
    ! (as usually is) then accumulation and forecast length are the same.
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_interval) :: fu_validity_length_of_field_id
    !
    ! Imported parameters with intent(in):
    TYPE(silja_field_id), INTENT(in) :: field_id

    fu_validity_length_of_field_id = field_id%length_of_validity

  END FUNCTION fu_validity_length_of_field_id


  ! ***************************************************************

  INTEGER FUNCTION fu_quantity_of_field_id(field_id)

    ! Description:
    ! Returns the quantity of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent(in):
    TYPE(silja_field_id), INTENT(in) :: field_id

    fu_quantity_of_field_id = field_id%quantity

  END FUNCTION fu_quantity_of_field_id


  ! ***************************************************************
  
  FUNCTION fu_species_of_field_id (id) result(species)
    !
    ! Returns the substance name of a field-id.
    !
    IMPLICIT NONE
    
    ! Return value
    type(silam_species) :: species

    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_field_id), intent(in) :: id

    if(defined(id))then
      species = id%species
    else
      call set_error('undefined field ID','fu_species_of_field_id')
      species = species_missing
    end if
    
  END FUNCTION fu_species_of_field_id


  !*****************************************************************

  function fu_cocktail_name_of_field_id(id) result(chCocktail)
    !
    ! Returns the name of the cocktail
    !
    implicit none

    character(len=substNmLen) :: chCocktail

    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_field_id), intent(in) :: id

    if(defined(id))then
      chCocktail = id%chCocktail
    else
      call set_error('undefined field ID','fu_cocktail_name_of_field_id')
      chCocktail = ''
    end if

  end function fu_cocktail_name_of_field_id


  ! ***************************************************************

  FUNCTION fu_chSubstNm_of_field_id(field_id)

    ! Returns the substance name of a field-id.
    !
    IMPLICIT NONE

    ! Return value
    character(len=clen) :: fu_chSubstNm_of_field_id
    ! Imported parameters with intent(in):
    TYPE(silja_field_id), INTENT(in) :: field_id
    
    fu_chSubstNm_of_field_id = fu_substance_name(field_id%species)

  END FUNCTION fu_chSubstNm_of_field_id


  ! ***************************************************************

  FUNCTION fu_mode_of_field_id(field_id) result(mode)

    ! Returns the mode (section) number of a field-id.
    !
    IMPLICIT NONE

    ! Imported parameters with intent(in):
    TYPE(silja_field_id), INTENT(in) :: field_id
    type(Taerosol_mode) :: mode

    mode = field_id%species%mode

  END FUNCTION fu_mode_of_field_id


  ! ***************************************************************

  real FUNCTION fu_optical_wave_len_field_id(field_id)

    ! Returns the mode (section) number of a field-id.
    !
    IMPLICIT NONE

    ! Imported parameters with intent(in):
    TYPE(silja_field_id), INTENT(in) :: field_id

    fu_optical_wave_len_field_id = field_id%species%wavelength

  END FUNCTION fu_optical_wave_len_field_id


  ! ***************************************************************

  INTEGER FUNCTION fu_field_kind(field_id)

    ! Description:
    ! Returns the kind of the field-id.
    !
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent(in):
    TYPE(silja_field_id), INTENT(in) :: field_id

    fu_field_kind= field_id%field_kind

  END FUNCTION fu_field_kind


  ! *****************************************************************

  LOGICAL FUNCTION fu_field_id_defined(id)

    ! Description:
    ! Returns a true value, if the identificaion is defined by
    ! setting correct values to it using the set-function.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silja_field_id), INTENT(in) :: id

    fu_field_id_defined = fu_true(id%defined)

  END FUNCTION fu_field_id_defined


  ! *****************************************************************

  LOGICAL FUNCTION fu_accumulated_id(id)

    ! Description:
    ! Returns a true value, if field is time-accumulated type.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silja_field_id), INTENT(in) :: id

    IF (defined(id)) THEN
      fu_accumulated_id = id%field_kind == accumulated_flag
    ELSE
      fu_accumulated_id = .false.
    END IF   

  END FUNCTION fu_accumulated_id

  !*************************************************************************

  subroutine set_field_kind_of_field_id(id, field_kind)
    !
    ! Encapsulation of field_kind
    !
    implicit none

    TYPE(silja_field_id), INTENT(inout) :: id
    integer, intent(in) :: field_kind

    if (any(field_kind == legal_field_kinds))then
       id%field_kind = field_kind
    else
       call set_error("Starnge field_kind requested:"+fu_str(field_kind), "set_field_kind_of_field_id")
    endif
  end subroutine set_field_kind_of_field_id


  !**********************************************************************
  
  subroutine set_grid_of_field_id(id, grid)
    !
    ! Just sets the grid to the field id
    !
    implicit none

    ! Imported parameters
    TYPE(silja_field_id), INTENT(inout) :: id
    TYPE(silja_grid), INTENT(in) :: grid

    id%grids(1) = grid

  end subroutine set_grid_of_field_id


  !**********************************************************************
  
  subroutine set_level_of_field_id(id, level)
    !
    ! Just sets the level to the field id
    !
    implicit none

    ! Imported parameters
    TYPE(silja_field_id), INTENT(inout) :: id
    TYPE(silja_level), INTENT(in) :: level

    id%level = level

  end subroutine set_level_of_field_id


  !***********************************************************************

  subroutine set_acc_len_of_field_id(id, acc_len)
    !
    ! Sets the accumulation length for the accumulated field
    !
    implicit none

    ! Imported parameters
    TYPE(silja_field_id), INTENT(inout) :: id
    type(silja_interval), intent(in) :: acc_len

    if(id%field_kind == accumulated_flag .or. &
     & id%field_kind == averaged_flag) then ! .or. &
!     & id%field_kind == accumulated_period_valid_flag)then
      id%length_of_accumulation = acc_len
    else
      if(acc_len > zero_interval)then
        call set_error('Not accumulated/averaged field','set_acc_len_of_field_id')
        return
      endif
    endif

  end subroutine set_acc_len_of_field_id


  !*************************************************************************

  subroutine set_valid_len_of_field_id(id, valid_len)
    !
    ! Sets the accumulation length for the accumulated field
    !
    implicit none

    ! Imported parameters
    TYPE(silja_field_id), INTENT(inout) :: id
    type(silja_interval), intent(in) :: valid_len
    !
    ! Whatever the field kind is, it can be valid for more than 1 nanosecond
    !
    id%length_of_validity = valid_len

  end subroutine set_valid_len_of_field_id

  
  !*************************************************************************

  subroutine set_valid_time_of_field_id(id, valid_time)
    !
    ! Encapsulation of valid_time
    !
    implicit none

    TYPE(silja_field_id), INTENT(inout) :: id
    type(silja_time), intent(in) :: valid_time

    id%valid_time = valid_time
    if(defined(valid_time) .and. defined(id%analysis_time))then
      id%forecast_length = id%valid_time - id%analysis_time
    else
      id%forecast_length = interval_missing
    endif
  end subroutine set_valid_time_of_field_id

  !*************************************************************************

  subroutine set_analysis_time_of_field_id(id, analysis_time)
    !
    ! Encapsulation of valid_time
    !
    implicit none

    TYPE(silja_field_id), INTENT(inout) :: id
    type(silja_time), intent(in) :: analysis_time

    id%analysis_time = analysis_time
    id%forecast_length = id%valid_time - id%analysis_time
  end subroutine set_analysis_time_of_field_id


  !*************************************************************************

  subroutine set_met_src_of_field_id(id, met_src)
    !
    ! Encapsulation of valid_time
    !
    implicit none

    TYPE(silja_field_id), INTENT(inout) :: id
    type(meteo_data_source), intent(in) :: met_src

    id%met_src = met_src
  end subroutine set_met_src_of_field_id

  !*************************************************************************

  subroutine set_species_of_field_id(id, species)
    implicit none

    TYPE(silja_field_id), INTENT(inout) :: id
    type(silam_species), intent(in) :: species
    
    id%species = species

  end subroutine set_species_of_field_id

  !*************************************************************************


  subroutine set_cocktail_name_of_field_id(id, chCocktail)
    implicit none

    TYPE(silja_field_id), INTENT(inout) :: id
    character(len=*), intent(in) :: chCocktail
    
    id%chCocktail = chCocktail

  end subroutine set_cocktail_name_of_field_id

  !*************************************************************************
  
  logical function fu_if_internal_silam_field(id)
    implicit none

    TYPE(silja_field_id), INTENT(in) :: id

    fu_if_internal_silam_field = (fu_met_src(id) == silam_internal_src)
  end function fu_if_internal_silam_field


  ! ***************************************************************
  ! ***************************************************************
  ! 
  !
  !                 TESTS
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  SUBROUTINE print_field_id_report(id)

    ! Description:
    ! Prints a report of a field identification section to screen.
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    TYPE(silja_field_id), INTENT(in) :: id
    !
    ! Local declarations:
    INTEGER :: i = 0

    IF (.NOT.fu_true(id%defined)) THEN
      call msg(' field id not defined ')
      RETURN
    END IF

    call msg('Quantity: ' // fu_quantity_string(id%quantity))

    if(defined(id%species)) call report(id%species)

    if(len_trim(id%chCocktail) > 0) call msg('Species mixture:' + id%chCocktail)

    CALL report(id%level)

    call msg("Field kind", id%field_kind)

    IF (defined(id%analysis_time)) &
      & call msg(fu_connect_strings('Analysis time:',fu_str(id%analysis_time)))

    IF (defined(id%forecast_length)) &
              & call msg(fu_connect_strings('Forecast length:', fu_str(id%forecast_length)))

    IF (defined(id%length_of_accumulation)) THEN
!      if(id%length_of_accumulation > zero_interval)then
        call msg(fu_connect_strings('Start time of accumulation: ',&
                       & fu_str(id%valid_time-id%length_of_accumulation)))
        call msg(fu_connect_strings('From which accumulated for ', fu_str(id%length_of_accumulation)))
!      endif
    END IF

    call msg('Valid time:' + fu_str(id%valid_time))

    if(defined(id%length_of_validity)) then 
         call msg('Length of validity: ' + fu_str(id%length_of_validity))
    else
         call msg('Length of validity: UNDEFINED')
    endif

    SELECT CASE (id%quantity)

      CASE (wind_flag)

      call msg(' U-grid: ')
      CALL report(id%grids(1))

      call msg(' V-grid: ')
      CALL report(id%grids(2))

      if(defined(id%grids(3)))then
        call msg(' W-grid: ')
        CALL report(id%grids(3))
      else
        call msg('No W component')
      endif

    CASE default

      CALL report (id%grids(1))

    END SELECT

  END SUBROUTINE print_field_id_report

END MODULE field_identifications
