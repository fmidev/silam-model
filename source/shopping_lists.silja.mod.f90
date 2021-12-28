MODULE shopping_lists

  ! Description:
  ! This modulue contains the type silja_shopping_list, which defined
  ! in several ways a set of data needed from an external source. A
  ! shopping list should typically contain data that is put to one
  ! stack in one call to a data-access routine. 
  !
  ! A new addition: type silam_shopping_variable. If an input file
  ! contains largely excessive information, for example, several grids
  ! for one quantity - the trouble is inevitable for the old shopping_list.
  ! This new type becomes the main entity of the shopping list. It
  ! contains exact combination of the quantity, grid, vertical and data
  ! source indentifications, which enables the selection of data in a 
  ! very tiny way.Of course, any element of the combination may have
  ! missing value, which would mean that any value is OK.
  !
  ! Original code: Mika Salonoja
  ! Author: Mikhail Sofiev, FMI mikhail.sofiev@fmi.fi
  !
  ! All units: SI
  ! 
  ! Language: ANSI standard Fortran 90
  ! 
  ! Modules used:
!  USE globals
!  USE times
!  USE levels
  !use netcdf_io
!  use grib_code_table
  USE field_identifications

  IMPLICIT NONE

  ! The public functions and subroutines available in this module:
  PUBLIC fu_set_shopping_list
  PUBLIC set_missing
  PUBLIC fix_shopping_time_boundaries
  PUBLIC fix_shopping_quantities
  PUBLIC fix_shopping_levels
  public fix_met_src
  PUBLIC add_shopping_quantities
  PUBLIC defined
  PUBLIC report
  PUBLIC fu_field_id_in_list
  PUBLIC fu_field_id_in_list_of_vars
  PUBLIC fu_met_src     ! from the main shopping list
  PUBLIC fu_quantity_in_list
  PUBLIC fu_start_time
  PUBLIC fu_end_time
  PUBLIC fu_quantities
  PUBLIC fu_nbr_of_quantities
  PUBLIC fu_requests
  public fu_quantity
  public fu_request
  public fu_if2D
  public fu_if2Ds
  public fu_mlev_quantity_in_list
  public set_request ! Changes the request value
  public set_list_time_indicator
  public fu_list_time_indicator
  public replace_quantity

  PUBLIC add_shopping_variable
  public fu_shopping_var
!  public fu_field_id  ! from shopping variable
  public fu_nbr_of_vars
  public fu_fld_corresponds_to_shop_var
  public clean_shop_var
  public fu_nbr_of_MDS  ! in shopping vars
  public set_quantity
  public replace_vertical_in_list

  public get_var_map_targets_and_factors
  public get_var_mapping
  public add_maplink_to_shopVar
  public fu_nrTargets
  public fu_time_in_list
  public fu_if_var_mapped


  ! The private functions and subroutines not to be used elsewhere:
  private set_shopping_list_missing
  private set_shopping_variable_missing
  PRIVATE fu_met_src_of_shopping
  PRIVATE fu_shopping_list_defined
  PRIVATE fu_shopping_variable_defined
  PRIVATE print_shopping_report
  private report_shopping_variable
  private fu_nbr_of_quantities_in_list
  private fu_quantity_from_shopping_list
  private fu_all_quantities_from_list
  private fu_request_from_shopping_list
  private fu_all_requests_from_list
  private fu_if2D_of_shopping_quantity
  private fu_all_if2D_from_shopping_list
  private set_request_in_shopping_list

  private fu_compare_variables_eq ! Interface for ==
  private add_shopping_var_by_params
  private add_shopping_var_by_var
  private add_shopping_var_by_fieldID
!  private fu_field_id_from_shopping_var
  private replace_quantity_in_shop_list
  private set_quantity_of_shop_var
  private fu_quantity_of_shopping_var
  private fu_shop_var_map_defined
  private set_var_mapping_missing
  private fu_set_var_mapping
  private fu_shlist_start_time
  private fu_shlist_end_time

  ! Generic names and operator-interfaces of some functions:
  INTERFACE defined
    MODULE PROCEDURE fu_shopping_list_defined
    MODULE PROCEDURE fu_shopping_variable_defined
    module procedure fu_shop_var_map_defined
  END INTERFACE

  interface set_missing
    module procedure set_shopping_list_missing
    module procedure set_shopping_variable_missing
    module procedure set_var_mapping_missing
  end interface

  INTERFACE fu_met_src 
    MODULE PROCEDURE fu_met_src_of_shopping
  END INTERFACE

  INTERFACE report
    MODULE PROCEDURE print_shopping_report
    module procedure report_shopping_variable
  END INTERFACE

  INTERFACE fu_start_time
    MODULE PROCEDURE fu_shlist_start_time
  END INTERFACE
  
  INTERFACE fu_end_time
    MODULE PROCEDURE fu_shlist_end_time
  END INTERFACE

  interface fu_nbr_of_quantities
    module procedure fu_nbr_of_quantities_in_list
  end interface

  interface fu_quantity
    module procedure fu_quantity_from_shopping_list
    module procedure fu_quantity_of_shopping_var
  end interface

  interface fu_quantities
    module procedure fu_all_quantities_from_list
  end interface

  interface  fu_request
    module procedure fu_request_from_shopping_list
  end interface

  interface  fu_requests
    module procedure fu_all_requests_from_list
  end interface

  interface  fu_if2D
    module procedure fu_if2D_of_shopping_quantity
  end interface

  interface  fu_if2Ds
    module procedure fu_all_if2D_from_shopping_list
  end interface

  interface  set_request
    module procedure set_request_in_shopping_list
  end interface

  interface replace_quantity
    module procedure replace_quantity_in_shop_list
  end interface

  interface set_quantity
    module procedure set_quantity_of_shop_var
  end interface

  interface operator(==)
    module procedure fu_compare_variables_eq
  end interface

  interface add_shopping_variable
    module procedure add_shopping_var_by_params
    module procedure add_shopping_var_by_var
    module procedure add_shopping_var_by_fieldID
  end interface
  
!  interface fu_field_id
!    module procedure fu_field_id_from_shopping_var
!  end interface

  ! Public types with private components defined in this module:


  type inVar2modVarsMap
    private
    type(silja_field_id), dimension(max_quantities):: modFieldIds
    real, dimension(max_quantities) :: in2modFactor
    integer :: nTargets = 0
    type(silja_logical) :: defined
  endtype inVar2modVarsMap 

  type silam_shopping_variable
    private
    integer :: quantity, request, vertLevNbr
    type(silam_species) :: species
    character(len=substNmLen) :: chCocktailNm ! single subst or mixture
    type(silam_vertical) :: vertical
    type(silja_grid) :: grid
    type(meteo_data_source) :: mds
    type(inVar2modVarsMap), pointer :: mapping ! one input variable might contribute to several model fields
    type(silja_logical) :: defined
  end type silam_shopping_variable


  TYPE silja_shopping_list
    PRIVATE 
    !    LOGICAL :: take_everything ! if true then accepts all data
    type(meteo_data_source) :: met_src
    INTEGER :: time_indicator
    TYPE(silja_time) :: earliest_valid_time
    TYPE(silja_time) :: latest_valid_time
    INTEGER :: level_indicator
    TYPE(silja_level) :: floor_level !boundaries in vertical
    TYPE(silja_level) :: ceiling_level
    type(silja_grid) :: grid
    INTEGER, DIMENSION(max_quantities) :: quantities, request
    type(silja_logical), DIMENSION(max_quantities) :: if2D
    integer :: nVars
    type(silam_shopping_variable), dimension(max_quantities) :: vars
    TYPE(silja_logical) :: defined
  END TYPE silja_shopping_list
  

  ! Possible values for level_indicator:
  INTEGER, PARAMETER, PRIVATE :: one_level_only = 1 ! take data only from level 1
  INTEGER, PARAMETER, PRIVATE :: between_levels = 2 ! take data from 
                                   ! all heights between level1 and level2
  INTEGER, PARAMETER, PRIVATE :: accept_all_levels = 3 ! take data from all heights

  ! Possible values for time_indicator:
  INTEGER, PARAMETER, public :: one_time_only = 621 ! take data only from one time 
                                           ! (earliest_valid_time and earliest_valid_time are equal)
  INTEGER, PARAMETER, public :: between_times = 622 ! take data from all times between t1 and t2
  INTEGER, PARAMETER, public :: accept_all_times = 623 ! take data from all times
  INTEGER, PARAMETER, public :: accept_same_month = 624  ! take data with same month,
                                                         ! ignoring year, day and time 


  ! Flag used in data retrieval, when all quantities are accepted,
  ! or when data is treated as black box and source is irrelevant:
  !MAS. Substituted with int_missing - in case of absent source it does
  !     not play any role anymore
!  INTEGER, PARAMETER, PUBLIC :: accept_all_met_sources = 229999


  type(silam_shopping_variable), parameter, public :: variable_missing = &
      & silam_shopping_variable(int_missing,  0,  int_missing, &
                             & species_missing, '',&
                             &  vertical_missing, &
                             &  grid_missing, &
                             &  met_src_missing, &
                             & null(), &
                             &  silja_false)

  type(inVar2modVarsMap), parameter, public :: varMap_missing = &
             & inVar2modVarsMap(field_id_missing, &
             & real_missing, &
             & 0, &
             & silja_false)

  TYPE (silja_shopping_list),  parameter, public :: shopping_list_missing = &
               &  silja_shopping_list(met_src_missing, int_missing, time_missing,  time_missing, &
               & int_missing, level_missing,  level_missing, grid_missing, int_missing, int_missing, &
               & silja_undefined, int_missing, variable_missing,  silja_false)

               
    
CONTAINS


  !********************************************************************

  subroutine set_shopping_list_missing(shList)
    !
    ! Nullifies the given shopping list - a substitution to shopping_list_missing constant
    !
    implicit none

    TYPE(silja_shopping_list), intent(out) :: shList

    ! Local variables
    integer :: i

    shList%met_src = met_src_missing
    shList%time_indicator = int_missing
    shList%earliest_valid_time = time_missing
    shList%latest_valid_time = time_missing
    shList%level_indicator = int_missing
    shList%floor_level = level_missing
    shList%ceiling_level = level_missing
    shList%grid = grid_missing
    shList%quantities = int_missing
    shList%request = int_missing
    shList%if2D = silja_undefined
    shList%nVars = 0
    do i=1, size(shlist%vars)
      call set_missing(shList%vars(i))
    end do
    shList%defined = silja_false

  end subroutine set_shopping_list_missing


  !***************************************************************************

  subroutine set_shopping_variable_missing(var)
    !
    ! Nullifies the variable
    !
    implicit none

    ! Imported parameter
    type(silam_shopping_variable), intent(out) :: var

    var%quantity=int_missing
    var%request=int_missing
    var%vertLevNbr=int_missing
    call set_missing(var%vertical, .true.)
    var%grid = grid_missing
    var%mds = met_src_missing
    var%species = species_missing
    var%chCocktailNm = ''
    var%defined = silja_false
    nullify(var%mapping)
  end subroutine set_shopping_variable_missing


  ! ***************************************************************
  
  FUNCTION fu_set_shopping_list (met_src, &
                               & quantities, &
                               & first_time_boundary, &
                               & second_time_boundary, &
                               & floor_level,&
                               & ceiling_level, &
                               & grid, &
                               & requests, &
                               & if2D) result(list)
    ! Sets a shopping list. 
    ! 
    ! If data retrieval is not determined by source (producer) but
    ! data from any source that matches other requirements is wanted,
    ! then set source = int_missing <=> accept_all_sources
    !
    ! If all available quantities from a given source (or file) are wanted
    ! then set quantities(1)= accept_all_quantities. In this case data may be
    ! searched by other qualities: levels, times and source.
    !
    ! Data is taken from all times between given boundaries.
    ! The order of boundaries doesn't matter.
    ! If data for only one time is wanted, then both boundaries
    ! should have this one value.
    ! If data from all possible times is
    ! wanted, then time-boundaries should be undefined (time_missing).
    ! 
    ! Data is taken from all levels between given boundaries.
    ! If data for only one height is wanted, then both upper
    ! and lower level should have this one value.
    ! If data from all possible levels is
    ! wanted, then level-boundaries should be undefined (level_missing).
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_shopping_list) :: list

    ! Imported parameters with intent(in):
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), INTENT(in) :: first_time_boundary,  second_time_boundary
    TYPE(silja_level), INTENT(in) :: floor_level, ceiling_level
    INTEGER, DIMENSION(:), INTENT(in) :: quantities
    INTEGER, DIMENSION(:), INTENT(in), optional :: requests
    type(silja_grid), intent(in), optional :: grid
    type(silja_logical), dimension(:), intent(in), optional :: if2D


    ! Local variables
    integer :: i, j
    !----------------------------------------
    !
    ! 1. Set empty first.
    !    ---------------

    call set_missing(list)
    list%nVars = 0

    !----------------------------------------
    !
    ! 2. Set time-boundaries
    !    -------------------

    IF ((.NOT.defined(first_time_boundary)).or.&
      & (.NOT.defined(second_time_boundary))) THEN
        list%time_indicator = accept_all_times
    ELSE
      IF (first_time_boundary < second_time_boundary) THEN
        list%earliest_valid_time = first_time_boundary
        list%latest_valid_time = second_time_boundary
      ELSE
        list%earliest_valid_time = second_time_boundary
        list%latest_valid_time = first_time_boundary
      END IF

      IF (first_time_boundary == second_time_boundary) THEN
        list%time_indicator = one_time_only
      ELSE
        list%time_indicator = between_times
      END IF

    END IF


    !----------------------------------------
    !
    ! 3. Set level-boundaries.
    !    --------------------

    IF ((.NOT.defined(floor_level)).or.(.NOT.defined(ceiling_level))) THEN
      list%level_indicator = accept_all_levels
    ELSE
      list%floor_level = floor_level
      list%ceiling_level = ceiling_level

      IF (fu_cmp_levs_eq(floor_level, ceiling_level)) THEN
        list%level_indicator = one_level_only
      ELSE
        list%level_indicator = between_levels
      END IF

    END IF


    !----------------------------------------
    !
    ! 4. Set quantites and their request types.
    !    If there is no quantities - issue warning but do not set error
    !

    IF (.NOT.(fu_known_quantity(quantities(1)).or.&
        & (quantities(1) == accept_all_quantities))) THEN
      CALL msg_warning('no quantities given','fu_set_shopping_list')
      list%quantities = int_missing
    else

      do i=1, size(quantities)
        if(quantities(i) == int_missing)exit
          list%quantities(i) = quantities(i)
          if(present(requests))then
            list%request(i) = requests(i)
          else
            if(fu_known_quantity(list%quantities(i)))then
              list%request(i) = 2
            else
              list%request(i) = 0
            endif
          endif
          if(present(if2D))list%if2D(i) = if2D(i)
      end do

    endif

    !----------------------------------------
    !
    ! 5. Set meteo source and grid
    !    -----------

    list%met_src = met_src

    if(present(grid))then
      list%grid = grid
    else
      list%grid = grid_missing
    end if

    list%defined = fu_set_true()

!    call report(list)

  END FUNCTION fu_set_shopping_list


! ***************************************************************

  SUBROUTINE fix_shopping_time_boundaries(list, new_earliest_time, new_latest_time)
    !
    ! Sets new values for time-boundaries of shopping list.
    !
    IMPLICIT NONE

    ! Imported parameters:
    TYPE(silja_time), INTENT(in) :: new_earliest_time, new_latest_time
    TYPE(silja_shopping_list), intent(inout) :: list

    IF ((.NOT.defined(new_earliest_time)) .or. (.NOT.defined(new_latest_time))) THEN
      CALL set_error('undefined times given','fix_time_boundaries')
      RETURN
    END IF

    IF (new_earliest_time < new_latest_time) THEN
      list%earliest_valid_time = new_earliest_time
      list%latest_valid_time = new_latest_time
    ELSE
      list%earliest_valid_time = new_latest_time
      list%latest_valid_time = new_earliest_time
    END IF

    IF (list%earliest_valid_time == list%latest_valid_time) THEN
      list%time_indicator = one_time_only
    ELSE
      list%time_indicator = between_times
    END IF
    
  END SUBROUTINE fix_shopping_time_boundaries


  ! ***************************************************************

  SUBROUTINE fix_shopping_quantities(shList, new_quantities, requests, if2D)

    ! Description:
    ! Sets new values for quantities of shopping shList.
    ! Does not check for duplicates
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    INTEGER, DIMENSION(:), INTENT(in) :: new_quantities
    INTEGER, DIMENSION(:), INTENT(in), optional :: requests
    type(silja_logical), DIMENSION(:), INTENT(in), optional :: if2D

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_shopping_list), INTENT(inout) :: shList

    INTEGER :: i, qc

    IF (.NOT.defined(shList)) THEN
      CALL set_error('undefined shList given','fix_shopping_quantities')
      RETURN
    END IF

    qc = 1
    shList%quantities = int_missing

    DO i = 1, SIZE(new_quantities)
      IF (fu_known_quantity(new_quantities(i))) THEN
        if(present(requests))then
          if(requests(i) > 0)then
            shList%quantities(qc) = new_quantities(i)
            shList%request(qc) = requests(i)
            if(present(if2D))then
              if(i > size(if2D))then
                call set_error('Too few if2D flags','fix_shopping_quantities')
                return
              endif
              shList%if2D(qc) = if2D(i)
            endif
            qc = qc + 1
          endif
        else
          shList%quantities(qc) = new_quantities(i)
          shList%request(qc) = 2 ! mandatory
          if(present(if2D))then
            if(i > size(if2D))then
              call set_error('Too few if2D flags','fix_shopping_quantities')
              return
            endif
            shList%if2D(qc) = if2D(i)
          endif
          qc = qc + 1
        endif
      END IF
    END DO

  END SUBROUTINE fix_shopping_quantities


  ! ***************************************************************

  SUBROUTINE fix_shopping_levels(shList, level_bottom, level_top)

    ! Description:
    ! Sets new values for levels of shopping shList.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silja_level), INTENT(in) :: level_bottom, level_top

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_shopping_list), INTENT(inout) :: shList

    IF (.NOT.defined(shList)) THEN
      CALL set_error('undefined shList given','fix_shopping_levels')
      RETURN
    END IF

    IF ((.NOT.defined(level_bottom)).or.(.NOT.defined(level_top))) THEN
      shList%level_indicator = accept_all_levels
    ELSE
      shList%floor_level = level_bottom
      shList%ceiling_level = level_top

      IF (fu_cmp_levs_eq(shList%floor_level, shList%ceiling_level)) THEN
        shList%level_indicator = one_level_only
      ELSE
        shList%level_indicator = between_levels
      END IF

    END IF
  END SUBROUTINE fix_shopping_levels



  ! ***************************************************************

  SUBROUTINE add_shopping_quantities(shList, new_quantities, requests, if2D)
    ! 
    ! Adds new values for quantities of shopping shList, including, if present,
    ! their request types
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    INTEGER, DIMENSION(:), INTENT(in) :: new_quantities
    INTEGER, DIMENSION(:), INTENT(in), optional :: requests
    type(silja_logical), DIMENSION(:), INTENT(in), optional :: if2D

    ! Imported parameters with intent INOUT
    TYPE(silja_shopping_list), INTENT(inout) :: shList

    INTEGER :: i, qc, iListSize

    IF (.NOT.defined(shList)) THEN
      CALL set_error('undefined shList given','add_shopping_quantities')
      RETURN
    END IF

    iListSize =SIZE(shList%quantities)

    !--------- Find the first empty quantity in original shList

    FindEndOfList: DO qc = 1, iListSize
      IF(shList%quantities(qc) == int_missing) exit FindEndOfList
    END DO FindEndOfList

    !--------- Now - go through the new shList and add what is new and known

    DO i = 1, SIZE(new_quantities)
      if (new_quantities(i) == int_missing) exit
      IF (fu_known_quantity(new_quantities(i))) THEN
        IF(ALL(shList%quantities(1:qc) /= new_quantities(i))) THEN
!          print *, 'added ', new_quantities(i)
          if(present(requests))then
            if(requests(i) > 0)then
              shList%quantities(qc) = new_quantities(i)
              shList%request(qc) = requests(i)
              if(present(if2D))then
                if(i > size(if2D))then
                  call set_error('Too few if2D flags','add_shopping_quantities')
                  return
                endif
                shList%if2D(qc) = if2D(i)
              endif
              qc = qc + 1
            endif
          else
            shList%quantities(qc) = new_quantities(i)
            shList%request(qc)=2 ! Mandatory one
            if(present(if2D))then
              if(i > size(if2D))then
                call set_error('Too few if2D flags','add_shopping_quantities')
                return
              endif
              shList%if2D(qc) = if2D(i)
            endif
            qc = qc + 1
          endif
          IF(qc > iListSize)THEN
            CALL set_error('Too many quantities','add_shopping_quantities')
            RETURN
          END IF
        END IF
      else
         call msg_warning("Unknown quantity","add_shopping_quantities")
         call msg("Quantities so far:",new_quantities(1:i))
         call msg("Requests   so far:",requests(1:i))
         call set_error("","add_shopping_quantities")
      END IF
    END DO

  END SUBROUTINE add_shopping_quantities


  !*****************************************************************

  subroutine add_shopping_var_by_params(shList,quantity,species,grid,vert,vertLevNbr,&
                                      & mds,request,mapping, chCocktail)
    !
    ! Adds one more shopping variable to the shopping shList.
    !
    implicit none

    ! Imported variables with intent IN
    integer, intent(in) :: quantity, vertLevNbr
    type(silam_species), intent(in) :: species
    integer, intent(in), optional :: request
    type(silja_grid), intent (in) :: grid
    type(silam_vertical), intent(in) :: vert
    type(meteo_data_source), intent(in) :: mds
    type(inVar2modVarsMap), pointer, optional :: mapping
    character(len=*), intent(in), optional :: chCocktail

    ! Imported variable with intent INOUT
    type(silja_shopping_list), intent(inout) :: shList

    ! Local valriables
    integer :: i

!    print *, 'adding variable'
    
    if(.not.defined(shList))then
      call msg_warning('Undefined shList given','add_shopping_var_by_params')
      shList%defined = silja_true
      shList%time_indicator = accept_all_times
      shList%nVars = 0
    end if
    !
    ! Find the first non-empty variable, also checking if the variable to add
    ! is already there.
    !
    do i = 1, shList%nVars
      if(.not.defined(shList%vars(i)))then
        call set_error('Empty variable in the shList','add_shopping_var_by_params')
        return
      end if
      if ( shList%vars(i)%quantity == quantity .and. &
        &  shList%vars(shList%nVars)%species == species .and.  &
        &  fu_cmp_verts_eq(shList%vars(i)%vertical, vert) .and. &
        &  shList%vars(i)%vertLevNbr == vertLevNbr .and. &
        &  shList%vars(i)%grid == grid .and. &
        &  shList%vars(i)%mds == mds)then
        if(present(chCocktail))then
          if(shList%vars(shList%nVars)%chCocktailNm == chCocktail) return
        else
          return
        endif
      endif
    end do

    if(shList%nVars == size(shList%vars))then
      call set_error('Too many variables','add_shopping_var_by_params')
      return
    else
      shList%nVars = shList%nVars + 1
    end if

    !
    ! Now - add variable
    !
    shList%vars(shList%nVars)%quantity = quantity
    shList%vars(shList%nVars)%species = species
    shList%vars(shList%nVars)%vertical = vert
    shList%vars(shList%nVars)%vertLevNbr = vertLevNbr
    shList%vars(shList%nVars)%grid = grid
    shList%vars(shList%nVars)%mds = mds
    if(present(mapping))then
      shList%vars(shList%nVars)%mapping => mapping 
    else
      nullify(shList%vars(shList%nVars)%mapping) 
    endif
    shList%vars(shList%nVars)%defined = silja_true
    if(present(request))then
      shList%vars(shList%nVars)%request = request
    else
      shList%vars(shList%nVars)%request = 2  ! Mandatory
    endif
    if(present(chCocktail))then
      shList%vars(shList%nVars)%chCocktailNm = chCocktail
    else
      shList%vars(shList%nVars)%chCocktailNm = ''
    endif
    
    call msg_test('Added to shopping shList (by params):' + &
                & fu_quantity_string(shList%vars(i)%quantity) + ',' + fu_str(shList%vars(i)%species))

  end subroutine add_shopping_var_by_params



  !*****************************************************************

  subroutine add_shopping_var_by_var(shList,varNew)
    !
    ! Adds one more shopping variable to the shopping shList.
    !
    implicit none

    ! Imported variables with intent IN
    type(silam_shopping_variable), intent(in) :: varNew

    ! Imported variable with intent INOUT
    type(silja_shopping_list), intent(inout) :: shList

    ! Local valriables
    integer :: i

!    print *, 'adding variable'
    
    if(.not.defined(shList))then
      call msg_warning('Undefined shList given','add_shopping_var_by_var')
      shList%defined = silja_true
      shList%time_indicator = accept_all_times
      shList%nVars = 0
    end if
    !
    ! Find the first non-empty variable, also checking if the variable to add
    ! is already there.
    !
    do i = 1, shList%nVars
      if(.not.defined(shList%vars(i)))then
        call set_error('Empty variable in the shList','add_shopping_var_by_var')
        return
      end if
      if ( shList%vars(i)%quantity == varNew%quantity .and. &
        &  shList%vars(shList%nVars)%species == varNew%species .and.  &
        &  shList%vars(shList%nVars)%chCocktailNm == varNew%chCocktailNm .and.  &
        &  fu_cmp_verts_eq(shList%vars(i)%vertical, varNew%vertical) .and. &
        &  shList%vars(i)%vertLevNbr == varNew%vertLevNbr .and. &
        &  shList%vars(i)%grid == varNew%grid .and. &
        &  shList%vars(i)%mds == varNew%mds) return
    end do

    if(shList%nVars == size(shList%vars))then
      call set_error('Too many variables','add_shopping_var_by_var')
      return
    else
      shList%nVars = shList%nVars + 1
    end if

    !
    ! Now - add variable
    !
    shList%vars(shList%nVars)%quantity = varNew%quantity
    shList%vars(shList%nVars)%species = varNew%species
    shList%vars(shList%nVars)%vertical = varNew%vertical
    shList%vars(shList%nVars)%vertLevNbr = varNew%vertLevNbr
    shList%vars(shList%nVars)%grid = varNew%grid
    shList%vars(shList%nVars)%mds = varNew%mds
    shList%vars(shList%nVars)%mapping = varNew%mapping 
    shList%vars(shList%nVars)%chCocktailNm = varNew%chCocktailNm
    shList%vars(shList%nVars)%defined = silja_true
    shList%vars(shList%nVars)%request = varNew%request

    call msg_test('Added to shopping shList (by var):' + fu_quantity_string(shList%vars(i)%quantity))

  end subroutine add_shopping_var_by_var


  !***********************************************************************

  subroutine add_shopping_var_by_fieldID(shList, fieldID, request, mapping)
    !
    ! Adds a new shopping variable using given field ID as a template
    !
    implicit none

    ! Imported parameters
    type(silja_shopping_list), intent(inout), target :: shList
    type(silja_field_id), intent(in) :: fieldID
    integer, intent(in) :: request
    type(inVar2modVarsMap), pointer, optional :: mapping

    ! Local variables
    type(silam_vertical) :: vertTmp
    !
    ! Just encapsulate adding via a set of parameters
    !
    if(fu_cmp_levs_eq(fu_level(fieldID),level_missing))then
      call set_missing(vertTmp, .true.)
    else
      call set_vertical(fu_level(fieldID), vertTmp)
    endif
    if(error)return
   
    if (present(mapping))then
      call add_shopping_var_by_params(shList, &
                                  & fu_quantity(fieldID), &
                                  & fu_species(fieldID), &
                                  & fu_grid(fieldID), &
                                  & vertTmp, &
                                  & 1, &  ! level number in the above vertical
                                  & fu_met_src(fieldID), &  
                                  & request, &
                                  & mapping, &
                                  & fu_cocktail_name(fieldID))
    else
      call add_shopping_var_by_params(shList, &
                                  & fu_quantity(fieldID), &
                                  & fu_species(fieldID), &
                                  & fu_grid(fieldID), &
                                  & vertTmp, &
                                  & 1, &  ! level number in the above vertical
                                  & fu_met_src(fieldID), &  
                                  & request, &
                                  & chCocktail = fu_cocktail_name(fieldID)) 
    endif

  end subroutine add_shopping_var_by_fieldID


  !******************************************************************

  function fu_shopping_var(shList, indexVar) result(varPtr)
    !
    ! Returns a pointer to a shopping variable, if idex is not outside the nVars
    !
    implicit none

    ! Return value
    type(silam_shopping_variable), pointer :: varPtr

    ! Imported parameters
    type(silja_shopping_list), pointer :: shList
    integer, intent(in) :: indexVar

    if(indexVar <= shList%nVars)then
      varPtr => shList%vars(indexVar)
    else
      call msg("shList%nVars, indexVar",shList%nVars, indexVar)
      call set_error('Index is bigger than the number of vars','fu_shopping_var')
      call msg("shList:")
      call report (shList)
      nullify(varPtr)
    endif

  end function fu_shopping_var


!  !***********************************************************************
!  
!  function fu_field_id_from_shopping_var(shlist, indVar) result(id)
!    !
!    ! Returns a field ID made of shopping variable, if idex is not outside the nVars
!    !
!    implicit none
!
!    ! Return value
!    type(silja_field_id) :: id
!
!    ! Imported parameters
!    type(silja_shopping_list), pointer :: shlist
!    integer, intent(in) :: indVar
!
!    if(indVar <= shlist%nVars)then
!      id = fu_set_field_id(met_src_missing, &                 ! met_src
!                         & shlist%vars(indexVar)%quantity, &  ! quantity
!                         & time_missing, &                    ! analysis_time
!                         & interval_missing, &                ! forecast_length
!                         &                                    ! grid,&
!                                     ! level,&
!                                     ! length_of_accumulation, & ! optional
!                                     ! length_of_validity, &     !optional
!                                     ! field_kind, &             ! optional
!                                     ! chSubstNm, &              ! optional
!                                     ! aerMode, &               ! optional
!                                     ! fWaveLen, &              ! optional
!                                     ! species, &               ! optional, alternative to the previous
!                                     ! chCocktail) &            ! optional, if a mixture of species
!!
!
!
!    type(silja_field_id), dimension(max_quantities):: modFieldIds
!    real, dimension(max_quantities) :: in2modFactor
!    integer :: nTargets = 0
!    type(silja_logical) :: defined
!  endtype inVar2modVarsMap 
!
!  type silam_shopping_variable
!    private
!    integer :: quantity, request, vertLevNbr
!    character(len=substNmLen) :: chSubstNm, chCocktailNm ! single subst or mixture
!    real :: fModeVal
!    type(silam_vertical) :: vertical
!    type(silja_grid) :: grid
!    type(meteo_data_source) :: mds
!    type(inVar2modVarsMap), pointer :: mapping ! one input variable might contribute to several model fields
!    type(silja_logical) :: defined
!  end type silam_shopping_variable
!
!
!
!    else
!      call set_error('Index is bigger than the number of vars','fu_shopping_var')
!      nullify(varPtr)
!    endif
! 
!  end function fu_field_id_from_shopping_var


  !*******************************************************************

  subroutine set_quantity_of_shop_var(var, quantity, species)

    implicit none

    type(silam_shopping_variable), intent(inout) :: var
    integer :: quantity
    type(silam_species), optional :: species

    var%quantity = quantity
    
    if(present(species))then
      var%species = species
    else
      var%species = species_missing
    endif

  end subroutine set_quantity_of_shop_var

  
  !*******************************************************************
  
  subroutine replace_vertical_in_list(list, vertical_from, vertical_to)
    !
    ! Replaces the vertical_from in all shopping variables, which have it, with vertical_to
    !
    implicit none
    
    ! Imported parameters
    type(silja_shopping_list), intent(inout) :: list
    type(silam_vertical), intent(in) :: vertical_from, vertical_to
    
    ! Local variables
    integer :: iVar
    
    do iVar = 1, list%nVars
      if(fu_cmp_verts_eq(list%vars(iVar)%vertical, vertical_from)) &
                                                & list%vars(iVar)%vertical = vertical_to
    end do
  end subroutine replace_vertical_in_list

  
  !*******************************************************************

  integer function fu_nbr_of_vars(list)
    !
    ! Returns the number of shopping variables
    !
    implicit none

    ! Imported parameters
    type(silja_shopping_list), intent(in), target :: list

    fu_nbr_of_vars = list%nVars

  end function fu_nbr_of_vars


  !*******************************************************************

  integer function fu_nbr_of_MDS(list)
    !
    ! Returns the number of shopping variables
    !
    implicit none

    ! Imported parameters
    type(silja_shopping_list), intent(in), target :: list

    ! Local variable
    integer :: i

    fu_nbr_of_MDS = 1
    
    do i=1,list%nVars
      if(.not. list%vars(i)%mds == list%vars(1)%mds) fu_nbr_of_MDS = fu_nbr_of_MDS + 1
    enddo

  end function fu_nbr_of_MDS


  !******************************************************************

  logical function fu_fld_corresponds_to_shop_var(fieldId, shopVar)
    !
    ! Checks of the field meets the limitations set by the shopping variable
    ! Important: should some parameter of the variable is missing, it is 
    ! considered as non-limiting
    !
    implicit none

    ! Imported parameters
    type(silja_field_id), intent(in) :: fieldId
    type(silam_shopping_variable), intent(in) :: shopVar

    fu_fld_corresponds_to_shop_var = .false.

    if(shopVar%quantity /= int_missing) then
      if(shopVar%quantity /= fu_quantity(fieldId))return
    endif

    if(defined(shopVar%species))then
      if(.not. (shopVar%species == fu_species(fieldID)))return
    else
      if(defined(fu_species(fieldID)))return
    endif

    if(len_trim(shopVar%chCocktailNm) > 0)then
      if(trim(shopVar%chCocktailNm) /= fu_cocktail_name(fieldID))return
    endif

    if(defined(shopVar%vertical))then
      if(fu_leveltype(shopVar%vertical) /= fu_leveltype(fu_level(fieldId))) return
      if(shopVar%vertLevNbr > 0)then
        if(.not. fu_cmp_levs_eq(fu_level(shopVar%vertical,shopVar%vertLevNbr), fu_level(fieldId)))return
      endif
    endif

    if(defined(shopVar%grid))then
      if(.not. shopVar%grid == fu_grid(fieldId))return
    endif

    if(.not. shopVar%mds == met_src_missing)then
      if(.not. shopVar%mds == fu_met_src(fieldId))return
    endif

    fu_fld_corresponds_to_shop_var = .true.

  end function fu_fld_corresponds_to_shop_var


  !***************************************************************************

  subroutine clean_shop_var(list, indexVar)
    !
    ! Cleans a variable pointed by index, if it is reasonable. The last
    ! variable is moved on-top of the current one
    !
    implicit none

    ! Imported parameters
    type(silja_shopping_list), intent(inout) :: list
    integer, intent(in) :: indexVar

    if(indexVar <= list%nVars)then
      if(indexVar <= list%nVars)then
        list%vars(indexVar)%quantity = list%vars(list%nVars)%quantity
        list%vars(indexVar)%vertical = list%vars(list%nVars)%vertical
        list%vars(indexVar)%vertLevNbr = list%vars(list%nVars)%vertLevNbr
        list%vars(indexVar)%grid = list%vars(list%nVars)%grid
        list%vars(indexVar)%mds  = list%vars(list%nVars)%mds
        list%vars(indexVar)%defined = list%vars(list%nVars)%defined
        list%vars(indexVar)%request = list%vars(list%nVars)%request
      endif

      list%vars(list%nVars)%quantity = int_missing
      list%vars(list%nVars)%species = species_missing
      list%vars(list%nVars)%chCocktailNm = ''
      call set_missing(list%vars(list%nVars)%vertical, .false.)
      list%vars(list%nVars)%vertLevNbr = int_missing
      list%vars(list%nVars)%grid = grid_missing
      list%vars(list%nVars)%mds = met_src_missing
      call set_missing(list%vars(list%nVars)%mapping)
      list%vars(list%nVars)%defined = silja_false
      list%vars(list%nVars)%request = int_missing

    else
      call set_error('Index is bigger than the number of vars','clean_shop_var')
    endif

  end subroutine clean_shop_var


  ! ****************************************************************

  subroutine fix_met_src(shopping_list, wdr, ifForceVars)
    !
    ! Checks the number of the sources in the wdr. So far the rule
    ! is straightforward: if there is only one source in wdr - set it
    ! in the shopping list. If there are several data sources - use
    ! met_src_missing to indicate that all should be accepted.
    !
    implicit none

    ! Imported parameters
    type(silja_wdr), intent(in) :: wdr
    type(silja_shopping_list), intent(inout) :: shopping_list
    logical, intent(in) :: ifForceVars

    ! Local variables
    integer :: iVar

    if(.not.defined(wdr))then
      call set_error('Undefined wdr given','fix_met_src')
      return
    end if

    if(.not.defined(shopping_list))then
      call set_error('Undefined shopping list given','fix_met_src')
      return
    end if

    if(fu_ifDisregardMDS(wdr) .or. fu_NbrOfMetSrcs(wdr) > 1)then
      shopping_list%met_src = met_src_missing
    else
      shopping_list%met_src = fu_met_src(wdr,1)
    end if

    !
    ! Reset own MDS of the variables only if all MDSs must be disregarded.
    !
    if(ifForceVars .and. fu_ifDisregardMDS(wdr))then
      do iVar = 1, shopping_list%nVars
        shopping_list%vars(iVar)%mds = met_src_missing
      enddo
    endif

  end subroutine fix_met_src

    !********************************************************

  logical function fu_shop_var_map_defined(mapping)
      
    implicit none

    type(inVar2modVarsMap), pointer :: mapping
      
    if(associated(mapping))then
      fu_shop_var_map_defined = fu_true(mapping%defined)
    else
      fu_shop_var_map_defined = .false.
    endif
 
  end function fu_shop_var_map_defined

  !********************************************************

  logical function fu_if_var_mapped(shopping_list, iVar)

    implicit none
    type(silja_shopping_list), intent(in) :: shopping_list
    integer, intent(in) :: iVar

    character (len=*), parameter :: sub_name="fu_if_var_mapped"

    fu_if_var_mapped = .false.
    if (.not. defined(shopping_list)) then
      call set_error("shopping_list undefined", sub_name)
    elseif (iVar < 1 .or. iVar > shopping_list%nvars ) then
      call msg("iVar ,shopping_list%nvars", iVar, shopping_list%nvars )
      call set_error("wrong ivar ", sub_name)
    elseif (.not. defined(shopping_list%vars(iVar))) then
      call set_error("undefined ivar ", sub_name)
    elseif ( associated(shopping_list%vars(iVar)%mapping)) then
      fu_if_var_mapped = defined(shopping_list%vars(iVar)%mapping)
    endif

  end function fu_if_var_mapped


  !********************************************************
  subroutine get_var_mapping(var, mapping)

    implicit none

    type(silam_shopping_variable), pointer :: var
    type(inVar2modVarsMap), pointer :: mapping

    if(associated(var%mapping))then
      mapping => var%mapping
      return
    else
!      mapping => varMap_missing      
      nullify(mapping)
    endif

  end subroutine get_var_mapping


  !********************************************************
 
  subroutine get_var_map_targets_and_factors(mapping, targets, factors)

    implicit none

    type(inVar2modVarsMap), pointer :: mapping
    type(silja_field_id), dimension(:), pointer :: targets
    real, dimension(:), pointer :: factors

    targets => mapping%modFieldIds
    factors => mapping%in2modFactor

   end subroutine get_var_map_targets_and_factors


  !********************************************************

  integer function fu_nrTargets(mapping)

    implicit none

    type(inVar2modVarsMap), intent(in) :: mapping
    
    fu_nrTargets = mapping%nTargets

  end function fu_nrTargets


  !********************************************************

  subroutine set_var_mapping_missing(mapping)

    implicit none

    type(inVar2modVarsMap) :: mapping

    mapping%nTargets = 0
    mapping%defined = silja_false

  end subroutine set_var_mapping_missing

  !********************************************************

  function fu_set_var_mapping(lstIds, factors, nTargets)

    implicit none

    type(inVar2modVarsMap) :: fu_set_var_mapping
    type(silja_field_id), dimension(:), intent(in) :: lstIds
    real, dimension(:), intent(in) :: factors
    integer, intent(in) :: nTargets
    integer :: iTmp

    fu_set_var_mapping%nTargets = nTargets

    do iTmp = 1, nTargets
      fu_set_var_mapping%modFieldIds(iTmp) = lstIds(iTmp)
      fu_set_var_mapping%in2modFactor(iTmp) = factors(iTmp)
    enddo

    fu_set_var_mapping%defined = silja_true


  end function fu_set_var_mapping
  

  !********************************************************

  subroutine add_maplink_to_shopVar(shpVar, idTarget, factor)

    implicit none

    type(silam_shopping_variable), pointer :: shpVar
    type(silja_field_id), intent(in) :: idTarget
    real, intent(in) :: factor


    if(.not. associated(shpVar%mapping))then
      allocate(shpVar%mapping)
    endif
    shpVar%mapping%nTargets = shpVar%mapping%nTargets + 1
    shpVar%mapping%modFieldIds(shpVar%mapping%nTargets) = idTarget
    shpVar%mapping%in2modFactor(shpVar%mapping%nTargets) = factor

    shpVar%mapping%defined = silja_true

  end subroutine add_maplink_to_shopVar

  ! ****************************************************************
  ! ****************************************************************
  ! 
  !
  !       HERE BEGINS USAGE OF SHOPPING LIST TO GET FIELDS
  !
  !
  ! ***************************************************************
  ! ***************************************************************
  
  LOGICAL FUNCTION fu_field_id_in_list(id, list, indexVar, data_time_features, time_to_force) &
                 & result(accept_field)

    ! Returns a true value, if the given field-id indicates a
    ! necessary field according to the shopping-list.
    !
    ! There are two ways of belonging to the list - belonging to the list of
    ! variables or satisfaction of the criteria like quantity, height range, 
    ! etc.
    ! The algorithm is the following. If the id belongs to the list of
    ! variables - only time range is checked. If the list of variables 
    ! does not accept the id - it is still checked for other general
    ! criteria. So, EITHER list of variables + time, OR general criteria
    ! + time may allow the field to come.
    !
    ! The reason for such tricky checking is purely historical. In principle,
    ! grib_analysis routines are deemed to describe complete set of 
    ! variables to shop. However, originally this technique was not used, so
    ! it will take time before it is eliminated from everywhere.
    !
    ! Post-processed/model level is not checked here, since it is not a
    ! quality of a single field.
    !
    ! The horizontal grid is not cheked here.
    !
    IMPLICIT NONE

    ! Imported parameters with intent(in):
    TYPE(silja_field_id), INTENT(in) :: id
    TYPE(silja_shopping_list), INTENT(in) :: list
    integer, intent(out) :: indexVar
    integer, optional, intent(in) :: data_time_features
    type(silja_time), optional, intent(out) :: time_to_force

    ! Local declarations:
    INTEGER :: i, quantity_in_field
    INTEGER :: lower_type, upper_type, level_type
    LOGICAL :: time_accept, quantity_accept
    LOGICAL :: lower_accept, upper_accept
    TYPE(silja_level) :: level
    type(silja_grid) :: gridTmp

    !----------------------------------------
    !
    ! 1. Some checking first.
    !
    IF (.NOT.defined(id)) THEN
      CALL set_error('field-id not defined','fu_field_id_in_list')
      accept_field = .false.
      RETURN
    END IF

    IF (.NOT.defined(list)) THEN
      CALL set_error('shopping-list not defined','fu_field_id_in_list')
      accept_field = .false.
      RETURN
    END IF

    if(present(time_to_force))time_to_force = time_missing

!call msg('Accept 1')
!call report(id)

    !---------------------------------------------
    !
    ! 2. Check the list of variables. If found there - check time and take it !
    ! If not - continue checking of more general rules
    !
!    if(present(indexVar))then
      accept_field = fu_field_id_in_list_of_vars(id, list, indexVar)
!    else
!      accept_field = fu_field_id_in_list_of_vars(id, list)
!    endif

!call msg('Accept 2')

    !
    ! If data features are given explicitly, they overwrite the shopping list
    !
    if(present(data_time_features))then
!      time_accept =  fu_time_in_list(fu_valid_time(id), zero_interval, list, data_time_features)
      time_accept =  fu_time_in_list(fu_valid_time(id), fu_validity_length(id), list, data_time_features)
    else
!      time_accept = fu_time_in_list(fu_valid_time(id), zero_interval, list)
      time_accept = fu_time_in_list(fu_valid_time(id), fu_validity_length(id), list)
    endif

!call msg('Accept 3')

    !
    ! If we decided to accept the field time but the input data are of some special temporal
    ! type, we might wish to reset the valid time of the input id to match the shopping list.
    ! However, we cannot put statis/monthly fields into the dynamic meteostack, i.e. into the
    ! multi-time stack. They must be in static single-time stack.
    !
    if(time_accept .and. present(time_to_force) .and. present(data_time_features))then
      select case(data_time_features)
        case(dynamic_map)
        case(monthly_climatology, static_climatology)
          if(list%earliest_valid_time == list%latest_valid_time)then
            time_to_force = list%earliest_valid_time
          else
            call set_error('Monthly/static field but two times to read. Dynamic meteo needed?', &
                         & 'fu_field_id_in_list')
            call set_error('Cannot put monthly/static field in dynamic stack','fu_field_id_in_list')
            return
          endif
        case default
          call set_error('Unknown data_time_features:' + fu_str(data_time_features), &
                       & 'fu_field_id_in_list')
          return
      end select
    endif
    !
    ! May be, we already accepted the variable ? If so, then time decides
    ! if the whole field is to be accepted
    !
    if(accept_field) then
      accept_field = time_accept
      return
    endif
    !
    ! If time is not acceptable, do not check anything else
    !
    if(.not. time_accept)then
      accept_field = .false.
      return
    endif
    !
    ! 3. Check source.
    !
    IF (.not.list%met_src == met_src_missing) THEN
      IF (.not.list%met_src == fu_met_src(id)) RETURN
    END IF
    !
    ! 4. Check quantity. Only known quantities are accepted.
    !
    IF(.not.fu_known_quantity(fu_quantity(id)))RETURN
    IF (list%quantities(1) /= accept_all_quantities) THEN
      IF (.NOT.fu_quantity_in_list(fu_quantity(id), list)) RETURN
    END IF

!call msg('Accept 5')


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
    if(.not.(list%grid == grid_missing))then
      !
      ! Get the grid and reposition it if needed
      !
      gridTmp = fu_grid(id)
      if(.not.fu_stdSilamGrid(gridTmp))then
        if(fu_ifLonGlobal(gridTmp))then
          if(.not.fu_stdSilamGrid(gridTmp)) call reposition_global_grid_lon(gridTmp)
        endif
        if(fu_ifLatGlobal(gridTmp))then
          if(.not.fu_stdSilamGrid(gridTmp)) call reposition_global_grid_lat(gridTmp)
        endif
        if(.not.fu_stdSilamGrid(gridTmp))return
      endif  ! input grid is non-standard

      if(.not.(fu_grids_arakawa_correspond(list%grid,gridTmp)))then
        if(.not.(fu_if_grid_covered(list%grid,gridTmp)))return
      endif
    end if


!call msg('Accept 6')


    !----------------------------------------
    !
    ! 5. Check levels.
    !
    SELECT CASE (list%level_indicator)

      ! --------------------------------
      !
      ! 5.1. One level case.

      CASE (one_level_only)

      IF (.NOT.fu_cmp_levs_eq(fu_level(id), list%floor_level)) RETURN

      ! --------------------------------
      !
      ! 5.2. Two level case.

      CASE (between_levels)

      lower_accept = .false.
      upper_accept = .false.

      lower_type = fu_leveltype(list%floor_level)  
      upper_type = fu_leveltype(list%ceiling_level)

      level = fu_level(id)
      level_type = fu_leveltype(level)


      lower: IF (level_type == lower_type) THEN
        IF (level >= list%floor_level) lower_accept = .true.
      ELSE IF ((lower_type == constant_altitude).or.&
	  & (lower_type == constant_height)) THEN
        lower_accept = .true.
      ELSE
        ! There's no way of telling which is higher, so we just take:
        CALL msg_warning('taking data from level of unknown height low')
        lower_accept = .true.
      END IF lower


      upper: IF (level_type == upper_type) THEN
        IF (level <= list%ceiling_level) upper_accept = .true.
      ELSE IF (fu_cmp_levs_eq(level, level_10m_above_ground).or.&
             & fu_cmp_levs_eq(level, level_2m_above_ground).or.&
             & fu_cmp_levs_eq(level, ground_level).or.&
             & fu_cmp_levs_eq(level, mean_sea_level)) THEN
        ! we'd better accept data from low ground levels, if quantity
        ! was on the list:
        upper_accept = .true.

      ELSE
        ! There's no way of telling which is higher, so we just take:
        CALL msg_warning('taking data from level of unknown height up')
        PRINT *, 'Level in check:'
        CALL report(level)
        PRINT *, 'Upper level in shoppig list:' 
        CALL report(list%ceiling_level)
        upper_accept = .true.
      END IF upper


      ! If either wrong, exit
      IF (.NOT.(lower_accept .and. upper_accept)) RETURN

      !
      ! 5.3. All levels case.

      CASE (accept_all_levels)   ! accepted in all cases

    END SELECT


!call msg('Accept 7')


    !
    ! 6. Accept this field-id and exit.
    !
    accept_field = .true.

  END FUNCTION fu_field_id_in_list


  !****************************************************************

  LOGICAL FUNCTION fu_time_in_list(time, validity_length_, list, data_time_features)

    ! Returns a true value, if the given field-id indicates a
    ! necessary field according to the shopping-list.
    !
    ! There are two ways of belonging to the list - belonging to the list of
    ! variables or satisfaction of the criteria like quantity, height range, 
    ! etc.
    ! The algorithm is the following. If the id belongs to the list of
    ! variables - only time range is checked. If the list of variables 
    ! does not accept the id - it is still checked for other general
    ! criteria. So, EITHER list of variables + time, OR general criteria
    ! + time may allow the field to come.
    !
    ! The reason for such tricky checking is purely historical. In principle,
    ! grib_analysis routines are deemed to describe complete set of 
    ! variables to shop. However, originally this technique was not used, so
    ! it will take time before it is eliminated from everywhere.
    !
    IMPLICIT NONE

    ! Imported parameters with intent(in):
    TYPE(silja_time), INTENT(in) :: time
    type(silja_interval), intent(in) :: validity_length_
    TYPE(silja_shopping_list), INTENT(in) :: list
    integer, intent(in), optional :: data_time_features

    ! local variable
    type(silja_interval) :: validity_length
    integer :: time_selector

    !
    ! 1. Some checking first.
    !
    IF (.NOT.defined(time)) THEN
      CALL set_error('time not defined','fu_time_in_list')
      RETURN
    END IF

    IF (.NOT.defined(list)) THEN
      CALL set_error('shopping-list not defined','fu_time_in_list')
      RETURN
    END IF
    !
    ! Validity_length can be undefined, then it is assumed zero
    !
    if(defined(validity_length_))then
      validity_length = validity_length_
    else
      validity_length = zero_interval
    endif
    !
    ! If data features are given explicitly, they overwrite the shopping list
    !
    if(present(data_time_features))then
      select case(data_time_features)
        case(dynamic_map)
          time_selector = list%time_indicator
        case(monthly_climatology)
          time_selector = accept_same_month
        case(static_climatology)
          time_selector = accept_all_times
        case default
          call set_error('Unknown data_time_features:' + fu_str(data_time_features), &
                       & 'fu_time_in_list')
          return
      end select
    else
      time_selector = list%time_indicator
    endif

    !
    ! Check times (if time_indicator == accept_all_times, then
    ! nothing is done here.
    !
    SELECT CASE (time_selector)   ! list%time_indicator)

      CASE (one_time_only)
        !
        ! earliest_valid_time should be inside the validity period
        !
!        fu_time_in_list = list%earliest_valid_time == time
        fu_time_in_list = fu_between_times(list%earliest_valid_time, & ! time
                                         & time, &                     ! lim 1
                                         & time+validity_length, &     ! lim2
                                         & .true.)             ! if accept boundaries

      CASE (between_times)
        !
        ! validity range should overlap with the earliest and latest range
        !
!        fu_time_in_list = time <= list%latest_valid_time .and. &
!                        & time >= list%earliest_valid_time
        fu_time_in_list = time <= list%latest_valid_time .and. &
                        & time + validity_length >= list%earliest_valid_time

      case(accept_all_times)
        fu_time_in_list = .true.

      case(accept_same_month)
        !
        ! months inside the validity period should be inside the months requries
        !
!        fu_time_in_list = fu_mon(time) <= fu_mon(list%latest_valid_time) .and. &
!                        & fu_mon(time) >= fu_mon(list%earliest_valid_time)
        fu_time_in_list = fu_mon(time) <= fu_mon(list%latest_valid_time) .and. &
                        & fu_mon(time+validity_length) >= fu_mon(list%earliest_valid_time)

      case default
        call set_error('Unknown time indicator','fu_time_in_list')
        fu_time_in_list = .false.
        return

    END SELECT

  end FUNCTION fu_time_in_list


  !****************************************************************

  logical function fu_field_id_in_list_of_vars(id, list, indexId)
    !
    ! Checks if the given field id is in the list of shopping
    ! variables. Should be called only if the logical switch
    ! useVariables is true. Otherwise - use the old way shown in
    ! fu_field_id_in_list
    !
    IMPLICIT NONE

    ! Imported parameters with intent(in):
    TYPE(silja_field_id), INTENT(in) :: id
    TYPE(silja_shopping_list), INTENT(in) :: list
    integer, intent(out), optional :: indexId

    ! Local variables
    integer :: iV
    type(silja_grid) :: gridTmp

    IF (.NOT.defined(id)) THEN
      CALL set_error('field-id not defined','fu_field_id_in_list_of_vars')
      RETURN
    END IF

    IF (.NOT.defined(list)) THEN
      CALL set_error('shopping-list not defined','fu_field_id_in_list_of_vars')
      RETURN
    END IF
    fu_field_id_in_list_of_vars = .false.
    if(present(indexId))then
      indexId = int_missing
    endif
    !
    ! The input grid might be non-standard for SILAM.
    ! In particular, it can be global with (0:360) coverage instead of -180,180.
    ! So far, a clumsy way to handle this is to use temporary standard grid,
    ! analogous to the input one. Later, the grid adjustment will be done
    ! if the field is accepted and stored into the stack.
    !
    ! So, get the grid and reposition it if needed
    !
    gridTmp = fu_grid(id)
    if(.not.fu_stdSilamGrid(gridTmp))then
      if(fu_ifLonGlobal(gridTmp))then
        if(.not.fu_stdSilamGrid(gridTmp)) call reposition_global_grid_lon(gridTmp)
      endif
      if(fu_ifLatGlobal(gridTmp))then
        if(.not.fu_stdSilamGrid(gridTmp)) call reposition_global_grid_lat(gridTmp)
      endif
      if(.not.fu_stdSilamGrid(gridTmp))return
    endif  ! input grid is non-standard

    !
    ! Variables are checked one-by-one
    !
    do iV = 1, size(list%vars)
      if(.not.defined(list%vars(iV))) return ! List is over...

      !  Quantity must correspond (unless accept all)
      !
      if(list%vars(iV)%quantity /= accept_all_quantities) then
        if(fu_quantity(id) /= list%vars(iV)%quantity) cycle
      end if

      ! Grid must Arakawa-correspond and have proper coverage unless accept all
      !
      if(.not.(list%vars(iV)%grid == grid_missing))then
        if(.not.fu_grids_arakawa_correspond(list%vars(iV)%grid, gridTmp))cycle
        if(.not.fu_if_grid_covered(list%vars(iV)%grid, gridTmp)) cycle
      end if

      !  Level must belong to vertial unless accept all verticals
      !
      if(defined(list%vars(iV)%vertical))then
        if(.not.fu_level_belongs_to_vertical(fu_level(id), list%vars(iV)%vertical)) cycle
      end if

      ! Meteo data source must correspond unless accept all sources
      !
      if(.not.(list%vars(iV)%mds == met_src_missing)) then 
         if(.not.(list%vars(iV)%mds == fu_met_src(id)))cycle
      end if

      !
      ! If substance name and mode value are present, they also must be checked.
      !
      if(trim(list%vars(iV)%chCocktailNm) /= '')then
         if(.not.(list%vars(iV)%chCocktailNm == fu_cocktail_name(id)))cycle        
      endif

      if(defined(list%vars(iV)%species))then
        if(.not. (list%vars(iV)%species == fu_species(id))) cycle
      else
        if(defined(fu_species(id))) cycle
      endif

      ! All is passed => id is in the list of variables
      !
      fu_field_id_in_list_of_vars = .true.
      if(present(indexId))then
        indexId = iV
      endif

      return

    end do   ! over variables

  end function fu_field_id_in_list_of_vars


  ! ***************************************************************

  LOGICAL FUNCTION fu_quantity_in_list(quantity, list)

    ! Description:
    ! Returns true value if the given quantity is on the list, or
    ! would be otherwise accepted.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_shopping_list), INTENT(in) :: list
    INTEGER, INTENT(in) :: quantity

    ! Local declarations:
    INTEGER :: i

    fu_quantity_in_list = fu_quantity_in_quantities(quantity, list%quantities)
    
    
  END FUNCTION fu_quantity_in_list
  

  !**************************************************************************

  function fu_mlev_quantity_in_list(list, quantity, ifVarsToo) result(ifMLev)
    !
    ! Checks if this quantity exists in the list and has set if2D flag
    !
    implicit none

    ! Return value
    type(silja_logical) :: ifMLev

    ! Imported parameters 
    TYPE(silja_shopping_list), INTENT(in) :: list
    INTEGER, INTENT(in) :: quantity
    logical, INTENT(in), optional :: ifVarsToo

    ! Local variables
    integer :: i

    ifMLev = silja_undefined

    IF (defined(list)) THEN
      do i=1,size(list%quantities)
        if(list%quantities(i) == int_missing)cycle
        if(list%quantities(i) /= quantity) cycle
        if(list%if2D(i) == silja_true)then
          ifMLev = silja_false
          return
        elseif(list%if2D(i) == silja_false)then
          ifMLev = silja_true
          return
        endif
      end do
    ELSE
      return
    END IF

    if(present(ifVarsToo))then
      if(ifVarsToo)then
        do i=1,size(list%vars)
          if(list%vars(i)%defined == silja_true)then
            if(list%vars(i)%quantity /= quantity)cycle
            if(list%vars(i)%vertLevNbr < int_missing)then
              ifMLev = silja_true
              return
            else
              ifMLev = silja_false
              return
            endif
          endif  ! if vars(i)%defined
        end do ! cycle through list%vars
      endif  ! ifVarsToo
    endif

  end function fu_mlev_quantity_in_list




  ! ****************************************************************
  ! ****************************************************************
  ! 
  !
  !       Private functions and subroutines.
  !
  !
  ! ***************************************************************
  ! ***************************************************************


  LOGICAL FUNCTION fu_shopping_list_defined(shopping_list)
    
    ! Description:
    ! Returns a true value, if the identificaion is defined by
    ! setting correct values to it using the set-function.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    
    ! Imported parameters with intent IN:
    TYPE(silja_shopping_list), INTENT(in) :: shopping_list
    
    fu_shopping_list_defined = fu_true(shopping_list%defined)
    
  END FUNCTION fu_shopping_list_defined


  ! ***************************************************************

  LOGICAL FUNCTION fu_shopping_variable_defined(shopping_var)
    
    ! Returns a true value, if the identificaion is defined by
    ! setting correct values to it using the set-function.
    ! 
    IMPLICIT NONE
    
    ! Imported parameters with intent IN:
    TYPE(silam_shopping_variable), INTENT(in) :: shopping_var
    
    fu_shopping_variable_defined = fu_true(shopping_var%defined)
    
  END FUNCTION fu_shopping_variable_defined


  ! ***************************************************************

  FUNCTION fu_met_src_of_shopping(shopping_list) result(met_src)
    ! 
    ! Returns the data met_src in the shopping list.
    !
    IMPLICIT NONE
 
    type(meteo_data_source) :: met_src
    ! Imported parameters with intent(in):
    TYPE(silja_shopping_list), INTENT(in) :: shopping_list
 
    IF (.NOT.defined(shopping_list)) THEN
      CALL set_error('undefined shopping list','fu_met_src_of_shopping')
      RETURN
    END IF

!    IF (shopping_list%met_src == met_src_missing) THEN
!      CALL msg_warning('no defined met_source, accepting all','fu_met_src_of_shopping')
!    END IF
    met_src = shopping_list%met_src

   END FUNCTION fu_met_src_of_shopping
  

  ! ***************************************************************

  FUNCTION fu_shlist_start_time(shopping_list)
    !
    ! Note that if all times accepted, the starting time is far in the past
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_time) :: fu_shlist_start_time
    
    ! Imported parameters with intent(in):
    TYPE(silja_shopping_list), INTENT(in) :: shopping_list

    if(shopping_list%time_indicator == accept_all_times)then
      fu_shlist_start_time = really_far_in_past
    else
      fu_shlist_start_time = shopping_list%earliest_valid_time
    endif

  END FUNCTION fu_shlist_start_time


  ! ***************************************************************

  FUNCTION fu_shlist_end_time(shopping_list)
    !
    ! Note that if all times are accepted, far-future time is set
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_time) :: fu_shlist_end_time
    
    ! Imported parameters with intent(in):
    TYPE(silja_shopping_list), INTENT(in) :: shopping_list
    
    if(shopping_list%time_indicator == accept_all_times)then
      fu_shlist_end_time = really_far_in_future
    else
      fu_shlist_end_time = shopping_list%latest_valid_time
    endif
    
  END FUNCTION fu_shlist_end_time


  ! ***************************************************************

  integer function fu_nbr_of_quantities_in_list(list, ifVarsToo)
    !
    ! Returns a number of reasonable quantities in the shopping list
    !
    IMPLICIT NONE
    
    ! Imported parameters with intent(in):
    TYPE(silja_shopping_list), INTENT(in) :: list
    logical, intent(in), optional :: ifVarsToo

    ! Local variables
    integer :: i

    fu_nbr_of_quantities_in_list = 0
    IF (defined(list)) THEN
      do i=1,size(list%quantities)
        if(list%quantities(i) == int_missing)cycle
        fu_nbr_of_quantities_in_list = fu_nbr_of_quantities_in_list + 1
      end do
    ELSE
      return
    END IF

    if(present(ifVarsToo))then
      if(ifVarsToo)then
        do i=1,size(list%vars)
          if(list%vars(i)%defined == silja_true)then
            fu_nbr_of_quantities_in_list = fu_nbr_of_quantities_in_list + 1
          endif
        end do ! cycle through list%vars
      endif  ! ifVarsToo
    endif

  end function fu_nbr_of_quantities_in_list


  ! ***************************************************************

  integer FUNCTION fu_quantity_from_shopping_list(shopping_list, indexQ)
    
    IMPLICIT NONE
    
    ! Imported parameters with intent(in):
    TYPE(silja_shopping_list), INTENT(in) :: shopping_list
    integer, intent(in) :: indexQ
    
    if(indexQ < 0 .or. indexQ > size(shopping_list%quantities))then
      print *, 'Given index: ', indexQ
      call set_error('Strange idnex','fu_quantity_from_shopping_list')
      fu_quantity_from_shopping_list = int_missing
      return
    endif
    fu_quantity_from_shopping_list = shopping_list%quantities(indexQ)

  END FUNCTION fu_quantity_from_shopping_list

  ! ***************************************************************

  integer FUNCTION fu_quantity_of_shopping_var(shopVar)
    
    IMPLICIT NONE
    
    ! Imported parameters with intent(in):
    type(silam_shopping_variable), intent(in) :: shopVar

  
    fu_quantity_of_shopping_var = shopVar%quantity

  END FUNCTION fu_quantity_of_shopping_var

  ! ***************************************************************

  FUNCTION fu_all_quantities_from_list(list, ifVarsToo) result(Qs)
    !
    ! Returns the complete list of quantities mentioned in the list. 
    ! If ifVarsToo == .true. then this list also includes all quantities
    ! mentioned in list%vars.
    ! This is not exactly correct procedure because list%vars are more limited
    ! than general list%quantities. But sometimes we need to know what in 
    ! principle we can get from the list
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    INTEGER, DIMENSION(max_quantities) :: Qs
    
    ! Imported parameters with intent(in):
    TYPE(silja_shopping_list), INTENT(in) :: list
    logical, intent(in), optional :: ifVarsToo

    ! Local variables
    integer :: i, iCount
  
    Qs = int_missing
    IF (defined(list)) THEN
      iCount=1
      do i=1,size(list%quantities)
        if(list%quantities(i) == int_missing)cycle
        Qs(iCount) = list%quantities(i)
        iCount = iCount + 1
      end do
    ELSE
      return
    END IF

    if(present(ifVarsToo))then
      if(ifVarsToo)then
        do i=1,size(list%vars)
          if(list%vars(i)%defined == silja_true)then
            Qs(iCount) = list%vars(i)%quantity
            iCount = iCount + 1
          endif  ! if vars(i)%defined
        end do ! cycle through list%vars
      endif  ! ifVarsToo
    endif

  END FUNCTION fu_all_quantities_from_list


  ! ***************************************************************

  integer FUNCTION fu_request_from_shopping_list(shopping_list, indexQ)
    
    IMPLICIT NONE
    
    ! Imported parameters with intent(in):
    TYPE(silja_shopping_list), INTENT(in) :: shopping_list
    integer, intent(in) :: indexQ
    
    if(indexQ < 0 .or. indexQ > size(shopping_list%quantities))then
      call set_error('Strange idnex','fu_quantity_from_shopping_list')
      fu_request_from_shopping_list = 0
      return
    endif
    fu_request_from_shopping_list = shopping_list%request(indexQ)

  END FUNCTION fu_request_from_shopping_list


  ! ***************************************************************

  FUNCTION fu_all_requests_from_list(list, ifVarsToo) result(RQs)
    !
    ! Returns the complete list of quantities mentioned in the list. 
    ! If ifVarsToo == .true. then this list also includes all quantities
    ! mentioned in list%vars.
    ! This is not exactly correct procedure because list%vars are more limited
    ! than general list%quantities. But sometimes we need to know what in 
    ! principle we can get from the list
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    INTEGER, DIMENSION(max_quantities) :: Rqs
    
    ! Imported parameters with intent(in):
    TYPE(silja_shopping_list), INTENT(in) :: list
    logical, intent(in), optional :: ifVarsToo

    ! Local variables
    integer :: i, iCount
  
    RQs = 0
    IF (defined(list)) THEN
      iCount=1
      do i=1,size(list%request)
!        if(list%request(i) == 0)cycle
        if(list%quantities(i) == int_missing)cycle ! synchronization with quantities
        Rqs(iCount) = list%request(i)
        iCount = iCount + 1
      end do
    ELSE
      return
    END IF

    if(present(ifVarsToo))then
      if(ifVarsToo)then
        do i=1,size(list%vars)
          if(list%vars(i)%defined == silja_true)then
            RQs(iCount) = list%vars(i)%request
            iCount = iCount + 1
          endif  ! if vars(i)%defined
        end do ! cycle through list%vars
      endif  ! ifVarsToo
    endif

  END FUNCTION fu_all_requests_from_list


  !*****************************************************************

  function fu_if2D_of_shopping_quantity(shopping_list, indexQ)

    IMPLICIT NONE
    
    ! return value
    type(silja_logical) :: fu_if2D_of_shopping_quantity

    ! Imported parameters with intent(in):
    TYPE(silja_shopping_list), INTENT(in) :: shopping_list
    integer, intent(in) :: indexQ
    
    if(indexQ < 0 .or. indexQ > size(shopping_list%quantities))then
      call set_error('Strange idnex','fu_quantity_from_shopping_list')
      fu_if2D_of_shopping_quantity = silja_undefined
      return
    endif
    fu_if2D_of_shopping_quantity= shopping_list%if2D(indexQ)

  end function fu_if2D_of_shopping_quantity


  ! ***************************************************************

  FUNCTION fu_all_if2D_from_shopping_list(list, ifVarsToo) result(if2Ds)
    !
    ! Returns the complete list of quantities mentioned in the list. 
    ! If ifVarsToo == .true. then this list also includes all quantities
    ! mentioned in list%vars.
    ! This is not exactly correct procedure because list%vars are more limited
    ! than general list%quantities. But sometimes we need to know what in 
    ! principle we can get from the list
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    type(silja_logical), DIMENSION(max_quantities) :: if2Ds
    
    ! Imported parameters with intent(in):
    TYPE(silja_shopping_list), INTENT(in) :: list
    logical, intent(in), optional :: ifVarsToo

    ! Local variables
    integer :: i, iCount
  
    if2Ds = silja_undefined
    IF (defined(list)) THEN
      iCount=1
      do i=1,size(list%request)
        if(list%quantities(i) == int_missing)cycle ! Synchronize with quantities
!        if(list%request(i) == 0)cycle
        if2Ds(iCount) = list%if2D(i)
        iCount = iCount + 1
      end do
    ELSE
      return
    END IF

    if(present(ifVarsToo))then
      if(ifVarsToo)then
        do i=1,list%nVars
          if(list%vars(i)%defined == silja_true)then
            if(list%vars(i)%vertLevNbr /= int_missing)then
              if2Ds(iCount) = silja_true
            else
              if2Ds(iCount) = silja_false
            endif
            iCount = iCount + 1
          endif  ! if vars(i)%defined
        end do ! cycle through list%vars
      endif  ! ifVarsToo
    endif

  END FUNCTION fu_all_if2D_from_shopping_list


  ! ***************************************************************

  subroutine set_request_in_shopping_list(shopping_list, qIndex, new_value)
    
    IMPLICIT NONE
    
    ! Imported parameters with intent(in):
    TYPE(silja_shopping_list), INTENT(inout) :: shopping_list
    integer, intent(in) :: qIndex, new_value

    
    if(qIndex < 0 .or. qIndex > size(shopping_list%quantities))then
      call set_error('Strange idnex','fu_quantity_from_shopping_list')
      return
    endif
    if(new_value < 0 .or. new_value > 2)then
      print *, 'New request value: ', new_value
      call set_error('Strange new request value','set_request_in_shopping_list')
      return
    endif

    shopping_list%request(qIndex) = new_value

  END subroutine set_request_in_shopping_list


  ! ***************************************************************

  subroutine set_list_time_indicator(shopping_list, fValue)

    IMPLICIT NONE
    ! Imported parameters with intent(in):
    TYPE(silja_shopping_list), INTENT(inout) :: shopping_list
    integer, intent(in) :: fValue

    shopping_list%time_indicator = fValue

  END subroutine set_list_time_indicator


  ! ***************************************************************

  integer function fu_list_time_indicator(shopping_list)

    IMPLICIT NONE

    TYPE(silja_shopping_list), INTENT(in) :: shopping_list

    fu_list_time_indicator = shopping_list%time_indicator

  end function fu_list_time_indicator


  !************************************************************************

  subroutine replace_quantity_in_shop_list(shop_lst, qOld, qNew, ifVarsToo)
    !
    ! Searches through the whole list and all cases of the given old quantity
    ! replaces it with the new one. No other parameters are touched
    !
    implicit none

    ! Imported parameters
    TYPE(silja_shopping_list), INTENT(inout) :: shop_lst
    integer, intent(in) :: qOld, qNew
    logical, intent(in) :: ifVarsToo

    ! Local variables
    integer :: i

    IF (defined(shop_lst)) THEN
      do i=1,size(shop_lst%quantities)
        if(shop_lst%quantities(i) == int_missing)cycle ! Synchronize with quantities
        if(shop_lst%quantities(i) == qOld)then
          shop_lst%quantities(i) = qNew
        endif
      end do
    ELSE
      return
    END IF

    if(ifVarsToo)then
      do i=1,shop_lst%nVars
        if(shop_lst%vars(i)%defined == silja_true)then
          if(shop_lst%vars(i)%quantity == qOld) shop_lst%vars(i)%quantity = qNew
        endif  ! if vars(i)%defined
      end do ! cycle through list%vars
    endif  ! ifVarsToo

  end subroutine replace_quantity_in_shop_list


  !********************************************************

  logical function fu_compare_variables_eq(var1, var2)
    !
    ! Checks if variables are equal.
    !
    implicit none

    ! Imported parameters with intent IN
    type(silam_shopping_variable), intent(in) :: var1, var2

    fu_compare_variables_eq = .false.
    !
    ! Undefined variables are always equal to each other...
    !
    if(defined(var1))then
      if(.not.defined(var2)) return
    else
      fu_compare_variables_eq = .not.defined(var2)
      return
    end if

    if(var1%quantity /= var2%quantity) return

    if(var1%request /= var2%request) return

    if(.not. fu_cmp_verts_eq(var1%vertical, var2%vertical)) return

    if(.not. var1%grid == var2%grid) return

    if(.not. var1%mds == var2%mds) return

    if(.not. (var1%species == var2%species))return

    if(.not. trim(var1%chCocktailNm) == trim(var2%chCocktailNm))return

    fu_compare_variables_eq = .true.

  end function fu_compare_variables_eq
 
 

  ! ***************************************************************
  ! ***************************************************************
  ! 
  !
  !                 TESTS
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  subroutine report_shopping_variable (var)
    !
    ! Prints the report of the shopping variable
    !
    implicit none

    ! Imported parameters with intent IN
    type(silam_shopping_variable), intent (in) :: var

    if(.not.fu_true(var%defined))then
      call msg('Undefined variable')
      return
    end if

    if(var%request == 0)then
      call msg(fu_connect_strings('Quantity (not needed)',fu_quantity_string(var%quantity)))
    elseif(var%request == 1)then
      call msg(fu_connect_strings('Quantity (desirable)',fu_quantity_string(var%quantity)))
    elseif(var%request == 2)then
      call msg(fu_connect_strings('Quantity (mandatory)',fu_quantity_string(var%quantity)))
    else
      call msg(fu_connect_strings('Quantity (necessity unknown)',fu_quantity_string(var%quantity)))
    endif

    call msg('Species:' + fu_str(var%species))

    call report(var%vertical)

    call report(var%grid)

    call report(var%mds)

  end subroutine report_shopping_variable



  ! ***************************************************************

  SUBROUTINE print_shopping_report(list)

    ! Description:
    ! Print contents of a shoppiog list for test purposes.
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silja_shopping_list), INTENT(in) :: list

    INTEGER :: iUnit, i
    INTEGER, DIMENSION(2) :: fUnits

    fUnits(1:2) = (/6, run_log_funit/)

    IF (.NOT.fu_shopping_list_defined(list)) THEN
      call msg(' Undefined shopping list. ')
      RETURN
    END IF

!!!$    IF (list%take_everything) THEN
!!!$      PRINT *, '*************************************'
!!!$      PRINT *, 'Take everything!'
!!!$      PRINT *, '*************************************'
!!!$      RETURN
!!!$    END IF

    call msg('************************shopping list ******')
    call msg('Variables in the list:')


   

    do i = 1, size(list%vars)
      if(defined(list%vars(i))) then 
         call report(list%vars(i))
!        PRINT *, fu_quantity_string(list%vars(i)%quantity), 'Necessity: ',list%vars(i)%request
!        write(run_log_funit,*) fu_quantity_string(list%vars(i)%quantity), 'Necessity: ',list%vars(i)%request
      endif
    end do


       SELECT CASE (list%time_indicator)

          CASE (one_time_only)
          call msg(fu_connect_strings('single time: ', fu_str(list%earliest_valid_time)))

          CASE (between_times)
          do iUnit = 1,2
             if (smpi_global_rank /= 0 .and. iUnit==1) cycle !Be quiet at stdut
            write(funits(iUnit),'(A,A,A)')'between times: ',&
                 & fu_str(list%earliest_valid_time),&
                 & fu_str(list%latest_valid_time)
           enddo

          CASE (accept_all_times)
          call msg('all times accepted')

       END SELECT

       IF ((.NOT.defined(list%floor_level)).and.&
           & (.NOT.defined(list%ceiling_level))) THEN
         call msg('no vertical level boundaries')
       ELSE
         call msg('lower and upper levels:')
         CALL report(list%floor_level)
         CALL report(list%ceiling_level)
       END IF

       call msg('Quantities we shop:')
       DO i = 1, SIZE(list%quantities)
         IF (fu_known_quantity(list%quantities(i))) THEN
          do iUnit = 1,2
             if (smpi_global_rank /= 0 .and. iUnit==1) cycle !Be quiet at stdut
           write(funits(iUnit),'(A,X,A,X,I2)')fu_quantity_string(list%quantities(i)), 'Necessity:',list%request(i)
          enddo !iUnit loop
         ELSE
           EXIT
         END IF
       END DO


    call msg('*************************************')

  END SUBROUTINE print_shopping_report

END MODULE shopping_lists


