MODULE fields_3d ! 3d scalar fields 

  ! Description:
  ! This module contains definition and tools for 
  ! 3-dimensional scalar fields (x,y,z) of meteorological data. For
  ! this type we define: silja_3d_field contains multilevel data on
  ! several vertical levels, whose leveltypes are same.
  !
  ! Author: Mika Salonoja, FMI email Mika.Salonoja@fmi.fi
  !
  ! All units: SI
  ! 
  ! All degrees: real numbers in degrees (degrees and parts per
  ! hundred/parts per thousand , NOT in degrees and minutes
  ! 
  ! Language: ANSI standard Fortran 90
  ! 
  ! Modules used:
!  USE globals
  USE fields
  !$ use OMP_LIB
  !USE levels

  IMPLICIT NONE

  ! The public functions and subroutines available in this module:

  ! Change contents of 3d-field:
  PUBLIC add_field_to_3d_field
  PUBLIC set_3d_field_empty
  PUBLIC organize_fields_vertically
  PUBLIC add_surf_pre_to_3d_field ! to a hybrid-level 3d-field

  ! Operations with 3d field
  PUBLIC fu_interpolate_to_position ! interpolate in x, y and z
  public ddz_of_field_3d      ! vertical metrical derivative d/dz
  public check_field_3d

  ! Return stuff from 3d-field:
  PUBLIC fu_field_from_3d_field ! i:th field from bottom
  public create_field_in_field_3d ! creates the next 2D field (note: CREATES)
  PUBLIC fu_grid_data_from_3d   ! i:th 2D grid-data
  public get_column_from_field_3d ! for the given grid index
  PUBLIC fu_find_level          ! find 2d-field by its level
  PUBLIC fu_lowest_field        ! field closest to ground
  PUBLIC fu_highest_field       ! field highest from ground
  PUBLIC fu_closest_field       ! up- or downwards
  PUBLIC fu_surface_pressure_field ! field in 3d-field
  PUBLIC vertical_levels        ! what 2D-levels there is in 3D-field?
  PUBLIC pressure_on_level      ! pr of a level of 3D field if every gridpoint

  ! Administrative data found in 3d-field:
  PUBLIC fu_grid
  PUBLIC fu_size
  PUBLIC fu_valid_time
  PUBLIC fu_analysis_time
  PUBLIC fu_forecast_length
  PUBLIC fu_met_src
  PUBLIC fu_quantity
  PUBLIC fu_leveltype
  PUBLIC fu_number_of_fields
  PUBLIC fu_pressures
  public fu_substance_name
  public fu_mode
  public fu_optical_wave_length
  public set_3dField_params_from_fieldId
  public set_defined
  PUBLIC report ! for testing
  public set_valid_time
  public defined
  public fu_species

  ! The private functions and subroutines not to be used elsewhere:
  PRIVATE find_closest_fields ! vertically from a position
  PRIVATE find_closest_fields_2 ! vertically from a gridpoint
  PRIVATE fu_field_3d_defined
  PRIVATE fu_met_src_of_3d_fi
  PRIVATE fu_grid_of_3d_fi
  PRIVATE fu_valid_time_of_3d_fi
  PRIVATE fu_quantity_of_3d_fi
  PRIVATE fu_leveltype_of_3d_fi
  PRIVATE fu_3d_fieldsize
  PRIVATE fu_find_level_scalar
  PRIVATE fu_lowest_field_scalar
  PRIVATE fu_highest_field_scalar
  PRIVATE fu_pressures_scalar
  PRIVATE fu_level_pressure_in_3d
  PRIVATE fu_surf_pre_field_scalar
  PRIVATE print_field_3d_report
  PRIVATE fu_analysis_time_of_3d
  PRIVATE fu_forecast_length_of_3d
  private fu_SubstNm_of_3d_fld
  private fu_mode_of_3d_fld
  private fu_optic_wavelen_of_3d_fld
  private set_field_3d_defined
  private set_valid_time_of_3d_fi
  private fu_closest_field_scalar
  private fu_species_of_3d_fld


  ! Generic names and operator-interfaces of some functions:
  INTERFACE defined
    MODULE PROCEDURE fu_field_3d_defined
  END INTERFACE

  interface set_defined
    module procedure set_field_3d_defined
  end interface

  INTERFACE report
    MODULE PROCEDURE print_field_3d_report
  END INTERFACE

  INTERFACE fu_size
    MODULE PROCEDURE fu_3d_fieldsize
  END INTERFACE

  INTERFACE fu_find_level
    MODULE PROCEDURE fu_find_level_scalar
  END INTERFACE

  INTERFACE fu_met_src  
    MODULE PROCEDURE fu_met_src_of_3d_fi
  END INTERFACE

  INTERFACE fu_grid 
    MODULE PROCEDURE fu_grid_of_3d_fi
  END INTERFACE

  INTERFACE fu_valid_time
    MODULE PROCEDURE fu_valid_time_of_3d_fi
  END INTERFACE

  INTERFACE fu_analysis_time
    MODULE PROCEDURE fu_analysis_time_of_3d
  END INTERFACE

  INTERFACE fu_forecast_length
    MODULE PROCEDURE fu_forecast_length_of_3d
  END INTERFACE

  INTERFACE fu_quantity
    MODULE PROCEDURE fu_quantity_of_3d_fi
  END INTERFACE

  INTERFACE fu_pressures
    MODULE PROCEDURE fu_pressures_scalar
  END INTERFACE

  INTERFACE fu_leveltype
    MODULE PROCEDURE fu_leveltype_of_3d_fi
  END INTERFACE

  INTERFACE fu_lowest_field
    MODULE PROCEDURE fu_lowest_field_scalar
  END INTERFACE

  INTERFACE fu_highest_field
    MODULE PROCEDURE fu_highest_field_scalar
  END INTERFACE

  INTERFACE fu_closest_field
    MODULE PROCEDURE fu_closest_field_scalar
  END INTERFACE

  INTERFACE fu_surface_pressure_field
    MODULE PROCEDURE fu_surf_pre_field_scalar
  END INTERFACE

  interface fu_substance_name
    module procedure fu_SubstNm_of_3d_fld
  end interface

  interface fu_mode
    module procedure fu_mode_of_3d_fld
  end interface

  interface fu_optical_wave_length
    module procedure fu_optic_wavelen_of_3d_fld
  end interface

  interface set_valid_time
    module procedure set_valid_time_of_3d_fi
  end interface

  interface fu_species
     module procedure fu_species_of_3d_fld
  end interface

  ! Public types with private components defined in this module:
  TYPE silja_3d_field
    PRIVATE
    INTEGER :: quantity
    type(silam_species) :: species
    TYPE(silja_time) :: valid_time
    type(meteo_data_source) :: met_src
    TYPE(silja_grid) :: grid
    INTEGER :: leveltype ! leveltype same for all fields in 3d-field!!
    TYPE(silja_fieldpointer), DIMENSION(max_levels) :: fields
    INTEGER, DIMENSION(max_levels) :: vertical_order
    LOGICAL :: in_order
    TYPE(silja_field), POINTER :: surface_pressure
    INTEGER :: number_of_fields
    TYPE(silja_logical) :: defined
  END TYPE silja_3d_field

  ! Public types with public components defined in this module:
  TYPE silja_field_3d_pointer
    TYPE(silja_3d_field), POINTER :: fp
  END TYPE silja_field_3d_pointer

CONTAINS

  ! ***************************************************************

  SUBROUTINE add_field_to_3d_field(field, field_3d)

    ! Description:
    ! Adds a 2D-field to field_3d. 
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_field), POINTER :: field
    TYPE(silja_3d_field), INTENT(inout) :: field_3d

    ! Local declarations:
    INTEGER ::i

    !----------------------------------------
    !
    ! 1. Check status of field_3d and the field that should be added
    !    ----------------------------------------------------------

    IF (.NOT.defined(field)) THEN
      CALL set_error('cannot put undefined field to field_3d', 'add_field_to_3d_field')
     RETURN
    END IF


    new_or_old: IF (.NOT.defined(field_3d)) THEN

      !----------------------------------------
      !
      ! 2. A new, untouched field_3d
      !    -----------------------

      field_3d%fields(1)%fp => field
      field_3d%number_of_fields = 1
      field_3d%quantity = fu_quantity(field)
      IF (error) RETURN
      
      field_3d%species = fu_species(field)
      IF (error) RETURN

      field_3d%valid_time = fu_valid_time(field)
      IF (error) RETURN

      field_3d%met_src = fu_met_src(field)
      IF (error) RETURN

      field_3d%leveltype = fu_leveltype(fu_level(field))
      IF (error) RETURN

      ! It is assumed, that is a vertical field_3d of a single
      ! quantity, the horizontal grid is the same for all levels:
      field_3d%grid = fu_grid(field)
      IF (error) RETURN

      ! No surface pressure given yet:
      NULLIFY(field_3d%surface_pressure)

      field_3d%in_order = .false.

      ! Everything ok:
      field_3d%defined = fu_set_true()


    ELSE

      !----------------------------------------
      !
      ! 3. Already stuff in field_3d, check that it isn't already
      ! there, and check that the new field belogs here.
      !    -----------------------

      belongs: IF (fu_field_belogs_to_this_3d()) THEN


        !----------------------------------------
        !
        ! 4. Add to bottom
        !    -----------------------

        IF (field_3d%number_of_fields < max_levels) THEN
          field_3d%number_of_fields = field_3d%number_of_fields + 1
          field_3d%fields(field_3d%number_of_fields)%fp => field
          field_3d%in_order = .false.
        ELSE
          CALL set_error('cannot add any more field to this field_3d','add_field_to_3d_field')
        END IF

      END IF belongs

    END IF new_or_old

  CONTAINS

    LOGICAL FUNCTION fu_field_belogs_to_this_3d()

      ! Description:
      ! Returns true value, if the current field should be added to
      ! the current 3d-field.
      !
      ! Method:
      ! 1. Check that the field is not already in 3d-field
      ! 2. Check that quantity, time, met_src and leveltype match.
      !
      IMPLICIT NONE

      ! Local declarations:
      INTEGER :: i

      fu_field_belogs_to_this_3d = .false.

      ! Already in 3d-field?
      DO i = 1, field_3d%number_of_fields
        IF (fu_id(field) == fu_id(field_3d%fields(i)%fp)) RETURN
      END DO

      IF (fu_leveltype(fu_level(field)) /= field_3d%leveltype) THEN
        call msg('')
        call msg('Field that comes:')
        CALL report(field)
        call msg('3D field:')
        call report(field_3d)
        CALL msg_warning('cannot mix different leveltypes','add_field_to_3d_field')
        RETURN
      END IF

      IF (.NOT.(fu_grid(field) == field_3d%grid)) THEN
        CALL report(field)
        CALL set_error('cannot mix different grids','add_field_to_3d_field')
        RETURN
      END IF

      IF (fu_valid_time(field) /= field_3d%valid_time) RETURN
      IF (.not.fu_met_src(field) == field_3d%met_src) RETURN
      IF (fu_quantity(field) /= field_3d%quantity) RETURN
      IF (.not. (fu_species(field) == field_3d%species)) RETURN

      fu_field_belogs_to_this_3d = .true.

    END FUNCTION fu_field_belogs_to_this_3d

  END SUBROUTINE add_field_to_3d_field



  ! ***************************************************************

  SUBROUTINE organize_fields_vertically(field_3d)

    ! Description:
    ! This routine sets the field in 3d-field into vertical order,
    ! from lowest upwards.
    !
    ! Method:
    ! For each level the number of lowel levels (closer to ground) is
    ! calculated.
    ! Comparison result are stored in vector vertical_order inside 3d
    ! -field.
    ! 
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_3d_field), INTENT(inout) :: field_3d

    ! Local declarations:
    INTEGER :: i, j
    TYPE(silja_field), POINTER :: field
    INTEGER :: lower_count
    TYPE(silja_level) :: level


    !----------------------------------------
    !
    ! 1. Check status of field_3d
    !    ------------------------

    IF (.NOT.defined(field_3d)) THEN
      CALL set_error('undefined field_3d given','organize_fields_vertically')
      RETURN
    END IF

    IF (field_3d%in_order) RETURN

    !
    ! 2. Make the counts.
    !    ----------------

    field_3d%vertical_order = int_missing

    outer: DO i = 1, field_3d%number_of_fields

      lower_count = 0
      level = fu_level(field_3d%fields(i)%fp)

      inner: DO j = 1, field_3d%number_of_fields
        IF (i == j) CYCLE inner

        IF (fu_level(field_3d%fields(j)%fp) < level) lower_count = lower_count + 1

      END DO inner

      field_3d%vertical_order(lower_count+1) = i

    END DO outer

    field_3d%in_order = .true.

  END SUBROUTINE organize_fields_vertically



  ! ***************************************************************

  SUBROUTINE add_surf_pre_to_3d_field(field_3d, surface_pressure)

    ! Description:
    ! Adds a correnponding ground surface pressure field to a 3d
    ! -field and checks it.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_3d_field), POINTER :: field_3d
    TYPE(silja_field), POINTER :: surface_pressure

    IF (.NOT.defined(field_3d)) THEN
      CALL set_error('undefined field_3d given','add_surf_pre_to_3d_field')
      RETURN
    END IF

    IF (.NOT.defined(surface_pressure)) THEN
      CALL set_error('undefined surface_pressure given'&
         & ,'add_surf_pre_to_3d_field')
      RETURN
    END IF

    matching_check: IF (&
!        & (fu_met_src(surface_pressure) == field_3d%met_src).and.&
        & (fu_valid_time(surface_pressure) == field_3d%valid_time)&
        & .and.(fu_grids_match_closely(&
        & fu_grid(surface_pressure),field_3d%grid))) THEN

      field_3d%surface_pressure => surface_pressure

    ELSE

      CALL set_error('fields do not match','add_surf_pre_to_3d_field')
      CALL report(field_3d)
      CALL report(surface_pressure)

    END IF matching_check

  END SUBROUTINE add_surf_pre_to_3d_field


  ! ***************************************************************

  SUBROUTINE set_3d_field_empty(field_3d)

    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_3d_field), INTENT(inout) :: field_3d

    INTEGER :: i

    field_3d%defined = silja_false
    field_3d%number_of_fields = 0

    DO i = 1, SIZE(field_3d%fields)
      NULLIFY(field_3d%fields(i)%fp)
    END DO

    nullify(field_3d%surface_pressure)

  END SUBROUTINE set_3d_field_empty


  ! ***************************************************************

  SUBROUTINE set_3dField_params_from_fieldId(field3d, id)
    implicit none
    type(silja_field_id), intent(in) :: id
    type(silja_3d_field), pointer :: field3d
    integer :: i
    field3d%quantity = fu_quantity(id)
    field3d%species = fu_species(id)
    field3d%valid_time = fu_valid_time(id)
    field3d%met_src = fu_met_src(id)
    field3d%grid = fu_grid(id)
    field3d%leveltype = fu_leveltype(fu_level(id))
     
!    call msg('Number of fields in f3d: ', field3d%number_of_fields)
    do i = 1, field3d%number_of_fields
      field3d%vertical_order(i) = i
    enddo
    field3d%in_order = .true.
    field3d%defined = silja_true
  
  end SUBROUTINE set_3dField_params_from_fieldId


  !************************************************************************

  subroutine create_field_in_field_3d(field_3d, indVert, id_2d, fValues)
    !
    ! Creates a requested 2D field at a given level and sets its grid to fValue
    !
    implicit none

    ! Imported parameters
    TYPE(silja_3d_field), POINTER :: field_3d
    integer, intent(in) :: indVert
    type(silja_field_id), intent(in) :: id_2d
    real, dimension(:), pointer :: fValues

    ! Local variables
    integer :: iStat

    field_3d%number_of_fields = field_3d%number_of_fields + 1
    if(field_3d%number_of_fields > max_levels)then
      call set_error('3D field is full','create_field_in_field_3d')
      return
    endif
    allocate(field_3d%fields(field_3d%number_of_fields)%fp, stat = iStat)
    if(iStat /= 0)then
      call msg('Failed to allocate field in the 3D field:',field_3d%number_of_fields)
      call set_error('Failed to allocate field in the 3D field','create_field_in_field_3d')
      return
    endif
!    call report(id_2d)
    call set_field(id_2d, fValues, field_3d%fields(field_3d%number_of_fields)%fp, .true.)
!    call msg('done')

  end subroutine create_field_in_field_3d


  ! ***************************************************************

  FUNCTION fu_field_from_3d_field(field_3d, number, ifSkipOrder) 

    ! Description:
    ! Returns a pointer to the n:th field in a field_3d. Number
    ! refers to the number in order inside field: 1 is lowest and the
    ! highest number highest from ground.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_field), POINTER :: fu_field_from_3d_field

    ! Imported parameters with intent(in):
    TYPE(silja_3d_field), POINTER :: field_3d
    INTEGER, INTENT(in) :: number
    logical, intent(in), optional :: ifSkipOrder

    IF ((number>=0).and.(number <= field_3d%number_of_fields)) THEN

      if(present(ifSkipOrder))then
        if(ifSkipOrder) then
          fu_field_from_3d_field => field_3d%fields(number)%fp
          return
        endif
      endif

      IF (field_3d%in_order) THEN
        fu_field_from_3d_field => field_3d%fields(field_3d%vertical_order(number))%fp
      ELSE
        CALL set_error('this 3dfield not in order','fu_field_from_3d_field')
        NULLIFY(fu_field_from_3d_field)
      END IF

    ELSE
      call msg('Number of field', number)
      call msg('Number of fields in 3d',field_3d%number_of_fields)
      call msg("####################################################################")
      call report(field_3d)
      call msg("####################################################################")

      CALL set_error('this number field not in this 3d-field','fu_field_from_3d_field')
      NULLIFY(fu_field_from_3d_field)
    END IF

  END FUNCTION fu_field_from_3d_field


  ! ***************************************************************

  FUNCTION fu_grid_data_from_3d(field_3d, number) 

    ! Description:
    ! Returns a pointer to the n:th grid-data in a field_3d. Number
    ! refers to the number in order inside field: 1 is lowest and the
    ! highest number highest from ground.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    REAL, DIMENSION(:), POINTER :: fu_grid_data_from_3d

    ! Imported parameters with intent(in):
    TYPE(silja_3d_field), POINTER :: field_3d
    INTEGER, INTENT(in) :: number

    IF ((number>=0).and.(number <= field_3d%number_of_fields)) THEN

      IF (field_3d%in_order) THEN

        fu_grid_data_from_3d => fu_grid_data(field_3d%fields(field_3d%vertical_order(number))%fp)
     
      ELSE
       CALL set_error('this 3dfield not in order', 'fu_grid_data_from_3d')
       NULLIFY(fu_grid_data_from_3d)
      END IF

    ELSE

      call msg('This number field not in this 3d-field', number)
      call msg("########### Field 3D report ##################")
      call report(field_3d)
      call msg("########### Field 3D report end ##################")
      CALL set_error('This number field not in this 3d-field', 'fu_grid_data_from_3d')
      NULLIFY(fu_grid_data_from_3d)

    END IF

  END FUNCTION fu_grid_data_from_3d


  !*****************************************************************

  subroutine get_column_from_field_3d(field_3d, grid_index, arValues)
    !
    ! Fills-in the given array with the elements from the given 3D field
    !
    implicit none

    ! Imported parameters
    TYPE(silja_3d_field), POINTER :: field_3d
    INTEGER, INTENT(in) :: grid_index
    real, dimension(:), pointer :: arValues

    ! Local variables
    integer :: iLev
    REAL, DIMENSION(:), POINTER :: grid_data_from_3d

    if(grid_index < 0 .or. grid_index > fu_number_of_gridpoints(field_3d%grid))then
      call msg('grid_index given does not fit into below grid:',grid_index)
      call report(field_3d%grid)
      call set_error('Strange grid_index given','get_column_from_field_3d')
      return
    endif
    IF (.not. field_3d%in_order) THEN
      CALL set_error('this 3dfield not in order', 'get_column_from_field_3d')
      return
    endif

    do iLev = 1, field_3d%number_of_fields
      grid_data_from_3d => fu_grid_data(field_3d%fields(field_3d%vertical_order(iLev))%fp)
      arValues(iLev) = grid_data_from_3d(grid_index)
    end do

  end subroutine get_column_from_field_3d


  ! ***************************************************************

  INTEGER FUNCTION fu_number_of_fields(field_3d)

    ! Description:
    ! Returns the current number of fields stored in field_3d
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_3d_field), POINTER :: field_3d

    IF (defined(field_3d)) THEN
      fu_number_of_fields = field_3d%number_of_fields
    ELSE
      fu_number_of_fields = 0
    END IF

  END FUNCTION fu_number_of_fields


  ! ***************************************************************

  SUBROUTINE vertical_levels(field_3d, levels, number_of_levels)

    ! Description:
    ! Returns the vertical level system of 3D field. That is: the 2D
    ! levels  that data is available inside 3D field. The levels are
    ! returned in order, so that levels(1) is the lowest (closest to
    ! the ground) and levels(number_of_levels) the highest.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent OUT:
    TYPE(silja_level), DIMENSION(:), INTENT(out) :: levels

    ! Optional parameters with intent OUT:
    INTEGER, INTENT(out), OPTIONAL :: number_of_levels

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_3d_field), POINTER :: field_3d

    ! Local:
    INTEGER :: i

    DO i = 1, field_3d%number_of_fields
      levels(i) = fu_level(&
          & field_3d%fields(field_3d%vertical_order(i))%fp)
    END DO

    IF (PRESENT(number_of_levels))&
       & number_of_levels = field_3d%number_of_fields

  END SUBROUTINE vertical_levels


  ! ***************************************************************
  
  SUBROUTINE pressure_on_level(field_3d, i, p)

    ! Description:
    ! Returns the pressure [Pa] of i:th level of 3D field
    ! in every gridpoint.
    ! The horizontal grid must be checked elsewhere.
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
 
    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: i
    TYPE(silja_3d_field), POINTER :: field_3d

    ! Imported parameters with intent(out):
    REAL, DIMENSION(:), INTENT(out) :: p

    ! Local declarations:
    REAL :: a, b
    INTEGER :: fs
    TYPE(silja_field), POINTER :: field, p_surf_field
    REAL, DIMENSION(:), POINTER :: p_surf

    IF (.NOT.defined(field_3d)) THEN
      CALL set_error('undefined 3D-field given','pressure_on_level')
      RETURN
    END IF

    field => fu_field_from_3d_field(field_3d, i)
    fs = fu_size(field)
    IF (error) RETURN
    p(1:fs) = real_missing

    SELECT CASE (field_3d%leveltype)

      CASE (constant_pressure)
         p(1:fs) = fu_pr_level_pressure(fu_level(field))

      CASE (hybrid, layer_btw_2_hybrid)
        p_surf_field => fu_surface_pressure_field(field_3d)
        IF (error.or.(.NOT.defined(p_surf_field))) THEN
            CALL set_error('no surface pressure in hybrid field 3D','fu_pressure_on_level')
            RETURN
        END IF
        p_surf => fu_grid_data(p_surf_field)
        a = fu_hybrid_level_coeff_a(fu_level(field))
        b = fu_hybrid_level_coeff_b(fu_level(field))
        IF (error) RETURN
        !$omp workshare
          p(1:fs) = a + (b * p_surf(1:fs))
        !$omp END workshare

      CASE (sigma_level)
        p_surf_field => fu_surface_pressure_field(field_3d)
        IF (error.or.(.NOT.defined(p_surf_field))) THEN
          CALL set_error('no surface pressure in sigma field 3D','fu_pressure_on_level')
          RETURN
        END IF
        p_surf => fu_grid_data(p_surf_field)
        a = fu_sigma_level_sigma(fu_level(field))
        IF (error) RETURN
        p(1:fs) = a * p_surf(1:fs)

    CASE default
      CALL report(field_3d)
      CALL set_error('cannot make pressure, wrong leveltype','fu_pressure_on_level')

    END SELECT

  END SUBROUTINE pressure_on_level



  ! ****************************************************************
  ! ****************************************************************
  ! 
  !
  !       Private functions and subroutines of this module.
  !
  !
  ! ***************************************************************
  ! ****************************************************************


  SUBROUTINE set_field_3d_defined(field_3d)

    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_3d_field), pointer :: field_3d

    field_3d%defined = silja_true

  END SUBROUTINE set_field_3d_defined

  ! ***************************************************************

  LOGICAL FUNCTION fu_field_3d_defined(field_3d)

    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_3d_field), INTENT(in) :: field_3d

    fu_field_3d_defined = fu_true(field_3d%defined)

  END FUNCTION fu_field_3d_defined



  ! ****************************************************************

  FUNCTION fu_grid_of_3d_fi(field_3d)

    ! Returns the grid of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_grid) :: fu_grid_of_3d_fi
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_3d_field), POINTER :: field_3d

    fu_grid_of_3d_fi = field_3d%grid

  END FUNCTION fu_grid_of_3d_fi



  ! ***************************************************************

  FUNCTION fu_met_src_of_3d_fi(field_3d) result(met_src)

    ! Returns the source of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    type(meteo_data_source) :: met_src
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_3d_field), POINTER :: field_3d

    met_src = field_3d%met_src

  END FUNCTION fu_met_src_of_3d_fi



  ! ***************************************************************

  FUNCTION fu_valid_time_of_3d_fi(field_3d)

    ! Description:
    ! Returns the valid time of a 3d-field
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_time) :: fu_valid_time_of_3d_fi
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_3d_field), POINTER :: field_3d

    fu_valid_time_of_3d_fi = field_3d%valid_time

  END FUNCTION fu_valid_time_of_3d_fi


  ! ***************************************************************

  subroutine set_valid_time_of_3d_fi(field_3d, time)

    ! Description:
    ! Returns the valid time of a 3d-field
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_time) :: time
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_3d_field), POINTER :: field_3d

    field_3d%valid_time = time

  END subroutine set_valid_time_of_3d_fi



  ! ***************************************************************

  FUNCTION fu_analysis_time_of_3d(field_3d)

    ! Description:
    ! Returns the analysis time of a 3d-field
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_time) :: fu_analysis_time_of_3d
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_3d_field), POINTER :: field_3d

    TYPE(silja_field), POINTER :: field

    field => field_3d%fields(1)%fp

    fu_analysis_time_of_3d = fu_analysis_time(field)

  END FUNCTION fu_analysis_time_of_3d



  ! ***************************************************************

  FUNCTION fu_forecast_length_of_3d(field_3d)

    ! Description:
    ! Returns the analysis time of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_interval) :: fu_forecast_length_of_3d
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_3d_field), POINTER :: field_3d

    TYPE(silja_field), POINTER :: field

    field => field_3d%fields(1)%fp

    fu_forecast_length_of_3d = fu_forecast_length(field)

  END FUNCTION fu_forecast_length_of_3d



  ! ***************************************************************

  INTEGER FUNCTION fu_quantity_of_3d_fi(field_3d)

    ! Description:
    ! Returns the quantity of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_3d_field), POINTER :: field_3d

    fu_quantity_of_3d_fi = field_3d%quantity

  END FUNCTION fu_quantity_of_3d_fi


  !*****************************************************************

  function fu_SubstNm_of_3d_fld(field_3d) result(chNm)
    !
    ! Encapsulation of substance name of the 3d field
    !
    IMPLICIT NONE

    ! Return value
    character(len=clen) :: chNm

    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_3d_field), POINTER :: field_3d

    if(defined(field_3d%species)) then 
      chNm = fu_substance_name(field_3d%species)
    else
      chNm = ''
    endif

  end function fu_SubstNm_of_3d_fld


  ! ***************************************************************
  
  FUNCTION fu_mode_of_3d_fld(field_3d) result(mode)

    IMPLICIT NONE
    TYPE(silja_3d_field), POINTER :: field_3d
    type(Taerosol_mode) :: mode

    if(defined(field_3d))then
      mode = fu_mode(field_3d%species)
    else
      call set_error('undefined field','fu_mode_of_3d_fld')
      mode = aerosol_mode_missing
    end if

  end FUNCTION fu_mode_of_3d_fld

  ! ***************************************************************

  FUNCTION fu_optic_wavelen_of_3d_fld(field_3d)

    IMPLICIT NONE
    TYPE(silja_3d_field), POINTER :: field_3d
    real :: fu_optic_wavelen_of_3d_fld

    if(defined(field_3d))then
      fu_optic_wavelen_of_3d_fld =  fu_optical_wave_length(field_3d%species)
    else
      call set_error('undefined field','fu_optic_wavelen_of_3d_fld')
      fu_optic_wavelen_of_3d_fld = real_missing
    end if

  end FUNCTION fu_optic_wavelen_of_3d_fld

  ! ***************************************************************

  function fu_species_of_3d_fld(fld3d) result(species)
    implicit none
    type(silja_3d_field), intent(in) :: fld3d
    
    type(silam_species) :: species
    
    species = fld3d%species
    
  end function fu_species_of_3d_fld


  INTEGER FUNCTION fu_leveltype_of_3d_fi(field_3d)

    ! Description:
    ! Returns the leveltype of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_3d_field), POINTER :: field_3d

    fu_leveltype_of_3d_fi = field_3d%leveltype

  END FUNCTION fu_leveltype_of_3d_fi



  ! ***************************************************************

  INTEGER FUNCTION fu_3d_fieldsize(field_3d)

    ! Description:
    ! Returns the number of real-values stores in 3d-field.
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
    TYPE(silja_3d_field), POINTER :: field_3d
    !
    ! Local declarations:
    INTEGER :: i

    fu_3d_fieldsize = 0

    IF (.NOT.defined(field_3d)) RETURN

    DO i = 1, field_3d%number_of_fields
      fu_3d_fieldsize = fu_3d_fieldsize +&
          & fu_size(field_3d%fields(i)%fp)
    END DO


  END FUNCTION fu_3d_fieldsize



  ! ***************************************************************

  FUNCTION fu_find_level_scalar(field_3d, level) result(field)

    ! Description:
    ! Finds the field inside 3d-field whose level is the given level.
    ! If the level is not found, an error is set. The levels are gone
    ! through in vertical order from ground up.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_field), POINTER :: field

    ! Imported parameters with intent(in):
    TYPE(silja_level), INTENT(in) :: level
    TYPE(silja_3d_field), POINTER :: field_3d

    ! Local declarations:
    INTEGER :: i

    IF (.NOT.defined(field_3d)) THEN
      CALL set_error('undefined field_3d given','fu_find_level_scalar')
      RETURN
    END IF

    IF (fu_leveltype(level) /= field_3d%leveltype) THEN
      CALL report(field_3d)
      CALL report(level)
      CALL set_error('this type of level not found in this 3d-field'&
          & ,'fu_find_level_scalar')
      RETURN
    END IF


    DO i = 1, field_3d%number_of_fields
      field => field_3d%fields(field_3d%vertical_order(i))%fp
      IF (fu_cmp_levs_eq(fu_level(field), level)) RETURN
    END DO

    CALL set_error('correct level not found in 3d-field','fu_find_level_scalar')

  END FUNCTION fu_find_level_scalar



  ! ***************************************************************

  FUNCTION fu_lowest_field_scalar(field_3d)

    ! Description:
    ! Finds the field inside 3d-field that is closest to ground.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_field), POINTER :: fu_lowest_field_scalar

    ! Imported parameters with intent(in):
    TYPE(silja_3d_field), POINTER :: field_3d

    IF (.NOT.defined(field_3d)) THEN
      CALL set_error('undefined field_3d given'&
          & ,'fu_lowest_field_scalar')
      RETURN
    END IF

    fu_lowest_field_scalar => &
        & field_3d%fields(field_3d%vertical_order(1))%fp

  END FUNCTION fu_lowest_field_scalar



  ! ***************************************************************

  FUNCTION fu_highest_field_scalar(field_3d)

    ! Description:
    ! Finds the field inside 3d-field that is highest from ground.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_field), POINTER :: fu_highest_field_scalar

    ! Imported parameters with intent(in):
    TYPE(silja_3d_field), POINTER :: field_3d

    IF (.NOT.defined(field_3d)) THEN
      CALL set_error('undefined field_3d given'&
          & ,'fu_highest_field_scalar')
      RETURN
    END IF

    fu_highest_field_scalar => &
        & field_3d%fields(&
        & field_3d%vertical_order(field_3d%number_of_fields))%fp

  END FUNCTION fu_highest_field_scalar


  ! ***************************************************************

  FUNCTION fu_closest_field_scalar(field_3d, gridpoint, pressure, direction)

    ! Description:
    ! Finds the field inside 3d-field that is closest to the given
    ! pressure in given gridpoint-index, and in the given direction.
    ! Works for hybrid and pressure level 3d fields.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_field), POINTER :: fu_closest_field_scalar

    ! Imported parameters with intent(in):
    TYPE(silja_3d_field), POINTER :: field_3d
    TYPE(silja_field), POINTER :: surf_pre_field
    INTEGER, INTENT(in) :: gridpoint, direction
    REAL, INTENT(in) :: pressure

    ! Local declarations:
    INTEGER :: i, closest
    REAL, DIMENSION(:), POINTER :: surface_pressure
    REAL :: p_now, p_ground
    TYPE(silja_level) :: level
    LOGICAL :: on_level

    on_level = .false.

    IF (.NOT.defined(field_3d)) THEN
      CALL set_error('undef. field_3d given','fu_closest_field_scalar')
      RETURN
    END IF

    IF (fu_leveltype(field_3d) == hybrid) THEN
      surface_pressure => fu_grid_data(fu_surface_pressure_field(field_3d))
      p_ground = surface_pressure(gridpoint)
      IF (error) RETURN
    END IF

    DO i = 1, field_3d%number_of_fields
      level = fu_level(fu_field_from_3d_field(field_3d, i))
      IF (error) RETURN

      SELECT CASE (fu_leveltype(level))
        CASE (constant_pressure)
        p_now = fu_pr_level_pressure(level)
        CASE (hybrid)
        p_now = fu_hybrid_level_pressure(level, p_ground)
      END SELECT

      on_level = (p_now.eps.pressure) 
      IF (p_now < pressure) EXIT
    END DO

    IF (on_level) THEN
      closest = i
    ELSE
      SELECT CASE (direction)
        CASE (upwards)
        closest = i

        CASE (downwards) 
        closest = i - 1

      CASE default
        CALL set_error('direction must be upwards or downwards',&
            & 'fu_closest_field_scalar')
      END SELECT
    END IF

    !    PRINT *, 'Closest:', closest
    closest = MAX(closest, 1)
    closest = MIN(closest, field_3d%number_of_fields)
    fu_closest_field_scalar => fu_field_from_3d_field(field_3d, closest)
        
  END FUNCTION fu_closest_field_scalar





  ! ***************************************************************

  FUNCTION fu_pressures_scalar(field_3d, position) result(pressures)

    ! Description:
    ! This function returns the pressures of it fields at position.
    ! Since no weather data is available to this module, calculation
    ! works on ly for constant pressure fields and hybrid level 3d
    ! -fields containing surface pressure. Values returned [Pa]
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
    REAL, DIMENSION(max_levels) :: pressures

    ! Imported parameters with intent(in):
    TYPE(silja_3d_field), POINTER :: field_3d
    TYPE(silam_grid_position), INTENT(in) :: position

    ! Local declarations:
    INTEGER :: i
    REAL :: surface_pressure

    IF (.NOT.defined(field_3d)) THEN
      CALL set_error('undefined field_3d given'&
          & ,'fu_pressures_scalar')
      RETURN
    END IF

    SELECT CASE (field_3d%leveltype)

      ! -----------------------------------------------------
      ! 
      ! 1. Pressures of constant pressure level system.
      !    -------------------------------------------

      CASE (constant_pressure)

      DO i = 1, field_3d%number_of_fields
        pressures(i) = fu_pr_level_pressure(&
           & fu_level(&
           & field_3d%fields(field_3d%vertical_order(i))%fp))
      END DO


      ! -----------------------------------------------------
      ! 
      ! 2. Pressures of hybrid level system.
      !    --------------------------------

      CASE (hybrid)

      IF (.NOT.defined(field_3d%surface_pressure)) THEN
        CALL set_error('field_3d not containig surface pressure'&
            & ,'fu_pressures_scalar')
        RETURN
      END IF

      surface_pressure = fu_interpolate_to_position(&
          & field_3d%surface_pressure, position, linear)
      IF (error) RETURN

      DO i = 1, field_3d%number_of_fields
        pressures(i) = fu_hybrid_level_pressure(&
            & fu_level(&
            & field_3d%fields(field_3d%vertical_order(i))%fp),&
            & surface_pressure)
      END DO


      ! -----------------------------------------------------
      ! 
      ! 3. Unknown level type.
      !    ------------------

    CASE default
      PRINT *, field_3d%leveltype
      CALL set_error('cannot handle this type of levels here'&
          & ,'fu_pressures_scalar')

    END SELECT

  END FUNCTION fu_pressures_scalar



  ! ***************************************************************

  REAL FUNCTION fu_level_pressure_in_3d(field_3d, i, position)

    ! Description:
    ! Returns the pressure of i:th level of 3d field at position. No
    ! weather data can be used here except what is already in 3d
    ! field.
    !
    ! All units: SI [Pa]
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silam_grid_position), INTENT(in) :: position
    INTEGER, INTENT(in) :: i ! the order from bottom: 1 is lowest

    TYPE(silja_3d_field), POINTER :: field_3d

    ! Local declarations:
    TYPE(silja_level) :: level
    REAL :: ground_pressure

    ! The level of i:th field from bottom:
    level = fu_level(field_3d%fields(field_3d%vertical_order(i))%fp)

    SELECT CASE (fu_leveltype(level))

      CASE (constant_pressure)

      fu_level_pressure_in_3d = fu_pr_level_pressure(level)

      CASE (hybrid)

      ground_pressure = fu_interpolate_to_position(&
          & field_3d%surface_pressure,&
          & position,&
          & nearest_point)
      IF (error) RETURN

      fu_level_pressure_in_3d = fu_hybrid_level_pressure(&
          & level,&
          & ground_pressure)
      !      PRINT *, 'scalar 3d level pressure', fu_level_pressure_in_3d

    CASE default

      CALL set_error('unknwon leveltype','fu_level_pressure_in_3d')

    END SELECT

  END FUNCTION fu_level_pressure_in_3d


  ! ***************************************************************

  FUNCTION fu_surf_pre_field_scalar(field_3d) 

    ! Description:
    ! Returns a pointer to th surface pressure field in a field_3d,
    ! in case it exists. If not pointer is nullified (undefined field)
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_field), POINTER :: fu_surf_pre_field_scalar

    ! Imported parameters with intent(in):
    TYPE(silja_3d_field), POINTER :: field_3d

    IF (ASSOCIATED(field_3d%surface_pressure)) THEN
      fu_surf_pre_field_scalar => field_3d%surface_pressure
    ELSE
      NULLIFY(fu_surf_pre_field_scalar)
    END IF


  END FUNCTION fu_surf_pre_field_scalar



  ! ***************************************************************

  SUBROUTINE find_closest_fields(&
      & position,&
      & field_3d,&
      & field_upwards,&
      & field_upwards_pressure,&
      & field_downwards,&
      & field_downwards_pressure,&
      & above_the_highest,&
      & below_the_lowest)

    ! Description:
    ! Inside a 3d-field finds the closest ones in vertical in both
    ! directions. If position if above the highest level, then
    ! field_upwards if nullified. Similarly for bottom.
    ! 
    ! Method:
    ! 
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent(in):
    TYPE(silam_grid_position), INTENT(in) :: position

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_3d_field), POINTER :: field_3d !in
    TYPE(silja_field), POINTER :: field_upwards,& !out
        & field_downwards

    ! Imported parameters with intent OUT:
    REAL, INTENT(out) :: field_upwards_pressure,&
        & field_downwards_pressure
    LOGICAL, INTENT(out) :: above_the_highest, below_the_lowest

    ! Local declarations: 
    REAL, DIMENSION(max_levels) :: field_pressures
    REAL, DIMENSION(max_levels) :: pr_differences
    REAL :: position_pressure
    TYPE(silja_field), POINTER :: field
    INTEGER :: nn
    INTEGER, DIMENSION(1) :: up_index, down_index


    !----------------------------------------
    !
    ! 1. Initial values.
    !    --------------

    IF (defined(position)) THEN
      position_pressure = fu_pressure(position)
    ELSE
      CALL set_error('undefined position given','find_closest_fields')
      RETURN
    END IF


    !----------------------------------------
    !
    ! 2. Get the pressure of fields at position.
    !    --------------------------------------

    field_pressures = fu_pressures_scalar(field_3d, position)
    IF (error) RETURN
    nn = field_3d%number_of_fields
    below_the_lowest = .false.
    above_the_highest = .false.


    !----------------------------------------------------------
    !
    ! 3. Check if we are below the lowest level
    !    --------------------------------------

    IF (position_pressure > field_pressures(1)) THEN   
      field_upwards => fu_lowest_field(field_3d)
      NULLIFY(field_downwards)
      field_upwards_pressure = field_pressures(1)
      field_downwards_pressure = real_missing
      below_the_lowest = .true.
      RETURN
    END IF


    !----------------------------------------------------------
    !
    ! 4. Check if we are above the highest level
    !    ---------------------------------------

    IF (position_pressure < field_pressures(nn)) THEN   
      NULLIFY(field_upwards)
      field_downwards => fu_highest_field(field_3d)
      field_upwards_pressure = real_missing
      field_downwards_pressure = field_pressures(nn)
      above_the_highest = .true.
      RETURN
    END IF


    !----------------------------------------------------------
    !
    ! 5. Between highest and lowest, find closest
    !    ---------------------------------------

    pr_differences(1:nn) = field_pressures(1:nn) - position_pressure

    ! Greatest negative difference:
    up_index = MAXLOC(pr_differences(1:nn),&
        & mask=pr_differences(1:nn) <= 0.)

    ! Smallest positive difference:
    down_index = MINLOC(pr_differences(1:nn),&
        & mask=pr_differences(1:nn) >= 0.)

    field_downwards => fu_field_from_3d_field(field_3d, down_index(1))
    field_downwards_pressure = field_pressures(down_index(1))

    field_upwards => fu_field_from_3d_field(field_3d, up_index(1))
    field_upwards_pressure = field_pressures(up_index(1))


  END SUBROUTINE find_closest_fields



  ! ***************************************************************

  SUBROUTINE find_closest_fields_2(gridpoint,&
                                 & pressure, &
                                 & field_3d,&
                                 & field_upwards,&
                                 & field_upwards_pressure,&
                                 & field_downwards,&
                                 & field_downwards_pressure,&
                                 & above_the_highest,&
                                 & below_the_lowest)

    ! Description:
    ! Inside a 3d-field finds the closest ones in vertical in both
    ! directions. If position if above the highest level, then
    ! field_upwards if nullified. Similarly for bottom.
    ! 
    ! Method:
    ! 
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: gridpoint
    REAL, INTENT(in) :: pressure

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_3d_field), POINTER :: field_3d !in
    TYPE(silja_field), POINTER :: field_upwards,& !out
        & field_downwards

    ! Imported parameters with intent OUT:
    REAL, INTENT(out) :: field_upwards_pressure,&
        & field_downwards_pressure
    LOGICAL, INTENT(out) :: above_the_highest, below_the_lowest

    ! Local declarations: 
    REAL :: p_ground, p_now, p_previous
    REAL, DIMENSION(:), POINTER :: surface_pressure
    INTEGER :: up_index, down_index
    INTEGER :: i
    TYPE(silja_level) :: level
    LOGICAL :: on_level, hybrid_level

    above_the_highest = .false.
    below_the_lowest = .false.
    on_level = .false.
    hybrid_level = .false.

    !----------------------------------------
    !
    ! 1. Find surface pressure for hybrid level
    !    ------------------

    SELECT CASE (fu_leveltype(field_3d))

      CASE (hybrid)
      hybrid_level= .true.
      surface_pressure => fu_grid_data(fu_surface_pressure_field(field_3d))
      p_ground = surface_pressure(gridpoint)

      CASE (constant_pressure)
      hybrid_level = .false.

    CASE default
      CALL report(field_3d)
      CALL set_error('sorry cannot handle','find_closest_fields_2')
      RETURN
    END SELECT


    !----------------------------------------
    !
    ! 2. Loop over  levels
    !    ------------------

    level_loop: DO i = 1, field_3d%number_of_fields
      level = fu_level(fu_field_from_3d_field(field_3d, i))
      IF (hybrid_level) THEN
        p_now = fu_hybrid_level_pressure(level, p_ground)
      ELSE
        p_now = fu_pr_level_pressure(level)
      END IF
      IF (error) RETURN

      ! On level:
      IF (pressure.eps.p_now) THEN 
        on_level= .true.
        up_index = i
        EXIT level_loop
      END IF

      ! Between levels or below the lowest (level with lower p found):
      IF (pressure > p_now) THEN
        below_the_lowest = (i == 1) 
        up_index = i
        down_index = i - 1
        EXIT level_loop
      ELSE
        p_previous = p_now
      END IF

      ! Above highest:
      IF (i == field_3d%number_of_fields) THEN
        IF (pressure > p_now)  THEN
          down_index = i
          above_the_highest = .true.
          EXIT
        END IF
      END IF
      
    END DO level_loop


    !----------------------------------------
    !
    ! 3. Set results
    !    -----------
    IF (on_level) THEN
      field_upwards => fu_field_from_3d_field(field_3d, up_index)
      field_downwards => field_upwards
      field_upwards_pressure = p_now
      field_downwards_pressure = p_now
      RETURN
    END IF

    IF (below_the_lowest) THEN
      NULLIFY(field_downwards)
      field_downwards_pressure = real_missing
      field_upwards => fu_lowest_field(field_3d)
      field_upwards_pressure = p_now
      RETURN
    END IF

    IF (above_the_highest) THEN
      NULLIFY(field_upwards)
      field_downwards => fu_highest_field(field_3d)
      field_upwards_pressure = real_missing
      field_downwards_pressure = p_now
      RETURN
    END IF

    field_upwards => fu_field_from_3d_field(field_3d, up_index)
    field_downwards => fu_field_from_3d_field(field_3d, down_index)
    field_upwards_pressure = p_now
    field_downwards_pressure = p_previous

  END SUBROUTINE find_closest_fields_2


  ! ***************************************************************



  !***********************************************************************

  subroutine ddz_of_field_3d(fld3d, height3d, grid, iLev, d_dz)

    !
    ! Derivative along the metrical vertical (or whatever is supplied as height 3d field)
    !
    implicit none

    ! Imported variables
    TYPE(silja_3d_field), intent(in) :: fld3d, height3d
    type(silja_grid), intent(in) :: grid
    integer, intent(in) :: iLev
    real, dimension(:), pointer :: d_dz

    ! Local variables
    integer :: iTmp, nx, ny
    real, dimension(:), pointer :: fld, fld_up, fld_down, h, h_up, h_down
    real :: fA, fB

    !
    ! Checking and preparation
    !
    if(.not. associated(d_dz))then
      call set_error('Output pointer is not associated','ddz_of_field_3d')
      return
    endif

    call grid_dimensions(grid, nx, ny)

    if(error .or. nx*ny > size(d_dz))then
      call set_error('Output pointer is not associated','ddz_of_field_3d')
      return
    endif

    if(.not. defined(grid))then
      call set_error('Grid is undefined','ddz_of_field_3d')
      return
    endif

    !
    ! Vertical derivative is taken via square-approximation for all levels except for
    ! the first and the last where the upward and downward first-order derivative is used
    !
    if(iLev == 1)then
      !
      ! Upward first-order
      !
      fld => fu_grid_data(fld3d%fields(iLev)%fp)
      fld_up => fu_grid_data(fld3d%fields(iLev+1)%fp)

      h => fu_grid_data(height3d%fields(iLev)%fp)
      h_up => fu_grid_data(height3d%fields(iLev+1)%fp)

      do iTmp = 1, nx*ny
        d_dz(iTmp) = (fld_up(iTmp) - fld(iTmp)) / (h_up(iTmp) - h(iTmp))
      end do

    elseif(iLev == height3d%number_of_fields)then
      !
      ! Downward first-order
      !
      fld => fu_grid_data(fld3d%fields(iLev)%fp)
      fld_down => fu_grid_data(fld3d%fields(iLev-1)%fp)

      h => fu_grid_data(height3d%fields(iLev)%fp)
      h_down => fu_grid_data(height3d%fields(iLev-1)%fp)

      do iTmp = 1, nx*ny
        d_dz(iTmp) = (fld(iTmp) - fld_down(iTmp)) / (h(iTmp) - h_down(iTmp))
      end do

    else
      !
      ! Second-order up-and downward approximation. See notebook 9, p.4
      !
      fld => fu_grid_data(fld3d%fields(iLev)%fp)
      fld_up => fu_grid_data(fld3d%fields(iLev+1)%fp)
      fld_down => fu_grid_data(fld3d%fields(iLev-1)%fp)

      h => fu_grid_data(height3d%fields(iLev)%fp)
      h_up => fu_grid_data(height3d%fields(iLev+1)%fp)
      h_down => fu_grid_data(height3d%fields(iLev-1)%fp)

      do iTmp = 1, nx*ny
        fA = ((fld(iTmp)-fld_up(iTmp)) / (h(iTmp)-h_up(iTmp)) - &
            & (fld_down(iTmp)-fld(iTmp)) / (h_down(iTmp)-h(iTmp))) / &
           & (h_up(iTmp)-h_down(iTmp))
        fB = (fld_down(iTmp)-fld(iTmp)) / (h_down(iTmp)-h(iTmp)) - fA * (h_down(iTmp)+h(iTmp))
        d_dz(iTmp) = 2. * fA * h(iTmp) + fB
      end do

    endif  ! value of iLev

  end subroutine ddz_of_field_3d


  !**********************************************************************

  subroutine check_field_3d(field_3d, fLowerLimit, fUpperLimit)
    !
    ! Checks validity of all the fields in the field_3d scanning the whole field
    !
    implicit none

    ! Improted parameters
    TYPE(silja_3d_field), intent(in) :: field_3d
    real, intent(in) :: fLowerLimit, fUpperLimit

    ! Local variables
    integer :: iLev, ix, iy, fs, nx, ny, iTmp
    logical :: ifError
    real, dimension(:), pointer :: f_data

    if(.not. defined(field_3d))then
      call set_error('Given field is undefined','check_field_3d')
      return
    endif

    call grid_dimensions(field_3d%grid, nx, ny)
    fs = nx * ny
    if(fs < 1 .or. fs > 1000000)then
      call msg('Strange grid size:',fs)
      call set_error('Strange grid size','check_field_3d')
      return
    endif
    
    ifError = .false.

    IF (.not.field_3d%in_order) THEN
      call msg_warning('Field_3d is not in vertical order')
    endif

    DO iLev = 1, field_3d%number_of_fields
      if(field_3d%in_order) THEN
        f_data => fu_grid_data(field_3d%fields(field_3d%vertical_order(iLev))%fp)
      else
        f_data => fu_grid_data(field_3d%fields(iLev)%fp)
      endif
      if(associated(f_data))then
        if(size(f_data) >= fs .and. size(f_data) < 10000000)then
          do iy = 1, ny
            do ix = 1, nx
              iTmp = ix + (iy -1)*nx
              if(f_data(iTmp) > fUpperLimit .or. f_data(iTmp) < fLowerLimit)then
                call msg('Data outside the limit: ix,iy',ix,iy)
                call msg('Level and value:',iLev,f_data(iTmp))
                ifError = .true.
              endif
            enddo
          enddo
        else
          call msg('Strange data array for level:',iLev)
          call msg('Size of data array and field_3d grid:',size(f_data),fs)
          ifError = .true.
          cycle
        endif
      else
        call msg('Not associated data array for level:',iLev)
        ifError = .true.
        cycle
      endif
    END DO   ! iLev

    if(ifError)then
      call msg_warning('Problems have been found in the following field','check_field_3d')
      call report(field_3d)
      call set_error('Problems have been found in the above field','check_field_3d')
    endif

  end subroutine check_field_3d


  ! ***************************************************************
  ! ***************************************************************
  ! 
  !
  !                 TESTS
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  SUBROUTINE print_field_3d_report(field_3d)

    ! Description:
    ! Print on standard output the contents of vertical field field_3d.
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silja_3d_field) :: field_3d

    ! Local declarations:
    INTEGER :: i
    TYPE(silja_field), POINTER :: field
    character(len=worksize_string) :: str, str1

    call msg(' *********** Field_3d report *************')


    IF (.NOT.defined(field_3d)) THEN
      call msg(' Undefined field_3d.')

    ELSE

      call msg(fu_connect_strings(' met_src and quantity: ', &
                                & fu_name(field_3d%met_src),&
                                & fu_quantity_string(field_3d%quantity)))


      if(field_3d%number_of_fields < 1)then
        call msg('No fields are stored in this 3d field')
      else
        call msg('Number of 2d fields: ', field_3d%number_of_fields)
        if(field_3d%vertical_order(1) > 0 .and. field_3d%vertical_order(1) < field_3d%number_of_fields)then
          call msg(' Full report of lowest field: ')
          CALL report(field_3d%fields(field_3d%vertical_order(1))%fp)
        else
          call msg(' Full report of one of fields:')
          CALL report(field_3d%fields(1)%fp)
        endif
      endif

      IF (field_3d%in_order) THEN
        call msg('Data available on following levels from down up:')
      ELSE
        call msg('not in order vertically:')
      END IF

      call msg("Data summary:")
        DO i = 1, field_3d%number_of_fields
          field => field_3d%fields(field_3d%vertical_order(i))%fp
        if(defined(field)) then 
           call level_to_short_string(fu_level(field), str1)
           WRITE (str, *) 'Level: ',i, ' ', trim(str1),' Field min, avg and max: ', &
                   & MINVAL(fu_grid_data(field)),&
                   & (SUM(fu_grid_data(field))/REAL(SIZE(fu_grid_data(field)))),&
                   & MAXVAL(fu_grid_data(field))
           call msg(str)
        endif
        END DO

    END IF

    IF (defined(field_3d%surface_pressure)) THEN
      CALL report(field_3d%surface_pressure)
    ELSE
      call msg(' No surface pressure field added to this 3d-field.')
    END IF

    call msg(' *********** End-of-field_3d report *************')

  END SUBROUTINE print_field_3d_report


END MODULE fields_3d

