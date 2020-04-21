MODULE windfields_3d 

  ! Description:
  ! This module contains definition and tools for 
  ! 3-dimensional scalar fields (x,y,z) of meteorological data.
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
!  USE vectors
  USE windfields
!  USE fields
!  USE levels

  IMPLICIT NONE

  ! The public functions and subroutines available in this module:

  ! Change contents of 3d-field:
  PUBLIC add_windfield_to_3d_windfield
  PUBLIC set_3d_windfield_empty
  PUBLIC organize_windfields_vertically
  PUBLIC add_surf_pre_to_3d_wind

  ! Return stuff from 3d-field:
  PUBLIC fu_windfield_from_3d_wind  ! i:th field from bottom
  PUBLIC fu_find_level ! find 2d-windfield by its level
  PUBLIC fu_lowest_field ! field closest to ground
  PUBLIC fu_highest_field ! field highest from ground
  PUBLIC fu_interpolate_to_position ! interpolate in x, y and z
  PUBLIC horizontal_interp_multiwind ! used to make soundings
  PUBLIC fu_surface_pressure_field ! field inside 3d-field
  PUBLIC fu_pressures

  ! Administrative data found in 3d-field:
  PUBLIC fu_number_of_windfields
  PUBLIC fu_size
  PUBLIC fu_valid_time
  PUBLIC fu_met_src
  PUBLIC fu_u_grid_from_3d
  PUBLIC fu_v_grid_from_3d
  PUBLIC fu_w_grid_from_3d
  PUBLIC fu_leveltype
  PUBLIC report ! for testing
  public defined

  ! The private functions and subroutines not to be used elsewhere:
  PRIVATE fu_windfield_3d_defined
  PRIVATE fu_met_src_of_3d_wind
  PRIVATE fu_valid_time_of_3d_wind
  PRIVATE fu_u_grid_of_3d_wind
  PRIVATE fu_v_grid_of_3d_wind
  PRIVATE fu_w_grid_of_3d_wind
  PRIVATE fu_3d_windfieldsize
  PRIVATE fu_find_level_wind
  PRIVATE fu_leveltype_of_3d_wind
  PRIVATE fu_pressures_wind
  PRIVATE fu_interpolate_wind_3d_to_pos
  PRIVATE find_closest_windfields ! vertically from a position
  PRIVATE print_windfield_3d_report ! for testing
  private fu_surf_pre_field_wind
  private fu_lowest_field_wind
  private fu_highest_field_wind
  private fu_level_pressure_in_3d

  ! Generic names and operator-interfaces of some functions:
  INTERFACE defined
    MODULE PROCEDURE fu_windfield_3d_defined
  END INTERFACE

  INTERFACE report
    MODULE PROCEDURE print_windfield_3d_report
  END INTERFACE

  INTERFACE fu_size
    MODULE PROCEDURE fu_3d_windfieldsize
  END INTERFACE

  INTERFACE fu_find_level
    MODULE PROCEDURE fu_find_level_wind
  END INTERFACE

  INTERFACE fu_met_src  
    MODULE PROCEDURE fu_met_src_of_3d_wind
  END INTERFACE

  INTERFACE fu_interpolate_to_position
    MODULE PROCEDURE fu_interpolate_wind_3d_to_pos
  END INTERFACE

  INTERFACE fu_valid_time
    MODULE PROCEDURE fu_valid_time_of_3d_wind
  END INTERFACE

  INTERFACE fu_u_grid_from_3d 
    MODULE PROCEDURE fu_u_grid_of_3d_wind
  END INTERFACE

  INTERFACE fu_v_grid_from_3d
    MODULE PROCEDURE fu_v_grid_of_3d_wind
  END INTERFACE

  INTERFACE fu_w_grid_from_3d
    MODULE PROCEDURE fu_w_grid_of_3d_wind
  END INTERFACE

  INTERFACE fu_lowest_field
    MODULE PROCEDURE fu_lowest_field_wind
  END INTERFACE

  INTERFACE fu_highest_field
    MODULE PROCEDURE fu_highest_field_wind
  END INTERFACE

  INTERFACE fu_surface_pressure_field
    MODULE PROCEDURE fu_surf_pre_field_wind
  END INTERFACE

  INTERFACE fu_leveltype
    MODULE PROCEDURE fu_leveltype_of_3d_wind
  END INTERFACE

  INTERFACE fu_pressures
    MODULE PROCEDURE fu_pressures_wind
  END INTERFACE



  ! Public types with private components defined in this module:
  TYPE silja_3d_windfield
    PRIVATE
    TYPE(silja_time) :: valid_time
    type(meteo_data_source) :: met_src
    TYPE(silja_grid), DIMENSION(3) :: grids ! u,v,w
    INTEGER :: leveltype ! leveltype same for all fields in 3d-field!!
    TYPE(silja_windfieldpointer), DIMENSION(max_levels) :: fields
    INTEGER, DIMENSION(max_levels) :: vertical_order
    LOGICAL :: in_order
    TYPE(silja_field), POINTER :: surface_pressure
    INTEGER :: number_of_fields
    TYPE(silja_logical) :: defined
  END TYPE silja_3d_windfield

  ! Public types with public components defined in this module:
  TYPE silja_windfield_3d_pointer
    TYPE(silja_3d_windfield), POINTER :: fp
  END TYPE silja_windfield_3d_pointer

CONTAINS

  ! ***************************************************************

  SUBROUTINE add_windfield_to_3d_windfield(wind, wind_3d)

    ! Description:
    ! Adds a field to vertical field-field_3d. 
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_windfield), POINTER :: wind
    TYPE(silja_3d_windfield), INTENT(inout) :: wind_3d

    ! Local declarations:
    INTEGER :: i

    !----------------------------------------
    !
    ! 1. Check status of wind_3d and the field that should be added
    !    ----------------------------------------------------------

    IF (.NOT.defined(wind)) THEN
      CALL set_error('cannot put undefined wind to wind_3d','add_windfield_to_3d_windfield')
      RETURN
    END IF

    ! For now:
!    IF (fu_cmp_levs_eq(fu_level(wind), level_10m_above_ground)) RETURN

    new_or_old: IF (.NOT.defined(wind_3d)) THEN

      !----------------------------------------
      !
      ! 2. A new, untouched wind_3d
      !    -----------------------

      wind_3d%fields(1)%fp => wind
      wind_3d%number_of_fields = 1
      IF (error) RETURN

      wind_3d%valid_time = fu_valid_time(wind)
      IF (error) RETURN

      wind_3d%met_src = fu_met_src(wind)
      IF (error) RETURN

      wind_3d%leveltype = fu_leveltype(fu_level(wind))
      IF (error) RETURN

      wind_3d%grids(1) = fu_u_grid(wind)
      wind_3d%grids(2) = fu_v_grid(wind)

      IF (fu_w_available(wind)) THEN
	wind_3d%grids(3) = fu_w_grid(wind)
      ELSE
	wind_3d%grids(3) = grid_missing
      END IF

      IF (error) RETURN

      ! No surface pressure given yet:
      NULLIFY(wind_3d%surface_pressure)

      wind_3d%in_order = .false.

      ! Everything ok:
      wind_3d%defined = fu_set_true()
    ELSE

      !----------------------------------------
      !
      ! 3. Already stuff in field_3d, check that it isn't already
      ! there, and check that the new field belogs here.
      !    -----------------------

      IF (.NOT.fu_field_belogs_to_this_3d()) RETURN


      !----------------------------------------
      !
      ! 4. Already stuff in wind_3d, add to bottom
      !    -----------------------


      IF (wind_3d%number_of_fields < max_levels) THEN

	wind_3d%number_of_fields = wind_3d%number_of_fields + 1

	wind_3d%fields(wind_3d%number_of_fields)%fp => wind

	wind_3d%in_order = .false.

      ELSE
        CALL set_error('cannot add any more field to this wind_3d','add_windfield_to_3d_windfield')
      END IF

    END IF new_or_old

  CONTAINS

    LOGICAL FUNCTION fu_field_belogs_to_this_3d()

      ! Description:
      ! Returns true value, if the current field should be added to
      ! the current 3d-field.
      !
      ! Method:
      ! 1. Check that the field is not already in 3d-field
      ! 2. Check that quantity, time met_src and leveltype match.
      !
      IMPLICIT NONE

      ! Local declarations:
      INTEGER :: i

      fu_field_belogs_to_this_3d = .false.

      if(.not. defined(wind))then
        call set_error('Undefined wind field','fu_field_belogs_to_this_3d')
        return
      endif

      ! Already in 3d-field?
      DO i = 1, wind_3d%number_of_fields
        IF (fu_id(wind) == fu_id(wind_3d%fields(i)%fp)) RETURN
      END DO

      IF (fu_leveltype(fu_level(wind)) /= wind_3d%leveltype) THEN
!        CALL report(wind)
!        CALL msg_warning('cannot mix different leveltypes','fu_field_belogs_to_this_3d')
        RETURN
      END IF


      IF (.NOT.(fu_u_grid(wind) == wind_3d%grids(1))) THEN
        call msg('input wind:')
        CALL report(wind)
        call msg('3D wind w-grid:')
        call report(wind_3d%grids(3))
        CALL set_error('cannot mix different u-grids','fu_field_belogs_to_this_3d')
        RETURN
      END IF

      IF (.NOT.(fu_v_grid(wind) == wind_3d%grids(2))) THEN
        call msg('input wind:')
        CALL report(wind)
        call msg('3D wind w-grid:')
        call report(wind_3d%grids(3))
        CALL set_error('cannot mix different v-grids','fu_field_belogs_to_this_3d')
        RETURN
      END IF

      IF (fu_w_available(wind)) THEN
        IF (.NOT.(fu_w_grid(wind) == wind_3d%grids(3))) THEN
          call msg('input wind:')
          CALL report(wind)
          call msg('3D wind w-grid:')
          call report(wind_3d%grids(3))
          CALL set_error('cannot mix different w-grids','fu_field_belogs_to_this_3d')
          RETURN
        END IF
      END IF

      IF (fu_valid_time(wind) /= wind_3d%valid_time) RETURN
      IF (.not.fu_met_src(wind) == wind_3d%met_src) RETURN

      fu_field_belogs_to_this_3d = .true.

    END FUNCTION fu_field_belogs_to_this_3d

  END SUBROUTINE add_windfield_to_3d_windfield



  ! ***************************************************************

  SUBROUTINE organize_windfields_vertically(field_3d)

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
    TYPE(silja_3d_windfield), INTENT(inout) :: field_3d

    ! Local declarations:
    INTEGER :: i, j
    TYPE(silja_windfield), POINTER :: field
    INTEGER :: lower_count
    TYPE(silja_level) :: level
    TYPE(silja_3d_windfield), TARGET :: fp
    TYPE(silja_3d_windfield), POINTER :: fpp

    !----------------------------------------
    !
    ! 1. Check status of field_3d
    !    ------------------------

    IF (.NOT.defined(field_3d)) THEN
      CALL set_error('undefined field_3d given','organize_windfields_vertically')
      RETURN
    END IF

    IF (field_3d%in_order) RETURN


    !----------------------------------------
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

!!!$    fp = field_3d
!!!$    fpp => fp
!!!$
!!!$    CALL print_windfield_3d_report(fpp)
!!!$    STOP

  END SUBROUTINE organize_windfields_vertically



  ! ***************************************************************

  SUBROUTINE add_surf_pre_to_3d_wind(wind_3d, surface_pressure)

    ! Description:
    ! Adds a correnponding ground surface pressure field to a 3d
    ! -windfield and checks it. Grid is not checked, since in hirlam
    ! windcomponents and scalar fields are in slightly different
    ! grids. A course check is done, though.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_3d_windfield), POINTER :: wind_3d
    TYPE(silja_field), POINTER :: surface_pressure

    IF (.NOT.defined(wind_3d)) THEN
      CALL set_error('undefined wind_3d given'&
	  & ,'add_surf_pre_to_3d_wind')
      RETURN
    END IF

    IF (.NOT.defined(surface_pressure)) THEN
      CALL set_error('undefined surface_pressure given'&
	  & ,'add_surf_pre_to_3d_wind')
      RETURN
    END IF

    matching_check: IF (&
!	& (fu_met_src(surface_pressure) == wind_3d%met_src).and.&
	& (fu_valid_time(surface_pressure) == wind_3d%valid_time)&
	& .and.(fu_grids_match_closely(&
	& fu_grid(surface_pressure),wind_3d%grids(1))))THEN

      wind_3d%surface_pressure => surface_pressure

    ELSE

      CALL set_error('fields do not match','add_surf_pre_to_3d_wind')
      CALL report(wind_3d)
      CALL report(surface_pressure)

    END IF matching_check

  END SUBROUTINE add_surf_pre_to_3d_wind


  ! ***************************************************************

  SUBROUTINE set_3d_windfield_empty(field_3d)

    ! Description:
    !  
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_3d_windfield), INTENT(inout) :: field_3d

    INTEGER :: i

    field_3d%defined = silja_false
    field_3d%number_of_fields = 0

    DO i = 1, SIZE(field_3d%fields)
      NULLIFY(field_3d%fields(i)%fp)
    END DO

  END SUBROUTINE set_3d_windfield_empty



  ! ***************************************************************

  SUBROUTINE horizontal_interp_multiwind(&
      & field_3d,&
      & position, &
      & method, &
      & values,&
      & levels,&
      & pressures,&
      & number_of_values)

    ! Description:
    ! Makes horizontal interpolation fast for all levels. Used to
    ! make soundings.
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

    ! Imported parameters with intent in:
    INTEGER, INTENT(in) :: method
    TYPE(silam_grid_position), INTENT(in) :: position

    ! Imported parameters with intent OUT:
    TYPE(silja_velocity), DIMENSION(max_levels), INTENT(out)&
	& :: values
    REAL, DIMENSION(max_levels), INTENT(out) :: pressures
    TYPE(silja_level), DIMENSION(max_levels), INTENT(out) ::&
	& levels
    INTEGER, INTENT(out) :: number_of_values

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_3d_windfield), POINTER :: field_3d

    ! Local declarations:
    INTEGER :: i

    !----------------------------------------
    !
    ! 1. Checkings.
    !    ---------

    IF (.NOT.defined(field_3d)) THEN
      CALL set_error('undefined field_3d given'&
	  & ,'fu_interp_scalar_3d_to_sounding')
      RETURN
    END IF

    IF (.NOT.field_3d%in_order) THEN
      CALL organize_windfields_vertically(field_3d)
      IF (error) RETURN
    END IF


    !----------------------------------------
    !
    ! 2. Interpolations. THIS SHOULD BE CHANGED TO MORE EFFICIENT!
    !    --------------

    DO i = 1, field_3d%number_of_fields

      values(i) = fu_interpolate_to_position(&
	  & field_3d%fields(field_3d%vertical_order(i))%fp,&
	  & position,&
	  & method)
      IF (error) RETURN

      levels(i) = fu_level(&
	  & field_3d%fields(field_3d%vertical_order(i))%fp)

      pressures(i) = fu_level_pressure_in_3d(&
	  & field_3d,&
	  & i,&
	  & position)

    END DO

    number_of_values = field_3d%number_of_fields


  END SUBROUTINE horizontal_interp_multiwind





  ! ***************************************************************

  FUNCTION fu_windfield_from_3d_wind(field_3d, number, ifSkipOrder) 

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
    TYPE(silja_windfield), POINTER :: fu_windfield_from_3d_wind

    ! Imported parameters with intent(in):
!    TYPE(silja_3d_windfield), POINTER :: field_3d
    TYPE(silja_3d_windfield), INTENT(in):: field_3d
    INTEGER, INTENT(in) :: number
    logical, intent(in), optional :: ifSkipOrder

    IF ((number>=0).and.(number <= field_3d%number_of_fields)) THEN

      if(present(ifSkipOrder))then
        if(ifSkipOrder)then
          fu_windfield_from_3d_wind => field_3d%fields(number)%fp
          return
        endif
      endif

      IF (field_3d%in_order) THEN
        fu_windfield_from_3d_wind => &
                        & field_3d%fields(field_3d%vertical_order(number))%fp
      ELSE
        CALL set_error('this 3dfield not in order'&
                    & ,'fu_windfield_from_3d_wind')
        NULLIFY(fu_windfield_from_3d_wind)
      END IF

    ELSE

      PRINT *, number
      CALL set_error('this number field not in this 3d-field','fu_windfield_from_3d_wind')
      NULLIFY(fu_windfield_from_3d_wind)

    END IF


  END FUNCTION fu_windfield_from_3d_wind






  ! ***************************************************************

  INTEGER FUNCTION fu_number_of_windfields(field_3d)

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
    TYPE(silja_3d_windfield), POINTER :: field_3d

    IF (defined(field_3d)) THEN
      fu_number_of_windfields = field_3d%number_of_fields
    ELSE
      fu_number_of_windfields = 0
    END IF

  END FUNCTION fu_number_of_windfields




  ! ****************************************************************
  ! ****************************************************************
  ! 
  !
  !       Private functions and subroutines of this module.
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  LOGICAL FUNCTION fu_windfield_3d_defined(field_3d)

    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_3d_windfield), INTENT(in) :: field_3d

    fu_windfield_3d_defined = fu_true(field_3d%defined)

  END FUNCTION fu_windfield_3d_defined




  ! ****************************************************************

  FUNCTION fu_u_grid_of_3d_wind(windfield)

    ! Returns the grid of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_grid) :: fu_u_grid_of_3d_wind
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_3d_windfield), POINTER :: windfield

    fu_u_grid_of_3d_wind = windfield%grids(1)

  END FUNCTION fu_u_grid_of_3d_wind


  ! ****************************************************************

  FUNCTION fu_v_grid_of_3d_wind(field)

    ! Returns the grid of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_grid) :: fu_v_grid_of_3d_wind
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_3d_windfield), POINTER :: field

    fu_v_grid_of_3d_wind = field%grids(2)

  END FUNCTION fu_v_grid_of_3d_wind


  ! ****************************************************************

  FUNCTION fu_w_grid_of_3d_wind(field)

    ! Returns the grid of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_grid) :: fu_w_grid_of_3d_wind
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_3d_windfield), POINTER :: field

    fu_w_grid_of_3d_wind = field%grids(3)

  END FUNCTION fu_w_grid_of_3d_wind



  ! ***************************************************************

  INTEGER FUNCTION fu_leveltype_of_3d_wind(field_3d)

    ! Description:
    ! Returns the leveltype of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_3d_windfield), POINTER :: field_3d

    fu_leveltype_of_3d_wind = field_3d%leveltype

  END FUNCTION fu_leveltype_of_3d_wind




  ! ***************************************************************

  FUNCTION fu_met_src_of_3d_wind(field_3d) result(met_src)

    ! Returns the met_src of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    type(meteo_data_source) :: met_src
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_3d_windfield), POINTER :: field_3d

    met_src = field_3d%met_src

  END FUNCTION fu_met_src_of_3d_wind


  ! ***************************************************************

  FUNCTION fu_valid_time_of_3d_wind(field_3d)

    ! Description:
    ! Returns the analysis time of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_time) :: fu_valid_time_of_3d_wind
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_3d_windfield), POINTER :: field_3d

    fu_valid_time_of_3d_wind = field_3d%valid_time

  END FUNCTION fu_valid_time_of_3d_wind




  ! ***************************************************************

  INTEGER FUNCTION fu_3d_windfieldsize(field_3d)

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
    TYPE(silja_3d_windfield), POINTER :: field_3d
    !
    ! Local declarations:
    INTEGER :: i

    fu_3d_windfieldsize = 0

    IF (.NOT.defined(field_3d)) RETURN

    DO i = 1, field_3d%number_of_fields
      fu_3d_windfieldsize = fu_3d_windfieldsize +&
	  & fu_size(field_3d%fields(i)%fp)
    END DO


  END FUNCTION fu_3d_windfieldsize


  ! ***************************************************************

  FUNCTION fu_surf_pre_field_wind(field_3d) 

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
    TYPE(silja_field), POINTER :: fu_surf_pre_field_wind

    ! Imported parameters with intent(in):
    TYPE(silja_3d_windfield), POINTER :: field_3d

    IF (ASSOCIATED(field_3d%surface_pressure)) THEN
      fu_surf_pre_field_wind => field_3d%surface_pressure
    ELSE
      NULLIFY(fu_surf_pre_field_wind)
    END IF


  END FUNCTION fu_surf_pre_field_wind



  ! ***************************************************************

  FUNCTION fu_find_level_wind(wind_3d, level) result(windfield)

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
    TYPE(silja_windfield), POINTER :: windfield

    ! Imported parameters with intent(in):
    TYPE(silja_level), INTENT(in) :: level
    TYPE(silja_3d_windfield), POINTER :: wind_3d

    ! Local declarations:
    INTEGER :: i

    IF (.NOT.defined(wind_3d)) THEN
      CALL set_error('undefined field_3d given'&
	  & ,'fu_find_level_wind')
      RETURN
    END IF

    DO i = 1, wind_3d%number_of_fields
      windfield => wind_3d%fields(wind_3d%vertical_order(i))%fp
      IF (fu_cmp_levs_eq(fu_level(windfield), level)) RETURN
    END DO

    CALL set_error('correct level not found in 3d-field'&
	& ,'fu_find_level_wind')


  END FUNCTION fu_find_level_wind





  ! ***************************************************************

  FUNCTION fu_lowest_field_wind(field_3d)

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
    TYPE(silja_windfield), POINTER :: fu_lowest_field_wind

    ! Imported parameters with intent(in):
    TYPE(silja_3d_windfield), POINTER :: field_3d

    IF (.NOT.defined(field_3d)) THEN
      CALL set_error('undefined field_3d given'&
	  & ,'fu_lowest_field_wind')
      RETURN
    END IF

    fu_lowest_field_wind => &
	& field_3d%fields(field_3d%vertical_order(1))%fp

  END FUNCTION fu_lowest_field_wind



  ! ***************************************************************

  FUNCTION fu_highest_field_wind(field_3d)

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
    TYPE(silja_windfield), POINTER :: fu_highest_field_wind

    ! Imported parameters with intent(in):
    TYPE(silja_3d_windfield), POINTER :: field_3d

    IF (.NOT.defined(field_3d)) THEN
      CALL set_error('undefined field_3d given'&
	  & ,'fu_highest_field_wind')
      RETURN
    END IF

    fu_highest_field_wind => &
	& field_3d%fields(&
	& field_3d%vertical_order(field_3d%number_of_fields))%fp

  END FUNCTION fu_highest_field_wind




  ! ***************************************************************

  FUNCTION fu_pressures_wind(field_3d, position) result(pressures)

    ! Description:
    ! This function returns the pressures of it fields at position.
    ! Since no weather data is available to this module, calculation
    ! works only for constant pressure fields and hybrid level 3d
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
    TYPE(silja_3d_windfield), POINTER :: field_3d
    TYPE(silam_grid_position), INTENT(in) :: position

    ! Local declarations:
    INTEGER :: i
    REAL :: surface_pressure

    IF (.NOT.defined(field_3d)) THEN
      CALL set_error('undefined field_3d given'&
	  & ,'fu_pressures_wind')
      RETURN
    END IF

    SELECT CASE (field_3d%leveltype)

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
             & ,'fu_pressures_wind')
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


    CASE default
      PRINT *, field_3d%leveltype
      CALL set_error('cannot handle this type of levels here'&
	  & ,'fu_pressures_wind')

    END SELECT

  END FUNCTION fu_pressures_wind







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

    TYPE(silja_3d_windfield), POINTER :: field_3d

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

  FUNCTION fu_interpolate_wind_3d_to_pos(&
      & field_3d,&
      & position,&
      & horizontal_method,&
      & vertical_method) result(vector)

    ! Description:
    ! Interpolates both horizontally and vertically the value from a
    ! 3d wind field. For sigma- and hybrid level data the 3d_field
    ! has to contain surface pressure field for vertical coordinate
    ! hassle.
    !
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

    ! Return value of this function:
    TYPE(silja_velocity) :: vector

    ! Imported parameters with intent IN:
    TYPE(silam_grid_position), INTENT(in) :: position
    INTEGER, INTENT(in) :: horizontal_method,&
	& vertical_method ! of interpolation

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_3d_windfield), POINTER :: field_3d


    ! Local declarations: 
    REAL :: field_upwards_pressure, field_downwards_pressure
    LOGICAL :: above_the_highest, below_the_lowest
    REAL, DIMENSION(max_levels) :: field_pressures
    REAL, DIMENSION(max_levels) :: pr_differences
    TYPE(silja_windfield), POINTER :: field_upwards, field_downwards
    TYPE(silja_velocity)  :: vector_downwards, vector_upwards
    REAL :: coeff_upwards, coeff_downwards


    !----------------------------------------
    !
    ! 1. Get closest fields in vertical.
    !    ------------------------------

    CALL find_closest_windfields(&
	& position,&
	& field_3d,&
	& field_upwards,&
	& field_upwards_pressure,&
	& field_downwards,&
	& field_downwards_pressure,&
	& above_the_highest,&
	& below_the_lowest)

    IF (error) RETURN

!!!$    CALL report(field_downwards)
!!!$    CALL report(field_upwards)
!!!$    PRINT *, field_downwards_pressure, field_upwards_pressure

    ! ------------------------------------------------------
    !
    ! 2. Above the highest.


    IF (above_the_highest) THEN
      vector = fu_interpolate_to_position (&
	  & field_downwards, &
	  & position,&
	  & horizontal_method)
      RETURN
    END IF


    ! ------------------------------------------------------
    !
    ! 3. Below the lowest

    IF (below_the_lowest) THEN
      vector = fu_interpolate_to_position (&
	  & field_upwards, &
	  & position,&
	  & horizontal_method)
      RETURN
    END IF


    ! ------------------------------------------------------
    !
    ! 4. Interpolate upper value to position.

    vector_upwards = fu_interpolate_to_position (&
	& field_upwards, &
	& position,&
	& horizontal_method)

    IF (error) RETURN


    ! ------------------------------------------------------
    !
    ! 5. Interpolate lower value to position.

    vector_downwards = fu_interpolate_to_position (&
	& field_downwards, &
	& position,&
	& horizontal_method)

    IF (error) RETURN


    ! ------------------------------------------------------
    !
    ! 6. Interpolate in vertical to position's height.


    CALL weight_coefficients(&
	& (fu_pressure(position) - field_upwards_pressure), &
	& (fu_pressure(position) - field_downwards_pressure), &
	& vertical_method, &
	& coeff_upwards, coeff_downwards)

    IF (error) RETURN

    vector = (coeff_upwards * vector_upwards) + &
	& (coeff_downwards * vector_downwards)



  END FUNCTION fu_interpolate_wind_3d_to_pos



  ! ***************************************************************

  SUBROUTINE find_closest_windfields(&
      & position,&
      & field_3d,&
      & field_upwards,&
      & field_upwards_pressure,&
      & field_downwards,&
      & field_downwards_pressure,&
      & above_the_highest,&
      & below_the_lowest)

    ! Description:
    ! Inside a 3d-windfield finds the closest ones in vertical in both
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
    TYPE(silja_3d_windfield), POINTER :: field_3d !in
    TYPE(silja_windfield), POINTER :: field_upwards,& !out
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

    IF (.NOT.defined(field_3d)) THEN
      CALL set_error('undefined windfield_3d given',&
          & 'find_closest_windfields')
      RETURN
    END IF


    IF (defined(position)) THEN
      position_pressure = fu_pressure(position)
    ELSE
      CALL set_error('undefined position given',&
          & 'find_closest_windfields')
      RETURN
    END IF


    !----------------------------------------
    !
    ! 2. Get the pressure of fields at position.
    !    --------------------------------------

    field_pressures = fu_pressures(field_3d, position)
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

    field_downwards =>&
	& fu_windfield_from_3d_wind(field_3d, down_index(1))
    field_downwards_pressure = field_pressures(down_index(1))

    field_upwards =>&
	& fu_windfield_from_3d_wind(field_3d, up_index(1))
    field_upwards_pressure = field_pressures(up_index(1))


  END SUBROUTINE find_closest_windfields





  ! ***************************************************************
  ! ***************************************************************
  ! 
  !
  !                 TESTS
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  SUBROUTINE print_windfield_3d_report(field_3d)

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
    TYPE(silja_3d_windfield), POINTER :: field_3d

    ! Local declarations:
    INTEGER :: i
    TYPE(silja_windfield), POINTER :: field

    call msg(' *********** Windfield_3d report *************')

    IF (.NOT.defined(field_3d)) THEN
      call msg(' Undefined field_3d.')

    ELSE

      if(field_3d%in_order)then
        call msg(' Full report of lowest field, index=: ',field_3d%vertical_order(1))
        CALL report(field_3d%fields(field_3d%vertical_order(1))%fp)

        call msg('Wind available at following levels from bottom up:')
      DO i = 1, field_3d%number_of_fields
	field => field_3d%fields(field_3d%vertical_order(i))%fp
	CALL report(fu_level(field))
      END DO
      else
        call msg(' Full report of first field')
        CALL report(field_3d%fields(1)%fp)

        call msg('Winds are stored at following levels:')
        DO i = 1, field_3d%number_of_fields
          field => field_3d%fields(i)%fp
          CALL report(fu_level(field))
        END DO

      endif

    END IF

    IF (defined(field_3d%surface_pressure)) THEN
      call msg(' Surface pressure:')
      CALL report(field_3d%surface_pressure)
    ELSE
      call msg('No surface pressure field added to this 3d-field')
    END IF

    call msg(' *********** End-of-windfield_3d report *********')

  END SUBROUTINE print_windfield_3d_report

END MODULE windfields_3d
