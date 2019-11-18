MODULE windfields ! 2-D wind-fields with identification

  ! Description:
  ! Contains the definition of silja_windfield -type which contains one
  ! fully definedfield. 
  ! Horizontal grid or vertical level can be of any type defined in
  ! modules grids and levels.
  !
  ! Author: Mika Salonoja, FMI, email Mika.Salonoja@fmi.fi
  !
  ! All units: SI
  ! 
  ! Language: ANSI standard Fortran 90
  !
  ! 
  ! Modules used:

!  USE globals
!  USE grids_v2
!  USE toolbox
!  USE areas
  USE vectors_v2
  USE fields

  IMPLICIT NONE

  ! The public functions and subroutines available in this module:
  PUBLIC set_windfield
  PUBLIC set_windfield_from_uvw
  PUBLIC add_w_to_windfield
  PUBLIC set_windfield_empty
  PUBLIC fu_id
  PUBLIC defined
  PUBLIC fu_size
  PUBLIC fu_interpolate_to_position
  PUBLIC fu_met_src
  PUBLIC fu_u_grid
  PUBLIC fu_v_grid
  PUBLIC fu_w_grid
  PUBLIC fu_get_u_data
  PUBLIC fu_get_v_data
  PUBLIC fu_get_w_data
  PUBLIC fu_get_u_id
  PUBLIC fu_get_v_id
  PUBLIC fu_get_w_id
  PUBLIC fu_valid_time
  PUBLIC fu_analysis_time
  PUBLIC fu_forecast_length
  PUBLIC fu_level
  PUBLIC fu_w_available
  PUBLIC fu_wind_interp_cpu_usage
  PUBLIC report

  ! The private functions and subroutines:
  PRIVATE fu_windfieldsize
  PRIVATE fu_windfield_id
  PRIVATE fu_met_src_of_windfield
  PRIVATE fu_u_grid_of_windfield
  PRIVATE fu_v_grid_of_windfield
  PRIVATE fu_w_grid_of_windfield
  PRIVATE fu_valid_time_of_windfield
  PRIVATE fu_analysis_time_of_windfield
  PRIVATE fu_forecast_length_of_windfield
  PRIVATE fu_level_of_windfield
  PRIVATE fu_windfieldpointer_defined
  PRIVATE fu_interpolate_wind_to_pos
  PRIVATE print_windfield_report

  INTERFACE defined
    MODULE PROCEDURE fu_windfieldpointer_defined
  END INTERFACE

  INTERFACE report
    MODULE PROCEDURE print_windfield_report
  END INTERFACE

  INTERFACE fu_id 
    MODULE PROCEDURE fu_windfield_id
  END INTERFACE

  INTERFACE fu_size
    MODULE PROCEDURE fu_windfieldsize
  END INTERFACE

  INTERFACE fu_interpolate_to_position
    MODULE PROCEDURE fu_interpolate_wind_to_pos
  END INTERFACE

  INTERFACE fu_met_src  
    MODULE PROCEDURE fu_met_src_of_windfield
  END INTERFACE

  INTERFACE fu_u_grid 
    MODULE PROCEDURE fu_u_grid_of_windfield
  END INTERFACE

  INTERFACE fu_v_grid 
    MODULE PROCEDURE fu_v_grid_of_windfield
  END INTERFACE

  INTERFACE fu_w_grid 
    MODULE PROCEDURE fu_w_grid_of_windfield
  END INTERFACE

  INTERFACE fu_valid_time
    MODULE PROCEDURE fu_valid_time_of_windfield
  END INTERFACE

  INTERFACE fu_analysis_time
    MODULE PROCEDURE fu_analysis_time_of_windfield
  END INTERFACE

  INTERFACE fu_forecast_length
    MODULE PROCEDURE fu_forecast_length_of_windfield
  END INTERFACE

  INTERFACE fu_level
    MODULE PROCEDURE fu_level_of_windfield
  END INTERFACE

  ! Public types with public components defined in this module:
  ! NOTE. In addition to the id - it is necessary to keep u_id,v_id,w_id
  ! Reasons: 1.It is borring to distinguish between e.g. wind, wind_10m;
  ! 2.You never know the unit of w - either Pa/s, or m/s, or layer/s. 

  TYPE silja_windfield ! 2D windfields (3D windvectors on a 2D level)
    PRIVATE
    TYPE(silja_field_id) :: id
    TYPE(silja_field_id), POINTER :: u_id, v_id, w_id
    REAL, DIMENSION(:), POINTER :: u, v, w ! windfield components
    TYPE(silja_logical) :: defined
  END TYPE silja_windfield

  TYPE silja_windfieldpointer
    TYPE(silja_windfield), POINTER :: fp
  END TYPE silja_windfieldpointer

  REAL, PRIVATE :: cpu_usage = 0.

CONTAINS 

  ! ***************************************************************

  SUBROUTINE set_windfield(id, u_grid_data, v_grid_data, w_grid_data, field)
    
    ! Description:
    ! Sets a value and status for a windfield. Allocates memory for
    ! grid-data. If there's alredy something in the field, it is
    ! overwritten. In this case memory area is not deallocated and
    ! then allocated again, if the sizes of the old and new data
    ! happen to be the same.
    !
    ! ADDITION. Since windfield has now pointers to u_id, v_id, w_id, they
    ! are also filled-in here. PROBLEM - since only u-,v-,w- data are given
    ! without id for the components - it is assumed that w- is always omega-
    ! and thus has Pa/s units. NEVER use this routine if this is not the case
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Original code : Mika Salonoja
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_field_id), INTENT(in) :: id
    
    ! Imported parameters with intent INOUT or POINTER:
    REAL, DIMENSION(:), POINTER :: u_grid_data, v_grid_data, w_grid_data
    TYPE(silja_windfield), INTENT(inout) :: field
 
    ! Local declarations:
    TYPE(silja_grid), DIMENSION(3) :: grids
    INTEGER :: gridpoints, u_q_tmp, v_q_tmp, w_q_tmp

    !--------------------------------------------------------------
    !
    ! 1. Check input parameters.
    !    ---------------------

    IF (.NOT.defined(id)) THEN
      CALL set_error('identification given is not ok','set_windfield')
      RETURN
    END IF


    ! Check also the grid, since in field_id it might be undefined:
    IF (.NOT.defined(fu_grid(id))) THEN
      CALL set_error('the grid of field-id not defined','set_windfield')
      RETURN
    END IF
    
 
    !--------------------------------------------------------------
    !
    ! 2. See that there's memory for the grid-data. 
    !    -----------------------------------------
  
    ! The number of gridpoints to be stored is found in the grid:
    grids = fu_grids_of_id(id)


    !----------------------------------------------------
    !
    ! 3. U-field:
    !    -------

    IF (.NOT.ASSOCIATED(u_grid_data)) THEN
      CALL set_error('no U-field given','set_windfield')
      RETURN
    END IF

    gridpoints = fu_number_of_gridpoints(grids(1))

    IF (gridpoints > SIZE(u_grid_data)) THEN
      CALL set_error('data-array smaller than number of gridpoints', &
                   & 'set_windfield')
      RETURN
    END IF

    CALL set_array_size(field%u, gridpoints)
 
    IF (error) RETURN

    field%u = u_grid_data(1:gridpoints)



    !----------------------------------------------------
    !
    ! 4. V-field:
    !    -------
    
    IF (.NOT.ASSOCIATED(v_grid_data)) THEN
      CALL set_error('no V-field given','set_windfield')
      RETURN
    END IF
    
    gridpoints = fu_number_of_gridpoints(grids(2))
    
    IF (gridpoints > SIZE(v_grid_data)) THEN
      CALL set_error('data-array smaller than number of gridpoints','set_windfield')
      RETURN
    END IF
    
    CALL set_array_size(field%v, gridpoints)
 
    IF (error) RETURN

    field%v = v_grid_data(1:gridpoints)


    !----------------------------------------------------
    !
    ! 5. W-field:
    
    w_available: IF (ASSOCIATED(w_grid_data)) THEN

      gridpoints = fu_number_of_gridpoints(grids(3))

      IF (gridpoints > SIZE(w_grid_data)) THEN
        CALL set_error('data-array smaller than number of gridpoints','set_windfield')
        RETURN
      END IF
    
      CALL set_array_size(field%w, gridpoints)
      IF (error) RETURN
    
      field%w = w_grid_data(1:gridpoints)

    ELSE
      NULLIFY(field%w)
    
    END IF w_available


    !------------------------------------------------------------
    !
    ! 6. Set the rest of stuff.
    !    ---------------------

    field%id = id

    SELECT CASE(fu_quantity(id))

      CASE (wind_flag)
        u_q_tmp = u_flag
        v_q_tmp = v_flag
        select case(fu_leveltype(fu_level(id)))
          case(constant_altitude)
            w_q_tmp = w_alt_msl_flag
          case(constant_height)
            w_q_tmp = w_height_srf_flag
          case(constant_pressure, hybrid, sigma_level)
            w_q_tmp = omega_flag
          case default
            call msg_warning('Non-supported level','set_windfield')
            call report(id)
            call set_error('Non-supported level','set_windfield')
            return
        end select

      CASE (wind_10m_flag)
        u_q_tmp = u_10m_flag
        v_q_tmp = v_10m_flag
        w_q_tmp = int_missing

      CASE DEFAULT
        CALL set_error('unknown wind flag','set_windfield')
        return
    END SELECT

    IF(.not.ASSOCIATED(field%u_id)) ALLOCATE(field%u_id)
    IF(.not.ASSOCIATED(field%v_id)) ALLOCATE(field%v_id)

    field%u_id = id
    CALL set_quantity(field%u_id,u_q_tmp)
    field%v_id = id
    CALL set_quantity(field%v_id,v_q_tmp)

    IF(ASSOCIATED(field%w))THEN
      IF(.not.ASSOCIATED(field%w_id)) ALLOCATE(field%w_id)
      field%w_id = id
      CALL set_quantity(field%w_id,w_q_tmp)
    ELSE
      IF(ASSOCIATED(field%w_id)) DEALLOCATE(field%w_id)
      NULLIFY(field%w_id)
    END IF

    field%defined = fu_set_true()

  END SUBROUTINE set_windfield



  ! ***************************************************************

  SUBROUTINE set_windfield_from_uvw(u, v, wind, w)
    
    ! Creates a windfield from scalar wind component fields. At least
    ! u and v have to be defined and filled with data. No memory is
    ! allocated here, only pointers are set. Ths functions is used in
    ! stack to arrange uv-feld to windfields.
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Original code: Mika Salonoja
    ! Corrections: M.Sofiev, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent(inout):
    TYPE(silja_field), POINTER :: u, v
    TYPE(silja_windfield), INTENT(inout) :: wind

    ! Optional parameters with intent(inout):
    TYPE(silja_field), POINTER, OPTIONAL :: w
   

    !------------------------------------------------------------
    !
    ! 1. Set identification.
    !    ------------------
    
    IF (PRESENT(w)) THEN 
      wind%id = fu_windfield_id_from_uvw_id(fu_id(u), fu_id(v), fu_id(w))
      wind%u_id => fu_id(u)
      wind%v_id => fu_id(v)
      wind%w_id => fu_id(w)
    ELSE
      wind%id = fu_windfield_id_from_uvw_id(fu_id(u), fu_id(v))
      wind%u_id => fu_id(u)
      wind%v_id => fu_id(v)
      NULLIFY(wind%w_id)
    END IF


    !------------------------------------------------------------
    !
    ! 2. Set wind-components.
    !    -------------------
    
    IF (defined(u)) THEN
      wind%u => fu_grid_data(u)
    ELSE
      CALL set_error('undefined u given','set_windfield_from_uvw')
      RETURN
    END IF

    IF (defined(v)) THEN
      wind%v => fu_grid_data(v)
    ELSE
      CALL set_error('undefined v given','set_windfield_from_uvw')
      RETURN
    END IF

    IF (PRESENT(w)) THEN
      wind%w => fu_grid_data(w)
    ELSE
      NULLIFY(wind%w)
    END IF


    !------------------------------------------------------------
    !
    ! 3. Everything ok.
    !    --------------
    
    wind%defined = fu_set_true()

  END SUBROUTINE set_windfield_from_uvw



  ! ***************************************************************

  SUBROUTINE add_w_to_windfield(w, windfield)
    
    ! Description:
    ! Adds w to 2d windfield in case windfield hes been set earlier,
    ! and w calculated in dq-module and added. Done for HIRLAM
    ! model level stuff.
    !
    ! No checking is done here!
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE


    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_windfield), POINTER :: windfield
    TYPE(silja_field), POINTER :: w

    ! Local declarations:
    !    TYPE(silja_windfield), POINTER :: wfp
    TYPE(silja_field_id) :: id
    TYPE(silja_grid), DIMENSION(3) :: grids

    grids = fu_grids_of_id(fu_id(windfield))

    id = fu_set_windfield_id(fu_quantity(fu_id(windfield)),&
                           & fu_met_src(windfield),&
                           & fu_field_kind(fu_id(windfield)),&
                           & fu_analysis_time(windfield),&
                           & fu_forecast_length(windfield), &
                           & grids(1),&
                           & grids(2),&
                           & fu_grid(w),&
                           & fu_level(w))
    IF (error) RETURN

    windfield%id = id

    windfield%w_id => fu_id(w)
    windfield%w => fu_grid_data(w)
    
  END SUBROUTINE add_w_to_windfield




  ! ***************************************************************

  SUBROUTINE set_windfield_empty(field)
    
    ! Description:
    ! Sets windfield empty. Nullifies the pointers to the real data.
    ! It is assumed here, that windfield contains only pointers to
    ! the field-components in corresponding uvw-fields. So the memory
    ! is not deallocated here.
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
    TYPE(silja_windfield), INTENT(inout) :: field
    
    call set_missing(field%id)
    
    field%defined = fu_set_false()
    
  END SUBROUTINE set_windfield_empty



  ! ***************************************************************
  
  LOGICAL FUNCTION fu_windfieldpointer_defined(field)
    !
    ! Description:
    ! Returns a true value, if the field has been filled with data and
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    
    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_windfield), POINTER :: field

    IF (ASSOCIATED(field)) THEN    
      fu_windfieldpointer_defined = fu_true(field%defined)
    ELSE
      fu_windfieldpointer_defined = .false.
    END IF

  END FUNCTION fu_windfieldpointer_defined
  


  ! ***************************************************************

  FUNCTION fu_get_u_data(windfield)
  !
  ! Returns a pointer to the u-field data
  !
  ! Code owner Mikhail Sofiev, FMI
  !
  IMPLICIT NONE

  ! return value:
  REAL, DIMENSION(:), POINTER :: fu_get_u_data

  ! Imported parameters with intent IN
  TYPE(silja_windfield), TARGET, INTENT(in) :: windfield

  ! Local stuff
  TYPE(silja_windfield), POINTER :: wf
  
  wf => windfield

    IF(.not.defined(wf))THEN
      CALL set_error('undefined field','fu_get_u_data')
      NULLIFY(fu_get_u_data)
      RETURN
    END IF

    fu_get_u_data => windfield%u

  END FUNCTION fu_get_u_data



  ! ***************************************************************

  FUNCTION fu_get_v_data(windfield)
  !
  ! Returns a pointer to the v-field data
  !
  ! Code owner Mikhail Sofiev, FMI
  !
  IMPLICIT NONE

  ! Return value
  REAL, DIMENSION(:), POINTER :: fu_get_v_data

  ! Imported parameters with intent IN
  TYPE(silja_windfield), POINTER :: windfield

    IF(.not.defined(windfield))THEN
      CALL set_error('undefined field','fu_get_v_data')
      NULLIFY(fu_get_v_data)
      RETURN
    END IF

    fu_get_v_data => windfield%v

  END FUNCTION fu_get_v_data



  ! ***************************************************************

  FUNCTION fu_get_w_data(windfield)
  !
  ! Returns a pointer to the w-field data
  !
  ! Code owner Mikhail Sofiev, FMI
  !
  IMPLICIT NONE

  TYPE(silja_windfield), POINTER :: windfield
  REAL, DIMENSION(:), POINTER :: fu_get_w_data

    IF(.not.defined(windfield))THEN
      CALL set_error('undefined field','fu_get_w_data')
      NULLIFY(fu_get_w_data)
      RETURN
    END IF

    IF(ASSOCIATED(windfield%w))THEN
      fu_get_w_data => windfield%w
    ELSE
      CALL set_error('undefined field','fu_get_w_data')
      NULLIFY(fu_get_w_data)
      RETURN
    END IF

  END FUNCTION fu_get_w_data




  ! ***************************************************************

  FUNCTION fu_get_u_id(windfield)
  !
  ! Returns a pointer to the u-field id
  !
  ! Code owner Mikhail Sofiev, FMI
  !
  IMPLICIT NONE

  ! return value:
  TYPE(silja_field_id), POINTER :: fu_get_u_id

  ! Imported parameters with intent IN
  TYPE(silja_windfield), TARGET, INTENT(in) :: windfield

  ! Local stuff
  TYPE(silja_windfield), POINTER :: wf
  
  wf => windfield

    IF(.not.defined(wf))THEN
      CALL set_error('undefined field','fu_get_u_id')
      NULLIFY(fu_get_u_id)
      RETURN
    END IF

    fu_get_u_id => windfield%u_id

  END FUNCTION fu_get_u_id



  ! ***************************************************************

  FUNCTION fu_get_v_id(windfield)
  !
  ! Returns a pointer to the v-field data
  !
  ! Code owner Mikhail Sofiev, FMI
  !
  IMPLICIT NONE

  ! Return value
  TYPE(silja_field_id), POINTER :: fu_get_v_id

  ! Imported parameters with intent IN
  TYPE(silja_windfield), POINTER :: windfield

    IF(.not.defined(windfield))THEN
      CALL set_error('undefined field','fu_get_v_id')
      NULLIFY(fu_get_v_id)
      RETURN
    END IF

    fu_get_v_id => windfield%v_id

  END FUNCTION fu_get_v_id



  ! ***************************************************************

  FUNCTION fu_get_w_id(windfield)
  !
  ! Returns a pointer to the w-field data
  !
  ! Code owner Mikhail Sofiev, FMI
  !
  IMPLICIT NONE

  TYPE(silja_windfield), POINTER :: windfield
  TYPE(silja_field_id),POINTER :: fu_get_w_id

    IF(.not.defined(windfield))THEN
      CALL set_error('undefined field','fu_get_w_id')
      NULLIFY(fu_get_w_id)
      RETURN
    END IF

    IF(ASSOCIATED(windfield%w))THEN
      fu_get_w_id => windfield%w_id
    ELSE
      CALL set_error('undefined field','fu_get_w_id')
      NULLIFY(fu_get_w_id)
      RETURN
    END IF

  END FUNCTION fu_get_w_id



  ! ***************************************************************
  ! ***************************************************************
  !
  !
  !      Return some basic stuff from identification
  !
  !
  ! ***************************************************************
  ! ***************************************************************
  
  FUNCTION fu_u_grid_of_windfield(field)

    ! Returns the grid of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_grid) :: fu_u_grid_of_windfield
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_windfield), POINTER :: field

    fu_u_grid_of_windfield = fu_u_grid(field%id)
    
  END FUNCTION fu_u_grid_of_windfield



  ! ***************************************************************
  
  FUNCTION fu_v_grid_of_windfield(field)

    ! Returns the grid of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_grid) :: fu_v_grid_of_windfield
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_windfield), POINTER :: field

    fu_v_grid_of_windfield = fu_v_grid(field%id)
    
  END FUNCTION fu_v_grid_of_windfield



  ! ***************************************************************
  
  FUNCTION fu_w_grid_of_windfield(field)

    ! Returns the grid of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_grid) :: fu_w_grid_of_windfield
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_windfield), POINTER :: field

    fu_w_grid_of_windfield = fu_w_grid(field%id)
    
  END FUNCTION fu_w_grid_of_windfield



  ! ***************************************************************
  
  FUNCTION fu_met_src_of_windfield(field) result(met_src)

    ! Returns the met_src of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    type(meteo_data_source) :: met_src
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_windfield), POINTER :: field

    met_src = fu_met_src(field%id)
    
  END FUNCTION fu_met_src_of_windfield



  ! ***************************************************************
  
  FUNCTION fu_level_of_windfield(field)
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
    TYPE(silja_level) :: fu_level_of_windfield
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_windfield), POINTER :: field
  
    fu_level_of_windfield = fu_level(field%id)
    
  END FUNCTION fu_level_of_windfield




  ! ***************************************************************
  
  FUNCTION fu_valid_time_of_windfield(field)
    
    ! Description:
    ! Returns the analysis time of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_time) :: fu_valid_time_of_windfield
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_windfield), POINTER :: field
    
    fu_valid_time_of_windfield = fu_valid_time(field%id)
    
  END FUNCTION fu_valid_time_of_windfield



  ! ***************************************************************
  
  FUNCTION fu_analysis_time_of_windfield(field)
    
    ! Description:
    ! Returns the analysis time of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_time) :: fu_analysis_time_of_windfield
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_windfield), POINTER :: field
    
    fu_analysis_time_of_windfield = fu_analysis_time(field%id)
    
  END FUNCTION fu_analysis_time_of_windfield



  ! ***************************************************************
  
  FUNCTION fu_forecast_length_of_windfield(field)
    
    ! Description:
    ! Returns the analysis time of a field-id.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_interval) :: fu_forecast_length_of_windfield
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_windfield), POINTER :: field
    
    fu_forecast_length_of_windfield = fu_forecast_length(field%id)
    
  END FUNCTION fu_forecast_length_of_windfield



  ! ***************************************************************

  FUNCTION fu_windfield_id(field) result(id)
    
    ! Description:
    ! Returns the identification section of a field.
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_field_id) :: id
    !
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_windfield), POINTER :: field
    
    id = field%id
    
  END FUNCTION fu_windfield_id

  

  ! ***************************************************************

  INTEGER FUNCTION fu_windfieldsize(field)
    
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
    TYPE(silja_windfield), POINTER :: field

    IF (defined(field)) THEN
      IF (ASSOCIATED(field%w)) THEN
	fu_windfieldsize =&
	    & SIZE(field%u) + SIZE(field%u) +SIZE(field%w)
      ELSE
	fu_windfieldsize =&
	    & SIZE(field%u) + SIZE(field%u)
      END IF

    ELSE
      fu_windfieldsize = 0
    END IF
    
  END FUNCTION fu_windfieldsize
  


  ! ***************************************************************

  LOGICAL FUNCTION fu_w_available(windfield)
    
    ! Description:
    ! Returns true value if windfield also contains vertical velocity
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_windfield), POINTER :: windfield
 
    fu_w_available = ASSOCIATED(windfield%w)

  END FUNCTION fu_w_available





  ! ****************************************************************
  ! 
  !
  !           INTERPOLATION
  !
  !
  ! ***************************************************************
  
  FUNCTION fu_interpolate_wind_to_pos (&
      & field, &
      & position,&
      & method) result(velocity)
    
    ! Description:
    ! Interpolates the field to given position and returns the wind
    ! as a velocity-vector.
    ! 
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Return value of this function:
    TYPE(silja_velocity) :: velocity

    ! Imported parameters with intent(in):
    TYPE(silam_grid_position), INTENT(in) :: position
    INTEGER, INTENT(in) :: method ! of interpolation
    
    ! Imported parameters with intent(inout) or POINTER:
    TYPE(silja_windfield), POINTER :: field

    ! Local declarations:
    REAL, DIMENSION(:), POINTER :: grid_data_pointer
    REAL :: u, v, w
    REAL :: before
    INTEGER :: local_method 

!    before = fu_total_cpu_time_sec()

    ! Interpolate the N-S component
    u = fu_2d_interpolation(field%u, fu_u_grid(field%id), position, method)
    IF (error) THEN
      PRINT *, ' % error in interpolating UU component %',u
!!!$      local_method = nearest_point
!!!$      u = fu_2d_interpolation(&
!!!$          & field%u, &
!!!$          & fu_u_grid(field%id),&
!!!$          & position,&
!!!$          & local_method)
!!!$      CALL unset_error('% Am trying once more % ') 
!!!$      PRINT *, ' % nearest point value =  %',u
      
      RETURN
      
    END IF
    
    ! Interpolate the E-W component
    v = fu_2d_interpolation(field%v, fu_v_grid(field%id), position, method)
    IF (error) THEN
      RETURN
    END IF
    
    ! Interpolate/fetch the vertical component
    IF (ASSOCIATED(field%w)) THEN
      w = fu_2d_interpolation(field%w, fu_w_grid(field%id), position, method)
      IF (error) THEN
        PRINT *, ' % error in interpolating WW component %'      
        CALL unset_error('%  set to zero ! % ') 
        w= 0.
      END IF
    ELSE
      w= 0.
    END IF
 
    velocity = fu_set_velocity(u, v, w, fu_u_grid(field%id), position)

!    cpu_usage = cpu_usage + fu_total_cpu_time_sec() - before
    
  END FUNCTION fu_interpolate_wind_to_pos

  
  ! ***************************************************************

  FUNCTION fu_wind_interp_cpu_usage()
    
    ! Description:
    ! Returns the total cpu-time used by the functions of this module
    ! during the run of program.
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
    TYPE(silja_interval) :: fu_wind_interp_cpu_usage

    fu_wind_interp_cpu_usage = fu_set_interval_sec(cpu_usage)

  END FUNCTION fu_wind_interp_cpu_usage




  ! *****************************************************************
  ! ****************************************************************
  ! 
  !
  !        Tests
  !
  !
  ! ***************************************************************
  ! ***************************************************************
  

  SUBROUTINE print_windfield_report(field)
    !
    ! Description:
    ! Prints a report of a field to screen. All the identifications,
    ! status values and max/min values are printed.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_windfield), POINTER :: field
    !
    ! Local declarations:
    INTEGER :: i, n = 0
    
    call msg(' %%%%%%%%%%%%%%%% Printing a windfield report %%%%%%%%%%%%%%%%')
    
    if(associated(field))then
    IF (.NOT.defined(field)) THEN
        call msg('Undefined field')
      RETURN
    END IF
    else
      call msg('Unassociated wind field')
      return
    endif

    CALL report(field%id)

    call msg(' Field size: ', fu_windfieldsize(field))

    PRINT *, ' U min and max: ', MINVAL(field%u), MAXVAL(field%u)

    PRINT *, ' V min and max: ', MINVAL(field%v), MAXVAL(field%v)

    IF (ASSOCIATED(field%w)) THEN    
      PRINT *, ' W min and max: ', MINVAL(field%w), MAXVAL(field%w)
    ELSE
      PRINT *, 'No W in this windfield.'
    END IF

    DO i = 1, SIZE(field%u)
      IF (field%u(i) .eps. real_missing) n = n+1
    END DO

    DO i = 1, SIZE(field%v)
      IF (field%v(i) .eps. real_missing) n = n+1
    END DO

    IF (ASSOCIATED(field%w)) THEN
      DO i = 1, SIZE(field%u)
	IF (field%w(i) .eps. real_missing) n = n+1
      END DO
    END IF

    PRINT *, ' number of missing values: ', n

    PRINT *, '%%%%%%%%% end of windfield report %%%%%%%%%%%%%%%%%%%%'
    

  END SUBROUTINE print_windfield_report

END MODULE windfields
