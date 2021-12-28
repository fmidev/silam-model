MODULE vectors_v2 ! velocity- and movement in x,y,z,t

  ! Description:
  ! Contains the definition of the public types silja_movement and
  ! silja_velocity.
  !
  ! Silja_movement defines a 3D + time-interval movement from a
  ! position to a new position. The movement is NOT a wind-vector
  ! (speed components) but a movement in relative grid space.
  !
  ! Silja_velocity defines a 3D-velocity vector in relative-units.
  ! 
  ! In this type we define:
  ! MOVEMENT FROM SOUTH TO NORTH IS POSITIVE
  ! MOVEMENT FROM WEST TO EAST IS POSITIVE
  !
  ! This module provides tools for setting a value for a movement and
  ! velocity, plus a toolkit for them.
  !
  ! MOVEMENTS AND VELOCITIES ARE IN RELATIVE UNITS IN THE SYSTEM_GRID.
  ! To be exact, horizontal movement is relative in the grid, while
  ! vertical is so far kept as Pa. Corresponding units for velocities:
  ! grid_cell/sec and Pa/sec.
  ! Actually it makes the most of below routines void, but in order 
  ! to keep the incapsulation they are not kicked out.
  !
  ! Original code: Mika Salonoja
  ! Author: Mikhail Sofiev, FMI
  ! 
  ! All units: RELATIVE, not SI
  ! 
  ! All degrees: real numbers in degrees; NOT in degrees and minutes
  ! 
  ! Language: ANSI standard Fortran 90
  ! 
  ! ***************************************************************
  !
  ! Modules used:
!  USE globals
!  USE times
!  USE poles
!  USE positions_v2
  use fields
!  USE grids

  IMPLICIT NONE

  ! The public functions and subroutines available in this module:
  PUBLIC fu_set_movement
  PUBLIC defined
  PUBLIC fu_compare_movements_eq ! interface ==
  PUBLIC fu_add_position_and_movement ! interface +
  PUBLIC fu_horizontal_length
  PUBLIC movement_tests
  PUBLIC fu_multiply_movement ! interface *
  PUBLIC fu_add_movements ! interface +
  PUBLIC report
  PUBLIC fu_subtract_positions ! interface - 
  PUBLIC fu_subs_position_and_movement ! interface - 
  PUBLIC fu_subtract_movements ! interface -
  PUBLIC fu_set_velocity
  PUBLIC fu_multiply_velocity1 ! velocity = real * velocity;
  ! interface *
  PUBLIC fu_multiply_velocity2 ! movement = velocity * interval;
  ! interface *
  PUBLIC fu_add_velocities ! interface +
  PUBLIC fu_speed
  PUBLIC fu_direction_deg
  PUBLIC scale_movement
  
  ! The private functions and subroutines not to be used elsewhere:
  PRIVATE turn_latlon_to_latlon
  PRIVATE turn_dxdy_to_latlon
  PRIVATE fu_velocity_defined
  PRIVATE fu_movement_defined
  PRIVATE fu_set_velocity_by_pole
  PRIVATE fu_set_velocity_by_grid
  PRIVATE fu_set_movement_by_pole
  PRIVATE fu_set_movement_by_grid
  private turndd
  PRIVATE print_velocity
  private print_movement


  ! Generic names and operator-interfaces of some functions:

  INTERFACE defined
    MODULE PROCEDURE fu_velocity_defined, fu_movement_defined
  END INTERFACE

  INTERFACE report
    MODULE PROCEDURE print_velocity, print_movement
  END INTERFACE

  INTERFACE fu_set_velocity
    MODULE PROCEDURE fu_set_velocity_by_pole, fu_set_velocity_by_grid
  END INTERFACE

  INTERFACE fu_set_movement
    MODULE PROCEDURE fu_set_movement_by_pole, fu_set_movement_by_grid
  END INTERFACE

  INTERFACE operator(==)
    MODULE PROCEDURE fu_compare_movements_eq
  END INTERFACE
  
  INTERFACE operator(+)
    MODULE PROCEDURE &
	& fu_add_position_and_movement,&
	& fu_add_movements,&
	& fu_add_velocities
  END INTERFACE
  
  INTERFACE operator(-)
    MODULE PROCEDURE &
	& fu_subs_position_and_movement,&
	& fu_subtract_movements,&
	& fu_subtract_positions
  END INTERFACE
  
  INTERFACE operator(*)
    MODULE PROCEDURE &
	& fu_multiply_movement,&
	& fu_multiply_velocity1,&
	& fu_multiply_velocity2
  END INTERFACE
  
   ! Public types with private components defined in this module:
  TYPE silja_movement 
    PRIVATE
    REAL :: dx, dy ! relative in system_grid
    REAL ::  dp ! in Pa
    TYPE(silja_interval) :: dt
    TYPE(silja_logical) :: defined
  END TYPE silja_movement
  
  TYPE silja_velocity
    PRIVATE
    TYPE(silja_movement) :: mm
  END TYPE silja_velocity
  
  ! Missing codes: 
  TYPE(silja_movement), PARAMETER, PUBLIC :: movement_missing = &
      & silja_movement(real_missing, real_missing, real_missing,&
      & interval_missing, silja_false)

  TYPE(silja_velocity), PARAMETER, PUBLIC :: velocity_missing = &
      & silja_velocity(movement_missing)
  
  ! Zero values: 
  TYPE(silja_movement), PARAMETER, PUBLIC :: zero_movement = &
      & silja_movement(0., 0., 0., zero_interval, silja_true)

  TYPE(silja_movement), PARAMETER, PRIVATE :: movement_zero_velo = &
      & silja_movement(0., 0., 0., one_second, silja_true)

  TYPE(silja_velocity), PARAMETER, PUBLIC :: zero_velocity = &
      & silja_velocity(movement_zero_velo)



CONTAINS

  ! ***************************************************************
     
  FUNCTION fu_set_movement_by_pole(dx, dy, pole, position, dp, dt) result(mov)
    !
    ! Sets a value for a movement. The horizontal components are
    ! along the axis defined by the pole. Be careful with the units!
    !
    ! Remember:
    ! MOVEMENT FROM SOUTH TO NORTH IS POSITIVE
    ! MOVEMENT FROM WEST TO EAST IS POSITIVE
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    
    ! The return value of this function:
    TYPE(silja_movement) :: mov
    !
    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: dx ! west-east, in metres
    REAL, INTENT(in) :: dy ! south-north, in metres
    TYPE(silam_pole), INTENT(in) :: pole ! pole to define the axis,
    ! along which the dx (west-east) and dy (south-north) are.
    TYPE(silam_grid_position), INTENT(in) :: position ! from which
    ! the dx and dy are
    REAL, INTENT(in) ::  dp ! Pa 
    TYPE(silja_interval), INTENT(in) :: dt

    IF (.NOT.defined(pole)) THEN
      CALL set_error('undefined pole', 'fu_set_movement_by_pole')
      RETURN
    END IF

    ! Set dx and dy, turn if necessary:
    IF (pole == dispersion_pole) THEN
      ! the movement components already in the right lat-lon-axis:
      !      PRINT *, ' System pole'
      mov%dx = dx
      mov%dy = dy
    ELSE ! the movement components have to be turned:

!!!$      PRINT *, ' NOT System pole *********'
!!!$      PRINT *,   fu_southpole_lat(pole),  fu_southpole_lon(pole)
!!!$      PRINT *,   fu_southpole_lat(pole_system),&
!!!$	  & fu_southpole_lon(pole_system)

      IF(.NOT.defined(position)) THEN
        CALL set_error('undefined position', 'fu_set_movement_by_pole')
        RETURN
      END IF
      
      CALL turn_latlon_to_latlon (&
           & dx, dy, position, &
           & pole,& ! the pole that defines the axis of dx and dy
           & dispersion_pole, & ! the pole that defines the axis which dx
                            ! and dy have to be converted to
           & mov%dx, mov%dy)
      
      IF (error) THEN
        mov = movement_missing
        RETURN
      END IF

    END IF
 
    ! set the rest:
    mov%dt = dt
    mov%dp = dp

    ! Everything ok:
    mov%defined = fu_set_true()

  END FUNCTION fu_set_movement_by_pole


  ! ***************************************************************
  
  FUNCTION fu_set_movement_by_grid (dx, dy, grid, position, dp, dt) &
      & result(movement) 
    
    ! Description:
    ! Sets a value for a movement-vector. The given u-v-vomponents
    ! have to be in the coordinate system defined by the grid.
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
    TYPE(silja_movement) :: movement
    
    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: dx, dy ! west-east, south-north, relative units
    TYPE(silja_grid), INTENT(in) :: grid ! to define the axis,
    ! along which the dx (west-east) and dy (south-north) are.
    TYPE(silam_grid_position), INTENT(in) :: position ! from which
    ! the dx and dy are
    REAL, INTENT(in) ::  dp ! Pa 
    TYPE(silja_interval), INTENT(in) :: dt
    
    IF (fu_position_inside_grid(position, grid) /= iPositionInside) THEN
      movement = movement_missing
      RETURN
    END IF

    SELECT CASE (fu_gridtype(grid))

    CASE (lonlat, anygrid)

      movement = fu_set_movement_by_pole(dx, dy, fu_pole(grid),&
	  & position, dp, dt)

    CASE default

      CALL set_error('sorry, can only handle latlon grids'&
                  & ,'fu_set_movement_by_grid')

    END SELECT

  END FUNCTION fu_set_movement_by_grid


  ! ***************************************************************

  LOGICAL FUNCTION fu_movement_defined(movement)
    
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_movement), INTENT(in) :: movement
    
    fu_movement_defined = fu_true(movement%defined)
    
  END FUNCTION fu_movement_defined
  
   

  ! ***************************************************************
  
  FUNCTION fu_subs_position_and_movement(pos, mov) result(pos_out)

    ! Subtracts a movement from a position, that is: adds the
    ! movement-vector of same length but opposite direction and time
    ! -interval.
    ! The handling of mapfactors is taken care of (the conversion of
    ! horizontal movement from meters to degrees).
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    
    TYPE(silam_grid_position) :: pos_out ! The return value of this
    ! function
    
    ! Imported parameters with intent IN:
    TYPE(silam_grid_position), INTENT(in) :: pos
    TYPE(silja_movement), INTENT(in) :: mov
  
    pos_out = pos + ((-1.) * mov)  

  END FUNCTION fu_subs_position_and_movement
  


  ! ***************************************************************
  
  FUNCTION fu_add_position_and_movement(pos, mov) result(pos_out)
    
    ! Adds a movement to a position. Also the time-interval is added.
    ! The handling of mapfactors is taken care of (the conversion of
    ! horizontal movement from meters to degrees).
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Original code: Mika Salonoja
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE
    
    TYPE(silam_grid_position) :: pos_out ! The return value 
    
    ! Imported parameters with intent IN:
    TYPE(silam_grid_position), INTENT(in) :: pos
    TYPE(silja_movement), INTENT(in) :: mov
    
    ! Local declarations:
    TYPE(silja_movement) :: piece
    TYPE(silam_grid_position) :: pos1, pos2
    REAL :: length
    INTEGER :: splitnumber, i

    ! Local parameters:
    REAL, PARAMETER :: too_long = 2. ! longer movements have to
    ! be split because of narrowing of longitude-axis in spherical
    ! coordinates
    
    !---------------------------------------------------------
    !
    ! 1. Add short movement.
    !    -------------------

    IF (.NOT.defined(pos)) THEN
      CALL set_error('undefined position given'&
	  & ,'fu_add_position_and_movement')
      RETURN
    END IF

    IF (.NOT.defined(mov)) THEN
      CALL set_error('undefined movement given'&
	  & ,'fu_add_position_and_movement')
      RETURN
    END IF

    length =  fu_horizontal_length(mov)

    IF (error) RETURN

    IF (length < too_long) THEN
      
      pos_out = fu_set_position(&
         & (fu_x(pos) + mov%dx),&
         & (fu_y(pos) + mov%dy),&
         & (fu_pressure(pos) + mov%dp), &
         & (fu_time(pos) + mov%dt))
      
    ELSE
      call report(mov)
      CALL msg_warning('Very long movement','add_position_and_movement')
      !
      ! Can not make a split because do not know the irregularity of the
      ! system grid

!      !---------------------------------------------------------
!      !
!      ! 2. Add long movement by splitting it into short movements.
!      !    ------------------------------------------------------
!
!      splitnumber = INT(length/too_long) + 2
!      piece = (1./REAL(splitnumber)) * mov
!      pos1 = pos
!     
!      DO i = 1, splitnumber
!        pos2 = fu_add_position_and_short(pos1, piece)
!        IF (error) RETURN
!        pos1 = pos2
!      END DO
!
!      pos_out = pos2
!      
    END IF

  END FUNCTION fu_add_position_and_movement


!  ! ***************************************************************
!  
!  FUNCTION fu_add_position_and_short(pos, mov)&
!      & result(pos_out)
!    
!    ! Description:
!    ! Adds a short movement to a position. Also the time-interval is
!    ! added.
!    ! The handling of mapfactors is taken care of (the conversion of
!    ! horizontal movement from meters to degrees).
!    !
!    ! All units: SI
!    !
!    ! Language: ANSI Fortran 90
!    !
!    ! Author: Mika Salonoja, FMI
!    ! 
!    IMPLICIT NONE
!    
!    TYPE(silja_position) :: pos_out ! The return value of this
!    ! function
!    
!    ! Imported parameters with intent IN:
!    TYPE(silja_position), INTENT(in) :: pos
!    TYPE(silja_movement), INTENT(in) :: mov
!    
!    ! Local declarations:
!    REAL :: dx_deg, dy_deg
!    REAL :: lat_system, lon_system
!    
!    !---------------------------------------------------------
!    !
!    ! 1. Convert the x-y-movements from meters to relative units
!    !    -------------------------------------------------
!    
!    CALL get_lonlat(pos, pole_system, lon_system, lat_system)
!    CALL distance_metres_to_degrees(mov%dx, mov%dy, lat_system,&
!                                  & dx_deg, dy_deg)
!    dx_relative = mov%dx
!
!    !---------------------------------------------------------
!    !
!    ! 2. Set the new position.
!    !    ---------------------
!    
!    pos_out = fu_set_position_modified(&
!       & (lat_system + dy_deg), north, &
!       & (lon_system + dx_deg), east, &
!       & pole_system, &
!       & (fu_pressure(pos) + mov%dp), &
!       & (fu_time(pos) + mov%dt))
!    
!  END FUNCTION fu_add_position_and_short


  ! ***************************************************************
  
  FUNCTION fu_add_movements(mov1, mov2) result(mov)
    !
    ! Description:
    ! Adds two movements. Also the time-interval is added.
    !
    IMPLICIT NONE
    
    TYPE(silja_movement) :: mov ! The return value of this
    ! function
    
    ! Imported parameters with intent IN:
    TYPE(silja_movement), INTENT(in) :: mov1, mov2

    mov%dx = mov1%dx + mov2%dx
    mov%dy = mov1%dy + mov2%dy
    mov%dp = mov1%dp + mov2%dp
    mov%dt = mov1%dt + mov2%dt
    mov%defined = fu_set_true()

  END FUNCTION fu_add_movements
  

  ! ***************************************************************
  
  FUNCTION fu_subtract_movements(mov1, mov2) result(mov)

    ! Description:
    ! Subtracts two movements. Also the time-interval is
    ! subtracted. If the direction of the vectors is same, and mov1
    ! is longer, then the result also points to the same direction as
    ! the originals (or the result is positive).
    !
    IMPLICIT NONE
    
    TYPE(silja_movement) :: mov ! The return value of this
    ! function
    
    ! Imported parameters with intent IN:
    TYPE(silja_movement), INTENT(in) :: mov1, mov2

    mov = mov1 + ((-1.) * mov2) 
    mov%defined = fu_set_true()
    
  END FUNCTION fu_subtract_movements
  

  ! ***************************************************************

  FUNCTION fu_multiply_movement(x, movement) result(new_movement)
    
    ! Description:
    ! Returns a vector multiplied by factor x. All the components are
    ! handled.
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_movement) :: new_movement
    !
    ! Imported parameters with intent(in):
    TYPE(silja_movement), INTENT(in) :: movement
    REAL, INTENT(in) :: x

    new_movement%dx = movement%dx * x
    new_movement%dy = movement%dy * x
    new_movement%dp = movement%dp * x
    new_movement%dt = movement%dt * x
    new_movement%defined = fu_set_true()
   
  END FUNCTION fu_multiply_movement


  ! ***************************************************************
  
  LOGICAL FUNCTION fu_compare_movements_eq(mov1, mov2)

    ! Description:
    ! Compares two movements and returns a ture value if equal. This
    ! function is used through an interface (mov1 == mov2)
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    TYPE(silja_movement), INTENT(in) :: mov1, mov2
    !
    IF (&
	&(mov1%dx == mov2%dx) .and. &
	&(mov1%dy == mov2%dy) .and. &
	&(mov1%dp == mov2%dp)) THEN
      fu_compare_movements_eq = .true.
    ELSE
      fu_compare_movements_eq = .false.
    END IF
    
  END FUNCTION fu_compare_movements_eq
  

  ! ***************************************************************

  FUNCTION fu_subtract_positions(pos2, pos1) result(mov)
    
    ! Description:
    ! Calculates the vector from position1 to position2. The result
    ! is of  type silja_movement, so that position1 + movement =
    ! position2.
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
    TYPE(silja_movement) :: mov
    !
    ! Imported parameters with intent(in):
    TYPE(silam_grid_position), INTENT(in) :: pos1, pos2
    
    ! Local declarations:
    ! 
    !----------------------------------------
    !
    ! 1. 
    !
    CALL set_error('sorry not working','fu_subtract_positions')

    mov = movement_missing

  END FUNCTION fu_subtract_positions


  ! ****************************************************************
  ! ****************************************************************
  ! 
  !
  !          LET'S TURN SOME VECTORS (OH SHIT!)
  !
  !
  ! ***************************************************************
  ! ***************************************************************
   
  SUBROUTINE turn_latlon_to_latlon & 
      & (dx1, dy1, position, pole1, &
      & pole2, dx2, dy2)

    ! Description:
    ! Turns a vector between two lat-lon systems. The two spherical
    ! coordinate systems are defined by the two poles given.
    !
    ! The turned vector can be a velocity- or movement-vector, and
    ! no unit conversion are made.
    !
    ! The position in the new lat-lon system is also calculated.
    !
    ! Northern latitudes and eastern longitudes are positive.
    ! Movement from south to north and from west to east is positive.
    !
    ! Method:
    ! TURNDD-subroutine from the FMI-HIRLAM HILA-library is used for
    ! rotation.
    ! 
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    REAL, INTENT(in) :: dx1, dy1 ! the original vector components
    TYPE(silam_grid_position), INTENT(in) :: position ! from where dx1 and
                                                      ! dx2 are
    TYPE(silam_pole), INTENT(in) :: pole1 ! the pole that defines the
                     ! lat-lon-axis along which the dx and dy are
    TYPE(silam_pole), INTENT(in) :: pole2  ! the pole that defines
    ! the lat-lon-system into which the dx and dy are converted
    
    ! Imported parameters with intent OUT:
    REAL, INTENT(inout) :: dx2, dy2 ! the new vector components
    ! in lat-lon-grid defined by the pole2

    ! Local declarations:
    REAL :: lat1, lon1, lat2, lon2
    REAL :: pa, pb, pc, pd ! the transformation coefficients
    !TYPE(silam_pole), SAVE :: pole1_previous, pole2_previous
    REAL :: dx_geo, dy_geo, geolat, geolon

    ! External functions used:
!    EXTERNAL turndd ! FMI HIRLAM HILA-library


    !----------------------------------------------------------------
    !
    ! 1. Get the coordinates of postion in relation to the two poles.
    !    ------------------------------------------------------------
    
    CALL get_lonlat(position, pole1, lon1, lat1)
    
    CALL get_lonlat(position, pole2, lon2, lat2)
      
    IF (error) RETURN

    !----------------------------------------
    !
    ! 2. Turn from geographical to modified.
    !    -----------------------------------
    
    IF (pole1 == pole_geographical) THEN
      
      ! 1.1. Get the lat-lon coordinates of position in respect to
      ! the two given poles:
      
      CALL turndd(dx1, dy1, &
                & dx2, dy2, &
                & pa, pb, pc, pd, &
                & lon1, lat1, & ! geographical
                & lon2, lat2, & ! modified
                & fu_southpole_lon(pole2), & ! the modified pole
                & fu_southpole_lat(pole2), &
                & 2) ! geographical to modified, preparation

      CALL turndd(dx1, dy1, &
                & dx2, dy2, &
                & pa, pb, pc, pd, &
                & lon1, lat1, & ! geographical
                & lon2, lat2, & ! modified
                & fu_southpole_lon(pole2), & ! the modified pole
                & fu_southpole_lat(pole2), &
                & 1) ! geographical to modified, calculation

!!!$      PRINT *, dx1, dy1, &
!!!$	  & dx2, dy2, &
!!!$	  & lon1, lat1, & ! geographical in
!!!$	  & lon2, lat2, & ! modified out
!!!$	  & fu_southpole_lon(pole2), & ! the modified pole
!!!$	  & fu_southpole_lat(pole2)

    ELSE IF (pole2 == pole_geographical) THEN

      !----------------------------------------
      !
      ! 3. Turn from modified to geographical
      !    ----------------------------------
      
      CALL turndd(dx1, dy1, &
                & dx2, dy2, &
                & pa, pb, pc, pd, &
                & lon2, lat2, & ! geographical
                & lon1, lat1, & ! modified
                & fu_southpole_lon(pole1), & ! the modified pole
                & fu_southpole_lat(pole1), &
                & -2) ! modified to geographical, preparation
  
      CALL turndd(dx1, dy1, &
                & dx2, dy2, &
                & pa, pb, pc, pd, &
                & lon2, lat2, & ! geographical
                & lon1, lat1, & ! modified
                & fu_southpole_lon(pole1), & ! the modified pole
                & fu_southpole_lat(pole1), &
                & -1) ! modified to geographical, calculation
    ELSE
      
      !----------------------------------------
      !
      ! 4. Turn from modified1 to modified2.
      !    ---------------------------------
      
      ! 4.1. The geographical coordnates of position.

      CALL get_lonlat(position, pole_geographical, geolon, geolat)

      IF (error) RETURN

      ! 4.1.  Turn from modified1 to geographical

      CALL turndd(dx1, dy1, & ! in 
                & dx_geo, dy_geo, & ! out
                & pa, pb, pc, pd, &
                & geolon, geolat, & ! geographical
                & lon1, lat1, & ! modified
                & fu_southpole_lon(pole1), & ! the modified pole
                & fu_southpole_lat(pole1), &
                & -2) ! modified to geographical, preparation
      
      CALL turndd(dx1, dy1, &
                & dx_geo, dy_geo, &
                & pa, pb, pc, pd, &
                & geolon, geolat, & ! geographical
                & lon1, lat1, & ! modified
                & fu_southpole_lon(pole1), & ! the modified pole
                & fu_southpole_lat(pole1), &
                & -1) ! modified to geographical, calculation


      ! 4.2.  Turn from geographical to modified2.

      CALL turndd(dx_geo, dy_geo, &
                & dx2, dy2, &
                & pa, pb, pc, pd, &
                & geolon, geolat, & ! geographical
                & lon2, lat2, & ! modified
                & fu_southpole_lon(pole2), & ! the modified pole
                & fu_southpole_lat(pole2), &
                & 2) ! geographical to modified, preparation

      CALL turndd(dx_geo, dy_geo, &
                & dx2, dy2, &
                & pa, pb, pc, pd, &
                & geolon, geolat, & ! geographical
                & lon2, lat2, & ! modified
                & fu_southpole_lon(pole2), & ! the modified pole
                & fu_southpole_lat(pole2), &
                & 1) ! geographical to modified, preparation
    END IF
   
  END SUBROUTINE turn_latlon_to_latlon




  ! ***************************************************************
  
  SUBROUTINE turn_dxdy_to_latlon(dx, dy, grid_x, grid_y, grid, new_pole, dx_new, dy_new)

    ! Description:
    ! Turns a movement in a given grid to
    ! a movement in a lat-lon system defined by the given pole. 
    !
    ! The movement and the gridcoordinates of its starting point have
    ! to be along the axis of the given grid. Also the starting point
    ! from which the dx and dy are, has to be provided in order the
    ! turning to be possible. 
    !
    ! The unit movement of movement is free, so it can be either
    ! meters, gridsquares or grid-degrees. The unit is not checked
    ! nor converted.
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
    !
    ! Imported parameters with intent IN:
    REAL, INTENT(in) ::dx, dy ! movement along grid axis
    REAL, INTENT(in) :: grid_x, grid_y ! The gridcoordinates in
    ! which the dx and dy are
    TYPE(silja_grid), INTENT(in) :: grid ! of the previous
    TYPE(silam_pole), INTENT(in) :: new_pole
    
    ! Imported parameters with intent OUT:
    REAL, INTENT(out) :: dx_new, dy_new
    
    ! Local declarations:

    !----------------------------------------
    !
    ! 1. 
    !
    dx_new = real_missing
    dy_new = real_missing

    CALL set_error('sorry, not working yet','turn_dxdy_to_latlon')
    
  END SUBROUTINE turn_dxdy_to_latlon
  

  ! ***************************************************************

  REAL FUNCTION fu_horizontal_length(mov)
    
    ! Description:
    ! Returns the horizontal length of a vector, in metres.
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
    TYPE(silja_movement), INTENT(in) :: mov

    fu_horizontal_length = SQRT((mov%dx**2) + (mov%dy**2))

  END FUNCTION fu_horizontal_length



  ! ***************************************************************
  ! ***************************************************************
  !
  !
  !        HERE BEGINS THE VELOCITY-STUFF!!!
  !
  !
  ! ***************************************************************
  ! ***************************************************************
  
  FUNCTION fu_set_velocity_by_pole(u, v, w, pole, position)&
      & result(velo) 
    
    ! Sets a value for a velocity-vector. The given u-v-vomponents
    ! have to be in the lat-lon-system defined by the given pole.
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
    TYPE(silja_velocity) :: velo
    
    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: u ! in m/s, positive from west to east
    REAL, INTENT(in) :: v ! in m/s, positive from south to north
    REAL, INTENT(in) :: w ! in Pa/s
    TYPE(silam_pole), INTENT(in) :: pole
    TYPE(silam_grid_position), INTENT(in) :: position
    
    velo%mm = fu_set_movement(u, v, pole, position, w, one_second)
  
  END FUNCTION fu_set_velocity_by_pole


  ! ***************************************************************
  
  FUNCTION fu_set_velocity_by_grid (u, v, w, grid, position)&
      & result(velocity) 
    
    ! Description:
    ! Sets a value for a velocity-vector. The given u-v-vomponents
    ! have to be in the coordinate system defined by the grid.
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
    TYPE(silja_velocity) :: velocity
    
    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: u ! in m/s, positive from west to east
    REAL, INTENT(in) :: v ! in m/s, positive from south to north
    REAL, INTENT(in) :: w ! in Pa/s
    TYPE(silja_grid), INTENT(in) :: grid
    TYPE(silam_grid_position), INTENT(in) :: position

    IF (fu_position_inside_grid(position, grid) /= iPositionInside) THEN
      velocity = velocity_missing
      RETURN
    END IF

    SELECT CASE (fu_gridtype(grid))

    CASE (lonlat, anygrid)

      velocity = fu_set_velocity_by_pole(u, v, w, fu_pole(grid),&
	  & position)

    CASE default

      CALL set_error('sorry, can only handle latlon grids'&
                  & ,'fu_set_velocity_by_grid')

    END SELECT

  
  END FUNCTION fu_set_velocity_by_grid


  ! ***************************************************************

  LOGICAL FUNCTION fu_velocity_defined(velocity)
    
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_velocity), INTENT(in) :: velocity
    
    fu_velocity_defined = defined(velocity%mm)
    
  END FUNCTION fu_velocity_defined
  


  ! ***************************************************************
  
  FUNCTION fu_multiply_velocity1(x, velocity) result(new_velocity)
    
    ! Description:
    ! Returns a vector multiplied by factor x. All the components are
    ! handled.
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
    TYPE(silja_velocity) :: new_velocity
    !
    ! Imported parameters with intent(in):
    TYPE(silja_velocity), INTENT(in) :: velocity
    REAL, INTENT(in) :: x

    new_velocity%mm = x * velocity%mm

    new_velocity%mm%dt = one_second

    new_velocity%mm%defined = fu_set_true()

  END FUNCTION fu_multiply_velocity1


  ! ***************************************************************
  
  FUNCTION fu_multiply_velocity2(interval, velocity) result(movement)
    
    ! Description:
    ! Multiplies a velocity-vector with time-interval to create a
    ! movement.
    !
    IMPLICIT NONE

    ! The return value of this function:
    TYPE(silja_movement) :: movement
    
    ! Imported parameters with intent(in):
    TYPE(silja_velocity), INTENT(in) :: velocity
    TYPE(silja_interval), INTENT(in) :: interval

    movement = fu_sec(interval) * velocity%mm

    movement%defined = fu_set_true()

  END FUNCTION fu_multiply_velocity2


  ! ***************************************************************
  
  FUNCTION fu_add_velocities(velo1, velo2) result(velo_sum)
    
    ! Description:
    ! Adds two velocity-vectors.
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
    TYPE(silja_velocity) :: velo_sum
    
    ! Imported parameters with intent(in):
    TYPE(silja_velocity), INTENT(in) :: velo1, velo2

    velo_sum%mm = velo1%mm + velo2%mm
    velo_sum%mm%dt = one_second ! always
    velo_sum%mm%defined = fu_set_true()

  END FUNCTION fu_add_velocities


  ! ***************************************************************

  REAL FUNCTION fu_speed(velocity)
    
    ! Description:
    ! Returns the absolute, directionless speed-value of a velocity
    ! -vector. Result in m/s, of course. The vertical velocity is not
    ! taken care of at the moment.
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    
    ! Imported parameters with intent(in):
    TYPE(silja_velocity), INTENT(in) :: velocity
      
    fu_speed = fu_horizontal_length(velocity%mm)
  
  END FUNCTION fu_speed


  ! ***************************************************************

  REAL FUNCTION fu_direction_deg(velocity, pole)
    
    ! Description:
    ! Returns the direction of the velocity-vector in degrees.
    ! Direction is returned in lat-lon system defined by the pole.
    ! Directions:
    ! 0 deg. = from north to south
    ! 90 deg = from east to west
    ! 180 deg = from south to north
    ! 270 deg = from west to east
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    
    ! Imported parameters with intent(in):
    TYPE(silja_velocity), INTENT(in) :: velocity
    TYPE(silam_pole), INTENT(in) :: pole

    fu_direction_deg = 0.
    CALL set_error('sorry, not working','fu_direction_deg')

  END FUNCTION fu_direction_deg



  ! ***************************************************************

  SUBROUTINE scale_movement(movement, x_scale, y_scale, p_scale, t_scale)
  !
  ! Scales the given movement with the given coefficients.
  ! Used for transferring the movement to the relative units
  !
  ! Author: Mikhail Sofiev, FMI
  !
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN
    REAL, INTENT(in) :: x_scale, y_scale, p_scale, t_scale

    ! Imported parameters with intent INOUT
    TYPE(silja_movement), INTENT(inout) :: movement

    movement%dx = movement%dx * x_scale
    movement%dy = movement%dy * y_scale
    movement%dp = movement%dp * p_scale
    movement%dt = movement%dt * t_scale

  END SUBROUTINE scale_movement






  ! ***************************************************************
  ! ***************************************************************
  ! 
  !
  !                 TESTS
  !
  !
  ! ***************************************************************
  ! ***************************************************************
  
  SUBROUTINE movement_tests

    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Local declarations:
    TYPE(silam_grid_position) :: position
    TYPE(silja_movement) :: wind
    REAL :: dx, dy, dx2, dy2

    position = fu_set_position(60., 0., 85000., fu_wallclock())

    CALL report(position)

    dx = 10000.
    dy = 0.

    PRINT *, 'Movement in: ',  dx, dy

    wind = fu_set_movement(dx, dy, pole_geographical, & 
	& position, 20., one_hour)

    CALL print_movement(wind)
    
    position = position + wind

    CALL report(position)

!!!$    CALL turn_latlon_to_latlon & 
!!!$	& (dx, dy, position, pole_geographical, &
!!!$	& pole_system, dx2, dy2)
    
    !PRINT*, ' dx and dy out: ', dx2, dy2

  END SUBROUTINE movement_tests


  ! ***************************************************************
  
  SUBROUTINE print_movement(mov)
    
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    TYPE(silja_movement), INTENT(in) :: mov
    
    IF (defined(mov)) THEN
      PRINT *, ' dx and dy in system pole axis: ', mov%dx, mov%dy
      PRINT *, ' dp and dt: ', mov%dp, fu_sec(mov%dt)
    ELSE
      PRINT *, 'Undefined velocity. '
    END IF

  END SUBROUTINE print_movement



  ! ***************************************************************
  
  SUBROUTINE print_velocity(velo, position)
    
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    TYPE(silja_velocity), INTENT(in) :: velo
    TYPE(silam_grid_position), INTENT(in), OPTIONAL :: position

    ! Local :
    REAL :: dx, dy

    IF (defined(velo)) THEN
!!!$      PRINT *, ' x- and y-velocities in system pole axis: ',&
!!!$	  & velo%mm%dx, velo%mm%dy, ' [m/s]'
!!!$      
!!!$      PRINT *, ' w: ', velo%mm%dp, ' [Pa/s]'
    
      IF (PRESENT(position)) THEN
	
        CALL turn_latlon_to_latlon & 
             & (velo%mm%dx, velo%mm%dy, &
             & position,&
             & dispersion_pole, &
             & pole_geographical,&
             & dx, dy)
	
        IF (error) RETURN
	
        PRINT *, ' x- and y-velocities in geographical coords: ',&
                 & dx, dy, ' [m/s]'
      END IF
      
    ELSE
      PRINT*, 'Undefined velocity.'
    END IF
    
  END SUBROUTINE print_velocity
  
  
  !*****************************************************************************************
  
        SUBROUTINE TURNDD(PUARG,PVARG,PURES,PVRES,PA,PB,PC,PD, &
                     &  PXREG,PYREG,PXROT,PYROT, &
                     &  PXCEN,PYCEN,KCALL)
!C
!C     $Id: turndd.F,v 1.1 1997/01/15 12:48:04 keerola Exp keerola $
!C
!C
!C-----------------------------------------------------------------------
!C
!C*    TURN HORIZONTAL VELOCITY COMPONENTS BETWEEN REGULAR AND
!C*    ROTATED SPHERICAL COORDINATES.
!C
!C*    PUARG(KXDIM,KYDIM) : INPUT U COMPONENTS
!C*    PVARG(KXDIM,KYDIM) : INPUT V COMPONENTS
!C*    PURES(KXDIM,KYDIM) : OUTPUT U COMPONENTS
!C*    PVRES(KXDIM,KYDIM) : OUTPUT V COMPONENTS
!C*    PA(KXDIM,KYDIM)    : TRANSFORMATION COEFFICIENTS
!C*    PB(KXDIM,KYDIM)    :    -``-
!C*    PC(KXDIM,KYDIM)    :    -``-
!C*    PD(KXDIM,KYDIM)    :    -``-
!C*    PXREG(KXDIM,KYDIM) : REGULAR LONGITUDES
!C*    PYREG(KXDIM,KYDIM) : REGULAR LATITUDES
!C*    PXROT(KXDIM,KYDIM) : ROTATED LONGITUDES
!C*    PYROT(KXDIM,KYDIM) : ROTATED LATITUDES
!C*    KXDIM              : DIMENSION IN THE X (LONGITUDE) DIRECTION
!C*    KYDIM              : DIMENSION IN THE Y (LATITUDE) DIRECTION
!C*    KX                 : NUMBER OF GRIDPOINTS IN THE X DIRECTION
!C*    KY                 : NUMBER OF GRIDPOINTS IN THE Y DIRECTION
!C*    PXCEN              : REGULAR LONGITUDE OF THE SOUTH POLE OF THE
!C*                         TRANSFORMED GRID
!C*    PYCEN              : REGULAR LATITUDE OF THE SOUTH POLE OF THE
!C*                         TRANSFORMED GRID
!C*    KCALL=-2 OR 2      : PREPARATION OF COEFFICIENTS.
!C*    KCALL=-1 OR 1      : MULTIPLICATION AFTER PREPARATIONS.
!C*
!C*    KCALL < 0          : FIND WIND COMPONENTS IN REGULAR COORDINATES
!C*                         FROM WIND COMPONENTS IN ROTATED COORDINATES
!C*    KCALL > 0          : FIND WIND COMPONENTS IN ROTATED COORDINATES
!C*                         FROM WIND COMPONENTS IN REGULAR COORDINATES
!C*    NOTE THAT ALL COORDINATES ARE GIVEN IN DEGREES N AND DEGREES E.
!C*       (NEGATIVE VALUES FOR S AND W)
!C
!C*    J.E. HAUGEN   HIRLAM   JUNE -92, MODIFIED BY K. EEROLA
!C
!C!-----------------------------------------------------------------------
!C

      IMPLICIT NONE

      INTEGER, intent(in) :: KCALL
!      INTEGER, intent(in) :: KXDIM,KYDIM,KX,KY,KCALL
!      REAL, intent(in) :: PUARG(KXDIM,KYDIM),PVARG(KXDIM,KYDIM),&
!        &  PURES(KXDIM,KYDIM),PVRES(KXDIM,KYDIM),&
!        &     PA(KXDIM,KYDIM),   PB(KXDIM,KYDIM),&
!        &     PC(KXDIM,KYDIM),   PD(KXDIM,KYDIM),&
!        &  PXREG(KXDIM,KYDIM),PYREG(KXDIM,KYDIM),&
!        &  PXROT(KXDIM,KYDIM),PYROT(KXDIM,KYDIM)
      REAL, intent(in) :: PUARG,PVARG, PXREG,PYREG, PXROT,PYROT
      REAL, intent(out) :: PURES,PVRES, PA,   PB,  PC,   PD

      REAL, intent (in) :: PXCEN,PYCEN

!!PUARG,PVARG,PURES,PVRES,PA,PB,PC,PD, &
!                     &  PXREG,PYREG,PXROT,PYROT,KXDIM,KYDIM,KX,KY, &
!                     &  PXCEN,PYCEN,KCALL

!-----------------------------------------------------------------------

      REAL :: PI,ZRAD,ZSYC,ZCYC,ZSXREG,ZCXREG,ZSYREG,ZCYREG,ZXMXC, &
            & ZSXMXC,ZCXMXC,ZSXROT,ZCXROT,ZSYROT,ZCYROT, &
            & ZPXCEN,ZPYCEN

!-----------------------------------------------------------------------

      IF(PYCEN .GT. 0.0) THEN
         ZPYCEN = -PYCEN
         ZPXCEN =  0.0
      ELSE
         ZPYCEN = PYCEN
         ZPXCEN = PXCEN
      END IF

      IF (KCALL.EQ.1 .OR. KCALL.EQ.-1) THEN

!*    MULTIPLICATION BETWEEN REGULAR AND ROTATED SPHERICAL GRID

      PURES = PA*PUARG + PB*PVARG
      PVRES = PC*PUARG + PD*PVARG

      ELSEIF (KCALL.EQ.2) THEN

!*    PRECALCULATE MATRIXES FROM REGULAR TO ROTATED GRID

      PI = 4.*ATAN(1.)
      ZRAD = PI/180.
      ZSYC = SIN(ZRAD*(ZPYCEN+90.))
      ZCYC = COS(ZRAD*(ZPYCEN+90.))

      ZSXREG = SIN(ZRAD*PXREG)
      ZCXREG = COS(ZRAD*PXREG)
      ZSYREG = SIN(ZRAD*PYREG)
      ZCYREG = COS(ZRAD*PYREG)

      ZXMXC  = ZRAD*(PXREG - ZPXCEN)
      ZSXMXC = SIN(ZXMXC)
      ZCXMXC = COS(ZXMXC)

      ZSXROT = SIN(ZRAD*PXROT)
      ZCXROT = COS(ZRAD*PXROT)
      ZSYROT = SIN(ZRAD*PYROT)
      ZCYROT = COS(ZRAD*PYROT)

      PA = ZCYC*ZSXMXC*ZSXROT + ZCXMXC*ZCXROT
      PB = ZCYC*ZCXMXC*ZSYREG*ZSXROT - ZSYC*ZCYREG*ZSXROT - ZSXMXC*ZSYREG*ZCXROT
      PC = ZSYC*ZSXMXC/ZCYROT
      PD = (ZSYC*ZCXMXC*ZSYREG + ZCYC*ZCYREG)/ZCYROT

      ELSEIF (KCALL.EQ.-2) THEN

!*    PRECALCULATE MATRIXES FROM ROTATED TO REGULAR GRID

      PI = 4.*ATAN(1.)
      ZRAD = PI/180.
      ZSYC = SIN(ZRAD*(ZPYCEN+90.))
      ZCYC = COS(ZRAD*(ZPYCEN+90.))

      ZSXREG = SIN(ZRAD*PXREG)
      ZCXREG = COS(ZRAD*PXREG)
      ZSYREG = SIN(ZRAD*PYREG)
      ZCYREG = COS(ZRAD*PYREG)

      ZXMXC  = ZRAD*(PXREG - ZPXCEN)
      ZSXMXC = SIN(ZXMXC)
      ZCXMXC = COS(ZXMXC)

      ZSXROT = SIN(ZRAD*PXROT)
      ZCXROT = COS(ZRAD*PXROT)
      ZSYROT = SIN(ZRAD*PYROT)
      ZCYROT = COS(ZRAD*PYROT)

      PA = ZCXMXC*ZCXROT + ZCYC*ZSXMXC*ZSXROT
      PB = ZCYC*ZSXMXC*ZCXROT*ZSYROT + ZSYC*ZSXMXC*ZCYROT - ZCXMXC*ZSXROT*ZSYROT
      PC =-ZSYC*ZSXROT/ZCYREG
      PD = (ZCYC*ZCYROT - ZSYC*ZCXROT*ZSYROT)/ZCYREG

      ELSE
      WRITE(6,'(1X,''INVALID KCALL IN TURNWI'')')
      STOP
      ENDIF

      RETURN
      END SUBROUTINE TURNDD



!!!$  ! ***************************************************************
!!!$  
!!!$  SUBROUTINE print_movement(velo, position)
!!!$    
!!!$    ! Language: ANSI Fortran 90
!!!$    !
!!!$    ! Author: Mika Salonoja, FMI
!!!$    ! 
!!!$    IMPLICIT NONE
!!!$    !
!!!$    ! Imported parameters with intent IN:
!!!$    TYPE(silja_movement), INTENT(in) :: mov
!!!$    TYPE(silja_position), INTENT(in), OPTIONAL :: position
!!!$
!!!$    ! Local :
!!!$    REAL :: dx, dy
!!!$
!!!$    IF (defined(mov)) THEN
!!!$      PRINT *, ' x- and y-movements in system pole axis: ',&
!!!$	  & mov%dx, mov%dy, ' [m]'
!!!$      
!!!$      PRINT *, ' w: ', mov%dp, ' [Pa/s]'
!!!$    
!!!$      IF (PRESENT(position)) THEN
!!!$	
!!!$	CALL turn_latlon_to_latlon & 
!!!$	    & (mov%dx, mov%dy, &
!!!$	    & position,&
!!!$	    & pole_system, &
!!!$	    & pole_geographical,&
!!!$	    & dx, dy)
!!!$	
!!!$	IF (error) RETURN
!!!$	
!!!$	PRINT *, ' x- and y-movements in geographical coords: ',&
!!!$	    & dx, dy, ' [m]'
!!!$	
!!!$      END IF
!!!$      
!!!$    ELSE
!!!$      PRINT*, 'Undefined movement'
!!!$      
!!!$      
!!!$    END IF
!!!$    
!!!$  END SUBROUTINE print_movement
!!!$  
END MODULE vectors_v2
