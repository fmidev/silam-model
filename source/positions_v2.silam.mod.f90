MODULE positions_v2 ! the position in x,y,z and time

  ! Description:
  ! This module contains the definition of the type silam_grid_position,
  ! which defines a (x,y,z,t)-location.
  !
  ! The main idea for this module (different from the old one):
  ! All motions in the model must be in dispersion_grid, which for Lagrangian 
  ! dynamics is defined by meteorology alone. The grid_position type
  ! makes this work - it contains x and y, not lon and lat co-ordinates
  ! and does no transformation in case if someone just wants to get 
  ! the co-ordinates. Purpose: keep the interpolation of meteodata during advection
  ! computations to a minimum level.
  !
  ! Routines to move from silam_grid_position to "normal" output grid
  ! are taken by appropriate routines below, where EASTERN LONGITUDES ARE 
  ! POSITIVE; NORTHERN LATITUDES ARE POSITIVE, unless explicit rule stated.
  !
  ! The vertical location can be also relative, or pressure (Pa), altitude
  ! (m above mean sea level) or height (m above ground level).  
  ! This module provides tools for returning the value position's z
  ! -coordinate in either relative height, altitude or pressure. Non-relative
  ! units are kept only temporarily and should be removed when vertical
  ! velocity tools cover all necessary stuff.
  !
  ! Author: Mikhail Sofiev, FMI
  ! Large portions of code are taken from original positions of Miko Salonoja
  ! 
  ! All units: SI
  ! 
  ! All degrees: real numbers in degrees; NOT in degrees and minutes
  ! 
  ! Language: ANSI Fortran 90
  ! 
  ! ***************************************************************
  !
  ! Modules used:

  USE silam_times  !optimisation
  use silam_levels
  USE grids_geo

  IMPLICIT NONE

  ! The public functions and subroutines available in this module:
  PUBLIC fu_set_position  ! With grid transformation, if needed
  PUBLIC set_position  ! Same, but subroutine
  public fu_set_pos_in_geo_global_grid
  PUBLIC update_height
  PUBLIC update_z
  PUBLIC update_time
  PUBLIC update_pressure
  PUBLIC fu_compare_positions_eq ! interface ==
  PUBLIC defined
  PUBLIC fu_pressure
  PUBLIC fu_height
  PUBLIC fu_height_up_to_date
  PUBLIC fu_geographical_latlon_string
!  PUBLIC fu_geographical_lat
!  PUBLIC fu_geographical_lon
  PUBLIC fu_time
  PUBLIC pos_tests
  PUBLIC fu_horizontal_distance ! Do not use ! SLOW.
  PUBLIC fu_pressure_difference
  PUBLIC fu_position_new_pressure
  PUBLIC position_elevator
  PUBLIC make_no_tunnels
  PUBLIC fu_2d_interpolation !Taken from grids. Fast.
!  PUBLIC fu_silja_4d_interpolation ! New. Uses 2d-routine. Fast
  PUBLIC fu_x  ! Relative grid units
  PUBLIC fu_y
  PUBLIC fu_z
  PUBLIC get_lonlat
  PUBLIC position_to_gridcoordinates
  PUBLIC fu_position_inside_grid
  PUBLIC fu_position_gridindex
  public fu_position_vert_index
  public position_to_new_grid_index  ! with a new grid


  PUBLIC report

  ! Private functions:
  PRIVATE fu_time_of_position
  PRIVATE fu_pr_diff_of_positions
  PRIVATE fu_position_defined
  PRIVATE fu_pressure_of_position
  PRIVATE fu_move_positions_time
  PRIVATE print_position
  PRIVATE fu_height_of_position
  PRIVATE fu_position_height_up_to_date
  private fu_2d_interp_to_position !Taken from grids. Fast.
  private fu_position_vert_index_levels
  private fu_position_vert_index_vert

  ! Generic names and operator-interfaces of some functions:

  INTERFACE defined
    MODULE PROCEDURE fu_position_defined
  END INTERFACE
  
  INTERFACE operator(==)
    MODULE PROCEDURE fu_compare_positions_eq
  END INTERFACE
  
  INTERFACE fu_time
    MODULE PROCEDURE fu_time_of_position
  END INTERFACE
  
!!!$  INTERFACE operator(+)
!!!$    MODULE PROCEDURE fu_move_positions_time
!!!$  END INTERFACE
  
  INTERFACE fu_pressure_difference
    MODULE PROCEDURE fu_pr_diff_of_positions
  END INTERFACE
    
  INTERFACE fu_pressure
    MODULE PROCEDURE fu_pressure_of_position
  END INTERFACE

  INTERFACE fu_height
    MODULE PROCEDURE fu_height_of_position
  END INTERFACE

  INTERFACE fu_height_up_to_date
    MODULE PROCEDURE fu_position_height_up_to_date
  END INTERFACE

  interface fu_2d_interpolation
    module procedure fu_2d_interp_to_position
  end interface

  INTERFACE report
    MODULE PROCEDURE print_position
  END INTERFACE

  interface fu_position_vert_index
    module procedure fu_position_vert_index_levels
    module procedure fu_position_vert_index_vert
  end interface


  ! Public types with private components defined in this module:
  TYPE silam_grid_position  
    PRIVATE
    REAL ::  x, y, z
    REAL :: pressure ! [Pa]
    REAL :: height ! [m] above ground
    LOGICAL :: height_up_to_date
    TYPE(silja_time) :: time
    TYPE(silja_logical) :: defined
  END TYPE silam_grid_position
  
  ! The missing codes for these types:
  TYPE(silam_grid_position), PARAMETER :: position_missing = &
      & silam_grid_position(real_missing,&
                          & real_missing,&
                          & real_missing,&
                          & real_missing,&
                          & real_missing,&
                          & .false.,&
                          & time_missing,&
                          & silja_false)

  integer, public, parameter :: iPositionInside = 4000 
  integer, public, parameter :: iPositionOutside_x0 = 4001  ! ix < 1
  integer, public, parameter :: iPositionOutside_y0 = 4002  ! iy < 1
  integer, public, parameter :: iPositionOutside_xM = 4003  ! ix > xMax
  integer, public, parameter :: iPositionOutside_yM = 4004  ! iy > yMax
  integer, public, parameter :: iPositionOutside_vert_0 = 4005  ! iLev < 0
  integer, public, parameter :: iPositionOutside_vert_M = 4006  ! iLev > nLev+1
  integer, public, parameter :: iPositionOutside_unknown = 4007


CONTAINS

  ! ***************************************************************
  ! ****************************************************************
  ! 
  !
  !    HERE BEGINS SETTING AND COMPARING OF POSITIONS
  !
  !
  ! ***************************************************************
  ! ***************************************************************
  
  FUNCTION fu_set_position(x,y,pressure,time,input_grid) RESULT(pos)
    !
    ! Sets a value for one position-variable, checks the legality
    ! of its components and converts all the components according to
    ! the internal position format of silam-program. 
    ! input_grid defines the way of treatment of other variables.
    !
    ! If input grid is defined, the x and y is, of course, co-ordinates
    ! in this grid.
    !
    ! NOTE. If present, the system_grid must be defined.
    ! TRICK. If no grid transformation is desireable - just skip input_grid
    ! 
    ! Author: Mikhail Sofiev, FMI
    
    IMPLICIT NONE
    
    ! The return value of this function
    TYPE(silam_grid_position) :: pos

    ! Imported parameters with intent IN:
    REAL, INTENT(in) :: x, y, pressure ! in Pa, NOT HPa!!!!
    TYPE(silja_time), INTENT(in) :: time
    TYPE(silja_grid), INTENT(in), OPTIONAL :: input_grid
    
    ! Local variables
    REAL :: meteo_corner_lon_E,meteo_corner_lat_N, input_corner_lon_E,input_corner_lat_N, &
          & southpole_lon_E, southpole_lat_N, & 
          & meteo_dx_deg, meteo_dy_deg, input_dx_deg, input_dy_deg
    LOGICAL :: corner_in_geographical_latlon, if_south_pole
    INTEGER :: number_of_points_x, number_of_points_y

    ! In theory - below each-to-each grid conversion should be located.
    ! So far just a part of the stuff is ready.

    IF(PRESENT(input_grid))THEN

!      IF(.not.defined(dispersion_grid))THEN
      IF(.not.defined(meteo_grid))THEN
        CALL set_error('Meteo grid is still undefined','set_position')
        RETURN
      END IF

      call project_point_to_grid(input_grid, x, y, meteo_grid, pos%x, pos%y)


    ELSE  ! No input_grid => simply accept the values given. NO CHECKING
      
      pos%x = x
      pos%y = y

    ENDIF

    !-------------------- Set time and pressure.
    !
    pos%pressure = pressure
    pos%time = time
    pos%height = real_missing
    pos%height_up_to_date = .false.

    pos%defined = fu_set_true()

  END FUNCTION fu_set_position


  !************************************************************

  subroutine set_position(x,y,pressure,time,input_grid, pos)
    !
    ! Sets a value for one position-variable, checks the legality
    ! of its components and converts all the components according to
    ! the internal position format of silam-program. 
    ! input_grid defines the way of treatment of other variables.
    !
    ! If input grid is defined, the x and y is, of course, co-ordinates
    ! in this grid.
    !
    ! NOTE. If present, the meteo_grid must be defined.
    ! TRICK. If no grid transformation is desireable - just skip input_grid
    ! 
    ! Author: Mikhail Sofiev, FMI
    
    IMPLICIT NONE
    
    ! The return value 
    TYPE(silam_grid_position), intent(out) :: pos

    ! Imported parameters with intent IN:
    REAL, INTENT(in) :: x, y, pressure ! in Pa, NOT HPa!!!!
    TYPE(silja_time), INTENT(in) :: time
    TYPE(silja_grid), INTENT(in), OPTIONAL :: input_grid
    
    ! Local variables
    REAL :: meteo_corner_lon_E, meteo_corner_lat_N, input_corner_lon_E,input_corner_lat_N, &
          & southpole_lon_E, southpole_lat_N, & 
          & meteo_dx_deg, meteo_dy_deg, input_dx_deg, input_dy_deg
    LOGICAL :: corner_in_geographical_latlon, if_south_pole
    INTEGER :: number_of_points_x, number_of_points_y

    ! In theory - below each-to-each grid conversion should be located.
    ! So far just a part of the stuff is ready.

    IF(PRESENT(input_grid))THEN

!      IF(.not.defined(dispersion_grid))THEN
!        CALL set_error('dispersion_grid is still undefined','set_position')
      IF(.not.defined(meteo_grid))THEN
        CALL set_error('meteo_grid is still undefined','set_position')
        RETURN
      END IF

      call project_point_to_grid(input_grid, x, y, meteo_grid, pos%x, pos%y)


    ELSE  ! No input_grid => simply accept the values given. NO CHECKING
      
      pos%x = x
      pos%y = y

    ENDIF

    !-------------------- Set time and pressure.
    !
    pos%pressure = pressure
    pos%time = time
    pos%height = real_missing
    pos%height_up_to_date = .false.

    pos%defined = fu_set_true()

  END subroutine set_position


  !****************************************************************

  function fu_set_pos_in_geo_global_grid(x,y,pressure,time) result(pos)
    !
    ! Sets a value for one position-variable, relative to the geo_global grid
    !
    ! Co-ordinates are supposed to be in normal geographical co-ordinates, north and east
    ! positive. all in degrees and decimals.
    !
    ! Author: Mikhail Sofiev, FMI
    
    IMPLICIT NONE
    
    ! The return value of this function
    TYPE(silam_grid_position) :: pos

    ! Imported parameters with intent IN:
    REAL, INTENT(in) :: x, y, pressure ! in Pa, NOT HPa!!!!
    TYPE(silja_time), INTENT(in) :: time
    

    !
    ! Call geo_global grid parameters and just scale the (x,y) position to it

    call project_point_to_grid(x, y, geo_global_grid, pos%x, pos%y)

    !-------------------- Set time and pressure.
    !
    pos%pressure = pressure
    pos%time = time
    pos%height = real_missing
    pos%height_up_to_date = .false.

    pos%defined = fu_set_true()

  end function fu_set_pos_in_geo_global_grid

  ! ***************************************************************

  SUBROUTINE update_height(position, height)
    
    ! Description:
    ! Adds metric height to a position. The value cannot be checked
    ! here, so it is up to user to given correct values.
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    REAL, INTENT(in) :: height ! above ground [m]

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silam_grid_position), INTENT(inout) :: position

    IF (.NOT.defined(position)) THEN
      CALL set_error('undefined position given', 'update_height')
    ELSE
      position%height = height
      position%height_up_to_date = .true.
    END IF

  END SUBROUTINE update_height


  ! ***************************************************************

  SUBROUTINE update_z(position, z) ! Relative grid units
    
    ! Adds relative height to a position. The value cannot be checked
    ! here, so it is up to user to given correct values. No checking to 
    ! speed up the process - so, do not use wrong or undefined position !
    !
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE

    REAL, INTENT(in) :: z ! relative units depending on grid

    TYPE(silam_grid_position), INTENT(inout) :: position

    position%z = z

  END SUBROUTINE update_z


  ! ***************************************************************

  SUBROUTINE update_time(position, time)
  !
  ! Sets new time of position
  !
  ! Author: Mikhail Sofiev, FMI
  ! 
    IMPLICIT NONE
 
    TYPE(silja_time), INTENT(in) :: time

    TYPE(silam_grid_position), INTENT(inout) :: position

    position%time = time

  END SUBROUTINE update_time



  ! ***************************************************************

  SUBROUTINE update_pressure(position, pressure)
  !
  ! Sets new time of position
  !
  ! Author: Mikhail Sofiev, FMI
  ! 
    IMPLICIT NONE
 
    REAL, INTENT(in) :: pressure

    TYPE(silam_grid_position), INTENT(inout) :: position

    position%pressure = pressure
    position%height_up_to_date = .false.

  END SUBROUTINE update_pressure


  ! ***************************************************************

  LOGICAL FUNCTION fu_compare_positions_eq(pos1, pos2)

    ! Compares two positions and returns a true value if equal. This
    ! function is used through an interface (pos1 == pos2)
    ! 
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    TYPE(silam_grid_position), INTENT(in) :: pos1, pos2

    fu_compare_positions_eq =(&
         &(pos1%x.eps.pos2%x) .and. &
         &(pos1%y.eps.pos2%y) .and. &
         &(pos1%pressure.eps.pos2%pressure) .and. &
         &(pos1%time == pos2%time))
    
  END FUNCTION fu_compare_positions_eq


  ! ***************************************************************

  SUBROUTINE get_lonlat(pos, pole, lon, lat)
    ! 
    ! Returns the coordinates of a position in lat-lon-system defined
    ! by the pole. 
    !
    ! Northern latitudes (N) and eastern longitudes (E) are positive.
    !
    ! Language: ANSI Fortran 90
    !
    ! Original code: Mika Salonoja
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    TYPE(silam_grid_position), INTENT(in) :: pos
    TYPE(silam_pole), INTENT(in) :: pole

    ! Imported parameters with intent OUT:
    REAL, INTENT(out) :: lat, lon

    ! Convert the coordinates:

    CALL modify_lonlat(fu_lat_native_from_grid(pos%x, pos%y, meteo_grid), &
                     & fu_lon_native_from_grid(pos%x, pos%y, meteo_grid), &
                     & meteo_pole, pole, lat, lon)
!    CALL modify_lonlat(fu_lat_native_from_grid(pos%x, pos%y, dispersion_grid), &
!                     & fu_lon_native_from_grid(pos%x, pos%y, dispersion_grid), &
!                     & dispersion_pole, pole, lat, lon)

  END SUBROUTINE get_lonlat
  
 

  ! ***************************************************************

  SUBROUTINE position_to_gridcoordinates(pos, grid, x, y, out_of_grid)

    ! Returns the given position as grid coordinates of given grid. Actually
    ! it is transformation of system_grid to another grid, because in
    ! both cases relative grid co-ordinates are used.
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    TYPE(silam_grid_position), INTENT(in) :: pos
    TYPE(silja_grid), INTENT(in) :: grid

    ! Imported parameters with intent OUT:
    REAL, INTENT(out) :: x,y ! the gridcoordinates of the position
    LOGICAL, INTENT(out) :: out_of_grid ! true, when position is
    ! outside the grid
    !
    ! Local declarations:
    REAL :: lat_grid, lon_grid, &
            & corner_lat_N, corner_lon_E, &
            & southpole_lat_N, southpole_lon_E, & 
            & dx_deg, dy_deg
    LOGICAL :: corner_in_geographical_latlon
    INTEGER :: nx, ny

!    IF(grid == dispersion_grid)THEN

    call grid_dimensions(grid, nx, ny)
    IF(grid == meteo_grid)THEN
      x = fu_x(pos)
      y = fu_y(pos)
      out_of_grid = x<0.5 .or. y<0.5 .or. x>nx+0.5 .or. y>ny+0.5
      RETURN
    END IF


    CALL get_lonlat(pos, fu_pole(grid), lon_grid, lat_grid)
    call project_point_to_grid(lon_grid, lat_grid, grid, x, y)
   
    IF ((x >= 1.) .and. (x <= REAL(nx)) .and. &
      & (y >= 1.) .and. (y <= REAL(ny))) THEN
      out_of_grid = .false.
    ELSE
      out_of_grid = .true.
    END IF
  

  END SUBROUTINE position_to_gridcoordinates


  ! ***************************************************************

  integer FUNCTION fu_position_inside_grid(position, grid, x_lag_, y_lag_)

    ! Returns true value if position is horizontally inside grid
    ! boundaries. Adapted to system_grid: position is related to system_grid.
    ! Ineffecient if grid /= system_grid.
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent(in):
    TYPE(silja_grid), INTENT(in) :: grid
    TYPE(silam_grid_position), INTENT(in) :: position
    REAL, INTENT(in), OPTIONAL :: x_lag_, y_lag_

    ! Local declarations:
    REAL :: x, y, x_lag, y_lag ! If Lags>0 => reserve band around grid border
    INTEGER :: nx,ny
    LOGICAL :: out

    IF(PRESENT(x_lag_)) THEN
      x_lag = MAX(0.,x_lag_)
    ELSE
      x_lag =0.
    END IF
    IF(PRESENT(y_lag_)) THEN
      y_lag = MAX(0.,y_lag_)
    ELSE
      y_lag =0.
    END IF
    fu_position_inside_grid = iPositionInside
    !
    ! All positions are represented in the system_grid relative indices. 
    ! Therefore, if the given grid is non-system, we have to apply transformation
    !
    CALL grid_dimensions(grid,nx,ny)

    IF(grid == meteo_grid)THEN
!    IF(grid == dispersion_grid)THEN
      if(ifPrintDump) call msg('from the position check: dispersion_gridgiven')
      x=fu_x(position)
      y=fu_y(position)
    ELSE
      if(ifPrintDump) call msg('from the position check: non-meteo_grid given')
      CALL position_to_gridcoordinates(position, grid, x, y, out)
      IF(out.or.error)THEN
        fu_position_inside_grid = iPositionOutside_unknown
        return
      endif
    END IF

    if(ifPrintDump)then
        call report(position)
        call msg('nx,x=',nx,x)
        call msg('ny,y=',ny,y)
    endif

    !
    ! Now the x, y coordinates are in the given grid and out-of-grid checking is simple
    !
    if(x < 1.+x_lag)then
      if(ifPrintDump) call msg('iPositionOutside_x0')
      fu_position_inside_grid = iPositionOutside_x0
      return
    endif
    if(y < 1.+y_lag)then
      if(ifPrintDump) call msg('iPositionOutside_y0')
      fu_position_inside_grid = iPositionOutside_y0
      return
    endif
    if(x > nx-x_lag)then
      if(ifPrintDump) call msg('iPositionOutside_xM')
      fu_position_inside_grid = iPositionOutside_xM
      return
    endif
    if(y > ny-y_lag)then
      if(ifPrintDump) call msg('iPositionOutside_yM')
      fu_position_inside_grid = iPositionOutside_yM
      return
    endif

  END FUNCTION fu_position_inside_grid


  ! ***************************************************************

  INTEGER FUNCTION fu_position_gridindex(position, grid)
    
    ! Returns the index of the gridsquare the position is in. If
    ! position is outside the grid, zero is returned.
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
    TYPE(silja_grid), INTENT(in) :: grid

    ! Local declarations:
    REAL :: x,y
    LOGICAL :: out_of_grid
    INTEGER :: nx,ny

    fu_position_gridindex = 0
    CALL position_to_gridcoordinates(position, grid, x, y, out_of_grid)
    IF (error.or.out_of_grid) RETURN

    CALL grid_dimensions(grid, nx,ny)

    fu_position_gridindex = ((INT(y+0.5)-1)*nx) + INT(x+0.5)
  
  END FUNCTION fu_position_gridindex


  !*******************************************************************

  real function fu_position_vert_index_levels(position, levels)
    !
    ! Finds relative index of the position in the vertical structure
    ! It does not compute anything - just sets the temporary level
    ! and calls the corresponding subroutine for comparison.
    ! A nice feature is that levels and thick layers will result in the 
    ! same index. So, we do not have to care about it
    !
    implicit none

    ! Imported parameters
    type(silam_grid_position), intent(in) :: position
    type(silja_level), dimension(:), intent(in) :: levels

    ! Local variables
    type(silja_level) :: levTmp

    fu_position_vert_index_levels = -1.
    select case(fu_leveltype(fu_central_level_of_layer(levels(1))))
      case(constant_pressure)
        levTmp = fu_set_pressure_level(position%pressure)
      case(constant_height)
        if(.not.position%height_up_to_date)then
          call report(position)
          call set_error('Height is not up to date','fu_position_vert_index_levels')
          return
        endif
        levTmp = fu_set_constant_height_level(position%height)
      case default
        call report(levels(1))
        call set_error('Unknown level type','fu_position_vert_index_levels')
        return
    end select

    fu_position_vert_index_levels = fu_level_index(levTmp,levels)

  end function fu_position_vert_index_levels

  !*******************************************************************

  real function fu_position_vert_index_vert(position, vertical)
    !
    ! Finds relative index of the position in the vertical structure
    ! It does not compute anything - just sets the temporary level
    ! and calls the corresponding subroutine for comparison.
    ! A nice feature is that levels and thick layers will result in the 
    ! same index. So, we do not have to care about it
    !
    implicit none

    ! Imported parameters
    type(silam_grid_position), intent(in) :: position
    type(silam_vertical), intent(in) :: vertical

    ! Local variables
    type(silja_level) :: levTmp

    fu_position_vert_index_vert = -1.
    select case(fu_leveltype(fu_central_level_of_layer(fu_level(vertical,1))))
      case(constant_pressure)
        levTmp = fu_set_pressure_level(position%pressure)
      case(constant_height)
        if(.not.position%height_up_to_date)then
          call report(position)
          call set_error('Height is not up to date','fu_position_vert_index_vert')
          return
        endif
        levTmp = fu_set_constant_height_level(position%height)
      case default
        call report(vertical)
        call set_error('Unknown level type','fu_position_vert_index_vert')
        return
    end select

    fu_position_vert_index_vert = fu_level_index(levTmp,vertical)

  end function fu_position_vert_index_vert


  ! ***************************************************************

  LOGICAL FUNCTION fu_position_defined(pos)
 
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silam_grid_position), INTENT(in) :: pos

    fu_position_defined = fu_true(pos%defined)
    
  END FUNCTION fu_position_defined



  ! ***************************************************************

  SUBROUTINE position_elevator(position, pressure_difference)

    ! Description:
    ! Moves a defined position in vertical the given amount.
    ! Position difference value = moving downwards (increasing
    ! pressure).
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    REAL, INTENT(in) :: pressure_difference

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silam_grid_position), INTENT(inout) :: position

    position%pressure = position%pressure + pressure_difference
    position%height_up_to_date = .false.


  END SUBROUTINE position_elevator


  ! ***************************************************************

  FUNCTION fu_position_new_pressure(position, new_pressure)
    
    ! Returns a position with new value for pressure. No height
    ! is updated.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silam_grid_position) :: fu_position_new_pressure
    !
    ! Imported parameters with intent(in):
    TYPE(silam_grid_position), INTENT(in) :: position
    REAL, INTENT(in) :: new_pressure

    fu_position_new_pressure = position
    fu_position_new_pressure%pressure = new_pressure
    fu_position_new_pressure%height_up_to_date = .false.
  
  END FUNCTION fu_position_new_pressure



  ! ***************************************************************

  SUBROUTINE make_no_tunnels(position, surface_pressure)
    
    ! Description:
    ! Lifts the position to ground level in case surface pressure is
    ! lower 
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    REAL, INTENT(in) :: surface_pressure

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silam_grid_position), INTENT(inout) :: position

    IF (defined(position)) THEN
      position%pressure = MIN(position%pressure, surface_pressure)
      position%height_up_to_date = .false.
    END IF
    
  END SUBROUTINE make_no_tunnels



  ! ***************************************************************

  FUNCTION fu_move_positions_time(position, interval) result(new_position)
  
    ! Description:
    ! Changes the position's time for the interval.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silam_grid_position) :: new_position
    !
    ! Imported parameters with intent(in):
    TYPE(silam_grid_position), INTENT(in) :: position
    TYPE(silja_interval), INTENT(in) :: interval

    new_position = position
    new_position%height_up_to_date = .false.
    new_position%time = position%time + interval
  
  END FUNCTION fu_move_positions_time



  ! ****************************************************************
  ! 
  !
  !          HORIZONTAL INTERPOLATION
  !
  !
  ! ***************************************************************


  integer function position_to_new_grid_index(position, &     ! position itself
                                            & grid_new, nx_new, ny_new) ! new grid and its dimensions
    !
    ! Position is always in meteo_grid relative co-ordinates. This function
    ! projects it to new grid. Made as fast as possible, thus the grid dimensions 
    ! are not called but rather imported to save a call of a sub.
    !
    implicit none

    ! Imported parameters
    TYPE(silja_grid), INTENT(in) :: grid_new
    TYPE(silam_grid_position), INTENT(in) :: position
    integer, intent(in) :: nx_new, ny_new

    ! Local variables
    real :: xNew, yNew

    call project_point_to_grid(meteo_grid, position%x, position%y, &
                             & grid_new, xNew, yNew)
    if(xNew < 0.5 .or. xNew > nx_new+0.5 .or. yNew < 0 .or. yNew > ny_new+0.5)then
      position_to_new_grid_index = int_missing
    else
      position_to_new_grid_index = int(xNew+0.5) + (int(yNew+0.5)-1)*nx_new
    endif

  end function position_to_new_grid_index


  !****************************************************************

  REAL FUNCTION fu_2d_interp_to_position (grid_data, grid, position, method) result(value)
    !
    ! Interpolates the field given in a grid to given position. This code
    ! version is adapted to meteo_grid conditions.
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Original code: Pilvi Siljamo, FMI
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parametersn with intent(in):
    TYPE(silja_grid), INTENT(in) :: grid
    TYPE(silam_grid_position), INTENT(in) :: position
    INTEGER, INTENT(in) :: method ! of interpolation

    ! Imported parametersn with intent(inout) or POINTER:
    REAL, DIMENSION(:), POINTER :: grid_data

    ! Local declarations:
    INTEGER :: southwest_corner, loop, xloop, yloop, ii
    INTEGER, DIMENSION(16)::corner_pointers 
    LOGICAL :: out_of_grid
    REAL :: x, y, zx, zy 
    REAL, DIMENSION(16)::ww
    INTEGER:: nx, ny
    REAL:: zxy1, zxy2, zxy3, zxy4, zmin
    REAL:: zx1,zx2,zx3,zx4,zy1,zy2,zy3,zy4     
    INTEGER:: ilat, ilon
    INTEGER :: method_local
    REAL, DIMENSION(4) :: distances
    INTEGER, DIMENSION(1) :: point

    ! -----------------------------------------------
    ! 
    ! 1. Get the gridpoint coordinates of position. Fast if grid==system_grid
    !    ------------------------------------------

    CALL position_to_gridcoordinates(position, grid, x, y, out_of_grid)
    IF (error) RETURN

    IF(out_of_grid)THEN
!      CALL report(position)
!      CALL report(grid)
      CALL set_error_nomsg('position out of grid!','fu_silja_2d_int (pos_to_grid)')
      RETURN
    END IF

    CALL grid_dimensions(grid, nx,ny)

    IF (method == cubic) THEN
      ! If cubic method and position too close to the border of the
      ! grid, use linear instead:
      IF ((x<2.).or.(y<2.).or.(x>(nx-1.)).or.(y>(ny-1.))) THEN
        method_local = linear
      ELSE
        method_local = cubic
      END IF
    ELSE
      method_local = method
    END IF

    method_of_interpolation: SELECT CASE (method_local)

      ! -----------------------------------------------
      ! 
      ! 2. Linear interpolation.
      !    --------------------

      CASE (linear)

      corner_pointers (1) = 0 
      corner_pointers (2) = 1 
      corner_pointers (3) = nx 
      corner_pointers (4) = nx + 1

      southwest_corner = ((INT(y) - 1)*nx) + INT(x)

      zx = x - REAL(INT(x))
      zy = y - REAL(INT(y))

      ww(1) = (1. - zx) * (1. - zy)
      ww(2) = zx * (1. - zy)
      ww(3) = (1. - zx) * zy
      ww(4) = zx * zy

      value = 0.

      DO loop = 1, 4
        value = value + ww(loop) * grid_data(southwest_corner + corner_pointers(loop))
      END DO

      ! -----------------------------------------------
      ! 
      ! 3. Bi-cubic interpolation.
      !    --------------------

      CASE (cubic)

      ii = 0
      ilat = -nx !      

      DO yloop = 1, 4 
        ilon = -1 !
        DO xloop = 1, 4 
          ii = ii + 1
          corner_pointers(ii) = ilon + ilat
          ilon = ilon + 1
        END DO
        ilat = ilat + nx
      END DO

      southwest_corner = ((INT(y) - 1)*nx) + INT(x)

      zx = x - REAL(INT(x))
      zy = y - REAL(INT(y))

      zx1 = ((-0.5*zx+1.0)*zx-0.5)*zx 
      zx2 = (( 1.5*zx-2.5)*zx    )*zx+1 
      zx3 = ((-1.5*zx+2.0)*zx+0.5)*zx 
      zx4 = (( 0.5*zx-0.5)*zx    )*zx 
      zy1 = ((-0.5*zy+1.0)*zy-0.5)*zy 
      zy2 = (( 1.5*zy-2.5)*zy    )*zy+1 
      zy3 = ((-1.5*zy+2.0)*zy+0.5)*zy 
      zy4 = (( 0.5*zy-0.5)*zy    )*zy

      ww( 1) = zx1*zy1 
      ww( 2) = zx2*zy1
      ww( 3) = zx3*zy1  
      ww( 4) = zx4*zy1 
      ww( 5) = zx1*zy2  
      ww( 6) = zx2*zy2  
      ww( 7) = zx3*zy2  
      ww( 8) = zx4*zy2  
      ww( 9) = zx1*zy3  
      ww(10) = zx2*zy3  
      ww(11) = zx3*zy3 
      ww(12) = zx4*zy3 
      ww(13) = zx1*zy4 
      ww(14) = zx2*zy4   
      ww(15) = zx3*zy4   
      ww(16) = zx4*zy4   

      value = 0.

      DO loop = 1, 16
        value = value + ww(loop) * grid_data(southwest_corner + corner_pointers(loop))
      END DO

      ! -----------------------------------------------
      ! 
      ! 4. Take the nearest point
      !    --------------------

      CASE (nearest_point)

      corner_pointers (1) = 0 
      corner_pointers (2) = 1 
      corner_pointers (3) = nx 
      corner_pointers (4) = nx + 1

      southwest_corner = ((INT(y) - 1)*nx) + INT(x)

      zx = x - REAL(INT(x))
      zy = y - REAL(INT(y))

      distances(1) = zx*zx + zy*zy
      distances(2) = (1.-zx)*(1.-zx) + zy*zy
      distances(3) = zx*zx + (1.-zy)*(1.-zy)
      distances(4) = (1.-zx)*(1.-zx) + (1.-zy)*(1.-zy)

      point = MINLOC(distances)

      value = grid_data(southwest_corner + corner_pointers(point(1)))




    CASE default

      CALL set_error('unknown interpolation method','fu_2d_interpolation')


    END SELECT method_of_interpolation

  END FUNCTION fu_2d_interp_to_position


  ! ***************************************************************
  ! ***************************************************************
  !
  !
  !        HERE BEGINS RETURNING THE VALUES OF POSITIONS
  !
  !
  ! ***************************************************************
  ! ***************************************************************
  
  REAL FUNCTION fu_pressure_of_position(position)
    
    ! Description:
    ! Returns the value of the position's z-coordinate in pressure
    ! (Pa, of course!).  
    !
    ! Author: Mika Salonoja, FMI
    ! 
    ! Imported parameters with intent IN:
    TYPE(silam_grid_position), INTENT(in) :: position

    IF (defined(position)) THEN
      fu_pressure_of_position = position%pressure
    ELSE
      CALL set_error('undefined position given'&
	  & ,'fu_pressure_of_position')
    END IF


  END FUNCTION fu_pressure_of_position
  

  ! ***************************************************************

  REAL FUNCTION fu_height_of_position(position)
    
    ! Description:
    ! Returns position's height above ground [m].
    !
    ! Author: Mika Salonoja, FMI
    ! 
    ! Imported parameters with intent IN:
    TYPE(silam_grid_position), INTENT(in) :: position

    IF (defined(position)) THEN
      IF (position%height_up_to_date) THEN
        fu_height_of_position = position%height
      ELSE
        CALL set_error('positon height not up to date',&
            & 'fu_height_of_position')
      END IF
    ELSE
      CALL set_error('undefined position given'&
	  & ,'fu_height_of_position')
    END IF

  END FUNCTION fu_height_of_position
  

  
  ! ***************************************************************

  LOGICAL FUNCTION fu_position_height_up_to_date(position)
    
    ! Returns height_up_to_date 
    !
    ! Author: Mikhail Sofiev, FMI
    ! 
    ! Imported parameters with intent IN:
    TYPE(silam_grid_position), INTENT(in) :: position

    IF (defined(position)) THEN
      fu_position_height_up_to_date = position%height_up_to_date
    ELSE
      CALL set_error('undefined position given'&
	  & ,'fu_height_of_position')
    END IF

  END FUNCTION fu_position_height_up_to_date

  ! ************************************************************

  FUNCTION fu_geographical_latlon_string(pos, grid) result(string)
    ! 
    ! Returns a string containing the geographical coordinates of a
    ! position. All values are positive and the letter markers S, N,
    ! W and E are used.
    ! grid is the grid in which the position is defined. If not - 
    ! the values will be returned without grid transformation
    ! 
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    CHARACTER (LEN=33) :: string
    
    ! Imported parameters with intent IN:
    TYPE(silam_grid_position), INTENT(in) :: pos
    type(silja_grid), INTENT(in) :: grid
    
    ! Local declarations:
    REAL :: geolat, geolon
    CHARACTER :: sn_char, we_char

    IF(.not.defined(grid))THEN
      CALL set_error('undefined grid given','fu_geographical_latlon_string')
      RETURN
    END IF

    geolat = fu_lat_geographical_from_grid(pos%x, pos%y, grid)
    geolon = fu_lon_geographical_from_grid(pos%x, pos%y, grid)

    IF (geolat>=0.) THEN
      sn_char = 'N'
    ELSE
      sn_char = 'S'
    END IF
    
    IF (geolon>=0.) THEN
      we_char = 'E'
    ELSE
      we_char = 'W'
    END IF
    
    WRITE(string,fmt='(F6.3,A,A,F7.3,A)') ABS(geolat),sn_char,' ',ABS(geolon),we_char

  END FUNCTION fu_geographical_latlon_string
  
  ! ***************************************************************

  REAL FUNCTION  fu_x(pos)
    ! Returns the x relative co-ordinate of a position 
    !
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE
    TYPE(silam_grid_position), INTENT(in) :: pos

    fu_x = pos%x
  END FUNCTION fu_x
  
  ! ***************************************************************
  REAL FUNCTION  fu_y(pos)
    ! Returns the y relative co-ordinate of a position 
    !
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE
    TYPE(silam_grid_position), INTENT(in) :: pos

    fu_y = pos%y
  END FUNCTION fu_y  
  
  ! ***************************************************************
  REAL FUNCTION  fu_z(pos)
    ! Returns the z relative co-ordinate of a position 
    !
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE
    TYPE(silam_grid_position), INTENT(in) :: pos

    fu_z = pos%z
  END FUNCTION fu_z
  
 
  ! ***************************************************************

  FUNCTION fu_time_of_position(position)
    !
    ! Description:
    ! Returns the time of a position.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_time) :: fu_time_of_position
    !
    ! Imported parameters with intent(in):
    TYPE(silam_grid_position) :: position
 
    fu_time_of_position = position%time
    
  END FUNCTION fu_time_of_position
  
  
  
  ! ***************************************************************
  
  REAL FUNCTION fu_horizontal_distance(pos1, pos2) result(dist)
    
    ! Returns the absolute value of the distance between two
    ! positions horizontally. 
    ! The height or time is not concerned here. Result in metres.
    ! The trick: as long as system_grid is defined, we also have 
    ! fields with the grid cell sizes in metres - both for x and y.
    ! Just use it. Attention: the routine is quite crude for large
    ! distances because the cell size is taken in middle point only.
    !
    ! NOTE. Since fields are not available here and thus permanent
    ! fields with cell size are not accessible, trigonometry is used.
    ! This is SLOW routine - do not use it.
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silam_grid_position), INTENT(in) :: pos1, pos2
    !
    ! Local declarations:
    LOGICAL :: corner_in_geographical_latlon, if_south_pole
    REAL :: dx_m, dy_m, dx, dy, &
            & corner_lat_N, corner_lon_E, &
            & southpole_lat_N, southpole_lon_E, & 
            & dx_deg, dy_deg

!    IF(fu_gridtype(meteo_grid) /= lonlat)THEN
!!    IF(fu_gridtype(dispersion_grid) /= lonlat)THEN
!      CALL set_error('So far only latlon grids possible','fu_horizontal_distance')
!      RETURN
!    END IF
!
!!    CALL grid_parameters(dispersion_grid,& ! Get parameters of the grid
!    CALL grid_parameters(meteo_grid,& ! Get parameters of the grid
!                       & corner_lon_E,&
!                       & corner_lat_N, &
!                       & corner_in_geographical_latlon, &
!                       & nx_meteo, ny_meteo, &
!!                       & nx_dispersion, ny_dispersion, &
!                       & southpole_lon_E, & 
!                       & southpole_lat_N, & 
!                       & if_south_pole, &
!                       & dx_deg, dy_deg)
!    IF(error)RETURN

    dx_m =fu_dx_cell_m(meteo_grid, int((pos1%x + pos2%x)*0.5+0.5), int((pos1%y + pos2%y)*0.5+0.5))
    dy_m =fu_dy_cell_m(meteo_grid, int((pos1%x + pos2%x)*0.5+0.5), int((pos1%y + pos2%y)*0.5+0.5))

    if(error)return
 
    dist = SQRT((pos2%x - pos1%x)*dx_m * (pos2%x - pos1%x)*dx_m + &
              & (pos2%y - pos1%y)*dy_m * (pos2%y - pos1%y)*dy_m)

  END FUNCTION fu_horizontal_distance
  


  ! ***************************************************************
  
  REAL FUNCTION fu_pr_diff_of_positions(pos1, pos2) result(diff)
    
    ! Description:
    ! Returns the difference between two positions' heights. Result
    ! in Pa. If the pressure of pos1 is greater, the result is
    ! positive.
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
    TYPE(silam_grid_position), INTENT(in) :: pos1, pos2

    diff = pos1%pressure - pos2%pressure
    
  END FUNCTION fu_pr_diff_of_positions








  ! ***************************************************************
  ! ***************************************************************
  ! 
  !
  !                 TESTS
  !
  !
  ! ***************************************************************
  ! ***************************************************************
     
  SUBROUTINE pos_tests()
    !
    ! Description:
    ! Test stuff for the positions and poles.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    !
    IMPLICIT NONE
    !
    ! Local declarations:
    TYPE(silam_grid_position) :: pos1, pos2
    
    TYPE(silam_pole) :: odd_pole
    
    REAL :: oddlat, oddlon

    PRINT *, ' The meteo south pole : ',&
             & fu_southpole_lat(meteo_pole), &
             & fu_southpole_lon(meteo_pole)
!             & fu_southpole_lat(dispersion_pole), &
!             & fu_southpole_lon(dispersion_pole)
    
    pos1 = fu_set_position(59., 0., 92500., fu_wallclock(), geo_global_grid)
        
    call msg(fu_connect_strings('Position1 geographical:', &
                         & fu_geographical_latlon_string(pos1,geo_global_grid)))
    
    ! this is not legal fortran:
    !PRINT *, 'Position 1 internal: ', pos1
    
    odd_pole = fu_set_pole(north_flag, 29., 180.)
    
    IF (error) RETURN
    
    CALL modify_lonlat(0., 0. , meteo_pole, odd_pole, oddlat, oddlon)
    
    PRINT*, ' The odd pole: ', fu_southpole_lat(odd_pole),&
                          & fu_southpole_lon(odd_pole)

    PRINT *, 'Position1 in odd coordinates: ', oddlat, oddlon

    PRINT *, '****************************************'

    pos1 = fu_set_position(59., 0., 92500., fu_wallclock(), geo_global_grid)

    pos2 = fu_set_position(60., 0., 92500., fu_wallclock(), geo_global_grid)

    CALL print_position(pos1)
    CALL print_position(pos2)

    PRINT *, 'Distance: ', fu_horizontal_distance(pos1, pos2)

  END SUBROUTINE pos_tests


  ! ***************************************************************
  
  SUBROUTINE print_position(pos, grid)

    ! Prints the position. Since the position has to be defined for some grid,
    ! which may not be known sometimes, we assume the following preference order:
    ! - if grid is present and defined - use it
    ! - if system_grid is defined - use it
    ! - use geo_global grid
    !
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    TYPE(silam_grid_position), INTENT(in) :: pos
    type(silja_grid), intent(in), optional :: grid

    IF (defined(pos)) THEN
      !
      ! If grid is present and defined
      !
      if(present(grid))then
        if(defined(grid))then
          call msg(fu_connect_strings('Lat&lon:',&
              & fu_geographical_latlon_string(pos,grid),&
              & ', Time:', fu_str(fu_time(pos)), &
              & ', Pres:'), fu_pressure(pos))
          if(pos%height_up_to_date) call msg('Height up to date:', pos%height)
          return
        endif
      endif
      !
      ! Either absent or undefined grid: try dispersion_grid
      !
      if(defined(meteo_grid))then
!      if(defined(dispersion_grid))then
        call msg(fu_connect_strings('Lat&lon:',&
            & fu_geographical_latlon_string(pos,meteo_grid),&
!            & fu_geographical_latlon_string(pos,dispersion_grid),&
            & ', Time:', fu_str(fu_time(pos)), &
            & ', Pres:'), fu_pressure(pos))
        if(pos%height_up_to_date) call msg('Height up to date:', pos%height)
        return
      endif
      !
      ! Neither is defined: use geo_global_grid
      !
      call msg(fu_connect_strings('Lat&lon:',&
          & fu_geographical_latlon_string(pos,geo_global_grid),&
          & ', Time:', fu_str(fu_time(pos)), &
          & ', Pres:'), fu_pressure(pos))
      if(pos%height_up_to_date) call msg('Height up to date:', pos%height)
    ELSE
      call msg('Undefined position')
    END IF
      
  END SUBROUTINE print_position
  
END MODULE positions_v2


