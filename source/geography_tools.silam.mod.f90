MODULE geography_tools
  !
  ! Provides the definitions and tools for geography-related activities, such
  ! as reprojection, coordinate transformation etc. 
  ! Poles and areas also migrated here now 
  ! Does NOT use any specific grid deficitions but does understand the areas, poles, 
  ! projectsions, etc.
  !
  ! The current person to blame: Julius Vira, FMI julius.vira@fmi.fi
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
  ! Modules used:
  !
  USE silam_namelist
  use silam_units

  IMPLICIT NONE

  ! Public functions and subroutines available in this module:
  ! poles
  PUBLIC fu_set_pole
  PUBLIC defined
  PUBLIC fu_southpole_lat
  PUBLIC fu_southpole_lon
  PUBLIC fu_northpole_lat
  PUBLIC fu_northpole_lon
  public report_as_namelist
  ! areas
  PUBLIC fu_set_area
  PUBLIC area_southwest_corner
  PUBLIC area_northwest_corner
  PUBLIC area_southeast_corner
  PUBLIC area_northeast_corner
  PUBLIC fu_name
  PUBLIC fu_pole
  PUBLIC report
  ! rotation
  public get_rotation
  public get_basis
  public transform_latlon
  PUBLIC modify_lonlat
!  public modify_lonlat_faster
  ! other
  public fu_gc_distance

  ! The private functions and subroutines not to be used elsewhere:
  ! poles
  private fu_set_pole_from_params
  private fu_set_pole_from_namelist
  PRIVATE fu_compare_poles_eq ! interface ==
  private fu_pole_relation
  private report_pole_as_namelist
  PRIVATE fu_pole_defined
  ! areas
  private fu_set_area_from_params
  private fu_set_area_from_nmlist
  PRIVATE fu_area_defined
  PRIVATE fu_compare_areas_eq ! interface ==
  PRIVATE print_area_report
  PRIVATE fu_pole_of_area
  PRIVATE fu_name_of_area

  ! Generic names and operator-interfaces of some functions:
  ! poles
  interface fu_set_pole
    module procedure fu_set_pole_from_params
    module procedure fu_set_pole_from_namelist
  end interface

  interface report_as_namelist
    module procedure report_pole_as_namelist
  end interface

  INTERFACE defined
    MODULE PROCEDURE fu_pole_defined
    MODULE PROCEDURE fu_area_defined
  END INTERFACE
      
  INTERFACE operator(==)
    MODULE PROCEDURE fu_compare_poles_eq
    MODULE PROCEDURE fu_compare_areas_eq
  END INTERFACE
  
  ! areas
  INTERFACE fu_set_area
    MODULE PROCEDURE fu_set_area_from_params
    MODULE PROCEDURE fu_set_area_from_nmlist
  END INTERFACE

  INTERFACE report
    MODULE PROCEDURE print_area_report
  END INTERFACE

  INTERFACE fu_pole
    MODULE PROCEDURE fu_pole_of_area
  END INTERFACE

  INTERFACE fu_name
    MODULE PROCEDURE fu_name_of_area
  END INTERFACE

  
  ! Public types with private components defined in this module:
  ! poles
  TYPE silam_pole 
    PRIVATE
    REAL :: lat, lon
    TYPE(silja_logical) :: defined
  END TYPE silam_pole

  ! The missing codes for these types:
  TYPE(silam_pole), PARAMETER :: pole_missing = silam_pole(real_missing, real_missing, silja_false)
  
  ! The pole for geographical doordinates:
  TYPE(silam_pole), PARAMETER :: pole_geographical = silam_pole(-90., 0., silja_true)
  
  ! just an example ..
  TYPE(silam_pole), PARAMETER, PUBLIC :: pole_hirlam = silam_pole(-30., 0., silja_true) 

  ! Poles of meteo and dispersion grid, set in field_buffer and io_server.
  ! These grids may be not of the lon-lat type, geographical pole is set for anygrid
  TYPE(silam_pole), PUBLIC, SAVE :: meteo_pole = pole_missing  ! 
  TYPE(silam_pole), PUBLIC, SAVE :: dispersion_pole = pole_missing

  ! areas
  ! Public types with private components defined in this module:
  TYPE silam_area
    PRIVATE
    CHARACTER (LEN = clen) :: name
    TYPE(silam_pole) :: pole
    REAL :: lat_south
    REAL :: lat_north
    REAL :: lon_west
    REAL :: lon_east
    TYPE(silja_logical) :: defined
  END TYPE silam_area
  
  ! Missing code:
  TYPE(silam_area), PARAMETER, PUBLIC :: area_missing =&
      & silam_area( 'undefined_area', pole_missing,&
      & real_missing, real_missing, real_missing, real_missing,&
      & silja_false)
  
  ! rotation
  ! Direction of the grid transformation
  integer, public, parameter :: world_to_rotated = 621
  integer, public, parameter :: rotated_to_world = 622

  integer, private, parameter :: iSamePoles = 630
  integer, private, parameter :: iOldPoleGeographical = 631
  integer, private, parameter :: iNewPoleGeographical = 632

  CONTAINS

  !
  ! Routines for the rotation geometry.
  ! 

  !**************************************************************************

      subroutine get_rotation(lat_pole, lon_pole, R)
        implicit none
        ! 
        ! Compute the 3D rotation matrix R corresponding to
        ! transformation from reference frame to frame positive z-axis
        ! is defined by given lat and long of the south pole.
        ! 
        ! Right multiplication -> world to rotated
        ! Left multiplication -> rotated to world 
        !
        ! reference (definition of the rotated lat/lon grid):
        ! http://dss.ucar.edu/docs/formats/grib2/grib2doc/
        !
        real, intent(in) :: lat_pole, lon_pole
        real, dimension(3,3), intent(out) :: R

        real :: alpha, beta, gamma, sa, sb, sg, ca, cb, cg
        real :: lat, lon

        lat = lat_pole * degrees_to_radians
        lon = lon_pole * degrees_to_radians
        
        ! rotation around the geographical polar axis
        alpha = lon
        ! rotate the pole to the correct latitude
        beta = pi/2 + lat
        ! possible rotation around the new polar axis
        gamma = 0.0

        ca = cos(alpha)
        cb = cos(beta)
        cg = cos(gamma)
        sa = sin(alpha)
        sb = sin(beta)
        sg = sin(gamma)

        R(1,1) = ca*cb*cg - sa*sg
        R(2,1) = sa*cb*cg + ca*sg
        R(3,1) = sb*cg
        R(1,2) = -ca*cb*sg - sa*cg
        R(2,2) = -sa*cb*sg + ca*cg
        R(3,2) = -sb*sg
        R(1,3) = -ca*sb
        R(2,3) = -sa*sb
        R(3,3) = cb

      end subroutine get_rotation
 
  !************************************************************************** 
      
      subroutine get_basis(latd, lond, B)
        ! 
        ! Return the 3d-cartesian components of the usual basis
        ! vectors for the tangent plane at latd, lond.
        !
        implicit none
        real, intent(in) :: latd, lond ! degrees
        real, dimension(3,2), intent(out) :: B 

        real :: lat, lon

        lat = latd * degrees_to_radians
        lon = lond * degrees_to_radians

        B(1,1) = -sin(lon)
        B(2,1) = cos(lon)
        B(3,1) = 0.
        B(1,2) = -cos(lon)*sin(lat)
        B(2,2) = -sin(lon)*sin(lat)
        B(3,2) = cos(lat)

      end subroutine get_basis


  !**************************************************************************

  subroutine transform_latlon(lon_from, lat_from, arRotation, direction, lon_new, lat_new)
    !
    ! Transforms the lon-lat *_from coordinates to/from the rotated grid system defined 
    ! by the arRotation rotation matrix - following the direction request.
    !
        ! Right multiplication -> world to rotated
        ! Left multiplication -> rotated to world 


    implicit none

    ! Imported parameters
    real, intent(in) :: lon_from, lat_from
    integer, intent(in) :: direction
    real, dimension(3,3) :: arRotation
    real, intent(out) :: lon_new, lat_new

    ! Local variables
    real, dimension(3) :: r1, r2
    real :: lon_rad, lat_rad

    lon_rad = lon_from * degrees_to_radians
    lat_rad = lat_from * degrees_to_radians

    r1(1) = cos(lon_rad)*cos(lat_rad)
    r1(2) = sin(lon_rad)*cos(lat_rad)
    r1(3) = sin(lat_rad)

    if(direction == rotated_to_world)then
      r2 = matmul(arRotation, r1)
    elseif(direction == world_to_rotated)then
      r2 = matmul(r1, arRotation)
    else
      call msg('Unknown grid transformation direction:', direction)
      call set_error('Unknown grid transformation direction', 'transform_latlon')
      return
    endif

    lon_new = atan2(r2(2), r2(1)) * radians_to_degrees
    lat_new = asin(r2(3)) * radians_to_degrees

  end subroutine transform_latlon

  
   ! ***************************************************************

    SUBROUTINE modify_lonlat(lat_from, lon_from, pole, pole_new, lat_new, lon_new)
    
    ! Description:
    ! Converts the lat-lon pair given in lat-lon
    ! -coordinates defined by the pole, to a new lat-lon pair
    ! defined by pole_new. Can be used to convert positions between
    ! modified and geograhical, or between two different modified
    ! coordinate systems.
    !
    ! This subroutine uses the parameter pole_geographical. No check
    ! of definition of poles are made, since this routine is private
    ! and so used in this module only.
    !
    ! Method:
    ! Uses the rotation-matrix based algorithm
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silam_pole), INTENT(in) :: pole, pole_new
    REAL, INTENT(in) :: lat_from, lon_from

    ! Imported parameters with intent OUT:
    REAL, INTENT(out) :: lat_new, lon_new

    ! Local declarations:
    REAL :: geolat,  geolon, lon_rad, lat_rad, alpha, beta, gamma, ca, cb, cg, sa, sb, sg, pole_lat, pole_lon
    real, dimension(3,3) :: R
    real, dimension(3) :: r1, r2


    !----------------------------------------
    !
    ! 1. Check if the poles happen to be the same already.
    select case(fu_pole_relation(pole, pole_new))
      case(iSamePoles)
!    IF (pole == pole_new) THEN
      lat_new = lat_from
      lon_new = lon_from
    
      !----------------------------------------
      !
      ! 2. Convert geographical to modified
      case(iOldPoleGeographical)
!    ELSE IF (pole == pole_geographical) THEN


      pole_lat = fu_southpole_lat(pole_new) * degrees_to_radians
      pole_lon = fu_southpole_lon(pole_new) * degrees_to_radians
 
      ! Euler angles. Some tricks needed, since this choice doesn't
      ! coincide with lat lon coordinates.
    
        ! rotation around the geographical polar axis
        alpha = pole_lon !_new
        ! rotate the pole to the correct latitude
        beta = pi/2 + pole_lat ! _new
        ! possible rotation around the new polar axis
        gamma = 0.0

        ca = cos(alpha)
        cb = cos(beta)
        cg = cos(gamma)
        sa = sin(alpha)
        sb = sin(beta)
        sg = sin(gamma)

        R(1,1) = ca*cb*cg - sa*sg
        R(2,1) = sa*cb*cg + ca*sg
        R(3,1) = sb*cg
        R(1,2) = -ca*cb*sg - sa*cg
        R(2,2) = -sa*cb*sg + ca*cg
        R(3,2) = -sb*sg
        R(1,3) = -ca*sb
        R(2,3) = -sa*sb
        R(3,3) = cb

      lon_rad = lon_from * degrees_to_radians
      lat_rad = lat_from * degrees_to_radians

      r1(1) = cos(lon_rad)*cos(lat_rad)
      r1(2) = sin(lon_rad)*cos(lat_rad)
      r1(3) = sin(lat_rad)

      r2 = matmul(r1, R)

      lon_new = atan2(r2(2), r2(1)) * radians_to_degrees
      lat_new = asin(r2(3)) * radians_to_degrees

      
      !----------------------------------------
      !
      ! 3. Convert modified to geographical
      case(iNewPoleGeographical)
!    ELSE IF (pole_new == pole_geographical) THEN

      pole_lat = fu_southpole_lat(pole) * degrees_to_radians
      pole_lon = fu_southpole_lon(pole) * degrees_to_radians
 
      ! Euler angles. Some tricks needed, since this choice doesn't
      ! coincide with lat lon coordinates.
    
        ! rotation around the geographical polar axis
        alpha = pole_lon
        ! rotate the pole to the correct latitude
        beta = pi/2 + pole_lat
        ! possible rotation around the new polar axis
        gamma = 0.0

        ca = cos(alpha)
        cb = cos(beta)
        cg = cos(gamma)
        sa = sin(alpha)
        sb = sin(beta)
        sg = sin(gamma)

        R(1,1) = ca*cb*cg - sa*sg
        R(2,1) = sa*cb*cg + ca*sg
        R(3,1) = sb*cg
        R(1,2) = -ca*cb*sg - sa*cg
        R(2,2) = -sa*cb*sg + ca*cg
        R(3,2) = -sb*sg
        R(1,3) = -ca*sb
        R(2,3) = -sa*sb
        R(3,3) = cb

      lon_rad = lon_from * degrees_to_radians
      lat_rad = lat_from * degrees_to_radians

      r1(1) = cos(lon_rad)*cos(lat_rad)
      r1(2) = sin(lon_rad)*cos(lat_rad)
      r1(3) = sin(lat_rad)

      r2 = matmul(R, r1)

      lon_new = atan2(r2(2), r2(1)) * radians_to_degrees
      lat_new = asin(r2(3)) * radians_to_degrees

      case default
!    ELSE
      

      pole_lat = fu_southpole_lat(pole) * degrees_to_radians
      pole_lon = fu_southpole_lon(pole) * degrees_to_radians
 
      ! Euler angles. Some tricks needed, since this choice doesn't
      ! coincide with lat lon coordinates.
    
        ! rotation around the geographical polar axis
        alpha = pole_lon
        ! rotate the pole to the correct latitude
        beta = pi/2 + pole_lat
        ! possible rotation around the new polar axis
        gamma = 0.0

        ca = cos(alpha)
        cb = cos(beta)
        cg = cos(gamma)
        sa = sin(alpha)
        sb = sin(beta)
        sg = sin(gamma)

        R(1,1) = ca*cb*cg - sa*sg
        R(2,1) = sa*cb*cg + ca*sg
        R(3,1) = sb*cg
        R(1,2) = -ca*cb*sg - sa*cg
        R(2,2) = -sa*cb*sg + ca*cg
        R(3,2) = -sb*sg
        R(1,3) = -ca*sb
        R(2,3) = -sa*sb
        R(3,3) = cb

      lon_rad = lon_from * degrees_to_radians
      lat_rad = lat_from * degrees_to_radians

      r1(1) = cos(lon_rad)*cos(lat_rad)
      r1(2) = sin(lon_rad)*cos(lat_rad)
      r1(3) = sin(lat_rad)

      r2 = matmul(R, r1)

      lon_rad = atan2(r2(2), r2(1))
      lat_rad = asin(r2(3))


      pole_lat = fu_southpole_lat(pole_new) * degrees_to_radians
      pole_lon = fu_southpole_lon(pole_new) * degrees_to_radians
 
      ! Euler angles. Some tricks needed, since this choice doesn't
      ! coincide with lat lon coordinates.
    
        ! rotation around the geographical polar axis
        alpha = pole_lon !_new
        ! rotate the pole to the correct latitude
        beta = pi/2 + pole_lat !_new
        ! possible rotation around the new polar axis
        gamma = 0.0

        ca = cos(alpha)
        cb = cos(beta)
        cg = cos(gamma)
        sa = sin(alpha)
        sb = sin(beta)
        sg = sin(gamma)

        R(1,1) = ca*cb*cg - sa*sg
        R(2,1) = sa*cb*cg + ca*sg
        R(3,1) = sb*cg
        R(1,2) = -ca*cb*sg - sa*cg
        R(2,2) = -sa*cb*sg + ca*cg
        R(3,2) = -sb*sg
        R(1,3) = -ca*sb
        R(2,3) = -sa*sb
        R(3,3) = cb

      r1(1) = cos(lon_rad)*cos(lat_rad)
      r1(2) = sin(lon_rad)*cos(lat_rad)
      r1(3) = sin(lat_rad)

      r2 = matmul(r1, R)

      lon_new = atan2(r2(2), r2(1)) * radians_to_degrees
      lat_new = asin(r2(3)) * radians_to_degrees

    END select !IF
    
  END SUBROUTINE modify_lonlat !_faster
  
  
  
  !==============================Poles=======================================
  
  
    FUNCTION fu_set_pole_from_params(pole_flag, lat, lon) result(pole)

    ! Description:
    ! Sets a value for a pole-system, NOT A SINGLE SOUTH- OR
    ! NORTHPOLE.
    ! 
    IMPLICIT NONE
    
    TYPE(silam_pole) :: pole ! The return value of this function 
    
    ! Imported parameters with intent IN:
    INTEGER, INTENT(in) :: pole_flag ! indicates, which
    ! one of the pole-coordinates are given. Allowed values are
    ! south and north.
    REAL, INTENT(in) :: lat, lon  ! geographical latitude and
    ! longitude of modified or unmodified pole
    
    !  Local declarations:
    REAL :: templat, templon

    ! -----------------------------------------------------------

    templat = fu_set_latitude(lat, north_flag)
    IF (error) RETURN
    
    templon = lon
    IF (error) RETURN
    
    SELECT CASE (pole_flag)

    CASE (north_flag)! Convert internally always to south pole:

      pole%lat = -1 * templat
      pole%lon = templon + 180.   
      IF (pole%lon > 180.) pole%lon = pole%lon - 360.

    CASE (south_flag) ! the internal format already
      
      pole%lat = templat
      pole%lon = templon
      
    CASE default
      
      call msg('pole_flag: ', pole_flag)
      CALL set_error('wrong pole_flag, must be south or north','fu_set_pole_from_params')
      pole = pole_missing
      RETURN
      
    END SELECT
    
    pole%defined = silja_true ! everything ok
    
  END FUNCTION fu_set_pole_from_params


  !********************************************************************* 

  FUNCTION fu_set_pole_from_namelist(nlSetup) result(pole)
    !
    ! Sets a value for a pole-system, NOT A SINGLE SOUTH- OR
    ! NORTHPOLE. Input data: silam_namelist, so standard SILAM
    ! conditions are applicable
    !
    IMPLICIT NONE
    
    TYPE(silam_pole) :: pole ! The return value of this function 
    
    type(Tsilam_namelist), intent(in) :: nlSetup

    pole%lon = fu_content_real(nlSetup,'lon_s_pole')
    pole%lat = fu_content_real(nlSetup,'lat_s_pole')

!    if((pole%lon .eps.0.) .and. (pole%lat.eps.0.))then
!      call report(nlSetup)
!      call set_error('Both lon and lat of the pole are zero','fu_set_pole_from_namelist')
!      return
!    endif
    if((pole%lon .eps.real_missing) .or. (pole%lat.eps.real_missing))then
      call report(nlSetup)
      call set_error('Both lon and lat of the pole are missing','fu_set_pole_from_namelist')
      return
    endif

    pole%defined = silja_true ! everything ok
    
  END FUNCTION fu_set_pole_from_namelist
 

  !***********************************************************************

  subroutine report_pole_as_namelist(pole, iUnit)
    !
    ! Writes the pole report in a form of namelist to the given unit
    !
    implicit none

    ! Imported parameters
    TYPE(silam_pole), intent(in) :: pole
    integer, intent(in) :: iUnit

    if(.not. (pole%defined == silja_true))then
      write(iUnit,*)'lon_s_pole = UNDEFINED_POLE'
      call set_error('Undefined pole given','report_pole_as_namelist')
      return
    endif

    write(iUnit, *)'lon_s_pole = ', pole%lon
    write(iUnit, *)'lat_s_pole = ', pole%lat

  end subroutine report_pole_as_namelist


  ! ***************************************************************
  
  LOGICAL FUNCTION fu_compare_poles_eq(pole1, pole2)
    !
    ! Description:
    ! Compares two poles and returns a ture value if equal. This
    ! function is used through an interface (pole1 == pole2)
    !  Note! pole_missing /= pole missing
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    TYPE(silam_pole), INTENT(in) :: pole1, pole2

    IF ((ABS(pole1%lat - pole2%lat) < 1.0e-6) .and. &
      & (ABS(pole1%lon - pole2%lon) < 1.0e-6) .and. &
      & pole1%defined == silja_true .and.  pole2%defined == silja_true) THEN
      fu_compare_poles_eq = .true.
    ELSE
      fu_compare_poles_eq = .false.
    END IF

  END FUNCTION fu_compare_poles_eq
  
  
  !****************************************************************
  
  integer function fu_pole_relation(pole_old, pole_new)
    !
    ! Checks whether the poles are the same or if one of them is geographical. 
    ! New and Old notations just refer to the function modify_lonlat, which is the only 
    ! user of this one
    ! Returns the 4-position switch
    !
    implicit none
    
    ! Imported parameters with intent IN:
    TYPE(silam_pole), INTENT(in) :: pole_old, pole_new

    if((ABS(pole_old%lat - pole_new%lat) < 1.0e-6) .and. &  ! poles are equal - see previous function
     & (ABS(pole_old%lon - pole_new%lon) < 1.0e-6))then
      fu_pole_relation = iSamePoles
      
    elseif((ABS(pole_old%lat - pole_geographical%lat) < 1.0e-6) .and. &  ! pole_old == pole_geographical
         & (ABS(pole_old%lon - pole_geographical%lon) < 1.0e-6))then
      fu_pole_relation = iOldPoleGeographical
      
    elseif((ABS(pole_new%lat - pole_geographical%lat) < 1.0e-6) .and. &  ! pole_new == pole_geographical
         & (ABS(pole_new%lon - pole_geographical%lon) < 1.0e-6))then
      fu_pole_relation = iNewPoleGeographical
    else
      fu_pole_relation = int_missing
    endif
    
  end function fu_pole_relation

  ! ***************************************************************

  LOGICAL FUNCTION fu_pole_defined(pole)
    !
    ! Description:
    ! Returns a true value, if the pole has been given a value using
    ! the function fu_set_pole.
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silam_pole), INTENT(in) :: pole

    fu_pole_defined = fu_true(pole%defined)
        
  END FUNCTION fu_pole_defined

  ! ***************************************************************
  
  REAL FUNCTION fu_southpole_lat(pole)
    
    ! Description:
    ! Returns the geographical latitude of
    ! the southpole of a pole-system. Northern latitudes are positive.
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    TYPE(silam_pole), INTENT(in) :: pole
    
    fu_southpole_lat = pole%lat
    
  END FUNCTION fu_southpole_lat
  

  ! ***************************************************************
  
  REAL FUNCTION fu_southpole_lon(pole)
    
    ! Description:
    ! Returns the geographical longitude of
    ! the southpole of a pole-system.
    ! Eastern longitudes are positive.
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    TYPE(silam_pole), INTENT(in) :: pole
    
    fu_southpole_lon = pole%lon
    
  END FUNCTION fu_southpole_lon
  
  ! ***************************************************************
  
  REAL FUNCTION fu_northpole_lat(pole)
    
    IMPLICIT NONE
    TYPE(silam_pole), INTENT(in) :: pole

    fu_northpole_lat = -1.0 * pole%lat
    
  END FUNCTION fu_northpole_lat
  


  ! ***************************************************************
  
  REAL FUNCTION fu_northpole_lon(pole)
    
    IMPLICIT NONE
    TYPE(silam_pole), INTENT(in) :: pole

    fu_northpole_lon = fu_scale_longitude(pole%lon + 180.)
    
  END FUNCTION fu_northpole_lon
 

  ! ***************************************************************

  FUNCTION fu_set_area_from_params( name_of_area, &
      & lat_south, sn_flag_south,  lat_north, sn_flag_north, &
      & lon_west, we_flag_west,  lon_east, we_flag_east, &
      & pole) result(area)
    
    ! Description:
    ! Sets a value for an area, which defines the bordeds of a latlon
    ! -box, but not the grid inside it.
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silam_area) :: area
    !
    ! Imported parameters with intent(in):
    CHARACTER (LEN=*), INTENT(in) :: name_of_area
    REAL, INTENT(in) :: lat_south, lat_north, lon_west, lon_east
    INTEGER, INTENT(in) :: sn_flag_south,& ! possible values given in module globals
                         & sn_flag_north,&
                         & we_flag_west,&
                         & we_flag_east
    TYPE(silam_pole), INTENT(in) :: pole !defined the lat-lon system
    ! of the coordinates above
    
    !----------------------------------------
    !
    ! 1. Set values and scale.
    !    ---------------------

    IF (defined(pole)) THEN
      area%pole = pole
    ELSE
      CALL set_error('undefined pole given','fu_set_latlon_area')
      RETURN
    END IF

    area%lat_south = fu_set_latitude(lat_south, sn_flag_south)
    area%lat_north = fu_set_latitude(lat_north, sn_flag_north)
    area%lon_west = lon_west
    area%lon_east = lon_east

    IF (error) RETURN

    IF (area%lat_south > area%lat_north) THEN
      call msg('area%lat_south, area%lat_north', area%lat_south, area%lat_north)
      CALL set_error('area lat_south bigger than lat_north', 'fu_set_lonlat_area_from_params')
      RETURN
    END IF

    area%name = TRIM(ADJUSTL(name_of_area))
    area%defined = fu_set_true()
    
  END FUNCTION fu_set_area_from_params


  !*****************************************************************

  FUNCTION fu_set_area_from_nmlist(nlSetup) result(area)
    !
    ! Sets a value for an area, which defines the bordeds of a latlon
    ! -box, but not the grid inside it. The parameters are set from 
    ! the namelist with a wide use of default values
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silam_area) :: area
    !
    ! Imported parameters with intent(in):
    type(Tsilam_namelist), intent(in) :: nlSetup

    !
    ! Local variables
    character(len=clen) :: strTmp
    integer :: status
    REAL :: lat_south, lat_north, lon_west, lon_east

    strTmp = fu_content(nlSetup,'area_borders')
    read(unit=strTmp,fmt=*,iostat=status) lat_south, lat_north, lon_west, lon_east
    if(status /= 0)then
      call set_error(fu_connect_strings('Failed to read area borders from:',strTmp), &
                   & 'fu_set_lonlat_area_from_nmlist')
      return
    endif

    !
    ! So far, in the namelists we allow only geographical areas with normal pole
    !
    area%pole = pole_geographical

    !
    ! In the namelists north and east are always positive
    !
    area%lat_south = fu_set_latitude(lat_south, north_flag) 
    area%lat_north = fu_set_latitude(lat_north, north_flag)
    area%lon_west = lon_west
    area%lon_east = lon_east

    IF (error) RETURN

    IF (area%lat_south > area%lat_north) THEN
      call msg('area%lat_south, area%lat_north', area%lat_south, area%lat_north)
      CALL set_error('area lat_south bigger than lat_north', 'fu_set_lonlat_area_from_nmlist')
      RETURN
    END IF
    IF (lon_east < lon_west) THEN
      call msg('lon_west lon_east', lon_west, lon_east)
      CALL set_error('area lon_west bigger than lon_east', 'fu_set_lonlat_area_from_nmlist')
      RETURN
    END IF

    area%name = fu_content(nlSetup,'area_title')
    area%defined = fu_set_true()
    
  END FUNCTION fu_set_area_from_nmlist


  ! ***************************************************************

  LOGICAL FUNCTION fu_area_defined(area)
 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silam_area), INTENT(in) :: area

    fu_area_defined = fu_true(area%defined)
    
  END FUNCTION fu_area_defined


  ! ***************************************************************

  FUNCTION fu_name_of_area(area)
    
    IMPLICIT NONE
    !
    ! Return valuye of this finction:
    CHARACTER (LEN=clen) :: fu_name_of_area

    ! Imported parameters with intent(in):
    TYPE(silam_area), INTENT(in) :: area

    fu_name_of_area = area%name
    
  END FUNCTION fu_name_of_area


  ! ***************************************************************

  FUNCTION fu_pole_of_area(area)
    
    IMPLICIT NONE
    !
    ! Return valuye of this finction:
    type(silam_pole) :: fu_pole_of_area

    ! Imported parameters with intent(in):
    TYPE(silam_area), INTENT(in) :: area

    IF (defined(area)) THEN
      fu_pole_of_area = area%pole
    ELSE
      fu_pole_of_area = pole_missing
    END IF

  END FUNCTION fu_pole_of_area

   ! ***************************************************************

  SUBROUTINE area_southwest_corner(area, pole, southwest_lat, southwest_lon)
    
    ! Returns southwest-corner coordinates in lat-lon-system defined
    ! by the pole. Northern latitudes and eastern longitudes are
    ! positive.

    IMPLICIT NONE
    
    ! Imported parameters with intent IN:
    TYPE(silam_area), INTENT(in) :: area
    TYPE(silam_pole) :: pole

    ! Imported parameters with intent OUT:
    REAL, INTENT(out) :: southwest_lat, southwest_lon
    
    CALL modify_lonlat(area%lat_south,&
                     & area%lon_west,&
                     & area%pole,&
                     & pole,&
                     & southwest_lat,&
                     & southwest_lon)
   
  END SUBROUTINE area_southwest_corner


  ! ***************************************************************

  SUBROUTINE area_northwest_corner(area, pole, northwest_lat, northwest_lon)
    
    ! Returns northwest-corner coordinates in lat-lon-system defined
    ! by the pole. Northern latitudes and eastern longitudes are
    ! positive.

    IMPLICIT NONE
    
    ! Imported parameters with intent IN:
    TYPE(silam_area), INTENT(in) :: area
    TYPE(silam_pole) :: pole

    ! Imported parameters with intent OUT:
    REAL, INTENT(out) :: northwest_lat, northwest_lon
    
    CALL modify_lonlat(area%lat_north,&
                     & area%lon_west,&
                     & area%pole,&
                     & pole,&
                     & northwest_lat,&
                     & northwest_lon)
   
  END SUBROUTINE area_northwest_corner


  ! ***************************************************************

  SUBROUTINE area_southeast_corner(area, pole, southeast_lat, southeast_lon)
    
    ! Returns southeast-corner coordinates in lat-lon-system defined
    ! by the pole. Northern latitudes and eastern longitudes are
    ! positive.
    IMPLICIT NONE
    
    ! Imported parameters with intent IN:
    TYPE(silam_area), INTENT(in) :: area
    TYPE(silam_pole) :: pole

    ! Imported parameters with intent OUT:
    REAL, INTENT(out) :: southeast_lat, southeast_lon
    
    CALL modify_lonlat(area%lat_south,&
                     & area%lon_east,&
                     & area%pole,&
                     & pole,&
                     & southeast_lat,&
                     & southeast_lon)
   
  END SUBROUTINE area_southeast_corner

 
  ! ***************************************************************

  SUBROUTINE area_northeast_corner(area, pole, northeast_lat, northeast_lon)
    
    ! Returns northeast-corner coordinates in lat-lon-system defined
    ! by the pole. Northern latitudes and eastern longitudes are
    ! positive.
    IMPLICIT NONE
    
    ! Imported parameters with intent IN:
    TYPE(silam_area), INTENT(in) :: area
    TYPE(silam_pole) :: pole

    ! Imported parameters with intent OUT:
    REAL, INTENT(out) :: northeast_lat, northeast_lon
    
    CALL modify_lonlat(area%lat_north,&
                     & area%lon_east,&
                     & area%pole,&
                     & pole,&
                     & northeast_lat,&
                     & northeast_lon)
   
  END SUBROUTINE area_northeast_corner

 
  ! ***************************************************************

  LOGICAL FUNCTION fu_compare_areas_eq(area1, area2) result(eq)
    
    ! Returns true if areas are equal.
    IMPLICIT NONE
    
    ! Imported parameters with intent(in):
    TYPE(silam_area), INTENT(in) :: area1, area2

    eq = .false.
    IF (.NOT.defined(area1)) THEN
      IF (.NOT.defined(area2)) THEN
	    eq = .true. ! Both undefined, so equal:
	    RETURN
      ELSE
	    RETURN
      END IF
    ELSE
      IF (.NOT.defined(area2))RETURN
    END IF
    
    IF (.NOT.(area1%pole == area2%pole)) RETURN

    IF ((area1%lat_south .eps. area2%lat_south).and.&
	  & (area1%lat_north .eps. area2%lat_north).and.&
	  & (area1%lon_west .eps. area2%lon_west).and.&
	  & (area1%lon_east .eps. area2%lon_east)) eq = .true.

  END FUNCTION fu_compare_areas_eq


  ! ***************************************************************

  SUBROUTINE print_area_report(area)
    
    ! Print contents of an area on screen for test purposes 
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silam_area), INTENT(in) :: area

    IF (.NOT.defined(area  )) THEN
      call msg(' **** Area report :undefined area *************')
      RETURN
    END IF
    
    call msg(fu_connect_strings('********* Area report of ',area%name,'************'))

    IF (area%pole == pole_geographical) THEN
        call msg('Area boundaries in geographical coordinates.')
    ELSE
        call msg('Area boundaries in modified coordinates.')
        WRITE(*, fmt='(A, 2F8.2)') 'lat and lon of southpole: ',&
                                   & fu_southpole_lat(area%pole),&
                                   & fu_southpole_lon(area%pole)
        WRITE(run_log_funit, fmt='(A, 2F8.2)') 'lat and lon of southpole: ',&
                                              & fu_southpole_lat(area%pole),&
                                              & fu_southpole_lon(area%pole)
    END IF

    WRITE(*, fmt='(F8.2, A, F8.2, A)') area%lat_south, ' ... ',&
                                       & area%lat_north, ' lat N boundaries '
    WRITE(*, fmt='(F8.2, A, F8.2, A)') area%lon_west, ' ... ',&
                                       & area%lon_east, ' lon E boundaries'
    WRITE(run_log_funit, fmt='(F8.2, A, F8.2, A)') area%lat_south, ' ... ',&
                                                   & area%lat_north,' lat N boundaries'
    WRITE(run_log_funit, fmt='(F8.2, A, F8.2, A)') area%lon_west, ' ... ',&
                                                   & area%lon_east, ' lon E boundaries'

    call msg(' ******** End-of-area report ****************')

  END SUBROUTINE print_area_report
  
  !*******************************************************************
  !*******************************************************************
  
  real function fu_gc_distance(x1, x2, y1, y2) result(dist)
    ! Returns great circle distance between two lon, lat points
    implicit none
    
    ! Imported parameters with intent(in):
    real, intent(in) :: x1, x2, y1, y2

    ! local variables
    real :: yr1, yr2, dx, dsg

    yr1 = y1 * pi / 180.
    yr2 = y2 * pi / 180.
    dx = abs(x1-x2) * pi / 180. 
    dsg = (sin((yr2-yr1)/2))**2 + cos(yr1)*cos(yr2)*(sin(dx/2))**2
    dsg = 2. * asin(sqrt(dsg))
    dist = earth_radius * dsg 

  end function fu_gc_distance
    
END MODULE geography_tools

