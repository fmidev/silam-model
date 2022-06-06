MODULE grids_geo

  ! Description:
  ! 1. Provides the definition an tools for SILJA_GRID, which
  ! defines a horizontal grid. This grid may be either a modified
  ! latitude-longitude or ps (polarstereographic) type. All the
  ! numerical
  ! weather model areas, their sub-areas, and the calculation
  ! grids of eulerian dispersion models is defined with this type,
  ! and ONLY this type in this system.
  !
  ! Provides function for area transformations of a field from one
  ! grid to another, when both are of the type SILJA_GRID.
  !
  ! Provides functions for horizontal interpolation to a point
  ! inside a grid of type SILJA_GRID.
  !
  !
  ! Original code: Mika Salonoja
  ! Author: Mikhail Sofiev, FMI mikhail.sofiev@fmi.fi
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
  USE geography_tools
  !!$use omp_lib

  IMPLICIT NONE

  ! Public functions and subroutines available in this module:
  
  ! Defining the grid
  !
  public fu_set_grid
  PUBLIC fu_set_lonlat_grid
  public fu_set_any_grid
  public setAnygridParam
  public setAgLonlatFlds
  PUBLIC fu_lonlat_grid_from_area
  public fu_area_from_grid
  public release_grid

  ! Information about the grid
  !
  PUBLIC defined
  PUBLIC report
  public report_as_namelist
  PUBLIC fu_number_of_gridpoints
!  PUBLIC grid_parameters
  PUBLIC lonlat_grid_parameters
  PUBLIC grid_dimensions
  public adjust_grid_dimensions
  PUBLIC fu_pole
  PUBLIC fu_gridtype
  PUBLIC fu_lon_native_from_grid
  PUBLIC fu_lat_native_from_grid
  PUBLIC fu_lon_geographical_from_grid
  PUBLIC fu_lat_geographical_from_grid
  public fu_geolons_fld
  public fu_geolats_fld
  public fu_cos_map_rot_fld
  public fu_sin_map_rot_fld
  public fu_dy_fld_m
  public fu_dx_fld_m
  public fu_dx_cell_m
  public fu_dy_cell_m
  public fu_dx_cell_deg
  public fu_dy_cell_deg
  public fu_cell_size
  public fu_ifLonGlobal
  public fu_ifLatGlobal
  public fu_stdSilamGrid
  public make_grid_lon_global
  public make_grid_lat_global
  public fu_boundary_grid
  public fu_ifPolarCapPossible
  public fu_ifPoleIncluded


  ! Grid transformations, re-projections and other 
  ! grid-related actions
  !
  public grid_data_hor_select_new_grid
  PUBLIC grid_data_horizontal_select
  public remap_field
!  public grid_transf_via_interp_struct
  PUBLIC fu_grid_containing_area
  public cut_grid_size     ! Reduce one grid to the size of another
  public cut_empty_lines   ! reduce the grid size using the given array
  public fu_pick_subgrid   ! A subdomain out of the big grid
  public fill_overlap_flag_array
  PUBLIC coriolis_parameters
  public grid_shift_indices
  public project_point_to_grid 
  public make_minimal_grid
  public extend_grid_to_coordinates
  public reposition_global_grid_lon
  public reposition_global_grid_lat
  public fu_point_inside_grid
  public set_grid_from_lores_templ

  public fu_if_mesh_intrpolatable

  !
  ! Interpolation structure
  !
  public fu_horiz_interp_struct  ! returns if exists or creates if not
  public fu_gridTo
  public fu_grid_index
  public fu_nCoefs
  public get_area_limits
  public fu_interpType
  public get_coefs
  public fu_ifWindRotationNeeded
  public fu_n_non_zero_output_cells
  public fu_rotation

  ! Numerical horizontal derivatives
  !
  PUBLIC ddx_of_field
  PUBLIC ddy_of_field
  public laplacian

  ! Grid comparison
  !
  PUBLIC fu_grids_match_closely ! Check for Arakawa shifts
  public fu_grids_correspond   ! Same grids, apart from size
  public fu_grids_arakawa_correspond ! Same or Arakawa-shifted grids, apart from size
  public fu_if_grid_covered ! If one grid covers the other one
  public fu_area_coverage ! Checks if some grid is covered by another grid
  public SubArea_Chk

  ! The private functions and subroutines
  !
  private fu_set_grid_from_namelist
  private fu_set_grid_from_params
  private make_global_lonlat_grid
  private fu_ag_param_index
  PRIVATE fu_compare_grids_eq ! interface ==
  PRIVATE fu_grid_defined
  PRIVATE fu_lonlat_grid_containing_area
  private fu_anygrid_containing_area
  PRIVATE adjust_grid_to_sample
  PRIVATE print_grid_report
  private report_grid_as_namelist
!  PRIVATE fu_grid_shift_index
  private write_anygrid_2_gradsfile
  private rd_grd_pars_from_gradsfile
  private fu_new_ag_param_index
  private cut_grid_size_grid
  private cut_grid_size_params

  PRIVATE fu_pole_of_grid
  private fu_cell_size_midpoint
  private fu_cell_size_n_th_cell  
  private fu_cell_size_x_y_cell  
  private fu_dx_x_y_cell_m
  private fu_dy_x_y_cell_m
  private fu_dx_x_y_cell_deg
  private fu_dy_x_y_cell_deg
  private project_point_to_grid_lonlat
  private project_point_to_grid_xy

  private fu_gridFrom_from_interp_struct
  private fu_gridTo_from_interp_struct
  private fu_nCoefs_interp_str_horiz
  private fu_interpType_interp_str_horiz
  private get_coefs_interp_str_horiz
  private get_coefs_interp_str_cell_horiz
  private fu_rotation_interp_str_horiz
  private fu_rotation_interp_str_cell_horiz
  private fu_rotation_interp_str_comp_horiz

  ! Generic names and operator-interfaces of some functions:

  INTERFACE fu_set_grid
    MODULE PROCEDURE fu_set_grid_from_namelist
    module procedure fu_set_grid_from_params
!    module procedure fu_set_any_grid
  END INTERFACE

  interface fu_cell_size
    module procedure fu_cell_size_midpoint
    module procedure fu_cell_size_n_th_cell
    module procedure fu_cell_size_x_y_cell
  end interface 

  interface fu_dx_cell_m
    module procedure fu_dx_x_y_cell_m
  end interface
    
  interface fu_dx_cell_deg
    module procedure fu_dx_x_y_cell_deg
  end interface

  interface fu_dy_cell_m
    module procedure fu_dy_x_y_cell_m
  end interface

  interface fu_dy_cell_deg
    module procedure fu_dy_x_y_cell_deg
  end interface

  interface project_point_to_grid
    module procedure project_point_to_grid_lonlat
    module procedure project_point_to_grid_xy
  end interface 

  interface cut_grid_size
    module procedure cut_grid_size_grid
    module procedure cut_grid_size_params
  end interface

  interface fu_gridFrom
    module procedure fu_gridFrom_from_interp_struct
  end interface 

  interface fu_gridTo
    module procedure fu_gridTo_from_interp_struct
  end interface 

  interface fu_nCoefs
    module procedure fu_nCoefs_interp_str_horiz
  end interface 

  interface fu_interpType
    module procedure fu_interpType_interp_str_horiz
  end interface 

  interface get_coefs
    module procedure get_coefs_interp_str_horiz
    module procedure get_coefs_interp_str_cell_horiz
  end interface 
  
  interface fu_rotation
    module procedure fu_rotation_interp_str_horiz
    module procedure fu_rotation_interp_str_cell_horiz
    module procedure fu_rotation_interp_str_comp_horiz
  end interface 

  INTERFACE defined
    MODULE PROCEDURE fu_grid_defined
  END INTERFACE

  INTERFACE report
    MODULE PROCEDURE print_grid_report
  END INTERFACE

  interface report_as_namelist
    module procedure report_grid_as_namelist
  end interface

  INTERFACE fu_pole
    MODULE PROCEDURE fu_pole_of_grid
  END INTERFACE

  INTERFACE operator(==)
    MODULE PROCEDURE fu_compare_grids_eq
  END INTERFACE

  private deallocate_horiz_interp_struct


  ! ***************************************************************

  ! ***************************************************************
  !
  !
  !        HERE'S THE LAT-LON GRID
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  TYPE silja_lonlat_grid
    REAL :: sw_corner_modlon, sw_corner_modlat  ! lon-lat coordinates of 
    ! the first gridpoint given in modified system defined by the
    ! grid's own pole
    INTEGER :: nx, ny ! number of points in x and y, respectively
    TYPE(silam_pole) :: pole 
    REAL :: dx_deg, dy_deg ! Grid distance in
    ! west-east and south-north, respectively.

  END TYPE silja_lonlat_grid

  PRIVATE silja_lonlat_grid

  TYPE(silja_lonlat_grid), PARAMETER, PRIVATE :: lonlat_grid_missing = &
            & silja_lonlat_grid(real_missing, real_missing, &
                              & int_missing, int_missing, pole_missing, &
                              & real_missing, real_missing)


  ! ***************************************************************

  ! ***************************************************************
  !
  !
  !        HERE'S THE POLARSTEREOGRAPHIC GRID
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  TYPE silja_ps_grid ! NOT FINISHED YET!!!
    INTEGER :: nx, ny  ! number of points in we and sn, respectively
    REAL :: vertical_longitude ! = the perpendicular longitude 
    REAL :: grid_dist_m ! grid distance in metres
    LOGICAL :: lonlat_of_sw_corner_given ! If true, the geographical
    ! (or normal) coordinates of grid's south-west corner are given.
    ! If false, the gridpoint coordinates of the north pole are given.
    REAL :: sw_geolat, sw_geolon ! The geographical
    ! coordinates of grid's south-west corner in case of the switch
    ! lonlat_of_sw_corner_given == .true.
    REAL :: x_of_northpole, y_of_northpole ! In case of
    ! polarstereographic grid's lonlat_of_sw_corner_given == .false.
  END TYPE silja_ps_grid

  PRIVATE silja_ps_grid

  TYPE(silja_ps_grid), PARAMETER, PRIVATE :: &
      & ps_grid_missing = &
      & silja_ps_grid(int_missing, int_missing, &
      & real_missing, real_missing, .false.,&
      & real_missing, real_missing, real_missing, real_missing)


  ! ***************************************************************

  ! ***************************************************************
  !
  !
  !        HERE'S THE GAUSS-KRUEGER GRID
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  REAL, PARAMETER, PRIVATE :: gk_origo_lat = 0.
  REAL, PARAMETER, PRIVATE :: gk_origo_lon = 0.

  TYPE silja_gauss_krueger_grid
    REAL :: sw_dist_x, sw_dist_y ! [m] south-west corner location in
    ! metres from the general origo
    REAL :: dx, dy ! [m] grid-dist inside grid
    INTEGER :: nx, ny ! number of gridpoints
  END TYPE silja_gauss_krueger_grid

  PRIVATE silja_gauss_krueger_grid

  TYPE(silja_gauss_krueger_grid), PARAMETER, PRIVATE ::&
      & gk_grid_missing =&
      & silja_gauss_krueger_grid(&
      & real_missing, real_missing, real_missing, real_missing,&
      & int_missing, int_missing)

  ! ****************************************************************
  ! ****************************************************************
  !
  !
  !       HERE's the ANY GRID defined via field of lons and lats
  !
  !
  ! ****************************************************************
  ! ****************************************************************

  type silam_any_grid_param
    !
    ! Type for storing the grid parameters in a global variable
    !
    real, dimension(:), allocatable :: xC,yC
    real, dimension(:), allocatable :: x3dC, y3dC, z3dC !!3d cartesian geocentric coordinates in units of earth_radius
    real, dimension(:), allocatable :: dx, dy, sin_map_rot, cos_map_rot ! sizes, rotation angle
    type(silja_logical) :: defined 
  end type silam_any_grid_param

  integer, public, parameter :: max_nbr_any_grids = 11 
 
  type(silam_any_grid_param), dimension(max_nbr_any_grids), public, target, save :: pAnyGrdParam 
  integer, dimension(max_nbr_any_grids), private, save :: anyGrdParamRefCount = int_missing

  type silam_any_grid
    integer :: indParam  ! index in the pAnyGrdParam array
    integer :: nx, ny
  end type silam_any_grid

  private silam_any_grid

  type(silam_any_grid), parameter, private :: ag_missing = &
                   & silam_any_grid(int_missing, int_missing, int_missing)


  ! ***************************************************************

  ! ***************************************************************
  !
  !
  !        HERE'S THE GENERAL DEFINITION OF ALL GRIDS
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  ! Public types with private components defined in this module:

  TYPE silja_grid ! defines one type of horizontal 2D-grid
    PRIVATE
    CHARACTER (LEN=80) :: name = "Undefined" ! Name of grid; propably the same as
    ! the variable defined ( TYPE(grid) :: huuhaa 
    ! huuhaa%name = huuhaa)
    INTEGER :: gridtype = int_missing ! the parameters below
    TYPE(silja_lonlat_grid) :: lonlat
    TYPE(silja_gauss_krueger_grid) :: gk
    TYPE(silja_ps_grid) :: ps 
    type(silam_any_grid) :: ag
  END TYPE silja_grid

  INTEGER, PARAMETER :: lonlat = 31
  INTEGER, PARAMETER :: polarstereographic = 32
  integer, parameter :: anygrid = 33

  TYPE(silja_grid), PARAMETER, PUBLIC :: grid_missing = silja_grid('missing',&
                                                                 & int_missing,&
                                                                 & lonlat_grid_missing,&
                                                                 & gk_grid_missing, &
                                                                 & ps_grid_missing, &
                                                                 & ag_missing)

  TYPE(silja_grid), PARAMETER, PUBLIC :: geo_global_grid = silja_grid( &
                        & 'global geographical grid',&
                        & lonlat,&
                        & silja_lonlat_grid(-180., -90.,360,180,pole_geographical,1.,1.), &
                        & gk_grid_missing, &
                        & ps_grid_missing, &
                        & ag_missing)
      
  TYPE(silja_grid), PARAMETER, PUBLIC :: coarse_geo_global_grid = silja_grid(&
                                               & 'coarse global geographical grid',&
                                               & lonlat,&
                          & silja_lonlat_grid(-180.,-90., 36,18,pole_geographical,10.,10.), &
                                               & gk_grid_missing, &
                                               & ps_grid_missing, &
                                               & ag_missing)
! Replaced by grib_scan_mode type in grib_io
!  INTEGER, PARAMETER, PUBLIC :: hirlam_grib_scanning_mode = 64 
!  INTEGER, PARAMETER, PUBLIC :: ecmwf_grib_scanning_mode = ??  

!  ! Grid shift indicators: (Arakawa A,C,D)
!  INTEGER, PARAMETER, PRIVATE :: same_corner = 1
!  INTEGER, PARAMETER, PRIVATE :: samelat_halfgrid_lon = 2 ! for two grids
!  
!  ! the latitudes (y) are the same, but longitudes (x) have a
!  ! shift of half grid square
!  INTEGER, PARAMETER, PRIVATE :: samelon_halfgrid_lat = 3 ! vice versa above
!  INTEGER, PARAMETER, PRIVATE :: halfgrid_lon_lat = 4 ! Both indices are shifted
!
!  ! or, if direction of shifting is important:
!  INTEGER, PARAMETER, PRIVATE :: plus_half = 1    ! i -> i+1/2
!  INTEGER, PARAMETER, PRIVATE :: minus_half = -1  ! i -> i-1/2
!  INTEGER, PARAMETER, PRIVATE :: no_shift = 99

  ! The grid is lat-global if the edges at south and north are closer to the poles then
  REAL, PARAMETER, PRIVATE :: lat_global_threshold = 10.0 ! degrees


  ! Operating with grids, one can be interested in two criteria: if grid covers
  ! the other one or is fully inside it
  integer, parameter, public :: cover_the_grid_area = 721
  integer, parameter, public :: inside_the_grid_area = 722


  ! ***************************************************************

  ! ***************************************************************

  !
  !
  !                    THE METEOROLOGICAL GRID (former SYSTEM_GRID)
  !                    =======================
  !
  !                    THE DISPERSION GRID
  !                    ===================
  !
  ! The meteo_grid is introduced, which makes SILAM flexible for 
  ! non-HIRLAM meteorological data. Meteorological pole now is just
  ! a derivative of the system grid (see module poles). It can be each 
  ! time taken from meteo_grid but this would take too much resources.
  !
  ! Now all meteofields will be forced to be in meteo_grid. Meteo_grid 
  ! is directly copied from the meteorological file with, possibly,
  ! some reductions of the covered area. 
  ! The same is true for the meteo_vertical (see module levels), which
  ! is selected from the meteo files as a sub-set of the best-fitting
  ! vertical structure availabel there.
  !
  ! Dispersion_grid will be used for positions, which are all in RELATIVE 
  ! coorinates of the dispersion_grid. It also serves the Eulerian 
  ! environment - together with the dispersion_vertical (see module levels)
  ! they form the main 3D spatial cube and its coordinates system.
  !
  TYPE(silja_grid),PUBLIC,TARGET,SAVE :: meteo_grid = grid_missing
  INTEGER,PUBLIC,SAVE :: nx_meteo = int_missing, ny_meteo = int_missing
  INTEGER,PUBLIC,SAVE :: fs_meteo = int_missing

  TYPE(silja_grid),PUBLIC,TARGET,SAVE :: dispersion_grid = grid_missing
  TYPE(silja_grid),PUBLIC,POINTER,SAVE :: dispersion_gridPtr
  INTEGER,PUBLIC,SAVE :: nx_dispersion = int_missing, ny_dispersion = int_missing, &
                       & fs_dispersion = int_missing

  real, public, dimension(:,:), allocatable  :: disp_grid_size_x, disp_grid_size_y

  TYPE(silja_grid),PUBLIC,TARGET,SAVE :: wholeMPIdispersion_grid = grid_missing
  TYPE(silja_grid),PUBLIC,POINTER,SAVE :: wholeMPIdispersion_gridPtr
  INTEGER,PUBLIC,SAVE :: nx_wholeMPIdispersion = int_missing, ny_wholeMPIdispersion = int_missing, &
                       & fs_wholeMPIdispersion = int_missing


  TYPE(silja_grid),PUBLIC,TARGET,SAVE :: output_grid = grid_missing
  TYPE(silja_grid),PUBLIC,POINTER,SAVE :: output_gridPtr
  INTEGER,PUBLIC,SAVE :: nx_output = int_missing, ny_output = int_missing, &
                       & fs_output = int_missing

  ! ***************************************************************


  !================================================================================
  !
  ! Horizontal interpolation structure. The interpoaltion from one grid to 
  ! another may happen thousands times.
  ! Below are the structures for the coefficiens of interpolation 
  ! between two grids of any type. Compute them once and then simply use
  ! the weighted-summing cycles to get the interpolated values.
  ! 
  type THorizInterpCells    ! All cells in the structure
    !
    ! Dimensions are all the same and depend on type of interpolation
    !
    integer, dimension(:,:,:), pointer :: indX, indY ! values in the gridFrom   (nCoefs, nxTo, nyTo)
    real, dimension(:,:,:), pointer :: weight, weight_X, weight_Y ! of the specific index  (nCoefs, nxTo, nyTo)
    logical, dimension(:,:), pointer :: ifValid => null() ! values in the gridFrom   (nxTo, nyTo)
  end type THorizInterpCells

  type THorizInterpOneCell    ! single cell in the structure
    !
    ! Dimensions are all the same and depend on type of interpolation
    !
    integer, dimension(:), pointer :: indX, indY ! values in the gridFrom   (nCoefs)
    real, dimension(:), pointer :: weight, weight_X, weight_Y ! of the specific index  (nCoefs)
    logical :: ifValid = .false.
  end type THorizInterpOneCell

  type THorizInterpStruct     ! The structure itself
    !  private
    logical :: ifInterpolate = .false., ifGridFromGlobal = .false.
    type(silja_grid) :: gridFrom = grid_missing, gridTo = grid_missing
    integer :: interp_type, nCoefs, ixStartFrom, iyStartFrom, ixEndFrom, iyEndFrom, &
                                  & ixStartTo, iyStartTo, ixEndTo, iyEndTo, iOutside
    integer, dimension(:,:,:), pointer :: indX => null() ! indices in the gridFrom   (nCoefs, nxTo, nyTo)
    integer, dimension(:,:,:), pointer :: indY => null() ! indices in the gridFrom   (nCoefs, nxTo, nyTo)
    logical, dimension(:,:), pointer :: ifValid => null() ! values in the gridFrom   (nxTo, nyTo)
    real, dimension(:,:,:), pointer :: weight => null(), weight_X  => null(), weight_Y => null() 
                                        ! of the specific index  (nCoefs, nxTo, nyTo)
    ! 
    ! Support for wind rotations
    !
    real, pointer, dimension(:,:,:,:) :: rotation ! (2, 2, nxTo, nyTo)
    logical :: ifRotation = .false.
  end type THorizInterpStruct
  type (THorizInterpStruct), public,parameter :: HorizInterpStruct_missing_par = &
             & THorizInterpStruct(.False., .false., grid_missing, grid_missing, int_missing,  int_missing,&
             & int_missing, int_missing, int_missing, int_missing, int_missing,&
             & int_missing, int_missing, int_missing, int_missing, null(),&
             & null(),null(), null(),null(), null(), null(),.false.)

  type (THorizInterpStruct), public, target :: HorizInterpStruct_missing = HorizInterpStruct_missing_par

  type THorizInterpStructPtr
    type(THorizInterpStruct), pointer :: ptr
  end type THorizInterpStructPtr
  private THorizInterpStructPtr
  type (THorizInterpStructPtr), public,parameter :: HorizInterpStructPtr_missing = &
                &  THorizInterpStructPtr(null())
  !
  ! The run will have a pool of interpolation structures, where everyone will be free to find the favorite one
  !
  integer, parameter :: horizInterpPoolSize = 60
  type THorizInterpPool
    ! horizontal interpolation structure pointer
    type(THorizInterpStructPtr), dimension(horizInterpPoolSize) :: pHIS  
    ! Reference count is incremented when pointer to a structure is returned
    integer, dimension(horizInterpPoolSize) :: refCount = 0
    integer :: nHIS = 0
  end type THorizInterpPool
  private THorizInterpPool
  
  type(THorizInterpPool), private, save :: HorizInterpPool  ! the bunch of interpolation structures
  

  
  
CONTAINS


  ! **************************************************************
  ! **************************************************************
  !
  !     Grid creation
    !
  ! **************************************************************
  ! **************************************************************


  function fu_set_grid_from_namelist(nlSetup)result(grid)
    ! 
    ! Sets a grid from the namelist. So far, the only grid type implemented 
    ! is the lonlat grid
    !
    IMPLICIT NONE

    ! A small value
    real, parameter :: eps = 1.0e-5 ! Small value for single presision

    ! The return value of this function:
    TYPE(silja_grid) :: grid

    ! Imported parameter with intent IN
    type(Tsilam_namelist), intent(in) :: nlSetup

    ! Local variables
    integer :: jTmp, nx,ny
    real :: lon_end, dx, lat_end, dy
    character(len = fnlen) :: grdfnm
    character(len = *), parameter :: sub_name = 'fu_set_grid_from_namelist'

    grid = grid_missing ! inital value 

    !-------------------------------------------------
    !
    ! General variables.
    !
    grid%name = fu_content(nlSetup,'grid_title')

    !
    ! The rest is grid-type specific
    !
    if(fu_str_u_case(fu_content(nlSetup,'grid_type')) == 'LON_LAT_GLOBAL')then
      !
      ! global grid is set in an automated way.
      !
      call make_global_lonlat_grid(fu_content_int(nlSetup,'nx'), fu_content_int(nlSetup,'ny'), &
                                 & fu_content_real(nlSetup,'dx'), fu_content_real(nlSetup,'dy'), &
                                 & grid)
      
    elseif(fu_str_u_case(fu_content(nlSetup,'grid_type')) == 'LON_LAT')then

      grid%lonlat%pole = fu_set_pole(nlSetup)
      if(.not.defined(grid%lonlat%pole)) then
         call report(nlSetup)
         call set_error('Cannot find the grid pole',sub_name)
         return
      endif

      !-------------------------------------------------
      !
      ! The first grid point is always in correct co-ordinates - no rotation is
      ! expected. In particular, we assume here that for rotated grid the corner
      ! will be given in the rotated co-ordinates
      !
      grid%lonlat%sw_corner_modlon = fu_content_real(nlSetup,'lon_start')
      grid%lonlat%sw_corner_modlat = fu_content_real(nlSetup,'lat_start')

      if((grid%lonlat%sw_corner_modlon == real_missing) .or. &
       & (grid%lonlat%sw_corner_modlat == real_missing)) then
         call report(nlSetup)
         call set_error('Cannot find the grid corner',sub_name)
         return
      endif
      IF (error) RETURN

      !-------------------------------------------------
      !
      ! Number of points and grid distance can be defined in several ways
      !
      !Longitude
      nx  = fu_content_int(nlSetup,'nx')
      dx  = fu_content_real(nlSetup,'dx')
      lon_end = fu_content_real(nlSetup,'lon_end')
      jTmp  = 0
      if (nx == int_missing) jTmp = jTmp + 1 
      if (dx == real_missing) jTmp = jTmp + 1
      if (lon_end == real_missing) jTmp = jTmp + 1
      if (jTmp > 1) then !More than one undefined
        call report(nlSetup)
        call set_error("at least two of nx dx and lon_end should be defined", sub_name)
        return
      endif
      if (nx == int_missing) then
          nx = nint((lon_end - grid%lonlat%sw_corner_modlon)/dx) + 1
      elseif (dx == real_missing) then
          dx = ( lon_end - grid%lonlat%sw_corner_modlon) / (nx-1)
      elseif (lon_end /= real_missing) then
        if (nint((lon_end - grid%lonlat%sw_corner_modlon)/dx) + 1 - nx /= 0) then
          call report(nlSetup)
          call set_error("inconsistent nx dx and lon_end", sub_name)
          return
        endif
      endif
      grid%lonlat%nx = nx
      grid%lonlat%dx_deg = dx
      !adjust the origin 
      grid%lonlat%sw_corner_modlon = fu_scale_longitude(grid%lonlat%sw_corner_modlon)

      !latitude
      ny  = fu_content_int(nlSetup,'ny')
      dy  = fu_content_real(nlSetup,'dy')
      lat_end = fu_content_real(nlSetup,'lat_end')
      jTmp  = 0
      if (ny == int_missing) jTmp = jTmp + 1 
      if (dy == real_missing) jTmp = jTmp + 1
      if (lat_end == real_missing) jTmp = jTmp + 1
      if (jTmp > 1) then
        call set_error("at least two of ny dy and lat_end should be defined", sub_name)
        return
      endif
      if (ny == int_missing) then
          ny = nint((lat_end - grid%lonlat%sw_corner_modlat)/dy) + 1
      elseif (dy == real_missing) then
          dy = ( lat_end - grid%lonlat%sw_corner_modlat) / (ny-1)
      elseif (lon_end /= real_missing) then
        if ( nint((lat_end - grid%lonlat%sw_corner_modlat)/dy) + 1 - ny /= 0) then
          call report(nlSetup)
          call set_error("inconsistent ny dy and lat_end", sub_name)
          return
        endif
      endif
      grid%lonlat%ny = ny
      grid%lonlat%dy_deg = dy


      !! Grid at least 2x2 and >100m step...
      if (dx< 1e-5 .or. dy < 1e-5 .or. nx < 2 .or. ny < 2 ) then
        call set_error("dx< 1e-5 .or. dy < 1e-5 .or. nx < 2 .or. ny < 2 ", sub_name)
        return
      endif

      grid%gridtype = lonlat

    elseif(fu_str_l_case(fu_content(nlSetup,'grid_type')) == 'anygrid')then

      nx = fu_content_int(nlSetup,'nx')
      if(nx <= 0)then
        call set_error('strange number of x gridcells',sub_name)
        return
      endif
      grid%ag%nx = nx
      ny = fu_content_int(nlSetup,'ny')
      if(ny <= 0)then
        call set_error('strange number of y gridcells',sub_name)
        return
      endif
      grid%ag%ny = ny
      grdFNm = fu_content(nlSetup,'grid_file')
      if(trim(grdFNm) == '')then
        call set_error('problem with filename',sub_name)
        return
      endif
      grid%gridtype = anygrid
      grid%ag%indParam = fu_new_ag_param_index()
      call rd_grd_pars_from_gradsFile(grdFNm, grid)
      if(error)return


    else
      call set_error(fu_connect_strings('Only lon_lat and any-grids available, not:', &
                                      & fu_content(nlSetup,'grid_type')), &
                   & sub_name)
      return
    endif

  end function fu_set_grid_from_namelist


  !***********************************************************************************************
  
  function fu_pick_subgrid(grid_in, iSubgrid, jSubgrid, &
                         & nXSubgrids, nYSubgrids, &
                         & division_hint, &
                         & xOffset, yOffset) &
    & result(grid_out)
    !
    ! Picks a subgrid from the big grid, which is split to nX x nY subgrids. Note that the subgrids
    ! can have different sizes - but not more than +- 1 cell in each direction
    !
    implicit none

    ! The return value of this function:
    TYPE(silja_grid) :: grid_out

    ! Imported parameter with intent IN
    type(silja_grid), intent(in) :: grid_in
    integer, intent(in) :: iSubgrid, jSubgrid, nXSubgrids, nYSubgrids
    real, dimension(:), intent(in) :: division_hint ! Relative fractions of sub-domains


    ! The index offsets are returned to the caller
    integer, intent(out) :: xOffset, yOffset

    ! Local variables
    integer :: nx_big, ny_big, nx_my, ny_my
    real ::  fTmp
    integer, dimension (max_divisions) :: ysizes


    call grid_dimensions(grid_in, nx_big, ny_big)

    
    xOffset = iSubgrid * nx_big / nxSubgrids
    nx_my =  (iSubgrid + 1 ) * nx_big / nxSubgrids - xOffset


    yOffset = 0
    if (all(division_hint(1:nYSubgrids)>0)) then !Use guided division for Y
      fTmp = sum(division_hint(1:nYSubgrids)) !Norm
      if (jSubgrid > 0) yOffset = int(sum(division_hint(1:jSubgrid)/fTmp*ny_big)) 
      ny_my = int(sum(division_hint(1:jSubgrid+1)/fTmp*ny_big)) - yOffset
    else
      yOffset = jSubgrid * ny_big / nYSubgrids
      ny_my =  (jSubgrid + 1 ) * ny_big / nYSubgrids - yOffset
    endif

    ! First, copy the input grid
    !
    grid_out = grid_in
    !
    ! then just cut and shift the grid dimension-wise
    !
    select case(grid_out%gridtype)
      case(lonlat)
        ! Reduce size and move the starting points of grids accordingly
        grid_out%lonlat%nx = nx_my
        grid_out%lonlat%sw_corner_modlon = grid_out%lonlat%sw_corner_modlon &
             & + grid_out%lonlat%dx_deg * xOffset

        grid_out%lonlat%ny = ny_my
        grid_out%lonlat%sw_corner_modlat = grid_out%lonlat%sw_corner_modlat &
             & + grid_out%lonlat%dy_deg * yOffset

      case default
        call set_error('Only lonlat grids are currently supported by the MPI version so far...', &
                     & 'fu_set_grid_mpi_from_nameslist')
        return
    end select  ! grid type

  end function fu_pick_subgrid


  !*************************************************************************

  function fu_set_grid_from_params(chName,iType, pole, sw_corner_modlon, sw_corner_modlat, &
                                 & nx, ny, dx, dy)result(grid)
    ! 
    ! Sets a grid from the namelist. So far, the only grid type implemented 
    ! is the lonlat grid
    !
    IMPLICIT NONE

    ! A small value
    real, parameter :: eps = 1.0e-5 ! Small value for single presision

    ! The return value of this function:
    TYPE(silja_grid) :: grid

    ! Imported parameters
    integer, intent(in) :: iType, nx, ny
    character(len=*), intent(in) :: chName
    type(silam_pole), intent(in) :: pole
    real, intent(in) :: sw_corner_modlon, sw_corner_modlat, dx, dy

    ! Local variables
    integer :: iTmp
    real :: fTmp, fTmp2

    grid = grid_missing ! inital value 

    !-------------------------------------------------
    !
    ! General variables.
    !
    grid%name = chName

    if(iType == lonlat)then

      grid%lonlat%pole = pole
      grid%lonlat%sw_corner_modlon = fu_scale_longitude(sw_corner_modlon)
      grid%lonlat%sw_corner_modlat = fu_scale_latitude(sw_corner_modlat)

      grid%lonlat%nx = nx
      grid%lonlat%ny = ny

      grid%lonlat%dx_deg = dx
      grid%lonlat%dy_deg = dy

      ! -------------------------------------------------
      !
      !  6. Check that all grid cells have longitude in [-180,180] range.
      !     May be, the starting point has to be rotated to one evolution.
      !     --------------
      if(grid%lonlat%dx_deg * (grid%lonlat%nx-1) > 360. + eps)then
        call report(grid)
        call set_error('More than 360 degrees in the grid','fu_set_grid_from_params')
        return
      end if
      if(grid%lonlat%sw_corner_modlon + grid%lonlat%dx_deg * (grid%lonlat%nx-1) > 180.+ eps)then
        grid%lonlat%sw_corner_modlon = grid%lonlat%sw_corner_modlon - 360.
      end if

      ! -------------------------------------------------
      !
      !  7. Everything ok.
      !     --------------

      grid%gridtype = lonlat

    else
      call set_error('Only lon_lat grids available', 'fu_set_grid_from_params')
      return
    endif

  end function fu_set_grid_from_params


  !*************************************************************************

  function fu_staggered_grid(chName, srcgrid, if_stag_x, if_stag_y) result(grid)
    ! 
    ! Makes staggered grid by increasing its size by one node 
    ! Works only for lonlat grids
    !
    IMPLICIT NONE

    ! A small value
    real, parameter :: eps = 1.0e-5 ! Small value for single presision

    ! The return value of this function:
    TYPE(silja_grid) :: grid

    ! Imported parameters
    TYPE(silja_grid), intent(in) :: srcgrid
    character(len=*), intent(in) :: chName
    logical, intent(in) :: if_stag_x, if_stag_y



    if(srcgrid%gridtype == lonlat)then
       grid = srcgrid ! inital value 
       grid%name = chName
       if (if_stag_x) then 
          grid%lonlat%sw_corner_modlon = grid%lonlat%sw_corner_modlon - 0.5*grid%lonlat%dx_deg
!          Extra cell in global x-staggered grids makes vector operations simpler
!          if (.not. fu_ifLonGlobal(grid)) then
                grid%lonlat%nx = grid%lonlat%nx + 1
!          endif
        endif
       if (if_stag_y) then 
          grid%lonlat%sw_corner_modlat = grid%lonlat%sw_corner_modlat - 0.5*grid%lonlat%dy_deg
          grid%lonlat%ny = grid%lonlat%ny + 1
        endif

      if(grid%lonlat%dx_deg * (grid%lonlat%nx-1) > 360. + 0.01*grid%lonlat%dx_deg)then
        call msg('Initial and staggered grids:')
        call report(srcgrid)
        call report(grid)
        call set_error('More than 360 degrees in the grid','fu_staggered_grid')
        return
      end if

    else
      grid = grid_missing
      call set_error('Only lon_lat grids available', 'fu_staggered_grid')
      
      return
    endif

  end function fu_staggered_grid


  !*********************************************************************

  subroutine report_grid_as_namelist(grid, iUnit, fName)
    !
    ! Stores the grid as the namelist to the given file, which must be open
    !
    implicit none

    type(silja_grid), intent(in) :: grid
    integer, intent(in) :: iUnit
    CHARACTER (LEN=*), optional :: fName
    CHARACTER (LEN=fnlen) :: fNm

    if(present(fName))then
      fNm = fu_connect_strings(fName, '.grid.grads')
    else
      fNm = 'grid.grads' 
    endif

    write(iUnit,*)'grid_title = ', adjustl(trim(grid%name))

    select case(grid%gridtype)
      case(lonlat)
        write(iUnit, *)'grid_type = LON_LAT'

        call report_as_namelist(grid%lonlat%pole,iUnit)
        if(error)return

        write(iUnit,'(A,F15.7)')'lon_start = ', grid%lonlat%sw_corner_modlon
        write(iUnit,'(A,F15.7)')'lat_start = ', grid%lonlat%sw_corner_modlat
        write(iUnit,'(A,I5)')'nx = ', grid%lonlat%nx
        write(iUnit,'(A,I5)')'ny = ', grid%lonlat%ny
        write(iUnit,'(A,F15.7)')'dx = ', grid%lonlat%dx_deg
        write(iUnit,'(A,F15.7)')'dy = ', grid%lonlat%dy_deg

        write(iUnit,'(A)')'resol_flag = 128'
        write(iUnit,'(A)')'ifReduced = 0'
        write(iUnit,'(A)')'earth_flag = 0'
        write(iUnit,'(A)')'wind_component = 0 '
        write(iUnit,'(A)')'reduced_nbr_str = 0'

        write(iUnit,'(A)')'lat_pole_stretch = 0. '
        write(iUnit,'(A)')'lon_pole_stretch = 0.'

      case(anygrid)
        !
        ! Have to write a binary file with fields lon, lat, dx, dy, sin_map_rot, cos_map_rot 
        !
        write(iUnit, *)'grid_type = anygrid'
        write(iUnit,'(A12, A)') 'grid_file = ', adjustl(trim(fu_connect_strings(fNm, '.ctl')))
        write(iUnit,'(A,F15.7)')'lon_start = ', pAnyGrdParam(grid%ag%indParam)%xc(1)
        write(iUnit,'(A,F15.7)')'lat_start = ', pAnyGrdParam(grid%ag%indParam)%yc(1)
        write(iUnit,'(A,I5)')'nx = ', grid%ag%nx
        write(iUnit,'(A,I5)')'ny = ', grid%ag%ny
        write(iUnit,'(A,F15.7)')'dx_m = ', pAnyGrdParam(grid%ag%indParam)%dx(1)
        write(iUnit,'(A,F15.7)')'dy_m = ', pAnyGrdParam(grid%ag%indParam)%dy(1)

        call write_anygrid_2_gradsfile(grid, fNm)

      case default
        call report(grid)
        call set_error('Non-supported grid type','report_grid_as_namelist')
        write(iUnit, '(A)')'grid_type = UNSUPPORTED'

    end select ! grid type

  end subroutine report_grid_as_namelist


  ! ************************************************************


  FUNCTION fu_set_lonlat_grid (name, &
                             & corner_lon, corner_lat, &
                             & corner_in_geographical_lonlat, &
                             & number_of_points_x, number_of_points_y, &
                             & pole, & 
                             & dx_deg, dy_deg) result(grid)

    ! Description:
    ! Sets a latitude-longitude grid. 
    ! See the desription of the type
    ! silja_grid in the beginning of this module.
    !
    IMPLICIT NONE

    ! A small value
    real, parameter :: eps = 1.0e-5 ! Small value for single presision

    ! The return value of this function:
    TYPE(silja_grid) :: grid

    ! Imported parameters with intent IN:
    CHARACTER (LEN=*), INTENT(in) :: name

    REAL, INTENT(in) :: corner_lat, corner_lon ! lat and lon of the
    ! first grid point. 

    LOGICAL, INTENT(in) :: corner_in_geographical_lonlat ! if true,
    ! the coordinates of the first point are given in geographical
    ! coordinates. If false, the coordinates are given in the
    ! modified coordinates defined by the pole.

    INTEGER, INTENT(in) :: number_of_points_x, number_of_points_y

    REAL, INTENT(in) ::  dx_deg, dy_deg

    TYPE(silam_pole) :: pole  ! Use fu_set_pole to give a
    ! value for this before setting the grid

    ! -------------------------------------------------
    !
    ! 1. General variables.
    !    ------------------

    grid = grid_missing ! inital value 
    grid%name = name


    ! -------------------------------------------------
    !
    ! 2. The pole of grid (if pole=pole_geographical, then grid is
    !     unmodified lat-lon).
    !    -----------------------------------------------------------

    IF (defined(pole)) THEN
      grid%lonlat%pole = pole
    ELSE
      CALL set_error ('no value for pole, use fu_set_pole','fu_set__lonlat_grid')
      grid = grid_missing
      RETURN
    END IF


    ! -------------------------------------------------
    !
    ! 3. The first grid point.
    !    ---------------------

    IF ((corner_in_geographical_lonlat).and.&
        & (.NOT.(pole == pole_geographical))) THEN
      ! the corner coordinates have to be modified:
      CALL modify_lonlat &
          & (fu_scale_latitude(corner_lat), &
          & fu_scale_longitude(corner_lon), &
          & pole_geographical, pole, &
          & grid%lonlat%sw_corner_modlat, grid%lonlat%sw_corner_modlon)
    ELSE 
      ! already in modified :
      grid%lonlat%sw_corner_modlat = fu_scale_latitude(corner_lat)
      grid%lonlat%sw_corner_modlon = fu_scale_longitude(corner_lon)
    END IF

    IF (error) RETURN


    ! -------------------------------------------------
    !
    ! 4. Number of points.
    !    -----------------

    IF (number_of_points_x > 0) THEN
      grid%lonlat%nx = number_of_points_x ! EW
    ELSE
      CALL set_error('number_of_points_x is not positive','fu_set_lonlat_grid')
      RETURN
    END IF

    IF (number_of_points_y > 0) THEN
      grid%lonlat%ny = number_of_points_y ! SN
    ELSE
      CALL set_error('number_of_points_y is not positive ','fu_set_lonlat_grid')
      RETURN
    END IF


    ! -------------------------------------------------
    !
    !  5. The grid distance
    !     -----------------

    IF (dx_deg.eps.0.) THEN
      CALL set_error('cannot have dx_deg=0' ,'fu_set_lonlat_grid')
    ELSE
      grid%lonlat%dx_deg = dx_deg
    END IF

    IF (dy_deg.eps.0.) THEN
      CALL set_error('cannot have dy_deg=0','fu_set_lonlat_grid')
    ELSE
      grid%lonlat%dy_deg = dy_deg
    END IF

!    ! -------------------------------------------------
!    !
!    !  6. Check that all grid cells have longitude in [-180,180] range.
!    !     May be, the starting point has to be rotated to one evolution.
!    !     --------------
    if(dx_deg * (number_of_points_x-1) > 360. + eps)then
      call set_error('More than 360 degrees in the grid','fu_set_lonlat_grid')
      return
    end if
!    if(grid%lonlat%sw_corner_modlon + dx_deg * (number_of_points_x-1) > 180.+ eps)then
!      grid%lonlat%sw_corner_modlon = grid%lonlat%sw_corner_modlon - 360.
!    end if


    ! -------------------------------------------------
    !
    !  7. Everything ok.
    !     --------------

    grid%gridtype = lonlat


  END FUNCTION fu_set_lonlat_grid


  ! ***************************************************************


  FUNCTION fu_lonlat_grid_from_area(area, ifMetres, &
                                        & gridsize_x, gridsize_y) result(grid)
    !
    ! Sets a horizontal grid inside given area havig given gridsize.
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_grid) :: grid
    !
    ! Imported parameters with intent(in):
    TYPE(silam_area), INTENT(in) :: area
    REAL, INTENT(in) :: gridsize_x, gridsize_y
    logical, intent(in) :: ifMetres  ! Unit of the grid cell size

    ! Local declarations:
    INTEGER :: number_of_points_x, number_of_points_y
    REAL :: dx_deg, dy_deg
    REAL :: south_lat, north_lat, west_lon, east_lon

    call area_southwest_corner(area, fu_pole(area), south_lat, west_lon)
    call area_northeast_corner(area, fu_pole(area), north_lat, east_lon)

    ! If resolution is given in km, select the right dx/dy in degrees to cover the area exactly; 
    ! If in degrees, keep dx/dy leaving the accuracy of coverage to the user
    if(ifMetres)then   ! if kilometres
      dx_deg = fu_dx_m_to_deg(gridsize_x, ((south_lat+north_lat)/2))
      dy_deg = fu_dy_m_to_deg(gridsize_y)
      number_of_points_x = NINT((east_lon-west_lon)/dx_deg)
      number_of_points_y = NINT((north_lat-south_lat)/dy_deg)
      dx_deg = (east_lon-west_lon)/number_of_points_x
      dy_deg = (north_lat-south_lat)/number_of_points_y
    else           ! degrees
      dx_deg = gridsize_x
      dy_deg = gridsize_y
      number_of_points_x = NINT((east_lon-west_lon)/dx_deg)
      number_of_points_y = NINT((north_lat-south_lat)/dy_deg)
    endif

    grid = fu_set_lonlat_grid (fu_name(area), &
                             & west_lon, south_lat, .true., &
                             & number_of_points_x, number_of_points_y, &
                             & fu_pole(area), & 
                             & dx_deg, dy_deg)

  END FUNCTION fu_lonlat_grid_from_area

  ! ***************************************************************



  integer function fu_ag_param_index(nc, lon, lat, ifNew) result(index)
    !
    ! Given lon & lat fields, finds either the matching anygrid or an index for storing a
    ! new one.
    !
    implicit none

    real, dimension(:), intent(in), optional :: lon,lat
    integer, intent(in) :: nc
    logical, intent(out) :: ifNew

    ! Local variables
    integer :: i, j

    do i=1, max_nbr_any_grids-1
      if (anyGrdParamRefCount(i) == int_missing) then
        ! never allocated
        pAnyGrdParam(i)%defined = silja_false
        anyGrdParamRefCount(i) = 0
      endif

      if(.not. pAnyGrdParam(i)%defined == silja_true)then
        ! no longer used
        if (anyGrdParamRefCount(i) /= 0) then
          call set_error('anygrid not defined but refCount not 0', 'fu_ag_param_index')
          return
        end if
        anyGrdParamRefCount(i) = 1
        pAnyGrdParam(i)%defined = silja_true
        index = i
        ifNew = .true.
        return
      else
        if(.not.(allocated(pAnyGrdParam(i)%xc) .and. allocated(pAnyGrdParam(i)%yc)))cycle
        if(size(pAnyGrdParam(i)%xc) /= nc)cycle
        ifNew = .false.   
        do j = 1, nc
          if(.not.((pAnyGrdParam(i)%xc(j) .eps. lon(j)).and.(pAnyGrdParam(i)%yc(j) .eps. lat(j))))then
            ifNew = .true.
            exit 
          endif
        enddo
        if(.not. ifNew)then
          index = i
          anyGrdParamRefCount(i) = anyGrdParamRefCount(i) + 1
          return
        endif
      endif
    end do

    call set_error('No free any grid structures left','fu_next_free_ag_param')

  end function fu_ag_param_index 

  
   ! ***************************************************************

  integer function fu_new_ag_param_index() result(index)
    !
    ! Finds the first free structure and returns its index
    !
    implicit none

    ! Local variables
    integer :: i

    do i=1, max_nbr_any_grids
      if(anyGrdParamRefCount(i) == int_missing)then
        index = i
        pAnyGrdParam(i)%defined = silja_true
        anyGrdParamRefCount(i) = 1
        return
      endif
    end do
    call set_error('No free any grid structures left','fu_next_free_ag_param')

  end function fu_new_ag_param_index 

  !*************************************************************************

  subroutine release_grid(grid)
    ! Decrement the reference counter for this anygrid, and if zero, deallocate.
    ! Does nothing for other grids
    implicit none
    type(silja_grid), intent(in) :: grid
    
    integer :: ind_ag
    type(silam_any_grid_param), pointer :: agptr

    if (grid%gridtype /= anygrid ) return
    
    ind_ag = grid%ag%indParam

    if (anyGrdParamRefCount(ind_ag) < 1) then
      call set_error('reference count < 1 for anygrid', 'release_anygrid')
    end if
    
    agptr => pAnyGrdParam(ind_ag)
    if (fu_fails(associated(agptr), 'Anygrid not associated', 'release_anygrid')) return
    if (fu_fails(fu_true(agptr%defined), 'Anygrid not defined', 'release_anygrid')) return
    
    anyGrdParamRefCount(ind_ag) = anyGrdParamRefCount(ind_ag) - 1
    if (anyGrdParamRefCount(ind_ag) == 0) then
      call msg('Deallocating anygrid, ind', ind_ag)
      if (allocated(agptr%xc)) deallocate(agptr%xc)
      if (allocated(agptr%yc)) deallocate(agptr%yc)
      if (allocated(agptr%dx)) deallocate(agptr%dx)
      if (allocated(agptr%dy)) deallocate(agptr%dy)
      if (allocated(agptr%sin_map_rot)) deallocate(agptr%sin_map_rot)
      if (allocated(agptr%cos_map_rot)) deallocate(agptr%cos_map_rot)
      agptr%defined = silja_false
    end if
    
  end subroutine release_grid

  !************************************************************************************
 
  FUNCTION fu_set_any_grid(name, nx, ny, lon, lat) result(grid)
    !
    ! Sets a horizontal grid inside given area havig given gridsize.
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_grid) :: grid
    !
    ! Imported parameters with intent(in):
    CHARACTER (LEN=*), INTENT(in) :: name
    real, dimension(:), intent(inout), optional :: lon, lat
    integer :: nx, ny

    ! Local declarations:
    logical :: ifnew, ifData
    integer :: i
    !
    ! Have to ensure that lons and lats follow silam standard (-180:180 and -90:90) 
    !
    if (nx < 5 .or. ny<5) call ooops("smallanygrid")
    ifData = .true. 
    if(present(lon))then
      do i = 1, nx*ny
        ! Maybe 0 to 360
        if(lon(i)>180.)lon(i)= lon(i)-360.
        ! if still abnormal, error
        if(lon(i)>180. .or. lon(i)<-180.)then
          call set_error('strange lon values given', 'fu_set_any_grid')
          return
        endif
      enddo
    else
      ifData = .false.
    endif
    if(present(lat))then
      if(.not. ifData)then
        call set_error('Only latfield given','fu_set_any_grid')
        return
      endif
      do i = 1, nx*ny
      if(lat(i)>90. .or. lat(i)<-90.)then
        call set_error('strange lat values given', 'fu_set_any_grid')
        return
      endif
    enddo
    else
      if(ifData)then
        call set_error('Only lonfield given','fu_set_any_grid')
        return
      endif
    endif

    grid = grid_missing ! inital value 
    grid%name = name

    grid%ag%nx = nx
    grid%ag%ny = ny
    
    if(ifData)then
      grid%ag%indParam = fu_ag_param_index(nx*ny, lon, lat, ifNew)
      if(error)return
      if(ifnew)then
          allocate(pAnyGrdParam(grid%ag%indParam)%xc(nx*ny))
          pAnyGrdParam(grid%ag%indParam)%xc(1:nx*ny) = lon(1:nx*ny)
          allocate(pAnyGrdParam(grid%ag%indParam)%yc(nx*ny))
          pAnyGrdParam(grid%ag%indParam)%yc(1:nx*ny) = lat(1:nx*ny)
          pAnyGrdParam(grid%ag%indParam)%defined = silja_true
      endif
    else
!      call msg_warning('Setting anygrid without lons or lats given')      
      grid%ag%indParam = int_missing
    endif
    grid%gridtype = anygrid

  END FUNCTION fu_set_any_grid

  
  
  !*****************************************************************
 
  subroutine setAgLonlatFlds(grid, lon, lat, ifnew)
    
    !
    ! Has to be used, if anygrid was set without lats and lons
    !

    IMPLICIT NONE
    TYPE(silja_grid), intent(inout) :: grid
    real, dimension(:), pointer  :: lon, lat
    logical, intent(out) :: ifNew
    integer :: nc


    if(grid%gridtype /= anygrid)then
      call set_error('Only for grids of anygrid type','setAnygridParam')
      return
    endif

    nc = grid%ag%nx*grid%ag%ny

    grid%ag%indParam = fu_ag_param_index(nc, lon, lat, ifNew)
    if(error)return
    if(ifnew)then
      allocate(pAnyGrdParam(grid%ag%indParam)%xc(nc))
      pAnyGrdParam(grid%ag%indParam)%xc(1:nc) = lon(1:nc)
      allocate(pAnyGrdParam(grid%ag%indParam)%yc(nc))
      pAnyGrdParam(grid%ag%indParam)%yc(1:nc) = lat(1:nc)
      pAnyGrdParam(grid%ag%indParam)%defined = silja_true
    endif

  end subroutine setAgLonlatFlds
 
 
  !*****************************************************************

  subroutine setAnygridParam(grid, gVarNm, values)

    ! If values given, will copy the given grid variable values to grid parameter field,
    ! else will try to compute from lat and lon fields.
    ! If exists, will be overwritten

    IMPLICIT NONE

    TYPE(silja_grid), intent(in) :: grid
    CHARACTER (LEN=*), intent(in) :: gVarNm
    real, dimension(:), intent(in) :: values

    !local variables
    real, dimension(:), allocatable :: arrTmp
    integer :: nc, iParam
    character (len=*), parameter :: sub_name = "setAnygridParam"
    
    if(grid%gridtype /= anygrid)then
      call set_error('Only for grids of anygrid type', sub_name)
      return
    endif

    nc = grid%ag%nx*grid%ag%ny
    iParam = grid%ag%indParam

    allocate(arrTmp(nc))
    arrTmp(1:nc) = values(1:nc)

    select case(gVarNm)
      case("lon")
        ! Maybe 0 to 360
        where(arrTmp(1:nc)>180.) arrTmp= arrTmp-360.
        if (any(abs(arrTmp)>180.)) then ! if still abnormal, error
          call set_error('strange lon values given', sub_name)
          return
        endif
        call move_alloc(arrTmp, pAnyGrdParam(iParam)%xc)
      case("lat")
        if (any(abs(arrTmp)>90.)) then
          call set_error('strange lat values given', sub_name)
          return
        endif
        call move_alloc(arrTmp, pAnyGrdParam(iParam)%yc)
      case("dx")
          call move_alloc(arrTmp, pAnyGrdParam(iParam)%dx)
      case("dy")
          call move_alloc(arrTmp, pAnyGrdParam(iParam)%dy)
      case("sin_map_rot")
          call move_alloc(arrTmp, pAnyGrdParam(iParam)%sin_map_rot)
      case("cos_map_rot")
          call move_alloc(arrTmp, pAnyGrdParam(iParam)%cos_map_rot)
      case("x3dc")
          call move_alloc(arrTmp, pAnyGrdParam(iParam)%x3dc)
      case("y3dc")
          call move_alloc(arrTmp, pAnyGrdParam(iParam)%y3dc)
      case("z3dc")
          call move_alloc(arrTmp, pAnyGrdParam(iParam)%z3dc)

      case default
        call set_error(fu_connect_strings('Unknown parameter name', gVarNm), sub_name)
    end select

  end subroutine setAnygridParam


  subroutine completeAnygrid3DParam(grid)

    !
    !  Fills 3D cartesian coordinates from xc and yc assuming that 
    !  thay are already set, but x3dc, xy3dc and zx3dc are not
    !
    IMPLICIT NONE

    TYPE(silja_grid), intent(in) :: grid
    integer :: nc, iGrid, nx, ny

    real, dimension(:), pointer :: x, y, fPtr
    real, dimension(:,:), pointer :: x3d, y3d, z3d, fPtr2

    real, dimension(:), allocatable :: arrTmp
    character (len=*), parameter :: sub_name = "completeAnygrid3DParam"

    if(grid%gridtype /= anygrid)then
      call set_error('Only for grids of anygrid type','setAnygridParam')
      return
    endif

    iGrid = grid%ag%indParam
    if (allocated( pAnyGrdParam(iGrid)%x3dc)) then
      call set_error(" pAnyGrdParam(iGrid)%x3dc already allocated", sub_name)
      return
    endif

    
    nx = grid%ag%nx
    ny = grid%ag%ny
    nc = nx*ny
    x(1:nc) =>  pAnyGrdParam(iGrid)%xc(1:nc)
    y(1:nc) =>  pAnyGrdParam(iGrid)%yc(1:nc)

    allocate(pAnyGrdParam(iGrid)%x3dc(nc))
    fPtr => pAnyGrdParam(iGrid)%x3dc(1:nc)
    fPtr(1:nc) = cos(x*degrees_to_radians)*cos(y*degrees_to_radians)


    allocate(pAnyGrdParam(iGrid)%y3dc(nc))
    fPtr => pAnyGrdParam(iGrid)%y3dc(1:nc)
    fPtr(1:nc) = sin(x*degrees_to_radians)*cos(y*degrees_to_radians)


    allocate(pAnyGrdParam(iGrid)%z3dc(nc))
    fPtr => pAnyGrdParam(iGrid)%z3dc(1:nc)
    fPtr(1:nc) = sin(y*degrees_to_radians)

    !
    !  Cell sizes
    ! 
    if (.not. allocated( pAnyGrdParam(iGrid)%dx) .and. &
       & .not. allocated( pAnyGrdParam(iGrid)%dy)) then
        
        !!! Just make those
        x3d(1:nx,1:ny) => pAnyGrdParam(iGrid)%x3dc(1:nc)
        y3d(1:nx,1:ny) => pAnyGrdParam(iGrid)%y3dc(1:nc)
        z3d(1:nx,1:ny) => pAnyGrdParam(iGrid)%z3dc(1:nc)
        allocate(pAnyGrdParam(iGrid)%dx(nc), pAnyGrdParam(iGrid)%dy(nc)) 
        !! dx
        fPtr2(1:nx,1:ny) => pAnyGrdParam(iGrid)%dx(1:nc)
        fPtr2(2:nx-1,2:ny-1) = 0.5 * sqrt( (x3d(3:nx,2:ny-1) - x3d(1:nx-2,2:ny-1)) ** 2 + &
                                         & (y3d(3:nx,2:ny-1) - y3d(1:nx-2,2:ny-1)) ** 2 + &
                                         & (z3d(3:nx,2:ny-1) - z3d(1:nx-2,2:ny-1)) ** 2 )
        fPtr2(1,2:ny-1) = fPtr2(2,2:ny-1)
        fPtr2(nx,2:ny-1) = fPtr2(nx-1,2:ny-1)
        fPtr2(1:nx,1) = fPtr2(1:nx,2)
        fPtr2(1:nx,ny) = fPtr2(1:nx,ny-1)
        fPtr2(1:nx,1:ny) = fPtr2(1:nx,1:ny) * earth_radius
        !! And dy
        fPtr2(1:nx,1:ny) => pAnyGrdParam(iGrid)%dy(1:nc)
        fPtr2(2:nx-1,2:ny-1) = 0.5 * sqrt( (x3d(2:nx-1,3:ny) - x3d(2:nx-1,1:ny-2)) ** 2 + &
                                         & (y3d(2:nx-1,3:ny) - y3d(2:nx-1,1:ny-2)) ** 2 + &
                                         & (z3d(2:nx-1,3:ny) - z3d(2:nx-1,1:ny-2)) ** 2 )
        fPtr2(1,2:ny-1) = fPtr2(2,2:ny-1)
        fPtr2(nx,2:ny-1) = fPtr2(nx-1,2:ny-1)
        fPtr2(1:nx,1) = fPtr2(1:nx,2)
        fPtr2(1:nx,ny) = fPtr2(1:nx,ny-1)
        fPtr2(1:nx,1:ny) = fPtr2(1:nx,1:ny) * earth_radius

        call msg("Generated dx and dy sizes for the grid:")
        call report(grid)

    elseif (.not. allocated( pAnyGrdParam(iGrid)%dx) .or. &
       & .not. allocated( pAnyGrdParam(iGrid)%dy)) then 
        call set_error("Partly initialised xy cell sizes", sub_name)
    endif


  end subroutine completeAnygrid3DParam

  !*****************************************************************


  function fu_boundary_grid(grid, bName) result(bGrid)

    type(silja_grid), intent(in) :: grid
    character :: bName
    type(silja_grid) :: bGrid
    real :: cornerY, cornerX, dx, dy, fTmp, fTmp1
    integer :: nX, nY, i
    logical :: corner_in_geo_lonlat
    real, dimension(:), pointer :: lats, lons, dx_ptr, dy_ptr, sin_map_rot, cos_map_rot
    type(silam_any_grid_param), pointer :: gr_param_ptr
                                  
    if (bName == 'T' .or. bName == 'B')then
        bGrid = grid
    else
      select case(grid%gridtype )
      case(lonlat)
         call lonlat_grid_parameters(grid,&
                                  & cornerX, cornerY, corner_in_geo_lonlat, &
                                  & nX, nY, &
                                  & fTmp, fTmp1, &
                                  & dx, dy)
         select case(bName)
           case('N')
             cornerY = cornerY + (nY - 0.5) * dy
             nY=1
           case('S')
             cornerY = cornerY - 0.5 * dy
             nY=1
           case('E')
             cornerX = cornerX + (nX - 0.5) * dx
             nX = 1
           case('W')
             cornerX = cornerX - 0.5 * dx
             nX = 1
           case default
             call set_error('Strange boundary name','fu_boundary_grid')
         end select
 
         bGrid = fu_set_lonlat_grid (bName,  cornerX, cornerY,  corner_in_geo_lonlat, &
                                 & nx, ny, fu_pole(grid), dx, dy)                                      
       case(anygrid)
         
         lons => fu_work_array() 
         lats => fu_work_array() 
         dx_ptr => fu_work_array() 
         dy_ptr => fu_work_array() 
         sin_map_rot => fu_work_array()
         cos_map_rot => fu_work_array()
         call grid_dimensions(grid, nx, ny)
         gr_param_ptr => pAnyGrdParam(grid%ag%indParam)

         select case(bName)

           case('N')
             do i = 1, nx
               lats(i) = gr_param_ptr%yc(nx*(ny-1)+i) + 0.5*fu_dy_m_to_deg(gr_param_ptr%dy(nx*(ny-1)+i)) 
               lons(i) = gr_param_ptr%xc(nx*(ny-1)+i)
               dx_ptr(i) = gr_param_ptr%dx(nx*(ny-1)+i)
               dy_ptr(i) = gr_param_ptr%dy(nx*(ny-1)+i)
               cos_map_rot = gr_param_ptr%cos_map_rot(nx*(ny-1)+i)
               sin_map_rot = gr_param_ptr%sin_map_rot(nx*(ny-1)+i)
             enddo
             nY=1
           case('S')
             do i = 1, nx
               lats(i) = gr_param_ptr%yc(i) - 0.5*fu_dy_m_to_deg(gr_param_ptr%dy(i)) 
               lons(i) = gr_param_ptr%xc(i)
               dx_ptr(i) = gr_param_ptr%dx(i)
               dy_ptr(i) = gr_param_ptr%dy(i)
               cos_map_rot(i) = gr_param_ptr%cos_map_rot(i)
               sin_map_rot(i) = gr_param_ptr%sin_map_rot(i)
             enddo
             nY=1
           case('E')
             do i = 1, ny
               lats(i) = gr_param_ptr%yc(i*nx)
               lons(i) = gr_param_ptr%xc(i*nx) + 0.5*fu_dx_m_to_deg(gr_param_ptr%dx(i*nx), gr_param_ptr%yc(i*nx))
               dx_ptr(i) = gr_param_ptr%dx(i*nx)
               dy_ptr(i) = gr_param_ptr%dy(i*nx)
               cos_map_rot(i) = gr_param_ptr%cos_map_rot(i*nx)
               sin_map_rot(i) = gr_param_ptr%sin_map_rot(i*nx)
             enddo
             nX = 1
           case('W')
             do i = 1, ny
               lats(i) = gr_param_ptr%yc((i-1)*nx+1)
               lons(i) = gr_param_ptr%xc((i-1)*nx+1) - 0.5*fu_dx_m_to_deg(gr_param_ptr%dx((i-1)*nx+1), gr_param_ptr%yc((i-1)*nx+1))
               dx_ptr(i) = gr_param_ptr%dx((i-1)*nx+1)
               dy_ptr(i) = gr_param_ptr%dy((i-1)*nx+1)
               cos_map_rot(i) = gr_param_ptr%cos_map_rot((i-1)*nx+1)
               sin_map_rot(i) = gr_param_ptr%sin_map_rot((i-1)*nx+1)
             enddo
             nX = 1
           case default
             call set_error('Strange boundary name','fu_boundary_grid')
         end select        
         bGrid = fu_set_any_grid(bName, nx, ny, lons, lats)
         call setAnygridParam(bGrid, "dx", dx_ptr)
         call setAnygridParam(bGrid, "dy", dy_ptr)
         call setAnygridParam(bGrid, "cos_map_rot", cos_map_rot)
         call setAnygridParam(bGrid, "sin_map_rot", sin_map_rot)

         call free_work_array(lons)
         call free_work_array(lats)
         call free_work_array(dx_ptr)
         call free_work_array(dy_ptr)
         call free_work_array(cos_map_rot)
         call free_work_array(sin_map_rot)


       case default
         call set_error('Wrong gridtype','fu_boundary_grid')
       end select
    endif

  end function fu_boundary_grid

  !*****************************************************************

  function fu_area_from_grid(grid) result(area)
    !
    ! Creates an area from the given grid. A problem is that fu_set_lonlat_area
    ! operates with absolute values of latitude and longitude, requiring
    ! extra information whther they are south/north latitude or west/east longitude.
    ! Grid always has northern latitude and eastern longitudes positive, the opposite
    ! ones - negative. So, this ambiguity has to be coped in the call.
    !
    implicit none

    ! Return value of the function
    type(silam_area) :: area

    ! Imported variables with intent IN
    type(silja_grid), intent(in) :: grid


    real :: ymax, ymin, xmax, xmin, y
    logical :: ifPole, ifCross180, ifCross0
    integer :: i

    if(.not.defined(grid))then
      area = area_missing
      return
    end if

    select case (grid%gridtype)
      case (lonlat)
        area= fu_set_area(grid%name, &
                & grid%lonlat%sw_corner_modlat, north_flag, & 
                & grid%lonlat%sw_corner_modlat+grid%lonlat%dy_deg*(grid%lonlat%ny-1), north_flag,&
                & grid%lonlat%sw_corner_modlon, east_flag, &
                & grid%lonlat%sw_corner_modlon+grid%lonlat%dx_deg*(grid%lonlat%nx-1), east_flag, &
                & grid%lonlat%pole)
      
      case(anygrid)
        ! As we have no idea about projection, the output area will be in geographical one,
        ! which makes this a quite strange thing to do in the first place - the area might be huge
        call msg_warning('Taking area of anygrid in geo coordinates - might be huge','fu_area_from_grid')
        ! Actually the area borders should always be on grid boundaries, unless pole is covered
        ! Out of lazyness, take maxval and minval of whole array anyway ..
        
        ! Latitudes are simple - no periodic boundary there
        ymax = maxval(pAnyGrdParam(grid%ag%indParam)%yc)
        ymin = minval(pAnyGrdParam(grid%ag%indParam)%yc)
        
        ! longitudes
        ! Check if any geo pole is located inside the grid, if so, have to cover all around
        ifPole = .false.
        if (fu_point_inside_grid(grid, 0., -90., pole_geographical))then
          ifPole = .true.
          ymin = -90.
          xmin = -180.
          xmax = 180. - fu_dx_cell_deg(grid, grid%ag%nx/2, grid%ag%ny/2)
        endif
        if (fu_point_inside_grid(grid, 0., 90., pole_geographical))then
          ifPole = .true.
          ymax = 90.
          xmin = -180.
          xmax = 180. - fu_dx_cell_deg(grid, grid%ag%nx/2, grid%ag%ny/2)
        endif
      
        if(.not. ifPole)then
          ! here have to consider that world is round - area start and end are not necessarily the max
          ! and min lon values in the grid
          ! Check if grid crosses the 180 meridian
          ifCross180 = .false.
          ifCross0 = .false.
          y = ymin
          do while(y <= ymax)
            if (fu_point_inside_grid(grid, 180., y, pole_geographical))ifcross180 = .true.
            if (fu_point_inside_grid(grid, 0., y, pole_geographical))ifcross0 = .true.
            y = y + fu_dy_x_y_cell_deg(grid, grid%ag%nx/2, grid%ag%ny/2)
          enddo 
        
          if(.not. ifCross180)then  ! simple
            xmax = maxval(pAnyGrdParam(grid%ag%indParam)%xc)
            xmin = minval(pAnyGrdParam(grid%ag%indParam)%xc)
          elseif(.not. ifCross0)then ! from smallest positive to largest negative
            ! cycle the boundaries
            xmin = 180.
            xmax = -180.
            do i = 1, grid%ag%nx
              if(pAnyGrdParam(grid%ag%indParam)%xc(i) > 0) &
                & xmin = min(xmin, pAnyGrdParam(grid%ag%indParam)%xc(i))
              if(pAnyGrdParam(grid%ag%indParam)%xc(grid%ag%nx*(grid%ag%ny-1)+i) > 0) &
                & xmin = min(xmin, pAnyGrdParam(grid%ag%indParam)%xc(grid%ag%nx*(grid%ag%ny-1)+i))     
              if(pAnyGrdParam(grid%ag%indParam)%xc(i) < 0) &
                & xmax = max(xmax, pAnyGrdParam(grid%ag%indParam)%xc(i))
              if(pAnyGrdParam(grid%ag%indParam)%xc(grid%ag%nx*(grid%ag%ny-1)+i) < 0) &
                & xmax = max(xmax, pAnyGrdParam(grid%ag%indParam)%xc(grid%ag%nx*(grid%ag%ny-1)+i)) 
            enddo
            do i = 1, grid%ag%ny
              if(pAnyGrdParam(grid%ag%indParam)%xc(grid%ag%nx*i) > 0) &
                & xmin = min(xmin, pAnyGrdParam(grid%ag%indParam)%xc(grid%ag%nx*i))
              if(pAnyGrdParam(grid%ag%indParam)%xc(grid%ag%nx*(i-1)+1) > 0) &
                & xmin = min(xmin, pAnyGrdParam(grid%ag%indParam)%xc(grid%ag%nx*(i-1)+1))     
              if(pAnyGrdParam(grid%ag%indParam)%xc(grid%ag%nx*i) < 0) &
                & xmax = max(xmax, pAnyGrdParam(grid%ag%indParam)%xc(grid%ag%nx*i))
              if(pAnyGrdParam(grid%ag%indParam)%xc(grid%ag%nx*(i-1)+1) < 0) &
                & xmax = max(xmax, pAnyGrdParam(grid%ag%indParam)%xc(grid%ag%nx*(i-1)+1)) 
            enddo
          else ! too much trouble, take the whole world, most of it has to be covered anyway
            xmin = -180.
            xmax = 180. - fu_dx_cell_deg(grid, grid%ag%nx/2, grid%ag%ny/2)
          endif
        endif
        area= fu_set_area(grid%name, ymin, north_flag, ymax, north_flag,&
                        & xmin, east_flag,  xmax, east_flag, pole_geographical)                  
      case default
        call set_error('Only lonlat and anygrids so far','fu_area_from_grid')
        return
    end select
    return
  end function fu_area_from_grid



  ! **************************************************************
  ! **************************************************************
  !
  !      Information about the grid
  !
  ! **************************************************************
  ! **************************************************************


  ! ***************************************************************


  LOGICAL FUNCTION fu_grid_defined(grid)
    !
    ! Description:
    ! Returns a true value, if the grid has been given a value using
    ! one of the the functions fu_set_xxx_grid.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silja_grid), INTENT(in) :: grid

    SELECT CASE (grid%gridtype)

      CASE (lonlat, polarstereographic, anygrid)

      fu_grid_defined = .true.

    CASE default 

      fu_grid_defined = .false.

    END SELECT

  END FUNCTION fu_grid_defined


  ! ***************************************************************


  SUBROUTINE print_grid_report(grid)

    ! Description:
    ! Prints a report to screen describing  the grid.
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    TYPE(silja_grid), INTENT(in) :: grid
    integer, dimension(2) :: fUnits
    integer :: iUnit
    
    logical :: speak
    
    fUnits(1:2) = (/run_log_funit, 6/)


    do iUnit = 1,2
      if (smpi_global_rank > 0 .and. iUnit > 1) exit
      SELECT CASE (grid%gridtype)

        CASE (lonlat) 

          WRITE (fUnits(iUnit), '(2A)') '  ------- Grid report of lonlat-grid--', grid%name


        IF (grid%lonlat%pole == pole_geographical) THEN
            WRITE(fUnits(iUnit), fmt = '(A)') '  Grid in geographical lon-lat'
            WRITE(fUnits(iUnit), fmt = '(A, 2(F15.7,1x), 2(I12,1x))') &
                 & '  SW-corner (lon,lat), number of points x,y: ', &
                 & grid%lonlat%sw_corner_modlon,&
                 & grid%lonlat%sw_corner_modlat, &
                 & grid%lonlat%nx, grid%lonlat%ny

        ELSE
          WRITE(fUnits(iUnit), fmt = '(A, 2(F15.7,1x), 2(I12,1x))') &
              & '  modified SW-corner(lon,lat),number of points x,y: ', &
              & grid%lonlat%sw_corner_modlon,&
              & grid%lonlat%sw_corner_modlat, &
              & grid%lonlat%nx, grid%lonlat%ny

          WRITE(fUnits(iUnit), fmt = '(A, 2(F15.7,1x))') &
              & '  Grid south pole lon(E) and lat (N):  ',&
              & fu_southpole_lon(grid%lonlat%pole),&
              & fu_southpole_lat(grid%lonlat%pole)
        END IF

        WRITE(fUnits(iUnit), fmt = '(A, 2(F15.7,1x))')'  Grid_distance x and y degrees: ',&
                                   & grid%lonlat%dx_deg, grid%lonlat%dy_deg

        case(anygrid)

        WRITE (fUnits(iUnit), '(A,I0,A,A,A)') '  ------- Grid report of anygrid (', &
                                &grid%ag%indParam, '): "', trim(grid%name), '"'
        WRITE(fUnits(iUnit), fmt = '(A)') '  Grid in geographical lon-lat'
        WRITE(fUnits(iUnit), fmt = '(A, 2(F15.7,1x), 2(I12,1x))') &
            & '  SW-corner (lon,lat), number of points x,y: ', &
            & pAnyGrdParam(grid%ag%indParam)%xc(1) ,&
            & pAnyGrdParam(grid%ag%indParam)%yc(1), &
            & grid%ag%nx, grid%ag%ny

        if (allocated(pAnyGrdParam(grid%ag%indParam)%dx)) then
          WRITE(fUnits(iUnit), fmt = '(A, 2(F15.7,1x))')  &
                      &' Grid_distance x and y in meters for first gridcell: ',&
                                     & pAnyGrdParam(grid%ag%indParam)%dx(1) , &
                                     & pAnyGrdParam(grid%ag%indParam)%dy(1)
        else
          WRITE(fUnits(iUnit), fmt = '(A, 2(F15.7,1x))')  &
                      &' Grid_distance x and y in meters for first gridcell: NOT DEFINED'
        endif


        CASE default
          WRITE (fUnits(iUnit), '(A)')  '  Grid report: unknown grid, type=', grid%gridtype
      END SELECT

      WRITE (fUnits(iUnit), '(A)') ' ------------ End-of-grid-report ------------'
    enddo

  END SUBROUTINE print_grid_report


  ! ***************************************************************


  INTEGER FUNCTION fu_number_of_gridpoints(grid)

    ! Description:
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silja_grid), INTENT(in) :: grid

    SELECT CASE (grid%gridtype)

      CASE (lonlat)

      fu_number_of_gridpoints = grid%lonlat%nx * grid%lonlat%ny

      CASE (polarstereographic)

      fu_number_of_gridpoints = grid%ps%nx * grid%ps%ny

      CASE (anygrid )

      fu_number_of_gridpoints = grid%ag%nx * grid%ag%ny

    CASE default

      CALL set_error('unknown gridtype', 'fu_number_of_gridpoints')

    END SELECT

  END FUNCTION fu_number_of_gridpoints


!  ! ***************************************************************


  SUBROUTINE lonlat_grid_parameters(grid,&
                                  & corner_lon_E, corner_lat_N,corner_in_geo_lonlat, &
                                  & number_of_points_x, number_of_points_y, &
                                  & southpole_lon_E, southpole_lat_N, &
                                  & dx_deg, dy_deg)

    ! Returns necessary parameters from a geographical or
    ! rotated lat-lon grid. Used for weatherserver overcoats, where
    ! only reals and integres can be used.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silja_grid), INTENT(in) :: grid

    ! Imported parameters with intent OUT:
    REAL, INTENT(out) :: corner_lat_N,corner_lon_E,southpole_lat_N,southpole_lon_E, dx_deg, dy_deg

    INTEGER, INTENT(out) :: number_of_points_x, number_of_points_y

    LOGICAL, INTENT(out) :: corner_in_geo_lonlat

    ! Local declarations:
    TYPE(silja_lonlat_grid) :: lalo

    IF (.NOT.defined(grid)) THEN
      CALL set_error('undefined grid given','lonlat_grid_parameters')
      RETURN
    END IF

    IF (grid%gridtype == lonlat) THEN
      lalo = grid%lonlat
    ELSE
      CALL report(grid)
      CALL set_error('other than lonlat grid given', 'lonlat_grid_parameters')
      RETURN
    END IF


    corner_lat_N = lalo%sw_corner_modlat
    corner_lon_E = lalo%sw_corner_modlon

    corner_in_geo_lonlat = .false.

    southpole_lat_N = fu_southpole_lat(lalo%pole)
    southpole_lon_E = fu_southpole_lon(lalo%pole)

    dx_deg = lalo%dx_deg
    dy_deg = lalo%dy_deg

    number_of_points_x = lalo%nx
    number_of_points_y = lalo%ny

  END SUBROUTINE lonlat_grid_parameters


  ! ***************************************************************


  SUBROUTINE grid_dimensions(grid, nx, ny)

    ! Description:
    ! Returns the dimensions of a grid in x and y. 
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silja_grid), INTENT(in) :: grid

    ! Imported parameters with intent OUT:
    INTEGER, INTENT(out) :: nx, ny

    SELECT CASE (grid%gridtype)

      CASE (lonlat)
      nx = grid%lonlat%nx
      ny = grid%lonlat%ny

      CASE (polarstereographic)
      nx = grid%ps%nx
      ny = grid%ps%ny

      CASE (anygrid)
      nx = grid%ag%nx
      ny = grid%ag%ny

    CASE default

      CALL set_error('unknown gridtype', 'grid_dimensions')

    END SELECT

  END SUBROUTINE grid_dimensions


  !******************************************************************
  
  subroutine adjust_grid_dimensions(grid, nxNew, nyNew)
    !
    ! Adjusts the grid dimensions to the given ones by changing its resolution.
    ! However, if the change is small in one of the dimensions, it is left intact
    !
    IMPLICIT NONE

    ! Imported parameters:
    TYPE(silja_grid), INTENT(inout) :: grid
    INTEGER, INTENT(in) :: nxNew, nyNew

    ! Local variables
    integer :: nx, ny

    if(nxNew <= 0 .or. nyNew <= 0)then
      call msg('Strange new grid dimensions given:',nxNew, nyNew)
      call set_error('Strange new grid dimensions given','adjust_grid_dimensions')
    endif
    !
    ! The actual work depends on the grid type
    !
    SELECT CASE (grid%gridtype)

      CASE (lonlat)
        if(grid%lonlat%nx < 0.9 * real(nxNew) .or. grid%lonlat%nx > 1.1 * real(nxNew))then
          grid%lonlat%dx_deg = grid%lonlat%dx_deg * (real(grid%lonlat%nx) / real(nxNew))
          grid%lonlat%nx = nxNew
        endif
        if(grid%lonlat%ny < 0.9 * real(nyNew) .or. grid%lonlat%ny > 1.1 * real(nyNew))then
          grid%lonlat%dy_deg = grid%lonlat%dy_deg * (real(grid%lonlat%ny) / real(nyNew))
          grid%lonlat%ny = nyNew
        endif

      CASE (polarstereographic)
        if(grid%ps%nx * grid%ps%ny < 0.8 * real(nxNew * nyNew) .or. &
         & grid%ps%nx * grid%ps%ny > 1.2 * real(nxNew * nyNew))then
          grid%ps%grid_dist_m = grid%ps%grid_dist_m * &
                              & sqrt(real(grid%ps%nx * grid%ps%ny) / real(nxNew * nyNew))
          grid%ps%grid_dist_m = nxNew * nyNew
        endif

      CASE (anygrid)
        call set_error('Dimensions of anygrid cannot be adjusted','adjust_grid_dimensions')
        call unset_error('adjust_grid_dimensions')

      CASE default
        CALL set_error('unknown gridtype', 'adjust_grid_dimensions')

    END SELECT

  end subroutine adjust_grid_dimensions


  ! ***************************************************************


  FUNCTION fu_pole_of_grid(grid)
    IMPLICIT NONE  
    ! The return value of this function:
    TYPE(silam_pole) :: fu_pole_of_grid
    ! Imported parameters with intent IN:
    TYPE(silja_grid), INTENT(in) :: grid

    IF (grid%gridtype == lonlat) THEN
      fu_pole_of_grid = grid%lonlat%pole
    elseif(grid%gridtype == anygrid)then
      fu_pole_of_grid = pole_missing
    ELSE
      CALL set_error('pole only avalilable for a lonlat-grid','fu_pole_of_grid')      
    END IF

  END FUNCTION fu_pole_of_grid


  ! ***************************************************************


  INTEGER FUNCTION fu_gridtype(grid)
    IMPLICIT NONE
    TYPE(silja_grid), INTENT(in) :: grid
    fu_gridtype = grid%gridtype
  END FUNCTION fu_gridtype


  ! ***************************************************************


  REAL FUNCTION fu_lon_native_from_grid(x, y, grid)
    !
    ! Returns geographical co-ordinates of the point (x,y) defined
    ! in relative grid co-ordinates. No pole transformation !
    !
    ! Author: Mikhail Sofiev, FMI
    !
    IMPLICIT NONE
    !
    ! Imported variables with intent IN:
    REAL, INTENT(in) :: x,y ! Relative co-ordinates in the grid
    TYPE(silja_grid), INTENT(in) :: grid

    SELECT CASE(grid%gridtype)
    CASE (lonlat)
      fu_lon_native_from_grid = (x-1.) * grid%lonlat%dx_deg + grid%lonlat%sw_corner_modlon
    case (anygrid)   ! return the geographical coordinate
      fu_lon_native_from_grid = fu_2d_interpolation(pAnyGrdParam(grid%ag%indParam)%xc, &
                                                  & x, y, grid%ag%nx, grid%ag%ny, linear, &
                                                  & notAllowed)
    CASE DEFAULT
      CALL set_error('Sorry, only lonlat grids','fu_lon_native_from_grid')
    END SELECT

  END FUNCTION fu_lon_native_from_grid


  ! ***************************************************************


  REAL FUNCTION fu_lat_native_from_grid(x, y, grid)
  !
  ! Returns geographical co-ordinates of the point (x,y) defined
  ! in relative grid co-ordinates. No pole transformation !
  !
  ! Author: Mikhail Sofiev, FMI
  !
  IMPLICIT NONE
  !
  ! Imported variables with intent IN:
  REAL, INTENT(in) :: x,y ! Relative co-ordinates in the grid
  TYPE(silja_grid), INTENT(in) :: grid

  SELECT CASE(fu_gridtype(grid))
  CASE (lonlat)
    fu_lat_native_from_grid = (y-1.) * grid%lonlat%dy_deg + grid%lonlat%sw_corner_modlat
  case (anygrid) ! return the geographical coordinate
     fu_lat_native_from_grid = fu_2d_interpolation(pAnyGrdParam(grid%ag%indParam)%yc, &
                                                 & x, y, grid%ag%nx, grid%ag%ny, linear, notAllowed)   
  CASE DEFAULT
    CALL set_error('Sorry, only lonlat grids','fu_lat_native_from_grid')
  END SELECT

  END FUNCTION fu_lat_native_from_grid


  ! ***************************************************************


  REAL FUNCTION fu_lon_geographical_from_grid(x, y, grid)
    !
    ! Returns geographical co-ordinates of the point (x,y) defined
    ! in relative grid co-ordinates. No pole transformation !
    !
    ! Author: Mikhail Sofiev, FMI
    !
    IMPLICIT NONE
    !
    ! Imported variables with intent IN:
    REAL, INTENT(in) :: x,y ! Relative co-ordinates in the grid
    TYPE(silja_grid), INTENT(in) :: grid

    ! Local variables
    real :: lon, lat

    SELECT CASE(grid%gridtype)
    CASE (lonlat)
      CALL modify_lonlat((y-1.) * grid%lonlat%dy_deg + grid%lonlat%sw_corner_modlat, &  ! latitude
                       & (x-1.) * grid%lonlat%dx_deg + grid%lonlat%sw_corner_modlon, &
                       & grid%lonlat%pole, pole_geographical, lat, lon)
      fu_lon_geographical_from_grid = lon
   
    case (anygrid)      
      fu_lon_geographical_from_grid =  fu_2d_interpolation( &
                                                     & pAnyGrdParam(grid%ag%indParam)%xc, &
                                                     & x, y, grid%ag%nx, grid%ag%ny, linear, &
                                                     & notAllowed) 
    CASE DEFAULT
      CALL set_error('Sorry, only lonlat grids','fu_lon_geographical_from_grid')
    END SELECT

  END FUNCTION fu_lon_geographical_from_grid


  ! ***************************************************************


  REAL FUNCTION fu_lat_geographical_from_grid(x, y, grid)
    !
    ! Returns geographical co-ordinates of the point (x,y) defined
    ! in relative grid co-ordinates. No pole transformation !
    !
    ! Author: Mikhail Sofiev, FMI
    !
    IMPLICIT NONE
    !
    ! Imported variables with intent IN:
    REAL, INTENT(in) :: x,y ! Relative co-ordinates in the grid
    TYPE(silja_grid), INTENT(in) :: grid

    ! Local variables
    real :: lon, lat

    SELECT CASE(grid%gridtype)
    CASE (lonlat)
      CALL modify_lonlat((y-1.) * grid%lonlat%dy_deg + grid%lonlat%sw_corner_modlat, &  ! latitude
                       & (x-1.) * grid%lonlat%dx_deg + grid%lonlat%sw_corner_modlon, &
                       & grid%lonlat%pole, pole_geographical, lat, lon)
      fu_lat_geographical_from_grid = lat
   
    case (anygrid)      
      fu_lat_geographical_from_grid =  fu_2d_interpolation(pAnyGrdParam(grid%ag%indParam)%yc, &
                                                               & x, y, grid%ag%nx, grid%ag%ny, linear, notAllowed)   
    CASE DEFAULT
      CALL set_error('Sorry, only lonlat grids','fu_lon_geographical_from_grid')
    END SELECT

  END FUNCTION fu_lat_geographical_from_grid
  
  !***************************************************************

  function fu_geolats_fld(grid)
  
    type(silja_grid), intent(in) :: grid
    real, dimension(:), pointer :: fu_geolats_fld
    integer :: ix, iy
    real :: lon, lat
    
    SELECT CASE(grid%gridtype)
    CASE (lonlat)
      fu_geolats_fld => fu_work_array()

      do ix = 1, grid%lonlat%nx
        do iy = 1, grid%lonlat%ny
          CALL modify_lonlat((iy-1.) * grid%lonlat%dy_deg + grid%lonlat%sw_corner_modlat, &  ! latitude
                       & (ix-1.) * grid%lonlat%dx_deg + grid%lonlat%sw_corner_modlon, &
                       & grid%lonlat%pole, pole_geographical, lat, lon)
          fu_geolats_fld(ix+grid%lonlat%nx*(iy-1)) = lat
        enddo
      enddo
    case (anygrid)
      fu_geolats_fld => pAnyGrdParam(grid%ag%indParam)%yc
    CASE DEFAULT
      CALL set_error('Sorry, only lonlat grids','fu_lon_geographical_from_grid')
    END SELECT
    
  end function fu_geolats_fld
  
  !***************************************************************

  function fu_geolons_fld(grid)
  
    type(silja_grid), intent(in) :: grid
    real, dimension(:), pointer :: fu_geolons_fld
    integer :: ix, iy
    real :: lon, lat
    
    SELECT CASE(grid%gridtype)
    CASE (lonlat)
      fu_geolons_fld => fu_work_array()

      do ix = 1, grid%lonlat%nx
        do iy = 1, grid%lonlat%ny
          CALL modify_lonlat((iy-1.) * grid%lonlat%dy_deg + grid%lonlat%sw_corner_modlat, &  ! latitude
                       & (ix-1.) * grid%lonlat%dx_deg + grid%lonlat%sw_corner_modlon, &
                       & grid%lonlat%pole, pole_geographical, lat, lon)
          fu_geolons_fld(ix+grid%lonlat%nx*(iy-1)) = lon
        enddo
      enddo
    case (anygrid)
      fu_geolons_fld => pAnyGrdParam(grid%ag%indParam)%xc
    CASE DEFAULT
      CALL set_error('Sorry, only lonlat grids','fu_lon_geographical_from_grid')
    END SELECT
    
  end function fu_geolons_fld

  !***************************************************************

  function fu_dx_fld_m(grid)
  
    type(silja_grid), intent(in) :: grid
    real, dimension(:), pointer :: fu_dx_fld_m
    integer :: ix, iy
    
    SELECT CASE(grid%gridtype)
    CASE (lonlat)
      fu_dx_fld_m => fu_work_array()

      do ix = 1, grid%lonlat%nx
        do iy = 1, grid%lonlat%ny
          fu_dx_fld_m(ix+grid%lonlat%nx*(iy-1)) = &
                         & fu_dx_deg_to_m(grid%lonlat%dx_deg, &
                                       & (real(iy)-1.)*grid%lonlat%dy_deg + &
                                               & grid%lonlat%sw_corner_modlat) 
        enddo
      enddo
    case (anygrid)
      fu_dx_fld_m => pAnyGrdParam(grid%ag%indParam)%dx
    CASE DEFAULT
      CALL set_error('Sorry, only lonlat grids','fu_lon_geographical_from_grid')
    END SELECT
    
  end function fu_dx_fld_m

  !***************************************************************
 
   function fu_dy_fld_m(grid)
  
    type(silja_grid), intent(in) :: grid
    real, dimension(:), pointer :: fu_dy_fld_m
    integer :: ix, iy
    
    SELECT CASE(grid%gridtype)
    CASE (lonlat)
      fu_dy_fld_m => fu_work_array()

      do ix = 1, grid%lonlat%nx
        do iy = 1, grid%lonlat%ny
          fu_dy_fld_m(ix+grid%lonlat%nx*(iy-1)) = fu_dy_deg_to_m(grid%lonlat%dy_deg)
        enddo
      enddo
    case (anygrid)
      fu_dy_fld_m => pAnyGrdParam(grid%ag%indParam)%dy
    CASE DEFAULT
      CALL set_error('Sorry, only lonlat grids','fu_lon_geographical_from_grid')
    END SELECT
    
  end function fu_dy_fld_m

  !***************************************************************
 
  function fu_sin_map_rot_fld(grid)
  
    type(silja_grid), intent(in) :: grid
    real, dimension(:), pointer :: fu_sin_map_rot_fld
    
    SELECT CASE(grid%gridtype)
    CASE (lonlat)
       call set_error('not implemented for lonlat grid','fu_sin_map_rot_fld')
    case (anygrid)
      fu_sin_map_rot_fld => pAnyGrdParam(grid%ag%indParam)%sin_map_rot
    CASE DEFAULT
      CALL set_error('Sorry, only any-grids','fu_sin_map_rot_fld')
    END SELECT
    
  end function fu_sin_map_rot_fld

  !***************************************************************
 
   function fu_cos_map_rot_fld(grid)
  
    type(silja_grid), intent(in) :: grid
    real, dimension(:), pointer :: fu_cos_map_rot_fld
    
    SELECT CASE(grid%gridtype)
    CASE (lonlat)
      call set_error('not implemented for lonlat grid','fu_cos_map_rot_fld')
    case (anygrid)
      fu_cos_map_rot_fld => pAnyGrdParam(grid%ag%indParam)%cos_map_rot
    CASE DEFAULT
      CALL set_error('Sorry, only any-grids','fu_cos_map_rot_fld')
    END SELECT
    
  end function fu_cos_map_rot_fld

  !***************************************************************


  real function fu_cell_size_midpoint(grid)
    !
    ! Computes the size of the first cell of the grid. Unit = [m2]
    !
    implicit none

    ! Imported parameters with intent IN
    type(silja_grid), intent(in) :: grid

    if(.not.defined(grid))then
      call set_error('Undefined grid','fu_cell_size_midpoint')
      return
    end if

    select case (fu_gridtype(grid))
      case(lonlat)
        fu_cell_size_midpoint = fu_dx_deg_to_m(grid%lonlat%dx_deg, &
                                             & grid%lonlat%sw_corner_modlat + grid%lonlat%dy_deg * int(grid%lonlat%ny * 0.5)) * &
                              & fu_dy_deg_to_m(grid%lonlat%dy_deg)
      case(anygrid)
        fu_cell_size_midpoint = pAnyGrdParam(grid%ag%indParam)%dx(nint(0.5*grid%ag%ny-1)*grid%ag%nx+nint(0.5*grid%ag%nx)) * &
                              & pAnyGrdParam(grid%ag%indParam)%dy(nint(0.5*grid%ag%ny-1)*grid%ag%nx+nint(0.5*grid%ag%nx))

      case default
        call set_error('Unsupported grid type','fu_cell_size_midpoint')
        return
    end select
  end function fu_cell_size_midpoint


  !***************************************************************

  real function fu_cell_size_n_th_cell(grid, iCell)
    !
    ! Computes the size of the cell iCell of the grid. Unit = [m2]
    ! The iCell is the co-ordinate in memory - from 1-D array. 
    ! Appropriate tranformation to 2-D is, of course, due.
    !
    implicit none

    ! Imported parameters with intent IN
    type(silja_grid), intent(in) :: grid
    integer, intent(in) :: iCell

    ! Local variables
    real :: fY

    if(.not.defined(grid))then
      call set_error('Undefined grid','fu_cell_size_n_th_cell')
      return
    end if

    select case (fu_gridtype(grid))
      case(lonlat)

        fY = int(iCell / grid%lonlat%nx) + 1
        
        fu_cell_size_n_th_cell = fu_dx_deg_to_m(grid%lonlat%dx_deg, &
                                              & (fY-1.)*grid%lonlat%dy_deg + &
                                              & grid%lonlat%sw_corner_modlat) * &
                               & fu_dy_deg_to_m(grid%lonlat%dy_deg)
      case(anygrid)
        fu_cell_size_n_th_cell = pAnyGrdParam(grid%ag%indParam)%dx(iCell) &
             & * pAnyGrdParam(grid%ag%indParam)%dy(iCell)


      case default
        call set_error('Unsupported grid type','fu_cell_size_n_th_cell')
        return
    end select
  end function fu_cell_size_n_th_cell


  !***************************************************************

  real function fu_cell_size_x_y_cell(grid, ix, iy)
    !
    ! Computes the size of the cell iCell of the grid. Unit = [m2]
    ! The iCell is the co-ordinate in memory - from 1-D array. 
    ! Appropriate tranformation to 2-D is, of course, due.
    !
    implicit none

    ! Imported parameters with intent IN
    type(silja_grid), intent(in) :: grid
    integer, intent(in) :: ix, iy

    if(.not.defined(grid))then
      call set_error('Undefined grid','fu_cell_size_x_y_cell')
      return
    end if

    select case (fu_gridtype(grid))
      case(lonlat)

        fu_cell_size_x_y_cell = fu_dx_deg_to_m(grid%lonlat%dx_deg, &
                                             & (real(iy)-1.)*grid%lonlat%dy_deg + &
                                             & grid%lonlat%sw_corner_modlat) * &
                              & fu_dy_deg_to_m(grid%lonlat%dy_deg)
      case(anygrid)
        fu_cell_size_x_y_cell = pAnyGrdParam(grid%ag%indParam)%dx(grid%ag%nx*(iy-1)+ix) * &
                              & pAnyGrdParam(grid%ag%indParam)%dy(grid%ag%nx*(iy-1)+ix)

      case default
        call set_error('Unsupported grid type','fu_cell_size_x_y_cell')
        return
    end select
  end function fu_cell_size_x_y_cell


  !****************************************************************

  real function fu_dx_x_y_cell_m(grid, ix, iy)
    !
    ! Computes the size of the cell iCell of the grid. Unit = [m2]
    ! The iCell is the co-ordinate in memory - from 1-D array. 
    ! Appropriate tranformation to 2-D is, of course, due.
    !
    implicit none

    ! Imported parameters with intent IN
    type(silja_grid), intent(in) :: grid
    integer, intent(in) :: ix, iy

    if(.not.defined(grid))then
      call set_error('Undefined grid','fu_dx_x_y_cell_m')
      return
    end if

    select case (fu_gridtype(grid))
      case(lonlat)

        fu_dx_x_y_cell_m = max(0., fu_dx_deg_to_m(grid%lonlat%dx_deg, &
                                             & (real(iy)-1.)*grid%lonlat%dy_deg + &
                                             & grid%lonlat%sw_corner_modlat) )
      case(anygrid)
        fu_dx_x_y_cell_m = pAnyGrdParam(grid%ag%indParam)%dx(grid%ag%nx*(iy-1)+ix) 

      case default
        call set_error('Unsupported grid type','fu_dx_x_y_cell_m')
        return
    end select
  end function fu_dx_x_y_cell_m

  !***************************************************************

  real function fu_dy_x_y_cell_m(grid, ix, iy)
    !
    ! Computes the size of the cell iCell of the grid. Unit = [m2]
    ! The iCell is the co-ordinate in memory - from 1-D array. 
    ! Appropriate tranformation to 2-D is, of course, due.
    !
    implicit none

    ! Imported parameters with intent IN
    type(silja_grid), intent(in) :: grid
    integer, intent(in) :: ix, iy

    if(.not.defined(grid))then
      call set_error('Undefined grid','fu_dy_x_y_cell_m')
      return
    end if

    select case (fu_gridtype(grid))
      case(lonlat)

        fu_dy_x_y_cell_m = fu_dy_deg_to_m(grid%lonlat%dy_deg)
      case(anygrid)
        fu_dy_x_y_cell_m = pAnyGrdParam(grid%ag%indParam)%dy(grid%ag%nx*(iy-1)+ix)

      case default
        call set_error('Unsupported grid type','fu_dy_x_y_cell_m')
        return
    end select
  end function fu_dy_x_y_cell_m
 
 !****************************************************************

  
  real function fu_dx_x_y_cell_deg(grid, ix, iy)
      implicit none

    ! Imported parameters with intent IN
    type(silja_grid), intent(in) :: grid
    integer, intent(in) :: ix, iy

    if(.not.defined(grid))then
      call set_error('Undefined grid','fu_dx_x_y_cell_deg')
      return
    end if

    select case (fu_gridtype(grid))
      case(lonlat)
        fu_dx_x_y_cell_deg = grid%lonlat%dx_deg
      case(anygrid)
        fu_dx_x_y_cell_deg = fu_dx_m_to_deg(pAnyGrdParam(grid%ag%indParam)%dx(grid%ag%nx*(iy-1)+ix), &
                                    &  pAnyGrdParam(grid%ag%indParam)%yc(grid%ag%nx*(iy-1)+ix))

      case default
        call set_error('Unsupported grid type','fu_dx_x_y_cell_deg')
        return
    end select
  end function fu_dx_x_y_cell_deg
  
  !****************************************************************
  
  real function fu_dy_x_y_cell_deg(grid, ix, iy)
      implicit none

    ! Imported parameters with intent IN
    type(silja_grid), intent(in) :: grid
    integer, intent(in) :: ix, iy

    if(.not.defined(grid))then
      call set_error('Undefined grid','fu_dy_x_y_cell_deg')
      return
    end if

    select case (fu_gridtype(grid))
      case(lonlat)
        fu_dy_x_y_cell_deg = grid%lonlat%dy_deg
      case(anygrid)
        fu_dy_x_y_cell_deg = fu_dy_m_to_deg(pAnyGrdParam(grid%ag%indParam)%dy(grid%ag%nx*(iy-1)+ix))

      case default
        call set_error('Unsupported grid type','fu_dy_x_y_cell_deg')
        return
    end select
  end function fu_dy_x_y_cell_deg

  !****************************************************************


  logical function fu_ifLonGlobal(grid)
    !
    ! Checks whether the given grid covers the whole 360 degrees along the longitude
    !
    implicit none

    ! Imported parameter
    type(silja_grid), intent(in) :: grid

    integer :: nx, ny, iy
    real :: dist, dist_max

    if(fu_fails(defined(grid), 'Undefined grid','fu_ifLonGlobal'))return

    select case (fu_gridtype(grid))
      case(lonlat)

        fu_ifLonGlobal = (grid%lonlat%dx_deg * (real(grid%lonlat%nx) + 0.1)) > 360.0 

      case(anygrid)   
      ! Global runs should be possible here too, but to ensure that the periodic boundaries 
      ! are applied correctly, grid must be reasonable 
      ! just check the lons of first and last gridcell - if reasonably close to each other  
      ! definition of reasonable - less than 2 dx in lon and lat
      ! have to scan all latitudes
        call grid_dimensions(grid, nx, ny)
        fu_ifLonGlobal = .true.
        do iy = 1, ny
          dist_max = fu_gc_distance(pAnyGrdParam(grid%ag%indParam)%xc(nx*iy),   &
                                  & pAnyGrdParam(grid%ag%indParam)%xc(nx*iy-1), &
                                  & pAnyGrdParam(grid%ag%indParam)%yc(nx*iy),   &
                                  & pAnyGrdParam(grid%ag%indParam)%yc(nx*iy-1)) + &
                   & fu_gc_distance(pAnyGrdParam(grid%ag%indParam)%xc(nx*(iy-1)+1), &
                                  & pAnyGrdParam(grid%ag%indParam)%xc(nx*(iy-1)+2), &
                                  & pAnyGrdParam(grid%ag%indParam)%yc(nx*(iy-1)+1), &
                                  & pAnyGrdParam(grid%ag%indParam)%yc(nx*(iy-1)+2))
         
          dist = fu_gc_distance(pAnyGrdParam(grid%ag%indParam)%xc(nx*iy), &
                              & pAnyGrdParam(grid%ag%indParam)%xc(nx*(iy-1)+1), &
                              & pAnyGrdParam(grid%ag%indParam)%yc(nx*iy), &
                              & pAnyGrdParam(grid%ag%indParam)%yc(nx*(iy-1)+1))
                              
          if(dist > dist_max)then
            fu_ifLonGlobal = .false.
            exit
          endif
        enddo

      case default
        call msg('Unsupported grid type:',fu_gridtype(grid))
        call set_error('Unsupported grid type','fu_ifLonGlobal')
        return
    end select

  end function fu_ifLonGlobal

  
  !****************************************************************

  logical function fu_ifLatGlobal(grid)
    !
    ! Checks whether the given grid covers the whole 180 degrees along
    ! the latitude -- this is checked with the tolerance lat_global_threshold
    !
    implicit none

    ! Imported parameter
    type(silja_grid), intent(in) :: grid

!    fu_ifLatGlobal =  (fu_ifNorthPole(grid) .and. fu_ifSouthPole(grid))
    fu_ifLatGlobal =  ( fu_ifPoleIncluded(grid, southern_boundary) .and. & 
                      & fu_ifPoleIncluded(grid, northern_boundary))

  end function fu_ifLatGlobal


  !****************************************************************

  subroutine dy_to_pole(grid, iBoundary, frac, dyDeg) 

    ! Calculates max distance of upper or lower mesh line to a grid pole. 
    ! Returns:
    !  frac -- distance in units of dy (averaged over the line)
    !  dyDeg -- dy (averaged over the line) in degrees
    implicit none

    ! Imported parameter
    type(silja_grid), intent(in) :: grid
    integer, intent(in) :: iBoundary
    real, intent(out) :: frac, dyDeg
    real, dimension(:), pointer :: xC, yC
    integer :: nx, ny, ix
    real :: fTmp, dy
        
    ! Some stupidty checks
    if (fu_fails(defined(grid), 'Undefined grid','fu_degrees_to_pole'))return
    if (all(iBoundary /= (/northern_boundary, southern_boundary/) )) then
       call msg("Strange pole  number", iBoundary)
       call set_error("Can make only southern and northern poles","fu_degrees_to_pole")
       return
    endif

    select case (fu_gridtype(grid))
      case(lonlat)
        dyDeg = grid%lonlat%dy_deg
        ! Distance from last point towards the pole (negative if behind the pole)
        if (iBoundary == northern_boundary) then ! Last grid point
           fTmp = 90. - (grid%lonlat%sw_corner_modlat + grid%lonlat%dy_deg * (grid%lonlat%ny-1))
        else
           fTmp = 90. + grid%lonlat%sw_corner_modlat
        endif
        frac =  fTmp / dyDeg
      case(anygrid)
        !
        ! First or last row
        !
        call grid_dimensions(grid, nx, ny)
        if (iBoundary == northern_boundary) then
           dy = sum(pAnyGrdParam(grid%ag%indParam)%dy(nx*(ny-1)+1:nx*ny))/nx
           xC => pAnyGrdParam(grid%ag%indParam)%xc(nx*(ny-1)+1:nx*ny)
           yC => pAnyGrdParam(grid%ag%indParam)%yc(nx*(ny-1)+1:nx*ny)
        else
           dy = sum(pAnyGrdParam(grid%ag%indParam)%dy(1:nx))/nx
           xC => pAnyGrdParam(grid%ag%indParam)%xc(1:nx)
           yC => pAnyGrdParam(grid%ag%indParam)%yc(1:nx)
        endif
      
        do ix = 2, nx
           fTmp = max(fTmp,fu_gc_distance(xC(ix), xC(1), yC(iX), yC(1)))
        enddo
        dyDeg =  fu_dy_m_to_deg(dy)
        frac =  fTmp*0.5 / dy
      case default
        call msg('Unsupported grid type:',fu_gridtype(grid))
        call set_error('Unsupported grid type','fu_ifNorthPole')
        return
    end select
   end subroutine dy_to_pole

  !****************************************************************

   logical function fu_ifPolarCapPossible(grid, iBoundary)
    !
    ! Checks if Polar cap can be done for the grid
    implicit none

    ! Imported parameter
    type(silja_grid), intent(in) :: grid
    integer, intent(in) :: iBoundary
    real  :: frac, dyDeg
        
    fu_ifPolarCapPossible = .false.
    if(fu_fails(defined(grid), 'Undefined grid','fu_ifNorthPolePossible'))return

    ! grid can only include pole if global in longitude direction ..
    if(.not. fu_ifLonGlobal(grid))return
    call dy_to_pole(grid, iBoundary, frac, dyDeg)

    if(dyDeg*frac > lat_global_threshold )then
       !Too big pole -- no cap
       return 
    else
       if(frac < 0.5)then
           call msg(boundary_name(iBoundary)//"polar cap radius in cells, dy in degrees", &
               & frac, dyDeg)
           call set_error(boundary_name(iBoundary)//' polar cap radius is smaller than half cell',&
                                 & 'fu_ifPolarCapPossible')
           return
       endif
       fu_ifPolarCapPossible = .true.
    endif

   end function fu_ifPolarCapPossible

  !****************************************************************

   logical function fu_ifPoleIncluded(grid, iBoundary)
    !
    ! Checks if the grid has a hole in the pole:
    ! i.e. field in the pole can be obtained with linear interpolation.
    ! For that it must be lon-global and the size of a 
    ! hole in the mesh must be not bigger than than 
    ! dy
    implicit none

    ! Imported parameter
    type(silja_grid), intent(in) :: grid
    integer, intent(in) :: iBoundary
    real  :: frac, dyDeg
        
    fu_ifPoleIncluded = .false.
    if(fu_fails(defined(grid), 'Undefined grid','fu_ifPoleIncluded'))return

    ! grid can only include pole if global in longitude direction ..
    if(.not. fu_ifLonGlobal(grid))return
    call dy_to_pole(grid, iBoundary, frac, dyDeg)

    fu_ifPoleIncluded =  (frac <= 0.5)

   end function fu_ifPoleIncluded




  !****************************************************************

  logical function fu_stdSilamGrid(grid)
    !
    ! Checks whether the given grid satisfies the SILAM standards
    !
    implicit none

    ! Imported parameter
    type(silja_grid), intent(in) :: grid

    ! Local variables
    logical :: ifLonOK, ifLatOK

    fu_stdSilamGrid = .false.
    ifLonOK = .false.
    ifLatOK = .false.

    if(.not.defined(grid))then
      call set_error('Undefined grid','fu_stdSilamGrid')
      return
    end if

    ! SILAM grid must have corners belonging to interval (-180,-90) - (180,90)
    ! Alternatively, it can cover more than a globe with starting point with less than 
    ! one grid cell from the (-180,-90) point.

    select case (fu_gridtype(grid))
      case(lonlat)
        !
        ! Check the globality of the grid. Global ones can eventually run over more than a globe. All
        ! what we require is to have starting point at -180,-90.
        !
        if(grid%lonlat%dx_deg * grid%lonlat%nx > 360.0 - 0.1 * grid%lonlat%dx_deg)then
          !
          ! The grid is global or over-the-globe along the longitude. In this case, the starting point
          ! must be close to -180 lon
          !
          ifLonOK = grid%lonlat%sw_corner_modlon + 0.501*grid%lonlat%dx_deg >= -180.0 .and. &
                  & grid%lonlat%sw_corner_modlon - grid%lonlat%dx_deg <= -180.0
        endif

        if(abs(grid%lonlat%dy_deg * grid%lonlat%ny - 180.0) / grid%lonlat%dy_deg < 0.1 .or. &
         & grid%lonlat%dy_deg * grid%lonlat%ny > 180.0)then
          !
          ! The grid is global or over-the-globe along the latitude. In this case, the starting point
          ! must be close to -90 lat
          !
          ifLatOK = grid%lonlat%sw_corner_modlat >= -90.0001 .and. &
                  & grid%lonlat%sw_corner_modlat - grid%lonlat%dy_deg <= -90.0
        endif
        !
        ! If the grid is either non-global or does not satusfy the above conditions, the last chance
        ! is to check it for limited-area conditions
        !
        if(.not. ifLonOK) ifLonOK = grid%lonlat%sw_corner_modlon >= -180.0 .and. &
                                  & grid%lonlat%sw_corner_modlon + grid%lonlat%dx_deg * grid%lonlat%nx <= 180.0
        if(.not. ifLonOK) ifLonOK = grid%lonlat%sw_corner_modlon >= 0.0 .and. &
                                  & grid%lonlat%sw_corner_modlon + grid%lonlat%dx_deg * grid%lonlat%nx <= 360.0
        if(.not. ifLonOK)then
          call msg('longitude not OK, corner, distance:', grid%lonlat%sw_corner_modlon, grid%lonlat%dx_deg)
          return
        endif

!        fu_stdSilamGrid = ifLatOK .or. (grid%lonlat%sw_corner_modlat >= -90.0 .and. &
!                                & grid%lonlat%sw_corner_modlat + grid%lonlat%dy_deg * (grid%lonlat%ny-1) <= 90.0)
        fu_stdSilamGrid = ifLatOK .or. (grid%lonlat%sw_corner_modlat > -90.01 .and. &
                                & grid%lonlat%sw_corner_modlat + grid%lonlat%dy_deg * (grid%lonlat%ny-1) < 90.01)

!        fu_stdSilamGrid = grid%lonlat%sw_corner_modlon >= -180.0 .and. &
!                        & grid%lonlat%sw_corner_modlat >= -90.0 .and. &
!                        & grid%lonlat%sw_corner_modlon + grid%lonlat%dx_deg * grid%lonlat%nx <= 180.0 .and. &
!                        & grid%lonlat%sw_corner_modlat + grid%lonlat%dy_deg * (grid%lonlat%ny-1) <= 90.0
        if(.not. fu_stdSilamGrid)then
          call msg('latitude not OK, corner:', grid%lonlat%sw_corner_modlat)
          call msg('latitude not OK, distance:', grid%lonlat%dy_deg)
        endif
      
      case(anygrid)
        !
        ! Should be some checking here that the grid is more or less normal, whatever that means (regular, monotonous, etc) 
        !
    !    call msg_warning('all anygrids are considered standard for now', 'fu_stdSilamGrid')
        fu_stdSilamGrid = .true. 
      
      case default
        call msg('Unsupported grid type:',fu_gridtype(grid))
        call set_error('Unsupported grid type','fu_stdSilamGrid')
        return
    end select

  end function fu_stdSilamGrid


  !**********************************************************************

  subroutine make_grid_lon_global(grid, ifAdjustmentAllowed)
    !
    ! Creates the global grid along longitude following the given template.
    ! Note that the only allowed changes are the starting grid coordinates and
    ! number of grid cells. Resolution and position of existing cells must stay.
    ! Some adjustment may be explicitly allowed only for longitude resolution
    !
    implicit none

    ! Imported parameters
    type(silja_grid), intent(inout) :: grid
    logical, intent(in) :: ifAdjustmentAllowed

    ! Local variables
    integer :: iTmp
    real :: fTmp

    select case (fu_gridtype(grid))
      case(lonlat)
        !
        ! Check the principal possibility to handle the case: 360 must be divided into 
        ! integer number of cells!
        !
        iTmp = int(360. / grid%lonlat%dx_deg + 0.5)
        if(abs(iTmp * grid%lonlat%dx_deg - 360.) < 0.01 * grid%lonlat%dx_deg)then
          grid%lonlat%nx = iTmp
        else
          if(ifAdjustmentAllowed)then
            grid%lonlat%dx_deg = 360. / real(iTmp)
            grid%lonlat%nx = iTmp
          else
            call set_error('Error in globality is more than 1% but adjustments forbidden', &
                         & 'make_grid_lon_global')
            return
          endif
        endif
        !
        ! Force starting longitude to be within -180:180 - and close to it
        !
        if(grid%lonlat%sw_corner_modlon > 180.)then
          grid%lonlat%sw_corner_modlon = grid%lonlat%sw_corner_modlon - 360.
        elseif(grid%lonlat%sw_corner_modlon < -180.)then
          iTmp = int((-180. - grid%lonlat%sw_corner_modlon) / grid%lonlat%dx_deg) + 1
          grid%lonlat%sw_corner_modlon = grid%lonlat%sw_corner_modlon + iTmp * grid%lonlat%dx_deg
        elseif(grid%lonlat%sw_corner_modlon - grid%lonlat%dx_deg > -180.0)then
          iTmp = int((-180. - grid%lonlat%sw_corner_modlon) / grid%lonlat%dx_deg)
          grid%lonlat%sw_corner_modlon = grid%lonlat%sw_corner_modlon + iTmp * grid%lonlat%dx_deg
        endif
      case (anygrid)
        call set_error('anygrids cannot be enlarged','make_grid_lon_global')
      case default
        call set_error('Cannot handle grids other than lon-lat so far','make_grid_lon_global ')
        return
    end select    

  end subroutine make_grid_lon_global


  !**********************************************************************

  subroutine make_grid_lat_global(grid, ifAdjustmentAllowed)
    !
    ! Creates the global grid along longitude following the given template.
    ! Note that the only allowed changes are the starting grid coordinates and
    ! number of grid cells. Resolution and position of existing cells must stay.
    ! Some adjustment may be explicitly allowed only for longitude resolution
    !
    implicit none

    ! Imported parameters
    type(silja_grid), intent(inout) :: grid
    logical, intent(in) :: ifAdjustmentAllowed

    ! Local variables
    integer :: iTmp
    real :: fTmp

    select case (fu_gridtype(grid))
      case(lonlat)
        !
        ! Check the principal possibility to handle the case: 360 must be divided into 
        ! integer number of cells!
        !
        iTmp = int(180. / grid%lonlat%dy_deg + 0.5)
        if(abs(iTmp * grid%lonlat%dy_deg - 180.) < 0.01 * grid%lonlat%dy_deg)then
          grid%lonlat%ny = iTmp
        else
          if(ifAdjustmentAllowed)then
            grid%lonlat%dy_deg = 180. / real(iTmp)
            grid%lonlat%ny = iTmp
          else
            call set_error('Error in globality is more than 1% but adjustments forbidded', &
                         & 'make_grid_lat_global')
            return
          endif
        endif
        !
        ! Force starting longitude to be within -180:180
        !
        if(grid%lonlat%sw_corner_modlat > 90.)then
          grid%lonlat%sw_corner_modlat = grid%lonlat%sw_corner_modlat - 180.
        elseif(grid%lonlat%sw_corner_modlat < -90.)then
          iTmp = int((-90. - grid%lonlat%sw_corner_modlat) / grid%lonlat%dy_deg) + 1
          grid%lonlat%sw_corner_modlat = grid%lonlat%sw_corner_modlat + iTmp * grid%lonlat%dy_deg
        endif

      case (anygrid)
        call set_error('anygrids cannot be enlarged','make_grid_lat_global')

      case default
        call set_error('Cannot handle grids other than lon-lat so far','make_grid_lat_global ')
        return
    end select    

  end subroutine make_grid_lat_global


  !*******************************************************************************************
  
  subroutine make_global_lonlat_grid(nx, ny, dx, dy, grid_out)
    !
    ! Makes a global grid with given resolution or given number of grid cells.
    ! If both cell size and number of cells are given, cell size is ignored and recalculated
    !
    implicit none
    
    ! imported parameters
    integer, intent(in) :: nx, ny
    real, intent(in) :: dx, dy
    type(silja_grid), intent(out) :: grid_out
    
    ! Local variables
    real :: dx_local, dy_local, lon_start, lat_start
    integer :: nx_local, ny_local
    !
    ! What is given - number of cells or resolution?
    !
    if(nx == int_missing .or. nx <= 1)then
      !
      ! No number of cells: resolution must exist
      !
      if(dx <= 0 .or. dy <= 0. .or. (dx .eps. real_missing) .or. (dy .eps. real_missing))then
        call msg('missing dx or dy with undefined nx',dx,dy)
        call set_error('missing dx or dy with undefined nx','make_global_lonlat_grid')
        return
      endif
      !
      ! Having the resolution given, find the number of grid cells, which will cover the whole globe
      !
      nx_local = nint(360. / dx)
      
      ny_local = nint((180. - 2. * lat_global_threshold) / dy)
      do while((ny_local + 2) * dy >= 180. .and. ny_local > 0)  ! ensure enough distance from poles
        ny_local = ny_local - 1
      end do
      dy_local = dy
      
    else
      !
      ! Number of cells must exist
      !
      if(ny == int_missing .or. ny <= 1)then
        call msg('Missing ny with reasonable nx:',nx,ny)
        call set_error('Missing ny with reasonable nx','make_global_lonlat_grid')
        return
      endif
      !
      ! Having the number of cells given, calculate the resolution. dx is simple, dy to be made
      !
      ny_local = ny
      dy_local = (180. - 2. * lat_global_threshold) / (real(ny_local - 2)) ! half-cell into threshold circle
      if(dy_local * real(ny_local) > 180. - lat_global_threshold)then      ! Coarse grid, half-cell is too much
        dy_local = (180. - lat_global_threshold) / (real(ny_local - 1))    ! half of threshold zone covered
      endif
    endif  ! what is given?
        
    dx_local = 360. / real(nx_local)            ! to get rid of rounding error
    lon_start = -180. + dx_local / 2.           ! left edge of the first cell is at 0 meridian
    lat_start = 90. - real(ny_local - 1) * dy_local / 2.
    
    grid_out =  fu_set_lonlat_grid ('autoglobal_grid', &
                                  & lon_start, lat_start, &
                                  & .true., & !                   corner_in_geographical_lonlat, &
                                  & nx_local, ny_local, &
                                  & pole_geographical, & 
                                  & dx_local, dy_local)
  end subroutine make_global_lonlat_grid

  
  !******************************************************************************************
  
  subroutine ensure_standard_grid(gridTmp, ifNonStadardGrid, grid_data, pDataTmp)
    !
    ! Checks the grid for being standard - or tries to reposition it if it is a global non-standard
    !
    implicit none
    
    ! Imported parameters
    type(silja_grid), intent(inout) :: gridTmp
    logical, intent(out) :: ifNonStadardGrid
    real, dimension(:), target, optional, intent(in) :: grid_data
    real, dimension(:), optional, pointer :: pDataTmp
    
    ! Local variables
    integer :: iShiftX, iShiftY, iXoriginal, iYoriginal, index_out, ix, iy
    
    ! Check the grid
    !
    if(fu_stdSilamGrid(gridTmp))then
      ifNonStadardGrid = .false. ! grid is OK, the working data and input one are the same
      if(present(grid_data) .and. present(pDataTmp)) pDataTmp => grid_data
    else
      !
      ! Non-standard global grid. Check the directions of repositioning and then repick the 
      ! data from input set to work array.
      !
      if(fu_gridtype(gridTmp) /= lonlat)then
        call set_error('Only lonlat grids can be non-standard','ensure_standard_grid')
        return
      endif
      ifNonStadardGrid = .true.
      if(fu_ifLonGlobal(gridTmp)) then
        call reposition_global_grid_lon(gridTmp, iShiftX)
      else
        iShiftX = 0
      endif
      if(fu_ifLatGlobal(gridTmp)) then
        call reposition_global_grid_lat(gridTmp, iShiftY)
      else
        iShiftY = 0
      endif
      if(error)return
      if(iShiftX /= 0 .or. iShiftY /= 0)then
        !
        ! Grid is repositioned. Check that now it is OK and proceed. Otherwise set the error 
        !
        if(fu_stdSilamGrid(gridTmp))then
          !
          ! Reposition has made the grid standard. Re-arrange the data accodingly
          ! Method: go through the whole gridTmp (the grid after repositioning) and 
          ! fill-in each array element from the grid_original (before repositioning)
          !
          if(present(grid_data) .and. present(pDataTmp))then
            pDataTmp => fu_work_array(fu_number_of_gridpoints(gridTmp))
            index_out = 0
            do iy = 1, gridTmp%lonlat%ny
              do ix = 1, gridTmp%lonlat%nx
                index_out = index_out + 1  ! We scan the array in correct order
                if(iShiftX /= 0)then
                  iXoriginal = ix + iShiftX ! + 1
                  if(iXoriginal < 1)then    ! We cross the grid border jumping to the other side!!
                    iXoriginal = iXoriginal + gridTmp%lonlat%nx
                  elseif(iXoriginal > gridTmp%lonlat%nx)then
                    iXoriginal = iXoriginal - gridTmp%lonlat%nx
                  endif
                else
                  iXoriginal = ix 
                endif
                if(iShiftY /= 0)then
                  call msg_warning('shiftY is not zero. Reposition along y axis is forbidden', 'ensure_standard_grid')
                  call report(gridTmp)
                  call msg('Suggested reposition:', iShiftX, iShiftY)
                  call set_error('shiftY is not zero. Reposition along y axis is forbidden', 'ensure_standard_grid')
                  return
         !         iYoriginal = iy + iShiftY !+ 1
         !         if(iYoriginal < 1)then     ! We cross the grid border jumping to the other side!!
         !           iYoriginal = iYoriginal + gridTmp%lonlat%ny
         !         elseif(iYoriginal > gridTmp%lonlat%ny)then
         !           iYoriginal = iYoriginal - gridTmp%lonlat%ny
         !         endif
                else
                  iYoriginal = iy
                endif
                pDataTmp(index_out) = grid_data(iXoriginal + gridTmp%lonlat%nx * (iYoriginal - 1))
              end do
            end do
          endif  ! if data pointers present

        else
          !
          ! Still problematic grid. Stop
          !
          call set_error('Reposition does not help non-standard grid','ensure_standard_grid')
          call msg('Original grid after reposition:')
          call report(gridTmp)
          return
        endif
      else
        !
        ! Grid was not repositioned, so it is still non-standard. Set the error
        !
        call set_error('Cannot do anything with non-standard grid','ensure_standard_grid')
        call report(gridTmp)
        return
      endif  ! If reposition leads to reaonsble shifts

    endif  ! if non-standard input grid
    
  end subroutine ensure_standard_grid

  
  !****************************************************************
  !****************************************************************
  !
  !  Grid transformations, re-projections and other 
  !  grid-related actions
  !
  !****************************************************************
  !****************************************************************

  subroutine grid_data_hor_select_new_grid(grid_original, &
                                         & storage_grid,&
                                         & storage_area, &
                                         & ifAdjust, &  ! if adjust dimensions to meteo_grid
                                         & grid_new)
    ! Imported parameters
    TYPE(silam_area), INTENT(in) :: storage_area
    TYPE(silja_grid), INTENT(in) :: storage_grid, grid_original
    logical, intent(in) :: ifAdjust
    type(silja_grid), intent(out) :: grid_new
    
    ! Local parameters
    type(silja_grid) :: grid_original_standardised, grid_new_intermed
    type(silja_logical) :: ifAdjusted
    logical :: ifNonStandardGrid
    !
    ! Prior to anything we have to check that the input grid is the standard SILAM one.
    ! ECMWF has a nasty habit to use global grids with longitude varying in 0:360, while
    ! SILAM standard is -180,180. This poses evident problems, so we correct them here as they
    ! are corrected in io_server: simply renumber the whole grid and rearrange its data array.
    ! This is a painful decision but all other ways are much more dangerous because involve
    ! modifications of HIRLAM grid routines, which are not suitable for handling these grids.
    !
    grid_original_standardised = grid_original ! Must save the original grid from possible changes

    call ensure_standard_grid(grid_original_standardised, ifNonStandardGrid)

    if(error)return

    !
    ! Check storage-definitons
    !
    IF (.NOT.defined(storage_grid)) THEN

      IF (.NOT.defined(storage_area)) THEN

        grid_new = grid_original_standardised
!        if(ifNonStadardGrid) call free_work_array(arTmp)
!        RETURN  ! original grid is needed

      ELSE
        !
        ! 2. Data from storage_area wanted.
        !
        ! Calculate the the new grid, that contains whole area.
        ! Picking to new grid (no interpolation expected !!) is done in hila library.
        !
        grid_new = fu_grid_containing_area(grid_original_standardised, storage_area)

        if(ifAdjust)then
          grid_new_intermed = grid_new
          call adjust_grid_to_sample(grid_new, meteo_grid, ifAdjusted)
          IF(fu_true(ifAdjusted))THEN
            IF(.not.fu_if_grid_covered(grid_new, grid_original_standardised))THEN
              call msg('Before adjustment, grid_new')
              call report(grid_new_intermed)
              call msg('Before adjustment, grid_template')
              call report(storage_grid)
              call msg('After adjustment, grid_new, grid_template')
              call report(grid_new)
              CALL set_error('Adjusted grid bigger than the data area 1','grid_data_hor_select_new_grid')
              RETURN
            END IF
          END IF
        end if

        IF (error) RETURN

      END IF ! area_given

    ELSE ! grid given

      !--------------------------------------------------------
      !
      ! 3. Data in storage_grid wanted - interpolation expected.
      !    However, if grids correspond to each other we can avoid costly
      !    interpolation. Check first...

      !    ifAdjust may force exact match of the grids or allow Arakawa shift
      !    However, adjustment to sample must always be done - plus-minus one
      !    grid cell due to Arakawa shift must be taken into account
      !
      !    ATTENTION. We recently allowed the ECMWF-type global grids covering
      !    the range from 0:360 and -90:90 lon-lat respectively. This is NOT the SILAM
      !    standard grid, which must be -180:180, -90:90. All below tricks are allowed 
      !    only for standard SILAM grids. Check and skip the whole story if the 
      !    grid is not the standard one
      !
      if(ifAdjust)then
        !
        ! The storage_grid is forced whether the input one corresponds to it or not, no checking 
        !
        grid_new = storage_grid  ! whatever happens, we force the storage_grid
      else
          !!
          !! Is it really a good idea to allow arakawa grids here ?
          !! If one wants shifted grids, better to tell it explicitly from calling routines

          !! This error is triggered to detect and eliminate such calls (RK)
          call set_error("ifAdjust == .False. came to grid_data_hor_select_new_grid", &
                                  & 'grid_data_hor_select_new_grid')
          return 

        !
        ! If the arakawa-correspondence is allowed, we may use this option to avoid
        ! costly interpolation. However, first check for exact correspondence to avoid 
        ! numerical problems arising from grid->area->grid cycle transformation
        !
        if(fu_grids_correspond(grid_original_standardised, storage_grid)) then
          grid_new = storage_grid
        else
          if(fu_grids_arakawa_correspond(grid_original_standardised, storage_grid))then
              if(storage_grid%gridtype == anygrid)then
                  if(fu_if_grid_covered(storage_grid, grid_original_standardised))then
                    grid_new = grid_original_standardised
                  else
                    call set_error('Storage anygrid not covered with arakawa shifted one', &
                                     & 'grid_data_hor_select_new_grid')
                    call report(storage_grid)
                    call report(grid_original_standardised)            
                    return
                  endif
              else

                grid_new = fu_grid_containing_area(grid_original_standardised, &
                                             & fu_area_from_grid(storage_grid))
                grid_new_intermed = grid_new
                call adjust_grid_to_sample(grid_new, storage_grid, ifAdjusted)
                IF(fu_true(ifAdjusted))THEN
                  IF(.not.fu_if_grid_covered(grid_new, grid_original_standardised))THEN
                    call msg('grid_original_standardised')
                    call report(grid_original_standardised)
                    call msg('Before adjustment, grid_new')
                    call report(grid_new_intermed)
                    call msg('Storage grid:')
                    call report(storage_grid)
                    call msg('After adjustment, grid_new, grid_template')
                    call report(grid_new)
                    CALL set_error('Adjusted grid bigger than the data area 2', &
                                 & 'grid_data_hor_select_new_grid')
                    RETURN
                  END IF
                END IF
              endif
          else
            grid_new = storage_grid ! no Arakawa-correspondence, force the storage_grid
          endif  ! grids Arakawa-correspond
        endif  ! grid correspond
      end if
      IF (error) RETURN

    END IF ! grid_given => area ignored

!    call msg('Select 3')
!open (21,file='before_grid_transformation.grads',access='direct',form='unformatted',recl=fu_number_of_gridpoints(gridTmp))
!write(21,rec=1) (arTmp(ix),ix=1,fu_number_of_gridpoints(gridTmp))
!close(21)
!call msg('')
!call msg('before_grid_transformation')
!call report(gridTmp)
!call msg('after_grid_transformation')
!call report(grid_new)

  
  end subroutine grid_data_hor_select_new_grid

  
  !********************************************************************************
  
  SUBROUTINE grid_data_horizontal_select(grid_original, &
                                       & grid_data,&
                                       & selection_done, &
                                       & grid_new,&
                                       & selected_grid_data, &
                                       & ifRandomise, iAccuracy, method, &
                                       & fMissingValue, iOutside)
    ! 
    ! Since often only a part of the total area covered by the field is
    ! wanted to store in memory, this routine selects the right area
    ! and makes the necessary changes to the field-id:s grid. The
    ! desired area is defined either by storage_area or storage_grid,
    ! but not both!
    !
    ! ifAdjust has different meanings depending on the storage_grids and
    ! storage_area.
    ! - If storage_area is defined then ifAdjust determines the nx and ny of the
    !   selected grid - they must correspond to system values
    ! - If storage_grid is defined - ifAdjust determines whether the Arakawa-grids
    !   must be interpolated to the system one.
    !
    ! Method:
    ! fu_grid_containing_area
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    REAL, DIMENSION(:), INTENT(in), target :: grid_data
    TYPE(silja_grid), INTENT(in) :: grid_original, grid_new
    LOGICAL, INTENT(in) :: ifRandomise
    real, intent(in) :: fMissingValue
    integer, intent(in) :: iOutside, iAccuracy, method

    ! Imported parameters with intent out:
    LOGICAL, INTENT(out) :: selection_done
    REAL, DIMENSION(:), pointer :: selected_grid_data

    ! Local declarations:
    TYPE(silja_logical) :: ifAdjusted
    real, dimension(:), pointer :: arTmp
    logical :: ifNonStadardGrid, ifSubArea
    type(silja_grid) :: gridTmp
    type(THorizInterpStruct), pointer :: pHIS
    integer :: nxSrc, nySrc, nxOut, nyOut, ix, iy, i, kOut, kSrc, ix_shift, iy_shift
    integer, dimension(:), pointer :: arFlag


!    call msg('Select 1')
!open (21,file='start_of_sub.grads',access='direct',form='unformatted',recl=fu_number_of_gridpoints(grid_original))
!write(21,rec=1) (grid_data(ix),ix=1,fu_number_of_gridpoints(grid_original))
!close(21)
!call msg('start_of_sub')
!call report(grid_original)

    !
    ! Prior to anything we have to check that the input grid is the standard SILAM one.
    ! ECMWF has a nasty habit to use global grids with longitude varying in 0:360, while
    ! SILAM standard is -180,180. This poses evident problems, so we correct them here as they
    ! are corrected in io_server: simply renumber the whole grid and rearrange its data array.
    ! This is a painful decision but all other ways are much more dangerous because involve
    ! modifications of HIRLAM grid routines, which are not suitable for handling these grids.
    !
    gridTmp = grid_original ! Must save the original grid from possible changes

    call ensure_standard_grid(gridTmp, ifNonStadardGrid, grid_data, arTmp)
    selection_done = ifNonStadardGrid  ! if something fails or further selection is not needed

    if(error)return

    call grid_dimensions(gridTmp, nxSrc, nySrc)
    call grid_dimensions(grid_new, nxOut, nyOut)

    !
    ! If grids are the same, just copy the data - or do nothing if the initial grid is standard
    !
    if(gridTmp == grid_new)then
      if(ifNonStadardGrid)then
        selected_grid_data(1 : nxOut*nyOut) = arTmp(1 : nxOut*nyOut)
        call free_work_array(arTmp)
      else
        call free_work_array(selected_grid_data)
      endif
      return
    endif
    !
    ! If the new grid is a sub-area of the original one, just pick the needed rectangular
    !
    call SubArea_Chk(gridTmp, grid_new, ix_shift, iy_shift, ifSubArea)
    if(ifSubArea)then
      do iy = 1, nyOut
       do ix = 1, nxOut
         selected_grid_data(ix+(iy-1)*nxOut) = arTmp(ix+ix_shift+(iy+iy_shift-1)*nxSrc)
       end do
      end do
      if(ifNonStadardGrid) call free_work_array(arTmp)
      selection_done = .true.
      RETURN
    end if
    !
    ! Final step - make actual data picking or interpolating
    !
    !
    ! Get the instrument for the interpolation and metadata
    !

    pHIS => fu_horiz_interp_struct(gridTmp, grid_new, method, ifRandomise, iAccuracy, iOutside, ifCreate=.false.) 
    if (.not. associated(pHIS)) then
        pHIS => fu_horiz_interp_struct(gridTmp, grid_new, method, ifRandomise, iAccuracy, iOutside)
    endif

    selected_grid_data(1:nxOut*nyOut) = 0.0
    arFlag => fu_work_int_array(nxOut*nyOut)
    if(error)return
    arFlag(1:nxOut*nyOut) = 0
    !
    ! The interpolation itself
    !
    do iy = pHIS%iyStartTo, pHIS%iyEndTo
      do ix = pHIS%ixStartTo, pHIS%ixEndTo
        if (associated(pHIS%ifValid)) then
                if (.not. pHIS%ifValid(ix,iy)) cycle
        endif
        kOut = ix + (iy-1)*nxOut
        arFlag(kOut) = 1  !Zero is a valid weight
        do i = 1, pHIS%nCoefs
          if(pHIS%indX(i,ix,iy) == 0) exit  ! all indices are used up
          kSrc = pHIS%indX(i,ix,iy) + (pHIS%indY(i,ix,iy)-1)*nxSrc
          if(abs(arTmp(kSrc) - real_missing) > 1e-5 * abs(real_missing)) then  ! if source value is valid
            selected_grid_data(kOut) = selected_grid_data(kOut) + arTmp(kSrc) * pHIS%weight(i,ix,iy)
          endif
          arFlag(kOut) = 1
        end do
      end do
    end do
    !
    ! Those cells that have not received any information should be set to missing value
    !
    where(arFlag(1:nxOut*nyOut) == 0) selected_grid_data(1:nxOut*nyOut) = fMissingValue
    call free_work_array(arFlag)
                           
    selection_done = .not. error

    if(ifNonStadardGrid) call free_work_array(arTmp)

!    call msg('Select 5')




  !===================================================================================                             

  END SUBROUTINE grid_data_horizontal_select

  !********************************************************************************
  
  SUBROUTINE remap_field( grid_in, &
                        & data_in,&
                        & grid_out,&
                        & data_out, &
                        & ifRandomise, iAccuracy, method, iOutside, omp_locks)
    !  Implementation grid_data_horizontal_select without any excessive 
    !  creativity: get THorizInterpStruct and apply it with minimal check
    !  Should be thread-safe
    !
    use omp_lib_kinds, only : OMP_lock_kind

    IMPLICIT NONE


    ! Imported parameters with intent IN:
    REAL, DIMENSION(:), INTENT(in) ::  data_in
    REAL, DIMENSION(:), INTENT(out), target :: data_out
    TYPE(silja_grid), INTENT(in) :: grid_in, grid_out

    ! These guys are passsed to interp_struct
    LOGICAL, INTENT(in) :: ifRandomise
    integer, intent(in) :: iOutside, iAccuracy, method
    integer(kind = OMP_lock_kind), dimension(:), intent(inout) :: omp_locks 
          !!Lock for OMP use of fu_horiz_interp_struct 
          !in multithread environment must be initialized

    ! Local declarations:
    type(THorizInterpStruct), pointer :: pHIS
    REAL, DIMENSION(:,:), pointer :: pOut !!2Dpointer to the target data
    integer :: nxSrc, nySrc, nxOut, nyOut, ix, iy, i, kSrc
    integer :: ithread, nthreads


    call grid_dimensions(grid_in, nxSrc, nySrc)

    ithread = 0
    nthreads = 1
    !$ ithread  = omp_get_thread_num()
    !$ nthreads = omp_get_num_threads()
    !
    ! If grids are the same, just copy the data - or do nothing if the initial grid is standard
    !
    if(grid_in == grid_out)then
      data_out(1:nxSrc*nySrc) = data_in(1:nxSrc*nySrc)
      return
    endif
    call grid_dimensions(grid_out, nxOut, nyOut)

!$    if (nthreads ==  1) then
      !No locks needed
      pHIS => fu_horiz_interp_struct(grid_in, grid_out, method, ifRandomise, iAccuracy, iOutside, ifCreate=.true.) 
!$    else
!$       call  OMP_SET_LOCK(omp_locks(ithread+1))
!$      !Try to get existing struct
!$       pHIS => fu_horiz_interp_struct(grid_in, grid_out, method, ifRandomise, iAccuracy, iOutside, ifCreate=.false.) 
!$       call  OMP_UNSET_LOCK(omp_locks(ithread+1)) !Unlock thread
!$      if (.not. associated(pHIS)) then
!$         do i=1, nthreads
!$             call  OMP_SET_LOCK(omp_locks(i)) !Lock all
!$         enddo
!$         !! Create the struct
!$         pHIS => fu_horiz_interp_struct(grid_in, grid_out, method, ifRandomise, iAccuracy, iOutside, ifCreate=.true.)
!$         do i=1, nthreads
!$            call  OMP_UNSET_LOCK(omp_locks(i)) !Unlock all
!$         enddo
!$      endif
!$    endif
    if (error) return

    
    ! The interpolation itself
    !
    pOut(1:nxOut, 1:nyOut) => data_out(1:nxOut*nyOut) !

    !set missing whatever is not covered by interp_struct
    if (pHIS%iOutside == setZero) then
      pOut(:,1:pHIS%iyStartTo-1) = 0.
      pOut(:,pHIS%iyEndTo+1:nyOut) = 0.
      pOut(1:pHIS%ixStartTo-1,:) = 0.
      pOut(  pHIS%ixEndTo+1:nxOut,:) = 0.
    else
      pOut(:,1:pHIS%iyStartTo-1) = real_missing
      pOut(:,pHIS%iyEndTo+1:nyOut) = real_missing
      pOut(1:pHIS%ixStartTo-1,:) = real_missing
      pOut(  pHIS%ixEndTo+1:nxOut,:) = real_missing
    endif


    !Loop over Target grid
    do iy = pHIS%iyStartTo, pHIS%iyEndTo
      do ix = pHIS%ixStartTo, pHIS%ixEndTo
        !missing interp
        if (pHIS%weight(1,ix,iy)  == real_missing) then
          pOut(ix,iy) = real_missing
          continue
        endif

        pOut(ix,iy) = 0.
        do i = 1, pHIS%nCoefs
          if(pHIS%indX(i,ix,iy) == 0) exit  ! all indices are used up

          kSrc = pHIS%indX(i,ix,iy) + (pHIS%indY(i,ix,iy)-1)*nxSrc  !!Index in the source grid
          if (data_in(kSrc) == real_missing) then !Any missing _in_ causes missing _out_
            pOut(ix,iy) = real_missing
            exit
          else
             pOut(ix,iy) = pOut(ix,iy) + data_in(kSrc) * pHIS%weight(i,ix,iy)
          endif
        end do
      end do!ix out
    end do!iy out

  end SUBROUTINE remap_field

  ! **********************************************************************************

  FUNCTION fu_grid_containing_area(grid, area) result(grid_new)
    ! 
    ! Calculates the grid, that is the smallest possible subgrib of
    ! the given grid_orig, that contains the whole given area.
    ! The gridpoints in grid_new are all exact gridpoints of
    ! grid_orig, only area if different.
    !
    ! If the area is totally outside the given grid, than
    ! grid_missing is returned and a warning is displayed. 
    !
    ! If the area is partially covered by the given grid, but is
    ! outside in some dimension, then that border of the new grid is
    ! determined by the original grid, instead of area. So the
    ! resulting grid never contains points that aren't in the
    ! original grid.
    !
    ! If ifAdjust is true than an attempt of adjustment to the dimensions
    ! of the system_grid is made. In case of failure - an error occurs (in
    ! particular, if grid is inconsistent or not close to the system_grid).
    ! Necessity of this adjustment: if the grid is shifted by half-cell than
    ! it may happen that plus-minus one cell in any dimension is enough to
    ! cover the requested area. In this case it is mandatory to requre 
    ! exactly the same size of the matrix as that of the system_grid. 
    ! Otherwise the derived fields calculations will be destroyed
    !
    ! Method:
    ! fu_lonlat_grid_containing_area
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_grid) :: grid_new
    !
    ! Imported parameters with intent(in):
    TYPE(silja_grid), INTENT(in) :: grid
    TYPE(silam_area), INTENT(in) :: area

    ! Local declarations:
    TYPE(silja_grid), SAVE :: grid_previous, grid_new_previous
    TYPE(silam_area), SAVE :: area_previous
    TYPE(silja_logical) :: ifAdjusted

    !----------------------------------------
    !
    ! 1. Check parameters.
    !    ----------------

    IF (.NOT.defined(grid)) THEN
      CALL set_error('undefined grid given','fu_grid_containing_area')
      RETURN
    END IF

    IF (.NOT.defined(area)) THEN
      CALL set_error('undefined area given','fu_grid_containing_area')
      RETURN
    END IF

    IF ((grid == grid_previous).and. (area == area_previous)) THEN
      grid_new = grid_new_previous
      RETURN
    END IF

    ! In case of error the original grid is returned:
    grid_new = grid

    SELECT CASE (grid%gridtype)

      CASE (lonlat)

      !---------------------------------------------
      !
      ! 2. Select area from lonlat grid.
      !    ----------------------------

      grid_new = fu_lonlat_grid_containing_area(grid%lonlat, area)
      IF (error) THEN
        CALL report(grid)
        RETURN
      END IF


      case(anygrid)
      
      grid_new = fu_anygrid_containing_area(grid, area)

      IF (error) THEN
        CALL report(grid)
        RETURN
      END IF

    CASE default

      CALL set_error('sorry, working only for lonlat grids','fu_grid_containing_area')
      RETURN

    END SELECT

    grid_new_previous = grid_new
    grid_previous = grid
    area_previous = area

  END FUNCTION fu_grid_containing_area




  ! ***************************************************************


  FUNCTION fu_lonlat_grid_containing_area(gr_orig, area) result(grid_new)

    ! Description:
    ! Calculates the lonlat grid, that is the smallest possible
    ! subgrib of
    ! the given grid_orig, that contains the whole given area.
    ! The gridpoints in grid_new are all exact gridpoints of
    ! grid_orig, only area is different.
    !
    ! If the area is totally outside the given grid, than
    ! grid_missing is returned and a warning is displayed. 
    !
    ! If the area is partially covered by the given grid, but is
    ! outside in some dimension, then that border of the new grid is
    ! determined by the original grid, instead of area. So the
    ! resulting grid never contains points that aren't in the
    ! original grid.
    !
    ! Method:
    ! If same pole, simply take the corners. 
    ! If projections are different, checking only corners or even boundaries is 
    ! not enough (poles can be located inside, can cross the 180 meridian)  
    ! Have to turn the area into a small grid (10x10 cells), project the cells
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_grid) :: grid_new
    !
    ! Imported parameters with intent(in):
    TYPE(silja_lonlat_grid), INTENT(in) :: gr_orig
    TYPE(silam_area), INTENT(in) :: area
    !
    ! Local declarations:
    REAL :: area_southest_lat, area_northest_lat,&
          & area_westest_lon, area_eastest_lon
    REAL :: orig_east_border_lon, orig_north_border_lat, &
          & sw_corner_modlon_new, sw_corner_modlat_new
    real :: fx, fy
    INTEGER :: ix, iy, nx_new, ny_new, mv_x, mv_y
    TYPE(silja_grid) :: area_grid
    
    ! area and grid boundaries
    call area_southwest_corner(area, fu_pole(area), area_southest_lat, area_westest_lon)
    call area_northeast_corner(area, fu_pole(area), area_northest_lat, area_eastest_lon)
    if(area_eastest_lon <= area_westest_lon) area_eastest_lon = area_eastest_lon + 360.
    orig_north_border_lat = gr_orig%sw_corner_modlat + (gr_orig%dy_deg * (gr_orig%ny-1))
    orig_east_border_lon = gr_orig%sw_corner_modlon + (gr_orig%dx_deg * (gr_orig%nx-1))
    
    if(gr_orig%pole == fu_pole(area))then ! just corners is enough

      IF (orig_north_border_lat <= area_southest_lat .or.  gr_orig%sw_corner_modlat >= area_northest_lat .or. &
        & orig_east_border_lon <= area_westest_lon .or.  gr_orig%sw_corner_modlon >= area_eastest_lon) THEN
        CALL set_error('Area totally outside the grid', 'fu_lonlat_grid_containing_area')
        CALL report(area)
        !call report(gr_orig)
        RETURN
      END IF   
      
      mv_x = NINT((area_westest_lon - gr_orig%sw_corner_modlon)/gr_orig%dx_deg)
      sw_corner_modlon_new = max(gr_orig%sw_corner_modlon + REAL(mv_x) * gr_orig%dx_deg, gr_orig%sw_corner_modlon)
      mv_y = NINT((area_southest_lat - gr_orig%sw_corner_modlat)/gr_orig%dy_deg)
      sw_corner_modlat_new = max(gr_orig%sw_corner_modlat + REAL(mv_y) * gr_orig%dy_deg, gr_orig%sw_corner_modlat)    
      nx_new = min(NINT((area_eastest_lon - sw_corner_modlon_new) / gr_orig%dx_deg) + 1, gr_orig%nx - mv_x)
      ny_new = min(NINT((area_northest_lat - sw_corner_modlat_new) / gr_orig%dy_deg) + 1, gr_orig%ny - mv_y)
    
    else ! different projections
      ! set a temporary small grid inside the area and reproject it to the original grid
      area_grid = fu_set_lonlat_grid(fu_name(area), area_westest_lon, area_southest_lat, &
                                  & .false., 11, 11, fu_pole(area), & 
                                  & (area_eastest_lon - area_westest_lon)/10., &
                                  & (area_northest_lat - area_southest_lat)/10.)
                                  
      ! project all gridcells to original grid, record the min and max grid coordinates 
      mv_x = gr_orig%nx
      mv_y = gr_orig%ny
      nx_new = 0
      ny_new = 0
      do ix = 1, 10
        do iy = 1, 10
          call modify_lonlat(area_grid%lonlat%sw_corner_modlat + (iy-1) * area_grid%lonlat%dy_deg, &
                           & area_grid%lonlat%sw_corner_modlon + (ix-1) * area_grid%lonlat%dx_deg, &
                           & area_grid%lonlat%pole, &
                           & gr_orig%pole, &
                           & fy, &
                           & fx)
          fx = (fx - gr_orig%sw_corner_modlon) / gr_orig%dx_deg
          fy = (fy - gr_orig%sw_corner_modlat) / gr_orig%dy_deg
          mv_x = min(mv_x, nint(fx))
          mv_y = min(mv_y, nint(fy))
          nx_new = max(nx_new, nint(fx))
          ny_new = max(ny_new, nint(fy))          
        enddo
      enddo
      
      if(mv_x == gr_orig%nx .or. mv_y == gr_orig%ny .or. nx_new == 0 .or. ny_new == 0)then
        CALL set_error('Area totally outside the grid', 'fu_lonlat_grid_containing_area')
        CALL report(area)
        !call report(gr_orig)
        RETURN
      endif
   
      sw_corner_modlon_new = max(gr_orig%sw_corner_modlon + REAL(mv_x) * gr_orig%dx_deg, gr_orig%sw_corner_modlon)
      sw_corner_modlat_new = max(gr_orig%sw_corner_modlat + REAL(mv_y) * gr_orig%dy_deg, gr_orig%sw_corner_modlat)    
      nx_new = min(nx_new-mv_x+1, gr_orig%nx - mv_x)
      ny_new = min(ny_new-mv_y+1, gr_orig%ny - mv_y)
    
    endif

    ! set the new grid:
    grid_new = fu_set_lonlat_grid(fu_connect_strings('grid_of_' ,fu_name(area)),&
                                & sw_corner_modlon_new, sw_corner_modlat_new, &
                                & .false.,&
                                & nx_new, ny_new,&
                                & gr_orig%pole,&
                                & gr_orig%dx_deg,&
                                & gr_orig%dy_deg)
    IF (error) RETURN

  END FUNCTION fu_lonlat_grid_containing_area


  ! ***************************************************************


  FUNCTION fu_anygrid_containing_area(gr_orig, area) result(grid_new)

    ! Description:
    ! Calculates the lonlat grid, that is the smallest possible
    ! subgrib of
    ! the given grid_orig, that contains the whole given area.
    ! The gridpoints in grid_new are all exact gridpoints of
    ! grid_orig, only area is different.
    !
    ! If the area is totally outside the given grid, than
    ! grid_missing is returned and a warning is displayed. 
    !
    ! If the area is partially covered by the given grid, but is
    ! outside in some dimension, then that border of the new grid is
    ! determined by the original grid, instead of area. So the
    ! resulting grid never contains points that aren't in the
    ! original grid.
    !
    ! turn the area into a small grid (10x10 cells), project the cells
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_grid) :: grid_new
    !
    ! Imported parameters with intent(in):
    TYPE(silja_grid), INTENT(in) :: gr_orig
    TYPE(silam_area), INTENT(in) :: area
    !
    ! Local declarations:
    integer :: ix, iy, i, nx, ny
    real :: xmax, xmin, ymax, ymin, x, y, flat, flon
    real, dimension(:), pointer :: lons, lats, dx, dy
    
    TYPE(silja_grid) :: area_grid

    call area_southwest_corner(area, fu_pole(area), ymin, xmin)
    call area_northeast_corner(area, fu_pole(area), ymax, xmax)
    if(xmax <= xmin) xmax = xmax + 360.
    area_grid = fu_set_lonlat_grid(fu_name(area), xmin, ymin, &
                                & .false., 11, 11, fu_pole(area), & 
                                & (xmax - xmin)/10., &
                                & (ymax - ymin)/10.)
                                  
    ! project all gridcells to original grid, record the min and max grid coordinates 
    xmin = gr_orig%ag%nx
    ymin = gr_orig%ag%ny
    xmax = 0
    ymax = 0
    x = real_missing
    y = real_missing
    do ix = 1, 10
      do iy = 1, 10
        call project_point_to_grid_xy(area_grid, real(ix), real(iy), gr_orig, x, y)
        xmin = min(xmin, x)
        ymin = min(ymin, y)
        xmax = max(xmax, x)
        ymax = max(ymax, y)          
      enddo
    enddo
      
    if(nint(xmin) == gr_orig%ag%nx .or. nint(ymin) == gr_orig%ag%ny .or. &
     & nint(xmax) == 0 .or. nint(ymax) == 0)then
      CALL set_error('Area totally outside the grid', 'fu_anygrid_containing_area')
      CALL report(area)
      call report(gr_orig)
      RETURN
    endif
   
    nx = nint(xmax)-nint(xmin) +1
    ny = nint(ymax)-nint(ymin) +1
  
    lons => fu_work_array()
    lats => fu_work_array()
    dx => fu_work_array()
    dy => fu_work_array()

    i = 1
    do iy = nint(ymin), nint(ymax)
      do ix = nint(xmin), nint(xmax)
      
        lons(i) = pAnyGrdParam(gr_orig%ag%indParam)%xc((iy-1)*gr_orig%ag%nx+ix)
        lats(i) = pAnyGrdParam(gr_orig%ag%indParam)%xc((iy-1)*gr_orig%ag%nx+ix)
        dx(i) = pAnyGrdParam(gr_orig%ag%indParam)%xc((iy-1)*gr_orig%ag%nx+ix)
        dy(i) = pAnyGrdParam(gr_orig%ag%indParam)%xc((iy-1)*gr_orig%ag%nx+ix)
      
        i = i+1
       
      enddo
    enddo

    grid_new = fu_set_any_grid(gr_orig%name, nx, ny, lons, lats)
    call setAnygridParam(grid_new, 'dx', dx)
    call setAnygridParam(grid_new, 'dy', dy)


   call free_work_array(lons)   
   call free_work_array(lats)   
   call free_work_array(dx)   
   call free_work_array(dy)

  end FUNCTION fu_anygrid_containing_area

  ! ***************************************************************


  subroutine cut_grid_size_grid(grid_to_cut, grid_area, iTypeOfCut, grid_to_preserve)
    !
    ! Depending on iTypeOfCut, cuts the grid_to_cut down to just-covering 
    ! or down to be fully inside the grid_area.
    ! The third grid, if given, is the one that must stay covered by the cut grid.
    !
    ! Algorithm contains several steps:
    ! 1. All points of the grid_to_cut are projected in grid_area 
    ! 2. Then the largest possible (in fact, just large enough) sub-grid
    !    is selected by iterational cutting off the non-covered lines and/or cells
    ! 3. During iterations, the grid_to_preserve is kept, essentially claiming very hgih
    !    costs for cutting the line that touches it
    !
    implicit none
    
    ! Imported variables
    type(silja_grid), intent(inout) :: grid_to_cut
    type(silja_grid), intent(in), target :: grid_area
    type(silja_grid), intent(in), target, optional :: grid_to_preserve
    integer, intent(in) :: iTypeOfCut

    ! Local variables
    type(silja_grid), target :: grdTmp
!    type(silja_lonlat_grid), pointer :: gtc,ga ! Just to shorten the notations
    integer, dimension(:,:), pointer :: arFlag
    integer, dimension(:), pointer :: iWork
    integer :: ix_start_count, ix_end_count, iy_start_count, iy_end_count, iCut
    integer :: i,j, ix_start, ix_end, iy_start, iy_end  !, iCnt !, iUnit
    real :: corner_x, corner_y, pole_x, pole_y, dx_deg, dy_deg
    real :: fXStart, fYStart, fXEnd, fYEnd 
    logical :: if_south_pole, if_corner_in_geo_coord, ifnew, ifAllocated
    integer :: nx, ny, grid_type
    type(silam_any_grid_param), pointer :: ptr_grid_pars
    real, dimension(:), pointer :: lon_tmp, lat_tmp, dx_tmp, dy_tmp


    if(fu_grids_arakawa_correspond(grid_to_cut, grid_area))return

    grdTmp = grid_to_cut ! If fail - the grid_to_cut will be returned without damage
    !
    ! No cut if covered area is global, whatever the request is
    !
!    if(fu_ifLonGlobal(grid_area)) then 
!            call msg("Only lonGlobal is checked by cut_grid_size")
!            call msg("Not cutting anything")
!            return
!    endif
    if(fu_ifLonGlobal(grid_area) .and. fu_ifLatGlobal(grid_area)) then
         call msg("Global grid, Not cutting anything")
         return
    endif

    SELECT CASE(fu_gridtype(grid_to_cut))
      CASE(lonlat)
        call lonlat_grid_parameters(grid_to_cut, &
                                  & corner_x, corner_y, if_corner_in_geo_coord, &
                                  & nx, ny, &
                                  & pole_x, pole_y, & 
                                  & dx_deg, dy_deg)
      case(anygrid)

        call grid_dimensions(grid_to_cut, nx, ny)  
        ptr_grid_pars => pAnyGrdParam(grid_to_cut%ag%indParam)
 
      case default
        call set_error('strange grid type','cut_grid_size')
    end select
    if(error)return

    iWork => fu_work_int_array(nx*ny)
    arFlag(1:nx,1:ny) => iWork(1:nx*ny)
    arFlag(1:nx, 1:ny) = -1

    !
    ! First, fill-in with 1 the cells of gtc, which belong to ga
    !
    call fill_overlap_flag_array(grdTmp, grid_area, arFlag, iTypeOfCut, 1)
    if(error)return

    if(present(grid_to_preserve))then
      !
      ! Stupidity check: grid to preserve must be preservable
      !
      if(.not. fu_if_grid_covered(grid_to_preserve, grid_area))then ! sml_grd, lrg_grd
        call set_error('Grid to preserve is not covered with grid_area','cut_grid_size')
        return
      endif
      if(.not. fu_if_grid_covered(grid_to_preserve, grid_to_cut))then ! sml_grd, lrg_grd
        call set_error('Grid to preserve is not covered with grid_to_cut','cut_grid_size')
        return
      endif
      !
      ! Now, fill-in with 2 those cells of gtc, which belong to grid_to_preserve
      !
      call fill_overlap_flag_array(grdTmp, grid_to_preserve, arFlag, cover_the_grid_area, 2)
      if(error)return
    endif

    !
    ! Remove empty lines of gtc - which entirely out of ga
    !
    ix_start = 1
    iy_start = 1
    ix_end = nx
    iy_end = ny

    call cut_empty_lines(arFlag, ix_start, iy_start, ix_end, iy_end, &
                       & -1, 1, 2) !value_to_cut, value_to_stay, value_to_preserve)
    !
    ! We cut out all the lines, which have no connection at all with the grid_area. 
    ! So, we got the grid just-covering the grid_area. Return if this was the task
    !
    select case(iTypeOfCut)
      
      case(cover_the_grid_area)
        !
        ! Only fully-emtpy line to be removed.
        ! Do nothing - all is above but below we should set the grid and free memory

      case(inside_the_grid_area)
        !
        ! Let's remove one-by-one verticals and horizontals not completely covered, each 
        ! time selecting one for cutting according to one of the following criteria:
        ! 
        ! Possible options are:
        ! - select that one with smallest number of good cells belonging to grid_area (BAD-tried)
        ! - select that one with largest number of cells not belonging to grid_area
        ! - select that one where the difference between the numbers of cells 
        !   not belonging and belonging is maximum
        ! - select the line with a max fraction of "bad" cells
        !
        do while(any(arFlag(ix_start:ix_end,iy_start:iy_end) == -1))
          !
          ! Count number of covered elements for each of boundary lines and
          ! add to this target function an extra penalty for non-square grid.
          ! However, do not be too strong - the user may want to have non-square
          ! grid
          !
          iCut = 0
          SELECT CASE(fu_gridtype(grid_to_cut)) ! Define the grid for ctl file
            CASE(lonlat)
               i = (sign(1,int(nx*dx_deg - ny*dy_deg))) ! penalty for non-aquare
            case(anygrid)
               ! don't know what we're doing here anyway
               i = 0
            case default
               call set_error('strange grid type','cut_grid_size')
          end select
          ix_start_count = count(arFlag(ix_start,iy_start:iy_end) == -1) + i
          ix_end_count = count(arFlag(ix_end,iy_start:iy_end) == -1) + i
          iy_start_count = count(arFlag(ix_start:ix_end,iy_start) == -1) - i
          iy_end_count = count(arFlag(ix_start:ix_end,iy_end) == -1) - i

!          write(iUnit, *)'ix_start_count', ix_start_count
!          write(iUnit, *)'ix_end_count', ix_end_count
!          write(iUnit, *)'iy_start_count', iy_start_count
!          write(iUnit, *)'iy_end_count', iy_end_count
!          write(iUnit, *)

          fXStart = real(ix_start_count) / real(nx)  ! a bad fraction of the line
          fYStart = real(iy_start_count) / real(ny)  ! a bad fraction of the line
          fXEnd = real(ix_end_count) / real(nx)  ! a bad fraction of the line
          fYEnd = real(iy_end_count) / real(ny)  ! a bad fraction of the line

          !
          ! Forbid cutting of the grid_to_preserve area
          !
          if(any(arFlag(ix_start,iy_start:iy_end) == 2)) fXStart = -1.0 ! negative!
          if(any(arFlag(ix_end,iy_start:iy_end) == 2))   fXEnd = -1.0
          if(any(arFlag(ix_start:ix_end,iy_start) == 2)) fYStart = -1.0
          if(any(arFlag(ix_start:ix_end,iy_end) == 2))   fYEnd = -1.0

          !
          ! Find and cut the worst line for each bad corner point that does not 
          ! For each point we just select the right side to cut
          !
          if(arFlag(ix_start,iy_start) == -1)then  ! fXStart vs fYStart
            if(fXStart < 0 .and. fYStart < 0)then
              call set_error('Failed to cut fXStart - fYStart corner','cut_grid_size')
              return
            endif
            if(fXStart > fYStart)then
              ix_start = ix_start + 1
            else
              iy_start = iy_start + 1
            endif
            iCut = iCut + 1
          end if

          if(arFlag(ix_start,iy_end) == -1)then
            if(fXStart < 0 .and. fYEnd < 0)then
              call set_error('Failed to cut fXStart - fYEnd corner','cut_grid_size')
              return
            endif
            if(fXStart > fYEnd)then
              ix_start = ix_start + 1
            else
              iy_end = iy_end - 1
            endif
            iCut = iCut + 1
          end if

          if(arFlag(ix_end,iy_start) == -1)then
            if(fXEnd < 0 .and. fYStart < 0)then
              call set_error('Failed to cut fXEnd - fYStart corner','cut_grid_size')
              return
            endif
            if(fXEnd > fYStart)then
              ix_end = ix_end - 1
            else
              iy_start = iy_start + 1
            endif
            iCut = iCut + 1
          end if

          if(arFlag(ix_end,iy_end) == -1)then
            if(fXEnd < 0 .and. fYEnd < 0)then
              call set_error('Failed to cut fXEnd - fYEnd corner','cut_grid_size')
              return
            endif
            if(fXEnd > fYEnd)then
              ix_end = ix_end - 1
            else
              iy_end = iy_end - 1
            endif
            iCut = iCut + 1
          end if

          call cut_empty_lines(arFlag, ix_start, iy_start, ix_end, iy_end, &
                             & -1, 1, 2)     ! value_to_cut, value_to_stay, value_to_preserve)
          !
          ! If we cannot cut anything but -1 is still with us - may be, some lines are 
          ! jumping out in the middle?
          !
          if(iCut < 1)then
            call set_error('Failed grid cutting for unknown reason','cut_grid_size')
            return
          endif  ! iCut < 1

          if(ix_end - ix_start < 2)then
            call set_error('X-axis: grids do not overlap','cut_grid_size')
            return
          endif
          if(iy_end - iy_start < 2)then
            call set_error('Y-axis: grids do not overlap','cut_grid_size')
            return
          endif

        end do ! cutting incomplete verticals and horizontals

      case default
        call msg('Strange type of the grid cut:',iTypeOfCut)
        call set_error('Strange type of the grid cut','cut_grid_size')
        return
    end select  ! type of cut

    !
    ! Cutting done. Record the results and set the output grid
    !
    call cut_grid_size_params(grid_to_cut, ix_start, iy_start, ix_end, iy_end)

    call free_work_array(iWork)

    call msg('cut grid')
    call report(grid_to_cut)

  end subroutine cut_grid_size_grid


  !*******************************************************************************

  subroutine fill_overlap_flag_array(grid_input, grid_template, arFlag, iTypeOfCut, val)
      !
      ! Fills-in the arFlag with the value reprojecting the input grid onto the grid_template
      !
      implicit none

      ! Imported parameters
      type(silja_grid), target, intent(in) :: grid_input, grid_template
      integer, dimension(:,:), pointer :: arFlag
      integer, intent(in) :: iTypeOfCut, val

      ! Local variables
      type(silja_lonlat_grid), pointer :: gtc,ga ! Just to shorten the notations
      real :: lon_ga, lat_ga  ! Co-ordinates in grid_area 
      integer :: i,j, iMax, jMax, ix_start, ix_end, iy_start, iy_end, iCnt !, iUnit
      integer :: nx_gtc, ny_gtc, nx_ga, ny_ga, ix, iy
      real :: xNew, yNew

        !
        ! Fill-in the array telling whether particular cell of grid_to_cut
        ! belongs to grid_area, or not
        !
        if(grid_input%gridtype == lonlat .and. grid_template%gridtype == lonlat)then

          gtc => grid_input%lonlat
          ga => grid_template%lonlat
          do j=gtc%ny,1,-1
            do i=1,gtc%nx
!              arFlag(i,j) = -1
              call modify_lonlat(gtc%sw_corner_modlat + (j-1.5)*gtc%dy_deg, &
                               & gtc%sw_corner_modlon + (i-1.5)*gtc%dx_deg, gtc%pole, &
                               & ga%pole, lat_ga, lon_ga)
              if(lon_ga >= ga%sw_corner_modlon - 0.5*ga%dx_deg .and. &
               & lon_ga <= ga%sw_corner_modlon + (ga%nx-0.5)*ga%dx_deg .and. &
               & lat_ga >= ga%sw_corner_modlat  - 0.5*ga%dy_deg .and. &
               & lat_ga <= ga%sw_corner_modlat + (ga%ny-0.5)*ga%dy_deg)then
                 if(iTypeOfCut == cover_the_grid_area)then
                   arFlag(i,j) = val
                   cycle
                 endif
              else
                 if(iTypeOfCut == inside_the_grid_area)then
                   cycle
                 endif
              endif
              call modify_lonlat(gtc%sw_corner_modlat + (j-0.5)*gtc%dy_deg, &
                               & gtc%sw_corner_modlon + (i-1.5)*gtc%dx_deg, gtc%pole, &
                               & ga%pole, lat_ga, lon_ga)
              if(lon_ga >= ga%sw_corner_modlon - 0.5*ga%dx_deg .and. &
               & lon_ga <= ga%sw_corner_modlon + (ga%nx-0.5)*ga%dx_deg .and. &
               & lat_ga >= ga%sw_corner_modlat - 0.5*ga%dy_deg .and. &
               & lat_ga <= ga%sw_corner_modlat + (ga%ny-0.5)*ga%dy_deg)then
                 if(iTypeOfCut == cover_the_grid_area)then
                   arFlag(i,j) = val
                   cycle
                 endif
              else
                if(iTypeOfCut == inside_the_grid_area)then
                  cycle
                endif
              endif
              call modify_lonlat(gtc%sw_corner_modlat + (j-1.5)*gtc%dy_deg, &
                               & gtc%sw_corner_modlon + (i-0.5)*gtc%dx_deg, gtc%pole, &
                               & ga%pole, lat_ga, lon_ga)
              if(lon_ga >= ga%sw_corner_modlon - 0.5*ga%dx_deg .and. &
               & lon_ga <= ga%sw_corner_modlon + (ga%nx-0.5)*ga%dx_deg .and. &
               & lat_ga >= ga%sw_corner_modlat - 0.5*ga%dy_deg .and. &
               & lat_ga <= ga%sw_corner_modlat + (ga%ny-0.5)*ga%dy_deg)then
                 if(iTypeOfCut == cover_the_grid_area)then
                   arFlag(i,j) = val
                   cycle
                 endif
              else
                 if(iTypeOfCut == inside_the_grid_area)then
                   cycle
                 endif
              endif
              call modify_lonlat(gtc%sw_corner_modlat + (j-0.5)*gtc%dy_deg, &
                               & gtc%sw_corner_modlon + (i-0.5)*gtc%dx_deg, gtc%pole, &
                               & ga%pole, lat_ga, lon_ga)
              if(lon_ga >= ga%sw_corner_modlon - 0.5*ga%dx_deg .and. &
               & lon_ga <= ga%sw_corner_modlon + (ga%nx-0.5)*ga%dx_deg .and. &
               & lat_ga >= ga%sw_corner_modlat - 0.5*ga%dy_deg .and. &
               & lat_ga <= ga%sw_corner_modlat + (ga%ny-0.5)*ga%dy_deg)then
                 arFlag(i,j) = val
              endif
            end do
          end do
          !
          ! If we want to cover the grid_area, we have to check the corners of that grid - 
          ! they do not necessarily touch the corners of the gtc.
          !
          if(iTypeOfCut == cover_the_grid_area)then
            call modify_lonlat(ga%sw_corner_modlat -1.5*ga%dy_deg, &
                             & ga%sw_corner_modlon -1.5*ga%dx_deg, ga%pole, &
                               & gtc%pole, lat_ga, lon_ga)
            i = nint((lon_ga - gtc%sw_corner_modlon) / gtc%dx_deg) + 1

            j = nint((lat_ga - gtc%sw_corner_modlat) / gtc%dy_deg) + 1

           if(i > 0 .and. i <= gtc%nx .and. j > 0 .and. j <= gtc%ny) arFlag(i,j)=val

            call modify_lonlat(ga%sw_corner_modlat -1.5*ga%dy_deg, &
                            & ga%sw_corner_modlon +(ga%nx+0.5)*ga%dx_deg, ga%pole, &
                            & gtc%pole, lat_ga, lon_ga)
            i = nint((lon_ga - gtc%sw_corner_modlon) / gtc%dx_deg) + 1

            j = nint((lat_ga - gtc%sw_corner_modlat) / gtc%dy_deg) + 1

           if(i > 0 .and. i <= gtc%nx .and. j > 0 .and. j <= gtc%ny) arFlag(i,j)=val

            call modify_lonlat(ga%sw_corner_modlat +(ga%ny+0.5)*ga%dy_deg, &
                            & ga%sw_corner_modlon -1.5*ga%dx_deg, ga%pole, &
                            & gtc%pole, lat_ga, lon_ga)
            i = nint((lon_ga - gtc%sw_corner_modlon) / gtc%dx_deg) + 1

            j = nint((lat_ga - gtc%sw_corner_modlat) / gtc%dy_deg) + 1

           if(i > 0 .and. i <= gtc%nx .and. j > 0 .and. j <= gtc%ny) arFlag(i,j)=val

           call modify_lonlat(ga%sw_corner_modlat +(ga%ny+0.5)*ga%dy_deg, &
                            & ga%sw_corner_modlon +(ga%nx+0.5)*ga%dx_deg, ga%pole, &
                            & gtc%pole, lat_ga, lon_ga)
            i = nint((lon_ga - gtc%sw_corner_modlon) / gtc%dx_deg) + 1

            j = nint((lat_ga - gtc%sw_corner_modlat) / gtc%dy_deg) + 1

           if(i > 0 .and. i <= gtc%nx .and. j > 0 .and. j <= gtc%ny) arFlag(i,j)=val
          endif


!        case(anygrid)
        else

          call grid_dimensions(grid_input, nx_gtc, ny_gtc)
          call grid_dimensions(grid_template, nx_ga, ny_ga)
          !$OMP PARALLEL DEFAULT(NONE) SHARED(ny_gtc, nx_gtc,grid_template, grid_input,&
          !$OMP & nx_ga, ny_ga, iTypeOfCut,val, arFlag) &
          !$OMP & PRIVATE(xNew, yNew, ix, iy)

          yNew  = real_missing
          !$OMP DO
          do iy = 1, ny_gtc
            if (modulo(iy, 100) == 0) call msg("fill_overlap_flag_array", iy, ny_gtc)
            !reset on discontinuity
            xNew = real_missing
            do ix = 1, nx_gtc
              call project_point_to_grid(grid_input, ix-0.5, iy-0.5, grid_template, xNew, yNew)
              if(xNew >= 0.5 .and. yNew >= 0.5 .and. xNew <= nx_ga+0.5 .and. yNew <= ny_ga+0.5)then
                 if(iTypeOfCut == cover_the_grid_area)then
                   arFlag(ix,iy) = val
                   cycle
                 endif
              else
                if(iTypeOfCut == inside_the_grid_area)then
                  cycle
                endif
              endif  
              call project_point_to_grid(grid_input, ix+0.5, iy+0.5, grid_template, xNew, yNew)
              if(xNew >= 0.5 .and. yNew >= 0.5 .and. xNew <= nx_ga+0.5 .and. yNew <= ny_ga+0.5)then
                 if(iTypeOfCut == cover_the_grid_area)then
                   arFlag(ix,iy) = val
                   cycle
                 endif
              else
                if(iTypeOfCut == inside_the_grid_area)then
                  cycle
                endif
              endif
              call project_point_to_grid(grid_input, ix-0.5, iy+0.5, grid_template, xNew, yNew)
              if(xNew >= 0.5 .and. yNew >= 0.5 .and. xNew <= nx_ga+0.5 .and. yNew <= ny_ga+0.5)then
                 if(iTypeOfCut == cover_the_grid_area)then
                   arFlag(ix,iy) = val
                   cycle
                 endif
              else
                if(iTypeOfCut == inside_the_grid_area)then
                  cycle
                endif
              endif
              call project_point_to_grid(grid_input, ix+0.5, iy-0.5, grid_template, xNew, yNew)
              if(xNew >= 0.5 .and. yNew >= 0.5 .and. xNew <= nx_ga+0.5 .and. yNew <= ny_ga+0.5)then
                 arFlag(ix,iy) = val
              endif

            end do
          end do
          !$OMP END DO
          !$OMP BARRIER


          if(iTypeOfCut == cover_the_grid_area)then

            ! probably should be grid borders here, not just corners ..

            !$OMP DO
            do iy = 1, ny_ga+1
               call project_point_to_grid(grid_template, 0.5, iy-0.5, grid_input, xNew, yNew)
               if(xNew > 0.5 .and. xNew <= nx_gtc .and. yNew > 0.5 .and. yNew <= ny_gtc) &
                                                                   & arFlag(nint(xNew),nint(yNew))=val
               call project_point_to_grid(grid_template, nx_ga+0.5, iy-0.5, grid_input, xNew, yNew)
               if(xNew > 0.5 .and. xNew <= nx_gtc .and. yNew > 0.5 .and. yNew <= ny_gtc) &
                                                                   & arFlag(nint(xNew),nint(yNew))=val
            enddo
            !$OMP END DO
            !$OMP DO
            do ix = 1, nx_ga+1
               call project_point_to_grid(grid_template, ix-0.5, 0.5, grid_input, xNew, yNew)
               if(xNew > 0.5 .and. xNew <= nx_gtc .and. yNew > 0.5 .and. yNew <= ny_gtc) &
                                                                   & arFlag(nint(xNew),nint(yNew))=val
               call project_point_to_grid(grid_template, ix-0.5, ny_ga+0.5, grid_input, xNew, yNew)
               if(xNew > 0.5 .and. xNew <= nx_gtc .and. yNew > 0.5 .and. yNew <= ny_gtc) &
                                                                   & arFlag(nint(xNew),nint(yNew))=val
            enddo
            !$OMP END DO
          endif
          !$OMP END PARALLEL 
        endif

  end subroutine fill_overlap_flag_array


  !*****************************************************************************

  subroutine cut_empty_lines(arFlag, ix_start, iy_start, ix_end, iy_end, &
                           & value_to_cut, value_to_stay, value_to_preserve)
      !
      ! Remove all totally empty verticals and horizontals
      ! and the lines with missing middle points 
      !
      implicit none
      
      ! Imported parameters
      integer, dimension(:,:), pointer :: arFlag
      integer, intent(inout) :: ix_start, iy_start, ix_end, iy_end
      integer, intent(in) :: value_to_cut, value_to_stay, value_to_preserve
      
      ! Local variables
      logical :: ifCut

      do while(all(arFlag(ix_start,iy_start:iy_end) == value_to_cut)) ! cut from x start
        ix_start = ix_start+1
        if(ix_end - ix_start < 1)then
          call set_error('X-axis: grids do not overlap','cut_empty_lines')
          return
        endif
      end do
      do while(all(arFlag(ix_end,iy_start:iy_end) == value_to_cut)) ! cut from x end
        ix_end = ix_end-1
        if(ix_end - ix_start < 1)then
          call set_error('X-axis: grids do not overlap','cut_empty_lines')
          return
        endif
      end do
      do while(all(arFlag(ix_start:ix_end,iy_start) == value_to_cut)) ! cut from y start
        iy_start = iy_start+1
        if(iy_end - iy_start < 1)then
          call set_error('Y-axis: grids do not overlap','cut_empty_lines')
          return
        endif
      end do
      do while(all(arFlag(ix_start:ix_end,iy_end) == value_to_cut)) ! cut from iy_end
        iy_end = iy_end-1
        if(iy_end - iy_start < 1)then
          call set_error('Y-axis: grids do not overlap','cut_empty_lines')
          return
        endif
      end do
      !
      ! Check and cut the lines with corners belonging to the grid and middle-points
      ! placed outside. These lines are self-responsible, so no tricks can be made for them
      !
      ifCut = .true.
      do while(ifCut)
        ifCut = .false.
        !
        ! Check the x-start line
        !
        if(arFlag(ix_start,iy_start) == value_to_stay .and. &
         & arFlag(ix_start,iy_end) == value_to_stay .and. &
         & any(arFlag(ix_start,iy_start:iy_end) == value_to_cut))then
          if(any(arFlag(ix_start,iy_start:iy_end) == value_to_preserve))then
            call set_error('Incompatible grids ix_start','cut_empty_lines')
            return
          endif
          ifCut = .true.
          ix_start = ix_start + 1
        endif
        if(arFlag(ix_start,iy_start) == value_to_stay .and. &
         & arFlag(ix_end,iy_start) == value_to_stay .and. &
         & any(arFlag(ix_start:ix_end,iy_start) == value_to_cut))then
          if(any(arFlag(ix_start:ix_end,iy_start) == value_to_preserve))then
            call set_error('Incompatible grids iy_start','cut_empty_lines')
            return
          endif
          ifCut = .true.
          iy_start = iy_start + 1
        endif
        if(arFlag(ix_end,iy_start) == value_to_stay .and. &
         & arFlag(ix_end,iy_end) == value_to_stay .and. &
         & any(arFlag(ix_end,iy_start:iy_end) == value_to_cut))then
          if(any(arFlag(ix_end,iy_start:iy_end) == value_to_preserve))then
            call set_error('Incompatible grids ix_end','cut_empty_lines')
            return
          endif
          ifCut = .true.
          ix_end = ix_end - 1
        endif
        if(arFlag(ix_start,iy_end) == value_to_stay .and. &
         & arFlag(ix_end,iy_end) == value_to_stay .and. &
         & any(arFlag(ix_start:ix_end,iy_end) == value_to_cut))then
          if(any(arFlag(ix_start:ix_end,iy_end) == value_to_preserve))then
            call set_error('Incompatible grids iy_end','cut_empty_lines')
            return
          endif
          ifCut = .true.
          iy_end = iy_end - 1
        endif
      end do
    end subroutine cut_empty_lines


  !*****************************************************************

  subroutine cut_grid_size_params(grid_to_cut, ix_start, iy_start, ix_end, iy_end)
    !
    ! Reduces the grid size down to the given start-end indices
    !
    implicit none
    
    ! Imported parameters
    type(silja_grid), intent(inout) :: grid_to_cut
    integer, intent(in) :: ix_start, iy_start, ix_end, iy_end

    ! Local variables
    integer :: iCut, i,j, nx, ny, nxnew, nynew
    type(silja_grid) :: grid_new
    type(silam_any_grid_param), pointer :: ptr_grid_pars, ptr_grid_pars_new
    real, dimension(:), pointer :: lon_tmp, lat_tmp
    real, dimension(:,:), pointer :: ptr2dSrc, ptr2dDst
    logical :: ifNew
    character(len = *), parameter :: sub_name = 'cut_grid_size_params'

    nxnew = ix_end - ix_start + 1
    nynew = iy_end - iy_start + 1

    select case (grid_to_cut%gridtype)
    case (lonlat)
      grid_to_cut%name = trim(grid_to_cut%name)//"_cut"
      grid_to_cut%lonlat%sw_corner_modlon = grid_to_cut%lonlat%sw_corner_modlon + &
                                                    & (ix_start-1) * grid_to_cut%lonlat%dx_deg
      grid_to_cut%lonlat%sw_corner_modlat = grid_to_cut%lonlat%sw_corner_modlat + &
                                                    & (iy_start-1) * grid_to_cut%lonlat%dy_deg
      grid_to_cut%lonlat%nx = nxnew
      grid_to_cut%lonlat%ny = nynew
 
    case(anygrid)
      !
      ! Have to change fields in grid params, including their sizes   
      !
      call grid_dimensions(grid_to_cut, nx, ny)

      ptr_grid_pars => pAnyGrdParam(grid_to_cut%ag%indParam)


      lon_tmp => fu_work_array(nxnew*nynew)
      lat_tmp => fu_work_array(nxnew*nynew)

      ! Lons
      ptr2dDst(1:nxnew,1:nynew) => lon_tmp(1:nxnew*nynew)
      ptr2dSrc(1:nx,1:ny) => ptr_grid_pars%xC(1:nx*ny)
      ptr2dDst(1:nxnew,1:nynew) = ptr2dSrc(ix_start:ix_end, iy_start:iy_end)

      ! Lats
      ptr2dDst(1:nxnew,1:nynew) => lat_tmp(1:nxnew*nynew)
      ptr2dSrc(1:nx,1:ny) => ptr_grid_pars%yC(1:nx*ny)
      ptr2dDst(1:nxnew,1:nynew) = ptr2dSrc(ix_start:ix_end, iy_start:iy_end)

      grid_new = fu_set_any_grid(trim(grid_to_cut%name)//"_cut", nxnew, nynew, lon_tmp, lat_tmp)
      call free_work_array(lon_tmp)
      call free_work_array(lat_tmp)

      ptr_grid_pars_new => pAnyGrdParam(grid_new%ag%indParam)

      if ( .not. allocated(ptr_grid_pars_new%dx)) then ! Cut grid is new

        ! Have to fill stuff for new grid
         allocate(ptr_grid_pars_new%dx(nxnew*nynew), ptr_grid_pars_new%dy(nxnew*nynew), &
           & ptr_grid_pars_new%x3dC(nxnew*nynew), ptr_grid_pars_new%y3dC(nxnew*nynew), &
           & ptr_grid_pars_new%z3dC(nxnew*nynew))

         !dx
         ptr2dSrc(1:nx,1:ny) => ptr_grid_pars%dx(1:nx*ny)
         ptr2dDst(1:nxnew,1:nynew) => ptr_grid_pars_new%dx(1:nxnew*nynew)
         ptr2dDst(1:nxnew,1:nynew) = ptr2dSrc(ix_start:ix_end, iy_start:iy_end)


         !dy
         ptr2dSrc(1:nx,1:ny) => ptr_grid_pars%dy(1:nx*ny)
         ptr2dDst(1:nxnew,1:nynew) => ptr_grid_pars_new%dy(1:nxnew*nynew)
         ptr2dDst(1:nxnew,1:nynew) = ptr2dSrc(ix_start:ix_end, iy_start:iy_end)


         !x3dC
         ptr2dSrc(1:nx,1:ny) => ptr_grid_pars%x3dC(1:nx*ny)
         ptr2dDst(1:nxnew,1:nynew) => ptr_grid_pars_new%x3dC(1:nxnew*nynew)
         ptr2dDst(1:nxnew,1:nynew) = ptr2dSrc(ix_start:ix_end, iy_start:iy_end)
         !y3dC
         ptr2dSrc(1:nx,1:ny) => ptr_grid_pars%y3dC(1:nx*ny)
         ptr2dDst(1:nxnew,1:nynew) => ptr_grid_pars_new%y3dC(1:nxnew*nynew)
         ptr2dDst(1:nxnew,1:nynew) = ptr2dSrc(ix_start:ix_end, iy_start:iy_end)
         !z3dC
         ptr2dSrc(1:nx,1:ny) => ptr_grid_pars%z3dC(1:nx*ny)
         ptr2dDst(1:nxnew,1:nynew) => ptr_grid_pars_new%z3dC(1:nxnew*nynew)
         ptr2dDst(1:nxnew,1:nynew) = ptr2dSrc(ix_start:ix_end, iy_start:iy_end)

         if(allocated(ptr_grid_pars%sin_map_rot))then
           allocate(ptr_grid_pars_new%sin_map_rot(nxnew*nynew), &
             & ptr_grid_pars_new%cos_map_rot(nxnew*nynew))
           !sin_map_rot
           ptr2dSrc(1:nx,1:ny) => ptr_grid_pars%sin_map_rot(1:nx*ny)
           ptr2dDst(1:nxnew,1:nynew) => ptr_grid_pars_new%sin_map_rot(1:nxnew*nynew)
           ptr2dDst(1:nxnew,1:nynew) = ptr2dSrc(ix_start:ix_end, iy_start:iy_end)
           !cos_map_rot
           ptr2dSrc(1:nx,1:ny) => ptr_grid_pars%cos_map_rot(1:nx*ny)
           ptr2dDst(1:nxnew,1:nynew) => ptr_grid_pars_new%cos_map_rot(1:nxnew*nynew)
           ptr2dDst(1:nxnew,1:nynew) = ptr2dSrc(ix_start:ix_end, iy_start:iy_end)
         endif
       endif

#ifdef DEBUG
       if (any(abs(ptr_grid_pars_new%x3dC(:)*ptr_grid_pars_new%x3dC(:) + &
                & ptr_grid_pars_new%y3dC(:)*ptr_grid_pars_new%y3dC(:) + &
                & ptr_grid_pars_new%z3dC(:)*ptr_grid_pars_new%z3dC(:) - 1.)>1e-6)) then

          call msg("ptr_grid_pars_new%xC(1:10)  ", ptr_grid_pars_new%xC(1:10)) 
          call msg("ptr_grid_pars_new%yC(1:10)  ", ptr_grid_pars_new%yC(1:10)) 
          call msg("ptr_grid_pars_new%x3dC(1:10)", ptr_grid_pars_new%x3dC(1:10)) 
          call msg("ptr_grid_pars_new%y3dC(1:10)", ptr_grid_pars_new%y3dC(1:10)) 
          call msg("ptr_grid_pars_new%z3dC(1:10)", ptr_grid_pars_new%z3dC(1:10)) 
          call msg("Should be ones:", ptr_grid_pars_new%x3dC(1:10)*ptr_grid_pars_new%x3dC(1:10) + &
                & ptr_grid_pars_new%y3dC(1:10)*ptr_grid_pars_new%y3dC(1:10) + &
                & ptr_grid_pars_new%z3dC(1:10)*ptr_grid_pars_new%z3dC(1:10))
          call set_error("Wrong 3D metrics of anygrid", sub_name)
       endif
#endif
      grid_to_cut = grid_new

       

    case default
      call set_error('Cannot update unknown grid type', sub_name)
      return
    end select

  end subroutine cut_grid_size_params


  ! ***************************************************************


  SUBROUTINE coriolis_parameters(grid, coriolis)
  
    ! Description:
    ! Calculates coriolis parameter for grid's every gridpoint.
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
    TYPE(silja_grid), INTENT(in) :: grid
    
    ! Imported parameters with intent(in):
    REAL, DIMENSION(:), INTENT(out) :: coriolis

    ! Local declarations:

    REAL :: geolat
    INTEGER :: ix, iy, nx, ny

    IF (SIZE(coriolis) < fu_number_of_gridpoints(grid)) THEN
      CALL set_error('result vector too small',&
          & 'coriolis_parameters')
      RETURN
    END IF
    
    coriolis(:) = real_missing
    call grid_dimensions(grid, nx, ny)

    do ix = 1, nx
      do iy = 1, ny
        geolat = fu_lat_geographical_from_grid(real(ix), real(iy), grid)
        if(error)return
        coriolis(nx*(iy-1)+ix) = 2.*earth_omega*SIN(geolat*degrees_to_radians)
      enddo
    enddo
      
  END SUBROUTINE coriolis_parameters


  !*****************************************************************

  subroutine project_point_to_grid_xy(gridIni, x, y, gridOut, xNew, yNew)
    !
    ! The routine projects a point given in RELATIVE co-ordinates x,y in the
    ! gridIni to the gridOut. Results are also relative co-ordinates ixNew, iyNew.
    ! Since it may be used very many times, NO CHECKING is included.
    !
    implicit none

    ! Imported parameters with intent IN
    type(silja_grid), intent(in) :: gridIni, gridOut
    real, intent(in) :: x,y

    ! Imported parameters with intent OUT
    real, intent (inout) :: xNew, yNew !! Needed for anygrid initiall guess


    ! Local variables
    integer :: ixC, iyC, ix, iy, ixD, iyD
    real, dimension(-1:1,-1:1) ::  distAr
    real :: fLon, fLat


    integer :: method
    method = linear


 !   call report(gridOut)

    ! Shortcut lonlat -> lonlat
    if (gridIni%gridtype == lonlat .and. gridOut%gridtype == lonlat) then

       call modify_lonlat(gridIni%lonlat%sw_corner_modlat + (y-1) * gridIni%lonlat%dy_deg, &
                        & gridIni%lonlat%sw_corner_modlon + (x-1) * gridIni%lonlat%dx_deg, &
                        & gridIni%lonlat%pole, &
                        & gridOut%lonlat%pole, &
                        & yNew, &
                        & xNew)
       xNew = 1+ (xNew - gridOut%lonlat%sw_corner_modlon) / gridOut%lonlat%dx_deg
       yNew = 1+ (yNew - gridOut%lonlat%sw_corner_modlat) / gridOut%lonlat%dy_deg

       ! If longitude out of grid, try rotate the earth once
       if(xNew < 0.5)then
         if(xNew + 360./gridOut%lonlat%dx_deg <= gridOut%lonlat%nx + 0.5) &
                                            & xNew = xNew + 360./gridOut%lonlat%dx_deg
       elseif(xNew > gridOut%lonlat%nx + 0.5)then
         if(xNew - 360./gridOut%lonlat%dx_deg >= 0.5) &
                                            & xNew = xNew - 360./gridOut%lonlat%dx_deg
       endif
          
     else ! Get georaphic  lon and lat
       if (gridIni%gridtype == anygrid) then 
            fLon = fu_2d_interpolation(pAnyGrdParam(gridIni%ag%indParam)%xc, &
                                             & x,y,gridIni%ag%nx,gridIni%ag%ny,method, notAllowed)
            fLat = fu_2d_interpolation(pAnyGrdParam(gridIni%ag%indParam)%yc, &
                                             & x,y,gridIni%ag%nx,gridIni%ag%ny,method, notAllowed)
       elseif (gridIni%gridtype == lonlat) then
            call modify_lonlat(gridIni%lonlat%sw_corner_modlat + (y-1) * gridIni%lonlat%dy_deg, &
                             & gridIni%lonlat%sw_corner_modlon + (x-1) * gridIni%lonlat%dx_deg, &
                             & gridIni%lonlat%pole, &
                             & pole_geographical, &
                             & fLat, &
                             & fLon)
       else
            call report (gridIni)
            call set_error('Non-supported input grid type','project_point_to_grid_xy')
       endif
       call  project_point_to_grid_lonlat(fLon, fLat, gridOut, xNew, yNew)
     endif

  end subroutine project_point_to_grid_xy


  !*****************************************************************

  subroutine project_point_to_grid_lonlat(fLon, fLat, gridOut, xNew, yNew)
    !
    ! The routine projects a point given in ABSOLUTE co-ordinates LON,LOT
    ! to the gridOut. Results are relative co-ordinates ixNew, iyNew.
    ! Since it may be used very many times, NO CHECKING is included.
    !
    implicit none

    ! Imported parameters with intent IN
    type(silja_grid), intent(in) :: gridOut
    real, intent(in) :: fLon, fLat

    ! Imported parameters with intent OUT
    real, intent (inout) :: xNew, yNew  !! Use as a hint for next search if possible

    ! Local variables
    integer :: ixC, iyC, ixD, iyD, nx, ny, itmp, jtmp, iAg
    real, dimension(:,:), pointer :: lats, lons, x3d, y3d, z3d
    real :: x0, x1, y0, y1, z0, z1, d
    real :: xp3d, yp3d, zp3d
    logical :: done, printit, ifGotBetter, ifBigError
    character(len = *), parameter :: sub_name = 'project_point_to_grid_lonlat'

    printit=.FALSE.
    if (.FALSE.) then
        printit=.true.
    endif


    select case(gridOut%gridtype)
      case(lonlat)
        call modify_lonlat(flat, &
                         & flon, &
                         & pole_geographical, &
                         & gridOut%lonlat%pole, &
                         & yNew, &
                         & xNew)
        xNew = 1+ (xNew - gridOut%lonlat%sw_corner_modlon) / gridOut%lonlat%dx_deg
        yNew = 1+ (yNew - gridOut%lonlat%sw_corner_modlat) / gridOut%lonlat%dy_deg
        
        ! If longitude out of grid, try rotate the earth once
        if(xNew < 0.5)then
          if(xNew + 360./gridOut%lonlat%dx_deg < gridOut%lonlat%nx + 0.5) &
                                               & xNew = xNew + 360./gridOut%lonlat%dx_deg
        elseif(xNew > gridOut%lonlat%nx + 0.5)then
          if(xNew - 360./gridOut%lonlat%dx_deg > 0.5) &
                                               & xNew = xNew - 360./gridOut%lonlat%dx_deg
        endif

      case (anygrid)
         nx = gridOut%ag%nx
         ny = gridOut%ag%ny
         iAg = gridOut%ag%indParam
         lons(1:nx,1:ny) => pAnyGrdParam(iAg)%xC(1:nx*ny)
         lats(1:nx,1:ny) => pAnyGrdParam(iAg)%yC(1:nx*ny)
         x3d(1:nx,1:ny) => pAnyGrdParam(iAg)%x3dc(1:nx*ny)
         y3d(1:nx,1:ny) => pAnyGrdParam(iAg)%y3dc(1:nx*ny)
         z3d(1:nx,1:ny) => pAnyGrdParam(iAg)%z3dc(1:nx*ny)

!         if (any(abs(cos(lons*degrees_to_radians)* cos(lats*degrees_to_radians) - x3d) > 1e-6 )) call ooops("X")
!         if (any(abs(sin(lons*degrees_to_radians)* cos(lats*degrees_to_radians) - y3d) > 1e-6 )) call ooops("Y")
!         if (any(abs(sin(lats*degrees_to_radians) - z3d) > 1e-6 )) call ooops("Y")

         xp3d = cos(flon*degrees_to_radians) * cos(flat*degrees_to_radians) 
         yp3d = sin(flon*degrees_to_radians) * cos(flat*degrees_to_radians) 
         zp3d = sin(flat*degrees_to_radians)

         
         ! Get initial guess: center, corners, previous value
         iXc = int_missing
         iYc = int_missing
         d   = real_missing !!Square euclidian distance in units of earth_raduis^2
         
          
         if (xNew > 0.5 .and. xNew < nx + 0.5 .and. yNew > 0.5 .and. yNew < ny + 0.5 ) then
           ! try previous
           call try(x3d,y3d,z3d, xp3d,yp3d,zp3d, iXc, iYc, nint(xNew), nint(yNew), d, ifGotBetter)
           if (printit) call msg("prev "//fu_str(ixC)//" "//fu_str(iYc), lons(ixC, iYC), lats(ixC, iYC))
         endif
         if (abs(d) > 1e-2) then  !! Might be grave wrong, check the grid center and corners
           call try(x3d,y3d,z3d, xp3d,yp3d,zp3d, iXc, iYc, nx/2, ny/2, d, ifGotBetter)
           if (printit) call msg("C "//fu_str(ixC)//" "//fu_str(iYc), lons(ixC, iYC), lats(ixC, iYC))
           call try(x3d,y3d,z3d, xp3d,yp3d,zp3d, iXc, iYc,    1,    1, d, ifGotBetter)
             if (printit) call msg("c1 "//fu_str(ixC)//" "//fu_str(iYc), lons(ixC, iYC), lats(ixC, iYC))
           call try(x3d,y3d,z3d, xp3d,yp3d,zp3d, iXc, iYc,    1,   ny, d, ifGotBetter)
             if (printit) call msg("c2 "//fu_str(ixC)//" "//fu_str(iYc), lons(ixC, iYC), lats(ixC, iYC))
           call try(x3d,y3d,z3d, xp3d,yp3d,zp3d, iXc, iYc,   nx,    1, d, ifGotBetter)
             if (printit) call msg("c3 "//fu_str(ixC)//" "//fu_str(iYc), lons(ixC, iYC), lats(ixC, iYC))
           call try(x3d,y3d,z3d, xp3d,yp3d,zp3d, iXc, iYc,   nx,   ny, d, ifGotBetter)
             if (printit) call msg("c4 "//fu_str(ixC)//" "//fu_str(iYc), lons(ixC, iYC), lats(ixC, iYC))
         endif


         do ixD = 1, max(nx,ny)
            done=.true.!
            ! iTmp is a dummy index, just to ensure that we do not gt out of the gtid
            do itmp = iyC, ny-1  ! Up
               call try(x3d,y3d,z3d, xp3d,yp3d,zp3d, iXc, iYc, iXc, iYc+1, d, ifGotBetter)
               if (.not. ifGotBetter) exit
               if (printit) call msg("Up "//fu_str(ixC)//" "//fu_str(iYc), lons(ixC, iYC), lats(ixC, iYC))
               done=.false.
            enddo
            do itmp = iyC, 2, -1              !Up
               call try(x3d,y3d,z3d, xp3d,yp3d,zp3d, iXc, iYc, iXc, iYc-1, d, ifGotBetter)
               if (.not. ifGotBetter) exit
               if (printit) call msg("Do "//fu_str(ixC)//" "//fu_str(iYc), lons(ixC, iYC), lats(ixC, iYC))
               done=.false.
            enddo

            ! For lon-global grids, wrapping should be allowed, for non-wrapped grids modulo/metrics 
            ! will ensure that we do not stick out of the grid 
            ! iTmp is dummy index, just to ensure that we do not get stuck in the loop if something goes wrong
            do itmp = 1, nx-1  ! Right
               call try(x3d,y3d,z3d, xp3d,yp3d,zp3d, iXc, iYc, modulo(iXc, nx)+1, iYc, d, ifGotBetter)
               if (.not. ifGotBetter) exit
               if (printit) call msg("Ri "//fu_str(ixC)//" "//fu_str(iYc), lons(ixC, iYC), lats(ixC, iYC))
               done=.false.
            enddo
            do itmp = nx, 2, -1              !Left
               call try(x3d,y3d,z3d, xp3d,yp3d,zp3d, iXc, iYc, modulo(iXc-2,nx)+1, iYc, d, ifGotBetter)
               if (.not. ifGotBetter) exit
               if (printit) call msg("Le "//fu_str(ixC)//" "//fu_str(iYc), lons(ixC, iYC), lats(ixC, iYC))
               done=.false.
            enddo
            if (done) exit
         enddo

         if (ixD >= max(nx,ny)) then
            !! Something went severely wrong: normally much less iterations needed
            call set_error("Could not project a point to anygrid", sub_name)
         endif
         !! 

         !Find next nearest in each direction from the nearest one
         if (ixC == nx) then
             ixD=-1
         elseif (ixC == 1) then
            ixD=1
         else
           d     = d3(x3d(iXc+1, iYc),xp3d , y3d(iXc+1, iYc),yp3d ,z3d(iXc+1, iYc), zp3d)
           if (d > d3(x3d(iXc-1, iYc),xp3d , y3d(iXc-1, iYc),yp3d ,z3d(iXc-1, iYc), zp3d)) then
             ixD=-1
            else
             ixD=1
            endif
         endif 

         if (iyC == ny) then
            iyD=-1
         elseif (iyC == 1) then
            iyD=1
         else
           d     = d3(x3d(iXc, iYc+1),xp3d , y3d(iXc, iYc+1),yp3d ,z3d(iXc, iYc+1), zp3d)
           if (d > d3(x3d(iXc, iYc-1),xp3d , y3d(iXc, iYc-1),yp3d ,z3d(iXc, iYc-1), zp3d)) then
             iyD=-1
            else
             iyD=1
            endif
         endif

         !! Make inter/extrapolation
         !! r0 Vector from nearest gridpoint to the target
         x0 = xp3d - x3d(iXc, iYc) 
         y0 = yp3d - y3d(iXc, iYc) 
         z0 = zp3d - z3d(iXc, iYc) 

         ! vector from nearest to second nearest in x grid direction
         x1 =  x3d(iXc + ixD, iYc) - x3d(iXc, iYc)
         y1 =  y3d(iXc + ixD, iYc) - y3d(iXc, iYc)
         z1 =  z3d(iXc + ixD, iYc) - z3d(iXc, iYc)

         !nearest point + projection of r0 onto the vector
         xNew = ixC + ixD * (x0*x1+y0*y1+z0*z1)/(x1*x1+y1*y1+z1*z1) 

         ! vector from nearest to second nearest in y grid direction
         x1 = x3d(iXc, iYc + iyD) - x3d(iXc, iYc) 
         y1 = y3d(iXc, iYc + iyD) - y3d(iXc, iYc) 
         z1 = z3d(iXc, iYc + iyD) - z3d(iXc, iYc) 

         !nearest point + projection of r0 onto the vector
         yNew = iyC + iyD * (x0*x1+y0*y1+z0*z1)/(x1*x1+y1*y1+z1*z1) 

#ifdef DEBUG
         !!! interpolated lon and lat
         x1=fu_2d_interpolation (pAnyGrdParam(iAg)%xC, xNew,yNew, nx,ny, linear, nearestPoint)
         y1=fu_2d_interpolation (pAnyGrdParam(iAg)%yC, xNew,yNew, nx,ny, linear, nearestPoint)

         x0 = abs(x1-flon) 
         if (x0 > 355.) x0 = abs(modulo(x1-flon+180., 360.) - 180.)

         x0 = x0 * cos(flat*degrees_to_radians) 
         y0 = abs(y1-flat)     

         ! fu_2d_interpolation does not extrapolate beyond the ceners mesh,
         ! while we continue the edge gradient infinitely, so the ceck 
         ! should work only within the centers mesh

         ! out of grid: do not force big error  
         if ((xnew < 1.     .and. ixc == 1)  .or. &
           & (xNew > nx  .and. ixc ==nx)  .or. &
           & (ynew < 1.     .and. iyc == 1)  .or. &
           & (yNew > ny  .and. iyc == ny) ) then
           x0 = -1.
           y0 = -1.
         endif

         ifBigError =  .not. (x0 < eps_degrees .and. y0 < eps_degrees ) 
         if (ifBigError) then
           call msg_warning("Too large error in reprojection to anygrid", sub_name)
           printit = .true.
         endif

         if (printit) then
           call msg("lon,              lat", fLon, fLat)
           call msg("lon-interp, latinterp", x1, y1)

           call msg("absolute error error in m", x0*1.1e5, y0*1.1e5)
           
           iTmp = iXc + (iYc-1)*nx !! 1D grid index
           call msg("relative error in grid cells", & 
                & x0*1.1e5 / pAnyGrdParam(iAg)%dx(iTmp),  y0*1.1e5/ pAnyGrdParam(iAg)%dy(iTmp))
           call msg("")

           call msg("Found nearest cell ", ixC, iyC)
           call msg("xNew, yNew", xNew, yNew)

           iTmp = max(1,ixC-1)
           jTmp = min(nx,ixC+1)

           call msg ("lons (ixmin,ixmax):", iTmp, jTmp)
           if (iyC>1) call msg("lons-1:", lons(itmp:jtmp, iyC-1))
             call msg("lons 0:", lons(itmp:jtmp, iyC))
           if (iyC<ny)  call msg("lons+1:", lons(itmp:jtmp, iyC+1))
           call msg ("lats (ixmin,ixmax):", iTmp, jTmp)
           if (iyC>1)call msg("lats-1:", lats(itmp:jtmp, iyC-1))
           call msg("lats 0:", lats(itmp:jtmp, iyC))
           if (iyC<ny) call msg("lats+1:", lats(itmp:jtmp, iyC+1))

           call ooops("")
         endif
         if (ifBigError) then
           call set_error("Too large error in reprojection to anygrid", sub_name)
         endif
#endif

      case default
        call report (gridOut)
        call set_error('Non-supported grid type', sub_name)
     end select ! gridIni%grid_type

     contains
       real function d3(x1,x2,y1,y2,z1,z2) 
         implicit none
         real, intent(in) :: x1,x2,y1,y2,z1,z2
         d3 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2)
  
       end function d3

       subroutine try(x3d,y3d,z3d, xp3d,yp3d,zp3d, ixC, iyC, ixtry, iytry, cost, ifGotBetter)
         real, dimension(:,:), intent(in) :: x3d,y3d,z3d !!grid
         real, intent(in) :: xp3d,yp3d,zp3d  !! point
         integer, intent(in) :: ixtry, iytry  !! Trial
         integer, intent(inout) :: ixC, iyC   !! Adjust if got better
         real, intent(inout) :: cost
         logical, intent(out) :: ifGotBetter

         real :: prevcost

         prevcost = cost
         
         cost = d3(x3d(ixtry, iytry),xp3d , y3d(ixtry, iytry),yp3d ,z3d(ixtry, iytry), zp3d)
         ifGotBetter = (prevcost == real_missing .or. cost < prevcost) 

         if (ifGotBetter) then
            ixC = ixTry
            iyC = iyTry
            if (printit) call msg("Better cost:"//fu_str(ixTry)//" "//fu_str(iyTry), cost) 
         else
           if (printit) call msg("Worse cost:"//fu_str(ixTry)//" "//fu_str(iyTry), cost) 
           cost = prevcost
         endif
       end subroutine try
         


  end subroutine project_point_to_grid_lonlat


  !****************************************************************************
  
  subroutine make_minimal_grid(grid, xIndex, yIndex)
    !
    ! Makes the 3x3 grid with ~1km cell size centered around the given point
    !
    implicit none
    
    type(silja_grid), intent(inout) :: grid
    real, intent(in) :: xIndex, yIndex                    ! grid co-ordinates
    
    ! Stupidity check
    !
    if(.not. defined(grid))then 
      call set_error('Undefined template grid','make_minimal_grid')
      return
    endif
    !
    ! Action depends on the grid type
    !
    select case(grid%gridtype)
      case (lonlat)
        grid%lonlat%nx = 3
        grid%lonlat%ny = 3
        grid%lonlat%sw_corner_modlon = grid%lonlat%sw_corner_modlon + (xIndex-1) * &
                                                                    & grid%lonlat%dx_deg - 0.01
        grid%lonlat%sw_corner_modlat = grid%lonlat%sw_corner_modlat + (yIndex-1) * &
                                                                    & grid%lonlat%dy_deg - 0.01
        grid%lonlat%dx_deg = 0.01
        grid%lonlat%dy_deg = 0.01

      case (anygrid)
        call msg_warning('Anygrid cannot be legally reduced','make_minimal_grid')
      case default
        call report(grid)
        call set_error('Unknown grid type','make_minilmal_grid')
    end select

  end subroutine make_minimal_grid


  !****************************************************************************
  
  subroutine extend_grid_to_coordinates(grid, x, y)
    !
    ! Extends the grid so that it covers the RELATIVE co-ordinates x, y
    ! So, x and y can be e.g. (-15,-9) and the grid has to be extended to cover them
    ! Evidently, in the new grid the values of x and y may be different, so
    ! they are recalculated if needed
    !
    implicit none

    ! Imported parameters
    type(silja_grid), intent(inout) :: grid
    real, intent(inout) :: x, y

    ! Stupidity check
    !
    if(.not. defined(grid))then 
      call set_error('Undefined grid given','extend_grid_to_coordinates')
      return
    endif
    !
    ! Checking depends on the grid type
    !
    select case(grid%gridtype)
      case (lonlat)
        !
        ! Enlarging longitude
        !
        if(x > real(grid%lonlat%nx))then  ! x is too large
          grid%lonlat%nx = int(x) + 1
        elseif(x < 1)then   ! x is negative or near zero
          grid%lonlat%nx = grid%lonlat%nx + int(2. - x)
          grid%lonlat%sw_corner_modlon = grid%lonlat%sw_corner_modlon - &
                                       & int(2. - x) * grid%lonlat%dx_deg 
          x = x + int(2. - x)
        else  ! x is inside - do nothing
        endif
        !
        ! Enlarging latitude
        !
        if(y > real(grid%lonlat%ny))then  ! y is too large
          grid%lonlat%ny = int(y) + 1
        elseif(y < 1)then   ! y is negative or near zero
          grid%lonlat%ny = grid%lonlat%ny + int(2. - y)
          grid%lonlat%sw_corner_modlat = grid%lonlat%sw_corner_modlat - &
                                       & int(2. - y) * grid%lonlat%dy_deg 
          y = y + int(2. - y)
        else  ! y is inside - do nothing
        endif
      case(anygrid)
        !
        ! Anygrid cannot legally be extended, so black magic ?
        !
        call msg_warning('Anygrid cannot legally be extended','extend_grid_to_coordinates')

      case default
        call set_error('Unknown grid type','extend_grid_to_coordinates')
        return
    end select
  end subroutine extend_grid_to_coordinates


  !*****************************************************************

  subroutine reposition_global_grid_lon(grid, iShift_)
    !
    ! Changes the starting corner of the global grid so that it 
    ! satisfies the SILAM standards: corners are within the rectangular
    ! (-180,-90) - (180,90)
    !
    implicit none

    ! Imported parameters
    type(silja_grid), intent(inout) :: grid
    integer, intent(out), optional :: iShift_

    ! Local variables
    integer :: iShift, iTmp

    !
    ! First, stupidity check. 
    !
    if(present(iShift_)) iShift_ = 0.
    if(.not. defined(grid))then 
      call set_error('Undefined grid given','reposition_global_grid_lon')
      return
    endif
    !
    ! Grids can be (i) over-the-global (funny but ECMWF keeps making such grids)
    ! (ii) global, (iii) smaller-than-global. In cases (i) and (ii) we just ensure that the 
    ! starting longitude is less than a cell above the -180. In the case (iii) the grid must
    ! satisfy the longitude standard range (-180,180).
    ! This all certainly depends on the grid type.
    !
    select case(grid%gridtype)
      case (lonlat)
        if(grid%lonlat%dx_deg * grid%lonlat%nx > 360.0 - 0.1 * grid%lonlat%dx_deg )then
          !
          ! Global grid with error of less than 10% of cell size or over-the-global grid. 
          ! Just rotate it to get (-180,180)
          !
          
          iShift = nint((-180. - grid%lonlat%sw_corner_modlon) / grid%lonlat%dx_deg)
          grid%lonlat%sw_corner_modlon = grid%lonlat%sw_corner_modlon + iShift * grid%lonlat%dx_deg
          !
          ! Due to stupidity with numerics and inaccuracy of GRIB definitions, 
          ! have to explicitly re-check the position an adjust with +-1 dx
          ! 
          ! Allow error of 0.1% of resolution
          if(grid%lonlat%sw_corner_modlon + 0.001*grid%lonlat%dx_deg < -180.0001) then
            grid%lonlat%sw_corner_modlon = grid%lonlat%sw_corner_modlon + grid%lonlat%dx_deg
            iShift = iShift + 1
          elseif(grid%lonlat%sw_corner_modlon - grid%lonlat%dx_deg >= -180.0) then
            grid%lonlat%sw_corner_modlon = grid%lonlat%sw_corner_modlon - grid%lonlat%dx_deg
            iShift = iShift - 1
          endif

        else
          !
          ! Less-than-global grid. Rotation is impossible, corners must satisfy the range
          !
          if((grid%lonlat%sw_corner_modlon < -180.0 .and. &
            & .not. (grid%lonlat%sw_corner_modlon .eps. -180.0)) .or. &
            & (grid%lonlat%sw_corner_modlon + grid%lonlat%dx_deg * grid%lonlat%nx > 180.0 .and. &
            & .not.(grid%lonlat%sw_corner_modlon + grid%lonlat%dx_deg * grid%lonlat%nx .eps. 180.0)))then
            !
            ! Grid is not standard and non-global along the longitude. Set the error
            !
            call msg_warning('Cannot modify the follwoing non-global grid to meet SILAM standard', &
                           & 'reposition_global_grid_lon')
            call report(grid)
            call set_error('Cannot modify the non-global grid to meet SILAM standard', &
                         & 'reposition_global_grid_lon')
            return
          else
            !
            ! The grid is standard SILAM one, do nothing
            !
            iShift = 0
          endif

        endif  ! global, over-global or less-than-global grid

      case(anygrid)

        call msg_warning('repositioning not done for anygrids','reposition_global_grid_lon')
      
      case default
        call set_error('Unknown grid type','reposition_global_grid_lon')
        return
    end select

    if(present(iShift_)) iShift_ = iShift

!        !
!        ! First, ensure that the starting corner so that it is above -180.
!        !
!        iShift = int(max(0., (-180. - grid%lonlat%sw_corner_modlon) / grid%lonlat%dx_deg + 0.99999))
!        grid%lonlat%sw_corner_modlon = grid%lonlat%sw_corner_modlon + iShift * grid%lonlat%dx_deg
!        !
!        ! Now, ensure that the end-corner is below 180.
!        !
!        iShift = -int(max(0., (grid%lonlat%sw_corner_modlon + grid%lonlat%dx_deg * grid%lonlat%nx - 180.0) / &
!                             & grid%lonlat%dx_deg + 0.99999))
!        grid%lonlat%sw_corner_modlon = grid%lonlat%sw_corner_modlon + iShift * grid%lonlat%dx_deg

  end subroutine reposition_global_grid_lon



  !*****************************************************************

  subroutine reposition_global_grid_lat(grid, iShift_)
    !
    ! Changes the startingh corner of the global grid so that it 
    ! satisfies the SILAM standards: corners are within the rectangular
    ! (-180,-90) - (180,90)
    !
    implicit none

    ! Imported parameters
    type(silja_grid), intent(inout) :: grid
    integer, intent(out), optional :: iShift_

    ! Local variables
    integer :: iShift

    !
    ! First, stupidity check. 
    !
    if(present(iShift_)) iShift_ = 0.
    if(.not. defined(grid))then 
      call set_error('Undefined grid given','reposition_global_grid_lat')
      return
    endif
    !
    ! Grids can be (i) over-the-global (funny but ECMWF keeps making such grids)
    ! (ii) global, (iii) smaller-than-global. In cases (i) and (ii) we just ensure that the 
    ! starting longitude is less than a cell above the -180. In the case (iii) the grid must
    ! satisfy the longitude standard range (-180,180).
    ! This all certainly depends on the grid type.
    !
    select case(grid%gridtype)
      case (lonlat)
        if(abs(grid%lonlat%dy_deg * grid%lonlat%ny - 180.0) / grid%lonlat%dy_deg < 0.1 .or. &
         & grid%lonlat%dy_deg * grid%lonlat%ny > 180.0)then
          !
          ! Global grid with error of less than 10% of cell size or over-the-global grid. 
          ! Just rotate it to get (-90,90)
          !
          iShift = int((-90. - grid%lonlat%sw_corner_modlat) / grid%lonlat%dy_deg)
          grid%lonlat%sw_corner_modlat = grid%lonlat%sw_corner_modlat + iShift * grid%lonlat%dy_deg
          !
          ! Due to stupidity with numerics and inaccuracy of GRIB definitions, 
          ! have to explicitly re-check the position an adjust with +-1 dx
          !
          if(grid%lonlat%sw_corner_modlat < -90.0) then
            grid%lonlat%sw_corner_modlat = grid%lonlat%sw_corner_modlat + grid%lonlat%dy_deg
            iShift = iShift + 1
          elseif(grid%lonlat%sw_corner_modlat - grid%lonlat%dy_deg >= -90.0) then
            grid%lonlat%sw_corner_modlat = grid%lonlat%sw_corner_modlat - grid%lonlat%dy_deg
            iShift = iShift - 1
          endif

        else
          !
          ! Less-than-global grid. Rotation is impossible, corners must satisfy the range
          !
          if((grid%lonlat%sw_corner_modlat < -90.0 .and. &
            & .not. (grid%lonlat%sw_corner_modlat .eps. -90.0)) .or. &
            & (grid%lonlat%sw_corner_modlat + grid%lonlat%dy_deg * (grid%lonlat%ny-1) > 90.0 .and. &
            & .not.(grid%lonlat%sw_corner_modlat + grid%lonlat%dy_deg * (grid%lonlat%ny-1) .eps. 90.0)))then
            !
            ! Grid is not standard and non-global along the latitude. Set the error
            !
            call msg_warning('Cannot modify the follwoing non-global grid to meet SILAM standard', &
                           & 'reposition_global_grid_lat')
            call report(grid)
            call set_error('Cannot modify the non-global grid to meet SILAM standard', &
                         & 'reposition_global_grid_lat')
            return
          else
            !
            ! The grid is standard SILAM one, do nothing
            !
            iShift = 0
          endif

        endif  ! global, over-global or less-than-global grid

      case(anygrid)

        call msg_warning('repositioning not done for anygrids','reposition_global_grid_lon')

      case default
        call set_error('Unknown grid type','reposition_global_grid_lat')
        return
    end select

    if(present(iShift_)) iShift_ = iShift

!        iShift = int(max(0., (-90. - grid%lonlat%sw_corner_modlat) / grid%lonlat%dy_deg + 0.99999))
!        grid%lonlat%sw_corner_modlat = grid%lonlat%sw_corner_modlat + iShift * grid%lonlat%dy_deg
!        !
!        ! Now, ensure that the end-corner is below 90.
!        !
!        iShift = -int(max(0., (grid%lonlat%sw_corner_modlat + grid%lonlat%dy_deg * (grid%lonlat%ny-1) - 90.0) / &
!                             & grid%lonlat%dy_deg + 0.99999))
!        grid%lonlat%sw_corner_modlat = grid%lonlat%sw_corner_modlat + iShift * grid%lonlat%ny 

  end subroutine reposition_global_grid_lat


  !******************************************************************
  
  subroutine set_grid_from_lores_templ(gridTemplate, gridLores, resolution_factor, &
                                     & gridOut, iCell_size_factor)
    !
    ! Makes are grid "almost" exactly as the gridTemplate but with resolution
    ! better than gridLores by about resolution_factor
    ! Needed when we want to keep easy indexing of lagrangian particles in low-resolution
    ! grid but still have sufficiently high-resolution underlying grid to use interpolation
    ! structures instead of trigonometry
    !
    implicit none

    ! Imported parameters
    type(silja_grid), intent(in) :: gridTemplate, gridLores
    real, intent(in) :: resolution_factor
    type(silja_grid), intent(out) :: gridOut
    integer, intent(out) :: iCell_size_factor

    ! Local variables
    integer :: nx, ny
    
    ! Copy the template grid...
    gridOut = gridTemplate
    !
    ! ... and change its resolution
    !
    iCell_size_factor = min(10, max(1, nint(resolution_factor * sqrt(fu_cell_size(gridTemplate) / &
                                                                   & fu_cell_size(gridLores)))))
    call grid_dimensions(gridOut, nx, ny)
    call adjust_grid_dimensions(gridOut, nx * iCell_size_factor, ny * iCell_size_factor)

  end subroutine set_grid_from_lores_templ

  !*****************************************************************
  !*****************************************************************
  !
  !  Horizontal interpolation structure
  !
  !*****************************************************************
  !*****************************************************************

  !*********************************************************************
  !
  !  Horizontal interpolation structure starts
  !
  !*********************************************************************
  !*********************************************************************

  
  subroutine deallocate_horiz_interp_struct(interp_struct)
    implicit none
    type(THorizInterpStruct), intent(inout) :: interp_struct
    
    if (interp_struct%interp_type == int_missing .or. .not. associated(interp_struct%weight)) then
      call set_error('Structure seems undefined', 'deallocate_horiz_interp_struct')
      return
    end if
    
    deallocate(interp_struct%indX, interp_struct%indY, interp_struct%weight)
    if (associated(interp_struct%weight_x)) deallocate(interp_struct%weight_x)
    if (associated(interp_struct%weight_y)) deallocate(interp_struct%weight_y)
    if (associated(interp_struct%ifValid)) deallocate(interp_struct%ifValid)
    
    interp_struct = HorizInterpStruct_missing
    
  end subroutine deallocate_horiz_interp_struct


  !*************************************************************************
  
  function fu_horiz_interp_struct(gridFrom, gridTo, interpolation_method, ifRandomise, &
      & iAccuracy, iOutside, ifMakeRotation, ifCreate) result(interpStructPtr)
    !
    ! Searches for an interpolation structure answering the request, creates it if
    ! not found
    !
    implicit none 

    ! return value
    type(THorizInterpStruct), pointer :: interpStructPtr

    ! Imported parameters
    type(silja_grid), intent(in) :: gridFrom, gridTo
    integer, intent(in) :: interpolation_method
    logical, intent(in) :: ifRandomise
    integer, intent(in), optional :: iAccuracy, iOutside
    logical, intent(in), optional :: ifMakeRotation
    logical, intent(in), optional :: ifCreate ! Can be set to false to 
                      !forbid creation and to make the routine thread-safe

    ! Local variables
    integer :: nxTo,nyTo, nxFrom, nyFrom, ix,iy, ixTo, iyTo, status, ii, iXLoop, iYLoop, maxCount, &
             & nSmall, iSml, jSml, iCount, maxExpectedCount, myiOutside, nWrap, iWrap
    real :: fX, fY, zx, zy, zx1, zy1, zx2, zy2, zx3, zy3, zx4, zy4, xFrom, yFrom, fTmp, ratio, &
          & fSmall_1, fSmall_2, xL, xR, yL, yU, fX_ll, fx_ur, fY_ll, fY_ur, fX_save
    real, dimension(4) :: distances
    integer, dimension(:,:,:), pointer :: ixTmp, iyTmp
    real, dimension(:,:,:), pointer :: weightTmp, weight_X_Tmp, weight_Y_Tmp
    integer, dimension(:,:), pointer :: arCount
    integer, dimension(:), pointer :: pIX, pIY, iWork
    real, dimension(:), pointer :: pWeight, pWeight_X, pWeight_Y, pWeight_Xmom, pWeight_Ymom,fWork
    logical :: ifFound, ifSimpleReprojection, ifLonGlobalTo, ifLonGlobalFrom
    real :: lat, lon, lat2, lon2
    real, dimension(3,2) :: B1, B2
    real, dimension(3,3) :: rotation_3d, rotation_3d_2
    real, dimension(2,2) :: rot_tmp
    type(silam_any_grid_param), pointer :: anyGridParamTo, anyGridParamFrom 

    !!!! same-projection conservative remap
    real :: y0From, y0To, x0From, x0To,  yeps, xeps, dxfrom, dxto, dyfrom, dyto, dxseg, dyseg, offXSeg, offYSeg
    integer :: nSegY, nSegX, nsegXFrom, nsegYFrom, iXFrom, iYfrom
    real,dimension(:), pointer :: bndXfrom, bndXto, bndX, bndYfrom, bndYto, bndY
    integer,dimension(:), pointer :: pIXfrom, pIXTo, pIYfrom, pIYTo 


    type(rng_type) :: rng !Random-number generator

    call rng_init(rng, 0) !Default seed

    interpStructPtr => null() !!Safe option

    !
    ! Stupidity check
    !
    if(.not. (defined(gridFrom) .and. defined(gridTo)))then
      call set_error('Undefined grids given for interpolation structure search','fu_horiz_interp_struct')
      return
    endif



    if (gridFrom == gridTo .and. HorizInterpPool%nHIS > 0) then !1:1 interpolation: Return missing structure
       interpStructPtr => HorizInterpPool%pHIS(1)%ptr
       HorizInterpPool%refCount(1) = HorizInterpPool%refCount(1) + 1
       return
    endif


    if (present(iOutside)) then
        myiOutside = iOutside
    else
        myiOutside = nearestPoint
    endif
    !
    ! First, try to find the structure
    !
    do iCount = 1, HorizInterpPool%nHIS
      if(HorizInterpPool%pHIS(iCount)%ptr%interp_type == interpolation_method)then
        if(HorizInterpPool%pHIS(iCount)%ptr%gridFrom == gridFrom)then
          if(HorizInterpPool%pHIS(iCount)%ptr%gridTo == gridTo)then
            if(HorizInterpPool%pHIS(iCount)%ptr%iOutside == myiOutside)then
            interpStructPtr => HorizInterpPool%pHIS(iCount)%ptr
            HorizInterpPool%refCount(iCount) = HorizInterpPool%refCount(iCount) + 1
!               call msg('Request horiz structure #, ref count', icount, HorizInterpPool%refCount(icount))
            return
          endif
        endif
      endif
      endif
    end do

    if (present(ifCreate)) then
      if (.not. ifCreate) return
    endif

    !
    ! Thread-unsafe stuff below
    !

    if (HorizInterpPool%nHIS == 0) then 
        allocate(HorizInterpPool%pHIS(1)%ptr, stat = iCount)
        if (iCount == 0)  then !Allocation ok
          HorizInterpPool%pHIS(1)%ptr = HorizInterpStruct_missing
          HorizInterpPool%refCount(1) =  1 ! Set counter to something, so it
                                           !  can not be released
          HorizInterpPool%nHIS = 1 !!This should be set only after the
                                   !  structure is fully initialized
        else
          call set_error('Failed allocation of identity interpolation strcutre',&
                    & 'fu_horiz_interp_struct')
          return
        endif
        
    endif
    !
    ! Strcutre is not found. Create new one!
    !
    if(HorizInterpPool%nHIS == size(HorizInterpPool%pHIS))then
      call set_error('No more space for horizonztal interpolation structures','fu_horiz_interp_struct')
      return
    endif
    HorizInterpPool%nHIS = HorizInterpPool%nHIS + 1
    allocate(HorizInterpPool%pHIS(HorizInterpPool%nHIS)%ptr, stat = iCount)
    if(fu_fails(iCount == 0,'Failed allocation of next horizontal interpolation strcutre', &
                          & 'fu_horiz_interp_struct'))return
    interpStructPtr => HorizInterpPool%pHIS(HorizInterpPool%nHIS)%ptr
    interpStructPtr%ifInterpolate = .True.
    interpStructPtr%ifGridFromGlobal = fu_ifLonGlobal(gridFrom)
    interpStructPtr%gridFrom = gridFrom
    interpStructPtr%gridTo = gridTo
    interpStructPtr%interp_type = interpolation_method
    horizInterpPool%refCount(HorizInterpPool%nHIS) = 1

    !
    ! Get the dimensions and prepare the indices for actually involved area
    !
    call grid_dimensions(interpStructPtr%gridFrom,nxFrom,nyFrom)
    call grid_dimensions(interpStructPtr%gridTo,nxTo,nyTo)
    interpStructPtr%ixStartFrom = nxFrom
    interpStructPtr%iyStartFrom = nyFrom
    interpStructPtr%ixEndFrom = 1
    interpStructPtr%iyEndFrom = 1
    interpStructPtr%ixStartTo = nxTo
    interpStructPtr%iyStartTo = nyTo
    interpStructPtr%ixEndTo = 1
    interpStructPtr%iyEndTo = 1
    interpStructPtr%iOutside = myiOutside

    ifLonGlobalTo = fu_ifLonGlobal(interpStructPtr%gridTo)
    ifLonGlobalFrom = fu_ifLonGlobal(interpStructPtr%gridFrom)
    !
    ! Allocate the necessary pointers and make up the reprojection coefficients
    !
    if (myiOutside == setMissVal) then
      allocate(interpStructPtr%ifValid(nxTo,nyTo),  stat = ii)
      if(fu_fails(ii == 0,'Failed interpolation-structure arrays allocation','fu_horiz_interp_struct'))return
    else 
        interpStructPtr%ifValid => null()
    endif

    select case(interpStructPtr%interp_type)

      case(average, summation)
        !
        ! These are the flux-preserving reprojection ways.
        !
        ifSimpleReprojection = gridFrom%gridtype == lonlat .and. &
                             & gridTo%gridtype == lonlat
        if(ifSimpleReprojection) ifSimpleReprojection = ifSimpleReprojection .and. &
                                                      & gridFrom%lonlat%pole == gridTo%lonlat%pole
        !
        ! For lonlat grids with the same projection the reprojection can be simple
        !
       !!!        if(.false.)then
       if(ifSimpleReprojection)then  
          !
          ! Poles are the same, i.e. the grid cell borders are parallel. Then a simple algebra can 
          ! serve the reprojection
          !
          
          call msg('Simplified same-pole reprojection')

          dxfrom = gridFrom%lonlat%dx_deg
          dxto   =  gridTo%lonlat%dx_deg
          x0From = gridFrom%lonlat%sw_corner_modlon - 0.5*dxfrom
          x0To = gridTo%lonlat%sw_corner_modlon - 0.5*dxto
          xeps = 1e-2 * dxto !! Epsilon to assume cell boundaries same
          if ( ifLonGlobalFrom ) then
              nSegXFrom =  nxFrom*2+1 !!Let wrap
              if ( x0From > x0To) x0From = x0From - 360. !
          else
             nSegXFrom =  nxFrom + 1
          endif
          nSegX = nSegXFrom + nxTo + 1
           
          dyfrom = gridFrom%lonlat%dy_deg
          dyto   = gridTo%lonlat%dy_deg
          y0From = gridFrom%lonlat%sw_corner_modlat - 0.5*dyfrom
          y0To = gridTo%lonlat%sw_corner_modlat - 0.5*dyto
          yeps = 1e-2 * gridTo%lonlat%dy_deg !! Epsilon to assume boundaries same
          nSegY = nyFrom + nYto + 2  !! Max no of overlapping segments

          ! Mazimum number of "From" cells to contribute to a "To" cell
          maxCount = (int(gridTo%lonlat%dx_deg / gridFrom%lonlat%dx_deg) + 2) * &
                    & (int(gridTo%lonlat%dy_deg / gridFrom%lonlat%dy_deg) + 2)

          iWork =>  fu_work_int_array(2*nSegX + 2*nSegY + nxTo*nYTo)
          ! Four slices for indices
          pIXfrom => iWork(1:nSegX)
          pIXTo   => iWork(nSegX+1:2*nSegX) 
          pIYfrom => iWork(2*nSegX+1 : 2*nSegX + nSegY)
          pIYTo   => iWork(2*nSegX + nSegY+1 : 2*nSegX + 2*nSegY)
          arCount(1:nxTo,1:nYTo)  => iWork(2*nSegX + 2*nSegY+1 : 2*nSegX + 2*nSegY + nxTo*nYTo)

          ! 
          fWork => fu_work_array(3*nSegX + 3*nSegY)
          ! six scices for
          bndXfrom => fWork(1:nSegX)
          bndXto   => fWork(nSegX +1   : 2*nSegX)
          bndX     => fWork(2*nSegX +1 : 3*nSegX)
          bndYfrom => fWork(3*nSegX +1 : 3*nSegX + nSegY)
          bndYto    => fWork(3*nSegX + nSegY +1 : 3*nSegX + 2*nSegY)
          bndY     => fWork(3*nSegX + 2*nSegY +1 : 3*nSegX + 3*nSegY)

          do ix = 0, nsegXFrom-1
            bndXfrom(ix+1) =  x0From + ix * dxFrom 
          enddo
          do ix = 0, nxTo
            bndXTo(ix+1) =  x0To + ix * dxTo 
          enddo

          do iy = 0, nyFrom
            bndYfrom(iy+1) =  y0From + iy * dyFrom
          enddo
          do iy = 0, nyTo
            bndYTo(iy+1) =  y0To + iy * dyTo
          enddo

          call remapcon_1D(bndXfrom(1:nsegXFrom), bndXto(1:nxTo+1), xeps, bndX, pIXfrom, pIXTo, nSegX)
          call remapcon_1D(bndYfrom(1:nyFrom+1), bndYto(1:nyTo+1), yeps, bndY, pIYfrom, pIYTo, nSegY)

          interpStructPtr%nCoefs = maxCount
          allocate(interpStructPtr%indX(maxCount,nxTo,nyTo), interpStructPtr%indY(maxCount,nxTo,nyTo), &
                 & interpStructPtr%weight(maxCount,nxTo,nyTo), &
                 & interpStructPtr%weight_X(maxCount,nxTo,nyTo), &
                 & interpStructPtr%weight_Y(maxCount,nxTo,nyTo), stat = ii)
          if(fu_fails(ii == 0,'Failed interpolation-structure arrays allocation','fu_horiz_interp_struct'))return

          interpStructPtr%indX(:,:,:) = 0
          interpStructPtr%indY(:,:,:) = 0
          interpStructPtr%weight(:,:,:) = 0.
          interpStructPtr%weight_X(:,:,:) = 0.
          interpStructPtr%weight_Y(:,:,:) = 0.
          arCount(:,:) = 0


          do iy = 1, nSegY
            iYTo = pIYTo(iy)
            iYFrom = pIYFrom(iy)
            if (iYTo < 1) cycle
            if (iYTo > nYto) cycle
            if (iYFrom < 1 .or. iYFrom > nYFrom ) then
              if (myioutside == setZero) cycle !!! Just assume all zeros outside the grid_from
              !! I could not figure out any consistent formulation for conservative reprojection 
              !! for nearestPoint or setMissVal.   RK
              call msg("grid from and to:")
              call report(gridFrom)
              call report(gridTo)
              call set_error("Cannot map Y conservatively", "fu_horiz_interp_struct")
              return
            endif
            dySeg = bndY(iY+1) - bndY(iY)

            !!Offset of center of segment from the trget cell (in units of target cell)
            offYSeg = 0.5* (bndY(iY+1) + bndY(iY+1) - bndYto(iYTo) - bndYto(iYTo+1)) / dyTo

            do ix = 1, nSegX
              iXTo = pIXTo(iX)
              if (iXTo < 1) cycle
              if (iXTo > nXto) cycle
              iXFrom = pIXFrom(iX)
              if (ifLonGlobalFrom) iXFrom = modulo(iXFrom-1, nxFrom) + 1
              if (iXFrom < 1 .or. iXFrom > nXFrom ) then
                if (myioutside == setZero) cycle  !!! Just assume all zeros outside the grid_from
                call msg("grid from and to:")
                call report(gridFrom)
                call report(gridTo)
                call set_error("Cannot map X conservatively", "fu_horiz_interp_struct")
                return
              endif

              dXSeg = bndX(iX+1) - bndX(iX)
              !!Offset of center of segment from the trget cell (in units of target cell)
              offXSeg = 0.5* (bndX(iX+1) + bndX(iX+1) - bndXto(iXTo) - bndXto(iXTo+1)) / dXTo

              arCount(iXTo, iYTo) = arCount(iXTo, iYTo) + 1
              iCount = arCount(iXTo, iYTo)

              interpStructPtr%indX(iCount,iXto,iYto) = iXFrom
              interpStructPtr%indY(iCount,iXto,iYto) = iYfrom
              if(interpStructPtr%interp_type == average)then
                fTmp = dXSeg * dySeg / (dxTo * dyTo) !!! Segment size as a fraction of output -- weights sum to one
              else
                fTmp = dXSeg * dySeg / (dxFrom * dyFrom) !!! Segment size as a fraction of input cell -- weights sum 
                                                         !! to ratio of cell areas 
              endif
              interpStructPtr%weight(iCount,iXto,iYto) =  fTmp
              interpStructPtr%weight_X(iCount,iXto,iYto) = fTmp * offXSeg
              interpStructPtr%weight_Y(iCount,iXto,iYto) = fTmp * offYSeg

              ! update the range of useful From and To cells 
              ! Only if not yet updated
              if (interpStructPtr%iyEndTo == 1) then !! Sic!  Y index is used for check
                interpStructPtr%ixStartTo = min(interpStructPtr%ixStartTo, ixTo)
                interpStructPtr%ixEndTo = max(interpStructPtr%ixEndTo, ixTo)
                interpStructPtr%ixStartFrom = min(interpStructPtr%ixStartFrom, ixFrom)
                interpStructPtr%ixEndFrom = max(interpStructPtr%ixEndFrom, ixFrom)
              endif
            enddo
            interpStructPtr%iyStartTo = min(interpStructPtr%iyStartTo, iYTo)
            interpStructPtr%iyEndTo = max(interpStructPtr%iyEndTo, iYTo)
            interpStructPtr%iyStartFrom = min(interpStructPtr%iyStartFrom, iYFrom)
            interpStructPtr%iyEndFrom = max(interpStructPtr%iyEndFrom, iYFrom)
          enddo
          call free_work_array(iWork)
          call free_work_array(fWork)
          
        else ! not same poles
          ! Note that the number of cells affecting the single output one is unknown. Have to use 
          ! temporaries as long as needed
          !
          ! Poles are different. No other way, have to go for tough work.
          ! Select the number of small cells, which the original grid cell will be broken into
          !
call msg('Full-scale different-poles reprojection')
        ! Let's get temporary array that will almost surely incorporate the neede coefs
        !
        ! FIXME black magic with ratios sometimes causes OOM errors
        ratio = fu_cell_size_midpoint(interpStructPtr%gridFrom) / &
              & fu_cell_size_midpoint(interpStructPtr%gridTo)
        
        call msg('Allocate average/sum horiz structure: #, resoltuion ratio:', HorizInterpPool%nHIS, ratio)
!          call msg('GridFrom:')
!          call report(gridFrom)
!          call msg('gridTo:')
!          call report(gridTo)
        if(ratio > 100. .or. ratio < 0.01)then
          call msg('Extreme ratio of resolution of gridFrom and gridTo:', ratio)
          call msg('GridFrom:')
          call report(gridFrom)
          call msg('gridTo:')
          call report(gridTo)
          call msg_warning('Extreme ratio of grid resolutions but continue','fu_horiz_interp_struct')
        endif
        
        maxExpectedCount = int(5.0 * ratio + 5.0)
        pIX => fu_work_int_array(maxExpectedCount * nxTo * nyTo)
        pIY => fu_work_int_array(maxExpectedCount * nxTo * nyTo)
        pWeight => fu_work_array(maxExpectedCount * nxTo * nyTo)
        pWeight_X => fu_work_array(maxExpectedCount * nxTo * nyTo)
        pWeight_Y => fu_work_array(maxExpectedCount * nxTo * nyTo)
        iWork => fu_work_int_array(nxTo * nyTo)
        
        if(maxExpectedCount * nxTo * nyTo < size(pIX)) maxExpectedCount = int(size(pIX) / (nxTo * nyTo))
call msg('maxExpectedCount =', maxExpectedCount)

        ixTmp(1:maxExpectedCount,1:nxTo,1:nyTo) => pIX(1:maxExpectedCount * nxTo * nyTo)
        iyTmp(1:maxExpectedCount,1:nxTo,1:nyTo) => pIY(1:maxExpectedCount * nxTo * nyTo)
        weightTmp(1:maxExpectedCount,1:nxTo,1:nyTo) => pWeight(1:maxExpectedCount * nxTo * nyTo)
        weight_X_Tmp(1:maxExpectedCount,1:nxTo,1:nyTo) => pWeight_X(1:maxExpectedCount * nxTo * nyTo)
        weight_Y_Tmp(1:maxExpectedCount,1:nxTo,1:nyTo) => pWeight_Y(1:maxExpectedCount * nxTo * nyTo)
        arCount(1:nxTo,1:nyTo) => iWork(1:nxTo * nyTo)

        if(error)return
        arCount = 0
        maxCount = 0
        ixTmp = 0
        iyTmp = 0
        weightTmp = 0.0
        weight_X_Tmp = 0.0
        weight_Y_Tmp = 0.0

          nSmall = int(sqrt(2000. * ratio)/2.) * 2 + 1

          nSmall = min(nSmall, int(sqrt(5.e4/nxFrom*nyFrom))+1) ! nSmall*nSmall*nx*ny < smth_big

          nSmall = min(1001, max(nSmall,9)) ! Hard limits

          if(present(iAccuracy)) nSmall = nSmall * iAccuracy / 10
          nSmall = max(nSmall, 3)  ! For reduced accuracy this will be the hard limit
        
          fSmall_1 = 1.0 / real(nSmall)
          if(interpStructPtr%interp_type == average)then
            fSmall_2 = 1.0
          else
            fSmall_2 = 1.0 / real(nSmall*nSmall)
          endif
  call msg('nSmall,,nxFrom,nyFrom,nxTo,nyTo =',(/nSmall,nxFrom,nyFrom,nxTo,nyTo/))
          !
          ! Actual work starts
          !
          fX = real_missing
          fY = real_missing
          do iy = 1, nyFrom
            do ix = 1, nxFrom
              !
              ! Check that corners of the cell are covered by gridTo: do we need to bother?
              !
              call project_point_to_grid(interpStructPtr%gridFrom,real(ix)-0.5,real(iy)-0.5, &
                                       & interpStructPtr%gridTo,fX,fY)
              if(fX<0.5 .or. fX>nxTo+0.5 .or. fY<0.5 .or. fY>nyTo+0.5)then
                call project_point_to_grid(interpStructPtr%gridFrom,real(ix)+0.5,real(iy)-0.5, &
                                         & interpStructPtr%gridTo,fX,fY)
                if(fX<0.5 .or. fX>nxTo+0.5 .or. fY<0.5 .or. fY>nyTo+0.5)then
                  call project_point_to_grid(interpStructPtr%gridFrom,real(ix)-0.5,real(iy)+0.5, &
                                           & interpStructPtr%gridTo,fX,fY)
                  if(fX<0.5 .or. fX>nxTo+0.5 .or. fY<0.5 .or. fY>nyTo+0.5)then
                    call project_point_to_grid(interpStructPtr%gridFrom,real(ix)+0.5,real(iy)+0.5, &
                                             & interpStructPtr%gridTo,fX,fY)
                    if(fX<0.5 .or. fX>nxTo+0.5 .or. fY<0.5 .or. fY>nyTo+0.5)cycle
                  endif
                endif
              endif    ! if all corners are outside the gridTo
              !
              ! OK, have to do the work: useful From-cell
              !
              interpStructPtr%ixStartFrom = min(interpStructPtr%ixStartFrom, ix)
              interpStructPtr%iyStartFrom = min(interpStructPtr%iyStartFrom, iy)
              interpStructPtr%ixEndFrom = max(interpStructPtr%ixEndFrom, ix)
              interpStructPtr%iyEndFrom = max(interpStructPtr%iyEndFrom, iy)
              do jSml = 1, nSmall
                do iSml = 1, nSmall
  if(ifRANDOMISE)then
                  call project_point_to_grid(interpStructPtr%gridFrom, &
                                           & ix-0.5+(iSml-0.5)*fSmall_1+fu_random_number_center_rng(rng, 0.,fSmall_1), &
                                           & iy-0.5+(jSml-0.5)*fSmall_1+fu_random_number_center_rng(rng, 0.,fSmall_1), &
                                           & interpStructPtr%gridTo, &
                                           & fX, fY)
  else
                  call project_point_to_grid(interpStructPtr%gridFrom, &
                                           & ix-0.5+(iSml-0.5)*fSmall_1, &
                                           & iy-0.5+(jSml-0.5)*fSmall_1, &
                                           & interpStructPtr%gridTo, &
                                           & fX, fY)
  endif
                  if(fX>=0.5 .and. fX<nxTo+0.5 .and. fY>=0.5 .and. fY<nyTo+0.5)then
                    !
                    ! Find out if this cellFrom has already contributed to this cellTo. Add weights 
                    ! if yes, make new if no
                    !
                    ixTo = nint(fX)
                    iyTo = nint(fY)
                    interpStructPtr%ixStartTo = min(ixTo, interpStructPtr%ixStartTo)
                    interpStructPtr%iyStartTo = min(iyTo, interpStructPtr%iyStartTo)
                    interpStructPtr%ixEndTo = max(ixTo, interpStructPtr%ixEndTo)
                    interpStructPtr%iyEndTo = max(iyTo, interpStructPtr%iyEndTo)
                    ifFound = .false.

  !if(ixTo == 41 .and. iyTo ==46)then
  !  call msg('ixFrom, iyFrom, m, fX, fY, pX, pY: (' + fu_str(ix) + ',' + fu_str(iy) + ')', &
  !                  & (/fSmall_2, fX, fY, fSmall_2 * (fX-ixTo), fSmall_2 * (fY-iyTo)/))
  !endif
                    do iCount = arCount(ixTo,iyTo), 1, -1
                      if(ixTmp(iCount,ixTo,iyTo) == ix .and. iyTmp(iCount,ixTo,iyTo) == iy)then
                        weightTmp(iCount,ixTo,iyTo) = weightTmp(iCount,ixTo,iyTo) + fSmall_2
                        weight_X_Tmp(iCount,ixTo,iyTo) = weight_X_Tmp(iCount,ixTo,iyTo) + &
                                                                                & fSmall_2 * (fX-ixTo)
                        weight_Y_Tmp(iCount,ixTo,iyTo) = weight_Y_Tmp(iCount,ixTo,iyTo) + &
                                                                                & fSmall_2 * (fY-iyTo)
                        ifFound = .true.
                        exit
                      endif
                    end do
                    if(ifFound)cycle  ! Finished with this small cell
                    !
                    ! This cellFrom has not yet contributed to the cellTo. Take the new element, 
                    ! beware on limited array size
                    !
                    if(arCount(ixTo,iyTo) == maxExpectedCount)then
  !                  if(arCount(ixTo,iyTo) == size(ixTmp,1))then
                      ii = max(10,(maxExpectedCount * 5) / 4) ! Enlarge by 20% or by 10
  call msg('Enlarging to/from:', ii, maxExpectedCount)
                      call enlarge_int_counter(ixTmp, pIX, ii, maxExpectedCount, nxTo, nyTo)
                      call enlarge_int_counter(iyTmp, pIY, ii, maxExpectedCount, nxTo, nyTo)
                      call enlarge_real_counter(weightTmp, pWeight, ii, maxExpectedCount, nxTo, nyTo)
                      call enlarge_real_counter(weight_X_Tmp, pWeight_X, ii, maxExpectedCount, nxTo, nyTo)
                      call enlarge_real_counter(weight_Y_Tmp, pWeight_Y, ii, maxExpectedCount, nxTo, nyTo)
                      maxExpectedCount = ii
                      if(error)return
                    endif
                    arCount(ixTo,iyTo) = arCount(ixTo,iyTo) + 1
                    maxCount = max(maxCount,arCount(ixTo,iyTo))
                    ixTmp(arCount(ixTo,iyTo),ixTo,iyTo) = ix    ! from
                    iyTmp(arCount(ixTo,iyTo),ixTo,iyTo) = iy    ! from
                    weightTmp(arCount(ixTo,iyTo),ixTo,iyTo) = fSmall_2
                    weight_X_Tmp(arCount(ixTo,iyTo),ixTo,iyTo) = fSmall_2 * (fX-ixTo)
                    weight_Y_Tmp(arCount(ixTo,iyTo),ixTo,iyTo) = fSmall_2 * (fY-iyTo)
                  endif
                end do  ! iSml
              end do  ! jSml
            end do   ! ix
          end do  ! iy
          !
          ! OK, we know the sizes and temporary arrays have been filled in. Copy stuff to the structure
          !
          interpStructPtr%nCoefs = maxCount
          allocate(interpStructPtr%indX(maxCount,nxTo,nyTo), interpStructPtr%indY(maxCount,nxTo,nyTo), &
                 & interpStructPtr%weight(maxCount,nxTo,nyTo), &
                 & interpStructPtr%weight_X(maxCount,nxTo,nyTo), &
                 & interpStructPtr%weight_Y(maxCount,nxTo,nyTo), stat = ii)
          if(fu_fails(ii == 0,'Failed interpolation-structure arrays allocation','fu_horiz_interp_struct'))return
          do iy = 1, nyTo
            do ix = 1, nxTo
              fX = sum(weightTmp(1:maxCount,ix,iy))
              if(fX > 0.0)then
                if(interpStructPtr%interp_type == average)then
                  do iCount = 1, maxCount
                    interpStructPtr%indX(iCount,ix,iy) = ixTmp(iCount,ix,iy)
                    interpStructPtr%indY(iCount,ix,iy) = iyTmp(iCount,ix,iy)
                    interpStructPtr%weight(iCount,ix,iy) = weightTmp(iCount,ix,iy) / fX
                    interpStructPtr%weight_X(iCount,ix,iy) = weight_X_Tmp(iCount,ix,iy) / fX
                    interpStructPtr%weight_Y(iCount,ix,iy) = weight_Y_Tmp(iCount,ix,iy) / fX
                  end do
                else
                  do iCount = 1, maxCount
  !if(ix == 41 .and. iy ==46)then
  !  call msg('ixFrom, iyFrom, weight, weight_X, weight_Y: (' + &
  !                & fu_str(ixTmp(iCount,ix,iy)) + ',' + fu_str(iyTmp(iCount,ix,iy)) + ')', &
  !                & (/weightTmp(iCount,ix,iy), weight_X_Tmp(iCount,ix,iy), weight_Y_Tmp(iCount,ix,iy)/))
  !endif
                    interpStructPtr%indX(iCount,ix,iy) = ixTmp(iCount,ix,iy)
                    interpStructPtr%indY(iCount,ix,iy) = iyTmp(iCount,ix,iy)
                    interpStructPtr%weight(iCount,ix,iy) = weightTmp(iCount,ix,iy)
                    interpStructPtr%weight_X(iCount,ix,iy) = weight_X_Tmp(iCount,ix,iy)
                    interpStructPtr%weight_Y(iCount,ix,iy) = weight_Y_Tmp(iCount,ix,iy)
                  end do
                endif  ! aver or sum
              else
                interpStructPtr%indX(1:maxCount,ix,iy) = 0
                interpStructPtr%indY(1:maxCount,ix,iy) = 0
                interpStructPtr%weight(1:maxCount,ix,iy) = 0.
                interpStructPtr%weight_X(1:maxCount,ix,iy) = 0.
                interpStructPtr%weight_Y(1:maxCount,ix,iy) = 0.
              endif  ! sum(1:maxCount)>0
            end do   ! ix
          end do  ! iy
          call free_work_array(iWork)
          call free_work_array(pIX)
          call free_work_array(pIY)
          call free_work_array(pWeight)
          call free_work_array(pWeight_X)
          call free_work_array(pWeight_Y)
  !        deallocate(ixTmp, iyTmp, weightTmp, weight_X_Tmp, weight_Y_Tmp)
        endif  ! if poles are same
          !
          ! We might need to deal with pieces of the output grid outside the input grid
          !
          select case(myiOutside)
              case(notAllowed, handleGlobalGrid)
                if(any(arCount(1:nxTo,1:nyTo) == 0))then
                      call set_error('Out-of-grid-from point in grid-to','fu_horiz_interp_struct')

                      call msg("arCount>0:")
                      do iy = nyTo , 1, -1 !Top-bottom, so map looks normal
                        call msg("",arCount(1:nxTo,iy) > 0)
                      enddo

                      return
                    endif
                
              case(nearestPoint)  ! find the first closest reasonable point and copy all its values
                call msg("Switch: iOutside = nearest")
                do iy = 1, nyTo
                  do ix = 1, nxTo
                    if(arCount(ix,iy) == 0)then
                      patch: do iCount = 1, max(nxTo, nyTo)
                        do jSml = max(iy-iCount, 1), min(iy+iCount, nyTo)
                          do iSml = max(ix-iCount, 1), min(ix+iCount, nxTo)
                            if(arCount(iSml,jSml) > 0)then
                              interpStructPtr%indX(1:maxCount,ix,iy) = interpStructPtr%indX(1:maxCount,iSml,jSml)
                              interpStructPtr%indY(1:maxCount,ix,iy) = interpStructPtr%indY(1:maxCount,iSml,jSml)
                              interpStructPtr%weight(1:maxCount,ix,iy) = interpStructPtr%weight(1:maxCount,iSml,jSml)
                              interpStructPtr%weight_X(1:maxCount,ix,iy) = interpStructPtr%weight_X(1:maxCount,iSml,jSml)
                              interpStructPtr%weight_Y(1:maxCount,ix,iy) = interpStructPtr%weight_Y(1:maxCount,iSml,jSml)
                              exit patch
                            endif
                          end do  ! iSml
                        end do  ! jSml
                      end do patch
                    endif  ! no contribution from input grid
                  end do  ! ix
                end do  ! iy
                
              case(setZero)  ! do nothing: zeroes are already there
                 call msg("Switch: iOutside = setZero")
              case(setMissVal)
                 call msg("Switch: iOutside = setMissVal")
                 where(arCount > 0) 
                      interpStructPtr%ifValid=.true.
                 elsewhere
                      interpStructPtr%weight(1,:,:) = real_missing
                      interpStructPtr%weight_X(1,:,:) = real_missing
                      interpStructPtr%weight_Y(1,:,:) = real_missing
                      interpStructPtr%ifValid=.false.
                 endwhere
              case default
                call set_error('Strange out-of-grid handle switch:'+fu_str(myiOutside),'fu_horiz_interp_struct')
          end select   ! iOutside

call msg('Finished averaging/sum structure')

      case(nearest_point)
        !
        ! coefs vectors are 1-element pointing to the nearest point. Out-of-grid is non-existent here
        !
        call msg('Allocate nearest-point hor structure: #:', HorizInterpPool%nHIS)
        interpStructPtr%nCoefs = 1
        allocate(interpStructPtr%indX(interpStructPtr%nCoefs,nxTo,nyTo), &
               & interpStructPtr%indY(interpStructPtr%nCoefs,nxTo,nyTo), &
               & interpStructPtr%weight(interpStructPtr%nCoefs,nxTo,nyTo), stat=status)
        if(status /= 0)then
          call msg('Coefs allocation status:',status)
          call set_error('Failed to allocate memory for coefs','fu_horiz_interp_struct')
          return
        endif
        nullify(interpStructPtr%ifValid) 
        nullify(interpStructPtr%weight_X)
        nullify(interpStructPtr%weight_Y)
        interpStructPtr%ixStartTo = 1   ! nearest point always defines the whole output grid
        interpStructPtr%iyStartTo = 1
        interpStructPtr%ixEndTo = nxTo
        interpStructPtr%iyEndTo = nyTo

        interpStructPtr%ixStartFrom = nxFrom
        interpStructPtr%iyStartFrom = nyFrom
        interpStructPtr%ixEndFrom = 1
        interpStructPtr%iyEndFrom = 1

        fX = real_missing
        fY = real_missing
        do iy=1,nyTo
          do ix=1,nxTo
            !
            ! A trick: back-projection gives the co-ordinates in the from-grid
            ! FIXME spherical distances in lat-lon must be used!!! (Roux)
            ! FIXME global grids should be handled hrere
            call project_point_to_grid(interpStructPtr%gridTo,real(ix),real(iy), &
                                     & interpStructPtr%gridFrom,fX,fY)
            
            zx = fX - REAL(INT(fX)) ! Exceedance of fX over lower border
            zy = fY - REAL(INT(fY))

            distances(1) = zx*zx + zy*zy
            distances(2) = (1.-zx)*(1.-zx) + zy*zy
            distances(3) = zx*zx + (1.-zy)*(1.-zy)
            distances(4) = (1.-zx)*(1.-zx) + (1.-zy)*(1.-zy)

            iCount = MINLOC(distances, 1)

            select case(iCount)
              case(1)
                interpStructPtr%indX(1,ix,iy) = min(nxFrom,max(1,int(fX)))
                interpStructPtr%indY(1,ix,iy) = min(nyFrom,max(1,int(fY)))
              case(2)
                interpStructPtr%indX(1,ix,iy) = min(nxFrom,max(1,int(fX)+1))
                interpStructPtr%indY(1,ix,iy) = min(nyFrom,max(1,int(fY)))
              case(3)
                interpStructPtr%indX(1,ix,iy) = min(nxFrom,max(1,int(fX)))
                interpStructPtr%indY(1,ix,iy) = min(nyFrom,max(1,int(fY)+1))
              case(4)
                interpStructPtr%indX(1,ix,iy) = min(nxFrom,max(1,int(fX)+1))
                interpStructPtr%indY(1,ix,iy) = min(nyFrom,max(1,int(fY)+1))
              case default
                call msg('Strange MINLOC output',iCount)
                call set_error('Strange MINLOC output','fu_horiz_interp_struct')
                return
            end select
            interpStructPtr%ixStartFrom = min(interpStructPtr%ixStartFrom, interpStructPtr%indX(1,ix,iy))
            interpStructPtr%iyStartFrom = min(interpStructPtr%iyStartFrom, interpStructPtr%indY(1,ix,iy))
            interpStructPtr%ixEndFrom = max(interpStructPtr%ixEndFrom, interpStructPtr%indX(1,ix,iy))
            interpStructPtr%iyEndFrom = max(interpStructPtr%iyEndFrom, interpStructPtr%indY(1,ix,iy))
            interpStructPtr%weight(1,ix,iy) = 1. ! Always
            if (myiOutside == setMissVal) then
                 if ((fX>1 .and. fX<nxFrom .and.  fY>1 .and. fY<nYfrom)) then
                    interpStructPtr%ifValid(ix,iy) = .True.
                 else
                    interpStructPtr%ifValid(ix,iy) = .False.
                    interpStructPtr%weight(1,ix,iy) = real_missing
                 endif
            elseif (myiOutside == setZero) then
                if (.not. (fX>1 .and. fX<nxFrom .and.  fY>1 .and. fY<nYfrom)) then
                    interpStructPtr%weight(1,ix,iy) = 0. 
                endif
            endif
          end do
        end do

      case(linear)
        !
        ! Four surrounding grid elements
        !
        call msg('Allocate linear hor structure: #:', HorizInterpPool%nHIS)
        interpStructPtr%nCoefs = 4
        allocate(interpStructPtr%indX(interpStructPtr%nCoefs,nxTo,nyTo), &
               & interpStructPtr%indY(interpStructPtr%nCoefs,nxTo,nyTo), &
               & interpStructPtr%weight(interpStructPtr%nCoefs,nxTo,nyTo), stat=status)
        if(status /= 0)then
          call msg('Coefs allocation status:',status)
          call set_error('Failed to allocate memory for coefs','fu_horiz_interp_struct')
          return
        endif
        nullify(interpStructPtr%weight_X)
        nullify(interpStructPtr%weight_Y)
        interpStructPtr%ixStartTo = 1
        interpStructPtr%iyStartTo = 1
        interpStructPtr%ixEndTo = nxTo
        interpStructPtr%iyEndTo = nyTo


        fX = real_missing
        fY = real_missing
        do iy=1,nyTo
          do ix=1,nxTo
            !
            ! A trick: back-projection gives the co-ordinates in the from-grid
            !
            call project_point_to_grid(interpStructPtr%gridTo,real(ix),real(iy), &
                                     & interpStructPtr%gridFrom,fX,fY)
            
            zx = fX - REAL(INT(fX)) ! Exceedance of fX over lower border
            zy = fY - REAL(INT(fY))

            interpStructPtr%indX(1,ix,iy) = min(nxFrom,max(1,int(fX)))
            interpStructPtr%indY(1,ix,iy) = min(nyFrom,max(1,int(fY)))
            interpStructPtr%weight(1,ix,iy) = (1. - zx) * (1. - zy)

            interpStructPtr%indX(2,ix,iy) = min(nxFrom,max(1,int(fX)+1))
            interpStructPtr%indY(2,ix,iy) = min(nyFrom,max(1,int(fY)))
            interpStructPtr%weight(2,ix,iy) = zx * (1. - zy)

            interpStructPtr%indX(3,ix,iy) = min(nxFrom,max(1,int(fX)))
            interpStructPtr%indY(3,ix,iy) = min(nyFrom,max(1,int(fY)+1))
            interpStructPtr%weight(3,ix,iy) = (1. - zx) * zy

            interpStructPtr%indX(4,ix,iy) = min(nxFrom,max(1,int(fX)+1))
            interpStructPtr%indY(4,ix,iy) = min(nyFrom,max(1,int(fY)+1))
            interpStructPtr%weight(4,ix,iy) = zx * zy

            interpStructPtr%ixStartFrom = min(interpStructPtr%ixStartFrom, interpStructPtr%indX(1,ix,iy))
            interpStructPtr%iyStartFrom = min(interpStructPtr%iyStartFrom, interpStructPtr%indY(1,ix,iy))
            interpStructPtr%ixEndFrom = max(interpStructPtr%ixEndFrom, interpStructPtr%indX(4,ix,iy))
            interpStructPtr%iyEndFrom = max(interpStructPtr%iyEndFrom, interpStructPtr%indY(4,ix,iy))

            ! Check for global grid running over 180                          
            if(fu_ifLonGlobal(interpStructPtr%gridFrom))then
              if(fX < 1)then
                  interpStructPtr%indX(1,ix,iy) = nxFrom
                  interpStructPtr%indX(3,ix,iy) = nxFrom
                  interpStructPtr%ixEndFrom = nxFrom
              elseif(fX > nxFrom)then
                  interpStructPtr%indX(2,ix,iy) = 1
                  interpStructPtr%indX(4,ix,iy) = 1
                  interpStructPtr%ixStartFrom = 1
              endif
              fX = 1. ! Pretend that the point is inside grid on X-direction for the
              ! out-of-grid  check below
            endif
            if (myiOutside == setMissVal) then
                if (fX>=1 .and. fX<=nxFrom .and.  fY>=1 .and. fY<=nYfrom) then
                  interpStructPtr%ifValid(ix,iy) = .True.
                else
                  interpStructPtr%ifValid(ix,iy) = .false.
                  interpStructPtr%weight(1:4,ix,iy) = real_missing
                endif
            elseif (myiOutside == setZero) then
                if ( (fX<1 .or. fX>nxFrom+1 .or.  fY<1 .or. fY>nYfrom+1)) then
                    interpStructPtr%weight(1:4,ix,iy) = 0. 
                endif
            endif
          end do
        end do

      case(cubic)
        ! Did not bother to make ifValid mask working (Roux)
        call msg('second_order and cubic interpolation methods were cinically killed')
        call set_error('Unimplemeted interpolation method','fu_horiz_interp_struct')
      case default
        call msg('Unknown interpolation method:',interpStructPtr%interp_type)
        call set_error('Unknown interpolation method','fu_horiz_interp_struct')
        return
    end select
    !
    ! Finally, check the values for the envelopes of meaningful areas
    !
    interpStructPtr%ixStartFrom = max(interpStructPtr%ixStartFrom, 1)
    interpStructPtr%iyStartFrom = max(interpStructPtr%iyStartFrom, 1)
    interpStructPtr%ixEndFrom = min(interpStructPtr%ixEndFrom, nxFrom)
    interpStructPtr%iyEndFrom = min(interpStructPtr%iyEndFrom, nyFrom)
    interpStructPtr%ixStartTo = max(interpStructPtr%ixStartTo, 1)
    interpStructPtr%iyStartTo = max(interpStructPtr%iyStartTo, 1)
    interpStructPtr%ixEndTo = min(interpStructPtr%ixEndTo, nxTo)
    interpStructPtr%iyEndTo = min(interpStructPtr%iyEndTo, nyTo)


    !-----------------------------------------------------------------------
    !
    ! Vector rotations, independent of the interpolation itself. 
    !
    if(present(ifMakeRotation))then
      if(.not. ifMakeRotation)return
    endif

    interpStructPtr%ifRotation = .false.   ! same poles, for instance

    select case (interpStructPtr%gridFrom%gridtype)
    case(lonlat) ! interpStructPtr%gridFrom%gridtype
      select case (interpStructPtr%gridTo%gridtype)
      case(lonlat) ! interpStructPtr%gridTo%gridtype
        !
        ! The method: 
        !
        ! 1. Calculate the 3D rotation matrix for rotating
        ! gridFrom -> gridTo 
        !
        ! 2. For each point in gridTo, compute the 2D basis vectors of the tangent 
        ! plane (in 3d cartesian components) wrt the native pole and the pole of gridFrom.
        !
        ! 3. Use these to define rotation matrix which is applied to the (u,v)-pair 
        ! after interpolation.

        if (.not. (fu_pole(interpStructPtr%gridFrom) == fu_pole(interpStructPtr%gridTo))) then
          call msg('Computing wind rotation lonlat to lonlat')
            
          allocate(interpStructPtr%rotation(2, 2, nxTo, nyTo), stat=status)
          if (status /= 0) then
            call set_error('Allocating rotation map failed', 'fu_horiz_interp_struct')
            return
          end if

          interpStructPtr%ifRotation = .true.

          call get_rotation(fu_southpole_lat(fu_pole(interpStructPtr%gridFrom)), &
                          & fu_southpole_lon(fu_pole(interpStructPtr%gridFrom)), &
                          & rotation_3d) ! gridFrom -> geo

          call get_rotation(fu_southpole_lat(fu_pole(interpStructPtr%gridTo)), &
                          & fu_southpole_lon(fu_pole(interpStructPtr%gridTo)), &
                          & rotation_3d_2) ! gridTo -> geo

          rotation_3d = matmul(transpose(rotation_3d_2), rotation_3d) ! gridFrom -> gridTo
          
!!$          call modify_lonlat(fu_southpole_lat(fu_pole(interpStructPtr%gridFrom)),&
!!$                          & fu_southpole_lon(fu_pole(interpStructPtr%gridFrom)),&
!!$                          & pole_geographical, fu_pole(interpStructPtr%gridTo) , &
!!$                          & lat, lon)
!!$
!!$          call get_rotation(lat, lon, rotation_3d) ! gridFrom -> gridTo

          do iy = 1, nyTo
            do ix = 1, nxTo
               
              if (myiOutside == setMissVal) then
                if (.not. interpStructPtr%ifValid(ix,iy)) then
                  interpStructPtr%rotation(:,:, ix,iy) = 0.
                  cycle
                endif
              endif
              !
              ! Here:
              ! lat, lon   : the modified coordinates of a point 
              ! lat2, lon2 : the same point in the reference system (gridTo) 
              ! B1         : the basis vectors at lat, lon in the modified system
              ! B2         : the basis vectors at lat, lon in the reference frame
              ! 
              lat2 = fu_lat_native_from_grid(real(ix), real(iy), interpStructPtr%gridTo)
              lon2 = fu_lon_native_from_grid(real(ix), real(iy), interpStructPtr%gridTo)

              call modify_lonlat(lat2, lon2, fu_pole(interpStructPtr%gridTo), &
                               & fu_pole(interpStructPtr%gridFrom), lat, lon)

              call get_basis(lat, lon, B1)
              call get_basis(lat2, lon2, B2)

              B1 = matmul(rotation_3d, B1) ! rotate B1 to the reference system
              interpStructPtr%rotation(:,:, ix,iy) = matmul(transpose(B2), B1)
              if(abs(interpStructPtr%rotation(1,1, ix,iy)*interpStructPtr%rotation(2,2, ix,iy)- &
               & interpStructPtr%rotation(1,2, ix,iy)*interpStructPtr%rotation(2,1, ix,iy) - 1.) > 0.05)then
                call set_error('rotation changing windspeed','fu_horiz_interp_struct')
                call unset_error('fu_horiz_interp_struct')
              endif
            end do
          end do
        end if

      case(anygrid) !interpStructPtr%gridTo%gridtype

        interpStructPtr%ifRotation = .true. 
        if( allocated(anyGridParamTo%cos_map_rot) .and. &
                    & allocated(anyGridParamTo%sin_map_rot))then

          call msg('Computing wind rotation lonlat to anygrid')
          anyGridParamTo => pAnyGrdParam(interpStructPtr%gridTo%ag%indParam)
          allocate(interpStructPtr%rotation(2, 2, nxTo, nyTo), stat=status)
          if (status /= 0) then
            call set_error('Allocating rotation map failed', 'fu_horiz_interp_struct')
            return
          end if


          call get_rotation(fu_southpole_lat(fu_pole(interpStructPtr%gridFrom)),&
                          & fu_southpole_lon(fu_pole(interpStructPtr%gridFrom)),&
                          & rotation_3d)

          do iy = 1, nyTo
            do ix = 1, nxTo
                
              if (myiOutside == setMissVal)then
                  if(.not. interpStructPtr%ifValid(ix,iy)) then
                     interpStructPtr%rotation(:,:, ix,iy) = 0.
                     cycle
                  endif
              endif

              lat2 = fu_lat_geographical_from_grid(real(ix), real(iy), interpStructPtr%gridTo)
              lon2 = fu_lon_geographical_from_grid(real(ix), real(iy), interpStructPtr%gridTo)

              call modify_lonlat(lat2, lon2, pole_geographical, &
                               & fu_pole(interpStructPtr%gridFrom), lat, lon)
              call get_basis(lat, lon, B1)
              call get_basis(lat2, lon2, B2)

              B1 = matmul(rotation_3d, B1) ! rotate B1 to the reference system
              interpStructPtr%rotation(:,:, ix,iy) = matmul(transpose(B2), B1)

              ! cos and sin of negative angle, for rotating from geographical to rotated
              rot_tmp(1, 1) =  anyGridParamTo%cos_map_rot(nxTo*(iy-1)+ix)
              rot_tmp(1, 2) =  anyGridParamTo%sin_map_rot(nxTo*(iy-1)+ix)        
              rot_tmp(2, 1) =  - rot_tmp(1, 2)
              rot_tmp(2, 2) =  rot_tmp(1, 1)

              interpStructPtr%rotation(:,:, ix,iy) = matmul(interpStructPtr%rotation(:,:, ix,iy), rot_tmp)
              if(abs(interpStructPtr%rotation(1,1, ix,iy)*interpStructPtr%rotation(2,2, ix,iy)- &
               & interpStructPtr%rotation(1,2, ix,iy)*interpStructPtr%rotation(2,1, ix,iy) - 1.) > 0.05)then
                call set_error('rotation changing windspeed','fu_horiz_interp_struct')
                call unset_error('fu_horiz_interp_struct')
              endif
            enddo
          enddo
        else
                call msg('Skipping interp-struct map rotation for gridTo')
                call report(interpStructPtr%gridTo)
                return
        endif


      case default !interpStructPtr%gridTo%gridtype
        call set_error('strange gridtype for gridTo','fu_horiz_interp_struct')
        return
      end select

    case(anygrid)! interpStructPtr%gridFrom%gridtype
      anyGridParamFrom => pAnyGrdParam(interpStructPtr%gridFrom%ag%indParam)
      select case(interpStructPtr%gridTo%gridtype)
      case(lonlat) 
        if(allocated(anyGridParamFrom%cos_map_rot) .and. &
                & allocated(anyGridParamFrom%sin_map_rot))then

          call msg('Computing wind rotation anygrid to lonlat')

          allocate(interpStructPtr%rotation(2, 2, nxTo, nyTo), stat=status)
          if (status /= 0) then
            call set_error('Allocating rotation map failed', 'fu_horiz_interp_struct')
            return
          end if

          interpStructPtr%ifRotation = .true.

          call get_rotation(fu_southpole_lat(fu_pole(interpStructPtr%gridTo)),&
                          & fu_southpole_lon(fu_pole(interpStructPtr%gridTo)),&
                          & rotation_3d)
          xFrom = real_missing
          yFrom = real_missing
          do iy = 1, nyTo
            do ix = 1, nxTo
                
              if (myiOutside == setMissVal)then
                if(.not. interpStructPtr%ifValid(ix,iy)) then
                  interpStructPtr%rotation(:,:, ix,iy) = 0.
                  cycle
                endif
              endif
              
              call project_point_to_grid_xy(interpStructPtr%gridTo, real(ix), real(iy), &
                                         & interpStructPtr%gridFrom, xFrom, yFrom)


              !Rotation is interpolated linearly, no matter what type of structure...
              rot_tmp(1, 1) = fu_2d_interpolation(anyGridParamFrom%cos_map_rot, &
                                              & xFrom, yFrom, nxFrom, nyFrom, linear, notAllowed) 
              rot_tmp(1, 2) = - fu_2d_interpolation(anyGridParamFrom%sin_map_rot, &
                                              & xFrom, yFrom, nxFrom, nyFrom, linear, notAllowed) 
              rot_tmp(2, 1) = - rot_tmp(1, 2)
              rot_tmp(2, 2) = rot_tmp(1, 1)

              lat2 = fu_lat_native_from_grid(real(ix), real(iy), interpStructPtr%gridTo)
              lon2 = fu_lon_native_from_grid(real(ix), real(iy), interpStructPtr%gridTo)

              call modify_lonlat(lat2, lon2, fu_pole(interpStructPtr%gridTo), &
                               & pole_geographical, lat, lon)
              call get_basis(lat, lon, B1)
              call get_basis(lat2, lon2, B2)

              B1 = matmul(transpose(rotation_3d), B1) ! rotate B1 to the reference system
              interpStructPtr%rotation(:,:, ix,iy) = matmul(transpose(B2), B1)
              interpStructPtr%rotation(:,:, ix,iy) = matmul(rot_tmp, interpStructPtr%rotation(:,:, ix,iy))
              if(abs(interpStructPtr%rotation(1,1, ix,iy)*interpStructPtr%rotation(2,2, ix,iy)- &
               & interpStructPtr%rotation(1,2, ix,iy)*interpStructPtr%rotation(2,1, ix,iy) - 1.) > 0.05)then
                call set_error('rotation changing windspeed','fu_horiz_interp_struct')
                call unset_error('fu_horiz_interp_struct')
              endif
            enddo
          enddo
        else
          call msg('Skipping wind rotation for HinterpStruct gridTo:')
          call report(interpStructPtr%gridTo)
        endif
      
      case(anygrid)
        anyGridParamTo   => pAnyGrdParam(interpStructPtr%gridTo%ag%indParam)
        interpStructPtr%ifRotation = .true.

        ! Assumes that gridTo is covered by gridFrom
        if(.not. fu_grids_correspond(interpStructPtr%gridFrom, interpStructPtr%gridTo))then  
          call msg('Computing wind rotation anygrid to anygrid')
          if(      allocated(anyGridParamFrom%cos_map_rot) .and. &
                  & allocated(anyGridParamFrom%sin_map_rot) .and. &
                  & allocated(anyGridParamTo%cos_map_rot) .and. &
                  & allocated(anyGridParamTo%sin_map_rot))then
        
            allocate(interpStructPtr%rotation(2, 2, nxTo, nyTo), stat=status)
            if (status /= 0) then
              call set_error('Allocating rotation map failed', 'fu_horiz_interp_struct')
              return
            end if

            xFrom = real_missing
            yFrom = real_missing
            do iy = 1, nyTo
              do ix = 1, nxTo
                  
                if (myiOutside == setMissVal)then
                  if(.not. interpStructPtr%ifValid(ix,iy)) then
                    interpStructPtr%rotation(:,:, ix,iy) = 0.
                    cycle
                  endif
                endif
                
                call project_point_to_grid(interpStructPtr%gridTo, real(ix), real(iy), &
                                         & interpStructPtr%gridFrom, xFrom, yFrom)


                interpStructPtr%rotation(1, 1, ix,iy) = fu_2d_interpolation(anyGridParamFrom%cos_map_rot, &
                                              & xFrom, yFrom, nxFrom, nyFrom, interpStructPtr%interp_type, notAllowed) 
                interpStructPtr%rotation(1, 2, ix,iy) = - fu_2d_interpolation(anyGridParamFrom%sin_map_rot, &
                                              & xFrom, yFrom, nxFrom, nyFrom, interpStructPtr%interp_type, notAllowed) 
                interpStructPtr%rotation(2, 1, ix,iy) = - interpStructPtr%rotation(1, 2, ix,iy)
                interpStructPtr%rotation(2, 2, ix,iy) = interpStructPtr%rotation(1, 1, ix,iy)

                ! cos and sin of negative angle, for rotating from geographical to rotated
                rot_tmp(1, 1) =  anyGridParamTo%cos_map_rot(nxTo*(iy-1)+ix)
                rot_tmp(1, 2) =  anyGridParamTo%sin_map_rot(nxTo*(iy-1)+ix)        
                rot_tmp(2, 1) =  - rot_tmp(1, 2)
                rot_tmp(2, 2) =  rot_tmp(1, 1)

                interpStructPtr%rotation(:,:, ix,iy) = matmul(interpStructPtr%rotation(:,:, ix,iy), rot_tmp)
                if(abs(interpStructPtr%rotation(1,1, ix,iy)*interpStructPtr%rotation(2,2, ix,iy)- &
                 & interpStructPtr%rotation(1,2, ix,iy)*interpStructPtr%rotation(2,1, ix,iy) - 1.) > 0.05)then
                  call msg("Norm", interpStructPtr%rotation(1,1, ix,iy)*interpStructPtr%rotation(2,2, ix,iy)- &
                                    & interpStructPtr%rotation(1,2, ix,iy)*interpStructPtr%rotation(2,1, ix,iy))
                  call set_error('rotation changing windspeed','fu_horiz_interp_struct')
                  call unset_error('fu_fu_horiz_interp_struct')
                endif
              end do
            end do
          else
            call msg('cannot turn winds for grids')
            call msg('From:')
            call report(interpStructPtr%gridFrom)
            call msg('To:')
            call report(interpStructPtr%gridTo)
          endif
        end if
     
      case default
        call set_error('strange gridtype for gridTo','fu_horiz_interp_struct')
        return
      end select
    case default
      call set_error('strange gridtype for gridFrom','fu_horiz_interp_struct')
      return
    end select

#ifdef DEBUG
   iCount = HorizInterpPool%nHIS
   call msg('Done horiz structure #, ref count', iCount, HorizInterpPool%refCount(iCount))
   call msg('Interp gridFrom: '//trim(gridFrom%name)//" gridTo: "//trim(gridTo%name))
   call msg("iOutside = "+fu_name(myiOutside))
   if (associated(interpStructPtr%ifValid)) then 
           call msg("ifValid mask nOfValid", count(interpStructPtr%ifValid))
           do iy = 1, nyTo
                   call msg(fu_str(iy),interpStructPtr%ifValid(:,iy))
           enddo
   endif
           
   call msg("ixToMin, ixtoMax", interpStructPtr%ixStartTo, interpStructPtr%ixEndTo)
   call msg("iytoMin, iYtoMax", interpStructPtr%iyStartTo, interpStructPtr%iyEndTo)
#endif

    CONTAINS
    
    !============================================================================

    subroutine enlarge_int_counter(arC,arC_work, nCountNew, nCountOld, nx, ny)
      implicit none
      integer, intent(in) :: nCountNew, nCountOld, nx, ny
      integer, dimension(:,:,:), pointer :: arC, arCTmp3d
      integer, dimension(:), pointer :: arC_work, arCTmp
      integer :: iTmp
      !
      ! The goal is to make a longer 1D array, type-cast appropriate dimensions
      ! and then copy the old stuff to the new structure. At the end, free old work array
      ! and redirect the pointers to the new place
      !
      arCTmp => fu_work_int_array(nCountNew*nx*ny)    ! new place
      if(error)return
      arCTmp3d(1:nCountNew,1:nx,1:ny) => arCTmp(1:nCountNew*nx*ny) ! new 3D structure, temporary
      do iTmp = 1, nCountOld
        arCTmp3d(iTmp,1:nx,1:ny) = arC(iTmp,1:nx,1:ny)     ! copy the information
      enddo
      arCTmp3d(1+nCountOld:nCountNew,1:nx,1:ny) = 0          ! zero the empty places

      call free_work_array(arC_work)       ! free old space
      arC_work => arCTmp                   ! old pointer looking at the new space
      arC(1:nCountNew, 1:nx, 1:ny) => arC_work(1:nCOuntNew*nx*ny)  ! 3D type cast
    end subroutine enlarge_int_counter
    
    !============================================================================
    
    subroutine enlarge_real_counter(arC,arC_work, nCountNew, nCountOld, nx, ny)
      implicit none
      integer, intent(in) :: nCountNew, nCountOld, nx, ny
      real, dimension(:,:,:), pointer :: arC, arCTmp3d
      real, dimension(:), pointer :: arC_work, arCTmp
      integer :: iTmp
      !
      ! The goal is to make a longer 1D array, type-cast appropriate dimensions
      ! and then copy the old stuff to the new structure. At the end, free old work array
      ! and redirect the pointers to the new place
      !
      arCTmp => fu_work_array(nCountNew*nx*ny)    ! new place
      if(error)return
      arCTmp3d(1:nCountNew,1:nx,1:ny) => arCTmp(1:nCountNew*nx*ny) ! new 3D structure, temporary
      do iTmp = 1, nCountOld
        arCTmp3d(iTmp,1:nx,1:ny) = arC(iTmp,1:nx,1:ny)     ! copy the information
      enddo
      arCTmp3d(1+nCountOld:nCountNew,1:nx,1:ny) = 0          ! zero the empty places

      call free_work_array(arC_work)       ! free old space
      arC_work => arCTmp                   ! old pointer looking at the new space
      arC(1:nCountNew, 1:nx, 1:ny) => arC_work(1:nCOuntNew*nx*ny)  ! 3D type cast
    end subroutine enlarge_real_counter

  end function fu_horiz_interp_struct

    
  !************************************************************************************
  
  subroutine release_horiz_interp_struct(strct_ptr)
    ! 
    ! Release the interpolation structure: if the structure has a positive reference
    ! count, decrement it, if zero, deallocate the structure. The last remaining
    ! strucuture is move to the freed slot, so rest of the subroutines need no
    ! changes.
    implicit none
    type(THorizInterpStruct), pointer :: strct_ptr
    
    integer :: ind_strct, ind_last_def
    
    if (fu_fails(associated(strct_ptr), 'Struct pointer not associated', 'release_horiz_interp_struct')) return

    ! Find the correct structure
    do ind_strct = 1, HorizInterpPool%nHIS
      if (associated(horizInterpPool%pHIS(ind_strct)%ptr, strct_ptr)) exit
    end do
    if (ind_strct > HorizInterpPool%nHIS) then
      call set_error('Interpolation structure not found in pool', 'release_horiz_interp_struct')
      return
    end if

    if (horizInterpPool%refCount(ind_strct) < 1) then
      call set_error('Cannot release structure: reference count already < 1', 'release_horiz_interp_struct')
      return
    end if

    horizInterpPool%refCount(ind_strct) = horizInterpPool%refCount(ind_strct) - 1
    call msg('Release structure #, ref count', ind_strct, HorizInterpPool%refCount(ind_strct))
    if (horizInterpPool%refCount(ind_strct) == 0) then
      ! Deallocate structure
      call deallocate_horiz_interp_struct(strct_ptr)
      ind_last_def = horizInterpPool%nHIS
      horizInterpPool%nHIS = horizInterpPool%nHIS - 1
      if (ind_last_def /= ind_strct) then
        ! if the last defined one was not just killed, move it to the vacated slot.
        horizInterpPool%pHIS(ind_strct)%ptr => horizInterpPool%pHIS(ind_last_def)%ptr 
        nullify(horizInterpPool%pHIS(ind_last_def)%ptr)
      end if
    end if
    
  end subroutine release_horiz_interp_struct

  !********************************************************************************

  integer function fu_grid_index(nxFrom, ixTo, iyTo, &
                                                & horizInterpStruct) result(iIndex)
    !
    ! The very basic function: returns an index in the From grid, which corresponds
    ! to the given ixTo, iyTo. It can later be used for direct addressing the values from the 
    ! gridFrom. However, it is faster then fu_get_value because that function
    ! makes horizontal interpolation rather than picks a single value
    !
    implicit none

    ! Imported parameters
    type(THorizInterpStruct), intent(in) :: horizInterpStruct
    integer, intent(in) :: nxFrom, ixTo, iyTo

    ! Local variables
    real :: fxFrom, fyFrom
    integer :: iCoef
    logical :: ifValid, ifDegradeMethod
    character(len = *), parameter :: sub_name = 'fu_grid_index'

    !
    !
    iIndex = int_missing
    if(horizInterpStruct%ifInterpolate)then
      !
      ! Full-blown horizontal interpolation
      !
      ifValid=.true.
      if (associated(horizInterpStruct%ifValid)) ifValid = horizInterpStruct%ifValid(iXTo, iYTo)
      if (ifValid) then 
          if ( horizInterpStruct%interp_type ==  summation) then
            call set_error("Called for summation structure", sub_name )
          elseif ( horizInterpStruct%interp_type == average) then
              !!for these structure we want a "center" of non-zero weights
              !! I failed to find any call in Silam that comes here
          
             if(horizInterpStruct%ifGridFromGlobal)then
               !
               ! In global grid cannot interpolate the indices across the longitude closure line
               ! Degrade interpolation along x-direction down to nearest-point if the coefficients are
               ! far from each other, i.e. their mean value strongly differs from the first of them.
               !
               ifDegradeMethod = abs ( sum(horizInterpStruct%indX(1:horizInterpStruct%nCoefs,ixTo,iyTo)) - &
                                     & horizInterpStruct%nCoefs * horizInterpStruct%indX(1,ixTo,iyTo)) > &
                               & horizInterpStruct%nCoefs
             else
               ifDegradeMethod = .false.   ! Non-global grid. Proceed!
             endif
#ifdef DEBUG_IS
if(ifDegradeMethod)then
call msg('Edge of global grid. Indices:',horizInterpStruct%indX(1:horizInterpStruct%nCoefs,ixTo,iyTo))
call msg('Criterion = ', abs ( sum(horizInterpStruct%indX(1:horizInterpStruct%nCoefs,ixTo,iyTo)) - &
                                & horizInterpStruct%nCoefs * horizInterpStruct%indX(1,ixTo,iyTo)))
endif
#endif
             !
             ! Now, do the actual interpolation with whatever method
             !
             if(ifDegradeMethod)then
               !
               ! Nearest-point is needed
               !
               fxFrom = horizInterpStruct%indX( &
                            & maxloc(horizInterpStruct%weight(1:horizInterpStruct%nCoefs, ixTo, iyTo), 1), &
                            & ixTo, iyTo)
             else
               !
               ! No need for nearest-point degraded interpolation. Sum the weighted indices 
               !
               if (horizInterpStruct%nCoefs < 5) then ! normal structures
                 fxFrom = sum(horizInterpStruct%weight(:,ixTo,iyTo) * horizInterpStruct%indX(:,ixTo,iyTo))
               else
                 fxFrom = 0.   ! Potentially long list of weights
                 do iCoef = 1, horizInterpStruct%nCoefs
                   if (horizInterpStruct%weight(iCoef,ixTo,iyTo) == 0.) exit
                   fxFrom = fxFrom + horizInterpStruct%weight(iCoef,ixTo,iyTo) * horizInterpStruct%indX(iCoef,ixTo,iyTo)
                 end do
               endif
             endif  ! whether nearest-point or full interpolation
             !
             ! For y-axis, no global closure can occur
             !
             if (horizInterpStruct%nCoefs < 5) then ! normal structures
               fyFrom = sum(horizInterpStruct%weight(:,ixTo,iyTo) * horizInterpStruct%indY(:,ixTo,iyTo))
             else
               fyFrom = 0.
               do iCoef = 1, horizInterpStruct%nCoefs   ! potentially long list of weights
                 if (horizInterpStruct%weight(iCoef,ixTo,iyTo) == 0.) exit
                 fyFrom = fyFrom + horizInterpStruct%weight(iCoef,ixTo,iyTo) * horizInterpStruct%indY(iCoef,ixTo,iyTo)
               end do
             endif
             !
             ! Finally, index
             iIndex = int(fxFrom+0.5)+(int(fyFrom+0.5)-1)*nxFrom
             if (iIndex < 1) iIndex = int_missing
         else
           !Structure is linear/nearest_point, so simple maxloc would do
          iCoef = maxloc(horizInterpStruct%weight(1:horizInterpStruct%nCoefs, ixTo, iyTo), 1)
          iIndex = horizInterpStruct%indX( iCoef, ixTo, iyTo) + &
                & nxFrom * (horizInterpStruct%indY( iCoef, ixTo, iyTo) - 1 )
         endif
      endif  ! ifValid
    else
      !
      ! No interpolation
      !
      iIndex = ixTo+(iyTo-1)*nxFrom
    endif  ! if horizontal interpolation needed

  end function fu_grid_index
  

  !******************************************************************
  !
  ! Encapsulation of the horizontal interpolation structure
  !
  !******************************************************************

  function fu_gridFrom_from_interp_struct(interpStructHoriz)result(gridFrom)
    implicit none
    type(silja_grid) :: gridFrom
    type(THorizInterpStruct), intent(in) :: interpStructHoriz
    gridFrom = interpStructHoriz%gridFrom
  end function fu_gridFrom_from_interp_struct

  function fu_gridTo_from_interp_struct(interpStructHoriz)result(gridTo)
    implicit none
    type(silja_grid) :: gridTo
    type(THorizInterpStruct), intent(in) :: interpStructHoriz
    gridTo = interpStructHoriz%gridTo
  end function fu_gridTo_from_interp_struct

  integer function fu_nCoefs_interp_str_horiz(interpStructHoriz)
    implicit none
    type(THorizInterpStruct), intent(in) :: interpStructHoriz
    fu_nCoefs_interp_str_horiz = interpStructHoriz%nCoefs
  end function fu_nCoefs_interp_str_horiz

  integer function fu_interpType_interp_str_horiz(interpStructHoriz)
    implicit none
    type(THorizInterpStruct), intent(in) :: interpStructHoriz
    fu_interpType_interp_str_horiz = interpStructHoriz%interp_type
  end function fu_interpType_interp_str_horiz

  subroutine get_coefs_interp_str_horiz(interpStructHoriz, coefsInterp)
    implicit none
    type(THorizInterpCells), intent(out) :: coefsInterp ! Just the size of gridTo
    type(THorizInterpStruct), intent(in) :: interpStructHoriz
    coefsInterp%indX => interpStructHoriz%indX
    coefsInterp%indY => interpStructHoriz%indY
    coefsInterp%weight => interpStructHoriz%weight
    coefsInterp%ifValid => interpStructHoriz%ifValid
    coefsInterp%weight_X => interpStructHoriz%weight_X
    coefsInterp%weight_Y => interpStructHoriz%weight_Y
  end subroutine get_coefs_interp_str_horiz

  subroutine get_coefs_interp_str_cell_horiz(interpStructHoriz, ix, iy, coefsInterpOneCell)
    implicit none
    integer, intent(in) :: ix,iy
    type(THorizInterpStruct), intent(in) :: interpStructHoriz
    type(THorizInterpOneCell), intent(out) :: coefsInterpOneCell ! Just the size of gridTo
    coefsInterpOneCell%indX => interpStructHoriz%indX(:,ix,iy)
    coefsInterpOneCell%indY => interpStructHoriz%indY(:,ix,iy)
    coefsInterpOneCell%weight => interpStructHoriz%weight(:,ix,iy)
    if(associated(interpStructHoriz%weight_X))then
      coefsInterpOneCell%weight_X => interpStructHoriz%weight_X(:,ix,iy)
    else
      nullify(coefsInterpOneCell%weight_X)
    endif
    if(associated(interpStructHoriz%weight_Y))then
      coefsInterpOneCell%weight_Y => interpStructHoriz%weight_Y(:,ix,iy)
    else
      nullify(coefsInterpOneCell%weight_Y)
    endif
    if(associated(interpStructHoriz%ifValid))then
        coefsInterpOneCell%ifValid = interpStructHoriz%ifValid(ix,iy)
    else
        coefsInterpOneCell%ifValid = .true.
    endif
  end subroutine get_coefs_interp_str_cell_horiz

  subroutine get_area_limits(interpStructHoriz, ixStartFrom, iyStartFrom, ixEndFrom, iyEndFrom, &
                                              & ixStartTo, iyStartTo, ixEndTo, iyEndTo)
    type(THorizInterpStruct), intent(in) :: interpStructHoriz
    integer, intent(out) :: ixStartFrom, iyStartFrom, ixEndFrom, iyEndFrom, &
                          & ixStartTo, iyStartTo, ixEndTo, iyEndTo
    ixStartFrom = interpStructHoriz%ixStartFrom
    iyStartFrom = interpStructHoriz%iyStartFrom
    ixEndFrom = interpStructHoriz%ixEndFrom
    iyEndFrom = interpStructHoriz%iyEndFrom
    ixStartTo = interpStructHoriz%ixStartTo
    iyStartTo = interpStructHoriz%iyStartTo
    ixEndTo = interpStructHoriz%ixEndTo
    iyEndTo = interpStructHoriz%iyEndTo
  end subroutine get_area_limits

  function fu_rotation_interp_str_horiz(interpStructHoriz)result(rotation)
    implicit none
    real, dimension(:,:,:,:), pointer :: rotation ! Just the size of gridTo
    type(THorizInterpStruct), intent(in) :: interpStructHoriz
    rotation => interpStructHoriz%rotation
  end function fu_rotation_interp_str_horiz

  function fu_rotation_interp_str_cell_horiz(interpStructHoriz, ix, iy)result(rotation)
    implicit none
    real, dimension(:,:), pointer :: rotation ! Just the size of gridTo
    integer, intent(in) :: ix,iy
    type(THorizInterpStruct), intent(in) :: interpStructHoriz
    rotation => interpStructHoriz%rotation(:,:,ix,iy)
  end function fu_rotation_interp_str_cell_horiz

  function fu_rotation_interp_str_comp_horiz(interpStructHoriz,i, j, ix, iy)result(rotation)
    implicit none
    real :: rotation ! Just the size of gridTo
    integer, intent(in) :: i, j, ix,iy
    type(THorizInterpStruct) :: interpStructHoriz
    rotation = interpStructHoriz%rotation(i,j,ix,iy)
  end function fu_rotation_interp_str_comp_horiz

  function fu_ifWindRotationNeeded(interpStructHoriz) result(ifNeeded)
    implicit none
    type(THorizInterpStruct), intent(in) :: interpStructHoriz
    logical :: ifNeeded
    ifNeeded = interpStructHoriz%ifRotation
  end function fu_ifWindRotationNeeded

  integer function fu_n_non_zero_output_cells(interpStructHoriz)result(iCount)
    ! Actually, returns potentially non-zero cells - those, which have non-zero
    ! contribution from input grid
    implicit none
    type(THorizInterpStruct), intent(in) :: interpStructHoriz
    integer :: ix, iy, iCoef, nxTo,nyTo
    call grid_dimensions(interpStructHoriz%gridTo,nxTo,nyTo)
    iCount = 0
    do iy = 1, nyTo
      do ix = 1, nxTo
        do iCoef = 1, interpStructHoriz%nCoefs
          if(interpStructHoriz%weight_Y(iCoef,ix,iy) > 0.0)then
            iCount = iCount + 1
            exit
          endif
        enddo
      enddo
    enddo
  end function fu_n_non_zero_output_cells

  ! ***************************************************************

  ! ***************************************************************

  !
  !  Numerical horizontal derivatives
  !
  ! ***************************************************************

  ! ***************************************************************



  SUBROUTINE ddx_of_field(grid_data, grd, grd_out, ddx)
    !
    ! Calculates d/dx of field (which has grid grd)
    ! in every gridpoint of grid grd_out.
    !
    ! The two grids must either be exactly the same or match one
    ! of the Arakawa grid systems (half a gridpoint difference in either
    ! direction).
    !
    ! Positive derivative to Silja's positive X-direction eastwards.
    ! 
    ! Unit is [unit-of-grid-data]/m.
    !
    ! Handles any +-1/2 cell shift.
    !
    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, DIMENSION(:) :: grid_data ! the data itself
    TYPE(silja_grid), INTENT(in) :: grd ! the grid the data is in
    TYPE(silja_grid), INTENT(in) :: grd_out ! the grid in which' points
    ! we want to have the values of derivatives

    ! Imported parameters with intent(out):
    REAL, DIMENSION(:), INTENT(out) :: ddx

    ! Local declarations:
    REAL :: lat, dist_x
    REAL :: lat1, dist_x1
    REAL :: lat2, dist_x2
    real :: fTmp, shift_lon, shift_lat
    INTEGER :: i, j, ind1, ind2, nx, ny, nxout, nyout, d, ixtr

    ! Check and initialize.
    
    IF (.NOT.fu_grids_match_closely(grd, grd_out)) THEN
      CALL report(grd)
      CALL report(grd_out)
      CALL set_error('cannot differentiate, no match','ddx_of_field')
      RETURN
    END IF

!    IF (grd%gridtype /= lonlat ) THEN
!      CALL set_error('sorry, can only handle lonlat grids','ddx_of_field')
!      RETURN
!    END IF

    IF (SIZE(grid_data) > SIZE(ddx)) THEN
      CALL set_error('ddx vector too small','ddx_of_field')
      RETURN
    END IF

    CALL grid_dimensions(grd, nx, ny)
    call grid_dimensions(grd_out, nxout, nyout)
    IF (error) RETURN
  
    if(nx /= nxout .or. ny /= nyout) then
      call set_error('Grids are of different size.', 'ddx_of_field')
      return
    end if
    
    call grid_shift_indices(grd_out, grd, shift_lon, shift_lat)

    ! First check the possible shift in x direction. This affects the
    ! finite differencing: if grids are shifted, the derivative is
    ! defined naturally in between two points in the "from" grid.
    
    !$OMP PARALLEL DEFAULT(NONE) SHARED(shift_lon, shift_lat, nx, ny, grd, grid_data, ddx) &
    !$OMP & PRIVATE(i, j, lat, dist_x, d)

    if(shift_lon .eps. 0.)then

      ! The grids are identical in x direction. Use central
      ! differences. Derivative is defined in the middle of the three
      ! points.
      
      !$OMP DO
      do j = 1, ny
        lat = grd%lonlat%sw_corner_modlat + (real(j-1) * grd%lonlat%dy_deg)
        dist_x = 2. * fu_dx_deg_to_m(grd%lonlat%dx_deg, lat)
        ! The interior of grid.
        do i = 2, nx-1
          dist_x = 2*fu_dx_cell_m(grd, i, j)
          ddx((j-1)*nx + i) = (grid_data(nx*(j-1) + i + 1) - grid_data(nx*(j-1) + i - 1)) &       
                              / dist_x
        end do
        
        ! Boundaries, use forward or backward differences (hence half of dist_x).
        ddx((j-1)*nx + 1) = 2. * (grid_data((j-1)*nx + 2) - grid_data((j-1)*nx + 1)) &
                            / dist_x
        ddx(j*nx) = 2. * (grid_data(j*nx) - grid_data(j*nx - 1)) / dist_x

      end do
      !$OMP END DO

    elseif(abs(shift_lon) .eps. 0.5) then
      ! Use central differences, the derivative is defined between two
      ! adjacent points in grd. Depending of the sign of shift, the
      ! first or last cell has to be extrapolated.
              
      if (shift_lon .eps. 0.5) then
        d = 1
      else
        d = 0
      end if

      !$OMP DO
      do j = 1, ny
        lat = grd%lonlat%sw_corner_modlat + (real(j-1+d) * grd%lonlat%dy_deg)
        dist_x = fu_dx_deg_to_m(grd%lonlat%dx_deg, lat)
        
        do i = 1, nx - 1
          dist_x = fu_dx_cell_m(grd, i, j)
          ddx((j-1)*nx + i + d) = (grid_data((j-1)*nx + i + 1) - grid_data((j-1)*nx + i)) &
                                  / dist_x
        end do

        ! Check which cell is missing:
        if (d == 1) then
          ddx((j-1)*nx + 1) = ddx((j-1)*nx + 2)
        else
          ddx(j*nx) = ddx(j*nx - 1)
        end if

      end do
      !$OMP END DO

    else
      call set_error('Grids do not match', 'ddx_of_field')
      !return

    endif

    ! Check the shift in y direction. If there is one, the derivatives
    ! are now defined between the gridpoints, and we need to
    ! interpolate. The last value is extrapolated from the two
    ! previous.
    !
    ! To avoid temporary variables, the grid is scanned from north to
    ! south in the case of positive shift, and south to north for the
    ! negative.

    if(shift_lat .eps. 0.)then
      continue
        
    elseif (shift_lat .eps. 0.5)then
      !$OMP DO
      do j = ny, 2, -1
        do i = 1, nx
          ddx((j-1) * nx + i) = 0.5 * (ddx((j-1) * nx + i) + ddx((j-2)*nx + i))
        end do
      end do
      !$OMP END DO
      !$OMP DO
      do i = 1, nx
        ddx(i) = 2. * ddx(nx + i) - ddx(2*nx + i)
      end do
      !$OMP END DO
    elseif (shift_lat .eps. -0.5)then
      !$OMP DO
      do j = 1, ny-1
        do i = 1, nx
          ddx((j-1) * nx + i) = 0.5 * (ddx((j-1) * nx + i) + ddx(j*nx + i))
        end do
      end do
      !$OMP END DO
      !$OMP DO
      do i = 1, nx
        ddx((ny-1)*nx + i) = 2. * ddx((ny-2)*nx + i) - ddx((ny-2)*nx + i)
      end do
      !$OMP END DO
    else
      call set_error('Grids do not match', 'ddx_of_field')
      !return

    endif
    
    !$OMP END PARALLEL
    
  END SUBROUTINE ddx_of_field

  !***************************************************************************************************


  SUBROUTINE ddy_of_field(grid_data, grd, grd_out, ddy)
    !
    ! Calculates d/dx of field (which has grid grd)
    ! in every gridpoint of grid grd_out.
    !
    ! The two grids must either be exactly the same or match one
    ! of the Arakawa grid systems (half a gridpoint difference in either
    ! direction).
    !
    ! Positive derivative to Silja's positive X-direction eastwards.
    ! 
    ! Unit is [unit-of-grid-data]/m.
    !
    ! Handles a +- half-cell shift in either direction.
    !
    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, DIMENSION(:) :: grid_data ! the data itself
    TYPE(silja_grid), INTENT(in) :: grd ! the grid the data is in
    TYPE(silja_grid), INTENT(in) :: grd_out ! the grid in which' points
    ! we want to have the values of derivatives

    ! Imported parameters with intent(out):
    REAL, DIMENSION(:), INTENT(out) :: ddy

    ! Local declarations:
    REAL :: lat, dist_y
    REAL :: lat1, dist_x1
    REAL :: lat2, dist_x2
    real :: fTmp, shift_lon, shift_lat
    INTEGER :: i, j, ind1, ind2, nx, ny, nxout, nyout, d, ixtr
    real :: metric_m, metric_s, metric_n

    ! Check and initialize.
    
    IF (.NOT.fu_grids_match_closely(grd, grd_out)) THEN
      CALL report(grd)
      CALL report(grd_out)
      CALL set_error('cannot differentiate, no match','ddy_of_field')
      RETURN
    END IF

!    IF (grd%gridtype /= lonlat ) THEN
!      CALL set_error('sorry, can only handle lonlat grids','ddy_of_field')
!      RETURN
!    END IF

    IF (SIZE(grid_data) > SIZE(ddy)) THEN
      CALL set_error('ddy vector too small','ddy_of_field')
      RETURN
    END IF

    CALL grid_dimensions(grd, nx, ny)
    call grid_dimensions(grd_out, nxout, nyout)
    IF (error) RETURN
  
    if(nx /= nxout .or. ny /= nyout) then
      call set_error('Grids are of different size.', 'ddy_of_field')
      return
    end if
    
    call grid_shift_indices(grd_out, grd, shift_lon, shift_lat)

    !$OMP PARALLEL DEFAULT(NONE) SHARED(shift_lon, shift_lat, nx, ny, grd, grid_data, ddy) &
    !$OMP & PRIVATE(i, j, dist_y, d, metric_n, metric_s, metric_m)

    if(shift_lat .eps. 0.)then
      ! The grids are identical in x direction. Use central
      ! differences. Derivative is defined in the middle of the three
      ! points.

!      dist_y = 2. * fu_dy_deg_to_m(grd%lonlat%dy_deg)

      !$OMP DO
      do j = 2, ny-1
        ! The interior of grid.
        do i = 1, nx
          metric_n = fu_dx_cell_m(grd, i, j+1)
          metric_s = fu_dx_cell_m(grd, i, j-1)
          metric_m = fu_dx_cell_m(grd, i, j)
          
          dist_y = 2*fu_dy_cell_m(grd, i, j)
          
          ddy((j-1)*nx + i) = (  grid_data(nx*(j) + i)*metric_n &
                             & - grid_data(nx*(j-2)+i)*metric_s) &       
                              / (dist_y*metric_m)
        end do
      end do
      !$OMP END DO

      ! Boundaries, use forward or backward differences (hence half of dist_y).
      !$OMP DO
      do i = 1, nx
        dist_y = fu_dy_cell_m(grd, i, 1)
        metric_n = fu_dx_cell_m(grd, i, 2)
        metric_s = fu_dx_cell_m(grd, i, 1)
        ddy(i) = (grid_data(nx+i)*metric_n - grid_data(i)*metric_s) / (dist_y*metric_s)
        
        dist_y = fu_dy_cell_m(grd, i, ny)
        metric_n = fu_dx_cell_m(grd, i, ny)
        metric_s = fu_dx_cell_m(grd, i, ny-1)
        ddy((ny-1)*nx + i) = (  grid_data((ny-1)*nx + i)*metric_n &
                            & - grid_data((ny-2)*nx + i)*metric_s) &
                            / (dist_y*metric_n)
      end do
      !$OMP END DO

    elseif (shift_lat .eps. 0.5)then

      ! Use central differences, the derivative is defined between two
      ! adjacent points in grd. Depending of the sign of shift, the
      ! first or last cell has to be extrapolated.
      ! 
      ! The column j=1 has to be extrapolated (actually, copied).

      !dist_y = fu_dy_deg_to_m(grd%lonlat%dy_deg)
      !$OMP DO
      do j = 1, ny - 1
        do i = 1, nx
          ! Here we take the meric coefs from grid of the differentiated data.
          !
          metric_n = fu_dx_cell_m(grd, i, j+1)
          metric_m = fu_dx_cell_m(grd, i, j)
          metric_s = fu_dx_cell_m(grd, i, j-1)          
          dist_y = fu_dy_cell_m(grd, i, j)
          ddy(j*nx + i) = (  grid_data(j*nx + i)*metric_n &
                         & - grid_data((j-1)*nx + i)*metric_s) &
                          / (dist_y*metric_m)
        end do
      end do
      !$OMP END DO
      !$OMP DO
      do i = 1, nx
        ddy(i) = ddy(nx + i)
      end do
      !$OMP END DO
    elseif (shift_lat .eps. -0.5)then
      ! As above, but copy the column j=ny.
      
!      dist_y = fu_dy_deg_to_m(grd%lonlat%dy_deg)
      !$OMP DO
      do j = 1, ny - 1
        do i = 1, nx
          dist_y = fu_dy_cell_m(grd, i, j)
          metric_n = fu_dx_cell_m(grd, i, j+1)
          metric_m = fu_dx_cell_m(grd, i, j)
          metric_s = fu_dx_cell_m(grd, i, j-1)          

          ddy((j-1)*nx + i) = (  grid_data(j*nx + i)*metric_n &
                             & - grid_data((j-1)*nx + i)*metric_s) &
                              / (dist_y*metric_m)
        end do
      end do
      !$OMP END DO
      !$OMP DO
      do i = 1, nx
        ddy((ny-1)*nx + i) = ddy((ny-2)*nx + i)
      end do
      !$OMP END DO

    else
      call set_error('Grids do not match', 'ddy_of_field')
      !return

    endif ! shift_lat

    ! Check the shift in x direction. If there is one, the derivatives
    ! are now defined between the gridpoints, and we need to
    ! interpolate. The last value is extrapolated from the two
    ! previous.
    !
    ! To avoid temporary variables, the grid is scanned from west to
    ! east in the case of negative shift, and east to west for the
    ! positive.

    if(shift_lon .eps. 0.)then
        continue
        
    elseif (shift_lon .eps. 0.5)then
      !$OMP DO
      do j = 1, ny
        do i = nx, 2, -1
          ddy((j-1) * nx + i) = 0.5 * (ddy((j-1)*nx + i) + ddy((j-1)*nx + i - 1))
        end do
      end do
      !$OMP END DO
      !$OMP DO
      do j = 1, ny
        ddy((j-1)*nx + 1) = 2. * ddy((j-1)*nx + 2) - ddy((j-1)*nx + 3)
      end do
      !$OMP END DO

    elseif (shift_lon .eps. -0.5)then
      !$OMP DO
      do j = 1, ny
        do i = 1, nx - 1
          ddy((j-1) * nx + i) = 0.5 * (ddy((j-1) * nx + i) + ddy((j-1)*nx + i + 1))
        end do
      end do
      !$OMP END DO
      !$OMP DO
      do j = 1, ny
        ddy(j*nx) = 2. * ddy(j*nx - 1) - ddy(j*nx - 2)
      end do
      !$OMP END DO

    else
      call set_error('Grids do not match', 'ddy_of_field')
      !return

    endif
    !$OMP END PARALLEL
    
  END SUBROUTINE ddy_of_field

  !************************************************************************************

  subroutine laplacian(grid_data, grid, grid_out, dd)
    implicit none
    
    real, dimension(:), intent(in) :: grid_data
    type(silja_grid) :: grid, grid_out
    real, dimension(:), intent(out) :: dd

    integer :: nx, ny, iy, ix, ii, iw, ie, is, in
    real :: hx,  hy, dd1, dd2, hxn, hxs, hxc

    if (.not. (grid == grid_out)) then
      call set_error('No staggered grids please', 'laplacian')
      return
    end if
    
    ! Laplacian == 1/(hx*hy) ( d/dx(hy/hx df/dx) + d/dy(hx/hy df/dy) )
    ! where hx, hy are the map factors - note that fu_dx_m gives
    ! essentially dx*hx, similar for y.
    
    call grid_dimensions(grid, nx, ny)

    dd(1:nx*ny) = 0.0
    do iy = 2, ny-1
      do ix = 2, nx-1
        ii = (iy-1)*nx + ix
        iw = (iy-1)*nx + ix - 1
        ie = (iy-1)*nx + ix + 1
        is = (iy-2)*nx + ix
        in = iy*nx + ix
        
        hxc = fu_dx_cell_m(grid, ix, iy)
        hy = fu_dy_cell_m(grid, ix, iy)
        ! df/dx in the left interface
        dd1 = (grid_data(ii) - grid_data(iw)) / hxc
        dd2 = (grid_data(ie) - grid_data(ii)) / hxc
        
        dd(ii) = (dd2 - dd1) / hxc

        ! width of the south/north interfaces
        hxn = 0.5 * (hxc + fu_dx_cell_m(grid, ix, iy+1))
        hxs = 0.5 * (hxc + fu_dx_cell_m(grid, ix, iy-1))
        
        ! derivatives at the s/n interfaces
        dd1 = (grid_data(ii) - grid_data(is)) / hy
        dd2 = (grid_data(in) - grid_data(ii)) / hy
        
        dd(ii) = dd(ii) + (dd2*hxn - dd1*hxs) / (hy*hxc)
        !dd(ii) = dd(ii) / (hy*hxc)
      end do
    end do
    
      

  end subroutine laplacian
  


  ! ***************************************************************

  ! ***************************************************************

  !
  !  Grid comparison
  !
  ! ***************************************************************

  ! ***************************************************************



  ! ***************************************************************


  LOGICAL FUNCTION fu_compare_grids_eq (grid1, grid2) result(eq)
    !
    ! Description:
    ! Compares two grids and return true value if they are the same.
    ! True is also returned if both are undefined.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    !
    ! Imported parameters with intent(in):
    TYPE(silja_grid), INTENT(in) :: grid1, grid2

    integer :: i


    IF (grid1%gridtype /= grid2%gridtype) THEN   
      eq = .false.
      RETURN
    END IF

    ! Same type:

    SELECT CASE (grid1%gridtype)

      CASE (lonlat)
        IF( (grid1%lonlat%sw_corner_modlat .eps. grid2%lonlat%sw_corner_modlat ) .and. &
          & (grid1%lonlat%sw_corner_modlon .eps. grid2%lonlat%sw_corner_modlon ) .and. &
          & (grid1%lonlat%nx == grid2%lonlat%nx ) .and. &
          & (grid1%lonlat%ny == grid2%lonlat%ny ) .and. &
          & (grid1%lonlat%pole == grid2%lonlat%pole ) .and. &
          & (grid1%lonlat%dx_deg .eps. grid2%lonlat%dx_deg ) .and. &
          & (grid1%lonlat%dy_deg .eps. grid2%lonlat%dy_deg )) THEN
          eq = .true.
        ELSE
          eq = .false.
        END IF


      case(anygrid)
        !! Same anygrid equals itself!
        eq = .true.
        if (grid1%ag%indParam == grid2%ag%indParam) return

        if((grid1%ag%nx /= grid2%ag%nx ) .or. (grid1%ag%ny /= grid2%ag%ny))then
          eq = .false.
          return
        endif

           

        do i = 1, grid1%ag%nx*grid1%ag%ny
          if(.not.( (pAnyGrdParam(grid1%ag%indParam)%xc(i) &
                   & .eps. pAnyGrdParam(grid2%ag%indParam)%xc(i)) .and. &
                  & (pAnyGrdParam(grid1%ag%indParam)%yc(i) &
                   & .eps. pAnyGrdParam(grid2%ag%indParam)%yc(i)) ))then
            eq = .false.
            exit
          endif
        enddo

      CASE (int_missing)
      eq = (grid2%gridtype == int_missing) 

    CASE default
      eq = .true.
      CALL msg_warning('comparing two unknown grids')

    END SELECT

  END FUNCTION fu_compare_grids_eq

  ! ***************************************************************


  subroutine SubArea_Chk(grid_large, grid_small, ix_shift, iy_shift, ifSubArea)

    !
    ! Checks if grid2 is a subarea of grid1 and if so, returns ix_shift and iy_shift as
    ! grid-coordinates of the first gridcell of grid1 in grid2 
    !

    implicit none

    type(silja_grid), intent(in) :: grid_large, grid_small
    integer, intent(out) :: ix_shift, iy_shift
    logical, intent(out) :: ifSubArea

    integer :: ix, iy, nx_large, ny_large, nx_small, ny_small
    real :: fx_shift, fy_shift

    ifSubArea = .false.
    ix_shift = int_missing
    iy_shift = int_missing

    call grid_dimensions(grid_large, nx_large, ny_large)
    call grid_dimensions(grid_small, nx_small, ny_small)

    if((nx_large < nx_small) .or. (ny_large < ny_small))return

    if(grid_large%gridtype == grid_small%gridtype)then
 
      select case(grid_large%gridtype)
      case(lonlat)
        if(grid_large%lonlat%dx_deg /= grid_small%lonlat%dx_deg)return
        if(grid_large%lonlat%dy_deg /= grid_small%lonlat%dy_deg)return
        if(.not. grid_large%lonlat%pole == grid_small%lonlat%pole)return

        fx_shift = (grid_small%lonlat%sw_corner_modlon - grid_large%lonlat%sw_corner_modlon) / &
                 & grid_large%lonlat%dx_deg
        fy_shift = (grid_small%lonlat%sw_corner_modlat - grid_large%lonlat%sw_corner_modlat) / &
                 & grid_large%lonlat%dy_deg
        ix_shift = nint(fx_shift)
        iy_shift = nint(fy_shift)

        if((ix_shift+nx_small > nx_large) .or. (iy_shift+ny_small > ny_large))return
        if((ix_shift < 0) .or. (iy_shift<0))return
        if(.not.(fx_shift  .eps. real(ix_shift))) return
        if(.not.(fy_shift  .eps. real(iy_shift))) return

        ifSubArea = .true.

      case(anygrid)

        ix_shift = -1
        iy_shift = -1
        do ix = 1, nx_large
          do iy = 1, ny_large
            if((pAnyGrdParam(grid_large%ag%indParam)%xc(nx_large*(iy-1)+ix) .eps. &
              & pAnyGrdParam(grid_small%ag%indParam)%xc(1)) .and. &
            &  (pAnyGrdParam(grid_large%ag%indParam)%yc(nx_large*(iy-1)+ix) .eps. &
              & pAnyGrdParam(grid_small%ag%indParam)%yc(1)))then
              ix_shift = ix-1
              iy_shift = iy-1
              exit
            endif
          enddo
          if(ix_shift>=0)exit
        enddo
        if(ix_shift>=0 .and. (ix_shift+nx_small)<=nx_large .and. (iy_shift+ny_small)<=ny_large)then
          ifSubArea = .true.
          do ix = 1, nx_small
            do iy = 1, ny_small
              if(.not. ((pAnyGrdParam(grid_small%ag%indParam)%xc(nx_small*(iy-1)+ix) .eps. &
                       & pAnyGrdParam(grid_large%ag%indParam)%xc(nx_large*(iy+iy_shift-1)+ix+ix_shift)) &
                       & .and. &
                      & (pAnyGrdParam(grid_small%ag%indParam)%yc(nx_small*(iy-1)+ix) .eps. &
                       & pAnyGrdParam(grid_large%ag%indParam)%yc(nx_large*(iy+iy_shift-1)+ix+ix_shift))))then
                ifSubArea = .false.
                return
              endif
            enddo
          enddo
        endif
      
      case default
        call set_error('input grid of strange type','SubArea_Chk')
      end select

    endif


  end subroutine SubArea_Chk 


  ! ***************************************************************


  logical function fu_grids_correspond(grd1, grd2)
    !
    ! Checks if the grids correspond to each other, meaning that they should
    ! have the same type, projection, pole, cell size etc. The only allowed
    ! difference is - the number of cells in x and y directions. In particular, 
    ! there must be no half-gridcell shift.
    !
    implicit none

    ! Imported parameters with intent IN
    type(silja_grid), intent(in) :: grd1, grd2

    ! Local variables
    real :: vTmp
    integer :: ix, iy, ix_shift, iy_shift, nx1, nx2, ny1, ny2

    fu_grids_correspond = .false.

    IF (grd1%gridtype /= grd2%gridtype) RETURN

!    if (grd1 == grd2) then ! Grid corresponds to itself
!           fu_grids_correspond = .true.
!           return
!    endif
!
    SELECT CASE (grd1%gridtype)

      CASE (lonlat)

        if (.not. grd1%lonlat%pole == grd2%lonlat%pole) return
        IF (.not. (grd1%lonlat%dx_deg .eps. grd2%lonlat%dx_deg)) RETURN
        IF (.not. (grd1%lonlat%dy_deg .eps. grd2%lonlat%dy_deg)) RETURN

        vTmp = modulo(ABS(grd1%lonlat%sw_corner_modlon - &
                         & grd2%lonlat%sw_corner_modlon), grd1%lonlat%dx_deg)
        if(.not.(vTmp.eps.0.))then
          if(.not.(vTmp.eps.grd1%lonlat%dx_deg))return
        end if

        vTmp = modulo(ABS(grd1%lonlat%sw_corner_modlat - &
                         & grd2%lonlat%sw_corner_modlat),grd1%lonlat%dy_deg)
        if(.not.(vTmp.eps.0.))then
          if(.not.(vTmp.eps.grd1%lonlat%dy_deg))return
        end if
       
        fu_grids_correspond = .true.

      case(anygrid)
        !
        ! Assumes that grd2 is covered by grd1
        !
        call grid_dimensions(grd1, nx1, ny1)
        call grid_dimensions(grd2, nx2, ny2)
        ix_shift = -1
        iy_shift = -1
        ! Find offset
        do ix = 1, nx1
          do iy = 1, ny1
            if((pAnyGrdParam(grd1%ag%indParam)%xc(nx1*(iy-1)+ix) .eps. &
              & pAnyGrdParam(grd2%ag%indParam)%xc(1)) .and. &
            &  (pAnyGrdParam(grd1%ag%indParam)%yc(nx1*(iy-1)+ix) .eps. &
              & pAnyGrdParam(grd2%ag%indParam)%yc(1)))then
              ix_shift = ix-1
              iy_shift = iy-1
              exit
            endif
          enddo
          if(ix_shift>=0)exit
        enddo
        ! Check that other cells are the same
        if(ix_shift>=0 .and. iy_shift>=0 .and. (ix_shift+nx2)<=nx1 .and. (iy_shift+ny2)<=ny1)then
          fu_grids_correspond = .true.
          do ix = 1, nx2
            do iy = 1, ny2
              if(.not. ((pAnyGrdParam(grd2%ag%indParam)%xc(nx2*(iy-1)+ix) .eps. &
                       & pAnyGrdParam(grd1%ag%indParam)%xc(nx1*(iy+iy_shift-1)+ix+ix_shift)) .and. &
                      & (pAnyGrdParam(grd2%ag%indParam)%yc(nx2*(iy-1)+ix) .eps. &
                       & pAnyGrdParam(grd1%ag%indParam)%yc(nx1*(iy+iy_shift-1)+ix+ix_shift))))then
                fu_grids_correspond = .false.
                return
              endif
            enddo
          enddo
        endif


    CASE default
      CALL set_error('input grid of strange type','fu_grids_correspond')
      CALL report(grd1)

    END SELECT

  end function fu_grids_correspond


  ! ***************************************************************


  logical function fu_grids_arakawa_correspond(grd1, grd2)   
    !
    ! Checks if the grids correspond to each other up to an Arakawa shift, meaning 
    ! that they must be of the same type, have the same pole and cell size.
    ! The allowed differences are - up to a half-gridcell shift (Arakawa grids)
    ! and different number of cells in x and y directions.
    !
    implicit none

    ! Imported parameters with intent IN
    type(silja_grid), intent(in) :: grd1, grd2

    ! Local variables
    real :: vTmp

    fu_grids_arakawa_correspond = .false.

    IF (grd1%gridtype /= grd2%gridtype) RETURN

    SELECT CASE (grd1%gridtype)

      CASE (lonlat)

        if (.not. (grd1%lonlat%pole == grd2%lonlat%pole)) return
        IF (.not. (grd1%lonlat%dx_deg .eps. grd2%lonlat%dx_deg)) RETURN
        IF (.not. (grd1%lonlat%dy_deg .eps. grd2%lonlat%dy_deg)) RETURN

        vTmp = modulo(ABS(grd1%lonlat%sw_corner_modlon - grd2%lonlat%sw_corner_modlon), grd1%lonlat%dx_deg)
        if(.not.(vTmp.eps.0.))then
          if(.not.(vTmp.eps.grd1%lonlat%dx_deg))then
            if(.not.(vTmp.eps.(grd1%lonlat%dx_deg*0.5)))return
          end if
        end if

        vTmp = modulo(ABS(grd1%lonlat%sw_corner_modlat - grd2%lonlat%sw_corner_modlat), grd1%lonlat%dy_deg)
        if(.not.(vTmp.eps.0.))then
          if(.not.(vTmp.eps.grd1%lonlat%dy_deg))then
            if(.not.(vTmp.eps.(grd1%lonlat%dy_deg*0.5)))return
          end if
        end if

        fu_grids_arakawa_correspond = .true.
    
    case(anygrid)
      !
      ! Allow +- 1 gridcells in each direction (can't allow more, checking would become complicated)
      !
      IF (ABS(grd1%ag%nx - grd2%ag%nx) > 1) RETURN 
      IF (ABS(grd1%ag%ny - grd2%ag%ny) > 1) RETURN
      !
      ! Check the 4 corners, if distance less than dx or dy
      !
      if(abs(pAnyGrdParam(grd1%ag%indParam)%xc(1) - &
           & pAnyGrdParam(grd2%ag%indParam)%xc(1)) > &
       & fu_dx_m_to_deg(pAnyGrdParam(grd1%ag%indParam)%dx(1), &
                      & pAnyGrdParam(grd1%ag%indParam)%yc(1)))return
      if(abs(pAnyGrdParam(grd1%ag%indParam)%yc(1) - &
           & pAnyGrdParam(grd2%ag%indParam)%yc(1)) > &
       & fu_dy_m_to_deg(pAnyGrdParam(grd1%ag%indParam)%dy(1)))return

      if(abs(pAnyGrdParam(grd1%ag%indParam)%xc(grd1%ag%nx) - &
           & pAnyGrdParam(grd2%ag%indParam)%xc(grd2%ag%nx)) > &
       & fu_dx_m_to_deg(pAnyGrdParam(grd1%ag%indParam)%dx(grd1%ag%nx), &
                      & pAnyGrdParam(grd1%ag%indParam)%yc(grd1%ag%nx)))return
      if(abs(pAnyGrdParam(grd1%ag%indParam)%yc(grd1%ag%nx) - &
           & pAnyGrdParam(grd2%ag%indParam)%yc(grd2%ag%nx)) > &
       & fu_dy_m_to_deg(pAnyGrdParam(grd1%ag%indParam)%dy(grd1%ag%nx)))return

      if(abs(pAnyGrdParam(grd1%ag%indParam)%xc(grd1%ag%nx*grd1%ag%ny) - &
           & pAnyGrdParam(grd2%ag%indParam)%xc(grd2%ag%nx*grd2%ag%ny)) > &
       & fu_dx_m_to_deg(pAnyGrdParam(grd1%ag%indParam)%dx(grd1%ag%nx*grd1%ag%ny), &
                      & pAnyGrdParam(grd1%ag%indParam)%yc(grd1%ag%nx*grd1%ag%ny)))return
      if(abs(pAnyGrdParam(grd1%ag%indParam)%yc(grd1%ag%nx*grd1%ag%ny) - &
           & pAnyGrdParam(grd2%ag%indParam)%yc(grd2%ag%nx*grd2%ag%ny)) > &
       & fu_dy_m_to_deg(pAnyGrdParam(grd1%ag%indParam)%dy(grd1%ag%nx*grd1%ag%ny)))return

      if(abs(pAnyGrdParam(grd1%ag%indParam)%xc(grd1%ag%nx*(grd1%ag%ny-1)+1) - &
           & pAnyGrdParam(grd2%ag%indParam)%xc(grd2%ag%nx*(grd2%ag%ny-1)+1)) > &
       & fu_dx_m_to_deg(pAnyGrdParam(grd1%ag%indParam)%dx(grd1%ag%nx*(grd1%ag%ny-1)+1), &
                      & pAnyGrdParam(grd1%ag%indParam)%yc(grd1%ag%nx*(grd1%ag%ny-1)+1)))return
      if(abs(pAnyGrdParam(grd1%ag%indParam)%yc(grd1%ag%nx*(grd1%ag%ny-1)+1) - &
           & pAnyGrdParam(grd2%ag%indParam)%yc(grd2%ag%nx*(grd2%ag%ny-1)+1)) > &
       & fu_dy_m_to_deg(pAnyGrdParam(grd1%ag%indParam)%dy(grd1%ag%nx*(grd1%ag%ny-1)+1)))return

      fu_grids_arakawa_correspond = .true.

    CASE default
      CALL set_error('grid of strange type given' ,'fu_grids_arakawa_correspond')
      CALL report(grd1)

    END SELECT

  end function fu_grids_arakawa_correspond

   
    ! ***************************************************************


    LOGICAL FUNCTION fu_grids_match_closely(grid1, grid2) result(eq)

    ! Description:
    ! Compares two grids and return true value if they match one of 
    ! the Arakawa grid systems.
    ! that is: same gridsize, number of gridpoints +-1, and
    ! corners less than a gridsquare away.
    ! 
    ! Used to check HIRLAM scalar-and vectorgrids.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    !
    ! Imported parameters with intent(in):
    TYPE(silja_grid), INTENT(in) :: grid1, grid2

    eq = .false.

    IF (grid1%gridtype /= grid2%gridtype) RETURN

    SELECT CASE (grid1%gridtype)

      CASE (lonlat)

       IF (.not. (grid1%lonlat%dx_deg .eps. grid2%lonlat%dx_deg)) RETURN

       IF (.not. (grid1%lonlat%dy_deg .eps. grid2%lonlat%dy_deg)) RETURN

       if (.not.(grid1%lonlat%pole == grid2%lonlat%pole))return

       IF (ABS(grid1%lonlat%nx - grid2%lonlat%nx) > 1) RETURN

       IF (ABS(grid1%lonlat%ny - grid2%lonlat%ny) > 1) RETURN

       IF (ABS(grid1%lonlat%sw_corner_modlat - grid2%lonlat%sw_corner_modlat) > &
         & grid1%lonlat%dy_deg) RETURN

       IF (ABS(grid1%lonlat%sw_corner_modlon - grid2%lonlat%sw_corner_modlon) > &
         & grid1%lonlat%dx_deg) RETURN

       eq = .true.


    case(anygrid)
      !
      ! Allow +- 1 gridcells in each direction
      !
      IF (ABS(grid1%ag%nx - grid2%ag%nx) > 1) RETURN 
      IF (ABS(grid1%ag%ny - grid2%ag%ny) > 1) RETURN
      !
      ! Check the 4 corners, if distance less than dx or dy
      !
      if(abs(pAnyGrdParam(grid1%ag%indParam)%xc(1) - &
           & pAnyGrdParam(grid2%ag%indParam)%xc(1)) > &
       & fu_dx_m_to_deg(pAnyGrdParam(grid1%ag%indParam)%dx(1), &
                      & pAnyGrdParam(grid1%ag%indParam)%yc(1)))return
      if(abs(pAnyGrdParam(grid1%ag%indParam)%yc(1) - &
           & pAnyGrdParam(grid2%ag%indParam)%yc(1)) > &
       & fu_dy_m_to_deg(pAnyGrdParam(grid1%ag%indParam)%dy(1)))return

      if(abs(pAnyGrdParam(grid1%ag%indParam)%xc(grid1%ag%nx) - &
           & pAnyGrdParam(grid2%ag%indParam)%xc(grid2%ag%nx)) > &
       & fu_dx_m_to_deg(pAnyGrdParam(grid1%ag%indParam)%dx(grid1%ag%nx), &
                      & pAnyGrdParam(grid1%ag%indParam)%yc(grid1%ag%nx)))return
      if(abs(pAnyGrdParam(grid1%ag%indParam)%yc(grid1%ag%nx) - &
           & pAnyGrdParam(grid2%ag%indParam)%yc(grid2%ag%nx)) > &
       & fu_dy_m_to_deg(pAnyGrdParam(grid1%ag%indParam)%dy(grid1%ag%nx)))return

      if(abs(pAnyGrdParam(grid1%ag%indParam)%xc(grid1%ag%nx*grid1%ag%ny) - &
           & pAnyGrdParam(grid2%ag%indParam)%xc(grid2%ag%nx*grid2%ag%ny)) > &
       & fu_dx_m_to_deg(pAnyGrdParam(grid1%ag%indParam)%dx(grid1%ag%nx*grid1%ag%ny), &
                      & pAnyGrdParam(grid1%ag%indParam)%yc(grid1%ag%nx*grid1%ag%ny)))return
      if(abs(pAnyGrdParam(grid1%ag%indParam)%yc(grid1%ag%nx*grid1%ag%ny) - &
           & pAnyGrdParam(grid2%ag%indParam)%yc(grid2%ag%nx*grid2%ag%ny)) > &
       & fu_dy_m_to_deg(pAnyGrdParam(grid1%ag%indParam)%dy(grid1%ag%nx*grid1%ag%ny)))return

      if(abs(pAnyGrdParam(grid1%ag%indParam)%xc(grid1%ag%nx*(grid1%ag%ny-1)+1) - &
           & pAnyGrdParam(grid2%ag%indParam)%xc(grid2%ag%nx*(grid2%ag%ny-1)+1)) > &
       & fu_dx_m_to_deg(pAnyGrdParam(grid1%ag%indParam)%dx(grid1%ag%nx*(grid1%ag%ny-1)+1), &
                      & pAnyGrdParam(grid1%ag%indParam)%yc(grid1%ag%nx*(grid1%ag%ny-1)+1)))return
      if(abs(pAnyGrdParam(grid1%ag%indParam)%yc(grid1%ag%nx*(grid1%ag%ny-1)+1) - &
           & pAnyGrdParam(grid2%ag%indParam)%yc(grid2%ag%nx*(grid2%ag%ny-1)+1)) > &
       & fu_dy_m_to_deg(pAnyGrdParam(grid1%ag%indParam)%dy(grid1%ag%nx*(grid1%ag%ny-1)+1)))return

      eq = .true.


      CASE default
        CALL set_error('grid of strange type given' ,'fu_grids_match_closely')
        CALL report(grid1)

    END SELECT

  END FUNCTION fu_grids_match_closely


  !**************************************************************************
  
  subroutine grid_shift_indices(grid1, grid2, shift_lon, shift_lat)
    implicit none
    
    ! Find the half-cell shift, if any, between two grids. Returned values (module constants) 
    ! are in relation to the first grid. Example: take x1 in grid1, and x2 in grid2. Then
    ! if x2_i == x1_{i+1/2} ==> shift_lon = plus_half.
    ! 
    ! As with the simpler version above, 0 is returned if grids match.

    type(silja_grid), intent(in) :: grid1, grid2
    real, intent(out) :: shift_lon, shift_lat

    real :: diff

    ! Longitude 

    if (grid1%gridtype /= grid2%gridtype)then
      call set_error('Grids of different types','grid_shift_indices')
      return
    endif

    SELECT CASE (grid1%gridtype)
    CASE (lonlat)

      diff = (grid1%lonlat%sw_corner_modlon &
            & - grid2%lonlat%sw_corner_modlon) / grid1%lonlat%dx_deg
    
      if (abs(diff) < 1.0e-4) then
        shift_lon = 0.
      else if (abs(diff - 0.5) < 1.0e-4) then
        shift_lon = -0.5
      else if (abs(diff + 0.5) < 1.0e-4) then
        shift_lon = 0.5
      else
        shift_lon = real_missing
      end if
    
      !  Latitude

      diff = (grid1%lonlat%sw_corner_modlat &
             - grid2%lonlat%sw_corner_modlat) / grid1%lonlat%dy_deg
    
      if (abs(diff) < 1.0e-4) then
        shift_lat = 0.
      else if (abs(diff - 0.5) < 1.0e-4) then
        shift_lat = -0.5
      else if (abs(diff + 0.5) < 1.0e-4) then
        shift_lat = 0.5
      else
        shift_lat = real_missing
      end if
    
    case(anygrid)

      if(fu_grids_arakawa_correspond(grid1, grid2))then
        !
        ! Now that we've checked the 4 corners and believe that grids belong 
        ! to arakawa system, compute the shift only for southwest corner
        !
        diff = (pAnyGrdParam(grid1%ag%indParam)%xc(1) -  &
              & pAnyGrdParam(grid2%ag%indParam)%xc(1)) / &
             & fu_dx_m_to_deg(pAnyGrdParam(grid1%ag%indParam)%dx(1), &
                            & pAnyGrdParam(grid2%ag%indParam)%yc(1))

        if (abs(diff) < 1.0e-4) then
          shift_lon = 0.
        else if (abs(diff - 0.5) < 1.0e-4) then
          shift_lon = -0.5
        else if (abs(diff + 0.5) < 1.0e-4) then
          shift_lon = 0.5
        else
          call msg_warning('anygrids check for arakawa but cannot compute lon shift','grid_shift_indices')
          call msg('diff = ', diff)
          shift_lon = real_missing
        end if

        diff = (pAnyGrdParam(grid1%ag%indParam)%yc(1) -  &
              & pAnyGrdParam(grid2%ag%indParam)%yc(1)) / &
             & fu_dy_m_to_deg(pAnyGrdParam(grid1%ag%indParam)%dy(1))

        if (abs(diff) < 1.0e-4) then
          shift_lat = 0.
        else if (abs(diff - 0.5) < 1.0e-4) then
          shift_lat = -0.5
        else if (abs(diff + 0.5) < 1.0e-4) then
          shift_lat = 0.5
        else
          call msg_warning('anygrids check for arakawa but cannot compute lat shift','grid_shift_indices')
          call msg('diff = ', diff)
          shift_lat = real_missing
        end if

      else
        shift_lon = real_missing
        shift_lat = real_missing
      endif

    case default

      CALL set_error('grid of strange type given' ,'grid_shift_indices')
      return

    end select


  end subroutine grid_shift_indices





  ! ***************************************************************

  ! ***************************************************************
  !
  !
  !       grid transformation
  !
  !
  ! ***************************************************************
  ! ***************************************************************


  ! ***************************************************************


  SUBROUTINE adjust_grid_to_sample(grid, grid_sample, ifAdjusted)
    !
    ! Compares grid with grid_sample and, in case of "closely matched" grids
    ! ensures the same dimensions of grid to be as of the grid_sample
    ! When to call:
    ! If some field belongs to Arakawa grid system, its half-cell shift may be
    ! enough for wrong selection of an area (see data_horizontal_select).
    ! For example, if the system_grid based on temperature or orography may have
    ! different size in x- and/or y- dimensions, which totally destroys all
    ! further computations.
    ! Problem if called for the system_grid adjustment: 
    ! This function does not care about availability of the data for ajusted grid.
    ! So, if changes were made, it may be needed to check if the new grid is
    ! filled with data. If not - the system_grid was incorectly selected and the 
    ! model has to be stopped.
    !
    ! Return value: silja_true if adjustment was done, silja_false is nothing to do
    ! and silja_undefined if the grid can not be adjusted - not close for example
    !
    IMPLICIT NONE

    ! Imported parameters with intent INOUT
    TYPE(silja_grid), INTENT (inout) :: grid
    type(silja_grid), intent(in) :: grid_sample

    ! Imported parameters with intent OUT
    TYPE(silja_logical), INTENT (out) :: ifAdjusted 
                         ! silja_true => adjustment made
                         ! silja_false => nothing to be done
                         ! silja_undefined => grids undefined or not close enough

    IF(.not.defined(grid_sample))THEN
      CALL msg_warning('grid_sample is undefined, no adjustment is done')
      ifAdjusted = silja_false
      RETURN
    END IF

    ifAdjusted = silja_undefined

    !------------------------------------
    !
    ! Check that grids are close enough (similar but not the same to fu_grid_match_closely)

    IF (grid%gridtype /= grid_sample%gridtype) RETURN

    SELECT CASE (grid_sample%gridtype)
    CASE (lonlat)
      IF (.not.(grid%lonlat%dx_deg .eps. grid_sample%lonlat%dx_deg)) RETURN

      IF (.not.(grid%lonlat%dy_deg .eps. grid_sample%lonlat%dy_deg)) RETURN

      IF (ABS(grid%lonlat%nx - grid_sample%lonlat%nx) > 3) RETURN !!!!!!

      IF (ABS(grid%lonlat%ny - grid_sample%lonlat%ny) > 3) RETURN !!!!!!

      IF (ABS(grid%lonlat%sw_corner_modlat - grid_sample%lonlat%sw_corner_modlat) &
           & > grid_sample%lonlat%dy_deg) RETURN

      IF (ABS(grid%lonlat%sw_corner_modlon - grid_sample%lonlat%sw_corner_modlon) &
          & > grid_sample%lonlat%dx_deg) RETURN



      !--------------------------------------
      ! May be, nothing to do ?

      IF(grid%lonlat%nx == grid_sample%lonlat%nx .and. &
       & grid%lonlat%ny == grid_sample%lonlat%ny) THEN

        ifAdjusted = fu_set_false()
        RETURN

      END IF

      !--------------------------------------
      ! Grids are close enough but not the same size => adjust dimensions
      ! It may be needed to adjust the corner, which is neglected so far
     
!     CALL report(grid)

      grid%lonlat%ny = grid_sample%lonlat%ny
      ! grid might stick beyond the north pole.
      ! Here is a dirty workaround for that
      if (grid%lonlat%sw_corner_modlat + (grid%lonlat%ny-1) * grid%lonlat%dy_deg > 90.) then
              call msg_warning("adjusted grid sticks beyond the North pole, removing one cell",&
              &"adjust_grid_to_sample")
              grid%lonlat%ny = grid%lonlat%ny - 1
              IF(grid%lonlat%ny == grid_sample%lonlat%ny .and. grid%lonlat%nx == grid_sample%lonlat%nx) THEN
                      ifAdjusted = silja_false
                      return
              endif
      endif
      grid%lonlat%nx = grid_sample%lonlat%nx


      ifAdjusted = fu_set_true()


    case(anygrid)

      IF (ABS(grid%ag%nx - grid_sample%ag%nx) > 3) RETURN 
      IF (ABS(grid%ag%ny - grid_sample%ag%ny) > 3) RETURN 
      ! can't enlarge
      IF (grid%ag%nx < grid_sample%ag%nx) RETURN 
      IF (grid%ag%ny < grid_sample%ag%ny) RETURN 
      ! check if corners close enough
      if(abs(pAnyGrdParam(grid%ag%indParam)%xc(1) - &
           & pAnyGrdParam(grid_sample%ag%indParam)%xc(1)) > &
       & fu_dx_m_to_deg(pAnyGrdParam(grid%ag%indParam)%dx(1), &
                      & pAnyGrdParam(grid%ag%indParam)%yc(1)))return
      if(abs(pAnyGrdParam(grid%ag%indParam)%yc(1) - &
           & pAnyGrdParam(grid_sample%ag%indParam)%yc(1)) > &
       & fu_dy_m_to_deg(pAnyGrdParam(grid%ag%indParam)%dy(1)))return

      if(abs(pAnyGrdParam(grid%ag%indParam)%xc(grid%ag%nx) - &
           & pAnyGrdParam(grid_sample%ag%indParam)%xc(grid_sample%ag%nx)) > &
       & fu_dx_m_to_deg(pAnyGrdParam(grid%ag%indParam)%dx(grid%ag%nx), &
                      & pAnyGrdParam(grid%ag%indParam)%yc(grid%ag%nx)))return
      if(abs(pAnyGrdParam(grid%ag%indParam)%yc(grid%ag%nx) - &
           & pAnyGrdParam(grid_sample%ag%indParam)%yc(grid_sample%ag%nx)) > &
       & fu_dy_m_to_deg(pAnyGrdParam(grid%ag%indParam)%dy(grid%ag%nx)))return

      if(abs(pAnyGrdParam(grid%ag%indParam)%xc(grid%ag%nx*grid%ag%ny) - &
           & pAnyGrdParam(grid_sample%ag%indParam)%xc(grid_sample%ag%nx*grid_sample%ag%ny)) > &
       & fu_dx_m_to_deg(pAnyGrdParam(grid%ag%indParam)%dx(grid%ag%nx*grid%ag%ny), &
                      & pAnyGrdParam(grid%ag%indParam)%yc(grid%ag%nx*grid%ag%ny)))return
      if(abs(pAnyGrdParam(grid%ag%indParam)%yc(grid%ag%nx*grid%ag%ny) - &
           & pAnyGrdParam(grid_sample%ag%indParam)%yc(grid_sample%ag%nx*grid_sample%ag%ny)) > &
       & fu_dy_m_to_deg(pAnyGrdParam(grid%ag%indParam)%dy(grid%ag%nx*grid%ag%ny)))return

      if(abs(pAnyGrdParam(grid%ag%indParam)%xc(grid%ag%nx*(grid%ag%ny-1)+1) - &
           & pAnyGrdParam(grid_sample%ag%indParam)%xc(grid_sample%ag%nx*(grid_sample%ag%ny-1)+1)) > &
       & fu_dx_m_to_deg(pAnyGrdParam(grid%ag%indParam)%dx(grid%ag%nx*(grid%ag%ny-1)+1), &
                      & pAnyGrdParam(grid%ag%indParam)%yc(grid%ag%nx*(grid%ag%ny-1)+1)))return
      if(abs(pAnyGrdParam(grid%ag%indParam)%yc(grid%ag%nx*(grid%ag%ny-1)+1) - &
           & pAnyGrdParam(grid_sample%ag%indParam)%yc(grid_sample%ag%nx*(grid_sample%ag%ny-1)+1)) > &
       & fu_dy_m_to_deg(pAnyGrdParam(grid%ag%indParam)%dy(grid%ag%nx*(grid%ag%ny-1)+1)))return

      !--------------------------------------
      ! May be, nothing to do ?

      IF(grid%ag%nx == grid_sample%ag%nx .and. &
       & grid%ag%ny == grid_sample%ag%ny) THEN

        ifAdjusted = fu_set_false()
        RETURN

      END IF
 
      !--------------------------------------
      ! Grids are close enough but not the same size => adjust dimensions
      ! Might need to change lons and lats fields, 
 
      grid%lonlat%nx = grid_sample%lonlat%nx
      grid%lonlat%ny = grid_sample%lonlat%ny
      ifAdjusted = fu_set_true()


    CASE default
      CALL set_error('cannot handle sample grid:','fu_adjust_grid_to_sample')
      CALL report(grid_sample)

    END SELECT

!    print *,'==============================================='
!    print *,'                GRID ADJUSTMENT: '

!    CALL report(grid)



  END SUBROUTINE adjust_grid_to_sample

  
  !*********************************************************************************************

  LOGICAL FUNCTION fu_if_mesh_intrpolatable(mesh_sml, mesh_lrg, covermask, ifSpk)
    !
    ! Checks all points of  mesh_sml can be obtained by interpolation
    ! from mesh_lrg
    ! Returns mask of large grid points needed for such interpolation
    ! No extraploation is allowed
    ! Stupid and bruteforcish method -- Not to be frequenlly used
    !
    ! Projects all points of small grid into the large one and fill the mask
    !

    IMPLICIT NONE

    TYPE(silja_grid), INTENT(in) :: mesh_sml, mesh_lrg
    logical, INTENT(in), optional :: ifSpk
    integer, dimension(:,:), optional, INTENT(out) :: covermask

    integer :: ix, iy, ix1, ix2, iy1, iy2, nx_lrg, ny_lrg, nx_sml, ny_sml
    real :: xNew, yNew
    logical :: ifLonGlobalLrg, ifLonGlobalSml,  ifNP, ifSP, ifSpeak, ifPointOK

    ifSpeak = .false.
    if (present( ifSpk)) ifSpeak = ifSpk

    call grid_dimensions(mesh_lrg, nx_lrg, ny_lrg)
    call grid_dimensions(mesh_sml, nx_sml, ny_sml)

    fu_if_mesh_intrpolatable = .true.
    !!!At voima/teho equal grids fail sometimes due to numerics
    !! Dirty hack
    if (mesh_sml == mesh_lrg) then  
       if (present(covermask)) covermask(1:nx_lrg,1:ny_lrg) = 1
       return
    endif



    ifLonGlobalSml = fu_ifLonGlobal(mesh_sml) ! Not to bother about X
    ifLonGlobalLrg = fu_ifLonGlobal(mesh_lrg) ! Not to bother about X
    ifNP = fu_ifPoleIncluded(mesh_lrg, northern_boundary) !Ignore stick out at North
    ifSP = fu_ifPoleIncluded(mesh_lrg,southern_boundary) !Ignore stick out at South


    if (present(covermask)) covermask(1:nx_lrg,1:ny_lrg) = 0

    !! Can be done OMP sharing covermask and fu_if_mesh_intrpolatable
    xNew = real_missing
    yNew = real_missing
    do ix = 1, nx_sml
      do iy = 1, ny_sml
        ifPointOK = .true.
        call project_point_to_grid_xy(mesh_sml, real(ix), real(iy), mesh_lrg, xNew, yNew)

        if( ((xNew < 1 .or. xNew > nx_lrg) .and. (.not. ifLonGlobalLrg)) .or. &
          & (yNew < 1 .and. (.not. ifSP)) .or. (yNew > ny_lrg .and.(.not. ifNP))) then
             !Point from mesh_sml sticks out -- drop the flag
           fu_if_mesh_intrpolatable = .false.
           ifPointOK = .false.
           
           ! if (ifSpeak) call msg("ix, iy, xNew, Ynew OOO", (/real(ix), real(iy), xNew, yNew/))
           if ( .not. present(covermask)) return
!        else 
!           if(ix==1 .and. iy==6) call msg("ix, iy, xNew, Ynew III", (/real(ix), real(iy), xNew, yNew/))

           !if (ifSpeak) call msg("ix, iy, xNew, Ynew III", (/real(ix), real(iy), xNew, yNew/))
        endif

        if ( present(covermask) .and. ifPointOK) then !Handle covermask

          if (ifLonGlobalLrg .and. (xNew <= 1.  .or. xNew > nx_lrg)) then
            ! Should handle the situation when xNew==1
            ix1 = nx_lrg   !Wrap around
            ix2 = 1
          else
            ix1 = floor(xNew)
            ix2 = ceiling(xNew)
          endif
        
          iy1 = floor(yNew)
          iy2 = ceiling(yNew)

          covermask(ix1,iy1:iy2) = 1 ! These cells are needed, Obs: possible ix1>ix2
          covermask(ix2,iy1:iy2) = 1 ! 
        endif

      enddo
    enddo


  end FUNCTION fu_if_mesh_intrpolatable



  ! ***************************************************************


  LOGICAL FUNCTION fu_if_grid_covered(grid_sml, grid_lrg)
    !
    ! Checks that grid_sml is completely covered by grid_lrg
    ! Method: all borders are projected onto grid_lrg.
    !

    IMPLICIT NONE

    TYPE(silja_grid), INTENT(in) :: grid_sml, grid_lrg

    integer :: ix, iy, nx_lrg, ny_lrg, nx_sml, ny_sml
    real :: xNew, yNew
    logical :: ifLonGlobal, ifLatGlobal

    fu_if_grid_covered = .false.
    ifLonGlobal = fu_ifLonGlobal(grid_lrg)
    ifLatGlobal = fu_ifLatGlobal(grid_lrg)

    if(ifLonGlobal .and. ifLatGlobal)then

        fu_if_grid_covered = .true.     ! fully global grid covers everything...

    elseif (grid_lrg == grid_sml) then
       fu_if_grid_covered = .true.
    elseif(grid_lrg%gridtype == lonlat .and. grid_sml%gridtype == lonlat .and. &
         & grid_sml%lonlat%pole == grid_lrg%lonlat%pole)then
        !
        ! Speed-up: no reprojection needed. Note possible lon and lat separate global closures
        !
        fu_if_grid_covered = ifLonGlobal .or. &
           & ((grid_sml%lonlat%sw_corner_modlon  - 0.5 * grid_sml%lonlat%dx_deg >= &
            & grid_lrg%lonlat%sw_corner_modlon - 0.5 * grid_lrg%lonlat%dx_deg - 1.e-6) .and. &
           & (grid_sml%lonlat%sw_corner_modlon + (grid_sml%lonlat%nx-0.5) * grid_sml%lonlat%dx_deg - 2.e-5  <= &
            & grid_lrg%lonlat%sw_corner_modlon + (grid_lrg%lonlat%nx-0.5) * grid_lrg%lonlat%dx_deg))
        if(.not. fu_if_grid_covered)return

        fu_if_grid_covered = ifLatGlobal .or. &
          & ((grid_sml%lonlat%sw_corner_modlat  - 0.5 * grid_sml%lonlat%dy_deg >= &
            & grid_lrg%lonlat%sw_corner_modlat - 0.5 * grid_lrg%lonlat%dy_deg - 1.e-6) .and. &
           & (grid_sml%lonlat%sw_corner_modlat + (grid_sml%lonlat%ny-0.5) * grid_sml%lonlat%dy_deg - 2.e-5  <= &
            & grid_lrg%lonlat%sw_corner_modlat + (grid_lrg%lonlat%ny-0.5) * grid_lrg%lonlat%dy_deg))

    else
        !
        ! Have to check the boundary line with reprojection
        ! 
      call grid_dimensions(grid_lrg, nx_lrg, ny_lrg)
      call grid_dimensions(grid_sml, nx_sml, ny_sml)
      xNew  = real_missing
      yNew  = real_missing

      do ix = 1, nx_sml+1
        call project_point_to_grid_xy(grid_sml, real(ix)-0.5, 0.5, grid_lrg, xNew, yNew)
        if(xNew < 0.5 .or. xNew > nx_lrg+0.5 .or. &
         & yNew < 0.5 .or. yNew > ny_lrg+0.5) return
        call project_point_to_grid_xy(grid_sml, real(ix)-0.5, real(ny_sml)+0.5, &
                                      & grid_lrg, xNew, yNew)
        if(xNew < 0.5 .or. xNew > nx_lrg+0.5 .or. &
         & yNew < 0.5 .or. yNew > ny_lrg+0.5) return
      end do
      do iy = 1, ny_sml+1
        call project_point_to_grid_xy(grid_sml, 0.5, real(iy)-0.5, grid_lrg, xNew, yNew)
        if(xNew < 0.5 .or. xNew > nx_lrg+0.5 .or. &
         & yNew < 0.5 .or. yNew > ny_lrg+0.5) return
        call project_point_to_grid_xy(grid_sml, real(nx_sml)+0.5, real(iy)-0.5, &
                                      & grid_lrg, xNew, yNew)
        if(xNew < 0.5 .or. xNew > nx_lrg+0.5 .or. &
         & yNew < 0.5 .or. yNew > ny_lrg+0.5) return
      end do
      
      
      if(ifLatGlobal)then ! Check for big grids poles inside the small grid
        call project_point_to_grid_lonlat(&
                           &fu_southpole_lon(fu_pole(grid_lrg)), &
                           &fu_southpole_lat(fu_pole(grid_lrg)), grid_sml, xNew, yNew)
        if (xNew >= 1.0 .and. xNew <= nx_sml .and. yNew >= 1.0 .and. yNew <= ny_sml) then
                if (grid_lrg%lonlat%sw_corner_modlat > -89.99999) then
                    return ! large grid does not cover the pole
                endif
        endif

        call project_point_to_grid_lonlat(&
                     & fu_northpole_lon(fu_pole(grid_lrg)), &
                     & fu_northpole_lat(fu_pole(grid_lrg)), grid_sml, xNew, yNew)
        if(xNew >= 1.0 .and. xNew <= nx_sml .and. yNew >= 1.0 .and. yNew <= ny_sml) then
                if (grid_lrg%lonlat%sw_corner_modlat + (grid_lrg%lonlat%ny-0.5) * grid_lrg%lonlat%dy_deg < 89.99999) then
                    return ! large grid does not cover the pole
                endif
        endif
      endif    
        
      fu_if_grid_covered = .true.

    endif

  END FUNCTION fu_if_grid_covered


  !*****************************************************************

  real function fu_area_coverage(grid_to_cover, studied_grid)
    !
    ! This function roughly computes the fraction of the grid_to_cover 
    ! covered by the studied_grid. The algorithm is quite slow in general
    ! case because each element of the grid_to_cover is checked 
    ! individually. So, DO NOT USE this function too often.
    !
    implicit none

    ! Imported parameters with intent IN
    type(silja_grid), intent(in), target :: grid_to_cover, studied_grid

    ! Local variables
    integer :: i,j, nx, ny, nCovered, nxS, nyS
    type(silja_grid), pointer :: GTC, SG ! Just abbreviations of input ones
    real :: xNew, yNew

    if(.not.defined(grid_to_cover))then
      call set_error('Undefined grid_to_cover','fu_area_coverage')
      return
    end if

    if(.not.defined(studied_grid))then
      call set_error('Undefined studied_grid','fu_area_coverage')
      return
    end if

    GTC => grid_to_cover
    SG => studied_grid
    
    ! Global grids cover everything 
    if(fu_ifLonGlobal(SG) .and. fu_ifLatGlobal(SG))then
      fu_area_coverage = 1.
      return
    endif
    
    ! If same rotation, can use speedup
    if(GTC%gridtype == lonlat .and. SG%gridtype == lonlat .and. GTC%lonlat%pole == SG%lonlat%pole)then
      if(fu_ifLatGlobal(SG) .or. &
       & GTC%lonlat%sw_corner_modlat - 0.5 * GTC%lonlat%dy_deg >= &
         & SG%lonlat%sw_corner_modlat - 0.5 * SG%lonlat%dy_deg - 1.e-6 .and. &
       & GTC%lonlat%sw_corner_modlat + GTC%lonlat%dy_deg * (GTC%lonlat%ny-1) <=  &
         & SG%lonlat%sw_corner_modlat + SG%lonlat%dy_deg * (SG%lonlat%ny-1) + 1.e-6)then 
        if(fu_ifLonGlobal(SG) .or. &
         & GTC%lonlat%sw_corner_modlon - 0.5 * GTC%lonlat%dx_deg >= &
           & SG%lonlat%sw_corner_modlon - 0.5 * SG%lonlat%dx_deg - 1.e-6 .and. &
         & GTC%lonlat%sw_corner_modlon + GTC%lonlat%dx_deg * (GTC%lonlat%nx-0.5) <=  &
           & SG%lonlat%sw_corner_modlon + SG%lonlat%dx_deg * (SG%lonlat%nx-0.5) + 1.e-6)then
          fu_area_coverage = 1.
          return

        end if
      end if
    end if

    ! Scan all cells of the studied_grid and check each for placement. 
    call grid_dimensions(GTC, nx,ny)
    call grid_dimensions(SG, nxS,nyS)
    nCovered = 0
    do j =1,ny
      do i=1, nx
        call project_point_to_grid(GTC, real(i), real(j), SG, xNew, yNew)
        if(xNew>=1 .and. yNew>=1 .and. xNew<=nxS .and. yNew<=nyS) NCovered=NCovered+1
      end do
    end do

    fu_area_coverage = real(nCovered) / (real(nx)*real(ny))

!    call msg("")
!    call msg("")
!    call msg("")
!    call msg("fu_area_coverage SG")
!    call report(SG)
!    call report(GTC)
!    call msg("Coverage:", fu_area_coverage)
!
!    call msg("fu_area_coverage")

  end function fu_area_coverage



  !***************************************************************

  logical function fu_point_inside_grid(grid, lon, lat, pole)
    !
    ! Checks if the point is inside the given grid. 
    !
    implicit none

    ! imported parameters
    type(silja_grid), intent(in) :: grid
    real, intent(in) :: lon, lat
    type(silam_pole), intent(in) :: pole
    integer :: nx, ny
    real :: lat_geo, lon_geo, xnew, ynew

    fu_point_inside_grid=.false.
    
    if(.not.defined(grid))then
      call set_error('Undefined grid given','fu_point_inside_grid')
      return
    end if
    if(.not.defined(pole))then
      call set_error('Undefined pole given','fu_point_inside_grid')
      return
    end if    
    
    call modify_lonlat(lat, lon, pole, pole_geographical, lat_geo, lon_geo)
    call grid_dimensions(grid, nx, ny)
    call project_point_to_grid(lat_geo, lon_geo, grid, xNew, yNew)
    if(error)return
    fu_point_inside_grid = (xNew>0.5 .and. yNew>0.5 .and. xNew<nx+0.5 .and. yNew<ny+0.5)

  end function fu_point_inside_grid


  !***************************************************************

  subroutine write_anygrid_2_gradsfile(grid, fNm)
    implicit none

    type(silja_grid), intent(in) :: grid
    integer :: nPoints, nVars, rec
    logical :: ifdx, ifdy, ifsin, ifcos
    integer :: u_grd_ctl, u_grd_bin
    CHARACTER (LEN=fnlen) :: fNm

    if(grid%gridtype /= anygrid)then
      call set_error('Only for grids of type anygrid','write_anygrid_2_gradsfile')
      return
    endif
    
    nPoints = grid%ag%nx*grid%ag%ny
    nVars = 0
    ifdx = .false.
    ifdy = .false.
    ifsin = .false.
    ifcos = .false.

    ! at least lon and lat have to exist; count the others
    if(allocated(pAnyGrdParam(grid%ag%indParam)%xc))then
      if(size(pAnyGrdParam(grid%ag%indParam)%xc) >= nPoints)then
        nVars = nVars + 1
     endif
    endif
    if(allocated(pAnyGrdParam(grid%ag%indParam)%yc))then
      if(size(pAnyGrdParam(grid%ag%indParam)%yc) >= nPoints)then
        nVars = nVars + 1
      endif
    endif
    if(nVars /= 2)then
      call set_error('problem with lats or lons','write_anygrid_2_gradsfile')
      return
    endif
    if(allocated(pAnyGrdParam(grid%ag%indParam)%dx))then
      if(size(pAnyGrdParam(grid%ag%indParam)%dx) >= nPoints)then
        nVars = nVars + 1
        ifdx = .true.
      endif
    endif
    if(allocated(pAnyGrdParam(grid%ag%indParam)%dy))then
      if(size(pAnyGrdParam(grid%ag%indParam)%dy) >= nPoints)then
        nVars = nVars + 1
        ifdy = .true.
      endif
    endif
    if(allocated(pAnyGrdParam(grid%ag%indParam)%sin_map_rot))then
     if(size(pAnyGrdParam(grid%ag%indParam)%sin_map_rot) >= nPoints)then
        nVars = nVars + 1
        ifsin = .true.
      endif
    endif
    if(allocated(pAnyGrdParam(grid%ag%indParam)%cos_map_rot))then
      if(size(pAnyGrdParam(grid%ag%indParam)%cos_map_rot) >= nPoints)then
        ifcos = .true.
        nVars = nVars + 1
      endif
    endif

    u_grd_ctl = fu_next_free_unit()
    open(u_grd_ctl, file=fu_connect_strings(fNm, '.ctl'))
    WRITE(u_grd_ctl,'(A5,A)') 'DSET ', fu_trim_grads_hat(fNm,fNm)
    WRITE(u_grd_ctl,'(A32)')  'TITLE SILAM GrADS anygrid output'
    WRITE(u_grd_ctl,'(A21)')  'OPTIONS LITTLE_ENDIAN'
    WRITE(u_grd_ctl,'(A6,1X,E15.6)') 'UNDEF ', real_missing
    WRITE(u_grd_ctl,'(A5,I5,A8,2(F15.7,1x))') 'XDEF ',grid%ag%nx,' LINEAR ', &
                                                & pAnyGrdParam(grid%ag%indParam)%xc(1), &
                                                & fu_dx_m_to_deg(pAnyGrdParam(grid%ag%indParam)%dx(1), &
                                                               & pAnyGrdParam(grid%ag%indParam)%yc(1))
    WRITE(u_grd_ctl,'(A5,I5,A8,2(F15.7,1x))') 'YDEF ',grid%ag%ny,' LINEAR ', &
                                                & pAnyGrdParam(grid%ag%indParam)%yc(1), &
                                                & fu_dy_m_to_deg(pAnyGrdParam(grid%ag%indParam)%dy(1))
    WRITE(u_grd_ctl,'(A15)') 'ZDEF 1 LEVELS 0'
    WRITE(u_grd_ctl,'(A16,A20,A4)') 'TDEF  1  LINEAR ', &
                             & '00:00Z01JAN2000', ' 1hr'
    WRITE(u_grd_ctl,'(A5,I5)')'VARS ',nVars

    u_grd_bin = fu_next_free_unit()
    call open_grads_binary_o(fNm, u_grd_bin, nPoints)
    rec = 1

    ! lon
    call write_grads_field(pAnyGrdParam(grid%ag%indParam)%xc, nPoints, rec, u_grd_bin)
    if(error)return
    WRITE(u_grd_ctl,'(A,A12,A,2x,A)') 'lon             ','0 99 99 0 ', 'longitude of the grid cell'
    rec = rec+1
   ! lat
    call write_grads_field(pAnyGrdParam(grid%ag%indParam)%yc, nPoints, rec, u_grd_bin)
    if(error)return
    WRITE(u_grd_ctl,'(A,A12,A,2x,A)') 'lat             ','0 99 99 0 ', 'latitude of the grid cell'
    rec = rec+1
   ! dx
    if(ifdx)then
      call write_grads_field(pAnyGrdParam(grid%ag%indParam)%dx, nPoints, rec, u_grd_bin)
      if(error)return
      WRITE(u_grd_ctl,'(A,A12,A,2x,A)') 'dx           ','0 99 99 0 ', 'x-size of the grid cell m'
      rec = rec+1
    endif
   ! dy
    if(ifdy)then
      call write_grads_field(pAnyGrdParam(grid%ag%indParam)%dy, nPoints, rec, u_grd_bin)
      if(error)return
      WRITE(u_grd_ctl,'(A,A12,A,2x,A)') 'dy           ','0 99 99 0 ', 'y-size of the grid cell m'
      rec = rec+1
    endif
   ! sin_map_rot
    if(ifsin)then
      call write_grads_field(pAnyGrdParam(grid%ag%indParam)%sin_map_rot, nPoints, rec, u_grd_bin)
      if(error)return
      WRITE(u_grd_ctl,'(A,A12,A,2x,A)') 'sin_map_rot  ','0 99 99 0 ', 'sine of local map rotation'
      rec = rec+1
    endif
   ! cos_map_rot
    if(ifcos)then
      call write_grads_field(pAnyGrdParam(grid%ag%indParam)%cos_map_rot, nPoints, rec, u_grd_bin)
      if(error)return
      WRITE(u_grd_ctl,'(A,A12,A,2x,A)') 'cos_map_rot  ','0 99 99 0 ', 'cosine of local map rotation'
      rec = rec+1
    endif

    WRITE(u_grd_ctl,'(A7)') 'ENDVARS'
    close(u_grd_ctl)
    close(u_grd_bin)

  end subroutine write_anygrid_2_gradsfile

  !***************************************************************


  subroutine rd_grd_pars_from_gradsFile(ctlFNm, grid)
    implicit none

    character(len=fnlen), intent(in) :: ctlFNm
    type(silja_grid), intent(in) :: grid

    type(silam_sp) :: sp, sp_u_case, spTmp
    integer :: nPoints, iUnit, unit_bin, iVar, nVars, iTmp
    logical :: eof
    real :: ftmp
    real, dimension(:), pointer :: dat
    character(len=fnlen) :: binFNm, chTmp
    character(len=11) :: varNm

    iUnit = fu_next_free_unit()
    IF(error)RETURN
    
    open(iUnit, file=ctlFNm, action='read', status='old', iostat=iTmp)
    if(iTmp /= 0)then
      call set_error(fu_connect_strings('Failed to open the ctl file:', ctlFNm), &
                   & 'rd_grd_pars_from_gradsFile')
      return
    endif
    nPoints = grid%ag%nx*grid%ag%ny
    sp%sp => fu_work_string()
    sp_u_case%sp => fu_work_string()
    spTmp%sp => fu_work_string()
    dat => fu_work_array()
    eof = .false.
   
    do while(.not.eof)
      call next_line_from_input_file(iUnit, sp%sp, eof)
      if(error.or.eof)exit
      sp_u_case%sp = fu_str_u_case(sp%sp)

      if(index(sp_u_case%sp,'DSET') == 1)then
        chTmp = trim(adjustl(sp%sp(6:)))
        binFNm = fu_extend_grads_hat(chTmp,ctlFNm)
        unit_bin = fu_next_free_unit()
        call open_grads_binary_i(binFNm, unit_bin, nPoints, .true.)

      elseif(index(sp_u_case%sp,'XDEF') == 1)then
        read(unit=sp%sp,fmt=*,iostat=iTmp) spTmp%sp, iVar, spTmp%sp, ftmp, ftmp ! to check
        if(iVar /= grid%ag%nx)then
          call set_error('X-dimension of grid is incorrect','rd_grd_pars_from_gradsFile')
          return
        endif     

      elseif(index(sp_u_case%sp,'YDEF') == 1)then
        read(unit=sp%sp,fmt=*,iostat=iTmp) spTmp%sp, iVar, spTmp%sp, ftmp, ftmp ! to check
        if(iVar /= grid%ag%ny)then
          call set_error('Y-dimension of grid is incorrect','rd_grd_pars_from_gradsFile')
          return
        endif     

      elseif(index(sp_u_case%sp,'VARS') == 1)then
        read(unit=sp%sp,fmt=*,iostat=iTmp) spTmp%sp, nvars
        do iVar = 1, nvars
          call next_line_from_input_file(iUnit, sp%sp, eof)
          if(error .or.eof)then
            call set_error('Failed to read all ctl variables','rd_grd_pars_from_gradsFile')
            return
          endif

          varNm = sp%sp(1:11)
          read(unit_bin, rec=iVar)(dat(iTmp),iTmp = 1, nPoints)
          if(error)then
            call set_error(fu_connect_strings('error reaading variable:', varNm),'rd_grd_pars_from_gradsFile')
          endif
          call setAnygridParam(grid, trim(adjustl(varNm)), dat)
          if(error)return

        end do  ! cycle over variables
      endif
    enddo

    call free_work_array(sp%sp)
    call free_work_array(sp_u_case%sp)
    call free_work_array(spTmp%sp)
    call free_work_array(dat)

  end subroutine rd_grd_pars_from_gradsFile

  !***************************************************************

  subroutine grid_tst()
    
    implicit none

    type(silam_area) :: area
    type(silja_grid) :: grid1, grid2

  
    
  
  
  
  end subroutine grid_tst


END MODULE grids_geo

