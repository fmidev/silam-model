MODULE silam_levels

  ! Description:
  ! Provides the definition and tools for vertical levels.
  !
  ! IN GENERAL WE DEFINE HERE:
  ! WORDS LOW AND LOWER MEANS CLOSER TO THE GROUND, NOT LOWER
  ! PRESSURE (which is higher from ground).
  ! THE OPPOSITE WORD FOR THESE ARE UP AND UPPER (not high and higher
  ! since they get more easily mixed with higher pressure)
  !
  ! Author: Mika Salonoja, email Mika.Salonoja@fmi.fi
  !
  ! All units: SI
  ! 
  ! All degrees: real numbers in degrees (degrees and parts per
  ! hundred/parts per thousand , NOT in degrees and minutes
  ! 
  ! Language: ANSI standard Fortran 90
  ! 
  ! Modules used:
!  USE silam_units
  use geography_tools
  use names_of_quantities
  use thermodynamic_tools
  
!  USE positions

  IMPLICIT NONE

  ! The public functions and subroutines available in this module:
  PUBLIC fu_set_pressure_level
  PUBLIC fu_set_constant_height_level
  PUBLIC fu_set_constant_altitude_level
  PUBLIC fu_set_sigma_level
  PUBLIC fu_set_hybrid_level
  public fu_set_depth_level
  public fu_set_layer_between_two
  public fu_set_level    ! General level. Requires the level type
  public set_named_level_with_fract ! From namelist item
  public set_vertical
  public fu_if_level_meteo_dependent
  public arrange_levels_in_vertical
  public set_missing
  public create_levels_from_namelist_v2
  public fu_central_level_of_layer
  PUBLIC defined
  PUBLIC report
  public report_as_namelist
  public level_to_short_string
  PUBLIC fu_leveltypes_eq
  PUBLIC fu_leveltype
  PUBLIC fu_pr_level_pressure
  PUBLIC fu_hybrid_level_number
  PUBLIC fu_hybrid_level_coeff_a
  PUBLIC fu_hybrid_level_coeff_b
  PUBLIC fu_hybrid_level_pressure
  PUBLIC fu_level_height
  PUBLIC fu_level_altitude
  public fu_sigma_level_sigma
  public fu_depth_level_depth
  public fu_std_pressure
  public fu_str2leveltype

  PUBLIC fu_level_count
  public fu_NbrOfLevels
  public fu_level
  public fu_level_nbr
  public set_level    ! in vertical
  public add_level
  public fu_verts_comparable
  public fu_range
  public fu_level_belongs_to_vertical
  public fu_level_index  ! Index of the level in some vertical structure
  public fu_level_to_vertical_crude  ! Crudely project level to vertical (interpolation possible !)
  public fu_level_to_vertical  ! Same but using actual meteodata
  public change_vertical_type
  public remove_last_level
  public fu_if_layer  ! thick layer of a single level ?
  public fu_lower_boundary_of_layer
  public fu_upper_boundary_of_layer
  public fu_layer_thickness_m         ! In metres !!
  public fu_layer_thickness_local_unit
  public fu_vert_overlap_fraction   ! Of one thick layer with another one
  public overlap_fraction_lyr_in_vert  ! overlap array for verticals
  public vertical_parameters
  public vertical_parameters_ab
  public make_vertical_of_levels
  public make_vertical_of_borders
  public make_vertical_of_layers 
  public fu_if_cut_vertical_size_from_above

  public fu_level_value           ! A quick function that delivers the VALUE of the level, local unit
  public fu_top_of_layer_value    ! A quick function that delivers the VALUE of the layer top, local unit
  public fu_bottom_of_layer_value
  public fu_layer_centre_value
  public reproject_verticals
  public hybrid_coefs

  public fu_project_level_crude
  public fu_project_level
  public projection_input_needs
  public vert_interp_data_crude
  public fu_height_for_press

  public fu_cmp_levs_eq
  public fu_cmp_verts_eq

  public vert_to_metric
  public test_vert_to_metric
  public vert_to_pressure
  public test_vert_to_pressure
  
  ! The private functions and subroutines not to be used elsewhere:
!  private vertical_layer_thicknesses_m

  private fu_set_level_general
  
  private set_vertical_single_level
  private set_vertical_multi_level
  private set_vertical_copy
  private set_vertical_from_namelist
  private set_level_in_vertical
  private set_missing_vertical

  private fu_first_lower_sametype
  private fu_first_upper_sametype
  private fu_first_lower_eq_sametype
  private fu_first_upper_eq_sametype

  PRIVATE fu_level_defined
  PRIVATE fu_vertical_defined
  PRIVATE fu_set_hybrid_level_no_coeff
  PRIVATE fu_set_hybrid_level_with_coeff
  PRIVATE fu_leveltype_of_level
  private fu_leveltype_of_vertical 
  private fu_set_layer_btw_two_ordinary
  private fu_set_layer_btw_two_general
  private fu_set_layer_between_two_hybrid

  PRIVATE print_level_report
  PRIVATE print_vertical_report
  private report_vertical_as_namelist
  private write_vert_and_fract_as_nmLst
  private write_level_and_fract_as_nmLst
  private fu_Nlevs_of_vertical
  private fu_level_from_vertical
  private fu_level_nbr_from_vertical
  private fu_range_of_vertical
  private fu_top_level_vert
  private fu_bottom_level_vert
  private fu_add_levels
  private fu_multiply_level
  private fu_level_index_in_levels
  private fu_level_index_in_vertical

  ! A quick function that delivers the VALUE of the layer top, local unit
  private fu_top_of_layer_value_lyr    
  ! A quick function that delivers the VALUE of the layer top, local unit
  private fu_top_of_layer_value_vert    
  private fu_bottom_of_layer_value_lyr
  private fu_bottom_of_layer_value_vert
  private fu_layer_centre_value_lyr
  private fu_layer_centre_value_vert
  private level_tests
  private test_overlap_fraction
  
  private fu_project_level_single
  private fu_project_level_from_vertical

  ! Generic names and operator-interfaces of some functions:

  INTERFACE fu_set_hybrid_level
    MODULE PROCEDURE fu_set_hybrid_level_no_coeff, fu_set_hybrid_level_with_coeff
  END INTERFACE

  interface fu_set_layer_between_two
    module procedure fu_set_layer_btw_two_ordinary
    module procedure fu_set_layer_btw_two_general
    module procedure fu_set_layer_between_two_hybrid
  end interface

  interface set_vertical
    module procedure set_vertical_single_level
    module procedure set_vertical_multi_level
!    module procedure set_vertical_copy
    module procedure set_vertical_from_namelist
  end interface

  interface set_missing
    module procedure set_missing_vertical
  end interface

  interface fu_set_level
    module procedure fu_set_level_general
  end interface

  interface set_level
    module procedure set_level_in_vertical
  end interface

  INTERFACE defined
    MODULE PROCEDURE fu_level_defined
    MODULE PROCEDURE fu_vertical_defined
  END INTERFACE

  INTERFACE report
    MODULE PROCEDURE print_level_report
    MODULE PROCEDURE print_vertical_report
  END INTERFACE

  INTERFACE report_as_namelist
    module procedure report_vertical_as_namelist     ! just vertical
    module procedure write_vert_and_fract_as_nmLst   ! vertical with fractions: different namelist!!
    module procedure write_level_and_fract_as_nmLst  ! level with fractions: different namelist!!
  END INTERFACE

!  INTERFACE operator(==)
!    MODULE PROCEDURE fu_cmp_levs_eq
!    MODULE PROCEDURE fu_cmp_verts_eq
!  END INTERFACE

  INTERFACE operator(<) ! same type levels only!!!!
    MODULE PROCEDURE fu_first_lower_sametype
  END INTERFACE

  INTERFACE operator(>) ! same type levels only!!!!
    MODULE PROCEDURE fu_first_upper_sametype
  END INTERFACE

  INTERFACE operator(<=) ! same type levels only!!!!
    MODULE PROCEDURE fu_first_lower_eq_sametype
  END INTERFACE

  INTERFACE operator(>=) ! same type levels only!!!!
    MODULE PROCEDURE fu_first_upper_eq_sametype
  END INTERFACE

  interface operator (+)
    module procedure fu_add_levels
  end interface

  interface operator (*)
    module procedure fu_multiply_level
  end interface

  interface assignment (=)
    module procedure set_vertical_copy
  end interface

  INTERFACE fu_leveltype
    MODULE PROCEDURE fu_leveltype_of_level
    MODULE PROCEDURE fu_leveltype_of_vertical
  END INTERFACE

  interface fu_NbrOfLevels
    module procedure fu_Nlevs_of_vertical 
  end interface

  interface fu_level
    module procedure fu_level_from_vertical
  end interface

  interface fu_level_nbr
    module procedure fu_level_nbr_from_vertical
  end interface

  interface fu_range
    module procedure fu_range_of_vertical
  end interface

  interface fu_level_index
    module procedure fu_level_index_in_levels
    module procedure fu_level_index_in_vertical
  end interface

  interface fu_top_of_layer_value
    module procedure fu_top_of_layer_value_lyr
    module procedure fu_top_of_layer_value_vert
  end interface

  interface fu_bottom_of_layer_value
    module procedure fu_bottom_of_layer_value_lyr
    module procedure fu_bottom_of_layer_value_vert
  end interface

  interface fu_layer_centre_value
    module procedure fu_layer_centre_value_lyr
    module procedure fu_layer_centre_value_vert
  end interface
  
  interface fu_project_level
    module procedure fu_project_level_single
    module procedure fu_project_level_from_vertical
  end interface


  ! The flags of the different level-types to be used everywhere:
  ! NOTE: these are WMO codes used in GRIB libraries !! Do not create
  ! own stuff - consult WMO GRIB manuals !
  !
  INTEGER, PARAMETER, PUBLIC :: surface = 001
  INTEGER, PARAMETER, PUBLIC :: top_of_the_atmosphere = 008
  INTEGER, PARAMETER, PUBLIC :: constant_pressure = 100
  INTEGER, PARAMETER, PUBLIC :: layer_btw_2_pressure = 101
  INTEGER, PARAMETER, PUBLIC :: mean_sea = 102
  INTEGER, PARAMETER, PUBLIC :: constant_altitude = 103 ! Above msl
  INTEGER, PARAMETER, PUBLIC :: layer_btw_2_altitude = 104
  INTEGER, PARAMETER, PUBLIC :: constant_height = 105
  INTEGER, PARAMETER, PUBLIC :: layer_btw_2_height = 106
  INTEGER, PARAMETER, PUBLIC :: sigma_level = 107
  INTEGER, PARAMETER, PUBLIC :: layer_btw_2_sigma = 108
  INTEGER, PARAMETER, PUBLIC :: hybrid = 109
  INTEGER, PARAMETER, PUBLIC :: layer_btw_2_hybrid = 110
  INTEGER, PARAMETER, PUBLIC :: depth_level = 111
  INTEGER, PARAMETER, PUBLIC :: layer_btw_2_depth = 112

  INTEGER, PARAMETER, PUBLIC :: entire_atmosphere_single_layer = 200
  INTEGER, PARAMETER, PUBLIC :: no_level = 0
  integer, parameter, public :: any_level = 333

  !==========================================================================
  !==========================================================================
  !
  ! The main vertical characteristic: silja_level and a set of levels
  ! silam_vertical
  !
  !==========================================================================
  !==========================================================================
  TYPE silja_level
    PRIVATE
    INTEGER :: leveltype  ! allowed values above
    LOGICAL :: hybrid_coeff_known
    INTEGER :: number, number2 ! for hybrid
    REAL :: a,b, a2, b2 ! Change meaning depending on the level type
                        ! for hybrid a point (x,y,z,t) the pressure at a hybrid
                        ! level n is defined as: P(n) = a(n) + b(n)*P(x,y,t,surf)
    TYPE(silja_logical) :: defined = silja_false
  END TYPE silja_level
  !
  ! Vertical is a sorted set of same-type levels
  !
  type silam_vertical
    private
    integer :: vert_type  = int_missing  ! Type of the vertical co-ordinate
    integer :: NLevs  =int_missing  ! Number of levels in the structure
    type(silja_level), dimension(:), pointer :: levs => null() ! Levels
    TYPE(silja_logical) :: defined   = silja_false
  end type silam_vertical

  !
  ! The missing constants
  !
  TYPE(silja_level), PARAMETER, PUBLIC :: level_missing = silja_level(&
                    & no_level,&
                    & .false.,&
                    & int_missing,int_missing,&
                    & real_missing, real_missing, real_missing, real_missing,&
                    & silja_false)

  type(silam_vertical), parameter,public :: vertical_missing = &
     & silam_vertical(int_missing, 0, null(), silja_false)
!
!  type(silam_vertical), parameter,public :: vertical_all_levels = &
!     & silam_vertical(any_level, 0, level_missing)

  !-----------------------------------------------------------------------
  !
  ! There are two verticals: meteo_vertical and dispersaion_vertical
  ! meteo_vertical is made to handle the meteorological fields, make derived ones
  ! and serves as the vertical coordinate for the Lagrangian environment
  ! while dispersion_vertical serves the Eulerian environment.
  !
  !-----------------------------------------------------------------------
  !
  type(silam_vertical),public,target,save :: meteo_vertical = vertical_missing
  type(silam_vertical),public,pointer,save :: meteo_verticalPtr => null()
  integer, public, save :: nz_meteo  ! Vertical dimension of the meteo_vertical
  real, dimension(:), public, pointer, save :: a_met, b_met !(0:nz_meteo+1)
  ! to be set by vertical_parameters_ab

  type(silam_vertical),public,target,save :: dispersion_vertical = vertical_missing
  type(silam_vertical),public,pointer,save :: dispersion_verticalPtr => null()
  integer, public, save :: nz_dispersion ! Vertical dimension of the dispersion_vertical
  real, dimension(:), public, pointer, save :: a_half_disp => null(), b_half_disp => null()  !(0:nz+2)
  real, dimension(:), public, pointer, save :: disp_layer_top_m => null() !(-1:nz+1)
  ! set by vertical_parameters below, called from io_init

  type(silam_vertical),public,target,save :: output_vertical = vertical_missing
  type(silam_vertical),public,pointer,save :: output_verticalPtr  => null()
  integer, public, save :: nz_output ! Vertical dimension of the dispersion_vertical
!  either  a_half_XXX and b_half_XXX, or out_layer_top_m  is null
  real, dimension(:), public, pointer, save :: a_half_out, b_half_out  !Set from io_server
  real, dimension(:), public, pointer, save :: out_layer_top_m
!  real, dimension(:), public, pointer, save :: output_layer_boundaries

  type(silam_vertical), private, save :: vertical_1000hpa = vertical_missing


  integer, public, save :: indexTemperature, indexSurfPressure, indexRelief

  !
  ! Parameters that control actions with the verticals: take as-is (may be, interpolate)
  ! integrate over it, or average.
  !
  integer, public, parameter :: do_nothing_flag  = 530, &
                              & integrate_column_flag = 531, &
                              & lowest_level_flag = 533, &
                              & level_3D_type_flag = 534, &
                              & level_2d_type_flag = 535

  ! Some parameters:

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_1000hpa = &
      & silja_level(constant_pressure,  .false.,&
      & int_missing,int_missing,&
      & 100000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_950hpa = &
      & silja_level(constant_pressure,  .false.,&
      & int_missing,int_missing,&
      & 95000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_925hpa = &
      & silja_level(constant_pressure,  .false.,&
      & int_missing,int_missing,&
      & 92500., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_900hpa = &
      & silja_level(constant_pressure, .false.,&
      & int_missing,int_missing,&
      & 90000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_850hpa = &
      & silja_level(constant_pressure, .false.,&
      & int_missing,int_missing,&
      & 85000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_800hpa = &
      & silja_level(constant_pressure, .false.,&
      & int_missing,int_missing,&
      & 80000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_750hpa = &
      & silja_level(constant_pressure, .false.,&
      & int_missing,int_missing,&
      & 75000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_700hpa = &
      & silja_level(constant_pressure, .false.,&
      & int_missing,int_missing,&
      & 70000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_650hpa = &
      & silja_level(constant_pressure, .false.,&
      & int_missing,int_missing,&
      & 65000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_600hpa = &
      & silja_level(constant_pressure, .false.,&
      & int_missing,int_missing,&
      & 60000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_550hpa = &
      & silja_level(constant_pressure, .false.,&
      & int_missing,int_missing,&
      & 55000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_500hpa = &
      & silja_level(constant_pressure, .false.,&
      & int_missing,int_missing,&
      & 50000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_450hpa = &
      & silja_level(constant_pressure, .false.,&
      & int_missing,int_missing,&
      & 45000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_400hpa = &
      & silja_level(constant_pressure, .false.,&
      & int_missing,int_missing,&
      & 40000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_350hpa = &
      & silja_level(constant_pressure, .false.,&
      & int_missing,int_missing,&
      & 35000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_300hpa = &
      & silja_level(constant_pressure, .false.,&
      & int_missing,int_missing,&
      & 30000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_250hpa = &
      & silja_level(constant_pressure, .false.,&
      & int_missing,int_missing,&
      & 25000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_200hpa = &
      & silja_level(constant_pressure, .false.,&
      & int_missing,int_missing,&
      & 20000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_150hpa = &
      & silja_level(constant_pressure, .false.,&
      & int_missing,int_missing,&
      & 15000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_100hpa = &
      & silja_level(constant_pressure, .false.,&
      & int_missing,int_missing,&
      & 10000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_90hpa = &
      & silja_level(constant_pressure, .false.,&
      & int_missing,int_missing,&
      & 9000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_80hpa = &
      & silja_level(constant_pressure, .false.,&
      & int_missing,int_missing,&
      & 8000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_70hpa = &
      & silja_level(constant_pressure, .false.,&
      & int_missing,int_missing,&
      & 7000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_60hpa = &
      & silja_level(constant_pressure, .false.,&
      & int_missing,int_missing,&
      & 6000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_50hpa = &
      & silja_level(constant_pressure, .false.,&
      & int_missing,int_missing,&
      & 5000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_40hpa = &
      & silja_level(constant_pressure, .false.,&
      & int_missing,int_missing,&
      & 4000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_30hpa = &
      & silja_level(constant_pressure, .false.,&
      & int_missing,int_missing,&
      & 3000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_20hpa = &
      & silja_level(constant_pressure, .false.,&
      & int_missing,int_missing,&
      & 2000., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: pr_level_10hpa = &
      & silja_level(constant_pressure, .false.,&
      & int_missing,int_missing,&
      & 1000., real_missing, real_missing, real_missing,&
      & silja_true)


  ! ***** Height and altitude levels

  TYPE(silja_level), PARAMETER, PUBLIC :: ground_level = &
      & silja_level(constant_height, .false.,&
      & int_missing,int_missing,&
      & 0., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: surface_level = &
      & silja_level(surface, .false.,&
      & int_missing,int_missing,&
      & real_missing, real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: level_2m_above_ground = &
      & silja_level(constant_height, .false.,&
      & int_missing,int_missing,&
      & 2., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: level_10m_above_ground = &
      & silja_level(constant_height, .false.,&
      & int_missing,int_missing,&
      & 10., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: mean_sea_level = &
      & silja_level(constant_altitude,  .false.,&
      & int_missing,int_missing,&
      & 0., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: level_999m_above_ground = &
      & silja_level(constant_height, .false.,&
      & int_missing,int_missing,&
      & 999., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: level_998m_above_ground = &
      & silja_level(constant_height, .false.,&
      & int_missing,int_missing,&
      & 998., real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: top_atmosphere_level = &
      & silja_level(top_of_the_atmosphere, .false.,&
      & int_missing,int_missing,&
      & real_missing, real_missing, real_missing, real_missing,&
      & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: entire_atmosphere_mean_level = &
      & silja_level(entire_atmosphere_single_layer, .false.,&
                  & 1001, int_missing,&
                  & real_missing, real_missing, real_missing, real_missing,&
                  & silja_true)

  TYPE(silja_level), PARAMETER, PUBLIC :: entire_atmosphere_integr_level = &
      & silja_level(entire_atmosphere_single_layer, .false.,&
                  & 1002, int_missing,&
                  & real_missing, real_missing, real_missing, real_missing,&
                  & silja_true)
  
  TYPE(silja_level), PARAMETER, PUBLIC :: lowest_atmosphere_level = &
      & silja_level(entire_atmosphere_single_layer, .false.,&  ! FIXME: Not really true
      & 1003,int_missing,&
      & real_missing, real_missing, real_missing, real_missing,&
      & silja_true)



CONTAINS

  !*****************************************************************

  function fu_set_level_general(level_type, fVal1, fVal2, fVal3, fVal4, iVal1, iVal2) result(level)
    ! 
    ! Sets one level of the given type. Meaning of all values depends on the level_type
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    type(silja_level) :: level

    ! Imported parameters with intent(in):
    integer, intent(in) :: level_type 
    real, intent(in), optional :: fVal1, fVal2, fVal3, fVal4
    integer, intent(in), optional :: iVal1, iVal2

    level = level_missing

    select case(level_type)
      case(surface)
        level = surface_level

      case(top_of_the_atmosphere)
        level = top_atmosphere_level

      case(constant_pressure)
        if(present(fVal1))then
          level = fu_set_pressure_level(fVal1)
        else
          call set_error('fVal1 is needed for pressure level','fu_set_level_general')
        endif

      case(layer_btw_2_pressure, layer_btw_2_altitude, &
         & layer_btw_2_height, layer_btw_2_sigma, layer_btw_2_depth)
        if(present(fVal1) .and. present(fVal2))then
          level = fu_set_layer_between_two(level_type, fVal1,fVal2)
        else
          call set_error('fVal1 and fVal2 are needed for layer between 2','fu_set_level_general')
        endif

      case(mean_sea)
        level = mean_sea_level

      case(constant_altitude)
        if(present(fVal1))then
          level = fu_set_constant_altitude_level(fVal1)
        else
          call set_error('fVal1 is needed for altitude level','fu_set_level_general')
        endif

      case(constant_height)
        if(present(fVal1))then
          level = fu_set_constant_height_level(fVal1)
        else
          call set_error('fVal1 is needed for height level','fu_set_level_general')
        endif

      case(sigma_level)
        if(present(fVal1))then
          level = fu_set_sigma_level(fVal1)
        else
          call set_error('fVal1 is needed for sigma level','fu_set_level_general')
        endif

      case(hybrid)
        if(present(iVal1))then
          if(present(fVal1) .and. present(fVal2))then
            level = fu_set_hybrid_level_with_coeff(iVal1, fVal1, fVal2)
          else
            level = fu_set_hybrid_level_no_coeff(iVal1)
          endif
        else
          if(present(fVal1) .and. present(fVal2))then
            level = fu_set_hybrid_level_with_coeff(0, fVal1, fVal2)
          else
            call set_error('Number or coefs are needed for hybrid level','fu_set_level_general')
          endif
        endif

      case(layer_btw_2_hybrid)
        if(present(fVal1) .and. present(fVal2) .and. &
         & present(fVal3) .and. present(fVal4) .and. &
         & present(iVal1) .and. present(iVal2))then
          level = fu_set_layer_between_two_hybrid(level_type,iVal2,iVal1,fVal1,fVal2, fVal3, fVal4)
        else
          call set_error('All values are needed for hybrid layer btw 2','fu_set_level_general')
        endif

      case(entire_atmosphere_single_layer)
        level = entire_atmosphere_mean_level

      case(no_level)  ! Do nothing

      case(any_level)
        call set_error('Level type must be defined','fu_set_level_general')

      case default
        call set_error('Unknown level type','fu_set_level_general')

    end select

  end function fu_set_level_general



  ! ***************************************************************


  FUNCTION fu_set_pressure_level(pressure) result(level)

    ! Description:
    ! Sets values for one constant pressure level.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_level) :: level

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: pressure ! pascals, not HPa

    level = level_missing
    level%leveltype = constant_pressure
    level%a = pressure
    level%defined = fu_set_true()

  END FUNCTION fu_set_pressure_level


  ! ***************************************************************


  FUNCTION fu_set_constant_height_level(height) result(level)

    ! Description:
    ! Sets values for one constant height above ground level.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_level) :: level

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: height

    level = level_missing
    level%leveltype = constant_height
    level%a = height
    level%defined = fu_set_true()

  END FUNCTION fu_set_constant_height_level



  ! ***************************************************************


  FUNCTION fu_set_constant_altitude_level(altitude) result(level)

    ! Description:
    ! Sets values for one constant altitude above mean sea level.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_level) :: level

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: altitude

    level = level_missing
    level%leveltype = constant_altitude
    level%a = altitude
    level%defined = fu_set_true()

  END FUNCTION fu_set_constant_altitude_level



  ! ***************************************************************


  FUNCTION fu_set_sigma_level(sigma_coefficient) result(level)

    ! Description:
    ! Sets values for one sigma-level, on which in all point pressure
    ! is defined as ground_pressure * sigma_coefficient.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_level) :: level

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: sigma_coefficient

    level = level_missing
    level%leveltype = sigma_level
    level%a = sigma_coefficient
    level%defined = fu_set_true()

  END FUNCTION fu_set_sigma_level



  ! ***************************************************************


  FUNCTION fu_set_hybrid_level_no_coeff(level_number) result(level)

    ! Description:
    ! Sets values for one hybrid-level, on which in all points pressure
    ! is defined as following:
    ! P(level) = a(level) + b(level)*ground_pressure
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_level) :: level

    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: level_number

    level = level_missing

    level%leveltype = hybrid
    level%number = level_number
    level%hybrid_coeff_known = .false.
    level%defined = fu_set_true()

  END FUNCTION fu_set_hybrid_level_no_coeff



  ! ***************************************************************


  FUNCTION fu_set_hybrid_level_with_coeff(level_number, a, b) result(level)

    ! Description:
    ! Sets values for one hybrid-level, on which in all points pressure
    ! is defined as following:
    ! P(level) = a + b*ground_pressure in which a and b depend only
    ! on level number
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_level) :: level

    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: level_number
    REAL, INTENT(in) :: a, b

    level = level_missing
    level%leveltype = hybrid
    level%number = level_number
    level%hybrid_coeff_known = .true.
    level%a = a
    level%b = b
    level%defined = fu_set_true()

  END FUNCTION fu_set_hybrid_level_with_coeff

  
  ! ***************************************************************


  FUNCTION fu_set_depth_level(depth) result(level)
    ! 
    ! Sets values for one depth level. Note that depth is a positive number, i.e. the z axis
    ! is directed downwards. Looks like a GRIB stadard
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_level) :: level

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: depth

    level = level_missing
    level%leveltype = depth_level
    level%a = depth
    level%defined = fu_set_true()

  END FUNCTION fu_set_depth_level


  ! ***************************************************************


  FUNCTION fu_set_layer_btw_two_general(top, bottom) result(level)

    ! Sets values for one thick layer between some two levels of any type.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_level) :: level

    ! Imported parameters with intent(in):
    type(silja_level), intent(in) :: top, bottom

    if(.not. defined(top) .or. .not. defined(bottom))then
      call set_error('Undefined top or bottom given', 'fu_set_layer_btw_two_general')
      return
    endif
    if(top%leveltype /= bottom%leveltype)then
      call set_error('Different level types of top and bottom', &
                   & 'fu_set_layer_btw_two_general')
      return
    endif
    select case(top%leveltype)
      case (constant_pressure)

        level%leveltype = layer_btw_2_pressure
        level%a = min(top%a,bottom%a)
        level%b = max(top%a,bottom%a)
        level%defined = silja_true
        level%hybrid_coeff_known = .false.
        level%number = int_missing
        level%number2 = int_missing

      case (constant_altitude)

        level%leveltype = layer_btw_2_altitude
        level%a = max(top%a,bottom%a)
        level%b = min(top%a,bottom%a)
        level%defined = silja_true
        level%hybrid_coeff_known = .false.
        level%number = int_missing
        level%number2 = int_missing

      case (constant_height)

        level%leveltype = layer_btw_2_height
        level%a = max(top%a,bottom%a)
        level%b = min(top%a,bottom%a)
        level%defined = silja_true
        level%hybrid_coeff_known = .false.
        level%number = int_missing
        level%number2 = int_missing

      case (sigma_level)

        level%leveltype = layer_btw_2_sigma
        level%a = min(top%a,bottom%a)
        level%b = max(top%a,bottom%a)
        level%defined = silja_true
        level%hybrid_coeff_known = .false.
        level%number = int_missing
        level%number2 = int_missing

      case (depth_level)

        level%leveltype = layer_btw_2_depth
        level%a = max(top%a,bottom%a)
        level%b = min(top%a,bottom%a)
        level%defined = silja_true
        level%hybrid_coeff_known = .false.
        level%number = int_missing
        level%number2 = int_missing

      case (hybrid)
        level%leveltype = layer_btw_2_hybrid
        
        !if(bottom < top)then
        level%a = top%a
        level%b = top%b
        level%a2 = bottom%a
        level%b2 = bottom%b
!!$        else
!!$          level%a2 = top%a
!!$          level%b2 = top%b
!!$          level%a = bottom%a
!!$          level%b = bottom%b
!!$        endif
        level%defined = silja_true
        level%hybrid_coeff_known = top%hybrid_coeff_known .and. bottom%hybrid_coeff_known
        level%number = top%number
        level%number2 = bottom%number
        
      case default
        call set_error('Strange level type','fu_set_layer_btw_two_general')
        return
    end select

  END FUNCTION fu_set_layer_btw_two_general



  ! ***************************************************************


  FUNCTION fu_set_layer_btw_two_ordinary(lev_typ, top, bottom) result(level)

    ! Sets values for one thick layer between some two levels of any type,
    ! except for the hybrid ones.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_level) :: level

    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: lev_typ
    real, intent(in) :: top, bottom

    select case(lev_typ)
    case (layer_btw_2_pressure, & ! = 101
        & layer_btw_2_sigma)    != 112

      level%leveltype = lev_typ
      level%a = min(top, bottom)
      level%b = max(top, bottom)
      level%defined = silja_true
      level%hybrid_coeff_known = .false.
      level%number = int_missing
      level%number2 = int_missing

    case (layer_btw_2_altitude, & ! = 104
        & layer_btw_2_height, & ! = 106
        & layer_btw_2_depth)    != 112

      level%leveltype = lev_typ
      level%a = max(top, bottom)
      level%b = min(top, bottom)
      level%defined = silja_true
      level%hybrid_coeff_known = .false.
      level%number = int_missing
      level%number2 = int_missing
    case default
      call set_error('Strange level type','fu_set_layer_btw_two_ordinary')
      return
    end select

  END FUNCTION fu_set_layer_btw_two_ordinary



  ! ***************************************************************


  FUNCTION fu_set_layer_between_two_hybrid(lev_typ,top,bottom,a,b, a2, b2)result(level)

    ! Sets values for one thick hybrid-level located between two others, 
    ! on those the pressure is defined as following:
    ! P(level) = a + b*ground_pressure in which a and b depend only
    ! on level number.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_level) :: level

    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: lev_typ, top, bottom
    REAL, INTENT(in) :: a, b, a2, b2

    level = level_missing
    if(lev_typ /= layer_btw_2_hybrid)then
      call set_error('Non-hybrid layer','fu_set_level_between_two_hybrid')
      return
    end if
    level%leveltype = layer_btw_2_hybrid
    level%number = top
    level%number2 = bottom
    level%hybrid_coeff_known = .true.
    level%a = a
    level%b = b
    level%a2 = a2
    level%b2 = b2
    level%defined = fu_set_true()
!    if ((a2-a) + std_pressure_sl*(b2-b) < 0)  call ooops("fu_set_layer_between_two_hybrid")

  END FUNCTION fu_set_layer_between_two_hybrid


  integer function fu_str2leveltype(str) result (lev_typ)

    IMPLICIT NONE
    character (len=*), intent(in) :: str

    select case (str) 

      case('*')
         lev_typ = any_level
      case( 'ALTITUDE_FROM_SEA_LEVEL')
         lev_typ = constant_altitude
      case( 'HEIGHT_FROM_SURF')
         lev_typ = constant_height
      case( 'PRESSURE')
         lev_typ = constant_pressure
      case( 'DEPTH')
         lev_typ = depth_level
      case( 'DEPTH_LAYER')
         lev_typ = layer_btw_2_depth
      case( 'ENTIRE_ATMOSPHERE_INTEGR_LAYER')
         lev_typ = entire_atmosphere_single_layer
      case( 'HYBRID')
         lev_typ = hybrid
      case( 'ALTITUDE_LAYER')
         lev_typ = layer_btw_2_altitude
      case( 'HEIGHT_LYR_FROM_SURF')
         lev_typ = layer_btw_2_height
      case( 'HYBRID_LAYER')  
         lev_typ = layer_btw_2_hybrid
      case( 'PRESSURE_LAYER')  
         lev_typ = layer_btw_2_pressure
      case( 'SIGMA_LAYER')  
         lev_typ = layer_btw_2_sigma
      case( 'MEAN_SEA_LEVEL')
         lev_typ = mean_sea
      case( 'SIGMA_LEVEL')
         lev_typ = sigma_level
      case( 'SURFACE_LEVEL')
         lev_typ = surface
      case( 'TOP_ATMOSPHERE_LEVEL')
         lev_typ = top_of_the_atmosphere
      case( 'XXX')
         lev_typ = no_level
      case default
         call msg_warning("Unknown leveltype string: '"//str//"'","fu_str2leveltype")
         lev_typ = int_missing
    end select
    end function fu_str2leveltype



  !*****************************************************************
  !*****************************************************************
  !
  !   Vertical-related routines
  !
  !*****************************************************************
  !*****************************************************************

  !*****************************************************************

  subroutine set_missing_vertical(vertical, ifNew)
    !
    ! Replacement for the vertical_missing parameter
    !
    implicit none

    ! Imported parameters
    type(silam_vertical), intent(inout) :: vertical
    logical, intent(in) :: ifNew

    ! Local variables
    integer :: iStat

    if(.not. ifNew)then
      if(vertical%defined == silja_true)then
        if(vertical%nLevs > 0 .and. vertical%nLevs < 10000)then
           if(associated(vertical%levs)) then
              deallocate(vertical%levs, stat=iStat)
           end if
        endif
      endif
    endif

    vertical%nLevs = int_missing
    vertical%vert_type = int_missing
    nullify(vertical%levs)
    vertical%defined = silja_false

  end subroutine set_missing_vertical


  !****************************************************************

  subroutine set_vertical_single_level(level, vertical)
    !
    ! Sets a new vertical structure copying its main parameters from 
    ! the given level
    !
    implicit none

    ! Imported parameters
    type(silja_level), intent(in) :: level
    type(silam_vertical), intent(out) :: vertical

    ! Local variables
    integer :: iStat

    vertical%defined = silja_false

    if(.not.defined(level))then
      call set_error('Undefined level given','set_vertical_single_level')
      return
    end if

    allocate(vertical%levs(1), stat=iStat)
    if(iStat /= 0)then
      call set_missing(vertical, .true.)
      call set_error('Failed to allocate level','set_vertical_single_level')
      return
    endif
    vertical%levs(1) = level
    vertical%vert_type = fu_leveltype(level)
    vertical%Nlevs = 1

    vertical%defined = silja_true

  end subroutine set_vertical_single_level


  !****************************************************************

  subroutine set_vertical_multi_level(levels, vertical)
    !
    ! Sets a new vertical structure copying its main parameters from 
    ! the given level
    !
    implicit none

    ! Imported values with intent IN
    type(silja_level), dimension(:), intent(in) :: levels
    type(silam_vertical) :: vertical

    ! local variables
    integer :: iLev

    vertical%defined = silja_false

    if(.not.defined(levels(1)))then
      call set_error('Undefined first level given','set_vertical_multi_level')
      return
    end if

    vertical%vert_type = fu_leveltype(levels(1))

    do iLev = 1, size(levels)
      if(.not.defined(levels(iLev)))exit
      vertical%Nlevs = iLev
    end do

    allocate(vertical%levs(vertical%Nlevs), stat=iLev)
    if(iLev /= 0)then
      vertical%Nlevs = int_missing
      call set_missing(vertical, .true.)
      call set_error('Failed to allocate space for levels','set_vertical_multi_level')
      return
    endif

    vertical%levs(1:vertical%nLevs) = levels(1:vertical%nLevs)

    vertical%defined = silja_true

!!    call report(vertical,.true.)

  end subroutine set_vertical_multi_level


  !****************************************************************

  subroutine set_vertical_copy(vertical, vertIn)
    !
    ! Sets a new vertical structure copying its main parameters from 
    ! the given level
    !
    implicit none

    ! Imported values with intent IN
    type(silam_vertical), intent(in) :: vertIn
    type(silam_vertical), intent(out) :: vertical

    ! local variables
    integer :: i

    vertical%defined = silja_false

    if(.not.defined(vertIn))then
!      call msg_warning('Undefined input vertical given','set_vertical_copy')
      call set_missing_vertical(vertical, .true.)
      return
    end if

    allocate(vertical%levs(vertIn%nLevs), stat=i)
    if(i /= 0)then
      call set_missing(vertical, .true.)
      call set_error('Failed to allocate the space for levels','set_vertical_copy')
      return
    endif

    vertical%vert_type = vertIn%vert_type
    vertical%Nlevs = vertIn%Nlevs
    vertical%levs(1:vertIn%nLevs) = vertIn%levs(1:vertIn%nLevs)

    vertical%defined = silja_true

  end subroutine set_vertical_copy


  !****************************************************************

  subroutine set_vertical_from_namelist(nlSetup, vertical)
    !
    ! Sets a new vertical structure getting its parameters from 
    ! the namelist and using create_levels_from_namelist
    !
    implicit none

    ! Imported values with intent IN
    type(Tsilam_namelist), pointer :: nlSetup
    type(silam_vertical), target :: vertical

    ! local variables
    integer :: iLev

    vertical%defined = silja_false

    nullify(vertical%levs)
    call create_levels_from_namelist_v2(vertical%levs,nlSetup)
    if(error)return

    if(associated(vertical%levs))then
      vertical%vert_type = fu_leveltype(vertical%levs(1))
      vertical%Nlevs = size(vertical%levs)
    else
      call set_error('Failed to create vertical','set_vertical_from_namelist')
      return
    endif

    vertical%defined = silja_true
    ! 
    ! It seems to make sense to sort the vertical, since namelist does
    ! not have any order.
! MAS: commented out since the arrangement is not needed since the namelist has the layer number
!
!    call arrange_levels_in_vertical(vertical)

  end subroutine set_vertical_from_namelist

  !************************************************************************************

  logical function fu_if_level_meteo_dependent(level_type)
    !
    ! True if the position of the level of this type changes when meteorology changes
    !
    implicit none

     ! Imported parameter
     integer, intent(in) :: level_type

     select case(level_type)
       case(surface, top_of_the_atmosphere, mean_sea, &
          & constant_altitude, layer_btw_2_altitude, constant_height, layer_btw_2_height, &
          & depth_level, layer_btw_2_depth, entire_atmosphere_single_layer)

         fu_if_level_meteo_dependent = .false.

       case(constant_pressure, layer_btw_2_pressure, sigma_level, layer_btw_2_sigma, &
          & hybrid, layer_btw_2_hybrid)

         fu_if_level_meteo_dependent = .true.

       case default
         call set_error('Unknown level type:' + fu_str(level_type),'fu_if_level_meteo_dependent')
         fu_if_level_meteo_dependent = .true. ! to minimise the possible trouble

     end select

  end function fu_if_level_meteo_dependent

  !************************************************************************************
  
  subroutine create_levels_from_namelist_v2(levels, nl)
    !
    ! Create levels from namelist. Allocates the levels array as needed.
    !
    ! The levels can be set either as thick or thin layers as
    ! determined by the vertical_method item. For non-hybrid levels,
    ! the level values are given either as layer_thickness, or as
    ! midpoints (the 'levels' item), or as a set of layers with given borders. 
    ! This is regardless of the thick/thin setting.
    ! 
    ! A different scheme is needed for the hybrid levels. Thin layers
    ! require entries hybrid_coefficients in the namelist, thick
    ! layers require hybrid_coefficients_top and
    ! hybrid_coefficients_bottom. Between these, no conversion is
    ! currently implemented.

    implicit none
    type(silja_level), dimension(:), pointer :: levels
    type(Tsilam_namelist), pointer :: nl
    
    ! Local variables
    character(len=worksize_string), pointer :: level_data, layer_data
    character(len=clen) :: vertical_method, level_type, chUnit, chSIunit
    logical :: have_levels, have_layers, need_levels, need_layers, done, have_layer_borders
    real, dimension(:), pointer :: level_values
    integer :: nlevs, iTmp, iStat
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: pLayerBorderSet
    real, dimension(:,:), pointer :: pBorderVals
    type(silam_sp) :: sp
    real :: fScale

    level_data => fu_work_string()
    layer_data => fu_work_string()
    level_values => fu_work_array()
    pBorderVals => fu_work_array_2d()
    sp%sp => fu_work_string()
    nullify(pLayerBorderSet)
    nLevs = 0
    !
    ! What do we have as desription of the levels?
    !
    level_data = fu_content(nl, 'levels')
    layer_data = fu_content(nl, 'layer_thickness')
    call get_items(nl, 'layer_borders', pLayerBorderSet, nLevs)
    if(nLevs > 0)then
      have_layer_borders = .true.
      do iTmp = 1, nLevs
        sp%sp = fu_content(pLayerBorderSet(iTmp))
        read(unit=sp%sp,fmt=*,iostat=iStat) pBorderVals(1,iTmp), pBorderVals(2,iTmp)
        if(fu_fails(iStat==0,'Failed to get layer borders from:'+fu_content(pLayerBorderSet(iTmp)), &
                           & 'create_levels_from_namelist_v2'))return
      end do
    else
      have_layer_borders = .false.
    endif
    chUnit = fu_content(nl,'vertical_unit')
    if(len_trim(chUnit) == 0)then
      fScale = 1.0
      chSIunit="???"
    else
      chSIunit = fu_SI_unit(chUnit)
      fScale = fu_factor_to_basic_unit(chUnit, chSIunit)
      if(error .or. fScale == real_missing)then
        call set_error('Failed to find factor to basic unit for:' + chUnit,'create_levels_from_namelist_v2')
        return
      endif
    endif
    
    ! See if we are given levels or layers
    
    vertical_method = fu_str_u_case(fu_content(nl, 'vertical_method'))
    level_type = fu_str_u_case(fu_content(nl, 'level_type'))

    need_layers = .false.
    need_levels = .false.
    done = .false.

    ! Check whether we have thickness or midpoints:
    !
    have_levels = level_data /= ''
    have_layers = layer_data /= ''

    ! Check whether we need thickness or midpoints:
    ! 
    select case (vertical_method)      ! We can also handle a few special cases here
    case ('CUSTOM_LEVELS')
      need_levels = .true.
    case ('CUSTOM_LAYERS')
      need_layers = .true.
    case ('SURFACE_LEVEL')
      call make_single_level(surface_level, levels) !, 1.0)
      done = .true.
    case ('TOP_ATMOSPHERE_LEVEL')
      call make_single_level(top_atmosphere_level, levels) !, 1.0)
      done =.true.
    case ('MEAN_SEA_LEVEL')
      call make_single_level(mean_sea_level, levels) !, 1.0)
      done = .true.
    case ('ENTIRE_ATMOSPHERE_MEAN_LAYER')
      call make_single_level(entire_atmosphere_mean_level, levels) !, 1.0)
      done = .true.
    case ('ENTIRE_ATMOSPHERE_LAYER')
      call make_single_level(entire_atmosphere_integr_level, levels) !, 1.0)
      done = .true.
    case ('ENTIRE_ATMOSPHERE_INTEGR_LAYER')
      call make_single_level(entire_atmosphere_integr_level, levels) !, 1.0)
      done = .true.
    case default
      call set_error('Strange vertical_method: '//trim(vertical_method), 'create_levels_from_namelist')
      call cleanup(have_layer_borders)
      return
    end select ! vertical_method
    
    if (done) then
      call cleanup(have_layer_borders)
      return
    end if
    !
    ! Simple cases have not materialised, have to do real work
    !
    if (.not. (need_levels .or. need_layers)) then
      call set_error('Something strange with creating levels', 'create_levels_from_namelist')
      call cleanup(have_layer_borders)
      return
    end if

    select case(level_type)
    case ('PRESSURE')
      if (.not. (have_layers .or. have_levels .or. have_layer_borders)) then
        call set_error('Pressure levels requested, but neither levels nor layer_thickness nor layer_borders given', &
                     & 'create_levels_from_namelist')
        call cleanup(have_layer_borders)
        return
      end if
      if(fu_fails(chSIunit=='Pa' .or. fScale==1.0,'For pressure level type, Pa must be vertical_unit','create_levels_from_namelist'))return
      if (have_levels) call split_string(level_data, ' ', level_values, nlevs)
      if (have_layers) call split_string(layer_data, ' ', level_values, nlevs)
      if (error) return

      if (need_levels) then
        if (.not. have_levels) call layers_to_levels(level_values, nlevs, constant_pressure, &
                                                   & have_layer_borders, pBorderVals)
        call make_levels(level_values, nlevs, levels, constant_pressure, fScale)
      else ! need_layers
        if(have_layer_borders)then
          call make_layers_from_borders(pBorderVals, nlevs, levels, layer_btw_2_altitude, fScale)
        else
          if (.not. have_layers) call levels_to_layers_pr(level_values, nlevs)
          call make_layers(level_values, nlevs, levels, layer_btw_2_pressure, fScale)
        endif
      end if
        
    case ('HEIGHT_FROM_SURFACE')
      if (.not. (have_layers .or. have_levels .or. have_layer_borders)) then
        call set_error('Height levels requested, but neither levels nor layer_thickness nor layer_borders given', &
                     & 'create_levels_from_namelist')
        call cleanup(have_layer_borders)
        return
      end if
      if(fu_fails(chSIunit=='m' .or. fScale==1.0,'For height level type, m must be vertical_unit','create_levels_from_namelist'))return
      if (have_levels) call split_string(level_data, ' ', level_values, nlevs)
      if (have_layers) call split_string(layer_data, ' ', level_values, nlevs)
      if (error) return

      if (need_levels) then
        if (.not. have_levels) call layers_to_levels(level_values, nlevs, constant_height, &
                                                   & have_layer_borders, pBorderVals)
        call make_levels(level_values, nlevs, levels, constant_height, fScale)
      else
        if(have_layer_borders)then
          call make_layers_from_borders(pBorderVals, nlevs, levels, layer_btw_2_altitude, fScale)
        else
          if (.not. have_layers) call levels_to_layers_hgt_alt(level_values, nlevs)
          call make_layers(level_values, nlevs, levels, layer_btw_2_height, fScale)
        endif
      end if

    case ('ALTITUDE_FROM_SEA')
      if (.not. (have_layers .or. have_levels)) then
        call set_error('Altitude levels requested, but neither levels nor layer_thickness nor layer_borders  given', &
                     & 'create_levels_from_namelist')
        call cleanup(have_layer_borders)
        return
      end if
      if(fu_fails(chSIunit=='m' .or. fScale==1.0,'For altitude level type, m must be vertical_unit','create_levels_from_namelist'))return
      if (have_levels) call split_string(level_data, ' ', level_values, nlevs)
      if (have_layers) call split_string(layer_data, ' ', level_values, nlevs)
      if (error) return

      if (need_levels) then
        if (.not. have_levels) call layers_to_levels(level_values, nlevs, constant_altitude, have_layer_borders, pBorderVals)
        call make_levels(level_values, nlevs, levels, constant_altitude, fScale)
      else
        if(have_layer_borders)then
          call make_layers_from_borders(pBorderVals, nlevs, levels, layer_btw_2_altitude, fScale)
        else
          if (.not. have_layers) call levels_to_layers_hgt_alt(level_values, nlevs)
          call make_layers(level_values, nlevs, levels, layer_btw_2_altitude, fScale)
        endif
      end if

    case ('HYBRID')
      !
      ! Note the different scheme:
!      nlevs = fu_content_int(nl, 'number_of_levels')
!      if (nlevs == int_missing) then
!        call set_error('Hybrid levels requires number_of_levels', 'create_levels_from_namelist')
!        call cleanup()
!        return
!      end if
      if(fu_fails(chSIunit=='Pa' .or. fScale==1.0,'For hybrid level type, Pa must be vertical_unit','create_levels_from_namelist'))return
      if (need_layers) then
        call make_layers_hybrid(nl, nlevs, levels, fScale)  ! nlevs and levels are output
      else
        call make_levels_hybrid(nl, nlevs, levels, fScale)
      end if
      
    case ('')
      call report(nl)
      call set_error('Missing level_type', 'create_levels_from_namelist')

    case default
      call set_error('Strange level_type: ' // trim(level_type), 'create_levels_from_namelist')
      
    end select ! level_type
    !
    ! Still have a chance to complain to the user:
    !
!    if (nlevs /= fu_content_int(nl, 'number_of_levels')) then
!      call msg_warning('number_of_levels does not match the level definition', 'create_levels_from_namelist')
!    end if

    call cleanup(have_layer_borders)
    
  contains

    !==========================================================================  

    subroutine cleanup(have_layer_borders)
      implicit none
      logical, intent(in) :: have_layer_borders
      call free_work_array(level_values)
      call free_work_array(layer_data)
      call free_work_array(level_data)
      call free_work_array(pBorderVals)
      call free_work_array(sp%sp)
      if(have_layer_borders)call destroy_items(pLayerBorderSet)
    end subroutine cleanup
    
    !========================================================================

    subroutine levels_to_layers_pr(values, nlevs)
      implicit none
      real, dimension(:), intent(inout) :: values
      integer, intent(in) :: nlevs

      real :: ftmp, gtmp
      logical :: ifOK
      integer :: i

      !
      ! Try rigorous method first. fTmp will keep the borders btw layers
      !
      fTmp = values(1)
      values(1) = std_pressure_sl ! Lower border
      ifOK = .true.
      do i=1,nLevs
        if(fTmp > values(i))then
          call msg_warning('Level centres are carelessly placed, use crude algorithm',&
                         & 'create_levels_from_namelist')
          ifOK = .false.
          exit
        endif
        gTmp = values(i+1)
        values(i+1)= 2.*fTmp - values(i)  !=values(i) - 2.*(values(i)-fTmp)
        fTmp = gTmp
      end do
      if(.not.ifOK)then
        !
        ! Have to use crude algorithm
        !
        fTmp = values(1)
        values(1) = std_pressure_sl
        do i=2, nLevs
          gTmp = values(i)
          values(i) = (fTmp + gTmp) * 0.5
          fTmp = gTmp
        enddo
        values(nLevs+1)= 2.*fTmp - values(nLevs)
      endif

    end subroutine levels_to_layers_pr
  
    !=======================================================================================
    
    subroutine levels_to_layers_hgt_alt(values, nlevs)
      implicit none
      real, dimension(:), intent(inout) :: values
      integer, intent(in) :: nlevs

      real :: ftmp, gtmp
      logical :: ifOK
      integer :: i

      !
      ! Non-pressure system. Try rigorous algorithm first
      !
      fTmp = values(1)
      values(1) = 0.
      ifOK = .true.
      do i=1, nLevs
        if(fTmp < values(i))then
          call msg_warning('Level centres are carelessly placed, use crude algorithm',&
               & 'create_levels_from_namelist')
          ifOK = .false.
          exit
        endif
        gTmp = values(i+1)
        values(i+1)= 2*fTmp - values(i)  !=values(i) + 2. * (fTmp - values(i))
        fTmp = gTmp
      end do
      if(.not.ifOK)then
        !
        ! Have to use crude algorithm
        !
        fTmp = values(1)
        values(1) = 0.
        do i=2, nLevs
          gTmp = values(i)
          values(i) = (fTmp + gTmp) * 0.5
          fTmp = gTmp
        enddo
        values(nLevs+1) = 2.* fTmp - values(nLevs)
      endif

    end subroutine levels_to_layers_hgt_alt

    !=======================================================================================

    subroutine layers_to_levels(values, nlevs, level_type_flag, ifLayerBorders, pBorders)
      implicit none
      real, dimension(:), intent(inout) :: values
      integer, intent(in) :: nlevs, level_type_flag
      logical, intent(in) :: ifLayerBorders
      real, dimension(:,:), pointer :: pBorders

      real :: ftmp, gtmp
      logical :: ifOK
      integer :: i

      if(ifLayerBorders)then   ! levels are just mid-points of layers
        do i = 1, nLevs
          values(i) = (pBorders(1,i) + pBorders(2,i)) / 2.
        enddo
        return
      endif
      
      if (level_type_flag == constant_pressure) then     ! here make them from thickness
        fTmp = std_pressure_sl - values(1)
        do i=1,nLevs
          values(i) = fTmp + values(i)*0.5 ! Upper border + half-thickness
          fTmp = fTmp - values(i+1)   ! values is big enough, last fTmp is void
        end do
      else
        fTmp=values(1)
        do i=1,nLevs
          values(i) = fTmp - values(i)*0.5
          fTmp = fTmp + values(i+1)
        end do
      end if

    end subroutine layers_to_levels

    !=======================================================================================

    subroutine make_single_level(level, levels)
      implicit none
      type(silja_level), intent(in) :: level
      type(silja_level), dimension(:), pointer :: levels
      
      integer :: stat

      allocate(levels(1), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', 'make_single_level')) return
      levels(1) = level
    end subroutine make_single_level

    !=======================================================================================

    subroutine make_levels(values, nlevs, levels, level_type_flag, fScale)
      implicit none
      real, dimension(:), intent(in) :: values
      integer, intent(in) :: nlevs, level_type_flag
      type(silja_level), dimension(:), pointer :: levels
      real, intent(in) :: fScale
      
      integer :: stat, ilev
      
      allocate(levels(nlevs), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', 'make_levels_hgt')) return
            
      do ilev = 1, nlevs
        levels(ilev) = fu_set_level_general(level_type_flag, fval1=values(ilev)*fScale)
      end do

    end subroutine make_levels
    
    !=======================================================================================

    subroutine make_layers(values, nlevs, levels, level_type_flag, fScale)
      implicit none
      real, dimension(:), intent(in) :: values
      integer, intent(in) :: nlevs, level_type_flag
      type(silja_level), dimension(:), pointer :: levels
      real, intent(in) :: fScale
      
      integer :: stat, ilev
      real :: top, bottom

      allocate(levels(nlevs), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', 'make_levels_hgt')) return
      
      if (level_type_flag == layer_btw_2_pressure) then
        bottom = std_pressure_sl
      else
        bottom = 0.0
      end if
      do ilev = 1, nlevs
        top = bottom + values(ilev)*fScale
        levels(ilev) = fu_set_level_general(level_type_flag, fval1=bottom, fval2=top)
        bottom = top
      end do

    end subroutine make_layers
    
    !=======================================================================
    
    subroutine make_layers_from_borders(pBorders, nlevs, levels, level_type_flag, fScale)
      implicit none
      integer, intent(in) :: nlevs, level_type_flag
      type(silja_level), dimension(:), pointer :: levels
      real, dimension(:,:), pointer :: pBorders
      real, intent(in) :: fScale

      integer :: stat, ilev

      allocate(levels(nlevs), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', 'make_levels_hgt')) return

      do ilev = 1, nlevs
        levels(ilev) = fu_set_level_general(level_type_flag, &
                                     & fval1=pBorders(1,iLev)*fScale, fval2=pBorders(2,iLev)*fScale)
      end do

    end subroutine make_layers_from_borders
    
    !=======================================================================================

    subroutine make_levels_hybrid(nlPtr, num_levels, levels, fScale)
      implicit none
!      type(Tsilam_namelist), intent(in), target :: nl
      type(silja_level), dimension(:), pointer :: levels
      integer, intent(out) :: num_Levels
      real, intent(in) :: fScale

      type(Tsilam_nl_item_ptr), dimension(:), pointer :: p_items
      type(Tsilam_namelist), pointer :: nlptr
      integer :: ind_item, lev_num, stat
      real :: a, b
      character(len=fnlen) :: content

      nullify(p_items)
      call get_items(nlptr, 'hybrid_coefficients', p_items, num_levels)
      if (fu_fails(num_levels > 0, 'Bad number of coefficients', 'set_hybrid_levels')) return

      allocate(levels(num_levels), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', 'make_levels_hybrid')) return

!      nlptr => nl
      do ind_item = 1, num_levels
        content = fu_content(p_items(ind_item))
        read(unit=content, fmt=*, iostat=stat) lev_num, a, b
        if (fu_fails(stat == 0, 'Failed to parse level:' // trim(content), 'set_hybrid_levels')) return
        levels(lev_num) = fu_set_hybrid_level(lev_num, a*fScale, b)
      end do

      deallocate(p_items)

    end subroutine make_levels_hybrid

    !=======================================================================================

    subroutine make_layers_hybrid(nlPtr, num_layers, layers, fScale)
      implicit none
      type(Tsilam_namelist), pointer :: nlPtr
      type(silja_level), dimension(:), pointer :: layers
      integer, intent(out) :: num_layers
      real, intent(in) :: fScale

      character(len=fnlen) :: content
      type(silja_level), dimension(:), allocatable :: tops, bottoms
      integer :: stat, lev_num, ilev, ind_item
      type(Tsilam_nl_item_ptr), dimension(:), pointer :: p_items
!      type(Tsilam_namelist), pointer :: nlptr
      real :: a, b

      ! Note the different numbering of the top and bottom levels!
      !
      nullify(p_items)
      call get_items(nlptr, 'hybrid_coefficients_bottom', p_items, num_layers)
      if (fu_fails(num_layers > 0, 'Bad number of bottom levels', 'set_hybrid_layers')) return

!      call get_items(nlptr, 'hybrid_coefficients_top', p_items, stat)
!      if (fu_fails(num_layers == stat + 1, 'Bad number of top levels', 'set_hybrid_layers')) return

      allocate(tops(num_layers), bottoms(num_layers), levels(num_layers), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', 'set_hybrid_layers')) return

      !
      ! Note: top of the level i is bottom of the level i+1
      !
      do ilev = 1, num_layers
        bottoms(ilev) = level_missing
      end do

      do ind_item = 1, num_layers
        content = fu_content(p_items(ind_item))
        read(unit=content, fmt=*, iostat=stat) lev_num, a, b
        if (fu_fails(stat == 0, 'Failed to parse level:' // trim(content), 'set_hybrid_layers')) return
        if(lev_num < 1 .or. lev_num > num_layers)then
          call set_error('Funny layer:'+content,'set_hybrid_layers')
          return
        endif
        bottoms(lev_num) = fu_set_hybrid_level(lev_num, a*fScale, b)
        if(lev_num > 1)tops(lev_num-1) = fu_set_hybrid_level(lev_num, a*fScale, b)
      end do
      
      content = fu_content(nlPtr,'hybrid_coefficients_domain_top')
      read(unit=content, fmt=*, iostat=stat) a, b
      if (fu_fails(stat == 0, 'Failed to parse domain_top:' + content, 'set_hybrid_layers')) return
      tops(num_layers) = fu_set_hybrid_level(num_layers, a*fScale, b)

!      do ind_item = 1, num_layers
!        content = fu_content(p_items(ind_item))
!        read(unit=content, fmt=*, iostat=stat) lev_num, a, b
!        if (fu_fails(stat == 0, 'Failed to parse level:' // trim(content), 'set_hybrid_layers')) return
!      end do
      
      do ilev = 1, num_layers
        if(defined(bottoms(iLev)))then
          layers(ilev) = fu_set_layer_between_two(tops(ilev), bottoms(ilev))
        else
          call set_error('Layer is missing:' + fu_str(ilev),'set_layers_hybrid')
          return
        endif
      end do
      
      deallocate(tops, bottoms)
      deallocate(p_items)

    end subroutine make_layers_hybrid

  end subroutine create_levels_from_namelist_v2


  !****************************************************************

  subroutine set_level_in_vertical(vertical, iLev, level)
    !
    ! Sets a new vertical structure copying its main parameters from 
    ! the given level
    !
    implicit none

    ! Imported values with intent IN
    type(silam_vertical), intent(inout) :: vertical
    type(silja_level), intent(in) :: level
    integer, intent(in) :: iLev

    if(.not.defined(vertical))then
      call set_error('Undefined vertical given','set_level_in_vertical')
      return
    endif
    if(iLev > vertical%nLevs)then
      call set_error('Too large index given','set_level_in_vertical')
      return
    endif
    if(.not.defined(level))then
      call set_error('Undefined level given','set_level_in_vertical')
      return
    end if

    vertical%levs(iLev) = level

  end subroutine set_level_in_vertical


  !******************************************************************

  subroutine arrange_levels_in_vertical(vert, ifChanged)
    !
    ! Arranges the vertical so that it satisfies the SILAM standard: the levels
    ! must start from the surface and go upwards one-by-one.
    ! Method used: a raising bulb
    !
    implicit none

    ! Imported parameter
    type(silam_vertical), intent(inout) :: vert
    logical, intent(out), optional :: ifChanged

    ! Local variables
    integer :: iLev, nChanges
    type(silja_level) :: levTmp
    logical :: ifAction
  
    if(.not.defined(vert))then
      call set_error('Undefined vertical given','arrange_levels_in_vertical')
      return
    endif
    
    if(vert%nLevs == 1)then
      if (present(ifChanged))then
        ifChanged = .false.
        return
      endif
    endif
    
    nChanges = 0
    ifAction = .true.
    do while(ifAction)
      ifAction = .false.
      do iLev = 1, vert%nLevs-1
        if(vert%levs(iLev) > vert%levs(iLev+1))then
          ifAction = .true.
          levTmp = vert%levs(iLev)
          vert%levs(iLev) = vert%levs(iLev+1)
          vert%levs(iLev+1) = levTmp
          nChanges = nChanges +1
        endif
      end do
    end do
    if (present(ifChanged))then
      if(nChanges > 0)then
        ifChanged = .true.
      else
        ifChanged = .false.
      endif
    endif

  end subroutine arrange_levels_in_vertical


  !*****************************************************************

  subroutine set_named_level_with_fract(nlItem, level, fFract, vertUnit)
    !
    ! Sets a level/layer that is defined as 
    ! vert_level = <type> <val1> <val2> <val3> <val4>
    ! The last value is ignored as it is by-default is the fraction os something
    ! that comes to this layer. Stupid agreement, of course but...
    !
    implicit none

    ! Imported parameters
    type(Tsilam_nl_item_ptr), intent(in) :: nlItem
    type(silja_level), intent(out) :: level
    real, intent(out), optional :: fFract
    character (len=*), optional, intent(in) :: vertUnit

    ! Local vars
    type(silam_sp) :: sp
    real, dimension(:), pointer :: fPtr
    real :: factor
    integer :: i, nVals, status

    factor = 1.
    sp%sp => fu_work_string()
    fPtr => fu_work_array()
    if(error)return


    sp%sp = fu_content(nlItem)
    if(present(fFract))then
      nVals = fu_nbrOfWords(sp%sp)-2  ! How many values => level or layer?
      read(unit = sp%sp(index(sp%sp,' '):),fmt=*, iostat=status) (fPtr(i),i=1,nVals),fFract
    else
      nVals = fu_nbrOfWords(sp%sp)-1  ! How many values => level or layer?
      read(unit = sp%sp(index(sp%sp,' '):),fmt=*, iostat=status) (fPtr(i),i=1,nVals)
    endif
    
    if(status /= 0)then
      call set_error('Failed to read:' + sp%sp,'set_named_level_with_fract')
      return
    endif
    
    !
    ! Set the level depending on the type and the number of values
    !
    if(index(sp%sp, 'PRESSURE') > 0)then
      if (present(vertUnit)) factor = fu_conversion_factor(vertUnit, "Pa")
      if(nVals == 1)then
        level = fu_set_pressure_level(fPtr(1)*factor)
      else
        level = fu_set_layer_between_two(layer_btw_2_pressure,fPtr(2)*factor,fPtr(1)*factor)
      endif

    elseif(index(sp%sp, 'HEIGHT_FROM_SURF') > 0)then
      if (present(vertUnit)) factor = fu_conversion_factor(vertUnit, "m")
      if(nVals == 1)then
        level = fu_set_constant_height_level(fPtr(1))
      else
        level = fu_set_layer_between_two(layer_btw_2_height,fPtr(2)*factor,fPtr(1)*factor)
      endif

    elseif(index(sp%sp, 'ALTITUDE_FROM_SEA_LEVEL') > 0)then
      if (present(vertUnit)) factor = fu_conversion_factor(vertUnit, "m")
      if(nVals == 1)then
        level = fu_set_constant_altitude_level(fPtr(1)*factor)
      else
        level = fu_set_layer_between_two(layer_btw_2_altitude,fPtr(2)*factor,fPtr(1)*factor)
      endif

    elseif(index(sp%sp, 'HYBRID') > 0)then
      if(nVals == 1)then
        level = fu_set_hybrid_level_no_coeff(int(fPtr(1)+0.5))
      else
        level = fu_set_layer_between_two(layer_btw_2_hybrid,fPtr(2),fPtr(1))
      endif

    elseif(index(sp%sp, 'SURFACE_LEVEL') > 0)then
      level = surface_level

    elseif(index(sp%sp, 'TOP_ATMOSPHERE_LEVEL') > 0)then
      level = top_atmosphere_level

    elseif(index(sp%sp, 'MEAN_SEA_LEVEL') > 0)then
      level = mean_sea_level

    elseif(index(sp%sp, 'ENTIRE_ATMOSPHERE_LAYER') > 0)then
      level = entire_atmosphere_integr_level

    elseif(index(sp%sp, 'ENTIRE_ATMOSPHERE_MEAN_LAYER') > 0)then
      level = entire_atmosphere_mean_level

    elseif(index(sp%sp, 'ENTIRE_ATMOSPHERE_INTEGR_LAYER') > 0)then
      level = entire_atmosphere_integr_level

    else
      call set_error(fu_connect_strings('Strange type of level:',sp%sp), &
                   & 'fu_set_named_level')
      return
    endif

    call free_work_array(sp%sp)
    call free_work_array(fPtr)

  end subroutine set_named_level_with_fract


  !*********************************************************************

  subroutine write_vert_and_fract_as_nmLst(uOut, vertical, fractions, ifIgnoreZeroes)
    !
    ! Stores the vertical levels/layers as a namelist in the style compatible with the
    ! source term definitions
    !
    implicit none

    ! Imported parameters
    integer, intent(in) :: uOut    ! file unit to write into
    type(silam_vertical), intent(in) :: vertical  ! vertical to report
    real, dimension(:), intent(in) :: fractions      ! fractions to be reported
    logical, intent(in) :: ifIgnoreZeroes         ! if skip levels with zero fractions

    !Local variables
    integer :: iLev

    ! Stupidity check
    !
    if(.not. (vertical%defined == silja_true))then
      call set_error('Cannot report undefined vertical','write_vert_and_fract_as_nmLst')
      return
    endif

    select case(vertical%levs(1)%levelType)
      case(constant_pressure, sigma_level, hybrid, layer_btw_2_pressure)
        write(uOut,fmt='(A)')'vertical_unit = Pa'

      case(constant_altitude, constant_height, depth_level, layer_btw_2_altitude, layer_btw_2_height)
        write(uOut,fmt='(A)')'vertical_unit = m'
      case default
        call report (vertical)
        call msg ('level type',vertical%levs(1)%levelType)
        call set_error('','')
    end select

    do iLev = 1, vertical%nLevs
      if(ifIgnoreZeroes)then
        if(fractions(iLev) < 1e-5)cycle
      endif
      call write_level_and_fract_as_nmLst(uOut, vertical%levs(iLev), fractions(iLev))
    end do

  end subroutine write_vert_and_fract_as_nmLst


  !************************************************************************

  subroutine write_level_and_fract_as_nmLst(uOut, level, fraction)
    !
    ! Writes level and its fraction into the given file in the namelist format
    !
    implicit none

    ! Imported parameters
    integer, intent(in) :: uOut
    type(silja_level), intent(in) :: level
    real, intent(in) :: fraction

    ! Local variable
    type(silam_sp) :: sp

    sp%sp => fu_work_string()
    if(error)return
    call level_to_short_string(level, sp%sp)
    write(uOut,fmt='(2A,2x,F15.7)') 'vert_level = ', trim(sp%sp), fraction
    call free_work_array(sp%sp)

  end subroutine write_level_and_fract_as_nmLst


  !************************************************************************

  subroutine level_to_short_string(level, strOut)
    !
    ! Writes level and its fraction into the given file in the namelist format
    !
    implicit none

    ! Imported parameters
    type(silja_level), intent(in) :: level
    character(len=*), intent(out) :: strOut

    !
    ! Writing depends on the level type
    !
    select case(level%levelType)
      case(constant_pressure)
        write(unit=strOut,fmt='(A,2x,F15.7)')'PRESSURE', level%a
      case(layer_btw_2_pressure)
        write(unit=strOut,fmt='(A,2x,2(F15.7,2x))')'PRESSURE', level%a, level%b

      case(constant_altitude)
        write(unit=strOut,fmt='(A,2x,F15.7)')'ALTITUDE_FROM_SEA_LEVEL', level%a
      case(layer_btw_2_altitude)
        write(unit=strOut,fmt='(A,2x,2(F15.7,2x))')'ALTITUDE_FROM_SEA_LEVEL', level%b, level%a

      case(constant_height)
        write(unit=strOut,fmt='(A,2x,F15.7)')'HEIGHT_FROM_SURF', level%a
      case(layer_btw_2_height)
        write(unit=strOut,fmt='(A,2x,2(F15.7,2x))')'HEIGHT_FROM_SURF', level%b, level%a

      case(sigma_level)
        write(unit=strOut,fmt='(A,2x,F15.7)')'SIGMA', level%a
      case(layer_btw_2_sigma)
        write(unit=strOut,fmt='(A,2x,2(F15.7,2x))')'SIGMA', level%b, level%a

      case(hybrid)
        write(unit=strOut,fmt='(A,2x,I3,2(2x,F15.7))')'HYBRID', level%number, level%a, level%b
      case(layer_btw_2_hybrid)
        write(unit=strOut,fmt='(A,2x,2(I3,2(2x,F15.7)))')'HYBRID', &
            & level%number2, level%a2, level%b2, level%number, level%a, level%b

      case(depth_level)
        write(unit=strOut,fmt='(A,2x,F15.7)')'DEPTH', level%a
      case(layer_btw_2_depth)
        write(unit=strOut,fmt='(A,2x,2(F15.7,2x))')'DEPTH', level%b, level%a

      case(entire_atmosphere_single_layer)
        write(unit=strOut,fmt='(A)')'ENTIRE_ATMOSPHERE_LAYER'

      case(surface)
        write(unit=strOut,fmt='(A)')'SURFACE_LEVEL'

      case(top_of_the_atmosphere)
        write(unit=strOut,fmt='(A)')'TOP_ATMOSPHERE_LEVEL'

      case(mean_sea)
        write(unit=strOut,fmt='(A)')'MEAN_SEA_LEVEL', level%a

      case default
        if(defined(level))then
          call msg_warning('Unsuppoted level','write_level_and_fract_as_nmLst')
          call report(level)
          call set_error('Unsuppoted level','write_level_and_fract_as_nmLst')
        else
          write(unit=strOut,fmt='(A)')'undefined'
        endif
    end select

  end subroutine level_to_short_string


  !*****************************************************************

  function fu_central_level_of_layer(layerIn) result(levelOut)
    !
    ! Creates the level placed in the central point of the input layer.
    !
    implicit none

    ! return value of the function
    type(silja_level) :: levelOut

    ! Imported parameter 
    type(silja_level), intent(in) :: layerIn

    select case (layerIn%leveltype)
      case(constant_pressure, &
         & constant_height, &
         & hybrid, &
         & surface, &
         & top_of_the_atmosphere, &
         & mean_sea, &
         & constant_altitude, &
         & sigma_level, &
         & depth_level, &
         & entire_atmosphere_single_layer)

        levelOut = layerIn

      case(layer_btw_2_pressure)
        levelOut%leveltype = constant_pressure
        levelOut%hybrid_coeff_known = .false.
        levelOut%number = int_missing
        levelOut%number2 = int_missing
        levelOut%a = (layerIn%a + layerIn%b) * 0.5
        levelOut%b = real_missing
        levelOut%defined = silja_true

      case(layer_btw_2_altitude)
        levelOut%leveltype = constant_altitude
        levelOut%hybrid_coeff_known = .false.
        levelOut%number = int_missing
        levelOut%number2 = int_missing
        levelOut%a = (layerIn%a + layerIn%b) * 0.5
        levelOut%b = real_missing
        levelOut%defined = silja_true

      case(layer_btw_2_height)
        levelOut%leveltype = constant_height
        levelOut%hybrid_coeff_known = .false.
        levelOut%number = int_missing
        levelOut%number2 = int_missing
        levelOut%a = (layerIn%a + layerIn%b) * 0.5
        levelOut%b = real_missing
        levelOut%defined = silja_true

      case(layer_btw_2_sigma)
        levelOut%leveltype = sigma_level
        levelOut%hybrid_coeff_known = .false.
        levelOut%number = int_missing
        levelOut%number2 = int_missing
        levelOut%a = (layerIn%a + layerIn%b) * 0.5
        levelOut%b = real_missing
        levelOut%defined = silja_true

      case(layer_btw_2_hybrid)
        levelOut%leveltype = hybrid
        levelOut%hybrid_coeff_known = layerIn%hybrid_coeff_known
        !levelOut%number = int((layerIn%number + layerIn%number2)*0.5 +0.5)
        levelOut%number = min(layerIn%number,layerIn%number2) ! Maps layers to 1-N levels 
        levelOut%number2 = int_missing
        levelOut%a = (layerIn%a + layerIn%a2) * 0.5
        levelOut%b = (layerIn%b + layerIn%b2) * 0.5
        levelOut%defined = silja_true

      case(layer_btw_2_depth)
        levelOut%leveltype = depth_level
        levelOut%hybrid_coeff_known = .false.
        levelOut%number = int_missing
        levelOut%number2 = int_missing
        levelOut%a = (layerIn%a + layerIn%b) * 0.5
        levelOut%b = real_missing
        levelOut%defined = silja_true

      case default
        call report(layerIn)
        call set_error('Unknown level type','fu_central_level_of_layer')
    end select

  end function fu_central_level_of_layer


  !*****************************************************************

  subroutine add_level(vertical, level)
    !
    ! Adds one more level to the existing vertical structure
    ! No level sorting is done.
    !
    implicit none

    ! Imported parameters
    type(silja_level), intent (in) :: level
    type(silam_vertical), intent (inout) :: vertical

    ! Local variables
    integer :: i
    type(silja_level), dimension(:), pointer :: levsTmp

    if(.not.defined(vertical))then
      call set_error('Undefined vertical given','add_level')
      return
    endif
    if(.not.defined(level))then
      call set_error('Undefined level is given','add_level')
      return
    end if
    !
    ! Vertical structure and the level must be compatible
    !
    if(vertical%vert_type /= fu_leveltype(level))then
      call set_error('Level type is not the same in vertical and level','add_level')
      return
    end if
    !
    ! Now we have to reallocate the space: a painful price for pointers
    !
    if(vertical%Nlevs > 0)then
      allocate(levsTmp(vertical%Nlevs), stat=i)
      if(i/= 0)then
        call set_error('Failed to allocate temporary levels','add_level')
        return
      endif
      levsTmp(1:vertical%Nlevs) = vertical%levs(1:vertical%Nlevs)
      deallocate(vertical%levs,stat=i)
      allocate(vertical%levs(vertical%Nlevs+1),stat=i)
      if(i/=0)then
        call set_error('Failed to reallocate memory','add_level')
        return
      endif
      vertical%levs(1:vertical%Nlevs) = levsTmp(1:vertical%Nlevs)
      deallocate(levsTmp, stat=i)
    else
      allocate(vertical%levs(vertical%Nlevs+1),stat=i)
      if(i/=0)then
        call set_error('Failed to reallocate memory','add_level')
        return
      endif
    endif
    !
    ! Now the last level can be added
    !
    vertical%Nlevs = vertical%Nlevs + 1
    vertical%levs(vertical%Nlevs) = level ! No sorting.

  end subroutine add_level


  ! ***************************************************************


  LOGICAL FUNCTION fu_level_defined(level)
    !
    ! Description:
    ! Returns a true value, if the level has been given a value using
    ! one of the setting functions.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silja_level), INTENT(in) :: level

    fu_level_defined = fu_true(level%defined)

  END FUNCTION fu_level_defined


  ! ***************************************************************


  LOGICAL FUNCTION fu_vertical_defined(vert)
    !
    ! Description:
    ! Returns a true value, if the vertical has been given a value using
    ! one of the setting functions.
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silam_vertical), INTENT(in) :: vert

    fu_vertical_defined = vert%defined == silja_true

  END FUNCTION fu_vertical_defined



  ! ***************************************************************


  LOGICAL FUNCTION fu_leveltypes_eq(lev1, lev2) result(eq)

    ! Description:
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    !
    ! Imported parameters with intent(in):
    TYPE(silja_level), INTENT(in) :: lev1, lev2

    ! 1. Check definitions.

    IF (.NOT.fu_level_defined(lev1)) THEN

      IF (.NOT.fu_level_defined(lev2)) THEN
        eq = .true. ! both undefined, so they're equal
        RETURN
      ELSE
        eq = .false. ! one defined, on not, so not equal
        RETURN
      END IF

    ELSE

      IF (.NOT.fu_level_defined(lev2)) THEN
        eq = .false. ! one defined, one not, so not equal
        RETURN

      ELSE

	! 2. Both defined.
	IF (lev1%leveltype == lev2%leveltype) THEN
	  eq = .true.
	ELSE
	  eq = .false.
	END IF
      END IF
    END IF

  END FUNCTION fu_leveltypes_eq


  ! ***************************************************************


  function fu_add_levels(lev1, lev2) result(level)
    !
    ! Overloading of the operator +
    !
    implicit none

    ! Return value of the function
    type(silja_level) :: level
    
    ! Imported parameters with intent IN
    type(silja_level), intent(in) :: lev1, lev2

    level = level_missing

    ! Stupidity check
    !
    if(.not.defined(lev1) .or. .not.defined(lev2))then
      call set_error('Undefiend level given','fu_add_levels')
      return
    endif

    ! Only same type levels can be summed
    !
    if(lev1%leveltype /= lev2%leveltype)then
      call set_error('Levels have different types','fu_add_levels')
      return
    endif

    ! Finally, sum-up the levels
    !
    select case(lev1%leveltype)
      case(surface,top_of_the_atmosphere,mean_sea, &
         & entire_atmosphere_single_layer,no_level)
        call set_error('2D levels can not be summed','fu_add_levels')
        return
      case(hybrid, layer_btw_2_hybrid)
        call set_error('Hybrid levels can not be summed','fu_add_levels')
        return
      case(constant_altitude,layer_btw_2_altitude, &
         & constant_height,layer_btw_2_height, & 
         & sigma_level,layer_btw_2_sigma, &
         & depth_level,layer_btw_2_depth, &
         & constant_pressure,layer_btw_2_pressure)

        level%leveltype = lev1%leveltype
        if((lev1%a .eps. real_missing) .or. (lev2%a .eps. real_missing)) then
          level%a = real_missing 
        else
          level%a = lev1%a + lev2%a
        endif
        if((lev1%b .eps. real_missing) .or. (lev2%b .eps. real_missing)) then
          level%b = real_missing 
        else
          level%b = lev1%b + lev2%b
        endif
        level%defined = silja_true

      case default
        call set_error('Unknown level type','fu_add_levels')
        return
    end select

  end function fu_add_levels


  ! ***************************************************************


  function fu_multiply_level(level_in, factor) result(level_out)
    !
    ! Overloading of the operator * as multiplication with the number
    !
    implicit none

    ! Return value of the function
    type(silja_level) :: level_out
    
    ! Imported parameters with intent IN
    type(silja_level), intent(in) :: level_in
    real, intent(in) :: factor

    level_out = level_missing

    ! Stupidity check
    !
    if(.not.defined(level_in) .or. (factor .eps. real_missing))then
      call set_error('Undefiend level or factor given','fu_multiply_level')
      return
    endif

    ! Multiply the level
    !
    select case(level_in%leveltype)
      case(surface,top_of_the_atmosphere,mean_sea, &
         & entire_atmosphere_single_layer,no_level)
        call set_error('2D level can not be multiplied with number','fu_multiply_level')
        return
      case(hybrid, layer_btw_2_hybrid)
        call set_error('Hybrid level can not be multiplied with number','fu_multiply_level')
        return
      case(constant_altitude,layer_btw_2_altitude, &
         & constant_height,layer_btw_2_height, & 
         & sigma_level,layer_btw_2_sigma, &
         & depth_level,layer_btw_2_depth, &
         & constant_pressure,layer_btw_2_pressure)

        level_out%leveltype = level_in%leveltype
        if(level_in%a .eps. real_missing) then
          level_out%a = real_missing 
        else
          level_out%a = level_in%a * factor
        endif
        if(level_in%b .eps. real_missing) then
          level_out%b = real_missing 
        else
          level_out%b = level_in%b * factor
        endif
        level_out%defined = silja_true

      case default
        call set_error('Unknown level type','fu_multiply_level')
        return
    end select

  end function fu_multiply_level


  ! ***************************************************************


  INTEGER FUNCTION fu_leveltype_of_level(level)

    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_level), INTENT(in) :: level
    !

    IF (fu_level_defined(level)) THEN
      fu_leveltype_of_level = level%leveltype
    ELSE
      fu_leveltype_of_level = no_level
    END IF

  END FUNCTION fu_leveltype_of_level


  ! ***************************************************************


  INTEGER FUNCTION fu_leveltype_of_vertical(vertical)

    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silam_vertical), INTENT(in) :: vertical
    !

    IF (fu_vertical_defined(vertical)) THEN
      fu_leveltype_of_vertical = vertical%vert_type
    ELSE
      fu_leveltype_of_vertical = no_level
    END IF

  END FUNCTION fu_leveltype_of_vertical



  ! ***************************************************************


  REAL FUNCTION fu_pr_level_pressure(level)

    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_level), INTENT(in) :: level

    IF (level%leveltype == constant_pressure) THEN
      fu_pr_level_pressure = level%a
    ELSE
      call msg('Strange level type and/or level:',level%leveltype, level%a)
      CALL set_error('this is not a pressure level','fu_pr_level_pressure')
    END IF

  END FUNCTION fu_pr_level_pressure



  ! ***************************************************************


  INTEGER FUNCTION fu_hybrid_level_number(level)

    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_level), INTENT(in) :: level

    IF (level%leveltype == hybrid) then
      fu_hybrid_level_number = level%number
    elseif (level%leveltype == layer_btw_2_hybrid) THEN
      fu_hybrid_level_number = level%number2
    ELSE
      CALL report(level)
      CALL set_error('this is not a defined hybird level'&
                   & ,'fu_hybrid_level_number')
    END IF

  END FUNCTION fu_hybrid_level_number


  ! ***************************************************************


  REAL FUNCTION fu_hybrid_level_coeff_a(level)

    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_level), INTENT(in) :: level

    fu_hybrid_level_coeff_a = real_missing

    IF (level%hybrid_coeff_known) THEN
      IF (level%leveltype == hybrid) THEN
        fu_hybrid_level_coeff_a = level%a
      ELSEIF (level%leveltype == layer_btw_2_hybrid) THEN
        fu_hybrid_level_coeff_a = 0.5*(level%a+level%a2)
      endif
    endif

    if (fu_hybrid_level_coeff_a == real_missing) then
      CALL report(level)
      CALL set_error('Failed' ,'fu_hybrid_level_coeff_a')
    END IF

  END FUNCTION fu_hybrid_level_coeff_a


  ! ***************************************************************


  REAL FUNCTION fu_hybrid_level_coeff_b(level)

    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_level), INTENT(in) :: level

    fu_hybrid_level_coeff_b = real_missing

    IF (level%hybrid_coeff_known) THEN
      IF (level%leveltype == hybrid) THEN
        fu_hybrid_level_coeff_b = level%b
      ELSEIF (level%leveltype == layer_btw_2_hybrid) THEN
        fu_hybrid_level_coeff_b = 0.5*(level%b+level%b2)
      endif
    endif

    if (fu_hybrid_level_coeff_b == real_missing) then
      CALL report(level)
      CALL set_error('Failed' ,'fu_hybrid_level_coeff_b')
    END IF

  END FUNCTION fu_hybrid_level_coeff_b

  !************************************************************************************
  
  subroutine hybrid_coefs(vertical, a_full, b_full, a_half, b_half)
    !
    ! Return the half or full-level hybrid coefficients in
    ! arrays. Half-level values are defined only for thick layers.
    implicit none
    type(silam_vertical), intent(in) :: vertical
    real, dimension(:), intent(out), optional :: a_full, b_full
    real, dimension(0:), intent(out), optional :: a_half, b_half
    
    integer :: nz, i

    nz = fu_nbrOfLevels(vertical)
        
    select case(fu_leveltype(vertical))
    case (hybrid)
      if (present(a_half) .or. present(b_half)) then
        call set_error('Requesting half-level parameters for thin layers', 'hybrid_coefs')
        return
      end if
      do i = 1, nz
        if (present(a_full)) then
          a_full(i) = fu_hybrid_level_coeff_a(fu_level(vertical, i))
        end if
        if (present(b_full)) then
          b_full(i) = fu_hybrid_level_coeff_b(fu_level(vertical, i))
        end if
      end do

    case (layer_btw_2_hybrid)
      do i = 1, nz
        if (present(a_full)) a_full(i) = fu_hybrid_level_coeff_a(fu_level(vertical, i, .true.))
        if (present(b_full)) b_full(i) = fu_hybrid_level_coeff_b(fu_level(vertical, i, .true.))
      end do

      if (present(a_half)) then
        a_half(0) = fu_hybrid_level_coeff_a(fu_lower_boundary_of_layer(fu_level(vertical, 1)))
        do i = 1, nz
          a_half(i) = fu_hybrid_level_coeff_a(fu_upper_boundary_of_layer(fu_level(vertical, i)))
        end do
      end if
      if (present(b_half)) then
        b_half(0) = fu_hybrid_level_coeff_b(fu_lower_boundary_of_layer(fu_level(vertical, 1)))
        do i = 1, nz
          b_half(i) = fu_hybrid_level_coeff_b(fu_upper_boundary_of_layer(fu_level(vertical, i)))
        end do
      end if

    case (sigma_level)
        if (present(a_full)) then
         do i = 1, nz
          a_full(i) = 0.
         enddo
        end if
        if (present(b_full)) then
         do i = 1, nz
          b_full(i) = vertical%levs(i)%a
         enddo
        end if


    case default
      call msg("Trouble getting hybrid_coefs from vertical")
      call report(vertical)
      do i = 1, nz
        call msg("Level", i)
        call report(vertical%levs(i))
          
      end do


      call set_error('Invalid leveltype for hybrid_coefs', 'hybrid_coefs')
    end select

  end subroutine hybrid_coefs


  ! ***************************************************************


  REAL FUNCTION fu_hybrid_level_pressure(level, ground_surface_pressure)

    ! Calculates the pressure on hybrid or sigma (model) level, when ground
    ! surface pressure is given.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_level), INTENT(in) :: level
    REAL, INTENT(in) :: ground_surface_pressure

    IF ((level%leveltype == hybrid).and.level%hybrid_coeff_known) THEN
      fu_hybrid_level_pressure = level%a + ground_surface_pressure * level%b
    else if (level%leveltype == sigma_level) then
      fu_hybrid_level_pressure = ground_surface_pressure * level%b
    else if (level%leveltype == layer_btw_2_hybrid .and. level%hybrid_coeff_known) then
      fu_hybrid_level_pressure = 0.5 * (level%a+level%a2) + 0.5 * (level%b + level%b2) * ground_surface_pressure
    ELSE
      CALL report(level)
      CALL set_error('cannot calculate hybrid level pressure','fu_hybrid_level_pressure')
    END IF

  END FUNCTION fu_hybrid_level_pressure



  ! ***************************************************************


  REAL FUNCTION fu_level_height(level)
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_level), INTENT(in) :: level

    IF (level%leveltype == constant_height .or. level%leveltype == constant_altitude ) THEN
      fu_level_height = level%a
    ELSE
      CALL report(level)
      fu_level_height = real_missing
      CALL set_error('this is not a constant height/altitude level','fu_level_height')
    END IF

  END FUNCTION fu_level_height



  ! ***************************************************************


  REAL FUNCTION fu_level_altitude(level)

    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_level), INTENT(in) :: level

    IF (level%leveltype == constant_altitude) THEN
      fu_level_altitude = level%a
    ELSE
      CALL report(level)
      CALL set_error('this is not a constant altitude level'&
	  & ,'fu_level_altitude')
    END IF

  END FUNCTION fu_level_altitude


  !*****************************************************************

  real function fu_sigma_level_sigma(level)
    !
    implicit none
    type(silja_level), intent(in) :: level

    IF (level%leveltype == sigma_level) THEN
      fu_sigma_level_sigma = level%a
    ELSE
      CALL report(level)
      CALL set_error('this is not a sigma level','fu_sigma_level_sigma')
    END IF

  end function fu_sigma_level_sigma


  !*******************************************************************

  real function fu_depth_level_depth(level)
    
    implicit none

    TYPE(silja_level), INTENT(in) :: level

    IF (level%leveltype == depth_level) THEN
      fu_depth_level_depth = level%a
    ELSE
      CALL report(level)
      CALL set_error('this is not a depth level','fu_level_altitude')
    END IF

  end function fu_depth_level_depth

 
  !*****************************************************************

  real function fu_level_value(level)
    !
    ! A universal function that returns the VALUE of the level, whatever it is.
    ! No unit conversion, no thick layers accepted, hybrid number is converted to real
    !
    implicit none 

    TYPE(silja_level), INTENT(in) :: level

    select case (level%leveltype)
      case(constant_pressure, &
         & constant_height, &
         & constant_altitude, &
         & sigma_level, &
         & depth_level)

        fu_level_value = level%a

      case(hybrid)

        fu_level_value = level%number

      case default
        call report(level)
        call set_error('Non-supported level type','fu_level_value')
    end select

  end function fu_level_value

  !*****************************************************************

  logical function fu_if_layer(level)
    !
    ! Checks whether the level represents a thick layer or it is a single level
    !
    implicit none

    ! Imported parameter 
    type(silja_level), intent(in) :: level

    select case (level%leveltype)
      case(constant_pressure, &
         & constant_height, &
         & hybrid, &
         & surface, &
         & top_of_the_atmosphere, &
         & mean_sea, &
         & constant_altitude, &
         & sigma_level, &
         & depth_level, &
         & entire_atmosphere_single_layer)

        fu_if_layer = .false.

      case(layer_btw_2_pressure, &
         & layer_btw_2_altitude, &
         & layer_btw_2_height, &
         & layer_btw_2_sigma, &
         & layer_btw_2_hybrid, &
         & layer_btw_2_depth)

        fu_if_layer = .true.

      case default
        call report(level)
        call set_error('Unknown level type','fu_if_layer')
    end select

  end function fu_if_layer


  !******************************************************************

  real function fu_layer_thickness_m(layer)
    !
    ! Computes the thickness of the given layer and translates it to metres.
    ! For the layers presented in the pressure or hybrid system, a 
    ! crude approximation is used.
    !
    implicit none

    ! Imported parameter 
    type(silja_level), intent(in) :: layer

    select case (layer%leveltype)

      case(layer_btw_2_pressure)
        fu_layer_thickness_m = fu_height_diff_between_pre(max(layer%a,layer%b), &
                                                        & min(layer%a,layer%b))
      case(layer_btw_2_altitude, layer_btw_2_height)
        fu_layer_thickness_m = abs(layer%a - layer%b)

      case(layer_btw_2_sigma)
        fu_layer_thickness_m = fu_height_diff_between_pre(std_pressure_sl* max(layer%a,layer%b), &
                                                        & std_pressure_sl* min(layer%a,layer%b))
      case(layer_btw_2_hybrid)
        fu_layer_thickness_m = fu_height_diff_between_pre( &
                     & fu_hybrid_level_pressure(fu_lower_boundary_of_layer(layer), std_pressure_sl), &
                     & fu_hybrid_level_pressure(fu_upper_boundary_of_layer(layer), std_pressure_sl))
      case(layer_btw_2_depth)

      case(constant_pressure, &
         & constant_height, &
         & hybrid, &
         & surface, &
         & top_of_the_atmosphere, &
         & mean_sea, &
         & constant_altitude, &
         & sigma_level, &
         & depth_level, &
         & entire_atmosphere_single_layer)
        call report(layer)
        call set_error('Not a thick layer','fu_layer_thickness_m')

      case default
        call report(layer)
        call set_error('Unknown level type','fu_layer_thickness_m')
    end select

  end function fu_layer_thickness_m


!  !************************************************************************************
!
!  subroutine vertical_layer_thicknesses_m(vert, thicknesses, nLevs)
!    !
!    ! Computes the thickness of the requested layer in the vertical. Importantly, 
!    ! does NOT use the above function  but also allows for levels in the vertical
!    ! In this case they are treated as middle-points
!    !
!    implicit none
!
!    ! Imported parameters    
!    type(silam_vertical), intent(in) :: vert
!    real, dimension(:), pointer :: thicknesses
!    integer, intent(out) :: nLevs
!
!    ! Local variables
!    real :: border
!    integer :: iLev
!
!    if(associated(thicknesses))then
!      if(size(thicknesses) < vert%nLevs)then
!        call msg('Too smal size of the thickness array',size(thicknesses))
!        call set_error('Too smal size of the thickness array','vertical_layer_thicknesses_m')
!        return
!      endif
!    else
!      allocate(thicknesses(vert%nLevs), stat = iLev)
!      if(iLev /= 0)then
!        call set_error('Failed to allocate thickness array','vertical_layer_thicknesses_m')
!        return
!      endif
!    endif
!    
!    nLevs = vert%nLevs
!
!    select case (vert%vert_type)
!
!      case(layer_btw_2_pressure)
!        do iLev = 1, vert%nLevs
!          thicknesses(iLev) = fu_height_diff_between_pre(max(vert%levs(iLev)%a, vert%levs(iLev)%b), &
!                                                        & min(vert%levs(iLev)%a, vert%levs(iLev)%b), &
!                                                        & temperature_ref)
!        end do
!
!      case(layer_btw_2_altitude, layer_btw_2_height, layer_btw_2_depth)
!        do iLev = 1, vert%nLevs
!          thicknesses(iLev) = abs(vert%levs(iLev)%a - vert%levs(iLev)%b)
!        end do
!
!      case(layer_btw_2_sigma)
!        do iLev = 1, vert%nLevs
!          thicknesses(iLev) = fu_height_diff_between_pre( &
!                                     & std_pressure_sl* max(vert%levs(iLev)%a,vert%levs(iLev)%b), &
!                                     & std_pressure_sl* min(vert%levs(iLev)%a,vert%levs(iLev)%b), &
!                                     & temperature_ref)
!        end do
!
!      case(layer_btw_2_hybrid)
!        do iLev = 1, vert%nLevs
!          thicknesses(iLev) = fu_height_diff_between_pre( &
!                     & fu_hybrid_level_pressure(fu_lower_boundary_of_layer(vert%levs(iLev)), std_pressure_sl), &
!                     & fu_hybrid_level_pressure(fu_upper_boundary_of_layer(vert%levs(iLev)), std_pressure_sl), &
!                     & temperature_ref)
!        end do
!
!      case(constant_pressure, &
!         & constant_height, &
!         & hybrid, &
!         & surface, &
!         & top_of_the_atmosphere, &
!         & mean_sea, &
!         & constant_altitude, &
!         & sigma_level, &
!         & depth_level, &
!         & entire_atmosphere_single_layer)
!         
!          thicknesses(1) = fu_level_height(vert%levs(1)) * 2.0
!          border = 0.
!          do iLev = 2, vert%nLevs
!            border = 2. * (fu_level_height(vert%levs(iLev))) - border
!            thicknesses(iLev) = (fu_level_height(vert%levs(iLev)) - border) * 2
!          end do
!
!      case default
!        call report(vert%levs(iLev))
!        call set_error('Unknown level type','vertical_layer_thicknesses_m')
!    end select
!
!  end subroutine vertical_layer_thicknesses_m


  !******************************************************************

  real function fu_layer_thickness_local_unit(layer)
    !
    ! Computes the thickness of the given layer and translates it to metres.
    ! For the layers presented in the pressure or hybrid system, a 
    ! crude approximation is used.
    !
    implicit none

    ! Imported parameter 
    type(silja_level), intent(in) :: layer

    select case (layer%leveltype)

      case(layer_btw_2_pressure)
        fu_layer_thickness_local_unit = abs(layer%b - layer%a) ! bottom - top pressure

      case(layer_btw_2_altitude, layer_btw_2_height)
        fu_layer_thickness_local_unit = abs(layer%a - layer%b) 

      case(layer_btw_2_sigma)
        fu_layer_thickness_local_unit = std_pressure_sl* abs(layer%a - layer%b)

      case(layer_btw_2_hybrid)
        fu_layer_thickness_local_unit = 1.0
        !fu_layer_thickness_local_unit = &
        !    & fu_hybrid_level_pressure(fu_lower_boundary_of_layer(layer), std_pressure_sl) - &
        !    & fu_hybrid_level_pressure(fu_upper_boundary_of_layer(layer), std_pressure_sl)

      case(layer_btw_2_depth)
        fu_layer_thickness_local_unit = abs(layer%a - layer%b)

      case(constant_pressure, &
         & constant_height, &
         & hybrid, &
         & surface, &
         & top_of_the_atmosphere, &
         & mean_sea, &
         & constant_altitude, &
         & sigma_level, &
         & depth_level, &
         & entire_atmosphere_single_layer)
        call report(layer)
        call set_error('Not a thick layer','fu_layer_thickness_local_unit')

      case default
        call report(layer)
        call set_error('Unknown level type','fu_layer_thickness_local_unit')
    end select

  end function fu_layer_thickness_local_unit

  !*****************************************************************

  real function fu_vert_overlap_fraction_old(layerThin, layerThick) result(fraction)
    !
    ! Returns a fraction of the thin layer, which is covered by the thick one
    ! Algorithm: if the layers are of the same type - just compute the overlap
    ! index, if they are of different types, first project the thin layer
    ! to the system of the thick one.
    !
    implicit none

    ! Imported parameters 
    type(silja_level), intent(in) :: layerThin, layerThick

    ! Local variables
    type(silja_level) :: layerTmp
    type(silam_vertical) :: vertTmp
    real :: top, bottom

    !
    ! Stupidity check: must be defined and must be layers, not levels
    !
    fraction = real_missing
    if(.not. defined(layerThin) .or. .not. defined(layerThick))then
      call set_error('One of layres is undefined','fu_vert_overlap_fraction')
      return
    endif
    if(.not. fu_if_layer(layerThin))then
      call set_error('Thin layer is not a layer','fu_vert_overlap_fraction')
      return
    endif
    if(.not. fu_if_layer(layerThick))then
      call set_error('Thick layer is not a layer','fu_vert_overlap_fraction')
      return
    endif
    
    !
    ! Same type ? If not - project layerThin to the layerThick system
    !
    if(layerThin%leveltype == layerThick%leveltype)then
      layerTmp = layerThin
    else
      call set_vertical(fu_lower_boundary_of_layer(layerThick), vertTmp)
      if(error)return
      layerTmp = fu_set_layer_between_two( &
         & fu_level_to_vertical_crude(fu_upper_boundary_of_layer(layerThin), vertTmp), &
         & fu_level_to_vertical_crude(fu_lower_boundary_of_layer(layerThin), vertTmp))
    endif
    if(error)return

    !
    ! Now they are of the same type - just compute the overlap. 
    !
    select case(layerThick%leveltype)
      case(layer_btw_2_altitude, &
         & layer_btw_2_height, &
         & layer_btw_2_depth)

        if(layerTmp%a <= layerThick%b .or. layerTmp%b >= layerThick%a)then
          fraction = 0.  ! No overlap
        else
          fraction = (min(layerTmp%a,layerThick%a) - &
                                    & max(layerTmp%b,layerThick%b)) /  &
                                   & (layerTmp%a-layerTmp%b)
        endif

      case(layer_btw_2_pressure, &
         & layer_btw_2_sigma)

        if(layerTmp%a >= layerThick%b .or. layerTmp%b <= layerThick%a)then
          fraction = 0.  ! No overlap
        else
          fraction = (max(layerTmp%a,layerThick%a) - &
                                    & min(layerTmp%b,layerThick%b)) /  &
                                   & (layerTmp%b-layerTmp%a)
        endif

      case(layer_btw_2_hybrid)

        if(layerTmp%hybrid_coeff_known)then
          if(layerTmp%a + std_pressure_sl*layerTmp%b >= &
           & layerThick%a2 + std_pressure_sl*layerThick%b2 &
           & .or. &
           & layerTmp%a2 + std_pressure_sl*layerTmp%b2 <= &
           & layerThick%a + std_pressure_sl*layerThick%b)then
            fraction = 0.                       ! No overlap
          else
            fraction = (max(layerTmp%a + std_pressure_sl*layerTmp%b, &
                                          & layerThick%a + std_pressure_sl*layerThick%b) - &
                                      & min(layerTmp%a2 + std_pressure_sl*layerTmp%b2, &
                                          & layerThick%a2 + std_pressure_sl*layerThick%b2)) &
                                     & / &
                                     & (layerTmp%a + std_pressure_sl*layerTmp%b - &
                                      & layerTmp%a2 - std_pressure_sl*layerTmp%b2)
          endif

        else ! If hybrid coeffs known

          if(layerTmp%number <= layerThick%number2 .or. &
           & layerTmp%number2 >= layerThick%number)then
            fraction = 0.  ! No overlap
          else
            fraction = (min(layerTmp%number,layerThick%number) - &
                                      & max(layerTmp%number2,layerThick%number2)) /  &
                                     & (layerTmp%number2-layerTmp%number)
          endif
          
        endif ! hybrid coeffs known
        
      case default
        call report(layerThick)
        call set_error('Non-supported layer type','fu_vert_overlap_fraction')
        return
    end select 

  end function fu_vert_overlap_fraction_old


  ! ***************************************************************


  INTEGER FUNCTION fu_level_count(levels)

    ! Description:
    ! For a vector of levels returns the number of defined levels.
    IMPLICIT NONE

    ! Imported parameters with intent(in):
    TYPE(silja_level), DIMENSION(:), INTENT(in) :: levels

    ! Local declarations:
    INTEGER :: i

    fu_level_count = 0

    DO i = 1, SIZE(levels)
      IF (defined(levels(i))) fu_level_count = fu_level_count + 1
    END DO

  END FUNCTION fu_level_count


  !*****************************************************************

  integer function fu_Nlevs_of_vertical(vert)
    !
    ! Returns the number of the levels in the vertical structure
    !
    IMPLICIT NONE

    !Imported parametsr with the intent IN
    type(silam_vertical), intent(in) :: vert

    if(defined(vert))then
      fu_Nlevs_of_vertical = vert%NLevs
    else
      fu_Nlevs_of_vertical = 0
    endif
  end function fu_Nlevs_of_vertical


  !*****************************************************************

  logical function fu_verts_comparable(vert1, vert2)
    !
    ! Comparison (apart from ==) of two verticals is possible if they are
    ! in similar types of the vertical co-ordinates, which is checked here.
    ! Term "comparable" about the vertical means that one can decide which
    ! of the vertical structures cover wider range and which one has better
    ! resolution. This should be possible without involvment of real data, 
    ! just looking at the list of levels.
    !
    implicit none

    ! Imported parameters with intent IN
    type(silam_vertical), intent (in) :: vert1, vert2

    !
    ! First, trivial part
    !
    if(.not.defined(vert1))then
      fu_verts_comparable = .not.defined(vert2)
      return
    endif
    if(.not.defined(vert2))then
      fu_verts_comparable = .false. ! undefined vert1 alread considered
      return
    endif
    if (vert1%vert_type == vert2%vert_type) then
      fu_verts_comparable = .true.
      return
    end if
    !
    ! Now what to do if the types are not exactly the same. Then comparability
    ! simply means that level heights can be expressed in similar units.
    ! Here we skip all conditions requiring equality of the types
    !
    select case (vert1%vert_type)
      case(surface, mean_sea)
        fu_verts_comparable = (vert2%vert_type == surface) .or. &
                           &  (vert2%vert_type == mean_sea)

      case(constant_pressure, &
         & layer_btw_2_pressure, &
         & hybrid, &
         & layer_btw_2_hybrid) 
        fu_verts_comparable = (vert2%vert_type == constant_pressure) .or. &
                            & (vert2%vert_type == hybrid) .or. &
                            & (vert2%vert_type == layer_btw_2_pressure) .or. &
                            & (vert2%vert_type == layer_btw_2_hybrid)

      case(constant_altitude, &
         & layer_btw_2_altitude, &
         & constant_height, &
         & layer_btw_2_height)
        fu_verts_comparable = (vert2%vert_type == constant_altitude) .or. &
                            & (vert2%vert_type == layer_btw_2_altitude).or. &
                            & (vert2%vert_type == constant_height).or. &
                            & (vert2%vert_type == layer_btw_2_height)

      case(sigma_level, layer_btw_2_sigma)
        fu_verts_comparable = (vert2%vert_type == sigma_level).or. &
                            & (vert2%vert_type == layer_btw_2_sigma)

      case(depth_level, layer_btw_2_depth)
        fu_verts_comparable = (vert2%vert_type == depth_level).or. &
                            & (vert2%vert_type == layer_btw_2_depth)
      case default
        call msg_warning('Unknown vertical type','fu_verts_comparable')
        fu_verts_comparable = .false. ! Unknown vertical
    end select


  end function fu_verts_comparable



  !*****************************************************************

  logical function fu_level_belongs_to_vertical(level, vert)
    !
    ! Checks if given level belongs to the vertical.
    !
    implicit none

    ! Imported parameters with intent IN
    type(silja_level), intent(in) :: level
    type(silam_vertical), intent(in) :: vert

    ! Local variables
    integer :: iL
!    logical :: equ
    
    fu_level_belongs_to_vertical = .false.

    if(.not.defined(level))then
      call set_error('Undefined level','fu_level_belongs_to_vertical')
      return
    end if
    if(.not.defined(vert))then
      call set_error('Undefined vertical','fu_level_belongs_to_vertical')
      return
    end if
    if(vert%vert_type == any_level)then    ! Accept all
      fu_level_belongs_to_vertical = .true.
      return
    end if

    if(vert%vert_type /= level%leveltype) return

    do iL = 1, vert%NLevs
  !  if(vert%levs(iL) == level)then
      if(fu_cmp_levs_eq(vert%levs(iL),level))then
        fu_level_belongs_to_vertical = .true.
        return
      end if
    end do

  end function fu_level_belongs_to_vertical


  !*********************************************************************

  real function fu_level_index_in_levels(levelIn, levs) result(fIndex)
    !
    ! Returns the relative index of the level in a vertical structure
    ! defined by a set of levels.
    ! Level must be of the same type as levs, so, no interpolation
    ! between the different types of levels is done here. 
    ! Important: levels are defined via their central points, so index
    ! should be computed for the correponding thick layers.
    ! Sorting of the levels is unknown
    ! If the level is outside the levs structure - it gets the first or last
    ! index correspondingly
    !
    ! All inuts: SI
    !
    implicit none
    !
    ! Imported parameters
    type(silja_level), intent(in) :: levelIn
    type(silja_level), dimension(:), intent(in) :: levs

    ! Local declarations
    integer :: i, nLevs
    real :: border, pressure, prTmp1, prTmp2
    type(silja_level) :: level

    fIndex = -1
    !
    ! Get the number of levels in the vertical structure
    !
    nLevs=size(levs)
    do i=1,size(levs)
      if(.not. levs(i)%defined == silja_true)then
        nLevs = i-1
        exit
      endif
    end do
    !
    ! If the level and the vertical structure are thick layers of similar type,
    ! we will compare their central points.
    !
    if(fu_if_layer(levelIn))then
      level = fu_central_level_of_layer(levelIn)
    else
      level = levelIn
    endif

    !
    ! Stupidity check. Levels must be the same OR be compatible. For example, 
    ! thick z-system based layers are compatible with z-level (not vise versa).
    ! Two thick systems are compatible too via their central points
    !
    if(level%levelType /= levs(1)%levelType .and. &
     & level%levelType /= fu_leveltype(fu_central_level_of_layer(levs(1))))then
      call msg_warning('Can not handle levels of different types, levelIn, levels(1):', &
                     & 'fu_level_index_in_levels')
      call report(levelIn)
      call report(levs(1))
      call set_error('Can not handle levels of substantially different types', &
                   & 'fu_level_index_in_levels')
      return
    endif

    !
    ! Actual reprojection starts
    !
    SELECT CASE (level%leveltype)

      CASE (constant_height) ! levels may be height or layers btw 2 heights

        border =0.
        if(nLevs == 1 .or. levs(1)%a < levs(nLevs)%a)then  ! First level is near surface
          if(levs(1)%leveltype == constant_height)then
            !
            ! Have to create borders, etc.
            !
            do i=1,nLevs
              border = 2.*levs(i)%a - border
              if(i < nLevs) then
                if(border > levs(i+1)%a) border = 0.5*(levs(i)%a + levs(i+1)%a)
              endif
              IF (level%a < border) THEN
                fIndex = real(i) - 0.5 + 0.5*(level%a + border - 2.*levs(i)%a) / &
                                           & (border - levs(i)%a)
                RETURN
              END IF
            end do

          elseif(levs(1)%leveltype == layer_btw_2_height)then
            !
            ! Borders are available immediately - bottom is b and top is a
            !
            do i=1,nLevs
              if(level%a <= levs(i)%a)then ! come to the layer for the first time
                fIndex = real(i) - 0.5 + (level%a - levs(i)%b) / (levs(i)%a - levs(i)%b)
                return
              end if
            end do
            
          else
            call set_error('Can not handle these levels','fu_level_index_in_levels')
            call report(level)
            call report(levs(1))
            return
          endif
          fIndex = real(nLevs) + 0.499999

        else           ! First level is the highest

          if(levs(1)%leveltype == constant_height)then
            !
            ! Have to create the borders
            !
            do i=nLevs,1,-1
              border = 2.*levs(i)%a - border
              if(i > 1) then
                if(border > levs(i-1)%a) border = 0.5*(levs(i)%a + levs(i-1)%a)
              endif
              IF (level%a < border) THEN
                fIndex = real(i) - 0.5 + 0.5*(level%a + border - 2.*levs(i)%a) / &
                                           & (border - levs(i)%a)
                RETURN
              END IF
            end do
          elseif(levs(1)%leveltype == layer_btw_2_height)then
            !
            ! Borders are available immediately - bottom is b and top is a
            !
            do i=nLevs, 1,-1
              if(level%a <= levs(i)%a)then ! come to the layer for the first time
                fIndex = real(i) - 0.5 + (level%a - levs(i)%b) / (levs(i)%a - levs(i)%b)
                return
              end if
            end do
            
          else
            call report(level)
            call report(levs(1))
            call set_error('Can not handle these levels','fu_level_index_in_levels')
            return
          endif
          fIndex = 0.5000001
        endif

      CASE (constant_pressure)

        if(nLevs == 1 .or. levs(1)%a > levs(nLevs)%a)then   ! Level 1 is close to surface
          if(levs(1)%leveltype == constant_pressure)then
            border =max(std_pressure_sl, levs(1)%a)
            do i=1,nLevs
              border = 2.*levs(i)%a - border
              if(i < nLevs) then
                if(border < levs(i+1)%a) border = 0.5*(levs(i)%a + levs(i+1)%a)
              endif
              IF (level%a > border) THEN
                fIndex = real(i) - 0.5 + 0.5*(level%a + border - 2.*levs(i)%a) / &
                                           & (border - levs(i)%a)
                RETURN
              END IF
            end do
          elseif(levs(1)%leveltype == layer_btw_2_pressure)then
            !
            ! Borders are available immediately - bottom is b and top is a
            !
            do i=1,nLevs
              if(level%a >= levs(i)%a)then ! come to the layer for the first time
                fIndex = real(i) - 0.5 + (level%a - levs(i)%b) / (levs(i)%a - levs(i)%b)
                return
              end if
            end do
            
          else
            call report(level)
            call report(levs(1))
            call set_error('Can not handle these levels','fu_level_index_in_levels')
            return
          endif
          fIndex = real(nLevs) + 0.499999

        else  ! Level 1 is the highest

          if(levs(1)%leveltype == constant_pressure)then
            border =max(std_pressure_sl, levs(nLevs)%a)
            do i=nLevs,1,-1
              border = 2.*levs(i)%a - border
              if(i > 1) then
                if(border < levs(i-1)%a) border = 0.5*(levs(i)%a + levs(i-1)%a)
              endif
              IF (level%a > border) THEN
                fIndex = real(i) - 0.5 + 0.5*(level%a + border - 2.*levs(i)%a) / &
                                           & (border - levs(i)%a)
                RETURN
              END IF
            end do
          elseif(levs(1)%leveltype == layer_btw_2_pressure)then
            !
            ! Borders are available immediately - bottom is b and top is a
            !
            do i=nLevs, 1, -1
              if(level%a >= levs(i)%a)then ! come to the layer for the first time
                fIndex = real(i) + 0.5 - (level%a - levs(i)%b) / (levs(i)%a - levs(i)%b)
                return
              end if
            end do
            
          else
            call report(level)
            call report(levs(1))
            call set_error('Can not handle these levels','fu_level_index_in_levels')
            return
          endif

          fIndex = 0.500001
        endif

      CASE (sigma_level)

        if(nLevs == 1 .or. levs(1)%a > levs(nLevs)%a)then   ! Level 1 is close to surface
          if(levs(1)%leveltype == sigma_level)then
            border = max(1.0, levs(1)%a)
            do i=1,nLevs
              border = 2.*levs(i)%a - border
              if(i < nLevs) then
                if(border < levs(i+1)%a) border = 0.5*(levs(i)%a + levs(i+1)%a)
              endif
              IF (level%a > border) THEN
                fIndex = real(i) - 0.5 + 0.5*(level%a + border - 2.*levs(i)%a) / &
                                           & (border - levs(i)%a)
                RETURN
              END IF
            end do
          elseif(levs(1)%leveltype == layer_btw_2_sigma)then
            !
            ! Borders are available immediately - bottom is b and top is a
            !
            do i=1,nLevs
              if(level%a >= levs(i)%a)then ! come to the layer for the first time
                fIndex = real(i) - 0.5 + (level%a - levs(i)%b) / (levs(i)%a - levs(i)%b)
                return
              end if
            end do
            
          else
            call report(level)
            call report(levs(1))
            call set_error('Can not handle these levels','fu_level_index_in_levels')
            return
          endif
          fIndex = real(nLevs) + 0.499999

        else  ! Level 1 is the highest

          if(levs(1)%leveltype == sigma_level)then
            border =max(1.0, levs(nLevs)%a)
            do i=nLevs,1,-1
              border = 2.*levs(i)%a - border
              if(i > 1) then
                if(border < levs(i-1)%a) border = 0.5*(levs(i)%a + levs(i-1)%a)
              endif
              IF (level%a > border) THEN
                fIndex = real(i) - 0.5 + 0.5*(level%a + border - 2.*levs(i)%a) / &
                                           & (border - levs(i)%a)
                RETURN
              END IF
            end do
          elseif(levs(1)%leveltype == layer_btw_2_sigma)then
            !
            ! Borders are available immediately - bottom is b and top is a
            !
            do i=nLevs, 1, -1
              if(level%a >= levs(i)%a)then ! come to the layer for the first time
                fIndex = real(i) + 0.5 - (level%a - levs(i)%b) / (levs(i)%a - levs(i)%b)
                return
              end if
            end do
            
          else
            call report(level)
            call report(levs(1))
            call set_error('Can not handle these levels','fu_level_index_in_levels')
            return
          endif

          fIndex = 0.500001
        endif

      case(hybrid)
        !
        ! If hybrid coefficients are known - use reference sea level pressure to 
        ! compute pressure at levels and then compare them as above. If not - just
        ! compare the level numbers. Since there can be holes in enumeration
        ! of levs - do not use the number itself as the index.
        !
        ! Be careful - nobody knows the sorting of the hybrid levels. It may happen
        ! that hybrid level n is below level n+1.
        !
        if(level%hybrid_coeff_known .and. all(levs(1:nLevs)%hybrid_coeff_known))then
          pressure = level%a + level%b * std_pressure_sl

          PrTmp1 = levs(1)%a+levs(1)%b*std_pressure_sl
          PrTmp2 = levs(nLevs)%a+levs(nLevs)%b*std_pressure_sl

          if(nLevs == 1 .or. PrTmp1 > PrTmp2)then  ! Level 1 starts from the surface

            if(levs(1)%leveltype == hybrid)then
              !
              ! Have to create borders
              !
              PrTmp2 = PrTmp1
              border = max(PrTmp2, std_pressure_sl)
              do i=1,nLevs
                prTmp1 = prTmp2
                border = 2.*prTmp1 - border
                if(i < nLevs) then
                  prTmp2 = levs(i+1)%a+levs(i+1)%b*std_pressure_sl
                  if(border < prTmp2) border = 0.5*(prTmp1 + prTmp2)
                endif
                if(pressure > border)then
                  fIndex = real(i) - 0.5 + 0.5 * &
                        & (pressure + border - 2.*prTmp1) / (border - prTmp1)
                  return
                endif
              end do
            elseif(levs(1)%leveltype == layer_btw_2_hybrid)then
              !
              ! Borders are available immediately - top is a, b, bottom is a2, b2.
              !
              do i=1,nLevs
                ! prTmp1 = the top pressure
                PrTmp1 = levs(i)%a + levs(i)%b * std_pressure_sl
                if(pressure >= PrTmp1)then   ! below layer top
                  fIndex = real(i) + 0.5 &
                       & - (pressure - PrTmp1) &
                       & / (levs(i)%a2 + levs(i)%b2 * std_pressure_sl - PrTmp1)
                  return
                end if
              end do
            
            else
              call report(level)
              call report(levs(1))
              call set_error('Can not handle these levels','fu_level_index_in_levels')
              return
            endif
            fIndex = real(nLevs) + 0.499999

          else     ! Level 1 is the up-most

            if(levs(1)%leveltype == hybrid)then
              !
              ! Have to create borders
              !
              border = max(PrTmp2, std_pressure_sl)
              do i=nLevs,1,-1
                prTmp1 = prTmp2
                border = 2.*prTmp1 - border
                if(i > 1) then
                  prTmp2 = levs(i-1)%a+levs(i-1)%b*std_pressure_sl
                  if(border < prTmp2) border = 0.5*(prTmp1 + prTmp2)
                endif
                if(pressure > border)then
                  fIndex = real(i) - 0.5 + 0.5 * &
                        & (pressure + border - 2.*prTmp1) / (border - prTmp1)
                  return
                endif
              end do
            elseif(levs(1)%leveltype == layer_btw_2_hybrid)then
              !
              ! Borders are available immediately - bottom is *2 and top is without index
              !
              do i=nLevs, 1, -1
                PrTmp1 = levs(i)%a2 + levs(i)%b2 * std_pressure_sl ! layer bottom
                if(pressure <= PrTmp1)then     ! come to the layer for the first time
                  fIndex = real(i) + 0.5 - (prTmp1 - pressure) / &
                         & (PrTmp1 - levs(i)%a + levs(i)%b * std_pressure_sl)
                  return
                end if
              end do
            
            else
              call report(level)
              call report(levs(1))
              call set_error('Can not handle these levels','fu_level_index_in_levels')
              return
            endif
            fIndex = 0.500001
          endif   ! Sorting

        else ! If coefficients are unknown - we have to assume that they are the same

          if(levs(1)%leveltype == hybrid)then
            !
            ! Have to create borders
            !
            border=0.
            do i=1,nLevs
              border = 2.*levs(i)%number - border
              if(i < nLevs) then
                if(border > levs(i+1)%number) border = 0.5*(levs(i)%number + levs(i+1)%number)
              endif
              if(level%number <= border)then
                fIndex = real(i) - 0.5 + 0.5 * &
                        & (level%number + border - 2.*levs(i)%number) / &
                        & (border - levs(i)%number)
                return
              endif
            end do
          elseif(levs(1)%leveltype == layer_btw_2_hybrid)then
            !
            ! Borders are available immediately - bottom is *2 and top is without index
            !
            do i=1, nLevs
              if(level%number <= levs(i)%number)then  ! come to the layer for the first time
                fIndex = real(i) - 0.5 + (level%number - levs(i)%number2) / &
                       & (levs(i)%number - levs(i)%number2)
                return
              end if
            end do
            
          else
            call report(level)
            call report(levs(1))
            call set_error('Can not handle these levels','fu_level_index_in_levels')
            return
          endif
        fIndex = real(nLevs) + 0.499999
        endif

      CASE default
        CALL msg_warning('sorry cannot handle the level type','fu_level_index_in_levels')
        call report(level)
        CALL set_error('sorry cannot handle the level type','fu_level_index_in_levels')
      END SELECT

  end function fu_level_index_in_levels


  !*********************************************************************

  real function fu_level_index_in_vertical(level, vert) result(fIndex)
    !
    ! Same as previous but the vertical structure is given via vertical
    ! type, not just a set of levels
    !
    implicit none

    ! Imported parameters
    type(silja_level), intent(in) :: level
    type(silam_vertical), intent(in) :: vert

    if(.not.defined(vert))then
      call set_error('Undefined vertical given','fu_level_index_in_vertical')
      fIndex = -1
      return
    endif
    fIndex = fu_level_index_in_levels(level, vert%levs)

  end function fu_level_index_in_vertical


  !***********************************************************************

  function fu_level_to_vertical_crude(levelIn, vert) result(levelOut)
    !
    ! Projects level to some vertical structure.
    ! If level type and verttical type are the same - nothing happens,
    ! the same level is returned. But if they are not the same, the
    ! level is interpolated to the vetical. Evidently, this may require
    ! knowledge about surface pressure and temperature, but for crude
    ! interpolation they are substituted with some standard values
    !
    implicit none

    ! return value
    type(silja_level) :: levelOut

    ! Imported parameters
    type(silja_level), intent(in) :: levelIn
    type(silam_vertical), intent(in)  :: vert

    ! Local variables
    integer :: iLev
    real :: fTmp, fTmp_1, fPressure
    logical :: start_from_surface

    !
    ! Stupidity check
    !
    !call report(levelIn)
    if(.not. defined(levelIn))then
      call set_error('Undefined level','level_to_vertical_crude')
      return
    endif
    if(.not. defined(vert))then
      call set_error('Undefined vertical','level_to_vertical_crude')
      return
    endif
    if(levelIn%leveltype == vert%vert_type)then
      levelOut = levelIn
      return
    endif
    !
    ! Now - project each-to-each
    !
    levelOut = level_missing ! By default they are non-comparable
    select case(vert%vert_type)

      case(surface) ! Require the levels to be close to the surface
        select case(levelIn%leveltype)
          case(mean_sea)
            levelOut = surface_level
          case(constant_pressure)
            if(levelIn%a >= 97500.) levelOut = surface_level
          case(constant_altitude)
            if(levelIn%a <= 50.) levelOut = surface_level
          case(constant_height)
            if(levelIn%a <= 50.) levelOut = surface_level
          case(sigma_level)
            if(levelIn%a > 0.975) levelOut = surface_level
          case(hybrid)
            if(levelIn%number == 1) levelOut = surface_level
          case default
            call report(vert)
            call report(levelIn)
            call set_error('Can not project above level to above vertical', &
                         & 'fu_level_to_vertical_crude')
        end select

      case(mean_sea)
        select case(levelIn%leveltype)
          case(surface)
            levelOut = mean_sea_level
          case(constant_pressure)
            if(levelIn%a >= 97500.) levelOut = mean_sea_level
          case(constant_altitude)
            if(levelIn%a <= 50.) levelOut = mean_sea_level
          case(constant_height)
            if(levelIn%a <= 50.) levelOut = mean_sea_level
          case(sigma_level)
            if(levelIn%a > 0.975) levelOut = mean_sea_level
          case(hybrid)
            if(levelIn%number == 1) levelOut = mean_sea_level
          case default
            call report(vert)
            call report(levelIn)
            call set_error('Can not project above level to above vertical', &
                         & 'level_to_vertical_crude')
        end select

      case(constant_altitude, layer_btw_2_altitude)
        select case(levelIn%leveltype)
          case(surface, mean_sea)
            levelOut = fu_set_constant_altitude_level(0.)
          case(constant_pressure)
            levelOut = fu_set_constant_altitude_level(&
                               & gas_constant_dryair * temperature_ref / g * &
                               & log(std_pressure_sl / fu_pr_level_pressure(levelIn)))
          case(constant_height, constant_altitude)
            levelOut = fu_set_constant_altitude_level(fu_level_height(levelIn))
          case(sigma_level)
            levelOut = fu_set_constant_altitude_level(&
                               & gas_constant_dryair * temperature_ref / g * &
                               & (-log(fu_sigma_level_sigma(levelIn))))
          case(hybrid)
            levelOut = fu_set_constant_altitude_level(&
                               & gas_constant_dryair * temperature_ref* &
                               & log(std_pressure_sl / &
                                   & fu_hybrid_level_pressure(levelIn, std_pressure_sl)))
          case(depth_level)
            levelOut = fu_set_constant_altitude_level(-fu_depth_level_depth(levelIn))
          case default
            call report(vert)
            call report(levelIn)
            call set_error('Can not project above level to above vertical', &
                         & 'level_to_vertical_crude')
        end select

      case(constant_height, layer_btw_2_height)
        select case(levelIn%leveltype)
          case(surface, mean_sea)
            levelOut = fu_set_constant_height_level(0.)
          case(constant_pressure)
            levelOut = fu_set_constant_height_level(&
                               & gas_constant_dryair * temperature_ref / g * &
                               & log(std_pressure_sl / fu_pr_level_pressure(levelIn)))
          case(constant_altitude, constant_height)
            levelOut = fu_set_constant_height_level(fu_level_height(levelIn))
          case(sigma_level)
            levelOut = fu_set_constant_height_level(&
                               & gas_constant_dryair * temperature_ref / g * &
                               & (- log(fu_sigma_level_sigma(levelIn))))
          case(hybrid)
            levelOut = fu_set_constant_height_level(&
                               & gas_constant_dryair * temperature_ref / g * &
                               & log(std_pressure_sl / &
                                   & fu_hybrid_level_pressure(levelIn, std_pressure_sl)))
          case(depth_level)
            levelOut = fu_set_constant_height_level(-fu_depth_level_depth(levelIn))
          case default
            call report(vert)
            call report(levelIn)
            call set_error('Can not project above level to above vertical', &
                         & 'fu_level_to_vertical_crude')
        end select

      case(constant_pressure, layer_btw_2_pressure)
        select case(levelIn%leveltype)
          case(surface, mean_sea) ! No topography...
            levelOut = fu_set_pressure_level(std_pressure_sl)
          case(constant_altitude)
            levelOut = fu_set_pressure_level(std_pressure_sl / &
                     & exp(fu_level_altitude(levelIn)/(gas_constant_dryair*temperature_ref/g)))
          case(constant_height)
            levelOut = fu_set_pressure_level(std_pressure_sl / &
                     & exp(fu_level_height(levelIn)/(gas_constant_dryair*temperature_ref/g)))
          case(constant_pressure)
            levelOut = levelIn
          case(sigma_level)
            levelOut = fu_set_pressure_level(std_pressure_sl * fu_sigma_level_sigma(levelIn))
          case(hybrid)
            levelOut = fu_set_pressure_level(fu_hybrid_level_pressure(levelIn, std_pressure_sl))
          case default
            call report(vert)
            call report(levelIn)
            call set_error('Can not project above level to above vertical', &
                         & 'fu_level_to_vertical_crude')
        end select

      case(hybrid, layer_btw_2_hybrid)
        !
        ! For the hybrid system, we get the pressure of the levelIn and then 
        ! find out which two hybrid levels surround it
        !
        select case(levelIn%leveltype)
          case(surface, mean_sea) ! No topography...
            fPressure = std_pressure_sl
          case(constant_altitude)
            fPressure = std_pressure_sl / &
                   & exp(fu_level_altitude(levelIn)/(gas_constant_dryair*temperature_ref/g))
          case(constant_height)
            fPressure = std_pressure_sl / &
                   & exp(fu_level_height(levelIn)/(gas_constant_dryair*temperature_ref/g))
          case(sigma_level)
            fPressure = std_pressure_sl * fu_sigma_level_sigma(levelIn)
          case(hybrid)
            fPressure = fu_hybrid_level_coeff_a(levelIn) + std_pressure_sl * fu_hybrid_level_coeff_b(levelIn)
          case(constant_pressure)
            fPressure = fu_pr_level_pressure(levelIn)
          case default
            call report(vert)
            call report(levelIn)
            call set_error('Can not project above level to above vertical', &
                         & 'fu_level_to_vertical_crude')
        end select

        !
        ! Hybrid levels can start from the surface and from the top of the domain,
        ! so there has to be two options here
        !
        fTmp_1 = fu_hybrid_level_pressure(vert%levs(1),std_pressure_sl)
        if (fu_nbrOfLevels(vert) == 1) then
          start_from_surface = .true. 
        else
          start_from_surface = fTmp_1 > fu_hybrid_level_pressure(vert%levs(2),std_pressure_sl)
        end if
        if(start_from_surface)then
          !
          ! Start from the surface
          !
          do iLev = 2,vert%nLevs
            if(.not. vert%levs(iLev)%hybrid_coeff_known)then
              call set_error('Cannot handle hybrid vertical with unknown coefficients', &
                           & 'fu_level_to_vertical_crude')
              return
            endif
            fTmp = fu_hybrid_level_pressure(vert%levs(iLev),std_pressure_sl)
            if(fPressure >= fTmp)then  ! our pressure in-between this and previous levels
              if(fPressure - fTmp > fTmp_1 - fPressure)then
                levelOut = vert%levs(iLev-1)
              else
                levelOut = vert%levs(iLev)
              endif
              return
            endif  ! Found the level higher than the levelIn
            fTmp_1 = fTmp  ! Store the current-level pressure for the next cycle
          end do ! levels over the vertical
          !
          ! Our level is closer to the surface than the lowest model level. Issue warning
          ! and force the last level to be the output one
          !
          call msg_warning('Level is higher than the highest one in the vertical', &
                         & 'fu_level_to_vertical_crude')
        else
          !
          ! Start from the top of the atmosphere
          !
          do iLev = 2,vert%nLevs
            if(.not. vert%levs(iLev)%hybrid_coeff_known)then
              call set_error('Cannot handle hybrid vertical with unknown coefficients', &
                           & 'fu_level_to_vertical_crude')
              return
            endif
            fTmp = fu_hybrid_level_pressure(vert%levs(iLev),std_pressure_sl)
            if(fPressure < fTmp)then  ! our pressure in-between this and previous levels
              if(fPressure - fTmp < fTmp_1 - fPressure)then
                levelOut = vert%levs(iLev-1)
              else
                levelOut = vert%levs(iLev)
              endif
              return
            endif
          end do ! levels over the vertical
          !
          ! Our level is closer to the surface than the lowest model level. Issue warning
          ! and force the last level to be the output one
          !
          call msg_warning('Level is closer to surface than the lowest one in the vertical', &
                         & 'fu_level_to_vertical_crude')

        endif ! sorting order in the vertical

        call report(vert)
        call report(levelIn)
        levelOut = vert%levs(vert%nLevs)

      case(sigma_level, layer_btw_2_sigma)
        select case(levelIn%leveltype)
          case(surface, mean_sea) ! No topography...
            levelOut = fu_set_sigma_level(1.)
           case(constant_altitude)
            levelOut = fu_set_sigma_level(1. / &
                     & exp(fu_level_altitude(levelIn)/(gas_constant_dryair*temperature_ref/g)))
          case(constant_height)
            levelOut = fu_set_sigma_level(1. / &
                     & exp(fu_level_height(levelIn)/(gas_constant_dryair*temperature_ref/g)))
          case(constant_pressure)
            levelOut = fu_set_sigma_level(fu_pr_level_pressure(levelIn) / std_pressure_sl)
          case(hybrid)
            levelOut = fu_set_sigma_level(fu_hybrid_level_pressure(levelIn, std_pressure_sl) / &
                                        & std_pressure_sl)
          case default
            call report(vert)
            call report(levelIn)
            call set_error('Can not project above level to above vertical', &
                         & 'fu_level_to_vertical_crude')
        end select

      case default
        call set_error('The following vertical is not handled yet','fu_level_to_vertical_crude')
        call report(vert)
        return
    end select

  end function fu_level_to_vertical_crude

  !************************************************************************************

  subroutine vert_interp_data_crude(vertical_from, vertical_to, data_column, data_surf, relief_height_m)
    !
    ! Return reference values for the meteorological quantities needed
    ! by fu_project_level. Can be used to perform vertical crude
    ! interpolation.
    implicit none
    type(silam_vertical), intent(in) :: vertical_from, & ! the levels are projected from this vertical
                                      & vertical_to ! ...to this
    real, dimension(:), intent(out) :: data_column
    real, intent(out) :: data_surf
    real, intent(in), optional :: relief_height_m

    integer :: ilev, nlevs_to, quantity_3d, quantity_2d
    real :: z, pres, density, temperature, press_ratio, pres_bottom, pres_top, z_top, z_bottom, relief
    type(silja_level) :: level, bottom, top
    logical :: have_layers

    if(present(relief_height_m))then
      relief = relief_height_m
    else
      relief = 0.
    endif
    nlevs_to = fu_nbrOfLevels(vertical_to)
    have_layers = fu_if_layer(fu_level(vertical_to,1))
    call projection_input_needs(vertical_from, vertical_to, quantity_3d, quantity_2d)
    if (error) return
    
    select case(quantity_3d)
      case (int_missing)
        !print *, 'nothing'
        continue

      case (pressure_flag)
        ! 
        ! Use US standard atmosphere to get presssure and barostatic equation to take care of relief
        !
        if (have_layers) then
          do ilev = 1, nlevs_to
            top = fu_upper_boundary_of_layer(fu_level(vertical_to, ilev))
            bottom = fu_lower_boundary_of_layer(fu_level(vertical_to, ilev))
            z_bottom = fu_level_height(bottom)
            z_top = fu_level_height(top)
            call us_standard_atmosphere(z_bottom + relief, density, press_ratio, temperature)
            data_column(ilev) = std_pressure_sl*press_ratio
            call us_standard_atmosphere(z_top + relief, density, press_ratio, temperature)
            data_column(ilev+1) = std_pressure_sl*press_ratio
          end do
        else
          do ilev = 1, nlevs_to
            z = fu_level_height(fu_level(vertical_to, ilev, .true.))
            call us_standard_atmosphere(z + relief, density, press_ratio, temperature)
            data_column(ilev) = std_pressure_sl*press_ratio
            !print *, z, data_column(ilev)
          end do
        end if

      case (height_flag)
        !
        ! Use binary search to invert the US standard atmosphere. Note that relief has to be subtracted
        ! from the altitude that falls out of MSL pressure.
        !
        if (have_layers) then
          do ilev = 1, nlevs_to
            top = fu_upper_boundary_of_layer(fu_level(vertical_to, ilev))
            bottom = fu_lower_boundary_of_layer(fu_level(vertical_to, ilev))
            if (top%leveltype == constant_pressure) then
              pres_top = fu_pr_level_pressure(top)
              pres_bottom = fu_pr_level_pressure(bottom)
            else if (top%leveltype == hybrid) then
              pres_top = fu_hybrid_level_pressure(top, &
                                                & std_pressure_sl * exp(- relief / std_height_scale))
              pres_bottom = fu_hybrid_level_pressure(bottom, &
                                                & std_pressure_sl * exp(- relief / std_height_scale))
            else
              call set_error('Cannot handle this level/parameter combination', 'vert_interp_data_crude')
            end if
            data_column(ilev) = fu_height_for_press(pres_bottom)
            data_column(ilev+1) = fu_height_for_press(pres_top)
          end do
        else
          do ilev = 1, nlevs_to
            level = fu_level(vertical_to, ilev, .true.)
            if (level%leveltype == constant_pressure) then
              pres = fu_pr_level_pressure(level)
            else if (level%leveltype == hybrid) then
              pres = fu_hybrid_level_pressure(level, std_pressure_sl * exp(- relief/std_height_scale))
            else
              call set_error('Cannot handle this level/parameter combination', 'vert_interp_data_crude')
            end if
            data_column(ilev) = fu_height_for_press(pres)
          end do
        end if

      case default
        call msg('Parameter requested: ' // fu_quantity_short_string(quantity_3d), quantity_3d) 
        call set_error('Cannot produce this parameter', 'vert_intepr_data_crude')
    end select
    
    select case(quantity_2d)
      case (int_missing)
        data_surf = real_missing
      case (surface_pressure_flag)
        data_surf = std_pressure_sl * exp(- relief / std_height_scale)
      case (relief_height_flag)
        data_surf = relief
      case default
        call msg('Parameter requested: ' // fu_quantity_short_string(quantity_2d), quantity_2d) 
        call set_error('Cannot produce this parameter', 'vert_intepr_data_crude')
    end select
  
  end subroutine vert_interp_data_crude

  
  !************************************************************************************
  
  real function fu_height_for_press(pressure) result(height)
    ! 
    ! Find the height corresponding to a given pressure according to
    ! the US standard atmosphere. Uses binary search to on the forward
    ! us_standard_atmosphere function -> can be slow.
    !
    implicit none
    real, intent(in) :: pressure
    real :: density_rel, height_up, height_down, press_ratio, &
         & tempr_ratio, dens_ratio, pressure_iter, pressure_up, pressure_down
    real, parameter :: max_height = 85e3, rel_tol = 0.0001, min_pressure = 0.36342, H_mono = 6700.0
    integer :: itercount
    integer, parameter :: max_iter = 50

    if (pressure < min_pressure) then
      height = max_height
      return
    end if
    if(pressure >= std_pressure_sl)then
      height = 0.
      return
    endif
    !
    ! First guess: single-layer 6700m atmopshrere and fork around it
    !
    height_up = - H_mono * log(pressure / std_pressure_sl)
    height_down = height_up
    call us_standard_atmosphere(height_up, dens_ratio, press_ratio, tempr_ratio)
    pressure_iter = std_pressure_sl * press_ratio
    pressure_up = pressure_iter
    pressure_down = pressure_iter
    
    do while(pressure_up > pressure)
      height_up = min(height_up * 1.1, max_height)
      call us_standard_atmosphere(height_up, dens_ratio, press_ratio, tempr_ratio)
      pressure_up = std_pressure_sl * press_ratio
    enddo

    do while(pressure_down < pressure)
      height_down = height_down * 0.9
      call us_standard_atmosphere(height_down, dens_ratio, press_ratio, tempr_ratio)
      pressure_down = std_pressure_sl * press_ratio
    enddo
    !
    ! Go into the binary iterations
    !
    do itercount = 1, max_iter
      height = 0.5 * (height_up+height_down)
      call us_standard_atmosphere(height, dens_ratio, press_ratio, tempr_ratio)
      pressure_iter = std_pressure_sl * press_ratio
      if (abs(pressure_iter - pressure) / pressure < rel_tol) exit
      !print *, height_up, height, height_down, pressure, pressure_iter
      if (pressure_iter < pressure) then
        height_up = height
        pressure_up = pressure_iter
      else
        height_down = height
        pressure_down = pressure_iter
      end if
    end do

    if(fu_fails(itercount < max_iter, 'binary search failed', 'find_height'))return
    !
    ! Remaining range should be narrow enough to allow for log approximation
    !
    height = height_down + &
           & (height_up - height_down) * log(pressure_down / pressure) / log(pressure_down / pressure_up)

!call us_standard_atmosphere(height, dens_ratio, press_ratio, tempr_ratio)
!pressure_iter = std_pressure_sl * press_ratio
!call msg('fu_height_for_press. Height, pressure, required pressure, accuracy:', &
!       & (/height, pressure_iter, pressure, (pressure_iter - pressure) / pressure/))
  end function fu_height_for_press

  !************************************************************************************

  real function fu_project_level_from_vertical(vertical_in, fIndex_in, &
                                             & vertical_out, &
                                             & met_data_column, met_data_surf, clip) result(f_ind)
    !
    ! Does exactly the same as fu_project_level_single, actually just calls that sub.
    ! The only difference is that input is not a given level but a given vertical and a relative 
    ! index in it. So, we make the level and then call that sub. Since the index is real,
    ! it is not too trivial: e.g. real index is hybrid vertical is a funny creature.
    !
    implicit none
    ! Imported parameters
    type(silam_vertical), intent(in) :: vertical_in, vertical_out  ! input and output vertical
    real, dimension(:), intent(in), target :: met_data_column ! 3d met data
    real, intent(in) :: fIndex_in, met_data_surf ! index in the input vertical and 2d met data
    logical, intent(in), optional :: clip ! clip mode, see above.

    ! Local variables
    real :: fract
    type(silja_level) :: level
    !
    ! Get the level: find the upper and lower ones in the vertical and make an average
    ! Note that if the vertical is made of thick layers, the index will actually point at the
    ! position inside the specific thick layer.
    !
    if(fu_if_layer(vertical_in%levs(1)))then
      fract = fIndex_in - nint(fIndex_in) + 0.5
    else
      fract = fIndex_in - int(fIndex_in)
    endif

    select case(vertical_in%vert_type)
      ! ---------------------------------------------------------- trivial levels
      case(surface)
        level = fu_set_level(constant_height, fval1=0.0)
      case(mean_sea)
        level = fu_set_level(constant_altitude, fval1=0.0)
      ! ---------------------------------------------------------- thin levels
      case(sigma_level)                                     ! met_data_surf must be surface pressure
        level = fu_set_pressure_level(vertical_in%levs(int(fIndex_in+1))%a * met_data_surf * fract + &
                                    & vertical_in%levs(int(fIndex_in))%a * met_data_surf * (1.-fract))
      case (constant_pressure)
        level = fu_set_pressure_level(vertical_in%levs(int(fIndex_in+1))%a * fract + &
                                    & vertical_in%levs(int(fIndex_in))%a * (1.- fract))
      case (constant_height)
        level = fu_set_constant_height_level(vertical_in%levs(int(fIndex_in+1))%a * fract + &
                                           & vertical_in%levs(int(fIndex_in))%a * (1.- fract))
      case (constant_altitude)
        level = fu_set_constant_altitude_level(vertical_in%levs(int(fIndex_in+1))%a * fract + &
                                             & vertical_in%levs(int(fIndex_in))%a * (1.- fract))
      case (hybrid)                                        ! met_data_surf must be surface pressure
        if(vertical_in%levs(1)%hybrid_coeff_known) then
          call report(vertical_in)
          call set_error('Unknown hybrid coefficients are not allowed', 'fu_project_level_from_vertical')
          return
        else
          level = fu_set_pressure_level( &
                           & (vertical_in%levs(int(fIndex_in+1))%a + &
                            & vertical_in%levs(int(fIndex_in+1))%a * met_data_surf)* fract + &
                           & (vertical_in%levs(int(fIndex_in))%a + &
                            & vertical_in%levs(int(fIndex_in+1))%a * met_data_surf)* (1.- fract))
        end if
      ! ---------------------------------------------------------- thick levels
      case (layer_btw_2_pressure)
        level = fu_set_pressure_level(vertical_in%levs(nint(fIndex_in))%a * fract + &
                                    & vertical_in%levs(nint(fIndex_in))%b * (1. - fract))
      case (layer_btw_2_altitude)
        level = fu_set_constant_altitude_level(vertical_in%levs(nint(fIndex_in))%a * fract + &
                                             & vertical_in%levs(nint(fIndex_in))%b * (1.- fract))
      case (layer_btw_2_height)
        level = fu_set_constant_height_level(vertical_in%levs(nint(fIndex_in))%a * fract + &
                                           & vertical_in%levs(nint(fIndex_in))%b * (1.- fract))
      case (layer_btw_2_sigma)
        level = fu_set_pressure_level(vertical_in%levs(nint(fIndex_in))%a * met_data_surf * fract + &
                                    & vertical_in%levs(nint(fIndex_in))%b * met_data_surf * (1.-fract))
      case (layer_btw_2_hybrid)
        if(vertical_in%levs(1)%hybrid_coeff_known) then
          call report(vertical_in)
          call set_error('Unknown hybrid coefficients are not allowed', 'fu_project_level_from_vertical')
          return
        else
          level = fu_set_pressure_level( &
                           & (vertical_in%levs(nint(fIndex_in))%a + &
                            & vertical_in%levs(nint(fIndex_in))%b * met_data_surf)* fract + &
                           & (vertical_in%levs(nint(fIndex_in))%a2 + &
                            & vertical_in%levs(nint(fIndex_in))%b2 * met_data_surf)* (1.- fract))
        endif
      case default
        call report(vertical_in)
        call set_error('Cannot handle input vertical', 'fu_project_level_from_vertical')
        return
    end select

    if(present(clip))then
      f_ind = fu_project_level_single(level, vertical_out, met_data_column, met_data_surf, clip)
    else
      f_ind = fu_project_level_single(level, vertical_out, met_data_column, met_data_surf)
    endif
    
  end function fu_project_level_from_vertical


  !************************************************************************************

  real function fu_project_level_crude(level_in, vertical, clip) result(f_ind)
    ! 
    ! A convenience function that combines fu_project_level and
    ! vert_interp_data_crude into a single call. Note that if this
    ! needs to be called frequently on the same vertical, it is better
    ! to use the two separately.
    ! 
    implicit none
    type(silja_level), intent(in) :: level_in ! level to project
    type(silam_vertical), intent(in) :: vertical ! vertical to project into
    logical, intent(in), optional :: clip
    
    real, dimension(max_levels) :: met_data_column
    real :: met_data_surf
    type(silam_vertical) :: vert_from_lev_in

    call set_vertical(level_in, vert_from_lev_in)
    call vert_interp_data_crude(vert_from_lev_in, vertical, met_data_column, met_data_surf)
    if ( .not. error) then
       if (present(clip)) then
         f_ind = fu_project_level(level_in, vertical, met_data_column, met_data_surf, clip)
       else
         f_ind = fu_project_level(level_in, vertical, met_data_column, met_data_surf)
       end if
    end if
    
  end function fu_project_level_crude


  !************************************************************************************

  real function fu_project_level_single(level_in, vertical, met_data_column, met_data_surf, clip) &
              & result(f_ind)
    !
    ! Project a level into a vertical. Returns a float index which can
    ! be used eg. for interpolation. Depending on the type of level
    ! and vertical, the projection requires up to one 3d quantity and
    ! one surface quantity, which are given as arguments. 
    !
    ! The above subroutine can be used to get 'standard' values for
    ! the meteo quantities, which provides a crude interpolation.
    ! 
    ! For thin ('full') level verticals, the index is defined so that
    ! the half indices are in the midpoints between levels (either by
    ! pressure or geometric height). If the target vertical consists
    ! of thick layers, the full values are handled similarly, but half
    ! levels are defined by the actual interfaces. 
    !
    ! If the target vertical consists of layers, the length of
    ! met_data_column must be nlevs + 1 with met_data_column(1)
    ! corresponding to the bottom of vertical.
    !
    ! If clip == .true., 1 <= ind <= nlevs for thin level verticals
    ! and 0.5 <= ind <= nlevs+0.5 for layers. This is the default. If
    ! clip == .false., the values are multiplied by -1 if the level is
    ! outside the range (ie. if the level is above the top for a thin
    ! level vertical, f_ind = -real(nlevs), etc.).
    
    implicit none
    type(silja_level), intent(in) :: level_in ! level to project
    type(silam_vertical), intent(in) :: vertical ! vertical to project into
    real, dimension(:), intent(in), target :: met_data_column ! 3d met data
    real, intent(in) :: met_data_surf ! 2d met data
    logical, intent(in), optional :: clip ! clip mode, see above.

    type(silja_level) :: level
    integer :: nlevs
    real :: pressure, clipfactor

    !print *, 'to_proj'

    if (fu_if_layer(level_in)) then
      level = fu_central_level_of_layer(level_in)
    else
      level = level_in
    end if
    nlevs = fu_nbrOflevels(vertical)
    
    if (present(clip)) then
      if (clip) then
        clipfactor = 1.0
      else
        clipfactor = -1.0
      end if
    else
      clipfactor = 1.0
    end if
      
    ! Simplify by reducing the special cases to the more generic types.
    !
    if (level%leveltype == surface)then
      level = fu_set_level(constant_height, fval1=0.0)
    ! The sigma level is hybrid level with a = 0, b = sigma_coef.
    elseif (level%leveltype == sigma_level)then
      level = fu_set_level(hybrid, fval1=0.0, fval2=level%a)
    elseif (level%leveltype == mean_sea)then
      level = fu_set_level(constant_altitude, fval1=0.0)
    elseif (level%leveltype == hybrid .and. .not. level%hybrid_coeff_known) then
      call report(level)
      call set_error('Unknown hybrid coefficients are not allowed', 'fu_project_level_single')
      return
    end if
    select case(level%leveltype)
    case (constant_pressure)
      ! met_data_column is the pressure if needed, met_data_surf is surface_pressure if needed
      f_ind = fu_find_pressure(level%a, vertical, met_data_column, met_data_surf)
    case (constant_height)
      ! met_data_column is the height if needed, met_data_surf is not needed.
      f_ind = fu_find_height(level%a, vertical, met_data_column, met_data_surf)
    case (constant_altitude)
      f_ind = fu_find_height(level%a-met_data_surf, vertical, met_data_column, met_data_surf)
    case (hybrid)
      pressure = level%a + met_data_surf*level%b ! surface pressure always needed for hybrid
      f_ind = fu_find_pressure(pressure, vertical, met_data_column, met_data_surf)
    case default
      call report(level)
      call set_error('Cannot handle input level', 'fu_project_level_single')
    end select
    
  contains
    
    !=================================================================================

    real function fu_find_pressure(pressure, vertical, pressure_column, pressure_surf) result(f_ind)
      ! Find the index for a given height in a vertical, using
      ! pressure_column and pressure_surf if needed.
      implicit none
      real, intent(in) :: pressure
      type(silam_vertical), intent(in) :: vertical
      real, dimension(:), intent(in) :: pressure_column
      real, intent(in) :: pressure_surf

      real level_press, next_level_press, pressure_firstlev, pressure_lastlev, top_press, bottom_press, a, b
      integer :: ilev, leveltype_full
      logical, parameter :: full_levels = .true.
      logical :: is_layer
      type(silja_level) :: top, bottom

      
      leveltype_full = fu_leveltype(fu_level(vertical, 1, full_levels))
      
      !print *, 'find_pressure:', pressure, pressure_column(1:nlevs), nlevs

      select case(fu_leveltype(fu_level(vertical, 1)))
      case (layer_btw_2_pressure)
        top = fu_upper_boundary_of_layer(fu_level(vertical, nlevs))
        bottom = fu_lower_boundary_of_layer(fu_level(vertical, 1))
        if (pressure > fu_pr_level_pressure(bottom)) then
          f_ind = 0.5*clipfactor
          return
        else if (pressure < fu_pr_level_pressure(top)) then
          f_ind = (real(nlevs) + 0.5)*clipfactor
          return
        end if
        do ilev = 1, nlevs
          bottom = fu_lower_boundary_of_layer(fu_level(vertical, ilev))
          top = fu_upper_boundary_of_layer(fu_level(vertical, ilev))
          top_press = fu_pr_level_pressure(top)
          bottom_press = fu_pr_level_pressure(bottom)
          level_press = 0.5*(top_press+bottom_press)
          if (pressure <= bottom_press .and. pressure >= top_press) then
            f_ind = ilev + (level_press - pressure) / (bottom_press - top_press)
            return
          end if
        end do
        call report(vertical, .true.)
        call msg('Failed to find pressure:', pressure)
        call set_error('Failed to find index for pressure:', 'fu_find_pressure')
        
      case (constant_pressure)
        top_press = fu_pr_level_pressure(fu_level(vertical, nlevs))
        bottom_press = fu_pr_level_pressure(fu_level(vertical, 1))
        if (pressure > bottom_press) then          
          f_ind = clipfactor
          return
        else if (pressure < top_press) then
          f_ind = real(nlevs)*clipfactor
          return
        end if

        do ilev = 1, nlevs-1
          level_press = fu_pr_level_pressure(fu_level(vertical, ilev, full_levels))
          next_level_press = fu_pr_level_pressure(fu_level(vertical, ilev+1, full_levels))
          if (pressure <= level_press .and. pressure >= next_level_press) then
            f_ind = ilev + (level_press - pressure) / (level_press - next_level_press)
            return
          end if
        end do
        
        call report(vertical, .true.)
        call msg('Failed to find pressure (in pressure vertical):', pressure)
        call set_error('Failed to find index for pressure:', 'fu_find_pressure')
        
      case (hybrid)
        pressure_firstlev = vertical%levs(1)%a  + vertical%levs(1)%b*pressure_surf

        if (pressure > pressure_firstlev) then
          f_ind = clipfactor
          return
        end if
        
        next_level_press = pressure_firstlev 
        do ilev = 1, nlevs-1
          level_press = next_level_press
          next_level_press = vertical%levs(ilev+1)%a  + vertical%levs(ilev+1)%b*pressure_surf
          if (pressure <= level_press .and. pressure >= next_level_press) then
            f_ind = ilev + (level_press - pressure) / (level_press - next_level_press)
            return
          end if
        end do
        
        f_ind = real(nlevs)*clipfactor

      case (sigma_level)
        pressure_firstlev = vertical%levs(1)%a*pressure_surf
        if (.not. (pressure_firstlev >= 0 .or. pressure_firstlev <= 1e6 )) then
           call report(vertical, .true.)
           call set_error("Gotcha!", "here")
         endif

        if (pressure > pressure_firstlev) then
          f_ind = clipfactor
          return
        end if
        next_level_press = pressure_firstlev 
        do ilev = 1, nlevs-1
          level_press =  next_level_press
          next_level_press =  vertical%levs(ilev+1)%a*pressure_surf
          if (pressure <= level_press .and. pressure >= next_level_press) then
            f_ind = ilev + (level_press - pressure) / (level_press - next_level_press)
            return
          end if
        end do
        
        f_ind = real(nlevs)*clipfactor

      case (layer_btw_2_hybrid, layer_btw_2_sigma)
        top = fu_upper_boundary_of_layer(fu_level(vertical, nlevs))
        bottom = fu_lower_boundary_of_layer(fu_level(vertical, 1))
        if (pressure > fu_hybrid_level_pressure(bottom, pressure_surf)) then
          f_ind = 0.5 * clipfactor
          return
        else if (pressure < fu_hybrid_level_pressure(top, pressure_surf)) then
          f_ind = (real(nlevs) + 0.5) * clipfactor
          return
        end if
        
        do ilev = 1, nlevs
          bottom = fu_lower_boundary_of_layer(fu_level(vertical, ilev))
          top = fu_upper_boundary_of_layer(fu_level(vertical, ilev))
          top_press = fu_hybrid_level_pressure(top, pressure_surf)
          bottom_press = fu_hybrid_level_pressure(bottom, pressure_surf)
          level_press = 0.5*(top_press+bottom_press)

          !print *, pressure, bottom_press, top_press
          if (pressure <= bottom_press .and. pressure >= top_press) then
            f_ind = ilev + (level_press - pressure) / (bottom_press - top_press)
            return
          end if
        end do
        
        call report(vertical, .true.)
        call msg('Failed to find pressure (in hybrid vertical):', pressure)
        call set_error('Failed to find index for pressure', 'fu_find_pressure')


      case (constant_height, constant_altitude)
        ! met_data_column is the level pressure.
        !
        if (pressure > pressure_column(1)) then
          f_ind = clipfactor
          return
        else if (pressure < pressure_column(nlevs)) then
          f_ind = real(nlevs)*clipfactor
          return
        end if
        
        do ilev = 1, nlevs-1
          level_press = pressure_column(ilev)
          next_level_press = pressure_column(ilev+1)
          if (pressure <= level_press .and. pressure >= next_level_press) then
            f_ind = ilev + (level_press - pressure) / (level_press - next_level_press)
            return
          end if
        end do
        
        call report(vertical, .true.)
        call msg('Failed to find pressure:', pressure)
        call set_error('Failed to find index for pressure', 'fu_find_pressure')

      case (layer_btw_2_height, layer_btw_2_altitude)
        ! met_data_column is the pressure at the layer interfaces.
        !
        if (pressure > pressure_column(1)) then
          f_ind = 0.5*clipfactor
          return
        else if (pressure < pressure_column(nlevs+1)) then
          f_ind = (real(nlevs) + 0.5)*clipfactor
          return
        end if
        
        do ilev = 1, nlevs
          top_press = pressure_column(ilev+1)
          bottom_press = pressure_column(ilev)
          level_press = 0.5*(top_press+bottom_press)
          if (pressure <= bottom_press .and. pressure >= top_press) then
            f_ind = ilev + (level_press - pressure) / (bottom_press - top_press)
            return
          end if
        end do
        
        call report(vertical, .true.)
        call msg('Failed to find pressure:', pressure)
        call set_error('Failed to find index for pressure', 'fu_find_pressure')

      case default
        call set_error('Cannot handle this type of vertical', 'fu_find_pressure')
        
      end select

    end function fu_find_pressure

    !==========================================================================
    
    real function fu_find_height(height_in, vertical, height_column, relief_height) result(f_ind)
      ! Find the index for a given height in a vertical, using height_column if needed.
      implicit none
      real, intent(in) :: height_in
      type(silam_vertical), intent(in) :: vertical
      real, dimension(:), intent(in) :: height_column
      real, intent(in) :: relief_height

      real :: height, level_height, next_level_height, bottom_height, top_height
      integer :: ilev, leveltype_full
      logical, parameter :: full_levels = .true.
      type(silja_level) :: top, bottom
      character(len=*), parameter :: sub_name = 'fu_find_height'

      height = height_in
      leveltype_full = fu_leveltype(fu_level(vertical, 1, full_levels)) 
      select case(vertical%vert_type)
      case (layer_btw_2_height, layer_btw_2_altitude)
        ! The case for constant altitude is confusing: now, height is
        ! indeed above ground, but fu_level_height returns level altitudes.
        if (leveltype_full == constant_altitude) height = height - relief_height
        top = fu_upper_boundary_of_layer(fu_level(vertical, nlevs))
        bottom = fu_lower_boundary_of_layer(fu_level(vertical, 1))
        if (height < fu_level_height(bottom)) then
          f_ind = 0.5*clipfactor
          return
        else if (height > fu_level_height(top)) then
          f_ind = (real(nlevs)+0.5)*clipfactor
          return
        end if
        do ilev = 1, nlevs
          top = fu_upper_boundary_of_layer(fu_level(vertical, ilev))
          bottom = fu_lower_boundary_of_layer(fu_level(vertical, ilev))
          top_height = fu_level_height(top)
          bottom_height = fu_level_height(bottom)
          level_height = fu_level_height(fu_level(vertical, ilev, full_levels))
          if (height >= bottom_height .and. height <= top_height) then
            f_ind = ilev + (height-level_height) / (top_height - bottom_height)
            return
          end if
        end do
        call report(vertical, .true.)
        call msg('Failed to find height:', height)
        call set_error('Failed to find index for height in z lyr:', 'fu_find_height')

      case (constant_height, constant_altitude)
        ! The case for constant altitude is confusing: now, height is
        ! indeed above ground, but fu_level_height returns level altitudes.
        if (vertical%vert_type == constant_altitude) height = height - relief_height
        if (height < fu_level_height(fu_level(vertical, 1, full_levels))) then
          f_ind = clipfactor
          return
        else if (height > fu_level_height(fu_level(vertical, nlevs, full_levels))) then
          f_ind = real(nlevs)*clipfactor
          return
        end if

        do ilev = 1, nlevs-1
          level_height = fu_level_height(fu_level(vertical, ilev, full_levels))
          next_level_height = fu_level_height(fu_level(vertical, ilev+1, full_levels))
          if (height >= level_height .and. height <= next_level_height) then
            f_ind = ilev + (height-level_height) / (next_level_height - level_height)
            return
          end if
        end do
        
        call report(vertical, .true.)
        call msg('Failed to find height:', height)
        call set_error('Failed to find index for height in z:', 'fu_find_height')
        
      case (hybrid, sigma_level, constant_pressure)
        ! Now height_column is the level height above ground.
        !
        if (height < height_column(1)) then
          f_ind = clipfactor
          return
        else if (height > height_column(nlevs)) then
          f_ind = real(nlevs)*clipfactor
          return
        end if
                
        do ilev = 1, nlevs-1
          level_height = height_column(ilev)
          next_level_height = height_column(ilev + 1)
          if (height >= level_height .and. height <= next_level_height) then
            f_ind = ilev + (height-level_height) / (next_level_height - level_height)
            return
          end if
        end do
        
        call report(vertical, .true.)
        call msg('Failed to find height:', height)
        call set_error('Failed to find index for height in hyb: ', 'fu_find_height')

      case (layer_btw_2_hybrid, layer_btw_2_sigma, layer_btw_2_pressure)
        !
        ! Now height_column is the level height above ground.
        ! Note thick layers: nLev+1 interfaces
        !
        if (height < height_column(1)) then
          f_ind = 0.5*clipfactor
          return
        else 
          if (height_column(nlevs+1) < height_column(nlevs) .or. height_column(nlevs+1) > 2*height_column(nlevs)) then
            !!! Less-dirty fix for meteolayers
            do ilev = 1, nlevs-1
              level_height = height_column(ilev)
              next_level_height = height_column(ilev + 1)
              if (height >= level_height .and. height <= next_level_height) then
                f_ind = ilev + (height-level_height) / (next_level_height - level_height)
                return
              end if
            end do
#ifdef DEBUG
            call ooops("Wrong height column in " // sub_name)
#endif
            return
          endif
          if (height > height_column(nlevs+1)) then
            f_ind = real(nlevs+0.5)*clipfactor
            return
          end if
        endif
                
        do ilev = 1, nlevs
          top_height = height_column(ilev+1)
          bottom_height = height_column(ilev)
          level_height = 0.5 * (top_height + bottom_height)
          if (height >= bottom_height .and. height <= top_height) then
            f_ind = ilev + (height-level_height) / (top_height - bottom_height)
            return
          end if
        end do
        
        call report(vertical, .true.)
        call msg('Failed to find height:', height)
        call set_error('Failed to find index for height in hybrid layers:', 'fu_find_height')

      case (surface)
        if (height == 0.) then 
          f_ind = 1
        else
          call report(vertical, .true.)
          call msg('Failed to find height in this vertical:', height)
          call set_error('Failed to find index for height in suface vert:', 'fu_find_height')
        endif


      case default
        call set_error('Cannot handle this type of vertical', 'fu_find_height')
        
      end select

    end function fu_find_height

  end function fu_project_level_single
  
  !************************************************************************************
  
  subroutine projection_input_needs(vert_from, vert_to, quantity_3d, quantity_2d)
    !
    ! Return input needs for the meteo-dependent projection of
    ! vertical levels. Depending on the type verticals, up to one 3d
    ! and one 2d quantity is required.
    !
    implicit none
    type(silam_vertical), intent(in) :: vert_to, vert_from
    integer, intent(out) :: quantity_3d, quantity_2d
    
    type(silja_level) :: level
    integer :: leveltype_from, leveltype_to

    ! The same quantities are required for thick and thin layers and
    ! we'll just see the type of the midpoints.
    level = fu_level(vert_from, 1, .true.)
    leveltype_from = level%leveltype
    level = fu_level(vert_to, 1, .true.)
    leveltype_to = level%leveltype

    select case (leveltype_from)
    case (constant_height)
      select case(leveltype_to)
      case (constant_height, surface)
        quantity_3d = int_missing
        quantity_2d = int_missing
      case (constant_altitude, mean_sea)
        quantity_3d = int_missing
        quantity_2d = relief_height_flag
      case (constant_pressure, hybrid, sigma_level)
        quantity_3d = height_flag
        quantity_2d = int_missing
      case default
        call msg('vert_from:')
        call report(vert_from)
        call msg('vert_to:')
        call report(vert_to)
        call set_error('Cannot project verticals', 'fu_vertInterpInputNeeds')
      end select

    case (constant_altitude)
      select case(leveltype_to)
      case (constant_height, surface)
        quantity_3d = int_missing
        quantity_2d = relief_height_flag
      case (constant_altitude, mean_sea)
        quantity_3d = int_missing
        quantity_2d = int_missing
      case (constant_pressure, hybrid, sigma_level)
        quantity_3d = height_flag
        quantity_2d = relief_height_flag
      case default
        call msg('vert_from')
        call report(vert_from)
        call msg('vert_to')
        call report(vert_to)
        call set_error('Cannot project verticals', 'fu_vertInterpInputNeeds')
      end select

    case (constant_pressure)
      select case(leveltype_to)
      case (constant_pressure)
        quantity_3d = int_missing
        quantity_2d = int_missing
      case (constant_altitude, mean_sea)
        quantity_3d = pressure_flag
        quantity_2d = relief_height_flag
      case (constant_height, surface)
        quantity_3d = pressure_flag
        quantity_2d = int_missing
      case (hybrid, sigma_level)
        quantity_3d = int_missing
        quantity_2d = surface_pressure_flag
      case default
        call msg('vert_from')
        call report(vert_from)
        call msg('vert_to')
        call report(vert_to)
        call set_error('Cannot project verticals', 'fu_vertInterpInputNeeds')
      end select

    case (sigma_level, hybrid)
      select case(leveltype_to)
      case (constant_pressure)
        quantity_3d = int_missing
        quantity_2d = surface_pressure_flag
      case (constant_altitude, mean_sea)
        quantity_3d = pressure_flag
        quantity_2d = surface_pressure_flag
      case (constant_height, surface)
        quantity_3d = pressure_flag
        quantity_2d = surface_pressure_flag
      case (hybrid, sigma_level)
        ! Projecting between hybrid levels always requires surface
        ! pressure. Between pure sigma levels, we could do without, but
        ! since surface pressure is usually available, we just handle
        ! sigma as a special case of hybrid.
        quantity_3d = int_missing
        quantity_2d = surface_pressure_flag
      case default
        call msg('vert_from')
        call report(vert_from)
        call msg('vert_to')
        call report(vert_to)
        call set_error('Cannot project verticals', 'fu_vertInterpInputNeeds')
      end select
    end select
    
  end subroutine projection_input_needs

  !************************************************************************************
  
  real function fu_vert_overlap_fraction(layer_fixed, layer_with, &
                                       & overlap_centre_idx_lyr_fixed, overlap_centre_idx_lyr_with) &
                                       & result(fraction)
    !
    ! Compute the fraction of layer_fixed covered by layer_with. The
    ! idea is that the sum of fractions is one when the first layer is
    ! held fixed and projected into levels of a vertical.
    !
    ! Algorithm: project the boundaries of layer_with into a vertical
    ! formed by layer_fixed.
    ! Note that for the overlap centre in the layer_with we have to do the opposite reprojection
    !
    implicit none
    type(silja_level), intent(in) :: layer_fixed, layer_with
    real, intent(out), optional :: overlap_centre_idx_lyr_fixed, overlap_centre_idx_lyr_with
    
    type(silam_vertical) :: vert_from_layer_fixed, vert_from_layer_with
    real, dimension(:), pointer :: met_data_column
    real :: met_data_surf, f_ind_bottom, f_ind_top
    type(silja_level) :: bottom, top

    if (fu_fails(fu_if_layer(layer_fixed), 'layer_fixed not a thick layer', 'fu_vert_verlap_fraction')) return
    if (fu_fails(fu_if_layer(layer_with), 'layer_with not a thick layer', 'fu_vert_verlap_fraction')) return

    call set_vertical(layer_fixed, vert_from_layer_fixed)
    call set_vertical(layer_with, vert_from_layer_with)

    met_data_column => fu_work_array()
    call vert_interp_data_crude(vert_from_layer_with, vert_from_layer_fixed, &
                              & met_data_column, met_data_surf)
    if (error) then
      call cleanup()
      return
    end if

    bottom = fu_lower_boundary_of_layer(layer_with)
    f_ind_bottom = fu_project_level(bottom, vert_from_layer_fixed, met_data_column, met_data_surf)
    top = fu_upper_boundary_of_layer(layer_with)
    f_ind_top = fu_project_level(top, vert_from_layer_fixed, met_data_column, met_data_surf)
    fraction = f_ind_top - f_ind_bottom
    !
    ! Indices of the overlap centre points in fixed and with- layers
    !
    if(present(overlap_centre_idx_lyr_fixed)) &
                                    & overlap_centre_idx_lyr_fixed = (f_ind_top + f_ind_bottom) / 2.
    !
    ! To get it the overlap centre in the with-layer, have to do the opposite projection of the 
    !
    if(present(overlap_centre_idx_lyr_with))then
      if(fraction .eps. 0.0)then
        overlap_centre_idx_lyr_with = real_missing
        call cleanup()
        return
      endif
      call vert_interp_data_crude(vert_from_layer_fixed, vert_from_layer_with, &
                                & met_data_column, met_data_surf)
      if (error) then
        call cleanup()
        return
      end if
      bottom = fu_lower_boundary_of_layer(layer_fixed)
      f_ind_bottom = fu_project_level(bottom, vert_from_layer_with, met_data_column, met_data_surf)
      top = fu_upper_boundary_of_layer(layer_fixed)
      f_ind_top = fu_project_level(top, vert_from_layer_with, met_data_column, met_data_surf)
      overlap_centre_idx_lyr_with = (f_ind_top + f_ind_bottom) / 2.
    endif
    call cleanup()

  contains
    subroutine cleanup()
      implicit none
      call free_work_array(met_data_column)
      call set_missing_vertical(vert_from_layer_fixed, .false.)
      call set_missing_vertical(vert_from_layer_with, .false.)
    end subroutine cleanup
    
  end function fu_vert_overlap_fraction

  
  !************************************************************************************

  subroutine overlap_fraction_lyr_in_vert(vert_fixed, vert_with, &
                                        & fractions, overlap_centre_point_index, &
                                        & met_data_column_fixed2with, met_data_srf_fixed2with, &
                                        & met_data_column_with2fixed, met_data_srf_with2fixed)
    !
    ! Computes the fraction of layers of vert_fixed covered by layers of vert_with. The
    ! sum of fractions along that vert_with is 1 for each layer of vert_with (if it has 
    ! sufficient span to cover the layer, of course).
    !
    ! Algorithm: project the boundaries of layer_with into a vertical formed by layer_fixed.
    ! We scan all the layers in both verticals. Since each layer_fix is affected by only
    ! few layers_with, we first do inverse projection to identify the important layers. 
    ! The others can be skipped.
    !
    ! ATTENTION.
    ! Has to be very quick since called for interpolation structures inside the grid cycle. 
    ! So, ABSOLUTELY NO CHECKING.
    !
    implicit none
    
    ! Imported parameters
    type(silam_vertical), intent(in) :: vert_fixed, vert_with
    real, dimension(:,:), intent(out) :: fractions, overlap_centre_point_index
    real, dimension(:), target, optional :: met_data_column_fixed2with, &
         & met_data_column_with2fixed
    real, intent(in), optional :: met_data_srf_fixed2with, met_data_srf_with2fixed
    
    ! Local variables
    real, dimension(:), pointer :: met_data_column_fix2with_loc, met_data_column_with2fix_loc
    real :: met_data_srf_fix2with_loc, met_data_srf_with2fix_loc, f_ind_bottom, f_ind_top
    type(silja_level) :: bottom, top
    integer :: iFixed, iWith, iWithStart, iWithEnd
    !
    ! Get the meteo data - either from the input or from default constants
    ! Note that, since meteo buffers are all for central points, we cannot
    ! handle layers here for default meteodata either. Fake vert_fixed__levels
    ! and vert_with__levels represent the corresponding verticals - but with levels.
    ! Useful only for reprojection.
    !
    if(present(met_data_column_fixed2with))then
      met_data_column_fix2with_loc => met_data_column_fixed2with
      met_data_column_with2fix_loc => met_data_column_with2fixed
      met_data_srf_fix2with_loc = met_data_srf_fixed2with
      met_data_srf_with2fix_loc = met_data_srf_with2fixed
    else
      met_data_column_fix2with_loc => fu_work_array()
      met_data_column_with2fix_loc => fu_work_array()
      call vert_interp_data_crude(vert_with, vert_fixed, &
                                & met_data_column_with2fix_loc, met_data_srf_with2fix_loc)
      call vert_interp_data_crude(vert_fixed, vert_with, &
                                & met_data_column_fix2with_loc, met_data_srf_fix2with_loc)
    endif
    if(error)return
    fractions(1:vert_fixed%nLevs, 1:vert_with%nLevs) = 0.0
    overlap_centre_point_index(1:vert_fixed%nLevs, 1:vert_with%nLevs) = 0.0
    !
    ! Scan all levels of vert_fixed, for each of them get the fractions covered by
    ! levels of vertical_with
    !
    do iFixed = 1, vert_fixed%nLevs
      !
      ! get the range of layers in vert_with that cover this layer: project this
      ! fixed layer to vert_with
      !
      bottom = fu_lower_boundary_of_layer(vert_fixed%levs(iFixed))
      top = fu_upper_boundary_of_layer(vert_fixed%levs(iFixed))
      f_ind_bottom = fu_project_level(bottom, vert_with, &
                                    & met_data_column_fix2with_loc, met_data_srf_fix2with_loc)
      f_ind_top = fu_project_level(top, vert_with, &
                                 & met_data_column_fix2with_loc, met_data_srf_fix2with_loc)
      iWithStart = max(1,int(f_ind_bottom))
      iWithEnd = min(vert_with%nLevs,int(f_ind_top)+1)
      !
      ! Scan the _with-range getting the overlaps
      !
      do iWith = iWithStart, iWithEnd
        bottom = fu_lower_boundary_of_layer(vert_with%levs(iWith))
        top = fu_upper_boundary_of_layer(vert_with%levs(iWith))
        f_ind_bottom = fu_project_level(bottom, vert_fixed, &
                                      & met_data_column_with2fix_loc, met_data_srf_with2fix_loc)
        f_ind_top = fu_project_level(top, vert_fixed, &
                                   & met_data_column_with2fix_loc, met_data_srf_with2fix_loc)
        fractions(iFixed,iWith) = f_ind_top - f_ind_bottom
        overlap_centre_point_index(iFixed,iWith) = (f_ind_top + f_ind_bottom) / 2.
      end do  ! iWith
    end do ! iFixed

    if(.not. present(met_data_column_fixed2with))then
      call free_work_array(met_data_column_fix2with_loc)
      call free_work_array(met_data_column_with2fix_loc)
    endif

  end subroutine overlap_fraction_lyr_in_vert


  !*****************************************************************************************

  subroutine vertical_parameters_ab(vertical, nz, a_full, b_full)
    !
    ! Computes the parameters of vertical and stores to the public parameters
    !
    implicit none
    
    ! Imported parameters
    type(silam_vertical), intent(in) :: vertical
    integer, intent(out) :: nz
    real, dimension(:), pointer :: a_full, b_full

    ! Local variables
    integer :: iTmp

    nz = fu_nbrOflevels(vertical)
    
    if (.not. associated(a_full))then 
      allocate(a_full(0:nz+1), b_full(0:nz+1), stat=iTmp)
      if (fu_fails(iTmp == 0, 'Allocate failed a_full  b_full', 'initEulerAdvectionFields'))return
    else
      if(size(a_full) /= nz+2)then
        deallocate(a_full, b_full)
        allocate(a_full(0:nz+1), b_full(0:nz+1), stat=iTmp)
        if (fu_fails(iTmp == 0, 'Allocate failed a_full  b_full', 'initEulerAdvectionFields'))return
      endif
    endif
    
    call hybrid_coefs(vertical, a_full=a_full(1:nz), b_full=b_full(1:nz))
    a_full(0) = 0. 
    b_full(0) = 1. !Surface
    a_full(nz+1) = 0. ! Pa
    b_full(nz+1) = 0. 
    
  end subroutine vertical_parameters_ab
  
  
  !*****************************************************************************************

  subroutine vertical_parameters(vertical, nz, a_half, b_half, layer_top_m, ifAllocate, ifExtend)
    !
    ! Computes the parameters of vertical of layers and stores to the (newly allocated) arrays
    ! Note that indexing of the pointers should be (-1:nz+1)
    !
    implicit none
    
    ! Imported parameters
    type(silam_vertical), intent(in) :: vertical
    logical, intent(in) :: ifAllocate, ifExtend
    integer, intent(out) :: nz
    real, dimension(:), pointer :: a_half, b_half, layer_top_m



    ! Local variables
    integer :: iTmp, iLev

    nz = fu_nbrOflevels(vertical)

    select case (fu_leveltype(vertical))

        case(layer_btw_2_height)
          if (ifAllocate) then 
              allocate(layer_top_m(-1:nz+1), stat=iTmp) ! Including underground and
                                        ! above-atmosphere
              if (iTmp /= 0) call set_error('Allocate failed layer_top_m', 'vertical_parameters')
              nullify(a_half,b_half)
          endif

          do iLev = 1, nz
            layer_top_m(iLev) = fu_top_of_layer_value(vertical, iLev)
          end do
          if (ifExtend) then
             ! The over-the-domain cover is faked by repeating twice the whole domain. Advection will 
             ! be never computed for it
             layer_top_m(-1) = -100000. ! -100 km  underground
             layer_top_m(0) = fu_bottom_of_layer_value(vertical, 1)
             layer_top_m(nz+1) = max(0., 3.* layer_top_m(nz) - 2.* layer_top_m(0))
           endif

        case (layer_btw_2_hybrid, layer_btw_2_pressure, layer_btw_2_sigma)
!          call msg('Preparing eta-level vertical advection')
          if (ifAllocate) then 
              allocate(a_half(0:nz+2), b_half(0:nz+2), stat=iTmp)
              if (iTmp /= 0) call set_error('Allocate failed a_half  b_half', 'vertical_parameters')
              nullify(layer_top_m)
           endif
          
          call hybrid_coefs(vertical, a_half=a_half(1:nz+1), b_half=b_half(1:nz+1))
          if (ifExtend) then
             a_half(0) = 0. 
             b_half(0) = 2. !whole atmosphere below
             a_half(nz+2) = 1. !Pa
                           !interstellar space above, still larger than meteo
             b_half(nz+2) = 0. 
          endif
          case default
          call set_error('Unknown dispersion vertical','vertical_parameters')
          return
    end select
    
  end subroutine vertical_parameters


  !***********************************************************************
  
  subroutine make_vertical_of_levels(vertIn, vertOfLevels)
    !
    ! Copies the input vertical, replacing each layer with its centre-point level
    !
    implicit none
    
    ! IMported parameters
    type(silam_vertical), intent(in) :: vertIn
    type(silam_vertical), intent(out) :: vertOfLevels

    ! Local variables
    integer :: iLev

    if(fu_if_layer(vertIn%levs(1)))then
      call set_vertical(fu_level(vertIn, 1, .true.), vertOfLevels)
      do iLev = 2, vertIn%nLevs
        call add_level(vertOfLevels, fu_level(vertIn, iLev, .true.))
      end do
    else
      vertOfLevels = vertIn
    endif

  end subroutine make_vertical_of_levels

  !***********************************************************************
  
  recursive subroutine make_vertical_of_borders(vertIn, vertOfBorders, ifConvertLevels2Layers)
    !
    ! Copies the input vertical, replacing each layer with its centre-point level
    !
    implicit none
    
    ! IMported parameters
    type(silam_vertical), intent(in) :: vertIn
    type(silam_vertical), intent(out) :: vertOfBorders
    logical, intent(in) :: ifConvertLevels2Layers

    ! Local variables
    integer :: iLev
    type(silam_vertical) :: vertTmp

    if(fu_if_layer(vertIn%levs(1)))then
      call set_vertical(fu_lower_boundary_of_layer(fu_level(vertIn, 1)), vertOfBorders)
      do iLev = 1, vertIn%nLevs
        call add_level(vertOfBorders, fu_upper_boundary_of_layer(fu_level(vertIn, iLev)))
      end do
    else
      if(ifConvertLevels2Layers)then
        call make_vertical_of_layers(vertIn, vertTmp)
        call make_vertical_of_borders(vertTmp, vertOfBorders, .false.)
      else
        call msg_warning('Original vertical of levels is returned instead of borders', &
                       & 'make_vertical_of_borders')
        vertOfBorders = vertIn
      endif
    endif

  end subroutine make_vertical_of_borders


  !***********************************************************************
  
  subroutine make_vertical_of_layers(vertIn, vertOfLayers)
    !
    ! Copies the input vertical, replacing each thin level with layer. Borders are to be
    ! generated...
    !
    implicit none
    
    ! IMported parameters
    type(silam_vertical), intent(in) :: vertIn
    type(silam_vertical), intent(out) :: vertOfLayers

    ! Local variables
    integer :: iLev

    if(fu_if_layer(vertIn%levs(1)))then
      vertOfLayers = vertIn
    else
      select case(vertIn%vert_type)
        case(constant_altitude, constant_height, depth_level)
          call turn_general_levels_2_layers(vertIn, vertOfLayers, 0.0, 1.0)
        case(constant_pressure, sigma_level)
          call turn_general_levels_2_layers(vertIn, vertOfLayers, std_pressure_sl, -1.0)
        case(hybrid)
          call turn_hybrid_levels_2_layers(vertIn, vertOfLayers)
        case default
          call report(vertIn)
          call set_error('Not a supported level type','make_vertical_of_layers')
      end select

    endif  ! if the input vertical is of layers

    CONTAINS
    
    !===========================================================
    integer function fu_layer_type_for_level(level_type)
      implicit none
      integer, intent(in) :: level_type
      select case(level_type)
        case(constant_altitude)
          fu_layer_type_for_level = layer_btw_2_altitude
        case(constant_height)
          fu_layer_type_for_level = layer_btw_2_height
        case(depth_level)
          fu_layer_type_for_level = layer_btw_2_depth
        case(constant_pressure)
          fu_layer_type_for_level = layer_btw_2_pressure
        case(sigma_level)
          fu_layer_type_for_level = layer_btw_2_sigma
        case(hybrid)
          fu_layer_type_for_level = layer_btw_2_hybrid
        case default
          call set_error('Not a supported level' + fu_str(level_type),'fu_layer_type_for_level')
          fu_layer_type_for_level = int_missing
      end select
    end function fu_layer_type_for_level

    !==========================================================================
    
    subroutine turn_general_levels_2_layers(vertIn, vertOfLayers, fLowBorder, fSign)
      implicit none
      type(silam_vertical), intent(in) :: vertIn
      type(silam_vertical), intent(out) :: vertOfLayers
      real, intent(in) :: fLowBorder, fSign

      integer :: iLayerType, iLev
      real :: fBorderUp, fBorderDown
      logical :: ifOK

      iLayerType = fu_layer_type_for_level(vertIn%vert_type)
      fBorderDown = fLowBorder
      !
      ! One level?
      !
      if(vertIn%nLevs < 2)then
        call set_vertical(fu_set_layer_between_two(iLayerType, vertIn%levs(1)%a * 2., fBorderDown), &
                        & vertOfLayers)
        return
      endif
      !
      ! Main work starts
      !
      ifOK = .true.
      fBorderUp = vertIn%levs(1)%a * 2.
      if(fBorderUp * fSign > vertIn%levs(2)%a * fSign)then
        ifOK = .false.
        fBorderUp = (vertIn%levs(1)%a + vertIn%levs(2)%a) / 2.
      else
        ifOK = .true.
      endif
      call set_vertical(fu_set_layer_between_two(iLayerType, fBorderUp, fBorderDown), vertOfLayers)

      if(ifOK)then
        do iLev = 2, vertIn%nLevs
          if(fBorderUp * fSign > vertIn%levs(iLev)%a * fSign)then
            ifOK = .false.
            exit
          endif
          fBorderDown = fBorderUp
          fBorderUp = 2. * vertIn%levs(iLev)%a - fBorderDown
          call add_level(vertOfLayers, fu_set_layer_between_two(iLayerType, fBorderUp, fBorderDown))
        end do
      endif

      if(.not. ifOK)then
        call msg_warning('Cannot create layers rigorously','turn_general_levels_2_layers')
        do iLev = 2, vertIn%nLevs-1
          fBorderDown = fBorderUp
          fBorderUp = (vertIn%levs(iLev)%a + vertIn%levs(iLev+1)%a) / 2.
          call add_level(vertOfLayers, fu_set_layer_between_two(iLayerType, fBorderUp, fBorderDown))
        end do
        call add_level(vertOfLayers, &
                     & fu_set_layer_between_two(iLayerType, &
                                              & vertIn%levs(vertIn%nLevs)%a * 2. - fBorderUp, &
                                              & fBorderUp))
      endif

    end subroutine turn_general_levels_2_layers

    !===============================================================================

    subroutine turn_hybrid_levels_2_layers(vertIn, vertOfLayers)
      implicit none
      type(silam_vertical), intent(in) :: vertIn
      type(silam_vertical), intent(out) :: vertOfLayers

      integer :: iLev
      real :: fBorderUp_a, fBorderUp_b, fBorderDown_a, fBorderDown_b

      fBorderDown_a = 0.
      fBorderDown_b = 1.
      !
      ! One level?
      !
      if(vertIn%nLevs < 2)then
        call set_vertical(fu_set_layer_between_two(layer_btw_2_hybrid, 1, 1, &
                                             & vertIn%levs(1)%a * 2., (1. - vertIn%levs(1)%b * 2.), &
                                             & fBorderDown_a, fBorderDown_b), &
                        & vertOfLayers)
        return
      endif
      !
      ! Main work starts. Coefs for half-levels are half-sums of the full-level values.
      !
      fBorderUp_a = (vertIn%levs(1)%a + vertIn%levs(2)%a) / 2.
      fBorderUp_b = (vertIn%levs(1)%b + vertIn%levs(2)%b) / 2.
      call set_vertical(fu_set_layer_between_two(layer_btw_2_hybrid, 1, 1, fBorderUp_a, fBorderUp_b, &
                                                                 & fBorderDown_a, fBorderDown_b), &
                      & vertOfLayers)
      do iLev = 2, vertIn%nLevs - 1
        fBorderDown_a = fBorderUp_a
        fBorderDown_b = fBorderUp_b
        fBorderUp_a = (vertIn%levs(iLev)%a + vertIn%levs(iLev+1)%a) / 2.
        fBorderUp_b = (vertIn%levs(iLev)%b + vertIn%levs(iLev+1)%b) / 2.
        call add_level(vertOfLayers, fu_set_layer_between_two(layer_btw_2_hybrid, iLev, iLev, &
                                                            & fBorderUp_a, fBorderUp_b, &
                                                            & fBorderDown_a, fBorderDown_b))
      end do

      call add_level(vertOfLayers, &
                   & fu_set_layer_between_two(layer_btw_2_hybrid, vertIn%nLevs, vertIn%nLevs, &
                                            & max(0., vertIn%levs(vertIn%nLevs)%a * 2. - fBorderUp_a), &
                                            & max(0., vertIn%levs(vertIn%nLevs)%b * 2. - fBorderUp_b), &
                                            & fBorderUp_a, fBorderUp_b))

    end subroutine turn_hybrid_levels_2_layers

  end subroutine make_vertical_of_layers


  !***********************************************************************

  recursive function fu_level_to_vertical(levIn, vert, ptrMeteo) result(levelOut)
    !
    ! Projects level to some vertical structure.
    ! If level type and verttical type are the same - nothing happens,
    ! the same level is returned. But if they are not the same, the
    ! level is interpolated to the vetical. Evidently, this may require
    ! knowledge about surface pressure and temperature, etc, which is
    ! provided via the ptrMeteo. 
    !
    ! Since this routine will be called inside grid cycle, it has to be made
    ! quick and with minimum checking
    !
    implicit none

    ! return value
    type(silja_level) :: levelOut

    ! Imported parameters
    type(silja_level), intent(in) :: levIn
    type(silam_vertical), intent(in)  :: vert
    real, dimension(:), pointer :: ptrMeteo

    ! Local variables
    type(silja_level) :: levTmp
    real :: pr1, pr2
    integer :: iLev

    !
    ! A speedup attempt
    !
    if(vert%vert_type == levIn%leveltype)then
      levelOut = levIn
      return
    else
      levelOut = level_missing
    endif

    if(fu_if_layer(levIn))then
      levTmp = fu_central_level_of_layer(levIn)
    else
      levTmp = levIn
    endif

    ! How about now?
    !
    if(vert%vert_type == levTmp%leveltype)then
      levelOut = levTmp
      return
    else
      levelOut = level_missing
    endif

    ! Have to project each-to-each
    !
    select case(vert%vert_type)
      case(surface) ! Require the levels to be close to the surface
        select case(levTmp%leveltype)
          case(mean_sea)
            levelOut = surface_level
          case(constant_pressure)
            if(levTmp%a >= 97500.) levelOut = surface_level
          case(constant_altitude)
            if(levTmp%a <= 50.) levelOut = surface_level
          case(constant_height)
            if(levTmp%a <= 50.) levelOut = surface_level
          case(sigma_level)
            if(levTmp%a > 0.975) levelOut = surface_level
          case(hybrid)
            if(levTmp%number == 1) levelOut = surface_level
          case default
            call report(vert)
            call report(levTmp)
            call set_error('Can not project above level to above vertical', &
                         & 'fu_level_to_vertical')
        end select

      case(mean_sea)
        select case(levTmp%leveltype)
          case(surface)
            levelOut = mean_sea_level
          case(constant_pressure)
            if(levTmp%a >= 97500.) levelOut = mean_sea_level
          case(constant_altitude)
            if(levTmp%a <= 50.) levelOut = mean_sea_level
          case(constant_height)
            if(levTmp%a <= 50.) levelOut = mean_sea_level
          case(sigma_level)
            if(levTmp%a > 0.975) levelOut = mean_sea_level
          case(hybrid)
            if(levTmp%number == 1) levelOut = mean_sea_level
          case default
            call report(vert)
            call report(levTmp)
            call set_error('Can not project above level to above vertical', &
                         & 'level_to_vertical')
        end select

      case(constant_altitude)
        select case(levTmp%leveltype)
          case(surface, mean_sea)
            levelOut = fu_set_constant_altitude_level(0.)
          case(constant_pressure)
            levelOut = fu_set_constant_altitude_level(ptrMeteo(indexRelief) + &
                               & gas_constant_dryair * ptrMeteo(indexTemperature) / g * &
                               & log(ptrMeteo(indexSurfPressure) / &
                                   & fu_pr_level_pressure(levTmp)))
          case(constant_height)
            levelOut = fu_set_constant_altitude_level(ptrMeteo(indexRelief) + &
                                                    & fu_level_height(levTmp))
          case(sigma_level)
            levelOut = fu_set_constant_altitude_level(ptrMeteo(indexRelief) + &
                               & gas_constant_dryair * ptrMeteo(indexTemperature) / g * &
                               & (-log(fu_sigma_level_sigma(levTmp))))
          case(hybrid)
            levelOut = fu_set_constant_altitude_level(ptrMeteo(indexRelief) + &
                               & gas_constant_dryair * ptrMeteo(indexTemperature) * &
                               & log(ptrMeteo(indexSurfPressure) / &
                                   & fu_hybrid_level_pressure(levTmp, &
                                                            & ptrMeteo(indexSurfPressure))))
          case(depth_level)
            levelOut = fu_set_constant_altitude_level(ptrMeteo(indexRelief) - &
                                                    & fu_depth_level_depth(levTmp))
          case default
            call report(vert)
            call report(levTmp)
            call set_error('Can not project above level to above vertical', &
                         & 'level_to_vertical')
        end select

      case(constant_height)
        select case(levTmp%leveltype)
          case(surface, mean_sea)
            levelOut = fu_set_constant_height_level(0.)
          case(constant_pressure)
            levelOut = fu_set_constant_height_level(&
                               & gas_constant_dryair * ptrMeteo(indexTemperature) / g * &
                               & log(ptrMeteo(indexSurfPressure) / fu_pr_level_pressure(levTmp)))
          case(constant_altitude)
            levelOut = fu_set_constant_height_level(fu_level_altitude(levTmp) - ptrMeteo(indexRelief))
          case(sigma_level)
            levelOut = fu_set_constant_height_level(&
                               & gas_constant_dryair * ptrMeteo(indexTemperature) / g * &
                               & (- log(fu_sigma_level_sigma(levTmp))))
          case(hybrid)
            levelOut = fu_set_constant_height_level(&
                               & gas_constant_dryair * ptrMeteo(indexTemperature) / g * &
                               & log(ptrMeteo(indexSurfPressure) / &
                                   & fu_hybrid_level_pressure(levTmp, ptrMeteo(indexSurfPressure))))
          case(depth_level)
            levelOut = fu_set_constant_height_level(-fu_depth_level_depth(levTmp))
          case default
            call report(vert)
            call report(levTmp)
            call set_error('Can not project above level to above vertical', &
                         & 'fu_level_to_vertical')
        end select

      case(constant_pressure)
        select case(levTmp%leveltype)
          case(surface) 
            levelOut = fu_set_pressure_level(ptrMeteo(indexSurfPressure))
          case(mean_sea)
            levelOut = fu_set_pressure_level(ptrMeteo(indexSurfPressure) / &
                     & exp(-ptrMeteo(indexRelief) / &
                         & (gas_constant_dryair*ptrMeteo(indexTemperature)/g)))
          case(constant_altitude)
            levelOut = fu_set_pressure_level(ptrMeteo(indexSurfPressure) / &
                     & exp((fu_level_altitude(levTmp) - ptrMeteo(indexRelief)) / &
                         & (gas_constant_dryair * ptrMeteo(indexTemperature) / g)))
          case(constant_height)
            levelOut = fu_set_pressure_level(ptrMeteo(indexSurfPressure) / &
                     & exp(fu_level_height(levTmp) / &
                         & (gas_constant_dryair * ptrMeteo(indexTemperature) / g)))
          case(sigma_level)
            levelOut = fu_set_pressure_level(ptrMeteo(indexSurfPressure) * &
                                           & fu_sigma_level_sigma(levTmp))
          case(hybrid)
            levelOut = fu_set_pressure_level(fu_hybrid_level_pressure(levTmp, &
                                                         & ptrMeteo(indexSurfPressure)))
          case default
            call report(vert)
            call report(levTmp)
            call set_error('Can not project above level to above vertical', &
                         & 'fu_level_to_vertical')
        end select

      case(sigma_level)
        select case(levTmp%leveltype)
          case(surface)
            levelOut = fu_set_sigma_level(1.)
          case(mean_sea)
            levelOut = fu_set_sigma_level(1. / &
                     & exp(-ptrMeteo(indexRelief) / &
                         & (gas_constant_dryair * ptrMeteo(indexTemperature) / g)))
          case(constant_altitude)
            levelOut = fu_set_sigma_level(1. / &
                     & exp((fu_level_altitude(levTmp) - ptrMeteo(indexRelief)) / &
                         & (gas_constant_dryair * ptrMeteo(indexTemperature) / g)))
          case(constant_height)
            levelOut = fu_set_sigma_level(1. / &
                     & exp(fu_level_height(levTmp) / &
                         & (gas_constant_dryair * ptrMeteo(indexTemperature) / g)))
          case(constant_pressure)
            levelOut = fu_set_sigma_level(fu_pr_level_pressure(levTmp) / ptrMeteo(indexSurfPressure))
          case(hybrid)
            levelOut = fu_set_sigma_level(fu_hybrid_level_pressure(levTmp, &
                                                                 & ptrMeteo(indexSurfPressure)) / &
                                                                 & ptrMeteo(indexSurfPressure))
          case default
            call report(vert)
            call report(levTmp)
            call set_error('Can not project above level to above vertical', &
                         & 'fu_level_to_vertical')
        end select

      case(hybrid)
        !
        ! For hybrid levels we have to follow two steps: make a pressure level and then 
        ! project it to hybrid one
        !
        if(error)return
        select case(levTmp%leveltype)
          case(surface)
            levelOut = fu_set_hybrid_level(int_missing, 0.0, 1.0)  ! The surface pressure: a+b*surfPr
          case(mean_sea)
            levelOut = fu_set_hybrid_level(int_missing, 0.0, 1.0 / &
                                         & exp(-ptrMeteo(indexRelief) / &
                                             & (gas_constant_dryair*ptrMeteo(indexTemperature)/g)))

          case(constant_altitude, constant_height, sigma_level, hybrid)
            levelOut = fu_level_to_vertical(fu_level_to_vertical(levTmp, vertical_1000hpa, ptrMeteo), &
                                          & vert, &
                                          & ptrMeteo)
          case(constant_pressure)
            !
            ! Have to scan the vertical looking for two hybrid levels surrounding levelIn
            ! and then interpolate the a,b coefficients
            ! Note: hybrid_level_pressure = level%a + ground_surface_pressure*level%b
            !
            pr1 = vert%levs(1)%a + vert%levs(1)%b * ptrMeteo(indexSurfPressure)
            do iLev = 2, vert%nLevs
              pr2 = vert%levs(iLev)%a + vert%levs(iLev)%b * ptrMeteo(indexSurfPressure)
              if((levTmp%a - pr1)*(levTmp%a - pr2) <= 0.0)then
                pr1 = (levTmp%a - pr1) / (pr2 - pr1)          ! Temporary use of the variable
                levelOut = fu_set_hybrid_level(int_missing, &
                                  & vert%levs(iLev-1)%a * (1.0-pr1) + vert%levs(iLev)%a * pr1, &
                                  & vert%levs(iLev-1)%b * (1.0-pr1) + vert%levs(iLev)%b * pr1)
                return
              endif
              pr1 = pr2
            enddo  ! over levels of the vert
            !
            ! No hybrid levels found. However, this can be just deficiency of the vertical.
            ! If the level itself is not stupid-looking, let's accept it without warning
            !
            levelOut = fu_set_hybrid_level(int_missing, levTmp%a, 0.0)
            if(levIn%a < 10000. .and. levIn%a > 102000.)then
              call msg_warning('Failed to find proper hybrid level for input one, fake some')
              call msg('levelIn:')
              call report(levIn)
              call msg('vertical:')
              call report(vert,.true.)
              call msg('levelOut:')
              call report(levelOut)
            endif

          case default
            call report(vert)
            call report(levTmp)
            call set_error('Can not project above level to above vertical', &
                         & 'fu_level_to_vertical')
        end select
        
      case default
        call report(vert)
        call set_error('Type of the above vertical is not handled yet', &
                     & 'fu_level_to_vertical')
        return
    end select

  end function fu_level_to_vertical


  !****************************************************************

  subroutine change_vertical_type(vertInOut, newType)
    !
    ! The sub changes the type of the given vertical
    !
    implicit none

    ! Imported parameters
    type(silam_vertical), intent(inout) :: vertInOut
    integer, intent(in) :: newType

    ! Local variables
    integer :: iLev

    ! If the new type is the same as previous one, exit
    !
    if(vertInOut%vert_type == newType)return
    
    ! Change the vertical type and project level-by-level to the new system
    !
    vertInOut%vert_type = newType
    do iLev = 1, vertInOut%nLevs
      vertInOut%levs(iLev) = fu_level_to_vertical_crude(vertInOut%levs(iLev), vertInOut)
    end do

  end subroutine change_vertical_type

  subroutine remove_last_level(vert)
    !
    ! The sub just removes last level (layer from vertical)
    ! st reduce nLevs without reallocating anything
    ! No check for arrangeent
    implicit none

    ! Imported parameters
    type(silam_vertical), intent(inout) :: vert

    ! Local variables
    integer :: iLev

    ! Dirty hack: just reduce nLevs without reallocating anything
    !
    if (vert%nLevs > 1) then 
         vert%nLevs = vert%nLevs-1
    else
       call set_error("Can't cut upper level from vertical","remove_last_level")
       call report(vert, .true.)
    endif

  end subroutine remove_last_level
  
  !*****************************************************************
  
  logical function fu_if_cut_vertical_size_from_above(vert_to_cut, vert_to_cover, max_relief_height_m)
    !
    ! Cuts the levels from above and below of the vert_to_cut, still leaving the safety
    ! margin to cover the given vert_to_cover with some reserve for actual meteodata
    !
    implicit none
    
    ! Imported parameters
    type(silam_vertical), intent(inout) :: vert_to_cut
    type(silam_vertical), intent(in) :: vert_to_cover
    real, intent(in) :: max_relief_height_m
    
    ! Local variables
    integer :: i, nLevsNew
    real :: fTmp, fTmp2
    real, dimension(:), pointer :: arTmp
    !
    ! Finally, if the out_vertical to cover by meteodata is provided, we may have a chance to reduce
    ! the number levels. Just project them one-by-one and see which one is not needed
    !
    fu_if_cut_vertical_size_from_above = .false.
    if(.not. defined(vert_to_cover))return
    arTmp => fu_work_array()
    if(error)return
    call msg('Attempting to cut the vertical size. Current one is:')
    call report(vert_to_cut)
    call msg('Max relief height [m] and vertical to be covered are:', max_relief_height_m)
    call report(vert_to_cover)
    !
    ! First, do it without bringing relief height into the game
    !
    nLevsNew = vert_to_cut%nLevs
    call vert_interp_data_crude(vert_to_cut, vert_to_cover, arTmp, fTmp)
    do i = vert_to_cut%nLevs, 1, -1
      fTmp2 = fu_project_level(vert_to_cut%levs(i),vert_to_cover,arTmp,fTmp)
      if(fu_project_level(vert_to_cut%levs(i), vert_to_cover, arTmp, fTmp) < &
                                                                  & vert_to_cover%nLevs+0.499)then
        nLevsNew = i
        call msg('Cut vertical from above without relief. Suggested nLevels =', nLevsNew)
        fu_if_cut_vertical_size_from_above = .true.
        exit
      endif
    end do
    !
    ! Now, repeat the same with relief height included. The bigger of the two will go, plus reserve
    !
    call vert_interp_data_crude(vert_to_cut, vert_to_cover, arTmp, fTmp, max_relief_height_m + 1500. )!+1.5km
    do i = vert_to_cut%nLevs, 1, -1
      fTmp2 = fu_project_level(vert_to_cut%levs(i),vert_to_cover,arTmp,fTmp)
      if(fu_project_level(vert_to_cut%levs(i), vert_to_cover, arTmp, fTmp) < &
                                                                  & vert_to_cover%nLevs+0.499)then
call msg('Suggestion with relief:',i)
        nLevsNew = max(nLevsNew, i)
        call msg('Cut vertical from above with relief. New nLevels =', nLevsNew)
        fu_if_cut_vertical_size_from_above = .true.
        exit
      endif
    end do
    !Dirty hack 
    vert_to_cut%nLevs = min(vert_to_cut%nLevs, nLevsNew + 3)
    call free_work_array(arTmp)

  end function fu_if_cut_vertical_size_from_above

  
  !*****************************************************************

  function fu_lower_boundary_of_layer(layer) result(level)
    !
    ! Returns the lower boundary of the layer between two levels (bottom)
    !
    implicit none

    ! Return value of the function
    type(silja_level) :: level

    ! Imported parameter
    type(silja_level), intent(in) :: layer

    select case(layer%leveltype)
      case (layer_btw_2_pressure)
        level = fu_set_pressure_level(layer%b)

      case (layer_btw_2_altitude)
        level = fu_set_constant_altitude_level(layer%b)

      case (layer_btw_2_height)
        level = fu_set_constant_height_level(layer%b)

      case (layer_btw_2_sigma)
        level = fu_set_sigma_level(layer%b)

      case (layer_btw_2_hybrid)
        level = fu_set_hybrid_level(layer%number2, layer%a2, layer%b2)

      case default
        call set_error('Inconsistent level type of layer','fu_lower_boundary_of_layer')
    end select

  end function fu_lower_boundary_of_layer


  !*****************************************************************

  function fu_upper_boundary_of_layer(layer) result(level)
    !
    ! Returns the lower boundary of the layer between two levels (top)
    !
    implicit none

    ! Return value of the function
    type(silja_level) :: level

    ! Imported parameter
    type(silja_level), intent(in) :: layer

    select case(layer%leveltype)
      case (layer_btw_2_pressure)
        level = fu_set_pressure_level(layer%a)

      case (layer_btw_2_altitude)
        level = fu_set_constant_altitude_level(layer%a)

      case (layer_btw_2_height)
        level = fu_set_constant_height_level(layer%a)

      case (layer_btw_2_sigma)
        level = fu_set_sigma_level(layer%a)

      case (layer_btw_2_hybrid)
        level = fu_set_hybrid_level(layer%number, layer%a, layer%b)

      case default
        level=level_missing
        call set_error('Inconsistent level type of layer','fu_lower_boundary_of_layer')
    end select

  end function fu_upper_boundary_of_layer


  !***************************************************************************************

  real function fu_top_of_layer_value_lyr(layer) result(top_value)
    !
    ! A quick function that delivers the VALUE of the layer top, local unit
    ! No checking, no conversion, nothing. Just selection of the layer type and returning
    ! its top.
    !
    implicit none

    ! Imported parameter
    type(silja_level), intent(in) :: layer

    select case(layer%leveltype)
      case (layer_btw_2_pressure, layer_btw_2_sigma)
        top_value = layer%b

      case (layer_btw_2_altitude, layer_btw_2_height)
        top_value = layer%a

      case (layer_btw_2_hybrid)
        top_value = real(layer%number)

      case default
        top_value = real_missing
        call set_error('Inconsistent level type of layer','fu_top_of_layer_value_lyr')
    end select

  end function fu_top_of_layer_value_lyr


  !***************************************************************************************

  real function fu_top_of_layer_value_vert(vert, ind) result(top_value)
    !
    ! A quick function that delivers the VALUE of the layer top, local unit
    ! No checking, no conversion, nothing. Just selection of the layer type and returning
    ! its top. Layer is taken from the vertical using the ind as index
    !
    implicit none

    ! Imported parameter
    type(silam_vertical), intent(in) :: vert
    integer, intent(in) :: ind

    select case(vert%vert_type)
      case (layer_btw_2_pressure, layer_btw_2_sigma)
        top_value = vert%levs(ind)%b

      case (layer_btw_2_altitude, layer_btw_2_height)
        top_value = vert%levs(ind)%a

      case (layer_btw_2_hybrid)
        top_value = real(vert%levs(ind)%number2)

      case default
        top_value = real_missing
        call set_error('Inconsistent level type of layer','fu_top_of_layer_value_vert')
    end select

  end function fu_top_of_layer_value_vert


  !***************************************************************************************

  real function fu_bottom_of_layer_value_lyr(layer) result(bottom_value)
    !
    ! A quick function that delivers the VALUE of the layer top, local unit
    ! No checking, no conversion, nothing. Just selection of the layer type and returning
    ! its top.
    !
    implicit none

    ! Imported parameter
    type(silja_level), intent(in) :: layer

    select case(layer%leveltype)
      case (layer_btw_2_pressure, layer_btw_2_sigma)
        bottom_value = layer%a

      case (layer_btw_2_altitude, layer_btw_2_height)
        bottom_value = layer%b

      case (layer_btw_2_hybrid)
        bottom_value = real(layer%number)

      case default
        bottom_value = real_missing
        call set_error('Inconsistent level type of layer','fu_bottom_of_layer_value_lyr')
    end select

  end function fu_bottom_of_layer_value_lyr


  !***************************************************************************************

  real function fu_bottom_of_layer_value_vert(vert, ind) result(bottom_value)
    !
    ! A quick function that delivers the VALUE of the layer top, local unit
    ! No checking, no conversion, nothing. Just selection of the layer type and returning
    ! its top. Layer is taken from the vertical using the ind as index
    !
    implicit none

    ! Imported parameter
    type(silam_vertical), intent(in) :: vert
    integer, intent(in) :: ind

    select case(vert%vert_type)
      case (layer_btw_2_pressure, layer_btw_2_sigma)
        bottom_value = vert%levs(ind)%a

      case (layer_btw_2_altitude, layer_btw_2_height)
        bottom_value = vert%levs(ind)%b

      case (layer_btw_2_hybrid)
        bottom_value = real(vert%levs(ind)%number)

      case default
        bottom_value = real_missing
        call set_error('Inconsistent level type of layer','fu_bottom_of_layer_value_vert')
    end select

  end function fu_bottom_of_layer_value_vert


  !***************************************************************************************

  real function fu_layer_centre_value_lyr(layer) result(centre_value)
    !
    ! A quick function that delivers the VALUE of the layer top, local unit
    ! No checking, no conversion, nothing. Just selection of the layer type and returning
    ! its top.
    !
    implicit none

    ! Imported parameter
    type(silja_level), intent(in) :: layer

    select case(layer%leveltype)
      case (layer_btw_2_pressure, layer_btw_2_altitude, &
          & layer_btw_2_height, layer_btw_2_sigma)
        centre_value = (layer%a + layer%b) * 0.5

      case default
        centre_value = real_missing
        call set_error('Unsupported level type of layer','fu_layer_centre_value_lyr')
    end select

  end function fu_layer_centre_value_lyr


  !***************************************************************************************

  real function fu_layer_centre_value_vert(vert, ind) result(centre_value)
    !
    ! A quick function that delivers the VALUE of the layer top, local unit
    ! No checking, no conversion, nothing. Just selection of the layer type and returning
    ! its top. Layer is taken from the vertical using the ind as index
    !
    implicit none

    ! Imported parameter
    type(silam_vertical), pointer :: vert
    integer, intent(in) :: ind

    select case(vert%vert_type)
      case (layer_btw_2_pressure, layer_btw_2_altitude, &
          & layer_btw_2_height, layer_btw_2_sigma)
        centre_value = (vert%levs(ind)%a + vert%levs(ind)%b) * 0.5

      case (layer_btw_2_hybrid)
        centre_value = (real(vert%levs(ind)%number) + real(vert%levs(ind)%number))*0.5

      case default
        centre_value = real_missing
        call set_error('Inconsistent level type of layer','fu_layer_centre_value_vert')
    end select

  end function fu_layer_centre_value_vert


  !**************************************************************************************


  ! ***************************************************************

  ! ***************************************************************
  !
  !
  !      Private functions and subroutines.
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  recursive LOGICAL FUNCTION fu_cmp_levs_eq(lev1, lev2) result(eq)
!  LOGICAL FUNCTION fu_cmp_levs_eq(lev1, lev2) 
    !
    ! Compares two levels to be the same, or represent the same
    ! elevation. For example, surface level and constant height, zero
    ! elevation are to be the same.
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    !
    ! Imported parameters with intent(in):
    TYPE(silja_level), INTENT(in) :: lev1, lev2

    
    IF (defined(lev1)) THEN
      IF (.NOT.defined(lev2)) THEN
        eq = .false. ! 1 defined, 2 undefined
        RETURN
      END IF
    ELSE
      IF (defined(lev2)) THEN
        eq = .false. ! 1 undefined, 2 defined
      ELSE
        ! Both undefined, then they are also equal:
        eq = .true.
        RETURN
      END IF
    END IF

!    ! Both defined, check level types:
!    IF  (.NOT.fu_leveltypes_eq(lev1, lev2)) THEN
!      eq = .false.
!      RETURN
!    END IF

    SELECT CASE (lev1%leveltype)

      CASE (hybrid)

        if(lev2%leveltype /= hybrid)then
          eq=.false.
          return
        endif

        IF (lev1%number == lev2%number) THEN

          ! If hybrid coefficients known, then the also must match.
          ! Otherwise levels are considered equal, if level number is the same.
          IF (lev1%hybrid_coeff_known .and.&  ! coeff_known
                         & lev2%hybrid_coeff_known) THEN
            IF ((lev1%a .eps. lev2%a).and.(lev1%b .eps. lev2%b)) THEN
              eq = .true.
            ELSE
              eq = .false.
            END IF
          ELSE
            ! Same level numbers and coefficients not known:
            eq = .true.
          END IF ! coeff_known

        ELSE
          ! Different level numbers:
          eq = .false.
        END IF

      CASE (surface) != 001
        eq = lev2%leveltype == surface .or. &
           & (lev2%leveltype == constant_height .and. (lev2%a .eps. 0.)).or. &
           & (lev2%leveltype == depth_level .and. (lev2%a .eps. 0.)).or. &
           & (lev2%leveltype == sigma_level .and. (lev2%a .eps. 1.))

      CASE (top_of_the_atmosphere) ! = 008
        eq = lev2%leveltype == top_of_the_atmosphere

      CASE (layer_btw_2_pressure, &
          & layer_btw_2_hybrid, &
          & layer_btw_2_altitude, &
          & layer_btw_2_depth, &
          & layer_btw_2_height, &
          & layer_btw_2_sigma) ! = 101
        if (.not. fu_if_layer(lev2)) then
          eq = .false.
        else
          eq = fu_cmp_levs_eq(fu_upper_boundary_of_layer(lev1),fu_upper_boundary_of_layer(lev2)) .and. &
               & fu_cmp_levs_eq(fu_lower_boundary_of_layer(lev1),fu_lower_boundary_of_layer(lev2)) 
        end if


      CASE (mean_sea) ! = 102
        eq = lev2%leveltype == mean_sea .or. &
           & (lev2%leveltype == constant_altitude .and. (lev2%a .eps. 0.))

      CASE (constant_pressure)
        eq = (lev2%leveltype == constant_pressure .and. (lev2%a .eps. lev1%a))

      case (constant_altitude)
        if(lev2%leveltype == constant_altitude)then
          eq = (lev2%a .eps. lev1%a)
        else
          eq = lev2%leveltype == mean_sea .and. (lev1%a .eps. 0.)
        endif

      CASE (constant_height) ! = 105
        if(lev2%leveltype == constant_height)then
          eq = (lev2%a .eps. lev1%a)
        else
          eq = lev2%leveltype == surface .and. (lev1%a .eps. 0.)
        endif

      CASE (sigma_level) ! = 107
        if(lev2%leveltype ==  sigma_level)then
          eq = (lev2%a .eps. lev1%a)
        else
          eq = lev2%leveltype == surface .and. (lev1%a .eps. 0.)
        endif

      CASE (depth_level) ! = 111
        if(lev2%leveltype == depth_level)then
          eq = (lev2%a .eps. lev1%a)
        else
          eq = lev2%leveltype == surface .and. (lev1%a .eps. 0.)
        endif

      CASE (entire_atmosphere_single_layer) ! = 200
        eq = lev2%leveltype == entire_atmosphere_single_layer .and. &
           & lev2%number == lev1%number

      CASE (no_level) ! = 0
        eq = lev2%leveltype == no_level

      CASE (any_level) ! = 333
        eq = lev2%leveltype == any_level

      CASE default
        call set_error('Unknown level type', 'fu_cmp_lev_eq')
        eq=.false.
	
    END SELECT

  END FUNCTION fu_cmp_levs_eq


  !****************************************************************

  LOGICAL FUNCTION fu_cmp_verts_eq(vert1, vert2) result(eq)
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silam_vertical), INTENT(in) :: vert1, vert2

    ! Local parameters
    integer :: i

    eq = .false.
    if((vert1%vert_type /= vert2%vert_type) .or. (vert1%NLevs /= vert2%NLevs))return
    
    do i=1, vert1%NLevs
      if(.not. fu_cmp_levs_eq(vert1%levs(i),vert2%levs(i)))return
    end do
    eq = .true.

  end function fu_cmp_verts_eq


  ! ***************************************************************


  recursive LOGICAL FUNCTION fu_first_lower_sametype(lev1, lev2)&
      & result(first_lower)

    ! Description: This function gives true value if the first of
    ! given levels is lower in height (closer to ground). Works only
    ! for levels of same type. If the levels are layers, the
    ! comparison is done with the midpoints.
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent(in):
    TYPE(silja_level), INTENT(in) :: lev1, lev2

    IF (.NOT.(fu_leveltypes_eq(lev1, lev2))) THEN
      CALL set_error('not same type levels',&
	  & 'fu_first_lower_sametype')
      RETURN
    END IF

    if (fu_if_layer(lev1)) then
      first_lower = fu_first_lower_sametype(fu_central_level_of_layer(lev1), &
                                          & fu_central_level_of_layer(lev2))
      return
    end if

    leveltype: SELECT CASE (lev1%leveltype)

      CASE(constant_pressure)
      IF (lev1%a > lev2%a) THEN
        first_lower = .true.
      ELSE
        first_lower = .false.
      END IF


      CASE(constant_altitude, constant_height)
      IF (lev1%a < lev2%a) THEN
        first_lower = .true.
      ELSE
        first_lower = .false.
      END IF


      CASE(sigma_level)
      IF (lev1%a > lev2%a) THEN
        first_lower = .true.
      ELSE
        first_lower = .false.
      END IF

    CASE(hybrid) ! for hybrid levels no. 1 is highest and last
      if (lev1%hybrid_coeff_known .and. lev2%hybrid_coeff_known) then
        first_lower = lev1%a + lev1%b*std_pressure_sl > lev2%a + lev2%b*std_pressure_sl
      else
        first_lower = lev1%number > lev2%number
      end if

    CASE default 

      CALL set_error('level1 not defined', &
	  & 'fu_first_lower_sametype')
      RETURN

    END SELECT leveltype


  END FUNCTION fu_first_lower_sametype



  ! ***************************************************************


  LOGICAL FUNCTION fu_first_lower_eq_sametype(lev1, lev2)

    ! Description:
    ! This function gives true value if the first of given levels is 
    ! upper in height or the levels are equal. Works only for levels
    ! of same type.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent(in):
    TYPE(silja_level), INTENT(in) :: lev1, lev2

!    IF (fu_first_lower_sametype(lev1, lev2).or.(lev1 == lev2)) THEN
    IF (fu_first_lower_sametype(lev1, lev2).or.fu_cmp_levs_eq(lev1,lev2)) THEN
      fu_first_lower_eq_sametype = .true.
    ELSE
      fu_first_lower_eq_sametype = .false.
    END IF

  END FUNCTION fu_first_lower_eq_sametype



  ! ***************************************************************


  LOGICAL FUNCTION fu_first_upper_sametype(lev1, lev2)

    ! Description:
    ! This function gives true value if the first of given levels is 
    ! upper in height. Works only for levels of same type.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent(in):
    TYPE(silja_level), INTENT(in) :: lev1, lev2

    fu_first_upper_sametype = .NOT.&
	& (fu_first_lower_eq_sametype(lev1, lev2))

  END FUNCTION fu_first_upper_sametype



  ! ***************************************************************


  LOGICAL FUNCTION fu_first_upper_eq_sametype(lev1, lev2)

    ! Description:
    ! This function gives true value if the first of given levels is 
    ! upper in height or the levels are equal. Works only for levels
    ! of same type.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent(in):
    TYPE(silja_level), INTENT(in) :: lev1, lev2

    fu_first_upper_eq_sametype = .NOT.&
	& (fu_first_lower_sametype(lev1, lev2))

  END FUNCTION fu_first_upper_eq_sametype


  !****************************************************************

  function fu_level_from_vertical(vert, levIndex, ifSingleLevelNeeded) result(level)
    !
    ! Picks the requested level from the vertical structure and
    ! computes the central-level for the thick layer, if needed
    !
    implicit none

    ! Returns level
    type(silja_level) :: level

    ! Imported parameters with intent IN
    type(silam_vertical), intent(in) :: vert
    integer, intent(in) :: levIndex
    logical, intent(in), optional :: ifSingleLevelNeeded

    if(levIndex > vert%NLevs .or. levIndex < 1) then
      call set_error('Incorrect index of the level','fu_level_from_vertical')
      call msg('Index', levIndex)
      call report(vert, .true.)
      level = level_missing
      return
    end if

    level = vert%levs(levIndex)

    if(present(ifSingleLevelNeeded))then
      if(ifSingleLevelNeeded)then
        level = fu_central_level_of_layer(level)
      endif
    endif

  end function fu_level_from_vertical


  !****************************************************************

  integer function fu_level_nbr_from_vertical(vert, level) 
    !
    ! Finds the requested level from the vertical structure
    !
    implicit none

    ! Imported parameters with intent IN
    type(silam_vertical), intent(in) :: vert
    type(silja_level), intent(in) :: level

    ! Local variables
    integer :: i

    fu_level_nbr_from_vertical = int_missing

    if(fu_level_belongs_to_vertical(level, vert))then
      do i = 1, vert%nLevs
 !      if(vert%levs(i) == level)then
	if(fu_cmp_levs_eq(vert%levs(i),level))then
          fu_level_nbr_from_vertical = i
          return
        endif
      end do
    endif

  end function fu_level_nbr_from_vertical


  !****************************************************************

  real function fu_range_of_vertical(vert)
    !
    ! Returns roughly estimated range covered by the vertical
    ! structure. IMPORTANT is that the unit of this range depends
    ! on the type of the layers.
    !
    implicit none

    ! Imported parameter with intent IN
    type(silam_vertical), intent(in) :: vert

    select case (vert%vert_type)
      case(surface, mean_sea)
        fu_range_of_vertical =0. ! No range - single-height levels

      case(constant_pressure, layer_btw_2_pressure)
        fu_range_of_vertical = fu_pr_level_pressure(fu_bottom_level_vert(vert)) - &
                              & fu_pr_level_pressure(fu_top_level_vert(vert))

      case(hybrid, layer_btw_2_hybrid)
        fu_range_of_vertical = fu_hybrid_level_pressure(fu_bottom_level_vert(vert),101000.) - &
                              & fu_hybrid_level_pressure(fu_top_level_vert(vert),101000.)

      case(constant_altitude, layer_btw_2_altitude)
        fu_range_of_vertical = fu_level_altitude(fu_top_level_vert(vert)) - &
                             & fu_level_altitude(fu_bottom_level_vert(vert))

      case(constant_height, layer_btw_2_height)
        fu_range_of_vertical = fu_level_height(fu_top_level_vert(vert)) - &
                             & fu_level_height(fu_bottom_level_vert(vert))

      case(sigma_level, layer_btw_2_sigma)
        fu_range_of_vertical = fu_sigma_level_sigma(fu_bottom_level_vert(vert)) - &
                             & fu_sigma_level_sigma(fu_top_level_vert(vert))

      case(depth_level, layer_btw_2_depth)
        fu_range_of_vertical = fu_depth_level_depth(fu_top_level_vert(vert)) - &
                             & fu_depth_level_depth(fu_bottom_level_vert(vert))

      case default
        call msg_warning('Unknown vertical type','fu_range_of_vertical')
        fu_range_of_vertical = real_missing
    end select

  end function fu_range_of_vertical


  !***********************************************************************

  function fu_top_level_vert(vert) result(level)
    !
    ! Returns the highest level in the vertical structure
    !
    implicit none

    ! Return value of the function
    type(silja_level) :: level

    ! Imported parameters with intent IN
    type(silam_vertical), intent(in) :: vert

    ! Local variables
    integer :: i

    if(.not.defined(vert) .or. vert%NLevs < 1)then
      level = level_missing
      call set_error('No levels or undefined vertical structure','fu_top_level_of_vert')
      return
    end if
    level = vert%levs(1)
    do i = 2, vert%NLevs
 !      if(.not.defined(vert%levs(i)) .or. vert%levs(i) == level_missing)return
       if(.not.defined(vert%levs(i)) .or. fu_cmp_levs_eq(vert%levs(i),level_missing))return
       if(vert%levs(i) > level)level = vert%levs(i)
    end do

  end function fu_top_level_vert


  !***********************************************************************

  function fu_bottom_level_vert(vert) result(level)
    !
    ! Returns the lowest level in the vertical structure
    !
    implicit none

    ! Return value of the function
    type(silja_level) :: level

    ! Imported parameters with intent IN
    type(silam_vertical), intent(in) :: vert

    ! Local variables
    integer :: i

    if(.not.defined(vert) .or. vert%NLevs < 1)then
      level = level_missing
      call set_error('No levels or undefined vertical structure','fu_bottom_level_of_vert')
      return
    end if
    level = vert%levs(1)
    do i = 2, vert%NLevs
       if(.not.defined(vert%levs(i)) .or. fu_cmp_levs_eq(vert%levs(i),level_missing))return
       if(vert%levs(i) < level)level = vert%levs(i)
    end do

  end function fu_bottom_level_vert


  !************************************************************************************

  subroutine reproject_verticals(vertFrom, fractFrom, vertTo, FractTo, fMassCentreTo, nLevsToActive, &
                               & ifMassCentreInRelUnit)

    !
    ! Scans the vertical TO layer by layer projecting each of them to the vertical FROM.
    ! Reason: each layer has to be distributed among all layers in the verticalFrom. 
    ! For that they must be projected onto it, not the other way round.
    !
    implicit none
    
    ! Imported parameters
    type(silam_vertical), intent(in) :: vertFrom, vertTo
    real, dimension(:), intent(in) :: fractFrom
    real, dimension(:), intent(out) :: fractTo, fMassCentreTo
    integer, intent(out) :: nLevsToActive
    logical, intent(in) :: IfMassCentreInRelUnit

    ! Local variables
    integer :: i, j
    real :: fTopInd, fBottomInd, fTopIndMom, fBottomIndMom, fOverlap, fCentreMom
    type(silja_level) :: lev_top, lev_bot
    real, dimension(max_levels) :: met_data_col, met_data_col_inv
    real :: met_data_srf, met_data_srf_inv

    nLevsToActive = -1
    call vert_interp_data_crude(vertTo, vertFrom, met_data_col, met_data_srf)
    call vert_interp_data_crude(vertFrom, vertTo, met_data_col_inv, met_data_srf_inv)
    if (error) return

    do i = 1,fu_NbrOfLevels(vertTo)
      !
      ! Project the layer of the given vertical to the source vertical 
      !
      lev_top = fu_upper_boundary_of_layer(fu_level(vertTo,i))
      lev_bot = fu_lower_boundary_of_layer(fu_level(vertTo,i))
      if (error) return
      fTopInd = fu_project_level(lev_top, vertFrom, met_data_col, met_data_srf)
      fBottomInd = fu_project_level(lev_bot, vertFrom, met_data_col, met_data_srf)
      if(error)return
      !
      ! Having upper and lower indices of the given vertical layer in the original 
      ! source vertical, sum-up the fractions that are within this layer
      ! To do that, scan one-by-one the source layers looking for the overlap
      !
      do j= max(1, int(fBottomInd)), min(fu_NbrOfLevels(vertFrom),int(fTopInd) + 1)
        ! Find out the overlap between the levels
        !
        fOverlap = max(0.,(min(fTopInd,real(j)+0.5) - max(fBottomInd,real(j)-0.5)))
        
        if(fOverlap > 0)then
          !
          ! Mass fraction is the overlap multiplied with mass fraction into the original layer
          !
          fractTo(i) = fractTo(i) + fractFrom(j) * fOverlap
          !
          ! Centre of mass in verticalTO requires projection of layer FROM to vertical TO:
          !
          lev_top = fu_upper_boundary_of_layer(fu_level(vertFrom,j))
          fTopIndMom = fu_project_level(lev_top, vertTo, met_data_col_inv, met_data_srf_inv) - i
          lev_bot = fu_lower_boundary_of_layer(fu_level(vertFrom,j))
          fBottomIndMom = fu_project_level(lev_bot, vertTo, met_data_col_inv, met_data_srf_inv) - i
          fCentreMom = 0.5*(max(fBottomIndMom,-0.5) + min(fTopIndMom,0.5))
          !
          ! The first moment of the emitted mass fraction is the centre of mass and the fraction itself
          ! Note: moment is RELATIVE with regard to vertTo layers. I.e., its centre of mass
          ! varies from -0.5 to 0.5.
          !
          fMassCentreTo(i) = fMassCentreTo(i) + fractFrom(j) * fOverlap * fCentreMom
        endif
      enddo
    end do
    !
    ! Having the vertical fractions filled-in, let's reduce, if possible, the number of layers
    ! to be considered: those with zero fractions are not interesting. 
    !
    do i = vertTo%nLevs, 1, -1
      if(fractTo(i) > 0.0001)then
        nLevsToActive = i
        exit
      endif
    end do
    if (nLevsToActive < 1) then
      call msg_warning("attempt to project beyond the vertical")
      do i = 1, vertTo%nLevs
        call msg("level fraction", fractTo(i), fMassCentreTo(i))
      end do
      nLevsToActive = 0
      call msg("levtop,fbot", fTopInd, fBottomInd)
      call msg("vertTo%nLevs",vertTo%nLevs)
      call msg("nLevsToActive",nLevsToActive)
      call msg("vertTo")
      do i = 1, vertTo%nLevs
        call report(fu_level(vertTo,i))
      end do
      call msg("vertFrom")
      do i = 1, vertFrom%nLevs
        call report(fu_level(vertFrom,i))
      end do
      call  set_error('Mass centre is outside the above layer 1234','reproject_verticals')
      return
    endif
!
! We got the vertical fractions and centres of masses in relative coordinates. This is enough
! in the current setup: the storage goes in relative indices [-0.5, 0.5] while calculations are
! in absolute units.
!
    !!
    !! Now we have vertical mass centre in relaive co-ordinates while the actual calculations 
    !! go on in absolute (but still related to the centre of the level). 
    !! Conversion is then straightforward:
    !!
    !do i = 1, nLevsToActive
    !  fMassCentreTo(i) = (fMassCentreTo(i) / (fractTo(i) + 1.e-20))
    !  if(abs(fMassCentreTo(i)) > 0.5)then
    !    call report(fu_level(vertTo,i))
    !    call msg('Mass centre is outside the above layer:',fMassCentreTo(i))
    !    call set_error('Mass centre is outside the above layer','reproject_verticals')
    !    return
    !  else if (.not. ifMassCentreInRelUnit) then
    !    fMassCentreTo(i) = fMassCentreTo(i) * fu_layer_thickness_local_unit(fu_level(vertTo,i))
    !  endif
    !enddo


  end subroutine reproject_verticals


  !********************************************************************************

  subroutine levels_to_layers(level_vertical, layer_vertical, first_level)
    !
    ! Attempt to create a vertical of thick layers based on the
    ! vertical defined by the midpoints and the first_level. The first
    ! level is the top of the highest layer for pressure vertical, and
    ! the bottom of the first layer for all others. 
    !
    ! By default, the first level is zero pressure for pressure levels,
    ! the ground (a=0, b=1 or height=0) for height or hybrid layers,
    ! and sea level for altitude levels.
    !
    ! If the levels cannot be transformed into layers, an error is set.
    !
    implicit none
    type(silam_vertical), intent(in) :: level_vertical
    type(silam_vertical), intent(out) :: layer_vertical
    type(silja_level), intent(in), optional :: first_level

    integer :: ilev, itmp, n_layers, stat
    type(silja_level) :: top, bottom, level
    real :: a, b, pres, height
    real, parameter :: ps_test = 101300
    type(silja_level), dimension(:), allocatable :: layers

    n_layers = fu_nbrOfLevels(level_vertical)

    allocate(layers(n_layers), stat=stat)
    if (stat /= 0) then
      call set_error('Allocate failed', 'levels_to_layers')
      return
    end if

    select case (fu_leveltype(level_vertical))
    case (hybrid)
      ! 
    
      if (present(first_level)) then
        bottom = first_level
      else
        bottom = fu_set_level(hybrid, fval1=0.0, fval2=1.0, ival1=n_layers)
      end if
      ! We assume that the levels are defined from the original layers by
      ! taking arithmetic averages of a and b.
      do ilev = 1, fu_nbrOfLevels(level_vertical)
        level = fu_level(level_vertical, ilev)
        a = level%a
        b = level%b
        if (a + b*ps_test > bottom%a + bottom%b*ps_test) then
          call msg("Problem:")
          call report(level)
          call msg("Bottom:")
          call report(bottom)
          call msg("Vertical:")
          do itmp = 1, fu_nbrOfLevels(level_vertical)
           call report(fu_level(level_vertical, itmp))
          enddo
          call set_error('Inconsistent hybrid levels', 'levels_to_layers')
          return
        end if
        ! The new b
        b = 2*b - bottom%b
        if (b < 0.0) b = 0.0 ! may happen due to lack of precision
        top = fu_set_level(hybrid, fval1=2*a - bottom%a, fval2=b, ival1=n_layers-ilev)
        layers(ilev) = fu_set_layer_between_two(top, bottom)
        bottom = top
      end do
            
    case (constant_pressure)
      if (present(first_level)) then
        top = first_level 
      else
        call msg_warning('Constructing a layer vertical assuming top pressure==0.0 Pa', 'levels_to_layers')
        top = fu_set_pressure_level(0.0)
      end if

      do ilev = fu_nbrOfLevels(level_vertical), 1, -1
        level = fu_level(level_vertical, ilev)
        pres = fu_pr_level_pressure(level)
        if (pres < fu_pr_level_pressure(top)) then
          call set_error('Inconsistent pressure levels', 'levels_to_layers')
          return
        end if
        bottom = fu_set_pressure_level(2*pres - fu_pr_level_pressure(top))
        layers(ilev) = fu_set_layer_between_two(top, bottom)
        top = bottom
      end do

    case (constant_height)
      if (present(first_level)) then
        bottom = first_level
      else
        bottom = fu_set_constant_height_level(0.0)
      end if

      do ilev = 1, fu_nbrOfLevels(level_vertical)
        level = fu_level(level_vertical, ilev)
        height = fu_level_height(level)
        if (height < fu_level_height(bottom)) then
          call set_error('Inconsistent constant height levels', 'levels_to_layers')
          return
        end if
        top = fu_set_constant_height_level(2*height - fu_level_height(bottom))
        layers(ilev) = fu_set_layer_between_two(top, bottom)
        bottom = top
      end do

    case (constant_altitude)
      if (present(first_level)) then
        bottom = first_level
      else
        bottom = fu_set_constant_altitude_level(0.0)
      end if

      do ilev = 1, fu_nbrOfLevels(level_vertical)
        level = fu_level(level_vertical, ilev)
        height = fu_level_altitude(level)
        if (height < fu_level_altitude(bottom)) then
          call set_error('Inconsistent constant altitude levels', 'levels_to_layers')
          return
        end if
        top = fu_set_constant_altitude_level(2*height - fu_level_altitude(bottom))
        layers(ilev) = fu_set_layer_between_two(top, bottom)
        bottom = top
      end do
      
    case default
      call set_error('Unsupported level type', 'levels_to_layers')
      return

    end select
    
    call set_vertical(layers, layer_vertical)
    deallocate(layers)

  end subroutine levels_to_layers

  subroutine vert_to_metric(vert_in, vert_out)
    ! Create a height vertical, which (in standard atmosphere) matches the given hybrid or
    ! pressure vertical.
    implicit none
    type(silam_vertical), intent(in) :: vert_in
    type(silam_vertical), intent(out) :: vert_out
    
    logical :: if_layer
    integer :: nlevs, allocstat, ilev
    type(silja_level), dimension(:), allocatable :: new_levels
    type(silja_level) :: top, bottom, top_m, bottom_m, new_level

    if (fu_fails(defined(vert_in), 'vert_in not defined', 'vert_to_metric')) return
    nlevs = fu_NbrOfLevels(vert_in)
    if (fu_fails(nlevs > 0, 'vert_in is empty', 'vert_to_metric')) return

    if (fu_leveltype(vert_in) == constant_altitude .or. fu_leveltype(vert_in) == constant_height) then
      vert_out = vert_in
      return
    end if

    allocate(new_levels(nlevs), stat=allocstat)
    if (fu_fails(allocstat == 0, 'Allocate failed', 'vert_to_metric')) return

    if_layer = fu_if_layer(fu_level(vert_in, 1))
    do ilev = 1, nlevs
      if (if_layer) then
        top = fu_upper_boundary_of_layer(fu_level(vert_in, ilev))
        bottom = fu_lower_boundary_of_layer(fu_level(vert_in, ilev))
        top_m = fu_set_level(constant_height, fu_std_height(top))
        bottom_m = fu_set_level(constant_height, fu_std_height(bottom))
        new_level = fu_set_layer_between_two(bottom_m, top_m)
      else
        new_level = fu_set_level(constant_height, fu_std_height(fu_level(vert_in, ilev)))
      end if
      new_levels(ilev) = new_level
    end do
    
    call set_vertical(new_levels, vert_out)
    call arrange_levels_in_vertical(vert_out)
    deallocate(new_levels)

  end subroutine vert_to_metric
  
  subroutine vert_to_pressure(vert_in, vert_out)
    ! Create a pressure vertical, which (in standard atmosphere) matches the given height
    ! altitude, or hybrid vertical.
    implicit none
    type(silam_vertical), intent(in) :: vert_in
    type(silam_vertical), intent(out) :: vert_out

    logical :: if_layer
    integer :: nlevs, allocstat, ilev
    type(silja_level), dimension(:), allocatable :: new_levels
    type(silja_level) :: top, bottom, top_Pa, bottom_Pa, new_level

    if (fu_fails(defined(vert_in), 'vert_in not defined', 'vert_to_metric')) return
    nlevs = fu_NbrOfLevels(vert_in)
    if (fu_fails(nlevs > 0, 'vert_in is empty', 'vert_to_metric')) return
    
    if (fu_leveltype(vert_in) == constant_pressure) then
      vert_out = vert_in
      return
    end if

    allocate(new_levels(nlevs), stat=allocstat)
    if (fu_fails(allocstat == 0, 'Allocate failed', 'vert_to_metric')) return

    if_layer = fu_if_layer(fu_level(vert_in, 1))

    do ilev = 1, nlevs
      if (if_layer) then
        top = fu_upper_boundary_of_layer(fu_level(vert_in, ilev))
        bottom = fu_lower_boundary_of_layer(fu_level(vert_in, ilev))
        top_Pa = fu_set_level(constant_pressure, fu_std_pressure(top))
        bottom_Pa = fu_set_level(constant_pressure, fu_std_pressure(bottom))
        new_level = fu_set_layer_between_two(bottom_Pa, top_Pa)
      else
        new_level = fu_set_level(constant_pressure, fu_std_pressure(fu_level(vert_in, ilev)))
      end if
      if (error) return
      new_levels(ilev) = new_level
    end do
    
    call set_vertical(new_levels, vert_out)
    call arrange_levels_in_vertical(vert_out)
    deallocate(new_levels)

  end subroutine vert_to_pressure
    
  real function fu_std_pressure(level) result(press)
    implicit none
    type(silja_level), intent(in) :: level
    
    real :: press_ratio, temperature, density
    
    select case (level%leveltype)
    case (constant_height, constant_altitude)
      call us_standard_atmosphere(fu_level_height(level), density, press_ratio, temperature)
      press = std_pressure_sl*press_ratio
    case (hybrid, sigma_level)
      press = fu_hybrid_level_pressure(level, std_pressure_sl)
    case default
      call set_error('Cannot handle this level', 'fu_std_pressure')
      return
    end select
  end function fu_std_pressure


  real function fu_std_height(level) result(height)
    implicit none
    type(silja_level), intent(in) :: level
    
    real :: press

    if (fu_fails(.not. fu_if_layer(level), 'fu_std_height works only for thin levels', 'fu_std_height')) return
    
    select case(fu_leveltype(level))
    case (surface, mean_sea)
      height = 0.0
    case (constant_altitude, constant_height)
      height = fu_level_height(level)
    case (constant_pressure)
      height = fu_height_for_press(fu_pr_level_pressure(level))
    case (hybrid, sigma_level)
      press = fu_hybrid_level_pressure(level, std_pressure_sl)
      height = fu_height_for_press(press)
    case default
      call report(level)
      call set_error('Cannot determine standard height for this type of level', 'fu_std_height')
    end select
  end function fu_std_height


  ! ***************************************************************

  ! ***************************************************************
  ! 
  !
  !                 TESTS
  !
  !
  ! ***************************************************************
  ! ***************************************************************

   function  fu_level_string(level) result (str)

    ! Description:
    ! Prints a report to screen describing  the level.
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
    TYPE(silja_level), INTENT(in) :: level
    character (len=fnlen) :: str

    integer, dimension(2) :: fUnits
    
    
    IF (.NOT. defined(level)) THEN
      WRITE(str,*)'Undefined level'
      RETURN
    END IF

    SELECT CASE (level%leveltype)

      case(surface)
        WRITE(str,*)'Surface level'

      case(mean_sea)
        WRITE(str,*)'Mean sea level'

      case(top_of_the_atmosphere)
        WRITE(str,*)'Top of the atmosphere'

      case(entire_atmosphere_single_layer)
        WRITE(str,*)'Entire atmosphere as a single layer'

      CASE (constant_pressure) 
        WRITE(str,*)' Constant pressure level ', NINT(level%a/100.),'HPa.'

      CASE (layer_btw_2_pressure) 
        WRITE(str,*)' Layer between two constant pressure levels', &
                & NINT(level%a/100.),'HPa.', NINT(level%b/100.),'HPa.'

      CASE (constant_altitude)
        WRITE(str,'(A,F10.2,A)')' Constant altitude level, height ',level%a,'m'

      CASE (layer_btw_2_altitude)
        WRITE(str,'(A,2(F10.2,A,1x))') &
                 & ' Layer between two constant altitude levels, height ',&
                 & level%a,'m', level%b,'m'

      CASE (constant_height)
        WRITE(str,'(A, F10.2, A)')' Constant height level, height ',level%a,'m'

      CASE (layer_btw_2_height)
        WRITE(str,'(A, 2(F10.2, A,1x))') &
              & ' Layer between two constant height levels, heights ', level%a,'m', level%b, 'm'

      CASE(hybrid)
        IF (level%hybrid_coeff_known) THEN
          WRITE(str,fmt = '(A, I4, F10.2, F12.7, F10.2)')&
                & ' Hybrid level number, coefficient a, b, stdpressure: ',&
                & level%number, level%a, level%b, level%a + level%b * std_pressure_sl 
        ELSE
          WRITE(str,fmt = '(A, I4, A)')&
                & ' Hybrid level number: ', level%number,' (a,b unknown)'
        END IF

      CASE(layer_btw_2_hybrid)
        IF (level%hybrid_coeff_known) THEN
          WRITE(str,fmt = '(A, 2(I4, F10.2, F12.7 F10.2))')&
                & ' Layer btw 2 hybrid: nbr, a, b, stdpress: ',&
                & level%number, level%a, level%b, level%a + level%b * std_pressure_sl, &
                & level%number2, level%a2, level%b2, level%a2 + level%b2 * std_pressure_sl
        ELSE
          WRITE(str,fmt = '(A, 2I4, A)')&
                & ' Layer between two hybrid levels: number: ', &
                & level%number,level%number2,' (a,b unknown)'
        END IF

      CASE (sigma_level)
        write(str, '(A, F12.7, F10.2)') 'Sigma level, coef, stdpress ', level%a, level%a * std_pressure_sl

      CASE (layer_btw_2_sigma)
        write(str,*)' Layer between two sigma levels:', level%a, level%b

      CASE (depth_level)
        write(str,*)' Depth level ', level%a

      CASE (layer_btw_2_depth)
        write(str,*)' Layer between two depth levels:', level%a, level%b

      CASE default
        write(str,*)'Unknown level. '  

    END SELECT

  END function  fu_level_string

  !***************************************************************** 
   
  SUBROUTINE print_level_report(level)

    ! Description:
    ! Prints a report to screen describing  the level.
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    TYPE(silja_level), INTENT(in) :: level

    integer, dimension(2) :: fUnits
    integer :: iUnit

    fUnits(1:2) = (/6, run_log_funit/)
   
    do iUnit = 1,2
       if (smpi_global_rank /= 0 .and. iUnit==1) cycle !Be quiet at stdut
       WRITE(funits(iUnit),fmt='(A)') trim(fu_level_string(level))
    enddo

  END SUBROUTINE print_level_report


  ! ***************************************************************

  SUBROUTINE print_vertical_report(vertical, ifFull)
    ! 
    ! Prints a report to screen describing  the level.
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    TYPE(silam_vertical), INTENT(in) :: vertical
    logical, intent(in), optional :: ifFull

    !Local variables
    integer :: iLev
    integer, dimension(2) :: fUnits
    integer :: iUnit

    IF (.NOT.defined(vertical)) THEN
      call msg('Undefined vertical')
      RETURN
    END IF



    fUnits(1:2) = (/6, run_log_funit/)
   
    do iUnit = 1,2
       if (smpi_global_rank /= 0 .and. iUnit==1) cycle !Be quiet at stdut
      
    SELECT CASE (vertical%vert_type)
      case(surface)
              WRITE(funits(iUnit),*)'Funny but surface vertical'

      case(mean_sea)
              WRITE(funits(iUnit),*)'Funny but mean sea vertical'

      case(top_of_the_atmosphere)
              WRITE(funits(iUnit),*)'Funny but top of the atmosphere vertical'

      CASE (constant_pressure) 
              WRITE(funits(iUnit),'(I4, A)') vertical%Nlevs, ' constant pressure levels'
      
      CASE (layer_btw_2_pressure) 
              WRITE(funits(iUnit),'(I4, A)') vertical%Nlevs, ' layers between constant pressure levels'

      CASE (constant_altitude)
              WRITE(funits(iUnit),'(I4, A)') vertical%Nlevs, ' constant altitude levels'

      CASE (layer_btw_2_altitude)
              WRITE(funits(iUnit),'(I4, A)') vertical%Nlevs, ' layers between constant altitude levels'

      CASE (constant_height)
              WRITE(funits(iUnit),'(I4, A)') vertical%Nlevs, ' constant height levels'

      CASE (layer_btw_2_height)
              WRITE(funits(iUnit),'(I4, A)') vertical%Nlevs, ' layers between constant height levels'

      CASE(hybrid)
              WRITE(funits(iUnit),'(I4, A)') vertical%Nlevs, ' hybrid levels'

      CASE(layer_btw_2_hybrid)
              WRITE(funits(iUnit),'(I4, A)') vertical%Nlevs, ' layers between hybrid levels'

      CASE (sigma_level)
              WRITE(funits(iUnit),'(I4, A)') vertical%Nlevs, ' sigma levels'

      CASE (layer_btw_2_sigma)
              WRITE(funits(iUnit),'(I4, A)') vertical%Nlevs, ' layers between sigma levels'

      CASE (depth_level)
              WRITE(funits(iUnit),'(I4, A)') vertical%Nlevs, ' depth levels'

      CASE (layer_btw_2_depth)
              WRITE(funits(iUnit),'(I4, A)') vertical%Nlevs, ' layers between depth levels'

      case(entire_atmosphere_single_layer)
        if(vertical%Nlevs /= 1)then
                WRITE(funits(iUnit),'(I4, A)') vertical%Nlevs, ' layers between depth levels'
        else
                WRITE(funits(iUnit),'(I4, A)') 'One layer of entire atmosphere'
        endif

      CASE default
              WRITE(funits(iUnit),'(I4, A)') vertical%Nlevs, ' levels of unknown type ',vertical%vert_type

    END SELECT
    enddo

    if(present(ifFull))then
      if(ifFull)then
        do iLev = 1, vertical%Nlevs
          call report(vertical%levs(iLev))
        enddo
      endif
    endif

  END SUBROUTINE print_vertical_report


  !**************************************************************************

  subroutine report_vertical_as_namelist(vert, iUnit)
    !
    ! Stores the level as the namelist to the given file, which must be open
    !
    implicit none

    TYPE(silam_vertical), intent(in) :: vert
    integer, intent(in) :: iUnit

    ! Local variables
    integer :: iLev
    !
    ! Stupidity check:
    !
    if(.not.defined(vert))then
      call set_error('Undefined vertical','report_vertical_as_namelist')
      write(iUnit,*)'vertical_method = UNDEFINED_VERTICAL'
      return
    endif

    !
    ! Depending on the vertical type, we print nesessary vertical_method
    !
    write(iUnit,*)
    select case (vert%vert_type)
      !
      ! Simple single-level verticals
      !
      case(surface)
        write(iUnit,*)'vertical_method = SURFACE_LEVEL'
      case(top_of_the_atmosphere)
        write(iUnit,*)'vertical_method = TOP_ATMOSPHERE_LEVEL'
      case(mean_sea)
        write(iUnit,*)'vertical_method = MEAN_SEA_LEVEL'
      case(entire_atmosphere_single_layer)
        write(iUnit,*)'vertical_method = ENTIRE_ATMOSPHERE_LAYER'
      !
      ! Multi-level verticals
      !
      case(constant_pressure)
        write(iUnit,*)'number_of_levels = ',vert%nLevs
        write(iUnit,*)'vertical_method = CUSTOM_LEVELS'
        write(iUnit,*)'level_type = PRESSURE'
        write(iUnit,'(1x,A9,200(1x,F9.2))')'levels = ',(vert%levs(iLev)%a,iLev=1,vert%nLevs)

      case(constant_height)
        write(iUnit,*)'number_of_levels = ',vert%nLevs
        write(iUnit,*)'vertical_method = CUSTOM_LEVELS'
        write(iUnit,*)'level_type = HEIGHT_FROM_SURFACE'
        write(iUnit,'(1x,A9,200(1x,F9.2))')'levels = ',(vert%levs(iLev)%a,iLev=1,vert%nLevs)

      case(constant_altitude)
        write(iUnit,*)'number_of_levels = ',vert%nLevs
        write(iUnit,*)'vertical_method = CUSTOM_LEVELS'
        write(iUnit,*)'level_type = ALTITUDE_FROM_SEA'
        write(iUnit,'(1x,A9,200(1x,F9.2))')'levels = ',(vert%levs(iLev)%a,iLev=1,vert%nLevs)

      case(sigma_level)
        write(iUnit,*)'number_of_levels = ',vert%nLevs
        write(iUnit,*)'vertical_method = CUSTOM_LEVELS'
        write(iUnit,*)'level_type = SIGMA'
        write(iUnit,'(1x,A9,200(1x,F9.4))')'levels = ',(vert%levs(iLev)%a,iLev=1,vert%nLevs)

      case(depth_level)
        write(iUnit,*)'number_of_levels = ',vert%nLevs
        write(iUnit,*)'vertical_method = CUSTOM_LEVELS'
        write(iUnit,*)'level_type = DEPTH'
        write(iUnit,'(1x,A9,200(1x,F9.4))')'levels = ',(vert%levs(iLev)%a,iLev=1,vert%nLevs)

      case(hybrid)
        write(iUnit,*)'number_of_levels = ',vert%nLevs
        write(iUnit,*)'vertical_method = CUSTOM_LEVELS'
        write(iUnit,*)'levels = HYBRID'
        write(iUnit,'(1x,A9,200(1x,I4))')'levels = ',(vert%levs(iLev)%number,iLev=1,vert%nLevs)
        do iLev = 1, vert%nLevs
          write(iUnit,*)'hybrid_coefficients = ',vert%levs(iLev)%number, &
                                               & vert%levs(iLev)%a, vert%levs(iLev)%b
        end do
      !
      ! Layers
      !
      case(layer_btw_2_pressure)
        write(iUnit,*)'number_of_levels = ',vert%nLevs
        write(iUnit,*)'vertical_method = CUSTOM_LAYERS'
        write(iUnit,*)'level_type = PRESSURE'
        write(iUnit,'(1x,A9,200(1x,F9.2))')'layer_thickness = ', &
                                   & (vert%levs(iLev)%a2 - vert%levs(iLev)%a,iLev=1,vert%nLevs)

      case(layer_btw_2_altitude)
        write(iUnit,*)'number_of_levels = ',vert%nLevs
        write(iUnit,*)'vertical_method = CUSTOM_LAYERS'
        write(iUnit,*)'level_type = ALTITUDE_FROM_SEA'
        write(iUnit,'(1x,A9,200(1x,F9.2))')'layer_thickness = ', &
                                   & (vert%levs(iLev)%a2 - vert%levs(iLev)%a,iLev=1,vert%nLevs)

      case(layer_btw_2_height)
        write(iUnit,*)'number_of_levels = ',vert%nLevs
        write(iUnit,*)'vertical_method = CUSTOM_LAYERS'
        write(iUnit,*)'level_type = HEIGHT_FROM_SURFACE'
        write(iUnit,'(1x,A,200(1x,F9.2))')'layer_thickness = ', &
                                   & (vert%levs(iLev)%a - vert%levs(iLev)%b,iLev=1,vert%nLevs)

      case(layer_btw_2_sigma)
        write(iUnit,*)'number_of_levels = ',vert%nLevs
        write(iUnit,*)'vertical_method = CUSTOM_LAYERS'
        write(iUnit,*)'level_type = SIGMA'
        write(iUnit,'(1x,A9,200(1x,F9.4))')'layer_thickness = ', &
                                   & (vert%levs(iLev)%a2 - vert%levs(iLev)%a,iLev=1,vert%nLevs)

      case(layer_btw_2_hybrid)
        write(iUnit,*)'number_of_levels = ',vert%nLevs
        write(iUnit,*)'vertical_method = CUSTOM_LAYERS'
        write(iUnit,*)'level_type = HYBRID'
        do iLev = 1, vert%nLevs
          write(iUnit,*)'hybrid_coefficients_bottom = ',vert%levs(iLev)%number2, &
                                                      & vert%levs(iLev)%a2, vert%levs(iLev)%b2
        end do
        write(iUnit,*)'hybrid_coefficients_domain_top = ', &
                                            & vert%levs(vert%nLevs)%a, vert%levs(vert%nLevs)%b

      case(layer_btw_2_depth)
        write(iUnit,*)'number_of_levels = ',vert%nLevs
        write(iUnit,*)'vertical_method = CUSTOM_LAYERS'
        write(iUnit,*)'level_type = DEPTH'
        write(iUnit,'(1x,A9,200(1x,F9.4))')'layer_thickness = ', &
                                   & (vert%levs(iLev)%a2 - vert%levs(iLev)%a,iLev=1,vert%nLevs)

      case default
        write(iUnit,*)'vertical_method = UNKNOWN METHOD'
        call report(vert)
        call set_error('Unknown vertical','report_vertical_as_namelist')

    end select

  end subroutine report_vertical_as_namelist


  ! ***************************************************************


  SUBROUTINE level_tests()

    IMPLICIT NONE
    TYPE(silja_level) :: level1, level2

    level1 = pr_level_700hpa
    CALL report(level1)

    level2 = fu_set_pressure_level(85000.)
    CALL report(level2)

    IF (level1 < level2) THEN
      call msg('2 korkeammalla')
    ELSE
      call msg('1 korkeammalla')
    END IF

  END SUBROUTINE level_tests

  subroutine test_overlap_fraction()
    implicit none
    real, dimension(5), parameter :: height_thin = (/0.0, 25.0, 50.0, 75., 100./)
    real, dimension(2), parameter :: height_thick = (/0.0, 1000./)
    real, dimension(2), parameter :: pressure = (/101324.89, 99298.39/)

    type(silja_level) :: top, bottom
    type(silja_level), dimension(5) :: levels
    type(silam_vertical) :: vert_hgt_thin, vert_hgt_thick, vert_press
    real :: fraction
    integer :: iz

    do iz = 1, 4
      bottom = fu_set_level(constant_height, fval1=height_thin(iz))
      top = fu_set_level(constant_height, fval1=height_thin(iz+1))
      levels(iz) = fu_set_layer_between_two(top, bottom)
    end do
    call set_vertical(levels, vert_hgt_thin)
    call arrange_levels_in_vertical(vert_hgt_thin)
    
    bottom = fu_set_level(constant_height, fval1=height_thick(1))
    top = fu_set_level(constant_height, fval1=height_thick(2))
    call set_vertical((/fu_set_layer_between_two(top, bottom)/), vert_hgt_thick)
    call arrange_levels_in_vertical(vert_hgt_thick)

    bottom = fu_set_level(constant_pressure, fval1=pressure(1))
    top = fu_set_level(constant_pressure, fval1=pressure(2))
    call set_vertical((/fu_set_layer_between_two(top, bottom)/), vert_press)
    call arrange_levels_in_vertical(vert_press)

    fraction = fu_vert_overlap_fraction(fu_level(vert_hgt_thin, 1), &
                                         & fu_level(vert_hgt_thick, 1))
    print *, 'fraction of hgt_thin(1) and hgt_thick(1)', fraction

    fraction = fu_vert_overlap_fraction(fu_level(vert_hgt_thick, 1), &
                                         & fu_level(vert_hgt_thin, 1))
    print *, 'fraction of hgt_thick(1) and hgt_thin(1)', fraction

    fraction = fu_vert_overlap_fraction(fu_level(vert_press, 1), &
                                         & fu_level(vert_hgt_thin, 1))
    print *, 'fraction of pressure(1) and hgt_thin(1)', fraction

    fraction = fu_vert_overlap_fraction(fu_level(vert_hgt_thick, 1), &
                                         & fu_level(vert_press, 1))
    print *, 'fraction of hgt_thick(1) and pressure(1)', fraction ! should ~170

  end subroutine test_overlap_fraction

  subroutine test_vert_to_metric()
    implicit none
    real, dimension(4) :: hyb_a = (/0.0, 0.0, 0.0, 783.0/), &
         & hyb_b = (/1.0, 0.992, 0.973, 0.929/), hyb_lev_press, lev_hgt
    type(silja_level), dimension(4) :: levels
    type(silja_level), dimension(3) :: layers
 
    character(len=*), parameter :: sub_name = 'test_vert_to_metric'
    type(silja_level) :: top, bottom
    type(silam_vertical) :: vert_hyb_thick, vert_m, vert_press
    real :: height_from_vert
    integer :: iz

    hyb_lev_press = hyb_a + std_pressure_sl*hyb_b
    lev_hgt = (/(fu_height_for_press(hyb_lev_press(iz)), iz=1, 4)/)

    ! hybrid to height, layers
    do iz = 1, 3
      bottom = fu_set_level(hybrid, hyb_a(iz), hyb_b(iz))
      top = fu_set_level(hybrid, hyb_a(iz+1), hyb_b(iz+1))
      layers(iz) = fu_set_layer_between_two(bottom, top)
    end do
    call set_vertical(layers, vert_hyb_thick)
    call arrange_levels_in_vertical(vert_hyb_thick)
    
    call msg('*** hybrid thick layers to metric ***')
    call vert_to_metric(vert_hyb_thick, vert_m)
    do iz = 1, 3
      bottom = fu_lower_boundary_of_layer(fu_level(vert_m, iz))
      height_from_vert = fu_level_height(bottom)
      if (fu_fails(height_from_vert .eps. lev_hgt(iz), 'Heights differ', sub_name)) continue
      call msg('heights (lower):', height_from_vert, lev_hgt(iz))
      height_from_vert = fu_level_height(fu_level(vert_m,iz,.true.))
      if (fu_fails(height_from_vert .eps. 0.5*(lev_hgt(iz) + lev_hgt(iz+1)), 'Heights differ', sub_name)) continue
      call msg('heights (mid):', height_from_vert, 0.5*(lev_hgt(iz) + lev_hgt(iz+1)))
    end do
    
    ! pressure to height, levels
    
    do iz = 1, 4
      levels(iz) = fu_set_level(constant_pressure, hyb_lev_press(iz))
    end do
    call set_vertical(levels, vert_press)
    call arrange_levels_in_vertical(vert_press)
    call msg('*** pressure levels to metric ***')
    call vert_to_metric(vert_press, vert_m)
    do iz = 1, 4
      height_from_vert = fu_level_height(fu_level(vert_m, iz))
      if (fu_fails(height_from_vert .eps. lev_hgt(iz), 'Heights differ', sub_name)) continue
      call msg('heights:', height_from_vert, lev_hgt(iz))
    end do
    
  end subroutine test_vert_to_metric

  subroutine test_vert_to_pressure()
    implicit none
    real, dimension(4) :: hyb_a = (/0.0, 0.0, 0.0, 783.0/), &
         & hyb_b = (/1.0, 0.992, 0.973, 0.929/), hyb_lev_press, lev_hgt
    type(silja_level), dimension(4) :: levels
    type(silja_level), dimension(3) :: layers
    type(silam_vertical) :: vert_hyb_thick, vert_m, vert_press
    type(silja_level) :: top, bottom
    character(len=*), parameter :: sub_name = 'test_vert_to_pressure'
    integer :: iz

    lev_hgt = (/0.0, 1500.0, 25000.0, 50000.0/)
    
    ! height to pressure, levels
    do iz = 1, 4
      levels(iz) = fu_set_level(constant_height, lev_hgt(iz))
    end do
    call set_vertical(levels, vert_m)
    call arrange_levels_in_vertical(vert_m)
    call msg('*** metric to pressure ***')
    call vert_to_pressure(vert_m, vert_press)
    call report(vert_press, .true.)
    
    ! hybrid to pressure, layers
    do iz = 1, 3
      bottom = fu_set_level(hybrid, hyb_a(iz), hyb_b(iz))
      top = fu_set_level(hybrid, hyb_a(iz+1), hyb_b(iz+1))
      call msg('Bottom', iz, fu_hybrid_level_pressure(bottom, std_pressure_sl))
      call msg('Top', iz, fu_hybrid_level_pressure(top, std_pressure_sl))
      layers(iz) = fu_set_layer_between_two(bottom, top)
    end do
    call set_vertical(layers, vert_hyb_thick)
    call arrange_levels_in_vertical(vert_hyb_thick)
    call msg('*** hybrid to pressure ***')
    call vert_to_pressure(vert_hyb_thick, vert_press)
    call report(vert_press, .true.)

  end subroutine test_vert_to_pressure

END MODULE silam_levels
