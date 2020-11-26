MODULE natural_constants

  ! Author: Mika Salonoja, email Mika.Salonoja@fmi.fi
  !
  ! All units: SI
  !
  ! Language: ANSI standard Fortran 90
  !
  USE globals, only : r4k, r8k
  IMPLICIT NONE

  ! general constants
  REAL, PARAMETER, PUBLIC :: g = 9.815
  REAL, PARAMETER, PUBLIC :: e = 2.718281828
  real, parameter, public :: e_1 = 1/2.718281828
  REAL, PARAMETER, PUBLIC :: pi = 3.141592654
  real, parameter, public :: ln_2 = 0.6931471805599453
  real, parameter, public :: ln_10 = 2.302585093
  real, parameter, public :: ln_10_1 = 0.43429448190
  real(r8k), parameter, public :: dsqrt_2 = 1.4142135623730950488016887242097

  real, parameter, public :: eps_degrees = 0.005 !!limit of single-precision lat-lon coordinates
         !! ~50m 
 
  REAL(8), PARAMETER, PUBLIC :: d_pi =  4*ATAN(1.d0)
  REAL(8), PARAMETER, PUBLIC :: drad_to_deg = 180.d0/d_pi
  REAL(8), PARAMETER, PUBLIC :: ddeg_to_rad = d_pi/180.d0

  real, parameter, public :: MAX_REAL = 3.4e38
  real, parameter, public :: LOG_MAX_REAL = log(MAX_REAL)

  REAL, PARAMETER, PUBLIC :: radians_to_degrees =  57.29577951
  REAL, PARAMETER, PUBLIC :: degrees_to_radians =  0.01745329252
  REAL, PARAMETER, PUBLIC :: karmann_c = 0.4 ! karmann-constant
  REAL, PARAMETER, PUBLIC :: zero_celcius = 273.15 ! 0C in Kelvins
  REAL, PARAMETER, PUBLIC :: temperature_ref = 293.15 ! reference temperature (=20C in K)
  REAL, PARAMETER, PUBLIC :: earth_radius = 6378000.0  ! meters
       ! CRC Handbook gives 2 values: 6378245 and 6378077.
  REAL, PARAMETER, PUBLIC :: earth_omega = 7.292E-5 ! earth rotation, rad/s
  REAL, PARAMETER, PUBLIC :: equator_length = 2.0 * pi * earth_radius
  
  ! basic values for constants that could be computed, too
  REAL, PARAMETER, PUBLIC :: basic_coriolis = 1.0E-4 ! J/kg
  REAL, PARAMETER, PUBLIC :: density_air_288K = 1.225 ! for dry air at 15C and 1013.25 hPa, kg/m**3   
  REAL, PARAMETER, PUBLIC :: density_std  = 1.29 ! for dry air at 0C and 1013.25 hPa, kg/m**3   
  REAL, PARAMETER, PUBLIC :: density_water = 1025. ! for moist air kg/m**3 
  REAL, PARAMETER, PUBLIC :: density_ice = 917. ! for ice at 0C in kg/m**3 
  REAL, PARAMETER, PUBLIC :: molecular_weight_air = 28.97e-3 ! ma, molecular weight for dry air, kg/mol
  REAL, PARAMETER, PUBLIC :: molecular_weight_water = 18.015e-3  ! [kg/mol]
  real, parameter, public :: srfTns_water = 0.073   ! surface tension of pure water @ ~ 293 K [J/m2]

  ! solar energy related constants
  REAL, PARAMETER, PUBLIC :: solar_const = 1370 !in W/m**2
  REAL, PARAMETER, PUBLIC :: solar_const_SI = 1.113 ! in Km/s
  REAL, PARAMETER, PUBLIC :: velocity_light = 2.99792458E8 ! in m/s
  REAL, PARAMETER, PUBLIC :: stefan_boltzmann = 5.67E-8 ! Stefan Boltzmann in W/(m**2K**4)
  real, parameter, public :: boltzmann_const = 1.38E-23 ! J/K
  
  ! atmospheric thermodynamics
  REAL, PARAMETER, PUBLIC :: std_pressure_sl = 101325. ! Pa - standard sea level pressure
  REAL, PARAMETER, PUBLIC :: std_height_scale = 7400. ! m, for barostatic eq
  REAL, PARAMETER, PUBLIC :: specific_heat_dryair = 1005.7 ! cpa in J/kgK
  REAL, PARAMETER, PUBLIC :: specific_heat_watervapour = 1875.0 ! cpw in J/kgK
  REAL, PARAMETER, PUBLIC :: specific_heat_liquidwater = 4200.0 ! cplw in J/kgK
  REAL, PARAMETER, PUBLIC :: specific_heat_ice = 2100.0 ! cpi in J/kgK
  
  REAL, PARAMETER, PUBLIC :: gas_constant_uni = 8.3143 ! universal R in J/molK
  REAL, PARAMETER, PUBLIC :: gas_constant_dryair = 287.04 ! Ra = R/ma in J/kg*K
  REAL, PARAMETER, PUBLIC :: gas_constant_watervapour = 461.5 ! Rw in J/kg*K
  REAL, PARAMETER, PUBLIC :: gas_constant_ratio = 0.622 ! Ra/Rw
  
  REAL, PARAMETER, PUBLIC :: fusion_latentheat = 3.34E5 ! J/kg
  REAL, PARAMETER, PUBLIC :: condensation_latentheat = 2.5E6 ! J/kg
  REAL, PARAMETER, PUBLIC :: sublimation_latentheat = 2.83E6 ! J/kg
  REAL, PARAMETER, PUBLIC :: vaporization_latentheat = 2.50E6 ! at 0C in J/kg

  REAL, PARAMETER, PUBLIC :: dry_adiabatic = -0.0098 ! dry adiabatic lapse rate, K/m
  REAL, PARAMETER, PUBLIC :: pseudo_adiabatic = -0.0065 ! pseudo adiabatic lapse rate, K/m

  ! aerodynamical properties of particles
  real, parameter, public :: avogadro = 6.025E23 ! Avogadro number
  REAL, PARAMETER, PUBLIC :: dynamic_viscosity = 1.81E-5 ! mean value in kg/ms**2
  REAL, PARAMETER, PUBLIC :: kinematic_viscosity = 1.461E-5 ! mean value in  m**2/s
  REAL, PARAMETER, PUBLIC :: free_path_mol = 6.53E-8 ! mean free path for molecules in air, m

  ! specials
  REAL, PARAMETER, PUBLIC :: R_per_c_dryair = 0.2854  ! = gas_constant_dryair/specific_heat_dryair
  REAL, PARAMETER, PUBLIC :: R_per_c_watervapour = 0.246  ! = gas_constant_watervapour/specific_heat_watervapour
  REAL, PARAMETER, PUBLIC :: heatcapacity_water = 4.295E6 ! Wm**2/(Kms)  ! = density_water*specific_heat_water 
  real, parameter, public :: prandtl_nbr = 0.72
  REAL, PARAMETER, PUBLIC :: small_epsilon  = 1.2E-38 ! a very small number, indeed
  REAL, PARAMETER, PUBLIC :: mev_in_joule = 1.6021773349e-13  ! CRC Handbook
  real, parameter, public :: max_wind_speed = 200.0           ! m/s
  real, parameter, public :: SWRad_2_PAR = 2.7e-3 * 0.45 *4.6 ! PAR=0.45*SWR(W/m2)*4.6(microphotons/s per W)
  real, parameter, public :: mode_diam_tolerance = 1e-8 ! 10nm

  real, parameter, public :: cloud_saturation = 1.5e-4 !!!kg/kg water content, above which clouds produce rain
                                                !! 3e-4 -- reported in literatire
  

  real, public, parameter :: Kz_FT_factor = 10. ! A ratio of Kz in ABL and free troposphere 
!  real, public, parameter :: Kz_FT_factor = 10000. ! A ratio of Kz in ABL and free troposphere 

END MODULE natural_constants

