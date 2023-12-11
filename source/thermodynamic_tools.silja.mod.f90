MODULE thermodynamic_tools
  !
  ! This module contains general tools for conversions between
  ! meteorological humidity and temperature quantites as well as
  ! computing of related thermodynamical parameters
  ! Some general functions (i.g. "constant" for gravitational
  ! acceleration ) are also included.
  !
  ! Author: Ilkka Valkama, FMI email Ilkka.Valkama@fmi.fi
  !
  ! All units: SI
  ! 
  ! All degrees: real numbers in degrees (degrees and parts per
  ! hundred/parts per thousand , NOT in degrees and minutes
  ! 
  ! Language: ANSI standard Fortran 90
  ! 
  !
  USE silam_namelist !toolbox !work_arrays
!  USE natural_constants

  !$use omp_lib
 
  IMPLICIT NONE

  ! The public functions and subroutines available in this module:
!!!$  PUBLIC fu_compute_g
  !      123456789012345678901234567890
  PUBLIC component_wind
  PUBLIC normal_wind
  PUBLIC fu_garratt_acceleration
  PUBLIC fu_latentheat_of_vaporization
  PUBLIC fu_relhum_from_dew_point
  PUBLIC fu_dew_point_from_rh
  PUBLIC fu_dew_point_from_temp
  PUBLIC fu_lifting_condensation_level
  PUBLIC fu_virtual_temp_from_mixratio  
  PUBLIC fu_virtual_temp_from_spechumi
  PUBLIC fu_virtual_temp_from_relhumi
  PUBLIC fu_gasconstant_of_water_vapour
  PUBLIC fu_specificheat_of_watervapour
  PUBLIC fu_specificheat_of_air
  PUBLIC fu_saturat_watervapourpressure
  PUBLIC fu_air_density_profile
  PUBLIC fu_watervapour_densityprofile
  PUBLIC fu_dynamic_viscosity
  PUBLIC fu_kinematic_viscosity
  PUBLIC fu_watervap_diffusivity
  public fu_air_diffusivity
  public fu_R_a
  public us_standard_atmosphere
  public fu_pr_diff_height
  public fu_height_diff_between_pre
  public std_atm_rho_t
  PUBLIC fu_dew_point_temperature
  PUBLIC fu_virtual_temperature
  PUBLIC fu_virtual_temp_from_spec_hum 
  PUBLIC fu_virtual_temp_from_dewpoint
  PUBLIC fu_potential_temperature
  PUBLIC fu_temperature_from_pot_temp
  PUBLIC fu_equivalent_temperature
  PUBLIC fu_equivalent_pot_temperature
  PUBLIC fu_wet_bulb_temperature
  PUBLIC fu_wet_bulb_pot_temperature
  PUBLIC fu_moist_potential_temperature 
  PUBLIC fu_lcl_temperature_from_rel_hum
  PUBLIC fu_lcl_temp_from_dew_point
  PUBLIC fu_lcl_pressure
  PUBLIC fu_pseudoadiab_equiv_pot_temp
  PUBLIC fu_watervapour_satur_p
  PUBLIC fu_vapor_pressure_dew_point
  PUBLIC spechum_to_relhum_field
  PUBLIC fu_spechum_to_relhum
  public fu_wind_gust_10m
  public fu_effective_wind_pwr_gust
  public fu_fade_in
  public fu_fade_out
  public set_fade_in_out_params
  public report

    ! Private types defined in this module:
  PRIVATE vapor_pressure_dewp
  private psim_prof
  private psih_prof
  private thermo_tests
  private temp_test
  private fu_wind_gust_default_10m
  private fu_wind_gust_general_10m
  private make_gust_lookup_table
  private report_fade_in_out_setup

  interface fu_wind_gust_10m
    module procedure fu_wind_gust_default_10m
    module procedure fu_wind_gust_general_10m
  end interface

  interface report
    module procedure  report_fade_in_out_setup
  end interface 
  

  !
  ! Wind gustinness seems to be needed for several applications. Let's make it properly.
  !
  integer, private, parameter :: nWinds_glob = 200, nZ_glob = 100, nWindPwr_glob = 10
  type Twind_gust_lookup
!    private
    integer :: nWinds=nWinds_glob, nZ=nZ_glob, nWindPwr=0
    real, dimension(:), pointer :: fWind, fZr, fWindPwr
    real, dimension(:,:,:), pointer :: pGustVal  !(nWinds, nZ_z0, nWindPower)
  end type Twind_gust_lookup
  public Twind_gust_lookup

  type(Twind_gust_lookup), private, save :: wgl
  !
  ! Tabulated US Standard Atmosphere
  !
  logical, private, save :: ifAtmosphereTabulated = .false.
  integer, private, parameter :: nStandardValues = 750
  real, private, parameter :: delta_standard_table = 100.
  real, dimension(nStandardValues), private :: arStandardDensity, arStandardPressure, &
                                             & arStandardTemperature
  !
  ! Two functions for smooth transition of whatever variable: 
  ! fade-in/out functions. Here are their features 
  !
  type Tfade_inout
    integer :: iTypeIn, iTypeOut
    real :: fAsymmetryIn, fExtraShiftIn, fAsymmetryOut, fExtraShiftOut
  end type Tfade_inout
 
  type (Tfade_inout), public, parameter ::  fade_inout_missing = &
       & Tfade_inout(int_missing, int_missing, &
       & real_missing, real_missing, real_missing, real_missing)

  
CONTAINS


  ! ****************************************************************

  REAL FUNCTION fu_garratt_acceleration(latitude,altitude)
    
    ! Description:
    ! Returns gravitational acceleration for a geographical position
    ! at or above mean sea level (altitude in metres). 
    !
    ! Method:
    ! The value for g is computed as a function  of latitude and
    ! altitude is computed as presented in Garratt (1992).
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Ilkka Valkama, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in) ::  altitude, latitude
    
    ! local declarations
    REAL :: lat_in_radians, ground_value
    REAL :: local_alfa, local_beta
    REAL :: altitude_correction

    ! The mean raius of the earth 
    ! r(equator)= 6378.137 km 
    ! r(pole)= 6356.75 km)
    !-------------------------    

    ! First trigonometrics
    !  

    local_alfa = (1.0 - COS(2.0*(degrees_to_radians*latitude)))/2.0
    ground_value = 0.0053*local_alfa

    ! Second trigonometrics

    local_alfa = (1.0 - COS(4.0*(degrees_to_radians*latitude)))/2.0
    ground_value = ground_value - 5.8E-6*local_alfa

    ! ground level value

    ground_value =  9.78033*(1.0 + ground_value)
 

    ! At higher levels when necessary

    IF(altitude > 0.0)THEN
	altitude_correction= 1.0 + altitude/earth_radius
	fu_garratt_acceleration=ground_value/(altitude_correction&
	    & *altitude_correction)
    ELSE
	fu_garratt_acceleration=ground_value
    END IF

    END FUNCTION fu_garratt_acceleration


 ! ****************************************************************

    SUBROUTINE component_wind(direction,force,&
        &  zonal_component, meridional_component)
    
      ! Description:
      ! Returns the zonal (u) and meridional (v) components for a given 
      ! wind vector (direction in rad and speed in m/s). 
      !
      ! Method: Based on the definition of vector wind.
      !
      ! All units: SI
      !
      ! Language: ANSI Fortran 90
      !
      ! Author: Ilkka Valkama, FMI
      
      IMPLICIT NONE
      
      ! Imported parameters with intent(in):
      REAL, INTENT(in) :: direction,force
      ! Imported parameters with intent(out):
      REAL, INTENT(out) :: zonal_component, meridional_component
      
      ! local declarations
      REAL :: vario
 
      IF(direction > 0.0) THEN
        ! testing with direction in degrees
!!!$        vario = direction/57.29578
        vario = direction
      ELSE
        vario = 0.0
      END IF
      IF(force > 0.0 )THEN
        meridional_component = -1.0*force
        zonal_component = meridional_component*SIN(vario)	
        meridional_component = meridional_component*COS(vario)
      ELSE
        meridional_component = 0.0
        zonal_component =  0.0
      END IF
      
    END SUBROUTINE component_wind


    ! ****************************************************************

    SUBROUTINE normal_wind(zonal_component, meridional_component,  direction,force)
    
      ! Description:
      ! Returns the synoptic wind  (direction in rad and speed in m/s)
      ! from the given wind vector components (u and v). 
      !
      ! Method: Based on the definition of vector wind.
      !
      ! All units: SI
      !
      ! Language: ANSI Fortran 90
      !
      ! Author: Ilkka Valkama, FMI
      
      IMPLICIT NONE
      
      ! Imported parameters with intent(in):
      REAL, INTENT(in) ::  zonal_component, meridional_component
      ! Imported parameters with intent(out):
      REAL, INTENT(out) :: direction,force
      
      ! local declarations
      REAL :: vario
      
      vario =  zonal_component*zonal_component +  meridional_component*meridional_component
 
      IF(vario > 0.0)THEN
        force = SQRT(vario)
        direction = ATAN2(zonal_component,meridional_component) + pi
      else
        force = 0.0    
        direction = 0.0
      END IF
      
    END SUBROUTINE normal_wind


  ! ****************************************************************

  REAL FUNCTION fu_saturat_watervapourpressure(temperature) 
    ! Description:
    !  Returns the saturation vapour pressure for a temperature
    !  valid for RANGE -40 C to +50 C
    !
    ! Method:
    ! The original Goff-Gratch-formulation
    ! ( in its  1946 modification). Function originally published in :
    ! 
    !
    ! All units: NOT SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Ilkka Valkama, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: temperature
 
    ! Local declarations :

    REAL local_temperature1, local_temperature2, local_temperature3
    REAL water_vapour,local_vapour1,local_vapour2

    !----------------------------------------
    !
    ! 1. Check the input  temperature.
    !    ------------------------------------
    
		
    IF(temperature.lt.100)THEN
      ! Only kelvin is acceptable
      local_temperature1=temperature+273.16
    ELSE
      local_temperature1=temperature
    END IF

    !----------------------------------------
    !
    ! 2. Start computing 
    !    ----------------
    
	
    local_temperature2= local_temperature1/273.16
    local_temperature3= 1./local_temperature2

    local_vapour1= 8.29692*(1.0 - local_temperature2)                   
    local_vapour2= 4.76955*(1.0 - local_temperature3)                 
    local_vapour1= 1.0 - 10.**local_vapour1
    local_vapour2= 10.**local_vapour2 - 1.

    water_vapour= 10.79574*(1.0 - local_temperature3) - 5.0280&
	&*(LOG10(local_temperature2))
    water_vapour=water_vapour + 1.50475E-4*local_vapour1 + 0.42873E-3&
	&*local_vapour2

    fu_saturat_watervapourpressure = 10.**(water_vapour + 0.78614)  ! in hPa

  END FUNCTION fu_saturat_watervapourpressure


 ! ****************************************************************

  REAL FUNCTION fu_latentheat_of_vaporization(temperature)
    
    ! Description:
    ! Returns the latent heat of vaporization (J/kg) for a given
    ! temperature 
    !
    ! Method:
    !   For temperatures given in Celsius the latent heat of 
    !    vaporization/sublimation is computed separately for 
    !    water and ice phases as in Andreas et al. (1996) :
    !   For temperture given ibn kelvin a common equation from 
    !    Bolton (1980) is used
    !
    ! If temperature is missing a constant value is returned
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Ilkka Valkama, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in), OPTIONAL :: temperature

    !----------------------------------------
    !
    ! 1. If no temperature given, use constant value at 0 degrees Celsiuss .
    !    ------------------------------------
    
    IF (.NOT.PRESENT(temperature)) THEN
      fu_latentheat_of_vaporization = 2.501E6
      RETURN
    END IF


    !----------------------------------------
    !
    ! 2. Calculate latent heat in J/kg
    !    -----------------------

    IF(temperature < 100.)THEN
      IF(temperature < 0.)THEN
      fu_latentheat_of_vaporization = (28.34 - 0.00149&
	  &*temperature)*1E5
      ELSE  
      fu_latentheat_of_vaporization = (25.0 - 0.002274&
	  &*temperature)*1E5
      END IF
    ELSE
      fu_latentheat_of_vaporization = (2.501 - 0.00237&
	  &*temperature)*1E6
    END IF
    
  END FUNCTION fu_latentheat_of_vaporization


 ! ****************************************************************

  REAL FUNCTION fu_relhum_from_dew_point(temperature,&
      & dewpoint_temperature) result(relative_humidity)
     

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: temperature,  dewpoint_temperature

    real :: T, Dp

    !----------------------------------------
    !
    ! 1. NB the rel_hum is given as a fraction (not %)
    !    ------------------------------------
    
    ! Iriginal code by Valkama, with reference to Babin 1995                                                                        |  --------------------------------------------------------------------------------------------------------------------
    ! ( adopted from Bras, 1990 )
    ! Seems to be wrong
    !ambient_relhum = 112. - 0.1*temperature +  dewpoint_temperature
    !ambient_relhum = ambient_relhum/(112. + 0.9*temperature)
    !relative_humidity = ambient_relhum**8.0
     
    !Expression from 
    ! https://www.omnicalculator.com/physics/relative-humidity
    ! Who refers to  Alduchov and Eskridge, 1996
    ! https://doi.org/10.1175/1520-0450(1996)035<0601:IMFAOS>2.0.CO;2

     T = temperature + zero_celcius
     Dp = dewpoint_temperature + zero_celcius

      relative_humidity = exp(17.625 * Dp/(243.04 + Dp))/exp(17.625 * T/(243.04 + T)) 
    
  END FUNCTION fu_relhum_from_dew_point


 ! ****************************************************************

  REAL FUNCTION fu_dew_point_from_temp(temperature) result(dewpoint)
    ! Description:
    ! Returns the dew point temperature for a given air temperature
    !
    ! Method: Subroutine copied from Sakari Kajosaari (c. 1980)
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Ilkka Valkama, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: temperature

    ! Local declarations
    REAL ::  sat_pres, ambient_temperature, ambient_relhum, ambient_dewp

    !----------------------------------------
    !
    ! 1. Check the input data &  Scale the units 
    !    ------------------------------------
    
    IF(temperature < 100.)THEN
      ambient_temperature = temperature + 273.16
    ELSE
      ambient_temperature = temperature
    END IF
    sat_pres = fu_saturat_watervapourpressure(temperature)
  
    ambient_relhum = LOG10((sat_pres/6.11))
 
    IF(sat_pres > 4.015) THEN
      IF(sat_pres < 6.108) THEN
        ambient_dewp = 1.282 - ambient_relhum
        ambient_relhum = 40.18*ambient_relhum
      ELSE
        ambient_dewp = 7.772 - ambient_relhum
        ambient_relhum = 244.65*ambient_relhum
      END IF
    ELSE
      ambient_dewp = 7.772 - ambient_relhum
      ambient_relhum = 209.06*ambient_relhum
    END IF

    dewpoint = ambient_relhum/(ambient_dewp + 273.16)
    
  END FUNCTION fu_dew_point_from_temp


 ! ****************************************************************

  REAL FUNCTION fu_dew_point_from_rh(temperature, relative_humidity) result(dewpoint)
     
    ! Description:
    ! Returns the dew point temperature for a given air temperature
    ! and relative humidity
    !
    ! Method: Subroutine copied from Kalle Eerola and NOT VERIFIED !!!
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Ilkka Valkama, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: temperature, relative_humidity

    ! Local declarations
    REAL ::  ambient_temperature, ambient_relhum, ambient_dewp

    !----------------------------------------
    !
    ! 1. Check the input data &  Scale the units 
    !    ------------------------------------
    

    IF(temperature < 100.)THEN
      ambient_temperature = temperature + 273.16
    ELSE
      ambient_temperature = temperature
    END IF

    IF(relative_humidity < 1.0)THEN
      ambient_relhum = MIN(relative_humidity,1.0)
    ELSE  
      ambient_relhum = MIN((relative_humidity/100.),1.0) 
    END IF

    !----------------------------------------
    !  ! 3. Calculate dew point 
    !    -----------------------

    ambient_relhum = MAX(ambient_relhum,0.01)
           
     ambient_dewp = ambient_temperature/&
	&(1. - ambient_temperature*1.8456591E-4&
	&*LOG(ambient_relhum))

    dewpoint = MIN(ambient_temperature, ambient_dewp)     
  
  END FUNCTION fu_dew_point_from_rh


 ! ****************************************************************

  REAL FUNCTION fu_lifting_condensation_level(temperature, relative_humidity,dew_point)
    
    ! Description:
    ! Returns the local lifting condensation level height
    !
    ! Method: Iribane & Godson (1981)
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Ilkka Valkama, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: temperature, relative_humidity
    REAL, INTENT(in), OPTIONAL :: dew_point

    ! Local declarations
    REAL ::  local_dewpoint, logarithmic_relhum, dewpoint_depression

    !----------------------------------------
    !
    ! 1. Check the input 
    !    ------------------------------------

    IF (PRESENT(dew_point)) THEN
      local_dewpoint = dew_point
    ELSE
      local_dewpoint = fu_dew_point_from_rh(temperature,&
	  & relative_humidity)
    END IF

    !----------------------------------------
    !
    ! 2. Compute the LCL
    !    ------------------------------------
    
    logarithmic_relhum = ABS(LOG10(relative_humidity))
    dewpoint_depression = 4.25E-4*temperature*local_dewpoint&
	&*logarithmic_relhum

    fu_lifting_condensation_level = 120.*dewpoint_depression

  END FUNCTION fu_lifting_condensation_level


  ! ****************************************************************

  REAL FUNCTION fu_virtual_temp_from_mixratio(temperature, mixing_ratio)                                           

    ! Description:
    ! Returns the virtual temperature from air temperature and mixing
    ! ratio
    !
    ! Method:   Using the definition of Tv
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Ilkka Valkama, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: temperature, mixing_ratio

    ! Local declarations
    REAL :: local_humidity

    !----------------------------------------
    !
    ! 1. Calculate the virtual temperature
    !    ------------------------------------

    local_humidity = (1.0 + 1.61*mixing_ratio)/(1.0 + mixing_ratio)

    fu_virtual_temp_from_mixratio = temperature&
	& *local_humidity

  END FUNCTION fu_virtual_temp_from_mixratio


 ! ****************************************************************

  REAL FUNCTION fu_virtual_temp_from_spechumi(temperature, specific_humidity)

    ! Description:
    ! Returns the virtual temperature from air temperature and
    ! specific humidity
    !
    ! Method:   Using the definition of Tv
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Ilkka Valkama, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: temperature, specific_humidity
 
    ! Local declarations

    REAL :: local_humidity

    !----------------------------------------
    !
    ! 1. Calculate the virtual temperature
    !    ------------------------------------

    local_humidity = (0.62198 + specific_humidity)/(0.62198*&
	& (1.0 + specific_humidity))

    fu_virtual_temp_from_spechumi = temperature*local_humidity

  END FUNCTION fu_virtual_temp_from_spechumi


  ! ****************************************************************

  REAL FUNCTION fu_virtual_temp_from_relhumi(pressure, temperature, relative_humidity)                                           

    ! Description:
    ! Returns the virtual temperature from air temperature and
    ! relative humidity
    !
    ! Method:   Using the definition of Tvirt from 
    !           WMO International Meteorological Tables
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Ilkka Valkama, FMI
    ! (original code by Sakari Kajosaari, c.1980)
    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: pressure, temperature,  relative_humidity
 
    ! Local declarations

    REAL :: local_humidity, sat_pres

    !----------------------------------------
    !
    ! 1. Calculate the virtual temperature
    !    ------------------------------------
    
    sat_pres =  fu_saturat_watervapourpressure(temperature)
    sat_pres = relative_humidity*sat_pres
    local_humidity = 0.62198*sat_pres/(pressure - sat_pres)

    fu_virtual_temp_from_relhumi = (temperature*(1.0 + &
	& local_humidity/0.62198))/(1.0 + local_humidity)

  END FUNCTION fu_virtual_temp_from_relhumi


  ! ****************************************************************

  REAL FUNCTION fu_gasconstant_of_water_vapour(mixing_ratio)                                           
    ! Description:
    !  Computes the gas constant for wet, but non-saturated air
    !  from the mixing ratio
    !
    ! Method:   If mixing ratio missing the value for saturated air
    ! is returned
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Ilkka Valkama, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in), OPTIONAL :: mixing_ratio

 
    !----------------------------------------
    !
    ! 1. Check the method
    !    ------------------------------------
    

    IF (.NOT.PRESENT(mixing_ratio)) THEN
      fu_gasconstant_of_water_vapour =  gas_constant_watervapour
      RETURN
    END IF


    !----------------------------------------
    !
    ! 2. Calculate the gas constant
    !    ------------------------------------

    IF(mixing_ratio < 1.0)THEN
    fu_gasconstant_of_water_vapour = gas_constant_dryair*(1.0 +&
	& 0.608E-3*mixing_ratio)
  ELSE
    fu_gasconstant_of_water_vapour =  gas_constant_watervapour
  END IF

  END FUNCTION fu_gasconstant_of_water_vapour

  ! ****************************************************************

  REAL FUNCTION fu_specificheat_of_watervapour(mixing_ratio)

    ! Description:
    !  Computes the value of specific heat for wet, but non-saturated air
    !  from the mixing ratio
    !
    ! Method:   If mixing ratio missing the value for saturated air
    ! is returned
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Ilkka Valkama, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in), OPTIONAL :: mixing_ratio

 
    !----------------------------------------
    !
    ! 1. Check the method
    !    ------------------------------------
    

    IF (.NOT.PRESENT(mixing_ratio)) THEN
      fu_specificheat_of_watervapour =   specific_heat_watervapour
      RETURN
    END IF


    !----------------------------------------
    !
    ! 2. Calculate the gas constant
    !    ------------------------------------

    IF(mixing_ratio < 1.0)THEN
      fu_specificheat_of_watervapour =  specific_heat_watervapour*(1.0 +&
	  & 0.877E-3*mixing_ratio)
    ELSE
      fu_specificheat_of_watervapour =   specific_heat_watervapour
    END IF

  END FUNCTION fu_specificheat_of_watervapour

  ! ****************************************************************

  REAL FUNCTION fu_specificheat_of_air(temperature)

    ! Description:
    !  Computes the value of specific heat for dry air (in J/kgK) 
    !  as a function of air temperature
    !
    ! Method:   If mixing ratio missing a constant value is returned
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Ilkka Valkama, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in), OPTIONAL :: temperature

    ! local declarations
    REAL :: local_constant

    !----------------------------------------
    !
    ! 1. Check the method
    !    ------------------------------------
    

    IF (.NOT.PRESENT(temperature)) THEN
      fu_specificheat_of_air =   specific_heat_dryair
      RETURN
    END IF


    !----------------------------------------
    !
    ! 2. Calculate the specific heat "constant"
    !    ------------------------------------

    local_constant = temperature - 250.
    fu_specificheat_of_air = 1005.7 +&
	  & ((local_constant*local_constant)/3364.)


  END FUNCTION fu_specificheat_of_air

  ! ****************************************************************

  REAL FUNCTION fu_air_density_profile(temperature, pressure)

    ! Description:
    !  Computes the "real" value of air density from  air temperature
    !  and pressure
    !
    ! Method:
    !   The general formula (see e.q. Andreas et al.1996) is used 
    !      density = 1.2923(273.156/T)(p/1013.23)
    !   where T is in K and p is in hPa
    !   To speed up the constants have been multiplied into one
    !   single  constant
    !
    !  Should either T or p is missing a constant value is returned
    !
    ! All units: SI, air density = kg/m**3
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Ilkka Valkama, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in), OPTIONAL :: temperature, pressure
    REAL :: local_pressure
    !----------------------------------------
    !
    ! 1. Check the method
    !    ------------------------------------
    

    IF (.NOT.PRESENT(temperature)) THEN
      fu_air_density_profile = density_air_288K
      RETURN
    END IF
    IF (.NOT.PRESENT(pressure)) THEN
      fu_air_density_profile = density_air_288K
      RETURN
    END IF

    !----------------------------------------
    !
    ! 2. Calculate the air density
    !    ------------------------------------

    local_pressure = pressure/100. ! Pascal to hectoPascal
    fu_air_density_profile = 0.34838*(local_pressure/temperature)


  END FUNCTION fu_air_density_profile

  ! ****************************************************************

  REAL FUNCTION fu_watervapour_densityprofile(temperature, saturation_pressure)

    ! Description:
    !  Computes the "real" value of density of water vapour
    !  present using ambient air temperature and saturation pressure 
    !  in atmosphere 
    !
    ! Method:
    !   The general formula (see e.q. Andreas et al.1996) is used 
    !      water vapour density = 100.(Mw*es(t))(R*T)
    !   where T is in K and es(t) is in hPa.
    !   To speed up the valuess for molecular weight of water
    !   (Mw = 18.016X10-3 kg/mol) and the universal gas constant 
    !   (R = 8.31441 J/K*mol) have been multiplied into one
    !   single  constant
    !
    !  Should the saturation pressure be missing it is computed
    !  If temperature is missing a constant value is returned
    !
    ! All units: SI, air density = kg/m**3
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Ilkka Valkama, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in), OPTIONAL :: temperature, saturation_pressure

    ! Local declarations
    REAL :: local_vapour

    !----------------------------------------
    !
    ! 1. Check the method
    !    ------------------------------------
    

    IF (.NOT.PRESENT(temperature)) THEN
      fu_watervapour_densityprofile = density_air_288K
      RETURN
    END IF
    IF (.NOT.PRESENT(saturation_pressure)) THEN
      local_vapour =  fu_saturat_watervapourpressure(temperature)
    ELSE
      local_vapour = saturation_pressure
    END IF

    !----------------------------------------
    !
    ! 2. Calculate the density of water vapour in air
    !    ------------------------------------

    fu_watervapour_densityprofile = 0.2167*(local_vapour/temperature)


  END FUNCTION fu_watervapour_densityprofile


 ! ****************************************************************
  
  REAL FUNCTION fu_dynamic_viscosity(temperature) result(dyn_viscosity)
    ! Description:
    !  Air dynamic viscisity
    ! Fit to the data 
!# https://www.engineersedge.com/physics/viscosity_of_air_dynamic_and_kinematic_14483.htm
!#T(C) mu(Pa s)
!-50  1.474e-5
!-40  1.527e-5
!-30  1.579e-5
!-20  1.630e-5
!-10  1.680e-5
!0    1.729e-5
!5    1.754e-5
!10   1.778e-5
!15   1.802e-5
!20   1.825e-5
!25   1.849e-5
!30   1.872e-5
!35   1.895e-5
!40   1.918e-5
!50   1.963e-5

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: temperature

    ! Local declarations
    REAL :: local_temp
    real, parameter :: A = 4.89032e-08, B = 1.72558e-05  
 
    dyn_viscosity = A*(temperature - zero_celcius) + B
    
  END FUNCTION fu_dynamic_viscosity


  ! ****************************************************************
  
  REAL FUNCTION fu_kinematic_viscosity(temperature, density) result(kin_viscosity)
    
    ! Description:
    !  Computes the "real" value for kinematic viscosity
    !  using the ambient air temperature &  density
    !
    ! Method:
    !   Approximation of the general formula
    ! 
    !  Should the air density be missing a cosntant r= 1.225 is used
    !  If temperature is missing a constant value is returned
    !
    ! Units: air density = kg/m**3
    !        temperature = Kelvins
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Ilkka Valkama, FMI
    
    IMPLICIT NONE
    
    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: temperature, density

       kin_viscosity =  fu_dynamic_viscosity(temperature)/density  ! typical value: 1.461E-5
    
  END FUNCTION fu_kinematic_viscosity


  ! ****************************************************************
  
  REAL FUNCTION  fu_watervap_diffusivity(temperature, pressure) result(wvap_diffusiv)
    
    ! Description:
    !  Computes the "real" value for kinematic viscosity
    !  using the ambient air temperature &  density
    !
    ! Method:
    !   Approximation of the general formula
    ! 
    !  Should the air density be missing a cosntant r= 1.225 is used
    !  If temperature is missing a constant value is returned
    !
    ! Units: air density = kg/m**3
    !        temperature = Kelvins
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Ilkka Valkama, FMI
    
    IMPLICIT NONE
    
    ! Imported parameters with intent(in):
    REAL, INTENT(in), OPTIONAL :: temperature, pressure

    ! Local declarations
    REAL :: local_temp, local_pres
 
    IF (.NOT.PRESENT(pressure)) THEN
      ! 1000 hPa assumed, 20C
      wvap_diffusiv = 2.42E-5
      RETURN
    END IF
    
    IF (.NOT.PRESENT(temperature)) THEN
      ! 1000 hpa, 20C assumed
      wvap_diffusiv = 2.42e-5
      RETURN
    END IF
    
    !  Calculate kinematic viscosity
    !  ------------------------------------
    !  
    local_temp = (temperature/273.15)**1.94
    wvap_diffusiv = 2.11E-5*local_temp*(101325.0/pressure)
    
  END FUNCTION fu_watervap_diffusivity



  ! ****************************************************************
  
  REAL FUNCTION  fu_air_diffusivity(temperature) result(air_diffusiv)
    
    ! Description:
    !  Computes the fit to air diffusivity using the ambient air temperature
    !
    ! Method:
    !   Approximation of the general fitting formula
    ! 
    !  If temperature is missing a constant value is returned
    !
    ! Units: temperature = Kelvins, diffusivity m**2/sec
    !
    ! Author: Mikhail Sofiev, FMI
    
    IMPLICIT NONE
    
    ! Imported parameters with intent(in):
    REAL, INTENT(in), OPTIONAL :: temperature

    IF (.NOT.PRESENT(temperature)) THEN
      ! 20C assumed
      air_diffusiv = 2.3e-5
      RETURN
    END IF
    
    !  Calculate the diffusivity
    !  
    air_diffusiv = (9.1018e-11 * temperature + 8.8197e-8) * temperature - 1.0654e-5
    
  END FUNCTION fu_air_diffusivity


  !************************************************************************

  real function fu_R_a(z, z0, L_inv, u_star)
    !
    ! Computes the Ra resistance for the given parameters
    !
    implicit none

    real, intent(in) :: z, L_inv, u_star
    real, intent(inout) :: z0

    if(z0 < 1e-10)then
!      call msg('Strange z0:', z0)
      z0=1.e-6
    endif
    if(z<z0*1.01)then ! within or very close to the roughness layer
      fu_R_a = 0.01  ! for sake of stability, do not set 0.
    else
      fu_R_a = (log(z/z0) - psih_prof(z*L_inv) + psih_prof(z0*L_inv)) / (karmann_c * u_star)
    endif

  end function fu_R_a



  ! ***************************************************************

  REAL FUNCTION vapor_pressure_dewp(dew_point) result(vapor_press)

    ! Description:
    ! 
    ! Computes the partial water vapor pressure  [Pa]
    ! for the given dew point temperature  [K] 
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Ilkka Valkama, FMI
    
    IMPLICIT NONE
    
    ! Imported parameters with intent(in/out):
    REAL, INTENT(in) :: dew_point
    
 
    !  local declarations
    REAL :: corr_a, corr_b, corr_c, corr_d
    
    !  some preliminaries

    IF(dew_point.LE.0.9) THEN
      vapor_press=0.0
      PRINT*,' SORRY: dew-point must be in Kelvins !'
      RETURN
    END IF

    ! let's start for real

    corr_a = 373.16/dew_point

    corr_b = -7.90298*(corr_a - 1.0)
    corr_b = corr_b + (5.02808*0.43429*ALOG(corr_a))

    corr_c = (1.0 - (1.0/dew_point))*11.344
    corr_c = -1.0 + (10.0**corr_c)
    corr_c= -1.3816*corr_c/(10.0**7)

    corr_d = (1.0 -  corr_a)*3.49149
    corr_d = -1.0 + (10.0**corr_d)
    corr_d = 8.1328*corr_d/(10.**3)
    
    corr_a = corr_b + corr_c + corr_d
    vapor_press = 101324.6*(10.**corr_a)
 
  END FUNCTION VAPOR_PRESSURE_DEWP


  !******************************************************************

  SUBROUTINE us_standard_atmosphere(alt_, sigma, delta, theta)
    !   -------------------------------------------------------------------------
    ! PURPOSE - Compute the properties of the 1976 standard atmosphere to 86 km.
    ! Author: - Ralph Carmichael, Public Domain Aeronautical Software
    ! NOTE - If alt > 86, the values returned will not be correct, but they will
    !   not be too far removed from the correct values for density.
    !   The reference document does not use the terms pressure and temperature
    !   above 86 km.
    IMPLICIT NONE
    !============================================================================
    !     A R G U M E N T S                                                     |
    !============================================================================
    REAL,INTENT(IN)::  alt_        ! geometric altitude, m.
    REAL,INTENT(OUT):: sigma      ! density/sea-level standard density
    REAL,INTENT(OUT):: delta      ! pressure/sea-level standard pressure
    REAL,INTENT(OUT):: theta      ! temperature/sea-level standard temperature
    !============================================================================
    !     L O C A L   C O N S T A N T S                                         |
    !============================================================================
    REAL,PARAMETER:: REARTH = 6369.0                 ! radius of the Earth (km)
    REAL,PARAMETER:: GMR = 34.163195                     ! hydrostatic constant
    INTEGER,PARAMETER:: NTAB=8       ! number of entries in the defining tables
    !============================================================================
    !     L O C A L   V A R I A B L E S                                         |
    !============================================================================
    INTEGER:: i,j,k                                                  ! counters
    REAL:: h                                       ! geopotential altitude (km)
    REAL:: tgrad, tbase      ! temperature gradient and base temp of this layer
    REAL:: tlocal                                           ! local temperature
    REAL:: deltah                             ! height above base of this layer
    !============================================================================
    !     L O C A L   A R R A Y S   ( 1 9 7 6   S T D.  A T M O S P H E R E )   |
    !============================================================================
    REAL,DIMENSION(NTAB),PARAMETER:: htab= &
                            (/0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852/)
    REAL,DIMENSION(NTAB),PARAMETER:: ttab= &
            (/288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946/)
    REAL,DIMENSION(NTAB),PARAMETER:: ptab= &
                 (/1.0, 2.233611E-1, 5.403295E-2, 8.5666784E-3, 1.0945601E-3, &
                                       6.6063531E-4, 3.9046834E-5, 3.68501E-6/)
    REAL,DIMENSION(NTAB),PARAMETER:: gtab= &
                                  (/-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0/)
    !----------------------------------------------------------------------------
    
    h=alt_*0.001*REARTH/(alt_*0.001+REARTH)      ! convert geometric to geopotential altitude

    i=1
    j=NTAB                                       ! setting up for binary search
    DO
      k=(i+j)/2                                              ! integer division
      IF (h < htab(k)) THEN
        j=k
      ELSE
        i=k
      END IF
      IF (j <= i+1) EXIT
    END DO

    tgrad=gtab(i)                                     ! i will be in 1...NTAB-1
    tbase=ttab(i)
    deltah=h-htab(i)
    tlocal=tbase+tgrad*deltah
    theta=tlocal/ttab(1)                                    ! temperature ratio

     IF (tgrad == 0.0) THEN                                     ! pressure ratio
      delta=ptab(i)*EXP(-GMR*deltah/tbase)
    ELSE
      delta=ptab(i)*(tbase/tlocal)**(GMR/tgrad)
    END IF

    sigma=delta/theta                                           ! density ratio
    RETURN
  END Subroutine us_standard_atmosphere   ! ----------------------------------------------- 


  ! ***************************************************************

  REAL FUNCTION fu_pr_diff_height(h1_m, h2_m) result (delta_pre)
    ! 
    ! Calculates the pressure difference, when moving from h1 to h2
    ! heights. The dumbest possible way is to use Standard Atmosphere routines.
    ! Since th function can be called inside some cycle, we shall use tabulated
    ! values for pressure and temperature
    !
    ! All units: SI
    !
    IMPLICIT NONE
    
    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: h1_m, h2_m
    
    ! Local variables
    real :: p1, p2, sigma, delta, theta

    call us_standard_atmosphere(h1_m, sigma, delta, theta)
    call us_standard_atmosphere(h2_m, sigma, delta, theta)
    delta_pre = (p2 - p1) * std_pressure_sl

  END FUNCTION fu_pr_diff_height


  ! ***************************************************************

  REAL FUNCTION fu_height_diff_between_pre(p1_Pa, p2_Pa) result(delta_height)
    !
    ! Calcuates height difference between two pressure-values using US stdandard atmosphere
    ! Since can be called many times, creates and uses tabulated array
    !
    ! All units: SI
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: p1_Pa, p2_Pa

    ! Loca variables
    real :: fInd1, fInd2
    integer :: iVal
    
    if(.not. ifAtmosphereTabulated)then
      do iVal = 1, nStandardValues
        call us_standard_atmosphere(real((iVal-1) * delta_standard_table), arStandardDensity(iVal), &
                                  & arStandardPressure(iVal), arStandardTemperature(iVal))
      enddo
      ifAtmosphereTabulated = .true.
    endif
    
    fInd1 = fu_value_index_in_array(p1_Pa / std_pressure_sl, arStandardPressure, nStandardValues)
    fInd2 = fu_value_index_in_array(p2_Pa / std_pressure_sl, arStandardPressure, nStandardValues)
    
    delta_height = (fInd2 - fInd1) * delta_standard_table

  END FUNCTION fu_height_diff_between_pre

  subroutine std_atm_rho_t(pressure, rho_over_rhoref, T_over_Tref) 
    !
    ! Calcuates standard density and temperature for given pressure
    ! using US stdandard atmosphere
    ! All units: SI
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: pressure
    REAL, INTENT(out) :: rho_over_rhoref, T_over_Tref


    ! Loca variables
    real :: fInd1, weight_up
    integer :: iVal
    
    if(.not. ifAtmosphereTabulated)then
      do iVal = 1, nStandardValues
        call us_standard_atmosphere(real((iVal-1) * delta_standard_table), arStandardDensity(iVal), &
                                  & arStandardPressure(iVal), arStandardTemperature(iVal))
      enddo
      ifAtmosphereTabulated = .true.
    endif
    
    fInd1 = fu_value_index_in_array(pressure / std_pressure_sl, arStandardPressure, nStandardValues)
    iVal = floor(fInd1)
    weight_up = fInd1 - floor(fInd1)
    if (iVal < 1 .or. iVal > nStandardValues-1) then
      call set_error("Can'T_over_Tref get us_standard_atmosphere for pessure: "//trim(fu_str(pressure)), &
                & "std_atm_rho_over_rhoref_t")
    endif

    rho_over_rhoref =     arStandardDensity(iVal)*(1.-weight_up) +     arStandardDensity(iVal+1)*weight_up
    T_over_Tref   = arStandardTemperature(iVal)*(1.-weight_up) + arStandardTemperature(iVal+1)*weight_up

  END subroutine std_atm_rho_t


  !****************************************************

  !  Private tools : Ilkka Valkama , Nov 2000 
  !     PRIVATE psim_prof
  !     PRIVATE psih_prof
  ! ***************************************************************


  !************************************************************************
  
  REAL FUNCTION psim_prof(z_per_mo_length) result(mom_gradient)
    ! 
    ! Evaluates the surface layer flux gradient function 
    ! for momentum  : the psim_function
    ! 
    ! Method : The classical approximation
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Ilkka Valkama, FMI
    
    IMPLICIT NONE
    
    ! Imported parameters with intent(in/out):
    REAL, INTENT(in) :: z_per_mo_length
    
    !  local declarations
    REAL :: local_grad,local_gr1,local_gr2
    
    IF(z_per_mo_length <= 0.0) THEN
      ! unstable CASE
      
      local_grad= (1.0 - 15.0*z_per_mo_length)**0.25
      local_gr1= (0.5 + 0.5*local_grad)*(0.5 + 0.5*local_grad)
      local_gr2= (1.0 + local_grad*local_grad)*0.5
      mom_gradient = LOG(local_gr1*local_gr2) - 2.0*ATAN(local_grad) + pi*0.5
    ELSE
      ! stable CASE
      
      mom_gradient = -4.7*z_per_mo_length
      
    ENDIF
  END FUNCTION psim_prof
  

  !************************************************************************

  REAL FUNCTION psih_prof(z_per_mo_length) result(temp_gradient)
    
    ! Description:
    ! 
    ! Evaluates the surface layer flux gradient function 
    ! for heat  : the psih_function
    ! 
    ! Method : The classical approximations of 
    !          Paulson (1970) and Beljaars & Holtslag (1991)
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Ilkka Valkama, FMI
    
    IMPLICIT NONE
    
    ! Imported parameters with intent(in/out):
    REAL, INTENT(in) :: z_per_mo_length
    
    !  constants for the vertical temperature gradients (psih-functions)
    !  Beljaars & Holtslag 1999
    !
    REAL, PARAMETER :: a_BH1991= 1.0
    REAL, PARAMETER :: b_BH1991= 0.667
    REAL, PARAMETER :: c_BH1991= 5.0
    REAL, PARAMETER :: d_BH1991= -0.35

    !  local variables
    REAL :: local_grad,local_gr1,local_gr2
    
    IF(z_per_mo_length <= 0.0) THEN
      ! unstable CASE
      
      local_grad= (1.0 - 16.0*z_per_mo_length)**0.25
      temp_gradient = 2.0*LOG((1.0 + local_grad*local_grad)/2.0)
      
    ELSE
      ! stable CASE
      local_gr1 = - (1.+0.667*a_BH1991*z_per_mo_length)**(1.5)
      local_gr2 = b_BH1991*(z_per_mo_length + c_BH1991/d_BH1991)*EXP(d_BH1991*z_per_mo_length)
      temp_gradient = local_gr1 - local_gr2 + b_BH1991*c_BH1991/d_BH1991 + 1.
      
    ENDIF
    
  END FUNCTION psih_prof



  ! ****************************************************************
  
  REAL FUNCTION fu_dew_point_temperature(humidity_mixing_ratio, pressure)
    
    ! Description:
    ! Returns the dew point temperature for a position 
    !
    ! Dew point temperature defined as the temperature to which moist
    ! air must be cooled, with pressure and humidity mixing ratio
    ! held constant, for it to reach saturation with respect to water.
    !
    ! Method: An anlytical approximation from Rogers & Yau (1991), A
    ! Short Course in Cloud Physics.
    !
    ! Humidity mixing ratio = water/dry_air [kg/kg]
    ! Note: specific humidity = water/(dry_air+water) [kg/kg]
    ! The difference is small for all reasonable casses
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Pilvi Siljamo, FMI
    ! Another formula is taken from Web by Marje Prank

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in):: humidity_mixing_ratio, pressure



!    fu_dew_point_temperature = 5.42E3/LOG((2.53E11*gas_constant_ratio) / &
!                                        & (humidity_mixing_ratio*pressure)) 

     !
     ! New formulation
     !
     fu_dew_point_temperature = 1./ &
              & (1./273. - log(humidity_mixing_ratio * pressure / (gas_constant_ratio * 611.)) * &
                         & gas_constant_watervapour / vaporization_latentheat)

  END FUNCTION fu_dew_point_temperature


  ! ****************************************************************

  REAL FUNCTION fu_virtual_temp_from_dewpoint(temperature,&
      & dew_point_temperature,pressure) result(virtual_temp)                                          

    ! Description:
    ! Returns the virtual temperature from  air temperature, 
    ! dew point temperature and pressure 
    !
    ! Virtual temperature is the temperature of dry air having the
    ! same  density as a sample of moist air at the same pressure.
    !
    ! Method:  Using the definition of T_virt
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Ilkka Valkama, FMI

    IMPLICIT NONE
    
    ! Imported parameters with intent(in):
    REAL, INTENT(in) ::  temperature, dew_point_temperature, pressure
    
    ! Local declarations
    REAL :: vapor_pres
    
    !----------------------------------------
    !
    ! 1. Calculate the virtual temperature
    !    ------------------------------------
    
    IF(pressure.GT.0.0) THEN
      vapor_pres = fu_vapor_pressure_dew_point(dew_point_temperature)
      virtual_temp = temperature*(1.0 + 0.378*vapor_pres/pressure)
    ELSE
      virtual_temp = 0.0
    END IF
    
  END FUNCTION fu_virtual_temp_from_dewpoint

  ! ****************************************************************

  REAL FUNCTION fu_virtual_temperature(&
      & temperature,&
      & humidity_mixing_ratio)                                           

    ! Description:
    ! Returns the virtual temperature from air temperature and
    ! humidity mixing ratio
    !
    ! Virtual temperature is the temperature of dry air having the
    ! same  density as a sample of moist air at the same pressure.
    !
    ! Method:  Using the definition of T_virt
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Ilkka Valkama, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: temperature, humidity_mixing_ratio

    ! Local declarations
    REAL:: local_humidity


    local_humidity = (1.0 + 1.61*humidity_mixing_ratio)/(1.0 + humidity_mixing_ratio)
    fu_virtual_temperature = temperature *local_humidity

  END FUNCTION fu_virtual_temperature


  ! ****************************************************************

  REAL FUNCTION fu_virtual_temp_from_spec_hum (temperature, specific_humidity)                                           

    ! Description:
    ! Returns the virtual temperature from air temperature and
    ! specific humidity
    !
    ! Virtual temperature is the temperature of dry air having the
    ! same  density as a sample of moist air at the same pressure.
    !
    ! Method:   Using the definition of T_virt
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Ilkka Valkama, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: temperature, specific_humidity

    ! Local declarations
    REAL:: local_humidity

    local_humidity = (0.622 + specific_humidity)/(0.622*(1.0 + specific_humidity))
    fu_virtual_temp_from_spec_hum = temperature *local_humidity

  END FUNCTION fu_virtual_temp_from_spec_hum

  ! ****************************************************************

  REAL FUNCTION fu_potential_temperature(temperature, pressure)                                           

    ! Description:
    ! Returns the potential temperature from air temperature and
    ! pressure
    !
    ! Method: Using the definition of potential temperature 
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Pilvi Siljamo, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: temperature, pressure

    fu_potential_temperature = temperature * (std_pressure_sl/pressure)**R_per_c_dryair

  END FUNCTION fu_potential_temperature



  ! ****************************************************************

  REAL FUNCTION fu_temperature_from_pot_temp(potential_temp, pressure)

    ! Description:
    ! Returns the air temperature from potential temperature and
    ! pressure
    !
    ! Method:  Using the definition of potential temperature 
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Pilvi Siljamo, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: potential_temp, pressure

    fu_temperature_from_pot_temp = potential_temp * (pressure / std_pressure_sl)**R_per_c_dryair

  END FUNCTION fu_temperature_from_pot_temp



  ! ****************************************************************

  REAL FUNCTION fu_equivalent_temperature(temperature, humidity_mixing_ratio)

    ! Description:
    ! Returns the equivalent temperature from air temperature and
    ! humidity mixing ratio
    !
    ! Equivalent temperature obtained by following up the
    ! pseudoadiabat from isentropic condensation point to very low
    ! pressure, thus condensing out all the water vapor, and then
    ! returning to the original pressure along a dry adiabat.
    !
    ! Method: From Rogers & Yau (1991), A Short  Course in Cloud Physics. 
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Pilvi Siljamo, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: temperature, humidity_mixing_ratio

    fu_equivalent_temperature = temperature + &
                              & (condensation_latentheat*humidity_mixing_ratio / &
                              &  specific_heat_dryair) 

  END FUNCTION fu_equivalent_temperature


  ! ****************************************************************

  REAL FUNCTION fu_equivalent_pot_temperature(&
      & temperature,&
      & pressure,&
      & humidity_mixing_ratio)

    ! Description:
    ! Returns the equivalent potential temperature from 
    ! potential temperature and humidity saturation mixing ratio
    !
    ! Equivalent potential temperature defined as the temperature a
    ! parcel of air would have if taken from its equivalent
    ! temperature to a pressure of 1000 hPa in a dry adiabatic process.
    !
    ! Method: Using the definition of potential temperature 
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Pilvi Siljamo, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: temperature, pressure, humidity_mixing_ratio

    fu_equivalent_pot_temperature =(temperature +&
                                  & (condensation_latentheat*humidity_mixing_ratio&
                                  &  /specific_heat_dryair))*&
                                  & (1E5/pressure)**R_per_c_dryair 


  END FUNCTION fu_equivalent_pot_temperature


  ! ****************************************************************

  REAL FUNCTION fu_wet_bulb_temperature(temperature,&
                                      & humidity_mixing_ratio,&
                                      & pressure)

    ! Description:
    ! Returns the wet bulb temperature from air temperature,
    ! humidity mixing ratio and humidity saturation mixing ratio
    !
    ! Wet bulb temrerature is the temperature to which air may be
    ! cooled by evaporating water into it at constant pressure, until
    ! saturation is reached.
    !
    ! Method:  From Rogers & Yau (1991), 
    ! A Short  Course in Cloud Physics  
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Pilvi Siljamo, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: temperature, humidity_mixing_ratio,&
	&pressure

    ! Local declarations
    REAL :: mixing_ratio_at_wet_bulb_temp, wet_bulb_temperature
    REAL::  wet_bulb_temperature2, dT_w
    INTEGER :: loop
    !----------------------------------------
    !
    ! 1. Calculate the wet bulb temperature
    !    ------------------------------------   


    wet_bulb_temperature = temperature - 3
    dT_w = 100
    wet_bulb_temperature2 = temperature
    loop = 0
    DO WHILE (dT_w > 0.01 .and. loop < 100)   

      mixing_ratio_at_wet_bulb_temp = 157.366E9/pressure *&
                                    & EXP((-5420)/wet_bulb_temperature)

      wet_bulb_temperature = temperature - &
       & (condensation_latentheat&
       & /specific_heat_dryair)*(mixing_ratio_at_wet_bulb_temp - &
       & humidity_mixing_ratio)

      dT_w = ABS( wet_bulb_temperature -  wet_bulb_temperature2)

      wet_bulb_temperature = 0.5*( wet_bulb_temperature + wet_bulb_temperature2)

      wet_bulb_temperature2  = wet_bulb_temperature

      loop = loop + 1 
    END DO

    fu_wet_bulb_temperature =  wet_bulb_temperature

  END FUNCTION fu_wet_bulb_temperature


  ! ****************************************************************

  REAL FUNCTION fu_wet_bulb_pot_temperature(&
      & temperature,&
      & pressure,&
      & humidity_mixing_ratio)


    ! Description:
    ! Returns the wet bulb potential temperature from air temperature,
    ! humidity mixing ratio and pressure
    !
    ! Wet bulb potential temperature defined by the intersection of
    ! the pseudoadiabat through the isentropic condensation point with the
    ! isobar  p=1000 hPa.
    !
    ! Method: From Rogers & Yau (1991), 
    ! A Short  Course in Cloud Physics and using the definition of
    ! potential temperature by pseudoadiabat 
    ! 
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Pilvi Siljamo, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: temperature, pressure,&
	& humidity_mixing_ratio

    ! Local declarations
    REAL :: mixing_ratio_at_wet_bulb_temp, wet_bulb_temperature,&
	& wet_bulb_temperature2, dT_w, pseudoadiabatic_lapse_rate
    INTEGER :: loop

    !-------------------------------------------------
    !
    ! 1. Calculate the wet bulb potential temperature
    !    ---------------------------------------------   

    wet_bulb_temperature = temperature - 3
    dT_w = 100
    wet_bulb_temperature2 = temperature
    loop = 0

    DO WHILE (loop < 100 .and. dT_w > 0.01)

      mixing_ratio_at_wet_bulb_temp = 157.366E9/pressure&
	  & *EXP((-5420)/wet_bulb_temperature)

      wet_bulb_temperature = temperature - &
	  & (condensation_latentheat&
 	  & /specific_heat_dryair)*(mixing_ratio_at_wet_bulb_temp - &
	  & humidity_mixing_ratio)

      dT_w = ABS( wet_bulb_temperature -  wet_bulb_temperature2)
      wet_bulb_temperature = 0.5*(wet_bulb_temperature + wet_bulb_temperature2)
      wet_bulb_temperature2 = wet_bulb_temperature
      loop = loop + 1

    END DO

    pseudoadiabatic_lapse_rate = R_per_c_dryair*(1+&
	& (condensation_latentheat* (157.366E9/pressure&
	& *EXP((-5420)/wet_bulb_temperature)))/&
	& (gas_constant_dryair*temperature))/(1+&
	& (((condensation_latentheat**2)*gas_constant_ratio*&
	&  (157.366E9/pressure&
	& *EXP((-5420)/wet_bulb_temperature)))/(gas_constant_dryair&
	& *specific_heat_dryair*(temperature**2)))) 

    fu_wet_bulb_pot_temperature = wet_bulb_temperature*&
	& (1.0E5/pressure)**&
	& pseudoadiabatic_lapse_rate

  END FUNCTION fu_wet_bulb_pot_temperature


  ! ****************************************************************^L 

  REAL FUNCTION fu_moist_potential_temperature(&
      & temperature,&
      & pressure,&
      & humidity_mixing_ratio)

    ! Description:
    ! Returns the moist  potential temperature from air
    ! temperature, humidity mixing ratio and pressure
    !    
    ! Method: 
    ! 
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Pilvi Siljamo, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: temperature, pressure, humidity_mixing_ratio

    fu_moist_potential_temperature = temperature*&
	& (1.0E5/pressure)**&
	&(R_per_c_dryair*&
	&((1 +0.608*humidity_mixing_ratio)/&
	&(1+0.887 *humidity_mixing_ratio)))

  END FUNCTION fu_moist_potential_temperature


  ! **************************************************************** 

  REAL FUNCTION fu_lcl_temperature_from_rel_hum(&
      & low_level_temperature,&
      & relative_humidity)

    ! Description:
    ! Returns the isentropic condensation temperature (temperature at
    ! lifting  condensation level) from air temperature and relative
    ! humidity
    !
    ! Isentropic condensation temperature defined as the temperature
    ! at which saturation is reached when moist air is cooled
    ! adiabatically with mixing ratio held constant.
    !
    ! Method:  
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Pilvi Siljamo, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: low_level_temperature ! for example lowest
    ! model level or 2m level
    REAL, INTENT(in) ::  relative_humidity

    !----------------------------------------
    !
    ! 1. Calculate the potential temperature
    !    ------------------------------------

    fu_lcl_temperature_from_rel_hum = 1/(1/(low_level_temperature-55)&
	& -(LOG(relative_humidity)/2840))+55


  END FUNCTION fu_lcl_temperature_from_rel_hum



  ! ****************************************************************

  REAL FUNCTION fu_lcl_temp_from_dew_point(&
      & low_level_temperature,&
      & dew_point_temperature)

    ! Description:
    ! Returns the isentropic condensation temperature (temperature at
    ! lifting  condensation level) from air temperature and dew point
    ! temperature. 
    !
    ! Isentropic condensation temperature defined as the temperature
    ! at which saturation is reached when moist air is cooled
    ! adiabatically with mixing ratio held constant.
    !
    ! Method:   
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Pilvi Siljamo, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in) ::low_level_temperature ! for example lowest
    ! model level or 2m level
    REAL, INTENT(in) ::dew_point_temperature
    !----------------------------------------
    !
    ! 1. Calculate the potential temperature
    !    ------------------------------------

    fu_lcl_temp_from_dew_point = ((1/(dew_point_temperature - 56))+&
                   &((LOG(low_level_temperature/dew_point_temperature))/800))**(-1)+56


  END FUNCTION fu_lcl_temp_from_dew_point


  ! ****************************************************************

  REAL FUNCTION fu_lcl_pressure(&
      & low_level_temperature,&
      & lcl_temperature,&
      & low_level_pressure)

    ! Description:
    ! Returns the isentropic condensation pressure from air
    ! temperature and isentropic condensation temperature and pressure.
    !
    ! Method:  
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Pilvi Siljamo, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: low_level_temperature ! for example lowest
    ! model level or 2m level
    REAL, INTENT(in) :: lcl_temperature
    REAL, INTENT(in) :: low_level_pressure
    !----------------------------------------
    !
    ! 1. Calculate the potential temperature
    !    ------------------------------------

    fu_lcl_pressure = low_level_pressure*(lcl_temperature&
	& /low_level_temperature)**&
	&(1/R_per_c_dryair) 

  END FUNCTION fu_lcl_pressure


  ! ****************************************************************

  REAL FUNCTION fu_pseudoadiab_equiv_pot_temp(wet_bulb_pot_temperature,&
                                            & lcl_temperature,&
                                            & humidity_mixing_ratio)

    ! Description:
    ! Returns the saturated, pseudoadiabatic equivalent potential
    ! temperature from wet bulb potential temperature,  isentropic
    ! condensation temperature and humidity mixing ratio.
    !
    ! Method:  
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Pilvi Siljamo, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: wet_bulb_pot_temperature, lcl_temperature,&
        & humidity_mixing_ratio

    !----------------------------------------
    !
    ! 1. Calculate the potential temperature
    !    ------------------------------------

    fu_pseudoadiab_equiv_pot_temp =  wet_bulb_pot_temperature * &
                & EXP(((3.376/lcl_temperature)-0.00254)*humidity_mixing_ratio&
                    & *(1+0.81*humidity_mixing_ratio))  

  END FUNCTION fu_pseudoadiab_equiv_pot_temp



  ! ****************************************************************

  REAL FUNCTION fu_spechum_to_relhum(q, t, p) result(rh)

    ! Description:
    ! Calculates relative humidity (0...1.) from specific humidity.
    !
    ! Method:
    ! Copied from HIRLAM's QTORH1
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
    REAL, INTENT(in) :: q, t, p
    !
    ! Local declarations:
    REAL :: e_sat

    e_sat = fu_watervapour_satur_p(t)

    rh = q*(p-(1.-gas_constant_ratio)*e_sat)/(gas_constant_ratio*e_sat)

    rh = MAX(rh, 0.)
    rh = MIN(rh, 1.)

  END FUNCTION fu_spechum_to_relhum


  ! ****************************************************************

  SUBROUTINE spechum_to_relhum_field(q, t, p, rh)

    ! Description:
    ! Calculates relative humidity (0...1.) field from specific humidity,
    ! t and p -fields. Sizes of q, t and p must be the same.
    !
    ! Method for RelHum:
    ! Copied from HIRLAM's QTORH1
    !
    ! Method for eSat:
    ! HIRLAM,  Gerard Cats  KNMI    24 September 1992
    ! Possible reference: Ivarson, K.I. Tests with separated tables for
    ! water vapor saturation pressure over ice and over water HIRLAM 
    ! technical report N 45, April 2000, Norrkoping Sweden, 22 pp.
    !
    ! All units: SI
    !
    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, DIMENSION(:), INTENT(in) :: q, t, p

    ! Imported parameters with intent(out):
    REAL, DIMENSION(:), INTENT(out) :: rh

    ! Local declarations:
!    REAL, DIMENSION(:), POINTER :: e_sat
    real e_sat
    INTEGER :: fs, i
    real, parameter :: one_m_gas_const_ratio = 1.0 - gas_constant_ratio

    REAL :: zc2, zc4, zes
    REAL, PARAMETER :: c1es = 610.78, &
                     & c2es =  17.269, &
                     & c2is =  21.875, &
                     & c3es = 273.16, &
                     & c4es =  35.86, &
                     & c4is =   7.66

    fs = SIZE(q)
    
    !$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE(zes, zc2, zc4, i, e_sat)
    do i = 1, fs
      zes = (t(i) - c3es + 15.)/15.
      zes = min(1.,max(zes,0.))   ! 0<= zes <= 1
      zc2 = c2is + (zes * (c2es-c2is))
      zc4 = c4is + (zes * (c4es-c4is))
      e_sat = c1es * EXP(zc2 *(t(i) - c3es)/(t(i)-zc4))
      rh(i) = q(i) * (p(i) - one_m_gas_const_ratio * e_sat) / (gas_constant_ratio * e_sat)
    end do
    !rh(1:fs) = q(1:fs)* (p(1:fs) - (1.-gas_constant_ratio) *e_sat(1:fs))/ &
    !         & (gas_constant_ratio*e_sat(1:fs))
    !$OMP END PARALLEL DO

  END SUBROUTINE spechum_to_relhum_field


  !***************************************************************************
  
  real function fu_wind_gust_default_10m(uStar)
    !
    ! Computes the max absolute wind gust for the given probability and 
    ! time period. Methodology: Schreur & Geertsema, HIRLAM Newsletter no 54, December 2008
    ! gustRelative = sqrt(2 * ln(refPeriod/(gustTime * sqrt(2*Pi) * ln(1./probability))))
    ! gustAbsolute = windSpeed + gustRelative * 2.185 * uStar
    ! WMO default values: probability = 0.5, gustTime = 3 sec, refPeriod = 1 hr
    ! We provide gust scaled with absolute 10m wind speed
    !
    implicit none
    
    ! Imported parameters
    real, intent(in) :: uStar

    !
    ! Default values lead to gustRelative = 3.616
    !
!    fu_wind_gust_default_10m = 1. + 7.9 * uStar / windSpeed 
    fu_wind_gust_default_10m = 7.9 * uStar   !/ windSpeed 
    
  end function fu_wind_gust_default_10m


  !***************************************************************************
  
  real function fu_wind_gust_general_10m(uStar, probability, gustTime, referencePeriod)
    !
    ! Computes the max absolute wind gust for the given probability and 
    ! time period. Methodology: Schreur & Geertsema, HIRLAM Newsletter no 54, December 2008
    ! gustRelative = sqrt(2 * ln(refPeriod/(gustTime * sqrt(2*Pi) * ln(1./probability))))
    ! gestAbsolute = windSpeed + gustRelative * 2.185 * uStar
    ! We provide gust scaled with absolute wind speed
    !
    implicit none
    
    ! Imported parameters
    real, intent(in) :: uStar, probability, gustTime, referencePeriod

    !
    ! The return value is actual extra 10m-wind, in m/s
    !
    fu_wind_gust_general_10m = 2.185 * uStar * &
                                         & sqrt(2.0 * log(-referencePeriod / &
                                                        & (gustTime * 2.50663 * log(probability))))
  end function fu_wind_gust_general_10m


  !**********************************************************************
  
  subroutine make_gust_lookup_table(fWindPower)
    !
    ! Creates new gust lookup table with appropriate space reservation
    !
    implicit none
    
    ! Imported parameters
    real, intent(in) :: fWindPower
    
    ! Local variables
    type(silja_rp_2d), dimension(:), pointer :: arTmp
    integer :: iWind, iZ, iAr
    !
    ! Get the space
    !
    if(wgl%nWindPwr > 0)then            
      !
      ! store old stuff
      call get_work_arrays_set(wgl%nWindPwr, arTmp)
      if(error)return
      do iAr = 1, wgl%nWindPwr
        do iZ = 1, wgl%nZ
          do iWind = 1, wgl%nWinds
            arTmp(iAr)%pp(iWind,iZ) = wgl%pGustVal(iWind,iZ,iAr)
          end do  ! wind-dim
        end do  ! z-dim
      end do  ! existing lookup tables
      deallocate(wgl%pGustVal)  ! kill the old table
    else
      !
      ! Make the scanning axes
      wgl%nZ = nZ_glob
      wgl%nWinds = nWinds_glob
      wgl%nWindPwr = 0
      allocate(wgl%fWind(wgl%nWinds), wgl%fZr(wgl%nZ), wgl%fWindPwr(1), stat=iWind)
      if(fu_fails(iWind==0,'Failed gust lookup axes allocation','make_gust_lookup_table'))return
      do iWind = 1, wgl%nWinds
        if(iWind == 1)then
          wgl%fWind(iWind) = 0.01   !m/s
        else
          wgl%fWind(iWind) = real(iWind-1) / 2.0  ! m/s
        endif
      end do
      do iZ = 1, wgl%nZ
        wgl%fZr(iZ) = 1.5**(iZ)  ! Z/z0, relative
      end do
    endif  ! if some lookup tables already exist
    !
    ! (Re)allocate and put old stuff back
    !
    allocate(wgl%pGustVal(wgl%nWinds, wgl%nZ, wgl%nWindPwr+1), stat=iZ)
    if(fu_fails(iZ==0,'Failed reallocation of windgust lookup table','make_gust_lookup_table'))return

    if(wgl%nWindPwr > 0)then
      do iAr = 1, wgl%nWindPwr    ! restore old stuff
        do iZ = 1, wgl%nZ
          do iWind = 1, wgl%nWinds
            wgl%pGustVal(iWind,iZ,iAr) = arTmp(iAr)%pp(iWind,iZ)
          end do  ! wind-dim
        end do  ! z-dim
      end do  ! existing lookup tables
      call free_work_array(arTmp)
    endif
    !
    ! Add the new lookup table
    !
    wgl%nWindPwr = wgl%nWindPwr + 1
    wgl%fWindPwr(wgl%nWindPwr) = fWindPower
    do iZ = 1, wgl%nZ
      do iWind = 1, wgl%nWinds
        wgl%pGustVal(iWind,iZ,wgl%nWindPwr) = fu_gust_integral(wgl%fWind(iWind), wgl%fZr(iZ), &
                                                             & wgl%fWindPwr(wgl%nWindPwr))
        if(error)return
      end do  ! wind-dim
    end do  ! z-dim

!do iZ = 1, wgl%nZ
!write(run_log_funit,fmt='(1000(F5.2,1x))')(wgl%pGustVal(iWind,iZ,wgl%nWindPwr),iWind=1,wgl%nWinds)
!end do

  CONTAINS
  
    !==============================================================================

    real function fu_gust_integral(fWind, fZ_z0, fPwr)
      !
      ! Just integrates the normal distribution with the given parameters
      ! gust_integral = S[(1+du/u)**a * p(du/u)] d du
      ! Here du is the deviation of wind from mean value (gust)
      ! u is mean wind itself (fWind)
      ! a is wind power (fPwr)
      ! p is Gaussian function
      ! integral S is taken from minus to plus infinity.
      ! See notebook 11, p.38
      !
      implicit none

      ! Imported parameters
      real, intent(in) :: fWind, fZ_z0, fPwr

      ! Local variables
      integer :: iTmp
      real :: fInt, s, s_2, du, fTmp
      real :: f1_sqrt_2Pi = 0.39894228040143288356171412555214
      !
      ! Relative standard deviation: sigma/u = 2.185 * k / ln(Z/z0)
      ! k is von Karman = 0.4
      !
      fInt = 0.0
      fTmp = 0.0
      s = 2.185 * 0.4 * fWind / log(fZ_z0)   ! sigma
      s_2 = 0.5 / (s * s)

      do iTmp = -1000, 1000
        du = s * real(iTmp) * 0.01          ! plus-minus 10 sigmas. Should be enough...
        if(fWind+du > 0) fInt = fInt + (fWind+du)**fPwr * exp(-du*du*s_2) * 0.01 * s
        fTmp = fTmp + exp(-du*du*s_2) * 0.01 * s ! wind speed cannot be negative
      end do  ! interation cycle
 
if(abs(fTmp*f1_sqrt_2Pi/s - 1.0) > 1e-5)then
 call msg('Pure integral not unity:',fTmp*f1_sqrt_2Pi/s)
endif
      fu_gust_integral = fInt * f1_sqrt_2Pi / s

! call msg('wind, z/z0, gust forcing factor:' + fu_str(fWind), fZ_z0, fu_gust_integral/fWind**fPwr)

    end function fu_gust_integral
  
  end subroutine make_gust_lookup_table


  !**********************************************************************
  
  real function fu_effective_wind_pwr_gust(fWind, fZ_relative, fWindPower)
    !
    ! Looks into the lookup table, returns the value if finds the table, otherwise
    ! first creates it and then still returns the value.
    ! 
    ! Has to be very fast, so somewhat bulky writing to minimise forking
    !
    implicit none
    
    real, intent(in) :: fWind, fZ_relative, fWindPower
    
    ! Local variables
    integer, save :: iPwr=int_missing, iWind=int_missing, iZ=int_missing
    logical :: ifFound
    real :: zx, zy, fSum
real(r8k), save :: fMinWind=1000., fMaxWind=0., fMinZ_z0=1.d20, fMaxZ_z0=0., fMeanWind=0., fMeanZ_z0=0.
integer :: iCount = 0
    !
    ! Find the table
    !
    if(iPwr == int_missing)then  ! nothing from previous case
      ifFound = .false.
      do iPwr = 1, wgl%nWindPwr
        if(fWindPower .eps. wgl%fWindPwr(iPwr))then
          ifFound = .true.
          exit
        endif
      end do
    else
      if(fWindPower .eps. wgl%fWindPwr(iPwr))then
        ifFound = .true.    ! previous request was for the same table
      else
        ifFound = .false.
        do iPwr = 1, wgl%nWindPwr
          if(fWindPower .eps. wgl%fWindPwr(iPwr))then
            ifFound = .true.
            exit
          endif
        end do
      endif ! not from previous request
    endif ! if there was a preious request
    !
    ! If not found, create the table
    !
    if(.not. ifFound)then
      call make_gust_lookup_table(fWindPower)
      if(error)return
    endif
    !
    ! Get the value
    !
    iWind = max(1,int(fWind * 2.0))+1    ! approx from below
    iZ = max(1,int(log(fZ_relative) / log(1.5)))   ! approx from below
    !
    ! Here we take bilinear interpolation. Not the most-accurate admittedly but still...
    !
    if(iWind > wgl%nWinds-1)then
      call msg('Very strong wind; u=',fWind)
      call set_error('Very strong wind','fu_effective_wind_pwr_gust')
      fu_effective_wind_pwr_gust = 0.
      return
    endif
    if(iZ > wgl%nZ-1)then
      call msg('Very high z; z/z0=',fZ_relative)
      call set_error('Very high z','fu_effective_wind_pwr_gust')
      fu_effective_wind_pwr_gust = 0.
      return
    endif
    zx = (fWind - wgl%fWind(iWind))/(wgl%fWind(iWind+1) - wgl%fWind(iWind))
    zy = (fZ_relative - wgl%fZr(iZ))/(wgl%fZr(iZ+1) - wgl%fZr(iZ))
    fSum = (1.-zx)*(1.-zy) + zx * (1. - zy) + (1. - zx) * zy + zx * zy

    fu_effective_wind_pwr_gust = ((1.-zx)*(1.-zy) * wgl%pGustVal(iWind, iZ, iPwr) + &
                                & zx * (1. - zy)  * wgl%pGustVal(iWind+1, iZ, iPwr) + &
                                & (1. - zx) * zy  * wgl%pGustVal(iWind, iZ+1, iPwr) + &
                                & zx * zy         * wgl%pGustVal(iWind+1, iZ+1, iPwr))/fSum

if(fWind > fMaxWind) fMaxWind = fWind
if(fWind < fMinWind) fMinWind = fWind
fMeanWind = fMeanWind + fWind
if(fZ_relative > fMaxZ_z0) fMaxZ_z0 = fZ_relative
if(fZ_relative < fMinZ_z0) fMinZ_z0 = fZ_relative
fMeanZ_z0 = fMeanZ_z0 + fZ_relative
iCount = iCount + 1

if(mod(iCount,100000) == 0)then
  call msg('Min and max wind:',fMinWind, fMaxWind)
  call msg('Min and max Z/z0:',fMinZ_z0, fMaxZ_z0)
  call msg('Mean wind and Z_z0:',fMeanWind/real(iCount),fMeanZ_z0/real(iCount))
endif

  end function fu_effective_wind_pwr_gust


  ! ***************************************************************

  REAL FUNCTION fu_watervapour_satur_p(t)

    ! Description:
    ! Calculates water vapour saturation partial pressure.
    !
    ! Method:
    ! HIRLAM,  Gerard Cats  KNMI    24 September 1992
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
    REAL, INTENT(in) :: t ! temperatre [K]

    ! Local declarations:
    REAL :: zc2,zc4,zes

    REAL, PARAMETER :: c1es = 610.78, &
                     & c2es =  17.269, &
                     & c2is =  21.875, &
                     & c3es = 273.16, &
                     & c4es =  35.86, &
                     & c4is =   7.66

    zes = (t - c3es + 15.)/15.
    if(zes < 0.) zes = 0.
    if(zes > 1.) zes = 1.
    zc2 = c2is + (zes * (c2es-c2is))
    zc4 = c4is + (zes * (c4es-c4is))
    fu_watervapour_satur_p = c1es * EXP(zc2 *(t - c3es) / (t - zc4))

    !
    ! Alternative, slightly less accurate (4-th digit) method
    !
!    fu_watervapour_satur_p = exp( -6094.4692 / t + 21.1249952 - 0.027245552 * t + &
!                                & 0.000016853396 * t * t + 2.4575506 * log(t))

  END FUNCTION fu_watervapour_satur_p


  ! ***************************************************************



 ! ***************************************************************
 
  REAL FUNCTION fu_vapor_pressure_dew_point(dew_point) result(vapor_press)
    
    ! Description:
    ! 
    ! Computes the partial water vapor pressure  [Pa]    ! 
    ! Computes the partial water vapor pressure  [Pa]
    ! for the given dew point temperature  [K] 
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Ilkka Valkama, FMI
    
    IMPLICIT NONE
    
    ! Imported parameters with intent(in/out):    
    REAL, INTENT(in) :: dew_point
    
    
    !  local declarations
    REAL :: corr_a, corr_b, corr_c, corr_d
    
    !  some preliminaries
    
    IF(dew_point.LE.0.9) THEN
      vapor_press=0.0
      PRINT*,' Is dew-point in Kelvins ??'
      RETURN
    END IF
    
    ! let's start for real
    
    corr_a = 373.16/dew_point
    
    corr_b = -7.90298*(corr_a - 1.0)
    corr_b = corr_b + (5.02808*0.43429*ALOG(corr_a))
    
    corr_c = (1.0 - (1.0/dew_point))*11.344
    corr_c = -1.0 + (10.0**corr_c)
    corr_c= -1.3816*corr_c/(10.0**7)
    
    corr_d = (1.0 -  corr_a)*3.49149
    corr_d = -1.0 + (10.0**corr_d)
    corr_d = 8.1328*corr_d/(10.**3)
    
    corr_a = corr_b + corr_c + corr_d
    vapor_press = 101324.6*(10.**corr_a)
    
    
  END FUNCTION fu_vapor_pressure_dew_point
  

 ! ***************************************************************
    

  !*******************************************************************
  !*******************************************************************
  !
  !
  !    FADE-IN, FADE-OUT functions
  !
  !
  !*******************************************************************
  !*******************************************************************


    !*******************************************************************

    real function fu_fade_in(fValue_relative, fUncertainty_relative, fade_in_params)
      !
      ! Computes the fade-in function. Linear is 0 at fValue_relative = 1.-fUncertainty_relative,
      ! Grows to 0.5 at fValue_relative = 1 and grows to 1 at fVaue_relative = 1.+fUncertainty_relative
      ! Sigmoid is similar but smoother, with width decided upon the uncertainty_relative
      !
      implicit none

      ! Imported parameters
      real, intent(in) :: fUncertainty_relative, fValue_relative
      type(Tfade_inout), intent(in) :: fade_in_params

      ! Local variables
      real :: fScale, fPower

      !
      ! If the uncertainty is very low, make it a jump
      !
      if(fUncertainty_relative < 1.e-5)then
        if(fValue_relative < 1.-max(0.,fUncertainty_relative))then
          fu_fade_in = 0.
        elseif(fValue_relative > 1.+max(0.,fUncertainty_relative))then
          fu_fade_in = 1.
        else
          fu_fade_in = 0.5  ! still allow a tiny corridor between the extremes
        endif
        return
      endif

      select case(fade_in_params%iTypeIn)
        case(linear)
          fu_fade_in = min(1.,max(0.,(fValue_relative - 1. + fUncertainty_relative) / &
                                   & (2. * fUncertainty_relative)))
        case(sigmoid)
          fu_fade_in = 1. / (1. + exp(2./fUncertainty_relative * (1. - fValue_relative)))

        case(double_sigmoid)
          fScale = 2. / fUncertainty_relative
          fPower = 1. / (1. + exp(fScale * (1. - fValue_relative)))
          fu_fade_in = 1. / (1. + exp(fScale * (fPower + fade_in_params%fAsymmetryIn) * &
                                    & (1. - fValue_relative + fade_in_params%fExtraShiftIn)))
        case default
          call set_error('Unknown fade_in_method:' + fu_str(fade_in_params%iTypeIn),'fu_fade_in')
      end select

    end function fu_fade_in


    !******************************************************************

    real function fu_fade_out(fValue_relative, fUncertainty_relative, fade_out_params)
      !
      ! Computes the linear fade-in function. It is 0 at fValue_relative = 1.-fUncertainty_relative,
      ! Grows to 0.5 at fValue_relative = 1 and grows to 1 at fVaue_relative = 1.+fUncertainty_relative
      !
      implicit none

      ! Imported parameters
      real, intent(in) :: fUncertainty_relative, fValue_relative
      type(Tfade_inout), intent(in) :: fade_out_params

      ! Local variables
      real :: fScale, fPower
      
      if(fUncertainty_relative < 1.e-5)then
        if(fValue_relative < 1.-max(0.,fUncertainty_relative))then
          fu_fade_out = 1.
        elseif(fValue_relative > 1.+max(0.,fUncertainty_relative))then
          fu_fade_out = 0.
        else
          fu_fade_out = 0.5  ! still allow a tiny corridor between the extremes
        endif
        return
      endif
      select case(fade_out_params%iTypeOut)
        case(linear)
          fu_fade_out = min(1.,max(0.,(1. + fUncertainty_relative - fValue_relative) / &
                                    & (2. * fUncertainty_relative)))
        case(sigmoid)
          fu_fade_out = 1. / (1. + exp(2./fUncertainty_relative * (fValue_relative - 1.)))
        
        case(double_sigmoid)
          fScale = 2. / fUncertainty_relative
          fPower = 1. / (1. + exp(fScale * (fValue_relative - 1.)))
          fu_fade_out = 1. / (1. + exp(fScale * (fPower + fade_out_params%fAsymmetryOut) * &
                                     & (fValue_relative - 1. - fade_out_params%fExtraShiftOut)))
        case default
          call set_error('Unknown fade_out_type:' + fu_str(fade_out_params%iTypeOut),'fu_fade_out')
      end select

    end function fu_fade_out

    
    !********************************************************************
    
    subroutine report_fade_in_out_setup(params)
      implicit none
      type(Tfade_inout), intent(in) :: params

      call msg('Type of the uncertainty shape of fade-IN:' + &
           & fu_str(Params%iTypeIn) + ', and fade=OUT:' + fu_str(Params%iTypeOut))
      if(Params%iTypeIn == double_sigmoid) & 
                      & call msg('Fade-IN asymmetry parameter and compensating shoft:', &
                                       & Params%fAsymmetryIn, Params%fExtraShiftIn)
      if(Params%iTypeOut == double_sigmoid) & 
                      & call msg('Fade-OUT asymmetry parameter and compensating shift:', &
                                       & Params%fAsymmetryOut, Params%fExtraShiftOut)
    end subroutine report_fade_in_out_setup

    
    !********************************************************************
    
    subroutine set_fade_in_out_params(nlSetup, params)
      implicit none
      type(Tsilam_namelist), pointer :: nlSetup
      type(Tfade_inout), intent(out) :: params
      
      if(fu_str_u_case(fu_content(nlSetup,'fade_in_function')) == 'LINEAR')then
        Params%iTypeIn = linear
      elseif(fu_str_u_case(fu_content(nlSetup,'fade_in_function')) == 'SIGMOID')then
        Params%iTypeIn = sigmoid
      elseif(fu_str_u_case(fu_content(nlSetup,'fade_in_function')) == 'DOUBLE_SIGMOID')then
        Params%iTypeIn = double_sigmoid
        Params%fAsymmetryIn = fu_content_real(nlSetup, 'fade_in_asymmetry')
        if(fu_fails(Params%fAsymmetryIn > 0.44999, &
         & 'Strange fade_in_asymmetry:' + fu_content(nlSetup,'fade_in_asymmetry') + &
         & ', must exceed 0.45 (largest asymmetry, smoothest shape)', &
         & 'set_fade_in_out_params'))return
        Params%fExtraShiftIn = min(0., -0.0093 / (Params%fAsymmetryIn-0.04)**1.6 + 1e-3)
      else
        call set_error('Cannot understand fade_in_function:' + &
                     & fu_content(nlSetup,'uncertainty_shape_function_start') + &
                     & '; must be LINEAR, SIGMOID or DOUBLE_SIGMOID','set_fade_in_out_params')
        return
      endif 
      if(fu_str_u_case(fu_content(nlSetup,'fade_out_function')) == 'LINEAR')then
        Params%iTypeOut = linear
      elseif(fu_str_u_case(fu_content(nlSetup,'fade_out_function')) == 'SIGMOID')then
        Params%iTypeOut = sigmoid
      elseif(fu_str_u_case(fu_content(nlSetup,'fade_out_function')) == 'DOUBLE_SIGMOID')then
        Params%iTypeOut = double_sigmoid
        Params%fAsymmetryOut = fu_content_real(nlSetup, 'fade_out_asymmetry')
        if(fu_fails(Params%fAsymmetryOut > 0.39999, &
         & 'Strange fade_out_asymmetry:' + fu_content(nlSetup,'fade_out_asymmetry') + &
         & ', must exceed 0.4 (largest asymmetry, smoothest shape)', &
         & 'set_fade_in_out_params'))return
        Params%fExtraShiftOut = min(0., -0.0093 / (Params%fAsymmetryOut-0.04)**1.6 + 1e-3)
      else
        call set_error('Cannot understand fade_out_function:' + &
                     & fu_content(nlSetup,'fade_out_function') + &
                     & '; must be LINEAR, SIGMOID or DOUBLE_SIGMOID','set_fade_in_out_params')
        return
      endif
    
    end subroutine set_fade_in_out_params
    



  ! ***************************************************************
  ! ***************************************************************
  ! 
  !
  !                 TESTS
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  SUBROUTINE thermo_tests()
    
    ! Description:
    ! Test stuff for diy-weather toolbox.
    !
    ! Author: Ilkka Valkama, FMI
    ! 
    !
    IMPLICIT NONE
    !
    ! Local declarations:
!!!$    TYPE(silja_position) :: position, point
!!!$    TYPE(silja_weather_data_rules) :: weather_data_rules
!!!$    TYPE(silja_time) :: time
    INTEGER :: i,j, ik, ii
    REAL :: ddd,ff, uu, vv, rankoe
    !   SET test place and time
!!!$    time = fu_round_closest_hour(fu_wallclock(),backwards)
!!!$    position = fu_helsinki()

!!!$    point = fu_set_position(60., north, 25., east, 100000.,&
!!!$	& fu_round_closest_hour(fu_wallclock(),backwards))
!!!$    CALL update_height(point, 100.)
!!!$    position = fu_set_position(60., north, 25., east, 100000., time)
!!!$    CALL update_height(position, 100.)
!!!$      IF (error) RETURN
!!!$
!!!$
!!!$
!!!$      position = fu_set_position(60.167, north, 24.967, east, 100000., time)
!!!$    CALL update_height(position, 100.)
!!!$      IF (error) RETURN
!!!$     
!!!$      PRINT*, fu_time_string_original(time)
!!!$      PRINT *, '  Solar Elevation at this time : ',&
!!!$	  & fu_solar_elevation_angle(position)
!!!$      PRINT*,' rel. hum 850 hpa ',&
!!!$	  & fu_relative_humidity_on_level(position&
!!!$	& ,pr_level_850hpa)
!!!$      PRINT*,' rel. hum 925 hpa ',&
!!!$	  & fu_relative_humidity_on_level(position&
!!!$	& ,pr_level_925hpa)

!!!$    DO i=1,10
!!!$      PRINT*, ' Anna paine  '
!!!$      READ(*,*)ddd
!!!$      IF(ddd.Lt.0.0) STOP
!!!$      ddd=100.0*ddd
!!!$      ik=20
!!!$      DO j=0,35,5       
!!!$        ff=float(j)+273.15
!!!$    
!!!$    uu=  fu_watervap_diffusivity(ff)
!!!$    PRINT*, ' diffusivity for (',(ff-273.15),',',(ddd/100.),') is ',uu
!!!$    
!!!$        PRINT*, '   '
!!!$        ik=ik+1
!!!$        WRITE(ik,*)(ff-273.15),(ddd/100.),uu
!!!$        rankoe = fu_random_number_center(0.2,0.2)
!!!$        PRINT*, ' PS. Current random = ', rankoe
!!!$        WRITE(ik,*)rankoe
!!!$      end DO
!!!$    END DO
!!!$    
!!!$    DO i=0,36,5
!!!$      ff = FLOAT(i)
     ff = 5.
!!!$      DO ii=0,361,45
!!!$        
!!!$        ddd = FLOAT(ii)/57.29578
!!!$        CALL component_wind(ddd,ff,uu,vv)
!!!$        PRINT*, ' '
!!!$        PRINT*, ' Wind vector : ',uu,vv,' (from ',(ddd*57.29578),ff,' )'
!!!$        PRINT*, ' '
!!!$        CALL normal_wind(uu,vv,ddd,ff)
!!!$        PRINT*, ' (and the other way round = ',(ddd*57.29578),ff,' )'
!!!$      END DO
!!!$    END DO
     ! ddd,ff, uu, vv, rankoe
     DO ii=1,10
       PRINT*, ' Anna lampotila : '
       READ(*,*)uu
       PRINT*, ' '
       uu=uu+273.15
       PRINT*, ' Temperature  is ',(uu-273.15)
       PRINT*, ' and satur.vapor pressure is ',fu_saturat_watervapourpressure(uu)

!!!$       ff = fu_dew_point_from_temp(uu)
!!!$       PRINT*, ' dew point for ',(uu-273.15),' is ',ff
!!!$       vv =  fu_relhum_from_dew_point(uu,ff)
!!!$       PRINT*, ' rel hum should be  ',vv
!!!$       vv =  fu_relhum_from_dew_point(uu,ff)
!!!$       PRINT*, ' (check: ',(fu_dew_point_from_rh(uu,vv)-273.15),') '
!!!$       
!!!$       vv =  fu_saturat_watervapourpressure(uu)
!!!$       ddd =  vapor_pressure_dewp(ff)  
!!!$       PRINT*, ' and satur.vapor pressure is ',vv
!!!$       PRINT*, ' ( flexpart gives : ',ddd,') '
!!!$        fu_saturat_watervapourpressure(uu)
       PRINT*, ' '
     END DO
     
    END SUBROUTINE thermo_tests


  !******************************************************************************
   
  SUBROUTINE temp_test()

    !Description:
    !Teststuff for temperature toolbox>
    !
    !Author:
    !
    !
    IMPLICIT NONE
    !
    ! 



    PRINT*,'263.15K:', fu_watervapour_satur_p(263.15)
    PRINT*,'293.15K:', fu_watervapour_satur_p(293.15)
!!!$    PRINT*,'300K:', fu_watervapour_satur_p(300.)
!!!$    PRINT*,'310K:', fu_watervapour_satur_p(310.)
!!!$    PRINT*,'320K:', fu_watervapour_satur_p(320.)



  END SUBROUTINE temp_test

END MODULE thermodynamic_tools

