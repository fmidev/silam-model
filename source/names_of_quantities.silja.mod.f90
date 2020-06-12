MODULE names_of_quantities

  ! This module contains the names of all physical and meteorological
  ! quantities and their units, including SILAM-owned quantities.
  !
  !*****************************************************************************
  ! A REQUEST. 
  ! There is a clumsy decision here - to keep quantity as an integer and then
  ! make all the rest via functions. However, it must be realised via some
  ! structure. Functions are very easy to forget / make mistake, etc. 
  ! It is important that integer quantity as it is now can still be a sort
  ! of index in the array of these structures. And in the most ambitious case
  ! they can even be read from an external file and then somehow linked to
  ! integer parameters "quantity"
  ! This should be implemented in SILAM v.2.3 or later
  !*****************************************************************************
  !
  ! Original code: Mika Salonoja
  ! Author: Mikhail Sofiev, FMI email Mikhail.Sofiev@fmi.fi
  !
  ! All units: SI
  ! 
  ! All degrees: real numbers in degrees (degrees and parts per
  ! hundred/parts per thousand , NOT in degrees and minutes
  ! 
  ! Language: ANSI standard Fortran 90
  ! 
!  USE work_arrays
  use toolbox

  IMPLICIT NONE

  ! The public functions and subroutines available in this module:
  PUBLIC fu_quantity_string
  PUBLIC fu_quantity_short_string
  PUBLIC fu_quantity_unit
  PUBLIC fu_realtime_quantity
  PUBLIC fu_vector_quantity
  PUBLIC fu_wind_quantity
  PUBLIC fu_known_quantity
  public fu_SILAM_dispersion_quantity ! Split meteorological and dispersion quantities
  public fu_multi_level_quantity
  PUBLIC quantity_feasible_range      ! physically meaningful range for quantity
  public check_quantity_range         ! forcing the feasible range in 1D array
  PUBLIC fu_accumulated_quantity
  PUBLIC fu_aver_quantity_for_acc
  PUBLIC fu_quantity_in_quantities
  PUBLIC fu_grib_missing_real
  public fu_real_missing_replacement
  public fu_regridding_method
  public fu_get_silam_quantity ! Finds quantity flag having its name
  public fu_quantity_standard_name
  public fu_if_staggerX
  public fu_if_staggerY
  public fu_outGridInterpType

  public reformat_names_of_quantities


  public report_list_of_quantities


  ! Flag used in data retrieval, when all quantities are accepted:
  INTEGER, PARAMETER, PUBLIC :: accept_all_quantities = 259999

  ! Flag used to show that some quantity is not mandatory from the input files
  ! Usage: if some value (e.g., surface roughness) can be substituted with
  ! a pre-defined constant - it should be claimed to be a derived quantity
  ! with silam_const_flag as an initial one.
  integer, parameter, public :: silam_const_flag = 259998

  ! ***************************************************************
  !
  !
  !    Quantites available on several levels [240xxx]
  !
  !
  ! ***************************************************************

  INTEGER, PARAMETER, PUBLIC :: temperature_flag = 240000 ! in K
  INTEGER, PARAMETER, PUBLIC :: day_temperature_acc_flag = 240001 ! in K-sec
  INTEGER, PARAMETER, PUBLIC :: eq_pot_temperature_flag = 240002 !K
  INTEGER, PARAMETER, PUBLIC :: potential_temperature_flag = 240003 !K
  INTEGER, PARAMETER, PUBLIC :: day_mean_temperature_flag = 240004
  integer, parameter, public :: perturb_pot_temperature_flag = 240005  ! pot_tempr-300

  ! A trick: omega_flag and w_flag constitute two types of vertical winds, Pa/sec and m/sec.
  ! They both can be defined in various vertical coordinates. 
  ! Omega is useful in pressure/sigma/hybrid systems, while w is good in altitude
  ! above msl and height above topography relief. We also need a pointer that is the 
  ! vertical-dependent vertical wind: a vertical_velocity_flag
  !
  INTEGER, PARAMETER, PUBLIC :: u_flag = 240006   ! m/s
  INTEGER, PARAMETER, PUBLIC :: v_flag = 240007   ! m/s
  INTEGER, PARAMETER, PUBLIC :: w_alt_msl_flag = 240008   ! m/s, z-system over mean sea level
  INTEGER, PARAMETER, PUBLIC :: w_height_srf_flag = 240009   ! m/s, z-system over topography
  INTEGER, PARAMETER, PUBLIC :: omega_flag = 240010  ! Pa /s, omega wind is always in pressure system
  INTEGER, PARAMETER, PUBLIC :: vertical_velocity_flag = 240011 ! itself, does not mean anything
  INTEGER, PARAMETER, PUBLIC :: u_mean_flag = 240012
  INTEGER, PARAMETER, PUBLIC :: v_mean_flag = 240013
  INTEGER, PARAMETER, PUBLIC :: mean_wind_flag = 240014
  INTEGER, PARAMETER, PUBLIC :: dxdt_flag = 240015 ! [grid cell / sec]
  INTEGER, PARAMETER, PUBLIC :: dydt_flag = 240016 ! [grid cell / sec]
  INTEGER, PARAMETER, PUBLIC :: dzdt_flag = 240017 ! [grid cell / sec]
  INTEGER, PARAMETER, PUBLIC :: wind_flag = 240018  ! a wind vector
  INTEGER, PARAMETER, PUBLIC :: windspeed_flag = 240019  ! a windspeed scalar
  INTEGER, PARAMETER, PUBLIC :: wind_divergence_flag = 240020 ! div(wind_vector), s-1
  INTEGER, PARAMETER, PUBLIC :: wind_vertical_shear_flag = 240021 ! dU/dz, m s-1 per m

  ! mass of vater vapour per unit mass of dry air
  INTEGER, PARAMETER, PUBLIC :: pressure_flag = 240022 ! in Pa
  INTEGER, PARAMETER, PUBLIC :: cloud_cover_flag = 240023 !%
  INTEGER, PARAMETER, PUBLIC :: cloud_water_flag = 240024 !kg/kg
  INTEGER, PARAMETER, PUBLIC :: cloud_ice_flag = 240025 !kg/kg
  INTEGER, PARAMETER, PUBLIC :: cloud_cond_water_flag = 240026 !kg/kg  water+ice
  INTEGER, PARAMETER, PUBLIC :: relative_humidity_flag = 240027 ! 0...1
  INTEGER, PARAMETER, PUBLIC :: specific_humidity_flag = 240028 ! kg/kg,
  INTEGER, PARAMETER, PUBLIC :: humidity_mixing_ratio_flag = 240029 ! kg/kg,
  INTEGER, PARAMETER, PUBLIC :: relative_vorticity_flag = 240030
  INTEGER, PARAMETER, PUBLIC :: absolute_vorticity_flag = 240031
  INTEGER, PARAMETER, PUBLIC :: abs_vorticity_advection_flag = 240032
  INTEGER, PARAMETER, PUBLIC :: ipv_flag = 240033
  INTEGER, PARAMETER, PUBLIC :: tfp_flag = 240034
  INTEGER, PARAMETER, PUBLIC :: bulk_richardson_nbr_flag = 240035
  INTEGER, PARAMETER, PUBLIC :: flux_richardson_nbr_flag = 240036
  INTEGER, PARAMETER, PUBLIC :: gradient_richardson_nbr_flag = 240037
  INTEGER, PARAMETER, PUBLIC :: bulk_tfp_flag = 240038
  INTEGER, PARAMETER, PUBLIC :: turb_kinetic_energy_NWP_flag = 240039 ! [m2/s2]
  INTEGER, PARAMETER, PUBLIC :: scavenging_coefficient_flag = 240040 ! [s-1]

  INTEGER, PARAMETER, PUBLIC :: Kz_momentum_3d_flag = 240041 ! [m2 s-1]
  INTEGER, PARAMETER, PUBLIC :: Kz_heat_3d_flag = 240042 ! [m2 s-1]
  INTEGER, PARAMETER, PUBLIC :: Kz_scalar_3d_flag = 240043 ! [m2 s-1]
  INTEGER, PARAMETER, PUBLIC :: turb_kinetic_energy_SILAM_flag = 240044 ! [m2/s2]
  INTEGER, PARAMETER, PUBLIC :: turb_length_scale_flag = 240045 ! [m2/s2]
  INTEGER, PARAMETER, PUBLIC :: brunt_vaisala_freq_flag = 240046 ! [m2/s2]

  INTEGER, PARAMETER, PUBLIC :: geopotential_flag = 240047      ! m2/s2
  INTEGER, PARAMETER, PUBLIC :: height_flag = 240048 !m

  integer, parameter, public :: dtheta_dz_flag = 240049    ! K/m
  integer, parameter, public :: dpressure_dz_flag = 240050 ! Pa /m 

  INTEGER, PARAMETER, PUBLIC :: layer_thickness_flag = 240051 !m
  
  integer, parameter, public :: eta_dot_flag = 240052 ! levels/s
  INTEGER, PARAMETER, PUBLIC :: cwcabove_3d_flag = 240053 ! kg/m2 !!!Cloud water 
  INTEGER, PARAMETER, PUBLIC :: cwcolumn_flag = 240054 ! kg/m2
  INTEGER, PARAMETER, PUBLIC :: pwcabove_3d_flag = 240055 ! kg/m2  !!Precipitable watwer -- cwc minus threshold
  INTEGER, PARAMETER, PUBLIC :: pwcolumn_flag = 240056 ! kg/m2
  INTEGER, PARAMETER, PUBLIC :: R_down_meteo_flag = 240057 ! [m-1 s] -- Aerodynamic
                                                      ! resistance to previous
                                                      ! meteo level - not
                                                      ! interpolatable in
                                                      ! vertical!!!!

  INTEGER, PARAMETER, PUBLIC :: lcwcabove_3d_flag = 240058 ! kg/m2 !!!Liquid cloud water 
  INTEGER, PARAMETER, PUBLIC :: lcwcolumn_flag = 240059 ! kg/m2

  integer, parameter, public :: first_multi_level_q = 240000
  integer, parameter, public :: last_multi_level_q = 240057


  ! ***************************************************************
  !
  !
  !    Quantites which also define their height [250xxx]
  !
  !
  ! ***************************************************************
  !
  ! Temperature related quantities
  !
  INTEGER, PARAMETER, PUBLIC :: temperature_2m_flag = 250001
  INTEGER, PARAMETER, PUBLIC :: day_temperature_2m_acc_flag = 250002
  INTEGER, PARAMETER, PUBLIC :: potential_temperature_2m_flag = 250003
  INTEGER, PARAMETER, PUBLIC :: dew_point_temp_2m_flag = 250004
  INTEGER, PARAMETER, PUBLIC :: temperature_1lyr_flag = 250005 ! K tempr at 1st layer
  INTEGER, PARAMETER, PUBLIC :: underground_temperature_flag = 250006
  integer, parameter, public :: day_mean_temperature_2m_flag = 250007
  !
  ! Wind related quantities
  !
  INTEGER, PARAMETER, PUBLIC :: u_10m_flag = 250008
  INTEGER, PARAMETER, PUBLIC :: v_10m_flag = 250009
  INTEGER, PARAMETER, PUBLIC :: windspeed_10m_flag = 250010  ! Scalar.
  INTEGER, PARAMETER, PUBLIC :: wind_10m_flag = 250011       ! Similar to wind_flag: vector
  INTEGER, PARAMETER, PUBLIC :: windspeed_1lyr_flag = 250012 ! m/s Windspeed at 1st layer
  !
  ! Physiography et al
  !
  integer, parameter, public :: relief_height_flag = 250013 ! m
  integer, parameter, public :: geopotential_sfc_flag = 250014 ! m2/s2
  integer, parameter, public :: longitude_flag = 250015 ! deg
  integer, parameter, public :: latitude_flag = 250016 ! deg
  INTEGER, PARAMETER, PUBLIC :: cell_size_x_flag = 250017 ! [m]
  INTEGER, PARAMETER, PUBLIC :: cell_size_y_flag = 250018 ! [m]
  INTEGER, PARAMETER, PUBLIC :: fraction_of_ice_flag = 250019    ! 1.0 = all water is ice, 0.0 = no ice
  INTEGER, PARAMETER, PUBLIC :: fraction_of_land_flag = 250020   ! 1.0 = all land, 0.0 = all water
  INTEGER, PARAMETER, PUBLIC :: fraction_of_water_flag = 250021  ! 1.0 = all land, 0.0 = all water
  INTEGER, PARAMETER, PUBLIC :: fraction_of_forest_flag = 250022  ! 1.0 = all forest, 0.= no forest
  INTEGER, PARAMETER, PUBLIC :: fraction_of_erodible_soil_flag = 250023  ! 1.0 = all forest, 0.= no forest
  INTEGER, PARAMETER, PUBLIC :: soiltype_flag = 250024            ! index
  INTEGER, PARAMETER, PUBLIC :: albedo_flag = 250025
  INTEGER, PARAMETER, PUBLIC :: climatological_albedo_flag = 250026
  INTEGER, PARAMETER, PUBLIC :: soil_moisture_vol_frac_nwp_flag = 250027 ! m3/m3, upper layer
  integer, parameter, public :: water_salinity_flag = 250028 ! fraction, NOT %
  INTEGER, PARAMETER, PUBLIC :: land_roughness_meteo_flag = 250029 ! = land, as in meteo model
  INTEGER, PARAMETER, PUBLIC :: land_roughness_disp_flag = 250030 ! = land, refined for dispersion model
  INTEGER, PARAMETER, PUBLIC :: water_roughness_flag = 250031
  INTEGER, PARAMETER, PUBLIC :: surface_roughness_meteo_flag = 250032  ! ground merged with water
  INTEGER, PARAMETER, PUBLIC :: surface_roughness_disp_flag = 250033  ! ground merged with water
  INTEGER, PARAMETER, PUBLIC :: ground_surface_temp_flag = 250034
  INTEGER, PARAMETER, PUBLIC :: water_surface_temp_flag = 250035
  INTEGER, PARAMETER, PUBLIC :: snow_depth_flag = 250036 ! m
  INTEGER, PARAMETER, PUBLIC :: water_eq_snow_depth_flag = 250037 !kg/m2
  INTEGER, PARAMETER, PUBLIC :: charnock_parameter_flag = 250038
  INTEGER, PARAMETER, PUBLIC :: leaf_area_index_flag = 250039
  INTEGER, PARAMETER, PUBLIC :: leaf_area_indexhv_flag = 250040  !High vegetation
  INTEGER, PARAMETER, PUBLIC :: leaf_area_indexlv_flag = 250041  !Low vegetation
  INTEGER, PARAMETER, PUBLIC :: soil_sand_mass_fraction_flag = 250042
  INTEGER, PARAMETER, PUBLIC :: soil_clay_mass_fraction_flag = 250043
  INTEGER, PARAMETER, PUBLIC :: alluvial_sedim_index_flag = 250044
  INTEGER, PARAMETER, PUBLIC :: fraction_hv_flag = 250045  !High vegetation
  INTEGER, PARAMETER, PUBLIC :: fraction_lv_flag = 250046  !Low vegetation
  INTEGER, PARAMETER, PUBLIC :: land_use_type_flag = 250047
  INTEGER, PARAMETER, PUBLIC :: ref_evapotranspiration_flag = 250048 ! reference, for standard grass field
  integer, parameter, public :: water_in_soil_srf_grav_flag = 250049      ! kg/m2 in the upper layer
  integer, parameter, public :: water_in_soil_deep_grav_flag = 250050     ! kg/m2 in the deep layer
  integer, parameter, public :: water_capac_soil_srf_grav_flag = 250051  ! kg/m2 in the upper layer
  integer, parameter, public :: water_capac_soil_deep_grav_flag = 250052 ! kg/m2 in the deep layer
  integer, parameter, public :: free_soil_water_grav_flag = 250053  ! underground water excess, kg/m2
  integer, parameter, public :: emission_mask_flag = 250054  ! total emission per run per m2
  integer, parameter, public :: dust_emis_0_flag = 250055
  !
  ! Precipitation and evaporation
  !
  ! Cumulative are moved to silam-owned quantities: Meteo stack cannot make
  ! them properly!
  INTEGER, PARAMETER, PUBLIC :: large_scale_accum_rain_flag = 250060 ! mm or kg/m2
  INTEGER, PARAMETER, PUBLIC :: convective_accum_rain_flag = 250061 ! mm or kg/m2
  INTEGER, PARAMETER, PUBLIC :: large_scale_rain_int_flag = 250062 ! mm/s or kg/m2s
  INTEGER, PARAMETER, PUBLIC :: convective_rain_int_flag = 250063! mm/s or kg/m2s
  INTEGER, PARAMETER, PUBLIC :: total_precipitation_acc_flag = 250064 !kg/m2
  INTEGER, PARAMETER, PUBLIC :: total_precipitation_rate_flag = 250065 !kg/m2s
  INTEGER, PARAMETER, PUBLIC :: evaporation_flag= 250066  ! Evaporation [kg/m2]
  INTEGER, PARAMETER, PUBLIC :: integr_cloud_water_flag = 250067 !kg/m2 integrated
  INTEGER, PARAMETER, PUBLIC :: snowfall_rate_weq_flag = 250068 !kg/m2s
  integer, parameter, public :: specific_humidity_2m_flag = 250069 ! kg/kg
  integer, parameter, public :: relative_humidity_2m_flag = 250070 ! kg/kg
  !
  ! Radiation budget.
  ! sw - short wave radiation, 
  ! lw = long-wave radiation; if downward then indirect (global) radiation of the atmosphere
  ! surf = at the surface, top - at the top of the atmopshere
  ! ac = accumulated, otherwise instant flux
  ! down - downward component, net - downward minus upward flux
  !
  INTEGER, PARAMETER, PUBLIC :: surf_sw_down_radiation_ac_flag = 250075      ! [J/m2]
  INTEGER, PARAMETER, PUBLIC :: surf_sw_down_radiation_flag = 250076    ! [W/m2]
  INTEGER, PARAMETER, PUBLIC :: surf_lw_down_radiation_ac_flag = 250077      ! [J/m2]
  INTEGER, PARAMETER, PUBLIC :: surf_lw_down_radiation_flag = 250078    ! [W/m2]
  INTEGER, PARAMETER, PUBLIC :: surf_sw_net_radiation_ac_flag= 250079        ! [J/m2]
  INTEGER, PARAMETER, PUBLIC :: surf_sw_net_radiation_flag= 250080      ! [W/m2]
  INTEGER, PARAMETER, PUBLIC :: surf_lw_net_radiation_ac_flag= 250081        ! [J/m2]
  INTEGER, PARAMETER, PUBLIC :: surf_lw_net_radiation_flag= 250082      ! [W/m2]
  INTEGER, PARAMETER, PUBLIC :: top_sw_net_radiation_ac_flag = 250083        ! [J/m2]
  INTEGER, PARAMETER, PUBLIC :: top_sw_net_radiation_flag = 250084      ! [W/m2]
  INTEGER, PARAMETER, PUBLIC :: top_lw_net_radiation_ac_flag = 250085        ! [J/m2]
  INTEGER, PARAMETER, PUBLIC :: top_lw_net_radiation_flag = 250086      ! [W/m2]
  INTEGER, PARAMETER, PUBLIC :: photosynth_active_rad_flag = 250087  ! Photosynthetically active rad. flux [W/m2]
  INTEGER, PARAMETER, PUBLIC :: photosynth_active_rad_ac_flag = 250088  ! Photosynt. active rad. cumul [J/m2]
  !
  ! Surface layer energy budget
  !
  INTEGER, PARAMETER, PUBLIC :: NWP_sensible_heatflux_ac_flag= 250090  ! Accumulated sensible heat flux [J/m2]
  INTEGER, PARAMETER, PUBLIC :: NWP_sensible_heatflux_flag= 250091  ! Sensible heat flux [W/m2]
  INTEGER, PARAMETER, PUBLIC :: NWP_latent_heatflux_ac_flag= 250092  ! Accumulated latent heat flux [J/m2]
  INTEGER, PARAMETER, PUBLIC :: NWP_latent_heatflux_flag= 250093  ! Latent heat flux [W/m2]
  INTEGER, PARAMETER, PUBLIC :: SILAM_sensible_heat_flux_flag = 250094 ! surface sensible heat flux [W/m2]
  INTEGER, PARAMETER, PUBLIC :: SILAM_latent_heat_flux_flag = 250095 ! Surf latent heat flux [W/m2]
  INTEGER, PARAMETER, PUBLIC :: u_momentum_flux_flag = 250096
  INTEGER, PARAMETER, PUBLIC :: v_momentum_flux_flag = 250097
  INTEGER, PARAMETER, PUBLIC :: cape_flag = 250098
  !
  ! Surface layer scaling parameters
  !
  INTEGER, PARAMETER, PUBLIC :: pasquill_class_flag = 250100
  INTEGER, PARAMETER, PUBLIC :: MO_length_inv_flag = 250101 ! m
  INTEGER, PARAMETER, PUBLIC :: friction_velocity_flag = 250102 ! m/s
  INTEGER, PARAMETER, PUBLIC :: convective_velocity_scale_flag = 250103 ! vertical, m/s
  INTEGER, PARAMETER, PUBLIC :: temperature_scale_flag = 250104 ! K
  integer, parameter, public :: humidity_scale_flag = 250105 ! kg/kg
  integer, parameter, public :: Prandtl_nbr_flag = 250106 ! kg/kg
  INTEGER, PARAMETER, PUBLIC :: Kz_scalar_1m_flag = 250107 ! [m2 s-1]
  !
  ! ABL 
  !
  INTEGER, PARAMETER, PUBLIC :: abl_height_m_flag = 250109
  INTEGER, PARAMETER, PUBLIC :: abl_top_pressure_flag = 250110
  INTEGER, PARAMETER, PUBLIC :: nwp_abl_height_m_flag = 250111

  ! A trick. Two below pressures form the low boundary for terrain-following 
  ! and for pure pressure vertical co-ordinate systems. So,
  ! we would need a universal pointer, which will point to either of
  ! them. It should be set by grib_analysis immediately after selection
  ! of the vertical structure of the model, similar to vertical_velocity_flag
  !
  INTEGER, PARAMETER, PUBLIC :: ground_pressure_flag = 250112
  INTEGER, PARAMETER, PUBLIC :: msl_pressure_flag = 250113
  integer, parameter, public :: surface_pressure_flag = 250114 ! In fact, just a pointer
  integer, parameter, public :: log_ground_pressure_flag = 250115 ! [log(Pa)]

  INTEGER, PARAMETER, PUBLIC :: low_level_cloud_cover_flag = 250116
  INTEGER, PARAMETER, PUBLIC :: medium_level_cloud_cover_flag = 250117
  INTEGER, PARAMETER, PUBLIC :: high_level_cloud_cover_flag = 250118
  INTEGER, PARAMETER, PUBLIC :: total_cloud_cover_flag = 250119
  INTEGER, PARAMETER, PUBLIC :: sub_grid_scale_snowfall_flag =250120
  INTEGER, PARAMETER, PUBLIC :: grid_scale_snowfall_flag = 250121
  INTEGER, PARAMETER, PUBLIC :: precipitable_water_flag = 250122  ! vertically-integrated [kg/m2]

  ! Some ISBA parameters. ATTENTION. They are introduced just temporarily
  ! for the needs of fix_odditites procedure. They can not be used for 
  ! real programming because each of variables covers several true-ISBA
  ! variables
  INTEGER, PARAMETER, PUBLIC :: ISBA_temperature = 250130
  INTEGER, PARAMETER, PUBLIC :: ISBA_u_wind = 250131
  INTEGER, PARAMETER, PUBLIC :: ISBA_v_wind = 250132 
  INTEGER, PARAMETER, PUBLIC :: ISBA_spec_humidity = 250133
  INTEGER, PARAMETER, PUBLIC :: ISBA_water_eq_of_snow = 250134
  INTEGER, PARAMETER, PUBLIC :: ISBA_land_fraction = 250135
  INTEGER, PARAMETER, PUBLIC :: ISBA_land_type = 250136
  INTEGER, PARAMETER, PUBLIC :: ISBA_moisture = 250137
  INTEGER, PARAMETER, PUBLIC :: ISBA_latent_hflux = 250138
  INTEGER, PARAMETER, PUBLIC :: ISBA_sensible_hflux = 250139
  INTEGER, PARAMETER, PUBLIC :: ISBA_roughness = 250140
  
  integer, parameter, public :: r_a_flag = 250141 ! For dry deposition resistive scheme
  integer, parameter, public :: r_b_flag = 250142 ! For dry deposition resistive scheme
  integer, parameter, public :: r_s_flag = 250143 ! For dry deposition resistive scheme
  !
  ! Variables storing the characteristics of the pollen source term
  ! and determining the thresholds for start of flowering 
  !
  integer, parameter, public :: start_calday_threshold_flag = 250144  ! [day]
  integer, parameter, public :: end_calday_threshold_flag = 250145    ! [day]
  integer, parameter, public :: calday_start_end_diff_flag = 250146   ! [day]
  integer, parameter, public :: start_heatsum_threshold_flag = 250147 ! [degday/deghr/bioday etc]
  integer, parameter, public :: end_heatsum_threshold_flag = 250148   ! [degday/deghr/bioday etc]
  integer, parameter, public :: heatsum_start_end_diff_flag = 250149  ! [degday/deghr/bioday etc]
  integer, parameter, public :: growth_season_start_day_flag = 250150       ! [day]
  integer, parameter, public :: heatsum_cutoff_tempr_flag = 250151    ! [K]
  integer, parameter, public :: temperature_threshold_flag = 250152   ! [K]
  integer, parameter, public :: daily_temp_threshold_flag = 250153    ! [K]
  integer, parameter, public :: soil_moisture_threshold_flag = 250154 ! [m3/m3]
  integer, parameter, public :: pollen_total_per_m2_flag = 250155     ! [grains/m2]
  integer, parameter, public :: pollen_left_relative_flag = 250156      ! [grains/m2]
  integer, parameter, public :: pollen_correction_flag = 250157       ! [relative]
  integer, parameter, public :: plant_growth_flag = 250158            ! [relative]
  integer, parameter, public :: pollen_potency_flag= 250159           ! [pg/grain]


  ! local timezone-related flags  
  integer, parameter, public :: timezone_index_flag= 250160           ! [integer]
  
  ! dry deposition related quantities
  integer, parameter, public :: canopy_height_flag = 250161
  integer, parameter, public :: stomatal_conductance_flag = 250162 ![m/s / per unit LAI]

  ! Surface pressure to use in single-time stack in hybrid coordinates
  integer, parameter, public :: srf_press_realtime_flag = 250163
                                                        


  integer, parameter, public :: first_single_level_q = 250000
  integer, parameter, public :: last_single_level_q = 250163

  ! ***************************************************************
  !
  !
  !    SILAM-owned quantities
  !
  !
  ! ***************************************************************

  INTEGER, PARAMETER, PUBLIC :: particle_counter_flag = 260000   ! In air
  INTEGER, PARAMETER, PUBLIC :: areas_of_risk_flag = 260001 ! Vertically integrated 
  INTEGER, PARAMETER, PUBLIC :: concentration_flag = 260002 !(moles, Bq...) /m3, in air
  INTEGER, PARAMETER, PUBLIC :: mass_in_air_flag = 260003 !moles or Bq, in air
  INTEGER, PARAMETER, PUBLIC :: drydep_flag = 260004 !moles or Bq /m2, dep
  INTEGER, PARAMETER, PUBLIC :: wetdep_flag = 260005 !moles or Bq /m2, dep
  INTEGER, PARAMETER, PUBLIC :: emission_intensity_flag = 260006 ! instant emission rate (kg/s)
  INTEGER, PARAMETER, PUBLIC :: emission_flux_flag = 260007 ! instant emission flus per unit area (kg/m2s)
  INTEGER, PARAMETER, PUBLIC :: advection_moment_X_flag = 260008  ! For Galperin's advection
  INTEGER, PARAMETER, PUBLIC :: advection_moment_Y_flag = 260009  ! For Galperin's advection
  INTEGER, PARAMETER, PUBLIC :: advection_moment_Z_flag = 260010  ! For Galperin's advection
  integer, parameter, public :: heatsum_flag = 260011             ! [degree day, degree_hour, bioday]
  integer, parameter, public :: chillsum_flag = 260012             ! [degree day, degree_hour, bioday]
  integer, parameter, public :: pollen_rdy_to_fly_flag = 260013  ! [grains]
  integer, parameter, public :: allergen_rdy_to_fly_flag= 260014

  integer, parameter, public :: physiography_field_set_flag = 260015 ! Not real variable but set
  !
  ! Interpolation coefficients between the grids. Meteogrid is the basic one,
  ! but met->dispersion, met->output and dispersion->output interpolations
  ! are possible
  integer, parameter, public :: interp_met2disp_coef_flag = 260016
  integer, parameter, public :: interp_met2out_coef_flag = 260017
  integer, parameter, public :: interp_disp2out_coef_flag = 260018

  !
  ! Diagnostic quantities
  !
  integer, parameter, public :: concentration_2m_flag = 260019
  integer, parameter, public :: optical_density_flag = 260020  ! 3D quantity
  integer, parameter, public :: optical_column_depth_flag = 260021  ! 2D column-integrated
  INTEGER, PARAMETER, PUBLIC :: absorption_coef_flag = 260022   ! 3D quantity, m-1
  INTEGER, PARAMETER, PUBLIC :: scattering_coef_flag = 260023   ! 3D quantity, m-1
  INTEGER, PARAMETER, PUBLIC :: back_scattering_coef_flag = 260024   ! 3D quantity, m-1


  integer, parameter, public :: volume_mixing_ratio_flag = 260025

  integer, parameter, public :: dispersion_u_flag = 260026
  integer, parameter, public :: dispersion_v_flag = 260027
  integer, parameter, public :: dispersion_w_flag = 260028

  integer, parameter, public :: aerosol_flag = 260029

  integer, parameter, public :: cell_size_z_flag = 260030 ! m 
  integer, parameter, public :: air_density_flag = 260031

  integer, parameter, public :: emission_scaling_flag = 260032

  integer, parameter, public :: emis_factor_fire_flame_flag = 260033
  integer, parameter, public :: emis_factor_fire_smold_flag = 260034

  integer, parameter, public :: z_layertop_flag = 260035 !Metric height of disp. layer top
  integer, parameter, public :: p_layertop_flag = 260036 ! Pressure at  disp. layer top
  integer, parameter, public :: disp_flux_cellt_rt_flag = 260037 !Mass flux upwards through cell top, diagnozed (kg/s)
  integer, parameter, public :: disp_flux_celle_rt_flag = 260038 !Mass flux rightwards through cell boarder, diagnozed (kg/s)
  integer, parameter, public :: disp_flux_celln_rt_flag = 260039 !Mass flux thirthwards through cell  boarder, (kg/s)
  integer, parameter, public :: disp_flux_celleast_flag = 260040 !Mass flux rightwards through the cell right (east) side
  integer, parameter, public :: disp_flux_cellnorth_flag = 260041 !Mass flux forward through the cell front (north) side 
  integer, parameter, public :: disp_flux_celltop_flag = 260042 !Mass flux forward through the cell front (north) side 
  integer, parameter, public :: disp_cell_airmass_flag = 260043 !Mass inside a dispersion cell

  integer, parameter, public :: cloud_cond_nucley_nbr_cnc_flag = 260044  ! cloud condensation nucley, [#/m3]
  integer, parameter, public :: ice_nucley_nbr_cnc_flag = 260045  ! ice nucley, [#/m3]

  integer, parameter, public :: reaction_rate_flag = 260046 ! Reaction rate (silam amount units / sec)
  integer, parameter, public :: Vd_correction_DMAT_flag = 260047

  integer, parameter, public :: first_silam_own_q = 260000
  integer, parameter, public :: last_silam_own_q = 260047

  !
  ! This is the pointer to one of two possible vertical velocities: omega or w.
  !
  integer, public, save :: vertical_velocity_pointer = int_missing




CONTAINS


  ! ****************************************************************

  FUNCTION fu_quantity_string(quantity) result(string)

    ! Description:
    ! Returns the name of quantity in a string.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    CHARACTER(LEN=45) :: string

    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: quantity

    SELECT CASE (quantity)

      CASE (temperature_flag)
      string = 'temperature [K]'

      CASE (day_temperature_acc_flag)
      string = 'daily temperature accum [K sec]'

      CASE (day_mean_temperature_flag)
      string = 'daily mean temperature [K]'

      CASE (day_mean_temperature_2m_flag)
      string = 'daily mean temperature 2m [K]'

      CASE (potential_temperature_flag)
      string = 'potential temperature [K]'

      CASE (potential_temperature_2m_flag)
      string = 'potential temperature 2m [K]'

      CASE (perturb_pot_temperature_flag)
      string = 'perturbation pot temperature [K]'

      CASE (pressure_flag)
      string = 'pressure [Pa]'

      CASE (u_flag)
      string = 'U-wind [m/s]'

      CASE (u_mean_flag)
      string = 'U-component of mean wind [m/s]'

      CASE (v_flag)
      string = 'V-wind [m/s]'

      CASE (v_mean_flag)
      string = 'V-component of mean wind [m/s]'

      CASE (omega_flag)
      string = 'Omega-wind [Pa/s]'

      CASE (w_alt_msl_flag)
      string = 'W-wind, above msl [m/s]'

      CASE (w_height_srf_flag)
      string = 'W-wind, above srf [m/s]'

      CASE (wind_flag)
      string = 'wind-velocity vector [m/s]'

      CASE (mean_wind_flag)
      string = 'mean wind velocity vector [m/s]'

      CASE (windspeed_flag)
      string = 'windspeed [m/s]'

      CASE (geopotential_flag)
      string = 'geopotential [m2/s2]'

      CASE (geopotential_sfc_flag)
      string = 'geopotential surface [m2/s2]'

      CASE (relief_height_flag)
      string = 'relief height [m]'

      CASE (longitude_flag)
      string = 'longitude [deg]'

      CASE (latitude_flag)
      string = 'latitude [deg]'

      CASE (relative_humidity_flag)
      string = 'relative humidity [frac.]'

      CASE (specific_humidity_flag)
      string = 'specific humidity [kg/kg]'

      CASE (humidity_mixing_ratio_flag)
      string = 'humidity mixing ratio[kg/kg]'

      CASE (wind_divergence_flag)
      string = 'wind divergence [sec-1]'

      CASE (wind_vertical_shear_flag)
      string = 'wind vertical shear [sec-1]'

      CASE (albedo_flag)
      string = 'albedo [frac.]'

      CASE (large_scale_accum_rain_flag)
      string = 'accumulated large scale rain [kg/m2]'

      CASE (convective_accum_rain_flag)
      string = 'accumulated convective rain [kg/m2]'

      CASE (large_scale_rain_int_flag)
      string = 'large scale rain intensity [kg/m2s]'

      CASE (convective_rain_int_flag)
      string = 'convective rain intensity [kg/m2s]'

      CASE (soil_moisture_vol_frac_nwp_flag)
      string = 'soil moisture content [m3/m3]'

      case(soil_sand_mass_fraction_flag)
      string = 'soil sand mass fraction [kg/kg]'

      case(soil_clay_mass_fraction_flag)
      string = 'soil clay mass fraction [kg/kg]'

      case(alluvial_sedim_index_flag)
      string = 'alluvial sediment index'

      CASE (fraction_of_ice_flag)
      string = 'fraction of ice[1.0=all ice, 0.0=no ice]'

      CASE (fraction_of_land_flag)
      string = 'fraction of land [1.0=all land, 0.0=all water]'

      CASE (fraction_of_water_flag)
      string = 'fraction of water [1.0=all water, 0.0=all land]'

      CASE (surface_roughness_meteo_flag)
      string = 'surface roughness meteo [m]'

      CASE (surface_roughness_disp_flag)
      string = 'surface roughness disp [m]'

      CASE (land_roughness_meteo_flag)
      string = 'land surface roughness meteo [m]'

      CASE (land_roughness_disp_flag)
      string = 'land surface roughness disp [m]'

      CASE (water_roughness_flag)
      string = 'water surface roughness [m]'

      CASE (ground_surface_temp_flag)
      string = 'ground surface temperature [K]'

      CASE (water_surface_temp_flag)
      string = 'water surface temperature [K]'

      CASE (snow_depth_flag)
      string = 'snow depth [m]'

      CASE (fraction_of_forest_flag)
      string = 'fraction of forest [1.0=all forest, 0.0=no forest]'

      CASE (fraction_of_erodible_soil_flag)
      string = 'fraction of erodible soil [1.0=all erodible]'

      CASE (underground_temperature_flag)
      string = 'soil temperature underground[K]'

      CASE (climatological_albedo_flag)
      string = 'climatological albedo [frac.]'

      CASE (temperature_2m_flag)
      string = '2m temperature [K]'

      CASE (day_temperature_2m_acc_flag)
      string = '2m temperature accum [K sec]'

      CASE (specific_humidity_2m_flag)
      string = '2m specific humidity [kg/kg]'

      CASE (relative_humidity_2m_flag)
      string = '2m relative humidity [frac.]'

      CASE (u_10m_flag)
      string = 'U-component of 10m wind [m/s]'

      CASE (v_10m_flag)
      string = 'V-component of 10m wind [m/s]'

      CASE (windspeed_10m_flag)
      string = '10m wind speed [m/s]'

      CASE (wind_10m_flag)
      string = '10m wind vector [m/s]'

      CASE(water_eq_snow_depth_flag)
      string = 'water equiv. of snow depth [m]'

      CASE(charnock_parameter_flag)
      string = 'Charnock parameter'

      CASE (snowfall_rate_weq_flag)
      string = 'snowfall rate water equiv. [kg/m2s]'

      CASE (soiltype_flag)
      string = 'soiltype index'

      CASE (land_use_type_flag)
      string = 'landuse type'

      case(ref_evapotranspiration_flag)
        string = 'reference evapotranspiration'
        
      case(water_in_soil_srf_grav_flag)
        string = 'soil water content,surf lyr'
        
      case(water_in_soil_deep_grav_flag)
        string = 'soil water content,deep lyr'
        
      case(water_capac_soil_srf_grav_flag)
        string = 'soil water capacity,srf lyr'
        
      case(water_capac_soil_deep_grav_flag)
        string = 'soil water capacity,deep lyr'
        
      case(free_soil_water_grav_flag)
        string = 'free soil water'

      case(emission_mask_flag)
        string = 'emission mask'

      case(dust_emis_0_flag)
        string = 'dust_emission'
      
      CASE (total_precipitation_acc_flag)
      string = 'tot. precipitation [kg/m2]'

      CASE (total_precipitation_rate_flag)
      string = 'precipitation rate [kg/m2s]'

      CASE (cloud_cover_flag)
      string = 'cloud cover [fract]'

      CASE (cloud_water_flag)
      string = 'cloud water [kg/kg]'

      case (integr_cloud_water_flag)
      string = 'integr. cloud water [kg/m2]'

      CASE (cloud_ice_flag)
      string = 'cloud ice [kg/kg]'

      CASE (cloud_cond_water_flag)
      string = 'cloud condensed water [kg/kg]'

      CASE (layer_thickness_flag)
      string = 'layer thickness from bottom [m]'

      CASE(height_flag)
      string = 'level height from ground [m]'

      CASE (surf_sw_down_radiation_ac_flag)
      string = 'SW_rad surface ac.down [J/m2]'

      CASE (surf_lw_down_radiation_ac_flag)
      string = 'LW_rad surface ac.down [J/m2]'

      CASE (surf_sw_net_radiation_ac_flag)
      string = 'SW_rad surface ac.netflux [J/m2]'

      CASE (surf_lw_net_radiation_ac_flag)
      string = 'LW_rad surface ac.netflux [J/m2]'

      CASE (NWP_latent_heatflux_ac_flag)
      string = 'NWP ac.latent heat flux [J/m2]'

      CASE (NWP_sensible_heatflux_ac_flag)
      string = 'NWP ac.sensible heat flux [J/m2]'

      CASE (surf_sw_down_radiation_flag)
      string = 'SW_rad surface down [W/m2]'

      CASE (surf_lw_down_radiation_flag)
      string = 'LW_rad surface down [W/m2]'

      CASE (surf_sw_net_radiation_flag)
      string = 'SW_rad surface netflux [W/m2]'

      CASE (surf_lw_net_radiation_flag)
      string = 'LW_rad surface netflux [W/m2]'

      CASE (NWP_latent_heatflux_flag)
      string = 'NWP latent heat flux [W/m2]'

      CASE (NWP_sensible_heatflux_flag)
      string = 'NWP sensible heat flux [W/m2]'

      CASE(evaporation_flag)
      string = 'evaporation [kg/m2]'

      CASE (eq_pot_temperature_flag)
      string = 'equivivalent potential temperature [K]'

      CASE (relative_vorticity_flag)
      string = 'relative vorticity [1/s]'

      CASE (absolute_vorticity_flag)
      string = 'absolute vorticity [1/s]'

      CASE (abs_vorticity_advection_flag)
      string = 'absolute vorticity advection [1/s2]'

      CASE (abl_height_m_flag)
      string = 'boundary layer height from ground [m]'

      CASE (nwp_abl_height_m_flag)
      string = 'NWP boundary layer height from ground [m]'

      CASE (abl_top_pressure_flag)
      string = 'boundary layer height top pressure [Pa]'
      
      CASE (pasquill_class_flag)
      string = 'pasquill stability class'

      CASE (srf_press_realtime_flag)
      string = 'surface pressure realtime [Pa]'

      CASE (ground_pressure_flag)
      string = 'ground pressure [Pa]'

      CASE (msl_pressure_flag)
      string = 'mean sea level pressure [Pa]'

      CASE (surface_pressure_flag)
      string = 'vertical type-dependent surface pressure [Pa]'

      CASE (ipv_flag)
      string = 'isentropic potential vorticity [K/Pa s]'

      CASE (tfp_flag)
      string = 'thermal front parameter [K/m2]'

      CASE (bulk_richardson_nbr_flag)
      string = 'bulk Richardson number'

      CASE (flux_richardson_nbr_flag)
      string = 'flux Richardson number'

      CASE (gradient_richardson_nbr_flag)
      string = 'gradient Richardson number'

      CASE (cwcabove_3d_flag)
      string = 'Cloud condensed water column above [kg/m2]'

      CASE (cwcolumn_flag)
      string = 'Condensed water column [kg/m2]'

      CASE (lcwcabove_3d_flag)
      string = 'Cloud liquid condensed water column above [kg/m2]'

      CASE (lcwcolumn_flag)
      string = 'Liquid condensed water column [kg/m2]'

      CASE (pwcabove_3d_flag)
      string = 'Precipitable column above [kg/m2]'

      CASE (pwcolumn_flag)
      string = 'Precipitable column [kg/m2]'

      CASE (scavenging_coefficient_flag)
      string = 'scavenging coefficient [1/s]'

      CASE (vertical_velocity_flag)
      string = 'Vertical-dependent W-wind [?]'

      CASE (areas_of_risk_flag)
      string = 'area of risk'

      CASE (bulk_tfp_flag)
      string = 'bulk thermal front parameter [K/m2]'

      CASE (dew_point_temp_2m_flag)
      string = 'dew point temperature 2m [K]'

      CASE (low_level_cloud_cover_flag)
      string = 'low level cloud cover'

      CASE (medium_level_cloud_cover_flag)
      string = 'medium level cloud cover'

      CASE (high_level_cloud_cover_flag)
      string = 'high level cloud cover'

      CASE (total_cloud_cover_flag)
      string = 'total cloud cover'

      CASE (sub_grid_scale_snowfall_flag)
      string = 'sub grid-scale snowfall [kg/m2]'

      CASE (grid_scale_snowfall_flag)
      string = 'grid-scale snowfall [kg/m2]'

      CASE (precipitable_water_flag)
      string = 'precipitable water [kg/m2]'

      CASE (top_sw_net_radiation_ac_flag)
      string = 'accumulated net short wave radiation at the top of the atmosphere [J/m2]'

      CASE (top_lw_net_radiation_ac_flag)
      string = 'accumulated net long wave radiation at the top of the atmosphere [J/m2]'

      CASE (u_momentum_flux_flag)
      string = 'momentum flux ( u-component) [N/m2 s]'

      CASE (v_momentum_flux_flag)
      string = 'momentum flux ( v-component) [N/m2 s]'

      CASE (cape_flag)
      string = 'convective available potential energy [J/kg]'

      CASE (photosynth_active_rad_ac_flag)
      string = 'Photosynthetically active radiation flux cumulative [J/m2]'

      CASE (photosynth_active_rad_flag)
      string = 'Photosynthetically active radiation flux [W/m2]'

      CASE (MO_length_inv_flag)
      string = 'Inverse Monin Obukhov legth scale [m-1]'
      
      CASE (friction_velocity_flag)
      string = 'friction velocity [m/s]'

      CASE (convective_velocity_scale_flag)
      string = 'convective velocity scale [m/s]'

      CASE (humidity_scale_flag)
      string = 'humidity scale [kg/kg]'

      CASE (Prandtl_nbr_flag)
      string = 'Prandtl number'

      CASE (temperature_scale_flag)
      string = 'turbulent temperature scale [K]'
      
      CASE (SILAM_sensible_heat_flux_flag)
      string = 'SILAM sensible heat flux [W/m2]'
 
      CASE (SILAM_latent_heat_flux_flag)
      string = 'SILAM latent heat flux [W/m2]'
 
      CASE (particle_counter_flag)
      string = 'particle count'

      CASE (emission_flux_flag)
      string = 'Emission flux [massunits/m2s]'

      CASE (emission_intensity_flag)
      string = 'Emission intensity [massunits/s]'

      CASE (concentration_flag)
      string = 'Concentration in air'

      case(volume_mixing_ratio_flag)
      string = 'Volume mixing ratio'

      CASE (mass_in_air_flag)
      string = 'Mass in air'

      CASE (advection_moment_X_flag)
      string = 'Galp. X-advection moment'

      CASE (advection_moment_Y_flag)
      string = 'Galp. Y-advection moment'

      CASE (advection_moment_Z_flag)
      string = 'Galp. Z-advection moment'

      CASE (drydep_flag)
      string = 'Cocktail dry deposition'

      CASE (wetdep_flag)
      string = 'Cocktail wet deposition'

      CASE (cell_size_x_flag)
      string = 'x-size of the grid cell'

      CASE (cell_size_y_flag)
      string = 'y-size of the grid cell'

      CASE (turb_kinetic_energy_NWP_flag)
      string = 'turbulence kinetic energy, NWP'

      CASE (Kz_momentum_3d_flag)
      string = 'Kz for momentum 3d'

      CASE (Kz_heat_3d_flag)
      string = 'Kz for heat 3d'

      CASE (Kz_scalar_3d_flag)
      string = 'Kz for scalar 3d'
      
      CASE (R_down_meteo_flag)
      string = 'Resistance to prev. meteo lvl[s/m]'
      
      CASE (turb_kinetic_energy_SILAM_flag)
      string = 'turbulence kinetic energy, SILAM [m2/s2]'

      CASE (Kz_scalar_1m_flag)
      string = 'kappaustar, aka Kz_scalar_1m [m/s]'
      
      CASE (turb_length_scale_flag)
      string = 'turbulence length scale [m]'

      CASE (brunt_vaisala_freq_flag)
      string = 'Brunt-Vaisala_frequency [1/s]'

      CASE (windspeed_1lyr_flag)
      string = 'windspeed at 1st layer'

      CASE (temperature_1lyr_flag)
      string = 'temperature at 1st layer'

      CASE (log_ground_pressure_flag)
      string = 'logarithm of surface pressure'

      case (silam_const_flag)
      string = 'SILAM constant'

      case (r_a_flag)
      string = 'Ra resistance to dry deposition [s/m]'

      case (r_b_flag)
      string = 'Rb resistance to dry deposition [s/m]'

      case (r_s_flag)
      string = 'Rs resistance to dry deposition [s/m]'

      case (dtheta_dz_flag)
      string = 'd(potential_tempr) over dz [K/m]'

      case (dpressure_dz_flag)
      string = 'd(pressure) over dz [Pa/m]'

      case(ISBA_temperature)
      string = 'ISBA temperature'
      case(ISBA_u_wind)
      string = 'ISBA u_wind'
      case(ISBA_v_wind)
      string = 'ISBA v_wind'
      case(ISBA_spec_humidity)
      string = 'ISBA specific humidity'
      case(ISBA_water_eq_of_snow)
      string = 'ISBA water equiv. of snow'
      case(ISBA_land_fraction)
      string = 'ISBA land fraction'
      case(ISBA_land_type)
      string = 'ISBA land type'
      case(ISBA_moisture)
      string = 'ISBA moisture'
      case(ISBA_latent_hflux)
      string = 'ISBA latent heat flux'
      case(ISBA_sensible_hflux)
      string = 'ISBA sensible heat flux'
      case(ISBA_roughness)
      string = 'ISBA roughness'

      case(heatsum_flag)
      string = 'heatsum'
      case(chillsum_flag)
      string = 'chillsum'

      case(pollen_rdy_to_fly_flag)
      string = 'Ready to fly pollen'
      case(allergen_rdy_to_fly_flag)
      string = 'Ready to fly allergen'

      case(start_calday_threshold_flag)
      string = 'start calendar-days threshold'
      case(end_calday_threshold_flag)
      string = 'end calendarday treshold'
      case(calday_start_end_diff_flag)
      string = 'Calendar days start-end diff'

      case(start_heatsum_threshold_flag)
      string = 'start_heatsum threshold'
      case(end_heatsum_threshold_flag)
      string = 'end heatsum threshold'
      case(temperature_threshold_flag)
      string = 'temperature threshold'
      case(daily_temp_threshold_flag)
      string = 'daily temperature threshold'
      case(soil_moisture_threshold_flag)
      string = 'soil moisture threshold'
      
      case(growth_season_start_day_flag)
      string = 'Growth season start day'
      case(heatsum_cutoff_tempr_flag)
      string = 'Cutoff tempr for heatsum'
      case(heatsum_start_end_diff_flag)
      string = 'Heatsum start-end difference'

      case(pollen_correction_flag)
      string = 'Climate correction of total pollen'

      case(pollen_left_relative_flag)
      string = 'Pollen left fraction'
      
      case(plant_growth_flag)
      string = 'Plant growth'
      
      case(pollen_total_per_m2_flag)
      string = 'Pollen total per m2'
      
      case(pollen_potency_flag)
      string = 'pollen potency'

      case(timezone_index_flag)
      string = 'timezone index'

      case(physiography_field_set_flag)
      string = 'physiopgraphy set of variables'

      case (canopy_height_flag)
      string = 'Canopy height [m]'

      case (stomatal_conductance_flag)
      string = 'Stomatal conductance [m/s]'

      case(Vd_correction_DMAT_flag)
      string = 'Vd correction DMAT'

      case(water_salinity_flag)
      string = 'water salinity'

      case(interp_met2disp_coef_flag)
      string = 'interpolation meteo -> dispersion grids'

      case(interp_met2out_coef_flag)
      string = 'interpolation meteo -> output grids'

      case(interp_disp2out_coef_flag)
      string = 'interpolation dispersion -> output grids'

      case(concentration_2m_flag)
      string = 'concentration at 2m height'

      case(optical_density_flag)
      string = 'optical density'

      case(optical_column_depth_flag)
      string = 'optical column depth'

      case(absorption_coef_flag)
        string = 'absorption coefficient'

      case(scattering_coef_flag)
        string = 'scattering coefficient'

      case(back_scattering_coef_flag)
        string = 'backscattering coeffisient'

      case(dispersion_u_flag)
      string = 'Dispersion U-wind [m/s]'

      case(dispersion_v_flag)
      string = 'Dispersion V-wind [m/s]'
        
      case(dispersion_w_flag)
      string = 'Dispersion W-wind [m/s]'
      
      case(cell_size_z_flag)
      string = 'Cell vertical extent [m]'

      case(z_layertop_flag)
      string = 'Dispersion celltop height [m]'

      case(p_layertop_flag)
      string = 'Dispersion celltop pressure [Pa]'

      case(disp_flux_cellt_rt_flag)
      string = 'Dispersion celltop flux corr. [kg/s]'

      case(disp_flux_celle_rt_flag)
      string = 'Dispersion celleast flux corr. [kg/s]'

      case(disp_flux_celln_rt_flag)
      string = 'Dispersion cellnorth flux corr. [kg/s]'

      case(disp_flux_celleast_flag)
      string = 'Dispersion celleast flux (raw) [kg/s]'

      case(disp_flux_cellnorth_flag)
      string = 'Dispersion cellnorth flux (raw) [kg/s]'

      case(disp_flux_celltop_flag)
      string = 'Dispersion celltop flux (raw) [kg/s]'

      case(disp_cell_airmass_flag)
      string = 'Dispersion cell air mass [kg]'

      case(reaction_rate_flag)
      string = 'reaction rate' !silam units/s

      case(cloud_cond_nucley_nbr_cnc_flag)
        string = 'cloud condensation nucley, [nbr/m3]'

      case(ice_nucley_nbr_cnc_flag)
        string = 'ice nucley, [nbr/m3]'

      case(leaf_area_index_flag)
      string = 'Leaf area index [m2/m2]'

      case(leaf_area_indexhv_flag)
      string = 'Leaf area index high veg [m2/m2]'

      case(leaf_area_indexlv_flag)
      string = 'Leaf area index low veg [m2/m2]'

      case(fraction_hv_flag)
      string = 'Fraction  high veg [frac.]'

      case(fraction_lv_flag)
      string = 'Fraction low veg [frac.]'

      case (eta_dot_flag)
      string = 'Hybrid vertical wind'

      case (air_density_flag)
      string = 'Air density [kg/m3]'

      case (emission_scaling_flag)
      string = 'Emission scaling'

      case (emis_factor_fire_flame_flag)
      string = 'Emission factor fire flames'

      case (emis_factor_fire_smold_flag)
      string = 'Emission factor fire smould'

    CASE default
      write(unit=string,fmt='(A,I12)')'UNKNOWN QUANTITY:',quantity

    END SELECT

  END FUNCTION fu_quantity_string



  ! ****************************************************************

  FUNCTION fu_quantity_short_string(quantity) result(string)

    ! Description:
    ! Returns the sbbreviation name of quantity in a string.
    !
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    CHARACTER(LEN=15) :: string

    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: quantity

    ! Local variables
    integer :: io_status

    SELECT CASE (quantity)

      CASE (temperature_flag)
      string = 'temperature'

      CASE (day_temperature_acc_flag)
      string = 'day_tempr_acc'

      CASE (perturb_pot_temperature_flag)
      string = 'pert_pot_temp'

      CASE (day_mean_temperature_flag)
      string = 'daymean_temp'

      CASE (day_mean_temperature_2m_flag)
      string = 'daymean_temp2m'

      CASE (potential_temperature_flag)
      string = 'pot_temp'

      CASE (potential_temperature_2m_flag)
      string = 'pot_temp_2m'

      CASE (pressure_flag)
      string = 'pressure'

      CASE (u_flag)
      string = 'U_wind'

      CASE (u_mean_flag)
      string = 'mean_U_wind'

      CASE (v_flag)
      string = 'V_wind'

      CASE (v_mean_flag)
      string = 'mean_V_wind'

      CASE (omega_flag)
      string = 'omega_wind'

      CASE (w_alt_msl_flag)
      string = 'w_msl'

      CASE (w_height_srf_flag)
      string = 'w_srf'

      CASE (wind_flag)
      string = 'windvel_vect'

      CASE (mean_wind_flag)
      string = 'mean_windvel_v'

      CASE (windspeed_flag)
      string = 'windspeed'

      CASE (geopotential_flag)
      string = 'geopotential'

      CASE (geopotential_sfc_flag)
      string = 'geopot_sfc'

      CASE (relief_height_flag)
      string = 'relief_height'

      CASE (longitude_flag)
      string = 'longitude'

      CASE (latitude_flag)
      string = 'latitude'

      CASE (relative_humidity_flag)
      string = 'rel_humid'

      CASE (specific_humidity_flag)
      string = 'spec_humid'

      CASE (humidity_mixing_ratio_flag)
      string = 'humid_mix_rat'

      CASE (wind_divergence_flag)
      string = 'wind_div'

      CASE (wind_vertical_shear_flag)
      string = 'dU_dz'

      CASE (albedo_flag)
      string = 'albedo'

      CASE (large_scale_accum_rain_flag)
      string = 'acc_ls_rain'

      CASE (convective_accum_rain_flag)
      string = 'acc_conv_rain'

      CASE (large_scale_rain_int_flag)
      string = 'ls_rain_inten'

      CASE (convective_rain_int_flag)
      string = 'cnv_rain_inten'

      CASE (soil_moisture_vol_frac_nwp_flag)
      string = 'soil_moisture'

      case(soil_sand_mass_fraction_flag)
      string = 'soil_sand_fract'

      case(soil_clay_mass_fraction_flag)
      string = 'soil_clay_fract'
      
      case(alluvial_sedim_index_flag)
      string = 'alluv_sedim_ind'
      
      CASE (fraction_of_ice_flag)
      string = 'ice_fract'

      CASE (fraction_of_land_flag)
      string = 'land_fract'

      CASE (fraction_of_water_flag)
      string = 'water_fract'

      CASE (surface_roughness_meteo_flag)
      string = 'srf_rough_met'

      CASE (surface_roughness_disp_flag)
      string = 'srf_rough_disp'

      CASE (land_roughness_meteo_flag)
      string = 'land_rough_met'

      CASE (land_roughness_disp_flag)
      string = 'land_rough_disp'

      CASE (water_roughness_flag)
      string = 'water_rough'

      CASE (ground_surface_temp_flag)
      string = 'gr_srf_temp'

      CASE (water_surface_temp_flag)
      string = 'water_srf_temp'

      CASE (snow_depth_flag)
      string = 'snow_depth'

      CASE (fraction_of_forest_flag)
      string = 'forest_fract'

      CASE (fraction_of_erodible_soil_flag)
      string = 'erodible_fract'

      CASE (underground_temperature_flag)
      string = 'soil_temp_u_gr'

      CASE (climatological_albedo_flag)
      string = 'clim_albedo'

      CASE (temperature_2m_flag)
      string = 'temp_2m'

      CASE (day_temperature_2m_acc_flag)
      string = 'temp_2m_acc'

      CASE (specific_humidity_2m_flag)
      string = 'spec_humid_2m'

      CASE (relative_humidity_2m_flag)
      string = 'relat_humid_2m'

      CASE (u_10m_flag)
      string = 'U_wind_10m'

      CASE (v_10m_flag)
      string = 'V_wind_10m'

      CASE (windspeed_10m_flag)
      string = 'windspeed_10m'

      CASE (wind_10m_flag)
      string = 'wind_10m'

      CASE(water_eq_snow_depth_flag)
      string = 'weq_snow_depth'

      CASE(charnock_parameter_flag)
      string = 'charnock'

      CASE (snowfall_rate_weq_flag)
      string = 'snowfall_weq'

      CASE (soiltype_flag)
      string = 'soiltype'

      CASE (land_use_type_flag)
      string = 'landuse'
      
      case(ref_evapotranspiration_flag)
        string = 'ref_evapotr'
        
      case(water_in_soil_srf_grav_flag)
        string = 'soil_water_surf'
        
      case(water_in_soil_deep_grav_flag)
        string = 'soil_water_deep'
        
      case(water_capac_soil_srf_grav_flag)
        string = 'soil_water_cap_srf'
        
      case(water_capac_soil_deep_grav_flag)
        string = 'soil_water_cap_deep'
        
      case(free_soil_water_grav_flag)
        string = 'free_soil_water'
        
      case(emission_mask_flag)
        string = 'emis_mask'

      case(dust_emis_0_flag)
        string = 'dust_emis_0'

      CASE (total_precipitation_acc_flag)
      string = 'tot_prec'
      
      CASE (total_precipitation_rate_flag)
      string = 'prec_rate'

      CASE (cloud_cover_flag)
      string = 'cloud_cover'

      CASE (cloud_water_flag)
      string = 'cloud_water'

      CASE (integr_cloud_water_flag)
      string = 'intg_cld_water'

      CASE (cloud_ice_flag)
      string = 'cloud_ice'

      CASE (cloud_cond_water_flag)
      string = 'cloud_cond_water'

      CASE (layer_thickness_flag)
      string = 'lyr_thick'

      CASE(height_flag)
      string = 'lev_height'

      CASE (surf_sw_down_radiation_ac_flag)
      string = 'SW_dn_srf_acc'

      CASE (surf_lw_down_radiation_ac_flag)
      string = 'LW_dn_srf_acc'

      CASE (surf_sw_net_radiation_ac_flag)
      string = 'SW_net_srf_acc'

      CASE (surf_lw_net_radiation_ac_flag)
      string = 'LW_net_srf_acc'

      CASE (NWP_latent_heatflux_ac_flag)
      string = 'nwp_acc_lat_hfl'

      CASE (NWP_sensible_heatflux_ac_flag)
      string = 'nwp_acc_sns_hfl'

      CASE (surf_sw_down_radiation_flag)
      string = 'SW_int_dn_srf'

      CASE (surf_lw_down_radiation_flag)
      string = 'LW_int_dn_srf'

      CASE (surf_sw_net_radiation_flag)
      string = 'SW_int_net_srf'

      CASE (surf_lw_net_radiation_flag)
      string = 'LW_int_net_srf'

      CASE (NWP_latent_heatflux_flag)
      string = 'nwp_lat_hflux'

      CASE (NWP_sensible_heatflux_flag)
      string = 'nwp_sens_hflux'

      CASE(evaporation_flag)
      string = 'evaporation'

      CASE (eq_pot_temperature_flag)
      string = 'equiv_pot_temp'

      CASE (relative_vorticity_flag)
      string = 'relat_vort'

      CASE (absolute_vorticity_flag)
      string = 'abs_vort'

      CASE (abs_vorticity_advection_flag)
      string = 'abs_vort_adv'

      CASE (abl_height_m_flag)
      string = 'BLH'

      CASE (nwp_abl_height_m_flag)
      string = 'nwp_BLH'

      CASE (abl_top_pressure_flag)
      string = 'BL_top_pr'
      
      CASE (pasquill_class_flag)
      string = 'pasquill_stab'

      CASE (ground_pressure_flag)
      string = 'gr_surf_pr'

      CASE (msl_pressure_flag)
      string = 'sea_lev_pr'

      CASE (surface_pressure_flag)
      string = 'vrt_dep_srf_pr'

      CASE (ipv_flag)
      string = 'isentr_pot_vrt'

      CASE (tfp_flag)
      string = 'therm_frnt_par'

      CASE (bulk_richardson_nbr_flag)
      string = 'bulk_Ri_nbr'

      CASE (flux_richardson_nbr_flag)
      string = 'flux_Ri_nbr'

      CASE (gradient_richardson_nbr_flag)
      string = 'grad_Ri_nbr'

      CASE (cwcabove_3d_flag)
      string = 'cwcabove_3d'
      
      CASE (cwcolumn_flag)
      string = 'cwcolumn'
      
      CASE (lcwcabove_3d_flag)
      string = 'liquid_cwcabove_3d'

      CASE (lcwcolumn_flag)
      string = 'liquid_cwcolumn'

      CASE (pwcabove_3d_flag)
      string = 'pwcabove_3d'
      
      CASE (pwcolumn_flag)
      string = 'pwcolumn'
      
      CASE (scavenging_coefficient_flag)
      string = 'scav_coef'

      CASE (vertical_velocity_flag)
      string = 'vert_dep_wind'

      CASE (turb_length_scale_flag)
      string = 'turb_l_scale'

      CASE (brunt_vaisala_freq_flag)
      string = 'brunt_vais_f'

      CASE (areas_of_risk_flag)
      string = 'area_of_risk'

      CASE (bulk_tfp_flag)
      string = 'bulk_tf_par'

      CASE (dew_point_temp_2m_flag)
      string = 'dew_p_temp2m'

      CASE (low_level_cloud_cover_flag)
      string = 'low_lev_cloud'

      CASE (medium_level_cloud_cover_flag)
      string = 'med_lev_cloud'

      CASE (high_level_cloud_cover_flag)
      string = 'high_lev_cloud'

      CASE (total_cloud_cover_flag)
      string = 'total_cloud'

      CASE (sub_grid_scale_snowfall_flag)
      string = 'subgr_snowfall'

      CASE (grid_scale_snowfall_flag)
      string = 'grid_snowfall'

      CASE (precipitable_water_flag)
      string = 'prec_water'

      CASE (top_sw_net_radiation_ac_flag)
      string = 'acc_net_SW_TOA'

      CASE (top_lw_net_radiation_ac_flag)
      string = 'acc_net_LW_TOA'

      CASE (u_momentum_flux_flag)
      string = 'mom_flux'

      CASE (v_momentum_flux_flag)
      string = 'mom_flux'

      CASE (cape_flag)
      string = 'cape'

      CASE (photosynth_active_rad_flag)
      string = 'photosyn_act_rad'

      CASE (photosynth_active_rad_ac_flag)
      string = 'photosyn_act_rad_ac'

      CASE (MO_length_inv_flag)
      string = 'MO_len_inv'
      
      CASE (friction_velocity_flag)
      string = 'fric_vel'

      CASE (convective_velocity_scale_flag)
      string = 'cnv_vel_scale'

      CASE (temperature_scale_flag)
      string = 'turb_temp'
      
      CASE (humidity_scale_flag)
      string = 'humid_scale'
      
      CASE (Prandtl_nbr_flag)
      string = 'Prandtl_nbr'
      
      CASE (SILAM_sensible_heat_flux_flag)
      string = 'SILAM_sens_hfl'
 
      CASE (SILAM_latent_heat_flux_flag)
      string = 'SILAM_lat_hfl'
 
      CASE (dtheta_dz_flag)
      string = 'dtheta_dz'
 
      CASE (dpressure_dz_flag)
      string = 'dpressure_dz'
 
      CASE (particle_counter_flag)
      string = 'part_count'

      CASE (emission_flux_flag)
      string = 'emf'

      CASE (emission_intensity_flag)
      string = 'ems'

      CASE (concentration_flag)
      string = 'cnc'

      case(volume_mixing_ratio_flag)
      string = 'vmr'

      CASE (mass_in_air_flag) 
      string = 'mass_in_air'

      CASE (advection_moment_X_flag) 
      string = 'adv_moment_x'

      CASE (advection_moment_Y_flag) 
      string = 'adv_moment_y'

      CASE (advection_moment_Z_flag) 
      string = 'adv_moment_z'

      CASE (drydep_flag)
      string = 'dd'

      CASE (wetdep_flag)
      string = 'wd'

      CASE (cell_size_x_flag)
      string = 'x_size_cell'

      CASE (cell_size_y_flag)
      string = 'y_size_cell'

      CASE (turb_kinetic_energy_NWP_flag)
      string = 'TKE_NWP'

      CASE (turb_kinetic_energy_SILAM_flag)
      string = 'TKE_SILAM'

      CASE (Kz_momentum_3d_flag)
      string = 'Kz_momentum'

      CASE (Kz_heat_3d_flag)
      string = 'Kz_heat'

      CASE (Kz_scalar_3d_flag)
      string = 'Kz_scalar'

      CASE (R_down_meteo_flag)
      string = 'R_down_met'

      CASE (Kz_scalar_1m_flag)
      string = 'Kz_1m'

      CASE (windspeed_1lyr_flag)
      string = 'windsp_1lyr'

      CASE (temperature_1lyr_flag)
      string = 'tempr_1lyr'

      CASE (log_ground_pressure_flag)
      string = 'log_surf_pr'

      case(silam_const_flag)
      string = 'silam_const'

      case (r_a_flag)
      string = 'Ra_res_drydep'

      case (r_b_flag)
      string = 'Rb_res_drydep'

      case (r_s_flag)
      string = 'Rs_res_drydep'

      case (canopy_height_flag)
      string = 'h_canopy'

      case (stomatal_conductance_flag)
      string = 'g_stomatal'

      case(ISBA_temperature)
      string = 'IS_tempr'
      case(ISBA_u_wind)
      string = 'IS_u_wind'
      case(ISBA_v_wind)
      string = 'IS_v_wind'
      case(ISBA_spec_humidity)
      string = 'IS_spec_humid'
      case(ISBA_water_eq_of_snow)
      string = 'IS_water_eq_sn'
      case(ISBA_land_fraction)
      string = 'IS_land_fr'
      case(ISBA_land_type)
      string = 'IS_land_type'
      case(ISBA_moisture)
      string = 'IS_moisture'
      case(ISBA_latent_hflux)
      string = 'IS_lat_hflux'
      case(ISBA_sensible_hflux)
      string = 'IS_sens_hflux'
      case(ISBA_roughness)
      string = 'IS_roughn'

      case(heatsum_flag)
      string = 'heatsum'
      case(chillsum_flag)
      string = 'chillsum'

      case(pollen_rdy_to_fly_flag)
      string = 'Poll_Rdy2fly'
      case(allergen_rdy_to_fly_flag)
      string = 'alrg_Rdy2fly'
      
      case(start_calday_threshold_flag)
      string = 'start_calday_th'
      case(end_calday_threshold_flag)
      string = 'end_calday_th'
      case(calday_start_end_diff_flag)
      string = 'Calday_beg_end'

      case(start_heatsum_threshold_flag)
      string = 'start_hsum_th'
      case(temperature_threshold_flag)
      string = 'temp_th'
      case(daily_temp_threshold_flag)
      string = 'daily_t_th'
      case(soil_moisture_threshold_flag)
      string = 'soilmoist_th'
      
      case(growth_season_start_day_flag)
      string = 'growth_start_dy'
      case(heatsum_cutoff_tempr_flag)
      string = 'hsum_cutoff_t'
      case(end_heatsum_threshold_flag)
      string = 'end_hsum_th'
      case(heatsum_start_end_diff_flag)
      string = 'hsum_start_end'

      case(pollen_correction_flag)
      string = 'pollen_corr'

      case(pollen_left_relative_flag)
      string = 'poll_left'
      
      case(pollen_total_per_m2_flag)
      string = 'poll_tot_m2'
      
      case(plant_growth_flag)
      string = 'plant_growth'
      
      case(pollen_potency_flag)
      string = 'poll_pot'

      case( timezone_index_flag)
      string = 'tz_index'

      case(physiography_field_set_flag)
      string = 'physiopgr_set'

      case(Vd_correction_DMAT_flag)
      string = 'Vd_corr_DMAT'

      case(water_salinity_flag)
      string = 'water_salinity'

      case(interp_met2disp_coef_flag)
      string = 'met2disp_intrp'

      case(interp_met2out_coef_flag)
      string = 'met2out_intrp'

      case(interp_disp2out_coef_flag)
      string = 'disp2out_intrp'

      case(concentration_2m_flag)
      string = 'cnc2m'

      case(optical_density_flag)
      string = 'od'

      case(optical_column_depth_flag)
      string = 'ocd'

      case(absorption_coef_flag)
        string = 'absorp_coef'

      case(scattering_coef_flag)
        string = 'scatt_coef'

      case(back_scattering_coef_flag)
        string = 'backscatt_coef'

      case(dispersion_u_flag)
      string = 'disp_u_wind'

      case(dispersion_v_flag)
      string = 'disp_v_wind'
        
      case(dispersion_w_flag)
      string = 'disp_w_wind'
      
      case(cell_size_z_flag)
      string = 'cellsize_z'

      case(z_layertop_flag)
      string = 'disp_ztop'

      case(p_layertop_flag)
      string = 'disp_ptop'

      case(disp_flux_cellt_rt_flag)
      string = 'd_flux_cellt_rt'

      case(disp_flux_celle_rt_flag)
      string = 'd_flux_celle_rt'

      case(disp_flux_celln_rt_flag)
      string = 'd_flux_celln_rt'

      case(disp_flux_celleast_flag)
      string = 'd_flux_celle'

      case(disp_flux_cellnorth_flag)
      string = 'd_flux_celln'

      case(disp_flux_celltop_flag)
      string = 'd_flux_cellt'

      case(disp_cell_airmass_flag)
      string = 'd_cellmass'

      case(reaction_rate_flag)
      string = 'react_rate'

      case(cloud_cond_nucley_nbr_cnc_flag)
        string = 'ccn'

      case(ice_nucley_nbr_cnc_flag)
        string = 'in'

      case(leaf_area_index_flag)
      string = 'lai'

      case(leaf_area_indexhv_flag)
      string = 'lai_hv'

      case(leaf_area_indexlv_flag)
      string = 'lai_lv'

      case(fraction_hv_flag)
      string = 'frac_hv'

      case(fraction_lv_flag)
      string = 'frac_lv'

      case(eta_dot_flag)
      string = 'eta_dot'

      case(air_density_flag)
      string = 'air_dens'

      case(emission_scaling_flag)
      string = 'ems_scale'

      case(emis_factor_fire_flame_flag)
      string = 'emfac_fire_flm'

      case(emis_factor_fire_smold_flag)
      string = 'emfac_fire_smld'

    CASE default
      write(unit=string,fmt='(i15)', iostat = io_status) quantity
      string = adjustl(string)
      string = fu_connect_strings(string, '-UNKNOWN')

    END SELECT

  END FUNCTION fu_quantity_short_string


  ! ****************************************************************

  FUNCTION fu_quantity_unit(quantity, massunit_) result(string)

    ! Description:
    ! Returns the unit of quantity in a readable string.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    CHARACTER(LEN=20) :: string

    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: quantity
    character( len=*), intent(in), optional :: massunit_

    character( len=10) :: massunit

    !Replace crazy massunit with something more reasonable
    if (present(massunit_)) then 
         massunit = massunit_
      else
         massunit = "mole,Bq,kg"
    endif


    SELECT CASE (quantity)

      CASE (temperature_flag, temperature_2m_flag, &
          & day_mean_temperature_flag, day_mean_temperature_2m_flag, &
          & dew_point_temp_2m_flag, &
          & potential_temperature_flag, potential_temperature_2m_flag, &
          & perturb_pot_temperature_flag, &
          & temperature_1lyr_flag, &
          & eq_pot_temperature_flag, &
          & ground_surface_temp_flag, &
          & water_surface_temp_flag)
        string = 'K'

      CASE (day_temperature_acc_flag, day_temperature_2m_acc_flag)
        string = 'K sec'

      CASE (pressure_flag)
        string = 'Pa'

      CASE (u_flag,&
          & v_flag,&
          & u_mean_flag,&
          & v_mean_flag,&
          & w_alt_msl_flag, w_height_srf_flag, &
          & wind_flag,&
          & mean_wind_flag, &
          & windspeed_flag, &
          & dispersion_u_flag,&
          & dispersion_v_flag)
       string = 'm/s' 

      CASE (dispersion_w_flag)
        string = 'may be m/s'

      CASE (omega_flag)
        string = 'Pa/s'

      CASE (geopotential_flag)
        string = 'm2/s2'

      CASE (geopotential_sfc_flag)
        string = 'm2/s2'

      CASE (relief_height_flag,&
          & cell_size_z_flag)
        string = 'm'

      CASE (longitude_flag)
        string = 'deg'

      CASE (latitude_flag)
        string = 'deg'

      CASE (relative_humidity_flag)
        string = 'frac.'

      CASE (specific_humidity_flag)
        string = 'kg/kg'

      CASE (humidity_mixing_ratio_flag)
        string = 'kg/kg'

      CASE (wind_divergence_flag)
        string = 'sec**-1'

      CASE (wind_vertical_shear_flag)
        string = 'sec**-1'

      CASE (albedo_flag, &
          & low_level_cloud_cover_flag, &
          & medium_level_cloud_cover_flag, &
          & high_level_cloud_cover_flag, &
          & total_cloud_cover_flag, & 
          & cloud_cover_flag) 
        string = 'frac.'

      CASE (large_scale_accum_rain_flag)
        string = 'kg/m2'

      CASE (convective_accum_rain_flag)
        string = 'kg/m2'

      CASE (large_scale_rain_int_flag)
        string = 'kg/m2s'

      CASE (convective_rain_int_flag)
        string = 'kg/m2s'

      CASE (soil_moisture_vol_frac_nwp_flag)
        string = 'm3/m3'

      case(soil_sand_mass_fraction_flag)
      string = 'kg/kg'

      case(soil_clay_mass_fraction_flag)
      string = 'kg/kg'

      case(alluvial_sedim_index_flag)
      string = ''

      CASE (fraction_of_ice_flag)
        string = 'frac.'

      CASE (fraction_of_land_flag, fraction_of_water_flag)
        string = 'frac.'

      CASE (surface_roughness_meteo_flag)
        string = 'm'

      CASE (surface_roughness_disp_flag)
        string = 'm'

      CASE (land_roughness_meteo_flag)
        string = 'm'

      CASE (land_roughness_disp_flag)
        string = 'm'

      CASE (water_roughness_flag)
        string = 'm'

      CASE (snow_depth_flag)
        string = 'm'

      CASE (fraction_of_forest_flag)
        string = 'frac.'

      CASE (fraction_of_erodible_soil_flag)
        string = 'frac.'

      CASE (underground_temperature_flag)
        string = 'K'

      CASE (climatological_albedo_flag)
        string = 'frac.'

      CASE (water_eq_snow_depth_flag)
        string = 'm'

      CASE (charnock_parameter_flag)
        string = ' '

      CASE (snowfall_rate_weq_flag)
        string = 'kg/m2s'

      CASE (soiltype_flag)
        string = ' '

      CASE (land_use_type_flag)
        string = ' '
        
      case(ref_evapotranspiration_flag)
        string = 'kg/m2sec'
        
      case(water_in_soil_srf_grav_flag)
        string = 'kg/m2'
        
      case(water_in_soil_deep_grav_flag)
        string = 'kg/m2'
        
      case(water_capac_soil_srf_grav_flag)
        string = 'kg/m2'
        
      case(water_capac_soil_deep_grav_flag)
        string = 'kg/m2'
        
      case(free_soil_water_grav_flag)
        string = 'kg/m2'
        
      case(emission_mask_flag)
        string = massunit+'/m2'

      case(dust_emis_0_flag)
        string = ''

      CASE (total_precipitation_acc_flag)
        string = 'kg/m2'
        
      CASE (total_precipitation_rate_flag)
        string = 'kg/m2s'

      CASE (integr_cloud_water_flag)
        string = 'kg/m2'

      CASE (cloud_ice_flag, cloud_water_flag, cloud_cond_water_flag)
        string = 'kg/kg'

      CASE (layer_thickness_flag, height_flag, cell_size_x_flag, cell_size_y_flag)
        string = 'm'

      CASE (evaporation_flag)
        string = 'kg/m2'

      CASE (relative_vorticity_flag,&
          & absolute_vorticity_flag)
        string = '1/s'

      CASE (abs_vorticity_advection_flag)
        string = 's**-2'

      CASE (abl_height_m_flag)
        string = 'm'

      CASE (nwp_abl_height_m_flag)
        string = 'm'

      CASE (abl_top_pressure_flag)
        string = 'Pa'

      CASE (surf_sw_down_radiation_ac_flag,&
          & surf_lw_down_radiation_ac_flag,&
          & surf_sw_net_radiation_ac_flag,&
          & surf_lw_net_radiation_ac_flag,&
          & NWP_latent_heatflux_ac_flag,&
          & NWP_sensible_heatflux_ac_flag)
        string = 'J/m2'

      CASE (surf_sw_down_radiation_flag,&
          & surf_lw_down_radiation_flag,&
          & surf_sw_net_radiation_flag,&
          & surf_lw_net_radiation_flag,&
          & NWP_latent_heatflux_flag,&
          & NWP_sensible_heatflux_flag, &
          & SILAM_sensible_heat_flux_flag, &
          & SILAM_latent_heat_flux_flag, &
          & photosynth_active_rad_flag)
        string = 'W/m2'

      CASE (ipv_flag)
        string = 'K/Pa*s'

      CASE (tfp_flag)
        string = 'K/m2 '

      CASE (ground_pressure_flag, msl_pressure_flag)
        string = 'Pa'

      CASE (cwcabove_3d_flag,cwcolumn_flag,pwcabove_3d_flag,pwcolumn_flag,lcwcabove_3d_flag,lcwcolumn_flag)
        string = 'kg/m2'

      CASE (scavenging_coefficient_flag)
        string = '1/s'

      CASE (vertical_velocity_flag)
        string = '??' 

      CASE (bulk_tfp_flag)
        string = 'K/m2 '

      CASE (sub_grid_scale_snowfall_flag, &
          & grid_scale_snowfall_flag, &
          & precipitable_water_flag)
        string = 'kg/m2'
 
      CASE (top_sw_net_radiation_ac_flag, &
          & top_lw_net_radiation_ac_flag, &
          & photosynth_active_rad_ac_flag)
        string = 'J/m2'

      CASE (u_momentum_flux_flag,& 
          & v_momentum_flux_flag )
        string = 'N/m2 s'

      CASE (cape_flag)
        string = 'J/kg'
 
      CASE (MO_length_inv_flag)
        string = ' m-1 '
      
      CASE (friction_velocity_flag, convective_velocity_scale_flag)
        string = 'm/s'
 
      CASE (humidity_scale_flag)
        string = 'kg/kg'
 
      CASE (dtheta_dz_flag)
        string = 'K/m'
 
      CASE (dpressure_dz_flag)
        string = 'Pa/m'
 
      CASE (particle_counter_flag)
        string = 'number'

      CASE (emission_flux_flag)
        string = massunit+'/m2s'

      CASE (emission_intensity_flag)
        string = massunit+'/sec'

      CASE (concentration_flag)
        string = massunit+'/m3'

      case(volume_mixing_ratio_flag)
        string = 'mole/mole'

      CASE (mass_in_air_flag)
        string = massunit

      CASE (advection_moment_X_flag, advection_moment_Y_flag, advection_moment_Z_flag)
        string = massunit

      CASE (drydep_flag)
        string = massunit+'/m2/s'

      CASE (wetdep_flag)
        string = massunit+'/m2/s'

      CASE (turb_kinetic_energy_NWP_flag, turb_kinetic_energy_SILAM_flag)
        string = 'J/m3'

      CASE (Kz_momentum_3d_flag, Kz_heat_3d_flag, Kz_scalar_3d_flag, Kz_scalar_1m_flag)
        string = 'm2/s'

      CASE (R_down_meteo_flag)
        string = 's/m'

      CASE (windspeed_1lyr_flag)
        string = 'm/s'

      CASE (log_ground_pressure_flag)
        string = 'log(Pa)'

      case (r_a_flag, r_b_flag, r_s_flag)
        string = 's/m'
      
      case (canopy_height_flag)
      string = 'm'

      case (stomatal_conductance_flag)
      string = 'm/s'

      case(heatsum_flag)
        string = 'deg_time'

      case(chillsum_flag)
        string = 'deg_time'

      case(pollen_rdy_to_fly_flag)
        string = 'number'

      case(allergen_rdy_to_fly_flag)
        string = 'ng'

      case(start_calday_threshold_flag)
        string = 'day'

      case(end_calday_threshold_flag)
        string = 'day'

      case(calday_start_end_diff_flag)
        string = 'day'

      case(start_heatsum_threshold_flag)
        string = 'deg_time'

      case(temperature_threshold_flag)
        string = 'K'
      case(daily_temp_threshold_flag)
        string = 'K'
        
      case(soil_moisture_threshold_flag)
        string = 'm3/m3'


      case(growth_season_start_day_flag)
        string = 'number'

      case(heatsum_cutoff_tempr_flag)
        string = 'K'

      case(end_heatsum_threshold_flag)
        string = 'deg_time'

      case(heatsum_start_end_diff_flag)
        string = 'deg_time'

      case(physiography_field_set_flag)
        string=''

      case(Vd_correction_DMAT_flag)
        string=''

      case(water_salinity_flag)
        string = 'kg/kg'

      case(concentration_2m_flag)
        string = massunit+'/m3'

      case(optical_density_flag)
        string = 'm2/m2'

      case(optical_column_depth_flag)
        string = ''

      case(absorption_coef_flag)
        string = '1/m'

      case(scattering_coef_flag)
        string = '1/m'

      case(back_scattering_coef_flag)
        string = '1/m'

      case(z_layertop_flag)
      string = 'm'

      case(p_layertop_flag)
      string = 'Pa'

      case(disp_flux_cellt_rt_flag,disp_flux_celle_rt_flag,disp_flux_celln_rt_flag)
      string = 'kg/s'

      case(disp_flux_celleast_flag, disp_flux_cellnorth_flag, disp_flux_celltop_flag)
      string = 'kg/s'

      case(disp_cell_airmass_flag)
      string = 'kg'

      case(reaction_rate_flag)
      string = 'massunit/s'

      case(cloud_cond_nucley_nbr_cnc_flag)
        string = '[nbr/m3]'

      case(ice_nucley_nbr_cnc_flag)
        string = '[nbr/m3]'

      case(eta_dot_flag)
        string = 'may be 1/s'
        
      case(air_density_flag)
        string = 'kg/m3'
  
      case(emission_scaling_flag)
        string = ''

      case(emis_factor_fire_flame_flag)
        string = 'kg/J'

      case(emis_factor_fire_smold_flag)
        string = 'kg/J'

    CASE default
        string = ' '

    END SELECT


  END FUNCTION fu_quantity_unit


  ! ***************************************************************

  LOGICAL FUNCTION fu_vector_quantity(quantity)

    ! Description:
    ! Returns ture value if the given quantity is of vector-type.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: quantity


    SELECT CASE (quantity)

      CASE (wind_flag, wind_10m_flag, mean_wind_flag)
      fu_vector_quantity = .true.

    CASE default
      fu_vector_quantity = .false.

    END SELECT

  END FUNCTION fu_vector_quantity


  ! ***************************************************************
  LOGICAL FUNCTION fu_realtime_quantity(quantity)
    !
    ! Returns true value if the given quantity has finite validity time and
    ! should be put to single_time meteo or dispersion stacks. 
    !
    IMPLICIT NONE
    INTEGER, INTENT(in) :: quantity


    SELECT CASE (quantity)

    CASE (large_scale_rain_int_flag, convective_rain_int_flag, &
            &  total_precipitation_rate_flag, scavenging_coefficient_flag, &
            &  disp_flux_cellt_rt_flag, disp_flux_celle_rt_flag, disp_flux_celln_rt_flag, &
            &  ground_pressure_flag, &
!            &  day_temperature_acc_flag, &
            & day_mean_temperature_flag, &
!            &  day_temperature_2m_acc_flag, &
            & day_mean_temperature_2m_flag)
      fu_realtime_quantity = .true.

    CASE default
      fu_realtime_quantity = .false.

    END SELECT

  END FUNCTION fu_realtime_quantity

  ! ***************************************************************

  LOGICAL FUNCTION fu_wind_quantity(quantity)

    ! Description:
    ! Returns true value if the given quantity is wind or a wind
    ! -component.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: quantity


    SELECT CASE (quantity)

      CASE (u_flag,&
          & v_flag,&
          & omega_flag,&
          & w_alt_msl_flag, w_height_srf_flag,&
          & wind_flag,&
          & u_10m_flag,&
          & v_10m_flag,&
          & wind_10m_flag,&
          & u_mean_flag,&
          & v_mean_flag,&
          & mean_wind_flag,&
          & dispersion_u_flag,&
          & dispersion_v_flag,&
          & dispersion_w_flag, &
          & eta_dot_flag )
      fu_wind_quantity = .true.

    CASE default
      fu_wind_quantity = .false.

    END SELECT

  END FUNCTION fu_wind_quantity



  ! ***************************************************************

  LOGICAL FUNCTION fu_known_quantity(quantity)

    ! Description:
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: quantity

    fu_known_quantity = ((quantity >= 240000).and.(quantity <=269999))

  END FUNCTION fu_known_quantity


  ! ***************************************************************

  LOGICAL FUNCTION fu_SILAM_dispersion_quantity(quantity)
    !
    ! Splits SILAM internal meteorological and dispersion quantities
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: quantity

    select case(quantity)
    case(particle_counter_flag, &
         & emission_intensity_flag, &
         & emission_flux_flag, &
         & areas_of_risk_flag, &
         & concentration_flag, &
         & mass_in_air_flag, &
         & advection_moment_X_flag, advection_moment_Y_flag, advection_moment_Z_flag, &
         & drydep_flag, &
         & wetdep_flag, &
         & heatsum_flag, chillsum_flag, start_heatsum_threshold_flag, &
         & daily_temp_threshold_flag, temperature_threshold_flag, soil_moisture_threshold_flag, &
         & start_calday_threshold_flag, end_calday_threshold_flag, &
         & growth_season_start_day_flag, heatsum_cutoff_tempr_flag, end_heatsum_threshold_flag, &
         & calday_start_end_diff_flag, heatsum_start_end_diff_flag, &
         & plant_growth_flag, pollen_left_relative_flag, pollen_total_per_m2_flag, &
         & pollen_correction_flag, pollen_rdy_to_fly_flag, &
         & pollen_potency_flag, allergen_rdy_to_fly_flag, &
         & day_temperature_acc_flag, day_temperature_2m_acc_flag, &
         & day_mean_temperature_flag, day_mean_temperature_2m_flag, &
         & physiography_field_set_flag, &
         & Vd_correction_DMAT_flag, &
         & interp_met2disp_coef_flag, &
         & interp_met2out_coef_flag, &
         & interp_disp2out_coef_flag, &
         & concentration_2m_flag, &
         & optical_density_flag, &
         & optical_column_depth_flag, &
         & absorption_coef_flag, scattering_coef_flag, back_scattering_coef_flag, &
         & volume_mixing_ratio_flag,&
         & cell_size_z_flag, &
         & cell_size_x_flag, &
         & cell_size_y_flag, &
         & dispersion_u_flag,&
         & dispersion_v_flag,&
         & dispersion_w_flag, &
         & air_density_flag, &
         & disp_cell_airmass_flag, &
         & cloud_cond_nucley_nbr_cnc_flag, ice_nucley_nbr_cnc_flag, &
         & disp_flux_cellt_rt_flag,disp_flux_celle_rt_flag,disp_flux_celln_rt_flag, &
         & disp_flux_celleast_flag, disp_flux_cellnorth_flag, disp_flux_celltop_flag, &
         & reaction_rate_flag,&
         & emission_scaling_flag, &
         & emis_factor_fire_flame_flag, emis_factor_fire_smold_flag)

        fu_SILAM_dispersion_quantity = .true.

      case default 
        fu_SILAM_dispersion_quantity = .false.
    end select

  END FUNCTION fu_SILAM_dispersion_quantity


  ! ***************************************************************

  LOGICAL FUNCTION fu_multi_level_quantity(quantity)
    ! 
    ! Returns false value, if a quantity is available on one level
    ! only, or has a meaning on one level only (like snow depth).
    ! This, however, does not mean that this multi-level quantity 
    ! can never appear at one level only.
    ! Be careful !!
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: quantity

    !---------------------------------------------------
    !
    ! 1. HIRLAM-quantities
    !    -----------------

    SELECT CASE (quantity)

      CASE (temperature_flag,&
          & day_temperature_acc_flag,&
          & perturb_pot_temperature_flag, &
          & u_flag,&
          & v_flag,&
          & u_mean_flag,&
          & v_mean_flag,&
          & omega_flag,&
          & w_alt_msl_flag, w_height_srf_flag,&
          & windspeed_flag, &
          & geopotential_flag,&
          & relative_humidity_flag,&
          & cloud_cover_flag,&
          & cloud_water_flag,&
          & cloud_cond_water_flag,&
          & cloud_ice_flag,&
          & specific_humidity_flag,&
          & layer_thickness_flag,&
          & height_flag,&
          & eq_pot_temperature_flag,&
          & potential_temperature_flag,&
          & relative_vorticity_flag,&
          & absolute_vorticity_flag,&
          & abs_vorticity_advection_flag,&
          & ipv_flag,&
          & tfp_flag,&
          & bulk_tfp_flag,&
          & pressure_flag,&
          & bulk_richardson_nbr_flag,&
          & flux_richardson_nbr_flag,&
          & gradient_richardson_nbr_flag,&
          & day_mean_temperature_flag,&
          & vertical_velocity_flag, &
          & Kz_momentum_3d_flag, &
          & Kz_heat_3d_flag, &
          & Kz_scalar_3d_flag, &
          & R_down_meteo_flag, &
          & turb_kinetic_energy_SILAM_flag, &
          & turb_kinetic_energy_NWP_flag, &
          & scavenging_coefficient_flag, &
          & cwcabove_3d_flag, &
          & cwcolumn_flag, &
          & lcwcabove_3d_flag, &
          & lcwcolumn_flag, &
          & pwcabove_3d_flag, &
          & pwcolumn_flag, &
          & wind_flag, &
          & humidity_mixing_ratio_flag, &
          & wind_divergence_flag, &
          & wind_vertical_shear_flag, & 
          & mean_wind_flag, &
          & turb_length_scale_flag, &
          & brunt_vaisala_freq_flag, &
          & dxdt_flag, &
          & dydt_flag, &
          & dzdt_flag, &
          & dtheta_dz_flag, &
          & dpressure_dz_flag, &
          !
          !  Dispersion quantities
          !
          & particle_counter_flag, &
          & concentration_flag, &
          & optical_density_flag, &
          & mass_in_air_flag, &
          & advection_moment_X_flag, advection_moment_Y_flag, advection_moment_Z_flag, &
          & emission_intensity_flag, &
          & emission_flux_flag, &
          & volume_mixing_ratio_flag, &
          & cell_size_z_flag,&
          & dispersion_u_flag,&
          & dispersion_v_flag,&
          & dispersion_w_flag, &
          & z_layertop_flag, p_layertop_flag, &
          & disp_flux_cellt_rt_flag, disp_flux_celle_rt_flag,disp_flux_celln_rt_flag, &
          & disp_flux_celleast_flag, disp_flux_cellnorth_flag, disp_flux_celltop_flag, &
          & disp_cell_airmass_flag, &
          & reaction_rate_flag, &
          & cloud_cond_nucley_nbr_cnc_flag, ice_nucley_nbr_cnc_flag, &
          & eta_dot_flag, &
          & air_density_flag)

        fu_multi_level_quantity = .true.

      CASE default
        fu_multi_level_quantity = .false.
      END SELECT

  END FUNCTION fu_multi_level_quantity


  ! ***************************************************************

  subroutine quantity_feasible_range(quantity, fMinAlert, fMinForce, fMaxForce, fMaxAlert)
    !
    ! Returns four values that define the range of the given quantity.
    ! fMinAlert < fMinForce <= fMaxForce < fMaxAlert
    ! Both fMin* and fMax* can be real_missing if the quantity can vary
    ! from plus to minus infinity, a rare occasion
    ! The quantity physical meaning requires it to be between
    ! fMinForce <= value <= fMaxForce.
    ! However, numerics, interpolation, etc may push it somewhat outside. This tolerance
    ! is defined via *alert limits, beyond which the situation is considered as
    ! completely impossible for normal operation.
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: quantity
    real, intent(out) :: fMinAlert, fMinForce, fMaxForce, fMaxAlert
    !
    ! Values are quantity-specific
    !
    SELECT CASE (quantity)

      CASE (temperature_flag)
      fMinAlert = 100.; fMinForce = 150.; fMaxForce = 400.; fMaxAlert = 500.

      CASE (day_temperature_acc_flag)
      fMinAlert = -10.; fMinForce = 0.; fMaxForce = 50000.; fMaxAlert = 100000.

      CASE (perturb_pot_temperature_flag)
      fMinAlert = -1000.; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = 1000.

      CASE (day_mean_temperature_flag)
      fMinAlert = 100.; fMinForce = 150.; fMaxForce = 400.; fMaxAlert = 500.

      CASE (day_mean_temperature_2m_flag)
      fMinAlert = 100.; fMinForce = 150.; fMaxForce = 400.; fMaxAlert = 500.

      CASE (potential_temperature_flag)
      fMinAlert = 100.; fMinForce = 150.; fMaxForce = 5000.; fMaxAlert = 10000.

      CASE (potential_temperature_2m_flag)
      fMinAlert = 100.; fMinForce = 150.; fMaxForce = 600.; fMaxAlert = 1000.

      CASE (pressure_flag)
      fMinAlert = -10.; fMinForce = 0.; fMaxForce = 200000.; fMaxAlert = 300000.

      CASE (u_flag)
      fMinAlert = -1000.; fMinForce = -500.; fMaxForce = 500.; fMaxAlert = 1000.

      CASE (u_mean_flag)
      fMinAlert = -1000.; fMinForce = -500.; fMaxForce = 500.; fMaxAlert = 1000.

      CASE (v_flag)
      fMinAlert = -1000.; fMinForce = -500.; fMaxForce = 500.; fMaxAlert = 1000.

      CASE (v_mean_flag)
      fMinAlert = -1000.; fMinForce = -500.; fMaxForce = 500.; fMaxAlert = 1000.

      CASE (omega_flag)
      fMinAlert = -1000.; fMinForce = -500.; fMaxForce = 500.; fMaxAlert = 1000.

      CASE (w_alt_msl_flag)
      fMinAlert = -1000.; fMinForce = -500.; fMaxForce = 500.; fMaxAlert = 1000.

      CASE (w_height_srf_flag)
      fMinAlert = -1000.; fMinForce = -500.; fMaxForce = 500.; fMaxAlert = 1000.

      CASE (wind_flag)
      fMinAlert = -1000.; fMinForce = -500.; fMaxForce = 500.; fMaxAlert = 1000.

      CASE (mean_wind_flag)
      fMinAlert = -1000.; fMinForce = -500.; fMaxForce = 500.; fMaxAlert = 1000.

      CASE (windspeed_flag)
      fMinAlert = -10.; fMinForce = 0.; fMaxForce = 500.; fMaxAlert = 1000.

      CASE (geopotential_flag)
      fMinAlert = -1.e4; fMinForce = -1.e4; fMaxForce = 100000.; fMaxAlert = 1000000.

      CASE (geopotential_sfc_flag)
      fMinAlert = -1.e4; fMinForce = -1.e4; fMaxForce = 100000.; fMaxAlert = 1000000.

      CASE (relief_height_flag)
      fMinAlert = -1.e3; fMinForce = -1.e3; fMaxForce = 10000.; fMaxAlert = 10000.

      CASE (longitude_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (latitude_flag)
      fMinAlert = -100.; fMinForce = -90.; fMaxForce = 90.; fMaxAlert = 100.

      CASE (relative_humidity_flag)
      fMinAlert = -1; fMinForce = 0; fMaxForce = 1; fMaxAlert = 10

      CASE (specific_humidity_flag)
      fMinAlert = -1; fMinForce = 0; fMaxForce = 1; fMaxAlert = 10

      CASE (humidity_mixing_ratio_flag)
      fMinAlert = -1; fMinForce = 0; fMaxForce = 1; fMaxAlert = 10

      CASE (wind_divergence_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (wind_vertical_shear_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (albedo_flag)
      fMinAlert = -1; fMinForce = 0.05; fMaxForce = 1; fMaxAlert = 10  !!! ECMWF has lowest limit of 0.06

      CASE (large_scale_accum_rain_flag)
      fMinAlert = -1; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (convective_accum_rain_flag)
      fMinAlert = -1; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (large_scale_rain_int_flag)
      fMinAlert = -1; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (convective_rain_int_flag)
      fMinAlert = -1; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (soil_moisture_vol_frac_nwp_flag)
      fMinAlert = -1; fMinForce = 0; fMaxForce = 1; fMaxAlert = 10

      case(soil_sand_mass_fraction_flag)
      fMinAlert = -1; fMinForce = 0; fMaxForce = 1; fMaxAlert = 10

      case(soil_clay_mass_fraction_flag)
      fMinAlert = -1; fMinForce = 0; fMaxForce = 1; fMaxAlert = 10
      
      case(alluvial_sedim_index_flag)
      fMinAlert = -1; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing
      
      CASE (fraction_of_ice_flag)
      fMinAlert = -1; fMinForce = 0; fMaxForce = 1; fMaxAlert = 10

      CASE (fraction_of_land_flag)
!      fMinAlert = -1; fMinForce = 0; fMaxForce = 1; fMaxAlert = 10
      fMinAlert = -1; fMinForce = 0; fMaxForce = 1000; fMaxAlert = 1000   ! pollen can show very funny numbers: it includes productivity scalng

      CASE (fraction_of_water_flag)
      fMinAlert = -1; fMinForce = 0; fMaxForce = 1; fMaxAlert = 10

      CASE (surface_roughness_meteo_flag)
      fMinAlert = -1; fMinForce = 1e-10; fMaxForce = 10; fMaxAlert = 100

      CASE (surface_roughness_disp_flag)
      fMinAlert = -1; fMinForce = 1e-10; fMaxForce = 10; fMaxAlert = 100

      CASE (land_roughness_meteo_flag)
      fMinAlert = -1; fMinForce = 1e-10; fMaxForce = 10; fMaxAlert = 200

      CASE (land_roughness_disp_flag)
      fMinAlert = -1; fMinForce = 1e-10; fMaxForce = 10; fMaxAlert = 100

      CASE (water_roughness_flag)
      fMinAlert = -1; fMinForce = 1e-10; fMaxForce = 1; fMaxAlert = 10

      CASE (ground_surface_temp_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 500; fMaxAlert = 1000

      CASE (water_surface_temp_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 500; fMaxAlert = 1000

      CASE (snow_depth_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (fraction_of_forest_flag)
      fMinAlert = -1; fMinForce = 0; fMaxForce = 1; fMaxAlert = 10

      CASE (fraction_of_erodible_soil_flag)
      fMinAlert = -1; fMinForce = 0; fMaxForce = 1; fMaxAlert = 10

      CASE (underground_temperature_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 500; fMaxAlert = 1000

      CASE (climatological_albedo_flag)
      fMinAlert = -1; fMinForce = 0; fMaxForce = 1; fMaxAlert = 10

      CASE (temperature_2m_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 500; fMaxAlert = 1000

      CASE (day_temperature_2m_acc_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (specific_humidity_2m_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1; fMaxAlert = 10

      CASE (relative_humidity_2m_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1; fMaxAlert = 10

      CASE (u_10m_flag)
      fMinAlert = -1000; fMinForce = -100; fMaxForce = 100; fMaxAlert = 1000

      CASE (v_10m_flag)
      fMinAlert = -1000; fMinForce = -100; fMaxForce = 100; fMaxAlert = 1000

      CASE (windspeed_10m_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 100; fMaxAlert = 1000

      CASE (wind_10m_flag)
      fMinAlert = -1000; fMinForce = -100; fMaxForce = 100; fMaxAlert = 1000

      CASE(water_eq_snow_depth_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1000; fMaxAlert = 10000

      CASE(charnock_parameter_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (snowfall_rate_weq_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 100; fMaxAlert = 1000

      CASE (soiltype_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (land_use_type_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (ref_evapotranspiration_flag)
      fMinAlert = 0.0; fMinForce = 0.0; fMaxForce = 1.0; fMaxAlert = 2.0

      case(water_in_soil_srf_grav_flag)
      fMinAlert = 0.0; fMinForce = 0.0; fMaxForce = 1000; fMaxAlert = 1.e5

      case(water_in_soil_deep_grav_flag)
      fMinAlert = 0.0; fMinForce = 0.0; fMaxForce = 1000; fMaxAlert = 1.e5
        
      case(water_capac_soil_srf_grav_flag)
      fMinAlert = 0.0; fMinForce = 0.0; fMaxForce = 1000; fMaxAlert = 1.e5
        
      case(water_capac_soil_deep_grav_flag)
      fMinAlert = 0.0; fMinForce = 0.0; fMaxForce = 1000; fMaxAlert = 1.e5
        
      case(free_soil_water_grav_flag)
      fMinAlert = 0.0; fMinForce = 0.0; fMaxForce = 1000; fMaxAlert = 1.e5
        
      CASE (emission_mask_flag)
      fMinAlert = -1.0; fMinForce = 0.0; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (dust_emis_0_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (total_precipitation_acc_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing
      
      CASE (total_precipitation_rate_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 100; fMaxAlert = 1000

      CASE (cloud_cover_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1; fMaxAlert = 10

      CASE (cloud_water_flag, cloud_cond_water_flag, cloud_ice_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1; fMaxAlert = 10

      CASE (integr_cloud_water_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (layer_thickness_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = 2e5

      CASE(height_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = 2e5

      CASE (surf_sw_down_radiation_ac_flag)
      fMinAlert = real_missing; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (surf_lw_down_radiation_ac_flag)
      fMinAlert = real_missing; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (surf_sw_net_radiation_ac_flag)
     fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (surf_lw_net_radiation_ac_flag)
     fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (NWP_latent_heatflux_ac_flag)
     fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (NWP_sensible_heatflux_ac_flag)
     fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (surf_sw_down_radiation_flag)
     fMinAlert = -10; fMinForce = 0; fMaxForce = 1e5; fMaxAlert = 1e6

      CASE (surf_lw_down_radiation_flag)
     fMinAlert = -10; fMinForce = 0; fMaxForce = 1e5; fMaxAlert = 1e6

      CASE (surf_sw_net_radiation_flag)
     fMinAlert = -1e6; fMinForce = -1e5; fMaxForce = 1e5; fMaxAlert = 1e6

      CASE (surf_lw_net_radiation_flag)
     fMinAlert = -1e6; fMinForce = -1e5; fMaxForce = 1e5; fMaxAlert = 1e6

      CASE (NWP_latent_heatflux_flag)
     fMinAlert = -1e6; fMinForce = -1e5; fMaxForce = 1e5; fMaxAlert = 1e6

      CASE (NWP_sensible_heatflux_flag)
     fMinAlert = -1e6; fMinForce = -1e5; fMaxForce = 1e5; fMaxAlert = 1e6

      CASE(evaporation_flag)
     fMinAlert = -10; fMinForce = 0; fMaxForce = 1e5; fMaxAlert = 1e6

      CASE (eq_pot_temperature_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1000; fMaxAlert = 10000

      CASE (relative_vorticity_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (absolute_vorticity_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (abs_vorticity_advection_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (abl_height_m_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1e4; fMaxAlert = 1e5

      CASE (nwp_abl_height_m_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1e4; fMaxAlert = 1e5

      CASE (abl_top_pressure_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1.1e5; fMaxAlert = 2e5
      
      CASE (pasquill_class_flag)
      fMinAlert = 0; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = 7

      CASE (ground_pressure_flag)
      fMinAlert = 10000; fMinForce = 10000; fMaxForce = 2e5; fMaxAlert = 2e5

      CASE (msl_pressure_flag)
      fMinAlert = 10000; fMinForce = 10000; fMaxForce = 2e5; fMaxAlert = 2e5

      CASE (surface_pressure_flag)
      fMinAlert = 10000; fMinForce = 10000; fMaxForce = 2e5; fMaxAlert = 2e5

      CASE (ipv_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (tfp_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (bulk_richardson_nbr_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (flux_richardson_nbr_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (gradient_richardson_nbr_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (cwcabove_3d_flag, cwcolumn_flag, pwcabove_3d_flag, pwcolumn_flag, lcwcabove_3d_flag, lcwcolumn_flag)
      fMinAlert = real_missing; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing
      
      CASE (scavenging_coefficient_flag)
      fMinAlert = -1; fMinForce = 0; fMaxForce = 1; fMaxAlert = 10

      CASE (vertical_velocity_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (turb_length_scale_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (brunt_vaisala_freq_flag)
      fMinAlert = -1; fMinForce = -1; fMaxForce = 1; fMaxAlert = 1

      CASE (areas_of_risk_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (bulk_tfp_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (dew_point_temp_2m_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 500; fMaxAlert = 500

      CASE (low_level_cloud_cover_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1; fMaxAlert = 10

      CASE (medium_level_cloud_cover_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1; fMaxAlert = 10

      CASE (high_level_cloud_cover_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1; fMaxAlert = 10

      CASE (total_cloud_cover_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1; fMaxAlert = 10

      CASE (sub_grid_scale_snowfall_flag)
      fMinAlert = -1; fMinForce = 0; fMaxForce = 1000; fMaxAlert = 1000

      CASE (grid_scale_snowfall_flag)
      fMinAlert = -1; fMinForce = 0; fMaxForce = 1000; fMaxAlert = 1000

      CASE (precipitable_water_flag)
      fMinAlert = -1; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (top_sw_net_radiation_ac_flag)
      fMinAlert = -1e4; fMinForce = -5e3; fMaxForce = 5e3; fMaxAlert = 1e4

      CASE (top_lw_net_radiation_ac_flag)
      fMinAlert = -1e4; fMinForce = -5e3; fMaxForce = 5e3; fMaxAlert = 1e4

      CASE (u_momentum_flux_flag)
      fMinAlert = -1e4; fMinForce = -5e3; fMaxForce = 5e3; fMaxAlert = 1e4

      CASE (v_momentum_flux_flag)
      fMinAlert = -1e4; fMinForce = -5e3; fMaxForce = 5e3; fMaxAlert = 1e4

      CASE (cape_flag)
      fMinAlert = -1e4; fMinForce = 0; fMaxForce = 2e4; fMaxAlert = 2e4

      CASE (photosynth_active_rad_flag)
      fMinAlert = -1e4; fMinForce = -5e3; fMaxForce = 5e3; fMaxAlert = 1e4

      CASE (photosynth_active_rad_ac_flag)
      fMinAlert = -1; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (MO_length_inv_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing
      
      CASE (friction_velocity_flag)
      fMinAlert = -1; fMinForce = 0; fMaxForce = 1e3; fMaxAlert = 1e3

      CASE (convective_velocity_scale_flag)
      fMinAlert = -1; fMinForce = 0; fMaxForce = 1e3; fMaxAlert = 1e3

      CASE (temperature_scale_flag)
      fMinAlert = -1000; fMinForce = -100; fMaxForce = 100; fMaxAlert = 1000
      
      CASE (humidity_scale_flag)
      fMinAlert = -1e3; fMinForce = -1; fMaxForce = 1; fMaxAlert = 1e3
      
      CASE (Prandtl_nbr_flag)
      fMinAlert = -1; fMinForce = 0; fMaxForce = 1e5; fMaxAlert = 1e5
      
      CASE (SILAM_sensible_heat_flux_flag) ! Make it work with UERRA
      fMinAlert = -1e5; fMinForce = -5e3; fMaxForce = 5e3; fMaxAlert = 1e5
 
      CASE (SILAM_latent_heat_flux_flag)
      fMinAlert = -1e5; fMinForce = -5e3; fMaxForce = 5e3; fMaxAlert = 1e5
 
      CASE (dtheta_dz_flag)
      fMinAlert = -100; fMinForce = -10; fMaxForce = 10; fMaxAlert = 100
 
      CASE (dpressure_dz_flag)
      fMinAlert = -1000; fMinForce = -100; fMaxForce = 100; fMaxAlert = 1000
 
      CASE (particle_counter_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (emission_flux_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (emission_intensity_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (concentration_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      case(volume_mixing_ratio_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (mass_in_air_flag) 
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (advection_moment_X_flag) 
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (advection_moment_Y_flag) 
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (advection_moment_Z_flag) 
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (drydep_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (wetdep_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (cell_size_x_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1e9; fMaxAlert = 1e9

      CASE (cell_size_y_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1e9; fMaxAlert = 1e9

      CASE (turb_kinetic_energy_NWP_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1e9; fMaxAlert = 1e9

      CASE (turb_kinetic_energy_SILAM_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1e9; fMaxAlert = 1e9

      CASE (Kz_momentum_3d_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1e5; fMaxAlert = 1e5

      CASE (Kz_heat_3d_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1e5; fMaxAlert = 1e5

      CASE (Kz_scalar_3d_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1e5; fMaxAlert = 1e5

      CASE (R_down_meteo_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE (Kz_scalar_1m_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1e5; fMaxAlert = 1e5

      CASE (windspeed_1lyr_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1e3; fMaxAlert = 1e3

      CASE (temperature_1lyr_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1e3; fMaxAlert = 1e3

      CASE (log_ground_pressure_flag)
      fMinAlert = 10; fMinForce = 10; fMaxForce = 20; fMaxAlert = 20

      case(silam_const_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      case (r_a_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing

      case (r_b_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing

      case (r_s_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing

      case (canopy_height_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1e3; fMaxAlert = 1e3

      case (stomatal_conductance_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing

      case(ISBA_temperature)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1e3; fMaxAlert = 1e3
      case(ISBA_u_wind)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1e3; fMaxAlert = 1e3
      case(ISBA_v_wind)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1e3; fMaxAlert = 1e3
      case(ISBA_spec_humidity)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1; fMaxAlert = 1
      case(ISBA_water_eq_of_snow)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1e5; fMaxAlert = 1e5
      case(ISBA_land_fraction)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1; fMaxAlert = 10
      case(ISBA_land_type)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing
      case(ISBA_moisture)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing
      case(ISBA_latent_hflux)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing
      case(ISBA_sensible_hflux)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing
      case(ISBA_roughness)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 100; fMaxAlert = 1e3

      case(heatsum_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1e10; fMaxAlert = 1e10
      case(chillsum_flag)
      fMinAlert = -1e10; fMinForce = -1e10; fMaxForce = 1e10; fMaxAlert = 1e10  ! can be negative

      case(pollen_rdy_to_fly_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing
      case(allergen_rdy_to_fly_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing

      case(start_calday_threshold_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 400; fMaxAlert = 500
      case(end_calday_threshold_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 400; fMaxAlert = 500
      case(calday_start_end_diff_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 400; fMaxAlert = 500

      case(start_heatsum_threshold_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing
      case(temperature_threshold_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 400; fMaxAlert = 500
      case(daily_temp_threshold_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 400; fMaxAlert = 500
      case(soil_moisture_threshold_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 400; fMaxAlert = 500
      
      case(growth_season_start_day_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 400; fMaxAlert = 400
      case(heatsum_cutoff_tempr_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 400; fMaxAlert = 400
      case(end_heatsum_threshold_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing
      case(heatsum_start_end_diff_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing

      case(pollen_correction_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing

      case(pollen_left_relative_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing
      
      case(pollen_total_per_m2_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing
      
      case(plant_growth_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing
      
      case(pollen_potency_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = real_missing; fMaxAlert = real_missing

      case( timezone_index_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      case(physiography_field_set_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      case(Vd_correction_DMAT_flag)
      fMinAlert = -1; fMinForce = 0; fMaxForce = 100; fMaxAlert = 100

      case(water_salinity_flag)
      fMinAlert = -1; fMinForce = 0; fMaxForce = 1; fMaxAlert = 10

      case(interp_met2disp_coef_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      case(interp_met2out_coef_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      case(interp_disp2out_coef_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      case(concentration_2m_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      case(optical_density_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1e5; fMaxAlert = 1e5

      case(optical_column_depth_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1e5; fMaxAlert = 1e5

      case(absorption_coef_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      case(scattering_coef_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      case(back_scattering_coef_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      case(dispersion_u_flag)
      fMinAlert = -1000.; fMinForce = -500.; fMaxForce = 500.; fMaxAlert = 1000.

      case(dispersion_v_flag)
      fMinAlert = -1000.; fMinForce = -500.; fMaxForce = 500.; fMaxAlert = 1000.
        
      case(dispersion_w_flag)
      fMinAlert = -1000.; fMinForce = -500.; fMaxForce = 500.; fMaxAlert = 1000.
      
      case(cell_size_z_flag)
      fMinAlert = -10.; fMinForce = 0.; fMaxForce = 1e10; fMaxAlert = 1e10

      case(z_layertop_flag)
      fMinAlert = -10.; fMinForce = 0.; fMaxForce = 1e10; fMaxAlert = 1e10

      case(p_layertop_flag)
      fMinAlert = -10.; fMinForce = 0.; fMaxForce = 1e10; fMaxAlert = 1e10

      case(disp_flux_cellt_rt_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing
      case(disp_flux_celle_rt_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing
      case(disp_flux_celln_rt_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing
      case(disp_flux_celleast_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing
      case(disp_flux_cellnorth_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing
      case(disp_flux_celltop_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing
      case(disp_cell_airmass_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      case(reaction_rate_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      case(cloud_cond_nucley_nbr_cnc_flag)
      fMinAlert = -1.0; fMinForce = 0.0; fMaxForce = real_missing; fMaxAlert = real_missing

      case(ice_nucley_nbr_cnc_flag)
      fMinAlert = -1; fMinForce = 0.0; fMaxForce = real_missing; fMaxAlert = real_missing

      case(leaf_area_index_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 100; fMaxAlert = 1e3

      case(leaf_area_indexhv_flag, leaf_area_indexlv_flag) !! routinely gets value of 9999 in MEPS
      fMinAlert = -10; fMinForce = 0; fMaxForce = 10; fMaxAlert = 1e4

      case(fraction_hv_flag, fraction_lv_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1.; fMaxAlert = 10.

      case(eta_dot_flag)
      fMinAlert = -1e3; fMinForce = -1e3; fMaxForce = 1e3; fMaxAlert = 1e3

      case(air_density_flag)
      fMinAlert = -10; fMinForce = 0; fMaxForce = 1e3; fMaxAlert = 1e3

      case(emission_scaling_flag)
      fMinAlert = real_missing; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      case(emis_factor_fire_flame_flag)
      fMinAlert = 0.0; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      case(emis_factor_fire_smold_flag)
      fMinAlert = 0.0; fMinForce = real_missing; fMaxForce = real_missing; fMaxAlert = real_missing

      CASE default
        call set_error('Unknown quantity:'+fu_quantity_string(quantity),'quantity_feasible_range')

    END SELECT
    
  END subroutine quantity_feasible_range

  
  !****************************************************************


  subroutine check_quantity_range(quantity, grid_data, grid_data_out, nPoints, nx, &
                                & ifRequireValidity, ifSilent, &
                                & n_out_of_range, npatched, nFailed)
    !
    ! Checks the 1D array for correspondence to the feasible range
    !
    implicit none
    
    ! Imported parameters
    integer, intent(in) :: quantity, nPoints, nx
    real, dimension(:), intent(inout) :: grid_data
    real, dimension(:), pointer ::  grid_data_out
    logical, intent(in) :: ifRequireValidity, ifSilent
    integer, intent(out) :: n_out_of_range, npatched, nFailed
    
    ! Internal variables
    integer :: ip, nPatch, ip1, iTmp, jTmp
    real :: fTmp, fMinAlert, fMinForce, fMaxForce, fMaxAlert, fMin, fMax
    logical :: ifMinAlert, ifMinForce, ifMaxForce, ifMaxAlert, ifCopy
    integer, dimension(:), pointer :: iPatch !Index to patch
    integer, dimension(4) :: ineighbour
    
    ! Get the range
    call quantity_feasible_range(quantity, fMinAlert, fMinForce, fMaxForce, fMaxAlert)
    if(error)return
    !
    ! Now careful. Any of the limits can be real_missing, then no checking happens
    ! Also, some elements of the array can be real_missing, those are skipped as well.
    ! Finally, note that this sub is going to be called plenty of times
    !
    ifMinAlert = fMinAlert /= real_missing
    ifMinForce = fMinForce /= real_missing
    ifMaxForce = fMaxForce /= real_missing
    ifMaxAlert = fMaxAlert /= real_missing
    ifCopy = associated(grid_data_out)
    n_out_of_range = 0
    npatched = 0
    nFailed = 0
    fMin = fMaxAlert
    fMax = fMinAlert
    iPatch => fu_work_int_array(nPoints)
    iPatch(1:nPoints) = 0
    !
    ! Scan the array, force min and max, and mark cells to patch.
    !
    do ip = 1, nPoints
      fTmp = grid_data(ip)
      if(.not. ifRequireValidity)then
        if(fTmp == real_missing)then
          if (ifCopy) grid_data_out(ip) = real_missing
          cycle  ! real_missing in the array
        endif
      endif     ! .not. ifRequireValidity
      
      if(ifMinAlert .and. fTmp < fMinAlert)then     ! much too low
        nPatched=npatched+1
        iPatch(ip) = -1
        fMin = min(fMin,fTmp)
#ifdef DEBUG            
        if(nPatched < 5) call make_warning('Will try to patch too low')
#endif        
        cycle
      endif
      if(ifMaxAlert .and. fTmp > fMaxAlert)then     ! much too high
        nPatched=npatched+1
        iPatch(ip) = 1
        fMax = max(fMax,fTmp)
#ifdef DEBUG            
        if(nPatched < 5) call make_warning('Will try to patch too high')
#endif        
        cycle
      endif
      if(ifMinForce .and. fTmp < fMinForce)then     ! somewhat too low
#ifdef DEBUG            
        if(n_out_of_range < 5 .and. fTmp < 0.999 * fMinForce) call make_warning('smaller-than-normal')
#endif        
        fMin = min(fMin,fTmp)
        fTmp =  fMinForce
        n_out_of_range = n_out_of_range + 1
      elseif(ifMaxForce .and. fTmp > fMaxForce)then ! somewhat too high
#ifdef DEBUG            
        if(n_out_of_range < 5 .and. fTmp > 1.001 * fMaxForce) call make_warning('bigger-than-normal')
#endif        
        fMax = max(fMax,fTmp)
        fTmp = fMaxForce
        n_out_of_range = n_out_of_range + 1
      endif
      if (ifCopy) then
        grid_data_out(ip) = fTmp
      else
        grid_data(ip) = fTmp
      endif
    end do  ! ip

    if (nPatched > 0) then
      do ip = 1, nPoints
        if(iPatch(ip) == 0) cycle
        fTmp=0
        iTmp=0
        ineighbour(:) = (/ip-nx, ip+nx, ip-1, ip+1/)
        do jTmp = 1,4 
          ip1 = ineighbour(jTmp)
          if (ip1 < 1 .or. ip1 > nPoints) cycle !Out of array
          if (iPatch(ip1) /= 0) cycle ! useless point
          fTmp = fTmp + grid_data(ip1)!
          iTmp = iTmp + 1
        enddo ! neighbours
        if (iTmp /= 0) then !Patch it
          if (ifCopy) then
            grid_data_out(ip) = fTmp / iTmp
          else
            grid_data(ip) = fTmp / iTmp
          endif
        else
          nFailed = nFailed + 1
          if(.not. ifSilent)then
            !Too bad!!
            call msg('Quantity:' + fu_quantity_string(quantity) + ', value and limits:', &
                            & (/grid_data(ip),fMinAlert,fMinForce, fMaxForce, fMaxAlert/))
            call msg('1D index, nearby values:'+fu_str(ip),grid_data(max(1,ip-5):min(nPoints,ip+5)))
            call set_error("Failed to patch the field", "check_quantity_range'")
            exit
          endif
        endif !Can patch
      end do  ! ip
      call msg('Quantity:' + fu_quantity_string(quantity) + ':' + &
                & "Total, out-of-range, patched:", (/nPoints, n_out_of_range, nPatched/))
    endif
 
    if(n_out_of_range > 0 )then
#ifdef VOIMA_GNU_BUG
!$OMP CRITICAL
#endif    
        call msg(trim('Quantity:' + fu_quantity_string(quantity) + ', out of range ')//' '&
                        &//trim(fu_str(n_out_of_range))//' times of '//trim(fu_str(nPoints))//' field min-max and limits:', &
               & (/fMin,fMax,fMinAlert,fMinForce, fMaxForce, fMaxAlert/))
#ifdef VOIMA_GNU_BUG
!$OMP END CRITICAL
#endif    
    endif
    call free_work_array(iPatch)
    
    contains
    
      subroutine make_warning(chMsg)
        character(len=*), intent(in) :: chMsg
        call msg('WARNING: Value is out of feasible range:' + chMsg + ', check_quantity_range')
        call msg('Quantity:' + fu_quantity_string(quantity) + ', value and limits:', &
               & (/fTmp,fMinAlert,fMinForce, fMaxForce, fMaxAlert/))
        call msg('1D index, nearby values:'+fu_str(ip),grid_data(max(1,ip-5):min(nPoints,ip+5)))
      end subroutine make_warning
  end subroutine check_quantity_range
  


  ! ***************************************************************

  LOGICAL FUNCTION fu_accumulated_quantity(quantity)

    ! Description:
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: quantity

    SELECT CASE (quantity)

      CASE (large_scale_accum_rain_flag,&
          & convective_accum_rain_flag,&
          & total_precipitation_acc_flag, & 
          & surf_sw_down_radiation_ac_flag,&
          & surf_lw_down_radiation_ac_flag,&
          & surf_sw_net_radiation_ac_flag,&
          & surf_lw_net_radiation_ac_flag,&
          & NWP_latent_heatflux_ac_flag,&
          & NWP_sensible_heatflux_ac_flag,&
          & sub_grid_scale_snowfall_flag,&
          & grid_scale_snowfall_flag,&
          & u_momentum_flux_flag,&
          & v_momentum_flux_flag,&
          & top_sw_net_radiation_ac_flag,&
          & top_lw_net_radiation_ac_flag, &
          & photosynth_active_rad_ac_flag, &
          & day_temperature_acc_flag, &
          & day_temperature_2m_acc_flag, &
          & drydep_flag, &
          & wetdep_flag, &
          & heatsum_flag, chillsum_flag)
    
      fu_accumulated_quantity = .true.

    CASE default

      fu_accumulated_quantity = .false.

    END SELECT

  END FUNCTION fu_accumulated_quantity



  ! ***************************************************************

  INTEGER FUNCTION fu_aver_quantity_for_acc(quantity)
    !
    ! Returns the "instantaneous" quantity flag corresponding to the 
    ! provided accumulated one. 
    !
    ! Code owner: Mikhail Sofiev

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: quantity

    IF(.not.fu_accumulated_quantity(quantity))THEN
      CAll set_error('Not an accumulated quantity','fu_aver_quantity_for_acc')
      fu_aver_quantity_for_acc = int_missing
      RETURN
    END IF

    SELECT CASE (quantity)

      CASE (large_scale_accum_rain_flag)
        fu_aver_quantity_for_acc = large_scale_rain_int_flag

      CASE (convective_accum_rain_flag)
        fu_aver_quantity_for_acc = convective_rain_int_flag

      case (total_precipitation_acc_flag)
        fu_aver_quantity_for_acc = total_precipitation_rate_flag

      CASE (surf_sw_down_radiation_ac_flag)
        fu_aver_quantity_for_acc = surf_sw_down_radiation_flag

      CASE (surf_lw_down_radiation_ac_flag)
        fu_aver_quantity_for_acc = surf_lw_down_radiation_flag

      CASE (surf_sw_net_radiation_ac_flag)
        fu_aver_quantity_for_acc = surf_sw_net_radiation_flag

      CASE (top_sw_net_radiation_ac_flag)
        fu_aver_quantity_for_acc = top_sw_net_radiation_flag

      CASE (surf_lw_net_radiation_ac_flag)
        fu_aver_quantity_for_acc = surf_lw_net_radiation_flag

      CASE (top_lw_net_radiation_ac_flag)
        fu_aver_quantity_for_acc = top_lw_net_radiation_flag

      CASE (NWP_latent_heatflux_ac_flag)
        fu_aver_quantity_for_acc = NWP_latent_heatflux_flag

      CASE (NWP_sensible_heatflux_ac_flag)
        fu_aver_quantity_for_acc = NWP_sensible_heatflux_flag

      case(photosynth_active_rad_ac_flag)
        fu_aver_quantity_for_acc = photosynth_active_rad_ac_flag

      case(day_temperature_acc_flag)
        fu_aver_quantity_for_acc = temperature_flag

      case(day_temperature_2m_acc_flag)
        fu_aver_quantity_for_acc = temperature_2m_flag

    CASE default
      CALL set_error('no match available for:' + fu_quantity_string(quantity),&
                   & 'fu_average_quantity')
      fu_aver_quantity_for_acc = int_missing

! Quantities without instantaneous analogies are: 
!
!                     sub_grid_scale_snowfall_flag,
!                     grid_scale_snowfall_flag,
!                     u_momentum_flux_flag,
!                     v_momentum_flux_flag,
!                     nuclear_cocktail_dep_flag

    END SELECT

  END FUNCTION fu_aver_quantity_for_acc


  ! ***************************************************************

  LOGICAL FUNCTION fu_quantity_in_quantities(quantity, quantities)
    
    ! Description:
    ! Returns true value if the given quantity is on the quantitylist, or
    ! would be otherwise accepted.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    INTEGER, DIMENSION(:), INTENT(in) :: quantities
    INTEGER, INTENT(in) :: quantity

    ! Local declarations:
    INTEGER :: i

    IF (quantities(1) == accept_all_quantities) THEN
      fu_quantity_in_quantities = .true.
      RETURN
    END IF

    do i=1,size(quantities)
      if(quantities(i) == int_missing)exit
      if(quantities(i) == quantity) then 
        fu_quantity_in_quantities = .true.
        return
      endif
    end do
    fu_quantity_in_quantities = .false.

  END FUNCTION fu_quantity_in_quantities


  ! ***************************************************************

  REAL FUNCTION fu_grib_missing_real(quantity)
    !
    ! Returns the value, which represents the missing data in output GRIB
    ! or GrADS files. 
    ! ATTENTION !!! 
    ! It is not necessarily really missing data. E.g. for particle_count 
    ! quantity it returns zero, while it has a clear physical meaning - 
    ! no particles. But there is another meaning for GrADS - leave the 
    ! grid cell empty, non-coloured. And there is another meaning for GRIB
    ! - mark this grid cell as no-particles in the bitmap and do not encode it.
    ! This very CONDITIONAL GRIB/GrADS MISSING is returned by this function.
    !
    ! Author: Mikhail Sofiev, FMI
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: quantity

    SELECT CASE (quantity)

    !------------------------------------------------------------------------------
    ! Quantities, for which 0.0 means "no event" and thus the grid cell should
    ! be left uncoloured
    !
    CASE ( large_scale_accum_rain_flag, &
         & convective_accum_rain_flag, &
         & large_scale_rain_int_flag, &
         & convective_rain_int_flag, &
         & soil_moisture_vol_frac_nwp_flag, &
         & fraction_of_ice_flag, &
         & fraction_of_land_flag, &
         & fraction_of_water_flag, &
         & snow_depth_flag, &
         & fraction_of_forest_flag, &
         & fraction_of_erodible_soil_flag, &
         & soil_sand_mass_fraction_flag, &
         & soil_clay_mass_fraction_flag, &
         & alluvial_sedim_index_flag, &
         & water_eq_snow_depth_flag, &
         & charnock_parameter_flag, &
         & snowfall_rate_weq_flag, &
         & total_precipitation_acc_flag, &
         & total_precipitation_rate_flag, &
         & low_level_cloud_cover_flag, &
         & medium_level_cloud_cover_flag, &
         & high_level_cloud_cover_flag, &
         & total_cloud_cover_flag, &
         & cloud_cover_flag, &
         & sub_grid_scale_snowfall_flag, &
         & grid_scale_snowfall_flag, &
         & precipitable_water_flag, &
         & particle_counter_flag, &
         & emission_intensity_flag, &
         & emission_flux_flag, &
         & concentration_flag, &
         & volume_mixing_ratio_flag, &
         & mass_in_air_flag, &
         & drydep_flag, &
         & wetdep_flag, &
         & areas_of_risk_flag, &
         & heatsum_flag, &

         & pollen_rdy_to_fly_flag, &
         & allergen_rdy_to_fly_flag, &

         & emission_mask_flag, &
         & concentration_2m_flag, &
         & optical_density_flag, &
         & optical_column_depth_flag, &
         & absorption_coef_flag, scattering_coef_flag, back_scattering_coef_flag)

      fu_grib_missing_real = 0.0

    CASE DEFAULT
      fu_grib_missing_real = real_missing
    END SELECT
    
  END FUNCTION fu_grib_missing_real


  !***************************************************************************
  
  real function fu_real_missing_replacement(quantity)
    !
    ! Returns the value, which has a funny meaning of "empty" value, e.g. zero for cloud cover
    ! In some cases, these values are replaced with grib-real missing in the input fields.
    ! Then they have to be replaced with this value.
    ! ATTENTION !!! 
    ! The function must be used only in connection with the fu_if_replace_real_missing from grib_api module.
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: quantity

    SELECT CASE (quantity)

    !------------------------------------------------------------------------------
    ! Quantities, for which 0.0 means "no event" and thus the grid cell should
    ! be left uncoloured: the clear case for stupid replacements of "empty"  with "missing"
    !
    CASE ( large_scale_accum_rain_flag, &
         & convective_accum_rain_flag, &
         & large_scale_rain_int_flag, &
         & convective_rain_int_flag, &
         & soil_moisture_vol_frac_nwp_flag, &
         & fraction_of_ice_flag, &
         & fraction_of_land_flag, &
         & fraction_of_water_flag, &
         & soil_sand_mass_fraction_flag, &
         & soil_clay_mass_fraction_flag, &
         & alluvial_sedim_index_flag, &
         & snow_depth_flag, &
         & fraction_of_forest_flag, &
         & fraction_of_erodible_soil_flag, &
         & emission_mask_flag, &
         & water_eq_snow_depth_flag, &
         & charnock_parameter_flag, &
         & snowfall_rate_weq_flag, &
         & total_precipitation_acc_flag, &
         & total_precipitation_rate_flag, &
         & low_level_cloud_cover_flag, &
         & medium_level_cloud_cover_flag, &
         & high_level_cloud_cover_flag, &
         & total_cloud_cover_flag, &
         & cloud_cover_flag, &
         & sub_grid_scale_snowfall_flag, &
         & grid_scale_snowfall_flag, &
         & precipitable_water_flag, &
         & relative_humidity_flag, &
         & specific_humidity_flag, &
         & particle_counter_flag, &
         & emission_intensity_flag, &
         & emission_flux_flag, &
         & concentration_flag, &
         & volume_mixing_ratio_flag, &
         & mass_in_air_flag, &
         & drydep_flag, &
         & canopy_height_flag, &
         & wetdep_flag, &
         & areas_of_risk_flag, &
         & heatsum_flag, &
         & pollen_rdy_to_fly_flag, &
         & allergen_rdy_to_fly_flag, &
         & concentration_2m_flag, &
         & optical_density_flag, &
         & optical_column_depth_flag, &
         & absorption_coef_flag, scattering_coef_flag, back_scattering_coef_flag, &
         & leaf_area_indexhv_flag, leaf_area_indexlv_flag, leaf_area_index_flag)

      fu_real_missing_replacement = 0.0

    CASE DEFAULT
      fu_real_missing_replacement = real_missing
    END SELECT
  
  end function fu_real_missing_replacement

  !***************************************************************************
  logical function fu_if_StaggerX(quantity)
    !
    ! Returns true for quantities that need X staggered grid
    implicit none
    
    ! Imported parameter
    integer, intent(in) :: quantity

    SELECT CASE (quantity)
      CASE (disp_flux_celleast_flag, disp_flux_celle_rt_flag)
          fu_if_StaggerX = .true.
      CASE DEFAULT
         fu_if_StaggerX = .false.
    END SELECT

   end function fu_if_StaggerX

  !***************************************************************************
  logical function fu_if_StaggerY(quantity)
    !
    ! Returns true for quantities that need X staggered grid
    implicit none
    
    ! Imported parameter
    integer, intent(in) :: quantity

    SELECT CASE (quantity)
     CASE (disp_flux_cellnorth_flag, disp_flux_celln_rt_flag)
         fu_if_StaggerY = .true.
    CASE DEFAULT
         fu_if_StaggerY = .false.
    END SELECT

   end function fu_if_StaggerY
  !***************************************************************************

  integer function fu_regridding_method(quantity)
    !
    ! Returns the default regridding method for the given quantity
    ! Essentially, we have three options:
    ! - interpolation of any kind, linear will be taken by default
    ! - averaging if some fraction or flux per m2
    ! - summing-up for total fluxes
    !
    implicit none
    
    ! Imported parameter
    integer, intent(in) :: quantity
    
    select case(quantity)
      case( surf_sw_down_radiation_ac_flag, &
         & surf_lw_down_radiation_ac_flag, &
         & surf_sw_net_radiation_ac_flag, &
         & surf_lw_net_radiation_ac_flag, &
         & NWP_latent_heatflux_ac_flag, &
         & NWP_sensible_heatflux_ac_flag, &
         & surf_sw_down_radiation_flag, &
         & surf_lw_down_radiation_flag, &
         & surf_sw_net_radiation_flag, &
         & surf_lw_net_radiation_flag, &
         & NWP_latent_heatflux_flag, &
         & NWP_sensible_heatflux_flag, &
         & evaporation_flag, &
         & sub_grid_scale_snowfall_flag, &
         & grid_scale_snowfall_flag, &
         & top_sw_net_radiation_ac_flag, &
         & top_lw_net_radiation_ac_flag, &
         & photosynth_active_rad_flag, &
         & photosynth_active_rad_ac_flag, &
         & SILAM_sensible_heat_flux_flag, &
         & SILAM_latent_heat_flux_flag, &
         & dust_emis_0_flag, &  !This should be conservative, but causes issues
                        ! when projecting coarse dust_emis to fine dispersion grids
         & leaf_area_indexhv_flag, &
         & leaf_area_indexlv_flag, &
         & concentration_flag, &        ! per m3, in theory should be averaged
         & drydep_flag, &
         & wetdep_flag)                 ! fluxes per m2 or fractions. In theory, should be averaged
         
         fu_regridding_method = linear  ! for the time being. Not too big error due to comparatively 
                                        ! small dynamic range of the fields but much faster

         case( large_scale_accum_rain_flag, &
         & convective_accum_rain_flag, &
         & large_scale_rain_int_flag, &
         & convective_rain_int_flag, &
         & snowfall_rate_weq_flag, &
         & soil_sand_mass_fraction_flag, &
         & soil_clay_mass_fraction_flag, &
         & alluvial_sedim_index_flag, &   ! Has to be averaged: high variability, wide dynamic range
         & leaf_area_index_flag, &
         & fraction_hv_flag, &
         & fraction_lv_flag, &
         & fraction_of_ice_flag, &
         & fraction_of_land_flag, &
         & fraction_of_water_flag, &
         & fraction_of_forest_flag, &
         & fraction_of_erodible_soil_flag, &
         & total_precipitation_acc_flag, &
         & total_precipitation_rate_flag, &
         & cwcabove_3d_flag, cwcolumn_flag, &
         & lcwcabove_3d_flag, lcwcolumn_flag, &
         & pwcabove_3d_flag, pwcolumn_flag, &
         & land_roughness_disp_flag, &
         & emission_flux_flag, &
         & emission_mask_flag, &
         & soil_moisture_vol_frac_nwp_flag)

        fu_regridding_method = average
         
      case(emission_intensity_flag, &
         & mass_in_air_flag)

        fu_regridding_method = summation  ! this is total-cell flux. Must be summed up

      case(timezone_index_flag, pollen_left_relative_flag) !! Sic!  

        fu_regridding_method = nearest_point 

      case(pasquill_class_flag, &
         & soiltype_flag, &
         & land_use_type_flag)
        call msg('Quantity:' + fu_quantity_string(quantity) + '- is not regriddable')

        fu_regridding_method = int_missing

      case default
        fu_regridding_method = linear

    end select
    
  end function fu_regridding_method


  !***************************************************************************

  integer function fu_get_silam_quantity(chQuantity) result(iQ)
  !
  ! Finds the integer SILAM quantity flag from its character name
  !
  implicit none

  ! Imported parameter
  character(len=*), intent(in) :: chQuantity

    if(trim(chQuantity) == "temperature")then
      iQ = temperature_flag

    elseif(trim(chQuantity) == "temperature_accum_day")then
      iQ = day_temperature_acc_flag

    elseif(trim(chQuantity) == "perturb_pot_temperature")then
      iQ = perturb_pot_temperature_flag 

    elseif(trim(chQuantity) == "temperature_2m")then
      iQ = temperature_2m_flag 

    elseif(trim(chQuantity) == "temperature_2m_acc")then
      iQ = day_temperature_2m_acc_flag 

    elseif(trim(chQuantity) == "potential_temperature")then
      iQ = potential_temperature_flag 

    elseif(trim(chQuantity) == "potential_temperature_2m")then
      iQ = potential_temperature_2m_flag 

    elseif(trim(chQuantity) == "eq_pot_temperature")then
      iQ = eq_pot_temperature_flag 

    elseif(trim(chQuantity) == "daily_mean_temperature")then
      iQ = day_mean_temperature_flag 

    elseif(trim(chQuantity) == "daily_mean_temperature_2m")then
      iQ = day_mean_temperature_2m_flag 

    elseif(trim(chQuantity) == "underground_temperature")then
      iQ = underground_temperature_flag 

    elseif(trim(chQuantity) == "pressure")then
      iQ = pressure_flag 

    elseif(trim(chQuantity) == "u")then
      iQ = u_flag 

    elseif(trim(chQuantity) == "v")then
      iQ = v_flag 

    elseif(trim(chQuantity) == "omega")then
      iQ = omega_flag

    elseif(trim(chQuantity) == "w_msl")then
      iQ = w_alt_msl_flag

    elseif(trim(chQuantity) == "w_srf")then
      iQ = w_height_srf_flag

    elseif(trim(chQuantity) == "wind")then
      iQ = wind_flag

    elseif(trim(chQuantity) == "windspeed")then
      iQ = windspeed_flag 

    elseif(trim(chQuantity) == "geopotential")then
      iQ = geopotential_flag 

    elseif(trim(chQuantity) == "geopotential_surface")then
      iQ = geopotential_sfc_flag 

    elseif(trim(chQuantity) == "relief_height")then
      iQ = relief_height_flag

    elseif(trim(chQuantity) == "longitude")then
      iQ = longitude_flag

    elseif(trim(chQuantity) == "latitude")then
      iQ = latitude_flag

    elseif(trim(chQuantity) == "relative_humidity")then
      iQ = relative_humidity_flag 

    elseif(trim(chQuantity) == "specific_humidity")then
      iQ = specific_humidity_flag 

    elseif(trim(chQuantity) == "humidity_mixing_ratio")then
      iQ = humidity_mixing_ratio_flag 

    elseif(trim(chQuantity) == "wind_divergence")then
      iQ = wind_divergence_flag 

    elseif(trim(chQuantity) == "wind_vertical_shear")then
      iQ = wind_vertical_shear_flag 

    elseif(trim(chQuantity) == "layer_thickness")then
      iQ = layer_thickness_flag 

    elseif(trim(chQuantity) == "cloud_cover")then
      iQ = cloud_cover_flag 

    elseif(trim(chQuantity) == "cloud_water")then
      iQ = cloud_water_flag 

    elseif(trim(chQuantity) == "cloud_cond_water")then
      iQ = cloud_cond_water_flag 

    elseif(trim(chQuantity) == "integr_cloud_water")then
      iQ = integr_cloud_water_flag 

    elseif(trim(chQuantity) == "cloud_ice")then
      iQ = cloud_ice_flag 

    elseif(trim(chQuantity) == "height")then
      iQ = height_flag 

    elseif(trim(chQuantity) == "relative_vorticity")then
      iQ = relative_vorticity_flag 

    elseif(trim(chQuantity) == "absolute_vorticity")then
      iQ = absolute_vorticity_flag 

    elseif(trim(chQuantity) == "abs_vorticity_advection")then
      iQ = abs_vorticity_advection_flag 

    elseif(trim(chQuantity) == "u_mean")then
      iQ = u_mean_flag 

    elseif(trim(chQuantity) == "v_mean")then
      iQ = v_mean_flag 

    elseif(trim(chQuantity) == "mean_wind")then
      iQ = mean_wind_flag 

    elseif(trim(chQuantity) == "ipv")then
      iQ = ipv_flag 

    elseif(trim(chQuantity) == "tfp")then
      iQ = tfp_flag 

    elseif(trim(chQuantity) == "bulk_richardson_nbr")then
      iQ = bulk_richardson_nbr_flag 

    elseif(trim(chQuantity) == "flux_richardson_nbr")then
      iQ = flux_richardson_nbr_flag 

    elseif(trim(chQuantity) == "gradient_richardson_nbr")then
      iQ = gradient_richardson_nbr_flag 

    elseif(trim(chQuantity) == "vertical_velocity")then
      iQ = vertical_velocity_flag

    elseif(trim(chQuantity) == "bulk_tfp")then
      iQ = bulk_tfp_flag 

    elseif(trim(chQuantity) == "dxdt")then
      iQ = dxdt_flag 

    elseif(trim(chQuantity) == "dydt")then
      iQ = dydt_flag 

    elseif(trim(chQuantity) == "dzdt")then
      iQ = dzdt_flag 

    elseif(trim(chQuantity) == "turb_kinetic_energy_NWP")then
      iQ = turb_kinetic_energy_NWP_flag 

    elseif(trim(chQuantity) == "turb_kinetic_energy_SILAM")then
      iQ = turb_kinetic_energy_SILAM_flag 

    elseif(trim(chQuantity) == "turb_length_scale")then
      iQ = turb_length_scale_flag 

    elseif(trim(chQuantity) == "brunt_vaisala_freq")then
      iQ = brunt_vaisala_freq_flag 

    elseif(trim(chQuantity) == "Kz_momentum_3d")then
      iQ = Kz_momentum_3d_flag

    elseif(trim(chQuantity) == "Kz_scalar_3d")then
      iQ = Kz_scalar_3d_flag

    elseif(trim(chQuantity) == "R_down_met")then
      iQ = R_down_meteo_flag

    elseif(trim(chQuantity) == "Kz_heat_3d")then
      iQ = Kz_heat_3d_flag 

    elseif(trim(chQuantity) == "cwcabove_3d")then
      iQ = cwcabove_3d_flag 

    elseif(trim(chQuantity) == "cwcolumn")then
      iQ = cwcolumn_flag

    elseif(trim(chQuantity) == "liquid_cwcabove_3d")then
      iQ = lcwcabove_3d_flag

    elseif(trim(chQuantity) == "liquid_cwcolumn")then
      iQ = lcwcolumn_flag

    elseif(trim(chQuantity) == "pwcabove_3d")then
      iQ = pwcabove_3d_flag 

    elseif(trim(chQuantity) == "pwcolumn")then
      iQ = pwcolumn_flag

    elseif(trim(chQuantity) == "scavenging_coefficient")then
      iQ = scavenging_coefficient_flag 

    elseif(trim(chQuantity) == "Kz_scalar_1m")then
      iQ = Kz_scalar_1m_flag 

    elseif(trim(chQuantity) == "albedo")then
      iQ = albedo_flag 

    elseif(trim(chQuantity) == "large_scale_accum_rain")then
      iQ = large_scale_accum_rain_flag 

    elseif(trim(chQuantity) == "convective_accum_rain")then
      iQ = convective_accum_rain_flag

    elseif(trim(chQuantity) == "large_scale_rain_int")then
      iQ = large_scale_rain_int_flag

    elseif(trim(chQuantity) == "convective_rain_int")then
      iQ = convective_rain_int_flag 

    elseif(trim(chQuantity) == "soil_moisture_content")then
      iQ = soil_moisture_vol_frac_nwp_flag 

    elseif(trim(chQuantity) == "soil_moisture")then !!! FIXME DUP to make test field working
      iQ = soil_moisture_vol_frac_nwp_flag 

    elseif(trim(chQuantity) == "soil_sand_mass_fraction")then
      iQ = soil_sand_mass_fraction_flag

    elseif(trim(chQuantity) == "soil_clay_mass_fraction")then
      iQ = soil_clay_mass_fraction_flag

    elseif(trim(chQuantity) == "alluv_sedim_ind")then
      iQ = alluvial_sedim_index_flag

    elseif(trim(chQuantity) == "fraction_of_ice")then
      iQ = fraction_of_ice_flag 

    elseif(trim(chQuantity) == "fraction_of_land")then
      iQ = fraction_of_land_flag 

    elseif(trim(chQuantity) == "fraction_of_water")then
      iQ = fraction_of_water_flag 

    elseif(trim(chQuantity) == "land_roughness_met")then
      iQ = land_roughness_meteo_flag 

    elseif(trim(chQuantity) == "land_roughness_disp")then
      iQ = land_roughness_disp_flag 

    elseif(trim(chQuantity) == "water_roughness")then
      iQ = water_roughness_flag 

    elseif(trim(chQuantity) == "surface_roughness_met")then
      iQ = surface_roughness_meteo_flag 

    elseif(trim(chQuantity) == "surface_roughness_disp")then
      iQ = surface_roughness_disp_flag 

    elseif(trim(chQuantity) == "ground_surface_temp")then
      iQ = ground_surface_temp_flag

    elseif(trim(chQuantity) == "water_surface_temp")then
      iQ = water_surface_temp_flag 

    elseif(trim(chQuantity) == "snow_depth")then
      iQ = snow_depth_flag 

    elseif(trim(chQuantity) == "fraction_of_forest")then
      iQ = fraction_of_forest_flag 

    elseif(trim(chQuantity) == "fraction_of_erodible_soil")then
      iQ = fraction_of_erodible_soil_flag

    elseif(trim(chQuantity) == "climatological_albedo")then
      iQ = climatological_albedo_flag 

    elseif(trim(chQuantity) == "specific_humidity_2m")then
      iQ = specific_humidity_2m_flag 

    elseif(trim(chQuantity) == "relative_humidity_2m")then
      iQ = relative_humidity_2m_flag 

    elseif(trim(chQuantity) == "u_10m" .or. trim(chQuantity) == "V_wind_10m")then
      iQ = u_10m_flag 

    elseif(trim(chQuantity) == "v_10m" .or. trim(chQuantity) == "V_wind_10m")then
      iQ = v_10m_flag 

    elseif(trim(chQuantity) == "windspeed_10m")then
      iQ = windspeed_10m_flag 

    elseif(trim(chQuantity) == "wind_10m")then
      iQ = wind_10m_flag 

    elseif(trim(chQuantity) == "snowfall_rate_weq")then
      iQ = snowfall_rate_weq_flag 

    elseif(trim(chQuantity) == "soiltype")then
      iQ = soiltype_flag 

    elseif(trim(chQuantity) == "landuse")then
      iQ = land_use_type_flag

    elseif(trim(chQuantity) == "reference_evapotranspiration")then
      iQ = ref_evapotranspiration_flag

    elseif(trim(chQuantity) == 'soil_water_surf')then
      iQ = water_in_soil_srf_grav_flag
        
    elseif(trim(chQuantity) == 'soil_water_deep')then
      iQ = water_in_soil_deep_grav_flag
        
    elseif(trim(chQuantity) == 'soil_water_cap_srf')then
      iQ = water_capac_soil_srf_grav_flag
        
    elseif(trim(chQuantity) == 'soil_water_cap_deep')then
      iQ = water_capac_soil_deep_grav_flag
        
    elseif(trim(chQuantity) == 'free_soil_water')then
      iQ = free_soil_water_grav_flag
        
    elseif(trim(chQuantity) == "emis_mask")then
      iQ = emission_mask_flag

    elseif(trim(chQuantity) == "dust_emis_0")then
      iQ = dust_emis_0_flag

    elseif(trim(chQuantity) == "total_precipitation")then
      iQ = total_precipitation_acc_flag 

    elseif(trim(chQuantity) == "total_precipitation_rate")then
      iQ = total_precipitation_rate_flag 

    elseif(trim(chQuantity) == "surf_sw_down_radiation_ac")then
      iQ = surf_sw_down_radiation_ac_flag 

    elseif(trim(chQuantity) == "surf_lw_down_radiation_ac")then
      iQ = surf_lw_down_radiation_ac_flag 

    elseif(trim(chQuantity) == "surf_sw_net_radiation_ac")then
      iQ = surf_sw_net_radiation_ac_flag

    elseif(trim(chQuantity) == "surf_lw_net_radiation_ac")then
      iQ = surf_lw_net_radiation_ac_flag

    elseif(trim(chQuantity) == "nwp_latent_heatflux_ac")then
      iQ = NWP_latent_heatflux_ac_flag

    elseif(trim(chQuantity) == "nwp_sensible_heatflux_ac")then
      iQ = NWP_sensible_heatflux_ac_flag

    elseif(trim(chQuantity) == "surf_sw_down_radiation")then
      iQ = surf_sw_down_radiation_flag 

    elseif(trim(chQuantity) == "surf_lw_down_radiation")then
      iQ = surf_lw_down_radiation_flag 

    elseif(trim(chQuantity) == "surf_sw_net_radiation")then
      iQ = surf_sw_net_radiation_flag

    elseif(trim(chQuantity) == "surf_lw_net_radiation")then
      iQ = surf_lw_net_radiation_flag

    elseif(trim(chQuantity) == "nwp_latent_heatflux")then
      iQ = NWP_latent_heatflux_flag

    elseif(trim(chQuantity) == "nwp_sensible_heatflux")then
      iQ = NWP_sensible_heatflux_flag

    elseif(trim(chQuantity) == "evaporation")then
      iQ = evaporation_flag

    elseif(trim(chQuantity) == "abl_height_m")then
      iQ = abl_height_m_flag 

    elseif(trim(chQuantity) == "abl_top_pressure")then
      iQ = abl_top_pressure_flag 

    elseif(trim(chQuantity) == "pasquill_class")then
      iQ = pasquill_class_flag 

    elseif(trim(chQuantity) == "ground_pressure")then
      iQ = ground_pressure_flag 

    elseif(trim(chQuantity) == "msl_pressure")then
      iQ = msl_pressure_flag 

    elseif(trim(chQuantity) == "surface_pressure")then
      iQ = surface_pressure_flag 

    elseif(trim(chQuantity) == "dew_point_temp_2m")then
      iQ = dew_point_temp_2m_flag 

    elseif(trim(chQuantity) == "low_level_cloud_cover")then
      iQ = low_level_cloud_cover_flag 

    elseif(trim(chQuantity) == "medium_level_cloud_cover")then
      iQ = medium_level_cloud_cover_flag

    elseif(trim(chQuantity) == "high_level_cloud_cover")then
      iQ = high_level_cloud_cover_flag 

    elseif(trim(chQuantity) == "total_cloud_cover")then
      iQ = total_cloud_cover_flag 

    elseif(trim(chQuantity) == "sub_grid_scale_snowfall")then
      iQ = sub_grid_scale_snowfall_flag 

    elseif(trim(chQuantity) == "grid_scale_snowfall")then
      iQ = grid_scale_snowfall_flag

    elseif(trim(chQuantity) == "precipitable_water")then
      iQ = precipitable_water_flag 

    elseif(trim(chQuantity) == "top_sw_net_radiation_ac")then
      iQ = top_sw_net_radiation_ac_flag

    elseif(trim(chQuantity) == "top_lw_net_radiation_ac")then
      iQ = top_lw_net_radiation_ac_flag

    elseif(trim(chQuantity) == "u_momentum_flux")then
      iQ = u_momentum_flux_flag 

    elseif(trim(chQuantity) == "v_momentum_flux")then
      iQ = v_momentum_flux_flag 

    elseif(trim(chQuantity) == "cape")then
      iQ = cape_flag

    elseif(trim(chQuantity) == "photosynth_active_rad")then
      iQ = photosynth_active_rad_flag 

    elseif(trim(chQuantity) == "photosynth_active_rad_ac")then
      iQ = photosynth_active_rad_ac_flag 

    elseif(trim(chQuantity) == "top_sw_net_radiation")then
      iQ = top_sw_net_radiation_flag 

    elseif(trim(chQuantity) == "top_lw_net_radiation")then
      iQ = top_lw_net_radiation_flag 

    elseif(trim(chQuantity) == "nwp_abl_height_m")then
      iQ = nwp_abl_height_m_flag 

    elseif(trim(chQuantity) == "silam_sensible_heat_flux")then
      iQ = SILAM_sensible_heat_flux_flag

    elseif(trim(chQuantity) == "silam_latent_heat_flux")then
      iQ = SILAM_latent_heat_flux_flag

    elseif(trim(chQuantity) == "cell_size_x")then
      iQ = cell_size_x_flag 

    elseif(trim(chQuantity) == "cell_size_y")then
      iQ = cell_size_y_flag 

    elseif(trim(chQuantity) == "weq_snow_depth")then
      iQ = water_eq_snow_depth_flag 

    elseif(trim(chQuantity) == "charnock")then
      iQ = charnock_parameter_flag 

    elseif(trim(chQuantity) == "MO_length_inv")then
      iQ = MO_length_inv_flag

    elseif(trim(chQuantity) == "friction_velocity")then
      iQ = friction_velocity_flag 

    elseif(trim(chQuantity) == "convective_velocity_scale")then
      iQ = convective_velocity_scale_flag 

    elseif(trim(chQuantity) == "cnv_vel_scale")then
      iQ = convective_velocity_scale_flag 

    elseif(trim(chQuantity) == "temperature_scale")then
      iQ = temperature_scale_flag 

    elseif(trim(chQuantity) == "humidity_scale")then
      iQ = humidity_scale_flag 

    elseif(trim(chQuantity) == "Prandtl_nbr")then
      iQ = Prandtl_nbr_flag 

    elseif(trim(chQuantity) == "windspeed_1lyr")then
      iQ = windspeed_1lyr_flag 

    elseif(trim(chQuantity) == "temperature_1lyr")then
      iQ = temperature_1lyr_flag 

    elseif(trim(chQuantity) == "ISBA_temperature")then
      iQ = ISBA_temperature 

    elseif(trim(chQuantity) == "ISBA_u_wind")then
      iQ = ISBA_u_wind 

    elseif(trim(chQuantity) == "ISBA_v_wind")then
      iQ = ISBA_v_wind 

    elseif(trim(chQuantity) == "ISBA_spec_humidity")then
      iQ = ISBA_spec_humidity 

    elseif(trim(chQuantity) == "ISBA_water_eq_of_snow")then
      iQ = ISBA_water_eq_of_snow 

    elseif(trim(chQuantity) == "ISBA_land_type")then
      iQ = ISBA_land_type 

    elseif(trim(chQuantity) == "ISBA_land_fraction")then
      iQ = ISBA_land_fraction 

    elseif(trim(chQuantity) == "ISBA_moisture")then
      iQ = ISBA_moisture 

    elseif(trim(chQuantity) == "ISBA_latent_hflux")then
      iQ = ISBA_latent_hflux

    elseif(trim(chQuantity) == "ISBA_sensible_hflux")then
      iQ = ISBA_sensible_hflux 

    elseif(trim(chQuantity) == "ISBA_roughness")then
      iQ = ISBA_roughness

    elseif(trim(chQuantity) == "log_ground_pressure")then
      iQ = log_ground_pressure_flag 

    elseif(trim(chQuantity) == "particle_counter")then
      iQ = particle_counter_flag 

    elseif(trim(chQuantity) == "emission_rate")then
      iQ = emission_intensity_flag 

    elseif(trim(chQuantity) == "emission_flux")then
      iQ = emission_flux_flag 

    elseif(trim(chQuantity) == "concentration")then
      iQ = concentration_flag 

    elseif(trim(chQuantity) == "volume_mixing_ratio")then
      iQ = volume_mixing_ratio_flag

    elseif(trim(chQuantity) == "mass_in_air")then
      iQ = mass_in_air_flag 

    elseif(trim(chQuantity) == "advection_moment_x")then
      iQ = advection_moment_X_flag 

    elseif(trim(chQuantity) == "advection_moment_y")then
      iQ = advection_moment_Y_flag 

    elseif(trim(chQuantity) == "advection_moment_z")then
      iQ = advection_moment_Z_flag 

    elseif(trim(chQuantity) == "drydep")then
      iQ = drydep_flag 

    elseif(trim(chQuantity) == "wetdep")then
      iQ = wetdep_flag 

    elseif(trim(chQuantity) == "areas_of_risk")then
      iQ = areas_of_risk_flag 

    elseif(trim(chQuantity) == "r_a_resistance")then
      iQ = r_a_flag 

    elseif(trim(chQuantity) == "r_b_resistance")then
      iQ = r_b_flag 

    elseif(trim(chQuantity) == "r_s_resistance")then
      iQ = r_s_flag 

    elseif(trim(chQuantity) == "h_canopy")then
      iQ = canopy_height_flag

    elseif(trim(chQuantity) == "g_stomatal")then
      iQ = stomatal_conductance_flag 

    elseif(trim(chQuantity) == "dtheta_dz")then
      iQ = dtheta_dz_flag 

    elseif(trim(chQuantity) == "dpressure_dz")then
      iQ = dpressure_dz_flag

    elseif(trim(chQuantity) == "heatsum")then
      iQ = heatsum_flag

    elseif(trim(chQuantity) == "chillsum")then
      iQ = chillsum_flag

    elseif(trim(chQuantity) == "pollen_rdy_to_fly")then
      iQ = pollen_rdy_to_fly_flag

    elseif(trim(chQuantity) == "allergen_rdy_to_fly")then
      iQ = allergen_rdy_to_fly_flag

    elseif(trim(chQuantity) == "start_calday_threshold")then
      iQ = start_calday_threshold_flag

    elseif(trim(chQuantity) == "end_calday_threshold")then
      iQ = end_calday_threshold_flag

    elseif(trim(chQuantity) == "calendar_day_start_end_diff")then
      iQ = calday_start_end_diff_flag

    elseif(trim(chQuantity) == "start_heatsum_threshold")then
      iQ = start_heatsum_threshold_flag

    elseif(trim(chQuantity) == "temperature_threshold")then
      iQ = temperature_threshold_flag

    elseif(trim(chQuantity) == "daily_temp_threshold")then
      iQ = daily_temp_threshold_flag
   
    elseif(trim(chQuantity) == "soil_moisture_threshold")then
      iQ = soil_moisture_threshold_flag

    elseif(trim(chQuantity) == "growth_season_start_day")then
      iQ = growth_season_start_day_flag

    elseif(trim(chQuantity) == "heatsum_cutoff_tempr")then
      iQ = heatsum_cutoff_tempr_flag

    elseif(trim(chQuantity) == "end_heatsum_threshold")then
      iQ = end_heatsum_threshold_flag

    elseif(trim(chQuantity) == "heatsum_start_end_diff")then
      iQ = heatsum_start_end_diff_flag

    elseif(trim(chQuantity) == "pollen_correction")then
      iQ = pollen_correction_flag

    elseif(trim(chQuantity) == "pollen_left")then
      iQ = pollen_left_relative_flag

    elseif(trim(chQuantity) == "pollen_total_m2")then
      iQ = pollen_total_per_m2_flag
    
    elseif(trim(chQuantity) == "plant_growth")then
      iQ = plant_growth_flag
    
    elseif(trim(chQuantity) == "pollen_potency")then
      iQ = pollen_potency_flag

    elseif(trim(chQuantity) == "tz_index")then
      iQ = timezone_index_flag

    elseif(trim(chQuantity) == "physiography_field_set")then
      iQ = physiography_field_set_flag

    elseif(trim(chQuantity) == "Vd_correction_DMAT")then
      iQ = Vd_correction_DMAT_flag

    elseif(trim(chQuantity) == "water_salinity")then
      iQ = water_salinity_flag

    elseif(trim(chQuantity) == "interp_met2disp_coef_flag")then
      iQ = interp_met2disp_coef_flag

    elseif(trim(chQuantity) == "interp_met2out_coef_flag")then
      iQ = interp_met2out_coef_flag

    elseif(trim(chQuantity) == "interp_disp2out_coef_flag")then
      iQ = interp_disp2out_coef_flag

    elseif(trim(chQuantity) == "concentration_2m")then
      iQ = concentration_2m_flag

    elseif(trim(chQuantity) == "optical_density")then
      iQ = optical_density_flag

    elseif(trim(chQuantity) == "optical_column_depth")then
      iQ = optical_column_depth_flag

    elseif(trim(chQuantity) == "absorp_coef")then
      iQ = absorption_coef_flag

    elseif(trim(chQuantity) == "scatt_coef")then
      iQ = scattering_coef_flag

    elseif(trim(chQuantity) == "backscatt_coef")then
      iQ = back_scattering_coef_flag

    elseif(trim(chQuantity) == "dispersion_u")then
      iQ = dispersion_u_flag

    elseif(trim(chQuantity) == "dispersion_v")then
      iQ = dispersion_v_flag

    elseif(trim(chQuantity) == "dispersion_w")then
      iQ = dispersion_w_flag

    elseif(trim(chQuantity) == "lai")then
      iQ = leaf_area_index_flag

    elseif(trim(chQuantity) == "lai_hv")then
      iQ = leaf_area_indexhv_flag

    elseif(trim(chQuantity) == "lai_lv")then
      iQ = leaf_area_indexlv_flag

    elseif(trim(chQuantity) == "frac_hv")then
      iQ = fraction_hv_flag

    elseif(trim(chQuantity) == "frac_lv")then
      iQ = fraction_lv_flag

    elseif(trim(chQuantity) == "eta_dot") then
      iQ = eta_dot_flag

    elseif(trim(chQuantity) == "disp_ztop") then
       iQ= z_layertop_flag

    elseif(trim(chQuantity) == "disp_ptop") then
       iQ= p_layertop_flag

    elseif(trim(chQuantity) == "d_flux_cellt_rt") then
       iQ= disp_flux_cellt_rt_flag

    elseif(trim(chQuantity) == "d_flux_celle_rt") then
       iQ= disp_flux_celle_rt_flag

    elseif(trim(chQuantity) == "d_flux_celln_rt") then
       iQ= disp_flux_celln_rt_flag

    elseif(trim(chQuantity) == "d_flux_celln") then
       iQ= disp_flux_cellnorth_flag

    elseif(trim(chQuantity) == "d_flux_celle") then
       iQ= disp_flux_celleast_flag

    elseif(trim(chQuantity) == "d_flux_cellt") then
       iQ= disp_flux_celltop_flag

    elseif(trim(chQuantity) == "d_cellmass") then
       iQ= disp_cell_airmass_flag

    elseif(trim(chQuantity) == "react_rate") then
       iQ= reaction_rate_flag

    elseif(trim(chQuantity) == "ccn") then
       iQ= cloud_cond_nucley_nbr_cnc_flag
       
    elseif(trim(chQuantity) == "in") then
       iQ= ice_nucley_nbr_cnc_flag

    elseif(trim(chQuantity) == "cellsize_z")then
      iQ = cell_size_z_flag

    elseif(trim(chQuantity) == 'air_dens') then
      iQ = air_density_flag

    elseif(trim(chQuantity) == 'ems_scale') then
      iQ = emission_scaling_flag

    elseif(trim(chQuantity) == 'emis_factor_fire_flame') then
      iQ = emis_factor_fire_flame_flag

    elseif(trim(chQuantity) == 'emis_factor_fire_smoulder') then
      iQ = emis_factor_fire_smold_flag

    else
      call set_error('Unknown quantity name:' + chQuantity, 'fu_get_silam_quantity')
    endif
  end function fu_get_silam_quantity

  !*************************************************************************


    FUNCTION fu_quantity_standard_name(quantity) result(string)

    ! Description:
    ! Returns the name of quantity in a string.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    CHARACTER(LEN=45) :: string

    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: quantity

    SELECT CASE (quantity)

      CASE (temperature_flag)
        string =  "air_temperature" 

      CASE (day_temperature_acc_flag)
        string =  "air_temperature_accum_daily" 

      CASE (perturb_pot_temperature_flag) 
        string =   "perturb_pot_temp" 

      CASE (pressure_flag) 
        string =   "air_pressure" 

      CASE (u_flag) 
        string =   "eastward_wind" 

      CASE (v_flag) 
        string =   "northward_wind" 

      CASE (omega_flag)
        string =   "lagrangian_tendency_of_air_pressure" 

      CASE (windspeed_flag) 
        string =   "wind_speed" 

      CASE (geopotential_flag) 
        string =   "geopotential" 

      CASE (geopotential_sfc_flag) 
        string =   "geopotential_surface" 

      CASE (longitude_flag)
        string =   "longitude" 

      CASE (latitude_flag)
        string =   "latitude" 

      CASE (relative_humidity_flag) 
        string =   "relative_humidity" 

      CASE (specific_humidity_flag) 
        string =   "specific_humidity" 

      CASE (humidity_mixing_ratio_flag) 
        string =   "humidity_mixing_ratio" 

      CASE (wind_divergence_flag) 
        string =   "divergence_of_wind" 

      CASE (potential_temperature_flag) 
        string =   "air_potential_temperature" 

      CASE (cloud_cover_flag) 
        string =   "cloud_area_fraction" 

      CASE (cloud_water_flag) 
        string =   "mass_fraction_of_cloud_liquid_water_in_air" 

      CASE (cloud_cond_water_flag) 
        string =   "mass_fraction_of_cloud_condensed_water_in_air" 

      CASE (integr_cloud_water_flag) 
        string =   "atmosphere_cloud_liquid_water_content"  

      CASE (cloud_ice_flag) 
        string =   "mass_fraction_of_cloud_ice_water_in_air" 

      CASE (height_flag) 
        string =   "height" 

      CASE (eq_pot_temperature_flag) 
        string =   "equivalent_potential_temperature" 

      CASE (relative_vorticity_flag) 
        string =   "atmosphere_relative_vorticity" 

      CASE (absolute_vorticity_flag) 
        string =   "atmosphere_absolute_vorticity" 

      CASE (albedo_flag) 
        string =   "surface_albedo" 

      CASE (large_scale_accum_rain_flag) 
        string =   "large_scale_rainfall_amount" 

      CASE (convective_accum_rain_flag)
        string =   "convective_rainfall_amount" 

      CASE (large_scale_rain_int_flag)
        string =   "large_scale_rainfall_flux"  

      CASE (convective_rain_int_flag) 
        string =   "convective_rainfall_flux" 

      CASE (soil_moisture_vol_frac_nwp_flag) 
        string =   "soil_moisture_content" 

      case(soil_sand_mass_fraction_flag)
        string = "soil_sand_mass_fraction"

      case(soil_clay_mass_fraction_flag)
        string = "soil_clay_mass_fraction"

      case(alluvial_sedim_index_flag)
        string = "alluv_sedim_ind"

      CASE (fraction_of_ice_flag) 
        string =   "sea_ice_area_fraction" 

      CASE (fraction_of_land_flag) 
        string =   "land_area_fraction" 

      CASE (fraction_of_water_flag) 
        string =   "water_area_fraction" 

      CASE (surface_roughness_meteo_flag) 
        string =   "surface_roughness_length_meteo" 

      CASE (surface_roughness_disp_flag) 
        string =   "surface_roughness_length_disp"

      CASE (water_surface_temp_flag) 
        string =   "sea_surface_temperature" 

      CASE (snow_depth_flag) 
        string =   "surface_snow_thickness" 

      CASE (underground_temperature_flag) 
        string =   "soil_temperature" 

      CASE (snowfall_rate_weq_flag) 
        string =   "lwe_snowfall_rate" 

      CASE (soiltype_flag) 
        string =   "soil_type" 

      CASE (land_use_type_flag) 
        string =   "landuse" 

      CASE (ref_evapotranspiration_flag) 
        string =   "ref_evapotr" 

      case(water_in_soil_srf_grav_flag)
        string = 'soil_water_surf'
        
      case(water_in_soil_deep_grav_flag)
        string = 'soil_water_deep'
        
      case(water_capac_soil_srf_grav_flag)
        string = 'soil_water_cap_srf'

      case(water_capac_soil_deep_grav_flag)
        string = 'soil_water_cap_deep'

      case(free_soil_water_grav_flag)
        string = 'free_soil_water'

      case(emission_mask_flag)
        string =   "emis_mask"
        
      case(dust_emis_0_flag)
        string =   "dust_emis"
    
      CASE (total_precipitation_acc_flag) 
        string =   "precipitation_amount" 

      CASE (total_precipitation_rate_flag) 
        string =   "precipitation_flux" 

      CASE (surf_sw_down_radiation_ac_flag) 
        string =   "surf_sw_down_radiation_ac" 

      CASE (surf_lw_down_radiation_ac_flag) 
        string =   "surf_lw_down_radiation_ac" 

      CASE (surf_sw_down_radiation_flag) 
        string =   "downwelling_shortwave_flux_in_air" 

      CASE (surf_lw_down_radiation_flag) 
        string =   "downwelling_longwave_flux_in_air" 

      CASE (surf_sw_net_radiation_ac_flag)
        string =   "integral_of_surface_net_downward_shortwave_flux_wrt_time" 

      CASE (surf_lw_net_radiation_ac_flag)
        string =   "integral_of_surface_net_downward_longwave_flux_wrt_time" 

      CASE (surf_sw_net_radiation_flag)
        string =   "surface_net_downward_shortwave_flux" 

      CASE (surf_lw_net_radiation_flag)
        string =   "surface_net_downward_longwave_flux" 

      CASE (top_sw_net_radiation_ac_flag)
        string =   "integral_of_toa_net_downward_shortwave_flux_wrt_time" 

      CASE (top_lw_net_radiation_ac_flag)
        string =   "integral_of_toa_outgoing_longwave_flux_wrt_time" 

      CASE (top_sw_net_radiation_flag) 
        string =   "toa_net_downward_shortwave_flux" 

      CASE (top_lw_net_radiation_flag) 
        string =   "toa_net_downward_longwave_flux" 

      CASE (NWP_latent_heatflux_ac_flag)
        string =   "integral_of_surface_downward_latent_heat_flux_wrt_time" 

      CASE (NWP_sensible_heatflux_ac_flag)
        string =   "integral_of_surface_downward_sensible_heat_flux_wrt_time" 

      CASE (NWP_latent_heatflux_flag)
        string =   "surface_downward_latent_heat_flux" 

      CASE (NWP_sensible_heatflux_flag)
        string =   "surface_downward_sensible_heat_flux" 

      CASE (evaporation_flag)
        string =   "water_evaporation_amount" 

      CASE (abl_height_m_flag) 
        string =   "atmosphere_boundary_layer_thickness" 

      CASE (ground_pressure_flag) 
        string =   "surface_air_pressure" 

      CASE (msl_pressure_flag) 
        string =   "air_pressure_at_sea_level" 

      CASE (dew_point_temp_2m_flag) 
        string =   "dew_point_temperature" 

      CASE (sub_grid_scale_snowfall_flag) 
        string =   "convective_snowfall_amount" 

      CASE (grid_scale_snowfall_flag)
        string =   "large_scale_snowfall_amount" 

      CASE (u_momentum_flux_flag) 
        string =   "downward_eastward_momentum_flux_in_air" 

      CASE (v_momentum_flux_flag) 
        string =   "downward_northward_momentum_flux_in_air" 

      CASE (cape_flag)
         string = "convective_available_potential_energy"

      CASE (water_eq_snow_depth_flag) 
        string =   "lwe_thickness_of_surface_snow_amount" 

      CASE (water_salinity_flag)
        string =   "sea_water_salinity" 

      CASE (emission_scaling_flag)
        string =   "emission_scaling" 

      case(emis_factor_fire_flame_flag)
        string = 'emis_factor_fire_flame'

      case(emis_factor_fire_smold_flag)
        string = 'emis_factor_fire_smoulder'



!      CASE (emission_intensity_flag)
!    string =   "tendency_of_atmosphere_mass_content_of_<SUBST>_due_to_emission" 
!      CASE (concentration_flag) 
!    string =   "mass_concentration_of_<SUBST>_in_air" , "mole_concentration_of_<SUBST>_in_air" 
!      CASE (drydep_flag)
!    string =   "tendency_of_atmosphere_mass_content_of_<SUBST>_due_to_dry_deposition" 
!      CASE (wetdep_flag)
!    string =   "tendency_of_atmosphere_mass_content_of_<SUBST>_due_to_wet_deposition" 
!      CASE (optical_column_depth_flag)
!    string =   "atmosphere_optical_thickness_due_to_<SUBST>" 

    CASE default
      string = ''

    END SELECT

  END FUNCTION fu_quantity_standard_name

!*************************************************************************

  integer function fu_outGridInterpType(iQ)

    implicit none

    integer, intent(in) :: iQ

    select case(iQ)

    case(large_scale_accum_rain_flag, &
       & convective_accum_rain_flag, &
       & large_scale_rain_int_flag, &
       & convective_rain_int_flag, &
       & snowfall_rate_weq_flag, &
       & total_precipitation_acc_flag, &
       & total_precipitation_rate_flag, &
       & sub_grid_scale_snowfall_flag, &
       & grid_scale_snowfall_flag, &
       & precipitable_water_flag, total_cloud_cover_flag, &
       & cloud_cond_water_flag, cloud_water_flag, cloud_ice_flag, cloud_cover_flag, &
       & cwcabove_3d_flag, cwcolumn_flag, &
       & lcwcabove_3d_flag, lcwcolumn_flag, &
       & pwcabove_3d_flag, pwcolumn_flag, &
       & scavenging_coefficient_flag, &
       & low_level_cloud_cover_flag, medium_level_cloud_cover_flag, high_level_cloud_cover_flag, &
       & integr_cloud_water_flag, fraction_of_forest_flag, plant_growth_flag, &
       & pollen_left_relative_flag, pollen_total_per_m2_flag, &
       & fraction_of_erodible_soil_flag, emission_mask_flag)

       fu_outGridInterpType = setZero

    case(fraction_of_land_flag, &
       & heatsum_flag, chillsum_flag, &
       & start_calday_threshold_flag, &
       & start_heatsum_threshold_flag, &
       & daily_temp_threshold_flag, &
       & temperature_threshold_flag, &
       & soil_moisture_threshold_flag, &
       & growth_season_start_day_flag, &
       & end_calday_threshold_flag, &
       & calday_start_end_diff_flag, &
       & heatsum_cutoff_tempr_flag, &
       & end_heatsum_threshold_flag, &
       & heatsum_start_end_diff_flag, &
       & pollen_correction_flag, day_temperature_2m_acc_flag, day_mean_temperature_2m_flag, &
       & water_salinity_flag, &
       & pollen_potency_flag)

       fu_outGridInterpType = nearestPoint
       
       case(u_flag, v_flag, w_alt_msl_flag, w_height_srf_flag, omega_flag, vertical_velocity_flag, &
       & relative_vorticity_flag, absolute_vorticity_flag, abs_vorticity_advection_flag, &
       & u_mean_flag, v_mean_flag, mean_wind_flag, wind_flag, windspeed_flag, &
       & wind_divergence_flag, wind_vertical_shear_flag, &
       & u_10m_flag, v_10m_flag, windspeed_10m_flag, wind_10m_flag, &
       & cell_size_x_flag, cell_size_y_flag, longitude_flag,  latitude_flag, &
       & dispersion_u_flag, dispersion_v_flag, dispersion_w_flag, eta_dot_flag, &
       & z_layertop_flag, p_layertop_flag, &
       & disp_flux_cellt_rt_flag,disp_flux_celle_rt_flag,disp_flux_celln_rt_flag, &
       & disp_flux_celleast_flag, &
       & disp_flux_cellnorth_flag, &
       & disp_flux_celltop_flag, &
       & disp_cell_airmass_flag, &
       & reaction_rate_flag, &
       & cloud_cond_nucley_nbr_cnc_flag, ice_nucley_nbr_cnc_flag, &
       & cell_size_z_flag, emission_scaling_flag, &
       & fraction_hv_flag, fraction_lv_flag,leaf_area_indexhv_flag, leaf_area_indexlv_flag,leaf_area_index_flag, &
       & timezone_index_flag,  &
       & emis_factor_fire_flame_flag, emis_factor_fire_smold_flag)

      fu_outGridInterpType = notAllowed

    case default
      ! no fu_outGridInterpType
      fu_outGridInterpType = notAllowed
    end select


  end function  fu_outGridInterpType



!*************************************************************************
subroutine report_list_of_quantities(list, title)
        !
        ! outputs array of quantities
        !
        integer, dimension(:), intent(in) :: list
        character(len=*), optional, intent(in) :: title
        integer :: i
        
        call msg("")
        if (present(title)) call msg("List:"+title)
        do i = 1,size(list)
                if (list(i) == int_missing) exit
                call msg(fu_quantity_string(list(i)))
        enddo
end subroutine report_list_of_quantities



!*************************************************************************

subroutine reformat_names_of_quantities()

  implicit none
  
  integer :: uOut, iQ
  
  uOut = fu_next_free_unit()

  open(uOut,file='d:\!model\2011\silam_v5_0\ini\quantity_summary.txt')

  do iQ = first_multi_level_q, last_multi_level_q
    if(fu_known_quantity(iQ)) call print_name_list_for_quantity(uOut, iQ)
  enddo
  
  do iQ = first_single_level_q, last_single_level_q
    if(fu_known_quantity(iQ)) call print_name_list_for_quantity(uOut, iQ)
  enddo
  
  do iQ = first_silam_own_q, last_silam_own_q
    if(fu_known_quantity(iQ)) call print_name_list_for_quantity(uOut, iQ)
  enddo
  

CONTAINS

subroutine print_name_list_for_quantity(uOut, quantity)
 implicit none
 integer, intent(in) :: uOut, quantity
 real :: fMinAlert, fMinForce, fMaxForce, fMaxAlert
 
 write(uOut,'(2A7)')'LIST = ',fu_str(quantity)
 write(unit=uOut,fmt='(A18,I6)')'  quantity_flag = ',quantity 
 write(unit=uOut,fmt='(2A)') '  short_name = ',fu_quantity_short_string(quantity)
 write(unit=uOut,fmt='(A18,I6)')'  long_name = ',fu_quantity_string(quantity)
 write(unit=uOut,fmt='(10A)') '  unit = ',fu_quantity_unit(quantity)
 call quantity_feasible_range(quantity, fMinAlert, fMinForce, fMaxForce, fMaxAlert)
 write(unit=uOut,fmt='(A,4(1x,E15.3))')' valid_ranges = ', fMinAlert, fMinForce, fMaxForce, fMaxAlert
 write(unit=uOut,fmt='(A,E15.3)')'  default_real_missing = ',fu_grib_missing_real(quantity)
 if(fu_accumulated_quantity(quantity))then
   write(unit=uOut,fmt='(A)')'  cumulative_quantity = YES'
   write(unit=uOut,fmt='(A,I15)')'  instant_quantity = ',fu_aver_quantity_for_acc(quantity)
 else
   write(unit=uOut,fmt='(A)')'  cumulative_quantity = NO'
 endif
 write(unit=uOut,fmt='(2A)') '  netcdf_cf_standard_name = ', fu_quantity_standard_name(quantity)
 write(uOut,'(A8)')'END_LIST = ',fu_str(quantity)
 write(uOut,*)''

! <flag> 
! <string_main> 
! <string_short> 
! <valid_ranges> 
! <real_missing> 
! <ifCumulative> 
! <instant for cumulative> 
! <string_standard> 
! <unit> 
! <full string>
end subroutine print_name_list_for_quantity


end subroutine reformat_names_of_quantities



END MODULE names_of_quantities



