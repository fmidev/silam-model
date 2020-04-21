MODULE derived_field_quantities

  ! Description:
  ! This module contains tools for calculating some new fields
  ! from the nwp weather data in supermarket. The new fields are
  ! stored to the supermarket along with the original, nwp's own fields.
  ! 
  ! Author: Mikhail Sofiev, FMI 
  ! Original code: Mika Salonoja, FMI 
  !
  ! All units: SI
  ! 
  ! All degrees: real numbers in degrees (degrees and parts per
  ! hundred/parts per thousand , NOT in degrees and minutes
  ! 
  ! Language: ANSI standard Fortran 90
  ! 
  ! Modules used:

  USE derived_field_quantities_2
!  use silam_times
!  USE temperature_tools

  IMPLICIT NONE

  ! The public functions and subroutines available in this module:
!  PUBLIC quantities_for_derived !FIXME: not single-time clean
  public quantities_for_derived_one
  public check_input_quantity
  public set_deriving_way
  PUBLIC make_derived_meteo_fields
  public ifCanBeSkippedInDerivation

  ! The private functions and subroutines not to be used elsewhere:
  private dq_surface_roughness
  private dq_surface_pressure
  PRIVATE dq_windspeed
  PRIVATE dq_layer_thickness
  PRIVATE dq_average_between_obstimes
  PRIVATE dq_height_from_z
  PRIVATE dq_height_from_t
  PRIVATE dq_rh_from_spechum
  PRIVATE dq_mean_vertically_from_bottom
  private dq_ground_pressure
  private dq_total_precipitation
  private dq_total_prec_rate
  private dq_ground_pr_from_msl_pr
  private dq_windspeed_10m
  private dq_3d_pressure
  private dq_rh_from_spechum_2m
  private dq_spechum_2m_from_dewp_2m
  private dq_total_precipitation_ls
  private dq_total_prec_rate_ls
  private dq_tests


CONTAINS


!       
!  ! ***************************************************************
!WARNINRG!!!  not single-time clean
!  SUBROUTINE quantities_for_derived (quantities, q_st, wdr) !, ifPP)
!  !
!  ! Takes the list of quantites. If it contains derived ones, whose
!  ! computing requires some other ones - it adds them.
!  !
!  IMPLICIT NONE
!
!  ! Imported with intent IN
!  type(silja_wdr), intent(in) :: wdr
!
!  ! Imported with intent inout
!  INTEGER, DIMENSION(:), INTENT(inout) :: quantities, q_st
!
!  ! Local stuff
!  INTEGER iTmp, iCount
!  integer, dimension(max_Q_4_derived) :: list, list_st
!
!  iCount=1
!
!  if(all(quantities == int_missing))then
!    call set_error('No quantities given','quantities_for_derived')
!    return
!  end if
!
!  ! A trick: any added quantity may be a derived one, so it must be a 
!  ! multi-pass cycle
!
!  DO WHILE (iCount > 0) ! Counter of added quantities
!!    print *, iCount
!    DO iTmp=1,SIZE(quantities)
!!      print *, iCount, iTmp
!      IF(quantities(iTmp) == int_missing) EXIT
!      call quantities_for_derived_one(quantities(iTmp), list, list_st, wdr)
!      iCount = max(fu_merge_int_arrays(list, quantities, .true.), &
!                 & fu_merge_int_arrays(list_st, q_st, .true.))
!!      print *, iCount
!    END DO
!  END DO
!
!  END SUBROUTINE quantities_for_derived
!

  !*************************************************************************
  logical function ifCanBeSkippedInDerivation(q)
    ! True for quantities that 
    !    1. might be needed only as intermediate ones in derivation chain
    !  and
    !    2. can come on their own from input or be 
    !  and 
    !    3. whose derivatives can be obtained from their precedors directly, 
    !      so q can be bypassed
    ! 
    IMPLICIT NONE
    integer :: q

    ifCanBeSkippedInDerivation =  any(q ==(/cloud_cond_water_flag/))

  end function ifCanBeSkippedInDerivation

  SUBROUTINE quantities_for_derived_one(dq, request, ifST, list, list_st, requests, requests_st, wdr)
    !
    ! This function returns a list of quantities needed for the
    ! calculation of a given derived quantity. For example for
    ! potential temperature this function returns a list containing
    ! only temperature, and for scavenging coefficient this function
    ! returns all rains and pressures.
    !
    ! This function must be updated every time new derived parameter introduced.
    !
    IMPLICIT NONE
    !
    ! Imported parameters 
    INTEGER, DIMENSION(:), intent(out) :: list, list_st, requests, requests_st
    INTEGER, INTENT(in) :: dq, request
    logical, intent(in) :: ifST 
    type(silja_wdr), intent(in) :: wdr

    integer :: iq
    
    list = int_missing
    list_st = int_missing

!    call msg('Looking for derived quantity input:' + fu_quantity_short_string(dq))

    SELECT CASE (dq)

      case(temperature_flag)
        list(1) = perturb_pot_temperature_flag
        requests(1) = request

      case(day_mean_temperature_flag)
        list_st(1) = day_temperature_acc_flag
        requests_st(1) = request

      case(day_temperature_acc_flag)
        list(1) = temperature_flag
        requests(1) = request

      case(day_mean_temperature_2m_flag)
        list_st(1) = day_temperature_2m_acc_flag
        requests_st(1) = request

      case(day_temperature_2m_acc_flag)
        list(1) = temperature_2m_flag
        requests(1) = request

      CASE (pasquill_class_flag)
        list(1:2) = (/MO_length_inv_flag, surface_roughness_meteo_flag/)
        requests(1:2) = request

      CASE (scavenging_coefficient_flag)
        if (ifST) then ! Only single-time now
          list(1:2) = (/pressure_flag, &
                      & temperature_flag/)
          requests(1:2) = request
          if (fu_number_of_precip_flds(wdr) == 2) then
            list_st(1:3) = (/large_scale_rain_int_flag, &
                           & convective_rain_int_flag, &
                           & ground_pressure_flag/)
            requests_st(1:3) = request
          else
            list_st(1:2) = (/large_scale_rain_int_flag,&
                           & ground_pressure_flag/)
            requests_st(1:2) = request
          endif
        else
            call msg('Requested quantity: ' + fu_quantity_short_string(dq))
            call set_error("Requesting realtime quantity as a dynamic one..", &
                        "quantities_for_derived_one")
        endif

      CASE (cloud_cond_water_flag) ! ...
        list(1:2) = (/cloud_water_flag, &
                          & cloud_ice_flag/)
        requests(1:2) = request

      CASE (cwcabove_3d_flag, cwcolumn_flag, pwcabove_3d_flag, pwcolumn_flag, &
           & lcwcabove_3d_flag, lcwcolumn_flag) ! ...
        list(1:2) = (/ground_pressure_flag, cloud_cond_water_flag/)
        list_st(1) = large_scale_rain_int_flag
        requests(1:2) = request
        requests_st(1) = 1 !!! optional, so large_scale_rain_int_flag is available 
               !! for scavenging if it can be done

      case (r_a_flag)  !, r_b_flag, r_s_flag)
        list(1:3) = (/MO_length_inv_flag, surface_roughness_meteo_flag, &
                    & friction_velocity_flag /)
        requests(1:3) = request

      CASE (pressure_flag, height_flag, cell_size_z_flag)
        list(1:2)=(/temperature_flag, ground_pressure_flag/)
        requests(1:2) = request

      CASE (windspeed_flag, mean_wind_flag)
        list(1:2)=(/u_flag, v_flag/)
        requests(1:2) = request

      CASE (windspeed_10m_flag)
        list(1:2)=(/u_10m_flag, v_10m_flag/)
        requests(1:2) = request

      CASE (layer_thickness_flag)
        list(1)=(geopotential_flag)
        requests(1) = request

      case(specific_humidity_2m_flag)
        list(1:2) = (/dew_point_temp_2m_flag, ground_pressure_flag/)
        requests(1:2) = request

      CASE (relative_humidity_flag)
        list(1:3)=(/specific_humidity_flag, &
                  & temperature_flag, pressure_flag/)
        requests(1:3) = request

      CASE (relative_humidity_2m_flag)
        list(1:3)=(/specific_humidity_2m_flag, &
                  & temperature_2m_flag, ground_pressure_flag/)
        requests(1:3) = request

      CASE (eq_pot_temperature_flag)
        list(1:2)=(/specific_humidity_flag, temperature_flag/)
        requests(1:2) = request

      CASE (potential_temperature_flag, tfp_flag) !tfp - thermal front parameter
        list(1:2)=(/temperature_flag, temperature_2m_flag/)  ! 2m is not needed for 3d one, but
                        ! dq_potential_temperature will use both in any case
        requests(1:2) = request

      CASE (relative_vorticity_flag, abs_vorticity_advection_flag)
        list(1:3)=(/temperature_flag, u_flag, v_flag/)
        requests(1:3) = request

      CASE (ipv_flag) ! isentropic potential vorticity
        list(1:4)=(/potential_temperature_flag, & 
                  & u_flag, v_flag, & 
                  & temperature_flag/)
        requests(1:4) = request

      CASE (bulk_richardson_nbr_flag)
        list(1:6)=(/height_flag, &
                  & u_flag, &
                  & v_flag, &
                  & friction_velocity_flag, &
                  & specific_humidity_flag, &
                  & potential_temperature_flag/)
        requests(1:6) = request

      CASE (gradient_richardson_nbr_flag)
        list(1:5)=(/height_flag, &
                  & windspeed_flag, &
                  & windspeed_10m_flag, &
                  & temperature_2m_flag, &
                  & potential_temperature_flag/)
        requests(1:5) = request

      case(brunt_vaisala_freq_flag)
        list(1:4) = (/potential_temperature_flag, &
                    & height_flag, &
                    & temperature_2m_flag, &
                    & abl_height_m_flag/)
        requests(1:4) = request

      CASE (w_alt_msl_flag, w_height_srf_flag) ! m/s vertical wind
        list(1:2) = (/omega_flag, temperature_flag/)
        requests(1:2) = request

      CASE (omega_flag)  ! Omega-vertical wind
        list(1:4)=(/u_flag, &
                  & v_flag, &
                  & temperature_flag, &
                  & ground_pressure_flag/)
        requests(1:4) = request

      case(vertical_velocity_flag) ! vertical-dependent z-component of wind
        list(1) = omega_flag   ! at least this ...
        requests(1) = request

      case(eta_dot_flag)
        list(1:3) = (/u_flag, v_flag, ground_pressure_flag/)
        requests(1:3) = request

      CASE(wind_divergence_flag)  ! div(wind_vector)
        list(1:4)=(/u_flag, &
                  & v_flag, &
                  & w_height_srf_flag, &
                  & height_flag/)
        requests(1:4) = request

      CASE(wind_vertical_shear_flag)  ! div(wind_vector)
        list(1:2)=(/u_flag, &
                  & v_flag/)
        requests(1:2) = request

      CASE (bulk_tfp_flag)
        list(1:2)=(/layer_thickness_flag, temperature_flag/)
        requests(1:2) = request

      CASE (friction_velocity_flag, & !They will be created by dq_ABL_params
          & temperature_scale_flag, &
          & SILAM_sensible_heat_flux_flag, &
          & MO_length_inv_flag, &
          & Kz_scalar_1m_flag, &
          & Prandtl_nbr_flag, &
          & convective_velocity_scale_flag, &
          & abl_height_m_flag)

        list(1:8)=(/ ground_pressure_flag,&
                   & temperature_2m_flag,&
                   & windspeed_10m_flag,&
                   & height_flag, &
                   & temperature_flag, &
                   & potential_temperature_flag, &
                   & u_flag, v_flag/) !, &
        iQ = 9
        list_st(1) = fraction_of_land_flag
        requests_st(1) = request

        select case(fu_abl_param(wdr))
          case(abl_full_param)
            list(iQ:iQ+1)=(/specific_humidity_flag, specific_humidity_2m_flag/)
            iQ=  iQ + 2
          case(abl_dry_param)
          case default
            call set_error('Only abl_full_param and abl_dry_param so far','quantities_for_derived_one')
        end select
        select case(fu_ablh_method(wdr))
          case(constant_abl_height)
          case(parcel_method)
          case(richardson_method)
          case(combination_method)
          case(coriolis_method)
          case(nwp_abl)
            list(iQ) = nwp_abl_height_m_flag
            iQ = iQ+1
          case default
            call set_error('Unknown ABL height method:'+fu_str(fu_ablh_method(wdr)),'')
        end select

        requests(1:iQ) = request
                
      case(SILAM_latent_heat_flux_flag, humidity_scale_flag)
        list(1:7)=(/ground_pressure_flag, &
                  & temperature_2m_flag, &
                  & height_flag, &
                  & specific_humidity_flag, &
                  & specific_humidity_2m_flag, &
                  & Kz_scalar_1m_flag, &
                  & Prandtl_nbr_flag/)
        requests(1:7) = request

      case(Kz_momentum_3d_flag, &
         & Kz_heat_3d_flag, &
         & Kz_scalar_3d_flag, &
         & turb_kinetic_energy_SILAM_flag)
        list(1:5) = (/gradient_richardson_nbr_flag, &
                    & windspeed_flag, &
                    & windspeed_10m_flag, &
                    & height_flag, &
                    & abl_height_m_flag/)
        requests(1:5) = request


    case(turb_length_scale_flag)
        list(1:6) = (/ ground_pressure_flag, &
                     & height_flag, &
                     & specific_humidity_flag, &
                     & u_flag, v_flag, & 
                     & temperature_flag/)
        requests(1:6) = request

    case(R_down_meteo_flag)
      select case(fu_kz_param(wdr))
        case(zero_kz)
          list(1:2) = (/ ground_pressure_flag, &
                             & height_flag/)
          requests(1:2) = request
        case(simple_abl_ec_ft_kz, silam_abl_ec_ft_kz, ec_kz,silam_kz_emulator, hunten_kz, simple_kz)

          list(1:8) = (/ ground_pressure_flag, &
                       & Kz_scalar_1m_flag, &
                       & abl_height_m_flag, &
                       & specific_humidity_flag, &
                       & height_flag, &
                       & u_flag, v_flag, & 
                       & temperature_flag/)
          requests(1:8) = request
        case default
          call msg('Unknown kz_method: ', fu_kz_param(wdr))
          call set_error("Unknown kz_method", "quantities_for_derived_one")
      endselect

    case(leaf_area_index_flag)
      select case(fu_LAIsrc(wdr))
        case(LAI_dynamic_2)
             if (.not. ifSt) then
                list(1:2) = (/leaf_area_indexhv_flag, &
                            & leaf_area_indexlv_flag /)
                 requests(1:2) = request

                 list_st(1:2) = (/fraction_hv_flag, &
                            & fraction_lv_flag /)
                 requests_st(1:2) = request
              else
                 call set_error("Tried to request lai static from dynamic ones", &
                     &  "quantities_for_derived_one")
                 return
              endif
        case(LAI_static_2)
             if (ifSt) then
                list_st(1:4) = (/leaf_area_indexhv_flag, &
                            & leaf_area_indexlv_flag, &
                            & fraction_hv_flag, &
                            & fraction_lv_flag /)
                 requests_st(1:4) = request
              else
                 call set_error("Tried to request lai dynamic from static ones", &
                                         &  "quantities_for_derived_one")
                 return
              endif
         case (LAI_static_1, LAI_dynamic_1)
               if (ifSt .eqv. (fu_LAIsrc(wdr) == LAI_static_1) ) then
                  call msg('Assumed initial quantity: ' + fu_quantity_short_string(dq))
               else
                  if (ifsT) call msg("Static LAI requested, while dynamic was set in control file")  
                 call set_error("Tried to request wrong LAI", &
                                         &  "quantities_for_derived_one")
                 return
               endif
         case(int_missing)
             call set_error('leaf area inxed is required but data source - use_lai - is missing','')
         case default
                 call set_error("unknnown wdr%LAIsrc", &
                                         &  "quantities_for_derived_one")
      end select
    
    case(stomatal_conductance_flag)
        select case(fu_LAIsrc(wdr))
           case (LAI_dynamic_1, LAI_dynamic_2)
             list(1:6) = (/  surface_pressure_flag,&
                       & total_cloud_cover_flag, &
                        & temperature_2m_flag, &
                        & relative_humidity_2m_flag,&
                        & leaf_area_index_flag, &
                        & soil_moisture_vol_frac_nwp_flag/)
             requests(1:6) = request
           case (LAI_static_1, LAI_static_2)
             list(1:5) = (/  surface_pressure_flag,&
                       & total_cloud_cover_flag, &
                        & temperature_2m_flag, &
                        & relative_humidity_2m_flag, &
                        & soil_moisture_vol_frac_nwp_flag/)
             requests(1:5) = request
             list_st(1) = leaf_area_index_flag 
             requests_st(1) = request
           case default
              call msg("LAIsrc", fu_LAIsrc(wdr))
              call set_error("Unknown LAIsrc", "quantities_for_derived_one")
        end select

      CASE (albedo_flag)
        list(1:2)=(/surf_sw_down_radiation_flag, surf_sw_net_radiation_flag/)
        requests(1:2) = request

      CASE (surf_lw_down_radiation_flag) ! W/m2 come from 2-time fields J/m2
        list(1) = surf_lw_down_radiation_ac_flag
        requests(1) = request

      CASE (surf_sw_down_radiation_flag) ! W/m2 come from 2-time fields J/m2
        list(1) = surf_sw_down_radiation_ac_flag
        requests(1) = request
        
      CASE (surf_sw_down_radiation_ac_flag) ! W/m2 come from 2-time fields J/m2
        list(1:2)=(/surf_sw_net_radiation_ac_flag, climatological_albedo_flag/)
        requests(1:2) = request
        
      CASE (surf_sw_net_radiation_flag) ! W/m2 come from 2-time fields J/m2
        list(1) = surf_sw_net_radiation_ac_flag
        requests(1) = request

      CASE (large_scale_rain_int_flag) ! kg/m2s come from 2-time kg/m2
        if (ifST) then
          list(1) = large_scale_accum_rain_flag
          requests(1) = request
        else
          call msg('Requested quantity: ' + fu_quantity_short_string(dq))
          call set_error("Requesting realtime quantity as a dynamic one..", &
                        "quantities_for_derived_one")
        endif

      CASE(convective_rain_int_flag) ! kg/m2s come from 2-time kg/m2
        if (ifST) then
          list(1) = convective_accum_rain_flag
          requests(1) = request
        else
            call msg('Requested quantity: ' + fu_quantity_short_string(dq))
            call set_error("Requesting realtime quantity as a dynamic one..", &
                        "quantities_for_derived_one")
        endif

      CASE(NWP_sensible_heatflux_flag) ! W/m2 come from 1-time J/m2
        list(1) = NWP_sensible_heatflux_ac_flag
        requests(1) = request

      CASE(NWP_latent_heatflux_flag) !  W/m2 come from 1-time J/m2
        list(1) = NWP_latent_heatflux_ac_flag
        requests(1) = request

      case (surface_roughness_meteo_flag)             ! Surface roughness has to be combined 
        list(1) = land_roughness_meteo_flag  ! from soil and water roughnesses
        requests(1) =  1 !!!request Can force default
        !list_st(1:2) = (/land_roughness_meteo_flag, fraction_of_land_flag/)  !DOES NOT WORK!!!
        !requests_st(1:2) =(/1, request/)  
        list_st(1) = fraction_of_land_flag
        requests_st(1) = request

      case (surface_roughness_disp_flag)             ! Surface roughness has to be combined 
         !Request both then priority is: Dynamic -> static -> default
        list(1) =  land_roughness_disp_flag  ! from soil and water roughnesses 
        list_st(1:2) = (/land_roughness_disp_flag, fraction_of_land_flag/)
        requests(1) = 1
        requests_st(1:2) =(/1, request/)  

!      case (surface_pressure_flag) ! Points to one of them depending on vertical type
!        list(1:2) = (/msl_pressure_flag, ground_pressure_flag/)
!        requests(1:2) = request

!      case (land_roughness_flag) ! Soil roughness can be substituted with     
!        list(1) = silam_const_flag ! a constant value without disasterous consequences  NOT ANY MORE
!        list_st(1) = fraction_of_land_flag

      case (water_roughness_flag) 
        list(1) = friction_velocity_flag
        list_st(1) = fraction_of_land_flag
        requests(1) = request
        requests_st(1) = request

      CASE (abl_top_pressure_flag)
        list(1:2)=(/abl_height_m_flag,temperature_flag/)
        requests(1:2) = request

      case (ground_pressure_flag)
        if (ifST) then ! Single-time is diagnosed from dynamic
          list(1) = ground_pressure_flag
        else  !taken from meteo
          list(1) = log_ground_pressure_flag ! Ohh, stupid, but forced by ECMWF files
        endif
        requests(1) = request

      case (log_ground_pressure_flag)
        list(1:2) = (/msl_pressure_flag, & ! Ohh, even worse but at least allows our standard way of branching
                    & temperature_2m_flag/)
        list_st(1) = relief_height_flag
        requests(1:2) = request
        requests_st(1) = request

      case (relief_height_flag)
        list_st(1) = geopotential_sfc_flag
        requests_st(1) = request

      case (total_precipitation_rate_flag)
        if (ifST)  then
          list(1) = total_precipitation_acc_flag
          requests(1) = request
        else
          call msg('Requested quantity: ' + fu_quantity_short_string(dq))
          call set_error("Requesting realtime quantity as a dynamic one..", &
                        "quantities_for_derived_one")
        endif

      case (total_precipitation_acc_flag)
        if (fu_number_of_precip_flds(wdr) == 2) then
          list(1:2) = (/large_scale_accum_rain_flag, convective_accum_rain_flag/)
          requests(1:2) = request
        else
          list(1) = large_scale_accum_rain_flag
          requests(1) = request
        end if
      
      case (heatsum_flag)
        list(1:2) = (/temperature_flag, temperature_2m_flag/)
        requests(1:2) = request

      case(concentration_flag)
        list(1:3) = (/volume_mixing_ratio_flag, temperature_flag, pressure_flag/)
        requests(1:3) = request

      case(ref_evapotranspiration_flag)
        list(1:7) = (/ temperature_2m_flag, surf_sw_net_radiation_flag, surf_lw_net_radiation_flag, &
                     & leaf_area_index_flag, surface_pressure_flag, &
!                     & air_density_flag, &
                     & relative_humidity_2m_flag, r_a_flag/)
        requests(1:7) = request
      
      case(photosynth_active_rad_flag)
!        list(1) = 
        call set_error('photosynth_active_rad_flag does not work yet','quantities_for_derived_one')
        
      case(total_cloud_cover_flag)
        list(1) = cloud_cover_flag
        requests(1) = request

      CASE DEFAULT
        call msg('Assumed initial quantity: ' + fu_quantity_short_string(dq))

    END SELECT
  
  END SUBROUTINE quantities_for_derived_one


  ! ***************************************************************

  recursive subroutine check_input_quantity(q_derived, ifST, &     ! Q to check
                                          & q_avail_dyn, &   ! list of available Qs
                                          & q_avail_st, &   ! list of available Qs
                                          & q_shop, &        ! Out: needed Qs from q_avail
                                          & q_shop_st, & ! Out: needed static to make q_derived
                                          & ifOK, &      ! If q_derived can be made
                                          & wdr, &  
                                          & q_intermediate, & ! Store the steps of deriving
                                          & q_intermediate_st, &
                                          & depthin)
    !
    ! The subroutine checks if q_derived can by ANY mean be obtained from 
    ! the q_avail set. ifOK is then true if it is possible. q_shop is a 
    ! sub-set of q_avail needed for the given q_derived
    !
    ! If present, the intermediate quantities are also returned. They appear when
    ! more than one deriving step is involved.
    !
    ! ATTENTION. Also returns the static quantities needed to derive the requested dynamic ones
    !            No checking of their availability is made, all is left for physiography
    !
    ! This a recursive procedure - be careful !
    !
    implicit none

    ! Imported parameters with intent IN
    integer, intent(in) :: q_derived
    logical, intent(in) :: ifST
    type(silja_wdr), intent(in) :: wdr
    integer, dimension(:), intent(in) :: q_avail_dyn,  q_avail_st
    integer, intent(in), optional ::  depthin


    ! Imported parameters with intent OUT
    integer, dimension(:), intent(out) :: q_shop, q_shop_st
    integer, dimension(:), intent(out), optional :: q_intermediate, q_intermediate_st
    logical, intent(out) :: ifOK

    ! Local variables
    integer :: i, NAvail, NAvailSt,depth
    integer, dimension(max_Q_4_derived) :: QTmp, q_shopTmp, q_shopStaticTmp, q_st, &
                                         & q_intermTmp, q_intermStaticTmp, req, req_st

    character(len=*), parameter :: lead= &
            & "--------------------------------------------------------------------------------------"

    if(.not. present(depthin))then
      depth = 1
    else
      depth = depthin
    endif
    !
    ! First, some preparation / speed-up
    !
!    if (ifST) then
!       call msg(lead(1:depth)+'>Checking singletime:' + fu_quantity_short_string(q_derived))
!    else
!       call msg(lead(1:depth)+'>Checking dyamic:' + fu_quantity_short_string(q_derived))
!    endif 

    do i=1, size(q_avail_st)
      if(q_avail_st(i) == int_missing)exit
    end do
    NAvailSt = i-1
    do i=1, size(q_avail_dyn)
      if(q_avail_dyn(i) == int_missing)exit
    end do
    NAvail = i-1
    q_shop = int_missing
    q_shop_st = int_missing
    ifOK = .false.
    QTmp = int_missing
    q_shopTmp = int_missing
    q_shopStaticTmp = int_missing
    q_intermTmp = int_missing
    q_intermStaticTmp = int_missing

    !
    ! If the quantity is silam_const_flag - it will be produced by the model,
    ! so nothing should be added into the shopping list
    !
    if(q_derived == silam_const_flag) then
      ifOK = .true.
 !      call msg(lead(1:depth)+'->Assumed constant:' + fu_quantity_short_string(q_derived))
      return
    end if

    !
    ! If the quantity is in q_avail - the task is solved. 
    !
    if(ifST) then
      if (any(q_avail_st(1:NAvailSt) == q_derived) .or. &
        & (.not. fu_realtime_quantity(q_derived) ))then
         ifOK = .true.
         q_shop_st(1) = q_derived
 !               call msg(lead(1:depth)+'->Found!:' + fu_quantity_short_string(q_derived))
         return
      endif
!      call msg("Checking static:"+fu_quantity_short_string(q_derived)) 
    else !dynamic
      if ( any(q_avail_dyn(1:NAvail) == q_derived)) then
         ifOK = .true.
         q_shop(1) = q_derived
 !               call msg(lead(1:depth)+'->Found!:' + fu_quantity_short_string(q_derived))
         return
      endif
!      call msg("Checking dynamic:"+fu_quantity_short_string(q_derived)) 
    end if

    !
    ! If the quantity is surface_pressure_flag - it will be EITHER msl_pressure_flag
    ! OR ground_pressure_flag. So, if any of these two flags is in the list or can be derived 
    ! - the task is solved. They are quite equal because can be obtained from each other using
    ! temperature field and orography.
    !
    if(q_derived == surface_pressure_flag) then

      if(any(q_avail_dyn(1:NAvail) == ground_pressure_flag)) then
        ifOK = .true.
        q_shop(1) = ground_pressure_flag
        return
      else
        !
        ! Try to derive the ground_pressure_flag. Here we may have recursive 
        ! procedure. Idea: make a one step in-depth and see if it is enough
        !
        q_intermTmp(count(q_intermTmp /= int_missing)+1) = ground_pressure_flag

        call quantities_for_derived_one(ground_pressure_flag, 2, .false., QTmp, q_st, req, req_st, wdr) ! One step only
        
        if(.not.all(QTmp(:) == int_missing)) then ! Can be derived, check if it helps
          do i=1, size(QTmp)
            if(QTmp(i) == int_missing) return ! all done
            ifOk = .True.
            if (req(i) == 2) call check_input_quantity(QTmp(i), .false.,&   ! q_derived
                                    & q_avail_dyn, &  ! available dynamic quantities
                                    & q_avail_st, &
                                    & q_shopTmp, &    ! needed to make the derivation
                                    & q_shopStaticTmp, &    ! needed for the derivation
                                    & ifOK, wdr, &
                                    & q_intermTmp, q_intermStaticTmp,depth+1)
            if(error)return

            if(ifOK) then
              NAvail = fu_merge_int_arrays(q_shopTmp, q_shop, .true.) ! Prepare arrays
              NAvailSt = fu_merge_int_arrays(q_shopStaticTmp, q_shop_st, .true.)
              if(present(q_intermediate))then
                NAvail = fu_merge_int_arrays(q_intermTmp, q_intermediate, .true.)
                NAvailSt = fu_merge_int_arrays(q_intermStaticTmp, q_intermediate_st, .true.)
              endif
            else
              return  ! faliure
            endif
          end do
        end if
      end if  ! if ground_pressure exists in list
      !
      ! Alternative - if ground is not availabe, take msl and convert through relief height
      !
      if(any(q_avail_dyn(1:NAvail) == msl_pressure_flag))then
        ifOK = .true.
        q_shop(1) = msl_pressure_flag
        q_shop_st(1) = relief_height_flag
        return
      else
        !
        ! Try to derive the msl_pressure_flag. Here we may have recursive 
        ! procedure. Idea: make a one step in-depth and see if it is enough
        !
        q_intermTmp(count(q_intermTmp /= int_missing)+1) = msl_pressure_flag
        call quantities_for_derived_one(msl_pressure_flag, 2, .false., QTmp, q_st, req, req_st, wdr) ! One step only

        if(.not.all(QTmp(:) == int_missing)) then ! Can be derived, check if it helps
          do i=1, size(QTmp)
            if(QTmp(i) == int_missing) return ! all done successfully. Arrays have been made
            ifOK = .True.
            if (req(i) == 2) call check_input_quantity(QTmp(i), .false., &   ! q_derived
                                    & q_avail_dyn, &  ! available dynamic quantities
                                    & q_avail_st, &
                                    & q_shopTmp, &    ! needed to make the derivation
                                    & q_shopStaticTmp, &    ! needed for the derivation
                                    & ifOK, wdr, &
                                    & q_intermTmp, q_intermStaticTmp, depth+1)
            if(error)return
            if(ifOK) then
              NAvail = fu_merge_int_arrays(q_shopTmp, q_shop, .true.)
              NAvailSt = fu_merge_int_arrays(q_shopStaticTmp, q_shop_st, .true.)
              if(present(q_intermediate))then
                NAvail = fu_merge_int_arrays(q_intermTmp, q_intermediate, .true.)
                NAvailSt = fu_merge_int_arrays(q_intermStaticTmp, q_intermediate_st, .true.)
              endif
            else
              return  ! failure
            endif
          end do
        end if
      end if  ! if msl_pressure exists in list

    end if  ! q_derived = surface_pressure

    !
    ! Try to derive the given q_derived quantity. Here we may have recursive 
    ! procedure. Idea: make a one step in-depth and see if it is enough
    !
    q_intermTmp(count(q_intermTmp /= int_missing)+1) = q_derived
    call quantities_for_derived_one(q_derived, 2, ifST, QTmp, q_st, req, req_st, wdr) ! One step only

    if((QTmp(1) == int_missing) .and. (q_st(1) == int_missing)) then ! No more deriving options.
        call msg(fu_connect_strings('*** FAILED to find initial  quantity for ***:', &
                                & fu_quantity_string(q_derived)))
      return 
    end if

    do i=1, size(q_st)
      if(q_st(i) == int_missing) exit ! all done
      ifOK = .True.
      if (req_st(i) == 2) call check_input_quantity(q_st(i), .true., &   ! q_derived static
                              & q_avail_dyn, &  ! available dynamic quantities
                              & q_avail_st, &
                              & q_shopTmp, &    ! needed to make the derivation
                              & q_shopStaticTmp, &    ! needed for the derivation
                              & ifOK, wdr, &
                              & q_intermTmp, q_intermStaticTmp, depth+1)

      if(error)return
      if( ifOK ) then
           NAvail = fu_merge_int_arrays(q_shopTmp, q_shop, .true.)  ! Fill-in arrays
           NAvailSt = fu_merge_int_arrays(q_shopStaticTmp, q_shop_st, .true.)
           if(present(q_intermediate))then
               NAvail = fu_merge_int_arrays(q_intermTmp, q_intermediate, .true.)
               NAvailSt = fu_merge_int_arrays(q_intermStaticTmp, q_intermediate_st, .true.)
           endif
       else
           call msg(fu_connect_strings('*** FAILED to find derived singletime quantity:', &
                                  & fu_quantity_string(q_st(i))))
           return
        endif
    end do

    do i=1, size(QTmp)
      if(QTmp(i) == int_missing) exit ! all done
      ifOk=.true.
      if (req(i) == 2) call check_input_quantity(QTmp(i), .false., &   ! q_derived dynamic
                              & q_avail_dyn, &  ! available dynamic quantities
                              & q_avail_st, &
                              & q_shopTmp, &    ! needed to make the derivation
                              & q_shopStaticTmp, &    ! needed for the derivation
                              & ifOK, wdr, &
                              & q_intermTmp, q_intermStaticTmp, depth+1)
      if(error)return

      if(ifOK)then
        NAvail = fu_merge_int_arrays(q_shopTmp, q_shop, .true.)  ! Fill-in arrays
        NAvailSt = fu_merge_int_arrays(q_shopStaticTmp, q_shop_st, .true.)
        if(present(q_intermediate))then
          NAvail = fu_merge_int_arrays(q_intermTmp, q_intermediate, .true.)
          NAvailSt = fu_merge_int_arrays(q_intermStaticTmp, q_intermediate_st, .true.)
        endif
      else  ! Failure
        call msg(fu_connect_strings('*** FAILED to find derived dynamic quantity:', &
                                  & fu_quantity_string(QTmp(i))))
        return
      end if
    end do

  end subroutine check_input_quantity


  !*********************************************************************

  subroutine set_deriving_way(available_list, & ! We do not know if it is st/dyn
                            & requested_dyn_list, final_dyn_list, &
                            & requested_static_list, final_static_list, &
                            & wdr)
    !
    ! Finds the set of quantities needed to derive the final set of 
    ! quantities from the input_list. Returns the overall list of 
    ! quantities, including both initial and final ones.
    ! Methodology is simple: for each quantity from final list we call 
    ! check_input_quantity, which will add necessary extra ones being limited
    ! by input_list.
    !
    implicit none

    ! Imported parameters
    type(silja_shopping_list), intent(in) :: available_list, requested_dyn_list, requested_static_list
    type(silja_shopping_list), intent(out) :: final_dyn_list, final_static_list
    type(silja_wdr), intent(in) :: wdr

    ! Local variables
    integer :: i
    integer, dimension(max_quantities) :: requests, q_shop, q_shop_static, &
                                        & q_intermediate, q_intermediate_static
    type(silja_logical), dimension(max_quantities) :: if2Ds
    logical :: found

    q_shop = int_missing
    q_intermediate = int_missing
    q_intermediate_static = int_missing
    !
    ! First, store both lists to the output.
    ! A trick. list = final_list sets all shopping parameters but does not
    ! checks the requests. So, after it one has to fix quantities using the
    ! the same final list once again. Looks stupid, but alternative is to
    ! call fix_everything routines, which is longer.
    !

    final_dyn_list = requested_dyn_list
    final_static_list = requested_static_list

    !
    ! Roux: fu_quantities(available_list,.true.) gives lots of duplicates,
    ! that are copied around. Not much harm if max_quantoties is large enough
    ! 
    !
    call fix_shopping_quantities(final_dyn_list, &
                               & fu_quantities(available_list,.true.), &
                               & fu_requests(available_list,.true.), &
                               & fu_if2Ds(available_list, .true.))
!    call msg("after fix_shopping_quantities: available")
!    call report(final_dyn_list)
!    call ooops("requested_dyn_list ")
    call add_shopping_quantities(final_dyn_list, &
                               & fu_quantities(requested_dyn_list,.true.), &
                               & fu_requests(requested_dyn_list,.true.), &
                               & fu_if2Ds(requested_dyn_list, .true.))
!    call msg("after fix_shopping_quantities: requested")
!    call report(final_dyn_list)
!    call ooops("requested_dyn_list ")
    !
    ! Cycle through the final_quantities and expanding each of its quantity to appropriate level. 
    ! Criterion is - the input_list must contain all needed quantities
    !
    do i=1, size(fu_quantities(requested_dyn_list))
      if(fu_quantity(requested_dyn_list,i) == int_missing .or. &
       & fu_request(requested_dyn_list,i) < 1) cycle ! There may be holes

      !
      ! At this stage we have no clue, if the quantity will be used as dynamic
      ! or single-time one. Provide available quantities as both....
      !
      call check_input_quantity(fu_quantity(requested_dyn_list,i), &
!                              & fu_realtime_quantity (fu_quantity(requested_dyn_list,i)),&
                              & .false.,& ! only dynamic list here !
                              & fu_quantities(available_list,.true.), & 
                              & fu_quantities(available_list,.true.), &
                              & q_shop, q_shop_static, &
                              & found, &
                              & wdr, &
                              & q_intermediate, q_intermediate_static)
      if(.not.found)then
        call set_error('Final quantity:' + &
                     & fu_quantity_short_string(fu_quantity(requested_dyn_list,i)) + &
                     & '- can not be obtained from input quantities', &
                     & 'set_deriving_way')
        return
      endif
      !
      ! Variable found - add it with the corresponding requests.
      ! Not to overlook - all intermediate quantities have to be added as well - they
      ! must be produced BEFORE the final ones can be made. So, this list must
      ! have them as well
      !
      requests = fu_request(requested_dyn_list,i)
      if2Ds = fu_if2D(requested_dyn_list,i)
      call add_shopping_quantities(final_dyn_list, q_shop, requests, if2Ds)
      call add_shopping_quantities(final_dyn_list, q_intermediate, requests, if2Ds)
      !
      if2Ds = fu_if2D(requested_static_list,i)
      call add_shopping_quantities(final_static_list, q_shop_static, requests, if2Ds)
      call add_shopping_quantities(final_static_list, q_intermediate_static, requests, if2Ds)

    end do ! Cycle through the final_quantities

    !
    ! Same for single-time (aca static) list
    ! Cycle through the final_quantities and expanding each of its quantity to appropriate level. 
    ! Criterion is - the input_list must contain all needed quantities
    !
    do i=1, size(fu_quantities(requested_static_list))
      if(fu_quantity(requested_static_list,i) == int_missing .or. &
       & fu_request(requested_static_list,i) < 1) cycle ! There may be holes
      !
      ! At this stage we have no clue, if the quantity will be used as dynamic
      ! or single-time one. Provide available quantities as both....
      !
      call check_input_quantity(fu_quantity(requested_static_list,i), &
!                              & fu_realtime_quantity (fu_quantity(requested_dyn_list,i)),&
                              & .true.,& ! static list here !
                              & fu_quantities(available_list,.true.), & 
                              & fu_quantities(available_list,.true.), &
                              & q_shop, q_shop_static, &
                              & found, &
                              & wdr, &
                              & q_intermediate, q_intermediate_static)
      if(.not.found)then
        call set_error('Final quantity:' + &
                     & fu_quantity_short_string(fu_quantity(requested_static_list,i)) + &
                     & '- can not be obtained from input quantities', &
                     & 'set_deriving_way')
        return
      endif 
      !
      ! Variable found - add it with the corresponding requests.
      ! Not to overlook - all intermediate quantities have to be added as well - they
      ! must be produced BEFORE the final ones can be made. So, this list must
      ! have them as well
      !
      requests = fu_request(requested_static_list,i)
      if2Ds = fu_if2D(requested_static_list,i)
      call add_shopping_quantities(final_dyn_list, q_shop, requests, if2Ds)
      call add_shopping_quantities(final_dyn_list, q_intermediate, requests, if2Ds)
      !
      if2Ds = fu_if2D(requested_static_list,i)
      call add_shopping_quantities(final_static_list, q_shop_static, requests, if2Ds)
      call add_shopping_quantities(final_static_list, q_intermediate_static, requests, if2Ds)

    end do ! Cycle through the final_quantities
    
  end subroutine set_deriving_way


  ! ***************************************************************

  SUBROUTINE make_derived_meteo_fields(meteoMarketPtr, list, wdr, if_high_stomatal_conductance, ifNewMeteoData)
    !
    ! This is the top interface to derived fields calculation.
    ! For all quantities with known calculation method,
    ! fields are calculated for all observation times that are found
    ! in supermarket.
    !
    ! If there is a quantity in the list, for which there is no
    ! calculation tool, it is simply skipped and no error occurs.
    ! So all the desired quantities (either from nwp external or
    ! calculated here) can be in the same list.
    ! 
    ! This routine checks that the desired quantity is not already in
    ! supermarket.
    !
    ! If some instantaneous quantities are requested, while the input
    ! is the accumulated one, proper extraction is performed
    !
    ! The supermarket is arranged after each 3d component added.
    ! So, the quantities computed later (physically later in the 
    ! routine) can already use any quantity computed in the same
    ! routine just a few lines above.
    ! Another trick: if a derived quantity requires some other ones to be made,
    ! they MUST be located above it.
    !
    ! IMPORTANT!!!. 
    ! 1. If a new quantity tool is added - do not forget to add corresponding lines 
    !    in quantities_for_derived.
    ! 2. Most derived fields are computaed only when new meteodata arrive but some
    !    are to be calculatd every time step. So, BE CAREFUL with ifNewMeteData switch
    !
    ! All units: SI
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silja_shopping_list), INTENT(IN) :: list
    type(silja_wdr), intent(in) :: wdr
    type(mini_market_of_stacks), pointer :: meteoMarketPtr
    logical, intent(in) :: if_high_stomatal_conductance
    logical, intent(in) :: ifNewMeteoData
    
    ! Local declarations:
    INTEGER :: i, t, q
    TYPE(silja_time), DIMENSION(max_times) :: valid_times
    INTEGER :: q_accum, q_average, number_of_times
    INTEGER, DIMENSION(max_quantities) :: sm_quantities_2d, sm_quantities_3d
    INTEGER :: method, NbrOf3dQ, NbrOf2dQ
    real :: fTmp
    type(meteo_data_source) :: met_src
    logical :: boolTmp1, boolTmp2, boolTmp3
    logical, save :: ifTemprMean = .false., ifFirst = .true.

    !
    ! Stupidity check
    !

    IF (.NOT.defined(list)) THEN
      CALL set_error('undefined shopping list','make_derived_meteo_fields')
      RETURN
    END IF

    met_src = fu_met_src(list)
    !
    ! Now, almost always there will be someone who will need these arrays
    !
    CALL supermarket_2d_quantities(meteoMarketPtr, met_src, multi_time_stack_flag, &
                                 & sm_quantities_2d, NbrOf2dQ)
    CALL supermarket_3d_quantities(meteoMarketPtr, met_src, multi_time_stack_flag, &
                                 & sm_quantities_3d, NbrOf3dQ)
    CALL supermarket_times(meteoMarketPtr, met_src, valid_times, number_of_times)
    IF (error) RETURN
!
!    call report(list)
!    do i = 1, NbrOf2dQ
!       call msg('sm_quantities_2d', sm_quantities_2d(i))
!    enddo
!    do i = 1, NbrOf3dQ
!       call msg('sm_quantities_3d', sm_quantities_3d(i))
!    enddo
!    do i = 1, number_of_times
!      call report(valid_times(i))
!    enddo

    !-----------------------------------------------------------------------------
    !
    ! The first grand block of derived quantities that are computed only when new meteodata show up
    ! Here we compute basics: pressure and temperature
    !
    if(ifNewMeteoData)then
      call msg('Making the following derived fields:')
      if (fu_LAIsrc(wdr) == LAI_DYNAMIC_2 ) then
         if(fu_quantity_in_list(leaf_area_index_flag, list))then ! .and. &
           if(supermarket_info) call msg_test('Making dynamic leaf-area index out of HV and LV ...')
           call dq_lai_dyn(meteoMarketPtr, met_src, valid_times, fu_LAIsrc(wdr), fu_if_randomise(wdr))
           if(error) return
           IF(.not.ANY(sm_quantities_2d == leaf_area_index_flag)) THEN
               NbrOf2dQ = NbrOf2dQ + 1
               sm_quantities_2d (NbrOf2dQ) = leaf_area_index_flag
               sm_quantities_2d (NbrOf2dQ+1) = int_missing
           END IF
         endif
      endif
      
      !
      ! Ground pressure from the logarithm of the ground pressure
      ! or mean-sea-level pressure (msl)
      !
      if(fu_quantity_in_list(ground_pressure_flag, list))then ! .and. &
        if(supermarket_info) call msg_test('Making ground pressure field...')
        call dq_ground_pressure(meteoMarketPtr, met_src, valid_times, sm_quantities_2d)
        if(error) THEN
           CALL unset_error('make_derived_meteo_fields')
           call msg("dq_ground_pressure failed. MeteoMarket contents:")
           call report(meteoMarketPtr)
           CALL set_error('dq_ground_pressure failed', 'make_derived_meteo_fields')
           return
        endif
        call arrange_supermarket_multitime(meteoMarketPtr) ! We need to set surface pressure to 3d fields
        IF(.not.ANY(sm_quantities_2d == ground_pressure_flag)) THEN
          NbrOf2dQ = NbrOf2dQ + 1
          sm_quantities_2d (NbrOf2dQ) = ground_pressure_flag
          sm_quantities_2d (NbrOf2dQ+1) = int_missing
        END IF
      endif
      !
      !  Make a surface_pressure_field using either msl_pressure or
      !  ground_pressure fields depending on the type of the system 
      !  vertical structure
      !
      if(fu_quantity_in_list(surface_pressure_flag, list))then
        if(supermarket_info) call msg_test('Selecting surface pressure field...')
        call dq_surface_pressure(meteoMarketPtr, met_src, valid_times)
        if(error) THEN
!          CALL unset_error('make_derived_meteo_fields')
           call msg("dq_surface_pressure failed. MeteoMarket contents:")
           call report(meteoMarketPtr)
           return
        ELSE
          IF(.not.ANY(sm_quantities_2d == surface_pressure_flag)) THEN
            NbrOf2dQ = NbrOf2dQ + 1
            sm_quantities_2d (NbrOf2dQ) = surface_pressure_flag
            sm_quantities_2d (NbrOf2dQ+1) = int_missing
          END IF
        END IF
      end if
      !
      ! 3D temperature field from perturbational potential temperature
      !
      IF (fu_quantity_in_list(temperature_flag, list)) then   ! .and. &
     !   & .not. fu_quantity_in_quantities(temperature_flag, sm_quantities_3d)) THEN

        IF (fu_quantity_in_quantities(perturb_pot_temperature_flag, sm_quantities_3d).and. &
            & fu_quantity_in_quantities(surface_pressure_flag, sm_quantities_2d)) THEN

          IF (supermarket_info) call msg_test('Making temperature from perturbed potential one...')

          CALL dq_temperature_from_perturb(meteoMarketPtr, met_src, valid_times)
          IF (error) THEN
            CALL unset_error('make_derived_meteo_fields')
          ELSE
            CALL arrange_supermarket_multitime(meteoMarketPtr)
            IF(.not.ANY(sm_quantities_3d == temperature_flag)) THEN
              NbrOf3dQ = NbrOf3dQ + 1
              sm_quantities_3d (NbrOf3dQ) = temperature_flag
              sm_quantities_3d (NbrOf3dQ+1) = int_missing
            END IF
          END IF
        endif
      END IF
    END IF   ! ifNewMeteodata

    ! ------------------------------------------------------------
    !
    ! 2. 3D pressure field of hybrid level data
    !
    IF (ifNewMeteoData)then
      IF(fu_quantity_in_list(pressure_flag, list) .and. &
       & (.NOT.fu_quantity_in_quantities(geopotential_flag, sm_quantities_3d))) THEN

        IF (supermarket_info) call msg_test('Making 3D hybrid level pressure...')

        CALL dq_3d_pressure(meteoMarketPtr, met_src, valid_times)
        IF (error) THEN
          CALL unset_error('make_derived_meteo_fields')
        ELSE
          CALL arrange_supermarket_multitime(meteoMarketPtr)
          IF(.not.ANY(sm_quantities_3d == pressure_flag)) THEN
            NbrOf3dQ = NbrOf3dQ + 1
            sm_quantities_3d (NbrOf3dQ) = pressure_flag
            sm_quantities_3d (NbrOf3dQ+1) = int_missing
          END IF
        END IF
      END IF

      !
      ! 3.1 Winspeed (scalar). U&V required.
      !
      IF (fu_quantity_in_list(windspeed_flag, list)) THEN

        IF (fu_quantity_in_quantities(u_flag, sm_quantities_3d).and. &
          & fu_quantity_in_quantities(v_flag, sm_quantities_3d)) THEN

          IF (supermarket_info) call msg_test('Making 3D windspeed....')

          CALL dq_windspeed(meteoMarketPtr, met_src, valid_times)
          IF (error) THEN
            CALL unset_error('make_derived_meteo_fields')
          ELSE
            CALL arrange_supermarket_multitime(meteoMarketPtr)
            IF(.not.ANY(sm_quantities_3d == windspeed_flag)) THEN
              NbrOf3dQ = NbrOf3dQ + 1
              sm_quantities_3d (NbrOf3dQ) = windspeed_flag
              sm_quantities_3d (NbrOf3dQ+1) = int_missing
            END IF
          END IF

        ELSE
          CALL msg_warning('windspeed required, but no U&V','make_derived_meteo_fields')
        END IF
      END IF

      !
      ! 3.2 Winspeed at 10m (scalar). U_10m & V_10m required.
      !
      IF (fu_quantity_in_list(windspeed_10m_flag, list)) THEN

        IF (fu_quantity_in_quantities(u_10m_flag, sm_quantities_2d).and. &
            & fu_quantity_in_quantities(v_10m_flag, sm_quantities_2d)) THEN

          IF (supermarket_info) call msg_test('Making 10m windspeed....')

          CALL dq_windspeed_10m(meteoMarketPtr, met_src, valid_times)
          IF (error) THEN
            CALL unset_error('make_derived_meteo_fields')
          ELSE
            CALL arrange_supermarket_multitime(meteoMarketPtr)
            IF(.not.ANY(sm_quantities_2d == windspeed_10m_flag)) THEN
              NbrOf2dQ = NbrOf2dQ + 1
              sm_quantities_2d (NbrOf2dQ) = windspeed_10m_flag
              sm_quantities_2d (NbrOf2dQ+1) = int_missing
            END IF
          END IF

        ELSE
          CALL msg_warning('10m windspeed requested, but no U_10m&V_10m',&
              & 'make_derived_meteo_fields')
        END IF
      END IF

      !
      ! 4. Mean wind U & V.
      !
      IF (fu_quantity_in_list(mean_wind_flag, list)) THEN

        IF (fu_quantity_in_quantities(u_flag, sm_quantities_3d).and. &
            & fu_quantity_in_quantities(v_flag, sm_quantities_3d)) THEN

          IF (supermarket_info) call msg_test('Making 3D mean windcomponents....')

          CALL dq_mean_vertically_from_bottom(meteoMarketPtr, u_flag,& 
                                            & u_mean_flag,& 
                                            & met_src,&
                                            & valid_times)

          CALL dq_mean_vertically_from_bottom(meteoMarketPtr, v_flag,& 
                                            & v_mean_flag,&
                                            & met_src,&
                                            & valid_times)
          IF (error) THEN
            CALL unset_error('make_derived_meteo_fields')
          ELSE
            CALL arrange_supermarket_multitime(meteoMarketPtr)
            IF(.not.ANY(sm_quantities_3d == u_mean_flag)) THEN
              NbrOf3dQ = NbrOf3dQ + 2
              sm_quantities_3d (NbrOf3dQ-1) = u_mean_flag
              sm_quantities_3d (NbrOf3dQ) = v_mean_flag
              sm_quantities_3d (NbrOf3dQ+1) = int_missing
            END IF
          END IF

        ELSE
          CALL msg_warning('mean wind required, but no U&V','make_derived_meteo_fields')
        END IF
      END IF

      !
      ! 5. Metric layer thickness between pressure levels and
      !    surface level. Z required.
      !
      IF (fu_quantity_in_list(layer_thickness_flag, list)) THEN

        IF (fu_quantity_in_quantities(geopotential_flag, sm_quantities_3d)) THEN

          IF (supermarket_info) call msg_test('Making metric layer thickness between pr-levels...')

          CALL dq_layer_thickness(meteoMarketPtr, met_src, valid_times)
          IF (error) THEN
            CALL unset_error('make_derived_meteo_fields')
          ELSE
            CALL arrange_supermarket_multitime(meteoMarketPtr)
            IF(.not.ANY(sm_quantities_3d == layer_thickness_flag)) THEN
              NbrOf3dQ = NbrOf3dQ + 1
              sm_quantities_3d (NbrOf3dQ) = layer_thickness_flag
              sm_quantities_3d (NbrOf3dQ+1) = int_missing
            END IF
          END IF

        ELSE
          CALL msg_warning('layer thickness required, but no Z','make_derived_meteo_fields')
        END IF
      END IF

      !
      ! 6. Metric height from ground surface.
      !
      IF (fu_quantity_in_list(height_flag, list)) THEN

        !
        ! 6.1. Metric height of pressure levels, Z and topo required
        !
        IF (fu_quantity_in_quantities(geopotential_flag, sm_quantities_3d)) THEN

          IF(supermarket_info)call msg_test('Making metric height of pressure levels...')

          CALL dq_height_from_z(meteoMarketPtr, met_src, valid_times)
          IF (error) THEN
            CALL unset_error('make_derived_meteo_fields')
          ELSE
            CALL arrange_supermarket_multitime(meteoMarketPtr)
            IF(.not.ANY(sm_quantities_3d == height_flag)) THEN
              NbrOf3dQ = NbrOf3dQ + 1
              sm_quantities_3d (NbrOf3dQ) = height_flag
              sm_quantities_3d (NbrOf3dQ+1) = int_missing
            END IF
          END IF

          !
          ! 6.2. Metric height of hybrid levels (model level data).
          !
        ELSE IF (fu_quantity_in_quantities(temperature_flag, sm_quantities_3d)) THEN
          IF(supermarket_info)call msg_test('Making metric height of hybrid levels...')

          CALL dq_height_from_t(meteoMarketPtr, met_src, valid_times)
          IF (error) THEN
            CALL unset_error('make_derived_meteo_fields')
          ELSE
            CALL arrange_supermarket_multitime(meteoMarketPtr)
            IF(.not.ANY(sm_quantities_3d == height_flag)) THEN
              NbrOf3dQ = NbrOf3dQ + 1
              sm_quantities_3d (NbrOf3dQ) = height_flag
              sm_quantities_3d (NbrOf3dQ+1) = int_missing
            END IF
          END IF

        ELSE
          CALL msg_warning('height required, but no Z of T','make_derived_meteo_fields')
        END IF
      END IF

      !
      ! Specific humidity from dew point temperature
      !
      IF (fu_quantity_in_list(specific_humidity_2m_flag, list)) THEN

        IF (fu_quantity_in_quantities(dew_point_temp_2m_flag,sm_quantities_2d)) THEN

          IF(supermarket_info)call msg_test('Making specific humidity 2m from dew point tempr 2m...')

          CALL dq_spechum_2m_from_dewp_2m(meteoMarketPtr, met_src, valid_times)
          IF (error) THEN
            CALL unset_error('make_derived_meteo_fields')
          ELSE
            CALL arrange_supermarket_multitime(meteoMarketPtr)
            IF(.not.ANY(sm_quantities_2d == specific_humidity_2m_flag)) THEN
              NbrOf2dQ = NbrOf2dQ + 1
              sm_quantities_2d (NbrOf2dQ) = specific_humidity_2m_flag
              sm_quantities_2d (NbrOf2dQ+1) = int_missing
            END IF
          END IF

        ELSE
          CALL msg_warning('Specific humidity required but no dew point temperature 2m', &
                         & 'make_derived_meteo_fields')
        END IF
      END IF

      !
      ! 7. Relative humidity from specific humidity
      !
      IF (fu_quantity_in_list(relative_humidity_flag, list)) THEN

        IF (fu_quantity_in_quantities(specific_humidity_flag,sm_quantities_3d)) THEN

          IF(supermarket_info)call msg_test('Making relative humidity from specific...')

          CALL dq_rh_from_spechum(meteoMarketPtr, met_src, valid_times)
          IF (error) THEN
            CALL unset_error('make_derived_meteo_fields')
          ELSE
            CALL arrange_supermarket_multitime(meteoMarketPtr)
            IF(.not.ANY(sm_quantities_3d == relative_humidity_flag)) THEN
              NbrOf3dQ = NbrOf3dQ + 1
              sm_quantities_3d (NbrOf3dQ) = relative_humidity_flag
              sm_quantities_3d (NbrOf3dQ+1) = int_missing
            END IF
          END IF
        ELSE
          CALL msg_warning('Rh required, but no spec.hum','make_derived_meteo_fields')
        END IF
      END IF

      !
      ! 7a. Relative humidity 2m from specific humidity 2m 
      !
      IF (fu_quantity_in_list(relative_humidity_2m_flag, list)) THEN

        IF (fu_quantity_in_quantities(specific_humidity_2m_flag,sm_quantities_2d)) THEN

          IF(supermarket_info)call msg_test('Making 2m relative humidity from specific...')

          CALL dq_rh_from_spechum_2m(meteoMarketPtr, met_src, valid_times)
          IF (error) THEN
            CALL unset_error('make_derived_meteo_fields')
          ELSE
            IF(.not.ANY(sm_quantities_2d == relative_humidity_2m_flag)) THEN
              NbrOf2dQ = NbrOf2dQ + 1
              sm_quantities_2d (NbrOf2dQ) = relative_humidity_2m_flag
              sm_quantities_2d (NbrOf2dQ+1) = int_missing
            END IF
          END IF

        ELSE
          CALL msg_warning('Rh_2m required, but no spec.hum_2m','make_derived_meteo_fields')
        END IF
      END IF
      !
      ! Total precipitation - accumulated field from large-scale and convective ones
      !
      IF (fu_quantity_in_list(total_precipitation_acc_flag, list)) THEN
        if (fu_quantity_in_quantities(large_scale_accum_rain_flag,sm_quantities_2d)) then

          IF(supermarket_info)call msg_test('Making total cumulative precipitation...')
        
          if (fu_number_of_precip_flds(wdr) == 1) then
            CALL dq_total_precipitation_ls(meteoMarketPtr, met_src, valid_times)
          else if (fu_quantity_in_quantities(convective_accum_rain_flag,sm_quantities_2d)) then
            ! Both convective and large-scale rain required and available
            CALL dq_total_precipitation(meteoMarketPtr, met_src, valid_times)
          else
            CALL msg_warning('Large-scale and convective prec required, but do not exist',&
                           & 'make_derived_meteo_fields')
          end if

          IF (error) THEN
            CALL unset_error('make_derived_meteo_fields')
          ELSE
            CALL arrange_supermarket_multitime(meteoMarketPtr)
            IF(.not.ANY(sm_quantities_2d == total_precipitation_acc_flag)) THEN
              NbrOf2dQ = NbrOf2dQ + 1
              sm_quantities_2d (NbrOf2dQ) = total_precipitation_acc_flag
              sm_quantities_2d (NbrOf2dQ+1) = int_missing
            END IF
          END IF

        ELSE
          CALL msg_warning('Large-scale prec required, but des not exist',&
                         & 'make_derived_meteo_fields')
        END IF
      END IF
    
      !
      ! 8. Equivivalent potential temperature (Pilvi's stuff).
      !
      IF (fu_quantity_in_list(eq_pot_temperature_flag, list)) THEN

        IF (supermarket_info) call msg_test('Making eq. pot. temperature...')

        CALL dq_eq_potential_temperature(meteoMarketPtr, met_src, valid_times)
        IF (error) THEN
          CALL unset_error('make_derived_meteo_fields')
        ELSE
          CALL arrange_supermarket_multitime(meteoMarketPtr)
          IF(.not.ANY(sm_quantities_3d == eq_pot_temperature_flag)) THEN
            NbrOf3dQ = NbrOf3dQ + 1
            sm_quantities_3d (NbrOf3dQ) = eq_pot_temperature_flag
            sm_quantities_3d (NbrOf3dQ+1) = int_missing
          END IF
        END IF

      END IF
      !
      ! 9. Potential temperature (Pilvi's stuff).
      !
      IF (fu_quantity_in_list(potential_temperature_flag, list)) THEN

        IF (fu_quantity_in_quantities(temperature_flag, sm_quantities_3d)) THEN

          IF (supermarket_info) call msg_test('Making potential temperature...')

          CALL dq_potential_temperature(meteoMarketPtr, met_src, valid_times)
          IF (error) THEN
            CALL unset_error('make_derived_meteo_fields')
          ELSE
            CALL arrange_supermarket_multitime(meteoMarketPtr)
            IF(.not.ANY(sm_quantities_3d == potential_temperature_flag)) THEN
              NbrOf3dQ = NbrOf3dQ + 1
              sm_quantities_3d (NbrOf3dQ) = potential_temperature_flag
              sm_quantities_3d (NbrOf3dQ+1) = int_missing
            END IF
          END IF

        ELSE
          CALL msg_warning('potential temperature required, but no T',&
              & 'make_derived_meteo_fields')
        END IF
      END IF

      !
      ! 10. Relative vorticity (Pilvi's stuff).
      !
      IF (fu_quantity_in_list(relative_vorticity_flag, list)) THEN

        IF (supermarket_info) call msg_test('Making relative vorticity....')

        CALL dq_relative_vorticity(meteoMarketPtr, met_src, valid_times)
        IF (error) THEN
          CALL unset_error('make_derived_meteo_fields')
        ELSE
          CALL arrange_supermarket_multitime(meteoMarketPtr)
          IF(.not.ANY(sm_quantities_3d == relative_vorticity_flag)) THEN
            NbrOf3dQ = NbrOf3dQ + 1
            sm_quantities_3d (NbrOf3dQ) = relative_vorticity_flag
            sm_quantities_3d (NbrOf3dQ+1) = int_missing
          END IF
        END IF

      END IF
      !
      ! 11. Absolute vorticity advection (Pilvi's stuff).
      !
      IF (fu_quantity_in_list(abs_vorticity_advection_flag, list)) THEN

        IF (supermarket_info) call msg_test('Making abs. vorticity advection...')

        CALL dq_abs_vorticity_advection(meteoMarketPtr, met_src, valid_times)
        IF (error) THEN
          CALL unset_error('make_derived_meteo_fields')
        ELSE
          CALL arrange_supermarket_multitime(meteoMarketPtr)
          IF(.not.ANY(sm_quantities_3d == abs_vorticity_advection_flag)) THEN
            NbrOf3dQ = NbrOf3dQ + 1
            sm_quantities_3d (NbrOf3dQ) = abs_vorticity_advection_flag
            sm_quantities_3d (NbrOf3dQ+1) = int_missing
          END IF
        END IF

      END IF

      !
      ! 13. Isentropic potential vorticity (Pilvi's stuff).
      !
      IF (fu_quantity_in_list(ipv_flag, list)) THEN

        IF (supermarket_info) call msg_test('Making IPV...')

        CALL dq_isentropic_pot_vorticity(meteoMarketPtr, met_src, valid_times)
        IF (error) THEN
          CALL unset_error('make_derived_meteo_fields')
        ELSE
          CALL arrange_supermarket_multitime(meteoMarketPtr)
          IF(.not.ANY(sm_quantities_3d == ipv_flag)) THEN
            NbrOf3dQ = NbrOf3dQ + 1
            sm_quantities_3d (NbrOf3dQ) = ipv_flag
            sm_quantities_3d (NbrOf3dQ+1) = int_missing
          END IF
        END IF

      END IF
      !
      ! 14b. Brunt-Vaisala frequency
      !
      IF (fu_quantity_in_list(brunt_vaisala_freq_flag, list)) THEN

        IF (supermarket_info) call msg_test('Making Brunt-Vaisala frequency...')

        CALL dq_brunt_vaisala_freq(meteoMarketPtr, met_src, valid_times)
        IF (error) THEN
          CALL unset_error('make_derived_meteo_fields')
        ELSE
          CALL arrange_supermarket_multitime(meteoMarketPtr)
          IF(.not.ANY(sm_quantities_3d == brunt_vaisala_freq_flag)) THEN
            NbrOf3dQ = NbrOf3dQ + 1
            sm_quantities_3d (NbrOf3dQ) = brunt_vaisala_freq_flag
            sm_quantities_3d (NbrOf3dQ+1) = int_missing
          END IF
        END IF

      END IF

      !
      ! 14c. Gradient Richardson number
      !
      IF (fu_quantity_in_list(gradient_richardson_nbr_flag, list)) THEN

        IF (supermarket_info) call msg_test('Making gradient Richardson number...')

        CALL dq_gradient_richardson_number(meteoMarketPtr, met_src, valid_times, fu_abl_param(wdr))
        IF (error) THEN
          CALL unset_error('make_derived_meteo_fields')
        ELSE
          CALL arrange_supermarket_multitime(meteoMarketPtr)
          IF(.not.ANY(sm_quantities_3d == gradient_richardson_nbr_flag)) THEN
            NbrOf3dQ = NbrOf3dQ + 1
            sm_quantities_3d (NbrOf3dQ) = gradient_richardson_nbr_flag
            sm_quantities_3d (NbrOf3dQ+1) = int_missing
          END IF
        END IF

      END IF

      !
      ! 15. Thermal front parameter (Pilvi's stuff).
      !
      IF (fu_quantity_in_list(tfp_flag, list)) THEN
      
        IF (supermarket_info) call msg_test('Making TFP...')

        CALL dq_thermal_front_parameter(meteoMarketPtr, met_src, valid_times)
        IF (error) THEN
          CALL unset_error('make_derived_meteo_fields')
        ELSE
          CALL arrange_supermarket_multitime(meteoMarketPtr)
          IF(.not.ANY(sm_quantities_3d == tfp_flag)) THEN
            NbrOf3dQ = NbrOf3dQ + 1
            sm_quantities_3d (NbrOf3dQ) = tfp_flag
            sm_quantities_3d (NbrOf3dQ+1) = int_missing
          END IF
        END IF

      END IF

      !
      !  vertical velocity omega (Pa/s)
      !
      IF (fu_quantity_in_list(omega_flag, list)) THEN
 
        IF (supermarket_info) call msg_test('Making vertical velocity omega...')

        CALL dq_omega(meteoMarketPtr, met_src, valid_times, .false.)  ! meteo source, time list, ifUpdate

        IF (error) THEN
          CALL unset_error('make_derived_meteo_fields')
        ELSE
          CALL arrange_supermarket_multitime(meteoMarketPtr)
          IF(.not.ANY(sm_quantities_3d == omega_flag)) THEN
            NbrOf3dQ = NbrOf3dQ + 1
            sm_quantities_3d (NbrOf3dQ) = omega_flag
            sm_quantities_3d (NbrOf3dQ+1) = int_missing
          END IF
        END IF

      END IF

      !
      !  vertical velocity w (m/s) in sea-level-based vertical coordinates
      !
      IF (fu_quantity_in_list(w_alt_msl_flag, list)) THEN
    
        IF (fu_quantity_in_quantities(omega_flag, sm_quantities_3d)) THEN

          IF (supermarket_info) call msg_test('Making vertical velocity w above msl...')

          CALL dq_vertical_velocity_msl(meteoMarketPtr, met_src, valid_times)
          IF (error) THEN
            CALL unset_error('make_derived_meteo_fields')
          ELSE
            CALL arrange_supermarket_multitime(meteoMarketPtr)
            IF(.not.ANY(sm_quantities_3d == w_alt_msl_flag)) THEN
              NbrOf3dQ = NbrOf3dQ + 1
              sm_quantities_3d (NbrOf3dQ) = w_alt_msl_flag
              sm_quantities_3d (NbrOf3dQ+1) = int_missing
            END IF
          END IF

        END IF

      END IF

      !
      !  vertical velocity w (m/s) in surface-based vertical coordinates
      !
      IF (fu_quantity_in_list(w_height_srf_flag, list)) THEN

        IF (fu_quantity_in_quantities(omega_flag, sm_quantities_3d)) THEN

          IF (supermarket_info) call msg_test('Making vertical velocity w above surface...')

          CALL dq_vertical_velocity_srf(meteoMarketPtr, met_src, valid_times)
          IF (error) THEN
            CALL unset_error('make_derived_meteo_fields')
          ELSE
            CALL arrange_supermarket_multitime(meteoMarketPtr)
            IF(.not.ANY(sm_quantities_3d == w_height_srf_flag)) THEN
              NbrOf3dQ = NbrOf3dQ + 1
              sm_quantities_3d (NbrOf3dQ) = w_height_srf_flag
              sm_quantities_3d (NbrOf3dQ+1) = int_missing
            END IF
          END IF

        END IF

      END IF

      !
      !  vertical-dependent velocity w (m/s or Pa/s)
      !  If it is omega, the above call has made it because this flag
      !  orders omega as a pre-requisite. But if it is w then it must be computed
      !
      IF (fu_quantity_in_list(vertical_velocity_flag, list)) THEN
    
        IF (supermarket_info) call msg_test('Making vertical-dependent w-wind ...')

        select case(vertical_velocity_pointer)
          case(omega_flag)
            CALL dq_omega(meteoMarketPtr, met_src, valid_times, .false.)  !meteo src, time list, ifUpdate
          case(w_alt_msl_flag)
            if(.not. fu_quantity_in_quantities(omega_flag, sm_quantities_3d)) &
                                & CALL dq_omega(meteoMarketPtr, met_src, valid_times, .false.)
            CALL dq_vertical_velocity_msl(meteoMarketPtr, met_src, valid_times)
          case(w_height_srf_flag)
            if(.not. fu_quantity_in_quantities(omega_flag, sm_quantities_3d)) &
                                & CALL dq_omega(meteoMarketPtr, met_src, valid_times, .false.)
            CALL dq_vertical_velocity_srf(meteoMarketPtr, met_src, valid_times)
          case default
            call msg('vertical_velocity_pointer',vertical_velocity_pointer)
            call set_error('Unknown vertical_velocity_pointer','make_derived_meteo_fields')
        end select

        IF (error) THEN
          CALL unset_error('make_derived_meteo_fields')
        ELSE
          CALL arrange_supermarket_multitime(meteoMarketPtr)
          CALL supermarket_3d_quantities(meteoMarketPtr, met_src, multi_time_stack_flag, &
                                       & sm_quantities_3d, NbrOf3dQ)
        END IF

      END IF

      if (fu_quantity_in_list(eta_dot_flag, list)) then
        if (supermarket_info) call msg_test('Making eta-dot...')
        call dq_hybrid_vertical_wind(meteoMarketPtr, met_src, valid_times)
      end if

      !
      !  Wind vector divergence
      !
      IF (fu_quantity_in_list(wind_divergence_flag, list)) THEN
   
        IF (supermarket_info) call msg_test('Making wind divergence...')
      
        CALL dq_wind_divergence(meteoMarketPtr, met_src, valid_times)
        IF (error) THEN
          CALL unset_error('make_derived_meteo_fields')
        ELSE
          IF(.not.ANY(sm_quantities_3d == wind_divergence_flag)) THEN
            NbrOf3dQ = NbrOf3dQ + 1
            sm_quantities_3d (NbrOf3dQ) = wind_divergence_flag
            sm_quantities_3d (NbrOf3dQ+1) = int_missing
          END IF
        END IF
      END IF

    !
    !  Bulk thermal front parameter (Pilvi's stuff).
    !
    IF (ifNewMeteoData .and. fu_quantity_in_list(bulk_tfp_flag, list)) THEN
   
      IF (.NOT.fu_quantity_in_quantities(layer_thickness_flag, sm_quantities_3d)) THEN
          CALL msg_warning('bulkTFP wanted but no layer_thickness',&
              & 'make_derived_meteo_fields')
        ELSE

          IF (supermarket_info) call msg_test('Making bulkTFP...')
        
          CALL dq_bulk_thermal_front_parameter(meteoMarketPtr, met_src, valid_times)
          IF (error) THEN
            CALL unset_error('make_derived_meteo_fields')
          ELSE
            CALL arrange_supermarket_multitime(meteoMarketPtr)
            IF(.not.ANY(sm_quantities_3d == bulk_tfp_flag)) THEN
              NbrOf3dQ = NbrOf3dQ + 1
              sm_quantities_3d (NbrOf3dQ) = bulk_tfp_flag
              sm_quantities_3d (NbrOf3dQ+1) = int_missing
            END IF
          END IF

        END IF
      END IF

      !
      !  Similarity parameters: Monin-obukhov length, friction velocity...
      !
      IF (fu_quantity_in_list(MO_length_inv_flag, list) .or. &
       & fu_quantity_in_list(friction_velocity_flag, list) .or. &
       & fu_quantity_in_list(temperature_scale_flag, list) .or. &
       & fu_quantity_in_list(humidity_scale_flag, list) .or. &
       & fu_quantity_in_list(SILAM_sensible_heat_flux_flag, list) .or. &
       & fu_quantity_in_list(SILAM_latent_heat_flux_flag, list) .or. &
       & fu_quantity_in_list(Kz_scalar_1m_flag, list) .or. &
       & fu_quantity_in_list(convective_velocity_scale_flag, list) .or.  &
       & fu_quantity_in_list(ABL_height_m_flag, list)) THEN
   
        ! Needed parameters: ground_pressure_flag, temperature_2m_flag, height_flag, 
        ! wind_10m_flag,temperature_flag,potential_temperature_flag,windspeed_flag
        ! They are not checked - in case of problem error is set in dq_monin_obukhov

          IF (supermarket_info) call msg_test('Making ABL parameters ...')

          CALL dq_ABL_params(meteoMarketPtr, met_src, valid_times, &
                                                     & fu_abl_param(wdr), fu_ablh_method(wdr))
          IF (error) THEN
            CALL unset_error('make_derived_meteo_fields')
          ELSE
            CALL arrange_supermarket_multitime(meteoMarketPtr) ! Several fields were added => full rebuild
            sm_quantities_2d(:)  = int_missing
            sm_quantities_3d(:)  = int_missing
            CALL supermarket_2d_quantities(meteoMarketPtr, met_src, multi_time_stack_flag, &
                                         & sm_quantities_2d, NbrOf2dQ)
            CALL supermarket_3d_quantities(meteoMarketPtr, met_src, multi_time_stack_flag, &
                                         & sm_quantities_3d, NbrOf3dQ)
          END IF

      END IF
      !
      ! Bulk Richardson number (Roux's stuff).
      !
      IF (fu_quantity_in_list(bulk_richardson_nbr_flag, list)) THEN

        IF (supermarket_info) call msg_test('Making bulk Richardson number...')

        CALL dq_bulk_richardson_number(meteoMarketPtr, met_src, valid_times, fu_abl_param(wdr))
        IF (error) THEN
          CALL unset_error('make_derived_meteo_fields')
        ELSE
          CALL arrange_supermarket_multitime(meteoMarketPtr)
          IF(.not.ANY(sm_quantities_3d == bulk_richardson_nbr_flag)) THEN
            NbrOf3dQ = NbrOf3dQ + 1
            sm_quantities_3d (NbrOf3dQ) = bulk_richardson_nbr_flag
            sm_quantities_3d (NbrOf3dQ+1) = int_missing
          END IF
        END IF
      END IF

      !
      !  R_a aerodynamic resistance
      !
      IF (fu_quantity_in_list(R_a_flag, list))THEN
   
        ! Needed parameters: basic surface scaling 
        ! They are not checked - in case of problem error is set in dq_monin_obukhov

        IF (supermarket_info) call msg_test('Making aerodynamic resistance...')
        
        CALL dq_aerodynamic_resistance(meteoMarketPtr, met_src, valid_times, 2.0)  ! R_a for 2m
        IF (error) THEN
          CALL unset_error('make_derived_meteo_fields')
        ELSE
          CALL arrange_supermarket_multitime(meteoMarketPtr)
          IF(.not.ANY(sm_quantities_2d == R_a_flag)) THEN
            NbrOf2dQ = NbrOf2dQ + 1
            sm_quantities_2d (NbrOf2dQ) = R_a_flag
            sm_quantities_2d (NbrOf2dQ+1) = int_missing
          END IF
        END IF

      END IF

      !
      !  ABL top pressure
      !
      IF (fu_quantity_in_list(abl_top_pressure_flag, list)) THEN

        IF (supermarket_info) call msg_test('Making ABL height [Pa]...')

        CALL dq_abl_top_pressure(meteoMarketPtr, met_src, valid_times)
        IF (error) THEN
          CALL unset_error('make_derived_meteo_fields')
        ELSE
          IF(.not.ANY(sm_quantities_2d == abl_top_pressure_flag)) THEN
            NbrOf2dQ = NbrOf2dQ + 1
            sm_quantities_2d (NbrOf2dQ) = abl_top_pressure_flag
            sm_quantities_2d (NbrOf2dQ+1) = int_missing
          END IF
        END IF

      END IF


      !
      ! Turbulence vertical profile: Kz for momentum, heat, scalar, TKE
      ! After Zilitinkevich et al, 2007
      !
      IF (fu_quantity_in_list(wind_vertical_shear_flag, list) .or. &
        & fu_quantity_in_list(Kz_momentum_3d_flag, list) .or. &
        & fu_quantity_in_list(Kz_heat_3d_flag, list) .or. &
        & fu_quantity_in_list(Kz_scalar_3d_flag, list) .or. &
!        & fu_quantity_in_list(turb_length_scale_flag, list) .or. & !Now has
!        different meaning used in R_down_meteo_flag
        & fu_quantity_in_list(turb_kinetic_energy_SILAM_flag, list)) THEN

        IF (.NOT.fu_quantity_in_quantities(gradient_richardson_nbr_flag, sm_quantities_3d)) THEN
          CALL msg_warning('Turbulence wanted but no Ri nbr','make_derived_meteo_fields')
        ELSE

          IF (supermarket_info) call msg_test('Making turbulence vertical profile...')
        
          CALL dq_turbulence_profile(meteoMarketPtr, met_src, valid_times, list)
          IF (error) THEN
            CALL unset_error('make_derived_meteo_fields')
          ELSE
            CALL arrange_supermarket_multitime(meteoMarketPtr)
            CALL supermarket_2d_quantities(meteoMarketPtr, met_src, multi_time_stack_flag, &
                                         & sm_quantities_2d, NbrOf2dQ)
            CALL supermarket_3d_quantities(meteoMarketPtr, met_src, multi_time_stack_flag, &
                                         & sm_quantities_3d, NbrOf3dQ)
          END IF

        END IF
      END IF


      IF (fu_quantity_in_list(R_down_meteo_flag, list) ) THEN
   
        ! Needed parameters: height_flag, temperature_flag, u_flag, v_flag,
        ! q_flag
        ! They are not checked - in case of problem error is set in
        ! dq_Rdown_profile

          IF (supermarket_info) call msg_test('Making dq_Rdown_profile ...')

          CALL dq_Rdown_profile(meteoMarketPtr, met_src, valid_times, fu_kz_param(wdr), fu_abl_min_m(wdr))
          IF (error) THEN
            CALL unset_error('make_derived_meteo_fields')
          ELSE
            CALL arrange_supermarket_multitime(meteoMarketPtr) 
            IF(.not.ANY(sm_quantities_3d == R_down_meteo_flag)) THEN
                 NbrOf3dQ = NbrOf3dQ + 1
                 sm_quantities_3d (NbrOf3dQ) = R_down_meteo_flag
                 sm_quantities_3d (NbrOf3dQ+1) = int_missing
            END IF
          END IF

      END IF

      IF (fu_quantity_in_list(turb_length_scale_flag, list) ) THEN
   
        ! Needed parameters: height_flag, temperature_flag, u_flag, v_flag,
        ! q_flag
        ! They are not checked - in case of problem error is set in
        ! dq_lscale_profile

          IF (supermarket_info) call msg_test('Making dq_lscale_profile ...')

          CALL dq_lscale_profile(meteoMarketPtr, met_src, valid_times, fu_kz_param(wdr), fu_abl_min_m(wdr))
          IF (error) THEN
            CALL unset_error('make_derived_meteo_fields')
          ELSE
            CALL arrange_supermarket_multitime(meteoMarketPtr) 
            IF(.not.ANY(sm_quantities_3d == turb_length_scale_flag)) THEN
                 NbrOf3dQ = NbrOf3dQ + 1
                 sm_quantities_3d (NbrOf3dQ) = turb_length_scale_flag
                 sm_quantities_3d (NbrOf3dQ+1) = int_missing
            END IF
          END IF

      END IF
            
      ! Downward radiation can be computed from net radiation and albedo
      ! use the climatological albedo here to indicate the meteo value
      ! and the cumulative quantities so the rest could work as it used to
      IF (fu_quantity_in_list(surf_sw_down_radiation_ac_flag, list)) THEN
    
        IF(fu_quantity_in_quantities(climatological_albedo_flag, sm_quantities_2d) .and. &
         & fu_quantity_in_quantities(surf_sw_net_radiation_ac_flag, sm_quantities_2d)) THEN

          IF (supermarket_info) call msg_test('Making solar down radiation...')

          CALL dq_rad_sw_down_sfc(meteoMarketPtr, met_src, valid_times)
          IF (error) THEN
            CALL unset_error('make_derived_meteo_fields')
          ELSE
            CALL arrange_supermarket_multitime(meteoMarketPtr)
            IF(.not.ANY(sm_quantities_2d == surf_sw_down_radiation_ac_flag)) THEN
              NbrOf2dQ = NbrOf2dQ + 1
              sm_quantities_2d (NbrOf2dQ) = surf_sw_down_radiation_ac_flag
              sm_quantities_2d (NbrOf2dQ+1) = int_missing
            END IF
          END IF

        ELSE
          CALL msg_warning('Missing albedo or net solar radiation for down radiation', &
                         & 'make_derived_meteo_fields')
        END IF

      END IF

      !
      !  Surface roughness field from soil and water ones
      !
      IF (fu_quantity_in_list(surface_roughness_meteo_flag, list))then

        IF (supermarket_info) call msg_test('Making surface roughness meteo ...')

        CALL dq_surface_roughness(meteoMarketPtr, .false., met_src, valid_times)
        IF (error) THEN
          CALL unset_error('make_derived_meteo_fields')
        ELSE
          IF(.not.ANY(sm_quantities_2d == water_roughness_flag)) THEN
            NbrOf2dQ = NbrOf2dQ + 1
            sm_quantities_2d (NbrOf2dQ) = water_roughness_flag
            sm_quantities_2d (NbrOf2dQ+1) = int_missing  
          END IF        
          IF(.not.ANY(sm_quantities_2d == surface_roughness_meteo_flag)) THEN
            NbrOf2dQ = NbrOf2dQ + 1
            sm_quantities_2d (NbrOf2dQ) = surface_roughness_meteo_flag
            sm_quantities_2d (NbrOf2dQ+1) = int_missing
          END IF
        END IF
      END IF
      IF (fu_quantity_in_list(surface_roughness_disp_flag, list))then

        IF (supermarket_info) call msg_test('Making surface roughness dispersion...')

        CALL dq_surface_roughness(meteoMarketPtr, .true., met_src, valid_times)
        IF (error) THEN
          CALL unset_error('make_derived_meteo_fields')
        ELSE
          IF(.not.ANY(sm_quantities_2d == surface_roughness_disp_flag)) THEN
            NbrOf2dQ = NbrOf2dQ + 1
            sm_quantities_2d (NbrOf2dQ) = surface_roughness_disp_flag
            sm_quantities_2d (NbrOf2dQ+1) = int_missing
          END IF
        END IF

      END IF

      !
      ! 23. Cumulative 2D-quantities => averages between obstimes.
      !     Temporary use of sm_quantities_*
      !     Method: if some accumulated quantity is in supermarket
      !     than corresponding average one can be created if needed
      !
      ! Create complete list of the supermarket
      sm_quantities_2d(NbrOf2dQ+1:NbrOf2dQ+NbrOf3dQ) =  sm_quantities_3d(1:NbrOf3dQ)
      NbrOf2dQ = NbrOf2dQ+NbrOf3dQ
      
      ! Check all available quantities
      DO q = 1, NbrOf2dQ
      
        IF(.not.fu_accumulated_quantity(sm_quantities_2d(q))) CYCLE
      
        q_accum = sm_quantities_2d(q)
        q_average = fu_aver_quantity_for_acc(q_accum)
        !
        ! Converting acc-av is wrong thing to do here. Skip it for rains...
        ! FIXME Should be disabled for all quantities.
        !
      
        if (q_average == large_scale_rain_int_flag .or. &
          & q_average ==  convective_rain_int_flag .or. &
          & q_average ==  total_precipitation_rate_flag .or. &
          & q_average == temperature_flag .or. &
          & q_average == temperature_2m_flag) cycle

        ! Is the available average quantity in shopping list ?
        IF (fu_quantity_in_list(q_average, list)) THEN
      
          IF(supermarket_info)THEN
            call msg_test(fu_connect_strings('Converting:', &
                                           & fu_quantity_string(q_accum),&
                                           & ', to:', &
                                           & fu_quantity_string(q_average)))
          END IF
      
          CALL dq_average_between_obstimes(meteoMarketPtr, q_accum, &
                                         & q_average, &
                                         & met_src,&
                                         & valid_times)
          IF (error) THEN
            CALL unset_error('make_derived_meteo_fields')
            CYCLE
          END IF
      
        END IF
      END DO
      
      ! Temporary variables must be restored
      
      sm_quantities_2d(:)  = int_missing
      sm_quantities_3d(:)  = int_missing
      
      CALL arrange_supermarket_multitime(meteoMarketPtr)
      CALL supermarket_2d_quantities(meteoMarketPtr, met_src, multi_time_stack_flag, &
                                   & sm_quantities_2d, NbrOf2dQ)
      CALL supermarket_3d_quantities(meteoMarketPtr, met_src, multi_time_stack_flag, &
                                   & sm_quantities_3d, NbrOf3dQ)
      !
      ! 24. Albedo can be calculated only after instantaneous radiation
      !     fluxes are known
      !
      IF (fu_quantity_in_list(albedo_flag, list)) THEN
    
        IF(fu_quantity_in_quantities(surf_sw_down_radiation_flag, sm_quantities_2d) .and. &
         & fu_quantity_in_quantities(surf_sw_net_radiation_flag, sm_quantities_2d)) THEN

          IF (supermarket_info) call msg_test('Making albedo...')

          CALL dq_albedo(meteoMarketPtr, met_src, valid_times)
          IF (error) THEN
            CALL unset_error('make_derived_meteo_fields')
          ELSE
            CALL arrange_supermarket_multitime(meteoMarketPtr)
            IF(.not.ANY(sm_quantities_2d == albedo_flag)) THEN
              NbrOf2dQ = NbrOf2dQ + 1
              sm_quantities_2d (NbrOf2dQ) = albedo_flag
              sm_quantities_2d (NbrOf2dQ+1) = int_missing
            END IF
          END IF

        ELSE
          CALL msg_warning('No instantaneous fluxes for albedo', &
                         & 'make_derived_meteo_fields')
        END IF

      END IF


      !
      ! 28. Pasquil stability class
      !
      IF (fu_quantity_in_list(pasquill_class_flag, list)) THEN
   
        IF (supermarket_info) call msg_test('Making Pasquill stability class...')
        
        CALL dq_pasquill(meteoMarketPtr, met_src, valid_times)
        IF (error) THEN
          CALL unset_error('make_derived_meteo_fields')
        ELSE
          CALL arrange_supermarket_multitime(meteoMarketPtr)
          IF(.not.ANY(sm_quantities_2d == pasquill_class_flag)) THEN
            NbrOf2dQ = NbrOf2dQ + 1
            sm_quantities_2d (NbrOf2dQ) = pasquill_class_flag
            sm_quantities_2d (NbrOf2dQ+1) = int_missing
          END IF
        END IF

      END IF

      !
      ! 29. 3D cloud water content column
      !
      boolTmp1 = fu_quantity_in_list(cwcabove_3d_flag, list) .or. fu_quantity_in_list(cwcolumn_flag, list)
      boolTmp2 = fu_quantity_in_list(pwcabove_3d_flag, list) .or. fu_quantity_in_list(pwcolumn_flag, list)
      boolTmp3 = fu_quantity_in_list(lcwcabove_3d_flag, list) .or. fu_quantity_in_list(lcwcolumn_flag, list)
      IF (boolTmp1 .or. boolTmp2 .or. boolTmp3) THEN
   
        ! Needed parameters: pressure_flag 
        ! (cloud_water_flag and cloud_ice_flag) or cloud_cond_water_flag
        ! Not checked - in case of problem error is set in dq_cwcabove_3d

        IF (supermarket_info) call msg_test('Making 3d cloud water column(s)...')

        CALL dq_cwcabove_3D(meteoMarketPtr, met_src, valid_times, boolTmp1, boolTmp2, boolTmp3)
        IF (error) THEN
          CALL set_error("Error after dq_cwcabove_3D", 'make_derived_meteo_fields')
          return
        ELSE
          CALL arrange_supermarket_multitime(meteoMarketPtr)
          IF(.not.ANY(sm_quantities_3d == cwcabove_3d_flag)) THEN
            NbrOf3dQ = NbrOf3dQ + 1
            sm_quantities_3d (NbrOf3dQ) = cwcabove_3d_flag
            sm_quantities_3d (NbrOf3dQ+1) = int_missing
          END IF
        END IF

      END IF

      !
     ! Total cloud cover from 3d
      IF (fu_quantity_in_list(total_cloud_cover_flag, list)) THEN
        IF (fu_quantity_in_quantities(cloud_cover_flag, sm_quantities_3d)) THEN

          IF(supermarket_info)call msg_test('Making total cloud coevr from 3d')

          CALL dq_cloudcvr_from_3d(meteoMarketPtr, met_src, valid_times)
          IF (error) THEN
            CALL unset_error('make_derived_meteo_fields')
          ELSE
            CALL arrange_supermarket_multitime(meteoMarketPtr)
            IF(.not.ANY(sm_quantities_2d == total_cloud_cover_flag)) THEN
              NbrOf2dQ = NbrOf2dQ + 1
              sm_quantities_2d (NbrOf2dQ) = total_cloud_cover_flag
              sm_quantities_2d (NbrOf2dQ+1) = int_missing
            END IF
          END IF

        ELSE
          CALL msg_warning('Total cloud cover required but no 3d cloud cover', &
                         & 'make_derived_meteo_fields')
        END IF
      END IF
      !
      IF (fu_quantity_in_list(stomatal_conductance_flag, list)) THEN

        IF (supermarket_info) call msg_test('Making stomatal conductance [m/s]...')

        CALL dq_stomatal_conductance(meteoMarketPtr, met_src, valid_times, fu_LAIsrc(wdr), &
             & if_high_stomatal_conductance)
        IF (error) THEN
          CALL unset_error('make_derived_meteo_fields')
        ELSE
          IF(.not.ANY(sm_quantities_2d == stomatal_conductance_flag)) THEN
            NbrOf2dQ = NbrOf2dQ + 1
            sm_quantities_2d (NbrOf2dQ) = stomatal_conductance_flag
            sm_quantities_2d (NbrOf2dQ+1) = int_missing
          END IF
        END IF
      END IF
      !
      !  Reference evapotranspiration
      !
      IF (fu_quantity_in_list(ref_evapotranspiration_flag, list)) THEN

        IF (supermarket_info) call msg_test('Making reference evapotranspiration [kg/m2sec]...')

        CALL dq_ref_evapotranspiration(meteoMarketPtr, met_src, valid_times)
        IF (error) THEN
          CALL unset_error('make_derived_meteo_fields')
        ELSE
          IF(.not.ANY(sm_quantities_2d == ref_evapotranspiration_flag)) THEN
            NbrOf2dQ = NbrOf2dQ + 1
            sm_quantities_2d (NbrOf2dQ) = ref_evapotranspiration_flag
            sm_quantities_2d (NbrOf2dQ+1) = int_missing
          END IF
        END IF
      END IF

!call report(meteoMarketPtr)
    else
       call msg("make_derived_meteo_fields: No new meteodata")
    ENDIF  ! ifNewMeteoData

!    IF(supermarket_info)THEN
!      sm_quantities_2d(:)  = int_missing
!      sm_quantities_3d(:)  = int_missing
!
!      CALL arrange_supermarket_multitime(meteoMarketPtr)
!      CALL supermarket_2d_quantities(meteoMarketPtr, met_src, sm_quantities_2d)
!      CALL supermarket_3d_quantities(meteoMarketPtr, met_src, sm_quantities_3d)
!
!      NbrOf2dQ = COUNT(sm_quantities_2d(1:SIZE(sm_quantities_2d)) /= int_missing)
!      NbrOf3dQ = COUNT(sm_quantities_3d(1:SIZE(sm_quantities_3d)) /= int_missing)
!
!      PRINT *,'Now the following 2D quantities are in SM:'
!      DO i= 1,NbrOf2dQ
!        print*,sm_quantities_2d(i),'  ',fu_quantity_string(sm_quantities_2d(i))
!      END DO
!      PRINT *,'And the following 3D quantities are in SM:'
!      DO i= 1,NbrOf3dQ
!        print*,sm_quantities_3d(i),'  ',fu_quantity_string(sm_quantities_3d(i))
!      END DO
!    END IF

 END SUBROUTINE make_derived_meteo_fields





  ! ****************************************************************
  ! ****************************************************************
  ! 
  !
  !       Private functions and subroutines of this module.
  !
  !
  ! ***************************************************************
  ! ***************************************************************




  ! ***************************************************************

  SUBROUTINE dq_surface_roughness(meteoMarketPtr, ifRoughnessDispersion, met_src, valid_times)

    ! Creates surface roughness from constant soil roughness and 
    ! dynamical water roughness

    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(mini_market_of_stacks), pointer :: meteoMarketPtr
    logical, intent(in) :: ifRoughnessDispersion
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times

    ! Local declarations:
    INTEGER ::  i,j
    TYPE(silja_field), POINTER :: land_roughness_field, water_roughness_field, &
                                & fraction_of_land_field, friction_velocity_field
    TYPE(silja_field_id) :: id
    REAL, DIMENSION(:), POINTER :: land_roughness, water_roughness, &
                                & land_roughnessTmp, &
                                & fraction_of_land, surf_r, u_star
    INTEGER :: t, iLandRoughnessFlag, iSurfaceRoughnessFlag
    TYPE(silja_time) :: time

!    surf_r => fu_work_array()
    land_roughnessTmp => fu_work_array()
    if(ifRoughnessDispersion)then
      iLandRoughnessFlag = land_roughness_disp_flag
      iSurfaceRoughnessFlag = surface_roughness_disp_flag
    else
      iLandRoughnessFlag = land_roughness_meteo_flag
      iSurfaceRoughnessFlag = surface_roughness_meteo_flag
    endif

    IF (error) RETURN

    loop_over_times: DO t = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(t))) EXIT loop_over_times
      time = valid_times(t)

      IF(fu_field_in_sm(meteoMarketPtr, met_src, &                     ! Already done before?
                      & iSurfaceRoughnessFlag,&
                      & time,&
                      & level_missing,&
                      & look_for_3d = .false.,&
                      & permanent = .false.)) CYCLE loop_over_times
      !
      ! Get or take default ground roughness
      !
      if(fu_field_in_sm(meteoMarketPtr, met_src, &
                      & iLandRoughnessFlag,&
                      & time,&
                      & level_missing,&
                      & .false.,&
                      & .false.))then
         land_roughness_field => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                                   & iLandRoughnessFlag, &
                                                   & level_missing, &
                                                   & time,&
                                                   & single_time)
        IF (error) RETURN
        land_roughness => fu_grid_data(land_roughness_field)

      elseIF(fu_field_in_sm(meteoMarketPtr, met_src, & 
                          & iLandRoughnessFlag, &
                          & time_missing,&
                          & level_missing,&
                          & .false.,&        ! look 3d
                          & .true.)) then   ! permanent
         call msg_warning('No dynamic land roughness, using static one','dq_surface_roughness')

        land_roughness_field => fu_sm_simple_field(meteoMarketPtr, met_src,&
                                                   & iLandRoughnessFlag, &
                                                   & level_missing, &
                                                   & single_time_stack_flag)
        IF (error) RETURN
        land_roughness => fu_grid_data(land_roughness_field)

      elseif(fu_field_in_sm(meteoMarketPtr, met_src, &
                      & iLandRoughnessFlag,&
                      & time,&
                      & level_missing,&
                      & .false.,&
                      & .false.))then
         call msg_warning('No static land roughness, using dynamic one','dq_surface_roughness')
         land_roughness_field => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                                   & iLandRoughnessFlag, &
                                                   & level_missing, &
                                                   & time,&
                                                   & single_time)
                      
        IF (error) RETURN
        land_roughness => fu_grid_data(land_roughness_field)

      else
        land_roughness => land_roughnessTmp
        land_roughness(1:fs_meteo) = 0.05
        call msg_warning('No land roughness, take default z0='+ fu_str(land_roughness(1))+'m', &
                       & 'dq_surface_roughness')
      endif
      !
      ! Get or create water roughness.
      !
      if(fu_field_in_sm(meteoMarketPtr, met_src, &
                      & water_roughness_flag,&
                      & time,&
                      & level_missing,&
                      & .true.,&
                      & .false.))then
        water_roughness_field => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                                   & water_roughness_flag, &
                                                   & level_missing, &
                                                   & time,&
                                                   & single_time)
        IF (error) RETURN
        water_roughness => fu_grid_data(water_roughness_field)

        ! Take a chance to copy id parameters for the new surface roughness
        !
        id = fu_set_field_id(fu_met_src(water_roughness_field),&
                           & iSurfaceRoughnessFlag,&
                           & fu_analysis_time(water_roughness_field),&
                           & fu_forecast_length(water_roughness_field), &
                           & meteo_grid,&
                           & surface_level)

      else
        !
        ! Create the water roughness from the Charnock formula
        ! z0 = Ny / (9.1 * u*) + 0.016 * (u*)**2 / g
        !
        friction_velocity_field => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                                     & friction_velocity_flag, &
                                                     & level_missing, &
                                                     & time,&
                                                     & single_time)
        if(error)return
        u_star => fu_grid_data(friction_velocity_field)
        if(error)return

        ! Store Water roughness to SM
        id = fu_set_field_id(fu_met_src(friction_velocity_field),&
                           & water_roughness_flag ,&
                           & fu_analysis_time(friction_velocity_field),&
                           & fu_forecast_length(friction_velocity_field), &
                           & meteo_grid,&
                           & surface_level)
        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, water_roughness)
        if(fu_fails(.not.error,'Failed water roughness field data pointer','dq_surface_roughness'))return

        water_roughness(1:fs_meteo) = kinematic_viscosity / (9.1 * u_star(1:fs_meteo)) + &
                                    & 0.016 *u_star(1:fs_meteo)*u_star(1:fs_meteo) / g

        ! Take a chance to copy id parameters for the new surface roughness
        !
        id = fu_set_field_id(fu_met_src(friction_velocity_field),&
                           & iSurfaceRoughnessFlag,&
                           & fu_analysis_time(friction_velocity_field),&
                           & fu_forecast_length(friction_velocity_field), &
                           & meteo_grid,&
                           & surface_level)

      endif ! Get or compute water roughness

      fraction_of_land_field => fu_sm_simple_field(meteoMarketPtr, met_src,&
                                                 & fraction_of_land_flag, &
                                                 & level_missing, &
                                                 & single_time_stack_flag)
      IF (error) RETURN

      fraction_of_land => fu_grid_data(fraction_of_land_field)

      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, surf_r)
      if(fu_fails(.not.error,'Failed roughness field data pointer','dq_surface_roughness'))return

      surf_r(1:fs_meteo) = land_roughness(1:fs_meteo) * fraction_of_land(1:fs_meteo) + &
                         & water_roughness(1:fs_meteo) * (1.- fraction_of_land(1:fs_meteo))

    END DO loop_over_times

    call free_work_array(land_roughnessTmp)

  END SUBROUTINE dq_surface_roughness



  ! ***************************************************************

  SUBROUTINE dq_ground_pressure(meteoMarketPtr, met_src, valid_times, sm_quantities_2d)

    ! Creates ground pressure from any input variable suitable for this task:
    ! logarithm of the ground pressure, msl pressure plus relief height
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr
    integer, dimension(:), intent(in) :: sm_quantities_2d

    ! Local declarations:
    INTEGER ::  i,j
    TYPE(silja_field), POINTER :: ground_pressure_fld, log_ground_pressure_fld, &
                                & msl_pressure_fld, temperature_2m_fld
    TYPE(silja_field_id) :: id
    REAL, DIMENSION(:), POINTER :: pr, log_ground_pressure, msl_pressure, t2m, topo
    INTEGER :: t
    TYPE(silja_time) :: time

    loop_over_times: DO t = 1, SIZE(valid_times)

      IF (.NOT.defined(valid_times(t))) EXIT loop_over_times
      time = valid_times(t)

      IF (fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
                       & ground_pressure_flag,&
                       & time,&
                       & level_missing,&
                       & .false.,&  ! If look for 3D
                       & .false.)) CYCLE loop_over_times  ! if Permanent
      !
      ! The deriving way depends on the available quantities
      !
      if(any(sm_quantities_2d == log_ground_pressure_flag))then
        !
        ! Make it from the logarithm
        !
        if(supermarket_info) then
          call msg_test('Making ground pressure from its logarithm for ' &
                        // trim(fu_str(time)))
        end if
        log_ground_pressure_fld => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                                     & log_ground_pressure_flag, &
                                                     & level_missing, &
                                                     & time,&
                                                     & single_time)
        IF (error) RETURN

        log_ground_pressure=> fu_grid_data(log_ground_pressure_fld)

        id = fu_set_field_id(fu_met_src(log_ground_pressure_fld),&
                           & ground_pressure_flag,&
                           & fu_analysis_time(log_ground_pressure_fld),&
                           & fu_forecast_length(log_ground_pressure_fld), &
                           & meteo_grid,&
                           & surface_level ) 

        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, pr)
        if(fu_fails(.not.error,'Failed pressure field data pointer','dq_ground_pressure'))return

        pr(1:fs_meteo) = exp(log_ground_pressure(1:fs_meteo))

      elseif(any(sm_quantities_2d == msl_pressure_flag))then
        !
        ! Make it from the msl pressure
        !
        if(supermarket_info) then
          call msg_test(fu_connect_strings('Making ground pressure from msl for:', &
                                         & fu_str(time)))
        end if
        msl_pressure_fld => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                              & msl_pressure_flag, &
                                              & level_missing, &
                                              & time,&
                                              & single_time)
        temperature_2m_fld => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                                & temperature_2m_flag, &
                                                & level_missing, &
                                                & time,&
                                                & single_time)
        IF (error) RETURN

        msl_pressure => fu_grid_data(msl_pressure_fld)
        t2m => fu_grid_data(temperature_2m_fld)
        if(associated(topography_fld))then
          topo => fu_grid_data(topography_fld)
        else
          call set_error('Topography field is not associated','dq_ground_pr_from_msl_pr')
          return
        endif
        id = fu_set_field_id(fu_met_src(temperature_2m_fld),&
                           & ground_pressure_flag,&
                           & fu_analysis_time(temperature_2m_fld),&
                           & fu_forecast_length(temperature_2m_fld), &
                           & meteo_grid,&
                           & surface_level ) 
        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, pr)
        if(fu_fails(.not.error,'Failed pressure field data pointer','dq_ground_pressure'))return

        pr(1:fs_meteo) = msl_pressure(1:fs_meteo) * exp(-topo(1:fs_meteo) * g / &
                                                      & (gas_constant_dryair * t2m(1:fs_meteo)))

      else
        call set_error('Can make ground pressure only from its logarithm or msl pressure plus relief',&
                     & 'dq_ground_pressure')
        return
      endif  ! available quantity

    END DO loop_over_times

  END SUBROUTINE dq_ground_pressure


  ! ***************************************************************

  SUBROUTINE dq_ground_pr_from_msl_pr(meteoMarketPtr, met_src, valid_times)

    ! Creates ground pressure from mean sea level one. Topography is evidently needed
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local declarations:
    INTEGER ::  i,j
    TYPE(silja_field), POINTER :: ground_pressure_field, msl_pressure_field, temperature_2m_field
    TYPE(silja_field_id) :: id
    REAL, DIMENSION(:), POINTER :: pr, msl_pressure, t2m, topo
    INTEGER :: t
    TYPE(silja_time) :: time

!    pr => fu_work_array()

    IF (error) RETURN

    loop_over_times: DO t = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(t))) EXIT loop_over_times
      time = valid_times(t)

      IF (fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
                       & ground_pressure_flag,&
                       & time,&
                       & level_missing,&
                       & .false.,&  ! If look for 3D
                       & .false.)) CYCLE loop_over_times  ! if Permanent
      if(supermarket_info) call msg_test('Making ground pressure from msl for:' + fu_str(time))

      msl_pressure_field => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                              & msl_pressure_flag, &
                                              & level_missing, &
                                              & time,&
                                              & single_time)
      temperature_2m_field => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                                & temperature_2m_flag, &
                                                & level_missing, &
                                                & time,&
                                                & single_time)
      IF (error) RETURN

      msl_pressure => fu_grid_data(msl_pressure_field)
      t2m => fu_grid_data(temperature_2m_field)
      if(associated(topography_fld))then
        topo => fu_grid_data(topography_fld)
      else
        call set_error('Topography field is not associated','dq_ground_pr_from_msl_pr')
        return
      endif

      id = fu_set_field_id(fu_met_src(msl_pressure_field),&
                         & ground_pressure_flag,&
                         & fu_analysis_time(msl_pressure_field),&
                         & fu_forecast_length(msl_pressure_field), &
                         & meteo_grid,&
                         & surface_level ) 
      if(error)return
      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, pr)
      if(fu_fails(.not.error,'Failed pressure field data pointer','dq_ground_pr_from_msl_pr'))return

      pr(1:fs_meteo) = msl_pressure(1:fs_meteo) &
           & * exp(-topo(1:fs_meteo) * g / (gas_constant_dryair * t2m(1:fs_meteo)))


!      CALL dq_store_2d(meteoMarketPtr, id, pr, multi_time_stack_flag )

    END DO loop_over_times

!    CALL free_work_array(pr)

  END SUBROUTINE dq_ground_pr_from_msl_pr



!*****************************************************************************

  subroutine dq_surface_pressure(meteoMarketPtr, met_src, valid_times)
    !
    ! Selects as a surface pressure either the mean sea level pressure (in
    ! case of the constant pressure vertical co-ordinate) or ground
    ! pressure (in case of terrain-following vertical co-ordinate).
    ! The field data are copied. Stupid looking but it may happen that
    ! some other vertical system may require more sophisticated treatment
    ! for its low-boundary field. The payment is anyway not too big.
    !
    ! Later there will be a routine for computing one pressure from another 
    ! one. So far - just selection of the correct one
    !
    implicit none

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local declarations:
    TYPE(silja_field), POINTER :: sp_field
    TYPE(silja_field_id) :: id
    REAL, DIMENSION(:), POINTER :: sp
    INTEGER :: t, fs, i, j
    TYPE(silja_time) :: time

!    sp => fu_work_array()

    loop_over_times: DO t = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(t))) EXIT loop_over_times
      time = valid_times(t)

      IF (fu_field_in_sm(meteoMarketPtr, & ! Already done before?
                       & met_src,&
                       & surface_pressure_flag,&
                       & time,&
                       & level_missing,&
                       & .false.,&
                       & .false.)) CYCLE loop_over_times

      select case(fu_leveltype(meteo_vertical))
        case(constant_pressure, constant_altitude)
          sp_field => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                        & msl_pressure_flag, &
                                        & level_missing, &
                                        & time,&
                                        & single_time)
        case(constant_height, hybrid, sigma_level, layer_btw_2_hybrid)
          sp_field => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                        & ground_pressure_flag, &
                                        & level_missing, &
                                        !& surface_level, &
                                        & time,&
                                        & single_time)
        case default
          call set_error('Unknown type of vertical structure','dq_surface_pressure')
          return
      end select

      IF (error) RETURN

      id = fu_id(sp_field)
      call set_quantity(id, surface_pressure_flag)

      fs = fu_number_of_gridpoints(fu_grid(sp_field))

      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, sp)
      if(fu_fails(.not.error,'Failed pressure field data pointer','dq_surface_pressure'))return

      sp(1:fs) = fu_grid_data(sp_field)
!      CALL dq_store_2d(meteoMarketPtr, id, sp, multi_time_stack_flag )

      if(error)return

    END DO loop_over_times

!    call free_work_array(sp)

  end subroutine dq_surface_pressure




!******************************************************************************

  SUBROUTINE dq_windspeed(meteoMarketPtr, met_src, valid_times)

    ! Description:
    ! Creates scalar windspeed 2D-fields for given one met_src and
    ! one given valid time, and stores them into the supermarket.
    ! The vertical level structure is defined by U and V -fields.
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Original code: Mika Salonoja
    ! Author: Mikhail Sofiev

    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local declarations:
    INTEGER ::  l, i,j
    TYPE(silja_3d_field), POINTER :: u3d, v3d, speed_3d, t3d 
    TYPE(silja_field_id) :: id
    REAL, DIMENSION(:), POINTER :: u, v, u_t_grid, v_t_grid, speed
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    INTEGER :: number_of_levels, t, nx, ny, i1, i2, j1, j2
    TYPE(silja_time) :: time
    real :: u_shift_lon, u_shift_lat, v_shift_lon, v_shift_lat
!    REAL, DIMENSION(2) :: u_grd_shift, v_grd_shift
!    LOGICAL :: OK
  !  TYPE(silja_grid) :: scalar_grid, u_grid, v_grid

!    speed => fu_work_array()
    u_t_grid => fu_work_array()
    v_t_grid => fu_work_array()

    IF (error) RETURN

    loop_over_times: DO t = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(t))) EXIT loop_over_times
      time = valid_times(t)

      IF (fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
                       & windspeed_flag,&
                       & time,&
                       & level_missing,&
                       & .true.,&
                       & .false.)) CYCLE loop_over_times

      u3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, u_flag, time, single_time)
      IF (error) RETURN

      v3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, v_flag, time, single_time)
      IF (error) RETURN

      t3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, temperature_flag, time, single_time)
      IF (error) RETURN

!      u => fu_grid_data_from_3d(u3d, 10)
!      v => fu_grid_data_from_3d(v3d, 10)
!      print *,'u,v before windspeed: ',u(200),v(200)

      IF (fu_number_of_fields(u3d) /=  fu_number_of_fields(v3d)) THEN
        call msg('U_ and real(V_) field sizes:',fu_number_of_fields(u3d), real(fu_number_of_fields(v3d)))
        CALL set_error('U and V sizes do not match','dq_windspeed')
        RETURN
      END IF

!!!$      IF (.NOT.(fu_grid(u3d) == fu_grid(v3d))) THEN
!!!$        CALL report(fu_grid(u3d))
!!!$        CALL report(fu_grid(v3d))
!!!$        CALL msg_warning('cannot make windspeed, U and V grid not matching',&
!!!$            & 'dq_windspeed')
!!!$        RETURN
!!!$      END IF

      CALL grid_dimensions(fu_grid(v3d), nx, ny)
      IF (error) RETURN 

      CALL vertical_levels(t3d, levels, number_of_levels)
      IF (error) RETURN

      loop_over_levels: DO l = 1, number_of_levels
       
        u => fu_grid_data_from_3d(u3d, l)
        v => fu_grid_data_from_3d(v3d, l)

        if(error)then 
          call unset_error('dq_windspeed')
          cycle
        end if

        !---------------------------------------------------------
        ! Calculations differ depending on mutual shift of the grids
        !
        u_t_grid(1:fs_meteo) = u  ! Temporal arrays for centralised wind velocity
        v_t_grid(1:fs_meteo) = v

!        CALL grid_chk(fu_grid(u3d), meteo_grid, u_grd_shift, OK)
        call grid_shift_indices(fu_grid(u3d), meteo_grid, u_shift_lon, u_shift_lat)
        IF((u_shift_lon .eps. real_missing) .or. (u_shift_lat .eps. real_missing))THEN
          call report(u3d)
          call report(meteo_grid)
          CALL set_error('Above u_grid is not close to the above meteo one','dq_speed')
          return
        END IF
!        CALL grid_chk(fu_grid(v3d), meteo_grid, v_grd_shift, OK)
        call grid_shift_indices(fu_grid(v3d), meteo_grid, v_shift_lon, v_shift_lat)
        IF((v_shift_lon .eps. real_missing) .or. (v_shift_lat .eps. real_missing))THEN
          call report(v3d)
          call report(meteo_grid)
          CALL set_error('Above v_grid is not close to the above meteo one','dq_speed')
          return
        END IF

        IF((u_shift_lon .eps. 0.).and.(u_shift_lat .eps. 0.))THEN 
            ! u grid corresponds to the system one - do nothing
        ELSEIF(u_shift_lat .eps. 0.)THEN
          DO j=1,ny   ! samelat_halfgrid_lon
          DO i=2,nx-1
            u_t_grid((j-1)*nx+i)= 0.5 * (u((j-1)*nx+i)+ u((j-1)*nx+i-NINT(SIGN(1.0,u_shift_lon))))
          END DO
          END DO
        ELSEIF(u_shift_lon .eps. 0.)THEN
          DO j=2,ny-1  ! samelon_halfgrid_lat
          DO i=1,nx
            u_t_grid((j-1)*nx+i)= 0.5*(u((j-1)*nx+i)+ u((j-1-NINT(SIGN(1.0,u_shift_lat)))*nx+i))
          END DO
          END DO
        ELSE
          DO j=2,ny-1  ! halfgrid_lon_lat - 4-point interpolation
          DO i=2,nx-1
            u_t_grid((j-1)*nx+i) = 0.25* (u((j-1)*nx+i)+ &
             & u((j-1-NINT(SIGN(1.0,u_shift_lat)))*nx+i)+ &
             & u((j-1)*nx+i-NINT(SIGN(1.0,u_shift_lon)))+ &
             & u((j-1+NINT(SIGN(1.0,u_shift_lat)))*nx+i+NINT(SIGN(1.0,u_shift_lon))))
          END DO
          END DO
        END IF

        IF((v_shift_lon .eps. 0.).and.(v_shift_lat .eps. 0.))THEN 
            ! v grid corresponds to the system one - do nothing
        ELSEIF(v_shift_lat .eps. 0.)THEN
          DO j=1,ny   ! samelat_halfgrid_lon
          DO i=2,nx-1
            v_t_grid((j-1)*nx+i)= 0.5 * (v((j-1)*nx+i)+ v((j-1)*nx+i-NINT(SIGN(1.0,v_shift_lon))))
          END DO
          END DO
        ELSEIF(v_shift_lon .eps. 0.)THEN
          DO j=2,ny-1  ! samelon_halfgrid_lat
          DO i=1,nx
            v_t_grid((j-1)*nx+i)= 0.5*(v((j-1)*nx+i)+ v((j-1-NINT(SIGN(1.0,v_shift_lat)))*nx+i))
          END DO
          END DO
        ELSE
          DO j=2,ny-1  ! halfgrid_lon_lat - 4-point interpolation
          DO i=2,nx-1
            v_t_grid((j-1)*nx+i) = 0.25* (v((j-1)*nx+i)+ &
             & v((j-1-NINT(SIGN(1.0,v_shift_lat)))*nx+i)+ &
             & v((j-1)*nx+i-NINT(SIGN(1.0,v_shift_lon)))+ &
             & v((j-1+NINT(SIGN(1.0,v_shift_lat)))*nx+i+NINT(SIGN(1.0,v_shift_lon))))
          END DO
          END DO
        END IF

        id = fu_set_field_id(fu_met_src(u3d),&
                           & windspeed_flag,&
                           & fu_analysis_time(u3d),&
                           & fu_forecast_length(u3d), &
                           & meteo_grid,&
                           & levels(l))

        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, speed)
        if(fu_fails(.not.error,'Failed pressure field data pointer','dq_windspeed'))return

        do i=1,fs_meteo
          speed(i) = SQRT(u_t_grid(i)*u_t_grid(i) + v_t_grid(i)*v_t_grid(i))
        end do

!        CALL dq_store_2d(meteoMarketPtr, id, speed, multi_time_stack_flag )

      END DO loop_over_levels

    END DO loop_over_times

!    CALL free_work_array(speed)
    CALL free_work_array(u_t_grid)
    CALL free_work_array(v_t_grid)

  END SUBROUTINE dq_windspeed



  ! ***************************************************************

  SUBROUTINE dq_windspeed_10m(meteoMarketPtr, met_src, valid_times)

    ! Creates scalar 10m windspeed field for given one met_src and
    ! one given valid time, and stores it into the supermarket.
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mikhail Sofiev

    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local declarations:
    INTEGER ::  i,j
    TYPE(silja_field), POINTER :: u_10m_field, v_10m_field, speed_10m_field
    TYPE(silja_field_id) :: id
    REAL, DIMENSION(:), POINTER :: u_10m, v_10m, u_t_grid, v_t_grid, speed_10m
    INTEGER :: t, nx, ny
    TYPE(silja_time) :: time
    real :: u_shift_lon, u_shift_lat, v_shift_lon, v_shift_lat
!    REAL, DIMENSION(2) :: u_grd_shift, v_grd_shift
!    LOGICAL :: OK

!    speed_10m => fu_work_array()
    u_t_grid => fu_work_array()
    v_t_grid => fu_work_array()

    IF (error) RETURN

    loop_over_times: DO t = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(t))) EXIT loop_over_times
      time = valid_times(t)

      IF (fu_field_in_sm(meteoMarketPtr, met_src, &           ! Already done before?
                       & windspeed_10m_flag,&
                       & time,&
                       & level_missing,&
                       & .false.,&    ! search for 3D
                       & .false.)) CYCLE loop_over_times

      u_10m_field => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                       & u_10m_flag, &
                                       & level_missing, &
                                       & time,&
                                       & single_time)
      IF (error) RETURN

      v_10m_field => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                       & v_10m_flag, &
                                       & level_missing, &
                                       & time,&
                                       & single_time)
      IF (error) RETURN

      u_10m => fu_grid_data(u_10m_field)
      v_10m => fu_grid_data(v_10m_field)
!      print *,'u,v before windspeed: ',u(200),v(200)

      CALL grid_dimensions(fu_grid(v_10m_field), nx, ny)
      IF (error) RETURN 

      !---------------------------------------------------------
      ! Calculations differ depending on mutual shift of the grids
      !
      u_t_grid(1:fs_meteo) = u_10m ! Temporal arrays for centralized wind velocity
      v_t_grid(1:fs_meteo) = v_10m

      !CALL grid_chk(fu_grid(u_10m_field), meteo_grid, u_grd_shift, OK)
      call grid_shift_indices(fu_grid(u_10m_field), meteo_grid, u_shift_lon, u_shift_lat)
      IF((u_shift_lon .eps. real_missing) .or. (u_shift_lat .eps. real_missing))THEN
      !IF(.not. OK)THEN
        CALL set_error('u_10m_grid is not close to the system one','dq_windspeed_10m')
        return
      END IF
      !CALL grid_chk(fu_grid(v_10m_field), meteo_grid, v_grd_shift, OK)
      call grid_shift_indices(fu_grid(v_10m_field), meteo_grid, v_shift_lon, v_shift_lat)
      IF((v_shift_lon .eps. real_missing) .or. (v_shift_lat .eps. real_missing))THEN
!      IF(.not. OK)THEN
        CALL set_error('v_10m_grid is not close to the system one','dq_windspeed_10m')
        return
      END IF

      IF((u_shift_lon.eps.0.).and.(u_shift_lat.eps.0.))THEN 
          ! u grid corresponds to the system one - do nothing
      ELSEIF(u_shift_lat.eps.0.)THEN
        DO j=1,ny   ! samelat_halfgrid_lon
        DO i=2,nx-1
          u_t_grid((j-1)*nx+i)= (u_10m((j-1)*nx+i+NINT(0.5+u_shift_lon))+ &
                               & u_10m((j-1)*nx+i+NINT(0.5-u_shift_lon)))* 0.5
        END DO
        END DO
      ELSEIF(u_shift_lon.eps.0.)THEN
        DO j=2,ny-1  ! samelon_halfgrid_lat
        DO i=1,nx
          u_t_grid((j-1)*nx+i)= (u_10m((j-1+NINT(0.5+u_shift_lat))*nx+i)+ &
                              &  u_10m((j-1+NINT(0.5-u_shift_lat))*nx+i))* 0.5
        END DO
        END DO
      ELSE
        DO j=2,ny-1  ! halfgrid_lon_lat - 4-point interpolation
        DO i=2,nx-1
         u_t_grid((j-1)*nx+i) = & 
          &(u_10m((j-1+NINT(0.5+u_shift_lat))*nx+i+NINT(0.5-u_shift_lon))+ &
          & u_10m((j-1+NINT(0.5-u_shift_lat))*nx+i+NINT(0.5+u_shift_lon))+ &
          & u_10m((j-1+NINT(0.5+u_shift_lat))*nx+i+NINT(0.5+u_shift_lon))+ &
          & u_10m((j-1+NINT(0.5-u_shift_lat))*nx+i+NINT(0.5-u_shift_lon)))*.25
        END DO
        END DO
      END IF

      IF((v_shift_lon.eps.0.).and.(v_shift_lat.eps.0.))THEN 
          ! u grid corresponds to the system one - do nothing
      ELSEIF(v_shift_lat.eps.0.)THEN
        DO j=1,ny   ! samelat_halfgrid_lon
        DO i=2,nx-1
          v_t_grid((j-1)*nx+i)= (v_10m((j-1)*nx+i+NINT(0.5+v_shift_lon))+ &
                               & v_10m((j-1)*nx+i+NINT(0.5-v_shift_lon)))* 0.5
        END DO
        END DO
      ELSEIF(v_shift_lon.eps.0.)THEN
        DO j=2,ny-1  ! samelon_halfgrid_lat
        DO i=1,nx
          v_t_grid((j-1)*nx+i)= (v_10m((j-1+NINT(0.5+v_shift_lat))*nx+i)+ &
                              &  v_10m((j-1+NINT(0.5-v_shift_lat))*nx+i))* 0.5
        END DO
        END DO
      ELSE
        DO j=2,ny-1  ! halfgrid_lon_lat - 4-point interpolation
        DO i=2,nx-1
         v_t_grid((j-1)*nx+i) = & 
          &(v_10m((j-1+NINT(0.5+v_shift_lat))*nx+i+NINT(0.5-v_shift_lon))+ &
          & v_10m((j-1+NINT(0.5-v_shift_lat))*nx+i+NINT(0.5+v_shift_lon))+ &
          & v_10m((j-1+NINT(0.5+v_shift_lat))*nx+i+NINT(0.5+v_shift_lon))+ &
          & v_10m((j-1+NINT(0.5-v_shift_lat))*nx+i+NINT(0.5-v_shift_lon)))*.25
        END DO
        END DO
      END IF
      !
      ! Set it up
      !
      id = fu_set_field_id(fu_met_src(u_10m_field),&
                         & windspeed_10m_flag,&
                         & fu_analysis_time(u_10m_field),&
                         & fu_forecast_length(u_10m_field), &
                         & meteo_grid,&
                         & level_10m_above_ground)
      if(error)return

      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, speed_10m)
      if(fu_fails(.not.error,'Failed pressure field data pointer','dq_windspeed'))return

      do i=1,fs_meteo
        speed_10m(i) = SQRT(u_t_grid(i)*u_t_grid(i) + v_t_grid(i)*v_t_grid(i))
      end do


!      CALL dq_store_2d(meteoMarketPtr, id, speed_10m, multi_time_stack_flag )

    END DO loop_over_times

!    CALL free_work_array(speed_10m)
    CALL free_work_array(u_t_grid)
    CALL free_work_array(v_t_grid)

  END SUBROUTINE dq_windspeed_10m



  ! ****************************************************************


  SUBROUTINE dq_layer_thickness(meteoMarketPtr, met_src, valid_times)

    ! Description:
    ! Creates scalar layer_thickness 2D-fields for given one met_src and
    ! valid times, and stores them into the supermarket.
    ! Layer thickness is calculated from the lowest level available,
    ! typically 1000HPa.
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local declarations:
    INTEGER :: i
    TYPE(silja_3d_field), POINTER :: z3d
    TYPE(silja_field_id) :: id
    REAL, DIMENSION(:), POINTER :: z, z_lowest, thickness
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    INTEGER :: number_of_levels, t
    TYPE(silja_time) :: time

!    thickness => fu_work_array()
    IF (error) RETURN

    loop_over_times: DO t = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(t))) EXIT loop_over_times
      time = valid_times(t)

      IF (fu_field_in_sm(meteoMarketPtr, & ! Already done before?
          & met_src,&
          & layer_thickness_flag,&
          & time,&
          & level_missing,&
          & .true.,&
          & .false.)) CYCLE loop_over_times

      z3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&
                                  & geopotential_flag, &
                                  & time,&
                                  & single_time)
      IF (error) RETURN

      CALL vertical_levels(z3d, levels, number_of_levels)
      IF (error) RETURN

      IF (.NOT.fu_cmp_levs_eq(levels(1), pr_level_1000hpa)) THEN
        CALL msg_warning('lowest z level not 1000HPa','dq_layer_thickness')
        CALL report(z3d)
      END IF

      z_lowest => fu_grid_data(fu_lowest_field(z3d))
      IF (error) RETURN

      loop_over_levels: DO i = 1, number_of_levels

        z => fu_grid_data_from_3d(z3d, i)

        id = fu_set_field_id(fu_met_src(z3d),&
                           & layer_thickness_flag,&
                           & fu_analysis_time(z3d),&
                           & fu_forecast_length(z3d), &
                           & fu_grid(z3d),&
                           & levels(i))
        IF (error) RETURN
        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, thickness)
        if(fu_fails(.not.error,'Failed pressure field data pointer','dq_windspeed'))return

        thickness(1:SIZE(z)) = (z - z_lowest)/g

!        CALL dq_store_2d(meteoMarketPtr, id, thickness, multi_time_stack_flag )

      END DO loop_over_levels
    END DO loop_over_times

!    CALL free_work_array(thickness)

  END SUBROUTINE dq_layer_thickness




  ! ***************************************************************


  SUBROUTINE dq_average_between_obstimes(meteoMarketPtr, quantity,& ! original accumulated quantity
                                       & average_quantity,& ! its obstime average counterpart
                                       & met_src,&
                                       & times)

    ! Description:
    ! For a quantity, which has been accumulated from the beginning
    ! of forecast, creates quantity accumulation average values
    ! between available observation times.
    ! The new quantity is also divided by its new
    ! cumulation period's length, so the  new values are like
    ! momentary ones.
    ! For example if given accumulated heat flux [J/m2] at 1h intervals,
    ! this routine creates average hourly heat flux [W/m2]
    !
    ! The given times must be in ascending order.
    !
    ! In time-series it is allowed to be fields from different model runs
    ! and thus different analysis times. In this case, however there is
    ! no way to tell the cumulation between them. In this case the latter
    ! one's average is calculated from its full length of cumulation.
    !
    ! 
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    INTEGER, INTENT(in) :: quantity, average_quantity
    type(meteo_data_source) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr


    ! Local declarations:
    TYPE(silja_field), POINTER :: latest, previous
    REAL :: aclen
    TYPE(silja_field_id) :: id
    INTEGER :: t
    REAL, DIMENSION(:), POINTER :: ave
    logical :: previous_same_start

    !----------------------------------------
    !
    ! 1. Loop over the times
    !    -------------------

    NULLIFY(latest)
    NULLIFY(previous)

!    ave => fu_work_array()
    IF (error) RETURN

    timeloop: DO t = 1, SIZE(times)

      IF (.NOT.defined(times(t))) EXIT

      latest => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                  & quantity,&
                                  & level_missing,&
                                  & times(t),&
                                  & single_time)
      IF (error) RETURN

      IF (fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
                       & average_quantity,&
                       & times(t),&
                       & level_missing,&
                       & fu_multi_level_quantity(average_quantity),&
                       & .false.)) then
        previous => latest
        CYCLE timeloop
      endif

      if(t > 1)then
        previous_same_start = fu_accumulation_start_time(latest) == &
                            & fu_accumulation_start_time(previous)
      else
        previous_same_start = .false.
      end if
      IF (error) RETURN

      id = fu_set_field_id(fu_met_src(latest),&
                         & average_quantity,&
                         & fu_analysis_time(latest),&
                         & fu_forecast_length(latest), &
                         & fu_grid(latest),&
                         & fu_level(latest))
      IF (error) RETURN

      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, ave)
      if(fu_fails(.not.error,'Failed pressure field data pointer','dq_windspeed'))return

      if(previous_same_start) then

        !----------------------------------------------------
        !
        ! 2. Previous field exists, and same start of accumulation.
        !    ------------------------------------------------

        aclen = fu_sec(fu_accumulation_length(latest) - fu_accumulation_length(previous))

        IF (aclen > 0.) THEN
          ave(1:(fu_size(latest))) = &
              & (fu_grid_data(latest) - fu_grid_data(previous))/aclen

        ELSE
          CALL report(latest)
          CALL report(previous)
          CALL set_error('zero accumulation between obstimes',&
              & 'dq_average_between_obstimes')
          RETURN
        END IF

      ELSE

        !----------------------------------------------------
        !
        ! 3. No previous field exists, or following fields
        !    from different analysis.
        !    ------------------------

        aclen = fu_sec(fu_accumulation_length(latest))

        IF (aclen > 0.) THEN

          ave(1:(fu_size(latest))) = fu_grid_data(latest)/aclen

        ELSE

          ave(1:(fu_size(latest))) = 0.

          CALL report(latest)
          CALL msg_warning('zero accumulation', 'dq_average_between_obstimes')
        END IF

      END IF


!      CALL dq_store_2d(meteoMarketPtr, id, ave, multi_time_stack_flag )
!      IF (error) RETURN

      previous => latest

    END DO timeloop

!    CALL free_work_array(ave)

  END SUBROUTINE dq_average_between_obstimes


  ! ****************************************************************


  SUBROUTINE dq_height_from_z(meteoMarketPtr, met_src, valid_times)

    ! Description:
    ! Calculates metric height from ground for constant pressure
    ! levels. Vertical coordinate defined by 3D geopotential (Z) field.
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr


    ! Local declarations:
    INTEGER :: i, t
    TYPE(silja_3d_field), POINTER :: z3d
    TYPE(silja_field_id) :: id
    REAL, DIMENSION(:), POINTER :: z, topo, height
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    INTEGER :: number_of_levels
    TYPE(silja_field), POINTER :: topofield
    TYPE(silja_time) :: time

!    height => fu_work_array()
    IF (error) RETURN

    loop_over_times: DO t = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(t))) EXIT loop_over_times
      time = valid_times(t)

      ! -------------------------------------------------------------
      !
      ! 1. Already done before?
      !
      IF (fu_field_in_sm(meteoMarketPtr, met_src,&
                       & height_flag,&
                       & time,&
                       & level_missing,&
                       & .true.,&
                       & .false.)) CYCLE loop_over_times

      ! -------------------------------------------------------------
      !
      ! 2. Get geopotential
      !
      z3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, geopotential_flag, time, single_time)
      IF (error) RETURN

      ! -------------------------------------------------------------
      !
      ! 3. For the first time find topo and check its grid.
      !
      first_time: IF (t == 1) THEN

        IF (fu_field_in_sm(meteoMarketPtr, met_src,&
                         & geopotential_flag,&
                         & time_missing,&
                         & ground_level,&
                         & .false.,&
                         & .true.)) THEN

          topofield => fu_sm_simple_field(meteoMarketPtr, met_src,&
                                        & geopotential_flag,&
                                        & ground_level, &
                                        & single_time_stack_flag) 
          IF (error) RETURN

        ELSE
          topofield => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                         & geopotential_flag,&
                                         & ground_level,&
                                         & time,&
                                         & back_and_forwards)
          IF (error) RETURN
        END IF

        IF (.NOT.(fu_grid(topofield) == fu_grid(z3d))) THEN
          CALL report(z3d)
          CALL report(topofield)
          CALL set_error('Z and topo grids not matching', 'dq_height_from_z')
          RETURN
        END IF

        topo => fu_grid_data(topofield)

      END IF first_time


      ! -------------------------------------------------------------
      !
      ! 4. Loop over levels
      !
      CALL vertical_levels(z3d, levels, number_of_levels)
      IF (error) RETURN

      loop_over_levels: DO i = 1, number_of_levels

        z => fu_grid_data_from_3d(z3d, i)

        id = fu_set_field_id(fu_met_src(z3d),&
                           & height_flag,&
                           & fu_analysis_time(z3d),&
                           & fu_forecast_length(z3d), &
                           & fu_grid(z3d),&
                           & levels(i))

        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, height)
        if(fu_fails(.not.error,'Failed pressure field data pointer','dq_windspeed'))return

        height(1:SIZE(z)) = (z-topo)/g


!        CALL dq_store_2d(meteoMarketPtr, id, height, multi_time_stack_flag )

      END DO loop_over_levels
    END DO loop_over_times

!    CALL free_work_array(height)

  END SUBROUTINE dq_height_from_z


  ! ****************************************************************


  SUBROUTINE dq_height_from_t(meteoMarketPtr, met_src, valid_times)
    !
    ! Calculates metric height from ground for hybrid (model)
    ! levels. Vertical coordinate defined by 3D temperature field.
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local declarations:
    INTEGER :: i, j, t
    TYPE(silja_3d_field), POINTER :: t3d
    TYPE(silja_field_id) :: id
    REAL, DIMENSION(:), POINTER :: t_up, &
                                 & t_down,&
                                 & p_surf,&
                                 & height, pHeightPrev, &
                                 & layer_thickness,&
                                 & p_level_up,&
                                 & p_level_down
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    INTEGER :: number_of_levels, fs, i2
    TYPE(silja_field), POINTER :: p_surf_field
    REAL :: a, b
    TYPE(silja_time) :: time

!    height => fu_work_array()
    layer_thickness => fu_work_array()
    p_level_up => fu_work_array()
    p_level_down => fu_work_array()
    IF (error) RETURN

    loop_over_times: DO t = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(t))) EXIT loop_over_times
      time = valid_times(t)

      ! -------------------------------------------------------------
      !
      ! 1. Already done before?
      !
      IF (fu_field_in_sm(meteoMarketPtr, met_src,&
                       & height_flag,&
                       & time,&
                       & level_missing,&
                       & .true.,&
                       & .false.)) CYCLE loop_over_times


      ! -------------------------------------------------------------
      !
      ! 2. Get temperature and surface pressure
      !
      t3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&
                                  & temperature_flag, &
                                  & time,&
                                  & single_time)
      IF (error) RETURN

      i = fu_leveltype(t3d)
      IF (all(i /=  (/hybrid,  sigma_level, layer_btw_2_hybrid/))) THEN
        CALL set_error('sorry, only works for hybrid level data', 'dq_height_from_t')
        RETURN
      END IF

      p_surf_field => fu_surface_pressure_field(t3d)
      IF (.NOT.defined(p_surf_field)) THEN
        CALL set_error('no surface pressure field in T 3D', 'dq_height_from_t')
        RETURN
      END IF

      p_surf => fu_grid_data(p_surf_field)
      fs = fu_size(p_surf_field)

      !
      ! 3. Loop over levels from bottom up
      !
      CALL vertical_levels(t3d, levels, number_of_levels)
      IF (error) RETURN

      DO i = 1, number_of_levels

        !      
        ! Temperatures of levels up and down:
        !
        IF (i == 1) THEN
          t_up => fu_grid_data_from_3d(t3d, 1)
          t_down => fu_grid_data_from_3d(t3d, 1)
        ELSE
          t_up => fu_grid_data_from_3d(t3d, i)
          if (any(t_up(1:fs) < 150.)) then
             call msg_warning("Temperatures  below 150K at meteo level:"+fu_str(i),'dq_height_from_t')
             call msg("Fixing then to 150K")
             WHERE(t_up(1:fs) < 150.) t_up = 150.
          endif
        END IF

        !
        ! Pressure of levels up and down in every gridpoint:
        !
        IF (i == 1) p_level_down(1:fs) = p_surf(1:fs)

        CALL pressure_on_level(t3d, i, p_level_up)
        IF (error) RETURN

#ifdef DEBUG                             
        if(any(p_level_down(1:fs) - p_level_up(1:fs) < 1.))then !1 pa increment
                j = minloc( p_level_down(1:fs) - p_level_up(1:fs) ,1) 
                call msg("Level: index of minimum thickness", i, j)
          call set_error('Non-positive pressure increment','dq_height_from_t')
          return
        end if
#endif

        !      
        ! Layer thickness between two hybrid levels or
        ! between lowest model level and ground:
        !
        !$omp workshare
        layer_thickness(1:fs) = gas_constant_dryair * (t_down + t_up) *&
                              & LOG(p_level_down(1:fs)/p_level_up(1:fs)) / (2.* g)

       !$omp end workshare
#ifdef DEBUG                             
        if (any(layer_thickness(1:fs)< 1.)) then ! Should be at least 1m
                j = minloc(layer_thickness(1:fs),1) 
                call msg("Level: index of minimum thickness", i, j)
                !                call msg("layer_thickness(1:fs)", layer_thickness(1:fs))
                call set_error("Gotcha! too thin layer", 'dq_height_from_t')
                return
        endif
#endif


        ! -------------------------------------------------
        !
        ! 3.4. Cumulate layer thickness to the height-field,
        ! and store its current status.

        id = fu_set_field_id(fu_met_src(t3d), height_flag, fu_analysis_time(t3d),&
                           & fu_forecast_length(t3d), fu_grid(t3d), levels(i))

        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, height)
        if(fu_fails(.not.error,'Failed pressure field data pointer','dq_windspeed'))return
        !
        ! Cumulate layer thickness to the height-field and store its current status.
        ! 
        if(i == 1)then
          !$omp workshare
          height(1:fs) = layer_thickness(1:fs)
          !$omp end workshare
        else
          !$omp workshare
          height(1:fs) = pHeightPrev(1:fs) + layer_thickness(1:fs)
          !$omp end workshare
        endif
        pHeightPrev => height

!        CALL dq_store_2d(meteoMarketPtr, id, height, multi_time_stack_flag )

        !
        ! Store some values for next loop
        !
        p_level_down(1:fs) = p_level_up(1:fs)
        t_down => t_up

      END DO
    END DO loop_over_times

!    call free_work_array(height)
    call free_work_array(layer_thickness)
    call free_work_array(p_level_up)
    call free_work_array(p_level_down)

  END SUBROUTINE dq_height_from_t



  ! ****************************************************************


  SUBROUTINE dq_3d_pressure(meteoMarketPtr, met_src, valid_times)

    ! Description:
    ! Calculates 3D pressure fields for hybrid (model)
    ! levels. Vertical coordinate defined by 3D temperature field.
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local declarations:
    INTEGER :: t, l
    TYPE(silja_3d_field), POINTER :: t3d
    TYPE(silja_field_id) :: id
    REAL, DIMENSION(:), POINTER :: p
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    INTEGER :: number_of_levels, fs
    TYPE(silja_time) :: time

!    p => fu_work_array()
    IF (error) RETURN

    loop_over_times: DO t = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(t))) EXIT loop_over_times
      time = valid_times(t)

      ! -------------------------------------------------------------
      !
      ! 1. Already done before?
      !    ---------------------

      IF (fu_field_in_sm(meteoMarketPtr, met_src,&
                       & pressure_flag,&
                       & time,&
                       & level_missing,&
                       & .true.,&
                       & .false.)) CYCLE loop_over_times

      ! -------------------------------------------------------------
      !
      ! 2. Get temperature and surface pressure
      !    ------------------------------------

      t3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&
                                  & temperature_flag, &
                                  & time,&
                                  & single_time)
      IF (error) RETURN

      !IF (fu_leveltype(t3d) == constant_pressure) THEN
      !  CALL set_error('do not calculate 3d pre for prlev data',& 
      !               & 'dq_3d_pressure')
      !  RETURN
      !END IF

      fs = fu_size(fu_lowest_field(t3d))
      IF (error) RETURN


      ! -------------------------------------------------------------
      !
      ! 3. Loop over levels from bottom up
      !    -------------------------------

      CALL vertical_levels(t3d, levels, number_of_levels)
      IF (error) RETURN

      DO l = 1, number_of_levels

        id = fu_set_field_id(fu_met_src(t3d),&
                           & pressure_flag,&
                           & fu_analysis_time(t3d),&
                           & fu_forecast_length(t3d), &
                           & fu_grid(t3d),&
                           & levels(l))
        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, p)
        if(fu_fails(.not.error,'Failed pressure field data pointer','dq_windspeed'))return

        CALL pressure_on_level(t3d, l, p)
        IF (error) RETURN

!        CALL dq_store_2d(meteoMarketPtr, id, p, multi_time_stack_flag )

      END DO
    END DO loop_over_times

!    call free_work_array(p)

  END SUBROUTINE dq_3d_pressure


  ! ****************************************************************


  SUBROUTINE dq_rh_from_spechum(meteoMarketPtr, met_src, valid_times)
    ! 
    ! Calculates relative humidity from specific humidity
    ! Vertical coordinate defined by 3D spec.hum field.
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local declarations:
    INTEGER :: i, tloop
    TYPE(silja_time) :: time
    TYPE(silja_3d_field), POINTER :: t3d, q3d
    TYPE(silja_field_id) :: id
    REAL, DIMENSION(:), POINTER :: t, p_surf, q, p, rh
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    INTEGER :: number_of_levels, fs

    p => fu_work_array()
    IF (error) RETURN

!    rh => fu_work_array()
!    IF (error) RETURN


    loop_over_times: DO tloop = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(tloop))) EXIT loop_over_times
      time = valid_times(tloop)

      ! -------------------------------------------------------------
      !
      ! 1. Already done before?
      !
      IF (fu_field_in_sm(meteoMarketPtr, met_src,&
                       & relative_humidity_flag,&
                       & time,&
                       & level_missing,&
                       & .true.,&
                       & .false.)) CYCLE loop_over_times

      ! -------------------------------------------------------------
      !
      ! 2. Get specific humidity, temperature and surface pressure
      !
      t3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, temperature_flag, time, single_time)
      IF (error) RETURN

      q3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, specific_humidity_flag, time, single_time)
      IF (error) RETURN

      fs = fu_size(fu_lowest_field(q3d))

      ! -------------------------------------------------------------
      !
      ! 3. Loop over levels from bottom up
      !
      CALL vertical_levels(t3d, levels, number_of_levels)
      IF (error) RETURN

      loop_over_levels: DO i = 1, number_of_levels

        t => fu_grid_data_from_3d(t3d, i)

        q => fu_grid_data_from_3d(q3d, i)

        CALL pressure_on_level(q3d, i, p)
        IF (error) RETURN

        id = fu_set_field_id(fu_met_src(t3d),&
                           & relative_humidity_flag,&
                           & fu_analysis_time(q3d),&
                           & fu_forecast_length(q3d), &
                           & fu_grid(q3d),&
                           & levels(i))

        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, rh)
        if(fu_fails(.not.error,'Failed pressure field data pointer','dq_windspeed'))return

        CALL spechum_to_relhum_field(q, t, p, rh)
        IF (error) RETURN

!        CALL dq_store_2d(meteoMarketPtr, id, rh, multi_time_stack_flag )

      END DO loop_over_levels
    END DO loop_over_times

    CALL free_work_array(p)
!    CALL free_work_array(rh)

  END SUBROUTINE dq_rh_from_spechum


  ! ****************************************************************


  SUBROUTINE dq_rh_from_spechum_2m(meteoMarketPtr, met_src, valid_times)
    ! 
    ! Calculates relative humidity at 2m from specific humidity
    !
    ! All units: SI
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local declarations:
    INTEGER :: i, tloop
    TYPE(silja_time) :: time
    TYPE(silja_field), POINTER :: t2m, p_srf, q2m
    TYPE(silja_field_id) :: id
    REAL, DIMENSION(:), POINTER :: t, p_surf, q, p, rh
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    INTEGER :: number_of_levels, fs

!    rh => fu_work_array()
    IF (error) RETURN

    loop_over_times: DO tloop = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(tloop))) EXIT loop_over_times
      time = valid_times(tloop)

      ! -------------------------------------------------------------
      !
      ! 1. Already done before?
      !
      IF (fu_field_in_sm(meteoMarketPtr, met_src,&
                       & relative_humidity_2m_flag,&
                       & time,&
                       & level_missing,&
                       & .true.,&
                       & .false.)) CYCLE loop_over_times

      ! -------------------------------------------------------------
      !
      ! 2. Get specific humidity, temperature and surface pressure
      !
      q2m => fu_sm_obstime_field(meteoMarketPtr, met_src, specific_humidity_2m_flag, &
                               & level_missing, time, single_time)
      t2m => fu_sm_obstime_field(meteoMarketPtr, met_src, temperature_2m_flag, &
                               & level_missing, time, single_time)
      p_srf => fu_sm_obstime_field(meteoMarketPtr, met_src, surface_pressure_flag, &
                                 & level_missing, time, single_time)
      IF (error) RETURN

      q => fu_grid_data(q2m)
      p => fu_grid_data(p_srf)
      t => fu_grid_data(t2m)

      id = fu_set_field_id(fu_met_src(q2m),&
                         & relative_humidity_2m_flag,&
                         & fu_analysis_time(q2m),&
                         & fu_forecast_length(q2m), &
                         & fu_grid(q2m),&
                         & fu_level(q2m))

      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, rh)
      if(fu_fails(.not.error,'Failed pressure field data pointer','dq_windspeed'))return

      CALL spechum_to_relhum_field(q, t, p, rh)
      IF (error) RETURN

!      CALL dq_store_2d(meteoMarketPtr, id, rh, multi_time_stack_flag )

    END DO loop_over_times

!    CALL free_work_array(rh)

  END SUBROUTINE dq_rh_from_spechum_2m



  ! ****************************************************************


  SUBROUTINE dq_spechum_2m_from_dewp_2m(meteoMarketPtr, met_src, valid_times)
    ! 
    ! Calculates specific humidity at 2m from dew point temperature
    !
    ! All units: SI
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local declarations:
    INTEGER :: i, tloop
    TYPE(silja_time) :: time
    TYPE(silja_field), POINTER :: td2m, p_srf, q2m
    TYPE(silja_field_id) :: id
    REAL, DIMENSION(:), POINTER :: td, p_surf, q2m_data
    INTEGER :: number_of_levels, fs

!    q2m_data => fu_work_array()
    IF (error) RETURN

    loop_over_times: DO tloop = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(tloop))) EXIT loop_over_times
      time = valid_times(tloop)

      ! -------------------------------------------------------------
      !
      ! 1. Already done before?
      !
      IF (fu_field_in_sm(meteoMarketPtr, met_src,&
                       & specific_humidity_2m_flag,&
                       & time,&
                       & level_missing,&
                       & .false.,&
                       & .false.)) CYCLE loop_over_times

      ! -------------------------------------------------------------
      !
      ! 2. Get dew point temperature and surface pressure
      !
      td2m => fu_sm_obstime_field(meteoMarketPtr, met_src, dew_point_temp_2m_flag, level_missing, time, single_time)
      p_srf => fu_sm_obstime_field(meteoMarketPtr, met_src, surface_pressure_flag, level_missing, time, single_time)
      IF (error) RETURN

      p_surf => fu_grid_data(p_srf)
      td => fu_grid_data(td2m)
      if(error)return

      fs = fu_number_of_gridpoints(fu_grid(td2m))
      if(error)return
      
      id = fu_set_field_id(fu_met_src(td2m),&
                         & specific_humidity_2m_flag,&
                         & fu_analysis_time(td2m),&
                         & fu_forecast_length(td2m), &
                         & fu_grid(td2m),&
                         & fu_level(td2m))

      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, q2m_data)
      if(fu_fails(.not.error,'Failed pressure field data pointer','dq_windspeed'))return

      do i = 1, fs
        q2m_data(i) = gas_constant_ratio * 611. / p_surf(i) * &
                             & exp(vaporization_latentheat / gas_constant_watervapour * &
                                 & (1./273. - 1./td(i)))
      end do

!      CALL dq_store_2d(meteoMarketPtr, id, q2m_data, multi_time_stack_flag )

    END DO loop_over_times

!    CALL free_work_array(q2m_data)

  END SUBROUTINE dq_spechum_2m_from_dewp_2m


  ! ***************************************************************


  SUBROUTINE dq_mean_vertically_from_bottom(meteoMarketPtr, quantity, mean_quantity, met_src, valid_times)
    !
    ! For a given scalar 3D-quantity calculates its mean value
    ! from lowest level up for every level available. The new quantity
    ! is named named as given parameter mean_quantity.
    ! 
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    INTEGER, INTENT(in) :: quantity   ! the original 3D scalar quantity
    INTEGER, INTENT(in) :: mean_quantity  ! its vertically average counterpart
    type(meteo_data_source), INTENT(in) ::  met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local declarations:
    INTEGER :: tloop, lev, fs
    REAL, DIMENSION(:), POINTER :: sum_from_bottom, mean_from_bottom, field_now
    TYPE(silja_field_id) :: id
    TYPE(silja_3d_field), POINTER :: field_3d
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    INTEGER :: number_of_levels
    TYPE(silja_time) :: time

    sum_from_bottom => fu_work_array()
!    mean_from_bottom => fu_work_array()
    field_now  => fu_work_array()
    IF (error) RETURN

    loop_over_times: DO tloop = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(tloop))) EXIT loop_over_times
      time = valid_times(tloop)


      ! -------------------------------------------------------------
      !
      ! 1. Already done before?
      !    ---------------------

      IF (fu_field_in_sm(meteoMarketPtr, & 
          & met_src,&
          & mean_quantity,&
          & time,&
          & level_missing,&
          & .true.,&
          & .false.)) CYCLE loop_over_times


      ! -------------------------------------------------------------
      !
      ! 2. Get the 3D scalar field
      !    ------------------------

      field_3d => fu_sm_obstime_3d_field(meteoMarketPtr, &
          & met_src,&
          & quantity, &
          & time,&
          & single_time)
      IF (error) RETURN

      fs = fu_size(fu_lowest_field(field_3d))


      ! -------------------------------------------------------------
      !
      ! 3. Loop over levels from bottom up
      !    -------------------------------

      CALL vertical_levels(field_3d, levels, number_of_levels)
      IF (error) RETURN

      sum_from_bottom = 0.

      loop_over_levels: DO lev = 1, number_of_levels

        field_now => fu_grid_data_from_3d(field_3d, lev)

        id = fu_set_field_id(fu_met_src(field_3d),&
                           & mean_quantity,&
                           & fu_analysis_time(field_3d),&
                           & fu_forecast_length(field_3d), &
                           & fu_grid(field_3d),&
                           & levels(lev))

        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, mean_from_bottom)
        if(fu_fails(.not.error,'Failed pressure field data pointer','dq_windspeed'))return

        sum_from_bottom(1:fs) = sum_from_bottom(1:fs) + field_now
        mean_from_bottom(1:fs) = sum_from_bottom(1:fs)/REAL(lev)


!        CALL dq_store_2d(meteoMarketPtr, id, mean_from_bottom, multi_time_stack_flag )

      END DO loop_over_levels
    END DO loop_over_times

    call free_work_array(sum_from_bottom)
!    call free_work_array(mean_from_bottom)
    call free_work_array(field_now)

  END SUBROUTINE dq_mean_vertically_from_bottom



  ! ***************************************************************


  subroutine dq_total_precipitation(meteoMarketPtr, met_src, valid_times)

    ! Creates total precipitation field for given one met_src and
    ! given valid times, and stores it into the supermarket.
    ! Not to overlook: it is an accumulated field
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mikhail Sofiev

    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local declarations:
    INTEGER ::  i,j
    TYPE(silja_field), POINTER :: ls_prec_field, conv_prec_field
    TYPE(silja_field_id) :: id
    REAL, DIMENSION(:), POINTER :: ls_prec, conv_prec, tot_prec
    INTEGER :: t, nx, ny
    TYPE(silja_time) :: time
    LOGICAL :: OK


    IF (error) RETURN

    loop_over_times: DO t = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(t))) EXIT loop_over_times
      time = valid_times(t)

      IF (fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
                       & total_precipitation_acc_flag,&
                       & time,&
                       & level_missing,&
                       & .false.,&           ! look for 3d
                       & .false.)) CYCLE loop_over_times

      ls_prec_field => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                         & large_scale_accum_rain_flag, &
                                         & level_missing, &
                                         & time,&
                                         & single_time)
      IF (error) RETURN

      conv_prec_field => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                           & convective_accum_rain_flag, &
                                           & level_missing, &
                                           & time,&
                                           & single_time)
      IF (error) RETURN

      ls_prec => fu_grid_data(ls_prec_field)
      conv_prec => fu_grid_data(conv_prec_field)

      id = fu_id(ls_prec_field)
      call set_quantity(id,total_precipitation_acc_flag)

      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, tot_prec)
      if(fu_fails(.not.error,'Failed pressure field data pointer','dq_total_precipitation'))return

      tot_prec(1:fs_meteo) = ls_prec(1:fs_meteo) + conv_prec(1:fs_meteo)

    END DO loop_over_times


  end subroutine dq_total_precipitation

  !*****************************************************************

  subroutine dq_total_precipitation_ls(meteoMarketPtr, met_src, valid_times)

    ! Creates total precipitation field for given one met_src and
    ! given valid times, and stores it into the supermarket. This
    ! subroutine uses only large-scale rain fields, so it should be
    ! called only if this is specified in wdr.
    !
    ! Not to overlook: it is an accumulated field
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mikhail Sofiev

    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local declarations:
    INTEGER ::  i,j
    TYPE(silja_field), POINTER :: ls_prec_field, conv_prec_field
    TYPE(silja_field_id) :: id
    REAL, DIMENSION(:), POINTER :: ls_prec, conv_prec, tot_prec
    INTEGER :: t, nx, ny
    TYPE(silja_time) :: time
    LOGICAL :: OK

!    tot_prec => fu_work_array()

    IF (error) RETURN

    loop_over_times: DO t = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(t))) EXIT loop_over_times
      time = valid_times(t)

      IF (fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
                       & total_precipitation_acc_flag,&
                       & time,&
                       & level_missing,&
                       & .false.,&                  ! Look for 3d
                       & .false.)) CYCLE loop_over_times

      ls_prec_field => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                         & large_scale_accum_rain_flag, &
                                         & level_missing, &
                                         & time,&
                                         & single_time)
      IF (error) RETURN

      ls_prec => fu_grid_data(ls_prec_field)

      id = fu_set_field_id(fu_met_src(ls_prec_field),&
                         & total_precipitation_acc_flag,&
                         & fu_analysis_time(ls_prec_field),&
                         & fu_forecast_length(ls_prec_field), &
                         & meteo_grid,&
                         & ground_level, &
                         & fu_accumulation_length(ls_prec_field), &
                         & field_kind = accumulated_flag)

      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, tot_prec)
      if(fu_fails(.not.error,'Failed pressure field data pointer','dq_total_precipitation_ls'))return

      tot_prec(1:fs_meteo) = ls_prec(1:fs_meteo)

!      CALL dq_store_2d(meteoMarketPtr, id, tot_prec, multi_time_stack_flag )

    END DO loop_over_times

!    CALL free_work_array(tot_prec)

  end subroutine dq_total_precipitation_ls


  !*****************************************************************

  subroutine dq_total_prec_rate(meteoMarketPtr, met_src, valid_times)

    ! Creates total precipitation field for given one met_src and
    ! given valid times, and stores it into the supermarket.
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mikhail Sofiev

    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local declarations:
    INTEGER ::  i,j
    TYPE(silja_field), POINTER :: ls_prec_field, conv_prec_field
    TYPE(silja_field_id) :: id
    REAL, DIMENSION(:), POINTER :: ls_prec, conv_prec, tot_prec
    INTEGER :: t, nx, ny
    TYPE(silja_time) :: time
    LOGICAL :: OK

!    tot_prec => fu_work_array()

    IF (error) RETURN

    loop_over_times: DO t = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(t))) EXIT loop_over_times
      time = valid_times(t)

      IF (fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
                       & total_precipitation_rate_flag,&
                       & time,&
                       & level_missing,&
                       & look_for_3d = .false.,&
                       & permanent = .false.)) CYCLE loop_over_times

      ls_prec_field => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                         & large_scale_rain_int_flag, &
                                         & level_missing, &
                                         & time,&
                                         & single_time)
      IF (error) RETURN

      conv_prec_field => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                           & convective_rain_int_flag, &
                                           & level_missing, &
                                           & time,&
                                           & single_time)
      IF (error) RETURN

      ls_prec => fu_grid_data(ls_prec_field)
      conv_prec => fu_grid_data(conv_prec_field)

      id = fu_set_field_id(fu_met_src(ls_prec_field),&
                         & total_precipitation_rate_flag,&
                         & fu_analysis_time(ls_prec_field),&
                         & fu_forecast_length(ls_prec_field), &
                         & meteo_grid,&
                         & ground_level)

      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, tot_prec)
      if(fu_fails(.not.error,'Failed pressure field data pointer','dq_total_prec_rate'))return

      tot_prec(1:fs_meteo) = ls_prec(1:fs_meteo) + conv_prec(1:fs_meteo)

!      CALL dq_store_2d(meteoMarketPtr, id, tot_prec, multi_time_stack_flag )

    END DO loop_over_times

!    CALL free_work_array(tot_prec)

  end subroutine dq_total_prec_rate

  !*****************************************************************

  subroutine dq_total_prec_rate_ls(meteoMarketPtr, met_src, valid_times)

    ! Creates total precipitation field for given one met_src and
    ! given valid times, and stores it into the supermarket. This
    ! subroutine uses only large-scale rain fields, so it should be
    ! called only if this is specified in wdr.
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mikhail Sofiev

    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local declarations:
    INTEGER ::  i,j
    TYPE(silja_field), POINTER :: ls_prec_field, conv_prec_field
    TYPE(silja_field_id) :: id
    REAL, DIMENSION(:), POINTER :: ls_prec, conv_prec, tot_prec
    INTEGER :: t, nx, ny
    TYPE(silja_time) :: time
    LOGICAL :: OK

!    tot_prec => fu_work_array()

    IF (error) RETURN

    loop_over_times: DO t = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(t))) EXIT loop_over_times
      time = valid_times(t)

      IF (fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
                       & total_precipitation_rate_flag,&
                       & time,&
                       & level_missing,&
                       & .true.,&
                       & .false.)) CYCLE loop_over_times

      ls_prec_field => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                         & large_scale_rain_int_flag, &
                                         & level_missing, &
                                         & time,&
                                         & single_time)
      IF (error) RETURN

      ls_prec => fu_grid_data(ls_prec_field)

      id = fu_set_field_id(fu_met_src(ls_prec_field),&
                         & total_precipitation_rate_flag,&
                         & fu_analysis_time(ls_prec_field),&
                         & fu_forecast_length(ls_prec_field), &
                         & meteo_grid,&
                         & ground_level)

      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, tot_prec)
      if(fu_fails(.not.error,'Failed pressure field data pointer','dq_total_prec_rate_ls'))return

      tot_prec(1:fs_meteo) = ls_prec(1:fs_meteo)

                         
!      CALL dq_store_2d(meteoMarketPtr, id, tot_prec, multi_time_stack_flag )

    END DO loop_over_times

!    CALL free_work_array(tot_prec)

  end subroutine dq_total_prec_rate_ls

  

  ! ***************************************************************

  ! ***************************************************************
  ! 
  !
  !                 TESTS
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  SUBROUTINE dq_tests(meteoMarketPtr, post_processed)

    IMPLICIT NONE
    LOGICAL, INTENT(in)  :: post_processed
    CHARACTER (LEN=200) :: fname
    type(meteo_data_source) :: met_src
    TYPE(silja_grid) :: grid
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

!!!$    CALL set_supermarket_storage_grid(&
!!!$        & fu_set_latlon_grid(&
!!!$        & 'harva_testi',&
!!!$        &  0., 0.,&
!!!$        & .false.,&
!!!$        & 40, 40,&
!!!$        & pole_system,&
!!!$        & 0.4, 0.4))
!!!$

    met_src = fmi_hirlam_src
    fname = 'test.v5d'

    CALL supermarket_test_fill(meteoMarketPtr) !post_processed)
    IF (error) RETURN

    CALL make_derived_meteo_fields(meteoMarketPtr, &
                                 & fu_set_shopping_list (met_src,&
                                                       & (/temperature_flag,&
                                                         & u_flag, v_flag,&
                                                         & pressure_flag,&
                                                         & windspeed_flag/),&
                                                       & time_missing, &
                                                       & time_missing, &
                                                       & level_missing,&
                                                       & level_missing), &
                                 & wdr_missing, &
                                 & .true., &
                                 & .true.)
    IF (error) RETURN

!!!$    CALL make_derived_meteo_fields(&
!!!$        & fu_set_shopping_list (&
!!!$        & met_src,&
!!!$        & post_processed,&
!!!$        & (/accept_all_quantities/),&
!!!$        & time_missing, &
!!!$        & time_missing, &
!!!$        & level_missing,&
!!!$        & level_missing))
!!!$    IF (error) RETURN


!    CALL report(meteoMarketPtr) 
!    IF (error) RETURN
    
    CALL arrange_supermarket_multitime(meteoMarketPtr)

!    CALL  report(meteoMarketPtr)

!    CALL supermarket_to_v5d_file(met_src, (/accept_all_quantities/), fname)


  END SUBROUTINE dq_tests

END MODULE derived_field_quantities


