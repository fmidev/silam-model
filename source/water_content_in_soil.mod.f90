module water_in_soil
  !
  ! A model (later, models) for water content in soil. 
  ! Model 1: Linkosale,T., Pasi Kolari,P., Pumpanen, J. (2013) New decomposition rate functions based
  ! on volumetric soil water ontent for the ROMUL soil organic matter dynamics model. Ecological 
  ! Modelling 263 (2013) 109–118.
  !
  ! The model integrates the soil water content using precipitation rates, temperature and humidity 
  ! in air. The model has two layers, the upper layer with organic content (may be) and lower 
  ! mineral layer. Precipitation fill-in the upper bucket, which drains the excess down to the 
  ! lower layer, which, in turn, drains the excess down to nowhere.
  ! Actual water content is stored as two single-time fields in the dispersion stack.
  !
  ! Code author: M.Sofiev
  !
  use field_buffer
 !  use silam_levels

  implicit none
  
  public set_rules_WIS
  public init_WIS
  public add_input_needs_WIS
  public update_water_in_soil
  public defined

  private defined_rulesWIS
  
  interface defined
    module procedure defined_rulesWIS
  end interface
  
  type Trules_water_in_soil
    private
    integer :: soil_water_model, evaporation_model, nRainFlds
    real :: fill_value, tauDrain_up2down, evapPlant_Scale, fract_roots_up
    type(silja_level), dimension(2) :: soil_levels
    type(Tsilam_namelist), pointer :: nlInputFiles  ! namelist for names of static meteo files
    type(silja_logical) :: defined
  end type Trules_water_in_soil
  public Trules_water_in_soil

  type(Trules_water_in_soil), public, parameter :: rulesWIS_missing = &
      & Trules_water_in_soil(int_missing, int_missing, int_missing, &
                           & real_missing, real_missing, real_missing, real_missing, &
                           & (/level_missing, level_missing/), null(), silja_false)
  !
  ! Evaporation and soil moisture model types
  !
  integer, private, parameter :: empirical_evap_rate_flag = 6100    
  integer, private, parameter :: latent_heat_proxy_evap_flag = 6101
  
  integer, private, parameter :: soil_water_nwp_flag = 6110
  integer, private, parameter :: soil_water_basic_mdl_flag = 6111
  
  
  CONTAINS
  
  
  
  !*****************************************************************************
  
  subroutine set_rules_WIS(nlSetup, rulesWIS)
    !
    ! Sets rules of making the soil water content. Note that by the moment these rules are set
    ! the necessity of these variables are not known. Therefore undefined rules as a result of 
    ! this sub is not an error.
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist), pointer :: nlSetup
    type(Trules_water_in_soil), intent(out) :: rulesWIS
    
    ! Local varaibles
    integer :: iFile, nFiles
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: pItems

    
    rulesWIS = rulesWIS_missing        ! precaution: some pieces may be undefined

    ! 
    ! Two possibilities for the soil water: take it from NWP or compute as suggested in this module
    !
    if(fu_str_l_case(fu_content(nlSetup,'soil_water_model')) == 'nwp')then
      rulesWIS%soil_water_model = soil_water_nwp_flag
    elseif(fu_str_l_case(fu_content(nlSetup,'soil_water_model')) == 'basic_dynamics')then
      rulesWIS%soil_water_model = soil_water_basic_mdl_flag
    else
      call msg('Unknown soil_water_model:>>' + fu_content(nlSetup,'soil_water_model') + '<<, need it?')
      call msg("Using NWP")
      rulesWIS%soil_water_model = soil_water_nwp_flag
      return
    endif
    ! 
    ! Two possibilities for the total water flux: through the latent heat and through the 
    ! semi-empirical relation of Linkosalo ea paper
    !
    if(fu_str_l_case(fu_content(nlSetup,'evaporation_model')) == 'latent_heat')then
      rulesWIS%evaporation_model = latent_heat_proxy_evap_flag
    elseif(fu_str_l_case(fu_content(nlSetup,'evaporation_model')) == 'empirical_rate')then
      rulesWIS%evaporation_model = empirical_evap_rate_flag
    else
      call msg('Unknown evaporation_model:>>' + fu_content(nlSetup,'evaporation_model') +'>>, need it?')
      return
    endif
    !
    ! Store the input fields for the rules
    !
    rulesWIS%nlInputFiles => fu_create_namelist('water_in_soil_rules_supplementary_files')
    call get_items(nlSetup, 'static_meteo_file', pItems, nFiles)
    if(error)return
    do iFile = 1, nFiles
      call add_namelist_item(rulesWIS%nlInputFiles, 'static_meteo_file', &
                           & fu_process_filepath(fu_content(pItems(iFile))))
    end do
    
    !
    ! soil water levels
    !
    rulesWIS%soil_levels(1:2) = &
                   & (/fu_set_layer_between_two(fu_set_depth_level(0.0), fu_set_depth_level(0.1)), &
                     & fu_set_layer_between_two(fu_set_depth_level(0.1), fu_set_depth_level(0.5)) /)


    rulesWIS%defined = silja_true
    
  end subroutine set_rules_WIS
  
    
  !*****************************************************************************
  
  subroutine init_WIS(dispMarketPtr, meteoMarketPtr, rulesWIS, start_time)
    !
    ! Initialises the Water_in_soil model from the setup namelist
    !
    implicit none

    ! Imported parameters
    type(mini_market_of_stacks), pointer :: dispMarketPtr, meteoMarketPtr
    type(Trules_water_in_soil), intent(in) :: rulesWIS
    type(silja_time), intent(in) :: start_time

    ! Local variables
    type(silja_shopping_list), target :: shop_list
    type(silja_shopping_list), pointer :: pShop_list
    type(silam_vertical) :: vertTmp
    integer, dimension(:), pointer ::  q_disp_dyn, q_disp_stat, stack_quantities
    integer, dimension(2) :: water_in_soil_flags = (/water_in_soil_srf_grav_flag, water_in_soil_deep_grav_flag/)
    integer :: iFlds
    logical :: ifOK
    type(silja_field_id) :: id
    type(silja_field), pointer :: pField
    real, dimension(:), pointer :: pVals
    !
    ! If we are here, the variables are requested, thus rules must be defined
    !
    if(fu_fails(defined(rulesWIS),'undefined water_in_soil rules but variables needed','init_WIS'))return
    
    pShop_list => shop_list
    !
    ! Request the water_in_soil fields: two are needed
    !
    q_disp_dyn => fu_work_int_array()
    q_disp_stat => fu_work_int_array()
    stack_quantities => fu_work_int_array()
    if(error)return
    q_disp_dyn(1:max_quantities) = int_missing
    q_disp_stat(1:max_quantities) = int_missing
    stack_quantities(1:max_quantities) = int_missing
    call set_missing(shop_list)
    call set_missing(vertTmp, .true.)
    if(error)return
    !
    ! Get the static quantities to the single-time stack. If they are already there, do nothing
    !
    call add_input_needs_WIS(rulesWIS, stack_quantities, stack_quantities, & ! meteo quantities, skip
                                     & q_disp_dyn, q_disp_stat)   ! dispersion-buffer quantities, use
    if(error)return
    
    do iFlds = 1, size(q_disp_stat)
      if(q_disp_stat(iFlds) == int_missing)exit
      call add_shopping_variable(shop_list, &
                               & q_disp_stat(iFlds), &     ! quantity
                               & species_missing, &
                               & grid_missing, &
                               & vertTmp, int_missing, &
                               & met_src_missing)
      if(error)return
    end do
    !
    ! Apart from the input fields, we need the water in soil fields themselves. They can be - 
    ! and better be - initialized by reading something from static_meteo_file. Climatology,
    ! for instance. But then we have to limit the validity time range
    !
    call set_vertical(fu_set_layer_between_two(fu_set_depth_level(0.0), &  ! upper layer, 0.1m
                                             & fu_set_depth_level(0.1)), vertTmp)
    if(error)return
    call add_level(vertTmp, fu_set_layer_between_two(fu_set_depth_level(0.1), &  ! lower layer, 0.1-0.5m
                                                   & fu_set_depth_level(0.5)))
    if(error)return
    do iFlds = 1, 2
      call add_shopping_variable(shop_list, &
                               & water_in_soil_flags(iFlds), &
                               & species_missing, &
                               & grid_missing, &
                               & vertTmp, iFlds, &
                               & met_src_missing)
      call fix_shopping_time_boundaries(pShop_list, start_time - one_minute, start_time + one_minute)
      if(error)return
    end do
    !
    ! The supplementary fields have to be taken from static_meteo_file items. In principle,
    ! it better be initialised: otherwise the field will be set to zero, i.e. the domain will 
    ! bone-dry.
    !
    call msg('Filling-in the water-in-soil input fields from static_meteo_file items')
    call fill_minimarket_from_namelist(meteoMarketPtr, &
                                     & rulesWIS%nlInputFiles, 'static_meteo_file', & ! namelist and item
                                     & shop_list, start_time, &
                                     & static_climatology, &  !make ever valid
                                     & create_if_absent, &        ! no error if a clash
                                     & wdr_missing, &
                                     & dispersion_gridPtr, &
                                     & 5, .true., & ! iAccuracy, ifAdjustGrid
                                     & ifOK)
    !
    ! Now, create the water fields in dispersion market. Two fields, one for near-surface, one for deep layer
    ! Their metrical depth is of no importance, it is measured in water holding capacity. 
    !
    call msg('creating/reading the water-in-soil fields')

    do iFlds = 1, 2
      id = fu_set_field_id_simple(met_src_missing,&
                                & water_in_soil_flags(iFlds), &
                                & time_missing, &        ! valid time
                                & fu_level(vertTmp,iFlds), &
                                & species_missing)
      pField => fu_get_field_from_mm_general(meteoMarketPtr, id, .false.)
      if(associated(pField))then
        ifOK= defined(pField)    ! exists already
      else
        ifOK = .false.   ! not found in the market
        !
        ! Make the place. If fill_value is given, make it up!
        !
        id = fu_set_field_id(met_src_missing,&
                             & water_in_soil_flags(iFlds), &
                             & start_time, &        ! analysis time
                             & zero_interval, &     ! forecast length
                             & dispersion_grid,&    ! grid
                             & fu_level(vertTmp,iFlds), &       ! level
                             & zero_interval, &     ! length of accumulation
                             & interval_missing, &     ! length of validity
                             & forecast_flag, &     ! field_kind
                             & species = species_missing) ! species
        !
        ! Actually makes the space
        !
        call find_field_data_storage_2d(meteoMarketPtr, id, single_time_stack_flag, pVals)
        if(error)return
        
        pVals(:) = rulesWIS%fill_value

      endif   ! if field is in market
    end do  ! Cycle over two depth layers to read
    
    call free_work_array(q_disp_dyn)
    call free_work_array(q_disp_stat)
    call free_work_array(stack_quantities)

  end subroutine init_WIS
  
  
  !*****************************************************************************
  
  subroutine add_input_needs_WIS(rulesWIS, q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st)
    !
    ! Returns input needs for the water_in_soil model. 
    !
    implicit none

    ! Imported parameters
    type(Trules_water_in_soil), intent(in) :: rulesWIS
    integer, dimension(:), intent(inout) :: q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st

    ! Local variables
    integer :: iTmp

    !
    ! If we use NWP input for soil water, therre is an appropriate flag
    !
    select case(rulesWIS%soil_water_model)
      case(soil_water_nwp_flag)
        iTmp = fu_merge_integer_to_array(soil_moisture_vol_frac_nwp_flag,    q_met_dynamic)
      case(soil_water_basic_mdl_flag)
        ! Meteo dynamic
        iTmp = fu_merge_integer_to_array(relative_humidity_flag,    q_met_dynamic)
        iTmp = fu_merge_integer_to_array(temperature_2m_flag,       q_met_dynamic)
        ! Meteo single-time
        iTmp = fu_merge_integer_to_array(total_precipitation_rate_flag,        q_met_st)
        iTmp = fu_merge_integer_to_array(soil_sand_mass_fraction_flag, q_met_st)
        iTmp = fu_merge_integer_to_array(soil_clay_mass_fraction_flag, q_met_st)
        iTmp = fu_merge_integer_to_array(fraction_of_land_flag, q_met_st)
        iTmp = fu_merge_integer_to_array(leaf_area_index_flag, q_met_st)
        iTmp = fu_merge_integer_to_array(water_capac_soil_srf_grav_flag, q_met_st)
        iTmp = fu_merge_integer_to_array(water_capac_soil_deep_grav_flag, q_met_st)
        if(rulesWIS%evaporation_model == latent_heat_proxy_evap_flag)then
          iTmp = fu_merge_integer_to_array(SILAM_latent_heat_flux_flag, q_met_st)
        endif
      case default
        call set_error('Unknown soil water model:'+fu_str(rulesWIS%soil_water_model),'add_input_needs_WIS')
    end select
  end subroutine add_input_needs_WIS
  
  
  !*****************************************************************************

  subroutine update_water_in_soil(met_buf, disp_buf, now, model_timestep, rulesWIS)
    !
    ! Takes the current state of water in soil fields, precipitation amount during the meteo time step
    ! and computes the next-model-time water amount. Since the rain rate is defined from now till 
    ! now + met_timestep ut filling-up and draining can have shorter time scales, we will stick to 
    ! model time step
    !
    implicit none
    
    ! Imported parameters
    type(Tfield_buffer), INTENT(in) :: met_buf, disp_buf
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: model_timestep
    type(Trules_water_in_soil), intent(in) :: rulesWIS
    
    ! Local variables
    integer, parameter :: indSrf = 1, indDeep = 2
    integer :: iMet, ixMet, iyMet, indWISsrf, indWISdeep
    real :: drain_up2down, evap, plant_uptake, capillar_uplift, dMup_dt, dMdown_dt, seconds
    real, dimension(:), pointer :: pTotalRain, pSWrad, pWISsrf, pWISdeep, pWISvolNWP, pLatentHeat, &
                                 & free_water, pWaterCapacitySrf, pWaterCapacityDeep, pSandMassFract
    !
    ! First, a simple case: if NWP data are used, just turn volumetric to gravimetric
    !
    if(rulesWIS%soil_water_model == soil_water_nwp_flag)then
      !
      ! Target field
      indWISsrf = fu_index(met_buf, water_in_soil_srf_grav_flag, pWISsrf)
      if(fu_fails(indWISsrf /= int_missing, 'Failed surface soil water','update_water_in_soil'))return  ! water in soil surface
      !
      ! Input data
      if(fu_fails(fu_index(met_buf, soil_moisture_vol_frac_nwp_flag, pWISvolNWP) /= int_missing, &     ! water in soil NWP
                            & 'Failed surface soil water','update_water_in_soil'))return
      if(fu_fails(fu_index(met_buf, soil_sand_mass_fraction_flag, pSandMassFract) /= int_missing, &  ! sand mass fraction
                                & 'Failed soil_sand_mass_fraction_flag','compute_emission_for_wb_dust'))return
      ! conversion: gravimentric -> volumetric soil moisture
      !
      pWISsrf(iMet) = pWISvolNWP(iMet) * 1000.0 / &
                    & (2600. * (0.51 + 0.126 * pSandMassFract(iMet)))  ! rho of dry soil 
      ! Update valid time
      !
      call set_valid_time(met_buf%p2d(indWISsrf)%present%idPtr, now)
      return

    elseif(rulesWIS%soil_water_model /= soil_water_nwp_flag)then
      !
      ! Basic dynamic soil water model 
      !
      ! Get the fields: 
      ! - total or cnv+ls precipitation rate; 
      ! - short-wave radiation
      ! - soil water capacity at two levels
      ! - latent heat flux
      ! - 2m temperature 
      ! - two current water_in_soil fields at two depths; 
      !
      if(fu_fails(fu_index(met_buf, total_precipitation_rate_flag, pTotalRain) /= int_missing, &    ! total precip
                              & 'Failed total rain','update_water_in_soil'))return
      if(fu_fails(fu_index(met_buf, surf_sw_down_radiation_flag, pSWrad) /= int_missing, &          ! short-w.rad
                              & 'Failed surface short-wave radiation','update_water_in_soil'))return
      if(fu_fails(fu_index(met_buf, water_capac_soil_srf_grav_flag, pWaterCapacitySrf) /= int_missing, &    ! water cap.srf
                              & 'Failed surface soil water capacity','update_water_in_soil'))return
      if(fu_fails(fu_index(met_buf, water_capac_soil_deep_grav_flag, pWaterCapacityDeep) /= int_missing, &  ! water cap.deep
                              & 'Failed deep soil water capacity','update_water_in_soil'))return
      if(fu_fails(fu_index(met_buf, SILAM_latent_heat_flux_flag, pLatentHeat) /= int_missing, &       ! latent heat flux
                              & 'Failed latent heat flux','update_water_in_soil'))return
      if(fu_fails(fu_index(met_buf, free_soil_water_grav_flag, free_water) /= int_missing, &       ! latent heat flux
                              & 'Failed rfee water in soil','update_water_in_soil'))return
      ! Target fields
      !
      indWISsrf = fu_index(met_buf, water_in_soil_srf_grav_flag, pWISsrf)
      if(fu_fails(indWISsrf /= int_missing, 'Failed surface soil water','update_water_in_soil'))return  ! water in soil surface
      indWISdeep = fu_index(met_buf, water_in_soil_deep_grav_flag, pWISdeep)
      if(fu_fails(indWISdeep /= int_missing, 'Failed deep soil water','update_water_in_soil'))return   ! water in soil deep
    
      seconds = fu_sec(model_timestep)
      !
      ! Cycle over the main grid
      !
      do iyMet = 1, ny_meteo
        do ixMet = 1, nx_meteo
  
          iMet = ixMet + (iyMet-1) * nx_meteo
          !
          ! Rates of changes in the layers. Note that the upper layer is almost independent from the 
          ! lower one except for tiny input through water capillar uplift when the upper layer gets dry.
          !
          ! Upper layer: +rain -drainage_from_upper_to_lower -plant_uptake_up -evaporation
          !
          drain_up2down = max(0., (pWISsrf(iMet) - pWaterCapacitySrf(iMet)) * &
                                & (1. - exp(- seconds / rulesWIS%tauDrain_up2down)))

          plant_uptake = rulesWIS%evapPlant_scale * SWRad_2_PAR * pSWrad(iMet)
          !
          ! Evaporation can be taken from several points of view.
          ! 1. Latent heat is total water release as seen by meteo model. Can be used like this:
          !
          if(plant_uptake > pLatentHeat(iMet) / vaporization_latentheat)then
            evap = pLatentHeat(iMet) / vaporization_latentheat - plant_uptake  ! W/m2 / J/kg = kg/m2sec <-> mm/sec
          else
            evap = 0.0
          endif
          !
          ! 2. Can be taken simply via dry deposition resistance assuming the water content is soil via capacity
          !    Pretty much the fugacity approach. Too many unknown constants so far.
          !
!        evap = (pWISup_past(iDisp) / pWaterFugacity(iMeteo) - q(iMeteo)) / (Ra(iMeteo) + Rb(iMeteo))  ! Rs=0 ??

          !
          ! Capillar uplist just needs water diffusivity in specific soil. So far, put zero
          !
          capillar_uplift = 0.
          !
          ! The water-in-soil change rate
          !
          dMup_dt = pTotalRain(iMet) - drain_up2down - rulesWIS%fract_roots_up * plant_uptake - evap + capillar_uplift

          if(dMup_dt < 0.)then
            !
            ! Discharge goes exponentially: processes taking water out will slow with reduced availability
            !
            pWISsrf(iMet) = pWISsrf(iMet) * exp(dMup_dt * seconds)
          else
            !
            ! The layer gets more water. Just make sure that it does not overfill
            !
            pWISsrf(iMet) = pWISsrf(iMet) + dMup_dt * seconds
        
            if(pWISsrf(iMet) > pWaterCapacitySrf(iMet))then
              drain_up2down = drain_up2down + pWISsrf(iMet) - pWaterCapacitySrf(iMet)
              pWISsrf(iMet) = pWaterCapacitySrf(iMet)                         ! fully charged
            endif
        
          endif  ! if upper layer gets filled or depleted
          !
          ! Lower layer: +drainage_from_upper_to_lower -drainage_from_lower_layer -plant_uptake_down -capillar uplift
          !
          dMDown_dt = drain_up2down - (1.-rulesWIS%fract_roots_up) * plant_uptake - capillar_uplift
      
          if(dMDown_dt < 0.)then
            !
            ! Discharge goes exponentially: processes taking water out will slow with reduced availability
            !
            pWISdeep(iMet) = pWISdeep(iMet) * exp(dMdown_dt * seconds)
          else
            !
            ! The layer gets water. Just get rid of excess: rivers will take it
            !
            pWISdeep(iMet) = pWISdeep(iMet) + dMdown_dt * seconds
        
            if(pWISdeep(iMet) > pWaterCapacityDeep(iMet))then
              free_water(iMet) = free_water(iMet) + (pWISdeep(iMet) - pWaterCapacityDeep(iMet))
              pWISdeep(iMet) = pWaterCapacityDeep(iMet)
            endif
        
          endif ! if lower layer gets filled or depleted
        end do ! x meteo grid
      end do ! y meteo grid
      !
      ! Update valid time
      !
      call set_valid_time(met_buf%p2d(indWISsrf)%present%idPtr, now)
      call set_valid_time(met_buf%p2d(indWISdeep)%present%idPtr, now)

    else
      !
      ! No such soil water model
      !
      call set_error('Unsupported soil water model:'+fu_str(rulesWIS%soil_water_model),'compute_emission_for_wb_dust')
      return

    endif  ! type of soil water model

  end subroutine update_water_in_soil


  !*****************************************************************************

  logical function defined_rulesWIS(rulesWIS)
    implicit none
    type(Trules_water_in_soil), intent(in) :: rulesWIS
    defined_rulesWIS = fu_true(rulesWIS%defined)
  end function defined_rulesWIS
                                
end module water_in_soil
