module fire_danger_indices
  !
  ! This module contains a few classes describing the fire danger indices
  !
  ! All units: SI, unless otherwise stated
  !
  ! Author Mikhail Sofiev mikhail.sofiev@fmi.fi
  !
  ! Language: ANSI FORTRAN-90 (or close to it)
  !
!  use globals
  use field_buffer

  implicit none
  private

  ! public functions of this module
  public set_rules_FDI
  public init_fire_danger_indices_FDI
  public add_input_needs_FDI
  public compute_fire_danger_indices
  public defined

  ! Private functions of this module
  private fu_FDI_4_quantity
  private init_KBDI_idx              ! KBDI
  private add_input_needs_KBDI
  private compute_FDI_KBDI_daily
  private get_parameters_KBDI
  private init_SDI_idx               ! SDI
  private add_input_needs_SDI
  private compute_FDI_SDI_daily
  private get_parameters_SDI
  private init_grass_idx             ! grass
  private add_input_needs_grass
  private compute_FDI_grass_daily
  private get_parameters_grass
  private init_fire_weather_idx      ! fire weather FWI
  private add_input_needs_fire_weather
  private compute_FDI_FWI_daily
  private get_parameters_FWI
  private init_fuel_moisture_idx     ! fuel moisture, not implemented
  private add_input_needs_fuel_moisture
!  private compute_FDI_fuel_moisture
  private get_parameters_fuel_moisture
  private rules_FDI_defined
  
  interface defined
    module procedure rules_FDI_defined
  end interface
  

  !
  ! Fire danger indices defined in this module
  !
  integer, private, parameter :: fdi_FireWeather = 6701
  integer, private, parameter :: fdi_Grass_mean = 6702
  integer, private, parameter :: fdi_Grass_max = 6703
  integer, private, parameter :: fdi_KBDI = 6704
  integer, private, parameter :: fdi_SDI = 6705
  integer, private, parameter :: fdi_FuelMoisture = 6706
  !
  ! Definitions of the FDIs
  ! The actual coefficients are given at the definition moment as parameters of the rules
  ! Each index however can be redefined from a given text file written in a form of a namelist
  ! 
  !---------------------------------------------------------------------------
  type Tdrought_factor
    private
    ! Drought factor
    ! Somewhat special class: it needs soil moisture deficit as an input, either from KBDI or SDI
    ! So, the actual computation is just a function that would return the drought_factor in response
    ! to a bunch of input variables
    !
    real :: y_rain = 3.0           ! 3.0 
    real :: dd_min = 0.8          ! 0.8
    real :: dd_power = 1.3           ! 1.3
    real :: sdi_a = 30           ! 30.0
    real :: sdi_b = 40          ! 40.0
    real :: sdi_c = 11          ! 10.5
    real :: y_a = 42          ! 42.0
    real :: y_b = 3           ! 3.0
    real :: y_c = 42          ! 42.0
    real :: sdi_max = 10      ! 10.0
    character(len=23) :: outName = 'drought_factor'
  end type Tdrought_factor
  
  !---------------------------------------------------------------------------
  type TKBDI
    private
    ! KBDI
    real :: PrIntercept = 5.0  ! 5mm of intercept
    real :: KBDI_max = 200     ! 203.2  # 203.2 mm of full-dry soil: 8 inches
    real :: Scale1 = 1.0       ! 0.968
    real :: exp1_a = 0.09      ! 0.0875
    real :: exp1_b = 1.56      ! 1.5552
    real :: c = 8.3            ! 8.30
    real :: Scale2 = 11.0      ! 10.88
    real :: exp2_a = 0.002     ! 0.00173
    type(Tdrought_factor) :: drought_factor
    character(len=23) :: outName = 'soil_moist_deficit_KBDI'
    type(silja_field), pointer :: pFld_moistDef => NULL(), pFld_DroughtFactor => NULL(), &
!                                & pFld_rainToday => NULL(), pFld_temprMax => NULL(), &
                                & pFld_rainYesterday => NULL()
  end type TKBDI

  !---------------------------------------------------------------------------
  type TSDI
    private
    real :: rain_a = 0.2
    real :: rain_b = 0.025
    real :: rain_max = 2.0
    real :: moist_a = 25.0
    real :: moist_b = 0.02
    real :: moist_scale = 6.0
    character(len=22) :: outName = 'soil_moist_deficit_SDI'
    type(silja_field), pointer :: pFld_SDI => NULL()
  end type TSDI

  !---------------------------------------------------------------------------
  type TGrass
    private
    real :: base = 10.0
    real :: intercept = 0.009   ! 0.009254
    real :: T_scale = 0.012     ! 0.01201
    real :: wind_scale = 0.28   ! 0.2789
    real :: cure_power = 1.5    ! 1.536
    real :: cure_scale = 0.004  ! 0.004096
    real :: RH_scale = 0.1      ! 0.09577
    real :: scale = 1.0         ! 1.021537
    integer :: FDI_type = FDI_Grass_max   ! mean or max
    character(len=22) :: outName = 'grass_fire_danger_mean'
    integer :: count_times = 0
    type(silja_field), pointer :: pFld_GFD => NULL(), pFld_temprMM => NULL(), pFld_windMM => NULL(), &
                                & pFld_RelHumMM => NULL()
  end type TGrass

  !---------------------------------------------------------------------------
  type TFWI_fine_fuel_moisture
    private
    real :: e1_1 = 150     ! 147.2
    real :: e1_2 = 100     ! 101
    real :: e1_3 = 60      ! 59.5
    real :: e3_1 = 42.5
    real :: e3_2 = 251     ! 251
    real :: e3_3 = 6.9     ! 6.93
    real :: e3_4 = 1.5e-3  ! 1.5e-3
    real :: e3_5 = 150     ! 150
    real :: e4_1 = 0.9     ! 0.941
    real :: e4_2 = 0.68    ! 0.679
    real :: e4_3 = 11      ! 11
    real :: e4_4 = 0.2     ! 0.18
    real :: e4_5 = 21      ! 21.1
    real :: e4_6 = 0.12    ! 0.115
    real :: e5_1 = 0.6     ! 0.118
    real :: e5_2 = 0.75    ! 0.753
    real :: e5_3 = 10      ! 10
    real :: ek_p1 = 1.7
    real :: ek_p2 = 8.
    real :: ek_1 = 0.6     ! 0.581
    real :: ek_2 = 0.04    ! 0.0365
    real :: ek_3 = 0.4     ! 0.424
    real :: ek_4 = 0.07    ! 0.0694
    real :: e10_1 = 60.    ! 59.5
    real :: e10_2 = 250.
    real :: e10_3 = 150.   ! 147.2
    character(len=19) :: outName = 'FWI_fine_fuel_moist'
    type(silja_field), pointer :: pFld_FFMoist => NULL()
  end type TFWI_fine_fuel_moisture

  !--------------------------------------------------------------
  type TFWI_duff_moisture
    real :: e11_1 = 0.9   ! 0.92
    real :: e11_2 = 1.3   ! 1.27
    real :: e12_0 = 20.
    real :: e12_1 = 5.6   ! 5.6348
    real :: e12_2 = 43.   ! 43.43
    real :: e13_1 = 33.
    real :: e13_2 = 0.5
    real :: e13_3 = 0.3
    real :: e13_4 = 65.
    real :: e13_5 = 14.
    real :: e13_6 = 1.3
    real :: e13_7 = 6.2
    real :: e13_8 = 17.2
    real :: e14 = 49      ! 48.77
    real :: e15_1 = 245   ! 244.72
    real :: e16_1 = 1.9   ! 1.894
    real :: e16_2 = 1.1
    real :: e16_3 = 3.
    character(len=14) :: outName = 'FWI_duff_moist'
    type(silja_field), pointer :: pFld_DuffMoist => NULL()
  end type TFWI_duff_moisture

  !--------------------------------------------------------------
  type TFWI_drought
    real :: e18_1 = 0.83
    real :: e18_2 = 1.3    ! 1.27
    real :: e19_1 = 800. 
    real :: e19_2 = 400.
    real :: e20_1 = 3.9    ! 3.937
    real :: e22_1 = 0.36
    real :: e22_2 = 2.8
    real :: lf_cutoff = -1.6
    real :: lf_corr = 12.
    character(len=11) :: outName = 'FWI_drought'
    type(silja_field), pointer :: pFld_drought => NULL()
  end type TFWI_drought

  !---------------------------------------------------------------------------
  type TFireWeather
    private
    type(TFWI_fine_fuel_moisture) :: FFM
    type(TFWI_duff_moisture) :: DM
    type(TFWI_drought) :: Drght
    real :: e24 = 0.05    ! 0.05039
    real :: e25_1 = 92    ! 91.9
    real :: e25_2 = 0.14  ! 0.1386
    real :: e25_3 = 5.3   ! 5.31
    real :: e25_4 = 1./4.9e7   ! 4.93e7
    real :: e26 = 0.2     ! 0.208
    real :: e27_1 = 0.8
    real :: e27_2 = 0.4
    real :: e27_3 = 0.9   ! 0.92
    real :: e27_4 = 0.01  ! 0.0114
    real :: e27_5 = 1.7
    real :: e28_1 = 80
    real :: e28_2 = 0.6   ! 0.626 
    real :: e28_3 = 0.8   ! 0.809
    real :: e28_4 = 2.
    real :: e28_5 = 1000.
    real :: e28_6 = 25.
    real :: e28_7 = 109   ! 108.65
    real :: e28_8 = 0.02  ! 0.023
    real :: e29 = 0.1
    real :: e30_1 = 2.7   ! 2.72 
    real :: e30_2 = 0.4   ! 0.34
    real :: e30_3 = 0.65  ! 0.647
    character(len=16) :: outName = 'FWI_fire_weather'
    integer :: count_times = 0
    type(silja_field), pointer :: pFld_tempr => NULL(), pFld_RH => NULL(), &
                                & pFld_rainToday => NULL(), pFld_wind => NULL(), &
                                & pFld_FireWeather => NULL()
  end type TFireWeather

  !---------------------------------------------------------------------------
  type TFuelMoisture
    private
    character(len=22) :: outName = 'fuel_moisture'
    type(silja_field), pointer :: pFld_FM => NULL()
  end type TFuelMoisture
  
  !
  ! The main type of the module: a collection of the fire danger indices
  ! Since they can contain some internal fields and variables, only those that are needed
  ! are activated
  ! Also serves as rules for fire danger indices
  !
  type Trules_fire_danger
    private
    integer :: nFDI
    integer, dimension(10) :: FDI_type
    type(TKBDI), allocatable :: FDI_KBDI
    type(TSDI), allocatable :: FDI_SDI
    type(TGrass), allocatable :: FDI_Grass_mean, FDI_Grass_max
    type(TFireWeather), allocatable :: FDI_FireWeather
    type(TFuelMoisture), allocatable :: FDI_Fuel_Moisture
    type(silja_logical) :: defined
  end type Trules_fire_danger
  public TRules_fire_danger
  
  type(Trules_fire_danger), parameter :: rulesFDI_missing = &
           & Trules_fire_danger(0, int_missing, NULL(), NULL(), NULL(), NULL(), &
                              & NULL(), NULL(), silja_false)
  public rulesFDI_missing


  CONTAINS

  
  !**************************************************************************************
  
  subroutine set_rules_FDI(nlFDI_ini, rules_FDI, out_quantities_st)
    !
    ! A wrapper for the initialization of a bunch of FDIs
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist), pointer :: nlFDI_ini
    type(Trules_fire_danger), intent(out) :: rules_FDI
    integer, dimension(:), intent(in) :: out_quantities_st
    
    ! local variables
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: pItems
    type(Tsilam_namelist_group), pointer :: nlFDI_grp
    type(Tsilam_namelist), pointer :: nlFDI
    integer :: iTmp, iFDI, iUnit
    character(len=clen) :: chFDI
    character(len=fnlen) :: chFNm, sp
    real, dimension(:), pointer :: pValues
    logical :: ifOK
    
    ! Get the indices from the control-file namelist
    !
    if(.not. associated(nlFDI_ini))then
      rules_FDI =  rulesFDI_missing
      return
    endif
    nullify(pItems)
    call get_items(nlFDI_ini, 'fire_danger_index', pItems, rules_FDI%nFDI)
    if(error)return
    !
    ! Set indices one by one
    !
    do iFDI = 1, rules_FDI%nFDI
      ! Each FDI is a name of the index and, may be, file name with parameters
      sp = fu_content(pItems(iFDI))
      read(unit=sp, iostat=iTmp, fmt=*) chFDI, chFNm  ! FDI type and coefficients file
      !
      ! Can be default coefficients - then the literature hard-coded values
      if(trim(fu_str_u_case(chFNm)) == 'DEFAULT')then
        nullify(nlFDI_grp)
        nullify(nlFDI)
      else
        iUnit = fu_next_free_unit()
        !! FIXME Someone should close it at some point
        OPEN(file = chFNm, unit = iUnit, action = 'read', status = 'old', iostat = iTmp)
        if(fu_fails(iTmp == 0, 'cannot open input file:' + chFNm, &
                             & 'set_rules_FDI'))return
        nlFDI_grp => fu_read_namelist_group(iUnit, .true.)
        if(fu_nbr_of_namelists(nlFDI_grp) == 1)then
          nlFDI => fu_namelist(nlFDI_grp, 1)
        elseif(fu_nbr_of_namelists(nlFDI_grp) == 0)then
          nullify(nlFDI)
        else
          call set_error('Strange number of namelists in:' + chFNm &
                     & + ', n=' + fu_str(fu_nbr_of_namelists(nlFDI_grp)), &
                     & 'set_rules_FDI')
          return
        endif
      endif
      !
      ! Initialise the one that is pointed out by the name
      !
      if(fu_str_u_case(chFDI) == fu_str_u_case('soil_moist_deficit_KBDI'))then   ! KBDI
        rules_FDI%FDI_type(iFDI) = fdi_KBDI
        allocate(rules_FDI%FDI_KBDI)
        call get_parameters_KBDI(rules_FDI%FDI_KBDI, nlFDI)

      elseif(fu_str_u_case(chFDI) == fu_str_u_case('soil_moist_deficit_SDI'))then    ! SDI
        rules_FDI%FDI_type(iFDI) = fdi_SDI
        allocate(rules_FDI%FDI_SDI)
        call get_parameters_SDI(rules_FDI%FDI_SDI, nlFDI)

      elseif(fu_str_u_case(chFDI) == fu_str_u_case('grass_fire_danger_mean'))then    ! GRASS-mean
        rules_FDI%FDI_type(iFDI) = fdi_Grass_mean
        allocate(rules_FDI%FDI_Grass_mean)
        call get_parameters_grass(rules_FDI%FDI_Grass_mean, nlFDI)
        
      elseif(fu_str_u_case(chFDI) == fu_str_u_case('grass_fire_danger_max'))then     ! GRASS-max
        rules_FDI%FDI_type(iFDI) = fdi_Grass_max
        allocate(rules_FDI%FDI_Grass_max)
        call get_parameters_grass(rules_FDI%FDI_Grass_max, nlFDI)

      elseif(fu_str_u_case(chFDI) == fu_str_u_case('fire_weather'))then                  ! FIRE WEATHER
        rules_FDI%FDI_type(iFDI) = fdi_FireWeather
        allocate(rules_FDI%FDI_FireWeather)
        call get_parameters_FWI(rules_FDI%FDI_FireWeather, nlFDI_grp)

      elseif(fu_str_u_case(chFDI) == fu_str_u_case('fuel_moisture'))then             ! FUEL_MOISTURE
        rules_FDI%FDI_type(iFDI) = fdi_FuelMoisture
        allocate(rules_FDI%FDI_Fuel_Moisture)
        call get_parameters_Fuel_Moisture(rules_FDI%FDI_Fuel_Moisture, nlFDI)

      else
        call set_error('Unknown FDI:' + chFDI,'set_rules_FDI')
        return
      endif  ! FDI name

      ! Prepare to the next round
      call reset_namelist_group(nlFDI_grp)
    end do  ! rules_FDI%nFDI
    
    ! Stupidity check: output must not require more than activated
    ! Opposite is not an issue: indices might be needed for IS4FIRES runs
    !
    do iTmp = 1, size(out_quantities_st)
      if(out_quantities_st(iTmp) == int_missing)exit   ! all fine
      iUnit = fu_FDI_4_quantity(out_quantities_st(iTmp))
      if(iUnit > 0)then
        ifOK = .false.
        do iFDI = 1, rules_FDI%nFDI
          ifOK = iUnit == rules_FDI%FDI_type(iFDI)  ! this FDI produces this quantity
          if(ifOK)exit
        end do
        if(.not. ifOK)then
          call set_error('Fire danger quantity ' + fu_quantity_string(out_quantities_st(iTmp)) &
                       & + '- in output but not in IS4FIRES namelist','set_rules_FDI')
          return
        endif
      endif   ! FDI quantity
    end do  ! output quantities

    if(associated(nlFDI_grp)) deallocate(nlFDI_grp)
    
    rules_FDI%defined = silja_true

  end subroutine set_rules_FDI
  
  
  !**************************************************************************************
  
  subroutine init_fire_danger_indices_FDI(rules_FDI, dispersionMarketPtr, start_time)
    !
    ! A wrapper for the initialization of a bunch of FDIs
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist), pointer :: nlFDI_ini
    type(Trules_fire_danger), intent(inout) :: rules_FDI
    type(mini_market_of_stacks), pointer :: dispersionMarketPtr
    type(silja_time), intent(in) :: start_time
    
    ! local variables
    integer :: iTmp, iFDI, iUnit
    logical :: ifOK
    !
    ! Set indices one by one
    !
    do iFDI = 1, rules_FDI%nFDI
      !
      ! Initialise the one that is pointed out by the name
      !
      select case(rules_FDI%FDI_type(iFDI))
        case(fdi_KBDI)
          call init_KBDI_idx(rules_FDI%FDI_KBDI, dispersionMarketPtr, start_time)

        case(fdi_SDI)
          call init_SDI_idx(rules_FDI%FDI_SDI, dispersionMarketPtr, start_time)

        case(fdi_Grass_mean)
          call init_grass_idx(rules_FDI%FDI_Grass_mean, dispersionMarketPtr, start_time, &
                            & fdi_Grass_mean)
        
        case(fdi_Grass_max)
          call init_grass_idx(rules_FDI%FDI_Grass_max, dispersionMarketPtr, start_time, &
                            & fdi_Grass_max)

        case(fdi_FireWeather)
          call init_fire_weather_idx(rules_FDI%FDI_FireWeather, dispersionMarketPtr, start_time)

        case(fdi_FuelMoisture)
          call init_fuel_moisture_idx(rules_FDI%FDI_fuel_moisture, dispersionMarketPtr, start_time)
        
        case default
          call set_error('Unknown FDI type:' + fu_str(rules_FDI%FDI_type(iFDI)),'init_fire_danger_indices_FDI')
          return
      end select  ! FDI type

    end do  ! rules_FDI%nFDI
    
  end subroutine init_fire_danger_indices_FDI


  !************************************************************************************
  
  subroutine add_input_needs_FDI(rules_FDI, q_met_dyn, q_met_st, q_disp_dyn, q_disp_st)
    !
    ! Returns input needs for the fire danger indices.
    ! Note that the indices themselves are not yet initialized, rules cannot be used
    !
    implicit none

    ! Imported parameters
    type(Trules_fire_danger), intent(in) :: rules_FDI
    integer, dimension(:), intent(inout) :: q_met_dyn, q_met_st, &
                                          & q_disp_dyn, q_disp_st
    ! Local variables
    integer :: iFDI
    !
    ! Add needed quantities depending on what indices are activated
    !
    do iFDI = 1, rules_FDI%nFDI
      select case(rules_FDI%FDI_type(iFDI))
        case(fdi_KBDI)
          call add_input_needs_KBDI(q_met_dyn, q_met_st, q_disp_dyn, q_disp_st)
        case(fdi_SDI)
          call add_input_needs_SDI(q_met_dyn, q_met_st, q_disp_dyn, q_disp_st)
        case(fdi_grass_mean, fdi_grass_max)
          call add_input_needs_grass(rules_FDI%FDI_type(iFDI), &
                                   & q_met_dyn, q_met_st, q_disp_dyn, q_disp_st)
        case(fdi_FireWeather)
          call add_input_needs_fire_weather(q_met_dyn, q_met_st, q_disp_dyn, q_disp_st)
        case(fdi_FuelMoisture)
          call add_input_needs_fuel_moisture(q_met_dyn, q_met_st, q_disp_dyn, q_disp_st)
        case default
          call set_error('Unknown fire danger type:' + fu_str(rules_FDI%FDI_type(iFDI)), &
                       & 'add_input_needs_FDI')
          return
      end select
    end do
  end subroutine add_input_needs_FDI

  
  !************************************************************************************
  
  subroutine compute_fire_danger_indices(rules_FDI, met_buf, now, &
!                                       & ifHorizInterp, pHorizInterpMet2DispStruct, &
                                       & dispMarketPtr)
    !
    ! This subroutine computes all requested FDIs
    ! It is a part of diagnostic_variables procedure, i.e. dispersion mini-market is ready
    ! but disp_buffer is not
    !
    ! Imported parameters
    type(Trules_fire_danger), target, intent(inout) :: rules_FDI
    type(Tfield_buffer), pointer ::  met_buf                ! meteo field buffer
    type(mini_market_of_stacks), pointer :: dispMarketPtr   ! mini-market to store data to
    type(silja_time), intent(in) :: now  ! time to diagnose for
!    type(THorizInterpStruct), pointer :: pHorizInterpMet2DispStruct  ! interpolation structure
!    logical, intent(in) :: ifHorizInterp                    ! whether interpolation is needed

    ! Local variables
    integer :: idx
    type(THorizInterpStruct), pointer :: pHorizInterpMet2DispStruct  ! interpolation structure
    logical :: ifHorizInterp                    ! whether interpolation is needed

    ifHorizInterp = .not. meteo_grid == dispersion_grid
    if(ifHorizInterp) &
          & pHorizInterpMet2DispStruct => fu_horiz_interp_struct(meteo_grid, dispersion_grid, &
                                                               & linear, .true.)
    if(error)return
    !
    ! Scan the whole set of defined FDIs in the rules
    ! Note that some indices are daily and some are hourly. Call appropriately
    ! Daily index can be computed at the end of the day and have a validity period from that day
    ! start till its end. The period does not cover the next day: that day can be rainy, for
    ! instance
    !
    do idx = 1, rules_FDI%nFDI
      select case(rules_FDI%FDI_type(idx))
        case(fdi_KBDI)
          if(now == fu_start_of_day_utc(now)) &
                  & call compute_FDI_KBDI_daily(rules_FDI%FDI_KBDI, met_buf, dispMarketPtr, &
                                              & pHorizInterpMet2DispStruct, ifHorizInterp, now)
        case(fdi_SDI)
          if(now == fu_start_of_day_utc(now)) &
                  & call compute_FDI_SDI_daily(rules_FDI%FDI_SDI, dispMarketPtr, &
                                             & pHorizInterpMet2DispStruct, ifHorizInterp, now)
        case(fdi_Grass_mean)
          if(now == fu_start_of_day_utc(now)) &
                  & call compute_FDI_grass_daily(rules_FDI%FDI_grass_mean, met_buf, dispMarketPtr, &
                                               & pHorizInterpMet2DispStruct, ifHorizInterp, now)
        case(fdi_Grass_max)
          if(now == fu_start_of_day_utc(now)) &
                  & call compute_FDI_grass_daily(rules_FDI%FDI_grass_max, met_buf, dispMarketPtr, &
                                               & pHorizInterpMet2DispStruct, ifHorizInterp, now)
        case(fdi_FireWeather)
          if(now == fu_start_of_day_utc(now)) &
                  & call compute_FDI_FWI_daily(rules_FDI%FDI_FireWeather, met_buf, dispMarketPtr, &
                                             & pHorizInterpMet2DispStruct, ifHorizInterp, now)
        case(fdi_FuelMoisture)
          call set_error('fuel moisture FDI is not implemented','compute_fire_danger_indices')
        case default
          call set_error('Unknown fire danger index:' + fu_str(idx),'compute_fire_danger_indices')
          return
      end select
    end do  ! rules_FDI%nFDI
  end subroutine compute_fire_danger_indices

  
  !**************************************************************************************

  integer function fu_FDI_4_quantity(quantity)
    !
    ! Checks whether the given quantity belongs here
    !
    implicit none
    
    ! Imported parameter
    integer, intent(in) :: quantity

    select case(quantity)
      case(FDI_KBDI_moisture_deficit_flag, FDI_KBDI_drought_factor_flag)
        fu_FDI_4_quantity = fdi_KBDI
      case(FDI_SDI_fire_danger_flag)
        fu_FDI_4_quantity = fdi_SDI
      case(FDI_grass_mean_fire_danger_flag)
        fu_FDI_4_quantity = fdi_Grass_mean
      case(FDI_grass_max_fire_danger_flag)
        fu_FDI_4_quantity = fdi_Grass_max
      case(FDI_FWI_fine_fuel_moist_flag, &
         & FDI_FWI_duff_moist_flag, FDI_FWI_drought_flag, &
         & FDI_fire_weather_index_flag, FDI_fuel_moisture_flag)
        fu_FDI_4_quantity = fdi_FireWeather
      case default
        fu_FDI_4_quantity = int_missing
    end select
    ! unused: fdi_FuelMoisture

  end function fu_FDI_4_quantity
  
  
  
  !**************************************************************************************
  !**************************************************************************************
  !
  ! Drought factor, based on either KBDI or SDI
  ! Keeth J.J. & Byram G.M., (1968) A drought factor index for forest fire control;. 
  !                              USDA Forest Service Research Paper SE-38, 32 pp.
  ! Soil Dryness Index (SDI) Mount (1972) The derivation and testing of a soil dryness index using
  !                              run-off data. Tasmanian Forestry Commission Bulletin No 4, 31 pp.
  !
  ! Formulas for DF can be found here:
  ! Finkele, K., Mills, G. A., Beard, G., and Jones, D. A.: National gridded drought factors 
  !                              and comparison of two soil moisture deficit formulations used
  !                              in prediction of Forest Fire Danger Index in Australia, Australian
  !                              Meteorological Magazine, 15, 2006.
  !
  ! Note that both indices remember their previous-day state
  !
  !**************************************************************************************
  !**************************************************************************************

  subroutine init_DF(self)
    !
    !
    implicit none
    
    ! Imported parameters
    type(Tdrought_factor), intent(inout) :: self
    
    !
    ! DF does not have variables and rather provides a function for computing drought factor
    ! if soil moisture deficit is known. So, the corresponding fields are included in KBDI and SDI
    !
  end subroutine init_DF

  
  !================================================================================

  subroutine add_input_needs_DF(q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st)
    !
    ! Returns input needs for the transformation routines. 
    !
    implicit none

    ! Imported parameters
    integer, dimension(:), intent(inout) :: q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st

    ! Local variables
    integer :: iTmp

    ! Add needed quantities Meteo
    iTmp = fu_merge_integer_to_array(day_sum_precipitation_flag, q_disp_st)
    iTmp = fu_merge_integer_to_array(dry_days_count_flag,        q_disp_st)

  end subroutine add_input_needs_DF
  

  !======================================================================================
  
  subroutine get_parameters_DF(self, nlFDI)
    !
    ! Reads the parameters of the KBDI index
    !
    implicit none
    
    ! Input parameter:
    type(Tdrought_factor), intent(out) :: self
    type(Tsilam_namelist), pointer :: nlFDI
    
    self%y_rain = fu_content_real(nlFDI,'y_rain')          ! 3.0 
    self%dd_min = fu_content_real(nlFDI,'dd_min')          ! 0.8
    self%dd_power = fu_content_real(nlFDI,'dd_power')      ! 1.3
    self%sdi_a = fu_content_real(nlFDI,'sdi_a')            ! 30.0
    self%sdi_b = fu_content_real(nlFDI,'sdi_b')            ! 40.0
    self%sdi_c = fu_content_real(nlFDI,'sdi_c')            ! 10.5
    self%y_a = fu_content_real(nlFDI,'y_a')            ! 42.0
    self%y_b = fu_content_real(nlFDI,'y_b')            ! 3.0
    self%y_c = fu_content_real(nlFDI,'y_c')            ! 42.0
    self%sdi_max = fu_content_real(nlFDI,'sdi_max')            ! 10.0

  end subroutine get_parameters_DF

  
  !================================================================================
  
  real function fu_compute_DF(self, moist_def, rain_yesterday, nDryDays) result(DF)
    !
    ! Computes the drought factor if a soil moisture deficit and some meteorology is given
    !
    implicit none

    ! Imported parameters
    type(Tdrought_factor), intent(inout) :: self
    real, intent(in) :: moist_def, rain_yesterday
    real, intent(in) :: nDryDays

    ! Local variables
    real :: y

    y = (rain_yesterday - self%y_rain) / (max(nDryDays, self%dd_min)) ** self%dd_power
    !
    ! exp can be of a large negative value. Take care of it
    if((moist_def + self%sdi_a) / self%sdi_b >= 10)then
      DF = max(min(self%sdi_c * (y + self%y_a) / (y*y + y * self%y_b + self%y_c), self%sdi_max), 0.0)
    else
      DF = max(min((1.0 - exp(-(moist_def + self%sdi_a) / self%sdi_b)) * self%sdi_c * &
                 & ((y + self%y_a) / (y*y + y * self%y_b + self%y_c)), self%sdi_max), 0.0)
    endif
    
  end function fu_compute_DF



  !**************************************************************************************
  !**************************************************************************************
  !
  ! Drought indices: KBDI 
  ! Keeth J.J. & Byram G.M., (1968) A drought factor index for forest fire control;. 
  !                                 USDA Forest Service Research Paper SE-38, 32 pp.
  ! Soil Dryness Index (SDI) Mount (1972) The derivation and testing of a soil dryness index using run-off data. 
  !                                       Tasmanian Forestry Commission Bulletin No 4, 31 pp.
  ! Note that both indices remember their previous-day state
  !
  !**************************************************************************************
  !**************************************************************************************

  subroutine init_KBDI_idx(self, dispersionMarketPtr, start_time)
    !
    !
    implicit none
    
    ! Imported parameters
    type(TKBDI), intent(inout) :: self
    type(mini_market_of_stacks), pointer :: dispersionMarketPtr
    type(silja_time), intent(in) :: start_time
    
    ! local variables
    logical :: ifFound
    integer :: status
    real, dimension(:), pointer :: pValues
    type(silja_field), pointer :: pFld
    
    call init_DF(self%drought_factor)

    !
    ! KBDI consists of two variables. They all require several fields to be created in 
    ! a dispersion buffer: diagnosed variables. They are diagnosed elsewhere, here we only 
    ! use them
    !
    call find_or_create_field(FDI_KBDI_moisture_deficit_flag, dispersionMarketPtr, species_missing, &
                            & start_time, surface_level, ifFound, self%pFld_moistDef, pValues)
    if(error)return
    if(.not. ifFound) call msg('KBDI moisture deficit created')
    pValues(:) = 0.0  ! for the case of cold start

    call find_or_create_field(FDI_KBDI_drought_factor_flag, dispersionMarketPtr, species_missing, &
                            & start_time, surface_level, ifFound, self%pFld_DroughtFactor, pValues)
    if(error)return
    if(.not. ifFound) call msg('KBDI drought factor created')
    pValues(:) = 0.0  ! for the case of cold start

    call find_or_create_field(yesterday_precipitation_flag, dispersionMarketPtr, species_missing, &
                            & start_time, surface_level, ifFound, self%pFld_rainYesterday, pValues)
    if(error)return
    if(.not. ifFound)call msg('yesterday rain created')
    pValues(:) = 0.0  ! for the case of cold start

    call find_or_create_field(dry_days_count_flag, dispersionMarketPtr, species_missing, &
                            & start_time, surface_level, ifFound, pFld, pValues)
    if(error)return
    if(.not. ifFound)call msg('nDryDays created')
    pValues(:) = 0.0  ! for the case of cold start
    
  end subroutine init_KBDI_idx

  
  !================================================================================

  subroutine add_input_needs_KBDI(q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st)
    !
    ! Returns input needs for the transformation routines. 
    !
    implicit none

    ! Imported parameters
    integer, dimension(:), intent(inout) :: q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st

    ! Local variables
    integer :: iTmp

    ! Add needed quantities
    ! Meteo
    iTmp = fu_merge_integer_to_array(day_sum_precipitation_flag,    q_disp_st)
    iTmp = fu_merge_integer_to_array(day_mean_temperature_2m_flag,  q_disp_st)
    iTmp = fu_merge_integer_to_array(day_max_temperature_2m_flag,   q_disp_st)
    iTmp = fu_merge_integer_to_array(dry_days_count_flag,           q_disp_st)
    iTmp = fu_merge_integer_to_array(mean_annual_precipitation_flag,    q_met_st)
    iTmp = fu_merge_integer_to_array(fraction_of_land_flag,             q_met_st)

    call add_input_needs_DF(q_met_dynamic, q_met_st, q_disp_dynamic,q_disp_st)

  end subroutine add_input_needs_KBDI
  

  !======================================================================================
  
  subroutine get_parameters_KBDI(self, nlFDI)
    !
    ! Reads the parameters of the KBDI index
    !
    implicit none
    
    ! Input parameter:
    type(TKBDI), intent(out) :: self
    type(Tsilam_namelist), pointer :: nlFDI
    
    if(.not. associated(nlFDI))return
    
    ! KBDI itself
    self%PrIntercept = fu_content_real(nlFDI,'PrIntercept')  ! 5mm of intercept
    self%KBDI_max = fu_content_real(nlFDI,'KBDI_max')     ! 203.2  # 203.2 mm of full-dry soil: 8 inches
    self%Scale1 = fu_content_real(nlFDI,'Scale1')  ! 0.968
    self%exp1_a = fu_content_real(nlFDI,'exp1_a')      ! 0.0875
    self%exp1_b = fu_content_real(nlFDI,'exp1_b')      ! 1.5552
    self%c = fu_content_real(nlFDI,'c')            ! 8.30
    self%Scale2 = fu_content_real(nlFDI,'Scale2')      ! 10.88
    self%exp2_a = fu_content_real(nlFDI,'exp2_a')     ! 0.00173

    ! And drought factor
    call get_parameters_DF(self%drought_factor, nlFDI)

  end subroutine get_parameters_KBDI

  
  !================================================================================
  
  subroutine compute_FDI_KBDI_daily(self, met_buf, dispMarketPtr, &
                                  & pHorizInterpMet2DispStruct, ifHorizInterp, now)
    !
    ! Computes the fire danger index KBDI
    ! Note that it uses the previous-day status, and the changes occur on a daily level
    ! This function is called once per day when now == midnight.
    ! It benefits from daily diagnostic variables, which are computed gradually.
    ! Its validity period is the next day.
    !
    implicit none

    ! Imported parameters
    type(TKBDI), intent(inout) :: self
    type(Tfield_buffer), pointer ::  met_buf  ! meteo and internal field buffers
    type(mini_market_of_stacks), pointer :: dispMarketPtr   ! working mini-market
    type(THorizInterpStruct), intent(in) :: pHorizInterpMet2DispStruct
    logical, intent(in) :: ifHorizInterp
    type(silja_time), intent(in) :: now

    ! Local variables
    integer :: iMeteo, ix, iy, ii, indAnnualPrecip
    real :: rainEffect, exp2, y, fTmp
    real, dimension(:), pointer :: fLandFractionMet, moistDef, DF
    logical :: ifFound
    type(silja_field), pointer :: pFld_RainYesterday, pFld_RainToday, pFld_TemprMax, pFld_NDryDays
    real, dimension(:), pointer :: pRainYesterday, pRainToday, pTemprMax, pNDryDays
!    type(silja_field_id), pointer :: idPtr

    fLandFractionMet => fu_grid_data(fraction_of_land_fld)
    !
    ! Fields to compute
    !
    moistDef => fu_grid_data(self%pFld_moistDef)
    DF => fu_grid_data(self%pFld_DroughtFactor)
    if(error)return

    ! meteo fields are provided by diagnostic_variables
    ! We are at the midnight, the just-comptued day-sum rain is TODAY (just-ended today, yes)
    ! The rain yesterday is the one that has been saved locally
    !
    pRainYesterday => fu_grid_data(self%pFld_rainYesterday)
    if(fu_fails(.not.error,'Failed local rain yesterday', 'computed_FDI_KBDI_daily'))return
    
    call find_or_create_field(day_sum_precipitation_flag, dispMarketPtr, species_missing, &
                            & now, surface_level, ifFound, pFld_rainToday, pRainToday)
    if(.not. ifFound)then
      call set_error('Cannot find in dispersion market:' + fu_quantity_string(day_sum_precipitation_flag), &
                   & 'compute_FDI_KBDI_daily')
      return
    endif

    call find_or_create_field(day_max_temperature_2m_flag, dispMarketPtr, species_missing, &
                            & now, level_2m_above_ground, ifFound, pFld_temprMax, pTemprMax)
    if(.not. ifFound)then
      call set_error('Cannot find in dispersion market:' + fu_quantity_string(day_max_temperature_2m_flag), &
                   & 'compute_FDI_KBDI_daily')
      return
    endif

    call find_or_create_field(dry_days_count_flag, dispMarketPtr, species_missing, &
                            & now, surface_level, ifFound, pFld_nDryDays, pNDryDays)
    if(.not. ifFound)then
      call set_error('Cannot find in dispersion market:' + fu_quantity_string(dry_days_count_flag), &
                   & 'compute_FDI_KBDI_daily')
      return
    endif

    indAnnualPrecip = fu_index(met_buf, mean_annual_precipitation_flag, .true.)
    if(error)return
    !
    ! Update the main daily field, will be valid through the next day
    !
    do iy = 1, ny_dispersion
      do ix = 1, nx_dispersion
        ii = ix + nx_dispersion * (iy -1)
        iMeteo =  fu_grid_index(nx_meteo, ix, iy, pHorizInterpMet2DispStruct)
        if(fLandFractionMet(iMeteo) < 0.01)cycle  ! speedup
        !
        ! Effective_rain = actual_rain - (interception,runoff), i.e. the first 5mm within 
        ! consequtive rainy days. So, if 2 days ago was a rainy day, yesterday rain goes as is
        ! otherwise subtract 5mm for initial watering of vegetation and runoff
        !
        ! Dry days
        if(pRainYesterday(ii) < 0.01 * self%PrIntercept)then
          pNDryDays(ii) = pNDryDays(ii) + 1
        else
          pNDryDays(ii) = 0
        endif
        ! Effective rain
        rainEffect = max(0., pRainToday(ii) - max(0., self%PrIntercept - pRainYesterday(ii)))
        !
        ! Evapotransporation increases the moisture deficit
        !
        if(self%exp2_a * met_buf%p2d(indAnnualPrecip)%present%ptr(iMeteo) < 10.)then
            exp2 = exp(-self%exp2_a * met_buf%p2d(indAnnualPrecip)%present%ptr(iMeteo))
        else
            exp2 = 0
        endif
        fTmp = moistDef(ii) + &
                       & ((self%KBDI_max - moistDef(ii)) * &
                        & (self%Scale1 * exp(min(10.,self%exp1_a * (pTemprMax(ii) - 273.15) + self%exp1_b)) - &
                        &  self%c) * 1e-3 / (1.0 + self%Scale2 * exp2)) - rainEffect
!          if(fu_fails(fTmp >= 0,'Moisture deficit<0','compute_FDI_KBDI_daily'))then
!            call msg('ix,iy,value:' + fu_str(ix) + ',' + fu_str(iy), fTmp)
!            call msg('moistDef(ii), self%KBDI_max - moistDef(ii)', moistDef(ii), self%KBDI_max - moistDef(ii))
!            call set_error('Non-positive KBDI moisture deficit','compute_FDI_KBDI_daily')
!            return
!          endif
        moistDef(ii) = max(fTmp, 0.)
        !
        ! Drought factor with KBDI soil moisture deficit
        !
        DF(ii) = fu_compute_DF(self%drought_factor, moistDef(ii), pRainYesterday(ii), pNDryDays(ii))
        !
        ! Prepare next 24 hours of rain accumulation
        !
        pRainYesterday(ii) = pRainToday(ii)
      end do  ! ix
    end do  ! iy
    !
    ! Set the parameters of the fields. They now describe _yesterday_ day as a whole.
    ! Validity is from its start till its end, same as accumulation
    ! But here is an ambiguity: valid_time and accumulation_length jointly decide
    ! start of accumulation, which is valid_time - accum_len. This contradicts to the definition
    ! of the mean value, which is valid from the start of accumulation till its end
    ! For now, I removed accumulation_length setup. Hopefully, can live without it. Also, nDryDays
    ! is not cumulative, this-far.
    !
!    call msg('')
!    call msg('Disperion market before the updates of times')
!    call report(dispMarketPtr)
!    call msg('') 
    
    call set_times(fu_id(self%pFld_moistDef), now)
    call set_times(fu_id(self%pFld_DroughtFactor), now)
    call set_times(fu_id(pFld_nDryDays), now)
    
!    idPtr => fu_id(self%pFld_moistDef)
!    do ix = 1,3
!      call set_valid_time(idPtr, now-one_day)
!      call set_analysis_time(idPtr, now-one_day)
!      call set_accumulation_length(idPtr, one_day)
!      call set_validity_length(idPtr, one_day)
!      idPtr => fu_id(self%pFld_DroughtFactor)
!    end do
    
!    call set_valid_time(self%pFld_DroughtFactor, now-one_day)
!    call set_analysis_time(, now-one_day)
!    call set_accumulation_length(fu_id(self%pFld_DroughtFactor), one_day)
!    call set_validity_length(fu_id(self%pFld_DroughtFactor), one_day)

!    call msg('')
!    call msg('Disperion market after the updates of times')
!    call report(dispMarketPtr)
!    call msg('')

    contains 
                                  
    subroutine set_times(idPtr, now)
      ! just to shorten the code
      implicit none
      type(silja_field_id), pointer :: idPtr
      type(silja_time), intent(in) :: now

      call set_valid_time(idPtr, now-one_day)
      call set_analysis_time(idPtr, now-one_day)
      call set_validity_length(idPtr, one_day)
    end subroutine set_times
    
  end subroutine compute_FDI_KBDI_daily



  !**************************************************************************************
  !**************************************************************************************
  !
  ! Drought indices: SDI, well, in principle.
  !
  ! This index is a major problem: Mount himself has made initial version for tabulated vegetaion
  ! with entirely unclear procedure and classes.
  !
  ! Soil Dryness Index (SDI) Mount (1972) The derivation and testing of a soil dryness index using run-off data. 
  !                                       Tasmanian Forestry Commission Bulletin No 4, 31 pp.
  !
  ! What is implemented below is _NOT_ the Mount's SDI but its interpretation here:
  ! Yeo, C. S., Kepert, J. D., and Hicks, R. (2015) Fire danger indices current  limitations and a 
  !                    pathway to better indices, Department of Industry, Innovation and Science,
  !                    Australian government
  ! 
  ! Note that index remembers its previous-day state
  !
  !**************************************************************************************
  !**************************************************************************************

  subroutine init_SDI_idx(self, dispersionMarketPtr, start_time)
    !
    ! Initialising the SDI index
    !
    implicit none
    
    ! Imported parameters
    type(TSDI), intent(inout) :: self
    type(mini_market_of_stacks), pointer :: dispersionMarketPtr
    type(silja_time), intent(in) :: start_time
    
    ! local variables
    logical :: ifFound
    integer :: status
    real, dimension(:), pointer :: pValues
    
    !
    ! SDI consists of one variable. It requires several fields to be created in 
    ! a dispersion buffer: diagnosed variables. They are diagnosed elsewhere, here we only 
    ! use them
    !
    call find_or_create_field(FDI_SDI_fire_danger_flag, dispersionMarketPtr, species_missing, &
                            & start_time, surface_level, ifFound, self%pFld_SDI, pValues)
    if(error)return
    if(.not. ifFound)call msg('SDI fire danger field created')
    pValues(:) = 0.0  ! for the case of cold start

  end subroutine init_SDI_idx

  
  !================================================================================

  subroutine add_input_needs_SDI(q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st)
    !
    ! Returns input needs for the transformation routines. 
    !
    implicit none

    ! Imported parameters
    integer, dimension(:), intent(inout) :: q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st

    ! Local variables
    integer :: iTmp

    ! Add needed quantities
    !
    iTmp = fu_merge_integer_to_array(day_sum_precipitation_flag,    q_disp_st)
    iTmp = fu_merge_integer_to_array(day_max_temperature_2m_flag,   q_disp_st)
    iTmp = fu_merge_integer_to_array(fraction_of_land_flag,             q_met_st)

  end subroutine add_input_needs_SDI

  
  !======================================================================================
  
  subroutine get_parameters_SDI(self, nlFDI)
    !
    ! Reads the parameters of the SDI index
    !
    implicit none
    
    ! Input parameter:
    type(TSDI), intent(out) :: self
    type(Tsilam_namelist), pointer :: nlFDI
    
    if(.not. associated(nlFDI))return

    self%rain_a = fu_content_real(nlFDI,'rain_a')
    self%rain_b = fu_content_real(nlFDI,'rain_b')
    self%rain_max = fu_content_real(nlFDI,'rain_max')
    self%moist_a = fu_content_real(nlFDI,'moist_a')
    self%moist_b = fu_content_real(nlFDI,'moist_b')
    self%moist_scale = fu_content_real(nlFDI,'moist_scale')

  end subroutine get_parameters_SDI
  
  
  !================================================================================
  
  subroutine compute_FDI_SDI_daily(self, dispMarketPtr, &
                                 & pHorizInterpMet2DispStruct, ifHorizInterp, now)
    !
    ! Computes the fire danger index SDI
    ! Note that it uses the previous-day status, and the changes occur on a daily level
    ! This function is called every timestep to accumulate meteorology but actual calcuations
    ! are done exclusively in the case of now == midnight
    !
    implicit none

    ! Imported parameters
    type(TSDI), intent(inout) :: self
    type(mini_market_of_stacks), pointer :: dispMarketPtr   ! working mini-market
    type(THorizInterpStruct), intent(in) :: pHorizInterpMet2DispStruct
    logical, intent(in) :: ifHorizInterp
    type(silja_time), intent(in) :: now

    ! Local variables
    logical :: ifFound
    integer :: iMeteo, ix, iy, ii, n_out_of_range, npatched, nFailed
    type(silja_field), pointer :: pFld_RainToday, pFld_TemprMax
    real, dimension(:), pointer :: fLandFractionMet, SDI, pRainToday, pTemprMax
    real :: fTmp

    fLandFractionMet => fu_grid_data(fraction_of_land_fld)
    !
    ! Field to compute
    !
    SDI => fu_grid_data(self%pFld_SDI)
    if(error)return
    !
    ! meteo fields are provided by diagnostic_variables
    ! We are at the midnight, the just-comptued day-sum rain is TODAY (just-ended today, yes)
    !
    call find_or_create_field(day_sum_precipitation_flag, dispMarketPtr, species_missing, &
                            & now, surface_level, ifFound, pFld_rainToday, pRainToday)
    if(.not. ifFound)then
      call set_error('Cannot find in dispersion market:' + fu_quantity_string(day_sum_precipitation_flag), &
                   & 'compute_FDI_SDI_daily')
      return
    endif

    call find_or_create_field(day_max_temperature_2m_flag, dispMarketPtr, species_missing, &
                            & now, level_2m_above_ground, ifFound, pFld_temprMax, pTemprMax)
    if(.not. ifFound)then
      call set_error('Cannot find in dispersion market:' + fu_quantity_string(day_max_temperature_2m_flag), &
                   & 'compute_FDI_SDI_daily')
      return
    endif
    
    do iy = 1, ny_dispersion
      do ix = 1, nx_dispersion
        ii = ix + nx_dispersion * (iy -1)
        iMeteo =  fu_grid_index(nx_meteo, ix, iy, pHorizInterpMet2DispStruct)
        if(fLandFractionMet(iMeteo) < 0.01)cycle  ! speedup
        !
        ! Actual calculations
        !
        fTmp = max(0., SDI(ii) - (pRainToday(ii) &
                                & - min(self%rain_max, pRainToday(ii) * self%rain_a) &
                                & - self%rain_b * pRainToday(ii))) &
             & + (max(0., pTemprMax(ii) - 273.15) &
                & / (self%moist_scale * 2.**((SDI(ii) - self%moist_a) * self%moist_b)))
        !
        ! Stupidity check
        !
        if(fu_fails(fTmp >= 0,'SDI<0','compute_FDI_SDI'))then
          call msg('ix,iy,value:' + fu_str(ix) + ',' + fu_str(iy), fTmp)
          call set_error('Negative SDI','compute_FDI_SDI')
          return
        endif
        SDI(ii) = fTmp
      end do  ! ix
    end do  ! iy
    
    call check_quantity_range(FDI_SDI_fire_danger_flag, SDI,  fs_dispersion, nx_dispersion, &
                            & .true., .false., n_out_of_range, npatched, nFailed)
    if(fu_fails(nFailed==0,'Failed patching the SDI field','compute_FDI_SDI_daily'))then
      call msg('Failed SDI')
    endif
    

  end subroutine compute_FDI_SDI_daily


  !**************************************************************************************
  !**************************************************************************************
  !
  ! Grass fire danger index
  !
  ! McArthur Grass Meter: Grass Fire Danger
  ! Importantly, for grass the fire danger does not have memory and is applicable at hourly basis
  ! So, we shall make three of them: daily-max, daily-mean and hourly
  ! Mary-Beth Schreck, Paul J. Howerton and Kenneth R.Cook (2010) Adapting Australia's Grassland
  ! Fire Danger Index for the United States' Central Plains
  !
  !**************************************************************************************
  !**************************************************************************************

  subroutine init_grass_idx(self, dispersionMarketPtr, start_time, FDI_type)
    !
    ! Initialize grass fire danger index. Can be based on mean and max/min daily values.
    !
    implicit none
    
    ! Imported parameters
    type(TGrass), intent(inout) :: self
    type(mini_market_of_stacks), pointer :: dispersionMarketPtr
    type(silja_time), intent(in) :: start_time
    integer, intent(in) :: FDI_type
    
    ! local variables
    logical :: ifFound
    integer :: status
    real, dimension(:), pointer :: pValues
    
    self%FDI_type = FDI_type
    !
    ! grass FDI is just one variable. It requires several fields to be created in 
    ! a dispersion buffer: forecasted variables
    !
    select case(self%FDI_type)
      case(FDI_grass_mean)
        call find_or_create_field(FDI_grass_mean_fire_danger_flag, dispersionMarketPtr, &
                                & species_missing, start_time, surface_level, ifFound, &
                                & self%pFld_GFD, pValues)
      case(FDI_grass_max)
        call find_or_create_field(FDI_grass_max_fire_danger_flag, dispersionMarketPtr, &
                                & species_missing, start_time, surface_level, ifFound, &
                                & self%pFld_GFD, pValues)
      case default
        call set_error('Unknown procedure for the grass FDI:' + fu_str(self%FDI_type), &
                     & 'init_grass_idx')
        return
    end select
    if(error)return
    if(.not. ifFound)call msg('grass FDI created for:' + fu_quantity_string(self%FDI_type))
    pValues(:) = 0.0  ! for the case of cold start

  end subroutine init_grass_idx

  
  !======================================================================================
  
  subroutine add_input_needs_grass(self_type, q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st)
    !
    ! Returns input needs for the transformation routines. 
    !
    implicit none

    ! Imported parameters
    integer, intent(in) :: self_type
    integer, dimension(:), intent(inout) :: q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st

    ! Local variables
    integer :: iTmp

    ! Add needed quantities
    ! Meteo
    select case(self_type)
      case(FDI_grass_mean)
        iTmp = fu_merge_integer_to_array(day_mean_temperature_2m_flag, q_disp_st)
        iTmp = fu_merge_integer_to_array(day_mean_windspeed_10m_flag,  q_disp_st)
        iTmp = fu_merge_integer_to_array(day_mean_relat_humid_2m_flag, q_disp_st)
      case(FDI_grass_max)
        iTmp = fu_merge_integer_to_array(day_max_temperature_2m_flag, q_disp_st)
        iTmp = fu_merge_integer_to_array(day_max_windspeed_10m_flag,  q_disp_st)
        iTmp = fu_merge_integer_to_array(day_min_relat_humid_2m_flag, q_disp_st)
      case default
        call set_error('Unknown procedure for the grass FDI:' + fu_str(self_type), &
                     & 'add_input_needs_grass')
        return
    end select
    iTmp = fu_merge_integer_to_array(leaf_area_indexlv_flag, q_met_dynamic)
    iTmp = fu_merge_integer_to_array(fraction_of_land_flag,         q_met_st)

  end subroutine add_input_needs_grass

  
 !======================================================================================
  
  subroutine get_parameters_grass(self, nlFDI)
    !
    ! Reads the parameters of the grass index
    !
    implicit none
    
    ! Input parameters
    type(Tgrass), intent(out) :: self
    type(Tsilam_namelist), pointer :: nlFDI
    
    if(.not. associated(nlFDI))return

    self%base = fu_content_real(nlFDI,'base')
    self%intercept = fu_content_real(nlFDI,'intercept')
    self%T_scale = fu_content_real(nlFDI,'T_scale')
    self%wind_scale = fu_content_real(nlFDI,'wind_scale')
    self%cure_power= fu_content_real(nlFDI,'cure_power')
    self%cure_scale = fu_content_real(nlFDI,'cure_scale')
    self%RH_scale = fu_content_real(nlFDI,'RH_scale')
    self%scale = fu_content_real(nlFDI,'scale')
    self%FDI_type = fu_content_real(nlFDI,'averaging_type')
   
  end subroutine get_parameters_grass
  

  !================================================================================
  
  subroutine compute_FDI_grass_daily(self, met_buf, dispMarketPtr, &
                                   & pHorizInterpMet2DispStruct, ifHorizInterp, now)
    !
    ! Computes the fire danger index grass
    !
    implicit none

    ! Imported parameters
    type(Tgrass), intent(inout) :: self
    type(Tfield_buffer), pointer ::  met_buf  ! meteo and internal field buffers
    type(mini_market_of_stacks), pointer :: dispMarketPtr   ! working mini-market
    type(THorizInterpStruct), intent(in) :: pHorizInterpMet2DispStruct
    logical, intent(in) :: ifHorizInterp
    type(silja_time), intent(in) :: now

    ! Local variables
    logical :: ifFound
    integer :: iMeteo, ix, iy, ii, qT, qW, qQ, indLAI
    real :: timestep_sec, curing, power
    real, dimension(:), pointer :: fLandFractionMet, danger_idx, pTempr, pWind, pRH, pLAI
    type(silja_field), pointer :: pFld_tempr, pFld_RH, pFld_wind, pFld_LAI

    fLandFractionMet => fu_grid_data(fraction_of_land_fld)
    !
    ! Field to compute
    !
    danger_idx => fu_grid_data(self%pFld_GFD)
    if(error)return
    !
    ! meteo fields
    !
    select case(self%FDI_type)
      case(FDI_grass_mean)
        qT = day_mean_temperature_2m_flag
        qW = day_mean_windspeed_10m_flag
        qQ = day_mean_relat_humid_2m_flag
      case(FDI_grass_max)
        qT = day_max_temperature_2m_flag
        qW = day_max_windspeed_10m_flag
        qQ = day_min_relat_humid_2m_flag
      case default
        call set_error('Unknown procedure for the grass FDI flag:' + fu_str(self%FDI_type), &
                     & 'compute_FDI_grass')
        return
    end select
    !
    ! Now find / create the appropriate fields
    !
    call find_or_create_field(qT, dispMarketPtr, species_missing, &
                            & now, level_2m_above_ground, ifFound, pFld_tempr, pTempr)
    if(.not. ifFound)then
      call set_error('Cannot find in dispersion market:' + fu_quantity_string(qT), &
                   & 'compute_FDI_grass_daily')
      return
    endif

    call find_or_create_field(qW, dispMarketPtr, species_missing, &
                            & now, level_10m_above_ground, ifFound, pFld_wind, pWind)
    if(.not. ifFound)then
      call set_error('Cannot find in dispersion market:' + fu_quantity_string(qW), &
                   & 'compute_FDI_grass_daily')
      return
    endif

    call find_or_create_field(qQ, dispMarketPtr, species_missing, &
                            & now, level_2m_above_ground, ifFound, pFld_RH, pRH)
    if(.not. ifFound)then
      call set_error('Cannot find in dispersion market:' + fu_quantity_string(qQ), &
                   & 'compute_FDI_grass_daily')
      return
    endif

    indLAI = fu_index(met_buf, leaf_area_indexlv_flag, .true.)
    if(error)then
      call set_error('Cannot find in meteo market:' + fu_quantity_string(leaf_area_indexlv_flag), &
                   & 'compute_FDI_grass_daily')
      return
    endif
    !
    ! At midnight, update the main daily field
    !
    !  
    do iy = 1, ny_dispersion
      do ix = 1, nx_dispersion
        ii = ix + nx_dispersion * (iy -1)
        iMeteo =  fu_grid_index(nx_meteo, ix, iy, pHorizInterpMet2DispStruct)
        if(fLandFractionMet(iMeteo) < 0.01)cycle  ! speedup
        !
        ! This is a problematic solution for curing. It is, indeed, zero when / where no grass
        ! but it is low at the growing / flowering season and only after that gets high.
        ! An ideal solution would be to multiply this with some bell shape based on phenology.
        ! But problems come with the hemispheres and grass growing seasons in equatorial regions.
        ! So, for now it is like this, then we'll see.
        !
        curing = met_buf%p2d(indLAI)%present%ptr(iMeteo) &
               & / (met_buf%p2d(indLAI)%present%ptr(iMeteo) + 1.0) * 100.   ! %
        !
        ! Actual computations
        !
        power = self%intercept - (100. - curing)**self%cure_power * self%cure_scale + &
              & (pTempr(ii) - 273.15) * self%T_scale + &      ! C
              & sqrt(pWind(ii) * 3.6) * self%wind_scale - &   ! km/hr
              & sqrt(pRH(ii) * 100.0) * self%RH_scale         ! %
        if(power > 10)then
          danger_idx(ii) = 100
        elseif(power < -10)then
          danger_idx(ii) = 0
        else
          danger_idx(ii) = min(100., self%base ** power * self%scale)
        endif
        !
        ! Stupidity check still
        !
        if(fu_fails(danger_idx(ii) >= 0,'Grass danger<0','compute_FDI_grass'))then
          call msg('ix,iy,meteo_type,value' + fu_str(ix) + ',' + fu_str(iy) + &
                 & ',' + fu_str(self%FDI_type), danger_idx(ii))
          call set_error('Negative grass fire danger','compute_FDI_grass')
          return
        endif
      end do  ! ix
    end do  ! iy

  end subroutine compute_FDI_grass_daily

  
  !**************************************************************************************
  !**************************************************************************************
  !
  ! Fire Weather Index
  !
  ! Canadian Forest Fire Weather Index
  ! Vam Wagner, C.E., (1974) Structure of the Canadian Forest Fire Weather Index
  !       Petawawa Forest Experiment Station, Chalk River, Ontario,
  !       Canadian Forestry Service, Dept. of the Environment, Publication 1333, 49 pp
  ! A better reference:
  ! Van Wagner,C.E., Picket,T.L. (1985) Equation and FORTRAN Program for the Canadian
  !       Forest Fire Weather Index System. Canadian Forestry Service, Government of
  !       Canada, Forestry Technical Report 33, Ottawa
  ! Note: English version misses symbol notations while French version is, well, in French
  !
  ! Another reasonable reference:
  ! Kumar, V. and Dharssi, I. (2015) Sources of soil dryness measures and forecasts for fire danger
  !                                  rating, Australian government Bureau of Meteorology
  !
  ! A series of classes describing it:
  ! Level 1
  ! - Fine Fuel Moisture Code
  ! - Duff Moisture Code
  ! - Drought Code
  ! Level 2
  ! - Initial Spread index (Fine Fuel Moisture)
  ! - Adjusted Duff Moisture (Duff Moisture, Drought Code)
  ! Level 3
  ! - Fire Weather Index (Initial Spread Index, Adjuted Duff Moisture)
  !
  !**************************************************************************************
  !**************************************************************************************

  subroutine init_fire_weather_idx(self, dispersionMarketPtr, start_time)
    !
    implicit none
    
    ! Imported parameters
    type(TFireWeather), intent(inout) :: self
    type(mini_market_of_stacks), pointer :: dispersionMarketPtr
    type(silja_time), intent(in) :: start_time
    
    ! local variables
    logical :: ifFound
    integer :: status
    real, dimension(:), pointer :: pValues
    
    !
    ! FWI contains three other indices. They all require several fields to be created in 
    ! a dispersion buffer: diagnosed variables. They are diagnosed elsewhere, here we only 
    ! use them
    
    call find_or_create_field(FDI_FWI_fine_fuel_moist_flag, dispersionMarketPtr, species_missing, &
                            & start_time, surface_level, ifFound, self%FFM%pFld_FFMoist, pValues)
    if(error)return
    if(.not. ifFound)call msg('FWI fine fuel moisture deficit created')
    pValues(:) = 0.0  ! for the case of cold start

    call find_or_create_field(FDI_FWI_duff_moist_flag, dispersionMarketPtr, species_missing, &
                            & start_time, surface_level, ifFound, self%DM%pFld_DuffMoist, pValues)
    if(error)return
    if(.not. ifFound)call msg('FWI duff moisture deficit created')
    pValues(:) = 0.0  ! for the case of cold start

    call find_or_create_field(FDI_FWI_drought_flag, dispersionMarketPtr, species_missing, &
                            & start_time, surface_level, ifFound, self%Drght%pFld_drought, pValues)
    if(error)return
    if(.not. ifFound)call msg('FWI drought created')
    pValues(:) = 0.0  ! for the case of cold start

    call find_or_create_field(FDI_fire_weather_index_flag, dispersionMarketPtr, species_missing, &
                            & start_time, surface_level, ifFound, self%pFld_FireWeather, pValues)
    if(error)return
    if(.not. ifFound)call msg('FWI fire danger index created')
    pValues(:) = 0.0  ! for the case of cold start

  end subroutine init_fire_weather_idx
  
  
  !======================================================================================
  
  subroutine add_input_needs_fire_weather(q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st)
    !
    ! Returns input needs for the transformation routines. 
    !
    implicit none

    ! Imported parameters
    integer, dimension(:), intent(inout) :: q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st

    ! Local variables
    integer :: iTmp

    ! Add needed quantities
    ! Meteo
    iTmp = fu_merge_integer_to_array(day_sum_precipitation_flag,    q_disp_st)
    iTmp = fu_merge_integer_to_array(day_mean_temperature_2m_flag,  q_disp_st)
    iTmp = fu_merge_integer_to_array(day_mean_relat_humid_2m_flag,  q_disp_st)
    iTmp = fu_merge_integer_to_array(fraction_of_land_flag,      q_met_st)

  end subroutine add_input_needs_fire_weather
  

  !======================================================================================
  
  subroutine get_parameters_FWI(self, nlFDI_grp)
    !
    ! Reads parameters of the Fire Weather Index - all its components
    !
    implicit none
    type(TFireWeather), intent(out) :: self
    type(Tsilam_namelist_group), pointer :: nlFDI_grp
    
    ! Local parameters
    type(Tsilam_namelist), pointer:: nlPtr
    
    if(.not. associated(nlFDI_grp))return

    nlPtr => fu_namelist(nlFDI_grp, 'FWI_fine_fuel_moisture')
    self%FFM%e1_1 = fu_content_real(nlPtr,'e1_1')
    self%FFM%e1_2 = fu_content_real(nlPtr,'e1_2')
    self%FFM%e1_3 = fu_content_real(nlPtr,'e1_3')
    self%FFM%e3_1 = fu_content_real(nlPtr,'e3_1')
    self%FFM%e3_2 = fu_content_real(nlPtr,'e3_2')
    self%FFM%e3_3 = fu_content_real(nlPtr,'e3_3')
    self%FFM%e3_4 = fu_content_real(nlPtr,'e3_4')
    self%FFM%e3_5 = fu_content_real(nlPtr,'e3_5')
    self%FFM%e4_1 = fu_content_real(nlPtr,'e4_1')
    self%FFM%e4_2 = fu_content_real(nlPtr,'e4_2')
    self%FFM%e4_3 = fu_content_real(nlPtr,'e4_3')
    self%FFM%e4_4 = fu_content_real(nlPtr,'e4_4')
    self%FFM%e4_5 = fu_content_real(nlPtr,'e4_5')
    self%FFM%e4_6 = fu_content_real(nlPtr,'e4_6')
    self%FFM%e5_1 = fu_content_real(nlPtr,'e5_1')
    self%FFM%e5_2 = fu_content_real(nlPtr,'e5_2')
    self%FFM%e5_3 = fu_content_real(nlPtr,'e5_3')
    self%FFM%ek_p1 = fu_content_real(nlPtr,'ek_p1')
    self%FFM%ek_p2 = fu_content_real(nlPtr,'ek_p2')
    self%FFM%ek_1 = fu_content_real(nlPtr,'ek_1')
    self%FFM%ek_2 = fu_content_real(nlPtr,'ek_2')
    self%FFM%ek_3 = fu_content_real(nlPtr,'ek_3')
    self%FFM%ek_4 = fu_content_real(nlPtr,'ek_4')
    self%FFM%e10_1 = fu_content_real(nlPtr,'e10_1')
    self%FFM%e10_2 = fu_content_real(nlPtr,'e10_2')
    self%FFM%e10_3 = fu_content_real(nlPtr,'e10_3')

  !--------------------------------------------------------------
    nlPtr => fu_namelist(nlFDI_grp, 'FWI_duff_moisture')
    self%DM%e11_1 = fu_content_real(nlPtr,'e11_1')
    self%DM%e11_2 = fu_content_real(nlPtr,'e11_2')
    self%DM%e12_0 = fu_content_real(nlPtr,'e12_0')
    self%DM%e12_1 = fu_content_real(nlPtr,'e12_1')
    self%DM%e12_2 = fu_content_real(nlPtr,'e12_2')
    self%DM%e13_1 = fu_content_real(nlPtr,'e13_1')
    self%DM%e13_2 = fu_content_real(nlPtr,'e13_2')
    self%DM%e13_3 = fu_content_real(nlPtr,'e13_3')
    self%DM%e13_4 = fu_content_real(nlPtr,'e13_4')
    self%DM%e13_5 = fu_content_real(nlPtr,'e13_5')
    self%DM%e13_6 = fu_content_real(nlPtr,'e13_6')
    self%DM%e13_7 = fu_content_real(nlPtr,'e13_7')
    self%DM%e13_8 = fu_content_real(nlPtr,'e13_8')
    self%DM%e14 = fu_content_real(nlPtr,'e14')
    self%DM%e15_1 = fu_content_real(nlPtr,'e15_1')
    self%DM%e16_1 = fu_content_real(nlPtr,'e16_1')
    self%DM%e16_2 = fu_content_real(nlPtr,'e16_2')
    self%DM%e16_3 = fu_content_real(nlPtr,'e16_3')

  !--------------------------------------------------------------
    nlPtr => fu_namelist(nlFDI_grp, 'FWI_drought')
    self%Drght%e18_1 = fu_content_real(nlPtr,'e18_1')
    self%Drght%e18_2 = fu_content_real(nlPtr,'e18_2')
    self%Drght%e19_1 = fu_content_real(nlPtr,'e19_1')
    self%Drght%e19_2 = fu_content_real(nlPtr,'e19_2')
    self%Drght%e20_1 = fu_content_real(nlPtr,'e20_1')
    self%Drght%e22_1 = fu_content_real(nlPtr,'e22_1')
    self%Drght%e22_2 = fu_content_real(nlPtr,'e22_2')
    self%Drght%lf_cutoff = fu_content_real(nlPtr,'lf_cutoff')
    self%Drght%lf_corr = fu_content_real(nlPtr,'lf_corr')

  !---------------------------------------------------------------------------
    nlPtr => fu_namelist(nlFDI_grp, 'FWI_fire?weather')
    self%e24 = fu_content_real(nlPtr,'e24')
    self%e25_1 = fu_content_real(nlPtr,'e25_1')
    self%e25_2 = fu_content_real(nlPtr,'e25_2')
    self%e25_3 = fu_content_real(nlPtr,'e25_3')
    self%e25_4 = fu_content_real(nlPtr,'e25_4')
    self%e26 = fu_content_real(nlPtr,'e26')
    self%e27_1 = fu_content_real(nlPtr,'e27_1')
    self%e27_2 = fu_content_real(nlPtr,'e27_2')
    self%e27_3 = fu_content_real(nlPtr,'e27_3')
    self%e27_4 = fu_content_real(nlPtr,'e27_4')
    self%e27_5 = fu_content_real(nlPtr,'e27_5')
    self%e28_1 = fu_content_real(nlPtr,'e28_1')
    self%e28_2 = fu_content_real(nlPtr,'e28_2')
    self%e28_3 = fu_content_real(nlPtr,'e28_3')
    self%e28_4 = fu_content_real(nlPtr,'e28_4')
    self%e28_5 = fu_content_real(nlPtr,'e28_5')
    self%e28_6 = fu_content_real(nlPtr,'e28_6')
    self%e28_7 = fu_content_real(nlPtr,'e28_7')
    self%e28_8 = fu_content_real(nlPtr,'e28_8')
    self%e29 = fu_content_real(nlPtr,'e29')
    self%e30_1 = fu_content_real(nlPtr,'e30_1')
    self%e30_2 = fu_content_real(nlPtr,'e30_2')
    self%e30_3 = fu_content_real(nlPtr,'e30_3')
    
  end subroutine get_parameters_FWI

  
  !================================================================================
  
  subroutine compute_FDI_FWI_daily(self, met_buf, dispMarketPtr, &
                                 & pHorizInterpMet2DispStruct, ifHorizInterp, now)
    !
    ! Computes the fire danger index Fire Weather
    ! Note that it uses the previous-day status, and the changes occur on a daily level
    !
    implicit none

    ! Imported parameters
    type(TFireWeather), intent(inout) :: self
    type(Tfield_buffer), pointer ::  met_buf  ! meteo and internal field buffers
    type(mini_market_of_stacks), pointer :: dispMarketPtr   ! working mini-market
    type(THorizInterpStruct), intent(in) :: pHorizInterpMet2DispStruct
    logical, intent(in) :: ifHorizInterp
    type(silja_time), intent(in) :: now

    ! Local variables
    logical :: ifFound
    integer :: iMeteo, ix, iy, ii
    real :: timestep_sec, r_f, E_d, E_w, arg1, arg2, k_dw, FF_m_0, D_m_0, b, K, r_e, M_r, P_r, &
          & R, U, daylen, tempr, RH, wind, fTmp
    real, dimension(:), pointer :: fLandFractionMet, FF_moist, Duff_moist, Drought, FireWeather, &
                                 & pRain, pRH, pWind, pTempr
    type(silja_field), pointer :: pFld_rainToday, pFld_tempr, pFld_RH, pFld_wind
    !
    ! Basic preparations
    !
    fLandFractionMet => fu_grid_data(fraction_of_land_fld)
    ! fields to compute
    FF_moist => fu_grid_data(self%FFM%pFld_FFMoist)
    Duff_moist => fu_grid_data(self%DM%pFld_DuffMoist)
    Drought => fu_grid_data(self%Drght%pFld_drought)
    FireWeather => fu_grid_data(self%pFld_FireWeather)
    if(error)return
    !
    ! meteo fields
    !
    call find_or_create_field(day_sum_precipitation_flag, dispMarketPtr, species_missing, &
                            & now, surface_level, ifFound, pFld_rainToday, pRain)
    if(.not. ifFound)then
      call set_error('Cannot find in dispersion market:' + fu_quantity_string(day_sum_precipitation_flag), &
                   & 'compute_FDI_FWI_daily')
      return
    endif

    call find_or_create_field(day_mean_relat_humid_2m_flag, dispMarketPtr, species_missing, &
                            & now, level_2m_above_ground, ifFound, pFld_RH, pRH)
    if(.not. ifFound)then
      call set_error('Cannot find in dispersion market:' + fu_quantity_string(day_mean_relat_humid_2m_flag), &
                   & 'compute_FDI_FWI_daily')
      return
    endif
    
    call find_or_create_field(day_mean_temperature_2m_flag, dispMarketPtr, species_missing, &
                            & now, level_2m_above_ground, ifFound, pFld_tempr, pTempr)
    if(.not. ifFound)then
      call set_error('Cannot find in dispersion market:' + fu_quantity_string(day_mean_temperature_2m_flag), &
                   & 'compute_FDI_FWI_daily')
      return
    endif
    
    call find_or_create_field(day_mean_windspeed_10m_flag, dispMarketPtr, species_missing, &
                            & now, level_10m_above_ground, ifFound, pFld_wind, pWind)
    if(.not. ifFound)then
      call set_error('Cannot find in dispersion market:' + fu_quantity_string(day_mean_windspeed_10m_flag), &
                   & 'compute_FDI_FWI_daily')
      return
    endif
    !
    ! Do the actual work
    !
    do iy = 1, ny_dispersion
      do ix = 1, nx_dispersion
        ii = ix + nx_dispersion * (iy -1)
        iMeteo =  fu_grid_index(nx_meteo, ix, iy, pHorizInterpMet2DispStruct)
        if(fLandFractionMet(iMeteo) < 0.01)cycle  ! speedup
        daylen = fu_hours(fu_daylength(fu_lat_geographical_from_grid(real(ix), real(iy), &
                                                                   & dispersion_grid), now))
        ! Adjust units
        tempr = pTempr(ii) - 273.15   ! to Centigrade
        RH = max(0., min(1., pRH(ii))) * 100    ! to %
        wind = pWind(ii) * 3.6        ! to km/hr
        !
        ! Actual calculations include several steps and produce several fields
        ! Equation numbers follow: Van Wagner, C.E., Pickett, T.L., 1985. Equations and 
        ! FORTRAN program for the Canadian Forest Fire Weather Index System (No. 33), 
        ! Forestry technical report. Canadian Forestry Service, Ottawa, Canada.
        !
        ! Fine fuel moisture
        !
        FF_m_0 = self%FFM%e1_1 * (self%FFM%e1_2 - FF_moist(ii)) / (self%FFM%e1_3 + FF_moist(ii))  ! (1)
        ! if rain, update FF_m_0
        r_f = pRain(ii) - 0.5                                                         ! (2)
        if(r_f > 0.01)then
          if(FF_m_0 <= 150)then
            FF_m_0 = FF_m_0 + self%FFM%e3_1 * r_f * exp(-100 / (self%FFM%e3_2 - FF_m_0)) * &      ! (3a)
                            & exp(-self%FFM%e3_3 / max(r_f, 0.1))              
          else
            FF_m_0 = FF_m_0 + self%FFM%e3_1 * r_f * exp(-100 / (self%FFM%e3_2 - FF_m_0)) * &
                            & exp(-self%FFM%e3_3 / max(r_f, 0.1)) +  &                            ! (3b)
                            & self%FFM%e3_4 * (FF_m_0 - self%FFM%e3_5) * (FF_m_0 - self%FFM%e3_5) * &
                            & sqrt(max(r_f,0.))
          endif
        endif  ! r_f > 0: effective rain

        E_d = self%FFM%e4_1 * RH**self%FFM%e4_2 + &                                   ! (4) 
              & (self%FFM%e4_3 * exp(RH * 0.1 - 10) + &
               & self%FFM%e4_4 * (self%FFM%e4_5 - tempr) * &
             & (1.-exp(-self%FFM%e4_6 * RH)))
        E_w = self%FFM%e5_1 * RH**self%FFM%e5_2 + &                                   ! (5)
            & (self%FFM%e5_3 * exp(RH * 0.1 - 10) + &
             & self%FFM%e4_4 * (self%FFM%e4_5 - tempr) * & ! yes, self%e4_* 
             & (1.-exp(-self%FFM%e4_6 * RH)))    ! # yes, self%e4_6
        if(FF_m_0 > E_d)then                                                              !  (6a)
          arg1 = 1.- (RH/100.)**self%FFM%ek_p1
        else
          arg1 = 1.- (1. - RH/100.)**self%FFM%ek_p1
        endif
        if(FF_m_0 < E_w)then
          arg2 = 1. - (RH/100.)**self%FFM%ek_p2
        else
          arg2 = 1.- (1. - RH/100.)**self%FFM%ek_p2                           ! (7a)
        endif
        k_dw = self%FFM%ek_1 * exp(self%FFM%ek_2 * tempr) * &                 ! (6b & 7b)
             & (self%FFM%ek_3 * arg1 + self%FFM%ek_4 * sqrt(wind) * arg2)
        if(FF_m_0 > E_d) FF_m_0 = E_d + (FF_m_0 - E_d) * 10.**(-k_dw)                     ! (8)
        if(FF_m_0 < E_w) FF_m_0 = E_w - (E_w - FF_m_0) * 10.**(-k_dw)                     ! (9)

        FF_moist(ii) = max(0., &
                         & self%FFM%e10_1 * (self%FFM%e10_2 - FF_m_0) / (self%FFM%e10_3 + FF_m_0)) ! (10)

        if(.not. (FF_moist(ii) >=0 .and. FF_m_0 >=0))then
          call msg('Strange fine-fuel moisture params: FF_m_0, FF_moist', FF_m_0, FF_moist(ii))
          call set_error('Problem in a grid cell (ix,iy)='+fu_str(ix)+','+fu_str(iy), 'compute_FDI_FWI')
          return
        endif
        !
        ! Duff moisture
        !
        D_m_0 = self%DM%e12_0 + exp(self%DM%e12_1 - Duff_moist(ii) / self%DM%e12_2)
        ! slope of rain effect
        if(Duff_moist(ii) <= self%DM%e13_1)then
          b = 100. / (self%DM%e13_2 + self%DM%e13_3 * Duff_moist(ii))
        else
          if(Duff_moist(ii) <= self%DM%e13_4)then
            b = self%DM%e13_5 - self%DM%e13_6 * log(Duff_moist(ii))
          else
            b = self%DM%e13_7 * log(Duff_moist(ii)) - self%DM%e13_8
          endif
        endif
        ! drying rate: in the ref paper, the scale is 1e-6 and then 100 at the next line. 
        K = (self%DM%e16_1 * max(0.,(tempr + self%DM%e16_2)) * (100. - RH) * &
                           & max(0., (daylen - self%DM%e16_3))) * 1e-4
        ! moistening the fuel only if a decent rain occurred
        ! effective rain
        r_e = max(0., pRain(ii) * self%DM%e11_1 - self%DM%e11_2)
        ! add moisture, subtract baseline minimum moisture and cut at 0.1% above it to calm the log
        M_r = max(0.1, D_m_0 + 1000 * r_e / (self%DM%e14 + b * r_e) - self%DM%e12_0)
        ! convert actual moisture to moisture code
        P_r = max(0.,self%DM%e15_1 - self%DM%e12_2 * log(M_r))   ! yes, connection to e12

        if(r_e > 0)then
          fTmp = P_r + K
        else
          fTmp = Duff_moist(ii) + K
        endif
        if( .not. fTmp >= 0.)then
          call msg('FWI duff moisture: Strange Duff_moist(ii).', fTmp)
          call set_error('Problem in a grid cell (ix,iy)='+fu_str(ix)+','+fu_str(iy), 'compute_FDI_FWI')
          return
        endif
        Duff_moist(ii) = fTmp
        !
        ! Drought
        !
        if(pRain(ii) > self%Drght%e18_2 / self%Drght%e18_1)then
          if(Drought(ii) / self%Drght%e19_2 < 50)then
            Drought(ii) = self%Drght%e19_2 * &
                        & log(self%Drght%e19_1 / &
                            & (self%Drght%e19_1 * exp(-Drought(ii) / self%Drght%e19_2) + &
                             & self%Drght%e20_1 * (max(self%Drght%e18_1 * pRain(ii) - &
                                                                               & self%Drght%e18_2, &
                                                     & 0.0001))))
          else
            Drought(ii) = self%Drght%e19_2 * &
                        & log(self%Drght%e19_1 / &
                            & (self%Drght%e20_1 * &
                             & (max(self%Drght%e18_1 * pRain(ii) - self%Drght%e18_2, 0.0001))))
          endif
        endif

        ! Note L_f table being just daylen-12hrs cut from below at -1.6
        Drought(ii) = max(0., Drought(ii) + &
                            & 0.5 * (self%Drght%e22_1 * (tempr + self%Drght%e22_2) + &
                                   & max(self%Drght%lf_cutoff, daylen - self%Drght%lf_corr)))
        if(.not. Drought(ii) >= 0.)then
          call msg('FWI drought: Strange Drought(ii)',Drought(ii))
          call set_error('Problem in a grid cell (ix,iy)='+fu_str(ix)+','+fu_str(iy), 'compute_FDI_FWI')
          return
        endif
        !
        ! Fire Weather itself
        !
        ! Initial spread index
        R = 0.0              ! by default, no spread: too wet
        if(self%e25_2 * FF_m_0 < 10)then
          R = self%e26 * exp(self%e24 * wind) * self%e25_1 * exp(-self%e25_2 * FF_m_0) * &
                        & (1.+ self%e25_4 * FF_m_0 ** self%e25_3)
        else
          R = 0.
        endif
        !
        ! Buildup index
        if(Duff_moist(ii) < self%e27_2 * Drought(ii))then
          U = self%e27_1 * Duff_moist(ii) * Drought(ii) / (Duff_moist(ii) + self%e27_2 * Drought(ii))
        else
          U = Duff_moist(ii) - ((1.- self%e27_1 * Drought(ii) / &
                                   & (Duff_moist(ii) + self%e27_2 * Drought(ii))) * &
                              & (self%e27_3 + (self%e27_4 * max(0.,Duff_moist(ii)))**self%e27_5))
        endif
        !
        ! Finally, fire weather index
        !
        if(U <= self%e28_1)then
          FireWeather(ii) = self%e28_2 * (max(U,0.))**self%e28_3 + self%e28_4
        else
          if(self%e28_8 * max(U, self%e28_1) < 50)then
            FireWeather(ii) = self%e29 * R * self%e28_5 / (self%e28_6 + self%e28_7 * &
                                                         & exp(-self%e28_8 * max(U,self%e28_1)))
          else
            FireWeather(ii) = self%e29 * R * self%e28_5 / self%e28_6
          endif
        endif
        if(FireWeather(ii) > 1)then
          FireWeather(ii) = exp(self%e30_1 * &
                              & (self%e30_2 * log(max(FireWeather(ii),1.0)))**self%e30_3)
        endif
        
        if(.not. FireWeather(ii) >= 0.)then
          call msg('FWI FireWeather: Strange FireWeather(ii)',FireWeather(ii))
          call set_error('Problem in a grid cell (ix,iy)='+fu_str(ix)+','+fu_str(iy), 'compute_FDI_FWI')
          return
        endif
      end do  ! ix
    end do  ! iy

  end subroutine compute_FDI_FWI_daily
  
  
        
  !**************************************************************************************
  !**************************************************************************************
  !
  ! Fuel moisture
  !
  !**************************************************************************************
  !**************************************************************************************

  subroutine init_fuel_moisture_idx(self, dispersionMarketPtr, start_time)
    !
    !
    implicit none
    
    ! Imported parameters
    type(TFuelMoisture), intent(inout) :: self
    type(mini_market_of_stacks), pointer :: dispersionMarketPtr
    type(silja_time), intent(in) :: start_time
    
    ! local variables
    logical :: ifFound
    real, dimension(:), pointer :: pValues
    
    ! Read ini namelist and assign the variables
    !
    ! If ini namelist is given, read it and assign the variables
    !
    ! FWI contains three other indices. They all requires several fields to be created in 
    ! a dispersion buffer: forecasted variables
    !
    call find_or_create_field(FDI_fuel_moisture_flag, dispersionMarketPtr, species_missing, &
                            & start_time, surface_level, ifFound, self%pFld_FM, pValues)
  end subroutine init_fuel_moisture_idx

  
  !======================================================================================
  
  subroutine add_input_needs_fuel_moisture(q_met_dynamic, q_met_st, &
                                              & q_disp_dynamic, q_disp_st)
    !
    ! Returns input needs for the transformation routines. 
    !
    implicit none

    ! Imported parameters
    integer, dimension(:), intent(inout) :: q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st

    ! Local variables
    integer :: iTmp

    ! Add needed quantities
    ! Meteo
    iTmp = fu_merge_integer_to_array(temperature_2m_flag,        q_met_dynamic)

  end subroutine add_input_needs_fuel_moisture
  
  
  !======================================================================================
  
  subroutine get_parameters_fuel_moisture(self, nlFDI)
    implicit none
    type(TFuelMoisture), intent(in) :: self
    type(Tsilam_namelist), pointer :: nlFDI
    
    call set_error('Not implemented','get_parameters_fuel_moist')
    
  end subroutine get_parameters_fuel_moisture

  
  !*******************************************************************************************
  
  LOGICAL FUNCTION rules_FDI_defined(rules)
    !
    ! Returns a true value, if the rules has been given a value using
    ! one of the setting functions.
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(Trules_fire_danger), INTENT(in) :: rules

    rules_FDI_defined = fu_true(rules%defined)

  END FUNCTION rules_FDI_defined
  
end module fire_danger_indices
