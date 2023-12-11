MODULE diagnostic_variables
  !
  ! This module contains a set of routines for deriving a diagnostic dispersion
  ! variables. 
  ! Do not mix with derived meteorological quantities: those are 
  ! input fields defined in meteo grid, while these are the outputs of SILAM 
  ! simulations of all kinds. By definition, diagnostic variables are those, which 
  ! are not predicted within the main dispersion cycle of the simualtions but 
  ! - either generated from input data and placed into dispersion environment, 
  !   including dispersion grid and vertical, 
  ! - or derived from SILAM products and do not affect the model simualtions, 
  !   thus being entirely output-oriented. The only reason why these routines are 
  !   not moved outside the SILAM code is that their computation might need internal 
  !   variables, which are either not easily available for output or bulky.
  ! There is a possibility for cross-diagnising one variable from another one. E.g., 
  ! mass can be diagnosed from concentration and the other way round.
  !
  ! Language: FORTRAN-90
  !
  ! Units: NOT NECESSARILY SI. 
  ! Warning is issued in each non-SI case
  !
  ! Authors Julius Vira, Rostislav Kouznetsov
  !         Marje Prank,  Mikhail Sofiev
  !         FMI
  !         <firstname>.<lastname>@fmi.fi
  !
  use cocktail_basic
  use vertical_motion
  use water_in_soil
  use fire_danger_indices
  use fishpack_silam
  !$use omp_lib

  implicit none

  public make_all_diagnostic_fields
  public fu_ifDiagnosticQuantity
  public fu_if_full_meteo_vertical_needed
  public init_wind_diag
  public set_rules_water_in_soil_model
  public set_rules_fire_danger_indices
  public init_water_in_soil_model
  public init_fire_danger_indices
  public ifUseMassfluxes
  public add_diagnostic_input_needs
  public disable_unavail_flds_4adjoint
  public enable_diag_fields

  private make_diagnostic_fields
  private fu_quantity_avail_4adjoint
  private df_DMAT_Vd_correction
  private df_make_cell_size_z_dyn
  private df_cumul_daily_variable

  private df_update_realtime_met_fields
  private df_daily_mean_variable
  private fu_accumulation_type

  private diag_cell_fluxes  ! Get massfluxes for past and future  ignores pressure tendency
  private df_update_cellfluxcorr_rt   ! realtime spt correction for meteo interval

  private exchange_wings
  private TBLKTRI
  private adjust_2D_fluxes
  
  !
  ! Private interpolation structures and variables needed for flux diagnostics
  !
  type(THorizInterpStruct), pointer, save, private :: pHorizInterpU => null(), pHorizInterpV => null(), &
       & pHorizInterpRho => null()
  type(TVertInterpStruct), pointer, save, private ::  pVertInterpRho => null()

  type(silam_vertical), private, save :: integration_vertical
  type(TVertInterpStruct), pointer, save, private :: wind_interp_struct
  integer, dimension(:), allocatable, private, save :: nsmall_in_layer
  real(vm_p), dimension(:,:), allocatable, private, save :: u3d, v3d, eta_dot_3d, w3d

  !
  ! Types of wind diagnostics
  !
  integer, parameter, private :: test_wind = 50000, incompressible = 50001, incompressible_v2 = 50002, &
                        & anelastic_v2 = 50003, from_nwp_omega = 50004, hybrid_top_down = 50005, & 
                        & hardtop = 50006, opentop = 50007, omegatop = 50008, topdown = 50009, hardtop_weighted = 50010
  !
  ! Types of averaging
  !
  integer, parameter, private :: av4mean = 50020, av4max = 50021, av4min = 50022, av4sum = 50023
  
  ! Local parameters: period-integrated quantities we deal with
  integer, dimension(8) :: day_mean_quantities = &
      & (/day_mean_temperature_flag, day_mean_temperature_2m_flag, day_max_temperature_2m_flag, &
        & day_mean_relat_humid_2m_flag, day_min_relat_humid_2m_flag, &
        & day_mean_windspeed_10m_flag, day_max_windspeed_10m_flag, &
        & day_sum_precipitation_flag/)
  integer, dimension(8) :: day_accum_quantities = &
      & (/day_temperature_acc_flag, day_temperature_2m_acc_flag, day_temperature_2m_acc_max_flag, &
        & day_relat_humid_2m_acc_flag, day_relat_humid_2m_acc_min_flag, &
        & day_windspeed_10m_acc_flag, day_windspeed_10m_acc_max_flag, & 
        & day_precipitation_acc_flag/)
  
  type Tdiagnostic_rules
     private
     integer :: continuity_equation = incompressible_v2
     integer :: wind_method = int_missing  ! Do not make any mass fluxes
     integer :: iTestWindType = 1
     integer :: a = int_missing, b = int_missing !Test Wind parameters
     logical :: w_on_half_levels, if_DoStatic
     type(silja_logical) :: ifAllowColdstartDailyVars = silja_undefined    ! used to be .false.
     logical :: if_high_stomatal_conductance=.True.
     type(Trules_water_in_soil) :: rulesWIS
     type(Trules_fire_danger) :: rules_FDI
  end type Tdiagnostic_rules

  type TWindMassBehind           
     real, dimension (:,:,:,:), pointer :: ptr => null()
     !! Intention is  to use (nx or ny, nz, 1:2 mass/flux,1:2 out/in),
     !! but it has to be transferrable as a flat array
  end type TWindMassBehind

  integer, parameter :: iflux = 1
  integer, parameter :: imass = 2
     


  type DispWindStruct
     private
     ! Arrays used for diagnostics with Poisson
     !Grid anisotropy:  Cos(lat) for v_grid with zeros at the edges
     real, dimension (:), allocatable :: CosLat  ! 1:ny_disp+1  
     real, dimension (:), allocatable :: LevelCorrWeight  ! 1:nz_disp
     !Stuff for BLKTRI
     real(fish_kind), dimension (:), allocatable :: AM, BM, CM ! 1:ny_disp 
     real(fish_kind), dimension (:), allocatable :: AN, BN, CN ! 1:nx_disp Stuff for BLKTRI
     ! Excess wind divergence / Velocity potential 
     real(fish_kind), dimension (:,:), allocatable :: Phi ! 1:nx_disp, 1:ny_disp



     ! Arrays used for test wind
     real(r8k), dimension(:,:), allocatable :: StreamFunction  !! 0:ny_disp, 0:nx_disp !Slice of above

     TYPE (fishworkspace) :: fish_work 
  end type DispWindStruct



   !Static stuff and temporary structures
  type(DispWindStruct), private, target :: Disp_wind


  CONTAINS

  !************************************************************************
  
  logical function ifUseMassfluxes(rules)
      type (Tdiagnostic_rules) :: rules
      ifUseMassfluxes = ( rules%wind_method > 0)
  end function


  !************************************************************************

  subroutine init_wind_diag(whole_disp_grid, disp_grid,  disp_vert, meteo_vert, rules)
    !
    ! allocates DispWindStruct and initialises other pieces needed for diagnostics
    !
    implicit none

    ! Imported parameters
    type (Silja_grid), intent(in) :: whole_disp_grid, disp_grid
    type (silam_vertical), intent(in) ::  disp_vert, meteo_vert
    type(Tdiagnostic_rules), intent(in) :: rules

    ! Local variables
    integer ::istat, iTmp, jTmp
    integer :: ny, nx, nz_d
    real a, delta, b, ftmp, std_pr_top
    type(silam_vertical) :: vertOfLayers
    REAL  :: corner_lat_N,corner_lon_E, southpole_lat_N, southpole_lon_E, dx_deg, dy_deg
    LOGICAL :: corner_in_geo_lonlat

    nz_d = fu_NbrOfLevels(disp_vert)
    iTmp = fu_leveltype(disp_vert)
    
    if ( .not. any(iTmp == (/layer_btw_2_height, layer_btw_2_sigma,  layer_btw_2_hybrid/)))then
      call set_error("Strange target  vertical","init_wind_diag")
    endif

    std_pr_top = fu_std_pressure(fu_upper_boundary_of_layer(fu_level(dispersion_vertical,nz_d)))
    call msg("Std pr at domain top (init_wind_diag)", std_pr_top)

     !
     ! Stuff needed for local domain diagnostics
     if (rules%wind_method == test_wind) then
       CALL grid_dimensions(disp_grid, nx, ny) 
       allocate( Disp_wind%StreamFunction((nx+1),(ny+1)), & ! Needed only for test winds....
             & stat=istat)
       if (istat /= 0) then
              call set_error("Failed to allocate DispWindStruct", "init_wind_diag")
       endif
     endif
     
     if (rules%wind_method == hardtop .or. rules%wind_method == hardtop_weighted) then
       call lonlat_grid_parameters(whole_disp_grid,&
                                  & corner_lon_E, corner_lat_N,corner_in_geo_lonlat, &
                                  & nx, ny, &
                                  & southpole_lon_E, southpole_lat_N, &
                                  & dx_deg, dy_deg)
    
       ! Global domain diagnostics
          ! Use poisson
      
          ! Mass-flux correction weight per level 
          allocate(Disp_wind%LevelCorrWeight(nz_d),  stat=istat)
          if (istat /= 0) then
             call set_error("Failed to allocate DispWindStruct", "init_wind_diag")
          endif
          
          if (rules%wind_method == hardtop_weighted) then
             fTmp=0
             do iTmp = 1, nz_d 
                if (fu_leveltype(disp_vert) == layer_btw_2_height) then 
                   call us_standard_atmosphere( &
                            & fu_level_height(fu_level(disp_vert,iTmp, .True.)), &
                            & a, delta, b)
                   Disp_wind%LevelCorrWeight(iTmp) = delta
                   fTmp = fTmp + delta
                else  !! Kind of hybrid level
                   Disp_wind%LevelCorrWeight(iTmp) = fu_hybrid_level_coeff_a(fu_central_level_of_layer(fu_level(disp_vert, iTmp))) + &
                                    & fu_hybrid_level_coeff_b(fu_central_level_of_layer(fu_level(disp_vert, iTmp))) * std_pressure_sl
                   fTmp = fTmp + Disp_wind%LevelCorrWeight(iTmp)
                endif
             enddo
             Disp_wind%LevelCorrWeight(1:nz_d) = Disp_wind%LevelCorrWeight(1:nz_d)/fTmp
           else ! Weight Poisson correction just by mass.
              Disp_wind%LevelCorrWeight(1:nz_d) = 10. ! Weights for levels in wind correction
           endif
           call msg("Colmass weights", Disp_wind%LevelCorrWeight(1:nz_d) )


          if (smpi_adv_rank == 0.) then
            allocate(&
                  & Disp_wind%CosLat(ny+1), &
               & Disp_wind%AM(nx), &
               & Disp_wind%BM(nx), &
               & Disp_wind%CM(nx), &
               & Disp_wind%AN(ny), &
               & Disp_wind%BN(ny), &
               & Disp_wind%CN(ny), &
                  & Disp_wind%Phi(nx,ny), & !Divergence/potential
               & stat=istat)
            if (istat /= 0) then
              call set_error("Failed to allocate DispWindStruct", "init_wind_diag")
           endif

             ! Fill the arrays
                     do iTmp = 2, ny
                !Warning!  Calculate cosine from parameters of non-staggered grid
                Disp_wind%CosLat(iTmp) = cos((corner_lat_N +  (iTmp - 1.5)*dy_deg) * degrees_to_radians )
                !Disp_wind%CosLat(iTmp) = 1.  ! Flat world
             enddo
                     ! It is not a "true" cosine it is rather a metrics
                     Disp_wind%CosLat(1) = 0
                     Disp_wind%CosLat(ny+1) = 0

             ! Metrics on a spherical grid
             Disp_wind%AN(1:ny) = Disp_wind%CosLat(1:ny)
                     Disp_wind%BN(1:ny) =  - Disp_wind%CosLat(1:ny) -  Disp_wind%CosLat(2:ny+1)
                     Disp_wind%CN(1:ny) = Disp_wind%CosLat(2:ny+1)

             Disp_wind%AM(:) = 1.
             Disp_wind%BM(:) = -2.
             Disp_wind%CM(:) = 1.

                     if (.not. fu_ifLonGlobal(whole_disp_grid)) then
                Disp_wind%AM(1) = 0.
                Disp_wind%BM(1) = -1.
                Disp_wind%BM(nx) = -1.
                Disp_wind%CM(nx) = 0.
             endif

             !FIXME TEST only Check if tblktri works at all....
             call msg("Calling tblktri to chek fishpack")
             if (fu_fails( tblktri(1), "Failed fishpack test","init_wind_diag")) return
             call msg("tblktri done")
                    
                    
             ! Prepare internal fishpack structures for our geometry....
             iTmp = 1  ! x boundary -- zero 
             if (fu_ifLonGlobal(whole_disp_grid)) iTmp = 0 ! periodic x boundaries
             Disp_wind%Phi(:,:) =  0 
             call  BLKTRI (0,  1,  ny, Disp_wind%AN, Disp_wind%BN, Disp_wind%CN, &
                           &  iTmp, nx, Disp_wind%AM, Disp_wind%BM, Disp_wind%CM, &
                             &   nx, Disp_wind%Phi, istat, Disp_wind%fish_work)
             if (istat /= 0) then
               call msg("BLKTRI returned, ny", istat, ny)
               call msg('An:',Disp_wind%AN)
               call msg('Bn:',Disp_wind%BN)
               call msg('Cn:',Disp_wind%CN)
               call msg('Am:',Disp_wind%AM)
               call msg('Bm:',Disp_wind%BM)
               call msg('Cm:',Disp_wind%CM)
               call set_error ("Failed poisson solver init", "init_wind_diag")
               return
             endif
          endif !Poisson-stuff only at rank==0
          if (std_pr_top > 2000) then
            call set_error("Attempt to usehardtop diagnostics for shallow (Top>2 hPa) domains","init_wind_diag")
            return 
          endif
        else
           if (std_pr_top < 2000 .and. rules%wind_method /= topdown) then
             call msg("Std pr at domain top", std_pr_top)
             call set_error("Attempt to use non-hardtop diagnostics for thick (Top<2 hPa) domains","init_wind_diag")
             return 
           endif
       endif ! Use poisson

  end subroutine init_wind_diag

  
  !************************************************************************************
  
  subroutine set_rules_water_in_soil_model(nlSetup, rules)
    !
    ! Checks the necessity and initialises prognostic model for water in soil 
    !
    implicit none  
    
    ! Imported parameters
    type(Tsilam_namelist), pointer :: nlSetup
    type(Tdiagnostic_rules), intent(inout) :: rules

    ! Set the rules
    !
    call set_rules_WIS(nlSetup, rules%rulesWIS)

  end subroutine set_rules_water_in_soil_model

  
  !************************************************************************************
  
  subroutine init_water_in_soil_model(dispMarketPtr, meteoMarketPtr, rules, start_time)
    !
    ! Checks the necessity and initialises prognostic model for water in soil 
    !
    implicit none  
    
    ! Imported parameters
    type(mini_market_of_stacks), pointer :: meteoMarketPtr, dispMarketPtr
    type(Tdiagnostic_rules), intent(in) :: rules
    type(silja_time), intent(in) :: start_time

    ! Calls the corresponding initialization, if needed
    !
    if(defined(rules%rulesWIS)) &
              & call init_WIS(dispMarketPtr, meteoMarketPtr, rules%rulesWIS, start_time)

  end subroutine init_water_in_soil_model


  !************************************************************************************

  subroutine set_rules_fire_danger_indices(nlIS4FIRES, rules, out_quantities_st)
    !
    ! sets rules for fire danger indices
    !
    implicit none
    type(Tsilam_namelist), pointer :: nlIS4FIRES
    type(Tdiagnostic_rules), intent(inout) :: rules
    integer, dimension(:), intent(in) :: out_quantities_st
    
    call set_rules_FDI(nlIS4FIRES, rules%rules_FDI, out_quantities_st)

  end subroutine set_rules_fire_danger_indices
  
  
  !************************************************************************************
  
  subroutine init_fire_danger_indices(rules, dispersionMarketPtr, start_time)
    !
    ! Checks the necessity and initialises prognostic model for fire danger indices
    !
    implicit none  
    
    ! Imported parameters
    type(mini_market_of_stacks), pointer :: dispersionMarketPtr
    type(Tdiagnostic_rules), intent(inout) :: rules
    type(silja_time), intent(in) :: start_time

    ! Send the request to the corresponding module, just open the rules
    !
    call init_fire_danger_indices_FDI(rules%rules_FDI, dispersionMarketPtr, start_time)

  end subroutine init_fire_danger_indices


  !************************************************************************************

  subroutine set_diagnostic_rules(nlMeteoSetup, nlStandardSetup, rules)
    !
    ! Set the vertical wind method and the half_levels option.
    ! 
    implicit none
    type(Tsilam_namelist), pointer :: nlMeteoSetup, nlStandardSetup
    type(Tdiagnostic_rules), intent(out) :: rules
    integer :: iTmp

    character(len=fnlen) :: strTmp
    character (len=*), parameter :: sub_name = "set_diagnostic_rules"

    rules%w_on_half_levels = .true.
    call msg('Vertical velocity computed at half-levels')
    call msg('')

    strTmp = fu_str_l_case(fu_content(nlStandardSetup, 'continuity_equation'))
    iTmp = INDEX(trim(strTmp),' ')
    if (iTmp > 0) then 
        strTmp = strTmp(1:iTmp-1)
    endif
    select case (strTmp)
      case ('incompressible')
        rules%continuity_equation = incompressible
        call msg('Incompressible continuity eqn. formulation')
      case ( 'incompressible_v2')
        rules%continuity_equation = incompressible_v2
        call msg('Revised incompressible continuity eqn. formulation')
      case ( 'anelastic_v2')
        rules%continuity_equation = anelastic_v2
        call msg('Revised anelastic continuity eqn. formulation')
      case ( 'nwp_omega')
        rules%continuity_equation = from_nwp_omega
        call msg('Vertical velocity derived from NWP omega field')
      case ( 'test_wind')
        rules%continuity_equation = test_wind
        if (iTmp > 0) then ! Requires a type of wind, may be, more parameters
          strTmp = fu_str_l_case(fu_content(nlStandardSetup, 'continuity_equation'))
          strTmp = trim(strTmp(iTmp+1:))
          read(unit=strTmp, fmt=*, iostat=iTmp)  rules%iTestWindType
          iTmp = INDEX(trim(strTmp),' ')
          if  (iTmp > 0)then
            read(unit=strTmp, fmt=*, iostat=iTmp) rules%a
            strTmp = trim(strTmp(iTmp+1:))
          endif
          iTmp = INDEX(trim(strTmp),' ')
          if  (iTmp > 0)then
            read(unit=strTmp(iTmp+1:), fmt=*, iostat=iTmp) rules%b
          endif
        else
          call set_error('test_wind requires at least 1 parameter: wind type index',sub_name)
          return
        endif
        call msg("Test wind parameters", rules%a, rules%b)
        call msg_warning('Fake wind wind')
      case ( 'hybrid_top_down')
        rules%continuity_equation = hybrid_top_down
        call msg('Top-down eta dot diagnostic')
      case default
        call msg('continuity_equation = ' // trim(strTmp))
        call msg_warning('Strange or no continuity_equation defined in standard setup,' + &
                       & 'will use differential form v2',sub_name)
        rules%continuity_equation = incompressible_v2
    end select
    

    ! To replace continuity_equation at some point
    strTmp = fu_str_l_case(fu_content(nlStandardSetup, 'wind_diagnostics'))
    iTmp = INDEX(trim(strTmp),' ')
    if (iTmp > 0) then 
        strTmp = strTmp(1:iTmp-1)
    endif
    select case (strTmp)
      case ('opentop')
        rules%wind_method = opentop
        call msg('Opentop vertical mass_flux')
      case ('topdown')
        rules%wind_method = topdown
        call msg('topdown vertical mass_flux')
      case ( 'hardtop')
        rules%wind_method = hardtop
        call msg('Hartdtop vertical mass_flux')
      case ( 'hardtop_weighted')
        rules%wind_method = hardtop_weighted
        call msg('Hartdtop  mass_fluxes with weighting')
      case ( 'none')
        rules%wind_method = int_missing
        call msg('No mass_flux diagnostics')
      case ( 'omegatop')
        rules%wind_method = omegatop
        call set_error("Wind diagnostics from omega is not implemented yet",&
                       & sub_name)
        return
      case ( 'test_wind')
        rules%wind_method = test_wind
        if (iTmp > 0) then ! Requires a type of wind, may be, more parameters
          strTmp = fu_str_l_case(fu_content(nlStandardSetup, 'wind_diagnostics'))
          strTmp = trim(strTmp(iTmp+1:))
          read(unit=strTmp, fmt=*, iostat=iTmp)  rules%iTestWindType
          iTmp = INDEX(trim(strTmp),' ')
          if  (iTmp > 0)then
            read(unit=strTmp, fmt=*, iostat=iTmp) rules%a
            strTmp = trim(strTmp(iTmp+1:))
          endif
          iTmp = INDEX(trim(strTmp),' ')
          if  (iTmp > 0)then
            read(unit=strTmp(iTmp+1:), fmt=*, iostat=iTmp) rules%b
          endif
        else
          call set_error('test_wind requires at least 1 parameter: wind type index',&
                           &sub_name)
          return
        endif
        call msg("Test wind parameters", rules%a, rules%b)
        call msg_warning('Fake massflux wind:'+fu_str(rules%iTestWindType))
      case default
        call msg('wind_diagnostics = ' // trim(strTmp))
        call msg("can be wind_diagnostics = opentop / hardtop / none")
        call msg_warning('Strange or no wind_diagnostics defined in standard setup,' + &
                       & 'will use winds opentop',sub_name)
        rules%wind_method = opentop
    end select

    strTmp = fu_str_l_case(fu_content(nlStandardSetup, 'stomatal_conductance'))
    select case (strTmp)
      case ('high')
        rules%if_high_stomatal_conductance = .true.
      case ('low')
        rules%if_high_stomatal_conductance = .false.
      case default
         call msg_warning('Strange or no stomatal conductance strength in standard setup,' + &
                       & 'will use high, suitable for CB5',sub_name)
        rules%if_high_stomatal_conductance = .true.
     end select

    rules%if_DoStatic = .true.  ! for the first set of diagnostic calls

    strTmp = fu_str_l_case(fu_content(nlStandardSetup, 'allow_coldstart_daily_variables'))
    rules%ifAllowColdstartDailyVars = silja_undefined
    if (strTmp=='yes')then
      rules%ifAllowColdstartDailyVars = silja_true
    elseif( strTmp=='no')then
      rules%ifAllowColdstartDailyVars = silja_false
    endif

    !! Backward compatibility
    strTmp = fu_str_l_case(fu_content(nlStandardSetup, 'allow_coldstart_day_temperature'))
    if (strTmp=='yes')then
      rules%ifAllowColdstartDailyVars = silja_true
    endif

    ! NO. 
    ! Seemingly a right place to set the soil water diagnostic rules
    ! Since here we have no clue if this is needed at all, just set them - this is always possible
    ! Later it will be decided if we want them
    !
!    call set_rules_WIS(nlMeteoSetup, rules%rulesWIS)
    
  end subroutine set_diagnostic_rules

  
  !****************************************************************************************
  
  subroutine add_diagnostic_input_needs(rules, ifVertical_MeteoDependent, &
                                      & q_met_dyn, q_met_st, q_disp_dyn, q_disp_st)
    !
    ! Includes the meteorological quantities needed for the diagnostics.
    ! The arrays must have all requestes from specific modules, so we can check what is
    ! needed to diagnose
    !
    implicit none
    
    ! Imported parameters
    type(Tdiagnostic_rules), intent(inout) :: rules
    logical, intent(in) :: ifVertical_MeteoDependent
    INTEGER, DIMENSION(:), INTENT(inout) :: q_met_dyn, q_met_st, q_disp_dyn, q_disp_st
    
    ! Local variables
    integer :: iTmp, iQ, iQ_rate

    ! The basic parameters that will be needed always
    !
    iTmp = fu_merge_integer_to_array(height_flag, q_met_dyn)
    iTmp = fu_merge_integer_to_array(surface_pressure_flag, q_met_dyn)
    !
    ! High-level diagnostic modules. May require some second-level diagnostics, so need to be first
    !
    ! Water in soil. Call it only if the corresponding quantities are needed
    !
    if (defined(rules%rulesWIS)) then
      if(fu_quantity_in_quantities(soil_moisture_vol_frac_nwp_flag, q_met_dyn) .or. &
       & fu_quantity_in_quantities(water_in_soil_srf_grav_flag, q_met_dyn) .or. &
       & fu_quantity_in_quantities(water_in_soil_deep_grav_flag, q_met_dyn))then
        call add_input_needs_WIS(rules%rulesWIS, q_met_dyn, q_met_st, q_disp_dyn, q_disp_st)
      endif
    endif
    !
    ! Fire danger indices. Should any of them be called, have to add needs
    !
    call add_input_needs_FDI(rules%rules_FDI, q_met_dyn, q_met_st, q_disp_dyn, q_disp_st)
    if(error)return
    !
    ! Basic meteorology in application to dispersion
    !
    !
    ! Vd correction
    !
    if(fu_quantity_in_quantities(Vd_correction_DMAT_flag, q_disp_dyn))then

      iTmp = fu_merge_integer_to_array(total_precipitation_int_flag, q_met_st)
      iTmp = fu_merge_integer_to_array(temperature_2m_flag, q_met_dyn)
      iTmp = fu_merge_integer_to_array(friction_velocity_flag, q_met_dyn)
    endif

    !
    ! z-size of grid cells
    !
    if(fu_quantity_in_quantities(cell_size_z_flag, q_disp_dyn)) &

      & iTmp = fu_merge_integer_to_array(surface_pressure_flag, q_met_dyn)
    
    !
    ! Fluxes at the edge of the grid cells
    !
    if(fu_quantity_in_quantities(disp_flux_celleast_flag, q_disp_dyn) .or. &
         & fu_quantity_in_quantities(disp_flux_cellnorth_flag, q_disp_dyn) .or. &
         & fu_quantity_in_quantities(disp_flux_celltop_flag, q_disp_dyn) .or. &
         & fu_quantity_in_quantities(disp_cell_airmass_flag, q_disp_dyn))then

      iTmp = fu_merge_integer_to_array(u_flag, q_met_dyn)
      iTmp = fu_merge_integer_to_array(v_flag, q_met_dyn)
      iTmp = fu_merge_integer_to_array(surface_pressure_flag, q_met_dyn)
      iTmp = fu_merge_integer_to_array(height_flag, q_met_dyn)
    endif

    !
    ! Mean quantities require corresponding cumulative ones (those are in dynamic minimarket)
    !
    do iQ = 1, size(day_mean_quantities)
      if(fu_quantity_in_quantities(day_mean_quantities(iQ), q_disp_st)) &
          & iTmp = fu_merge_integer_to_array(fu_cumul_quantity_for_mean_one(day_mean_quantities(iQ)), &
                                           & q_disp_dyn)
    end do
    !
    ! Cumulative quantities require corresponding rates (may or may not be in dynamic minimarket)
    !
    do iQ = 1, size(day_accum_quantities)
      if(fu_quantity_in_quantities(day_accum_quantities(iQ), q_disp_dyn))then
        iQ_rate = fu_rate_quantity_for_acc(day_accum_quantities(iQ))    ! rate quantity
        if(fu_SILAM_disp_grid_quantity(iQ_rate))then                ! dispersion grid or meteo?
          if(fu_realtime_quantity(iQ_rate))then                     ! multi-time or single-time?
            iTmp = fu_merge_integer_to_array(iQ_rate, q_disp_st)
          else
            iTmp = fu_merge_integer_to_array(iQ_rate, q_disp_dyn)
          endif
        else
          if(fu_realtime_quantity(iQ_rate))then                  ! multi-time or single-time?
            iTmp = fu_merge_integer_to_array(iQ_rate, q_met_st)
          else
            iTmp = fu_merge_integer_to_array(iQ_rate, q_met_dyn)
          endif
        endif
      endif
    end do

  end subroutine add_diagnostic_input_needs

  !****************************************************************************************

  subroutine disable_unavail_flds_4adjoint(buf_ptr, timeStart, timeEnd)
    !
    ! In case of 4D-VAR, adjoint run may not need the full list of meteo quantities
    ! For instance, cumulative temeprature, which drives pollen sources, is not needed
    ! for adjoint.
    ! Removing these fields from existing meteo buffer and minimarkets is a huge task,
    ! entirely unnecessary because these fields will be needed back in the forward run
    ! Therefore, we just declare them always valid, thus turning their update off.
    !
    implicit none
    
    ! Imported parameters
    type(Tfield_buffer), pointer :: buf_ptr
    type(silja_time), intent(in) :: timeStart, timeEnd  ! interval of validity
    
    ! Local variables
    integer :: iQbuf, Qidx
    type(silja_field_id), pointer :: idPtr
    
    ! Stupidity check
    !
    if(fu_fails(timeStart < timeEnd,'timeStart is not earlierthan timeEnd','disable_unavail_flds_4adjoint'))return
    !
    ! Scan the buffer and reset validity time for each field
    !
    do Qidx = 1, size(buf_ptr%buffer_quantities)
      ! Finished?
      if(buf_ptr%buffer_quantities(Qidx) == int_missing)return ! all done
      ! available?
      if(fu_quantity_avail_4adjoint(buf_ptr%buffer_quantities(Qidx)))cycle
      ! disable
      if (fu_dimension(buf_ptr, Qidx) == 4) then 
        idPtr => buf_ptr%p4d(Qidx)%present%p2d(1)%idPtr
      else
        idPtr => buf_ptr%p2d(Qidx)%present%idPtr
      endif
      call set_valid_time(idPtr, timeStart)
      if(error)return
      call set_validity_length(idPtr, timeEnd - timeStart)
      if(error)return
    end do ! Qidx
  
  end subroutine disable_unavail_flds_4adjoint

  
  !****************************************************************************************

  subroutine enable_diag_fields(buf_ptr, now)
    !
    ! In case of 4D-VAR, adjoint run may not need the full list of meteo quantities
    ! For instance, cumulative temeprature, which drives pollen sources, is not needed
    ! for adjoint.
    ! Removing these fields from existing meteo buffer and minimarkets is a huge task,
    ! entirely unnecessary because these fields will be needed back in the forward run
    ! Therefore, we just declare them always valid, thus turning their update off.
    !
    implicit none
    
    ! Imported parameters
    type(Tfield_buffer), pointer :: buf_ptr
    type(silja_time), intent(in) :: now
    
    ! Local variables
    integer :: Qidx
    type(silja_field_id), pointer :: idPtr
    
    !
    ! Scan the buffer and reset validity time to zero for each field
    !
    do Qidx = 1, size(buf_ptr%buffer_quantities)
      ! finished?
      if(buf_ptr%buffer_quantities(Qidx) == int_missing)return ! all done
      ! Need to reset?
      if(fu_quantity_avail_4adjoint(buf_ptr%buffer_quantities(Qidx)))cycle
      ! Enable
      if (fu_dimension(buf_ptr, Qidx) == 4 ) then 
        idPtr => buf_ptr%p4d(Qidx)%present%p2d(1)%idPtr
      else
        idPtr => buf_ptr%p2d(Qidx)%present%idPtr
      endif
      call set_valid_time(idPtr, now)
      if(error)return
      call set_validity_length(idPtr, zero_interval)
      if(error)return
    end do ! Qidx
  
  end subroutine enable_diag_fields
  
  
  !************************************************************************************
  
  logical function fu_quantity_avail_4adjoint(quantity)
    !
    ! Some quantities cannot be computed in adjoint mode if the are in single-time stack
    ! They cannot be just skipped - it can be an error in setup - but explicitly
    ! disable them upon request is OK
    !
    implicit none
    
    ! Imported parameter
    integer, intent(in) :: quantity
    
    fu_quantity_avail_4adjoint = .not. (any(day_mean_quantities(:) == quantity) &
                                      & .or. any(day_accum_quantities(:) == quantity))
    
  end function fu_quantity_avail_4adjoint


  !****************************************************************************************

  subroutine make_all_diagnostic_fields(meteoMarketPtr, meteo_buf_ptr, &         ! meteo market & buffer
                                      & dispersionMarketPtr, disp_buf_ptr, & ! disp market & buffer
                                      & outputMarketPtr, output_buf_ptr, & ! output market & buffer
                                      & meteo_dyn_shopping_list, & ! full list of meteo dynamic fields
                                      & disp_dyn_shopping_list, &  ! full list of disp dynamic fields
                                      & output_dyn_shopping_list, &  ! full list of output dynamic fields
                                      & disp_stat_shopping_list, & ! full list of disp static fields
                                      & output_stat_shopping_list, & ! full list of output static fields
                                      & wdr, &
                                      & diagnostic_rules, &
                                      & ifNewMeteoData, &
                                      & now, timestep)
    !
    ! This is the main umbrella for deriving the diagnostic variables of all types.
    ! For meteo-only diagnosing it uses the existing derived_field_quantities modules
    ! but it is mainly historical. 
    ! Dispersion/output markets is the next point: diagnosing there may need the meteo information
    ! and the fields from dispersion/output market. Note that meteo static full list is handled by 
    ! the physiography module while dispersion static list is treated here
    !
    ! Boundary stuff should also be here but so far it is not.
    !
    implicit none

    ! Imported parameters
    type(mini_market_of_stacks), pointer :: meteoMarketPtr, dispersionMarketPtr, outputMarketPtr
    type(Tfield_buffer), pointer :: meteo_buf_ptr, disp_buf_ptr, output_buf_ptr
    type(silja_shopping_list), intent(inout) :: meteo_dyn_shopping_list
    type(silja_shopping_list), intent(in) :: disp_dyn_shopping_list, disp_stat_shopping_list, &
                                           & output_dyn_shopping_list, output_stat_shopping_list
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep
    type(silja_wdr), intent(in) :: wdr
    type(Tdiagnostic_rules), intent(inout) :: diagnostic_rules
    logical, intent(in) :: ifNewMeteoData
    
    ! Local variables
    type(silja_time), save :: valid_border_time_1 = really_far_in_past, &
                            & valid_border_time_2 = really_far_in_past
    type(silja_interval) :: meteo_time_shift
    integer :: iTmp, nDiagnosed
    type(silja_field), pointer :: fieldPtr
    logical :: ifOK, if_update_realtime_met_fields
    !
    ! If meteo is shifted, take this into account
    !
    meteo_time_shift = fu_meteo_time_shift(wdr)

    !
    ! Should be called before make_derived_meteo_fields to get tz_offset field properly
    !
    call update_timezone_offsets(now)
    call SolarSetup(now)  ! Solar zenith angle, sunset, sunrize etc would be right
    !
    ! Part 1: deal with the meteo data. 
    !
    CALL fix_shopping_time_boundaries(meteo_dyn_shopping_list, now + meteo_time_shift, &
                                                             & now + timestep + meteo_time_shift)
    if(error)return
    !
    ! Derived meteo quantities
    !
    CALL make_derived_meteo_fields(meteoMarketPtr, meteo_dyn_shopping_list, wdr, &
                                 & diagnostic_rules%if_high_stomatal_conductance, ifNewMeteoData)
    if(error)return

    !
    ! Set meteo buffer pointers
    !
    call msg("make_all_diagnostic_fields sets met_buf pointers")

    CALL arrange_buffer(meteoMarketPtr, fu_met_src(wdr,1), now + meteo_time_shift, timestep,&
                      & meteo_buf_ptr, meteo_verticalPtr)
    if(error)return
    !
    ! Preparatory: refine all vertical interpolation coefficients in the whole pool
    ! of the vertical interpolation structures.
    !
    call start_count('Refine vertical coefs v2')
    call refine_all_vert_interp_coefs(meteo_buf_ptr, ifNewMeteoData, &
                                    & now + meteo_time_shift + timestep * 0.5)

!    if(ifMeteo2DispVertInterp) call refine_interp_vert_coefs_v2(pMeteo2DispVertInterp, & 
!                                                              & meteo_buf_ptr, &
!                                                              & now, &
!                                                              & pMeteo2DispHorizInterp)
    call stop_count('Refine vertical coefs v2')
    if(error)return

    !
    ! Need to update real-time meteo fields prior to start making dispersion ones.
    ! valid_border_time_1, valid_border_time_2 are remembered through the model run
    ! do not reset them
    !
    if(fu_between_times(now + meteo_time_shift + timestep * 0.5, &
                      & valid_border_time_1, valid_border_time_2, .true.))then
      !
      ! There is nothing to do, all is still valid
      !
      if_update_realtime_met_fields = .false.
    else
      !
      ! work it out
!      call msg("Updating realime for time:"+fu_str(now + timestep * 0.5))
!      call msg("Valid range = "+fu_str(valid_border_time_1)+ ":" +fu_str(valid_border_time_2))

      ! Reset_times is not used since at least v.5.7, so commented it here
      !
      if_update_realtime_met_fields = .true.

      call df_update_realtime_met_fields(meteoMarketPtr, meteo_buf_ptr, &
                                       & now + meteo_time_shift + timestep * 0.5, wdr, diagnostic_rules,&
                                       & valid_border_time_1, valid_border_time_2) ! , .true.)   ! interval of validity, output
      if(error)return
    endif
!    !                        
!    ! Probably, here is the right point to take care of water in soil
! Wrong solution: the update must use dispersion market rather than buffer. For now, commented it all out
!    !
!    if(defined(diagnostic_rules%rulesWIS))then
!      call update_water_in_soil(meteo_buf_ptr, disp_buf_ptr, now + meteo_time_shift, timestep, diagnostic_rules%rulesWIS)
!      if(error)return
!    endif
    !
    ! Finally, derive the dispersion diagnostic fields using their shopping lists
    !
    call make_diagnostic_fields(meteo_buf_ptr, & ! meteo buffer
                              & dispersionMarketPtr, & ! disp market
                              & disp_dyn_shopping_list, & ! full list of disp dynamic fields
                              & disp_stat_shopping_list, & ! full list of disp static fields
                              & diagnostic_rules, wdr, nDiagnosed)
    if(error)return
    !
    ! In case we completely failed the request, something must be wrong
    !
    iTmp = fu_nbr_of_quantities(disp_dyn_shopping_list,.false.) &
       & + fu_nbr_of_quantities(disp_stat_shopping_list,.false.)
    if(nDiagnosed == 0 .and. iTmp > 0)then
      call set_error('Diagnostic quantities requested but none diagnosed out of:' + fu_str(iTmp), &
                   & 'make_all_diagnostic_fields for dispersionMarket')
      return
    endif

    call msg_test('Setting dispersion pointers...')
    call arrange_buffer(dispersionMarketPtr, met_src_missing, now + meteo_time_shift, timestep, & 
                      & disp_buf_ptr,  dispersion_verticalPtr)
    if(error)return
    call msg_test('Done dispersion pointers...')
    !
    ! fields must be valid for the middle of step
    ! The sub prepares the fields with certain validity period. As soon as out target time is outside
    ! we call it again. 
    ! Note that its call must be synchronised with the above update_realtime_met_fields for meteo market
    ! This sub also handles daily-mean quantities
    !
    if(if_update_realtime_met_fields)then
#ifdef DEBUG
      call msg("Updating realime for time:" + fu_str(now + timestep * 0.5))
      call msg("Valid range =" + fu_str(valid_border_time_1) + ":" + fu_str(valid_border_time_2))
#endif
      ! Reset_times is not used since at least v.5.7, so commented it here
      !
!      call df_update_realtime_met_fields(meteoMarketPtr, meteo_buf_ptr, &
!                                       & now + meteo_time_shift + timestep * 0.5, wdr, diagnostic_rules,&
!                                       & valid_border_time_1, valid_border_time_2) ! , .true.)   ! interval of validity, output
!      if(error)return
!
#ifdef DEBUG_DIAG
      call msg('')
      call msg('DISPERSION market')
      call report(dispersionMarketPtr, .false.)
      call msg('')
#endif
      call df_update_realtime_met_fields(dispersionMarketPtr, disp_buf_ptr, &
                                       & now + meteo_time_shift + timestep * 0.5, wdr, diagnostic_rules, &
                                       & valid_border_time_1, valid_border_time_2)  !, .false.)   ! interval of validity, output. Do not reset it
      if(error)return
#ifdef DEBUG_DIAG
      call msg('')
      call msg('DISPERSION market after update')
      call report(dispersionMarketPtr, .false.)
      call msg('')
#endif
    endif  ! outside time interval

    !
    ! The meteorology is over, including dispersion-specific quantities. 
    ! Something very specific: fire danger indices and water in soil
    !
    ! Fire danger indices. Determined not by their quantity but whether rules are defined or not
    !
    if(defined(diagnostic_rules%rules_FDI))then
      call compute_fire_danger_indices(diagnostic_rules%rules_FDI, meteo_buf_ptr, now + meteo_time_shift, &
!                                     & ifHorizInterp, pHorizInterpStruct, &
                                     & dispersionMarketPtr)
      if(error)return
!      call arrange_supermarket(dispersionMarketPtr, .true.)
!      if(error)return
      nDiagnosed = nDiagnosed + 1
    endif

    !
    ! Water in soil
    !
    if(defined(diagnostic_rules%rulesWIS))then
      !
      ! Water in soil appears here but it has to deal with dispersion market, not buffer since
      ! buffer pointers are not yet set
      !        call update_water_in_soil(met_buf, dispMarketPtr, obstimes, rules%rulesWIS)
      call set_error('Water in soil is not implemented to the end','make_all_diagnostic_fields')
      return
!      if(error)return
!      call arrange_supermarket(dispMarketPtr, .true.)
!      if(error)return
!      nDiagnosed = nDiagnosed + 1
    endif
    

!      call msg("**********************Diagnozing output market with shopping list:")
!      call report(output_dyn_shopping_list)
!      call msg("********************")
    !
    ! Finally, the output market: make the fields if needed and then updated them if needed
    !
    call make_diagnostic_fields(meteo_buf_ptr,   &   ! meteo buffer
                              & outputMarketPtr, &   ! output market
                              & output_dyn_shopping_list,  &  ! full list of output dynamic fields
                              & output_stat_shopping_list, & ! full list of output static fields
                            !  & now, &
                              & diagnostic_rules, wdr, nDiagnosed)
    if(error)return
    !
    ! In case we completely failed the request, something must be wrong
    !
    iTmp = fu_nbr_of_quantities(output_dyn_shopping_list,.false.) &
       & + fu_nbr_of_quantities(output_stat_shopping_list,.false.)
    if(nDiagnosed == 0 .and. iTmp > 0)then
      call set_error('Diagnostic quantities requested but none diagnosed out of:' + fu_str(iTmp), &
                   & 'make_all_diagnostic_fields for outputMarket')
      return
    endif
    !
    ! After all derivations are done, arrange the minimarket
    !
    call arrange_supermarket(dispersionMarketPtr, .true., .true.)
    if(error)return


!    call msg_test('Setting output pointers...')
!    call arrange_buffer(outputMarketPtr, met_src_missing, now, & 
!                      & weight_past, output_buf_ptr, output_gridPtr, output_verticalPtr)
!    if(error)return
    
    diagnostic_rules%if_DoStatic = .false.  ! no need to repeat it at each time step

  end subroutine make_all_diagnostic_fields


  !***************************************************************************************************

  subroutine make_diagnostic_fields(met_buf, dispMarketPtr, &
                                  & dqListDyn, dqListStat,&
!                                  & ifMeteo2DispHorizInterp, ifMeteo2DispVertInterp, &
!                                  & pMeteo2DispHorizInterp, pMeteo2DispVertInterp, &
!                                  & now, &
                                  & diagnostic_rules, wdr, nDiagnosed)
    !
    ! Creates the quantities in list in every stack in the minimarket, where possible
    ! Most of arrange_minimarket calls have been commented out. Only possible chains
    ! are processed at once: if winds or 3d fields are needed for later diagnostics
    ! A final arranging is done at the end.
    !
    implicit none

    ! Imported variables
    type(silja_shopping_list), intent(in) ::  dqListDyn, dqListStat
    type(mini_market_of_stacks), pointer :: dispMarketPtr
    type(Tfield_buffer), pointer :: met_buf
!    type(THorizInterpStruct), pointer :: pMeteo2DispHorizInterp
!    type(TVertInterpStruct), pointer :: pMeteo2DispVertInterp
!    logical :: ifMeteo2DispHorizInterp, ifMeteo2DispVertInterp
!    type(silja_time), intent(in) :: now
    type(Tdiagnostic_rules), intent(inout) :: diagnostic_rules
    type(silja_wdr), intent(in) :: wdr
    integer, intent(out) :: nDiagnosed

    ! Local variables
    type(meteo_data_source), dimension(max_met_srcs) :: metSrcsMultitime, metSrcsSingletime
    integer :: nMetSrcsMultitime, nMetSrcsSingletime, ntimes, iTmp
    TYPE(silja_time), DIMENSION(2) :: obstimes
    type(silja_stack), pointer :: stackPtr 
    type(silja_field), pointer :: fieldPtr
    type(silja_field_id), pointer :: idIn
    type(silja_field_id) :: idRequest
    integer :: iVarLst, iMetSrc, iT, iFldStack, shopQ, iVar, nVars, iQ
    integer, dimension(max_quantities) :: iArr, lst, lst_st
    logical :: iffound,  ifvalid, ifHorizInterp, ifVertInterp      ! if interpolation needed
    type(THorizInterpStruct), pointer :: pHorizInterpStruct ! meteo 2 dispersion horizontal
    type(TVertInterpStruct), pointer :: pVertInterpStruct   ! meteo 2 dispersion vertical
    !
    ! Static market is set by physiography
    ! no realtime fields here
    !
    IF (.NOT.defined(dqListStat)) cALL msg_warning('undefined static list','make_diagnostic_fields')
    IF (.NOT.defined(dqListDyn)) THEN
      CALL msg_warning('undefined dynamic list','make_diagnostic_fields')
    END IF
      
!    call msg("make_diagnostic_fields got Dyn list")
!    call report(dqListDyn)
!    call msg("make_diagnostic_fields got ST list")
!    call report(dqListDyn)
!    call msg("oops.")

    ! Reorder quantities in the dynamic market
    ! cell_size_z_flag must be the first
    !
    nVars = fu_nbr_of_quantities(dqListDyn,.false.)
    nDiagnosed = 0
    !
    ! Get interpolation structures
    !
    ifHorizInterp = .not. meteo_grid == dispersion_grid
    ifVertInterp = .not. fu_cmp_verts_eq(meteo_vertical, dispersion_vertical)
    if(ifHorizInterp) &
          & pHorizInterpStruct => fu_horiz_interp_struct(meteo_grid, dispersion_grid, linear, .true.)
    if(error)return
    if(ifVertInterp) &
          & pVertInterpStruct => fu_vertical_interp_struct(meteo_vertical, dispersion_vertical, &
                                                         & dispersion_grid, linear, &
                                                         & one_hour, 'main_meteo_to_disp')
    !
    ! Do the main job: call diagnostic routines one by one
    !
    ! Cell size-z
    !
    if(fu_quantity_in_list(cell_size_z_flag, dqListDyn) .or. &
     & fu_quantity_in_list(cell_size_z_flag, dqListStat))then
      call msg("Diagnosing dynamic dispersion vertical cell size..")
      call check_obstimes(cell_size_z_flag)
      call df_make_cell_size_z_dyn(met_src_missing, met_buf, obstimes, &
                                 & pHorizInterpStruct, ifHorizInterp, &
                                 & pVertInterpStruct, ifVertInterp, &
                                 & dispMarketPtr, dispersion_grid, dispersion_vertical)
      if(error)return
      call arrange_supermarket(dispMarketPtr, .true., .true.)
      if(error)return
      nDiagnosed = nDiagnosed + 1
    endif  ! cell_size_z_flag
    
    !
    ! Vd correction DMAT
    !
    if(fu_quantity_in_list(Vd_correction_DMAT_flag, dqListDyn))then
      call msg('Diagnosing Vd correction')
      call check_obstimes(Vd_correction_DMAT_flag)
      call df_DMAT_Vd_correction(met_src_missing, met_buf, obstimes, &
                               & ifHorizInterp, pHorizInterpStruct, &
                               & dispMarketPtr)
      if(error)return
!      call arrange_supermarket(dispMarketPtr, .false., .true.)
!      if(error)return
      nDiagnosed = nDiagnosed + 1
    endif  ! Vd_correction_DMAT_flag
    
    !
    ! dispersion cell flux. Smart routine, has all interpolations inside, etc
    ! Diagnoses several quantities but presence of just one is enough for calling
    !
    if(fu_quantity_in_list(disp_flux_celleast_flag, dqListDyn))then !, disp_flux_cellnorth_flag, disp_flux_celltop_flag,  &
                                                                    !       & disp_cell_airmass_flag)
      call check_obstimes(disp_flux_celleast_flag)
      call diag_cell_fluxes(met_src_missing, met_buf, obstimes, &
                          & dispMarketPtr, dispersion_grid, dispersion_vertical, diagnostic_rules)
      if(error)return
      call arrange_supermarket(dispMarketPtr, .true., .true.)
      if(error)return
      nDiagnosed = nDiagnosed + 4
    endif   ! cell fluxes

    !
    ! Make daily cumulative quantities. Each cumulative quantity must have one and only one 
    ! initial quantity. We know a few, check those. If more is needed, include it in the array
    !
    do iQ = 1, size(day_accum_quantities)
      shopQ = day_accum_quantities(iQ)
      if(fu_quantity_in_list(shopQ, dqListDyn))then
        call df_cumul_daily_variable(shopQ, fu_rate_quantity_for_acc(shopQ), fu_accumulation_type(shopQ), &
                                     & met_buf, dispMarketPtr, &
                                     & diagnostic_rules%ifAllowColdstartDailyVars)
        if(error)return
        call arrange_supermarket(dispMarketPtr, .true., .true.)
        if(error)return
        nDiagnosed = nDiagnosed + 1
      endif
    end do  ! known cumulative quantities

    !    call msg('End of diagnosing:' + fu_name(MarketPtr))

    contains

                                  
      !======================================================================

      subroutine check_obstimes(quantity)
        !
        ! Returns the times to be diagnosed
        !
        implicit none
      
        ! Imported parameters
        integer, intent(in), optional :: quantity
        
        ! Local variables
        integer :: iTmp
      
        iTmp = 0
        if(present(quantity))then
          !
          ! Quantity is given, can check if it has already been diagnosed
          !
          obstimes(1:2) = time_missing
          !
          ! day_temperature_2m_acc_flag and other cumulative ones have all the intelligence inside:
          ! future has to be rediagnozed after past initialized
          !
          ! past
          if(.not. fu_field_in_sm(dispMarketPtr, & ! Past already done before?
                                & met_src_missing, & ! metSrcsMultitime(iMetSrc), &
                                & quantity, & !fu_quantity(shopVar), &
                                & met_buf%time_past, &
                                & level_missing, &
                                & fu_multi_level_quantity(quantity), & !.true., &
                                & permanent=.false.))then
            iTmp = iTmp + 1
            obstimes(iTmp) = met_buf%time_past
          endif
          !
          ! future
          !
          if(.not. fu_field_in_sm(dispMarketPtr, & ! Future already done before?
                                & met_src_missing, & ! metSrcsMultitime(iMetSrc), &
                                & quantity, & !fu_quantity(shopVar), &
                                & met_buf%time_future, &
                                & level_missing, &
                                & fu_multi_level_quantity(quantity), & !.true., &
                                & permanent=.false.))then
            iTmp = iTmp + 1
            obstimes(iTmp) = met_buf%time_future
          endif
        else
          !
          ! No quantity - take all times in the meteo_ptr
          !
          obstimes(1) = met_buf%time_past
          obstimes(2) = met_buf%time_future
        endif  ! if quantity is given
        
        if(error)then
          call msg_warning('Something went wrong with detection of already-done dynamic fields but continue')
          call unset_error('make_diagnostic_fields')
        endif
        if(iTmp==0) nDiagnosed = nDiagnosed + 1  ! even if doing nothing, claim the credits: quantity is diagnosed
    
      end subroutine check_obstimes
    
  end subroutine make_diagnostic_fields


  !******************************************************************

  logical function fu_ifDiagnosticQuantity(qIn, qDiag)
    !
    ! Checks if qDiag can be diagnosed directly from qIn
    !
    implicit none

    integer, intent(in) :: qIn, qDiag

    if(qIn == qDiag)then        ! Simple things first
      fu_ifDiagnosticQuantity = .true.
      return
    endif
    !
    ! Remember that it may happen so that one quantity can be diagnosed from the other but
    ! not necessarily vise versa. 
    !
    select case(qIn)
      case(concentration_flag)
        fu_ifDiagnosticQuantity = (qDiag == mass_in_air_flag)
      case(mass_in_air_flag)
        fu_ifDiagnosticQuantity = (qDiag == concentration_flag)
      case(volume_mixing_ratio_flag)
        fu_ifDiagnosticQuantity = (qDiag == concentration_flag .or. qDiag == mass_in_air_flag)
      case default
        fu_ifDiagnosticQuantity = .false.
    end select

  end function fu_ifDiagnosticQuantity

  
  !***************************************************************************************************
  
  logical function fu_if_full_meteo_vertical_needed(diagnostic_rules)
    !
    ! If wind disgnostic goes from bottom to top and does not involve the global Poison-type closure
    ! we might not need the full vertical to be read from meteo file
    !
    implicit none

    type(Tdiagnostic_rules), intent(in) :: diagnostic_rules

    select case(diagnostic_rules%continuity_equation)
      case(test_wind, incompressible, incompressible_v2, anelastic_v2, from_nwp_omega)
        fu_if_full_meteo_vertical_needed = .false.
      case(hybrid_top_down)
        fu_if_full_meteo_vertical_needed = .true.
      case default
        call msg("diagnostic_rules%continuity_equation", diagnostic_rules%continuity_equation)
        call set_error('Unknown wind diagnostic rule','fu_if_full_meteo_vertical_needed')
        fu_if_full_meteo_vertical_needed = .true.
    end select
  end function fu_if_full_meteo_vertical_needed

  !***************************************************************************************************

  subroutine df_DMAT_Vd_correction(met_src, met_buf, &
                                 & obstimes, &
                                 & ifHorizInterp, interpCoefMet2DispHoriz, &
                                 & dispMarketPtr)
    !
    ! Computes the Vd correction factor following the DMAT procedure. 
    !
    implicit none

    ! Imported parameters
    TYPE(Tfield_buffer), POINTER :: met_buf
    type(meteo_data_source), intent(in) :: met_src
    type(silja_time), dimension(:), intent(in) :: obstimes
    type(THorizInterpStruct), pointer :: interpCoefMet2DispHoriz
    logical, intent(in) :: ifHorizInterp
    type(mini_market_of_stacks), pointer :: dispMarketPtr

    ! Local variables
    INTEGER :: iQ, ix, iy, iDisp, iMeteo, iTime
    INTEGER, DIMENSION(:), POINTER :: mdl_in_q
    real, DIMENSION(:), POINTER :: pLandFraction, pTotPrecipPast, pTotPrecipFuture, pVdCorrection, &
                                 & pTempr2mPast, pTempr2mFuture, pFricVelPast, pFricVelFuture
    real :: fPrec, fU_star, weight_past
    TYPE(Tfield_buffer), POINTER :: met_bufPtr
    type(silja_field_id) :: idTmp

    mdl_in_q => met_buf%buffer_quantities
    met_bufPtr => met_buf

    !
    ! Basic: land fraction in Meteo grid
    !
    pLandFraction => fu_grid_data(fraction_of_land_fld)

    !
    ! Dry deposition needs also precipitation field(s) and 2m temperature
    !
    ! First, precipitation field(s)
    !
    iQ = fu_index(mdl_in_q, total_precipitation_int_flag)
    if(iQ <= 0)then
      call set_error('No total precipitation rate','df_DMAT_Vd_correction')
      return
    endif
    pTotPrecipPast => met_buf%p2d(iQ)%past%ptr
    pTotPrecipFuture => met_buf%p2d(iQ)%future%ptr

    ! Now - 2m temperature and specific humidity
    !
    iQ = fu_index(mdl_in_q, temperature_2m_flag)
    if(iQ <= 0)then
      call set_error('No 2m temperature','df_DMAT_Vd_correction')
      return
    endif
    pTempr2mPast => met_buf%p2d(iQ)%past%ptr
    pTempr2mFuture => met_buf%p2d(iQ)%future%ptr

    !
    ! Dry deposition requires a correction function via the u* 
    !
    iQ = fu_index(mdl_in_q, friction_velocity_flag)
    if(iQ <= 0)then
      call set_error('No friction velocity in buffer','df_DMAT_Vd_correction')
      return
    endif
    pFricVelPast => met_buf%p2d(iQ)%past%ptr
    pFricVelFuture => met_buf%p2d(iQ)%future%ptr

    !
    ! Finally, reserve the space for the Vd correction itself
    !
    pVdCorrection => fu_work_array()

    !
    ! Now we make the computations of the Vd correction
    ! Scavenging will be taken as it is in SILAM standard routine
    !
    do iTime = 1, size(obstimes)
      if(.not.defined(obstimes(iTime)))exit
      weight_past = (fu_valid_time(met_buf%p2d(iQ)%future%idPtr) - obstimes(iTime)) / &
                  & (fu_valid_time(met_buf%p2d(iQ)%future%idPtr) - &
                                       & fu_valid_time(met_buf%p2d(iQ)%past%idPtr))
      if(weight_past < 0. .or. weight_past > 1.)then
        call set_error('Incompatible times (obstimes(i), meteo_past, meteo_future):' + &
                     & fu_str(obstimes(iTime)) + ',' + &
                     & fu_str(fu_valid_time(met_buf%p2d(iQ)%past%idPtr)) + ',' + &
                     & fu_str(fu_valid_time(met_buf%p2d(iQ)%future%idPtr)), &
                     & 'df_DMAT_Vd_correction')
        return
      endif

      do iy = 1, ny_dispersion
        do ix = 1, nx_dispersion
         iDisp = ix + (iy-1)*nx_dispersion
         iMeteo = fu_grid_index(nx_meteo, ix, iy, interpCoefMet2DispHoriz)
         if(iMeteo < 1 .or. iMeteo > fs_meteo)then
           call msg('ix,real(iy) in dispersion grid',ix,real(iy))
           call msg('iMeteo: ',iMeteo)
           call set_error('Strange meteo index','df_DMAT_Vd_correction')
           return
         endif
         if(error)return
         !
         ! Precipitation amount
         !
         fPrec = pTotPrecipPast(iMeteo) * weight_past + pTotPrecipFuture(iMeteo) * (1.-weight_past)
         fU_star = pFricVelPast(iMeteo) * weight_past + pFricVelFuture(iMeteo) * (1.-weight_past)
         !
         ! Scavenging and Vd correction depend on temperature 
         !
         if(pTempr2mPast(iMeteo) * weight_past + pTempr2mFuture(iMeteo) * (1.-weight_past) > 271.)then  ! Above -2 C
           !
           ! Vd correction over land and sea is different
           !
           if(pLandFraction(iMeteo) < 0.1)then
             !
             ! Over sea the Charnock formula is OK
             !
             pVdCorrection(iDisp) = 3.+ 5. * fU_star * fU_star
           else
             !
             ! Over land the surface moisture plays its role
             !
             if(fPrec < 1.e-10)then                     ! No precip.
               pVdCorrection(iDisp) = 1.
             elseif(fPrec < 3.e-4)then                  ! less than 1mm/hr
               pVdCorrection(iDisp) = 2.
             else
               pVdCorrection(iDisp) = 3.
             endif    ! Precipitation amount
           endif  ! fraction of land

         else   ! Below -2 C - all frozen

           pVdCorrection(iDisp) = 1.

         endif  ! Tempr at 2m is higher/lower than -2 C

         if(pVdCorrection(iDisp) < 0.)then
           call msg('Negative VdCorrection,iDisp,VdCorr:', iDisp, pVdCorrection(iDisp))
           call set_error('Negative VdCorrection','df_DMAT_Vd_correction')
         endif

        end do  ! ix_dispersion
      end do  ! iy_dispersion

    end do  ! obstimes

  end subroutine df_DMAT_Vd_correction


  !************************************************************************************

  subroutine df_make_cell_size_z_dyn(met_src, met_buf, obstimes, &
                                   & p_horiz_interp_struct, if_horiz_interp, &
                                   & p_vert_interp_struct, if_vert_interp, &
                                   & miniMarket, gridTarget, vertTarget)
    !
    ! Generic subroutine. Makes the cell thickness dz, whether dynamic or static and
    ! puts it to the given miniMarket. 
    ! Warning!  Should be called only for output market until proper
    ! interpolation of dz from dispersion to output vertical is made. 
    ! Then diag_cell_fluxes should be used and only for dispersion market 
    !
    implicit none
    type(meteo_data_source), intent(in) :: met_src
    type(tfield_buffer), intent(in) :: met_buf
    type(silja_time), dimension(:), intent(in) :: obstimes
    type(THorizInterpStruct), pointer :: p_horiz_interp_struct
    type(TVertInterpStruct), pointer :: p_vert_interp_struct   ! contains target grid and vertical
    logical, intent(in) :: if_horiz_interp, if_vert_interp
    type(mini_market_of_stacks), intent(inout) :: miniMarket
    type(silam_vertical), intent(in) :: vertTarget
    type(silja_grid), intent(in) :: gridTarget

    ! Local variables
    integer :: ind_time, ilev, i1d, ind_z, ind_ps, nlevs, ind_pres, nx, ny, iCellTo, iCoef, iLevMet, iTmp, ixTo, iyTo, nx_met
    type(silja_time) :: now 
    type(silja_field_id) :: idtmp
    real, dimension(:), pointer ::  surf_pres, z_met, ps_met
    type (field_3d_data_ptr), pointer :: z_met3d
    real, dimension(max_levels) :: a_half, b_half
    type (TrealPtr), dimension (1:max_levels) ::  dz_3d_ptr
    real :: pTop, prev_p_c_met, prev_z_c_met, zbot, fTmp, psTmp, ztmp, p_c_met

    ind_ps = fu_index(met_buf, surface_pressure_flag)
    ind_z  = fu_index(met_buf, height_flag)
    if (any((/ind_ps,ind_z/) < 1 )) then
       call msg ("ind_ps,ind_z", (/ind_ps,ind_z/))
      call set_error('Surface pressure or height not found', 'df_make_cell_size_z_dyn')
      return
    end if

    nlevs = fu_NbrOfLevels(vertTarget)
    call grid_dimensions(gridTarget, nx, ny)
!    call msg("Grid size requested for dz", nx*ny)
!    call report(gridTarget)
    call grid_dimensions( fu_grid(met_buf%p2d(ind_ps)%past%idPtr),nx_met,iTmp)



    if (fu_leveltype(vertTarget) == layer_btw_2_height) then
      do ind_time = 1, size(obstimes)
        now = obstimes(ind_time)
        if (.not. defined(now)) exit
        do ilev = 1, nlevs
           !ever valid
          idTmp = fu_set_field_id(met_src, &                                            
                                & cell_size_z_flag, &                                  
                                & really_far_in_past, &                                      
                                & very_long_interval, &                                    
                                & gridTarget, &
                                & fu_level(vertTarget, iLev))
          call find_field_data_storage_2d(miniMarket, idTmp, single_time_stack_flag, dz_3d_ptr(iLev)%ptr)
          dz_3d_ptr(iLev)%ptr(1:nx*ny) = fu_layer_thickness_m(fu_level(vertTarget, iLev))
        enddo
        return
      enddo
    end if

    !For hybrid levels -- more complicated

    call hybrid_coefs(vertTarget, a_half=a_half, b_half=b_half)

    !$OMP PARALLEL default(none),  shared(obstimes, nlevs, idTmp, &
    !$OMP    &    met_src,gridTarget, vertTarget, miniMarket, dz_3d_ptr, surf_pres, &
    !$OMP    &    ps_met, z_met3d, ny,nx, nx_met, met_buf, ind_ps, ind_z, error, &
    !$OMP    &  if_horiz_interp, p_horiz_interp_struct, nz_meteo, a_half_disp, &
    !$OMP    &  b_half_disp, a_met, b_met ) &
    !$OMP & private(ind_time, now,  iyTo, ixTo, iCellTo, psTmp, iCoef, iLev,&
    !$OMP & zbot, prev_z_c_met,prev_p_c_met,iLevMet,z_met,zTmp, itmp, p_c_met, &
    !$OMP & ptop, ftmp)

    do ind_time = 1, size(obstimes)
      now = obstimes(ind_time)
      if (.not. defined(now)) exit
      !$OMP MASTER
      do ilev = 1, nlevs
        idTmp = fu_set_field_id(met_src, &                                            
                              & cell_size_z_flag, &                                  
                              & now, &                                      
                              & zero_interval, &                                    
                              & gridTarget, &
                              & fu_level(vertTarget, iLev))
        call find_field_data_storage_2d(miniMarket, idTmp, multi_time_stack_flag, dz_3d_ptr(iLev)%ptr)
      enddo

      ! For hybrid we need some meteo
      if (now == met_buf%time_past) then
        ps_met  => met_buf%p2d(ind_ps)%past%ptr
        z_met3d => met_buf%p4d(ind_z)%past  
      else if (now == met_buf%time_future) then
        ps_met  => met_buf%p2d(ind_ps)%future%ptr
        z_met3d => met_buf%p4d(ind_z)%future  
      else
        call set_error('Strange time:' // fu_str(obstimes(ind_time)), &
                     & 'df_make_cell_size_z_dyn')
      end if

      idTmp = fu_set_field_id(met_src, &                                            
                            & surface_pressure_flag, &
                            & now, &                                      
                            & zero_interval, &                                    
                            & gridTarget, &                                    
                            & surface_level)
      call find_field_data_storage_2d(miniMarket, idTmp, multi_time_stack_flag, surf_pres)
      !$OMP END MASTER
      !$OMP BARRIER
      if (error) cycle

      !$OMP DO COLLAPSE(2)
      do iyTo = 1, ny  ! 
        do ixTo = 1, nx
          if (error) cycle
          iCellTo = ixTo + nx*(iyTo-1)
          
          ! Surf pressure for grid centers
          if (if_horiz_interp) then
                  psTmp = 0.
                  do iCoef = 1, p_horiz_interp_struct%nCoefs
                    iTmp =  p_horiz_interp_struct%indX(iCoef,ixTo,iYto) + (p_horiz_interp_struct%indY(iCoef,ixTo,iYto)-1)*nx_met
                    psTmp = psTmp + p_horiz_interp_struct%weight(iCoef,ixTo,iYto) * ps_met(iTmp)
                  enddo
          else 
                  psTmp = ps_met(iCellTo)
          endif
          surf_pres(iCellTo) = psTmp !Save as surface pressure

          ! Cellsize
          iLev = 1
          zbot=0
          zTmp = 0. ! pretend-to be prevois meteo height
          p_c_met = psTmp ! and corresponding pressure

          do iLevMet = 1, nz_meteo + 1
            if (iLevMet <= nz_meteo) then !Can interpolate
               prev_z_c_met = zTmp
               prev_p_c_met = p_c_met
               p_c_met = a_met(iLevMet)+b_met(iLevMet)*psTmp
               if (if_horiz_interp) then
                    z_met => z_met3d%p2d(iLevMet)%ptr
                    zTmp = 0.
                    do iCoef = 1, p_horiz_interp_struct%nCoefs
                        iTmp =  p_horiz_interp_struct%indX(iCoef,ixTo,iYto) + (p_horiz_interp_struct%indY(iCoef,ixTo,iYto)-1)*nx_met
                        zTmp = zTmp + p_horiz_interp_struct%weight(iCoef,ixTo,iYto) *  z_met(iTmp)
                    enddo
               else
                   zTmp =  z_met3d%p2d(iLevMet)%ptr(iCellTo)
               endif
            else!! If we are at the top of the atmosphere
               if (prev_p_c_met  > 2.5 * p_c_met) then !! TOA: 
                 fTmp = zTmp !Save it
                 zTmp = 2*zTmp - prev_z_c_met ! Add full previous layer (approximately)
                 prev_z_c_met = fTmp
                 prev_p_c_met = p_c_met
                 p_c_met = 0. ! Top of the atmosphere
               else
                 call set_error("Meteo vertical too shallow: upper met level below the dispersion top!",&
                 & "df_make_cell_size_z_dyn")
                 exit
               endif   
            endif
            do iLev = iLev, nlevs
                pTop = a_half_disp(iLev+1)+b_half_disp(iLev+1)*psTmp !Disp. layer top perssure
                if (pTop <= p_c_met) exit ! Can't interpolate to this level 
                        fTmp = (pTop-p_c_met)/(prev_p_c_met - p_c_met) !Weight-prev
                        fTmp = zTmp*(1.-fTmp)+ prev_z_c_met*fTmp !zTop
                        dz_3d_ptr(iLev)%ptr(iCellTo) = fTmp - zbot
                        zbot = fTmp !Previous
            enddo! iLevDisp
            if (iLev ==  nlevs + 1) exit ! All disp. levels done
           enddo!iLevMet
          if (iLevMet > nz_meteo + 1) then
                 call set_error("Dispersion top beyond top of the atmosphere!",&
                 & "df_make_cell_size_z_dyn")
          endif   
        enddo !ix
      enddo !iy
      !$OMP END DO


    end do ! time
    !$OMP END PARALLEL

  end subroutine df_make_cell_size_z_dyn


  !**************************************************************************************************

  subroutine diag_cell_fluxes(met_src, met_buf, obstimes, &
                            & dispMarketPtr, gridTarget, vertTarget, diag_rules)
    !
    ! Diagnozes mass fluxes through cell borders, as they appear
    ! from meteo, diagnoses vertical wind assuming zero surface pressure tendency
    ! if hardtop -- Poisson is solved to adjust horizontal winds
    !
    implicit none
    type(meteo_data_source), intent(in) :: met_src
    type(tfield_buffer), intent(in) :: met_buf
    type(silja_time), dimension(:), intent(in) :: obstimes
    type(mini_market_of_stacks), intent(inout) :: dispMarketPtr
    type(silam_vertical), intent(in) :: vertTarget
    type(silja_grid), intent(in) :: gridTarget
    type(Tdiagnostic_rules), intent(in) :: diag_rules
    character (len=*), parameter :: sub_name = "diag_cell_fluxes"

    ! Local variables
    integer :: ind_time, i1d, ind_ps, ind_u, ind_v, ind_z, sendcnt, rcvcnt
    integer :: ixTo, iyTo, iLev, iLevMet, iCell, iTmp, jTmp, iBound, iCellTo, iCoef, itime, last_vert_index
    type(silja_time) :: now, analysis_time
    real ::  da, db,  psTmp, zTmp, uTmp, vTmp, mTmp, cellsizeX, cellsizeY
    real :: p_top_met, p_c_met, dp_met, prevtop
    real(r8k) :: fTmp, fTmp1, inflow, absinflow, domainmass, inflowfactor
    type(silja_field_id) :: idtmp
    type (field_3d_data_ptr), pointer :: u_met3d, v_met3d,  w_met3d, z_met3d
    TYPE(silja_3d_field), POINTER ::  dz_disp3d
    REAL, DIMENSION(:), POINTER :: ps_met, u_met, v_met, z_met 
    type(silja_grid)  :: u_grid, v_grid  ! Properly staggered grids for U and V
    type(silja_grid)  :: u_grid_met, v_grid_met, grid_met
    type(THorizInterpStruct), pointer ::  pHorizInterpUu, pHorizInterpUc, pHorizInterpVv, pHorizInterpVc, &
                                          & pHorizInterpC, pHorizInterpUv, pHorizInterpVu
    real, dimension(:), pointer :: uPtr, vPtr, uLptr, vLptr, wLptr, mLptr,  colmass, fPtr, fPtr1, surf_pres    
    real, dimension(:), pointer ::  fWAptr, fBoundWAptr  !Work arrays
    real, dimension(:,:), pointer :: uPtr2D, vPtr2D, mPtr2D, tmpPtr2D, uLPtr2D, vLPtr2D, wLPtr2D,mLPtr2D, fPtr2D

    integer :: nx_met_u, nx_met, nx_disp, ny_disp,  fs_disp, nz_disp ! Grids dimensions
    integer :: offx, offy, gnx, gny
!    type(TVertInterpStruct), pointer :: celltop_interp_struct ! Not really needed yet...

!    TYPE(silja_3d_field), POINTER :: u_flux_3d, v_flux_3d, w_flux_3d, cell_mass_3d
    TYPE(silja_field), POINTER :: fldPtr
    type (TrealPtr), dimension (1:max_levels) :: uFlux_3d_ptr, vFlux_3d_ptr, wFlux_3d_ptr, cellMass_3d_ptr, airdens_3d_ptr, dz_3d_ptr
    logical :: ifRotate, ifZlevels
    logical :: ifLonGlobal, ifPoisson, ifBottomUp

    ! Fluxes and masses at the grid boarder stuff
    type (TRealPtr), dimension(4) :: FMBflat  !Slices of a work array
    real, dimension(:,:,:,:), pointer :: FMBN,FMBE,FMBW,FMBS !Fluxes and masses at the grid boarder, 
    !just views of above (nxy, nz, iMass/iFlux, our/their)
    integer, dimension(4) :: boundLength  !length of boundaries
    character(len=20) ::  EndDiagStr ! "top" or "bottom"
    logical :: ifGetWind  ! Do we need to extract wind from meteo? (not needed for test winds)
    integer :: my_y_coord, my_x_coord


    call msg ("Diagnozing cells mass fluxes from meteo for:"+fu_name(dispMarketPtr))

    !
    ! Verticals
    !
    nz_disp = fu_NbrOfLevels(vertTarget)
    select case (fu_leveltype(vertTarget))
        case (layer_btw_2_height) 
             ifZlevels = .True.
        case (layer_btw_2_sigma, layer_btw_2_hybrid)
             ifZlevels = .False.
        case default
             call set_error("Strange target  vertical","here")
    end select

    ! Met-buf indices

    ind_u = fu_index(met_buf%buffer_quantities, u_flag)
    if (fu_fails(ind_u > 0, 'u-wind not available', sub_name)) return

    ind_v = fu_index(met_buf%buffer_quantities, v_flag)
    if (fu_fails(ind_v > 0, 'v-wind not available', sub_name)) return
    
    ind_ps = fu_index(met_buf%buffer_quantities, surface_pressure_flag)
    if (fu_fails(ind_ps > 0, 'surface pressure not available', sub_name)) return

!    ind_z = int_missing
!    if (ifZlevels) then
    ind_z = fu_index(met_buf%buffer_quantities, height_flag)
    if (fu_fails(ind_z > 0, 'met height not available', sub_name)) return
!    endif

    if(error)return



    !
    ! Grids and horizontal interpolations
    !
    
    ! Need also dispersion_grid etc grids. If change this, change also cell size defs below.
    if (fu_fails(gridTarget == dispersion_grid, &
               & 'Not implememted for non-dispersion grid', sub_name)) return
    
    grid_met   = fu_grid(met_buf%p2d(ind_ps)%past%idPtr)
    u_grid_met = fu_grid(met_buf%p4d(ind_u)%past%p2d(1)%idPtr)
    v_grid_met = fu_grid(met_buf%p4d(ind_v)%past%p2d(1)%idPtr)
    call grid_dimensions(u_grid_met,nx_met_u,iTmp)
    call grid_dimensions(grid_met,nx_met,iTmp)

    u_grid =  fu_staggered_grid("u_grid_stgrd", gridTarget, .true. , .false.)
    v_grid =  fu_staggered_grid("v_grid_stgrd", gridTarget, .false., .true. )

    call  smpi_get_decomposition(nx_disp, ny_disp, offx, offy, gnx, gny)
    call smpi_get_process_xy_topology(my_y_coord,my_x_coord, iTmp, jTmp)

    fs_disp=nx_disp*ny_disp
    ifLonGlobal = fu_ifLonGlobal(gridTarget)

    pHorizInterpUu => fu_horiz_interp_struct(u_grid_met, u_grid, linear, .true.)
    pHorizInterpUc => fu_horiz_interp_struct(grid_met, u_grid, linear, .true.)
    pHorizInterpVv => fu_horiz_interp_struct(v_grid_met, v_grid, linear, .true.)
    pHorizInterpVc => fu_horiz_interp_struct(grid_met, v_grid, linear, .true.)
    pHorizInterpC => fu_horiz_interp_struct(grid_met, dispersion_grid, linear, .true.)

    if (fu_pole(grid_met) == fu_pole(gridTarget)) then
       ifRotate = .False.
       call msg("NO rotation for massfluxes")
       pHorizInterpUv => null()
       pHorizInterpVu => null()
    else
       call msg("Using rotation for massfluxes")
       ifRotate = .true.
       pHorizInterpUv => fu_horiz_interp_struct(v_grid_met, u_grid, linear, .true.)
       pHorizInterpVu => fu_horiz_interp_struct(u_grid_met, v_grid, linear, .true.)
    endif


    do itime = 1, size(obstimes)
      now = obstimes(itime)
      if (.not. defined(now)) exit
     
      if (now == met_buf%time_past) then
          u_met3d => met_buf%p4d(ind_u)%past
          v_met3d => met_buf%p4d(ind_v)%past
          ps_met =>  met_buf%p2d(ind_ps)%past%ptr
          z_met3d => met_buf%p4d(ind_z)%past
      else if (now == met_buf%time_future) then
          u_met3d => met_buf%p4d(ind_u)%future
          v_met3d => met_buf%p4d(ind_v)%future
          ps_met =>  met_buf%p2d(ind_ps)%future%ptr
          z_met3d => met_buf%p4d(ind_z)%future  
      else
        call set_error('Strange time:'+fu_str(obstimes(itime)),sub_name)
        return
      end if
      if(error)return

     !create fields
     !Pressure
     idTmp = fu_set_field_id(met_src, &                                            
             & surface_pressure_flag, &
             & now, &                                      
             & zero_interval, &                                    
             & dispersion_grid, &                                    
             & surface_level)
      call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, surf_pres)
      if (error) return

      !dZ Must be already there
      if (ifZlevels) then
          dz_disp3d => fu_sm_simple_3d_field(dispMarketPtr, met_src,&
                                       & cell_size_z_flag, &
                                       & single_time_stack_flag)
      else
          dz_disp3d => fu_sm_obstime_3d_field(dispMarketPtr, met_src,&
                                       & cell_size_z_flag, &
                                       & now,&
                                       & single_time)
      endif
      if (error) return
      do ilev = 1, nz_disp

        dz_3d_ptr(iLev)%ptr => fu_grid_data_from_3d(dz_disp3d, iLev)
         

        !rho
        idTmp = fu_set_field_id(met_src, &                                            
             & air_density_flag, &
             & now, &                                      
             & zero_interval, &                                    
             & dispersion_grid, &                                    
             & fu_level(vertTarget, iLev))
        call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, airdens_3d_ptr(iLev)%ptr)
        if (error) return

        !m
        call set_quantity(idTmp, disp_cell_airmass_flag)
        call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, cellMass_3d_ptr(iLev)%ptr)
        if (error) return
        cellMass_3d_ptr(iLev)%ptr = 0.

        !w
        call set_quantity(idTmp,disp_flux_celltop_flag)
        call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, wFlux_3d_ptr(iLev)%ptr)
        if (error) return
         wFlux_3d_ptr(iLev)%ptr = 0.

        !u
        call set_quantity(idTmp,disp_flux_celleast_flag)
        call set_grid(idTmp, u_grid)
        call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, uFlux_3d_ptr(iLev)%ptr)
        if (error) return
        uFlux_3d_ptr(iLev)%ptr = 0.
        
        !v
        call set_quantity(idTmp,disp_flux_cellnorth_flag)
        call set_grid(idTmp, v_grid)
        call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, vFlux_3d_ptr(iLev)%ptr)
        if (error) return
        vFlux_3d_ptr(iLev)%ptr = 0.
      enddo



      ! Prepare for opentop, hardtop diagnostics
      ! or just set cellmass for test_wind and return
      ifGetWind = .True. !Will be reset if no wind needed from meteo, i.e. fot test_wind
      select case (diag_rules%wind_method)
         case  (opentop)
            ifPoisson = .false.
            ifBottomUp = .true.
            EndDiagStr = "top"
            call msg(fu_str(now)+": Diagnozing cell fluxes with opentop, aka bottom-up....")

         case  (hardtop, hardtop_weighted)
            ifPoisson = .True.
            ifBottomUp = .False.
            call msg(fu_str(now)+": Diagnozing cell fluxes with Poisson....")
            EndDiagStr = "bottom"
         case (topdown)
            ifPoisson = .false.
            ifBottomUp = .false.
            EndDiagStr = "bottom"
            call msg(fu_str(now)+": Diagnozing cell fluxes with topdown....")
            
         case (test_wind)
            ! All non-realtime massfluxes are set to zero already
            ! Set  cells airmasses to cell area in m^2 at all levels
            if (diag_rules%iTestWindType == 0) then 
               !real air density is needed for settling 
               call msg(fu_str(now)+":Forcing zero cell fluxes and cell masses from air density....")
               ifGetWind = .False.
            else
               if (fu_fails(ifZlevels,"Only z-coordiantes are allowed for this wind dype", sub_name )) return
               call msg(fu_str(now)+":Forcing zero cell fluxes and 1kg/m2 cell masses....")
               fPtr =>  cellMass_3d_ptr(1)%ptr
               do iyTo = 1, ny_disp  ! 
                   iCellTo = 1 + nx_disp*(iyTo-1)
                   fPtr(iCellTo:iCellTo+nx_disp-1) = &
                     & disp_grid_size_x(1:nx_disp, iyto) * disp_grid_size_y(1:nx_disp, iyto) !FIXME
               enddo
               do iTmp=2,nz_disp
                  cellMass_3d_ptr(iTmp)%ptr(1:fs_disp) = fPtr(1:fs_disp)
               enddo
               do iTmp=1,nz_disp
                  fTmp = disp_layer_top_m(iTmp)-disp_layer_top_m(iTmp-1)
                  dz_3d_ptr(iTmp)%ptr = fTmp
                  airdens_3d_ptr(iTmp)%ptr(:) = 1.
               enddo
               return !Done for now
            endif

         case default
          call msg ("diag_rules%wind_method=",diag_rules%wind_method)
          call set_error("Unknown massflux diagnostic method",sub_name)
          return
         end select
      
      !$OMP PARALLEL default(none), shared(my_x_coord,my_y_coord,&
      !$OMP & ny_disp, nx_disp,nx_met, pHorizInterpUc, pHorizInterpUu,pHorizInterpUv, pHorizInterpVc, &
      !$OMP & pHorizInterpVu, pHorizInterpVv, phorizinterpc, u_met3d,v_met3d, z_met3d, ifZlevels, uFlux_3d_ptr, vFlux_3d_ptr, &
      !$OMP & a_half_disp, b_half_disp, nz_disp, disp_layer_top_m, airdens_3d_ptr, cellMass_3d_ptr, dz_3d_ptr, &
      !$OMP & ifGetWind, u_grid, v_grid, a_met, b_met, ifrotate, nz_meteo, nx_met_u, & 
      !$OMP & nx_dispersion, disp_grid_size_x, disp_grid_size_y, ps_met, surf_pres, error, &
      !$OMP & meteo_vertical, dispersion_vertical &
      !$OMP & ), private( &
      !$OMP & cellSizeY, cellSizeX, iCellTo, iXto, iyTo, psTmp, iCoef,&
      !$OMP & u_met, v_met, z_met, p_top_met, prevtop, last_vert_index,p_c_met, dp_met, uTmp, vTmp,&
      !$OMP & iLevMet, iLev, iTmp,  ztmp,  mTmp &
      !$OMP & &
      !$OMP &)

      if (ifGetWind) then
         ! U stuff:  Interpolate all to x-staggered grid
         !$OMP DO COLLAPSE(2)
         do iyTo = 1, ny_disp
           do ixTo = 1, nx_disp+1
             ! note: for anyGrids both y and x size can depend on y and x.
             cellSizeY = fu_dy_cell_m(u_grid,ixTo,iyTo)
             iCellTo = iXto + (nx_disp+1)*(iyTo-1)
             
             ! Surf pressure for U
             psTmp = 0.
             do iCoef = 1, pHorizInterpUc%nCoefs
               iTmp =  pHorizInterpUc%indX(iCoef,ixTo,iYto) + (pHorizInterpUc%indY(iCoef,ixTo,iYto)-1)*nx_met
               psTmp = psTmp + pHorizInterpUc%weight(iCoef,ixTo,iYto) * ps_met(iTmp)
             enddo
             
             p_top_met = psTmp
             prevtop = 0.
             last_vert_index = 1
             do iLevMet = 1, nz_meteo
               p_c_met = a_met(iLevMet)+b_met(iLevMet)*psTmp
               dp_met =  2*(p_top_met - p_c_met)
               p_top_met = p_top_met - dp_met
               
               u_met => u_met3d%p2d(iLevMet)%ptr
               uTmp = 0.
               do iCoef = 1, pHorizInterpUu%nCoefs !Umet-> u_disp
                 iTmp =  pHorizInterpUu%indX(iCoef,ixTo,iYto) + (pHorizInterpUu%indY(iCoef,ixTo,iYto)-1)*nx_met_u
                 uTmp = uTmp + pHorizInterpUu%weight(iCoef, ixTo, iYto) *  u_met(iTmp)
               enddo
               
               if (ifRotate) then !Vmet-> u_disp 
                 v_met => v_met3d%p2d(iLevMet)%ptr
                 vTmp = 0.
                 do iCoef = 1, pHorizInterpUv%nCoefs
                   iTmp=pHorizInterpUv%indX(iCoef,ixTo,iYto) + (pHorizInterpUv%indY(iCoef,ixTo,iYto)-1)*nx_met
                   vTmp = vTmp + pHorizInterpUv%weight(iCoef,ixTo,iYto) * v_met(iTmp)
                 enddo
                 uTmp = uTmp*pHorizInterpUu%rotation(1,1,ixTo,iyTo) + vTmp*pHorizInterpUu%rotation(1,2,ixTo,iyTo)
               endif
               
               ! Actual mass flux in kg/s in meteo layer
               uTmp =  uTmp * dp_met /g * CellSizeY
               
               ! Distribute U mass-flux between dispersion layers
               if (ifZlevels)  then 
                 z_met => z_met3d%p2d(iLevMet)%ptr
                 zTmp = 0.
                 do iCoef = 1, pHorizInterpUc%nCoefs
                   iTmp =  pHorizInterpUc%indX(iCoef,ixTo,iYto) + (pHorizInterpUc%indY(iCoef,ixTo,iYto)-1)*nx_met
                   zTmp = zTmp + pHorizInterpUc%weight(iCoef,ixTo,iYto) *  z_met(iTmp)
                 enddo
                 call distributePieceZ(uTmp, uFlux_3d_ptr, iCellTo, last_vert_index, &
                                     & prevtop,  zTmp, disp_layer_top_m, nz_disp)
               else
                 !         call msg("Before redist", p_top_met, last_vert_index)
                 call distributePieceHyb(uTmp, uFlux_3d_ptr, iCellTo, last_vert_index,  &
                                       & dp_met, p_top_met, a_half_disp, b_half_disp, nz_disp, psTmp)
               endif
             end do ! iLevMet
           end do ! ixTo
         end do ! iyTo
         !$OMP END DO

         ! V stuff:  Interpolate all to y-staggered grid
         !$OMP DO COLLAPSE(2)
         do iyTo = 1, ny_disp+1  ! 
           do ixTo = 1, nx_disp
             cellsizeX = fu_dx_cell_m(v_grid,ixTo,iyTo)
             iCellTo = ixTo + nx_dispersion*(iyTo-1)
             
             ! Surf pressure for U
             psTmp = 0.
             do iCoef = 1, pHorizInterpVc%nCoefs
               iTmp =  pHorizInterpVc%indX(iCoef,ixTo,iYto) + (pHorizInterpVc%indY(iCoef,ixTo,iYto)-1)*nx_met
               psTmp = psTmp + pHorizInterpVc%weight(iCoef,ixTo,iYto) * ps_met(iTmp)
             enddo
             
             p_top_met = psTmp
             prevtop = 0.
             last_vert_index = 1
             do iLevMet = 1, nz_meteo
               p_c_met = a_met(iLevMet)+b_met(iLevMet)*psTmp
               dp_met =  2*(p_top_met - p_c_met)
               p_top_met = p_top_met - dp_met
               
               v_met => v_met3d%p2d(iLevMet)%ptr
               vTmp = 0.
               do iCoef = 1, pHorizInterpVv%nCoefs !Umet-> u_disp
                 iTmp =  pHorizInterpVv%indX(iCoef,ixTo,iYto) + (pHorizInterpVv%indY(iCoef,ixTo,iYto)-1)*nx_met
                 vTmp = vTmp + pHorizInterpVv%weight(iCoef, ixTo, iYto) *  v_met(iTmp)
               enddo
               
               if (ifRotate) then !Vmet-> u_disp 
                 u_met => u_met3d%p2d(iLevMet)%ptr
                 uTmp = 0.
                 do iCoef = 1, pHorizInterpVu%nCoefs
                   iTmp=pHorizInterpVu%indX(iCoef,ixTo,iYto) + (pHorizInterpVu%indY(iCoef,ixTo,iYto)-1)*nx_met_u
                   uTmp = uTmp + pHorizInterpVv%weight(iCoef,ixTo,iYto) * u_met(iTmp)
                 enddo
                 vTmp = uTmp*pHorizInterpVv%rotation(2,1,ixTo,iyTo) + vTmp*pHorizInterpVv%rotation(2,2,ixTo,iyTo)
               endif
               
               ! Actual mass flux in kg/s in meteo layer
               vTmp =  vTmp * dp_met /g * CellSizeX
               !                 vTmp = 0
               !                 if (modulo(iYto,10) ==0) vTmp = 1.
               
               ! Distribute U mass-flux between dispersion layers
               if (ifZlevels)  then 
                 z_met => z_met3d%p2d(iLevMet)%ptr
                 zTmp = 0.
                 do iCoef = 1, pHorizInterpVc%nCoefs
                   iTmp =  pHorizInterpVc%indX(iCoef,ixTo,iYto) + (pHorizInterpVc%indY(iCoef,ixTo,iYto)-1)*nx_met
                   zTmp = zTmp + pHorizInterpVc%weight(iCoef,ixTo,iYto) *  z_met(iTmp)
                 enddo
                 call distributePieceZ(vTmp, vFlux_3d_ptr, iCellTo, last_vert_index, &
                                     & prevtop,  zTmp, disp_layer_top_m, nz_disp)
               else
                 call distributePieceHyb(vTmp, vFlux_3d_ptr, iCellTo, last_vert_index,  &
                                       & dp_met, p_top_met, a_half_disp, b_half_disp, nz_disp, psTmp)
               endif
             end do ! iLevTo
           end do ! ixTo
         end do ! iyTo
        !$OMP END DO
      endif


      ! Mass, density and surface pressure stuff:  Interpolate all to non-staggered grid

      !$OMP DO COLLAPSE(2)
      do iyTo = 1, ny_disp  ! 
ix:      do ixTo = 1, nx_disp
          if (error) cycle
          cellsizeX = disp_grid_size_x(ixto, iyto)
          cellsizeY = disp_grid_size_y(ixto, iyto)
          iCellTo = ixTo + nx_disp*(iyTo-1)
          
          ! Surf pressure for grid centers
          psTmp = 0.
          do iCoef = 1, pHorizInterpC%nCoefs
            iTmp =  pHorizInterpC%indX(iCoef,ixTo,iYto) + (pHorizInterpC%indY(iCoef,ixTo,iYto)-1)*nx_met
            psTmp = psTmp + pHorizInterpC%weight(iCoef,ixTo,iYto) * ps_met(iTmp)
          enddo
          surf_pres(iCellTo) = psTmp !Save as surface pressure
          
          p_top_met = psTmp
          prevtop = 0.
          last_vert_index = 1
          do iLevMet = 1, nz_meteo
            p_c_met = a_met(iLevMet)+b_met(iLevMet)*psTmp
            dp_met =  2*(p_top_met - p_c_met)
            p_top_met = p_top_met - dp_met
            
            ! Actual mass in kg in meteo layer
            mTmp =  dp_met * CellSizeX * CellSizeY / g
            
            ! Distribute mass between dispersion layers
            if (ifZlevels)  then 
              z_met => z_met3d%p2d(iLevMet)%ptr
              zTmp = 0.
              do iCoef = 1, pHorizInterpC%nCoefs
                iTmp =  pHorizInterpC%indX(iCoef,ixTo,iYto) + (pHorizInterpC%indY(iCoef,ixTo,iYto)-1)*nx_met
                zTmp = zTmp + pHorizInterpC%weight(iCoef,ixTo,iYto) *  z_met(iTmp)
              enddo
              if (iLevMet == nz_meteo) then
                 if (zTmp <  disp_layer_top_m(nz_disp)) then
                   !$OMP CRITICAL (bark_diag_cell_fluxes)
                   call msg("")
                   call msg("")
                   call msg("Dispersion vertical sticks out of meteo")
                   call msg("Zmet_top, zdisp_top",  zTmp,  disp_layer_top_m(nz_disp))
                   call msg("Dispersion ix, iy", ixTo, iyTo)
                   call msg('Meteo vertical')
                   call report(meteo_vertical, .true.)
                   call msg('dispersion vertical')
                   call report(dispersion_vertical, .true.)
                   do jTmp = 1, nz_meteo
                     call msg('Meteo level and its height:',jTmp,z_met3d%p2d(jTmp)%ptr(iTmp))
                   end do
                   call set_error("Dispersion  z-vertical is not covered by meteo", sub_name)
                   !$OMP END CRITICAL (bark_diag_cell_fluxes)
                   cycle ix
                endif
              endif
              call distributePieceZ(mTmp, cellMass_3d_ptr, iCellTo, last_vert_index, &
                                  & prevtop,  zTmp, disp_layer_top_m, nz_disp)
            else
              ! This could be done through a and b coefficients
              call distributePieceHyb(mTmp, cellMass_3d_ptr, iCellTo, last_vert_index,  &
                                    & dp_met, p_top_met, a_half_disp, b_half_disp, nz_disp, psTmp)
            endif
          end do ! iLevMet

          !Air density Could have been done with vector operations
          do iLev = 1, nz_disp
             airdens_3d_ptr(iLev)%ptr(iCellTo) = cellMass_3d_ptr(iLev)%ptr(iCellTo) / (dz_3d_ptr(iLev)%ptr(iCellTo) * CellSizeX * CellSizeY )
          enddo
        end do ix! ixTo
      end do ! iyTo
      !$OMP END DO
      !$OMP END PARALLEL


      if (.not. ifGetWind .or. error) return ! No further diagnostics is needed
    
      if (ifPoisson .or. smpi_is_mpi_version()) then
         ! Need to ensure that domains have exactly the same idea on their boundary vlocities
         ! and behind-boundary masses. Same stuff is used for boundary masses in poisson.
         boundLength(1:4) = (/nx_disp, nx_disp, ny_disp, ny_disp/)
         fBoundWAptr => fu_work_array((2*nx_disp+ 2*ny_disp)*4*nz_disp)
         do iTmp = 1,4
            jTmp = sum(boundLength(1:iTmp-1))*4*nz_disp
            FMBflat(iTmp)%ptr(1:boundLength(iTmp)*4*nz_disp) => &
               & fBoundWAptr(jTmp+1:jTmp+boundLength(iTmp)*4*nz_disp)
         enddo
         FMBN(1:nx_disp, 1:nz_disp, 1:2, 1:2) =>  FMBflat(northern_boundary)%ptr(:)
         FMBS(1:nx_disp, 1:nz_disp, 1:2, 1:2) =>  FMBflat(southern_boundary)%ptr(:)
         FMBE(1:ny_disp, 1:nz_disp, 1:2, 1:2) =>  FMBflat(eastern_boundary)%ptr(:)
         FMBW(1:ny_disp, 1:nz_disp, 1:2, 1:2) =>  FMBflat(western_boundary)%ptr(:)

         ! Fill stuff from behind-boundaries for transmission
!         if (ifLonGlobal) call msg("Global")
         do iLev = 1,nz_disp
               uPtr2D(1:nx_disp+1,1:ny_disp)   => uFlux_3d_ptr(iLev)%ptr(1:fs_disp+ny_disp)
               vPtr2D(1:nx_disp, 1:ny_disp+1) => vFlux_3d_ptr(iLev)%ptr(1:fs_disp+nx_disp)
               mPtr2D(1:nx_disp, 1:ny_disp)    => cellMass_3d_ptr(iLev)%ptr(1:fs_disp) 
               FMBN(:,iLev,iMass,our) = mPtr2D(:,ny_disp)
               FMBN(:,iLev,iFlux,our) = vPtr2D(:,ny_disp+1)
               FMBS(:,iLev,iMass,our) = mPtr2D(:,1)
               FMBS(:,iLev,iFlux,our) = vPtr2D(:,1)

               if (ifLonGlobal) then 
                  !average with opposite side
                  FMBE(:,iLev,iMass,our) = 0.5 * (mPtr2D(nx_disp,:) + mPtr2D(1,:))
                  FMBE(:,iLev,iFlux,our) = 0.5 * (uPtr2D(nx_disp+1,:) + uPtr2D(1,:))
                  FMBW(:,iLev,iMass,our) =  FMBE(:,iLev,iMass,our) 
                  FMBW(:,iLev,iFlux,our) =  FMBE(:,iLev,iFlux,our) 
               else !from this side
                  FMBW(:,iLev,iMass,our) = mPtr2D(1,:)
                  FMBW(:,iLev,iFlux,our) = uPtr2D(1,:)
                  FMBE(:,iLev,iMass,our) = mPtr2D(nx_disp,:)
                  FMBE(:,iLev,iFlux,our) = uPtr2D(nx_disp+1,:)
               endif
         enddo

         ! Exchange Flux and Mass behind  and average with neighbour 
         do iTmp=1,4 !Boundaries
            if (mod(my_y_coord+my_x_coord,2) == 0) then !Should be different order for the neighbour
                !2,1,4,3
                iBound = iTmp + 2*mod(iTmp,2) - 1 ! iTmp-1 for even, iTmp+1 for odd iTmp
            else
                iBound = iTmp ! 1,2,3,4
            endif
            if(adv_mpi_neighbours(iBound)>=0) then
               call msg(fu_str(smpi_adv_rank)+"Exchanging FMB with neighbor, boundary:"+fu_str(iBound))
               call msg("boundLength(iTmp)", boundLength(iTmp), nz_disp)
               sendcnt = boundLength(iTmp)*nz_disp*2
               fPtr => FMBflat(iBound)%ptr(1:sendcnt)  !!! Our
               fPtr1 => FMBflat(iBound)%ptr(sendcnt+1:2*sendcnt) !!!Their
               flush(run_log_funit)

               call smpi_exchange_wings(adv_mpi_neighbours(iBound), &
                       & fPtr, sendcnt,  fPtr1, rcvcnt)
               if (sendcnt /= rcvcnt) then
                 call msg("after smpi_exchange_wings sendcnt,recvcnt", sendcnt,rcvcnt )
                 call set_error("sendcnt and recvcnt do not match", sub_name)
                 return
               endif

               fPtr =  0.5*(fPtr + fPtr1)  !! Take a compromise as valid for this boundary
            endif
         enddo


         !Put mean boundary massfluxes back
         do iLev = 1,nz_disp
               uPtr2D(1:nx_disp+1,1:ny_disp)   => uFlux_3d_ptr(iLev)%ptr(1:fs_disp+ny_disp)
               vPtr2D(1:nx_disp,1:ny_disp+1)   => vFlux_3d_ptr(iLev)%ptr(1:fs_disp+nx_disp)

               uPtr2D(nx_disp+1,:) = FMBE(:,iLev,iFlux,our)
               uPtr2D(1,:)         = FMBW(:,iLev,iFlux,our)

               vPtr2D(:,ny_disp+1) = FMBN(:,iLev,iFlux,our) 
               vPtr2D(:,1)         = FMBS(:,iLev,iFlux,our)
         enddo
      endif
      !! Masses in FMBX(:,:,iMass,our) are ready to use as masses at the domain boundary



      !Done with getting U_flux, v_flux and mass from meteo. Now let us adjust U_flux and v_flux,
      ! so the flux is non-divergent, diagnoze w_flux. 
      ! Pressure tendency (mass change) is ignored here completely. realtime correction will account for that.
      ! Now we deal strictly with one time step. Goal -- non-divergent massfluxin the domain

      if (ifPoisson)  then
         fWAptr => fu_work_array(3*fs_disp+nx_disp+ny_disp) 
         uPtr(1:fs_disp+ny_disp) => fWAptr(1:fs_disp+ny_disp)
         iTmp = fs_disp+ny_disp+1
         vPtr(1:fs_disp+nx_disp) => fWAptr(iTmp:iTmp+fs_disp+nx_disp-1)
         iTmp = iTmp + fs_disp+nx_disp
         colmass(1:fs_disp) => fWAptr(iTmp:iTmp+fs_disp-1)

         uPtr(1:fs_disp+ny_disp) = 0.
         vPtr(1:fs_disp+nx_disp) = 0.
         colmass(1:fs_disp) = 0.
        
         !integtate masses and fluxes over vertical
         do iLev = 1,nz_disp
           colmass(1:fs_disp) = colmass(1:fs_disp) + cellMass_3d_ptr(iLev)%ptr(1:fs_disp)  * Disp_Wind%LevelCorrWeight(iLev)
           iTmp = ny_disp*(nx_disp+1) 
           uPtr(1:iTmp) =       uPtr(1:iTmp) +     uFlux_3d_ptr(iLev)%ptr(1:iTmp)
           iTmp = nx_disp*(ny_disp+1) 
           vPtr(1:iTmp) =       vPtr(1:iTmp) +     vFlux_3d_ptr(iLev)%ptr(1:iTmp)
         enddo

         if (gnx == nx_disp .and.  gny == ny_disp) then ! All in one domain
            uPtr2D(1:nx_disp+1,1:ny_disp)   => uPtr(1:fs_disp+ny_disp)
            vPtr2D(1:nx_disp,1:ny_disp+1)   => vPtr(1:fs_disp+nx_disp)
            mPtr2D(1:nx_disp,1:ny_disp)     => colmass(1:fs_disp) 
         else
            uLptr => fu_work_array((gnx+1)*gny)
            uPtr2D(1:gnx+1,1:gny) => uLptr(1:(gnx+1)*gny)
            call smpi_gather_field(uPtr, uPtr2D , .true., .false.)

            vLptr => fu_work_array((gny+1)*gnx)
            vPtr2D(1:gnx,1:gny+1) => vLptr(1:gnx*(gny+1))
            call smpi_gather_field(vPtr, vPtr2D, .false., .true.)
            
            mLptr => fu_work_array(gnx*gny)
            mPtr2D(1:gnx,1:gny) => mLptr(1:gnx*gny)
            call smpi_gather_field(colmass, mPtr2D , .false., .false.)
         endif

         if (smpi_adv_rank == 0) then ! Master
           call adjust_2D_fluxes(uPtr2D, vPtr2D, null(), mPtr2D, &
                      & fu_ifLonGlobal(wholeMPIdispersion_grid), &
                      & fu_ifPolarCapPossible(wholeMPIdispersion_grid, southern_boundary), &
                      & fu_ifPolarCapPossible(wholeMPIdispersion_grid, northern_boundary), &
                      & gnx, gny, fTmp)  
         endif

         !Make uPtr2D, vPtr2D, and mPtr2D look to local domain
         if (gnx /= nx_disp .or.  gny /= ny_disp) then
            call smpi_bcast_aray(uLptr, (gnx+1)*gny, smpi_adv_comm,0)
            fPtr2D(1:gnx+1,1:gny) => uLptr(1:(gnx+1)*gny)
            uPtr2D(1:nx_disp+1,1:ny_disp)   => uPtr(1:fs_disp+ny_disp)
            uPtr2D(1:nx_disp+1,1:ny_disp) = fPtr2D(offx+1:offx+nx_disp+1,offy+1:offy+ny_disp)

            call smpi_bcast_aray(vLptr, gnx*(gny+1),smpi_adv_comm,0)
            fPtr2D(1:gnx,1:gny+1) => vLptr(1:(gny+1)*gnx)
            vPtr2D(1:nx_disp,1:ny_disp+1) => vPtr(1:fs_disp+nx_disp)
            vPtr2D(1:nx_disp,1:ny_disp+1) = fPtr2D(offx+1:offx+nx_disp,offy+1:offy+ny_disp+1)

            ! This did not change 
            mPtr2D(1:nx_disp,1:ny_disp)  => colmass(1:fs_disp)

            !Cleanup
            call free_work_array(uLptr)
            call free_work_array(vLptr)
            call free_work_array(mLptr)
         endif 



         !Apply corrections to  U and V fluxes at levels
         do iLev = 1,nz_disp
               fTmp =  Disp_Wind%LevelCorrWeight(iLev)
               uLptr2D(1:nx_disp+1,1:ny_disp) => uFlux_3d_ptr(iLev)%ptr(:)
               vLptr2D(1:nx_disp,1:ny_disp+1) => vFlux_3d_ptr(iLev)%ptr(:)
               mLptr2D(1:nx_disp,1:ny_disp) => cellMass_3d_ptr(iLev)%ptr(:)

               uLptr2D(2:nx_disp,1:ny_disp) = uLptr2D(2:nx_disp,1:ny_disp) - &
                   & uPtr2D(2:nx_disp,1:ny_disp)* &
                   & 0.5 * fTmp * (mLptr2D(2:nx_disp,1:ny_disp)+ mLptr2D(1:nx_disp-1,1:ny_disp))
               uLptr2D(1,1:ny_disp) = uLptr2D(1,1:ny_disp) - &
                     & uPtr2D(1,1:ny_disp) * fTmp * FMBW(:,iLev,iMass,our) 
               uLptr2D(nx_disp+1,1:ny_disp) = uLptr2D(nx_disp+1,1:ny_disp) - &
                     & uPtr2D(nx_disp+1,1:ny_disp)  * fTmp * FMBE(:,iLev,iMass,our)

               vLptr2D(1:nx_disp,2:ny_disp) =  vLptr2D(1:nx_disp,2:ny_disp) - &
                     & vPtr2D(1:nx_disp,2:ny_disp) * &
                     & 0.5 * fTmp * (mLptr2D(1:nx_disp,2:ny_disp)+ mLptr2D(1:nx_disp,1:ny_disp-1))  
               vLptr2D(1:nx_disp,1) =  vLptr2D(1:nx_disp,1) - &
                     & vPtr2D(1:nx_disp,1) * fTmp * FMBS(:,iLev,iMass,our)
               vLptr2D(1:nx_disp,ny_disp+1) =  vLptr2D(1:nx_disp,ny_disp+1) - &
                     & vPtr2D(1:nx_disp,ny_disp+1) * fTmp * FMBN(:,iLev,iMass,our)
         enddo
         ! Check it
!         call msg("*************************CHECKUP******************")
!         colmass(1:fs_disp) = 0.
!         uPtr(1:fs_disp+ny_disp) = 0.
!         vPtr(1:fs_disp+nx_disp) = 0.
!        
!         !integtate masses and fluxes over vertical
!         do iLev = 1,nz_disp
!           colmass(1:fs_disp) = colmass(1:fs_disp) + cellMass_3d_ptr(iLev)%ptr(1:fs_disp)  
!           iTmp = ny_disp*(nx_disp+1) 
!           uPtr(1:iTmp) =  uPtr(1:iTmp) + uFlux_3d_ptr(iLev)%ptr(1:iTmp)
!           iTmp = nx_disp*(ny_disp+1) 
!           vPtr(1:iTmp) =  vPtr(1:iTmp) + vFlux_3d_ptr(iLev)%ptr(1:iTmp)
!         enddo
!
!            uPtr2D(1:nx_disp+1,1:ny_disp)   => uPtr(1:fs_disp+ny_disp)
!            vPtr2D(1:nx_disp,1:ny_disp+1)   => vPtr(1:fs_disp+nx_disp)
!            mPtr2D(1:nx_disp,1:ny_disp)     => colmass(1:fs_disp) 
!           call adjust_2D_fluxes(uPtr2D, vPtr2D, null(), mPtr2D, norm2D, &
!                      & fu_ifLonGlobal(wholeMPIdispersion_grid), &
!                       & fu_ifPolarCapPossible(wholeMPIdispersion_grid, southern_boundary), &
!                       & fu_ifPolarCapPossible(wholeMPIdispersion_grid, northern_boundary), &
!                      & gnx, gny, fTmp)  
!         call msg("*************************CHECKUP OVER******************")
         
        call free_work_array(fWaPtr)

      endif ! Poisson


      if (ifPoisson .or. smpi_is_mpi_version()) call free_work_array(fBoundWAptr)

      fTmp = 0.
      fTmp1 = 0.
      ! Bottom-up diagnostics for w_flux
      if (ifBottomUp) then
         call msg("Doing bottom-up")
         tmpPtr2D(1:nx_disp,1:ny_disp) => wFlux_3d_ptr(1)%ptr(:)
         tmpPtr2D(1:nx_disp,1:ny_disp) = 0. !No flux from below
         do iLev = 1,nz_disp
             uLptr2D(1:nx_disp+1,1:ny_disp) => uFlux_3d_ptr(iLev)%ptr(:)
             vLptr2D(1:nx_disp,1:ny_disp+1) => vFlux_3d_ptr(iLev)%ptr(:)
             wLptr2D(1:nx_disp,1:ny_disp)   => wFlux_3d_ptr(iLev)%ptr(:)

             wLPtr2D(1:nx_disp,1:ny_disp) = tmpPtr2D(1:nx_disp,1:ny_disp)  & !Flux from below
                  & +  vLptr2D(1:nx_disp,1:ny_disp)  -  vLptr2D(1:nx_disp,2:ny_disp+1) & !Y-convergence
                  & +  uLptr2D(1:nx_disp,1:ny_disp)  -  uLptr2D(2:nx_disp+1,1:ny_disp) !X-convergence

             tmpPtr2D(1:nx_disp,1:ny_disp) => wFlux_3d_ptr(iLev)%ptr(:)
             fTmp = ftmp - sum(vLptr2D(1:nx_disp,1)) !flux off Norhern boundary
             fTmp1 = ftmp1 + sum(vLptr2D(1:nx_disp,ny_disp+1)) !flux off Southern boundary

         enddo
      else !TopDown
         call msg("Doing top-down")
         tmpPtr2D(1:nx_disp,1:ny_disp) => wFlux_3d_ptr(nz_disp)%ptr(:)
         tmpPtr2D(1:nx_disp,1:ny_disp) = 0. !No flux from above
         do iLev = nz_disp, 2, -1 !No w field below the domain
             uLptr2D(1:nx_disp+1,1:ny_disp) => uFlux_3d_ptr(iLev)%ptr(:)
             vLptr2D(1:nx_disp,1:ny_disp+1) => vFlux_3d_ptr(iLev)%ptr(:)
             wLptr2D(1:nx_disp,1:ny_disp)   => wFlux_3d_ptr(iLev-1)%ptr(:)

             wLPtr2D(1:nx_disp,1:ny_disp) = tmpPtr2D(1:nx_disp,1:ny_disp)  & !Flux from below
                  & -  vLptr2D(1:nx_disp,1:ny_disp)  +  vLptr2D(1:nx_disp,2:ny_disp+1) & !Y-convergence
                  & -  uLptr2D(1:nx_disp,1:ny_disp)  +  uLptr2D(2:nx_disp+1,1:ny_disp) !X-convergence

             tmpPtr2D(1:nx_disp,1:ny_disp) => wFlux_3d_ptr(iLev-1)%ptr(:)
             fTmp = ftmp - sum(vLptr2D(1:nx_disp,1)) !flux off Norhern boundary
             fTmp1 = ftmp1 + sum(vLptr2D(1:nx_disp,ny_disp+1)) !flux off Southern boundary

         enddo
             uLptr2D(1:nx_disp+1,1:ny_disp) => uFlux_3d_ptr(1)%ptr(:)
             vLptr2D(1:nx_disp,1:ny_disp+1) => vFlux_3d_ptr(1)%ptr(:)
             fTmp = ftmp - sum(vLptr2D(1:nx_disp,1)) !flux off Norhern boundary
             fTmp1 = ftmp1 + sum(vLptr2D(1:nx_disp,ny_disp+1)) !flux off Southern boundary

      endif

      call msg("Sum of vertical flux at the domain_"+EndDiagStr+"_kg/s",  sum(tmpPtr2D))
      call msg("Sum of abs vertical flux at the domain_"+EndDiagStr+"_kg/s",  sum(abs(tmpPtr2D)))
      call msg("Sum of upward flux at domain_"+EndDiagStr, sum(tmpPtr2D, mask=(tmpPtr2D.gt.0))) 
      call msg("Max vertical cell flux at the domain_"+EndDiagStr+"_kg/s",  maxval(tmpPtr2D))
      call msg("Min vertical cell flux at the domain_"+EndDiagStr+"_kg/s",  minval(tmpPtr2D))
      call msg("Mean abs vertical mass flux at the domain_"+EndDiagStr+"_kg/(m2*s)", &
            &  sum(abs(tmpPtr2D))/(4*pi*earth_radius*earth_radius))
      call msg("Total massflux at S boundary (pole divergence)",  fTmp)
      call msg("Total massflux at N boundary (pole divergence)",  fTmp1)

 
  enddo ! Time loop
 

  contains

    !====================================================================

    subroutine distributePieceZ(Piece, targetPtr, iCellTo, lastiZ, metltop, metZ, z_top_d, nz_disp)
            ! Distribute piece between Z dispersion layers
            implicit none
            real, intent(in) :: Piece, metZ   
            type (Trealptr), dimension(1:), intent(inout) :: targetPtr
            real, dimension(-1:), intent(in) :: z_top_d
            integer, intent(inout) :: lastiZ
            real,    intent(inout) :: metltop  ! Top of previous met level
            integer, intent(in) :: nz_disp, iCellTo

            real :: bl_tmp, SS, br_abs, zTmp, fTmp, frachere

            bl_tmp = metltop !Slab bottom height
            SS = 2*(metZ - metltop)    
            metltop = bl_tmp + SS !Slab top height
            fTmp = 1.0
            do lastiZ = lastiZ, nz_disp
               if (z_top_d(lastiZ) < metltop ) then
                  frachere = (z_top_d(lastiZ) - bl_tmp) / SS
                  fTmp = fTmp - frachere
                  bl_tmp = z_top_d(lastiZ)
                  targetPtr(lastiZ)%ptr(iCellTo) = &
                   & targetPtr(lastiZ)%ptr(iCellTo) + frachere * Piece  
                else
                  targetPtr(lastiZ)%ptr(iCellTo) = &
                     & targetPtr(lastiZ)%ptr(iCellTo) + fTmp * Piece 
                  exit
                endif
            enddo
    end  subroutine distributePieceZ
    !====================================================================

    subroutine distributePieceHyb(Piece, targetPtr, iCellTo, lastiZ, dp, p_top, a_top_d, b_top_d, nz_disp, ps)
            implicit none
            real, intent(in) :: Piece, dp, p_top, ps 
            type (Trealptr), dimension(1:), intent(inout) :: targetPtr
            real, dimension(-1:), intent(in) :: a_top_d, b_top_d
            integer, intent(inout) :: lastiZ
            integer, intent(in) :: nz_disp, iCellTo

            real(8) :: bl_tmp, fTmp, frachere, pTmp

            bl_tmp = p_top + dp
            ! call msg("Inslab  top_p, bottp", real(p_top), real(bl_tmp))
            fTmp = 1.0
            do lastiZ = lastiZ, nz_disp
               pTmp = a_top_d(lastiZ)  + ps * b_top_d(lastiZ) !p Disp. layertop
               if (pTmp > p_top) then
                  frachere = (bl_tmp - pTmp) / dp
                  fTmp = fTmp - frachere
                  bl_tmp = pTmp
                  targetPtr(lastiZ)%ptr(iCellTo) = &
                   & targetPtr(lastiZ)%ptr(iCellTo) + frachere * Piece
            !       call msg("iLev, frachere", lastiZ, frachere)
                else
                  targetPtr(lastiZ)%ptr(iCellTo) = &
                   & targetPtr(lastiZ)%ptr(iCellTo) + fTmp * Piece
             !      call msg("iLev, frachere", lastiZ, fTmp)
                   exit
                endif
            enddo
     end subroutine distributePieceHyb

  end subroutine diag_cell_fluxes

  
  !***************************************************************************************************
  
  subroutine df_cumul_daily_variable(quantity_cumul, quantity_inst, iAvType, &
                                   & met_buf, dispMarketPtr, ifColdstartAllowed)
    !
    ! Computes the cumulative daily value for a given quantity from its instant counterpart
    ! Instant quantity _must_ be in meteo buffer, whereas the cumulative one will be put
    ! to dispersion market. 
    ! Note that allowColdStart flag should be handled on the _second_ call
    ! As a way to detect the initialization we store the zero-accumualted
    ! past-field if allowColdStart==false.
    !
    implicit none

    ! Imported parameters
    integer, intent(in) :: quantity_cumul, quantity_inst, iAvType
    TYPE(Tfield_buffer), POINTER :: met_buf
    type(mini_market_of_stacks), intent(inout) :: dispMarketPtr
    type(silja_logical), intent(in) :: ifColdstartAllowed

    ! Local variables
    INTEGER :: iQ, ix, iy, iTo, iMeteo, iTime, fsTo, iLev, nLevs, nxTo, nyTo, iCoef, iTmp, indQinst
    real, DIMENSION(:), POINTER :: pVarPast, pVarFuture, pTCum, pTCum_past
    TYPE(Tfield_buffer), POINTER :: met_bufPtr
    type(silja_field_id) :: idTmp
    type(silja_field_id), pointer :: idVarPast, idVarFuture, idPastPtr, idFuturePtr
    type(silja_time) :: start_of_day, timeTmp
    type(silja_interval) :: intervalTmp, accLen
    logical :: ifAddPAst, ifPastMadeHere, ifQInstant_is_singletime
    real :: seconds, vPast, vFuture
    real, dimension(:), pointer :: pData
    type(silja_field), pointer :: fldPtr
    type(THorizInterpStruct), pointer :: pHorizInterpStruct
    character(len=20), parameter :: sub_name = 'df_cumul_daily_variable'

    pHorizInterpStruct => fu_horiz_interp_struct(meteo_grid, dispersion_grid, linear, .true.)
    met_bufPtr => met_buf
    call grid_dimensions(dispersion_grid, nxTo, nyTo)
    fsTo = nxTo * nyTo
    !
    indQinst =  fu_index(met_buf%buffer_quantities, quantity_inst)
    if(indQinst == int_missing)then
      call msg_warning('Needed instant quantity is not in met_buf:' + fu_quantity_string(quantity_inst), &
                   & sub_name)
      call msg('Quantities in the buffer:')
      do iQ = 1, size(met_buf%buffer_quantities)
        if(met_buf%buffer_quantities(iQ) == int_missing)exit
        call msg(fu_quantity_string(met_buf%buffer_quantities(iQ)))
      end do
      call set_error('Needed instant quantity is not in met_buf:' + fu_quantity_string(quantity_inst), &
                   & sub_name)
      return
    endif   ! indQinst is missing
    !
    ! The instant quantity can be single-time, e.g., precipitation rate. Then
    ! past and future times have to be taken from met_buf central timing and the 
    ! field value itself - from the present pointer
    !
    ifQInstant_is_singletime = fu_realtime_quantity(quantity_inst)
    if(ifQInstant_is_singletime)then
      pVarPast => met_buf%p2d(indQinst)%present%ptr
      pVarFuture => met_buf%p2d(indQinst)%present%ptr
      idVarPast => met_buf%p2d(indQinst)%present%idPtr
      idVarFuture => met_buf%p2d(indQinst)%present%idPtr
    else
      pVarPast => met_buf%p2d(indQinst)%past%ptr
      pVarFuture => met_buf%p2d(indQinst)%future%ptr
      idVarPast => met_buf%p2d(indQinst)%past%idPtr
      idVarFuture => met_buf%p2d(indQinst)%future%idPtr
    endif
    !
    !  Detect cold/warm start, otherwise past should be ready: 
    !  Cumulative quantity must never have zero accumulation length
    !
    if(ifQInstant_is_singletime)then
      timeTmp = met_buf%time_past
    else
      timeTmp = fu_valid_time(idVarPast)
    endif
   ! Should be previous day for midnight, otherwise today morning
    start_of_day = fu_start_of_day_utc(timeTmp - one_second)
    accLen = timeTmp - start_of_day !! Should never be zero!

    idTmp = fu_set_field_id_simple(fu_met_src(idVarPast), &
                                 & quantity_cumul, timeTmp, level_missing) 
    !
    ! Dangerous call - idTmp time is not guaranteed and will be reset to a close time 
    ! that is directly available from stacks. It better be checked
    !
    call get_field_from_mm_general(dispMarketPtr, idTmp, fldPtr, .false.)

    if (.not. associated(fldPtr)) then ! Coldstart
      !
      ! Create something here as mean from past and future temperature accumulated since day start
      !
      ifPastMadeHere = .true.

      idTmp = fu_set_field_id(fu_met_src(idVarPast),&
                            & quantity_cumul, &
                            & start_of_day, &      !analysis_time,&
                            & accLen, &            !forecast_length, &
                            & dispersion_grid, &
                            & fu_level(idVarPast), &
                            & accLen, & !        & !length_of_accumulation, optional
                            & zero_interval, &     !length_of_validity, optional
                            & accumulated_flag)    !field_kind, optional

      call msg_warning("Making surrogate past for:" + fu_quantity_string(quantity_cumul),sub_name)
      call msg("ID for the surrogatre field")
      call report(idTmp)

      if(error)return
      call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, pTCum_past)
      if(fu_fails(.not.error,'Failed field data pointer', sub_name))return
         
      seconds = fu_sec(accLen) 
      if (pHorizInterpStruct%ifInterpolate) then
        do iy = 1, nyTo
          do ix = 1, nxTo
            iTo = ix + (iy-1) * nxTo
            vPast = 0.0
            vFuture = 0.0
            do iCoef = 1, pHorizInterpStruct%nCoefs
              iTmp = pHorizInterpStruct%indX(iCoef,ix,iy) + &
                        & (pHorizInterpStruct%indY(iCoef,ix,iy) - 1) * nx_meteo
              vPast = vPast + pHorizInterpStruct%weight(iCoef,ix,iy) * pVarPast(iTmp)
              vFuture = vFuture + pHorizInterpStruct%weight(iCoef,ix,iy) * pVarFuture(iTmp)
            end do
            !
            ! Past and future of the instant quantity computed. Accumulate them now.
            ! Note that there was nothing in pTCum_past till now
            !
            select case(iAvType)
              case(av4mean, av4sum)
                pTCum_past(iTo) = (vPast + vFuture) * seconds * 0.5
              case(av4min)
                pTCum_past(iTo) = min(vPast, vFuture)
              case(av4max)
                pTCum_past(iTo) = max(vPast, vFuture)
              case default
                call set_error('Unknown averaging type 1:' + fu_str(iAvType), sub_name)
                return
            end select
          end do  ! ix
        end do  ! iy
      else
        select case(iAvType)
          case(av4mean, av4sum)
            pTCum_past(1:fsTo) = (pVarPast(1:fsTo) + pVarFuture(1:fsTo)) * seconds * 0.5
          case(av4min)
            pTCum_past(1:fsTo) = min(pVarPast(1:fsTo),pVarFuture(1:fsTo))
          case(av4max)
            pTCum_past(1:fsTo) = max(pVarPast(1:fsTo),pVarFuture(1:fsTo))
          case default
            call set_error('Unknown averaging type 2:' + fu_str(iAvType), sub_name)
            return
        end select
      endif  ! if horizontal interpolation
    else
      !
      !  fldPtr exists. Pick it up
      !
      idPastPtr => fu_id(fldPtr) 
      pTCum_past => fu_grid_data(fldPtr)
      intervalTmp = fu_accumulation_length(idPastPtr)
      ifPastMadeHere = .false.
      ! Should be right time
      if(.not. fu_valid_time(idPastPtr) == timeTmp)then
        call msg_warning('Wrong time found in stack',sub_name)
        call msg('Requested and provided times:' + fu_str(timeTmp) + '<==  ==>' + &
                                                 & fu_str(fu_valid_time(idPastPtr)))
        call set_error('Wrong time found in stack',sub_name)
        return
      endif

      !!Should be valid field
      if (pTCum_past(1) == real_missing) then
        if (fu_true(ifColdStartAllowed)) then
          call  set_error("This should not happen", sub_name)
          return
        elseif(fu_false(ifColdStartAllowed)) then
          call  set_error("This should not happen either", sub_name)
          return
        else
          call msg_warning("Looks like:" + fu_quantity_string(quantity_cumul) &
                         & + "-was not initialised. crashing..", sub_name)
          call msg("This behaviour can be overridden by 'allow_coldstart_daily_variables = yes' in standard_setup")
          call set_error("Missing values for past", sub_name)
          return
        endif
      endif

      ! Correct accumulation that might have been broken from the initialization
      if (.not. intervalTmp == accLen) then
          call msg_warning("Correcting acc_length for past:" + fu_quantity_string(quantity_cumul), sub_name)
          call msg("Initial ID")
          call report(idPastPtr)
          call set_accumulation_length(idPastPtr, accLen)
          call msg("Corrected ID")
          call report(idPastPtr)
      endif
      idTmp = idPastPtr
    endif  ! past cumulative field exists
    if(error)return
    !
    ! Past should be okay. Future now
    ! We should make it in any case: past might have been reset with ini
    
    ! Should be previous day for midnight
    if(ifQInstant_is_singletime)then
      timeTmp = met_buf%time_future
    else
      timeTmp = fu_valid_time(idVarFuture)
    endif
    start_of_day = fu_start_of_day_utc(timeTmp - one_second)
    intervalTmp = timeTmp - start_of_day !! Should never be zero!

    ifAddPast = (fu_hour(fu_valid_time(idTmp)) > 0) 

    call set_valid_time(idTmp, timeTmp)
    call set_accumulation_length(idTmp,intervalTmp)

    ! On the second diagnostics we should get the existig field
    !
!
!     Dangerous call - idTmp time is not guaranteed and will be reset to something
!     close and directly available, e.g. previous time step
!
!    fldPtr => fu_get_field_from_mm_general(dispMarketPtr, idTmp, .false.)
!    if (associated(fldPtr)) then
!      idFuturePtr => fu_id(fldPtr)
!      idFuturePtr = idTmp
!      pTCum => fu_grid_data(fldPtr)
!    else
!      !This guy always creates a field 
!      call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, pTCum)
!      if(fu_fails(.not.error,'Failed field data pointer', sub_name))return
!    endif

    ! safe call is this one (also creates the field if it does not exist yet):
    call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, pTCum)
    if(fu_fails(.not.error,'Failed field data pointer', sub_name))return

    pTCum(1:fsTo) =  0.
    if (pHorizInterpStruct%ifInterpolate) then
      do iy = 1, nyTo
        do ix = 1, nxTo
          iTo = ix + (iy-1) * nxTo
          vPast = 0.0
          vFuture = 0.0
          do iCoef = 1, pHorizInterpStruct%nCoefs
              iTmp =  pHorizInterpStruct%indX(iCoef,ix,iy) + &
                   & (pHorizInterpStruct%indY(iCoef,ix,iy)-1) * nx_meteo
              vPast = vPast + pHorizInterpStruct%weight(iCoef,ix,iy) * pVarPast(iTmp)
              vFuture = vFuture + pHorizInterpStruct%weight(iCoef,ix,iy) * pVarFuture(iTmp)
!              pTCum(iTo) = pTCum(iTo) + pHorizInterpStruct%weight(iCoef,ix,iy) * &
!                                      & (pVarPast(iTmp) + pVarFuture(iTmp))
          enddo  ! iCoef
          select case(iAvType)
            case(av4mean, av4sum)
              pTCum(iTo) = vPast + vFuture
            case(av4min)
              pTCum_past(iTo) = min(vPast, vFuture)
            case(av4max)
              pTCum_past(iTo) = max(vPast,vFuture)
            case default
              call set_error('Unknown averaging type 3:' + fu_str(iAvType), sub_name)
              return
          end select
        end do  ! ix
      end do  ! iy
    else
      select case(iAvType)
        case(av4mean, av4sum)
          pTCum(1:fsTo) = pVarPast(1:fsTo) + pVarFuture(1:fsTo)
        case(av4min)
          pTCum(1:fsTo) = min(pVarPast(1:fsTo),pVarFuture(1:fsTo))
        case(av4max)
          pTCum(1:fsTo) = max(pVarPast(1:fsTo),pVarFuture(1:fsTo))
        case default
          call set_error('Unknown averaging type 4:' + fu_str(iAvType), sub_name)
          return
      end select
    endif

    !! Finalize accumulation and add earlier accumulation if needed
    seconds = fu_sec(fu_valid_time(idVarFuture) - fu_valid_time(idVarPast))
    if (ifAddPast) then
      select case(iAvType)
        case(av4mean, av4sum)
          pTCum(1:fsTo) = pTCum(1:fsTo)*0.5*seconds + pTCum_past(1:fsTo) 
        case(av4min)
          pTCum(1:fsTo) = min(pTCum(1:fsTo), pTCum_past(1:fsTo))
        case(av4max)
          pTCum(1:fsTo) = max(pTCum(1:fsTo), pTCum_past(1:fsTo))
        case default
          call set_error('Unknown averaging type 5:' + fu_str(iAvType), sub_name)
          return
      end select
    else
      select case(iAvType)
        case(av4mean, av4sum)
          pTCum(1:fsTo) = pTCum(1:fsTo) * 0.5 * seconds
        case(av4min, av4max)
        case default
          call set_error('Unknown averaging type 6:' + fu_str(iAvType), sub_name)
          return
      end select
    endif
    !
    ! Guard for missing initialisation
    if (ifPastMadeHere .and. fu_false(ifColdStartAllowed))then
       pTCum_past(1:fsTo) = real_missing
    endif

  end subroutine df_cumul_daily_variable
  
  
  !***************************************************************************************************

  subroutine df_update_realtime_met_fields(pMMarket, buf_ptr, now, wdr, diag_rules, &
                                         & valid_time_border_1, valid_time_border_2) !, ifResetTimes)
    !
    ! Updates realtime fields in meteo market/buffer
    ! Uses buffer wherever possible
    !
    implicit none

    ! Imported variables
    type(mini_market_of_stacks), pointer :: pMMarket
    type(Tfield_buffer), pointer :: buf_ptr
    type(silja_time), intent(in) :: now
    type(silja_wdr), intent(in) :: wdr
    type(Tdiagnostic_rules), intent(in) :: diag_rules
    type(silja_time), intent(inout) :: valid_time_border_1, valid_time_border_2
!    logical, intent(in) :: ifResetTimes

    !Local 
    type(silja_field), pointer :: field
    type(silja_field_id), pointer :: idIn
    integer :: iVar, shopQ, nQstat
    integer, dimension(max_quantities) :: quantities
    logical :: ifFound
    integer, parameter, dimension(3) :: rain_int = (/large_scale_rain_int_flag, &
                                                   & convective_rain_int_flag, &
                                                   & total_precipitation_int_flag/)
    valid_time_border_1 = really_far_in_past
    valid_time_border_2 = really_far_in_future

    !
    ! Go over the available diagnostics checking if the specific one has been ordered
    ! Order matters!
    !
    ! First, get the single-time quantities
    !
    call supermarket_quantities(pMMarket, fu_met_src(wdr,1), single_time_stack_flag, quantities, nQstat)
    if(error .or. nQstat < 1)return
    call msg("Updating single-time fields for now =" + fu_str(now))

    !
    ! Update the ground pressure field
    !
    if(fu_quantity_in_quantities(ground_pressure_flag, buf_ptr%buffer_quantities))then
      ! only multitime ground_pressure_field is accesible from the buffer.
      ! Should use stack directly
      call find_field_from_stack(met_src_missing, ground_pressure_flag, time_missing,& ! any time 
                                 & fu_stack(pMMarket, 1),& !ST stack
                                 & field, ifFound)
      if (ifFound)then
        idIn => fu_id(field)
        !
        ! check validity and update if needed
        if(fu_do_work(buf_ptr, int_missing, idIn))then
          call msg('Updating ground pressure')
          CALL df_update_ground_pressure(field, buf_ptr, valid_time_border_1, valid_time_border_2)
        endif
      else
        call msg("Not found in stack:"+ fu_name(fu_stack(pMMarket, 1)) )
      endif   ! found in stack
    endif  ! ground pressure

    !
    ! Precipitation intensity
    !
    do iVar = 1, size(rain_int)
      shopQ = rain_int(iVar)
      if(fu_quantity_in_quantities(shopQ, buf_ptr%buffer_quantities))then
        ! update
        if(fu_do_work(buf_ptr, shopQ))then
          call msg('Updating rain:' + fu_quantity_string(shopQ))
          CALL df_update_rain_intensity(buf_ptr, shopQ, wdr, valid_time_border_1, valid_time_border_2)
        endif
      endif ! needed?
    end do  ! rain vars

    !
    ! Scavenging coefficient
    !
    if(fu_quantity_in_quantities(scavenging_coefficient_flag, buf_ptr%buffer_quantities))then
      ! update
      if(fu_do_work(buf_ptr, scavenging_coefficient_flag))then
        call msg('Updating scavenging coefficient')
        CALL df_update_scav_coefficient(buf_ptr, valid_time_border_1, valid_time_border_2)
      endif
    endif
    
    !
    ! Dispersion flux. Needed just one quantity in buffer to activate
    !
    if(fu_quantity_in_quantities(disp_flux_cellt_rt_flag, buf_ptr%buffer_quantities))then
      ! update
      if(fu_do_work(buf_ptr, disp_flux_cellt_rt_flag))then
        call msg('Updating dispersion fluxes')
        CALL df_update_cellfluxcorr_rt(buf_ptr, valid_time_border_1, valid_time_border_2, &
                                     & diag_rules, now)
        call exchange_wings(buf_ptr)
      endif
    endif

    !
    ! daily mean parameters
    !
    do iVar = 1, size(day_mean_quantities)
      shopQ = day_mean_quantities(iVar)
      if(fu_quantity_in_quantities(shopQ, buf_ptr%buffer_quantities))then
        ! update
        if(fu_do_work(buf_ptr, shopQ))then
          call msg('Updating daily mean:' + fu_quantity_string(shopQ))
          if(fu_fails(fu_get_time_direction_sm(pMMarket) == forwards, &
                    & 'Mean daily and cumulative variables cannot be made in adjoint runs', &
                    & 'df_update_realtime_met_fields'))return
          call df_daily_mean_variable(fu_cumul_quantity_for_mean_one(shopQ), shopQ, fu_accumulation_type(shopQ), &
                                    & buf_ptr, now, valid_time_border_1, valid_time_border_2, &
                                    & diag_rules%ifAllowColdstartDailyVars)
        endif
      endif ! needed?
    end do  ! daily mean parameters


    CONTAINS
    
    !=============================================================================================
    
    logical function fu_do_work(buf_ptr, shopQ, idPtr_)
      !
      ! Checks if a field is still valid
      !
      implicit none
      ! Imported parameters
      type(Tfield_buffer), pointer :: buf_ptr
      integer, intent(in) :: shopQ
      type(silja_field_id), pointer, optional :: idPtr_
      
      ! Local parameters
      integer :: iTmp, Qidx
      type(silja_field_id) :: idRequest
      type(silja_field_id), pointer :: idPtr
      
      fu_do_work = .false.

      ! get the input id
      !
      if(present(idPtr_))then
        idPtr => idPtr_
      else
        Qidx  = fu_index(buf_ptr, shopQ) 
        if (Qidx < 1) return
        if (fu_dimension(buf_ptr, Qidx) == 4 ) then 
          idPtr => buf_ptr%p4d(Qidx)%present%p2d(1)%idPtr
        else
          idPtr => buf_ptr%p2d(Qidx)%present%idPtr
        endif
      endif
      !
      ! ID found. Does the field cover the required time?
      !
      idRequest = fu_set_field_id_simple(met_src_missing, shopQ, now, level_missing)
      if (error) return

      if (fu_field_id_covers_request(idPtr, idRequest, .true.)) then
        call msg("Still valid:"+ fu_quantity_string(fu_quantity(idPtr)))
        if(valid_time_border_1 < fu_valid_time(idPtr)) valid_time_border_1 = fu_valid_time(idPtr)
        if(valid_time_border_2 > fu_valid_time(idPtr) + fu_validity_length(idPtr)) &
                        & valid_time_border_2 = fu_valid_time(idPtr) + fu_validity_length(idPtr)
        return
      endif
      fu_do_work = .true.
    end function fu_do_work

  end subroutine df_update_realtime_met_fields


  !*********************************************************************************************

  subroutine df_update_rain_intensity(met_buf,  desiredQ, wdr, &
                                    & valid_time_border_1, valid_time_border_2)
    !
    ! Updates rain intensity in met_buf. No checks. it just uses
    ! past and future accumulated rains to create present rain.
    ! The resulting rain is set  valid for the accumulation interval
    !
    type(Tfield_buffer), pointer :: met_buf    ! input meteo data
    integer, intent(in) :: desiredQ            ! rain quantitity wanted
    type(silja_wdr), intent(in) :: wdr
    type(silja_time), intent(inout) :: valid_time_border_1, valid_time_border_2

    integer :: accQ, indAcc, indAv, iVar
    real, dimension(:), pointer :: instant_rain
    real ::  aclen
    type(silja_field_id), pointer :: idFuture, idPast, idAv
    logical :: previous_same_start
  
    select case (desiredQ)
      case (large_scale_rain_int_flag)
        accQ =  large_scale_accum_rain_flag
      case (convective_rain_int_flag)
        accQ = convective_accum_rain_flag
      case (total_precipitation_int_flag)
        accQ = total_precipitation_acc_flag
      case default
        call set_error("Can't make rain out of:"+ fu_quantity_string(desiredQ), &
                       "df_update_rain_intensity")
        return
    end select
    
!call msg("Filling field: "//trim(fu_quantity_string(desiredQ)))

    indAv = fu_index(met_buf%buffer_quantities, desiredQ)
    instant_rain => met_buf%p2d(indAv)%present%ptr
    idAv => met_buf%p2d(indAv)%present%IDptr

    indAcc = fu_index(met_buf%buffer_quantities, accQ)
    if (indAcc < 1) then
      call msg("")
      call msg("Could not find accumulated rain in meteo buf...")
      call report_list_of_quantities(met_buf%buffer_quantities,"Buffer contents")
      call set_error("Failed to make rain intensity","df_update_rain_intensity")
      return
    endif
    idFuture => met_buf%p2d(indAcc)%future%idPtr
    idPast  => met_buf%p2d(indAcc)%past%idPtr

!call report(idPast)
!call report(idFuture)

    previous_same_start = fu_accumulation_start_time(idFuture) == &
                            & fu_accumulation_start_time(idPast)

    if(previous_same_start) then 
      ! Previous field has same start of accumulation.
      aclen = fu_sec(fu_accumulation_length(idFuture) - fu_accumulation_length(idPast))

      IF (aclen > 0.) THEN
        instant_rain(1:fs_meteo) =  ( met_buf%p2d(indAcc)%future%ptr(1:fs_meteo) &
                                  &  - met_buf%p2d(indAcc)%past%ptr(1:fs_meteo) ) / aclen
      ELSE
        instant_rain(1:fs_meteo) = 0.
        call msg("**********FUTURE RAIN*************")
        CALL report(met_buf%p2d(indAcc)%future%idPtr)
        call msg("**********PAST RAIN*************")
        CALL report(met_buf%p2d(indAcc)%past%idPtr)
        CALL set_error('zero accumulation between obstimes','df_update_rain_intensity')
        RETURN
      END IF
    ELSE ! New accumulation interval

      aclen = fu_sec(fu_accumulation_length(idFuture))

      IF (aclen > 0.) THEN
        instant_rain(1:fs_meteo) = met_buf%p2d(indAcc)%future%ptr(1:fs_meteo) /aclen
      ELSE
        instant_rain(1:fs_dispersion) = 0.
        call msg("**********FUTURE RAIN*************")
        CALL report(met_buf%p2d(indAcc)%future%idPtr)
          call msg("**********PAST RAIN (not used...)*************")
        CALL report(met_buf%p2d(indAcc)%past%idPtr)
        CALL msg_warning('zero accumulation time for future', 'df_update_rain_intensity')
      END IF
    END IF  ! previous exist
!call msg("Accumulation length:", aclen)
!call msg("average rain past", sum(met_buf%p2d(indAcc)%past%ptr(1:fs_meteo))/fs_meteo)
!call msg("average rain future", sum(met_buf%p2d(indAcc)%future%ptr(1:fs_meteo))/fs_meteo)
!call msg("average rain now *aclen", sum(instant_rain(1:fs_meteo)*aclen)/fs_meteo)

!     apply minimum rain rate threshold

    where (instant_rain(1:fs_meteo) < fu_precipitation_low_limit(wdr)) instant_rain(1:fs_meteo) = 0.

!call msg("Cutoff & average rain now after cutoff, *aclen", fu_precipitation_low_limit(wdr), sum(instant_rain(1:fs_meteo)*aclen)/fs_meteo)
!call msg('')

    !
    ! Almost done. Now set proper validity in field ID
    !
    call set_valid_time(idAv, fu_valid_time(idPast))
    call set_analysis_time(idAv, fu_analysis_time(idPast))
    call set_validity_length(idAv, fu_set_interval_sec(aclen))

    !
    ! Adjust the borders of the global validity period of all single-time fields
    !
    if(valid_time_border_1 < fu_valid_time(idAv))valid_time_border_1 = fu_valid_time(idAv)
    if(valid_time_border_2 > fu_valid_time(idAv) + fu_validity_length(idAv)) &
                    & valid_time_border_2 = fu_valid_time(idAv) + fu_validity_length(idAv)

!call msg('Reporting the field')
!call report(idAv)
!call msg("df_update_rain_intensity Done!")
!call msg("")
!call msg("")
  end subroutine df_update_rain_intensity


  !***********************************************************************************
  
  subroutine df_update_scav_coefficient( met_buf,  valid_time_border_1, valid_time_border_2)
    !
    ! Updates scav_coefficient in met_buf. No checks. it just uses
    ! real-time rain intensities.
    ! The resulting rain is set  valid for the rain validity time,
    ! as dependency on instant meteo fields is assumed minor.
    ! Their values are taken for the middle of validity interval
    ! rewritten from dq_scavenging_coef (derived_field_quantities_2)
    ! Roux

    !
    ! First-guess rough estimate of the scavenging coefficient regardless the
    ! reatures of the scavenged substance. 
    !
    ! The in-cloud/sub-cloud division is determined for a 'reference' cloud,
    ! with constant cloud top and bottom heights. The water/snow
    ! phase is chosen according to the air temperature at the particle
    ! elevation.
    !     
    ! "Instatenous" HIRLAM precipitation values are used
    ! whenever they are available. Otherwise the accumulated values are taken
    !
    ! The water/snow phase are split regarding the temperature at the particle height
    ! The scavenging coefficients are computed from the NAME-II type methods
    !
    ! Author Mikhail Sofiev
    ! Implemented algorithm is taken from diy_weather_tools of I.Valkama
    !
    ! All units: SI

    type(Tfield_buffer), pointer :: met_buf    ! input meteo data
    type(silja_time), intent(inout) :: valid_time_border_1, valid_time_border_2

    integer ::   iScav,iLS, iCnv, iP, iT
    real, dimension(:), pointer :: instant_rain
    real ::  aclen
    type(silja_field_id), pointer :: idRain, idScav
    logical :: previous_same_start
    integer :: iTmp, iLev, ix
    REAL, DIMENSION(:), POINTER :: scav, pr_ls, pr_conv, t, p
    integer, parameter :: lookup_length = 1000 !index in 0.1 mm/h
    REAL, DIMENSION(0:lookup_length), save :: sc_rain, sc_snow ! lookup tables
    logical, save :: ifLookupReady = .false.
    integer, dimension(:), pointer :: ind_ls, ind_cnv
    real :: fMaxScav

    !
    ! Fill lookup table
    if (.not. ifLookupReady) then
      call msg("making lookup table for scav. coefficient....")
      do iTmp = 0, lookup_length
         !Subcloud-rain, incloud-largscale-rain, 0.25*incloud_convective-any
          sc_rain(iTmp) = (iTmp/36000.)**0.79 * 0.05417     !8.40E-5*(3600**0.79)
         !Snow, except for incloud_convective
          sc_snow(iTmp) = (iTmp/36000.)**0.305 * 9.722e-04  ! 8.00E-5* 3600.**0.305
      enddo 
      ifLookupReady = .true.
    endif

    iScav = fu_index(met_buf%buffer_quantities, scavenging_coefficient_flag)
    iLS   = fu_index(met_buf%buffer_quantities, large_scale_rain_int_flag)
    iCnv  = fu_index(met_buf%buffer_quantities, convective_rain_int_flag)

    idRain => met_buf%p2d(iLS)%present%IDptr

    !rain  as integer in 0.1 mm/h
    ind_ls =>  fu_work_int_array()
    ind_cnv =>  fu_work_int_array()

    !   Note precipitation unit: mm s-1
    pr_ls => met_buf%p2d(iLS)%present%ptr
    DO ix = 1,fs_meteo
      ind_ls(ix)  = min( nint(pr_ls(ix)*36000.),lookup_length)
      if(ind_ls(ix) < 0)then
        call set_error('Negative LS precipitation at (' + fu_str(ix) + ')' + fu_str(pr_ls(ix)),'')
        call msg('LS Precipitation:',pr_ls)
        call msg('CONV precip:',pr_conv)
      endif
    enddo
    !print *,ind_ls(1:fs_meteo)

    if (iCnv > 0) then
      pr_conv => met_buf%p2d(iCnv)%present%ptr
      DO ix = 1,fs_meteo
        ind_cnv(ix) = min( nint(pr_conv(ix)*36000.),lookup_length)
        if(ind_cnv(ix) < 0)then
          call set_error('Negative CONV precipitation at (' + fu_str(ix) + ')' + fu_str(pr_conv(ix)),'')
          call msg('LS Precipitation:',pr_ls)
          call msg('CONV precip:',pr_conv)
        endif
      enddo
    else
      ind_cnv(1:fs_meteo) = 0
    endif
   

    iT = fu_index(met_buf%buffer_quantities, temperature_flag)
    iP = fu_index(met_buf%buffer_quantities, pressure_flag)

    t => fu_work_array(fs_meteo)
    p => fu_work_array(fs_meteo)
    if(error)return

    do iLev = 1,fu_NbrOfLevels(meteo_vertical)
        t(1:fs_meteo) = 0.5* ( met_buf%p4d(iT)%past%p2d(iLev)%ptr(1:fs_meteo) &
                              + met_buf%p4d(iT)%future%p2d(iLev)%ptr(1:fs_meteo))
        p(1:fs_meteo) = 0.5* ( met_buf%p4d(iP)%past%p2d(iLev)%ptr(1:fs_meteo) &
                              + met_buf%p4d(iP)%future%p2d(iLev)%ptr(1:fs_meteo))

        scav => met_buf%p4d(iScav)%present%p2d(iLev)%ptr
        fMaxScav = 0.

!        !FIXME Fixed value!!
!        scav(1:fs_meteo) = 0.001
!        cycle


        loop_grid: DO ix = 1,fs_meteo

          !---------- Ensure no scavenging for fog or above the cloud top
          IF(ind_ls(ix)+ind_cnv(ix) < 1 .or. p(ix) < 70000.)THEN 
            scav(ix)=0.0
            CYCLE loop_grid
          END IF

          !---------- If something exists - check the in/below cloud options
          IF(p(ix) > 90000.)THEN 
            !---  below the cloud base convective and dynamic efficiencies equal
            IF (t(ix) > 273.16)THEN  !---rain
              scav(ix)= sc_rain(ind_ls(ix)) +  sc_rain(ind_cnv(ix))
              !conv_scav = pr_conv(ix)**0.79 * 0.05417  !8.40E-5*(3600**0.79)
              !ls_scav = pr_ls(ix)**0.79 * 0.05417
            ELSE   !-----snow
              scav(ix)= sc_snow(ind_ls(ix)) +  sc_snow(ind_cnv(ix))
!              conv_scav = pr_conv(ix)**0.305 * 9.722e-04  ! 8.00E-5* 3600.**0.305
!              ls_scav = pr_ls(ix)**0.305 * 9.722e-04
            END IF

          ELSE   
            !
            !------- for in-cloud scavenging (700 hPa < height < 900 hpa) 
            !
            ! 1. Convective rain efficiencies of water and snow are equal:
            ! Four times normal rain scavenging
!            conv_scav = pr_conv(ix)**0.79 * 0.2167 ! 3.36E-4*(3600**0.79)

            ! 2. Dynamic (large-scale) rain/snow efficiencies
            !
            IF(t(ix) > 273.16)THEN  !------ rain
              scav(ix)= sc_rain(ind_ls(ix)) +  4. * sc_rain(ind_cnv(ix))
 !             ls_scav = pr_ls(ix)**0.79 * 0.05417  ! 8.40E-5*(3600**0.79)

            ELSE !----------- snow
              scav(ix)=  sc_snow(ind_ls(ix)) +  4.* sc_rain(ind_cnv(ix))
!              ls_scav = pr_ls(ix)**0.305 * 9.722e-04  ! 8.00E-5* 3600.**0.305
            END IF

          END IF ! In/sub cloud scavenging

          fMaxScav = max(fMaxScav, scav(ix))

        END DO loop_grid

        !
        ! Almost done. Now set proper validity in field ID
        ! Note, that only level 1 is checked above....
        !
        idScav => met_buf%p4d(iScav)%present%p2d(iLev)%IDptr
        call set_valid_time(idScav, fu_valid_time(idRain))
        call set_analysis_time(idScav, fu_analysis_time(idRain))
        call set_validity_length(idScav, fu_validity_length(idRain))
!      call msg("")
!call msg("level, max Scavenging", iLev, fMaxScav)
!       call msg("iScav", iScav)
!      call report(met_buf%p4d(iScav)%past%p2d(iLev)%IDptr)
!      call report(met_buf%p4d(iScav)%present%p2d(iLev)%IDptr)
!      call report(met_buf%p4d(iScav)%future%p2d(iLev)%IDptr)
      
     enddo !levels

    ! fieldID3d times should be kept coherently -- otherwise arrange stack gets
    ! mad...
    call set_valid_time(met_buf%p4d(iScav)%present%field3d, fu_valid_time(idRain))
    !
    ! Adjust the borders of the global validity period of all single-time fields
    !
    if(valid_time_border_1 < fu_valid_time(idScav))valid_time_border_1 = fu_valid_time(idScav)
    if(valid_time_border_2 > fu_valid_time(idScav) + fu_validity_length(idScav)) &
                    & valid_time_border_2 = fu_valid_time(idScav) + fu_validity_length(idScav)
   
    call free_work_array(t)
    call free_work_array(p)
    call free_work_array(ind_ls)
    call free_work_array(ind_cnv)

  end subroutine df_update_scav_coefficient


  !**********************************************************************************************

  subroutine df_update_ground_pressure(field, met_buf, valid_time_border_1, valid_time_border_2) 
    !
    ! Updates ground pressure field. No checks.
    ! The resulting pressure is set valid between past and future
    ! Roux
    implicit none
    
    ! Imported parameters
    type(Tfield_buffer), pointer :: met_buf    ! input meteo data
    type(silja_time), intent(inout) :: valid_time_border_1, valid_time_border_2
    
    ! Local variables
    REAL, DIMENSION(:), POINTER ::  p
    TYPE(silja_field), POINTER :: field
    TYPE(silja_field_id), POINTER :: idFuture, idPast, idPrST
    integer :: iTmp, iPrDyn

    iPrDyn = fu_index(met_buf%buffer_quantities, ground_pressure_flag)
    if (iPrDyn < 1) then
       call msg("Failed to use "+ &
               & fu_quantity_string(ground_pressure_flag) + &
               & ". Trying " + fu_quantity_string(surface_pressure_flag))
       iPrDyn = fu_index(met_buf%buffer_quantities, surface_pressure_flag)
       if (iPrDyn < 1) then

          call msg("Buffer contains:")
          do iTmp = 1,size(met_buf%buffer_quantities)
               call msg(fu_quantity_string(met_buf%buffer_quantities(iTmp)))
          enddo
          call set_error("Could not find find_pressure_flag in met buffer", &
                  "df_update_ground_pressure")
          return
       endif
       call msg("Using " + fu_quantity_string(surface_pressure_flag))
    endif
  
    idFuture => met_buf%p2d(iPrDyn)%future%idPtr
    idPast  => met_buf%p2d(iPrDyn)%past%idPtr

    p => fu_grid_data(field)

    p(1:fs_meteo) = 0.5* ( met_buf%p2d(iPrDyn)%past%ptr(1:fs_meteo) &
                         + met_buf%p2d(iPrDyn)%future%ptr(1:fs_meteo))

    idPrST => fu_id(field)
       !
       ! Almost done. Now set proper validity in field ID
       ! Note, that only level 1 is checked above....
       !
    call set_valid_time(idPrST, fu_valid_time(idPast))
    call set_analysis_time(idPrST, fu_analysis_time(idPast))
    call set_validity_length(idPrST, fu_valid_time(idFuture) - fu_valid_time(idPast))

    !
    ! Adjust the borders of the global validity period of all single-time fields
    !
    if(valid_time_border_1 < fu_valid_time(idPrST)) valid_time_border_1 = fu_valid_time(idPrST)
    if(valid_time_border_2 > fu_valid_time(idPrST) + fu_validity_length(idPrST)) &
                    & valid_time_border_2 = fu_valid_time(idPrST) + fu_validity_length(idPrST)

!  call msg("Singletime ground_perssure after update by df_update_ground_pressure")   
!  call report(field)
!  call msg("")
!  call msg("")
!  call msg("")

  end subroutine df_update_ground_pressure


  !**********************************************************************************************

 subroutine df_update_cellfluxcorr_rt(disp_buf, valid_time_border_1, valid_time_border_2, diag_rules, now) 
   !
   ! Updates the correction for  cell fluxes according to given continuity equation
   ! Hardly relies on staggered grids
   ! The resulting fields are set valid between past and future

   ! For hardtop diagnostics -- zero mass flux at the top is aimed.
   ! For non-global domain -- by adjusting inflow from boundaries:
   ! discrepancy is distributed proportionally to mass in cells
   ! For global domain mass change rate in cells is adjusted propotionally to masses
   ! to get mass-conservative domain.

   ! The fluxes in multitime stack are assumed to be non-divergent
   ! for each time step. Here we calculate a correction for these fluxes
   ! due to mass change

   ! 
   ! Roux
   implicit none
   
   ! Imported parameters
   type(Tfield_buffer), pointer :: disp_buf  ! buffer in which we diagnoze
   type(silja_time), intent(inout) :: valid_time_border_1, valid_time_border_2
   type(silja_time), intent(in) :: now
   type(Tdiagnostic_rules), intent(in) :: diag_rules
   
   ! Local variables
   REAL, DIMENSION(:), POINTER ::  p, fPtr, fPtr1, mbPtr, maPtr, ulPtr, vlPtr, wLptr, colmass, colmassdelta!Pointers

   REAL, DIMENSION(:), POINTER ::  uPtr,vPtr, fWAptr, mLptr, dmLptr,fBoundWAptr !Work arays


   TYPE(silja_field), POINTER :: field
   integer :: iTmp, iPrDyan, nx_disp, ny_disp, iLev, jTmp, fs_disp, iY, iX, offx, offy, gnx, gny,iBound, nz_disp
   integer :: ind_m, ind_furt, ind_fvrt, ind_fwrt, i1d,  sendcnt, rcvcnt
   type (silja_time) :: valid, analysis 
   type (silja_interval) :: validity
   logical :: ifLonGlobal, ifPoisson, ifBottomUp
   real :: seconds
   real(8) ::  inflowfactor, fTmp, fTmp1, masscorrfactor
   real(8) :: domainmassdelta, deltamassfactor, domainmass
   real, dimension(:,:), pointer :: uPtr2D, vPtr2D, mPtr2D, dmPtr2d, tmpPtr2D, uLPtr2D, vLPtr2D, wLPtr2D, mbPtr2D, maPtr2d



   integer, parameter ::  T =  12*3600 ! seconds -- period 
   REAL  :: corner_lat_N,corner_lon_E, southpole_lat_N, southpole_lon_E, dx_deg, dy_deg
   INTEGER :: number_of_points_x, number_of_points_y
   LOGICAL :: corner_in_geo_lonlat
   real :: pi_t_over_T, tfactor, lon, lat, cosTheta, sin2Theta, &
             & lambda_prime, sinlprime2, sin2lprime, cosTheta2, sinTheta, cosLambda, sinLambda

   ! Fluxes and masses at the grid boarder stuff
   type (TRealPtr), dimension(4) :: MBflat  !Slices of a work array
   real, dimension(:,:,:), pointer :: MBN,MBE,MBW,MBS !Fluxes and masses at the grid boarder, 
                                       !just views of above
   integer, dimension(4) :: boundLength  !length of boundaries
   integer :: my_y_coord, my_x_coord

   ! For test winds
    TYPE(silam_pole) :: pole45 
    real :: rlon, rlat  

    character(len=19), parameter :: sub_name = 'df_update_cellfluxcorr_rt'

   !Realtime (output)
   ind_furt = fu_index(disp_buf%buffer_quantities,disp_flux_celle_rt_flag) 
   ind_fvrt = fu_index(disp_buf%buffer_quantities,disp_flux_celln_rt_flag) 
   ind_fwrt = fu_index(disp_buf%buffer_quantities,disp_flux_cellt_rt_flag) 
   ! Multitime (input)
   ind_m = fu_index(disp_buf%buffer_quantities,disp_cell_airmass_flag)
   if (any((/ind_furt, ind_fvrt, ind_m, ind_fwrt/) < 1)) then
      call msg("Failed to find quantities in buffer:")
      call msg("ind_furt, ind_fvrt, ind_m, ind_fwrt", &
          (/ind_furt, ind_fvrt, ind_m, ind_fwrt/))

      call msg("Buffer contains:")
      do iTmp = 1,size(disp_buf%buffer_quantities)
              call msg(fu_quantity_string(disp_buf%buffer_quantities(iTmp)))
      enddo
      call set_error("Could not find input/output fields in disp buffer", &
                 "df_update_cellfluxes_rt")
         return
   endif

   call smpi_get_decomposition(nx_disp, ny_disp, offx, offy, gnx, gny)
   call smpi_get_process_xy_topology(my_y_coord,my_x_coord, iTmp, jTmp)
   fs_disp = nx_disp*ny_disp

   ! For this domain, not for the whole!!!
   ifLonGlobal = fu_ifLonGlobal( fu_grid(disp_buf%p4d(ind_m)%past%p2d(1)%idPtr))

   nz_disp = fu_number_of_fields(disp_buf%p4d(ind_m)%past%field3d)
   boundLength(1:4) = (/nx_disp, nx_disp, ny_disp, ny_disp/)




   select case (diag_rules%wind_method)
            
      case (opentop, hardtop, topdown, hardtop_weighted)
         ifPoisson = .false.

         if (any(diag_rules%wind_method == (/hardtop, hardtop_weighted/) )) then
            call msg("Diagnozing cell fluxcorrection with poisson....")
            ifPoisson = .True.
            ifBottomUp = .False.
         else
            if (diag_rules%wind_method == opentop) then
               call msg("Diagnozing cell fluxcorrection with opentop....")
               ifBottomUp = .True.
            else
               call msg("Diagnozing cell fluxcorrection with topdown....")
               ifBottomUp = .false.
            endif

         endif

         disp_buf%ifMassFluxBottomUp = ifBottomUp ! Quick and dirty
         
         valid  =      fu_valid_time(disp_buf%p4d(ind_m)%past%p2d(1)%idptr)
         analysis = fu_analysis_time(disp_buf%p4d(ind_m)%future%p2d(1)%idptr)
         validity =   fu_valid_time(disp_buf%p4d(ind_m)%future%p2d(1)%idptr) &
                  & - fu_valid_time(disp_buf%p4d(ind_m)%past%p2d(1)%idptr) 
         seconds = fu_sec(validity)
         
         if (ifPoisson)  then 
            ! Need to ensure that domains have exactly the same idea on 
            ! boundary masses. Same stuff is used for boundary masses in poisson.
            boundLength(1:4) = (/nx_disp, nx_disp, ny_disp, ny_disp/)
            fBoundWAptr => fu_work_array((2*nx_disp+ 2*ny_disp)*2*nz_disp)
            do iTmp = 1,4
               jTmp = sum(boundLength(1:iTmp-1))*2*nz_disp
               MBflat(iTmp)%ptr(1:boundLength(iTmp)*2*nz_disp) => &
                    & fBoundWAptr(jTmp+1:jTmp+boundLength(iTmp)*2*nz_disp)
            enddo
            MBN(1:nx_disp, 1:nz_disp, 1:2) =>  MBflat(northern_boundary)%ptr(:)
            MBS(1:nx_disp, 1:nz_disp, 1:2) =>  MBflat(southern_boundary)%ptr(:)
            MBE(1:ny_disp, 1:nz_disp, 1:2) =>  MBflat(eastern_boundary)%ptr(:)
            MBW(1:ny_disp, 1:nz_disp, 1:2) =>  MBflat(western_boundary)%ptr(:)

            ! Fill stuff from behind-boundaries for transmission Only mass is needed
            do iLev = 1,nz_disp
               mbptr2d(1:nx_disp,1:ny_disp)     => disp_buf%p4d(ind_m)%past%p2d(iLev)%ptr
               maptr2d(1:nx_disp,1:ny_disp)     => disp_buf%p4d(ind_m)%future%p2d(iLev)%ptr
               MBN(:,iLev,our) = 0.5*(mbPtr2D(:,ny_disp) + maPtr2D(:,ny_disp))
               MBS(:,iLev,our) = 0.5*(mbPtr2D(:,1)       + maPtr2D(:,1))
               if (ifLonGlobal) then 
                  !average with opposite side
                  MBE(:,iLev,our) = 0.25 * (mbPtr2D(1,:) + maPtr2D(1,:) + &
                                      & mbPtr2D(nx_disp,:) + maPtr2D(nx_disp,:))
                  MBW(:,iLev,our) = MBE(:,iLev,our) 
               else !from this side
                  MBW(:,iLev,our) =  0.5*(mbPtr2D(1,:) + maPtr2D(1,:) )
                  MBE(:,iLev,our) =  0.5*(mbPtr2D(nx_disp,:) + maPtr2D(nx_disp,:))
               endif
            enddo

            ! Exchange boundaries and average with neighbour 
            do iTmp=1,4 !Boundaries
               if (mod(my_y_coord+my_x_coord,2) == 0) then !Should be different order for the neighbour
                   !2,1,4,3
                   iBound = iTmp + 2*mod(iTmp,2) - 1 ! iTmp-1 for even, iTmp+1 for odd iTmp
               else
                   iBound = iTmp ! 1,2,3,4
               endif
               if(adv_mpi_neighbours(iBound)>=0) then
                  !call msg(fu_str(smpi_adv_rank)+"Exchanging MassBehind with neighbor, boundary:"+fu_str(iBound))
                  !call flush()
                   call msg(fu_str(smpi_adv_rank)+"Exchanging FMB with neighbor, boundary:"+fu_str(iBound))
                   call msg("boundLength(iTmp)", boundLength(iTmp), nz_disp)
                   sendcnt = boundLength(iTmp)*nz_disp
                  fPtr  => MBflat(iBound)%ptr(1:sendcnt)  !!! Our
                  fPtr1 => MBflat(iBound)%ptr(sendcnt+1:2*sendcnt) !!!Their
                  flush(run_log_funit)

                   call smpi_exchange_wings(adv_mpi_neighbours(iBound), &
                           & fPtr, sendcnt,  fPtr1, rcvcnt)
                   if (sendcnt /= rcvcnt) then
                     call msg("after smpi_exchange_wings mass sendcnt,recvcnt", sendcnt,rcvcnt )
                     call set_error("sendcnt and recvcnt do not match", sub_name)
                     return
                   endif
                  fPtr =  0.5*(fPtr + fPtr1)  !! Take a compromise as valid for this boundary
               endif
            enddo
            ! MBX(:,:,our) now has the right mass for boundary fluxes adjustment


            fWAptr => fu_work_array(2*fs_disp) 
            colmass(1:fs_disp) => fWAptr(1:fs_disp)
            colmassdelta(1:fs_disp) => fWAptr(fs_disp+1:2*fs_disp)

           
            !integtate masses over vertical
            colmass(1:fs_disp) = 0.
            colmassdelta(1:fs_disp) = 0.
            do iLev = 1,fu_NbrOfLevels(dispersion_vertical)
              colmass(1:fs_disp) = colmass(1:fs_disp) + & 
                   & Disp_Wind%LevelCorrWeight(iLev) * &
                   & 0.5 * (disp_buf%p4d(ind_m)%past%p2d(iLev)%ptr(:) + disp_buf%p4d(ind_m)%future%p2d(iLev)%ptr(:))
              colmassdelta(1:fs_disp) = colmassdelta(1:fs_disp) + &
                   & ((disp_buf%p4d(ind_m)%future%p2d(iLev)%ptr(:) - disp_buf%p4d(ind_m)%past%p2d(iLev)%ptr(:)))
            enddo
            colmassdelta(1:fs_disp) = colmassdelta(1:fs_disp)/seconds ! To kg/s

            ! Place to get the correction
            uPtr => fu_work_array((gnx+1)*gny) 
            vPtr => fu_work_array(gnx*(gny+1)) 
            uPtr2D(1:gnx+1,1:gny) => uPtr(1:(gnx+1)*gny)
            vPtr2D(1:gnx,1:gny+1) => vPtr(1:gnx*(gny+1))



            if (gnx == nx_disp .and.  gny == ny_disp) then ! All in one domain
               mPtr2D(1:gnx,1:gny)  =>  colmass(1:fs_disp) 
               dmPtr2D(1:gnx,1:gny)  =>  colmassdelta(1:fs_disp) 
            else

               mLptr => fu_work_array(gnx*gny)
               mPtr2D(1:gnx,1:gny) => mLptr(1:gnx*gny)
              call smpi_gather_field(colmass, mPtr2D , .false., .false.)

               dmLptr => fu_work_array(gnx*gny)
               dmPtr2D(1:gnx,1:gny) => dmLptr(1:gnx*gny)
               call smpi_gather_field(colmassdelta, dmPtr2D , .false., .false.)

            endif

            masscorrfactor = 0 !! Just make it initialized for non-master
            if (smpi_adv_rank == 0) then ! Master
              call adjust_2D_fluxes(uPtr2D, vPtr2D, dmPtr2D, mPtr2D, &
                         & fu_ifLonGlobal(wholeMPIdispersion_grid), &
                         & fu_ifPolarCapPossible(wholeMPIdispersion_grid, southern_boundary), &
                         & fu_ifPolarCapPossible(wholeMPIdispersion_grid, northern_boundary), &
                         & gnx, gny, masscorrfactor)  
            endif
            
            if (gnx /= nx_disp .or.  gny /= ny_disp) then
              ! Make uPtr2D and vPtr2D to look to local arrays,
              ! Get masscorrfactor
               fPtr => fu_work_array((gnx+1)*(gny+1))
               fPtr(1:(gnx+1)*gny) = uPtr(1:(gnx+1)*gny)
               !Fit masscorrfactor to the same transaction
               fPtr((gnx+1)*gny+1) = real( masscorrfactor)
               call smpi_bcast_aray(fptr, (gnx+1)*gny+1, smpi_adv_comm,0)
               masscorrfactor =  fPtr((gnx+1)*gny+1)

               tmpptr2d(1:gnx+1,1:gny) => fptr(1:(gnx+1)*gny)
               uPtr2D(1:nx_disp+1,1:ny_disp)   => uPtr(1:fs_disp+ny_disp)
               uPtr2D(1:nx_disp+1,1:ny_disp) = tmpptr2d(offx+1:offx+nx_disp+1,offy+1:offy+ny_disp)

               fPtr(1:(gny+1)*gnx) = vPtr(1:(gny+1)*gnx)
               call smpi_bcast_aray(fPtr, gnx*(gny+1),smpi_adv_comm,0)
               tmpptr2d(1:gnx,1:gny+1) => fPtr(1:(gny+1)*gnx)
               vPtr2D(1:nx_disp,1:ny_disp+1) => vPtr(1:fs_disp+nx_disp)
               vPtr2D(1:nx_disp,1:ny_disp+1) = tmpptr2d(offx+1:offx+nx_disp,offy+1:offy+ny_disp+1)

               ! These did not change 
               mPtr2D(1:nx_disp,1:ny_disp)  => colmass(1:fs_disp)
               dmPtr2D(1:nx_disp,1:ny_disp)  => colmassdelta(1:fs_disp)

               !Cleanup
               call free_work_array(fPtr)
               call free_work_array(mLptr)
               call free_work_array(dmLptr)
            endif 

         endif ! Poisson


          if (ifBottomUp) then
            ! Apply horizontal wind  corrections and diagnoze vertical wind bottom-up
            if (ifPoisson) then
                 call set_error("No Poisson with bottom-Up diagnostics",sub_name)
                 return
            endif
            !Pretend it was flux correction from the underground level
            tmpPtr2D(1:nx_disp,1:ny_disp) => disp_buf%p4d(ind_fwrt)%present%p2d(1)%ptr
            tmpPtr2D(1:nx_disp,1:ny_disp) = 0. !No flux from below
            do iLev = 1,fu_NbrOfLevels(dispersion_vertical)
               uLptr2D(1:nx_disp+1,1:ny_disp) => disp_buf%p4d(ind_furt)%present%p2d(iLev)%ptr 
               vLptr2D(1:nx_disp,1:ny_disp+1) => disp_buf%p4d(ind_fvrt)%present%p2d(iLev)%ptr 
               wLptr2D(1:nx_disp,1:ny_disp)   => disp_buf%p4d(ind_fwrt)%present%p2d(iLev)%ptr 

               mbptr2d(1:nx_disp,1:ny_disp)   => disp_buf%p4d(ind_m)%past%p2d(ilev)%ptr
               maptr2d(1:nx_disp,1:ny_disp)   => disp_buf%p4d(ind_m)%future%p2d(ilev)%ptr

               if(ifPoisson) then 
                 fTmp = Disp_Wind%LevelCorrWeight(iLev)
                 ! U
                 uLptr2D(2:nx_disp,1:ny_disp) =  - uPtr2D(2:nx_disp,1:ny_disp)* &
                     & 0.25 *fTmp* (mbptr2D(2:nx_disp,1:ny_disp)+ mbptr2D(1:nx_disp-1,1:ny_disp) + &
                     &         maptr2D(2:nx_disp,1:ny_disp)+ maptr2D(1:nx_disp-1,1:ny_disp))
                 uLptr2D(1,1:ny_disp) =   - uPtr2D(1,1:ny_disp) *fTmp* MBW(:,iLev,our) 
                 uLptr2D(nx_disp+1,1:ny_disp) = - uPtr2D(nx_disp+1,1:ny_disp) *fTmp* MBE(:,iLev,our)
                 ! V
                 vLptr2D(1:nx_disp,2:ny_disp) =  -vPtr2D(1:nx_disp,2:ny_disp) * &
                       & 0.25 *fTmp* (mbptr2D(1:nx_disp,2:ny_disp)+ mbptr2D(1:nx_disp,1:ny_disp-1) + &
                       &        maptr2D(1:nx_disp,2:ny_disp)+ maptr2D(1:nx_disp,1:ny_disp-1) )  
                 vLptr2D(1:nx_disp,1) =   - vPtr2D(1:nx_disp,1) *fTmp* MBS(:,iLev,our)
                 vLptr2D(1:nx_disp,ny_disp+1) =  - vPtr2D(1:nx_disp,ny_disp+1) *fTmp* MBN(:,iLev,our)
                 ! W
                 wLPtr2D(1:nx_disp,1:ny_disp) = tmpPtr2D(1:nx_disp,1:ny_disp)  & !Flux from below
                     & +  vLptr2D(1:nx_disp,1:ny_disp)  -  vLptr2D(1:nx_disp,2:ny_disp+1) & !Y-convergence
                     & +  uLptr2D(1:nx_disp,1:ny_disp)  -  uLptr2D(2:nx_disp+1,1:ny_disp) & !X-convergence
                     & +  (mbptr2d(1:nx_disp,1:ny_disp) -  maptr2d(1:nx_disp,1:ny_disp))/seconds   ! Masschange
                 !  

                 ! If masschange correction used (global, hardtop)
                 if (masscorrfactor /= 0.) then
                    wLPtr2D(1:nx_disp,1:ny_disp)  = wLPtr2D(1:nx_disp,1:ny_disp) &
                       & - 0.5*fTmp*(mbptr2d(1:nx_disp,1:ny_disp) + maptr2d(1:nx_disp,1:ny_disp))*masscorrfactor
                 endif
                 !Point to bottom of the next level == top of this one
                 tmpPtr2D(1:nx_disp,1:ny_disp) => disp_buf%p4d(ind_fwrt)%present%p2d(iLev)%ptr

               else ! No corrections for U and V
               ! correction to w just due to mass change
               uLptr2D(1:nx_disp+1,1:ny_disp) = 0.
               vLptr2D(1:nx_disp,1:ny_disp+1) = 0.
               wLptr2D(1:nx_disp,1:ny_disp) =  tmpPtr2D(1:nx_disp,1:ny_disp) + &
                  & (mbptr2d(1:nx_disp,1:ny_disp) -  maptr2d(1:nx_disp,1:ny_disp))/seconds
               tmpPtr2D(1:nx_disp,1:ny_disp) => disp_buf%p4d(ind_fwrt)%present%p2d(iLev)%ptr

   !            wLptr(1:fs_disp) = p(1:fs_disp) +  (mbptr(1:fs_disp) - maptr(1:fs_disp)) / seconds
   !            p => wLptr !switch to current level
               endif ! poisson
                
             enddo !Levels
             call msg("Sum of flux at domain top", sum(tmpPtr2D)) 
             call msg("Sum of upward flux at domain top", sum(tmpPtr2D, mask=(tmpPtr2D.gt.0))) 
           else 
              !Top-Down
              ! Apply horizontal wind  corrections and diagnoze vertical wind correction Top-Down
              tmpPtr2D(1:nx_disp,1:ny_disp) => disp_buf%p4d(ind_fwrt)%present%p2d(fu_NbrOfLevels(dispersion_vertical))%ptr
              tmpPtr2D(1:nx_disp,1:ny_disp) = 0. !No flux from above
              do iLev = fu_NbrOfLevels(dispersion_vertical),1,-1
                  uLptr2D(1:nx_disp+1,1:ny_disp) => disp_buf%p4d(ind_furt)%present%p2d(iLev)%ptr 
                  vLptr2D(1:nx_disp,1:ny_disp+1) => disp_buf%p4d(ind_fvrt)%present%p2d(iLev)%ptr
                  if (iLev /= 1) then 
                     wLptr2D(1:nx_disp,1:ny_disp) =>  disp_buf%p4d(ind_fwrt)%present%p2d(iLev-1)%ptr !Below this layer
                  else
                     wLptr2D => null()  !Should not be used
                  endif

                  mbptr2d(1:nx_disp,1:ny_disp)     => disp_buf%p4d(ind_m)%past%p2d(ilev)%ptr
                  maptr2d(1:nx_disp,1:ny_disp)     => disp_buf%p4d(ind_m)%future%p2d(ilev)%ptr

                  if(ifPoisson) then 
                    fTmp = Disp_Wind%LevelCorrWeight(iLev)
                    ! U
                    uLptr2D(2:nx_disp,1:ny_disp) =  - uPtr2D(2:nx_disp,1:ny_disp)* &
                        & 0.25 *fTmp* (mbptr2D(2:nx_disp,1:ny_disp)+ mbptr2D(1:nx_disp-1,1:ny_disp) + &
                        &         maptr2D(2:nx_disp,1:ny_disp)+ maptr2D(1:nx_disp-1,1:ny_disp))
                    uLptr2D(1,1:ny_disp) =   - uPtr2D(1,1:ny_disp) *fTmp* MBW(:,iLev,our) 
                    uLptr2D(nx_disp+1,1:ny_disp) = - uPtr2D(nx_disp+1,1:ny_disp) *fTmp* MBE(:,iLev,our)
                    ! V
                    vLptr2D(1:nx_disp,2:ny_disp) =  -vPtr2D(1:nx_disp,2:ny_disp) * &
                          & 0.25 *fTmp* (mbptr2D(1:nx_disp,2:ny_disp)+ mbptr2D(1:nx_disp,1:ny_disp-1) + &
                          &        maptr2D(1:nx_disp,2:ny_disp)+ maptr2D(1:nx_disp,1:ny_disp-1) )  
                    vLptr2D(1:nx_disp,1) =   - vPtr2D(1:nx_disp,1) *fTmp* MBS(:,iLev,our)
                    vLptr2D(1:nx_disp,ny_disp+1) =  - vPtr2D(1:nx_disp,ny_disp+1) *fTmp* MBN(:,iLev,our)

                     if (iLev /= 1) then !Need to adjust below-layer W
                       ! W
                       wLPtr2D(1:nx_disp,1:ny_disp) = tmpPtr2D(1:nx_disp,1:ny_disp)  & !Flux from above
                           & -  vLptr2D(1:nx_disp,1:ny_disp)  +  vLptr2D(1:nx_disp,2:ny_disp+1) & !Y-convergence
                           & -  uLptr2D(1:nx_disp,1:ny_disp)  +  uLptr2D(2:nx_disp+1,1:ny_disp) & !X-convergence
                           & -  (mbptr2d(1:nx_disp,1:ny_disp) -  maptr2d(1:nx_disp,1:ny_disp))/seconds   ! Masschange
                       !  
                       ! If masschange correction used (global, hardtop)
                       if (masscorrfactor /= 0.) then
                          wLPtr2D(1:nx_disp,1:ny_disp)  = wLPtr2D(1:nx_disp,1:ny_disp) &
                             & + 0.5*fTmp*(mbptr2d(1:nx_disp,1:ny_disp) + maptr2d(1:nx_disp,1:ny_disp))*masscorrfactor
                       endif
                       tmpPtr2D(1:nx_disp,1:ny_disp) => disp_buf%p4d(ind_fwrt)%present%p2d(iLev-1)%ptr
                     endif

                  else !No Poisson -- no corrections for U and V

                     ! correction due to mass change to W only, u and v untouched 
                     uLptr2D(1:nx_disp+1,1:ny_disp) = 0.
                     vLptr2D(1:nx_disp,1:ny_disp+1) = 0.
                     if (iLev /= 1) then
                        wLptr2D(1:nx_disp,1:ny_disp) =  tmpPtr2D(1:nx_disp,1:ny_disp)  &
                            & - (mbptr2d(1:nx_disp,1:ny_disp) -  maptr2d(1:nx_disp,1:ny_disp))/seconds
                        !Point to next-layers top == this layers bottom
                        tmpPtr2D(1:nx_disp,1:ny_disp) => disp_buf%p4d(ind_fwrt)%present%p2d(iLev-1)%ptr
                     else

                     endif
                  endif
                
              enddo !Levels

              call msg("Sum of remaining flux at 1st (or last) layer top", sum(tmpPtr2D)) 
              call msg("Sum of remaining upward flux at 1st layer top", sum(tmpPtr2D, mask=(tmpPtr2D.gt.0)))
              call msg("Through-ground flux is assumed zero anyhow...")
             ! Here we have 2D divergence as a vertical velocity at the top of the upper layer
             ! Should be close to zero in case of hardtop


           endif  !Top-down/bottom-up   

           if (ifPoisson) then
             call free_work_array(uPtr)
             call free_work_array(vPtr)
             call free_work_array(fWAptr)
             call free_work_array(fBoundWAptr)
           endif

         
       case(test_wind) 
          !Valid forever by default
          valid = really_far_in_past
          analysis = valid
          validity = very_long_interval

          !Should use whole MPI grid and offsets to avoid inconsistences in MPI runs
          call  lonlat_grid_parameters(wholempidispersion_grid,&
                           & corner_lon_E, corner_lat_N,corner_in_geo_lonlat, &
                           & number_of_points_x, number_of_points_y, &
                           & southpole_lon_E, southpole_lat_N, &
                           & dx_deg, dy_deg)
            ! Only our size and offsets are needed

             if (error) return
          select case (diag_rules%iTestWindType)
            case(0) ! 0: zero wind
               do iLev = 1,fu_NbrOfLevels(dispersion_vertical)
                  disp_buf%p4d(ind_furt)%present%p2d(iLev)%ptr(:) = 0.
                  disp_buf%p4d(ind_fvrt)%present%p2d(iLev)%ptr(:) = 0.
                  disp_buf%p4d(ind_fwrt)%present%p2d(iLev)%ptr(:) = 0.
               enddo
               call msg("Zero wind, valid forever")
               valid = really_far_in_past
               analysis = valid
               validity = very_long_interval
            case(2) !Varying Courant, divergent
                p => disp_buf%p4d(ind_furt)%present%p2d(1)%ptr(:)
                maptr => disp_buf%p4d(ind_m)%future%p2d(1)%ptr

                do iy = 1, ny_disp
                  do ix = 1, nx_disp
                    iTmp = ix + (iy-1) * nx_disp !non_staggered
                    i1d = ix + (iy-1) * (nx_disp+1)!Staggered index
                    p(i1d) = maptr(iTmp) * 1e-3 *(1.1 + sin(real(ix+offx)/15.))
                  end do
                  i1d = ix + (iy-1) * (nx_disp+1)
                  p(i1d) = maptr(iTmp) * 1e-3 *(1.1 + sin(real(ix+offx)/15.))
                end do

                do iLev = 1,fu_NbrOfLevels(dispersion_vertical)
                     disp_buf%p4d(ind_furt)%present%p2d(iLev)%ptr(:) =  p(:)
                     disp_buf%p4d(ind_fvrt)%present%p2d(iLev)%ptr(:) = 0.
                     disp_buf%p4d(ind_fwrt)%present%p2d(iLev)%ptr(:) = 0.
                enddo
                call msg('Linear varying-Courant test (.1~2.1*1e-3, 0, 0) cells/sec:')


             case(6)   ! 6: Pole-to-pole: to north at bottom, to south at top
                do iLev = 1,fu_NbrOfLevels(dispersion_vertical)
                     disp_buf%p4d(ind_furt)%present%p2d(iLev)%ptr(:) = 0.
                     disp_buf%p4d(ind_fvrt)%present%p2d(iLev)%ptr(:) = 0.
                     disp_buf%p4d(ind_fwrt)%present%p2d(iLev)%ptr(:) = 0.
                enddo
                p => disp_buf%p4d(ind_fvrt)%present%p2d(1)%ptr(:)
                fTmp = pi*earth_radius/T*earth_radius*dx_deg*degrees_to_radians !Half-turn per T
                p(1:fs_disp+nx_disp) = fTmp
                iTmp = fu_NbrOfLevels(dispersion_vertical)
                do iLev = 1,iTmp/2-1
                     disp_buf%p4d(ind_fvrt)%present%p2d(iLev)%ptr(:) =  p(:)
                     disp_buf%p4d(ind_fvrt)%present%p2d(iTmp-iLev+1)%ptr(:) =  -p(:)
                enddo
                call msg('Pole-to-pole circulation')


             !Stream-function based diagnostics 2D non-divergent winds
             case(7,8,9,10)     
               
               select case (diag_rules%iTestWindType)

                case (7) !Lauritzen
                  pi_t_over_T =( pi * modulo(nint(fu_sec(now  - ref_time_01012000)), 2*T)) / T 
                  call msg("Lauritzen phase (seconds, of total):", (modulo(nint(fu_sec(now - ref_time_01012000)), 2*T)  ), T)
                  valid = now
                  analysis = now
                  validity = one_second

                  tfactor = earth_radius/T*cos(pi_t_over_T)
                  do iy = 1, ny_disp+1
                     fTmp = (corner_lat_N + (iY+offy-1.5)* dy_deg)* degrees_to_radians
                     cosTheta2 = cos(fTmp)
                     cosTheta2 = cosTheta2*cosTheta2
                     sinTheta = sin (ftmp)
                     do ix = 1, nx_disp+1
                       lon = corner_lon_E + dx_deg *(ix+offx-1.5)
                       lambda_prime = lon * degrees_to_radians - 2 * pi_t_over_T 
                       sinlprime2 = sin(lambda_prime) ! sin^2
                       sinlprime2 = sinlprime2 * sinlprime2
                       Disp_wind%StreamFunction(ix,iy) = (-10.* tfactor * sinlprime2*cosTheta2 + 2*pi*earth_radius/T*sinTheta)*earth_radius
                     end do
                  enddo
                case(8) !Solid body rotation pole-to-polea
                  call msg('Solid body rotation Pole-to-pole')
                  do iy = 1, ny_disp+1
                     lat = corner_lat_N + dy_deg *((iy+offy)-1.5)
                     fTmp = 2*pi*earth_radius/T*earth_radius*cos(lat * degrees_to_radians) 
                     do ix = 1, nx_disp+1
                        lon = corner_lon_E + dx_deg *((ix+offx)-1.5)
                           Disp_wind%StreamFunction(ix,iy) = fTmp * sin(lon * degrees_to_radians)
                     enddo
                  enddo
                case(9) !Solid body rotation Along equator
                  call msg('Solid body rotation Along equator')
                  do iy = 1, ny_disp+1
                     fTmp = -(corner_lat_N + (iY+offy-1.5)* dy_deg)* degrees_to_radians
                     Disp_wind%StreamFunction(:,iy) = 2*pi*earth_radius*earth_radius/T*sin(fTmp) !Full turn / period
                  enddo

                case(10) !Solid body rotation 45 deg to equator
                  call msg('Solid body rotation around 45S 90W')
                  pole45 = fu_set_pole(south_flag, -45., -90.)
                  do iy = 1, ny_disp+1
                     lat = corner_lat_N + dy_deg *((iy+offy)-1.5)
                     do ix = 1, nx_disp+1
                        lon = corner_lon_E + dx_deg *((ix+offx)-1.5)
                        call modify_lonlat(lat, lon, pole_geographical, pole45, rlat, rlon)
                        Disp_wind%StreamFunction(ix,iy) = 0.1 * 2*pi*earth_radius*earth_radius/T*sin(rlat* degrees_to_radians)
                     enddo
                  enddo

                end select
                !Make streamfunction exactly matching at boundaries
                if (fu_ifLonGlobal(dispersion_grid)) then
                   Disp_wind%StreamFunction(1,:)=0.5*(Disp_wind%StreamFunction(1,:) + Disp_wind%StreamFunction(nx_disp+1,:))
                   Disp_wind%StreamFunction(nx_disp+1,:) = Disp_wind%StreamFunction(1,:)
                endif
                  
                !Make first level
                uLptr => disp_buf%p4d(ind_furt)%present%p2d(1)%ptr(:)
                vLptr => disp_buf%p4d(ind_fvrt)%present%p2d(1)%ptr(:)
                disp_buf%p4d(ind_fwrt)%present%p2d(1)%ptr(:) = 0.

                do iy = 1, ny_disp
                    i1d = 1+(nx_disp+1)*(iy-1) !x-staggered
                    uLptr(i1d:i1d+nx_disp) = Disp_wind%StreamFunction(:,iy+1) - &
                                            & Disp_wind%StreamFunction(:,iy)
                enddo
                do iy = 1, ny_disp+1
                    i1d = 1+nx_disp*(iy-1) !y-staggered
                    vLptr(i1d:i1d+nx_disp-1) = Disp_wind%StreamFunction(1:nx_disp,iy) - &
                                            & Disp_wind%StreamFunction(2:nx_disp+1,iy)
                enddo
                do iLev = 2,fu_NbrOfLevels(dispersion_vertical)
                     disp_buf%p4d(ind_furt)%present%p2d(iLev)%ptr(:) = uLptr(:)
                     disp_buf%p4d(ind_fvrt)%present%p2d(iLev)%ptr(:) = vLptr(:)
                     disp_buf%p4d(ind_fwrt)%present%p2d(iLev)%ptr(:) = 0.
                enddo


            case default
                call msg ("rules%iTestWindType =", diag_rules%iTestWindType)
                call set_error("Not implemented (yet?):",sub_name)
            end select ! Test wind method


       case default
          call msg ("diag_rules%wind_method=",diag_rules%wind_method)
          call set_error("Unknown massflux diagnostic method",sub_name)
    end select

    !
    ! Almost done. Now set proper validity in field IDs
    ! Note, that only level 1 is checked above....
    !

     do iLev = 1,fu_NbrOfLevels(dispersion_vertical)
       call set_valid_time     (disp_buf%p4d(ind_furt)%present%p2d(iLev)%idptr, valid)
       call set_analysis_time  (disp_buf%p4d(ind_furt)%present%p2d(iLev)%idptr, analysis)
       call set_validity_length(disp_buf%p4d(ind_furt)%present%p2d(iLev)%idptr, validity)
       call set_valid_time(     disp_buf%p4d(ind_fvrt)%present%p2d(iLev)%idptr, valid)
       call set_analysis_time(  disp_buf%p4d(ind_fvrt)%present%p2d(iLev)%idptr, analysis)
       call set_validity_length(disp_buf%p4d(ind_fvrt)%present%p2d(iLev)%idptr, validity)
       call set_valid_time(     disp_buf%p4d(ind_fwrt)%present%p2d(iLev)%idptr, valid)
       call set_analysis_time(  disp_buf%p4d(ind_fwrt)%present%p2d(iLev)%idptr, analysis)
       call set_validity_length(disp_buf%p4d(ind_fwrt)%present%p2d(iLev)%idptr, validity)
     enddo
     ! fieldID3d times should be kept coherently -- otherwise arrange stack gets
     ! mad...
     call set_valid_time(disp_buf%p4d(ind_furt)%present%field3d, valid)
     call set_valid_time(disp_buf%p4d(ind_fvrt)%present%field3d, valid)
     call set_valid_time(disp_buf%p4d(ind_fwrt)%present%field3d, valid)
    !
    ! Adjust the borders of the global validity period of all single-time fields
    !
    if(valid_time_border_1 < valid) valid_time_border_1 = valid

    analysis = valid + validity ! Recycle fariable
    if(valid_time_border_2 > analysis)  valid_time_border_2 =  analysis

  end subroutine  df_update_cellfluxcorr_rt


  !***********************************************************************************************

  subroutine df_daily_mean_variable(quantity_accum, quantity_mean, iAvType, buf, now, &
                                  & valid_time_border_1, valid_time_border_2, ifColdstartAllowed)
    ! 
    ! Mean temperature over _previous_ UTC day (00Z-00Z)
    ! Field set as average from 00Z-00Z for today
    !
    implicit none

    ! Imported parameters
    integer, intent(in) :: quantity_accum, quantity_mean, iAvType
    type(Tfield_buffer), INTENT(in) :: buf   ! actually, dispersion buffer
    type(silja_time) :: now                  ! Used only for temporary solution
    type(silja_time), intent(inout) :: valid_time_border_1, valid_time_border_2
    type(silja_logical), intent(in) :: ifColdstartAllowed

    ! Local variables
    TYPE(silja_field_id), pointer :: pIdCum_past, pIdCum_future, pIdMean
    REAL, DIMENSION(:), POINTER :: cum_T_past, cum_T_future, mean_T
    type(silja_time) :: time_past, time_future, start_of_day
    type(silja_time) :: time_tmp
    integer :: iT_mean, iT_cum, iLev, nLevs, nx, ny, fs
    real :: seconds
    logical :: ifMakeIt, ifMakeSurrogate, stillValid
    character(len=22), parameter :: sub_name = 'df_daily_mean_variable'
    
    integer, save :: callCnt = 0

    callCnt = callCnt + 1
    !
    ! Do we need to do anything?
    ! is mean field valid for the time range between pst and future?
    !
    ! Preparatory steps
    !
    iT_cum  = fu_index(buf%buffer_quantities, quantity_accum)
    iT_mean = fu_index(buf%buffer_quantities, quantity_mean)
    if (any((/iT_cum, iT_mean/) < 1 )) then
      call msg("(/iT_cum, iT_mean/)",(/iT_cum, iT_mean/))
      call set_error("Could not find indices in a buffer", sub_name)
      return
    endif

    time_past = buf%time_past
    time_future = buf%time_future

    pIdCum_past => buf%p2d(iT_cum)%past%IDptr
    pIdCum_future => buf%p2d(iT_cum)%future%IDptr
    pIdMean => buf%p2d(iT_mean)%present%IDptr

    cum_T_past => buf%p2d(iT_cum)%past%ptr
    cum_T_future => buf%p2d(iT_cum)%future%ptr
    mean_T => buf%p2d(iT_mean)%present%ptr   ! single-time field, only present

    stillValid = (fu_hour(time_past) /= 0 .and. callCnt > 2) 
      !Always rediagnoze for after-midnight interval
      ! or on the second call
    if (stillValid) stillValid =  (fu_valid_time(pIdMean) >= time_future) 
    if (stillValid) stillValid =  &
        & (fu_valid_time(pIdMean) - fu_accumulation_length(pIdMean) <= time_future) 

    if( stillValid) then
      call msg("Still valid: "//fu_quantity_string(quantity_mean))
      ! av_end (end_of "applicable time")
      time_tmp = fu_valid_time(pIdMean)
      if(valid_time_border_2 > time_tmp) valid_time_border_2 = time_tmp
      ! av_start (start of "applicable time")
      time_tmp = fu_valid_time(pIdMean) - fu_accumulation_length(pIdMean)
      if(valid_time_border_1 < time_tmp) valid_time_border_1 = time_tmp
      return
    endif

    fs = fu_number_of_gridpoints(fu_grid(pIdMean))

    ! Few options:
    ! 1. Very first diagnostics (pIdMean valid really_far_in_past):  
    !     make an instant field valid for "now" with analysis  really_far_in_past
    !     to satisfy the first output (Make Instant)
    ! 2. Second diagnostic, not initialised (instant with analysis really_far_in_past):
    !     If possible, Rediagnose from cum_T_past for 24h
    !       otherwise make a surrogate if allowed, or crash (Make )
    ! 3. Second diagnostic, initialised (instant field with analysis not really_far_in_past ):
    !     If it is time, Rediagnose from cum_T_past for 24h
    !     otherwise just make it 24-h valid.
    ! 4. Regular diagnostics: 24-hour average field:
    !      Just make it     

    !! FIXME
    ! For now, just make it if we can or reset the averaging time

    ifMakeIt = fu_accumulation_length(pIdCum_past) == one_day .and. fu_hour(time_past) == 0
    if (ifMakeIt) then 
      !
      ! Just reset it: next day came
      !
      call msg("Recalculating daymean T")
      call msg("Before mean_T(1:10)", mean_T(1:10))
      select case(iAvType)
        case(av4mean)
          mean_T(1:fs) = cum_T_past(1:fs) / seconds_in_day
        case(av4sum, av4min, av4max)
          mean_T(1:fs) = cum_T_past(1:fs)
        case default
          call set_error('Unknown averaging type 1:' + fu_str(iAvType), 'df_daily_mean_variable')
          return
      end select
      call msg("After mean_T(1:10)", mean_T(1:10))

    elseif (mean_T(1) == real_missing) then
      !
      ! fingerprint of gobal_io_init is recognised
      !
      if (CallCnt == 1) then
        !
        ! Field was just created by gobal_io_init and not reset. Allowed?
        !
        if (fu_true(ifColdstartAllowed)) then
          !
          ! Put coldstart values here, so next time will not crash on missing ini
          !
          call msg("Filling daymean T with something meaningful for the first out")
          select case(iAvType)
            case(av4mean)
              mean_T(1:fs) = cum_T_past(1:fs) / fu_sec(fu_accumulation_length(pIdCum_past))
            case(av4sum, av4min, av4max)
              mean_T(1:fs) = cum_T_past(1:fs)
            case default
              call set_error('Unknown averaging type 2:' + fu_str(iAvType), 'df_daily_mean_variable')
              return
          end select

        elseif (fu_false(ifColdstartAllowed)) then
          ! Forbidden thing happened
          call set_error('Cold start is not allowed but happened','df_daily_mean_variable')
          return
        else
          ! Do not know what to do
!          call msg("Leaving real_missing daymean T ")
          call set_error('Cold start permission undefined but it happened','df_daily_mean_variable')
          return
        endif
      else
        call set_error("Uninitialized daymean T!",  sub_name)
        return
      endif
    else
      call msg("Keeping mean_T(1:10)", mean_T(1:10))
      call msg_warning("Assuming daymean T was initialized or recalcualated") 
    endif
    call msg("resetting ID for daymean T")
    !(Re)set the field ID, that might have got outdated and/or came broken from initialization
     !! It is just an average field for the whole day
    start_of_day = fu_start_of_day_utc(time_past)
    pIdMean =  fu_set_field_id(met_src_missing,&
                             & quantity_mean, &
                             & start_of_day,& ! Analysis
                             & one_day, &             ! Forecast length 
                             & fu_grid(pIdCum_past),&
                             & fu_level(pIdCum_past),&
                             & one_day, & ! optional accumulation interval
                             & zero_interval, & ! !! We should get rid of validity_length at some point
                             & averaged_flag) !field_kind

    if (callCnt < 2) then !! Force recalculation on the second call
      time_tmp = start_of_day
    else
      time_tmp = start_of_day + one_day
    endif
      if(valid_time_border_1 < start_of_day) valid_time_border_1 = start_of_day
      if(valid_time_border_2 > time_tmp) valid_time_border_2 = time_tmp

  end subroutine df_daily_mean_variable

  
  !**********************************************************************************************
  
  integer function fu_accumulation_type(quantity)
    !
    ! Averaged quantities can be period-mean, period-max and period-min. This function
    ! checks the action judjing from the quantity
    !
    implicit none
    
    ! Imported parameter
    integer, intent(in) :: quantity
    
    select case(quantity)
      case(day_temperature_acc_flag, day_mean_temperature_flag, &
         & day_temperature_2m_acc_flag, day_mean_temperature_2m_flag, &
         & day_windspeed_10m_acc_flag, day_mean_windspeed_10m_flag, &
         & day_relat_humid_2m_acc_flag, day_mean_relat_humid_2m_flag)

        fu_accumulation_type = av4Mean

      case(day_temperature_2m_acc_max_flag, day_max_temperature_2m_flag, &
         & day_windspeed_10m_acc_max_flag, day_max_windspeed_10m_flag)

        fu_accumulation_type = av4Max

      case(large_scale_accum_rain_flag, convective_accum_rain_flag, total_precipitation_acc_flag, &
         & day_precipitation_acc_flag, day_sum_precipitation_flag)

        fu_accumulation_type = av4sum

      case(day_relat_humid_2m_acc_min_flag, day_min_relat_humid_2m_flag)

        fu_accumulation_type = av4Min

      case default
        call set_error('Not an accumulation quantity:' + fu_quantity_string(quantity),'fu_accumulation_type')
        fu_accumulation_type = int_missing
    end select
    
  end function fu_accumulation_type
  
  
  !**********************************************************************************************
  
 subroutine adjust_2D_fluxes(U2D, V2D, deltaM2Drate, m2D, ifLonGlobal, ifSpole, ifNpole, nx, ny,&
                          &  masscorrfactor)
    !
    !Calls Poisson solver and gets corrections for column-integrated fluxes.
    ! The corrections are normalized by cell masses: mean of two cells for 
    ! inner cellboundaries and one cell for outer boundaries

    implicit none
    real, dimension(:,:), pointer, intent(inout) :: U2D, V2D !Fluxes on input, Corrections per unit mass on output
    real, dimension(:,:), pointer, intent(in) :: deltaM2Drate, m2D !DeltaM (can be null), colamss
    logical, intent(in) :: ifLonGlobal, ifSpole, ifNpole 
    integer, intent(in) ::  nx, ny !Dimensions of Phi
    real(r8k), intent(out) :: masscorrfactor  ! masschange correction per unit mass 
                   ! zero for non-global domains or if deltaM2D == NULL


    real(fish_kind) :: masschangerate, domainmass, inflowfactor, bmass, fTmp, absinflow, inflow
    integer :: fs_phi, iTmp, jTmp, iY

      
    masscorrfactor = 0.

    ! Adjust massfluxes at boundaries and prepare divergence for Poisson solver

    if (associated(deltaM2Drate)) then  ! Using DeltaM to correct, ignore Fluxes
      masschangerate = sum(deltaM2Drate(1:nx,1:ny)) 
      if (ifLonGlobal .and. ifNpole .and. ifSpole) then
         !adjust  masschange to force zero over domain
         domainmass      = sum(m2D(1:nx,1:ny)) 
         call msg("Global domain.. Adjusting masschange rate by",  -masschangerate)
         masscorrfactor = -masschangerate/domainmass !Non-conservative mass in closed domain
         Disp_wind%Phi(:,:) = deltaM2Drate(1:nx,1:ny) + m2D(1:nx,1:ny)*masscorrfactor
                        ! Will have to tune mass change rate in each cell for W calculations
         V2D(1:nx,1) = 0. 
         V2D(1:nx,ny+1) = 0. 
      else ! Need to adjust the boarder flow -- mass-weighted
         Disp_wind%Phi(:,:) = deltaM2Drate(1:nx,1:ny)
         ! Normalize by boarder cell masses
         bmass = 0.
         if (ifSpole) then 
                V2D(1:nx,1) = 0.    ! No inflow adjustment for poles
         else
                bmass = bmass + sum(m2D(1:nx,1))
         endif
         if (ifNpole) then
            V2D(1:nx,ny+1) = 0.
         else
            bmass = bmass + sum(m2D(1:nx,ny))
         endif 
         if (.not. ifLonGlobal) bmass = bmass + sum(m2D(1,1:ny)) +  sum(m2D(nx,1:ny))

         ! Adjust boundary fluxes and divergence
         inflowfactor = -masschangerate/bmass

         call   msg("Sum of masschange inside the domain before adj  (kg/s)",sum(Disp_wind%Phi(:,:)))
         if (.not. ifSpole) then
            V2D(1:nx,1) = inflowfactor *  m2D(1:nx,1)
            Disp_wind%Phi(:,1) = Disp_wind%Phi(:,1) + V2D(1:nx,1)
            call   msg("Sum of masschange inside the domain after S adj  (kg/s)",sum(Disp_wind%Phi(:,:)))
          endif
          if ( .not. ifNpole)  then 
            V2D(1:nx,ny+1) = - inflowfactor *  m2D(1:nx,ny)
            Disp_wind%Phi(:,ny) = Disp_wind%Phi(:,ny) - V2D(1:nx,ny+1)
            call   msg("Sum of masschange inside the domain after N adj  (kg/s)",sum(Disp_wind%Phi(:,:)))
         endif
      if (.not. ifLonGlobal) then
            U2D(1,1:ny)    =  inflowfactor * m2D(1,1:ny)
            U2D(nx+1,1:ny) =  - inflowfactor * m2D(nx,1:ny)
            Disp_wind%Phi(1,:) = Disp_wind%Phi(1,:) + U2D(1,1:ny)
            call   msg("Sum of masschange inside the domain after W adj  (kg/s)",sum(Disp_wind%Phi(:,:)))
            Disp_wind%Phi(nx,:) = Disp_wind%Phi(nx,:) - U2D(nx+1,1:ny)
            call   msg("Sum of masschange inside the domain after E adj  (kg/s)",sum(Disp_wind%Phi(:,:)))
      endif

      endif
   else ! Ignore  DeltaM and use MassFluxes
        ! Divergence
        Disp_wind%Phi(:,:) =  - V2D(1:nx,1:ny) + V2D(1:nx,2:ny+1) &
                            & - U2D(1:nx,1:ny) + U2D(2:nx+1,1:ny)

      !      call msg("DeltaU sum", sum(- U2D(1:nx,1:ny) + U2D(2:nx+1,1:ny), 1))
!        call msg("latconvergence", sum(Disp_wind%Phi(:,:), 1))
     !   call msg("domainmass", sum(m2D))
     !   call msg("into NP", sum(V2D(1:nx,ny+1)) )
     !   call msg("into SP", sum(V2D(1:nx,1)) )
        call   msg("Sum of divergence inside the domain before adj      (kg/s)",sum(Disp_wind%Phi(:,:)))
        absinflow = 0.
        inflow = 0.
        if (ifSpole) then
         V2D(1:nx,1) = sum(V2D(1:nx,1)) / nx
         Disp_wind%Phi(:,1) = Disp_wind%Phi(:,1) + V2D(1:nx,1) 
         call   msg("Sum of divergence inside the domain after SP adj   (kg/s)",sum(Disp_wind%Phi(:,:)))
      else
          inflow =  inflow + sum(V2D(1:nx,1))
          absinflow = absinflow + sum(abs( V2D(1:nx,1) ))
      endif
      !  call msg("latconvergence", sum(Disp_wind%Phi(:,:), 1))

        if (ifNpole) then
         V2D(1:nx,ny+1) = sum(V2D(1:nx,ny+1)) / nx
         Disp_wind%Phi(:,ny) = Disp_wind%Phi(:,ny) - V2D(1:nx,ny+1) 
         call   msg("Sum of divergence inside the domain after NP adj   (kg/s)",sum(Disp_wind%Phi(:,:)))
      else
          inflow =  inflow - sum(V2D(1:nx,ny+1))
          absinflow = absinflow + sum(abs( V2D(1:nx,ny+1) ))
      endif
     !   call msg("latconvergence", sum(Disp_wind%Phi(:,:), 1))

        if (.not. ifLonGlobal) then 
          absinflow = absinflow + sum(abs( U2D(1,1:ny) )) +  sum(abs( U2D(nx+1,1:ny) ))  
          inflow =  inflow + sum(U2D(1,1:ny)) - sum(U2D(nx+1,1:ny))
       endif

       call msg("inflow", inflow)
      
       if (.not. ifSpole) then
            inflowfactor = inflow/absinflow
            V2D(1:nx,1) = inflowfactor * abs(V2D(1:nx,1))
            Disp_wind%Phi(:,1) = Disp_wind%Phi(:,1) + V2D(1:nx,1)
            call   msg("Sum of divergence inside the domain after S adj  (kg/s)",sum(Disp_wind%Phi(:,:)))
       endif
       if ( .not. ifNpole)  then 
          inflowfactor = inflow/absinflow
          V2D(1:nx,ny+1) = - inflowfactor *  abs( V2D(1:nx,ny+1) )
          Disp_wind%Phi(:,ny) = Disp_wind%Phi(:,ny) - V2D(1:nx,ny+1)
          call   msg("Sum of masschange inside the domain after N adj  (kg/s)",sum(Disp_wind%Phi(:,:)))
       endif
      if (.not. ifLonGlobal) then
          inflowfactor = inflow/absinflow
          U2D(1,1:ny)    =  inflowfactor * abs( U2D(1,1:ny) )
          U2D(nx+1,1:ny) =  - inflowfactor * abs( U2D(nx+1,1:ny) )
          Disp_wind%Phi(1,:) = Disp_wind%Phi(1,:) + U2D(1,1:ny)
          call   msg("Sum of masschange inside the domain after W adj  (kg/s)",sum(Disp_wind%Phi(:,:)))
          Disp_wind%Phi(nx,:) = Disp_wind%Phi(nx,:) - U2D(nx+1,1:ny)
          call   msg("Sum of masschange inside the domain after E adj  (kg/s)",sum(Disp_wind%Phi(:,:)))
       endif
    endif

      !
      ! Poison solver
      !
      iTmp = 1  ! x boundary -- zero (BLKTRI)
      if (ifLonGlobal) iTmp=0 ! Periodic x boundaries
     
      jTmp = 0 !Return status 
    call  BLKTRI (1,     1, ny, Disp_wind%AN, Disp_wind%BN, Disp_wind%CN, &
                      & iTmp , nx, Disp_wind%AM, Disp_wind%BM, Disp_wind%CM, &
                      & nx, Disp_wind%Phi, jTmp, Disp_wind%fish_work)
      if (jTmp /= 0) then
        call msg("BLKTRI returned", jTmp)
        call set_error ("Failed poisson solver", "df_update_cellfluxes_rt")
        return
      endif
     
      !get Velocity compensation for in-domain from potential
     
    U2D(2:nx,1:ny) = Disp_wind%Phi(2:nx,1:ny) - Disp_wind%Phi(1:nx-1,1:ny)
         if (ifLonGlobal) then ! Put correction to both  east and west side
          U2D(1,:) =  Disp_wind%Phi(1,:) -  Disp_wind%Phi(nx,:)
          U2D(nx+1,1:ny) =  U2D(1,1:ny)
         endif

    V2D(1:nx,2:ny) = (Disp_wind%Phi(1:nx,2:ny) - Disp_wind%Phi(1:nx,1:ny-1))
    do iY = 2,ny 
      V2D(1:nx,iY) = V2D(1:nx,iY) *  Disp_wind%CosLat(iY)
    enddo

    !
    !FIXME Checkup
!     call msg("Correction")
!
!    Disp_wind%Phi(:,:) =  - V2D(1:nx,1:ny) + V2D(1:nx,2:ny+1) &
!                        & - U2D(1:nx,1:ny) + U2D(2:nx+1,1:ny)
!    call msg("latconvergence", sum(Disp_wind%Phi(:,:), 1))
!    if (associated(deltaM2Drate)) then
!     !    call msg("Deltam       :", sum(deltaM2Drate(1:nx,1:ny) + m2D(1:nx,1:ny)*masscorrfactor, 1))
     !   call msg("sumconvergence", sum(Disp_wind%Phi(:,:)))
!      do iy=ny,1,-1
!          print '(I2,X,A3,X,10F10.0)', iy, "div",   Disp_wind%Phi(:,iy) 
!          print '(I2,X,A3,X,10F10.0)', iy, "d_m",   deltaM2Drate(:,iy) 
!      enddo
!   endif

   !Now u2D and v2D have corrections to column-integrated mass flux
   !normalize by colmass
    U2D(2:nx,1:ny) = U2D(2:nx,1:ny) * 2. / (m2d(1:nx-1,1:ny)+ m2d(2:nx,1:ny))
    if (ifLonGlobal) then 
       U2D(1,1:ny) = U2D(1,1:ny) * 2. / (m2d(1,1:ny)+ m2d(nx,1:ny))
       U2D(nx+1,1:ny) = U2D(1,1:ny)  
    else
       U2D(1,1:ny) = U2D(1,1:ny) / m2d(1,1:ny)  !Use just boundary mass
       U2D(nx+1,1:ny) = U2D(nx+1,1:ny) /  m2d(nx,1:ny)
    endif
    V2D(1:nx,2:ny) =  V2D(1:nx,2:ny) * 2. /  (m2d(1:nx,1:ny-1)+ m2d(1:nx,2:ny))
    V2D(1:nx,1) =  V2D(1:nx,1)  /  m2d(1:nx,1)
    V2D(1:nx,ny+1) =  V2D(1:nx,ny+1)  / m2d(1:nx,ny)

 end subroutine adjust_2D_fluxes


  !***********************************************************************************************

subroutine  exchange_wings(disp_buf)
   !
   ! Allocates disp_buf%wingsMPI if needed
     ! Fills disp_buf%wingsMPI with neighbors winds and masses
   ! iDepth, Mass/Flux, ixy, iLev, itime, iBound
   !
    IMPLICIT NONE
   
   ! Imported parameters
    TYPE(Tfield_buffer), POINTER :: disp_buf  ! buffer in which we diagnose
   
   ! Local variables
    ! Work arrays for sends
    REAL, DIMENSION(:), POINTER :: sendN_WAptr, sendS_WAptr, sendE_WAptr, sendW_WAptr
    ! Work arrays for receives
    REAL, DIMENSION(:), POINTER :: recvN_WAptr, recvS_WAptr, recvE_WAptr, recvW_WAptr
    ! Size of slices
    INTEGER :: nSlice, nxSlice, nySlice

   TYPE(silja_field), POINTER :: field
    INTEGER :: ind_m, ind_furt, ind_fvrt, ind_fu, ind_fv
    REAL, DIMENSION(:,:), POINTER :: uPtr2D, vPtr2D, mPtr2D
    INTEGER :: recv_count
    INTEGER :: iTmp, nx_disp, ny_disp, iLev, jTmp, iBound, &
         & nz_disp, nDepth, xbuf_size, ybuf_size

   ! Fluxes and masses at the grid boarder stuff
    TYPE(TRealPtr), DIMENSION(4) :: flat  !Slices of a work array
    ! Fluxes and masses at the grid boarder
    REAL, DIMENSION(:,:,:,:,:), POINTER :: wingsSendN, wingsSendS, wingsSendE, wingsSendW
    ! Length of boundaries
    INTEGER, DIMENSION(4) :: boundLength
    integer :: my_x_coord, my_y_coord, iDirection

   ! No neigbours -- no wings
   if (all(adv_mpi_neighbours(:) < 0)) return

   call smpi_get_process_xy_topology(my_x_coord, my_y_coord, iTmp, jTmp)

   !Realtime 
   ind_furt = fu_index(disp_buf%buffer_quantities,disp_flux_celle_rt_flag) 
   ind_fvrt = fu_index(disp_buf%buffer_quantities,disp_flux_celln_rt_flag) 
   ! Multitime 
   ind_m = fu_index(disp_buf%buffer_quantities,disp_cell_airmass_flag)
     ind_fu = fu_index(disp_buf%buffer_quantities,disp_flux_celleast_flag)
     ind_fv = fu_index(disp_buf%buffer_quantities,disp_flux_cellnorth_flag)
   if (any((/ind_furt, ind_fvrt, ind_m, ind_fu, ind_fv/) < 1)) then
      call msg("Failed to find quantities in buffer:")
       call msg("ind_furt, ind_fvrt, ind_m, ind_fu, ind_fv", (/ind_furt, ind_fvrt, ind_m, ind_fu, ind_fv/))

      call msg("Buffer contains:")
      do iTmp = 1,size(disp_buf%buffer_quantities)
              call msg(fu_quantity_string(disp_buf%buffer_quantities(iTmp)))
      enddo
       call set_error("Could not find input/output fields in disp buffer", "exchange_wings")
         return
   endif

   call grid_dimensions(fu_grid(disp_buf%p4d(ind_m)%past%p2d(1)%idPtr), nx_disp, ny_disp)
    nDepth = wing_depth ! Use the global maximum as buffer depth
   if (wing_depth > nx_disp .or. wing_depth > ny_disp) then
         call msg("Depth, nx, ny", (/wing_depth, nx_disp, ny_disp/))
         call set_error("Wing depth exceeds subdomain dimensions.","exchange_wings")
         return
   endif
 
   nz_disp = fu_number_of_fields(disp_buf%p4d(ind_m)%past%field3d)
   boundLength(1:4) = (/nx_disp, nx_disp, ny_disp, ny_disp/)
    if (.not. associated(disp_buf%wings_N)) then
      call msg("Allocating MPI wings for North:", (/3,2,nz_disp,nx_disp,nDepth/))
      allocate(disp_buf%wings_N(3,2,nz_disp,nx_disp,nDepth), stat=iTmp)
      if (iTmp/=0) then 
         call set_error("failed to allocate disp_buf%wingsMPI","exchange_wings")
         return
      endif
   endif
    if (.not. associated(disp_buf%wings_S)) then
      call msg("Allocating MPI wings for South:", (/3,2,nz_disp,nx_disp,nDepth/))
      allocate(disp_buf%wings_S(3,2,nz_disp,nx_disp,nDepth), stat=iTmp)
      if (iTmp/=0) then 
        call set_error("failed to allocate disp_buf%wingsMPI","exchange_wings")
        return
      endif
    endif
    if (.not. associated(disp_buf%wings_E)) then
      call msg("Allocating MPI wings for East:", (/3,2,nz_disp,nDepth,ny_disp/))
      allocate(disp_buf%wings_E(3,2,nz_disp,nDepth,ny_disp), stat=iTmp)
      if (iTmp/=0) then
        call set_error("failed to allocate disp_buf%wingsMPI", "exchange_wings")
        return
      endif
    endif
    if (.not. associated(disp_buf%wings_W)) then
      call msg("Allocating MPI wings for West:", (/3,2,nz_disp,nDepth,ny_disp/))
      allocate(disp_buf%wings_W(3,2,nz_disp,nDepth,ny_disp), stat=iTmp)
       if (iTmp/=0) then
         call set_error("failed to allocate disp_buf%wingsMPI", "exchange_wings")
         return
       endif
     endif

     nxSlice = 3 * nz_disp * nx_disp * 2 * nDepth
     nySlice = 3 * nz_disp * nDepth * 2 * ny_disp

    sendN_WAptr => fu_work_array(nxSlice)
    sendS_WAptr => fu_work_array(nxSlice)
    sendE_WAptr => fu_work_array(nySlice)
    sendW_WAptr => fu_work_array(nySlice)
    recvN_WAptr => fu_work_array(nxSlice)
    recvS_WAptr => fu_work_array(nxSlice)
    recvE_WAptr => fu_work_array(nySlice)
    recvW_WAptr => fu_work_array(nySlice)

    wingsSendN(1:3,1:2,1:nz_disp,1:nx_disp,1:nDepth) => sendN_WAptr
    wingsSendS(1:3,1:2,1:nz_disp,1:nx_disp,1:nDepth) => sendS_WAptr
    wingsSendE(1:3,1:2,1:nz_disp,1:nDepth,1:ny_disp) => sendE_WAptr
    wingsSendW(1:3,1:2,1:nz_disp,1:nDepth,1:ny_disp) => sendW_WAptr

   ! Stuff send buffer with wings... 
   do iLev = 1,nz_disp
      ! Past
       uPtr2D(1:nx_disp+1,1:ny_disp) => disp_buf%p4d(ind_fu)%past%p2d(iLev)%ptr
         vPtr2D(1:nx_disp,1:ny_disp+1) => disp_buf%p4d(ind_fv)%past%p2d(iLev)%ptr 
         mPtr2D(1:nx_disp, 1:ny_disp)  => disp_buf%p4d(ind_m)%past%p2d(iLev)%ptr

      ! North
      wingsSendN(iPast,iMass,iLev,1:nx_disp,1:nDepth) = mPtr2D(1:nx_disp, ny_disp-nDepth+1:ny_disp)
      wingsSendN(iPast,iFlux,iLev,1:nx_disp,1:nDepth) = vPtr2D(1:nx_disp, ny_disp-nDepth+1:ny_disp)
      ! South
      wingsSendS(iPast,iMass,iLev,1:nx_disp,1:nDepth) = mPtr2D(1:nx_disp, 1:nDepth  )
      wingsSendS(iPast,iFlux,iLev,1:nx_disp,1:nDepth) = vPtr2D(1:nx_disp, 2:nDepth+1)
      ! East
      wingsSendE(iPast,iMass,iLev,1:nDepth,1:ny_disp) = mPtr2D(nx_disp-nDepth+1:nx_disp, 1:ny_disp)
      wingsSendE(iPast,iFlux,iLev,1:nDepth,1:ny_disp) = uPtr2D(nx_disp-nDepth+1:nx_disp, 1:ny_disp)
      ! West
      wingsSendW(iPast,iMass,iLev,1:nDepth,1:ny_disp) = mPtr2D(1:nDepth,   1:ny_disp)
      wingsSendW(iPast,iFlux,iLev,1:nDepth,1:ny_disp) = uPtr2D(2:nDepth+1, 1:ny_disp)

         !Future
       uPtr2D(1:nx_disp+1,1:ny_disp) => disp_buf%p4d(ind_fu)%Future%p2d(iLev)%ptr
         vPtr2D(1:nx_disp,1:ny_disp+1) => disp_buf%p4d(ind_fv)%Future%p2d(iLev)%ptr 
         mPtr2D(1:nx_disp, 1:ny_disp)  => disp_buf%p4d(ind_m)%Future%p2d(iLev)%ptr

      ! North
      wingsSendN(iFuture,iMass,iLev,1:nx_disp,1:nDepth) = mPtr2D(1:nx_disp, ny_disp-nDepth+1:ny_disp)
      wingsSendN(iFuture,iFlux,iLev,1:nx_disp,1:nDepth) = vPtr2D(1:nx_disp, ny_disp-nDepth+1:ny_disp)
      ! South
      wingsSendS(iFuture,iMass,iLev,1:nx_disp,1:nDepth) = mPtr2D(1:nx_disp, 1:nDepth  )
      wingsSendS(iFuture,iFlux,iLev,1:nx_disp,1:nDepth) = vPtr2D(1:nx_disp, 2:nDepth+1)
      ! East
      wingsSendE(iFuture,iMass,iLev,1:nDepth,1:ny_disp) = mPtr2D(nx_disp-nDepth+1:nx_disp, 1:ny_disp)
      wingsSendE(iFuture,iFlux,iLev,1:nDepth,1:ny_disp) = uPtr2D(nx_disp-nDepth+1:nx_disp, 1:ny_disp)
      ! West
      wingsSendW(iFuture,iMass,iLev,1:nDepth,1:ny_disp) = mPtr2D(1:nDepth,   1:ny_disp)
      wingsSendW(iFuture,iFlux,iLev,1:nDepth,1:ny_disp) = uPtr2D(2:nDepth+1, 1:ny_disp)

      ! Realtime
      wingsSendN(iRealTime,iFlux,iLev,:,:) = CONST_NAN
      wingsSendS(iRealTime,iFlux,iLev,:,:) = CONST_NAN
      wingsSendE(iRealTime,iFlux,iLev,:,:) = CONST_NAN
      wingsSendW(iRealTime,iFlux,iLev,:,:) = CONST_NAN
       uPtr2D(1:nx_disp+1,1:ny_disp) => disp_buf%p4d(ind_furt)%present%p2d(iLev)%ptr
         vPtr2D(1:nx_disp,1:ny_disp+1) => disp_buf%p4d(ind_fvrt)%present%p2d(iLev)%ptr 

      ! North
      wingsSendN(iRealTime,iFlux,iLev,1:nx_disp,1:nDepth) = vPtr2D(1:nx_disp, ny_disp-nDepth+1:ny_disp)
      ! South
      wingsSendS(iRealTime,iFlux,iLev,1:nx_disp,1:nDepth) = vPtr2D(1:nx_disp, 2:nDepth+1)
      ! East
      wingsSendE(iRealTime,iFlux,iLev,1:nDepth,1:ny_disp) = uPtr2D(nx_disp-nDepth+1:nx_disp, 1:ny_disp)
      ! West
      wingsSendW(iRealTime,iFlux,iLev,1:nDepth,1:ny_disp) = uPtr2D(2:nDepth+1, 1:ny_disp)
     enddo

     ! Exchange the boundaries in North-South direction
    do iDirection = 0,1
      if (mod(my_y_coord, 2) == iDirection) then
        if(adv_mpi_neighbours(northern_boundary)>=0)then
          call smpi_exchange_wings(adv_mpi_neighbours(northern_boundary), sendN_WAptr, nxSlice, &
               & recvN_WAptr, recv_count)
          disp_buf%wings_N(1:3,1:2,1:nz_disp,1:nx_disp,1:nDepth) = &
               & reshape(recvN_WAptr(1:nxSlice), (/3,2,nz_disp,nx_disp,nDepth/))
        else
          disp_buf%wings_N(1:3,1:2,1:nz_disp,1:nx_disp,1:nDepth) = -8888
        endif
      else
        if(adv_mpi_neighbours(southern_boundary)>=0)then
          call smpi_exchange_wings(adv_mpi_neighbours(southern_boundary), sendS_WAptr, nxSlice, &
               & recvS_WAptr, recv_count)
          disp_buf%wings_S(1:3,1:2,1:nz_disp,1:nx_disp,1:nDepth) = &
               & reshape(recvS_WAptr(1:nxSlice), (/3,2,nz_disp,nx_disp,nDepth/))
        else
          disp_buf%wings_S(1:3,1:2,1:nz_disp,1:nx_disp,1:nDepth) = -8888
        endif
      endif
    enddo

     ! Exchange the boundaries in East-West direction
    do iDirection = 0,1 
      if (mod(my_x_coord, 2) == iDirection) then
        if(adv_mpi_neighbours(eastern_boundary)>=0)then
          call smpi_exchange_wings(adv_mpi_neighbours(eastern_boundary), sendE_WAptr, nySlice, &
               & recvE_WAptr, recv_count)
          disp_buf%wings_E(1:3,1:2,1:nz_disp,1:nDepth,1:ny_disp) = &
               & reshape(recvE_WAptr(1:nySlice), (/3,2,nz_disp,nDepth,ny_disp/))
        else
          disp_buf%wings_E(1:3,1:2,1:nz_disp,1:nDepth,1:ny_disp) = -8888
        endif
      else
        if(adv_mpi_neighbours(western_boundary)>=0)then
          call smpi_exchange_wings(adv_mpi_neighbours(western_boundary), sendW_WAptr, nySlice, &
               & recvW_WAptr, recv_count)
          disp_buf%wings_W(1:3,1:2,1:nz_disp,1:nDepth,1:ny_disp) = &
               & reshape(recvW_WAptr(1:nySlice), (/3,2,nz_disp,nDepth,ny_disp/))
        else
          disp_buf%wings_W(1:3,1:2,1:nz_disp,1:nDepth,1:ny_disp) = -8888
        endif
      endif
    enddo

    call free_work_array(sendN_WAptr)
    call free_work_array(sendS_WAptr)
    call free_work_array(sendE_WAptr)
    call free_work_array(sendW_WAptr)

    call free_work_array(recvN_WAptr)
    call free_work_array(recvS_WAptr)
    call free_work_array(recvE_WAptr)
    call free_work_array(recvW_WAptr)

end subroutine  exchange_wings


END MODULE diagnostic_variables
