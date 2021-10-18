MODULE diagnostic_variables
  !
  ! This module contains a set of routines for deriving a diagnostic dispersion
  ! variables. 
  ! Do not mix with derived meteorological quantities: those are 
  ! input fields, while these are the outputs of dispersion simulations. By definition,
  ! diagnostic variables are those, which are not predicted within the main cycle of
  ! the simualtions but rather derived from its products. They also do not affect the
  ! further model simualtions, thus being entirely output-oriented.
  ! The only reason why these routines are not moved outside the SILAM code is that
  ! their computation might need internal variables, which are either
  ! not easily available for output or bulky.
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
  use fishpack_silam
  !$use omp_lib

  implicit none

  public make_all_diagnostic_fields
  public fu_ifDiagnosticQuantity
  public init_wind_diag
  public init_water_in_soil_model
  public ifUseMassfluxes
  public diagnostic_input_needs

  private make_diagnostic_fields
  private diagnose_single_dyn_field
  private diag_vertical_wind_incompr
  private diagnostic_wind_test
  private diag_vertical_wind_incompr_v2
  private diag_vertical_wind_anelastic
  private diag_vertical_wind_eta
  private df_cnc_from_vmr
  private df_DMAT_Vd_correction
  private df_make_cell_size_z_dyn
  private df_cumul_daily_tempr_2m

  private df_update_realtime_met_fields
  private df_daily_mean_tempr_2m

  private diag_cell_fluxes  ! Get massfluxes for past and future  ignores pressure tendency
  private df_update_cellfluxcorr_rt   ! realtime spt correction for meteo interval

  private exchange_wings
  private TBLKTRI
  private adjust_2D_fluxes
  
  


  type(THorizInterpStruct), pointer, save, private :: pHorizInterpU => null(), pHorizInterpV => null(), &
       & pHorizInterpRho => null()
  type(TVertInterpStruct), pointer, save, private ::  pVertInterpRho => null()

  type(silam_vertical), private, save :: integration_vertical
  type(TVertInterpStruct), pointer, save, private :: wind_interp_struct
  integer, dimension(:), allocatable, private, save :: nsmall_in_layer

  integer, parameter, private :: test_wind = 50000, incompressible = 50001, incompressible_v2 = 50002, &
                        & anelastic_v2 = 50003, from_nwp_omega = 50004, hybrid_top_down = 50005, & 
                        & hardtop = 50006, opentop = 50007, omegatop = 50008, topdown = 50009, hardtop_weighted = 50010

  real(vm_p), dimension(:,:), allocatable, private, save :: u3d, v3d, eta_dot_3d, w3d
  
  type Tdiagnostic_rules
     private
     integer :: continuity_equation = incompressible_v2
     integer :: wind_method = int_missing  ! Do not make any mass fluxes
     integer :: iTestWindType = 1
     integer :: a = int_missing, b = int_missing !Test Wind parameters
     logical :: w_on_half_levels, if_DoStatic
     logical :: ifAllowColdstartDayTemperature=.False.
     logical :: if_high_stomatal_conductance=.True.
     type(Trules_water_in_soil) :: rulesWIS
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
     real, dimension (:), pointer :: CosLat => Null() ! 1:ny_disp+1  
     real, dimension (:), pointer :: LevelCorrWeight => Null() ! 1:nz_disp
     !Stuff for BLKTRI
     real(fish_kind), dimension (:), pointer :: AM => Null(), BM => Null(), CM => Null() ! 1:ny_disp 
     real(fish_kind), dimension (:), pointer :: AN => Null(), BN => Null(), CN => Null() ! 1:nx_disp Stuff for BLKTRI
     ! Excess wind divergence / Velocity potential 
     real(fish_kind), dimension (:,:), pointer :: Phi=> Null()  ! 1:nx_disp, 1:ny_disp



     ! Arrays used for test wind
     real(r8k), dimension(:), pointer :: StreamFunction1D => Null() !! 0:ny_disp*nx_disp
     real(r8k), dimension(:,:), pointer :: StreamFunction => Null() !! 0:ny_disp, 0:nx_disp !Slice of above


     real, dimension(:,:), pointer :: WholeMapU => Null(), WholeMapV => Null(), WholeMapM => Null()
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


        
    if (.not. associated(Disp_wind%StreamFunction1D)) then
       !
       ! Stuff needed for local domain diagnostics
       call lonlat_grid_parameters(whole_disp_grid,&
                                  & corner_lon_E, corner_lat_N,corner_in_geo_lonlat, &
                                  & nx, ny, &
                                  & southpole_lon_E, southpole_lat_N, &
                                  & dx_deg, dy_deg)

       allocate( Disp_wind%StreamFunction1D((nx+1)*(ny+1)), & ! Needed only for test winds....
               & stat=istat)
       if (istat /= 0) then
              call set_error("Failed to allocate DispWindStruct", "init_wind_diag")
       endif
       Disp_wind%StreamFunction(1:nx+1,1:ny+1) => Disp_wind%StreamFunction1D(1:(nx+1)*(ny+1))

       ! Global domain diagnostics
       if (rules%wind_method == hardtop .or. rules%wind_method == hardtop_weighted) then
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
    endif !allocation needed

  end subroutine init_wind_diag

  
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

  subroutine set_diagnostic_rules(nlMeteoSetup, nlStandardSetup, rules)
    !
    ! Set the vertical wind method and the half_levels option.
    ! 
    implicit none
    type(Tsilam_namelist), pointer :: nlMeteoSetup, nlStandardSetup
    type(Tdiagnostic_rules), intent(out) :: rules
    integer :: iTmp

    character(len=255) :: method
    character (len=*), parameter :: sub_name = "set_diagnostic_rules"

    rules%w_on_half_levels = .true.
    call msg('Vertical velocity computed at half-levels')
    call msg('')

    method = fu_str_l_case(fu_content(nlStandardSetup, 'continuity_equation'))
    iTmp = INDEX(trim(method),' ')
    if (iTmp > 0) then 
        method = method(1:iTmp-1)
    endif
    select case (method)

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
          method = fu_str_l_case(fu_content(nlStandardSetup, 'continuity_equation'))
          method = trim(method(iTmp+1:))
          read(unit=method, fmt=*, iostat=iTmp)  rules%iTestWindType
          iTmp = INDEX(trim(method),' ')
          if  (iTmp > 0)then
            read(unit=method, fmt=*, iostat=iTmp) rules%a
            method = trim(method(iTmp+1:))
          endif
          iTmp = INDEX(trim(method),' ')
          if  (iTmp > 0)then
            read(unit=method(iTmp+1:), fmt=*, iostat=iTmp) rules%b
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
        call msg('continuity_equation = ' // trim(method))
        call msg_warning('Strange or no continuity_equation defined in standard setup,' + &
                       & 'will use differential form v2',sub_name)
        rules%continuity_equation = incompressible_v2
    end select
    

    ! To replace continuity_equation at some point
    method = fu_str_l_case(fu_content(nlStandardSetup, 'wind_diagnostics'))
    iTmp = INDEX(trim(method),' ')
    if (iTmp > 0) then 
        method = method(1:iTmp-1)
    endif
    select case (method)

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
          method = fu_str_l_case(fu_content(nlStandardSetup, 'wind_diagnostics'))
          method = trim(method(iTmp+1:))
          read(unit=method, fmt=*, iostat=iTmp)  rules%iTestWindType
          iTmp = INDEX(trim(method),' ')
          if  (iTmp > 0)then
            read(unit=method, fmt=*, iostat=iTmp) rules%a
            method = trim(method(iTmp+1:))
          endif
          iTmp = INDEX(trim(method),' ')
          if  (iTmp > 0)then
            read(unit=method(iTmp+1:), fmt=*, iostat=iTmp) rules%b
          endif
        else
          call set_error('test_wind requires at least 1 parameter: wind type index',&
                           &sub_name)
          return
        endif
        call msg("Test wind parameters", rules%a, rules%b)
        call msg_warning('Fake massflux wind:'+fu_str(rules%iTestWindType))
      case default
        call msg('wind_diagnostics = ' // trim(method))
        call msg("can be wind_diagnostics = opentop / hardtop / none")
        call msg_warning('Strange or no wind_diagnostics defined in standard setup,' + &
                       & 'will use winds opentop',sub_name)
        rules%wind_method = opentop
    end select

    method = fu_str_l_case(fu_content(nlStandardSetup, 'stomatal_conductance'))
    select case (method)
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

    if (fu_str_l_case(fu_content(nlStandardSetup, 'allow_coldstart_day_temperature'))=='yes') &
      rules%ifAllowColdstartDayTemperature=.True.

    !
    ! Seemingly a right place to set the soil water diagnostic rules
    ! Since here we have no clue if this is needed at all, just set them - this is always possible
    ! Later it will be decided if we want them
    !
    call set_rules_WIS(nlMeteoSetup, rules%rulesWIS)
    
  end subroutine set_diagnostic_rules

  
  !****************************************************************************************
  
  subroutine diagnostic_input_needs(rules, q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st)
    !
    ! Includes the meteorological quantities needed for the diagnostics.
    ! The arrays already have all requestes from specific modules, so we can check what is
    ! needed to diagnose
    !
    implicit none
    
    ! Imported parameters
    type(Tdiagnostic_rules), intent(inout) :: rules
    INTEGER, DIMENSION(:), INTENT(inout) :: q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st
    
    ! Local variables
    integer :: iTmp

    ! The basic parameters that will be needed always
    !
    iTmp = fu_merge_integer_to_array(height_flag, q_met_dynamic)
    iTmp = fu_merge_integer_to_array(surface_pressure_flag, q_met_dynamic)

    !
    ! Part of the game: water in soil. Not exactly diagnostic but still quite nice here
    ! Call it only if the corresponding quantities are needed
    !
    if(fu_quantity_in_quantities(soil_moisture_vol_frac_nwp_flag, q_met_dynamic) .or. &
     & fu_quantity_in_quantities(water_in_soil_srf_grav_flag, q_met_dynamic) .or. &
     & fu_quantity_in_quantities(water_in_soil_deep_grav_flag, q_met_dynamic))then
      call add_input_needs_WIS(rules%rulesWIS, q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st)
    else
      rules%rulesWIS = rulesWIS_missing  ! a precaution. Actually, should be set elsewhere too
    endif
    
    !
    ! For each sub here, there is a requested quantities and input ones. Then:
    !
    if(fu_quantity_in_quantities(day_mean_temperature_flag, q_disp_st)) &
                              & iTmp = fu_merge_int_to_array(day_temperature_acc_flag, q_disp_dynamic)
    if(fu_quantity_in_quantities(day_mean_temperature_2m_flag, q_disp_st)) &
                              & iTmp = fu_merge_int_to_array(day_temperature_2m_acc_flag, q_disp_dynamic)
    if(fu_quantity_in_quantities(day_temperature_acc_flag, q_disp_dynamic)) &
                              & iTmp = fu_merge_int_to_array(temperature_flag, q_met_dynamic)
    if(fu_quantity_in_quantities(day_temperature_2m_acc_flag, q_disp_dynamic)) &
                              & iTmp = fu_merge_int_to_array(temperature_2m_flag, q_met_dynamic)

  end subroutine diagnostic_input_needs


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
                                      & now, timestep, &
                                      & ifRunDispersion)
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
    logical, intent(in) :: ifNewMeteoData, ifRunDispersion
!    !FIXME delete this field
!    type(silja_field), pointer :: field
    logical :: iffound
    
    ! Local variables
    type(silja_time), save :: valid_border_time_1 = really_far_in_past, &
                            & valid_border_time_2 = really_far_in_past
    type(silja_interval) :: meteo_time_shift

    !
    ! If meteo is shifted, take this into account
    !
    meteo_time_shift = fu_meteo_time_shift(wdr)
    !
    ! Should be called before make_derived_meteo_fields to get tz_offset field
    ! properly
    !
    call update_timezone_offsets(now)
    call SolarSetup(now)  ! Solar zenith angle, sunset, sunrize etc would be right
   
    !-----------------------------------------------------------------------------
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


    !------------------------------------------------------------------------------
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

    ! Probably, here is the right point to take care of water in soil
    !
    if(defined(diagnostic_rules%rulesWIS))then
      call update_water_in_soil(meteo_buf_ptr, disp_buf_ptr, now + meteo_time_shift, timestep, diagnostic_rules%rulesWIS)
      if(error)return
    endif
    
    if(ifRunDispersion)then
      !
      ! Finally, derive the dispersion diagnostic fields using their shopping lists
      !
      call make_diagnostic_fields(meteo_buf_ptr, & ! meteo buffer
                                & dispersionMarketPtr, & ! disp market
                                & disp_dyn_shopping_list, & ! full list of disp dynamic fields
                                & disp_stat_shopping_list, & ! full list of disp static fields
           !                     & now, &
                                & diagnostic_rules)
      if(error)return
      call msg_test('Setting dispersion pointers...')
      call arrange_buffer(dispersionMarketPtr, met_src_missing, now + meteo_time_shift, timestep, & 
                        & disp_buf_ptr,  dispersion_verticalPtr)
      if(error)return
    endif   ! ifRunDispersion

    ! fields must be valid for the middle of step
    ! The sub prepares the fields with certain validity period. As soon as out target time is outside
    ! we call it again.
    !
    if(.not. fu_between_times(now + meteo_time_shift + timestep * 0.5, &
                            & valid_border_time_1, valid_border_time_2, .true.))then
!      call msg("Updating realime for time:"+fu_str(now + timestep * 0.5))
!      call msg("Valid range = "+fu_str(valid_border_time_1)+ ":" +fu_str(valid_border_time_2))

      call df_update_realtime_met_fields(meteoMarketPtr, meteo_buf_ptr, &
                                       & now + meteo_time_shift + timestep * 0.5, wdr, diagnostic_rules,&
                                       & valid_border_time_1, valid_border_time_2, .true.)   ! interval of validity, output
      if(error)return

      if(ifRunDispersion)then 
        call df_update_realtime_met_fields(dispersionMarketPtr, disp_buf_ptr, &
                                       & now + meteo_time_shift + timestep * 0.5, wdr, diagnostic_rules, &
                                       & valid_border_time_1, valid_border_time_2, .false.)   
                                    ! interval of validity, output. Do not reset it
      endif

      if(error)return

    endif
    
!      call msg("**********************Diagnozing output market with shopping list:")
!      call report(output_dyn_shopping_list)
!      call msg("********************")

    call make_diagnostic_fields(meteo_buf_ptr,   &   ! meteo buffer
                              & outputMarketPtr, &   ! output market
                              & output_dyn_shopping_list,  &  ! full list of output dynamic fields
                              & output_stat_shopping_list, & ! full list of output static fields
                            !  & now, &
                              & diagnostic_rules)
    if(error)return

    
!    call msg_test('Setting output pointers...')
!    call arrange_buffer(outputMarketPtr, met_src_missing, now, & 
!                      & weight_past, output_buf_ptr, output_gridPtr, output_verticalPtr)
!    if(error)return

    diagnostic_rules%if_DoStatic = .false.  ! no need to repeat it at each time step

  end subroutine make_all_diagnostic_fields


  !***************************************************************************************************

  subroutine make_diagnostic_fields(meteo_ptr, MarketPtr, &
                                  & dqListDyn, dqListStat,&
!                                  & ifMeteo2DispHorizInterp, ifMeteo2DispVertInterp, &
!                                  & pMeteo2DispHorizInterp, pMeteo2DispVertInterp, &
!                                  & now, &
                                  & diagnostic_rules)
    !
    ! Creates the quantities in list in every stack in the minimarket, where possible
    !
    implicit none

    ! Imported variables
    type(silja_shopping_list), intent(in) ::  dqListDyn, dqListStat
    type(mini_market_of_stacks), pointer :: MarketPtr
    type(Tfield_buffer), pointer :: meteo_ptr
!    type(THorizInterpStruct), pointer :: pMeteo2DispHorizInterp
!    type(TVertInterpStruct), pointer :: pMeteo2DispVertInterp
!    logical :: ifMeteo2DispHorizInterp, ifMeteo2DispVertInterp
!    type(silja_time), intent(in) :: now
    type(Tdiagnostic_rules), intent(in) :: diagnostic_rules

    ! Local variables
    type(meteo_data_source), dimension(max_met_srcs) :: metSrcsMultitime, metSrcsSingletime
    TYPE(silja_time), DIMENSION(2) :: obstimes
    integer :: nMetSrcsMultitime, nMetSrcsSingletime, ntimes, iTmp
!    type(silam_shopping_variable) :: shopVar, varTmp
    type(silja_stack), pointer :: stackPtr 
    type(silja_field), pointer :: fieldPtr
    type(silja_field_id), pointer :: idIn
    type(silja_field_id) :: idRequest
    integer :: iVarLst, iMetSrc, iT, iFldStack, shopQ, iVar, nVars
    integer, dimension(max_quantities) :: iArr
    logical :: iffound,  ifvalid

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

    !
    ! Dynamic market
    !
    !cell_size_z_flag must be diagnozed before all others
    nVars = fu_nbr_of_quantities(dqListDyn,.false.)
    do iVar = 1, nVars
      iArr(iVar) = fu_quantity(dqListDyn, iVar)
      if (iArr(iVar) == cell_size_z_flag) then
              iArr(iVar) = iArr(1) 
              iArr(1) = cell_size_z_flag
      endif
    enddo
    !
    ! The main cycle over the diagnostic variables
    !
    do iVarLst = 1, nVars

      obstimes(1:2) = time_missing
      shopQ = iArr(iVarLst)
      ! day_temperature_2m_acc_flag has all the intelligence inside:
      ! future has to be rediagnozed after past initialized
      ! 
      if ( shopQ /=  day_temperature_2m_acc_flag) then
        iTmp = 0
        if(.not. fu_field_in_sm(MarketPtr, & ! Past already done before?
                                & met_src_missing, & ! metSrcsMultitime(iMetSrc), &
                                & shopQ, & !fu_quantity(shopVar), &
                                & meteo_ptr%time_past, &
                                & level_missing, &
                                & fu_multi_level_quantity(shopQ), & !.true., &
                                & permanent=.false.))then
            iTmp = iTmp + 1
            obstimes(iTmp) = meteo_ptr%time_past
        endif

        if(.not. fu_field_in_sm(MarketPtr, & ! Future already done before?
                                & met_src_missing, & ! metSrcsMultitime(iMetSrc), &
                                & shopQ, & !fu_quantity(shopVar), &
                                & meteo_ptr%time_future, &
                                & level_missing, &
                                & fu_multi_level_quantity(shopQ), & !.true., &
                                & permanent=.false.))then
            iTmp = iTmp + 1
            obstimes(iTmp) = meteo_ptr%time_future
        endif

        if(error)then
          call msg_warning('Something went wrong with detection of already-done dynamic fields but continue')
          call unset_error('make_diagnostic_fields')
        endif
        
        ! Dirty hack: force diagnostics on a second time (right after initalisation)
        ! Obs! day_temperature_2m_acc_flag does not care about obstimes
        if(iTmp==0) cycle ! Nothing to diagnoze

      endif
      ! Diagnose the quantity
      !
      call diagnose_single_dyn_field(meteo_ptr, &
                                   & MarketPtr, &
                                   & obstimes, &
                                   & shopQ, &  ! quantity to be diagnosed
                                   & diagnostic_rules)
      if(error)return
      call arrange_supermarket_multitime(MarketPtr)
      if(error)return
    end do  ! varLstDyn


!    call msg('End of diagnosing:' + fu_name(MarketPtr))

  end subroutine make_diagnostic_fields


  !*********************************************************************

  subroutine diagnose_single_dyn_field(met_buf, &
                                     & dispMarketPtr, &
                                     & obstimes, &
                                     & desiredQ, &
                                     & rules)
    !
    ! An actual dispatcher for diagnosing the fields in the whatever buffer.
    ! Gets an desired quantity form the shopping list and, if needed, the input
    ! quantity from where it needs to be diagnosed and calls the corresponding sub
    !
    implicit none

    ! Imported parameters
    type(Tfield_buffer), pointer :: met_buf              ! input meteo data: buffer!
    type(mini_market_of_stacks), pointer :: dispMarketPtr   ! working mini-market
    type(silja_time), dimension(:), intent(in) :: obstimes
    integer, intent(in) :: desiredQ                 ! input and output quantities wanted
    type(Tdiagnostic_rules), intent(in) :: rules

    ! Local variables
    type(THorizInterpStruct), pointer :: pHorizInterpStruct ! meteo 2 dispersion horizontal
    type(TVertInterpStruct), pointer :: pVertInterpStruct   ! meteo 2 dispersion vertical
    logical :: ifHorizInterp, ifVertInterp      ! if interpolation needed

    select case (desiredQ)

!      case(concentration_flag)
!        !
!        ! Concentration from mixing ratio
!        !
!        if(inputQ == volume_mixing_ratio_flag)then
!          call df_cnc_from_vmr(met_buf, weight_past, &
!                             & pHorizInterpStruct, pVertInterpStruct, ifHorizInterp, ifVertInterp, &
!                             & in_fld_3d, out_fld_3d)
!        else
!          call set_error(fu_connect_strings('Cannot diagnose concentration from: ', &
!                                          & fu_quantity_string(inputQ)), &
!                       & 'diagnose_single_field')
!        endif

      case(dispersion_u_flag, dispersion_v_flag, dispersion_w_flag)
        !
        ! Solenoidal wind field in dispersion grid
        !
        call msg('Diagnosing dispersion vertical wind')
        !
        ! Get the interpolation structures
        !
        if(rules%continuity_equation == test_wind)then
          
          call diagnostic_wind_test(met_src_missing, met_buf, obstimes, dispMarketPtr, rules%iTestWindType)
          
        else
          ifHorizInterp = .not. (meteo_grid == dispersion_grid)
          ifVertInterp = .not. fu_cmp_verts_eq(meteo_vertical, dispersion_vertical)
          if(ifHorizInterp) &
             & pHorizInterpStruct => fu_horiz_interp_struct(meteo_grid, dispersion_grid, linear, .true.)
          if(ifVertInterp) &
             & pVertInterpStruct => fu_vertical_interp_struct(meteo_vertical, dispersion_vertical, &
                                                            & dispersion_grid, linear, &
                                                            & one_hour, 'wind_meteo_to_disp')
          select case(fu_leveltype(dispersion_vertical))
            case(layer_btw_2_hybrid)                                         ! hybrid levels
              call diag_vertical_wind_eta(met_src_missing, &
                                        & met_buf, &
                                        & obstimes, &
                                        & ifHorizInterp, pHorizInterpStruct, &
                                        & ifVertInterp, pVertInterpStruct, &
                                        & dispMarketPtr,&
                                        & rules%w_on_half_levels, &
                                        & rules%continuity_equation == from_nwp_omega, &
                                        & rules%continuity_equation == hybrid_top_down)
            case(layer_btw_2_height)                                         ! height levels
              select case(rules%continuity_equation)
                case(incompressible_v2)
                  call diag_vertical_wind_incompr_v2(met_src_missing, &
                                                   & met_buf, &
                                                   & obstimes, &
                                                   & ifHorizInterp, pHorizInterpStruct, &
                                                   & ifVertInterp, pVertInterpStruct, &
                                                   & dispMarketPtr,&
                                                   & rules%w_on_half_levels)
                case(incompressible)
                  call diag_vertical_wind_incompr(met_src_missing, &
                                                & met_buf, &
                                                & obstimes, &
                                                & ifHorizInterp, pHorizInterpStruct, &
                                                & ifVertInterp, pVertInterpStruct, &
                                                & dispMarketPtr,&
                                                & rules%w_on_half_levels)
                case(anelastic_v2)
                  call diag_vertical_wind_anelastic(met_src_missing, &
                                                  & met_buf, &
                                                  & obstimes, &
                                                  & ifHorizInterp, pHorizInterpStruct, &
                                                  & ifVertInterp, pVertInterpStruct, &
                                                  & dispMarketPtr,&
                                                  & rules%w_on_half_levels)
                case default
                  call set_error('Unsupported continuity equation' + fu_str(rules%continuity_equation), &
                               & 'diagnose_single_dyn_field')
                  return
              end select  ! continuity equation
            case default
              call set_error('Unsupported dispersion vertical type', 'diagnose_single_dyn_field')
          end select ! leveltype
        endif

      case(Vd_correction_DMAT_flag)
        !
        ! V_d correction for DMAT/acid chem_dep.
        !
        call msg('Diagnosing Vd correction')
        !
        ! Get the interpolation structures
        !
        ifHorizInterp = .not. meteo_grid == dispersion_grid
        if(ifHorizInterp) &
                & pHorizInterpStruct => fu_horiz_interp_struct(meteo_grid, dispersion_grid, linear, .true.)
        if(error)return
        call df_DMAT_Vd_correction(met_src_missing, &
                                 & met_buf, &
                                 & obstimes, &
                                 & ifHorizInterp, pHorizInterpStruct, &
                                 & dispMarketPtr)

      case(cell_size_z_flag, surface_pressure_flag)
        !
        ! Dynamic thickness of the layers.
        !
        call msg("Diagnosing dynamic dispersion vertical cell size..")
        !
        ! Get the interpolation structures
        !
        ifHorizInterp = .not. meteo_grid == dispersion_grid
        ifVertInterp = .not. fu_cmp_verts_eq(meteo_vertical, dispersion_vertical)
        if(ifHorizInterp) &
          & pHorizInterpStruct => fu_horiz_interp_struct(meteo_grid, dispersion_grid, linear, .true.)

        if(ifVertInterp) &
             & pVertInterpStruct => fu_vertical_interp_struct(meteo_vertical, dispersion_vertical, &
                                                            & dispersion_grid, linear, &
                                                            & one_hour, 'main_meteo_to_disp')
        call df_make_cell_size_z_dyn(met_src_missing, &
                                   & met_buf, &
                                   & obstimes, &
                                   & pHorizInterpStruct, ifHorizInterp, &
                                   & pVertInterpStruct, ifVertInterp, &
                                   & dispMarketPtr, dispersion_grid, dispersion_vertical)

        if(error)return

      case(disp_flux_celleast_flag,disp_flux_cellnorth_flag, disp_flux_celltop_flag,  &
              & disp_cell_airmass_flag, air_density_flag)

        call diag_cell_fluxes(met_src_missing, &
                                   & met_buf, &
                                   & obstimes, &
                                   & dispMarketPtr, dispersion_grid, dispersion_vertical, rules)

      case(day_temperature_2m_acc_flag)
        !
        ! Make daily cumulative temperature, 2m
        !
        call df_cumul_daily_tempr_2m(  met_buf,  dispMarketPtr, rules%ifAllowColdstartDayTemperature)
                
      case default
        call msg_warning('Cannot diagnose dyn:' + fu_quantity_string(desiredQ) + ', skipping...', &
                       & 'diagnose_single_dyn_field')

      

    end select

  end subroutine diagnose_single_dyn_field


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
  
  subroutine diagnostic_wind_test(met_src, metBufPtr, obstimes, dispMarketPtr, iWindType)

    ! Fake diagnostics sets u v and w to zero
    !
    ! The fields are stored into dispersion market.
    !
    implicit none

    ! Units: SI 

    type(meteo_data_source), intent(in) :: met_src
    type(Tfield_buffer), pointer :: metBufPtr
    type(silja_time), dimension(:), intent(in) :: obstimes
    integer, intent(in) :: iWindType
    type(mini_market_of_stacks), intent(inout) :: dispMarketPtr

    ! Local variables
    type(silja_time) :: now
    type(silja_time) :: time, analysis_time
    type(silja_time), dimension(max_times) :: dispMarketTimes
    real, dimension(:), pointer :: zero
    integer :: iTime, iLev, i,  it, n_dispMarketTimes, ix, iy, i1d
    real :: weight_past, fTmp
    type(silja_field_id) :: idTmp !, meteoWindID
    type(silja_interval) :: forecast_length
    
    real :: pi_t_over_T, tfactor, lon, lat, cosTheta, sin2Theta, &
             & lambda_prime, sinlprime2, sin2lprime, sinTheta, cosLambda, sinLambda

    integer, parameter ::  T =  12*3600 ! seconds -- period 
    real, dimension(:), pointer :: pXCellSize, pYCellSize, u, v, w

    REAL  :: corner_lat_N,corner_lon_E,southpole_lat_N,southpole_lon_E, dx_deg, dy_deg
    INTEGER :: number_of_points_x, number_of_points_y
    LOGICAL :: corner_in_geo_lonlat
    real ::  lambda0, dlambda, theta0, dtheta, lambda, theta, dx_m, dy_m, x_centre, y_centre, speedScale

    u => fu_work_array()
    v => fu_work_array()
    w => fu_work_array()
    
    call supermarket_times(dispMarketPtr, met_src_missing, dispMarketTimes, n_dispMarketTimes)
    
    forecast_length = zero_interval

    do itime = 1, size(obstimes)
      now = obstimes(itime)
      if (.not. defined(now)) exit
      if (now == metBufPtr%time_past) then
        weight_past = 1.0
      else if (now == metBufPtr%time_future) then
        weight_past = 0.0
      else
        call set_error('Strange time:'+fu_str(obstimes(itime)),'diagnostic_wind_test')
        return
      end if
      
      analysis_time = now 

      if (error) return

      pi_t_over_T =( pi * modulo(nint(fu_sec(now  - ref_time_01012000)), 2*T)) / T

      call msg("Phase (seconds, of total):", (modulo(nint(fu_sec(now - ref_time_01012000)), 2*T)  ), T)
      tfactor = earth_radius/T*cos(pi_t_over_T)

      pXCellSize => fu_grid_data(dispersion_cell_x_size_fld)
      pYCellSize => fu_grid_data(dispersion_cell_y_size_fld)

      call msg("Producing fake wind: iWindType, param", iWindType, T/24./3600. )
      do iLev = 1, nz_dispersion
         select case(iWindType)
           case(0)                               ! 0: zero wind
             w = 0.0
             u = 0.0
             v = 0.0
             call msg('Linear test wind (u,v,w):',(/u(1),v(1),w(1)/))

           case(1)                               ! 1: linear const wind, m/s
             w = 0.0
             u = 100.0
             v = -75.
             call msg('Linear test wind (u,v,w):',(/u(1),v(1),w(1)/))

           case(2)                               ! 2: linear varying-Courant wind
             w = 0.
             do iy = 1, ny_dispersion
               do ix = 1, nx_dispersion
                 i1d = ix + (iy-1) * nx_dispersion
                 u(i1d) = 50. * (1.1 + sin(real(ix) / 15.))
               end do
             end do
             v = 0.
             call msg('Linear varying-Courant test wind (~u_max,v,w):',(/u(23),v(1),w(1)/))
             
           case(31)                                ! 3: circular wind. Unit: m/s
             w = 0.
             do iy = 1, ny_dispersion
               do ix = 1, nx_dispersion
                 i1d = ix + (iy-1) * nx_dispersion
                 u(i1d) = real(iy-(ny_dispersion/2)) / (min(nx_dispersion, ny_dispersion) / 4)
                 v(i1d) = real((nx_dispersion/2)-ix) / (min(nx_dispersion, ny_dispersion) / 4)
   !              u(i1d) = real(iy-200) / ny_dispersion * pXCellSize(i1d) /100.
   !              v(i1d) = real(200-ix) / nx_dispersion * pYCellSize(i1d) /100.
               end do
             end do

           case(32)                                ! 3: 2-cell divergent circular wind. Unit: m/s
             w = 0.
             if(fu_fails(mod(nx_dispersion,2)==0,'2-vortex wind needs even nx','diagnostic_wind_test'))return
             x_centre = real(nx_dispersion / 4) + 0.5 + mod(nx_dispersion/2, 2) / 2.
             y_centre = real(ny_dispersion / 2) + 0.5 + mod(ny_dispersion, 2) / 2.
!             speedScale = 4. / min(nx_dispersion, ny_dispersion)
             speedScale = pXCellSize(fs_dispersion/2) * 2. * Pi / T * 2.
             do iy = 1, ny_dispersion
               do ix = 1, nx_dispersion
                 i1d = ix + (iy-1) * nx_dispersion
                 if(ix <= nx_dispersion/2)then
                   u(i1d) = - real(iy-y_centre) * speedScale
                   v(i1d) = - real(x_centre+3-ix) * speedScale
                 else
                   u(i1d) = real(iy-y_centre) * speedScale
                   v(i1d) = real(x_centre-3+nx_dispersion/2 - ix) * speedScale
                 endif
               end do
             end do

           case(4)                              ! 5: double-vortex non-divergent rotation
             w = 0.                        ! Stream function: sin(x*2*Pi/nx)*cos(y*Pi/ny)
             do iy = 1, ny_dispersion
               do ix = 1, nx_dispersion
                 i1d = ix + (iy-1) * nx_dispersion
   ! 1x2 vortices
                 u(i1d) = 4.* sin(real(ix)*2.*Pi/nx_dispersion) * cos(real(iy)*Pi/ny_dispersion)
                 v(i1d) = -8.* cos(real(ix)*2.*Pi/nx_dispersion) * sin(real(iy)*Pi/ny_dispersion)
               end do
             end do
             
           case(5)                              ! 5: 8-vortex non-divergent rotation
             w = 0.                        ! Stream function: sin(x*2*Pi/nx)*cos(y*Pi/ny)
             do iy = 1, ny_dispersion
               do ix = 1, nx_dispersion
                 i1d = ix + (iy-1) * nx_dispersion
   ! 2x4 vortices
                 u(i1d) = 40.* sin(real(ix)*4.*Pi/nx_dispersion) * cos(real(iy)*2.*Pi/ny_dispersion)
                 v(i1d) = -80.* cos(real(ix)*4.*Pi/nx_dispersion) * sin(real(iy)*2.*Pi/ny_dispersion)
               end do
             end do
             
           case(6)                           ! 6: LauritzenGMD 2012 non-divergent Eq 18 and 19
                                             ! No rotation        
             w = 0.
             call  lonlat_grid_parameters(dispersion_grid,&
                                     & corner_lon_E, corner_lat_N,corner_in_geo_lonlat, &
                                     & number_of_points_x, number_of_points_y, &
                                     & southpole_lon_E, southpole_lat_N, &
                                     & dx_deg, dy_deg)
             do iy = 1, ny_dispersion
                lat = corner_lat_N + dy_deg *(iy-1)
                cosTheta = cos(lat * degrees_to_radians)
                sin2Theta = sin (2 * lat * degrees_to_radians)

               do ix = 1, nx_dispersion
                 lon = corner_lon_E + dx_deg *(ix-1)
                 lambda_prime = lon * degrees_to_radians !- 2 * pi_t_over_T 
                 sinlprime2 = sin(lambda_prime) ! sin^2
                 sinlprime2 = sinlprime2 * sinlprime2 
                 sin2lprime = sin(2*lambda_prime) 
                   
                 i1d = ix + (iy-1) * nx_dispersion

                 u(i1d) = 10 * sinlprime2 * sin2Theta * tfactor ! + 2* pi * earth_radius/T * cosTheta
                 v(i1d) = 10 * sin2lprime * cosTheta  * tfactor 
               end do
             end do

         case(7)                                   ! 7: Solid body rotation 
             w = 0.
             call  lonlat_grid_parameters(dispersion_grid,&
                                     & corner_lon_E, corner_lat_N,corner_in_geo_lonlat, &
                                     & number_of_points_x, number_of_points_y, &
                                     & southpole_lon_E, southpole_lat_N, &
                                     & dx_deg, dy_deg)
             do iy = 1, ny_dispersion
                lat = corner_lat_N + dy_deg *(iy-1)
                cosTheta = cos(lat * degrees_to_radians) 

               do ix = 1, nx_dispersion
                 i1d = ix + (iy-1) * nx_dispersion
                 u(i1d) = 2* pi * earth_radius/T * cosTheta
                 v(i1d) = 0.
               end do
             end do

           case(8)                           ! 6: LauritzenGMD 2012 non-divergent Eq 18 and 19
                                             ! With rotation        
             w = 0.
             call  lonlat_grid_parameters(dispersion_grid,&
                                     & corner_lon_E, corner_lat_N,corner_in_geo_lonlat, &
                                     & number_of_points_x, number_of_points_y, &
                                     & southpole_lon_E, southpole_lat_N, &
                                     & dx_deg, dy_deg)
             do iy = 1, ny_dispersion
                lat = corner_lat_N + dy_deg *(iy-1)
                cosTheta = cos(lat * degrees_to_radians)
                sin2Theta = sin (2 * lat * degrees_to_radians)

               do ix = 1, nx_dispersion
                 lon = corner_lon_E + dx_deg *(ix-1)
                 lambda_prime = lon * degrees_to_radians - 2 * pi_t_over_T 
                 sinlprime2 = sin(lambda_prime) ! sin^2
                 sinlprime2 = sinlprime2 * sinlprime2 
                 sin2lprime = sin(2*lambda_prime) 
                   
                 i1d = ix + (iy-1) * nx_dispersion

                 u(i1d) = 10 * sinlprime2 * sin2Theta * tfactor  + 2* pi * earth_radius/T * cosTheta
                 v(i1d) = 10 * sin2lprime * cosTheta  * tfactor 
               end do
             end do
!                                             ! With rotation        
!             w = 0.
!             call  lonlat_grid_parameters(dispersion_grid,&
!                                     & corner_lon_E, corner_lat_N,corner_in_geo_lonlat, &
!                                     & number_of_points_x, number_of_points_y, &
!                                     & southpole_lon_E, southpole_lat_N, &
!                                     & dx_deg, dy_deg)
!             lambda0 = corner_lon_E * degrees_to_radians
!             dlambda = dx_deg * degrees_to_radians
!             theta0 = corner_lat_N * degrees_to_radians
!             dtheta = dy_deg * degrees_to_radians
!             tfactor = cos(pi_t_over_T)
!             fTmp = 1.  / (T * dlambda) !Common factor
!             dy_m =  pYCellSize(1) ! All y sizes are same on sphere
!
!!             call msg("dy, dy", pYCellSize(1), dtheta*earth_radius)
!
!             do iy = 1, ny_dispersion
!               theta = theta0 + dtheta *(iy-1) 
!               sinTheta = sin(theta) ! divided by cos(theta) to account for
!               cosTheta = cos(theta)
!               sin2Theta = 2 * sinTheta * cosTheta
!               dx_m = pXCellSize(1 + (iy-1) * nx_dispersion) !Same at same latitude
!
!               !               call msg("dx, dx", pXCellSize(1 + (iy-1) * nx_dispersion), cos(theta)* dlambda *earth_radius)
!
!               do ix = 1, nx_dispersion
!                  lambda = lambda0 +  dlambda*(ix-1) 
!                  lambda_prime = lambda  - 2 * pi_t_over_T 
!                  sinlprime2 = sin(lambda_prime) ! sin^2
!                  sinlprime2 = sinlprime2 * sinlprime2
!                  sin2lprime = sin(2*lambda_prime) 
!                   
!                  i1d = ix + (iy-1) * nx_dispersion
!
!                  u(i1d) = dx_m * fTmp *  (5. * sinlprime2 * 2* sinTheta * tfactor  + 2* pi)
!                  v(i1d) = dy_m * fTmp *   5. * sin2lprime *    cosTheta * tfactor
!               end do
!             end do

           case(9)                           ! 6: solid body rotation around horizontal (0,0) line
             w = 0.
             call  lonlat_grid_parameters(dispersion_grid,&
                                     & corner_lon_E, corner_lat_N,corner_in_geo_lonlat, &
                                     & number_of_points_x, number_of_points_y, &
                                     & southpole_lon_E, southpole_lat_N, &
                                     & dx_deg, dy_deg)
             do iy = 1, ny_dispersion
                lat = corner_lat_N + dy_deg *(iy-1)
                cosTheta = cos(lat * degrees_to_radians)
                sinTheta = sin (lat * degrees_to_radians)

               do ix = 1, nx_dispersion
                 lon = corner_lon_E + dx_deg *(ix-1)
                 cosLambda = cos(lon * degrees_to_radians)
                 sinLambda = sin (lon * degrees_to_radians)

                 i1d = ix + (iy-1) * nx_dispersion

                 u(i1d) = Pi * earth_radius/T * sinTheta * cosLambda
                 v(i1d) = -Pi * earth_radius/T * sinLambda
               end do
             end do
         end select

        idTmp = fu_set_field_id(met_src, &
                              & dispersion_u_flag, &
                              & analysis_time, &
                              & forecast_length, &
                              & dispersion_grid, &
                              & fu_level(dispersion_vertical, iLev, .false.))
        call dq_store_2d(dispMarketPtr, idTmp, u, multi_time_stack_flag, .false., &
                         iUpdateType = overwrite_field, storage_grid=dispersion_grid)

        idTmp = fu_set_field_id(met_src, &
                              & dispersion_v_flag, &
                              & analysis_time, &
                              & forecast_length, &
                              & dispersion_grid, &
                              & fu_level(dispersion_vertical, iLev, .false.))

        call dq_store_2d(dispMarketPtr, idTmp, v, multi_time_stack_flag, .false., &
                         iUpdateType = overwrite_field, storage_grid=dispersion_grid)

        idTmp = fu_set_field_id(met_src, &
                              & dispersion_w_flag, &
                              & analysis_time, &
                              & forecast_length, &
                              & dispersion_grid, &
                              & fu_level(dispersion_vertical, iLev, .false.))

        call dq_store_2d(dispMarketPtr, idTmp, w, multi_time_stack_flag, .false., &
                         iUpdateType = overwrite_field, storage_grid=dispersion_grid)

      end do ! level
    end do ! time

    call free_work_array(u)
    call free_work_array(v)
    call free_work_array(w)

    call arrange_supermarket(dispMarketPtr)
    
  end subroutine diagnostic_wind_test

  
  !***************************************************************************************************

  subroutine diag_vertical_wind_incompr(met_src, metBufPtr, &
                                    & obstimes, &
                                    & ifHorizInterp, pHorizInterpStruct, &
                                    & ifVertInterp, pVertInterpStruct, &
                                    & dispMarketPtr, &
                                    & if_half_levels)
    !
    ! Determine the vertical wind in the midpoints of each grid cell
    ! from the incompressible continuity equation du/dx + dv/dy +
    ! dw/dz = 0, where z is height and w is the height-system vertical
    ! velocity.
    ! 
    ! Method: du/dx + dv/dy is computed in centres of grid cells and
    ! integrated from the bottom of the layer to the top. The value at
    ! the vertical midpoint ("full levels", w_full) is interpolated linearly. 
    ! 
    ! Vanishing w-wind is assumed on the bottom of first layer, while
    ! the wind at the top of last leyer can be nonvanishing.
    ! 
    ! This subroutine includes the possible interpolation from meteo
    ! to dispersion grid. The interpolated horizontal and the
    ! resulting vertical wind are stored into dispersion miniMarket.
    !
    ! The fields are stored into dispersion market.
    !
    implicit none


    ! Units: SI 
    ! Author: J. Vira

    type(meteo_data_source), intent(in) :: met_src
    type(Tfield_buffer), pointer :: metBufPtr
    type(silja_time), dimension(:), intent(in) :: obstimes
    
    ! meteo -> dispersion interpolation
    type(TVertInterpStruct), pointer :: pVertInterpStruct
    type(THorizInterpStruct), pointer :: pHorizInterpStruct
    ! See also module variables for interpolating winds.
    logical, intent(in) :: ifHorizInterp, ifVertInterp
    type(mini_market_of_stacks), intent(inout) :: dispMarketPtr
    logical, intent(in) :: if_half_levels

    ! Local variables
    
    ! Grids of the horizontal winds as in meteo.
    type(silja_grid) :: u_grid, v_grid
    type(silja_time) :: now

    type(silja_time), dimension(max_times) :: dispMarketTimes

    ! fields needed from the meteo stack
    real, dimension(:), pointer ::  u, v
   
    ! fields computed
    real, dimension(:), pointer :: dudx, dvdy, w, w_full 
    real :: dz       ! == h_{i+1/2} - h_{i-1/2}
    real :: rho_adv
    integer :: iTime, iLev, i, iu, iv, ip, it
    logical :: ifInterpolateUV
    ! Temperature and pressure for calculating time derivative of density. 
    real, dimension(:), pointer :: drho_dt, windfield, rho, drhodx, drhody
    real :: weight_past
    ! Check if dispersion stack must be arranged
    logical :: ifFirst = .true.
    ! Debug output
    logical :: ifWrite = .false.
    integer :: funit, iTmp, ix, iy, n_sm_times, n_dispMarketTimes
    
    ! Pressure and temperature fields to calculate d(rho) / dt.
    type(silja_field) :: p_2d, t_2d
    type(silja_field_id) :: idTmp !, meteoWindID
    type(silja_time) :: time, analysis_time
    type(silja_interval) :: forecast_length, dt

    rho => fu_work_array()
    drhodx => fu_work_array()
    drhody => fu_work_array()

    drho_dt => fu_work_array()
    dudx => fu_work_array()
    dvdy => fu_work_array()
    
    w => fu_work_array()

    iu = fu_index(metBufPtr%buffer_quantities, u_flag)
    iv = fu_index(metBufPtr%buffer_quantities, v_flag)
    ip = fu_index(metBufPtr%buffer_quantities, pressure_flag)
    it = fu_index(metBufPtr%buffer_quantities, temperature_flag)

    !
    ! Before starting anything, let's check whether wind has already been diagnosed
    !
    call supermarket_times(dispMarketPtr, met_src_missing, dispMarketTimes, n_dispMarketTimes)


    ! Interpolation of horizontal winds. The possibility of Arakawa
    ! shifts complicates matters, and currently it is only considered
    ! if the dispersion and meteo grids are otherwise equal. In this
    ! case, if ifHorizInterp is false but u_grid /= dispersion_grid,
    ! and separate interpolation structures are set for u and v.

    ! If horiz. interpolated is needed globally, always interpolate u, v.
    ifInterpolateUV = ifHorizInterp

    if (dispersion_grid == meteo_grid) then
      ! The u/v grids might still be shifted by a half cell.
      u_grid = fu_grid(metBufPtr%p4d(iu)%past%p2d(1)%idPtr)
      v_grid = fu_grid(metBufPtr%p4d(iv)%past%p2d(1)%idPtr)
!      u_grid = fu_u_grid(meteoWindID)
!      v_grid = fu_v_grid(meteoWindID)
    else
      !call set_error('Meteo to dispersion interpolation is broken', 'diag_vertical_wind_incompr')
      !return

      ! We have global interpolation, which have been defined before.
      u_grid = dispersion_grid
      v_grid = dispersion_grid
      pHorizInterpU => pHorizInterpStruct
      pHorizInterpV => pHorizInterpStruct

    end if
    
    if (.not. (u_grid == dispersion_grid) .or. .not. (v_grid == dispersion_grid)) then
      ! The case of shifted grids.

      ifInterpolateUV = .true.
      
      ! Only allocate once:
      if (.not. associated(pHorizInterpU)) then
        call msg('Allocating interpolation for shifted horizontal winds')
        allocate(pHorizInterpU, pHorizInterpV, stat=iTmp)

        if (iTmp /= 0) then
          call set_error('Cannot allocate', 'diag_vertical_wind_incompr')
          return
        end if

        pHorizInterpU => fu_horiz_interp_struct(u_grid, dispersion_grid, linear, .true.)
        pHorizInterpV => fu_horiz_interp_struct(v_grid, dispersion_grid, linear, .true.)
        
      end if
    end if
      
!!$    if (ifWrite) then
!!$      funit = fu_next_free_unit()
!!$      open(funit, file='diagn_wind', form='binary', iostat=iTmp)
!!$      if (iTmp /= 0) then
!!$        call msg('iostat = ', itmp)
!!$        call set_error('Failed to write debug output', 'diag_vertical_wind_incompr')
!!$        return
!!$      end if
!!$    end if
  
    !call msg('Making dispersion vertical w-wind')
    
    forecast_length = zero_interval
    !
    ! Start computations: loop over times. In principle, obstimes can
    ! contain any times, but since values are taken from meteo buffer,
    ! only times defined in buffer are actually considered. Normally
    ! these are the same, of course.
    ! 
    ! Attention. It is assumed that no other part of the program deals
    ! with the time-dependent part of dispersion buffer.

    do itime = 1, size(obstimes)
      now = obstimes(itime)
      if (.not. defined(now)) exit
      if (now == metBufPtr%time_past) then
        weight_past = 1.0
      else if (now == metBufPtr%time_future) then
        weight_past = 0.0
      else
        call set_error('Strange time:'+fu_str(obstimes(itime)),'diag_vertical_wind_incompr')
        return
      end if
      
!      if (n_DispMarketTimes > 0) then
!        !
!        ! If the market is not empty, check if this time has been processed already.
!        !
!        if (fu_closest_time(now, dispMarketTimes(1:n_dispMarketTimes), single_time, .true.) > 0) cycle
!      end if

      analysis_time = now 

      !
      ! Prepare all pointers
      !
      idTmp = fu_set_field_id(met_src, &
                            & dispersion_u_flag, &
                            & analysis_time, &
                            & forecast_length, &
                            & dispersion_grid, &
                            & fu_level(dispersion_vertical, 1, .false.))
      call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, u)       ! wind u
      if(fu_fails(.not.error,'Failed u field data pointer','diag_vertical_wind_incompr'))return

      call set_quantity(idTmp,dispersion_v_flag)
      call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, v)       ! wind v
      if(fu_fails(.not.error,'Failed v field data pointer','diag_vertical_wind_incompr'))return

      call set_quantity(idTmp, dispersion_w_flag)
      call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, w_full)  ! wind w
      if(fu_fails(.not.error,'Failed w_full field data pointer','diag_vertical_wid_incompr'))return

      call get_drho_dt(1)
      call get_winds(1, u, v, ifInterpolateUV)

      call ddx_of_field(u, dispersion_grid, dispersion_grid, dudx)
      call ddy_of_field(v, dispersion_grid, dispersion_grid, dvdy)

      if (ifWrite) write(funit) u(1:fs_dispersion), v(1:fs_dispersion)
      if (error) return
      !
      ! The first level. Assume that at ground level, w == 0.
      !
      call get_rho(1)
      call ddx_of_field(rho, dispersion_grid, dispersion_grid, drhodx)
      call ddy_of_field(rho, dispersion_grid, dispersion_grid, drhody)

      dz = fu_layer_thickness_m(fu_level(dispersion_vertical, 1))
!!$      do i = 1, fs_dispersion
!!$        rho_adv = (u(i)*drhodx(i) + v(i)*drhody(i)) / rho(i)
!!$        w(i) = -(dudx(i) + dvdy(i) + drho_dt(i)/rho(i) + rho_adv) * dz
!!$        w_full(i) = w(i) / 2.
!!$      end do

      if (if_half_levels) then
          do i = 1, fs_dispersion
            w_full(i) = -(dudx(i) + dvdy(i)) * dz
          end do
        else
          do i = 1, fs_dispersion
            rho_adv = (u(i)*drhodx(i) + v(i)*drhody(i)) / rho(i)
            w_full(i) = w(i) ! temporary
            !w(i) = w(i) - (dudx(i) + dvdy(i) + drho_dt(i)/rho(i) + rho_adv) * dz
            w(i) = w(i) - (dudx(i) + dvdy(i)) * dz
            w_full(i) = w(i) / 2.
          end do
        end if

      if (ifWrite) write(funit) w_full(1:fs_dispersion)
      if (error) return
      !
      ! Level 2 and up:
      !
      do iLev = 2, nz_dispersion

        call set_level(idTmp, fu_level(dispersion_vertical, iLev, .false.))
        call set_quantity(idTmp,dispersion_u_flag)
        call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, u)       ! wind u
        if(fu_fails(.not.error,'Failed u field data pointer','diag_vertical_wind_incompr'))return

        call set_quantity(idTmp,dispersion_v_flag)
        call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, v)       ! wind v
        if(fu_fails(.not.error,'Failed v field data pointer','diag_vertical_wind_incompr'))return

        call set_quantity(idTmp, dispersion_w_flag)
        call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, w_full)  ! wind w
        if(fu_fails(.not.error,'Failed w_full field data pointer','diag_vertical_wid_incompr'))return

        call get_winds(ilev, u, v, ifInterpolateUV)
        call get_drho_dt(iLev)

        if (ifWrite) write(funit) u(1:fs_dispersion), v(1:fs_dispersion)

        call ddx_of_field(u, dispersion_grid, dispersion_grid, dudx)
        call ddy_of_field(v, dispersion_grid, dispersion_grid, dvdy)

        call get_rho(iLev)
        call ddx_of_field(rho, dispersion_grid, dispersion_grid, drhodx)
        call ddy_of_field(rho, dispersion_grid, dispersion_grid, drhody)

        if (error) return

        dz = fu_layer_thickness_m(fu_level(dispersion_vertical, iLev))
        if (if_half_levels) then
          do i = 1, fs_dispersion
            w_full(i) = w_full(i) - (dudx(i) + dvdy(i)) * dz
          end do
        else
          do i = 1, fs_dispersion
            rho_adv = (u(i)*drhodx(i) + v(i)*drhody(i)) / rho(i)
            w_full(i) = w(i) ! temporary
            !w(i) = w(i) - (dudx(i) + dvdy(i) + drho_dt(i)/rho(i) + rho_adv) * dz
            w(i) = w(i) - (dudx(i) + dvdy(i)) * dz
            w_full(i) = (w_full(i) + w(i)) / 2.
          end do
        end if

        if (ifWrite) write(funit) w_full(1:fs_dispersion)

      end do ! level
    end do ! time

    call free_work_array(dudx)
    call free_work_array(dvdy)
    call free_work_array(w)
    call free_work_array(drho_dt)
    
    call free_work_array(rho)
    call free_work_array(drhodx)
    call free_work_array(drhody)

    
    call arrange_supermarket(dispMarketPtr)
   

    if (ifWrite) then
      close(funit)
      ifWrite = .false.
    end if
    
    !==========================================================================================
    
    contains 
      
      subroutine get_rho(ilev)
        implicit none
        integer :: ilev
        ! Local variables
        integer :: ix, iy
        real :: pr, t

        
        do iy = 1, ny_dispersion
          do ix = 1, nx_dispersion
            
            pr = fu_get_value(metBufPtr%p4d(ip), nx_meteo, ix, iy, iLev, weight_past, &
                 pHorizInterpStruct, pVertInterpStruct, ifHorizInterp, ifVertInterp)
            
            t = fu_get_value(metBufPtr%p4d(it), nx_meteo, ix, iy, iLev, weight_past, &
                 pHorizInterpStruct, pVertInterpStruct, ifHorizInterp, ifVertInterp)
            
            rho(ix + (iy-1)*nx_dispersion) = pr / (gas_constant_dryair*t)
            
          end do
        end do
        
      end subroutine get_rho

      subroutine get_drho_dt(ilev)
        ! Compute d rho / dt on level iLev
        implicit none
        integer :: iLev
        ! Local declarations
        integer :: ix, iy
        real :: p_now, p_next, p_prev, t_now, t_next, t_prev
        real :: dt_sec

        !dt_sec = fu_sec(fu_valid_time(metBufPtr%p4d(ip)%future%p2d(ilev)%idptr) &
        !     & - fu_valid_time(metBufPtr%p4d(ip)%past%p2d(ilev)%idptr))
        dt_sec = fu_sec(metBufPtr%time_future - metBufPtr%time_past)
        
        if (dt_sec .eps. 0.0) then
        !if (.true.) then
          call msg_warning('dt_sec is zero', 'get_drho_dt @ diag_vertical_wind_incompr')
          drho_dt(1:fs_dispersion) = 0.
          return
        end if

        if (ifHorizInterp .or. ifVertInterp) then
          do iy = 1, ny_dispersion
            do ix = 1, nx_dispersion

              p_now = fu_get_value(metBufPtr%p4d(ip), nx_meteo, ix, iy, iLev, weight_past, &
                   pHorizInterpStruct, pVertInterpStruct, ifHorizInterp, ifVertInterp)
              p_next = fu_get_value(metBufPtr%p4d(ip), nx_meteo, ix, iy, iLev, 0., &
                   pHorizInterpStruct, pVertInterpStruct, ifHorizInterp, ifVertInterp)
              p_prev = fu_get_value(metBufPtr%p4d(ip), nx_meteo, ix, iy, iLev, 1., &
                   pHorizInterpStruct, pVertInterpStruct, ifHorizInterp, ifVertInterp)

              t_now = fu_get_value(metBufPtr%p4d(it), nx_meteo, ix, iy, iLev, weight_past, &
                   pHorizInterpStruct, pVertInterpStruct, ifHorizInterp, ifVertInterp)
              t_next = fu_get_value(metBufPtr%p4d(it), nx_meteo, ix, iy, iLev, 0., &
                   pHorizInterpStruct, pVertInterpStruct, ifHorizInterp, ifVertInterp)
              t_prev = fu_get_value(metBufPtr%p4d(it), nx_meteo, ix, iy, iLev, 1., &
                   pHorizInterpStruct, pVertInterpStruct, ifHorizInterp, ifVertInterp)

              drho_dt(ix + (iy-1)*nx_dispersion) = (t_now * (p_next-p_prev)/dt_sec &
                   & + p_now * (t_next-t_prev)/dt_sec) / (gas_constant_dryair*t_now**2)

            end do
          end do

        else
          do iy = 1, ny_dispersion
            do ix = 1, nx_dispersion

              p_now = metBufPtr%p4d(ip)%past%p2d(iLev)%ptr(ix + (iy-1)*nx_meteo) * weight_past &
                   & + metBufPtr%p4d(ip)%future%p2d(iLev)%ptr(ix + (iy-1)*nx_meteo) * (1.-weight_past)
              p_prev = metBufPtr%p4d(ip)%past%p2d(iLev)%ptr(ix + (iy-1)*nx_meteo)
              p_next = metBufPtr%p4d(ip)%future%p2d(iLev)%ptr(ix + (iy-1)*nx_meteo)

              t_now = metBufPtr%p4d(it)%past%p2d(iLev)%ptr(ix + (iy-1)*nx_meteo) * weight_past &
                   & + metBufPtr%p4d(it)%future%p2d(iLev)%ptr(ix + (iy-1)*nx_meteo) * (1.-weight_past)
              t_prev = metBufPtr%p4d(it)%past%p2d(iLev)%ptr(ix + (iy-1)*nx_meteo)
              t_next = metBufPtr%p4d(it)%future%p2d(iLev)%ptr(ix + (iy-1)*nx_meteo)

              drho_dt(ix + (iy-1)*nx_dispersion) = &
                   (t_now * (p_next-p_prev)/dt_sec &
                  & + p_now * (t_next-t_prev)/dt_sec) / (gas_constant_dryair*t_now**2)

            end do
          end do

        end if

      end subroutine get_drho_dt

!    subroutine field_get_value(iFld, iLev, pHorizInterpStructNow, grid_data, ifInterpolate)
!      ! Interpolate from meteo buffer into vector grid_data (with
!      ! the size of dispersion grid) using a given interpolation
!      ! structure and the standard meteo -> dispersion vertical interpolation.
!
!      implicit none
!      integer, intent(in) :: iFld, iLev
!      type(THorizInterpStruct), pointer :: pHorizInterpStructNow
!      real, dimension(:), pointer :: grid_data
!      logical, intent(in) :: ifInterpolate
!
!      integer :: ix, iy
!
!      do iy = 1, ny_dispersion
!       do ix = 1, nx_dispersion
!          grid_data(ix + (iy-1)*nx_dispersion) = &
!               & fu_get_value(metBufPtr%p4d(iFld), nx_meteo, ix, iy, iLev, weight_past, &
!               &              pHorizInterpStructNow, pVertInterpStruct, ifInterpolate, ifVertInterp)
!        end do
!      end do
!
!    end subroutine field_get_value

      subroutine get_winds(iLev, data_u, data_v, ifInterpolate)
        ! 
        ! Get fields full of winds using wind_from_buffer.
        !
        implicit none
        integer, intent(in) :: iLev
        real, dimension(:), pointer :: data_u, data_v
        logical, intent(in) :: ifInterpolate

        integer :: ix, iy
        
        do iy = 1, ny_dispersion
          do ix = 1, nx_dispersion
            call wind_from_buffer_4d(metBufPtr%p4d(iu), metBufPtr%p4d(iv), &
                                   & nx_meteo, ix, iy, ilev, &
                                   & weight_past, &
                                   & pHorizInterpU, pVertInterpStruct, &
                                   & PHorizInterpV, pVertInterpStruct, &
                                   & ifInterpolate, ifVertInterp, &
                                   & data_u(ix + (iy-1)*nx_dispersion), &
                                   & data_v(ix + (iy-1)*nx_dispersion))
          end do
        end do

      end subroutine get_winds

      subroutine get_winds_l1(iLev, data_u, data_v, ifInterpolate)
        ! 
        ! Get fields full of winds using wind_from_buffer.
        !
        implicit none
        integer, intent(in) :: iLev
        real, dimension(:), pointer :: data_u, data_v
        logical, intent(in) :: ifInterpolate

        integer :: ix, iy
        
        call msg('Special l1 winds')

        do iy = 1, ny_dispersion
          do ix = 1, nx_dispersion
            call wind_from_buffer_4d(metBufPtr%p4d(iu), metBufPtr%p4d(iv), &
                                   & nx_meteo, ix, iy, 1, &
                                   & weight_past, &
                                   & pHorizInterpU, pVertInterpStruct, &
                                   & PHorizInterpV, pVertInterpStruct, &
                                   & ifInterpolate, .false., &
                                   & data_u(ix + (iy-1)*nx_dispersion), &
                                   & data_v(ix + (iy-1)*nx_dispersion))
          end do
        end do

      end subroutine get_winds_l1

  end subroutine diag_vertical_wind_incompr

  !***************************************************************************************************

  subroutine diag_vertical_wind_incompr_v2(met_src, metBufPtr, &
                                         & obstimes, &
                                         & ifHorizInterp, pHorizInterpStruct, &
                                         & ifVertInterp, pVertInterpStruct, &
                                         & dispMarketPtr, &
                                         & if_half_levels)
    !
    ! Determine the vertical wind in the midpoints of each grid cell
    ! from the incompressible continuity equation du/dx + dv/dy +
    ! dw/dz = 0, where z is height and w is the height-system vertical
    ! velocity.
    ! 
    ! Method: du/dx + dv/dy is computed in centres of grid cells and
    ! integrated from the bottom of the layer to the top. The value at
    ! the vertical midpoint ("full levels", w_full) is interpolated linearly. 
    ! 
    ! Vanishing w-wind is assumed on the bottom of first layer, while
    ! the wind at the top of last leyer can be nonvanishing.
    ! 
    ! This subroutine includes the possible interpolation from meteo
    ! to dispersion grid. The interpolated horizontal and the
    ! resulting vertical wind are stored into dispersion miniMarket.
    !
    ! The fields are stored into dispersion market.
    !
    implicit none


    ! Units: SI 
    ! Author: J. Vira

    type(meteo_data_source), intent(in) :: met_src
    type(Tfield_buffer), pointer :: metBufPtr
    type(silja_time), dimension(:), intent(in) :: obstimes
    
    ! meteo -> dispersion interpolation
    type(TVertInterpStruct), pointer :: pVertInterpStruct
    type(THorizInterpStruct), pointer :: pHorizInterpStruct
    ! See also module variables for interpolating winds.
    logical, intent(in) :: ifHorizInterp, ifVertInterp
    type(mini_market_of_stacks), intent(inout) :: dispMarketPtr
    logical, intent(in) :: if_half_levels

    ! Local variables
    
    ! Grids of the horizontal winds as in meteo.
    type(silja_grid) :: u_grid, v_grid
    type(silja_time) :: now

    type(silja_time), dimension(max_times) :: dispMarketTimes

    ! fields needed from the meteo stack
    real, dimension(:), pointer ::  u, v
   
    ! fields computed
    real, dimension(:), pointer :: dudx, dvdy, w, w_full, u_mean, v_mean 
    real :: dz       ! == h_{i+1/2} - h_{i-1/2}
    real :: rho_adv, dz_full
    integer :: iTime, iLev, i, iu, iv, ip, it, ismall, ismall_prev
    logical :: ifInterpolateUV
    ! Temperature and pressure for calculating time derivative of density. 
    real, dimension(:), pointer :: drho_dt, windfield, rho, drhodx, drhody
    real :: weight_past
    ! Check if dispersion stack must be arranged
    logical :: ifFirst = .true.
    ! Debug output
    logical :: ifWrite = .false.
    integer :: funit, iTmp, ix, iy, n_sm_times, n_dispMarketTimes
    integer :: ind_disp, ind_u, ind_v, nx_u, ny_v, fs_uv

    ! Pressure and temperature fields to calculate d(rho) / dt.
    type(silja_field) :: p_2d, t_2d
    type(silja_field_id) :: idTmp !, meteoWindID
    type(silja_time) :: time, analysis_time
    type(silja_interval) :: forecast_length, dt

    !rho => fu_work_array()
    !drhodx => fu_work_array()
    !drhody => fu_work_array()

    !drho_dt => fu_work_array()
    dudx => fu_work_array()
    dvdy => fu_work_array()
    
    w => fu_work_array()
    u => fu_work_array()
    v => fu_work_array()

    iu = fu_index(metBufPtr%buffer_quantities, u_flag)
    iv = fu_index(metBufPtr%buffer_quantities, v_flag)
    ip = fu_index(metBufPtr%buffer_quantities, pressure_flag)
    it = fu_index(metBufPtr%buffer_quantities, temperature_flag)

    call msg('Thin-layer vertical wind diagnosis')

    !
    ! Before starting anything, let's check whether wind has already been diagnosed
    !
    call supermarket_times(dispMarketPtr, met_src_missing, dispMarketTimes, n_dispMarketTimes)


    ! Interpolation of horizontal winds. The possibility of Arakawa
    ! shifts complicates matters, and currently it is only considered
    ! if the dispersion and meteo grids are otherwise equal. In this
    ! case, if ifHorizInterp is false but u_grid /= dispersion_grid,
    ! and separate interpolation structures are set for u and v.


    ! Interpolating the winds
    ! The u/v grids might still be shifted by a half cell.
    u_grid = fu_grid(metBufPtr%p4d(iu)%past%p2d(1)%idPtr)
    v_grid = fu_grid(metBufPtr%p4d(iv)%past%p2d(1)%idPtr)

    if (.not. (u_grid == meteo_grid) .or. .not. (v_grid == meteo_grid) &
      & .or. .not. (dispersion_grid == dispersion_grid) &
      & .or. .not. (dispersion_grid == dispersion_grid)) then
      ! Need separate u/v interpolation because either dispersion or meteo has different
      ! grids for the components.  
      pHorizInterpU => fu_horiz_interp_struct(u_grid, dispersion_grid, linear, .true.)
      pHorizInterpV => fu_horiz_interp_struct(v_grid, dispersion_grid, linear, .true.)
      ifInterpolateUV = .true.
    else
      ! The winds are not shifted from the main meteo grid, which may or 
      ! may not be the dispersion grid. If not, use interpolation as usual.
      pHorizInterpU => pHorizInterpStruct
      pHorizInterpV => pHorizInterpStruct
      ifInterpolateUV = ifHorizInterp
    end if

    if (ifFirst) call get_integration_vertical()
    
    forecast_length = zero_interval
    !
    ! Start computations: loop over times. In principle, obstimes can
    ! contain any times, but since values are taken from meteo buffer,
    ! only times defined in buffer are actually considered. Normally
    ! these are the same, of course.
    !
    do itime = 1, size(obstimes)
      now = obstimes(itime)
      if (.not. defined(now)) exit
      if (now == metBufPtr%time_past) then
        weight_past = 1.0
      else if (now == metBufPtr%time_future) then
        weight_past = 0.0
      else
        call set_error('Strange time:'+fu_str(obstimes(itime)),'diag_vertical_wind_incompr_v2')
        return
      end if
      
!      print *, "wind_interp_struct, Before:"
!      print *, "Ilev", wind_interp_struct%indLev(1:2,10,10,4)
!      print *, "Coef", wind_interp_struct%weight(1:2,10,10,4)
      call refine_interp_vert_coefs_v2(wind_interp_struct, metBufPtr, now)
!      print *, "wind_interp_struct, After:"
!      print *, "Ilev", wind_interp_struct%indLev(1:2,10,10,4)
!      print *, "Coef", wind_interp_struct%weight(1:2,10,10,4)
     
            
!      if (n_DispMarketTimes > 0) then
!        !
!        ! If the market is not empty, check if this time has been processed already.
!        !
!        if (fu_closest_time(now, dispMarketTimes(1:n_dispMarketTimes), single_time, .true.) > 0) cycle
!      end if

      analysis_time = now 
      ismall_prev = 1
      
      if (.not. if_half_levels) then
        call set_error('very unhappy: not if_half_levels', 'diag_vertical_wind_incompr_v2')
        return
      end if

      idTmp = fu_set_field_id(met_src, &
                            & dispersion_u_flag, &
                            & analysis_time, &
                            & forecast_length, &
                            & dispersion_grid, &
                            & fu_level(dispersion_vertical, 1, .false.))
      call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, u_mean)       ! wind u
      if(fu_fails(.not.error,'Failed u field data pointer','diag_vertical_wind_incompr_v2'))return

      call set_quantity(idTmp,dispersion_v_flag)
      call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, v_mean)       ! wind v
      if(fu_fails(.not.error,'Failed v field data pointer','diag_vertical_wind_incompr_v2'))return

      call set_quantity(idTmp, dispersion_w_flag)
      call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, w_full)  ! wind w
      if(fu_fails(.not.error,'Failed w_full field data pointer','diag_vertical_wind_incompr_v2'))return

      w_full(1:fs_dispersion) = 0.0
      
      fs_uv = nx_dispersion * ny_dispersion
      u_mean(1:fs_uv) = 0.0
      v_mean(1:fs_uv) = 0.0
      dz_full = 0.0

      nx_u = nx_dispersion
      ny_v = ny_dispersion
      
      do ismall = ismall_prev, ismall_prev + nsmall_in_layer(1) - 1
        dz = fu_layer_thickness_m(fu_level(integration_vertical, ismall))
        call get_winds(ismall, u, v, ifInterpolateUV)
        
        call ddx_of_field(u, dispersion_grid, dispersion_grid, dudx)
        call ddy_of_field(v, dispersion_grid, dispersion_grid, dvdy)
        u_mean(1:fs_uv) = u_mean(1:fs_uv) + u(1:fs_uv)*dz
        v_mean(1:fs_uv) = v_mean(1:fs_uv) + v(1:fs_uv)*dz
        
        do iy = 1, ny_dispersion
          do ix = 1, nx_dispersion
            ind_disp = ix + (iy-1)*nx_dispersion
            ind_u = wing_depth_w + ix + (iy-1)*(nx_u)
            ind_v = wing_depth_s*nx_dispersion + ix + (iy-1)*nx_dispersion
            w_full(ind_disp) = w_full(ind_disp) - (dudx(ind_u) + dvdy(ind_v))*dz
          end do
        end do
        dz_full = dz_full + dz
      end do

      u_mean(1:fs_uv) = u_mean(1:fs_uv) / dz_full
      v_mean(1:fs_uv) = v_mean(1:fs_uv) / dz_full
      
      ismall_prev = ismall
      if (error) return

      !call msg('Means of v_mean inside dispersion area:')
!!$      call msg('Means of v:')
!!$      do iy = 1, ny_dispersion_mpi
!!$        !windfield => v_mean(  1 +(iy-1)*nx_dispersion &
!!$        !                  & : (iy)*nx_dispersion)
!!$        !call msg('iy, mean', iy, sum(windfield)/nx_dispersion)
!!$        call msg('iy, mean', iy, v(1 + (iy-1)*nx_dispersion))
!!$      end do
      
!!$      if (ifFirst) then
!!$        open(999, access='stream', form='unformatted', file='vfld_' // fu_str(global_rank))
!!$        write(999) v(1:fs_v)
!!$        close(999)
!!$      end if
      

      if (ifWrite) write(funit) w_full(1:fs_dispersion)
      if (error) return
      !
      ! Level 2 and up:
      !
      do iLev = 2, nz_dispersion

        call set_level(idTmp, fu_level(dispersion_vertical, iLev, .false.))
        call set_quantity(idTmp,dispersion_u_flag)
        call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, u_mean)       ! wind u
        if(fu_fails(.not.error,'Failed u field data pointer','diag_vertical_wind_incompr_v2'))return

        call set_quantity(idTmp,dispersion_v_flag)
        call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, v_mean)       ! wind v
        if(fu_fails(.not.error,'Failed v field data pointer','diag_vertical_wind_incompr_v2'))return

        call set_quantity(idTmp, dispersion_w_flag)
        call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, w_full)  ! wind w
        if(fu_fails(.not.error,'Failed w_full field data pointer','diag_vertical_wind_incompr_v2'))return

        u_mean(1:fs_uv) = 0.0
        v_mean(1:fs_uv) = 0.0
        dz_full = 0.0
        
        do ismall = ismall_prev, ismall_prev + nsmall_in_layer(ilev) - 1
          dz = fu_layer_thickness_m(fu_level(integration_vertical, ismall))
          call get_winds(ismall, u, v, ifInterpolateUV)
          call ddx_of_field(u, dispersion_grid, dispersion_grid, dudx)
          call ddy_of_field(v, dispersion_grid, dispersion_grid, dvdy)
          u_mean(1:fs_uv) = u_mean(1:fs_uv) + u(1:fs_uv)*dz
          v_mean(1:fs_uv) = v_mean(1:fs_uv) + v(1:fs_uv)*dz
          
          do iy = 1, ny_dispersion
            do ix = 1, nx_dispersion
              ind_disp = ix + (iy-1)*nx_dispersion
              ind_u = wing_depth_w + ix + (iy-1)*(nx_u)
              ind_v = wing_depth_s*nx_dispersion + ix + (iy-1)*nx_dispersion
              w_full(ind_disp) = w_full(ind_disp) - (dudx(ind_u) + dvdy(ind_v))*dz
            end do
          end do
          dz_full = dz_full + dz
        end do
        
        u_mean(1:fs_uv) = u_mean(1:fs_uv) / dz_full
        v_mean(1:fs_uv) = v_mean(1:fs_uv) / dz_full

        ismall_prev = ismall
        
      end do ! level
    end do ! time

    call free_work_array(u)
    call free_work_array(v)
    call free_work_array(dudx)
    call free_work_array(dvdy)
    call free_work_array(w)
    
    call arrange_supermarket(dispMarketPtr)
   
    ifFirst = .false.

    if (ifWrite) then
      close(funit)
      ifWrite = .false.
    end if
    
    !==========================================================================================
    
  contains 
    

    subroutine get_winds(iLev, data_u, data_v, ifInterpolate)
      ! 
      ! Pick winds from buffer.  Will go thru grid_uv but store only grid_u for u and
      ! grid_v for v. The vertical and horizontal interpolations must be defined with
      ! grid_uv as target.
      !
      implicit none
      integer, intent(in) :: iLev
      real, dimension(:), pointer :: data_u, data_v
      logical, intent(in) :: ifInterpolate
      
      integer :: ix, iy, ind_u, ind_v
      real :: u, v
      
      !if (.not. ifInterpolate) then
      !  if (fu_fails(nx_dispersion_mpi == nx_dispersion, 'Need interpolation for this (x)', 'get_winds')) return
      !  if (fu_fails(ny_dispersion_mpi == ny_dispersion, 'Need interpolation for this (y)', 'get_winds')) return
      !end if
      do iy = 1, ny_dispersion
        do ix = 1, nx_dispersion
          call wind_from_buffer_4d(metBufPtr%p4d(iu), metBufPtr%p4d(iv), &
                                 & nx_meteo, ix, iy, ilev, &
                                 & weight_past, &
                                 & pHorizInterpU, wind_interp_struct, &
                                 & PHorizInterpV, wind_interp_struct, &
                                 & ifInterpolate, .true., &
                                 & u, v)
          
          ! in u_grid, iy' = iy-wing_depth_s, ix' = ix, i1d = ix' + (iy'-1)*nx_dispersion_mpi
          ind_u = ix + (iy-wing_depth_s-1)*nx_dispersion
          ! in v_grid, iy' = iy, ix' = ix-wing_depth_w, i1d = ix' + (iy'-1)*nx_dispersion
          ind_v = ix-wing_depth_w + (iy-1)*nx_dispersion

          ! fill u for 1 <= ix <= nx_dispersion_mpi, shift_ympi < iy <= ny_dispersion + shift_ympi
          ! fill v for 1 <= iy <= ny_dispersion_mpi, shift_xmpi < ix <= nx_dispersion + shift_xmpi
          if (wing_depth_s < iy .and. iy <= ny_dispersion + wing_depth_s) then
            data_u(ind_u) = u
          end if
          if (wing_depth_w < ix .and. ix <= nx_dispersion + wing_depth_w) then
            data_v(ind_v) = v
          end if

!!$          call wind_from_buffer_4d(metBufPtr%p4d(iu), metBufPtr%p4d(iv), &
!!$                                 & nx_meteo, ix, iy, ilev, &
!!$                                 & weight_past, &
!!$                                 & pHorizInterpU, wind_interp_struct, &
!!$                                 & PHorizInterpV, wind_interp_struct, &
!!$                                 & ifInterpolate, .true., &
!!$                                 & data_u(ix + (iy-1)*nx_dispersion), &
!!$                                 & data_v(ix + (iy-1)*nx_dispersion))
        end do
      end do
      
    end subroutine get_winds

    subroutine get_integration_vertical()
      implicit none
      
      integer :: nsublevs, ilev, lev_ind, nsmall, ilev_new, status
      type(silja_level) :: metlev, displev, new_level, top_small, bottom_big, bottom_small
      integer, dimension(:), pointer :: metlev_counter
      real :: lev_ind_float, thickness, height, bottom_height
      real, dimension(:), pointer :: met_data_col
      real :: met_data_srf
      
      metlev_counter => fu_work_int_array()
      metlev_counter(1:nz_dispersion) = 0

      met_data_col => fu_work_array()
      call vert_interp_data_crude(meteo_vertical, dispersion_vertical, met_data_col, met_data_srf)
      if (error) return

      do ilev = 1, nz_meteo
        metlev = fu_level(meteo_vertical, ilev, .true.)
        lev_ind_float = fu_project_level(metlev, dispersion_vertical, met_data_col, met_data_srf)
        lev_ind = int(lev_ind_float + 0.5)
        if (lev_ind > nz_dispersion) cycle
        metlev_counter(lev_ind) = metlev_counter(lev_ind) + 1
      end do
      call free_work_array(met_data_col)
      
      allocate(nsmall_in_layer(nz_dispersion), stat=status)
      if (status /= 0) then
        call set_error('Allocate failed', 'get_winds')
        return
      end if
      
      height = 0.0
      
      do ilev = 1, nz_dispersion
        nsmall = max(2, metlev_counter(ilev))
        nsmall_in_layer(ilev) = nsmall
        thickness = fu_layer_thickness_m(fu_level(dispersion_vertical, ilev))
        bottom_big = fu_lower_boundary_of_layer(fu_level(dispersion_vertical, ilev))
        bottom_height = fu_bottom_of_layer_value(dispersion_vertical, ilev)
        bottom_small = bottom_big
        do ilev_new = 1, nsmall
          height = height + thickness/nsmall
          top_small = fu_set_constant_height_level(height)
          new_level = fu_set_layer_between_two(bottom_small, top_small)
          bottom_small = top_small
          if (ilev == 1 .and. ilev_new == 1) then
            call set_vertical(new_level, integration_vertical)
          else
            call add_level(integration_vertical, new_level)
          end if
        end do
      end do

      call arrange_levels_in_vertical(integration_vertical)

      wind_interp_struct => fu_vertical_interp_struct(meteo_vertical, &
                                                    & integration_vertical, &
                                                    & dispersion_grid, &
                                                    & linear, &
                                                    & very_long_interval, 'wind_interp_incompr_2')

      call msg('Vertical for integrating continuity equation')

      call report(integration_vertical, ifFull=.true.)

      ! set inter
      ! refine
      call free_work_array(metlev_counter)

    end subroutine get_integration_vertical

  end subroutine diag_vertical_wind_incompr_v2
    
  !***************************************************************************************************
  
  subroutine diag_vertical_wind_anelastic(met_src, metBufPtr, &
                                        & obstimes, &
                                        & ifHorizInterp, pHorizInterpStruct, &
                                        & ifVertInterp, pVertInterpStruct, &
                                        & dispMarketPtr, &
                                        & if_half_levels)
    !
    ! Determine the vertical wind in the midpoints of each grid cell
    ! from the incompressible continuity equation du/dx + dv/dy +
    ! dw/dz = 0, where z is height and w is the height-system vertical
    ! velocity.
    ! 
    ! Method: du/dx + dv/dy is computed in centres of grid cells and
    ! integrated from the bottom of the layer to the top. The value at
    ! the vertical midpoint ("full levels", w_full) is interpolated linearly. 
    ! 
    ! Vanishing w-wind is assumed on the bottom of first layer, while
    ! the wind at the top of last leyer can be nonvanishing.
    ! 
    ! This subroutine includes the possible interpolation from meteo
    ! to dispersion grid. The interpolated horizontal and the
    ! resulting vertical wind are stored into dispersion miniMarket.
    !
    ! The fields are stored into dispersion market.
    !
    implicit none
    

    ! Units: SI 
    ! Author: TBA

    type(meteo_data_source), intent(in) :: met_src
    type(Tfield_buffer), pointer :: metBufPtr
    type(silja_time), dimension(:), intent(in) :: obstimes
    
    ! meteo -> dispersion interpolation
    type(TVertInterpStruct), pointer :: pVertInterpStruct
    type(THorizInterpStruct), pointer :: pHorizInterpStruct
    ! See also module variables for interpolating winds.
    logical, intent(in) :: ifHorizInterp, ifVertInterp
    type(mini_market_of_stacks), intent(inout) :: dispMarketPtr
    logical, intent(in) :: if_half_levels

    ! Local variables
    
    ! Grids of the horizontal winds as in meteo.
    type(silja_grid) :: u_grid, v_grid
    type(silja_time) :: now

    type(silja_time), dimension(max_times) :: dispMarketTimes

    ! fields needed from the meteo stack
    real, dimension(:), pointer ::  u, v
   
    ! fields computed
    real, dimension(:), pointer :: w_full, vertFluxb, vertFluxt
    real :: dz       ! == h_{i+1/2} - h_{i-1/2}
    real :: rho_adv
    integer :: iTime, iLev, i, iu, iv, ip, it, ismall_prev, ismall
    logical :: ifInterpolateUV
    ! Temperature and pressure for calculating time derivative of density. 
    real, dimension(:), pointer :: drho_dt, windfield, rho, drhodx, drhody, rho_c, rho_t
    real :: weight_past
    ! Check if dispersion stack must be arranged
    logical :: ifFirst = .true.
    ! Debug output
    logical :: ifWrite = .false.
    integer :: funit, iTmp, ix, iy, n_sm_times, n_dispMarketTimes
    character (len=*), parameter :: sub_name = "diag_vertical_wind_anelastic"
    
    ! Pressure and temperature fields to calculate d(rho) / dt.
    type(silja_field) :: p_2d, t_2d
    type(silja_field_id) :: idTmp !, meteoWindID
    type(silja_time) :: time, analysis_time
    type(silja_interval) :: forecast_length, dt
    
    type(silja_level), dimension(:), pointer :: rhoLevels
    type(silam_vertical) :: vertRho
    real, dimension(:), pointer ::  xSizePtr, ySizePtr, zSizePtr, u_mean, v_mean, total_mass
    real :: uw, ue, vn, vs, rhoe, rhow, rhon, rhos, dxn, dxs, dyw, dye, area
    integer :: fs_uv
    integer :: ind_rho, ind_disp, ind_u, ind_v, ix_u, iy_v, nx_u, ny_v
    integer :: nx_dispersion_mpi, ny_dispersion_mpi

    nx_dispersion_mpi = nx_dispersion + wing_depth_e + wing_depth_w
    ny_dispersion_mpi = ny_dispersion + wing_depth_s + wing_depth_n    

    rho_t => fu_work_array()
    rho_c => fu_work_array()

    drho_dt => fu_work_array()

    vertFluxb => fu_work_array()
    vertFluxt => fu_work_array()

    u_mean => fu_work_array()
    v_mean => fu_work_array()

    total_mass => fu_work_array()

    iu = fu_index(metBufPtr%buffer_quantities, u_flag)
    iv = fu_index(metBufPtr%buffer_quantities, v_flag)
    ip = fu_index(metBufPtr%buffer_quantities, pressure_flag)
    it = fu_index(metBufPtr%buffer_quantities, temperature_flag)
    if (.not. all((/iu, iv, ip, it/) > 0 )) then
          call msg ("iu, iv, ip, it", (/iu, iv, ip, it/))
          call set_error("Meteo input is not available", sub_name) 
          return
    endif

    xSizePtr => fu_grid_data(dispersion_cell_x_size_fld)
    ySizePtr => fu_grid_data(dispersion_cell_y_size_fld)
    zSizePtr => fu_work_array()
    do iLev = 1, nz_dispersion
      zSizePtr(iLev) = fu_layer_thickness_m(fu_level(dispersion_vertical, iLev)) 
    enddo

    !
    ! Before starting anything, let's check whether wind has already been diagnosed
    !
    call supermarket_times(dispMarketPtr, met_src_missing, dispMarketTimes, n_dispMarketTimes)


    ! Interpolation of horizontal winds. The possibility of Arakawa
    ! shifts complicates matters, and currently it is only considered
    ! if the dispersion and meteo grids are otherwise equal. In this
    ! case, if ifHorizInterp is false but u_grid /= dispersion_grid,
    ! and separate interpolation structures are set for u and v.

    ! Interpolating the winds
    ! The u/v grids might still be shifted by a half cell.
    u_grid = fu_grid(metBufPtr%p4d(iu)%past%p2d(1)%idPtr)
    v_grid = fu_grid(metBufPtr%p4d(iv)%past%p2d(1)%idPtr)

    if (.not. (u_grid == meteo_grid) .or. .not. (v_grid == meteo_grid) ) then

      ! Need separate u/v interpolation because either dispersion or meteo has different
      ! grids for the components.  
      pHorizInterpU => fu_horiz_interp_struct(u_grid, dispersion_grid, linear, .true.)
      pHorizInterpV => fu_horiz_interp_struct(v_grid, dispersion_grid, linear, .true.)
      pHorizInterpRho => fu_horiz_interp_struct(meteo_grid, dispersion_grid, linear, .true.)
      ifInterpolateUV = .true.
    else
      ! The winds are not shifted from the main meteo grid, which may or 
      ! may not be the dispersion grid. If not, use interpolation as usual.
      pHorizInterpU => pHorizInterpStruct
      pHorizInterpV => pHorizInterpStruct
      pHorizInterpRho => pHorizInterpStruct
      ifInterpolateUV = ifHorizInterp
    end if
    ! 
    ! Now u, v and rho will be interpolated into whole dispersion-uv grid. The diagnostics
    ! will be done for this area. This means some extra work, but avoids massive
    ! complexity with the grid indexing inside the diagnostic part.
    
    ! And vertical interpolation structure for level top
    if(.not. associated(pVertInterpRho))then
      call msg('Allocating vertical interpolation for vertical winds')
      allocate(rhoLevels(nz_dispersion))
      do iLev = 1, nz_dispersion
        rhoLevels(iLev) = fu_upper_boundary_of_layer(fu_level(dispersion_vertical, iLev))
      enddo
      call set_vertical(rhoLevels, vertRho)
      pVertInterpRho => fu_vertical_interp_struct(meteo_vertical, vertRho, dispersion_grid, linear, &
                                                & one_hour, 'meteo_to_disp_leveltop')
      if(error)return
      call set_missing(vertRho, .false.)
      deallocate(rhoLevels)
    endif


    !call msg('Making dispersion vertical w-wind')
    
    forecast_length = zero_interval
    !
    ! Start computations: loop over times. In principle, obstimes can
    ! contain any times, but since values are taken from meteo buffer,
    ! only times defined in buffer are actually considered. Normally
    ! these are the same, of course.
    ! 
    
    if (ifFirst) call get_integration_vertical()
    ifFirst = .false.
    
    ! The grids used by this sub, with mpi:
    ! u, v -> interpolated into dispersion u, v grids
    ! rho -> interpolated into joint area of u, v but not corners
    fs_uv = nx_dispersion*ny_dispersion

    do itime = 1, size(obstimes)
      now = obstimes(itime)
      if (.not. defined(now)) exit
      if (now == metBufPtr%time_past) then
        weight_past = 1.0
      else if (now == metBufPtr%time_future) then
        weight_past = 0.0
      else
        call set_error('Strange meteo time:'+fu_str(obstimes(itime)), &
                     & 'diag_vertical_wind_anelastic')
      end if

      analysis_time = now 

      vertFluxt(1:fs_uv) = 0.
      
      call start_count('refine in wind')
      call refine_interp_vert_coefs_v2(wind_interp_struct, metBufPtr, now)
      call refine_interp_vert_coefs_v2(pVertInterpRho,  metBufPtr,  now)
      call stop_count('refine in wind')
      if(error)return

      ismall_prev = 1

      nx_u = nx_dispersion
      ny_v = ny_dispersion

      idTmp = fu_set_field_id(met_src, &                                            
                            & dispersion_u_flag, &                                  
                            & now, &                                      
                            & zero_interval, &                                    
                            & dispersion_grid, &                                    
                            & fu_level(dispersion_vertical, 1, .false.)) ! Some level

      do iLev = 1, nz_dispersion
        
        call set_level(idTmp, fu_level(dispersion_vertical, iLev, .false.))
        call set_quantity(idTmp,dispersion_u_flag)
        call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, u)       ! wind u
        if(fu_fails(.not.error,'Failed u field data pointer','diag_vertical_wind_anelastic'))return

        call set_quantity(idTmp,dispersion_v_flag)
        call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, v)       ! wind v
        if(fu_fails(.not.error,'Failed v field data pointer','diag_vertical_wind_anelastic'))return

        call set_quantity(idTmp, dispersion_w_flag)
        call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, w_full)  ! wind w
        if(fu_fails(.not.error,'Failed w_full field data pointer','diag_vertical_wind_anelastic'))return

        call get_rho(ilev, pVertInterpRho, rho_t)
        u_mean(1:fs_uv) = 0.0
        v_mean(1:fs_uv) = 0.0
        total_mass(1:fs_uv) = 0.0

        do ismall = ismall_prev, ismall_prev + nsmall_in_layer(ilev) - 1
          !$OMP PARALLEL DEFAULT(NONE) &
          !$OMP & SHARED(u, v, integration_vertical, ismall, ifInterpolateUV, wind_interp_struct, rho_c, &
          !$OMP        & error, vertFluxb, vertFluxt, nx_dispersion, ny_dispersion, xSizePtr, ySizePtr, &
          !$OMP        & u_mean, v_mean, total_mass, fs_dispersion, drho_dt, fs_uv, &
          !$OMP        & nx_dispersion_mpi, ny_dispersion_mpi, disp_grid_size_x, disp_grid_size_y) &
          !$OMP & PRIVATE(dz, ix, iy, i, area, uw, ue, rhow, rhoe, dye, dyw, dxs, dxn, rhon, rhos, vs, vn)
          
          dz = fu_layer_thickness_m(fu_level(integration_vertical, ismall))
          ! In this sub, the diagnostic will be done on full dispersion_grid.  Some
          ! work could be avoided by restricting into dispersion_grid, but u, v are still
          ! needed in wider area. Dealing with 1D indices of 4 different grids just gets
          ! too technical.
          call get_winds(ismall, u, v, ifInterpolateUV)
          call get_drho_dt(ismall, wind_interp_struct)
          call get_rho(ismall, wind_interp_struct, rho_c)
          
          !rho_c = 1.
          !rho_t = 1.0
          !drho_dt = 0.0
          !if (error) return

          !$OMP SINGLE
          vertFluxb(1:fs_dispersion) = vertFluxt(1:fs_dispersion)
          !$OMP END SINGLE
          !$OMP BARRIER
          i = 0
          !$OMP DO
          do iy = 1, ny_dispersion
            do ix = 1, nx_dispersion
              i = ix + (iy-1) * nx_dispersion

              if(ix == 1)then 
                !on borders just take the values in the middle of the gridcell and cell size as half
                uw   = u(i)
                rhow = rho_c(i)
                dyw  = disp_grid_size_y(ix,iy) ! ySizePtr(i)
              else
                uw   = 0.5 * (u(i-1)+u(i))
                rhow = 0.5 * (rho_c(i-1)+rho_c(i))
                dyw  = 0.5 * (disp_grid_size_y(ix-1,iy) + disp_grid_size_y(ix,iy))!(ySizePtr(i-1) + ySizePtr(i))
              endif
              if(ix == nx_dispersion_mpi)then
                ue   = u(i)
                rhoe = rho_c(i)
                dye  = disp_grid_size_y(ix,iy)!ySizePtr(i)
              else
                ue   = 0.5 * (u(i+1)+u(i))
                rhoe = 0.5 * (rho_c(i+1)+rho_c(i))
                dye  = 0.5 * (disp_grid_size_y(ix+1,iy) + disp_grid_size_y(ix,iy))!(ySizePtr(i+1)+ ySizePtr(i))
              endif
              if(iy == 1)then
                vs   = v(i)
                rhos = rho_c(i)
                dxs  = disp_grid_size_x(ix,iy)!xSizePtr(i)
              else
                vs   = 0.5 * (v(i-nx_dispersion) + v(i))
                rhos = 0.5 * (rho_c(i-nx_dispersion) + rho_c(i))
                dxs  = 0.5 * (disp_grid_size_x(ix,iy-1) + disp_grid_size_x(ix,iy))
                !(xSizePtr(i-nx_dispersion) + xSizePtr(i))
              endif
              if(iy == ny_dispersion)then
                vn   = v(i)
                rhon = rho_c(i)
                dxn  = disp_grid_size_x(ix,iy) !xSizePtr(i)
              else
                vn   = 0.5 * (v(i+nx_dispersion) + v(i)) 
                rhon = 0.5 * (rho_c(i+nx_dispersion) + rho_c(i))
                dxn  = 0.5 * (disp_grid_size_x(ix,iy+1) + disp_grid_size_x(ix,iy))
                !(xSizePtr(i+nx_dispersion) + xSizePtr(i))
              endif

              if(iy == ny_dispersion .or. iy == 1)then
                dyw =  0.5 * dyw
                dye =  0.5 * dye
              endif
              if(ix == nx_dispersion .or. ix == 1)then
                dxs =  0.5 * dxs
                dxn =  0.5 * dxn
              endif
              area = 0.5*(dyw+dye) * 0.5*(dxs+dxn)

              ! vertFluxt(i)= ((uw*rhow -ue*rhoe)*ySizePtr(i) + (vs*rhos - vn*rhon)*xSizePtr(i))*zSizePtr(iLev) + &
              ! & vertFluxb(i) - drho_dt(i)*xSizePtr(i)*ySizePtr(i)*zSizePtr(iLev)

              vertFluxt(i) = vertFluxb(i)*area + &
                   & (uw*rhow*dyw - ue*rhoe*dye + vs*rhos*dxs - vn*rhon*dxn) * dz - &
                   &  drho_dt(i) * area * dz
              ! Convert to flux density
              vertFluxt(i) = vertFluxt(i) / area

              u_mean(i) = u_mean(i) + rho_c(i)*dz*u(i)
              v_mean(i) = v_mean(i) + rho_c(i)*dz*v(i)
              total_mass(i) = total_mass(i) + rho_c(i) * dz
            enddo
          enddo
          !$OMP END DO
          !$OMP END PARALLEL
        end do ! ismall
        
       

        ismall_prev = ismall

        do iy = 1, ny_dispersion
          iy_v = iy + wing_depth_s
          do ix = 1, nx_dispersion
            ! Computing 1d indices for u, v and uv grids:
            ! u -> ix_u, iy, nx_u, (ny)
            ! v -> ix,   iy_v, nx, (ny_v)
            ! uv -> ix_u, iy_v, nx_u, (ny_v)
            ! where nx, ny are for dispersion grid
            ix_u = ix + wing_depth_w

            ind_disp = ix + (iy-1) * nx_dispersion
            ind_u = ix_u + (iy-1)*nx_u
            ind_v =  ix + (iy_v-1)*nx_dispersion
            ind_rho = ix_u + (iy_v-1)*nx_u
            
            if (if_half_levels) then
              w_full(ind_disp) = vertFluxt(ind_rho) / rho_t(ind_rho)
            else
              w_full(ind_disp) = 0.5 * (vertFluxb(ind_rho) + vertFluxt(ind_rho)) &
                   & / rho_t(ind_rho)
            end if
          end do
          ! mean u wind in u grid
          do ix_u = 1, nx_u
            ind_rho = ix_u + (iy_v-1)*nx_u
            ind_u = ix_u + (iy-1)*nx_u
            u(ind_u) = u_mean(ind_rho) / total_mass(ind_rho)
          end do
        end do
        ! mean v in v grid
        do iy_v = 1, ny_v
          do ix = 1, nx_dispersion
            ix_u = ix + wing_depth_w
            ind_rho = ix_u + (iy_v-1)*nx_u
            ind_v =  ix + (iy_v-1)*nx_dispersion
            v(ind_v) = v_mean(ind_rho) / total_mass(ind_rho)
          end do
        end do
        
        
!!$        if(if_half_levels)then ! Now actually at the top of the layer
!!$          w_full(1:fs_dispersion) = vertFluxt(1:fs_dispersion) / (rho_t(1:fs_dispersion))
!!$        else
!!$          w_full(i:fs_dispersion) = 0.5*(vertFluxb(1:fs_dispersion) &
!!$                                       & + vertFluxt(1:fs_dispersion)) &
!!$                                     & /(rho_t(1:fs_dispersion))
!!$        endif
!!$
!!$        u(1:fs_dispersion) = u_mean(1:fs_dispersion) / total_mass(1:fs_dispersion)
!!$        v(1:fs_dispersion) = v_mean(1:fs_dispersion) / total_mass(1:fs_dispersion)

        !call msg('Setting test winds...')
        !u(1:fs_dispersion) = 0.0
        !v(1:fs_dispersion) = 0.0
        !w_full(1:fs_dispersion) = 0.0

      end do ! level

    end do ! time

    call free_work_array(drho_dt)
    
    call free_work_array(rho_t)
    call free_work_array(rho_c)
    call free_work_array(vertFluxb)
    call free_work_array(vertFluxt)

    call free_work_array(zSizePtr)
    
    call free_work_array(u_mean)
    call free_work_array(v_mean)
    call free_work_array(total_mass)

    call arrange_supermarket(dispMarketPtr)
   

    if (ifWrite) then
      close(funit)
      ifWrite = .false.
    end if
    
    !==========================================================================================
    
  contains 

    subroutine get_rho(ilev, VertInterp, data_rho)
      ! interpolate rho to rho grid == uv_grid
      implicit none
      integer :: ilev
      type(TVertInterpStruct), pointer  ::  VertInterp
      real, dimension(:), pointer :: data_rho
      ! Local variables
      integer :: ix, iy
      real :: pr, t
      
      !$OMP DO
      do iy = 1, ny_dispersion
        do ix = 1, nx_dispersion
          
          pr = fu_get_value(metBufPtr%p4d(ip), nx_meteo, ix, iy, iLev, weight_past, &
               pHorizInterpRho, VertInterp, ifHorizInterp, ifVertInterp)
          
          t = fu_get_value(metBufPtr%p4d(it), nx_meteo, ix, iy, iLev, weight_past, &
               pHorizInterpRho, VertInterp, ifHorizInterp, ifVertInterp)
!          print *, smpi_adv_rank, ': ', pr, gas_constant_dryair, t 
          data_rho(ix + (iy-1)*nx_dispersion_mpi) = pr / (gas_constant_dryair*t)
          !data_rho(ix + (iy-1)*nx_dispersion) = 1.0
          
        end do
      end do
      !$OMP END DO
    end subroutine get_rho

    subroutine get_drho_dt(ilev, vertInterp)
      ! Compute d rho / dt on level iLev
      implicit none
      integer :: iLev
      type(TVertInterpStruct), pointer  ::  VertInterp

      ! Local declarations
      integer :: ix, iy
      real :: p_now, p_next, p_prev, t_now, t_next, t_prev
      real :: dt_sec

      logical, parameter :: no_drho_dt = .false.

      !dt_sec = fu_sec(fu_valid_time(metBufPtr%p4d(ip)%future%p2d(ilev)%idptr) &
      !     & - fu_valid_time(metBufPtr%p4d(ip)%past%p2d(ilev)%idptr))
      dt_sec = fu_sec(metBufPtr%time_future - metBufPtr%time_past)
      
      if ((dt_sec .eps. 0.0) .or. no_drho_dt) then
        !if (.true.) then
        !call msg_warning('dt_sec is zero', 'get_drho_dt @ diag_vertical_wind_anelastic')
        !$OMP SINGLE
        drho_dt(1:fs_dispersion) = 0.
        !$OMP END SINGLE
        return
      end if
      
      if (ifHorizInterp .or. ifVertInterp) then
      !$OMP DO
        do iy = 1, ny_dispersion_mpi
          do ix = 1, nx_dispersion_mpi

            p_now = fu_get_value(metBufPtr%p4d(ip), nx_meteo, ix, iy, iLev, weight_past, &
                 pHorizInterpRho, vertInterp, ifHorizInterp, ifVertInterp)
            p_next = fu_get_value(metBufPtr%p4d(ip), nx_meteo, ix, iy, iLev, 0., &
                 pHorizInterpRho, vertInterp, ifHorizInterp, ifVertInterp)
            p_prev = fu_get_value(metBufPtr%p4d(ip), nx_meteo, ix, iy, iLev, 1., &
                 pHorizInterpRho, vertInterp, ifHorizInterp, ifVertInterp)

            t_now = fu_get_value(metBufPtr%p4d(it), nx_meteo, ix, iy, iLev, weight_past, &
                 pHorizInterpRho, vertInterp, ifHorizInterp, ifVertInterp)
            t_next = fu_get_value(metBufPtr%p4d(it), nx_meteo, ix, iy, iLev, 0., &
                 pHorizInterpRho, vertInterp, ifHorizInterp, ifVertInterp)
            t_prev = fu_get_value(metBufPtr%p4d(it), nx_meteo, ix, iy, iLev, 1., &
                 pHorizInterpRho, vertInterp, ifHorizInterp, ifVertInterp)

            drho_dt(ix + (iy-1)*nx_dispersion_mpi) = (t_now * (p_next-p_prev)/dt_sec &
                 & + p_now * (t_next-t_prev)/dt_sec) / (gas_constant_dryair*t_now**2)

          end do
        end do
        !$OMP END DO
      else
        if (fu_fails(nx_dispersion*ny_dispersion == nx_dispersion_mpi*ny_dispersion_mpi, &
                   & 'Should not be here', 'get_drho_dt')) return
        !$OMP DO
        do iy = 1, ny_dispersion
          do ix = 1, nx_dispersion

            p_now = metBufPtr%p4d(ip)%past%p2d(iLev)%ptr(ix + (iy-1)*nx_meteo) * weight_past &
                 & + metBufPtr%p4d(ip)%future%p2d(iLev)%ptr(ix + (iy-1)*nx_meteo) * (1.-weight_past)
            p_prev = metBufPtr%p4d(ip)%past%p2d(iLev)%ptr(ix + (iy-1)*nx_meteo)
            p_next = metBufPtr%p4d(ip)%future%p2d(iLev)%ptr(ix + (iy-1)*nx_meteo)

            t_now = metBufPtr%p4d(it)%past%p2d(iLev)%ptr(ix + (iy-1)*nx_meteo) * weight_past &
                 & + metBufPtr%p4d(it)%future%p2d(iLev)%ptr(ix + (iy-1)*nx_meteo) * (1.-weight_past)
            t_prev = metBufPtr%p4d(it)%past%p2d(iLev)%ptr(ix + (iy-1)*nx_meteo)
            t_next = metBufPtr%p4d(it)%future%p2d(iLev)%ptr(ix + (iy-1)*nx_meteo)

            drho_dt(ix + (iy-1)*nx_dispersion) = &
                 (t_now * (p_next-p_prev)/dt_sec &
                 & + p_now * (t_next-t_prev)/dt_sec) / (gas_constant_dryair*t_now**2)

          end do
        end do
        !$OMP END DO
      end if
      
    end subroutine get_drho_dt

    subroutine get_winds(iLev, data_u, data_v, ifInterpolate)
      ! 
      ! Get fields full of winds using wind_from_buffer.
      !
      implicit none
      integer, intent(in) :: iLev
      real, dimension(:), intent(inout) :: data_u, data_v
      logical, intent(in) :: ifInterpolate
      
      integer :: ix, iy
      
      !$OMP DO
      do iy = 1, ny_dispersion_mpi
        do ix = 1, nx_dispersion_mpi
          call wind_from_buffer_4d(metBufPtr%p4d(iu), metBufPtr%p4d(iv), &
                                 & nx_meteo, ix, iy, ilev, &
                                 & weight_past, &
                                 & pHorizInterpU, wind_interp_struct, &
                                 & PHorizInterpV, wind_interp_struct, &
                                 & ifInterpolate, .true., &
                                 & data_u(ix + (iy-1)*nx_dispersion_mpi), &
                                 & data_v(ix + (iy-1)*nx_dispersion_mpi))
        end do
      end do
      !$OMP END DO
    end subroutine get_winds

    subroutine get_integration_vertical()
      implicit none
      
      integer :: nsublevs, ilev, lev_ind, nsmall, ilev_new, status
      type(silja_level) :: metlev, displev, new_level, top_small, bottom_big, bottom_small
      integer, dimension(:), pointer :: metlev_counter
      real :: lev_ind_float, thickness, height, bottom_height
      real, dimension(:), pointer :: met_data_col
      real :: met_data_srf

      metlev_counter => fu_work_int_array()
      metlev_counter(1:nz_dispersion) = 0

      met_data_col => fu_work_array()
      call vert_interp_data_crude(meteo_vertical, dispersion_vertical, met_data_col, met_data_srf)
      if (error) return

      do ilev = 1, nz_meteo
        metlev = fu_level(meteo_vertical, ilev, .true.)
        !displev = fu_level_to_vertical_crude(metlev, dispersion_vertical)
        lev_ind_float = fu_project_level(metlev, dispersion_vertical, met_data_col, met_data_srf)
        !lev_ind_float = fu_level_index(displev, dispersion_vertical)
        lev_ind = int(lev_ind_float + 0.5)
        if (lev_ind > nz_dispersion) cycle
        metlev_counter(lev_ind) = metlev_counter(lev_ind) + 1
      end do
      call free_work_array(met_data_col)

      allocate(nsmall_in_layer(nz_dispersion), stat=status)
      if (status /= 0) then
        call set_error('Allocate failed', 'get_integration_vertical')
        return
      end if
      
      height = 0.0
      
      do ilev = 1, nz_dispersion
        nsmall = max(2, metlev_counter(ilev))
        !nsmall = 1
        nsmall_in_layer(ilev) = nsmall
        thickness = fu_layer_thickness_m(fu_level(dispersion_vertical, ilev))
        bottom_big = fu_lower_boundary_of_layer(fu_level(dispersion_vertical, ilev))
        bottom_height = fu_bottom_of_layer_value(dispersion_vertical, ilev)
        bottom_small = bottom_big
        do ilev_new = 1, nsmall
          height = height + thickness/nsmall
          top_small = fu_set_constant_height_level(height)
          new_level = fu_set_layer_between_two(bottom_small, top_small)
          bottom_small = top_small
          if (ilev == 1 .and. ilev_new == 1) then
            call set_vertical(new_level, integration_vertical)
          else
            call add_level(integration_vertical, new_level)
          end if
        end do
      end do

      call arrange_levels_in_vertical(integration_vertical)

      wind_interp_struct => fu_vertical_interp_struct(meteo_vertical, integration_vertical, &
                                                    & dispersion_grid, linear, &
                                                    & one_hour, 'wind_itnerp_anelastic')
      call msg('Vertical for integrating anelastic continuity equation')

      call report(integration_vertical, ifFull=.true.)

      ! set inter
      ! refine
      call free_work_array(metlev_counter)

    end subroutine get_integration_vertical

  end subroutine diag_vertical_wind_anelastic


  !************************************************************************************
  
  subroutine diag_vertical_wind_eta(met_src, metBufPtr, &
                                  & obstimes, &
                                  & ifHorizInterp, pHorizInterpStruct, &
                                  & ifVertInterp, pVertInterpStruct, &
                                  & dispMarketPtr, &
                                  & if_half_levels, if_from_omega, if_top_down)
    !
    ! Determine the vertical wind in the midpoints of each grid cell
    ! from the incompressible continuity equation du/dx + dv/dy +
    ! dw/dz = 0, where z is height and w is the eta-system vertical
    ! velocity.
    ! 
    ! Method: du/dx + dv/dy is computed in centres of grid cells and
    ! integrated from the bottom of the layer to the top. The value at
    ! the vertical midpoint ("full levels", w_full) is interpolated linearly. 
    ! 
    ! Vanishing w-wind is assumed on the bottom of first layer, while
    ! the wind at the top of last leyer can be nonvanishing.
    ! 
    ! This subroutine includes the possible interpolation from meteo
    ! to dispersion grid. The interpolated horizontal and the
    ! resulting vertical wind are stored into dispersion miniMarket.
    !
    ! The fields are stored into dispersion market.
    !
    ! Units: SI 
    ! Author: J. Vira
    implicit none
    type(meteo_data_source), intent(in) :: met_src
    type(Tfield_buffer), pointer :: metBufPtr
    type(silja_time), dimension(:), intent(in) :: obstimes
    ! meteo -> dispersion interpolation
    type(TVertInterpStruct), pointer :: pVertInterpStruct
    type(THorizInterpStruct), pointer :: pHorizInterpStruct
    ! See also module variables for interpolating winds.
    logical, intent(in) :: ifHorizInterp, ifVertInterp
    type(mini_market_of_stacks), intent(inout) :: dispMarketPtr
    logical, intent(in) :: if_half_levels, if_from_omega, if_top_down

    ! Local variables
    !
    ! Grids of the horizontal winds as in meteo.
    type(silja_grid) :: u_grid, v_grid
    type(silja_time) :: now
    type(silja_time), dimension(max_times) :: dispMarketTimes
    ! fields computed
    real, dimension(:), pointer :: dudx, dvdy, w2d
    real :: dz       ! == h_{i+1/2} - h_{i-1/2}
    real :: rho_adv, dz_full
    integer :: iTime, iLev, i, iu, iv, ip, it, ismall, ismall_prev
    logical :: ifInterpolateUV
    ! Temperature and pressure for calculating time derivative of density. 
    real, dimension(:), pointer :: u2d, v2d, ps, spt, a_half, b_half, a_half_small, b_half_small
    real :: weight_past, da_full, db_full, da_small, db_small, dp
    ! Check if dispersion stack must be arranged
    logical :: ifFirst = .true.
    ! Debug output
    logical :: ifWrite = .false.
    integer :: funit, iTmp, ix, iy, n_sm_times, n_dispMarketTimes, stat, ips, iw
    ! Pressure and temperature fields to calculate d(rho) / dt.
    type(silja_field_id) :: idTmp !, meteoWindID
    type(silja_time) :: time, analysis_time
    type(silja_interval) :: forecast_length, dt

    ps => null()
    spt => fu_work_array()
    a_half => fu_work_array()
    b_half => fu_work_array()
    a_half_small => fu_work_array()
    b_half_small => fu_work_array()
    
    iu = fu_index(metBufPtr%buffer_quantities, u_flag)
    if (fu_fails(iu /= int_missing, 'u-wind not available', 'diagnostic_vertical_wind_eta')) return
    iv = fu_index(metBufPtr%buffer_quantities, v_flag)
    if (fu_fails(iv /= int_missing, 'v-wind not available', 'diagnostic_vertical_wind_eta')) return
    ips = fu_index(metBufPtr%buffer_quantities, surface_pressure_flag)
    if (fu_fails(ips /= int_missing, 'surface pressure not available', 'diagnostic_vertical_wind_eta')) return

    if (if_from_omega) then
      iw = fu_index(metBufPtr%buffer_quantities, omega_flag)
      if (fu_fails(iw /= int_missing, 'omega not available', 'diagnostic_vertical_wind_eta')) return
    end if
    
    call msg('Eta vertical wind diagnosis')

    !
    ! Before starting anything, let's check whether wind has already been diagnosed
    !
    call supermarket_times(dispMarketPtr, met_src_missing, dispMarketTimes, n_dispMarketTimes)
    ! Interpolation of horizontal winds. The possibility of Arakawa
    ! shifts complicates matters, and currently it is only considered
    ! if the dispersion and meteo grids are otherwise equal. In this
    ! case, if ifHorizInterp is false but u_grid /= dispersion_grid,
    ! and separate interpolation structures are set for u and v.

    ! If horiz. interpolated is needed globally, always interpolate u, v.
    ifInterpolateUV = ifHorizInterp

    ! Interpolating the winds
    ! The u/v grids might still be shifted by a half cell.
    u_grid = fu_grid(metBufPtr%p4d(iu)%past%p2d(1)%idPtr)
    v_grid = fu_grid(metBufPtr%p4d(iv)%past%p2d(1)%idPtr)

    if (.not. (u_grid == meteo_grid) .or. .not. (v_grid == meteo_grid)) then
      ! The meteo grid of winds differs from the main one. Create separate 
      ! interpolateion structures.
      if (.not. associated(pHorizInterpU)) then
        call msg('Allocating interpolation for shifted horizontal winds')
        allocate(pHorizInterpU, pHorizInterpV, stat=iTmp)
        if (iTmp /= 0) then
          call set_error('Cannot allocate', 'diag_vertical_wind_eta')
          return
        end if
        pHorizInterpU => fu_horiz_interp_struct(u_grid, dispersion_grid, linear, .true.)
        pHorizInterpV => fu_horiz_interp_struct(v_grid, dispersion_grid, linear, .true.)
      end if
      ifInterpolateUV = .true.
    else
      ! The winds are not shifted from the main meteo grid, which may or 
      ! may not be the dispersion grid. If not, use interpolation as usual.
      pHorizInterpU => pHorizInterpStruct
      pHorizInterpV => pHorizInterpStruct
      ifInterpolateUV = ifHorizInterp
    end if

    if (ifFirst) then 
      call get_integration_vertical()
      allocate(u3d(fs_dispersion, nz_dispersion), v3d(fs_dispersion, nz_dispersion), &
             & eta_dot_3d(fs_dispersion, nz_dispersion+1), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', 'diagnostic_vertical_wind_eta')) return
      if (if_from_omega) allocate(w3d(fs_dispersion, nz_dispersion), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed (omega)', 'diagnostic_vertical_wind_eta')) return
      ifFirst = .false.
    end if



    call hybrid_coefs(dispersion_vertical, a_half=a_half, b_half=b_half)
    call hybrid_coefs(integration_vertical, a_half=a_half_small, b_half=b_half_small)
    forecast_length = zero_interval
    !
    ! Start computations: loop over times. In principle, obstimes can
    ! contain any times, but since values are taken from meteo buffer,
    ! only times defined in buffer are actually considered. Normally
    ! these are the same, of course.
    !
    do itime = 1, size(obstimes)
      now = obstimes(itime)
      if (.not. defined(now)) exit
      if (now == metBufPtr%time_past) then
        weight_past = 1.0
      else if (now == metBufPtr%time_future) then
        weight_past = 0.0
      else
        call set_error('Strange time:'+fu_str(obstimes(itime)),'diag_vertical_wind_eta')
        return
      end if

      analysis_time = now 
      
      call start_count('refine in wind')
      call refine_interp_vert_coefs_v2(wind_interp_struct, metBufPtr, now)
      call stop_count('refine in wind')
      if(error)return
      
      if (.not. if_half_levels) then
        call set_error('very unhappy: not if_half_levels', 'diag_vertical_wind_eta')
        return
      end if

      ! Make the ps field
      !
      idTmp = fu_set_field_id(met_src, &
                            & surface_pressure_flag, &
                            & analysis_time, &
                            & forecast_length, &
                            & dispersion_grid, &
                            & surface_level)
      call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, ps)
      if(fu_fails(.not.error,'Failed u field data pointer','diag_vertical_wind_anelastic'))return
      call get_ps_and_spt(ps, spt, weight_past)
      
      ! Get the winds into 3d cubes. 
      ! 
      u3d = 0.0
      v3d = 0.0
      ismall_prev = 1


      if (if_from_omega) w3d = 0.0

      do ilev = 1, nz_dispersion
        !
        ! Get the fields to update
        !
        idTmp = fu_set_field_id(met_src, &
                              & dispersion_u_flag, &
                              & analysis_time, &
                              & forecast_length, &
                              & dispersion_grid, &
                              & fu_level(dispersion_vertical, iLev, .false.))
        call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, u2d)
        if(fu_fails(.not.error,'Failed u2d field data pointer','diag_vertical_wind_eta'))return
        call set_quantity(idTmp, dispersion_v_flag)
        call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, v2d)
        if(fu_fails(.not.error,'Failed v2d field data pointer','diag_vertical_wind_eta'))return
        call set_quantity(idTmp, dispersion_w_flag)
        call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, w2d)
        if(fu_fails(.not.error,'Failed w2d field data pointer','diag_vertical_wind_eta'))return

        da_full = a_half(ilev) - a_half(ilev+1)
        db_full = b_half(ilev) - b_half(ilev+1)
        do ismall = ismall_prev, ismall_prev + nsmall_in_layer(ilev) - 1
          da_small = a_half_small(ismall) - a_half_small(ismall+1)
          db_small = b_half_small(ismall) - b_half_small(ismall+1)
          call get_winds(ismall, u2d, v2d, ifInterpolateUV)
          do i = 1, fs_dispersion
            dp = (da_small+db_small*ps(i)) / (da_full + db_full*ps(i))
            u3d(i,ilev) = u3d(i,ilev) + u2d(i)*dp
            v3d(i,ilev) = v3d(i,ilev) + v2d(i)*dp
          end do
          if (if_from_omega) then
            call get_omega(ismall, w2d, ifHorizInterp)
            do i = 1, fs_dispersion
              dp = (da_small+db_small*ps(i)) / (da_full + db_full*ps(i))
              w3d(i,ilev) = w3d(i,ilev) + w2d(i)*dp
            end do
          end if
        end do ! small levels
        ismall_prev = ismall
      end do

      ! Diagnose the vertical motion. 
      ! 
      if (.not. if_from_omega) then
        call diagnose_eta_dot(u3d, v3d, ps, spt, &
                            & dispersion_vertical, &
                            & dispersion_grid, dispersion_grid, dispersion_grid, &
                            & .false., if_top_down, eta_dot_3d)
      else
        call eta_from_omega(u3d, v3d, ps, spt, w3d, &
                            & dispersion_vertical, &
                            & dispersion_grid, dispersion_grid, dispersion_grid, &
                            & eta_dot_3d)
      end if

      do iLev = 1, nz_dispersion
        u2d(1:fs_dispersion) = u3d(1:fs_dispersion, ilev)
        v2d(1:fs_dispersion) = v3d(1:fs_dispersion, ilev)
        w2d(1:fs_dispersion) = eta_dot_3d(1:fs_dispersion, ilev+1)
        
        !call msg('Setting test winds...')
        !u(1:fs_dispersion) = 0.0
        !v(1:fs_dispersion) = 0.0!v3d(1:fs_dispersion, ilev)
        !w(1:fs_dispersion) = 0.0!eta_dot_3d(1:fs_dispersion, ilev+1)

      end do ! level
      
!!$      ! To avoid potential inconsistencies and to save time, the
!!$      ! surface pressure is also stored in dispersion grid.
!!$      !
!!$      idtmp = fu_set_field_id(met_src, &
!!$                            & surface_pressure_flag, &
!!$                            & analysis_time, forecast_length, &
!!$                            & dispersion_grid, surface_level)
!!$      call dq_store_2d(dispMarketPtr, idtmp, ps, multi_time_stack_flag, &
!!$                     & iupdateType = overwrite_field, storage_grid=dispersion_grid)
      ! OK gave up, too difficult to make it a dispersion buffer quantity.

    end do ! time

    call free_work_array(spt)
    call free_work_array(a_half)
    call free_work_array(b_half)
    call free_work_array(a_half_small)
    call free_work_array(b_half_small)
    
    call arrange_supermarket(dispMarketPtr)
    
    !==========================================================================================
    
  contains 

    subroutine get_ps_and_spt(ps, spt, weight_past)
      ! Compute d rho / dt on level iLev
      implicit none
      real, dimension(:), intent(out) :: ps, spt
      real, intent(in) :: weight_past 
      ! Local declarations
      integer :: ix, iy
      real :: ps_next, ps_prev
      real :: dt_sec
      
      !dt_sec = fu_sec(fu_valid_time(metBufPtr%p4d(ip)%future%p2d(ilev)%idptr) &
      !     & - fu_valid_time(metBufPtr%p4d(ip)%past%p2d(ilev)%idptr))
      dt_sec = fu_sec(metBufPtr%time_future - metBufPtr%time_past)
      
      if (dt_sec .eps. 0.0) then
        !if (.true.) then
        call msg_warning('dt_sec is zero', 'get_ps_and_spt @ diag_vertical_wind_eta')
        spt(1:fs_dispersion) = 0.0
      end if
      
      do iy = 1, ny_dispersion
        do ix = 1, nx_dispersion

          ps_next = fu_get_value(metBufPtr%p2d(ips), &
                               & nx_meteo, ix, iy, &
                               & 0.0, &
                               & pHorizInterpStruct, ifHorizInterp, .true.)
          ps_prev = fu_get_value(metBufPtr%p2d(ips), &
                               & nx_meteo, ix, iy, &
                               & 1.0, &
                               & pHorizInterpStruct, ifHorizInterp, .true.)

          if (dt_sec > 0.0) spt(ix + (iy-1)*nx_dispersion) = (ps_next - ps_prev) / dt_sec
          ps(ix + (iy-1)*nx_dispersion) = ps_next * (1.0-weight_past) + ps_prev*weight_past
          ! SET CONST PRES
          !ps(ix + (iy-1)*nx_dispersion) = std_pressure_sl
          !spt(ix + (iy-1)*nx_dispersion) = 0.0
        end do
      end do
        
      
    end subroutine get_ps_and_spt

    subroutine get_winds(iLev, data_u, data_v, if_interp_horiz)
      ! 
      ! Get fields full of winds using wind_from_buffer.
      !
      implicit none
      integer, intent(in) :: iLev
      real, dimension(:), pointer :: data_u, data_v
      logical, intent(in) :: if_interp_horiz
      
      integer :: ix, iy
      real :: rand

      !data_u(1:fs_dispersion) = 10.0 !* ilev*0.2
      !data_v(1:fs_dispersion) = 15.0

      do iy = 1, ny_dispersion
        do ix = 1, nx_dispersion
          call wind_from_buffer_4d(metBufPtr%p4d(iu), metBufPtr%p4d(iv), &
                                 & nx_meteo, ix, iy, ilev, &
                                 & weight_past, &
                                 & pHorizInterpU, wind_interp_struct, &
                                 & pHorizInterpV, wind_interp_struct, &
                                 & if_interp_horiz, ifVertInterp, &
                                 & data_u(ix + (iy-1)*nx_dispersion), &
                                 & data_v(ix + (iy-1)*nx_dispersion))

          !call random_number(rand)
          !data_v(ix + (iy-1)*nx_dispersion) = data_v(ix + (iy-1)*nx_dispersion) + rand
          !if (ix > 1) then
          !  data_u(ix + (iy-1)*nx_dispersion) = data_u(ix-1 + (iy-1)*nx_dispersion) - 0.5*(rand-0.5)
          !end if
        end do
      end do
      
    end subroutine get_winds

    subroutine get_omega(ilev, omega, if_interp_horiz)
      implicit none
      integer, intent(in) :: ilev
      real, dimension(:), intent(out) :: omega
      logical, intent(in) :: if_interp_horiz

      integer :: ix, iy

      do iy = 1, ny_dispersion
        do ix = 1, nx_dispersion
          omega(ix + (iy-1)*nx_dispersion) = fu_get_value(metBufPtr%p4d(iw), &
                                                        & nx_meteo, ix, iy, ilev, &
                                                        & weight_past, &
                                                        & pHorizInterpStruct, wind_interp_struct, &
                                                        & if_interp_horiz, ifVertInterp)
        end do
      end do
      
    end subroutine get_omega

    subroutine get_omega_no_small(ilev, omega, if_interp_horiz)
      implicit none
      integer, intent(in) :: ilev
      real, dimension(:), intent(out) :: omega
      logical, intent(in) :: if_interp_horiz

      integer :: ix, iy

      do iy = 1, ny_dispersion
        do ix = 1, nx_dispersion
          omega(ix + (iy-1)*nx_dispersion) = fu_get_value(metBufPtr%p4d(iw), &
                                                        & nx_meteo, ix, iy, ilev, &
                                                        & weight_past, &
                                                        & pHorizInterpStruct, pVertInterpStruct, &
                                                        & if_interp_horiz, ifVertInterp)
        end do
      end do
      
    end subroutine get_omega_no_small

    subroutine get_integration_vertical()
      implicit none
      
      integer :: nsublevs, ilev, lev_ind, nsmall, ilev_new, status, ilev_small_counter
      type(silja_level) :: metlev, displev, new_level, top_small, top_big, bottom_small
      integer, dimension(:), pointer :: metlev_counter
      real :: lev_ind_float, a, b, da_disp, db_disp, interp_data_surf, lev_ind_float_2
      real, dimension(:), pointer :: a_disp, b_disp, interp_data_col

      metlev_counter => fu_work_int_array()
      metlev_counter(1:nz_dispersion) = 0

      interp_data_col => fu_work_array()
      call vert_interp_data_crude(meteo_vertical, dispersion_vertical, interp_data_col, interp_data_surf)
      !print *, interp_data_col(1:nz_dispersion), interp_data_surf
      
      do ilev = 1, nz_meteo
        metlev = fu_level(meteo_vertical, ilev, .true.)
        !displev = fu_level_to_vertical_crude(metlev, dispersion_vertical)
        !lev_ind_float = fu_level_index(displev, dispersion_vertical)
        lev_ind_float = fu_project_level(metlev, dispersion_vertical, interp_data_col, interp_data_surf)
        !print *, lev_ind_float, lev_ind_float_2
        lev_ind = int(lev_ind_float + 0.5)
        if (lev_ind > nz_dispersion) cycle
        metlev_counter(lev_ind) = metlev_counter(lev_ind) + 1
      end do
      call free_work_array(interp_data_col)
      
      if (.not. allocated(nsmall_in_layer)) then
        allocate(nsmall_in_layer(nz_dispersion), stat=status)
        if (status /= 0) then
          call set_error('Allocate failed', 'get_winds')
          return
        end if
      end if

      a_disp => fu_work_array()
      b_disp => fu_work_array()
      call hybrid_coefs(dispersion_vertical, a_half=a_disp, b_half=b_disp)
      
      a = 0.0
      b = 0.0
      ilev_small_counter = 1
      do ilev = nz_dispersion, 1, -1
        nsmall = max(1, metlev_counter(ilev))
        nsmall_in_layer(ilev) = nsmall
        da_disp = a_disp(ilev) - a_disp(ilev+1)
        db_disp = b_disp(ilev) - b_disp(ilev+1)
        top_big = fu_upper_boundary_of_layer(fu_level(dispersion_vertical, ilev))
        top_small = top_big
        a = fu_hybrid_level_coeff_a(top_big)
        b = fu_hybrid_level_coeff_b(top_big)
        do ilev_new = 1, nsmall
          !height = height + thickness/nsmall
          a = a + da_disp/nsmall
          b = b + db_disp/nsmall
          bottom_small = fu_set_hybrid_level(ilev_small_counter, a, b)
          !top_small = fu_set_constant_height_level(height)
          new_level = fu_set_layer_between_two(top_small, bottom_small)
          top_small = bottom_small
          ilev_small_counter = ilev_small_counter + 1
          if (ilev == nz_dispersion .and. ilev_new == 1) then
            call set_vertical(new_level, integration_vertical)
          else
            call add_level(integration_vertical, new_level)
          end if
        end do
      end do

      call arrange_levels_in_vertical(integration_vertical)
      

      wind_interp_struct => fu_vertical_interp_struct(meteo_vertical, integration_vertical, &
                                                    & dispersion_grid, linear, &
                                                    & one_hour, 'wind_interp_eta')

      call msg('Vertical for integrating continuity equation in eta vert')

      call report(integration_vertical, ifFull=.true.)

      call msg('Dispersion vertical:')
      call report(dispersion_vertical, ifFull=.true.)
      
      ! set inter
      ! refine
      call free_work_array(metlev_counter)
      call free_work_array(a_disp)
      call free_work_array(b_disp)
    end subroutine get_integration_vertical

  end subroutine diag_vertical_wind_eta


  !************************************************************************************

  subroutine df_cnc_from_vmr(met_buf, weight_past, &
                           & pHorizInterpStruct, pVertInterpStruct, ifHorizInterp, ifVertInterp, &
                           & vmr_fld_3d, cnc_fld_3d)

    ! Creates concentrations from volume mixing ratio
    ! It is assumed, that the concentration field 3d exists and whatever is there, is overwritten
    ! 

    implicit none

    ! Imported parameters with intent IN:

    type(Tfield_buffer), intent(in) :: met_buf
    real, intent(in) :: weight_past
    type(THorizInterpStruct), pointer :: pHorizInterpStruct
    type(TVertInterpStruct), pointer :: pVertInterpStruct
    logical, intent(in) :: ifHorizInterp, ifVertInterp
    type(silja_3d_field), pointer :: vmr_fld_3d, cnc_fld_3d
    type(silja_field_id), pointer :: id

    ! Local declarations:

!    TYPE(silja_field), POINTER :: cnc_fld
    REAL, DIMENSION(:), POINTER :: volume_mixing_ratio, concentration
    real :: temperature, pressure
!    TYPE(silja_field_id) :: id



    INTEGER :: tIndex, pIndex, ix, iy, iz, iCell, nx, ny 
    type(silam_vertical) :: vertTmp


    call grid_dimensions(fu_grid(vmr_fld_3d), nx, ny)

    tIndex = fu_index(met_buf%buffer_quantities, temperature_flag)
    pIndex = fu_index(met_buf%buffer_quantities, pressure_flag)

    if (.not. all((/tIndex, pIndex/) > 0 )) then
          call msg ("tIndex, pIndex", (/tIndex, pIndex/))
          call set_error("Meteo input is not available", "df_cnc_from_vmr") 
          return
    endif

!    concentration => fu_work_array()

!    call msg_warning('VMR-to-concentration with surface density')

    do iz = 1, fu_number_of_fields(vmr_fld_3d)

      volume_mixing_ratio => fu_grid_data_from_3d(vmr_fld_3d, iz)
      concentration => fu_grid_data_from_3d(cnc_fld_3d, iz)

      do iy = 1, ny
        do ix = 1, nx
          iCell = ix + (iy-1)*nx

          ! Get temperature and pressure
          ! 
          temperature = fu_get_value(met_buf%p4d(tIndex), nx_meteo, ix, iy, iz, &
                                    & weight_past, &
                                    & pHorizInterpStruct, pVertInterpStruct, &
                                    & ifHorizInterp, ifVertInterp)
          IF (error) RETURN
          pressure = fu_get_value(met_buf%p4d(pIndex), nx_meteo, ix, iy, iz, &
                                & weight_past, &
                                & pHorizInterpStruct, pVertInterpStruct, &
                                & ifHorizInterp, ifVertInterp)
          IF (error) RETURN
          
          !Compute the concentration

          concentration(iCell) = volume_mixing_ratio(iCell) *pressure /( gas_constant_uni * temperature)

!          concentration(iCell) = volume_mixing_ratio(iCell) / 22.4e-3
 
        enddo
      enddo
      

!      id => fu_id(fu_field_from_3d_field(vmr_fld_3d, iz))
!      call set_quantity(id, concentration_flag)
!      call set_quantity(id, concentration_flag)
!      call set_field(id, concentration, cnc_fld)
!      call add_field_to_3d_field(cnc_fld, field_3d)

    enddo

!    CALL free_work_array(concentration)

  END SUBROUTINE df_cnc_from_vmr


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
    iQ = fu_index(mdl_in_q, total_precipitation_rate_flag)
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
  
  subroutine df_cumul_daily_tempr_2m( met_buf,  dispMarketPtr, ifColdstartAllowed)
    !
    ! Computes the cumulative daily temperature from instant one
    ! Note that allowColdStart flag should be handled on the _second_ call
    ! As a way to detect the initialization we store the zero-accumualted
    !  past-field if allowColdStart==false.  The accumulation length   
    !
    implicit none

    ! Imported parameters
    TYPE(Tfield_buffer), POINTER :: met_buf
    type(mini_market_of_stacks), intent(inout) :: dispMarketPtr
    logical, intent(in) :: ifColdstartAllowed

    ! Local variables
    INTEGER :: iQ, ix, iy, iTo, iMeteo, iTime, fsTo, iLev, nLevs, nxTo, nyTo, iCoef, iTmp, indT2m
    real, DIMENSION(:), POINTER :: pTemprPast, pTemprFuture, pTCum, pTCum_past
    TYPE(Tfield_buffer), POINTER :: met_bufPtr
    type(silja_field_id) :: idTmp
    type(silja_field_id), pointer :: idTemprPast, idTemprFuture, idPastPtr, idFuturePtr
    type(silja_time) :: start_of_day, timeTmp
    type(silja_interval) :: intervalTmp, accLen
    logical :: ifAddPAst, ifPastMadeHere
    real :: seconds
    real, dimension(:), pointer :: pData
    type(silja_field), pointer :: fldPtr
    type(THorizInterpStruct), pointer :: pHorizInterpStruct
    character(len=20), parameter :: sub_name = 'df_cumul_daily_tempr'


    pHorizInterpStruct => fu_horiz_interp_struct(meteo_grid, dispersion_grid, linear, .true.)
    met_bufPtr => met_buf
    call grid_dimensions(dispersion_grid, nxTo, nyTo)
    fsTo = nxTo * nyTo
    !

    indT2m =  fu_index(met_buf%buffer_quantities, temperature_2m_flag)
    pTemprPast => met_buf%p2d(indT2m)%past%ptr
    pTemprFuture => met_buf%p2d(indT2m)%future%ptr
    idTemprPast => met_buf%p2d(indT2m)%past%idPtr
    idTemprFuture => met_buf%p2d(indT2m)%future%idPtr


    !
    !   Detect cold/warm start, otherwise past should be ready: 
    !  day_temperature_2m_acc_flag must never have zero accumulation length
    !
    timeTmp = fu_valid_time(idTemprPast)
    idTmp = fu_set_field_id_simple(fu_met_src(idTemprPast), &
                                 & day_temperature_2m_acc_flag, timeTmp, level_missing) 
    fldPtr => fu_get_field_from_mm_general(dispMarketPtr, idTmp, .false.)

   ! Should be previous day for midnight
    start_of_day = fu_start_of_day_utc(timeTmp - one_second)
    accLen = timeTmp - start_of_day !! Should never be zero!

    if (.not. associated(fldPtr)) then !! Coldstart
          ! Create something here as mean from past and future temperature
          ! accumulated since day start
          ifPastMadeHere = .true.

          idTmp = fu_set_field_id(fu_met_src(idTemprPast),&
                                & day_temperature_2m_acc_flag, &
                                & start_of_day, &  !analysis_time,&
                                & accLen, & !forecast_length, &
                                & dispersion_grid, &
                                & fu_level(idTemprPast), &
                                & accLen, & ! & !length_of_accumulation, optional
                                & zero_interval, &          !length_of_validity, optional
                                & accumulated_flag)         !field_kind, optional

          call msg_warning("Making surrogate past day_temperature_2m_acc_flag",sub_name)
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
                  pTCum_past(iTo) = 0.
                  do iCoef = 1, pHorizInterpStruct%nCoefs
                    iTmp = pHorizInterpStruct%indX(iCoef,ix,iy) + &
                                          & (pHorizInterpStruct%indY(iCoef,ix,iy) - 1) * nx_meteo
                    pTCum_past(iTo) = pTCum_past(iTo) + pHorizInterpStruct%weight(iCoef,ix,iy) * &
                                                  & (pTemprPast(iTmp) + pTemprFuture(iTmp))
                  enddo
                  pTCum_past(iTo) = pTCum_past(iTo) * seconds * 0.5
              end do  ! ix
            end do  ! iy
          else
            pTCum_past(1:fsTo) = (pTemprPast(1:fsTo) + pTemprFuture(1:fsTo)) * 0.5 * seconds
          endif  ! if horizontal interpolation
    else
      idPastPtr => fu_id(fldPtr) 
      pTCum_past => fu_grid_data(fldPtr)
      intervalTmp = fu_accumulation_length(idPastPtr)
      ifPastMadeHere = .false.

      !!Should be valid field
      if (pTCum_past(1) == real_missing) then
         if (ifColdStartAllowed) then
            call  set_error("This should not happen", sub_name)
            return
         else
            call msg_warning("Looks like day_temperature_2m_acc_flag was not initialised. crashing..", sub_name)
            call msg("This behaviour can be overridden by 'allow_coldstart_day_temperature = yes' in standard_setup")
            call set_error("Missing values for past", sub_name)
            return
         endif

      endif

      ! Correct accumulation that might have been broken from the initialization
      if (.not. intervalTmp == accLen) then
          call msg_warning("Correcting acc_length for past day_temperature_2m_acc_flag", sub_name)
          call msg("Initial ID")
          call report(idPastPtr)
          call set_accumulation_length(idPastPtr, accLen)
          call msg("Corrected ID")
          call report(idPastPtr)
      endif
      idTmp = idPastPtr
    endif

    !
    ! Past should be okay. Future now
    ! We should make it in any case: past might have been reset with ini
    
    ! Should be previous day for midnight
    timeTmp = fu_valid_time(idTemprFuture)
    start_of_day = fu_start_of_day_utc(timeTmp - one_second)
    intervalTmp = timeTmp - start_of_day !! Should never be zero!


    ifAddPast = (fu_hour(fu_valid_time(idTmp)) > 0) 

    call set_valid_time(idTmp, timeTmp)
    call set_accumulation_length(idTmp,intervalTmp)

    ! On the second diagnostics we should get the existig field
    fldPtr => fu_get_field_from_mm_general(dispMarketPtr, idTmp, .false.)
    if (associated(fldPtr)) then
      idFuturePtr => fu_id(fldPtr)
      idFuturePtr = idTmp
      pTCum => fu_grid_data(fldPtr)
    else
      !This guy always creates a field 
      call find_field_data_storage_2d(dispMarketPtr, idTmp, multi_time_stack_flag, pTCum)
      if(fu_fails(.not.error,'Failed field data pointer', sub_name))return
    endif

    pTCum(1:fsTo) =  0.
    if (pHorizInterpStruct%ifInterpolate) then
      do iy = 1, nyTo
        do ix = 1, nxTo
          iTo = ix + (iy-1) * nxTo
            do iCoef = 1, pHorizInterpStruct%nCoefs
              iTmp =  pHorizInterpStruct%indX(iCoef,ix,iy) + &
                   & (pHorizInterpStruct%indY(iCoef,ix,iy)-1) * nx_meteo
              pTCum(iTo) = pTCum(iTo) + pHorizInterpStruct%weight(iCoef,ix,iy) * &
                                      & (pTemprPast(iTmp) + pTemprFuture(iTmp))
            enddo
        end do
      end do
    else
        pTCum(1:fsTo) = pTemprPast(1:fsTo) + pTemprFuture(1:fsTo)
    endif

    !! Finalize accumulation and add earlier accumulation if needed
    seconds = fu_sec(fu_valid_time(idTemprFuture) - fu_valid_time(idTemprPast))
    if (ifAddPast) then
        pTCum(1:fsTo) = pTCum(1:fsTo)*0.5*seconds + pTCum_past(1:fsTo) 
    else
        pTCum(1:fsTo) = pTCum(1:fsTo)*0.5*seconds
    endif


    !
    ! Guard for missing initialisation
    if (ifPastMadeHere .and. (.not. ifColdStartAllowed)) then
       pTCum_past(1:fsTo) = real_missing
    endif


  end subroutine df_cumul_daily_tempr_2m
  
  
  !***************************************************************************************************

  subroutine df_update_realtime_met_fields(pMMarket, buf_ptr, now, wdr, diag_rules, &
                                         & valid_time_border_1, valid_time_border_2, ifResetTimes)
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
    logical, intent(in) :: ifResetTimes

    !Local 
    type(silja_field), pointer :: field
    type(silja_field_id), pointer :: idIn
    type(silja_field_id) :: idRequest
    logical :: ifvalid, ifFound
    integer :: iVar, shopQ, Qidx

    !
    ! List of realtime quantities that we know how to update
    ! Huom! Order matters:  No validity checks in the deriving routines.
    !
    integer, parameter, dimension(1:10) :: st_update_list = &
       & (/ground_pressure_flag, large_scale_rain_int_flag, convective_rain_int_flag, &
            & total_precipitation_rate_flag, scavenging_coefficient_flag,& 
            & disp_flux_cellt_rt_flag, disp_flux_celle_rt_flag, disp_flux_celln_rt_flag, &
            & day_mean_temperature_flag, day_mean_temperature_2m_flag/)

    call msg("Updating single-time fields for now=" + fu_str(now))
    valid_time_border_1 = really_far_in_past
    valid_time_border_2 = really_far_in_future

    do iVar = 1, size(st_update_list)
      shopQ = st_update_list(iVar)
      if (shopQ == ground_pressure_flag) then
        ! only multitime ground_pressure_field is accesible from the buffer.
        ! Should use stack directly
        call find_field_from_stack(met_src_missing, &
                                 & ground_pressure_flag,&
                                 & time_missing,& ! any time 
                                 & fu_stack(pMMarket, 1),& !ST stack
                                 & field,&
                                 & ifFound)
        if (.not. ifFound) then
          call msg("Not found in stack:"+ fu_name(fu_stack(pMMarket, 1)) )
          cycle
        endif
        ! call msg("Found in stack:"+ fu_name(fu_stack(metMarket, 1)) )
        ! call report(fu_stack(metMarket, 1))
        idIn => fu_id(field)
      else
        Qidx  = fu_index(buf_ptr, shopQ) 
        if (Qidx < 1) cycle 
             
        if (fu_dimension(buf_ptr, Qidx) == 4 ) then 
          idIn => buf_ptr%p4d(Qidx)%present%p2d(1)%idPtr
        else
          idIn => buf_ptr%p2d(Qidx)%present%idPtr
        endif
      endif
!call msg(fu_quantity_string(shopQ))
      idRequest = fu_set_field_id_simple(met_src_missing, shopQ, now, level_missing)
      ifvalid = fu_field_id_covers_request(idIn, idRequest, .true.)
      if (error) return
!        call msg("---------------")
!        call msg("requested ID:")
!        call report(idRequest)
!        call msg("found ID:")
!        call report(idIn)
!        call ooops("Achtung!")


      if (ifvalid) then
        call msg("Still valid:"+ fu_quantity_string(shopQ))
        if(valid_time_border_1 < fu_valid_time(idIn)) valid_time_border_1 = fu_valid_time(idIn)
        if(valid_time_border_2 > fu_valid_time(idIn) + fu_validity_length(idIn)) &
                        & valid_time_border_2 = fu_valid_time(idIn) + fu_validity_length(idIn)
        cycle ! no need to update
      endif

call msg("Updating realtime meteo field:"+ fu_quantity_string(shopQ))

      select case (shopQ)
        case (ground_pressure_flag)
            CALL df_update_ground_pressure(field, buf_ptr, &
                                         & valid_time_border_1, valid_time_border_2) 
! call report(fu_id(field))
        case (large_scale_rain_int_flag, convective_rain_int_flag, &
                    &  total_precipitation_rate_flag)
            CALL df_update_rain_intensity(buf_ptr, shopQ, wdr, &
                                        & valid_time_border_1, valid_time_border_2) 

        case (scavenging_coefficient_flag)
            CALL df_update_scav_coefficient(buf_ptr, valid_time_border_1, valid_time_border_2)

        case (disp_flux_cellt_rt_flag, disp_flux_celle_rt_flag, disp_flux_celln_rt_flag)
             ! Should be called for dispersion buffer only
            CALL df_update_cellfluxcorr_rt(buf_ptr, valid_time_border_1, valid_time_border_2, &
                                         & diag_rules, now)
            call exchange_wings(buf_ptr)

        case(day_mean_temperature_2m_flag)
            if(fu_get_time_direction_sm(pMMarket) == backwards)then
              call set_error('Mean daily and cumulative temperature cannot be made in adjoint runs', &
                           & 'df_update_realtime_met_fields')
              return
            endif
            call df_daily_mean_tempr_2m(buf_ptr, now, valid_time_border_1, valid_time_border_2, &
               & diag_rules%ifAllowColdstartDayTemperature)

        case default
          call msg_warning('Do not know how to update:' + fu_quantity_string(shopQ) + ', skipping...', &
                                   & 'df_update_realtime_met_fields')

      end select
    enddo  ! Updatable meteo Var list
    
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
      case (total_precipitation_rate_flag)
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
             case(7,8,9)     
               
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

  subroutine df_daily_mean_tempr_2m(buf, now, valid_time_border_1, valid_time_border_2, ifColdstartAllowed)
    ! 
    ! Mean temperature over _previous_ UTC day (00Z-00Z)
    ! Field set as average from 00Z-00Z for today
    !
    implicit none

    ! Imported parameters
    type(Tfield_buffer), INTENT(in) :: buf   ! actually, dispersion buffer
    type(silja_time) :: now !! Used only for temporary solution
    type(silja_time), intent(inout) :: valid_time_border_1, valid_time_border_2
    logical, intent(in) :: ifColdstartAllowed

    ! Local variables
    TYPE(silja_field_id), pointer :: pIdCum_past, pIdCum_future, pIdMean
    REAL, DIMENSION(:), POINTER :: cum_T_past, cum_T_future, mean_T
    type(silja_time) :: time_past, time_future, start_of_day
    type(silja_time) :: time_tmp
    integer :: iT_mean, iT_cum, iLev, nLevs, nx, ny, fs
    real :: seconds
    logical :: ifMakeIt, ifMakeSurrogate, stillValid
    character(len=19), parameter :: sub_name = 'df_daily_mean_tempr_2m'
    
    integer, save :: callCnt = 0


    callCnt = callCnt + 1
    !
    ! Do we need to do anything?
    ! is mean temperature field valid for the time range between pst and future?
    !
    ! Preparatory steps
    !
    iT_cum  = fu_index(buf%buffer_quantities, day_temperature_2m_acc_flag)
    iT_mean = fu_index(buf%buffer_quantities, day_mean_temperature_2m_flag)
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
      call msg("Still valid: "//fu_quantity_string(day_mean_temperature_2m_flag))
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
    if (ifMakeIt) then ! Just reset it
          call msg("Recalculating daymean T")
          call msg("Before mean_T(1:10)", mean_T(1:10))
          mean_T(1:fs) = cum_T_past(1:fs) / seconds_in_day
          call msg("After mean_T(1:10)", mean_T(1:10))
    elseif (mean_T(1) == real_missing) then  !fingerprint of gobal_io_init
      if (CallCnt == 1) then  !Field was just created by gobal_io_init
        if (ifColdstartAllowed) then
          ! Put coldstart values here, so next time will not crash on missing ini
          call msg("Filling daymean T with something meaningful for the first out")
          mean_T(1:fs) = cum_T_past(1:fs) / fu_sec(fu_accumulation_length(pIdCum_past))
        else
          call msg("Leaving real_missing daymean T ")
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
    call msg("before")
    call report(pIdMean)

     !! It is just an average field for the whole day
    start_of_day = fu_start_of_day_utc(time_past)
    pIdMean =  fu_set_field_id(met_src_missing,&
                                     & day_mean_temperature_2m_flag, &
                                     & start_of_day,& ! Analysis
                                     & one_day, &             ! Forecast length 
                                     & fu_grid(pIdCum_past),&
                                     & fu_level(pIdCum_past),&
                                     & one_day, & ! optional accumulation interval
                                     & zero_interval, & ! !! We should get rid of validity_length at some point
                                     & averaged_flag) !field_kind
    call msg ("After")                                                     
    call report(pIdMean)
    !call ooops("")

    if (callCnt < 2) then !! Forece recalculation on the second call
      time_tmp = start_of_day
    else
      time_tmp = start_of_day + one_day
    endif
      if(valid_time_border_1 < start_of_day) valid_time_border_1 = start_of_day
      if(valid_time_border_2 > time_tmp) valid_time_border_2 = time_tmp




  end subroutine df_daily_mean_tempr_2m


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

   ! No neigbours -- no wings
   if (all(adv_mpi_neighbours(:) < 0)) return

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
    if(adv_mpi_neighbours(northern_boundary)>=0)then
      call smpi_exchange_wings(adv_mpi_neighbours(northern_boundary), sendN_WAptr, nxSlice, &
           & recvN_WAptr, recv_count)
      disp_buf%wings_N(1:3,1:2,1:nz_disp,1:nx_disp,1:nDepth) = &
           & reshape(recvN_WAptr(1:nxSlice), (/3,2,nz_disp,nx_disp,nDepth/))
        else
      disp_buf%wings_N(1:3,1:2,1:nz_disp,1:nx_disp,1:nDepth) = -8888
        endif
    if(adv_mpi_neighbours(southern_boundary)>=0)then
      call smpi_exchange_wings(adv_mpi_neighbours(southern_boundary), sendS_WAptr, nxSlice, &
           & recvS_WAptr, recv_count)
      disp_buf%wings_S(1:3,1:2,1:nz_disp,1:nx_disp,1:nDepth) = &
           & reshape(recvS_WAptr(1:nxSlice), (/3,2,nz_disp,nx_disp,nDepth/))
    else
      disp_buf%wings_S(1:3,1:2,1:nz_disp,nx_disp,nDepth) = -8888
    endif

     ! Exchange the boundaries in East-West direction
    if(adv_mpi_neighbours(eastern_boundary)>=0)then
      call smpi_exchange_wings(adv_mpi_neighbours(eastern_boundary), sendE_WAptr, nySlice, &
           & recvE_WAptr, recv_count)
      disp_buf%wings_E(1:3,1:2,1:nz_disp,1:nDepth,1:ny_disp) = &
           & reshape(recvE_WAptr(1:nySlice), (/3,2,nz_disp,nDepth,ny_disp/))
        else
      disp_buf%wings_E(1:3,1:2,1:nz_disp,1:nDepth,1:ny_disp) = -8888
        endif
    if(adv_mpi_neighbours(western_boundary)>=0)then
      call smpi_exchange_wings(adv_mpi_neighbours(western_boundary), sendW_WAptr, nySlice, &
           & recvW_WAptr, recv_count)
      disp_buf%wings_W(1:3,1:2,1:nz_disp,1:nDepth,1:ny_disp) = &
           & reshape(recvW_WAptr(1:nySlice), (/3,2,nz_disp,nDepth,ny_disp/))
    else
      disp_buf%wings_W(1:3,1:2,1:nz_disp,1:nDepth,1:ny_disp) = -8888
    endif

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
