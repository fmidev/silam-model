module chem_dep_cbm5_strato_SOA
  !
  ! This module contains full description of the DMAT-based SOx chemistry and deposition.
  ! To be exact, there is a description of the SO2 and SO4=, removal computation - both
  ! dry and wet.
  !
  ! One type, which explains how this stuff has to be treated is Tchem_rules_S_DMAT,
  ! which contains all necessary rules for the operation, including the types
  ! of the functions used for the computations. This type is used 
  ! only here and is not seen from outside, except for the command "create".
  !
  ! All units: SI, unless otherwise stated
  !
  ! Code owner Mikhail Sofiev mikhail.sofiev@fmi.fi
  !
  ! Language: ANSI FORTRAN-90 (or close to it)
  !
  use cbm5_strato_SOA_interface
  use cocktail_basic
  use cbm5_strato_SOA_Parameters !R.H. 10Jan2018 for adjusting the accuracy for different species!
  !$use omp_lib

  implicit none
  private

  public set_chem_rules_cbm5_strato_SOA
  public set_missing
  
  public fu_lifetime
  public check_low_mass_threshold_cbm5_strato_SOA
  public registerSpeciescbm5_strato_SOA
  public inventory_cbm5_strato_SOA
  public transform_cbm5_strato_SOA
  public transform_cbm5_strato_SOA_adj

  public cbm5_strato_SOA_input_needs
  public fu_if_specific_deposition
  public prepare_step_cbm5_strato_SOA
  
  ! inherited from the interface
  public init_chemicals_cbm5_strato_SOA

  !
  ! Private routines of the sulphur DMAT cocktail and chemical rules
  !
  private set_miss_chem_rules_cbm5_strato_SOA

  private fu_lifetime_cbm5_strato_SOA
  private fu_if_specific_deposition_cbm5_strato_SOA

  interface set_missing
    module procedure set_miss_chem_rules_cbm5_strato_SOA
  end interface

  interface fu_lifetime
    module procedure fu_lifetime_cbm5_strato_SOA
  end interface

  interface fu_if_specific_deposition
    module procedure fu_if_specific_deposition_cbm5_strato_SOA
  end interface

  public fu_if_tla_required
  private fu_if_tla_required_cbm5_strato_SOA
  interface fu_if_tla_required
     module procedure fu_if_tla_required_cbm5_strato_SOA
  end interface


  type Tchem_rules_cbm5_strato_SOA
     private
     logical :: defined = .false.
     integer, dimension(20) :: kpp_icntrl
     ! solver relative tolerance
     real(precision_cbm5_strato_SOA), dimension(num_species_cbm5_strato_SOA) :: kpp_rtol
     ! solver absolute tolerance, solver units
     real(precision_cbm5_strato_SOA), dimension(num_species_cbm5_strato_SOA) :: kpp_atol
     ! the rcntrl doesn't work here: it is updated per cell
     !real, dimension(num_species_cbm5_strato_SOA) :: kpp_rcntrl
     real :: ASOA_aging, BSOA_aging, IVOC_aging, ifMTprods !NOTE: SOA specific
     real :: MCF_factor !Extra factor for CH3CCl3+OH reaction rate!
  end type Tchem_rules_cbm5_strato_SOA
  
  public Tchem_rules_cbm5_strato_SOA


  !
  ! The mark of this transformation module:
  !
  integer, parameter, public :: transformation_cbm5_strato_SOA = 5014
  ! an argument to the solver. Set initially except for the parameter which defined
  ! initial internal step: this is set for each cell according to the previous value.
  real(precision_cbm5_strato_SOA), dimension(20), private, save :: kpp_rcntrl
  !$OMP THREADPRIVATE(kpp_rcntrl)

  real, parameter, private :: mol2molec = avogadro*1e-6

  integer, pointer :: ind_tempr, ind_press, ind_cloud_covr, ind_spec_hum
  integer, private, save :: iNO, iXO2, iHO2 !NOTE: SOA specific


CONTAINS



  !************************************************************************************

  subroutine prepare_step_cbm5_strato_SOA(if_report_rates, photoarr)
    implicit none
    logical, optional, intent(in) :: if_report_rates
    real, dimension(:), optional, intent(in) :: photoarr
    character(len=*), parameter :: sub_name = 'prepare_step_cbm5_strato_SOA'
    real :: cAir
    
    call set_const_rates_cbm5_strato_SOA()
    kpp_rcntrl = 0.0
    kpp_rcntrl(1) = 1.0e-3 ! 1 ms minimum step length
    if (present(if_report_rates)) then
      if (fu_fails(present(photoarr), 'Need photoarr to report_rates', sub_name)) return
      if (if_report_rates) then
        ! arguments below need to be matched with those in kpp_interface!
        cAir = 101325.0 / (300.0 * gas_constant_uni) * mol2molec 
        call set_rates_cbm5_strato_SOA(photoarr, temp=300.0, cAir=cAir, press=101325.0, cH2O=0.01&
                             &, soa_b=0.0, ASOA_aging=0.0, BSOA_aging=0.0, IVOC_aging=0.0, ifMTprods=0.0&
                             &, MCF_factor=1.0&
	                       &)	
        call report_rates_cbm5_strato_SOA()
      end if
    end if
  end subroutine prepare_step_cbm5_strato_SOA


  subroutine inventory_cbm5_strato_SOA(rules, &
                               & speciesEmis, speciesTransp, speciesShortlived, speciesAerosol,&
                               & nSpeciesEmis, nSpeciesTransp, nSpeciesShortlived, nSpeciesAerosol, &
                               & iClaimedSpecies)
    implicit none
    !
    ! The list of species + claiming the species. Call the routine from interface, but
    ! without rules which is not known there.
    type(Tchem_rules_cbm5_strato_SOA), intent(in) :: rules
    type(silam_species), dimension(:), pointer :: speciesEmis, speciesTransp, speciesShortlived, &
                                                & speciesAerosol
    integer, intent(in) :: nSpeciesEmis
    integer, intent(inout) :: nspeciesTransp, nspeciesShortlived, nspeciesAerosol
    integer, dimension(:), intent(inout) :: iClaimedSpecies

    ! Local variables
    integer :: iEmis

    call inventory_int(transformation_cbm5_strato_SOA, speciesEmis, &
                     & speciesTransp, speciesShortlived, speciesAerosol,&
                     & nSpeciesEmis, nSpeciesTransp, nSpeciesShortlived, nSpeciesAerosol, &
                     & iClaimedSpecies)

  end subroutine inventory_cbm5_strato_SOA


  !************************************************************************************

  subroutine registerSpeciescbm5_strato_SOA(rules, &
                                    & speciesTransp, speciesShortlived, speciesAerosol,&
                                    & nspeciesTransp, nspeciesShortlived, nspeciesAerosol)
    !
    ! The integrator requires fixed indices for the species. At this
    ! point, we just check that the transport species follow this ordering.
    implicit none
    type(Tchem_rules_cbm5_strato_SOA), intent(in) :: rules
    type(silam_species), dimension(:), pointer :: speciesTransp, speciesShortlived, &
                                                & speciesAerosol
    integer, intent(in) :: nspeciesTransp, nspeciesShortlived, nspeciesAerosol
    
    ! Local variable
    integer :: i

    call register_species_int(speciesTransp, speciesShortlived, speciesAerosol,&
    !                     & nspeciesTransp, nspeciesShortlived, nspeciesAerosol)                !Without SOA
                         & nspeciesTransp, nspeciesShortlived, nspeciesAerosol,iNO, iXO2, iHO2) !With SOA
    
  end subroutine registerSpeciescbm5_strato_SOA


  !*****************************************************************

  subroutine cbm5_strato_SOA_input_needs(rules, metdat)
    !
    ! Returns input needs for the transformation routine. Change if needed.  
    ! 
    implicit none

    ! Imported parameters
    type(Tchem_rules_cbm5_strato_SOA), intent(in) :: rules
    type(Tmeteo_input), intent(out), target :: metdat

    ! Local variables
    integer :: iStat, nq

    if (.not. rules%defined) then
      call set_error('Undefined cbm5_strato_SOA rules', 'cbm5_strato_SOA_input_needs')
      return
    end if

    metdat = meteo_input_empty
    metdat%nQuantities = 4

    metdat%quantity(1) = temperature_flag
    metdat%q_type(1) = meteo_dynamic_flag
    ind_tempr => metdat%idx(1)

    metdat%quantity(2) = pressure_flag
    metdat%q_type(2) = meteo_dynamic_flag
    ind_press => metdat%idx(2)

    metdat%quantity(3) = total_cloud_cover_flag
    metdat%q_type(3) = meteo_dynamic_flag
    ind_cloud_covr => metdat%idx(3)
    
    metdat%quantity(4) = specific_humidity_flag
    metdat%q_type(4) = meteo_dynamic_flag
    ind_spec_hum => metdat%idx(4)

    metdat%defined = silja_true

  end subroutine cbm5_strato_SOA_input_needs

  !***********************************************************************

  subroutine set_miss_chem_rules_cbm5_strato_SOA(rules)
    implicit none
    type(Tchem_rules_cbm5_strato_SOA), intent(out) :: rules

    rules%defined = .false.

  end subroutine set_miss_chem_rules_cbm5_strato_SOA
  
  !***********************************************************************

  subroutine set_chem_rules_cbm5_strato_SOA(nlSetup, rules)
    !
    ! Creates the set of rules, in accordance with which the computations are to
    ! be made. Input - SILAM namelist
    ! Taking this chance, let's also set once-and-forever a few chemical coefficients
    !
    implicit none

    ! Imported parameter
    type(Tsilam_namelist), pointer :: nlSetup 
    type(Tchem_rules_cbm5_strato_SOA) :: rules

    ! Local variables
    integer :: i, jTmp, ierr
    real :: fTempr, z, ratio_K_prod, ratio_K_loss, f300_T, fRelatDensity, fRelatPress, fRelatTempr
    CHARACTER(len=*), PARAMETER :: sub_name = 'set_chem_rules_cbm5_strato_SOA'

    !
    ! Stupidity checking
    !
    rules%defined = .false.
    if(.not.defined(nlSetup))then
      call set_error('Undefined namelist given', sub_name)
      return
    endif
    
    rules%kpp_icntrl = 0 ! all default, except:
    rules%kpp_icntrl(1) = 1 ! autonomous system
    !rules%kpp_icntrl(4) = 900*1000 ! max number of sub-steps: enough for 15 minutes

    ! The error tolerances. Inside the KPP intergrator, the error is computed something
    ! like err = |y_i+1 - y_i| / scale, where scale = atol + rtol*y.  We can set rtol as
    ! high as we tolerate, but it seems that atol has better be small in order to avoid
    ! negative concentrations.
    rules%kpp_atol = 3e-15*mol2molec  !R.H. New more accurate values. Otherwise negative garbage may appear.
    rules%kpp_rtol = 0.0001           !R.H. New more accurate values. Otherwise negative garbage may appear.
    !R.H. Improving the atol for specific species. NOTE: These are specific for strato/strato_SOA.
    rules%kpp_atol(ind_Br2)  = 1e-18*mol2molec  !Only for strato/strato_SOA
    rules%kpp_atol(ind_Br)   = 1e-17*mol2molec  !Only for strato/strato_SOA
    rules%kpp_atol(ind_BrCl) = 1e-17*mol2molec  !Only for strato/strato_SOA
    rules%kpp_atol(ind_Cl2)  = 1e-16*mol2molec  !Only for strato/strato_SOA
    rules%kpp_atol(ind_OClO) = 2e-17*mol2molec  !Only for strato/strato_SOA
    rules%kpp_atol(ind_O1D)  = 1e-16*mol2molec
    rules%kpp_atol(ind_NO)   = 1e-16*mol2molec
    rules%kpp_atol(ind_NO2)  = 1e-16*mol2molec
    rules%kpp_atol(ind_C2O3) = 5e-17*mol2molec
    rules%kpp_atol(ind_BrO)  = 1e-16*mol2molec   !Only for strato/strato_SOA
    rules%kpp_atol(ind_Cl)   = 1e-17*mol2molec   !Only for strato/strato_SOA
    rules%kpp_atol(ind_O)    = 1e-15*mol2molec
    !!!rules%kpp_atol(ind_BrNO2)= 1e-16*mol2molec
    rules%kpp_atol(ind_Cl2O2)= 1e-15*mol2molec   !Only for strato/strato_SOA
    select case (fu_content(nlSetup, 'cbm_tolerance'))
         case ('FAST')
            call msg("Using fast (aka old v55 tolerances) for cbm5_strato_SOA")
            rules%kpp_atol = 1e-12*mol2molec !old ver5.5
            rules%kpp_rtol = 0.01            !old ver5.5
         case ('LOW')
            call msg("Using low tolerances for cbm5_strato_SOA")
            rules%kpp_atol = 3e-13*mol2molec
            rules%kpp_rtol = 0.003
         case ('MODERATE')
            call msg("Using moderate tolerances for cbm5_strato_SOA")
            rules%kpp_atol = 1e-13*mol2molec
            rules%kpp_rtol = 0.001
         case ('HIGH')
            call msg("Using high tolerances for cbm5_strato_SOA")
            rules%kpp_atol = 2e-14*mol2molec
            rules%kpp_rtol = 0.0003
         case ('STRICT')
            call msg("Using strict tolerances for cbm5_strato_SOA")
         case ('')
           call msg_warning("cbm_tolerance can be FAST or STRICT")
           call msg("STRICT tolerances are used")
         case default
           call set_error("Strange cbm_tolerance", sub_name)
    end select

    call init_solver(rules%kpp_icntrl, ierr)
    if (fu_fails(ierr == 0, 'Error initializing solver', sub_name)) return
    call init_solver_adj(rules%kpp_icntrl, ierr)
    if (fu_fails(ierr == 0, 'Error initializing adjoint solver', sub_name)) return

    !Next are specific to SOA only:
    rules%BSOA_aging = fu_content_real(nlSetup, 'biogenic_SOA_aging_rate')
    if(rules%BSOA_aging .eps. real_missing)then
      call set_error('Failed to find biogenic_SOA_aging_rate', sub_name)
      return
    endif
    call msg('Aging rate for biogenic SOA', rules%BSOA_aging)

    rules%ASOA_aging = fu_content_real(nlSetup, 'anthropogenic_SOA_aging_rate')
    if(rules%ASOA_aging .eps. real_missing)then
      call set_error('Failed to find anthropogenic_SOA_aging_rate', sub_name)
      return
    endif
    call msg('Aging rate for anthropogenic SOA', rules%ASOA_aging)

    rules%IVOC_aging = fu_content_real(nlSetup, 'intermediate_volatility_OC_aging_rate')
    if(rules%IVOC_aging .eps. real_missing)then
      call set_error('Failed to find intermediate_volatility_OC_aging_rate', sub_name)
      return
    endif
    call msg('Aging rate for intermediate volatility OC', rules%IVOC_aging)

    rules%ifMTprods = fu_content_real(nlSetup, 'if_monoterpene_products')
    if(rules%ifMTprods .eps. real_missing)then
      call set_error('Failed to find if_monoterpene_products', sub_name)
      return
    endif
    call msg('Monoterpenes oxydation products', rules%ifMTprods)
    !END: Next are specific to SOA only.

    !Since the OH levels are typically too low in SILAM we make it possible to have an extra factor
    !for the rate of the reaction OH+CH3CCl3 -> products. Needed only for strato.
    rules%MCF_factor = fu_content_real(nlSetup, 'methylchloroform_OH_rate_factor')
    if(rules%MCF_factor .eps. real_missing)then
      call msg('Failed to find methylchloroform_OH_rate_factor! Using default value 1.0.')
      rules%MCF_factor = 1.0
    endif
    call msg('Extra rate factor for methylchloroform + OH reaction:', rules%MCF_factor)    

    rules%defined = .true.

  end subroutine set_chem_rules_cbm5_strato_SOA


  !***********************************************************************
  subroutine transform_cbm5_strato_SOA(vSp, photo, rules, metdat, timestep_sec, garbage, &
                               & zenith_cos, lat, h_start, now, print_it)

    implicit none
    real(precision_cbm5_strato_SOA), dimension(:), intent(inout) :: vSp
    real(precision_cbm5_strato_SOA), dimension(:), intent(in) :: photo
    type(Tchem_rules_cbm5_strato_SOA), intent(in) :: rules
    real, dimension(:), intent(in) :: metdat
    real, intent(in) :: timestep_sec
    real, dimension(:), intent(inout) :: garbage
    real, intent(in) :: zenith_cos, lat
    real, intent(inout) :: h_start
    type(silja_time), intent(in) :: now
    logical,  intent(out) :: print_it

    ! Local variables
    real :: tempr, sun, sunD, &
         & fPressure, fLandFr, &
         & fMixing_ratio_2_cnc, &
         & cH2O, cH2, cO2, cCH4, cCO2, fYear, cAir
    real :: SOA_b !With SOA
    integer :: indT, indZ, iSubst, iMode, nSteps, iStep, iHourTmp, jTmp, i
    real, dimension(num_species_cbm5_strato_SOA) :: garbage_local
    integer, save :: iCount = 0
    
    character(len=clen) :: preparation = 'preparation', integration = 'integration'
    ! Moles / m3 -> molecules / cm3
    !
    ! KPP integrator arguments
    !
    integer, dimension(20):: istatus
    real(precision_cbm5_strato_SOA), dimension(20) :: rstatus
    !real, dimension(num_species_cbm5_strato_SOA) :: atol, rtol

    integer :: ierr
    !real, dimension(num_react_cbm5_strato_SOA) :: rc
    !real, dimension(num_fixed_cbm5_strato_SOA) :: fix

    ! The initial size for the substep taken by the integrator. This is saved from
    ! the previous integration in this gridcell. Saving h_start gives about 50% cut
    ! in computation time for the cb4 transformation.
    if (h_start > 0.0) then
      kpp_rcntrl(3) = h_start
    else
      kpp_rcntrl(3) = 0.0
    end if
    print_it = .false.
    !--------------------------------------------------------------------------------
    !
    ! Get temperature from the meteo buffer
    !
    fPressure = metdat(ind_press)
    tempr = metdat(ind_tempr)
    
    !fLandFr = ptrLandFr(fu_grid_index(nx_meteo, ix, iy, interpCoefMet2DispHoriz, &
    !                                & ifInterpMet2DispHoriz))

    fMixing_ratio_2_cnc = fPressure / (gas_constant_uni * tempr)
    
    ! Fixed species' concentrations
    !
    ! Methane:
    !
    ! Old altitude variation that give somewhat worse fit than the new one below:
    ! It omits the latitude dependence and gives bigger values at high latitudes but underestimates
    ! the mixing ratio at 20-30km altitudes near the equator. 
    !if (fPressure > 0.3*std_pressure_sl) then ! CH4 in ppm, fit for vertical profile
    !  cCH4 = 1.
    !else
    !  cCH4 = 1. - 3.25e4 * (0.35 - fPressure/std_pressure_sl)**10   ! 1 -> 0.1 in stratosphere
    !end if
    ! Apply latitudinal distribution and refine absolute value. R.H. NOTE: THIS IS BUGGY AND OMITTED NOW!
    !cCH4 = cCH4 * (0.65 + 0.52 * (1.-fPressure/std_pressure_sl) * cos(lat * Pi / 180.0))
    !
    ! New altitude profile based on E.M. Buzan et al., Atmos. Meas. Tech. 9, 1095-1111 (2016).
    ! The surface mixing ratio extends higher in statosphere for latitudes near the equator.
    cCH4 = tanh((fPressure/std_pressure_sl)/(0.01+0.09*sin(Pi*lat/180.0)**2))**0.5;    
    !
    ! New methane temporal variation, 31July2018. Data from http://www.methanelevels.org/ (Years 1983->2017 are NOOA data):    
    fYear = fu_year(now)
    if (fYear > 2019.0) then !Old fit which MIGHT be a reasonable extrapolation in the future!?!? 
       cCH4 = cCH4 * (4.09e-5 / 46.1) * (0.83 + 1.5e-3*(fYear-1800) + 0.012 * (fYear-1930) / (1.+exp((1920-fYear)/50.)))
    elseif (fYear > 2006.0) then  !linear fit for years 2006->2017
       cCH4 = cCH4*(6.668006993*fYear-11603.725)*1e-9
    elseif (fYear > 1980.0) then  !This region between 1980->2016 is difficult to fit, but the following is a reasonable estimate.
       cCH4 = cCH4*(1491.64+283.0*sin((fYear-1980)*Pi/52.0)**0.42)*1e-9
    elseif (fYear > 1900.0) then  !fit for years 1900->1980.
       cCH4 = cCH4*exp(6.06038980e-07*fYear**3-3.48031329e-03*fYear**2+6.66681891*fYear-4253.08570)*1e-9  
    elseif (fYear > 1750.0) then  !fit for years 1750->1900.
       cCH4 = cCH4*exp(4.13844474e-08*fYear**3-2.20634143e-04*fYear**2+0.393026442*fYear-227.357323)*1e-9  
    else !Take the value 694 ppb at 1750 for earlier years. 
       cCH4 = cCH4*0.694e-6
    end if
    !Make some difference between North and South hemispheres (recent decades the max difference is about 140 ppb, about 8 percent):
    cCH4 = cCH4*(1.0+0.04*sin(Pi*lat/180.0))
    !Convert from mixing ratio -> molecules/cm3
    cCH4 = cCH4 * fPressure / (tempr * gas_constant_uni) !mixing ratio -> mole/m3
    cCH4 = cCH4 * mol2molec                              !mole/m3 -> molec/cm3 (needed for kpp)

    !Carbon dioxide, CO2:
    !The following fits are taken yearly averaged values from NASA and NOOA web-pages:
    !Volume mixing ratio in ppm
    if (fYear > 1958.5) then !Using the fit for years 1950 to 2017.
       !A reasonable extrapolation estimate up to 2050 if we compare the NASA 2 Degree C Scenario 
       cCO2 = 0.013462512*(fYear-1988.7312)*(fYear-1873.6165)+350.0
       cCO2 = min(cCO2,500.0) !Limit to 500 ppm
    elseif (fYear > 1939.525) then  !fit for years 1940 to 1960.
       cCO2 = 0.027739160*(fYear-1957.8590)*(fYear-1932.2721)+315.0
    elseif (fYear > 1850.0) then  !fit for years 1850->1940.
       cCO2=0.0015308339*(fYear-1910.8875)*(fYear-1681.4547)+300.0
    else !Take the value 284.29 ppm at 1850 for earlier years. 
       cCO2 = 284.29
    end if
    cCO2 = cCO2*1e-6 !since the above values are in ppm
    cCO2 = cCO2*mol2molec*fPressure / (tempr * gas_constant_uni)
    !END: Carbon dioxide, CO2:

    cAir = fPressure / (tempr * gas_constant_uni) * mol2molec
    cH2 = 5.0e-07 * cAir
    cO2 = 0.2095 * cAir
    cH2O = metdat(ind_spec_hum) * fMixing_ratio_2_cnc * molecular_weight_air / molecular_weight_water
    cH2O = cH2O * mol2molec ! to molec/cm3
    
    !Next are specific to SOA only:
    ! Branching ratio for SOA formation reactions for low and hi NOx cases
    if (vSp(iXO2) .eps. 0.0) then
      SOA_b = 1.0 
    else
      SOA_b = vSp(iNO)* 2.6* EXP(365.0/tempr) / &
            & (vSp(iNO)* 2.6* EXP(365.0/tempr) +  &
            &  vSp(iXO2)*6.8E-2 + &
            &  vSp(iHO2)*7.737e-2 * EXP(-1300.0/tempr))
    endif
    !END: Next are specific to SOA only.

    call set_rates_cbm5_strato_SOA(photo, tempr, cAir, fPressure, cH2O&
                             &, SOA_b, rules%ASOA_aging, rules%BSOA_aging, rules%IVOC_aging, rules%ifMTprods&
                             &, rules%MCF_factor&
                             &)

! Check which argumets are needed, use named arguments to avoid mismatches. For cbm4_SOA remove H2 and CO2 species.
    call set_fixed_cbm5_strato_SOA(cnc_h2o=cH2O, cnc_o2=cO2, cnc_m=cAir, cnc_ch4=cCH4, cnc_h2=cH2, cnc_co2=cCO2) !NOTE: Arguments specific to cbm42_strato/strato_SOA!  
    !call set_fixed_cbm5_strato_SOA(cnc_h2o=cH2O, cnc_o2=cO2, cnc_m=cAir, cnc_ch4=cCH4, cnc_h2=cH2) !NOTE: Arguments specific to cbm5_SOA!  
    !call set_fixed_cbm4_SOA(cnc_h2o=cH2O, cnc_o2=cO2, cnc_m=cAir, cnc_ch4=cCH4) !NOTE: Arguments specific to cbm4_SOA! 

    vSp(:) = vSp(:) * mol2molec
!    garbage(1:num_species_cbm5_strato_SOA) = garbage(1:num_species_cbm5_strato_SOA)*mol2molec
    garbage_local = 0.0

    call rosenbrock(num_species_cbm5_strato_SOA, vSp(1:num_species_cbm5_strato_SOA), &
                  & 0.0_precision_cbm5_strato_SOA, real(timestep_sec, precision_cbm5_strato_SOA),&
                  & rules%kpp_atol, rules%kpp_rtol, garbage_local, &
!                  & rules%kpp_atol, rules%kpp_rtol, garbage(1:num_species_cbm5_strato_SOA), &
		  & kpp_rcntrl, rules%kpp_icntrl, &
                  & rstatus, istatus, ierr) ! , rc, fix)
      

    ! The negativity limits as of 22.6.2009:
    ! < -1e-7 -> error
    ! < -1e-8 -> warning, add to garbage
    ! < -1e-9 -> to garbage
    ! < -1e-11 -> neglect

    vSp(:) = vSp(:) / mol2molec
    garbage_local(:) = garbage_local(:) / mol2molec

    do i = 1, num_species_cbm5_strato_SOA

!      if(garbage_local(i) < 0.0)then
!        iCount = iCount + 1
!        if(mod(iCount,100) == 0) call msg('After Rosenbrock garbage_local <0. Species:' + &
!                                        & fu_str(species_cbm5_strato_SOA(i)) + ', count, value =' + &
!                                        & fu_str(iCount), i, garbage_local(i))
!        !  garbage_local(i) = 0.
!      endif
      
!!$ if (isnan(vsp(i))) call msg('nan in cocktail', omp_get_thread_num())
      if (vSp(i) < 0.) then
        if (vSp(i) > -1e-11) then
          call msg('Concentration between -1e-11 and 0.0 after chemistry'+ &
                 & fu_str(species_cbm5_strato_SOA(i)),vSp(i))
        else if (vSp(i) > -1e-9) then
          call msg('Concentration between -1e-9 and -1e-11 after chemistry'+ &
                 & fu_str(species_cbm5_strato_SOA(i)),vSp(i))
        else if (vSp(i) > -1e-8) then
          call msg_warning('Concentration between -1e-8 and -1e-9 after chemistry, species:' + &
                         & fu_str(species_cbm5_strato_SOA(i)) + '_' + fu_str(vSp(i)), &
                         & 'transform_cbm5_strato_SOA')
          !do jTmp = 1, num_species_cbm5_strato_SOA
          !  call msg('Concentration for:' + fu_str(species_cbm5_strato_SOA(jTmp)), vSp(jTmp))
          !enddo
        else
          call msg_warning('Something somewhere went terribly wrong...(iSp,vSp)'+ &
                         & fu_str(species_cbm5_strato_SOA(i)) + '_' + fu_str(vSp(i)), &
                         & 'transform_cbm5_strato_SOA')
        end if
        print_it = .true.
        garbage_local(i) = garbage_local(i) + vSp(i)
        vSp(i) = 0.0
        ! Otherwise leave as is - will be trapped later.
      end if
    end do

    !call stop_count(chCounterNm = integration)
    garbage(1:num_species_cbm5_strato_SOA) = garbage(1:num_species_cbm5_strato_SOA) + &
                                        & garbage_local(1:num_species_cbm5_strato_SOA)

    if (ierr < 0)then
      call msg_warning('The integrator had an error:' + fu_str(ierr), 'transf_cbm5_strato_SOA_single_val')
      !do i = 1, num_species_cbm5_strato_SOA
      !  call msg('Concentration for:' + fu_str(species_cbm5_strato_SOA(i)), vSp(i))
      !enddo
      call msg('h_start on entry:', h_start)
      call set_error('The integrator had an error', 'transf_cbm5_strato_SOA_single_val')
      print_it = .true.
      
    endif
    h_start = rstatus(3)
    
  end subroutine transform_cbm5_strato_SOA
  
  !********************************************************************************

  subroutine transform_cbm5_strato_SOA_adj(vSp, vSp_lin, photo, rules, metdat, timestep_sec, garbage, &
                                   & zenith_cos, lat, h_start, now, print_it)

    implicit none
    real(precision_cbm5_strato_SOA), dimension(:), intent(inout) :: vSp
    real(precision_cbm5_strato_SOA), dimension(:), intent(inout) :: vSp_lin ! linearization point
    real(precision_cbm5_strato_SOA), dimension(:), intent(in) :: photo
    type(Tchem_rules_cbm5_strato_SOA), intent(in) :: rules
    real, dimension(:), intent(in) :: metdat
    real, intent(in) :: timestep_sec
    real, dimension(:), intent(inout) :: garbage
    real, intent(in) :: zenith_cos, lat
    real, intent(inout) :: h_start
    type(silja_time), intent(in) :: now
    logical, intent(out) :: print_it

    ! Local variables
    real :: tempr, sun, sunD, &
         & fPressure, fLandFr, &
         & fMixing_ratio_2_cnc, &
         & cH2O, cH2, cO2, cCH4, cCO2, fYear, cAir
    integer :: indT, indZ, iSubst, iMode, nSteps, iStep, iHourTmp, jTmp, i
    
    character(len=clen) :: preparation = 'preparation', integration = 'integration'
    ! Moles / m3 -> molecules / cm3
    !
    ! KPP integrator arguments
    !
    integer, dimension(20):: istatus
    real(precision_cbm5_strato_SOA), dimension(20) :: rstatus
    !real, dimension(num_species_cbm5_strato_SOA) :: atol, rtol
    !real, dimension(32) :: cTmp
    integer :: ierr
    !real, dimension(num_react_cbm5_strato_SOA) :: rc
    !real, dimension(num_fixed_cbm5_strato_SOA) :: fix
    
    ! The initial size for the substep taken by the integrator. This is saved from
    ! the previous integration in this gridcell. Saving h_start gives about 50% cut
    ! in computation time for the cb4 transformation.
    if (h_start > 0.0) then
      kpp_rcntrl(3) = h_start
    else
      kpp_rcntrl(3) = 0.0
    end if
    print_it = .false.
    !----------------------------------	----------------------------------------------
    !
    ! Get temperature from the meteo buffer
    !
    fPressure = metdat(ind_press)
    tempr = metdat(ind_tempr)
    
    !fLandFr = ptrLandFr(fu_grid_index(nx_meteo, ix, iy, interpCoefMet2DispHoriz, &
    !                                & ifInterpMet2DispHoriz))

    fMixing_ratio_2_cnc = fPressure / (gas_constant_uni * tempr)
    
    ! Fixed species' concentrations
    !
    ! Methane:
    !
    ! New altitude profile based on E.M. Buzan et al., Atmos. Meas. Tech. 9, 1095-1111 (2016).
    ! The surface mixing ratio extends higher in statosphere for latitudes near the equator.
    cCH4 = tanh((fPressure/std_pressure_sl)/(0.01+0.09*sin(Pi*lat/180.0)**2))**0.5;    
    !
    ! New methane temporal variation, 31July2018. Data from http://www.methanelevels.org/ (Years 1983->2017 are NOOA data):    
    fYear = fu_year(now)
    if (fYear > 2019.0) then !Old fit which MIGHT be a reasonable extrapolation in the future!?!? 
       cCH4 = cCH4 * (4.09e-5 / 46.1) * (0.83 + 1.5e-3*(fYear-1800) + 0.012 * (fYear-1930) / (1.+exp((1920-fYear)/50.)))
    elseif (fYear > 2006.0) then  !linear fit for years 2006->2017
       cCH4 = cCH4*(6.668006993*fYear-11603.725)*1e-9
    elseif (fYear > 1980.0) then  !This region between 1980->2016 is difficult to fit, but the following is a reasonable estimate.
       cCH4 = cCH4*(1491.64+283.0*sin((fYear-1980)*Pi/52.0)**0.42)*1e-9
    elseif (fYear > 1900.0) then  !fit for years 1900->1980.
       cCH4 = cCH4*exp(6.06038980e-07*fYear**3-3.48031329e-03*fYear**2+6.66681891*fYear-4253.08570)*1e-9  
    elseif (fYear > 1750.0) then  !fit for years 1750->1900.
       cCH4 = cCH4*exp(4.13844474e-08*fYear**3-2.20634143e-04*fYear**2+0.393026442*fYear-227.357323)*1e-9  
    else !Take the value 694 ppb at 1750 for earlier years. 
       cCH4 = cCH4*0.694e-6
    end if
    !Make some difference between North and South hemispheres (recent decades the max difference is about 140 ppb, about 8 percent):
    cCH4 = cCH4*(1.0+0.04*sin(Pi*lat/180.0))
    !Convert from mixing ratio -> molecules/cm3
    cCH4 = cCH4 * fPressure / (tempr * gas_constant_uni) !mixing ratio -> mole/m3
    cCH4 = cCH4 * mol2molec                              !mole/m3 -> molec/cm3 (needed for kpp)

    !Carbon dioxide, CO2:
    !The following fits are taken yearly averaged values from NASA and NOOA web-pages:
    !Volume mixing ratio in ppm
    if (fYear > 1958.5) then !Using the fit for years 1950 to 2017.
       !A reasonable extrapolation estimate up to 2050 if we compare the NASA 2 Degree C Scenario 
       cCO2 = 0.013462512*(fYear-1988.7312)*(fYear-1873.6165)+350.0
       cCO2 = min(cCO2,500.0) !Limit to 500 ppm
    elseif (fYear > 1939.525) then  !fit for years 1940 to 1960.
       cCO2 = 0.027739160*(fYear-1957.8590)*(fYear-1932.2721)+315.0
    elseif (fYear > 1850.0) then  !fit for years 1850->1940.
       cCO2=0.0015308339*(fYear-1910.8875)*(fYear-1681.4547)+300.0
    else !Take the value 284.29 ppm at 1850 for earlier years. 
       cCO2 = 284.29
    end if
    cCO2 = cCO2*1e-6 !since the above values are in ppm
    cCO2 = cCO2*mol2molec*fPressure / (tempr * gas_constant_uni)
    !END: Carbon dioxide, CO2:

    cAir = fPressure / (tempr * gas_constant_uni) * mol2molec
    cH2 = 5.0e-07 * cAir
    cO2 = 0.2095 * cAir
    cH2O = metdat(ind_spec_hum) * fMixing_ratio_2_cnc * molecular_weight_air / molecular_weight_water
    cH2O = cH2O * mol2molec ! to molec/cm3

    call set_rates_cbm5_strato_SOA(photo, tempr, cAir, fPressure, cH2O&
                             &, soa_b=0.0, ASOA_aging=rules%ASOA_aging, BSOA_aging=rules%BSOA_aging&
                             &, IVOC_aging=rules%IVOC_aging, ifMTprods=rules%ifMTprods&
                             &, MCF_factor=rules%MCF_factor&
                             &)
    ! Check which arguments are needed, use named arguments to avoid mismatches. For cbm4_SOA remove H2 and CO2 species.
    call set_fixed_cbm5_strato_SOA(cnc_h2o=cH2O, cnc_o2=cO2, cnc_m=cAir, cnc_ch4=cCH4, cnc_h2=cH2, cnc_co2=cCO2) !NOTE: Arguments specific to cbm42_strato/strato_SOA!
    !call set_fixed_cbm5_strato_SOA(cnc_h2o=cH2O, cnc_o2=cO2, cnc_m=cAir, cnc_ch4=cCH4, cnc_h2=cH2) !NOTE: Arguments specific to cbm5_SOA!
    !call set_fixed_cbm4_SOA(cnc_h2o=cH2O, cnc_o2=cO2, cnc_m=cAir, cnc_ch4=cCH4) !NOTE: Arguments specific to cbm4_SOA!
    
    vSp_lin(:) = vSp_lin(:) * mol2molec
    ! One could consider scaling also the perturbations for better numerics, but there
    ! seems to be no practical benefit.

    !call set_error('Adjoint not available', 'transform_cbm5_strato_SOA')
    call rosenbrockAdj(vSp_lin(1:num_species_cbm5_strato_SOA), 1, & 
                     & vSp(1:num_species_cbm5_strato_SOA), &
                     & 0.0_precision_cbm5_strato_SOA, -real(timestep_sec, precision_cbm5_strato_SOA), &
                     & rules%kpp_atol, rules%kpp_rtol, & 
                     & rules%kpp_atol, rules%kpp_rtol, kpp_rcntrl, &
                     & rules%kpp_icntrl, rstatus, istatus, ierr) ! , rc, fix)
    vSp_lin(:) = vSp_lin(:) / mol2molec
    h_start = rstatus(3)
    if (ierr < 0) then
      call msg('ierr', ierr)
      call set_error('The integrator had an error', 'transf_cbm5_strato_SOA_adj')
      print_it = .true.
    end if
    
  end subroutine transform_cbm5_strato_SOA_adj
  
  !************************************************************************************

  ! Functions for correcting photolysis rates based on cloud cover, from RADM (Chang et
  ! al., 1987). The functions require a value for the cloud transmissivity; I have used
  ! 0.5 as representative for mid-level clouds; Savijarvi's lecture notes give a range of
  ! 0.2 (covective or low-level clouds) to 0.8 (high-level clouds). The transmissivity
  ! could be evaluated from LWP; a formula supposedly used in some HIRLAM version is
  ! T = 40 * (0.5 + m) / (LWP + 40 * (0.5+m)) 
  ! where m is zenith angle cosine and LWP is in g/m2. The Chang paper gives another formula.
  ! 
  ! Different functions are for below, within and above clouds.  To use, consider for
  ! example: 
  ! call set_rates_cbm4(tempr, fPressure, sun, fu_cloudcorr_below(fCloudCover, 0.5, sun))

  real function fu_cloudcorr_below(cloud_fract, cloud_trans, zenith_cos) result(jscale)
    implicit none
    real, intent(in) :: cloud_fract
    real, intent(in) :: cloud_trans ! cloud transmittance, 0...1
    real, intent(in) :: zenith_cos
    
    real :: fcld

    fcld = 1.6 * cloud_trans * zenith_cos
    jscale = 1.0 + cloud_fract * (fcld - 1.0)

  end function fu_cloudcorr_below

  real function fu_cloudcorr_incld(cloud_fract, cloud_trans, zenith_cos) result(jscale)
    implicit none
    real, intent(in) :: cloud_fract
    real, intent(in) :: cloud_trans ! cloud transmittance, 0...1
    real, intent(in) :: zenith_cos
    
    real :: fcld

    fcld = 1.4 * zenith_cos
    jscale = 1.0 + cloud_fract * (fcld - 1.0)

  end function fu_cloudcorr_incld

  real function fu_cloudcorr_above(cloud_fract, cloud_trans, zenith_cos) result(jscale)
    implicit none
    real, intent(in) :: cloud_fract
    real, intent(in) :: cloud_trans ! cloud transmittance, 0...1
    real, intent(in) :: zenith_cos
    
    real :: fcld
    real, parameter :: alpha = 0.7 ! for ozone

    fcld = 1.0 + alpha*(1.0 - cloud_trans) * zenith_cos
    jscale = 1.0 + cloud_fract * (fcld - 1.0)

  end function fu_cloudcorr_above


  !********************************************************************************

  
  function fu_lifetime_cbm5_strato_SOA(rules) result(lifetime)
    !
    ! A typical life time due to degradation and deposition for SOx materials
    ! So far, whatever the species are, we set the "typical" lifetime to be two days
    !
    implicit none

    type(silja_interval) :: lifetime
    type(Tchem_rules_cbm5_strato_SOA), intent(in) :: rules

    lifetime = one_day * 2.0

  end function fu_lifetime_cbm5_strato_SOA

  logical function fu_if_specific_deposition_cbm5_strato_SOA(rules_cbm5_strato_SOA) result(if_specific)
    implicit none
    type(Tchem_rules_cbm5_strato_SOA), intent(in) :: rules_cbm5_strato_SOA
    if_specific = .false.     ! no specific deposition whatsoever
  end function fu_if_specific_deposition_cbm5_strato_SOA

  !************************************************************************************

  !The following subroutine is not needed after PAR-OLE tweaks are all implemented.
  !This routine is called from chemistry_manager.mod.f90 where the calls are disabled
  !for all the chemistries where the PAR-OLE tweaks are applied.
  subroutine check_low_mass_threshold_cbm5_strato_SOA(rules_cbm5_strato_SOA, speciesTrn, nSpeciesTrn, ifForwardRun, &
                                        & low_cnc_threshold, low_mass_threshold)
    !
    ! This subroutine ensures that, in case of strong chemical connection, one substance can never
    ! be dropped to garbage while the other one still stays alive. In-essence, for cb4 we want 
    ! PAR to be always present if OLE are in
    !
    implicit none
    
    ! Imported parameters
    type(Tchem_rules_cbm5_strato_SOA), intent(in) :: rules_cbm5_strato_SOA
    logical, intent(in) :: ifForwardRun
    type(silam_species), dimension(:), pointer :: speciesTrn
    real, dimension(:), pointer :: low_cnc_threshold, low_mass_threshold
    integer, intent(in) :: nSpeciesTrn
    
    ! Local variables
    type(silam_species) :: spOLE, spPAR
    integer :: indOLE, indPAR
    !
    ! Find the species
    !
    call set_species(spOLE, fu_get_material_ptr('OLE '), in_gas_phase)
    call set_species(spPAR, fu_get_material_ptr('PAR '), in_gas_phase)
    indOLE = fu_index(spOLE, speciesTrn, nSpeciesTrn)
    indPAR = fu_index(spPAR, speciesTrn, nSpeciesTrn)
    if(fu_fails(indOLE /= int_missing,'No OLE in CB4 list','check_low_mass_threshold_cbm5_strato_SOA'))return
    if(fu_fails(indPAR /= int_missing,'No PAR in CB4 list','check_low_mass_threshold_cbm5_strato_SOA'))return
    !
    ! Check the relation
    !
    if(low_cnc_threshold(indPAR) > low_cnc_threshold(indOLE))then
      call msg('cbm5_strato_SOA: bad relation between OLE and PAR thresholds. Take OLE (the second):', &
                                            & low_cnc_threshold(indPAR), low_cnc_threshold(indOLE))
      low_cnc_threshold(indPAR) = low_cnc_threshold(indOLE)
      low_mass_threshold(indPAR) = low_mass_threshold(indOLE)
    endif
    
  end subroutine check_low_mass_threshold_cbm5_strato_SOA
  
  
  

  !************************************************************************************

  logical function fu_if_tla_required_cbm5_strato_SOA(rules) result(required)
    implicit none
    type(Tchem_rules_cbm5_strato_SOA), intent(in) :: rules
    
    ! call set_error('TLA not available for cbm5_strato_SOA', 'fu_if_tla_required_cbm5_strato_SOA')
     required = .true.
    
  end function fu_if_tla_required_cbm5_strato_SOA


END MODULE chem_dep_cbm5_strato_SOA

