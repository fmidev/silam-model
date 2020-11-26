module chem_dep_cbm4
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
  use cbm4_interface
  use cocktail_basic
  !$use omp_lib

  implicit none
  private

  public set_chem_rules_cbm4
  public set_missing
  
  public fu_lifetime
  public check_low_mass_threshold_cb4
  public registerSpecies
  public transform_cbm4
  public transform_cbm4_adj

  public cbm4_input_needs
  public fu_if_specific_deposition
  public prepare_step_cbm4
  
  ! inherited from the interface
  public init_chemicals_cbm4
  public register_reaction_rates_cbm4

  !
  ! Private routines of the sulphur DMAT cocktail and chemical rules
  !
  private set_miss_chem_rules_cbm4

  private fu_lifetime_cbm4
  private registerSpeciescbm4
  public  inventory_cbm4
  private fu_if_specific_deposition_cbm4

  interface registerSpecies
     module procedure registerSpeciescbm4
  end interface

  interface set_missing
    module procedure set_miss_chem_rules_cbm4
  end interface

  interface fu_lifetime
    module procedure fu_lifetime_cbm4
  end interface

  interface fu_if_specific_deposition
    module procedure fu_if_specific_deposition_cbm4
  end interface

  public fu_if_tla_required
  private fu_if_tla_required_cbm4
  interface fu_if_tla_required
     module procedure fu_if_tla_required_cbm4
  end interface


  type Tchem_rules_cbm4
     private
     logical :: defined = .false.
     integer, dimension(20) :: kpp_icntrl
     ! solver relative tolerance
     real(precision_cbm4), dimension(num_species_cbm4) :: kpp_rtol
     ! solver absolute tolerance, solver units
     real(precision_cbm4), dimension(num_species_cbm4) :: kpp_atol
     ! the rtol doesn't work here: it is updated per cell
     !real, dimension(num_species_cbm4) :: kpp_rcntrl
     integer :: react_rates_idx0, react_rates_idx1 !Indices of our reactions in reaction rates massmap 
  end type Tchem_rules_cbm4
  
  public Tchem_rules_cbm4


  !
  ! The mark of this transformation module:
  !
  integer, parameter, public :: transformation_cbm4 = 5009
  ! an argument to the solver. Set initially except for the parameter which defined
  ! initial internal step: this is set for each cell according to the previous value.
  real(precision_cbm4), dimension(20), private, save :: kpp_rcntrl
  !$OMP THREADPRIVATE(kpp_rcntrl)

  real, parameter, private :: mol2molec = avogadro*1e-6

  integer, pointer :: ind_tempr, ind_press, ind_cloud_covr, ind_spec_hum


CONTAINS



  !************************************************************************************

  subroutine prepare_step_cbm4()
    implicit none

    call set_const_rates_cbm4()
    kpp_rcntrl = 0.0
    kpp_rcntrl(1) = 1.0e-3
  end subroutine prepare_step_cbm4


  !*************************************************************************************

  subroutine inventory_cbm4(rules, &
                               & speciesEmis, speciesTransp, speciesShortlived, speciesAerosol,&
                               & nSpeciesEmis, nSpeciesTransp, nSpeciesShortlived, nSpeciesAerosol, &
                               & iClaimedSpecies)
    implicit none
    !
    ! The list of species + claiming the species. Call the routine from interface, but
    ! without rules which is not known there.
    type(Tchem_rules_cbm4), intent(in) :: rules
    type(silam_species), dimension(:), pointer :: speciesEmis, speciesTransp, speciesShortlived, &
                                                & speciesAerosol
    integer, intent(in) :: nSpeciesEmis
    integer, intent(inout) :: nspeciesTransp, nspeciesShortlived, nspeciesAerosol
    integer, dimension(:), intent(inout) :: iClaimedSpecies

    ! Local variables
    integer :: iEmis

    call inventory_int(transformation_cbm4, speciesEmis, &
                     & speciesTransp, speciesShortlived, speciesAerosol,&
                     & nSpeciesEmis, nSpeciesTransp, nSpeciesShortlived, nSpeciesAerosol, &
                     & iClaimedSpecies)

  end subroutine inventory_cbm4


  !************************************************************************************

  subroutine registerSpeciescbm4(rules, &
                                    & speciesTransp, speciesShortlived, speciesAerosol,&
                                    & nspeciesTransp, nspeciesShortlived, nspeciesAerosol)
    !
    ! The integrator requires fixed indices for the species. At this
    ! point, we just check that the transport species follow this ordering.
    implicit none
    type(Tchem_rules_cbm4), intent(in) :: rules
    type(silam_species), dimension(:), pointer :: speciesTransp, speciesShortlived, &
                                                & speciesAerosol
    integer, intent(in) :: nspeciesTransp, nspeciesShortlived, nspeciesAerosol
    
    ! Local variable
    integer :: i

    call register_species_int(speciesTransp, speciesShortlived, speciesAerosol,&
                            & nspeciesTransp, nspeciesShortlived, nspeciesAerosol)
    
  end subroutine registerSpeciescbm4

  subroutine  register_reaction_rates_cbm4(rules,nReactRates)

    implicit none
    type(Tchem_rules_cbm4), intent(inout) :: rules
    integer, intent (inout) :: nReactRates

    rules%react_rates_idx0 = nReactRates + 1 
    nReactRates = nReactRates + num_react_cbm4
    rules%react_rates_idx1 = nReactRates
    call msg ("CBM4 reaction rates from, to", rules%react_rates_idx0, rules%react_rates_idx1 ) 

  end subroutine  register_reaction_rates_cbm4


  !*****************************************************************

  subroutine cbm4_input_needs(rules, metdat)
    !
    ! Returns input needs for the transformation routine. Change if needed.  
    ! 
    implicit none

    ! Imported parameters
    type(Tchem_rules_cbm4), intent(in) :: rules
    type(Tmeteo_input), intent(out), target :: metdat


    if (.not. rules%defined) then
      call set_error('Undefined cbm4 rules', 'cbm4_input_needs')
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


  end subroutine cbm4_input_needs

  !***********************************************************************

  subroutine set_miss_chem_rules_cbm4(rules)
    implicit none
    type(Tchem_rules_cbm4), intent(out) :: rules

    rules%defined = .false.

  end subroutine set_miss_chem_rules_cbm4
  
  !***********************************************************************

  subroutine set_chem_rules_cbm4(nlSetup, rules)
    !
    ! Creates the set of rules, in accordance with which the computations are to
    ! be made. Input - SILAM namelist
    ! Taking this chance, let's also set once-and-forever a few chemical coefficients
    !
    implicit none

    ! Imported parameter
    type(Tsilam_namelist), intent(in) :: nlSetup 
    type(Tchem_rules_cbm4) :: rules

    ! Local variables
    integer :: i, jTmp, ierr
    real :: fTempr, z, ratio_K_prod, ratio_K_loss, f300_T, fRelatDensity, fRelatPress, fRelatTempr

    !
    ! Stupidity checking
    !
    rules%defined = .false.
    if(.not.defined(nlSetup))then
      call set_error('Undefined namelist given','set_chem_rules_cbm4')
      return
    endif
    
    rules%kpp_icntrl = 0 ! all default, except:
    rules%kpp_icntrl(1) = 1 ! autonomous system

    ! The error tolerances. Inside the KPP intergrator, the error is computed something
    ! like err = |y_i+1 - y_i| / scale, where scale = atol + rtol*y.  We can set rtol as
    ! high as we tolerate, but it seems that atol has better be small in order to avoid
    ! negative concentrations.
    rules%kpp_atol = 1e-15*mol2molec
    !rules%kpp_atol = 1.0
    rules%kpp_rtol = 0.01

    call init_solver(rules%kpp_icntrl, ierr)
    if (fu_fails(ierr == 0, 'Error initializing solver', 'set_chem_rules_cbm4')) return
    call init_solver_adj(rules%kpp_icntrl, ierr)
    if (fu_fails(ierr == 0, 'Error initializing adjoint solver', 'set_chem_rules_cbm4')) return
    rules%defined = .true.

  end subroutine set_chem_rules_cbm4


  !***********************************************************************
  subroutine transform_cbm4(vSp, rules, metdat, timestep_sec, garbage, &
                               & zenith_cos, h_start, print_it, reactrates)

    implicit none
    real(precision_cbm4), dimension(:), intent(inout) :: vSp
    type(Tchem_rules_cbm4), intent(in) :: rules
    real, dimension(:), intent(in) :: metdat
    real, intent(in) :: timestep_sec
    real, dimension(:), intent(inout) :: garbage
    real, intent(in) :: zenith_cos
    real, intent(inout) :: h_start
    logical, intent(out) :: print_it
    real, dimension(:), optional, intent(inout) :: reactrates

    ! Local variables
    real :: tempr, sun, sunD, &
         & fPressure, fCloudCover, fLandFr, &
         & fMixing_ratio_2_cnc, &
         & cH2O
    integer :: indT, indZ, iSubst, iMode, nSteps, iStep, iHourTmp, jTmp, i
    
    character(len=clen) :: preparation = 'preparation', integration = 'integration'
    ! Moles / m3 -> molecules / cm3
    !
    ! KPP integrator arguments
    !
    integer, dimension(20):: istatus
    real(precision_cbm4), dimension(20) :: rstatus
    !real, dimension(num_species_cbm4) :: atol, rtol
    !real, dimension(32) :: cTmp
    integer :: ierr
    !real, dimension(num_react_cbm4) :: rc
    !real, dimension(num_fixed_cbm4) :: fix
    
    print_it = .false.  ! set to true and chemistry manager will give complete dump for this cell

    ! The initial size for the substep taken by the integrator. This is saved from
    ! the previous integration in this gridcell. Saving h_start gives about 50% cut
    ! in computation time for the cb4 transformation.
    if (h_start > 0.0) then
      kpp_rcntrl(3) = h_start
    else
      kpp_rcntrl(3) = 0.0
    end if
    
    !--------------------------------------------------------------------------------
    !
    ! Get temperature from the meteo buffer
    !
    fPressure = metdat(ind_press)
    tempr = metdat(ind_tempr)
    
    !fLandFr = ptrLandFr(fu_grid_index(nx_meteo, ix, iy, interpCoefMet2DispHoriz, &
    !                                & ifInterpMet2DispHoriz))

    fMixing_ratio_2_cnc = fPressure / (gas_constant_uni * tempr)

    cH2O = metdat(ind_spec_hum) * fMixing_ratio_2_cnc * molecular_weight_air / molecular_weight_water
    cH2O = cH2O * mol2molec ! to molec/cm3
    
    fCloudCover = metdat(ind_cloud_covr)
   
    sun = zenith_cos
    if(sun < 0.) sun = 0.  ! night
    
    call set_rates_cbm4(tempr, fPressure, sun, 1.0 - 0.5*fCloudCover**1.0)
    call set_fixed_cbm4(cnc_h2o=cH2O) ! check which argumets are needed!

    vSp(:) = vSp(:) * mol2molec

    call rosenbrock(num_species_cbm4, vSp(1:num_species_cbm4), &
                  & 0.0_precision_cbm4, real(timestep_sec, precision_cbm4),&
                  & rules%kpp_atol, rules%kpp_rtol, kpp_rcntrl, rules%kpp_icntrl, &
                  & rstatus, istatus, ierr) ! , rc, fix)
      
    h_start = rstatus(3)

    ! The negativity limits as of 22.6.2009:
    ! < -1e-7 -> error
    ! < -1e-8 -> warning, add to garbage
    ! < -1e-9 -> to garbage
    ! < -1e-11 -> neglect

    do i = 1, num_species_cbm4
!!$ if (isnan(vsp(i))) call msg('nan in cocktail', omp_get_thread_num())
      if (vSp(i) < 0.) then
        if (vSp(i) > -1e-11 * mol2molec) then
          !call msg_warning('Concentration between -1e-12 and 0.0 after chemistry')
          vSp(i) = 0.0
        else if (vSp(i) > -1e-9 * mol2molec) then
          garbage(i) = garbage(i) + vSp(i) / mol2molec
          vSp(i) = 0.0
        else if (vSp(i) > -1e-8 * mol2molec) then
          call msg_warning('Concentration between -1e-8 and -1e-11 after chemistry, species:' + &
               & fu_str(species_cbm4(i)))
          do jTmp = 1, num_species_cbm4
            call msg('Concentration for:' + fu_str(species_cbm4(jTmp)), &
                 & vSp(jTmp) / mol2molec)
          enddo
          garbage(i) = garbage(i) + vSp(i) / mol2molec
          vSp(i) = 0.0
        end if
        ! Otherwise leave as is - will be trapped later.
      end if
    end do

    !call stop_count(chCounterNm = integration)

    vSp(:) = vSp(:) / mol2molec

    if (ierr < 0)then
      call msg_warning('The integrator had an error:' + fu_str(ierr), 'transf_cbm4_single_val')
      do i = 1, num_species_cbm4
        call msg('Concentration for:' + fu_str(species_cbm4(i)), &
             & vSp(i))
      enddo
      call set_error('The integrator had an error', 'transf_cbm4_single_val')
    endif

    ! Dump reaction rates in units of mole/m3 sec
    if (present(reactrates)) then
      reactrates(rules%react_rates_idx0:rules%react_rates_idx1) = a(:) / mol2molec
    endif

    
  end subroutine transform_cbm4
  
  !********************************************************************************

  subroutine transform_cbm4_adj(vSp, vSp_lin, rules, metdat, timestep_sec, garbage, &
                                   & zenith_cos, h_start, print_it)

    implicit none
    real(precision_cbm4), dimension(:), intent(inout) :: vSp
    real(precision_cbm4), dimension(:), intent(inout) :: vSp_lin ! linearization point
    type(Tchem_rules_cbm4), intent(in) :: rules
    real, dimension(:), intent(in) :: metdat
    real, intent(in) :: timestep_sec
    real, dimension(:), intent(inout) :: garbage
    real, intent(in) :: zenith_cos
    real, intent(inout) :: h_start
    logical, intent(out) :: print_it

    ! Local variables
    real :: tempr, sun, sunD, &
         & fPressure, fCloudCover, fLandFr, &
         & fMixing_ratio_2_cnc, &
         & cH2O
    integer :: indT, indZ, iSubst, iMode, nSteps, iStep, iHourTmp, jTmp, i
    
    character(len=clen) :: preparation = 'preparation', integration = 'integration'
    ! Moles / m3 -> molecules / cm3
    !
    ! KPP integrator arguments
    !
    integer, dimension(20):: istatus
    real(precision_cbm4), dimension(20) :: rstatus
    !real, dimension(num_species_cbm4) :: atol, rtol
    !real, dimension(32) :: cTmp
    integer :: ierr
    !real, dimension(num_react_cbm4) :: rc
    !real, dimension(num_fixed_cbm4) :: fix
    
    print_it = .false.  ! set to true and chemistry manager will give complete dump for this cell

    ! The initial size for the substep taken by the integrator. This is saved from
    ! the previous integration in this gridcell. Saving h_start gives about 50% cut
    ! in computation time for the cb4 transformation.
    if (h_start > 0.0) then
      kpp_rcntrl(3) = h_start
    else
      kpp_rcntrl(3) = 0.0
    end if
    
    !--------------------------------------------------------------------------------
    !
    ! Get temperature from the meteo buffer
    !
    fPressure = metdat(ind_press)
    tempr = metdat(ind_tempr)
    
    !fLandFr = ptrLandFr(fu_grid_index(nx_meteo, ix, iy, interpCoefMet2DispHoriz, &
    !                                & ifInterpMet2DispHoriz))

    fMixing_ratio_2_cnc = fPressure / (gas_constant_uni * tempr)

    cH2O = metdat(ind_spec_hum) * fMixing_ratio_2_cnc * molecular_weight_air / molecular_weight_water
    cH2O = cH2O * mol2molec ! to molec/cm3
    
    fCloudCover = metdat(ind_cloud_covr)
   
    sun = zenith_cos
    if(sun < 0.) sun = 0.  ! night
    
    call set_rates_cbm4(tempr, fPressure, sun, 1.0 - 0.5*fCloudCover**1.0)
    call set_fixed_cbm4(cnc_h2o=cH2O) ! check which argumets are needed!

    vSp_lin(:) = vSp_lin(:) * mol2molec

    ! One could consider scaling the perturbations as well for better numerics, but there
    ! seems to be no practical benefit.
    vSp(:) = vSp(:) * mol2molec 
    !'transform_cbm4')
    call rosenbrockAdj(vSp_lin(1:num_species_cbm4), 1, & 
                     & vSp(1:num_species_cbm4), &
                     & 0.0_precision_cbm4, -real(timestep_sec, precision_cbm4), &
                     & rules%kpp_atol, rules%kpp_rtol, & 
                     & rules%kpp_atol, rules%kpp_rtol, kpp_rcntrl, &
                     & rules%kpp_icntrl, rstatus, istatus, ierr) ! , rc, fix)
    vSp_lin(:) = vSp_lin(:) / mol2molec
    vSp(:) = vSp(:) / mol2molec

    h_start = rstatus(3)
    if (ierr < 0) then
      call msg('ierr', ierr)
      call set_error('The integrator had an error', 'transf_cbm4_adj')
    end if
    
  end subroutine transform_cbm4_adj
  
  !********************************************************************************

  
  function fu_lifetime_cbm4(rules) result(lifetime)
    !
    ! A typical life time due to degradation and deposition for SOx materials
    ! So far, whatever the species are, we set the "typical" lifetime to be two days
    !
    implicit none

    type(silja_interval) :: lifetime
    type(Tchem_rules_cbm4), intent(in) :: rules

    lifetime = one_day * 2.0

  end function fu_lifetime_cbm4

  logical function fu_if_specific_deposition_cbm4(rulescbm4) result(if_specific)
    implicit none
    type(Tchem_rules_cbm4), intent(in) :: rulescbm4
    if_specific = .false.     ! no specific deposition whatsoever
  end function fu_if_specific_deposition_cbm4

  !************************************************************************************
  
  subroutine check_low_mass_threshold_cb4(rulesCBM4, speciesTrn, nSpeciesTrn, ifForwardRun, &
                                        & low_cnc_threshold, low_mass_threshold)
    !
    ! This subroutine ensures that, in case of strong chemical connection, one substance can never
    ! be dropped to garbage while the other one still stays alive. In-essence, for cb4 we want 
    ! PAR to be always present if OLE are in
    !
    implicit none
    
    ! Imported parameters
    type(Tchem_rules_cbm4), intent(in) :: rulesCBM4
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
    if(fu_fails(indOLE /= int_missing,'No OLE in CB4 list','check_low_mass_threshold_cb4'))return
    if(fu_fails(indPAR /= int_missing,'No PAR in CB4 list','check_low_mass_threshold_cb4'))return
    !
    ! Check the relation
    !
    if(low_cnc_threshold(indPAR) > low_cnc_threshold(indOLE))then
      call msg('CB4: bad relation between OLE and PAR thresholds. Take OLE (the second):', &
                                            & low_cnc_threshold(indPAR), low_cnc_threshold(indOLE))
      low_cnc_threshold(indPAR) = low_cnc_threshold(indOLE)
      low_mass_threshold(indPAR) = low_mass_threshold(indOLE)
    endif
    
  end subroutine check_low_mass_threshold_cb4
  
  !************************************************************************************

  logical function fu_if_tla_required_cbm4(rules) result(required)
    ! Collect transformations' requests for tangent linearization. If tangent linear
    ! is not allowed, corresponding subroutine sets error. 
    implicit none
    type(Tchem_rules_cbm4), intent(in) :: rules
    
    ! call set_error('TLA not available for cbm4', 'fu_if_tla_required_cbm4')
     required = .true.
    
  end function fu_if_tla_required_cbm4


END MODULE chem_dep_cbm4

