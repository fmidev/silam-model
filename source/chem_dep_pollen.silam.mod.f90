MODULE chem_dep_pollen
  !
  ! This module contains a set of routines that describe the transformations
  ! of pollen and allergen. 
  ! Possible processes include pollen potency change in time,
  ! pollen break-up and allergen release due to humidity and pollution,
  ! allergen interactions with pollution and radiation, etc
  ! Detailed aerosol-dynamic interactions between the particles and pollen
  ! should be handled by aerosol dynamics module!
  !
  
  use cocktail_basic

  implicit none
  private
  

  ! public routines
  public inventoryPollen
  public registerSpeciesPollen
  public pollen_proc_input_needs
  public transform_pollen
  public set_chem_rules_pollen
  public fu_lifetime
  public fu_if_specific_deposition

  ! Private routines
  private fu_lifetime_pollen
  private fu_if_specific_dep_pollen


  interface fu_lifetime
    module procedure fu_lifetime_pollen
  end interface

  interface fu_if_specific_deposition
    module procedure fu_if_specific_dep_pollen
  end interface

  public fu_if_tla_required
  private fu_if_tla_required_pollen
  interface fu_if_tla_required
     module procedure fu_if_tla_required_pollen
  end interface


  !--------------------------------------------------------------------
  !
  !  Tchem_rules_pollen type definition
  !  Determines the rules of operating with pollen and allergen
  !
  type Tchem_rules_pollen
    private
    logical :: ifRHbreak, ifChemBreak, ifPotChange, ifBreak
    integer, dimension(:), pointer :: nReacts
    real, dimension(:), pointer :: fRatePotChange, fBreakRate, fThresRH, fBreakRateRH
    integer, dimension(:,:), pointer :: iPolBreaker
    real, dimension(:,:), pointer :: fPolBreakeRate
    type(silja_logical) :: defined
  end type Tchem_rules_pollen
  public Tchem_rules_pollen

  !-------------------------------------------------------------------------
  !
  ! Pointers to the meteorological fields, which can be requested
  ! by the corresponding transformation or deposition routines
  !
  real, dimension(:), private, pointer, save :: ptrRH
  type(field_4d_data_ptr), private, pointer, save :: fldPtrRH

  integer, private, pointer, save :: ind_RH  ! Used for addressing the meteo data

  !-------------------------------------------------------------------------
  !
  ! The main arrays pointing at the pollen and allergen species in the transport cocktail
  !
  integer, dimension(:), pointer, private, save :: indPoll_transp 
  integer, dimension(:), pointer, private, save :: indAlrgPoll_transp
  integer, dimension(:), pointer, private, save :: indAlrgFree_transp

  
  integer, private, save :: nPol = 0, nPolAlrg = 0, nFreeAlrg = 0

  ! The mark of this transformation type
  INTEGER, PARAMETER, public :: transformation_pollen = 5003


  CONTAINS


  !************************************************************************************

  subroutine inventoryPollen(rules, &
                           & speciesEmis, speciesTransp, speciesShortlived, speciesAerosol,&
                           & nSpeciesEmis, nspeciesTransp, nspeciesShortlived, nspeciesAerosol, &
                           & iClaimedSpecies)
    implicit none
    ! 
    ! Claim all pollen and allergen species. 
    ! It is assumed that all necessary species are emitted and no new ones are created here.
    ! It might be an option to request the free allergen field here, but then this place would
    ! have to know the particle size which might be species dependent and generally messy, 
    ! so it's easier to let the source worry about it and here just deal with the conciquences.
    !
    type(Tchem_rules_pollen), intent(in) :: rules
    type(silam_species), dimension(:), pointer :: speciesEmis, speciesTransp, speciesShortlived, &
                                                & speciesAerosol
    integer, intent(in) :: nspeciesEmis
    integer, intent(inout) :: nspeciesTransp, nspeciesShortlived, nspeciesAerosol
    integer, dimension(:), intent(inout) :: iClaimedSpecies

    ! Local variables
    integer :: i
    !
    ! Scan the emission species, find pollen material and put it to transport mass map, also 
    ! claiming ownership. 
    !
    do i = 1, nspeciesEmis
      if(fu_pollen_material_type(speciesEmis(i)%material) == pollenGrains .or. &
      &  fu_pollen_material_type(speciesEmis(i)%material) == pollenAllergen .or. &
      &  fu_pollen_material_type(speciesEmis(i)%material) == freeAllergen )then
        call addSpecies(speciesTransp, nSpeciesTransp, (/speciesEmis(i)/), 1, .true.)
        if(iClaimedSpecies(i) < 0)then
          call msg('Pollen owns:' + fu_substance_name(speciesEmis(i)))
          iClaimedSpecies(i) = transformation_pollen
        else
          call msg('Cannot claim ownership because he claimed it already:',iClaimedSpecies(i))
          call set_error('Cannot claim ownership because someone claimed it already','inventoryPollen')
          return
        endif
        if(fu_pollen_material_type(speciesEmis(i)%material) == pollenGrains)then
          nPol = nPol + 1      
        elseif(fu_pollen_material_type(speciesEmis(i)%material) == pollenAllergen)then
          nPolAlrg = nPolAlrg + 1 
        elseif(fu_pollen_material_type(speciesEmis(i)%material) == freeAllergen)then
          nFreeAlrg = nFreeAlrg + 1 
        endif
      endif 
    end do  ! emission species
   
    ! There can't be more allergens than there is pollens
    if(nPolAlrg>nPol .or. nFreeAlrg>nPol)then
      call set_error('More allergen species found than pollen','inventoryPollen')
    endif

  end subroutine inventoryPollen


  !***********************************************************************

  subroutine registerSpeciesPollen(rules, &
                                 & speciesTransp, speciesShortlived, speciesAerosol,&
                                 & nspeciesTransp, nspeciesShortlived, nspeciesAerosol)
    ! 
    ! Find the indices of pollen and allergen species and store them to the module
    ! variable, colocating the ones for same species and checking for consistency
    ! for existing species and requested transformations. 
    ! Also have to catch the chemical species requested for the transformations and
    ! remember them in iPolBreaker and fPolBreakeRate matrices 
    !
    implicit none
    type(Tchem_rules_pollen), intent(inOut) :: rules
    type(silam_species), dimension(:), pointer :: speciesTransp, speciesShortlived, speciesAerosol
    integer, intent(in) :: nspeciesTransp, nspeciesShortlived, nspeciesAerosol
    
    integer :: iSpTrans, iSpPol, iSpAlrg, nReacts, iTmp
    character(len=substNmLen), dimension(:), pointer :: agents
    real, dimension(:), pointer :: rates
    type(silam_material), pointer :: materialPtr

    if(nPol > 0)then
      allocate(indPoll_transp(nPol))
    else
      nullify(indPoll_transp)
      return
    endif
    if(nPolAlrg > 0)then
      allocate(indAlrgPoll_transp(nPol))
      indAlrgPoll_transp(:) = int_missing
    else
      nullify(indAlrgPoll_transp)
    endif
    if(nFreeAlrg > 0)then
      allocate(indAlrgFree_transp(nPol))
      indAlrgFree_transp(:) = int_missing
    else
      nullify(indAlrgFree_transp)
    endif
    
    iSpPol = 0
    nReacts = 0
    do iSpTrans = 1, nSpeciesTransp
      if(fu_pollen_material_type(speciesTransp(iSpTrans)%material) == pollenGrains)then
        iSpPol = iSpPol + 1
        indPoll_transp(iSpPol) = iSpTrans
        ! Check for reactions
        if(fu_pollen_break_rate(speciesTransp(iSpTrans)%material) > 0.0 .or. &
         & fu_pollen_RH_break_rate(speciesTransp(iSpTrans)%material) > 0.0 .or. &
         & fu_n_pollen_break_chem_reacts(speciesTransp(iSpTrans)%material) > 0)then
        ! both free and pollen allergens must exist
          if(.not. (associated(indAlrgPoll_transp) .and. associated(indAlrgFree_transp)))then
            call set_error('pollen reactions requested but no allergen fields given','registerSpeciesPollen')
            return
          endif
          do iSpAlrg = 1, nSpeciesTransp
            if(fu_pollen_material_type(speciesTransp(iSpAlrg)%material) == pollenAllergen .and. &
               & fu_pollen_taxon_name(speciesTransp(iSpTrans)%material) == &
               & fu_pollen_taxon_name(speciesTransp(iSpAlrg)%material) .and. &
               & fu_mode(speciesTransp(iSpTrans)) == fu_mode(speciesTransp(iSpAlrg)))then
              indAlrgPoll_transp(iSpPol) = iSpAlrg
            elseif(fu_pollen_material_type(speciesTransp(iSpAlrg)%material) == freeAllergen .and. &
               & fu_pollen_taxon_name(speciesTransp(iSpTrans)%material) == &
               & fu_pollen_taxon_name(speciesTransp(iSpAlrg)%material))then
              indAlrgFree_transp(iSpPol) = iSpAlrg
            endif
          enddo
          if(indAlrgPoll_transp(iSpPol) == int_missing .or. indAlrgFree_transp(iSpPol) == int_missing)then
            call set_error('pollen reactions requested but no allergen fields given','registerSpeciesPollen')
            return
          endif
        elseif(fu_pollen_potency_change_rate(speciesTransp(iSpTrans)%material) > 0.0)then
        ! at least pollen allergen must exist
          if(.not. associated(indAlrgPoll_transp))then
            call set_error('pollen reactions requested but no allergen field given','registerSpeciesPollen')
            return
          endif
          do iSpAlrg = 1, nSpeciesTransp
            if(fu_pollen_material_type(speciesTransp(iSpAlrg)%material) == pollenAllergen .and. &
               & fu_pollen_taxon_name(speciesTransp(iSpTrans)%material) == &
               & fu_pollen_taxon_name(speciesTransp(iSpAlrg)%material) .and. &
               & fu_mode(speciesTransp(iSpTrans)) == fu_mode(speciesTransp(iSpAlrg)))then
              indAlrgPoll_transp(iSpPol) = iSpAlrg
            endif
          enddo
          if(indAlrgPoll_transp(iSpPol) == int_missing)then
            call set_error('pollen reactions requested but no allergen field given','registerSpeciesPollen')
            return
          endif
        else
          call msg('No pollen transformations requested')
        endif
        ! count some processes ..
        if(fu_pollen_RH_break_rate(speciesTransp(iSpTrans)%material) > 0.0) rules%ifRHbreak = .true.
        if(fu_pollen_break_rate(speciesTransp(iSpTrans)%material) > 0.0)rules%ifBreak = .true.
        if(.not. (fu_pollen_potency_change_rate(speciesTransp(iSpTrans)%material) .eps. 0.0))rules%ifPotChange = .true.
        if(fu_n_pollen_break_chem_reacts(speciesTransp(iSpTrans)%material) > nReacts) &
            & nReacts = fu_n_pollen_break_chem_reacts(speciesTransp(iSpTrans)%material)
        
      endif
    end do
    
    ! Allocate and fill the arrays for requested transformations
    if(rules%ifRHbreak)then 
      allocate(rules%fBreakRate(nPol),rules%fThresRH(nPol))
      do iSpPol = 1, nPol
        rules%fBreakRate(iSpPol) = fu_pollen_RH_break_rate(speciesTransp(indPoll_transp(iSpPol))%material)
        rules%fThresRH(iSpPol) = fu_pollen_RH_break_thres(speciesTransp(indPoll_transp(iSpPol))%material)
      enddo
    endif
    
    if(rules%ifPotChange)then 
      allocate(rules%fRatePotChange(nPol))
      do iSpPol = 1, nPol
        rules%fRatePotChange(iSpPol) = fu_pollen_potency_change_rate(speciesTransp(indPoll_transp(iSpPol))%material)
      enddo
    endif
    
    if(rules%ifBreak)then
      allocate(rules%fBreakRate(nPol))
      do iSpPol = 1, nPol
        rules%fBreakRate(iSpPol) = fu_pollen_break_rate(speciesTransp(indPoll_transp(iSpPol))%material)
      enddo
    endif
    
    ! Reactions with other chemicals stored in iPolBreaker and fPolBreakeRate matrices
    if(nReacts > 0)then
      rules%ifChemBreak = .true.
      allocate(rules%iPolBreaker(nPol,nReacts), rules%fPolBreakeRate(nPol,nReacts), rules%nReacts(nPol))
      rules%iPolBreaker(:,:) = int_missing
      rules%fPolBreakeRate(:,:) = 0.0
      do iSpPol = 1, nPol
        call pollen_break_chem_reacts_ptr(speciesTransp(indPoll_transp(iSpPol))%material, agents, rates, nReacts)
        rules%nReacts(iSpPol) = nReacts
        do iTmp = 1, nReacts
          materialPtr => fu_get_material_ptr(agents(iTmp))
          do iSpTrans = 1, nSpeciesTransp
            if(associated(materialPtr, speciesTransp(iSpTrans)%material))then
              rules%iPolBreaker(iSpPol,iTmp) = iSpTrans
              exit
            endif
          enddo
          if(rules%iPolBreaker(iSpPol,iTmp) == int_missing)then
            call set_error('Species requested for pollen transformation not found in the run: ' + &
                                                              & agents(iTmp),'registerSpeciesPollen')
            return
          endif
          rules%fPolBreakeRate(iSpPol,iTmp) = rates(iTmp)
        enddo
      enddo
    endif 
    
    ! Cycle the transport species once again to catch the allergen species that were not requested
    ! for pollen transformations. No idea why ..
    if(nPolAlrg > 0)then
      do iSpAlrg = 1, nSpeciesTransp
        if(fu_pollen_material_type(speciesTransp(iSpAlrg)%material) == pollenAllergen)then
          if(any(indAlrgPoll_transp == iSpAlrg))cycle ! already there
          do iSpPol = 1, nPol
            if(fu_pollen_taxon_name(speciesTransp(indPoll_transp(iSpPol))%material) == &
               & fu_pollen_taxon_name(speciesTransp(iSpAlrg)%material) .and. &
               & fu_mode(speciesTransp(indPoll_transp(iSpPol))) == fu_mode(speciesTransp(iSpAlrg)))then
              indAlrgPoll_transp(iSpPol) = iSpAlrg
              exit
            endif
          enddo
        endif
      enddo
    endif    
    if(nFreeAlrg > 0)then
      do iSpAlrg = 1, nSpeciesTransp        
        if(fu_pollen_material_type(speciesTransp(iSpAlrg)%material) == freeAllergen)then
          if(any(indAlrgFree_transp == iSpAlrg))cycle ! already there
          do iSpPol = 1, nPol
            if(fu_pollen_taxon_name(speciesTransp(indPoll_transp(iSpPol))%material) == &
             & fu_pollen_taxon_name(speciesTransp(iSpAlrg)%material))then
              indAlrgFree_transp(iSpPol) = iSpAlrg
              exit
            endif
          enddo
        endif
      enddo
    endif 

  end subroutine registerSpeciesPollen

  !*****************************************************************

  subroutine pollen_proc_input_needs(rulesPollen, meteo_input_local)
    !
    ! Returns input needs for the transformation routines. 
    !
    implicit none

    ! Imported parameters
    type(Tchem_rules_pollen), intent(in) :: rulesPollen
    type(Tmeteo_input), intent(out), target :: meteo_input_local

    ! Local variables
    integer :: iTmp

    meteo_input_local = meteo_input_empty 

    if(.not.(rulesPollen%defined == silja_true))return

    ! Prepare the local mapping for the meteo input
    !
    if(rulesPollen%ifRHbreak)then
      meteo_input_local%nQuantities = 1

      ! 
      ! .. and store it into the local mapping. Note that q_static is not needed for pollen
      !
      meteo_input_local%quantity(1) = relative_humidity_flag    ! For self-degradation
      meteo_input_local%q_type(1) = meteo_dynamic_flag
      !
      ! Finally, prepare the named indices, which will be used in the actual calculations
      !
      ind_RH => meteo_input_local%idx(1)
    endif
  end subroutine pollen_proc_input_needs


  !*************************************************************************

  subroutine get_Rs_pollen_crd(ptrRs, ifRs_meteo_depend)
    !
    ! Returns the surface resistance for the pollen substance
    !
    implicit none

    ! Imported parameters
    real, dimension(:), intent(out) :: ptrRs
    logical, intent(out) :: ifRs_meteo_depend

    ! Local variables
    integer :: iSpecies

    do iSpecies = 1, nPol
      ptrRs(indPoll_transp(iSpecies)) = 0.
    enddo
    ifRs_meteo_depend = .false.

  end subroutine get_Rs_pollen_crd


  !*******************************************************************

  subroutine transform_pollen(vMass, rulesPollen, metdat, timestep_sec, print_it) 
    !
    ! Computes the transformation for the pollen species in the given array of masses 
    ! using the provided meteodata
    !
    implicit none

    ! Imported parameters
    real, dimension(:), intent(inout) :: vMass
    real, dimension(:), intent(in) :: metdat
    type(Tchem_rules_pollen), intent(in) :: rulesPollen
    real, intent(in) :: timestep_sec
    logical, intent(out) :: print_it

    ! Local parameters
    integer :: iSpecies, iReact
    real :: fBreak

    print_it = .false.  ! set to true and chemistry manager will give complete dump for this cell

    if(.not. rulesPollen%defined == silja_true)return

    ! cycle over pollen species
    do iSpecies = 1, nPol
      ! linear growth of allergen in pollen grain
      if(rulesPollen%ifPotChange)then
        vMass(indAlrgPoll_transp(iSpecies)) = vMass(indAlrgPoll_transp(iSpecies)) + &
            & vMass(indPoll_transp(iSpecies)) * rulesPollen%fRatePotChange(iSpecies) * timestep_sec
      endif
      ! Pollen break processes (simple, RH and chemistry dependent)
      ! As all of them result in loss of pollen and allergen mofing from pollen allergen to 
      ! free allergen, the rates are added up 
      
      fBreak = 0.0
      if(rulesPollen%ifBreak)then
        fBreak = fBreak + rulesPollen%fBreakRate(iSpecies)
      endif
      if(rulesPollen%ifRHbreak)then
        if(metdat(ind_RH) > rulesPollen%fThresRH(iSpecies))then
          fBreak = fBreak + (metdat(ind_RH) - rulesPollen%fThresRH(iSpecies)) &
                                                          & * rulesPollen%fBreakRateRH(iSpecies)
        endif
      endif
      if(rulesPollen%ifChemBreak)then
        do iReact = 1, rulesPollen%nReacts(iSpecies)
          fBreak = fBreak + vMass(rulesPollen%iPolBreaker(iSpecies,iReact)) * &
                                                   & rulesPollen%fPolBreakeRate(iSpecies,iReact)
        enddo
      endif
      if(fBreak > 0.)then
        fBreak = 1.0 - exp(- fBreak * timestep_sec)
        vMass(indPoll_transp(iSpecies)) = vMass(indPoll_transp(iSpecies)) *  (1.0 - fBreak)
        vMass(indAlrgFree_transp(iSpecies)) = vMass(indAlrgFree_transp(iSpecies)) + &
                                                   & vMass(indAlrgPoll_transp(iSpecies)) * fBreak
        vMass(indAlrgPoll_transp(iSpecies)) = vMass(indAlrgPoll_transp(iSpecies)) * (1.0 - fBreak)
      endif
    end do

  end subroutine transform_pollen



  !*******************************************************************************
  !*******************************************************************************
  !
  !  Passive rules
  !
  !*******************************************************************************
  !*******************************************************************************

  !***********************************************************************

  subroutine set_chem_rules_pollen(nlSetup, rulesPollen)
    !
    ! Creates the set of rules, in accordance with which the computations are to
    ! be made. Input - SILAM namelist
    !
    implicit none

    ! Imported parameter
    type(Tsilam_namelist), pointer :: nlSetup 
    type(Tchem_rules_pollen), intent(out) :: rulesPollen

    rulesPollen%defined = silja_false
    
    rulesPollen%ifRHbreak = .false.
    rulesPollen%ifChemBreak = .false.
    rulesPollen%ifPotChange = .false.
    rulesPollen%ifBreak = .false.
    nullify(rulesPollen%fRatePotChange, rulesPollen%fBreakRate, rulesPollen%fThresRH, rulesPollen%fBreakRateRH)
    nullify(rulesPollen%iPolBreaker, rulesPollen%fPolBreakeRate, rulesPollen%nReacts)

    rulesPollen%defined = silja_true

  end subroutine set_chem_rules_pollen


  !********************************************************************************

  function fu_lifetime_pollen(rules)result(lifetime)
    !
    ! A typical life time due to degradation and deposition for aerosol
    ! So far, whatever the species are, we set the "typical" lifetime to be two days
    !
    implicit none

    type(silja_interval) :: lifetime
    type(Tchem_rules_pollen), intent(in) :: rules

    lifetime = one_day * 2.0

  end function fu_lifetime_pollen

  !**********************************************************************************

  logical function fu_if_specific_dep_pollen(rulesPassive)
    implicit none
    type(Tchem_rules_pollen), intent(in) :: rulesPassive
    fu_if_specific_dep_pollen = .false.     ! no specific deposition whatsoever
  end function fu_if_specific_dep_pollen

  !************************************************************************************

  logical function fu_if_tla_required_pollen(rules) result(required)
    ! Collect transformations' requests for tangent linearization. If tangent linear
    ! is not allowed, corresponding subroutine sets error. 
    implicit none
    type(Tchem_rules_pollen), intent(in) :: rules
    
    call set_error('TLA not available for pollen transormation', 'fu_if_tla_required_pollen')
    required = .true.

  end function fu_if_tla_required_pollen


END MODULE chem_dep_pollen

