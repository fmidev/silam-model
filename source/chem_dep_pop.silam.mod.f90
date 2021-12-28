MODULE cocktail_persistent_organics
  !
  ! This module contains description of the Persistent Organic Pollutants (POP) chemistry and deposition.
  ! To be exact, POP can be anything, which has gas-particle partitioning and degradation in the 
  ! atmosphere. Later, re-emission from soil/vegetation and water will be added, with appropriate
  ! treatment of degradation in these media.
  !
  ! All units: SI, unless otherwise stated
  !
  ! Code owner Mikhail Sofiev mikhail.sofiev@fmi.fi
  !
  ! Language: ANSI FORTRAN-90 (or close to it)
  !
  use cocktail_basic

  implicit none

  !
  ! PUBLIC routines of POP cocktail 
  !
  public init_POP_emis_chemicals
  public full_subst_lst_4ckt_persis_org
  public POP_cocktail_input_needs
  public POP_EmisFlds ! For emission computation
  public POP_ChemistryFlds ! For transformation computation
  public InitPOP_EmissionFields
  public InitPOP_ChemistryFields
  public Init_POP_CocktailMap
  public FinalisePOP_EmissionFields
  public FinalisePOP_ChemistryFields
  public FinalisePOP_CocktailMap

  public fu_NModes_POP
  public fu_gas_phase_idx_persist_org
  public set_ModeMap_POP

  public prepare_POP_transdep
  public make_dry_deposition_POP
  public get_Rb_persistent_organics
  public get_Rs_persistent_organics_crd
  public refine_Rs_persistent_organics
  public make_scavenging_POP
  public make_transformation_POP
  public link_ems_trn_subst_pop

  ! PUBLIC routines of Tchem_rules_POP
  public set_chem_rules_POP  ! returns the pointer to the rules
  public set_missing
  public set_low_mass_threshold
  public fu_low_mass_threshold
  public fu_zero_mass_threshold
  public fu_lifetime
  public fu_if_tla_if_required

  !
  ! Private routines of the POP cocktail and chemical rules
  !
  private transform_POP_single_val
  private transform_POP_map
  private scavenge_POP_single_val
  private scavenge_POP_map

  private set_missing_chem_rules_POP
  private fu_low_mass_threshold_POP
  private fu_zero_mass_threshold_POP
  private fu_lifetime_POP
  private set_low_mass_thrsh_pop
  private fu_tla_if_required_pop

  interface make_transformation_POP
    module procedure transform_POP_single_val
    module procedure transform_POP_map
  end interface

  interface make_scavenging_POP
    module procedure scavenge_POP_single_val
    module procedure scavenge_POP_map
  end interface

  interface set_missing
    module procedure set_missing_chem_rules_POP
  end interface

  interface set_low_mass_threshold
    module procedure set_low_mass_thrsh_POP
  end interface

  interface fu_low_mass_threshold
    module procedure fu_low_mass_threshold_POP
  end interface

  interface fu_zero_mass_threshold
    module procedure fu_zero_mass_threshold_POP
  end interface

  interface fu_lifetime
    module procedure fu_lifetime_POP
  end interface

  interface init_tangent_linear
     module procedure init_tlin_pop
  end interface

  interface fu_if_tla_required
     module procedure fu_if_tla_required_pop
  end interface


  !--------------------------------------------------------------------
  !
  !  Tchem_rules_POP type definition
  !  Determines the rules of operating with the POPs
  !
  type Tchem_rules_POP
    private
    integer :: nPrecFlds
!    real, dimension(:), pointer :: degradation_rate
    real :: massLowThreshold
    type(silja_logical) :: defined
  end type Tchem_rules_POP

  !-------------------------------------------------------------------------
  !
  ! There can many POPs in one run, so the index is an array
  ! The main array pointing at the passive species in teh transport cocktail
  !
  integer, dimension(:), pointer, private, save :: indPOP_1d, indPOPTranspSubst
  integer, private, save :: nPOPTranspSubst, nPOPTranspSpecies

  !
  ! Stuff needed for computations of transformation and deposition.
  ! Pointers are set during prepare_POP_transdep
  !
  integer, private, save :: indTempr2m
  REAL, DIMENSION(:), private, save, POINTER :: ptrTempr2m, ptrCloudCover, ptrFricVel
  type(field_4d_data_ptr), private, save, pointer :: ptrTempr, ptrScavStd, ptrRelHumidity

  real, dimension(:,:), private, save, pointer :: ptrDegradation_rate


CONTAINS


  !*********************************************************************

  subroutine init_POP_emis_chemicals(nlStdSetup, chSubstNames, nSubst, SubstanceList, fDefaultDensity)
    !
    ! Initializes the acid basic chemicals. 
    ! Actually, we just request for the global chemical database to be read.
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist), pointer :: nlStdSetup
    character(len=*), dimension(:), intent(in) :: chSubstNames
    integer, intent(in) :: nSubst
    type(TSubstanceList), intent(inout) :: SubstanceList
    real, intent(in) :: fDefaultDensity

    if(.not. fu_if_chemical_data_loaded()) call init_chemical_materials(nlStdSetup, SubstanceList)
    if(error)return

    !
    ! Now, mark the given species as used in teh emission cocktail. Add if something is missing
    !
    call include_subst_into_run(chSubstNames, nSubst, species_persistent_organics, &
                              & emission_substance_set, SubstanceList, fDefaultDensity )
    if(error)then
      call set_error('Failed to include the pollen into run','init_POP_emis_chemicals')
    endif

  end subroutine init_POP_emis_chemicals


  !*******************************************************************

  subroutine full_subst_lst_4ckt_persis_org(lstSubstNm, nSubst, SubstanceList, fDefaultDensity)
    !
    ! From a given chemical list of emission, explores the full set of species for transport
    ! For POPs cocktail it is the same - there is no production of new species
    !
    implicit none

    ! Imported parameters
    character(len=*), dimension(:), pointer :: lstSubstNm
    integer, intent(in) :: nSubst
    type(TSubstanceList), intent(inout) :: SubstanceList
    real, intent(in) :: fDefaultDensity

    ! Simply send the list of materials to include them into transport cocktail
    !
    call include_subst_into_run(lstSubstNm, nSubst, species_persistent_organics, &
                              & transport_substance_set, SubstanceList, fDefaultDensity)

  end subroutine full_subst_lst_4ckt_persis_org


  !**********************************************************************

  subroutine link_ems_trn_subst_POP(cocktEmis, cocktTransp, SubstLst, ChemRunSetup)
    !
    ! Links emission and transport cocktails and sets the local indices for all pollen
    ! species (NOT substances!!) in the transport cocktail
    !
    implicit none

    ! Imported parameters
    type(silam_cocktail), intent(in) :: cocktEmis, cocktTransp
    type(TsubstanceList), intent(in) :: SubstLst
    type(TChemicalRunSetup), intent(inout) :: ChemRunSetup

    ! Local variables
    integer :: iS, iMode, iTmp, iTmp_1d, iSubstTr

    !
    ! Link the transport cocktail and the local transformation routines
    ! Fill-in the local arrays of indices for sea salt species and local counter nSeaSaltTranspSubst
    ! Note that for sea salt substances, one substance may mean several species due to size modes
    !
    nPOPTranspSubst = 0
    nPOPTranspSpecies = 0
    do iTmp = 1, SubstLst%nSubstSet
      if(SubstLst%iSubstTypes(iTmp) == species_persistent_organics)then
        nPOPTranspSubst = nPOPTranspSubst + 1
        nPOPTranspSpecies = nPOPTranspSpecies + cocktTransp%nModesTotSubst(iTmp)
      endif
    end do
    allocate(indPOP_1d(nPOPTranspSpecies), indPOPTranspSubst(nPOPTranspSubst), stat=iTmp)
    if(iTmp /= 0)then
      call set_error('Failed to allocate local indices','link_ems_trn_subst_POP')
      call report(cocktTransp)
      return
    endif
    !
    ! Now, fill-in the local array of pointers. Note the difference between the 1d coding
    ! array iSpecies and just substance number  in the cocktail
    !
    iTmp = 0
    iTmp_1d = 0
    do iS = 1, SubstLst%nSubstSet
      if(SubstLst%iSubstTypes(iS) == species_persistent_organics)then
        iTmp = iTmp + 1
        indPOPTranspSubst(iTmp) = iS
        do iMode = 1, cocktTransp%nModesTotSubst(iS)
          iTmp_1d = iTmp_1d + 1
          indPOP_1d(iTmp_1d) = cocktTransp%iSpecies(iS,iMode,1)
        end do
      endif
    end do
    !
    ! Go through emission cocktail step-by-step. Whenever find POP species, 
    ! find out the counterpart in transport cocktail and set the link
    !
    do iS = 1, ChemRunSetup%lstEmis%nSpeciesEms
      if(ChemRunSetup%lstEmis%iSpeciesTypes(iS) == species_persistent_organics)then
        allocate(ChemRunSetup%lstEmis%refsToTranspSpecies(iS)%ptr, stat=iTmp)
        if(iTmp /= 0)then
          call msg('Failed ref ptr for emission species:',iS)
          call set_error('Failed to allocate reference ptr','link_ems_trn_subst_POP')
          call report(cocktEmis)
          return
        endif
        allocate(ChemRunSetup%lstEmis%refsToTranspSpecies(iS)%ptr%fract(1), &
               & ChemRunSetup%lstEmis%refsToTranspSpecies(iS)%ptr%indSpeciesTo(1), stat=iTmp)
        if(iTmp /= 0)then
          call msg('Failed ref coef for emission species:',iS)
          call set_error('Failed to allocate reference coefficients','link_ems_trn_subst_POP')
          call report(cocktEmis)
          return
        endif
        ChemRunSetup%lstEmis%refsToTranspSpecies(iS)%ptr%nRefSpecies = 1  ! Link 1:1
        ChemRunSetup%lstEmis%refsToTranspSpecies(iS)%ptr%fract(1) = 1.
        !
        ! Sea salt species simply goes from emission to transport - get the counterpart!
        !
        iSubstTr = 1
        do 
          if(cocktTransp%iSrcInventory(iSubstTr) == cocktEmis%iSubst(iS)) exit
          iSubstTr = iSubstTr +1
          if(iSubstTr > cocktTransp%nSubstTot)then
            call msg('No transport substance for emission species:',iS)
            call set_error('No transport substance for emission','link_ems_trn_subst_POP')
            return
          endif
        end do
        ChemRunSetup%lstEmis%refsToTranspSpecies(iS)%ptr%indSpeciesTo(1) = &
                       & cocktTransp%iSpecies(iSubstTr, cocktEmis%iMode(iS), cocktEmis%iWave(iS))
      endif  ! pollen_species
    end do  ! emision species

  end subroutine link_ems_trn_subst_POP


  !*****************************************************************

  subroutine POP_cocktail_input_needs(rulesPOP, rulesAerosol, q_dynamic, q_st)
    !
    ! Returns input needs for the transformation routines. 
    !
    implicit none

    ! Imported parameters
    type(Tchem_rules_POP), intent(in) :: rulesPOP
    type(Taerosol_rules),  intent(in) :: rulesAerosol

    integer, dimension(:), intent(out) :: q_dynamic, q_st

    ! Local variables
    integer :: iDyn, iST

    q_dynamic = int_missing
    q_st = int_missing

    !
    ! Each POP can be in two phases, so we have to take this into account
    !
    call aerosol_input_needs(rulesAerosol, q_dynamic, q_st)

    !
    ! Find free space in both dynamic and static arrays
    !
    iDyn = 1
    do while (iDyn < size(q_dynamic))
      if(.not. fu_known_quantity(q_dynamic(iDyn))) exit
      iDyn = iDyn + 1
    end do
    if(iDyn > size(q_dynamic)-5)then
      call set_error('Too small dynamic array','POP_cocktail_input_needs')
      return
    endif

    iST = 1
    do while (iST < size(q_st))
      if(.not. fu_known_quantity(q_st(iST))) exit
      iST = iST + 1
    end do
    if(iST > size(q_st)-1)then
      call set_error('Too small static array','POP_cocktail_input_needs')
      return
    endif
    !
    ! Add needed dynamic quantities
    !
    if(.not. any(q_dynamic(1:iDyn-1) == temperature_flag))then
      q_dynamic(iDyn) = temperature_flag   ! For convertion
      iDyn = iDyn + 1
    endif

    if(.not. any(q_dynamic(1:iDyn-1) == temperature_2m_flag))then
      q_dynamic(iDyn) = temperature_2m_flag   ! for dry deposition and conversion
      iDyn = iDyn + 1
    endif

    if(.not. any(q_dynamic(1:iDyn-1) == total_cloud_cover_flag))then
      q_dynamic(iDyn) = total_cloud_cover_flag   ! for OH concentration
      iDyn = iDyn + 1
    endif


    if(.not. any(q_dynamic(1:iDyn-1) == friction_velocity_flag))then
      q_dynamic(iDyn) = friction_velocity_flag
      iDyn = iDyn + 1
    endif

    if(rulesAerosol%ifHumidityDependent)then
      q_dynamic(iDyn) = relative_humidity_flag
      iDyn = iDyn + 1
    else
      nullify(ptrRelHumidity)
    endif

    !
    ! Single-time quantities
    !
    if(.not. any(q_st(1:iST-1) == fraction_of_land_flag))then
      q_st(iST) = fraction_of_land_flag
      iST = iST + 1
    endif

    if(.not. any(q_st(1:iST-1) == large_scale_rain_int_flag))then
      q_st(iST) = large_scale_rain_int_flag   ! for wet deposition 
      iST = iST + 1
    endif

    if(.not. any(q_st(1:iST-1) == convective_rain_int_flag) .and. &
     & rulesPOP%nPrecFlds == 2)then
      q_st(iST) = convective_rain_int_flag   ! for wet deposition 
      iST = iST + 1
    endif
   
    if(.not. any(q_st(1:iST-1) == scavenging_coefficient_flag))then
      q_st(iST) = scavenging_coefficient_flag
      iST = iST + 1
    endif

  end subroutine POP_cocktail_input_needs

  
  !**********************************************************************

  subroutine POP_EmisFlds(Rules, q_field, q_cocktail_map, iNbrOfFields, iNbrOfCocktailMaps)
    !
    ! Returns the number of fields needed for POP emission computations
    !
    implicit none

    ! Imported parameter
    type(Tchem_rules_POP), intent(in) :: Rules
    integer, dimension(:), pointer :: q_field, q_cocktail_map
    integer, intent(inout) :: iNbrOfFields, iNbrOfCocktailMaps

  end subroutine POP_EmisFlds


  !**********************************************************************

  subroutine POP_ChemistryFlds(Rules, arQuantities, iNbrOfFlds)
    !
    ! Returns the number of fields needed for POP transformation and
    ! deposition computations. 
    !
    implicit none

    ! Imported parameter
    type(Tchem_rules_POP), intent(in) :: Rules
    integer, dimension(:), pointer :: arQuantities
    integer, intent(out) :: iNbrOfFlds

    iNbrOfFlds = 0
    arQuantities(1) = int_missing

  end subroutine POP_ChemistryFlds


  !**********************************************************************

  subroutine InitPOP_EmissionFields(Rules, wdr)
    !
    ! Actually creates the fields serving the emission computations
    ! and drops them to the internal dispersion buffer
    !
    implicit none

    ! Imported parameter
    type(Tchem_rules_POP), intent(in) :: Rules
    type(silja_wdr), intent(in) :: wdr

    ! No specific fields are needed

  end subroutine InitPOP_EmissionFields


  !**********************************************************************

  subroutine InitPOP_ChemistryFields(Rules, wdr, cocktail_template)
    !
    ! Actually creates the fields serving the chemistry and deposition
    ! computations and drops them to the internal dispersion buffer, if needed.
    !
    ! For POPs, here we create and precompute the degradation coefficient??????????
    !
    implicit none

    ! Imported parameter
    type(Tchem_rules_POP), intent(in) :: Rules
    type(silja_wdr), intent(in) :: wdr
    TYPE(silam_cocktail), intent(in) :: cocktail_template

    ! Local variables
    integer :: iTmp, iSubst, iTempr
    real :: Ea_NPT,alpha,beta

    allocate(ptrDegradation_rate(fu_nbr_of_subst(cocktail_template),201),stat=iTmp)
    if(iTMp /= 0)then
      call set_error('Failed to allocate space for degradation coefficient','InitPOP_ChemistryFields')
      return
    endif
    !
    ! Degradation rate description valid for 200-400 K
    ! T.H.Lay, J.W.Bozzelli, J.H.Seinfield (1996), J.Phys.Chem., 100, 6543-6554
    !
    do iTmp = 1, nPOPTranspSubst
      iSubst = indPOPTranspSubst(iTmp)
      !
      ! Get reaction parameters for the specific substance
      !
      call get_reaction_parameters(cocktail_template%materials(iSubst)%ptr,'OH','',Ea_NPT,alpha,beta)
      if(error)return
      !
      ! Array index runs from 1 to 201, which corresponds to 200 to 400 K, so the index is T-199
      !
      do iTempr = 1, 201
        ptrDegradation_rate(iSubst,iTempr) = beta * (199.+iTempr)**alpha * &
                                         & exp(-Ea_NPT/(gas_constant_dryair * (199.+iTempr)))
      end do   ! temperatures
    end do  ! species


  end subroutine InitPOP_ChemistryFields


  !***************************************************************************
  
  subroutine Init_POP_CocktailMap(Rules, ptrCocktailMap)
    !
    ! Should POP operations request some internal field,
    ! here it must store it for local pointer. This function is called once.
    ! Map pointers never change, so they are never restored
    !
    implicit none

    ! Imported parameters
    type(Tchem_rules_POP), intent(in) :: rules
    type(Tmass_map), pointer :: ptrCocktailMap

  end subroutine Init_POP_CocktailMap


  !**************************************************************************

  subroutine FinalisePOP_EmissionFields(Rules, disp_buf)
    !
    ! Finalises chemistry-related internal POP fields.
    ! Can store them for the next run or alike. Called once at
    ! the very end of the model run
    !
    implicit none

    ! Imported parameters
    type(Tchem_rules_POP), intent(in) :: Rules
    type(Tfield_buffer), pointer :: disp_buf

  end subroutine FinalisePOP_EmissionFields


  !***************************************************************************

  subroutine FinalisePOP_ChemistryFields(Rules, disp_buf, wdr)
    !
    ! Finalises chemistry-related internal POP fields.
    ! Can stoe them for the next run or alike. Called once at
    ! the very end of the model run
    !
    implicit none

    ! Imported parameters
    type(Tchem_rules_POP), intent(in) :: Rules
    type(Tfield_buffer), pointer :: disp_buf
    type(silja_wdr), intent(in) :: wdr

    ! There are no fields and nothing to finalise

  end subroutine FinalisePOP_ChemistryFields


  !**************************************************************************

  subroutine FinalisePOP_CocktailMap(Rules, ptrCocktailMap)
    !
    ! Finalises chemistry-related internal POP fields.
    ! Can stoe them for the next run or alike. Called once at
    ! the very end of the model run
    !
    implicit none

    ! Imported parameters
    type(Tchem_rules_POP), intent(in) :: Rules
    type(Tmass_map), pointer :: ptrCocktailMap

    ! There are no fields and nothing to finalise

  end subroutine FinalisePOP_CocktailMap



  !***********************************************************************

  integer function fu_NModes_POP(chName, aerosol)
    !
    ! Determines the number of modes for the given substance, which is always 2:
    ! all POPs can fluctuate between the fine-mode particles and gas phase
    !
    implicit none

    ! Imported parameters
    character(len=*), intent(in) :: chName
    type(Taerosol), intent(in) :: aerosol

    fu_NModes_POP = 2

  end function fu_NModes_POP


  !***********************************************************************

  integer function fu_gas_phase_idx_persist_org(aerosol)
    !
    ! Determines the number of modes for the given passive substance, which is always 1
    !
    implicit none

    ! Imported parameters
    type(Taerosol), intent(in) :: aerosol

    fu_gas_phase_idx_persist_org = fu_n_modes(aerosol) + 1

  end function fu_gas_phase_idx_persist_org


  !************************************************************************

  subroutine set_ModeMap_POP(chName, aerosol, iGasPhase, nModesTotSubst, iModeMapSubst)
    !
    ! Selects the correct mapping of size classes and phases for the species.
    ! Pointers must be already allocated but the cocktail is evidently not yet defined
    !
    implicit none

    ! Imported parameters
    character(len=*), intent(in) :: chName
    type(Taerosol), intent(in) :: aerosol
    integer, intent(in) :: nModesTotSubst, iGasPhase
    integer, dimension(:), intent(out) :: iModeMapSubst   ! (nModes)

    !
    ! Acid basic chemistry assumes that all aerosol species are sub-micron. Let's find such mode
    !
    if(fu_n_modes(aerosol) < 1)then
      call report(aerosol)
      call set_error('Strange aerosol','set_ModeMap_POP')
      return
    endif

    iModeMapSubst(1) = 1           ! the finest aerosol mode
    iModeMapSubst(2) = iGasPhase

  end subroutine set_ModeMap_POP


  !***********************************************************************

  subroutine set_chem_rules_POP(nlSetup, rulesPOP)
    !
    ! Creates the set of rules, in accordance with which the computations are to
    ! be made. Input - SILAM namelist
    ! Taking this chance, let's also set once-and-forever a few chemical coefficients
    !
    implicit none

    ! Imported parameter
    type(Tsilam_namelist), pointer :: nlSetup 
    type(Tchem_rules_POP) :: rulesPOP

    ! Local variables
    integer :: iTmp

    !
    ! Stupidity checking
    !
    rulesPOP%defined = silja_false
    if(.not.defined(nlSetup))then
      call set_error('Undefined namelist given','set_chem_rules_POP')
      return
    endif

    rulesPOP%nPrecFlds = fu_content_int(nlSetup,'number_of_precipitation_fields')

    !
    ! For Eulerian scheme we will need the minimum threshold of amount of species in the grid cell
    ! The cell will not be processed by ANY Eulerian-pool routine should it have smaller amount 
    ! than that threshold. A data-based setting is possible only when the total emitted amount
    ! is known, as well as a grid size. So far let's put it to something small enough.
    !
    rulesPOP%massLowThreshold = 1.0e-3

    rulesPOP%defined = silja_true

  end subroutine set_chem_rules_POP


  !***********************************************************************

  subroutine set_missing_chem_rules_POP(rulesPOP)
    implicit none
    type(Tchem_rules_POP), intent(out) :: rulesPOP

    rulesPOP%defined = silja_false
    rulesPOP%nPrecFlds = int_missing
    rulesPOP%massLowThreshold = real_missing
  end subroutine set_missing_chem_rules_POP
  

  !**********************************************************************

  subroutine set_low_mass_thrsh_POP(rules, fThreshold)
    !
    ! Sets the minimum mass per grid cell for Eulerian scheme. Cells with smaller mass
    ! are skipped from the computations
    !
    implicit none

    ! Imported parameters
    type(Tchem_rules_POP), intent(inout) :: rules
    real, intent(in) :: fThreshold

    rules%massLowThreshold = fThreshold

  end subroutine set_low_mass_thrsh_POP


  !********************************************************************************

  real function fu_low_mass_threshold_POP(rules)
    !
    ! Returns the low-mass threshold stored in some specific chemical rules
    !
    implicit none

    type(Tchem_rules_POP), intent(in) :: rules

    fu_low_mass_threshold_POP = rules%massLowThreshold

  end function fu_low_mass_threshold_POP


  !********************************************************************************

  real function fu_zero_mass_threshold_POP(rules)
    !
    ! Returns the zero-mass threshold below which the mass is considered to be
    ! practically zero
    !
    implicit none

    type(Tchem_rules_POP), intent(in) :: rules

    fu_zero_mass_threshold_POP = 1.0e-20

  end function fu_zero_mass_threshold_POP


  !********************************************************************************

  function fu_lifetime_POP(rules, cocktail_template)result(lifetime)
    !
    ! A typical life time due to degradation and deposition for POPs
    ! So far, whatever the species are, we set the "typical" lifetime to be two days
    !
    implicit none

    type(silja_interval) :: lifetime

	type(silam_cocktail), intent(in) :: cocktail_template
    type(Tchem_rules_POP), intent(in) :: rules

    ! Internal variables
	integer :: iSubst, iTmp
    real :: Ea_NPT, alpha, beta, fK

    !
	! There can be several POPs, we take the longest-living and give its lifetime for
	! typical conditions
    !
    do iTmp = 1, nPOPTranspSubst
      iSubst = indPOPTranspSubst(iTmp)
      call get_reaction_parameters(cocktail_template%materials(iSubst)%ptr,'OH','',Ea_NPT,alpha,beta)
      fK = beta * 298.**alpha * exp(-Ea_NPT/(gas_constant_dryair * 298.)) * &
         & 2.0e-12   ! typical OH concentration mole/m3: EMEP night 1.5e-14, day max 5e-12
      if(fu_sec(lifetime) < 1./fK) lifetime = fu_set_interval_sec(1./fK)
	end do

    call msg('Life time due to degradation of POP, hrs:',fu_hour(lifetime))

  end function fu_lifetime_POP


  !***********************************************************************

  subroutine prepare_POP_transdep(met_buf, disp_buf, RulesPOP, aerosolRules, &
                                & interpCoefMet2DispHoriz, ifHorizInterp)
    !
    ! Prepares to computation of the POP deposition and transofrmation
    ! Dry deposition follows resistive analogy, while the wet deposition for aerosols
    ! goes as it goes, gases take their scaveniging coefficients.
    !
    implicit none

    ! Imported parameters
    TYPE(Tfield_buffer), POINTER :: met_buf, disp_buf
    type(Tchem_rules_POP), intent(inout) :: RulesPOP
    type(Taerosol_rules), intent(in) :: aerosolRules
    type(THorizInterpStruct), pointer :: interpCoefMet2DispHoriz
    logical, intent(in) :: ifHorizInterp

    ! Local declarations
    TYPE(silja_field_id) :: idTmp
    INTEGER :: iQ, ix, iy, iDisp, iMeteo
    INTEGER, DIMENSION(:), POINTER :: mdl_in_q
    real, DIMENSION(:), POINTER :: ptrFrOfLand, ptrLsPrecip, ptrCnvPrecip
    real :: fPrec
    TYPE(Tfield_buffer), POINTER :: met_bufPtr, disp_bufPtr
    type(silja_time) :: now

    mdl_in_q => met_buf%buffer_quantities
    disp_bufPtr => disp_buf
    met_bufPtr => met_buf

    !
    ! For safety have to prepare the aerosol processes as well
    !
    call prepare_aerosol_transdep(met_buf, aerosolRules)
    if(error)return

    !
    ! For conversion - 3D temperature is enough. 
    !
    iQ = fu_index(mdl_in_q, temperature_flag)
    if(iQ <= 0)then
      call set_error('No temperature','prepare_POP_transdep')
      return
    endif
    ptrTempr => met_buf%p4d(iQ)

    ! Now - 2m temperature
    !
    iQ = fu_index(mdl_in_q, temperature_2m_flag)
    if(iQ <= 0)then
      call set_error('No 2m temperature','prepare_POP_transdep')
      return
    endif
    ptrTempr2m => met_buf%p2d(iQ)%present%ptr
    now = fu_valid_time(met_buf%p2d(iQ)%present%IdPtr)

    ! Now - cloud cover
    !
    iQ = fu_index(mdl_in_q, total_cloud_cover_flag)
    if(iQ <= 0)then
      call set_error('No cloud cover','prepare_POP_transdep')
      return
    endif
    ptrCloudCover => met_buf%p2d(iQ)%present%ptr

    ! Now - friction velocity
    !
    iQ = fu_index(mdl_in_q, friction_velocity_flag)
    if(iQ <= 0)then
      call set_error('No friction velocity','prepare_POP_transdep')
      return
    endif
    ptrFricVel => met_buf%p2d(iQ)%present%ptr

    !
    ! Deposition needs also precipitation field(s)
    !
    ! First, precipitation field(s)
    !
    iQ = fu_index(mdl_in_q, large_scale_rain_int_flag)
    if(iQ <= 0)then
      call set_error('No large-scale precipitation','prepare_POP_transdep')
      return
    endif
    ptrLsPrecip => met_buf%p2d(iQ)%present%ptr

    if(RulesPOP%nPrecFlds == 2)then
      iQ = fu_index(mdl_in_q, convective_rain_int_flag)
      if(iQ <= 0)then
        call set_error('No convective precipitation','prepare_POP_transdep')
        call unset_error('prepare_POP_transdep')
        RulesPOP%nPrecFlds = 1
      else
        ptrCnvPrecip => met_buf%p2d(iQ)%present%ptr
      endif
    endif

    !
    ! Standard SILAM scavenging
    !
    iQ = fu_index(mdl_in_q, scavenging_coefficient_flag)
    if(iQ <= 0)then
      call set_error('No default scavenging coefficient','prepare_POP_transdep')
      return
    endif
    ptrScavStd => met_buf%p4d(iQ) ! 3D coefficient

    !
    ! Store the dependence on relative humidity switcher
    !
    if(aerosolRules%ifHumidityDependent)then
      iQ = fu_index(mdl_in_q, relative_humidity_flag)
      if(iQ < 0)then
        call set_error('No relative humidity','prepare_POP_transdep')
        return
      endif
      ptrRelHumidity => met_buf%p4d(iQ)
    endif

    !
    ! Finally, sea-land mask
    !
    ptrFrOfLand => fu_grid_data(fraction_of_land_fld)
    if(error)return

  end subroutine prepare_POP_transdep



  !***********************************************************************

  subroutine transform_POP_single_val(mass_vector_1d, now, timestep_sec, weight_past, &
                                    & index_horiz_meteo, index_vert_meteo, &
                                    & rules_POP, cocktail_template)
    !
    ! Implements equilibrium between aerosol and gas phases and degradation due to OH
    !
    implicit none

    ! Imported parameters
    integer, intent(in) :: index_horiz_meteo, index_vert_meteo
    type(silja_time), intent(in) :: now
    real, intent(in) :: timestep_sec, weight_past
    type(Tchem_rules_POP), intent(in) :: rules_POP
    real, dimension(:), intent(inout) :: mass_vector_1d
    type(silam_cocktail), intent(in) :: cocktail_template

    ! Local variables
    real :: tempr, dSO2, fTmp, fOH, fK_OH
    real, dimension(:), pointer :: lonPtr, latPtr
    integer :: iSubst, iMode, iTmp

    ! Get temperature from the meteo buffer
    !
    tempr = weight_past * ptrTempr%past%p2d(index_vert_meteo)%ptr(index_horiz_meteo) + &
          & (1.-weight_past) * ptrTempr%future%p2d(index_vert_meteo)%ptr(index_horiz_meteo)
    if(tempr > 399. .or. tempr < 200.)then
      call msg_warning('Strange temperature:','transform_POP_single_val')
      call msg('Strange temperature at the level:',index_vert_meteo,tempr)
      call set_error('Strange temperature','transform_POP_single_val')
      return
    endif
    !
    ! Compute the actual transformation rate. Unit = moles, so no problem with mass
    !
    lonPtr => fu_grid_data(meteo_longitude_fld)
    latPtr => fu_grid_data(meteo_latitude_fld)

    fOH = max(1.0e4, 1.0e4+4.0e6 * exp(-0.25 / fu_solar_zenith_angle_cos(lonPtr(index_horiz_meteo), &
                                                                       & latPtr(index_horiz_meteo), &
                                                                       & fu_julian_date(now), &
                                                                       & fu_year(now), &
                                                                       & fu_mon(now), &
                                                                       & fu_day(now), &
                                                                       & fu_hour(now), &
                                                                       & fu_min(now), &
                                                                       & int(fu_sec(now))))) * &
        & (1. - ptrCloudCover(index_horiz_meteo) * 0.5)
    if(error)return

    do iTmp = 1, nPOPTranspSubst
      iSubst = indPOPTranspSubst(iTmp)
      do iMode = 1, cocktail_template%nModesTotSubst(iSubst)
        !
        ! Second order non-biased implicit
        !
        fTMp = ptrDegradation_rate(iSubst,int(tempr-198.5)) * fOH * timestep_sec
        fTmp = fTmp + fTmp*fTmp * 0.5
        fTmp = mass_vector_1d(cocktail_template%iSpecies(iSubst,iMode,1)) * fTmp/(1.+fTmp)
        mass_vector_1d(cocktail_template%iSpecies(iSubst,iMode,1)) = &
                        & mass_vector_1d(cocktail_template%iSpecies(iSubst,iMode,1)) - fTmp

      end do !  iMode
    end do ! iSubst
      
    !
    ! Gas-particle partitioning
    !
    call set_error('Gas-particle partitioning does not work yet','transform_POP_single_val')
    return

  end subroutine transform_POP_single_val


  !***********************************************************************

  subroutine make_dry_deposition_POP(src, dest, timestep_sec, weight_past, index_horiz, lyr_thickness, &
                                   & RulesPOP, aerosolRules, cocktail_template)
    !
    ! Dry deposition is treated in accordance with the resistance rules.
    !
    implicit none

    ! Imported parameters
    real, dimension(:), pointer :: src, dest
    integer, intent(in) :: index_horiz
    real, intent(in) :: lyr_thickness, timestep_sec, weight_past
    type(Tchem_rules_POP), intent(in) :: rulesPOP
    type(Taerosol_rules), intent(in) :: aerosolRules
    type(silam_cocktail), pointer :: cocktail_template

    ! Local variables
    real, dimension(:), pointer :: srcPtr, destPtr
    type(TwetParticle) :: wetParticle
    integer :: iSubst, iTmp

    srcPtr => src
    destPtr => dest
    !
    ! Make the dry deposition for each substance one-by-one. But, if there is no humidity dependence,
    ! can make the wetParticle once
    !
    do iTmp = 1, nPOPTranspSubst
      iSubst = indPOPTranspSubst(iTmp)

      if(aerosolRules%ifHumidityDependent)then
        call make_dry_deposition_aerosol(srcPtr, destPtr, &
                        & iSubst, cocktail_template%nModesTotSubst, &
                        & fu_wet_particle_features( &
                                   & cocktail_template%materials(iSubst)%ptr, &
                                      & ptrRelHumidity%past%p2d(1)%ptr(index_horiz)* weight_past + &
                                      & ptrRelHumidity%future%p2d(1)%ptr(index_horiz)* (1.-weight_past)), &
                        & timestep_sec, index_horiz, lyr_thickness, &
                        & cocktail_template%aerosol, aerosolRules, &
                        & cocktail_template%iSpecies, &
                        & cocktail_template%iModeMap)
      else
        call make_dry_deposition_aerosol(srcPtr, destPtr, &
                        & iSubst, cocktail_template%nModesTotSubst, &
                        & fu_wet_particle_features(cocktail_template%materials(iSubst)%ptr, &
                                                 & aerosolRules%fDefaultRelHumidity), &
                        & timestep_sec, index_horiz, lyr_thickness, &
                        & cocktail_template%aerosol, aerosolRules, &
                        & cocktail_template%iSpecies, &
                        & cocktail_template%iModeMap)
      endif
    end do  ! iSubst
    !
    ! Dry deposition for gas phase: so far put to zero.
    ! Later, if ever needed for Lagrangian dynamics, can be made through resistive scheme
    !
  end subroutine make_dry_deposition_POP


  !*****************************************************************************

  subroutine get_Rb_persistent_organics(cocktail_template, indexMeteo, weight_past, &
                                      & RulesPOP, rulesAerosol, ptrRb)
    !
    ! Returns the Rb 
    !
    implicit none

    ! Imported parameters
    type(silam_cocktail), intent(in) :: cocktail_template
    integer, intent(in) :: indexMeteo
    type(Taerosol_rules), intent(in) :: rulesAerosol
    type(Tchem_rules_POP), intent(in) :: RulesPOP
    type(Tmasses), intent(out) :: ptrRb
    real, intent(in) :: weight_past

    ! Local variables
    integer :: iSubst, iMode, iTmp
    type(TwetParticle) :: wetParticle

    do iTmp = 1, nPOPTranspSubst
      iSubst = indPOPTranspSubst(iTmp)
      if(rulesAerosol%ifHumidityDependent)then
        wetParticle = fu_wet_particle_features(cocktail_template%materials(iSubst)%ptr, &
                                             & ptrRelHumidity%past%p2d(1)%ptr(indexMeteo) * weight_past + &
                                             & ptrRelHumidity%future%p2d(1)%ptr(indexMeteo) * (1.-weight_past))
      else
        wetParticle = fu_wet_particle_features(cocktail_template%materials(iSubst)%ptr, &
                                            & rulesAerosol%fDefaultRelHumidity)
      endif
      if(error)return
      do iMode = 1, cocktail_template%nModesTotSubst(iSubst)
        if(cocktail_template%iModeMap(iSubst,iMode) == cocktail_template%iGasPhase)then
          !
          ! Rb = 5 * (Sc)**2/3 / u*
          ! Sc=nu/D, typical nu=1.5e-5 m2/s, D = kT Cun/(3pi d mu) - molecular diffusivity, from tables
          !
          ptrRb%subst(iSubst)%mass(cocktail_template%iModeMap(iSubst,iMode),1) = &
              & 5. / ptrFricVel(indexMeteo) * &
              & (fu_kinematic_viscosity(ptrTempr2m(indexMeteo), density_air_288K) / &
               & fu_molecular_diffusivity_air(cocktail_template%materials(iSubst)%ptr))**0.66666667
        else
          ptrRb%subst(iSubst)%mass(cocktail_template%iModeMap(iSubst,iMode),1) = &
                               & fu_R_b_aerosol(cocktail_template%aerosol, &
                                              & cocktail_template%iModeMap(iSubst,iMode), &
                                              & wetParticle, &
                                              & indexMeteo, &
                                              & rulesAerosol)
        endif
      end do
    enddo

  end subroutine get_Rb_persistent_organics


  !*************************************************************************

  subroutine get_Rs_persistent_organics_crd(cocktail_template, Rules_POP, ptrRs, ifRs_meteo_depend)
    !
    ! Returns the surface resistances for the POP cocktail components
    ! For aerosols zero surface resistance, while gases have Rs dependent 
    ! on meteorology and surface type.
    !
    implicit none

    ! Imported parameters
    type(silam_cocktail), intent(in) :: cocktail_template
    type(Tmasses), intent(out) :: ptrRs
    logical, intent(out) :: ifRs_meteo_depend
    type(Tchem_rules_POP), intent(in) :: rules_POP

    ! Local variables
    integer :: iSubst, iMode, iTmp

    !
    ! For crude Rs, take the default parameters
    !
    do iTmp = 1, nPOPTranspSubst
      iSubst = indPOPTranspSubst(iTmp)
      do iMode = 1, fu_n_modes(cocktail_template, iSubst)
        if(cocktail_template%iModeMap(iSubst,iMode) == cocktail_template%iGasPhase)then
          ptrRs%subst(iSubst)%mass(iMode,1) = &
                       & 0.5 * (fu_Rs_default_soil(cocktail_template%materials(iSubst)%ptr) + &
                              & fu_Rs_default_water(cocktail_template%materials(iSubst)%ptr))
        else
          ptrRs%subst(iSubst)%mass(iMode,1) = 0.
        endif
      enddo
    end do
    ifRs_meteo_depend = .true.

  end subroutine get_Rs_persistent_organics_crd


  !***********************************************************************

  subroutine refine_Rs_persistent_organics(cocktail_template, Rules_POP, indexMeteo, indexDisp, ptrRs)
    !
    ! Surface resistance for POP cocktail is meteo-dependent, so we have to 
    ! refine it deeply within the cycle over grid. This routine makes this trick
    !
    implicit none

    ! Imported parameters
    type(silam_cocktail), intent(in) :: cocktail_template
    integer, intent(in) :: indexMeteo, indexDisp
    type(Tmasses), intent(out) :: ptrRs
    type(Tchem_rules_POP), intent(in) :: rules_POP

    ! Local variables
    real, dimension(:), pointer :: ptrFractionOfLand
    integer :: iSubst, iMode, iTmp

    ptrFractionOfLand => fu_grid_data(fraction_of_land_fld)

    do iTmp = 1, nPOPTranspSubst
      iSubst = indPOPTranspSubst(iTmp)
      do iMode = 1, fu_n_modes(cocktail_template, iSubst)
        if(cocktail_template%iModeMap(iSubst,iMode) == cocktail_template%iGasPhase)then
          ptrRs%subst(iSubst)%mass(iMode,1) = &
                 & ptrFractionOfLand(indexMeteo) * &
                                 & fu_Rs_default_soil(cocktail_template%materials(iSubst)%ptr) + &
                 & (1.-ptrFractionOfLand(indexMeteo)) * &
                                 & fu_Rs_default_water(cocktail_template%materials(iSubst)%ptr) 
        endif
      enddo
    end do

  end subroutine refine_Rs_persistent_organics


  !***********************************************************************

  subroutine scavenge_POP_single_val(src, dest, timestep_sec, weight_past, &
                                   & index_horiz, index_vert, &
                                   & rulesPOP, cocktail_template)
    !
    ! Computes wet deposition of the SO2 and SO4. Since the Galperin's scheme of
    ! scavenging is (i) non-linear and (ii) involves a ratio of concentrations in 
    ! air and in precipitation, it cannot be used here. Solutions are: (i) linearise
    ! Galperin's scheme and somehow translate it to masses from concentrations, 
    ! (ii) use SILAM standard scavenging.
    ! Linearization and conc->mass conversion are not very straightforward and will
    ! eventually lead to somewhat outdated scheme without its essence of saturation,
    ! sub- and in-cloud split, etc. Therefore, it seems to be better to utilize
    ! the SILAM standard scheme for SO4 and then scale it by a factor of ... for SO2.
    ! This is what is done below.
    ! The method employed is: (m(t+dt)-m(t))/dt = -L*0.5*(m(t+dt)+m(t))
    ! It seems to give the second-order accuracy.
    ! The other way is to put exponent up to the second-order Taylor polynome but
    ! it is somewhat longer and involves more strange stuff
    !
    implicit none

    ! Imported parameters
    integer, intent(in) :: index_horiz, index_vert
    real, intent(in) :: weight_past, timestep_sec
    type(Tchem_rules_POP), intent(in) :: rulesPOP
    real, dimension(:), intent(inout) :: src, dest
    type(silam_cocktail), intent(in) :: cocktail_template

    ! Local variables
    real :: scav_coef_std, scav_coef, ScavAmount
    integer :: iSubst, iMode, iTmp
    !
    ! Gases will follow the scaling from standard scavenging. Aerosols just accept it
    ! The initial coefficient is quite comprehensive: it is 3D, rain/snow, etc.
    !
    scav_coef_std = weight_past * ptrScavStd%past%p2d(index_vert)%ptr(index_horiz) + &
              & (1.-weight_past) * ptrScavStd%future%p2d(index_vert)%ptr(index_horiz) * timestep_sec  

    do iTmp = 1, nPOPTranspSubst
      iSubst = indPOPTranspSubst(iTmp)
      do iMode = 1, fu_n_modes(cocktail_template, iSubst)

        if(cocktail_template%iModeMap(iSubst,iMode) == cocktail_template%iGasPhase)then
          scav_coef = scav_coef_std * fu_scavenging_scale(cocktail_template%materials(iSubst)%ptr)
        else
          scav_coef = scav_coef_std
        endif

        scav_coef = scav_coef * (1. + scav_coef * (0.5 + scav_coef / 3.))

        ScavAmount = src(cocktail_template%iSpecies(iSubst,iMode,1)) * scav_coef/(1.+scav_coef)
        src(cocktail_template%iSpecies(iSubst,iMode,1)) = &
                         & src(cocktail_template%iSpecies(iSubst,iMode,1)) - ScavAmount
        dest(cocktail_template%iSpecies(iSubst,iMode,1)) = &
                         & dest(cocktail_template%iSpecies(iSubst,iMode,1)) + ScavAmount
      end do
    end do        

  end subroutine scavenge_POP_single_val



  !****************************************************************

  subroutine scavenge_POP_map(mapConc, mapWetDep, &
                            & timestep_sec, weight_past, &
                            & rulesPOP, &
                            & interpCoefMet2DispHoriz, interpCoefMet2DispVert, &
                            & ifInterpMet2DispHoriz, ifInterpMet2DispVert)
    !
    ! Implements self-degradation of a passive substance.
    ! Handles the whole 3D map of concentrations
    !
    implicit none

    ! Imported parameters
    type(TMass_map), pointer :: mapConc, mapWetDep
    real, intent(in) :: timestep_sec, weight_past
    type(Tchem_rules_POP), intent(in) :: rulesPOP
    type(THorizInterpStruct), pointer :: interpCoefMet2DispHoriz
    type(TVertInterpStruct), pointer :: interpCoefMet2DispVert
    logical, intent(in) :: ifInterpMet2DispHoriz, ifInterpMet2DispVert

    ! Local variables
    real :: scav_coef_std, scav_coef, scavAmount, fTmp
    integer :: ix, iy, iLev, iSrc, iSubst, iMode, iTmp

    do iSrc = 1, mapConc%nSrc
      do iLev = 1, mapConc%n3D
        do iy = 1, mapConc%ny
          do ix = 1, mapConc%nx
            !
            ! Mass check
            !
            if(sum(mapConc%arM(1:mapConc%nSpecies,iSrc,iLev,ix,iy)) > rulesPOP%massLowThreshold)then

              scav_coef_std = fu_get_value(ptrScavStd, &  ! field_ptr_4d
                                         & nx_meteo, &  ! nxFrom
                                         & ix, iy, iLev, & 
                                         & weight_past, &
                                         & interpCoefMet2DispHoriz, interpCoefMet2DispVert, &
                                         & ifInterpMet2DispHoriz, ifInterpMet2DispVert)


              do iTmp = 1, nPOPTranspSubst
                iSubst = indPOPTranspSubst(iTmp)
                do iMode = 1, fu_n_modes(mapConc%cocktail_template, iSubst)

                  if(mapConc%cocktail_template%iModeMap(iSubst,iMode) == mapConc%cocktail_template%iGasPhase)then
                    scav_coef = scav_coef_std * &
                              & fu_scavenging_scale(mapConc%cocktail_template%materials(iSubst)%ptr)
                  else
                    scav_coef = scav_coef_std
                  endif

                  scav_coef = scav_coef * (1. + scav_coef * (0.5 + scav_coef / 3.))

                  ScavAmount = mapConc%arM(mapConc%iSp(iSubst,iMode,1),iSrc,iLev,ix,iy) * &
                             & scav_coef/(1.+scav_coef)
                  mapConc%arM(mapConc%iSp(iSubst,iMode,1),iSrc,iLev,ix,iy) = &
                                   & mapConc%arM(mapConc%iSp(iSubst,iMode,1),iSrc,iLev,ix,iy) - ScavAmount
                  mapWetDep%arM(mapWetDep%iSp(iSubst,iMode,1),iSrc,1,ix,iy) = &
                                   & mapWetDep%arM(mapWetDep%iSp(iSubst,iMode,1),iSrc,1,ix,iy) + ScavAmount
                end do
              end do        

            endif    ! mass above threshold

          end do   ! ix
        enddo   ! iy 
      enddo   ! iLev
    enddo   ! iSrc

  end subroutine scavenge_POP_map



  !***********************************************************************************
  !
  ! As a speed-up for Eulerian scheme, a few map-based routines
  !
  !***********************************************************************************

  !****************************************************************

  subroutine transform_POP_map(mapConc, mapP, &
                             & now, timestep_sec, weight_past, &
                             & rulesPOP, &
                             & interpCoefMet2DispHoriz, interpCoefMet2DispVert, &
                             & ifInterpMet2DispHoriz, ifInterpMet2DispVert)
    !
    ! Implements POP chemistry for the whole 3D map of concentrations
    !
    implicit none

    ! Imported parameters
    type(TMass_map), pointer :: mapConc, mapP
    type(silja_time), intent(in) :: now
    real, intent(in) :: timestep_sec, weight_past
    type(Tchem_rules_POP), intent(in) :: rulesPOP
    type(THorizInterpStruct), pointer :: interpCoefMet2DispHoriz
    type(TVertInterpStruct), pointer :: interpCoefMet2DispVert
    logical, intent(in) :: ifInterpMet2DispHoriz, ifInterpMet2DispVert

    ! Local variables
    real :: tempr, fTmp, fOH
    integer :: ix, iy, iLev, iSrc, iSubst, iMode, iMetIndex, iTmp
    real, dimension(:), pointer :: lonPtr, latPtr

    do iSrc = 1, mapConc%nSrc
      do iLev = 1, mapConc%n3D
        do iy = 1, mapConc%ny
          do ix = 1, mapConc%nx
            if(sum(mapConc%arM(1:mapConc%nSpecies,iSrc,iLev,ix,iy)) > rulesPOP%massLowThreshold)then
              tempr = fu_get_value(ptrTempr, &  ! field_prt_4d
                                 & nx_meteo, &  ! nxFrom
                                 & ix, iy, iLev, & 
                                 & weight_past, &
                                 & interpCoefMet2DispHoriz, interpCoefMet2DispVert, &
                                 & ifInterpMet2DispHoriz, ifInterpMet2DispVert)
              if(tempr > 400. .or. tempr < 200.)then
                call msg_warning('Strange temperature:','transform_POP_map')
                call msg('Strange temperature at the level:',iLev,tempr)
                call set_error('Strange temperature','transform_POP_map')
                return
              endif
              !
              ! Compute the actual transformation rate. Unit = moles, so no problem with mass 
              !
              lonPtr => fu_grid_data(meteo_longitude_fld)
              latPtr => fu_grid_data(meteo_latitude_fld)
              iMetIndex = fu_grid_index(nx_meteo, ix, iy, interpCoefMet2DispHoriz, ifInterpMet2DispHoriz)

              fOH = max(1.0e4, 1.0e4+4.0e6 * exp(-0.25 / fu_solar_zenith_angle_cos(lonPtr(iMetIndex), &
                                                                       & latPtr(iMetIndex), &
                                                                       & fu_julian_date(now), &
                                                                       & fu_year(now), &
                                                                       & fu_mon(now), &
                                                                       & fu_day(now), &
                                                                       & fu_hour(now), &
                                                                       & fu_min(now), &
                                                                       & int(fu_sec(now))))) * &
                   & (1. - ptrCloudCover(iMetIndex) * 0.5)
              if(error)return

              do iTmp = 1, nPOPTranspSubst
                iSubst = indPOPTranspSubst(iTmp)
                do iMode = 1, mapConc%nModes(iSubst)
                  !
                  !  First-order implicit
!                 mass_vector%subst(iSubst)%mass(iMode) = mass_vector%subst(iSubst)%mass(iMode) / &
!                                  & (1. + ptrDegradation_rate(int(tempr-198.5)) * fOH * timestep_sec)
                  !
                  ! Second order non-biased implicit
                  !
                  fTMp = ptrDegradation_rate(iSubst,int(tempr-198.5)) * fOH * timestep_sec
                  fTmp = fTmp + fTmp*fTmp * 0.5
                  fTmp = mapConc%arM(mapConc%iSp(iSubst,iMode,1),iSrc,iLev,ix,iy) * fTmp/(1.+fTmp)
                  mapConc%arM(mapConc%iSp(iSubst,iMode,1),iSrc,iLev,ix,iy) = &
                                   & mapConc%arM(mapConc%iSp(iSubst,iMode,1),iSrc,iLev,ix,iy) - fTmp

                end do !  iMode
              end do ! iSubst

            endif    ! mass above threshold
          end do   ! ix
        enddo   ! iy 
      enddo   ! iLev
    enddo   ! iSrc

  end subroutine transform_POP_map

  logical function fu_if_tla_required_pop(rules) result(required)
    implicit none
    type(Tchem_rules_pop), intent(in) :: rules
    
    call set_error('TLA not available for POP', 'fu_if_tla_required_pop')

  end function tla_required


END MODULE cocktail_persistent_organics

