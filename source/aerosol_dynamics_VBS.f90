MODULE aer_dyn_VBS

  use cocktail_basic

  implicit none
  
  !
  ! public routines for aer_dyn_VBS
  !
  public init_AerDynVBS
  public transform_AerDynVBS
  public set_rules_AerDynVBS
  public AerDynVBS_input_needs
  public fu_if_tla_required

  public  full_spec_lst_4_AerDynVBS
  public  registerSpecies_4_AerDynVBS
  private fu_if_tla_required_VBS

  interface fu_if_tla_required
     module procedure fu_if_tla_required_VBS
  end interface




  !--------------------------------------------------------------------
  !
  !  Tchem_rules_AeroDynVBS type definition
  !
  type Tchem_rules_AerDynVBS
    private
    type(Taerosol_mode) :: modeSOA
    real, dimension(4) :: Csat_298, dH, MolMass
    
    type(silja_logical) :: defined
  end type Tchem_rules_AerDynVBS
  !
  ! Used for addressing the meteo data
  !
  integer, private, pointer, save :: ind_tempr, ind_rh, ind_pres
  !
  ! Used for addressing the species
  !
  integer, private, save :: iAVB0_gas, iAVB1e0_gas, iAVB1e1_gas, iAVB1e2_gas, iAVB1e3_gas, & ! iAVB1e4_gas, iAVB1e5_gas, iAVB1e6_gas, &
                          & iBVB0_gas, iBVB1e0_gas, iBVB1e1_gas, iBVB1e2_gas, iBVB1e3_gas, &
                          & iAVB0_aer, iAVB1e0_aer, iAVB1e1_aer, iAVB1e2_aer, iAVB1e3_aer, & ! iAVB1e4_aer, iAVB1e5_aer, iAVB1e6_aer, &
                          & iBVB0_aer, iBVB1e0_aer, iBVB1e1_aer, iBVB1e2_aer, iBVB1e3_aer
  !
  ! Label of this module
  !
  integer, parameter, public :: aerosol_dynamics_VBS = 5023


  CONTAINS


  !***********************************************************************

  subroutine init_AerDynVBS(rulesAerDynVBS, species_lst_initial)

    ! Set modes, transport- and short-living species arrays 

    implicit none

    ! Imported parameter
    type(Tchem_rules_AerDynVBS), intent(inout) :: rulesAerDynVBS
    type(silam_species), dimension(:), allocatable, intent(out) :: species_lst_initial

  end subroutine init_AerDynVBS


  !***********************************************************************

  subroutine AerDynVBS_input_needs(meteo_input_local)
    !
    ! Returns input needs for the transformation routines. 
    !
    implicit none

    ! Imported parameters
    type(Tmeteo_input), intent(out), target :: meteo_input_local

    ! Local variables
    integer :: iTmp

    meteo_input_local = meteo_input_empty

    meteo_input_local%quantity(1) = temperature_flag    
    meteo_input_local%q_type(1) = meteo_dynamic_flag
    ind_tempr => meteo_input_local%idx(1)

    meteo_input_local%quantity(2) = relative_humidity_flag
    meteo_input_local%q_type(2) = meteo_dynamic_flag
    ind_rh => meteo_input_local%idx(2)
    
    meteo_input_local%quantity(3) = pressure_flag
    meteo_input_local%q_type(3) = meteo_dynamic_flag
    ind_pres => meteo_input_local%idx(3)
    meteo_input_local%nQuantities = 3

  end subroutine AerDynVBS_input_needs


  !***********************************************************************

  subroutine full_spec_lst_4_AerDynVBS(rulesAerDynVBS, &
                                        & speciesEmis, speciesTransp, speciesSL, speciesAerosol, &
                                        & nSpeciesEmis, nSpeciesTransp, nSpeciesSL, nSpeciesAerosol, &
                                        & iClaimedSpecies)
    ! VBS bins in aerosol phase
    ! Gas phase chemistry owns the gas phase bins
    ! IVOCs can be ignored here, as they mainly stay in gas phase
    
    implicit none

    ! Imported parameters
    type(Tchem_rules_AerDynVBS), intent(inout) :: rulesAerDynVBS
    type(silam_species), dimension(:), pointer :: speciesTransp, speciesEmis, speciesSL, speciesAerosol
    integer, intent(in) :: nSpeciesEmis
    integer, intent(inout) :: nSpeciesTransp, nSpeciesSL, nSpeciesAerosol
    integer, dimension(:), intent(inout) :: iClaimedSpecies

    ! Local variable
    integer :: nSelected, iTmp, iEmis
    type(silam_species) :: speciesTmp
    integer, dimension(:), pointer :: indices

    ! Add species to transport (might be already added from emission, but add_species takes care of duplicates)
    call set_species(speciesTmp, fu_get_material_ptr('AVB0'), rulesAerDynVBS%modeSOA)
    call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
    call set_species(speciesTmp, fu_get_material_ptr('AVB1e0'), rulesAerDynVBS%modeSOA)
    call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
    call set_species(speciesTmp, fu_get_material_ptr('AVB1e1'), rulesAerDynVBS%modeSOA)
    call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
    call set_species(speciesTmp, fu_get_material_ptr('AVB1e2'), rulesAerDynVBS%modeSOA)
    call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
    call set_species(speciesTmp, fu_get_material_ptr('AVB1e3'), rulesAerDynVBS%modeSOA)
    call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
!    call set_species(speciesTmp, fu_get_material_ptr('AVB1e4'), rulesAerDynVBS%modeSOA)
!    call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
!    call set_species(speciesTmp, fu_get_material_ptr('AVB1e5'), rulesAerDynVBS%modeSOA)
!    call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
!    call set_species(speciesTmp, fu_get_material_ptr('AVB1e6'), rulesAerDynVBS%modeSOA)
!    call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
    call set_species(speciesTmp, fu_get_material_ptr('BVB0'), rulesAerDynVBS%modeSOA)
    call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
    call set_species(speciesTmp, fu_get_material_ptr('BVB1e0'), rulesAerDynVBS%modeSOA)
    call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
    call set_species(speciesTmp, fu_get_material_ptr('BVB1e1'), rulesAerDynVBS%modeSOA)
    call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
    call set_species(speciesTmp, fu_get_material_ptr('BVB1e2'), rulesAerDynVBS%modeSOA)
    call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
    call set_species(speciesTmp, fu_get_material_ptr('BVB1e3'), rulesAerDynVBS%modeSOA)
    call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
    
  end subroutine full_spec_lst_4_AerDynVBS


  !*******************************************************************

  subroutine registerSpecies_4_AerDynVBS(rulesAerDynVBS, speciesTrans, speciesSL, &
                                           & nSpeciesTrans, nSpeciesSL)
    !
    ! Set the species indices arrays
    !
    implicit none
    
    ! Imported parameters
    type(Tchem_rules_AerDynVBS), intent(in) :: rulesAerDynVBS
    type(silam_species), dimension(:), intent(in) :: speciesTrans, speciesSL
    integer, intent(in) :: nSpeciesTrans, nSpeciesSL

    ! Local variables
    type(silam_species) :: speciesTmp
    
    call set_species(speciesTmp, fu_get_material_ptr('AVB0'), in_gas_phase)
    iAVB0_gas = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
    call set_species(speciesTmp, fu_get_material_ptr('AVB1e0'), in_gas_phase)
    iAVB1e0_gas = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
    call set_species(speciesTmp, fu_get_material_ptr('AVB1e1'), in_gas_phase)
    iAVB1e1_gas = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
    call set_species(speciesTmp, fu_get_material_ptr('AVB1e2'), in_gas_phase)
    iAVB1e2_gas = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
    call set_species(speciesTmp, fu_get_material_ptr('AVB1e3'), in_gas_phase)
    iAVB1e3_gas = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
!    call set_species(speciesTmp, fu_get_material_ptr('AVB1e4'), in_gas_phase)
!    iAVB1e4_gas = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
!    call set_species(speciesTmp, fu_get_material_ptr('AVB1e5'), in_gas_phase)
!    iAVB1e5_gas = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
!    call set_species(speciesTmp, fu_get_material_ptr('AVB1e6'), in_gas_phase)
!    iAVB1e6_gas = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
    call set_species(speciesTmp, fu_get_material_ptr('BVB0'), in_gas_phase)
    iBVB0_gas = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
    call set_species(speciesTmp, fu_get_material_ptr('BVB1e0'), in_gas_phase)
    iBVB1e0_gas = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
    call set_species(speciesTmp, fu_get_material_ptr('BVB1e1'), in_gas_phase)
    iBVB1e1_gas = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
    call set_species(speciesTmp, fu_get_material_ptr('BVB1e2'), in_gas_phase)
    iBVB1e2_gas = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
    call set_species(speciesTmp, fu_get_material_ptr('BVB1e3'), in_gas_phase)
    iBVB1e3_gas = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
    call set_species(speciesTmp, fu_get_material_ptr('AVB0'), rulesAerDynVBS%modeSOA)
    iAVB0_aer = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
    call set_species(speciesTmp, fu_get_material_ptr('AVB1e0'), rulesAerDynVBS%modeSOA)
    iAVB1e0_aer = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
    call set_species(speciesTmp, fu_get_material_ptr('AVB1e1'), rulesAerDynVBS%modeSOA)
    iAVB1e1_aer = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
    call set_species(speciesTmp, fu_get_material_ptr('AVB1e2'), rulesAerDynVBS%modeSOA)
    iAVB1e2_aer = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
    call set_species(speciesTmp, fu_get_material_ptr('AVB1e3'), rulesAerDynVBS%modeSOA)
    iAVB1e3_aer = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
!    call set_species(speciesTmp, fu_get_material_ptr('AVB1e4'), rulesAerDynVBS%modeSOA)
!    iAVB1e4_aer = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
!    call set_species(speciesTmp, fu_get_material_ptr('AVB1e5'), rulesAerDynVBS%modeSOA)
!    iAVB1e5_aer = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
!    call set_species(speciesTmp, fu_get_material_ptr('AVB1e6'), rulesAerDynVBS%modeSOA)
!    iAVB1e6_aer = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
    call set_species(speciesTmp, fu_get_material_ptr('BVB0'), rulesAerDynVBS%modeSOA)
    iBVB0_aer = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
    call set_species(speciesTmp, fu_get_material_ptr('BVB1e0'), rulesAerDynVBS%modeSOA)
    iBVB1e0_aer = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
    call set_species(speciesTmp, fu_get_material_ptr('BVB1e1'), rulesAerDynVBS%modeSOA)
    iBVB1e1_aer = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
    call set_species(speciesTmp, fu_get_material_ptr('BVB1e2'), rulesAerDynVBS%modeSOA)
    iBVB1e2_aer = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
    call set_species(speciesTmp, fu_get_material_ptr('BVB1e3'), rulesAerDynVBS%modeSOA)
    iBVB1e3_aer = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)

    if (any(  (/iAVB0_gas, iAVB1e0_gas, iAVB1e1_gas, iAVB1e2_gas, iAVB1e3_gas, & ! iAVB1e4_gas, iAVB1e5_gas, iAVB1e6_gas, &
              & iBVB0_gas, iBVB1e0_gas, iBVB1e1_gas, iBVB1e2_gas, iBVB1e3_gas, &
              & iAVB0_aer, iAVB1e0_aer, iAVB1e1_aer, iAVB1e2_aer, iAVB1e3_aer, & ! iAVB1e4_aer, iAVB1e5_aer, iAVB1e6_aer, &
              & iBVB0_aer, iBVB1e0_aer, iBVB1e1_aer, iBVB1e2_aer, iBVB1e3_aer/) < 1)) then
          call report(speciesTrans)
          call set_error("Clould not find all transport species for VBS", &
                    & "registerSpecies_4_AerDynVBS")
   endif
              
              

  end subroutine registerSpecies_4_AerDynVBS

  !*******************************************************************

  subroutine transform_AerDynVBS(vMassTrn, vMassSL, garbage, rulesAerDynVBS, fLowCncThresh, &
                                  & metdat, seconds, print_it)
    !
    ! Gas-Particle partitioning
    !
    implicit none

    ! Imported parameters
    real, dimension(:), intent(inout) :: vMassTrn, vMassSL
    type(Tchem_rules_AerDynVBS), intent(in) :: rulesAerDynVBS
    real, dimension(:), intent(inout) :: garbage
    real, dimension(:), intent(in) :: metdat, fLowCncThresh
    real, intent(in) :: seconds
    logical, intent(out) :: print_it

    ! Local variables
    integer :: iBin, nIter
    real :: cOAer, cOGas, cOAer_new, OAchange, cnc_bin
    real, dimension(4) :: c_sat, aer_fract
    !type(TwetParticle) :: wetParticle
    
    print_it = .false.  ! set to true and chemistry manager will give complete dump for this cell

    ! Move the nonvolatile stuff that chemistry created to aerosol bin
    vMassTrn(iAVB0_aer) = vMassTrn(iAVB0_aer) + vMassTrn(iAVB0_gas)
    vMassTrn(iAVB0_gas) = 0.0
    vMassTrn(iBVB0_aer) = vMassTrn(iBVB0_aer) + vMassTrn(iBVB0_gas)
    vMassTrn(iBVB0_gas) = 0.0 
    
    ! Totals
    cOAer =  vMassTrn(iAVB0_aer)+vMassTrn(iAVB1e0_aer)+vMassTrn(iAVB1e1_aer)+ &
         & vMassTrn(iAVB1e2_aer)+vMassTrn(iAVB1e3_aer)+ &
         !& vMassTrn(iAVB1e4_aer)+ vMassTrn(iAVB1e5_aer)+vMassTrn(iAVB1e6_aer)+ &
         & vMassTrn(iBVB0_aer)+vMassTrn(iBVB1e0_aer)+vMassTrn(iBVB1e1_aer)+ &
         & vMassTrn(iBVB1e2_aer)+vMassTrn(iBVB1e3_aer)
         
    cOGas =  vMassTrn(iAVB0_gas)+vMassTrn(iAVB1e0_gas)+vMassTrn(iAVB1e1_gas)+ &
         & vMassTrn(iAVB1e2_gas)+vMassTrn(iAVB1e3_gas)+ &
         !& vMassTrn(iAVB1e4_gas)+ vMassTrn(iAVB1e5_gas)+vMassTrn(iAVB1e6_gas)+ &
         & vMassTrn(iBVB0_gas)+vMassTrn(iBVB1e0_gas)+vMassTrn(iBVB1e1_gas)+ &
         & vMassTrn(iBVB1e2_gas)+vMassTrn(iBVB1e3_gas)
    
!    call msg('Beginning cOAer, cOGas', cOAer, cOGas)     
    if(cOAer+cOGas < fLowCncThresh(iAVB1e0_aer))then ! threshold of whichever bin 
        garbage(iAVB0_gas) = garbage(iAVB0_gas) + vMassTrn(iAVB0_gas)
        garbage(iAVB1e0_gas) = garbage(iAVB1e0_gas) + vMassTrn(iAVB1e0_gas)
        garbage(iAVB1e1_gas) = garbage(iAVB1e1_gas) + vMassTrn(iAVB1e1_gas)
        garbage(iAVB1e2_gas) = garbage(iAVB1e2_gas) + vMassTrn(iAVB1e2_gas)
        garbage(iAVB1e3_gas) = garbage(iAVB1e3_gas) + vMassTrn(iAVB1e3_gas)
        !garbage(iAVB1e4_gas) = garbage(iAVB1e4_gas) + vMassTrn(iAVB1e4_gas)
        !garbage(iAVB1e5_gas) = garbage(iAVB1e5_gas) + vMassTrn(iAVB1e5_gas)
        !garbage(iAVB1e6_gas) = garbage(iAVB1e6_gas) + vMassTrn(iAVB1e6_gas)
        garbage(iBVB0_gas) = garbage(iBVB0_gas) + vMassTrn(iBVB0_gas)
        garbage(iBVB1e0_gas) = garbage(iBVB1e0_gas) + vMassTrn(iBVB1e0_gas)
        garbage(iBVB1e1_gas) = garbage(iBVB1e1_gas) + vMassTrn(iBVB1e1_gas)
        garbage(iBVB1e2_gas) = garbage(iBVB1e2_gas) + vMassTrn(iBVB1e2_gas)
        garbage(iBVB1e3_gas) = garbage(iBVB1e3_gas) + vMassTrn(iBVB1e3_gas)
        garbage(iAVB0_aer) = garbage(iAVB0_aer) + vMassTrn(iAVB0_aer)
        garbage(iAVB1e0_aer) = garbage(iAVB1e0_aer) + vMassTrn(iAVB1e0_aer)
        garbage(iAVB1e1_aer) = garbage(iAVB1e1_aer) + vMassTrn(iAVB1e1_aer)
        garbage(iAVB1e2_aer) = garbage(iAVB1e2_aer) + vMassTrn(iAVB1e2_aer)
        garbage(iAVB1e3_aer) = garbage(iAVB1e3_aer) + vMassTrn(iAVB1e3_aer)
        !garbage(iAVB1e4_aer) = garbage(iAVB1e4_aer) + vMassTrn(iAVB1e4_aer)
        !garbage(iAVB1e5_aer) = garbage(iAVB1e5_aer) + vMassTrn(iAVB1e5_aer)
        !garbage(iAVB1e6_aer) = garbage(iAVB1e6_aer) + vMassTrn(iAVB1e6_aer)
        garbage(iBVB0_aer) = garbage(iBVB0_aer) + vMassTrn(iBVB0_aer)
        garbage(iBVB1e0_aer) = garbage(iBVB1e0_aer) + vMassTrn(iBVB1e0_aer)
        garbage(iBVB1e1_aer) = garbage(iBVB1e1_aer) + vMassTrn(iBVB1e1_aer)
        garbage(iBVB1e2_aer) = garbage(iBVB1e2_aer) + vMassTrn(iBVB1e2_aer)
        garbage(iBVB1e3_aer) = garbage(iBVB1e3_aer) + vMassTrn(iBVB1e3_aer)
    
        vMassTrn(iAVB0_gas) = 0.0
        vMassTrn(iAVB1e0_gas) = 0.0
        vMassTrn(iAVB1e1_gas) = 0.0
        vMassTrn(iAVB1e2_gas) = 0.0
        vMassTrn(iAVB1e3_gas) = 0.0
        !vMassTrn(iAVB1e4_gas) = 0.0
        !vMassTrn(iAVB1e5_gas) = 0.0
        !vMassTrn(iAVB1e6_gas) = 0.0
        vMassTrn(iBVB0_gas) = 0.0
        vMassTrn(iBVB1e0_gas) = 0.0
        vMassTrn(iBVB1e1_gas) = 0.0
        vMassTrn(iBVB1e2_gas) = 0.0
        vMassTrn(iBVB1e3_gas) = 0.0
        vMassTrn(iAVB0_aer) = 0.0
        vMassTrn(iAVB1e0_aer) = 0.0
        vMassTrn(iAVB1e1_aer) = 0.0
        vMassTrn(iAVB1e2_aer) = 0.0
        vMassTrn(iAVB1e3_aer) = 0.0
        !vMassTrn(iAVB1e4_aer) = 0.0
        !vMassTrn(iAVB1e5_aer) = 0.0
        !vMassTrn(iAVB1e6_aer) = 0.0
        vMassTrn(iBVB0_aer) = 0.0
        vMassTrn(iBVB1e0_aer) = 0.0
        vMassTrn(iBVB1e1_aer) = 0.0
        vMassTrn(iBVB1e2_aer) = 0.0
        vMassTrn(iBVB1e3_aer) = 0.0
        return
    endif
     
    ! claculate the saturation concentrations
    do iBin = 1,4  !7
        c_sat(iBin) = rulesAerDynVBS%Csat_298(iBin) / metdat(ind_tempr) * 298.0 * &
            & exp(rulesAerDynVBS%dH(iBin) / gas_constant_uni * &
                & (1.0/298.0 - 1.0/metdat(ind_tempr)))
    enddo
    
    ! iterate until aerosol changes less than 5%
    OAchange = 1.0
    nIter = 0
    do while(OAchange > 0.05 .and. nIter <= 5)
      nIter = nIter + 1
      if (cOAer > 0.0)then
        do iBin = 1, 4 !7
          aer_fract(iBin) = 1.0 / (1.0 + c_sat(iBin)/cOAer)
        enddo
      else
!call msg('0 cOAer, cOGas', cOAer, cOGas)
        ! Hope that iterations take care of this
        do iBin = 1, 4  !7
          aer_fract(iBin) = 0.5
        enddo
      endif
      cnc_bin = vMassTrn(iAVB1e0_aer) + vMassTrn(iAVB1e0_gas)
      vMassTrn(iAVB1e0_aer) = cnc_bin * aer_fract(1) 
      vMassTrn(iAVB1e0_gas) = cnc_bin - vMassTrn(iAVB1e0_aer)
      cnc_bin = vMassTrn(iAVB1e1_aer) + vMassTrn(iAVB1e1_gas)
      vMassTrn(iAVB1e1_aer) = cnc_bin * aer_fract(2)
      vMassTrn(iAVB1e1_gas) = cnc_bin - vMassTrn(iAVB1e1_aer)
      cnc_bin = vMassTrn(iAVB1e2_aer) + vMassTrn(iAVB1e2_gas)
      vMassTrn(iAVB1e2_aer) = cnc_bin * aer_fract(3)
      vMassTrn(iAVB1e2_gas) = cnc_bin - vMassTrn(iAVB1e2_aer)
      cnc_bin = vMassTrn(iAVB1e3_aer) + vMassTrn(iAVB1e3_gas)
      vMassTrn(iAVB1e3_aer) = cnc_bin * aer_fract(4)
      vMassTrn(iAVB1e3_gas) = cnc_bin - vMassTrn(iAVB1e3_aer)
!      cnc_bin = vMassTrn(iAVB1e4_aer) + vMassTrn(iAVB1e4_gas)
!      vMassTrn(iAVB1e4_aer) = cnc_bin * aer_fract(5)
!      vMassTrn(iAVB1e4_gas) = cnc_bin - vMassTrn(iAVB1e4_aer)
!      cnc_bin = vMassTrn(iAVB1e5_aer) + vMassTrn(iAVB1e5_gas)
!      vMassTrn(iAVB1e5_aer) = cnc_bin * aer_fract(6)
!      vMassTrn(iAVB1e5_gas) = cnc_bin - vMassTrn(iAVB1e5_aer)
!      cnc_bin = vMassTrn(iAVB1e6_aer) + vMassTrn(iAVB1e6_gas)
!      vMassTrn(iAVB1e6_aer) = cnc_bin * aer_fract(7)
!      vMassTrn(iAVB1e6_gas) = cnc_bin - vMassTrn(iAVB1e6_aer)
      cnc_bin = vMassTrn(iBVB1e0_aer) + vMassTrn(iBVB1e0_gas)
      vMassTrn(iBVB1e0_aer) = cnc_bin * aer_fract(1)
      vMassTrn(iBVB1e0_gas) = cnc_bin - vMassTrn(iBVB1e0_aer)
      cnc_bin = vMassTrn(iBVB1e1_aer) + vMassTrn(iBVB1e1_gas)
      vMassTrn(iBVB1e1_aer) = cnc_bin * aer_fract(2)
      vMassTrn(iBVB1e1_gas) = cnc_bin - vMassTrn(iBVB1e1_aer)
      cnc_bin = vMassTrn(iBVB1e2_aer) + vMassTrn(iBVB1e2_gas)
      vMassTrn(iBVB1e2_aer) =  cnc_bin * aer_fract(3)
      vMassTrn(iBVB1e2_gas) = cnc_bin - vMassTrn(iBVB1e2_aer)
      cnc_bin = vMassTrn(iBVB1e3_aer) + vMassTrn(iBVB1e3_gas)
      vMassTrn(iBVB1e3_aer) = cnc_bin * aer_fract(4)
      vMassTrn(iBVB1e3_gas) = cnc_bin - vMassTrn(iBVB1e3_aer)
      
      cOAer_new = vMassTrn(iAVB0_aer)+vMassTrn(iAVB1e0_aer)+vMassTrn(iAVB1e1_aer)+ &
                & vMassTrn(iAVB1e2_aer)+vMassTrn(iAVB1e3_aer)+ &
                !& vMassTrn(iAVB1e4_aer)+vMassTrn(iAVB1e5_aer)+vMassTrn(iAVB1e6_aer)+ &
                & vMassTrn(iBVB0_aer)+ vMassTrn(iBVB1e0_aer)+vMassTrn(iBVB1e1_aer)+ &
                & vMassTrn(iBVB1e2_aer)+vMassTrn(iBVB1e3_aer)
      if (cOAer > 0.0)then
          OAchange = abs(cOAer_new - cOAer) / cOAer  
      else
          if(cOAer_new > 0.0)then
              OAchange = 1.0
          else
              OAchange = 0.0
          endif
      endif
      cOAer = cOAer_new 
    enddo
!    cOGas =  vMassTrn(iAVB0_gas)+vMassTrn(iAVB1e0_gas)+vMassTrn(iAVB1e1_gas)+ &
!         & vMassTrn(iAVB1e2_gas)+vMassTrn(iAVB1e3_gas)+ &
!         !& vMassTrn(iAVB1e4_gas)+vMassTrn(iAVB1e5_gas)+vMassTrn(iAVB1e6_gas)+ &
!         & vMassTrn(iBVB0_gas)+vMassTrn(iBVB1e0_gas)+vMassTrn(iBVB1e1_gas)+ &
!         & vMassTrn(iBVB1e2_gas)+vMassTrn(iBVB1e3_gas)
!if(nIter == 5)then
!call msg('End cOAer, cOGas', cOAer, cOGas)
!endif
 
  end subroutine transform_AerDynVBS
  !*******************************************************************


  subroutine set_rules_AerDynVBS(nlSetup, rulesAerDynVBS)
    !
    ! Creates the set of rules, in accordance with which the computations are to
    ! be made. Input - SILAM namelist
    !
    implicit none

    ! Imported parameter
    type(Tsilam_namelist), intent (in) :: nlSetup 
    type(Tchem_rules_AerDynVBS), intent(out) :: rulesAerDynVBS

    ! local variables
    rulesAerDynVBS%modeSOA = fu_set_mode(fixed_diameter_flag, 0.1e-6, 1.0e-6, 0.5e-6) !whatever
    ! parameters for the bins from Tsimpidi ea 2010
    ! First 2 bins assumed to be always aerosol (lumped to VB0) and last 3 always gas(do not exist in particulate phase)
    !rulesAerDynVBS%Csat_298 = /0.01, 0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0, 1000000.0/  !ug/m3
    rulesAerDynVBS%Csat_298 = (/1.0e-9, 10.0e-9, 100.0e-9, 1000.0e-9/) !kg/m3
    !rulesAerDynVBS%dH = /112.0, 106.0, 100.0, 94.0, 88.0, 82.0, 76.0, 70.0, 64.0/  !kJ/mol
    rulesAerDynVBS%dH = (/100.0e3, 94.0e3, 88.0e3, 82.0e3/) !J/mol
!    rulesAerDynVBS%MolMass = (/250.0e-3, 250.0e-3, 250.0e-3, 250.0e-3/) !kg/mol
    rulesAerDynVBS%defined = silja_true

  end subroutine set_rules_AerDynVBS
 
  !************************************************************************************

  logical function fu_if_tla_required_VBS(rules) result(required)
    ! Collect transformations' requests for tangent linearization. If tangent linear
    ! is not allowed, corresponding subroutine sets error. 
    implicit none
    type(Tchem_rules_AerDynVBS), intent(out) :: rules
    
    call set_error('Aerosol dynamics VBS does not have TLA/ADJOINT','fu_if_tla_required_VBS')
    required = .false.
    
  end function fu_if_tla_required_VBS

END MODULE aer_dyn_VBS
