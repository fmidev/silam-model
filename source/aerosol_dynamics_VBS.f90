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

  public  full_spec_lst_4_AerDynVBS
  public  registerSpecies_4_AerDynVBS
  public fu_tla_size_VBS

  !--------------------------------------------------------------------
  !
  !  Tchem_rules_AeroDynVBS type definition
  !
  type Tchem_rules_AerDynVBS
    private
    type(Taerosol_mode) :: modeSOA
    real, dimension(0:4) :: Csat_298, dH
    
    type(silja_logical) :: defined
  end type Tchem_rules_AerDynVBS
  !
  ! Used for addressing the meteo data
  !
  integer, private, pointer, save :: ind_tempr, ind_rh, ind_pres
  !
  ! Used for addressing the species
  !
  
  !! Same as matrix 
  integer, dimension(1:2,0:4,1:2) :: iTr !(gas:part, ivMode, anthro:bio) 
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
    integer  :: iAVB0_gas, iAVB1e0_gas, iAVB1e1_gas, iAVB1e2_gas, iAVB1e3_gas, & ! iAVB1e4_gas, iAVB1e5_gas, iAVB1e6_gas, &
              & iAVB0_aer, iAVB1e0_aer, iAVB1e1_aer, iAVB1e2_aer, iAVB1e3_aer, & ! iAVB1e4_aer, iAVB1e5_aer, iAVB1e6_aer, &
              & iBVB0_gas, iBVB1e0_gas, iBVB1e1_gas, iBVB1e2_gas, iBVB1e3_gas, &
              & iBVB0_aer, iBVB1e0_aer, iBVB1e1_aer, iBVB1e2_aer, iBVB1e3_aer
    
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

    !!!  integer, dimension(1:2,0:4,1:2) :: iTr !(gas:part, ivMode, anthro:bio) 
    itr(1,:,1) = (/iAVB0_gas, iAVB1e0_gas, iAVB1e1_gas, iAVB1e2_gas, iAVB1e3_gas/)
    itr(1,:,2) = (/iBVB0_gas, iBVB1e0_gas, iBVB1e1_gas, iBVB1e2_gas, iBVB1e3_gas/)
    itr(2,:,1) = (/iAVB0_aer, iAVB1e0_aer, iAVB1e1_aer, iAVB1e2_aer, iAVB1e3_aer/)
    itr(2,:,2) = (/iBVB0_aer, iBVB1e0_aer, iBVB1e1_aer, iBVB1e2_aer, iBVB1e3_aer/)


    if (any( itr < 1)) then
          call report(speciesTrans)
          call set_error("Clould not find all transport species for VBS", &
                    & "registerSpecies_4_AerDynVBS")
   endif
              
              

  end subroutine registerSpecies_4_AerDynVBS

  !*******************************************************************

  subroutine transform_AerDynVBS(vMassTrn, rulesAerDynVBS, metdat, seconds, vtla)
    !
    ! Gas-Particle partitioning
    !
    implicit none

    ! Imported parameters
    real, dimension(:), intent(inout) :: vMassTrn
    type(Tchem_rules_AerDynVBS), intent(in) :: rulesAerDynVBS
    real, dimension(:), intent(in) :: metdat
    real, dimension(:), pointer :: vtla !!!Single value
    real, intent(in) :: seconds

    ! Local variables
    integer :: iBin, iAB, iaer, igas 
    real :: cOAer, cOGas, cOAer_new, OAchange, cnc_bin
    real :: csat, ctot
    integer :: bitmask 
    !type(TwetParticle) :: wetParticle
    character(len=*), parameter :: sub_name = 'transform_AerDynVBS'
    
     if (associated(vtla)) then
       if (fu_fails(size(vtla) == 1, "Wrong size TLA", sub_name)) return
     endif

    if (seconds > 0) then
      bitmask = 0
      do iBin = 0,4  
         ! claculate the saturation concentration
         csat = rulesAerDynVBS%Csat_298(iBin) / metdat(ind_tempr) * 298.0 * &
              & exp(rulesAerDynVBS%dH(iBin) / gas_constant_uni * &
                  & (1.0/298.0 - 1.0/metdat(ind_tempr)))
         do iAB = 1,2
            !!integer, dimension(1:2,0:4,1:2) :: iTr (gas:part, ivMode, anthro:bio) 
            bitmask = bitmask * 2 ! <<1
            igas = iTr(1,iBin,iAB)
            iaer = iTr(2,iBin,iAB)
            ctot =  vMassTrn(igas) + vMassTrn(iaer)
            bitmask = bitmask * 2 ! <<1
            if (ctot > csat) then                     ! / 0   0 \                        !
               vMassTrn(igas) = csat                  !|         |   TL  for saturation  !
               vMassTrn(iaer) = ctot - csat           ! \ 1   1 /                        !
               bitmask = bitmask + 1
            else
               vMassTrn(igas) = ctot                  ! / 1   1 \                               !
               vMassTrn(iaer) = 0.                    !|         |   TL  for no saturation      !
            endif                                     ! \ 0   0 /                               !
         enddo
      enddo
      if (associated(vtla)) vtla(1) = bitmask !!! WARNING: only 24 bit available
    else 
      if (fu_fails( associated(vtla) .and. size(vtla) == 1, & !!! Adjoint, single value needed
              & "Unassociated TLA in adjoint", sub_name)) return
      !!10 bit allowed for now
      if (fu_fails( vtla(1) >= 0 .and. vtla(1) < 2**11, "Strane TLA in adjoint", sub_name)) return 
      bitmask = int(vtla(1))
      do iBin = 4, 0, -1  !! reverse order to pop corresponding bits from bitmask
         do iAB = 2,1, -1
            !!integer, dimension(1:2,0:4,1:2) :: iTr (gas:part, ivMode, anthro:bio) 
            igas = iTr(1,iBin,iAB)
            iaer = iTr(2,iBin,iAB)
            if (modulo(bitmask,1) == 1) then                           ! / 0   1 \                           !
              !! Saturation, aerosol only could be affected by either  !|         |   TLA  for saturation    !
               vMassTrn(igas) = vMassTrn(iaer)                         ! \ 0   1 /                           !
               !!! vMassTrn(iaer) =  !!!aerosol the same
            else        
              !! No saturation gas only could be affected by either    ! / 1   0 \                            !
              !vMassTrn(igas) = ctot  !!! gas the same                 !|         |   TLA  for no saturation  !
               vMassTrn(iaer) = vMassTrn(igas)                         ! \ 1   0 /                            !
            endif
            bitmask = bitmask / 2 !>>1
         enddo
      enddo
    endif
 
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
    rulesAerDynVBS%Csat_298 = (/0., 1.0e-9, 10.0e-9, 100.0e-9, 1000.0e-9/) !kg/m3
    rulesAerDynVBS%dH = (/0., 100.0e3, 94.0e3, 88.0e3, 82.0e3/) !J/mol
    rulesAerDynVBS%defined = silja_true

  end subroutine set_rules_AerDynVBS
 
  !************************************************************************************

  integer function fu_tla_size_VBS(rules) result(n)
    implicit none
    type(Tchem_rules_AerDynVBS), intent(in) :: rules
    
    n  = 1
    
  end function fu_tla_size_VBS

END MODULE aer_dyn_VBS
