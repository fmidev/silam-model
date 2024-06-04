MODULE aer_dyn_simple

  use cocktail_basic

  implicit none
  
  !
  ! public routines for aer_dyn_simple
  !
  public init_AerDynSimple
  public transform_AerDynSimple
  public set_rules_AerDynSimple
  public full_spec_lst_4_AerDynSimple
  public registerSpecies_4_AerDynSimple
  public AerDynSimple_input_needs
  public fu_tla_size_ADS

  !--------------------------------------------------------------------
  !
  !  Tchem_rules_AeroDynSimple type definition
  !
  type Tchem_rules_AerDynSimple
    private
    type(Taerosol_mode) :: modeCondensation, modeDroplet, modeCoarse
    character(len=substNmLen) :: ssltName4NO3c, ssltName4SO4c
    integer :: nOHaging, maxAgingModes
    character(len=substNmLen), dimension(:), pointer :: OHaging_freshSps, OHaging_agedSps
    real, dimension(:), pointer :: OHaging_rates
    logical :: SSLTmodes4NO3c, SSLTmodes4SO4c
    real :: fParticleSrfRatio, no3_StickingCoeff, so4_StickingCoeff
    type(silja_logical) :: defined
  end type Tchem_rules_AerDynSimple
  !
  ! Used for addressing the meteo data
  !
  integer, private, pointer, save :: ind_tempr=>null(), ind_rh=>null(), &
                                     & ind_pres=>null(), ind_hgt=>null(), ind_cloud_cvr=>null()
  !
  ! Used for addressing the species
  !
  integer, private, save :: iSO4f_aer, iSO4c_aer, iSO4w_aer_sl, iH2SO4_gas_sl, iNH415SO4f_aer, iNH415SO4c_aer, &
                          & iNH3, iHNO3, iNH4NO3_aer, iNO3c_aer
  integer, dimension(:), pointer, private, save :: iSslt_aer => null(), iSsltNO3  => null()
  integer, dimension(:,:), pointer, private, save :: iFresh  => null(), iAged  => null()
  real, dimension(:), pointer, private, save :: sslt_D => null()

  !
  ! Label of this module
  !
  integer, parameter, public :: aerosol_dynamics_simple = 5021

  !Sizes for lookup table
  integer, private, parameter :: NH3_HNO3_vs_NH4NO3_eq_nT = 211, &
                               & NH3_HNO3_vs_NH4NO3_eq_nRH = 22, &
                               & NH3_HNO3_vs_NH4NO3_eq_Tmin = 149.0

  real, dimension(NH3_HNO3_vs_NH4NO3_eq_nRH), parameter :: RH = (/0., 0.50, 0.54, 0.58, 0.62, 0.66, 0.70, &
       & 0.73, 0.76, 0.79, 0.82,  0.85, 0.87, 0.89, 0.91, 0.93, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0/)

  !!!  Lookup table
  real, dimension(NH3_HNO3_vs_NH4NO3_eq_nRH,NH3_HNO3_vs_NH4NO3_eq_nT), &
                                 & private :: NH3_HNO3_vs_NH4NO3_eq = F_NAN

  CONTAINS


  !***********************************************************************

  subroutine init_AerDynSimple(species_lst_initial)

    ! Set modes, transport- and short-living species arrays 
    ! Create  NH3_HNO3_vs_NH4NO3 LUT

    implicit none

    ! Imported parameter
    type(silam_species), dimension(:), allocatable, intent(out) :: species_lst_initial

    real, dimension(7) :: a, b, c
    integer :: nT, iRh, nRh, iT, i
    real:: dryVFract, m, T
    type(TwetParticle) :: wetParticle

    nT = NH3_HNO3_vs_NH4NO3_eq_nT
    nRH = NH3_HNO3_vs_NH4NO3_eq_nRH
    
    ! RH dependent equilibrium from Mozurkewich 1993 
    a(1:7) = (/ -9.293,   33.983,  45.4557, -28.5816,  5.9599, -0.44257, 0.003962 /) 
    b(1:7) = (/  15876.0, 1259.4,  2867.7,  -1704.56,  378.94, -36.144,  1.2071   /)
    c(1:7) = (/  11.208, -5.223,  -6.318,    4.013,   -0.822,   0.0564,  0.0      /)
    
    do iRH = 1, nRh
      if (RH(iRH) > fu_deliquescence_humidity(fu_get_material_ptr('NH4NO3')))then
        wetParticle = fu_wet_particle_features(fu_get_material_ptr('NH4NO3'), RH(iRH))
        dryVFract = 1./(wetParticle%fGrowthFactor * wetParticle%fGrowthFactor * wetParticle%fGrowthFactor)
        m = fu_dry_part_density(fu_get_material_ptr('NH4NO3')) * dryVFract / &
             & (fu_mole_mass(fu_get_material_ptr('NH4NO3')) * 1000.0 * (1. - dryVFract)) ! molality [mol/kgH2O]
        
        do iT = 1, nT
          T = NH3_HNO3_vs_NH4NO3_eq_Tmin + real(iT)
          NH3_HNO3_vs_NH4NO3_eq(iRH, iT) = 2*log(m) - 2.3523 * m**0.5 / (1.+ 0.925 * m**0.5)
          do i = 0, 6
            NH3_HNO3_vs_NH4NO3_eq(iRH, iT) = &
                          & NH3_HNO3_vs_NH4NO3_eq(iRH, iT) + &
                          & (a(i+1) - b(i+1) / T + c(i+1) * log(T)) * m**(i/2.)
          enddo
          NH3_HNO3_vs_NH4NO3_eq(iRH, iT) = &
                     & exp(NH3_HNO3_vs_NH4NO3_eq(iRH, iT)) * &
                     & 1.0e-8 / (gas_constant_uni*gas_constant_uni*T*T)   ! unit conversion from nb**2 partial pressure
        enddo
      else
        do iT = 1, nT
          T = NH3_HNO3_vs_NH4NO3_eq_Tmin + real(iT) 
          NH3_HNO3_vs_NH4NO3_eq(iRH, iT) = exp(118.87-24084/T-6.025*log(T)) * &
                     & 1.0e-8 / (gas_constant_uni*gas_constant_uni*T*T)   ! unit conversion from nb**2 partial pressure
        enddo
      endif
    enddo

!!!Old!!!
!    allocate(NH3_HNO3_vs_NH4NO3_eq(nIndT))
!    do iTmp = 1, nIndT
!      fTempr = real(iTmp) + 72.15
! Wrong way (loss of accuracy, zero for already -10C):
!      NH3_HNO3_vs_NH4NO3_eq(iTmp) = 1.12e22*(298./fTempr)**6.1 * exp(-24220./fTempr)
!
! Right way: all goes OK down to -100C
!      NH3_HNO3_vs_NH4NO3_eq(iTmp) = exp(50.7702 + 6.1*log(298./fTempr) - (24220./fTempr))
!    enddo    
!!!Old!!!

    
  end subroutine init_AerDynSimple


  !***********************************************************************

  subroutine AerDynSimple_input_needs(meteo_input_local)
    !
    ! Returns input needs for the transformation routines. 
    !
    implicit none

    ! Imported parameters
    type(Tmeteo_input), intent(out), target :: meteo_input_local

    ! Local variables
    integer :: iTmp

    ! Prepare the local mapping structure for the meteo input
    !
    meteo_input_local =  meteo_input_empty

    meteo_input_local%nQuantities = 5
    meteo_input_local%quantity(1) = temperature_flag    
    meteo_input_local%q_type(1) = meteo_dynamic_flag
    ind_tempr => meteo_input_local%idx(1)

    meteo_input_local%quantity(2) = relative_humidity_flag
    meteo_input_local%q_type(2) = meteo_dynamic_flag
    ind_rh => meteo_input_local%idx(2)
    
    meteo_input_local%quantity(3) = pressure_flag
    meteo_input_local%q_type(3) = meteo_dynamic_flag
    ind_pres => meteo_input_local%idx(3)
    
    meteo_input_local%quantity(4) = height_flag
    meteo_input_local%q_type(4) = meteo_dynamic_flag
    ind_hgt => meteo_input_local%idx(4)
    
    meteo_input_local%quantity(5) = total_cloud_cover_flag
    meteo_input_local%q_type(5) = meteo_dynamic_flag
    ind_cloud_cvr => meteo_input_local%idx(5)

  end subroutine AerDynSimple_input_needs


  !***********************************************************************

  subroutine full_spec_lst_4_AerDynSimple(rulesAerDynSimple, &
                                        & speciesEmis, speciesTransp, speciesSL, speciesAerosol, &
                                        & nSpeciesEmis, nSpeciesTransp, nSpeciesSL, nSpeciesAerosol, &
                                        & iClaimedSpecies)
    !
    ! Here we select the species to be transformed (some are just copied from gaseous to aerosol)
    ! Species that can be considered are: 
    ! sulphates: this stuff is emitted together with SO2 with unknown cations
    ! ammonium sulphate - reaction of H2SO4 and NH3, as well as of SO4_general with NH3
    ! ammonium nitrate - reaction with HNO3
    ! H2SO4 in aerosol - essentially, sulphate
    !
    implicit none

    ! Imported parameters
    type(Tchem_rules_AerDynSimple), intent(inout) :: rulesAerDynSimple
    type(silam_species), dimension(:), pointer :: speciesTransp, speciesEmis, speciesSL, speciesAerosol
    integer, intent(in) :: nSpeciesEmis
    integer, intent(inout) :: nSpeciesTransp, nSpeciesSL, nSpeciesAerosol
    integer, dimension(:), intent(inout) :: iClaimedSpecies

    ! Local variable
    integer :: nSelected, iTmp, iSO4Emis, iNH3Emis, iEmis, maxModes
    type(silam_species) :: speciesTmp
    logical :: ifSulphateProduction
    integer, dimension(:), pointer :: indices
    indices => null()
    indices => fu_work_int_array()
    ! SO4
    ! Chemistry puts SO4 and H2SO4 to short living, here make transport species
    call set_species(speciesTmp, fu_get_material_ptr('SO4'), in_gas_phase)
    if(fu_index(speciesTmp, speciesSL, nSpeciesSL) /= int_missing)then
      ifSulphateProduction = .true.
    else
      call set_species(speciesTmp, fu_get_material_ptr('H2SO4'), in_gas_phase)
      ifSulphateProduction = (fu_index(speciesTmp, speciesSL, nSpeciesSL) /= int_missing)
    endif
    if(ifSulphateProduction)then
      !
      ! If somebody already added SO4 aerosol to transport species, potentially trouble
      ! However, it can be Mid-Atmopshere dynamics, then it is OK. So, just check that
      ! the sizes are correct
      !
      call select_species(speciesTransp, nSpeciesTransp, 'SO4', aerosol_mode_missing, real_missing, &
                        & indices, nSelected)
      if(nSelected > 0)then
        call msg_warning('Aerosol SO4 found in transport species','full_spec_lst_4_AerDynSimple')
        do iTmp = 1, nSelected
          call report(speciesTransp(indices(iTmp)))
          if(.not. (fu_mode(speciesTransp(indices(iTmp))) == rulesAerDynSimple%modeCondensation .or. &
                &   fu_mode(speciesTransp(indices(iTmp))) == rulesAerDynSimple%modeDroplet .or. &
                &   fu_mode(speciesTransp(indices(iTmp))) == rulesAerDynSimple%modeCoarse))then
            call set_error('Transported SO4 has to be in one of the following modes', 'full_spec_lst_4_AerDynSimple')
            call report(rulesAerDynSimple%modeCondensation)
            call report(rulesAerDynSimple%modeDroplet)
            call report(rulesAerDynSimple%modeCoarse)
            return
          endif
        end do
      endif
    
      ! Check the emitted species. For the time being, demand them to be the same mode as the ones here
      call select_species(speciesEmis, nSpeciesEmis, 'SO4', aerosol_mode_missing, real_missing, indices, nSelected)
      if(error)return
      if(nSelected > 0)then
        do iTmp = 1, nSelected 
          if(.not. (fu_mode(speciesEmis(indices(iTmp))) == rulesAerDynSimple%modeCondensation .or. &
                &   fu_mode(speciesEmis(indices(iTmp))) == rulesAerDynSimple%modeDroplet .or. &
                &   fu_mode(speciesEmis(indices(iTmp))) == rulesAerDynSimple%modeCoarse))then
            call set_error('Emitted SO4 has to be in one of the following modes', 'full_spec_lst_4_AerDynSimple')
            call report(rulesAerDynSimple%modeCondensation)
            call report(rulesAerDynSimple%modeDroplet)
            call report(rulesAerDynSimple%modeCoarse)
            return
          else   ! claim & add to transport
            if(iClaimedSpecies(indices(iTmp)) < 0)then
              iClaimedSpecies(indices(iTmp)) = aerosol_dynamics_simple
              call msg('Simple aerosol dynamics owns emitted:' + fu_substance_name(speciesEmis(indices(iTmp))))
              call report(fu_mode(speciesEmis(indices(iTmp))))
              call set_species(speciesTmp, fu_get_material_ptr('SO4'), fu_mode(speciesEmis(indices(iTmp))))
              call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1) 
            else
              call msg('Sulphates in aerosol are claimed in emission by:', iClaimedSpecies(indices(iTmp)))
              call msg_warning('Sulphates are claimed in emission by:', 'full_spec_lst_4_AerDynSimple')
            endif
          endif
        end do
      endif
      
      ! Add SO4 species to transport (might be already added from emmission, but add_species takes care of duplicates)
      call set_species(speciesTmp, fu_get_material_ptr('SO4'), rulesAerDynSimple%modeCondensation)
      call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1) !duplicates not added
      call set_species(speciesTmp, fu_get_material_ptr('SO4'), rulesAerDynSimple%modeDroplet)
      call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1) !duplicates not added
!      call set_species(speciesTmp, fu_get_material_ptr('SO4'), rulesAerDynSimple%modeCoarse)
!      call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1) !duplicates not added
!      call set_species(speciesTmp, fu_get_material_ptr('H2SO4'), in_gas_phase)
!      call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1) !duplicates not added
      if(error)return
       rulesAerDynSimple%fParticleSrfRatio = (fu_massmean_D(rulesAerDynSimple%modeCondensation) / &
                            & fu_massmean_D(rulesAerDynSimple%modeDroplet)) ** 0.6666666666667

    endif  ! sulphate production
    
    
    ! NH3 (can be oxidised by H2SO4/SO4 and HNO3)
    call set_species(speciesTmp, fu_get_material_ptr('NH3'), in_gas_phase)
    iNH3Emis = fu_index(speciesTmp, speciesEmis, nSpeciesEmis)
    if(iNH3Emis /= int_missing)then
      ! If not claimed by anyone else, claim and add to transport
      if(iClaimedSpecies(iNH3Emis) < 0)then
        iClaimedSpecies(iNH3Emis) = aerosol_dynamics_simple
        call msg('Simple aerosol dynamics owns:' + fu_substance_name(speciesEmis(iNH3Emis)))
        call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
      endif
      if(error)return

      if(ifSulphateProduction)then ! (NH4)1.5SO4
        call set_species(speciesTmp, fu_get_material_ptr('NH415SO4'), rulesAerDynSimple%modeCondensation)
        call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
        call set_species(speciesTmp, fu_get_material_ptr('NH415SO4'), rulesAerDynSimple%modeDroplet)
        call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
!        call set_species(speciesTmp, fu_get_material_ptr('NH415SO4'), rulesAerDynSimple%modeCoarse)
!        call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
        if(error)return
      endif
      
      ! HNO3 
      call set_species(speciesTmp, fu_get_material_ptr('HNO3'), in_gas_phase)
      if(fu_index(speciesTmp, speciesTransp, nSpeciesTransp) /= int_missing)then ! NH4NO3
        !call set_species(speciesTmp, fu_get_material_ptr('NH4NO3'), rulesAerDynSimple%modeCondensation)
        !call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
        call set_species(speciesTmp, fu_get_material_ptr('NH4NO3'), rulesAerDynSimple%modeDroplet)
        call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
!        call set_species(speciesTmp, fu_get_material_ptr('NH4NO3'), rulesAerDynSimple%modeCoarse)
!        call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
        if(error)return
      endif
    endif  ! NH3 present in emission

    ! Coarse nitrates on seasalt
    if(rulesAerDynSimple%no3_StickingCoeff > 0.0)then
      call select_species(speciesEmis, nSpeciesEmis, rulesAerDynSimple%ssltName4NO3c, &
                        & aerosol_mode_missing, real_missing, indices, nSelected)
      if(nSelected > 0)then
        ! claim the seasalt and add to transport
        do iTmp = 1, nSelected
          speciesTmp = speciesEmis(indices(iTmp))
          call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
          if(iClaimedSpecies(indices(iTmp)) < 0)then
            iClaimedSpecies(indices(iTmp)) = aerosol_dynamics_simple
            call msg('Simple aerosol dynamics owns:' + fu_substance_name(speciesEmis(indices(iTmp))))
          endif
        enddo
      if(error)return
      ! NaNO3  
      if(rulesAerDynSimple%SSLTmodes4NO3c)then
        !do iTmp = 1, nSelected
        !  call set_species(speciesTmp, fu_get_material_ptr('NO3_c'), fu_mode(speciesEmis(indices(iTmp))))
        !  call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
        !enddo
        call set_error('Seasalt modes for coarse NO3 do not work yet', 'full_spec_lst_4_AerDynSimple')
      else     
        call set_species(speciesTmp, fu_get_material_ptr('NO3_c'), rulesAerDynSimple%modeCoarse)
        call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
      endif
      else
        !
        ! If sea salt has not been found in emission species, check the transport - may be,
        ! they are already there, owing to boundaries or alike
        !
        call select_species(speciesTransp, nSpeciesTransp, rulesAerDynSimple%ssltName4NO3c, &
                          & aerosol_mode_missing, real_missing, indices, nSelected)
        if(nSelected == 0)then
          call set_error('Cannot find my seasalt:' + rulesAerDynSimple%ssltName4NO3c, &
                       & 'full_spec_lst_4_AerDynSimple')
          return
        endif
      endif  ! nSelected
    endif   ! NO3 sticking coef is > 0
        
    ! Aging reactions with OH (black carbon or anything else - materials of fresh and aged species 
    ! and rate with OH taken from control file). Can be several.
    maxModes = 0
    do iTmp = 1, rulesAerDynSimple%nOHaging
      call select_species(speciesEmis, nSpeciesEmis, rulesAerDynSimple%OHaging_freshSps(iTmp), &
                          & aerosol_mode_missing, real_missing, indices, nSelected)
      if(nSelected == 0)then
        call set_error(fu_connect_strings('Cannot find in emission material requested for aging:', rulesAerDynSimple%OHaging_freshSps(iTmp)), &
                       & 'full_spec_lst_4_AerDynSimple')
        return
      else 
        maxModes = max(maxModes, nSelected)
        do iEmis = 1, nSelected
          speciesTmp = speciesEmis(indices(iEmis))
          call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
          if(iClaimedSpecies(indices(iEmis)) < 0)then
            iClaimedSpecies(indices(iEmis)) = aerosol_dynamics_simple
            call msg('Simple aerosol dynamics owns:' + fu_substance_name(speciesEmis(indices(iEmis))))
          endif
          call set_species(speciesTmp, fu_get_material_ptr(rulesAerDynSimple%OHaging_agedSps(iTmp)), fu_mode(speciesEmis(indices(iEmis))))
          call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
          if(error)return
        enddo
      endif
    enddo
    rulesAerDynSimple%maxAgingModes = maxModes
    call msg('Max nr modes for aging species:', maxModes)
    allocate(iFresh(rulesAerDynSimple%nOHaging, maxModes))
    allocate(iAged(rulesAerDynSimple%nOHaging, maxModes))    
    iFresh(:,:) = -1
    iAged(:,:) = -1   

    
    
    call free_work_array(indices)
    
  end subroutine full_spec_lst_4_AerDynSimple


  !*******************************************************************

  subroutine registerSpecies_4_AerDynSimple(rulesAerDynSimple, speciesTrans, speciesSL, &
                                                             & nSpeciesTrans, nSpeciesSL)
    !
    ! Set the species indices arrays
    !
    implicit none
    
    ! Imported parameters
    type(Tchem_rules_AerDynSimple), intent(in) :: rulesAerDynSimple
    type(silam_species), dimension(:), intent(in) :: speciesTrans, speciesSL
    integer, intent(in) :: nSpeciesTrans, nSpeciesSL

    ! Local variables
    type(silam_species) :: speciesTmp
    integer :: nSelected, iMod, iTmp, nAged
    integer, dimension(:), pointer :: indicesTmp
    
    indicesTmp => fu_work_int_array()
    !
    ! Short-lived species
    !
    call set_species(speciesTmp, fu_get_material_ptr('SO4'), in_gas_phase)  ! rulesAerDynSimple%modeCloudProcessedAerosol)  ! SO4 short-living
    iSO4w_aer_sl = fu_index(speciesTmp, speciesSL, nSpeciesSL)

    call set_species(speciesTmp, fu_get_material_ptr('H2SO4'), in_gas_phase)
    iH2SO4_gas_sl = fu_index(speciesTmp, speciesSL, nSpeciesSL)
    !
    ! Transport species
    !
    call set_species(speciesTmp, fu_get_material_ptr('SO4'), rulesAerDynSimple%modeCondensation)
    iSO4f_aer = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)

    call set_species(speciesTmp, fu_get_material_ptr('SO4'), rulesAerDynSimple%modeDroplet) !rulesAerDynSimple%modeSecondaryAerosol)
    iSO4c_aer = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)

    call set_species(speciesTmp, fu_get_material_ptr('NH415SO4'), rulesAerDynSimple%modeCondensation)
    iNH415SO4f_aer = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)

    call set_species(speciesTmp, fu_get_material_ptr('NH415SO4'), rulesAerDynSimple%modeDroplet)
    iNH415SO4c_aer = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)

    call set_species(speciesTmp, fu_get_material_ptr('NH3'), in_gas_phase)
    iNH3 = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)

    call set_species(speciesTmp, fu_get_material_ptr('HNO3'), in_gas_phase)
    iHNO3 = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)

    call set_species(speciesTmp, fu_get_material_ptr('NH4NO3'), rulesAerDynSimple%modeDroplet)
    iNH4NO3_aer = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)

    if(rulesAerDynSimple%no3_StickingCoeff > 0.0)then
      call select_species(speciesTrans, nSpeciesTrans, rulesAerDynSimple%ssltName4NO3c, &
                        & aerosol_mode_missing, real_missing, &
                        & indicesTmp, nSelected)
      if(nSelected > 0)then
        allocate(iSslt_aer(nSelected))
        allocate(sslt_D(nSelected))
        iSslt_aer(1:nSelected) = indicesTmp(1:nSelected)
        do iMod = 1, nSelected
            sslt_D(iMod) = fu_massmean_D(speciesTrans(indicesTmp(iMod)))
        enddo
      else
        call msg_warning('Lost my seasalt:' + rulesAerDynSimple%ssltName4NO3c, &
                       & 'registerSpecies_4_AerDynSimple')
        do iMod = 1, nSpeciesTrans
          call report(speciesTrans(iMod))
        end do
        call set_error('Lost my seasalt in the above list:' + rulesAerDynSimple%ssltName4NO3c, &
                     & 'registerSpecies_4_AerDynSimple')
        return
      endif
      
      call select_species(speciesTrans, nSpeciesTrans, 'NO3_c',  aerosol_mode_missing, real_missing, &
                          & indicesTmp, nSelected)
      if(nSelected > 0)then
        if(rulesAerDynSimple%SSLTmodes4NO3c .and. nSelected /= size(iSslt_aer))then
          call set_error('NO3c should exist in all seasalt modes','registerSpecies_4_AerDynSimple')
          return 
        elseif((.not. rulesAerDynSimple%SSLTmodes4NO3c) .and. nSelected /= 1)then
          call set_error('NO3c should exist only in coarse mode','registerSpecies_4_AerDynSimple')
          return
        endif 
        if(rulesAerDynSimple%SSLTmodes4NO3c)then
          allocate(iSsltNO3(nSelected))
          iSsltNO3(1:nSelected) = indicesTmp(1:nSelected)
          iNO3c_aer = int_missing
        else
          iNO3c_aer = indicesTmp(1)
          nullify(iSsltNO3)
        endif
      else
        call set_error('Lost my NO3c', 'registerSpecies_4_AerDynSimple')
      endif

    else
      nullify(iSslt_aer)
      nullify(iSsltNO3)
    endif
    ! Aging reactions with OH (black carbon or anything else - materials of fresh and aged species 
    ! and rate with OH taken from control file). Can be several.
    do iTmp = 1, rulesAerDynSimple%nOHaging
      call select_species(speciesTrans, nSpeciesTrans, rulesAerDynSimple%OHaging_freshSps(iTmp), &
                        & aerosol_mode_missing, real_missing, indicesTmp, nSelected)
      if(nSelected > 0)then
        iFresh(iTmp, 1:nSelected) = indicesTmp(1:nSelected)
      else
        call set_error('Lost aging fresh material ' + rulesAerDynSimple%OHaging_freshSps(iTmp),&
                     & 'registerSpecies_4_AerDynSimple')
        return
      endif
      do iMod = 1, nSelected
        call select_species(speciesTrans, nSpeciesTrans, rulesAerDynSimple%OHaging_agedSps(iTmp), &
                        & fu_mode(speciesTrans(iFresh(iTmp,iMod))), real_missing, indicesTmp, nAged)
        if(nAged == 1)then       
          iAged(iTmp,iMod) = indicesTmp(1)
        else
          ! This should never happen
          call set_error(fu_connect_strings('Strange nr of aged species found for', rulesAerDynSimple%OHaging_agedSps(iTmp)),&
                     & 'registerSpecies_4_AerDynSimple')
          call msg('Found:', nAged)
          return
        endif
      enddo
    enddo
    call free_work_array(indicesTmp)
    

  end subroutine registerSpecies_4_AerDynSimple

  !*******************************************************************

  subroutine transform_AerDynSimple(vMassTrn, vMassSL, TL, rulesAerDynSimple, &
                                  & metdat, seconds, zenith_cos)
    !
    ! Closes-up the loop for the gas-phase chemistry that stops at creating the aerosols
    !
    implicit none

    ! Imported parameters
    real, dimension(:), intent(inout) :: vMassTrn, vMassSL
    real, dimension(:), pointer :: TL ! Save/restore 
    type(Tchem_rules_AerDynSimple), intent(in) :: rulesAerDynSimple
    real, dimension(:), intent(in) :: metdat
    real, dimension(:), pointer :: TLams, TLamn, TLcn ! Save/restore 
    real, intent(in) :: seconds, zenith_cos
    character(len = *), parameter :: sub_name = 'transform_AerDynSimple'

    if (associated(TL)) then
       if (fu_fails(size(TL) == 6, "Wrong size TLA", sub_name)) return
       TLams => TL(1:3) !!size 3
       TLamn => TL(4:5) !! size 2
       TLcn  => TL(6:6)  !! size 1
    else
      TLams => null()
      TLamn => null()
      TLcn  => null()
    endif
 
    !! Aging not related to others. Can be done at any time
    call ocAging(vMassTrn, rulesAerDynSimple, metdat, zenith_cos, seconds)
    
    !! These should be ordered
    if (seconds > 0) then
        ! Move h2so4 from sort-lived

!        if (any( (/vMassTrn(iSO4f_aer), vMassTrn(iSO4c_aer)/) < 0)) call ooops("1")
        call so4_h2so4(vMassTrn, vMassSL, seconds)
!        if (any( (/vMassTrn(iSO4f_aer), vMassTrn(iSO4c_aer)/) < 0)) call ooops("2")


        if (iNH3 > 0) then
          ! Break ammonium nitrate, so sulphates can consume NH3
          call break_nh4no3(vMassTrn,  seconds)

          !SO4 consumes free ammonia (3 TL vlues)
          call so4_nh3(vMassTrn, TLams,rulesAerDynSimple, metdat, seconds)
  !        if (any( (/vMassTrn(iSO4f_aer), vMassTrn(iSO4c_aer)/) < 0)) call ooops("3")
          
          !restore nh4no3 to equilibrium (2 TL values)
          call nh3_hno3_to_nh4no3(vMassTrn, TLamn, metdat, seconds )

          !create coarse nitrates (1 TL value)
          call no3_to_no3c(vMassTrn, TLcn, rulesAerDynSimple, metdat, seconds)
        endif
    else
        if (iNH3 > 0) then
          call no3_to_no3c(vMassTrn, TLcn, rulesAerDynSimple, metdat, seconds)
          !create coarse nitrates (1 TL value)

          call nh3_hno3_to_nh4no3(vMassTrn, TLamn, metdat, seconds )
          !restore nh4no3 to equilibrium (2 TL values)
          
          call so4_nh3(vMassTrn, TLams, rulesAerDynSimple, metdat, seconds)
          !SO4 consumes free ammonia (3 TL vlues)

          call break_nh4no3(vMassTrn,  seconds)
          ! Break ammonium nitrate, so sulphates can consume NH3
        endif

        call so4_h2so4(vMassTrn, vMassSL, seconds)
        ! Move h2so4 from sort-lived
    endif

  end subroutine transform_AerDynSimple

  !***********************************************************************

  subroutine set_rules_AerDynSimple(nlSetup, rulesAerDynSimple)
    !
    ! Creates the set of rules, in accordance with which the computations are to
    ! be made. Input - SILAM namelist
    !
    implicit none

    ! Imported parameter
    type(Tsilam_namelist), intent(in) :: nlSetup 
    type(Tchem_rules_AerDynSimple), intent(out) :: rulesAerDynSimple

    ! local variables
    integer :: iT, nT, iStat
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: ptrItems
    character(len=fnlen) :: spContent
    character(len = *), parameter :: sub_name = 'nh3_hno3_nh4no3'


    ! 3 modes for secondary stuff (roughly from Seinfeld & Pandis 1997, page 441)
    ! Secondary particles from grown from nucleation
   !     rulesAerDynSimple%modeCondensation = fu_set_mode(fixed_diameter_flag, 0.1e-6, 0.3e-6, 0.2e-6)

   rulesAerDynSimple%modeCondensation = fu_set_mode(fixed_diameter_flag, 0.01e-6, 0.3e-6, 0.2e-6)
    ! Pre-existing / cloud processed particles
    rulesAerDynSimple%modeDroplet = fu_set_mode(fixed_diameter_flag, 0.3e-6, 1.5e-6, 0.7e-6)
    ! Sea salt
    rulesAerDynSimple%modeCoarse = fu_set_mode(fixed_diameter_flag, 1.5e-6, 5.0e-6, 3.0e-6)

!    !18.09.2016. Considerations of AOD and more specific spectrum analysis sugested the following adjustments
!    ! 3 modes for secondary stuff (roughly from Seinfeld & Pandis 1997, page 441)
!    ! Secondary particles from grown from nucleation
!    rulesAerDynSimple%modeCondensation = fu_set_mode(fixed_diameter_flag, 0.01e-6, 0.2e-6, 0.1e-6)
!    ! Pre-existing / cloud processed particles
!    rulesAerDynSimple%modeDroplet = fu_set_mode(fixed_diameter_flag, 0.2e-6, 1.e-6, 0.5e-6)
!    ! Sea salt
!    rulesAerDynSimple%modeCoarse = fu_set_mode(fixed_diameter_flag, 1.5e-6, 5.0e-6, 3.0e-6)

    ! Actually everything should go to all modes, but for now: 
    ! gas phase produced SO4 (H2SO4) to condensation
    ! heterogeneous SO4 to droplet
    ! (NH4)1.5SO4 follows SO4
    ! NH4NO3 to droplet 
    ! NaNO3 to coarse / all seasalt modes if requested (to be done)
    ! NaSO4 & NH4Cl to be done:
    ! Pio, C.A. and Harrison, R.M. (1987). The Equilibrium of Ammonium Chloride Aerosol with Gaseous 
    ! Hydrochloric Acid and Ammonia under Tropospheric Conditions. Atmos. Environ. 21: 1243–1246
    ! supposedly claims 1.5-2 times higher eq constant for HCl than HNO3 for same amount of NH3.
    

    ! Coarse NO3 & SO4
    nullify(ptrItems)
    call get_items(nlSetup, 'make_coarse_no3', ptrItems, nT)
    if (nT == 1)then
      spContent = fu_content(ptrItems(1))
      read(unit=spContent, iostat=iStat, fmt=*) rulesAerDynSimple%ssltName4NO3c, rulesAerDynSimple%no3_StickingCoeff
      if(iStat /= 0)then
        call set_error('Invalid make_coarse_no3 line in transformation rules','set_rules_AerDynSimple')
        call msg('Required format:   make_coarse_no3 = ssltName alpha')
        return
      else
        call msg(fu_connect_strings('Making coarse NO3 on seasalt: ', rulesAerDynSimple%ssltName4NO3c), &
                                   & rulesAerDynSimple%no3_StickingCoeff)
      endif
      if(index(spContent,'SSLTmodes') > 0)then
        rulesAerDynSimple%SSLTmodes4NO3c = .true.
        call msg('Using seasalt modes for coarse NO3')
      else
        rulesAerDynSimple%SSLTmodes4NO3c = .false.
      endif
    elseif(nT == 0)then
      call msg('No coarse NO3 production on seasalt')
      rulesAerDynSimple%ssltName4NO3c = ''
      rulesAerDynSimple%no3_StickingCoeff = -1.0
    else
      ! In future might allow several, to form also on dust ..
      call set_error('Only one make_coarse_no3 line allowed in transformation rules','set_rules_AerDynSimple')
    endif
    if(error)return
    
    call get_items(nlSetup, 'make_coarse_so4', ptrItems, nT)
    if (nT == 1)then
      spContent = fu_content(ptrItems(1))
      read(unit=spContent, iostat=iStat, fmt=*) rulesAerDynSimple%ssltName4SO4c, rulesAerDynSimple%so4_StickingCoeff
      if(iStat /= 0)then
        call set_error('Invalid make_coarse_so4 line in transformation rules','set_rules_AerDynSimple')
        call msg('Required format:   make_coarse_so4 = ssltName alpha')
        return
      else
        call msg(fu_connect_strings('Making coarse SO4 on seasalt: ', rulesAerDynSimple%ssltName4SO4c), &
                                   & rulesAerDynSimple%so4_StickingCoeff)
      endif
      if(index(spContent,'SSLTmodes') > 0)then
        rulesAerDynSimple%SSLTmodes4SO4c = .true.
        call msg('Using seasalt modes for coarse SO4')
      else
        rulesAerDynSimple%SSLTmodes4SO4c = .false.
      endif
    elseif(nT == 0)then
      call msg('No coarse SO4 production on seasalt')
      rulesAerDynSimple%ssltName4So4c = ''
      rulesAerDynSimple%so4_StickingCoeff = -1.0
    else
      ! In future might allow several ...
      call set_error('Only one make_coarse_so4 line allowed in transformation rules','set_rules_AerDynSimple')
    endif
    if(error)return
    
    ! Aging reactions with OH (black carbon or anything else - materials of fresh and aged species 
    ! and rate with OH taken from control file). Can be several.
    nullify(iFresh) ! Allocated when emitted modes are known
    nullify(iAged)  ! Allocated when emitted modes are known
    rulesAerDynSimple%maxAgingModes = -1
    call get_items(nlSetup, 'OH_aging', ptrItems, nT)
    rulesAerDynSimple%nOHaging = nT
    call msg('Number of OH aging reactions:',nT)
    if(nT == 0)then
      nullify(rulesAerDynSimple%OHaging_freshSps)
      nullify(rulesAerDynSimple%OHaging_agedSps)
      nullify(rulesAerDynSimple%OHaging_rates)
    else
      allocate(rulesAerDynSimple%OHaging_freshSps(nT))
      allocate(rulesAerDynSimple%OHaging_agedSps(nT))
      allocate(rulesAerDynSimple%OHaging_rates(nT))
      do iT = 1, nT
        spContent = fu_content(ptrItems(iT))
        read(unit=spContent, iostat=iStat, fmt=*) rulesAerDynSimple%OHaging_freshSps(iT), &
               & rulesAerDynSimple%OHaging_agedSps(iT), rulesAerDynSimple%OHaging_rates(iT)
        if (iStat /= 0) then
          call msg("parsing failed '"//trim(spContent)//"'")
          call set_error("Failed to parse OH_aging", sub_name)
        endif
        call msg(rulesAerDynSimple%OHaging_freshSps(iT) + ' -> ' + rulesAerDynSimple%OHaging_agedSps(iT), rulesAerDynSimple%OHaging_rates(iT))
      enddo
    endif

    

    rulesAerDynSimple%defined = silja_true
  end subroutine set_rules_AerDynSimple

  
  !!!!!************************************************************************************
  !!!!!
  !!!!! Effects on clouds
  !!!!!
  !!!!!************************************************************************************
  !!!!!************************************************************************************
  !!!!
  !!!!subroutine make_cld_droplet_nbr_cnc(mapCnc, met_buf, disp_buf)
  !!!!  !
  !!!!  ! Computes the number concentration of cloud droplets.
  !!!!  ! Sometimes, similar relations are used for condensation / ice nuclei
  !!!!  ! See: Menon et al, J.Atmosph.Sci., 2002, v.59, 692-714
  !!!!  !
  !!!!  implicit none
  !!!!  
  !!!!  ! Imported parameters
  !!!!  type(Tmass_map), pointer :: mapCnc
  !!!!  type(Tfield_buffer), pointer :: met_buf, disp_buf
  !!!!  
  !!!!  ! Local variables
  !!!!  integer :: ix, iy, iz
  !!!!  real, parameter :: def_CDC = 10 ** (2.41 - 6)  ! To make it per m3
  !!!!
  !!!!  do iy = 1, mapCnc%ny
  !!!!    do ix = 1, mapCnc%nx
  !!!!      i1d_disp = (iy-1)*mapCnc%nx + ix
  !!!!      do iz = 1, mapCnc%n3d
  !!!!        dz_past => disp_buf%p4d(ind_dz)%past%p2d(iz)%ptr
  !!!!        dz_future => disp_buf%p4d(ind_dz)%future%p2d(iz)%ptr
  !!!!        CDC_now => disp_buf%p4d(ind_CDC)%past%p2d(iz)%ptr
  !!!!        dz = met_buf%weight_past*dz_past(i1d_disp) + (1.0-met_buf%weight_past)*dz_future(i1d_disp)
  !!!!        cell_volume = dz * cell_size_x(i1d_disp) * cell_size_y(i1d_disp)
  !!!!        !
  !!!!        ! Get the concentrations of the components: sulphates, organic matter and sea salt.
  !!!!        ! Depending on the run setup, these can be different species. Careful
  !!!!        !
  !!!!        cncSO4_ugm3 = 0.0
  !!!!        do iSp = 1, max_index_value_CDC
  !!!!          if(indSO4(iSp) < 0) cncSO4_ugm3 = cncSO4_ugm3 + mole2ug_SO4 / cell_volume * &
  !!!!                                            & sum(mapCnc%arM(indSO4(iSp),1:mapCnc%nSrc,iz,ix,iy))
  !!!!          if(indOM(iSp) < 0) cncOM_ugm3 = cncOM_ugm3 + kg2ug_OM / cell_volume * &
  !!!!                                            & sum(mapCnc%arM(indOM(iSp),1:mapCnc%nSrc,iz,ix,iy))
  !!!!          if(indSSLT(iSp) < 0) cncSSLT_ugm3 = cncSSLT_ugm3 + kg2ug_SSLT  / cell_volume * &
  !!!!                                            & sum(mapCnc%arM(indSSLT(iSp),1:mapCnc%nSrc,iz,ix,iy))
  !!!!        end do
  !!!!        !
  !!!!        ! The main formula here is modified from Menon et al.
  !!!!        ! In the paper, zeroying or very small value of any of the components
  !!!!        ! destroys the whole thing by putting cloud droplets cnc to zero.
  !!!!        ! To avoid that, the formula is changed by adding 1 ug/m3 to each of the species
  !!!!        ! concentrations
  !!!!        ! Initial formula: N = 10^ (2.41 + 0.5 * lg(SO4) + 0.13 * lg(OM) + 0.05 * lg(sslt))
  !!!!        ! Modified: N = 10^ (2.41 + 0.5 * lg(SO4 + 1) + 0.13 * lg(OM + 1) + 0.05 * lg(sslt + 1))
  !!!!        ! All concentrations are in ug/m3 of the species, N is in #/cm3.
  !!!!        !
  !!!!        ! Below, CDC is in SI, i.e. #/m3.
  !!!!        !
  !!!!        CDC_now(i1d_disp) = def_CDC * ((cncSO4_ugm3 + 1.0) ** 0.5) * &
  !!!!                                    & ((cncOM_ugm3 + 1.0) ** 0.13) * &
  !!!!                                    & ((cnc_SSLT_ugm3 + 1.0) ** 0.05)
  !!!!      end do   ! iz
  !!!!    end do   ! ix
  !!!!  end do   ! iy
  !!!!  
  !!!!  
  !!!!  call set_error('Does not work yet','make_cld_droplet_nbr_cnc')
  !!!!
  !!!!end subroutine make_cld_droplet_nbr_cnc
  !!!!
  
  
  !************************************************************************************

  integer function fu_tla_size_ADS(rules) result(n)
    ! Collect transformations' requests for tangent linearization. If tangent linear
    ! is not allowed, corresponding subroutine sets error. 
    implicit none
    type(Tchem_rules_AerDynSimple), intent(in) :: rules
    
    n = 6
    
  end function fu_tla_size_ADS

   !*******************************************************
  
  subroutine ocAging(vMassTrn, rulesAerDynSimple, metdat, zenith_cos, seconds)
    ! Aging reactions with OH (black carbon or anything else - materials of fresh and aged species 
    ! and rate with OH taken from control file). Can be several.
    !
      implicit none
      ! Imported parameters
      real, dimension(:), intent(inout) :: vMassTrn
      type(Tchem_rules_AerDynSimple), intent(in) :: rulesAerDynSimple
      real, dimension(:), intent(in) :: metdat
      real, intent(in) :: zenith_cos, seconds

      ! Local variables
      integer :: iMod, iTmp
      real :: differ, sun, cOH_forced, aging
      character(len = *), parameter :: sub_name = 'ocAging'

      if(rulesAerDynSimple%nOHaging > 0)then
        sun = max(zenith_cos, 0.0) * (1. - metdat(ind_cloud_cvr)*0.5)
        cOH_forced = fu_OH_cnc_forced(sun, metdat(ind_hgt))
        do iTmp = 1, rulesAerDynSimple%nOHaging
          aging = exp(-abs(seconds) * rulesAerDynSimple%OHaging_rates(iTmp) * cOH_forced)
          do iMod = 1, rulesAerDynSimple%maxAgingModes
            if(iFresh(iTmp, iMod) < 1)cycle
            if (seconds>0) then 
              differ = vMassTrn(iFresh(iTmp, iMod)) * aging                        ! / 1 - aging     0 \  !
              vMassTrn(iFresh(iTmp, iMod)) = vMassTrn(iFresh(iTmp, iMod)) - differ !|                   | !
              vMassTrn(iAged(iTmp, iMod)) = vMassTrn(iAged(iTmp, iMod)) + differ   ! \ aging         1 /  !
            else
              differ = vMassTrn(iFresh(iTmp, iMod)) * aging
              vMassTrn(iFresh(iTmp, iMod)) = vMassTrn(iFresh(iTmp, iMod)) * (1. - aging) + vMassTrn(iAged(iTmp, iMod)) * aging
              ! vMassTrn(iAged(iTmp, iMod)) stays as is
              ! / 1 - aging     aging \                     !
              !|                       |                    !
              ! \  0                1 /                     !
            endif
          enddo
        enddo 
      endif
  end subroutine ocAging
  
  !************************************************************************************

  subroutine so4_h2so4(vMassTrn, vMassSL, seconds)
    ! Two short-living species are sent in: H2SO4 and SO4 - from gas-phase and heterogeneous
    ! oxidations, respectively.
    ! First move sulphates created this timestep to sulphuric acid aerosol
    !
    implicit none

    ! Imported parameters
    real, dimension(:), intent(inout) :: vMassTrn, vMassSL
    real, intent(in) :: seconds
    character(len = *), parameter :: sub_name = 'so4_h2so4'


    if (seconds > 0) then
      if(iH2SO4_gas_sl /= int_missing)then                                 ! fine mode
        vMassTrn(iSO4f_aer) = vMassTrn(iSO4f_aer) + vMassSL(iH2SO4_gas_sl)
        vMassSL(iH2SO4_gas_sl) = 0.0 
      endif
      if(iSO4w_aer_sl /= int_missing)then                                   ! coarse mode
        vMassTrn(iSO4c_aer) = vMassTrn(iSO4c_aer) + vMassSL(iSO4w_aer_sl)
        vMassSL(iSO4w_aer_sl) = 0.0 
      endif
    else
      ! Adjoint
      if(iH2SO4_gas_sl /= int_missing)then                                 ! fine mode
        vMassSL(iH2SO4_gas_sl) = vMassTrn(iSO4f_aer)
      endif
      if(iSO4w_aer_sl /= int_missing)then                                   ! coarse mode
        vMassSL(iSO4w_aer_sl) = vMassTrn(iSO4c_aer)
      endif
      return
    end if
  end subroutine so4_h2so4

  !************************************************************************************

  subroutine break_nh4no3(vMassTrn,  seconds)
    ! Break nh4no3: it will simplify uptake of NH3 by sulphates
     !  needed nh4no3 will be restored later by nh3_hno3_nh4no3
    !
    implicit none
    ! Imported parameters
    real, dimension(:), intent(inout) :: vMassTrn
    real, intent(in) :: seconds

    character(len = *), parameter :: sub_name = 'break_nh4no3'
 
    if (seconds >0) then
        vMassTrn(iNH3) = vMassTrn(iNH3) + vMassTrn(iNH4NO3_aer)      ! /  1  0  1 \ !  
        vMassTrn(iHNO3) = vMassTrn(iHNO3) + vMassTrn(iNH4NO3_aer)    !|   0  1  1  |!
        vMassTrn(iNH4NO3_aer) = 0                                    ! \  0  0  0 / !
     else 
        !vMassTrn(iNH3)  asis
        !vMassTrn(iHNO3) asis
        vMassTrn(iNH4NO3_aer) = vMassTrn(iNH3) + vMassTrn(iHNO3) 
    endif
  end subroutine break_nh4no3

  !************************************************************************************

  subroutine so4_nh3(vMassTrn, TL, rulesAerDynSimple, metdat, seconds)
    ! Uptake of ammonia by sulphates. needs 3 values for TL
    ! This subroutine not equivalent to v5_8:
    ! 1) Uptake can't be mass-wise and surface-wise at the same time
    !    We might end up destroying more SO4 than we have. 
    !    It was creating negative masses, see workarounds/hacks in v5_8
    !    Here we do everything mass-wise, since we try to keep mass budget
    ! 2) in v5_8 sulphates took ammonia only from fine-mode NH4NO3,
    !    taking it from the coarse-mode was done implicitly via 
    !    NH4NO3 saturation over several timesteps
    !    Here sulphates take any ammonia explicitly
    
    ! Uptake of amonuim from amonuim sulphate is also handled here 
    ! due to break_nh4no3 called before
    implicit none
    ! Imported parameters
    real, dimension(:), intent(inout) :: vMassTrn 
    real, intent(in) :: seconds
    type(Tchem_rules_AerDynSimple), intent(in) :: rulesAerDynSimple
    real, dimension(:), intent(in) :: metdat
    real, dimension(:), pointer :: TL ! Save/restore 
    character(len = *), parameter :: sub_name = 'so4_nh3'

    ! Local variables
    real SO4f, SO4c, NH3, SO4tot, excNH3, fracc, fracf, deltaSc
 
    if (seconds >0) then
      SO4c = vMassTrn(iSO4c_aer)
      SO4f = vMassTrn(iSO4f_aer)
      NH3  =  vMassTrn(iNH3) 
      if (associated(TL)) TL(1:3) = (/SO4c, SO4f, NH3/) !!Store TL
    else
      if (fu_fails(associated(TL), "so4_nh3 adjoint without TL point", sub_name )) return
      SO4c = TL(1)
      SO4f = TL(2)
      NH3  = TL(3)
    endif

    ! Valent sulphates take all what they need out of NH3. Ammonium sulphate receives the product

    SO4tot  = SO4c + SO4f
    excNH3  = NH3 - SO4tot * 1.5

    if(excNH3 <= 0) then ! 1 mole of NH4_1.5_SO4 takes 1.5 moles of NH3
        !MAS Too little NH3, sulpathes take it all 
        fracc = SO4c / SO4tot ! 
        fracf = SO4f / SO4tot ! !!!NOT (1 - fracc)!!

      if (seconds > 0) then
        vMassTrn(iSO4c_aer) =  - excNH3 * fracc / 1.5
        vMassTrn(iSO4f_aer) =  - excNH3 * fracf / 1.5

        vMassTrn(iNH415SO4c_aer) = vMassTrn(iNH415SO4c_aer) + vMassTrn(iNH3) * fracc / 1.5
        vMassTrn(iNH415SO4f_aer) = vMassTrn(iNH415SO4f_aer) + vMassTrn(iNH3) * fracf / 1.5
        vMassTrn(iNH3) = 0.
      else
        deltaSc =  - vMassTrn(iSO4c_aer) + vMassTrn(iSO4f_aer) + vMassTrn(iNH415SO4c_aer) - vMassTrn(iNH415SO4f_aer)

        vMassTrn(iSO4c_aer) = vMassTrn(iSO4c_aer) + NH3 * fracf * deltaSc / (1.5 * SO4tot)
        vMassTrn(iSO4f_aer) = vMassTrn(iSO4f_aer) - NH3 * fracc * deltaSc / (1.5 * SO4tot)
        !!vMassTrn(iNH415SO4c_aer) AS IS
        !!vMassTrn(iNH415SO4f_aer) AS IS
        vMassTrn(iNH3) = 1./1.5 * (fracc*(-vMassTrn(iSO4c_aer) + vMassTrn(iNH415SO4c_aer)) + &
                            &      fracf*(-vMassTrn(iSO4f_aer) + vMassTrn(iNH415SO4f_aer)) )
      endif
    else 
        ! abundance of ammonia, engage all SO4, free ammonia left
      if (seconds > 0) then
        vMassTrn(iSO4c_aer) = 0.0
        vMassTrn(iSO4f_aer) = 0.0
        vMassTrn(iNH415SO4c_aer) = vMassTrn(iNH415SO4c_aer) + SO4c
        vMassTrn(iNH415SO4f_aer) = vMassTrn(iNH415SO4f_aer) + SO4f
        vMassTrn(iNH3) = excNH3
      else
        vMassTrn(iSO4c_aer) = vMassTrn(iNH415SO4c_aer)  - 1.5 * vMassTrn(iNH3)
        vMassTrn(iSO4f_aer) = vMassTrn(iNH415SO4f_aer)  - 1.5 * vMassTrn(iNH3)
        ! vMassTrn(iNH415SO4f_aer) ASIS
        ! vMassTrn(iNH415SO4c_aer) ASIS
        ! vMassTrn(iNH3) 
      endif
    endif  ! NH3-SO4 ratio
  end subroutine so4_nh3

  !************************************************************************************

  real function fu_cequil(T,rh) result(CEQUIL)
    implicit none
    real, intent(in) :: T, rh
    integer indT, indRH

    ! NH3-HNO3-NH4NO3 equilibrium
    indT = nint(T-NH3_HNO3_vs_NH4NO3_eq_Tmin)
    if (indT < 1) then
      call msg("Temperature below NH4NO3eq lookup table range", T)
      indT = 1
    elseif (indT > NH3_HNO3_vs_NH4NO3_eq_nT) then
      call msg("Temperature above NH4NO3eq lookup table range", T)
      indT = NH3_HNO3_vs_NH4NO3_eq_nT
    endif

    if (rh >= 1.0) then
        indRH = 22
    elseif (rh >= 0.95) then
        indRH = int((rh - 0.775)*100) 
    elseif (rh >= 0.85) then
        indRH = int((rh - 0.60)*50) 
    elseif (rh >= 0.70) then
        indRH = int((rh - 0.475)*33.333) 
    elseif (rh >= 0.50) then
        indRH = int((rh - 0.40)*25)
    else
        indRH = 1
    endif
    
    !Species here in concentrations; CEQUIL given for nb**2 partial pressures
    CEQUIL = NH3_HNO3_vs_NH4NO3_eq(indRH, indT)
      
  end function fu_cequil

  !************************************************************************************

  subroutine nh3_hno3_to_nh4no3(vMassTrn, TL, metdat, seconds )
    !
    ! NH3-HNO3-NH4NO3 equilibrium.
    ! WARNIG! nh4no3  _Must_ have been set to zero before by brak_nh4o3 in forward run
    !
    implicit none

    ! Imported parameters
    real, dimension(:), intent(inout) :: vMassTrn
    real, dimension(:), intent(in) :: metdat
    real, dimension(:), pointer :: TL
    real, intent(in) :: seconds 
    character(len = *), parameter :: sub_name = 'nh3_hno3_nh4no3'

    ! Local variables
    real :: NH3, HNO3, NH4NO3, CEQUIL, CEQUIL2, dNH3, dHNO3, dNH4NO3, fTmp, delta, eps
    
    !Product of gas concentrations that equilibrates the aerosol fraction
    !Species here in concentrations; CEQUIL is geometric mean of concentrations
    CEQUIL2 = fu_cequil(metdat(ind_tempr), metdat(ind_rh))
    
    if (seconds > 0.) then
      NH3 = vMassTrn(iNH3) 
      HNO3 = vMassTrn(iHNO3)
      NH4NO3 = vMassTrn(iNH4NO3_aer)
      if (associated(TL)) TL(1:2) = (/NH3, HNO3/) !!Store TL
    else
      NH3  = TL(1)
      HNO3 = TL(2)
      NH4NO3 = 0.
      dNH3 = vMassTrn(iNH3)
      dHNO3 = vMassTrn(iHNO3)
      dNH4NO3 = vMassTrn(iNH4NO3_aer)
    endif
    !  vMassTrn(iNH4NO3_aer) _Must_ have been set to zero before by brak_nh4o3
    if (.not. ( NH3 >= 0. .and. HNO3 >=0 .and. NH4NO3 == 0.)) then
        call msg("NH3, HNO3, NH4NO3", (/NH3, HNO3, NH4NO3/))
        call set_error('Negative HNO3 before HNO3-NH3-NO3 eq', sub_name)
        return
    endif

    eps = NH3 * HNO3 - CEQUIL2
    if  (eps <= 0.  ) return !! No aerosol, do nothing. NH4NO3 should be already broken
                          !! Also handles both-zero case


    CEQUIL = sqrt(CEQUIL2)
    delta = 0.5 * (NH3 - HNO3) / CEQUIL

    if (delta > 1e3) then  !! sqrt(delta*delta) + 1 \simeq delta
      !NH3 >> iNO3
      if (seconds > 0) then 
        vMassTrn(iNH4NO3_aer) = HNO3       ! / 0    1    0 \     /  0  \  !
        vMassTrn(iHNO3) = 0.               !|  0    0    0  | * | HNO3  | !
        vMassTrn(iNH3) = -HNO3  + NH3      ! \ 0   -1    1 /     \ NH3 /  !
      else
        vMassTrn(iNH4NO3_aer) = 0.         ! / 0    0    0 \      / dNH4NO3 \   !
        vMassTrn(iHNO3)  =  dNH4NO3 - dNH3 !|  1    0   -1  |  * |  dHNO3    |  !
        !vMassTrn(iNH3) =  dNH3 !! as is   ! \ 0    0    1 /      \ dNH3    /   !
      endif
    elseif (delta < -1e3) then !!! sqrt(delta*delta) + 1 \simeq  (- delta)
      !NH3 << iNO3
      if (seconds > 0) then 
        vMassTrn(iNH4NO3_aer) = NH3        ! / 0    0    1 \     /  0  \  !
        vMassTrn(iHNO3) = HNO3  - NH3      !|  0    1   -1  | * | HNO3  | !
        vMassTrn(iNH3) =  0                ! \ 0    0    0 /     \ NH3 /  !
      else
        vMassTrn(iNH4NO3_aer) = 0          ! / 0    0    0 \      / dNH4NO3 \   !
        !!vMassTrn(iHNO3) = dNH3 !as is    !|  0    1    0  |  * |  dHNO3    |  !
        vMassTrn(iNH3) = dNH4NO3 - dHNO3   ! \ 1   -1    0 /      \ dNH3    /   !
      endif
    else
      eps = eps /  ((NH3 + HNO3)**2) !!eps is strictly positive dimensionless now

      if (eps < 1e-4 ) then  !! Tiny fraction goes to aerosol
        if (seconds > 0) then
          fTmp =  (NH3 + HNO3) * eps    !! Final cnc of NH4NO3 
          vMassTrn(iNH4NO3_aer) = fTmp   ! / 0    bcTmp   acTmp \     /  0  \  !
          vMassTrn(iHNO3) = HNO3 - fTmp  !|  0  1-bcTmp  -acTmp  | * | HNO3  | !
          vMassTrn(iNH3) = NH3 - fTmp    ! \ 0   -bcTmp 1-acTmp /     \ NH3 /  !
        else
          !!  acTmp =  (HNO3**2 + CEQUIL2) / ((NH3 + HNO3)**2)  !! d NH4NO3 / d NH3
          !!  bcTmp =  (NH3**2  + CEQUIL2) / ((NH3 + HNO3)**2)  !! d NH4NO3 / d HNO3
           ! / 0        0         0   \     / dNH4NO3 \  !
           !|  bcTmp  1-bcTmp  -bcTmp  | * |    dHNO3  | !
           ! \ acTmp   -acTmp 1-acTmp /     \    dNH3 /  !

          fTmp = dNH4NO3 - dHNO3 - dNH3
          vMassTrn(iNH4NO3_aer) = 0
          vMassTrn(iHNO3) = dHNO3 + (NH3**2  + CEQUIL2) / ((NH3 + HNO3)**2) * fTmp
          vMassTrn(iNH3)  = dNH3  + (HNO3**2 + CEQUIL2) / ((NH3 + HNO3)**2) * fTmp
        endif
      else !! substantial fracion to aerosol
        if (seconds > 0) then
          fTmp =  CEQUIL * sqrt(1. + delta*delta) !!Dimension of cnc
          vMassTrn(iNH4NO3_aer) = 0.5 * (NH3 + HNO3) - fTmp
          vMassTrn(iHNO3) = 0.5 * (HNO3 -  NH3)  + fTmp
          vMassTrn(iNH3)  = 0.5 * (NH3 -  HNO3)  + fTmp
        else
          fTmp =  0.5 * delta / sqrt(1. + delta*delta) ! dimensionless, range (-0.5: +0.5)
          vMassTrn(iNH4NO3_aer) = 0.
          vMassTrn(iHNO3) = (0.5 + fTmp) * dNH4NO3 + (0.5 - fTmp) * dHNO3 - (0.5 + fTmp) * dNH3 
          vMassTrn(iNH3)  = (0.5 - fTmp) * dNH4NO3 - (0.5 - fTmp) * dHNO3 + (0.5 + fTmp) * dNH3
        endif
      endif 
    endif

    if (seconds > 0) then
      if (.not. all((/vMassTrn(iNH4NO3_aer),vMassTrn(iHNO3),vMassTrn(iNH3)/)>= 0.))then
          call msg("nh3_hno3_to_nh4no3 failed output (HNO3, NH3, NH4NO3):", &
                          & (/vMassTrn(iHNO3),vMassTrn(iNH3),vMassTrn(iNH4NO3_aer)/))
          call msg("nh3_hno3_to_nh4no3 input (HNO3, NH3, CEQUIL2):", (/HNO3, NH3, CEQUIL2/))
          call set_error("Negative mass on the output", sub_name)
      endif
    endif
  end subroutine nh3_hno3_to_nh4no3

  !*******************************************************************

  subroutine no3_to_no3c(vMassTrn, TL, rulesAerDynSimple, metdat, seconds)
    !
    ! And finaly try to make some more NO3 in seasalt particles (NaNO3) 
    ! MP: Rather retarded parmeterization assuming diffusion of gas to particle srf as main limit, 
    ! MP: amount of sodium as sharp limit and a tunable sticking coefficient to be interpreted however we wish.
    !
    ! RK: Note that the thing does not deplete seasalt at all, so eventually _all_ HNO3 will be converted
    ! RK: at _any_ presence of sea salt
    !
    implicit none

    ! Imported parameters
    real, dimension(:), intent(inout) :: vMassTrn
    real, dimension(:), pointer :: TL
    type(Tchem_rules_AerDynSimple), intent(in) :: rulesAerDynSimple
    real, dimension(:), intent(in) :: metdat
    real, intent(in) :: seconds

    ! Local variables
    integer :: iMod, iTmp
    real ::  D_HNO3, mfp, sink, Kn, na, fracconv
    type(TwetParticle) :: wetParticle
    
    if (seconds > 0.) then
      if(rulesAerDynSimple%no3_StickingCoeff > 0.0)then
        if(vMassTrn(iHNO3) > 0.0)then
          na = 0.0
          do iMod = 1, size(iSslt_aer)
            na = na + vMassTrn(iSslt_aer(iMod)) * 0.4 / 23.0e-3 ! moles of sodium (0.3 ->0.4 includes other cations)
          enddo
          sink = 0.0
          if(na > 0)then
              wetParticle = fu_wet_particle_features(fu_get_material_ptr('sslt'), metdat(ind_rh))

              D_HNO3 = 13.13e-6 * (metdat(ind_tempr)/273.15)**1.75 * std_pressure_sl/metdat(ind_pres)  ! diffusion coefficient [m2/s]
              mfp = 2. * D_HNO3 * sqrt(pi * fu_mole_mass(fu_get_material_ptr('HNO3')) / &
                                            & (8. * gas_constant_uni * metdat(ind_tempr)))! mean free path [m]
              do iMod = 1, size(iSslt_aer)
                  Kn = 2. * mfp / wetParticle%fGrowthFactor / sslt_D(iMod)
                  sink = sink + 12. * wetParticle%fGrowthFactor * D_HNO3 *  vMassTrn(iSslt_aer(iMod)) *  &
                      & (1. + Kn) / (1. + (2. * Kn * (1. + Kn)) / rulesAerDynSimple%no3_StickingCoeff) / &  
                      & (fu_dry_part_density(fu_get_material_ptr('sslt')) * &
                      & sslt_D(iMod) * sslt_D(iMod))
              enddo
          endif
          fracconv = exp(- seconds * sink) !! Not-converted fraction
        else
          fracconv = 1.
        endif
      else
        fracconv = 1.
      endif
      if (associated(TL)) TL(1) = fracconv !!Store TL

      vMassTrn(iNO3c_aer) = vMassTrn(iNO3c_aer) + (1 - fracconv) * vMassTrn(iHNO3)  !   1          (1 - fracconv)
      vMassTrn(iHNO3)     = fracconv * vMassTrn(iHNO3)                              !   0          fracconv          
    else
      fracconv = TL(1)
      !vMassTrn(iNO3c_aer) = vMassTrn(iNO3c_aer)   ! ASIS
      vMassTrn(iHNO3)  = (1. - fracconv) * vMassTrn(iNO3c_aer) + fracconv * vMassTrn(iHNO3)
    endif
  end subroutine no3_to_no3c

END MODULE aer_dyn_simple
