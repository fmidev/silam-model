MODULE chem_dep_acid_basic
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
  ! Author Mikhail Sofiev mikhail.sofiev@fmi.fi
  !
  ! Language: ANSI FORTRAN-90 (or close to it)
  !
  use cocktail_basic
  use hydroxyl

  implicit none
  private

  public init_chemicals_acid_basic
  public inventory_acid_basic
  public registerSpeciesAcidBasic
  public acidBasic_input_needs
  public transform_acid_basic

  ! PUBLIC routines of Tchem_rules_DMAT_S 
  !
  public set_chem_rules_acidBasic  ! returns the pointer to the rules
  public set_missing
  public fu_lifetime
  
  public fu_acidBasic_name
  public fu_OH_cnc_forced
  public fu_if_specific_deposition
  !
  ! Private routines of the sulphur DMAT cocktail and chemical rules
  !
!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***
!
! ATTENTION. All what is between the above !!!*** lines is to be swapped for the actual
!            SILAM v.4.1.2.
!            They reflect the only important thing: for box model meteo data have to be given
!            from outside, not requested from meteobuffer.
!
  !private transform_acid_basic
!  public transf_acidBasic_single_val
!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***
  
  private set_miss_chem_rules_acidBasic
  private fu_lifetime_acidBasic
  private fu_if_specific_dep_acid_basic

  private int_ptr


  interface set_missing
    module procedure set_miss_chem_rules_acidBasic
  end interface

  interface fu_lifetime
    module procedure fu_lifetime_acidBasic
  end interface

  interface fu_if_specific_deposition
    module procedure fu_if_specific_dep_acid_basic
  end interface

  public fu_if_tla_required
  private fu_if_tla_required_acid_basic
  interface fu_if_tla_required
     module procedure fu_if_tla_required_acid_basic
  end interface

  !--------------------------------------------------------------------
  !
  ! Tchem_rules_acidBasic type definition
  ! Determines the rules of operating with the NOx-NHx-O3 following the basic acid chemistry
  !
  type Tchem_rules_acidBasic
    private
    real :: vdSO2, vdSO4, vdNO2, vdPAN, vdNO3_h, vdNH3, vdHNO3
    !integer :: emisBioVOCMethod
    logical :: ifNitrogen
    !
    ! Temperature-dependent coefficients
    !
    real, dimension(:), pointer :: NO_2_NO2_with_O3, NO2_2_NO3_dark, &
                                 & PAN_radical_prod, NH3_HNO3_vs_NH4NO3_eq, PAN_2_NO2_sun, &
                                 & NO_2_NO2_direct, NO_NO3_2_NO2, NO3_NO2_N2O5_eq, &
                                 & HO2_NO_2_NO2,  H2O2_OH_2_HO2, HO2_O3_2_OH, OH_O3_2_HO2, &
                                 & CH4_OH_2_CH3O2, CH3O2_NO_2_HCHO_NO2, CH3O2_HO2_2_CH3OOH, &
                                 & CH3OOH_OH_2_CH3O2, HCHO_OH_2_CO, HO2_HO2_H2O_2_H2O2, &
                                 & HO2_PANrad_2_CH3O2, HONO_OH_2_NO2, CH4_cnc_global, CH4_vert_prof, &
                                 & CH3O2_CH3O2_2_HCHO, CH3O2_CH3O2_2_CH3O
    !
    ! Temperature and pressure / altitude
    !
    real, dimension(:,:), pointer :: SO2_OH_2_SO3, PANrad_2_PAN_eq, NO2_OH_2_HNO3, CO_OH_2_CO2, &
                                   & HO2_HO2_2_H2O2_basic, O1d_N2, O1d_O2, &
                                   & HO2NO2_2_HO2_NO2, HO2_NO2_2_HO2NO2, NO_OH_2_HONO, &
                                   & CH3O2_NO2_2_CH3O2NO2, CH3O2NO2_2_CH3O2_NO2
    real :: massLowThreshold
    real :: O3_decay_dark
    real :: conv_SO2_2_SO4_aqua_basic, &
          & NO2_2_PAN_rate, &
          & N2O5_2_HNO3, &
          & HCHO_NO3_2_HNO3, &
          & NO2_decay_sun_srf, NO2_decay_sun_25km, &
          & NO3_2_NO2_sun_srf, NO3_2_NO2_sun_25km, &
          & NO3_2_NO_sun_srf, NO3_2_NO_sun_25km, &
          & HONO_sun_srf, HONO_sun_25km, &
          & H2O2_sun_srf, H2O2_sun_25km, &
          & HCHO_sun_2_CO_HO2_srf, HCHO_sun_2_CO_HO2_25km,  &
          & HCHO_sun_2_CO_H2_srf, HCHO_sun_2_CO_H2_25km,  &
          & CH3OOH_sun_2_HCHO_srf, CH3OOH_sun_2_HCHO_25km, &
          & OH_NO3_2_NO2, CH3OOH_OH_2_HCHO, &
          & O3_sun_2_O_srf, O3_sun_2_O_25km, O3_sun_2_O1d_srf, O3_sun_2_O1d_25km, O1d_H2, O1d_H2O, &
          & HNO3_sun_2_OH_NO2_25km, HNO3_sun_2_OH_NO2_srf,&
          & PAN_radical_NO_loss, HO2_NO3_2_HNO3
    ! Load climatologicaly OH / used a fixed parameterisation
    logical :: use_clim_oh = .false.
    type(silja_logical) :: defined
  end type Tchem_rules_acidBasic
  PUBLIC Tchem_rules_acidBasic

  type int_ptr
     integer, pointer :: ptr
  end type int_ptr


  !--------------------------------------------------------------------------
  !
  ! Acid basic chemistry knows:
  !    SO2, SO4, H2SO4, 
  !    NH3, NH4_S, NH4_N, 
  !    NO, NO2, NOx for emission, NO3, HNO3, PAN
  !    O3
  !    C5H8-isoprene, C5H8_2 - monoterpenes
  ! We will point them by ind<name>
   integer, target, private, save :: iSO2_1d, & ! iH2SO4_1d, iNH3_1d, iNH4_S_1d, iNH4NO3_1d, &
                          & iNO_1d, iNO2_1d, iNO3rad_1d, iHNO3_1d, iPAN_1d, &
!                          & iPANrad_1d, &
                          & iNOP_1d, iCO_1d, iCH3OOH_1d, iHCHO_1d, iH2O2_1d, iCH3O2_1d, iHONO_1d, &
!                          & iC5H8_1d, iC5H8_2_1d, &
                          & nAcidBasicTranspSpecies !, &
!                          & iSO2_subst, iSO4_subst, iH2SO4_subst, & !iNH3_subst, iNH4_S_subst, iNH4NO3_subst, &
!                          & iNO_subst, iNO2_subst, iNO3rad_subst, iHNO3_subst, &
!                          & iPAN_subst, iPANrad_subst, iO3_subst, iCO_subst, iCH3OOH_subst, &
!                          & iHCHO_subst, iH2O2_subst, iCH3O2_subst, iHONO_subst, iOH_subst, &
!                          & iHO2_subst, iC5H8_subst, iC5H8_2_subst, nAcidBasicTranspSubst
  integer, target, private, save :: iH2SO4_1d_sl, iSO4w_1d_sl, iOH_1d_sl, iHO2_1d_sl

  integer, parameter, private :: nSpeciesAcidBasic = 13 !17 !23
  integer, parameter, private :: nSpeciesAcidBasic_SL = 4
  
  ! The list of species in acid_basic. Filled by init_chemicals_acid_basic.
  type(silam_species), dimension(:), allocatable, target, private, save :: species_acid_basic, &
                                                                         & species_acid_basic_SL

  ! List of integer pointers set to point at the index variable
  ! (above) corresponding to each species in the species_acid_basic array.
  type(int_ptr), dimension(:), allocatable, private, save :: species_indices, species_indices_SL


  ! Indices for meteorological quantities needed by the
  ! transformation. These point to the values in the meteo_input structure.
  integer, private, save, pointer :: ind_tempr, ind_tempr_2m, ind_press, ind_cloud_covr, &
                                   & ind_hgt, ind_spec_hum, ind_rel_hum, ind_ablh

  !integer, private, save :: indTempr2m, indFricVel, indVdCorrection
  
  !
  ! Stuff needed for computations of transformation and deposition.
  ! Pointers are set during prepare_DMAT_S_transdep
  !
  !  type(field_2d_data_ptr), private, save, POINTER :: fldTempr2m, fldSpecHumid2m, fldFricVel, &
  type(field_2d_data_ptr), private, save, POINTER :: fldTempr2m, fldFricVel, &
                                                   & fldVdCorrection, &
                                                   & fldTotalCloudCover, fldPrecLs, fldPrecCnv
  type(field_4d_data_ptr), private, save, pointer :: fldTempr, fldScavStd, &
                                                   & fldSpecHumid, fldRelHumidity, fldPressure

  real, dimension(:,:), private, pointer :: ptrLon, ptrLat, ptrCellSize
  real, DIMENSION(:), private,  POINTER :: ptrTotPrecip, ptrABLH, ptrLandFr

  integer, private, parameter :: nIndZ = 51
  integer, private, parameter :: nIndT = 301

  !
  ! Label of this module
  !  
  integer, parameter, public :: transformation_acid_basic = 5006


  integer, private, save :: barkcount = 0

CONTAINS


  !*******************************************************************

  subroutine init_chemicals_acid_basic()
    !
    ! Create the list of acid basic species. A trick is used for
    ! setting the indices: for each entry in species_acid_basic, an
    ! element in species_indices is set to point at the corresponding
    ! index variable (i*_1d). Then at registerSpecies, the indices are
    ! set to their actual values through the pointers. This avoids
    ! many lines of handwritten boilerplate code.
    !
    implicit none

    integer :: i, isp

    allocate(species_acid_basic(nSpeciesAcidBasic), species_indices(nSpeciesAcidBasic), &
           & species_acid_basic_SL(nSpeciesAcidBasic_SL), species_indices_SL(nSpeciesAcidBasic_SL), stat=i)
    if (i /= 0) then
      call set_error('Allocate failed', 'init_chemicals_acid_basic')
      return
    end if
    
    ! Transport species
    
    isp = 0
    call put_species('SO2', iSO2_1d, isp)
    call put_species('NO', iNO_1d, isp)
    call put_species('NO2', iNO2_1d, isp)
    call put_species('NO3rad', iNO3rad_1d, isp)
    call put_species('HNO3', iHNO3_1d, isp)
    call put_species('PAN', iPAN_1d, isp)
!    call put_species('PANrad', iPANrad_1d)
    call put_species('NOP', iNOP_1d, isp)
    call put_species('CO', iCO_1d, isp)
    call put_species('CH3O2', iCH3O2_1d, isp)
    call put_species('CH3OOH',iCH3OOH_1d, isp)
    call put_species('HCHO', iHCHO_1d, isp)
    call put_species('H2O2', iH2O2_1d, isp)
    call put_species('HONO', iHONO_1d, isp)
!    call put_species('NH3', iNH3_1d)
!    call put_species('H2SO4', iH2SO4_1d)
!    call put_species('NH4_S', iNH4_S_1d)
!    call put_species('NH4NO3', iNH4NO3_1d)
!    call put_species('C5H8', iC5H8_1d)
!    call put_species('C5H8_2', iC5H8_2_1d)

    ! Shortlived species
    isp = 0
    call put_species_SL('OH', iOH_1d_sl, isp, in_gas_phase)
    call put_species_SL('HO2', iHO2_1d_sl, isp, in_gas_phase)
    call put_species_SL('SO4', iSO4w_1d_sl, isp, in_gas_phase)   !fu_set_mode(fixed_diameter_flag, 0.5e-6, 2.0e-6, 1.11e-6))
    call put_species_SL('H2SO4', iH2SO4_1d_sl, isp, in_gas_phase)

    
  contains

    !====================================================================

    subroutine put_species(substname, ind_var, isp)
      implicit none
      character(len=*) :: substname
      integer, intent(inout) :: isp
      integer, target :: ind_var

      isp = isp + 1
      if (isp > nSpeciesAcidBasic) then
        call set_error('isp > nspeciesAcidBasic', 'put_species')
        return
      end if
      call set_species(species_acid_basic(isp), fu_get_material_ptr(substname), in_gas_phase)
      species_indices(isp)%ptr => ind_var
      
    end subroutine put_species

    !====================================================================

    subroutine put_species_SL(substname, ind_var, isp, mode)
      implicit none
      character(len=*) :: substname
      integer, intent(inout) :: isp
      integer, target :: ind_var
      type(Taerosol_mode), intent(in) :: mode

      isp = isp + 1
      if (isp > nSpeciesAcidBasic_SL) then
        call set_error('isp > nspeciesAcidBasic_SL', 'put_species_SL')
        return
      end if
      call set_species(species_acid_basic_SL(isp), fu_get_material_ptr(substname), mode)
      species_indices_SL(isp)%ptr => ind_var
      
    end subroutine put_species_SL

  end subroutine init_chemicals_acid_basic


  !************************************************************************************

  subroutine inventory_acid_basic(rules, &
                                & speciesEmis, speciesTransp, speciesShortlived, speciesAerosol,&
                                & nSpeciesEmis, nspeciesTransp, nspeciesShortlived, nspeciesAerosol, &
                                & iClaimedSpecies)
    !
    ! Essentially, return a pointer to the acid_basic chemicals. Take
    ! care not to change it!
    implicit none
    type(Tchem_rules_acidBasic), intent(in) :: rules
    type(silam_species), dimension(:), pointer :: speciesEmis, speciesTransp, speciesShortlived, &
                                                & speciesAerosol
    integer, intent(in) :: nSpeciesEmis
    integer, intent(inout) :: nSpeciesTransp, nSpeciesShortlived, nSpeciesAerosol
    integer, dimension(:), intent(inout) :: iClaimedSpecies

    ! Local variables
    integer :: iEmis

!    speciesTransp => species_acid_basic
!    nSpeciesTransp = nSpeciesAcidBasic
!    speciesShortLived => species_acid_basic_SL
!    nSpeciesShortlived = nSpeciesAcidBasic_SL
!    nullify(speciesAerosol)
!    nSpeciesAerosol = 0

    call addSpecies(speciesTransp, nSpeciesTransp, species_acid_basic, nSpeciesAcidBasic, .true.)
    call addSpecies(speciesShortlived, nSpeciesShortlived, species_acid_basic_SL, nSpeciesAcidBasic_SL, &
                  & .true.)
    if(error)return

    !
    ! Presence of emission species in the list of transported acid basic species means ownership
    ! Presence of emission species in short-living acid basic species means error
    !
    do iEmis = 1, nSpeciesEmis
      if(fu_index(speciesEmis(iEmis), species_acid_basic, nSpeciesAcidBasic) > 0)then
        if(iClaimedSpecies(iEmis) < 0)then
          call msg('Acid basic owns:' + fu_substance_name(speciesEmis(iEmis)))
          iClaimedSpecies(iEmis) = transformation_acid_basic
        else
          call msg('Cannot claim ownership because he claimed it already:',iClaimedSpecies(iEmis))
          call set_error('Cannot claim ownership because someone claimed it already','inventory_acid_basic')
          return
        endif
      endif  ! emission species belongs to acid basic transport one

      if(fu_index(speciesEmis(iEmis), species_acid_basic_SL, nSpeciesAcidBasic_SL) > 0)then
        call msg('This species are emitted to short-lived mass map:')
        call report(speciesEmis(iEmis))
        call set_error('This species are emitted to short-lived mass map:','inventory_acid_basic')
        return
      endif
    end do  ! emission species

  end subroutine inventory_acid_basic


  !************************************************************************************

  subroutine registerSpeciesAcidBasic(rules, &
                                    & speciesTransp, speciesShortlived, speciesAerosol,&
                                    & nspeciesTransp, nspeciesShortlived, nspeciesAerosol)
    !
    ! The final step in chemical init: collect the indices of acid
    ! basic species in the run. The species_indices pointer array is used.
    implicit none
    type(Tchem_rules_acidBasic), intent(in) :: rules
    type(silam_species), dimension(:), pointer :: speciesTransp, speciesShortlived, speciesAerosol
    integer, intent(in) :: nspeciesTransp, nspeciesShortlived, nspeciesAerosol
    type(chemical_adaptor) :: adaptor

    integer :: i, ind
    
    
    print *, "Before adaptor", iNO_1d, iNO2_1d, iNO3rad_1d
    
    !
    ! Transport species
    !
    call create_adaptor(species_acid_basic, speciesTransp, adaptor)
    if(error)return
    do i = 1, nSpeciesAcidBasic
      species_indices(i)%ptr = adaptor%iSp(i)
    end do

    print *, "After adaptor", iNO_1d, iNO2_1d, iNO3rad_1d

    !
    ! Short lived species
    !
    call create_adaptor(species_acid_basic_SL, speciesShortLived, adaptor)
    if(error)return
    do i = 1, nSpeciesAcidBasic_SL
      species_indices_SL(i)%ptr = adaptor%iSp(i)
    end do
    !
    ! Aerosol belongs to the aerosol dynamics, not to the acid basic
    !
!    call create_adaptor(species_acid_basic_AER, speciesAerosol, adaptor)
!    if(error)return
!    do i = 1, nSpeciesAcidBasic_AER
!      species_indices_AER(i)%ptr = adaptor%iSp(i)
!    end do
    

!    do i = 1, nSpeciesTransp
!      ind = fu_index(species_acid_basic(i), speciesTransp)
!      if (ind < 0) then
!        call msg('Looking for:')
!        call report(species_acid_basic(i))
!        call set_error('Species lost in transport', 'registerSpeciesAcidBasic')
!        return
!      end if
!      species_indices(i)%ptr = ind
!    end do
!
!    do i = 1, nSpeciesShortLived
!      ind = fu_index(species_acid_basic_SL(i), speciesShortLived)
!      if (ind < 0) then
!        call msg('Looking for:')
!        call report(species_acid_basic_SL(i))
!        call set_error('Species lost in short lived', 'registerSpeciesAcidBasic')
!        return
!      end if
!      species_indices_SL(i)%ptr = ind
!    end do

    call destroy_adaptor(adaptor)
    
  end subroutine registerSpeciesAcidBasic


  !*****************************************************************

  subroutine acidBasic_input_needs(rules, metdat)
    !
    ! Returns input needs for the transformation routines. 
    !
    implicit none

    ! Imported parameters
    type(Tchem_rules_acidBasic), intent(in) :: rules
    type(Tmeteo_input), intent(out), target :: metdat

    ! Local variables
    integer :: iStat, nq

    if (.not. rules%defined == silja_true) then
      call set_error('Undefined acid basic rules', 'acidBasic_input_needs')
      return
    end if

    metdat%nQuantities = 8

    metdat%quantity(1) = temperature_flag
    metdat%q_type(1) = meteo_dynamic_flag
    ind_tempr => metdat%idx(1)

    metdat%quantity(2) = temperature_2m_flag
    metdat%q_type(2) = meteo_dynamic_flag
    ind_tempr_2m => metdat%idx(2)

    metdat%quantity(3) = pressure_flag
    metdat%q_type(3) = meteo_dynamic_flag
    ind_press => metdat%idx(3)

    metdat%quantity(4) = total_cloud_cover_flag
    metdat%q_type(4) = meteo_dynamic_flag
    ind_cloud_covr => metdat%idx(4)

    metdat%quantity(5) = height_flag
    metdat%q_type(5) = meteo_dynamic_flag
    ind_hgt => metdat%idx(5)
    
    metdat%quantity(6) = specific_humidity_flag
    metdat%q_type(6) = meteo_dynamic_flag
    ind_spec_hum => metdat%idx(6)

    metdat%quantity(7) = relative_humidity_flag
    metdat%q_type(7) = meteo_dynamic_flag
    ind_rel_hum => metdat%idx(7)

    metdat%quantity(8) = abl_height_m_flag
    metdat%q_type(8) = meteo_dynamic_flag
    ind_ablh => metdat%idx(8)

  end subroutine acidBasic_input_needs


  !***********************************************************************

  subroutine set_miss_chem_rules_acidBasic(rules)
    implicit none
    type(Tchem_rules_acidBasic), intent(out) :: rules

    rules%defined = silja_false
    rules%vdSO2 = real_missing
    rules%vdSO4 = real_missing
    rules%vdNO2 = real_missing  
    rules%vdPAN = real_missing
    rules%vdNO3_h = real_missing
    rules%vdNH3 = real_missing
    rules%vdHNO3 = real_missing
    rules%massLowThreshold = real_missing

    nullify(rules%SO2_OH_2_SO3, &  !conv_SO2_2_SO4, &
          & rules%NO_2_NO2_with_O3, rules%NO2_2_NO3_dark, &
          & rules%NO_NO3_2_NO2, &
          & rules%PANrad_2_PAN_eq, rules%NH3_HNO3_vs_NH4NO3_eq, rules%HO2_PANrad_2_CH3O2, &
          & rules%NO3_NO2_N2O5_eq, rules%NO2_OH_2_HNO3, rules%HONO_OH_2_NO2, &
          & rules%CH4_OH_2_CH3O2, rules%CH3O2_HO2_2_CH3OOH, &
          & rules%CH3OOH_OH_2_CH3O2, rules%CH3O2_NO_2_HCHO_NO2, &
          & rules%CH3O2_NO2_2_CH3O2NO2, rules%CH3O2NO2_2_CH3O2_NO2, &
          & rules%CH3O2_CH3O2_2_HCHO, rules%CH3O2_CH3O2_2_CH3O, &
          & rules%CO_OH_2_CO2, rules%HCHO_OH_2_CO, &
          & rules%HO2_NO_2_NO2, rules%HO2_HO2_2_H2O2_basic, &
          & rules%HO2NO2_2_HO2_NO2, rules%HO2_NO2_2_HO2NO2, rules%NO_OH_2_HONO, &
          & rules%HO2_HO2_H2O_2_H2O2, rules%H2O2_OH_2_HO2, rules%HO2_O3_2_OH, rules%OH_O3_2_HO2, &
          & rules%O1d_N2, rules%O1d_O2, &
          & rules%CH4_cnc_global, rules%CH4_vert_prof)
    
    rules%CH3OOH_sun_2_HCHO_srf = real_missing
    rules%CH3OOH_sun_2_HCHO_25km = real_missing
    rules%CH3OOH_OH_2_HCHO = real_missing

    rules%HCHO_sun_2_CO_HO2_srf = real_missing
    rules%HCHO_sun_2_CO_HO2_25km = real_missing
    rules%HCHO_sun_2_CO_H2_srf = real_missing
    rules%HCHO_sun_2_CO_H2_25km = real_missing

    rules%H2O2_sun_srf = real_missing
    rules%H2O2_sun_25km = real_missing

    rules%HONO_sun_srf = real_missing
    rules%HONO_sun_25km = real_missing

    rules%O3_decay_dark = real_missing
    rules%O3_sun_2_O_srf = real_missing
    rules%O3_sun_2_O_25km = real_missing
    rules%O3_sun_2_O1d_srf = real_missing
    rules%O3_sun_2_O1d_25km = real_missing
    rules%O1d_H2O = real_missing
    rules%O1d_H2 = real_missing
    !
    ! Basic NOx
    !
    rules%NO_2_NO2_direct = real_missing
    rules%NO2_decay_sun_25km = real_missing
    rules%HNO3_sun_2_OH_NO2_srf = real_missing
    rules%HNO3_sun_2_OH_NO2_25km = real_missing
    rules%NO2_decay_sun_srf = real_missing
    rules%NO2_2_PAN_rate = real_missing
    rules%N2O5_2_HNO3 = real_missing
    rules%HCHO_NO3_2_HNO3 = real_missing
    rules%NO3_2_NO2_sun_srf = real_missing
    rules%NO3_2_NO2_sun_25km = real_missing
    rules%NO3_2_NO_sun_srf = real_missing
    rules%NO3_2_NO_sun_25km = real_missing
    rules%OH_NO3_2_NO2 = real_missing
    rules%massLowThreshold = real_missing
    rules%PAN_radical_prod = real_missing
    rules%PAN_radical_NO_loss = real_missing
    rules%HO2_NO3_2_HNO3 = real_missing

  end subroutine set_miss_chem_rules_acidBasic


  !***********************************************************************

  subroutine set_chem_rules_acidBasic(nlSetup, nlStdSetup, rules)
    !
    ! Creates the set of rules, in accordance with which the computations are to
    ! be made. Input - SILAM namelist
    ! Taking this chance, let's also set once-and-forever a few chemical coefficients
    !
    implicit none

    ! Imported parameter
    type(Tsilam_namelist), intent(in) :: nlSetup, nlStdSetup
    type(Tchem_rules_acidBasic), intent(out) :: rules

    ! Local variables
    integer :: iTmp, jTmp, ilev
    real :: fTempr, z, ratio_K_prod, ratio_K_loss, f300_T, fRelatDensity, fRelatPress, fRelatTempr,&
         & fprec
    real, parameter :: vres = 250.0 ! vertical resolution the OH lookup table is interpolated into
    character(len=fnlen) :: filename
    type(silam_vertical) :: vert_lut
    character(len=*), parameter :: sub_name = 'set_chem_rules_acidBasic'
    character(len=fnlen) :: nl_content
    !
    ! Stupidity checking
    !
    rules%defined = silja_false
    if(.not.defined(nlSetup))then
      call set_error('Undefined namelist given', sub_name)
      return
    endif

    !
    ! Dry deposition for SO2 and SO4 - so far we just take them as they are in DMAT.
    ! Later one can take something more sophisticated - e.g. link the whole stuff
    ! to resistance laws, etc. So far it is not supported.
    ! Another problem is that SILAM does not have the near-surface concentration profile,
    ! which essentially limits the dry deposition velocity. Therefore, I might need to
    ! cut it manually - and risk to over-state the absolute concentration levels.
    !
    rules%vdSO2 =  0.002   ! in DMAT: 0.002    [m/s]
    rules%vdSO4 =  0.0005  ! in DMAT : 0.0002  [m/s]
    rules%vdNO2 =  0.00083  
    rules%vdPAN =  0.00056
    rules%vdNO3_h =  0.0005  ! same as SO4 but again twice the DMAT value
    rules%vdNH3 =  0.00033 ! should be high but compensation point forces this
    rules%vdHNO3 = 0.0083

    if(fu_str_u_case(fu_content(nlSetup,'if_full_acid_chemistry')) == 'YES')then
      rules%ifNitrogen = .true.
    elseif(fu_str_u_case(fu_content(nlSetup,'if_full_acid_chemistry')) == 'NO')then
      rules%ifNitrogen = .false.
    else
      call set_error('Strange value of namelist item if_full_acid_chemistry:' + &
                   & fu_content(nlSetup,'if_full_acid_chemistry'), &
                   & sub_name)
      return
    endif

    !
    ! SO2 to SO4 conversion as a fraction of the SO2 mass per second
    ! It is temperature-dependent and in DMAT is set as a tabulated array following 
    ! the formula:
    ! TLA = 4e-4 * 1. * 1.1419**(T-261K) + 3e-3
    ! However, according to paper, it is: TLA = 138500 * exp(-4517/T)
    !
    ! It is not too far for low temperatures but MUCH smaller than that in the code for 
    ! higher temperatures. Therefore, I replace it with a similar-type exponent, which
    ! meets the power-law-based coef at lowest and highest temperatures. In-between it 
    ! is slightly larger (up to a factor of 2).
    ! 
    ! TLA = 4.2e7 * exp(-6000/T)  [1/hour] = 1.17e4 * exp(-6000/T)  [1/sec]
    !
    ! We will set nIndT=301 value from -200 to +100 degrees Celsius, which should be enough
    !
    ! Tuning the scheme:
    ! reduced conversion rate: 1e4 * exp(-6000/T)
    !
    ! TLA_N = 0.65 * T - 143, according to DMAT code. Paper has something different and strange
    !
    allocate(rules%NO_2_NO2_with_O3(nIndT), &
           & rules%SO2_OH_2_SO3(nIndT,nIndZ),  &
           & rules%CH4_cnc_global(nIndZ),      rules%CH4_vert_prof(nIndZ), &
           & rules%NO2_2_NO3_dark(nIndT),      rules%PAN_2_NO2_sun(nIndT), &
           & rules%PANrad_2_PAN_eq(nIndT,nIndZ),  rules%PAN_radical_prod(nIndT), &
           & rules%NH3_HNO3_vs_NH4NO3_eq(nIndT),  rules%NO_2_NO2_direct(nIndT), &
           & rules%NO_NO3_2_NO2(nIndT),        rules%NO3_NO2_N2O5_eq(nIndT), &
           & rules%CO_OH_2_CO2(nIndT,nIndZ),   rules%CH3O2_NO_2_HCHO_NO2(nIndT), &
           & rules%HO2_NO_2_NO2(nIndT),        rules%H2O2_OH_2_HO2(nIndT), &
           & rules%O1d_N2(nIndT,nIndZ),        rules%O1d_O2(nIndT,nIndZ), &
           & rules%HO2_O3_2_OH(nIndT),         rules%OH_O3_2_HO2(nIndT), &
           & rules%NO2_OH_2_HNO3(nIndT,nIndZ),     rules%NO_OH_2_HONO(nIndT,nIndZ), &
           & rules%HO2_NO2_2_HO2NO2(nIndT,nIndZ),  rules%HO2NO2_2_HO2_NO2(nIndT,nIndZ), &
           & rules%HO2_HO2_2_H2O2_basic(nIndT,nIndZ), rules%HO2_HO2_H2O_2_H2O2(nIndT), &
           & rules%CH3O2_HO2_2_CH3OOH(nIndT),         rules%CH4_OH_2_CH3O2(nIndT), &
           & rules%CH3O2_NO2_2_CH3O2NO2(nIndT,nIndZ), rules%CH3O2NO2_2_CH3O2_NO2(nIndT,nIndZ), &
           & rules%CH3O2_CH3O2_2_HCHO(nIndT),  rules%CH3O2_CH3O2_2_CH3O(nIndT), &
           & rules%CH3OOH_OH_2_CH3O2(nIndT),   rules%HCHO_OH_2_CO(nIndT), & 
           & rules%HO2_PANrad_2_CH3O2(nIndT),  rules%HONO_OH_2_NO2(nIndT), &
           & stat=iTmp)
    if(iTmp /= 0)then
      call set_error('Failed to allocate memory for conversion rates', sub_name)
      return
    endif

    !
    ! Cycle over temperature
    !
    do iTmp = 1, nIndT

      fTempr = real(iTmp) + 72.15
      f300_T = 300. / fTempr

!call msg('800. * exp(-6000. / fTempr',800. * exp(-6000. / fTempr))
!      rules%conv_SO2_2_SO4_aqua(iTmp) = 800. * exp(-6000. / fTempr)
!      rules%conv_SO2_2_SO4_aqua(iTmp) = 400. * exp(-6000. / fTempr)

      !  After Sander(2003), Also used by Seinfield & Pandis (2006):
      !
      rules%NO_2_NO2_with_O3(iTmp) = 1.8e6 * exp(-1500./fTempr)

      ! Direct NO->NO2 is the main O3 driving factor (seemingly)
      ! for CO it is ~5e-4, a replacement for VOC-driven direct oxidation
      ! Tunning:
      ! 5.7.2007 : 0.0005 as for CO - far too small
      ! 6.7.2007 morning : 0.002*(1.+sunD) - too big
      ! 6.7.2007 evening: 0.002*(sunD+0.01) day and 0.002*0.01 at night
      ! 13.7.2007: temperature dependence is introduced as K*exp(-E/RT) with K=0.01, E/R=300
      !            for T=288 that would correspond to previous 0.002
      !
      ! __.7.2007 : 0.002 was substituted by 0.006, tunning 
      ! 27.7.2007 : 0.006 was substituted by 0.009, keep tunning is still too low 
      ! 01.8.2007 : 0.006 is again used when the ozone night production is being tunned
      ! 29.08.2007: NEW SILAM version, 0.006 seams to be the best coefficient for now
      ! 30.10.2007: 4-fold increase. Further makes no sense: NOx-limited case
      ! 13.12.2007: roll back the increase. CO and CH4 seem to take care of large
      !             part of the transformation. Ozone is high even in December. Run 27
      ! 19.12.2007: 10-fold cut. Ozone production is tremendous - in December there must be
      !             titration everywhere in industrial regions. Background ozone is a problem

      rules%NO_2_NO2_direct(iTmp) = 4. * 0.0006 * exp(-f300_T)

!      rules%NO2_2_NO3_dark(iTmp) = 0.1557*(fTempr - 215.)     ! DMAT code
!      rules%NO2_2_NO3_dark(iTmp) = 0.157*(fTempr - 228.7)     ! DMAT paper
!      rules%NO2_2_NO3_dark(iTmp) = 7.2e4 * exp(-2450./fTempr) ! Atkinson 1997
      rules%NO2_2_NO3_dark(iTmp) = 1200. * exp(-1350./fTempr)  ! Atkinson 1997 fit to code

      rules%PAN_2_NO2_sun(iTmp) = 7.9444e13*exp(-12530./fTempr)  !DMAT code~paper~Simpson
! Wrong way: for 32 bit machines all zero at already -10C
!      rules%NH3_HNO3_vs_NH4NO3_eq(iTmp) = 1.12e22*(298./fTempr)**6.1 * exp(-24220./fTempr)
! Right way: all goes OK down to -100C
      rules%NH3_HNO3_vs_NH4NO3_eq(iTmp) = exp(50.7702 + 6.1*log(298./fTempr) - (24220./fTempr))
      rules%NO_NO3_2_NO2(iTmp) = 9.e6*exp(170./fTempr)            ! S&P

      !No overflows!
      rules%NO3_NO2_N2O5_eq(iTmp) = 1.81e-9 * exp(min(10990./fTempr, 60.))  ! Atkinson 1997
      rules%HO2_NO_2_NO2(iTmp) = 2.23e6 * exp(240./fTempr)        ! n.9, p.18
      rules%H2O2_OH_2_HO2(iTmp) = 1.747e6 * exp(-160./fTempr)     ! n.9,p.18
      rules%CH3O2_HO2_2_CH3OOH(iTmp) = 2.47e5 * exp(750./fTempr)  ! n.9,p.23
      rules%CH4_OH_2_CH3O2(iTmp) = 1.48e6 * exp(-1775./fTempr)    ! n.9,p.23
      rules%CH3OOH_OH_2_CH3O2(iTmp) = 2.17e6 * exp(200./fTempr)   ! n.9,p.23
      rules%CH3O2_NO_2_HCHO_NO2(iTmp) = 1.69e6 * exp(f300_T)      ! n.9,p.23
      rules%CH3O2_CH3O2_2_CH3O(iTmp) = 5.9e-13*exp(-509./fTempr) * 6.025e17 !n.10,p3
      rules%CH3O2_CH3O2_2_HCHO(iTmp) = 7.04e-14*exp(365./fTempr) * 6.025e17 !n.10,p3
      rules%HCHO_OH_2_CO(iTmp) = 5.18e6 * exp(20./fTempr)         ! Jakobson
      rules%HO2_HO2_H2O_2_H2O2(iTmp) = 8.44e-4 * exp(2200./fTempr) ! n.9.p.18
      rules%HO2_O3_2_OH(iTmp) = 8.43e3 * exp(-600./fTempr)         ! Jakobson
      rules%OH_O3_2_HO2(iTmp) = 1.14e6 * exp(-1000./fTempr)        ! Jakobson
      rules%HONO_OH_2_NO2(iTmp) = 1.08e7 * exp(-390./fTempr)       ! Jakobson
      rules%HO2_PANrad_2_CH3O2(iTmp) = 1.9e5 * exp(1040./fTempr)   ! Jakobson

      !
      ! PAN is always taken in equilibrium as a ratio of production and loss (notebook 8 p.77)
      ! Altitude dependence for [M] is taken from: [M]=[M0]*exp(-H_km/7.4). Error >10% for H>15km
      ! For N2 and O2 separately, the magic 7.4km is replaced with 7.66 and 6.7 km, respectively.
      ! See notebook 9, p.28 for details.
      !
      ! ATTENTION. UNITS are NOT SI
      !
      rules%PAN_radical_prod(iTmp) = 0.9e-14*exp(-f300_T)

      do jTmp = 1, nIndZ
        z = ((jTmp - 1) * 500.) / 7400.            ! a ratio to 7.4km scale
        
        call us_standard_atmosphere(real(jTmp-1)*500., fRelatDensity, fRelatPress, fRelatTempr)
        if(error)return

        ratio_K_prod = 575. * (f300_T)**6.2 * fRelatDensity     ! K0_3*M / Kinf_3
        ratio_K_loss = 2.3 * exp(1730./fTempr) * fRelatDensity             ! K0_4*M / Kinf_4
!        ratio_K_prod = 575. * (f300_T)**6.2 * exp(-z)     ! K0_3*M / Kinf_3
!        ratio_K_loss = 2.3 * exp(1730./fTempr - z)             ! K0_4*M / Kinf_4
        rules%PANrad_2_PAN_eq(iTmp, jTmp) = &
              & (6.9e-9 * (f300_T)**7.1 * fRelatDensity / (1.+ratio_K_prod) * &      ! PAN production
!              & (6.9e-9 * (f300_T)**7.1 * exp(-z) / (1.+ratio_K_prod) * &      ! PAN production
               & 0.3**(1./(1.+(log(ratio_K_prod)*ln_10_1)**2))) / &
              & (1.25e17 * exp(-12100./fTempr) * fRelatDensity / (1.+ratio_K_loss) * & ! PAN destroying
!              & (1.25e17 * exp(-12100./fTempr - z) / (1.+ratio_K_loss) * &     ! PAN destroying
                     & 0.3**(1./(1.+(log(ratio_K_loss)*ln_10_1)**2)) + &
               & (7.23e-8*(1.-real(jTmp-1)/50.) + 3.54e-6*(real(jTmp-1)/50.))) ! photo destroying

        rules%NO2_OH_2_HNO3(iTmp,jTmp) = (1.75e8*fRelatDensity*(f300_T)**3.5) / &     ! n.9,p.31,nomin
              & (4.4*(f300_T)**0.6 + 3.98*fRelatDensity*(f300_T)**2.9) * &            ! denominator
              & 0.43 ** (1./(1.+(log(0.905*fRelatDensity*(f300_T)**2.3)*ln_10_1)**2)) ! broadening

        rules%SO2_OH_2_SO3(iTmp,jTmp) = (7.42e6*fRelatDensity*(f300_T)**3.3) / &          ! n.9,p.32: nom
                   & (1.21 + 6.13*fRelatDensity*(f300_T)**3.3) * &                        ! denominator
                   & 0.45 ** (1./(1.+(log(5.07*fRelatDensity*(f300_T)**3.3)*ln_10_1)**2)) ! broadening

!        rules%NO2_OH_2_HNO3(iTmp,jTmp) = (1.75e8*exp(-z)*(f300_T)**3.5) / &       ! nominator
!              & (4.4*(f300_T)**0.6 + 3.98*exp(-z)*(f300_T)**2.9) * &              ! denominator
!              & 0.43 ** (1./(1.+(log(0.905*exp(-z)*(f300_T)**2.3)*ln_10_1)**2))   ! broadening

!        rules%SO2_OH_2_SO3(iTmp,jTmp) = (7.42e6*exp(-z)*(f300_T)**3.3) / &          ! n.9,p.32: nom
!                   & (1.21 + 6.13*exp(-z)*(f300_T)**3.3) * &                        ! denominator
!                   & 0.45 ** (1./(1.+(log(5.07*exp(-z)*(f300_T)**3.3)*ln_10_1)**2)) ! broadening

        rules%CO_OH_2_CO2(iTmp,jTmp) = 7.83e4 * (1+0.6*fRelatPress) * exp(f300_T) !n9,p18
        rules%O1d_N2(iTmp,jTmp) = 4.59e8 * 0.8 * exp(107./fTempr - (((jTmp-1) * 500.) / 7660.)) ! sec-1
        rules%O1d_O2(iTmp,jTmp) = 8.16e8 * 0.16 * exp(67./fTempr - (((jTmp-1) * 500.) / 6700.)) ! sec-1

        rules%HO2_HO2_2_H2O2_basic(iTmp,jTmp) = 1.39e5 * exp(600./fTempr) + &
                                              & 2.61e4 * exp(1000./fTempr) * fRelatDensity
!                                              & 2.61e4 * exp(1000./fTempr - z)

        rules%HO2NO2_2_HO2_NO2(iTmp,jTmp) = 42.3 * fRelatDensity * (2.83e33 * exp(-10900/fTempr)) / & !n.9,p.38
                        & (1.8*42.3*fRelatDensity + 1.57e3 * exp(-900/fTempr)) * &
                        & 0.6**(1./(1.+(log(1.15e-3*42.3*fRelatDensity*exp(900/fTempr))*ln_10_1)**2))

        rules%HO2_NO2_2_HO2NO2(iTmp,jTmp) = 1.85e7 * 42.3 * fRelatDensity * (f300_T)**1.4 / & !n.9,p.39
                        & (6.53 * 42.3 * fRelatDensity + 2.83e2 * (f300_T)**0.2) * &
                        & 0.6**(1./(1.+(log(2.31e-2*42.3*fRelatDensity*(fTempr/300.)**0.2)*ln_10_1)**2))

        rules%NO_OH_2_HONO(iTmp,jTmp) = 2.7e7 * 42.3 * fRelatDensity * (f300_T)**2.4 / & !n.9,p.39. NOTE DAYTIME PRESENCE
                        & (42.3*fRelatDensity*(f300_T)**2.4 + 100) * &
                        & 0.9**(1./(1.+(log(0.423*fRelatDensity*(f300_T)**2.4)*ln_10_1)**2))


!call msg('Drop HONO production')
!rules%NO_OH_2_HONO(iTmp,jTmp) = 0.0




        rules%CH3O2_NO2_2_CH3O2NO2(iTmp,jTmp) = 4.52e6 * 42.3*fRelatDensity * f300_T**5.5 / & !n.10,p.5
                        & (5. + 42.3*fRelatDensity * f300_T**5.5) * &
                        & 0.36**(1./(1.+(log(8.46*fRelatDensity * f300_T**5.5)*ln_10_1)**2))

        rules%CH3O2NO2_2_CH3O2_NO2(iTmp,jTmp) = 42.3*fRelatDensity * (1.1e16 * exp(-10560/fTempr)) / &
                        & (42.3*fRelatDensity + 203. * exp(-870./fTempr)) * &
                        & 0.36**(1./(1.+(log(4.93e-3 * 42.3*fRelatDensity * exp(870./fTempr))*ln_10_1)**2))
        !
        ! Here we set the relative vertical profile for CH4 using the year 2000 surface value
        ! Adjustment to a specific year is going on in prepare_transdep. 
        ! Both formulas are brute-force fitting to IPCC data and sounding profiles
        !
        if(iTmp == 1)then  ! methane we need only once
          if(fRelatPress > 0.3)then
            rules%CH4_vert_prof(jTmp) = 1.
          else
            rules%CH4_vert_prof(jTmp) = 1. - 3.25e4 * (0.35-fRelatPress)**10
          endif
          rules%CH4_vert_prof(jTmp) = rules%CH4_vert_prof(jTmp) * 4.09e-5 * fRelatDensity ! ppm -> mole/m3
        endif
      end do  ! height cycle

    end do  ! tempr cycle
    !fPrec = real(fu_year(fu_valid_time(fldTempr2m%present%idPtr))) + &
    !      & real(fu_mon(fu_valid_time(fldTempr2m%present%idPtr))) / 12.
    fprec = 2009.0
    rules%CH4_cnc_global(1:nIndZ) = rules%CH4_vert_prof(1:nIndZ) * &
              & (0.83 + 1.5e-3*(fPrec-1800) + 0.012 * (fPrec - 1930) / (1.+exp((1920-fPrec)/50.)))
    !
    ! Sulphur aqueous-phase conversion basic coefficient.
    !
!!!!!    rules%conv_SO2_2_SO4_aqua_basic = 2.e-6  ! =1000*exp(-6000/300), see notebook 11, pp.24-26
    rules%conv_SO2_2_SO4_aqua_basic = 1.e-6
    
    
    !
    ! Ozone.
    !
    rules%O3_decay_dark = 0. !4.e-8    !sec-1. 1.39e-6 for rate: 0.005 hr-1
    rules%O3_sun_2_O_srf = 4.17e-4   ! sec-1, to be scaled with sunD
    rules%O3_sun_2_O_25km = 4.95e-4  ! sec-1, to be scaled with sunD
    rules%O3_sun_2_O1d_srf = 5.08e-5 ! sec-1, to be scaled with sunD
    rules%O3_sun_2_O1d_25km = 1.13e-4  ! sec-1, to be scaled with sunD

    rules%O1d_H2O = 1.32e8   ! m3 mole-1 sec-1
    rules%O1d_H2  = 6.63e7   ! m3 mole-1 sec-1

    !
    ! Nitrogen oxides
    !
    rules%NO2_2_PAN_rate = 2.583333e-5 ! 0.093/3600
    rules%N2O5_2_HNO3 = 1.2e-3         ! hydrolysis, m3 mole-1 sec-1
    rules%OH_NO3_2_NO2 = 1.2e7         ! Jacobson
    rules%HO2_NO3_2_HNO3 = 2.41e6      ! Jakobson

    rules%NO2_decay_sun_25km = 0.0124  ! Jacobson, sun-related NO2 decay stratosph (paper: 12 hr-1)
    rules%NO2_decay_sun_srf = 0.00882  ! Jacobson, sun-related NO2 decay surface (paper: 12 hr-1)
    rules%NO3_2_NO2_sun_srf = 0.284    ! sec-1
    rules%NO3_2_NO2_sun_25km = 0.304   ! sec-1
    rules%NO3_2_NO_sun_srf = 0.0249    ! sec-1
    rules%NO3_2_NO_sun_25km = 0.0265   ! sec-1

    !
    ! methane, formaldehyde, and CO
    !
    rules%CH3OOH_OH_2_HCHO = 1.145e6       ! m3 mole-1 sec-1
    rules%HCHO_NO3_2_HNO3 = 348.0          ! m3 mole-1 sec-1

    rules%H2O2_sun_25km = 1.18e-5          ! sec-1
    rules%H2O2_sun_srf = 7.72e-6           ! sec-1
    rules%HONO_sun_25km = 2.75e-3          ! sec-1
    rules%HONO_sun_srf = 1.94e-3           ! sec-1
    rules%HNO3_sun_2_OH_NO2_25km = 5.61e-6 ! sec-1
    rules%HNO3_sun_2_OH_NO2_srf = 8.23e-7  ! sec-1
    rules%CH3OOH_sun_2_HCHO_25km = 9.92e-6 ! sec-1
    rules%CH3OOH_sun_2_HCHO_srf = 5.73e-6  ! sec-1
    rules%HCHO_sun_2_CO_HO2_25km = 6.17e-5 ! sec-1
    rules%HCHO_sun_2_CO_HO2_srf = 3.29e-5  ! sec-1
    rules%HCHO_sun_2_CO_H2_25km = 7.64e-5  ! sec-1
    rules%HCHO_sun_2_CO_H2_srf = 4.4e-5    ! sec-1

    rules%PAN_radical_NO_loss = 1.2e6  ! m3 mole-1 sec-1. Reaction with NO
    !
    ! For Eulerian scheme we will need the minimum threshold of amount of species in the grid cell
    ! The cell will not be processed by ANY Eulerian-pool routine should it have smaller amount 
    ! than that threshold. A data-based setting is possible only when the total emitted amount
    ! is known, as well as a grid size. So far let's put it to something small enough.
    !
    rules%massLowThreshold = 1.0e-12

    nl_content = fu_str_l_case(fu_content(nlsetup, 'oh_param_method'))
    select case(nl_content)
    case ('climatology')
      rules%use_clim_oh = .true.
    case ('parametrization')
      rules%use_clim_oh = .false.
    case default
      call msg_warning('Missing oh_param_method, take parameterization')
      rules%use_clim_oh = .false.
    end select
    
    if (rules%use_clim_oh) then
      filename = fu_content(nlStdSetup, 'oh_climatology_file')
      if (fu_fails(filename /= '', 'Missing oh_climatology_file', sub_name)) return
      call set_vertical((/(fu_set_level(constant_height, ilev*vres), ilev = 0, int(20e3/vres))/), vert_lut)
      call setup_oh_lut(fu_process_filepath(filename), vert_lut)
      call set_missing(vert_lut, ifNew=.false.)
      if (error) return
    end if

    rules%defined = silja_true

  end subroutine set_chem_rules_acidBasic


  !***********************************************************************

!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***
!!$  subroutine transf_acidBasic_single_val(vSp, garbage, now, timestep_sec, weight_past, &
!!$                                       & ix, iy, iLev, fHeight, fCellVolume, &
!!$                                       & rules, cocktail_template, &
!!$                                       & interpCoefMet2DispHoriz, interpCoefMet2DispVert, &
!!$                                       & ifInterpMet2DispHoriz, ifInterpMet2DispVert)

!  subroutine transf_acidBasic_single_val(vSp, garbage, now, timestep_sec, weight_past, &
!                                       & ix, iy, iLev, fHeight, fCellVolume, &
!                                       & fTemprRef, fTemprSpan, fCloudCover, fLon, fLat, fSpecHumid, uFileOH, &
!                                      & rules, cocktail_template, &
!                                       & interpCoefMet2DispHoriz, interpCoefMet2DispVert, &
!                                       & ifInterpMet2DispHoriz, ifInterpMet2DispVert)
!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***

  subroutine transform_acid_basic(vSp, vSpSL, rules, metdat, timestep_sec, garbage, &
                               & zenith_cos, species_list, fLowCncThresh, ifReport, ix, iy, iz, &
                               & lat, lon, now, print_it)
    !
    ! Implements chemical conversion for SOx, NOx and NHx as it is done in DMAT.
    ! Note that ozone is an internal substance, which is computed regardless the NOx/VOC
    ! concentrations and thus must not be given to the output (well, it may but must not be 
    ! shown widely). The only use for ozone is to keep the oxidation capacity and NO-NO2
    ! ratio.
    ! Later, this equilibrium (probably the largest hole in the mechanism) must be replaced
    ! with standard NO-NO2-VOC-O3 cycle. To start, VOC can be treated as ozone creation potential
    !
    implicit none

    ! Imported parameters
    real, dimension(:), intent(inout) :: vSp, vSpSL
    type(Tchem_rules_acidBasic), intent(in) :: rules
    real, dimension(:), intent(in) :: metdat
    real, intent(in) :: timestep_sec
    real, dimension(:), intent(inout) :: garbage
    real, intent(in) :: zenith_cos
    type(silam_species), dimension(:), intent(in) :: species_list
    real, dimension(:)  :: fLowCncThresh
    logical, intent(in) :: ifReport   ! a request of report from the calling sub
    integer, intent(in) :: ix, iy, iz
    real, intent(in) :: lat, lon ! for OH concentration lookup
    type(silja_time), intent(in) :: now
    logical, intent(out) :: print_it  ! a request of report to the calling sub
    
    ! Local variables
!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***
!    real :: tempr, sun, sunD, dO3, dNO, fNO_eq, fNO2_eq, fO3_eq, z, fTmp, fRate, fRelatDensity, fRelatPress, fRelatTempr, &
    real :: tempr, fCloudCover, sun, sunD, dO3, dNO, fNO_eq, fNO2_eq, fO3_eq, z, fTmp, fRate, fPressure, &
!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***
          & rate_air, rate_water, scaling, dHNO3, DNO2_2_PAN, dSO2, dH2SO4_SO4, dNO2, dPAN, &
          & DHNO3_NO3, CEQUIL, sqrtmp, differ, TNH3, THNO3, dNH3_4, fMixing_ratio_2_cnc, &
          & fHeight, dSO2_water, dSO2_air, &
          & fPANprod, fEqRatio, cH2, cH2O, cN2O5, &
          & cOH_forced, cHO2NO2, fK_HO2_HO2,&
          & fJ_NO2, fJ_NO3_2_NO,  fJ_NO3_2_NO2, fJ_HONO, fJ_O3_O, fJ_O3_O1d, fJ_H2O2, &
          & fJ_CH3OOH, fJ_HCHO_2_HO2, fJ_HCHO_2_H2, fJ_HNO3, &
          & pCH3O2, dCO, a, b, aP_OH, bP_OH, aL_OH, aP_HO2, bP_HO2, aL_HO2, dCH3x, dHCHO, dH2O2, &
          & detBB, det4AC

    integer :: indT, indZ, iSubst, iMode, nSteps, iStep, iHourTmp, jTmp, nSpecies
    real, parameter :: fOH_night_vs_day = 0.01
    real, parameter :: fSunMin = 0.01 !Minimum cos of Solar Zenith angle for day

    nSpecies = size(vSp)


   vSpSL(iH2SO4_1d_sl) = 0.
   vSpSL(iSO4w_1d_sl) = 0.
   vSpSL(iOH_1d_sl) = 0. 
   vSpSL(iHO2_1d_sl) = 0.

!  call msg('Acid transform start')
!  call report_vector('before')
    !--------------------------------------------------------------------------------
    !
    ! Get temperature from the meteo buffer
    !
!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***
!    sunD = fu_solar_zenith_angle_cos(fLon, fLat, now) * (1. - fCloudCover*0.5)
!
!    tempr = fTemprRef + fTemprSpan * (1+max(0.,sunD))
!
!    call us_standard_atmosphere(fHeight, fRelatDensity, fRelatPress, fRelatTempr)
!
!    fMixing_ratio_2_cnc = fRelatPress * std_pressure_sl / (gas_constant_uni * tempr)
!
!
    tempr = metdat(ind_tempr)
    fCloudCover = metdat(ind_cloud_covr)
    fHeight = metdat(ind_hgt)

    sun = zenith_cos
    if(sun < 0.) sun = 0.  ! night
    sunD = sun * (1. - fCloudCover*0.5)
    if (.not. fCloudCover <= 1.)then
      call set_error('Negative sunshine', 'transform_acid_basic') 
      call msg('sun, fCloudCover', sun, fCloudCover)
    endif
    if(error)return
!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***

    indT = int(tempr-71.65)

    indZ = int(fHeight / 500.) + 1

    !iHourTmp = fu_hour(now)

    !
    ! ATTENTION.
    ! All computations must go in concentrations while advection et al require masses.
    ! Thus, we divide the mass vector with the cell volume. At the end we return it back
    !
!    vSp(:) = vSp(:) / fCellVolume
!    vSpSL(:) = vSpSL(:) / fCellVolume

    !
    ! So far, we do not believe the equilibrium OH-HO2 model, thus will
    ! use OH_forced mass in the grid cell. For sulphur-only chemistry, it is the only 
    ! available option.
    ! 
    ! Forced OH can be either from the tabulated climatology or the EMEP parametrisation.
    if (rules%use_clim_oh) then
      cOH_forced = fu_OH_cnc_clim(lat, lon, fHeight, now)
    else
      cOH_forced = fu_OH_cnc_forced(sunD, fHeight)
    end if


    if(rules%ifNitrogen)then
      !
      ! Nitrogen chemistry. 
      ! Solve two budget equations for masses M_OH and M_HO2 (notebook 9, p.35)
      ! aP_OH  + bP_OH  * M_HO2 = M_OH * aL_OH
      ! aP_HO2 + bP_HO2 * M_OH  = K_HO2_HO2 * M_HO2**2 + aL_HO2 * M_HO2
      !
      ! Photodissociation rates
      !
      if(sunD > fSunMin)then
        z = min(fHeight / 25000., 1.)
        fTmp = sun**0.2 * exp(-0.25/sun) * (1. - fCloudCover*0.5)

        fJ_NO2 = fTmp * (rules%NO2_decay_sun_srf * (1.-z) + rules%NO2_decay_sun_25km * z)
        !
        ! NO3 is copied from NO2
        !
        fJ_NO3_2_NO = fTmp  * (rules%NO3_2_NO_sun_srf * (1.-z) + rules%NO3_2_NO_sun_25km * z)
        fJ_NO3_2_NO2 = fTmp * (rules%NO3_2_NO2_sun_srf * (1.-z) + rules%NO3_2_NO2_sun_25km * z)
        fJ_HONO = fTmp * (rules%HONO_sun_srf * (1.-z) + rules%HONO_sun_25km * z)
        fJ_HNO3 = fTmp * (rules%HNO3_sun_2_OH_NO2_srf * (1.-z) + rules%HNO3_sun_2_OH_NO2_25km * z)
        !
        ! O3_O is copied from O3_O1d
        !
        fTmp = sun**2.7245 * (1. - fCloudCover*0.5)

        fJ_O3_O = fTmp * (rules%O3_sun_2_O_srf * (1.-z) + rules%O3_sun_2_O_25km * z)
        fJ_O3_O1d = fTmp * (rules%O3_sun_2_O1d_srf * (1.-z) + rules%O3_sun_2_O1d_25km * z)
        fJ_CH3OOH = sun**0.764 * exp(-0.249/sun)  * (1. - fCloudCover*0.5) * &
                               & (rules%CH3OOH_sun_2_HCHO_srf * (1.-z) + rules%CH3OOH_sun_2_HCHO_25km * z)
        !
        ! H2O2 is copied from O3_o1d
        !
        fJ_H2O2 = fTmp * (rules%H2O2_sun_srf * (1.-z) + rules%H2O2_sun_25km * z )
        fJ_HCHO_2_HO2 = sun**0.781 * exp(-0.349/sun) * (1. - fCloudCover*0.5) * &
                         & (rules%HCHO_sun_2_CO_HO2_srf * (1.-z) + rules%HCHO_sun_2_CO_HO2_25km * z)
        fJ_HCHO_2_H2 = sun**0.565 * exp(-0.275/sun) * (1. - fCloudCover*0.5) * &
                           & (rules%HCHO_sun_2_CO_H2_srf * (1.-z) + rules%HCHO_sun_2_CO_H2_25km * z)
      else
        fJ_no2 = 0.
        fJ_O3_O = 0.
        fJ_O3_O1d = 0.
        fJ_CH3OOH = 0.
        fJ_H2O2 = 0.
        fJ_HONO = 0.
        fJ_HCHO_2_HO2 = 0.
        fJ_HCHO_2_H2 = 0.
        fJ_HNO3 = 0.
      endif

      !
      ! Other chemical coefficients
      !
!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***
      !
      ! specific humidity is in kg/kg, have to convert to mole/m3
      !
      fPressure = metdat(ind_press)
      fMixing_ratio_2_cnc = fPressure / (gas_constant_uni * tempr)
      
      cH2O = metdat(ind_spec_hum)

!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***

      cH2 = fMixing_ratio_2_cnc * 6.e-4    ! H2 concentration is 0.6 ppm

      fK_HO2_HO2 = rules%HO2_HO2_2_H2O2_basic(indT,indZ) * (1. + rules%HO2_HO2_H2O_2_H2O2(indT) * cH2O)

      cHO2NO2 = 0.

!      call report_vector('Start')
      call check_nan('Start')


      !
      ! Preparing the coefficients for the squared equation, do not miss the point: a-s and B-s 
      ! have different units: mole/m3sec and 1/sec, respectively
      ! See notebook 9, p.35
      !
      aP_OH = 2. * fJ_H2O2   * vSp(iH2O2_1d) + &                            ! H2O2
            & fJ_HONO * vSp(iHONO_1d) + &                                   ! HONO
            & fJ_HNO3 * vSp(iHNO3_1d) + &                                   ! HNO3
            & fJ_CH3OOH * vSp(iCH3OOH_1d) + &
            & fJ_O3_O1d * vSp(iNOP_1d) * (2. * rules%O1d_H2O * cH2O + rules%O1d_H2 * cH2) / &
                                      & (rules%O1d_N2(indT,indZ) + rules%O1d_O2(indT,indZ) + &
                                       & rules%O1d_H2O * cH2O + rules%O1d_H2 * cH2)

      bP_OH = rules%HO2_NO_2_NO2(indT) * vSp(iNO_1d) + rules%HO2_O3_2_OH(indT) * vSp(iNOP_1d)

      aL_OH = rules%CO_OH_2_CO2(indT,indZ)  * vSp(iCO_1d) + &
            & rules%CH4_OH_2_CH3O2(indT)    * rules%CH4_cnc_global(indZ) + &
            & rules%OH_O3_2_HO2(indT)       * vSp(iNOP_1d) + &
            & rules%NO_OH_2_HONO(indT,indZ) * vSp(iNO_1d) + &
            & rules%HONO_OH_2_NO2(indT)     * vSp(iHONO_1d) + &
            & rules%NO2_OH_2_HNO3(indT,indZ)* vSp(iNO2_1d) + &
            & rules%H2O2_OH_2_HO2(indT)     * vSp(iH2O2_1d) + &
            & rules%CH3OOH_OH_2_CH3O2(indT) * vSp(iCH3OOH_1d) + &
            & rules%CH3O2_CH3O2_2_HCHO(indT)* vSp(iCH3O2_1d) *  vSp(iCH3O2_1d)   ! see n.10,p.3

      aP_HO2 = 2. * fJ_HCHO_2_HO2 * vSp(iHCHO_1d) + &
             & fJ_CH3OOH * vSp(iCH3OOH_1d) + &
             & rules%HO2NO2_2_HO2_NO2(indT,indZ) * cHO2NO2 + &
             & rules%CH3O2_NO_2_HCHO_NO2(indT)   * vSp(iCH3O2_1d) * vSp(iNO_1d) + &
             & (2.*rules%CH3O2_CH3O2_2_CH3O(indT) + rules%CH3O2_CH3O2_2_HCHO(indT)) * &
                                                                  & vSp(iCH3O2_1d)*vSp(iCH3O2_1d)

      bP_HO2 = rules%OH_O3_2_HO2(indT)      * vSp(iNOP_1d) + &
             & rules%H2O2_OH_2_HO2(indT)    * vSp(iH2O2_1d) + &
             & rules%OH_NO3_2_NO2           * vSp(iNO3rad_1d) + &
             & rules%CO_OH_2_CO2(indT,indZ) * vSp(iCO_1d) + &
             & rules%HCHO_OH_2_CO(indT)     * vSp(iHCHO_1d) + &
             & rules%HCHO_NO3_2_HNO3 * vSp(iNO3rad_1d) * vSp(iHCHO_1d) + &
             & rules%SO2_OH_2_SO3(indT,indZ)* vSp(iSO2_1d)

      aL_HO2 = rules%HO2_NO_2_NO2(indT)         * vSp(iNO_1d) + &
             & rules%HO2_NO2_2_HO2NO2(indT,indZ)* vSp(iNO2_1d) + &
             & rules%HO2_NO3_2_HNO3             * vSp(iNO3rad_1d) + &
             & rules%CH3O2_HO2_2_CH3OOH(indT)   * vSp(iCH3O2_1d) !+ &
!             & rules%HO2_PANrad_2_CH3O2(indT)   * vSp(iPANrad_1d)

      fTmp = bP_OH *bP_HO2 / aL_OH - aL_HO2

      !
      ! Attention!
      ! HO2 production terms can be zero during night if some components disappear.
      ! Thus, HO2 will go to zero by itself. Corresponding reactions for CH3O2 etc have to be
      ! left for day-time only!!
      !

      detBB = fTmp*fTmp
      det4AC = 4. * (aP_HO2 + aP_OH*bP_HO2/aL_OH) * fK_HO2_HO2 ! -4AC, indeed
      if (detBB >= 1e4*det4AC) then ! Linear equation
           !Should be here also for detBB == det4AC == 0
           vSpSL(iHO2_1d_sl) = - (aP_HO2 + aP_OH*bP_HO2/aL_OH)/fTmp
      else
           vSpSL(iHO2_1d_sl) = (fTmp + sqrt(detBB + det4AC)) / (2.*fK_HO2_HO2)
      endif
!           vSpSL(iHO2_1d_sl) = (fTmp + sqrt(fTmp * fTmp + 4. * (aP_HO2 + aP_OH*bP_HO2/aL_OH) * fK_HO2_HO2)) / &
!                            & (2.*fK_HO2_HO2)
      !  NEVER allow zero cHO2
      vSpSL(iHO2_1d_sl) = max(vSpSL(iHO2_1d_sl),  1e-13 * fMixing_ratio_2_cnc)

      !if (vSp(iHO2_1d) < 1e-13 * fMixing_ratio_2_cnc) vSp(iHO2_1d) = 1e-13 * fMixing_ratio_2_cnc  ! NEVER allow zero cHO2
!      if (vSpSL(iHO2_1d_sl) < 0.) then
!              ! NEVER allow zero cHO2
!              call msg("Fixing negative HO2, was", vSpSL(iHO2_1d_sl))
!              call msg("aL_OH aL_HO2",aL_OH, aL_HO2)
!              call msg("aL_OH - aL_HO2, fTmp",aL_OH - aL_HO2, fTmp)
!              call msg("BB, 4AC", detBB, det4AC)
!              call msg("Quadratic eq gives", &
!                        & (fTmp + sqrt(detBB + det4AC)) / (2.*fK_HO2_HO2))
!              call msg("Linear gives ", - (aP_HO2 + aP_OH*bP_HO2/aL_OH)/fTmp)
!              call msg("Forcing HO2 to ", 1e-13 * fMixing_ratio_2_cnc)
!              vSpSL(iHO2_1d_sl) =  1e-13 * fMixing_ratio_2_cnc
!      endif

      vSpSL(iOH_1d_sl) = (aP_OH + bP_OH * vSpSL(iHO2_1d_sl)) / aL_OH

    endif   ! ifNitrogen

if(ifReport)then
  call msg('Sun,J_HONO:',sun, fJ_HONO)
  call msg('SunD,OH:',sunD, cOH_forced)
  call report_vector('before')
endif

!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***
!
!    write(uFileOH,'(I4,2i2,a1,i2,a1,i2,2x,100(E8.2,1x))')fu_year(now),fu_mon(now),fu_day(now), &
!                                                     & ',',fu_hour(now),':',fu_min(now), &
!                                                     & cOH, cHO2, cNO3_rad, cOH_forced
!
!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***!!!***


    call check_nan('start after OH-HO2')

    !-------------------------------------------------------------------------------
    !
    ! Sulphur and ammonium sulphate
    ! SOx = SO2 + SO4 + H2SO4
    ! NHx = NH3 + NH4_S + NH4_N
    !
    if(vSp(iSO2_1d) < fLowCncThresh(iSO2_1d))then
      !
      ! No or little sulphur in cell. Put everything to garbage
      !
      garbage(iSO2_1d) = garbage(iSO2_1d) + vSp(iSO2_1d)
      vSp(iSO2_1d) = 0.
      vSpSL(iH2SO4_1d_sl) = 0.
    else
      !
      ! SO2 gas phase oxidation
      !
      rate_air = timestep_sec * rules%SO2_OH_2_SO3(indT,indZ) * cOH_forced
      !
      ! SO2 aquesous-phase oxidation - see n.11,pp24-26
      !
!!!!      if(metdat(ind_rel_hum) > 0.4 .and. metdat(ind_tempr) > 265.)then ! not too dry, neither too cold
!!!!        if(metdat(ind_rel_hum) > 0.95)then
!!!!!          scaling = 4.3267487  !9.0^(2/3)
!!!!!          scaling = 2.08  !9.0^(1/3)
!!!!          scaling = 2.42  !9.0^(0.3) * 1.25
!!!!        else
!!!!!          scaling = ((1.4 - metdat(ind_rel_hum)) / (1. - metdat(ind_rel_hum)))**0.666666666666667
!!!!!          scaling = ((1.4 - metdat(ind_rel_hum)) / (1. - metdat(ind_rel_hum)))**0.3333333333333333
!!!!          scaling = ((1.4 - metdat(ind_rel_hum)) / (1. - metdat(ind_rel_hum)))**0.3 * 1.25
!!!!        endif
!!!!        if(metdat(ind_tempr) < 275.)then  ! decrease the rate for freezing droplets
!!!!          scaling = scaling * (1.-(275-metdat(ind_tempr)) / 10.0)
!!!!        endif
!!!!        rate_water = timestep_sec * rules%conv_SO2_2_SO4_aqua_basic * scaling
!!!!      else
!!!!        rate_water = 0.0
!!!!      endif
!!!!      !
!!!!      ! This is for the aboo-o-ve-ABL cloud-cover dependent in-cloud processing
!!!!      ! Possibly, is to be added to the below-ABL mechanism.
!!!!      !
!!!!!      if(metdat(ind_hgt) > metdat(ind_ablh) .and. metdat(ind_tempr) > 255.)then ! above ABL, not too cold
!!!!!        scaling = 10. * metdat(ind_cloud_covr)  ! the more clouds the faster conversion
!!!!!        if(metdat(ind_tempr) < 275.)then  ! decrease the rate for freezing droplets
!!!!!          scaling = scaling * (1.-(275-metdat(ind_tempr)) / 20.0)
!!!!!        endif
!!!!!        rate_water = rate_water + timestep_sec * rules%conv_SO2_2_SO4_aqua_basic * scaling
!!!!!!      else
!!!!!!        rate_water = 0.0
!!!!!      endif



      ! Relative humidity, cloud cover and temperature dependent scaling 
      ! for the aqueous conversion rate
      scaling = 0.
   
      ! option 1: linear to rh above 0.5
      if(metdat(ind_rel_hum) > 0.5)then 
        scaling = metdat(ind_rel_hum) - 0.5
      endif
    
      !    ! option 2: water in the wet particle volume
      !    if(metdat(ind_rh) > 0.95)then
      !      scaling = 0.5
      !    elseif(metdat(ind_rh) > 0.4)then
      !      scaling = (((1.4 - metdat(ind_rh)) / (1. - metdat(ind_rh))) - 1.) * 0.5 / 8.
      !    endif
      !  
      ! in cloud (above abl)
      if(metdat(ind_hgt) > metdat(ind_ablh))then 
        scaling = scaling + metdat(ind_cloud_covr)
      endif
      
      ! Temperature dependence and freezing out
      if(metdat(ind_tempr) < 270.)then
        scaling = 0.
      else   
        scaling = scaling * (metdat(ind_tempr) - 270.)
        !    elseif(metdat(ind_tempr) < 275.)then  
        !      scaling = scaling * (1.-(275. - metdat(ind_tempr)) / 5.0)
        !    else
        !      scaling = scaling * (1.-(275. - metdat(ind_tempr)) / 5.0)
      endif
      
      ! Constant conversion rate
      !    scaling = 0.5
      
      !   ! Lotos-Euros parameterisation. Seems to be enormous + illogical, but for a try .. 
      !   ! factor 8.3 replaced with 0.5 * rules%conv_SO2_2_SO4_aqua_basic !
      !    scaling = 0.5 + metdat(ind_cloud_cvr)
      !    if(metdat(ind_rh) >= 0.9)then
      !      scaling = scaling * (1. + 10.*(metdat(ind_rh)-0.9))
      !    endif
      
      
      ! aqueous conversion rate
      rate_water = rules%conv_SO2_2_SO4_aqua_basic * scaling * timestep_sec
      
      !
      ! Take possible fast reduction of the SO2 mass
      !
      ftmp = rate_air + rate_water
      if (ftmp < 1e-3 ) then
        ftmp = 1. - 0.5*ftmp ! Accurate within single-precision
      else
        fTmp = (1.-exp(-fTmp)) / fTmp
      end if
      
      dSO2_air = vSp(iSO2_1d) * rate_air * fTmp
      dSO2_water = vSp(iSO2_1d) * rate_water * fTmp
      
      vSp(iSO2_1d) = vSp(iSO2_1d) - dSO2_air - dSO2_water
      vSpSL(iH2SO4_1d_sl) = dSO2_air
      vSpSL(iSO4w_1d_sl) = dSO2_water

    endif  ! If enough sulphur components


    !---------------------------------------------------------------------------------
    !
    ! Check ifNitrogen switch and compute the rest if needed
    !
    if(rules%ifNitrogen)then

      call check_nan('After sulphur')

      !
      ! Nitrogen chemistry is computed only in case of sufficient mass in cell
      !
      if(   (vSp(iNO_1d)     < fLowCncThresh(iNO2_1d)     ) .and. &
          & (vSp(iNO2_1d)    <  fLowCncThresh(iNO_1d)     ) .and. &
          & (vSp(iPAN_1d)    <  fLowCncThresh(iPAN_1d)   ) .and. &
          & (vSp(iHNO3_1d)   <  fLowCncThresh(iHNO3_1d) ) .and. &
          & (vSp(iNO3rad_1d) <  fLowCncThresh(iNO3rad_1d)    ) .and. &
          & (vSp(iCH3OOH_1d) <  fLowCncThresh(iCH3OOH_1d) ) .and. &
          & (vSp(iH2O2_1d)   <  fLowCncThresh(iH2O2_1d)   ) .and. &
          & (vSp(iCO_1d)     <  fLowCncThresh(iCO_1d)     ) .and. &
          & (vSp(iCH3O2_1d)  <  fLowCncThresh(iCH3O2_1d)  ) .and. &
          & (vSp(iHONO_1d)   <  fLowCncThresh(iHONO_1d)) )then
        !
        ! Small masses. Put everything to garbage
        !
        garbage(iNO_1d) = garbage(iNO_1d) + vSp(iNO_1d)
        garbage(iNO2_1d) = garbage(iNO2_1d) + vSp(iNO2_1d)
        garbage(iPAN_1d) = garbage(iPAN_1d) + vSp(iPAN_1d)
        garbage(iNO3rad_1d) = garbage(iNO3rad_1d) + vSp(iNO3rad_1d)
        garbage(iHNO3_1d) = garbage(iHNO3_1d) + vSp(iHNO3_1d)
        garbage(iH2O2_1d) = garbage(iH2O2_1d) + vSp(iH2O2_1d)
        garbage(iCO_1d) = garbage(iCO_1d) + vSp(iCO_1d)
        garbage(iCH3OOH_1d) = garbage(iCH3OOH_1d) + vSp(iCH3OOH_1d)
        garbage(iCH3O2_1d) = garbage(iCH3O2_1d) + vSp(iCH3O2_1d)
        garbage(iHONO_1d) = garbage(iHONO_1d) + vSp(iHONO_1d)

        vSp(iNO_1d) = 0.
        vSp(iNO2_1d) = 0.
        vSp(iPAN_1d) = 0.
        vSp(iNO3rad_1d) = 0.
        vSp(iHNO3_1d) = 0.
        vSp(iH2O2_1d) = 0.
        vSp(iCO_1d) = 0.
        vSp(iCH3OOH_1d) = 0.
        vSp(iCH3O2_1d) = 0.
        vSp(iHONO_1d) = 0.

      else  ! non-negligible mass. Make-up chemistry
        !
        ! Nitrogen chemistry depens on sun. 
        !
!call report_vector('Start nitrogen')

        if(vSp(iNO3rad_1d) > 0.)then  ! decompose it
          !
          ! Seemingly the fastest set of reactions is degradation of NO3 radical to NO and NO2.
          ! Notebook 9, pg 5-7.
          !
!call msg(fu_connect_strings(fu_time_to_io_string(now),', ix,real(iy)'),ix,real(iy))
!call msg('iLev, sunD',iLev,sunD)

          if(sunD > fSunMin) then
            !
            ! If NO is small enough and NO3 radical is in abundance, assume NO in equilibrium
            ! In another extreme, NO changes just a bit
            !
            if (fJ_NO3_2_NO2 + 1.5 * fJ_NO3_2_NO + 0.5 * rules%NO_NO3_2_NO2(indT) * vSp(iNO_1d) < 1e-20) then
               !This thing causes nans when all of fJ_NO3_2_NO2, fJ_NO3_2_NO and vSp(iNO_1d) == 0

               if (barkcount<101) then
                 if (barkcount<100) then
                    call msg_warning("Implementation gotcha Daytime NO3_rad chemistry in transform_acid_basic")
                 else
                   call msg("Gotcha Daytime NO3_rad chemistry in transform_acid_basic. Giving up....")
                 endif
               endif
               !$OMP ATOMIC
               barkcount = barkcount + 1 

               fEqRatio =  1e6 !Whatver
            else 



              fEqRatio = (rules%NO_NO3_2_NO2(indT) * vSp(iNO3rad_1d)) / &
                     & (fJ_NO3_2_NO2 + 1.5 * fJ_NO3_2_NO + 0.5 * rules%NO_NO3_2_NO2(indT) * vSp(iNO_1d))
            endif


            fTmp = fJ_NO3_2_NO / rules%NO_NO3_2_NO2(indT) ! equilibrium point for NO

            if(fEqRatio > 2.0)then
              !
              ! Enough NO3: equilibrium for NO, zero for NO3 radical, NO2 takes the rest
              !
              fTmp = fJ_NO3_2_NO / rules%NO_NO3_2_NO2(indT) ! equilibrium point for NO
            else
              !
              ! Lack of NO3, equilibrium is approached but not entirely. Update fTmp - it is now
              ! the equilibrium point plus shift.
              !
              fTmp = fTmp + (vSp(iNO_1d) - fTmp) / &
                          & (1. + fEqRatio * (1. + fEqRatio * (0.5 + fEqRatio * 0.33)))
            endif ! if NO equil
            dNO = fTmp - vSp(iNO_1d)

            vSp(iNO_1d) = max(0.,fTmp)
            vSp(iNO2_1d) = max(0.,vSp(iNO2_1d) + vSp(iNO3rad_1d) - dNO)
            vSp(iNO3rad_1d) = 0.

            call check_nan('Daytime NO3_rad destroyed')

          else
            !
            ! At night, the fast reaction annihilates NO and NO3_rad: NO + NO3_rad -> 2*NO2
            ! Daytime, the same reaction competes with sun, thus a ratio was taken
            !
            if(vSp(iNO_1d) > vSp(iNO3rad_1d))then
              vSp(iNO_1d) =  vSp(iNO_1d) - vSp(iNO3rad_1d)
              vSp(iNO2_1d) = vSp(iNO2_1d) + 2. * vSp(iNO3rad_1d)
              vSp(iNO3rad_1d) = 0.
            else
              vSp(iNO3rad_1d) = vSp(iNO3rad_1d) - vSp(iNO_1d)
              vSp(iNO2_1d) = vSp(iNO2_1d) + 2. * vSp(iNO_1d)
              vSp(iNO_1d) = 0.
            endif

          endif   ! day/night time for NO3 radical

        endif  ! if NO3 radical present

        !
        ! If there is some HONO, it has comparatively fast photolysis making a blow to OH and
        ! NO in the morning. Still, due to sun the lifetime can be longer than the time step.
        ! All-in-all, three reactions at present: 
        ! NO+OH->HONO
        ! HONO+OH->NO2
        ! HONO+sun->NO+OH
        ! Several options depending on the insolation: 
        ! - fast decay, equilibrium
        ! - moderate decay, equilibrium with relaxation to it
        ! - slow decay, incremental
        !
        call check_mass_vector(vSp, garbage, species_list, 'Before HONO', fLowCncThresh, nSpecies, ix, iy, iz, print_it)

        if(vSp(iHONO_1d) + vSp(iNO_1d) > fLowCncThresh(iHONO_1d) + fLowCncThresh(iNO_1d))then
          if (fJ_HONO*timestep_sec > tiny(fJ_HONO) .or. cOH_forced > 1e-15) then
            ! below fails if fJ_HONO and cOH_forced are both zero (in which case none of
            ! the reactions happen). 
            dNO = vSp(iNO_1d) - (vSp(iNO_1d) + vSp(iHONO_1d)) * &
                 & fJ_HONO / (fJ_HONO + rules%NO_OH_2_HONO(indT,indZ) * cOH_forced)
            
!          if(fJ_HONO * timestep_sec < 10.)then
!            !
!            ! Moderate or slow photolysis rate. Equilibrium will not be reached within one time step
!            ! Also, a question who is faster arises
!            !
            if(fJ_HONO > rules%NO_OH_2_HONO(indT,indZ) * cOH_forced)then
              dNO = dNO * (1. - exp(-fJ_HONO * timestep_sec))
            else
              dNO = dNO * (1. - exp(-rules%NO_OH_2_HONO(indT,indZ) * cOH_forced * timestep_sec))
            endif
            !          endif
            vSp(iNO_1d) = vSp(iNO_1d) - dNO
            vSp(iHONO_1d) = vSp(iHONO_1d) + dNO
            !
            ! The final kick: HONO + OH -> NO2. Since it is slower than the above reactions
            ! can be taken after them.
            !
            dNO2 = vSp(iHONO_1d) * (1. - exp(-rules%HONO_OH_2_NO2(indT) * cOH_forced * timestep_sec))
            vSp(iHONO_1d) = vSp(iHONO_1d) - dNO2
            vSp(iNO2_1d) = vSp(iNO2_1d) + dNO2
          end if
        else
          vSp(iNO_1d) = vSp(iNO_1d) + vSp(iHONO_1d)
          vSp(iHONO_1d) = 0.
        endif

        call check_mass_vector(vSp, garbage, species_list, 'After HONO', fLowCncThresh, nSpecies, ix, iy, iz, print_it)


        ! ATTENTION. Below TLAM3=rules%NO_2_NO2_direct is a fake coefficient made in order to
        !            exclude VOC. It then corresponds to the assumption that everywhere the
        !            ozone situation is NOx-limited. In original DMAT scheme it is so.
        !            Whenever possible it should be replaced with NO_2_NO2_by_VOC thus including
        !            the ozone-creation potential approach
        !

        !---------------------------------------------------------------------------
        !
        ! The main photochemical lab NO-NO2-O3 during day:
        ! Equilibrium of  NO +O3 =NO2 in presence of light. Here we apply the standard photostationary
        ! state relation with dynamic relaxation to it.
        !
        if(sunD > fSunMin)then
        
          if(vSp(iNOP_1d) < 1.e3 * (vSp(iNO_1d) + vSp(iNO2_1d)) .and. &
           & (vSp(iNO_1d) + vSp(iNO2_1d) +  vSp(iNOP_1d)) > &
            & fLowCncThresh(iNO_1d) + fLowCncThresh(iNO2_1d) +  fLowCncThresh(iNOP_1d))then
            !
            ! Sufficient NOx mass to make equilibrium
            !
            ! Compute the photostationary equilibrium values from initial ones (notebook 8, p.73)
            ! This is just a squared algebraic equation. However, componets can be very different,
            ! so have to be careful cutting out wrong stuff
            !
!fTmp = vSp(iNOP_1d) - vSp(iNO_1d) - (fJ_no2 / rules%NO_2_NO2_with_O3(indT))
!fO3_eq = (fTmp + sqrt(fTmp*fTmp + 4 * fJ_no2 / rules%NO_2_NO2_with_O3(indT) * &
!                                                  & (vSp(iNOP_1d) + vSp(iNO2_1d)))) * 0.5
!call msg('Compare the eq O3 (short, long):', &
!       & (vSp(iNOP_1d) + vSp(iNO2_1d)), &
!       & fO3_eq)


            if(vSp(iNOP_1d) + vSp(iNO2_1d) + vSp(iNO_1d) < 0.01 * fJ_no2 / rules%NO_2_NO2_with_O3(indT))then
              fO3_eq = vSp(iNOP_1d) + vSp(iNO2_1d)
              fNO_eq = vSp(iNO_1d) + vSp(iNO2_1d)
              fNO2_eq = 0.0
            else
              fTmp = vSp(iNOP_1d) - vSp(iNO_1d) - (fJ_no2 / rules%NO_2_NO2_with_O3(indT))
              fO3_eq = max((fTmp + sqrt(fTmp*fTmp + 4 * fJ_no2 / rules%NO_2_NO2_with_O3(indT) * &
                                                      & (vSp(iNOP_1d) + vSp(iNO2_1d)))) * 0.5, 0.)
              fNO_eq = max(vSp(iNO_1d) + fO3_eq - vSp(iNOP_1d), 0.)
              fNO2_eq = vSp(iNO2_1d) + vSp(iNO_1d) - fNO_eq ! cannot put max(.,0): keep the budget!
              if(fNO2_eq < 0.0)then   ! full consumption
                fNO2_eq = 0.0
                fNO_eq = vSp(iNO2_1d) + vSp(iNO_1d)
              endif
            endif

            !
            ! Having the equilibrium values, let's take the dynamic relaxation in order to avoid too fast
            ! forcing of this equilibrium. The relaxation rate depends on all components and thus dynamic.
            ! However, we hope that our time step is small enough to allow for linearity assumption.
            !
            fRate = max(rules%NO_2_NO2_with_O3(indT) * max((vSp(iNO_1d) + fNO_eq), &
                                                         & (vSp(iNOP_1d) + fO3_eq)) * 0.5, &
                      & fJ_NO2) * timestep_sec
            if(fRate  < 2)then
              !
              ! Slow relaxation, be careful. Compute the dynamic approximation to the equilibrium
              !
              fO3_eq = fO3_eq + (vSp(iNOP_1d)-fO3_eq) / (1.+fRate*(1.+fRate*(0.5+fRate/3.)))
              fNO_eq = max(vSp(iNO_1d) + fO3_eq - vSp(iNOP_1d), 0.)
              fNO2_eq = max(vSp(iNO_1d) + vSp(iNO2_1d) - fNO_eq, 0.)
            endif

if(fO3_eq < 0. .or. fNO_eq < 0. .or. fNO2_eq < 0.)then
call report_vector('photost.eq')
call msg('',fRate)
call msg('fJ_no2, rules%NO_2_NO2_with_O3(indT)',fJ_no2,  rules%NO_2_NO2_with_O3(indT))
call msg('NO and NO2 eq: ', fNO_eq, fNO2_eq)
call msg('fJ_no2/rules%NO_2_NO2_with_O3(indT) O3 eq: ', fJ_no2/rules%NO_2_NO2_with_O3(indT), fO3_eq)
endif

            vSp(iNOP_1d) = fO3_eq
            vSp(iNO_1d) = fNO_eq
            vSp(iNO2_1d) = fNO2_eq

          endif  ! if mass is enough for equilibrium

call check_mass_vector(vSp, garbage, species_list, 'After photostationary eq', &
                               & fLowCncThresh, nSpecies, ix, iy, iz, print_it)

          !
          ! Ozone photo-destruction is quite slow, still take implicit scheme
          ! Its destruction is the key production of OH: do it!
          !
          dO3 = (fJ_O3_O + fJ_O3_O1d) * timestep_sec
!          dO3 = dO3 * (1 + dO3 * (0.5 + 0.333 * dO3)) !exp(dO3) - 1
!          dO3 = vSp(iNOP_1d) * dO3 / (1. + dO3)   !* 
!          vSp(iNOP_1d) = vSp(iNOP_1d) - dO3
          vSp(iNOP_1d) = vSp(iNOP_1d) * exp(-dO3)

call check_mass_vector(vSp, garbage, species_list, 'After ozone photolysis', &
                               & fLowCncThresh, nSpecies, ix, iy, iz, print_it)

!        !
!        ! PAN family halnding consists of 3 parts: production, internal equilibrium and loss
!        ! Production comes from formaldehyde by OH, loss via NO, equilibrium with NO2
!        ! Since formaldehyde is not yet available explicitly, have to fake its vertical profile to
!        ! avoid over-producing PAN in the upper troposphere where there is no NO to kill it.
!        !
!        ! See notebook 8 p.79
!        !
!        if(sunD > 0.)then
!          fPANprod = rules%PAN_radical_prod(indT) * (1.+100.*sunD) ! OH
!        else
!          fPANprod = rules%PAN_radical_prod(indT) ! OH
!        endif
!
!        if(vSp(iNO2_1d) < 1.e-20)then 
!          !
!          ! if no NO2 - assume all new stuff destroyed by NO. Actually, it may be different in case of 
!          ! strong losses but so far this is a minor headache.
!          !
!          dNO = min(rules%PAN_radical_prod(indT) * exp(-z*5.) * timestep_sec, &
!                  & vSp(iNO_1d))
!          vSp(iNO_1d) = vSp(iNO_1d) - dNO
!
!        else  ! some NO2 exists
!
!          dPAN = (rules%PAN_radical_prod(indT) * exp(-z*5.) * &
!                & rules%PAN_radical_2_PAN_eq(indT, indZ) * vSp(iNO2_1d) - &
!                & rules%PAN_radical_NO_loss * vSp(iNO_1d) * vSp(iPAN_1d)) / &
!               & (1. + vSp(iPAN_1d) / vSp(iNO2_1d))
!
!          if(dPAN > 0)then  ! Sucks NO2
!            nSteps = max(1., 2.*dPAN*timestep_sec / vSp(iNO2_1d))
!          else              ! Produces NO2, destroys PAN
!            nSteps = max(1., -2.*dPAN*timestep_sec / vSp(iPAN_1d))
!          endif
!
!          dNO = rules%PAN_radical_NO_loss * vSp(iNO_1d) * vSp(iPAN_1d) / &
!              & (1. + vSp(iPAN_1d) / vSp(iNO2_1d))
!          nSteps = max(nSteps, int(2.*dNO*timestep_sec / vSp(iNO_1d)))
!
!          do iStep = 1, nSteps
!            dPAN = (rules%PAN_radical_prod((indT)) * exp(-z*5.) * &
!                  & rules%PAN_radical_2_PAN_eq(indT, indZ) * vSp(iNO2_1d) - &
!                  & rules%PAN_radical_NO_loss * vSp(iNO_1d) * vSp(iPAN_1d)) / &
!                 & (1. + vSp(iPAN_1d) / vSp(iNO2_1d))  * timestep_sec / nSteps
!            dNO = rules%PAN_radical_NO_loss * vSp(iNO_1d) * vSp(iPAN_1d) / &
!                & (1. + vSp(iPAN_1d) / vSp(iNO2_1d)) * timestep_sec /nSteps
!            vSp(iPAN_1d) = vSp(iPAN_1d) + dPAN
!            vSp(iNO2_1d) = vSp(iNO2_1d) - dPAN
!            vSp(iNO_1d) = vSp(iNO_1d) - dNO
!          end do
!
!        endif  ! if there is some NO2 to store in PAN
!
!
!        if(vSp(iPAN_1d) < 0.)then
!          call msg('PAN<0')
!        endif
!        if(vSp(iNO2_1d) < 0.)then
!          call msg('after PAN: NO2<0')
!        endif
!        if(vSp(iNO_1d) < 0.)then
!          call msg('after PAN: NO2<0')
!        endif

          !
          ! So far restore the original DMAT PAN equilibrium
          !
          ! Equilibriums: NO2 <--> PAN. Since the representation is crude, no reason to 
          ! go after funcy exponents or so. Just ensure positive masses. Note that for very hot air
          ! and long time step (e.g. 20 min and T>320K) more than 100% of PAN is destroyed
          !
          dNO2_2_PAN = min(vSp(iNO2_1d), max(-vSp(iPAN_1d), &
                     & (rules%NO2_2_PAN_rate * vSp(iNO2_1d) * sunD - &
                      & vSp(iPAN_1d) * rules%PAN_2_NO2_sun(indT)) * timestep_sec))
          vSp(iPAN_1d) = vSp(iPAN_1d) + dNO2_2_PAN
          vSp(iNO2_1d) = vSp(iNO2_1d) - dNO2_2_PAN
call check_mass_vector(vSp, garbage, species_list, 'After PAN', &
                               & fLowCncThresh, nSpecies, ix, iy, iz, print_it)

          !
          ! A slow photolysis of HNO3, no need in copmplex implicit form
          !
          dNO2 = fJ_HNO3 * timestep_sec
          dNO2 = vSp(iHNO3_1d) * dNO2 / (1. + dNO2)
          vSp(iNO2_1d) = vSp(iNO2_1d) + dNO2
          vSp(iHNO3_1d) = vSp(iHNO3_1d) - dNO2

          !
          ! CHECK MASSES
          !
          call check_nan('End of primary photo-lab')
call check_mass_vector(vSp, garbage, species_list, 'End of primary photo-lab', &
                               & fLowCncThresh, nSpecies, ix, iy, iz, print_it)
          if(error)return

        else
          !---------------------------------------------------------------------------
          !---------------------------------------------------------------------------
          !---------------------------------------------------------------------------
          !
          ! Darkness
          !
call check_mass_vector(vSp, garbage, species_list, 'Start of dark chem', &
                               & fLowCncThresh, nSpecies, ix, iy, iz, print_it)

! ATTENTION. Remnants of DMAT code for ozone and NO-NO2. No justification exist. Commented so far
!          !
!          ! Ozone self-decay exists but it is no longer pushed by sun
!          !
!          vSp(iNOP_1d) = vSp(iNOP_1d) / (1. + rules%O3_decay_dark * timestep_sec)
!
!          !
!          ! Oxidation NO->NO2 by VOCs continues at low pace: stuff like CO still remains important
!          !
!          dNO = rules%NO_2_NO2_direct(indT) * fOH_night_vs_day * timestep_sec
!          dNO = dNO * (1 + dNO * (0.5 + 0.333 * dNO))
!          dNO = dNO / (1. + dNO) * vSp(iNO_1d)
!          vSp(iNO2_1d) = vSp(iNO2_1d) + dNO
!          vSp(iNO_1d) = vSp(iNO_1d) - dNO

          !
          ! NO goes to NO2 eating out ozone.
          !
          if(vSp(iNO_1d) < vSp(iNOP_1d))then
            dNO = vSp(iNO_1d) * (1. - exp(-rules%NO_2_NO2_with_O3(indT) * vSp(iNOP_1d) * timestep_sec))
          else
            dNO = vSp(iNOP_1d) * (1. - exp(-rules%NO_2_NO2_with_O3(indT) * vSp(iNO_1d) * timestep_sec))
          endif

          vSp(iNO_1d) = vSp(iNO_1d) - dNO
          vSp(iNOP_1d) = vSp(iNOP_1d) - dNO
          vSp(iNO2_1d) = vSp(iNO2_1d) + dNO
          !
          ! convertion NO2+O3 --> NO3_rad in darkness. Can be fast! Can eat out all O3 or all NO2.
          ! Have to be very careful. The agent of smaller quantities will be the fast-moving
          ! one with regard to which we solve the problem via implicit scheme. The other one
          ! will be considered slowly changing.
          !
!          DNO2_2_no3 = vSp(iNO2_1d) * vSp(iNOP_1d) * rules%NO2_2_NO3_dark(indT) * timestep_sec

          if(vSp(iNO2_1d) > vSp(iNOP_1d))then
            dNO2 = vSp(iNO2_1d) * rules%NO2_2_NO3_dark(indT) * timestep_sec
!            dNO2 = dNO2 * (1. + dNO2 * (0.5 + dNO2 * 0.333))  !!exp(x)-1
!            dNO2 = vSp(iNOP_1d) * dNO2 / (1. + dNO2)  ! (exp(x)-1) / exp(x) = 1 - exp(-x)
             dNO2 = vSp(iNOP_1d) *(1. - exp(-dNO2))
          else
            dNO2 = vSp(iNOP_1d) * rules%NO2_2_NO3_dark(indT) * timestep_sec
!            dNO2 = dNO2 * (1. + dNO2 * (0.5 + dNO2 * 0.333))
!            dNO2 = vSp(iNO2_1d) * dNO2 / (1. + dNO2)
            dNO2 = vSp(iNO2_1d) * (1. - exp(-dNO2))
          endif

          vSp(iNO2_1d) = vSp(iNO2_1d) - dNO2
          vSp(iNOP_1d) = vSp(iNOP_1d) - dNO2
          vSp(iNO3rad_1d) = vSp(iNO3rad_1d) + dNO2

          if(.not.(vSp(iNOP_1d) >= 0.0 .and. vSp(iNO2_1d)>=0) )then
            call msg('Negative O3 in dark after comp,NO2=', vSp(iNO2_1d))
            call msg('Negative O3 in dark after comp,O3=', vSp(iNOP_1d))
            call msg('Negative O3 in dark after comp,indT, NO2_2_NO3_', indT, rules%NO2_2_NO3_dark(indT))
            call set_error("Gotcha  NO2+O3 --> NO3_rad!", "here")
          endif

          call check_nan('Darkness: before N2O5')

          if(vSp(iNO3rad_1d)  > fLowCncThresh(iNO3rad_1d) .and. &
           & vSp(iNO2_1d) > fLowCncThresh(iNO2_1d))then
            !
            ! Conversion of NO3_rad + NO2 <-> N2O5 -> HNO3
            ! The first step is equilibrium for getting N2O5, which is then hydrolised to HNO3. 
            ! The N2O5 production takes one N molecule from each. Having the equilibrium constant, 
            ! one can get the square equation to be solved. See notebook 9 p.12-14
            ! Here we assume N2O5 start = 0
            !
            fTmp = vSp(iNO3rad_1d) + vSp(iNO2_1d) + 1. / rules%NO3_NO2_N2O5_eq(indT)
            if(isNaN(fTmp))then
              call msg('fTmp is NaN: fTmp',fTmp)
              call msg('fTmp is NaN: NO3_rad',vSp(iNO3rad_1d))
              call msg('fTmp is NaN: NO2',vSp(iNO2_1d))
              call msg('fTmp is NaN: iTempr,eq_coef',indT,rules%NO3_NO2_N2O5_eq(indT))
            endif
            !
            ! Have to take care about numerics:
            !
            if(fTmp * fTmp > 400. * vSp(iNO3rad_1d) * vSp(iNO2_1d))then
              cN2O5 = vSp(iNO3rad_1d) * vSp(iNO2_1d) / fTmp
            else
              cN2O5 = 0.5 * (fTmp - sqrt(max(0., fTmp * fTmp - 4. * vSp(iNO3rad_1d) * vSp(iNO2_1d))))
            endif
            if(cN2O5 > vSp(iNO2_1d))   cN2O5 = vSp(iNO2_1d)
            if(cN2O5 > vSp(iNO3rad_1d))cN2O5 = vSp(iNO3rad_1d)

            if(.not. cN2O5 >= 0)then
              call msg('N2O5 < 0, N2O5',cN2O5)
              call msg('N2O5 < 0, NO2',vSp(iNO2_1d))
              call msg('N2O5 < 0, NO3_rad',vSp(iNO3rad_1d))
              call msg('N2O5 < 0, ,tempr',tempr)
              call msg('N2O5 < 0, eq_coef', rules%NO3_NO2_N2O5_eq(indT))
            endif

            ! ... and hydrolysis. Here fH2O comes from specific humidity, converted to mole/m3
            ! Implicit scheme, of course, although the process seems to be slow.
            !
            dNO2 = rules%N2O5_2_HNO3 * cH2O * timestep_sec  
!            dNO2 = dNO2 * (1 + dNO2 * (0.5 + 0.333 * dNO2)) !!exp(x)-1
!            dNO2 = cN2O5 * dNO2 / (1. + dNO2)  ! (exp(x)-1) / exp(x) = 1 - exp(-x) 
            dNO2 = cN2O5 * (1. - exp(-dNO2)) 

            vSp(iNO2_1d) = vSp(iNO2_1d) - dNO2
            vSp(iNO3rad_1d) = vSp(iNO3rad_1d) - dNO2
            vSp(iHNO3_1d) = vSp(iHNO3_1d) + 2. * dNO2

            if(.not. all( (/vSp(iNO3rad_1d), vSp(iNO2_1d), vSp(iHNO3_1d)/) >= 0.0 ))then
              !call msg('Negative NO3 after NO2+NO3->HNO3, coord:',ix,real(iy))
              call msg('Negative NO3 after NO2+NO3->HNO3, height:',fHeight)
              call msg('Negative NO3 after NO2+NO3->HNO3, N2O5:',  cN2O5)
              call msg('Negative NO3 after NO2+NO3->HNO3, NO3_rad:',  vSp(iNO3rad_1d))
              call msg('Negative NO3 after NO2+NO3->HNO3, dNO2:',  dNO2)
              call msg('Negative NO3 after NO2+NO3->HNO3, water:',  cH2O)
            endif
          endif

          !
          ! Slow reaction to push-up HNO3 from HCHO and NO3 radical
          !
          dNO2 = rules%HCHO_NO3_2_HNO3 * vSp(iHCHO_1d) * vSp(iNO3rad_1d) * timestep_sec
          if(vSp(iHCHO_1d) < dNO2)dNO2 = vSp(iHCHO_1d)
          if(vSp(iNO3rad_1d) < dNO2)dNO2 = vSp(iNO3rad_1d)
          vSp(iHCHO_1d) = vSp(iHCHO_1d) - dNO2
          vSp(iNO3rad_1d) = vSp(iNO3rad_1d) - dNO2
          vSp(iHNO3_1d) = vSp(iHNO3_1d) + dNO2
          vSp(iCO_1d) = vSp(iCO_1d) + dNO2

          call check_mass_vector(vSp, garbage, species_list, 'Before HONO night', &
                               & fLowCncThresh, nSpecies, ix, iy, iz, print_it)
          if(error) call report_vector('Before HONO night')
!          !
!          ! A comparatively slow reaction creating HONO from NO and OH (noticeable also daytime - n.9.p.39))
!          ! and even slower one converting HONO to NO2
!          !
!          dNO = vSp(iNO_1d) * (1. - exp(-rules%NO_OH_2_HONO(indT,indZ) * cOH_forced * timestep_sec))
!          vSp(iHONO_1d) = vSp(iHONO_1d) + dNO
!          vSp(iNO_1d) = vSp(iNO_1d) - dNO
!
!          dNO2 = vSp(iHONO_1d) * (1. - exp(-rules%HONO_OH_2_NO2(indT) * cOH_forced * timestep_sec))
!          vSp(iHONO_1d) = vSp(iHONO_1d) - dNO2
!          vSp(iNO2_1d) = vSp(iNO2_1d) + dNO2
!
!          call check_mass_vector(vSp, garbage, species_list, 'After HONO, darkness', &
!                               & fCncLowThresh, nSpecies, ix, iy, iz, print_it)
!          if(error)then
!            call report_vector('After HONO night')
!            call msg('dNO, dNO2', dNO, dNO2)
!          endif

        end if  ! light or darkness

!        call report_vector('End split')
        call check_mass_vector(vSp, garbage, species_list, 'End split', &
                             & fLowCncThresh, nSpecies, ix, iy, iz, print_it)

if(ifReport)then
  call msg('Sun,J_HONO:',sun, fJ_HONO)
  call msg('SunD,OH:',sunD, cOH_forced)
  call report_vector('after')
endif

        !---------------------------------------------------------------------------
        !---------------------------------------------------------------------------
        !---------------------------------------------------------------------------
        !---------------------------------------------------------------------------
        !
        ! Having initial NO-NO2-O3 computed, let's make up the production of NO2 from NO
        ! without O3: by CO, CH4 and VOC. Note that here some of these guys are represented by 
        ! TLAM3 - direct oxidation of NO by VOC. Can be fast, use implicit scheme (notebook 8 p.46).
        !
        ! Methane cycle starts from CH4+OH -> CH3O2 -> HCHO -> CO
        ! with NO->NO2 happening in a few occasions (n.9 p.30 for the scheme)
        ! The CH3O2 is an intermediate creature, with lifetime dependent on concentrations 
        ! of oxidising agents
        !
!        call report_vector('Before CH3O2, 1')

        pCH3O2 = (rules%CH4_cnc_global(indZ) * rules%CH4_OH_2_CH3O2(indT) + &
                & rules%CH3OOH_OH_2_CH3O2(indT) * vSp(iCH3OOH_1d)) * cOH_forced
        !
        ! Start CH3O2 cycle from getting the lifetimes and ratios for the key players
        ! The life time of CH3O2 is determined from HO2 and NO 
        !
!        call report_vector('Before CH3O2, 2')
!        call msg('pCH3O2:', pCH3O2)

        fRate = rules%CH3O2_HO2_2_CH3OOH(indT) * vSpSL(iHO2_1d_sl) + &
              & rules%CH3O2_NO_2_HCHO_NO2(indT) * vSp(iNO_1d)
!        call msg('fRate:', fRate)

        !
        ! Current equilibrium is made only under linear assumption, so have to check that
        ! squared term is small.
        !
!        if(fRate > 0.)then
!          if((rules%CH3O2_CH3O2_2_HCHO(indT) + rules%CH3O2_CH3O2_2_CH3O(indT)) * &
!                          & vSp(iCH3O2_1d) / fRate > 0.5) &
!          & call msg('Squared term vs linear:',  &
!                         & (rules%CH3O2_CH3O2_2_HCHO(indT) + rules%CH3O2_CH3O2_2_CH3O(indT)) * &
!                          & vSp(iCH3O2_1d) / fRate )
!        endif

!        if((rules%CH3O2_CH3O2_2_HCHO(indT) + rules%CH3O2_CH3O2_2_CH3O(indT)) * &
!          & vSp(iCH3O2_1d) * vSp(iCH3O2_1d) < 0.1 * fRate .and. fRate * timestep_sec > 4)then
        if(fRate * timestep_sec > 4)then
          !
          ! Short lifetime for CH3O2. Equilibrium between charge-discharge
          !
          vSp(iCH3O2_1d) = pCH3O2 / fRate

          if(vSp(iCH3O2_1d) < 0) call report_vector('CH3O2 equilibrium, 1')
          call check_nan('CH3O2 equilibrium, 1')
          call check_mass_vector(vSp, garbage, species_list, 'CH3O2 equilibrium, 1', &
                               & fLowCncThresh, nSpecies, ix, iy, iz, print_it)

          !
          ! NO -> NO2 conversion. Careful: if CH3O2 is in abundance, NO may change on the way,
          ! then the above equilibrium is to be adjusted. Strictly speaking, there must be
          ! an iterative procedure but we limit here with the first iteration because the 
          ! NO (and/or HO2) are in abundance anyway, so the correction must be small.
          !
          dNO = rules%CH3O2_NO_2_HCHO_NO2(indT) * vSp(iCH3O2_1d) * timestep_sec
          if(dNO > 0.5)then
            dNO = dNO * (1 + dNO * (0.5 + 0.333 * dNO))
            dNO = dNO / (1. + dNO)
!            call msg('Adjusted CH3O2-NO dNO*1000, CH3O2:',int(dNO*1000.),vSp(iCH3O2_1d))
            dNO = vSp(iNO_1d) * dNO
            fRate = rules%CH3O2_HO2_2_CH3OOH(indT) * vSpSL(iHO2_1d_sl) + dNO
            vSp(iCH3O2_1d) = pCH3O2 / fRate
          else
            dNO = vSp(iNO_1d) * dNO
          endif

          if(vSp(iCH3O2_1d) < 0) call report_vector('CH3O2 equilibrium, 2')

          vSp(iNO_1d) = vSp(iNO_1d) - dNO    ! CH3O2 + NO -> HCHO + NO2
          vSp(iNO2_1d) = vSp(iNO2_1d) + dNO
          vSp(iHCHO_1d) = vSp(iHCHO_1d) + dNO

          call check_nan('CH3O2 equilibrium, 2')
          call check_mass_vector(vSp, garbage, species_list, 'CH3O2 equilibrium, 2', &
                               & fLowCncThresh, nSpecies, ix, iy, iz, print_it)

          !
          ! Here we make a first part of CH3OOH budget - due to exchange with CH3O2
          !
          vSp(iCH3OOH_1d) = vSp(iCH3OOH_1d) + timestep_sec * &
                          & (rules%CH3O2_HO2_2_CH3OOH(indT) * vSpSL(iHO2_1d_sl) * vSp(iCH3O2_1d) - &
                           & rules%CH3OOH_OH_2_CH3O2(indT) * cOH_forced  * vSp(iCH3OOH_1d))

          call check_nan('CH3O2 equilibrium, 3')
          call check_mass_vector(vSp, garbage, species_list, 'CH3O2 equilibrium, 3', &
                               & fLowCncThresh, nSpecies, ix, iy, iz, print_it)

        else
          !
          ! Oxidants are in lack, long life time for CH3O2. The temporary reservoir CH3O2NO2 
          ! starts to play a role - being itself always in equilibrium with CH3O2 + NO2 (n.10,p.6)
          ! Split the time step to small sub-steps so that explicit scheme is OK.
          !
          fEqRatio = rules%CH3O2_NO2_2_CH3O2NO2(indT,indZ) * vSp(iNO2_1d) / &
                   & rules%CH3O2NO2_2_CH3O2_NO2(indT,indZ)

          nSteps = int(fRate * timestep_sec * 20.) + 1 ! <20, otherwise above equilibrium is forced
          fTmp = timestep_sec / real(nSteps)

          if(isNAN(fEqRatio) .or. isNAN(fTmp))then
           call msg('fEqRatio:',fEqRatio)
           call msg('nSteps:',real(nSteps))
           call msg('fTmp:',fTmp)
           call msg('indT,CH3O2NO2_2_CH3O2_NO2:',indT,rules%CH3O2NO2_2_CH3O2_NO2(indT,indZ))
           call msg('indZ,CH3O2_NO2_2_CH3O2NO2:',indZ,rules%CH3O2_NO2_2_CH3O2NO2(indT,indZ))
           call msg('iNO2, iNO2:',iNO2_1d,vSp(iNO2_1d))
          endif
          call check_nan('CH3O2 equilibrium, 4')

          do iStep = 1, nSteps
            dNO = rules%CH3O2_NO_2_HCHO_NO2(indT) * vSp(iCH3O2_1d) * fTmp
            dNO = dNO * (1 + dNO * (0.5 + 0.333 * dNO))
            dNO = vSp(iNO_1d) * dNO / (1. + dNO)

            vSp(iNO_1d) = vSp(iNO_1d) - dNO
            vSp(iNO2_1d) = vSp(iNO2_1d) + dNO
            vSp(iHCHO_1d) = vSp(iHCHO_1d) + dNO
            vSp(iCH3O2_1d) = vSp(iCH3O2_1d) - dNO



            call check_nan('CH3O2 equilibrium, 5')

            dCH3x = rules%CH3OOH_OH_2_CH3O2(indT) * cOH_forced  * vSp(iCH3OOH_1d) * fTmp
            vSp(iCH3OOH_1d) = vSp(iCH3OOH_1d) - dCH3x + &
                   & fTmp * rules%CH3O2_HO2_2_CH3OOH(indT) * vSpSL(iHO2_1d_sl) * vSp(iCH3O2_1d)

            vSp(iCH3O2_1d) = vSp(iCH3O2_1d) + (dCH3x + &
                 & rules%CH4_cnc_global(indZ) * rules%CH4_OH_2_CH3O2(indT) * cOH_forced * fTmp) / &
                 & (1. + fEqRatio)

            call check_nan('CH3O2 equilibrium, 6')

!            vSp(iCH3O2_1d) = vSp(iCH3O2_1d) - dNO + fTmp * &
!                         & (  (rules%CH4_cnc_global(indZ) * rules%CH4_OH_2_CH3O2(indT) + &
!                             & rules%CH3OOH_OH_2_CH3O2(indT) * vSp(iCH3OOH_1d)) * cOH_forced / &
!                            & (1. + fEqRatio) - &                          ! part goes to CH3O2NO2
!                            & (rules%CH3O2_HO2_2_CH3OOH(indT) * vSpSL(iHO2_1d_sl) + &
!                             & (rules%CH3O2_CH3O2_2_HCHO(indT) + rules%CH3O2_CH3O2_2_CH3O(indT)) * &
!                             & vSp(iCH3O2_1d)) * &
!                            & vSp(iCH3O2_1d) &
!                         & )
            dCH3x = fTmp * (rules%CH3O2_HO2_2_CH3OOH(indT) * vSpSL(iHO2_1d_sl) + &
                & (rules%CH3O2_CH3O2_2_HCHO(indT) + rules%CH3O2_CH3O2_2_CH3O(indT)) * vSp(iCH3O2_1d))
            dCH3x = dCH3x * (1 + dCH3x * (0.5 + 0.333 * dCH3x))
            dCH3x = vSp(iCH3O2_1d) * dCH3x / (1. + dCH3x)
            vSp(iCH3O2_1d) = vSp(iCH3O2_1d) - dCH3x

            if(vSp(iCH3O2_1d) < 0)then
              call report_vector('Inner CH3O2 cycle')
              call msg('iStep,fTmp',iStep,fTmp)
              call msg('indZ,CH4',indZ,rules%CH4_cnc_global(indZ))
              call msg('indT,rules%CH4_OH_2_CH3O2(indT)',indT,rules%CH4_OH_2_CH3O2(indT))
              call msg('rules%CH3OOH_OH_2_CH3O2(indT)', rules%CH3OOH_OH_2_CH3O2(indT))
              call msg('rules%CH3O2_HO2_2_CH3OOH(indT)', rules%CH3O2_HO2_2_CH3OOH(indT))
              call msg('rules%CH3O2_CH3O2_2_HCHO(indT)', rules%CH3O2_CH3O2_2_HCHO(indT))
              call msg('rules%CH3O2_CH3O2_2_CH3O(indT)', rules%CH3O2_CH3O2_2_CH3O(indT))
              call msg('fEqRatio', fEqRatio)
              call msg('dNO', dNO)
              call msg('cOH_forced', cOH_forced)
            endif

            call check_nan('CH3O2 equilibrium, 7')
            call check_mass_vector(vSp, garbage, species_list, 'CH3O2 equilibrium, 7', &
                                 & fLowCncThresh, nspecies, ix, iy, iz, print_it)

          end do  ! smaller time steps for CH3O2

        endif  ! if equilibrium for CH3O2

!        call report_vector('After CH3O2')
        call check_nan('After CH3O2')

        if(vSp(iCH3O2_1d) < 0) call report_vector('after the CH3O2 equilibrium')

!        if(iHourTmp == 5 .and. iLev == 1)then
!          call msg('indT, timestep',indT, timestep_sec)
!!          call msg('indT, pCH3O2',indT, pCH3O2)
!          call msg('iCH3OOH, cCH3O2',iCH3OOH_1d,vSp(iCH3O2_1d))
!!          call msg('size of CH3OOH_OH_2_CH3O2',size(rules%CH3OOH_OH_2_CH3O2))
!          call msg('rules%CH3OOH_OH_2_CH3O2(indT)', rules%CH3OOH_OH_2_CH3O2(indT))
!          call msg('vSp(iCH3OOH_1d)', vSp(iCH3OOH_1d))
!          call msg('dCH3OOH', dCH3OOH)
!        endif

        !
        ! Remaining CH3OOH -> HCHO are to be made with implicit scheme: OH and sun
        !
        dCH3x = (rules%CH3OOH_OH_2_HCHO * cOH_forced + fJ_CH3OOH) * timestep_sec
        dCH3x = dCH3x * (1 + dCH3x * (0.5 + 0.333 * dCH3x))
        dCH3x = vSp(iCH3OOH_1d) * dCH3x / (1. + dCH3x)

        vSp(iCH3OOH_1d) = vSp(iCH3OOH_1d) - dCH3x
        vSp(iHCHO_1d) = vSp(iHCHO_1d) + dCH3x

        !
        ! Finally, HCHO -> CO: three slow reactions
        !
        dHCHO = (fJ_HCHO_2_HO2 + fJ_HCHO_2_H2 + rules%HCHO_OH_2_CO(indT) * cOH_forced) * &
              & timestep_sec
        dHCHO = dHCHO * (1 + dHCHO * (0.5 + 0.333 * dHCHO))
        dHCHO = vSp(iHCHO_1d) * dHCHO / (1. + dHCHO)
        vSp(iHCHO_1d) = vSp(iHCHO_1d) - dHCHO
        vSp(iCO_1d) = vSp(iCO_1d) + dHCHO

        call check_nan('After HCHO')

        !
        ! CO oxidation with OH
        !
        dCO = rules%CO_OH_2_CO2(indT,indZ) * cOH_forced * timestep_sec ! CO oxidation
        dCO = dCO * (1 + dCO * (0.5 + 0.333 * dCO))
        dCO = vSp(iCO_1d) * dCO / (1. + dCO)
        vSp(iCO_1d) = vSp(iCO_1d) - dCO

        !
        ! Flux from CO oxidation is split to H2O2 and NO2 production, depending on the 
        ! concentrations of HO2 and NO. Careful with NO: it may be in lack
        !
        if(vSp(iNO_1d) == 0.)then ! speedup
          dNO = 0.
        else 
          if(rules%HO2_NO_2_NO2(indT) * vSp(iNO_1d) / (fK_HO2_HO2 * vSpSL(iHO2_1d_sl)) > 1000.) then 
            dNO = dCO
          else
            dNO = dCO * (rules%HO2_NO_2_NO2(indT) * vSp(iNO_1d)) / &
                    & (rules%HO2_NO_2_NO2(indT) * vSp(iNO_1d) + fK_HO2_HO2 * vSpSL(iHO2_1d_sl))
          endif
        endif
        if(dNO > vSp(iNO_1d)) dNO = vSp(iNO_1d)
        vSp(iNO_1d) = vSp(iNO_1d) - dNO
        vSp(iNO2_1d) = vSp(iNO2_1d) + dNO
        vSp(iH2O2_1d) = vSp(iH2O2_1d) + (dCO - dNO)
         

        call check_nan('After CO')

        !
        ! A slower reaction of H2O2 destruction via OH and sun
        ! Interestingly, since we do not have OH and HO2 detailed budget, these fluxes are
        ! of no interest, apart from H2O2 desruction itself
        !
        dH2O2 = (fJ_H2O2 + rules%H2O2_OH_2_HO2(indT) * cOH_forced) * timestep_sec
        dH2O2 = dH2O2 * (1. + dH2O2 * (0.5 + 0.333 * dH2O2))
        dH2O2 = vSp(iH2O2_1d) * dH2O2 / (1. + dH2O2)
        vSp(iH2O2_1d) = vSp(iH2O2_1d) - dH2O2

        if (vSp(iH2O2_1d) < 0. ) then
          if (abs(vSp(iH2O2_1d)) < fLowCncThresh(iH2O2_1d))then
            garbage(iH2O2_1d) = garbage(iH2O2_1d) + vSp(iH2O2_1d)
            vSp(iH2O2_1d) = 0.
          else
            call msg ('Negative H2O2 in the beginning, H2O2',   vSp(iH2O2_1d))
            call msg ('Negative H2O2 in the beginning, dNO, dCO',   dNO, dCO)
            call msg ('Negative H2O2 in the beginning, NO2',   vSp(iNO2_1d)) 
            call msg ('Negative H2O2 in the beginning, NO',   vSp(iNO_1d))
            call msg ('Negative H2O2 in the beginning, CO',   vSp(iCO_1d))
          endif
        endif 

        call check_nan('After H2O2')

        !
        ! This is the remnants of direct VOC oxidation. All what is left is 
        ! the anthropogenic and non-methane biogenic VOCs. So, the coefficient can be 
        ! significantly reduced
        !
        dNO = rules%NO_2_NO2_direct(indT) * (sunD + fOH_night_vs_day) * timestep_sec   
        dNO = dNO * (1 + dNO * (0.5 + 0.333 * dNO))
        dNO = vSp(iNO_1d) * dNO / (1. + dNO)
        vSp(iNO2_1d) = vSp(iNO2_1d) + dNO
        vSp(iNO_1d) = vSp(iNO_1d) - dNO

        if(vSp(iNO_1d) < 0 .or. vSp(iNO2_1d) < 0.0)then
          call msg("indT,indZ", indT,indZ)
          call msg("sunD,  fOH_night_vs_day", sunD,  fOH_night_vs_day )
          call msg("rules%NO_2_NO2_direct(indT) timestep_sec", rules%NO_2_NO2_direct(indT), timestep_sec)
          call msg("dNO", dNO)
          call msg("Now NO2, NO:", vSp(iNO2_1d) , vSp(iNO_1d) )
          call msg("Initially (apprx) NO2, NO:", vSp(iNO2_1d) - dNO, vSp(iNO_1d) + dNO)
          call check_nan('Before PAN equilibrium')
          call check_mass_vector(vSp, garbage, species_list, 'Before PAN equilibrium', &
                               & fLowCncThresh, nspecies, ix, iy, iz, print_it)
          if(error)return
        endif

        !------------------------------------------------------------------------
        !
        ! convertion NO2 + OH -> HNO3. About 0.05-0.1 hr-1. Slow.
        ! In DMAT paper: 0.08 * sunD**2, in code: 0.05 * sunD
        !
        dHNO3 = rules%NO2_OH_2_HNO3(indT,indZ) * vSp(iNO2_1d) * cOH_forced * timestep_sec
        vSp(iNO2_1d) = vSp(iNO2_1d) - dHNO3
        vSp(iHNO3_1d) = vSp(iHNO3_1d) + dHNO3

        if(vSp(iHNO3_1d) < 0.0 .or. vSp(iNO2_1d) < 0.0)then
          call msg("indT,indZ", indT,indZ)
          call msg("cOH_forced,  timestep_sec", cOH_forced, timestep_sec)
          call msg("dHNO3", dHNO3)
          call msg("Now  NO2, HNO3:", vSp(iNO2_1d), vSp(iHNO3_1d) )
          call msg("Initially (apprx) NO2, HNO3:", vSp(iNO2_1d) + dHNO3, vSp(iHNO3_1d) - dHNO3)

          call check_nan('After NO2->HNO3 conversion')
          call check_mass_vector(vSp, garbage, species_list, 'After NO2->HNO3 conversion', &
                               & fLowCncThresh, nSpecies, ix, iy, iz, print_it)
          if(error)return
        endif

        call check_nan('NH3-NH4_N')

!        !-----------------------------------------------------------------------------------
!        !
!        ! Here we start dealing with ammonium nitrates. So far 
!        ! we were dealing with NO3 radicals only
!        !
!        if(vSp(iNH3_1d) < 0.0)then
!          call msg('Negative NH3 before NH3')
!        endif
!        if(vSp(iNH4NO3_1d) < 0.0)then
!          call msg('Negative NH4NO3 before NH3')
!        endif
!
!        !
!        ! Equilibrium HNO3+NH3<=>NH4NO3
!        ! In case of low concentrations, ammonium nitrate desintegrates to ammonia and nitric acid,
!        ! otherwise they stay in equilibrium.
!        !
!        if(vSp(iHNO3_1d) < 0.0)then
!          call msg('Negative HNO3 before tNH3-tHNO3 eq')
!        endif
!
!        call check_nan('Before NH3-NH4')
!
!        tNH3 = vSp(iNH3_1d) + vSp(iNH4NO3_1d)
!        tHNO3 = vSp(iHNO3_1d) + vSp(iNH4NO3_1d)
!
!        if(tNH3 + tHNO3 > fCncLowThresh)then
!          !
!          ! Enough masses, make the equilibrium
!          !
!          CEQUIL = rules%NH3_HNO3_vs_NH4NO3_eq(indT) * (3. -2. * fCloudCover)
!
!          if(vSp(iHNO3_1d) < 0.0)then
!            call msg('Negative HNO3 before HNO3-NH3-NO3 eq')
!          endif
!
!          IF(TNH3 * THNO3 <= CEQUIL)THEN
!          !
!          ! Desintegration of ammonium nitrate, i.e. NO3 and NH4_N go to zero
!          !
!            vSp(iNH3_1d) = tNH3
!            vSp(iHNO3_1d) = tHNO3
!            vSp(iNH4NO3_1d) = 0.
!
!            if(vSp(iHNO3_1d) < 0.0)then
!              call msg('Negative HNO3 after HNO3-NH3-NO3 eq branch 1, TNH3=',  tNH3)
!              call msg('Negative HNO3 after HNO3-NH3-NO3 eq branch 1, THNO3=',  tHNO3)
!            endif
!
!            call check_nan('After NH4 desintegration')
!
!          ELSE
!            !
!            ! Equilibrium. Note possible full consumption of NH3 or HNO3
!            !
!            sqrtmp = sqrt((tHNO3 - tNH3) * (tHNO3 - tNH3) * 0.25 + CEQUIL)
!            differ = (vSp(iHNO3_1d) - vSp(iNH3_1d))*0.5
!            vSp(iHNO3_1d)= sqrtmp + differ
!            vSp(iNH3_1d) = sqrtmp - differ
!            vSp(iNH4NO3_1d) = tNH3 - vSp(iNH3_1d)
!
!            if(vSp(iHNO3_1d) < 0.0)then ! Check for full consumption
!              if(-1.* vSp(iHNO3_1d) < 1.e-5 * tHNO3 .or. (sqrtmp .eps. -1*differ))then
!                vSp(iHNO3_1d) = 0.
!              else
!                call msg('Negative HNO3 after HNO3-NH3-NO3 eq branch 2, sqrtmp=',   sqrtmp)
!                call msg('Negative HNO3 after HNO3-NH3-NO3 eq branch 2, differ=',   differ)
!                call msg('Negative HNO3 after HNO3-NH3-NO3 eq branch 2, fNH3',   tNH3)
!                call msg('Negative HNO3 after HNO3-NH3-NO3 eq branch 2, HNO3:',   vSp(iHNO3_1d))
!                call msg('Negative HNO3 after HNO3-NH3-NO3 eq branch 2, NH4NO3:',   vSp(iNH4NO3_1d))
!                call msg('Negative HNO3 after HNO3-NH3-NO3 eq branch 2, NH3:',   vSp(iNH3_1d))
!              endif
!            endif
!
!            if(vSp(iNH3_1d) < 0.0)then ! Check for full consumption
!              if(-1.* vSp(iNH3_1d) < 1.e-5 * tNH3 .or. (sqrtmp .eps. differ))then
!                vSp(iNH3_1d) = 0.
!              else
!                call msg('Negative NH3 after HNO3-NH3-NO3 eq branch 2, sqrtmp=',   sqrtmp)
!                call msg('Negative NH3 after HNO3-NH3-NO3 eq branch 2, differ=',   differ)
!                call msg('Negative NH3 after HNO3-NH3-NO3 eq branch 2, tNH3',   tNH3)
!                call msg('Negative NH3 after HNO3-NH3-NO3 eq branch 2, HNO3:',   vSp(iHNO3_1d))
!                call msg('Negative NH3 after HNO3-NH3-NO3 eq branch 2, NH4NO3:',   vSp(iNH4NO3_1d))
!                call msg('Negative NH3 after HNO3-NH3-NO3 eq branch 2, NH3:',   vSp(iNH3_1d))
!              endif
!            endif
!
!            if(vSp(iNH4NO3_1d) < 0.0)then ! Check for full consumption
!              if(-1.* vSp(iNH4NO3_1d) < 1.e-6 * tNH3+tHNO3 .or. (tNH3 .eps. vSp(iNH3_1d)))then
!                vSp(iNH4NO3_1d) = 0.
!              else
!                call msg('Negative NH4_N after HNO3-NH3-NO3 eq branch 2, sqrtmp=',   sqrtmp)
!                call msg('Negative NH4_N after HNO3-NH3-NO3 eq branch 2, differ=',   differ)
!                call msg('Negative NH4_N after HNO3-NH3-NO3 eq branch 2, HNO3:',   vSp(iHNO3_1d))
!                call msg('Negative NH4_N after HNO3-NH3-NO3 eq branch 2, NH4NO3:',   vSp(iNH4NO3_1d))
!                call msg('Negative NH4_N after HNO3-NH3-NO3 eq branch 2, NH3:',   vSp(iNH3_1d))
!              endif
!            endif
!
!            call check_nan('After NH3-NH4 equilibrium')
!
!          END IF  ! NH3 + HNO3 <-> NH4NO3 equilibrium
!
!          if(vSp(iHNO3_1d) < 0.0)then
!            call msg('Negative HNO3 after HNO3-NH3-NO3, indT, CEQUIL=',indT, cequil)
!          endif
!        else
!          !
!          ! Not enough masses => all desintegrates to primary pollutants
!          !
!          vSp(iNH3_1d) = vSp(iNH3_1d) + vSp(iNH4NO3_1d)
!          vSp(iHNO3_1d) = vSp(iHNO3_1d) + vSp(iNH4NO3_1d)
!          vSp(iNH4NO3_1d) = 0.
!        endif  ! if enough masses for NH3 <-> NH4_N equilibrium

      endif ! sufficient NOx+NHx mass

    endif  ! ifNitrogen

    !
    ! ATTENTION.
    ! All computations were in concentrations while advection et al require masses.
    ! Thus, we return them to masses here.
    !
!    vSp(:) = vSp(:) * fCellVolume
!    vSpSL(:) = vSpSL(:) * fCellVolume

!    call check_nan('End of chemistry')
!    call check_mass_vector(vSp, garbage, species_list, 'End of chemistry', fLowCncThresh, nSpecies, ix, iy, iz, print_it)

    CONTAINS

    subroutine check_nan(chPlace)
      implicit none
      character(len=*), intent(in) :: chPlace
      integer :: iTmp, jTmp

!      if(sunD > 0. .and. iHourTmp == 5 .and. iLev == 1) call msg(chPlace)
      do iTmp = 1, nSpecies
        !call msg('NaN species:' + fu_str(species_list(iTmp)),iTmp,vSp(iTmp))
        if(isNaN(vSp(iTmp)))then
          do jTmp = 1, nSpecies
            call msg(chPlace + ',' + fu_str(species_list(jTmp)), vSp(jTmp))
          end do
          call msg('NaN found at N2O5:',cN2O5)
          call msg('NaN found at cloud:', fCloudCover)
          call msg('NaN found at iTempr, tempr:',indT,tempr)
          call msg('NaN found at sunD, H2O:', sunD,cH2O)
          stop
        endif
      enddo
    end subroutine check_nan

    subroutine report_vector(chPlace)
      implicit none
      character(len=*), intent(in) :: chPlace
      integer :: iTmp, jTmp

      call msg(chPlace)
      do jTmp = 1, nSpecies
        call msg(fu_connect_strings(chPlace,':',&
                                  & fu_substance_name(species_list(jTmp)),&
                                  & '='), &
               & vSp(jTmp))
      end do
    end subroutine report_vector

  end subroutine transform_acid_basic



  !***********************************************************************************
  !
  ! Acid basic chemistry rules
  !
  !***********************************************************************************


  !********************************************************************************

  function fu_lifetime_acidBasic(rules)result(lifetime)
    !
    ! A typical life time due to degradation and deposition for SOx materials
    ! So far, whatever the species are, we set the "typical" lifetime to be two days
    !
    implicit none

    type(silja_interval) :: lifetime
    type(Tchem_rules_acidBasic), intent(in) :: rules

    lifetime = one_day * 2.0

  end function fu_lifetime_acidBasic


  !*********************************************************************************

  function fu_acidBasic_name(indexSubst) result(chNm)
    !
    ! Returns the name of the substance using its index. Materials must be initialised by 
    ! the time of calling this funciton
    !
    implicit none

    ! return value
    character(len=10) :: chNm

    ! Imported parameter
    integer, intent(in) :: indexSubst

    if(indexSubst == iSO2_1d)then
        chNm = 'SO2'
    elseif(indexSubst == iSO4w_1d_sl)then
        chNm = 'SO4'
    elseif(indexSubst == iH2SO4_1d_sl)then
        chNm = 'H2SO4'
!    elseif(indexSubst == iNH3_1d)then
!        chNm = 'NH3'
!    elseif(indexSubst == iNH4_S_1d)then
!        chNm = 'NH4_S'
!    elseif(indexSubst == iNH4NO3_1d)then
!        chNm = 'NH4NO3'
    elseif(indexSubst == iNO_1d)then
        chNm = 'NO'
    elseif(indexSubst == iNO2_1d)then
        chNm = 'NO2'
    elseif(indexSubst == iNO3rad_1d)then
        chNm = 'NO3rad'
    elseif(indexSubst == iHNO3_1d)then
        chNm = 'HNO3'
    elseif(indexSubst == iPAN_1d)then
        chNm = 'PAN'
!    elseif(indexSubst == iPANrad_1d)then
!        chNm = 'PANrad'
    elseif(indexSubst == iNOP_1d)then
        chNm = 'NOP'
    elseif(indexSubst == iCO_1d)then
        chNm = 'CO'
    elseif(indexSubst == iCH3OOH_1d)then
        chNm = 'CH3OOH'
    elseif(indexSubst == iHCHO_1d)then
        chNm = 'HCHO'
    elseif(indexSubst == iH2O2_1d)then
        chNm = 'H2O2'
    elseif(indexSubst == iCH3O2_1d)then
        chNm = 'CH3O2'
    elseif(indexSubst == iHONO_1d)then
        chNm = 'HONO'
    elseif(indexSubst == iOH_1d_sl)then
        chNm = 'OH'
    elseif(indexSubst == iHO2_1d_sl)then
        chNm = 'HO2'
!    elseif(indexSubst == iC5H8_1d)then
!        chNm = 'C5H8'
!    elseif(indexSubst == iC5H8_2_1d)then
!        chNm = 'C5H8_2'
    else
      call msg('Unknown index:',indexSubst)
      call set_error('Unknown index','fu_acidBasic_name')
      chNm = ''
    endif

  end function fu_acidBasic_name


  !***************************************************************************************
  
  logical function fu_if_specific_dep_acid_basic(rulesAcidBasic)
    implicit none
    type(Tchem_rules_acidBasic), intent(in) :: rulesAcidBasic
    fu_if_specific_dep_acid_basic = .false.     ! no specific deposition whatsoever
  end function fu_if_specific_dep_acid_basic

  !************************************************************************************
  
  logical function fu_if_tla_required_acid_basic(rules) result(required)
    ! Collect transformations' requests for tangent linearization. If tangent linear
    ! is not allowed, corresponding subroutine sets error. 
    implicit none
    type(Tchem_rules_acidBasic), intent(in) :: rules
    
    call set_error('TLA not available for acid basic', 'fu_if_tla_required_acid_basic')
    required = .true.

  end function fu_if_tla_required_acid_basic


END MODULE chem_dep_acid_basic
