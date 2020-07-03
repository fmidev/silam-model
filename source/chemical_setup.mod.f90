module chemical_setup

  use materials
  use names_of_quantities

  implicit none

  ! public routines
  
  ! Basic routines for species class
  !
  public set_species
  public select_species
  public addSpecies
!  public fu_if_species_lists_overlap
  public fu_set_mode
  public fu_species_match ! defined, have same material and modes match
  public fu_modes_match ! both gas or both no_mode or have same nominal diameter 
  public fu_str      ! fu_str
  public fu_species_from_short_string
  public decode_id_params_from_io_str
  public fu_optical_wave_length
  public set_optical_wave_length
  public fu_mode
  public fu_massmean_d
  public fu_nominal_d
  public fu_max_d
  public fu_min_d
  public fu_index   ! item in the list of items
  public put_species_to_namelist
  public fu_material

  !
  ! mapping
  !
  public create_chemical_mapper
  public copy_chemical_mapper
  public free_chemical_mapper
  public create_adaptor
  public chemical_adaptor_missing
  public destroy_adaptor

  public fu_nSubst
  public fu_name
  public fu_iSp
  public fu_iSubst

  !
  ! Cocktail descriptor
  !
  public set_cocktail_description
  public set_descriptor_from_species  
  public copy_descr_to_cocktail_link
  public copy_cocktail_descriptor
  public normalise_cocktail_description
  public fu_cocktail_descriptior_ptr
  public fu_conversion_factor
  public read_standard_cocktail_descrs
!  public set_standard_cocktails_fname
  public get_inventory
  public collect_descriptor_emis_rates
  public link_src_species_to_given_list
  public fu_basic_mass_unit
  public fu_species

  !
  ! Aerosol
  !
  public fu_set_aerosol_rules
  public fu_if_humidity_dependent
  public fu_default_humidity
  public defined
  public report
  public fu_ifMeanDiamDynamic
  public fu_mode_type
  public fu_substance_name
  public set_missing
  public fu_nModes
  public fu_nWaves
  public fu_nSpecies
  public fu_emis_fractions
  public set_aerosol
  public fu_lifetime

  public set_aerosol_mode
  public set_mode_name
  public set_massmean_d
  public fu_aerosol_mode_size_to_str
  public fu_optical_wave_length_to_str

  public fu_index_of_variable
  public fu_emis_cocktail_mapping
  public fu_emis_cocktail_scaling
  public set_speciesReferences
  public kill_speciesReferences
  public create_mode_projection
  public fu_integrate_volume
  public fu_integrate_number
  public fu_mass_mean_diam
  public destroy

  ! private routines
  private report_species
  private report_species_as_namelist
  private report_species_list
  private fu_species_to_short_string
  private set_species_from_params
  private fu_optical_wave_length_species
  private set_optical_wave_length_sp
  private fu_massmean_d_species
  private fu_min_d_species
  private fu_max_d_species
  private fu_mode_type_species
  private fu_substance_name_species
  private fu_index_of_species
  private fu_species_equal
  private fu_species_defined
  private fu_adaptor_defined
  private fu_iSubst_name
  private fu_iSubst_ind
  private fu_material_species
  private copy_species
  private enlarge_species_list
  private fu_set_mode_from_params
  private fu_set_mode_from_namelist

  private fu_nwaves_mapper
  private fu_nModes_mapper
  private fu_substance_name_mapper
  private fu_iSp_wl
  private fu_iSp_short
  private fu_iSp_from_values
  private fu_iSp_from_species
  private fu_descr_defined
  private cocktail_descriptor_report

  !
  ! Private cocktail descriptor
  !

  private read_cocktail_description
  private read_one_description_v4_6
  private reset_used_cocktail_descr   ! actually, deallocate
  private reset_new_cocktail_descr
  private fu_nSpecies_descr
  private fu_emis_fractions_of_descr
  private fu_descriptors_equal
  private fu_descr_index_in_descr_by_name
  private fu_conversion_factor_descr
  private get_inventory_descr
  private fu_basic_mass_unit_descr
  private fu_name_of_descriptor
  private fu_species_of_descriptor
  private link_emis_descr_to_species
  !
  ! Private aerosol
  !  
  private set_aerosol_missing
  private set_aerosol_rules_missing
  private set_aerosol_from_copy
  private set_aerosol_mode_from_params
  private set_aerosol_mode_from_d
  private fu_compare_modes_bigger
  private fu_compare_modes_eq
  private fu_compare_modes_mean_diam
  private fu_compare_modes_general
  private set_aerosol_from_namelist
  private fu_n_modes_of_spectrum
  private fu_mode_of_aerosol
  private fu_aerosol_defined
  private report_aerosol
  private report_aerosol_mode
  private fu_lifetime_aerosol
  private set_chem_mapper_missing
  private fu_mode_from_species
  private fu_min_d_mode
  private fu_max_d_mode
  private fu_mode_name
  private fu_mode_type_mode
  private fu_massmean_d_mode
  private fu_mode_defined
  private sort_modes

  interface operator(==)
     module procedure fu_species_equal
     module procedure fu_compare_modes_eq
     module procedure fu_compare_aerosols_eq
     module procedure fu_descriptors_equal
  end interface

  interface assignment (=)
    module procedure copy_species
    module procedure copy_adaptor
  end interface
  
  interface operator (>)
     module procedure fu_compare_modes_bigger
  end interface

  interface set_species
    module procedure set_species_from_params
  end interface

  interface fu_str
    module procedure fu_species_to_short_string
  end interface
  
  interface fu_set_mode
    module procedure fu_set_mode_from_params
    module procedure fu_set_mode_from_namelist
  end interface

  interface defined
     module procedure fu_species_defined
     module procedure fu_aerosol_defined
     module procedure fu_descr_defined
     module procedure fu_adaptor_defined
     module procedure fu_mode_defined
  end interface

  interface report
     module procedure report_aerosol
     module procedure report_aerosol_mode
     module procedure report_species
     module procedure report_species_list
     module procedure report_species_as_namelist
     module procedure cocktail_descriptor_report
  end interface
  

  interface fu_index
     module procedure fu_index_of_species
     module procedure fu_descr_index_in_descr_by_name
  end interface

  interface fu_material
    module procedure fu_material_species
  end interface

  interface get_inventory
     module procedure get_inventory_descr
  end interface

  interface fu_optical_wave_length
     module procedure fu_optical_wave_length_species
  end interface

  interface set_optical_wave_length
     module procedure set_optical_wave_length_sp
  end interface

  interface fu_mode
    module procedure fu_mode_from_species
    module procedure fu_mode_of_aerosol
  end interface
  
  interface fu_mode_type
    module procedure fu_mode_type_species
    module procedure fu_mode_type_mode
  end interface
  
  interface fu_massmean_d
     module procedure fu_massmean_d_species
     module procedure fu_massmean_d_mode
  end interface

  interface fu_nominal_d
     module procedure fu_nominal_d_species
     module procedure fu_nominal_d_mode
  end interface

  interface fu_max_d
     module procedure fu_max_d_species
     module procedure fu_max_d_mode
  end interface

  interface fu_min_d
     module procedure fu_min_d_species
     module procedure fu_min_d_mode
  end interface


  interface fu_density
    module procedure fu_density_species
  end interface


  interface fu_iSp
     module procedure fu_iSp_wl
     module procedure fu_iSp_short
     module procedure fu_iSp_from_values
     module procedure fu_iSp_from_species
  end interface

  interface fu_iSubst
     module procedure fu_iSubst_name
     module procedure fu_iSubst_ind
  end interface

  interface fu_substance_name
     module procedure fu_substance_name_mapper
     module procedure fu_substance_name_species
  end interface

  interface set_missing
     module procedure set_aerosol_missing
     module procedure set_aerosol_rules_missing
     module procedure reset_new_cocktail_descr
     module procedure set_chem_mapper_missing
  end interface

  interface set_aerosol_mode
    module procedure set_aerosol_mode_from_params
    module procedure set_aerosol_mode_from_d
  end interface

  interface fu_nModes
     module procedure fu_nModes_mapper
     module procedure fu_n_modes_of_spectrum  ! aerosol
  end interface
  
  interface fu_nWaves
     module procedure fu_nWaves_mapper
  end interface
  
  interface fu_nSpecies
     module procedure fu_nSpecies_descr
  end interface

  interface fu_conversion_factor
    module procedure fu_conversion_factor_descr
  end interface

  interface fu_basic_mass_unit
    module procedure fu_basic_mass_unit_descr
  end interface

  interface fu_name
    module procedure fu_name_of_descriptor
    module procedure fu_mode_name
  end interface

  interface fu_species
    module procedure fu_species_of_descriptor
  end interface

  interface fu_emis_fractions
     module procedure fu_emis_fractions_of_descr
  end interface

  interface set_aerosol
     module procedure set_aerosol_from_namelist
     module procedure set_aerosol_from_copy
  end interface
  
  interface fu_lifetime
     module procedure fu_lifetime_aerosol
  end interface

  interface destroy
     module procedure reset_used_cocktail_descr
  end interface

  ! 
  ! A structure for defining the aerosol used as a sub-element of the silam_species.
  ! Note on solubility: 0-insoluble, 1..100 refers to solubility. May be useful for wet deposition.
  !
  type Taerosol_mode
     private
     character(len=substNmLen) :: name = ''
     integer :: distr_type =int_missing       ! distribution type or gas_phase_flag
     integer :: solubility = int_missing ! 0=insoluble, 1=soluble, later: 1..100 reflects the solubility
     ! Shape parameters:
     !
     real :: fp1 = real_missing, fp2 = real_missing
     real :: mass_mean_d = real_missing  !To be used in calculations
     real :: nominal_d = int_missing ! to match species and modes (within 1nm)
     logical :: ifMeanDiamDynamic = .false.
     ! Bins: 
     ! fp1 = lower limit, fp2 = max limit
     ! Lognormal modes:
     ! fp1 = (number) median, fp2 = sigma, 
     ! 
     ! ip : reservation for gamma distribution
     integer :: ip1 = int_missing
     type(silja_logical) :: defined = silja_false
  end type Taerosol_mode
  
  ! Shapes of the aerosol size distribution density
  integer, public, parameter :: gas_phase_flag = 7000
  integer, public, parameter :: fixed_diameter_flag = 7001
  integer, public, parameter :: gamma_function_flag = 7002
  integer, public, parameter :: moving_diameter_flag = 7003
  integer, public, parameter :: lognormal_flag = 7004
  integer, public, parameter :: sea_salt_mode_flag = 7005
  integer, public, parameter :: fire_mode_flag = 7006
  integer, public, parameter :: no_mode_flag = 7007

  type(Taerosol_mode), parameter, public :: aerosol_mode_missing = Taerosol_mode( '', &
                                                         & int_missing, int_missing, &
                                                         & real_missing, &
                                                         & real_missing, real_missing, real_missing,&
                                                         & .false., &
                                                         & int_missing, silja_false)
  !
  ! The gas-phase "mode" is currently implemented using a special
  ! type-flag. The below constant can be used for constructing
  ! gas-phase species.
  !
  type(Taerosol_mode), parameter, public :: in_gas_phase = Taerosol_mode( '', &
                                                      & gas_phase_flag, 0, &
                                                      & real_missing, &
                                                      & real_missing, real_missing, real_missing, &
                                                      & .false., &
                                                      & int_missing, silja_true)

  ! For strange quantities like PM10, PM25  total dust, total seasalt etc...
  type(Taerosol_mode), parameter, public :: no_mode = Taerosol_mode( '', &
                                                      & no_mode_flag, 0, &
                                                      & real_missing, &
                                                      & real_missing, real_missing, real_missing, &
                                                      & .false., &
                                                      & int_missing, silja_true)
  ! 
  ! The silam species. Contains complete chemical and physcial
  ! information of one consituent. Also optical species are possible
  ! with the wavelength field.
  !
  ! The species are typically handled in an array (eg. in descriptors
  ! and cocktails). In almost all cases the species array arguments
  ! for the subroutines in this module are assumed to contain only
  ! unique species, but this is not checked explicitly. Therefore, the
  ! preferred way of composing such lists is using the addSpecies
  ! subroutine which does the checking.
  !
  type silam_species
     type(Taerosol_mode) :: mode = aerosol_mode_missing
     real :: wavelength = real_missing
     type(silam_material), pointer :: material => null()
     type(silja_logical) :: defined = silja_false
  end type silam_species
  type silam_species_arr_ptr
    type(silam_species), dimension(:), pointer :: pArSp
  end type silam_species_arr_ptr
  
  type (silam_species), parameter, public :: species_missing = &
      & silam_species(aerosol_mode_missing, real_missing, null(), silja_false)

  ! 
  ! The type for an overarching aerosol definition. I did not find use
  ! for it in cocktails, but it remains in use for passing information
  ! to some routines, as well as for reading descriptors from a file.
  !
  TYPE Taerosol
     !PRIVATE
     integer :: n_modes =int_missing !, shape_p
     type(Taerosol_mode), dimension(:), pointer :: modes => null()
     type(silja_logical) :: defined = silja_false
  END TYPE Taerosol
  TYPE (Taerosol),  parameter, public :: aerosol_missing = &
      & Taerosol(int_missing, null(), silja_false)
  !
  ! Aerosol has own rules for computation - first of all, scavenging
  ! Note that it is public since nearly all cocktails need to know
  ! them. The content now overlaps with deposition_rules defined in
  ! the respective module, however, this might turn out to be
  ! necessary.
  !
  type Taerosol_rules
    !    private
    logical :: ifHumidityDependent     ! or just take the defualt humidity
    real :: fDefaultRelHumidity        ! the very default humidity
    logical :: ifNumberCnc
    type(silja_logical) :: defined
  end type Taerosol_rules

  type(Taerosol_rules), parameter, public :: aerosolRules_missing = Taerosol_rules( &
                                                                       & .false., &
                                                                       & real_missing, &
                                                                       & .false., &
                                                                       & silja_false)
  !
  ! Chemical mapper. The preferred way (in my opinion, JV 9/2009) of
  ! accessing the species in the run should by using the 1D species
  ! index, since 1) this is used in mass_maps and 2) the list of
  ! species, accessed with the species index, usually contains the all
  ! the necessary information. However, in some situations the
  ! subtance-mode-wavelength-triplet (or substance and mode) is still
  ! necessary and this is now handled using a generic structure
  ! defined below.
  !
  ! The mapper is created by giving an array of silam species. After
  ! this, eg. following is possible:
  ! do iSubst = 1, fu_nSubst(mapper)
  !  do iMode = 1, fu_nModes(iSubst, mapper)
  !   do iWave = 1, fu_nWaves(iSubst, iMode, mapper)
  !   iSpecies = fu_iSp(iSubst, iMode, iWave, mapper)
  !   ....
  !
  ! The structure allows arbitrary number of modes for each substance
  ! and an arbitrary number of wavelengths for each
  ! substance-wavelength combination, though this is not used
  ! currently (the number if same for all modes). Consequently, the
  ! structure is somewhat complicated, which mandates using the the
  ! accessor functions whenever possible. This also ensures that
  ! nonexistent subst-mode-wave combinations are not sent to fu_isp,
  ! which for speed reasons does not contain such checks.
  !
  type chemical_mapper
     character(len=substNmLen), dimension(:), pointer :: names
     type(iarrptr), dimension(:), pointer :: iSp    ! for substance, mode and wave
     integer, dimension(:), pointer :: iSubst       ! (iSpecies)
     integer, dimension(:), pointer :: nmodes       ! for substances
     type(iarrptr), dimension(:), pointer :: nwaves ! for substance and mode
     integer :: nSubstances
     integer :: nSpecies
     type(silja_logical) :: defined != .false.
  end type chemical_mapper

  type (chemical_mapper), public, parameter :: chemical_mapper_missing = &
        & chemical_mapper(null(), null(), null(), null(), null(),&
                         & int_missing, int_missing, silja_false)

  ! The cocktail descriptor. Largely unchanged from v4.5 except for
  ! one thing: all arrays now have dimension of descr%nspecies, and
  ! they are indexed by the species index of the descriptor. For
  ! convenience, the descriptor has a chemical_mapper for traditional
  ! subst-mode indexing.
  !
  type Tcocktail_descr
     private
     integer :: nSpecies
     character(len=substNmLen) :: cocktail_name                 ! Name of the mixture
     character(len=unitNmLen) :: chFractionMassUnit !!! _BASIC_ unit of one specified in the cocktail file
     real, dimension(:), pointer :: mass_fractions              ! (nSpecies)
     type(silam_species), dimension(:), pointer :: species      ! (nSpecies)
     integer, dimension(:), pointer :: iEmisCocktSpeciesMapping  ! (nspecies)
     real, dimension(:), pointer :: f_descr_unit_2_species_unit  ! (nSpecies)
     type(silja_logical) :: defined 
     logical :: if_normalise
  end type Tcocktail_descr

  ! 
  ! The chemical adaptor is actually just array with dimension
  ! adaptor%nspecies, where every element contains the index of a
  ! species of a smaller list in a bigger species list.
  !
  type chemical_adaptor
     integer, dimension(:), allocatable :: iSp
     integer :: nSpecies
     type(silja_logical) :: defined != .false.
  end type chemical_adaptor

  type iarrptr
     integer, dimension(:), pointer :: ptr1d
     type(iarrptr), dimension(:), pointer :: ptr2d
  end type iarrptr

 !----------------------------------------------------------------------
  !
  ! The type for the chemical setup of the run. It consists of four lists: 
  ! - Native, non-sorted full list of species with reference to materials and to all 
  !   three main cocktails: Emis, Transp, OptDns
  ! - Emission, which contains only emission species sorted as in the emission
  !   cocktail, with references to materials and to transport cocktail
  ! - Transport, which contains only transported species, in corresponding order and
  !   with references to both emission and OptDns cocktails
  ! - OptDns, which contains only OptDns species, with references to materials and to
  !   transport cocktail
  !
  ! A reference from a substance to another cocktail means a list of indices in that 
  ! cocktail and corresponding fractions (summing-up to 1), so that a step from this
  ! substance to another cocktail is just a single cycle over the list
  !
  type TspeciesReference      ! basic reference type
    integer :: nRefSpecies
    integer, dimension(:), pointer :: indSpeciesTo  => null()
    real, dimension(:), pointer :: fract  => null()
  end type TspeciesReference

  type TspeciesReferencePtr
     type(TspeciesReference), pointer :: ptr
  end type TspeciesReferencePtr

  !
  ! List of standard cocktails
  !
  type(Tcocktail_descr), dimension(:), pointer, private, save :: StandardCocktailPtr => null() 
!  !
!  ! The PATH of the cocktail description files
!  !
!  character(len=fnlen),private,save :: chCocktailDescrFile





CONTAINS
   
  !**********************************************************************************
  !
  ! Encapsulation of the silam_species and aerosol_mode types
  !
  !**********************************************************************************

  function fu_optical_wave_length_species(species) result(wavelength)
    type(silam_species) :: species
  
    real :: wavelength
    wavelength = species%wavelength
  end function fu_optical_wave_length_species

  !********************************************************************************

  subroutine set_optical_wave_length_sp(species, fWaveLength)
    type(silam_species) :: species
    real, intent(in) :: fWaveLength
  
    species%wavelength = fWaveLength
  end subroutine set_optical_wave_length_sp

  !************************************************************************************

  function fu_material_species(species) result(material)
    type(silam_species) :: species
  
    type(silam_material), pointer  :: material

    material => species%material
  end function fu_material_species


  !************************************************************************************

  subroutine set_species_from_params(species, material, aerosol_mode, wavelength_)
    ! 
    ! Defined a new species structure with given aerosol mode,
    ! material, and optionally, wavelength.
    ! 
    implicit none
    type(silam_species), intent(out) :: species
    type(Taerosol_mode), intent(in) :: aerosol_mode
    type(silam_material), pointer :: material
    real, intent(in), optional :: wavelength_

    if (.not. associated(material)) then
      call set_error('Material not associated','set_species_from_params')
      return
    end if
!    if (.not. defined(aerosol_mode)) then
!      call set_error('Undefined aerosol mode', 'set_species_from_params')
!      return
!    end if

    species%material => material
    species%mode = aerosol_mode

    if(present(wavelength_))then
      species%wavelength = wavelength_
    else
      species%wavelength = real_missing
    end if

    species%defined = silja_true

  end subroutine set_species_from_params


  !---------------------------------------------------------------------------
  !
  ! Accessrors for the aerosol properties
  !
  !==================================================================================

  function fu_mode_from_species(species) result(aerosol_mode)
    implicit none
    type(silam_species), intent(in) :: species
    type(Taerosol_mode) :: aerosol_mode
    aerosol_mode = species%mode
  end function fu_mode_from_species

  !================================================================================
  function fu_massmean_d_species(species) result(massmean_d)
    implicit none
    type(silam_species), intent(in) :: species

    real :: massmean_d

    if (.not. species%mode%defined) then
      call set_error('Undefined mode', 'fu_massmean_d_species')
      return
    end if
    
    massmean_d = fu_massmean_d(species%mode)
    
  end function fu_massmean_d_species
  !================================================================================
  function fu_nominal_d_species(species) result(nominal_d)
    implicit none
    type(silam_species), intent(in) :: species

    real :: nominal_d

    if (.not. species%mode%defined) then
      call set_error('Undefined mode', 'fu_nominal_d_species')
      return
    end if
    
    nominal_d = species%mode%nominal_d
    
  end function fu_nominal_d_species

  !==================================================================================
  real function fu_min_d_species(species) result(min_d)
    implicit none
    type(silam_species), intent(in) :: species

    if (.not. species%mode%defined) then
      call set_error('Undefined mode', 'fu_min_d_species')
      return
    end if

    min_d = fu_min_d(species%mode)

  end function fu_min_d_species

  !==================================================================================
  real function fu_max_d_species(species) result(max_d)
    implicit none
    type(silam_species), intent(in) :: species

    if (.not. species%mode%defined) then
      call set_error('Undefined mode', 'fu_max_d_species')
      return
    end if
    
    max_d = fu_max_d(species%mode)

  end function fu_max_d_species

  !==================================================================================
  integer function fu_mode_type_species(species) result(mType)
    implicit none
    type(silam_species), intent(in) :: species
    mType = species%mode%distr_type
  end function fu_mode_type_species

  !==================================================================================
  real function fu_density_species(species) result(density)
    implicit none
    type(silam_species), intent(in) :: species

    if (.not. species%defined == silja_true) then
      call set_error('Undefined species', 'fu_density_species')
      return
    end if

    if (species%mode == in_gas_phase) then
      density = real_missing
      return
    end if
    
    density = fu_dry_part_density(species%material)

  end function fu_density_species

  !==================================================================================

  function fu_substance_name_species(species) result(chName)
    implicit none
    character(len=substNmLen) :: chName
    type(silam_species), intent(in) :: species

    if (defined(species)) then
      chName = fu_name(species%material)
    else
      chName = ''
    end if

  end function fu_substance_name_species

  !==================================================================================

  subroutine select_species(species_in, nSpecies_in, substNm, mode, wavelength, indices, nSelected)
    !
    ! Select species from a list using given criteria and return their
    ! indices. Any of the conditions (susbtNm, mode, wavelength) can
    ! be set to missing, in which case any species is taken to match.
    !
    implicit none
    type(silam_species), dimension(:), intent(in) :: species_in
    integer, intent(in) :: nSpecies_in
    character(len=*) :: substNm
    type(Taerosol_mode), intent(in) :: mode
    real, intent(in) :: wavelength
    
    integer, dimension(:), intent(inout) :: indices
    integer, intent(out) :: nSelected
    
    integer :: i
    logical :: ifselect
    
    nSelected = 0
    
    do i = 1, nSpecies_in
      if(.not.associated(species_in(i)%material))cycle
      ifselect = .true.
      ifselect = ifselect .and. (substNm == char_missing .or. &
                               & trim(fu_name(species_in(i)%material)) == trim(substNm))
      ifselect = ifselect .and. (mode == aerosol_mode_missing .or. &
                & species_in(i)%mode == mode .or. &
                &( species_in(i)%mode%nominal_d > 0. .and. &
                   & abs(species_in(i)%mode%nominal_d - mode%nominal_d) < mode_diam_tolerance) )
      if (.not. (wavelength .eps. real_missing)) then
        ifselect = ifselect .and. (wavelength .eps. species_in(i)%wavelength)
      end if

      if (ifselect) then
        nselected = nselected + 1
        if(nselected > size(indices))then
          call msg('Too small size of the indices array:',size(indices))
          call set_error('Too small size of the indices array','select_species')
          return
        endif
        indices(nselected) = i
      end if

    end do

  end subroutine select_species


  !==================================================================================

  subroutine report_species_list(species_list)
    implicit none
    type(silam_species), dimension(:), intent(in) :: species_list

    integer :: ind_species
    
    do ind_species = 1, size(species_list)
      call report(species_list(ind_species))
    end do

  end subroutine report_species_list
  
  !************************************************************************************

  subroutine report_species(species)
    implicit none
    type(silam_species), intent(in) :: species
    
    type(silam_sp):: line

    if (fu_undefined(species%defined)  ) then ! neither true nor false
            call set_error("Can't tell if species defined", "report_species")
    endif

    if (.not. defined(species)) then
      call msg('Undefined species')
      return
    end if

    line%sp => fu_work_string()
!    if(error)return

    if (species%wavelength .eps. real_missing) then
      if (.not. defined(species%mode)) then
        write(line%sp,fmt='("Material: ", A, ", mode ", A, " UNDEFINED mode")') trim(fu_name(species%material))
      elseif (species%mode == in_gas_phase) then
        write(line%sp,fmt='("Material: ", A, ", mode ", A, " gas phase, ", A)') trim(fu_name(species%material)), trim(species%mode%name)
      elseif (species%mode%distr_type == fixed_diameter_flag .or. &
             & species%mode%distr_type == moving_diameter_flag) then
        write(line%sp,fmt='("Material: ", A, ", mode ", A, " fixed/moving diam mode (um): ",3F6.2, " Nominal diam (um): ",F6.2,", solubility:", I7)') &
             & trim(fu_name(species%material)), &
             & trim(species%mode%name), &
             & species%mode%fp1*1e6, species%mode%fp2*1e6, species%mode%mass_mean_d*1e6, &
             & species%mode%nominal_d*1e6, &
             & species%mode%solubility
      elseif (species%mode%distr_type == lognormal_flag) then
        write(line%sp,fmt='("Material: ", A, ", mode ", A, " lognormal mode (um): ",F6.2, " Nominal diam (um): ",F6.2, ", solubility:", I7)') &
             & trim(fu_name(species%material)), &
             & trim(species%mode%name), &
             & species%mode%mass_mean_d*1e6, &
             & species%mode%nominal_d*1e6, &
             & species%mode%solubility
      else
        write(line%sp,fmt='("Material: ", A, ", mode ", A, " UNKNOWN, solubility:", I7)') &
             & trim(fu_name(species%material)), &
             & trim(species%mode%name), &
             & species%mode%solubility
      end if

    else

      if (species%mode == in_gas_phase) then
        write(line%sp,fmt='("Material: ", A, ", mode ", A, " gas phase, wavelength (nm): ", F6.1)') &
             & trim(fu_name(species%material)), &
             & trim(species%mode%name), &
             & species%wavelength*1e9
      else if (species%mode%distr_type == fixed_diameter_flag .or. &
             & species%mode%distr_type == moving_diameter_flag) then
        write(line%sp,fmt='("Material: ",A, ", mode ", A, " fixed/moving diam mode:",4F6.2, ", solubility:",I7, ", wavelength (nm):" F6.1)') &
             & trim(fu_name(species%material)),&
             & trim(species%mode%name),  &
             & species%mode%fp1*1e6,  &
             & species%mode%fp2*1e6, &
             & species%mode%mass_mean_d*1e6, &
             & species%mode%nominal_d*1e6, species%mode%solubility, species%wavelength*1e9
      else if (species%mode%distr_type == lognormal_flag) then
        write(line%sp, fmt='(5A, F6.2, A, I3, A, F6.1)') &
             & 'Material ', &
             & trim(fu_name(species%material)), &
             & ' ', trim(species%mode%name), &
             & ', lognormal mode (um): ', &
             & species%mode%mass_mean_d*1e6, &
             & ', nominal d (um): ', &
             & species%mode%nominal_d*1e6, &
             & ', solubility: ', &
             & species%mode%solubility, &
             & ', wavelength: ', &
             & species%wavelength*1e9
      elseif (species%mode%distr_type == no_mode_flag) then
        write(line%sp,fmt='("Material: ", A, ", mode ", A, " NO_MODE")') &
             & trim(fu_name(species%material)), &
             & trim(species%mode%name)
      else
        write(line%sp,fmt='("Material: ", A, ", mode ", A, " UNKNOWN mode, solubility:", I7, ", wavelength (nm):" F6.1)') &
             & trim(fu_name(species%material)), &
             & trim(species%mode%name), &
             & species%mode%solubility, &
             & species%wavelength*1e9
      end if

    end if
    if(species%mode%ifMeanDiamDynamic) then 
      line%sp = line%sp + ', dynamic mean diameter,'
    endif
    
    line%sp = line%sp + fu_str(species)
    
    call msg(line%sp)
    call free_work_array(line%sp)

  end subroutine report_species


  !==================================================================================
  
  subroutine put_species_to_namelist(species, nl) 
    !
    ! Reports species in the form of namelist to the given unit. NOT duplicated to the
    ! log file.
    !
    implicit none
    
    type (Tsilam_namelist), pointer :: nl
    ! Imported parameters
    character (len=*), parameter :: subname="put_species_to_namelist"
    type(silam_species), intent(in) :: species
    real :: fTmp

    if (.not. associated(nl)) then
       call set_error("not associated namelist", subname)
    endif

    if(defined(species))then
      call add_namelist_item(nl, 'substance_name', trim(fu_name(species%material)))

      ! If this thing is set in a file on input, silam units to be assumed: 
      ! no netcdf scaling applied. the one from a nametabe is still possible though
      call add_namelist_item(nl,'silam_amount_unit', fu_basic_mass_unit(species%material))

      fTmp = fu_mole_mass(species%material)
      if (fTmp < 9999  .and. fTmp  >0) then
          call add_namelist_item(nl, 'molar_mass', fu_str(fTmp)//' kg/mole' )
      endif

      if(species%WaveLength /= real_missing) &
           call add_namelist_item(nl, 'optical_wavelength', fu_str(species%WaveLength*1.e6)//' um')

      if(species%mode%defined == silja_true) then
        call add_namelist_item(nl, 'mode_name', species%mode%name)
        if(species%mode%solubility /= int_missing) &
            & call add_namelist_item(nl, 'mode_solubility', fu_str(species%mode%solubility))

        select case(species%mode%distr_type)
         case(gas_phase_flag)
           call add_namelist_item(nl, 'mode_distribution_type', 'GAS_PHASE' )
           
         case(fixed_diameter_flag)
           call add_namelist_item(nl, 'mode_distribution_type', 'FIXED_DIAMETER' )
           call add_namelist_item(nl, 'mode_nominal_diameter', fu_str(species%mode%nominal_d*1e6)//' um')
           call add_namelist_item(nl, 'fix_diam_mode_min_diameter',  fu_str(species%mode%fp1*1.e6)//' um')
           call add_namelist_item(nl, 'fix_diam_mode_max_diameter',  fu_str(species%mode%fp2*1.e6)//' um')
           if( species%mode%mass_mean_d /= real_missing) then
             call add_namelist_item(nl, 'fix_diam_mode_mean_diameter',  &
                     & fu_str(species%mode%mass_mean_d*1.e6)//' um')
           end if

         case (moving_diameter_flag)
           call add_namelist_item(nl, 'mode_distribution_type', 'MOVING_DIAMETER' )
           call add_namelist_item(nl, 'mode_nominal_diameter', fu_str(species%mode%nominal_d*1e6)//' um')
           call add_namelist_item(nl, 'moving_diam_mode_min_diameter',  fu_str(species%mode%fp1*1.e6)//' um')
           call add_namelist_item(nl, 'moving_diam_mode_max_diameter',  fu_str(species%mode%fp2*1.e6)//' um')
           if( species%mode%mass_mean_d /= real_missing) then
             call add_namelist_item(nl, 'moving_diam_mode_mean_diameter',  &
                     & fu_str(species%mode%mass_mean_d*1.e6)//' um')
           end if

         case(gamma_function_flag)
           call add_namelist_item(nl, 'mode_distribution_type', 'GAMMA_FUNCTION' )
           call add_namelist_item(nl, 'mode_nominal_diameter', fu_str(species%mode%nominal_d*1e6)//' um')
           call add_namelist_item(nl, 'gamma_mode_gamma_number', fu_str(species%mode%ip1) )

         case(lognormal_flag)
           call add_namelist_item(nl, 'mode_distribution_type', 'LOGNORMAL' )
           call add_namelist_item(nl, 'mode_nominal_diameter', fu_str(species%mode%nominal_d*1e6)//' um')
           call add_namelist_item(nl, 'lognormal_mode_gamma_number', fu_str(species%mode%fp1*1e6)//' um' )
           call add_namelist_item(nl, 'lognormal_mode_gamma_number', fu_str(species%mode%fp2*1e6)//' um' )

         case(sea_salt_mode_flag)
           call add_namelist_item(nl, 'mode_distribution_type', 'SEA_SALT' )
           call add_namelist_item(nl, 'mode_nominal_diameter', fu_str(species%mode%nominal_d*1e6)//' um')
           call add_namelist_item(nl, 'sea_salt_mode_gamma_number', fu_str(species%mode%fp1*1e6)//' um' )
           call add_namelist_item(nl, 'sea_salt_mode_gamma_number', fu_str(species%mode%fp2*1e6)//' um' )

         case(no_mode_flag)
           call add_namelist_item(nl, 'mode_distribution_type', 'NO_MODE' )

         case default
           call msg('Unknown mode type:',species%mode%distr_type)
           call set_error('Unknown mode type',subname)
           return
       end select
      endif
   else
      call set_error("Trying to make a list of undefined species",subname)
   endif

  end subroutine put_species_to_namelist

  !==================================================================================
  
  subroutine report_species_as_namelist(species, iUnit)
    !
    ! Reports species in the form of namelist to the given unit. NOT duplicated to the
    ! log file.
    !
    implicit none
    
    ! Imported parameters
    type(silam_species), intent(in) :: species
    integer, intent(in) :: iUnit
    
    if(defined(species))then
      write(iUnit,fmt='(3x,2A)')'substance_name = ', trim(fu_name(species%material))
      if(defined(fu_mode(species))) call report_as_namelist_aer_mode(fu_mode(species), iUnit)
      if(.not. (species%WaveLength .eps. real_missing)) &
           & write(iUnit,fmt='(3x,A,E15.7,A)')'optical_wavelength = ',species%WaveLength*1.e6,' mkm'
    else
      write(iUnit,*)'substance_name = '
    endif

  end subroutine report_species_as_namelist


  !==================================================================================

  integer function fu_index_of_species(species, speciesList, nspecies_, ifByName_) result(ind)
    ! 
    ! Return the index of species in list speciesList optionally
    ! scanning only nspecies_ first species. The species in the list
    ! are assumed to be unique.
    !
    implicit none
    type(silam_species), intent(in) :: species
    type(silam_species), dimension(:), intent(in) :: speciesList
    integer, intent(in), optional :: nspecies_
    logical, intent(in), optional :: ifByName_

    integer :: i, nspecies
    logical ifByName
    
    if (present(nspecies_)) then
      nspecies = nspecies_
    else
      nspecies = size(speciesList)
    end if
    if(present(ifByName_))then
      ifByName = ifByName_
    else
      ifByName = .false.
    endif

    ind = int_missing

    if(ifByName)then
      do i = 1, nspecies
        if (fu_str(speciesList(i)) == fu_str(species)) then
          ind = i
          return
        end if
      end do
    else
      do i = 1, nspecies
        if (speciesList(i) == species) then
          ind = i
          return
        end if
      end do
    endif

  end function fu_index_of_species


  !==================================================================================

  function fu_index_of_variable(quantity, quantityList, species, speciesList, nVars_) result(ind)
    ! 
    ! Return the index of the given variable (= quantity + species) in the list of 
    ! variables (= quantityList + speciesList), optionally scanning only nspecies_ first 
    ! species. The variables in the list are assumed to be unique.
    ! 
    ! To use without speciesList fu_index_of_integer should be called
    !
    implicit none
    integer, intent(in) :: quantity
    integer, dimension(:), intent(in) :: quantityList
    type(silam_species), intent(in) :: species
    type(silam_species), dimension(:), intent(in) :: speciesList
    integer, intent(in), optional :: nVars_
    
    integer :: ind, i, nVars
    
    if (present(nVars_)) then
      nVars = nVars_
    else
      nVars = size(quantityList)
    end if
    if (nvars > size(speciesList)) call ooops("Achtung")

    ind = int_missing

    do i = 1, nVars
      if(quantity == quantityList(i))then
        if (species == speciesList(i)) then
          ind = i
          return
        end if
      endif
    end do

  end function fu_index_of_variable


  !==================================================================================

  logical function fu_species_equal(species1, species2) result(equal)
    !
    ! Compare equality of two species. This requires the same material
    ! (name) and the same mode (as defined by compare_modes), and if
    ! wavelength is not missing, the same wavelength.
    !
    implicit none
    type(silam_species), intent(in) :: species1, species2

    equal = .false.
    !
    ! Comparison is possible only if both species are defined 
    !
    if(species1%defined == silja_true)then
      if(.not. species2%defined == silja_true)then
        equal = .false.
        return
      endif
    else
      equal = species1%defined == species2%defined
      return
    endif
    
    if (.not. associated(species1%material, species2%material)) return
    !if(.not. fu_name(species1%material) == fu_name(species2%material)) return
    if(.not. (species1%mode == species2%mode)) return

    if (species1%wavelength .eps. real_missing) then
      if(.not. (species2%wavelength .eps. real_missing))return
    else
      if (.not. (species1%wavelength*1.e7 .eps. species2%wavelength*1.e7)) return
    end if
    
    equal = .true.

  end function fu_species_equal
  
  !==================================================================================

  logical function fu_modes_match(mode1,mode2) result(match)
    !Check if modes nominally match
    implicit none
    type(Taerosol_mode), intent(in) :: mode1,mode2

    match = .false.
    if ((mode1 == in_gas_phase .and. mode1 == in_gas_phase) .or. &
      & (mode1 == no_mode      .and. mode2 == no_mode)  ) then
      match = .true.
    !Reasonable number within 10nm... 
    elseif(mode1%nominal_d /= real_missing .and. &
          &  abs(mode1%nominal_d - mode2%nominal_d) < mode_diam_tolerance ) then
      match = .true.
    endif

  end function fu_modes_match

  !==================================================================================

  logical function fu_species_match(species1, species2) result(match)
    !Check if species nominally match
    implicit none
    type(silam_species), intent(in) :: species1, species2

    match = .false.
    if (defined(species1) .and. defined(species2)) then
      if (associated(species1%material, species2%material)) then
        match =  fu_modes_match(species1%mode, species2%mode)
      endif
    endif

  end function fu_species_match

  !==================================================================================

  logical function fu_species_defined(species)
    type(silam_species), intent(in) :: species

    fu_species_defined = fu_true(species%defined)  ! note change from silja_logical to logical

  end function fu_species_defined

  !==================================================================================


  subroutine addSpecies(speciesList, nSpecies,  moreSpecies, nMoreSpecies, ifAllNew_)
    !
    ! Add nMoreSpecies species from array moreSpecies into array
    ! speciesList, ensuring that no species are duplicated. On entry,
    ! the argument nspecies must contain the number of defined
    ! elements in speciesList. On exit, it holds the number of species
    ! defined after adding.
    !
    implicit none
    type(silam_species), dimension(:), intent(in) :: moreSpecies
    type(silam_species), dimension(:), pointer :: speciesList
    integer, intent(in) :: nMoreSpecies
    integer, intent(inout) :: nspecies
    logical, intent(in), optional :: ifAllNew_ ! Aall moreSpecies
                              !  must be new in speciesList 

    integer :: iTmp, jTmp
    logical :: ifAllNew

    ifAllNew = .false.  !Default
    if(present(ifAllNew_)) ifAllNew = ifAllNew_ 

    if(nSpecies == 0)then
      !
      ! No old species, copy the new ones to the list
      !
      if(.not.associated(speciesList))then
        call enlarge_species_list(speciesList, 0, nMoreSpecies)
        if(error)return
      else
        if(size(speciesList) < nMoreSpecies)then
          call enlarge_species_list(speciesList, size(speciesList), nMoreSpecies)
          if(error)return
      endif
      endif
      do iTmp = 1, nMoreSpecies
        speciesList(iTmp) = moreSpecies(iTmp)
      end do
      nspecies = nMoreSpecies
      return
    
    else
      !
      ! There are old species, have to check for absence of the duplicates
      !
      NEW_SPECIES : do iTmp = 1, nMoreSpecies
        do jTmp = 1, nspecies
          if (speciesList(jTmp) == moreSpecies(iTmp))then
              if(ifAllNew)then
                call msg_warning('Duplicate species faced, which is strictly forbidden', &
                               & 'addSpecies')
                call report(speciesList(jTmp))
                call set_error('Duplicate species faced, which is strictly forbidden', &
                             & 'addSpecies')
                return
              endif
              cycle NEW_SPECIES
          endif
        end do
        if (nSpecies == size(speciesList)) then
          call enlarge_species_list(speciesList, nSpecies, nSpecies+1) 
          if(error)return
        end if

        nspecies = nspecies+1
        speciesList(nspecies) = moreSpecies(iTmp)
  
      end do NEW_SPECIES
    endif
  
  end subroutine addSpecies


  !************************************************************************************

  function fu_species_to_short_string(species)result(chString)
    !
    ! Returns a string describing the species in a condenced form. 
    !
    implicit none

    ! Imported parameter
    type(silam_species), intent(in) :: species

    ! Return value
    character(len=clen) :: chString

    if(.not.defined(species))then
      chString = 'undefined'
      return
    else
      chString = ''
    endif
    !
    ! Species description consist of substance name, mode outlook, and wave length
    !
    if(associated(species%material))then
      chString = fu_name(species%material)
    endif
    if(defined(species%mode))then
      if(species%mode == in_gas_phase)then
        chString = chString + '_gas'
      elseif (.not. species%mode == no_mode) then
          chString = chString + '_m' + fu_aerosol_mode_size_to_str(species%mode%nominal_d) + &
                              & species%mode%name
      endif
    else
      ! Would be preferable to use an actual mode number. 
      call msg_warning('Undefined Mode of the species:', 'fu_species_to_short_string')
      call report(species)
    endif
    if(.not. (species%waveLength .eps. real_missing))then
      chString = chString + '_w' + fu_optical_wave_length_to_str(species%waveLength)
    endif
    
  end function fu_species_to_short_string


  !************************************************************************************
  
  function fu_species_from_short_string(chInString) result(species)
    !
    ! Decodes the three species parameters: substance, mode size and wave length from the 
    ! short string written by the above funciotn. The string looks like 
    ! <substance_name>_m<mode_value>_w<value>' at the end, which means aerosol mode mean 
    ! diameter in micrometres and optical wavelength in nanometres. In this case, the 
    ! pattern must resemble the real number with the decimal dot replaced with underscores. 
    !
    implicit none
    
    ! Imported parameter
    character(len=*), intent(in) :: chInString
    
    ! Return value
    type(silam_species) :: species
    
    ! Internal variables
    real :: fModeValue, fWaveLength
    integer :: iQ, iEnd
    character(len=clen) :: chInTmp, chTmp
    type(Taerosol_mode) :: aer_mode
    type(silam_material), pointer :: material
    
    species = species_missing
    !
    ! Stupidity check
    !
    if(len_trim(chInString) == 0 .or. fu_str_l_case(chInString) == 'undefined')return
    !
    ! Wave length is at the very end - find and cut it out.
    !
    fModeValue = real_missing
    fWaveLength = real_missing
    chInTmp = chInString

    iQ = index(chInTmp,'_w')
    if(iQ > 0)then              ! pattern found. There must be real value after it
      chTmp = chInTmp(iQ+2:)
      if(len_trim(chTmp) >0)then
        do iQ = 1, len_trim(chTmp)
          if(chTmp(iQ:iQ) == '_') chTmp(iQ:iQ) = '.'  ! prepare to reading
        end do
        read(unit=chTmp,fmt=*,iostat=iQ)fWaveLength
        if(iQ == 0)then
          chInTmp = chInTmp(1:index(chInTmp,'_w')-1)  ! cut the wave length part
          fWaveLength = fWaveLength * 1.e-9  ! to SI units [m]
        endif
      endif
    endif  ! '_w' pattern found
    !
    ! Having the wavelength cut out, check for aerosol mode
    !
    iQ = index(chInTmp,'_m')
    if(iQ > 0)then              
      !
      ! pattern found. There must be real value after it
      ! plus, possibly, name of the mode written immediately after the last digit.
      !
      chTmp = chInTmp(iQ+2:)
      if(len_trim(chTmp) >0)then
        iEnd = len_trim(chTmp)
        do iQ = 1, len_trim(chTmp)
          if(chTmp(iQ:iQ) == '_')then
            chTmp(iQ:iQ) = '.'  ! digital dot; prepare for reading
          elseif(.not. fu_if_digit(chTmp(iQ:iQ)))then
            iEnd = iQ-1
            exit
          endif
        end do
        read(unit=chTmp(1:iEnd), fmt=*, iostat=iQ) fModeValue
        if(iQ == 0)then    ! iQ is now read status
          chInTmp = chInTmp(1:index(chInTmp,'_m')-1)  ! cut the mode part from the main string
          fModeValue = fModeValue * 1.e-6  ! to SI units
          call set_aerosol_mode_from_d(aer_mode, fModeValue)
          call set_mode_name(aer_mode, chTmp(iEnd+1:len_trim(chTmp)))
        endif
      else
        call set_error('Strange mode: got _m but no diameter in:' + chInTmp,'fu_species_from_short_string')
        return
      endif
!call msg('Aerosol mode from fu_species_from_short_string:')
!call report(aer_mode)
    else
      !
      ! No _m => substance exists, no mode. Has to be gas
      ! For gaseous species the mode value goes as _gas
      !
      iQ = index(chInTmp,'_gas')
      if(iQ > 0)then              ! pattern found. There must be real value after it
        chInTmp = chInTmp(1:index(chInTmp,'_gas')-1)  ! cut the mode part
        aer_mode = in_gas_phase
      else
        aer_mode = no_mode
      endif
    endif  ! _m is present
    !
    ! What is left must be a material.
    ! All is then ready, can set the species
    !
    material => fu_get_material_ptr(chInTmp)
    if (associated(material)) then
       call set_species_from_params(species, &
                               & material, &
                               & aer_mode, &
                               & fWaveLength)
    else
       species = species_missing
    endif
    
  end function fu_species_from_short_string


  !*************************************************************************

  subroutine decode_id_params_from_io_str(chInputString, ifMultiLevel, iSilamQ, species, ifPush)
    !
    ! Searches the SILAM quantity from the given GrADS name of variable
    ! The main problem that GrADS has variable as quantity_short_string + '_' + substance,
    ! if any + '_m' + iAerosolMode, if any, and + '_w' + wave length, if any. Their separation 
    ! is non-trivial because it is easy to misread tempr_2m as quantity "tempr" and subst "2m". 
    ! Non-trivial combinations have to be checked before the simple ones to avoid misinterpreting.
    !
    implicit none

    ! Imported parameters
    character(len=*), intent(in) :: chInputString
    logical, intent(in) :: ifMultiLevel
    integer, intent(out) :: iSilamQ
    type(silam_species), intent(out) :: species
    logical, intent(in) :: ifPush    ! whether the quantity is really expected to be found

    ! Local variables
    integer, dimension(max_quantities) :: arQ
    integer :: iQ, iFound, iTmp, jTmp
    !
    ! Algorithm: check the whole currently active range of SILAM quantities and
    ! find out the suitable one from quantity_short_string. Shoudl several appear 
    ! to satisfy the comparison, the one with the longest name is taken.
    !

    arQ(1) = int_missing
    iFound = 1
    iSilamQ = int_missing

    !
    ! Have to check three ranges - multi-level, single-level and silam-own quantities
    ! Silam-own quantities can both require substance name and free from it. Be careful!
    !
    do iQ = first_silam_own_q, last_silam_own_q
      if(index(chInputString, trim(fu_quantity_short_string(iQ))) == 1)then
        arQ(iFound) = iQ
        iFound = iFound + 1
        arQ(iFound) = int_missing
      endif
    end do

    do iQ = first_multi_level_q, last_multi_level_q
      if(index(chInputString, trim(fu_quantity_short_string(iQ))) == 1)then
        arQ(iFound) = iQ
        iFound = iFound + 1
        arQ(iFound) = int_missing
      endif
    end do

    if(.not. ifMultiLevel)then    ! Do not check this if n_levs > 1
      do iQ = first_single_level_q, last_single_level_q
        if(index(chInputString, trim(fu_quantity_short_string(iQ))) == 1)then
          arQ(iFound) = iQ
          iFound = iFound + 1
          arQ(iFound) = int_missing
        endif
      end do
    endif
    iFound = iFound -1

    !
    ! All currently existing quantities have been checked. Let's now analyse the outcome
    !
    if(iFound == 0)then  
      !
      ! nothing is found
      !
      if(ifPush) call msg_warning('Failed to find the SILAM quantity for variable:' + chInputString, &
                                & 'decode_id_params_from_io_str')
      iSilamQ = int_missing
      species = species_missing

    else
      !
      ! One or more found -- take the longest one
      !
      iSilamQ = arQ(1)
      iTmp = len_trim(fu_quantity_short_string(iSilamQ))
      do iQ = 2, iFound
        jTmp = len_trim(fu_quantity_short_string(arQ(iQ))) ! New string length
        if( iTmp > jTmp ) continue  
        !Longer than old one
        iSilamQ = arQ(iQ)  
        iTmp = jTmp
      end do

      if (len_trim(chInputString) > iTmp + 1 ) then  !!!Some leftover in the name
         species = fu_species_from_short_string( chInputString(iTmp+2:))
         if( .not. defined(species))then
           if(ifPush) call set_error('Failed to find species for the string:' + chInputString + &
                             & ', resetting the quantity', 'decode_id_params_from_io_str')
           iSilamQ = int_missing
         endif
      else
          species = species_missing
      endif
    endif

  end subroutine decode_id_params_from_io_str


  !************************************************************************************

  subroutine enlarge_species_list(speciesLst, sizeOld, sizeNew)
    !
    ! Increases the size of the species list
    !
    implicit none 

    ! Imported parameters
    type(silam_species), dimension(:), pointer :: speciesLst
    integer, intent(in) :: sizeOld, sizeNew

    ! Local variables
    integer :: iStatus, i
    type(silam_species), dimension(:), pointer :: speciesTmp

    
    allocate(speciesTmp(sizeNew),stat=iStatus)
    if(iStatus /= 0)then
       call set_error('Failed to allocate memory for values.1','enlarge_species_list')
       return
     endif
     do i = 1, sizeOld
        speciesTmp(i) = speciesLst(i)
     end do
     do i = sizeOld+1, sizeNew
        speciesTmp(i) = species_missing
     end do

     if(sizeOld /= 0 .or. associated(speciesLst))then
        deallocate(speciesLst)
     endif 

     speciesLst => speciesTmp

  end subroutine enlarge_species_list


  !********************************************************************

  function fu_optical_wave_length_to_str(fWaveLen) result(strWaveLen)
    !
    ! Prints the mode size to the string and makes it short enough
    ! Rule: aerosol mode size goes in micrometres with max 2 digits after decimal point
    ! Total number of digits is not more than 3. Leading zero is omitted, dot is replaced
    ! with the underscore
    !
    implicit none

    ! Return value
    character(len=10) :: strWaveLen

    ! Imported parameter
    real, intent(in) :: fWaveLen  ! in m = SI unit

    ! Local variables
    integer :: iTmp

    iTmp = int(fWaveLen * 1.e9 + 0.5)  ! to nanometres
    if(iTmp < 1)then            ! Less than a nanometer
      write(unit=strWaveLen,fmt='(A1)')'0'
    elseif(iTmp < 10)then
      write(unit=strWaveLen,fmt='(I1)')iTmp
    elseif(iTmp < 100)then
      write(unit=strWaveLen,fmt='(I2)')iTmp
    elseif(iTmp < 1000)then
      write(unit=strWaveLen,fmt='(I3)')iTmp
    elseif(iTmp < 10000)then
      write(unit=strWaveLen,fmt='(I4)')iTmp
    else
      call msg('Too large wave length:', fWaveLen)
      call set_error('Too large wave length','fu_optical_wave_length_to_str')
    endif

  end function fu_optical_wave_length_to_str


  !************************************************************************************
  !
  ! Stuff for aerosol_modes
  !
  !************************************************************************************
  
  function fu_compare_modes_eq(mode1, mode2) result(equal)
    !
    ! Compare equality of aerosol modes, in practice, their max, min
    ! and mean diameters. Modes of different type are never considered
    ! equal, but no error is set.
    !
    type(Taerosol_mode), intent(in) :: mode1, mode2
    logical :: equal
    
    equal = .false.

    if (mode1%distr_type /= mode2%distr_type) return
    equal = (fu_compare_modes_general(mode1, mode2) == 0)
!    equal = (fu_compare_modes_mean_diam(mode1, mode2) == 0)

  end function fu_compare_modes_eq
  
  !==================================================================================

  logical function fu_mode_defined(mode) result(defined)
    implicit none
    type(Taerosol_mode), intent(in) :: mode
    defined = .not. (.not. mode%defined)    ! change from silja_logical to logical
  end function fu_mode_defined


  !==================================================================================
  
  function fu_set_mode_from_params(mode_type, fp1, fp2, fp3,  ip1, solubility, label) result(mode)
    !
    ! Create an aerosol mode of given type. Allowed options: 
    ! 1. fixed_diameter_flag, moving_diameter_flag -- fp1 and fp2 define
    ! limits. If known, volume mean diameter can be given as fp3.
    ! 2. lognormal_flag -- fp1 = median, fp2 = sigma
    ! 
    implicit none
    integer, intent(in) :: mode_type
    real, intent(in), optional :: fp1, fp2, fp3
    integer, intent(in), optional :: ip1, solubility
    character(len=*), optional :: label

    ! Return value of the function
    type(Taerosol_mode) :: mode

    ! Local variables
    real :: fTmp

    ! The creation procedure depends on the mode type
    !
    select case(mode_type)

      case(gas_phase_flag)
        mode = in_gas_phase

      case(no_mode_flag)
        mode = no_mode

      case (moving_diameter_flag, fixed_diameter_flag)
        if (.not. present(fp1) .or. .not. present(fp2)) then
          call set_error('Two parameters are needed for setting a size bin', 'fu_set_mode_from_params')
          return
        end if
        ! 
        ! fp1 = min_diam, fp2 = max_diam, fp3 = volume_mean_diameter
        ! if present, otherwise uniform distribution assumed.
        !
        if (fp1 < 0.0 .or. fp2 < 0.0 .or. fp2 < fp1) then
          call msg('MinD', fp1)
          call msg('MaxD', fp2)
          call set_error('Strange parameters for size bin', 'fu_set_mode_from_params')
          return
        end if
        mode%fp1 = fp1
        mode%fp2 = fp2

        if (present(fp3)) then
          if (fp3 > fp2 .or. fp3 < fp1 .or. fp3 < 0.0) then
            call msg('MinD', fp1)
            call msg('MaxD', fp2)
            call msg('MeanD', fp3)
            call set_error('Strange parameters for size bin', 'fu_set_mode_from_params')
            return
          end if
          mode%mass_mean_d = fp3
        else
          mode%mass_mean_d = ((fp1**2 + fp2**2)*(fp1+fp2) / 4.0) ** (1.0/3.0)
        end if

      case (lognormal_flag)
        ! 
        ! fp1 = (number) median diameter, fp2 = sigma. mode%fp3 set to
        ! volume mean diameter.
        !
        if (.not. present(fp1) .or. .not. present(fp2)) then
          call set_error('Two parameters are needed for setting a lognormal mode', &
                       & 'fu_set_mode_from_params')
          return
        end if

        if (fp1 < 0.0 .or. fp2 < 0.0) then
          call msg('Median:', fp1)
          call msg('Sigma', fp2)
          call set_error('Strange parameters for lognormal mode', 'fu_set_mode_from_params')
        end if

        mode%fp1 = fp1 !/ exp(log(fp2)**2 / 2.0)  ! nbr median
        mode%fp2 = fp2  ! std dev

        mode%mass_mean_d = mode%fp1 * exp(3. * log(fp2)**2) * exp(log(fp2)**2 / 2.0)
!      mode%mass_mean_d = fp1 * exp(log(fp2)**2 / 2.0)

      case(sea_salt_mode_flag)
        !
        ! Sea salt mode is set in the sea salt source term, from here we need only the range
        !
        if (.not. present(fp1) .or. .not. present(fp2)) then
          call set_error('Two parameters are needed for setting a sea salt mode', &
                       & 'fu_set_mode_from_params')
          return
        end if

        if (fp1 < 0.0 .or. fp2 < 0.0) then
          call msg('Min, max:', fp1, fp2)
          call set_error('Strange parameters for sea salt mode', 'fu_set_mode_from_params')
        end if

        mode%fp1 = fp1  ! min
        mode%fp2 = fp2  ! max

        mode%mass_mean_d = real_missing

      case default
        call msg('mode_type = ', mode_type)
        call set_error('Unsupported mode type', 'fu_set_mode_from_params')
    end select

    mode%distr_type = mode_type
    
    if (present(solubility)) then
      mode%solubility = solubility
    else
      mode%solubility = int_missing
    end if

    if (present(label)) then
      mode%name = label
    else
      mode%name = ''
    end if

    mode%nominal_d = mode%mass_mean_d

    mode%defined = silja_true

  end function fu_set_mode_from_params


  !*******************************************************************************

  function fu_set_mode_from_namelist(nlSetup) result(mode)
    !
    ! Create an aerosol mode of given type. Allowed options: 
    ! 1. fixed_diameter_flag, moving_diameter_flag -- fp1 and fp2 define
    ! limits. If known, volume mean diameter can be given as fp3.
    ! 2. lognormal_flag -- fp1 = median, fp2 = sigma
    ! 
    implicit none

    ! Imported parameter
    type(Tsilam_namelist), intent(in) :: nlSetup

    ! Return value
    type(Taerosol_mode) :: mode

    ! Local declarations
    type(silam_sp) :: spTmp
    real :: dens

    mode = aerosol_mode_missing
    spTmp%sp => fu_work_string()
    if(error)return
    !
    ! The mode definitions depend on the mode type
    !
    spTmp%sp = fu_str_u_case(fu_content(nlSetup,'mode_distribution_type'))
    if(spTmp%sp == '')then
      call set_error('No mode_distribution_type in the mode namelist','fu_set_mode_from_namelist')
      return
    endif
    !
    ! Start the selection itself
    !
    if(fu_content(nlSetup,'mode_density')== '')then
      dens = real_missing
    else
      dens = fu_set_named_value(fu_content(nlSetup,'mode_density'))
    endif
    if(spTmp%sp == 'GAS_PHASE')then
      mode = fu_set_mode_from_params(gas_phase_flag)

    elseif(spTmp%sp == 'NO_MODE')then
      mode = fu_set_mode_from_params(no_mode_flag)

    elseif(spTmp%sp == 'FIXED_DIAMETER')then
      if(len_trim(fu_content(nlSetup,'fix_diam_mode_mean_diameter')) > 0)then
        mode = fu_set_mode_from_params( &
                         & fixed_diameter_flag, &
                         & fu_set_named_value(fu_content(nlSetup,'fix_diam_mode_min_diameter')), &
                         & fu_set_named_value(fu_content(nlSetup,'fix_diam_mode_max_diameter')), &
                         & fu_set_named_value(fu_content(nlSetup,'fix_diam_mode_mean_diameter')), &
                         & solubility = fu_content_int(nlSetup,'mode_solubility'), &
                         & label=fu_content(nlSetup,'mode_name'))
      else
        mode = fu_set_mode_from_params( &
                         & fixed_diameter_flag, &
                         & fu_set_named_value(fu_content(nlSetup,'fix_diam_mode_min_diameter')), &
                         & fu_set_named_value(fu_content(nlSetup,'fix_diam_mode_max_diameter')), &
                         & solubility = fu_content_int(nlSetup,'mode_solubility'), &
                         & label = fu_content(nlSetup,'mode_name'))
      endif

    elseif(spTmp%sp == 'GAMMA_FUNCTION')then
      mode = fu_set_mode_from_params( &
                         & gamma_function_flag, &
                         & ip1 = fu_content_int(nlSetup,'gamma_mode_gamma_number'), &
                         & solubility = fu_content_int(nlSetup,'mode_solubility'), &
                         & label=fu_content(nlSetup,'mode_name'))

    elseif(spTmp%sp == 'MOVING_DIAMETER')then
      if(len_trim(fu_content(nlSetup,'mode_mean_diameter')) > 0)then
        mode = fu_set_mode_from_params( &
                         & moving_diameter_flag, &
                         & fu_set_named_value(fu_content(nlSetup,'moving_diam_mode_min_diameter')), &
                         & fu_set_named_value(fu_content(nlSetup,'moving_diam_mode_max_diameter')), &
                         & fu_set_named_value(fu_content(nlSetup,'moving_diam_mode_mean_diameter')), &
                         & solubility = fu_content_int(nlSetup,'mode_solubility'), &
                         & label=fu_content(nlSetup,'mode_name'))
      else
        mode = fu_set_mode_from_params( &
                         & moving_diameter_flag, &
                         & fu_set_named_value(fu_content(nlSetup,'moving_diam_mode_min_diameter')), &
                         & fu_set_named_value(fu_content(nlSetup,'moving_diam_mode_max_diameter')), &
                         & solubility = fu_content_int(nlSetup,'mode_solubility'), &
                         & label = fu_content(nlSetup,'mode_name'))
      endif

    elseif(spTmp%sp == 'LOGNORMAL')then
      mode = fu_set_mode_from_params( &
                         & lognormal_flag, &
                         & fu_set_named_value(fu_content(nlSetup,'lognormal_mode_mean_diameter')), &
                         & fu_set_named_value(fu_content(nlSetup,'lognormal_mode_sigma')), &
                         & solubility = fu_content_int(nlSetup,'mode_solubility'), &
                         & label = fu_content(nlSetup,'mode_name'))

    elseif(spTmp%sp == 'SEA_SALT')then
      mode = fu_set_mode_from_params( &
                         & sea_salt_mode_flag, &
                         & fu_set_named_value(fu_content(nlSetup,'sea_salt_mode_min_diameter')), &
                         & fu_set_named_value(fu_content(nlSetup,'sea_salt_mode_max_diameter')), &
                         & solubility = fu_content_int(nlSetup,'mode_solubility'), &
                         & label = fu_content(nlSetup,'mode_name'))
    else
      call set_error('Unknown distribution type:' + spTmp%sp,'fu_set_mode_from_namelist')
      return
    endif
    
    ! Set nominal diameter if available
    spTmp%sp = fu_content(nlSetup,'mode_nominal_diameter')
    if(spTmp%sp /= '')then
      mode%nominal_d =  fu_set_named_value(spTmp%sp)
    else
      mode%nominal_d = mode%mass_mean_d
    endif


    if(error)call set_error('Failed to set the mode','fu_set_mode_from_namelist')

    call free_work_array(spTmp%sp)

  end function fu_set_mode_from_namelist


  !==================================================================================

  subroutine set_mode_name(mode, mode_name)

    implicit none
    
    type(Taerosol_mode), intent(inout) :: mode
    character(len=*), intent(in) :: mode_name

    mode%name = mode_name
   
  end subroutine set_mode_name

  !==================================================================================

  function fu_mode_name(mode)

    implicit none
    
    type(Taerosol_mode), intent(in) :: mode
    character(len=substNmLen) :: fu_mode_name

    fu_mode_name = mode%name
   
  end function fu_mode_name

  !==================================================================================


  subroutine sort_modes(modes, indices, nmodes)
    !
    ! Sort an array of modes and an array of integers associated to
    ! them (in practice, their indices in a species list). Smallest
    ! goes first, gas phase is assumed smaller than any aerosol. same-diameter
    ! insoluble (or less soluble) is smaller than more soluble
    ! Only works on modes of type fixed diameter.
    !
    implicit none

    type(Taerosol_mode), dimension(:), intent(inout) :: modes
    integer, dimension(:), intent(inout) :: indices
    integer, intent(in) :: nModes

    type(Taerosol_mode) :: modeTmp
    integer :: iMode, indexTmp, iCompare
    logical :: ifChanged

    ifChanged = .true.
    do while(ifChanged)
      ifChanged = .false.
      do iMode = 1, nModes-1
#ifdef DEBUG
call msg('Mode i, i+1:')
call report(modes(iMode))
call report(modes(iMode+1))
#endif
        iCompare = fu_compare_modes_mean_diam(modes(iMode), modes(iMode+1))
        if(iCompare == 1)then
#ifdef DEBUG
call msg('Flip...')
#endif
          modeTmp = modes(iMode)
          modes(iMode) = modes(iMode+1)
          modes(iMode+1) = modeTmp
          indexTmp = indices(iMode)
          indices(iMode) = indices(iMode+1)
          indices(iMode+1) = indexTmp
          ifChanged = .true.
        elseif(iCompare == 0)then
          call msg('The following overlapping modes found')
          call report(modes(iMode))
          call report(modes(iMode+1))
          call set_error('Overlapping modes found','sort_modes')
          return
        endif
      end do
    end do



!    type(Taerosol_mode), dimension(maxSizeModes) :: modesTmp
!    integer, dimension(maxSizeModes) :: indicesTmp
!    integer :: i, j, inserted
!
!    modesTmp(1:nmodes) = aerosol_mode_missing
!    indicesTmp(1:nmodes) = int_missing
!    inserted = 0
!
!    OUTER: do i = 1, nModes
!      INNER: do j = 1, nModes
!
!        if (indicesTmp(j) == int_missing) then
!          modesTmp(j) = modes(i)
!          indicesTmp(j) = indices(i)
!          inserted = inserted + 1
!          cycle OUTER
!        end if
!
!        select case(fu_compare_modes_mean_diam(modes(i), modesTmp(j)))
!        case(0)
!          call set_error('Overlapping modes found','sort_modes')
!        case (-1)
!          modesTmp(j+1:j+inserted) = modesTmp(j:j+inserted-1)
!          modesTmp(j) = modes(i)
!          indicesTmp(j+1:j+inserted) = indicesTmp(j:j+inserted-1)
!          indicesTmp(j) = indices(i)
!          inserted = inserted + 1
!          cycle OUTER
!        case (1)
!          cycle INNER
!          
!        end select
!        
!        if (error) return
!
!      end do INNER
!    end do OUTER
!       
!    modes(1:nmodes) = modesTmp(1:nmodes)
!    indices(1:nmodes) = indicesTmp(1:nmodes)

  end subroutine sort_modes

  !==================================================================================

  real function fu_massmean_d_mode(mode) result(meanD)
    implicit none
    type(Taerosol_mode), intent(in) :: mode
    meanD = mode%mass_mean_d
  end function fu_massmean_d_mode

  !==================================================================================
  real function fu_nominal_d_mode(mode) result(meanD)
    implicit none
    type(Taerosol_mode), intent(in) :: mode
    meanD = mode%nominal_d
  end function fu_nominal_d_mode

  !==================================================================================
  real function fu_min_d_mode(mode) result(minD)
    implicit none
    type(Taerosol_mode), intent(in) :: mode

    if (mode%distr_type == fixed_diameter_flag .or. mode%distr_type == moving_diameter_flag) then
      minD = mode%fp1
    else
      call msg('Mode type = ', mode%distr_type)
      call set_error('min diameter not defined for this type of mode', 'fu_min_d_mode')
    end if
  end function fu_min_d_mode

  !================================================================================== 
  real function fu_max_d_mode(mode) result(maxD)
    implicit none
    type(Taerosol_mode), intent(in) :: mode
    
    if (mode%distr_type == fixed_diameter_flag .or. mode%distr_type == moving_diameter_flag) then
      maxD = mode%fp2
    else
      call msg('Mode type = ', mode%distr_type)
      call set_error('max diameter not defined for this type of mode', 'fu_max_d_mode')
    end if

  end function fu_max_d_mode

  !================================================================================== 
  integer function fu_solubility_mode(mode) result(solubility)
    implicit none
    type(Taerosol_mode), intent(in) :: mode
    solubility = mode%solubility
  end function fu_solubility_mode

  !==================================================================================
  subroutine set_massmean_d(mode, meanD, ifDynamic)
    implicit none
    type(Taerosol_mode), intent(inout) :: mode
    real, intent(in) :: meanD
    logical, intent(in) :: ifDynamic
    if(.not. (mode%defined == silja_true))then
      call set_error('Attempting to set mean diameter for undefined mode','set_massmean_d')
      return
    endif
    if (mode%distr_type /= fixed_diameter_flag) then
      call set_error('Attempting to change mode for non-bin type aerosol mode', &
                   & 'set_mean_d_mode')
      return
    end if
    mode%mass_mean_d = meanD
    mode%ifMeanDiamDynamic = ifDynamic
  end subroutine set_massmean_d

  !==================================================================================
  integer function fu_mode_type_mode(mode) result(mType)
    implicit none
    type(Taerosol_mode), intent(in) :: mode
    mType = mode%distr_type
  end function fu_mode_type_mode

  !==================================================================================
  logical function fu_ifMeanDiamDynamic(mode) result(ifDyn)
    implicit none
    type(Taerosol_mode), intent(in) :: mode
    ifDyn = mode%ifMeanDiamDynamic
  end function fu_ifMeanDiamDynamic




  !********************************************************************

  function fu_aerosol_mode_size_to_str(fAerosolModeSize) result(strModeSize)
    !
    ! Prints the mode size to the string and makes it short enough
    ! Rule: aerosol mode size goes in micrometres with max 2 digits after decimal point
    ! Total number of digits is not more than 3. Leading zero is omitted, dot is replaced
    ! with the underscore
    !
    implicit none

    ! Return value
    character(len=10) :: strModeSize

    ! Imported parameter
    real, intent(in) :: fAerosolModeSize  ! in m = SI unit

    ! Local variables
    integer :: iTmp
    real :: fTmp

    fTmp = fAerosolModeSize * 1.e6
    if(fTmp < 5.e-4)then            ! Less than half-a nanometer is actually a gas
      write(unit=strModeSize,fmt='(A1)')'0'
    elseif(fTmp < 0.01)then
      write(unit=strModeSize,fmt='(f4.3)')fTmp  ! ".00x"
    elseif(fTmp < 0.99)then
      write(unit=strModeSize,fmt='(f3.2)')fTmp  ! ".0x"
    elseif(fTmp < 10.)then
      write(unit=strModeSize,fmt='(f3.1)')fTmp  ! "x.x"
    elseif(fTmp < 99.5)then
      write(unit=strModeSize,fmt='(i2)')int(fTmp+0.5)  ! "xx"
    elseif(fTmp < 999.5)then
      write(unit=strModeSize,fmt='(i3)')int(fTmp+0.5)  ! "xxx"
    else
      call msg('Too large mode size:', fAerosolModeSize)
      call set_error('Too large mode size','fu_aerosol_mode_size_to_str')
    endif

    do iTmp = 1, len_trim(strModeSize)
      if(strModeSize(iTmp:iTmp) == '.')strModeSize(iTmp:iTmp) = '_'
    end do

  end function fu_aerosol_mode_size_to_str



  !**********************************************************************************
  !
  ! Encapsulation of the chemical_mapper and chemical_adaptor types
  !
  !**********************************************************************************

  subroutine create_adaptor(species_in, species_out, adaptor, allow_missing)
    ! 
    ! Allocate and define a adaptor between to species lists. After
    ! created, it holds that species_in(i) = species_out(adaptor%isp(i)).
    ! It is an error if species_out does not contain all of species_in. 
    ! 
    implicit none
    type(silam_species), dimension(:), intent(in) :: species_in
    type(silam_species), dimension(:), intent(in) :: species_out
    type(chemical_adaptor), intent(out) :: adaptor
    logical, intent(in), optional :: allow_missing

    integer :: i, j, k, stat
    logical :: allow_missing_loc

    ! For each species in species_in, find corresponding species_out.

    adaptor%defined = silja_false

    if (present(allow_missing)) then
      allow_missing_loc = allow_missing
    else
      allow_missing_loc = .false.
    end if
    
    allocate(adaptor%iSp(size(species_in)), stat = i)
    if(fu_fails(i==0,'Failed to allocate the chemical adaptor indices','create_adaptor'))return
    
    do i = 1, size(species_in)
      if(.not. defined(species_in(i)))exit   ! all done, exit
      adaptor%nSpecies = i
      j = fu_index(species_in(i), species_out, size(species_out), .true.) !! by name
      if (j == int_missing .and. .not. allow_missing_loc) then
        call msg('Species not found in list: '+fu_str(species_in(i)))
        call msg('The list includes:')
        do k=1, size(species_out)
          call msg("species_out("+fu_str(k)+")="+fu_str(species_out(k)))
        end do
        call set_error('Species not found', 'create_adaptor')
        return
      end if
      adaptor%iSp(i) = j
    end do
    
    adaptor%defined = silja_true
    
  end subroutine create_adaptor

  
  !************************************************************************************
  
  function chemical_adaptor_missing()result(adaptor_missing)
    type (chemical_adaptor) :: adaptor_missing
!    nullify(adaptor_missing%iSp)
    adaptor_missing%nSpecies =  int_missing
    adaptor_missing%defined = silja_false
  end function chemical_adaptor_missing
  
  
  !************************************************************************************
  
  subroutine destroy_adaptor(adaptor)
    type (chemical_adaptor), intent(inout) :: adaptor
    if(fu_true(adaptor%defined)) deallocate(adaptor%iSp)
    adaptor%defined = silja_false
  end subroutine destroy_adaptor
  
  
  !************************************************************************************
  
  subroutine copy_adaptor(adaptor_out, adaptor_in)
    !
    ! Copies one adaptor to another one
    !
    implicit none

    ! Imported parameters
    type(chemical_adaptor), intent(in) :: adaptor_in
    type(chemical_adaptor), intent(out) :: adaptor_out

    if(defined(adaptor_in))then
      allocate(adaptor_out%iSp(size(adaptor_in%iSp)))
      adaptor_out%iSp(:) = adaptor_in%iSp(:)
      adaptor_out%nSpecies = adaptor_in%nSpecies
      adaptor_out%defined = silja_true
    else
      call destroy_adaptor(adaptor_out)
    endif
  end subroutine copy_adaptor
    

  !************************************************************************************

  subroutine set_speciesReferences(species_in, species_out, speciesReferences, &
                                 & ignoreMode, ignoreSubst, allocateFractions)
    
    implicit none

    type(silam_species), dimension(:), intent(in) :: species_in, species_out
    type(TspeciesReference), dimension(:) :: speciesReferences
    logical, intent(in), optional :: ignoreMode, ignoreSubst, allocateFractions

    integer :: i, j, nrefs, allocstat
    integer, dimension(:), pointer :: indices

    type(Taerosol_mode) :: mode_in, mode_out
    real :: wave_in, wave_out
    logical :: ifMode, ifSubst

    if(present(ignoreMode))then
      ifMode = .not. ignoreMode
    else
      ifMode = .true.
    endif
    if(present(ignoreSubst))then
      ifSubst = .not. ignoreSubst
    else
      ifSubst = .true.
    endif

    indices => fu_work_int_array()

    do i = 1, size(species_in)
      nrefs = 0

      do j = 1, size(species_out)
        if (ifSubst .and. (fu_substance_name(species_in(i)) /= fu_substance_name(species_out(j)))) cycle
        if(ifMode)then
          mode_in = fu_mode(species_in(i))
          mode_out = fu_mode(species_out(j))
          if (defined(mode_in) .and. defined(mode_out) .and. .not.(mode_in == mode_out)) cycle
        endif
        wave_in = fu_optical_wave_length(species_in(i))
        wave_out = fu_optical_wave_length(species_out(j))
        if (.not.(wave_in .eps. real_missing) .and. .not.(wave_out .eps. real_missing) &
          & .and. .not.(wave_in*1e9 .eps. wave_out*1e9)) cycle
        nrefs = nrefs + 1
        indices(nrefs) = j
        
      end do
        
      if (nrefs > 0) then
        allocate(speciesReferences(i)%indSpeciesTo(nrefs), stat=allocstat)
        if (allocstat /= 0) then
          call set_error('Allocate failed', 'set_speciesReferences')
          return
        end if
        if (present(allocateFractions)) then
          if (allocateFractions) then
            allocate(speciesReferences(i)%fract(nrefs), stat=allocstat)
            if (allocstat /= 0) then
              call set_error('Allocate failed', 'set_speciesReferences')
              return
            end if
            speciesReferences(i)%fract = 1.0
          end if
        else
          nullify(speciesReferences(i)%fract)
        end if

        speciesReferences(i)%indSpeciesTo(1:nrefs) = indices(1:nrefs)
        speciesReferences(i)%nrefSpecies = nrefs
      else
        nullify(speciesReferences(i)%indSpeciesTo)
        speciesReferences(i)%nrefSpecies = 0
      end if

    end do
    
    call free_work_array(indices)


  end subroutine set_speciesReferences

  !==================================================================================

  subroutine kill_speciesReferences(speciesReferences)
    implicit none

    type(TspeciesReference), dimension(:) :: speciesReferences
    integer :: i
 
    do i = 1, size(speciesReferences)
      if(associated(speciesReferences(i)%indSpeciesTo))deallocate(speciesReferences(i)%indSpeciesTo)
      if(associated(speciesReferences(i)%fract))deallocate(speciesReferences(i)%fract)
    enddo
 
  end subroutine kill_speciesReferences

  !==================================================================================
  function fu_adaptor_defined(adaptor) result(defined)
    implicit none
    type(chemical_adaptor), intent(in) :: adaptor
    logical :: defined
    
    defined = .not. (.not. adaptor%defined)

  end function fu_adaptor_defined



  !************************************************************************************

  subroutine create_chemical_mapper(species_list, mapper, nspecies_)
    !
    ! Create and allocate a chemical mapper (see general description
    ! on the top of the module). After completed, the species of
    ! species_list can be indexed by substance, mode and wavelength
    ! using the functions below. Optionally, it is possible to map
    ! only the first nspecies_ species in the list (the mapper should
    ! then work with no problems for the subset).
    !
    implicit none
    type(silam_species), dimension(:), intent(in) :: species_list
    integer, optional, intent(in) :: nspecies_
    type(chemical_mapper), intent(out) :: mapper

    character(len=substNmLen), dimension(:), pointer :: names
    character(len=substNmLen) :: name
    integer :: stat, iTmp,  jTmp, kTmp, nsubstances, nmodes, nspecies, nwaves
    integer, dimension(:), pointer :: indices
    type(Taerosol_mode), dimension(:), allocatable :: modelist
    real, dimension(:), pointer :: waves
    type(silam_species), dimension(:), allocatable :: species_tmp
   
    waves => fu_work_array()
    indices => fu_work_int_array()

    if (.not. present(nspecies_)) then
      nspecies = size(species_list)
    else
      nspecies = nspecies_
    end if

    allocate(names(nspecies), modelist(nspecies), stat=stat)
    if (stat /= 0) then
      call set_error('Allocate failed', 'create_chemical_mapper')
      return
    end if
   
    call substances(species_list(1:nspecies), names, nsubstances)

    allocate(mapper%names(nsubstances), mapper%nmodes(nsubstances), mapper%isp(nsubstances),&
           & mapper%iSubst(nspecies), stat=stat)
    if (stat /= 0) then
      call set_error('Allocate failed', 'create_chemical_mapper')
      return
    end if

    mapper%nSubstances = nsubstances    
    mapper%nSpecies = nSpecies
    
    allocate(mapper%nwaves(nsubstances), stat=stat)
    if (stat /= 0) then
      call set_error('Allocate failed', 'create_chemical_mapper')
      return
    end if
      
    ! In the followins, we 1) get the modes for substances (they will
    ! be sorted by size) and 2) for each mode, get the
    ! wavelengths. The elements of iSp are allocated accordingly.

    do iTmp = 1, nsubstances
      mapper%names(iTmp) = names(iTmp)
      call modes_for_subst(species_list(1:nspecies), names(iTmp), modelist, indices, nmodes)
      mapper%nmodes(iTmp) = nmodes
      
      if (stat /= 0) then
        call set_error('Allocate failed', 'create_chemical_mapper')
        return
      end if
            
      allocate(mapper%nwaves(iTmp)%ptr1d(nmodes), mapper%iSp(iTmp)%ptr2d(nmodes), stat=stat)
      if (stat /= 0) then
        call set_error('Allocate failed', 'create_chemical_mapper')
        return
      end if
      nullify(mapper%nwaves(iTmp)%ptr2d)
      
      do jTmp = 1, nmodes
        call waves_for_subst_mode(species_list(1:nspecies), names(iTmp), modelist(jTmp),&
                                & waves, indices, nwaves)
        allocate(mapper%iSp(iTmp)%ptr2d(jTmp)%ptr1d(nwaves), stat=stat)
        if (stat /= 0) then
          call set_error('Allocate failed', 'create_chemical_mapper')
          return
        end if
        nullify(mapper%nwaves(iTmp)%ptr2d)
        mapper%nwaves(iTmp)%ptr1d(jTmp) = nwaves
        do kTmp = 1, nwaves
          mapper%iSp(iTmp)%ptr2d(jTmp)%ptr1d(kTmp) = indices(kTmp)
          mapper%iSubst(indices(kTmp)) = iTmp
        end do
      end do
         
    end do

    mapper%defined = silja_true
    
    call free_work_array(waves)
    call free_work_array(indices)

    contains 

      !----------------------------------------------------------------------------

      subroutine substances(species, names, nsubstances)
        type(silam_species), dimension(:), intent(in) :: species
        character(len=substNmLen), dimension(:), intent(out) :: names
        integer, intent(out) :: nsubstances
    
        integer :: iSp, nspecies
        
        nspecies = size(species)
        
        nsubstances = 0
        names(1:nspecies) = char_missing
        do iSp = 1, nspecies
          if (.not. defined(species_list(iSp)) ) cycle
          if (fu_index(fu_name(species_list(iSp)%material), names) < 0) then
            nsubstances = nsubstances + 1
            names(nsubstances) = fu_name(species_list(iSp)%material)
          end if
        end do
        
      end subroutine substances

      !----------------------------------------------------------------------------

      subroutine modes_for_subst(species, subst, modes, indices, nmodes)
        type(silam_species), dimension(:), intent(in) :: species
        character(len=substNmLen), intent(in) :: subst

        type(Taerosol_mode), dimension(:), intent(out) :: modes
        integer, dimension(:), intent(inout) :: indices
        integer, intent(out) :: nmodes

        integer :: iSp2, iMode
        
        !modes(1:size(species)) = aerosol_mode_missing
        nmodes = 0

        SPECIES_LOOP : do iSp2 = 1, size(species)
          if (fu_name(species(iSp2)%material) == subst) then
            ! check if the mode already exists (for a different wave)
            if (.not. defined(species(iSp2)) ) cycle
            do iMode = 1, nmodes
              if (modes(iMode) == species(iSp2)%mode) cycle SPECIES_LOOP
            end do
            nmodes = nmodes+1
            modes(nmodes) = species(iSp2)%mode
            indices(nmodes) = iSp2
          end if
        end do SPECIES_LOOP
        
        if (nmodes > 1) call sort_modes(modes, indices, nmodes)
        if(error)call msg('Mode sorting failed for: ' +subst)

      end subroutine modes_for_subst

      !----------------------------------------------------------------------------

      subroutine waves_for_subst_mode(species, subst, mode, waves, indices, nwaves)
        type(silam_species), dimension(:), intent(in) :: species
        character(len=substNmLen), intent(in) :: subst
        type(Taerosol_mode), intent(in) :: mode

        real, dimension(:), intent(out) :: waves
        integer, dimension(:), intent(out) :: indices
        integer, intent(out) :: nwaves

        integer :: iSp3, iWave
        real, dimension(:), pointer :: wavesTmp
        integer, dimension(:), pointer :: arTmp
        logical :: temp
        
        wavesTmp => fu_work_array()
        arTmp => fu_work_int_array()
        
        nwaves = 0
        !waves(1:size(species)) = real_missing

        do iSp3 = 1, size(species)
          if (.not. defined(species(iSp3)) ) cycle
          if (fu_name(species(iSp3)%material) == subst .and. species(iSp3)%mode == mode) then
            nwaves = nwaves + 1
            waves(nwaves) = species(iSp3)%wavelength
            arTmp(nwaves) = iSp3
          end if
        end do
        
        wavesTmp(1:nwaves) = waves(1:nwaves)

        if (nwaves > 1) then 
          temp = fu_sort_real_array(waves(1:nwaves), ascending, .false.)
          do iSp3 = 1, nwaves
            do iWave = 1, nwaves
              if (waves(iSp3) .eps. wavesTmp(iWave)) indices(iSp3) = arTmp(iWave)
            end do
          end do
        else
          indices(1) = arTmp(1)
        end if

        call free_work_array(arTmp)
        call free_work_array(wavesTmp)

      end subroutine waves_for_subst_mode

  end subroutine create_chemical_mapper


  !************************************************************************************

  subroutine copy_chemical_mapper(mapper1, mapper2)
    !
    ! Allocate and define mapper 2 as copy of mapper1.
    !
    implicit none
    type(chemical_mapper), intent(in) :: mapper1
    type(chemical_mapper), intent(out) :: mapper2

    integer :: i, j, stat, nSubst, nSp
    
    mapper2%defined = silja_false

    if (.not. mapper1%defined) then
      call set_error('Undefined mapper', 'copy_chemical_mapper')
    end if

    nSubst = mapper1%nSubstances
    nSp = mapper1%nSpecies

    allocate(mapper2%names(nSubst), mapper2%iSubst(nSp), mapper2%nModes(nSubst),&
           & mapper2%iSp(nSubst), mapper2%nWaves(nSubst), stat=stat)
    if (stat /= 0) then
      call set_error('Allocate failed', 'copy_chemical_mapper')
      return
    end if

    mapper2%names(1:nsubst) = mapper1%names(1:nsubst)
    mapper2%iSubst(1:nSp) = mapper1%iSubst(1:nSp)
    mapper2%nModes(1:nSubst) = mapper1%nModes(1:nSubst)

    mapper2%nSubstances = mapper1%nSubstances
    mapper2%nSpecies = mapper1%nSpecies
    
    
    do i = 1, nSubst
      allocate(mapper2%iSp(i)%ptr2d(mapper1%nModes(i)), &
             & mapper2%nWaves(i)%ptr1d(mapper1%nModes(i)), stat=stat)
      if (stat /= 0) then
        call set_error('Allocate failed', 'copy_chemical_mapper')
        return
      end if
      nullify(mapper2%nwaves(i)%ptr2d)
      mapper2%nWaves(i)%ptr1d(:) = mapper1%nWaves(i)%ptr1d(:)
      do j = 1, mapper1%nModes(i)
        allocate(mapper2%iSp(i)%ptr2d(j)%ptr1d(mapper1%nWaves(i)%ptr1d(j)),&
               & stat=stat)
        if (stat /= 0) then
          call set_error('Allocate failed', 'copy_chemical_mapper')
          return
        end if
        mapper2%iSp(i)%ptr2d(j)%ptr1d(:) = mapper1%iSp(i)%ptr2d(j)%ptr1d(:)
      end do
    end do
    
    mapper2%defined = silja_true

  end subroutine copy_chemical_mapper


  !************************************************************************************

  subroutine free_chemical_mapper(mapper)
    ! 
    ! Release the mapper structure.
    implicit none
    type(chemical_mapper), intent(inout) :: mapper

    integer :: i, j

    do i = 1, fu_nSubst(mapper)
      do j = 1, fu_nModes(i, mapper)
        deallocate(mapper%iSp(i)%ptr2d(j)%ptr1d)
      end do
      deallocate(mapper%iSp(i)%ptr2d, mapper%nWaves(i)%ptr1d)
    end do

    deallocate(mapper%names, mapper%iSp, mapper%iSubst, mapper%nModes, mapper%nWaves)

  end subroutine free_chemical_mapper

  !************************************************************************************

  subroutine set_chem_mapper_missing(mapper)
    implicit none
    type(chemical_mapper), intent(out) :: mapper

    nullify(mapper%names, mapper%isp, mapper%iSubst, mapper%nmodes, mapper%nwaves)
    mapper%nsubstances = int_missing
    mapper%nspecies = int_missing
    mapper%defined = silja_false

  end subroutine set_chem_mapper_missing
    

  !************************************************************************************

  function fu_iSubst_name(substNm, mapper) result(isubst)
    !
    ! Find the index of a substance name mapped by the mapper.
    !
    implicit none
    character(len=*), intent(in) :: substNm
    type(chemical_mapper), intent(in) :: mapper
    
    integer :: isubst

    isubst = fu_index(substNm, mapper%names)
    if (isubst < 1) then
      call set_error('No such substance: ' + substNm,'fu_iSubst_name')
      return
    end if

  end function fu_iSubst_name

  !************************************************************************************

  function fu_iSubst_ind(iSpecies, mapper) result(iSubst)
    !
    ! Find the substance index in the mapper of species iSpecies.
    !
    implicit none
    integer, intent(in) :: iSpecies
    type(chemical_mapper), intent(in) :: mapper
    
    integer :: iSubst

    iSubst = mapper%iSubst(iSpecies)
    
  end function fu_iSubst_ind


  !************************************************************************************

  function fu_nwaves_mapper(iSubst, iMode, mapper) result(nwaves)
    ! Return the number of wavelengths for given substance and
    ! mode. At least one is always returned, even if the wavelength is
    ! not defined.
    ! 
    implicit none
    type(chemical_mapper), intent(in) :: mapper
    integer, intent(in) :: iSubst, iMode

    integer :: nwaves

    nwaves = mapper%nwaves(isubst)%ptr1d(imode)
  end function fu_nwaves_mapper

  !************************************************************************************

  function fu_nModes_mapper(iSubst, mapper) result(nmodes)
    ! Return the number of modes for a given substance.
    !
    implicit none
    type(chemical_mapper), intent(in) :: mapper
    integer, intent(in) :: iSubst

    integer :: nmodes
    
    nmodes = mapper%nmodes(isubst)
  end function fu_nModes_mapper

  !************************************************************************************

  function fu_nSubst(mapper) result(nsubst)
    ! Return the number of substances mapped by the mapper.
    !
    implicit none
    type(chemical_mapper), intent(in) :: mapper
    
    integer :: nSubst
    
    if (.not. mapper%defined) then
      call set_error('Undefined mapper', 'fu_nsubst')
      return
    end if

    nSubst = mapper%nSubstances
  end function fu_nSubst

  !************************************************************************************

  function fu_substance_name_mapper(iSubst, mapper) result(name)
    !
    ! Return the name of substance isubst. 
    ! Same as fu_name(species_list(fu_iSp(iSubst,1,1),mapper)).
    !
    implicit none
    integer, intent(in) :: iSubst
    type(chemical_mapper), intent(in) :: mapper
    character(len=substNmLen) :: name

    name = mapper%names(iSubst)
  end function fu_substance_name_mapper

  !************************************************************************************

  function fu_iSp_wl(iSubst, iMode, iWave, mapper) result(isp)
    !
    ! Using a chemical_mapper, return the 1D species index corresponding
    ! to given substance, mode and wavelength indices.
    !
    implicit none
    integer, intent(in) :: isubst, iMode, iWave
    type(chemical_mapper), intent(in) :: mapper

    integer :: isp

    if (.not. mapper%defined) then
      call set_error('Undefined subst-mode-mapper','fu_iSp')
      return
    end if
    
    isp = mapper%isp(isubst)%ptr2d(imode)%ptr1d(iWave)

  end function fu_iSp_wl

  !************************************************************************************

  function fu_iSp_short(iSubst, iMode, mapper) result(isp)
    implicit none
    integer, intent(in) :: isubst, iMode
    type(chemical_mapper), intent(in) :: mapper

    integer :: isp

    if (.not. mapper%defined) then
      call set_error('Undefined subst-mode-mapper','fu_iSp')
      return
    end if
    
    isp = mapper%isp(isubst)%ptr2d(imode)%ptr1d(1)

  end function fu_iSp_short


  !************************************************************************************

  integer function fu_iSp_from_values(species_list, nSpecies, chSubstNm, aerMode, fWaveLength)
    !
    ! Finds the given species in the given list
    !
    implicit none

    ! Imported parameters
    type(silam_species), dimension(:), pointer :: species_list
    integer, intent(in) :: nSpecies
    character(len=*), intent(in) :: chSubstNm
    type(TAerosol_mode), intent(in) :: aerMode
    real, intent(in) :: fWaveLength

    ! Local variables
    integer :: iSp

    do iSp = 1, nSpecies
      if(fu_substance_name(species_list(iSp)) == chSubstNm)then
        !
        ! Substance name must exist and now it is found
        !
        if(defined(aerMode))then
          !
          ! Defined mode - check for equivalence
          !
          if(fu_mode(species_list(iSp)) == aerMode)then
            if(fWaveLength .eps. real_missing)then   ! if undefined wave length then nothing to check
              fu_iSp_from_values = iSp
              return
            else
              !
              ! Feasible wave length - check it too
              !
              if(fu_optical_wave_length(species_list(iSp)) .eps. fWaveLength)then
                fu_iSp_from_values = iSp
                return
              endif
            endif
          endif
        else
          !
          ! Mode is undefined but wave length can still be both defined and not
          !
          if(fWaveLength .eps. real_missing)then
            fu_iSp_from_values = iSp
            return
          else
            if(fu_optical_wave_length(species_list(iSp)) .eps. fWaveLength)then
              fu_iSp_from_values = iSp
              return
            endif
          endif
        endif
      endif
    end do

    fu_iSp_from_values = int_missing

  end function fu_iSp_from_values


  !***********************************************************************************
  
  integer function fu_iSp_from_species(species_list, nSpecies, speciesIn)
    !
    ! Finds the given species in the given list
    !
    implicit none

    ! Imported parameters
    type(silam_species), dimension(:), pointer :: species_list
    integer, intent(in) :: nSpecies
    type(silam_species), intent(in) :: speciesIn

    ! Local variables
    integer :: iSp

    do iSp = 1, nSpecies
      if(species_list(iSp) == speciesIn)then
        fu_iSp_from_species = iSp
        return
      endif
    end do

    fu_iSp_from_species = int_missing

  end function fu_iSp_from_species

  !**********************************************************************************
  !
  ! Encapsulation of the cocktail descriptors
  !
  !**********************************************************************************

  function fu_descr_defined(descriptor) result(defined)
    implicit none
    type(Tcocktail_descr), intent(in) :: descriptor
    logical :: defined

    defined = (descriptor%defined == silja_true) .and. &
            & (descriptor%nSpecies > 0) .and. &
            & (descriptor%nSpecies < 1000)
    
  end function fu_descr_defined


  !**************************************************************

  subroutine get_inventory_descr(descriptor, species, nspecies)
    !
    ! Return a pointer to the list of species included by this desriptor.
    !
    implicit none
    type(Tcocktail_descr), intent(in) :: descriptor
    type(silam_species), dimension(:), pointer :: species
    integer, intent(out) :: nspecies
    
    if (.not. defined(descriptor)) then
      nullify(species)
      nspecies = 0
      return
    end if
    
    species => descriptor%species
    nspecies = descriptor%nspecies
    
  end subroutine get_inventory_descr

  !*****************************************************************************

  subroutine cocktail_descriptor_report(descr)
    !
    ! Prints the report of the cocktail descriptor
    !
    implicit none

    ! Imported parameters
    type(Tcocktail_descr), intent(in) :: descr

    ! Local variables
    integer :: iSpecies

    call msg('---------- cocktail descriptor report -------------')
    if(descr%defined == silja_true)then
      call msg("cocktail_name = "//trim(descr%cocktail_name))
      do iSpecies = 1, descr%nSpecies
        call msg(fu_name(descr%species(iSpecies)%material), descr%mass_fractions(iSpecies))
      end do
    else
      call msg('Undefined descriptor')
    endif  ! if the descriptor defined
    call msg('---------------------------------------------------')

  end subroutine cocktail_descriptor_report


  !*********************************************************************

  subroutine link_emis_descr_to_species(species_list, descriptor)
    !
    ! Create the indexing from descriptor species to species_list
    ! (normally emission species). Uses adaptors, which includes
    ! checking that descriptor species are contained by species_list.
    ! 
    implicit none
    type(silam_species), dimension(:), intent(in) :: species_list
    type(Tcocktail_descr), intent(inout) :: descriptor

    type(chemical_adaptor) :: adaptor
    integer :: i

    call create_adaptor(descriptor%species, species_list, adaptor)
    if (error) return
    
    do i = 1, descriptor%nSpecies
      descriptor%iEmisCocktSpeciesMapping(i) = adaptor%iSp(i)
    end do

    
  end subroutine link_emis_descr_to_species


  !***********************************************************************************

  subroutine collect_descriptor_emis_rates(descr, rate_descr_unit, arRate)
    !
    ! Returns the arRate of rates of emission species - in the order of the emission mass map
    ! and in the unit of the emission mass map, i.e., with the unit dictated by the species
    ! The idea is that the descriptor upon definition is explored into the species list - see 
    ! link_emis_descr_to_species. Here we use this link to explore the rates
    !
    ! ATTENTION.
    ! The sub is called many times, so no checking whatsoever
    !
    implicit none

    ! Imported parameters
    type(Tcocktail_descr), intent(in) :: descr
    real, intent(in) :: rate_descr_unit
    real, dimension(:), pointer :: arRate

    ! Local variables
    integer :: iSp
    real :: factor

    do iSp = 1, descr%nSpecies
      arRate(descr%iEmisCocktSpeciesMapping(iSp)) = arRate(descr%iEmisCocktSpeciesMapping(iSp)) + &
                                                  & rate_descr_unit * &
                                                           & descr%f_descr_unit_2_species_unit(iSp)
    end do

  end subroutine collect_descriptor_emis_rates


  !***************************************************************************

  subroutine link_src_species_to_given_list(nDescriptors, &
                                          & nSpeciesInDescr, &
                                          & cocktail_descr_lst, &
                                          & fDescr2SpeciesUnit, &
                                          & pEmisSpeciesMapping, &
                                          & species_list)
    ! Imported parameters
    integer, intent(in) :: nDescriptors
    integer, dimension(:), pointer :: nSpeciesInDescr
    type(Tcocktail_descr), dimension(:), pointer :: cocktail_descr_lst  ! (nCocktails)
    real, dimension(:,:), pointer :: fDescr2SpeciesUnit
    integer, dimension(:,:), pointer :: pEmisSpeciesMapping
    type(silam_species), dimension(:), pointer :: species_list

    ! Local variables
    integer :: iTmp, iDescr, iSpecies
    integer, dimension(:), pointer :: pMapping
    real, dimension(:), pointer :: factorTmp, fractions

    !
    ! To start with, allocate the space
    !
    iTmp = fu_nSpecies(cocktail_descr_lst(1))
    do iDescr = 2, nDescriptors
      if(iTmp < fu_nSpecies(cocktail_descr_lst(iDescr))) &
                                          & iTmp = fu_nSpecies(cocktail_descr_lst(iDescr))
    end do
    allocate(nSpeciesInDescr(nDescriptors), &
           & fDescr2SpeciesUnit(iTmp, nDescriptors), &
           & pEmisSpeciesMapping(iTmp, nDescriptors), stat = iDescr)
    if(fu_fails(iDescr == 0, 'Failed to allocate the space for descriptor mapping', &
                           & 'link_src_species_to_given_list'))return
    !
    ! Now go descriptors one by one making the link
    !
    do iDescr = 1, nDescriptors
      call link_emis_descr_to_species(species_list, cocktail_descr_lst(iDescr))
      if(error)return
      nSpeciesInDescr(iDescr) = fu_nSpecies(cocktail_descr_lst(iDescr))
      pMapping => fu_emis_cocktail_mapping(cocktail_descr_lst(iDescr))
      factorTmp => fu_emis_cocktail_scaling(cocktail_descr_lst(iDescr))
      fractions => fu_emis_fractions(cocktail_descr_lst(iDescr))
      do iSpecies = 1, nSpeciesInDescr(iDescr)
        pEmisSpeciesMapping(iSpecies,iDescr) = pMapping(iSpecies)
        fDescr2SpeciesUnit(iSpecies,iDescr) = factorTmp(iSpecies) * fractions(iSpecies)
      end do
    end do
  
  end subroutine link_src_species_to_given_list


  !************************************************************************************

  function fu_nSpecies_descr(descriptor) result(nSpecies)
    implicit none
    type(Tcocktail_descr), intent(in) :: descriptor
    integer :: nSpecies
    
    nSpecies = descriptor%nSpecies

  end function fu_nSpecies_descr


  !************************************************************************************

  function fu_emis_cocktail_mapping(descriptor) result(indices)
    implicit none
    type(Tcocktail_descr), intent(in) :: descriptor
    
    integer, dimension(:), pointer :: indices

    indices => descriptor%iEmisCocktSpeciesMapping
    
  end function fu_emis_cocktail_mapping


  !************************************************************************************

  function fu_emis_cocktail_scaling(descriptor) result(factors)
    implicit none
    type(Tcocktail_descr), intent(in) :: descriptor
    
    real, dimension(:), pointer :: factors

    factors => descriptor%f_descr_unit_2_species_unit
    
  end function fu_emis_cocktail_scaling


  !************************************************************************************

  subroutine set_cocktail_description(chCocktailNm, cocktail_descr, ifSpecies)
    !
    ! Initialises all types of cocktails. Actually, just reads the 
    ! substance and material data from the files to memory.
    ! Corresponding database names have to be already loaded from the
    ! setup file.
    !
    implicit none

    ! Imported parameters 
    character(len=*), intent(in) :: chCocktailNm
    type(Tcocktail_descr), intent(out) :: cocktail_descr
    logical, intent(out) :: ifSpecies    ! if there is no cocktail but chemical exists, 
                                         ! try to use that one
    ! Local variables
    integer :: i, j
    real :: fTmp
    logical , save :: iFirst = .true.
    logical :: ifFound
    type(silam_material), pointer :: substTmp
    type(silam_species) :: speciesTmp
    
    cocktail_descr%defined = silja_false
    ifSpecies = .false.

!    ! Cocktails can be different: some prescribed names, or just
!    ! the name of the describing file.
!    !
!    if(iFirst) then
!      call read_standard_cocktail_descrs()
!      iFirst = .false.
!      if(error)return
!    endif

!    print *, 'Standard cocktails for checking: '
    ifFound = .false.
    do i=1,size(StandardCocktailPtr)
!      call msg(StandardCocktailPtr(i)%cocktail_name + ',' + chCocktailNm)
      if(fu_str_u_case(StandardCocktailPtr(i)%cocktail_name) == fu_str_u_case(chCocktailNm))then
        call copy_cocktail_descriptor(StandardCocktailPtr(i), cocktail_descr)
        ifFound = .true.
        exit
      endif
    end do
    !
    ! No standard cocktails found - read from chCocktail if it is a file name
    !
    if(.not. ifFound)then
      INQUIRE(file=chCocktailNm, exist=ifFound)
      if(.not. ifFound)then
        !
        ! Last chance: may be, it is species or "almost species"?
        ! Here we allow: (i) decodable species string, (ii) gas-phase chemical name
        ! Note that even if complete the species, the name of the cocktail must stay as given
        !
        call msg_warning('Strange string:' + chCocktailNm + &
                       & '- neither cocktail nor file name. Trying species...', &
                       & 'set_cocktail_description')
        speciesTmp = fu_species_from_short_string(chCocktailNm)  ! decodable species string?
        if(error)call unset_error('set_cocktail_description')
        !
        ! Now, what did we get?
        ! Attention: there can be substance in chCocktailNm, which erroneously
        ! recognised as no_mode species. This case must be understood and eliminated
        !
        if(defined(speciesTmp))then
          if(fu_mode(speciesTmp) == no_mode)then
            substTmp => fu_material(speciesTmp)
            if(fu_true(fu_if_gas(substTmp))) call set_species(speciesTmp, substTmp, in_gas_phase)
          endif
        else               ! if not, may be, gaseous substance?
          substTmp => fu_get_material_ptr(chCocktailNm)
          if(associated(substTmp))then
            if(fu_true(fu_if_gas(substTmp))) call set_species(speciesTmp, substTmp, in_gas_phase)
          else
            call set_error('Strange string:' + chCocktailNm + &
                       & '- not cocktail, not species, not file name', &
                       & 'set_cocktail_description')
            return
          endif   ! if substance name
        endif   ! if species can be deduced
        ifSpecies = .true.
        call set_descriptor_from_species(speciesTmp, cocktail_descr, chCocktailNm)
        return
      endif   ! if no file is found
      call read_cocktail_description(chCocktailNm, cocktail_descr)
    endif   ! if no standard cocktail found
    if(error)return

    !
    ! The cocktail fractions might be not scaled to 1.
    !
    call normalise_cocktail_description(cocktail_descr)
    
    cocktail_descr%defined = silja_true

  end subroutine set_cocktail_description


  !************************************************************************************
  
  subroutine set_descriptor_from_species(species, descr, force_name)
    !
    ! For the single species given, creates a cocktail with this species
    !
    type(silam_species), intent(in) :: species
    type(Tcocktail_descr), intent (out) :: descr
    character(len=*), intent(in) :: force_name
    
    ! Local variables
    integer :: iStat
    
    descr%defined = silja_false
    descr%nSpecies = 1
    descr%cocktail_name = trim(force_name)   !fu_str(species)
    descr%chFractionMassUnit = fu_basic_mass_unit(species%material)
    allocate(descr%mass_fractions(1), descr%species(1), descr%iEmisCocktSpeciesMapping(1), &
           & descr%f_descr_unit_2_species_unit(1), stat = iStat)
    if(iStat /= 0)then
      call set_error('Failed descriptor allocations for:' + fu_str(species), 'set_descriptor_from_species')
      return
    endif
    descr%mass_fractions(1) = 1.0
    descr%species(1) = species
    descr%iEmisCocktSpeciesMapping(1) = int_missing
    descr%f_descr_unit_2_species_unit(1) = 1.0
    descr%if_normalise = .true.
    if(error)return
    descr%defined = silja_true
    
  end subroutine set_descriptor_from_species

  
  !************************************************************************************

  subroutine read_cocktail_description(chCocktailFNm, descr)
    !
    ! Reads the cocktail description - either from the standard file
    ! or from the file pointed by the chCocktail itself.
    !
    implicit none

    ! Imported parameters
    character(len=*), intent(in) :: chCocktailFNm
    type(Tcocktail_descr), intent (out) :: descr

    ! Local variables
    integer :: status, fUnit
    logical :: eof
    character(len=clen) :: line

    descr%defined = silja_false

    ! First, let's try to find the standard file and see if the cocktail is there
    !
    fUnit = fu_next_free_unit()

    !
    ! Standard file is not found. Try to open the chCocktail like it is a file name
    !
    open(file=chCocktailFNm, unit=fUnit,action='read',status='old',iostat=status)
    if(status /= 0)then
      CALL set_error('Failed to find:' + chCocktailFNm, 'read_cocktail_description')
      RETURN
    endif
    !
    ! Find the first line of the description and read the first cocktail
    !
    eof = .false.
    do while(.not.eof)
      CALL next_line_from_input_file(fUnit, line, eof)
      if(error .or. eof)then
        call set_error('File:' + chCocktailFNm + ': has no cocktails', &
                     & 'read_cocktail_description')
        close(fUnit)
        return
      endif

      if(line =='COCKTAIL_DESCRIPTION_V3_2')then
        call read_one_description_v4_6(fUnit, descr) ! Namelist-based description
      else
        call msg_warning('Strange cocktail version:' + line + ', ignore', &
                       & 'read_cocktail_description')
        cycle
      endif

      close(fUnit)
      call normalise_cocktail_description(descr)
      descr%defined = silja_true
      return

    end do  ! while not eof
    !
    ! We must never come here
    !
    call set_error('File:' + chCocktailFNm + '- has no cocktails', &
                 & 'read_cocktail_description')
    close(fUnit)
    return

  end subroutine read_cocktail_description


  !*******************************************************************

  subroutine copy_descr_to_cocktail_link(descr1, descr2)
    !
    ! Copy (but not allocate) the mapping to emission species from
    ! descr1 to descr2.
    ! 
    implicit none
    type(Tcocktail_descr), intent(in) :: descr1
    type(Tcocktail_descr), intent(inout) :: descr2

    integer :: i

    if (.not. defined(descr1) .or. .not. defined(descr2)) then
      call set_error('Undefined descriptor', 'copy_descr_to_cocktail_link')
      return
    end if
    
    do i = 1, descr1%nspecies
      descr2%iEmisCocktSpeciesMapping(i) = descr1%iEmisCocktSpeciesMapping(i)
    end do
  end subroutine copy_descr_to_cocktail_link

  !*******************************************************************

  subroutine read_standard_cocktail_descrs(pItems, nItems)
    !
    ! Opens the standard file and reads all descriptions from there
    ! Fills in StandardCocktailPtr array.
    !
    implicit none

    ! Imported parameters
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: pItems
    integer, intent(in) :: nItems

    ! Local variables
    integer :: iCount, status, fUnit, iItem
    logical :: eof
    character(len=clen) :: line

    ! Stupidity check
    !
    if(nItems < 1)then
      call set_error('Strange number of standard cocktail files:' + fu_str(nItems), &
                   & 'read_standard_cocktail_descrs')
      return
    endif
    
    fUnit = fu_next_free_unit()
    if(error)return
    iCount = 0

    do iItem = 1, nItems
      open(file=fu_content(pItems(iItem)), unit=fUnit,action='read',status='old',iostat=status)
      IF (fu_fails(status == 0,'Failed to open the standard cocktails file:' + &
                            & fu_content(pItems(iItem)), 'read_standard_cocktail_descrs'))then
        nullify(StandardCocktailPtr)
        return
      endif
      !
      ! Standard file exists - count the number of cocktails
      !
      eof = .false.
      do while(.not.eof)
        CALL next_line_from_input_file(fUnit, line, eof)
        if(index(line,'COCKTAIL_DESCRIPTION_V3') == 1) iCount = iCount + 1
        if(error .or. eof)exit
      end do
      close(fUnit)
    end do  ! standard cocktail files

    !
    ! Having the full number of cocktails, allocate the structure and read the whole set
    !
    if(fu_fails(iCount >= 1,'No standard cocktails found','read_standard_cocktail_descrs')) then
      do iItem = 1, nItems
        call msg('Tried cocktail file:' + fu_content(pItems(iItem)))
      end do
      nullify(StandardCocktailPtr)
      return
    end if

    allocate(StandardCocktailPtr(iCount), stat = status)
    if(fu_fails(status == 0, 'Failed to allocate memory for cocktails','read_standard_cocktail_descrs'))return
    !
    ! Read the cocktails
    !
    iCount = 1
    do iItem = 1, nItems

      open(file=fu_content(pItems(iItem)), unit=fUnit,action='read',status='old',iostat=status)
      IF (fu_fails(status == 0,'Descriptions of standard cocktails cannot be red from:' + &
                            & fu_content(pItems(iItem)), 'read_standard_cocktail_descrs'))then
        nullify(StandardCocktailPtr)
        return
      endif
      !
      ! Go through the file and read all cocktails
      !
      eof=.false.
      do while(.not.eof)
        CALL next_line_from_input_file(fUnit, line, eof)
        if(error .or. eof)exit
        if(line == 'COCKTAIL_DESCRIPTION_V3_2')then
          call read_one_description_v4_6(fUnit, StandardCocktailPtr(iCount))
          if (error) return
        else
          call msg_warning('Strange cocktatil starting line:' + line + ', ignored', &
                         & 'read_standard_cocktail_descrs')
          cycle
        endif
        !
        ! Check uniqueness
        !
        do status = 1, iCount-1
          if(StandardCocktailPtr(status)%cocktail_name == StandardCocktailPtr(iCount)%cocktail_name)then
            call set_error('Duplicated cocktail name:' +  StandardCocktailPtr(status)%cocktail_name + &
                         & ', file=' + fu_content(pItems(iItem)), &
                         & 'read_standard_cocktail_descrs')
            return
          endif
        end do
        !
        ! When cocktail is read, its masses have to be normalised to fractions
        !
        call normalise_cocktail_description(StandardCocktailPtr(iCount))
        StandardCocktailPtr(iCount)%defined = silja_true
        iCount = iCount + 1
      end do  ! file
      close(fUnit)
    end do  ! standard cocktail files

  end subroutine read_standard_cocktail_descrs


  !***************************************************************

  subroutine read_one_description_v4_6(unit, descr)
    !
    ! Read one descriptor record from file unit. On succesful exit,
    ! the descriptor is completely allocated, but the arrays for unit
    ! conversion and emission cocktail mapping are not yet
    ! defined. The chemical mapper is created.
    !
    implicit none

    integer, intent(in) :: unit
    type(Tcocktail_descr), intent(out) :: descr

    
    type(Tsilam_namelist), pointer :: nlTmp
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: pItems ! output array
    type(silam_sp) :: spContent
    character(len=20) :: chTmp
    type(Taerosol) :: aerosol
    logical :: gas_phase
    integer :: nSubst, iSubst, iMode, gas_shift, iSpecies, nModes, nSpecies, iTMp, stat
    type(silam_species), dimension(max_species), save :: species_tmp
    type(Taerosol_mode) :: aero_mode
    real, dimension(:), pointer :: fractions, species_fractions
    real :: factor
    character(len=substNmLen) :: substNm
    character (len=*), parameter :: subname="read_one_description_v4_6"
   
    aerosol = aerosol_missing 
    fractions => fu_work_array()
    species_fractions => fu_work_array()

    nlTmp => fu_read_namelist(unit, .false., 'END_COCKTAIL_DESCRIPTION')
    spContent%sp => fu_work_string()

    descr%cocktail_name = fu_content(nlTmp, 'cocktail_name')

    call set_missing(aerosol, .true.)

    !
    ! Cocktail fractions may be mass fractions, number fractions,
    ! mole, volume, etc.
    !
    descr%chFractionMassUnit = fu_SI_unit(fu_content(nlTmp, 'mass_unit'))
    if (error .or. descr%chFractionMassUnit == '') then
      call set_error('No mass_unit for ' + descr%cocktail_name, subname)
      return
    end if
    !
    ! Set the aerosol from the namelist
    !
    nullify(pitems)
    call get_items(nlTmp, 'aerosol_mode', pItems, iTmp) ! returns iTmp= %n_modes
    if(iTmp > 0) then
      call set_aerosol(nlTmp, aerosol)
    endif
    if(error) return
    !
    ! Is the gas phase present ?
    !
    if (fu_str_u_case(fu_content(nlTmp,'gas_phase')) == 'YES') then
      gas_shift = 1
    else
      gas_shift = 0
    end if

    ! There used to be checking for radioactive species having a
    ! Junge-mode. So far there is no way to tell whether we have
    ! something radioactive, so we leave that for the radioactive
    ! transformation.
    
    !
    ! The component fraction lines. The first field is the name of the
    ! element, the second one and all others - its mass fractions in
    ! the mixture for each size mode and then for gaseous component,
    ! if any.
    !
    call get_items(nlTmp, 'component_fraction', pItems, nsubst)
    if(error .or. nSubst < 1)then
      call set_error('Strange component_fraction lines in the namelist:', subname)
      call report(nlTmp)
      return
    endif

    nSpecies = 0
    nModes = aerosol%n_modes
    do iSubst = 1, nSubst
      spContent%sp = fu_content(pItems(iSubst))
      read(unit=spContent%sp, iostat=stat, fmt=*) substNm, fractions(1:nModes+gas_shift)
      if (stat /= 0) then
        call set_error('Failed to parse fraction string:' + spContent%sp, subname)
        return
      end if
      
      do iMode = 1, nModes + gas_shift
        if (.not. fractions(iMode) > 0.0) cycle

        if (gas_shift > 0 .and. iMode == nModes + gas_shift) then
          aero_mode = in_gas_phase
        else
          aero_mode = aerosol%modes(iMode) 
        end if

        nSpecies = nSpecies + 1
        if (nSpecies > max_species) then
          call set_error('Too many species', subname)
          return
        end if
        call set_species(species_tmp(nSpecies), fu_get_material_ptr(substNm), aero_mode)
        do iSpecies = 1, nSpecies-1
          if (species_tmp(nSpecies) == species_tmp(iSpecies)) then
            call msg("")
            call msg("Namelist:")
            call report(nlTmp)
            call msg_warning("Duplicate species", subname)
            call report(species_tmp(iSpecies))
            call set_error("Duplicate species", subname)
            exit
          endif
        enddo
        if (error) return
        species_fractions(nSpecies) = fractions(iMode)
       end do
    end do

    descr%nSpecies = nSpecies

    allocate(descr%mass_fractions(nSpecies), &
           & descr%f_descr_unit_2_species_unit(nSpecies), &
           & descr%iEmisCocktSpeciesMapping(nSpecies),&
           & descr%species(nSpecies), stat=stat)
    if (stat /= 0) then
      call set_error('Allocate failed', subname)
      return
    end if

    do iSpecies = 1, nSpecies
      descr%species(iSpecies) = species_tmp(iSpecies)
      descr%mass_fractions(iSpecies) = species_fractions(iSpecies)
    end do
    !
    ! Having all species and the descriptor unit, can set up the conversion factors
    !
    do iSpecies = 1, nSpecies
      factor = fu_conversion_factor(fu_basic_mass_unit(descr%species(iSpecies)%material), &
                                  & descr%chFractionMassUnit, &
                                  & descr%species(iSpecies)%material)
      if (factor <= 0.0 .or. error) then
        call set_error('Invalid conversion factor:' + &
                     & fu_basic_mass_unit(descr%species(iSpecies)%material) + '->' + &
                     & descr%chFractionMassUnit + &
                     & ', for:' + fu_name(descr%species(ispecies)%material), subname)
        return
      end if
      descr%f_descr_unit_2_species_unit(iSpecies) = 1.0 / factor
    end do

    ! Is normalisation required?
    !
    ! The line if_normalize was quietly ignored
    if (fu_content(nlTmp,'if_normalize') /= '') then
      call msg("use 'if_normalise' instead of 'if_normalize'")
      call set_error("Too american spelling in the cocktails file!", subname)
      return
    endif

    if (fu_str_u_case(fu_content(nlTmp,'if_normalise')) == 'NO') then
      descr%if_normalise = .false.
    else
      descr%if_normalise = .true.
    end if

    descr%defined = silja_true
              
    call free_work_array(fractions)
    call free_work_array(species_fractions)
    call free_work_array(spContent%sp)
    call set_missing(aerosol, .false.) ! deallocate the possibly allocated memory

  end subroutine read_one_description_v4_6


  !***************************************************************

  subroutine copy_cocktail_descriptor(descr_in, descr_out)
    !
    ! Allocate descr_out and copy descr_in into it.
    !
    implicit none

    type(Tcocktail_descr), intent(in) :: descr_in
    type(Tcocktail_descr), intent(out) :: descr_out

    ! Local variables
    integer :: i, nModes

    ! JV 9/09: this subroutine no longer tries to deallocate descr_out
    ! - it is the caller's duty.
    
    if(.not.defined(descr_in))then
      call set_error('Cannot copy the undefined descriptor','copy_cocktail_descriptor')
      return
    endif

!    if(defined(descr_out))then
!      deallocate(descr_out%mass_fractions, descr_out%species, descr_out%iEmisCocktSpeciesMapping, &
!               & descr_out%f_descr_unit_2_species_unit, stat=i)
!      if(i /= 0)then
!        call set_error('Failed deallocation of the output defined descriptor', &
!                     & 'copy_cocktail_descriptor')
!        call unset_error('copy_cocktail_descriptor')
!      endif
!    endif

    descr_out%defined = silja_false

    allocate(descr_out%mass_fractions(descr_in%nSpecies), &
           & descr_out%species(descr_in%nSpecies), &
           & descr_out%iEmisCocktSpeciesMapping(descr_in%nSpecies), &
           & descr_out%f_descr_unit_2_species_unit(descr_in%nSpecies), stat=i)
    if (i /= 0) then
      call set_error('Failed to allocate memory','copy_cocktail_description')
      return
    endif

    descr_out%nSpecies = descr_in%nSpecies 
    descr_out%cocktail_name = descr_in%cocktail_name
    descr_out%chFractionMassUnit = descr_in%chFractionMassUnit

    do i = 1, descr_out%nSpecies
      descr_out%mass_fractions(i) = descr_in%mass_fractions(i)
      descr_out%iEmisCocktSpeciesMapping(i) = descr_in%iEmisCocktSpeciesMapping(i)
      descr_out%species(i) = descr_in%species(i)
      descr_out%f_descr_unit_2_species_unit(i) = descr_in%f_descr_unit_2_species_unit(i)
    end do

    descr_out%if_normalise = descr_in%if_normalise
    descr_out%defined = descr_in%defined

  end subroutine copy_cocktail_descriptor
  

  !***************************************************************

  subroutine normalise_cocktail_description(descr, fScale)
    !
    ! File is processed, but user is not forced to sum-up the mole fractions to 1
    ! So, let's do it for him. The scaling to be enforced: the sum of all mass
    ! fractions in all size modes is equal to 1.
    !
    implicit none

    ! Imported parameters 
    type(Tcocktail_descr), intent(inout) :: descr
    real, intent(out), optional :: fScale

    ! Local variables
    real :: vTmp
    integer :: i,j

    if (.not. descr%if_normalise) return

    vTmp =0.
    do i = 1, descr%nspecies
      if(descr%mass_fractions(i) < 0.)then
        call msg('Negative fraction for following species:')
        call report(descr%species(i))
        call set_error('Negative fraction','normalise_cocktail_description')
        return
      endif
      vTmp = vTmp + descr%mass_fractions(i)
    end do

    ! Scale the fractions, if needed
    !
    if(vTmp < 0. .or. ((vTmp*1.e10) .eps. 0.))then
      !
      ! If the total mass is unusable, we theoretically must return with error.
      ! However, one can try to zero all fractions and still continue. Risky, of course
      ! becaue theoretically sum of fractions must be 1 but we will see the outcome.
      !
      call msg_warning('Very small or non-positive fractions','normalise_cocktail_description')
      descr%mass_fractions(:) = 0.0

    elseif(.not.(vTmp .eps. 1.0))then
      
      descr%mass_fractions(:) = descr%mass_fractions(:) / vTmp

    endif   ! total mass <> 1 or 0

    if(present(fScale)) fScale = vTmp

  end subroutine normalise_cocktail_description


  !****************************************************************

  subroutine reset_used_cocktail_descr(descr)
    !
    ! Deallocates the cocktail descriptors, and set it missing. To set missing an
    ! undefined descriptor see below.
    !
    implicit none

    ! Imported paramneter
    type(Tcocktail_descr), intent(inout) :: descr

    ! Local variable
    integer :: i

    if (fu_fails(fu_true(descr%defined), 'Descriptor must be defined', 'reset_used_cocktail_descr')) return

    deallocate(descr%mass_fractions, descr%species, descr%iEmisCocktSpeciesMapping)
    call set_missing(descr)

  end subroutine reset_used_cocktail_descr
  

  !***************************************************************
  
  subroutine reset_new_cocktail_descr(descr)
    !
    ! Sets the given descriptor to a missing value. NO MEMORY DEALLOCATION.
    ! Do NOT call this function if the descriptor has allocated pointers.
    !
    implicit none

    ! Imported paramneter
    type(Tcocktail_descr), intent(inout) :: descr

    ! Local variable
    integer :: i

    descr%defined = silja_false
    descr%nSpecies = 0
    descr%cocktail_name = ''
    
    descr%chFractionMassUnit = ''
  
    nullify(descr%mass_fractions)
    nullify(descr%iEmisCocktSpeciesMapping)

  end subroutine reset_new_cocktail_descr


  !=======================================================================================
  function fu_emis_fractions_of_descr(descr)result(fractions)    ! encapsulation
    implicit none
    real, dimension(:), pointer:: fractions
    type(Tcocktail_descr), intent(in) :: descr
    fractions => descr%mass_fractions
  end function fu_emis_fractions_of_descr

  !=======================================================================================
  function fu_basic_mass_unit_descr(descr)result(chUnit)    ! encapsulation
    implicit none
    character(len=unitNmLen) :: chUnit
    type(Tcocktail_descr), intent(in) :: descr
    chUnit = descr%chFractionMassUnit
  end function fu_basic_mass_unit_descr

  !=======================================================================================
  function fu_name_of_descriptor(descr)result(chName)    ! encapsulation
    implicit none
    character(len=clen) :: chName
    type(Tcocktail_descr), intent(in) :: descr
    chName = descr%cocktail_name
  end function fu_name_of_descriptor

  !=======================================================================================
  function fu_species_of_descriptor(descr)result(pSpecies)    ! encapsulation
    implicit none
    type(silam_species), dimension(:), pointer :: pSpecies
    type(Tcocktail_descr), intent(in) :: descr
    pSpecies => descr%species
  end function fu_species_of_descriptor


  !************************************************************************************

  real function fu_conversion_factor_descr(descr, chUnitFrom, chUnitTo) result(factor)
    !
    ! Finds the factor converting the mass fractions or rates from one unit to another one.
    ! Time conversion is not a problem but mass units can require tricks: the
    ! descriptor fractions are given in the descriptor mass units.
    ! Procedure: split the descriptor rate to species rates in the original unitFrom, then
    ! convert rates species-wise to the unitTo, then sum-up.
    !
    implicit none

    ! Imported parameters
    type(Tcocktail_descr), intent(in) :: descr
    character(len=*), intent(in) :: chUnitFrom, chUnitTo

    ! Local variables
    character(len=unitNmLen) :: chMassUnitFrom, chMassUnitTo
    integer :: iSp
    real, dimension(:,:), pointer :: arTmp
    real :: factor_to_descr, factor_to_out

    ! Masses or rates?
    !
    if(index(chUnitFrom,'/') > 0)then
      if(index(chUnitTo,'/') == 0)THEN
        call set_error('Unit FROM is rate, unit TO is not:' + chUnitFrom + ',' + chUnitTo, &
                     & 'fu_conversion_factor_descr')
      endif
      !
      ! We are dealing with rates. Convert masses and denominators separately
      !
      chMassUnitFrom = chUnitFrom(1:index(chUnitFrom,'/')-1)
      chMassUnitTo = chUnitTo(1:index(chUnitTo,'/')-1)
    else
      !
      ! No denominators - masses only
      !
      if(index(chUnitTo,'/') > 0)THEN
        call set_error('Unit TO is rate, unit FROM is not:' + chUnitFrom + ',' + chUnitTo, &
                     & 'fu_conversion_factor_descr')
      endif
      chMassUnitFrom = chUnitFrom
      chMassUnitTo = chUnitTo
    endif

    !
    ! A hack to make the same source to emit cocktails with incompatible mass units
    if (chMassUnitFrom == "CocktailBasicUnit") then
      chMassUnitFrom = descr%chFractionMassUnit
    endif

    !
    ! Start from the mass conversion.  First of all, get the mass
    ! fractions of descriptor to the FROM unit.
    !
    arTmp => fu_work_array_2D()
    if(error)return

    if (descr%if_normalise) then
      ! 
      ! If the descriptor is normalised, we allow an arbitrary conversion chain.
      !
      do iSp = 1, descr%nSpecies
        arTmp(iSp,1) = descr%mass_fractions(iSp) * &
             & fu_conversion_factor(descr%chFractionMassUnit, chMassUnitFrom, &
                                  & descr%species(iSp)%material)
      end do
      !
      ! Fractions have to be re-normalised
      !
      arTmp(descr%nSpecies+1,1) = sum(arTmp(1:descr%nSpecies,1))
      arTmp(1:descr%nSpecies,2) = arTmp(1:descr%nSpecies,1) / arTmp(descr%nSpecies+1,1) ! fractions in FROM unit
      
      !
      ! The unit mass of the descriptor in the FROM unit is now the sum of the above fractions.
      ! Then this mass in the TO unit will be the sum of the converted fractions
      !
      do iSp = 1, descr%nSpecies
        arTmp(iSp,3) = arTmp(iSp,2) * fu_conversion_factor(chMassUnitFrom, chMassUnitTo, &
                                                         & descr%species(iSp)%material)
      end do
    
      factor = sum(arTmp(1:descr%nSpecies,3))  ! mass only so far

    else
      ! If the descriptor unit is not normalisable, we handle only the
      ! case where massUnitIn, descriptor mass unit and massUnitOut
      ! are basically the same. 
      !
      ! A material-dependent conversion from from the input unit to
      ! the descriptor one cannot be defined for non-normalised
      ! descriptors. From the descriptor unit to a different unitOut
      ! it could, but that would lead into trouble, because this is
      ! currently done elsewhere, and we would then need to distinguish
      ! between normalised and non-normalised fractions.
      ! 
      factor_to_descr = fu_conversion_factor(chMassUnitFrom, descr%chFractionMassUnit)
      if (error) then
        call set_error('Illegal unit conversion for non-normalised descriptor:' &
                     & // trim(chUnitFrom) // ' -> ' // trim(descr%chFractionMassUnit), &
                     & 'fu_conversion_factor_descr')
      end if
      factor_to_out = fu_conversion_factor(descr%chFractionMassUnit, chMassUnitTo)
      if (error) then
        call set_error('Illegal unit conversion for non-normalised descriptor:' &
                     & // trim(chUnitFrom) // ' -> ' // trim(descr%chFractionMassUnit), &
                     & 'fu_conversion_factor_descr')
      end if
      factor = factor_to_descr * factor_to_out
    end if
    
    !
    ! Now convert time
    !
    if(index(chUnitFrom,'/') > 0)then
      factor = factor / fu_conversion_factor(chUnitFrom(index(chUnitFrom,'/')+1:), &
                                           & chUnitTo(index(chUnitTo,'/')+1:))
    endif
    call free_work_array(arTmp)
  end function fu_conversion_factor_descr


  !************************************************************************************
  
  logical function fu_descriptors_equal(descr1, descr2) result(OK)
    !
    ! Overloads the comparison operator
    !
    implicit none

    ! Imported parameters
    type(Tcocktail_descr), intent(in) :: descr1, descr2

    ! Local variables
    integer :: i,j

    if(.not. descr1%defined == silja_true)then
     OK = (.not. (descr2%defined == silja_true))
     return
    endif

    OK = descr1%nSpecies == descr2%nSpecies &
         & .and. descr1%chFractionMassUnit == descr2%chFractionMassUnit

    if (.not. OK) return
      

    ! The species. It would seem that allowing different orders for
    ! the lists would be to unsafe already.

    do i = 1, descr1%nSpecies
      if (.not. (descr1%species(i) == descr2%species(i)) .or. &
        & .not. (descr1%mass_fractions(i) .eps. descr2%mass_fractions(i))) then
        OK = .false.
        return
      end if
    end do

  end function fu_descriptors_equal

  
  !****************************************************************
  
  integer function fu_descr_index_in_descr_by_name(chDescr, arDescrLst) result(iC)
    !
    ! Finds a descriptor given by name in the ilst of descriptors
    !
    implicit none
    ! Imported parameters
    type(Tcocktail_descr), dimension(:), intent(in) :: arDescrLst
    character(len=*), intent(in) :: chDescr
    
    do iC=1, size(arDescrlst)
      if(trim(fu_str_u_case(fu_name(arDescrlst(iC)))) == trim(fu_str_u_case(chDescr)))return
    end do
    iC = int_missing
  end function fu_descr_index_in_descr_by_name

  
  !***************************************************************

  function fu_cocktail_descriptior_ptr(chDescr) result(pDescr)
    !
    ! Returns the pointer to the cocktail descriptor, providing it is found in 
    ! StandardCocktailPtr.
    !
    implicit none

    ! Imported parameters
    character(len=*), intent(in) :: chDescr

    ! result
    type(Tcocktail_descr), pointer :: pDescr

    ! Local variables
    integer :: iDescr

    do iDescr = 1, size(StandardCocktailPtr)
!      print *, trim(StandardCocktailPtr(i)%cocktail_name)
      if(fu_str_u_case(StandardCocktailPtr(iDescr)%cocktail_name) == fu_str_u_case(chDescr))then
        pDescr => StandardCocktailPtr(iDescr)
        return
      endif
    end do
  end function fu_cocktail_descriptior_ptr

  
!  !*****************************************************************************
!
!  subroutine set_standard_cocktails_fname(chFNm)
!    implicit none
!    character(len=*), intent(in) :: chFNm
!    chCocktailDescrFile = trim(chFNm)
!  end subroutine set_standard_cocktails_fname


  !*******************************************************************
  !*******************************************************************
  !
  !  AEROSOL RULES
  !
  !*******************************************************************
  !*******************************************************************

  
  subroutine set_aerosol_rules_missing(rules)
    implicit none
    type(Taerosol_rules) :: rules
    
    rules%fDefaultRelHumidity = real_missing
    rules%defined = silja_false
  end subroutine set_aerosol_rules_missing

  !*******************************************************************

  function fu_set_aerosol_rules(nlSetup) result(rulesAerosol)
    !
    ! Creates the set of rules, in accordance with which the computations are to
    ! be made. Input - SILAM namelist
    !
    implicit none

    ! Return value
    !
    type(Taerosol_rules) :: rulesAerosol

    ! Imported parameter
    type(Tsilam_namelist), pointer :: nlSetup 

    if(.not.defined(nlSetup))then
      call set_error('Undefined namelist given','fu_set_aerosol_rules')
      rulesAerosol%defined = silja_false
      return
    endif
    
    !
    ! In theory, aerosol features depend on humidity. However, in simple cases we can ignore them
    ! just setting the default relative humidity value = 80%
    !
    if(fu_str_u_case(fu_content(nlSetup,'if_actual_humidity_for_particle_size')) == 'YES')then
      rulesAerosol%ifHumidityDependent = .true.
    elseif(fu_str_u_case(fu_content(nlSetup,'if_actual_humidity_for_particle_size')) == 'NO')then
      rulesAerosol%ifHumidityDependent = .false.
    else
      call report(nlSetup)
      call msg(fu_content(nlSetup,'if_actual_humidity_for_particle_size'))
      call set_error('if_actual_humidity_for_particle_size must be YES or NO','fu_set_aerosol_rules')
      return
    endif

    if(.not. rulesAerosol%ifHumidityDependent)then
      rulesAerosol%fDefaultRelHumidity = fu_content_real(nlSetup,'default_relative_humidity')
      if(error .or. rulesAerosol%fDefaultRelHumidity < 0 .or. rulesAerosol%fDefaultRelHumidity > 1.)then
        call msg('if_actual_humidity_for_particle_size = NO but strange default_relative_humidity')
        call report(nlSetup)
        call set_error('default_relative_humidity must be real [0..1]','fu_set_aerosol_rules')
        return
      endif
    endif

    rulesAerosol%ifNumberCnc = .false.  ! if needed, set to true later on

    !
    ! For Eulerian scheme we will need the minimum threshold of amount of species in the grid cell
    ! The cell will not be processed by ANY Eulerian-pool routine should it have smaller amount 
    ! than that threshold. A data-based setting is possible only when the total emitted amount
    ! is known, as well as a grid size. So far let's put it to something small enough.
    !
    rulesAerosol%defined = silja_true

  end function fu_set_aerosol_rules


  !********************************************************************************

  logical function fu_if_humidity_dependent(rules)
    !
    ! Returns the universal zero-mass threshold, below which the aerosol mass
    ! is considered to be practically zero
    !
    implicit none

    type(Taerosol_rules), intent(in) :: rules

    fu_if_humidity_dependent = rules%ifHumidityDependent

  end function fu_if_humidity_dependent


  !********************************************************************************

  real function fu_default_humidity(rules)
    !
    ! Returns the universal zero-mass threshold, below which the aerosol mass
    ! is considered to be practically zero
    !
    implicit none

    type(Taerosol_rules), intent(in) :: rules

    fu_default_humidity = rules%fDefaultRelHumidity

  end function fu_default_humidity


  !********************************************************************************

  function fu_lifetime_aerosol(rules)result(lifetime)
    !
    ! A typical life time due to degradation and deposition for aerosol
    ! So far, whatever the species are, we set the "typical" lifetime to be two days
    !
    implicit none

    type(silja_interval) :: lifetime
    type(Taerosol_rules), intent(in) :: rules

    lifetime = one_day * 2.0

  end function fu_lifetime_aerosol



  !**********************************************************************************
  !
  ! SILAM_AEROSOLS
  !
  !**********************************************************************************

!!$  subroutine init_aerosol_parameters(aerosol, nModes, shape_p, min_d, max_d, aver_d, &
!!$                                   & density, solubility, mass)
!!$    !
!!$    ! Initialises the aerosol parameters
!!$    !
!!$    implicit none
!!$    
!!$    ! Imported variables
!!$    type(Taerosol), intent(inout) :: aerosol
!!$    integer, intent(in) :: nModes
!!$    real, intent(in) :: shape_p
!!$    real, dimension(:), intent(in) :: min_d, max_d, aver_d, mass, density
!!$    integer, dimension(:), intent(in) :: solubility
!!$    ! local variables
!!$    integer :: i
!!$
!!$    aerosol%n_modes = nModes
!!$
!!$    if(size(min_d) < nModes .or. size(max_d) < nModes .or. size(aver_d) < nModes)then
!!$      call set_error('Inconsistent sizes of arrays','init_aerosol_parameters')
!!$      return
!!$    endif
!!$
!!$    if(nModes > 0)then
!!$      allocate(aerosol%modes(nModes), stat=i) 
!!$      if(i /= 0)then
!!$        call set_error('Failed to allocate memory','init_aerosol_parameters')
!!$        return
!!$      endif
!!$      do i = 1, nModes
!!$        
!!$        aerosol%modes(i)%mode_type = shape_p
!!$        aerosol%modes(i)%minD = min_d(i)
!!$        aerosol%modes(i)%maxD = max_d(i)
!!$        aerosol%modes(i)%meanD = aver_d(i)
!!$        aerosol%modes(i)%density = density(i)
!!$        aerosol%modes(i)%solubility = solubility(i)
!!$      end do
!!$    else
!!$      nullify(aerosol%modes)
!!$    endif
!!$
!!$  end subroutine init_aerosol_parameters


  !*****************************************************************

  subroutine set_aerosol_missing(aerosol, ifNew)
    !
    ! Substitutes the aerosol_missing parameter. If it was somehow defined
    ! the allocatable parts are carefully destroyed.
    ! To avoid problems - use iStat to block the runtime error generation. 
    !
    implicit none

    type(Taerosol), intent(inout) :: aerosol
    logical, intent(in) :: ifNew  ! if this is a newly allocated aerosol alias

    ! Local variables
    integer :: iStat

!!$    if(.not. ifNew)then
!!$      if(aerosol%defined == silja_true)then
!!$        if(aerosol%n_modes > 0 .and. aerosol%n_modes < 100)then  ! Something reasonable
!!$          if(associated(aerosol%modes))deallocate(aerosol%modes, stat=iStat)
!!$          if(iStat /= 0)then
!!$            call msg('Memory deallocation failed, iStat=',iStat)
!!$            call msg('Number of aerosol modes=',aerosol%n_modes)
!!$            call msg_warning('Memory deallocation failed','set_aerosol_missing')
!!$            call msg_warning('set_aerosol_missing')
!!$          endif
!!$        endif
!!$      endif
!!$    endif
    aerosol%n_modes = 0
    nullify(aerosol%modes)
    aerosol%defined = silja_false

  end subroutine set_aerosol_missing


  !******************************************************************

  subroutine set_aerosol_from_namelist(nlInput, aerosol)
    !
    ! Gets the aerosol size spectrum from namelist. It contains 
    !    integer :: n_modes, shape_p
    !    real :: density
    !    REAL, DIMENSION(:), POINTER :: min_d, max_d, aver_d
    !
    ! Density and shape_p are optional, while modes are described
    ! with min, max and average diameter and mass fraction
    ! In the description, each aerosol mode
    ! occupies one string where the mean diameter, a distribution function and
    ! fraction of the cocktail mass placed in this mode are defined. 
    ! Density has to be defined separately and evidently the same
    ! for all modes - same as chemical composition
    !
    ! This function creates a new aerosol alias and returns pointer on it.
    !
    implicit none

    ! Imported variables
    type(Tsilam_namelist), intent(in) :: nlInput
    type(Taerosol), intent(out) :: aerosol

    ! Local variables
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: pItems ! output array
    type(silam_sp) :: spContent
    character(len=20) :: chTmp
    integer :: status, iTmp, iMode, iShape

    aerosol%defined = silja_false
    !
    ! Get the number of modes in aerosol and set it
    !
    nullify(pItems)
    call get_items(nlInput, 'aerosol_mode', pItems, iTmp) ! iTmp= %n_modes
    if(iTmp < 1)then
      call set_missing(aerosol, .true.)
      return
    endif
    !
    ! Allocate necessary memory for aerosol description
    !
    !if(aerosol%n_modes > 0 .and. aerosol%n_modes < 100) deallocate(aerosol%modes, stat=status)
    allocate(aerosol%modes(iTmp), stat=status)
    if(status /= 0)then
      call msg('Requested number of aerosol modes:', iTmp)
      call set_error('Failed to allocate aerosol modes','set_aerosol_from_namelist')
      return
    endif
    aerosol%n_modes = iTmp
    do iMode = 1, aerosol%n_modes
      aerosol%modes(iMode)%defined = silja_false
    end do

    spContent%sp => fu_work_string()
    if(error)return

      ! Shape of the mass distribution inside the nodes
      !
      select case(fu_str_u_case(fu_content(nlInput,'mode_distribution_type')))
        case('FIXED_DIAMETER')
          iShape = fixed_diameter_flag
        case('LOGNORMAL')
          iShape = lognormal_flag
        !case('GAMMA_FUNCTION')
        !  iShape = gamma_function_flag
        case default
          call set_error('mode_distribution_type is only fixed_diameter and lognormal so far', &
                       & 'set_aerosol_from_namelist')
          return
      end select

      !
      ! Aerosol basic parameters
      !
      ! Scan modes and read the diameters: min, max, mean, and unit
      !
      do iTmp = 1, aerosol%n_modes
        spContent%sp = fu_content(pItems(iTmp))
        read(unit=spContent%sp, iostat=status,fmt=*) iMode
        if(iMode <=0 .or. iMode > aerosol%n_modes)then
          call msg('Strange mode number in aerosol namelist:',iMode)
          call report(nlInput)
          call set_error('Strange mode number in aerosol namelist','set_aerosol_from_namelist')
          return
        else
          read(unit=spContent%sp, iostat=status,fmt=*) iMode, &
                                                     & aerosol%modes(iMode)%fp1, &
                                                     & aerosol%modes(iMode)%fp2, &
                                                     & aerosol%modes(iMode)%mass_mean_d, &
                                                     & chTmp
        endif
        if(status /= 0)then
          call set_error(fu_connect_strings('Cannot read fraction string:',spContent%sp), &
                       & 'set_aerosol_from_namelist')
          return
        endif

        aerosol%modes(iMode)%fp1 = aerosol%modes(iMode)%fp1 * fu_conversion_factor(chTmp,'m')
        aerosol%modes(iMode)%fp2 = aerosol%modes(iMode)%fp2 * fu_conversion_factor(chTmp,'m')
        aerosol%modes(iMode)%mass_mean_d = aerosol%modes(iMode)%mass_mean_d &
                                         & * fu_conversion_factor(chTmp,'m')
        aerosol%modes(iMode)%nominal_d = aerosol%modes(iMode)%mass_mean_d ! Set it now and forever
        aerosol%modes(iMode)%distr_type = iShape
        if(error)return
        aerosol%modes(iMode)%defined = silja_true

      end do  ! over the namelist mode items

      !
      ! Now set the mode names, which are in the separate set of lines
      !
      call get_items(nlInput, 'aerosol_mode_name', pItems, status)
      if(status == 0)then
        do iMode = 1, aerosol%n_modes
          aerosol%modes(iMode)%name = ''
        end do
      else
        do iTmp = 1, status
          spContent%sp = fu_content(pItems(iTmp))
          read(unit=spContent%sp, iostat=status,fmt=*) iMode
          if(iMode <=0 .or. iMode > aerosol%n_modes)then
            call msg('Strange mode number in aerosol namelist:',iMode)
            call report(nlInput)
            call set_error('Strange mode number in aerosol namelist','set_aerosol_from_namelist')
            return
          endif
          read(unit=spContent%sp, iostat=status,fmt=*) iMode, aerosol%modes(iMode)%name 
        end do
      endif

      !
      ! Now set the mode solubility, which are in the separate set of lines
      !
      call get_items(nlInput, 'aerosol_mode_relative_solubility', pItems, status)
      if(status == 0)then
        do iMode = 1, aerosol%n_modes
          aerosol%modes(iMode)%solubility = int_missing
        end do
      else
        do iTmp = 1, status
          spContent%sp = fu_content(pItems(iTmp))
          read(unit=spContent%sp, iostat=status,fmt=*) iMode
          if(iMode <=0 .or. iMode > aerosol%n_modes)then
            call msg('Strange mode number in aerosol namelist:',iMode)
            call report(nlInput)
            call set_error('Strange mode number in aerosol namelist','set_aerosol_from_namelist')
            return
          endif
          read(unit=spContent%sp, iostat=status,fmt=*) iMode, aerosol%modes(iMode)%solubility
        end do
      endif

      call free_work_array(spContent%sp)

      aerosol%defined = silja_true

  end subroutine set_aerosol_from_namelist


  !***************************************************************************

  subroutine set_aerosol_mode_from_params(mode, chNm, fp1, fp2, &
                                        & mass_mean_d, density, distr_type, solubility)
    !
    ! Just sets the aerosol mode with no checking
    !
    implicit none
    type(Taerosol_mode), intent(out) :: mode
    character(len=*), intent(in) :: chNm
    integer, intent(in) :: distr_type, solubility
    real, intent(in) :: fp1, fp2, mass_mean_d, density

    mode%name = chNm
    mode%fp1 = fp1
    mode%fp2 = fp2
    mode%mass_mean_d = mass_mean_d
    mode%nominal_d = mode%mass_mean_d ! Set it now and forever
    mode%distr_type = distr_type
    mode%solubility = solubility
    mode%defined = silja_true

  end subroutine set_aerosol_mode_from_params


  !***************************************************************************

  subroutine set_aerosol_mode_from_d(mode, d)
    !
    ! Just sets the simple aerosol mode with no checking
    !
    implicit none
    type(Taerosol_mode), intent(out) :: mode
    real, intent(in) :: d

    mode%name = ''
    mode%fp1 = d
    mode%fp2 = d
    mode%mass_mean_d = d
    mode%nominal_d = d
    mode%distr_type = fixed_diameter_flag
    mode%solubility = int_missing
    mode%defined = silja_true

  end subroutine set_aerosol_mode_from_d



  !***************************************************************************

  function fu_mode_of_aerosol(aerosol, iMode)result(aerMode)
    !
    ! Encapsumation
    !
    implicit none
    ! Imported parameter
    type(Taerosol), intent(in) :: aerosol
    integer, intent(in) :: iMode

    ! Return value
    type(Taerosol_mode) :: aerMode

    if(iMode < 1 .or. imode > aerosol%n_modes)then
      call msg('iMode is outside the allowed range (1:nModes):',iMode, aerosol%n_modes)
      call set_error('Strange iMode','fu_mode_of_aerosol')
    else
      aerMode = aerosol%modes(iMode)
    endif

  end function fu_mode_of_aerosol


  !***************************************************************************

  subroutine copy_species(species_out, species_in)
    !
    ! Copies one species into another one
    !
    implicit none

    ! Imported parameters
    type(silam_species), intent(in) :: species_in
    type(silam_species), intent(out) :: species_out

    species_out%mode = species_in%mode 
    species_out%wavelength = species_in%wavelength
    species_out%material => species_in%material
    species_out%defined = species_in%defined

  end subroutine copy_species

  !***************************************************************************

  logical function fu_compare_modes_bigger(mode1, mode2)
    !
    ! Compares the two modes for one being bigger than the other
    ! Strict comparison!
    !
    implicit none

    type(Taerosol_mode), intent(in) :: mode1, mode2

    fu_compare_modes_bigger = (fu_compare_modes_mean_diam(mode1, mode2) == 1)

  end function fu_compare_modes_bigger


  !*************************************************************************

  function fu_compare_modes_mean_diam(mode1, mode2) result(cmp)
    !
    ! Compares two modes. Rules: gas pahse is always smaller than aerosol phase,
    ! smaller mean diameters mean smaller modes, same-diameters less soluble mode is smaller.
    !
    implicit none

    ! Imported parameters
    type(Taerosol_mode), intent(in) :: mode1, mode2

    ! Return value
    integer :: cmp ! -1, 0, 1


    cmp = 0.

    ! We may later think what to do for modes with different
    ! type. Right now, I'll only allow fixed diameter and gas phase
    ! (which is considered smaller than any aerosol).
    if (mode1%distr_type < mode2%distr_type) then
      cmp = -1
      return
    else if (mode1%distr_type > mode2%distr_type) then
      cmp = 1
      return
    end if

    ! From now on, we hase same type of modes to compare.
    !
    if (mode1%distr_type == gas_phase_flag) then  ! two gas phases
      cmp = 0
      return
    end if

!    if (mode1%mass_mean_d .eps. mode2%mass_mean_d) then
!      cmp = 0
!      return
!    end if

    if (mode1%mass_mean_d < mode2%mass_mean_d) then
      cmp = -1 
      return
    else if (mode1%mass_mean_d > mode2%mass_mean_d) then
      cmp = 1
      return
    end if

        
    if(mode1%solubility < mode2%solubility)then  ! diameters are same. Solubility then plays
      cmp = -1
      return
    elseif(mode1%solubility > mode2%solubility)then
      cmp = 1
      return
    end if
          
    ! Same mean, same type, same solubility...
    ! The last option is mode name. Empty name is always smaller, 
    ! otherwise alphabetical comparison just goes.
    !
!    return
!    !
!    ! Below comparison of the mode names has been disabled since
!    ! otherwise recognision of the modes from strings becomes impossible.
!    !
    if(trim(mode1%name) == trim(mode2%name))return
    
    if(len_trim(mode1%name) == 0)then   ! first zero-length
      if(len_trim(mode2%name) == 0)then   ! second zero-length
        return             ! same modes
      else
        cmp = -1           ! first smaller (zero-len)
        return
      endif
    elseif(len_trim(mode2%name) == 0)then   ! second smaller (zero-len)
      cmp = 1
      return
    else
      if(mode1%name < mode2%name)then  ! both names make sense
        cmp = -1
        return
      else
        cmp = 1
        return
      end if
    endif
    
  end function fu_compare_modes_mean_diam


  !*************************************************************************

  function fu_compare_modes_general(mode1, mode2) result(cmp)
    !
    ! Compares two modes. Rules: gas phase is always smaller than aerosol phase,
    ! smaller parameter fp1 means smaller mode, then fp2 and ip, then less soluble mode is smaller.
    !
    implicit none

    ! Imported parameters
    type(Taerosol_mode), intent(in) :: mode1, mode2

    ! Return value
    integer :: cmp ! -1, 0, 1


    cmp = 0.

    !
    ! Two undefined modes are identical, otherwise the undefined mode is smaller
    !
    if(.not. (mode1%defined == silja_true))then
      if(mode2%defined == silja_true) cmp = -1    ! undefined mode 1 and defined mode 2
      return
    endif

    ! We may later think what to do for modes with different
    ! type. Right now, I'll only allow fixed diameter and gas phase
    ! (which is considered smaller than any aerosol).
    if (mode1%distr_type < mode2%distr_type) then
      cmp = -1
      return
    else if (mode1%distr_type > mode2%distr_type) then
      cmp = 1
      return
    end if

    ! From now on, we have same type of modes to compare.
    !
    if (mode1%distr_type == no_mode_flag) then  ! two no_modes
      cmp = 0
      return
    end if
    if (mode1%distr_type == gas_phase_flag) then  ! two gas phases
      cmp = 0
      return
    end if
    !
    ! fp1 means different things for different distribution types but can be the diameter
    ! then it has to be multiplied with 1e6 at least to allow .eps. comparison
    !
    if(.not. (mode1%fp1*1e6 .eps. mode2%fp1*1e6))then
      if (mode1%fp1 < mode2%fp1) then
        cmp = -1 
        return
      else
        cmp = 1
        return
      end if
    endif
    !
    ! fp1 means different things for different distribution types but can be the diameter
    ! then it has to be multiplied with 1e6 at least to allow .eps. comparison
    !
    if(.not. (mode1%fp2*1e6 .eps. mode2%fp2*1e6))then
      if (mode1%fp2 < mode2%fp2) then
        cmp = -1 
        return
      else
        cmp = 1
        return
      end if
    endif
    !
    ! fp* are the same. The last-but-one option is solubility
    !
    if(mode1%solubility < mode2%solubility)then  ! diameters are same. Solubility then plays
      cmp = -1
      return
    elseif(mode1%solubility > mode2%solubility)then
      cmp = 1
      return
    end if
    !
    ! All numbers are the same. The last option is mode name. Empty name is always smaller, 
    ! otherwise alphabetical comparison just goes.
    !
    if(trim(mode1%name) == trim(mode2%name))return

    if(len_trim(mode1%name) == 0)then   ! first zero-length
      if(len_trim(mode2%name) == 0)then   ! second zero-length
        return             ! same modes
      else
        cmp = -1           ! first smaller (zero-len)
        return
      endif
    elseif(len_trim(mode2%name) == 0)then   ! second smaller (zero-len)
      cmp = 1
      return
    else
      if(mode1%name < mode2%name)then  ! both names make sense
        cmp = -1
        return
      else
        cmp = 1
        return
      end if
    endif

  end function fu_compare_modes_general


  !*****************************************************************

  subroutine set_aerosol_from_copy(aerosol_out, aerosol_in)
    !
    ! Since aerosol involves memory allocation, one has to ensure proper copying
    !
    implicit none

    ! Imported parameters
    type(Taerosol), intent(in) :: aerosol_in
    type(Taerosol), intent(inout) :: aerosol_out

    ! Local variables
    integer :: status
    logical :: ifAllocate

    if(aerosol_in%n_modes <= 0 .or. .not.defined(aerosol_in)) then
      call set_missing(aerosol_out, .true.)
      return
    endif
    !
    ! UNFORTUNATELY, with loosy memory functions of FORTRAN, this checking is 
    ! dangerous. Maximum what is possible is to deallocate memory but even this is 
    ! not always safe. What is below is a compromise
    !
    ! May be, aerosol_out is already something reasonable, then
    ! let's handle it properly: we might not need to re-allocate memory
    ! if the size of the _out arrays is correct
    !

!    call msg('aer_1')
    if(aerosol_out%defined == silja_true)then
!      call msg('aer_2')
      if(aerosol_out%n_modes > 0 .and. aerosol_out%n_modes < 100) then
!        call msg('aer_3')
        if(associated(aerosol_out%modes))then
          !
          ! aerosol_out is not empty. May be, it is the same as the aerosol_in?
          !
!          call msg('aer_4')
          if(aerosol_out%n_modes /= aerosol_in%n_modes)then
            deallocate(aerosol_out%modes, stat=status)
            ifAllocate = .true.
!            call msg('aer_5, status:', status)
          else
            ifAllocate = .false.
          endif
        else
          ifAllocate = .true.
        endif
      else
        ifAllocate = .true.
      endif
    else
      ifAllocate = .true.
    endif

!    call msg('aer_6, allocate nModes:',aerosol_in%n_modes)
    if(ifAllocate)then
      allocate(aerosol_out%modes(aerosol_in%n_modes), stat=status)
      if(status /= 0)then
        call msg('Failed memory allocation: status=',status)
        call set_error('Failed to allocate memory','set_aerosol_from_copy')
        return
      endif
    endif
!    call msg('aer_7')
    !
    ! Now reset the _out variables. It is faster than comparing and, if needed, resetting
    !
    aerosol_out%n_modes = aerosol_in%n_modes
!    call msg('aer_8')
    aerosol_out%modes(1:aerosol_out%n_modes) = aerosol_in%modes(:)

    aerosol_out%defined = silja_true
!    call msg('aer_9')

  end subroutine set_aerosol_from_copy


  !*****************************************************************

  integer function fu_n_modes_of_spectrum(aerosol)
    implicit none
    type(Taerosol), intent(in) :: aerosol
    fu_n_modes_of_spectrum = aerosol%n_modes
  end function fu_n_modes_of_spectrum


  !*****************************************************************

  logical function fu_aerosol_defined(aerosol)
    implicit none
    type(Taerosol), intent(in) :: aerosol
!    fu_aerosol_defined = aerosol%n_modes > 0 .and. aerosol%n_modes < 100
    fu_aerosol_defined = aerosol%defined == silja_true

  end function fu_aerosol_defined


  !******************************************************************

  subroutine report_aerosol(aerosol)
    implicit none
    type(Taerosol), intent(in) :: aerosol
    integer :: iTmp

    if(aerosol%defined == silja_true)then
      call msg('Number of aerosol modes:',aerosol%n_modes)
      do iTmp=1,aerosol%n_modes
        call msg('Mean mode diameter and mode number:', iTmp, aerosol%modes(iTmp)%mass_mean_d)
      end do
    else
      call msg('Undefined aerosol')
    endif

  end subroutine report_aerosol


  !******************************************************************

  subroutine report_aerosol_mode(aerosol_mode)
    implicit none
    type(Taerosol_mode), intent(in) :: aerosol_mode

    call msg('Mode name, distribution type, solubility:' + aerosol_mode%name + ',' + &
           & fu_name(aerosol_mode%distr_type), aerosol_mode%solubility)
    if(aerosol_mode%distr_type /= gas_phase_flag)then
      call msg('fp1, fp2, nom.diam.:', (/aerosol_mode%fp1, aerosol_mode%fp2, aerosol_mode%nominal_d/))
      if(aerosol_mode%ifMeanDiamDynamic)then
        call msg('Mean diameter is dynamic, value now:',aerosol_mode%mass_mean_d)
      else
        call msg('Mean diameter is static:',aerosol_mode%mass_mean_d)
      endif
      if(aerosol_mode%distr_type /= gamma_function_flag) &
                                  & call msg('Gamma function parameter:',aerosol_mode%ip1)
    endif

  end subroutine report_aerosol_mode


  !******************************************************************

  subroutine report_as_namelist_aer_mode(aer_mode, iUnit)
    !
    ! Writes down the aerosol mode as a namelist, so that later it can be read by fu_set_mode_from_namelist
    !
    implicit none
    
    ! Imported parameters
    type(Taerosol_mode), intent(in) :: aer_mode
    integer, intent(in) :: iUnit

    if(.not. (aer_mode%defined == silja_true))then
      call msg_warning('Undefined mode given','report_as_namelist_aer_mode')
      write(iUnit,fmt='(A)')'mode_distribution_type = UNDEFINED'
      return
    endif

    write(iUnit,fmt='(2A)')'  mode_name = ', aer_mode%name
    if(aer_mode%solubility /= int_missing) write(iUnit,fmt='(A,I4)')'  mode_solubility = ', &
                                                                  & aer_mode%solubility
!    if(.not.(aer_mode%density .eps. real_missing)) then 
!      write(iUnit,fmt='(A,F15.7,A)')'  mode_density = ', aer_mode%density,' kg/m3'
!    end if
    select case(aer_mode%distr_type)
      case(gas_phase_flag)
        write(iUnit,fmt='(A)')'  mode_distribution_type = GAS_PHASE'
        
      case(no_mode_flag)
        write(iUnit,fmt='(A)')'  mode_distribution_type = NO_MODE'
        
      case(fixed_diameter_flag)
        write(iUnit,fmt='(A)')'  mode_distribution_type = FIXED_DIAMETER'
        write(iUnit,fmt='(A,F15.7,A)')'  mode_nominal_diameter =', aer_mode%nominal_d*1e6,' mkm'
        write(iUnit,fmt='(A,F15.7,A)')'  fix_diam_mode_min_diameter = ',aer_mode%fp1*1.e6,' mkm'
        write(iUnit,fmt='(A,F15.7,A)')'  fix_diam_mode_max_diameter = ',aer_mode%fp2*1.e6,' mkm'
        if(.not.(aer_mode%mass_mean_d .eps. real_missing)) then
          write(iUnit,fmt='(A,F15.7,A)')'  fix_diam_mode_mean_diameter = ', &
               & aer_mode%mass_mean_d*1.e6,' mkm'
        end if
      case(moving_diameter_flag)
        write(iUnit,fmt='(A)')'  mode_distribution_type = MOVING_DIAMETER'
        write(iUnit,fmt='(A,F15.7,A)')'  mode_nominal_diameter =', aer_mode%nominal_d*1e6,' mkm'
        write(iUnit,fmt='(A,F15.7,A)')'  moving_diam_mode_min_diameter = ',aer_mode%fp1*1.e6,' mkm'
        write(iUnit,fmt='(A,F15.7,A)')'  moving_diam_mode_max_diameter = ',aer_mode%fp2*1.e6,' mkm'
        if(.not.(aer_mode%mass_mean_d .eps. real_missing)) then 
          write(iUnit,fmt='(A,F15.7,A)')'  moving_diam_mode_mean_diameter = ', &
               & aer_mode%mass_mean_d*1.e6,' mkm'
        end if
      case(gamma_function_flag)
        write(iUnit,fmt='(A)')'  mode_distribution_type = GAMMA_FUNCTION'
        write(iUnit,fmt='(A,F15.7,A)')'  mode_nominal_diameter =', aer_mode%nominal_d*1e6,' mkm'
        write(iUnit,fmt='(A,I15)')'  gamma_mode_gamma_number = ',aer_mode%ip1

      case(lognormal_flag)
        write(iUnit,fmt='(A)')'  mode_distribution_type = LOGNORMAL'
        write(iUnit,fmt='(A,F15.7,A)')'  mode_nominal_diameter =', aer_mode%nominal_d*1e6,' mkm'
        write(iUnit,fmt='(A,F15.7,A)') &
             & '  lognormal_mode_mean_diameter = ',aer_mode%fp1*1.e6,' mkm'
        write(iUnit,fmt='(A,F15.7,A)') &
             & '  lognormal_mode_sigma = ',aer_mode%fp2*1.e6,' mkm'

      case(sea_salt_mode_flag)
        write(iUnit,fmt='(A)')'  mode_distribution_type = SEA_SALT'
        write(iUnit,fmt='(A,F15.7,A)')'  mode_nominal_diameter =', aer_mode%nominal_d*1e6,' mkm'
        write(iUnit,fmt='(A,F15.7,A)')' &
             & sea_salt_mode_min_diameter = ',aer_mode%fp1*1.e6,' mkm'
        write(iUnit,fmt='(A,F15.7,A)') &
             & '  sea_salt_mode_max_diameter = ',aer_mode%fp2*1.e6,' mkm'

      case default
        call msg('Unknown mode type:',aer_mode%distr_type)
        call set_error('Unknown mode type','report_as_namelist_aer_mode')
        return
    end select

  end subroutine report_as_namelist_aer_mode


  !********************************************************************

  real function fu_mean_d_aerosol(aerosol, iMode)
    !
    ! Just encapsulation fo the mean diameter
    !
    implicit none
    type(Taerosol), intent(in) :: aerosol
    integer, intent(in) :: iMode

    if(aerosol%n_modes < 1 .or. aerosol%n_modes > 100 .or. iMode > aerosol%n_modes)then
      call set_error('Undefined aerosol or iMode>n','fu_mean_d_aerosol')
      fu_mean_d_aerosol = real_missing
    else
      fu_mean_d_aerosol = aerosol%modes(iMode)%mass_mean_d
    endif


  end function fu_mean_d_aerosol


  !********************************************************************


  real function fu_min_d_aerosol(aerosol, iMode)
    !
    ! Just encapsulation fo the mean diameter
    !
    implicit none
    type(Taerosol), intent(in) :: aerosol
    integer, intent(in) :: iMode

    if(aerosol%n_modes < 1 .or. aerosol%n_modes > 100 .or. iMode > aerosol%n_modes)then
      call set_error('Undefined aerosol','fu_min_d_aerosol')
      return
    endif

    fu_min_d_aerosol = fu_min_d(aerosol%modes(imode))

  end function fu_min_d_aerosol


  !********************************************************************

  real function fu_max_d_aerosol(aerosol, iMode)
    !
    ! Just encapsulation fo the mean diameter
    !
    implicit none
    type(Taerosol), intent(in) :: aerosol
    integer, intent(in) :: iMode

    if(aerosol%n_modes < 1 .or. aerosol%n_modes > 100 .or. iMode > aerosol%n_modes)then
      call set_error('Undefined aerosol','fu_max_d_aerosol')
      return
    endif

    fu_max_d_aerosol = fu_max_d(aerosol%modes(imode))

  end function fu_max_d_aerosol


  !************************************************************************

  integer function fu_aerosol_mode_index(aerosol, fDiam)
    !
    ! Finds out the mode index into which the given diameter falls
    !
    implicit none

    ! Imported parameters 
    type(Taerosol), intent(in) :: aerosol
    real, intent(in) :: fDiam

    ! Local variables
    integer :: iMode

    if(aerosol%n_modes < 1 .or. aerosol%n_modes > 100)then
      call set_error('Undefined aerosol or wrong mode number','fu_mean_d_aerosol')
      return
    endif
    !
    ! If a given diam is inside one of the modes - get it and return
    !
    do iMode = 1, aerosol%n_modes
      if(fu_min_d(aerosol%modes(iMode)) <= fDiam &
       & .and. fDiam < fu_max_d(aerosol%modes(iMode)))then
        fu_aerosol_mode_index = iMode
        return
      endif
    end do
    if(fDiam == fu_max_d(aerosol%modes(aerosol%n_modes)))then ! at the upper edge??
      fu_aerosol_mode_index = aerosol%n_modes
      return
    endif
    !
    ! Given diameter is outside the size spectrum of aerosol
    !
    call msg('Given aerosol size is outside the range', fDiam)
    call report(aerosol)
    call set_error('Given aerosol size is outside the range','fu_aerosol_mode_index')
    fu_aerosol_mode_index = int_missing

  end function fu_aerosol_mode_index


  !******************************************************************

  logical function fu_compare_aerosols_eq (aerosol1, aerosol2) result(OK)
    !
    ! Comparison operator encapsulation
    !
    implicit none
    type(Taerosol), intent(in) :: aerosol1, aerosol2
    integer :: i

     OK = aerosol1%n_modes == aerosol2%n_modes
     if(.not.OK) return
     !
     ! If no aerosol modes - skip the rest
     !
     if(aerosol1%n_modes > 0)then
       do i=1,aerosol1%n_modes
         OK = (aerosol1%modes(i) == aerosol2%modes(i))
         if(.not.OK) return
       end do
     endif

  end function fu_compare_aerosols_eq


  !**********************************************************************************
  !
  ! Unit test
  !
  !**********************************************************************************

  subroutine chemical_setup_tests(nlSetup)
    implicit none
    type(Tsilam_namelist), pointer :: nlSetup
    type(silam_species) :: species1, species2, species3
    type(silam_species), dimension(500), target :: species_list
    type(Taerosol_mode) :: small, big, strange

    integer :: nspecies = 0, i, j, k, bignum = 100
    type(chemical_mapper) :: mapper
    logical :: ifSpecies
    type(Tcocktail_descr) :: descr1

    type(silam_species), dimension(:), pointer :: species_list_ptr

    species_list_ptr => species_list

    call init_chemical_materials(fu_expand_environment(fu_content(nlSetup,'chemical_database_fnm')), &
                               & fu_expand_environment(fu_content(nlSetup,'nuclide_database_fnm')))

    !small = fu_set_mode(.5e-6, .9e-6, .1e-6, 1.e3)
    !big = fu_set_mode(5.0e-6, 10.0e-6, 2.5e-6, 1.e3)
   
    call set_species(species1, fu_get_material_ptr('passive'), small)
    call set_species(species2, fu_get_material_ptr('passive'), in_gas_phase)
    call set_species(species3, fu_get_material_ptr('generic_aerosol'), big)
    call report(species1)
    call report(species2)
    call report(species3)

    call addSpecies(species_list_ptr, nspecies, (/species1/), 1)
    print *, nspecies
    call addSpecies(species_list_ptr, nspecies, (/species2, species3, species3/), 3)
    print *, nspecies
    print *, ''
    do i = 1, nspecies
      call report(species_list(i))
    end do
    
    call create_chemical_mapper(species_list, mapper, nspecies)
    
    do i = 1, fu_nsubst(mapper)
      do j = 1, fu_nmodes(i, mapper)
        print *, i, j, fu_isp(i,j,1,mapper), fu_nwaves(1,j,mapper)
      end do
    end do

    call msg('1 OK')
    
    ! a stress test
    do i = bignum, 1, -1
      call set_species(species_list(i), &
                     fu_get_material_ptr('generic_aerosol'), &
                     & fu_set_mode(fixed_diameter_flag, &
                                 & i*1.0e-8+1e-9, i*1.0e-8-1e-9, 1.e3))
      if (error) return
    end do
    nspecies = bignum
    call addSpecies(species_list_ptr, nspecies, (/species2/), 1)
    
    call free_chemical_mapper(mapper)
    do i = 1, 100
      !call msg('i = ', i)
      call create_chemical_mapper(species_list, mapper, nSpecies)
      call free_chemical_mapper(mapper)
    end do

    call msg('2 OK')

    call create_chemical_mapper(species_list, mapper, nSpecies)

    do i = 1, fu_nsubst(mapper)
      do j = 1, fu_nmodes(i, mapper)
        call report(species_list(fu_isp(i,j,1,mapper)))
      end do
    end do

    call msg('3 OK')
    
    ! try some wavelengths
    call set_species(species_list(1), fu_get_material_ptr('passive'), small, 550.0e-9)
    call set_species(species_list(2), fu_get_material_ptr('passive'), in_gas_phase, 550.0e-9)
    call set_species(species_list(3), fu_get_material_ptr('generic_aerosol'), big, 550.0e-9)
    call set_species(species_list(4), fu_get_material_ptr('generic_aerosol'), big, 650.0e-9)
    call set_species(species_list(5), fu_get_material_ptr('generic_aerosol'), big, 750.0e-9)
    nspecies = 5
    call set_species(species1, fu_get_material_ptr('passive'), in_gas_phase, 650.0e-9)
    call set_species(species2, fu_get_material_ptr('passive'), small, 950.0e-9)
    call addspecies(species_list_ptr, nspecies, (/species1,species2,species1/), 3)
    
    do i = 1, nspecies
      call report(species_list(i))
    end do
    
    call create_chemical_mapper(species_list, mapper, nspecies)
    
    call msg('4 OK')

    do i = 1, fu_nsubst(mapper)
      do j = 1, fu_nmodes(i, mapper)
        do k = 1, fu_nwaves(i,j,mapper)
          call report(species_list(fu_isp(i,j,k,mapper)))
        end do
      end do
    end do

    ! descriptors

!    chCocktailDescrFile = '/home/vira/silam/silam_v4_6/ini/v4_6_test.cocktails'

    call set_cocktail_description('TEST_1', descr1, ifSpecies)
    call set_cocktail_description('TEST_2', descr1, ifSpecies)
    call set_cocktail_description('TEST_2', descr1, ifSpecies)

  end subroutine chemical_setup_tests


  !********************************************************************************
  
  subroutine create_mode_projection(species_in, nSpecies_in, &
                                  & species_mass, nspecies_mass, &
                                  & species_nbr, nspecies_nbr, &
                                  & ref_list_mass, ref_list_nbr, &
                                  & ifFullSpectrum)
    !
    ! Reprojects the input aerosol size distribution to the given template. 
    ! Creates the "extended adaptor" - the SpeciesReference structure
    !
    implicit none

    ! Imported parameters
    type(silam_species), dimension(:), pointer :: species_in, species_mass, species_nbr
    integer, intent(in) :: nspecies_in, nspecies_mass, nspecies_nbr
    type(TspeciesReference), dimension(:), pointer :: ref_list_mass, ref_list_nbr
    logical, intent(in) :: ifFullSpectrum
!    type(Tmoment_mapping), intent(out) :: mapping_nbr_2_mass

    ! Local variables
    integer :: stat, iSpIn, iSpOut, isubst, isp, modetype, nrefs, iTmp, iSkipCounter
    type(chemical_mapper) :: mapper
    type(Taerosol_mode) :: input_mode, output_mode
    logical :: species_ok, ifSkip
    integer, dimension(:), pointer  :: indices, indSkippedSpeciesModes
    real, dimension(:), pointer :: fractions
    real :: fraction, fraction_total, mass_to_nbr, d1, d2, dens

    real :: fnbr, fmass
    type(Taerosol_mode) :: test_mode

    fractions => fu_work_array()
    indices => fu_work_int_array()
    indSkippedSpeciesModes => fu_work_int_array()
    if(error)return

    call msg('Create_mode_projection: species in')
    do iSp = 1, nSpecies_in
      call report(species_in(iSp))
    end do
    call msg('')
    call msg('Create_mode_projection: species mass')
    do iSp = 1, nspecies_mass
      call report(species_mass(iSp))
    end do
    call msg('')
    call msg('Create_mode_projection: species number')
    do iSp = 1, nspecies_nbr
      call report(species_nbr(iSp))
    end do

    allocate(ref_list_mass(nspecies_in), ref_list_nbr(nspecies_in), stat=stat)
    if (stat /= 0) then
      call set_error('Allocate failed', 'create_mode_projection')
      return
    end if
    
    ! Mass mapping from species_mass to subst-mode pairs
    !
    call create_chemical_mapper(species_mass, mapper, nspecies_mass)
    if(error)return
    !
    ! Grand cycle over the input species mapping them to the mass and number
    !
    do iSpIn = 1, nspecies_in

      input_mode = fu_mode(species_in(iSpIn))
      iSkipCounter = 0  ! will count the output modes skipped from the reprojection

      !
      ! No mapping of number species to mass
      !
      if(fu_name(fu_material(species_in(iSpIn))) == 'nbr_aer')then
        nullify(ref_list_mass(iSpIn)%indSpeciesTo, ref_list_mass(iSpIn)%fract)
        ref_list_mass(iSpIn)%nRefSpecies = 0
      else
!        cycle
!      endif

        !
        ! Detect the way of projection
        !
        if (fu_index(species_in(iSpIn), species_mass, nspecies_mass) > 0) then
          !
          ! Got 1-1 mapping - easy.
          !
          allocate(ref_list_mass(iSpIn)%indSpeciesTo(1), ref_list_mass(iSpIn)%fract(1), stat=stat)
          if (stat /= 0) then
            call set_error('Allocate failed', 'create_mode_projection')
            return
          end if
          ref_list_mass(iSpIn)%nRefSpecies = 1
          ref_list_mass(iSpIn)%indSpeciesTo(1) = fu_index(species_in(iSpIn), species_mass, nspecies_mass)
          ref_list_mass(iSpIn)%fract(1) = 1.
        else
          !
          ! no 1:1 mapping
          !
!        cycle
!      end if

          ! Gas phase stuff must go 1-1, so if not - set error
          !
          if (input_mode == in_gas_phase) then
            call msg('Species to project:')
            call report(species_in(iSpIn))
            call set_error('Cannot project mode for gas phase species', 'create_mode_projection')
            return
          end if
          ! 
          ! No free lunch, need to loop over the relevant modes in species_mass and
          ! find out the fraction of species_in that goes into them
          !
          isubst = fu_isubst(fu_substance_name(species_in(iSpIn)), mapper)
          nrefs = 0
          fraction_total = 0.0

          do iSpOut = 1, fu_nmodes(isubst, mapper)    ! over modes of the iSubst
            isp = fu_isp(isubst, iSpOut, 1, mapper)   ! j-th mode, iSubst-th substance -> mass species index
            output_mode = fu_mode(species_mass(isp))  ! the mass mode
            !
            ! Know how to deal with only some distributions
            !
            if (.not. (output_mode%distr_type == moving_diameter_flag .or. &
                     & output_mode%distr_type == fixed_diameter_flag)) cycle
            !
            ! Solubilities must be the same
            !
            if(.not. (input_mode%solubility == int_missing .or. &
                    & input_mode%solubility == output_mode%solubility))cycle

            d1 = fu_min_d(output_mode)
            d2 = fu_max_d(output_mode)
            !
            ! Fraction of the input mode falling between d1 and d2
            !
            fraction = fu_integrate_volume(d1, d2, input_mode)

            if (fraction > 0.0) then
              nrefs = nrefs + 1 
              fraction_total = fraction_total + fraction
              fractions(nrefs) = fraction
              indices(nrefs) = isp
            else
              iSkipCounter = iSkipCounter + 1
              indSkippedSpeciesModes(iSkipCounter) = iSp
              call msg('Zero mass fraction for mode:')
              call report(output_mode)
            end if

            if (error) return

          end do  ! nModes for specific subst in species_in

          !
          ! If full spectrum is given as projection target, check the fractionation: must be ~1
          !
          if(ifFullSpectrum)then
            if (fraction_total < 0.98 .or. fraction_total > 1.02) then
              call report(species_in(iSpIn))
              call msg('Failed to project mode, total mass fraction = ', fraction_total)
              call set_error('Failed to project mode', 'create_mode_projection')
              return
            end if
          end if

          allocate(ref_list_mass(iSpIn)%indSpeciesTo(nrefs), ref_list_mass(iSpIn)%fract(nrefs), &
                 & stat=stat)
          if (stat /= 0) then
            call set_error('Allocate failed', 'create_mode_projection')
            return
          end if

          ref_list_mass(iSpIn)%nRefSpecies = nrefs
          ref_list_mass(iSpIn)%indSpeciesTo(1:nrefs) = indices(1:nrefs)
          ref_list_mass(iSpIn)%fract(1:nrefs) = fractions(1:nrefs)

        endif  ! if 1:1 mapping found
      endif  ! if species_in of number type

!    end do  ! species_in

!    call free_chemical_mapper(mapper)

      !
      ! Number mapping. Only mass species with static mean mass diameter are taken
      !
      if(nspecies_nbr > 0)then
!        do i = 1, nspecies_in
        !
        ! 1-1 mapping cannot exist: input species must have material,
        ! while species_nbr must not have one. Gas phase species are of
        ! course ignored.
        ! Also, skip the dynamic-diameter species, such as sea salt
        !
        if(input_mode == in_gas_phase .or. fu_ifMeanDiamDynamic(input_mode)) then
          nullify(ref_list_nbr(iSpIn)%indSpeciesTo, ref_list_nbr(iSpIn)%fract)
          ref_list_nbr(iSpIn)%nRefSpecies = 0
          cycle
        end if

        dens = fu_density(species_in(iSpIn)) / &
             & fu_conversion_factor(fu_basic_mass_unit(fu_material(species_in(iSpIn))), &
                                  & 'kg', &
                                  & fu_material(species_in(iSpIn)))
        select case(input_mode%distr_type)
          case (lognormal_flag)
            ! Formula from some strange place ..
            mass_to_nbr = 6.0 / (dens * pi * input_mode%fp1**3.0) * exp(-4.5 * log(input_mode%fp2)**2)
          case (fixed_diameter_flag, moving_diameter_flag)
            mass_to_nbr = 6.0 / (dens * pi * input_mode%mass_mean_d**3.0) 
          case default
            call set_error('Unknown distribution','create_mode_projection')
            return
        end select

        nrefs = 0
        fraction_total = 0.0
        do iSpOut = 1, nspecies_nbr
          output_mode = fu_mode(species_nbr(iSpOut))
          !
          ! Search whether this mode has been zero above
          !
          ifSkip = .false.
          do iTmp = 1, iSkipCounter
            if(output_mode == fu_mode(species_mass(indSkippedSpeciesModes(iTmp))))then
              call msg('Needs skipping the mode:')
              call report(output_mode)
              ifSkip = .true.
              exit
            endif
          end do

          if (.not. (output_mode%distr_type == moving_diameter_flag .or. &
                   & output_mode%distr_type == fixed_diameter_flag)) cycle
          
          if(.not. (input_mode%solubility == int_missing .or. &
                  & input_mode%solubility == output_mode%solubility))cycle

          d1 = fu_min_d(output_mode)
          d2 = fu_max_d(output_mode)

          fraction = fu_integrate_number(d1, d2, input_mode)

          if (fraction > 0.0) then
            !
            ! Mode got something from species_in(iSpIn). But above it
            ! can be skipped for it. Check!
            !
            if(ifSkip)then
              !
              ! Skipped mass-mode has to be forced here
              !
              call msg_warning('Inconsistency: skipped mass mode has non-zero number fraction', &
                             & 'create_mode_projection')
              call report(output_mode)
              call msg('Received fraction to be forced to zero:',fraction)
              fraction = 0.0
            else
              !
              ! Success: both mass and number projections are valid
              !
              nrefs = nrefs + 1 
              fraction_total = fraction_total + fraction
              fractions(nrefs) = fraction
              indices(nrefs) = iSpOut
              call msg('SUCCESS: mode is valid')
              call report(output_mode)
            endif

          else
            !
            ! This mode is skipped here - but it may be not skipped above. Check!
            !
            if(.not. ifSkip)then
              !
              ! Indeed, mass-projection fraction is not zero. It has to be nullified
              !
              call msg_warning('Inconsistency: reasonable mass mode has zero number fraction', &
                             & 'create_mode_projection')
              call report(output_mode)
              !
              ! Find all references to this mode and zero the fractions
              !
              do iTmp = 1, ref_list_mass(iSpIn)%nRefSpecies
                if(fu_mode(species_mass(ref_list_mass(iSpIn)%indSpeciesTo(iTmp))) == output_mode)then
                  call msg('Nullifying reference:',iTmp)
                  ref_list_mass(iSpIn)%fract(iTmp) = 0.0
                endif
              end do
            endif
          end if  ! if fraction > 0

        end do  ! species out (nbr projection)

        if (fraction_total < 0.98 .or. fraction_total > 1.02) then
          call report(species_in(iSpIn))
          call msg('Failed to project mode (number density), total number fraction = ', &
                 & fraction_total)
          call set_error('Failed to project mode', 'create_mode_projection')
          return
        end if

        allocate(ref_list_nbr(iSpIn)%indSpeciesTo(nrefs), ref_list_nbr(iSpIn)%fract(nrefs), stat=stat)
        if (stat /= 0) then
          call set_error('Allocate failed', 'create_mode_projection')
          return
        end if

        ref_list_nbr(iSpIn)%nRefSpecies = nrefs
        ref_list_nbr(iSpIn)%indSpeciesTo(1:nrefs) = indices(1:nrefs)
        ref_list_nbr(iSpIn)%fract(1:nrefs) = mass_to_nbr*fractions(1:nrefs)

      else
        ref_list_nbr(iSpIn)%nRefSpecies = 0
      endif    ! if any number species are to be mapped
    end do  ! species_in

    call free_work_array(fractions)
    call free_work_array(indices)
    call free_work_array(indSkippedSpeciesModes)

 !   test_mode = fu_set_mode(lognormal_flag, 0.038e-6, 1.8)
 !   d1 = 0.3000000E-08
 !   d2 = 0.7663094E-08
 !   fmass = fu_integrate_volume(d1,d2,test_mode)
 !   fnbr = fu_integrate_number(d1,d2,test_mode)
 !   call msg('', (fmass * fu_mean_d(test_mode)**3 /fnbr)**(1./3.))
 !   call msg ('',fnbr, fmass )
 !   call msg('md', fu_mean_d(test_mode))
!
  !  stop
  end subroutine create_mode_projection


  !**************************************************************************************

  function fu_integrate_volume(d1, d2, mode) result(fract)
    !
    ! Computes the fraction of the given aerosol mode that falls between the given boundaries
    !
    implicit none
    real, intent(in) :: d1, d2
    type(Taerosol_mode), intent(in) :: mode

    real :: fract

    real(r8k) :: cumu1, cumu2, mass_median_log, a, b, mean_vol_input, nbr_fract, mfp1, mfp2, dd1, dd2

    select case(mode%distr_type)
      case (lognormal_flag)
        !
        ! Seinfeld & Pandis, 1998: pp.421->, 2006: 358->
        !
        mfp1 = mode%fp1
        mfp2 = mode%fp2
        dd1 = d1
        dd2 = d2
        mass_median_log = log(mfp1) + 3.0*(log(mfp2))**2
!        call msg('Erfect argument:',(log(d1) - mass_median_log) / (sqrt2*log(mode%fp2)))
!        cumu1 = (1.0 + ERF((log(dd1) - mass_median_log) / (dsqrt_2*log(mfp2)))) * 0.5
!        cumu2 = (1.0 + ERF((log(dd2) - mass_median_log) / (dsqrt_2*log(mfp2)))) * 0.5
        cumu1 = (1.0 + erf((log(dd1) - mass_median_log) / (dsqrt_2*log(mfp2)))) * 0.5
        cumu2 = (1.0 + erf((log(dd2) - mass_median_log) / (dsqrt_2*log(mfp2)))) * 0.5
        fract = cumu2 - cumu1

      case (fixed_diameter_flag, moving_diameter_flag)
        ! 
        ! Here...sorry. The bin-type mode can have nonuniform density
        ! (partly) described by the mass mean diameter. However
        ! projecting it correctly becomes rather tedious and the
        ! custom mass-mean diameter is ignored except in the trivial
        ! case were the bins have equal borders (this covers sea salt
        ! emission with ADB).
        !
        if ((d1 .eps. mode%fp1) .and. (d2 .eps. mode%fp2)) then
          fract = 1.0
          return
        end if

        a = max(mode%fp1, d1)
        b = min(mode%fp2, d2)
        
        if (a > b) then
          fract = 0.0
          return
        end if

        nbr_fract = (b-a) / (mode%fp2-mode%fp1)
        !
        ! Mass fraction == nbr fraction * ratio of average volumes
        !
        mean_vol_input = (mode%fp1**2 + mode%fp2**2) * (mode%fp1 + mode%fp2)
        fract = nbr_fract * ((a**2+b**2) * (a+b)) / (mean_vol_input)
        
      case default
        call set_error('Only lognormal and bin modes can be projected', 'fu_integrate_volume')
        
    end select

  end function fu_integrate_volume


  !*******************************************************************************************

  function fu_integrate_number(d1, d2, mode) result(fract)
    !
    ! Computes the fraction of the given aerosol mode that falls between the given boundaries
    !
    implicit none
    real, intent(in) :: d1, d2
    type(Taerosol_mode), intent(in) :: mode

    real(r8k) :: fract
    real(r8k) :: cumu1, cumu2, a, b, median_log, mfp1, mfp2, dd1, dd2, fTmp

    select case(mode%distr_type)
      case (lognormal_flag)
        ! Seinfeld & Pandis, pp. 421 ->
        !
        mfp1 = mode%fp1
        mfp2 = mode%fp2
        dd1 = d1
        dd2 = d2
        median_log = log(mfp1)
!        fTmp = (dlog(dd1) - median_log)
!        fTmp = (dlog(dd1) - median_log) / (dsqrt_2*dlog(mfp2))
!        fTmp = ERF((dlog(dd1) - median_log) / (dsqrt_2*dlog(mfp2)))
!        fTmp = (dlog(dd2) - median_log)
!        fTmp = (dlog(dd2) - median_log) / (dsqrt_2*dlog(mfp2))
!        fTmp = ERF((dlog(dd2) - median_log) / (dsqrt_2*dlog(mfp2)))
!        cumu1 = (1.0 + ERF((log(dd1) - median_log) / (dsqrt_2*log(mfp2)))) * 0.5
!        cumu2 = (1.0 + ERF((log(dd2) - median_log) / (dsqrt_2*log(mfp2)))) * 0.5
        cumu1 = (1.0 + erf((log(dd1) - median_log) / (dsqrt_2*log(mfp2)))) * 0.5
        cumu2 = (1.0 + erf((log(dd2) - median_log) / (dsqrt_2*log(mfp2)))) * 0.5
        fract = cumu2 - cumu1
        
      case (fixed_diameter_flag, moving_diameter_flag)
        if ((d1 .eps. mode%fp1) .and. (d2 .eps. mode%fp2)) then
          fract = 1.0
          return
        end if

        a = max(mode%fp1, d1)
        b = min(mode%fp2, d2)
        
        if (a > b) then
          fract = 0.0
          return
        end if

        fract = (b-a) / (mode%fp2-mode%fp1)

      case default
        call set_error('Only bins and lognormal modes can be projected', 'fu_integrate_number')
        
    end select

  end function fu_integrate_number


  !**************************************************************************************

  real function fu_mass_mean_diam(d1, d2, mode) result(diam)
    !
    ! Computes the fraction of the given aerosol mode that falls between the given boundaries
    !
    implicit none
    
    ! Imported parameters
    real, intent(in) :: d1, d2
    type(Taerosol_mode), intent(in) :: mode

    ! Local variables
    real(r8k) :: cumu1, cumu2, a, b, mean_vol_input, nbr_fract, mfp1, mfp2, dd1, dd2, fTmp, &
            & step, step_2, lgD1, lgD2, fNSteps
    integer :: iStep, nSteps
    real :: DStart, DEnd

    select case(mode%distr_type)
      case (lognormal_flag)
        !
        ! Seinfeld & Pandis, 1998, pp.421->  ; 2006: pp362->
        ! Note: mass distribution is also logarithimic with this median diameters
        !
        ! The trouble is that I have not found easy form for the first moment of the bounded lognormal 
        ! distribution. To save the effort, brute-force integration seems to be the easiest option
        !
        DStart = min(d1,d2)
        lgD1 = log10(DStart)
        lgD2 = log10(max(d1,d2))
        nSteps = max(100, nint(100. * (lgD2 - lgD1)))  ! 100 step per decimal order of diameter size
        fNSteps = nSteps  ! turn to double
        step = 10** ((lgD2 - lgD1) / fNSteps)
        step_2 = sqrt(step)
        cumu1 = 0.
        cumu2 = 0.
        do iStep = 0, nSteps-1
          DEnd = Dstart * step
          fTmp = fu_integrate_volume(Dstart, DEnd, mode)
          cumu1 = cumu1 + fTmp * Dstart * step_2  ! weighted diameter of the step
          cumu2 = cumu2 + fTmp                    ! amount of the step
          Dstart = Dend
        end do
        !
        ! Now careful: if we are far from mean diameter, there can be no contribution from this range:
        ! exponents fall fast.
        !
        if(cumu2 > 1e-20)then
          diam = cumu1 / cumu2
!call msg('Given range, mode diam and sigma, and mass-mean diameter:',(/d1,d2,mode%fp1,mode%fp2,diam/))
        else
          if(min(d1,d2) > mode%fp1)then  ! If no meaningful contribution, can take virtually anything...
            diam = 0.75 * d1 + 0.25 * d2   ! still let's preserve the monotonicity of the distribution
          else
            diam = 0.25 * d1 + 0.75 * d2
          endif
!call msg('CORRECTION, CORRECTION CORRECTION and mass-mean diameter:',(/d1,d2,mode%fp1,mode%fp2,diam/))
        endif
!if(.not. diam > 0)then
!  call msg('*')
!endif
!        !
!        ! volume distribution is lognormal with the same sigma and above median diam:
!        !
!        diam = exp(mass_median_log + 0.5*fTmp)   ! this is for the whole distribution.
        

      case (fixed_diameter_flag, moving_diameter_flag)
        ! 
        ! Here...sorry. The bin-type mode can have nonuniform density
        ! (partly) described by the mass mean diameter. However
        ! projecting it correctly becomes rather tedious and the
        ! custom mass-mean diameter is ignored except in the trivial
        ! case were the bins have equal borders (this covers sea salt
        ! emission with ADB).
        !
        if ((d1 .eps. mode%fp1) .and. (d2 .eps. mode%fp2)) then
          diam = 1.0
          return
        end if

        a = max(mode%fp1, d1)
        b = min(mode%fp2, d2)
        
        if (a > b) then
          diam = 0.0
          return
        end if

        nbr_fract = (b-a) / (mode%fp2-mode%fp1)
        !
        ! Mass fraction == nbr fraction * ratio of average volumes
        !
        mean_vol_input = (mode%fp1**2 + mode%fp2**2) * (mode%fp1 + mode%fp2)
        diam = nbr_fract * ((a**2+b**2) * (a+b)) / (mean_vol_input)
        
      case default
        call set_error('Only lognormal and bin modes can be projected', 'fu_mass_mean_diam')
        
    end select

  end function fu_mass_mean_diam

end module chemical_setup
