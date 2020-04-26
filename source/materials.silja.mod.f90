MODULE materials
  !
  ! This module contains the highest-level handling of materials in the SILAM
  ! model, and should, by itself, only define data structures for general
  ! properties, i.e. such properties that are common to all materials. These
  ! include aerodynamic and deposition properties.
  !
  ! Each physical material in the SILAM model (corresponding to one single
  ! chemical compound) MAY, in addition, have one set of chemistry-specific
  ! properties and contain any number of radionuclides. These additional data
  ! are not contained within the highest-level structure ('silja_material'),
  ! but are rather made targets of pointer variables.

!????????
  ! Should one physical material have alternative general properties, i.e.
  ! appear as particles of differing sizes, it is intended that these
  ! alternatives be defined as different silja_materials, however sharing the
  ! same chemical data and radionuclide pointers, if such exist.
!????????
  ! 
  ! Current (05.02.2003) idea of using materials is the following. Material
  ! is a component of the silam_cocktail containing the key features like 
  ! chemical, radioactive. In radioactive simulations, 
  ! a list of materials in the cocktail corresponds to the source inventory. 
  ! However, radioactive decay creates extra stuff, which is stored in 
  ! decay_chain - a global list of everything available in the model run. 
  ! The same can be done for various chemical materials.
  !
  ! Author: Mikhail Sofiev, FMI mikhail.sofiev@fmi.fi
  !
  ! All units: SI
  ! 
  ! Language: ANSI standard Fortran 90
  ! 
  ! Modules used:
  USE nuclides
!  USE silam_units

  IMPLICIT NONE

  ! The public functions and subroutines available in this module:

  public fu_get_material_ptr
  PUBLIC set_material_missing
  public init_chemical_materials   ! Reads the chemical data file
  public fu_if_chemical_data_loaded  ! check if they are read
  public fu_if_radioactive           ! whether valid radioactive data exist
 ! public fu_type
!  public fu_SILAM_StandardUnit
  public fu_conversion_factor
  public fu_factor_to_basic_unit
  public get_reaction_parameters
  public fu_gas_deposition_param
  public fu_aerosol_param
  public defined
  public report

  public fu_name
  public fu_mole_mass
  public fu_dry_part_density
  public fu_wet_particle_features
  public fu_half_life
  public fu_nuclide
  public fu_mole_to_kg
  public fu_if_gas
  public fu_if_gas_depositing
  public fu_gas_molecular_diffusivity_air
  public fu_gas_molecular_diffusivity_water
  public fu_gas_Rs_default_soil
  public fu_gas_Rs_default_water
  public fu_scavenging_scale_gas_rain
  public fu_scavenging_scale_gas_snow
  public fu_optical_features
  public fu_n_opt_subst
  public fu_names_opt_subst
  public fu_fractions_opt_subst
  public fu_low_mass_threshold
  public fu_basic_mass_unit
  public fu_set_named_amount
  public fu_pollen_taxon_name
  public fu_pollen_material_type
  public fu_pollen_potency_change_rate
  public fu_pollen_break_rate
  public fu_pollen_RH_break_rate
  public fu_pollen_RH_break_thres
  public pollen_break_chem_reacts_ptr
  public set_pollen_taxon_name
  public set_pollen_material_type 

  public fu_deliquescence_humidity
!  ! temporary fix to make this work with passive substance
!  public init_passive_chemicals

  ! The private functions and subroutines not to be used elsewhere:
  ! Encapsulation

  private fu_get_material_ptr_via_name
!  private fu_get_material_ptr_via_sl
  private fu_get_material_ptr_via_index
  private fu_material_defined
  private fu_name_of_material
  private fu_mole_mass_of_material
  private fu_dry_part_dens_of_material
  private fu_half_life_of_material    ! In fact, for nuclide material
  private fu_nuclide_from_material
  private fu_if_gas_material
  ! private fu_material_type
  private get_reaction_param_material
  private copy_material
  private fu_compare_materials_eq
  private fu_molec_diffus_water_gas
  private fu_Rs_default_soil_gas
  private fu_Rs_default_water_gas
  private fu_basic_mass_unit_material
  private report_material
  private report_gas_deposition_data
  private report_aerosol_data
  private report_reaction_data
  private report_pollen_data
  private report_optical_data

  private fu_conversion_factor_material
  private fu_mole_to_kg_by_material
!  private fu_unit_type
  private fu_factor_to_basic_unit_mater
  private fu_check_new_material
  private fu_low_mass_threshold_material
  private fu_del_humid_of_material
  private fu_wet_particle_simple
  private fu_wet_particle_Kelvin
  private test_amount_units

  interface fu_get_material_ptr
!    module procedure fu_get_material_ptr_via_sl
    module procedure fu_get_material_ptr_via_index
    module procedure fu_get_material_ptr_via_name
  end interface

  interface fu_name
    module procedure fu_name_of_material
  end interface

  interface defined
    module procedure fu_material_defined
  end interface
  
  interface fu_mole_mass
    module procedure fu_mole_mass_of_material
  end interface

  interface fu_dry_part_density
    module procedure fu_dry_part_dens_of_material
  end interface

  interface fu_deliquescence_humidity
    module procedure fu_del_humid_of_material
  end interface

  interface fu_half_life
    module procedure fu_half_life_of_material
  end interface

  interface fu_nuclide
    module procedure fu_nuclide_from_material
  end interface

  interface fu_mole_to_kg
    module procedure fu_mole_to_kg_by_material
  end interface

  interface fu_if_gas
    module procedure fu_if_gas_material
  end interface

  interface get_reaction_parameters
    module procedure get_reaction_param_material
  end interface

  interface fu_conversion_factor
    module procedure fu_conversion_factor_material
  end interface

  interface fu_factor_to_basic_unit
    module procedure fu_factor_to_basic_unit_mater
  end interface

  interface fu_basic_mass_unit
    module procedure fu_basic_mass_unit_material
  end interface

  interface assignment(=)
    module procedure copy_material
  end interface

  INTERFACE operator(==)
    MODULE PROCEDURE fu_compare_materials_eq
  END INTERFACE

  interface fu_gas_molecular_diffusivity_water
    module procedure fu_molec_diffus_water_gas
  end interface

  interface fu_gas_Rs_default_soil
    module procedure fu_Rs_default_soil_gas
  end interface

  interface fu_gas_Rs_default_water
    module procedure fu_Rs_default_water_gas
  end interface

  interface report
    module procedure report_material
    module procedure report_gas_deposition_data
    module procedure report_aerosol_data
    module procedure report_reaction_data
    module procedure report_pollen_data
    module procedure report_optical_data
  end interface

  interface fu_low_mass_threshold
     module procedure fu_low_mass_threshold_material
  end interface

  interface fu_wet_particle_features
     module procedure fu_wet_particle_simple
     module procedure fu_wet_particle_Kelvin
  end interface



  !
  ! Some information for deposition for gas-phase species.
  ! Aerosol deposition is treated in cocktail_basic module in a standard way
  !
  TYPE Tgas_deposition_param
    type(silja_logical) :: ifGasDeposited
    real :: Rs_gas_water, &
          & Rs_gas_land_default, &  ! Something suitable as default value if landuse is loosy
          & Rs_gas_bare_land, Rs_gas_urban, Rs_gas_mountain, &
          & Rs_gas_agriculture, Rs_gas_range_land, Rs_gas_range_and_agric, &
          & Rs_gas_decid_forest, Rs_gas_conif_forest, Rs_gas_mixed_forest, &
          & Rs_gas_wetland, &
          & gas_molecular_diffus_in_air, gas_molecular_diffus_in_water, Wesely_f0 
    real :: Henry_const_298K, Henry_const_T_dep, &  ! H(T)=Hmain*exp(H_T_dep*((1/T)-(1/298)))
          & washout_scale_rain_to_standard, washout_scale_snow_to_standard  ! scale from standard SILAM scavenging
  END TYPE Tgas_deposition_param
  type(Tgas_deposition_param) :: gas_depositon_missing = Tgas_deposition_param(silja_undefined, &
                          & real_missing, real_missing, real_missing, real_missing, real_missing, &
                          & real_missing, real_missing, real_missing, real_missing, real_missing, &
                          & real_missing, real_missing, real_missing, real_missing, real_missing, &
                          & real_missing, real_missing, real_missing, real_missing)
  type(Tgas_deposition_param) :: gas_depositon_default = Tgas_deposition_param(silja_true, &
                                                        & 100., 100., 100., 100., 100., & 
                                                        & 100., 100., 100., 100., 100., 100., 100., &
                                                        & 1e-5, 1e-9, 0., 1.2e5, 0., 1., 1.)
  type(Tgas_deposition_param) :: gas_depositon_nodep = Tgas_deposition_param(silja_false, &
                                                  & 1e10, 1e10, 1e10, 1e10, 1e10, & 
                                                  & 1e10, 1e10, 1e10, 1e10, 1e10, 1e10, 1e10,&
                                                  & 1e-5, 1e-9, 0., 0., 0., 0., 0.)

  !
  ! Water soluble particles will grow in high-humidity conditions, which is described
  ! via wet particle density and growth factor for the particle diameter.
  ! SILAM operates with dry particles but deposition, aerosol dynamics and optical stuff
  ! are to be made for wet particles
  !
  type TwetParticle
    real :: fGrowthFactor       ! a ratio with regard to dry diameter
    real :: fWetParticleDensity  ! absolute density, kg/m3
  end type TwetParticle

  type Taerosol_param
    private
    real :: dry_part_density, humidity_growth_alpha, humidity_growth_beta, deliquescence_humidity
  end type Taerosol_param
  !
  ! Default density of the radioactive aerosol
  !
  real, private, parameter :: default_density_nuclide_aeros = 2000.0  ! kg/m3
  !
  ! Standard values for the aerosol parameters
  ! 
  type(Taerosol_param) :: aerosol_param_missing = Taerosol_param( &
                                  & real_missing, real_missing, real_missing, real_missing)
  type(Taerosol_param) :: aerosol_param_soluble = Taerosol_param( &
                                  & default_density_nuclide_aeros, 1.0, 1.4, 0.5)   ! like in sulphates
  type(Taerosol_param) :: aerosol_param_insoluble = Taerosol_param( &
                                  & default_density_nuclide_aeros, 1.0, 1.0, 1.0)  ! deliquescence = 1
  !
  ! Chemical feature of the substance
  !
  ! Reactions will be written after modified Arrhenius equation:
  !
  ! k = beta * T ** alpha * exp(- activation_energy / (gas_constant * T))
  !
  ! Activation energy is measured for reaction of current substance
  ! with one or two agents, so have to make space for them
  !
  type Treaction_param
    character(len=substNmLen) :: reaction_agent1, reaction_agent2
    real :: activation_energy_NPT, alpha, beta          ! NPT = normal pressure & temperature
  end type Treaction_param

  type Treactions
    type(Treaction_param), dimension(:), pointer :: reaction  => null()
  end type Treactions

  
  !  pollen_material
  type Tpollen_material
    ! private
    character(len=substNmLen) :: taxon_nm=''
    integer :: material_type=int_missing
    real :: fRatePotChange=0.0, fBreakRate = 0.0,  fThresRH=real_missing, fBreakRateRH=0.0
    character(len=substNmLen), dimension(:), pointer :: break_agents
    real, dimension(:), pointer :: break_rates
  end type Tpollen_material

  integer, public, parameter :: pollenGrains = 11001  
  integer, public, parameter :: pollenAllergen = 11002  
  integer, public, parameter :: freeAllergen = 11003  

  !
  ! Optical features of the material
  ! If it is a gas, the gas name determines the features
  ! If it is an aerosol, its features are determined using a basis of N known so-far
  ! aerosols. Their names are prescribed. Each aerosol is assumed to be an external mixture
  ! of up to N (so far, 5, jan2008) known substances
  !
  type Toptical_features
    integer :: nOptSubst ! how many standard optical substances are in the material
    character(len=substNmLen), dimension(:), pointer :: chOpticStdSubst   => null() ! (nOptSubst), refers to known "optical" subst
    real, dimension(:), pointer :: fractionOptStdSubst   => null() ! (nOptSubst), composition fraction
  end type Toptical_features


  !----------------------------------------------------------------------
  !
  ! The SILAM main chemical metadata element: silam_material
  !
  TYPE silam_material
    PRIVATE 
    CHARACTER(LEN=substNmLen) :: name
   ! integer                   :: material_type ! the major step from 4.5.1!!
    TYPE(silja_logical)       :: defined
    real                      :: mole_mass
    logical :: ifGasPossible, ifAerosolPossible
    type(Taerosol_param)    :: aerosol_param
    TYPE(Tgas_deposition_param) :: gas_dep_param
    TYPE(Treactions), pointer       :: reaction_data  => null()
    TYPE(silam_nuclide), pointer         :: nuc_data  => null()
    type(Tpollen_material), pointer :: pollen  => null()
    type(Toptical_features), pointer:: optic_data  => null()
    real :: low_mass_thrsh = real_missing
    character(len=substNmLen) :: chLow_mass_thrsh_unit = ''
  END TYPE silam_material

!  type(silam_material), public, parameter :: material_missing = silam_material( &
!                         & '',silja_false,real_missing,.false.,.false., &
!                         & null(),null(),null(),null(),null(),null())


  type silam_material_ptr
    type(silam_material), pointer :: ptr ! => null()
  end type silam_material_ptr

  !---------------------------------------------------------------------
  !----------------------------------------------------------------------
  !
  ! The database type for materials
  !
  type Tchemical_data
    private
    type(silam_material_ptr), dimension(:), pointer :: pChem ! => null()
    integer :: nChemical
    type(silja_logical) :: defined
  end type Tchemical_data
  !
  ! The grand chemical structure of materials
  ! It is complementary to the nuclide datafile nuc_array defined in the nuclides module
  !
  type(Tchemical_data), private, target, save :: ChemicalData

  logical, private, save :: ifChemicalDataLoaded = .false.


  integer, public, parameter :: emission_substance_set = 12001
  integer, public, parameter :: transport_substance_set = 12002

  logical, private, save :: ifRadioactiveSpecies


CONTAINS

  !*****************************************************************

  function fu_get_material_ptr_via_index(indexSubst) result(materialPtr)
    !
    ! All materials are now stored in a single place of the grand dataset 
    ! ChemicalData, this function returns a pointer to the material requested by its
    ! name and type
    !
    implicit none

    ! Return value of this function
    type(silam_material), pointer :: materialPtr

    ! Imported parameters
    integer, intent(in) :: indexSubst

    if(indexSubst >= 1 .and. indexSubst <= ChemicalData%nChemical)then
      materialPtr => ChemicalData%pChem(indexSubst)%ptr
    else
      nullify(materialPtr)
      call msg('Strange substance index:',indexSubst)
      call msg('Number of available substances in the database:',ChemicalData%nChemical)
      call set_error('Strange substance index','fu_get_material_ptr_via_index')
    endif

  end function fu_get_material_ptr_via_index


  !******************************************************************

  function fu_get_material_ptr_via_name(substNm) result(materialPtr)
    implicit none
    character(len=*) :: substNm
    type(silam_material), pointer :: materialPtr
    
    integer :: i
    character(len=len(substNm)) :: substNmUpper

    substNmUpper = fu_str_u_case(substNm)
    do i = 1, chemicalData%nChemical
      !if (chemicalData%pChem(i)%ptr%name == substNmUpper) then
      if (chemicalData%pChem(i)%ptr%name == substNm) then
        materialPtr => chemicalData%pChem(i)%ptr
        return
      end if
    end do
    nullify(materialPtr)

    !!! No errors here 
    call msg_warning('No such chemical:' + substNm, 'get_material_ptr_via_name')

  end function fu_get_material_ptr_via_name
  

  ! ***************************************************************
  
  subroutine set_material_missing(material)
    !
    ! Returns a material with no properties set, but only field 'defined'
    ! set to 'silja_false'.
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE (silam_material) :: material
    !
    material%defined = silja_false
    material%ifAerosolPossible = .false.
    material%ifGasPossible = .false.
    material%low_mass_thrsh = real_missing
    material%aerosol_param = aerosol_param_missing
    material%gas_dep_param = gas_depositon_missing
    nullify(material%reaction_data)
    nullify(material%nuc_data)
    nullify(material%pollen)
  
  END subroutine set_material_missing


  !*****************************************************************

  subroutine init_chemical_materials(chChemicalDataFNm, chNuclideDataFNm) !nlStdSetup) !, nNuclides) !SubstanceList, nNuclides)
    !) !
    ! Reads the chemical database and fills in the appropriate structure.
    ! Note that the radiological features of the materials are taken here from 
    ! the corresponding database. There are just too many of them, so we separated
    ! radioactive and non-radioactive species.
    !
    implicit none

    ! Imported parameter
!    type(Tsilam_namelist), pointer :: nlStdSetup
!    type(TSubstanceList), intent(out) :: SubstanceList
!    integer, intent(in), optional :: nNuclides
    character(len=*), intent(in) :: chChemicalDataFNm, chNuclideDataFNm

    ! Local variables
    type(Tsilam_namelist_group), pointer :: nlGroupSubst, nlGroupNuclides
    type(Tsilam_namelist), pointer :: nlSubst_
    integer :: iTmp, iStat, iSubst, nItems, jTmp
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: pItems
    type(silam_material), pointer :: pSubst_
    logical :: file_exists


!      call msg('======================================TEMPORARY==========================================')
!
!    TYPE(silam_material) :: materialTmp
!    character(len=nuc_name_len) :: chTmp1, chTmp2, chTmp3
!    TYPE (silja_branching), DIMENSION (max_nuclides_used) :: branching
!    INTEGER :: no_of_nucs_used
!    CHARACTER (LEN = nuc_name_len), DIMENSION (max_nuclides_used) :: nuclides_used
!      call msg('======================================TEMPORARY==========================================')



    ! Stupidity check
    !
!    chChemicalDataFNm = fu_expand_environment(fu_content(nlStdSetup,'chemical_database_fnm'))
!    chNuclideDataFNm =  fu_expand_environment(fu_content(nlStdSetup,'nuclide_database_fnm'))

    INQUIRE(file=chChemicalDataFNm, exist=file_exists)
    if(.not. file_exists)then
      call set_error('Chemical data file does not exist:' + chChemicalDataFNm, &
                   & 'init_chemical_materials')
      return
    endif
    call msg("Reading chemical database from :>>>" + chChemicalDataFNm + '<<<')
    
    INQUIRE(file=chNuclideDataFNm, exist=file_exists)

    ! Somehow the Intel runtime thinks that blank-named files exist here:
    if(file_exists .and. chNuclideDataFnm /= '')then
      ifRadioactiveSpecies = .true.
    else
      call msg_warning('Radiological data file does not exist:' + chNuclideDataFNm, &
                     & 'init_chemical_materials')
      ifRadioactiveSpecies = .false.
    endif

    ! Chemical data file is a set of namelists, each describing a single substance
    !
    ! Read the group
    !
    iTmp = fu_next_free_unit()
    open(unit=iTmp, file=chChemicalDataFNm, iostat=iStat, status='old')
    if(iStat == 0)then
      nlGroupSubst => fu_read_namelist_group (iTmp, .false.)
      close(iTmp)
    else
      call set_error('Failed to open chemical data file:' + chChemicalDataFNm, &
                   & 'init_chemical_materials')
      return
    endif

!call msg('')
!call msg('')
!call msg('Material namelist group')
!call write_namelist_group(run_log_funit, nlGroupSubst, .true.)
!call msg('')
!call msg('')

    ! Create the grand chemical structure
    !
!    ifRadioactiveSpecies = .false.
    if(ifRadioactiveSpecies)then
      !
      ! Radioactive species are present: read the whole dataset
      !
      iTmp = fu_next_free_unit()
      open(unit=iTmp, file=chNuclideDataFNm, iostat=iStat, status='old')
      if(iStat == 0)then
        nlGroupNuclides => fu_read_namelist_group (iTmp, .false.)
        close(iTmp)
      else
        call set_error('Failed to open radiological data file:' + chNuclideDataFNm, &
                     & 'init_chemical_materials')
        return
      endif

      allocate(ChemicalData%pChem(fu_nbr_of_namelists(nlGroupSubst) + 30 + &
                                & fu_nbr_of_namelists(nlGroupNuclides)), stat=iStat)
    else
      allocate(ChemicalData%pChem(fu_nbr_of_namelists(nlGroupSubst) + 30), stat=iStat)
    endif
    if(fu_fails(iStat == 0,'Failed to allocate chemical structure','init_chemical_materials'))return

    do iSubst = 1, size(ChemicalData%pChem)
      allocate(ChemicalData%pChem(iSubst)%ptr)
      call set_material_missing(ChemicalData%pChem(iSubst)%ptr)
    end do

    !
    ! Set the chemical materials on-by-one
    !
    do iSubst = 1, fu_nbr_of_namelists(nlGroupSubst)

!call msg('namelist:',iSubst)

      pSubst_ => ChemicalData%pChem(iSubst)%ptr
      nlSubst_ => fu_namelist(nlGroupSubst,iSubst)
!call report(nlSubst_)
      if(.not. associated(nlSubst_))cycle
      if(.not. defined(nlSubst_))cycle

      call set_substance(nlSubst_, pSubst_)
      if(error)return

    end do  ! cycle over substances

    !
    ! Having done the chemicals, let's do the nuclides. They follow the same template but also
    ! include the radiological properties.
    !
    if(ifRadioactiveSpecies)then
      do iSubst = 1, fu_nbr_of_namelists(nlGroupNuclides)

!call msg('nuclide namelist:',iSubst)

        pSubst_ => ChemicalData%pChem(iSubst + fu_nbr_of_namelists(nlGroupSubst))%ptr
        nlSubst_ => fu_namelist(nlGroupNuclides,iSubst)
        if(.not. associated(nlSubst_))cycle
        if(.not. defined(nlSubst_))cycle

        call set_substance(nlSubst_, pSubst_)
        if(error)return

      end do  ! cycle over nuclides
      
      ChemicalData%nChemical = fu_nbr_of_namelists(nlGroupSubst) + fu_nbr_of_namelists(nlGroupNuclides)
      call destroy_namelist_group(nlGroupNuclides)
    else
      ChemicalData%nChemical = fu_nbr_of_namelists(nlGroupSubst)
    endif
    call destroy_namelist_group(nlGroupSubst)

    ChemicalData%defined = silja_true
    ifChemicalDataLoaded = .true.

    call msg('Number of chemical materials set:',ChemicalData%nChemical)
    !
    ! Check for duplicates. Break the run if found something.
    !
    do iSubst = 1, ChemicalData%nChemical
      do iTmp = iSubst+1, ChemicalData%nChemical
        if(ChemicalData%pChem(iSubst)%ptr%name == ChemicalData%pChem(iTmp)%ptr%name)then
          call msg("Substance numbers:", iSubst, iTmp)      
          call set_error('Duplicated substance: '// trim(ChemicalData%pChem(iSubst)%ptr%name), &
                       & 'init_chemical_materials')
          return
        endif
      enddo
!      call report(ChemicalData%pChem(iSubst))
    end do

!    do iSubst = 1, ChemicalData%nChemical
!      call report(ChemicalData%pChem(iSubst)%ptr,.true.)
!    end do

    !call test_amount_units()
    !stop

    CONTAINS

    !==================================================================================================
    
    subroutine set_substance(nlSubst, pSubst)
      !
      ! Sets a single chemical substance from the given namelist
      !
      implicit none
      
      ! Imported parameters
      type(Tsilam_namelist), pointer :: nlSubst
      type(silam_material), pointer :: pSubst
      !
      ! Internal variables
      !
      integer :: iStat, iTmp, nItems
      type(Tsilam_nl_item_ptr), dimension(:), pointer :: pItems
      type(silam_sp) :: sp

      sp%sp => fu_work_string()

      !
      ! Basic parameters
      !
      pSubst%name = fu_content(nlSubst,'chemical_name')
      pSubst%mole_mass = fu_set_named_value(fu_content(nlSubst,'molar_mass'))
      call set_named_value_and_unit(fu_content(nlSubst,'low_concentration_threshold'), &
                                  & pSubst%low_mass_thrsh, &
                                  & pSubst%chLow_mass_thrsh_unit)
      if(error)return

      !RISTO: To get the low_concentration_thresholds to moles/m3 in module midAtm something like this is needed. But this does not work:
      !pSubst%low_mass_thrsh = pSubst%low_mass_thrsh * &
      !                      & fu_factor_to_basic_unit_mater(pSubst%chLow_mass_thrsh_unit, fu_basic_mass_unit(pSubst), pSubst)
      !
      !RISTO TEST: The lines below just illustrates that the low concentration threshold does not work in proper units, but
      !            rather takes the numerical value from silam_chemical.dat and uses the basic unit for the species.
      !chDescrUnit = fu_basic_mass_unit(cocktail_descr) + '/sec'
      !write(*,"(A,3X,A10,3X,F8.6,3X,E12.4,3X,A8,3X,A8,3X,A11)") 'RISTOTEST  : ', &
      !     & trim(pSubst%name), pSubst%mole_mass, pSubst%low_mass_thrsh, trim(pSubst%chLow_mass_thrsh_unit), &
      !     & trim(fu_basic_mass_unit(pSubst)), trim(fu_basic_mass_unit(pSubst)+'/m3')
      !
      !The following attemp to change to correct units does not work either!!!
      !pSubst%low_mass_thrsh = pSubst%low_mass_thrsh * fu_conversion_factor_material(pSubst%chLow_mass_thrsh_unit, fu_basic_mass_unit(pSubst) + '/m3', pSubst)
      !pSubst%chLow_mass_thrsh_unit = fu_basic_mass_unit(pSubst) + '/m3'
      !write(*,"(A,3X,A10,3X,F8.6,3X,E12.4,3X,A8,3X,A8,3X,A11)") 'RISTOTEST2 : ', &
      !     & trim(pSubst%name), pSubst%mole_mass, pSubst%low_mass_thrsh, trim(pSubst%chLow_mass_thrsh_unit), &
      !     & trim(fu_basic_mass_unit(pSubst)), trim(fu_basic_mass_unit(pSubst)+'/m3')
      !fu_conversion_factor_material(chUnitFrom, chUnitTo, material)
      !END RISTO TEST

      !
      ! Somewhat clumsy: gas or aerosol
      !
      if(fu_str_u_case(fu_content(nlSubst,'can_be_gas')) == 'YES')then
        pSubst%ifGasPossible = .true.
      elseif(fu_str_u_case(fu_content(nlSubst,'can_be_gas')) == 'NO')then
        pSubst%ifGasPossible = .false.
      else
        call set_error('the item can_be_gas is not recognised:' + fu_content(nlSubst,'can_be_gas'), &
                     & 'init_chemical_materials')
        return
      endif
      !
      ! Aerosol features
      !
      if(fu_str_u_case(fu_content(nlSubst,'can_be_aerosol')) == 'YES')then
        !
        ! Material can be aerosol (possibly, also a gas)
        !
        pSubst%ifAerosolPossible = .true.
        pSubst%aerosol_param%dry_part_density = fu_content_real(nlSubst,'dry_particle_density_kg_per_m3')
        pSubst%aerosol_param%humidity_growth_alpha = fu_content_real(nlSubst,'humidity_growth_scale_alpha')
        pSubst%aerosol_param%humidity_growth_beta = fu_content_real(nlSubst,'humidity_growth_curvature_beta')
        pSubst%aerosol_param%deliquescence_humidity = fu_content_real(nlSubst,'deliquescence_humidity_0_to_1')

        if(pSubst%aerosol_param%dry_part_density > 1.0e6 .or. &
                                         & pSubst%aerosol_param%dry_part_density <= 0.)then
          call msg('Problem with dry particle density:', pSubst%aerosol_param%dry_part_density)
          call set_error('Problem with dry particle density for substance:' + pSubst%name, &
                       & 'init_chemical_materials')
          return
        endif
        if(pSubst%aerosol_param%humidity_growth_alpha > 1.0e6 .or. &
                                         & pSubst%aerosol_param%humidity_growth_alpha <= 0.)then
          call msg('Problem with humidity growth scale:', pSubst%aerosol_param%humidity_growth_alpha)
          call set_error('Problem with humidity growth scale for substance:' + pSubst%name,  &
                       & 'init_chemical_materials')
          return
        endif
        if(pSubst%aerosol_param%humidity_growth_beta > 1.0e6 .or. &
                                         & pSubst%aerosol_param%humidity_growth_beta <= 0.)then
          call msg('Problem with humidity growth scale:', pSubst%aerosol_param%humidity_growth_alpha)
          call set_error(fu_connect_strings('Problem with dry particle density for substance:', &
                                          & pSubst%name),'init_chemical_materials')
          return
        endif
        if(pSubst%aerosol_param%deliquescence_humidity > 10.0 .or. &
                                         & pSubst%aerosol_param%deliquescence_humidity < 0.)then
          call msg('Problem with deliquescence humidity:', pSubst%aerosol_param%deliquescence_humidity)
          call set_error(fu_connect_strings('Problem with deliquescence humidity for substance:', &
                                         & pSubst%name),'init_chemical_materials')
          return
        endif

      elseif(fu_str_u_case(fu_content(nlSubst,'can_be_aerosol')) == 'NO')then
        !
        ! Gas only
        !
        pSubst%ifAerosolPossible = .false.
        pSubst%aerosol_param = aerosol_param_missing
      else
        !
        ! Information is missing
        !
        call set_error('Item can_be_aerosol must be YES or NO in suibst' + pSubst%name, &
                     & 'init_chemical_materials')
        return
      endif

      !-------------------------------------------------------------------------
      !
      ! Chemical reactions: if none is defined, nullify the pointer
      ! iStat will work as a counter
      !
      nullify(pItems)
      call get_items(nlSubst, 'reaction_first_order', pItems, nItems)
      call get_items(nlSubst, 'reaction_second_order', pItems, iTmp)
      nItems = nItems + iTmp
      call get_items(nlSubst, 'reaction_third_order', pItems, iTmp)
      nItems = nItems + iTmp

      if(nItems > 0)then
        allocate(pSubst%reaction_data, stat=iStat)
        if(iStat /= 0)then
          call set_error('Failed to allocate chemical reaction structure','init_chemical_materials')
          return
        endif
        allocate(pSubst%reaction_data%reaction(nItems), stat=iStat)
        if(iStat /= 0)then
          call set_error('Failed to allocate chemical reactions','init_chemical_materials')
          return
        endif
        !
        ! Now fill-in the reaction parameters
        ! All reactions are treated via k = alpha * T**beta * exp(-Eact/(RT))
        ! First-order reaction has no agents, just these three coefs
        !
        call get_items(nlSubst, 'reaction_first_order', pItems, nItems)
        do iTmp = 1, nItems
          sp%sp = fu_content(pItems(iTmp))
          if(len_trim(sp%sp) < 1)then
            pSubst%reaction_data%reaction(iTmp)%alpha = real_missing
            pSubst%reaction_data%reaction(iTmp)%beta = real_missing
            pSubst%reaction_data%reaction(iTmp)%activation_energy_NPT = real_missing
          else
            read(unit=sp%sp,fmt=*,iostat=iStat)pSubst%reaction_data%reaction(iTmp)%alpha, &
                                             & pSubst%reaction_data%reaction(iTmp)%beta, &
                                             & pSubst%reaction_data%reaction(iTmp)%activation_energy_NPT
          endif
          pSubst%reaction_data%reaction(iTmp)%reaction_agent1 = ''
          pSubst%reaction_data%reaction(iTmp)%reaction_agent2 = ''
        end do ! 1st-order reactions

        ! Second-order has one agent (apart from the current material itself, of course)
        !
        call get_items(nlSubst, 'reaction_second_order', pItems, nItems)
        do iTmp = 1, nItems
          sp%sp = fu_content(pItems(iTmp))
          if(len_trim(sp%sp) < 1)then
            pSubst%reaction_data%reaction(iTmp)%alpha = real_missing
            pSubst%reaction_data%reaction(iTmp)%beta = real_missing
            pSubst%reaction_data%reaction(iTmp)%activation_energy_NPT = real_missing
            pSubst%reaction_data%reaction(iTmp)%reaction_agent1 = ''
          else
            read(unit=sp%sp,fmt=*,iostat=iStat)pSubst%reaction_data%reaction(iTmp)%reaction_agent1, &
                                             & pSubst%reaction_data%reaction(iTmp)%alpha, &
                                             & pSubst%reaction_data%reaction(iTmp)%beta, &
                                             & pSubst%reaction_data%reaction(iTmp)%activation_energy_NPT
          endif
          pSubst%reaction_data%reaction(iTmp)%reaction_agent2 = ''
        end do ! 2d order reactions

        ! Third-order has two agents
        !
        call get_items(nlSubst, 'reaction_third_order', pItems, nItems)
        do iTmp = 1, nItems
          sp%sp = fu_content(pItems(iTmp))
          if(len_trim(sp%sp) < 1)then
            pSubst%reaction_data%reaction(iTmp)%alpha = real_missing
            pSubst%reaction_data%reaction(iTmp)%beta = real_missing
            pSubst%reaction_data%reaction(iTmp)%activation_energy_NPT = real_missing
            pSubst%reaction_data%reaction(iTmp)%reaction_agent1 = ''
            pSubst%reaction_data%reaction(iTmp)%reaction_agent2 = ''
          else
            read(unit=sp%sp,fmt=*,iostat=iStat)pSubst%reaction_data%reaction(iTmp)%reaction_agent1, &
                                             & pSubst%reaction_data%reaction(iTmp)%reaction_agent2, &
                                             & pSubst%reaction_data%reaction(iTmp)%alpha, &
                                             & pSubst%reaction_data%reaction(iTmp)%beta, &
                                             & pSubst%reaction_data%reaction(iTmp)%activation_energy_NPT
          endif
        end do  ! 3-d order reactions

      else
        nullify(pSubst%reaction_data)
      endif ! chemical reactions are defined


      !----------------------------------------------------------------------------------
      !
      ! Deposition features. If none defined, nullify the pointer
      !
      if(pSubst%ifGasPossible)then
        pSubst%gas_dep_param%Rs_gas_water = fu_content_real(nlSubst,'gas_surf_resistance_over_water')
        pSubst%gas_dep_param%Rs_gas_land_default = fu_content_real(nlSubst,'gas_surf_resistance_land_default')  !if landuse is loosy
        pSubst%gas_dep_param%Rs_gas_bare_land= fu_content_real(nlSubst,'gas_surf_resistance_bare_land')
        pSubst%gas_dep_param%Rs_gas_urban = fu_content_real(nlSubst,'gas_surf_resistance_over_urban')
        pSubst%gas_dep_param%Rs_gas_mountain = fu_content_real(nlSubst,'gas_surf_resistance_over_mountain')
        pSubst%gas_dep_param%Rs_gas_agriculture = fu_content_real(nlSubst,'gas_surf_resistance_over_agriculture')
        pSubst%gas_dep_param%Rs_gas_range_land = fu_content_real(nlSubst,'gas_surf_resistance_over_range_land')
        pSubst%gas_dep_param%Rs_gas_range_and_agric = fu_content_real(nlSubst,'gas_surf_resistance_over_range_and_agriculture')
        pSubst%gas_dep_param%Rs_gas_decid_forest = fu_content_real(nlSubst,'gas_surf_resistance_over_deciduous_forest')
        pSubst%gas_dep_param%Rs_gas_conif_forest = fu_content_real(nlSubst,'gas_surf_resistance_over_coniferous_forest')
        pSubst%gas_dep_param%Rs_gas_mixed_forest = fu_content_real(nlSubst,'gas_surf_resistance_over_mixed_forest')
        pSubst%gas_dep_param%Rs_gas_wetland = fu_content_real(nlSubst,'gas_surf_resistance_over_wetland')
        pSubst%gas_dep_param%gas_molecular_diffus_in_air = fu_content_real(nlSubst,'gas_molecular_diffusivity_in_air')
        pSubst%gas_dep_param%gas_molecular_diffus_in_water = fu_content_real(nlSubst,'gas_molecular_diffusivity_in_water')

          ! H(T)=Hmain*exp(H_T_dep*((1/T)-(1/298)))
          !
        pSubst%gas_dep_param%Henry_const_298K = fu_content_real(nlSubst,'Henry_constant_at_298K')
        pSubst%gas_dep_param%Wesely_f0 = fu_content_real(nlSubst,'Wesely_f0')
        pSubst%gas_dep_param%Henry_const_T_dep = fu_content_real(nlSubst,'Henry_constant_T_dependence_exponent_scale')
        pSubst%gas_dep_param%washout_scale_rain_to_standard = &
                         & fu_content_real(nlSubst,'washout_coef_rain_ratio_to_SILAM_standard')
        pSubst%gas_dep_param%washout_scale_snow_to_standard = &
                             & fu_content_real(nlSubst,'washout_coef_snow_ratio_to_SILAM_standard')
        !
        ! Is there anything useful? If not, this is a non-depositing gas
        !
        if(all ((/pSubst%gas_dep_param%Rs_gas_water, pSubst%gas_dep_param%Rs_gas_land_default, &
                & pSubst%gas_dep_param%Rs_gas_bare_land, pSubst%gas_dep_param%Rs_gas_urban, &
                & pSubst%gas_dep_param%Rs_gas_mountain, pSubst%gas_dep_param%Rs_gas_agriculture, &
                & pSubst%gas_dep_param%Rs_gas_range_land, pSubst%gas_dep_param%Rs_gas_range_and_agric, &
                & pSubst%gas_dep_param%Rs_gas_decid_forest, pSubst%gas_dep_param%Rs_gas_conif_forest, &
                & pSubst%gas_dep_param%Rs_gas_mixed_forest, pSubst%gas_dep_param%Rs_gas_wetland, &
                & pSubst%gas_dep_param%gas_molecular_diffus_in_air, &
                & pSubst%gas_dep_param%gas_molecular_diffus_in_water, &
                & pSubst%gas_dep_param%Henry_const_298K, pSubst%gas_dep_param%Wesely_f0, &
                & pSubst%gas_dep_param%Henry_const_T_dep, &
                & pSubst%gas_dep_param%washout_scale_rain_to_standard, &
                & pSubst%gas_dep_param%washout_scale_snow_to_standard /) == real_missing))then
          pSubst%gas_dep_param%ifGasDeposited = silja_false
        else
          pSubst%gas_dep_param%ifGasDeposited = silja_true
        endif
      else
          pSubst%gas_dep_param = gas_depositon_missing
      endif   ! if gas possible

      !-------------------------------------------------------------------------
      !
      ! Radiological data - for radioactive nuclides only
      !
      sp%sp = fu_content(nlSubst,'half_life_period')
      if(len_trim(sp%sp) > 0)then
        allocate(pSubst%nuc_data, stat=iStat)
        if(fu_fails(iStat == 0, 'Failed to allocate nuclide data','init_chemical_materials'))return
        
        call set_nuclide_from_namelist(nlSubst, pSubst%nuc_data) !, sp%sp, sp%sp, sp%sp)
        !
        ! For now, we shall no distinguish the details of deposition features of nuclides
        ! Just check if they are deposible or not
        !
        call msg('Setting default deposition parameters for nuclide: '//trim(fu_nuc_name(pSubst%nuc_data)) &
                        &//' substance: '// fu_name(pSubst))
        ! gas

        if( pSubst%ifGasPossible ) then
          if(fu_content(nlSubst,'gas_deposible') == 'YES')then
            pSubst%gas_dep_param = gas_depositon_default
          elseif(fu_content(nlSubst,'gas_deposible') == 'NO')then
            pSubst%gas_dep_param = gas_depositon_nodep
          else
            call set_error('gas_deposible must be YES or NO for nuclides','init_chemical_materials')
            return
          endif ! can_be_gas
        endif

        ! aerosol
        if(pSubst%ifAerosolPossible)then
          pSubst%aerosol_param = aerosol_param_soluble
        endif 
      else
        nullify(pSubst%nuc_data)
      endif

      !-------------------------------------------------------------------------
      !
      ! Pollen data. Nullify pointer if absent
      !
      sp%sp = fu_content(nlSubst,'pollen_material_type')
      if(len_trim(sp%sp) > 0)then
        allocate(pSubst%pollen, stat=iStat)
        if(iStat /= 0)then
          call set_error('Failed to allocate pollen data','init_chemical_materials')
          return
        endif
        if(fu_str_u_case(sp%sp) == 'POLLEN_GRAINS')then
          pSubst%pollen%material_type = pollenGrains
        elseif(fu_str_u_case(sp%sp) == 'POLLEN_ALLERGEN')then
          pSubst%pollen%material_type = pollenAllergen
        elseif(fu_str_u_case(sp%sp) == 'FREE_ALLERGEN')then
          pSubst%pollen%material_type = freeAllergen
        else
          pSubst%pollen%material_type = int_missing
          call set_error(fu_connect_strings('Unknown pollen material type:',sp%sp),'init_chemical_materials')
          return
        endif
        
        pSubst%pollen%taxon_nm = fu_content(nlSubst,'taxon_name')
        
        if(pSubst%pollen%material_type == pollenGrains)then
          
          if(.not. fu_content(nlSubst,'potency_change_rate')=='') &
            & call setNamedValue(fu_content(nlSubst,'potency_change_rate'), 'kg/sec',  &
                         & pSubst%pollen%fRatePotChange)
          if(.not. fu_content(nlSubst,'pollen_break_rate')=='') &
            & call setNamedValue(fu_content(nlSubst,'pollen_break_rate'), 'fraction/sec',  &
                         & pSubst%pollen%fBreakRate)
          if(.not. fu_content(nlSubst,'pollen_break_RH_rate')=='') &
                      & pSubst%pollen%fBreakRateRH = fu_content_real(nlSubst,'pollen_break_RH_rate')
          if(pSubst%pollen%fBreakRateRH > 0.0)then
            if(.not. fu_content(nlSubst,'pollen_break_RH_threshold')=='')then
              call setNamedValue(fu_content(nlSubst,'pollen_break_RH_threshold'), 'fraction',  &
                         & pSubst%pollen%fThresRH)
            else
              call set_error('Pollen RH break rate given but no RH threshold','init_chemical_materials')
            endif
          endif
          
          call get_items(nlSubst, 'pollen_break_reaction', pItems, nItems)
          if(nItems > 0)then
            allocate(pSubst%pollen%break_agents(nItems), pSubst%pollen%break_rates(nItems), stat=iStat)
            if(iStat /= 0)then
              call set_error('Failed to allocate pollen breaking','init_chemical_materials')
              return
            endif
            do iTmp = 1, nItems
              sp%sp = fu_content(pItems(iTmp))
              read(unit=sp%sp,fmt=*,iostat=iStat)pSubst%pollen%break_agents(nItems), &
                                               & pSubst%pollen%break_rates(nItems)
            end do
          else
            nullify(pSubst%pollen%break_agents, pSubst%pollen%break_rates)
          endif
        endif     
      else
        nullify(pSubst%pollen)
      endif

      !--------------------------------------------------------------------------
      !
      ! Optical data
      !
      sp%sp = fu_content(nlSubst,'number_of_reference_subst_in_mixture')
      if(len_trim(sp%sp) > 0)then
        iTmp = fu_content_int(nlSubst,'number_of_reference_subst_in_mixture')
        allocate(pSubst%optic_data, stat=iStat)
        if(iStat /= 0)then
          call set_error('Failed to allocate optical features','init_chemical_materials')
          return
        endif
        pSubst%optic_data%nOptSubst = iTmp
        allocate(pSubst%optic_data%chOpticStdSubst(pSubst%optic_data%nOptSubst), &
               & pSubst%optic_data%fractionOptStdSubst(pSubst%optic_data%nOptSubst), stat=iStat)
        if(iStat /= 0)then
          call set_error('Failed to allocate optical features arrays','init_chemical_materials')
          return
        endif
        sp%sp = fu_content(nlSubst,'names_of_ref_subst_in_mixture')
        read(unit=sp%sp,fmt=*,iostat=iStat) &
                & (pSubst%optic_data%chOpticStdSubst(jTmp),jTmp = 1,pSubst%optic_data%nOptSubst)
        sp%sp = fu_content(nlSubst,'fractions_of_ref_subst_in_mixture')
        read(unit=sp%sp,fmt=*,iostat=iStat) &
                & (pSubst%optic_data%fractionOptStdSubst(jTmp),jTmp = 1,pSubst%optic_data%nOptSubst)
      else
        nullify(pSubst%optic_data)
      endif

      pSubst%defined = silja_true
      
      call free_work_array(sp%sp)
      
    end subroutine set_substance

  end subroutine init_chemical_materials


  !*********************************************************************************

  logical function fu_check_new_material(material)result(materialOK)
    !
    ! Justs checks the reasonable ranges for the material features
    !
    implicit none

    ! Imported parameter
    type(silam_material), pointer :: material

    materialOK = len_trim(material%name) > 0 
    if(.not.materialOK)then
      call report(material)
      call set_error('Failed name and type check','fu_check_new_material')
      return
    endif

    materialOK = material%aerosol_param%dry_part_density >=0. .and. &
               & material%aerosol_param%dry_part_density <= 1.e10
    if(.not.materialOK)then
      call report(material)
      call msg('Failed aerosol check: funny density', material%aerosol_param%dry_part_density)
      call set_error('Failed aerosol check: funny density','fu_check_new_material')
      return
    endif

    call msg_warning('Material checking is incomplete')

!  TYPE silam_material
!    PRIVATE 
!    CHARACTER(LEN=substNmLen) :: name
!    TYPE(silja_logical)       :: defined
!    real                      :: mole_mass
!    logical :: ifAerosolPossible
!    type(Taerosol_data), pointer         :: aerosol_data => null()
!    TYPE(silam_deposition_data), pointer :: dep_data  => null()
!    TYPE(silam_reactions), pointer       :: reaction_data  => null()
!    TYPE(silam_nuclide), pointer         :: nuc_data  => null()
!    type(silam_pollen_material), pointer :: pollen  => null()
!    type(silam_optical_features), pointer :: optic_data  => null()
!  END TYPE silam_material

  end function fu_check_new_material


  !********************************************************************
  logical function fu_if_chemical_data_loaded()
    implicit none
    fu_if_chemical_data_loaded = ifChemicalDataLoaded
  end function fu_if_chemical_data_loaded


  !*******************************************************************
  logical function fu_if_radioactive(material)
    implicit none
    type(silam_material), intent(in) :: material
    if(associated(material%nuc_data))then
      fu_if_radioactive = defined(material%nuc_data)
    else
      fu_if_radioactive = .false.
    endif
  end function fu_if_radioactive


  !*******************************************************************

  subroutine set_nuc_ptr(material, nuclide_ptr)
    !
    ! Sets the pointer nuc_data to the nuclide. No data is copied, just pointers
    !
    implicit none
    !
    ! Impoerted parameters 
    type(silam_material), intent(inout) :: material
    type(silam_nuclide), pointer :: nuclide_ptr

    material%nuc_data => nuclide_ptr ! Even if it is NULL

  end subroutine set_nuc_ptr


  !***********************************************************************

  function fu_if_gas_material(material) result(ifGas)
    !
    ! Returns silja_true if the given material is a gas and silja_false if it is aerosol, and
    ! silja_undefined if can be both.
    ! Note: it is different from the nuclide_material module function fu_is_it_noble_gas, which 
    ! answers true if the nuclide component is an inert gas.
    ! Also note that the radioactive database is separate from the chemical one, have to check both
    !
    implicit none

    ! result value
    type(silja_logical) :: ifGas
    
    ! Imported parameters 
    type(silam_material), intent(in) :: material

    if(fu_if_noble_gas(material%name))then   ! if known to be a noble gas - all clear
      ifGas = silja_true
      return
    endif

    if(material%ifAerosolPossible .and. material%ifGasPossible)then
      ifGas = silja_undefined
    elseif(material%ifAerosolPossible)then
      ifGas = silja_false
    elseif(material%ifGasPossible)then
      ifGas = silja_true
    else
      call set_error('Material is neither gas nor aerosol:'+material%name,'fu_if_gas_material')
    endif

  end function fu_if_gas_material


  !*******************************************************************

  subroutine get_reaction_param_material(material, agent1, agent2, activation_energy, alpha, beta)
    !
    ! Returns three parameters for the requested reaction
    !
    implicit none

    ! Imported parameters
    type(silam_material), intent(in) :: material
    character(len=*), intent(in) :: agent1, agent2
    real, intent(out) :: activation_energy, alpha, beta

    ! Local variables
    integer :: iReaction

    ! Scans the reactions looking for the right one. Check the agents one-by-one
    !
    do iReaction = 1, size(material%reaction_data%reaction)
      if(index(material%reaction_data%reaction(iReaction)%reaction_agent1, agent1) == 1)then
        if(len_trim(agent2) > 0)then
          if(index(material%reaction_data%reaction(iReaction)%reaction_agent2, agent2) /= 1)cycle
        endif  ! if second agent is needed
        activation_energy = material%reaction_data%reaction(iReaction)%activation_energy_NPT
        alpha = material%reaction_data%reaction(iReaction)%alpha
        beta = material%reaction_data%reaction(iReaction)%beta
        return
      endif  ! if first agent is found
    end do  ! cycle over reactions

    call set_error('Cannot find reaction for agents:'+agent1+','+agent2,'get_reaction_param_material')

  end subroutine get_reaction_param_material


  !***********************************************************************

  function fu_wet_particle_simple(material, fRelHumidity)result(wetParticle)
    !
    ! SILAM operates with dry particles to avoid dynamic size. But for deposition, 
    ! aerosol dynamics, and optical computations wet parameter is needed.
    ! We take into account the deliquescence relative humidity (without hysteresis)
    ! and a possibility of RH >= 1. See Notebook 10, pp.50-51 for details.
    ! For simple fast computation, Kelvin effect is ignored
    !
    implicit none

    ! return structure
    type(TwetParticle) :: wetParticle

    ! Imported parameters
    type(silam_material), intent(in) :: material
    real, intent(in) :: fRelHumidity

    ! Local variable
    real :: fTmp

    !
    ! Growth is defined via three parameters: deliquescence RH point [0..1], 
    ! growth scaling alpha [unitless] and growth curvature beta [unitless]
    !
    if(material%ifAerosolPossible)then
      !
      ! Aerosol, possibly, soluble
      !
      if(fRelHumidity < material%aerosol_param%Deliquescence_Humidity)then
        wetParticle%fGrowthFactor = 1.0
        wetParticle%fWetParticleDensity = material%aerosol_param%dry_part_density
      elseif(fRelHumidity <= 0.98)then
        wetParticle%fGrowthFactor = material%aerosol_param%humidity_growth_alpha * &
                                  & ((material%aerosol_param%humidity_growth_beta - fRelHumidity) / &
                                   & (1.0 - fRelHumidity))** 0.333333333333
        fTmp = 1 / (wetParticle%fGrowthFactor * wetParticle%fGrowthFactor * wetParticle%fGrowthFactor)
        wetParticle%fWetParticleDensity = material%aerosol_param%dry_part_density * fTmp + &  ! core of the particle
                                        & 1000.0 * (1. - fTmp)   ! water fraction in particle volume
      elseif(fRelHumidity <= 1.03)then
!        wetParticle%fGrowthFactor = material%aerosol_param%humidity_growth_alpha * &
!         & (4.64 * (material%aerosol_param%humidity_growth_beta - 0.99)**0.3333333333333 + &
!         & 155 * (fRelHumidity - 0.99) * (material%aerosol_param%humidity_growth_beta - 1) / &
!                                       & ((material%aerosol_param%humidity_growth_beta - 0.99)**0.6666666667))
       
        ! after 98% relative humidity take just linear growth
        wetParticle%fGrowthFactor = material%aerosol_param%humidity_growth_alpha * &
         & (3.684 * (material%aerosol_param%humidity_growth_beta - 0.98)**0.3333333333333 + &
         & 61.4 * (fRelHumidity - 0.98) * (material%aerosol_param%humidity_growth_beta - 1) / &
                                       & ((material%aerosol_param%humidity_growth_beta - 0.98)**0.6666666667))
        wetParticle%fWetParticleDensity = 1000.0  ! just water density
      else
       ! Do not allow relative humidities bigger than 103%
        wetParticle%fGrowthFactor = material%aerosol_param%humidity_growth_alpha * &
         & (3.684 * (material%aerosol_param%humidity_growth_beta - 0.98)**0.3333333333333 + &
         & 8.6 * (material%aerosol_param%humidity_growth_beta - 1) / &
                                       & ((material%aerosol_param%humidity_growth_beta - 0.98)**0.6666666667))
        wetParticle%fWetParticleDensity = 1000.0  ! just water density

      endif
    else
      !
      ! Gas
      !
      wetParticle%fGrowthFactor = 1.0
      wetParticle%fWetParticleDensity = 0.0
    endif

  end function fu_wet_particle_simple

  
  !****************************************************************************

  function fu_wet_particle_Kelvin(material, Ddry, fTemperature, fRelHumidity)result(wetParticle)
    !
    ! SILAM operates with dry particles to avoid dynamic size. But for deposition, 
    ! aerosol dynamics, and optical computations wet parameter is needed.
    ! We take into account the deliquescence relative humidity (without hysteresis)
    ! and a possibility of RH >= 1. See Notebook 10, pp.50-51 for details.
    ! Kelvin effect is computed using simple bisectional algorithm.
    !
    ! !!! Uses iterational algorithm, so probably slow and not to be used too often !!!
    ! 
    implicit none

    ! return structure
    type(TwetParticle) :: wetParticle

    ! Imported parameters
    type(silam_material), intent(in) :: material
    real, intent(in) :: Ddry, fTemperature, fRelHumidity

    ! Local variable
    real :: fTmp1, fTmp2

    !
    ! Growth is defined via three parameters: deliquescence RH point [0..1], 
    ! growth scaling alpha [unitless] and growth curvature beta [unitless]
    !
    if(.not. material%ifAerosolPossible)then
      !
      ! Gas
      !
      wetParticle%fGrowthFactor = 1.0
      wetParticle%fWetParticleDensity = 0.0
      return
    else
      !
      ! Aerosol, possibly, soluble
      !
      if(fRelHumidity < material%aerosol_param%Deliquescence_Humidity)then
        wetParticle%fGrowthFactor = 1.0
        wetParticle%fWetParticleDensity = material%aerosol_param%dry_part_density
        return
      elseif(fRelHumidity <= 0.99)then

        wetParticle%fGrowthFactor = fu_wetSize(Ddry, &
                 & material%aerosol_param%humidity_growth_alpha, &
                 & material%aerosol_param%humidity_growth_beta, &
                 & material%mole_mass / (1000. * avogadro * material%aerosol_param%dry_part_density), &
                 & fTemperature, fRelHumidity) / Ddry

      else
        ! after 98% relative humidity take just linear growth

        fTmp1 = fu_wetSize(Ddry, &
                 & material%aerosol_param%humidity_growth_alpha, &
                 & material%aerosol_param%humidity_growth_beta, &
                 & material%mole_mass / (1000. * avogadro * material%aerosol_param%dry_part_density), &
                 & fTemperature, 0.98) / Ddry
        fTmp2 = fu_wetSize(Ddry, &
                 & material%aerosol_param%humidity_growth_alpha, &
                 & material%aerosol_param%humidity_growth_beta, &
                 & material%mole_mass / (1000. * avogadro * material%aerosol_param%dry_part_density), &
                 & fTemperature, 0.99) / Ddry

        if(fRelHumidity <= 1.03)then 
          wetParticle%fGrowthFactor = fTmp2 + (fTmp2 - fTmp1) * (fRelHumidity - 0.99) / 0.01
        else ! Do not allow relative humidities bigger than 103%
          wetParticle%fGrowthFactor = fTmp2 + (fTmp2 - fTmp1) * (1.03 - 0.99) / 0.01
        endif
   
      endif
      fTmp1 = 1 / (wetParticle%fGrowthFactor * wetParticle%fGrowthFactor * wetParticle%fGrowthFactor)
      wetParticle%fWetParticleDensity = material%aerosol_param%dry_part_density * fTmp1 + &  ! core of the particle
                                        & 1000.0 * (1. - fTmp1)   ! water fraction in particle volume

    endif


    contains

      real function fu_wetSize(Ddry, alpha, beta, Vmolec, T, RH)

        !
        ! Solve the equation f(x) = a*((b*exp(x)-RH)/(exp(x)-RH))*Ddry - Dcrit/x = 0 
        ! where x = Dcrit/Dwet
        !
        implicit none
        
        real :: Ddry, alpha, beta, Vmolec, T, RH
        real :: Dcrit, hi, lo, mid, epsilon
        
        !
	    ! Critical diameter
		Dcrit = 2. * srfTns_water * Vmolec / (boltzmann_const * T)

		! Bisection algorithm:
        ! Set fu(lo) < 0 and fu(hi) > 0; epsilon small

        epsilon = 0.00001 * (Dcrit / Ddry) ! a small number
        lo = epsilon                       ! a small number
        hi = min(1. / lo, 700.)            ! a large number
	 
	    if(fu_f(lo, alpha, beta, RH, Ddry, Dcrit) > 0. .or. fu_f(hi, alpha, beta, RH, Ddry, Dcrit) < 0.)then
		  call set_error('Starting values fail','fu_wet_particle_Kelvin')
		  call msg('hi, lo: ', hi, lo)
		  return
        endif

        mid = lo + (hi-lo)/2
        do while ((mid - lo > epsilon) .and. (hi - mid > epsilon)) 
          if (fu_f(mid, alpha, beta, RH, Ddry, Dcrit) <= 0) then
            lo = mid
          else
            hi = mid
          endif
          mid = lo + (hi-lo)/2
        enddo
 
        ! Dwet
        fu_wetSize = Dcrit / mid

      end function fu_wetSize


	  real function fu_f(x, a, b, RH, Dd, Dc)

	    implicit none
		real :: x, b, a, RH, Dd, Dc

		fu_f = a * Dd / Dc * ((b * exp(x) - RH) / (exp(x) - RH))**(1./3.)  -  1. / x

	  end function fu_f


  end function fu_wet_particle_Kelvin


  !*********************************************************************

  subroutine report_material(material, ifDetailed)
    !
    ! Prints the report about the material - in a form of namelist
    !
    implicit none
    type(silam_material), intent(in) :: material
    logical, intent(in), optional :: ifDetailed
    type(silam_sp) :: sp

    sp%sp => fu_work_string()
    if(error)return

    if(material%defined == silja_true)then
      call msg('#------------------- SILAM material report ---------------------')
      call write_namelist_item('LIST',material%name)
      call write_namelist_item('chemical_name',material%name)
      call write_namelist_item('molar_mass',material%mole_mass)
      write(unit=sp%sp,fmt='(E9.3,1x,A10)')material%low_mass_thrsh, trim(material%chLow_mass_thrsh_unit)
      call write_namelist_item('low_concentration_threshold',sp%sp)
      
      if(material%ifAerosolPossible)then
        call write_namelist_item('can_be_aerosol','YES')
        call write_namelist_item('dry_particle_density_kg_per_m3', &
                               & material%aerosol_param%dry_part_density)
      else
        call write_namelist_item('can_be_aerosol','NO')
      endif

      if(material%ifGasPossible)then
        call write_namelist_item('can_be_gas','YES')
      else
        call write_namelist_item('can_be_gas','NO')
      endif

      if(present(ifDetailed))then 
        if(ifDetailed)then
          call report(material%aerosol_param)
          call report(material%gas_dep_param)
          if(associated(material%reaction_data)) call report(material%reaction_data)
          if(associated(material%nuc_data)) call report(material%nuc_data)
          if(associated(material%pollen)) call report(material%pollen)
          if(associated(material%optic_data)) call report(material%optic_data)
        endif
      endif
    else
      call msg('#------------------- Undefined material ---------------------')
    endif

    call write_namelist_item('END_LIST',material%name)
    call msg('#------------------- End material report ---------------------')

    call free_work_array(sp%sp)

  end subroutine report_material


  !******************************************************************

  subroutine report_gas_deposition_data(gdd)
    !
    ! Writes the parameters of the gaseous deposition parameterization
    !
    implicit none
    
    ! Imported parameters
    type(Tgas_deposition_param) :: gdd

    if(fu_true(gdd%ifGasDeposited))then
        call write_namelist_item('gas_deposible','NO')
    else
        call write_namelist_item('gas_deposible','YES')
    endif
    call write_namelist_item('gas_surf_resistance_over_water', gdd%Rs_gas_water)
    call write_namelist_item('gas_surf_resistance_land_default', gdd%Rs_gas_land_default)
    call write_namelist_item('gas_surf_resistance_bare_land', gdd%Rs_gas_bare_land)
    call write_namelist_item('gas_surf_resistance_over_urban', gdd%Rs_gas_urban)
    call write_namelist_item('gas_surf_resistance_over_mountain', gdd%Rs_gas_mountain)
    call write_namelist_item('gas_surf_resistance_over_agriculture', gdd%Rs_gas_agriculture)
    call write_namelist_item('gas_surf_resistance_over_range_land', gdd%Rs_gas_range_land) 
    call write_namelist_item('gas_surf_resistance_over_range_and_agriculture', gdd%Rs_gas_range_and_agric)
    call write_namelist_item('gas_surf_resistance_over_deciduous_forest', gdd%Rs_gas_decid_forest)
    call write_namelist_item('gas_surf_resistance_over_coniferous_forest', gdd%Rs_gas_conif_forest)
    call write_namelist_item('gas_surf_resistance_over_mixed_forest', gdd%Rs_gas_mixed_forest)
    call write_namelist_item('gas_surf_resistance_over_wetland', gdd%Rs_gas_wetland)
    call write_namelist_item('gas_molecular_diffusivity_in_air', gdd%gas_molecular_diffus_in_air)
    call write_namelist_item('gas_molecular_diffusivity_in_water', gdd%gas_molecular_diffus_in_water)
    call write_namelist_item('Henry_constant_at_298K', gdd%Henry_const_298K)
    call write_namelist_item('Henry_constant_T_dependence_exponent_scale', gdd%Henry_const_T_dep)
    call write_namelist_item('Wesely_f0', gdd%Wesely_f0)
    call write_namelist_item('washout_coef_rain_ratio_to_SILAM_standard', gdd%washout_scale_rain_to_standard)
    call write_namelist_item('washout_coef_snow_ratio_to_SILAM_standard', gdd%washout_scale_snow_to_standard)  

  end subroutine report_gas_deposition_data


  !******************************************************************

  subroutine report_aerosol_data(ad)
    !
    ! Writes the parameters of the aerosol parameterization
    !
    implicit none
    
    ! Imported parameters
    type(Taerosol_param) :: ad
    
    call write_namelist_item('dry_particle_density_kg_per_m3', ad%dry_part_density)
    call write_namelist_item('humidity_growth_scale_alpha', ad%humidity_growth_alpha)
    call write_namelist_item('humidity_growth_curvature_beta', ad%humidity_growth_beta)
    call write_namelist_item('deliquescence_humidity_0_to_1', ad%deliquescence_humidity)
    
  end subroutine report_aerosol_data


  !******************************************************************

  subroutine report_reaction_data(rd)
    !
    ! Writes the parameters of the gaseous deposition parameterization
    !
    implicit none
    
    ! Imported parameters
    type(Treactions), pointer :: rd

    ! Local variables
    integer :: iR
    type(silam_sp) :: sp

    call write_namelist_item('number_of_reactions',size(rd%reaction))
    if(size(rd%reaction) > 0)then
      sp%sp => fu_work_string()
      if(error)return
      do iR = 1, size(rd%reaction)
        if(len_trim(rd%reaction(iR)%reaction_agent1) < 1)then ! No agents: 1-st order
          write(unit=sp%sp,fmt='(3(E10.5,1x))')rd%reaction(iR)%alpha, &
                                             & rd%reaction(iR)%beta, &
                                             & rd%reaction(iR)%activation_energy_NPT
          call write_namelist_item('reaction_first_order',sp%sp)
        elseif(len_trim(rd%reaction(iR)%reaction_agent2) < 1)then ! only one agent: 2-nd order
          write(unit=sp%sp,fmt='(A,1x,3(E10.5,1x))')rd%reaction(iR)%reaction_agent1, &
                                                  & rd%reaction(iR)%alpha, &
                                                  & rd%reaction(iR)%beta, &
                                                  & rd%reaction(iR)%activation_energy_NPT
          call write_namelist_item('reaction_second_order',sp%sp)
        else                                                      ! all exists: 3rd order
          write(unit=sp%sp,fmt='(2(A,1x),3(E10.5,1x))')rd%reaction(iR)%reaction_agent1, &
                                                     & rd%reaction(iR)%reaction_agent2, &
                                                     & rd%reaction(iR)%alpha, &
                                                     & rd%reaction(iR)%beta, &
                                                     & rd%reaction(iR)%activation_energy_NPT
          call write_namelist_item('reaction_third_order',sp%sp)
        endif        
      end do
      call free_work_array(sp%sp)
    end if

  end subroutine report_reaction_data


  !******************************************************************

  subroutine report_pollen_data(pd)
    !
    ! Writes the pollen features of the material
    !
    implicit none
    
    ! Imported parameters
    type(Tpollen_material), pointer :: pd

    call write_namelist_item('taxon', pd%taxon_nm)
    call write_namelist_item('material_type',pd%material_type)
    call write_namelist_item('potency_change_rate', pd%fRatePotChange)
    call write_namelist_item('break_rate', pd%fBreakRate)
    call write_namelist_item('RH_threshold', pd%fThresRH)
    call write_namelist_item('pollen_break_rate_RH',pd%fBreakRateRH)
    
  end subroutine report_pollen_data


  !******************************************************************

  subroutine report_optical_data(od)
    !
    ! Writes the parameters of the gaseous deposition parameterization
    !
    implicit none
    
    ! Imported parameters
    type(Toptical_features), pointer :: od

    ! Local variables
    type(silam_sp) :: sp
    integer :: iS

    sp%sp => fu_work_string()
    if(error)return
    
    call write_namelist_item('number_of_reference_subst_in_mixture',od%nOptSubst)
    
    write(unit=sp%sp,fmt='(100(A,1x))')(trim(od%chOpticStdSubst(iS)), iS=1,od%nOptSubst)
    call write_namelist_item('names_of_ref_subst_in_mixture',sp%sp)

    write(unit=sp%sp,fmt='(100(F6.4,1x))')(od%fractionOptStdSubst(iS), iS=1,od%nOptSubst)
    call write_namelist_item('fractions_of_ref_subst_in_mixture',sp%sp)
    
    call free_work_array(sp%sp)

  end subroutine report_optical_data
  

  !*****************************************************************

  function fu_basic_mass_unit_material(material) result(unit)
    !
    ! This function will return the basic unit (as used in transport mass map) for 
    ! the given material. Have to invove some guesses:
    ! - radioactive materials are in Bq
    ! - known chemicals with useful mole mass are in moles
    ! - known pollen is in numbers
    ! - the rest is in kg
    !
    implicit none

    type(silam_material), intent(in) :: material
    
    character(len=10) :: unit

    if(associated(material%nuc_data))then  ! A known radioactive material
      if(defined(material%nuc_data))then
        unit = 'Bq'
      else
        call set_error('Nuclide data are associated but undefined for:'+material%name, &
                     & 'fu_basic_mass_unit_material')
      endif
    elseif(material%mole_mass > 0.0)then  ! known chemical
      unit = 'mole'
    elseif(associated(material%pollen))then   ! only pollen is considered as numbers
      if(material%pollen%material_type == pollenGrains)then
        unit = 'number'
      else
        unit = 'kg'
      endif
    else
      unit = 'kg'
    endif

  end function fu_basic_mass_unit_material


  !***********************************************************************

  real function fu_set_named_amount(chInputLine, material)
    !
     ! Get a line, such as "10 kg" and a new unit, such as "g" and returns
     ! the given amount in the given unit. Should material is mandatory for
     ! the unit conversion, it has to be given as well
     !
     implicit none
 
     ! Imported parameters
     character(len=*), intent(in) :: chInputLine  ! number and unit
     type(silam_material), intent(in), optional :: material
 
     ! Local variables
     real :: fIn
     integer :: iStat
     character(len=10) :: chGivenUnit
     !
     ! Just read the line and, find out the conversion factor and return the result
     !
     fu_set_named_amount = real_missing
     read(unit=chInputLine,fmt=*,iostat=iStat) fIn, chGivenUnit
     if(iStat /= 0)then
       call set_error(fu_connect_strings('Failed the line:',chInputLine),'fu_set_named_amount')
       return
     endif
 
     if(present(material))then
       fu_set_named_amount = fIn *  fu_conversion_factor(chGivenUnit, fu_basic_mass_unit(material))
     else
       fu_set_named_amount = fIn
     end if
 
   end function fu_set_named_amount



  !***********************************************************************
  !
  ! ENCAPSULATION OF MATERIAL FEATURES
  !
  !***********************************************************************

  ! Pollen material
  !=======================================================================
  function fu_pollen_taxon_name(material) result(name)
    implicit none
    character(len = substNmLen) :: name
    type(silam_material), intent(in) :: material
    if(.not. associated(material%pollen))then
      name = ''
    else
      name = material%pollen%taxon_nm
    endif
  end function fu_pollen_taxon_name
  !=======================================================================
  integer function fu_pollen_material_type(material)
    implicit none
    type(silam_material), intent(in) :: material
    if(.not. associated(material%pollen))then
      fu_pollen_material_type = int_missing
    else
      fu_pollen_material_type = material%pollen%material_type
    endif
  end function fu_pollen_material_type
  !=======================================================================  
  real function fu_pollen_potency_change_rate(material)
    implicit none
    type(silam_material), intent(in) :: material
    if(.not. associated(material%pollen))then
      fu_pollen_potency_change_rate = 0.0
    else
      fu_pollen_potency_change_rate = material%pollen%fRatePotChange
    endif
  end function fu_pollen_potency_change_rate
  !=======================================================================  
  real function fu_pollen_break_rate(material)
    implicit none
    type(silam_material), intent(in) :: material
    if(.not. associated(material%pollen))then
      fu_pollen_break_rate = 0.0
    else
      fu_pollen_break_rate = material%pollen%fBreakRate
    endif
  end function fu_pollen_break_rate
  !=======================================================================  
  real function fu_pollen_RH_break_rate(material)
    implicit none
    type(silam_material), intent(in) :: material
    if(.not. associated(material%pollen))then
      fu_pollen_RH_break_rate = 0.0
    else
      fu_pollen_RH_break_rate = material%pollen%fBreakRateRH
    endif
  end function fu_pollen_RH_break_rate
  !=======================================================================  
  real function fu_pollen_RH_break_thres(material)
    implicit none
    type(silam_material), intent(in) :: material
    if(.not. associated(material%pollen))then
      fu_pollen_RH_break_thres = real_missing
    else
      fu_pollen_RH_break_thres = material%pollen%fThresRH
    endif
  end function fu_pollen_RH_break_thres
  !=======================================================================  
  integer function fu_n_pollen_break_chem_reacts(material)
    implicit none
    type(silam_material), intent(in) :: material
    if(.not. associated(material%pollen))then
      fu_n_pollen_break_chem_reacts = 0.0
      return
    endif
    if(.not. associated(material%pollen%break_agents))then
      fu_n_pollen_break_chem_reacts = 0.0
      return
    endif
    fu_n_pollen_break_chem_reacts = len(material%pollen%break_agents)
  end function fu_n_pollen_break_chem_reacts
  !=======================================================================  
  subroutine pollen_break_chem_reacts_ptr(material, agents, rates, nReacts)
    implicit none
    type(silam_material), intent(in) :: material
    character(len=substNmLen), dimension(:), pointer :: agents
    real, dimension(:), pointer :: rates
    integer, intent(out) :: nReacts
    if(.not. associated(material%pollen))then
      nReacts = 0
      return
    endif
    if(.not. associated(material%pollen%break_agents))then
      nReacts = 0
      return
    endif
    nReacts = len(material%pollen%break_agents)
    agents => material%pollen%break_agents
    rates => material%pollen%break_rates
  end subroutine pollen_break_chem_reacts_ptr
  
  !=======================================================================
  subroutine set_pollen_taxon_name(material, name)
    type(silam_material), intent(inout) :: material
    character(len = substNmLen), intent(in) :: name
    if(.not. associated(material%pollen))then
        call set_error('pollen not associated','set_pollen_taxon_name')
    endif
    material%pollen%taxon_nm = name
  end subroutine set_pollen_taxon_name
  !=======================================================================
  subroutine set_pollen_material_type(material, material_type)
    type(silam_material), intent(inout) :: material
    integer, intent(in) :: material_type
    if(.not. associated(material%pollen))then
        call set_error('pollen not associated','set_pollen_material_type')
    endif
    select case(material_type)
      case(pollenGrains, pollenAllergen, freeAllergen)
        material%pollen%material_type = material_type
      case default
        call set_error('Unknown pollen material type','set_pollen_material_type')
    end select
  end subroutine set_pollen_material_type
  !=======================================================================
  
  
  !
  !=======================================================================
  real function fu_mole_mass_of_material(material) result (mass)
    implicit none
    type(silam_material), intent(in) :: material
    mass = material%mole_mass
  end function fu_mole_mass_of_material

  !=======================================================================
  real function fu_dry_part_dens_of_material(material) result (dens)
    implicit none
    type(silam_material), intent(in) :: material
    dens = material%aerosol_param%dry_part_density
  end function fu_dry_part_dens_of_material
  !=======================================================================
  real function fu_del_humid_of_material(material) result (del_humid)
    implicit none
    type(silam_material), intent(in) :: material
    del_humid = material%aerosol_param%deliquescence_humidity
  end function fu_del_humid_of_material
  !=======================================================================
  real function fu_half_life_of_material(material) 
    implicit none
    type(silam_material), intent(in) :: material
    if(.not.associated(material%nuc_data))then
      call set_error('Non-radioactive material','fu_half_life_of_material')
    else
      fu_half_life_of_material = fu_half_life(material%nuc_data)
    endif
  end function fu_half_life_of_material

  !=========================================================================
  function fu_optical_features(material) result(features)
    implicit none
    type(Toptical_features), pointer :: features
    type(silam_material), intent(in) :: material
    features => material%optic_data
  end function fu_optical_features

  !=======================================================================
  integer function fu_n_opt_subst(material)
    implicit none
    type(silam_material), intent(in) :: material
    if (associated(material%optic_data)) then
      fu_n_opt_subst = material%optic_data%nOptSubst
    else
      fu_n_opt_subst = 0
    end if
  end function fu_n_opt_subst

  !=======================================================================
  function fu_names_opt_subst(material) result(names)
    implicit none
    type(Toptical_features), pointer :: features
    type(silam_material), intent(in) :: material
    character(len=substNmLen), dimension(material%optic_data%nOptSubst) :: names
    features => material%optic_data
    names = features%chOpticStdSubst
  end function fu_names_opt_subst

  !=======================================================================
  function fu_fractions_opt_subst(material) result(fractions)
    implicit none
    type(Toptical_features), pointer :: features
    type(silam_material), intent(in) :: material
    real, dimension(material%optic_data%nOptSubst) :: fractions
    features => material%optic_data
    fractions = features%fractionOptStdSubst
  end function fu_fractions_opt_subst

  !=======================================================================


  subroutine copy_material(material_out, material_in)
    !
    ! This subroutine overloads the assignment operator = . Necessity is,
    ! of course, dictated by the presence of pointers in the 
    ! material structure definition
    !
    implicit none

    ! Imported parameters
    type(silam_material), intent(in) :: material_in
    type(silam_material), intent(out) :: material_out

    ! Local variables
    integer :: iTmp

    material_out%name = material_in%name
!    material_out%material_type = material_in%material_type
    material_out%defined = material_in%defined
    material_out%mole_mass = material_in%mole_mass
    material_out%ifAerosolPossible = material_in%ifAerosolPossible
    material_out%aerosol_param = material_in%aerosol_param
    material_out%gas_dep_param = material_in%gas_dep_param  ! inside are just pointers
    if(associated(material_in%reaction_data))then
      material_out%reaction_data => material_in%reaction_data
    else
      nullify(material_out%reaction_data)
    endif
    if(associated(material_in%pollen))then
      material_out%pollen => material_in%pollen
    else
      nullify(material_out%pollen)
    endif
    if(associated(material_in%nuc_data))then
      material_out%nuc_data => material_in%nuc_data  ! inside are just pointers
    else
      nullify(material_out%nuc_data )
    endif
  end subroutine copy_material


  !=======================================================================

  logical function fu_compare_materials_eq(material1, material2)
    !
    ! Compares two materials
    !
    implicit none
    type(silam_material), intent(in) :: material1, material2

    if(material1%defined == silja_true)then
      if(material2%defined == silja_true)then
        fu_compare_materials_eq = trim(material1%name) == trim(material2%name) .and. &
!                                & material1%material_type == material2%material_type .and. &
                                & material1%mole_mass == material2%mole_mass
      else
        fu_compare_materials_eq = .false.
      endif
    else
      fu_compare_materials_eq = material1%defined == material2%defined
    endif
  end function fu_compare_materials_eq

  !=======================================================================
  logical FUNCTION fu_material_defined (material)
    IMPLICIT NONE
    TYPE (silam_material), intent(in) :: material
    fu_material_defined = material%defined == silja_true
  END FUNCTION fu_material_defined
  
  !=======================================================================
  function fu_name_of_material(material) result(name)
    implicit none
    character(len = substNmLen) :: name
    type(silam_material), intent(in) :: material
    if(material%defined == silja_true)then
      name = trim(material%name)
    else
      name = "missing_material"
    end if
  end function fu_name_of_material

  !=======================================================================
  function fu_nuclide_from_material(material) result(nuclide)
    !
    ! If the material is radioactive, returns the nuclide pointer
    !
    implicit none
    type(silam_nuclide), pointer :: nuclide
    type(silam_material), intent(in) :: material

    nullify(nuclide)
    if(.not. (material%defined == silja_true))then
      call set_error('Undefined material','fu_nuclide_from_material')
      return
    endif
    if(.not.associated(material%nuc_data))then
      call set_error('Non-radioactive material:' + fu_name(material), 'fu_nuclide_from_material')
      return
    endif
    nuclide => material%nuc_data
  end function fu_nuclide_from_material


  !*********************************************************************
  !
  ! ENCAPSULATION of chemistry and deposition
  !
  !*********************************************************************

  logical function fu_if_gas_depositing(material)
    implicit none
    type(silam_material), intent(in) :: material
    fu_if_gas_depositing = fu_true(material%gas_dep_param%ifGasDeposited)
  end function fu_if_gas_depositing

  !====================================================================
  real function fu_gas_molecular_diffusivity_air(material)
    implicit none
    type(silam_material), intent(in) :: material
    fu_gas_molecular_diffusivity_air = material%gas_dep_param%gas_molecular_diffus_in_air
  end function fu_gas_molecular_diffusivity_air

  !====================================================================
  real function fu_molec_diffus_water_gas(material)
    implicit none
    type(silam_material), intent(in) :: material
    fu_molec_diffus_water_gas = material%gas_dep_param%gas_molecular_diffus_in_water
  end function fu_molec_diffus_water_gas

  !====================================================================
  real function fu_Wesely_f0(material)
    implicit none
    type(silam_material), intent(in) :: material
    fu_Wesely_f0 = material%gas_dep_param%Wesely_f0
  end function fu_Wesely_f0

  !====================================================================
  real function fu_Henry_const_298K(material)
    implicit none
    type(silam_material), intent(in) :: material
    fu_Henry_const_298K = material%gas_dep_param%Henry_const_298K
  end function fu_Henry_const_298K

  !====================================================================
  real function fu_Rs_default_soil_gas(material)
    implicit none
    type(silam_material), intent(in) :: material
    fu_Rs_default_soil_gas = material%gas_dep_param%Rs_gas_land_default
  end function fu_Rs_default_soil_gas

  !====================================================================
  real function fu_Rs_default_water_gas(material)
    implicit none
    type(silam_material), intent(in) :: material
    fu_Rs_default_water_gas = material%gas_dep_param%Rs_gas_water
  end function fu_Rs_default_water_gas

  !====================================================================
  function fu_gas_deposition_param(material) result(DepData)     !!!!!!!! pointer
    implicit none
    type(Tgas_deposition_param), pointer :: DepData
    type(silam_material), target, intent(in) :: material
    DepData => material%gas_dep_param
  end function fu_gas_deposition_param

  !====================================================================
  function fu_aerosol_param(material) result(AerData)            !!!!!!!!! pointer
    implicit none
    type(silam_material), target, intent(in) :: material
    type(Taerosol_param), pointer:: AerData
    AerData => material%aerosol_param
  end function fu_aerosol_param

  !====================================================================
  real function fu_scavenging_scale_gas_rain(material)
    implicit none
    type(silam_material), intent(in) :: material
    fu_scavenging_scale_gas_rain = material%gas_dep_param%washout_scale_rain_to_standard
  end function fu_scavenging_scale_gas_rain

  !====================================================================
  real function fu_scavenging_scale_gas_snow(material)
    implicit none
    type(silam_material), intent(in) :: material
    fu_scavenging_scale_gas_snow = material%gas_dep_param%washout_scale_snow_to_standard
  end function fu_scavenging_scale_gas_snow


  !*******************************************************************
  !*******************************************************************
  !
  !
  !    UNIT CONVERSION ROUTINES FOR MATERIALS
  !
  !
  !*******************************************************************
  !*******************************************************************

  !*******************************************************************

  real function fu_mole_to_kg_by_material(material, moles) result(amount_kg)
    !
    ! Transfers the mount in moles of some material to kg.
    !
    implicit none

    type(silam_material), intent(in) :: material
    real, intent(in) :: moles

    if(material%mole_mass < 1.0e-10)then
      call msg('Material:' + material%name, material%mole_mass)
      call set_error('Material without mole mass','fu_mole_to_kg_by_material')
      amount_kg = real_missing
    else
      amount_kg = moles * material%mole_mass
    endif

  end function fu_mole_to_kg_by_material


  ! ******************************************************************

  real function fu_conversion_factor_material(chUnitFrom, chUnitTo, material)
    !
    ! This function makes the conversion from one unit to another. The idea
    ! is to enable the user to enter e.g. emission flux in some unit, e.g., tons
    ! per day instead of the model standard moles per second.
    ! Since moles to kg and Bq can be translated only knowing the specific 
    ! substance, it is also required as an input parameter
    !
    ! Idea: conversion between the units is performed in two steps. First, 
    ! unitFrom is converted to some basic SI unit, second, the UnitTo is 
    ! obtained from that basic SI unit.
    ! Material may be not needed for some conversions, so let it be optional
    !
    implicit none

    ! Imported parameters with intent IN
    !
    character(len=*), intent(in) :: chUnitFrom, chUnitTo
    type(silam_material), intent(in) :: material

    ! Local variables
    real :: tmpFactor
    character(len=clen) :: chSubUnit1, chSubUnit2, chSubUnit3, chSubUnit4, chBasicUnit
    integer :: unitType, iSlashFrom, iSlashTo

    if(chUnitFrom == "" .or. chUnitTo == "")then
      call set_error('Empty unit given','fu_conversion_factor_material')
      return
    endif
    tmpFactor = 1.
    iSlashFrom = index(chUnitFrom,"/")
    iSlashTo = index(chUnitTo,"/")
    if(iSlashFrom == 0)then ! [unit1] -> [unit3]
      if(iSlashTo == 0)then
        chSubUnit1 = chUnitFrom
        chSubUnit3 = chUnitTo
        chSubUnit2 =''  ! No denominators
        chSubUnit4 =''
      else
        call set_error(fu_connect_strings('Incompatible units:',chUnitFrom,'<->',chUnitTo), &
                     & 'fu_conversion_factor_material')
        return
      endif
    else   ! [unit_1]/[unit_2] -> [unit_3]/[unit_4]
      if(iSlashTo == 0)then
        call set_error(fu_connect_strings('Incompatible units:',chUnitFrom,'<->',chUnitTo), &
                     & 'fu_conversion_factor_material')
        return
      else
        chSubUnit1 = chUnitFrom(1:iSlashFrom-1)
        chSubUnit3 = chUnitTo(1:iSlashTo-1)
        chSubUnit2 = chUnitFrom(iSlashFrom+1 : )
        chSubUnit4 = chUnitTo(iSlashTo+1 : )
      endif
    endif

    !------------------------------------------------------
    ! 
    ! Unit nominators: select basic unit and transform each given unit to basic one
    ! Reason for that is: sometimes properly selected basic unit will make the task 
    ! much simpler and, e.g., independent from the availability of material
    !
!    call msg('Select basic unit')
    call select_basic_unit(chSubUnit1, chSubUnit3, chBasicUnit)
!    call msg(fu_connect_strings('Selected:',chBasicUnit))

    tmpFactor = tmpFactor * fu_factor_to_basic_unit(chSubUnit1, chBasicUnit, material)

!    call msg('Calling unit type...')

    unitType = fu_unit_type(chSubUnit1)

!    call msg('Unit1 to basic - type and factor:',unitType,tmpFactor)

    if(unitType /= fu_unit_type(chSubUnit3))then
      call set_error(fu_connect_strings('Incompatible units:',chUnitFrom,'_and_',chUnitTo), &
                   & 'fu_conversion_factor_material')
      return
    endif

    tmpFactor = tmpFactor / fu_factor_to_basic_unit(chSubUnit3, chBasicUnit, material)

!    call msg('Unit 3 to basic: new factor',real_value=tmpFactor)

    !------------------------------------------------------
    !
    ! Unit denominators:
    !
    if(iSlashFrom /= 0)then

      call select_basic_unit(chSubUnit2, chSubUnit4, chBasicUnit)

      tmpFactor = tmpFactor / fu_factor_to_basic_unit(chSubUnit2, chBasicUnit, material)

      unitType = fu_unit_type(chSubUnit2)

      if(unitType /= fu_unit_type(chSubUnit4))then
        call set_error('Incompatible units:' + chUnitFrom + '_and_' + chUnitTo, &
                     & 'fu_conversion_factor_material')
        return
      endif

      tmpFactor = tmpFactor * fu_factor_to_basic_unit(chSubUnit4, chBasicUnit, material)

    endif  ! Denominators exist

    if(error)then
      fu_conversion_factor_material = -1.
    else
      fu_conversion_factor_material = tmpFactor
    endif

  end function fu_conversion_factor_material


  !****************************************************************************

  recursive real function fu_factor_to_basic_unit_mater(chUnit, chBasicUnit, material) &
                                         & result(factor_to_basic_unit)
    !
    ! Finds the conversion factor from some unit to the appropriate basic unit
    !
    implicit none
    character(len=*), intent(in) :: chUnit, chBasicUnit
    type(silam_material), intent(in) :: material

    ! Local variables
    integer :: iPower, status, lTmp
    real :: fTmp
    !
    ! Check first if the unit is in some power - 2, 3 or 4 are allowed. 
    !
    lTmp = len_trim(chUnit)
!    call msg(fu_connect_strings('Unit before power and len_trim:',chUnit), lTmp)

    if(fu_if_digit(chUnit(lTmp:lTmp)))then
      read(unit=chUnit(lTmp:lTmp),fmt=*,iostat=status) iPower
!      call msg(fu_connect_strings('Unit and len_trim after power:',chUnit),lTmp, real(iPower))
      if(status == 0 .and. iPower > 1 .and. iPower < 5)then
!        call msg('recursing...')
        factor_to_basic_unit = fu_factor_to_basic_unit(chUnit(1:lTmp-1), &
                                                     & chBasicUnit, &
                                                     & material)
        !
        ! Take care of the power
        !
        fTmp = factor_to_basic_unit
        do status = 1, iPower-1
          factor_to_basic_unit = factor_to_basic_unit * fTmp
        end do
        return
      endif
    endif
    !
    ! Now we should have a basic unit, which is defined by the substance type.
    !
!    call msg ('case')
    select case(trim(chUnit))
        !
        ! Amount units -> basic units. Here mole is default, the others are case-specific
        !
      case('mkg')
        if(chBasicUnit == 'kg')then
          factor_to_basic_unit = 0.000000001
        elseif(chBasicUnit == 'mole')then
          factor_to_basic_unit = 0.000000001/fu_mole_to_kg_by_material(material,1.)
        elseif(chBasicUnit == 'Bq')then
          factor_to_basic_unit = 0.000000001 * fu_mole_to_bq(fu_nuclide(material),1.) / &
                                             & fu_mole_to_kg_by_material(material,1.)
        else
          call set_error(fu_connect_strings('Failed factor to basic unit:',chUnit,',',chBasicUnit), &
                       & 'fu_factor_to_basic_unit_mater')
          factor_to_basic_unit = 1. ! Zero is forbidden!!
        endif

      case('mg')
        if(chBasicUnit == 'kg')then
          factor_to_basic_unit = 0.000001
        elseif(chBasicUnit == 'mole')then
          factor_to_basic_unit = 0.000001/fu_mole_to_kg_by_material(material,1.)
        elseif(chBasicUnit == 'Bq')then
          factor_to_basic_unit = 0.000001 * fu_mole_to_bq(fu_nuclide(material),1.) / &
                                          & fu_mole_to_kg_by_material(material,1.)
        else
          call set_error(fu_connect_strings('Failed factor to basic unit:',chUnit,',',chBasicUnit), &
                       & 'fu_factor_to_basic_unit_mater')
          factor_to_basic_unit = 1. ! Zero is forbidden!!
        endif

      case('g')
        if(chBasicUnit == 'kg')then
          factor_to_basic_unit = 0.001
        elseif(chBasicUnit == 'mole')then
          factor_to_basic_unit = 0.001/fu_mole_to_kg_by_material(material,1.)
        elseif(chBasicUnit == 'Bq')then
          factor_to_basic_unit = 0.001 * fu_mole_to_bq(fu_nuclide(material),1.) / &
                                       & fu_mole_to_kg_by_material(material,1.)
        else
          call set_error(fu_connect_strings('Failed factor to basic unit:',chUnit,',',chBasicUnit), &
                       & 'fu_factor_to_basic_unit_mater')
          factor_to_basic_unit = 1. ! Zero is forbidden!!
        endif

      case('kg')
        if(chBasicUnit == 'kg')then
          factor_to_basic_unit = 1.
        elseif(chBasicUnit == 'mole')then
          factor_to_basic_unit = 1./fu_mole_to_kg_by_material(material,1.)
        elseif(chBasicUnit == 'Bq')then
          factor_to_basic_unit = fu_mole_to_bq(fu_nuclide(material),1.) / &
                               & fu_mole_to_kg_by_material(material,1.)
        else
          call set_error(fu_connect_strings('Failed factor to basic unit:',chUnit,',',chBasicUnit), &
                       & 'fu_factor_to_basic_unit_mater')
          factor_to_basic_unit = 1. ! Zero is forbidden!!
        endif

      case('Mg','ton')
        if(chBasicUnit == 'kg')then
          factor_to_basic_unit = 1000.
        elseif(chBasicUnit == 'mole')then
          factor_to_basic_unit = 1000./fu_mole_to_kg_by_material(material,1.)
        elseif(chBasicUnit == 'Bq')then
          factor_to_basic_unit = 1000. * fu_mole_to_bq(fu_nuclide(material),1.) / &
                                       & fu_mole_to_kg_by_material(material,1.)
        else
          call set_error(fu_connect_strings('Failed factor to basic unit:',chUnit,',',chBasicUnit), &
                       & 'fu_factor_to_basic_unit_mater')
          factor_to_basic_unit = 1. ! Zero is forbidden!!
        endif

      case('Gg','kton')
        if(chBasicUnit == 'kg')then
          factor_to_basic_unit = 1000000.
        elseif(chBasicUnit == 'mole')then
          factor_to_basic_unit = 1000000./fu_mole_to_kg_by_material(material,1.)
        elseif(chBasicUnit == 'Bq')then
          factor_to_basic_unit = 1000000. * fu_mole_to_bq(fu_nuclide(material),1.) / &
                                          & fu_mole_to_kg_by_material(material,1.)
        else
          call set_error('Failed factor to basic unit:' + chUnit + ',' + chBasicUnit, &
                       & 'fu_factor_to_basic_unit_mater')
          factor_to_basic_unit = 1. ! Zero is forbidden!!
        endif

      case('Mton')
        if(chBasicUnit == 'kg')then
          factor_to_basic_unit = 1000000000.
        elseif(chBasicUnit == 'mole')then
          factor_to_basic_unit = 1000000000./fu_mole_to_kg_by_material(material,1.)
        elseif(chBasicUnit == 'Bq')then
          factor_to_basic_unit = 1000000000. * fu_mole_to_bq(fu_nuclide(material),1.) / &
                                             & fu_mole_to_kg_by_material(material,1.)
        else
          call set_error('Failed factor to basic unit:' + chUnit + ',' + chBasicUnit, &
                       & 'fu_factor_to_basic_unit_mater')
          factor_to_basic_unit = 1. ! Zero is forbidden!!
        endif

      case('mole')
        if(chBasicUnit == 'mole')then
          factor_to_basic_unit = 1.
        elseif(chBasicUnit == 'Bq')then
          factor_to_basic_unit = fu_mole_to_bq(fu_nuclide(material),1.)
        else
          call set_error('Failed factor to basic unit:' + chUnit + ',' + chBasicUnit, &
                       & 'fu_factor_to_basic_unit_mater')
          factor_to_basic_unit = 1. ! Zero is forbidden!!
        endif

      case('kmole')
        if(chBasicUnit == 'mole')then
          factor_to_basic_unit = 1000.
        elseif(chBasicUnit == 'Bq')then
          factor_to_basic_unit = 1000. * fu_mole_to_bq(fu_nuclide(material),1.)
        else
          call set_error('Failed factor to basic unit:' + chUnit + ',' + chBasicUnit, &
                       & 'fu_factor_to_basic_unit_mater')
          factor_to_basic_unit = 1. ! Zero is forbidden!!
        endif

      case('Mmole')
        if(chBasicUnit == 'mole')then
           factor_to_basic_unit = 1000000.
        elseif(chBasicUnit == 'Bq')then
          factor_to_basic_unit = 1000000. * fu_mole_to_bq(fu_nuclide(material),1.)
        else
          call set_error('Failed factor to basic unit:' + chUnit + ',' + chBasicUnit, &
                       & 'fu_factor_to_basic_unit_mater')
          factor_to_basic_unit = 1. ! Zero is forbidden!!
        endif

      case('Bq')
        if(chBasicUnit == 'Bq')then
          factor_to_basic_unit = 1.
        elseif(chBasicUnit == 'mole')then
          factor_to_basic_unit = 1./fu_mole_to_bq(fu_nuclide(material),1.)
        else
          call set_error('Failed factor to basic unit:' + chUnit + ',' + chBasicUnit, &
                       & 'fu_factor_to_basic_unit_mater')
          factor_to_basic_unit = 1. ! Zero is forbidden!!
        endif

      case('kBq')
        if(chBasicUnit == 'Bq')then
          factor_to_basic_unit = 1000.
        elseif(chBasicUnit == 'mole')then
          factor_to_basic_unit = 1000./fu_mole_to_bq(fu_nuclide(material),1.)
        else
          call set_error('Failed factor to basic unit:' + chUnit + ',' + chBasicUnit, &
                       & 'fu_factor_to_basic_unit_mater')
          factor_to_basic_unit = 1. ! Zero is forbidden!!
        endif

      case('MBq')
        if(chBasicUnit == 'Bq')then
          factor_to_basic_unit = 1000000.
        elseif(chBasicUnit == 'mole')then
          factor_to_basic_unit = 1000000./fu_mole_to_bq(fu_nuclide(material),1.)
        else
          call set_error('Failed factor to basic unit:' + chUnit + ',' + chBasicUnit, &
                       & 'fu_factor_to_basic_unit_mater')
          factor_to_basic_unit = 1. ! Zero is forbidden!!
        endif

      case default
        !
        ! All units not related to the mass transformation are handled in toolbox
        !
        factor_to_basic_unit = fu_factor_to_basic_unit(chUnit, chBasicUnit)

    end select

!    call msg('Done case, factor:',real_value=factor_to_basic_unit)

  end function fu_factor_to_basic_unit_mater


  !******************************************************************************

  function fu_low_mass_threshold_material(material) result(thrsh)
    ! 
    ! Return the low mass threshold for the material. This is the last
    ! resort for defining one for the species.
    !
    implicit none
    type(silam_material), intent(in) :: material
    real :: thrsh
    thrsh = material%low_mass_thrsh
  end function fu_low_mass_threshold_material

  !************************************************************************************

  subroutine test_amount_units()
    implicit none
    character(len=2), parameter, dimension(5) :: mod_names = (/'mk', ' m', ' k', ' M', ' G'/)
    real, parameter, dimension(5) :: mod_values = (/1e-6, 1e-3, 1e3, 1e6, 1e9/)
    type(silam_material), pointer :: no2, i_131, pm
    real :: mole_to_bq
    character(len=*), parameter :: sub_name = 'test_amount_units'

    no2 => fu_get_material_ptr('NO2')
    if (error) return

    call msg('Testing unit conversions...')

    if (fu_fails(test_conversion('mole', 'g', no2, 46.0, .true.), 'Failed conversion', sub_name)) return
    if (fu_fails(test_conversion('mole', 'mole', no2, 1.0, .true.), 'Failed conversion', sub_name)) return
    if (fu_fails(test_conversion('g', 'g', no2, 1.0, .true.), 'Failed conversion', sub_name)) return
    
    if (fu_fails(test_conversion('ton', 'g', no2, 1e6, .false.), 'Failed conversion', sub_name)) return
    if (fu_fails(test_conversion('ton', 'kg', no2, 1e3, .false.), 'Failed conversion', sub_name)) return
    if (fu_fails(test_conversion('Mton', 'g', no2, 1e12, .false.), 'Failed conversion', sub_name)) return
    if (fu_fails(test_conversion('kton', 'g', no2, 1e9, .false.), 'Failed conversion', sub_name)) return

    pm => fu_get_material_ptr('PM')
    if (error) return
    call msg('Will test conversion for PM')
    if (fu_fails(test_conversion('g', 'g', pm, 1.0, .true.), 'Failed conversion', sub_name)) return
    
    ! Don't have molar mass or nuclide data, but mole->mole or Bq->Bq should still work!
    if (fu_fails(test_conversion('mole', 'mole', pm, 1.0, .true.), 'Failed conversion', sub_name)) return
    if (fu_fails(test_conversion('Bq', 'kBq', pm, 1e-3, .false.), 'Failed conversion', sub_name)) return
    
    call msg('** Expected failures:')
    if (fu_fails(.not.test_conversion('mole', 'g', pm, 46.0, .true.), 'Error not set', sub_name)) return
    if (fu_fails(.not.test_conversion('number', 'g', pm, 46.0, .true.), 'Error not set', sub_name)) return
    if (fu_fails(.not.test_conversion('hour', 'g', pm, 46.0, .true.), 'Error not set', sub_name)) return

    call unset_error(sub_name)

    i_131 => fu_get_material_ptr('I_131')
    if (error) return
    call msg('Will test conversion for I-131')
    if (fu_fails(test_conversion('mole', 'g', i_131, 131.0, .true.), 'Failed conversion', sub_name)) return
    mole_to_bq = fu_mole_to_bq(fu_nuclide(i_131), 1.0)
    ! some modifiers for Bq are not supported, test them separately:
    if (fu_fails(test_conversion('mole', 'Bq', i_131, mole_to_bq, .false.), 'Failed conversion', sub_name)) return
    if (fu_fails(test_conversion('mole', 'kBq', i_131, mole_to_bq*1e-3, .false.), &
               & 'Failed conversion', sub_name)) then
      return
    end if
    if (fu_fails(test_conversion('mmole', 'kBq', i_131, mole_to_bq*1e-6, .false.), &
               & 'Failed conversion', sub_name)) then
      return
    end if

  contains

    logical function test_conversion(unit1, unit2, material, expected, test_modifiers) result(ok)
      character(len=*), intent(in) :: unit1, unit2
      real, intent(in) :: expected
      type(silam_material), intent(in) :: material
      logical, intent(in) :: test_modifiers

      character(len=clen) :: unit1_mod, unit2_mod, mod1_name, mod2_name
      integer :: ind_mod1, ind_mod2, num_modifiers
      real :: mod_total_val, factor

      ok = .false.

      factor = fu_conversion_factor(unit1, unit2, material)
      if (.not. (factor .eps. expected)) then
        call msg('Expected:', expected)
        call msg('Obtained:', factor)
        call msg('Failed conversion: ' + unit1 + ' -> ' + unit2)
        return
      else
        call msg('OK: ' + unit1 + ' -> ' + unit2)
      end if
        
      factor  = fu_conversion_factor(unit2, unit1, material)
      if (.not. (factor .eps. 1/expected)) then
        call msg('Expected:', 1 / expected)
        call msg('Obtained:', factor)
        call msg('Failed conversion: ' + unit2 + ' -> ' + unit1)
        return
      else
        call msg('OK: ' + unit2 + ' -> ' + unit1)
      end if

      if (.not. test_modifiers) then
        ok = .true.
        return
      end if

      num_modifiers = size(mod_names)
      do ind_mod1 = 1, num_modifiers
        do ind_mod2 = 1, num_modifiers
          mod1_name = mod_names(ind_mod1)
          mod2_name = mod_names(ind_mod2)
          
          unit1_mod = trim(adjustl(mod1_name)) // trim(unit1)
          unit2_mod = trim(adjustl(mod2_name)) // trim(unit2)
          
          mod_total_val = mod_values(ind_mod1) / mod_values(ind_mod2)
      
          factor = fu_conversion_factor(unit1_mod, unit2_mod, material)
          if (.not. (factor .eps. expected * mod_total_val)) then
            call msg('Expected:', expected * mod_total_val)
            call msg('Obtained:', factor)
            call msg('Failed conversion: ' + unit1_mod + ' -> ' + unit2_mod)
            return
          else
            call msg('OK: ' + unit1_mod + ' -> ' + unit2_mod)
          end if
          
          factor  = fu_conversion_factor(unit2_mod, unit1_mod, material)
          if (.not. (factor .eps. 1 / (expected*mod_total_val))) then
            call msg('Expected:', 1 / (expected * mod_total_val))
            call msg('Obtained:', factor)
            call msg('Failed conversion: ' + unit2_mod + ' -> ' + unit1_mod)
            return
          else
            call msg('OK: ' + unit2_mod + ' -> ' + unit1_mod)
          end if
          
        end do
      end do

      ok = .true.
      
    end function test_conversion

  end subroutine test_amount_units

END MODULE materials
