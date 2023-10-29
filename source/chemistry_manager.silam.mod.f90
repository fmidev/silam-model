module chemistry_manager
  !
  ! The module controls the set of chemical and physical transformations
  ! in SILAM. Its task is to initialise the chemico-physical modules, collect
  ! their input needs and then dispatch the transformations during the simulations
  ! At present, SILAM has several independent transformation lines:
  ! - chemical and radioactive
  ! - aerosol dynamics
  ! - deposition
  ! These lines are parallel to each other, i.e., each of them can be defined independently
  ! from the other lines. Some compatibility limitations exist but made as small as possible.
  !
  ! Code owners: Mikhail Sofiev, Julius Vira, Marje Prank, email <firstname>.<lastname>@fmi.fi
  !
  ! Units: SI
  !
  ! Code: ANSI FORTRAN-90 (or close to it)
  !

  use source_terms_general
  use optical_density
  use depositions
  use chem_dep_passive
  use chem_dep_pollen
  use chem_dep_pm_general
  use chem_dep_sulphur_dmat
  use chem_dep_acid_basic
  use chem_dep_cbm4
  use chem_dep_cbm4_SOA
  use chem_dep_cbm42_strato
  use chem_dep_cbm42_strato_SOA
  use chem_dep_cbm5_SOA
  use chem_dep_cbm5_strato_SOA
  use chem_dep_radioactive
  use aer_dyn_basic
  use aer_dyn_simple
  use photolysis
  use aer_dyn_middle_atmosph
  use aer_dyn_VBS
  use ini_boundary_conditions
  
  !$use omp_lib

  implicit none

  private


!!#define DUMP
!!#define FORCING_TOP

#ifdef DEBUG_TRANSFORMATIONS
  logical, parameter :: if_OMP_chemistry = .false.
#else
  logical, parameter :: if_OMP_chemistry = .true.
#endif

  
  !
  ! Public functions 
  !
  public global_chemical_init
  public add_transformation_input_needs   ! for general transofrmations
  public transform_maps
  public transform_lagrangian_part
  public transform_deposited
  public make_wet_deposition_lagr

  public set_chemistry_rules     ! for Global chemical rules
  public set_missing
  public set_low_mass_threshold
  public fu_mass_lagrangian_particle
  public fu_number_of_lagr_particles
  public init_chemicalRunSetup
  public fu_low_mass_threshold
  public fu_aerosol_rules
  public fu_deposition_rules
  public fu_chemRunSetup
  public fu_optical_rules
  public point_low_mass_thresholds
  public fu_low_mass_thres_defined
  public get_tla_chemistry
  ! Private functions
  !
  private init_chemisryStuff
  private set_chemrules_missing
  private fu_low_mass_thresholds
  private fu_compute_thresh_for_release
  private fu_lagr_part_mass_for_release

  ! Interface for the generic functions

  interface set_missing
     module procedure set_chemrules_missing
  end interface

  interface fu_low_mass_threshold
    module procedure fu_low_mass_thresholds
  end interface
  
  !-----------------------------------------------------------------------
  !
  ! Maximum number of transformations used in one run
  !
  integer, private, parameter :: max_transformations = 32
  integer, parameter, public :: max_aux_cocktails = 10

  ! Couple more parameters
  integer, private, parameter :: LowMassThreshUseEmission = 11000 ! From emission
  integer, private, parameter :: LowMassThreshUseDefault  = 11001 ! from Chemicals
  integer, private, parameter :: LowMassThreshUseConst = 11002 ! one_size_fits_all
  real, public, parameter :: fLowMassThreshConstant = 1e-20

  type Tchem_rules
    !private
    integer, dimension(max_transformations) :: iTransformTypes, iAerosolDynTypes, &
                                             & iWhomToApplyTransform, iWhomToApplyAerDyn
    integer :: nTransformations=int_missing, nAerosolDynamics=int_missing
    character(len=fnlen) :: filename_photo_lut = char_missing
    logical :: need_photo_lut = .false.
    logical :: useDynamicAlbedo = .false.        ! By default use static albedo for photolysis
    real :: defaultStaticAlbedo = 0.3            ! Some default value   
    integer :: PhotoO3col = int_missing          ! Use O3 column to calculate photolysis 
                                          !can be standard_atmosphere, massmap or meteo_column 
    logical :: ifPhotoAOD  = .false.             ! Use AOD to calculate photolysis
    integer :: cloud_model_for_photolysis = simple_cloud
    real :: photoAODwavelength = real_missing    ! Wavelength to use for photolysis attenuation by AOD
    logical :: ifOnesAdjust  = .false.             ! Use ONES tracer to adjust concentrations
    type(Tchem_rules_passive) :: rulesPassive    ! Passive self-degrading substance
    type(Tchem_rules_pm_general) :: rulesPM      ! inert particulate matter
    type(Tchem_rules_radioactive) :: rulesRadioactive
    type(Tchem_rules_pollen) :: rulesPollen
    type(Tchem_rules_DMAT_S) :: rulesSulphurDMAT
    type(Tchem_rules_acidBasic) :: rulesAcidBasic
!    !type(Tchem_rules_sea_salt) :: rulesSeaSalt
!    !type(Tchem_rules_POP) :: rulesPOP
    type(Tchem_rules_cbm4) :: rulesCBM4
    type(Tchem_rules_cbm4_SOA) :: rulesCBM4_SOA
    type(Tchem_rules_cbm42_strato) :: rulesCBM42_strato
    type(Tchem_rules_cbm42_strato_SOA) :: rulesCBM42_strato_SOA
    type(Tchem_rules_cbm5_SOA) :: rulesCBM5_SOA
    type(Tchem_rules_cbm5_strato_SOA) :: rulesCBM5_strato_SOA
    type(Tchem_rules_AerDynBasic) :: rulesAerDynBasic
    type(Tchem_rules_AerDynSimple) :: rulesAerDynSimple
    type(Tchem_rules_AerDynMidAtm) :: rulesAerDynMidAtmosph
    type(Tchem_rules_AerDynVBS) :: rulesAerDynVBS
    type(silam_species), dimension(:), allocatable :: expected_species   ! if aerosol dynamics wants it strict
    type(Taerosol_rules) :: rulesAerosol
    type(Toptical_density_rules), pointer :: rulesOpticDens  => null()
    type(TchemicalRunSetup), pointer :: ChemRunSetup  => null()
    type(Tdeposition_rules), pointer :: rulesDeposition  => null()
    ! The below must be nullified, if not here then in set_chem_rules.
    real, dimension(:), pointer :: low_cnc_trsh => null(), low_mass_trsh => null(), &
                                 & mass_lagr_particle => null(), &
                                 & low_cnc_trsh_frw => null(), low_mass_trsh_frw => null(), &
                                 & low_cnc_trsh_adj => null(), low_mass_trsh_adj => null()
    integer :: LowMassThresh  = LowMassThreshUseEmission
    character(len=clen), dimension(max_aux_cocktails) :: auxCocktName = ''
    integer :: nAuxCocktails=int_missing
    type(silja_logical) :: defined = silja_false

  end type Tchem_rules
  public Tchem_rules



  ! Allowed types of transformations. INFO ONLY HERE. The real declarations are
  ! in the corresponding modules
  ! 
  ! need to make these public for the box model:
  ! 
  public transformation_passive
  !public species_radioactive
  !public species_pollen
  public transformation_inert_PM
  public transformation_sulphur_dmat
  public transformation_acid_basic
  public transformation_cbm4
  public transformation_cbm4_SOA
  public transformation_radioactive
  public transformation_cbm42_strato
  public transformation_cbm42_strato_SOA
  public transformation_cbm5_SOA
  public transformation_cbm5_strato_SOA

  public aerosol_dynamics_basic
  public aerosol_dynamics_simple
  public aerosol_dynamics_Mid_Atmosph
  public aerosol_dynamics_VBS
  
  public transform_passive
  public transform_cbm4
  public transform_cbm4_SOA
  public transform_cbm4_adj
  public transform_aerdynbasic
  public transform_acid_basic
  public transform_dmat
  public transform_radioactive
  public transform_aerdynsimple
  public transform_cbm42_strato
  public transform_cbm42_strato_adj
  public transform_cbm42_strato_SOA
  public transform_cbm42_strato_SOA_adj
  public transform_cbm5_SOA
  public transform_cbm5_SOA_adj
  public transform_cbm5_strato_SOA
  public transform_cbm5_strato_SOA_adj
  public transform_AerDynMidAtm
  public transform_AerDynVBS
  public prepare_step_cbm4
  public prepare_step_cbm4_SOA
  public prepare_step_cbm42_strato
  public prepare_step_cbm42_strato_SOA
  public prepare_step_cbm5_SOA
  public prepare_step_cbm5_strato_SOA
  !--------------------------------------------------------------------
  !
  ! The global meteo_input data structure and an array of the local meteo_input ones
  !
  type(Tmeteo_input), dimension(:), allocatable, private, target, save :: pMeteo_input_local

  real, dimension(:,:,:), allocatable, save, private :: cb4_h_start_array
  !
  ! Thread-owned scratch stuff
  !
  type Tchemical_thread_stuff
    real, dimension(:,:),  allocatable :: garb_array  ! nSrc, nSpecies
    real, dimension(:,:),  allocatable :: cncTrn ! nSpecies,1:mapTransport%nSrc
    real, dimension(:),  allocatable :: cncAer ! 1:mapAerosol%nSpecies
    real, dimension(:),  allocatable :: cncSL  ! 1:mapShortLived%nSpecies
    real, dimension(:,:),  allocatable :: photorates ! 1:num_react, 1:nz_dispersion
    real, dimension(:),  allocatable :: aodext ! 1:nz_disperson
    real, dimension(:),  allocatable :: aodscat ! 1:nz_disperson
    real, dimension(:,:),  allocatable :: o3column ! 1:nz_disperson, 1:nSrc     !Ozone column above
    real, dimension(:),  allocatable :: soot_col ! 1:nz_disperson
    real, dimension(:),  allocatable :: pwc_col ! 1:nz_disperson
    real, dimension(:),  allocatable :: tau_above_bott ! 1:nz_disperson
    real, dimension(:,:),  allocatable :: metdat ! meteo_input%nQuantities
    real, dimension(:,:),  allocatable :: reactRates ! 1:nSrc, 1:nRates
  end type Tchemical_thread_stuff

  type Tchemical_stuff
     type (Tchemical_thread_stuff), dimension(:), allocatable :: arrStuff
     type (Silja_Logical) :: defined = silja_false
  end type Tchemical_stuff

  type (Tchemical_stuff), private,  target,save :: ChemStuff  ! Container for temporary things

#ifdef FORCING_TOP
integer, private, save :: indOzoneToForce = int_missing
#endif




CONTAINS

  !***********************************************************************************


  subroutine global_chemical_init(chemRules, timestep, timestep_output, &  ! input
                                & speciesEmission, nSpeciesEmission, &     ! input
                                & speciesTransport, nSpeciesTransport, &      ! output
                                & speciesShortlived, nSpeciesShortlived, &    ! output
                                & speciesAerosol, nSpeciesAerosol, nReactRates)  ! output
    ! 
    ! This subroutine (eventually) collects the emission species from
    ! the source, the required transport species from transformations,
    ! and sets the chemistry rules for individual modules. It should
    ! be called when the source is chemically initialized.
    ! 
    ! The algorithm is as follows:
    ! 
    ! 1. Get a list of emission species from the source. These for the
    ! emission cocktail. 
    !
    ! 2. The transformations requested in chemical
    ! rules can now add the species they require to the *transport*
    ! cocktail. Transformations (currently radioactive and cb4) that
    ! require specific species indexing are given priority. Actually,
    ! the emission species could be handed to the transformations at
    ! this point (for instance, radioactive will then choose the
    ! nuclide chains). 
    !
    ! 3. Finally, if some emission species are not in
    ! the transport cocktail, they are added.
    !
    implicit none
    
    !type(Tsilam_namelist), intent(in) :: nlChem ! Chemistry namelist
!    type(silam_source), intent(in) :: source

    type(Tchem_rules), intent(inout) :: chemRules
    type(silja_interval), intent(in) :: timestep, timestep_output
    type(silam_species), dimension(:),  pointer :: speciesEmission                     ! input
    type(silam_species), dimension(:),  pointer :: speciesTransport, speciesShortlived, &
                                                 & speciesAerosol ! output
    integer, intent(in) :: nspeciesEmission
    integer, intent(out) :: nspeciesTransport, nspeciesShortlived, nSpeciesAerosol, nReactRates
    integer, dimension(max_species) :: indDepositionType
    integer, dimension(max_species) :: iClaimedSpecies

    ! Local declarations
    integer :: iTmp, stat, nspecies_aux
    type(silam_species), dimension(:), pointer :: speciesTmp, speciesTmp2, speciesTmp3, species_aux
    logical :: ifActive
    type(Tcocktail_descr) :: descr_aux_cockt
    logical :: ifSpecies


   
!    integer ::  nspeciesTmpShort, nspeciesTmpTransp, nspeciesTmpAero

    nSpeciesTransport = 0
    nSpeciesShortlived = 0
    nSpeciesAerosol = 0
    nReactRates = 0
    nullify(speciesTransport, speciesShortlived, speciesAerosol)
    !
    ! This array allows the transformation procedures to claim specific emission species
    ! to their ownership. This means that the corresponding transformation knows this
    ! species, knows what to do with it and how to treat the corresponding transport species.
    ! Note that the aerosol dynamics is a separate line, so the chemical part is (:,1) and 
    ! aerosol dynamics is in (:,2).
    ! The species unclaimed by anyone will be at the end added "as-is". That would mean that 
    ! they are just transported without any trasnformation - but with deposition as written 
    ! in the chemical database.
    !
    if(error)return
    iClaimedSpecies(1:nSpeciesEmission) = -1

#ifdef DEBUG 
        call msg("global_chemical_init got species emission", nSpeciesEmission)
        call report(speciesEmission)
        call msg("END species emission")
#endif



    !
    ! 1. if a specific ordering is needed, it is handled first.
    !
    do iTmp = 1, chemRules%nTransformations
      select case(chemRules%iTransformTypes(iTmp))
!!$      case (radioactive_flag)
!!$        call inventory(chemRules%rulesRadioactive, speciesTmpTransp, speciesTmpShort,&
!!$                     & nspeciesTmpShort, nspeciesTmp)
!!$        call addSpecies(speciesTmp, speciesTmp)
!!$        call addSpecies(speciesTmp2, speciesTmpShortlived)
!!$        exit

      case(transformation_cbm4)
        call init_chemicals_cbm4()
        call register_reaction_rates_cbm4(chemRules%rulesCBM4, nReactRates)
        call inventory_cbm4(chemRules%rulesCBM4, &
                     & speciesEmission, speciesTransport, speciesShortlived, speciesAerosol, &
                     & nSpeciesEmission, nSpeciesTransport, nSpeciesShortlived, nSpeciesAerosol, &
                     & iClaimedSpecies)
!        call addSpecies(speciesTransport, nSpeciesTransport, speciesTmp, nspeciesTmpTransp)
        call registerSpecies(chemRules%rulesCBM4, &
                           & speciesTransport, speciesShortlived, speciesAerosol, &
                           & nspeciesTransport, nspeciesShortlived, nspeciesAerosol)
      case(transformation_cbm4_SOA)
        call init_chemicals_cbm4_SOA()
        call inventory_cbm4_SOA(chemRules%rulesCBM4_SOA, &
                     & speciesEmission, speciesTransport, speciesShortlived, speciesAerosol, &
                     & nSpeciesEmission, nSpeciesTransport, nSpeciesShortlived, nSpeciesAerosol, &
                     & iClaimedSpecies)
!        call addSpecies(speciesTransport, nSpeciesTransport, speciesTmp, nspeciesTmpTransp)
        call registerSpeciescbm4_SOA(chemRules%rulesCBM4_SOA, &
                           & speciesTransport, speciesShortlived, speciesAerosol, &
                           & nspeciesTransport, nspeciesShortlived, nspeciesAerosol)                     

      case(transformation_cbm42_strato)
        call init_chemicals_cbm42_strato()
        call inventory_cbm42_strato(chemRules%rulesCBM42_strato, &
                     & speciesEmission, speciesTransport, speciesShortlived, speciesAerosol, &
                     & nSpeciesEmission, nSpeciesTransport, nSpeciesShortlived, nSpeciesAerosol, &
                     & iClaimedSpecies)
!        call addSpecies(speciesTransport, nSpeciesTransport, speciesTmp, nspeciesTmpTransp)
        call registerSpeciescbm42_strato(chemRules%rulesCBM42_strato, &
                           & speciesTransport, speciesShortlived, speciesAerosol, &
                           & nspeciesTransport, nspeciesShortlived, nspeciesAerosol)
      case(transformation_cbm42_strato_SOA)
        call init_chemicals_cbm42_strato_SOA()
        call inventory_cbm42_strato_SOA(chemRules%rulesCBM42_strato_SOA, &
                     & speciesEmission, speciesTransport, speciesShortlived, speciesAerosol, &
                     & nSpeciesEmission, nSpeciesTransport, nSpeciesShortlived, nSpeciesAerosol, &
                     & iClaimedSpecies)
!        call addSpecies(speciesTransport, nSpeciesTransport, speciesTmp, nspeciesTmpTransp)
        call registerSpeciescbm42_strato_SOA(chemRules%rulesCBM42_strato_SOA, &
                           & speciesTransport, speciesShortlived, speciesAerosol, &
                           & nspeciesTransport, nspeciesShortlived, nspeciesAerosol)
      case(transformation_cbm5_SOA)
        call init_chemicals_cbm5_SOA()
        call inventory_cbm5_SOA(chemRules%rulesCBM5_SOA, &
                     & speciesEmission, speciesTransport, speciesShortlived, speciesAerosol, &
                     & nSpeciesEmission, nSpeciesTransport, nSpeciesShortlived, nSpeciesAerosol, &
                     & iClaimedSpecies)
!        call addSpecies(speciesTransport, nSpeciesTransport, speciesTmp, nspeciesTmpTransp)
        call registerSpeciescbm5_SOA(chemRules%rulesCBM5_SOA, &
                           & speciesTransport, speciesShortlived, speciesAerosol, &
                           & nspeciesTransport, nspeciesShortlived, nspeciesAerosol)
      case(transformation_cbm5_strato_SOA)
        call init_chemicals_cbm5_strato_SOA()
        call inventory_cbm5_strato_SOA(chemRules%rulesCBM5_strato_SOA, &
                     & speciesEmission, speciesTransport, speciesShortlived, speciesAerosol, &
                     & nSpeciesEmission, nSpeciesTransport, nSpeciesShortlived, nSpeciesAerosol, &
                     & iClaimedSpecies)
!        call addSpecies(speciesTransport, nSpeciesTransport, speciesTmp, nspeciesTmpTransp)
        call registerSpeciescbm5_strato_SOA(chemRules%rulesCBM5_strato_SOA, &
                           & speciesTransport, speciesShortlived, speciesAerosol, &
                           & nspeciesTransport, nspeciesShortlived, nspeciesAerosol)
      case(int_missing)
        exit
      case default
        continue
      end select
  
    end do
    !
    ! 2. Then all the rest.
    !
    do iTmp = 1, chemRules%nTransformations

      ifActive = .true.
      select case(chemRules%iTransformTypes(iTmp))

      case(transformation_passive)
        call inventoryPassive(chemRules%rulesPassive, &
                     & speciesEmission, speciesTransport, speciesShortlived, speciesAerosol, &
                     & nSpeciesEmission, nSpeciesTransport, nSpeciesShortlived, nSpeciesAerosol, &
                     & iClaimedSpecies)
!        call addSpecies(speciesTransport, nspeciesTransport, speciesTmp, nspeciesTmpTransp)
        call registerSpeciesPassive(chemRules%rulesPassive, &
                           & speciesTransport, speciesShortlived, speciesAerosol, &
                           & nspeciesTransport, nspeciesShortlived, nspeciesAerosol)
      
      case(transformation_pollen)
        call inventoryPollen(chemRules%rulesPollen, &
                     & speciesEmission, speciesTransport, speciesShortlived, speciesAerosol, &
                     & nSpeciesEmission, nSpeciesTransport, nSpeciesShortlived, nSpeciesAerosol, &
                     & iClaimedSpecies)
        call registerSpeciesPollen(chemRules%rulesPollen, &
                           & speciesTransport, speciesShortlived, speciesAerosol, &
                           & nspeciesTransport, nspeciesShortlived, nspeciesAerosol)

      case(transformation_inert_PM)
        call inventoryPM(chemRules%rulesPM, &
                     & speciesEmission, speciesTransport, speciesShortlived, speciesAerosol, &
                     & nSpeciesEmission, nSpeciesTransport, nSpeciesShortlived, nSpeciesAerosol, &
                     & iClaimedSpecies)
        call registerSpeciesPM(chemRules%rulesPM, &
                           & speciesTransport, speciesShortlived, speciesAerosol, &
                           & nspeciesTransport, nspeciesShortlived, nspeciesAerosol)

      case(transformation_sulphur_dmat)
        call init_chemicals_dmat()
        call inventory_dmat(chemRules%rulesSulphurDMAT, &
                     & speciesEmission, speciesTransport, speciesShortlived, speciesAerosol, &
                     & nSpeciesEmission, nSpeciesTransport, nSpeciesShortlived, nSpeciesAerosol, &
                     & iClaimedSpecies, ifActive)
        call register_species_dmat(chemRules%rulesSulphurDMAT, &
                            & speciesTransport, speciesShortlived, speciesAerosol, &
                            & nspeciesTransport, nspeciesShortlived, nspeciesAerosol)

      case(transformation_acid_basic)
        call init_chemicals_acid_basic()
        call inventory_acid_basic(chemRules%rulesAcidBasic, &
                     & speciesEmission, speciesTransport, speciesShortlived, speciesAerosol, &
                     & nSpeciesEmission, nSpeciesTransport, nSpeciesShortlived, nSpeciesAerosol, &
                     & iClaimedSpecies)
        call registerSpeciesAcidBasic(chemRules%rulesAcidBasic, &
                           & speciesTransport, speciesShortlived, speciesAerosol, &
                           & nspeciesTransport, nspeciesShortlived, nspeciesAerosol)

      case(transformation_cbm4, transformation_cbm4_SOA, transformation_cbm42_strato, transformation_cbm42_strato_SOA, &
         & transformation_cbm5_SOA, transformation_cbm5_strato_SOA)
        ! Handled above.
        continue

      case(transformation_radioactive)
        call init_radioactive(chemRules%rulesRadioactive, &
                     & speciesEmission, speciesTransport,&
                     & nSpeciesEmission, nSpeciesTransport, &
                     & iClaimedSpecies, ifActive, timestep, timestep_output)
      case(int_missing)
        exit

      case default
        call msg('All transformations in the run:',chemRules%iTransformTypes(1:chemRules%nTransformations))
        call msg('Transformation type: ', chemRules%iTransformTypes(iTmp))
        call set_error('Strange transformation', 'global_chemical_init')
        return
      end select
      if (error) return
      !
      ! If the chemical module decides that it has nothing to do, it can switch itself out
      !
      if(.not. ifActive) chemRules%iTransformTypes(iTmp) = int_missing
#ifdef DEBUG
call msg('Transport species after chemistry: ',nSpeciesTransport, chemRules%iTransformTypes(iTmp))
call report(speciesTransport)
call msg('End transport species')
#endif
    end do  ! iTransformations

!    !
!    ! If for whatever reasons the chemical transformation modules decide that they cannot do anything
!    ! they can switch themselves off.
!    !
!    call compress_int_array(chemRules%iTransformTypes, int_missing)
!    if(error)return

    !
    ! Finally, the aerosol dynamics
    !
    do iTmp = 1, chemRules%nAerosolDynamics

      select case(chemRules%iAerosolDynTypes(iTmp))
        case(aerosol_dynamics_basic)
          call full_species_lst_4_ADB(chemRules%rulesAerDynBasic, &
                       & speciesEmission, speciesTransport, speciesShortlived, speciesAerosol, &
                       & nSpeciesEmission, nSpeciesTransport, nSpeciesShortlived, nSpeciesAerosol, &
                       & iClaimedSpecies)
!          call addSpecies(speciesTransport, nSpeciesTransport,  speciesTmp, nspeciesTmpTransp)
!          call addSpecies(speciesAerosol, nSpeciesAerosol, speciesTmp3, nspeciesTmpAero)
!          call registerSpecies(chemRules%rulesAerDynBasic, &
!                             & speciesEmission, speciesTransport, speciesShortlived, speciesAerosol, &
!                             & nSpeciesEmission, nSpeciesTransport, nSpeciesShortlived, nSpeciesAerosol, &
!                             & iClaimedSpecies)

        case(aerosol_dynamics_simple)
          call full_spec_lst_4_AerDynSimple(chemRules%rulesAerDynSimple, &
                       & speciesEmission, speciesTransport, speciesShortlived, speciesAerosol, &
                       & nSpeciesEmission, nSpeciesTransport, nSpeciesShortlived, nSpeciesAerosol, &
                       & iClaimedSpecies)
!          call addSpecies(speciesTransport, nSpeciesTransport,  speciesTmp, nspeciesTmpTransp)
          call registerSpecies_4_AerDynSimple(chemRules%rulesAerDynSimple, &
                             & speciesTransport, speciesShortlived, &
                             & nSpeciesTransport, nSpeciesShortlived)

        case(aerosol_dynamics_Mid_Atmosph)
          call full_spec_lst_4_AerDynMidAtm(chemRules%rulesAerDynMidAtmosph, &
                        & speciesEmission, speciesTransport, speciesShortlived, speciesAerosol, &
                        & nSpeciesEmission, nSpeciesTransport, nSpeciesShortlived, nSpeciesAerosol, &
                        & iClaimedSpecies)
          call registerSpecies(chemRules%rulesAerDynMidAtmosph, &
                             & speciesTransport, speciesShortlived, speciesAerosol, &
                             & nSpeciesTransport, nSpeciesShortlived, nSpeciesAerosol) !, &
!                             & passenger_tickets)
        
        case (aerosol_dynamics_VBS)
          call full_spec_lst_4_AerDynVBS(chemRules%rulesAerDynVBS, &
                        & speciesEmission, speciesTransport, speciesShortlived, speciesAerosol, &
                        & nSpeciesEmission, nSpeciesTransport, nSpeciesShortlived, nSpeciesAerosol, &
                        & iClaimedSpecies)
          call registerSpecies_4_AerDynVBS(chemRules%rulesAerDynVBS, &
                             & speciesTransport, speciesShortlived, &
                             & nSpeciesTransport, nSpeciesShortlived)
        
        
      case(int_missing)
        exit

      case default
        call msg('All aerosol dynamics in the run:', chemRules%iAerosolDynTypes(1: chemRules%nAerosolDynamics))
        call msg('Strange aerosol dynamic type: ', chemRules%iAerosolDynTypes(iTmp))
        call set_error('Strange aerosol dynamics', 'global_chemical_init')
        return
      end select
      if (error) return
#ifdef DEBUG
call msg('Transport species after aerosol dynamics: ',nSpeciesTransport, chemRules%iAerosolDynTypes(iTmp))
do stat = 1, nSpeciesTransport
  call report(speciesTransport(stat))
end do
#endif
    end do   ! aerosol dynamics

    !
    ! 3. Include emitted species not required by any transformation.
    !
    do iTmp = 1, nSpeciesEmission
      if(iClaimedSpecies(iTmp) < 0) &
           & call addSpecies(speciesTransport, nSpeciesTransport, (/speciesEmission(iTmp)/), 1)
      if(error)return
    end do
!!#ifdef DEBUG
!!call msg('Transport species after inclusion emission species: ',nSpeciesTransport, chemRules%iAerosolDynTypes(iTmp))
!!do stat = 1, nSpeciesTransport
!!  call report(speciesTransport(stat))
!!end do
!!#endif

    ! 
    ! 3b. Include auxiliary species after user input
    do iTmp=1,chemRules%nAuxCocktails
      call msg("Adding AUX cocktail: "//trim(chemRules%auxCocktName(iTmp)))
      call set_cocktail_description(chemRules%auxCocktName(iTmp), descr_aux_cockt, ifSpecies)
      if (error) return
      call get_inventory(descr_aux_cockt, species_aux, nspecies_aux)
      if (error) return
      if (fu_fails(nspecies_aux > 0, 'Aux Cocktail given but no species found', 'global_chemical_init')) return 
      call addSpecies(speciesTransport, nSpeciesTransport, species_aux, nspecies_aux)
      if (error) return
      deallocate(species_aux)
    enddo

    ! 4. We now know the final number of transport and shortlived
    ! species and we can allocate the lists.

    if (nSpeciesTransport == 0) then
      call set_error('No transport species', 'global_chemical_init')
      return
    else
      call msg('')
      call msg('Global chemical init. Transport species: ',nSpeciesTransport)
      do iTmp = 1, nSpeciesTransport
        call report(speciesTransport(iTmp))
      end do
    end if

    !
    ! In case the aerosol species are present, we have to raise the flag in chemical rules
    !
    chemRules%rulesAerosol%ifNumberCnc = (nSpeciesAerosol > 0)
    
!!$    if (nSpeciesShortlived > 0) then
!!$      allocate(speciesShortlived(nSpeciesShortlived), stat=stat)
!!$      if (stat /= 0) then
!!$        call set_error('Allocate failed', 'global_chemical_init')
!!$        return
!!$      end if
!!$      speciesShortlived(1:nSpeciesShortlived) = speciesTmpShortlived(1:nSpeciesShortlived)
!!$    end if

    ! Knowing emission and transport species, we can set the according
    ! adaptor in chemRunSetup
    ! 
    call msg('Mapping emission species to transport ones')

    call create_mode_projection(speciesEmission, nSpeciesEmission, &    ! emission
                              & speciesTransport,nSpeciesTransport, &   ! transport
                              & speciesAerosol, nSpeciesAerosol, &      ! aerosol
                              & chemRules%ChemRunSetup%refEmis2Transp_mass, &   ! mapping for mass emission
                              & chemRules%ChemRunSetup%refEmis2Transp_nbr, &    ! mapping for number emission
                              & .true.)        ! if full-spectrum check
    if(error)return
    call msg('Done emission to transport mapping')

    if (chemRules%ifPhotoAOD) &
      call init_photoatt_lut(speciesTransport, nSpeciesTransport, chemRules%photoAODwavelength)




    !----------------------------------------------------------------------------------
    !
    ! Upon creating the transformation structures, we can initialise the deposition ones
    !
    indDepositionType(1:nSpeciesTransport+1) = int_missing
    if(error)return

    do iTmp = 1, chemRules%rulesDeposition%nDepositionTypes
    
      select case(chemRules%rulesDeposition%iDepositionType(iTmp))

        case(transformation_passive, transformation_inert_PM, transformation_sulphur_dmat)
!          call init_deposition_<whatever>(speciesTransport, indDepositionType, nSpeciesTransport, &
!                                     & chemRules%rulesPassive)


        case(aerosol_dynamics_basic, aerosol_dynamics_simple, aerosol_dynamics_VBS)
          call set_error('Aerosol dynamics does not have own deposition','global_chemical_init')

        case(deposition_standard)
          call init_standard_deposition(speciesTransport, indDepositionType, nSpeciesTransport, &
                                      & chemRules%rulesDeposition)

        case default
          call msg('Transformation not supported:',chemRules%rulesDeposition%iDepositionType(iTmp))
          call set_error('Transformation not supported','global_chemical_init')
      end select
      if(error)return
    end do  ! nDepositions
    
  end subroutine global_chemical_init


  !***************************************************************************

  subroutine add_transformation_input_needs(chemRules, &
                                          & q_met_dyn, q_met_stat, q_disp_dyn, q_disp_stat, &
                                          & meteo_input_global)
    !
    ! Collects the input data requests from all chemical modules and
    ! returns the meteorological part in q_static and q_dynamic arrays.
    ! This switch is for historical reasons and in general the Tmeteo_input
    ! structures are more universal. 
    !
    implicit none

    ! Imported parameters
    type(Tchem_rules), intent(in) :: chemRules
    integer, dimension(:), intent(inout) :: q_met_dyn, q_met_stat, q_disp_dyn, q_disp_stat
    type(Tmeteo_input), intent(out) :: meteo_input_global
    
    ! Local variables
    integer :: iTmp, iType, nMeteoInputs, iMetInput
    logical :: ifNeedPhotolysis

    if(.not.(chemRules%defined == silja_true))then
      call set_error('Undefined chemistry rules','transformation_input_needs')
      return
    endif

    iMetInput = 0

    ! The last one for wet deposition
    nMeteoInputs = chemRules%nTransformations + chemRules%nAerosolDynamics + 1
    if (chemRules%need_photo_lut) then
        nMeteoInputs = nMeteoInputs + 1
      if (chemRules%ifPhotoAOD) nMeteoInputs = nMeteoInputs + 1
    endif

    allocate(pMeteo_input_local(nMeteoInputs), stat=iTmp)

    if(iTmp /= 0)then
      call set_error('Failed memory allocation for meteo_input_local','transformation_input_needs')
      return
    endif

    do iType = 1, chemRules%nTransformations

      select case(chemRules%iTransformTypes(iType))

        case(transformation_passive)
          call passive_proc_input_needs(chemRules%rulesPassive, pMeteo_input_local(iType))
                
        case(transformation_pollen)
          call pollen_proc_input_needs(chemRules%rulesPollen, pMeteo_input_local(iType))

        case(transformation_inert_PM)
          call pm_proc_input_needs(chemRules%rulesPM, pMeteo_input_local(iType))

        case(transformation_sulphur_dmat)
          call sulphur_dmat_input_needs(chemrules%rulesSulphurDMAT, pMeteo_input_local(iType))

        case(transformation_acid_basic)
          call acidBasic_input_needs(chemRules%rulesAcidBasic, pMeteo_input_local(iType))

        case(transformation_cbm4)
          call cbm4_input_needs(chemRules%rulesCBM4, pMeteo_input_local(iType))
          
        case(transformation_cbm4_SOA)
          call cbm4_SOA_input_needs(chemRules%rulesCBM4_SOA, pMeteo_input_local(iType))
 

        case(transformation_cbm42_strato)
          call cbm42_strato_input_needs(chemRules%rulesCBM42_strato, pMeteo_input_local(iType))

        case(transformation_cbm42_strato_SOA)
          call cbm42_strato_SOA_input_needs(chemRules%rulesCBM42_strato_SOA, pMeteo_input_local(iType))

        case(transformation_cbm5_SOA)
          call cbm5_SOA_input_needs(chemRules%rulesCBM5_SOA, pMeteo_input_local(iType))
          
        case(transformation_cbm5_strato_SOA)
          call cbm5_strato_SOA_input_needs(chemRules%rulesCBM5_strato_SOA, pMeteo_input_local(iType))
          
        case(transformation_radioactive)
          call radioactive_input_needs(chemRules%rulesRadioactive, pMeteo_input_local(iType))

!!$        case(species_persistent_organics)
!!$          call POP_cocktail_input_needs(chemRules%rulesPOP, chemRules%rulesAerosol, &
!!$                                      & q_dynamicTmp, q_staticTmp)
        case(int_missing)
          ! if transformation was claimed inactive by init_XXX (e.g. init_radioactive) 
          ! chemRules%iTransformTypes(iType) is set to int_missing
          cycle

        case default
          call msg('Unknown transformationi: no, type:',iType, chemRules%iTransformTypes(iType))
          call set_error('Unknown transformation type','transformation_input_needs')
          return
      end select

    end do  ! active transformations

    !
    ! The same for aerosol dynamics
    !
    iMetInput = chemRules%nTransformations  !! Assign it in case nAerosolDynamics = 0
    do iType = 1, chemRules%nAerosolDynamics
      iMetInput = iMetInput + 1
      select case(chemrules%iAerosolDynTypes(iType))

        case(aerosol_dynamics_basic)
          call AerDynBasic_input_needs(chemRules%rulesAerDynBasic, pMeteo_input_local(iMetInput))

        case(aerosol_dynamics_simple)
          call AerDynSimple_input_needs( pMeteo_input_local(iMetInput))

        case(aerosol_dynamics_Mid_Atmosph)
          call AerDynMidAtm_input_needs(pMeteo_input_local(iMetInput))
          
                          
        case(aerosol_dynamics_VBS)
            call AerDynVBS_input_needs(pMeteo_input_local(iMetInput))
            
        case default
          call msg('Unknown aerosol dynamics type:',iType)
          call set_error('Unknown aerosol dynamics type','transformation_input_needs')
          return
      end select

    end do  ! Active aerosol dynamics
   
    !
    ! Check if photolysis is needed
    !
    if (chemRules%need_photo_lut) then
      iMetInput = iMetInput +1
      call photolysis_input_needs(chemRules%useDynamicAlbedo, chemRules%photoO3col == meteo_column, &
            &chemRules%cloud_model_for_photolysis,  pMeteo_input_local(iMetInput))

      if (chemRules%ifPhotoAOD) then
        iMetInput = iMetInput +1
        call  optical_column_input_needs(chemRules%rulesOpticDens, pMeteo_input_local(iMetInput))
      endif
    endif

    iMetInput = iMetInput +1
    call wet_deposition_input_needs(chemRules%rulesDeposition, pMeteo_input_local(iMetInput))

    !
    ! Having all the local requests collected, we can define the global meteo_input structure
    !
    call define_meteo_input(meteo_input_global, pMeteo_input_local, nMeteoInputs)
    if(error)return

    !
    ! Now we have to turn the nice structure to the arrays, as requested from the dispersion_models
    !
    do iType = 1, meteo_input_global%nQuantities

      if(meteo_input_global%q_type(iType) == meteo_dynamic_flag)then
        iTmp = fu_merge_integer_to_array(meteo_input_global%quantity(iType), q_met_dyn)

      elseif(meteo_input_global%q_type(iType) == meteo_single_time_flag)then
        iTmp = fu_merge_integer_to_array(meteo_input_global%quantity(iType), q_met_stat)

      elseif(meteo_input_global%q_type(iType) == dispersion_dynamic_flag)then
        iTmp = fu_merge_integer_to_array(meteo_input_global%quantity(iType), q_disp_dyn)

      elseif(meteo_input_global%q_type(iType) == dispersion_single_time_flag)then
        iTmp = fu_merge_integer_to_array(meteo_input_global%quantity(iType), q_disp_stat)

      else
        call msg('Unknown type of input data:',meteo_input_global%q_type(iType))
        call set_error('Unknown type of input data','add_transformation_input_needs')
        return
      endif

    end do  ! cycle over requested quantities

  end subroutine add_transformation_input_needs


  !************************************************************************************
  
  subroutine get_tla_chemistry(chemrules, traj)
    !
    ! Collect transformations' requests for tangent linearization. If tangent linear
    ! is not allowed, corresponding subroutine sets error. 
    !
    implicit none
    type(Tchem_rules), intent(in) :: chemrules
    type(t_tla_trajectory), intent(inout) :: traj

    integer :: iType
    logical :: required, need_hstart

    need_hstart = .false.
    !
    ! Scan transformation types - some of them need TLA and other stuff stored
    !
    do iType = 1, chemrules%nTransformations
      select case(chemrules%iTransformTypes(iType))
        case(transformation_passive)
          required = fu_if_tla_required(chemrules%rulesPassive)
        case(transformation_inert_pm)
          required = fu_if_tla_required(chemrules%rulesPM)
        case(transformation_sulphur_dmat)
          required = fu_if_tla_required(chemrules%rulesSulphurDMAT)
        case(transformation_acid_basic)
          required = fu_if_tla_required(chemrules%rulesAcidBasic)
        case(transformation_cbm4)
          required = fu_if_tla_required(chemrules%rulesCbm4)
          need_hstart = .true.
        case(transformation_cbm4_SOA)
          required = fu_if_tla_required(chemrules%rulesCbm4_SOA)
          need_hstart = .true.
        case(transformation_cbm42_strato)
          required = fu_if_tla_required(chemrules%rulesCbm42_strato)
          need_hstart = .true.
        case(transformation_cbm42_strato_SOA)
          required = fu_if_tla_required(chemrules%rulesCBM42_strato_SOA)
          need_hstart = .true.
        case(transformation_cbm5_SOA)
          required = fu_if_tla_required(chemrules%rulesCBM5_SOA)
          need_hstart = .true.
        case(transformation_cbm5_strato_SOA)
          required = fu_if_tla_required(chemrules%rulesCBM5_strato_SOA)
          need_hstart = .true.
        case(transformation_radioactive)
          required = fu_if_tla_required(chemrules%rulesRadioactive)
        case(int_missing)
          exit ! no more transformations
        case default
          call set_error('Unsupported transformation', 'get_tladj_chemistry')
      end select      ! transform type
      if (error) return
      if (required) call add_tla_traj(chemrules%iTransformTypes(iType), traj)
      ! the storage for the h_start array
      if (need_hstart) call add_tla_traj(-chemrules%iTransformTypes(iType), traj, dimensions=3)
    end do  ! transformation types
    !
    ! The same for aerosol dynamics
    !
    required = .false.
    do iType = 1, chemRules%nAerosolDynamics
      select case(chemrules%iAerosolDynTypes(iType))
        case(aerosol_dynamics_basic)
!          required = fu_if_tla_required(chemRules%rulesAerDynBasic)
        case(aerosol_dynamics_simple)
          required = fu_if_tla_required(chemRules%rulesAerDynSimple)
        case(aerosol_dynamics_Mid_Atmosph)
!          required = fu_if_tla_required(chemRules%rulesAerDynMAAD)
        case(aerosol_dynamics_VBS)
!          required = fu_if_tla_required(chemRules%rulesAerDynVBS)
        case default
          call msg('Unknown aerosol dynamics type:',iType)
          call set_error('Unknown aerosol dynamics type','get_tla_chemistry')
          return
      end select
      if (error) return
      if (required) call add_tla_traj(chemrules%iAerosolDynTypes(iType), traj)
    end do  ! Active aerosol dynamics
    !
    ! Same for deposition: need TLA to store?
    !
    required = fu_if_TLA_required(chemrules%rulesDeposition)
    if (error) return
    if (required) call add_tla_traj(chemrules%iTransformTypes(iType), traj)

  end subroutine get_tla_chemistry

  !************************************************************************************

  subroutine transform_maps(mapTransport, mapShortLived, mapAerosol, mapDryDep, mapWetDep, &
                          & mapReactRates, &
                          & pBoundaryBuffer, &
                          & garbage, tla_step, &
                          & met_buf, disp_buf, meteo_input, &
                          & pHorizInterpStruct, pVertInterpStruct, &
                          & ifHorizInterp, ifVertInterp, &
                          & ifDryDep_cumulative_in_output, ifWetDep_cumulative_in_output, &
                          & chemRules, seconds, now)
    !
    ! Main chemistry and aerosol transformation routine fo Eulerian environment
    !
    implicit none
    type(Tmass_map), intent(inout) :: mapTransport, mapShortlived, mapAerosol, mapDryDep, mapWetDep, mapReactRates
    type(TboundaryBuffer), pointer :: pBoundaryBuffer
    real, dimension(:,:), intent(inout) :: garbage
    type(t_tla_step), intent(inout) :: tla_step
    type(Tfield_buffer), pointer :: met_buf, disp_buf
    type(Tmeteo_input), intent(in) :: meteo_input
    type(THorizInterpStruct), intent(in) :: pHorizInterpStruct
    type(TVertInterpStruct), intent(in) :: pVertInterpStruct
    logical, intent(in) :: ifHorizInterp, ifVertInterp
    type(silja_logical), intent(in) :: ifDryDep_cumulative_in_output, ifWetDep_cumulative_in_output
    type(Tchem_rules), intent(inout) :: chemRules
    type(silja_time), intent(in) :: now
    real, intent(in) :: seconds

    integer :: ix, iy, i3d, iSrc, iTransf, iSpecies, status, ind_dz, i1d, cbm_type, itransf_cbm, &
             & nReact
    real, dimension(:,:), pointer :: cncTrn
    real, dimension(:,:), pointer :: garb_array, photorates, reactRates, metdat_col, o3column
    real, dimension(:), pointer :: soot_col, pwc_col, tau_above_bott
    real :: cell_volume, zenith_cos, lat, lon, dz
    real, dimension(:), pointer :: cell_size_x, cell_size_y, dz_past, dz_future, cncAer, cncSL, &
         & aodext,  aodscat, metdat, cell_mass_past, cell_mass_future
    logical :: have_cb4, have_tla, print_it, ifNeeded
    real, dimension(:,:,:,:), pointer :: tla_point_cb4 => null(), tla_point_hstart => null()
    integer :: ithread, iLev, iSoot, istat
    real :: fixed_albedo, ssa 
    logical, save :: if_first = .true.
    character(len=*), parameter :: sub_name = 'transform_maps'
    integer, save :: indO3, n_soot
    real :: mass_air, mass_ones !!Masses of ir and of ones for the cell
    integer, dimension(max_species), save :: ind_soot
    integer, save :: ind_air_mass, iSpOnes
    real, dimension(max_species), save :: frac_soot
    real, dimension(max_species), save :: diam, density

#ifdef DEBUG_MORE
     !Temporary arrays to save the state before chemistry
    real, dimension (max_species) :: cnctrn_tmp,  cnctrn_tmp1,  cnctrn_tmp2,  cnctrn_tmp3 !FIXME
    real :: fTmp
#endif
    character(len=180) :: chTmp

    if (if_first) then
      call allocate_scav_amount(mapTransport)
      call data_for_photolysis_and_cld_model(mapTransport, chemRules, indO3, ind_soot, &
           & frac_soot, n_soot, density, diam)
      if (chemRules%ifOnesAdjust) then
        call msg("Chemistry will use ones tracer as cell size")
        iSpOnes = select_single_species(mapTransport%species, mapTransport%nSpecies, &
                          & 'ones', in_gas_phase, real_missing)
        if (error .or. (iSpOnes < 1)) then
          call set_error("Couldn't find 'ones' species", sub_name)
          return
        endif
      endif
      if_first = .false.
    end if

    !call msg_warning('Skipping chemistry')
    !return

    if (chemRules%nTransformations + chemRules%nAerosolDynamics == 0 .and. &
         & chemRules%rulesDeposition%scavengingType == scavNoScav) return

    if (.not. ChemStuff%defined == silja_true) then
      if (defined(mapReactRates)) then
        nReact = mapReactRates%nSpecies
      else
        nReact = 100   !Whatever size
                             ! Have no clue wher to get num_reactions
                                               ! for photorates, took 100
      endif
       call init_chemisryStuff(mapTransport%nSrc, mapTransport%nSpecies, &
                             & mapAerosol%nSpecies, mapShortLived%nSpecies, &
                             & nReact,  meteo_input%nQuantities, & 
                             & mapTransport%n3d)
       if (error) return
    endif


    ind_dz = fu_index(disp_buf, cell_size_z_flag)
    if (ind_dz < 1) then
      call set_error('Failed to find the dz field', sub_name)
      return
    end if
    
    ind_air_mass = fu_index(disp_buf, disp_cell_airmass_flag)
    if (ind_air_mass < 1 .and. chemRules%ifOnesAdjust) then
      call set_error('Failed to find the disp_cell_airmass field', sub_name)
      return
    end if
    
    cell_size_x => fu_grid_data(dispersion_cell_x_size_fld)
    cell_size_y => fu_grid_data(dispersion_cell_y_size_fld)

    if (fu_index(transformation_cbm4, chemrules%iTransformTypes) > 0) then
      cbm_type = transformation_cbm4
    elseif (fu_index(transformation_cbm4_SOA, chemrules%iTransformTypes) > 0) then
      cbm_type = transformation_cbm4_SOA
    else if (fu_index(transformation_cbm42_strato, chemrules%iTransformTypes) > 0) then
      cbm_type = transformation_cbm42_strato
    else if (fu_index(transformation_cbm42_strato_SOA, chemrules%iTransformTypes) > 0) then
      cbm_type = transformation_cbm42_strato_SOA
    else if (fu_index(transformation_cbm5_SOA, chemrules%iTransformTypes) > 0) then
      cbm_type = transformation_cbm5_SOA
    else if (fu_index(transformation_cbm5_strato_SOA, chemrules%iTransformTypes) > 0) then
      cbm_type = transformation_cbm5_strato_SOA
    else
      cbm_type = int_missing
    end if
    itransf_cbm = fu_index(cbm_type, chemrules%iTransformTypes)
    have_cb4 = cbm_type /= int_missing .and. chemrules%iWhomToApplyTransform(itransf_cbm) == eulerian_flag

    if (have_cb4 .and. .not. allocated(cb4_h_start_array)) then
      ! Speedup for (or other KPP schemes): use the initial step length stored at previous
      ! step for each grid cell.
      allocate(cb4_h_start_array(mapTransport%n3d, mapTransport%nx, mapTransport%ny), &
             & stat=status)
      if(fu_fails(status == 0, 'Allocate failed', sub_name))return
      cb4_h_start_array(:,:,:) = -1.0
    end if

    if (have_cb4 .and. defined(tla_step)) then
      call msg('Have tla step')
      ! More speedup: instead of requiesting the linearization point inside the loop,
      ! store it here if needed.
      tla_point_cb4 => fu_get_tla_point(tla_step, cbm_type)
      tla_point_hstart => fu_get_tla_point(tla_step, -cbm_type)
    else
      nullify(tla_point_cb4)
    end if
    have_tla = associated(tla_point_cb4) .and. associated(tla_point_hstart)
    if (seconds < 0 .and. have_cb4 .and. .not. have_tla) then
      call set_error('Have cb4, timestep < 0 but no linearization trajectory', sub_name)
      return
    end if


    if (chemRules%useDynamicAlbedo) then 
       fixed_albedo = real_missing  !!Use it from meteo input
    else
       fixed_albedo = chemRules%defaultStaticAlbedo
    end if
    
    !-------------------------------------------------------------------------
    !
    ! The main cycle over the grid
    !

    !call msg('Total N before transformation:', get_total_n(mapTransport))

    !$OMP PARALLEL if (if_OMP_chemistry) DEFAULT(SHARED) PRIVATE(metdat, metdat_col, reactRates, ix, iy, i3d, itransf, isrc, print_it, &
    !$OMP & cell_volume, zenith_cos, lat, lon, dz_past, dz_future, dz, i1d, garb_array, cncTrn, mass_air, mass_ones, &
    !$OMP & photorates, aodext, aodscat, o3column, cncAer, cncSL, ithread, soot_col, pwc_col, tau_above_bott, ssa, iSoot)

    iThread = 0
    !$ iThread = OMP_GET_THREAD_NUM()
     
    metdat_col       => ChemStuff%arrStuff(iThread)%metdat(:,:) !!nMetInput, nz
    garb_array   => ChemStuff%arrStuff(iThread)%garb_array(:,:)
    cncTrn       => ChemStuff%arrStuff(iThread)%cncTrn(:,:)
    cncAer       => ChemStuff%arrStuff(iThread)%cncAer(:)
    cncSL        => ChemStuff%arrStuff(iThread)%cncSL(:)
    photorates   => ChemStuff%arrStuff(iThread)%photorates(:,:)
    aodext       => ChemStuff%arrStuff(iThread)%aodext(:)
    aodscat      => ChemStuff%arrStuff(iThread)%aodscat(:)
    o3column     => ChemStuff%arrStuff(iThread)%o3column(:,:)
    soot_col     => ChemStuff%arrStuff(iThread)%soot_col(:)
    pwc_col      => ChemStuff%arrStuff(iThread)%pwc_col(:)
    tau_above_bott      => ChemStuff%arrStuff(iThread)%tau_above_bott(:)
    reactRates   => ChemStuff%arrStuff(iThread)%reactRates(:,:)
    
    garb_array(1:mapTransport%nSrc, 1:mapTransport%nSpecies) = 0.0
    o3column(:,:) = real_missing !! (1:nLev,1:nSrc) should not be used uninitialized
    
    if (cbm_type == transformation_cbm4) call prepare_step_cbm4()
    if (cbm_type == transformation_cbm4_SOA) call prepare_step_cbm4_SOA()
    if (cbm_type == transformation_cbm42_strato) call prepare_step_cbm42_strato()
    if (cbm_type == transformation_cbm42_strato_SOA) call prepare_step_cbm42_strato_SOA()
    if (cbm_type == transformation_cbm5_SOA) call prepare_step_cbm5_SOA()
    if (cbm_type == transformation_cbm5_strato_SOA) call prepare_step_cbm5_strato_SOA()

    !$OMP DO collapse (2) schedule (guided)
    do iy = 1, mapTransport%ny
      do ix = 1, mapTransport%nx
        if (error) cycle

        i1d = (iy-1)*mapTransport%nx + ix

        call fill_in_meteo_input_column(meteo_input, metdat_col, &
                                 & met_buf, disp_buf, &
                                 & ix, iy, mapTransport%n3d, &
                                 & pHorizInterpStruct, pVertInterpStruct, &  ! meteo to disp
                                 & ifHorizInterp, ifVertInterp)

        lat = fu_lat_geographical_from_grid(real(ix), real(iy), dispersion_grid)
        lon = fu_lon_geographical_from_grid(real(ix), real(iy), dispersion_grid)
        if (chemRules%need_photo_lut .or. chemRules%rulesDeposition%scavengingType == scav2020) then

          if (chemRules%ifPhotoAOD) then
            call get_photoatt_aod(metdat_col, & 
                    & mapTransport%arM(1:mapTransport%nSpecies,1:mapTransport%nSrc,1:mapTransport%n3d,ix,iy),&
                    & aodext, aodscat)
            if (error) cycle
          endif

          if (chemRules%cloud_model_for_photolysis == detailed_cloud) then
            soot_col = 0.0
            do iSoot = 1, n_soot
              soot_col = soot_col + frac_soot(iSoot)*mapTransport%arM(ind_soot(iSoot), 1, 1:mapTransport%n3d, ix,iy) &
                   & / (cell_size_x(i1d)*cell_size_y(i1d))
            end do
          end if

          if (chemRules%cloud_model_for_photolysis == detailed_cloud .or. chemRules%rulesDeposition%scavengingType == scav2020) then
            call compute_water_cloud_properties(mapTransport%arM(1:mapTransport%nSpecies, 1:mapTransport%nSrc, 1:mapTransport%n3d, ix,iy), &
                 & metdat_col, aodext, aodscat, density, diam, soot_col, ix, iy, pwc_col, tau_above_bott, ssa)
          end if

          if (chemRules%need_photo_lut) then 
            !calculate the ozone column above (in Dobson units)
            if (chemRules%PhotoO3col == mass_map) then
              call get_o3column(metdat_col, mapTransport%arM(indO3,1:mapTransport%nSrc,1:mapTransport%n3d,ix,iy), o3column)
              !call msg('RISTO TEST O3: ix,iy,Ozone column: ', (/real(ix), real(iy), o3column(1,1) /) )
            endif
            if (error) cycle
            zenith_cos = fu_solar_zenith_angle_cos(lon, lat, now)
            if (error) cycle

            !NOTE: Only the first source is used for o3column when calculating the photorates!!!!
            call get_photorates_column(metdat_col, zenith_cos, fixed_albedo, now, aodext, &
                 & aodscat, o3column(1:mapTransport%n3d,1), photorates, chemRules%ifPhotoAOD, &
                 & chemRules%PhotoO3col, tau_above_bott, ssa, chemRules%cloud_model_for_photolysis)
          end if
          if (error) cycle
        end if

        if (seconds < 0) then
          ! WARNING! adjoint scavenging makes some noncense. Check for pollen leads to substantially
          ! different mass in air with the below hack and without. Unless the check passes
          ! the hack should be kept...
          !  HACK here: 
          ! Apply forward scavenging here (just flip he sign of "seconds")
          ! Should be fine for particles and very soluble gases, and for non-soluble gases
          ! Medium-solubility gasas and SO2 did not work anyway...
           call scavenge_column(mapTransport, mapWetDep, garb_array, tla_step, -seconds, &
                & chemRules%rulesDeposition, chemRules%low_mass_trsh, metdat_col, ix, iy, pwc_col, &
                & pHorizInterpStruct, pVertInterpStruct, ifHorizInterp, ifVertInterp, met_buf)                                                                                                                                                                                                                 
           if(error)call set_error('Trouble with scavenging', sub_name)
        end if

        do i3d = 1, mapTransport%n3d
          
          if (error) cycle

          print_it = .false. !ix == 30 .and. iy == 55 .and. now >= fu_set_time_utc(2009,12,19,20,0,0.)


          metdat => metdat_col(:,i3d)
          dz_past => disp_buf%p4d(ind_dz)%past%p2d(i3d)%ptr
          dz_future => disp_buf%p4d(ind_dz)%future%p2d(i3d)%ptr
          dz = met_buf%weight_past*dz_past(i1d) + (1.0-met_buf%weight_past)*dz_future(i1d)
          cell_volume = dz * cell_size_x(i1d) * cell_size_y(i1d)

          if (chemRules%ifOnesAdjust) then
            mass_air = disp_buf%p4d(ind_air_mass)%past%p2d(i3d)%ptr(i1d) * met_buf%weight_past + &
                     & disp_buf%p4d(ind_air_mass)%future%p2d(i3d)%ptr(i1d) * (1.0-met_buf%weight_past)

            mass_ones = mapTransport%arM(iSpOnes,1,i3d,ix,iy)
            if (.not. ( mass_ones > 0 )) then
              call set_error("Strane mass_ones="//fu_str(mass_ones)//" at (i3d,ix,iy)=("//&
                 &fu_str(i3d)//","//fu_str(ix)//","//fu_str(iy)//")", sub_name)
              cycle
            endif
            cell_volume = cell_volume * mass_ones / mass_air !! more ones -- biger cell
          endif

          !
          ! Masses are to be turned into concentrations
          !
          cncTrn(1:mapTransport%nSpecies,1:mapTransport%nSrc) = &
                         & mapTransport%arM(1:mapTransport%nSpecies,1:mapTransport%nSrc,i3d,ix,iy) / &
                         & cell_volume
          cncSL(1:mapShortLived%nSpecies) = 0.0
!          mapTransport%arM(1:mapTransport%nSpecies,1:mapTransport%nSrc,i3d,ix,iy) = &
!                         & mapTransport%arM(1:mapTransport%nSpecies,1:mapTransport%nSrc,i3d,ix,iy) / &
!                         & cell_volume
          garb_array(1:mapTransport%nSrc, 1:mapTransport%nSpecies) = &
                                & garb_array(1:mapTransport%nSrc, 1:mapTransport%nSpecies) / cell_volume

          do iSrc = 1, mapTransport%nSrc

if (print_it) then
  !   if (print_it) call ooops("here")
call msg('Dump, x,y,z,Src=('+fu_str(ix)+','+fu_str(iy)+','+fu_str(i3d), isrc)
do iSpecies = 1, meteo_input%nQuantities
  call msg(fu_quantity_string(meteo_input%quantity(iSpecies)))
enddo
call msg('metdat:',metdat(1:meteo_input%nQuantities))
do iSpecies = 1, mapTransport%nSpecies
  call msg('Mass before all calls:' + &
         & fu_str(mapTransport%species(iSpecies)), &
         & cncTrn(iSpecies,isrc))
end do
endif
            
            ! Aerosol dynamics, adjoint ===============================================================
            
            if (seconds < 0) then
              do iTransf = 1, chemRules%nAerosolDynamics
                ! skip lagrangian-only transform:
                if(.not. fu_if_eulerian_present(chemrules%iWhomToApplyAerDyn(iTransf)))cycle 
                
#ifdef DEBUG_MORE
                call check_mass_vector(cnctrn(1:maptransport%nspecies,isrc), &
                                     & garb_array(isrc,:), &
                                     & maptransport%species, 'chem manager, before adj ad:'+fu_str(itransf), &
                                     & chemrules%low_cnc_trsh, &
                                     & maptransport%nspecies, ix, iy, i3d, print_it)
                if (error) continue
#endif
                
                select case(chemrules%iAerosolDynTypes(iTransf))
                case(aerosol_dynamics_basic)
                  call set_error('aerosol_dynamics_basic is not available for adjoint', sub_name)
                  
                  
!                  cncAer(1:mapAerosol%nSpecies) = &   ! used only here
!                               & mapAerosol%arM(1:mapAerosol%nSpecies,isrc,i3d,ix,iy) / cell_volume
!!                  cncSL(1:mapShortLived%nSpecies) = 0.
!                  call transform_AerDynBasic(chemRules%rulesAerDynBasic, &
!                                           & cncTrn(1:mapTransport%nSpecies,isrc), &
!                                           & cncSL(1:mapShortLived%nSpecies), &
!                                           & cncAer(1:mapAerosol%nSpecies), &
!                                           & chemRules%low_cnc_trsh, &
!                                           & garb_array(isrc,:), &
!                                           & chemRules%ChemRunSetup%mapVolume2NumberCnc, &
!                                           & metdat, &
!                                           & seconds, &
!                                           & print_it)
!                  mapAerosol%arM(1:mapAerosol%nSpecies,isrc,i3d,ix,iy) = &
!                                                          & cncAer(1:mapAerosol%nSpecies) * cell_volume
!
                case(aerosol_dynamics_simple)
                  call set_error('aerosol_dynamics_simple is not available for adjoint', sub_name)
                  
                  !zenith_cos = fu_solar_zenith_angle_cos(lon, lat, now)
                  !call transform_AerDynSimple(cncTrn(1:mapTransport%nSpecies,isrc), &
                  !                          & cncSL(1:mapShortLived%nSpecies), &
                  !                          & garb_array(isrc,:), &
                  !                          & chemRules%rulesAerDynSimple, chemRules%low_cnc_trsh, &
                  !                          & metdat, seconds, zenith_cos, print_it)

                case(aerosol_dynamics_Mid_Atmosph)
                  call set_error('aerosol_dynamics_Mid_Atmosph is not available for adjoint', sub_name)
                  
                  !call transform_AerDynMidAtm(cncTrn(1:mapTransport%nSpecies,isrc), &
                  !                          & cncSL(1:mapShortLived%nSpecies), &
                  !                          & cncAer(1:mapAerosol%nSpecies), &
                  !                          & garb_array(isrc,:), &
                  !                          & chemRules%rulesAerDynMidAtmosph, &
                  !                          & chemRules%low_cnc_trsh, &
                  !                          & metdat, &
                  !                          & seconds, &
                  !                          & mapTransport%nSpecies, &
                  !                          & mapTransport%Species, ix, iy, i3d, &
                  !                          & print_it)                         

                case(aerosol_dynamics_VBS)
                  call set_error('aerosol_dynamics_VBS is not available for adjoint', sub_name)
                  
                  !call transform_AerDynVBS(cncTrn(1:mapTransport%nSpecies,isrc), &
                  !                       & cncSL(1:mapShortLived%nSpecies), &
                  !                       & garb_array(isrc,:), &
                  !                       & chemRules%rulesAerDynVBS, chemRules%low_cnc_trsh, &
                  !                       & metdat, seconds, print_it)
                case(int_missing)
                  continue

                case default
                  call msg('Unknown aerosol dynamics:',chemrules%iAerosolDynTypes(itransf))
                  call set_error('Unknown aerosol dynamics', sub_name)
                  cycle
                end select
              enddo    ! nAerosolDynamics
            else
              ! forward run
              if (debug_level > 0) then
                call check_mass_vector(cncTrn(1:mapTransport%nSpecies,isrc), &
                                     & garb_array(iSrc,:), &
                                     & mapTransport%species, 'chem manager, before chemistry', &
                                     & chemRules%low_cnc_trsh, &
                                     & mapTransport%nSpecies, ix, iy, i3d, print_it)
              end if
            end if   ! if seconds < 0

            !
            ! Chemical transformations, forward and adjoint =================================================
            !
            do iTransf = 1, chemRules%nTransformations
              ! skip lagrangian-only transfomr:
              if(.not. fu_if_eulerian_present(chemrules%iWhomToApplyTransform(iTransf)))cycle 
#ifdef DUMP_MORE
call msg('before the chemistry call:' + fu_str(chemrules%iTransformTypes(itransf)))
do iSpecies = 1, mapTransport%nSpecies
  call msg('Mass before the next call:' + &
         & fu_str(mapTransport%species(iSpecies)), &
         & cncTrn(iSpecies,isrc))
end do
#endif
              select case(chemrules%iTransformTypes(itransf))
                case(transformation_passive)
                  call transform_passive(cncTrn(1:mapTransport%nSpecies,isrc), &
                                       & chemRules%rulesPassive, &
                                       & metdat, i3d, &
                                       & seconds, print_it)
                                     
                case(transformation_pollen)
                  call transform_pollen(cncTrn(1:mapTransport%nSpecies,isrc), &
                                       & chemRules%rulesPollen, &
                                       & metdat, &
                                       & seconds, &
                                       & print_it)

                case (transformation_sulphur_dmat)
                  
                  zenith_cos = fu_solar_zenith_angle_cos(lon, lat, now)
                  call transform_dmat(cncTrn(1:mapTransport%nSpecies,isrc), &  ! handles also adjoint
                                    & cncSL(1:mapShortLived%nSpecies), &
                                    & chemRules%rulesSulphurDMAT, &
                                    & metdat,&
                                    & zenith_cos, &
                                    & now, &
                                    & lat, &
                                    & lon, &
                                    & seconds, &
                                    & chemRules%low_cnc_trsh, &
                                    & garb_array(isrc,:), &
                                    & print_it)

                case (transformation_acid_basic)
                  zenith_cos = fu_solar_zenith_angle_cos(lon, lat, now)
                  call transform_acid_Basic(cncTrn(1:mapTransport%nSpecies,isrc), &  ! no adjoint
                                          & cncSL(1:mapShortLived%nSpecies), &
                                          & chemRules%rulesAcidBasic, &
                                          & metdat, &
                                          & seconds, &
                                          & garb_array(isrc,:), &
                                          & zenith_cos, &
                                          & mapTransport%species, &
                                          & chemRules%low_cnc_trsh, &
                                          & .false., & ! if report (ix == 100 .and. iy == 100 .and. i3d == 1))
                                          & ix, iy, i3d, &
                                          & lat, lon, now, &
                                          & print_it)

                case (transformation_cbm4)
                  zenith_cos = fu_solar_zenith_angle_cos(lon, lat, now)
                  if (seconds > 0) then
                    !
                    ! forward and adjoint are different subroutines
                    ! Forward
                    !
                    if (have_tla) then
                      tla_point_cb4(:,i3d,ix,iy) = mapTransport%arM(:, isrc, i3d, ix, iy) / cell_volume
                      tla_point_hstart(1,i3d,ix,iy) = cb4_h_start_array(i3d,ix,iy)
                    end if

#ifdef DEBUG_MORE
call check_mass_vector(cncTrn(1:mapTransport%nSpecies,isrc), &
                     & garb_array(iSrc,:), &
                     & mapTransport%species, 'chem manager, before transform_cbm4', &
                     & chemRules%low_cnc_trsh, &
                     & mapTransport%nSpecies, ix, iy, i3d, print_it)
if (error) cycle
#endif
                    ! commented arguments to be included when chemistry regenerated...

                    call transform_cbm4(cncTrn(1:mapTransport%nSpecies,isrc), &
!                                      & photorates(:, i3d), &
                                      & chemRules%rulesCBM4, &
                                      & metdat, &
                                      & seconds, &
                                      & garb_array(isrc,:), &
                                      & zenith_cos, &
                                      & cb4_h_start_array(i3d,ix,iy), &
!                                      & now, &
                                      & print_it, & 
                                      & reactrates = reactRates(:,iSrc))
#ifdef DEBUG_MORE
call check_mass_vector(cncTrn(1:mapTransport%nSpecies,isrc), &
                     & garb_array(iSrc,:), &
                     & mapTransport%species, 'chem manager, after transform_cbm4', &
                     & chemRules%low_cnc_trsh, &
                     & mapTransport%nSpecies, ix, iy, i3d, print_it)
if (error) cycle
#endif

                  else
                    !
                    ! Adjoint
#ifdef DEBUG_MORE
                     cnctrn_tmp(1:mapTransport%nSpecies) = cncTrn(1:mapTransport%nSpecies,isrc)
                     cnctrn_tmp1(1:mapTransport%nSpecies) = tla_point_cb4(1:mapTransport%nSpecies,i3d,ix,iy)
                     cnctrn_tmp2(1:mapTransport%nSpecies) = garb_array(isrc,1:mapTransport%nSpecies)
                     fTmp = tla_point_hstart(1,i3d,ix,iy)
#endif
                    call transform_cbm4_adj(cncTrn(1:mapTransport%nSpecies,isrc), &
                                          & tla_point_cb4(:,i3d,ix,iy), &
!                                          & photorates(:, i3d), &
                                          & chemRules%rulesCBM4, &
                                          & metdat, &
                                          & seconds, &
                                          & garb_array(isrc,:), &
                                          & zenith_cos, &
                                          & tla_point_hstart(1,i3d,ix,iy), &
!                                          & now, &
                                          & print_it)
#ifdef DEBUG
                  if (.not. sum(abs(cncTrn(1:mapTransport%nSpecies,isrc))) >=0. ) then
#ifdef DEBUG_MORE
                     call msg("Before:")
                     call msg("cnctrn_tmp", cnctrn_tmp(1:mapTransport%nSpecies))
                     call msg("tla_point_cb4(:,i3d,ix,iy)", cnctrn_tmp1(1:mapTransport%nSpecies))
                     call msg(" garb_array(isrc,1:mapTransport%nSpecies)", cnctrn_tmp2(1:mapTransport%nSpecies))
                     call msg("hstart: before, after", fTmp, tla_point_hstart(1,i3d,ix,iy))
#endif
                     call msg("After:")
                     call msg(" garb_array(isrc,1:mapTransport%nSpecies)",garb_array(isrc,1:mapTransport%nSpecies))
                     call msg("tla_point_cb4(:,i3d,ix,iy)", tla_point_cb4(1:mapTransport%nSpecies,i3d,ix,iy))
                     call msg("cncTrn(1:mapTransport%nSpecies,isrc)", cncTrn(1:mapTransport%nSpecies,isrc))
                     call set_error("Gotcha cb4 adj", "Here")
                  endif
#endif
                  end if

                  
                  
                case (transformation_cbm4_SOA)
                  zenith_cos = fu_solar_zenith_angle_cos(lon, lat, now)
                  if (seconds > 0) then
                    if (have_tla) then
                      tla_point_cb4(:,i3d,ix,iy) = mapTransport%arM(:, isrc, i3d, ix, iy) / cell_volume
                      tla_point_hstart(1,i3d,ix,iy) = cb4_h_start_array(i3d,ix,iy)
                    end if

                    call transform_cbm4_SOA(cncTrn(1:mapTransport%nSpecies,isrc), &
                                              & photorates(:, i3d), &
                                              & chemRules%rulesCBM4_SOA, &
                                              & metdat, &
                                              & seconds, &
                                              & garb_array(isrc,:), &
                                              & zenith_cos, lat, &
                                              & cb4_h_start_array(i3d,ix,iy), &
                                              & now, &
                                              & print_it)
                  else
                    call transform_cbm4_SOA_adj(cncTrn(1:mapTransport%nSpecies,isrc), &
                                                  & tla_point_cb4(:,i3d,ix,iy), &
                                                  & photorates(:, i3d), &
                                                  & chemRules%rulesCBM4_SOA, &
                                                  & metdat, &
                                                  & seconds, &
                                                  & garb_array(isrc,:), &
                                                  & zenith_cos, lat, &
                                                  & tla_point_hstart(1,i3d,ix,iy), &
                                                  & now, &
                                                  & print_it)

                  end if

                  
                case (transformation_cbm42_strato)
                  if (seconds > 0) then
                    !
                    ! forward and adjoint are different subroutines
                    ! Forward
                    !
                    if (have_tla) then
                      tla_point_cb4(:,i3d,ix,iy) = mapTransport%arM(:, isrc, i3d, ix, iy) / cell_volume
                      tla_point_hstart(1,i3d,ix,iy) = cb4_h_start_array(i3d,ix,iy)
                    end if

                    call transform_cbm42_strato(cncTrn(1:mapTransport%nSpecies,isrc), &
                                              & photorates(:, i3d), &
                                              & chemRules%rulesCBM42_strato, &
                                              & metdat, &
                                              & seconds, &
                                              & garb_array(isrc,:), &
                                              & zenith_cos, lat, &
                                              & cb4_h_start_array(i3d,ix,iy), &
                                              & now, &
                                              & print_it)
#ifdef DEBUG_MORE
call check_mass_vector(cncTrn(1:mapTransport%nSpecies,isrc), &
                     & garb_array(iSrc,:), &
                     & mapTransport%species, 'chem manager, after transform_cbm42_strato', &
                     & chemRules%low_cnc_trsh, &
                     & mapTransport%nSpecies, ix, iy, i3d, print_it)
#endif

                  else
                    !
                    ! adjoint
                    !
                    call transform_cbm42_strato_adj(cncTrn(1:mapTransport%nSpecies,isrc), &
                                                  & tla_point_cb4(:,i3d,ix,iy), &
                                                  & photorates(:, i3d), &
                                                  & chemRules%rulesCBM42_strato, &
                                                  & metdat, &
                                                  & seconds, &
                                                  & garb_array(isrc,:), &
                                                  & zenith_cos, lat, &
                                                  & tla_point_hstart(1,i3d,ix,iy), &
                                                  & now, &
                                                  & print_it)

                  end if

                case (transformation_cbm42_strato_SOA)
                  if (seconds > 0) then
                    !
                    ! forward and adjoint are different subroutines
                    ! Forward
                    !
                    if (have_tla) then
                      tla_point_cb4(:,i3d,ix,iy) = mapTransport%arM(:, isrc, i3d, ix, iy) / cell_volume
                      tla_point_hstart(1,i3d,ix,iy) = cb4_h_start_array(i3d,ix,iy)
                    end if

                    call transform_cbm42_strato_SOA(cncTrn(1:mapTransport%nSpecies,isrc), &
                                              & photorates(:, i3d), &
                                              & chemRules%rulesCBM42_strato_SOA, &
                                              & metdat, &
                                              & seconds, &
                                              & garb_array(isrc,:), &
                                              & zenith_cos, lat, &
                                              & cb4_h_start_array(i3d,ix,iy), &
                                              & now, &
                                              & print_it)
#ifdef FORCING_TOP
if(i3d == mapTransport%n3d)then
if(ix == 1 .and. iy == 1) call msg_warning('Forcing top O3 to ~0.1ppm', sub_name)
! Concentration at 62km is around 2e-10 mole/m3, roughly corresponds to 0.1ppm.
! Note that daytime it is essentially zero, night-time is fine. But it looks like MOZART also has problems with this

if(indOzoneToForce /= int_missing) cncTrn(indOzoneToForce,isrc) = 2e-10 + (1.-zenith_cos) * 3e-10

endif
#endif




#ifdef DEBUG_MORE
call check_mass_vector(cncTrn(1:mapTransport%nSpecies,isrc), &
                     & garb_array(iSrc,:), &
                     & mapTransport%species, 'chem manager, after transform_cbm42_strato_SOA', &
                     & chemRules%low_cnc_trsh, &
                     & mapTransport%nSpecies, ix, iy, i3d, print_it)
#endif

                  else
                    !
                    ! adjoint
                    !
                    call transform_cbm42_strato_SOA_adj(cncTrn(1:mapTransport%nSpecies,isrc), &
                                                  & tla_point_cb4(:,i3d,ix,iy), &
                                                  & photorates(:, i3d), &
                                                  & chemRules%rulesCBM42_strato_SOA, &
                                                  & metdat, &
                                                  & seconds, &
                                                  & garb_array(isrc,:), &
                                                  & zenith_cos, lat, &
                                                  & tla_point_hstart(1,i3d,ix,iy), &
                                                  & now, &
                                                  & print_it)

                  end if

                case (transformation_cbm5_SOA)
                  zenith_cos = fu_solar_zenith_angle_cos(lon, lat, now)
                  if (seconds > 0) then
                    if (have_tla) then
                      tla_point_cb4(:,i3d,ix,iy) = mapTransport%arM(:, isrc, i3d, ix, iy) / cell_volume
                      tla_point_hstart(1,i3d,ix,iy) = cb4_h_start_array(i3d,ix,iy)
                    end if

                    call transform_cbm5_SOA(cncTrn(1:mapTransport%nSpecies,isrc), &
                                              & photorates(:, i3d), &
                                              & chemRules%rulesCBM5_SOA, &
                                              & metdat, &
                                              & seconds, &
                                              & garb_array(isrc,:), &
                                              & zenith_cos, lat, &
                                              & cb4_h_start_array(i3d,ix,iy), &
                                              & now, &
                                              & print_it)
                  else
                    call transform_cbm5_SOA_adj(cncTrn(1:mapTransport%nSpecies,isrc), &
                                                  & tla_point_cb4(:,i3d,ix,iy), &
                                                  & photorates(:, i3d), &
                                                  & chemRules%rulesCBM5_SOA, &
                                                  & metdat, &
                                                  & seconds, &
                                                  & garb_array(isrc,:), &
                                                  & zenith_cos, lat, &
                                                  & tla_point_hstart(1,i3d,ix,iy), &
                                                  & now, &
                                                  & print_it)

                  end if

                  
                case (transformation_cbm5_strato_SOA)
                  zenith_cos = fu_solar_zenith_angle_cos(lon, lat, now)
                  if (seconds > 0) then
                    if (have_tla) then
                      tla_point_cb4(:,i3d,ix,iy) = mapTransport%arM(:, isrc, i3d, ix, iy) / cell_volume
                      tla_point_hstart(1,i3d,ix,iy) = cb4_h_start_array(i3d,ix,iy)
                    end if

                    call transform_cbm5_strato_SOA(cncTrn(1:mapTransport%nSpecies,isrc), &
                                              & photorates(:, i3d), &
                                              & chemRules%rulesCBM5_strato_SOA, &
                                              & metdat, &
                                              & seconds, &
                                              & garb_array(isrc,:), &
                                              & zenith_cos, lat, &
                                              & cb4_h_start_array(i3d,ix,iy), &
                                              & now, &
                                              & print_it)
                  else
                    call transform_cbm5_strato_SOA_adj(cncTrn(1:mapTransport%nSpecies,isrc), &
                                                  & tla_point_cb4(:,i3d,ix,iy), &
                                                  & photorates(:, i3d), &
                                                  & chemRules%rulesCBM5_strato_SOA, &
                                                  & metdat, &
                                                  & seconds, &
                                                  & garb_array(isrc,:), &
                                                  & zenith_cos, lat, &
                                                  & tla_point_hstart(1,i3d,ix,iy), &
                                                  & now, &
                                                  & print_it)

                  end if

                  
                case (transformation_radioactive)
                  call transform_radioactive(cncTrn(1:mapTransport%nSpecies,isrc), &   ! no adjoint
!                                            & mapTransport%arM(:,isrc,i3d,ix,iy), &
                                           & chemRules%rulesRadioactive, &
                                           & metdat,&
                                           & seconds, &
                                           & print_it)
                  if(i3d == 1)then
                    !
                    ! Radioactive decay is defined also for the deposited species but
                    ! the actual decay must be computed only if cumulative output is requested
                    ! For rates, we need just that - rates, i.e. no decay is allowed.
                    ! Since species-wise definition is technically allowed, have to be careful here
                    !
                    if(fu_true(ifDryDep_cumulative_in_output))then
                    call transform_radioactive(mapDryDep%arM(:,isrc,1,ix,iy), &
                                             & chemRules%rulesRadioactive, &
                                             & metdat,&
                                             & seconds, &
                                             & print_it)
                    elseif(fu_false(ifDryDep_cumulative_in_output))then
                    else
                      call set_error('Dry deposition averaging type silja_undefined is not allowed in radioactive runs', &
                                   & sub_name)
                      cycle
                    endif
                    if(fu_true(ifWetDep_cumulative_in_output))then
                    call transform_radioactive(mapWetDep%arM(:,isrc,1,ix,iy), &
                                             & chemRules%rulesRadioactive, &
                                             & metdat,&
                                             & seconds, &
                                             & print_it)
                    elseif(fu_false(ifWetDep_cumulative_in_output))then
                    else
                      call set_error('Wet deposition averaging type silja_undefined is not allowed in radioactive runs', &
                                   & sub_name)
                      cycle
                    endif
                  endif  ! i3d ==1
                case(int_missing, transformation_inert_PM)
                  continue
                case default
                  call msg('Unknown transformation:',chemrules%iTransformTypes(itransf))
                  call set_error('Unknown transformation', sub_name)
                  cycle
              end select
              if (error) cycle
            end do      ! iTransf
            
            ! Aerosol dynamics, forward ===============================================================
            
            if (seconds > 0 .and. chemRules%nAerosolDynamics > 0) then
             

#ifdef DEBUG_MORE
call check_mass_vector(cncTrn(1:mapTransport%nSpecies,isrc), &
                     & garb_array(iSrc,:), &
                     & mapTransport%species, 'chem manager, before AD', chemRules%low_cnc_trsh, &
                     & mapTransport%nSpecies, ix, iy, i3d, print_it)
#endif

              do iTransf = 1, chemRules%nAerosolDynamics
                ! skip lagrangian-only transform
                if(.not. fu_if_eulerian_present(chemrules%iWhomToApplyAerDyn(iTransf)))cycle 
#ifdef DUMP
call msg('before the aerosol dynamics call:' + fu_str(chemrules%iAerosolDynTypes(itransf)))
do iSpecies = 1, mapTransport%nSpecies
  call msg('Mass before the next call:' + &
         & fu_str(mapTransport%species(iSpecies)), &
         & cncTrn(iSpecies,isrc))
end do
#endif
                select case(chemrules%iAerosolDynTypes(iTransf))
                case(aerosol_dynamics_basic)
                  cncAer(1:mapAerosol%nSpecies) = &
                       & mapAerosol%arM(1:mapAerosol%nSpecies,isrc,i3d,ix,iy) / cell_volume
                  call transform_AerDynBasic(chemRules%rulesAerDynBasic, &
                                           & cncTrn(1:mapTransport%nSpecies,isrc), &
                                           & cncSL(1:mapShortLived%nSpecies), &
                                           & cncAer(1:mapAerosol%nSpecies), &
!                                           & mapTransport%arM(:,isrc,i3d,ix,iy), &
!                                           & mapShortlived%arM(:,isrc,i3d,ix,iy), &
!                                           & mapAerosol%arM(:,isrc,i3d,ix,iy), &
                                           & chemRules%low_cnc_trsh, &
                                           & garb_array(isrc,:), &
                                           & chemRules%ChemRunSetup%mapVolume2NumberCnc, &
                                           & metdat, &
                                           & seconds, print_it)
                  mapAerosol%arM(1:mapAerosol%nSpecies,isrc,i3d,ix,iy) = &
                                                    & cncAer(1:mapAerosol%nSpecies) * cell_volume
!                       & mapAerosol%arM(1:mapAerosol%nSpecies,isrc,i3d,ix,iy) * cell_volume

                case(aerosol_dynamics_simple)
                  zenith_cos = fu_solar_zenith_angle_cos(lon, lat, now)
                  call transform_AerDynSimple(cncTrn(1:mapTransport%nSpecies,isrc), &
                                            & cncSL(1:mapShortLived%nSpecies), &
!                                            & mapTransport%arM(:,isrc,i3d,ix,iy), &
!                                            & mapShortLived%arM(:,isrc,i3d,ix,iy), &
                                            & garb_array(isrc,:), &
                                            & chemRules%rulesAerDynSimple, chemRules%low_cnc_trsh, &
                                            & metdat, seconds, zenith_cos, print_it)

                case(aerosol_dynamics_Mid_Atmosph)
!$OMP CRITICAL(aerdyn)  !This is only to order the error etc messages. Can be removed if one this slows the code too much!
                  call transform_AerDynMidAtm(cncTrn(1:mapTransport%nSpecies,isrc), &
                                            & cncSL(1:mapShortLived%nSpecies), &
                                            & cncAer(1:mapAerosol%nSpecies), &
                                            & garb_array(isrc,:), &
                                            & chemRules%rulesAerDynMidAtmosph, &
                                            & chemRules%low_cnc_trsh, &
                                            & metdat, &
                                            & seconds, &
                                            & mapTransport%nSpecies, &
                                            & mapTransport%Species, ix, iy, i3d, &
                                            & print_it)
!$OMP END CRITICAL(aerdyn)
                                
                case(aerosol_dynamics_VBS)
                  call transform_AerDynVBS(cncTrn(1:mapTransport%nSpecies,isrc), &
                                         & cncSL(1:mapShortLived%nSpecies), &
                                         & garb_array(isrc,:), &
                                         & chemRules%rulesAerDynVBS, chemRules%low_cnc_trsh, &
                                         & metdat, seconds, print_it)

                case(int_missing)
                  continue

                case default
                  call msg('Unknown aerosol dynamics:',chemrules%iAerosolDynTypes(itransf))
                  call set_error('Unknown aerosol dynamics', sub_name)
                  cycle
                end select
              enddo    ! nAerosolDynamics
              call check_mass_vector(cncTrn(1:mapTransport%nSpecies,isrc), &
                     & garb_array(iSrc,:), &
                     & mapTransport%species, 'chem manager, after AD', chemRules%low_cnc_trsh, &
                     & mapTransport%nSpecies, ix, iy, i3d, print_it)
            end if   ! seconds > 0 , forward aerosol dynamics

            ! Any trouble?
            !
            if(error .or. print_it)then
               !$OMP  CRITICAL (barkMassVector)
 
              if(error)call set_error('Trouble with chemistry start report', sub_name)
              call msg('ix,iy,i3d,iSrc:' + fu_str(ix) + ',' + fu_str(iy),i3d,iSrc)
              call msg('longitude =',fu_lon_geographical_from_grid(real(ix),real(iy),dispersion_grid))
              call msg('latitude =', fu_lat_geographical_from_grid(real(ix),real(iy),dispersion_grid))
              call msg('start_time =',(/fu_year(now), fu_mon(now), fu_day(now), fu_hour(now), &
                                      & fu_min(now), nint(fu_sec(now))/))
              call msg('end_time =',(/fu_year(now+fu_set_interval_sec(seconds)), &
                                    & fu_mon(now+fu_set_interval_sec(seconds)), &
                                    & fu_day(now+fu_set_interval_sec(seconds)), &
                                    & fu_hour(now+fu_set_interval_sec(seconds)), &
                                    & fu_min(now+fu_set_interval_sec(seconds)), &
                                    & nint(fu_sec(now+fu_set_interval_sec(seconds)))/))
              do iSpecies = 1, meteo_input%nQuantities
                write(unit=chTmp,fmt=*)fu_quantity_short_string(meteo_input%quantity(iSpecies)), &
                                     & ' = ', metdat(iSpecies)
                call msg(chTmp)
              enddo
              call msg("photorates(:, i3d)", photorates(:, i3d))
              call msg('Dump, x,y,z,Src=('+fu_str(ix)+','+fu_str(iy)+','+fu_str(i3d), isrc)
              call msg('Mass before & after:')
              do iSpecies = 1, mapTransport%nSpecies
                write(unit=chTmp,fmt=*)trim(fu_str(mapTransport%species(iSpecies))), ' = ' , &
                                     & mapTransport%arM(iSpecies,iSrc,i3d,ix,iy) / cell_volume, &
                                     & ' mole   ! ', &
                                     & cncTrn(iSpecies,isrc)
                call msg(chTmp)
              end do
              do iSpecies = 1, mapShortLived%nSpecies
                write(unit=chTmp,fmt=*)trim(fu_str(mapShortLived%species(iSpecies))), ' = ' , &
                                     & mapShortLived%arM(iSpecies,iSrc,i3d,ix,iy) / cell_volume, &
                                     & ' mole   ! ', &
                                     & cncSL(iSpecies)
                call msg(chTmp)
              end do
              do iSpecies = 1, mapAerosol%nSpecies
                write(unit=chTmp,fmt=*)trim(fu_str(mapAerosol%species(iSpecies))), ' = ' , &
                                     & mapAerosol%arM(iSpecies,iSrc,i3d,ix,iy) / cell_volume, &
                                     & ' mole   ! ', &
                                     & cncAer(iSpecies)
                call msg(chTmp)
              end do
              call msg("")
              call msg("")
              call msg("Trouble with chemistry end_report")
              print_it = .false.
              !$OMP END CRITICAL (barkMassVector)
              if(error)call set_error('Trouble with chemistry end_report', sub_name)
            endif

#ifdef DUMP
call msg('After all calls')
do iSpecies = 1, mapTransport%nSpecies
  call msg('Mass after all calls:' + &
         & fu_str(mapTransport%species(iSpecies)), &
         & cncTrn(iSpecies,isrc))
end do
#endif
          end do     ! iSrc
          
          if (defined(mapReactRates)) then
             mapReactRates%arM(:,1:mapReactRates%nSrc,i3d,ix,iy) = reactRates(:,1:mapReactRates%nSrc)
          endif

          mapTransport%arM(1:mapTransport%nSpecies,1:mapTransport%nSrc,i3d,ix,iy) = &
               & cncTrn(1:mapTransport%nSpecies,1:mapTransport%nSrc) * &
!               & mapTransport%arM(1:mapTransport%nSpecies,1:mapTransport%nSrc,i3d,ix,iy) * &
               & cell_volume
         if (mapShortLived%nSpecies > 0) then
          mapShortLived%arM(1:mapShortLived%nSpecies,1,i3d,ix,iy) = cncSL(1:mapShortLived%nSpecies) * &
                                                                  & cell_volume
         endif                                                          
          garb_array(1:mapTransport%nSrc, 1:mapTransport%nSpecies) = &
                                & garb_array(1:mapTransport%nSrc, 1:mapTransport%nSpecies) * cell_volume
        end do     !n3d
        
        if (seconds > 0) then
           call scavenge_column(mapTransport, mapWetDep, garb_array, tla_step, seconds, &
                & chemRules%rulesDeposition, chemRules%low_mass_trsh, metdat_col, ix, iy, pwc_col, &
                & pHorizInterpStruct, pVertInterpStruct, ifHorizInterp, ifVertInterp, met_buf)
           if(error)call set_error('Trouble with scavenging', sub_name)
        end if
      end do    ! nx
    end do   ! ny
    !$OMP END DO

    !$OMP CRITICAL(transf_maps_work)
    ! Collect the garbage
    garbage(1:mapTransport%nSrc, 1:mapTransport%nSpecies) &
         & = garbage(:,:) + garb_array(1:mapTransport%nSrc,1:mapTransport%nSpecies)
    !$OMP END CRITICAL(transf_maps_work)

    !$OMP END PARALLEL
!call msg("Transform_maps done")


 
  end subroutine transform_maps
  

  !**************************************************************************************

  subroutine transform_lagrangian_part(lpSet, mapDryDep, mapWetDep, &
                                     & garbage, tla_step, &
                                     & met_buf, disp_buf, meteo_input, &
                                     & pHorizInterpStruct, pVertInterpStruct, &
                                     & ifHorizInterp, ifVertInterp,&
                                     & chemRules, seconds, now)
    !
    ! Main chemistry and aerosol transformation routine for Lagrnagian environment
    !
    implicit none
    
    ! Imported parameters
    type(Tlagrange_particles_set), intent(inout) :: lpSet
    type(Tmass_map), intent(inout) :: mapDryDep, mapWetDep
    real, dimension(:,:), intent(inout) :: garbage
    type(t_tla_step), intent(inout) :: tla_step
    type(Tfield_buffer), pointer :: met_buf, disp_buf
    type(Tmeteo_input), pointer :: meteo_input
    type(THorizInterpStruct), pointer :: pHorizInterpStruct
    type(TVertInterpStruct), pointer :: pVertInterpStruct
    logical, intent(in) :: ifHorizInterp, ifVertInterp
    type(Tchem_rules), intent(inout) :: chemRules
    type(silja_time), intent(in) :: now
    real, intent(in) :: seconds

    ! Local variables
    integer :: iTransf, iP, iSrc, iSpecies
    real :: cell_volume, lon, lat, zenith_cos
    real, dimension(:), pointer :: metdat, cncTrn, cncAer, cncSL
    real, dimension(:,:), pointer :: garb_array
    logical :: print_it
    character(len=100) :: chTmp

    if (chemRules%nTransformations + chemRules%nAerosolDynamics == 0) return

    !
    ! The main cycle over the particles
    !
    ! Garbage comes as an array. For parallelisation, we temporarily use 
    ! a real array which will be private for everyone.
    !

    !$OMP PARALLEL if (if_OMP_chemistry) DEFAULT(SHARED) PRIVATE(metdat, iP, itransf, isrc, iSpecies, &
    !$OMP & cell_volume, zenith_cos, lat, lon, garb_array, print_it, cncTrn, cncAer, cncSL)
    metdat => fu_work_array(meteo_input%nQuantities)
    garb_array => fu_work_array_2D(lpSet%nSrcs,lpSet%nSpeciesTrn)
    cncTrn => fu_work_array(max_species)
    cncAer => fu_work_array(max_species)
    cncSL => fu_work_array(max_species)
    garb_array(1:lpSet%nSrcs, 1:lpSet%nSpeciesTrn) = 0.0

    !$OMP DO
    do iP = 1, lpSet%nop
      if (error) cycle
      if(lpSet%lpStatus(iP) == int_missing)cycle ! skip void particles
      iSrc = mod(lpSet%lpStatus(iP), 100)
      
      call fill_in_meteo_input(meteo_input, metdat, &
                             & met_buf, disp_buf, &
                             & nint(lpSet%lpDyn(lp_x,iP)), nint(lpSet%lpDyn(lp_y,iP)), &
                             & nint(lpSet%lpDyn(lp_z,iP)), &
                             & pHorizInterpStruct, pVertInterpStruct, &  ! meteo to disp
!                            & ifHorizInterp, ifVertInterp)
                             & .false., .false.) ! LP fly in meteo!!!
      cell_volume = lpSet%lpDyn(lp_dx,iP) * lpSet%lpDyn(lp_dy,iP) * lpSet%lpDyn(lp_dz,iP)

      cncTrn(1:lpSet%nSpeciesTrn) = lpSet%lpMassTrn(1:lpSet%nSpeciesTrn,iP) / cell_volume
      cncSL(1:lpSet%nSpeciesSL) = 0.
!      lpSet%lpMassTrn(1:lpSet%nSpeciesTrn,iP) = lpSet%lpMassTrn(:,iP) / cell_volume
      garb_array(1:lpSet%nSrcs,1:lpSet%nSpeciesTrn) = garb_array(1:lpSet%nSrcs,1:lpSet%nSpeciesTrn) / cell_volume

      ! Aerosol dynamics, adjoint ===============================================================
            
      if (seconds < 0) then

        do iTransf = 1, chemRules%nAerosolDynamics
          if(.not. fu_if_lagrangian_present(chemrules%iWhomToApplyAerDyn(iTransf)))cycle ! skip eulerian-only AD
          select case(chemrules%iAerosolDynTypes(iTransf))
            case(aerosol_dynamics_basic)
              cncAer(1:lpSet%nSpeciesAer) = lpSet%lpMassAer(:,iP) / cell_volume
!              lpSet%lpMassAer(1:lpSet%nSpeciesAer,iP) = lpSet%lpMassAer(:,iP) / cell_volume
              call transform_AerDynBasic(chemRules%rulesAerDynBasic, &
                                       & cncTrn(1:lpSet%nSpeciesTrn), &
                                        & cncSL(1:lpSet%nSpeciesSL), &
                                        & cncAer(1:lpSet%nSpeciesAer), &
!                                        &lpSet%lpMassTrn(:,iP), &
!                                       & lpSet%lpMassSL(:,iP), &
!                                       & lpSet%lpMassAer(:,iP), &
                                       & chemRules%low_cnc_trsh, &
                                       & garb_array(iSrc,:), &
                                       & chemRules%ChemRunSetup%mapVolume2NumberCnc, &
                                       & metdat, &
                                       & seconds, print_it)
              lpSet%lpMassAer(1:lpSet%nSpeciesAer,iP) = cncAer(:) * cell_volume
!              lpSet%lpMassAer(1:lpSet%nSpeciesAer,iP) = lpSet%lpMassAer(:,iP) * cell_volume

            case(aerosol_dynamics_simple)
              zenith_cos = fu_solar_zenith_angle_cos(lon, lat, now)
              call transform_AerDynSimple(cncTrn(1:lpSet%nSpeciesTrn), &
                                        & cncSL(1:lpSet%nSpeciesSL), &
!                                        &lpSet%lpMassTrn(:,iP), &
!                                        & lpSet%lpMassSL(:,iP), &
                                        & garb_array(iSrc,:), &
                                        & chemRules%rulesAerDynSimple, chemRules%low_cnc_trsh, &
                                        & metdat, seconds,zenith_cos, print_it)
            case(aerosol_dynamics_Mid_Atmosph)
              call transform_AerDynMidAtm(cncTrn(1:lpSet%nSpeciesTrn), &
                                        & cncSL(1:lpSet%nSpeciesSL), &
                                        & cncAer(1:lpSet%nSpeciesAer), &
                                        & garb_array(isrc,:), &
                                        & chemRules%rulesAerDynMidAtmosph, &
                                        & chemRules%low_cnc_trsh, &
                                        & metdat, &
                                        & seconds, &
                                        & lpSet%nSpeciesTrn, &
                                        & lpSet%spTransp, iP, int_missing, int_missing, &
                                        & print_it)
            case(aerosol_dynamics_VBS)
              call transform_AerDynVBS(cncTrn(1:lpSet%nSpeciesTrn), &
                                         & cncSL(1:lpSet%nSpeciesSL), &
                                         & garb_array(isrc,:), &
                                         & chemRules%rulesAerDynVBS, chemRules%low_cnc_trsh, &
                                         & metdat, seconds, print_it)
            case(int_missing)

              continue

            case default
              call msg('Unknown aerosol dynamics:',chemrules%iAerosolDynTypes(itransf))
              call set_error('Unknown aerosol dynamics','transform_lagrangian_part')
              cycle
          end select  ! AerDyn type
        enddo    ! nAerosolDynamics
      end if  ! seconds < 0

      if (debug_level > 0 .and. seconds > 0) then
        call check_mass_vector(cncTrn(1:lpSet%nSpeciesTrn), &
!                                        & lpSet%lpMassTrn(:,iP), &
                             & garb_array(iSrc,:), &
                             & lpSet%spTransp, 'chem manager, before Lagr chemistry', &
                             & chemRules%low_cnc_trsh, &
                             & lpSet%nSpeciesTrn, iP, int_missing, int_missing, print_it)
      end if
            
      do iTransf = 1, chemRules%nTransformations
        if(.not. fu_if_eulerian_present(chemrules%iWhomToApplyTransform(iTransf)))cycle ! skip lagrangian-only transform
        select case(chemrules%iTransformTypes(itransf))
          case(transformation_passive)
            call transform_passive(cncTrn(1:lpSet%nSpeciesTrn), &
!                                        & lpSet%lpMassTrn(:,iP), &
                                 & chemRules%rulesPassive, &
                                 & metdat, nint(lpSet%lpDyn(lp_z,iP)), &
                                 & seconds,  print_it)

            case(transformation_pollen)
              call transform_pollen(cncTrn(1:lpSet%nSpeciesTrn), &
!                                        & lpSet%lpMassTrn(:,iP), &
                                   & chemRules%rulesPollen, &
                                   & metdat, &
                                   & seconds, print_it)

            case (transformation_sulphur_dmat)
              lat = fu_lat_geographical_from_grid(lpSet%lpDyn(lp_x,iP), lpSet%lpDyn(lp_y,iP), &
                                                & dispersion_grid)
              lon = fu_lon_geographical_from_grid(lpSet%lpDyn(lp_x,iP), lpSet%lpDyn(lp_y,iP), &
                                                & dispersion_grid)
              zenith_cos = fu_solar_zenith_angle_cos(lon, lat, now)
              call transform_dmat(cncTrn(1:lpSet%nSpeciesTrn), &
                                & cncSL(1:lpSet%nSpeciesSL), &
!                                        &lpSet%lpMassTrn(:,iP), &
!                                & lpSet%lpMassSL(:,iP), &
                                & chemRules%rulesSulphurDMAT, &
                                & metdat,&
                                & zenith_cos, &
                                & now, &
                                & lat, &
                                & lon, &
                                & seconds, &
                                & chemRules%low_cnc_trsh, &
                                & garb_array(isrc,:), &
                                & print_it)

            case (transformation_acid_basic)
              lat = fu_lat_geographical_from_grid(lpSet%lpDyn(lp_x,iP), lpSet%lpDyn(lp_y,iP), &
                                                & dispersion_grid)
              lon = fu_lon_geographical_from_grid(lpSet%lpDyn(lp_x,iP), lpSet%lpDyn(lp_y,iP), &
                                                & dispersion_grid)
              zenith_cos = fu_solar_zenith_angle_cos(lon, lat, now)

              call transform_acid_Basic(cncTrn(1:lpSet%nSpeciesTrn), &
                                      & cncSL(1:lpSet%nSpeciesSL), &
!                                        & lpSet%lpMassTrn(:,iP), &
!                                      & lpSet%lpMassSL(:,iP), &
                                      & chemRules%rulesAcidBasic, &
                                      & metdat, &
                                      & seconds, &
                                      & garb_array(isrc,:), &
                                      & zenith_cos, &
                                      & lpSet%spTransp, &
                                      & chemRules%low_cnc_trsh, &
                                      & .false., & ! if report (ix == 100 .and. iy == 100 .and. i3d == 1))
                                      & iP, int_missing, int_missing, &
                                      & lat, lon, now, &
                                      & print_it)

            case (transformation_cbm4)
              call set_error('CB4 is not applicable in Lagrangian environment','')
                  
            case (transformation_radioactive)
              call transform_radioactive(cncTrn(1:lpSet%nSpeciesTrn), &
!                                        & lpSet%lpMassTrn(:,iP), &
                                       & chemRules%rulesRadioactive, &
                                       & metdat,&
                                       & seconds, print_it)
            case(int_missing, transformation_inert_PM)
              continue
            case default
              call msg('Unknown transformation:',chemrules%iTransformTypes(itransf))
              call set_error('Unknown transformation','transform_lagrangian_part')
              cycle
        end select
        if (error) cycle
      end do      ! iTransf
            
      if (debug_level > 0 .and. seconds > 0) then
        call check_mass_vector(cncTrn(1:lpSet%nSpeciesTrn), &
!                                        & lpSet%lpMassTrn(:,iP), &
                             & garb_array(iSrc,:), &
                             & lpSet%spTransp, 'chem manager, before Lagr AD', chemRules%low_cnc_trsh, &
                             & lpSet%nSpeciesTrn, iP, int_missing, int_missing, print_it)
      end if

      ! Aerosol dynamics, forward ===============================================================

      if (seconds > 0) then
        do iTransf = 1, chemRules%nAerosolDynamics
          if(.not. fu_if_lagrangian_present(chemrules%iWhomToApplyAerDyn(iTransf)))cycle ! skip eulerian-only transform
          select case(chemrules%iAerosolDynTypes(iTransf))
            case(aerosol_dynamics_basic)
              cncAer(1:lpSet%nSpeciesAer) = lpSet%lpMassAer(:,iP) /cell_volume
!              lpSet%lpMassAer(1:lpSet%nSpeciesAer,iP) = lpSet%lpMassAer(:,iP) /cell_volume
              call transform_AerDynBasic(chemRules%rulesAerDynBasic, &
                                       & cncTrn(1:lpSet%nSpeciesTrn), &
                                       & cncSL(1:lpSet%nSpeciesSL), &
                                       & cncAer(1:lpSet%nSpeciesAer), &
!                                        &lpSet%lpMassTrn(:,iP), &
!                                       & lpSet%lpMassSL(:,iP), &
!                                       & lpSet%lpMassAer(:,iP), &
                                       & chemRules%low_cnc_trsh, &
                                       & garb_array(isrc,:), &
                                       & chemRules%ChemRunSetup%mapVolume2NumberCnc, &
                                       & metdat, &
                                       & seconds, print_it)
              lpSet%lpMassAer(1:lpSet%nSpeciesAer,iP) = cncAer(:) * cell_volume
!              lpSet%lpMassAer(1:lpSet%nSpeciesAer,iP) = lpSet%lpMassAer(:,iP) * cell_volume

            case(aerosol_dynamics_simple)
              zenith_cos = fu_solar_zenith_angle_cos(lon, lat, now)
              call transform_AerDynSimple(cncTrn(1:lpSet%nSpeciesTrn), &
                                        & cncSL(1:lpSet%nSpeciesSL), &
!                                        & lpSet%lpMassTrn(:,iP), &
!                                        & lpSet%lpMassSL(:,iP), &
                                        & garb_array(isrc,:), &
                                        & chemRules%rulesAerDynSimple, chemRules%low_cnc_trsh, &
                                        & metdat, seconds, zenith_cos, print_it)

            case(aerosol_dynamics_Mid_Atmosph)
              call transform_AerDynMidAtm(cncTrn(1:lpSet%nSpeciesTrn), &
                                        & cncSL(1:lpSet%nSpeciesSL), &
                                        & cncAer(1:lpSet%nSpeciesAer), &
                                        & garb_array(isrc,:), &
                                        & chemRules%rulesAerDynMidAtmosph, &
                                        & chemRules%low_cnc_trsh, &
                                        & metdat, &
                                        & seconds, &
                                        & lpSet%nSpeciesTrn, lpSet%spTransp, &
                                        & iP, int_missing, int_missing, &
                                        & print_it)
                          
            case(aerosol_dynamics_VBS)
              call transform_AerDynVBS(cncTrn(1:lpSet%nSpeciesTrn), &
                                         & cncSL(1:lpSet%nSpeciesTrn), &
                                         & garb_array(isrc,:), &
                                         & chemRules%rulesAerDynVBS, chemRules%low_cnc_trsh, &
                                         & metdat, seconds, print_it)

            case(int_missing)
              continue

            case default
              call msg('Unknown aerosol dynamics:',chemrules%iAerosolDynTypes(itransf))
              call set_error('Unknown aerosol dynamics','transform_lagrangian_part')
              cycle
          end select
        enddo    ! nAerosolDynamics
      end if  ! seconds > 0

      if (debug_level > 0 .and. seconds > 0) then
        call check_mass_vector(cncTrn(1:lpSet%nSpeciesTrn), &
!                             & lpSet%lpMassTrn(:,iP), &
                             & garb_array(iSrc,:), &
                             & lpSet%spTransp, 'chem manager, after Lagr AD', chemRules%low_cnc_trsh, &
                             & lpSet%nSpeciesTrn, iP, int_missing, int_missing, print_it)
      end if

      if(error)then
        !$OMP CRITICAL(barkMassVector)
        call set_error('Trouble with chemistry','transform_lagrangian_part')
        call msg('iParticle,iSrc:',iP,iSrc)
!      endif
!      if(print_it)then
        call msg('Dump, x,y,z =', (/lpSet%lpDyn(lp_x,iP), lpSet%lpDyn(lp_y,iP), lpSet%lpDyn(lp_z,iP)/))
        call msg('longitude =',fu_lon_geographical_from_grid(lpSet%lpDyn(lp_x,iP), &
                                                           & lpSet%lpDyn(lp_y,iP), dispersion_grid))
        call msg('latitude =', fu_lat_geographical_from_grid(lpSet%lpDyn(lp_x,iP), &
                                                           & lpSet%lpDyn(lp_y,iP),dispersion_grid))
        call msg('start_time =',(/fu_year(now), fu_mon(now), fu_day(now), fu_hour(now), &
                                & fu_min(now), nint(fu_sec(now))/))
        call msg('end_time =',(/fu_year(now+fu_set_interval_sec(seconds)), &
                              & fu_mon(now+fu_set_interval_sec(seconds)), &
                              & fu_day(now+fu_set_interval_sec(seconds)), &
                              & fu_hour(now+fu_set_interval_sec(seconds)), &
                              & fu_min(now+fu_set_interval_sec(seconds)), &
                              & nint(fu_sec(now+fu_set_interval_sec(seconds)))/))
        do iSpecies = 1, meteo_input%nQuantities
          write(unit=chTmp,fmt=*)fu_quantity_short_string(meteo_input%quantity(iSpecies)), &
                                     & ' = ', metdat(iSpecies)
          call msg(chTmp)
        enddo
        call msg('Mass before & after:')
        do iSpecies = 1, lpSet%nSpeciesTrn
          write(unit=chTmp,fmt=*)fu_str(lpSet%spTransp(iSpecies)), ' = ' , &
                       & lpSet%lpMassTrn(iSpecies,iP) / cell_volume, &
                       & ' mole   ! ', &
                       & cncTrn(iSpecies)
                call msg(chTmp)
!        call msg('Mass before & after:' + fu_str(lpSet%spTransp(iSpecies)), &
!                                        & lpSet%lpMassTrn(iSpecies,iP) / cell_volume, &
!                                        & cncTrn(iSpecies))
        end do
        print_it = .false.
        !$OMP END CRITICAL(barkMassVector)
      endif

      lpSet%lpMassTrn(1:lpSet%nSpeciesTrn,iP) = cncTrn(1:lpSet%nSpeciesTrn) * cell_volume
!      lpSet%lpMassTrn(1:lpSet%nSpeciesTrn,iP) = lpSet%lpMassTrn(:,iP) * cell_volume
      garb_array(1:lpSet%nSrcs, 1:lpSet%nSpeciesTrn) = garb_array(1:lpSet%nSrcs, 1:lpSet%nSpeciesTrn) * cell_volume
    end do   ! iP
    !$OMP END DO

    call free_work_array(metdat)
    ! Collect the garbage
    !
    !$OMP CRITICAL(transf_maps_work)
    garbage(1:lpSet%nSrcs, 1:lpSet%nSpeciesTrn) = garbage(:,:) + &
                                                & garb_array(1:lpSet%nSrcs,1:lpSet%nSpeciesTrn)
    !$OMP END CRITICAL(transf_maps_work)
    call free_work_array(garb_array)
    call free_work_array(cncTrn)
    call free_work_array(cncAer)
    call free_work_array(cncSL)


    !$OMP END PARALLEL
    
  end subroutine transform_lagrangian_part




  !************************************************************************************

  subroutine transform_deposited(mapDryDep, mapWetDep, &
                               & met_buf, disp_buf, meteo_input, &
                               & pHorizInterpStruct, pVertInterpStruct, &
                               & ifHorizInterp, ifVertInterp,&
                               & chemRules, seconds, now)
    !
    ! Transformations of deposited stuff: so far, only for radioactive 
    !
    implicit none
    
    ! Improted parameters
    type(Tmass_map), intent(inout) :: mapDryDep, mapWetDep
    type(Tfield_buffer), pointer :: met_buf, disp_buf
    type(Tmeteo_input), intent(in) :: meteo_input
    type(THorizInterpStruct), pointer :: pHorizInterpStruct
    type(TVertInterpStruct), pointer :: pVertInterpStruct
    logical, intent(in) :: ifHorizInterp, ifVertInterp
    type(Tchem_rules), intent(in) :: chemRules
    type(silja_time), intent(in) :: now
    real, intent(in) :: seconds

    ! Local variables
    integer :: ix, iy, iSrc, iTransf, iSpecies, status
    real, dimension(1:max_quantities) :: metdat
    logical :: print_it

!call msg_warning('Skipping chemistry')
!return

    ! Garbage comes as an array. For parallelisation, we temporarily use 
    ! a real array which will be private for everyone.

    if (chemRules%nTransformations + chemRules%nAerosolDynamics == 0) return

    !
    ! The main cycle over the grid
    !

    !call msg('Total N before transformation:', get_total_n(mapTransport))

    !$OMP PARALLEL if (if_OMP_chemistry) DEFAULT(SHARED) PRIVATE(metdat, ix, iy, itransf, isrc, print_it)
    
    !$OMP DO
    do iy = 1, mapDryDep%ny
      do ix = 1, mapDryDep%nx
        if (error) cycle
        call fill_in_meteo_input(meteo_input, metdat, &
                               & met_buf, disp_buf, &
                               & ix, iy, 1, &
                               & pHorizInterpStruct, pVertInterpStruct, &  ! meteo to disp
                               & ifHorizInterp, ifVertInterp)
        do iSrc = 1, mapDryDep%nSrc
          do iTransf = 1, chemRules%nTransformations
            select case(chemrules%iTransformTypes(itransf))
              case(transformation_passive)
              case(transformation_pollen)
              case (transformation_sulphur_dmat)
              case (transformation_acid_basic)
              case (transformation_cbm4)
              case(int_missing, transformation_inert_PM)
                  
              case (transformation_radioactive)
                !
                ! Radioactive decay is defined also for the deposited species
                !
                call transform_radioactive(mapDryDep%arM(:,isrc,1,ix,iy), &
                                         & chemRules%rulesRadioactive, &
                                         & metdat(:),&
                                         & seconds, print_it)
                call transform_radioactive(mapWetDep%arM(:,isrc,1,ix,iy), &
                                         & chemRules%rulesRadioactive, &
                                         & metdat(:),&
                                         & seconds, print_it)
              case default
                call msg('Unknown transformation:',chemrules%iTransformTypes(itransf))
                call set_error('Unknown transformation','transform_deposited')
                cycle
            end select
            if (error) cycle
          end do      ! iTransf
          if(error)then
            call set_error('Trouble with deposited-stuff chemistry','transform_deposited')
            call msg('ix,iy,iSrc:' + fu_str(ix) + ',' + fu_str(iy),iSrc)
          endif

        end do     ! iSrc
      end do    ! nx
    end do   ! ny
    !$OMP END DO


    !$OMP END PARALLEL

  end subroutine transform_deposited
    

  !***********************************************************************************

  subroutine make_wet_deposition_lagr()
    implicit none
    
    call set_error('Wet depostion not implemented', 'make_wet_deposition_lagr')
    return

  end subroutine make_wet_deposition_lagr


  !************************************************************************************
  !************************************************************************************
  !
  ! Global chemical rules
  !
  !************************************************************************************
  !************************************************************************************

  !*********************************************************************************

  subroutine set_chemistry_rules(nlTransf, nlStdSetup, rulesChemistry, ifOldFileFormat, if_lagr_present)
    !
    ! Initialises the basic chemical rules for the run
    !
    implicit none
    !
    ! Imported parameters
    type(Tsilam_namelist), intent(in) :: nlTransf, nlStdSetup
    type(Tchem_rules), intent(out) :: rulesChemistry
    logical, intent(in)  :: ifOldFileFormat
    logical, intent(in)  :: if_lagr_present
    !
    ! Local variables
    integer :: iTmp, jTmp, kTmp
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: items
    logical :: ifOK
    !
    !ATTENTION. Order in Dynamics calls is fixed and follows these arrays:
    character(len=17), dimension(5), parameter :: chDynamics = &
                                                & (/'NONE             ', 'SIMPLE           ', &
                                                  & 'MIDDLE_ATMOSPHERE', 'BASIC            ', &
                                                  & 'VBS              '/)
    integer, dimension(5), parameter :: iDynamics = &
                                                & (/int_missing, aerosol_dynamics_simple, &
                                                  & aerosol_dynamics_Mid_Atmosph, aerosol_dynamics_basic, &
                                                  & aerosol_dynamics_VBS/)
    character(len=fnlen) :: strTmp
    character(len=*), parameter :: sub_name='set_chemistry_rules'

    nullify(items)

    rulesChemistry%rulesAerosol = fu_set_aerosol_rules(nlTransf)
    if(error)return
    
    allocate(rulesChemistry%rulesOpticDens, &
           & rulesChemistry%ChemRunSetup, &
           & rulesChemistry%rulesDeposition, stat=iTmp)
    if (iTmp /= 0) then
      call set_error('Allocate failed', sub_name)
      return
    end if

    rulesChemistry%iTransformTypes(:) = int_missing
    rulesChemistry%iWhomToApplyTransform(:) = int_missing
    rulesChemistry%iWhomToApplyAerDyn(:) = int_missing

    !
    ! 1. Include chemical transformations.
    !
    call get_items(nlTransf, 'transformation', items, rulesChemistry%nTransformations)
    do iTmp = 1, rulesChemistry%nTransformations
      
      if(index(fu_str_u_case(fu_content(items(iTmp))),'PASSIVE') == 1)then
        rulesChemistry%iTransformTypes(iTmp) = transformation_passive
        call set_chem_rules_passive(nlTransf, rulesChemistry%rulesPassive)
        
      elseif(index(fu_str_u_case(fu_content(items(iTmp))),'POLLEN') == 1)then
        rulesChemistry%iTransformTypes(iTmp) = transformation_pollen
        call set_chem_rules_pollen(nlTransf, rulesChemistry%rulesPollen)

      elseif(index(fu_str_u_case(fu_content(items(iTmp))),'PM_GENERAL') == 1)then
        rulesChemistry%iTransformTypes(iTmp) = transformation_inert_PM
        call set_chemRules_PM_general(nlTransf, rulesChemistry%rulesPM)

      elseif(index(fu_str_u_case(fu_content(items(iTmp))),'DMAT_SULPHUR') == 1)then
        rulesChemistry%iTransformTypes(iTmp) = transformation_sulphur_dmat
        call set_chem_rules_sulphur_dmat(nlTransf, nlStdSetup, rulesChemistry%rulesSulphurDMAT)

      elseif(index(fu_str_u_case(fu_content(items(iTmp))),'ACID_BASIC') == 1)then
        rulesChemistry%iTransformTypes(iTmp) = transformation_acid_basic
        call set_chem_rules_acidBasic(nlTransf, nlStdSetup, rulesChemistry%rulesAcidBasic)

     ! Should be checked before CB4 
      elseif(index(fu_str_u_case(fu_content(items(iTmp))),'CB4_STRATO_SOA') == 1)then
        rulesChemistry%iTransformTypes(iTmp) = transformation_cbm42_strato_SOA
        call set_chem_rules_CBM42_strato_SOA(nlTransf, rulesChemistry%rulesCBM42_strato_SOA)
      elseif(index(fu_str_u_case(fu_content(items(iTmp))),'CB4_STRATO') == 1)then
        rulesChemistry%iTransformTypes(iTmp) = transformation_cbm42_strato
        call set_chem_rules_CBM42_strato(nlTransf, rulesChemistry%rulesCBM42_strato)
      elseif(index(fu_str_u_case(fu_content(items(iTmp))),'CB4_SOA') == 1)then
        rulesChemistry%iTransformTypes(iTmp) = transformation_cbm4_SOA
        call set_chem_rules_CBM4_SOA(nlTransf, rulesChemistry%rulesCBM4_SOA)
      elseif(index(fu_str_u_case(fu_content(items(iTmp))),'CB5_SOA') == 1)then
        rulesChemistry%iTransformTypes(iTmp) = transformation_cbm5_SOA
        call set_chem_rules_CBM5_SOA(nlTransf, rulesChemistry%rulesCBM5_SOA)
      elseif(index(fu_str_u_case(fu_content(items(iTmp))),'CB5_STRATO_SOA') == 1)then
        rulesChemistry%iTransformTypes(iTmp) = transformation_cbm5_strato_SOA
        call set_chem_rules_CBM5_strato_SOA(nlTransf, rulesChemistry%rulesCBM5_strato_SOA)

      elseif(index(fu_str_u_case(fu_content(items(iTmp))),'CB4') == 1)then
        rulesChemistry%iTransformTypes(iTmp) = transformation_cbm4
        call set_chem_rules_CBM4(nlTransf, rulesChemistry%rulesCBM4)
        

      elseif(index(fu_str_u_case(fu_content(items(iTmp))),'RADIOACTIVE') == 1)then
        rulesChemistry%iTransformTypes(iTmp) = transformation_radioactive
        call set_chem_rules_radioactive(nlTransf, rulesChemistry%rulesRadioactive)

      else
        call set_error('Transfomation not supported:' + fu_content(items(iTmp)), &
                     & sub_name)
        return

      endif   ! transformation
      !
      ! Now set the type of dynamics this transformation should be applied
      !
      if(ifOldFileFormat)then
        rulesChemistry%iWhomToApplyTransform(iTmp) = eulerian_flag
      else
        if(index(fu_str_u_case(fu_content(items(iTmp))),'LAGRANGIAN') > 0)then
          rulesChemistry%iWhomToApplyTransform(iTmp) = lagrangian_flag
        elseif(index(fu_str_u_case(fu_content(items(iTmp))),'EULERIAN') > 0)then
          rulesChemistry%iWhomToApplyTransform(iTmp) = eulerian_flag
        elseif(index(fu_str_u_case(fu_content(items(iTmp))),'HYBRID') > 0)then
          rulesChemistry%iWhomToApplyTransform(iTmp) = hybrid_flag
        else
          call set_error('Unknown type of dynamics in:' + fu_content(items(iTmp)),sub_name)
          call msg_warning('Dynamics must be LAGRANGIAN or EULERIAN or HYBRID')
          return
        endif
      endif
    end do  ! list of transformations

    !
    ! Aerosol dynamics. There can be nothing or something - but here we check for impossible combinations
    ! Note that the selected dynamics will be applied to both Lagrangian and Eulerian environments
    ! Later this limitation can be removed
    !
    rulesChemistry%iAerosolDynTypes(:) = int_missing
    call get_items(nlTransf, 'aerosol_dynamics', items, rulesChemistry%nAerosolDynamics)

    if(rulesChemistry%nAerosolDynamics <= 0)then
      !
      ! User does not want any aerosol dynamics. However, acid_basic and sulphur_dmat
      ! put some of their stuff to short-living species, so at least simple dynamics
      ! must be involved if these are in game
      !
      if(fu_quantity_in_quantities(transformation_sulphur_dmat, rulesChemistry%iTransformTypes) .or. &
       & fu_quantity_in_quantities(transformation_acid_basic, rulesChemistry%iTransformTypes))then
        call set_error('Acid_basic and sulphur_dmat transformations require aerosol dynamics', &
                     & sub_name)
        return
      endif
    else
      !
      ! Check the supported dynamics types: the namelist must have only them
      !
      ifOK = .false.
      do jTmp = 1, rulesChemistry%nAerosolDynamics
        do iTmp = 1, size(chDynamics)
          if(index(fu_str_u_case(fu_content(items(jTmp))),trim(chDynamics(iTmp))) > 0)then
            ifOK = .true.
            exit
          endif
        end do
        if(.not. ifOK)then
          call set_error('Unsupported aerosol dynamics:'+fu_content(items(jTmp)),sub_name)
          return
        endif
      end do
      
      !
      ! Find the right ones - for Lagrangian and Eulerian parts
      ! Note that the dynamics can order the species by filling-in the list of expected_species.
      ! If they are flexible - nothing to be written there.
      !
      !
      ! ATTENTION.
      ! Simple and Mid-Atmosphere dynamics are compatible in principle but their order is prescribed:
      ! Mid-Atm dynamics must be called first.
      !
      kTmp = 1
      do iTmp = 1, size(chDynamics)
        do jTmp = 1, rulesChemistry%nAerosolDynamics
          if(index(fu_str_u_case(fu_content(items(jTmp))),trim(chDynamics(iTmp))) > 0)then
            rulesChemistry%iAerosolDynTypes(kTmp) = iDynamics(iTmp)
            kTmp = kTmp + 1
          endif
        end do
      end do  ! cycle over supported aerosol dynamics
      !
      ! Now, having the dynamics calls in right order, initialise the stuff
      !
      do iTmp = 1, rulesChemistry%nAerosolDynamics
        select case (rulesChemistry%iAerosolDynTypes(iTmp))
          case(aerosol_dynamics_basic)
            call set_rules_AerDynBasic(nlTransf, rulesChemistry%rulesAerDynBasic)
            call init_AerDynBasic(rulesChemistry%rulesAerDynBasic, rulesChemistry%expected_species)
          case(aerosol_dynamics_VBS)
            call set_rules_AerDynVBS(nlTransf, rulesChemistry%rulesAerDynVBS)
            call init_AerDynVBS(rulesChemistry%rulesAerDynVBS, rulesChemistry%expected_species)
          case(aerosol_dynamics_simple)
            call set_rules_aerDynSimple(nlTransf, rulesChemistry%rulesAerDynSimple)
            call init_aerDynSimple(rulesChemistry%expected_species)
          case(aerosol_dynamics_Mid_Atmosph)
            call set_rules_AerDynMidAtm(nlTransf, rulesChemistry%rulesAerDynMidAtmosph)
            call init_AerDynMidAtm(rulesChemistry%rulesAerDynMidAtmosph, rulesChemistry%expected_species)
          case default
        end select  ! types of aerosol dynamics
        !
        ! Now set the type of dynamics this aerosol dynamics should be applied
        !
        if(ifOldFileFormat)then
          rulesChemistry%iWhomToApplyAerDyn(iTmp) = eulerian_flag
        else
          if(index(fu_str_u_case(fu_content(items(iTmp))),'LAGRANGIAN') > 0)then
            rulesChemistry%iWhomToApplyAerDyn(iTmp) = lagrangian_flag
          elseif(index(fu_str_u_case(fu_content(items(iTmp))),'EULERIAN') > 0)then
            rulesChemistry%iWhomToApplyAerDyn(iTmp) = eulerian_flag
          elseif(index(fu_str_u_case(fu_content(items(iTmp))),'HYBRID') > 0)then
            rulesChemistry%iWhomToApplyAerDyn(iTmp) = hybrid_flag
          else
            call set_error('Unknown type of dynamics in:' + fu_content(items(iTmp)),sub_name)
            call msg_warning('Dynamics must be LAGRANGIAN or EULERIAN or HYBRID')
            return
          endif
        endif   ! old or new format
      enddo  ! 1:nAerosolDynamics
    endif ! how many aerosol dynamics?

    if(rulesChemistry%nAerosolDynamics < 0) rulesChemistry%nAerosolDynamics = 0

    ! The optical rules are set later
    !
    call set_missing(rulesChemistry%rulesOpticDens)

    !-----------------------------------------------------------------------------------
    !
    ! Deposition rules are set as part of chemical ones
    !
    ! The standard deposition will be applied only for those species, which are not covered by the 
    ! specific chemistry-related procedures
    !
    allocate(rulesChemistry%rulesDeposition%iDepositionType(rulesChemistry%nTransformations+1), &
           & stat = iTmp)
    if(iTmp /= 0)then
      call set_error('Failed to allocate memory for deposition types',sub_name)
      return
    endif
    rulesChemistry%rulesDeposition%iDepositionType = int_missing
    rulesChemistry%rulesDeposition%nDepositionTypes = 0

    do iTmp = 1, rulesChemistry%nTransformations
      select case(rulesChemistry%iTransformTypes(iTmp))

        case(transformation_passive)
          if(fu_if_specific_deposition(rulesChemistry%rulesPassive))then
            rulesChemistry%rulesDeposition%nDepositionTypes = &
                                 & fu_merge_integer_to_array(transformation_passive, &
                                                   & rulesChemistry%rulesDeposition%iDepositionType)
          endif
          
        case(transformation_pollen)
          if(fu_if_specific_deposition(rulesChemistry%rulesPollen))then
            rulesChemistry%rulesDeposition%nDepositionTypes = &
                                 & fu_merge_integer_to_array(transformation_pollen, &
                                                   & rulesChemistry%rulesDeposition%iDepositionType)
          endif

        case(transformation_inert_PM)
          if(fu_if_specific_deposition(rulesChemistry%rulesPM))then
            rulesChemistry%rulesDeposition%nDepositionTypes = &
                                 & fu_merge_integer_to_array(transformation_inert_PM, &
                                                   & rulesChemistry%rulesDeposition%iDepositionType)
          endif

        case (transformation_sulphur_dmat)
          if(fu_if_specific_deposition(rulesChemistry%rulesSulphurDMAT))then
            rulesChemistry%rulesDeposition%nDepositionTypes = &
                                 & fu_merge_integer_to_array(transformation_sulphur_dmat, &
                                                   & rulesChemistry%rulesDeposition%iDepositionType)
          endif

        case (transformation_acid_basic)
          if(fu_if_specific_deposition(rulesChemistry%rulesAcidBasic))then
            rulesChemistry%rulesDeposition%nDepositionTypes = &
                                 & fu_merge_integer_to_array(transformation_acid_basic, &
                                                   & rulesChemistry%rulesDeposition%iDepositionType)
          endif

        case (transformation_radioactive)
          if(fu_if_specific_deposition(rulesChemistry%rulesRadioactive))then
            rulesChemistry%rulesDeposition%nDepositionTypes = &
                                 & fu_merge_integer_to_array(transformation_radioactive, &
                                                   & rulesChemistry%rulesDeposition%iDepositionType)
          endif

        case (transformation_cbm4)
          if(fu_if_specific_deposition(rulesChemistry%rulesCBM4))then
            rulesChemistry%rulesDeposition%nDepositionTypes = &
                                 & fu_merge_integer_to_array(transformation_cbm4, &
                                                   & rulesChemistry%rulesDeposition%iDepositionType)
          endif
                  
        case (transformation_cbm4_SOA)
          if(fu_if_specific_deposition(rulesChemistry%rulesCBM4_SOA))then
            rulesChemistry%rulesDeposition%nDepositionTypes = &
                                 & fu_merge_integer_to_array(transformation_cbm4_SOA, &
                                                   & rulesChemistry%rulesDeposition%iDepositionType)
          endif
          

        case (transformation_cbm42_strato)
          if(fu_if_specific_deposition(rulesChemistry%rulesCBM42_strato))then
            rulesChemistry%rulesDeposition%nDepositionTypes = &
                                 & fu_merge_integer_to_array(transformation_cbm42_strato, &
                                                           & rulesChemistry%rulesDeposition%iDepositionType)
          endif
         
        case (transformation_cbm42_strato_SOA)
          if(fu_if_specific_deposition(rulesChemistry%rulesCBM42_strato_SOA))then
            rulesChemistry%rulesDeposition%nDepositionTypes = &
                                 & fu_merge_integer_to_array(transformation_cbm42_strato_SOA, &
                                                           & rulesChemistry%rulesDeposition%iDepositionType)
          endif

        case (transformation_cbm5_SOA)
          if(fu_if_specific_deposition(rulesChemistry%rulesCBM5_SOA))then
            rulesChemistry%rulesDeposition%nDepositionTypes = &
                                 & fu_merge_integer_to_array(transformation_cbm5_SOA, &
                                                           & rulesChemistry%rulesDeposition%iDepositionType)
          endif

        case (transformation_cbm5_strato_SOA)
          if(fu_if_specific_deposition(rulesChemistry%rulesCBM5_strato_SOA))then
            rulesChemistry%rulesDeposition%nDepositionTypes = &
                                 & fu_merge_integer_to_array(transformation_cbm5_strato_SOA, &
                                                           & rulesChemistry%rulesDeposition%iDepositionType)
          endif

        case default
          call msg('Transfomation not supported:',rulesChemistry%iTransformTypes(iTmp))
          call set_error('Transfomation not supported',sub_name)
      end select
      !if(error)return
    end do  ! transformation types
!    call unset_error(sub_name)
    !
    ! Standard deposition is always computed
    !
    rulesChemistry%rulesDeposition%nDepositionTypes = &
                     & fu_merge_integer_to_array(deposition_standard, &
                               & rulesChemistry%rulesDeposition%iDepositionType)
    call set_deposition_rules(nlTransf, rulesChemistry%rulesDeposition, if_lagr_present)
    if(error)return


    ! Extra species. If not defined, no species added.
    call get_items(nlTransf, 'auxiliary_cocktail', items, rulesChemistry%nAuxCocktails)
    do iTmp = 1, rulesChemistry%nAuxCocktails
        rulesChemistry%auxCocktName(iTmp) = fu_content(items(iTmp))
    enddo
    !
    ! 
    !
    strTmp = fu_str_u_case(fu_content(nlStdSetup,'reference_4_low_mass_threshold'))
    select case(strTmp)
      case('EMISSION')
        rulesChemistry%LowMassThresh = LowMassThreshUseEmission
        call msg("LowMassThresh UseEmission")
      case('DEFAULT')
        rulesChemistry%LowMassThresh = LowMassThreshUseDefault
        call msg("LowMassThresh UseDefault")
      case('CONST')
        rulesChemistry%LowMassThresh = LowMassThreshUseConst
        call msg("LowMassThresh UseConst")
      case('')
        call msg("Default behaviour: LowMassThresh UseEmission") !!Backward-compatible
        rulesChemistry%LowMassThresh = LowMassThreshUseEmission
      case default
        call msg_warning('empty reference_4_low_mass_threshold, must be EMISSION or DEFAULT or CONST', &
                       & sub_name)
        call set_error("Strange reference_4_low_mass_threshold: "//trim(strTmp), sub_name)
        return
    end select

    ! For CB4 strato, initialize photolysis
    !
!!$    if (fu_index(transformation_cbm42_strato, chemRules%iTransformTypes) > 0 &
!!$      & .or. fu_index(transformation_cbm4, chemRules%iTransformTypes) > 0 ) then
    if (fu_index(transformation_cbm42_strato, rulesChemistry%iTransformTypes) > 0 &
      & .or. fu_index(transformation_cbm42_strato_SOA, rulesChemistry%iTransformTypes) > 0 &        
      & .or. fu_index(transformation_cbm4_SOA, rulesChemistry%iTransformTypes) > 0 &
      & .or. fu_index(transformation_cbm5_SOA, rulesChemistry%iTransformTypes) > 0 &
      & .or. fu_index(transformation_cbm5_strato_SOA, rulesChemistry%iTransformTypes) > 0) then
      rulesChemistry%need_photo_lut = .true.


      !Check if one requires a dynamic albedo or a static one
      if (fu_str_u_case(fu_content(nlTransf,'use_dynamic_albedo')) == 'YES') then
         rulesChemistry%useDynamicAlbedo = .true.
         call msg('Photolysis rates will be calculated using dynamic (forecast) albedo from meteo!')
      else
         rulesChemistry%useDynamicAlbedo = .false.
         rulesChemistry%defaultStaticAlbedo = fu_content_real(nlTransf, 'default_static_albedo')
         if (rulesChemistry%defaultStaticAlbedo .eps. real_missing) then
            rulesChemistry%defaultStaticAlbedo = 0.3
            call msg('No default albedo found, using predefined value 0.3!') 
         elseif (.not. abs(rulesChemistry%defaultStaticAlbedo-0.5) <= 0.5) then
            call msg('Strange default_static_albedo, must be real [0..1]')
            call report(nlTransf)
            call set_error('default_static_albedo, must be real [0..1]',sub_name)
            return
         end if
         call msg('Photolysis rates will be calculated using albedo = ', rulesChemistry%defaultStaticAlbedo)
      end if 

      ! 
      strTmp = fu_content(nlStdSetup, 'photolysis_data_file')
      rulesChemistry%filename_photo_lut = fu_process_filepath(strTmp, must_exist=.true.)
      if (error) then
        call set_error("Error after photolysis_data_file = '"//trim(strTmp)//"'", sub_name)
        return
      endif
        
      call init_photolysis_lut(rulesChemistry%filename_photo_lut)
      if (error) then
        call set_error("Error after init_photolysis_lut", sub_name)
        return
      endif
      

      strTmp = fu_str_u_case(fu_content(nlTransf,'photolysis_affected_by_aod'))
      if( strTmp == 'YES') then
        strTmp = fu_content(nlTransf,'photolysis_AOD_wavelength')
        if (trim(strTmp) == "") then
          call msg("Using default photolysis_AOD_wavelength")
          strTmp = "550 nm"
        endif
        rulesChemistry%photoAODwavelength = fu_set_named_value(strTmp)
        call msg("photolysis_AOD_wavelength, nm", rulesChemistry%photoAODwavelength*1e9)
        rulesChemistry%ifPhotoAOD = .True.
        call msg("photolysis_affected_by_aod = YES")
      elseif (strTmp == 'NO') then
        rulesChemistry%ifPhotoAOD = .False.
        rulesChemistry%photoAODwavelength = real_missing
        call msg("Photolysis ignores aerosols: photolysis_affected_by_aod = NO")
    else
         call msg("Usage: photolysis_affected_by_aod = (YES|NO)")
         call set_error("Unknown photolysis_affected_by_aod value: '"//trim(strtmp)//"'", sub_name)
      endif

      strTmp = fu_str_u_case(fu_content(nlTransf,'photolysis_affected_by_o3col'))
      if(  strTmp == 'YES' .or.  strTmp == 'MASS_MAP') then !! Backward compatible
        rulesChemistry%PhotoO3col = mass_map
        call msg("photolysis_affected_by_o3col = MASS_MAP") 
      elseif (strTmp == 'METEO') then
        rulesChemistry%PhotoO3col = meteo_column
        call msg("photolysis_affected_by_o3col = METEO") 
      elseif  ((strTmp == 'STANDARD_ATMOSPHERE') .or. (strTmp == '') .or. (strTmp == 'NO')) then
        rulesChemistry%PhotoO3col = standard_atmosphere
        call msg("photolysis_affected_by_o3col = STANDARD_ATMOSPHERE") 
      else
        call msg("Usage: photolysis_affected_by_o3col = (MASS_MAP|METEO|STANDARD_ATMOSPHERE)")
        call set_error("Unknown photolysis_affected_by_o3col value: '"//trim(strtmp)//"'", sub_name)
      end if

      if (fu_str_u_case(fu_content(nlTransf,'cloud_model_for_photolysis')) == 'SIMPLE_CLOUD') then
        rulesChemistry%cloud_model_for_photolysis = simple_cloud
        call msg("Simple water cloud model for photolysis")
      else if (fu_str_u_case(fu_content(nlTransf,'cloud_model_for_photolysis')) == 'DETAILED_CLOUD') then
        rulesChemistry%cloud_model_for_photolysis = detailed_cloud
        call msg("Detailed water cloud model for photolysis")
      else if (fu_str_u_case(fu_content(nlTransf,'cloud_model_for_photolysis')) == 'FAKE_CLOUD') then
        rulesChemistry%cloud_model_for_photolysis = fake_cloud
        call msg("Fake water cloud model for photolysis")
      end if
    else
      rulesChemistry%cloud_model_for_photolysis = int_missing
      rulesChemistry%PhotoO3col = int_missing
      rulesChemistry%ifPhotoAOD = .False.
      rulesChemistry%photoAODwavelength = real_missing

      rulesChemistry%need_photo_lut = .false.
    end if

    strTmp = fu_str_u_case(fu_content(nlTransf,'adjust_cell_volume_for_ones'))
    if(  strTmp == 'YES') then !! Backward compatible
      rulesChemistry%ifOnesAdjust = .True.
      call msg("Using ONES to adjust air density for chemistry")
    endif

    rulesChemistry%defined = silja_true

  end subroutine set_chemistry_rules


  !*********************************************************************************

  subroutine set_chemrules_missing(chemrules)
    implicit none
    type(Tchem_rules), intent(out) :: chemrules
    
    call set_missing(chemrules%rulesAerosol)
    nullify(chemrules%rulesOpticDens, chemrules%chemRunSetup, chemrules%rulesDeposition,&
          & chemrules%low_mass_trsh)
    chemrules%defined = silja_false

  end subroutine set_chemrules_missing


  !*********************************************************************************

  subroutine set_low_mass_threshold(chemRules, &
                                  & src, &
                                  & transport_species, emission_species, &
                                  & nspecies_transport, &
                                  & residence_interval, &
                                  & period_to_compute, &
                                  & timestep, &
                                  & iComputationAccuracy, &
                                  & earliest_start, &
                                  & ifForwardRun, ifFromEmission, & ! frw/adj and emis/observations
                                  & DA_RMSE_species, nDA_ObservationPoints, &
                                  & ifDefineLagrangianParticle)
    type(Tchem_rules), intent(inout) :: chemRules
    type(silam_source), intent(in) :: src
    type(silam_species), dimension(:), pointer :: transport_species, emission_species
    integer, intent(in) :: nspecies_transport, iComputationAccuracy
    integer, intent(in) :: nDA_ObservationPoints
    type(silja_interval), intent(in) :: residence_interval, period_to_compute, timestep
    type(silja_time), intent(in) :: earliest_start
    real, dimension(:), intent(in) :: DA_RMSE_species
    logical, intent(in) :: ifForwardRun, ifFromEmission, ifDefineLagrangianParticle

    ! Local variables
    integer :: i, istat
    real :: fTmp
    real(8), dimension(max_species) :: amounts_r8
    real, dimension(max_species) :: amounts_my, amounts_all
    logical, save :: ifFirst = .true.
    !
    ! Attention! 
    ! The sub can be called several times in case of data assimilation. Careful with memory
    ! allocation and unnecesary repeating of some calculations.
    !
    if(ifFirst)then
      !
      ! Allocate memory, calculate thresholds, store them to thresholds_saved array.
      !
      ! But first, a stupidity test
      !
      allocate(chemrules%low_mass_trsh_frw(nspecies_transport), &
             & chemrules%low_cnc_trsh_frw(nspecies_transport), &
             & chemrules%low_mass_trsh_adj(nspecies_transport), &
             & chemrules%low_cnc_trsh_adj(nspecies_transport), &
             & stat=istat)
      if(fu_fails(istat == 0, 'Allocate failed', 'set_low_mass_threshold'))return
      chemrules%low_cnc_trsh_frw = -1.0
      chemrules%low_mass_trsh_frw = -1.0
      chemrules%low_cnc_trsh_adj = -1.0
      chemrules%low_mass_trsh_adj = -1.0
      amounts_r8(:) = 0.0
      amounts_my(:) = 0.0
      amounts_all(:) = 0.0

      !      chemrules%low_mass_trsh_frw(:) = 1e-23
      !chemrules%low_cnc_trsh_frw(:) = 1e-37

      if(ifDefineLagrangianParticle)then
        allocate(chemrules%mass_lagr_particle(nspecies_transport), stat=istat)
        if(fu_fails(istat == 0, 'Lagrangian allocation failed', 'set_low_mass_threshold'))return
        chemrules%mass_lagr_particle = -1.0
      end if
      ifFirst = .false.
    endif  ! ifFirst

    !
    ! Having space allocated, get the required threshold
    !
    if(ifForwardRun)then
      if(chemrules%low_mass_trsh_frw(1) < 0)then   ! not available, calculate
        call calculate_thresholds(chemrules%low_cnc_trsh_frw, chemrules%low_mass_trsh_frw)
        if(error)return
      endif
    else
      if(chemrules%low_mass_trsh_adj(1) < 0)then   ! not available, calculate

        call calculate_thresholds(chemrules%low_cnc_trsh_adj, chemrules%low_mass_trsh_adj)
        if(error)return
        
        ! Extrapolation: in da-runs, some thresholds aree not available since species are not observed
        ! get the scaling for observed species and extrapolate it
        fTmp = 0.0
        do i = 1, nSpecies_transport
          if(chemrules%low_cnc_trsh_adj(i) > 0.0) then
            fTmp = max(fTmp, chemrules%low_cnc_trsh_frw(i) / chemrules%low_cnc_trsh_adj(i))
          end if
        end do
        ! and now use it for missing components
        do i = 1, nSpecies_transport
          if(chemrules%low_cnc_trsh_adj(i) == 0.0) then
              chemrules%low_cnc_trsh_adj(i) = chemrules%low_cnc_trsh_frw(i) / fTmp
              chemrules%low_mass_trsh_adj(i) = chemrules%low_mass_trsh_frw(i) / fTmp
            call msg('Low-concentration threshold by scaling, ' // fu_substance_name(transport_species(i)), &
                   & chemrules%low_cnc_trsh_adj(i))
           endif
        end do
      endif
      
    endif  ! ifForwardRun

    call point_low_mass_thresholds(chemrules, ifForwardRun)

    ! No matter how the thresholds were calculated, there must be no zero or negative thresholds
    !
    do i = 1, nSpecies_transport
      if(.not. chemrules%low_cnc_trsh(i) > 0.0)then
        call msg('Non-positive low-cnc threshold for species:'+fu_str(i), &
                                                     & chemrules%low_cnc_trsh(1:nSpecies_transport))
        call set_error('Not all low-concentration thresholds are positive','set_low_mass_threshold')
        return
      endif
    end do
    !
    ! Some species must exist if another one exists, so we will ensure that by reducing its threshold
    !
    call msg('Checking the consistency of the low-mass thresholds')
    do i = 1, chemrules%nTransformations
      select case(chemrules%iTransformTypes(i))
        case(transformation_cbm4)
          call check_low_mass_threshold_cb4(chemrules%rulesCBM4, transport_species, &
                                          & nSpecies_transport, ifForwardRun, &
                                          & chemrules%low_cnc_trsh, chemrules%low_mass_trsh)
!The following is not needed after PAR-OLE tweaks:
!        case(transformation_cbm4_SOA)
!          call check_low_mass_threshold_cbm4_SOA(chemrules%rulesCBM4_SOA, transport_species, &
!                                          & nSpecies_transport, ifForwardRun, &
!                                          & chemrules%low_cnc_trsh, chemrules%low_mass_trsh)          
!        case(transformation_cbm42_strato)
!          call check_low_mass_threshold_cbm42_strato(chemrules%rulesCBM42_strato, transport_species, &
!                                             & nSpecies_transport, ifForwardRun, &
!                                             & chemrules%low_cnc_trsh, chemrules%low_mass_trsh)
!        case(transformation_cbm42_strato_SOA)
!          call check_low_mass_threshold_cbm42_strato_SOA(chemrules%rulesCBM42_strato_SOA, transport_species, &
!                                             & nSpecies_transport, ifForwardRun, &
!                                             & chemrules%low_cnc_trsh, chemrules%low_mass_trsh)
!        case(transformation_cbm5_SOA)
!          call check_low_mass_threshold_cbm5_SOA(chemrules%rulesCBM5_SOA, transport_species, &
!                                             & nSpecies_transport, ifForwardRun, &
!                                             & chemrules%low_cnc_trsh, chemrules%low_mass_trsh)
!        case(transformation_cbm5_strato_SOA)
!          call check_low_mass_threshold_cbm5_strato_SOA(chemrules%rulesCBM5_strato_SOA, transport_species, &
!                                             & nSpecies_transport, ifForwardRun, &
!                                             & chemrules%low_cnc_trsh, chemrules%low_mass_trsh)
        case(transformation_cbm4_SOA,transformation_cbm42_strato,transformation_cbm42_strato_SOA, &
           & transformation_cbm5_SOA, transformation_cbm5_strato_SOA, &
           & transformation_passive, transformation_pollen, &
           & transformation_inert_PM, transformation_sulphur_dmat, &
           & transformation_acid_basic, transformation_radioactive, &
           & int_missing ) !int_missing is perfectly valid here, as global_chemical_init might set it
          ! Do nothing - these folks are fine with any combinations of any thresholds
        case default
          call set_error('Unknown transformation type:'+fu_str(chemrules%iTransformTypes(i)),&
                       & 'set_low_mass_threshold')
      end select
    end do   ! iTrn
    !
    ! Report the thresholds
    !
    if(ifForwardRun)then
      call msg('Low-concentration thresholds from emission data:')
    else
      call msg('Low-concentration thresholds from model-measurement discrepancy:')
    endif
    do i = 1, nspecies_transport
      call msg(fu_substance_name(transport_species(i)) + ':', chemrules%low_cnc_trsh(i))
    end do
    


    CONTAINS

    !=======================================================================================

    subroutine calculate_thresholds(pThresholdsCnc, pThresholdsMass)
      !
      ! Does the actual threshold calculation
    
      implicit none
      
      ! Imported parameter
      real, dimension(:), intent(inout) :: pThresholdsCnc, pThresholdsMass
      
      ! Local variables
      integer :: i
      real :: fCellSize
      real, parameter :: typical_cnc_factor = 1e-3
      logical :: ifOk
      
      fCellSize = fu_cell_size(wholeMPIdispersion_grid, fs_wholeMPIdispersion/2)
      ! No need to touch sources or do exchange at all
      if (chemrules%LowMassThresh == LowMassThreshUseConst) then
        pThresholdsCnc(1:nspecies_transport) = fLowMassThreshConstant / (fCellSize * 100.)

        pThresholdsMass(1:nspecies_transport) = pThresholdsCnc(1:nspecies_transport) * fCellSize * 100.0
        return
      endif
      !
      ! Stuff an be injected either from emission sources or from model-measurement discrepancy
      !
      if(ifFromEmission)then
        !
        ! For each transport species:
        ! 1. try to find the amount released from inventory sources, compute threshold, OR
        ! 2. try to get typical concentration from dynamic sources, compute threshold  OR
        ! 3. if both are available, use the smaller one, OR
        ! 4. if neither are available, use one from materials
        !
        ! 1. try to find the amount released from inventory sources, compute threshold
        !    note that it can be forbidden if we want default values
        !
        amounts_r8 = 0.0
        
        if(chemrules%LowMassThresh == LowMassThreshUseEmission)then
          call msg ("calling amounts_from_src_species_unit" )
          call amounts_from_src_species_unit(src, transport_species, nspecies_transport, &
                                           & amounts_r8, &
                                           & earliest_start, fu_abs(period_to_compute), &
                                           & emission_species, &
                                           & chemRules%ChemRunSetup%refEmis2Transp_mass)
          if(error)return

          if (smpi_global_tasks > 1) then
            call msg("Exchanging emitted amounts")
            amounts_my(nspecies_transport+1) = real(fu_nbr_of_disp_grd_cells(src))!  exchange it as well
            amounts_my(1:nspecies_transport) = real(amounts_r8(1:nspecies_transport))
            call  smpi_allreduce_add(amounts_my(1:nspecies_transport+1), &
                                   & amounts_all(1:nspecies_transport+1), &
                                   & smpi_adv_comm, ifOk)
            if (fu_fails(ifOk, "Amonts exchange failed","calculate_thresholds" )) return
          else
            amounts_all(1:nspecies_transport) = real(amounts_r8(1:nspecies_transport))
            amounts_all(nspecies_transport+1) = real(fu_nbr_of_disp_grd_cells(src))
          endif
          call msg("amounts_my(i)", amounts_my(1:nspecies_transport+1))
          call msg("amounts_all(i)", amounts_all(1:nspecies_transport+1))
          call msg("amounts_from_src_species_unit earliest_start: "+fu_str(earliest_start) )
          call msg("amounts_from_src_species_unit fu_abs(period_to_compute):"+ fu_str(fu_abs(period_to_compute)))
          call msg("amounts_from_src_species_unit end: "+fu_str(earliest_start + fu_abs(period_to_compute)) )
          call msg("nx_wholeMPIdispersion, ny_wholeMPIdispersion, fs_wholeMPIdispersion", (/nx_wholeMPIdispersion, ny_wholeMPIdispersion, fs_wholeMPIdispersion/))
          call msg("fCellSize,", fCellSize)

          do i = 1, nspecies_transport
            if(amounts_all(i) > 0) pThresholdsCnc(i) = &
                       & fu_compute_thresh_for_release(amounts_all(i), &
                                                & nint(amounts_all(nspecies_transport+1)), & !& fu_nbr_of_disp_grd_cells(src),&
                                                     & fu_abs(residence_interval), &
                                                     & fu_abs(period_to_compute), &
                                                     & fu_abs(timestep), &
                                                     & nx_wholeMPIdispersion, ny_wholeMPIdispersion, nz_dispersion, &
                                                     & fCellSize, &
                                                     & 1000.0, & ! layer thickness
                                                     & iComputationAccuracy)
          end do

        endif  ! rulesChemistry%ifLowMassThreshUseEmission

        ! 2. try to get typical concentration from dynamic sources, compute threshold  OR
        ! 3. if both are available, use the smaller one, OR
        !
        amounts_my(:) = 0.0
        call typical_cnc_from_src_species_unit(src, transport_species, nspecies_transport, amounts_my)
!        ! Do we need an exchange here????
!        if (smpi_global_tasks > 1) then
!          amounts_my(1:nspecies_transport) = real(amounts_r8(1:nspecies_transport))
!          call  smpi_allreduce_add(amounts_my(1:nspecies_transport), &
!                                 & amounts_all(1:nspecies_transport), &
!                                 & smpi_adv_comm, ifOk)
!          if (fu_fails(ifOk, "cnc exchange failed","calculate_thresholds" )) return
!        else
          amounts_all(1:nspecies_transport) = amounts_my(1:nspecies_transport)
!        endif

        do i = 1, nspecies_transport
          if (.not. (amounts_all(i) .eps. real_missing)) then
            if (pThresholdsCnc(i) > 0) then
              pThresholdsCnc(i) = min(pThresholdsCnc(i), typical_cnc_factor * amounts_all(i))
            else
              pThresholdsCnc(i) = amounts_all(i) * typical_cnc_factor
            end if
          end if
        end do

        ! 4. if neither are available, use one from materials
        !
        do i = 1, nspecies_transport
          if (pThresholdsCnc(i) < 0) then
            pThresholdsCnc(i) = fu_low_mass_threshold(fu_material(transport_species(i)))
            if (pThresholdsCnc(i) .eps. real_missing) pThresholdsCnc(i) = -1.0
            if (error) exit
          end if
        end do
        

        pThresholdsMass(1:nspecies_transport) = pThresholdsCnc(1:nspecies_transport) * fCellSize * 100.0

        


        if (.not. all(pThresholdsCnc > 0)) then
          call set_error('Not all low mass thresholds set', 'set_low_mass_threshold')
        end if
      else
        !
        ! Threshold comes from model-measurement discrepancy
        !
        do i = 1, nSpecies_transport
          pThresholdsCnc(i) = &
                       & fu_compute_thresh_for_release(DA_RMSE_species(i), &
                                                     & nDA_ObservationPoints, &
                                                     & fu_abs(residence_interval), &
                                                     & fu_abs(period_to_compute), &
                                                     & fu_abs(timestep), &
                                                     & nx_dispersion, ny_dispersion, nz_dispersion, &
                                                     & fCellSize, &
                                                     & 1000.0, & ! layer thickness
                                                     & iComputationAccuracy)
          amounts_r8(i) = DA_RMSE_species(i)   ! for the case of Lagrangian particles
          pThresholdsMass(i) = pThresholdsCnc(i) * fCellSize * 100.0
        end do
      endif  ! if threshold comes from emission or mdl-meas discrepancy

      ! Final step: having the low mass threshold, determine the initial mass of each species in 
      ! the lagrangian particle.
      ! Algorithm to try:
      ! Low-mass threshold uniformly distributed over the whole grid gives global minimal 
      ! noticeable mass. It is to be distributed over lagrangian particles taking into account
      ! typical cnc non-uniformity factor
      !
      if(ifDefineLagrangianParticle)then
        iStat = fu_number_of_lagr_particles(iComputationAccuracy)
        call msg('Approximate Mass in a single Lagrangian particle:')
        do i = 1, nspecies_transport
call msg('chemrules%low_mass_trsh(i)',pThresholdsMass(i))
call msg('nx*ny*nz', nx_dispersion * ny_dispersion * nz_dispersion)
call msg('min(1.,(fu_abs(residence_interval) / fu_abs(period_to_compute)))', min(1.,(fu_abs(residence_interval) / fu_abs(period_to_compute))))
call msg('accuracy, nop',iComputationAccuracy, iStat )
          if(amounts_r8(i) > 0)then
            chemRules%mass_lagr_particle(i) = fu_lagr_part_mass_for_release(real(amounts_r8(i)), &
                                                   & fu_nbr_of_disp_grd_cells(src),&
                                                   & fu_abs(residence_interval), &
                                                   & fu_abs(period_to_compute), &
                                                   & fu_abs(timestep), &
                                                   & iStat)
          else
            chemRules%mass_lagr_particle(i) = 1.0  ! anything not too stupid
          endif
          call msg(fu_substance_name(transport_species(i)) + ':', chemRules%mass_lagr_particle(i))
        end do
      endif  ! lagrangian defined
    
    end subroutine calculate_thresholds
    
  end subroutine set_low_mass_threshold


  !  ****************************************************************************** 
   
  subroutine point_low_mass_thresholds(chemrules, ifForward)
    ! 
    ! Point the low mass thresholds to adjoint/forward values, asseming these are already
    ! calculated.
    implicit none
    type(Tchem_rules), intent(inout) :: chemrules
    logical, intent(in) :: ifForward
    character(len=*), parameter :: subname = 'point_low_mass_thresholds'

    if (.not. (chemrules%defined == silja_true)) then
      call set_error('Chemistry rules not defined', subname)
      return
    end if

    if (ifforward) then
      if (fu_fails(associated(chemrules%low_cnc_trsh_frw), 'low_cnc_trsh_frw not associated', subname)) return
      if (fu_fails(associated(chemrules%low_mass_trsh_frw), 'low_mass_trsh_frw not associated', subname)) return
      if (fu_fails(all(chemrules%low_cnc_trsh_frw > 0), 'low_cnc_trsh_frw < 0', subname)) return
      if (fu_fails(all(chemrules%low_mass_trsh_frw > 0), 'low_mass_trsh_frw < 0', subname)) return
      chemrules%low_cnc_trsh => chemrules%low_cnc_trsh_frw
      chemrules%low_mass_trsh => chemrules%low_mass_trsh_frw
    else
      if (fu_fails(associated(chemrules%low_cnc_trsh_adj), 'low_cnc_trsh_adj not associated', subname)) return
      if (fu_fails(associated(chemrules%low_mass_trsh_adj), 'low_mass_trsh_adj not associated', subname)) return
      if (fu_fails(all(chemrules%low_cnc_trsh_adj > 0), 'low_cnc_trsh_adj < 0', subname)) return
      if (fu_fails(all(chemrules%low_mass_trsh_adj > 0), 'low_cnc_trsh_adj < 0', subname)) return
      chemrules%low_cnc_trsh => chemrules%low_cnc_trsh_adj
      chemrules%low_mass_trsh => chemrules%low_mass_trsh_adj
    end if

  end subroutine point_low_mass_thresholds

  
  !*******************************************************************

  logical function fu_low_mass_thres_defined(chemrules, ifforward) result(defined)
    implicit none
    type(Tchem_rules), intent(inout) :: chemrules
    logical, intent(in) :: ifForward
    character(len=*), parameter :: subname = 'point_low_mass_thresholds'

    if (.not. (chemrules%defined == silja_true)) then
      call set_error('Chemistry rules not defined', subname)
      return
    end if

    if (ifforward) then
      if (.not. associated(chemrules%low_mass_trsh_frw)) then
        defined = .false.
      else
        defined = chemrules%low_mass_trsh_frw(1) > 0
      end if
    else
      if (.not. associated(chemrules%low_mass_trsh_adj)) then
        defined = .false.
      else
        defined = chemrules%low_mass_trsh_adj(1) > 0
      end if
    end if
  end function fu_low_mass_thres_defined

  !*******************************************************************
  
  integer function fu_number_of_lagr_particles(iAccuracy)
    !
    ! Determines the number of particles used in the Lagrangian simulations
    !
    implicit none
    
    ! Imported parameter
    integer, intent(in) :: iAccuracy
    
    if(iAccuracy > -10 .and. iAccuracy < 20)then
      fu_number_of_lagr_particles = nint(100000 * 3.0**(iAccuracy - 5))
    else
      call msg('Strange accuracy given: ', iAccuracy)
      call set_error('Strange accuracy given: ','fu_number_of_lagr_particles')
      fu_number_of_lagr_particles = int_missing
    endif
    
  end function fu_number_of_lagr_particles


  !*******************************************************************
  subroutine  init_chemisryStuff (nSrc, nSpTr, nSpAer,nSpSl,nReact,nMetdat,nz )
    ! Just allocates necessary space
    implicit none

    ! Imported parameter
    integer, intent(in) :: nSrc, nSpTr, nSpAer,nSpSl,nReact,nMetdat,nz

    ! Local variables
    integer :: iStat, iStatP
    integer :: nthreads, ithread
    character(len=*), parameter :: subname = 'init_chemistryStuff'

    if(ChemStuff%defined == silja_true)then
      call set_error('pChemicalRunSetup structure is already allocated',subname)
      return
    endif
    nthreads = 1
    ithread = 0

    call msg("Allocating ChemStuff%arrStuff, memusage (kB)", fu_system_mem_usage())
    !
    ! Allocate chemical stuff
    !
    !$OMP PARALLEL if (if_OMP_chemistry) DEFAULT(NONE) &
    !$OMP SHARED(nthreads,iStat, ChemStuff, &
    !$OMP & nSrc, nSpTr, nSpAer,nSpSl,nReact,nMetdat,nz) &
    !$OMP & PRIVATE(ithread, istatP)

    !$OMP MASTER
    !$ nthreads = omp_get_num_threads()
      allocate(ChemStuff%arrStuff(0:nthreads-1), stat=iStat)
      if (iStat/=0) call set_error('Failed ChemStuff%arrStuff allocation', subname)
    !$OMP END MASTER
    !$OMP BARRIER
    !$ ithread = omp_get_thread_num()
    if (iStat==0) then
       allocate( &
         & ChemStuff%arrStuff(ithread)%garb_array(1:nSrc, 1:nSpTr), &
         & ChemStuff%arrStuff(ithread)%cncTrn(1:nSpTr, 1:nSrc), &
         & ChemStuff%arrStuff(ithread)%cncAer(1:nSpAer), &
         & ChemStuff%arrStuff(ithread)%cncSL(1:nSpSl), &
         & ChemStuff%arrStuff(ithread)%photorates(1:maxPhotoIndex,1:nz), &
         & ChemStuff%arrStuff(ithread)%aodext(1:nz), &
         & ChemStuff%arrStuff(ithread)%aodscat(1:nz), &
         & ChemStuff%arrStuff(ithread)%o3column(1:nz,1:nSrc) , &
         & ChemStuff%arrStuff(ithread)%pwc_col(1:nz) , &
         & ChemStuff%arrStuff(ithread)%soot_col(1:nz) , &
         & ChemStuff%arrStuff(ithread)%tau_above_bott(1:nz) , &
         & ChemStuff%arrStuff(ithread)%metdat(1:nMetdat, 1:nz), &!Should be :meteo_input%nQuantities
         & ChemStuff%arrStuff(ithread)%reactRates(1:nReact,1:nSrc), &
         & stat=iStatP)
       if (iStatP ==0) then
          ChemStuff%arrStuff(ithread)%photorates(1:maxPhotoIndex,1:nz) = F_NAN  !Should break if not initialized
          ChemStuff%arrStuff(ithread)%reactRates(1:nReact,1:nSrc) = F_NAN
       else
          !$OMP CRITICAL (BARK_ALLOC)
          call set_error("Failed ChemStuff%arrStuff("//trim(fu_str(iThread))//") allocation", subname)
          !$OMP END CRITICAL (BARK_ALLOC)
       endif
    endif
    !$OMP END PARALLEL
    if (.not. error) ChemStuff%defined = silja_true
    call msg("Allocating ChemStuff%arrStuff done for numthreads, memusage (kB)", nthreads, fu_system_mem_usage())
    
  end subroutine init_chemisryStuff

  !*************************************************************************************************

  subroutine init_chemicalRunSetup(nSpeciesSet, ChemicalRunSetup)
    !
    ! Just allocates necessary space
    !
    implicit none

    ! Imported parameter
    integer, intent(in) :: nSpeciesSet
    type(TChemicalRunSetup), intent(out) :: ChemicalRunSetup

    ! Local variables
    integer :: iStat

    if(ChemicalRunSetup%defined == silja_true)then
      call set_error('pChemicalRunSetup structure is already allocated','init_chemicalRunSetup')
      return
    endif

!    allocate(ChemicalRunSetup%iSpeciesTypes(nSpeciesSet), &
!           & ChemicalRunSetup%indChemical(nSpeciesSet), &
!           & ChemicalRunSetup%indEmisSpecies(nSpeciesSet), &
!           & ChemicalRunSetup%indTranspSpecies(nSpeciesSet), &
!           & ChemicalRunSetup%indOptDnsSpecies(nSpeciesSet), stat=iStat)
!    if(iStat /= 0)then
!      call set_error('Failed to allocate pChemicalRunSetup','init_chemicalRunSetup')
!      return
!    endif
!
!    ChemicalRunSetup%iSpeciesTypes(1:nSpeciesSet) = int_missing
!    ChemicalRunSetup%indChemical(1:nSpeciesSet) = int_missing
!    ChemicalRunSetup%indEmisSpecies  (1:nSpeciesSet) = int_missing
!    ChemicalRunSetup%indTranspSpecies(1:nSpeciesSet) = int_missing
!    ChemicalRunSetup%indOptDnsSpecies(1:nSpeciesSet) = int_missing

    ChemicalRunSetup%defined = silja_true
    
  end subroutine init_chemicalRunSetup

  
  !*******************************************************************************
  
  subroutine data_for_photolysis_and_cld_model(mapTransport, chemRules, indO3, ind_soot, &
                                             & frac_soot, n_soot, density, diam)
    !
    ! Catches the index of ozone in the transport species and counts soot-related aerosols
    !
    implicit none

    ! Imported parameters
    type(Tmass_map), intent(in) :: mapTransport
    type(Tchem_rules), intent(in) :: chemRules
    integer, intent(out) :: indO3
    integer, dimension(max_species), intent(out) :: ind_soot
    real, dimension(max_species), intent(out) :: frac_soot, density, diam
    integer, intent(out) :: n_soot

    ! Local parameter
    logical, parameter :: if_soot_affects_photolysis = .true.
    ! Local variables
    real, dimension(:), pointer :: fractions_opt_subst
    character(len=substNmLen), dimension(max_species) :: names_opt_subst
    integer :: iSpecies, iSubst, is, n_opt_subst
    character(len=*), parameter :: sub_name='data_for_photolysis_and_cld_model'
    
    ! ozone needed?
    if (chemRules%need_photo_lut) then
      do iSpecies = 1, mapTransport%nSpecies
        if(trim(fu_str_u_case(fu_name(fu_material(mapTransport%species(iSpecies))))) == 'O3')then
          indO3 = iSpecies
          exit
        endif
      end do
      if(fu_fails(indO3 /= int_missing,'Undefined ozone index: Needed for photolysis LUT', sub_name))return
    else
      indO3 = int_missing
    end if
    !
    ! Store aerosol particle features
    !
    diam(:) = 0.0
    density(:) = real_missing

    do iSpecies = 1, mapTransport%nSpecies
      diam(iSpecies) = fu_nominal_d(mapTransport%species(iSpecies))
      density(iSpecies) = fu_density(mapTransport%species(iSpecies))
    end do
    !
    ! Find soot-referring particles (those, which have soot in referece optical species)
    !
    ind_soot(:) = int_missing
    frac_soot(:) = 0.0
    n_soot = 0
    if (if_soot_affects_photolysis) then
      fractions_opt_subst => fu_work_array()
      is = 1
      do iSpecies = 1, mapTransport%nSpecies
        n_opt_subst = fu_n_opt_subst(fu_material(mapTransport%species(iSpecies)))
        if (n_opt_subst > 0) then
          names_opt_subst(1:n_opt_subst) = fu_names_opt_subst(fu_material(mapTransport%species(iSpecies)))
          fractions_opt_subst(1:n_opt_subst) = fu_fractions_opt_subst(fu_material(mapTransport%species(iSpecies)))
          do iSubst = 1, n_opt_subst
            if (trim(names_opt_subst(iSubst)) == 'soot') then
              ind_soot(is) = iSpecies
              frac_soot(is) = fractions_opt_subst(iSubst)
              is = is + 1
            end if
          end do  ! n_opt_subst
        end if  ! n_opt_subst > 0
      end do  ! transport species
      call free_work_array(fractions_opt_subst)
      n_soot = is - 1
    end if   ! if_soot_affects_photolysis

  end subroutine data_for_photolysis_and_cld_model


  !=====================================================================================

  function fu_low_mass_thresholds(chemRules) result(thrsh)
    !
    ! Returns the low-mass threshold stored in some specific chemical rules
    !
    implicit none
    type(Tchem_rules), intent(in) :: chemRules
    real, dimension(:), pointer :: thrsh
    
    if (associated(chemrules%low_mass_trsh)) then
       thrsh => chemrules%low_mass_trsh
    else
       call msg_warning("Low mass threshold is not set, returning null!", "fu_low_mass_thresholds")
       thrsh => null()
    endif

  end function fu_low_mass_thresholds

  !=====================================================================================

  function fu_mass_lagrangian_particle(chemRules) result(mass)
    !
    ! Returns the low-mass threshold stored in some specific chemical rules
    !
    implicit none
    type(Tchem_rules), intent(in) :: chemRules
    real, dimension(:), pointer :: mass

    mass => chemrules%mass_lagr_particle

  end function fu_mass_lagrangian_particle


  !=====================================================================================

  function fu_aerosol_rules(chemRules) result(rulesAerosol)
    ! Encapsulates the optical-density computation rules
    implicit none
    type(Taerosol_rules), pointer :: rulesAerosol
    type(Tchem_rules), intent(in), target :: chemRules
    rulesAerosol => chemRules%rulesAerosol
  end function fu_aerosol_rules


  !======================================================================================  

  function fu_deposition_rules(chemrules) result(deprules)
    implicit none
    type(Tchem_rules), intent(in) :: chemrules
    
    type(Tdeposition_rules), pointer :: deprules

    deprules => chemrules%rulesDeposition
  end function fu_deposition_rules


  !=====================================================================================

  function fu_chemRunSetup(chemrules) result(setup)
    type(Tchem_rules), intent(in) :: chemrules
    type(TchemicalRunSetup), pointer :: setup

    setup => chemrules%chemrunsetup
  end function fu_chemRunSetup


  !=====================================================================================

  function fu_optical_rules(chemrules) result(rules)
    implicit none
    type(Tchem_rules), intent(in) :: chemrules

    type(Toptical_density_rules), pointer :: rules

    rules => chemrules%rulesOpticDens
  end function fu_optical_rules


  !***************************************************************************
  !
  !  PRIVATE-only functions of this module
  !
  !***************************************************************************

  !***************************************************************************

  real function fu_compute_thresh_for_release(total_amt_src, &
                                            & nEmittingDispGridCells, &
                                            & residence_time_adv, computed_period, timestep, &
                                            & nx, ny, nz, &
                                            & dx_dy, dz, &
                                            & iComputationAccuracy) result(massLowThreshold)
    implicit none
    real, intent(in) :: total_amt_src, dx_dy, dz
    integer, intent(in) :: nx, ny, nz, nEmittingDispGridCells, iComputationAccuracy
    type(silja_interval), intent(in) :: computed_period, residence_time_adv, timestep
!    type(Tsilam_namelist), pointer :: nlChemSetup
        
    ! Local parameters
    real, parameter :: fNonUniformity = 2.0e-5

    ! Local variables
    type(silja_interval) :: lifetime
    real :: fMassUniform, fEmisMassQuantum
    integer :: iType

!    !
!    ! A bit of checking and distributing the tasks to species-specific rules
!    !
!    if(.not. associated(nlChemSetup))then
!      call set_error('Not associated chemistry namelist','fu_compute_thresh_for_release')
!      return
!    endif
    !
    ! iComputationAccuracy varies from 0 to 10, neutral 5. Then, let's make the threshold
    ! accuracy varying roughly +- 2 orders of magnitude
    !
    call msg('Applying standard threshold with accuracy factor:', 3**real(5-iComputationAccuracy))
    !
    ! Get the lifetime due to chemistry; for dynamically-computed emissions
    ! we use the "typical" concentrations rather then the total amnt from the emission
    ! source
    !
    if (residence_time_adv > computed_period) then
            lifetime = computed_period
    else
            lifetime = residence_time_adv
    endif
    !
    ! Mean mass in a grid cell (uniform distribution). That concludes the step 1. Note extremely 
    ! crude assumptions: uniform grid and uniform mass distribution along the vertical.
    ! The second limitation is the mean mass amount emitted per tmie step per cell, for which 
    ! the emission into 2 layers is assumed
    !
    fMassUniform = total_amt_src * (lifetime/computed_period) / (nx*ny*nz*dz*dx_dy)
    fEmisMassQuantum = total_amt_src * (timestep/computed_period) / (nEmittingDispGridCells*2)
    !
    ! The mass threshold is the min of these two
    !
    massLowThreshold = min(fEmisMassQuantum, fMassUniform) * fNonUniformity * &
                     & 3**real(5-iComputationAccuracy)

  end function fu_compute_thresh_for_release


  !*****************************************************************************************

  real function fu_lagr_part_mass_for_release(total_amt_src, &
                                       & nEmittingDispGridCells, &
                                       & residence_time_adv, computed_period, timestep, &
                                       & nLPs) result(lpMass)
    implicit none
    
    real, intent(in) :: total_amt_src
    integer, intent(in) :: nEmittingDispGridCells, nLPs
    type(silja_interval), intent(in) :: computed_period, residence_time_adv, timestep

    ! Local variables
    type(silja_interval) :: lifetime
    real :: fMassUniform, fEmisMassQuantum
    integer :: iType

!    !
!    ! A bit of checking and distributing the tasks to species-specific rules
!    !
!    if(.not. associated(nlChemSetup))then
!      call set_error('Not associated chemistry namelist','fu_compute_thresh_for_release')
!      return
!    endif
    !
    ! iComputationAccuracy varies from 0 to 10, neutral 5. Then, let's make the threshold
    ! accuracy varying roughly +- 2 orders of magnitude
    !
    call msg('Applying standard threshold for this number of LPs:', nLPs)
    !
    ! Get the lifetime due to chemistry; for dynamically-computed emissions
    ! we use the "typical" concentrations rather then the total amnt from the emission
    ! source
    !
!    lifetime = residence_time_adv
    !
    ! Mean mass in a grid cell (uniform distribution). That concludes the step 1. Note extremely 
    ! crude assumptions: uniform grid and uniform mass distribution along the vertical.
    ! The second limitation is the mean mass amount emitted per tmie step per cell, for which 
    ! the emission into 2 layers is assumed
    !
!    fEmisMassQuantum = total_amt_src * (timestep/computed_period) / nEmittingDispGridCells
    !
    ! The mass threshold is the min of these two
    !
!    lpMass = fEmisMassQuantum / nLPs * (residence_time_adv/computed_period)

    lpMass = 1.01 * total_amt_src / nLPs * min(1., (residence_time_adv/computed_period))

  end function fu_lagr_part_mass_for_release

end module chemistry_manager
