
MODULE nuclides

  ! Description:
  ! Module 'nuclides' contains the implementation of a radionuclide database
  ! for SILAM. The properties defined in this module should be such that their
  ! only use is in the context of radioactive decay and dose assessment. Other
  ! properties should be defined in module 'materials' (general properties) or
  ! module 'chemicals' (properties needed for calculation of chemical
  ! reactions). Modules 'nuclides' and 'chemicals' are used by module
  ! 'materials'.
  ! Note that the information about nuclides is in the main material database
  ! whereas here we only set pointers to these nuclides. This pointer-wise approach
  ! allows for an easy treatment of decay chains: the pointers are ordered in
  ! a proper way.
  ! In this module, a data structure (type) suitable for the representation of
  ! one nuclide's properties is defined. Thereafter, nuclides are defined as
  ! public constants of this type for easy individual reference. For sequential
  ! reference, a pointer array is defined and its elements are set to point to
  ! the individual nuclides. Functions are provided e.g. to retrieve nuclides
  ! by their name, values by nuclide and property name, and units by property
  ! name.
  !
  ! Original code: Mikko Ilvonen, VTT/ENE, email mikko.ilvonen@vtt.fi
  ! Corrected and extended by M.Sofiev, FMI, 2004
  ! New code for SILAM v.5: M.Sofiev, 2010
  !
  ! All units: SI
  !
  ! Language: ANSI standard Fortran 90
  !
  ! Modules used:
  !
  use silam_namelist
  use silam_times

  IMPLICIT NONE

  ! The public functions and subroutines available in this module:

!  PUBLIC initialize_nuclides
!  public fu_if_nuclides_initialized
  public set_nuclide_from_namelist
  public fu_mole_to_bq

  ! Public encapsulation functions for the silam_nuclide
  PUBLIC fu_nuc_name
  PUBLIC fu_half_life
  public fu_n_daughters
  public fu_daughter
  public fu_decay_branch_fraction
  PUBLIC fu_cloud_dc
  PUBLIC fu_fallout_dc
  PUBLIC fu_if_noble_gas
  public fu_decay_enegries
  public report
  public defined

  ! The private functions and subroutines not to be used elsewhere:
!  public read_nuclide_new
!  private set_nuclide_from_namelist
  private fu_mole_to_bq_mole_amount
  private fu_mole_to_bq_nuclide
  private report_nuclide_data

  ! Encapsulation functions by name
  PRIVATE fu_if_noble_gas_by_name

  ! Encapsulation functions by nuclide reference
  PRIVATE fu_nuc_name_of_nuclide
  PRIVATE fu_half_life_of_nuclide
  PRIVATE fu_cloud_dc_of_nuclide
  PRIVATE fu_fallout_dc_of_nuclide
  PRIVATE fu_if_noble_gas_of_nuclide
  private set_missing_nuclide
  private fu_nuclide_defined
  private copy_nuclide


  ! Public encapsulation functions for the silam_nuclide are all interfaced
  ! to the by_index, by_name and of_nuclide functions. The idea that nuclides
  ! are stored in the PUBLIC array nuc_array, so that each nuclide can be
  ! addressed from there either through its index in that array, or
  ! directly via silam_nuclide pointer. The first way is kept only for
  ! compatibility with the historical code. Correct way is, of course, a
  ! pointer addressing. Direct indexing is VERY DANGERROUS and should not
  ! be used in newly created routines. By_name addressing is done just
  ! for the conveniency - it is quite slow because the name has to be
  ! searched in the array of ~500 elements.
  !
  interface fu_nuc_name
    module procedure fu_nuc_name_of_nuclide
  end interface

  interface report
    module procedure report_nuclide_data
  end interface

  interface fu_half_life
    module procedure fu_half_life_of_nuclide
  end interface

  interface fu_cloud_dc
    module procedure fu_cloud_dc_of_nuclide
  end interface

  interface fu_fallout_dc
    module procedure fu_fallout_dc_of_nuclide
  end interface

  interface fu_if_noble_gas
    module procedure fu_if_noble_gas_by_name
    module procedure fu_if_noble_gas_of_nuclide
  end interface

  interface fu_mole_to_bq
    module procedure fu_mole_to_bq_mole_amount
    module procedure fu_mole_to_bq_nuclide
  end interface
  
  interface defined
    module procedure fu_nuclide_defined
  end interface
  
  interface assignment (=)
     module procedure copy_nuclide
  end interface

  ! References:
  !
  ! The dose-rate conversion factors for external exposure are based on the
  ! following sources of information:
  !
  ! D. C. Kocher, Dose-rate conversion factors for external exposure to photon
  ! and electron radiation from radionuclides occurring in routine releases
  ! from nuclear fuel cycle facilities, Health Physics, Vol. 38 (April),
  ! 543-621, 1980
  !
  ! D. C. Kocher, Dose-rate conversion factors for external exposure to photons
  ! and electrons, Health Physics, Vol. 45 (September), 665-686, 1983
  !
  ! RSIC (Radiation Shielding Information Center), DOSFACTER II: Calculation of
  ! dose-rate conversion factors for exposure to photons and electrons, RSIC
  ! computer code collection, CCC-400, Oak Ridge National Laboratory, Oak
  ! Ridge, Tennessee, 1981
  !
  ! RSIC (Radiation Shielding Information Center), DOSFACTER-DOE: Dose-rate
  ! conversion factors for external exposure to photons and electrons, RSIC
  ! computer code collection, CCC-536, Oak Ridge National Laboratory, Oak
  ! Ridge, Tennessee, 1988
  !
  ! RSIC (Radiation Shielding Information Center), DOSDAT-DOE: Dose-rate
  ! conversion factors for external exposure to photons and electrons, RSIC
  ! data library collection, DLC-144, Oak Ridge National Laboratory, Oak
  ! Ridge, Tennessee, 1989

  ! Integer constants to facilitate the handling of structured nuclide data:

  INTEGER, PARAMETER :: max1_external_dc   = 4
  INTEGER, PARAMETER :: max1_inhalation_dc = 7
  INTEGER, PARAMETER :: max1_ingestion_dc  = 6  ! if this changes, code too !!

  ! The exposure modes for gamma radiation are immersion in contaminated air,
  ! immersion in contaminated water, and standing on contaminated ground
  ! surface.

  INTEGER, PARAMETER :: no_of_exp_modes = 3
  INTEGER, PARAMETER :: air    = 1
  INTEGER, PARAMETER :: water  = 2
  INTEGER, PARAMETER :: ground = 3

  ! The target organs for gamma radiation are adrenals, bladder, brain, breast,
  ! heart, small intestine, upper large intestine, lower large intestine,
  ! kidneys, liver, lungs, marrow, red marrow, ovaries, pancreas, skeleton,
  ! spleen, stomach, testes, thymus, thyroid, uterus and skin. In addition, the
  ! factor for effective dose is given as the last 'organ'.

  INTEGER, PARAMETER :: no_of_organs_ext = 24
  INTEGER, PARAMETER :: adrenals     = 1
  INTEGER, PARAMETER :: bladder      = 2
  INTEGER, PARAMETER :: brain        = 3
  INTEGER, PARAMETER :: breast       = 4
  INTEGER, PARAMETER :: heart        = 5
  INTEGER, PARAMETER :: s_intestine  = 6
  INTEGER, PARAMETER :: ul_intestine = 7
  INTEGER, PARAMETER :: ll_intestine = 8
  INTEGER, PARAMETER :: kidneys      = 9
  INTEGER, PARAMETER :: liver        = 10
  INTEGER, PARAMETER :: lungs        = 11
  INTEGER, PARAMETER :: marrow       = 12
  INTEGER, PARAMETER :: r_marrow     = 13
  INTEGER, PARAMETER :: ovaries      = 14
  INTEGER, PARAMETER :: pancreas     = 15
  INTEGER, PARAMETER :: skeleton     = 16
  INTEGER, PARAMETER :: spleen       = 17
  INTEGER, PARAMETER :: stomach      = 18
  INTEGER, PARAMETER :: testes       = 19
  INTEGER, PARAMETER :: thymus       = 20
  INTEGER, PARAMETER :: thyroid      = 21
  INTEGER, PARAMETER :: uterus       = 22
  INTEGER, PARAMETER :: skin         = 23
  INTEGER, PARAMETER :: effective    = 24

  ! Factors for beta radiation are given for four depths in skin.

  INTEGER, PARAMETER :: no_of_depths = 4
  INTEGER, PARAMETER :: depth1 = 1
  INTEGER, PARAMETER :: depth2 = 2
  INTEGER, PARAMETER :: depth3 = 3
  INTEGER, PARAMETER :: depth4 = 4

  INTEGER, PARAMETER :: max1_internal_dc        = 2
  INTEGER, PARAMETER :: no_of_integration_times = 9
  INTEGER, PARAMETER :: no_of_organs_int        = 10

  ! Public types with private components defined in this module:

  TYPE silja_branching  ! the coding of all daughters for one single nuclide
    INTEGER :: n_daughters
    CHARACTER (LEN = nuc_name_len), DIMENSION (max_daughters) :: daughter
    REAL, DIMENSION (max_daughters) :: branch_fraction  ! ~ probability
!    INTEGER :: chain_start  ! length of chain at progenitress, otherwise -1
  END TYPE silja_branching

  ! When using type 'silam_nuclide', it is meant that the pointer fields
  ! relative1 and relative2 are set to point to the right nuclides by
  ! 'initialize_nuclides', based on the nuclide names returned by
  ! 'read_nuclide'. The purpose of these pointers is to speed up chain decay
  ! calculation.

  TYPE silam_nuclide
    PRIVATE

    TYPE (silja_logical)  :: defined
    integer               :: atomic_no
    logical               :: metastable
    CHARACTER (LEN = nuc_name_len) :: name
    type(silja_interval)           :: half_life

    ! The following fields contain the information needed for the old ARANO
    ! style chain decay scheme:

!    CHARACTER (LEN=3)                       :: chain_type
!    CHARACTER (LEN=1)                       :: role_in_chain
!    TYPE (silam_nuclide), POINTER           :: relative1
!    TYPE (silam_nuclide), POINTER           :: relative2
!    character(len=nuc_name_len) :: relative1, relative2
!    REAL (KIND = SELECTED_REAL_KIND (4, 0)) :: probability1
!    REAL (KIND = SELECTED_REAL_KIND (4, 0)) :: probability2

    ! The following field contains the information for the new coding system of
    ! arbitrary decay chains:

    TYPE (silja_branching) :: branching

    REAL(r4k), DIMENSION (max1_external_dc)   :: external_dc
    REAL(r4k), DIMENSION (max1_inhalation_dc) :: inhalation_dc
    REAL(r4k), DIMENSION (max1_ingestion_dc)  :: ingestion_dc

    REAL(r4k), DIMENSION (no_of_exp_modes, no_of_organs_ext) :: ext_gamma_dc
    REAL(r4k), DIMENSION (no_of_exp_modes, no_of_depths)     :: ext_beta_dc

  ! In the next two arrays, first index means: 1 = inhalation, 2 = ingestion

    REAL(r4k), DIMENSION (max1_internal_dc, no_of_integration_times)  :: int_eff_dc
    REAL(r4k), DIMENSION (max1_internal_dc, no_of_organs_int)         :: int_organ_dc
    
    ! Energy of the decay and its probability (2, nEnergies)
    real(r4k), dimension(2,max_nuclide_energies) :: pDecayEnergy

  END TYPE silam_nuclide

  ! Because of pointers, 'nuclide_missing' can't be initialized. Instead of a
  ! parameter, a function should be used. (Until that, this is a variable...)
  !
  TYPE (silam_nuclide), PUBLIC, SAVE :: nuclide_missing

  type silam_nuclide_ptr
    TYPE (silam_nuclide), POINTER :: ptr
  end type silam_nuclide_ptr
  public silam_nuclide_ptr


CONTAINS

  ! ***************************************************************

  subroutine set_missing_nuclide()
    !
    ! Sets a nuclide_missing variable to truly missing value, nullifies pointers
    !
    ! Author: Mikhail Sofiev, FMI
    !
    implicit none

    nuclide_missing%defined = silja_false
    nuclide_missing%atomic_no = 0
    nuclide_missing%metastable = .false.
    nuclide_missing%name = ''
    nuclide_missing%half_life = interval_missing
!    nuclide_missing%chain_type = ''
!    nuclide_missing%role_in_chain = ''

!    nullify(nuclide_missing%relative1)
!    nullify(nuclide_missing%relative2)
!    nullify(nuclide_missing%relative1)
!    nullify(nuclide_missing%relative2)
!
!    nuclide_missing%probability1 = 0.
!    nuclide_missing%probability2 = 0.

    nuclide_missing%branching%n_daughters =0 ! This is enough to remove branching

    nuclide_missing%external_dc = 0.
    nuclide_missing%inhalation_dc = 0.
    nuclide_missing%ingestion_dc = 0.

    nuclide_missing%ext_gamma_dc = 0.
    nuclide_missing%ext_beta_dc = 0.

    nuclide_missing%int_eff_dc = 0.
    nuclide_missing%int_organ_dc = 0.
    nuclide_missing%pDecayEnergy = 0.0

  end subroutine set_missing_nuclide


!******************************************************************************

  real function fu_mole_to_Bq_mole_amount(amount_mole, half_life_sec) result(amount_Bq)
    !
    ! Just converts the amount of the stuff from amount of moles to becquerrels
    !
    ! Imported parameters with the intent IN
    real, intent(in) :: amount_mole, half_life_sec ![mole],[sec] respectively

    amount_Bq = ln_2 * avogadro * amount_mole / half_life_sec

    !    A = lambda * N, or
    !activity = decay rate constant (1/s) * number of atoms,
    !and lambda * half_life = ln (2)

  end function fu_mole_to_Bq_mole_amount


!******************************************************************************

  real function fu_mole_to_Bq_nuclide(nuclide, amount_mole) result(amount_Bq)
    !
    ! Just converts the amount of the stuff from amount of moles to becquerrels
    ! Half-life is taken from the nuclide ID
    !
    ! Imported parameters with the intent IN
    real, intent(in) :: amount_mole  ![mole]
    type(silam_nuclide), intent(in) :: nuclide
    
    if(defined(nuclide%half_life))then
      amount_Bq = ln_2 * avogadro * amount_mole / fu_sec(nuclide%half_life)
    else
      call msg('Nuclide:' + nuclide%name + ', half-life:' + fu_interval_to_named_string(nuclide%half_life))
      call set_error('Nuclide with strange half-life time','fu_mole_to_Bq_nuclide')
      amount_Bq = 0.0
    end if
    
    !    A = lambda * N, or
    !activity = decay rate constant (1/s) * number of atoms,
    !and lambda * half_life = ln (2)

  end function fu_mole_to_Bq_nuclide


  !******************************************************************

  subroutine set_nuclide_from_namelist(nlN, pN) !, nuc1, nuc2, nuc3)
    !
    ! Writes the parameters of the gaseous deposition parameterization
    !
    implicit none
    
    ! Imported parameters
    type(Tsilam_namelist), intent(in) :: nlN
    type(silam_nuclide), intent(out) :: pN
!    character(len=*), intent(out) :: nuc1, nuc2, nuc3

    ! Local variables
    type(silam_sp) :: sp
    integer :: iB, iMode, iItem, nItems
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: pItems
    
    sp%sp => fu_work_string()
    if(error)return
    !
    ! Basic information
    !
    pN%name = fu_content(nlN,'nuclide_name')
    pN%atomic_no = fu_content_int(nlN,'nuclide_atomic_nbr')
!    pN%mass_no = fu_content_int(nlN,'nuclide_mass_nbr')
    pN%metastable = (fu_str_u_case(fu_content(nlN,'nuclide_metastable')) == 'YES')
    pN%half_life = fu_set_named_interval(fu_content(nlN,'half_life_period'))
    if(error)return
    !
    ! Decay branching
    !
    pN%branching%n_daughters = fu_content_int(nlN,'nuclide_branching_nbr_of_daughters')
!    pN%branching%chain_start = fu_content_int(nlN,'nuclide_branching_chain_start')
    do iB = 1, pN%branching%n_daughters
      write(unit=sp%sp,fmt='(A27,I1)')'nuclide_branching_daughter_',iB
      pN%branching%daughter(iB) = fu_content(nlN,sp%sp)
      write(unit=sp%sp,fmt='(A36,I1)')'nuclide_branching_daughter_fraction_',iB
      pN%branching%branch_fraction(iB) = fu_content_real(nlN,sp%sp)
    end do
    !
    ! Doses
    !
    sp%sp = fu_content(nlN,'nuclide_external_dose_rates')
    read(unit=sp%sp, fmt=*) (pN%external_dc(iB),iB=1,max1_external_dc)

    sp%sp = fu_content(nlN,'nuclide_inhalation_dose_rates')
    read(unit=sp%sp, fmt=*) (pN%inhalation_dc(iB),iB=1,max1_inhalation_dc)

    sp%sp = fu_content(nlN,'nuclide_ingestion_dose_rates')
    read(unit=sp%sp, fmt=*) (pN%ingestion_dc(iB),iB=1,max1_ingestion_dc)

    nullify(pItems)
    call get_items(nlN,'nuclide_external_gamma_dose_rates_for_mode',pItems,nItems)
    if(nItems /= no_of_exp_modes)then
      call msg('Number of nuclide_external_gamma_dose_rates_for_mode lines is not equal to number of modes:', &
             & nItems, no_of_exp_modes)
      call set_error('Strange number of nuclide_external_gamma_dose_rates_for_mode lines', &
                   & 'set_nuclide_from_namelist')
      return
    endif
    do iItem = 1, nItems
      iMode = iItem       ! just to be on a safe side
      sp%sp = fu_content(pItems(iMode))
      read(unit=sp%sp, fmt=*) iMode, (pN%ext_gamma_dc(iMode,iB), iB = 1, no_of_organs_ext)
    end do
    
    call get_items(nlN,'nuclide_external_beta_dose_rates_for_mode',pItems,nItems)
    if(nItems /= no_of_exp_modes)then
      call msg('Number of nuclide_external_beta_dose_rates_for_mode lines is not equal to number of modes:', &
             & nItems, no_of_exp_modes)
      call set_error('Strange nuclide_external_beta_dose_rates_for_mode of nuclide_external_gamma_dose_rates lines', &
                   & 'set_nuclide_from_namelist')
      return
    endif
    do iItem = 1, nItems
      iMode = iItem       ! just to be on a safe side
      sp%sp = fu_content(pItems(iMode))
      read(unit=sp%sp, fmt=*) iMode, (pN%ext_beta_dc(iMode,iB), iB = 1, no_of_depths)
    end do
    
    call destroy_items(pItems)
    
    ! In the first index, 1 = inhalation, 2 = ingestion

    sp%sp = fu_content(nlN,'nuclide_internal_effective_inhalation_dose')
    read(unit=sp%sp, fmt=*)(pN%int_eff_dc(1,iB),iB=1,no_of_integration_times)
    
    sp%sp = fu_content(nlN,'nuclide_internal_effective_ingestion_dose')
    read(unit=sp%sp, fmt=*)(pN%int_eff_dc(2,iB),iB=1,no_of_integration_times)
    
    sp%sp = fu_content(nlN,'nuclide_internal_per_organ_inhalation_dose')
    read(unit=sp%sp, fmt=*)(pN%int_organ_dc(1,iB),iB=1,no_of_organs_int)
    
    sp%sp = fu_content(nlN,'nuclide_internal_per_organ_ingestion_dose')
    read(unit=sp%sp, fmt=*)(pN%int_organ_dc(2,iB),iB=1,no_of_organs_int)

    ! Have to read the decay energy and probability of it - as many as entered. 
    ! Format: <energy> <probability> [<energy> <probability> [...]]
    ! If unknown, the first couple is -1 -1. 
    !
    sp%sp = fu_content(nlN,'energy_and_decay_prob')
    nItems = fu_nbrOfWords(sp%sp)   ! how many items?
    if (nItems/2 > max_nuclide_energies) then
       call set_error('max number of energies set too low -> just set it larger','set_nuclide_from_namelist')
       !This error will happen if the nuclides.dat file is updated and
       !some nuclide gets larger list of energies and probabilities
    end if
    pN%pDecayEnergy = 0.0 !Set all to zero, then fill the actual values
    read(unit=sp%sp, fmt=*)((pN%pDecayEnergy(iMode,iItem),iMode=1,2),iItem=1,nItems/2)
    
    pN%defined = silja_true

    call free_work_array(sp%sp)
    !call report_nuclide_data(pN)

  end subroutine set_nuclide_from_namelist


  !******************************************************************

  subroutine report_nuclide_data(pN)
    !
    ! Writes the parameters of the gaseous deposition parameterization
    !
    implicit none
    
    ! Imported parameters
    type(silam_nuclide), pointer :: pN

    ! Local variables
    type(silam_sp) :: sp
    integer :: iB, iMode, iItem, nItems, i
    
    sp%sp => fu_work_string()
    if(error)return
    !
    ! Basic information
    !
    call write_namelist_item('nuclide_name',pN%name)
    iB = pN%atomic_no
    call write_namelist_item('nuclide_atomic_nbr',iB)
    if(pN%metastable)then
      call write_namelist_item('nuclide_metastable','YES')
    else
      call write_namelist_item('nuclide_metastable','NO')
    endif
    call write_namelist_item('half_life_period',fu_interval_to_named_string(pN%half_life))
    sp%sp = fu_interval_to_named_string(pN%half_life)
    if(index(sp%sp,'*') > 0)then
      call msg('')
    endif
    !
    ! Decay branching
    !
    call write_namelist_item('nuclide_branching_nbr_of_daughters',pN%branching%n_daughters)
    do iB = 1, pN%branching%n_daughters
      write(unit=sp%sp,fmt='(A27,I1)')'nuclide_branching_daughter_',iB
      call write_namelist_item(sp%sp,pN%branching%daughter(iB))
      write(unit=sp%sp,fmt='(A36,I1)')'nuclide_branching_daughter_fraction_',iB
      call write_namelist_item(sp%sp,pN%branching%branch_fraction(iB))
    end do
    !
    ! Doses
    !
    write(unit=sp%sp, fmt='(100(E9.3,1x))') (pN%external_dc(iB),iB=1,max1_external_dc)
    if(index(sp%sp,'*') > 0)then
      call msg('')
    endif

    call write_namelist_item('nuclide_external_dose_rates',sp%sp)

    write(unit=sp%sp, fmt='(100(E9.3,1x))') (pN%inhalation_dc(iB),iB=1,max1_inhalation_dc)
    call write_namelist_item('nuclide_inhalation_dose_rates',sp%sp)

    write(unit=sp%sp, fmt='(100(E9.3,1x))') (pN%ingestion_dc(iB),iB=1,max1_ingestion_dc)
    call write_namelist_item('nuclide_ingestion_dose_rates',sp%sp)

    call write_namelist_item('nuclide_external_gamma_dose_number_of_modes',no_of_exp_modes)
    do iMode = 1, no_of_exp_modes
      write(unit=sp%sp, fmt='(I2,1x,100(E9.3,1x))') iMode, (pN%ext_gamma_dc(iMode,iB), iB = 1, no_of_organs_ext)
      call write_namelist_item('nuclide_external_gamma_dose_rates_for_mode',sp%sp)
    end do

    call write_namelist_item('nuclide_external_beta_dose_number_of_modes',no_of_exp_modes)
    do iMode = 1, no_of_exp_modes
      write(unit=sp%sp, fmt='(I2,1x,100(E9.3,1x))') iMode, (pN%ext_beta_dc(iMode,iB), iB = 1, no_of_depths)
      call write_namelist_item('nuclide_external_beta_dose_rates_for_mode',sp%sp)
    enddo

    ! In the first index, 1 = inhalation, 2 = ingestion

    write(unit=sp%sp, fmt='(100(E9.3,1x))')(pN%int_eff_dc(1,iB),iB=1,no_of_integration_times)
    call write_namelist_item('nuclide_internal_effective_inhalation_dose',sp%sp)
    
    write(unit=sp%sp, fmt='(100(E9.3,1x))')(pN%int_eff_dc(2,iB),iB=1,no_of_integration_times)
    call write_namelist_item('nuclide_internal_effective_ingestion_dose',sp%sp)
    
    write(unit=sp%sp, fmt='(100(E9.3,1x))')(pN%int_organ_dc(1,iB),iB=1,no_of_organs_int)
    call write_namelist_item('nuclide_internal_per_organ_inhalation_dose',sp%sp)
    
    write(unit=sp%sp, fmt='(100(E9.3,1x))')(pN%int_organ_dc(2,iB),iB=1,no_of_organs_int)
    call write_namelist_item('nuclide_internal_per_organ_ingestion_dose',sp%sp)
    
    !Energies and probabilities
    do i = 1,size(pN%pDecayEnergy(1,:))
       if (pN%pDecayEnergy(1,i) == 0.0) then
          nItems = size(pN%pDecayEnergy(1,:i-1))*2
          exit
       end if
    end do
    write(unit=sp%sp, fmt='(100(E9.3,1x))')((pN%pDecayEnergy(iMode,iItem),iMode=1,2),iItem=1,nItems/2)
    call write_namelist_item('energy_and_decay_prob',sp%sp)
       
    
    call free_work_array(sp%sp)

  end subroutine report_nuclide_data


  !***********************************************************************

  logical function fu_nuclide_defined(nuc)
    implicit none
    
    ! Imported parameters
    type(silam_nuclide), intent(in) :: nuc

    fu_nuclide_defined = (nuc%defined == silja_true)

  end function fu_nuclide_defined


  !******************************************************************************

  logical FUNCTION fu_if_noble_gas_by_name (name)
    CHARACTER (LEN = *), INTENT (in) :: name
    IF (name (1 : 2) == 'H_' .or. &
        name (1 : 2) == 'HE' .or. &
        name (1 : 2) == 'N_' .or. &
        name (1 : 2) == 'O_' .or. &
        name (1 : 2) == 'F_' .or. &
        name (1 : 2) == 'NE' .or. &
        name (1 : 2) == 'CL' .or. &
        name (1 : 2) == 'AR' .or. &
        name (1 : 2) == 'KR' .or. &
        name (1 : 2) == 'XE' .or. &
        name (1 : 2) == 'RN') THEN
      fu_if_noble_gas_by_name = .true.
    ELSE
      fu_if_noble_gas_by_name = .false.
    END IF
  END FUNCTION fu_if_noble_gas_by_name


  ! ***************************************************************

  FUNCTION fu_nuc_name_of_nuclide (nuclide)
    CHARACTER (LEN = nuc_name_len) :: fu_nuc_name_of_nuclide
    TYPE(silam_nuclide), INTENT (in) :: nuclide
    fu_nuc_name_of_nuclide = nuclide%name
  END FUNCTION fu_nuc_name_of_nuclide

  REAL FUNCTION fu_half_life_of_nuclide (nuclide)
    type(silam_nuclide), INTENT (in) :: nuclide
    fu_half_life_of_nuclide = fu_sec(nuclide%half_life)
  END FUNCTION fu_half_life_of_nuclide

  REAL FUNCTION fu_cloud_dc_of_nuclide (nuclide)
    type(silam_nuclide), INTENT (in) :: nuclide
    fu_cloud_dc_of_nuclide = nuclide%ext_gamma_dc (air, effective)
  END FUNCTION fu_cloud_dc_of_nuclide

  REAL FUNCTION fu_fallout_dc_of_nuclide (nuclide)
    type(silam_nuclide), INTENT (in) :: nuclide
    fu_fallout_dc_of_nuclide = nuclide%ext_gamma_dc (ground, effective)
  END FUNCTION fu_fallout_dc_of_nuclide

  logical FUNCTION fu_if_noble_gas_of_nuclide (nuclide)
    type(silam_nuclide), INTENT (in) :: nuclide

    IF (nuclide%name (1 : 2) == 'H_' .or. &
        nuclide%name (1 : 2) == 'HE' .or. &
        nuclide%name (1 : 2) == 'N_' .or. &
        nuclide%name (1 : 2) == 'O_' .or. &
        nuclide%name (1 : 2) == 'F_' .or. &
        nuclide%name (1 : 2) == 'NE' .or. &
        nuclide%name (1 : 2) == 'CL' .or. &
        nuclide%name (1 : 2) == 'AR' .or. &
        nuclide%name (1 : 2) == 'KR' .or. &
        nuclide%name (1 : 2) == 'XE' .or. &
        nuclide%name (1 : 2) == 'RN') THEN
      fu_if_noble_gas_of_nuclide = .true.
    ELSE
      fu_if_noble_gas_of_nuclide = .false.
    END IF
  END FUNCTION fu_if_noble_gas_of_nuclide

  integer function fu_n_daughters(nuclide)
    type(silam_nuclide), INTENT (in) :: nuclide
    fu_n_daughters = nuclide%branching%n_daughters
  end function fu_n_daughters

  function fu_daughter(nuclide, iDaughter)
    character(len=nuc_name_len) :: fu_daughter
    type(silam_nuclide), INTENT (in) :: nuclide
    integer, intent(in) :: iDaughter
    if(iDaughter > 0 .and. iDaughter <= nuclide%branching%n_daughters)then
      fu_daughter = nuclide%branching%daughter(iDaughter)
    else
      call msg('Wrong daughter index:',iDaughter)
      call set_error('Wrong daughter index','fu_daughter')
      fu_daughter = ''
    endif
  end function fu_daughter

  real function fu_decay_branch_fraction(nuclide, iDaughter)
    type(silam_nuclide), INTENT (in) :: nuclide
    integer, intent(in) :: iDaughter
    if(iDaughter > 0 .and. iDaughter <= nuclide%branching%n_daughters)then
      fu_decay_branch_fraction= nuclide%branching%branch_fraction(iDaughter)
    else
      call msg('Wrong daughter index:',iDaughter)
      call set_error('Wrong daughter index','fu_decay_fraction')
      fu_decay_branch_fraction = 0.0
    endif
  end function fu_decay_branch_fraction

  function fu_decay_enegries(nuclide) result(pEnergies)
    !
    ! Returns energies and emission probabilities of a given nuclide
    !
    implicit none
    type(silam_nuclide), target, intent(in) :: nuclide
    real(r4k), dimension(:,:), pointer :: pEnergies
    integer :: i,n
    nullify(pEnergies)
    if(.not. (nuclide%defined == silja_true))then
      call set_error('Undefined nuclide','fu_decay_enegries')
      return
    endif
    do i = 1,size(nuclide%pDecayEnergy(1,:))
       if (nuclide%pDecayEnergy(1,i) == 0.0) then
          n = size(nuclide%pDecayEnergy(1,:i-1))
          exit
       end if
    end do
    pEnergies => nuclide%pDecayEnergy(:,:n)
  end function fu_decay_enegries
  
  subroutine copy_nuclide(nuc_out,nuc_in)
    ! Just to test where nuclide is copied
    ! error invoked if it is done somewhere
    implicit none
    type(silam_nuclide), intent(in) :: nuc_in
    type(silam_nuclide), intent(out) :: nuc_out
    
    call set_error('nuclide copied','copy_nuclide')
    
  end subroutine copy_nuclide
  
END MODULE nuclides
