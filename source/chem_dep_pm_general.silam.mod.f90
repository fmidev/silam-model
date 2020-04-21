MODULE chem_dep_pm_general
  !
  ! This module contains a set of routines that describe the processes related to 
  ! transformation and some (non-default) features of deposition of general aerosol species. 
  ! Essentially, this is the turbulent thermal diffusion stuff. Note that this is NOT
  ! the aerosol dynamics but rather a "general-aerosol chemistry", i.e. belongs to the
  ! same line as all other chemical modules.
  !
  ! All units: SI, unless otherwise stated
  !
  ! Code owner Mikhail Sofiev mikhail.sofiev@fmi.fi
  !
  ! Language: ANSI FORTRAN-90 (or close to it)
  !
  use cocktail_basic

  implicit none
  private

  !
  ! Specific functions for the general-aerosol transformation etc
  !
  public pm_proc_input_needs
  public init_inert_pm_emis_chem
  public inventoryPM
  public registerSpeciesPM
  public set_chemRules_PM_general
  public fu_if_specific_deposition

  private fu_if_specific_dep_PM

  interface fu_if_specific_deposition
    module procedure fu_if_specific_dep_PM
  end interface

  public fu_if_tla_required
  private fu_if_tla_required_pm
  interface fu_if_tla_required
     module procedure fu_if_tla_required_pm
  end interface


  !--------------------------------------------------------------------
  !
  ! Tchem_rules_aerosol_gen type definition
  ! Determines the rules of operating with passive species with defined lifetime
  ! So far, define only way the depostiion is handled and whether thermodiffusion
  ! is included into the processes
  !
  type Tchem_rules_pm_general
    private
    integer :: dryDepositionType, wetDepositionType
    real :: ratio_to_standardScavenging
!    real :: massLowThreshold
    type(silja_logical) :: defined
  end type Tchem_rules_pm_general
  public Tchem_rules_pm_general

  integer, private, parameter :: dry_dep_standard_full = 11001   ! resistances and sedimentation, ref to deposition module
  integer, private, parameter :: wet_dep_standard_scaled = 11011 ! 3D scavenging coef, possibly scaled, ref to dep module

  !
  ! Each moduel has to know own species it treats.
  ! The indices of these species are different in emission, transport, and short-lived mass maps
  ! The corresponding indices have to be stored in each module. For inert PM the situation is more
  ! difficult: there can be arbitrary number fo them due to many size modes. Hense, here we solve 
  ! the problem via an allocatable array
  !
  integer, dimension(:), pointer, private, save :: indPM_emis, indPM_transp, indPM_short_lived, indPM_optic
  integer, private, save :: nPM_emis=0, nPM_transp=0, nPM_short_lived=0, nPM_optic=0

  !
  ! Label of this module
  !
  integer, parameter, public :: transformation_inert_PM = 5004

  CONTAINS

  !***********************************************************************************
  !
  !   Input data request and handling
  !
  !***********************************************************************************

  subroutine pm_proc_input_needs(rulesPM, meteo_input_local)
    !
    ! Returns the required fields for the thermodiffusion computations
    !
    implicit none

    ! Imported parameters
    type(Tchem_rules_pm_general), intent(in) :: rulesPM
    type(Tmeteo_input), intent(out), target :: meteo_input_local

    ! Local variable
    integer :: iTmp

    meteo_input_local = meteo_input_empty

    !
    ! In fact, nothing special: inert PM is not transformed. Additional items
    ! can come from depositon part if it goes in some specific way but so far nothing 
    ! is in mind (10.2009)

  end subroutine pm_proc_input_needs



  !**************************************************************************************
  !**************************************************************************************
  !
  !    SILAM general-PM specific routines
  !
  !**************************************************************************************
  !**************************************************************************************
  !
  ! Since the general-aerosol does not have any chemical transformation, this part
  ! is practically void here. 
  ! The default-aerosol deposition processes are handled in the deposition module. 
  ! Here we can overrule them but so far there is nothing in mind (Oct 2009)
  ! Hense, the functions below are almost all empty.
  !

  !***********************************************************************

  subroutine init_inert_PM_emis_chem(speciesEmis, nSpecies, species_picked)
    !
    ! Stores the indices of the inert-PM species in the emission species list. 
    ! The chemical database is supposed to be read already
    !
    implicit none

    ! Imported parameters
    type(silam_species), dimension(:), intent(in) :: speciesEmis    ! all species from emission cocktails
    integer, intent(in) :: nSpecies                ! number of thes species in the list
    integer, dimension(:), intent(inout) :: species_picked ! the recognised ones are marked here

    ! Local variables
    integer :: iSpecies
    integer, dimension(:), pointer :: arTmp

    arTmp => fu_work_int_array()
    if(error)return

    !
    ! Go through the list and store the species, which are recognised to be the inert PM.
    ! Mark the corresponding place in the specie_picked
    !
    nPM_emis = 0
    do iSpecies = 1, nSpecies
      if(fu_str_l_case(fu_name(fu_material(speciesEmis(iSpecies)))) == 'PM')then
        species_picked(iSpecies) = species_picked(iSpecies) + 1
        nPM_emis = nPM_emis + 1
        arTmp(nPM_emis) = iSpecies
      endif
    end do
    !
    ! Aloocate the space for the indices of inert PM
    !
    allocate(indPM_emis(nPM_emis), stat = iSpecies)
    if(iSpecies /= 0)then
      call set_error('Failed memory allocation','init_inert_PM_emis_chem')
      return
    endif
    !
    ! and store the indices of inert PM
    !
    indPM_emis(1:nPM_emis) = arTmp (1:nPM_emis)

    call free_work_array(arTmp)

  end subroutine init_inert_PM_emis_chem


  !************************************************************************************

  subroutine inventoryPM(rules, &
                       & speciesEmis, speciesTransp, speciesShortlived, speciesAerosol,&
                       & nSpeciesEmis, nSpeciesTransp, nSpeciesShortlived, nSpeciesAerosol, &
                       & iClaimedSpecies)
    ! 
    ! The pm transformation does not actually have an inventory:
    ! instead, it just adds the emission species with the pm
    ! material. We will assume that speciesEmission contains no
    ! duplicates!
    !
    implicit none

    ! Imported parameters
    type(Tchem_rules_pm_general), intent(in) :: rules
    type(silam_species), dimension(:), pointer :: speciesEmis, speciesTransp, &
                                                & speciesShortlived, speciesAerosol
    integer, intent(in) :: nSpeciesEmis
    integer, intent(inout) :: nSpeciesTransp, nspeciesShortlived, nspeciesAerosol
    integer, dimension(:), intent(inout) :: iClaimedSpecies
    
    ! Local variables
    integer :: i
    type(silam_material), pointer :: pm_ptr

    pm_ptr => fu_get_material_ptr('PM')
    if (error) return
    nPM_emis = 0

!    nullify(speciesTransp)
!    nSpeciesTransp = 0

    !
    ! Scan the emission looking for PM material, claim ownership if found
    !
    do i = 1, nSpeciesEmis
      if (associated(speciesEmis(i)%material, pm_ptr)) then
        nPM_emis = nPM_emis + 1
        call addSpecies(speciesTransp, nSpeciesTransp, (/speciesEmis(i)/), 1, .true.)

        if(iClaimedSpecies(i) < 0)then
          call msg('Inert PM owns:' + fu_substance_name(speciesEmis(i)))
          iClaimedSpecies(i) = transformation_inert_pm
        else
          call msg('Cannot claim ownership because he claimed it already:',iClaimedSpecies(i))
          call set_error('Cannot claim ownership because someone claimed it already','inventoryPM')
          return
        endif
      end if
    end do

!    nSpeciesTransp = nPM_emis
!    nullify(speciesShortlived)
!    nSpeciesShortlived = 0
!    nullify(speciesAerosol)
!    nSpeciesAerosol = 0

  end subroutine inventoryPM


  !***********************************************************************

  subroutine registerSpeciesPM(rules, &
                             & speciesTransp, speciesShortlived, speciesAerosol,&
                             & nspeciesTransp, nspeciesShortlived, nspeciesAerosol)
    ! 
    ! Find the indices of passive species and store them to the module
    ! variable. Number of passive species in transport must be equal
    ! to their number of in emission species, since emission species
    ! is a subset of transport species.
    !
    implicit none
    type(Tchem_rules_pm_general), intent(in) :: rules
    type(silam_species), dimension(:), pointer :: speciesTransp, speciesShortlived, speciesAerosol
    integer, intent(in) :: nspeciesTransp, nspeciesShortlived, nspeciesAerosol

    integer :: i
    type(silam_material), pointer :: pm_ptr

    pm_ptr => fu_get_material_ptr('PM')

    allocate(indPM_transp(nPM_emis), stat=i)
    if (i /= 0) then
      call set_error('Allocate failed', 'registerSpeciesPM')
      return
    end if

    nPM_transp = 0

    do i = 1, nSpeciesTransp
      if (associated(speciesTransp(i)%material, pm_ptr)) then
        nPM_transp = nPM_transp + 1
        if (nPM_transp > nPM_emis) then
          ! Not that critical, actually, but currently I haven't
          ! bothered to accommodate this so far.
          call set_error('More PM species in transport than in emission, won''t work',&
                       & 'registerSpeciesPM')
          return
        end if
        indPM_transp(nPM_transp) = i

      end if
    end do
    
  end subroutine registerSpeciesPM


!  !***********************************************************************
!
!  subroutine full_species_lst_4ckt_inert_pm(speciesEmis, materialsAvailable, &
!                                          & speciesTransp, nSpeciesTransp, &
!                                          & speciesShortLived, nSpeciesShortLived, &
!                                          & ifFixedNumbering)
!    !
!    ! Here we select the species to be transported and short-lived having the 
!    ! list of emitted species and the list of available ones (the SILAM database).
!    ! Also, an option for fixed order can be executed here.
!    ! The inert PM species go one-to-one from emission to transport cocktail.
!    !
!    implicit none
!
!    ! Imported parameters
!    type(silam_species), dimension(:), pointer :: speciesEmis  ! all species from emission cocktails
!    type(silam_material), dimension(:), pointer :: materialsAvailable  ! the SILAM chemical database
!    type(silam_species), dimension(:), pointer :: speciesTransp, speciesShortLived
!    integer, intent(out) :: nSpeciesTransp, nSpeciesShortLived
!    logical, intent(out) :: ifFixedNumbering
!
!    ! Local variable
!    integer :: iSpecies
!
!    ifFixedNumbering = .false.  ! order does not matter
!    nSpeciesTransp = nPM_emis   ! no new species, just the emitted stuff
!    nSpeciesShortLived = 0      ! no short-lived PM speceis
!
!    ! Emitted and transported species are the same, which simplifies the story.
!    ! However, we follow the general procedure here, which requires the reference in the
!    ! speceisAvailable list, not from the speciesEmis list. This is because other 
!    ! transformation modules may have non-emitted species as transported and short-lived 
!    !
!    do iSpecies = 1, nPM_emis
!      speciesTransp(iSpecies) = speciesEmis(indPM_emis(iSpecies)) ! no need to make new, just copy
!    end do
!
!  end subroutine full_species_lst_4ckt_inert_pm


  !***************************************************************************************
  !
  !  Aerosol-general rules
  !
  !***************************************************************************************

  subroutine set_chemRules_PM_general(nlSetup, rulesPM)
    !
    ! Sets the aerosol chemical rules from the given namelist
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist), pointer :: nlSetup
    type(Tchem_rules_pm_general), intent(out) :: rulesPM

    rulesPM%defined = silja_false

    if(.not.defined(nlSetup))then
      call set_error('Undefined namelist given','set_chemRules_PM_general')
      return
    endif
    
    !
    ! Here is the place to request any non-standard deposition, if they are going to be computed
    ! in this module. So far, we do not have any, so nothing to worry about
    !
    rulesPM%dryDepositionType = int_missing
    rulesPM%wetDepositionType = int_missing
    rulesPM%ratio_to_standardScavenging = real_missing

!    if(fu_str_u_case(fu_content(nlSetup,'dry_deposition_scheme')) == 'STANDARD')then
!      rulesPM%dryDepositionType = dry_dep_standard_full
!    else
!      call set_error('dry_deposition_scheme can be only STANDARD so far','set_chemRules_PM_general')
!      return
!    endif
!
!    if(fu_str_u_case(fu_content(nlSetup,'wet_deposition_scheme')) == 'STANDARD')then
!      rulesPM%wetDepositionType = wet_dep_standard_scaled
!      rulesPM%ratio_to_standardScavenging = 1.0
!
!    elseif(fu_str_u_case(fu_content(nlSetup,'wet_deposition_scheme')) == 'STANDARD_SCALED')then
!      rulesPM%wetDepositionType = wet_dep_standard_scaled
!      rulesPM%ratio_to_standardScavenging = fu_content_real(nlSetup,'ratio_to_standard_scavenging')
!      if(rulesPM%ratio_to_standardScavenging .eps. real_missing)then
!        call set_error('Scaling wet deposition scheme requested but no ratio_to_standard_scavenging', &
!                     & 'set_chemRules_PM_general')
!      endif
!
!    else
!      call set_error('wet_deposition_scheme can be only STANDARD or STANDARD_SCALED', &
!                   & 'set_chemRules_PM_general')
!      return
!    endif

    !
    ! For Eulerian scheme we will need the minimum threshold of amount of species in the grid cell
    ! The cell will not be processed by ANY Eulerian-pool routine should it have smaller amount 
    ! than that threshold. A data-based setting is possible only when the total emitted amount
    ! is known, as well as a grid size. So far let's put it to something small enough.
    !
!    rulesPM%massLowThreshold = 1.0e-13

    rulesPM%defined = silja_true

  end subroutine set_chemRules_PM_general


  !**********************************************************************************

  logical function fu_if_specific_dep_PM(rulesPM)
    implicit none
    type(Tchem_rules_pm_general), intent(in) :: rulesPM
    fu_if_specific_dep_PM = .false.     ! no specific deposition whatsoever
  end function fu_if_specific_dep_PM

  !************************************************************************************
  

  logical function fu_if_tla_required_pm(rules) result(required)
    ! Collect transformations' requests for tangent linearization. If tangent linear
    ! is not allowed, corresponding subroutine sets error. 
    implicit none
    type(Tchem_rules_pm_general), intent(in) :: rules
    
    required = .false.

  end function fu_if_tla_required_pm
  

end module chem_dep_pm_general


