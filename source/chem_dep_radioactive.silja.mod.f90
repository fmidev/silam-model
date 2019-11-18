
MODULE chem_dep_radioactive
  !
  ! The handling of radioactive transformations. Treatment of the radioactive 
  ! species is done via interface to the Mikko Ilvonen's dose assessment, for
  ! which the special transfer functions were made. Neither very effective
  ! nor elegant method, but seems to be mandatory for keeping the split
  ! between the dispersion and radioactive parts of the code.
  ! May be, later something more efficient will be generated.
  !
  ! Author: Mikhail Sofiev, FMI, mikhail.sofiev@fmi.fi
  !
  ! All units: SI
  ! 
  ! All degrees: real numbers in degrees (degrees and parts per
  ! hundred/parts per thousand , NOT in degrees and minutes)
  ! 
  ! Language: ANSI standard Fortran 90
  !
  use cocktail_basic
  !$use omp_lib

  IMPLICIT NONE

  ! The public functions and subroutines available in this module:
  !
  public init_radioactive
  public fu_if_specific_deposition

  public radioactive_input_needs
  public transform_radioactive
  public report

  public defined
  public set_chem_rules_radioactive
  public set_missing
  public fu_lifetime

  ! The private functions and subroutines not to be used elsewhere:
  !
  private add_radioact_decay_sp
  private fu_defined_radioactive_rules
  private set_missing_radioact_rules
  private fu_lifetime_radioact

  interface fu_if_specific_deposition
    module procedure fu_if_specific_dep_radioact
  end interface

  interface defined
    module procedure fu_defined_radioactive_rules
  end interface

  interface set_missing
    module procedure set_missing_radioact_rules
  end interface
  
  public fu_if_tla_required
  private fu_if_tla_required_radioact
  interface fu_if_tla_required
     module procedure fu_if_tla_required_radioact
  end interface


  !
  ! The decay can, in reality, be pre-computed following:
  ! dA/dt=-D*A => A(t+delta_t)=A(t)exp(-D*delta_t). Note that activity A can 
  ! be vector and the decay rate D will then be a matrix. This matrix or
  ! scalar exponent can be pre-computed since delta_t is constant throughout the run.
  ! A structure for that:
  !
  type Tprecomputed_decay
    private
    integer :: nuclide_starting_index        ! the first transport species, which is radioactive
    integer :: nNuclides                             ! number of nuclides == size of the matrix
    real*8, dimension(:,:), allocatable :: decay_matrix    ! precomputed matrix exponent
    type(silam_species), dimension(:), allocatable :: species  ! Corresponding nuclide-containing species
    type(silja_logical) :: defined
  end type Tprecomputed_decay

  !--------------------------------------------------------------------
  !
  !  Tchem_rules_passive type definition
  !  Determines the rules of operating with passive species with defined lifetime
  !
  type Tchem_rules_radioactive
    private
!    real :: massLowThreshold   ! Lower threshold for activity transported by Eulerian advection
    real :: vD_gas  ! A dry deposition velocity for gaseous species for dry smooth surface
    type(Tprecomputed_decay) :: precomputed_decay
    type(silja_logical) :: defined
  end type Tchem_rules_radioactive


  INTEGER, PARAMETER, public :: transformation_radioactive = 5002

!  REAL, DIMENSION(max_species),  private :: ptrDecayOutputTmp
!  !$OMP THREADPRIVATE(ptrDecayOutputTmp)

CONTAINS


  !************************************************************************************

  subroutine init_radioactive(rules, &
                         & speciesEmis, speciesTransp, &
                         & nSpeciesEmis, nSpTrn, &
                         & iClaimedSpecies, ifActive, timestep)
    implicit none
    !
    ! Adds decay products to speciesTransp,
    ! Claims radioactive spiecies to iClaimedSpecies
    !
    type(Tchem_rules_radioactive), intent(inout) :: rules
    type(silam_species), dimension(:), intent(in) :: speciesEmis
    type(silam_species), dimension(:), pointer :: speciesTransp
    integer, intent(in) :: nSpeciesEmis
    integer, intent(inout) :: nSpTrn
    integer, dimension(:), intent(inout) :: iClaimedSpecies !Claimed emission species
    logical, intent(out) :: ifActive  ! whether the radioactive transformation is active
    type(silja_interval), intent(in) :: timestep

    ! Local variables
    integer :: iEmis, nNucsEmis, nNucsAdded,  iSp
!    real, dimension(:), pointer :: fWork
    real, dimension(:,:), pointer :: reactTmp 
    character (len=15) :: strTmp
    character (len=*), parameter :: subname = "init_radioactive"

    !
    ! Each emission species has to be checked for valid radioactive features of the material.
    ! If yes, ownership is claimed, and species is added to speciesTransp
    !

    nNucsEmis = 0
    do iEmis = 1, nSpeciesEmis
      if(fu_if_radioactive(fu_material(speciesEmis(iEmis))))then
        if(iClaimedSpecies(iEmis) < 0)then
          call msg('Radioactive owns:' + fu_str(speciesEmis(iEmis)))
          iClaimedSpecies(iEmis) = transformation_radioactive
          nNucsEmis = nNucsEmis + 1
          call addSpecies(speciesTransp, nSpTrn, speciesEmis(iEmis:iEmis), 1)
        else
          call msg('Cannot claim ownership for:' + fu_str(speciesEmis(iEmis)) + &
                 & ', because he claimed it already:',iClaimedSpecies(iEmis))
          call set_error('Cannot claim ownership because someone claimed it already',subname)
          return
        endif  ! if can claim the emission species
      endif  ! emission species has valid radiological description
    end do  ! emission species

    if(nNucsEmis == 0)then
      ifActive = .false.
      return  ! If nothing to do - exit
    end if
    !
    ! Now we can proceed with the transformations species.
    !
!    fWork => fu_work_array(max_species*max_species)
    reactTmp => fu_work_array_2d()      !(1:max_species,1:max_species) => fWork(1:max_species*max_species)
    ! dv/dt = reactTmp * v
    reactTmp(:,:)=0.0

    nNucsAdded = nNucsEmis
    do iEmis = 1, nSpeciesEmis
      if(iClaimedSpecies(iEmis) == transformation_radioactive)then
!        call msg('Trying Species:' + fu_str(speciesEmis(iEmis)))
        call add_radioact_decay_sp(speciesEmis(iEmis), speciesTransp, nSpTrn, reactTmp, nNucsAdded)
      end if
      if(error)return
    end do
    if(nNucsAdded > 0)then
      ifActive = .not. error
    else
      call set_error('No nuclides in decay chain: nothing useful in the emission list?', subname)
      ifActive = .false.
    endif

    rules%precomputed_decay%nuclide_starting_index = nSpTrn - nNucsAdded + 1
    rules%precomputed_decay%nNuclides = nNucsAdded

    allocate(rules%precomputed_decay%decay_matrix(nNucsAdded,nNucsAdded), &
           & rules%precomputed_decay%species(nNucsAdded), stat=iSp)
    if(fu_fails(iSp == 0, 'Cannot allocate the precomputed decay matrix',subname))return
    
    ! A*t
    rules%precomputed_decay%decay_matrix(:,:) = fu_sec(timestep)*reactTmp(1:nNucsAdded,1:nNucsAdded)
    ! Save species
    rules%precomputed_decay%species(1:nNucsAdded) = speciesTransp(nSpTrn-nNucsAdded+1:nSpTrn)

    if (debug_level > 0) then
      call msg('Decay matrixd before the exponent:')
      do iSp = 1, nNucsAdded
        write(unit=strTmp,fmt='(A10,X,I3,A)') fu_str(rules%precomputed_decay%species(iSp)), iSp, ':'
        call msg(strTmp, rules%precomputed_decay%decay_matrix(1:nNucsAdded,iSp))
      end do
      call msg('')
    end if
    ! exp(A*t)
    call  matrix_exponent(rules%precomputed_decay%decay_matrix, nNucsAdded)
    if (debug_level > 0) then
      call msg('Decay matrixd after the exponent:')
      do iSp = 1, nNucsAdded
        write(unit=strTmp,fmt='(A10,X,I3,A)') fu_str(rules%precomputed_decay%species(iSp)), iSp, ':'
        call msg(strTmp, rules%precomputed_decay%decay_matrix(1:nNucsAdded,iSp))
      end do
      call msg('')
    end if
    call free_work_array(reactTmp)  !fWork)

  end subroutine init_radioactive


  !************************************************************************************

  recursive subroutine add_radioact_decay_sp(spIn, speciesTransp, nSpTrn, reactTmp, nNucsAdded)
    !
    ! Adds daughters to speciesTransp
    ! Fills reactTmp matrix with decay rates (1/tau)
    ! Do not care about order the species added
    !
    ! Particle-to-particle keeps mode
    ! Particle-to-gas forgets mode
    ! gas-to-particle choses smallest existing mode of reasonable size or default mode
    !

    implicit none
    
    ! Imported parameters
    type(silam_species), intent(in) :: spIn !Decaying species 
    type(silam_species), dimension(:), pointer :: speciesTransp !Transport species to which we add stuff
    integer, intent(inout) :: nSpTrn ! number of filled species in speciesTransp
    real, dimension (:,:), intent(inout) :: reactTmp ! decay rates 
    integer, intent(inout) :: nNucsAdded ! number of nucleides added to speciesTransp
!                Nucleides have indices (nSpTrn-nNucsAdded+1):nSpTrn

    ! Local variables
    type(silam_species), dimension(1) :: spDaughter
    type(silam_nuclide), pointer :: pNuc
    type(silam_material), pointer :: pSubst, pDaughterSubst
    integer :: iNucStart, iDec, iMother, iDaughter, iSp
    type(Taerosol_mode) :: mode, modeOut
    real :: d, decay_rate
    logical :: ifAerosol ! spIn is aerosol

    !
    ! If there are daughters, call recursion for each of them. 
    !
    pSubst => fu_material(spIn)
    if(.not.fu_if_radioactive(pSubst))then
      call set_error('Cannot find nuclide data in:'+fu_name(pSubst),'add_radioact_decay_sp')
      return
    endif


    mode  = fu_mode(spIn)
    pNuc => fu_nuclide(pSubst)
    ifAerosol = .not. (fu_if_gas(pSubst) == silja_true)
    iNucStart = nSpTrn - nNucsAdded + 1


    iMother = fu_index(spIn, speciesTransp, nSpTrn)
    if (iMother < 0) then
      call msg("Could not find species "//trim(fu_str(spIn))//" in speciesTransp")
      call report(speciesTransp(1:nSpTrn))
      call set_error("No mother species included","add_radioact_decay_sp")
          return
        endif
    
    iMother = iMother - iNucStart + 1 !! Now -- index in radioactive subset

    decay_rate =  ln_2 / fu_half_life(pNuc)
    reactTmp(iMother,iMother) = - decay_rate 


    do iDec = 1, fu_n_daughters(pNuc) ! Does not execute if no daughters

      pDaughterSubst => fu_get_material_ptr(fu_daughter(pNuc,iDec))

      !Get or guess daughter mode
      if (fu_if_gas(pDaughterSubst) == silja_true) then
        modeOut = in_gas_phase
      elseif (ifAerosol) then 
        !Inherit mode from mother
        modeOut=mode
      else
        ! Try to find existing species with suitable mode
        modeOut = aerosol_mode_missing
        do iSp = iNucStart, nSpTrn
         if (associated(fu_material(speciesTransp(iSp)), pDaughterSubst)) then
           d = fu_nominal_d(speciesTransp(iSp))
           if (d < 5.e-8 .or. d > 2.e-6) cycle !! Wrong-sized mode
           if (defined(modeOut)) then
             if (fu_nominal_d(modeOut) < d) modeOut = fu_mode(speciesTransp(iSp))
           else
             modeOut = fu_mode(speciesTransp(iSp))
    endif
        endif
    enddo
       ! Last resort: just invent a mode 
       if (.not. defined(modeOut)) call set_aerosol_mode(modeOut, 0.5e-6)
      endif

      call set_species(spDaughter(1), pDaughterSubst, modeOut)
      call addSpecies(speciesTransp, nSpTrn, spDaughter, 1)
      nNucsAdded = nSpTrn - iNucStart + 1

      iDaughter = fu_index(spDaughter(1), speciesTransp, nSpTrn) - iNucStart + 1

      reactTmp(iMother,iDaughter) = decay_rate * fu_decay_branch_fraction(pNuc, iDec)

      !! Handle further generations
      call add_radioact_decay_sp(spDaughter(1), speciesTransp, nSpTrn, reactTmp, nNucsAdded)
    end do


  end subroutine add_radioact_decay_sp

  !************************************************************************************

  logical function fu_if_specific_dep_radioact(rules)
    !
    ! Turn into TRUE if this module is to control the deposition computations
    !
    implicit none
    
    type(Tchem_rules_radioactive), intent(in) :: rules

    fu_if_specific_dep_radioact = .false.
    
  end function fu_if_specific_dep_radioact

  !************************************************************************************

  subroutine radioactive_input_needs(rules, meteo_input_local)
    implicit none
    type(Tchem_rules_radioactive), intent(in) :: rules
    type(Tmeteo_input), intent(out) :: meteo_input_local

    meteo_input_local = meteo_input_empty

  end subroutine radioactive_input_needs


  !***********************************************************************

  subroutine transform_radioactive(mass_vector_1d_tr, rules, metdat, timestep_sec, print_it)
    !
    ! Implements radioactive decay
    !
    implicit none

    ! Imported parameters
    real, dimension(:), intent(inout) :: mass_vector_1d_tr
    type(Tchem_rules_radioactive), intent(in) :: rules
    real, dimension(:), intent(in) :: metdat
    real, intent(in) :: timestep_sec
    logical, intent(out) :: print_it

    ! Local variables
    integer :: iSp, iSp2, iTMp, jTmp, nuc_start_ind, nNuclides
    type(silam_sp) :: sp
    real, dimension(:), pointer :: ptrDecayOutputTmp
    
    print_it = .false.  ! set to true and chemistry manager will give complete dump for this cell

    nuc_start_ind = rules%precomputed_decay%nuclide_starting_index

    if(any(mass_vector_1d_tr(nuc_start_ind : &
                           & nuc_start_ind &
                           & + rules%precomputed_decay%nNuclides-1) > 0.0)) then

      if(debug_level > 1)then
        sp%sp => fu_work_string()
        call msg('Decay matrix at the transformation:')
        do iSp = 1, rules%precomputed_decay%nNuclides
    write(unit=sp%sp,fmt='(I2,200(1x,F9.7))')iSp,(rules%precomputed_decay%decay_matrix(iSp2,iSp),iSp2=1,rules%precomputed_decay%nNuclides)
          call msg(fu_str(rules%precomputed_decay%species(iSp)) + ':' + sp%sp)
        end do
        call msg('Mass vector before:')
        write(unit=sp%sp,fmt='(2x,200(1x,E9.3))')&
             & (mass_vector_1d_tr(iSp),iSp=rules%precomputed_decay%nuclide_starting_index, &
              & rules%precomputed_decay%nuclide_starting_index + &
              & rules%precomputed_decay%nNuclides - 1)
        call msg(sp%sp)
      endif

      ptrDecayOutputTmp => fu_work_array()
      do jTmp = 1, rules%precomputed_decay%nNuclides
        ptrDecayOutputTmp(jTmp) = 0
        do iTmp = 1, rules%precomputed_decay%nNuclides
          ptrDecayOutputTmp(jTmp) = ptrDecayOutputTmp(jTmp) + rules%precomputed_decay%decay_matrix(iTmp,jTmp) * &
               & mass_vector_1d_tr(iTmp + nuc_start_ind - 1)
        enddo
      end do
      do iTmp = 1, rules%precomputed_decay%nNuclides
        mass_vector_1d_tr(nuc_start_ind + iTmp - 1) = ptrDecayOutputTmp(iTmp) 
      end do
      call free_work_array(ptrDecayOutputTmp)
      
      !     mass_vector_1d_tr(rules%precomputed_decay%nuclide_starting_index : &
      !                    & rules%precomputed_decay%nuclide_starting_index + &
      !                                               & rules%precomputed_decay%nNuclides-1) = &
      !           & MATMUL(rules%precomputed_decay%decay_matrix, &
      !                  & mass_vector_1d_tr(rules%precomputed_decay%nuclide_starting_index : &
      !                                    & rules%precomputed_decay%nuclide_starting_index + &
      !                                                       & rules%precomputed_decay%nNuclides - 1))
      if(debug_level > 1)then
        call msg('Mass vector after:')
        write(unit=sp%sp,fmt='(2x,200(1x,E9.3))') &
             & (mass_vector_1d_tr(iSp),iSp=rules%precomputed_decay%nuclide_starting_index, &
              & rules%precomputed_decay%nuclide_starting_index + &
              & rules%precomputed_decay%nNuclides - 1)
        call msg(sp%sp)
        call free_work_array(sp%sp)
      endif
    endif
  end subroutine transform_radioactive


  !*****************************************************************
  !
  ! Chemistry rules for radioactive cocktails
  !
  !*****************************************************************

  
  !***********************************************************************

  subroutine set_chem_rules_radioactive(nlSetup, rulesRadioactive)
    !
    ! Creates the set of rules, in accordance with which the computations are to
    ! be made. Input - SILAM namelist
    !
    implicit none

    ! Imported parameter
    type(Tsilam_namelist), pointer :: nlSetup !, nlStdSetup
    type(Tchem_rules_radioactive), intent(out) :: rulesRadioactive

    rulesRadioactive%defined = silja_false

    if(.not.defined(nlSetup))then
      call set_error('Undefined namelist given','set_chem_rules_radioactive')
      return
    endif
    
    rulesRadioactive%defined = silja_true

  end subroutine set_chem_rules_radioactive


  !***********************************************************************

  subroutine set_missing_radioact_rules(rulesRadioactive)
    !
    ! Creates the set of rules, in accordance with which the computations are to
    ! be made. Input - SILAM namelist
    !
    implicit none

    ! Imported parameter
    type(Tchem_rules_radioactive), intent(out) :: rulesRadioactive

    rulesRadioactive%defined = silja_false

!    rulesRadioactive%massLowThreshold = real_missing

  end subroutine set_missing_radioact_rules


  !**********************************************************************
  
  logical function fu_defined_radioactive_rules(rulesRadioactive)
    !
    ! Checks whether the rules are reasonable
    !
    implicit none
    
    type(Tchem_rules_radioactive), intent(in) :: rulesRadioactive
    
    fu_defined_radioactive_rules = fu_true(rulesRadioactive%defined)
    
  end function fu_defined_radioactive_rules


  !********************************************************************************

  function fu_lifetime_radioact(rules)result(lifetime)
    !
    ! A typical life time due to degradation and deposition for radioactive materials
    ! So far, whatever the species are, we set the "typical" lifetime to be two days
    !
    implicit none

    type(silja_interval) :: lifetime
    type(Tchem_rules_radioactive), intent(in) :: rules

    lifetime = one_day * 2.0

  end function fu_lifetime_radioact

  !************************************************************************************

  logical function fu_if_tla_required_radioact(rules) result(required)
    ! Collect transformations' requests for tangent linearization. If tangent linear
    ! is not allowed, corresponding subroutine sets error. 
    implicit none
    type(Tchem_rules_radioactive), intent(in) :: rules
    
    call set_error('TLA not available', 'fu_if_tla_required')
    required = .false.

  end function fu_if_tla_required_radioact


END MODULE chem_dep_radioactive

