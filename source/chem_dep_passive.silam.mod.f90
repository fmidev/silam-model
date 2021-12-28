MODULE chem_dep_passive
  !
  ! This module contains a set of routines that describe the transformation and
  ! some (non-default) features of deposition of passive species. It
  ! does not react with anything but can degrade by itself reproducing a 
  ! substance with tricky/unknown features but defined lifetime.
  !
  ! All units: SI, unless otherwise stated
  !
  ! Code owner Mikhail Sofiev mikhail.sofiev@fmi.fi
  !
  ! Language: ANSI FORTRAN-90 (or close to it)
  !
  use cocktail_basic

  implicit none
  
  !
  ! public routines for PASSIVE cocktail
  !
  !public init_passive_emis_chemicals
  !public full_species_lst_4ckt_passive
  public inventoryPassive
  public registerSpeciesPassive
  public passive_proc_input_needs
  public transform_passive

  ! public routines for passive chemistry rules
  !
  public set_chem_rules_passive
  public fu_lifetime
  public fu_if_specific_deposition
  public if_need_init_time_tag
  public if_need_init_ones


  ! Private routines for passive chemistry rules
  !
  private fu_lifetime_passive
  private fu_if_specific_dep_passive

  interface fu_lifetime
    module procedure fu_lifetime_passive
  end interface

  interface fu_if_specific_deposition
    module procedure fu_if_specific_dep_passive
  end interface

  public fu_if_tla_required
  private fu_if_tla_required_passive
  interface fu_if_tla_required
     module procedure fu_if_tla_required_passive
  end interface

  !--------------------------------------------------------------------
  !
  !  Tchem_rules_passive type definition
  !  Determines the rules of operating with passive species with defined lifetime
  !
  type Tchem_rules_passive
    private
    real :: basicTempr, basicLifeTime, DLifeTime_DT
    logical :: ifTimeTracer = .false., ifInitialisedTimeTag =.false.
    logical :: ifOnesTracer = .false., ifInitialisedOnes =.false.
    logical :: ifSF6destroy = .false.
    real :: SF6LifetimeUpperLayer = real_missing
    type(silja_logical) :: defined
  end type Tchem_rules_passive

  !-------------------------------------------------------------------------
  !
  ! Pointers to the meteorological fields, which can be requested
  ! by the corresponding transformation or deposition routines
  !
  real, dimension(:), private, pointer, save :: ptrTempr2m
  type(field_4d_data_ptr), private, pointer, save :: ptrTempr

  integer, private, pointer, save :: ind_tempr  ! Used for addressing the meteo data
  integer, private, pointer, save :: ind_pressure  ! Used for addressing the meteo data

  !-------------------------------------------------------------------------
  !
  ! The main array pointing at the passive species in the transport cocktail
  !
  integer, dimension(:), pointer, private, save :: indPass_transp
  !, indPass_short_lived, indPass_optic
  integer, private, save :: indPassTime_transp = int_missing ! A dedicated index for time tag
  integer, private, save :: indPassOnes_transp = int_missing ! A dedicated index for ones

  ! SF6 kinds that can be destroyed
  integer, private, save :: indSF6_transp = int_missing ! A dedicated index for SF6
  integer, private, save :: indSF6ng_transp = int_missing ! A dedicated index for SF6 with no moecilar diff
  integer, private, save :: nPass_emis = 0, nPass_transp = 0
  !, nPass_short_lived = 0, nPass_optic = 0


  ! The mark of this transformation type
  !
  INTEGER, PARAMETER, public :: transformation_passive = 5001


  CONTAINS


  !************************************************************************************

  subroutine inventoryPassive(rules, &
                            & speciesEmis, speciesTransp, speciesShortlived, speciesAerosol,&
                            & nSpeciesEmis, nspeciesTransp, nspeciesShortlived, nspeciesAerosol, &
                            & iClaimedSpecies)
    implicit none
    ! 
    ! The passive transformation adds the emission species with the passive
    ! material. We will assume that speciesEmission contains no duplicates
    ! The only material that might need to be added is the passive time tracer: a special
    ! material, which does not need to be emitted and which concentration grows linearly with time
    !
    type(Tchem_rules_passive), intent(in) :: rules
    type(silam_species), dimension(:), pointer :: speciesEmis, speciesTransp, speciesShortlived, &
                                                & speciesAerosol
    integer, intent(in) :: nspeciesEmis
    integer, intent(inout) :: nspeciesTransp, nspeciesShortlived, nspeciesAerosol
    integer, dimension(:), intent(inout) :: iClaimedSpecies

    ! Local variables
    integer :: i
    type(silam_material), pointer :: passive_ptr
    type(silam_species) :: speciesTmp

    passive_ptr => fu_get_material_ptr('passive')
    if (error) return
    nPass_emis = 0

    !
    ! Scan the emission species, find passive material and put it to transport mass map, also claiming
    ! ownership
    !
    do i = 1, nspeciesEmis
      if (associated(speciesEmis(i)%material, passive_ptr)) then
        nPass_emis = nPass_emis + 1
        call addSpecies(speciesTransp, nSpeciesTransp, (/speciesEmis(i)/), 1)

        if(iClaimedSpecies(i) < 0)then
          call msg('Passive owns:' + fu_substance_name(speciesEmis(i)))
          iClaimedSpecies(i) = transformation_passive
        else
          call msg('Cannot claim ownership because he claimed it already:',iClaimedSpecies(i))
          call set_error('Cannot claim ownership because someone claimed it already','inventory_passive')
          return
        endif
      endif  ! emission species belongs to passive transport one
    end do  ! emission species
    !
    ! Taking care of the time tracer
    !


    if(rules%ifTimeTracer)then
      call set_species(speciesTmp, fu_get_material_ptr('timetag'), in_gas_phase)
      if (fu_fails( fu_index(speciesTmp, speciesTransp, nSpeciesTransp) <0 ,&
                   & "Can't add timetag twice",'inventory_passive')) return
      call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
    endif

    if(rules%ifOnesTracer)then
      call set_species(speciesTmp, fu_get_material_ptr('ones'), in_gas_phase)
      if (fu_fails( fu_index(speciesTmp, speciesTransp, nSpeciesTransp) <0 ,&
                   & "Can't add ones twice",'inventory_passive')) return
      call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
    endif

    if(rules%ifSF6destroy)then
      call set_species(speciesTmp, fu_get_material_ptr('SF6'), in_gas_phase)
      if (fu_fails( fu_index(speciesTmp, speciesTransp, nSpeciesTransp) <0 ,&
                   & "Can't add SF6 twice",'inventory_passive')) return
      call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
      call set_species(speciesTmp, fu_get_material_ptr('SF6ng'), in_gas_phase)
      if (fu_fails( fu_index(speciesTmp, speciesTransp, nSpeciesTransp) <0 ,&
                   & "Can't add SF6ng twice",'inventory_passive')) return
      call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1)
    endif

  end subroutine inventoryPassive


  !***********************************************************************

  subroutine registerSpeciesPassive(rules, &
                                  & speciesTransp, speciesShortlived, speciesAerosol,&
                                  & nspeciesTransp, nspeciesShortlived, nspeciesAerosol)
    ! 
    ! Find the indices of passive species and store them to the module
    ! variable. Number of passive species in transport must be equal
    ! to their number of in emission species, since emission species
    ! is a subset of transport species.
    !
    implicit none
    type(Tchem_rules_passive), intent(in) :: rules
    type(silam_species), dimension(:), pointer :: speciesTransp, speciesShortlived, speciesAerosol
    integer, intent(in) :: nspeciesTransp, nspeciesShortlived, nspeciesAerosol
    
    integer :: i
    type(silam_material), pointer :: passive_ptr, time_ptr, ones_ptr, sf6_ptr, sf6ng_ptr

    passive_ptr => fu_get_material_ptr('passive')
    time_ptr  => null()
    ones_ptr  => null()
    sf6_ptr   => null()
    sf6ng_ptr => null()
    if(rules%ifTimeTracer) time_ptr => fu_get_material_ptr('timetag')
    if(rules%ifOnesTracer) ones_ptr => fu_get_material_ptr('ones')
    if(rules%ifSF6destroy) sf6_ptr => fu_get_material_ptr('SF6')
    if(rules%ifSF6destroy) sf6ng_ptr => fu_get_material_ptr('SF6ng')
    
    allocate(indPass_transp(nPass_emis), stat=i)
    if(fu_fails(i == 0,'Allocate failed', 'registerSpeciesPassive'))return

    npass_transp = 0
    
    do i = 1, nSpeciesTransp
      if (associated(speciesTransp(i)%material, passive_ptr)) then
        npass_transp = npass_transp + 1
        if (npass_transp > nPass_emis) then
          ! Not that critical, actually, but currently I haven't
          ! bothered to accommodate this so far.
          call set_error('More passive species in transport than in emission, won''t work',&
                       & 'registerSpeciesPassive')
          return
        end if
        indPass_transp(npass_transp) = i
      end if
      !
      ! A separate timetag species
      !
      call msg("Species", i, nSpeciesTransp)
      call report(speciesTransp(i)%material)
      if(rules%ifTimeTracer)then
        if (associated(speciesTransp(i)%material, time_ptr)) indPassTime_transp = i
      endif
      if(rules%ifOnesTracer)then
        if (associated(speciesTransp(i)%material, ones_ptr)) indPassOnes_transp = i
      endif
      if(rules%ifSF6destroy) then
        if (associated(speciesTransp(i)%material, sf6_ptr)) indSF6_transp = i
        if (associated(speciesTransp(i)%material, sf6ng_ptr)) indSF6ng_transp = i
      endif
    end do  ! transport species

    if ( rules%ifSF6destroy .and. indSF6_transp <= 0 .and. indSF6ng_transp <= 0) then
      call msg("indSF6_transp indSF6ng_transp", indSF6_transp, indSF6ng_transp)
        call set_error("Destroy SF6 requested, but no destructable SF6 in the run..", &
           & 'registerSpeciesPassive')
        return
    endif
    
  end subroutine registerSpeciesPassive


  !*****************************************************************

  subroutine passive_proc_input_needs(rulesPassive, meteo_input_local)
    !
    ! Returns input needs for the transformation routines. 
    !
    implicit none

    ! Imported parameters
    type(Tchem_rules_passive), intent(in) :: rulesPassive
    type(Tmeteo_input), intent(out), target :: meteo_input_local

    ! Local variables
    integer :: iTmp

    meteo_input_local = meteo_input_empty

    if(.not.(rulesPassive%defined == silja_true))return

    iTmp = 0 !Free index in meteo input
    if (rulesPassive%dLifeTime_DT /= 0.0) then 
      iTmp = iTmp + 1
      meteo_input_local%quantity(iTmp) = temperature_flag    ! For self-degradation
      meteo_input_local%q_type(iTmp) = meteo_dynamic_flag
      ind_tempr => meteo_input_local%idx(iTmp)
      call msg("rulesPassive%dLifeTime_DT", rulesPassive%dLifeTime_DT)
    endif
    if(rulesPassive%ifSF6destroy) then
      iTmp = iTmp + 1
      meteo_input_local%quantity(iTmp) = pressure_flag    ! For self-degradation
      meteo_input_local%q_type(iTmp) = meteo_dynamic_flag
      ind_pressure => meteo_input_local%idx(iTmp)
      call msg("rulesPassive%ifDestroySf6 enabled")
    endif

    meteo_input_local%nQuantities = iTmp
    meteo_input_local%defined = silja_true

  end subroutine passive_proc_input_needs



  !*******************************************************************

  subroutine transform_passive(vMass, rulesPassive, metdat, iz,  timestep_sec, print_it) 
    !
    ! Computes the transformation for the passive species in the given array of masses 
    ! using the provided meteodata
    !
    implicit none

    ! Imported parameters
    real, dimension(:), intent(inout) :: vMass
    real, dimension(:), intent(in) :: metdat
    integer, intent(in) :: iz
    type(Tchem_rules_passive), intent(in) :: rulesPassive
    real, intent(in) :: timestep_sec
    logical, intent(out) :: print_it

    ! Local parameters
    integer :: iSpecies
    real :: fLifeTime, fTmp
    integer, save :: iCount=0

    print_it = .false.  ! set to true and chemistry manager will give complete dump for this cell

    if(rulesPassive%defined == silja_true)then
      !
      ! Just do the transformation
      !
      if(nPass_transp > 0)then
        if ( rulesPassive%DLifeTime_DT == 0) then 
            fLifeTime = rulesPassive%basicLifeTime
        else
           fLifeTime = rulesPassive%basicLifeTime + &
                    & (metdat(ind_tempr) - rulesPassive%basicTempr) * rulesPassive%DLifeTime_DT
           if(fLifeTime < 0.01)then
             call msg('Warning. Very short lifetime. temperature and lifetime:', metdat(ind_tempr), fLifeTime)
             if(iCount < 100)then
               call msg('basic tempr and dtau/dT:',rulesPassive%basicTempr, rulesPassive%DLifeTime_DT)
               iCount = iCount + 1
             endif
             fLifeTime = 0.01
           endif
        endif
        !
        ! We cycle over the passive spesies only picking them by means of the indPass_transp mapping
        !
        if (fLifeTime > 0) then
           do iSpecies = 1, nPass_transp
             vMass(indPass_transp(iSpecies)) = vMass(indPass_transp(iSpecies))  / (1. + timestep_sec / fLifeTime) 
           end do
        endif
      endif  ! self-decaying species exist
      !
      ! Time tag treatment
      !
      if(rulesPassive%ifTimeTracer)then
        if(iz == 1) then
          vMass(indPassTime_transp) = vMass(indPassOnes_transp)*timestep_sec
        else
          vMass(indPassTime_transp) = vMass(indPassTime_transp) + vMass(indPassOnes_transp)*timestep_sec
        endif
      endif

      if(rulesPassive%ifSF6destroy .and. iz == nz_dispersion) then
          ! Too straightforward -- no destruction here, but rather diffusion+ destruction
          ! should be very similar to deposition velocity 
          !    ! 1yr at 20 pa, 100 yrs at 100 pa 
          !    ! Very vague reading of Fig 9 in Totterdill et al. Journal of Physical Chemistry A 2015 v 119 pp 2016-2025
          !    ! doi: 10.1021/jp5123344
          !    fLifeTime = 3e7 * (metdat(ind_pressure) / 20 )**3

          ! The destruction is driven by resistance to the layer of ~0.8 Pa, thus relying on normal diffusion
          ! except for the uppermost level
          ! The assuming profile of Kz from Hunten 1974 (very vague table 3 from Massie&Hunten 1981) Kz = 80*(100Pa/p)**.7
          ! i.e.. 15 m/s2 @ 1000Pa and 80 @ 100 Pa
          ! Assuming for mesosphere T=275 K
          ! R/rho = 2.3e5 * int (dp/p**1.25) 

          fTmp = 1./(1. + timestep_sec/rulesPassive%SF6LifetimeUpperLayer)

          if (indSF6_transp > 0) &
            & vMass(indSF6_transp) = vMass(indSF6_transp) * fTmp
          if (indSF6ng_transp > 0) &
            & vMass(indSF6ng_transp) = vMass(indSF6ng_transp) * fTmp
      endif

    endif ! defined rules

  end subroutine transform_passive



  !*******************************************************************************
  !*******************************************************************************
  !
  !  Passive rules
  !
  !*******************************************************************************
  !*******************************************************************************

  !***********************************************************************

  subroutine set_chem_rules_passive(nlSetup, rulesPassive)
    !
    ! Creates the set of rules, in accordance with which the computations are to
    ! be made. Input - SILAM namelist
    !
    implicit none

    ! Imported parameter
    type(Tsilam_namelist), intent(in) :: nlSetup 
    type(Tchem_rules_passive), intent(out) :: rulesPassive

    rulesPassive%defined = silja_false

    if(.not.defined(nlSetup))then
      call set_error('Undefined namelist given','set_chem_rules_passive')
      return
    endif
    !
    ! If basic temperature and corresponding life time are not defined - set the whole story 
    ! to non-existent.
    !
    if(fu_content(nlSetup,'passive_subst_ref_lifetime') == '' .or. &
     & fu_content(nlSetup,'passive_subst_ref_tempr') == '')then
      rulesPassive%basicLifeTime = real_missing
      rulesPassive%basicTempr = 273.15
      rulesPassive%DLifeTime_DT = 0.0
    endif
    !
    ! Set the reference interval and temperature
    !
    rulesPassive%basicLifeTime = &
             & fu_sec(fu_set_named_interval(fu_content(nlSetup,'passive_subst_ref_lifetime')))
    if(error .or. rulesPassive%basicLifeTime < 0.001)return

    rulesPassive%basicTempr = fu_content_real(nlSetup,'passive_subst_ref_tempr')
    if(error .or. rulesPassive%basicTempr < 0.001 .or. (rulesPassive%basicTempr .eps. real_missing))return

    !
    ! Derivative of lifetime over temperature
    ! Note that we know that this derivative is time/temperature intervals. The fu_set_named_value
    ! will return the value in the SI units, i.e., sec/degree_K.
    !
    if(fu_content(nlSetup,'passive_subst_dLifeTime_dT') /= '')then
      rulesPassive%DLifeTime_DT = &
           & fu_sec(fu_set_interval_sec(fu_set_named_value(fu_content(nlSetup,'passive_subst_dLifeTime_dT'))))
      if(error .or. rulesPassive%DLifeTime_DT < 0.001) rulesPassive%DLifeTime_DT = 0.0
    else
      rulesPassive%DLifeTime_DT = 0.0
    endif

    rulesPassive%ifInitialisedOnes = .false.
    if(fu_str_u_case(fu_content(nlSetup,'passive_ones_tracer')) == 'YES')then
      rulesPassive%ifOnesTracer = .true.
    else
      rulesPassive%ifOnesTracer = .false.
    endif
    !
    ! A special type of the passive tracer: time tracer
    !
    rulesPassive%ifInitialisedTimeTag = .false.
    if(fu_str_u_case(fu_content(nlSetup,'passive_time_tracer')) == 'YES')then
      rulesPassive%ifTimeTracer = .true.
      if (.not. rulesPassive%ifOnesTracer ) then
        call msg('passive_time_tracer requested.  please add "passive_ones_tracer = YES" to the tranformation namelist')
        call set_error("Timetag requires ones tracer to work properly", "set_chem_rules_passive")
        return
      endif
    else
      rulesPassive%ifTimeTracer = .false.
    endif


    if(fu_str_u_case(fu_content(nlSetup,'destroy_SF6')) == 'YES')then
      rulesPassive%ifSF6destroy = .true.

      !! For Hunten profile and rate of destruction from ***
          ! The assuming profile of Kz from Hunten 1974 (very vague table 3 from Massie&Hunten 1981) Kz = 80*(100Pa/p)**.7
          ! i.e.. 15 m/s2 @ 1000Pa and 80 @ 100 Pa
          ! Assuming for mesosphere T=275 K
          ! R/rho = 2.3e5 * int (dp/p**1.25) 
      !! R/rho from upper layer to destruction layer (0.8 Pa) == 9.8e5 (1.05 - p_upper_layer**(-0.25))
      !! Then lifetime = R/rho * dp_upper_layer/g
      !! For p=15Pa and dp=10Pa(aoa setup) R/rho dp/g ~= 4.9e5 s ~= 6 days.... 
      rulesPassive%SF6LifetimeUpperLayer = &
           & fu_sec(fu_set_interval_sec(fu_set_named_value(fu_content(nlSetup,'SF6LifetilemeTop'))))
      if (.not. rulesPassive%SF6LifetimeUpperLayer > 0.) then
        call set_error("destroy_SF6 requested, but no SF6LifetilemeTop specified", "set_chem_rules_passive") 
        return
      endif
    else
      rulesPassive%ifSF6destroy = .false.
    endif

    rulesPassive%defined = silja_true

  end subroutine set_chem_rules_passive


  !********************************************************************************

  function fu_lifetime_passive(rules)result(lifetime)
    !
    ! A typical life time due to degradation and deposition for radioactive materials
    ! So far, whatever the species are, we set the "typical" lifetime to be two days
    !
    implicit none

    type(silja_interval) :: lifetime
    type(Tchem_rules_passive), intent(in) :: rules

    if(rules%defined == silja_true)then
      if(rules%basicLifeTime > 0 .and. rules%basicLifeTime < 3.2e8)then
        lifetime = fu_set_interval_sec(rules%basicLifeTime)
      else
        lifetime = one_day * 365.0
      endif
    else
      lifetime = interval_missing
    endif

  end function fu_lifetime_passive


  !**********************************************************************************

  logical function fu_if_specific_dep_passive(rulesPassive)
    implicit none
    type(Tchem_rules_passive), intent(in) :: rulesPassive
    fu_if_specific_dep_passive = .false.     ! no specific deposition whatsoever
  end function fu_if_specific_dep_passive

  !************************************************************************************
  

  logical function fu_if_tla_required_passive(rules) result(required)
    ! Collect transformations' requests for tangent linearization. If tangent linear
    ! is not allowed, corresponding subroutine sets error. 
    implicit none
    type(Tchem_rules_passive), intent(in) :: rules
    
    required = .false.

  end function fu_if_tla_required_passive

  
  !************************************************************************************
  
  subroutine if_need_init_time_tag(rules, indTimeTag, ifNeedInit, ifReset)
    !
    ! Just says whether the time tracer field has been initialised - and returns the index of its species
    !
    implicit none
    
    ! Imported parameters
    type(Tchem_rules_passive), intent(inout) :: rules
    logical, intent(in) :: ifReset
    integer, intent(out) :: indTimeTag
    logical, intent(out) :: ifNeedInit
    
    ifNeedInit = rules%ifTimeTracer .and. (.not. rules%ifInitialisedTimeTag)
    
    indTimeTag = indPassTime_transp
    if(ifReset) rules%ifInitialisedTimeTag = .true.  ! somewhat dangerous but saves a function call to reset this switch
    
  end subroutine if_need_init_time_tag

  subroutine if_need_init_ones(rules, indOnes, ifNeedInit, ifReset)
    !
    ! Just says whether the time tracer field has been initialised - and returns the index of its species
    !
    implicit none
    
    ! Imported parameters
    type(Tchem_rules_passive), intent(inout) :: rules
    logical, intent(in) :: ifReset
    integer, intent(out) :: indOnes
    logical, intent(out) :: ifNeedInit
    
    ifNeedInit = rules%ifOnesTracer .and. (.not. rules%ifInitialisedOnes)
    
    indOnes = indPassOnes_transp
    if(ifReset) rules%ifInitialisedOnes = .true.  ! somewhat dangerous but saves a function call to reset this switch
    
  end subroutine if_need_init_ones
  
END MODULE chem_dep_passive

