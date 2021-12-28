MODULE aer_dyn_basic

  use cocktail_basic
  use globals, only : r4k, r8k

  implicit none
  
  !
  ! public routines for ADB
  !
  public set_rules_AerDynBasic
  public AerDynBasic_input_needs
  public init_AerDynBasic
  public transform_AerDynBasic
!  public registerSpecies
  public fu_if_tla_required

  public full_species_lst_4_ADB
!  private registerSpecies_4_ADB

  private equilibration
  private condensation
  private coagulation
  private update_distributions
  private advect_nbr_distr
  private distribute_SO4_aq

  private fu_if_tla_required_ADB
  
  interface fu_if_tla_required
     module procedure fu_if_tla_required_ADB
  end interface


!  interface registerSpecies
!    module procedure registerSpecies_4_ADB
!  end interface


  !----------------------------------------------------------------------
  !
  ! Dimensions of vols 3d matrix
  integer, private, save :: nMaterials, nBins, nParallelBins
  !
  ! Named indices for vols matrix:
  ! Named indices of regime start & end bins
  INTEGER, private, save :: in1, in2, in3, fn1, fn2, fn3
  ! Named indices of aerosol materials
  integer, private, save :: iSO4, iOC, iBC, iSslt, iDust, iNO3, iNH4
  ! Named indices of parallel bins
  integer, private, save :: soluble, insoluble
  !
  ! Material properties  
  character(len = substNmLen), dimension(:), pointer, private, save :: chADBMaterialNms
  real, dimension(:), pointer, private, save :: factor2volume, density, molarMass, volMolec
  !
  ! Aerosol modes
  type(Taerosol_mode), dimension(:,:), pointer, private, save :: modesADB
  ! Mode properties  
  real, dimension(:), pointer, save :: vhilim, vlolim, dpmid   ! as they are used, lets keep them for now
  !
  ! Aerosol species (mass and number concentrations)
  type(silam_species), dimension(:), pointer, private, save :: speciesADB_aerTransp, speciesADB_aerosol
  type(silam_material), pointer, save :: materialADB_aerosol

  integer, private, save :: nADB_aerTransp = 0, nADB_aerosol = 0
  !
  ! Indices for aerosol species in mass and number maps
  integer, dimension(:,:,:), pointer, private, save :: indADB_transp
  integer, dimension(:,:), pointer, private, save :: indADB_aerosol
  !
  ! Named indices for gaseous species in transport & shortliving maps
  INTEGER, private :: iSO4_gas_sl, iH2SO4_gas, iNH3, iHNO3, iOC_nv, iOC_sv
  !
  ! Named indices for meteo input
  integer, private, pointer, save :: ind_t, ind_p, ind_rh  

  !
  !-------------------------------------------------------------------------
  !
  !  Tchem_rules_AerDynBasic type definition
  !
  type Tchem_rules_AerDynBasic
    private

    logical :: nucl_on, coag_on, cond_on, cloud_on, recalc
    integer :: nucl_scheme ! 1 = binary nucleation
                           ! 2 = ternary nucleation
                           ! 3 = kinetic nucleation
                           ! 4 = activation nucleation
   
    INTEGER :: nreg        ! number of main size regimes
    REAL, dimension(:), pointer ::    regLims !  low/high diameter limits of main size regimes [m]      
    INTEGER, dimension(:), pointer :: regNbin ! number of bins in each main regime
    REAL(r8k), dimension(:,:), pointer :: massacc
    REAL(r8k), dimension(:,:,:), pointer :: coagtable
    real, dimension(:,:,:,:), pointer :: wetParticleGrowth 
    REAL(r8k) :: epsv    !critical dry volume ratio of soluble matter
    REAL(r8k) :: epsoc   ! water uptake of organic material
    
    type(silam_species), dimension(:), pointer :: species_SO4aq  
    type(TspeciesReference), dimension(:), pointer :: mapMass_SO4aq,  mapNbr_SO4aq   

    real :: nbrLowThreshold !  low number concentration treshold [#/m3] 
    logical :: ifSO4, ifNH4, ifNO3, ifOC, ifBC, ifSslt, ifDust  ! Materials in run

    type(silja_logical) :: defined
  end type Tchem_rules_AerDynBasic

  !
  ! Label of this module
  !
  integer, parameter, public :: aerosol_dynamics_basic = 5020



  !-------------------------------------------------------------------------

  CONTAINS

  !***********************************************************************

  subroutine set_rules_AerDynBasic(nlSetup, rulesADB)
    !
    ! Creates the set of rules, in accordance with which the computations are to
    ! be made. Input - SILAM namelist
    !
    implicit none

    ! Imported parameter
    type(Tsilam_namelist), intent(in) :: nlSetup 
    type(Tchem_rules_AerDynBasic), intent(out) :: rulesADB

    ! Stuff to be later read from ini files !!!!!!!!!!!!!!!!
    rulesADB%nreg = 3            ! number of main size regimes

    allocate(rulesADB%reglims(rulesADB%nreg+1))
    allocate(rulesADB%regNbin(rulesADB%nreg))

    rulesADB%reglims = (/ 3.e-9, 4.e-8, 1.e-6, 5.e-4 /)   ! low/high diameter limits of main size regimes [m]      
    rulesADB%regNbin = (/ 20, 200, 1 /)  ! number of bins in each main regime


    nMaterials = 7
    nBins = sum(rulesADB%regNbin(:))
    nParallelBins = 2

    rulesADB%defined = silja_false

    if(.not.defined(nlSetup))then
      call set_error('Undefined namelist given','set_chem_rules_ADB')
      return
    endif
    
    if(fu_str_u_case(fu_content(nlSetup,'ADB_if_compute_nucleation')) == 'YES')then
      rulesADB%nucl_on = .true.
      if(fu_str_u_case(fu_content(nlSetup,'ADB_nucleation_scheme')) == 'BINARY')then
        rulesADB%nucl_scheme = 1
      elseif(fu_str_u_case(fu_content(nlSetup,'ADB_nucleation_scheme')) == 'TERNARY')then
        rulesADB%nucl_scheme = 2
      elseif(fu_str_u_case(fu_content(nlSetup,'ADB_nucleation_scheme')) == 'KINETIC')then
        rulesADB%nucl_scheme = 3
      elseif(fu_str_u_case(fu_content(nlSetup,'ADB_nucleation_scheme')) == 'ACTIVATION')then
        rulesADB%nucl_scheme = 4
      else
        call set_error('Strange content for ADB_nucleation_scheme','set_chem_rules_ADB')
        return
      endif
    elseif(fu_str_u_case(fu_content(nlSetup,'ADB_if_compute_nucleation')) == 'NO')then
      rulesADB%nucl_on = .false.
    else
      call set_error('Strange content for ADB_if_compute_nucleation','set_chem_rules_ADB')
      return
    endif

    if(fu_str_u_case(fu_content(nlSetup,'ADB_if_compute_coagulation')) == 'YES')then
      rulesADB%coag_on = .true.
    elseif(fu_str_u_case(fu_content(nlSetup,'ADB_if_compute_coagulation')) == 'NO')then
      rulesADB%coag_on = .false.
    else
      call set_error('Strange content for ADB_if_compute_coagulation','set_chem_rules_ADB')
      return
    endif

    if(fu_str_u_case(fu_content(nlSetup,'ADB_if_compute_condensation')) == 'YES')then
      rulesADB%cond_on = .true.
    elseif(fu_str_u_case(fu_content(nlSetup,'ADB_if_compute_condensation')) == 'NO')then
      rulesADB%cond_on = .false.
    else
      call set_error('Strange content for ADB_if_compute_condensation','set_chem_rules_ADB')
      return
    endif

    if(fu_str_u_case(fu_content(nlSetup,'ADB_if_compute_cloud_activation')) == 'YES')then
      rulesADB%cloud_on = .true.
    elseif(fu_str_u_case(fu_content(nlSetup,'ADB_if_compute_cloud_activation')) == 'NO')then
      rulesADB%cloud_on = .false.
    else
      call set_error('Strange content for ADB_if_compute_cloud_activation','set_chem_rules_ADB')
      return
    endif

    if(fu_str_u_case(fu_content(nlSetup,'ADB_if_recalc_wet_d')) == 'YES')then
      rulesADB%recalc = .true.
    elseif(fu_str_u_case(fu_content(nlSetup,'ADB_if_recalc_wet_d')) == 'NO')then
      rulesADB%recalc = .false.
    else
      call set_error('Strange content for ADB_if_recalc_wet_d','set_chem_rules_ADB')
      return
    endif

    !
    ! For the time being ..
    !
    rulesADB%nbrLowThreshold = 0.001  ! #/m3

    rulesADB%ifSO4 = .false.
    rulesADB%ifNH4 = .false.
    rulesADB%ifNO3 = .false.
    rulesADB%ifOC = .false.
    rulesADB%ifBC = .false.
    rulesADB%ifSslt = .false.
    rulesADB%ifDust = .false.

    allocate(rulesADB%massacc(nbins,nParallelBins))    
    rulesADB%massacc(:,:) = 1.


    rulesADB%epsv = real_missing !critical dry volume ratio of soluble matter
    rulesADB%epsoc = 0.15_r8k    ! water uptake of organic material

    rulesADB%defined = silja_true


  end subroutine set_rules_AerDynBasic


  !***********************************************************************

  subroutine AerDynBasic_input_needs(rulesADB, meteo_input_local)
    !
    ! Returns input needs for the transformation routines. 
    !
    implicit none

    ! Imported parameters
    type(Tchem_rules_AerDynBasic), intent(in) :: rulesADB
    type(Tmeteo_input), intent(out), target :: meteo_input_local

    ! Local variables
    integer :: iTmp

    meteo_input_local = meteo_input_empty 

    if(.not.(rulesADB%defined == silja_true))return

    meteo_input_local%nQuantities = 3

    meteo_input_local%quantity(1) = temperature_flag    
    meteo_input_local%q_type(1) = meteo_dynamic_flag

    meteo_input_local%quantity(2) = pressure_flag    
    meteo_input_local%q_type(2) = meteo_dynamic_flag

    meteo_input_local%quantity(3) = relative_humidity_flag    
    meteo_input_local%q_type(3) = meteo_dynamic_flag

    !
    ! Finally, prepare the named indices, which will be used in the actual calculations
    !
    ind_t => meteo_input_local%idx(1)
    ind_p => meteo_input_local%idx(2)
    ind_rh => meteo_input_local%idx(3)

  end subroutine AerDynBasic_input_needs


  !***********************************************************************

  subroutine init_AerDynBasic(rules, species_lst_initial)

    ! Set ADB modes, transport- and aerosol species arrays 

    implicit none

    ! Imported parameter
    type(silam_species), dimension(:), allocatable, intent(out) :: species_lst_initial
    type(Tchem_rules_AerDynBasic), intent(inout) :: rules
    ! local variables 
    integer :: i, iMode
    real ::  ratio, dlo, dhi, dmid
    character(len=substNmLen) :: modeNm
    type(silam_material), pointer :: materialPtrTmp

    !-- local variables ------
    REAL(r8k) :: aa, bb   ! constants in Kohler equation
    REAL(r8k), parameter :: ions = 3.0_r8k,   & ! van't Hoff factor (ions produced upon dissociation)
                           slim = 1.005_r8k    ! water saturation used as limit                               

    integer :: nRH, nT, iRH, iT, klev
    real, dimension(:), pointer :: rh, t
    type(TwetParticle) :: wetP

    REAL, dimension(:), pointer :: pres0, & ! pressure at each vertical level [Pa]
                                   temp0   ! temperature at each vertical level [K]

    INTEGER :: ii, jj, kk  ! loop indices

    REAL(r8k) :: mpart(nBins)   ! approximate mass of particles [kg]

    real :: ftmp


    ! Named indices
    in1 = 1                           ! regime 1
    in2 = in1 + rules%regNbin(1)      ! regime 2
    in3 = in2 + rules%regNbin(2)      ! regime 3
                        
    ! number/radius: last index
    fn1 = in2 - 1                    ! regime 1a
    fn2 = in3 - 1                    ! regime 2b
    fn3 = in3 + rules%regNbin(3) - 1 ! regime 3b

    allocate(chADBMaterialNms(nMaterials))
    allocate(density(nMaterials))
    allocate(molarMass(nMaterials))
    allocate(volMolec(nMaterials))
    allocate(factor2volume(nMaterials))
    chADBMaterialNms = (/'SO4 ', 'OC  ', 'BC  ', 'sslt', 'dust', '_NO3', 'NH4 '/)
  
    iSO4 = 1
    iOC = 2   
    iBC = 3   
    iSslt = 4 
    iDust = 5 
    iNO3 = 6  
    iNH4 = 7

    soluble = 1
    insoluble = 2


    ! Compute the factors for converting silam mass map to ADB vols matrix
    ! Have to take into account, that some things are in moles and some in kg-s
    do i = 1, nMaterials
      materialPtrTmp => fu_get_material_ptr(chADBMaterialNms(i))
      if(error)return
      density(i) = fu_dry_part_density(materialPtrTmp)
      molarMass(i) = fu_mole_mass(materialPtrTmp)
      volMolec(i) = molarMass(i)/avogadro/density(i)
      factor2volume(i) = fu_conversion_factor(fu_basic_mass_unit(materialPtrTmp), &
                                            & 'kg', materialPtrTmp) / &
                       & fu_dry_part_density(materialPtrTmp)
    enddo


    !
    ! First fill the modes array (20 bins, sometimes parallel)
    !-------------------------------------------------------------------------------
    !-- 1) size regime 1: --------------------------------------
    ! 1 parallel
    !

    allocate(modesADB(nBins,nParallelBins))
    allocate(vlolim(nBins))
    allocate(vhilim(nBins))
    allocate(dpmid(nBins))
!    allocate(speciesADB_aerosol(nBins*2-fn1))


    ratio = rules%reglims(2)/rules%reglims(1)   ! section spacing

    DO i = in1,fn1
       dlo = rules%reglims(1)*ratio**(real(i-1)/rules%regNbin(1))
       dhi = rules%reglims(1)*ratio**(real(i)/rules%regNbin(1))
       dpmid(i) = (0.5 * (dhi**3 + dlo**3))**(1./3.)
       modesADB(i,soluble) = fu_set_mode(moving_diameter_flag, dlo, dhi, &
                                 & label = fu_connect_strings('1a', fu_aerosol_mode_size_to_str(dpmid(i))), &
                                 & solubility = 1)
       ! ADB wants these ..
       vlolim(i) = pi / 6. * dlo**3
       vhilim(i) = pi / 6. * dhi**3

    END DO


    !-- 2) size regime 2: --------------------------------------
    ! 2 parallel
    !
    ratio = rules%reglims(3)/rules%reglims(2)   ! section spacing

    DO i = in2,fn2
       dlo = rules%reglims(2)*ratio**(real(i - in2)/rules%regNbin(2))
       dhi = rules%reglims(2)*ratio**(real(i - in2 + 1)/rules%regNbin(2))
       dmid = (0.5*(dhi**3 + dlo**3))**(1./3.)
       modesADB(i,soluble) = fu_set_mode(moving_diameter_flag, dlo, dhi, &
                                 & label = fu_connect_strings('2a', fu_aerosol_mode_size_to_str(dpmid(i))), &
                                 & solubility = 1)
       modesADB(i,insoluble) = fu_set_mode(moving_diameter_flag, dlo, dhi, &
                                 & label = fu_connect_strings('2b', fu_aerosol_mode_size_to_str(dpmid(i))), &
                                 & solubility = 0)
       ! ADB wants these ..
       vlolim(i) = pi/6 * dlo**3
       vhilim(i) = pi/6 * dhi**3
       dpmid(i) = (0.5*(dhi**3 + dlo**3))**(1./3.)

    END DO

    !-- 3) size regime 3: --------------------------------------
    ! 2 parallel
    !
    ratio = rules%reglims(4)/rules%reglims(3)   ! section spacing

    DO i = in3, fn3
       dlo = rules%reglims(3)*ratio**(real(i - in3)/rules%regNbin(3))
       dhi = rules%reglims(3)*ratio**(real(i - in3 + 1)/rules%regNbin(3))
       dpmid(i) = (0.5*(dhi**3 + dlo**3))**(1./3.)
       modesADB(i,soluble) = fu_set_mode(moving_diameter_flag, dlo, dhi, &
                                 & label = fu_connect_strings('3a', fu_aerosol_mode_size_to_str(dpmid(i))), &
                                 & solubility = 1)
       modesADB(i,insoluble) = fu_set_mode(moving_diameter_flag, dlo, dhi, &
                                 & label = fu_connect_strings('3b', fu_aerosol_mode_size_to_str(dpmid(i))), &
                                 & solubility = 0)
       ! ADB wants these ..
       vlolim(i) = pi/6 * dlo**3
       vhilim(i) = pi/6 * dhi**3
    END DO

    !
    ! Since ADB wants the aerosol features to be strict, it has to fill them in into the 
    ! chemical_aerosol in the chemical rules. That aerosol will be later used for mapping, 
    ! setting source terms with continuous size spectrum, etc.
    !

    allocate(species_lst_initial(nMaterials*nBins*nParallelBins))
    species_lst_initial(:) = species_missing
   
    i = 0
    materialPtrTmp => fu_get_material_ptr('SO4')
    do iMode = in1,fn3
      i = i + 1
      call set_species(species_lst_initial(i), materialPtrTmp, modesADB(iMode,soluble))
    enddo  
    do iMode = in2,fn3
      i = i + 1
      call set_species(species_lst_initial(i), materialPtrTmp, modesADB(iMode,insoluble))
    enddo 

    materialPtrTmp => fu_get_material_ptr('_NO3')
    do iMode = in1,fn3
      i = i + 1
      call set_species(species_lst_initial(i), materialPtrTmp, modesADB(iMode,soluble))
    enddo  
    do iMode = in2,fn3
      i = i + 1
      call set_species(species_lst_initial(i), materialPtrTmp, modesADB(iMode,insoluble))
    enddo  
 
    materialPtrTmp => fu_get_material_ptr('NH4')
    do iMode = in1,fn3
      i = i + 1
      call set_species(species_lst_initial(i), materialPtrTmp, modesADB(iMode,soluble))
    enddo  
    do iMode = in2,fn3
      i = i + 1
      call set_species(species_lst_initial(i), materialPtrTmp, modesADB(iMode,insoluble))
    enddo  

    materialPtrTmp => fu_get_material_ptr('OC')
    do iMode = in1,fn3
      i = i + 1
      call set_species(species_lst_initial(i), materialPtrTmp, modesADB(iMode,soluble))
    enddo 
    do iMode = in2,fn3
      i = i + 1
      call set_species(species_lst_initial(i), materialPtrTmp, modesADB(iMode,insoluble))
    enddo 

    materialPtrTmp => fu_get_material_ptr('BC')
    do iMode = in2,fn3
      i = i + 1
      call set_species(species_lst_initial(i), materialPtrTmp, modesADB(iMode,soluble))
    enddo 
    do iMode = in2,fn3
      i = i + 1
      call set_species(species_lst_initial(i), materialPtrTmp, modesADB(iMode,insoluble))
    enddo 

    materialPtrTmp => fu_get_material_ptr('sslt')
    do iMode = in2,fn2
      i = i + 1
      call set_species(species_lst_initial(i), materialPtrTmp, modesADB(iMode,soluble))
    enddo 
    do iMode = in2,fn3
      i = i + 1
      call set_species(species_lst_initial(i), materialPtrTmp, modesADB(iMode,insoluble))
    enddo 

    materialPtrTmp => fu_get_material_ptr('dust')
    do iMode = in2,fn3
      i = i + 1
      call set_species(species_lst_initial(i), materialPtrTmp, modesADB(iMode,soluble))
    enddo 
    do iMode = in2,fn3
      i = i + 1
      call set_species(species_lst_initial(i), materialPtrTmp, modesADB(iMode,insoluble))
    enddo 

    allocate(speciesADB_aerTransp(i))
    speciesADB_aerTransp(1:i) = species_lst_initial(1:i)

call set_error('The number concentrations now in transported mass map, linked aerosol masses in aerosol mass map', &
             & 'init_AerDynBasic')
! ATTENTION.
! The idea is: all that is transported is in transported mass map. 
! Short-lived and aerosol species are not advected and thus have no moments. Instead,
! they are linked to specific "hosts" in the transport mass map, which control their advection.
! In particular, number concentrations are advected, whereas the masses constituting them are
! in aerosol mass map, linked and advected as "passengers".
! Below code seems to follow different approach: number concentrations are in aerosol mass map.
! Should be changed.
!

    allocate(materialADB_aerosol)
    call set_material_missing(materialADB_aerosol)
    materialPtrTmp => materialADB_aerosol
    allocate(speciesADB_aerosol(nBins*nParallelBins-fn1))

    do i = 1, nBins
      call set_species(speciesADB_aerosol(i), materialPtrTmp, modesADB(i,soluble))
    enddo
    do i = in2, nBins
      call set_species(speciesADB_aerosol(nBins+i-fn1), materialPtrTmp, modesADB(i,insoluble))
    enddo

    !-- Calculation of critical soluble material 
    !   (must be divided with temp**3 when used) ------------------
    aa = 4. * molecular_weight_water * srfTns_water / (gas_constant_uni * density_water)   ! curvature (Kelvin) effect
    bb = 6. * molecular_weight_water * ions * density(iSO4) / (pi * density_water * molarMass(iSO4)) ! solute (Raoult) effect
    rules%epsv = 4._r8k * aa**3 / (27._r8k * bb * (log(slim))**2)

! Fill the particle RH growth table
    RH => fu_work_array()
    T => fu_work_array()
    nRH = 34
    nT = 11
    RH(1:nRH) = (/0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.54, 0.58, 0.62, 0.66, &
               & 0.70, 0.73, 0.76, 0.79, 0.82, 0.85, 0.87, 0.89, 0.91, 0.93, 0.95, 0.96, 0.97, 0.98, 0.99, &
               & 100., 101., 102., 103./)
    do iT = 1, nT
      T(iT) = 183+15*iT
    enddo
    allocate(rules%wetParticleGrowth(nMaterials, nBins, nRH, nT))
    do i = 1, nMaterials
      do iMode = 1, nBins 
        do iRH = 1, nRH ! Relative humidity 
          do iT = 1, nT
            wetP = fu_wet_particle_features(fu_get_material_ptr(chADBMaterialNms(i)), fu_massmean_D(modesADB(iMode, 1)), T(iT), RH(iRH))
            if(error)return
            rules%wetParticleGrowth(i, iMode, iRH, iT) = wetP%fGrowthFactor
          enddo
        enddo
      enddo
    enddo

    call free_work_array(RH)
    call free_work_array(T)

!    ! Fill the coagulation coeffitients table 
!    if(.not.(fu_leveltype(dispersion_vertical) == layer_btw_2_height .or. &
!           & fu_leveltype(dispersion_vertical) == layer_btw_2_h_altit_above_msl))then
!      call msg('Cannot compute coagulation coeffitients for dispersion vertical type', &
!              & fu_leveltype(dispersion_vertical))
!      call set_error('Cannot compute coagulation coeffitients','set_coagc')
!      return
!    endif
!    klev = fu_NbrOfLevels(dispersion_vertical)
!    pres0 => fu_work_array()
!    temp0 => fu_work_array()
!
!    allocate(rules%coagtable(nBins,nBins,klev))
!    rules%coagtable(:,:,:) = 0.
!
!    do ii = 1, klev
!      call us_standard_atmosphere(fu_layer_centre_value(dispersion_verticalPtr, ii), &
!                                & ftmp, pres0(ii), temp0(ii)) 
!      pres0(ii) = pres0(ii) * 288.15
!      temp0(ii) = temp0(ii) * std_pressure_sl
!    enddo
!  
!    !-- particle mass; density of 1500 kg/m3 assumed [kg] 
!    mpart = (pi/6)*dpmid(1:nBins)**3 * 1500.
!
!    !-- calculation of coagulation coefficients [m3/s]
!    DO kk = 1,nBins
!      DO jj = 1,nBins
!          DO ii = 1,klev 
!             rules%coagtable(ii,jj,kk) = coagc(dpmid(jj),dpmid(kk),mpart(jj),mpart(kk),temp0(ii),pres0(ii))
!          END DO
!       END DO
!    END DO
!
!    call free_work_array(pres0)
!    call free_work_array(temp0)


    contains

!====================================================================================================

      FUNCTION coagc(diam1,diam2,mass1,mass2,temp,pres)

      !
      ! Purpose:
      ! --------
      ! Calculates the coagulation coefficient for two
      ! colliding particles. 
      !
      !
      ! Method: 
      ! -------  
      ! Only Brownian coagulation taken into account.
      ! Transition regime correction is done with Fuchs  
      ! flux matching.
      !
      !
      ! Interface:
      ! ----------
      ! Called from subroutine SET_COAG
      ! (only at the beginning of simulation)
      !
      !
      ! Coded by:
      ! ---------
      ! Hannele Korhonen (FMI) 2005 
      !

        IMPLICIT NONE

        !-- Input variables ----------
        REAL, INTENT(IN) :: diam1, diam2   ! diameters of colliding particles [m]
        REAL(r8k), INTENT(IN) :: mass1, mass2   ! masses -"- [kg]
        
        real, intent(in) ::     temp,   &   ! ambient temperature [K]
                                pres        ! ambient pressure [fxm]

        !-- Output variables ---------
        REAL(r8k) :: coagc       ! coagulation coefficient of particles [m3/s]

        !-- Local variables ----------  
        REAL(r8k) :: visc,   &   ! viscosity of air [kg/(m s)]
                    mfp,    &   ! mean free path of air molecules [m]
                    mdiam,  &   ! mean diameter of colliding particles [m]
                    fmdist      ! distance of flux matching [m]

        REAL(r8k), DIMENSION (2) :: diam,   &   ! diameters of particles [m]
                                   mpart,  &   ! masses of particles [kg]
                                   knud,   &   ! particle knudsen number [1]
                                   beta,   &   ! Cunningham correction factor [1]
                                   dfpart, &   ! particle diffusion coefficient [m2/s]
                                   mtvel,  &   ! particle mean thermal velocity [m/s]
                                   omega,  &   !
                                   tva,    &   ! temporary variable [m]
                                   flux        ! flux in continuum and free molec. regime [m/s]

        !-- 0) Initializing particle and ambient air variables --------------------
        diam = (/ diam1, diam2 /)       ! particle diameters [m]
        mpart = (/ mass1, mass2 /)       ! particle masses [kg]
        visc = (7.44523e-3_r8k*temp**1.5_r8k)/(5093._r8k*(temp+110.4_r8k)) ! viscosity of air [kg/(m s)]
        mfp = (1.656e-10_r8k*temp+1.828e-8_r8k)*std_pressure_sl/pres ! mean free path of air [m]

        !-- 2) Slip correction factor for small particles -------------------------
        knud = 2._r8k*mfp/diam                                    ! Knudsen number
        beta = 1._r8k+knud*(1.142_r8k+0.558_r8k*exp(-0.999_r8k/knud))! Cunningham correction factor
        ! (Allen and Raabe, Aerosol Sci. Tech. 4, 269)

        !-- 3) Particle properties ------------------------------------------------
        dfpart = beta*boltzmann_const*temp/(3._r8k*pi*visc*diam)  ! diffusion coefficient [m2/s]
        mtvel = sqrt((8._r8k*boltzmann_const*temp)/(pi*mpart))    ! mean thermal velocity [m/s]
        omega = 8._r8k*dfpart/(pi*mtvel)
        mdiam = 0.5_r8k*(diam(1)+diam(2))               ! mean diameter [m]

        !-- 4) Calculation of fluxes and flux matching ----------------------------
        flux(1) = 4._r8k*pi*mdiam*(dfpart(1)+dfpart(2)  )    ! flux in continuum regime [m3/s]
        flux(2) = pi*sqrt(mtvel(1)**2+mtvel(2)**2)*mdiam**2 !  -"- in free molec. regime [m3/s]
        ! temporary variable [m]
        tva(1) = ((mdiam+omega(1))**3 - (mdiam**2+omega(1)**2)* &
                  sqrt((mdiam**2+omega(1)**2)))/(3._r8k*mdiam*omega(1)) - mdiam
        ! temporary variable [m]
        tva(2) = ((mdiam+omega(2))**3 - (mdiam**2+omega(2)**2)* &
                 sqrt((mdiam**2+omega(2)**2)))/(3._r8k*mdiam*omega(2)) - mdiam
        ! flux matching distance [m]
        fmdist = sqrt(tva(1)**2+tva(2)**2)             

        !-- 5) Coagulation coefficient [m3/s] -------------------------------------
        coagc = flux(1) / (mdiam/(mdiam+fmdist) + flux(1)/flux(2)) 

      END FUNCTION coagc


  end subroutine init_AerDynBasic


  !***********************************************************************

  subroutine full_species_lst_4_ADB(rulesADB, &
                                    & speciesEmis, speciesTransp, speciesShortLived, speciesAerosol, &
                                    & nSpeciesEmis, nSpeciesTransp, nSpeciesShortLived, nSpeciesAerosol, &
                                    & iClaimedSpecies)
    !
    ! Here we select the species to be transported and short-lived having the 
    ! list of emitted species and the list of available ones (the SILAM database).
    ! Also, an option for fixed order can be executed here.
    !
    implicit none

    ! Imported parameters
    type(Tchem_rules_AerDynBasic), intent(inout) :: rulesADB
    type(silam_species), dimension(:), pointer :: speciesEmis, speciesTransp, speciesShortLived, &
                                                & speciesAerosol
    integer, intent(in) :: nSpeciesEmis
    integer, intent(inout) :: nSpeciesTransp, nSpeciesShortLived, nSpeciesAerosol
    integer, dimension(:), intent(inout) :: iClaimedSpecies

    ! Local variables
    type(silam_species) :: speciesTmp
    integer :: iTmp, nSelected, iMat, iBin, iS
    integer, dimension(:), pointer :: indices
    type(silam_material), pointer :: materialPtrTmp

    !
    ! ADB has no right to request any gaseous species to be added to the run. It will just check 
    ! for the available ones and decide, if that's enough. If not, error, else compile the 
    ! necessary subset of its own species.
    !
    ! Not all species are in transport list - NH3 has to be checked in emission and SO4 in 
    ! shortliving list, primary particles in emission
    !
    ! So now our feeble attempt to figure out, what's around ..
    ! First gaseous part: 
    !
    call set_species(speciesTmp, fu_get_material_ptr('SO4'), in_gas_phase)
    if(fu_index(speciesTmp, speciesShortLived, nSpeciesShortLived) /= int_missing) rulesADB%ifSO4 = .true.

    call set_species(speciesTmp, fu_get_material_ptr('H2SO4'), in_gas_phase)
    if(fu_index(speciesTmp, speciesTransp, nSpeciesTransp) /= int_missing) rulesADB%ifSO4 = .true.

    call set_species(speciesTmp, fu_get_material_ptr('NH3'), in_gas_phase)
    if(fu_index(speciesTmp, speciesEmis, nSpeciesEmis) /= int_missing)then
      call claim_emission(fu_index(speciesTmp, speciesEmis, nSpeciesEmis))
      if(error)return
      rulesADB%ifNH4 = .true.
    endif

    call set_species(speciesTmp, fu_get_material_ptr('HNO3'), in_gas_phase)
    if(fu_index(speciesTmp, speciesTransp, nSpeciesTransp) /= int_missing) rulesADB%ifNO3 = .true.

    call set_species(speciesTmp, fu_get_material_ptr('OC_nv'), in_gas_phase)
    if(fu_index(speciesTmp, speciesTransp, nSpeciesTransp) /= int_missing) rulesADB%ifOC = .true.

    call set_species(speciesTmp, fu_get_material_ptr('OC_sv'), in_gas_phase)
    if(fu_index(speciesTmp, speciesTransp, nSpeciesTransp) /= int_missing) rulesADB%ifOC = .true.

    !
    ! And primary aerosol emissions in emission map: sea-salt, dust, black carbon, organic carbon, SO4)
    ! Here we don't exactly know, what we're looking for (emission mode can be anything)
    do iTmp = 1, nSpeciesEmis
     
      if(fu_substance_name(speciesEmis(iTmp)) == 'sslt')then
        call claim_emission(iTmp)
        rulesADB%ifSslt = .true.
      endif
      if(fu_substance_name(speciesEmis(iTmp)) == 'dust')then
        call claim_emission(iTmp)
        rulesADB%ifDust = .true.
      endif
      if(fu_substance_name(speciesEmis(iTmp)) == 'BC')then
        call claim_emission(iTmp)
        rulesADB%ifBC = .true.
      endif
      if(fu_substance_name(speciesEmis(iTmp)) == 'OC')then
        call claim_emission(iTmp)
        rulesADB%ifOC = .true.
      endif
      if(fu_substance_name(speciesEmis(iTmp)) == 'SO4')then
        call claim_emission(iTmp)
        rulesADB%ifSO4 = .true.
      endif
      if(error)return
    enddo

    !
    ! So now, knowing what's around, set the requested species
    !

    indices => fu_work_int_array()

    if(rulesADB%ifSO4)then
      call select_species(speciesADB_aerTransp, size(speciesADB_aerTransp), 'SO4', aerosol_mode_missing, real_missing, indices, nSelected)
      do iTmp = 1, nSelected 
        call addSpecies(speciesTransp, nSpeciesTransp,  (/speciesADB_aerTransp(indices(iTmp))/), 1, .true.)
      enddo
    endif

    if(rulesADB%ifOC)then
      call select_species(speciesADB_aerTransp, size(speciesADB_aerTransp), 'OC', aerosol_mode_missing, real_missing, indices, nSelected)
      do iTmp = 1, nSelected 
        call addSpecies(speciesTransp, nSpeciesTransp,  (/speciesADB_aerTransp(indices(iTmp))/), 1, .true.)
      enddo
    endif

    if(rulesADB%ifBC)then
      call select_species(speciesADB_aerTransp, size(speciesADB_aerTransp), 'BC', aerosol_mode_missing, real_missing, indices, nSelected)
      do iTmp = 1, nSelected 
        call addSpecies(speciesTransp, nSpeciesTransp,  (/speciesADB_aerTransp(indices(iTmp))/), 1, .true.)
      enddo
    endif

    if(rulesADB%ifSslt)then
      call select_species(speciesADB_aerTransp, size(speciesADB_aerTransp), 'sslt', aerosol_mode_missing, real_missing, indices, nSelected)
      do iTmp = 1, nSelected 
        call addSpecies(speciesTransp, nSpeciesTransp,  (/speciesADB_aerTransp(indices(iTmp))/), 1, .true.)
      enddo
    endif

    if(rulesADB%ifDust)then
      call select_species(speciesADB_aerTransp, size(speciesADB_aerTransp), 'dust', aerosol_mode_missing, real_missing, indices, nSelected)
      do iTmp = 1, nSelected 
        call addSpecies(speciesTransp, nSpeciesTransp,  (/speciesADB_aerTransp(indices(iTmp))/), 1, .true.)
      enddo
    endif
    
    if(rulesADB%ifNO3)then
      call select_species(speciesADB_aerTransp, size(speciesADB_aerTransp), '_NO3', aerosol_mode_missing, real_missing, indices, nSelected)
      do iTmp = 1, nSelected 
        call addSpecies(speciesTransp, nSpeciesTransp,  (/speciesADB_aerTransp(indices(iTmp))/), 1, .true.)
      enddo
    endif

    if(rulesADB%ifNH4)then
      call select_species(speciesADB_aerTransp, size(speciesADB_aerTransp), 'NH4', &
                        & aerosol_mode_missing, real_missing, indices, nSelected)
      do iTmp = 1, nSelected 
        call addSpecies(speciesTransp, nSpeciesTransp,  (/speciesADB_aerTransp(indices(iTmp))/), 1, .true.)
      enddo
    endif

    call free_work_array(indices)

    ! If NH3 in emission
    if(rulesADB%ifNH4)then
      call set_species(speciesTmp, fu_get_material_ptr('NH3'), in_gas_phase)
      call addSpecies(speciesTransp, nSpeciesTransp,  (/speciesTmp/), 1)
    endif

    !
    ! Request the number concentrations
    !
    call addSpecies(speciesAerosol, nSpeciesAerosol,  speciesADB_aerosol, size(speciesADB_aerosol))



    !
    ! Registration.
    ! This is stupid of course, but as ADB wants everything in its own dp arrays
    ! indexed in its prefered way, it ADB still deals with the big full list of its species,
    ! not just the ones available in the current run. The ones not existing are kept as zeroes,
    ! indices to mass maps are set to int_missing
    !

    ! Gaseous species
    call set_species(speciesTmp, fu_get_material_ptr('SO4'), in_gas_phase)
    iSO4_gas_sl = fu_index(speciesTmp, speciesShortLived, nSpeciesShortLived) 

    call set_species(speciesTmp, fu_get_material_ptr('H2SO4'), in_gas_phase)
    iH2SO4_gas = fu_index(speciesTmp, speciesTransp, nSpeciesTransp) 

    call set_species(speciesTmp, fu_get_material_ptr('NH3'), in_gas_phase)
    iNH3 = fu_index(speciesTmp, speciesTransp, nSpeciesTransp) 

    call set_species(speciesTmp, fu_get_material_ptr('HNO3'), in_gas_phase)
    iHNO3 = fu_index(speciesTmp, speciesTransp, nSpeciesTransp)

    call set_species(speciesTmp, fu_get_material_ptr('OC_nv'), in_gas_phase)
    iOC_nv = fu_index(speciesTmp, speciesTransp, nSpeciesTransp)

    call set_species(speciesTmp, fu_get_material_ptr('OC_sv'), in_gas_phase)
    iOC_sv = fu_index(speciesTmp, speciesTransp, nSpeciesTransp)

    ! Aerosols in transport map
    ! 2d indices for vols matrix 
    allocate(indADB_transp(nMaterials,nBins,nParallelBins))
    do iMat = 1, nMaterials
      do iBin = 1, nBins
        do iS = 1, nParallelBins
          call set_species(speciesTmp, fu_get_material_ptr(chADBMaterialNms(iMat)), modesADB(iBin, iS))
          call report(speciesTmp)
          indADB_transp(iMat,iBin,iS) = fu_index(speciesTmp, speciesTransp, nSpeciesTransp)
          call msg('i', indADB_transp(iMat,iBin,iS))
        enddo
      enddo
    enddo

    ! Aerosol number concentrations
    allocate(indADB_aerosol(nBins,nParallelBins))
    materialPtrTmp => materialADB_aerosol
    do iBin = 1, nBins
      do iS = 1, nParallelBins
        call set_species(speciesTmp, materialPtrTmp, modesADB(iBin,iS))
        indADB_aerosol(iBin,iS) = fu_index(speciesTmp, speciesAerosol, nSpeciesAerosol)
        print *, indADB_aerosol(iBin,iS) 
      enddo
    enddo

    ! Species that get the SO4 created in aqueous phase
    allocate(rulesADB%species_SO4aq(1))
    call set_species(rulesADB%species_SO4aq(1), fu_get_material_ptr('SO4'), &
                   & fu_set_mode(lognormal_flag, 1.0e-6, 1.5, solubility = 1))

    call create_mode_projection(rulesADB%species_SO4aq, 1, &   
                              & speciesTransp, nSpeciesTransp, &   
                              & speciesAerosol, nSpeciesAerosol, &     
                              & rulesADB%mapMass_SO4aq, &   ! mapping for mass 
                              & rulesADB%mapNbr_SO4aq, &    ! mapping for number 
                              & .true.)    ! if full spectrum


    CONTAINS

      !===================================================================================
      subroutine claim_emission(iEmis)
        implicit none
        integer, intent(in) :: iEmis
        if(iClaimedSpecies(iEmis) < 0)then
          iClaimedSpecies(iEmis) = aerosol_dynamics_basic
        else
          call msg('This species is already claimed:', iEmis, iClaimedSpecies(iEmis))
          call set_error('The emission species is already claimed', 'claim_emission')
        endif
      end subroutine claim_emission

  endsubroutine full_species_lst_4_ADB


  !*******************************************************************

  subroutine transform_AerDynBasic(rulesADB, vMass, vMassShortliving, vNumber, massLowTrsh, garbage, &
                           & mapVolume2NumberCnc, metdat, timestep_sec, print_it) 
    !
    ! Calls all aerosol dynamics processes (coagulation, condensation + nucleation); 
    ! equilibration & distribution update. Handles copying the data between mass map
    ! and its own structures
    !
    implicit none

    ! Imported parameters
    real, dimension(:), intent(inout) :: vMass, vMassShortliving, vNumber, garbage
    real, dimension(:), intent(in) :: massLowTrsh
    real, dimension(:), intent(in) :: metdat
    type(Tchem_rules_AerDynBasic), intent(in) :: rulesADB
    type(Tmoment_mapping), pointer :: mapVolume2NumberCnc
    real, intent(in) :: timestep_sec
    logical, intent(out) :: print_it

    ! local variables
    integer :: iMat, iBin, iS
    integer :: iTmp, iTmp2, nSmallSteps
    real :: d

    ! ADB 
    
    REAL(r8k) :: vols(nMaterials, nBins, nParallelBins), &
                naero(nBins, nParallelBins),   &
                core(nBins, nParallelBins),   &
                dwet(nBins, nParallelBins),  &
                v0(nBins, nParallelBins)

    print_it = .false.  ! set to true and chemistry manager will give complete dump for this cell

    ! Organize mass concentrations of aerosol chem. compounds into a matrix of volumes.
    ! Also some checking here, that the volume-number ratio is reasonable, as noone outside this module
    ! should be able to change the particle size to be outside the bin.
    !
    vols = 0._r8k
    core = 0._r8k
    naero = 0._r8k

    !
    ! Below all goes per-cubic-metre, so we have to convert all masses to concentrations
    !
!    vMass(:) = vMass(:) / fCellVolume
!    vMassShortliving(:) = vMassShortliving(:) / fCellVolume
!    vNumber(:) = vNumber(:) / fCellVolume

    do iBin = 1,nBins 
call msg('d', fu_massmean_D(modesADB(iBin, 1)))
!call msg('d', fu_massmean_D(modesADB(iBin, 2)))
!call msg('log_d_bin', log10(fu_max_D(modesADB(iBin, 1))*1e6) - log10(fu_min_D(modesADB(iBin, 1))*1e6))
      do iS = 1, nParallelBins 
        if(indADB_aerosol(iBin,iS) == int_missing)cycle
        if(vNumber(indADB_aerosol(iBin,iS)) > rulesADB%nbrLowThreshold)then
          !
          ! Reasonable number cnc. Check for mass and proceed
          !
          naero(iBin,iS) = vNumber(indADB_aerosol(iBin,iS))
          do iMat = 1, nMaterials
            if(indADB_transp(iMat,iBin,iS) /= int_missing)then
              if(vMass(indADB_transp(iMat,iBin,iS)) > massLowTrsh(indADB_transp(iMat,iBin,iS)))then
                vols(iMat,iBin,iS) =  vMass(indADB_transp(iMat,iBin,iS)) * factor2volume(iMat)
              else
                call msg('Garbage, iBin, vMass:',iBin,vMass(indADB_transp(iMat,iBin,iS)))
                garbage(indADB_transp(iMat,iBin,iS)) = garbage(indADB_transp(iMat,iBin,iS)) + &
                                                       & vMass(indADB_transp(iMat,iBin,iS))
                vMass(indADB_transp(iMat,iBin,iS)) = 0.0
              endif
            endif  ! material-bin-solubility combination exists
          enddo   ! materials
          !
          ! Mean volume of a single particle in the bin
          !
          core(iBin,iS) = sum(vols(1:nMaterials,iBin,iS)) / naero(iBin, iS)
          !
          ! Check for the its diameter
          !
          d =(6. * core(iBin,iS) / pi)**(1./3.)
          if(d < fu_min_D(modesADB(iBin, iS)) .or. d > fu_max_D(modesADB(iBin, iS)))then
            call set_error('Particle diameter at the beginning is outside bin limits', 'transform_AerDynBasic')
            call msg('iBin, iS: ', iBin, iS)
            call msg('d, naero', d, real(naero(iBin, iS)))
            call msg('core, sum(vols):',core(iBin, iS),sum(vols(1:nMaterials,iBin,iS)))
            call msg('min, max', fu_min_D(modesADB(iBin, iS)), fu_max_D(modesADB(iBin, iS)))
          endif
        else
          !
          ! Too low number cnc. Send to garbage but check for mass to be sure in consistency
          !
          do iMat = 1, nMaterials
            if(indADB_transp(iMat,iBin,iS) /= int_missing)then
              !
              ! If mass is sufficient, there has to be an error but since the mass and number 
              ! thresholds are largely independent, can still proceed. Just make some noise.
              !
              if(vMass(indADB_transp(iMat,iBin,iS)) > massLowTrsh(indADB_transp(iMat,iBin,iS)))then
                call set_error('Low number but enough mass cnc at the beginning', 'transform_AerDynBasic')
                call msg('iBin, nbr', iBin, vNumber(indADB_aerosol(iBin,iS)))
                call msg('Material, iS, mass:' + chADBMaterialNms(iMat), &
                       & iS, vMass(indADB_transp(iMat,iBin,iS)))
                call msg('Mass and number thresholds:', massLowTrsh(indADB_transp(iMat,iBin,iS)), &
                                                      & rulesADB%nbrLowThreshold)
                call unset_error('transform_AerDynBasic')
              endif
              garbage(indADB_transp(iMat,iBin,iS)) = garbage(indADB_transp(iMat,iBin,iS)) + &
                                                     & vMass(indADB_transp(iMat,iBin,iS))

            endif   ! material-bin-solubility combination exists
          enddo  ! materials
          vNumber(indADB_aerosol(iBin,iS)) = 0.
        endif
      enddo
    enddo
    if(error)return

    v0 = core   ! Remember the initial particle size

    ! And now calling processes ..
    !
    call msg('equilibration')
    CALL equilibration(rulesADB, naero, vols, metdat(ind_rh), metdat(ind_t), core, dwet)
    if(error)return

!    call msg('nbrs')
!    do iS = 1, nParallelBins
!      do iBin = 1,nBins 
!        if(iS == 2 .and. iBin < in2)cycle
!        call msg(fu_name(modesADB(iBin,iS)), naero(iBin,iS))
!      enddo
!    enddo

!    call msg('vols')
!    do iMat = 1, nMaterials
!      call msg(chADBMaterialNms(iMat), sum(vols(iMat,:,:)))
!    enddo

    call msg('coagulation')
    IF (rulesADB%coag_on) CALL coagulation(rulesADB, naero, vols, dwet, core, timestep_sec, &
                                         & metdat(ind_t), metdat(ind_p))
    if(error)return
  


    !- For more accurate condensation sink, dry & wet diameter must be recalculated  
    !
    IF(rulesADB%recalc)then
      do iBin = 1,nBins 
        do iS = 1, nParallelBins 
          if(naero(iBin,iS) > rulesADB%nbrLowThreshold)then 
!!!            if(all(vols(1:nMaterials,iBin,iS) / factor2volume(1:nMaterials) < massLowTrsh(indADB_transp(1:nMaterials,iBin,iS))))then
!!!              call set_error('Strange number and mass concentrations','distr_update')
!!!              return
!!!            endif
            core(iBin, iS) = sum(vols(1:nMaterials,iBin,iS)) / naero(iBin, iS)
          else
            naero(iBin,iS) = 0.
            do iMat = 1, nMaterials
!!!              if(vols(iMat,iBin,iS) / factor2volume(iMat) < massLowTrsh(indADB_transp(iMat,iBin,iS)))then
                garbage(indADB_transp(iMat,iBin,iS)) = garbage(indADB_transp(iMat,iBin,iS)) + &
                                                       & vols(iMat,iBin,iS) / factor2volume(iMat) 
                vols(iMat,iBin,iS) = 0.
                core(iBin, iS) = fu_massmean_D(modesADB(iBin, soluble))
!!!              else
!!!                call set_error('Strange number and mass concentrations', 'distr_update')
!!!                return
!!!              endif
            enddo
          endif
        enddo
      enddo
      
      call msg('equilibration')
      call equilibration(rulesADB, naero, vols, metdat(ind_rh), metdat(ind_t), core, dwet)
      if(error)return
    endif
     
    call msg('nbrs 1')
!    do iS = 1, nParallelBins
is = 1 
!      fTmp = 0.0
      do iBin = 1,nBins 
        if(iS == 2 .and. iBin < in2)cycle
        call msg(fu_name(modesADB(iBin,iS)), naero(iBin,iS))
!        fTmp = fTmp + Pi/6.0 * d**3 * rho
      enddo
!    enddo

call msg('sum ', sum(naero(:,iS)))


!    call msg('vols 1')
    do iMat = 1, nMaterials
      call msg(chADBMaterialNms(iMat), sum(vols(iMat,:,:)))
exit
    enddo

    call msg('condensation')
    IF (rulesADB%cond_on) CALL condensation(rulesADB, naero, vols, dwet, vMass, &
                                            & metdat(ind_rh), metdat(ind_t), metdat(ind_p), &
                                            & massLowTrsh, garbage, timestep_sec)
    if(error)return


    call msg('nbrs 2')
!    do iS = 1, nParallelBins
      do iBin = 1,nBins
is = 1  
        if(iS == 2 .and. iBin < in2)cycle
        call msg(fu_name(modesADB(iBin,iS)), naero(iBin,iS))
      enddo
!    enddo
call msg('sum ', sum(naero(:,iS)))
!    call msg('vols 2')
    do iMat = 1, nMaterials
      call msg(chADBMaterialNms(iMat), sum(vols(iMat,:,:)))
exit
    enddo


    !
    ! Update the aerosol distribution & return everything to its right place
    !
    call msg('distr_update')
    CALL update_distributions(rulesADB, naero, vols, v0, core, metdat(ind_t), massLowTrsh, garbage)
    if(error)return


    call msg('nbrs 3')
!    do iS = 1, nParallelBins
is = 1 
      do iBin = 1,nBins 
        if(iS == 2 .and. iBin < in2)cycle
        call msg(fu_name(modesADB(iBin,iS)), naero(iBin,iS))
      enddo
!    enddo
call msg('sum ', sum(naero(:,iS)))
!    call msg('vols 3')
    do iMat = 1, nMaterials
      call msg(chADBMaterialNms(iMat), sum(vols(iMat,:,:)))
exit
    enddo


    do iBin = 1,nBins 
      do iS = 1, nParallelBins 
        if(indADB_aerosol(iBin,iS) == int_missing)then
          do iMat = 1, nMaterials
            if(indADB_transp(iMat,iBin,iS) /= int_missing)then
              call msg('iBin, iMat',iBin,iMat)
              call set_error('Aerosol index absent but material index exists','transform_AerDynBasic')
            endif
          end do
        else
          vNumber(indADB_aerosol(iBin,iS)) = naero(iBin,iS)
          do iMat = 1, nMaterials
            if(indADB_transp(iMat,iBin,iS) /= int_missing)then
              vMass(indADB_transp(iMat,iBin,iS)) = vols(iMat,iBin,iS) / factor2volume(iMat)
            endif
          enddo
        endif
      enddo
    enddo

    ! Finally take care of sulphates created in aqueous phase and stored in short living map
!    if(iSO4_gas_sl /= int_missing)then
if(.false.)then
      call msg('distribute_SO4_aq')
      if(vMassShortliving(iSO4_gas_sl) > 0.) call distribute_SO4_aq(rulesADB, vMassShortliving, vMass, vNumber)
      if(error)return
    endif

!    vMass(:) = vMass(:) * fCellVolume
!    vMassShortliving(:) = vMassShortliving(:) * fCellVolume
!    vNumber(:) = vNumber(:) * fCellVolume

  end subroutine transform_AerDynBasic


  !***********************************************************************

  subroutine distribute_SO4_aq(rules, vMassShortliving, vMass, vNumber)
    !
    ! As ADB soesn't deal with so4 created in aqueous phase, 
    ! we have to distribute it somehow 
    !
    ! rules%mapMass_SO4aq & rules%mapNbr_SO4aq are precomputed 1 element arrays 
    !
    implicit none

    real, dimension(:), intent(inout) :: vMass, vMassShortliving, vNumber
    type(Tchem_rules_AerDynBasic), intent(in) :: rules

    integer :: i

    ! Mass
    do i = 1, rules%mapMass_SO4aq(1)%nRefSpecies
      vMass(rules%mapMass_SO4aq(1)%indSpeciesTo(i)) = vMass(rules%mapMass_SO4aq(1)%indSpeciesTo(i)) + &
                                   & vMassShortliving(iSO4_gas_sl) * rules%mapMass_SO4aq(1)%fract(i)
    enddo

    ! Number
    do i = 1, rules%mapNbr_SO4aq(1)%nRefSpecies
      vNumber(rules%mapNbr_SO4aq(1)%indSpeciesTo(i)) = vNumber(rules%mapNbr_SO4aq(1)%indSpeciesTo(i)) + &
                                   & vMassShortliving(iSO4_gas_sl) * rules%mapNbr_SO4aq(1)%fract(i)
    enddo

    vMassShortliving(iSO4_gas_sl) = 0.

  end subroutine distribute_SO4_aq


  !***********************************************************************

  SUBROUTINE equilibration(rules, n_aero, vols, rh, temp, core, dwet)

    !
    ! Calculates ambient sizes of particles by equilibrating soluble fraction of particles with water
	! 

    IMPLICIT NONE

    !-- input variables -------------
    type(Tchem_rules_AerDynBasic), intent(in) :: rules
    REAL(r8k), dimension(:,:,:), intent(inout) :: vols       ! vol cnc (nMaterials, nBins, soluble/insoluble)
    REAL(r8k), dimension(:,:), intent(inout) :: n_aero       ! nbr cnc (nBins, soluble/insoluble)
    REAL, INTENT(in) :: rh, temp                            ! relative humidity [0-1], temperature
    REAL(r8k), INTENT(in) :: core(nBins,nParallelBins)       ! particle dry volume [fxm]
    REAL(r8k), INTENT(out) :: dwet(nBins,nParallelBins)      ! particle ambient diameter [m]

    !-- local variables --------------
    integer :: iBin, iMaterial, iT, iRH
    real :: wetV, d
    type(TwetParticle) :: wetP

    ! Insoluble
    do iBin = 1, nBins
      if(core(iBin, insoluble) > 0.)then
        dwet(iBin, insoluble) = (6. * core(iBin, insoluble) / pi)**(1./3.)
      else
        dwet(iBin, insoluble) = fu_massmean_D(modesADB(iBin, insoluble))
      endif
    enddo

    ! Soluble
    ! Pick the values from precomputed table
    if (temp <= 198.) then
      iT = 1
    elseif (temp >= 348.) then
      iT = 11
	else
      iT = nint(((temp - 198.) / 15.)) + 1
    endif

    if(rh < 0.)then
      call msg_warning('negative rh')
      iRH = 1
    elseif(rh < 0.5)then
      iRH = nint(rh / 0.05) + 1
    elseif(rh < 0.7)then
      iRH = nint((rh - 0.5) / 0.04) + 11
    elseif(rh < 0.85)then
      iRH = nint((rh - 0.7) / 0.03) + 16
    elseif(rh < 0.95)then
      iRH = nint((rh - 0.85) / 0.02) + 21
    elseif(rh < 1.03)then
      iRH = nint((rh - 0.95) / 0.01) + 26
    else
      iRH = 34
    endif

    do iBin = 1, nBins
      if(n_aero(iBin, soluble) > rules%nbrLowThreshold)then
        d = (6. * core(iBin, soluble) / pi)**(1./3.)
        wetV = 0. 
        do iMaterial = 1, nMaterials
          wetV = wetV + vols(iMaterial,iBin,soluble) * rules%wetParticleGrowth(iMaterial, iBin, iRH, iT)**3 / n_aero(iBin, soluble)
        enddo
        dwet(iBin, soluble) = (6. * wetV / pi)**(1./3.)
      else
        dwet(iBin, soluble) = fu_massmean_D(modesADB(iBin, soluble))
      endif
    enddo


!    ! Soluble
!    ! Call fu_wet_particle_features to compute the wet size
!    do iBin = 1, nBins
!      core(iBin, soluble) = sum(vols(1:nMaterials,iBin,soluble)) / n_aero(iBin, soluble)
!      wetV = 0. 
!      d = (6. * core(iBin, soluble) / pi)**(1./3.)
!      do iMaterial = 1, nMaterials
!        wetP = fu_wet_particle_features(fu_get_material_ptr(chADBMaterialNms(iMaterial)), &
!            &  fTmp, temp, rh)
!        wetV = wetV + vols(iMaterial,iBin,soluble) * wetP%fGrowthFactor**3 / n_aero(iBin, soluble)
!      enddo
!      dwet(iBin, soluble) = (6. * wetV / pi)**(1./3.)
!    enddo

  END SUBROUTINE equilibration


  !***********************************************************************

  SUBROUTINE update_distributions(rules, naero, vols, v0, core, temp, massLowTrsh, garbage)

    !
    ! Uses Galperin 1d advection for absolute coordinates.
    ! If the soluble fraction of an insoluble bin exceeds the limit, 
    ! all the particles there are moved to soluble bin.
    !

    IMPLICIT NONE

    !-- Input and output variables ----------
    type(Tchem_rules_AerDynBasic) :: rules
    REAL, INTENT(IN) :: temp         ! ambient temperature [K]
    REAL(r8k), dimension(:,:,:), intent(inout) :: vols   ! vol cnc (nMaterials, nBins, soluble/insoluble)
    REAL(r8k), dimension(:,:), intent(inout) :: naero, core   ! nbr cnc (nBins, soluble/insoluble)
    real(r8k), dimension(:,:), intent(in) :: v0   ! particle size
    real, dimension(:), intent(in) :: massLowTrsh
    real, dimension(:), intent(inout) :: garbage

	!-- Local variables ----------------------
    INTEGER ::iBin, iMat, iS
    real :: volsTmp(nMaterials, nBins), v0tmp(nBins)
    real :: N(nBins), dv(nBins), binVBorder(nBins+1)
    REAL(r8k) :: epspart

    !
    ! Recompute core
    ! As advection returns garbage as only one number, send masses below low number threshold to garbage 
    ! by species and call advection with 0. threshold.
    !
    do iBin = 1,nBins 
      do iS = 1, nParallelBins 
        call msg('iBin, iS', iBin, iS)
        if(naero(iBin,iS) > rules%nbrLowThreshold)then 
          do iMat = 1, nMaterials
            if(indADB_transp(iMat,iBin,iS) /= int_missing)then
              call msg('iMat, index',iMat,indADB_transp(iMat,iBin,iS))
              call msg('vols, factor2volume', real(vols(iMat,iBin,iS)), factor2volume(iMat))
              call msg('threshold',massLowTrsh(indADB_transp(iMat,iBin,iS)))
              if(vols(iMat,iBin,iS) / factor2volume(iMat) < massLowTrsh(indADB_transp(iMat,iBin,iS)))then
                call set_error('Strange number and mass concentrations','distr_update')
                return
              endif
            endif
          end do
          core(iBin, iS) = sum(vols(1:nMaterials,iBin,iS)) / naero(iBin, iS)
        else
          naero(iBin,iS) = 0.
          do iMat = 1, nMaterials
!!!            if(vols(iMat,iBin,iS) / factor2volume(iMat) < massLowTrsh(indADB_transp(iMat,iBin,iS)))then
            if(indADB_transp(iMat,iBin,iS) /= int_missing)then
              garbage(indADB_transp(iMat,iBin,iS)) = garbage(indADB_transp(iMat,iBin,iS)) + &
                                                       & vols(iMat,iBin,iS) / factor2volume(iMat) 
            endif
            vols(iMat,iBin,iS) = 0.
            core(iBin, iS) = pi / 6. * fu_massmean_D(modesADB(iBin, iS))**3
!!!            else
!!!              call set_error('Strange number and mass concentrations', 'distr_update')
!!!              return
!!!            endif
          enddo
        endif
      enddo
    enddo

    do iBin = 1, nBins
      binVBorder(iBin) = vlolim(iBin)
    enddo
    binVBorder(nBins+1) = vhilim(nBins)


    ! Soluble bins
    N(1:nBins) = naero(:,soluble)
    volsTmp(1:nMaterials,1:nBins) = vols(:,:,soluble)

    v0tmp(1:nBins) = 0.

    ! Have to compute the initial volumes and volume change
    do iBin = 1, nBins
      if(indADB_aerosol(iBin,soluble) > 0)then
!        if(vNbr(indADB_aerosol(iBin,soluble))>rules%nbrLowThreshold)then   
!          v0(iBin) = 0.
!          do iMat = 1, nMaterials     
!            v0(iBin) =  v0(iBin) + vMass(indADB_transp(iMat,iBin,soluble)) * factor2volume(iMat)
!          enddo
!          v0(iBin) = v0(iBin) / vNbr(indADB_aerosol(iBin,soluble))
! !         v0(iBin) = sum(factor2volume(1:nMaterials) * vMass(indADB_transp(1:nMaterials,iBin,soluble))) / &
! !                 & vNbr(indADB_aerosol(iBin,soluble))
!          
          v0tmp(iBin) = v0(iBin, soluble) 
          dv(iBin) = core(iBin, soluble) - v0tmp(iBin)
!        else
!          v0tmp(iBin) = pi / 6. * fu_massmean_D(modesADB(iBin, soluble))**3
!          dv(iBin) = 0. 
!        endif
      else
        v0tmp(iBin) = pi / 6. * fu_massmean_D(modesADB(iBin, soluble))**3
        dv(iBin) = 0.
      endif
    enddo

    call advect_nbr_distr(rules, N, v0tmp, dv, binVBorder, volsTmp)
    if(error)return


    naero(1:nBins,soluble) = N(1:nBins) 
    vols(1:nMaterials,1:nBins,soluble) = volsTmp(1:nMaterials,1:nBins)


    ! Insoluble bins
    N(1:nBins) = naero(:,insoluble)
    volsTmp(1:nMaterials,1:nBins) = vols(:,:,insoluble)


    ! Have to compute the initial dry diameters and diameter growth
    do iBin = 1, nBins
      if(indADB_aerosol(iBin,insoluble) > 0)then
!        if(vNbr(indADB_aerosol(iBin,insoluble))>rules%nbrLowThreshold)then   
!          v0(iBin) = (sum(factor2volume(1:nMaterials) * vMass(indADB_transp(1:nMaterials,iBin,insoluble))) / &
!                  & vNbr(indADB_aerosol(iBin,insoluble)) * 6. / pi) ** (1./3.)
!
          v0tmp(iBin) = v0(iBin, insoluble) 
          dv(iBin) = core(iBin, insoluble) - v0tmp(iBin)

!        else
!          v0tmp(iBin) = fu_massmean_D(modesADB(iBin, insoluble))
!          dv(iBin) = 0. 
!        endif
      else
        v0tmp(iBin) = fu_massmean_D(modesADB(iBin, soluble))
        dv(iBin) = 0.
      endif
    enddo

    call advect_nbr_distr(rules, N, v0tmp, dv, binVBorder, volsTmp)
    if(error)return

    naero(1:nBins,insoluble) = N(1:nBins) 
    vols(1:nMaterials,1:nBins,insoluble) = volsTmp(1:nMaterials,1:nBins)


    !
    ! Solubility update
    !
    do iBin = in2, nBins
      if (naero(iBin, inSoluble) < rules%nbrLowThreshold)cycle
      
      !-- soluble volume fraction of particles
      epspart = (vols(iSO4,iBin,insoluble) + rules%epsoc*vols(iOC,iBin,insoluble) + &
               & vols(iSslt,iBin,insoluble) + vols(iNO3,iBin,insoluble) + vols(iNH4,iBin,insoluble)) / &
              & sum(vols(:,iBin,insoluble))

      IF (rules%epsv * vhilim(iBin) / temp**3 < epspart) THEN

        naero(iBin,soluble) = naero(iBin,soluble) + naero(iBin,insoluble)
        naero(iBin,insoluble) = 0._r8k

        vols(:,iBin,soluble) = vols(:,iBin,soluble) + vols(:,iBin,insoluble)   
        vols(:,iBin,insoluble) = 0._r8k

      END IF
    END DO

  END SUBROUTINE update_distributions 

  !*************************************************************************

  subroutine advect_nbr_distr(rules, N, v0, dv, binVLims, vols)
    
    ! Uses different distributions inside the bins, depending on the slope of the whole distribution  
    
    implicit none
     
    ! Imported
    type(Tchem_rules_AerDynBasic) :: rules
    real, dimension(:,:), intent(inout) :: vols
    real, dimension(:), intent(inout) :: N, v0, dv
 
    real :: binVLims(nBins+1)
    real :: v1, v2
    real :: nTmp(nBins)
    real :: vTmp(nMaterials,nBins)
    integer :: iBinFrom, iBinTo, iTmp
    
    real :: vTmp1, vTmp2, x, a, b
    real :: fractN(nBins), fractV(nBins)

    logical :: ifConstNbrDistr

    nTmp(1:nBins) = 0.
    vTmp(1:nMaterials, 1:nBins) = 0.

    do iBinFrom = 1, nBins

      ! Some checking
      if(n(iBinFrom) < rules%nbrLowThreshold)cycle ! Don't know the transport indices, so someone else has to garbage

      if(v0(iBinFrom) > binVLims(iBinFrom+1) .or. v0(iBinFrom) < binVLims(iBinFrom))then
        call set_error('Strange particles before aerosol dynamics','advect_nbr_distr')
        call msg('Bin, particle size: ', iBinFrom, (v0(1)* 6 / pi)**(1./3.))
        return
      endif

      if(((v0(iBinFrom) + dv(iBinFrom))* 6 / pi)**(1./3.) < 1e-9 .or. &
       & ((v0(iBinFrom) + dv(iBinFrom))* 6 / pi)**(1./3.) > 1e-4)then
        call msg_warning('Strange particles after aerosol dynamics','advect_nbr_distr')
        call msg('Bin, particle size: ', iBinFrom, ((v0(1) + dv(1))* 6 / pi)**(1./3.))
        return
      endif

      ! Ends of the distribution
      if((iBinFrom == 1 .and. dv(iBinFrom) < 0.) .or. (iBinFrom == nBins .and. dv(iBinFrom) > 0.))then
        nTmp(iBinFrom) = nTmp(iBinFrom) + N(iBinFrom)
        vTmp(1:nMaterials, iBinFrom) = vTmp(1:nMaterials, iBinFrom) + vols(1:nMaterials, (iBinFrom))
        cycle
      endif

      ! Advection
      fractN(1:nBins) = 0.      
      fractV(1:nBins) = 0.
    
      ! First have to set V1 and V2
      if((n(iBinFrom) / (fu_max_d(modesADB(iBinFrom,1)) - fu_min_d(modesADB(iBinFrom,1)))) < &
       & (n(iBinFrom+1) / (fu_max_d(modesADB(iBinFrom+1,1)) - fu_min_d(modesADB(iBinFrom+1,1)))))then
        ifConstNbrDistr = .true.  
      else
        ifConstNbrDistr = .false.
        v2 = binVLims(iBinFrom+1)

        if(3. * (4. * v0(iBinFrom) - v2) >= 0.)then
          x = (- sqrt(v2) + sqrt(3. * (4. * v0(iBinFrom) - v2))) / 2. 
          v1 = x**2
        else
          x = -1.
          v1 = -1.
        endif
        if(x < 0. .or. v1 > binVLims(iBinFrom+1) .or. v1 < binVLims(iBinFrom))then
          v1 = binVLims(iBinFrom)
          if(3. * (4. * v0(iBinFrom) - v1) >= 0.)then
            x = (- sqrt(v1) + sqrt(3. * (4. * v0(iBinFrom) - v1))) / 2. 
            v2 = x**2
          else
            x = -1.
            v1 = -1.
          endif
          if(x < 0. .or. v2 > binVLims(iBinFrom+1) .or. v2 < binVLims(iBinFrom))then
            call msg('Failed to set square root nbr distr')
            ifConstNbrDistr = .true.
          endif
        endif
      endif

!ifConstNbrDistr = .true.

      if(ifConstNbrDistr)then

        if((binVLims(iBinFrom+1) - v0(iBinFrom)) < (v0(iBinFrom) - binVLims(iBinFrom)))then 
          v2 = binVLims(iBinFrom+1) 
          v1 = v2 - 2. * (v2 - v0(iBinFrom))
        else
          v1 = binVLims(iBinFrom)
          v2 = v1 + 2. * (v0(iBinFrom) - v1)
        endif
      endif
!call msg('')
!call msg('Bin from: ', iBinFrom)
!if(iBinFrom == 100)then
!  call msg('Reached 100-th bin')
!endif
!if(.not. ifConstNbrDistr) call msg('SQRT DISTR')
!call msg('Bin from limits: ', binVLims(iBinFrom), binVLims(iBinFrom+1))
!call msg('Advected from limits: ', v1, v2)
!call msg('v0, dv: ', v0(iBinFrom), dv(iBinFrom))

      ! Uniform shift, assuming that all the particles in the bin gain the same volume
      v1 = v1 + dv(iBinFrom)
      v2 = v2 + dv(iBinFrom)

!call msg('Advected to limits: ', v1, v2)

      ! And now we have to figure out, where it ended up and divide the particles and volumes to appropriate bins 
      do iBinTo = 1, nBins
        if(v1 < binVLims(iBinTo+1))exit
      enddo
            
      if(v2 <= binVLims(iBinTo+1))then ! All to one bin
!call msg('All to bin: ', iBinTo)
        nTmp(iBinTo) = nTmp(iBinTo) + N((iBinFrom))
        vTmp(1:nMaterials, iBinTo) = vTmp(1:nMaterials, iBinTo) + vols(1:nMaterials, (iBinFrom))

      else

        iTmp = iBinTo
        do while (iTmp <= nBins)

!call msg('Bin to: ', iTmp)
!call msg('Bin to limits: ', binVLims(iTmp), binVLims(iTmp+1))

          vTmp1 = max(v1, binVLims(iTmp))
          vTmp2 = min(v2, binVLims(iTmp+1))

!call msg('Integration limits: ', vTmp1, vTmp2)

          if(ifConstNbrDistr)then  

            fractN(iTmp) = (vTmp2 - vTmp1) / (v2 - v1)
            fractV(iTmp) = (vTmp2**2 - vTmp1**2) / (v2**2 - v1**2)
!            fractV(iTmp) = (vTmp1 + vTmp2) / 2. * fractN(iTmp) / (v0(iBinFrom) + dv(iBinFrom)) 
          else

            x = 1. / (sqrt(v2 - dv(iBinFrom)) - sqrt(v1 - dv(iBinFrom)))
            b = sqrt(vTmp2 - dv(iBinFrom))
            a = sqrt(vTmp1 - dv(iBinFrom))
            fractN(iTmp) = (b - a) * x
            fractV(iTmp) = ((b**3. - a**3.) / 3. +  dv(iBinFrom) * (b - a)) * N(iBinFrom) * x &
                         & / sum(vols(1:nMaterials, iBinFrom)) 
          endif
          
          nTmp(iTmp) = nTmp(iTmp) + N(iBinFrom) * fractN(iTmp)
          vTmp(1:nMaterials, iTmp) = vTmp(1:nMaterials, iTmp) + vols(1:nMaterials, iBinFrom) * fractV(iTmp) 

          if((sum(vTmp(1:nMaterials, iTmp)) / nTmp(iTmp)) < binVLims(iTmp) .or. &
           & (sum(vTmp(1:nMaterials, iTmp)) / nTmp(iTmp)) > binVLims(iTmp+1))then
            call set_error('Strange advected volume','advect_nbr_distr')
            call msg('v', sum(vTmp(1:nMaterials, iTmp)) / nTmp(iTmp))
            call msg('limits', binVLims(iTmp), binVLims(iTmp+1))
            call msg('Fractions for n and v: ', fractN(iTmp), fractV(iTmp))
            call unset_error('advect_nbr_distr')
          endif


          if(v2 < binVLims(iTmp+1))exit
          iTmp = iTmp + 1
        enddo
        
        ! Check here, that fractions sum up to one
        if(abs(sum(fractN)-1.) > 0.01)then
          call set_error('Numbers not conserved','advect_nbr_distr')
          call msg('Fractions sum', sum(fractN))
        endif
        if(abs(sum(fractV)-1.) > 0.01)then
          call set_error('Masses not conserved','advect_nbr_distr')
          call msg('Fractions sum', sum(fractV))
        endif

      endif
    enddo ! Bins from

    ! And copy the temporary arrays back to the main ones
    N(1:nBins) = nTmp(1:nBins)
    vols(1:nMaterials,1:nBins) = vTmp(1:nMaterials,1:nBins)

  end subroutine advect_nbr_distr

  !***********************************************************************

  SUBROUTINE coagulation(rules, naero, vols, dwet, core, tstep, t, p)

  !
  ! Purpose:
  ! --------
  ! Calculates particle loss and change in size distribution
  !  due to (Brownian) coagulation
  !
  ! Method:
  ! -------  
  !  implicit, non-iterative method:
  !  Volume concentrations of the smaller colliding particles
  !  added to the bin of the larger colliding particles.
  !
  !!  Exact coagulation coefficients for each pressure level
  !!  are calculated in subroutine SET_COAGC (in init_ADB) 
  !!  which is called once at the beginning of the simulation 
  !!  In subroutine COAGULATION, these exact 
  !!  coefficients are scaled according to current particle wet size
  !!  (linear scaling).
  !
  !---------------------------------------------------------------------
  !

    IMPLICIT NONE

    type(Tchem_rules_AerDynBasic), intent(in) :: rules
   
    REAL(r8k), dimension(:,:,:), intent(inout) :: vols   ! vol cnc (nMaterials, nBins, soluble/insoluble)
    REAL(r8k), dimension(:,:), intent(inout) :: naero    ! nbr cnc (nBins, soluble/insoluble)

    real, intent(in) :: t, p                          ! index for coagtable 
    REAL(r8k), dimension(:,:), INTENT(IN) :: dwet, core  ! wet diameter & dry volume (nBins, soluble/insoluble)
    REAL, INTENT(IN) :: tstep                           ! timestep [s]
!    REAL(r8k), dimension(:,:,:), INTENT(IN) :: coagtable ! precomputed coefficients for bin mid try diameter (lev,s,l) 
!                                                        ! no soluble/insoluble dimension necessary

    !-- Local variables ------------------------
    real(r8k), dimension(nBins,nBins,nParallelBins,nParallelBins) :: cct ! updated coefs (smaller_particle, larger_particle, soluble/insoluble, soluble/insoluble) 
    integer :: s, l, k                 ! indices 
    REAL(r8k), dimension(nBins,nParallelBins) :: nTmp    ! temporary array for numbers     
    REAL(r8k), dimension(nMaterials,nBins,nParallelBins) :: vTmp    ! temporary array for volume changes 
    real(r8k) :: ftmp, fTmp2
    REAL(r8k), dimension(nMaterials) :: dif    ! temporary array for volume changes 
   integer :: iMat


    ! Need the number change later for volumes, so have to store this temporarily
    nTmp = naero
    vTmp = 0. 
    cct(:,:,:,:) = 0.

    !
    ! If horribly slow, have to go back to precomputation, but for now:
    DO l = 1,nBins  ! larger colliding particle 
      DO s = 1,nBins  ! smaller colliding particle
        fTmp = (pi/6)*dwet(s,soluble)**3 * 1500.
        fTmp2 = (pi/6)*dwet(l,soluble)**3 * 1500.
        cct(s,l,soluble,soluble) = tstep * coagc(dwet(s,soluble),dwet(l,soluble), fTmp, fTmp2, t, p)        
        if(l > fn1)then
          fTmp2 = (pi/6)*dwet(l,insoluble)**3 * 1500.
          cct(s,l,soluble, insoluble) = tstep * coagc(dwet(s,soluble),dwet(l,insoluble), fTmp, fTmp2, t, p)
        endif
        if(s > fn1 .and. l > fn1)then
          fTmp = (pi/6)*dwet(s,insoluble)**3 * 1500.
          cct(s,l,insoluble, insoluble) = tstep * coagc(dwet(s,insoluble),dwet(l,insoluble), fTmp, fTmp2, t, p)
        endif
        if(s > fn1)then
          fTmp2 = (pi/6)*dwet(l,soluble)**3 * 1500.
          cct(s,l,insoluble,soluble) = tstep * coagc(dwet(s,insoluble),dwet(l,soluble), fTmp, fTmp2, t, p)
        endif

!        if(l == s)then
!          cct(s,l,soluble,soluble) = 0.5 * cct(s,l,soluble,soluble)
!          cct(s,l,insoluble, insoluble) = 0.5 * cct(s,l,insoluble, insoluble)
!        endif

      END DO
    END DO

    !-- 1) Updating coagulation coefficients (the whole horrible 4d matrix for 2d bins)
    !      timestep and factor 0.5 for self coagulation also go to cc

!    DO l = 1,nBins  ! larger colliding particle 
!      DO s = 1,nBins  ! smaller colliding particle
!        cct(s,l,soluble,soluble) = tstep * rules%coagtable(lev,s,l) * dpmid(s) * dwet(l,soluble) / (dpmid(l) * dwet(s,soluble)) 
!        cct(s,l,soluble, insoluble) = tstep * rules%coagtable(lev,s,l) * dpmid(s) * dwet(l,insoluble) / (dpmid(l) * dwet(s,soluble)) 
!        cct(s,l,insoluble,soluble) = tstep * rules%coagtable(lev,s,l) * dpmid(s) * dwet(l,soluble) / (dpmid(l) * dwet(s,insoluble)) 
!        cct(s,l,insoluble, insoluble) = tstep * rules%coagtable(lev,s,l) * dpmid(s) * dwet(l,insoluble) / (dpmid(l) * dwet(s,insoluble)) 
!        if(l == s)then
!           cct(s,l,soluble,soluble) = 0.5 * cct(s,l,soluble,soluble)
!           cct(s,l,insoluble, insoluble) = 0.5 * cct(s,l,insoluble, insoluble)
!        endif
!      END DO
!    END DO


    ! The idea is, that the coagulated creature stays in the bin of the bigger
    ! parent. In case of same size soluble/insoluble parents it goes to soluble bin. 
    ! So .. we have to compute every bins losses to bigger soluble and insoluble bins 
    ! (for insoluble guys also to soluble same size) + selfcoagulation
    ! 
    do k = nBins, 1, -1
      ! Soluble
      if(naero(k,soluble) > 0.)then
        ftmp = 1. + sum(cct(k,k+1:nBins,soluble,soluble)*naero(k+1:nBins,soluble))            ! soluble
        ftmp = ftmp + sum(cct(k,k+1:nBins,soluble,insoluble)*naero(k+1:nBins,insoluble)) ! insoluble
      
        naero(k,soluble) = ( - ftmp + sqrt(ftmp**2. + 4. * cct(k,k,soluble,soluble) * naero(k,soluble))) / &
                         & (2. * cct(k,k,soluble,soluble))
      endif
      ! Insoluble
      if(k > fn1)then
        if(naero(k,insoluble) > 0.)then
          ftmp = 1. + sum(cct(k,k:nBins,insoluble,soluble)*naero(k:nBins,soluble))            ! soluble
          ftmp = ftmp + sum(cct(k,k+1:nBins,insoluble,insoluble)*naero(k+1:nBins,insoluble)) ! insoluble       
    
          naero(k,insoluble) = ( - ftmp + sqrt(ftmp**2 + 4. * cct(k,k,insoluble,insoluble) * naero(k,insoluble))) / &
                              & (2. * cct(k,k,insoluble,insoluble))
        endif
      endif
    enddo

    ! Distribute the volumes
    ! First bin 
    ! normalization factor  
    ftmp = sum(cct(1,1:nBins,soluble,soluble) * 0.5 * (naero(1:nBins,soluble) + nTmp(1:nBins,soluble)))               ! soluble
    ftmp = ftmp + sum(cct(1,in2:nBins,soluble,insoluble) * 0.5 * (naero(in2:nBins,insoluble) + nTmp(in2:nBins,insoluble))) ! insoluble
 
    ! fraction of original particles removed
    ftmp2 = (nTmp(1,soluble) - naero(1,soluble)) / nTmp(1,soluble)
    if(ftmp2 > 0.)then
      do l = 1, nBins  ! some "removed" to the same bin by selfcoagulation
        vTmp(:,l,soluble) = vols(:,1,soluble) * ftmp2 * cct(1,l,soluble,soluble) * 0.5 * (naero(l,soluble) + nTmp(l,soluble)) / fTmp

        vols(:,1,soluble) = vols(:,1,soluble) - vTmp(:,l,soluble)
        if(l > fn1)then
          vTmp(:,l,insoluble) = vols(:,1,soluble) * fTmp2 * cct(k,l,soluble,insoluble) * 0.5 * (naero(l,insoluble) + nTmp(l,insoluble)) / fTmp
          vols(:,1,soluble) = vols(:,1,soluble) - vTmp(:,l,insoluble)

        endif
      enddo
    endif

    ! Now cycle the rest
    do k = 2, nBins
      ! Soluble part
      ! Normalization factor 
      ftmp = sum(cct(k,k:nBins,soluble,soluble) * 0.5 * (naero(k:nBins,soluble) + nTmp(k:nBins,soluble)))            ! soluble
      ftmp = ftmp + sum(cct(k,k+1:nBins,soluble,insoluble) * 0.5 * (naero(k+1:nBins,insoluble) + nTmp(k+1:nBins,insoluble))) ! insoluble

      ! fraction of particles removed (some "removed" to the same bin by selfcoagulation)
      ftmp2 = (nTmp(k,soluble) - naero(k,soluble)) / nTmp(k,soluble)
      if(ftmp2 > 0.)then
        do l = k, nBins 

          dif(:) = (vols(:,k,soluble) + 0.5 * vTmp(:,k,soluble)) * ftmp2 * &
               & cct(k,l,soluble,soluble) *  0.5 * (naero(l,soluble) + nTmp(l,soluble)) / fTmp 
                     
          vTmp(:,l,soluble) = vTmp(:,l,soluble) + dif(:)
                      
          vols(:,k,soluble) = vols(:,k,soluble) - dif(:)
     
          if(l > fn1 .and. l > k)then
            dif(:) = fTmp2 / fTmp *  cct(k,l,soluble,insoluble) * &
                     & 0.5 * (naero(l,insoluble) + nTmp(l,insoluble)) * &
                     & (vols(:,k,soluble) + 0.5 * vTmp(:,k,soluble))

            vTmp(:,l,insoluble) = vTmp(:,l,insoluble) + dif(:) 
                      
            vols(:,k,soluble) = vols(:,k,soluble) - dif(:)
          endif     
        enddo
      endif
      ! Insoluble part
      if(k > fn1)then
        ! Normalization factor 
        ftmp = sum(cct(k,k:nBins,insoluble,soluble) * 0.5 * (naero(k:nBins,soluble) + nTmp(k:nBins,soluble)))           ! soluble
        ftmp = ftmp + sum(cct(k,k:nBins,insoluble,insoluble) * 0.5 * (naero(k:nBins,insoluble) + nTmp(k:nBins,insoluble))) ! insoluble
          
        ! fraction of particles removed
        fTmp2 = (nTmp(k,insoluble)-naero(k,insoluble)) / nTmp(k,insoluble)
        if(fTmp2 > 0.)then
          do l = k, nBins  
            dif(:) =  fTmp2 / fTmp *  cct(k,l,insoluble,soluble) * &
                       & 0.5 * (naero(l,soluble) + nTmp(l,soluble)) * &
                       & (vols(:,k,insoluble) + 0.5 * vTmp(:,k,insoluble))
           
            vTmp(:,l,soluble) = vTmp(:,l,soluble) +  dif(:)
       
            vols(:,k,insoluble) = vols(:,k,insoluble) - dif(:)

     
            dif(:) = fTmp2 / fTmp *  cct(k,l,insoluble,insoluble) * &
                       & 0.5 * (naero(l,insoluble) + nTmp(l,insoluble)) * &
                       & (vols(:,k,insoluble) + 0.5 * vTmp(:,k,insoluble))
            
            vTmp(:,l,insoluble) = vTmp(:,l,insoluble) + dif(:)
       
            vols(:,k,insoluble) = vols(:,k,insoluble) - dif(:)
           
          enddo
        endif
      endif
    enddo
    vols = vols + vTmp

    contains

!====================================================================================================

      FUNCTION coagc(diam1,diam2,mass1,mass2,temp,pres)

      !
      ! Purpose:
      ! --------
      ! Calculates the coagulation coefficient for two
      ! colliding particles. 
      !
      !
      ! Method: 
      ! -------  
      ! Only Brownian coagulation taken into account.
      ! Transition regime correction is done with Fuchs  
      ! flux matching.
      !
      !
      ! Interface:
      ! ----------
      ! Called from subroutine SET_COAG
      ! (only at the beginning of simulation)
      !
      !
      ! Coded by:
      ! ---------
      ! Hannele Korhonen (FMI) 2005 
      !

        IMPLICIT NONE

        !-- Input variables ----------
        REAL(r8k), INTENT(IN) :: diam1, diam2   ! diameters of colliding particles [m]
        REAL(r8k), INTENT(IN) :: mass1, mass2   ! masses -"- [kg]
        
        real, intent(in) ::     temp,   &   ! ambient temperature [K]
                                pres        ! ambient pressure [fxm]

        !-- Output variables ---------
        REAL(r8k) :: coagc       ! coagulation coefficient of particles [m3/s]

        !-- Local variables ----------  
        REAL(r8k) :: visc,   &   ! viscosity of air [kg/(m s)]
                    mfp,    &   ! mean free path of air molecules [m]
                    mdiam,  &   ! mean diameter of colliding particles [m]
                    fmdist      ! distance of flux matching [m]

        REAL(r8k), DIMENSION (2) :: diam,   &   ! diameters of particles [m]
                                   mpart,  &   ! masses of particles [kg]
                                   knud,   &   ! particle knudsen number [1]
                                   beta,   &   ! Cunningham correction factor [1]
                                   dfpart, &   ! particle diffusion coefficient [m2/s]
                                   mtvel,  &   ! particle mean thermal velocity [m/s]
                                   omega,  &   !
                                   tva,    &   ! temporary variable [m]
                                   flux        ! flux in continuum and free molec. regime [m/s]


        if((diam1 <= 0._r8k) .or. (diam2 <= 0._r8k) .or. (mass1 <= 0._r8k) .or. (mass2 <= 0._r8k))then
    !      call msg_warning('Coagulating zeroes', 'coagc')
          coagc = 0.
          return
        endif

        !-- 0) Initializing particle and ambient air variables --------------------
        diam = (/ diam1, diam2 /)       ! particle diameters [m]
        mpart = (/ mass1, mass2 /)       ! particle masses [kg]
        visc = (7.44523e-3_r8k*temp**1.5_r8k)/(5093._r8k*(temp+110.4_r8k)) ! viscosity of air [kg/(m s)]
        mfp = (1.656e-10_r8k*temp+1.828e-8_r8k)*std_pressure_sl/pres ! mean free path of air [m]

        !-- 2) Slip correction factor for small particles -------------------------
        knud = 2._r8k*mfp/diam                                    ! Knudsen number
        beta = 1._r8k+knud*(1.142_r8k+0.558_r8k*exp(-0.999_r8k/knud))! Cunningham correction factor
        ! (Allen and Raabe, Aerosol Sci. Tech. 4, 269)

        !-- 3) Particle properties ------------------------------------------------
        dfpart = beta*boltzmann_const*temp/(3._r8k*pi*visc*diam)  ! diffusion coefficient [m2/s]
        mtvel = sqrt((8._r8k*boltzmann_const*temp)/(pi*mpart))    ! mean thermal velocity [m/s]
        omega = 8._r8k*dfpart/(pi*mtvel)
        mdiam = 0.5_r8k*(diam(1)+diam(2))               ! mean diameter [m]

        !-- 4) Calculation of fluxes and flux matching ----------------------------
        flux(1) = 4._r8k*pi*mdiam*(dfpart(1)+dfpart(2)  )    ! flux in continuum regime [m3/s]
        flux(2) = pi*sqrt(mtvel(1)**2+mtvel(2)**2)*mdiam**2 !  -"- in free molec. regime [m3/s]
        ! temporary variable [m]
        tva(1) = ((mdiam+omega(1))**3 - (mdiam**2+omega(1)**2)* &
                  sqrt((mdiam**2+omega(1)**2)))/(3._r8k*mdiam*omega(1)) - mdiam
        ! temporary variable [m]
        tva(2) = ((mdiam+omega(2))**3 - (mdiam**2+omega(2)**2)* &
                 sqrt((mdiam**2+omega(2)**2)))/(3._r8k*mdiam*omega(2)) - mdiam
        ! flux matching distance [m]
        fmdist = sqrt(tva(1)**2+tva(2)**2)             

        !-- 5) Coagulation coefficient [m3/s] -------------------------------------
        coagc = flux(1) / (mdiam/(mdiam+fmdist) + flux(1)/flux(2)) 

      END FUNCTION coagc


  END SUBROUTINE coagulation


  !***********************************************************************

  SUBROUTINE condensation(rules, naero, vols, dwet, vMass, rh, temp, pres, massLowTrsh, garbage, tstep)

  !
  !  Calculates the increase in particle volume and decrease in gas phase concentrations 
  !  due to condensation of sulphuric acid, ammonia, nitric acid and non-volatile
  !  and semivolatile organic compounds.
  !  Sulphuric acid and non-volatile organic compounds condense onto all sized particles. 
  !  Semivolatile organic compounds condense only onto particles in regimes 2 & 3
  !  Condensation of ammonia and nitrates is computed using equilibrium gas concentrations 
  !  for particles from ISORROPIA and taking into account Kelvin effect.
  !  Sulphuric acid condensation is computed together with nucleation  
  !
  !  New gas and aerosol phase concentrations calculated according to Jacobson (1997): 
  !  Numerical techniques to solve condensational and dissolutional growth equations 
  !  when growth is coupled to reversible reactions, Aerosol Sci. Tech., 27, pp 491-498.
  !
  ! Following parameterization has been used:
  ! ------------------------------------------
  !
  ! Molecular diffusion coefficient of condensing vapour [m2/s] 
  ! (Reid et al. (1987): Properties of gases and liquids, McGraw-Hill, New York.)
  !
  ! D = {1.d-7*sqrt(1/M_air + 1/M_gas)*T^1.75} / {p_atm/p_stand * (d_air^(1/3) + d_gas^(1/3))^2 }
  !   
  ! M_air = 28.965 : molar mass of air [g/mol]
  ! d_air = 19.70  : diffusion volume of air
  ! M_h2so4 = 98.08  : molar mass of h2so4 [g/mol]
  ! d_h2so4 = 51.96  : diffusion volume of h2so4
  !
    IMPLICIT NONE

    type(Tchem_rules_AerDynBasic), intent(in) :: rules

    real, dimension(:), intent(inout) ::  vMass, garbage
    real, dimension(:), intent(in) :: massLowTrsh


    !-- Input and output variables ----------

    REAL(r8k), dimension(:,:), INTENT(IN) :: dwet              ! wet diameter of particles in each bin [m]
    REAL, INTENT(IN) ::     rh,                             & ! relative humidity [0-1]
                            temp,                           & ! ambient temperature [K]
                            pres,                           & ! ambient pressure [Pa]
                            tstep                             ! timestep [s]

    REAL(r8k), dimension(:,:,:), intent(inout) :: vols   ! vol cnc (nMaterials, nBins, soluble/insoluble)
    REAL(r8k), dimension(:,:), intent(inout) :: naero    ! nbr cnc (nBins, soluble/insoluble)

    !-- Local variables ----------------------
    INTEGER :: class, iS     ! loop indices

    REAL(r8k) :: visc,                               & ! viscosity of air [kg/(m s)]
                dfvap,                              & ! air diffusion coefficient [m2/s]
                mfp,                                & ! mean free path of condensing vapour [m]              
                dfpart(in1:fn1),                    & ! particle diffusion coefficient
                knud(nBins, nParallelBins),         & ! particle Knudsen number
                beta(nBins, nParallelBins),         & ! transitional correction factor
                colrate(nBins, nParallelBins),      & ! collision rate of molecules to particles [1/s]
                cs_tot,                             & ! total condensation sink [1/s] 
                cvap_new,                           & ! vapour concentration after time step [#/m3]
                dvap,                               & ! change in vapour concentration [#/m3]
                dvol(nBins, nParallelBins),         & ! change of volume in each bin                 
                ceqno(nBins, nParallelBins),        & ! equilibrium concentration of HNO3
                ceqnh(nBins, nParallelBins),        & ! equilibrium concentration of NH3
                kelvin(nBins, nParallelBins)          ! kelvin effect for all sizes
   
    real :: aNH3, bNH3, aHNO3, bHNO3, fTmp

    ! Nucleation
    REAL(r8k) :: jnuc,       & ! nucleation rate at ~ 1 nm [#/(
                nsa,        & ! number of H2SO4 molecules in critical cluster [1]
                dcrit,      & ! diameter of critical cluster [m]
                grclust,    & ! growth rate of formed clusters [nm/h]
                csink,      & ! condensational sink [#/m2]
                dmean,      & ! mean diameter of existing particles [nm]
                gamma,      & ! proportionality factor [(nm2 m2)/h]
                eta,        & ! constant; proportional to ratio of CS/GR [m]
                j3,         & ! number conc. of formed 3 nm particles [#/m3]
                j3n3,       & ! nucleation rate at ~ 1 nm [#/(m3 s)]
                n_vs_c        ! ratio of nucleation of all mass transfer in the smallest bin

    real(r8k), parameter :: n3 = 158.79_r8k        ! number of H2SO4 molecules in 3 nm cluster assuming d_sa = 5.54 �     
    REAL(r8k), PARAMETER :: d_sa   = 5.539376964394570e-10   ! diameter of condensing molecule [m]


    ! FROM/TO ISOROPIA
    REAL(r8k) :: WI_S(5)                ! Total species concentrations in moles/m**3 air
    REAL(r8k) :: CNTRL_S(2)             ! nug for different types of problems solved
                                       ! different state of aerosols (deliquescent or metastable)
    REAL(r8k) :: WT_S(1)                ! nucleation rate at ~ 1 nm [#/(m3 s)]
    REAL(r8k) :: AERLIQ_S(12)           ! Aqueous-phase concentration array in moles/m**3air
    REAL(r8k) :: AERSLD_S(9)            ! Solid-phase concentration array in moles/m**3 air  
    REAL(r8k) :: OTHER_S(6)             ! Solution information array
    REAL(r8k) :: GAS_S(3)               ! Gas-phase concentration array in moles/m**3 air
    CHARACTER*15 SCASI_S               ! Returns the subcase which the input corresponds to

    !
    ! Properties of air and condensing gases --------------------
    !
!    visc  = (7.44523e-3_r8k * temp**1.5_r8k ) / (5093._r8k * (temp + 110.4_r8k)) ! viscosity of air [kg/(m s)] 
    visc = fu_dynamic_viscosity(temp)
    dfvap = 5.1111e-10_r8k * temp**1.75_r8k * std_pressure_sl / pres          ! diffusion coefficient [m2/s]
    mfp   = 3._r8k * dfvap * sqrt(pi * molarMass(iSO4) / (8._r8k * gas_constant_uni * temp))    ! mean free path [m]

    !
    ! Transition regime correction factor for particles
    ! Fuchs and Sutugin (1971), In: Hidy et al. (ed.) Topics in current aerosol research, Pergamon.  
    ! Size of condensing molecule considered only for nucleation mode (3 - 20 nm) 
    !
    ! particle Knudsen number
    !
    knud(in1:fn1, :) = 2._r8k * mfp / (dwet(in1:fn1,:) + d_sa) 
    knud(in2:nBins,:) = 2._r8k * mfp / dwet(in2:nBins,:)
    !
    ! transitional correction factor
    !
    beta = (knud + 1._r8k) / (0.377_r8k * knud + 1._r8k + 4._r8k / (3._r8k * rules%massacc) * (knud + knud**2))  

    !
    ! Collision rate of molecules to particles
    !
    !  Particle diffusion coefficient [m2/s] considered only for nucleation mode (3 - 20 nm) 
    !
    dfpart = boltzmann_const * temp * beta(in1:fn1,soluble) / (3._r8k * pi * visc * dwet(in1:fn1,soluble))  

    !
    ! collision rate [1/s]
    !
    colrate(in1:fn1,insoluble) = 0._r8k     !empty
    colrate(in1:fn1,soluble) =  2._r8k * pi * (dfvap + dfpart) * (dwet(in1:fn1,soluble) + d_sa) * &
                              & beta(in1:fn1,soluble) * naero(in1:fn1,soluble)     
    colrate(in2:nBins,:) = 2._r8k * pi * dfvap * dwet(in2:nBins,:) * beta(in2:nBins,:) * naero(in2:nBins,:)  

!colrate = 2.4*colrate

    cs_tot = sum(colrate)                              ! total sink 


vMass(iH2SO4_gas) = 4.078e-8
call msg('H2SO4 before', vMass(iH2SO4_gas)) 

    !-- 5) Changes in gas-phase concentrations and particle volume -----
    !--- 5.1) Sulphuric acid -------------------------
    if(iH2SO4_gas /= int_missing)then
      if(vMass(iH2SO4_gas) > massLowTrsh(iH2SO4_gas))then

        ! Nucleation 
        j3n3 = 0._r8k
        n_vs_c = 0._r8k
        csink = cs_tot
        IF (rules%nucl_on)then         
          SELECT CASE (rules%nucl_scheme)
          CASE(1) ! binary H2SO4-H2O nucleation
            CALL binnucl(vMass(iH2SO4_gas), temp, rh, jnuc, nsa, dcrit)   
          CASE(2) ! ternary H2SO4-H2O-NH3 nucleation
            call set_error('Ternary nucleation crashes the compiler','')
!            CALL ternucl(vMass(iH2SO4_gas), vMass(iNH3), temp, rh, jnuc, nsa, dcrit)             
          CASE(3) ! kinetically limited nucleation of (NH4)HSO4 clusters
            CALL kinnucl(vMass(iH2SO4_gas), temp, jnuc, dcrit)
            nsa = 2._r8k
          CASE(4) ! activation type nucleation
            call actnucl(vMass(iH2SO4_gas), jnuc, dcrit)
            nsa = 2._r8k
          END SELECT
  
          !-- 2) Change of particle and gas concentrations ------------------------------------
          IF (jnuc > 1.e3_r8k) then  !-- very small nucleation rates neglected

            !-- 2.1) Check that there is enough H2SO4 to produce the nucleation ---------
           
            dvap = tstep * jnuc * nsa / avogadro 
            if(dvap >= vMass(iH2SO4_gas))then
              call set_error('Strange nucleation rate', 'condensation')
              call unset_error('condensation')
              dvap = vMass(iH2SO4_gas)
            endif

            cvap_new = vMass(iH2SO4_gas) - dvap     ! H2SO4 concentration after nucleation [#/m3]

            !-- 2.2) Formation rate of 3 nm particles (Kerminen & Kulmala, 2002) -----------------
            !--- 2.2.1) Growth rate of formed clusters by H2SO4 (Kerminen & Kulmala, 2002)
            !  GR = 3d-15 / cluster_density * molecular_speed * molarmass * concentration 
            grclust = 2.3623e-15_r8k * sqrt(temp) * cvap_new * avogadro  ! [nm/h]

            !--- 2.2.2) Condensational sink of pre-existing particle population    
            !--- condensational sink [#/m2]
            csink = sum(dwet(:,:) * beta * naero(:,:))   
    
            !--- 2.2.3) Parameterized formation rate of detectable 3 nm particles
            !--- Constants needed for the parameterization
            ! The following values have been used to simplify the parameterization presented in
            ! Kerminen & Kulmala (2002):
            ! dapp = 3 nm
            ! dens_nuc = 1830 kg/m3 

            IF(csink < 1.e-30_r8k) THEN 
              eta = 0._r8k
            ELSE
              dmean = sum(naero(:,:) * dwet(:,:)) * 1.e9_r8k / sum(naero(:,:)) ! mean particle diameter [nm]
              gamma = 0.18842_r8k * dcrit**0.2_r8k * (dmean / 150._r8k)**0.048_r8k * (temp / 293._r8k)**(-0.75_r8k) ! [nm2*m2/h] 
              eta = gamma * csink / grclust      ! [nm]
            END IF

            !-- Number conc. of clusters surviving to 3 nm in a time step
            j3 = jnuc * exp(eta / 3._r8k - eta / (dcrit * 1.e9_r8k)) * tstep ! [#/m3]
  
            !-- If J3 very small (< 1 #/cm3), neglect particle formation
            !  In real atmosphere this would mean that clusters form but coagulate to pre-existing
            !  particles who gain sulphate.Since CoagS ~ CS, we do *not* update H2SO4 concentration 
            !  here but let condensation take care of it.
            IF (j3 > 1.e6_r8k)then
              ! change in H2SO4 concentration [#/m3]
              j3n3 = j3*n3 / tstep      
              !-- Ratio of mass transfer between nucleation and condensation
              n_vs_c = j3n3 / (j3n3 + vMass(iH2SO4_gas) * avogadro * colrate(in1,soluble))
              !   collision rate in the smallest bin, including nucleation and condensation
              !   see Mark Z. Jacobson, Fundamentals of Atmospheric Modeling, Second Edition (2005)
              !   equation (16.73)  
              colrate(in1,soluble) = colrate(in1,soluble) + j3n3 / (vMass(iH2SO4_gas)* avogadro)
              ! total sink for sulfate
              csink = cs_tot + j3n3 / (vMass(iH2SO4_gas)*avogadro)
            endif
          endif
        endif
        IF(csink > 1.e-30_r8k ) THEN
          cvap_new = vMass(iH2SO4_gas) / (1._r8k + tstep * csink)    ! new gas phase concentration [#/m3]
          dvap = vMass(iH2SO4_gas) - cvap_new                       ! change in gas concentration [#/m3]
          vMass(iH2SO4_gas) = cvap_new                              ! updating vapour concentration [mol/m3]
       
          ! volume change of particles [m3(SO4)/m3(air)] by condensation
          dvol = dvap * avogadro * volMolec(iSO4) * colrate(in1:nBins,soluble:insoluble) / csink 

          !-- Change of volume concentration of sulphate in aerosol [fxm]
          vols(iSO4,in1:nBins,soluble:insoluble) = vols(iSO4,in1:nBins,soluble:insoluble) + dvol

          !-- Change of number concentration in the smallest bin caused by nucleation
          !   Jacobson (2005), equation (16.75)
          naero(in1,soluble) = naero(in1,soluble) + n_vs_c * dvol(in1,soluble) / (volMolec(iSO4) * n3)
        endif
      else
        garbage(iH2SO4_gas) = garbage(iH2SO4_gas) + vMass(iH2SO4_gas)
        vMass(iH2SO4_gas) = 0.
      endif
    endif
call msg('H2SO4 after', vMass(iH2SO4_gas)) 
return


    !--- 5.1) Organic vapours ------------------------        
    !---- 5.1.1) Semivolatile organic compound: regimes 2 and 3
    if(iOC_sv /= int_missing)then
      if(vMass(iOC_sv) > massLowTrsh(iOC_sv))then
     
        csink = sum(colrate(in2:nBins,soluble:insoluble))      ! sink for semivolatile organics
        cvap_new = vMass(iOC_sv) / (1._r8k + tstep * csink)    ! new gas phase concentration [#/m3]
        dvap = vMass(iOC_sv) - cvap_new                        ! change in gas concentration [#/m3]
        vMass(iOC_sv) = cvap_new                               ! updating gas concentration [#/m3]
      
        ! volume change of particles 
        dvol(in2:nBins,soluble:insoluble) = dvap * avogadro * volMolec(iOC) * colrate(in2:nBins,soluble:insoluble) / csink     !  [m3(OC)/m3(air)]
        !-- change of volume due to condensation
        vols(iOC,in2:nBins,soluble:insoluble) = vols(iOC,in2:nBins,soluble:insoluble) + dvol(in2:nBins,soluble:insoluble)                      

      else
        garbage(iOC_sv) = garbage(iOC_sv) + vMass(iOC_sv)
        vMass(iOC_sv) = 0.
      endif
    endif

    !---- 5.1.1) Non-volatile organic compound: condenses onto all bins 
    if(iOC_nv /= int_missing)then
      if(vMass(iOC_nv) > massLowTrsh(iOC_nv))then
   
        cvap_new = vMass(iOC_nv) / (1._r8k + tstep * cs_tot) ! new gas phase concentration [#/m3]
        dvap = vMass(iOC_nv) - cvap_new                   ! change in gas concentration [#/m3]
        vMass(iOC_nv) = cvap_new                            ! updating vapour concentration [#/m3]
     
        ! volume change of particles 
        dvol = dVap * avogadro * volMolec(iOC) * colrate(in1:nBins,soluble:insoluble) / cs_tot !  [m3(OC)/m3(air)]
        !-- change of volume due to condensation
        vols(iOC,in1:nBins,soluble:insoluble) = vols(iOC,in1:nBins,soluble:insoluble) + dvol    
         
      else
        garbage(iOC_nv) = garbage(iOC_nv) + vMass(iOC_nv)
        vMass(iOC_nv) = 0.
      endif
    endif

	!--- 5.3) Ammoniae and nitrates -------------------------
    if(iNH3 /= int_missing .and. iHNO3 /= int_missing)then

      if(vMass(iNH3) > massLowTrsh(iNH3) .or. vMass(iHNO3) > massLowTrsh(iHNO3) .or. &
     & any(vols(iNH4,:,:) > 0.) .or. any(vols(iNO3,:,:) > 0.))then 

        ! Sums for new gas concentrations 
        aHNO3 = 0.
        bHNO3 = 0.
        aNH3 = 0.
        bNH3 = 0.
        Kelvin = 1.
     
        !--- Calling ISORROPIA for all sizebins -----------------
        WI_S = 0._r8k 
        ceqno = 0.
        ceqnh = 0.

        DO class = in1,nBins
          do iS = 1, nParallelBins
            if(naero(class,iS) < rules%nbrLowThreshold)cycle

            Kelvin(class,iS) = exp(2._r8k *srfTns_water*volMolec(iSO4)/(boltzmann_const*temp*dwet(class,iS)))

            ! WI moles/m3, aerosol phase, NA, SO4, HN4, NO3, CL
            !
            ! Seasalt: 55.03% Cl, 30.59% Na, 7.68% SO4, 6.7% everything else (by weight)
            !
            WI_S(1) = 0.3059 * vols(iSslt,class,iS) * density(iSslt) / 0.023            !Na
            WI_S(2) = vols(iSO4,class,iS) * density(iSO4) / molarMass(iSO4) + &
                    & 0.0768 * vols(iSslt,class,iS) * density(iSslt) / molarMass(iSO4)  !SO4
            WI_S(3) = vols(iNH4,class,iS) * density(iNH4) / molarMass(iNH4)             !NH4
            WI_S(4) = vols(iNO3,class,iS) * density(iNO3) / molarMass(iNO3)             !NO3
            WI_S(5) = 0.5503 * vols(iSslt,class,iS) * density(iSslt) / 0.0355           !CL
        
            CNTRL_S(1)=1  ! reverse problem solved with WI only aerosol phase 
            CNTRL_S(2)=0  ! The aerosol can have both solid+liquid phases (deliquescent)               

!            CALL ISOROPIA (WI_S, rh, temp, CNTRL_S, WT_S, GAS_S, AERLIQ_S, AERSLD_S, SCASI_S, OTHER_S)

            ! The only thing we want to use is the output GAS_S. Gaseous
            ! species concentrations, expressed in moles/m3. .GAS_S(1) - NH3
            ! GAS_S(2) - HNO3,  GAS_S(3) - HCl
            ceqno(class,iS)=GAS_S(2) * Kelvin(class,iS)
            ceqnh(class,iS)=GAS_S(1) * Kelvin(class,iS)
  
           
            ! Sums for new gas concentrations
            aHNO3 = aHNO3 - colrate(class, iS) + &
                   & colrate(class, iS) * colrate(class, iS) * ceqno(class,iS) * tstep / &
                   & (vols(iNO3,class,iS) * density(iNO3) / molarMass(iNO3) + &
                    & tstep * colrate(class, iS) * ceqno(class,iS))
         
            bHNO3 = bHNO3 + colrate(class, iS) * ceqno(class,iS) / &
                          & (1. + tstep * colrate(class, iS) * ceqno(class,iS) / &
                           & (vols(iNO3,class,iS) * density(iNO3) / molarMass(iNO3)))
        
          
            aNH3 = aNH3 - colrate(class, iS) + &
                 & colrate(class, iS) * colrate(class, iS) * ceqnh(class,iS) * tstep / &
                 & (vols(iNH4,class,iS) * density(iNH4) / molarMass(iNH4) + &
                  & tstep * colrate(class, iS) * ceqnh(class,iS))
          
            bNH3 = bNH3 + colrate(class, iS) * ceqnh(class,iS) / &
                        & (1. + tstep * colrate(class, iS) * ceqnh(class,iS) / &
                         & (vols(iNH4,class,iS) * density(iNH4) / molarMass(iNH4)))


          enddo
        END DO

        vMass(iHNO3) = (vMass(iHNO3) + bHNO3 * tstep) / (1. - aHNO3 * tstep)
        vMass(iNH3) = (vMass(iNH3) + bNH3 * tstep) / (1. - aNH3 * tstep)

        DO class = in1,nBins
          do iS = 1, nParallelBins

            vols(iNO3,class,iS) = (tstep * colrate(class,iS) * vMass(iHNO3) * avogadro * volMolec(iNO3) + vols(iNO3,class,iS)) / &
                                & (1. + tstep * colrate(class,iS) * ceqno(class,iS) / &
                                      & (vols(iNO3,class,iS) * density(iNO3) / molarMass(iNO3)))


            vols(iNH4,class,iS) = (tstep * colrate(class,iS) * vMass(iNH3) * avogadro * volMolec(iNO3) + vols(iNH4,class,iS)) / &
                                & (1. + tstep * colrate(class,iS) * ceqnh(class,iS) / &
                                      & (vols(iNH4,class,iS) * density(iNH4) / molarMass(iNH4)))
          enddo
        enddo
      else
        garbage(iNH3) = garbage(iNH3) + vMass(iNH3)
        vMass(iNH3) = 0.
        garbage(iHNO3) = garbage(iHNO3) + vMass(iHNO3)
        vMass(iHNO3) = 0.

      END IF
    endif


  contains
 
!=======================================================================================================

    !---------------------------------------------------------------------
    !
    ! subroutine BINNUCL
    ! subroutine TERNUCL
    ! subroutine ACTNUCL
    ! subroutine KINNUCL
    !
    !---------------------------------------------------------------------
    !
    ! Purpose:
    ! --------
    ! Calculate the nucleation rate and the size of critical clusters
    !
    !
    ! Method:
    ! -------  
    ! 1) Binary nucleation calculated according to 
    !  
    !  Vehkam�ki et al. (2002): An improved parameterization for
    !  sulphuric acid/water nucleation rates for tropospheric
    !  and stratospheric conditions, JGR, 107, D22, 4622.
    !
    ! 2) Ternary nucleation calculated according to
    !
    !  Napari et al. (2002), An improved model for ternary nucleation of
    !  sulfuric acid - ammonia- water, J. Chem. Phys., 116, 4221-4227
    !  
    !  Napari et al. (2002), Parameterization of ternary nucleation rates for
    !  H2SO4 - NH3 - H2O vapors, J. Geophys. Res., 107(D19), AAC 6-1
    !
    ! 3) Kinetic nucleation calculated assuming that each sulphuric acid
    !  molecule forms an (NH4)HSO4 molecule in the atmosphere and that
    !  two colliding (NH4)HSO4 molecules form a stable cluster
    !
    !
    ! Interface:
    ! ----------
    ! Called from subroutine nucleation
    ! 
    !
    ! Coded by:
    ! ---------
    ! Hannele Korhonen (FMI)      2005
    ! Hanna Vehkam�ki (University of Helsinki) 2002
    ! Ismo Napari (University of Helsinki)  2002 
    ! Sanna-Liisa Sihto (University of Helsinki) 2004
    !
    !---------------------------------------------------------------------

    SUBROUTINE binnucl(c_sa, temp, rh, nuc_rate, n_crit_sa, d_crit)

      IMPLICIT NONE

      !-- Input variables -------------------

      REAL, INTENT(IN) :: c_sa       ! sulphuric acid concentration [#/cm3]
      REAL, INTENT(IN) :: temp,   &  ! ambient temperature [K]
                          rh         ! relative humidity [0-1]

      !-- Output variables ------------------
      REAL(r8k), INTENT(OUT) :: nuc_rate,  &  ! nucleation rate [#/(m3 s)]
                               n_crit_sa, &  ! number of sulphuric acid molecules in cluster [1]
                               d_crit        ! diameter of critical cluster [m]

      !-- Local variables -------------------

      REAL(r8k) :: x,   &   ! mole fraction of sulphate in critical cluster
                  ntot,&   ! number of molecules in critical cluster
                  zt,  &   ! temperature
                  csa, &   ! sulfuric acid concentration
                  zrh, &   ! relative humidity
                  ma,mw,xmass, a, b, c, roo, m1, m2, v1, v2, zcoll


      nuc_rate  = 0._r8k
      d_crit    = 1.e-9_r8k

      !-- 1) Checking that we are in the validity range of the parameterization -----------

      zt = max(temp, 190.15_r8k)
      zt = min(zt, 300.15_r8k)

      csa = max(c_sa * avogadro * 1.e-6, 1.e4_r8k)
      csa = min(csa, 1.e11_r8k)

      zrh = max(rh, 0.0001_r8k)
      zrh = min(zrh, 1._r8k)

      ! validity of parametrization : DO NOT REMOVE!
!      IF(temp < 190.15_r8k) STOP '  INVALID INPUT VALUE (bin. nucleation): temperature < 190.15 K'
!      IF(temp > 300.15_r8k) STOP '  INVALID INPUT VALUE (bin. nucleation): temperature > 300.15 K'
!      IF(rh < 0.0001_r8k) STOP '  INVALID INPUT VALUE (bin. nucleation): relative humidity < 0.01 %'
!      IF(rh > 1._r8k) STOP '  INVALID INPUT VALUE (bin. nucleation): relative humidity > 100 %'
      !IF(c_sa < 1.e4_r8k) CYCLE! STOP '  INVALID INPUT VALUE (bin. nucleation): H2SO4 concentration < 10^4 1/cm3'
!      IF(c_sa > 1.e11_r8k) STOP '  INVALID INPUT VALUE (bin. nucleation): H2SO4 concentration > 10^11 1/cm3'

      !-- 2) Mole fraction of sulphate in a critical cluster ------------------------------
    
      x = 0.7409967177282139_r8k                    &
          - 0.002663785665140117_r8k*zt         &
          + 0.002010478847383187_r8k*log(zrh)      &
          - 0.0001832894131464668_r8k*zt*log(zrh)     &
          + 0.001574072538464286_r8k*log(zrh)**2         &  
          - 0.00001790589121766952_r8k*zt*log(zrh)**2 &
          + 0.0001844027436573778_r8k*log(zrh)**3        & 
          - 1.503452308794887e-6_r8k*zt*log(zrh)**3   &
          - 0.003499978417957668_r8k*log(csa)    &
          + 0.0000504021689382576_r8k*zt*log(csa)

      !-- 3) Nucleation rate -----------------------------------------------------------

      nuc_rate = 0.1430901615568665_r8k        &
                 + 2.219563673425199_r8k*zt        &
                 - 0.02739106114964264_r8k*zt**2         &
                 + 0.00007228107239317088_r8k*zt**3            &
                 + 5.91822263375044_r8k/x                   &
                 + 0.1174886643003278_r8k*log(zrh)                &
                 + 0.4625315047693772_r8k*zt*log(zrh)                &
                 - 0.01180591129059253_r8k*zt**2*log(zrh)            &
                 + 0.0000404196487152575_r8k*zt**3*log(zrh)          &
                 + (15.79628615047088_r8k*log(zrh))/x                   &
                 - 0.215553951893509_r8k*log(zrh)**2                    &
                 - 0.0810269192332194_r8k*zt*log(zrh)**2             &
                 + 0.001435808434184642_r8k*zt**2*log(zrh)**2        &
                 - 4.775796947178588e-6_r8k*zt**3*log(zrh)**2        &
                 - (2.912974063702185_r8k*log(zrh)**2)/x                &
                 - 3.588557942822751_r8k*log(zrh)**3                    &
                 + 0.04950795302831703_r8k*zt*log(zrh)**3            &
                 - 0.0002138195118737068_r8k*zt**2*log(zrh)**3             &
                 + 3.108005107949533e-7_r8k*zt**3*log(zrh)**3        &
                 - (0.02933332747098296_r8k*log(zrh)**3)/x              &
                 + 1.145983818561277_r8k*log(csa)                     &
                 - 0.6007956227856778_r8k*zt*log(csa)              &
                 + 0.00864244733283759_r8k*zt**2*log(csa)          &
                 - 0.00002289467254710888_r8k*zt**3*log(csa)             &
                 - (8.44984513869014_r8k*log(csa))/x                  &
                 + 2.158548369286559_r8k*log(zrh)*log(csa)            &
                 + 0.0808121412840917_r8k*zt*log(zrh)*log(csa)     &
                 - 0.0004073815255395214_r8k*zt**2*log(zrh)*log(csa)     &
                 - 4.019572560156515e-7_r8k*zt**3*log(zrh)*log(csa)      & 
                 + (0.7213255852557236_r8k*log(zrh)*log(csa))/x             &
                 + 1.62409850488771_r8k*log(zrh)**2*log(csa)          &
                 - 0.01601062035325362_r8k*zt*log(zrh)**2*log(csa)       &
                 + 0.00003771238979714162_r8k*zt**2*log(zrh)**2*log(csa) &
                 + 3.217942606371182e-8_r8k*zt**3*log(zrh)**2*log(csa)   &
                 - (0.01132550810022116_r8k*log(zrh)**2*log(csa))/x         &
                 + 9.71681713056504_r8k*log(csa)**2                   &
                 - 0.1150478558347306_r8k*zt*log(csa)**2           &
                 + 0.0001570982486038294_r8k*zt**2*log(csa)**2           &
                 + 4.009144680125015e-7_r8k*zt**3*log(csa)**2            &
                 + (0.7118597859976135_r8k*log(csa)**2)/x             &
                 - 1.056105824379897_r8k*log(zrh)*log(csa)**2         &
                 + 0.00903377584628419_r8k*zt*log(zrh)*log(csa)**2       &
                 - 0.00001984167387090606_r8k*zt**2*log(zrh)*log(csa)**2 &
                 + 2.460478196482179e-8_r8k*zt**3*log(zrh)*log(csa)**2   &
                 - (0.05790872906645181_r8k*log(zrh)*log(csa)**2)/x         &
                 - 0.1487119673397459_r8k*log(csa)**3                 &
                 + 0.002835082097822667_r8k*zt*log(csa)**3         &
                 - 9.24618825471694e-6_r8k*zt**2*log(csa)**3             &
                 + 5.004267665960894e-9_r8k*zt**3*log(csa)**3            &
                 - (0.01270805101481648_r8k*log(csa)**3)/x

      nuc_rate = exp(nuc_rate) ! [#/(cm3 s)]

      IF (nuc_rate < 1.e-7_r8k) THEN ! validity of parameterization
        nuc_rate = 0._r8k
        d_crit = 1.e-9_r8k
      END IF

      !-- 4) Total number of molecules in the critical cluster -------------------------
      ntot = - 0.002954125078716302_r8k                     &
             - 0.0976834264241286_r8k*zt                   &
             + 0.001024847927067835_r8k*zt**2                    &
             - 2.186459697726116e-6_r8k*zt**3                    &
             - 0.1017165718716887_r8k/x                       &
             - 0.002050640345231486_r8k*log(zrh)                    &
             - 0.007585041382707174_r8k*zt*log(zrh)              &
             + 0.0001926539658089536_r8k*zt**2*log(zrh)          &
             - 6.70429719683894e-7_r8k*zt**3*log(zrh)            &
             - (0.2557744774673163_r8k*log(zrh))/x                  &
             + 0.003223076552477191_r8k*log(zrh)**2                 &
             + 0.000852636632240633_r8k*zt*log(zrh)**2           &
             - 0.00001547571354871789_r8k*zt**2*log(zrh)**2            &
             + 5.666608424980593e-8_r8k*zt**3*log(zrh)**2        &
             + (0.03384437400744206_r8k*log(zrh)**2)/x              &
             + 0.04743226764572505_r8k*log(zrh)**3                  &
             - 0.0006251042204583412_r8k*zt*log(zrh)**3          &
             + 2.650663328519478e-6_r8k*zt**2*log(zrh)**3        &
             - 3.674710848763778e-9_r8k*zt**3*log(zrh)**3        &
             - (0.0002672510825259393_r8k*log(zrh)**3)/x            &
             - 0.01252108546759328_r8k*log(csa)                   &
             + 0.005806550506277202_r8k*zt*log(csa)            &
             - 0.0001016735312443444_r8k*zt**2*log(csa)        &
             + 2.881946187214505e-7_r8k*zt**3*log(csa)         &
             + (0.0942243379396279_r8k*log(csa))/x                &
             - 0.0385459592773097_r8k*log(zrh)*log(csa)           &
             - 0.0006723156277391984_r8k*zt*log(zrh)*log(csa)        &
             + 2.602884877659698e-6_r8k*zt**2*log(zrh)*log(csa)      &
             + 1.194163699688297e-8_r8k*zt**3*log(zrh)*log(csa)      &
             - (0.00851515345806281_r8k*log(zrh)*log(csa))/x            &
             - 0.01837488495738111_r8k*log(zrh)**2*log(csa)             &
             + 0.0001720723574407498_r8k*zt*log(zrh)**2*log(csa)     &
             - 3.717657974086814e-7_r8k*zt**2*log(zrh)**2*log(csa)   &
             - 5.148746022615196e-10_r8k*zt**3*log(zrh)**2*log(csa)  &
             + (0.0002686602132926594_r8k*log(zrh)**2*log(csa))/x       &
             - 0.06199739728812199_r8k*log(csa)**2                &
             + 0.000906958053583576_r8k*zt*log(csa)**2         &
             - 9.11727926129757e-7_r8k*zt**2*log(csa)**2             &
             - 5.367963396508457e-9_r8k*zt**3*log(csa)**2            &
             - (0.007742343393937707_r8k*log(csa)**2)/x                 &
             + 0.0121827103101659_r8k*log(zrh)*log(csa)**2        &
             - 0.0001066499571188091_r8k*zt*log(zrh)*log(csa)**2     &
             + 2.534598655067518e-7_r8k*zt**2*log(zrh)*log(csa)**2   &
             - 3.635186504599571e-10_r8k*zt**3*log(zrh)*log(csa)**2  &
             + (0.0006100650851863252_r8k*log(zrh)*log(csa)**2)/x       &
             + 0.0003201836700403512_r8k*log(csa)**3              &
             - 0.0000174761713262546_r8k*zt*log(csa)**3        &
             + 6.065037668052182e-8_r8k*zt**2*log(csa)**3            &
             - 1.421771723004557e-11_r8k*zt**3*log(csa)**3           &
             + (0.0001357509859501723_r8k*log(csa)**3)/x

      ntot = exp(ntot)

      !-- 5) Size of the critical cluster ----------------------------------
      n_crit_sa = x*ntot
      d_crit = 2.e-9_r8k*exp(-1.6524245_r8k+0.42316402_r8k*x+0.33466487_r8k*log(ntot)) ! [m]

      IF (n_crit_sa < 4._r8k) THEN
        !set nucleation rate to collision rate

        !volumes of the colliding objects
        ma = 96.   ![g/mol]
        mw = 18.   ![g/mol]
        xmass = 1.
        !    xmass = zxmole*ma/((1.0-zxmole)*mw + zxmole*ma) !mass fraction of h2so4
        a = 0.7681724  +xmass*(2.1847140   +xmass*(7.1630022 + &  
            xmass*(-44.31447 + xmass * (88.75606  + xmass*(-75.73729+ &
            xmass*23.43228)))))
        b = 1.808225e-3+xmass*(-9.294656e-3+xmass*(-0.03742148+& 
            xmass*(0.2565321 + xmass * (-0.5362872+ xmass*(0.4857736- &
            xmass*0.1629592)))))

        c = -3.478524e-6+xmass*(1.335867e-5+xmass*(5.195706e-5+ &
            xmass*(-3.717636e-4+xmass * (7.990811e-4+xmass*(-7.458060e-4+ &
            xmass*2.58139e-4)))))

        roo = a+zt*(b+c*zt) ! g/cm^3
        roo = roo*1.e+3_r8k !kg/m^3

        m1 = 0.098 ! [kg/mol]/ Na !kg
        m2 = m1
        v1 = m1/avogadro/roo
        v2 = v1

        zcoll = csa*csa* &
                (3._r8k*pi/4._r8k)**(1._r8k/6._r8k) * &
                SQRT(6._r8k*gas_constant_uni*zt/m1+6._r8k*gas_constant_uni*zt/m2) * &
                (v1**(1._r8k/3._r8k) + v2**(1._r8k/3._r8k))**2. * &
                1.e+6_r8k        ! m3 -> cm3      

        zcoll=MIN(zcoll,1.e10_r8k)

        nuc_rate=zcoll

      ELSE

        nuc_rate = min(nuc_rate,1.e10_r8k)

      END IF

      nuc_rate = nuc_rate*1.e6_r8k ! [#/(m3 s)]

    END SUBROUTINE binnucl


!    !===================================================================================
!
!    SUBROUTINE ternucl(csa, c_nh3, temp, rh, nuc_rate, n_crit_sa, d_crit)     
!
!      IMPLICIT NONE
!
!      !-- Input variables -------------------
! 
!      REAL, INTENT(IN) :: csa,      & ! sulphuric acid concentration [#/cm3]
!                          c_nh3        ! ammonia mixing ratio [ppt]              !???????????????????????????????????????
!      REAL, INTENT(IN) :: temp,      & ! ambient temperature [K]
!                          rh     ! relative humidity
!
!      !-- Output variables -------------------
!      REAL(r8k), INTENT(OUT) :: nuc_rate,  & ! nucleation rate [#/(m3 s)]
!                               n_crit_sa, & ! number of H2SO4 molecules in cluster [1]
!                               d_crit       ! diameter of critical cluster [m]
!      REAL(r8k) :: lnj     ! logarithm of nucleation rate
!      REAL(r8k) :: c_sa
!
!      c_sa = csa * avogadro * 1.e-6
!
!      call set_error('ammonia mixing ratio [ppt]???','ternucl')
!      call set_error('A-a-a-ahhhh, this crap crushes the compiler !!!!','ternucl')
!      return
!
!      !-- 1) Checking that we are in the validity range of the parameterization -----------
!
!      ! validity of parameterization : DO NOT REMOVE!
!      IF (temp < 240._r8k)  STOP '  INVALID INPUT VALUE (ter. nucleation): temperature < 240 K'
!      IF (temp > 300._r8k)  STOP '  INVALID INPUT VALUE (ter. nucleation) temperature > 300 K'
!      IF (rh < 0.05_r8k)    STOP '  INVALID INPUT VALUE (ter. nucleation) relative humidity < 5 %'
!      IF (rh > 0.95_r8k)    STOP '  INVALID INPUT VALUE (ter. nucleation) relative humidity > 95 %'
!      IF (c_sa < 1.e4_r8k)  STOP '  INVALID INPUT VALUE (ter. nucleation) H2SO4 concentration < 10^4 1/cm3'
!      IF (c_sa > 1.e9_r8k)  STOP '  INVALID INPUT VALUE (ter. nucleation) H2SO4 concentration > 10^9 1/cm3'
!      IF (c_nh3 < 0.1_r8k)  STOP '  INVALID INPUT VALUE (ter. nucleation) ammonia mixing ratio < 0.1 ppt'
!      IF (c_nh3 > 100._r8k) STOP '  INVALID INPUT VALUE (ter. nucleation) ammonia mixing ratio > 100 ppt'
!
!  
!      !-- 2) Nucleation rate ---------------------------------------------------
!
!      lnj = -84.7551114741543_r8k + 0.3117595133628944_r8k*rh                &
!            + 1.640089605712946_r8k*rh*temp                  &
!            - 0.003438516933381083_r8k*rh*temp**2                  &
!            - 0.00001097530402419113_r8k*rh*temp**3                &
!            - 0.3552967070274677_r8k/log(c_sa)                      &
!            - (0.06651397829765026_r8k*rh)/log(c_sa)               &
!            - (33.84493989762471_r8k*temp)/log(c_sa)               &
!            - (7.823815852128623_r8k*rh*temp)/log(c_sa)           &
!            + (0.3453602302090915_r8k*temp**2)/log(c_sa)                 &
!            + (0.01229375748100015_r8k*rh*temp**2)/log(c_sa)            &
!            - (0.000824007160514956_r8k*temp**3)/log(c_sa)               &
!            + (0.00006185539100670249_r8k*rh*temp**3)/log(c_sa)         & 
!            + 3.137345238574998_r8k*log(c_sa)                       &
!            + 3.680240980277051_r8k*rh*log(c_sa)                   &
!            - 0.7728606202085936_r8k*temp*log(c_sa)                &
!            - 0.204098217156962_r8k*rh*temp*log(c_sa)             &
!            + 0.005612037586790018_r8k*temp**2*log(c_sa)                 &
!            + 0.001062588391907444_r8k*rh*temp**2*log(c_sa)             &
!            - 9.74575691760229e-6_r8k*temp**3*log(c_sa)            &
!            - 1.265595265137352e-6_r8k*rh*temp**3*log(c_sa)             &
!            + 19.03593713032114_r8k*log(c_sa)**2                    &
!            - 0.1709570721236754_r8k*temp*log(c_sa)**2             &
!            + 0.000479808018162089_r8k*temp**2*log(c_sa)**2              &
!            - 4.146989369117246e-7_r8k*temp**3*log(c_sa)**2              &
!            + 1.076046750412183_r8k*log(c_nh3)                      &
!            + 0.6587399318567337_r8k*rh*log(c_nh3)                 &
!            + 1.48932164750748_r8k*temp*log(c_nh3)                 & 
!            + 0.1905424394695381_r8k*rh*temp*log(c_nh3)           &
!            - 0.007960522921316015_r8k*temp**2*log(c_nh3)                &
!            - 0.001657184248661241_r8k*rh*temp**2*log(c_nh3)            &
!            + 7.612287245047392e-6_r8k*temp**3*log(c_nh3)                &
!            + 3.417436525881869e-6_r8k*rh*temp**3*log(c_nh3)            & 
!            + (0.1655358260404061_r8k*log(c_nh3))/log(c_sa)              & 
!            + (0.05301667612522116_r8k*rh*log(c_nh3))/log(c_sa)         &
!            + (3.26622914116752_r8k*temp*log(c_nh3))/log(c_sa)          &
!            - (1.988145079742164_r8k*rh*temp*log(c_nh3))/log(c_sa)     &
!            - (0.04897027401984064_r8k*temp**2*log(c_nh3))/log(c_sa)    &
!            + (0.01578269253599732_r8k*rh*temp**2*log(c_nh3))/log(c_sa)&
!            + (0.0001469672236351303_r8k*temp**3*log(c_nh3))/log(c_sa)        &
!            - (0.00002935642836387197_r8k*rh*temp**3*log(c_nh3))/log(c_sa) &
!            + 6.526451177887659_r8k*log(c_sa)*log(c_nh3)                 & 
!            - 0.2580021816722099_r8k*temp*log(c_sa)*log(c_nh3)          &
!            + 0.001434563104474292_r8k*temp**2*log(c_sa)*log(c_nh3)     &
!            - 2.020361939304473e-6_r8k*temp**3*log(c_sa)*log(c_nh3)     &
!            - 0.160335824596627_r8k*log(c_sa)**2*log(c_nh3)              &
!            + 0.00889880721460806_r8k*temp*log(c_sa)**2*log(c_nh3)      &
!            - 0.00005395139051155007_r8k*temp**2*log(c_sa)**2*log(c_nh3)      &
!            + 8.39521718689596e-8_r8k*temp**3*log(c_sa)**2*log(c_nh3)         &
!            + 6.091597586754857_r8k*log(c_nh3)**2                   &
!            + 8.5786763679309_r8k*rh*log(c_nh3)**2                 &
!            - 1.253783854872055_r8k*temp*log(c_nh3)**2             &
!            - 0.1123577232346848_r8k*rh*temp*log(c_nh3)**2        &
!            + 0.00939835595219825_r8k*temp**2*log(c_nh3)**2              &
!            + 0.0004726256283031513_r8k*rh*temp**2*log(c_nh3)**2        &
!            - 0.00001749269360523252_r8k*temp**3*log(c_nh3)**2           &
!            - 6.483647863710339e-7_r8k*rh*temp**3*log(c_nh3)**2         &
!            +(0.7284285726576598_r8k*log(c_nh3)**2)/log(c_sa)            &
!            + (3.647355600846383_r8k*temp*log(c_nh3)**2)/log(c_sa)      &
!            - (0.02742195276078021_r8k*temp**2*log(c_nh3)**2)/log(c_sa)       &
!            + (0.00004934777934047135_r8k*temp**3*log(c_nh3)**2)/log(c_sa)    &
!            + 41.30162491567873_r8k*log(c_sa)*log(c_nh3)**2              &
!            - 0.357520416800604_r8k*temp*log(c_sa)*log(c_nh3)**2        &
!            + 0.000904383005178356_r8k*temp**2*log(c_sa)*log(c_nh3)**2        &
!            - 5.737876676408978e-7_r8k*temp**3*log(c_sa)*log(c_nh3)**2        &
!            - 2.327363918851818_r8k*log(c_sa)**2*log(c_nh3)**2           &
!            + 0.02346464261919324_r8k*temp*log(c_sa)**2*log(c_nh3)**2         &
!            - 0.000076518969516405_r8k*temp**2*log(c_sa)**2*log(c_nh3)**2     &
!            + 8.04589834836395e-8_r8k*temp**3*log(c_sa)**2*log(c_nh3)**2      &
!            - 0.02007379204248076_r8k*log(rh)                       &
!            - 0.7521152446208771_r8k*temp*log(rh)                  &
!            + 0.005258130151226247_r8k*temp**2*log(rh)             &
!            - 8.98037634284419e-6_r8k*temp**3*log(rh)              &
!            + (0.05993213079516759_r8k*log(rh))/log(c_sa)                &
!            + (5.964746463184173_r8k*temp*log(rh))/log(c_sa)            &
!            - (0.03624322255690942_r8k*temp**2*log(rh))/log(c_sa)       &
!            + (0.00004933369382462509_r8k*temp**3*log(rh))/log(c_sa)    &
!            - 0.7327310805365114_r8k*log(c_nh3)*log(rh)            &
!            - 0.01841792282958795_r8k*temp*log(c_nh3)*log(rh)           &
!            + 0.0001471855981005184_r8k*temp**2*log(c_nh3)*log(rh)      &
!            - 2.377113195631848e-7_r8k*temp**3*log(c_nh3)*log(rh)
!
!      nuc_rate = exp(lnj) ! [#/(cm3 s)]
!
!      IF (nuc_rate < 1.e-5_r8k) THEN ! validity of parametrization
!        nuc_rate = 0._r8k
!        d_crit = 1.e-9_r8k
!        RETURN
!      END IF
!
!      ! validity of parametrization
!      IF (nuc_rate > 1.e6_r8k) STOP '  INVALID OUTPUT VALUE for ternary nucleation: nucleation rate > 10^6 1/cm3s'
!      nuc_rate = nuc_rate*1.e6_r8k ! [#/(m3 s)]
!
!      !-- 3) Number of H2SO4 molecules in a critical cluster -----------------------
!      n_crit_sa = 38.16448247950508_r8k + 0.7741058259731187_r8k*lnj +             &
!                  0.002988789927230632_r8k*lnj**2 - 0.3576046920535017_r8k*temp -          &
!                  0.003663583011953248_r8k*lnj*temp + 0.000855300153372776_r8k*temp**2
!      n_crit_sa = MAX(n_crit_sa,2.e0_r8k) ! kinetic limit: at least 2 H2SO4 molecules in a cluster
!
!      !-- 4) Size of the critical cluster ------------------------------------------
!      d_crit = 0.1410271086638381_r8k - 0.001226253898894878_r8k*lnj - &
!               7.822111731550752e-6_r8k*lnj**2 - 0.001567273351921166_r8k*temp - &
!               0.00003075996088273962_r8k*lnj*temp + 0.00001083754117202233_r8k*temp**2  ! [nm]
!      d_crit = d_crit*2.e-9_r8k  ! [m]
!
!    END SUBROUTINE ternucl

    !===================================================================================


    SUBROUTINE actnucl(sa_conc,  nuc_rate, d_crit)

      IMPLICIT NONE

      REAL, INTENT(IN) ::  sa_conc     ! sulphuric acid concentration [#/m3]

      REAL(r8k), INTENT(OUT) :: nuc_rate, & ! nucleation rate [#/(m3 s)]
                               d_crit      ! critical diameter of clusters [m]

      nuc_rate = 1.e-7_r8k * sa_conc * avogadro ! [#/(m3 s)]
      d_crit = 7.9375e-10_r8k ! [m]

    END SUBROUTINE actnucl

    !===================================================================================


    SUBROUTINE kinnucl(sa_conc, temp, nuc_rate, d_crit)

      ! 
      ! Below the following assumption have been made:W
      !
      !  nucrate = coagcoeff*csa**2
      !  coagcoeff = 8*sqrt(3*boltz*temp*r_abs/dens_abs)
      !  r_abs = 0.315d-9 radius of bisulphate molecule [m]
      !  dens_abs = 1465  density of - " - [kg/m3]
      !

      IMPLICIT NONE

      !-- Input variables -------------------
        
      REAL, INTENT(IN) :: sa_conc ! sulphuric acid concentration [#/m3]
      REAL, INTENT(IN) :: temp        ! ambient temperature [K]

      REAL(r8k), INTENT(OUT) ::  nuc_rate, &  ! nucleation rate [#/(m3 s)]
                                d_crit   ! critical diameter of clusters [m]
      real(r8k) :: csa

      nuc_rate = 1.19386e-17_r8k * sqrt(temp) * (sa_conc * avogadro)**2 ! [#/(m3 s)]
      d_crit = 7.9375e-10_r8k ! [m]

    END SUBROUTINE kinnucl

  END SUBROUTINE condensation


  !***********************************************************************

  logical function fu_if_tla_required_ADB(rules) result(required)
    ! Collect transformations' requests for tangent linearization. If tangent linear
    ! is not allowed, corresponding subroutine sets error. 
    implicit none
    type(Tchem_rules_AerDynBasic), intent(out) :: rules

    call set_error('Aerosol dynamics basic does not have TLA/ADJOINT','fu_if_tla_required_ADB')
    
     required = .false.
    
  end function fu_if_tla_required_ADB

END MODULE aer_dyn_basic



