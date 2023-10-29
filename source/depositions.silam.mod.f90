module depositions
  !
  ! A module for standard deposition calculations. At least dry deposition
  ! rates (Rs,Rb,etc) and settling should be defined here, as well as 
  ! the wet deposition. All the necessary information is included in a 
  ! species structure, which can be given as arguments to the subroutines.
  !
  ! The deposition setup should come from two places: material should
  ! include material-specific data, while general options can be
  ! stored in the deposition_rules structure.
  !
  ! Should the non-standard deposition be used for some species, it has to 
  ! be either defined here and used via the deposition_rules or defined in
  ! the corresponding transformation modules. The later decision, however,
  ! should have a reason, such as the compatibility requirement with the
  ! corresponding chemistry.
  !
  ! Author: M.Sofiev, FMI, R.Kouznetsov, FMI
  !
  ! Language: ANSI FORTRAN-90 (or close to it)
  !
  ! Units: SI
  !
  use lagrange_particles

  ! to be able to share the metdat column indices:
  use photolysis

  implicit none

!  private

!!!!!!!!!!!!!!!!#define DEBUG_TLA
  
  public set_deposition_rules
  public init_standard_deposition
  public add_deposition_input_needs
  public wet_deposition_input_needs
  public prepare_deposition

!  public get_Rb
  public get_Rs
  public get_Rs_2013
  public fu_get_vd
!    public check_deposition_rules
  public get_settling_motion_species
  public make_deposition_lagr
  public scavenge_column
  public verify_deposition
  public allocate_scav_amount
  public fu_if_TLA_required

  private scavenge_lagr
  private scavenge_puff_2011
  private scavenge_column_2011_adjoint
  private fu_if_TLA_required_scav
  private material_properties
  private test_TL_ADJ
  private make_rain_in_cell
  private make_rain_in_cell_col
  private report_rain
  private fu_vdplus_slinnSP98
  private get_vd_species
  private get_settling_velocity_species
  private iVrk4
  private DERIV
  private DERIV2
  private fu_Psi
  private fu_vdplus_rough
  
  interface fu_if_TLA_required
    module procedure fu_if_TLA_required_scav
  end interface

  ! max_scav_rate_depends on one of these
  integer, public, parameter :: horiz_wind = 7001
  integer, public, parameter :: cape = 7002
  integer, public, parameter :: cape_and_wind = 7003

  type Tdeposition_rules
    integer, dimension(:), pointer :: iDepositionType
    integer :: nDepositionTypes
    type(silam_species), dimension(:), pointer :: pSpecies
    integer, dimension(:), pointer :: indGasesDepositing, indAerosolDepositing ! grouping the species
    integer :: nGasesDepositing, nAerosolsDepositing
    integer :: scavengingType ! Method for computing the scavenging
    integer :: DryDepType     ! Method for computing the dry deposition
    integer :: RsType         ! Method for surface resistance for gases
    integer :: max_scav_rate_depends_on = int_missing
    real :: max_scav_rate_cape_scaling = 0.5, max_scav_rate_wind_scaling = 0.3
    logical :: ifHumidityDependent ! or just take the defualt humidity
    real :: fDefaultRelHumidity       ! the very default humidity
    logical :: ifSulphurSaturation
    logical :: if_lagr = .false.
    real :: fSulphurSaturation
    type(silja_logical) :: defined
  end type Tdeposition_rules
  public Tdeposition_rules
  
  ! Private placeholders for rain features in the cell and column
  !
  type Train_in_cell ! type to hold  wet-dep rlated parameters of a grid-cell
                     ! and params needed to diagnose Train_in_cell below 
    logical :: ifRainyCell, ifRainSfc, ifRealCloud, ifValid=.FALSE.
    integer :: iMeteo, ix, iy, iLev
    real :: fCloudCover, fPressure, fTemperature, fRainRateSfc, cwcColumn, pwcColumn ! kg/m2
                !Generic meteo
    real :: max_scav_rate  !!!1./sec
    real :: mu_air, rho_air, nu_air, lambda_air !air parameters
    real :: rain_rate_in, deltarain, & !rain rate from above 
                   cwcColumnLayer, pwcColumnLayer !  kg/m2
    real :: is_liquid, drop_size, drop_mass, drop_vel, drop_Re !droplets/snowflakes
    real :: cell_zsize, cell_volume
  end type Train_in_cell
  private Train_in_cell
  !
  ! Private palceholders for internal scavenging coefficients in cell and column
  !
  type Tscav_coefs
    real, dimension(max_species) :: scav_coef_ic, scav_coef_sc, scav_coef_tau, A, tauf_over_RdCd
    real :: fB, fC, Ceq, dLsc_dWc, dLic_dWc, dCeq_dMs, dWC_dMs, dCeq_dMa, dWC_dMa
  end type Tscav_coefs
  private Tscav_coefs
  type Tscav_coefs_1d
    type(Tscav_coefs), dimension(max_levels) :: c
  end type Tscav_coefs_1d
  private Tscav_coefs_1d
  type(Tscav_coefs), private, parameter :: coefs_zero = Tscav_coefs(0.,0.,0.,0.,0., &
                                                                  & 0.,0.,0.,0.,0.,0.,0.,0.,0.)

  ! For each grid column, we will have to collect the scavenging amount into a single array 
  !
  real, dimension(:,:,:,:), private, pointer :: scavAmount  ! (nSpecies, nSrc, nLevels, nthreads)

  ! Standard deposition flag
  !
  integer, public, parameter :: deposition_standard = 5100

  real, public, parameter :: freezing_point_of_water = 263

  ! Scavenging computation allowed so far
  integer, public, parameter :: scavStandard = 6101 ! Standard 3D scavenging ratio

  integer, public, parameter :: scav2011 = 6112     ! Scavenge accounting
                                                     ! for rain profile,
                                                     ! saturations etc..
  integer, public, parameter :: scav2011fc = 6113   ! Same as above, but with fake cloud profile
  integer, public, parameter :: scavNoScav = 6114     ! No scavenging at all
  integer, public, parameter :: scav2018 = 6115     ! #Limited washout
  integer, public, parameter :: scav2018entr = 6116     ! #Limited washout with proper entrainment
  integer, public, parameter :: scav2020 = 6117    !the precipitable water is based on a cloud droplet model

  ! Dry deposition computation allowed so far
!  integer, public, parameter :: DryD_none = 6200 ! No dry deposition at all.
  ! Resistance-based schemes
!  integer, private, parameter :: DryD_Rb_only = 6201 ! Only Rb of resistive scheme 
!  integer, private, parameter :: DryD_resist_only = 6202 ! Only resistive scheme, all R-s 
!  integer, private, parameter :: DryD_grav_only = 6203 ! Only settling
!  integer, private, parameter :: DryD_Rb_grav_combine = 6204 ! Rb and settling
!  integer, private, parameter :: DryD_resist_grav_combine = 6205 ! Both resist and settling

! for schemes >= DryD_KS2011  Vd is calculated  
  integer, public, parameter :: DryD_KS2011 = 6266 ! Our deposition scheme
  integer, public, parameter :: DryD_KS2011_TF = 6267 ! Our deposition scheme
                                                      ! with thermophoresis
  integer, public, parameter :: DryD_Vd_difsed = 6268 ! Diffusion-settling
  integer, public, parameter :: DryD_Vd_SP98 = 6269 ! Slinn -- after SP98
                                                      ! as it was in Silam
  integer, public, parameter :: DryD_Vd_Zhang = 6271 ! Zhang -- after SP2006
  integer, public, parameter :: DryD_Vd_Zero = 6272 ! Vd=0


  integer, public, parameter :: DryD_Rs_standard = 6300 ! Standard land-sea Rs
  integer, public, parameter :: DryD_Rs_2013     = 6301 ! Fancy Rs                                                      

! Types of Zhang parameterization 
! Stub. To be removed at some point..
  integer, public, parameter :: ZhangGeneric = 10000
  integer, public, parameter :: ZhangSmooth = 10001
  integer, public, parameter :: ZhangGrass = 10002
  integer, public, parameter :: ZhangSnow = 10003
  integer, public, parameter :: ZhangBlforest = 10004



!****************************************************
  !-------------------------------------------------------------------------
  !
  ! Pointers to the meteorological fields, which can be requested by the deposition routines
  !
  ! Surface fields, meteorological grid: via 1-D pointers
  !
  real , dimension(:), private, pointer, save :: pMetLandFr => null(), pMetTempr2m => null(), pMetRh2m => null(), pMetSrfPressure => null(), &
                                              & pMetSrfRoughMeteo => null(), pMetSrfRoughDisp => null(), &
                                              & pMetPrecTot => null(), pMetPrecLs => null(), pMetPrecCnv => null(), pMetTotCloud => null(), &
                                              & pMetABLHeight => null(), pMetFricVel => null(),  pMetMO_Len_inv => null(), &
                                              & pMetConvVel => null(), pMetSensHF => null(), pMetSnowDepth => null(), pMetLAI => null(), &
                                              & pMetGsto => null(), pMetIcefr => null(), pCAPE => null()

  integer, private, pointer, save :: imet_cwc3d, imet_pwc3d, imet_press, imet_landfrac,   &
       & imet_temp, imet_u, imet_v, imet_cape, imet_tcc, imet_scav,  imet_dx_size, imet_dy_size, &
       & imet_dz_size, imet_abl, imet_airdens, imet_airmass, imet_cc3d, imet_cic, imet_cwc, imet_prec

  real, private :: meteo_cell_size !! Linear size of a meteo cell, m
  !
  ! Surface fields, dispersion grid: via 1-D pointers
  !
!  real, dimension(:), private, pointer, save :: pDispVdCorrection
  !
  ! 3-D fields have to be given via field_4d_data_ptr structures because they may not have
  ! present value
  !
  type(field_4d_data_ptr), private, pointer, save :: fldScavStd, fldPrecip, fldTempr, fldRelHumidity, fldPressure, &
                                & fldCellSizeX, fldCellSizeY, fldCellSizeZ, fldCwcabove,  fldPwcabove, fldUmet, fldVmet

  !  type(field_4d_data_ptr), private, pointer, save :: fldDispZsize

  integer, dimension(:), private, pointer :: indSulphurSpecies, indStrongAcidSpecies, &
                                           & indAcidityAffectingSp, iAciditySign  ! acid and alcaline species
  integer, private :: nSulphurSpecies, indSO2, nStrongAcidSpecies, nAcidityAffectingSpecies, &
                    & indNH3, indHNO3, indO3



CONTAINS

  !*****************************************************************************
  !
  !  THE DEPOSITION MANAGER
  !
  !*****************************************************************************

  subroutine set_deposition_rules(nlSetup, rulesDeposition, if_lagr_present)
    !
    ! Creates the set of standard deposition rules. The input is taken
    ! from the transformation_parameters namelist in the control
    ! file. The option for non-standard deposition taken by each
    ! module is left for the future initialisation
    !
    implicit none

    ! Imported parameter
    type(Tsilam_namelist), intent(in) :: nlSetup 
    type(Tdeposition_rules) :: rulesDeposition
    logical, intent(in) :: if_lagr_present

    ! Local parameters
    integer :: iTmp, iTransf
    character (len=*), parameter :: sub_name="set_deposition_rules"

    if(.not.defined(nlSetup))then
      call set_error('Undefined namelist given',sub_name)
      rulesDeposition%defined = silja_false
      return
    endif

    rulesDeposition%if_lagr = if_lagr_present 
    !
    ! Scavenging type. May be undefined - then take the default one
    !
    select case(trim(fu_str_u_case(fu_content(nlSetup,'wet_deposition_scheme'))))
      case('STANDARD_3D_SCAVENGING')
        rulesDeposition%scavengingType = scavStandard

      case('NEW2011_SCAVENGING')
        rulesDeposition%scavengingType = scav2011

      case('2018_SCAVENGING')
        rulesDeposition%scavengingType = scav2018

      case('2018_SCAVENGING_ENTRAINMENT')
        rulesDeposition%scavengingType = scav2018entr

      case('2020_SCAVENGING')
        rulesDeposition%scavengingType = scav2020

      case('NEW2011_SCAVENGING_FAKECLOUD')
        rulesDeposition%scavengingType = scav2011fc

      case('NO_SCAVENGING')
        rulesDeposition%scavengingType = scavNoScav

      case default
        call set_error('Unknown scavenging type: "'// trim(fu_content(nlSetup,'wet_deposition_scheme')) // '"',&
                     & sub_name)
        call msg("Valid values for wet_deposition_scheme are:")
        call msg(" (STANDARD_3D_SCAVENGING|NEW2011_SCAVENGING|NEW2011_SCAVENGING_FAKECLOUD|NO_SCAVENGING|2020_SCAVENGING|2018_SCAVENGING)")
        call msg(" 2019_SCAVENGING depricated since r588934")
        rulesDeposition%defined = silja_false
        return
    end select


    if (fu_str_u_case(fu_content(nlSetup,'max_scav_rate_depends_on')) == 'HORIZONTAL_WIND') then
      call msg('Max scavenging rate depends on the horizontal wind')
      rulesDeposition%max_scav_rate_depends_on = horiz_wind
    else if (fu_str_u_case(fu_content(nlSetup,'max_scav_rate_depends_on')) == 'CAPE') then
      call msg('Max scavenging rate depends on the convective available potential energy')
      rulesDeposition%max_scav_rate_depends_on = cape
    else if (fu_str_u_case(fu_content(nlSetup,'max_scav_rate_depends_on')) == 'CAPE_AND_WIND') then
      call msg('Max scavenging rate depends on the convective available potential energy and the horizontal wind')
      rulesDeposition%max_scav_rate_depends_on = cape_and_wind
    end if
    
    if (rulesDeposition%max_scav_rate_depends_on == int_missing .neqv. &
       & any(rulesDeposition%scavengingType == (/scavStandard, scav2011, scav2011fc, scavNoScav/))) then
       call set_error("Incompatible max_scav_rate_depends_on and wet_deposition_scheme", sub_name)
       return
    endif

    if (rulesDeposition%max_scav_rate_depends_on == cape_and_wind) then
      rulesDeposition%max_scav_rate_wind_scaling = fu_content_real(nlSetup,'max_scav_rate_wind_scaling')
      if(error .or. rulesDeposition%max_scav_rate_wind_scaling < 0.0)then
        call msg('strange max_scav_rate_wind_scaling')
        call set_error('Strange or missing max_scav_rate_wind_scaling. It could be 1.0 for a global ERA5 run or 0.3 for &
             &a European run with EC oper',sub_name)
        return
      endif
      rulesDeposition%max_scav_rate_cape_scaling = fu_content_real(nlSetup,'max_scav_rate_cape_scaling')
      if(error .or. rulesDeposition%max_scav_rate_wind_scaling < 0.0)then
        call msg('strange max_scav_rate_cape_scaling')
        call set_error('Strange or missing max_scav_rate_cape_scaling. It could be 1.0 for a global ERA5 run or 0.4 &
             &for a European run with EC oper.',sub_name)
        return
      endif
    end if
    
    !
    ! Saturation of sulphur scavenging is so far the only non-linear species-specific process
    !
    rulesDeposition%ifSulphurSaturation = &
               & fu_str_u_case(fu_content(nlSetup,'if_sulphur_scavenging_saturation')) == 'YES'
    rulesDeposition%fSulphurSaturation = 2.3e-2  ! about 750 mg S / l
 
    !
    ! Dry deposition type. May be undefined - then take the default one.
    !
!call report(nlSetup)
    select case(fu_str_u_case(fu_content(nlSetup,'dry_deposition_scheme')))
!      case('SIMPLE_DIFFUSION_ONLY')
!        rulesDeposition%DryDepType = DryD_Rb_only
!
!      case('FULL_DIFFUSION_ONLY')
!        rulesDeposition%DryDepType = DryD_resist_only
!
!      case('GRAVITATIONAL_ONLY')
!        rulesDeposition%DryDepType = DryD_grav_only
!
!      case('GRAVITATIONAL_AND_SIMPLE_DIFFUSION')
!        rulesDeposition%DryDepType = DryD_Rb_grav_combine
!
!      case('GRAVITATIONAL_AND_FULL_DIFFUSION')
!        rulesDeposition%DryDepType = DryD_resist_grav_combine

      case('KS2011')
        rulesDeposition%DryDepType = DryD_KS2011

      case('KS2011_TF')
        rulesDeposition%DryDepType = DryD_KS2011_TF
      
      case('VD_DIFSED')
        rulesDeposition%DryDepType = DryD_Vd_difsed

      case('NO_DD') 
         rulesDeposition%DryDepType = DryD_Vd_Zero

      case('VD_SP98') 
         rulesDeposition%DryDepType = DryD_Vd_SP98


      case('') ! No such line at all
        call msg_warning('Dry deposition is not defined. Take default',sub_name)
        rulesDeposition%DryDepType = DryD_KS2011_TF

      case default
        call set_error('Unknown dry deposition type:' + &
                                        & fu_content(nlSetup,'dry_deposition_scheme'), &
                     & sub_name)
        rulesDeposition%defined = silja_false
        return
    end select

    select case(fu_str_u_case(fu_content(nlSetup,'surface_resistance_method')))

      case('STANDARD')
        rulesDeposition%RsType = DryD_Rs_standard
        call msg('Using SILAM standard Rs for dry depo')

      case('WES2013')
        rulesDeposition%RsType = DryD_Rs_2013
        call msg('Using Wesely-type Rs for dry depo')

      case('') ! No such line at all
        call msg_warning('Rs for dry deposition is not defined. Using standard',sub_name)
        rulesDeposition%RsType = DryD_Rs_standard

      case default
              call set_error('Unknown Rs methos dry deposition type:' + &
                                        & fu_content(nlSetup,'surface_resistance_method'), &
                     & sub_name)
        rulesDeposition%defined = silja_false
        return
    end select

    !
    ! In theory, aerosol features depend on humidity. However, in simple cases we can ignore them
    ! just setting the default relative humidity value = 80%
    !
    if(fu_str_u_case(fu_content(nlSetup,'if_actual_humidity_for_particle_size')) == 'YES')then
      rulesDeposition%ifHumidityDependent = .true.
    elseif(fu_str_u_case(fu_content(nlSetup,'if_actual_humidity_for_particle_size')) == 'NO')then
      rulesDeposition%ifHumidityDependent = .false.
    else
      call set_error('if_actual_humidity_for_particle_size must be YES or NO',sub_name)
      return
    endif

    if(.not. rulesDeposition%ifHumidityDependent)then
      rulesDeposition%fDefaultRelHumidity = fu_content_real(nlSetup,'default_relative_humidity')
      if(error .or. rulesDeposition%fDefaultRelHumidity < 0 .or. rulesDeposition%fDefaultRelHumidity > 1.)then
        call msg('if_actual_humidity_for_particle_size = NO but strange default_relative_humidity')
        call report(nlSetup)
        call set_error('default_relative_humidity must be real [0..1]',sub_name)
        return
      endif
    endif

    rulesDeposition%defined = silja_true

  end subroutine set_deposition_rules


  !**********************************************************************************

  subroutine init_standard_deposition(speciesTransport, indDepositionType, nSpeciesTransport, &
                                    & rulesDeposition)
    !
    ! Essentially, distributes the species among the deposition procedures,
    ! Sets indices for some special species and checks if they have needed parameters set.
    ! For the standard deposition computed here, we have two groups: aerosols and gases
    ! They need somewhat different input parameters, etc, so it is easier to treat them
    ! together.
    !
    implicit none

    ! Imported parameters
    type(silam_species), dimension(:), pointer :: speciesTransport
    integer, dimension(:), intent(inout) :: indDepositionType
    integer, intent(in) :: nSpeciesTransport
    type(Tdeposition_rules), intent(inout) :: rulesDeposition

    ! Local variables
    type(Tgas_deposition_param) :: DepData
    integer :: iSpecies, iSpTr
    integer :: iTmp
    integer, dimension(:), pointer :: iarTmp, iarTmpSa, iarTmpAAS
    !
    ! Initialisation means that we have to mark the species that will be treated via the standard
    ! deposition
    !
    rulesDeposition%nGasesDepositing = 0
    rulesDeposition%nAerosolsDepositing = 0
    allocate(rulesDeposition%indGasesDepositing(nSpeciesTransport), &
           & rulesDeposition%indAerosolDepositing(nSpeciesTransport), &
           & stat=iTmp)
    if(iTmp /= 0)then
      call set_error('Failed to allocate the depositing species arrays','init_standard_deposition')
      return
    endif
    
    iarTmp => fu_work_int_array()
    iarTmpSa => fu_work_int_array()
    iarTmpAAS => fu_work_int_array()
    if(error)return
    nSulphurSpecies = 0
    nStrongAcidSpecies = 0
    nAcidityAffectingSpecies = 0

    indSO2 = -1 
    indNH3 = -1
    indHNO3 = -1
    indO3 = -1
    do iTmp = 1, nSpeciesTransport
      if(indDepositionType(iTmp) == int_missing)then
        indDepositionType(iTmp) = deposition_standard
        if(fu_mode(speciesTransport(iTmp)) == in_gas_phase)then
          if(fu_if_gas_depositing(speciesTransport(iTmp)%material))then
            rulesDeposition%nGasesDepositing = rulesDeposition%nGasesDepositing + 1
            rulesDeposition%indGasesDepositing(rulesDeposition%nGasesDepositing) = iTmp
          endif
        else
          rulesDeposition%nAerosolsDepositing = rulesDeposition%nAerosolsDepositing + 1
          rulesDeposition%indAerosolDepositing(rulesDeposition%nAerosolsDepositing) = iTmp
        endif
      endif
      select case(fu_str_u_case(fu_name(fu_material(speciesTransport(iTmp)))))
        case('NH3')
          indNH3 = iTmp
          nAcidityAffectingSpecies = nAcidityAffectingSpecies + 1
          iarTmpAAS(nAcidityAffectingSpecies) = -iTmp   ! it is alcaline, mark with the sign
        case('O3','NOP')
          indO3 = iTmp
        case ('HNO3')
          nStrongAcidSpecies = nStrongAcidSpecies + 1
          nAcidityAffectingSpecies = nAcidityAffectingSpecies + 1
          iarTmpSa(nStrongAcidSpecies) = iTmp
          iarTmpAAS(nAcidityAffectingSpecies) = iTmp
          indHNO3 = iTmp
        case('SO2')
          nSulphurSpecies = nSulphurSpecies + 1
          nAcidityAffectingSpecies = nAcidityAffectingSpecies + 1
          iarTmp(nSulphurSpecies) = iTmp
          iarTmpAAS(nAcidityAffectingSpecies) = iTmp
          indSO2 = iTmp
        case('H2SO4')
          nSulphurSpecies = nSulphurSpecies + 1
          nStrongAcidSpecies = nStrongAcidSpecies + 1
          nAcidityAffectingSpecies = nAcidityAffectingSpecies + 1
          iarTmp(nSulphurSpecies) = iTmp
          iarTmpSa(nStrongAcidSpecies) = iTmp
          iarTmpAAS(nAcidityAffectingSpecies) = iTmp
        case('SO4','SOX','NH415SO4')
          nSulphurSpecies = nSulphurSpecies + 1
          nAcidityAffectingSpecies = nAcidityAffectingSpecies + 1
          iarTmp(nSulphurSpecies) = iTmp
          iarTmpAAS(nAcidityAffectingSpecies) = iTmp
        case default
      end select
    end do
    !
    ! Store the locations of the sulphur species. Needed for old standard scavenging with saturation
    !
    if(nSulphurSpecies > 0)then
      allocate(indSulphurSpecies(nSulphurSpecies), stat=iTmp)
      if(fu_fails(iTmp == 0, 'Failed to allocate sulphur species array','init_standard_deposition'))return
      indSulphurSpecies(1:nSulphurSpecies) = iarTmp(1:nSulphurSpecies)
    else
      nullify(indSulphurSpecies)
    endif
    !
    ! Store the locations of the strong-acid species
    !
    if(nStrongAcidSpecies > 0)then
      allocate(indStrongAcidSpecies(nStrongAcidSpecies), stat=iTmp)
      if(fu_fails(iTmp == 0,'Failed to allocate strong acid species array','init_standard_deposition'))return
      indStrongAcidSpecies(1:nStrongAcidSpecies) = iarTmpSa(1:nStrongAcidSpecies)
    else
      nullify(indStrongAcidSpecies)
    endif
    !
    ! Store the locations of the acidity-affecting species
    !
    if(nAcidityAffectingSpecies > 0)then
      allocate(indAcidityAffectingSp(nAcidityAffectingSpecies), stat=iTmp)
      if(fu_fails(iTmp == 0,'Failed to allocate acidity-affecting species array','init_standard_deposition'))return
      indAcidityAffectingSp(1:nAcidityAffectingSpecies) = iarTmpAAS(1:nAcidityAffectingSpecies)
    else
      nullify(indAcidityAffectingSp)
    endif

    rulesDeposition%pSpecies => speciesTransport
    call free_work_array(iarTmp)
    call free_work_array(iarTmpSa)

    !
    ! Now some checks: species must have needed parameters set. 
    ! No need to repeat them everytime...
    !
    !
    select case (rulesDeposition%scavengingType)
      case(scavStandard)
                !Gaseous species must have scavenging coeff scaling
      case(scavNoScav)
                ! Nothing is needded
      case(scav2011, scav2011fc, scav2018, scav2018entr, scav2020)
        do iSpecies = 1,rulesDeposition%nGasesDepositing
          iSpTr = rulesDeposition%indGasesDepositing(iSpecies)
          ! call msg ("Gas Depositing " +  fu_name(species(iSpTr)%material))
          DepData = fu_gas_deposition_param(speciesTransport(iSpTr)%material)
          if(fu_true(DepData%ifGasDeposited))then
            !Negative Henry constant now omits scavenging. So, NO need for the following:
            !if    (ptrDepData%Henry_const_298K < 0.)    then
            !  ptrDepData%Henry_const_298K = 0.
            !  call msg ("Henry constant is missing for:" + &
            !                 & fu_name(speciesTransport(iSpTr)%material))
            !  call msg_warning("Henry constant set to zero", "init_standard_deposition")
            !endif
            if (DepData%Henry_const_T_dep == real_missing) then
              call msg ("Assuming zero Henry_const_T_dep for:" + &
                               & fu_name(speciesTransport(iSpTr)%material))
              DepData%Henry_const_T_dep = 0.
            endif
          endif  ! ifGasDeposited
       enddo
      case default
       call msg('Unknown scavenging type',rulesDeposition%scavengingType)
       call set_error('Unknown scavenging type','init_standard_deposition')
    end select ! type of scavenging


    select case (rulesDeposition%RsType)
      case(DryD_Rs_standard)
              ! Default Rs
                !Gaseous species must have scavenging coeff scaling
      case(DryD_Rs_2013)
       do iSpecies = 1,rulesDeposition%nGasesDepositing
              iSpTr = rulesDeposition%indGasesDepositing(iSpecies)
              ! call msg ("Gas Depositing " +  fu_name(species(iSpTr)%material))
          DepData = fu_gas_deposition_param(speciesTransport(iSpTr)%material)
          if(fu_true(DepData%ifGasDeposited))then
              !Negative Henry constant now omits dry deposition. So, NO need for the following:
              !if    (ptrDepData%Henry_const_298K < 0.)    then
              !  ptrDepData%Henry_const_298K = 0.
              !  call msg ("Henry constant is missing for:" + &
              !                 & fu_name(speciesTransport(iSpTr)%material))
              !  call msg_warning("Henry constant set to zero", "init_standard_deposition")
              !endif
              if  (DepData%Henry_const_T_dep == real_missing) then
                call msg ("Assuming zero Henry_const_T_dep for:" + &
                               & fu_name(speciesTransport(iSpTr)%material))
                      DepData%Henry_const_T_dep = 0.
              endif
            if(DepData%Wesely_f0 < 0) then
              call msg ("Assuming zero Wesely_f0 for:" + fu_name(speciesTransport(iSpTr)%material))
              DepData%Wesely_f0 = 0.
            endif
          endif  ! ifGasDeposited
        enddo  ! iSpecies
      case default
       call msg('Unknown Rs method',rulesDeposition%scavengingType)
       call set_error('Unknown Rs method','init_standard_deposition')

    end select ! type of Rs


  end subroutine init_standard_deposition


  !**********************************************************************************

  subroutine add_deposition_input_needs(q_met_dyn, q_met_st, q_disp_dyn, q_disp_st, &
                                      & rulesDeposition, wdr)
    !
    ! Adds the list of quantities needed for the deposition computations
    !
    implicit none

    ! Imported parameters
    integer, dimension(:), intent(inout) :: q_met_dyn, q_met_st, q_disp_dyn, q_disp_st
    type(Tdeposition_rules), intent(in) :: rulesDeposition
    type(silja_wdr), intent(in) :: wdr

    ! Local variables
    integer :: iTmp

    !
    ! The information needed depends on the bulkiness and complexity of the procedures involved.
    ! They, in turn, are decided by the list of species to be treated and rules
    !
    ! Aerosol settling: so far here. Later, possibly, will be moved to special module.
    !
    iTmp = fu_merge_integer_to_array(temperature_flag, q_met_dyn)
    iTmp = fu_merge_integer_to_array(pressure_flag, q_met_dyn)

    if(rulesDeposition%ifHumidityDependent)then
      iTmp = fu_merge_integer_to_array(relative_humidity_flag, q_met_dyn)
      iTmp = fu_merge_integer_to_array(relative_humidity_2m_flag, q_met_dyn)
    endif

    if (rulesDeposition%DryDepType /= DryD_Vd_Zero) then
      ! the schemes that use vd are below
      ! they do not distinguish between gas and particle deposition
      ! they also have same requirements
      !
      iTmp = fu_merge_integer_to_array(relative_humidity_flag, q_met_dyn)
      iTmp = fu_merge_integer_to_array(friction_velocity_flag, q_met_dyn)
      iTmp = fu_merge_integer_to_array(temperature_2m_flag, q_met_dyn)
      iTmp = fu_merge_integer_to_array(surface_pressure_flag, q_met_dyn)
      iTmp = fu_merge_integer_to_array(MO_length_inv_flag, q_met_dyn)
      !iTmp = fu_merge_integer_to_array(surface_roughness_disp_flag, q_met_dyn)
      iTmp = fu_merge_integer_to_array(surface_roughness_meteo_flag, q_met_dyn)
      iTmp = fu_merge_integer_to_array(fraction_of_land_flag, q_met_st)
      iTmp = fu_merge_integer_to_array(total_precipitation_rate_flag, q_met_st)
      iTmp = fu_merge_integer_to_array(SILAM_sensible_heat_flux_flag, q_met_dyn)
!      iTmp = fu_merge_integer_to_array(Vd_correction_DMAT_flag, q_disp_dyn) ! Needed for Rs
    endif


      select case(rulesDeposition%DryDepType)
        case(DryD_KS2011)
#ifdef DEBUG
         call msg('**DEP**: Preparing KS2011 depo')
#endif
        
        case(DryD_KS2011_TF)
#ifdef DEBUG
          call msg('**DEP**: Preparing KS2011-TF depo')
#endif

        case(DryD_Vd_Zero)
#ifdef DEBUG
           call msg('**DEP**: Preparing for NO_DD depo')
#endif

        case(DryD_Vd_difsed)
#ifdef DEBUG
        call msg('**DEP**: Preparing Vd_difsed depo')
#endif

        case(DryD_Vd_SP98)
#ifdef DEBUG
        call msg('**DEP**: Preparing Vd_SP98 depo')
#endif

        case default
          call msg('Unknown dry deposition type:',rulesDeposition%DryDepType)
          call set_error('Unknown deposition type','add_deposition_input_needs')
          return
      end select ! rulesDeposition%DryDepType

!    endif 


    !
    ! Rs
    !
    if (rulesDeposition%DryDepType /= DryD_Vd_Zero) then
       iTmp = fu_merge_integer_to_array(total_precipitation_rate_flag, q_met_st)
       iTmp = fu_merge_integer_to_array(fraction_of_land_flag, q_met_st)
    endif
    select case (rulesDeposition%RsType)
      case(DryD_Rs_standard)
        ! Default Rs
        !Gaseous species must have scavenging coeff scaling
      case(DryD_Rs_2013)
        iTmp = fu_merge_integer_to_array(relative_humidity_2m_flag, q_met_dyn)
        iTmp = fu_merge_integer_to_array(water_eq_snow_depth_flag, q_met_dyn)
        iTmp = fu_merge_integer_to_array(stomatal_conductance_flag, q_met_dyn)
!        iTmp = fu_merge_integer_to_array(fraction_of_ice_flag, q_met_dyn)
        if ( wdr == wdr_missing) then  ! Can be undefined
           iTmp = fu_merge_integer_to_array(leaf_area_index_flag, q_met_dyn)
        elseif (any (fu_LAIsrc(wdr) == (/LAI_dynamic_1, LAI_dynamic_2/))) then
          iTmp = fu_merge_integer_to_array(leaf_area_index_flag, q_met_dyn)
        elseif(any (fu_LAIsrc(wdr) == (/LAI_static_1, LAI_static_2/)))then
          iTmp = fu_merge_integer_to_array(leaf_area_index_flag, q_met_st)
        else
          call set_error('Dry deposition Rs_2013 requires leaf area index','add_deposition_input_needs')
          return
        endif
      case default
        call msg('Unknown Rs method',rulesDeposition%RsType)
        call set_error('Unknown Rs method','add_deposition_input_needs')

    end select ! rulesDeposition%RsType

    !
    ! Now scavenging. Might need the ABL height as well to position the cloud bottoma

    !! Needed for Lagrangian transport only, the input needs for scavenging in Eulerian transport 
    !! are set in the subroutine wet_deposition_input_needs

    if (rulesDeposition%if_lagr) then

     select case(rulesDeposition%scavengingType)
       case(scavStandard)
         iTmp = fu_merge_integer_to_array(total_precipitation_rate_flag, q_met_st)
         iTmp = fu_merge_integer_to_array(scavenging_coefficient_flag, q_met_st)
 !        iTmp = fu_merge_integer_to_array(abl_height_m_flag, q_met_dyn)

       case(scav2011,scav2011fc)
         iTmp = fu_merge_integer_to_array(total_precipitation_rate_flag, q_met_st)
         iTmp = fu_merge_integer_to_array(total_cloud_cover_flag, q_met_dyn)
         iTmp = fu_merge_integer_to_array(pressure_flag, q_met_dyn)
         iTmp = fu_merge_integer_to_array(temperature_flag, q_met_dyn)
         iTmp = fu_merge_integer_to_array(cell_size_z_flag, q_disp_dyn)
         if (rulesDeposition%scavengingType == scav2011) &
                & iTmp = fu_merge_integer_to_array(cwcabove_3d_flag, q_met_dyn)
      
       case(scavNoScav)

       case default
         call set_error('Unsupported scavenging type','add_deposition_input_needs')
     end select ! rulesDeposition%scavengingType
  end if

  end subroutine add_deposition_input_needs

  subroutine wet_deposition_input_needs(rulesDeposition, meteo_input_local)
    implicit none
    ! Imported parameters                                                                                                        
    type (Tdeposition_rules), intent(in) :: rulesDeposition
    type(Tmeteo_input), intent(out), target :: meteo_input_local
    character(len=*), parameter :: subname="wet_deposition_input_needs"

    ! Local variables                                                                        
    integer :: iQ, iTmp, nq

    meteo_input_local =  meteo_input_empty
    if(.not. rulesDeposition%defined == silja_true)then
      call set_error('Deposition rules are not defined',subname)
      return
    endif

    if (rulesDeposition%scavengingType == scavNoScav) return
      
    nq=2
    meteo_input_local%quantity(1:2) = (/cell_size_x_flag, cell_size_y_flag/)
    meteo_input_local%q_type(1:2) = dispersion_single_time_flag
    imet_dx_size =>  meteo_input_local%idx(1)
    imet_dy_size =>  meteo_input_local%idx(2)

    nq = nq +1
    meteo_input_local%quantity(nq) = total_precipitation_rate_flag
    meteo_input_local%q_type(nq) = meteo_single_time_flag
    imet_prec => meteo_input_local%idx(nq)

    select case(rulesDeposition%scavengingType)                 
      case(scavStandard)        
        nq = nq +1
        meteo_input_local%quantity(nq) = scavenging_coefficient_flag
        meteo_input_local%q_type(nq) = meteo_single_time_flag
        imet_scav => meteo_input_local%idx(nq)

        nq = nq +1
        meteo_input_local%quantity(nq) = temperature_flag
        meteo_input_local%q_type(nq) = meteo_dynamic_flag
        imet_temp => meteo_input_local%idx(nq)

      case(scav2011,scav2011fc)
        nq = nq +1
        meteo_input_local%quantity(nq) = total_cloud_cover_flag
        meteo_input_local%q_type(nq) = meteo_dynamic_flag
        imet_tcc => meteo_input_local%idx(nq)

        nq = nq +1
        meteo_input_local%quantity(nq) = pressure_flag
        meteo_input_local%q_type(nq) = meteo_dynamic_flag
        imet_press => meteo_input_local%idx(nq)

        nq = nq +1
        meteo_input_local%quantity(nq) = temperature_flag
        meteo_input_local%q_type(nq) = meteo_dynamic_flag
        imet_temp => meteo_input_local%idx(nq)

        nq = nq +1
        meteo_input_local%quantity(nq) = cell_size_z_flag
        meteo_input_local%q_type(nq) = dispersion_dynamic_flag
        imet_dz_size => meteo_input_local%idx(nq)

        if (rulesDeposition%scavengingType == scav2011) then
          nq = nq +1
          meteo_input_local%quantity(nq) = cwcabove_3d_flag
          meteo_input_local%q_type(nq) = meteo_dynamic_flag
          imet_cwc3d => meteo_input_local%idx(nq)
        end if

      case(scav2018, scav2018entr, scav2020)
        nq = nq +1
        meteo_input_local%quantity(nq) = total_cloud_cover_flag
        meteo_input_local%q_type(nq) = meteo_dynamic_flag
        imet_tcc => meteo_input_local%idx(nq)

        nq = nq +1
        meteo_input_local%quantity(nq) = pressure_flag
        meteo_input_local%q_type(nq) = meteo_dynamic_flag
        imet_press => meteo_input_local%idx(nq)

        nq = nq +1
        meteo_input_local%quantity(nq) = temperature_flag
        meteo_input_local%q_type(nq) = meteo_dynamic_flag
        imet_temp => meteo_input_local%idx(nq)

        nq = nq +1
        meteo_input_local%quantity(nq) = cell_size_z_flag
        meteo_input_local%q_type(nq) = dispersion_dynamic_flag
        imet_dz_size => meteo_input_local%idx(nq)

        nq = nq +1
        meteo_input_local%quantity(nq) = cwcabove_3d_flag
        meteo_input_local%q_type(nq) = meteo_dynamic_flag
        imet_cwc3d => meteo_input_local%idx(nq)
        
        if(rulesDeposition%scavengingType == scav2018 .or. &
             & rulesDeposition%scavengingType == scav2018entr) then
          nq = nq +1
          meteo_input_local%quantity(nq) = pwcabove_3d_flag
          meteo_input_local%q_type(nq) = meteo_dynamic_flag
          imet_pwc3d => meteo_input_local%idx(nq)
        end if

        if(rulesDeposition%scavengingType == scav2020) then
           nq = nq +1
           meteo_input_local%quantity(nq) = cloud_water_flag
           meteo_input_local%q_type(nq) = meteo_dynamic_flag
           imet_cwc => meteo_input_local%idx(nq)

           nq = nq +1
           meteo_input_local%quantity(nq) = cloud_ice_flag
           meteo_input_local%q_type(nq) = meteo_dynamic_flag
           imet_cic => meteo_input_local%idx(nq)

           nq = nq +1
           meteo_input_local%quantity(nq) = cloud_cover_flag
           meteo_input_local%q_type(nq) = meteo_dynamic_flag
           imet_cc3d => meteo_input_local%idx(nq)

           nq = nq +1
           meteo_input_local%quantity(nq) = air_density_flag
           meteo_input_local%q_type(nq) = dispersion_dynamic_flag
           imet_airdens => meteo_input_local%idx(nq)

           nq = nq +1
           meteo_input_local%quantity(nq) = disp_cell_airmass_flag
           meteo_input_local%q_type(nq) = dispersion_dynamic_flag
           imet_airmass => meteo_input_local%idx(nq)

        end if

        if(rulesDeposition%max_scav_rate_depends_on == horiz_wind .or. &
             & rulesDeposition%max_scav_rate_depends_on == CAPE_AND_WIND) then
          nq = nq +1
          meteo_input_local%quantity(nq) = u_flag
          meteo_input_local%q_type(nq) = meteo_dynamic_flag
          imet_u => meteo_input_local%idx(nq)

          nq = nq +1
          meteo_input_local%quantity(nq) = v_flag
          meteo_input_local%q_type(nq) = meteo_dynamic_flag
          imet_v => meteo_input_local%idx(nq)
       else if (rulesDeposition%max_scav_rate_depends_on == CAPE) then
          nq = nq +1
          meteo_input_local%quantity(nq) = cape_flag
          meteo_input_local%q_type(nq) = meteo_dynamic_flag
          imet_cape => meteo_input_local%idx(nq)

          nq = nq +1
          meteo_input_local%quantity(nq) = abl_height_m_flag
          meteo_input_local%q_type(nq) = meteo_dynamic_flag
          imet_abl => meteo_input_local%idx(nq)
          
        end if
        if (rulesDeposition%max_scav_rate_depends_on == CAPE_AND_WIND) then
          nq = nq +1
          meteo_input_local%quantity(nq) = cape_flag
          meteo_input_local%q_type(nq) = meteo_dynamic_flag
          imet_cape => meteo_input_local%idx(nq)

          nq = nq +1
          meteo_input_local%quantity(nq) = abl_height_m_flag
          meteo_input_local%q_type(nq) = meteo_dynamic_flag
          imet_abl => meteo_input_local%idx(nq)

          nq = nq +1
          meteo_input_local%quantity(nq) = fraction_of_land_flag
          meteo_input_local%q_type(nq) = meteo_single_time_flag
          imet_landfrac => meteo_input_local%idx(nq)
        end if

      case(scavNoScav)

      case default                                                   
        call set_error('Unknown scavenging type', subname)                

    end select ! rulesDeposition%scavengingType

    meteo_input_local%nQuantities = nq

  end subroutine wet_deposition_input_needs


  !**********************************************************************************

  subroutine prepare_deposition(met_buf, disp_buf, rulesDeposition)
    !
    ! The subroutine prepares the private module pointers to the fields requested in the 
    ! deposition input needs
    !
    implicit none

    ! Imported parameters
    type(Tfield_buffer), pointer :: met_buf, disp_buf
    type(Tdeposition_rules), intent(in) :: rulesDeposition

    ! Local variables
    integer, dimension(max_quantities) :: q_met_dyn, q_met_st, q_disp_dyn, q_disp_st
    integer, dimension(:), pointer :: mdl_in_q
    integer :: iTmp, iQ

    !
    ! Start from the collection of what actually is needed
    !
    q_met_dyn(:) = int_missing
    q_met_st(:) = int_missing
    q_disp_dyn(:) = int_missing
    q_disp_st(:) = int_missing

    ! We do not really care here if patricular quantity is singletime or dynamic,
    ! thus wdr_missing is a good option....
    call add_deposition_input_needs(q_met_dyn, q_met_st, q_disp_dyn, q_disp_st, &
                                  & rulesDeposition, wdr_missing)
    
    if(error)return
    !
    ! Scan the meteo buffer
    !
    mdl_in_q => met_buf%buffer_quantities

    do iQ = 1, size(q_met_dyn)
      if(q_met_dyn(iQ) == int_missing)exit
      iTmp = fu_index(mdl_in_q, q_met_dyn(iQ))
      if(iTmp < 0)then
        call set_error('No dynamic meteo quantity:'+fu_quantity_string(q_met_dyn(iQ)), &
                     & 'prepare_deposition')
        return
      endif
      call set_dep_pointer(met_buf, iTmp, q_met_dyn(iQ))
    end do
    do iQ = 1, size(q_met_st)
      if(q_met_st(iQ) == int_missing)exit
      iTmp = fu_index(mdl_in_q, q_met_st(iQ))
      if(iTmp < 0)then
        call set_error('No static meteo quantity:'+fu_quantity_string(q_met_st(iQ)), &
                     & 'prepare_deposition')
        return
      endif
      call set_dep_pointer(met_buf, iTmp, q_met_st(iQ))
    end do
    !
    ! Scan the dispersion buffer
    !
    mdl_in_q => disp_buf%buffer_quantities

    do iQ = 1, size(q_disp_dyn)
      if(q_disp_dyn(iQ) == int_missing)exit
      iTmp = fu_index(mdl_in_q, q_disp_dyn(iQ))
      if(iTmp < 0)then
        call set_error('No dynamic dispersion quantity:'+fu_quantity_string(q_disp_dyn(iQ)), &
                     & 'prepare_deposition')
        return
      endif
      call set_dep_pointer(disp_buf, iTmp, q_disp_dyn(iQ))
    end do
    do iQ = 1, size(q_disp_st)
      if(q_disp_st(iQ) == int_missing)exit
      iTmp = fu_index(mdl_in_q, q_disp_st(iQ))
      if(iTmp < 0)then
        call set_error('No static dispersion quantity:'+fu_quantity_string(q_disp_st(iQ)), &
                     & 'prepare_deposition')
        return
      endif
      call set_dep_pointer(disp_buf, iTmp, q_disp_st(iQ))
    end do

    meteo_cell_size = sqrt(fu_cell_size(meteo_grid)) !!!! horizontal scale of meteo grid


    CONTAINS

    !==================================================================================

    subroutine set_dep_pointer(buf, indexQ, quantity)
      !
      ! Actually sets the requested pointer. First it checks whether the given index is 2D or 4D, then
      ! decides what to set
      !
      implicit none

      ! Imported parameters
      type(Tfield_buffer), pointer :: buf
      integer, intent(in) :: indexQ, quantity

      ! Local variables
      integer, pointer :: iP
      !
      ! 2D or 4D ?
      !
      if(fu_dimension(buf, indexQ) == 4)then
        !
        ! Choose the pointer to 4D field
        !
        select case(quantity)
          case(scavenging_coefficient_flag)
            fldScavStd => buf%p4d(indexQ)
! call msg('**DEP**: scavenging_coefficient_flag')

         case(u_flag)
           fldUmet => buf%p4d(indexQ)
! call msg('**DEP**: u_flag')

         case(v_flag)
           fldVmet => buf%p4d(indexQ)
! call msg('**DEP**: v_flag')

         case(temperature_flag)
           fldTempr => buf%p4d(indexQ)
! call msg('**DEP**: temperature_flag')

         case(relative_humidity_flag)
           fldRelHumidity => buf%p4d(indexQ)
! call msg('**DEP**: relative_humidity_flag')

         case(pressure_flag)
           fldPressure => buf%p4d(indexQ)
! call msg('**DEP**: pressure_flag')

         case(cell_size_z_flag)
           fldCellSizeZ => buf%p4d(indexQ)
! call msg('**DEP**: cell_size_z_flag')

         case(cwcabove_3d_flag)
           fldCwcabove => buf%p4d(indexQ)
! call msg('**DEP**: cwcabove_3d_flag')

         case(pwcabove_3d_flag)
           fldPwcabove => buf%p4d(indexQ)
! call msg('**DEP**: pwcabove_3d_flag')

         case default
           write (*,*) "Unknown 4D quantity dimension: ", (fu_dimension(buf, indexQ))
           call set_error('Unknown quantity:'+fu_quantity_string(quantity),'set_dep_pointer')
           return
       end select
     else
       !
       ! Choose the pointer to 2D array of reals
       !
       select case(quantity)
         case(fraction_of_land_flag)
           pMetLandFr => buf%p2d(indexQ)%present%ptr
! call msg('**DEP**: fraction_of_land_flag')

         case(fraction_of_ice_flag)
           pMetIceFr => buf%p2d(indexQ)%present%ptr
! call msg('**DEP**: fraction_of_ice_flag')

         case(temperature_2m_flag)
           pmettempr2m => buf%p2d(indexq)%present%ptr
! call msg('**DEP**: temperature_2m_flag')

         case(relative_humidity_2m_flag)
           pMetRh2m => buf%p2d(indexq)%present%ptr
! call msg('**DEP**: temperature_2m_flag')

         case(friction_velocity_flag)
           pMetFricVel => buf%p2d(indexQ)%present%ptr
! call msg('**DEP**: friction_velocity_flag')

         case(total_precipitation_rate_flag)
           pMetPrecTot => buf%p2d(indexQ)%present%ptr
! call msg('**DEP**: total_precipitation_rate_flag')

         case(large_scale_rain_int_flag)
           pMetPrecLs => buf%p2d(indexQ)%present%ptr
! call msg('**DEP**: total_precipitation_rate_flag')

         case(convective_rain_int_flag)
           pMetPrecCnv => buf%p2d(indexQ)%present%ptr
! call msg('**DEP**: total_precipitation_rate_flag')

         case(abl_height_m_flag)
           pMetABLHeight => buf%p2d(indexQ)%present%ptr
! call msg('**DEP**: abl_height_m_flag')

         case(MO_length_inv_flag)
           pMetMO_Len_inv => buf%p2d(indexQ)%present%ptr
! call msg('**DEP**: MO_length_inv_flag')

         case(surface_pressure_flag)
           pMetSrfPressure => buf%p2d(indexQ)%present%ptr
! call msg('**DEP**: surface_pressure_flag')

         case(surface_roughness_meteo_flag)
           pMetSrfRoughMeteo => buf%p2d(indexQ)%present%ptr
! call msg('**DEP**: surface_roughness_flag')
         
         case(surface_roughness_disp_flag)
           pMetSrfRoughDisp => buf%p2d(indexQ)%present%ptr
! call msg('**DEP**: surface_roughness_flag')
         
       case(convective_velocity_scale_flag)
           pMetConvVel => buf%p2d(indexQ)%present%ptr
! call msg('**DEP**: convective_velocity_scale_flag')

       case(cape_flag)
           pCAPE => buf%p2d(indexQ)%present%ptr

       case(total_cloud_cover_flag)
           pMetTotCloud => buf%p2d(indexQ)%present%ptr
! call msg('**DEP**: total_cloud_cover_flag')

       case(SILAM_sensible_heat_flux_flag)
           pMetSensHF => buf%p2d(indexQ)%present%ptr
! call msg('**DEP**: SILAM_sensible_heat_flux_flag')

!          case(Vd_correction_DMAT_flag)
!            pDispVdCorrection => buf%p2d(indexQ)%present%ptr
!call msg('**DEP**: Vd_correction_DMAT_flag')

        case(leaf_area_index_flag)
                pMetLAI => buf%p2d(indexQ)%present%ptr

        case(stomatal_conductance_flag)
                pMetGsto => buf%p2d(indexQ)%present%ptr

        case(water_eq_snow_depth_flag)
                pMetSnowDepth => buf%p2d(indexQ)%present%ptr

          case default
            call set_error('Unknown 2D quantity:'+fu_quantity_string(quantity),'set_dep_pointer')
            return
        end select
      endif  ! 2D or 4D

    end subroutine set_dep_pointer

  end subroutine prepare_deposition


  !********************************************************************************************
  !
  !    Dry deposition
  !
  !********************************************************************************************

  real function  fu_get_vd(zref,  speciesTransport, &  ! Reference height, what to deposit
                         & indexMeteo, weight_past, &  ! position in space and time
                         & Rs, deptype, timeSign, invVd2m)     ! surface resistance
    !
    ! Returns the deposition velocity for given species and given height
    ! with new KS2011 deposition scheme plus thermpohoresis.
    ! For gases the return value is an inverse resistance
    !
    ! ATTENTION. Uses the precomputed arrays and parameters, so the species must be the 
    !            very transport species, which were used for the precomputation step. This is 
    !            checked upon the first call
    ! Units: SI
    ! Author: Roux
    !
    implicit none

    ! Imported parameters
    real, intent(in) :: zref                         ! at what height?
    type(silam_species), intent(in) :: speciesTransport !what depositing
    integer, intent(in) :: indexMeteo                   ! where
    real, intent(in) :: weight_past                     ! when
    real, intent(in) :: Rs                              ! surface resistance
    integer, intent(in) :: deptype                      ! deposition type
    integer, intent(in) :: timeSign                     ! 1 for forward, -1 for backward
    real, intent(inout) :: invVd2m                            ! resistance from zref to 2m 
    
    ! Local variables
    real :: fLambda, fCun, fWetDiam, fSc, fRoughfr, fStickRatio, fKin_visc, fDyn_visc, fDiffusivity, &
          & fSetlVelocity, fTau, fHumidity, invL, u_star, mol_diff, fZ0, fTemperature, &
          & fVTFplus, fSensHeatVelplus, vsplusCorrTF, Cp_ro, fTmp, R2m, V, fIvd
    real, parameter :: ascale = 2e-3  ! Collection scale is prescribed for
                                          ! land. Should be taken from landuse
    real, parameter :: fPrm = 0.7  ! Molecular prandtl number
    real, parameter :: fkgkp = 0.02  ! Ratio of heat conductivities
                                         ! for gas and particle
    real, parameter :: Z2m = 2.  ! Height fo 2m

    real :: tauplus, vsplus, rplus, Rsplus, invLplus, zplus,  Re, Sc, St, dpa, fKn !  dimensionless
    logical :: fwdOnly ! Can calculate vd only for forward time
                                                       !  parameters
    real :: vdsmooth, vdrough       ! vd for two landuse types. Should be more
    type(TwetParticle) :: wetParticle
    integer, save :: iCount=0

    vdsmooth = 0.
    vdrough = 0.
    ! Meteo
    fDyn_visc = fu_dynamic_viscosity(pMetTempr2m(indexMeteo))
    fKin_visc = fDyn_visc * pMetSrfPressure(indexMeteo) / &
              & (gas_constant_dryair * pMetTempr2m(indexMeteo))
    u_star =   pMetFricVel(indexMeteo)

    invL =     pMetMO_Len_inv(indexMeteo)
    fRoughfr = pMetLandFr(indexMeteo)
    !fZ0      =  pMetSrfRoughDisp(indexMeteo)
    fZ0      =  pMetSrfRoughMeteo(indexMeteo)
        
    Cp_ro = specific_heat_dryair*pMetSrfPressure(indexMeteo) / &
          & (gas_constant_dryair*pMetTempr2m(indexMeteo))
    fSensHeatVelplus = pMetSensHF(indexMeteo)/(Cp_ro * pMetTempr2m(indexMeteo) * u_star )
    !        - u_star*u_star*invL/(0.4*g) 
    ! <w'T'>/T/ u_*


    if(fZ0 < 1e-10)then
      if(mod(iCount, 1000) == 1) call msg('Strange z0 many times:', fZ0, iCount)
      iCount = iCount + 1
      fz0=1.e-6
    endif

    if (fZ0 < 3. * fKin_visc/u_star) then !FIXME: May be some other coefficient needed here
      fRoughfr = 0.
    endif
                                               
    if(pMetPrecTot(indexMeteo) > 0.0) then ! Dimensionless meteoswitch for impaction 
      fStickRatio = 1.0 
    else
      fStickRatio = 0.0
    endif 
        
    !
    ! Dimensionless species properties
    if(fu_mode(speciesTransport) == in_gas_phase)then    !gas
      ! Sc = nu / D (kinematic viscosity over molecular diffusivity)
      fSc = fKin_visc/fu_gas_molecular_diffusivity_air(speciesTransport%material)
      fWetDiam = real_missing 
      tauplus = 0.
      vsplus = 0.
      rplus = 0.
      Rsplus = Rs*u_star
      fTau = 0.
      fVTFplus = 0.
    else    
      ! particles
      ! No deposition rules here....
      !if(rulesDeposition%ifHumidityDependent)then
      ! FIXME: Surface-layer humidity should be here!
      fHumidity = fldRelHumidity%past%p2d(1)%ptr(indexMeteo) * weight_past + &
           & fldRelHumidity%future%p2d(1)%ptr(indexMeteo) * (1.-weight_past)
      fTemperature = fldTempr%past%p2d(1)%ptr(indexMeteo) * weight_past + &
           & fldTempr%future%p2d(1)%ptr(indexMeteo) * (1.-weight_past)
      !else
      !  fHumidity = rulesDeposition%fDefaultRelHumidity
      !endif
      wetParticle = fu_wet_particle_features( speciesTransport%material, fHumidity)
      fWetDiam = wetParticle%fGrowthFactor * fu_massmean_D( speciesTransport%mode )
      if(fWetDiam > 25e-6 .and. wetParticle%fGrowthFactor > 5)then
        fWetDiam = 25e-6
      endif
      fLambda = 2.37e-5 * pMetTempr2m(indexMeteo) / pMetSrfPressure(indexMeteo)
      fKn = 2*fLambda / fWetDiam
      fCun = 1. + fKn * (1.257 + 0.4 * exp(-1.1*fKn))

      fDiffusivity = boltzmann_const * pMetTempr2m(indexMeteo) * fCun / &
                   & (3. * Pi * fDyn_visc * fWetDiam)
      fTau =  fWetDiam * fWetDiam * wetParticle%fWetParticleDensity * &
                       & fCun / (18. * fDyn_visc)
      fVTFplus = - fSensHeatVelplus * fPrm  & !Thermophoretic velocity
          &        * 2*fCun *1.17 *(fkgkp+2.18*fKn) & ! (Friedlander 2000)
          &        / (1+3*1.14*fKn)/(1.+2.*fkgkp +2.*2.18*fKn) * timeSign
      fSc = fKin_visc / fDiffusivity
      tauplus = fTau * u_star * u_star / ( fKin_visc)
      vsplus = fTau * g / u_star * timeSign
      rplus  = 0.5 * fWetDiam * u_star / fKin_visc
      Rsplus = 0. ! No surface resistance for particles
    endif ! gas or particles
        
!! Temporary hack
!fu_get_vd = & !1e-8 +  fTau * g
!  & u_star * fu_vdplus_DS(tauplus,fSc,vsplus, Rsplus, invL*fZ0, zref/fZ0)
!return
    SELECT CASE(deptype)
      CASE(DryD_KS2011, DryD_KS2011_TF) ! Original KS scheme
        vsplusCorrTF=0. ! correction of VsPlus for thermophoresis
        if (deptype == DryD_KS2011_TF)  vsplusCorrTF = fvtfplus !  thermophoresis
        if (fRoughfr < 1.) then ! smooth fraction exist
          fTmp = fKin_visc/u_star ! z0smooth
          vdsmooth = (1. - fRoughfr) * u_star * fu_vdplus_smooth(tauplus, fSc, vsplus + vsplusCorrTF, &
                                                               & rplus, Rsplus, invL*fTmp , zref/fTmp)
          if (.not. (vdsmooth .ge. 0.)) then
            if (zref < fTmp) then
              ! Could be just too-low center of masses
              vdsmooth =  (1. - fRoughfr) * u_star ! whatever
            else
              call msg ("**DEPO** -fu_get_vd-------------------------")
              call msg ("**DEPO** Negative vdsmooth!, WetDp=",fWetDiam)
              call msg ("**DEPO** Negative vdsmooth!, tauplus=",tauplus)
              call msg ("**DEPO** Negative vdsmooth!, fsc=",fsc)
              call msg ("**DEPO** Negative vdsmooth!, vsplus=",vsplus+vsplusCorrTF)
              call msg ("**DEPO** Negative vdsmooth!, rplus=",rplus)
              call msg ("**DEPO** Negative vdsmooth!, Rsplus=",Rsplus)
              call msg ("**DEPO** Negative vdsmooth!, Rs=",Rs)
              call msg ("**DEPO** Negative vdsmooth!, vdsmooth=",vdsmooth)
              call msg ("**DEPO** Negative vdsmooth!, zref=",zref)
              call msg ("**DEPO** Negative vdsmooth!, z0smooth=",fTmp)
              call msg ("**DEPO** Negative vdsmooth!, invL=",invL)
              vdsmooth = 0. ! Should not happen, but happens. 
                            ! sometimes returns NAN due to over-/under- flow
            endif
          endif 
       endif
        vsplusCorrTF = 0.5*(vsplusCorrTF+abs(vsplusCorrTF)) ! Can not be negative for rough surfaces. 
                                                            ! Could be treated in some better way
        if (fRoughfr > 0.) then ! rough fraction exist
           Re =  ascale * u_star / fKin_visc
           St = 2*tauplus/Re * fStickRatio ! No impaction for non-sticky
                                          ! surfaces
           dpa = 2.*rplus/Re
           vdrough = fRoughfr * u_star * &
                fu_vdplus_rough(St, fSc, vsplus + vsplusCorrTF, dpa, Rsplus, Re, invL*fZ0 , max(zref/fZ0, 1.))
            if (.not. (vdrough .ge. 0.)) then
                  call msg("**DEPO** Bad vdrough in fu_get_vd! vdrough,fz0:", vdrough, fZ0)
                  call msg('fRoughfr, u_star:',fRoughfr, u_star)
                  call msg('St, fSc',St, fSc)
                  call msg('vsplus, vsplusCorrTF',vsplus, vsplusCorrTF)
                  call msg('dpa, Rsplus',dpa, Rsplus)
                  call msg('Re, invL',Re, invL)
                  call msg('fZ0 , zref',fZ0 , zref)
                  vdrough = 0. ! FIXME Should never happen
            endif
        endif          
        fu_get_vd = vdsmooth + vdrough 
        fwdOnly = .false.
        


      CASE (DryD_Vd_difsed) ! Vd ==  diffusion + settling
         fu_get_vd = u_star * &
               & fu_vdplus_DS(tauplus,fSc,abs(vsplus), Rsplus, invL*fZ0, zref/fZ0)
        fwdOnly = .true.

      CASE (DryD_Vd_SP98) ! Vd Slinn (Accoding to Seinfeld Pandis 2006)
         fu_get_vd = u_star * & 
               & fu_vdplus_slinnSP98(tauplus,fSc,abs(vsplus), Rsplus, invL*fZ0, zref/fZ0)
        fwdOnly = .true.

      CASE (DryD_VD_Zhang)! Generic Zhang... Grass for all
        fu_get_vd = fu_vd_Zhang(fTau,fSc,fWetDiam,u_star, Zref, Rs, invL, ZhangGrass)
        fwdOnly = .true.

      CASE (DryD_VD_Zero)! Generic Zhang... Grass for all
        fu_get_vd = 0
        fwdOnly = .false.

      CASE DEFAULT ! Old scheme
            call msg("Unknown dry depo  rules! rules:", deptype) 
            call set_error("Unknown dry depo  rules", "fu_get_vd") 
    END SELECT

    ! correct Vd for inverse time
    if (fwdOnly .and. (timeSign < 0)) then
      fu_get_vd = fu_get_vd - u_star*abs(vsplus)
      fu_get_vd = 0.5*(fu_get_vd + abs(fu_get_vd)) + 1e-10 ! 1/Vd is still meaningful
    endif

    if (invVd2m == real_missing) return  !! No need in invVd2m

    !! Remaining stuff -- Vd2m

    if (fu_get_vd == 0.) then
      invVd2m = real_missing
    elseif (Z2m > fZ0) then  ! reasonable profile to 2m is possible
      !
      ! R2m > 0 for zref > Z2m and < 0  zref < Z2m. Can be wrong sign in strong stratifications
      !
      R2m = 2.5 * (log(Zref/Z2m) + fu_Psi(Zref*invL) - fu_Psi(Z2m*invL))
      V = g * fTau * timeSign
      invVd2m = 1./ fu_get_vd
      if (zref>Z2m) then 
        R2m = max(0., R2m)   ! Could become negative in unstable
        fTmp = V*R2m
        if (-fTmp < LOG_MAX_REAL ) then  !exp defined
          if (abs(fTmp) .gt. 0.001) then
            fTmp = exp(-fTmp)
            invVd2m = invVd2m*fTmp + (1.-fTmp)/V
          else
            invVd2m = invVd2m - R2m
          endif
        endif
      else
        R2m = -min(R2m, 0.0)   ! Could become wrong sign. In the end, should become positive
        fTmp = V*R2m
        if (fTmp < 0.1*LOG_MAX_REAL ) then  !exp defined Plus hack.. 
          if (abs(fTmp) .gt. 0.001) then
            fTmp = exp(fTmp)
            invVd2m = invVd2m*fTmp + (fTmp - 1.)/V !!Hack needed for this line
          else
            invVd2m = invVd2m + R2m
          endif
        endif
      endif
      if (.not.(invVd2m > 0)) then
        if (.not. fRoughfr > 0) then !! Initialize them to something
          Re = -1
          St = -1
          dpa = -1
        endif
        call msg_warning("Negative Vd2m", "fu_get_vd")
        call msg("Resistance invVd2m", invVd2m)
        call msg("Resistance zRef-2m", R2m)
        call msg("v_d", fu_get_vd)
        call msg('fRoughfr, u_star:',fRoughfr, u_star)
        call msg('St, fSc',St, fSc)
        call msg('fTau', fTau)
        call msg('vsplus, vsplusCorrTF',vsplus, vsplusCorrTF)
        call msg('dpa, Rsplus',dpa, Rsplus)
        call msg('Re, invL',Re, invL)
        call msg('fZ0 , zref',fZ0 , zref)
        call msg("vdsmooth, vdrough", vdsmooth, vdrough)
        invVd2m = fu_get_vd !Safe choice
  
!                  !FIXME This should not happen, but happened in apta run....
!   Negative Vd2m     -22.5560760
!   fRoughfr, u_star:       0.8557873       1.4377022
!   St, fSc   0.0000000E+00   0.2825128E+01
!   vsplus, vsplusCorrTF   0.0000000E+00   0.0000000E+00
!   dpa, Rsplus   0.0000000E+00   0.3594255E+02
!   Re, invL     127.2245331       0.2871192
!   fZ0 , zref       0.0015802      71.2181854
!   Strange Cnc2m:  pCnc2m%arM(iSpecies,iSrc,1,ix,iy):      -0.0442355
!   Resulting concentration (/m3)  -0.5795835E-11
!   weightUp  -0.9999990E+15
!   dh1, dh2      79.4501953     210.7581329
!   zCM(1),zCM(2)       0.3963878      -0.0166495
!   M(1), M(2)      34.2621117    3755.0971680
!   vd, vd2m       0.0045477      -0.0443340
!   Rs       25.0000000
      endif
    else
       !Inside roughness 
       invVd2m = 1./ fu_get_vd
    endif

  end function  fu_get_vd

  !****************************************************************************************
  
  subroutine  get_vd_species(zref,  speciesTransport, &  ! Reference height, what to deposit
                         & indexMeteo, weight_past, &  ! position in space and time
                         & Rs, deptype, timeSign, fVd, invVd2m)     ! surface resistance
    !
    ! Returns the deposition velocity for given species and given height
    ! with new KS2011 deposition scheme plus thermpohoresis.
    ! For gases the return value is an inverse resistance
    !
    ! ATTENTION. Uses the precomputed arrays and parameters, so the species must be the 
    !            very transport species, which were used for the precomputation step. This is 
    !            checked upon the first call
    ! Units: SI
    ! Author: Roux
    !
    implicit none

    ! Imported parameters
    real, dimension(:,:),  intent(in) :: zref    ! (nSp, nSrc)  at what height?
    type(silam_species), dimension(:), intent(in) :: speciesTransport ! (nSp) what depositing
    integer, intent(in) :: indexMeteo            ! where
    real, intent(in) :: weight_past              ! when
    real, dimension(:), intent(in) :: Rs         ! (nSp) surface resistance
    integer, intent(in) :: deptype               ! deposition type
    integer, intent(in) :: timeSign              ! 1 for forward, -1 for backward
    real, dimension(:,:), intent(out) :: fVd, invVd2m ! resistance from zref to 2m 

    
    ! Local variables
    real :: fLambda, fCun, fWetDiam, fSc, fRoughfr, fStickRatio, fKin_visc, fDyn_visc, fDiffusivity, &
          & fSetlVelocity, fTau, fHumidity, invL, u_star, mol_diff, fZ0, fTemperature, &
          & fVTFplus, fSensHeatVelplus, vsplusCorrTF, Cp_ro, fTmp, R2m, V, fIvd
        real, parameter :: ascale = 2e-3  ! Collection scale is prescribed for
                                          ! land. Should be taken from landuse
        real, parameter :: fPrm = 0.7  ! Molecular prandtl number
        real, parameter :: fkgkp = 0.02  ! Ratio of heat conductivities
                                         ! for gas and particle
     real, parameter :: Z2m = 2.  ! Height fo 2m

    real :: tauplus, vsplus, rplus, Rsplus, invLplus, zplus,  Re, Sc, St, dpa, fKn !  dimensionless
    logical :: fwdOnly ! Can calculate vd only for forward time
                                                       !  parameters
    real :: vdsmooth, vdrough       ! vd for two landuse types. Should be more
    type(TwetParticle) :: wetParticle
    integer :: iSp, iSrc
    integer, save :: iCount=0

    call set_error("Not done yet!!","get_vd_species")
    return

        vdsmooth = 0.
        vdrough = 0.
        ! Meteo
        fDyn_visc = fu_dynamic_viscosity(pMetTempr2m(indexMeteo))
        fKin_visc = fDyn_visc * pMetSrfPressure(indexMeteo) / &
                  & (gas_constant_dryair * pMetTempr2m(indexMeteo))
        u_star =   pMetFricVel(indexMeteo)

        invL =     pMetMO_Len_inv(indexMeteo)
        fRoughfr = pMetLandFr(indexMeteo)
        !fZ0      =  pMetSrfRoughDisp(indexMeteo)
        fZ0      =  pMetSrfRoughMeteo(indexMeteo)
        
        Cp_ro = specific_heat_dryair*pMetSrfPressure(indexMeteo) / &
              & (gas_constant_dryair*pMetTempr2m(indexMeteo))
        fSensHeatVelplus = pMetSensHF(indexMeteo)/(Cp_ro * pMetTempr2m(indexMeteo) * u_star )
        !        - u_star*u_star*invL/(0.4*g) 
        ! <w'T'>/T/ u_*


            if(fZ0 < 1e-10)then
              if(iCount < 1000)then
                call msg('Strange z0:', fZ0)
                iCount = iCount + 1
              endif
              fz0=1.e-6
            endif

        if (fZ0 < 3. * fKin_visc/u_star) then !FIXME: May be some other coefficient 
                                              !needed here
                fRoughfr = 0.
        endif
                                               
        if(pMetPrecTot(indexMeteo) > 0.0) then ! Dimensionless meteo
                                                !switch for impaction 
                 fStickRatio = 1.0 
        else
                 fStickRatio = 0.0
        endif 
        
       !
       ! Dimensionless species properties
        if(fu_mode(speciesTransport(iSp)) == in_gas_phase)then    !gas
        ! Sc = nu / D (kinematic viscosity over molecular diffusivity)
            fSc = fKin_visc/fu_gas_molecular_diffusivity_air(speciesTransport(iSp)%material)
            fWetDiam = real_missing 
            tauplus = 0.
            vsplus = 0.
            rplus = 0.
            Rsplus = Rs(iSp)*u_star
            fTau = 0.
            fVTFplus = 0.
        else    ! particles
          ! No deposition rules here....
          !if(rulesDeposition%ifHumidityDependent)then
          ! FIXME: Surface-layer humidity should be here!
            fHumidity = fldRelHumidity%past%p2d(1)%ptr(indexMeteo) * weight_past + &
                 & fldRelHumidity%future%p2d(1)%ptr(indexMeteo) * (1.-weight_past)
            fTemperature = fldTempr%past%p2d(1)%ptr(indexMeteo) * weight_past + &
                 & fldTempr%future%p2d(1)%ptr(indexMeteo) * (1.-weight_past)
          !else
          !  fHumidity = rulesDeposition%fDefaultRelHumidity
          !endif
            wetParticle = fu_wet_particle_features( speciesTransport(iSp)%material, fHumidity)
            fWetDiam = wetParticle%fGrowthFactor * fu_massmean_D( speciesTransport(iSp)%mode )
            if(fWetDiam > 25e-6 .and. wetParticle%fGrowthFactor > 5)then
              fWetDiam = 25e-6
            endif
            fLambda = 2.37e-5 * pMetTempr2m(indexMeteo) / pMetSrfPressure(indexMeteo)
            fKn = 2*fLambda / fWetDiam
            fCun = 1. + fKn * (1.257 + 0.4 * exp(-1.1*fKn))
            
            fDiffusivity = boltzmann_const * pMetTempr2m(indexMeteo) * fCun / &
                         & (3. * Pi * fDyn_visc * fWetDiam)
            fTau =  fWetDiam * fWetDiam * wetParticle%fWetParticleDensity * &
                          & fCun / (18. * fDyn_visc)

            
            fVTFplus = - fSensHeatVelplus * fPrm  & !Thermophoretic velocity 
                &        * 2*fCun *1.17 *(fkgkp+2.18*fKn) & ! (Friedlander 2000)
                &        / (1+3*1.14*fKn)/(1.+2.*fkgkp +2.*2.18*fKn) * timeSign
            fSc = fKin_visc / fDiffusivity
            tauplus = fTau * u_star * u_star / ( fKin_visc)
            vsplus = fTau * g / u_star * timeSign
            rplus  = 0.5 * fWetDiam * u_star / fKin_visc   
            Rsplus = 0. ! No surface resistance for particles
          endif 




        
!! Temporary hack
!fu_get_vd = & !1e-8 +  fTau * g
!  & u_star * fu_vdplus_DS(tauplus,fSc,vsplus, Rsplus, invL*fZ0, zref/fZ0)
!return

    SELECT CASE(deptype)
      CASE(DryD_KS2011, DryD_KS2011_TF) ! Original KS scheme
        vsplusCorrTF=0. ! vorrection of VsPlus for thermophoresis
        if (deptype == DryD_KS2011_TF)  vsplusCorrTF = fvtfplus 
                                !  thermophoresis
        if  (fRoughfr < 1.) then ! smooth fraction exist
           fTmp = fKin_visc/u_star ! z0smooth
           vdsmooth = (1. - fRoughfr) * u_star * &
                fu_vdplus_smooth(tauplus, fSc, vsplus + vsplusCorrTF, rplus, Rsplus, invL*fTmp , zref(iSp, iSrc)/fTmp)
              if (.not. (vdsmooth .ge. 0.)) then
                 if (zref(iSp, iSrc) < fTmp) then
                         ! Could be just too-low center of masses
                         vdsmooth =  (1. - fRoughfr) * u_star ! whatever
                 else
                  call msg ("**DEPO** --get_vd_species------------------------")
                  call msg ("**DEPO** Negative vdsmooth!, WetDp=",fWetDiam)
                  call msg ("**DEPO** Negative vdsmooth!, tauplus=",tauplus)
                  call msg ("**DEPO** Negative vdsmooth!, fsc=",fsc)
                  call msg ("**DEPO** Negative vdsmooth!, vsplus=",vsplus)
                  call msg ("**DEPO** Negative vdsmooth!, rplus=",rplus)
                  call msg ("**DEPO** Negative vdsmooth!, Rsplus=",Rsplus)
                  call msg ("**DEPO** Negative vdsmooth!, Rs=",Rs)
                  call msg ("**DEPO** Negative vdsmooth!, vdsmooth=",vdsmooth)
                  call msg ("**DEPO** Negative vdsmooth!, zref=",zref(iSp, iSrc))
                  call msg ("**DEPO** Negative vdsmooth!, z0smooth=",fTmp)
                  vdsmooth = 0. ! Should not happen, but happens. 
                                ! sometimes returns NAN due to over-/under- flow
                endif
              endif 
       endif
        vsplusCorrTF = 0.5*(vsplusCorrTF+abs(vsplusCorrTF)) ! Can not be
                       ! negative for rough surfaces. 
                       ! Could be treated in some better way
        
        if (fRoughfr > 0.) then ! rough fraction exist
           Re =  ascale * u_star / fKin_visc
           St = 2*tauplus/Re * fStickRatio ! No impaction for non-sticky
                                          ! surfaces
           dpa = 2.*rplus/Re
           vdrough = fRoughfr * u_star * &
                fu_vdplus_rough(St, fSc, vsplus + vsplusCorrTF, dpa, Rsplus, Re, invL*fZ0 , max(zref(iSp, iSrc)/fZ0, 1.))
            if (.not. (vdrough .ge. 0.)) then
                  call msg("**DEPO** Bad vdrough in get_vd_species! vdrough,fz0:", vdrough, fZ0)
                  call msg('fRoughfr, u_star:',fRoughfr, u_star)
                  call msg('St, fSc',St, fSc)
                  call msg('vsplus, vsplusCorrTF',vsplus, vsplusCorrTF)
                  call msg('dpa, Rsplus',dpa, Rsplus)
                  call msg('Re, invL',Re, invL)
                  call msg('fZ0 , zref',fZ0 , zref(iSp, iSrc))
                  vdrough = 0. ! FIXME Should never happen
            endif
        endif          
        fVd = vdsmooth + vdrough 
        fwdOnly = .false.
        


      CASE (DryD_Vd_difsed) ! Vd ==  diffusion + settling
         fVd = u_star * &
               & fu_vdplus_DS(tauplus,fSc,abs(vsplus), Rsplus, invL*fZ0, zref(iSp, iSrc)/fZ0)
        fwdOnly = .true.

      CASE (DryD_Vd_SP98) ! Vd Slinn (Accoding to Seinfeld Pandis 2006)
         fVd = u_star * & 
               & fu_vdplus_slinnSP98(tauplus,fSc,abs(vsplus), Rsplus, invL*fZ0, zref(iSp, iSrc)/fZ0)
        fwdOnly = .true.

      CASE (DryD_VD_Zhang)! Generic Zhang... Grass for all
        fVd = fu_vd_Zhang(fTau,fSc,fWetDiam,u_star, Zref(iSp, iSrc), Rs(iSp), invL, ZhangGrass)
        fwdOnly = .true.

      CASE (DryD_VD_Zero)! 
        fVd = 0
        fwdOnly = .false.

      CASE DEFAULT ! Old scheme
            call msg("Unknown dry depo  rules! rules:", deptype) 
            call set_error("Unknown dry depo  rules", "fVd") 
    END SELECT

    ! correct Vd for inverse time
    if (fwdOnly .and. (timeSign < 0)) then
            fVd = fVd - u_star*abs(vsplus)
            fVd = 0.5 * (fVd + abs(fVd)) + 1e-10 ! 1/Vd is still
                                                                 ! meaningful
    endif

    !! Remaining stuff -- Vd2m

    invVd2m = 1./ fVd
    if (Z2m > fZ0) then  ! reasonable profile to 2m is possible
            R2m  = 2.5 * (log(Zref(iSp, iSrc)/Z2m) + fu_Psi(Zref(iSp, iSrc)*invL) - fu_Psi(Z2m*invL))
            ! R2m is positive for zref > Z2m and negative for  zref < Z2m
            V = g * fTau * timeSign
            invVd2m = 1./ fVd
            if (zref(iSp, iSrc)>Z2m) then 
                   R2m  = 0.5*(R2m+abs(R2m))   ! Could become negative in
                                             ! unstable
                   fTmp = V*R2m
                   if (abs(fTmp) .gt. 0.001) then
                        fTmp = exp(-fTmp)
                        invVd2m = invVd2m*fTmp +  (1.-fTmp)/V
                   else
                        invVd2m = invVd2m - R2m
                   endif
            else
                   R2m  = -0.5*(R2m-abs(R2m))   ! Could become wrong sign
                               ! Should be    positive here

                   fTmp = V*R2m
                   if (abs(fTmp) .gt. 0.001) then
                        fTmp = exp(fTmp)
                        invVd2m = invVd2m*fTmp +  (fTmp - 1.)/V
                   else
                        invVd2m = invVd2m + R2m
                   endif
            endif
    endif
    if (invVd2m(iSp, iSrc) < 0) then
            call msg("Negative Vd2m", invVd2m(iSp, iSrc))
    endif
    return
  end subroutine get_vd_species


!**********************************************************************************
! New deposition DryD_Vd_SP98 does the same at revision 57750  FIXME 
! To be removed?
!
!  subroutine get_Rb(speciesTransport, nSpecies, &                ! resistance for them
!                  & indexMeteo, weight_past, &  ! position in space and time
!                  & rulesDeposition, &           ! rules for standard deposition
!                  & arRb)                        ! output array for Rb
!    !
!    ! Returns an array of Rs for the subset of the given species - those, which are handled 
!    ! by the standard deposition procedure and thus referred to via the module private 
!    ! index arrays indGasesDepositing and indAerosolDepositing.
!    !
!    ! ATTENTION. Uses the precomputed arrays and parameters, so the species must be the 
!    !            very transport species, which were used for the precomputation step. This is 
!    !            checked upon the first call
!    !
!    implicit none
!
!    ! Imported parameters
!    type(silam_species), dimension(:), intent(in) :: speciesTransport
!    integer, intent(in) :: nSpecies, indexMeteo
!    type(Tdeposition_rules), intent(in) :: rulesDeposition
!    real, intent(in) :: weight_past
!    real, dimension(:), intent(out) :: arRb
!    
!    ! Local variables
!    integer :: iTmp, index
!    ! Local variables
!    real :: fLambda, fCun, fWetDiam, fSc, fSt, fStickRatio, fKin_visc, fDyn_visc, fDiffusivity, &
!          & fSetlVelocity, fHumidity, u_star, prandtl, mol_diff
!    type(TwetParticle) :: wetParticle
!
!    !
!    ! Actual work starts
!    !
!    select case(rulesDeposition%DryDepType)
!
!      case(DryD_grav_only)
!        do iTmp = 1, rulesDeposition%nGasesDepositing
!          arRb(rulesDeposition%indGasesDepositing(iTmp)) = -1     ! no deposition for gases
!        end do
!        do iTmp = 1, rulesDeposition%nAerosolsDepositing
!          arRb(rulesDeposition%indAerosolDepositing(iTmp)) = 1.e10  ! only settling is considered for aerosols
!        end do
!
!      case(DryD_Rb_only, DryD_resist_only, DryD_Rb_grav_combine, DryD_resist_grav_combine)
!
!        fLambda = 2.37e-5 * pMetTempr2m(indexMeteo) / pMetSrfPressure(indexMeteo)
!
!        fDyn_visc = fu_dynamic_viscosity(pMetTempr2m(indexMeteo))
!
!        fKin_visc = fDyn_visc * pMetSrfPressure(indexMeteo) / &
!                  & (gas_constant_dryair * pMetTempr2m(indexMeteo))
!
!        !
!        ! For gases:
!        ! Hicks et al. (1987): Rb = 2 / (u_star * karmann) * (Sc/Pr)**2/3
!        ! Alternatively in Seinfeld & Pandis: Rb = 5*Sc**2/3 / u_star
!        ! (attributed to Wesely (1989) but mentioned there).
!        ! Sc = nu / D (kinematic viscosity over molecular diffusivity)
!        !
!        do iTmp = 1, rulesDeposition%nGasesDepositing
!!          call msg('iTmp, indGasesDepositing(iTMp)',iTmp,rulesDeposition%indGasesDepositing(iTmp))
!!          call msg('Reporting species:')
!!          call report(speciesTransport(rulesDeposition%indGasesDepositing(iTmp)))
!          
!          if(fu_if_gas_depositing(speciesTransport(rulesDeposition%indGasesDepositing(iTmp))%material))then
!            u_star = pMetFricVel(indexMeteo)
!            prandtl = pMetPrandtl(indexMeteo)
!            index = rulesDeposition%indGasesDepositing(iTmp)
!            mol_diff = fu_gas_molecular_diffusivity_air(speciesTransport(index)%material)
!            arRb(index) = 2.0 / (u_star*karmann_c) * (fKin_visc / (mol_diff*prandtl))**0.666667            
!          else
!            arRb(rulesDeposition%indGasesDepositing(iTmp)) = -1.
!          endif
!        end do
!        !
!        ! For aerosols:
!        ! 1/Rb = u* (Sc**(-2/3) + 10**(-3/St))
!        ! Sc=nu/D, nu=1.46e-5 m2/s, D = kT Cun/(3pi d mu), d is 
!        ! particle diameter, mu is dynamic viscosity=1.8e-5, k=1.38e-23 J/K - 
!        ! Boltzman's constant, Cun is no-slip correction
!        ! St is Stokes number (inertial passing through layer), St=v_sedim * u* * u* /(g*nu)
!        !
!        ! Alternative from Seinfield & Pandis:
!        ! 1/Rb = 3 * u* (Sc**-gamma + (St/alfa+St)**2 + 0.5(Dp/A)**2)
!        ! In principle, gamma, alfa and A are land-use dependent but some mean values are:
!        ! gamma = 0.55, alpha=1, A=5mm
!        !
!        if(rulesDeposition%nAerosolsDepositing > 0)then
!          if(rulesDeposition%ifHumidityDependent)then
!            fHumidity = fldRelHumidity%past%p2d(1)%ptr(indexMeteo) * weight_past + &
!                      & fldRelHumidity%future%p2d(1)%ptr(indexMeteo) * (1.-weight_past)
!          else
!            fHumidity = rulesDeposition%fDefaultRelHumidity
!          endif
!
!          do iTmp = 1, rulesDeposition%nAerosolsDepositing
!
!            wetParticle = fu_wet_particle_features( &
!                           & speciesTransport(rulesDeposition%indAerosolDepositing(iTmp))%material, &
!                           & fHumidity)
!
!            fWetDiam = wetParticle%fGrowthFactor * &
!                     & fu_massmean_D(speciesTransport(rulesDeposition%indAerosolDepositing(iTmp))%mode)
!            if(fWetDiam > 25e-6 .and. wetParticle%fGrowthFactor > 5)then
!              fWetDiam = 25e-6
!            endif
!            fCun = 1. + 2*fLambda / fWetDiam * (1.257 + 0.4 * exp(-0.55 * fWetDiam / fLambda))
!            fDiffusivity = boltzmann_const * pMetTempr2m(indexMeteo) * fCun / &
!                         & (3. * Pi * fDyn_visc * fWetDiam)
!            fSc = fKin_visc / fDiffusivity
!            fSetlVelocity = g * fWetDiam * fWetDiam * wetParticle%fWetParticleDensity * &
!                          & fCun / (18. * fDyn_visc)
!
!            fSt = fSetlVelocity * pMetFricVel(indexMeteo) * pMetFricVel(indexMeteo) / (g * fKin_visc)
!
!            !
!            ! Stick-ratio decides whether the particle stays in air after hitting the surface
!            ! exp(-sqrt(fSt)) for dry 1 for wet surfaces
!            !
!!            call msg('Large scale and convective precipitation:',ptrLrgScalePrec(indMeteo), ptrConvPrec(indMeteo))
!            if(pMetPrecTot(indexMeteo) > 0.0)then
!              fStickRatio = 1.0
!            else
!              fStickRatio = 1.0 * (1.0-pMetLandFr(indexMeteo)) + exp(-sqrt(fSt)) * pMetLandFr(indexMeteo)
!            endif
!            !
!            ! Now - two options. This is from Slinn.
!            !
!            arRb(rulesDeposition%indAerosolDepositing(iTmp)) = 1. / &
!                              & (pMetFricVel(indexMeteo) * (fSc**(-0.66666667) + 10**(-3.0/fSt)))
!            !
!            ! Option 2 seems to be better but the interception and stokes number require detailed land use, 
!            ! which is still in the future. So far, let's just distinguish between the land and sea
!            !
!! 14.1.2011.
!! This option is from Zhang(2001) and, in fact, is rubbish. The deposition velocity is much too high.
!!
!!            arRb(rulesDeposition%indAerosolDepositing(iTmp)) = & 
!!                          &   1. / (3. * pMetFricVel(indexMeteo) * fStickRatio * &
!!                                  & (fSc**(-0.55) + &                                ! Brownian diffusion
!!                                   & fSt*fSt/((1+fSt)*(1+fSt)) + &                 ! Impaction
!!                                   & pMetLandFr(indexMeteo)*fWetDiam*fWetDiam/5e-5))  ! Interception
!
!          end do  ! depositing aerosols
!        endif  ! if any aerosol
!      case default
!        call set_error('Unknown dry deposition type','fu_R_b_aerosol')
!    end select
!
!  end subroutine get_Rb

  !**********************************************************************************

  subroutine get_Rs_2013(speciesTransport, mmr_lowest_lev, nSpecies, &        ! resistance for them
                  & indexMeteo, weight_past, &           ! position in space and time
                  & rulesDeposition, &                   ! rules for standard deposition
                  & arRs, ifTuned)                                ! output array for Rs
                  !& arRs, massair,massDep)                                ! output array for Rs
    !
    ! Returns an array of Rs for the subset of the given species - those, which are handled 
    ! by the standard deposition procedure and thus referred to via the module private 
    ! index arrays indGasesDepositing and indAerosolDepositing.
    !
    ! ATTENTION. Uses the precomputed arrays and parameters, so the species must be the 
    !            very transport species, which were used for the precomputation step
    !
    !  Losely adapted from EMEP unified model
    !

    implicit none

    ! Imported parameters
    type(silam_species), dimension(:), intent(in) :: speciesTransport
    real, dimension(:), intent(in) :: mmr_lowest_lev 
    integer, intent(in) :: nSpecies, indexMeteo
    type(Tdeposition_rules), intent(in) :: rulesDeposition
    real, intent(in) :: weight_past
    !    real, dimension(:), intent(in) :: massAir, massDep
    real, dimension(:), intent(out) :: arRs
    logical, intent(out) :: ifTuned

    ! Local variables
    type(Tgas_deposition_param), pointer :: DepData
    real :: alphaSNdep, alphaSNair
    integer :: iTmp, iSpecies
    real :: rWater, rSoil, fWet, fHumidity
    real :: fice, fsnow,lowTcorr,fZ0,sdepth,hveg,pressure,t2m,t2c,rh2m,LAI,SAI,landfr,u_star, g_sto
    real :: Sdmax,RsnowS,RsnowO, F1, F2, Hstar, drx, fTmp
    real :: GigsO, RigsO, Gns, GnsO, RnsO, Gmesophyl
    real, parameter :: BETA = 1.0/22.0
    logical :: canopy, leafy_canopy, is_veg
    real :: GnsS,Rinc,Rns_NH3,Rns_SO2
    real :: RgsOsfc, RgsSsfc
!! external resistance for Ozone
  real, parameter :: RextO =  2500.0   ! gives Gext=0.2 cm/s for LAI=5
!!  real, parameter :: RextO =  10000.0   ! Try to adjust

    real, parameter :: mmrSO2sat = 10.e-6 /1.2  / .064 !10 ug/m3 / 1.2kg/m3 / .064 kg/mole

!! Here, "Gext=0.2cm/s" refers to the external conductance, G_ext, where 
!! G_ext=LAI/R_ext. In many studies, it has been assumed 
!! that G_ext should be low, particularly relative to stomatal conductance g_s.
!! Results from a variety of experiments, however, have made the above 
!! estimates  Rext0 and RextS plausible.  The above equation for G_ext has been
!! designed on the basis of these experimental results. 
!
!! Notice also that given the equations for the canopy resistance R_sur and the 
!! deposition velocity V_g, V_g>=LAI/R_ext. The value of G_ext can therefore be
!! interpreted as the minimum value for V_g.
    
    !
    ! Initially, all is negative. Those, which deposit, will be set to something reasonable
    arRs(1:nSpecies) = -1.0

    if (rulesDeposition%DryDepType == DryD_Vd_Zero) return


    is_veg = .true.
!    alphaSNdep = 100. ! ~ alphaSNair = 1
    alphaSNdep = 1000 !! Increase Non-stomatal resistance

    alphaSNair = 1.
!    if (indSO2 > 0) then
!            if (indNH3 > 0) then
!                    alphaSNdep = min(max(massdep(indSO2)/massdep(indNH3),3.), 2e5) ![3:2e5]
!                    alphaSNair = min(massair(indSO2)/massair(indNH3), 3.) ! [0:3]
!            else ! no ammonia at all...
!                    alphaSNdep = 2e5
!                    alphaSNair = 3.
!            endif
!    endif



    !
    ! Rs for aerosols is zero
    do iTmp = 1, rulesDeposition%nAerosolsDepositing
      arRs(rulesDeposition%indAerosolDepositing(iTmp)) = 0
    end do

    if  (rulesDeposition%nGasesDepositing < 1) return


    fZ0      =  pMetSrfRoughMeteo(indexMeteo)
    hveg     = 10 * fZ0 ! Vegetation height used for snow cover, incanopy
                        ! resistance etc
    u_star = pMetFricVel(indexMeteo)
    pressure = pMetSrfPressure(indexMeteo) 
    t2m      = pMetTempr2m(indexMeteo)
    t2c      = t2m - 273.15 !In Celsius
    rh2m     = pMetRh2m(indexMeteo)
    sdepth   = pMetSnowDepth(indexMeteo)* 5.  !! Water eq. to normal depth
    LAI      = pMetLAI(indexMeteo)
    if (LAI < 0.) LAI=0  !Can be missing or something...
    g_sto    = pMetGsto(indexMeteo)
    landfr   = pMetLandFr(indexMeteo)
    !fice     = pMetIceFr(indexMeteo) 

    
    if (hveg > 2.) then ! high vegetation
        SAI = LAI + 1 !! Surface area index of canopy
        canopy = .true.
    else
        SAI = LAI
        canopy = .false.
    endif
    if (LAI > 0.1) then
            leafy_canopy = .true.
    else
            leafy_canopy = .false.
    endif

    landfr = pMetLandFr(indexMeteo)

    ! Get snow
    Sdmax = max( hveg/10.0, 0.01) !meters
    fsnow = 2.0 *sdepth/Sdmax

!Treat ice in the same way as snow
   ! fsnow = max(fsnow,fice) !if snow_flag, ice_nwp probably has snow_flag
                            !but it might be ice without snow..
    fsnow = min(fsnow,1.0)
    fsnow = max(fsnow,0.0)
 
  !===========================================================================
  !  Adjustment for low temperatures (Wesely, 1989, p.1296, left column)

    lowTcorr = exp(0.2*(-1 - t2c))!Zhang,2003 & Erisman 1994
    lowTcorr = min(2.0,lowTcorr)   !Zhang,2003 & Erisman 1994
    lowTcorr = max(1.0,lowTcorr)   !Zhang,2003 & Erisman 1994
                                   !effectivelluy means that it will only
                                   !kick in for T<-1

! Rsnow for sulphur and O3, Erisman, 1994 + Zhang, 2003. Also used for ice. 
    RsnowS = 70.0*(2.0 - t2c) !Used for snow_flag and ice_nwp
    if (t2C < -1.0) RsnowS = 700.0 !700 from Cadle,1985
    RsnowS = min(700.0,RsnowS) !Erisman 1994=500,very low.. Puts to 2000
    RsnowS = max(70.0,RsnowS)  !Erisman 1994. 70 above 1 degree

    RsnowO = 2000.0 !same for snow_flag, ice_nwp, water. Later corrected with lowTcorr
                    !as recommended by Juha-Pekka


  !===========================================================================
  ! Get Rns_SO2

!       call CoDep_factors(G%so2nh3ratio24hr,G%so2nh3ratio,&
!              L%t2C,L%rh,L%is_forest, debug_flag)
    if (T2c > 0 ) then    ! Use "rh" - now in fraction 0..1.0
        
      F1 = 10.0 * log10(T2C+2.0) * exp(100.0*(1.0-rh2m)/7.0) ! tab_exp_rh(IRH)
      F2 = 48 * exp(-2.53 * alphaSNair)
      !a_SN =  ia_SN * MAX_SN /real(NTAB)
      !ia_SN = nint( NTAB * a_SN/MAX_SN )   ! Spread values from 0-3 to 0:100
      !tab_F2 (ia_SN)            = 10.0**( (-1.1099 * a_SN)+1.6769 )
      !120 * alphaSNdep**(-0.7)  
      !!tab_F2( iaSN  )
      Rns_NH3 = BETA * F1 * F2
      Rns_NH3 = min( 200.0, Rns_NH3)  ! After discussion with Ron
      Rns_NH3 = max(  10.0,Rns_NH3)   ! Sic! FIXME Does not depend on LAI/SAI 
        
      !Rns_SO2_dry = 11.84  * exp(1.1*so2nh3ratio24hr) * ( frh**(-1.67) )
      Rns_SO2 = 8.0 * alphaSNdep **0.3 * ( rh2m**(-1.67) )
      Rns_SO2 = min( 1000.0, Rns_SO2)  ! Set because rh goes almost 
      ! to zero occationally over the Alps, Greenland and Svalbard
      ! 1000 chosen sort of random
      Rns_SO2 = max(  10.0,Rns_SO2) !hf CoDep SHOULD WE LIMIT IT to 10??

      !Attempt to reduce winter deposition of SO2

!    else if (  T2c > -5 ) then
!      ! In a future version, we might test for  Ts_C > -2 instead of
!      ! Ts_C > 0  as we have now.
!
!      Rns_NH3 = 100.0 ! will be modified by low-T factor later
!      Rns_SO2 = 100.0
    else
      Rns_NH3= 500.0 ! will be modified by low-T factor later
      Rns_SO2= 500.0

    end if !Ts_C

  !** Calculate Rinc, Gext 

    if(  canopy ) then

      Rinc = 14.0 * SAI * hveg  / u_star    ! Erisman's b.LAI.h/u*

      ! for now, use CEH stuff for canopies,and soils (canopies ouside 
      ! growing season)
      ! keep Ggs for non-canopy

      GnsS = (1.-fsnow)/(Rns_SO2 * lowTcorr) + fsnow/RsnowS

    elseif  ( is_veg ) then ! vegetation outside growing season

      Rinc = 0.0
      GnsS = (1.-fsnow)/(Rns_SO2 * lowTcorr) + fsnow/RsnowS


    else   ! No canopy or soil present

      Rinc = 0.0

      !/ Here we preserve the values from the ukdep_gfac table
      !  giving higher deposition to water, less to deserts

      !SILAM value, map is used in EMEP model
      
      !Again error in EMEP
      !RgsSsfc = 50 * (1. - landfr) + 500 * landfr  ! Default for bare surface


      RgsSsfc = 1./(  (1. - landfr)/50. +  landfr / 200.)  ! Default for bare surface
                !gas_surf_resistance_over_water = 0.5 # Replaced with 50
                !gas_surf_resistance_land_default = 100. # 500 as typical EMEP

      GnsS = (1.-fsnow)/(RgsSsfc * lowTcorr) + fsnow/RsnowS

    end if !  canopy

    !!snow treated as in Zhang 2003
    !!But Zhang wse 2*fsnow for ground surface because Sdmax(snow depth when total coverage is assumed)
    !!for soils under vegetation is assumed to stay snow covered longer than 'the leafs'
    !!but - we have underlying surfaces only for O3 and for simplicity we treat them equally
    !!RECONSIDER THIS ESPECIALLY BASED ON SATELITTES

    !!no snow corrections (or low temperature) for Rinc 
    !!RgsO 'corrected for snow' and low temp
    !!as adviced by Juha-Pekka

    !SILAM value, map is used in EMEP model
    !!! FIXME BULLSHIT from EMEP!! Can't add resistances this way!!!
    !!RgsOsfc = 2000.0 * (1. - landfr) + 200. * landfr ! Default for bare surface
    !! More correct way
!    RgsOsfc = 1./( (1. - landfr)/2000. +  landfr / 200.) ! Default for bare surface
!ifTuned = .False.
    RgsOsfc = 1./( (1. - landfr)/2000. +  landfr / 1000.) ! Default for bare surface
ifTuned = .true.
!        Not particulary sensitive
!  gas_surf_resistance_over_water = 2000.0
!  #gas_surf_resistance_land_default = 250.0 # 200 as typical EMEP
!  gas_surf_resistance_land_default = 1500.0

    GigsO=  (1.-fsnow)/RgsOsfc   + fsnow/RsnowO
    RigsO = lowTcorr/GigsO +  Rinc

!!####   2. Calculate Surface Resistance, Rsur, for HNO3 and Ground Surface 
!!####      Resistance, Rgs, for the remaining Gases of Interest                                

   !!/ Ozone values....

        !!RextO corrected for low temp
        !!as adviced by Juha-Pekka
    GnsO   = SAI/(RextO * lowTcorr) + 1.0/ RigsO     ! (SAI=0 if no canopy)

    do iTmp = 1, rulesDeposition%nGasesDepositing

      iSpecies = rulesDeposition%indGasesDepositing(iTmp)
      DepData => fu_gas_deposition_param(speciesTransport(iSpecies)%material)
      if (DepData%Henry_const_298K < 0.) cycle ! No deposition for the gas
     
      if (iSpecies == indHNO3) then
        arRs(iSpecies) = max(10., -2.*T2c) ! Simpson 2012, eq.56
        if (fu_fails(arRs(iSpecies) >= 0., "Negative surface resistance HNO3","get_Rs_2013"))return
        cycle
      endif

      Hstar = DepData%Henry_const_298K * &
            & exp (DepData%Henry_const_T_dep * ((1./t2m)-(1./298.)) )
      
      !Diffusivity ratio D(x)/D(H20)
      ! In EMEP more diffusive stuff deposites slower!!! Corrected here. Roux
      !DRx =  fu_gas_molecular_diffusivity_air(speciesTransport(iSpecies)%material) / 2e-5
      DRx =  fu_gas_molecular_diffusivity_air(speciesTransport(iSpecies)%material) / 2.178e-5 !RH2018 update: Using the Massman 1998 value for H2O!

      if (iSpecies == indNH3) then
        Gns = (1.-fsnow)/(Rns_NH3 * lowTcorr) + fsnow/RsnowS
      else
        Gns = Hstar * 0.67e-5 * GnsS + DepData%Wesely_f0 * GnsO 
      endif
      !
      ! Stomatal resistance is a sum of the leaf stomatal resistance and mesophyll resistance
      ! See Seinfeld & Pandis 2006, p.921.
      ! Note that in Eq. (19.52) the Henry law contant should be in units of M/atm but
      ! SILAM has it in units of mol/kg/Pa. Therefore an extra factor of 1.01e5 is needed. 
      !
      Gmesophyl = 3.3e-4 * 1.01e5*Hstar + 100.0 * DepData%Wesely_f0

!      fTmp = max( LAI * DRx * G_sto + Gns, 1e-10 ) !1/Rs, misses mesophyl resistor
!      fTmp = max( LAI / (1./(DRx * G_sto) + 1. / Gmesophyl) + Gns, 1e-10 ) !1/Rs
      fTmp = max( LAI * DRx * G_sto * Gmesophyl / (Gmesophyl + DRx*G_sto) + Gns, 1e-10 ) !1/Rs
            
      arRs(iSpecies) = 1.0 / fTmp
      if (.not. (arRs(iSpecies) >= 0.)) then
        call msg("Non-positive surface resistance for:" + fu_name(speciesTransport(iSpecies)%material), arRs(iSpecies))
        call set_error("Negative surface resistance","get_Rs_2013")
      endif

      if (iSpecies == indSO2) then
        ! Reduce deposition of very dense plumes
        arRs(iSpecies) = arRs(iSpecies) * ( mmr_lowest_lev(iSpecies) + mmrSO2sat ) / mmrSO2sat
      endif

    end do  ! Gaseous species

  end subroutine get_Rs_2013


  !**********************************************************************************

  subroutine get_Rs(speciesTransport, nSpecies, &        ! resistance for them
                  & indexMeteo, weight_past, &           ! position in space and time
                  & rulesDeposition, &                   ! rules for standard deposition
                  & arRs)                                ! output array for Rs
    !
    ! Returns an array of Rs for the subset of the given species - those, which are handled 
    ! by the standard deposition procedure and thus referred to via the module private 
    ! index arrays indGasesDepositing and indAerosolDepositing.
    !
    ! ATTENTION. Uses the precomputed arrays and parameters, so the species must be the 
    !            very transport species, which were used for the precomputation step
    !
    implicit none

    ! Imported parameters
    type(silam_species), dimension(:), intent(in) :: speciesTransport
    integer, intent(in) :: nSpecies, indexMeteo
    type(Tdeposition_rules), intent(in) :: rulesDeposition
    real, intent(in) :: weight_past
    real, dimension(:), intent(out) :: arRs

    ! Local variables
    integer :: iTmp
    real :: rWater, rSoil, fWet, fHumidity

    !
    ! Initially, all is negative. Those, which deposit, will be set to something reasonable
    !
    arRs(1:nSpecies) = -1.0

    if (rulesDeposition%DryDepType == DryD_Vd_Zero) return

!call msg('Thread reset completed',OMP_GET_THREAD_NUM())

    !
    ! Rs for aerosols is zero
    !
    do iTmp = 1, rulesDeposition%nAerosolsDepositing
      arRs(rulesDeposition%indAerosolDepositing(iTmp)) = 0
    end do

    do iTmp = 1, rulesDeposition%nGasesDepositing

      if(pMetLandFr(indexMeteo) > 0.01)then
        !
        ! Land exists. Apply the standard soil resistance (so far!) and the correction coefficient for the 
        ! soil wetness. The correction can be non-integer but in general, it is 1 for dry conditions,
        ! 2 for slight rain and 3 for strong rain. Algorithm: for dry case take the soil parameter, for
        ! slight rain a mixture of soil and water, for strong rain - water parameter.
        ! Note that resistance onto water can be negative, which means no deposition. Take care of it.
        !
        fHumidity = fldRelHumidity%past%p2d(1)%ptr(indexMeteo) * weight_past + &
                  & fldRelHumidity%future%p2d(1)%ptr(indexMeteo) * (1.-weight_past)
        
        
        !
        ! Procedure: 
        ! 1. Strong rain => water coefficeint
        ! 2. Mild rain OR high humidity with not too high temperature => intermediate value
        ! 3. Else: (No rain, low humidity or high temperature) => dry value
        !
        if(pMetTempr2m(indexMeteo) > 270.)then
          if(pMetPrecTot(indexMeteo) > 3.0e-4)then
            !
            ! Above 1mm/hr, all wet
            !
            fWet = 1.
            rWater = fu_gas_Rs_default_water( &
                          & speciesTransport(rulesDeposition%indGasesDepositing(iTmp))%material)
            rSoil = rWater
          else
            !
            ! No strong rain. Mild one?
            !
            fWet = pMetPrecTot(indexMeteo)/3.e-4
            if(fHumidity > 0.7)then
              !
              ! High humidity, for low temperature probably means wet surface
              !
              if(pMetTempr2m(indexMeteo) > 270. .and. pMetTempr2m(indexMeteo) < 285.) &
                                                         & fWet = max(fWet,(fHumidity - 0.7) / 0.3)
            endif
            !
            ! Using the wettness parameter, get the Rs keeping in mind hydrophobic/filic species
            !
            if(fWet > 1.e-5)then
              rWater = fu_gas_Rs_default_water( &
                          & speciesTransport(rulesDeposition%indGasesDepositing(iTmp))%material)
              if(rWater < 0)then
                !
                ! For species not depositing on water, increase the resitance with stronger rain
                !
                rSoil = fu_gas_Rs_default_soil( &
                            & speciesTransport(rulesDeposition%indGasesDepositing(iTmp))%material) * &
                      & exp(4.*fWet)
              else
                rSoil = fWet * rWater + (1.-fWet) * &
                                  & fu_gas_Rs_default_soil(speciesTransport( &
                                               & rulesDeposition%indGasesDepositing(iTmp))%material)
              endif
            else
              rSoil = fu_gas_Rs_default_soil( &
                            & speciesTransport(rulesDeposition%indGasesDepositing(iTmp))%material)
            endif  ! fWet
          endif  ! precipitation
        else
          !
          ! All frozen
          !
          fWet = 0.
          rSoil = fu_gas_Rs_default_soil( &
                          & speciesTransport(rulesDeposition%indGasesDepositing(iTmp))%material)
        endif
        !
        ! Now the fraction over water
        !
        if(pMetLandFr(indexMeteo) < 0.99 .and. pMetTempr2m(indexMeteo) > 270.)then
          !
          ! Water exists. The standard coefficient is there
          !
          if(fWet <= 1e-5)rWater = fu_gas_Rs_default_water( &
                             & speciesTransport(rulesDeposition%indGasesDepositing(iTmp))%material)
          arRs(rulesDeposition%indGasesDepositing(iTmp)) = 0.0  ! turn into conductivity

          if(rWater > 0) arRs(rulesDeposition%indGasesDepositing(iTmp)) = &
                               & arRs(rulesDeposition%indGasesDepositing(iTmp)) + &
                                                        & (1. / rWater) * (1.-pMetLandFr(indexMeteo))
          if(rSoil > 0) arRs(rulesDeposition%indGasesDepositing(iTmp)) = &
                               & arRs(rulesDeposition%indGasesDepositing(iTmp)) + &
                                                       & (1. / rSoil) * pMetLandFr(indexMeteo)
          if(arRs(rulesDeposition%indGasesDepositing(iTmp)) > 0.0001)then
            arRs(rulesDeposition%indGasesDepositing(iTmp)) = 1. / arRs(rulesDeposition%indGasesDepositing(iTmp))
          else
            arRs(rulesDeposition%indGasesDepositing(iTmp)) = -1.0  ! no deposition
          endif
        else
          !
          ! No water, just soil
          !
          arRs(rulesDeposition%indGasesDepositing(iTmp)) = rSoil

        endif   ! water exists

      else
        !
        ! No land, just water
        !
        arRs(rulesDeposition%indGasesDepositing(iTmp)) = &
               & fu_gas_Rs_default_water(speciesTransport(rulesDeposition%indGasesDepositing(iTmp))%material)

      endif   ! land exists
    end do  ! Gaseous species

  end subroutine get_Rs


  !*************************************************************************************

  subroutine get_settling_velocity_species(speciesTransport, nSpecies, &
                                         & indexMeteo, iLevMeteo, weight_past, &
                                         & rulesDeposition, &
                                         & velocity, &
                                         & forced_humidity)
    !
    ! Computes the general settling velocity for the given substance and given mode.
    ! Searches the right metedata by itself
    !
    !
    ! ATTENTION. NO CHECKING to make the story fast
    !
    implicit none
    type(silam_species), dimension(:), intent(in) :: speciesTransport
    integer, intent(in) :: nSpecies, indexMeteo, iLevMeteo ! horizontal and vertical coordinates
    real, intent(in) :: weight_past
    type(Tdeposition_rules), intent(in) :: rulesDeposition
    real, dimension(:), intent(inout):: velocity
    real, intent(in), optional :: forced_humidity
    
    ! Local variables
    real :: fTempr, fPressure, fWetDiam, mu_air
    type(TwetParticle) :: WetParticle
    integer :: iTmp

    ! The settling velocity is taken as a function of air dynamic viscosity and aerosol particle
    ! characteristics. Note that for gases it is zero, of course, as well as for non-depositing species
    !
    velocity(1:nSpecies) = 0.0
!    do iTmp = 1, rulesDeposition%nGasesDepositing
!      velocity(rulesDeposition%indGasesDepositing(iTmp)) = 0
!    end do

    if(rulesDeposition%nAerosolsDepositing == 0)return

    !
    ! Aerosols exist. Some preparation and proceed
    !
    if(iLevMeteo == 0)then
      fTempr = pMetTempr2m(indexMeteo)
      fPressure = pMetSrfPressure(indexMeteo)
    else
      fTempr = fldTempr%past%p2d(iLevMeteo)%ptr(indexMeteo) * weight_past + &
             & fldTempr%future%p2d(iLevMeteo)%ptr(indexMeteo) * (1.-weight_past)
      fPressure = fldPressure%past%p2d(iLevMeteo)%ptr(indexMeteo) * weight_past + &
                & fldPressure%future%p2d(iLevMeteo)%ptr(indexMeteo) * (1.-weight_past)
    endif

    mu_air = fu_dynamic_viscosity(fTempr)
    do iTmp = 1, rulesDeposition%nAerosolsDepositing
      !
      ! Take the aerosol growth into account
      !
      if(rulesDeposition%ifHumidityDependent)then
        if(present(forced_humidity))then
          wetParticle = fu_wet_particle_features( &
                             & speciesTransport(rulesDeposition%indAerosolDepositing(iTmp))%material, &
                                              & forced_humidity)
        else
          if(iLevMeteo == 0)then
            wetParticle = fu_wet_particle_features( &
                               & speciesTransport(rulesDeposition%indAerosolDepositing(iTmp))%material, &
                               & fldRelHumidity%past%p2d(1)%ptr(indexMeteo) * weight_past + &
                               & fldRelHumidity%future%p2d(1)%ptr(indexMeteo) * (1.-weight_past))
          else
            wetParticle = fu_wet_particle_features( &
                            & speciesTransport(rulesDeposition%indAerosolDepositing(iTmp))%material, &
                               & fldRelHumidity%past%p2d(iLevMeteo)%ptr(indexMeteo) * weight_past + &
                               & fldRelHumidity%future%p2d(iLevMeteo)%ptr(indexMeteo) * (1.-weight_past))
          endif
        endif  ! forced_humidity present
      else
        !
        ! No relative humidity dependence
        !
        wetParticle = fu_wet_particle_features( &
                             & speciesTransport(rulesDeposition%indAerosolDepositing(iTmp))%material, &
                                              & rulesDeposition%fDefaultRelHumidity)
      endif

      fWetDiam = fu_massmean_d(speciesTransport(rulesDeposition%indAerosolDepositing(iTmp))%mode) * &
               & WetParticle%fGrowthFactor

      velocity(rulesDeposition%indAerosolDepositing(iTmp)) = &
              & fu_settling_vel(fWetDiam, WetParticle%fWetParticleDensity, mu_air, fTempr, fPressure)
    end do ! aerosol species

  end subroutine get_settling_velocity_species

 !*******************************************************

real function fu_settling_vel(dp, rhop, mu_air, T, P) result(vs)
  !
  ! Implements settling velocity for the whole range of particles
  !
  ! Implements Solution for equation taht F_g + F_stokes + F_drag = 0
  ! F_g = - pi/18  rho_p g  dp**3    
  ! F_stokes = 3 pi mu_air dp vs / Cun
  ! where Cun = 1 + 2*lambda/d * (1.26 + 0.4 * exp(-0.55 * d/lambda))
  ! lambda = 2.18e-5 * T / P ## Constant to match eq 9.7 of SeinfeldPandis2006, p.399
            !! Earlier was 2.37e-5 from unknown source
  ! F_drag = pi/8 dp**2 rho_air C_d
 
    implicit none
    real, intent(in) :: dp, rhop, mu_air, T, P
    real :: A, B, C, Det, lambda_over_dp, fCun, rho_air
    character(len = *), parameter :: sub_name = 'fu_settling_vel'
  
      rho_air = P / (gas_constant_dryair * T)
    ! Equation: A * vs * vs + B * vs + C =0 
      A = rho_air * dp * 0.45  !! A = F_t * 6 /pi / U^2 /dp  Turbulent drag
                                 !! Drag coefficient for a sphre Cd=0.45
      B = 18 * mu_air   ! B = F_s *6 /pi  / U /dp ##Stokes drag without Cunningham
      C = - 0.75 *  rhop  * g * dp * dp !! C = Fg *6 /pi /dp ##Gravity

      if (abs(A*C) < 1e-4 * B*B) then !! Stokes dominates -- linear case
        lambda_over_dp = 2.18e-5 * T / (P*dp) !!! lambda/d_p
        fCun = 1.+2*lambda_over_dp * (1.26 + 0.4 * exp(-0.55 / lambda_over_dp))
        vs = - fCun * C / B
      else  !! quadractic can be solved without numerical issues
        Det = B*B - 4*A*C !! No cunningham for large particles (hopefully)
        vs = 0.5*(-B + sqrt(Det)) / A
      endif
  
end function fu_settling_vel


  !*************************************************************************************

  subroutine get_settling_motion_species(speciesTransport, nSpecies, &
                                         & fTempr, fPressure, fRhoAir, fRelHum, &
                                         & rulesDeposition, &
                                         & velocity, ifPressure)
    !
    ! Computes the general settling velocity in m/s or in Pa/s 
    ! for the given substance and given mode.
    !
    ! Positive -- down!!!!
    !
    ! Formula: v = g * ro_part * d_part**2 * Cun / (18 * dyn_visc)
    ! where Cun = 1 + 2*lambda/d * (1.26 + 0.4 * exp(-0.55 * d/lambda))
    ! lambda = 2.37e-5 * T / P
    !
    !
    implicit none
    type(silam_species), dimension(:), intent(in) :: speciesTransport
    integer, intent(in) :: nSpecies
    real, intent(in) :: fTempr, fPressure, fRhoAir, fRelHum
    logical, intent(in) :: ifPressure
    type(Tdeposition_rules), intent(in) :: rulesDeposition
    real, dimension(:), intent(out):: velocity
    
    ! Local variables
    real :: mu_air, fWetDiam
    type(TwetParticle) :: WetParticle
    integer :: iTmp, iSp

    ! The settling velocity is taken as a function of air dynamic viscosity and aerosol particle
    ! characteristics. Note that for gases it is zero, of course, as well as for non-depositing species
    !
    velocity(1:nSpecies) = 0.0
    do iTmp = 1, rulesDeposition%nGasesDepositing
      velocity(rulesDeposition%indGasesDepositing(iTmp)) = 0
    end do

    if(rulesDeposition%nAerosolsDepositing == 0)return

    !
    ! Aerosols exist. Some preparation and proceed

    mu_air = fu_dynamic_viscosity(fTempr)
    do iTmp = 1, rulesDeposition%nAerosolsDepositing
      iSp = rulesDeposition%indAerosolDepositing(iTmp)
      !
      ! Take the aerosol growth into account
      !
      if(rulesDeposition%ifHumidityDependent)then
#ifdef DEBUG
          if (fu_fails(fRelHum /= real_missing, "Missing relaive humidity", "get_settling_motion_species")) return
#endif        
         wetParticle = fu_wet_particle_features(speciesTransport(iSp)%material, fRelHum)
      else
         wetParticle = fu_wet_particle_features(speciesTransport(iSp)%material, rulesDeposition%fDefaultRelHumidity)
      endif

      fWetDiam = fu_massmean_d(speciesTransport(iSp)%mode) * WetParticle%fGrowthFactor

      velocity(iSp) = fu_settling_vel(fWetDiam, WetParticle%fWetParticleDensity, mu_air, fTempr, fPressure)
    end do ! aerosol species

    !Turn it to Pa/s Positive -- down
    if (ifPressure) then
       velocity(1:nSpecies)  = g * fRhoAir *velocity(1:nSpecies)
    endif

  end subroutine get_settling_motion_species


  !*******************************************************************

  subroutine verify_deposition(speciesTransport, nspecies, rulesDeposition)
    !
    ! Check for the consistency of the deposition setup. Used to be
    ! inside resistance-computing subroutines.
    implicit none
    type(silam_species), dimension(:), intent(in) :: speciesTransport
    integer, intent(in) :: nSpecies
    type(Tdeposition_rules), intent(in) :: rulesDeposition

    integer :: i

    do i = 1, nspecies
      if(.not. rulesDeposition%pSpecies(i) == speciesTransport(i))then
        call msg('Species in the deposition list:')
        call report(rulesDeposition%pSpecies(i))
        call msg('Species in the input transport list:')
        call report(speciesTransport(i))
        call set_error('Species are not the same','verify_deposition')
        return
      endif
    end do

    if(rulesDeposition%nAerosolsDepositing < 0 .or. rulesDeposition%nGasesDepositing < 0 .or. &
     & .not. rulesDeposition%nAerosolsDepositing + rulesDeposition%nGasesDepositing <= nSpecies)then
      call msg('Strange number of aerosol or gas depositing species:', &
             & rulesDeposition%nAerosolsDepositing, rulesDeposition%nGasesDepositing)
      call msg('.. or, may be, input transport species:',nSpecies)
      call set_error('Inconsistent number of species','verify_deposition')
    endif
    
  end subroutine verify_deposition

  !************************************************************************************

  subroutine make_deposition_lagr(lpMassTrn, lpDyn, lpStatus, nop, &
                                & mapDryDep, mapWetDep, garbage_mass, &
                                & timestep_sec, &  ! seconds !!
                                & weight_past, &   ! meteo time interpolation
!                                & pMeteo2DispInterpHoriz, pMeteo2DispInterpVert, &
                                & rulesDeposition)
    !
    ! Computes the dry depositon of species in the lagrangian particle.
    ! Idea is that the particle appearing sufficiently close to surface will start
    ! loosing its mass due to dry deposition to the surface. Each particle
    ! then represents the mass in the layer it appeared to be.
    !
    implicit none

    ! Imported parameters
    real, dimension(:,:), intent(inout) :: lpMassTrn, lpDyn
    integer, dimension(:), intent(inout) :: lpStatus
    integer, intent(in) :: nop
    type(Tmass_map),  intent(inout)  :: mapDryDep, mapWetDep
    real, dimension(:,:),  intent(inout)  :: garbage_mass
    real, intent(in) :: timestep_sec, weight_past
    type(Tdeposition_rules), intent(in) :: rulesDeposition

    ! Local variables
    integer :: iSp, iSrc, iP, indexMeteoHoriz, iTmp, nDispLev
    real, dimension(max_species) :: arRs
    type(THorizInterpStruct), pointer :: pMeteo2DispInterpHoriz
    type(TVertInterpStruct), pointer :: pMeteo2DispInterpVert
    logical :: ifVertInterp, ifHorizInterp
    real :: dz1, dM, fVd, fiVd2m

    !
    ! Get the vertical and horizontal interpolation structures
    !
    pMeteo2DispInterpHoriz => fu_horiz_interp_struct(meteo_grid, dispersion_grid, linear, .true.)
    ifHorizInterp = meteo_grid == dispersion_grid

    pMeteo2DispInterpVert => fu_vertical_interp_struct(meteo_vertical, dispersion_vertical, &
                                                     & dispersion_grid, linear, &
                                                     & one_hour, 'meteo_to_dispersion')
    ifVertInterp = fu_cmp_verts_eq(meteo_vertical, dispersion_vertical)

    dz1 = fu_layer_thickness_m(fu_level(dispersion_vertical, 1))
    if(error)return

    nDispLev = fu_NbrOfLevels(dispersion_vertical)
    
    !
    ! Grand loop over particles
    !
    do iP = 1, nop

      if(lpStatus(iP) == int_missing)cycle

      indexMeteoHoriz = fu_grid_index(nx_Meteo, nint(lpDyn(lp_x,iP)),nint(lpDyn(lp_y,iP)), &
                                    & pMeteo2DispInterpHoriz)
      call get_rs(mapDryDep%species, mapDryDep%nSpecies, indexMeteoHoriz, &
                      & weight_past, rulesDeposition, arRs)
      iSrc = mod(lpStatus(iP),100)

      if(nint(lpDyn(lp_z,iP)) == 1)then       ! Dry deposition is done only in the first layer
        do iSp = 1, mapDryDep%nSpecies
          if(arRs(iSp) > -0.5)then
            fVd = fu_get_vd(dz1 * (lpDyn(lp_z,iP)-0.4999),  &    ! Reference height
                          & mapDryDep%species(iSp), &         ! what to deposit
                          & indexMeteoHoriz, weight_past, &   ! position in space and time
                          & arRs(iSp), &                      ! surface resistance 
                          & rulesDeposition%DryDepType, 1, fiVd2m)       ! Depo rules
            dM = lpMassTrn(iSp,iP) * (1. - exp(-fVd * timestep_sec / lpDyn(lp_dz,iP))) 
                                ! ?????? Should it be something like (z+dz)/2
                                ! instead of dz here??? 
                                ! Now point particles would deposit instantly
                                ! FIXME R.
            lpMassTrn(iSp,iP) = lpMassTrn(iSp,iP) - dM
            mapDryDep%arM(iSp, iSrc, 1, nint(lpDyn(lp_x,iP)), nint(lpDyn(lp_y,iP))) = &
                      & mapDryDep%arM(iSp, iSrc, 1, nint(lpDyn(lp_x,iP)), nint(lpDyn(lp_y,iP))) + dM
          endif  ! Rs > 0, i.e. there is dry deposition
        enddo  ! species
      endif  ! first layer

      !
      ! Do the scavenging
      !
      call scavenge_lagr(lpMassTrn(:,iP:iP), &
                       & mapWetDep%arM(:, iSrc:iSrc, 1, &
                                     & nint(lpDyn(lp_x,iP)), nint(lpDyn(lp_y,iP))), &
                       & timestep_sec, weight_past, &
                       & lpDyn(:,iP) ,&
                       & nDispLev, &
                       & mapDryDep%species, mapDryDep%nSpecies, &
                       & rulesDeposition, &
                       & pMeteo2DispInterpHoriz,  pMeteo2DispInterpVert, &
                       & ifVertInterp, ifHorizInterp)

    end do  ! Particle

  end subroutine make_deposition_lagr


  !********************************************************************************************
  !
  !    Wet deposition
  !
  !********************************************************************************************

  !************************************************************************************

  subroutine scavenge_column(mapConc, mapWetDep, arGarbage, tla_step, timestep_sec, &
                        & rulesDeposition, arLowMassThreshold, metdat_col, ix, iy, pwc_col, &
                        & pHorizInterpStruct, pVertInterpStruct, ifHorizInterp, ifVertInterp, met_buf)

    type(TMass_map), intent(inout)  :: mapConc, mapWetDep
    real, dimension(:,:), intent(inout)  :: arGarbage
    type(t_tla_step), intent(inout) :: tla_step
    real, dimension(:), pointer :: arLowMassThreshold
    real, dimension(:,:), pointer, intent(in) :: metdat_col
    real, intent(in) :: timestep_sec !, tcc
    type(Tdeposition_rules), intent(in) :: rulesDeposition
    Integer, intent(in) :: ix, iy
    real, dimension(:), intent(in) :: pwc_col
    type(THorizInterpStruct), intent(in) :: pHorizInterpStruct
    type(TVertInterpStruct), intent(in) :: pVertInterpStruct
    logical, intent(in) :: ifHorizInterp, ifVertInterp
    type(Tfield_buffer), intent(in) :: met_buf

    ! Local variables                                                                                                         
    real :: scav_coef_std, scav_coef,  timestep_sec_abs, fSO2, fS_rest, fMaxScav, precip_rate, &
         & cape_met, tcc, fTempr_1, cell_size_x, cell_size_y, abl_height, landfrac
    integer :: iLev, iSrc, iSpecies, iSpTr, ithread

    real, dimension(max_species) :: rSettling

    logical :: not_adjoint
    real :: pwc_above_top, pwc_above_bottom, pwc_column, cwc_column_layer

    real, dimension(:,:), pointer :: precipContent
    real, dimension(:,:,:), pointer :: precipContent3d
    real, dimension(:,:,:,:), pointer :: TL_store_mass_in_air => null()
    real, dimension(mapConc%nSpecies,mapConc%n3D) :: tl_array

    integer :: istat, iMeteo, num_levs
    logical :: have_tla
    character(len=*), parameter :: sub_name = 'scavenge_column'
    real, dimension(:), pointer ::  pMetPrecTot, pCAPE, pMetTotCloud
    integer, dimension(:), pointer :: mdl_in_q
    type(field_4d_data_ptr), pointer ::  fldCwcabove, fldScavCoefStd
    real, dimension(mapConc%n3D) :: cwc_col
    type (Train_in_cell), dimension(mapConc%n3D) :: r1d

    if (rulesDeposition%scavengingType == scavNoScav) return

    num_levs = mapConc%n3D
    mdl_in_q => met_buf%buffer_quantities

    pMetPrecTot => met_buf%p2d(fu_index(mdl_in_q, total_precipitation_rate_flag))%present%ptr

    iMeteo = fu_grid_index(nx_meteo, ix, iy, pHorizInterpStruct)
    precip_rate = pMetPrecTot(iMeteo)


    if (precip_rate < 1e-30/3600.) return

    iThread = 0
    !$ iThread = OMP_GET_THREAD_NUM() 

    Timestep_sec_abs = Abs(Timestep_sec)
    not_adjoint = timestep_sec > 0

    select case (rulesDeposition%scavengingType)

    case(scavStandard)
      !                                                                                                          
      ! Standard 3D scavenging coefficient. Each material has own scaling applied as a                           
      ! correction factor: material-wise for gases, size-wise for aerosols.                                      
      !                                                                                                          
      ! For aerosols, we will want to have a settling velocity as a scaling parameter                            
      ! with regard to some reference value of a "default" aerosol.                                              
      ! Note that settling is ro*d^2 - proportional. Hence, here we can simply take                              
      ! diameter as the proportionality factor                                                                   
      !                                                                                                          
!      call msg('Standard scavenging...',rulesDeposition%scavengingType)

      fldScavCoefStd => met_buf%p4d(fu_index(mdl_in_q, scavenging_coefficient_flag))

      if(rulesDeposition%nAerosolsDepositing > 0)then
        call get_settling_velocity_species(mapConc%species, mapConc%nSpecies, &
             & 1, 1, 1.0, &      ! indexMeteo, iLev, weight_past                    
             & rulesDeposition, &
             & rSettling, &
             & 0.9)  ! forced relative humidity                                     
        !                                                                                                        
        ! Default settling velocity is for 1 micrometre aerosol with 1000 kg/m3 density,                         
        ! which is about 3.3e-5 m/sec. We take 4.e-5m/s just voluntarily. The sub-cloud scavenging               
        ! is directly proportional to sedimentation velocity but in-cloud is entirely independent                
        ! from it. So, a geometric average is a crude option.                                                    
        ! Note that for small particles, the proportionality to sedimentation disappears                         
        !                                                                                                        
        do iSpecies = 1, mapConc%nSpecies
          if(rSettling(iSpecies) > 4.0e-5)then  ! a formal settling of slightly over-1um particle                
            rSettling(iSpecies) = sqrt(min(1000.0, rSettling(iSpecies) * 2.5e4))
            !             rSettling(iSpecies) = min(1000.0, rSettling(iSpecies) * 2.5e4)  ! sqrt may be too much, try linear  
          else
            rSettling(iSpecies) = 1.0 ! gases or whatever else with small deposition                             
          endif
        end do
      endif  ! if there are aerosols

      cell_size_x = metdat_col(imet_dx_size,1)
      cell_size_y = metdat_col(imet_dy_size,1)

      do iLev = 1, num_levs
        do iSrc = 1, mapConc%nSrc
 
          if(.not. mapConc%ifColumnValid(iSrc,ix,iy))cycle
           
          scav_coef_std = metdat_col(imet_scav, iLev)
          fTempr_1 = metdat_col(imet_temp, 1)

          if(scav_coef_std <= 0.0)cycle
          scav_coef_std = scav_coef_std * timestep_sec_abs

          !                                                                                                     
          ! The scaling with regard to standard 3D coefficient depends on the phase.                            
          ! Gases keep their scaling in the chemical database                                                   
          !
          do iSpecies = 1, rulesDeposition%nGasesDepositing
            iSpTr = rulesDeposition%indGasesDepositing(iSpecies)
            if(mapConc%arM(iSpTr,iSrc,iLev,ix,iy) < arLowMassThreshold(iSpTr) .and. not_adjoint)then
              arGarbage(iSrc, iSpTr) = arGarbage(iSrc, iSpTr) + mapConc%arM(iSpTr,iSrc,iLev,ix,iy)
              mapConc%arM(iSpTr,iSrc,iLev,ix,iy) = 0.0
            else
              if (fTempr_1 > freezing_point_of_water) then
                scav_coef = scav_coef_std * &
                     & fu_scavenging_scale_gas_rain(mapConc%species(iSpTr)%material)
              else
                scav_coef = scav_coef_std * &
                     & fu_scavenging_scale_gas_snow(mapConc%species(iSpTr)%material)
              endif
              ! Negatives for non-depositing stuff                                                              
              if (scav_coef < 0.0) cycle  ! no wet deposition
              ScavAmount(iSpTr,iSrc,iLev,iThread) = mapConc%arM(iSpTr,iSrc,iLev,ix,iy) * (1.-exp(-scav_coef))
            endif
          end do  !  gases                                                                                      

          !                                                                                                     
          ! Aerosols scaling can be considered either via the settling velocity or via                          
          ! the material properties, just like gases. The default option is probably the size.                  
          !
          do iSpecies = 1, rulesDeposition%nAerosolsDepositing
            iSpTr = rulesDeposition%indAerosolDepositing(iSpecies)
            if(mapConc%arM(iSpTr,iSrc,iLev,ix,iy) < arLowMassThreshold(iSpTr) .and. not_adjoint)then
              arGarbage(iSrc, iSpTr) = arGarbage(iSrc, iSpTr) + mapConc%arM(iSpTr,iSrc,iLev,ix,iy)
              mapConc%arM(iSpTr,iSrc,iLev,ix,iy) = 0.0
            else
              ScavAmount(iSpTr,iSrc,iLev,iThread) = mapConc%arM(iSpTr,iSrc,iLev,ix,iy) * &
                   & (1.-exp(-scav_coef_std * rSettling(iSpTr)))
              !                mapConc%arM(iSpTr,iSrc,iLev,ix,iy) = mapConc%arM(iSpTr,iSrc,iLev,ix,iy) - ScavAmount             
              !                mapWetDep%arM(iSpTr,iSrc,1,ix,iy) = mapWetDep%arM(iSpTr,iSrc,1,ix,iy) + ScavAmount               
            endif
          end do  !  aerosols
        end do   ! iSrc                                                                                         
      enddo   ! iLev
      !                                                                                                        
      ! Saturation of the sulphur deposition and other non-linear effects can be evaluated now.                
      !                                                                                                        
      if(rulesDeposition%ifSulphurSaturation)then
        !                                                                                                      
        ! Sum-up the total amount and find out the excess over the max scavengeable value                      
        !                                                                                                      
        fSO2 = 0.
        fS_rest = 0.
        do iLev = 1, mapConc%n3D
          do iSrc = 1, mapConc%nSrc
            do iSpecies = 1, nSulphurSpecies
              if(indSulphurSpecies(iSpecies) == indSO2)then
                fSO2 = fSO2 + ScavAmount(indSulphurSpecies(iSpecies),iSrc,iLev,iThread)
              else
                fS_rest = fS_rest + ScavAmount(indSulphurSpecies(iSpecies),iSrc,iLev,iThread)
              endif
            end do
          end do  ! iSrc                                                                                       
        end do  ! iLev                                                                                         
        if (fTempr_1  > freezing_point_of_water) then    ! Rain                                                  
          fMaxScav = precip_rate * cell_size_x * cell_size_y * timestep_sec_abs * &
               & rulesDeposition%fSulphurSaturation
        else                                   ! Snow                                                          
          fMaxScav = precip_rate * cell_size_x * cell_size_y * timestep_sec_abs * &
               & rulesDeposition%fSulphurSaturation / 1.5
        endif
        !                                                                                                      
        ! We assume that SO2 cannot go into droplet if there is enough sulfur. For                             
        ! sulphates that are already there, nothing changes - they can even go higher                          
        ! than the saturation level                                                                            
        !                                                                                                      
        if(fSO2 + fS_rest > fMaxScav)then
          call msg('SOx saturation (ix,iy), fSO2, fS_rest:' + fu_str(ix) + ',' + fu_str(iy),fSO2, fS_rest)
          fMaxScav = max(fMaxScav - fS_rest, 0.) / fSO2
          do iLev = 1, mapConc%n3D
            do iSrc = 1, mapConc%nSrc
              ScavAmount(indSO2,iSrc,iLev,iThread) = ScavAmount(indSO2,iSrc,iLev,iThread) * fMaxScav
            end do  ! iSrc                                                                                     
          end do  ! iLev                                                                                       
        end if
      endif   ! ifSulphur saturation                                                                           
      !                                                                                                        
      ! Having all magic done with the scavenged amount, introduce it into the main arrays                     
      !                                                                                                        
      do iLev = 1, mapConc%n3D
        do iSrc = 1, mapConc%nSrc
          do iSpTr = 1, mapConc%nSpecies
            mapConc%arM(iSpTr,iSrc,iLev,ix,iy) = mapConc%arM(iSpTr,iSrc,iLev,ix,iy) - ScavAmount(iSpTr,iSrc,iLev,iThread)
            mapWetDep%arM(iSpTr,iSrc,1,ix,iy) = mapWetDep%arM(iSpTr,iSrc,1,ix,iy) + ScavAmount(iSpTr,iSrc,iLev,iThread)
            ScavAmount(iSpTr,iSrc,iLev,iThread) = 0.0
          end do
        end do
      end do

      
      !case(scavNoScav) 
      !Should crash if we got here
   !   call msg('Dummy scavenging -- no scavenging at all!')
   !   return
      
    case(scav2011fc, scav2011, scav2018, scav2018entr, scav2020)    
      ! Should be checked by init_standard_deposition                                                                                                                              
      !       call  check_deposition_rules(rulesDeposition, mapConc%species)                                                                                                       
      !                                                                                                                                                                    
      ! If adjoint is needed then TLA structure is defined and we need to pick or create                                                                                   
      ! the 4D array for TL storage. Its dimensions correspond to massMap with eliminated nSrc=1                                                                           
      !                                                                                                                                                                    
      if (defined(tla_step)) then
        call msg('Have tla step, SCAV20XX')
        have_tla = .true.
        TL_store_mass_in_air => fu_get_tla_point(tla_step, rulesDeposition%scavengingType)
        if(fu_fails(mapConc%n3d == 1,'TLA allow only one source in mass map', sub_name))return
      else
        nullify(TL_store_mass_in_air)
        have_tla = .false.
      end if
      
      if (.not. any(mapConc%ifColumnValid(:,ix, iy)))then
        if (have_tla .and. not_adjoint) TL_store_mass_in_air(:,:,ix,iy) = 0.
        return
      endif

      ! Store TLA on forward run

      if ( not_adjoint) then
        if (have_tla) TL_store_mass_in_air(:,:,ix,iy) = mapConc%arM(:, 1, :, ix, iy)
      else !! Adjoint
        if (have_tla) then
           tl_array(:,:) = TL_store_mass_in_air(:,:,ix,iy)
        else
           tl_array(:,:)  = 0. !! No TLA, pretend they were all zeroes
        endif
      end if

                  
      if (rulesDeposition%scavengingType == scav2020) then
        pwc_column = sum(pwc_col(1:num_levs))
      else
        pwc_column = real_missing
      end if

      ! Contents of precipitation at the bottom of each layer. Essentialy, reuse ScavAmount array                                                                          
      !                                                                                                                                                                    
      if(timestep_sec > 0.)then
        ! forward: no z-dimension needed, cumulative is good enough                                                                                                        
        precipContent => scavAmount(1:mapConc%nSpecies,1:mapConc%nSrc,1,ithread)
      else
        ! adjoint: need all levels                                                                                                                                         
        precipContent3d => scavAmount(1:mapConc%nSpecies, 1:mapConc%nSrc, 1:mapConc%n3D, ithread)
      endif
      
    
      tcc = metdat_col(imet_tcc,1)
      if (rulesDeposition%max_scav_rate_depends_on == cape .or. &
           & rulesDeposition%max_scav_rate_depends_on == cape_and_wind) then
        cape_met = metdat_col(imet_cape,1)
        abl_height = metdat_col(imet_abl,1)
      end if
      if (rulesDeposition%max_scav_rate_depends_on == cape_and_wind) then
        landfrac = metdat_col(imet_landfrac,1)
      end if
      if (associated(imet_cwc3d)) then ! Enable fakecloud branch
        cwc_col(1:num_levs) = metdat_col(imet_cwc3d,1:num_levs)
      else
        cwc_col(1:num_levs) = 0.
      endif
      
      do iLev = num_levs,1,-1
          if (rulesDeposition%scavengingType == scav2020) then
            if (iLev < num_levs) then
              pwc_above_top = sum(pwc_col(iLev+1:num_levs))
            else
              pwc_above_top = 0.0
            end if

            if (iLev <= num_levs) then
               pwc_above_bottom = sum(pwc_col(iLev:num_levs))
            else
               pwc_above_bottom = 0.0
            end if
          else
            pwc_above_top = real_missing
            pwc_above_bottom = real_missing
          end if

          if (iLev < num_levs) then
            cwc_column_layer = cwc_col(iLev) - cwc_col(iLev+1)
          else
            cwc_column_layer = cwc_col(iLev)
          end if

          cwc_column_layer = max(0.0, cwc_column_layer)
          !
          call make_rain_in_cell_col(iLev, num_levs, r1d(iLev),  &
               & rulesDeposition%scavengingType, rulesDeposition%max_scav_rate_depends_on, &
               & metdat_col, pwc_col(iLev), pwc_column, pwc_above_top, pwc_above_bottom, &
               & cwc_column_layer, cwc_col(1), cwc_col(iLev), landfrac, precip_rate, tcc, cape_met, &
               & abl_height, rulesDeposition%max_scav_rate_wind_scaling, &
               & rulesDeposition%max_scav_rate_cape_scaling)

          r1d(iLev)%ix = ix
          r1d(iLev)%iy = iy
          r1d(iLev)%iLev = iLev
      enddo

      if(not_adjoint)then
                
        precipContent(1:mapConc%nSpecies, 1:mapConc%nSrc) = 0. !Reset
        do iLev = num_levs,1,-1
          if (.not. r1d(iLev)%ifRainyCell) cycle  ! no rain at this level
          
          call scavenge_puff_2011(mapConc%arM(:,:,iLev,ix,iy), &    ! stuff in the air                                                                                               
               & precipContent(:,:), & ! stuff in precipitation                                                                                                     
               & arLowMassThreshold, &
               & mapConc%species,      &
               & mapConc%nSpecies,&
               & mapConc%nSrc,&
               & timestep_sec,&
               & rulesDeposition, &
               & r1d(iLev), r1d(iLev)%cell_volume, r1d(iLev)%cell_zsize)
        end do
        mapWetDep%arM(1:mapConc%nSpecies,1:mapConc%nSrc,1,ix,iy) = &
             & mapWetDep%arM(1:mapConc%nSpecies,1:mapConc%nSrc,1,ix,iy) + &
             & precipContent(1:mapConc%nSpecies,1:mapConc%nSrc)

      else
        !                                                                                                                                      
        ! Adjoint                                                                                                                                                                     
        ! Now, call the adjoint scavenging keeping in mind that we need the whole column                                                                                              
        ! There will be several passes along the vertical handled inside                                                                                                              
        !                                                                                                                                                                             
        call scavenge_column_2011_adjoint(mapConc%arM(:,:,:,ix,iy), & ! sensitivity in the air, column                                                                                
             & precipContent3d(:,:,:), &   ! sensitivity mass in precip, column                                                                            
             & arLowMassThreshold, &
             & tl_array, &     ! TL- stored mass in the air                                                                         
             & mapConc%species,  &
             & mapConc%nSpecies, mapConc%nSrc, mapConc%n3D, &
             & - timestep_sec, &   ! invert the sign to have it positive                                                                                   
             & rulesDeposition, &
             & r1d)
        ! Finally deposit rain contents....      
        do iLev = 1, mapConc%n3D
          mapWetDep%arM(1:mapConc%nSpecies,1:mapConc%nSrc,1,ix,iy) = &
               & mapWetDep%arM(1:mapConc%nSpecies,1:mapConc%nSrc,1,ix,iy) + &
               & precipContent3d(1:mapConc%nSpecies,1:mapConc%nSrc, iLev)
        end do
      endif  ! forward or adjoint
      
    case default
      call msg('Unknown scavenging type',rulesDeposition%scavengingType)
      call set_error('Unknown scavenging type', sub_name)
    end select ! type of scavenging

  end subroutine scavenge_column


  !*************************************************************************************                      
                        
  subroutine report_rain(r)
    type (Train_in_cell), intent(in) :: r
    type(silam_sp) :: l
    
    l%sp => fu_work_string()
    if(error)return
    
    call msg("***Reporting rain in cell")
    WRITE (l%sp, '(A,L1,A,L1,A,L1,A,L1)')  "ifRainyCell = ", r%ifRainyCell, & 
                        &", fRainSfc = ", r%ifRainSfc, &
                        &", ifRealCloud = ",r%ifRealCloud, &
                        &", ifValid = ", r%ifValid
    call msg(l%sp)

    WRITE (l%sp, '(A,F5.3,A,F7.0,A,F5.1,A,F5.2,A,F5.2)') &
                        &"fCloudCover = ",r%fCloudCover , &
                        &", fPressure = ", r%fPressure, &
                        &", fTemperature = ", r%fTemperature, &
                        &", fRainRateSfc(mm/h) = ", r%fRainRateSfc*3600,&
                        &", cwcColumn(mm) = ", r%cwcColumn
    call msg(l%sp)

    WRITE (l%sp, '(A,F10.3,A,F10.3,A,F10.3,A,F10.3)') &
                        &"rain_rate_in(mm/h) = ", r%rain_rate_in*3600, &
                        &", deltarain (mm/h) = ", r%deltarain*3600, &
                        &", cwcColumnLayer(mm) = ", r%cwcColumnLayer, &
                        &", pwcColumnLayer(mm) = ", r%pwcColumnLayer 
    call msg(l%sp)

    WRITE (l%sp, '(A,F3.1,A,F5.3,A,F5.0,A,F5.2,A,F10.2)') &
                        &", is_liquid = ", r%is_liquid, &
                        &", drop_size(mm) = ", r%drop_size*1000, &
                        &", drop_mass(mg) = ", r%drop_mass*1e9, &
                        &", drop_vel = ", r%drop_vel, &
                        &", drop_Re = ", r%drop_Re
    call msg(l%sp)

    WRITE (l%sp, '(A,F5.0,A,F10.0)') &
                        &", cell_zsize = ", r%cell_zsize, &
                        &", cell_volume(km3) = ", r%cell_volume*1e-9
    call msg(l%sp)
    call free_work_array(l%sp)
    
  end subroutine report_rain

!subroutine check_deposition_rules(rulesDeposition, species)
!    type(Tdeposition_rules), intent(in) :: rulesDeposition
!    type(silam_species), dimension(:), pointer, intent(in) :: species
!    type(Tgas_deposition_param), pointer :: DepData
!    integer :: iSpecies,iSpTr
!       ! check if needed parameters available
!
!       do iSpecies = 1,rulesDeposition%nGasesDepositing
!              iSpTr = rulesDeposition%indGasesDepositing(iSpecies)
!              if(fu_scavenging_scale_gas_snow(species(iSpTr)%material) <= 0.)cycle
!                ! Should be some more inteligent way to figure out that the
!                ! stuff does not deposit....
!
!             ! call msg ("Gas Depositing " +  fu_name(species(iSpTr)%material))
!              DepData = fu_gas_deposition_param(species(iSpTr)%material)
!              if    (DepData%Henry_const_298K .eps. real_missing)    then
!                call msg ("Henry constant is missing for:" + &
!                               & fu_name(species(iSpTr)%material))
!                call set_error("Henry constant missing", "check_deposition_rules")
!              endif
!              if  (DepData%Henry_const_T_dep .eps. real_missing) then
!                call msg ("Assuming zero Henry_const_T_dep for:" + &
!                               & fu_name(species(iSpTr)%material))
!                      DepData%Henry_const_T_dep = 0.
!              endif
!              if  (DepData%Wesely_f0 .eps. real_missing) then
!                call msg ("Assuming zero Wesely_f0 for:" + &
!                               & fu_name(species(iSpTr)%material))
!                      DepData%Wesely_f0 = 0.
!              endif
!       enddo
!end subroutine check_deposition_rules

  
  !*****************************************************************************************************

  
  subroutine scavenge_puff_2011(mass_in_air, &    ! stuff in the air
                              & mass_in_water, &  ! stuff in precipitation
                              & arLowMassThreshold, &
                              & species,      &   ! array of 
                              & Nspecies,&
                              & Nsources,&
                              & timestep_sec,&
                              & rulesDeposition, &
                              & r, cell_volume, cell_zsize, &  ! volume and zsize differ in Lagrange environ
                              & coefs)  ! internal coefficients, if requested
    !
    ! Scavenges stuff from a puff into precipitation
    ! caller has to handle the stuff from mass_in_water.
    ! In Eulerian case -- puff is a dispersion cell
    ! In Lagrangian case -- puff is a particle
    ! R
    ! Imported parameters
    type (Train_in_cell), intent(in) :: r
    type(Tdeposition_rules), intent(in) :: rulesDeposition
    real, dimension(:), intent(in) :: arLowMassThreshold
!        real, dimension(:,:), pointer, intent(in) :: arGarbage
    type(silam_species), dimension(:), intent(in) :: species
    real, dimension(:,:), intent(inout) :: mass_in_air, mass_in_water
    integer , intent (in) :: Nspecies, Nsources
    real, intent(in) :: timestep_sec, cell_volume, cell_zsize
    type(Tscav_coefs), intent(out), optional :: coefs  ! if need to store internal values
    
    ! Local variables
    integer :: iSpTr,iSrc ! indices for mass_in_air, precipContent,  arLowMassThreshold, arGarbage
    integer :: iSp, iSpeciesSA
    real :: efficiency, acid_mol
    real :: equilibrium_aq_SIV, fB, fC, fD ! stuff for sulphur equilibrium
    real, parameter :: H0_SO2    = 1.2e-5  !True Henry const for SO2 @ 298.15 K  [Mol/kg/Pa]
   !!!    real, parameter :: fKs10_SO2 = 1.3e-2  !Dissociation constant for H2SO3   @ 298.15 K [Mol/kg]
                                   ! I have not found any reasonable data on temperature-dependence. R.
    real, parameter :: H0_SO2_Tconst = 3000 ! [K] https://webbook.nist.gov/cgi/cbook.cgi?ID=C7446095&Mask=10
    real ::  H_SO2 ! Henry const for SO2
    real ::  fKs1_SO2 ! First dissociation of H2SO3
    real :: precip_rate, scav_coef,scav_coef_ic,scav_coef_sc, tauf_over_RdCd, &
          & total_SIV_in_scav_zone, scav_zone_volume, total_water_in_scav_zone, &
          & unscav_mass_in_scav_zone, & ! Mass in equilibrium with incoming precipitation
          & water_capacitance_factor, scav_amt, scav_coef_tau, A, fTmp

    scav_zone_volume = cell_volume * r%fCloudCover
    precip_rate = r%rain_rate_in /  r%fCloudCover
    tauf_over_RdCd = 0.0
    
    !
    ! Sources are considered indeoendently from each other - for the case if we run many
    ! runs with varying source configuration
    !    
    do iSrc = 1, nSources
      do iSp = 1, nSpecies
        !
        ! Note that SO2 has to be scavenged last because it is the only one that currently depends
        ! on the rain acidity, i.e. it is affected by other acidic and alcaline species
        !
        iSptr=iSp
        if (iSp == indSO2) iSpTr = nSpecies ! scavenge instead of SO2
        if ((indSO2 > 0) .and. (iSp == nSpecies)) iSpTr = indSO2   ! Finally SO2

        !
        ! Below, in case of small mass etc, we skip the scavenging. Prepare for this
        ! by setting TL coefs to zero
        !
        if(present(coefs))then
          coefs%scav_coef_ic(iSptr) = 0
          coefs%scav_coef_sc(iSptr) = 0
          coefs%scav_coef_tau(iSptr) = 0
          coefs%A(iSptr) = 0
          coefs%tauf_over_RdCd(iSptr) = 0
          if(iSptr== indSO2)then
            coefs%fB = 0
            coefs%fC = 0
            coefs%Ceq = 0
            coefs%dLsc_dWc = 0
            coefs%dLic_dWc = 0
            coefs%dCeq_dMs = 0
            coefs%dWC_dMs = 0
            coefs%dCeq_dMa = 0
            coefs%dWC_dMa = 0
          endif  ! if SO2
        endif   ! if internal coefs are requested

        ! puting directly to garbage here is risky....
        if(mass_in_air(iSpTr,iSrc) +  mass_in_water(iSpTr,iSrc) < arLowMassThreshold(iSpTr))cycle
        !
        ! Get material-dependent properties. 
        !
        call material_properties(species(iSpTr), r, water_capacitance_factor, efficiency)

        if(error)return

        if (water_capacitance_factor <0) cycle  !! not scavenged


        if(iSpTr == indSO2) then       !adjust water_capacitance_factor for SO2
          !
          ! Capacitance depends on acids and SO2 itself...
          ! Note that SO2 should be the last species to scavenge.
          ! here it is defined formally as above mentioned ratio
          ! for total rain amount + clouds and total SO2 (air+water)
          !
          total_water_in_scav_zone = (r%rain_rate_in * timestep_sec + r%cwcColumnLayer) * &  ! kg/m2
                                   & scav_zone_volume / cell_zsize       ! cell area, m2
          ! Total amount of SO2 in air and, in a dissolved form, in water
          total_SIV_in_scav_zone = mass_in_air(iSpTr,iSrc)*r%fCloudCover + mass_in_water(iSpTr,iSrc)
          ! Calculate very roughly amount of strong acids and alkalines in rain
          acid_mol=0.
          do iSpeciesSA = 1, nStrongAcidSpecies
            acid_mol = acid_mol + mass_in_water(indStrongAcidSpecies(iSpeciesSA),iSrc)  !acid is in moles
          enddo
          ! Subtract alcalines....
          if (indNH3 > 0) acid_mol = acid_mol - mass_in_water(indNH3,iSrc)
          if (.not. total_SIV_in_scav_zone > 0.) then
            if (.not. total_SIV_in_scav_zone >= 0.) then
!$OMP CRITICAL(WC_SO2_1)
              call msg(".not. total_SIV_in_scav_zone > 0., iSpTr, iSrc",iSpTr,iSrc)
              call msg("mass_in_air(iSpTr,iSrc)",mass_in_air(iSpTr,iSrc))
              call msg("mass_in_water(iSpTr,iSrc)",mass_in_water(iSpTr,iSrc))
              call msg("total_water_in_scav_zone", total_water_in_scav_zone)
              
              call msg("r%cwcColumnLayer", r%cwcColumnLayer)
              call msg("cell_zsize", cell_zsize)
              call msg("scav_zone_volume", scav_zone_volume)
              call msg("r%fCloudCover", r%fCloudCover)
              call msg("r%rain_rate_in", r%rain_rate_in)
              call msg("precip_rate", precip_rate)
!$OMP END CRITICAL(WC_SO2_1)
              cycle
            endif
          else
            if (total_water_in_scav_zone > 0.) then
              acid_mol = acid_mol / total_water_in_scav_zone ! Mole/kg
            else
!$OMP CRITICAL(WC_SO2_1)
              call msg("Rain with no water!")
              call msg(".not. total_water_in_scav_zone > 0., iSpTr, iSrc",iSpTr,iSrc)
              call msg("mass_in_air(iSpTr,iSrc)",mass_in_air(iSpTr,iSrc))
              call msg("mass_in_water(iSpTr,iSrc)",mass_in_water(iSpTr,iSrc))
              call msg("total_water_in_scav_zone", total_water_in_scav_zone)
              call msg("acid_mol (moles)", acid_mol)
              
              call msg("r%cwcColumnLayer", r%cwcColumnLayer)
              call msg("cell_zsize", cell_zsize)
              call msg("scav_zone_volume", scav_zone_volume)
              call msg("r%fCloudCover", r%fCloudCover)
              call msg("r%rain_rate_in", r%rain_rate_in)
              call msg("precip_rate", precip_rate)
!$OMP END CRITICAL(WC_SO2_1)
              cycle
            endif
          endif  ! total SIV in scavenging zone >< 0
 
          if (acid_mol < 2.5e-6) acid_mol = 2.5e-6 ! pH5.6

          ! Solubility 

          ! Equilibrium molar fraction of S(IV)  [mol/kg]
          ! stuff_in_... supposed to be moles...
          ! quadratic equation for equilibrium concentration in water
          ! [mol/kg]  fA [S(IV)]^2 + fB [S(IV)] + C = 0
          ! See notebook R1 p. 37.
          ! fKs1_SO2 - Dissociation constant for HSO3 [mol/kg]
          !  fA = -1./(fKs1_SO2 * water_capacitance_factor) !True Henry constant needed here
          !
          
          ! https://webbook.nist.gov/cgi/cbook.cgi?ID=C7446095&Mask=10#Solubility
          !
          H_SO2 = H0_SO2 * exp(H0_SO2_Tconst*(1./r%fTemperature - 1./298.5)) ![Mol/kg]

          ! Sulfur‐Dioxide/Water Equilibria Between 0° and 50°C. An Examination of Data at Low Concentrations
          ! Howard G. Maahs (1982) https://doi.org/10.1029/GM026p0187 
          ! Cited after https://www.slideserve.com/asasia/b-krasovitov-t-elperin-and-a-fominykh-department-of-mechanical-engineering
          !  10**(853/T - 4.74) is there
          !
          fKs1_SO2 = exp(1964.1/r%fTemperature - 10.91) ![mol/kg]  

          fTmp = (fKs1_SO2 * gas_constant_uni*r%fTemperature * H_SO2)/ scav_zone_volume
          fB =  acid_mol + fTmp * total_water_in_scav_zone
          fC = - fTmp * total_SIV_in_scav_zone
          if (fB*fB > 10000*abs(fC)  ) then ! No need for quadratic equation, Taylor series sufficient
            equilibrium_aq_SIV = - fC/fB - fC*fC / (fB*fB*fB)
          else
            fD = fB*fB - 4*fC
            equilibrium_aq_SIV =  0.5*(-fB + sqrt(fD))
          endif
          !
          ! In case of no strong acids (acid_mol=0) in the droplet and quadratic equation reduced to 
          ! linear approximation, capasitance will be infinity. Should prevent this from happening
          !
          fTmp =  1 -  equilibrium_aq_SIV*total_water_in_scav_zone / total_SIV_in_scav_zone 
                  !! fraction of S_IV that wants to stay in air
          if (abs(fTmp) < 1e-3 ) fTmp = 1e-3  !!Some numeric can be

          water_capacitance_factor =  scav_zone_volume / total_water_in_scav_zone / fTmp

          if (.not. water_capacitance_factor >= 0) then 
!$OMP CRITICAL(WC_SO2)
            fD = fB*fB - 4*fC
            call msg(".not. water_capacitance_factor > 0., iSpTr, iSrc",iSpTr,iSrc)
            call msg("mass_in_air(iSpTr,iSrc)",mass_in_air(iSpTr,iSrc))
            call msg("mass_in_water(iSpTr,iSrc)",mass_in_water(iSpTr,iSrc))
            call msg("total_water_in_scav_zone", total_water_in_scav_zone)
            
            call msg("r%cwcColumnLayer", r%cwcColumnLayer)
            call msg("cell_zsize", cell_zsize)
            call msg("scav_zone_volume", scav_zone_volume)
            call msg("r%fCloudCover", r%fCloudCover)
            call msg("r%rain_rate_in", r%rain_rate_in)
            call msg("precip_rate", precip_rate)

            call msg("*** Warning! Negative SO2 capacitance factor:",water_capacitance_factor)
            call msg(" fraction of S_IV that wants to stay in air", fTmp)
            call msg("fB, fC, fD, fB^2, 4*fC, -fC/fB, quadr_slv", &
                                & (/fB, fC, fD, fB*fB, 4.*fC, -fC/fB, 0.5*(-fB + sqrt(fD))/))
            call msg('fKs1_SO2, gas_constant_uni, r%fTemperature, H_SO2',(/fKs1_SO2, gas_constant_uni, r%fTemperature, H_SO2/))
            call msg("total_water_in_scav_zone, total_SIV_in_scav_zone, scav_zone_volume, cell_zsize, acid_mol", &
                   & (/total_water_in_scav_zone, total_SIV_in_scav_zone, scav_zone_volume, cell_zsize, acid_mol/))
            call msg('fC, equilibrium_aq_SIV', fC, equilibrium_aq_SIV)
            call msg('their ratios:', fC/equilibrium_aq_SIV, total_water_in_scav_zone/scav_zone_volume)
            water_capacitance_factor = 1e-3 !something small 
            call msg("Rain report:")
            call report_rain(r)
            call set_error("WC_SO2","scavenge_puff_2011")

!$OMP END CRITICAL(WC_SO2)
          endif

        endif ! SO2

        if (water_capacitance_factor < 1e-10) cycle       ! less than 1l of gas in 1 l of water

        if (precip_rate > 1e-5) then 
          !See notebook R1 p35. bottom
          ! Flight_time / Droplet_equilibration_time
          tauf_over_RdCd =  1.5* cell_zsize * efficiency/(1000.* r%drop_size *water_capacitance_factor) 
                                                       ! water density 
          ! 1.-exp(-tauf_over_RdCd) -- rain-air equilibration factor 
          if (tauf_over_RdCd > 0.01) then  ! drop approaches saturation, exp is needed
            scav_coef_sc = (1.0 - exp(-tauf_over_RdCd)) * &
                                            & precip_rate * water_capacitance_factor / cell_zsize
          elseif  (tauf_over_RdCd >= 0.) then ! drop is far from saturation -- linear scavenging
            scav_coef_sc = 3. / (2.*1000) * precip_rate / r%drop_size * efficiency 
          else ! too bad....
            call msg_warning("Non-positive scavenging coefficient...", "scavenging")
            call msg("tauf_over_RdCd", tauf_over_RdCd)
            call msg("Drop size", r%drop_size)
            call msg("Re", r%drop_Re)
            call msg("water path", cell_zsize)
            call msg("efficiency", efficiency )
            call msg("water_capacitance_factor", water_capacitance_factor)
            call set_error("Non-positive scavenging coefficient...", "scavenging")
          endif
        else
          scav_coef_sc = 0.
        endif ! precip rate

        ! in-cloud -- fraction of rain increase. (p.38)
!        if (r%deltarain < 0.001 * r%rain_rate_in) then ! no in-cloud
!          scav_coef_ic = 0.
!        else
          scav_coef_ic = r%deltarain * water_capacitance_factor /  &
                           &  (r%cwcColumnLayer * water_capacitance_factor + cell_zsize)
!        endif

        if (scav_coef_ic < 0.0) then
          call msg_warning("Negative in-cloud  scavenging coefficient...", "scavenging")
          call msg("scav_coef_ic", scav_coef_ic)
          call msg("scav_coef_sc", scav_coef_sc)
          call msg("cell_zsize", cell_zsize)
          call msg("r%deltarain", r%deltarain)
          call msg("r%pwcColumn", r%pwcColumn)
          call msg("r%pwcColumnLayer", r%pwcColumnLayer)
          call msg("r%cwcColumn", r%cwcColumn)
          call msg("r%cwcColumnLayer", r%cwcColumnLayer)
          call msg("Water capacitance_factor", water_capacitance_factor)
          call msg("r%rain_rate_in", r%rain_rate_in)
          call set_error("Negative in-cloud  scavenging coefficient...", "scavenging")
        end if

        scav_coef = scav_coef_ic + scav_coef_sc

        if ((scav_coef - scav_coef) /= 0.) then
          call msg_warning("NaN scavenging coefficient...", "scavenging")
          call msg("**Scav** NaN scavenging_coef...", scav_coef)
          call msg("In_cloud", scav_coef_ic)
          call msg("Sub_cloud", scav_coef_sc)
          call msg("Water path", cell_zsize)
          call msg("Water capacitance_factor", water_capacitance_factor)
          call msg("tauf_over_RdCd", tauf_over_RdCd )
          call msg("Precipitation rate", r%rain_rate_in*3600.)
          call msg("CwC rate", r%cwcColumnLayer )
          call msg("Water timestep equiv column", r%rain_rate_in*water_capacitance_factor *timestep_sec)
          call set_error("NaN scavenging coefficient...", "scavenging")
        endif

!        if (scav_coef < 1e-6) cycle !Scav. time more than ~10 days
        !
        ! Stuff can not be scvenged faster than cell turnover rate (rate=wind/cellsize)
        !
        select case(rulesDeposition%scavengingType)
          case(scav2018, scav2020)
            scav_coef = scav_coef * r%max_scav_rate / (scav_coef + r%max_scav_rate)
          case(scav2018entr)
            ! logistic function, seems to be a bit sharper transition
            scav_coef = 2.0 * r%max_scav_rate * (1. / (1. + exp(-scav_coef/r%max_scav_rate)) - 0.5)
          case default
        end select

        scav_coef_tau = scav_coef * timestep_sec
        if (scav_coef_tau > 0.01) scav_coef_tau = 1. - exp(-scav_coef_tau) ! equilibration factor to incoming precipitation

        ! scavenge from all sources       
        ! p36 
        ! in-air mass that would equilibrate the precipitation content
        !
        if (scav_coef_sc > 0) then
          A = cell_zsize * scav_coef_sc / (scav_coef * precip_rate * timestep_sec * water_capacitance_factor)
          unscav_mass_in_scav_zone = mass_in_water(iSpTr, iSrc) * A
        else 
          A = 0
          unscav_mass_in_scav_zone = 0
        endif

        ! stuff to be transferred  air -> water (happens only in cloudy part) 
        ! Yes. It makes the whole thing timestep-dependent....
        !
        scav_amt = ( mass_in_air(iSpTr,iSrc) * r%fCloudCover - unscav_mass_in_scav_zone) * scav_coef_tau
        
        !
        ! Deposit!!!

        mass_in_air(iSpTr,iSrc) =  mass_in_air(iSpTr,iSrc) - scav_amt
        mass_in_water(iSpTr,iSrc) =  mass_in_water(iSpTr,iSrc) + scav_amt

        ! Check if it is still reasonable. Should trigger NaNs alert as well..
        !
        if ( .not. ((mass_in_water(iSpTr,iSrc) >= 0.) .and. &
                                              &  (mass_in_air(iSpTr,iSrc) >= 0.))  ) then
          !
          ! For forward runs, this may be an alert
          !
          ! Could be just numerics
          if( (mass_in_water(iSpTr,iSrc) < 0.) .and. &
                                    & (mass_in_water(iSpTr,iSrc) >  -1e-4*abs(scav_amt)) ) then
            mass_in_water(iSpTr,iSrc) = 0.
            mass_in_air(iSpTr,iSrc) = mass_in_air(iSpTr,iSrc) + mass_in_water(iSpTr,iSrc)
          elseif (  (mass_in_air(iSpTr,iSrc) < 0.) .and. &
                                     & (mass_in_air(iSpTr,iSrc) > -1e-4*abs(scav_amt))) then
            mass_in_air(iSpTr,iSrc) = 0.
            mass_in_water(iSpTr,iSrc) = mass_in_air(iSpTr,iSrc) + mass_in_water(iSpTr,iSrc)
          else ! too bad...
            call set_error("**Scav forward**  Wrong mass scavenged..", "scavenge_puff_2011")
          endif
          !
          ! More reporting in case of problem
          if(error)then
            call msg("**Scavenging subroutine** Wrong mass scavenged..", scav_amt)
            call report_rain(r)
            call msg("ScavFactor", scav_coef_tau)
            call msg("ScavCoef incl, subcl:", scav_coef_ic, scav_coef_sc)
            call msg("water cap. factor [m3/kg]",water_capacitance_factor)
            call msg("mass in air in scav.zone   before", (mass_in_air(iSpTr,iSrc) + scav_amt)*r%fCloudCover)
            call msg("mass to stay in scav.zone   before", unscav_mass_in_scav_zone)
            call msg("mass in water before", (mass_in_water(iSpTr,iSrc) - scav_amt))
            call msg("mass in air in scav.zone equilibrium", unscav_mass_in_scav_zone)
            call msg("mass in air scav. zone   after", mass_in_air(iSpTr,iSrc)*r%fCloudCover)
            call msg("mass in water after", mass_in_water(iSpTr,iSrc) )
          endif
        endif  ! potentially mass problem
 
        !
        ! If the internal coefficients are needed, make them available
        ! See notebook mas 10a, pp.27-36. From these formulas, we get the needed variables
        !
        if(present(coefs))then
          coefs%scav_coef_ic(iSpTr) = scav_coef_ic
          coefs%scav_coef_sc(iSpTr) = scav_coef_sc
          coefs%scav_coef_tau(iSpTr) = scav_coef_tau 
          coefs%A(iSpTr) = A               ! factor for non-scavenging mass, notebook mas 10a p.27-28
          coefs%tauf_over_RdCd(iSpTr) = tauf_over_RdCd   ! rest is needed for SO2, n. mas 10a, pp.34-36
          if(iSpTr == indSO2)then
            coefs%fB = fB  ! coef from quadratic equation
            coefs%fC = fC  ! coef from quadratic equation
            coefs%Ceq = equilibrium_aq_SIV  ! solution of quadratic equation
            coefs%dLsc_dWc = precip_rate / cell_zsize * &           ! d lambda_subcloud / d WC, p.35
                           & (1. - exp(-scav_coef_tau) * (1. + scav_coef_tau))
            coefs%dLic_dWc = r%deltarain * cell_zsize /  &          ! d lambda_incloud / d WC, p.35
                           & (r%cwcColumnLayer * water_capacitance_factor + cell_zsize) ** 2
            coefs%dCeq_dMs = -1. / (scav_zone_volume * (2 * coefs%Ceq + fB)) ! d Ceq / d M_SO2, p.37
            coefs%dWC_dMs = scav_zone_volume * &                  ! d Wc / d M_SO2, p.37
                          & (total_SIV_in_scav_zone * coefs%dCeq_dMs - coefs%Ceq) / &
                          & (total_SIV_in_scav_zone - total_water_in_scav_zone * coefs%Ceq) ** 2
            coefs%dCeq_dMa = - coefs%Ceq / (2. * coefs%Ceq  + fB)  ! d Ceq / d Ma - non-SO2 acids, p.37
            coefs%dWC_dMa = scav_zone_volume * &                  ! d Wc / d M_a, p.37
                          & total_SIV_in_scav_zone * coefs%dCeq_dMa / &
                          & (total_SIV_in_scav_zone - total_water_in_scav_zone * coefs%Ceq) ** 2
          endif  ! if SO2
        endif   ! if internal coefs are requested
        
      enddo !species-wise scavenging 
    enddo  ! sources
              
  end subroutine scavenge_puff_2011


  !**************************************************************************************

  subroutine scavenge_column_2011_adjoint(sens_mass_in_air, &    ! sensitivity in the air   (nSpecies, nSrc, n3d)
                                        & sens_mass_in_water, &  ! sensit mass in precipitation (sSpecies, nSrc, nz)
                                        & sens_arLowMassThreshold, &  ! for sensitivity 
                                        & TL_store_mass_in_air, &   ! mass in air, from TL storage (nSpecies,1,n3d)
                                        & species, &
                                        & Nspecies, Nsources, n3d, &
                                        & timestep_sec_abs, &
                                        & rulesDeposition, &
                                        & r1d, pTL_L0_in, pTL_dLdm_in, pAdj_in)
    !
    ! Scavenges sensitivity from a puff into precipitation, adjoint run
    ! caller has to handle the stuff from mass_in_water.
    ! In Eulerian case -- puff is a dispersion cell
    ! In Lagrangian case -- puff is a particle
    ! The trouble is in saturation: that one comes out of total mass stored as a part of TLA
    ! while the scavenged fraction is applied to the sensitivity
    ! Since the upper levels affect the saturation and scavenging efficiency of lower layers,
    ! the problem is non-local and tranformation matrix has a dimension (nSpecies*nLayers) x (nSpecies*nLayers)
    ! It consists mostly of zeroes but storing and computing even zeroes is expensive.
    ! However: (i) lower layers do not affect saturation of the upper ones => 
    !                     run top-down with accumulation (forwards)
    !         (ii) in adjoint, the transposed matrix has no propagation from uppwer levels down => 
    !                     bottom-up run with accumulation
    !        (iii) TL matrix consists of two parts: L(m0) and dL/dm, to be multiplied with (m-m0)
    !                     but differential adjoint has only dL/dm. m0 is mass stored as TLA point
    !         (iv) non-scavengeable, linearly scavengeable and perfectly scavengeable species are trivial
    ! See notebook 10a, pp.22-35 and notebook 11a pp. 51-59
    ! mas
  
    ! Imported parameters
    real, dimension(:,:), intent(in) :: TL_store_mass_in_air   ! (nSpecies, n3d)
    real, dimension(:,:,:), intent(inout) :: sens_mass_in_air, sens_mass_in_water ! (nSpecies, nSrc, n3d)
    type (Train_in_cell), dimension(n3d),  intent(in) :: r1d                 ! rain info for the whole column
    type(Tdeposition_rules), intent(in) :: rulesDeposition
    real, dimension(:), intent(in) :: sens_arLowMassThreshold
    type(silam_species), dimension(:), intent(in) :: species
    integer , intent (in) :: Nspecies, Nsources, n3d
    real, intent(in) :: timestep_sec_abs
    type(silja_rp_2d), dimension(:), optional, pointer :: pTL_L0_in, pTL_dLdm_in, pAdj_in

    ! Local parameters
    real, dimension(max_species), parameter :: TL_arLowMassThreshold = 0.0  ! cannot skip cells
    
!    ! Local variables
    integer :: iLev, iSp, iTmp, jTmp, iAAS, iSign, iSrc
    real :: fTmp
    real, dimension(Nspecies, 1, n3d) :: TL_mass_in_water
    real, dimension(Nspecies, 1) :: TL_store_mass_in_air_local   ! (nSpecies, nSrc=1)
    type(Tscav_coefs_1d) :: c1d  ! internal coefficients stored from forward model
    real, dimension(:,:), pointer :: pTL_L0, pTL_dLdm, pAdj
    logical :: ifBusy_pTL_L0 , ifBusy_pTL_dLdm, ifBusy_pAdj  !!

    !
    ! Get the coefficients that depend on the system state. Essentially, we run the 
    ! forward non-linear model with TL- stored masses extracting the needed parameters.
    ! For TL, we need only the main mass, i.e. nSrc=1 in the TL-store. Note top-down scanning
    !
    TL_mass_in_water(1:nSpecies,1,1:n3d) = 0.
    do iLev = n3d, 1, -1
      ! nullify the internal coefficients, hide and reformat input array
      c1d%c(iLev) = coefs_zero
      TL_store_mass_in_air_local(1:nSpecies, 1) = TL_store_mass_in_air(1:nSpecies, iLev)
      !
      ! run full forward model storing its coefficients
      if (r1d(iLev)%ifRainyCell) then
        call scavenge_puff_2011(TL_store_mass_in_air_local, &    ! stuff in the air, need (nSpecies,nSrc=1)
                              & TL_mass_in_water(:,:,iLev), & ! stuff in precipitation, need (nSpecies, nSrc=1)
                              & TL_arLowMassThreshold, &  ! what to do if forward model did not have mass ???
                              & species, &
                              & Nspecies, &
                              & Nsources, &
                              & timestep_sec_abs, &  ! positive
                              & rulesDeposition, &
                              & r1d(iLev), r1d(iLev)%cell_volume, r1d(iLev)%cell_zsize, &
                              & c1d%c(iLev))     ! internal coefficients
      endif
      if(iLev>1)TL_mass_in_water(:,:,iLev-1) = TL_mass_in_water(:,:,iLev)
    end do  ! iLev, inverse order
    !
    ! Having the basic coefficients stored, can compute the TL matrix. For each species but SO2 it is
    ! a squared matrix n3d x n3d with above-diagonal and diagonal elemnents non-zero (notebook 10a, p.28)
    ! We shall process species one by one, which allows for processing nAcidityAffectingSpecies matrices 
    ! sized n3d x n3d instead of one huge nSpecies*n3d x nSpecies*n3d.
    ! Note that for non-SO2 species the scavenging problem is linear, i.e. this very triangle matrix is
    ! both the exact solution and the tangent-linear.
    ! Exception is SO2, which is affected by all acidity affecting species. Its TL matrix is then
    ! rectangular n3d x nAcidityAffectingSpecies*n3d. 
    !
    do iSp = 1, nSpecies
      !
      ! For the time being, SO2 is treated in exactly the same way as all others. See Notebook 11 p.49
      ! for explanations.
      !
      ifBusy_pTL_dLdm = .false.
      ifBusy_pAdj = .false. 
      ifBusy_pTL_L0 = .false.
      if (iSp == indSO2)then
        !
        ! SO2 has complicated matrix dependent on other species and own concentrations
        ! Its matrix consists of nAcidityAffectingSpecies squared matrices n3d x n3d, notebook 11a p.32
        !
        if(present(pTL_L0_in))then
          pTL_L0 => pTL_L0_in(iSp)%pp
        else
          pTL_L0 => fu_work_array_2d(n3d, n3d*nAcidityAffectingSpecies)   ! L0=L(m0) for sulphur is same n3d x n3d matrix
          ifBusy_pTL_L0 = .true.
        endif
        if(present(pTL_dLdm_in))then
          pTL_dLdm => pTL_dLdm_in(iSp)%pp
        else
          pTL_dLdm => fu_work_array_2d(n3d, n3d*nAcidityAffectingSpecies)   ! Ldiff for SO2 is affected by other species
          ifBusy_pTL_dLdm = .true.
        endif
        pTL_L0(1:n3d, 1:n3d*nAcidityAffectingSpecies) = 0.0
        pTL_dLdm(1:n3d, 1:n3d*nAcidityAffectingSpecies) = 0.0
        do iAAS = 1, nAcidityAffectingSpecies  ! non-zero matrices are only for acidity-affecting species
          iSign = sign(1, indAcidityAffectingSp(iAAS))  ! if alcaline, it will help SO2 solubility
          !
          ! L(m0) the same as n3d x n3d (all info is stored in scav_coeff) but Ldiff is 
          ! is sensitive to the above-cumulated amount of other acidity-affecting species
          !
          if (indAcidityAffectingSp(iAAS) == indSO2)then
            call fill_triangle_L0(pTL_L0, (iAAS-1)*n3d, indAcidityAffectingSp(iAAS)*iSign, n3d, r1d, c1d)
          endif
!          call fill_triangle_dLdm(pTL_dLdm, (iAAS-1)*n3d, indAcidityAffectingSp(iAAS)*iSign, n3d, r1d, c1d)
          !!!
          !!! Now, scale the matrix with the corresponding partial derivatives, layer by layer
          !!!
          !!if(indAcidityAffectingSp(iAAS) == indSO2)then
          !!  do iTmp = 1, n3d
          !!    pTL_dLdm(iTmp,1+iAAS:n3d+iAAS) = pTL_dLdm(iTmp,1+iAAS:n3d+iAAS) * &   ! dMs / dmSO2, non-local 
          !!                   & (c1d%c(iTmp)%dLsc_dWc + c1d%c(iTmp)%dLic_dWc) * &  ! d lambda / d Wc
          !!                   & c1d%c(iTmp)%dWC_dMs         ! d Wc / dMs  - sensitivity to itself
          !!  end do
          !!else
          !!  do iTmp = 1, n3d
          !!    pTL_dLdm(iTmp,1+iAAS:n3d+iAAS) = pTL_dLdm(iTmp,1+iAAS:n3d+iAAS) * &   ! dMs / dmSO2, non-local 
          !!                   & (c1d%c(iTmp)%dLsc_dWc + c1d%c(iTmp)%dLic_dWc) * &  ! d lambda / d Wc
          !!                   & c1d%c(iTmp)%dWC_dMa * iSign       ! d Wc / dMa - sensitivity to others
          !!  end do
          !!endif  ! acid or SO2 itself
        end do  ! all species
        !
        ! Having TL matrix made, transpose it to adjoint
        !
        if(present(pAdj_in))then
          pAdj => pAdj_in(iSp)%pp
        else
          pAdj => fu_work_array_2d(n3d*nAcidityAffectingSpecies, n3d)
          ifBusy_pAdj = .true.
        endif
        do iTmp = 1,n3d
          pAdj(1:n3d*nAcidityAffectingSpecies, iTmp) = pTL_dLdm(iTmp, 1:n3d*nAcidityAffectingSpecies)
        enddo

      else
        !
        ! non-SO2 species depend only on own overlaying and current-level concentrations.
        ! See notebook mas 10a, pp. 27-28
        !
        ! Write down explicit tangent-linear matricx for this species
        !
        if(present(pTL_L0_in))then
          pTL_L0 => pTL_L0_in(iSp)%pp
        else
          pTL_L0 => fu_work_array_2d(n3d, n3d)
          ifBusy_pTL_L0 = .true.
        endif
        if(present(pTL_dLdm_in))then
          pTL_dLdm => pTL_dLdm_in(iSp)%pp
        else
          pTL_dLdm => pTL_L0
        endif
        call fill_triangle_L0(pTL_L0, 0, iSp, n3d, r1d, c1d)        ! Basic case: non-SO2 gas
        !
        ! Having TL matrix made, transpose it to adjoint
        !
        if(present(pAdj_in))then
          pAdj => pAdj_in(iSp)%pp
        else
          pAdj => fu_work_array_2d(n3d, n3d)
          ifBusy_pAdj = .true.
        endif
        do iTmp = 1,n3d
          pAdj(1:n3d,iTmp) = pTL_dLdm(iTmp,1:n3d)
        enddo
      endif  ! if SO2
      !
      ! Having done the TL and Adjoint matrices, now can make the adjoint scavenging
      !
      do iLev = 1, n3d
!        do iSrc = 1, nSources
!          sens_mass_in_air(iSp, iSrc, ilev) = sens_mass_in_air(iSp, iSrc, ilev) + &
!                                            & pAdj(iLev,iLev) * sens_mass_in_air(iSp, iSrc, ilev)
!        end do
        
          fTmp = 0.0
          if(iSp == indSO2)then  
            ! SO2 depends on several species
            do iAAS = 1, nAcidityAffectingSpecies
!              call msg('SO2 Diagonal only')
!              fTmp = fTmp + pTL_dLdm(iSp)%pp(iLev, iLev+(iAAS-1)*n3d) * &
!                          & sens_mass_in_air(abs(indAcidityAffectingSp(iAAS)), 1, iLev)
!              call msg('SO2 Triangle full')
              do iTmp = 1, n3d
                fTmp = fTmp + pTL_dLdm(iLev, iTmp+(iAAS-1)*n3d) * &
                            & sens_mass_in_air(abs(indAcidityAffectingSp(iAAS)), 1, iTmp)
              end do  ! iTmp 1:n3d
            end do
          else
            ! non-SO2, squared TL matrix
            fTmp = fTmp + pAdj(iLev,iLev) * sens_mass_in_air(iSp,1,iLev)
!            call msg('All Triangle full 1')
            do iTmp = 1, n3d
              fTmp = fTmp + pTL_dLdm(iLev,iTmp) * sens_mass_in_air(iSp,1,iTmp)
            end do  ! iTmp 1:n3d
          endif  ! iSp == indSO2
          sens_mass_in_air(iSp, 1, iLev) = sens_mass_in_air(iSp, 1, iLev) + fTmp
          sens_mass_in_water(iSp, 1, iLev) = sens_mass_in_water(iSp, 1, iLev) - fTmp
      end do  ! iLev

      if(ifBusy_pTL_L0) call free_work_array(pTL_L0)
      if(ifBusy_pTL_dLdm) call free_work_array(pTL_dLdm)
      if(ifBusy_pAdj) call free_work_array(pAdj)

    end do  ! iSp=1:nSpecoes

    contains
    
      !================================================================================

      subroutine fill_triangle_L0(pTL, iShift, iSp, n3d, r1d, c1d)
        !
        ! Triagonal matrix of this type (p.28 notebook 10a mas) is presented in several incarnations
        ! Its physical meaning is the impact of upper layers to the amoung droplet can pick at the current one.
        ! Specific coefficients are different for non-SO2 acidic gases and SO2, which has water_capacity
        ! dependent on already dissolved SO2.
        !
        implicit none
        !
        ! Imported parameters
        type (Train_in_cell), dimension(n3d), intent(in) :: r1d                 ! rain info for the whole column
        type(Tscav_coefs_1d), intent(in) :: c1d  ! internal coefficients stored from forward model
        real, dimension(:,:), pointer :: pTL
        integer, intent(in) :: iSp, n3d, iShift
        ! Local varisbles
        real :: fTmp
        integer iTmp, jTmp

        pTL(1:n3d,1+iShift:n3d+iShift) = 0.0     ! nullify the region to fill-in
        pTL(n3d,n3d+iShift) = -c1d%c(n3d)%scav_coef_tau(iSp)
        pTL(n3d-1,n3d-1+iShift) = -c1d%c(n3d-1)%scav_coef_tau(iSp)
        pTL(n3d-1,n3d+iShift) = c1d%c(n3d)%scav_coef_tau(iSp) * c1d%c(n3d-1)%scav_coef_tau(iSp) * &
                                                                     & c1d%c(n3d-1)%A(iSp)
        do iTmp = n3d-2, 1, -1
          pTL(iTmp,iTmp+iShift) = -c1d%c(iTmp)%scav_coef_tau(iSp)
          fTmp = c1d%c(iTmp)%scav_coef_tau(iSp) * c1d%c(iTmp)%A(iSp)
          pTL(iTmp,iTmp+1+iShift) = c1d%c(iTmp+1)%scav_coef_tau(iSp) * fTmp
          do jTmp = iTmp+2, n3d
            fTmp = fTmp * (1. - c1d%c(jTmp-1)%A(iSp) * c1d%c(jTmp-1)%scav_coef_tau(iSp))
            pTL(iTmp,jTmp+iShift) = c1d%c(jTmp)%scav_coef_tau(iSp) * fTmp
          end do
        enddo  ! iTmp=1:n3d-2
        pTL(1:n3d,1+iShift:n3d+iShift) = pTL(1:n3d,1+iShift:n3d+iShift) * r1d(1)%fCloudCover
      end subroutine fill_triangle_L0
      
      !====================================================================
      
      subroutine fill_triangle_dLdm(pTL_diff, iShift, iSp, n3d, r1d, c1d)
        !
        ! Triagonal matrix of this type (p.28 notebook 10a mas) is presented in several incarnations
        ! Its physical meaning is the impact of upper layers to the amoung droplet can pick at the current one.
        ! Specific coefficients are different for non-SO2 acidic gases and SO2, which has water_capacity
        ! dependent on already dissolved SO2.
        !
        implicit none
        !
        ! Imported parameters
        type (Train_in_cell), dimension(n3d), intent(in) :: r1d  ! rain info for the whole column
        type(Tscav_coefs_1d), intent(in) :: c1d  ! internal coefficients stored from forward model
        real, dimension(:,:), pointer :: pTL_diff
        integer, intent(in) :: iSp, n3d, iShift
        ! Local varisbles
        real :: fTmp
        integer iTmp, jTmp

        pTL_diff(1:n3d,1+iShift:n3d+iShift) = 0.0     ! nullify the region to fill-in

        call set_error('Does not work yet','fill_triangle_dLdm')

      end subroutine fill_triangle_dLdm
      
  end subroutine scavenge_column_2011_adjoint

  
  !****************************************************************************************      
      
  subroutine material_properties(species, r, water_capacitance_factor, efficiency)
    !
    ! Material properties with respect to rain properties
    !
    implicit none
    
    ! Imported paramters
    type(silam_species), intent(in) :: species
    real, intent(out) :: water_capacitance_factor, efficiency
    type (Train_in_cell), intent(in) :: r
    
    ! Local variables
    type(Tgas_deposition_param) :: DepData
    type(TwetParticle) :: wetParticle
    real :: fSc, fWetDiam, fKn, fCun, fDiffusivity, HenryConst
    real :: f_omegainv, f_phi, f_StEff, fTau
    
    water_capacitance_factor = real_missing
    if ( species%mode == in_gas_phase) then

        DepData = fu_gas_deposition_param(species%material)

        if (DepData%Henry_const_298K < 0.) return

        ! Should be some more inteligent way to figure out that the stuff does not deposit....

        fSc = r%nu_air / DepData%gas_molecular_diffus_in_air 
        ! drop  collection efficiency 
        efficiency =  4/(r%drop_Re*fSc) * &
             & (2 + 0.5657*sqrt(r%drop_Re)*(fSc**0.33333+0.4*sqrt(fSc)))*r%is_liquid
                               ! no gas scavenging by snow
        HenryConst = DepData%Henry_const_298K &
                        * exp (DepData%Henry_const_T_dep * ((1./r%fTemperature)-(1./298.)) )

      else ! particles
        ! debug temperature!!!
        wetParticle = fu_wet_particle_features(species%material, 1.03) 
                                              !As high as it can be
        fWetDiam = wetParticle%fGrowthFactor * fu_massmean_D(species%mode)
          !         if(fWetDiam > 25e-6 .and. wetParticle%fGrowthFactor > 5)then
          !            fWetDiam = 25e-6
          !         endif
        fKn = 2*r%lambda_air / fWetDiam
        fCun = 1. + fKn * (1.257 + 0.4 * exp(-1.1*fKn))
                      
        fDiffusivity = boltzmann_const * r%fTemperature * fCun / (3. * Pi * r%mu_air * fWetDiam)
        fSc = r%nu_air / fDiffusivity 
        fTau =  fWetDiam * fWetDiam * wetParticle%fWetParticleDensity * &
                                    & fCun / (18. * r%mu_air)
        f_phi=fWetDiam/r%drop_size
        f_omegainv = 0.012 ! omega = mu_w/mu_air 
                                 !   mu_water(0C) = 1.8 mPa s
                                 ! mu_water(25C)  = 0.9 mPa s
                                 ! mu_air  = 1.8e-5 Pa s 

        ! Effective Stokes minus critical stokes
        f_StEff = 2. * fTau * r%drop_vel / r%drop_size - 1./sqrt(r%drop_Re) - 0.08
        f_StEff = max(f_StEff, 1e-5) ! always positive 

            ! fraction of drop frontal area 
            ! As in SP98, Note that Re there is based on drop radius
            ! Here Re is based on drop diameter
        efficiency = 4. / (r%drop_Re*fSc) * &
                                & (2 + 0.5657*sqrt(r%drop_Re)*(fSc**0.33333+0.4*sqrt(fSc))) &
                                & + 4*f_phi*(f_omegainv*r%is_liquid + (1+2*sqrt(0.5*r%drop_Re)*f_phi)) &
                                & + exp(-1.8/(f_StEff+sqrt(f_StEff)))
                                ! diffusion + interception + impaction

        HenryConst = 1e10 ! As much as you can imagine

    endif !gas/particle 
                   
      ! see  notebook R1 p. 34.
      ! (equivalent m^3 air)/(kg water):
             
    water_capacitance_factor =  gas_constant_uni * r%fTemperature * HenryConst ! According to Henry const

  end subroutine material_properties


  !********************************************************************************************

  subroutine make_rain_in_cell_col(iLev, nLevels, r, scavType, max_scav_rate_depends_on, &
       & metdat_col, pwcColumnLayer, pwcColumn, pwcAboveT, pwcAboveB, cwcColumnLayer, cwcColumn, cwcAbove, landfrac, &
       & rain_rate_scf, tcc, cape_met, abl_height, max_scav_rate_wind_scaling, max_scav_rate_cape_scaling)
    ! Prepares precipitation parameters for a dispersion cell,
    ! for which there is ready-made columnar meteo data
                                                
    implicit none

    ! Imported parameters             
    integer, intent(in) :: iLev, nLevels  ! Dispersion indices
    integer, intent(in) :: scavType, max_scav_rate_depends_on
    real, dimension(:,:), intent(in) :: metdat_col
    real, intent(in) :: rain_rate_scf, tcc, cape_met, abl_height, landfrac
    real, intent(in) :: pwcColumnLayer, pwcColumn, pwcAboveT, pwcAboveB, cwcColumnLayer, cwcColumn, cwcAbove
    real, intent(in) :: max_scav_rate_wind_scaling, max_scav_rate_cape_scaling
    
    type (Train_in_cell), intent(out) :: r

    ! Local variables 
    real  :: fTmp, fTmp1, precip_rate, u, v, wCld, dxDisp, fEntrain, Phi
    real :: cwcAboveBottom, cwcAboveTop, pwcAboveTop, pwcAboveBottom
    !Integer :: iDisp, iMeteo

    R%ifRealCloud = .TRUE.  ! Use meteo clouds by default                                                                                                                                   
    r%ifRainSfc = .TRUE.
    r%ifRainyCell = .TRUE.
    r%ifValid = .True.

    r%fRainRateSfc = rain_rate_scf

    r%fCloudCover = min(0.999,max(tcc, 1e-2)) !tcc
    r%fPressure = metdat_col(imet_press, iLev)
    r%fTemperature = metdat_col(imet_temp, iLev)
    r%pwcColumn = 0. !!! Not sure if it is right choice, but if scavType is not in scav2018, scav2020, 
                     !! it stays undefined and is used in if, which is wrong in any case

    select case(scavType)
      case(scav2018, scav2020)

        if (max_scav_rate_depends_on == horiz_wind) then
          u = metdat_col(imet_u, iLev)
          v = metdat_col(imet_v, iLev)
          r%max_scav_rate = max_scav_rate_wind_scaling * sqrt(u*u+v*v+1.)/10000.0 !!! 1m/s minimum cell turnover speed
        else if (max_scav_rate_depends_on == cape) then
          ! "optimal" max_scav_rate_cape_scaling:
          ! 1.0 for global 0.5 x 0.5 deg run with ERA5
          ! 0.5 for European run with ECMWF operational meteodata
          r%max_scav_rate = max_scav_rate_cape_scaling * 0.01 * (1. + sqrt(cape_met))/max(500.0, abl_height)
        else if (max_scav_rate_depends_on == cape_and_wind) then
          ! "good" max_scav_rate_cape_scaling:
          ! 1.0 for global 0.5 x 0.5 deg run with ERA5
          ! 0.4 for European run with ECMWF operational meteodata
          ! "good" max_scav_rate_wind_scaling:
          ! 1.0 for global 0.5 x 0.5 deg run with ERA5
          ! 0.3 for European run with ECMWF operational meteodata
          u = metdat_col(imet_u, iLev)
          v = metdat_col(imet_v, iLev) 
          r%max_scav_rate = max_scav_rate_cape_scaling * 0.01 * (1. + sqrt(cape_met))/max(500.0, abl_height)
          r%max_scav_rate = landfrac * r%max_scav_rate + &
                & (1.0-landfrac) * max_scav_rate_wind_scaling * sqrt(u*u+v*v+1.)/10000.0
        end if

        if (scavType == scav2018) then
          r%pwcColumn = metdat_col(imet_pwc3d, 1)
        else if (scavType == scav2020) then
          r%pwcColumn = pwcColumn
        end if

      case(scav2018entr)
        r%pwcColumn = metdat_col(imet_pwc3d, 1)
        !                                                                                                              
        ! Entrainment is considered as a function of the cloud updraft velocity. See notebook 11a, pp.63-69            
        !                                                                                                              
        wCld = 3.0                ! m/s, quite stable but depends on condensation rate, range 1-10 m/s                 
        dxDisp = 0.5 * metdat_col(imet_dx_size,1)  + metdat_col(imet_dy_size,1)  ! m, S_entrain / V_cell                    
        fEntrain = 0.1            ! unitless, u_entr_horiz / wCld, 0.1 from observations, subject for tuning           
        Phi = wCld * fEntrain / dxDisp
        r%max_scav_rate = Phi / max(1e-3, r%fCloudCover * (1.- r%fCloudCover))

      case default
        !                                                                                                              
        ! No limitation due to entrainment                                                                             
        !                                                                                                              
        r%max_scav_rate = 1e3  ! pretty much infinity for this quantity                                                
     end select

    if (r%pwcColumn > 1e-25) then ! this may be put to a higher value, but not before testing that
                                  ! the global species budgets are not altered
      
      r%cwcColumnLayer = cwcColumnLayer
      cwcAboveBottom = cwcAbove
      r%cwcColumn = cwcColumn

      if (scavType /= scav2020) then
        pwcAboveBottom = metdat_col(imet_pwc3d, iLev)
        if (iLev < nLevels) then
          pwcAboveTop = metdat_col(imet_pwc3d, iLev+1)
        else
          pwcAboveTop = 0.0
        end if
        r%pwcColumnLayer = pwcAboveBottom - pwcAboveTop
      else
        r%pwcColumnLayer = pwcColumnLayer
        pwcAboveBottom = pwcAboveB
        pwcAboveTop = pwcAboveT
      end if
    else
      r%cwcColumn = r%fRainRateSfc * 3600. ! hourly rain                                                                                                                                                          
      r%ifRealCloud = .False.
      cwcAboveBottom =  r%cwcColumn * fu_frac_cwc_above(r%fPressure)
      if (cwcAboveBottom < 1e-20 * r%cwcColumn) then
        r%ifRainyCell = .FALSE.
        return
      endif
      if (iLev < nLevels) then
        ftmp = metdat_col(imet_press, iLev+1)
      else   ! last level
        ftmp = metdat_col(imet_press, iLev-1)
        fTmp = 2*r%fPressure - fTmp  !!extrapolate pressure                 
      endif
      cwcAboveTop =  r%cwcColumn * fu_frac_cwc_above(fTmp)
      r%cwcColumnLayer = cwcAboveBottom - cwcAboveTop

      pwcAboveTop = cwcAboveTop
      r%pwcColumnLayer = r%cwcColumnLayer !!This is fake cloud -- copy stuff   
      r%pwcColumn = r%cwcColumn
      pwcAboveBottom = cwcAboveBottom
   end if

    if (r%cwcColumnLayer < 1e-7 * cwcAboveBottom) r%cwcColumnLayer = 0.  !! Avoid numerics
    if (r%pwcColumnLayer < 1e-7 * pwcAboveBottom) r%pwcColumnLayer = 0.  !! Avoid numerics

    if (pwcAboveBottom < 1e-7*r%pwcColumn) then ! could be done earlier...
      r%ifRainyCell = .FALSE.
      return
   endif

    ! how much rain get formed in the cell?                                                                            
    fTmp = r%fRainRateSfc / r%pwcColumn
    r%deltarain = fTmp * r%pwcColumnLayer  ! rain int. for in-cloud scavenging                                         
    r%rain_rate_in = fTmp * pwcAboveTop ! how much rain comes from above?

    r%cell_zsize = metdat_col(imet_dz_size, iLev) !fu_layer_thickness_m(fu_level(dispersion_vertical, ilev))
    r%cell_volume = metdat_col(imet_dx_size, 1) * metdat_col(imet_dy_size, 1) * r%cell_zsize

    r%mu_air  = fu_dynamic_viscosity(r%fTemperature)
    r%rho_air = r%fPressure / (gas_constant_dryair * r%fTemperature)
    r%nu_air  =  r%mu_air / r%rho_air

    ! Precipitation properties                                                                                                
    precip_rate = r%rain_rate_in / r%fCloudCover ! actual rain intensity where it rains                                       

    if(r%fTemperature > freezing_point_of_water)then ! rain                                                                   
      r%is_liquid = 1.

      r%drop_size = max(1e-5, min(1e-2, 7e-4 * sqrt(sqrt(precip_rate*3600.))))  ! Slinn 1977 Eq. 14 for steady fronal rain    
      ! limits added to allow for very strong / small rain                                                                    

      if (r%drop_size - r%drop_size .ne. 0.) then
        call msg("Drop size", r%drop_size)
        call msg("precip_rate", precip_rate)
        call msg("fTempr", r%fTemperature )
        call msg("rho_air", r%rho_air )
        call set_error("NaN scavenging coefficient...", "scavenging")
      endif

      r%drop_mass =1000.*0.166667*3.1415936*r%drop_size*r%drop_size*r%drop_size  ! 1/6 rho_w * Pi *d^3                        

      ! Gossard (1992, JTECH) Saturation due to drop shape change Density correction could be                                 
      ! power 0.4 (Foote, 1969,JAM)                                                                                           
      r%drop_vel  = 10. * (1.-exp(-5e2*r%drop_size))*sqrt(1.2/r%rho_air) ! m/s                                                
    else
      r%is_liquid = 0. ! melted_size = 1.4e-3 * sqrt(precip_rate*3600.)  median volume MELTED diameter                        
                       ! Power 0.45 in Cekhon&Srivastava                                                                      
    !drop_mass = 1000.*0.166667*3.1415926*melted_size*melted_size*melted_size                                                 

      r%drop_size = 1e-3 ! m                                                                                                  
      r%drop_mass = 3e-8 ! kg                                                                                                 
      r%drop_vel  = 0.8  ! m/s   As done in ECHHAM-5. (Croft 2009, ACP) Looks ugly, but...                                    
                    ! Some dependencies from Matrosov 2007, Gossard 1992 and other radar papers could be better               
    endif

    r%drop_Re = r%drop_vel*r%drop_size/r%nu_air

    ! Scavenging rate in cloud for particles                                                                                  
    r%lambda_air = 2.37e-5 * r%fTemperature / r%fPressure

    contains
      real function fu_frac_cwc_above(pressure)   
        real, intent(in) :: pressure

        !fake cloud  structure: constant density cloud,    
        if (pressure>9e4) then ! subcloud       
          fu_frac_cwc_above = 1
        elseif (pressure > 7e4) then !incloud          
          fu_frac_cwc_above =  (pressure-7e4)/2e4 ! 900hPa-700hPa   
        else
          fu_frac_cwc_above = 0.
        endif
      end function fu_frac_cwc_above

    end subroutine make_rain_in_cell_col

  !*********************************************************************************

  subroutine make_rain_in_cell(iX, iY, iLev, nLevels, r, &
                             & weight_past, &
                             & pHorizInterpStruct, pVertInterpStruct,&
                             & ifHorizInterp, ifVertInterp, scavType)

    ! Used only for Lagrangian runs!!!!!!
    ! Prepares parameters of precipitation within a dispersion cell: fills Train_in_cell structure               
    implicit none

    ! Imported parameters                                                                                        
    integer, intent(in) :: iX, iY, iLev, nLevels  ! Dispersion indices                                           
    type (Train_in_cell), intent(out) :: r
    real, intent(in) ::  weight_past
    type(THorizInterpStruct), intent(in) :: pHorizInterpStruct
    type(TVertInterpStruct), intent(in) :: pVertInterpStruct
    logical, intent(in) :: ifHorizInterp, ifVertInterp
    integer, intent(in) :: scavType
    real, dimension(:), pointer :: xSizePtr, ySizePtr

    ! Local variables                                                                                            
    real ::  cwcAboveBottom, cwcAboveTop, pwcAboveBottom, pwcAboveTop
    real  :: fTmp, fTmp1, precip_rate, u, v, wCld, dxDisp, fEntrain, Phi
    integer :: iDisp, iMeteo

    r%ifRealCloud = .TRUE.  ! Use meteo clouds by default                                                        
    r%ifRainSfc = .TRUE.
    r%ifRainyCell = .TRUE.
    r%ifValid = .TRUE.
    !                                                                                                            
    ! Meteo                                                                                                      
    !                                                                                                            
    iMeteo = fu_grid_index(nx_meteo, iX, iY, pHorizInterpStruct)

    r%fRainRateSfc = pMetPrecTot(iMeteo)

    if (r%fRainRateSfc < 1e-7/3600.) then  ! <1e-7 mm/h is not a rain                                            
      r%ifRainSfc = .FALSE.
      r%ifRainyCell = .FALSE.
      return
   endif

    r%fCloudCover = min(0.999,max(pMetTotCloud(iMeteo), 1e-2))  ! No rain without clouds ->                      
                                                                ! Force 1% cloud cover                           
    r%fPressure  = fu_get_value(fldPressure, &   ! field_ptr_4d                                                  
                              & nx_meteo, &     ! nxFrom                                                         
                              & ix, iy, iLev, &
                              & weight_past, &
                              & pHorizInterpStruct, pVertInterpStruct, ifHorizInterp, ifVertInterp)
    r%fTemperature  = fu_get_value(fldTempr, &   ! field_ptr_4d                                                  
                                 & nx_meteo, &     ! nxFrom                                                      
                                 & ix, iy, iLev, &
                                 & weight_past, &
                                 & pHorizInterpStruct, pVertInterpStruct, ifHorizInterp, ifVertInterp)
    
    r%max_scav_rate = 1e3

    ! Cloud water content                                                                                        
    !                                                                                                            
    if (scavType == scav2011) then
      r%cwcColumn =  fu_get_value(fldCwcabove, nx_meteo, ix, iy, 1 , 0.0, &
           & pHorizInterpStruct, pVertInterpStruct, ifHorizInterp, .false.)
      r%pwcColumn = r%cwcColumn
    else
       r%cwcColumn = real_missing
       r%pwcColumn = real_missing
    endif

    if (r%pwcColumn > 1e-4) then ! 0.1 mm cloud water ...                                                        

      cwcAboveBottom =  fu_get_value(fldCwcabove, nx_meteo, ix, iy, iLev , 0.0, &
                                   & pHorizInterpStruct, pVertInterpStruct, ifHorizInterp, ifVertInterp)

      pwcAboveBottom = cwcAboveBottom

      if (iLev < nLevels) then
        cwcAboveTop   =  fu_get_value(fldCwcabove, nx_meteo, ix, iy, iLev + 1, 0.0, &
                                    & pHorizInterpStruct, pVertInterpStruct, ifHorizInterp, ifVertInterp)
        pwcAboveTop = cwcAboveTop

      else
        cwcAboveTop = 0.
        pwcAboveTop = 0.
      endif

      r%cwcColumnLayer = cwcAboveBottom - cwcAboveTop
      r%pwcColumnLayer = pwcAboveBottom - pwcAboveTop
      if (pwcAboveBottom < 1e-10 * r%pwcColumn) then
        r%ifRainyCell = .FALSE.
        return
      endif

    else  ! Some cloud structure to be invented                                                                  

      r%cwcColumn = r%fRainRateSfc * 3600. ! hourly rain                                                         
      r%ifRealCloud = .False.
      cwcAboveBottom =  r%cwcColumn * fu_frac_cwc_above(r%fPressure)
      if (cwcAboveBottom < 1e-10 * r%cwcColumn) then
        r%ifRainyCell = .FALSE.
        return
      endif
      if (iLev < nLevels) then
        fTmp = fu_get_value(fldPressure, &   ! Pressure above                                                    
                          & nx_meteo, ix, iy, iLev + 1, weight_past, &
                          & pHorizInterpStruct, pVertInterpStruct, ifHorizInterp, ifVertInterp)
      else   ! last level                                                                                        
        fTmp = fu_get_value(fldPressure, &   ! Pressure below                                                    
                          & nx_meteo, ix, iy, iLev - 1, weight_past, &
                          & pHorizInterpStruct, pVertInterpStruct, ifHorizInterp, ifVertInterp)
        fTmp = 2*r%fPressure - fTmp  !!extrapolate pressure                                                      
      endif
      cwcAboveTop =  r%cwcColumn * fu_frac_cwc_above(fTmp)
      r%cwcColumnLayer = cwcAboveBottom - cwcAboveTop

      pwcAboveTop = cwcAboveTop
      r%pwcColumnLayer = r%cwcColumnLayer !!This is fake cloud -- copy stuff                                     
      r%pwcColumn = r%cwcColumn
      pwcAboveBottom = cwcAboveBottom
    endif

    if (r%cwcColumnLayer < 1e-6 * cwcAboveBottom) r%cwcColumnLayer = 0.  !! Avoid numerics                       
    if (r%pwcColumnLayer < 1e-6 * pwcAboveBottom) r%pwcColumnLayer = 0.  !! Avoid numerics                       

    ! how much rain get formed in the cell?                                                                      
    fTmp = r%fRainRateSfc / r%pwcColumn
    r%deltarain = fTmp * r%pwcColumnLayer  ! rain int. for in-cloud scavenging                                   
    r%rain_rate_in = fTmp * pwcAboveTop ! how much rain comes from above?                                        

#ifdef DEBUG
    if (.not. (r%deltarain < 1. .and. r%deltarain >= 0.)) call ooops("Strange r%deltarain")
    if (.not. (r%rain_rate_in < 1. .and. r%rain_rate_in >= 0.)) call ooops("Strange r%rain_rate_in")
    if (.not. (r%cwcColumnLayer >= 0.)) call ooops("Strange r%cwcColumnLayer")
    if (.not. (r%pwcColumnLayer >= 0.)) call ooops("Strange r%pwcColumnLayer")
#endif

    if (pwcAboveBottom < 1e-10*r%pwcColumn) then ! could be done earlier...                                      
      ! less than 5% of cloud above the bottom of this cell all rain is below                                    
      r%ifRainyCell = .FALSE.
      return
    endif

    ! cell parameters                                                                                            

    iDisp = ix+(iy-1)*nx_dispersion
    xSizePtr => fu_grid_data(dispersion_cell_x_size_fld)
    ySizePtr => fu_grid_data(dispersion_cell_y_size_fld)
    r%cell_zsize = fldCellSizeZ%past%p2d(iLev)%ptr(iDisp)   * weight_past &
             & + fldCellSizeZ%future%p2d(iLev)%ptr(iDisp) * (1. - weight_past)
    r%cell_volume = xSizePtr(iDisp) * ySizePtr(iDisp) * r%cell_zsize

    r%mu_air  = fu_dynamic_viscosity(r%fTemperature)
    r%rho_air = r%fPressure / (gas_constant_dryair * r%fTemperature)
    r%nu_air  =  r%mu_air / r%rho_air

    ! Precipitation properties                                                                                   
    precip_rate = r%rain_rate_in / r%fCloudCover ! actual rain intensity where it rains                          

    if(r%fTemperature > freezing_point_of_water)then ! rain                                                      
      r%is_liquid = 1.

      r%drop_size = max(1e-5, min(1e-2, 7e-4 * sqrt(sqrt(precip_rate*3600.))))  ! Slinn 1977 Eq. 14 for steady fronal rain                                                                                                       
      ! limits added to allow for very strong / small rain                                                       

      if (r%drop_size - r%drop_size .ne. 0.) then
        call msg("Drop size", r%drop_size)
        call msg("precip_rate", precip_rate)
        call msg("fTempr", r%fTemperature )
        call msg("rho_air", r%rho_air )
        call set_error("NaN scavenging coefficient...", "scavenging")
      endif

      r%drop_mass =1000.*0.166667*3.1415936*r%drop_size*r%drop_size*r%drop_size  ! 1/6 rho_w * Pi *d^3           

      ! Gossard (1992, JTECH) Saturation due to drop shape change Density correction could be                    
      ! power 0.4 (Foote, 1969,JAM)                                                                              
      r%drop_vel  = 10. * (1.-exp(-5e2*r%drop_size))*sqrt(1.2/r%rho_air) ! m/s                                   
    else
      r%is_liquid = 0. ! melted_size = 1.4e-3 * sqrt(precip_rate*3600.)  median volume MELTED diameter           
                       ! Power 0.45 in Cekhon&Srivastava                                                         
    !drop_mass = 1000.*0.166667*3.1415926*melted_size*melted_size*melted_size                                    

      r%drop_size = 1e-3 ! m                                                                                     
      r%drop_mass = 3e-8 ! kg                                                                                    
      r%drop_vel  = 0.8  ! m/s   As done in ECHHAM-5. (Croft 2009, ACP) Looks ugly, but...                       
                    ! Some dependencies from Matrosov 2007, Gossard 1992 and other radar papers could be better  
    endif

    r%drop_Re = r%drop_vel*r%drop_size/r%nu_air

    ! Scavenging rate in cloud for particles                                                                     

    r%lambda_air = 2.37e-5 * r%fTemperature / r%fPressure

    contains

      real function fu_frac_cwc_above(pressure)   
        real, intent(in) :: pressure

        !fake cloud  structure: constant density cloud,                                                          
        if (pressure>9e4) then ! subcloud                                                                        
          fu_frac_cwc_above = 1
        elseif (pressure > 7e4) then !incloud                                                                    
          fu_frac_cwc_above =  (pressure-7e4)/2e4 ! 900hPa-700hPa                                                
        else
          fu_frac_cwc_above = 0.
        endif
      end function fu_frac_cwc_above

  end subroutine make_rain_in_cell


  !*********************************************************************************


  subroutine scavenge_lagr(src_masses, dest_masses, &
                         & timestep_sec, weight_past, &
                         & p_dyn, &  ! slice of lp_dyn
                         & DispNevels, &
                         & speciesTransport, nSpecies, &
                         & rulesDeposition, &
                         & pMeteo2DispInterpHoriz,  pMeteo2DispInterpVert, &
                         & ifVertInterp, ifHorizInterp &
                        )
    !
    ! Computes the wet depositon of species in the lagrangian particle.
    !
    implicit none

    ! Imported parameters
    real, dimension(:,:), intent(inout) :: src_masses, dest_masses  ! (nSpecies)
    real, dimension(:),  intent(in) :: p_dyn
    real, intent(in) :: timestep_sec, weight_past
    integer, intent(in) :: nSpecies, DispNevels
    type(Tdeposition_rules), intent(in) :: rulesDeposition
    type(THorizInterpStruct), pointer, intent(in) :: pMeteo2DispInterpHoriz
    type(TVertInterpStruct), pointer, intent(in) :: pMeteo2DispInterpVert
    logical, intent(in)  :: ifVertInterp, ifHorizInterp
    integer :: indexMeteoHoriz, indexMeteoVert
    type(silam_species), dimension(:), pointer :: speciesTransport

    ! Local variables
    real :: scav_coef_std, scav_coef, fTempr_1, scav_amt
    integer :: iSpecies
    real, dimension(:), pointer :: fptrTmp
    real, dimension(max_species), target :: rSettling
    real, dimension(max_species,1) :: tmpDepo
    type (Train_in_cell) :: r

    !
    ! Simply go through the map applying the requested type of deposition
    !
    select case (rulesDeposition%scavengingType)

      case(scavStandard)
       indexMeteoHoriz = fu_grid_index(nx_Meteo, nint(p_dyn(lp_x)),nint(p_dyn(lp_y)), &
                                    & pMeteo2DispInterpHoriz)
       indexMeteoVert = fu_vert_index(nint(p_dyn(lp_x)),nint(p_Dyn(lp_y)), &
                                & nint(p_dyn(lp_z)), pMeteo2DispInterpVert, ifVertInterp)
       !
       ! Standard 3D scavenging coefficient. Each material has own scaling applied as a 
       ! correction factor: material-wise for gases, size-wise for aerosols.
       !
       ! For aerosols, we will want to have a settling velocity as a scaling parameter
       ! with regard to some reference value of a "default" aerosol.
       !
       scav_coef_std = fldScavStd%past%p2d(indexMeteoVert)%ptr(indexMeteoHoriz) * weight_past + &
                     & fldScavStd%future%p2d(indexMeteoVert)%ptr(indexMeteoHoriz) * (1. - weight_past) * &
                     & timestep_sec

       if(scav_coef_std .eps. 0.0)return  ! speedup

       if(rulesDeposition%nAerosolsDepositing > 0)then
         if(error)return
         call get_settling_velocity_species(speciesTransport, nSpecies, &
                                          & 1, 1, 1.0, &      ! indexMeteo, iLev, weight_past
                                          & rulesDeposition, &
                                          & rSettling)
         !
         ! Default settling velocity is for 1 micrometre aerosol with 1000 kg/m3 density,
         ! which is about 3.3e-5 m/sec
         !
         do iSpecies = 1, nSpecies
           if(rSettling(iSpecies) > 1.0e-10)then  ! 3e-10 is a formal settling of 1nm particle
             rSettling(iSpecies) = min(100.0, rSettling(iSpecies) * 3.0e4)
           else
             rSettling(iSpecies) = 1.0 ! gases or whatever non-depositing for whatever reasons
           endif
         end do

       endif  ! if there are aerosols
       !
       !
       fTempr_1 = fldTempr%past%p2d(1)%ptr(indexMeteoHoriz) * weight_past + &
                & fldScavStd%future%p2d(1)%ptr(indexMeteoHoriz) * (1. - weight_past)
       !
       ! The scaling with regard to standard 3D coefficient depends on the phase. 
       ! Gases keep their scaling in the chemical database
       !
       do iSpecies = 1, rulesDeposition%nGasesDepositing

!         if(.not. fu_if_gas_depositing(speciesTransport( &
!                                      & rulesDeposition%indGasesDepositing(iSpecies))%material))cycle

         if(fTempr_1 > freezing_point_of_water)then
           scav_coef = scav_coef_std * &
                   & fu_scavenging_scale_gas_rain( &
                          & speciesTransport(rulesDeposition%indGasesDepositing(iSpecies))%material)
         else
           scav_coef = scav_coef_std * &
                   & fu_scavenging_scale_gas_snow( &
                          & speciesTransport(rulesDeposition%indGasesDepositing(iSpecies))%material)
         endif
         scav_coef = scav_coef * (1. + scav_coef * (0.5 + scav_coef / 3.))  ! implicit 3rd order

         scav_amt = src_masses(rulesDeposition%indGasesDepositing(iSpecies),1) * scav_coef/(1. + scav_coef)
         src_masses(rulesDeposition%indGasesDepositing(iSpecies),1) = &
                                & src_masses(rulesDeposition%indGasesDepositing(iSpecies),1) - scav_amt
         dest_masses(rulesDeposition%indGasesDepositing(iSpecies),1) = &
                                & dest_masses(rulesDeposition%indGasesDepositing(iSpecies),1) + scav_amt
       end do  !  gases

       !
       ! Aerosols scaling can be considered either via the settling velocity or via 
       ! the material properties, just like gases. The default option is probably the size.
       !
       do iSpecies = 1, rulesDeposition%nAerosolsDepositing
         scav_coef = scav_coef_std * rSettling(rulesDeposition%indAerosolDepositing(iSpecies))
         scav_coef = scav_coef * (1. + scav_coef * (0.5 + scav_coef / 3.))  ! implicit 3rd order

         scav_amt = src_masses(rulesDeposition%indAerosolDepositing(iSpecies),1) * scav_coef/(1. + scav_coef)
         src_masses(rulesDeposition%indAerosolDepositing(iSpecies),1) = &
                            & src_masses(rulesDeposition%indAerosolDepositing(iSpecies),1) - scav_amt
         dest_masses(rulesDeposition%indAerosolDepositing(iSpecies),1) = &
                            & dest_masses(rulesDeposition%indAerosolDepositing(iSpecies),1) + scav_amt
       end do  !  aerosols


       !
       ! New Scavenging scheme
       !
      case(scavNoScav)
!        call msg('Dummy scavenging -- no scavenging at all!')

      case(scav2011,scav2011fc)
         !call msg('New 2011 scavenging scheme...',rulesDeposition%scavengingType)
         
        call  make_rain_in_cell(nint(p_dyn(lp_x)), nint(p_Dyn(lp_y)), nint(p_dyn(lp_z)), &
                       & DispNevels, r,  weight_past, &
                       & pMeteo2DispInterpHoriz,  pMeteo2DispInterpVert, &
                       & ifHorizInterp, ifVertInterp, rulesDeposition%scavengingType)

        if (.not. r%ifRainyCell) return
        fptrTmp => rSettling  ! recycle array
        fptrTmp(1:nSpecies) = 0. !Dummy arLowMassThreshold
        tmpdepo(1:nSpecies,1) = 0. ! clean rain from above
             ! call msg("***Making rain, cell "+fu_str(ix)+","+fu_str(iY)+', level= '+fu_str(iLev))
             !             call report_rain(r)
        call  scavenge_puff_2011( src_masses, &    ! stuff in the air
                     &  tmpdepo, & ! stuff in precipitation
                     &  fptrTmp, &   !   arLowMassThreshold, &
                     &  speciesTransport,      & 
                     &  nSpecies,&
                     &  1, & ! nSources
                     &  timestep_sec,&
                     &  rulesDeposition, &
                     &  r,  &
                     &  p_Dyn(lp_dx) * p_Dyn(lp_dy) * p_Dyn(lp_dz) , &  !puff volume
                     &  p_Dyn(lp_dz))
             
         ! Finally deposit rain contents....
         dest_masses(1:nSpecies,1) = dest_masses(1:nSpecies,1) + tmpdepo(1:nSpecies,1)

      case default
       call msg('Unknown scavenging type',rulesDeposition%scavengingType)
       call set_error('Unknown scavenging type','scavenge_lagr')

    end select ! type of scavenging

  end subroutine scavenge_lagr


  subroutine allocate_scav_amount(mapTransport)
    implicit none
    type(Tmass_map),  intent(in) :: mapTransport 
    integer :: nthreads, istat
    character(len=*), parameter :: sub_name = 'allocate_scav_amount'

    nthreads = 1
    !$OMP PARALLEL DEFAULT(NONE) &    
    !$OMP SHARED(nthreads)     
    !$OMP MASTER    
    !$ nthreads = omp_get_num_threads()    
    !$OMP END MASTER                
    !$OMP END PARALLEL  

    allocate(scavAmount(mapTransport%nSpecies, mapTransport%nSrc, mapTransport%n3D, 0:nthreads), stat=istat)
    if(fu_fails(istat == 0, 'Failed to allocate the scavenging amount array', sub_name))return
    scavAmount(1:mapTransport%nSpecies, 1:mapTransport%nSrc, 1:mapTransport%n3D, 0:nthreads) = 0.

  end subroutine allocate_scav_amount

!************************************************************************************

  ! Ugly implementation of ugly parameterization
  real function  fu_vd_Zhang(tau,Sc,dp,ustar, Zref, RsIn, MOs, sfc_typ)
!       tau particle        tau_p , seconds
!       Schmidt number      Sc   = \nu / D  
!       particle size       dp, meters
!       friction velocity   ustar, m/s
!       reference height    Zref, m
!       surface resistance  R_s, s/m
!       Stability parameter MOs = 1/L_{MO}, 1/m
!       reference height    Zref, m
!       surface type        
!       
!       returns             v_d, m/s
        implicit none
        real, INTENT(IN) :: tau,Sc,dp,ustar, Zref, RsIn, MOs
        integer, intent(IN) :: sfc_typ
        ! implemented according to SP06 with approximation for unstable
        ! stratification
        ! Note that different Rb are used there for gases and particles
        ! Particles are distinguished by non-zero taup or vsplus
        real :: ra,rb,rs
        real :: gama, alpha ! Surface parameters
        real :: A, z0       !  Surface length scales
        real :: Vdif, Vimp, Vint, Vs ! dimensionless components in denominator
                                     !of R_b
        real :: R ! gebound
        real :: StZ ! Stokes number in Zhang interpretation
        real, parameter :: g = 9.8 !m/s -- gravity
        real, parameter :: nu = 1.5e-5 !m^2/s -- kinematic viscosity
        

        select case (sfc_typ)
                case (ZhangSmooth)
                        z0 = nu/ustar
                        gama = 0.5
                        alpha = 100.
                        StZ  = tau*ustar/z0
                        Vint = 0
                        R=1
                case (ZhangSnow)
                        z0 = 2e-4
                        gama =0.54
                        alpha = 1.2
                        A = 5e-3 
                        StZ = tau*ustar/A
                        Vint = 0.5*dp*dp/(A*A)
                case (ZhangGrass)
                        z0 = 1e-2
                        gama = 0.54 
                        alpha = 1.2
                        A = 2e-3
                        StZ = tau*ustar/A
                        Vint = 0.5*dp*dp/(A*A)
                case (ZhangBlForest)
                        z0 = 1.
                        gama = 0.56
                        alpha = 0.8
                        A = 5e-3
                        StZ = tau*ustar/A
                        Vint = 0.5*dp*dp/(A*A)
                case default
                 fu_vd_Zhang = 0.
                 return
        end select
        
        Rb = 1./(3.*ustar) / &
               & (Sc**(-gama) + (StZ/(StZ+alpha))**2 + Vint )
        

        ra = 2.5*  (log(Zref/z0) + fu_Psi(Zref*MOs)) / ustar
        ra = 0.5 * (ra + abs(ra))

        vs = g * tau

        if (tau + abs(vs) .gt. 1e-10) then !particles
                rs = ra*rb*vs ! virtual resistance
        else
!            !   rb = 5*Sc**.6666667
                rs = RsIn
        endif
        
        fu_vd_Zhang = vs + 1./(ra+rb+rs)
        return
  end function fu_vd_Zhang


  !************************************************************************************

  real function  fu_vdplus_slinnSP98(taup,Sc,vsplus, Rsplus, MOplus, Zpmax)
        implicit none
        real, INTENT(IN) :: taup ,Sc ,vsplus, Rsplus, MOplus, Zpmax
        ! implemented according to SP98 with approximation for unstable
        ! stratification
        ! Note that different Rb are used there for gases and particles
        ! Particles are distinguished by non-zero taup or vsplus
        real ra,rb,rs
        
        ra = 2.5 * (log(Zpmax) + fu_Psi(Zpmax*MOplus))
        ra = 0.5 * (ra + abs(ra))

        if (taup + abs(vsplus) .gt. 1e-10) then !particles
                rb = 1./ ( Sc**(-.6666667) + 10**(-3./taup))
                rs = ra*rb*vsplus ! virtual resistance
        else
               rb = 5*Sc**.6666667
               rs = Rsplus
        endif
        
        fu_vdplus_slinnSP98 = vsplus + 1./(ra+rb+rs)
        return
  end function fu_vdplus_slinnSP98

  !************************************************************************************

  real function  fu_vdplus_DS(taup,Sc,vsplus, Rsplus, MOplus, Zpmax)
        implicit none
        real, INTENT(IN) :: taup ,Sc ,vsplus, Rsplus, MOplus, Zpmax
       ! Simple vd for smooth surfaces with diffusion and settling.
       ! Minimalistic deposition scheme that deposites any sizes.
       ! Coefficient 0.1 chousen to fit smooth-surface diffusion at  Sc >> 1
        fu_vdplus_DS = vsplus + 0.1*Sc**(-.666667)
        return
  end function fu_vdplus_DS

  !***************************************************************************


  real function fu_vdplus_smooth(taup,Sc,vsplus,rplus, Rsplus, MOplus, Zpmax)
        implicit none 
        real, INTENT(IN) :: taup ,Sc ,vsplus ,rplus ,Rsplus, MOplus, Zpmax
        ! Copyleft Rostislav Kouznetsov, FMI (firstname.lastname@fmi.fi)
        ! See: Kouznetsov, R., and M. Sofiev (2012), A methodology for evaluation of vertical dispersion and dry deposition
        ! of atmospheric aerosols, J. Geophys. Res., 117, D01202, doi:10.1029/2011JD016366.

!       tau plus            taup = tau_p  u_*^2 / \nu   
!       Schmidt number      Sc   = \nu / D  
!       settling velocity   vsplus = v_s / u_*; usually v_s = g * tau_p 
!       particle size       rplus = d_p u_* / (2\nu)
!       surface resistance  Rsplus  = u_* R_s
!       reference height    zpmax  = z u_* / \nu
!       Stability parameter MOplus = \nu/(u_* L)?
!       
!       returns            v_d / u_*
!
!	(Dimensionless) deposition velocity through split layer...
!	Analytics for resistances, Lagrangian time is infinite except for 
!       turbophoretic layer, where \tau_L=5, picewize fit for turbophoretic velocity. 
!       Structure of the layers Bottom->top:
!               zplus=rplus; particle radius  
!	           laminar layer: Vs,  Full diffusivity
!	        zplus=30*Sc^{-1/3}; top of laminar layer (\nu_t >> D)
!	           buffer layer:  Vs, Eddy diff
!	        zplus=4; 	
!	            Turbophoretic layer  Vs+Vtf, Eddy diff
!	        zplus=25; 
!	             turbulent layer Vs, Eddy diff 
!               zplus=Zpmax
!	
        !local parameters
        real :: Zl, Ztf2
        real, parameter :: Zbuf=3.
        real, parameter :: Ztf = 18.   ! turbophoretic sublayer heights
        real, parameter :: taultf = 5. ! Lagrangian time in turbophoretic layer
        real :: V, R, R1, S,fTmp, fTmp1, fIvd, x ! Velocity, resistance, Sc**0.33333
        real :: Il, It, Nutp                ! integrals (functions)
        
        Nutp(x)=.4*x*x*x/(x*x+200)
        ! integrals
        !\int (1+x^3)^{-1}dx
        ! 3**(-0.5) = 0.57735; 
        Il(x) = -0.16667 * log(x*x-x+1) + 0.57735 * atan((2*x-1)*0.57735) + 0.3333*log(x+1)

        ! \int (x^2+200)/(0.4*x^3) dx 
        It(x) = 2.5*log(x) - 100./(x*x) 

!	Laminar+buffer
        S=Sc**0.3333
        Zl=20. / S ! laminar layer height

        if (Zl .gt. Zbuf) then ! Very diffusive thing -- no turbophoresis
                               ! only diffusion and turbulent diffusion 
                               ! stiched in ad-hoc manner.
                               ! Regular velocity
                               ! (electricity, thermophoresis etc.)
                               !  is not accounted here
                if (Zl .gt. Zpmax) then ! Whole layer is diffusive
                                        !should never happen normally 
                       fu_vdplus_smooth=1./(Rsplus + (Zpmax)*Sc)
                       return
                endif
                !solving equation nu_t(Zl) = 0.4*Zl^3/(Zl^2+200) = 1/Sc
                ! to find a point where the turbulent diffusion equals to molecular
               fTmp= 2.5/Sc
               fTmp1 = (fTmp*fTmp*fTmp/27. + &
                        fTmp*(100. +5.*sqrt(8.*fTmp/27.+400.)) )** 0.3333333
               Zl= fTmp1 + fTmp*fTmp/(9.*fTmp1) + 0.333333*fTmp
               ! linear approximation of nu_t around Zl
                !  (d nu / dz)(Zl)
                fTmp = Zl*Zl/(Zl*Zl+200)
                fTmp = 1.2* fTmp - 0.8*fTmp*fTmp;
                ! half-width of transitional layer to calculate with Simpsons parabola
                 fTmp1 = 1./(Sc*fTmp) 
                 fTmp = 1./Sc;
                
                fIvd = Rsplus + (Zl-fTmp1)*Sc  &
                        + 0.333333*fTmp1*( &
                          1./(fTmp + Nutp(Zl-fTmp1))+ 4./(fTmp + Nutp(Zl))+ 1./(fTmp + Nutp(Zl+fTmp1))&
                         ) ! up to the top of transitional layer
                fTmp =  It(Zpmax) - It(Zl+fTmp1) + 2.5*fu_Psi(Zpmax*MOplus) !R_a
                fIvd = fIvd + 0.5*(fTmp+abs(fTmp)) ! no negative R_a allowed
               fu_vdplus_smooth=1./fIvd
               return
        endif
        


        ! 500^0.333 = 7.92
        R = 7.92 * S*S * (  Il(Zl*S/7.92) - Il(rplus*S/7.92) ) !laminar resistance
        R = 0.5*(R+abs(R))  ! should be zero if rplus>Zl;
        ! Buffer layer
        Zl = 0.5*(rplus+Zl + abs(rplus-Zl)) ! max(rplus,Zl) - lower limit of buffer layer
        R1 = It(Zbuf) - It(Zl);             ! buffer resistance
        R = R + 0.5*(R1+abs(R1));           ! total resistance
        
        fTmp = vsplus*R
        if (abs(fTmp) .gt. 0.001) then
                if (fTmp < -30) then 
                        !fIvd turns to infinity at large negative vsplus
                        ! backward-time large particles
                       fu_vdplus_smooth = 0
                       return
                endif
                fTmp = exp(-fTmp)
                fIvd = Rsplus*fTmp +  (1-fTmp)/vsplus
        else
                fIvd = Rsplus + R
        endif
        if (fIvd /= fIvd) then
                call msg("Strange fu_vdplus_smooth after laminar",fu_vdplus_smooth)
        endif
        

        !turbophoretic layer
        V = 0.81 * taup / (Ztf-Zbuf)/(1+taup/taultf)*sign(1.,vsplus) + vsplus
        R = (It(Ztf) - It(Zbuf))*(1+taup/taultf);
        fTmp = V*R
        if (abs(fTmp) .gt. 0.001) then
                if (fTmp < -30) then 
                        !fIvd turns to infinity at large negative vsplus
                        ! backward-time large particles
                       fu_vdplus_smooth = 0
                       return
                endif
                fTmp = exp(-fTmp)
                fIvd = fIvd*fTmp +  (1-fTmp)/V
        else
                fIvd = fIvd + R
        endif

        ! Not exactly accurate but saves when the center of mass is at the
        ! surface....               
        if (Zpmax .lt. Ztf) then
                fu_vdplus_smooth=1./fIvd
                return
        endif
        

        ! Lagrangian turbophoretic layer 
        Ztf2 = 2*taup
        if (Ztf2 .gt. Zpmax) then !shold not really happen, 
                                  !but just in case
                Ztf2=Zpmax
        endif

!        if (1 .eq. 0) then
        if (Ztf2 .gt. Ztf) then ! Lagrangian turbophoretic layer present
                                ! tau_L = 0.5 zplus
                V = 0.4 + vsplus

                R =  0.16667 * ( &
                            (1 + taup/(0.5*Ztf))        / Nutp(Ztf) + &
                       4. * (1 +taup/(0.25*(Ztf+Ztf2))) / Nutp((0.5*(Ztf+Ztf2)))+&
                            (1 + taup/(0.5*Ztf2))       / Nutp(Ztf2) &
                            ) * (Ztf2 - Ztf)
 !               write (*,*) "#Rtf2-parabolic", R
!                R = -5.*taup*(1./Ztf2 - 1./Ztf)
 !               write (*,*) "#Rtf2-quad", R
                fTmp = V*R
                if (abs(fTmp) .gt. 0.001) then
                        fTmp = exp(-fTmp)
                        fIvd = fIvd*fTmp +  (1-fTmp)/V
                else
                        fIvd = fIvd + R
                endif
        else
                
                Ztf2 = Ztf
        endif
        if (fIvd /= fIvd) then
          call msg('taup,Sc',taup,Sc)
          call msg('vsplus,rplus',vsplus,rplus)
          call msg('Rsplus, MOplus', Rsplus, MOplus)
          call msg('Zpmax', Zpmax)
          call msg("Strange fu_vdplus_smooth1",fu_vdplus_smooth)
        endif
        !aerodynamic layer

        R = It(Zpmax) - It(Ztf2) + 2.5*fu_Psi(Zpmax*MOplus) 
                                ! Stability correction applied to
                                ! the upper limit of integration
        R = 0.5* (R + abs(R)) ! Could happen to be negative
                              ! in very unstable case
        fTmp = vsplus*R
        if (abs(fTmp) .gt. 0.001) then
                fTmp = exp(-fTmp)
                fIvd = fIvd*fTmp +  (1-fTmp)/vsplus
        else
                fIvd = fIvd + R
        endif

        fu_vdplus_smooth=1./fIvd
        if (fIvd /= fIvd) then
          call msg('taup,Sc',taup,Sc)
          call msg('vsplus,rplus',vsplus,rplus)
          call msg('Rsplus, MOplus', Rsplus, MOplus)
          call msg('Zpmax', Zpmax)
          call msg("Strange fu_vdplus_smooth2",fu_vdplus_smooth)
        endif

        return
  end function fu_vdplus_smooth


!************************************************************************************


  ! adaptive step runge-kutta for smooth surface
  ! Not to be actually used. Just to check the 
  ! approximation used in fu_vdplus_smooth
  real function fu_vdplus_smooth_rkad(taup,Sc,vsplus,rplus, Rsplus, Zpmax,tol)
        implicit none 
        real, INTENT(IN) :: taup ,Sc ,vsplus ,rplus ,Rsplus ,Zpmax,tol
!       tau plus           taup = tau_p  u_*^2 / \nu   
!       Schmidt number     Sc   = \nu / D  
!       settling velocity  vsplus = v_s / u_*; usually v_s = g * tau_p 
!       particle size      rplus = d_p u_* / (2\nu)
!       surface resistance Rsplus  = u_* R_s
!       reference height   zpmax  = z u_* / \nu
!       
!       returns            v_d / u_*
!	(Dimensionless) deposition velocity through constant flux layer...
        double precision  :: R,R1,R2,Rprev,z,zstep,dz, steperr
        integer :: ngood, nbad
        ngood=0
        nbad=0
        z=rplus 
        zstep=0.5*rplus;
        R=Rsplus
        steperr = tol /(Zpmax-rplus)
        do 
                R1=R;
                R2=R;
                call iVrk4(z,R1,zstep,taup ,Sc ,vsplus)
                call iVrk4(z,R2,0.5*zstep,taup ,Sc ,vsplus)
                call iVrk4(z+0.5*zstep,R2,0.5*zstep,taup ,Sc ,vsplus)
                
                if ((abs(R1-R2) .gt. steperr*zstep) .and. (zstep .gt. 1e-4)) then
                                zstep = 0.1*zstep
                                nbad = nbad + 1 
!                                write (*,*) z, zstep, R1, R2
                else
                        ngood =ngood + 1
                        R=R1
                        zstep = zstep * 1.123456789
                        z=z+zstep  
                        if  ( z .gt. Zpmax) then ! almost finished
                                if  ( z .gt. Zpmax + 0.999*zstep ) then
                                               exit
                                endif
                                zstep = zstep - Zpmax + z
                                z = Zpmax
                        endif
                endif
        enddo
!write (*,*) '# ngood=', ngood,  '# nbad=', nbad
                fu_vdplus_smooth_rkad = 1./R
! 
        return
  end function fu_vdplus_smooth_rkad

  !************************************************************************************

  !   fourth-order Runge-Kutta step subroutine  for smooth surfaces
  SUBROUTINE iVrk4(t, y, tstep,taup ,Sc ,vsplus)
       IMPLICIT none
       double precision  ::  h, t, tstep, y, x
        real :: taup ,Sc ,vsplus
        double precision  ::  k1, k2,k3, k4, temp1, temp2, temp3
        
        h=tstep/2.0
      k1 = tstep * DERIV(t, y,taup ,Sc ,vsplus)
      temp1 = y + 0.5*k1
      k2 = tstep * DERIV(t+h, temp1,taup ,Sc ,vsplus)
      temp2 = y + 0.5*k2
      k3 = tstep * DERIV(t+h, temp2,taup ,Sc ,vsplus)
      temp3 = y + k3
      k4 = tstep * DERIV(t+tstep, temp3,taup ,Sc ,vsplus)
      y = y + (k1 + (2.*(k2 + k3)) + k4)/6.0
       RETURN
  END SUBROUTINE iVrk4

!************************************************************************************


  double precision function DERIV(x,y,taup ,Sc ,vsplus)
        implicit none 
        double precision :: x,y
        real :: vsplus, Sc, taup
        real :: vtf, nu
 
! Translated with Maxima       
!        nu(x):=2/5 * x^3 /(x^2+200);
!        w(x):=9/10 * x^2/ (x^2 + 90);
!        tauL(x):=nu(x)/(w(x)^2);
!        fortran(diff(w(x)*w(x)/(1+ taup/tauL(x)),x));
!        fortran  (nu(x)/(1+ taup/tauL(x)));

        vtf = vsplus + taup *( &
              &  8.1E+1*x**3/(2.5E+1*(x**2+90)**2*(8.1E+1*taup*x*(x**2+200)/&
              &  (4.0E+1*(x**2+90)**2)+1))+(-8.1E+1)*x**5/(2.5E+1*(x**2+90)**3&
              &  *(8.1E+1*taup*x*(x**2+200)/(4.0E+1*(x**2+90)**2)+1))+(-8.1E+1)&
              &  *x**4*(8.1E+1*taup*(x**2+200)/(4.0E+1*(x**2+90)**2)+(-8.1E+1)&
              &  *taup*x**2*(x**2+200)/(1.0E+1*(x**2+90)**3)+8.1E+1*taup*x**2/&
              &  (2.0E+1*(x**2+90)**2))/(1.0E+2*(x**2+90)**2*(8.1E+1*taup*x*&
              &  (x**2+200)/(4.0E+1*(x**2+90)**2)+1)**2)&
              &  )
        nu = 1./Sc + &
                2.0E+0*x**3/(5.0E+0*(x**2+200)*(8.1E+1*taup*x*(x**2+200)/(4.0E+1*(x**2+90)**2)+1))
        DERIV = (1-vtf*y)/nu
        return
  end  function DERIV

!************************************************************************************

  double precision function DERIV2(x,y,taup ,Sc ,vsplus)
        implicit none 
        double precision :: x,y
        real :: vsplus, Sc, taup
        real :: vtf, nu
 
! Translated with Maxima from 
!        nu(x):=2/5 * x^3 /(x^2+10*x+200);
!        w(x):=9/10 * x^2/ (x^2+ 2*x + 91);
!        tauL(x):=nu(x)/(w(x)^2);
!        fortran(diff(w(x)*w(x)/(1+ taup/tauL(x)),x));
!        fortran  (nu(x)/(1+ taup/tauL(x)));

        vtf = vsplus + taup *( &
             &   8.1E+1*x**3/(2.5E+1*(x**2+2*x+91)**2*(8.1E+1*taup*x*(x**2+10*x+200)&
             &   /(4.0E+1*(x**2+2*x+91)**2)+1))+(-8.1E+1)*x**4*(2*x+2)/(5.0E+1*&
             &   (x**2+2*x+91)**3*(8.1E+1*taup*x*(x**2+10*x+200)/(4.0E+1*(x**2+2&
             &   *x+91)**2)+1))+(-8.1E+1)*x**4*(8.1E+1*taup*(x**2+10*x+200)/(4.0E+1&
             &   *(x**2+2*x+91)**2)+(-8.1E+1)*taup*x*(2*x+2)*(x**2+10*x+200)/&
             &   (2.0E+1*(x**2+2*x+91)**3)+8.1E+1*taup*x*(2*x+10)/(4.0E+1*(x**2+&
             &   2*x+91)**2))/(1.0E+2*(x**2+2*x+91)**2*(8.1E+1*taup*x*(x**2+10*x&
             &   +200)/(4.0E+1*(x**2+2*x+91)**2)+1)**2)&
             &   )
        nu = 1./Sc + &
                2.0E+0*x**3/(5.0E+0*(x**2+10*x+200)*(8.1E+1*taup*x*(x**2+10*x+200) &
                           /(4.0E+1*(x**2+2*x+91)**2)+1))
        DERIV2 = (1-vtf*y)/nu
        return
  end  function DERIV2

!************************************************************************************


  real function fu_Psi(zL)
        implicit none 
             real :: zL
             !Stability correction to neutral resistance
             ! based on SP98, eq. (19.14) on p. 963
             ! Note, that it is Psi for momentum!
             ! Claimed in SP limits are -1<z/L<1.
             ! For stable stratification used AS IS.
             ! for unstable -- simpler approximation is used.
             ! Approximation agrees to original expressions
             ! within the range -10<z/L<0 (p.25 in Roux2009 notebook)
             ! The correction should be added to log(z/z0) 
             ! WARNING: the check should be made for negative resistances!
             real :: s, u ! stable and unstable part of correction
             s = 2.35 * (zL+abs(zL))  !2.35=4.7/2
                                      ! zero for unstable
             u = 0.5*(abs(zL)-zL)     ! zero for stable
             u = -4.*u/( 2.65*sqrt(u*sqrt(u)) + 1.) ! 4x/(2.65*x^(3/4)+1)
             fu_Psi = s+u
    return
  end  function fu_Psi

!************************************************************************************


  real function fu_vdplus_rough(St,Sc,vsplus,dpa, Rsplus, Restar, MOplus, Zz0) 
        implicit none 
        real, INTENT(IN) :: St,Sc,vsplus,dpa, Rsplus, Restar, MOplus, Zz0
!       Stokes number       St = 2 \tau_p u* / a
!       Schmidt number     Sc   = \nu / D  
!       settling velocity  vsplus = v_s / u_*; usually v_s = g * \tau_p 
!       particle size      dpa = d_p / a
!       surface resistance Rsplus  = u_* R_s
!       Reynolds* number   Restar = u_* a / \nu
!       Stability parameter MOplus = z0 / L
!       reference height   Zz0  = z / z_0
!       
!       returns            v_d / u_*
!	(Dimensionless) deposition velocity 

        real :: Vdif, Vint, Vimp, Impar, etaimp
        real ::  Re12,  Ra, fIvd, fTmp ! 
        
        real, parameter :: Stcr=0.15 ! critical stokes number
        real, parameter :: Utus=3.0 ! U_{top}/u_*
        
        Re12 = sqrt(Restar)

        Vdif = 1. / (Rsplus + 0.5*Re12 * Sc**.6666667)
        Vint = 80. * dpa * dpa *Re12
        Vimp = 0;
        Impar = St - 1./(Utus * Re12) - Stcr ! Impaction parameter
        if (Impar .gt. 0) then 
                Vimp = 2./Utus * exp(-0.1/Impar- 1./sqrt(Impar)) 
        endif
        
        fTmp = Vdif+Vint+Vimp+vsplus
        if (fTmp .le. 0.) then !  no deposition
                fu_vdplus_rough = 1e-6
                return
        endif

        fIvd = 1./fTmp
       

        Ra  = 2.5 * (log(Zz0) + fu_Psi(Zz0*MOplus)) ! no stability correction
                                                      ! for lower limit
        Ra = 0.5*(Ra+abs(Ra))                   ! Could become negative in
                                                !unstable
        fTmp = vsplus*Ra
        if (abs(fTmp) .gt. 0.001) then
                fTmp = exp(-fTmp)
                fIvd = fIvd*fTmp +  (1-fTmp)/vsplus
        else
                fIvd = fIvd + Ra
        endif
        fu_vdplus_rough = 1./fIvd
        return
  end function fu_vdplus_rough

  
  !************************************************************************************

  logical function fu_if_tla_required_scav(rules) result(required)
    !
    ! Collect transformations' requests for tangent linearization. If tangent linear
    ! is not allowed, corresponding subroutine sets error. 
    !
    implicit none
    type(Tdeposition_rules), intent(in) :: rules
    
    integer :: iType
    
    ! call set_error('TLA not available for cbm4', 'fu_if_tla_required_cbm4')
    required = .false.
    do iType = 1, rules%nDepositionTypes
      select case(rules%iDepositionType(iType))
        case(scav2011, scav2011fc)
          required = .true.
          return
        case default
      end select
    end do
    
  end function fu_if_tla_required_scav

  
  !************************************************************************************
  
  subroutine test_TL_ADJ(mass_in_air, &  ! sensit mass in precipitation (sSpecies, nSrc, nz)
                       & arLowMassThreshold, &  ! for sensitivity 
                       & species, &
                       & Nspecies, Nsources, n3d, &
                       & timestep_sec_abs, &
                       & rulesDeposition, &
                       & r1d)
    !
    ! Tests tangent linear and adjoint matrices by computing non-linear forward and then
    ! repeating it with the TL. The Adjoint matrix is checked by transposition
    !
    implicit none
    
    ! Imported parameters
    real, dimension(:,:), intent(in) :: mass_in_air ! (nSpecies, n3d)
    type (Train_in_cell), dimension(n3d), intent(in) :: r1d  ! rain info for the whole column
    type(Tdeposition_rules), intent(in) :: rulesDeposition
    real, dimension(:), intent(in) :: arLowMassThreshold
    type(silam_species), dimension(:), intent(in) :: species
    integer , intent (in) :: Nspecies, Nsources, n3d
    real, intent(in) :: timestep_sec_abs

    ! Local variables
    integer, parameter :: nTimes=10
    integer :: iLev, iTime, iAAS, iSp, iTmp, indSpecies
    real, dimension(Nspecies,1,n3d,nTimes) :: Full_mass_in_air, TL_mass_in_air, &
                          & Full_mass_in_water, TL_mass_in_water   ! (nSpecies, nSrc, n3d, nTimes)
    type(Tscav_coefs_1d) :: c1d  ! internal coefficients stored from forward model
    type(silja_rp_2d), dimension(:), pointer :: pTL_L0, pTL_dLdm, pAdj
    real :: fTmp
    real, dimension(Nspecies) :: mass_factor

    mass_factor(:) = 1
    !
    ! Get the vectors of species and masses, then compute its forward non-linear scavenging
    ! over several time steps, then do the same with tangent linear matrix, which is made only 
    ! once, of course. The first time step must show identical results, then they can deviate
    !
    !
!    ! Store mass in air - the initial conditions
!    mass_in_air_local(1:Nspecies,1,1:n3d,1) = mass_in_air(1:Nspecies,1:n3d)
!    mass_in_water(1:Nspecies,1,1:n3d, 1:nTimes) = 0.0
!    ! Store the coefficients, if someone needs them
!    call msg('getting coefficients')
!    c1d%c(iLev) = coefs_zero
!    do iLev = n3d, 1, -1
!      call scavenge_puff_2011(mass_in_air_local(1:Nspecies,:,iLev,iTime), &    ! stuff in the air, need (nSpecies,nSrc=1)
!                            & mass_in_water(1:Nspecies,:,iLev,iTime), & ! stuff in precipitation, need (nSpecies, nSrc=1)
!                           & arLowMassThreshold, &  ! what to do if forward model did not have mass ???
!                           & species, &
!                            & Nspecies, &
!                            & Nsources, &
!                            & timestep_sec_abs, &  ! positive
!                            & rulesDeposition, &
!                            & r1d(iLev), r1d(iLev)%cell_volume, r1d(iLev)%cell_zsize, &
!                            & c1d%c(iLev))     ! internal coefficients
!    end do
!    ! (re)store intial masses
    
    !
    ! Make the matrices - the only thing needed now from the adjoint sub
    !
    call msg('Getting TL matrix')
    call get_work_arrays_set(nSpecies, pTL_L0)
    call msg('done pTL_L0')
    call get_work_arrays_set(nSpecies, pTL_dLdm)
    call msg('done pTL_dLdm')
    call get_work_arrays_set(nSpecies, pAdj)
    call msg('done pAdj')
    if(error)then
      call set_error('Failed some of work array sets','test_TL_ADJ')
      return
    endif
    call scavenge_column_2011_adjoint(TL_mass_in_air(1:Nspecies,:,1:n3d,1), &    ! mass in the air   (nSpecies, nSrc, n3d)
                                    & TL_mass_in_water(1:Nspecies,:,1:n3d,1), &  !  mass in precipitation (sSpecies, nSrc, nz)
                                    & arLowMassThreshold, &  ! for sensitivity 
                                    & mass_in_air, &   ! mass in air, from full-model storage (nSpecies,1,n3d)
                                    & species, &
                                    & Nspecies, Nsources, n3d, &
                                    & timestep_sec_abs, &
                                    & rulesDeposition, &
                                    & r1d, pTL_L0, pTL_dLdm, pAdj)
    ! set working masses
    !
!    mass_factor(indHNO3) = 100000
!    mass_factor(indSO2) = 100000
    do iSp = 1, nSpecies
      if(mass_factor(iSp) /= 1) &
         & call msg('===========>>>>>>>>> Species:' + fu_str(species(iSp)) + ', multiplied with', &
                  & mass_factor(iSp))
      Full_mass_in_air(iSp,1,1:n3d,1) = mass_in_air(iSp,1:n3d) * mass_factor(iSp)
      Full_mass_in_air(iSp,1,1:n3d,2) = mass_in_air(iSp,1:n3d) * mass_factor(iSp)
      TL_mass_in_air(iSp,1,1:n3d,1) = mass_in_air(iSp,1:n3d) * mass_factor(iSp)
      TL_mass_in_air(iSp,1,1:n3d,2) = mass_in_air(iSp,1:n3d) * mass_factor(iSp)
    end do
    Full_mass_in_water(1:Nspecies,1,1:n3d, 1:nTimes) = 0.0
    TL_mass_in_water(1:Nspecies,1,1:n3d, 1:nTimes) = 0.0
    
    !
    ! run full forward model
    !
    do iTime = 2, nTimes-1
      do iLev = n3d, 1, -1
        if (r1d(iLev)%ifRainyCell)then  ! there is rain at this level
!           call msg('Repository scavenging')
!          call scavenge_puff_2011_repo(mass_in_air_local(1:Nspecies,:,iLev,iTime), &    ! stuff in the air, need (nSpecies,nSrc=1)
          call msg('Adjusted scavenging')
          call scavenge_puff_2011(Full_mass_in_air(1:Nspecies,:,iLev,iTime), &    ! stuff in the air, need (nSpecies,nSrc=1)
                                & Full_mass_in_water(1:Nspecies,:,iLev,iTime), & ! stuff in precipitation, need (nSpecies, nSrc=1)
                                & arLowMassThreshold, &  ! what to do if forward model did not have mass ???
                                & species, &
                                & Nspecies, &
                                & Nsources, &
                                & timestep_sec_abs, &  ! positive
                                & rulesDeposition, &
                                & r1d(iLev), r1d(iLev)%cell_volume, r1d(iLev)%cell_zsize)
        endif  ! if there is rain in the cell
        ! Cumulate the scavenged amount
        if(iLev> 1) Full_mass_in_water(1:Nspecies,:,iLev-1,iTime) = &
                                         & Full_mass_in_water(1:Nspecies,:,iLev,iTime)
      end do  ! iLev, inverse order
      !
      ! Copy the current value to next time step
      !
      Full_mass_in_air(1:Nspecies,:,1:n3d,iTime+1) = Full_mass_in_air(1:Nspecies,:,1:n3d,iTime)
!      Full_mass_in_water(1:Nspecies,:,1:n3d,iTime+1) = Full_mass_in_water(1:Nspecies,:,1:n3d,iTime)
    end do ! iTime
    !
    ! Have nTimes-step evolution of the masses in air and in water as non-linear model computed
    ! Do nTimes steps with TL matrix
    !
    do iTime = 2, nTimes-1
      do iSp = 1, nSpecies
        do iLev = n3d, 1, -1
          fTmp = 0.0
          if(iSp == indSO2)then  
            ! SO2 depends on several species
            do iAAS = 1, nAcidityAffectingSpecies
!              call msg('SO2 Diagonal only')
!              fTmp = fTmp + pTL(iSp)%pp(iLev, iLev+(iAAS-1)*n3d) * &
!                          & TL_mass_in_air(abs(indAcidityAffectingSp(iAAS)), 1, iLev, iTime)
              call msg('SO2 Triangle full')
              do iTmp = 1, n3d
                fTmp = fTmp + pTL_L0(iSp)%pp(iLev, iTmp+(iAAS-1)*n3d) * &
                            & TL_mass_in_air(abs(indAcidityAffectingSp(iAAS)), 1, iTmp, iTime)
              end do  ! iTmp 1:n3d
            end do
          else              
            ! non-SO2, squared TL matrix
!            call msg('All Diagonal only')
!            fTmp = fTmp + pTL(iSp)%pp(iLev,iLev) * TL_mass_in_air(iSp,1,iLev,iTime)
            call msg('All Triangle full')
            do iTmp = 1, n3d
              fTmp = fTmp + pTL_L0(iSp)%pp(iLev,iTmp) * TL_mass_in_air(iSp,1,iTmp,iTime)
            end do  ! iTmp 1:n3d
          endif  ! iSp == indSO2
          ! Make the scavenging...
          TL_mass_in_air(iSp, 1, iLev, iTime) = TL_mass_in_air(iSp, 1, iLev, iTime) + fTmp
          TL_mass_in_water(iSp, 1, iLev, iTime) = TL_mass_in_water(iSp, 1, iLev, iTime) - fTmp
          ! Cumulate scavengned amount downwards
          if(iLev> 1)then
            TL_mass_in_water(iSp,1,iLev-1,iTime) = TL_mass_in_water(iSp,1,iLev,iTime)
          endif
        end do  ! iLev
      end do  ! iSp
      !
      ! Copy the current value to next time step
      !
      TL_mass_in_air(1:Nspecies,:,1:n3d,iTime+1) = TL_mass_in_air(1:Nspecies,:,1:n3d,iTime)
!      TL_mass_in_water(1:Nspecies,:,1:n3d,iTime+1) = TL_mass_in_water(1:Nspecies,:,1:n3d,iTime)
    end do  ! iTime
    !
    ! Now, just write the results down to the log file
    !
    call msg('')
    call msg('>>>>>>>>>>>>>>>>>>>>>>>>>> TEST TL SCAVENGING <<<<<<<<<<<<<<<<<<<<<<<<<<')
    call msg('')
    
    call msg('Rain parameters for iMeteo=', r1d(1)%iMeteo)
    call msg('iLev fCloudCover  fPressure  fTemperature  fRainRateSfc  cwcColumn    mu_air     rho_air     nu_air     lambda_air   rain_rate_in   deltarain   cwcColumnLayer   is_liquid    drop_size    drop_mass    drop_vel   drop_Re   cell_zsize  cell_volume')
    do iLev = 1, n3d
      call msg(fu_str(iLev),(/r1d(iLev)%fCloudCover, r1d(iLev)%fPressure, r1d(iLev)%fTemperature, &
                            & r1d(iLev)%fRainRateSfc, r1d(iLev)%cwcColumn, &
                            & r1d(iLev)%mu_air, r1d(iLev)%rho_air, r1d(iLev)%nu_air, r1d(iLev)%lambda_air, &
                            & r1d(iLev)%rain_rate_in, r1d(iLev)%deltarain, r1d(iLev)%cwcColumnLayer, &
                            & r1d(iLev)%is_liquid, r1d(iLev)%drop_size, r1d(iLev)%drop_mass, &
                            & r1d(iLev)%drop_vel, r1d(iLev)%drop_Re, &
                            & r1d(iLev)%cell_zsize, r1d(iLev)%cell_volume/))
    end do

    do iSp = 1, nSpecies
      if(sum(Full_mass_in_air(iSp,1,1:n3d,1)) < 1e-10)then
        call msg(fu_str(species(iSp)) + ', ZERO MASS')
        cycle
      else
        call msg(fu_str(species(iSp)))
      endif
      call msg('TL matrix:' + fu_str(species(iSp)) + ', L0')
      do iLev = 1, n3d
        if(iSp == indSO2)then
          call msg(fu_str(iLev) + ':', (/pTL_L0(iSp)%pp(iLev,1:n3d*nAcidityAffectingSpecies)/))
        else
          call msg(fu_str(iLev) + ':', (/pTL_L0(iSp)%pp(iLev,1:n3d)/))
        endif
      end do
      call msg('Full scheme, mass in air:' + fu_str(species(iSp)))
      call msg('time', (/(iTime,iTime=1,nTimes)/))
      do iLev = 1, n3d
        call msg(fu_str(iLev) + ',full', (/(Full_mass_in_air(iSp,1,iLev,iTime),iTime=1,nTimes)/))
      end do
      call msg('TL scheme, mass in air:' + fu_str(species(iSp)))
      call msg('time', (/(iTime,iTime=1,nTimes)/))
      do iLev = 1, n3d
        call msg(fu_str(iLev) + ',TL', (/(TL_mass_in_air(iSp,1,iLev,iTime),iTime=1,nTimes)/))
      end do
      call msg('(TL - full) / (TL + full), mass in air:' + fu_str(species(iSp)))
      call msg('time', (/(iTime,iTime=1,nTimes)/))
      do iLev = 1, n3d
        call msg(fu_str(iLev) + ',(-)/(+)', (/((TL_mass_in_air(iSp,1,iLev,iTime) - &
                                           & Full_mass_in_air(iSp,1,iLev,iTime)) / &
                                      & (TL_mass_in_air(iSp,1,iLev,iTime) + &
                                           & Full_mass_in_air(iSp,1,iLev,iTime) + 1e-20),iTime=1,nTimes)/))
      end do
      call msg('Full scheme, mass in water:' + fu_str(species(iSp)))
      call msg('time', (/(iTime,iTime=1,nTimes)/))
      do iLev = 1, n3d
        call msg(fu_str(iLev) + ',full', (/(Full_mass_in_water(iSp,1,iLev,iTime),iTime=1,nTimes)/))
      end do
      call msg('TL scheme, mass in water:' + fu_str(species(iSp)))
      call msg('time', (/(iTime,iTime=1,nTimes)/))
      do iLev = 1, n3d
        call msg(fu_str(iLev) + ',TL', (/(TL_mass_in_water(iSp,1,iLev,iTime),iTime=1,nTimes)/))
      end do
      call msg('(TL - full) / (TL + full), mass in water:' + fu_str(species(iSp)))
      call msg('time', (/(iTime,iTime=1,nTimes)/))
      do iLev = 1, n3d
        call msg(fu_str(iLev) + '(-)/(+)', (/((TL_mass_in_water(iSp,1,iLev,iTime) - &
                                           & Full_mass_in_water(iSp,1,iLev,iTime)) / &
                                      & (TL_mass_in_water(iSp,1,iLev,iTime) + &
                                           & Full_mass_in_water(iSp,1,iLev,iTime) + 1e-20),iTime=1,nTimes)/))
      end do
    end do
    call msg('')
    call msg('>>>>>>>>>>>>>>>>>>>>>>>>>> END TEST TL SCAVENGING <<<<<<<<<<<<<<<<<<<<<<<<<<')
    call msg('')
 
    call free_work_array(pTL_L0)
    call free_work_array(pTL_dLdm)
    call free_work_array(pAdj)

    
  end subroutine test_TL_ADJ

end module depositions


