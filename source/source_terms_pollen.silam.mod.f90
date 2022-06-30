MODULE source_terms_pollen
  !
 ! This module contains full description of the pollen source.
  !
  ! All units: SI, unless otherwise stated (pollen source is indeed in own units)
  !
  ! Author Mikhail Sofiev mikhail.sofiev@fmi.fi
  !
  ! Language: ANSI FORTRAN-90 (or close to it)
  !
  use cocktail_basic
  use globals, only : DEFAULT_REAL_KIND, r4k, r8k

  implicit none
  private

  public fill_pollen_src_from_namelist
  public reserve_pollen_source
  public add_input_needs
  public add_source_species_pollen_src
  public init_emission_pollen
  public source_2_second_grid
  public create_source_containing_grid
  public compute_emission_for_pollen ! Dynamic emission, including heat indices
  public fu_pollen_emis_owned_quantity
  public link_source_to_species
  public fu_name
  public fu_source_nbr
  public fu_source_id_nbr
  public typical_species_conc
  public report

  !
  ! Private routines of the pollen source
  !
  private fu_species_src
  private add_input_needs_pollen_source
  private link_pollen_src_to_species
  private project_pollen_src_second_grd
  private create_src_cont_grd_pollen_src
  private compute_HS_sigmoid_parameters
  private compute_HS_sigmoid_parameters2
  private fu_source_id_nbr_of_pollen_src
  private fu_source_nbr_of_pollen_src
  private fu_pollen_source_name
  private typical_species_cnc_pollen
  private report_pollen_src

  interface add_input_needs
    module procedure add_input_needs_pollen_source
  end interface

  interface link_source_to_species
    module procedure link_pollen_src_to_species
  end interface

  interface source_2_second_grid
    module procedure project_pollen_src_second_grd
  end interface

  interface create_source_containing_grid
    module procedure create_src_cont_grd_pollen_src
  end interface
  interface fu_source_nbr
    module procedure fu_source_nbr_of_pollen_src
  end interface

  interface fu_source_id_nbr
    module procedure fu_source_id_nbr_of_pollen_src
  end interface

  interface fu_name
    module procedure fu_pollen_source_name
  end interface

  interface typical_species_conc
    module procedure typical_species_cnc_pollen
  end interface

  interface report
    module procedure report_pollen_src
  end interface

  ! 
  ! Modified sigmoid function for heat-sum accumulation
  !
  type Theat_sum
    private
    integer :: iHeatSumType              ! degreeDay, degreeHour, bioDay, sigmoidPeriodUnits
    integer :: iChillSumType             ! degreeDay

    ! cut-linear parameters
!    integer :: indHSTemprCutOff, indHS_StartDay  ! just indices of the fields in the buffer
!    real :: fChillCutOff, fChill_StartDay, fChill_impact_HS ! parameters and effect of chill sum
    real :: chill_Tmax, chill_Topt, chill_Tmin
    integer :: chill_startDay
    logical :: chill_ifNegativeAboveMax
    
    ! sigmoid parameters
    real :: maxVal   ! = 273.15 + (fHSMaxRate-273.15)*1.1  ! +1-2 DD/day to handle T > fHSTemprSatur
    real :: TemprCutOff    ! above freesing temperature but close to it. Significantly positive only for
                           ! very warmth-loving species
    real :: smoother       ! in K: ((T-Tco) / (T-Tco+smoother)) ^ sm_power
    real :: sm_power
    real :: TemprMidPoint    ! in the middle between the cut-off and saturation, a bit to the left
                             ! to handle the asymmetric fade-in curve
    real :: exp_power        ! Slope should be 1 around the mid-point to meet the linear heatsum curve
  end type Theat_sum
  private Theat_sum
  !--------------------------------------------------------------------
  !
  !  The pollen source term
  !
  ! The source can be controlled by a series of thresholds: calendar day, heat sum, plant growth, etc.
  ! The corresponding uncertainties serve as fade-in/out of the season. Since in most cases the
  ! open-pocket principle is followed, uncertainty in the total pollen amount released is the uncertainty
  ! of the season duration, i.e. its end date
  !
  TYPE silam_pollen_source
    PRIVATE
    
    ! Source
    CHARACTER(len=clen) :: src_nm, sector_nm        ! Name of the area source and sector
    integer :: src_nbr, id_nbr                      ! A source and id numbers in a WHOLE source list
    type(silja_field), pointer :: mask_src_grid, mask_disp_grid
    TYPE(silja_grid) :: grid
    type(Tsilam_namelist), pointer :: nlInputFiles, nlSourceMask  ! namelist for names of supplementary files

    ! Pollen species
    type(silam_species), dimension(:), pointer :: species
    integer:: nSpecies
    integer :: indPolAlrg = int_missing, indFreeAlrg = int_missing, indPol = int_missing
    
    type(chemical_adaptor) :: adaptor
    
    ! Pollen amount
    real :: standardPollenTotal, &           ! magic 3.7e8 or whatever total pollen #/m2
          & fUncertainty_tot_pollen_relat    ! uncertainty of the above number: for open-pocket
                                             ! it translates into unceratinty of the season length
    ! Allergen
    real :: alrgGrowthRate, freeAlrgFract
    
    ! Flowering thresholds
    logical :: ifStartHSThr=.false., ifStartCDThr=.false., &  !start
             & ifEndHSThr=.false., ifEndCDThr=.false., & ! end
             & ifTempThr=.false., ifDayTempThr=.false., ifSWthr=.false., & ! release type dependent
             & ifStartEndGammaThr = .false., &
             & ifChillSum = .false., &                ! if chill sum is involved
             & ifStartDLThr = .false., ifEndDLThr = .false.  ! if day length is a trigger

    ! Calendar day
    type(silja_interval) :: timeMapShift ! time difference of this-year flowering from climate
    real :: fUncertainty_CD_days_start   ! uncertainty of the calendar day threshold dates
    ! Daylength
    real :: fDayLenStart_hrs, fDaylenEnd_hrs  ! day length of the start and end 
    logical :: ifShortDayTriggers             ! whether long-enough or short-enough day is the trigger
    real :: fUncertainty_DL_start_hrs, fUncertainty_DL_end_hrs  ! uncertainty
    ! Heatsum 
    real :: fUncertainty_HS_relative_start   ! uncertainty of heat sum threshold
    type(Theat_sum) :: heatsum_params
    ! Bioday 
    real :: LoTemp, HiTemp, OptTemp, &  ! Cardinal temperatures for temperatre response function
          & photoperiod                 ! photoperiod responce
    integer :: iTRfunc                  ! Temperature response function (beta, linear triangle)

    ! Gamma-type season: exp(-x/beta)*sum(sc_i * (max(x-loc_i,0))^p_i), i=1,n
    real :: gf_beta                     ! rate in dumping exponent
    real, dimension(:), pointer :: gf_scales, gf_start_times, gf_powers
    integer :: gf_nTerms

    ! Pre-season plant growth
    logical :: ifHSGrowth=.false., ifSWGrowth=.false.
    real :: SWGrowthMid, SWGrowthSigma, SWGrowthDeath

    ! Pollen ripening
    integer :: iRipeningType                 ! heatsum, time; normal, linear
    real :: dayMid1, dayMid2, daySgm1, daySgm2, dayFrac1  ! Diurnal release pattern (2 normal modes)
    
    ! Pollen release
    integer :: iReleaseType
    logical :: ems2wholeABL=.false.
    real :: LowHumidThresh, HighHumidThresh, & ! full emis if lower than low, no emis if higher than high
          & precipThreshold, &                 ! no emis if precip is above
          & windSaturation, windMaxScale, &    ! saturation of wind scaling and max impact value
          & fBufReleaseLimit, &                ! how fast the ready-to-fly-pollen can be thrown out
          & fEmissionBottom, fEmissionTop      ! emission layer    

    type(Tfade_inout) :: UncertaintyParams       ! linear or sigmoid
    
    integer :: nLevsDispVert
    real, dimension(:), pointer :: levFractDispVert, fzDisp
    type(silam_vertical) :: vertLevsDispVert
    
    type(silja_logical) :: defined
  END TYPE silam_pollen_source

  type pollen_src_ptr
    type(silam_pollen_source) :: pollen_src
  end type pollen_src_ptr
  public pollen_src_ptr

  ! Chillsum types
  integer, private, parameter :: csSarvasGeneralised  = 6011
  
  ! Heatsum types
  integer, private, parameter :: hsDegreeDay = 6015
  integer, private, parameter :: hsDegreeHour  = 6016
  integer, private, parameter :: hsBioDay = 6017
  integer, private, parameter :: hsSigmoidPeriodUnits = 6018
  
  ! Pollen ripening types
  integer, private, parameter :: prCDLinear = 6021
  integer, private, parameter :: prHSLinear = 6022
  integer, private, parameter :: prCDNormal = 6023
  integer, private, parameter :: prHSNormal = 6024
  integer, private, parameter :: prCDNormDrn = 6025
  integer, private, parameter :: prCDGammaWTails = 6026
!  integer, private, parameter :: prDLNormal = 6023

  ! Pollen release types
  integer, private, parameter :: relExp = 6041
  integer, private, parameter :: relLim = 6042
  integer, private, parameter :: relInst = 6043
  
  ! Temperature response function
  integer, private, parameter :: trfBeta = 6051
  integer, private, parameter :: trfLinear = 6052

  
!!!!!!!!!!!!!!!!!!#define CAMS_DUMP

#ifdef CAMS_DUMP
  integer, private :: iDayInYear_, iCAMSdisp, iCAMSmet
  integer, private, parameter :: ixCAMSdisp  = 95, iyCAMSdisp = 40
  real :: now_sec_since_sunrise_, &
             & mdl_timestep_sec_, &
             & dayLength_hours_, &
             & T2m_, &
             & DailyTempr_, &
             & StartCDThr_, &
             & HS_in_, HS_out_, &
             & PollenLeft_in_, PollenLeft_out_, &
             & PollenRdyToFly_in_, PollenRdyToFly_out_, &
             & fMassInjected_out_
#endif


  
CONTAINS


  !**************************************************************************

  subroutine reserve_pollen_source(pollen_src, &   ! Src to initialise
                                 & iSrcNbr, &      ! Src number in the grand list
                                 & iSrcIdNbr)     ! SrcID number
    !
    ! Initialises the source:
    ! - stores the reference information: source number and source ID number
    ! - stores the total number of chemical descriptors that will be stored in the source
    ! - nullifies the source dynamic arrays 
    ! - set a few basic internal variables of the source
    !
    implicit none

    ! Imported parameters
    type(silam_pollen_source), intent(inout) :: pollen_src
    integer, intent(in) :: iSrcNbr, iSrcIdNbr
    !
    ! Nullify the basic variables
    !
    pollen_src%src_nm = ''
    pollen_src%sector_nm = ''
    !
    ! A bit of other stuff
    !
    nullify(pollen_src%fZDisp)
    nullify(pollen_src%levFractDispVert)
    !
    ! Main source parameters - enough to identify it in the global information list
    !
    pollen_src%src_nbr = iSrcNbr
    pollen_src%id_nbr = iSrcIdNbr
    !
    ! Finally, mark the source as incomplete
    !
    pollen_src%defined = silja_false

  end subroutine reserve_pollen_source

  
  !***********************************************************************

  subroutine fill_pollen_src_from_namelist(nlSetup, srcPollen, chDataDir)

    !
    ! Creates the set of rules, in accordance with which the computations are to
    ! be made. Input - SILAM namelist
    !
    implicit none

    !Imported parameters
    type(silam_pollen_source), intent(inout) :: srcPollen
    type(Tsilam_namelist), pointer :: nlSetup 
    character(len=*), intent(in) :: chDataDir

    ! Local variables
    integer :: nFiles, iFile, iTmp, nMaterials, iMode, iModePollen
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: pItems
    type(silam_sp) :: sp
    type(Taerosol) :: aerosoltmp
    type(silam_material), pointer :: materialTmp
    real :: fTmp

    if(.not.defined(nlSetup))then
      call set_error('Undefined namelist given','fill_pollen_src_from_namelist')
      return
    endif

    sp%sp => fu_work_string()
    nullify(pItems)
    
    srcPollen%defined = silja_false

    ! Source 
    srcPollen%src_nm = fu_content(nlSetup,'source_name')
    srcPollen%sector_nm = fu_content(nlSetup,'source_sector_name')
    
    ! input fields for the source features.
    srcPollen%nlInputFiles => fu_create_namelist('pollen_src_supplementary_files')
    if(error)return
    call get_items(nlSetup, 'supplementary_file', pItems, nFiles)
    if(error)return
    do iFile = 1, nFiles
      call add_namelist_item(srcPollen%nlInputFiles, 'supplementary_file', &
                           & fu_extend_grads_hat(fu_expand_environment(fu_content(pItems(iFile))), &
                                               & chDataDir))
    end do
!    ! source mask - no more. It is now a supplementary file
!    call add_namelist_item(srcPollen%nlInputFiles, 'source_area_mask', &
!                         & fu_extend_grads_hat(fu_expand_environment(fu_content(nlSetup,'source_area_mask')), &
!                                             & chIniFNm))
!    if(error)return

    ! Pollen species
    srcPollen%nSpecies = 0
    nullify(srcPollen%species)
    call set_aerosol(nlSetup, aerosolTmp)  ! Modes
    call get_items(nlSetup, 'emitted_material', pItems, nMaterials) ! Materials
    allocate(srcPollen%species(nMaterials))
    srcPollen%nSpecies = nMaterials
    if(error)return
    call msg('The following species are set for the pollen source:' + srcPollen%src_nm)
    iModePollen = -1
    do iTmp = 1, nMaterials
      sp%sp = fu_content(pItems(iTmp))
      read(unit=sp%sp,fmt=*) iMode
      sp%sp(1:) = sp%sp(index(sp%sp,' ')+1:)
      sp%sp = adjustl(sp%sp)
      materialTmp => fu_get_material_ptr(sp%sp)
      if(.not. associated(materialTmp))then
        call set_error('No pollen material found','fill_pollen_src_from_namelist')
        return
      endif
      if(defined(aerosolTmp) .and. iMode > 0)then
        call set_species(srcPollen%species(iTmp), materialTmp, fu_mode(aerosolTmp, iMode))
      elseif(iMode == 0)then
        call set_species(srcPollen%species(iTmp), materialTmp, in_gas_phase)
      else
        call set_error('If mode number is not 0 aerosol has to be given','fill_pollen_src_from_namelist')
        return
      endif
      !
      ! Checking the compatibility of emitted species:
      ! Only 1 aerosol mode is allowed for pollen, allergen in pollen and free allergen
      ! Other materials are allowed, pollen emission value will be emitted for them
      !
      select case(fu_pollen_material_type(materialTmp))
      case(pollenGrains)
        if(srcPollen%indPol > 0)then
          call set_error('Only one pollen grain species allowed for single source','fill_pollen_src_from_namelist') 
        else
          srcPollen%indPol = iTmp
        endif
      case(pollenAllergen) 
        if(srcPollen%indPolAlrg > 0)then
          call set_error('Only one pollen allergen species allowed for single source','fill_pollen_src_from_namelist') 
        else
          srcPollen%indPolAlrg = iTmp
        endif
      case(freeAllergen) 
        if(srcPollen%indFreeAlrg > 0)then
          call set_error('Only one free allergen species allowed for single source','fill_pollen_src_from_namelist') 
        else
          srcPollen%indFreeAlrg = iTmp
        endif
      if(error)return
      case default
        ! Other materials are allowed, pollen emission value will be emitted for them
      end select
      call report(srcPollen%species(iTmp))
    enddo
    if(error)return

    ! Checking the compatibility of emitted species:
    ! Allergen cannot be emitted without pollen
    !
    if(srcPollen%indPol < 0 .and. (srcPollen%indPolAlrg > 0 .or. srcPollen%indPolAlrg > 0))then
      call set_error('Allergen requested but no pollen','fill_pollen_src_from_namelist')
      return
    endif
    ! pollen grains and pollen allergen must have same mode
    if(srcPollen%indPolAlrg > 0)then
      if(.not. fu_mode(srcPollen%species(srcPollen%indPolAlrg)) == fu_mode(srcPollen%species(srcPollen%indPol)))then
        call set_error('Something strange with allergen emission modes','fill_pollen_src_from_namelist')
      endif
    endif

    ! Pollen amount
    !
    srcPollen%standardPollenTotal = fu_content_real(nlSetup,'climate_pollen_amt_per_m2')
    if(fu_fails(.not. (srcPollen%standardPollenTotal .eps. real_missing), &
              & 'Failed line climate_pollen_amt_per_m2','fill_pollen_src_from_namelist'))return

    call setNamedValue(fu_content(nlSetup,'uncertainty_of_total_pollen_amt'), 'fraction', &
                     & srcPollen%fUncertainty_tot_pollen_relat)
    if(fu_fails(.not. (srcPollen%fUncertainty_tot_pollen_relat .eps. real_missing), &
              & 'Failed line uncertainty_of_total_pollen_amt','fill_pollen_src_from_namelist'))return
    if(error)return

    ! Flowering thresholds
    sp%sp = fu_content(nlSetup,'flowering_thresholds')
    
    ! for start of flowering
    srcPollen%ifChillSum = index(sp%sp, 'chillsum') > 0
    srcPollen%ifStartHSThr = index(sp%sp, 'start_heatsum') > 0
    srcPollen%ifStartCDThr = index(sp%sp, 'start_calendar_day') > 0
    srcPollen%ifStartDLThr = index(sp%sp, 'start_daylength') > 0

    ! for end of flowering
    srcPollen%ifEndHSThr = index(sp%sp, 'end_heatsum') > 0
    srcPollen%ifEndCDThr = index(sp%sp, 'end_calendar_day') > 0
    srcPollen%ifEndDLThr = index(sp%sp, 'end_daylength') > 0

    ! Gamma season sets both start and end dates. 
    srcPollen%ifStartEndGammaThr = index(sp%sp, 'start_end_gamma_percentiles') > 0

    ! Meteo
    srcPollen%ifTempThr = index(sp%sp, 'instant_temperature') > 0
    srcPollen%ifDayTempThr = index(sp%sp, 'daily_temperature') > 0
    srcPollen%ifSWthr = index(sp%sp, 'soilwater') > 0
    
    ! dynamic plant size
    sp%sp = fu_content(nlSetup,'plant_growth')
    if(index(sp%sp, 'heatsum')>0)then
      if(.not. srcPollen%ifStartHSThr)then
        call set_error('heatsum dependent plant growth requires start heatsum threshold','fill_pollen_src_from_namelist')
        return
      endif
      srcPollen%ifHSGrowth = .true. 
    else
      srcPollen%ifHSGrowth = .false. 
    endif
    if(index(sp%sp, 'soilwater')>0)then
      srcPollen%SWGrowthMid = fu_content_real(nlSetup, 'soilwater_mid')
      srcPollen%SWGrowthSigma = fu_content_real(nlSetup, 'soilwater_sigma')
      srcPollen%SWGrowthDeath = fu_content_real(nlSetup, 'soilwater_death')
      srcPollen%ifSWGrowth = .true.
      if(fu_fails((.not. (srcPollen%SWGrowthMid .eps. real_missing)).and. &
                & (.not. (srcPollen%SWGrowthSigma .eps. real_missing)).and. &
                & (.not. (srcPollen%SWGrowthDeath .eps. real_missing)), &
       & 'soilwater_mid/sigma/death failed','fill_pollen_src_from_namelist'))return
    else
      srcPollen%ifSWGrowth = .false.
    endif
    
    ! Specific parameters for each threshold
    
    ! Calendar day and gamma
    if(srcPollen%ifStartCDThr .or. srcPollen%ifEndCDThr .or. srcPollen%ifStartEndGammaThr)then 
      ! flowering map shift
        sp%sp = fu_content(nlSetup,'flowering_map_shift')
        if(len_trim(sp%sp) > 0)then
          srcPollen%timeMapShift = fu_set_named_interval(sp%sp)
        else
          call msg('flowering_map_shift is absent from control file, take zero interval')
          srcPollen%timeMapShift = zero_interval
        endif
        if(error)then
          call set_error('Failed to set flowering_map_shift','fill_pollen_src_from_namelist')
          return
      endif
    endif

    ! Calendar day and daylength require uncertainty in days
    if(srcPollen%ifStartCDThr .or. srcPollen%ifEndCDThr .or. srcPollen%ifStartDLThr .or. srcPollen%ifEndDLThr)then 
      ! uncertainty of the threshold dates
      call setNamedValue(fu_content(nlSetup,'uncertainty_of_calendar_day_threshold_start'), 'day',  &
                         & srcPollen%fUncertainty_CD_days_start)
      if(fu_fails(.not. (srcPollen%fUncertainty_CD_days_start .eps. real_missing), &
        & 'failed line uncertainty_of_calendar_day_threshold_start','fill_pollen_src_from_namelist'))return
    endif
    !
    ! daylength
    !
    if(srcPollen%ifStartDLThr .or. srcPollen%ifEndDLThr)then 
      call setNamedValue(fu_content(nlSetup,'daylength_start'), 'hr',  &
                       & srcPollen%fDayLenStart_hrs)
      call setNamedValue(fu_content(nlSetup,'daylength_end'), 'hr',  &
                       & srcPollen%fDaylenEnd_hrs)
      if(fu_str_u_case(fu_content(nlSetup,'if_short_day_flowering')) =='YES')then
        srcPollen%ifShortDayTriggers = .true.
      elseif(fu_str_u_case(fu_content(nlSetup,'if_short_day_flowering')) =='NO')then
        srcPollen%ifShortDayTriggers = .false.
      else
        call set_error('Missing if_short_day_flowering = YES/NO','fill_pollen_src_from_namelist')
        return
      endif
      call setNamedValue(fu_content(nlSetup,'uncertainty_of_daylength_start'), 'hr',  &
                       & srcPollen%fUncertainty_DL_start_hrs)
      call setNamedValue(fu_content(nlSetup,'uncertainty_of_daylength_end'), 'hr',  &
                       & srcPollen%fUncertainty_DL_end_hrs)
      ! calendar day variables
      srcPollen%timeMapShift = zero_interval
    endif  ! daylength
    !
    ! Chillsum
    !
    if(srcPollen%ifChillSum)then
      select case(fu_str_l_case(fu_content(nlSetup,'chillsum_type')))
        case('sarvas_generalised')
          srcPollen%heatsum_params%iChillSumType = csSarvasGeneralised
          srcPollen%heatsum_params%chill_Tmax = fu_content_real(nlSetup, 'chill_Tmax_K')
          srcPollen%heatsum_params%chill_Topt = fu_content_real(nlSetup, 'chill_Toptimal_K')
          srcPollen%heatsum_params%chill_Tmin = fu_content_real(nlSetup, 'chill_Tmin_K')
          srcPollen%heatsum_params%chill_startDay = fu_content_real(nlSetup, 'chill_start_julian_day')
          srcPollen%heatsum_params%chill_ifNegativeAboveMax = &
                          & fu_str_u_case(fu_content(nlSetup, 'chill_if_negative_above_Tmax')) == 'YES'
        case default
          call set_error('Unknown chillsum type:' + fu_content(nlSetup,'chillsum_type'), &
                       & 'fill_pollen_src_from_namelist')
          return
      end select
    endif
    
    ! Heatsum
    !
    if(srcPollen%ifStartHSThr .or. srcPollen%ifEndHSThr)then ! T and day length functions for bioday accumulation
      select case(fu_str_l_case(fu_content(nlSetup,'heatsum_type')))
        case('degree_day')
          srcPollen%heatsum_params%iHeatSumType = hsDegreeDay
        case('degree_hour')
          srcPollen%heatsum_params%iHeatSumType = hsDegreeHour
        case('bioday')
          srcPollen%heatsum_params%iHeatSumType = hsBioDay
          select case(fu_str_l_case(fu_content(nlSetup,'temperature_response_function')))
            case('linear')
              srcPollen%iTRfunc = trfLinear
            case('beta')
              srcPollen%iTRfunc = trfBeta
            case default
              call set_error('Unknown temperature response function:' + fu_content(nlSetup,'temperature_response_function'), &
                       & 'fill_pollen_src_from_namelist')
              return
          end select
          srcPollen%loTemp = fu_content_real(nlSetup, 'loTemp')
          srcPollen%hiTemp = fu_content_real(nlSetup,'hiTemp')
          srcPollen%optTemp = fu_content_real(nlSetup,'optTemp')
          srcPollen%photoperiod = fu_content_real(nlSetup,'photoperiod')
          if(fu_fails((.not. (srcPollen%loTemp .eps. real_missing)).and. &
                    & (.not. (srcPollen%hiTemp .eps. real_missing)) .and. &
                    & (.not. (srcPollen%optTemp .eps. real_missing)) .and. &
                    & (.not. (srcPollen%photoperiod .eps. real_missing)), &
           & 'Failed loTemp/hiTemp/optTemp lines','fill_pollen_src_from_namelist'))return
        case('degree_day_sigmoid')
          srcPollen%heatsum_params%iHeatSumType = hsSigmoidPeriodUnits
          call msg('HS_midpoint_temp =', fu_content_real(nlSetup, 'HS_midpoint_temp'))
          call msg('HS_max_rate_per_day =', fu_content_real(nlSetup, 'HS_max_rate_per_day'))
          !
          ! A temporary solution: compute_HS_sigmoid_parameters2 is the new and better one
          ! but I already have plenty of source term files with the initial formulation.
          ! Will remove it anyway when finish the fitting
          !
          if(len_trim(fu_content(nlSetup,'HS_midpoint_temp')) > 0)then
            call compute_HS_sigmoid_parameters2(fu_content_real(nlSetup, 'HS_midpoint_temp'), &
                                              & fu_content_real(nlSetup, 'HS_cutoff_temp'), &
                                              & fu_content_real(nlSetup, 'HS_smoother'), &
                                              & fu_content_real(nlSetup, 'HS_smoother_power'), &
                                              & fu_content_real(nlSetup, 'HS_max_rate_per_day'), &
                                              & fu_content_real(nlSetup, 'HS_exp_temperature_response'), &
                                              & srcPollen%heatsum_params)
          else
            call compute_HS_sigmoid_parameters(fu_content_real(nlSetup, 'HS_saturation_temp'), &
                                             & fu_content_real(nlSetup, 'HS_max_rate_per_day'), &
                                             & fu_content_real(nlSetup, 'HS_exp_temperature_response'), &
                                             & srcPollen%heatsum_params)
          endif
        case default
          call set_error('Unknown heatsum type:' + fu_content(nlSetup,'heatsum_type'), &
                       & 'fill_pollen_src_from_namelist')
          return
      end select
      call setNamedValue(fu_content(nlSetup,'uncertainty_of_heat_sum_threshold_start'), 'fraction',  &
                       & srcPollen%fUncertainty_HS_relative_start)
      if(fu_fails(.not. (srcPollen%fUncertainty_HS_relative_start .eps. real_missing), &
                & 'Failed line uncertainty_of_heat_sum_threshold_start','fill_pollen_src_from_namelist'))return
    endif  ! if heatsum start-end-threshods
 
    ! Gamma season
    if(srcPollen%ifStartEndGammaThr)then
        srcPollen%gf_nTerms = fu_content_int(nlSetup,'gf_n_terms')
        if(fu_fails(srcPollen%gf_nTerms /= int_missing,'did not findnumber of Gamma-distribution terms', &
                                                                      & 'fill_pollen_src_from_namelist'))return
        srcPollen%gf_beta = fu_content_real(nlSetup,'gf_beta')
        if(fu_fails(.not. (srcPollen%gf_beta .eps. real_missing),'did not find beta for Gamma-distribution', &
                                                                      & 'fill_pollen_src_from_namelist'))return
        allocate(srcPollen%gf_powers(srcPollen%gf_nTerms), srcPollen%gf_scales(srcPollen%gf_nTerms), &
               & srcPollen%gf_start_times(srcPollen%gf_nTerms), stat = iTmp)
        if(fu_fails(iTmp==0,'Failed Gamma-distribution terms allocation','fill_pollen_src_from_namelist'))return
        sp%sp = fu_content(nlSetup,'gf_scales')
        read(unit=sp%sp,fmt=*,iostat=iMode) (srcPollen%gf_scales(iTmp),iTmp=1,srcPollen%gf_nTerms)
        if(fu_fails(iMode==0,'Failed reading gamma-distr. scales:'+sp%sp,'fill_pollen_src_from_namelist'))return
        sp%sp = fu_content(nlSetup,'gf_powers')
        read(unit=sp%sp,fmt=*,iostat=iMode) (srcPollen%gf_powers(iTmp),iTmp=1,srcPollen%gf_nTerms)
        if(fu_fails(iMode==0,'Failed reading gamma-distr. powers:'+sp%sp,'fill_pollen_src_from_namelist'))return
        sp%sp = fu_content(nlSetup,'gf_start_times')
        read(unit=sp%sp,fmt=*,iostat=iMode) (srcPollen%gf_start_times(iTmp),iTmp=1,srcPollen%gf_nTerms)
        if(fu_fails(iMode==0,'Failed reading gamma-distr. start times:'+sp%sp,'fill_pollen_src_from_namelist'))return
    endif
    
    if(error)return
    !
    ! Pollen ripening
    !
    select case(fu_str_l_case(fu_content(nlSetup,'ripening_method')))
      case('calday_linear')
        if(fu_fails(srcPollen%ifStartCDThr .and. srcPollen%ifEndCDThr, &
                  & 'For calday based ripening need to have both start and end calday thresholds', &
                  & 'fill_pollen_src_from_namelist'))return
        srcPollen%iRipeningType = prCDLinear

      case('calday_gamma_tails')
        if(fu_fails(srcPollen%ifStartEndGammaThr, &
                  & 'Gamma-type season ripening requires start_end_gamma_percentiles threshold', &
                  & 'fill_pollen_src_from_namelist'))return
        srcPollen%iRipeningType = prCDGammaWTails
      
      case('heatsum_linear')
        if(fu_fails(srcPollen%ifStartHSThr .and. srcPollen%ifEndHSThr, &
                  & 'For heatsum based ripening need to have both start and end heatsum thresholds', &
                  & 'fill_pollen_src_from_namelist'))return
        srcPollen%iRipeningType = prHSLinear
      
      case('calday_normal')
        if(fu_fails(srcPollen%ifStartCDThr .and. srcPollen%ifEndCDThr, &
                  & 'For calday based ripening need to have both start and end calday thresholds', &
                  & 'fill_pollen_src_from_namelist'))return
        srcPollen%iRipeningType = prCDNormal

      case('daylen_normal')
        if(fu_fails(srcPollen%ifStartDLThr .and. srcPollen%ifEndDLThr, &
                  & 'For daylen based ripening need to have both start and end daylen thresholds', &
                  & 'fill_pollen_src_from_namelist'))return
        srcPollen%iRipeningType = prCDNormal  ! we later turn daylaength into calendar day

      case('heatsum_normal')
        if(fu_fails(srcPollen%ifStartHSThr .and. srcPollen%ifEndHSThr, &
                  & 'For heatsum based ripening need to have both start and end heatsum thresholds', &
                  & 'fill_pollen_src_from_namelist'))return
        srcPollen%iRipeningType = prHSNormal 
        
        case('calday_normal_diurnal')
        if(fu_fails(srcPollen%ifStartCDThr .and. srcPollen%ifEndCDThr, &
                  & 'For calday based ripening need to have both start and end calday thresholds', &
                  & 'fill_pollen_src_from_namelist'))return
        srcPollen%iRipeningType = prCDNormDrn  
        call setNamedValue(fu_content(nlSetup,'diurnal_peak1'), 'sec', srcPollen%dayMid1)
        call setNamedValue(fu_content(nlSetup,'diurnal_peak2'), 'sec', srcPollen%dayMid2)
        call setNamedValue(fu_content(nlSetup,'diurnal_sigma1'), 'sec', srcPollen%daySgm1)
        call setNamedValue(fu_content(nlSetup,'diurnal_sigma2'), 'sec', srcPollen%daySgm2)
        call setNamedValue(fu_content(nlSetup,'diurnal_fraction1'), 'fraction', srcPollen%dayFrac1)         
        if(fu_fails((.not. (srcPollen%dayMid1 .eps. real_missing)) .and. &
                  & (.not. (srcPollen%dayMid2 .eps. real_missing)) .and. &
                  & (.not. (srcPollen%daySgm1 .eps. real_missing)) .and. &
                  & (.not. (srcPollen%daySgm2 .eps. real_missing)) .and. &
                  & (.not. (srcPollen%dayFrac1 .eps. real_missing)), &
         & 'Failed diurnal_peak1/2 or diurnal sigma1/2 of diurnal_fraction','fill_pollen_src_from_namelist'))return 
      case default
        call set_error('Unknown ripening_method:' + fu_content(nlSetup,'ripening_method'), &
                     & 'fill_pollen_src_from_namelist')
        return
    end select
    
    !
    ! For some types of the season curves, get the uncertainty parameters of start and end
    !
    if(srcPollen%iRipeningType == prCDLinear .or. &
     & srcPollen%iRipeningType == prHSLinear)then
      call set_fade_in_out_params(nlSetup, srcPollen%UncertaintyParams)
      if(error)return
    else 
      srcPollen%UncertaintyParams = fade_inout_missing
    endif
 
    ! Pollen release from buffer

    select case(fu_str_l_case(fu_content(nlSetup,'release_type')))
     case('exponential') 
       srcPollen%iReleaseType = relExp
       call setNamedValue(fu_content(nlSetup,'shortest_full_release_period'), 'sec',  &
                     & srcPollen%fBufReleaseLimit)
       if(fu_fails(.not. (srcPollen%fBufReleaseLimit .eps. real_missing), &
                 & 'Failed shortest_full_release_period','fill_pollen_src_from_namelist'))return
       srcPollen%fBufReleaseLimit = 1./srcPollen%fBufReleaseLimit
     case('limited')
       srcPollen%iReleaseType = relLim
       call setNamedValue(fu_content(nlSetup,'max_relative_release_rate'), 'fraction/sec',  &
                     & srcPollen%fBufReleaseLimit)
       if(fu_fails(.not. (srcPollen%fBufReleaseLimit .eps. real_missing), &
                 & 'Failed max_relative_release_rate','fill_pollen_src_from_namelist'))return
       srcPollen%fBufReleaseLimit = srcPollen%fBufReleaseLimit * srcPollen%standardPollenTotal
     case('instant')
       srcPollen%iReleaseType = relInst
     case default
       call set_error('release_type not given','fill_pollen_src_from_namelist')
       return
    end select
    
    call setNamedValue(fu_content(nlSetup,'low_humidity_threshold'), 'fraction', srcPollen%LowHumidThresh)
    if(fu_fails(.not. (srcPollen%LowHumidThresh .eps. real_missing), &
              & 'Failed low_humidity_threshold','fill_pollen_src_from_namelist'))return

    call setNamedValue(fu_content(nlSetup,'high_humidity_threshold'), 'fraction', srcPollen%HighHumidThresh)
    if(fu_fails(.not. (srcPollen%HighHumidThresh .eps. real_missing), &
              & 'Failed high_humidity_threshold','fill_pollen_src_from_namelist'))return

    call setNamedValue(fu_content(nlSetup,'precipitation_threshold'), 'mm/sec', srcPollen%precipThreshold)
    if(fu_fails(.not. (srcPollen%precipThreshold .eps. real_missing), &
              & 'Failed precipitation_threshold','fill_pollen_src_from_namelist'))return

    call setNamedValue(fu_content(nlSetup,'wind_speed_saturation_level'), 'm/sec', srcPollen%windSaturation)
    if(fu_fails(.not. (srcPollen%windSaturation .eps. real_missing), &
              & 'Failed wind_speed_saturation_level','fill_pollen_src_from_namelist'))return

    srcPollen%windMaxScale = fu_content_real(nlSetup,'wind_speed_max_impact')
    if(fu_fails(.not. (srcPollen%windMaxScale .eps. real_missing), &
              & 'Failed wind_speed_max_impact','fill_pollen_src_from_namelist'))return

    ! Emission height
    call setNamedValue(fu_content(nlSetup,'lowest_emission_injection_height'), 'm', &
                     & srcPollen%fEmissionBottom)
    if(fu_fails(.not. (srcPollen%fEmissionBottom .eps. real_missing), &
              & 'Failed lowest_emission_injection_height','fill_pollen_src_from_namelist'))return
    if(fu_content(nlSetup,'highest_emission_injection_height') == 'BLH')then
      srcPollen%ems2wholeABL = .true.
      srcPollen%fEmissionTop = real_missing
    else
      call setNamedValue(fu_content(nlSetup,'highest_emission_injection_height'), 'm', &
                       & srcPollen%fEmissionTop)
      if(fu_fails(.not. (srcPollen%fEmissionTop .eps. real_missing), &
                & 'Failed highest_emission_injection_height','fill_pollen_src_from_namelist'))return
    endif

    srcPollen%defined = silja_true   ! all done

    call free_work_array(sp%sp)
    call report(srcPollen)
    
  end subroutine fill_pollen_src_from_namelist


  !*****************************************************************

  subroutine add_input_needs_pollen_source(srcPollen, q_met_dynamic, q_met_st, &
                                                    & q_disp_dynamic, q_disp_st)
    !
    ! Returns input needs for the transformation routines. 
    !
    implicit none

    ! Imported parameters
    type(silam_pollen_source), intent(in) :: srcPollen
    integer, dimension(:), intent(inout) :: q_met_dynamic, q_met_st, &
                                          & q_disp_dynamic, q_disp_st
    ! Local variables
    integer :: iTmp

    ! Add needed quantities
    ! Meteo
    iTmp = fu_merge_integer_to_array(relative_humidity_2m_flag,        q_met_dynamic)
    iTmp = fu_merge_integer_to_array(windspeed_10m_flag,               q_met_dynamic)
    iTmp = fu_merge_integer_to_array(convective_velocity_scale_flag,   q_met_dynamic)
    iTmp = fu_merge_integer_to_array(total_precipitation_rate_flag,    q_met_st)
    if(srcPollen%ems2wholeABL)then
      iTmp = fu_merge_integer_to_array(abl_height_m_flag,              q_met_dynamic)
    endif
    ! pollen
    iTmp = fu_merge_integer_to_array(emission_mask_flag,               q_disp_st)
    iTmp = fu_merge_integer_to_array(pollen_correction_flag,           q_disp_st)
    iTmp = fu_merge_integer_to_array(pollen_left_relative_flag,        q_disp_st)
    iTmp = fu_merge_integer_to_array(pollen_total_per_m2_flag,         q_disp_st)
    iTmp = fu_merge_integer_to_array(pollen_rdy_to_fly_flag,           q_disp_st)
    
    ! allergen
    if(.not. srcPollen%indPolAlrg == int_missing)then
      iTmp = fu_merge_integer_to_array(pollen_potency_flag,           q_disp_st)
      iTmp = fu_merge_integer_to_array(allergen_rdy_to_fly_flag,      q_disp_st)
    endif
    
    ! fields depending on the emission computation principle
    if(srcPollen%ifStartHSThr .or. srcPollen%ifEndHSThr)then
        iTmp = fu_merge_integer_to_array(temperature_2m_flag,          q_met_dynamic)
        iTmp = fu_merge_integer_to_array(heatsum_flag,                 q_disp_st)
        iTmp = fu_merge_integer_to_array(growth_season_start_day_flag, q_disp_st)
        if(srcPollen%heatsum_params%iHeatSumType == hsDegreeDay .or. &
         & srcPollen%heatsum_params%iHeatSumType == hsBioDay)then
          iTmp = fu_merge_integer_to_array(day_mean_temperature_2m_flag, q_disp_st)
        endif
        if(srcPollen%heatsum_params%iHeatSumType == hsDegreeDay .or. &
         & srcPollen%heatsum_params%iHeatSumType == hsDegreeHour)then
          iTmp = fu_merge_integer_to_array(heatsum_cutoff_tempr_flag,    q_disp_st)
        endif
    endif
    if(srcPollen%ifStartHSThr)then
        iTmp = fu_merge_integer_to_array(start_heatsum_threshold_flag, q_disp_st)
    endif
    if(srcPollen%ifEndHSThr)then  
        iTmp = fu_merge_integer_to_array(end_heatsum_threshold_flag,   q_disp_st)
        if(srcPollen%ifStartHSThr)then 
            ! actually either end or length and start but request both here just in case ..
            iTmp = fu_merge_integer_to_array(heatsum_start_end_diff_flag,  q_disp_st)
        endif
    endif
    if(srcPollen%ifChillSum)then
        iTmp = fu_merge_integer_to_array(chillsum_flag, q_disp_st)
    endif
    if(srcPollen%ifStartCDThr .or. srcPollen%ifStartDLThr .or. srcPollen%ifStartEndGammaThr)then 
        iTmp = fu_merge_integer_to_array(start_calday_threshold_flag,  q_disp_st)
    endif
    if(srcPollen%ifEndCDThr .or. srcPollen%ifEndDLThr .or. srcPollen%ifStartEndGammaThr)then 
        iTmp = fu_merge_integer_to_array(end_calday_threshold_flag,    q_disp_st)
        iTmp = fu_merge_integer_to_array(calday_start_end_diff_flag,   q_disp_st)
    endif
    if(srcPollen%ifTempThr)then
        iTmp = fu_merge_integer_to_array(temperature_2m_flag,         q_met_dynamic)
        iTmp = fu_merge_integer_to_array(temperature_threshold_flag,  q_disp_st)
    endif
    if(srcPollen%ifDayTempThr)then
        iTmp = fu_merge_integer_to_array(daily_temp_threshold_flag,     q_disp_st)
        iTmp = fu_merge_integer_to_array(day_mean_temperature_2m_flag,  q_disp_st)
    endif
    
    if(srcPollen%ifSWthr)then
        iTmp = fu_merge_integer_to_array(soil_moisture_vol_frac_nwp_flag,   q_met_dynamic)
        iTmp = fu_merge_integer_to_array(soil_moisture_threshold_flag, q_disp_st)
    endif
    
    if(srcPollen%ifSWGrowth)then
        iTmp = fu_merge_integer_to_array(soil_moisture_vol_frac_nwp_flag,   q_met_dynamic)
        iTmp = fu_merge_integer_to_array(plant_growth_flag,            q_disp_st)
        iTmp = fu_merge_integer_to_array(growth_season_start_day_flag, q_disp_st)
    endif
     
  end subroutine add_input_needs_pollen_source


  !**********************************************************************

  subroutine init_emission_pollen(srcPollen, dispersionMarketPtr, start_time)
    !
    ! Actually creates the fields serving the pollen emission computations
    ! and drops them to the internal dispersion buffer
    !
    implicit none

    ! Imported parameters
    type(silam_pollen_source), intent(inout) :: srcPollen
    type(mini_market_of_stacks), pointer :: dispersionMarketPtr
    type(silja_time), intent(in) :: start_time

    ! Local variables
    type(silja_shopping_list) :: shop_list
    integer, dimension(:), pointer ::  q_disp_dyn, q_disp_stat, stack_quantities
    real, dimension(:), pointer :: pTmp1, pValues
    type(silja_field_id) :: id, idTmp
    logical :: ifOK
    type(silam_sp) :: strTmp
    integer :: iFlds, iTmp, nQuantities, ix, iy
    type(Tsilam_namelist), pointer :: nlPtr
    type(silja_field), pointer :: fieldPtr
    type(silam_vertical) :: vertTmp
    type(silja_field), pointer :: field
    real :: fSum
    
    q_disp_dyn => fu_work_int_array()
    q_disp_stat => fu_work_int_array()
    stack_quantities => fu_work_int_array()
    q_disp_dyn(1:max_quantities) = int_missing
    q_disp_stat(1:max_quantities) = int_missing
    stack_quantities(1:max_quantities) = int_missing
    strTmp%sp => fu_work_string()

    !
    ! Start from storing the quantities needed in the dispersion stack
    ! The meteo stack has already been requested
    !
    call add_input_needs_pollen_source(srcPollen, &
                                     & stack_quantities, stack_quantities, & ! meteo quantities, skip
                                     & q_disp_dyn, q_disp_stat)   ! dispersion-buffer quantities, use
    if(error)return
    !
    ! Make the shopping list with all needed quantities. Note that we can have several 
    ! pollen sources with the same or different species emitted - e.g. grass and birch pollen.
    ! The only way to allow them in one run is to use the species names for the dispersion stack 
    ! fields. If the species names are same for some sources, no problem, it is just the same source
    ! written in several parts.
    !
    call set_missing(shop_list)
    call set_missing(vertTmp, .true.)
    if(error)return

    do iFlds = 1, size(q_disp_stat)
      if(q_disp_stat(iFlds) == int_missing)exit
      call add_shopping_variable(shop_list, &
                               & q_disp_stat(iFlds), &     ! quantity
                               & fu_species_src(srcPollen, q_disp_stat(iFlds), .false.), &   !srcPollen%species(iSp), &
                               & grid_missing, &
                               & vertTmp, int_missing, &
                               & met_src_missing)
      if(error)return
    end do

    !
    ! The supplementary fields have to be taken from supplementary_info files. Note
    ! that due to species names used explicitly in the shopping list, they must be in the files.
    !
    call msg('Filling-in the pollen emission info from supplementary fields:' + srcPollen%src_nm)
    call fill_minimarket_from_namelist(dispersionMarketPtr, &
                                     & srcPollen%nlInputFiles, 'supplementary_file', & ! namelist and item
                                     & shop_list, start_time, &
                                     & static_climatology, &  ! ever valid
                                     & create_field, &            ! error if a clash
                                     & wdr_missing, &
                                     & dispersion_gridPtr, &
                                     & 5, .true., & ! iAccuracy, ifAdjustGrid
                                     & ifOK)
    !
    ! Find the source mask and store its pointer for the future use
    ! Check for existence (might be initialized and thus should not be overwritten) 
    !
    id = fu_set_field_id_simple(met_src_missing,&
                              & emission_mask_flag, &
                              & time_missing, &        ! valid time
                              & surface_level, &
                              & fu_species_src(srcPollen, emission_mask_flag, .true.))
    !
    ! Store the pointer to the source mask in dispersion grid
    !
    srcPollen%mask_disp_grid => fu_get_field_from_mm_general(dispersionMarketPtr, id, .false.)
    if(fu_fails(defined(srcPollen%mask_disp_grid),'Failed to find the source mask in stack','init_emission_pollen'))return
    
    ! Check that all the requested quantities are in the stack.
    ! Some of them might be initialized but if not, can still be added: 
    ! heatsum , pollen ready to fly and the amount of pollen still available for emission.
    ! Note that some of the quantities are specific only for heat-sum type of emission.
    
    call find_or_create_field(pollen_total_per_m2_flag, dispersionMarketPtr, ifOK, field, pValues)
    
    if(.not. ifOK)then  ! Not there yet

      if(srcPollen%ifHSGrowth .or. srcPollen%ifSWGrowth)then
        pValues(1:fs_dispersion) = 0.0
      
      else  ! No dynamic plant growth
        call find_field_from_stack(met_src_missing, &
                                 & pollen_correction_flag,&
                                 & time_missing,&
                                 & fu_stack(dispersionMarketPtr, 1),&
                                 & fieldPtr, &
                                 & ifOK, &
                                 & fu_species_src(srcPollen, pollen_correction_flag, .true.))
        if(ifOK)then
          pTmp1 => fu_grid_data(fieldPtr)
          if(error)return

          ! The amount of pollen per m2 is zero outside the source area and is a year-corrected 
          ! standard constant inside
          fSum = 0.0
          do iTmp = 1, fs_dispersion
              pValues(iTmp) = srcPollen%standardPollenTotal * pTmp1(iTmp)
              fSum = fSum + srcPollen%standardPollenTotal * pTmp1(iTmp)
          end do
          call msg('Total pollen without source mask:', fSum*fu_cell_size(dispersion_grid,iTmp))
        else
          call set_error('No pollen-amount correction field in dispersion stack','init_emission_pollen')
          return
        endif
      endif   ! if HSGrowth or SWGrowth
    endif  ! pollen_total_per_m2_flag is in market
 
    ! Adding the left pollen 
    call find_or_create_field(pollen_left_relative_flag, dispersionMarketPtr, ifOK, field, pValues, 1.0)
    
!    ArTmp(1:fs_dispersion) = 0.0
    ! Temperature-related stuff is to be taken only for heat/chill-sum emission type
    if(srcPollen%ifStartHSThr .or. srcPollen%ifEndHSThr)then
      call find_or_create_field(heatsum_flag, dispersionMarketPtr, ifOK, field, pValues, 0.0)
    endif
    
    ! Chill sum, make and zero it
    if(srcPollen%ifChillSum)then
      call find_or_create_field(chillsum_flag, dispersionMarketPtr, ifOK, field, pValues, 0.0)
    endif
    
    ! Plant growth field
    if(srcPollen%ifSWGrowth)then
      call find_or_create_field(plant_growth_flag, dispersionMarketPtr, ifOK, field, pValues, 0.0)
    endif
    
    ! Adding the pollen ready to fly amount 
    call find_or_create_field(pollen_rdy_to_fly_flag, dispersionMarketPtr, ifOK, field, pValues, 0.0)

    ! Adding the allergen ready to fly amount 
    if(.not. srcPollen%indPolAlrg == int_missing)then
      call find_or_create_field(allergen_rdy_to_fly_flag, dispersionMarketPtr, ifOK, field, pValues, 0.0)
    endif
     
    ! For backward compatibility, season end thresholds (CD, HS) can be computed from start and length 
    ! Alternatively, season length from start and end
    ! To avoid complex logic analoguous to meteo input, make both in any case
    if(srcPollen%ifEndHSThr)then
      call find_or_create_field(end_heatsum_threshold_flag, dispersionMarketPtr, ifOK, field, pValues)
  
      if(.not. ifOK)then  ! make it from start and length

        call find_field_from_stack(met_src_missing, &
                                 & start_heatsum_threshold_flag,&
                                 & time_missing,&
                                 & fu_stack(dispersionMarketPtr, 1),&
                                 & fieldPtr, &
                                 & ifOK, &
                                 & fu_species_src(srcPollen, start_heatsum_threshold_flag, .true.))
        if(ifOK)then
          pTmp1 => fu_grid_data(fieldPtr)
          if(error)return
          pValues(1:fs_dispersion) = pTmp1(1:fs_dispersion)
        else
          call set_error('No start heatsum field in dispersion stack','init_emission_pollen')
          return
        endif
        call find_field_from_stack(met_src_missing, &
                                 & heatsum_start_end_diff_flag,&
                                 & time_missing,&
                                 & fu_stack(dispersionMarketPtr, 1),&
                                 & fieldPtr, &
                                 & ifOK, &
                                 & fu_species_src(srcPollen, heatsum_start_end_diff_flag, .true.))
        if(ifOK)then
          pTmp1 => fu_grid_data(fieldPtr)
          if(error)return
          do iTmp = 1, fs_dispersion
            pValues(iTmp) = pValues(iTmp) + pTmp1(iTmp)
          end do
        else
          call set_error('No heatsum diff field in dispersion stack','init_emission_pollen')
          return
        endif

      else
        ! end heatsum exists, make duration
        if(srcPollen%ifStartHSThr)then ! Make the season length (no reason apart from filling up the stack)
          call find_or_create_field(heatsum_start_end_diff_flag, dispersionMarketPtr, ifOK, field, pValues)

          if(ifOK)then
            ! If exists, weird, as all start, end and length are given and might not be consistent
            call set_error('Start end and length are given and might not be consistent','init_emission_pollen')
            return
          else

            call find_field_from_stack(met_src_missing, &
                                     & start_heatsum_threshold_flag,&
                                     & time_missing,&
                                     & fu_stack(dispersionMarketPtr, 1),&
                                     & fieldPtr, &
                                     & ifOK, &
                                     & fu_species_src(srcPollen, start_heatsum_threshold_flag, .true.))
            if(ifOK)then
              pTmp1 => fu_grid_data(fieldPtr)
              if(error)return
              pValues(1:fs_dispersion) = pTmp1(1:fs_dispersion)
            else
              call set_error('No start heatsum field in dispersion stack','init_emission_pollen')
              return
            endif
            call find_field_from_stack(met_src_missing, &
                                     & end_heatsum_threshold_flag,&
                                     & time_missing,&
                                     & fu_stack(dispersionMarketPtr, 1),&
                                     & fieldPtr, &
                                     & ifOK, &
                                     & fu_species_src(srcPollen, end_heatsum_threshold_flag, .true.))
            if(ifOK)then
              pTmp1 => fu_grid_data(fieldPtr)
              if(error)return
              do iTmp = 1, fs_dispersion
                pValues(iTmp) = pTmp1(iTmp) - pValues(iTmp)
              end do
            else
              call set_error('No end heatsum field in dispersion stack','init_emission_pollen')
              return
            endif

          endif   ! if duration exists
        endif  ! if start HS threshold
      endif ! end HS exists
    endif  ! if end HS threshold
    
    ! End of season might need to be created still
    if(srcPollen%ifEndCDThr)then
      call find_or_create_field(end_calday_threshold_flag, dispersionMarketPtr, ifOK, field, pValues)

      if(.not. ifOK)then  ! make it from start and length

        call find_field_from_stack(met_src_missing, &
                                 & start_calday_threshold_flag,&
                                 & time_missing,&
                                 & fu_stack(dispersionMarketPtr, 1),&
                                 & fieldPtr, &
                                 & ifOK, &
                                 & fu_species_src(srcPollen, start_calday_threshold_flag, .true.))
        if(ifOK)then
          pTmp1 => fu_grid_data(fieldPtr)
          if(error)return
          pValues(1:fs_dispersion) = pTmp1(1:fs_dispersion)
        else
          call set_error('No start calday field in dispersion stack','init_emission_pollen')
          return
        endif
        call find_field_from_stack(met_src_missing, &
                                 & calday_start_end_diff_flag,&
                                 & time_missing,&
                                 & fu_stack(dispersionMarketPtr, 1),&
                                 & fieldPtr, &
                                 & ifOK, &
                                 & fu_species_src(srcPollen, calday_start_end_diff_flag, .true.))
        if(ifOK)then
          pTmp1 => fu_grid_data(fieldPtr)
          if(error)return
          do iTmp = 1, fs_dispersion
            pValues(iTmp) = pValues(iTmp) + pTmp1(iTmp)
          end do
        else
          call set_error('No calday diff field in dispersion stack','init_emission_pollen')
          return
        endif
      endif  ! EndCDThr found
    endif  ! ifEndCDThr
 
    ! Day length number is converted here to calendar day map
    if(srcPollen%ifStartDLThr)then
      ! Start day field
      call find_or_create_field(start_calday_threshold_flag, dispersionMarketPtr, ifOK, field, pValues)
      if(ifOK)then
        ! If exists, weird, day length and calendar day must never be together
        call set_error('Start day length and calendar day cannot be both present','init_emission_pollen')
        return
      else
        ! convert the day length to canedar day
        do iy = 1, ny_dispersion
          do ix = 1, nx_dispersion
            pValues(ix + (iy -1) * nx_dispersion) = fu_get_day_from_daylen( &
                                    & fu_lat_geographical_from_grid(real(ix),real(iy),dispersion_grid), &
                                    & srcPollen%fDayLenStart_hrs, srcPollen%ifShortDayTriggers)
          end do
        end do  
      endif  ! if start day field present
      !
!call msg('>>>>>>> Dump of CD for the daylength,', srcPollen%fDayLenStart_hrs, srcPollen%fDayLenEnd_hrs)
!iTmp = open_gradsfile_o('', &  ! directory not used 
!                                   & 'CD_maps_from_daylen.grads', fu_grid(fu_id(field)))
!call write_next_field_to_gradsfile(iTmp, fu_id(field), pValues(1:nx_dispersion*ny_dispersion))
      !
      ! End day field
      call find_or_create_field(end_calday_threshold_flag, dispersionMarketPtr, ifOK, field, pValues)
      if(ifOK)then
        ! If exists, weird, day length and calendar day must never be together
        call set_error('End day length and calendar day cannot be both present','init_emission_pollen')
        return
      else
        ! convert the day length to canedar day
        do iy = 1, ny_dispersion
          do ix = 1, nx_dispersion
            pValues(ix + (iy -1) * nx_dispersion) = fu_get_day_from_daylen( &
                                    & fu_lat_geographical_from_grid(real(ix),real(iy),dispersion_grid), &
                                    & srcPollen%fDayLenEnd_hrs, srcPollen%ifShortDayTriggers)
          end do
        end do  
      endif  ! if end day field present
      !
!call write_next_field_to_gradsfile(iTmp, fu_id(field), pValues(1:nx_dispersion*ny_dispersion))
!call close_gradsfile_o(iTmp,"")
      !
      ! Now, switch daylength to calendar day
      srcPollen%ifStartDLThr = .false.
      srcPollen%ifEndDLThr = .false.
      srcPollen%ifStartCDThr = .true.
      srcPollen%ifEndCDThr = .true.
    endif  ! day length
    !
    ! If CD thresholds are involved, make/check the season length. Not really needed: start and end 
    ! are forced above but procedure below uses all three fields
    !
    if((srcPollen%ifStartCDThr .and. srcPollen%ifEndCDThr) .or. srcPollen%ifStartEndGammaThr)then
      ! Make the season length
      call find_or_create_field(calday_start_end_diff_flag, dispersionMarketPtr, ifOK, field, pValues)
      if(.not. ifOK)then
        ! make it
        call find_field_from_stack(met_src_missing, &
                                 & start_calday_threshold_flag,&
                                 & time_missing,&
                                 & fu_stack(dispersionMarketPtr, 1),&
                                 & fieldPtr, &
                                 & ifOK, &
                                 & fu_species_src(srcPollen, start_calday_threshold_flag, .true.))
        if(ifOK)then
          pTmp1 => fu_grid_data(fieldPtr)
          if(error)return
          pValues(1:fs_dispersion) = pTmp1(1:fs_dispersion)
        else
          call set_error('No start calday field in dispersion stack','init_emission_pollen')
          return
        endif
        call find_field_from_stack(met_src_missing, &
                                 & end_calday_threshold_flag,&
                                 & time_missing,&
                                 & fu_stack(dispersionMarketPtr, 1),&
                                 & fieldPtr, &
                                 & ifOK, &
                                 & fu_species_src(srcPollen, end_calday_threshold_flag, .true.))
        if(ifOK)then
          pTmp1 => fu_grid_data(fieldPtr)
          if(error)return
          do iTmp = 1, fs_dispersion
            pvalues(iTmp) = pTmp1(iTmp) - pValues(iTmp)
          end do
        else
          call set_error('No end calday field in dispersion stack','init_emission_pollen')
          return
        endif
      endif   ! if found 
    endif  ! if field added/created

    ! Now check all
    ! dynamic dispersion stack
    CALL supermarket_2d_quantities(dispersionMarketPtr, &
                                 & met_src_missing, multi_time_stack_flag, &
                                 & stack_quantities, nQuantities)
    if(error)return
    do iFlds = 1, size(q_disp_dyn)
      if(q_disp_dyn(iFlds) == int_missing) exit
      if(.not. fu_quantity_in_quantities(q_disp_dyn(iFlds), stack_quantities))then
        call set_error('Missing dynamic quantity:' + fu_quantity_string(q_disp_dyn(iFlds)), &
                     & 'init_emission_pollen')
      endif
    end do ! checking the quantities
    if(error)return
    
    ! static dispersion stack
    CALL supermarket_2d_quantities(dispersionMarketPtr, &
                                 & met_src_missing, single_time_stack_flag, &
                                 & stack_quantities, nQuantities)
    call msg('Quantities in the static dispersion stack')
    do iFlds = 1, nQuantities
      if(stack_quantities(iFlds) == int_missing)exit
      call msg(fu_quantity_string(stack_quantities(iFlds)))
    enddo
    call msg('')
    do iFlds = 1, size(q_disp_stat)
      if(q_disp_stat(iFlds) == int_missing) exit
      if(.not. fu_quantity_in_quantities(q_disp_stat(iFlds), stack_quantities))then
        ! Except season lengths (HS, CD) as they are not necessary, just might be used for computing 
        ! the end threshold if it is not given
        if(q_disp_stat(iFlds) == heatsum_start_end_diff_flag .or. &
                                      & q_disp_stat(iFlds) == calday_start_end_diff_flag)cycle
        call set_error('Missing static quantity:' + fu_quantity_string(q_disp_stat(iFlds)), &
                     & 'init_emission_pollen')
      endif
    end do ! checking the quantities
    if(error)return

    ! Since the threshold may be computed for sub-areas of the domain, clean the field
    ifOK = .true.
    if(srcPollen%ifStartCDThr .or. srcPollen%ifStartEndGammaThr)then
      if(fu_fails(clean(start_calday_threshold_flag),'Failed cleaning start_calday_threshold_flag','init_emission_pollen'))return
    endif
    if(srcPollen%ifStartHSThr)then
      if(fu_fails(clean(start_heatsum_threshold_flag),'Failed cleaning start_heatsum_threshold_flag','init_emission_pollen'))return
    endif
    if(srcPollen%ifEndCDThr .or. srcPollen%ifStartEndGammaThr)then
      if(fu_fails(clean(end_calday_threshold_flag),'Failed cleaning end_calday_threshold_flag','init_emission_pollen'))return
    endif
    if(srcPollen%ifEndHSThr)then
      if(fu_fails(clean(end_heatsum_threshold_flag),'Failed cleaning end_heatsum_threshold_flag','init_emission_pollen'))return
    endif
    if(srcPollen%ifTempThr)then
      if(fu_fails(clean(temperature_threshold_flag),'Failed cleaning temperature_threshold_flag','init_emission_pollen'))return
    endif
    if(srcPollen%ifDayTempThr)then
      if(fu_fails(clean(daily_temp_threshold_flag),'Failed cleaning daily_temp_threshold_flag','init_emission_pollen'))return
    endif
    if(srcPollen%ifSWthr)then
      if(fu_fails(clean(soil_moisture_threshold_flag),'Failed cleaning soil_moisture_threshold_flag','init_emission_pollen'))return
    endif
    
    call free_work_array(strTmp%sp)
    call free_work_array(stack_quantities)
    call free_work_array(q_disp_dyn)
    call free_work_array(q_disp_stat)

#ifdef CAMS_DUMP
open(55,file = 'd:\model\silam_v5_5\ragweed_src_extract\output\src_dump.txt',mode="WRITE")
write(55,'(A)') "## year mon day hour, min, iDayInYear_  now_sec_since_sunrise_ mdl_timestep_sec_ dayLength_hours_  T2m_  DailyTempr_ StartCDThr_ HS_in_ HS_out_ PollenLeft_in_ PollenLeft_out_ PollenRdyToFly_in_ PollenRdyToFly_out_ fMassInjected_out_"
#endif

  contains
  
    !=================================================================================
    
    subroutine find_or_create_field(quantity, ptrMarket, ifFound, pField, pVals, fill_value)
      ! 
      ! Adding the total pollen amount as the climatologic amount and year-specific correction
    
      implicit none 
      
      ! Imported parameters
      integer :: quantity
      type(mini_market_of_stacks), pointer :: ptrMarket
      logical, intent(out) :: ifFound
      type(silja_field), pointer :: pField
      real, dimension(:), pointer :: pVals
      real, intent(in), optional :: fill_value
      
      ! Local variables
      type(silja_field_id) :: id
      
      ! Check for existence (might be initialized and thus should not be overwritten) 
      id = fu_set_field_id_simple(met_src_missing,&
                                & quantity, &
                                & time_missing, &        ! valid time
                                & surface_level, &
                                & fu_species_src(srcPollen, quantity, .true.))
      pField => fu_get_field_from_mm_general(ptrMarket, id, .false.)
      if(associated(pField))then
        ifFound= defined(pField)    ! exists already
      else
        ifFound = .false.   ! not found in the market
        !
        ! If fill_value is given, make it up!
        !
        id = fu_set_field_id(met_src_missing,&
                           & quantity, &
                           & start_time, &        ! analysis time
                           & zero_interval, &     ! forecast length
                           & dispersion_grid,&    ! grid
                           & surface_level, &     ! level
                           & zero_interval, &     ! length of accumulation
                           & zero_interval, &     ! length of validity
                           & accumulated_flag, &  ! field_kind
                           & species = fu_species_src(srcPollen, quantity, .true.)) ! species
        call find_field_data_storage_2d(ptrMarket, id, single_time_stack_flag, pVals)
        if(error)return

        pField => fu_get_field_from_mm_general(ptrMarket, id, .false.)
        
        if(present(fill_value)) pVals(:) = fill_value          

      endif   ! if field is in market

    end subroutine find_or_create_field
    
    !================================================================================
    
    logical function clean(flag)
      implicit none
      integer, intent(in) :: flag
      logical :: ifOK
      
      clean = .false.
      call find_field_from_stack(met_src_missing, &
                               & flag,&
                               & time_missing,&
                               & fu_stack(dispersionMarketPtr, 1),&
                               & fieldPtr, &
                               & ifOK)
      if(.not. ifOK)then
        call set_error('Failed to find field in stack','init_emission_pollen')
        call msg('quantity missing:', flag)
        return
      else
        pTmp1 => fu_grid_data(fieldPtr)
        do iFlds = 1, fs_dispersion
          if(pTmp1(iFlds) < 0.) pTmp1(iFlds) = real_missing
        end do
      endif
      clean = .true.
    end function clean
  
  end subroutine init_emission_pollen


  !************************************************************************

  function fu_species_src(pollen_src, quantity, ifSpeciesMandatory)result(species_src)
    !
    ! Selects the proper name of the substance for the given quantity - just to be able to 
    ! choose the right input and dispersion-stack fields for each source.
    ! Indices of main pollen species, pollen and free allergen are known by the source
    !
    implicit none

    ! Return value
    type(silam_species) :: species_src
    
    ! Imported parameters
    type(silam_pollen_source), intent(in) :: pollen_src
    integer, intent(in) :: quantity
    logical, intent(in) :: ifSpeciesMandatory

    select case(quantity)
      
      ! We have two types of quantities so far: those, which are universal for all pollen sources,
      ! and those, which are related to the pollen or allergen release of the specific taxon. 
      case (pollen_rdy_to_fly_flag, pollen_correction_flag, &
          & pollen_left_relative_flag, pollen_total_per_m2_flag, &
          & start_calday_threshold_flag, end_calday_threshold_flag, &
          & start_heatsum_threshold_flag, end_heatsum_threshold_flag, &
          & temperature_threshold_flag, daily_temp_threshold_flag, soil_moisture_threshold_flag, &
          & heatsum_flag, chillsum_flag, growth_season_start_day_flag, heatsum_cutoff_tempr_flag, &
          & emission_mask_flag, & !fraction_of_land_flag, &
          & pollen_potency_flag, &
          & heatsum_start_end_diff_flag, calday_start_end_diff_flag, plant_growth_flag)
        species_src = pollen_src%species(pollen_src%indPol)
        
      case(allergen_rdy_to_fly_flag)
        species_src = pollen_src%species(pollen_src%indPolAlrg)

      case default
        species_src = species_missing
        
    end select
      
    if(ifSpeciesMandatory)then
      if(species_src == species_missing)then
        call set_error('Species are required for quantity:>>' + fu_quantity_string(quantity)+ '<<', &
                     & 'fu_species_src')
      endif
    endif
                 
  end function fu_species_src


  !****************************************************************************

  subroutine add_source_species_pollen_src(pollen_src, species_list, nSpecies)
    !
    ! Get the species emitted by this source and make a list of
    ! them. The source must be initialized, of course.
    !
    implicit none
    type(silam_pollen_source), intent(in) :: pollen_src
    type(silam_species), dimension(:), pointer :: species_list
    integer, intent(inout) :: nSpecies

    call addSpecies(species_list, nSpecies, pollen_src%species, pollen_src%nSpecies)

  end subroutine add_source_species_pollen_src


  !**************************************************************************

  subroutine link_pollen_src_to_species(species_list, pollen_src)
    !
    ! Having the cocktails created, we should establish the shortcut links between the 
    ! source descriptors and species list.
    ! Note that cocktail is "anything" and the only requirement is that it has all the species
    ! existing in the source inventory.
    !
    implicit none

    ! Imported parameters
    type(silam_pollen_source), intent(inout) :: pollen_src
    type(silam_species), dimension(:), pointer :: species_list
    !
    ! Linkage is actually just creation of the chemical adaptor
    !
    call create_adaptor(pollen_src%species, species_list, pollen_src%adaptor)
    
  end subroutine link_pollen_src_to_species


  !**********************************************************************
  
  subroutine create_src_cont_grd_pollen_src(pollen_src, grid_template, ifVerbose, ifMinimal, ifExtended)
    !
    ! Creates the grid that covers the area with active pollen emission: source mask
    !
    implicit none

    ! Imported parameters
    type(silam_pollen_source), intent(in) :: pollen_src
    type(silja_grid), intent(inout) :: grid_template
    logical, intent(in) :: ifVerbose, ifMinimal
    logical, intent(out) :: ifExtended

    ! Local variables
    integer, dimension(:,:), pointer :: arFlag
    integer :: nx, ny, ix_start, iy_start, ix_end, iy_end

    ifExtended = .false.

    ! If template is underfined, source mask is the right grid, if defined - it can be cut
    !
    if(defined(grid_template))then
      !
      ! Try to cut the grid_template leaving only the area covered by the source
      !
      if(ifMinimal)then
        ix_start = 1
        iy_start = 1
        call grid_dimensions(grid_template, ix_end, iy_end)
        if(error)return

        arFlag => fu_work_int_array_2d()
        if(error)return
        arFlag(ix_start:ix_end, iy_start:iy_end) = -1
        !
        ! Cut non-overlapping parts of the grid_template
        !
        call fill_overlap_flag_array(grid_template, &             ! grid to cut
                                   & fu_grid(pollen_src%mask_src_grid), & ! grid area to preserve
                                   & arFlag, &                    ! overlap flag array
                                   & cover_the_grid_area, 2)      ! filling algorithm, overlap mark
        if(error)then
          call free_work_array(arFlag)
          call unset_error('create_src_cont_grd_a_src')
          return
        endif

        call cut_empty_lines(arFlag, ix_start, iy_start, ix_end, iy_end, &
                           & -1, 1, 2)  !value_to_cut, value_to_stay, value_to_preserve)
        if(error)return

        call cut_grid_size(grid_template, ix_start, iy_start, ix_end, iy_end)

        call free_work_array(arFlag)
      endif  ! ifMinimal
    else
      !
      ! Undefined grid_template. Just return the source grid
      !
      grid_template = fu_grid(pollen_src%mask_src_grid)

    endif  ! if defined grid_template
    
  end subroutine create_src_cont_grd_pollen_src


  !*****************************************************************

  subroutine project_pollen_src_second_grd(pollen_src, grid, vert_disp, vert_proj, iAccuracy)
    !
    ! Make the horizontal and vertical reprojection
    ! Horizontal part stolen from area source
    ! 
    implicit none

    ! Imported parameters
    type(silam_pollen_source), intent(inout) :: pollen_src
    type(silja_grid), intent(in) :: grid
    type(silam_vertical), intent(in) :: vert_disp, vert_proj
    integer, intent(in) :: iAccuracy

    ! Local variables
    integer :: i,j, ixSrc, iySrc, iSrc, iNew, iDisp, nSmall, nxSrc, nySrc, nxNew, nyNew, nSrc, nDisp
    real :: fxNew, fyNew, val, fArea
    type(silam_vertical) :: vertTmp
    real, dimension(:), pointer ::pMaskSrcGrd, pMaskDispGrd
    
    ! Now re-project the vertical of the source to the given vertical
    pollen_src%vertLevsDispVert = vert_disp
    pollen_src%nLevsDispVert = fu_NbrOfLevels(vert_disp)
    if(.not. pollen_src%ems2wholeABL)then
      allocate(pollen_src%levFractDispVert(pollen_src%nLevsDispVert), &
             & pollen_src%fzDisp(pollen_src%nLevsDispVert), stat=i)
      if(i /= 0)then
        call set_error('Failed to allocate dispersion-vertical level fractions','project_pollen_src_second_grd')
        return
      endif
      pollen_src%levFractDispVert(:) = 0.0
      pollen_src%fzDisp(:) = 0.0

      ! Create the vertical, which can depend on the emission method
 
      call set_vertical(fu_set_layer_between_two(layer_btw_2_height, pollen_src%fEmissionBottom, &
                                                                 & pollen_src%fEmissionTop), vertTmp)
      if(error)return

      call reproject_verticals(vertTmp, (/1.0/), &                   ! vertical from, fractions from
                             & vert_proj, pollen_src%levFractDispVert, &   ! vertical to, fractions to
                             ! mass centres, number of non-zero levels:
                             & pollen_src%fzDisp, pollen_src%nLevsDispVert, &
                             & ifMassCentreInRelUnit=.true.)
      call msg("pollen_src%levFractDispVert",pollen_src%levFractDispVert)
      call set_missing(vertTmp, .false.)
    else
         nullify(pollen_src%levFractDispVert, pollen_src%fzDisp)
    endif

  end subroutine project_pollen_src_second_grd


  !**************************************************************************
  !**************************************************************************
  !**************************************************************************

  subroutine compute_emission_for_pollen(srcPollen, &
                                       & met_buf, disp_buf, & 
                                       & now, timestep, &            
                                       & pHorizInterpMet2DispStruct, ifHorizInterp, &
                                       & ifSpeciesMoment, &
                                       & emisMap, mapCoordX, mapCoordY, mapCoordZ, &
                                       & fMassInjected, &
                                       & arSourceIdMapping, ifUseSourceIdMapping)
    !
    ! Computes the intermediate heat sums and final emission fields for pollen.
    ! Note: computations are done regardless the land use and forest maps.
    ! Corresponding map is supposed to come as an emission map itself, so that
    ! it determines the forest fraction, while here we compute the mass flux 
    ! (a number of grains) emitted per unit forest area.
    ! Consequently, the computations are tightly connected to meteorological rather
    ! than to dispersion buffer.
    !
    implicit none

    ! Imported parameters
    type(Tfield_buffer), pointer ::  met_buf, disp_buf  ! meteo and internal field buffers
    type(silam_pollen_source), target, intent(in) :: srcPollen  ! Also driving rules for pollen
    type(silja_time), intent(in) :: now                 ! current time
    type(silja_interval), intent(in) :: timestep        ! model time step
    type(THorizInterpStruct), intent(in) :: pHorizInterpMet2DispStruct
    logical, intent(in) :: ifHorizInterp, ifSpeciesMoment
    type(Tmass_map), intent(inout) :: emisMap, mapCoordX, mapCoordY, mapCoordZ  ! pointer to dynamic emission masses
    real(8), dimension(:), intent(inout) :: fMassInjected
    integer, dimension(:,:), intent(in) :: arSourceIdMapping   ! if we want different source indices over the grid
    logical, intent(in) :: ifUseSourceIdMapping             ! if we want different source indices over the grid
    ! Local variables
    integer :: iMeteo, iDisp, ix, iy, iSp, iLev, iSrcNbr, &
             & indStartCDThr,indEndCDThr,indStartHSThr,indEndHSThr,indTempThr,indDayTempThr,indSoilWaterThr, &
             & iDayInYear, indT2m, indDailyTempr, indSoilWater, &
             & indHS, indStartDay, indTemprCutOff, indCS, &
             & indHumid, indWind10m, indConvVelocity, indPrecip, &
             & indAnnualTotalCorr, indPollenLeft, indPollenRdyToFly, &
             & indAlrgRdyToFly, indPotency, indBLH, indCDDiff, indHSDiff, indPlantGrowth, indPollenTotal
    real :: fTmp, fHS, fThresh, fCellTotal, fPolAlrg, emsAmntLev, rdyPollen, rdyAllergen, leftPollen, &
          & overlapBottom, overlapTop, dHSdt, fzDisp, levFractDispVert, check, timestep_sec
    real, dimension(:), pointer :: xSize, ySize, fLandFractionMet, pSrcMask
    real, dimension(max_species) :: emsAmnt
    type(silja_interval) :: HS_aver_period
    logical :: ifEms, ifSeasonStarted
    type(silam_vertical) :: vertTmp


    fLandFractionMet => fu_grid_data(fraction_of_land_fld)
    pSrcMask => fu_grid_data(srcPollen%mask_disp_grid)
    xSize => fu_grid_data(dispersion_cell_x_size_fld)
    ySize => fu_grid_data(dispersion_cell_y_size_fld)
    emsAmnt(1:srcPollen%nSpecies) = 0.0
    rdyPollen = 0.0
    rdyAllergen = 0.0
    leftPollen = 0.0
    timestep_sec = fu_sec(timestep)
    
    if(error)return

    if(.not. emisMap%gridTemplate == dispersion_grid)then
      call msg('Emission map grid:')
      call report(emisMap%gridTemplate)
      call msg('Dispersion grid:')
      call report(dispersion_grid)
      call set_error('Grid of pollen emission map differs from dispersion grid', &
                   & 'compute_emission_for_pollen')
      return
    endif
    
    ! General fields, always needed no matter what method is used
    ! meteo (relative humidity(3D), windspeed, precipitation, convective velocity scale)
    indHumid =       fu_index(met_buf, relative_humidity_2m_flag, .true.)
    indWind10m =     fu_index(met_buf, windspeed_10m_flag, .true.)
    indConvVelocity =fu_index(met_buf, convective_velocity_scale_flag, .true.)
    indPrecip =      fu_index(met_buf, total_precipitation_rate_flag, .true.)
    if(srcPollen%ems2wholeABL) indBLH = fu_index(met_buf, abl_height_m_flag, .true.)
    
    ! pollen (amount left, annual correction, ready to fly)
    indPollenLeft   =   fu_get_buffer_index(srcPollen, disp_buf, pollen_left_relative_flag, .true.)
    indPollenTotal  =   fu_get_buffer_index(srcPollen, disp_buf, pollen_total_per_m2_flag, .true.)
    indAnnualTotalCorr =fu_get_buffer_index(srcPollen, disp_buf, pollen_correction_flag, .true.)
    indPollenRdyToFly = fu_get_buffer_index(srcPollen, disp_buf, pollen_rdy_to_fly_flag, .true.) 
    
    ! Allergen (pollen potency, allergen ready to fly , nothing for free allergen this far)
    if(.not. srcPollen%indPolAlrg == int_missing)then 
      indAlrgRdyToFly = fu_get_buffer_index(srcPollen, disp_buf, allergen_rdy_to_fly_flag, .true.) 
      indPotency      = fu_get_buffer_index(srcPollen, disp_buf, pollen_potency_flag, .true.) 
    endif
    
    ! Fields depending on the emission computation principle
    ! Calendar day
    if(srcPollen%ifStartCDThr .or. srcPollen%ifEndCDThr .or. srcPollen%ifStartEndGammaThr)then 
      ! day and start/end thresholds
      iDayInYear = fu_julian_date(now) - nint(fu_days(srcPollen%timeMapShift))
      if(srcPollen%ifStartCDThr .or. srcPollen%ifStartEndGammaThr)then 
        indStartCDThr = fu_get_buffer_index(srcPollen, disp_buf, start_calday_threshold_flag, .true.) 
      endif
      if(srcPollen%ifEndCDThr .or. srcPollen%ifStartEndGammaThr)then
        indEndCDThr = fu_get_buffer_index(srcPollen, disp_buf, end_calday_threshold_flag, .true.) 
      endif
      if(srcPollen%ifStartCDThr .and. srcPollen%ifEndCDThr)then 
        !season length
        indCDDiff = fu_get_buffer_index(srcPollen, disp_buf, calday_start_end_diff_flag, .true.) 
      endif
    endif  
    
    ! Chillsum
    if(srcPollen%ifChillSum)then
      indCS = fu_get_buffer_index(srcPollen, disp_buf, chillsum_flag, .true.) 
      if(srcPollen%heatsum_params%iChillSumType == csSarvasGeneralised)then
        indT2m = fu_index(met_buf, temperature_2m_flag, .true.)
      endif
    endif  ! chillsum
    
    ! Heatsum
    if(srcPollen%ifStartHSThr .or. srcPollen%ifEndHSThr)then ! temperature, heatsum stuff, thresholds
      indT2m = fu_index(met_buf, temperature_2m_flag, .true.)
      indHS = fu_get_buffer_index(srcPollen, disp_buf, heatsum_flag, .true.) 
      indStartDay = fu_get_buffer_index(srcPollen, disp_buf, growth_season_start_day_flag, .true.) 
      if(srcPollen%heatsum_params%iHeatSumType == hsDegreeDay .or. &
       & srcPollen%heatsum_params%iHeatSumType == hsBioDay)then  ! daily mean
        indDailyTempr = fu_index(disp_buf, day_mean_temperature_2m_flag, .true.)
      endif
      if(srcPollen%heatsum_params%iHeatSumType == hsDegreeDay .or. &
       & srcPollen%heatsum_params%iHeatSumType == hsDegreeHour)then
        indTemprCutOff = fu_get_buffer_index(srcPollen, disp_buf, heatsum_cutoff_tempr_flag, .true.) 
      endif
      if(srcPollen%ifStartHSThr)then !threshold
        indStartHSThr = fu_get_buffer_index(srcPollen, disp_buf, start_heatsum_threshold_flag, .true.) 
      endif
      if(srcPollen%ifEndHSThr)then ! threshold
        indEndHSThr = fu_get_buffer_index(srcPollen, disp_buf, end_heatsum_threshold_flag, .true.) 
      endif
      if(srcPollen%ifStartHSThr .and. srcPollen%ifEndHSThr)then !season length
        indHSDiff = fu_get_buffer_index(srcPollen, disp_buf, heatsum_start_end_diff_flag, .true.) 
      endif
    endif   ! heatsum
    
    ! Low temperature to end the season
    if(srcPollen%ifTempThr)then   !T2m and threshold
      indT2m = fu_index(met_buf, temperature_2m_flag, .true.)
      indTempThr = fu_get_buffer_index(srcPollen, disp_buf, temperature_threshold_flag, .true.) 
    endif
    if(srcPollen%ifDayTempThr)then   !daily T2m and threshold
      indDailyTempr = fu_index(disp_buf, day_mean_temperature_2m_flag,.true.)
      indDayTempThr = fu_get_buffer_index(srcPollen, disp_buf, daily_temp_threshold_flag, .true.) 
    endif
    ! Soil moisture 
    if(srcPollen%ifSWthr)then ! soil moisture and threshold
      indSoilWater = fu_index(met_buf, soil_moisture_vol_frac_nwp_flag, .true.)
      indSoilWaterThr = fu_get_buffer_index(srcPollen, disp_buf, soil_moisture_threshold_flag, .true.) 
    endif
    
    ! Soil moisture for plant growth
    if(srcPollen%ifSWGrowth)then 
      indSoilWater =  fu_index(met_buf, soil_moisture_vol_frac_nwp_flag, .true.)
      indPlantGrowth =fu_get_buffer_index(srcPollen, disp_buf, plant_growth_flag, .true.)
      indStartDay =   fu_get_buffer_index(srcPollen, disp_buf, growth_season_start_day_flag, .true.) 
    endif

    if(error)return

#ifdef CAMS_DUMP
iDayInYear_ = fu_julian_date(now)
mdl_timestep_sec_ = fu_sec(timestep)
iCAMSdisp = ixCAMSdisp + (iyCAMSdisp-1)*nx_dispersion
iCAMSmet = fu_grid_index(nx_meteo, ixCAMSdisp, iyCAMSdisp, pHorizInterpMet2DispStruct)
now_sec_since_sunrise_ = fu_sec(now - fu_sunrise_utc(fu_lon_geographical_from_grid(real(ixCAMSdisp),  real(iyCAMSdisp), dispersion_grid), &
                                                   & fu_lat_geographical_from_grid(real(ixCAMSdisp),  real(iyCAMSdisp), dispersion_grid), now)) ! seconds after sunrise
dayLength_hours_ = fu_hours(fu_daylength(fu_lat_geographical_from_grid(real(ixCAMSdisp), &
                                                            & real(iyCAMSdisp), dispersion_grid), now))

T2m_ = met_buf%p2d(indT2m)%present%ptr(iCAMSmet)
DailyTempr_ = min(400., max(250., disp_buf%p2d(indDailyTempr)%present%ptr(iCAMSdisp)))
StartCDThr_ = disp_buf%p2d(indStartCDThr)%present%ptr(iCAMSdisp)
HS_in_ = disp_buf%p2d(indHS)%present%ptr(iCAMSdisp)
PollenLeft_in_ = disp_buf%p2d(indPollenLeft)%present%ptr(iCAMSdisp)
PollenRdyToFly_in_ = disp_buf%p2d(indPollenRdyToFly)%present%ptr(iCAMSdisp)
#endif    
    
    !
    ! If necessary, update chill sum
    !
    if(srcPollen%ifChillSum)then
      select case(srcPollen%heatsum_params%iChillSumType)
        case(csSarvasGeneralised)
          fTmp = fu_min(now) * 60.0 + fu_sec(now)
          if(fTmp < timestep_sec .or. (fTmp + timestep_sec > 3600))then ! less than timestep to round hour
            call update_chill_sum(met_buf, disp_buf, indCS, indT2m, &
                                & one_hour, now, timestep, srcPollen%heatsum_params)
          endif
        case default
          call set_error('Unknown chill sum type: >>' + &
                & fu_str(srcPollen%heatsum_params%iChillSumType) + '<<','compute_emission_for pollen')
          return
      end select
    endif
    !
    ! If necessary, update the heatsum, possibly, accounting for the chill sum
    !
    if(srcPollen%ifStartHSThr .or. srcPollen%ifEndHSThr)then
      select case(srcPollen%heatsum_params%iHeatSumType)
        case(hsDegreeDay)
          if(fu_abs(now - fu_start_of_day_utc(now)) < timestep)then ! less than timestep to midnight
            call update_heat_sum(met_buf, disp_buf, indHS, indStartDay, indDailyTempr, &
                               & one_day, now, timestep, srcPollen%heatsum_params)
          endif

        case(hsDegreeHour)
          fTmp = fu_min(now) * 60.0 + fu_sec(now)
          if(fTmp < timestep_sec .or. (fTmp + timestep_sec > 3600))then ! less than timestep to round hour
            call update_heat_sum(met_buf, disp_buf, indHS, indStartDay, indT2m, &
                               & one_hour, now, timestep, srcPollen%heatsum_params)
          endif
      
        case(hsBioDay)
          call update_heat_sum(met_buf, disp_buf, indHS, indStartDay,indT2m, timestep, now, timestep, &
                             & srcPollen%heatsum_params)

        case(hsSigmoidPeriodUnits)
          fTmp = fu_min(now) * 60.0 + fu_sec(now)
          if(fTmp < timestep_sec .or. (fTmp + timestep_sec > 3600))then ! less than timestep to round hour
            call update_heat_sum(met_buf, disp_buf, &
                               & indHS, indStartDay, indT2m, &
                               & one_hour, now, timestep, srcPollen%heatsum_params)
          endif
      end select
    endif

#ifdef CAMS_DUMP
HS_out_ = disp_buf%p2d(indHS)%present%ptr(iCAMSdisp)
fMassInjected_out_ = 0
#endif

    !
    ! Main emission cycle over the whole grid.
    !
!    call report(pHorizInterpMet2DispStruct%gridfrom)
!    call report(pHorizInterpMet2DispStruct%gridto)

    do iy = 1, ny_dispersion
      do ix = 1, nx_dispersion
!if(ix == 50 .and. iy == 30)then
!  call msg('Start (50,40) emission')
!endif
        iDisp = ix + (iy-1) * nx_dispersion  
        if(pSrcMask(iDisp) < 1e-10)then
!  Due to imprecise grids interpolation PollenLeft diffused and cut off on every MPI restart
! here and in few places above. Made it nearest-neighbour to avoid diffusion
!              if (disp_buf%p2d(indPollenLeft)%future%ptr(iDisp) /= 0) then 
!                      call msg("Cutting pollen left pSrcMask"//trim(srcPollen%src_nm)//", cell, value, maskval", &
!                       & (/real(iDisp), disp_buf%p2d(indPollenLeft)%future%ptr(iDisp), pSrcMask(iDisp)/) )
                      disp_buf%p2d(indPollenLeft)%future%ptr(iDisp) = 0.0
!              endif
          cycle
        endif
        iMeteo =  fu_grid_index(nx_meteo, ix, iy, pHorizInterpMet2DispStruct)
        if(fLandFractionMet(iMeteo) < 0.01)cycle  ! speedup but cannot zero pollen_left
        
        ! Check against flowering thresholds
        ! First the meteo thresholds, that can zero the heatsum!!!
        ifEms = .true.
        if(srcPollen%ifTempThr)then
          if(met_buf%p2d(indT2m)%present%ptr(iMeteo) < disp_buf%p2d(indTempThr)%present%ptr(iDisp))then
            ifEms = .false.
            if(srcPollen%ifStartHSThr .or. srcPollen%ifHSGrowth)then ! Zero the heatsum 
              disp_buf%p2d(indHS)%future%ptr(iDisp) = 0.0               
            endif     
            ! If flowering started, zero the pollen left to end the emission in gridcell
            if(disp_buf%p2d(indPollenLeft)%present%ptr(iDisp) > 1.0e-3 .and. &
             & 1.0 - disp_buf%p2d(indPollenLeft)%present%ptr(iDisp) > 1.0e-3) then 
!              if (disp_buf%p2d(indPollenLeft)%future%ptr(iDisp) /= 0) then 
!                      call msg("Cutting pollen left ifTempThr"//trim(srcPollen%src_nm)//", cell, value", &
!                       & real(iDisp), disp_buf%p2d(indPollenLeft)%future%ptr(iDisp))
                      disp_buf%p2d(indPollenLeft)%future%ptr(iDisp) = 0.0
!              endif
            endif
          endif
        endif
        if(srcPollen%ifDayTempThr)then
          if(disp_buf%p2d(indDailyTempr)%present%ptr(iDisp) < disp_buf%p2d(indDayTempThr)%present%ptr(iDisp))then
            ifEms = .false.
            if(srcPollen%ifStartHSThr .or. srcPollen%ifHSGrowth)then ! Zero the heatsum 
              disp_buf%p2d(indHS)%future%ptr(iDisp) = 0.0               
            endif     
            ! If flowering started, zero the pollen left to end the emission in gridcell
            if(disp_buf%p2d(indPollenLeft)%present%ptr(iDisp) > 1.0e-3 .and. &
             & 1.0 - disp_buf%p2d(indPollenLeft)%present%ptr(iDisp) > 1.0e-3) then
!              if (disp_buf%p2d(indPollenLeft)%future%ptr(iDisp) /= 0) then 
!                      call msg("Cutting pollen left ifDayTempThr"//trim(srcPollen%src_nm)//", cell, value", &
!                       & real(iDisp), disp_buf%p2d(indPollenLeft)%future%ptr(iDisp))
                      disp_buf%p2d(indPollenLeft)%future%ptr(iDisp) = 0.0
!              endif
            endif
          endif
        endif
        if(srcPollen%ifSWThr)then
          if(met_buf%p2d(indSoilWater)%present%ptr(iMeteo) < disp_buf%p2d(indSoilWaterThr)%present%ptr(iDisp))then
            ifEms = .false.  
            ! If flowering started, zero the pollen left to end the emission in gridcell
            if(disp_buf%p2d(indPollenLeft)%present%ptr(iDisp) > 1.0e-3 .and. &
             & 1.0 - disp_buf%p2d(indPollenLeft)%present%ptr(iDisp) > 1.0e-3) then 
!              if (disp_buf%p2d(indPollenLeft)%future%ptr(iDisp) /= 0) then 
!                      call msg("Cutting pollen left ifSWThr"//trim(srcPollen%src_nm)//", cell, value", &
!                       & real(iDisp), disp_buf%p2d(indPollenLeft)%future%ptr(iDisp))
                      disp_buf%p2d(indPollenLeft)%future%ptr(iDisp) = 0.0
!             endif
            endif
          endif
        endif

#ifdef CAMS_DUMP
if(ix == ixCAMSdisp .and. iy == iyCAMSdisp)then
 HS_out_ = disp_buf%p2d(indHS)%present%ptr(iCAMSdisp)
endif
#endif        

        ifSeasonStarted = .true.  ! formal, without uncertainty
        if(srcPollen%ifStartHSThr)then
          ifSeasonStarted = disp_buf%p2d(indHS)%present%ptr(iDisp) > &
                          & disp_buf%p2d(indStartHSThr)%present%ptr(iDisp)
          if(disp_buf%p2d(indHS)%present%ptr(iDisp) < &
                                      & (1. - srcPollen%fUncertainty_HS_relative_start) * &
                                      & disp_buf%p2d(indStartHSThr)%present%ptr(iDisp)) then
            ifEms = .false.  ! outside uncertainty too, so no emission
            ! delaying past calendar day threshold - change its value to preserve distribution function
            if(srcPollen%ifStartCDThr)then
              if(iDayInYear > disp_buf%p2d(indStartCDThr)%present%ptr(iDisp))then
                disp_buf%p2d(indStartCDThr)%future%ptr(iDisp) = fu_julian_date_real(now)
                if(srcPollen%ifEndCDThr)disp_buf%p2d(indCDDiff)%future%ptr(iDisp) = &
                                      & max(0.0, disp_buf%p2d(indEndCDThr)%future%ptr(iDisp) - &
                                               & disp_buf%p2d(indStartCDThr)%future%ptr(iDisp))
              endif
            endif  ! ifStart CD threshold
          endif  ! if not yet enough heat sum
        endif  ! if start HS threshold
        if(srcPollen%ifStartCDThr)then
          ifSeasonStarted = ifSeasonStarted .and. &
                          & (iDayInYear > disp_buf%p2d(indStartCDThr)%present%ptr(iDisp))
          if(iDayInYear < disp_buf%p2d(indStartCDThr)%present%ptr(iDisp) - &
                        & srcPollen%fUncertainty_CD_days_start) &
           & ifEms = .false.
        endif
        if(srcPollen%ifStartEndGammaThr)then
          !
          ! So far, do nothing. In principle, starting tail can be checked
          !
        endif
!if(ix == 20 .and. iy == 25)then
!  if (ifSeasonStarted .and. ifEms)then
!    call msg('Season check 198-121 season started and emission,ix,iy:',ix,iy)
!  elseif (ifSeasonStarted) then
!    call msg('Season check 198-121 season started and but not emission')
!  elseif (ifEms)then
!    call msg('Season check 198-121 season NOT started but emission true')
!  else
!    call msg('Season check 198-121 season not started and emission is false')
!  endif
!endif
! No need for end checking here: the open-pocket principle serves us here: 
! uncertainty in the total pollen amount translates into the uncertainty of the season end
!
!        if(srcPollen%ifEndHSThr)then
!          if(disp_buf%p2d(indHS)%present%ptr(iDisp) >= &
!                                & (1 + srcPollen%fUncertainty_HS_relative) * &
!                                 & disp_buf%p2d(indEndHSThr)%present%ptr(iDisp)) &
!           & ifEms = .false.
!        endif
!        if(srcPollen%ifEndCDThr)then
!          if(iDayInYear >= &
!           & disp_buf%p2d(indEndCDThr)%present%ptr(iDisp) + srcPollen%fUncertainty_CD_days) &
!           & ifEms = .false.
!        endif
        
        ! If there's pollen waiting in buffer and allergen is requested, 
        ! potency growth in buffer has to happen before new pollen comes
        fPolAlrg = 0.0
        if(srcPollen%indPolAlrg /= int_missing)then
          if(disp_buf%p2d(indPollenRdyToFly)%future%ptr(iDisp) > 0.0)then
            fPolAlrg = fu_make_allergen_in_rdy_pollen(&
                             & disp_buf%p2d(indPollenRdyToFly)%future%ptr(iDisp), &
                             & srcPollen%alrgGrowthRate)
          endif
        endif

        if((srcPollen%iRipeningType == prHSLinear .or. srcPollen%iRipeningType == prHSNormal))then
          !
          ! If ripening is linear with regard to heatsum derivative, have to compute that one
          !
          select case(srcPollen%heatsum_params%iHeatSumType)
          case(hsDegreeDay)
            dHSdt = max(0.0,(met_buf%p2d(indT2m)%present%ptr(iMeteo) - &
                                          & disp_buf%p2d(indTemprCutOff)%present%ptr(iDisp))) / 24. ! scale to hour
          case(hsDegreeHour)
            dHSdt = max(0.0,(met_buf%p2d(indT2m)%present%ptr(iMeteo) - &
                                          & disp_buf%p2d(indTemprCutOff)%present%ptr(iDisp)))
          case(hsSigmoidPeriodUnits)
            if(met_buf%p2d(indT2m)%present%ptr(iMeteo) >= &
                                                     & srcPollen%heatsum_params%TemprCutOff)then
              dHSdt = fu_sigmoid_hsum_response_tempr(met_buf%p2d(indT2m)%present%ptr(iMeteo), &
                                                   & srcPollen%heatsum_params)
            else
              dHSdt = 0.0
            endif
          case(hsBioDay)
            dHSdt = fu_btime_tempr_response(met_buf%p2d(indT2m)%present%ptr(iMeteo)) *  &
                  & fu_btime_daylen_response(fu_hours(fu_daylength(fu_lat_geographical_from_grid( &
                                                                              & real(ix), real(iy), &
                                                                              & dispersion_grid), &
                                                                & now)), &
                                           & disp_buf%p2d(indHS)%present%ptr(iDisp))
          end select
          dHSdt = dHSdt * fu_hours(timestep)
        endif 
        
        ! Total amount of pollen determined by the pre-season conditions
        if((srcPollen%ifSWGrowth .or. srcPollen%ifHSGrowth) .and. (.not. ifSeasonStarted))then 
          disp_buf%p2d(indPollenTotal)%future%ptr(iDisp) = &
                              & disp_buf%p2d(indAnnualTotalCorr)%present%ptr(iDisp) * srcPollen%standardPollenTotal
          if(srcPollen%ifSWGrowth)then  
            ! Plant response to water in soil
            if(met_buf%p2d(indSoilWater)%present%ptr(iMeteo) < srcPollen%SWGrowthDeath)then
              disp_buf%p2d(indPlantGrowth)%future%ptr(iDisp) = 0.0   
            else 
              fTmp = fu_plant_growth_SW(met_buf%p2d(indSoilWater)%present%ptr(iMeteo))
              ! Normalization depending on whether the main driver is time or thermal time
              if(srcPollen%iRipeningType == prHSLinear .or. srcPollen%iRipeningType == prHSNormal )then
                fTmp = fTmp * dHSdT / disp_buf%p2d(indStartHSThr)%present%ptr(iDisp)
              elseif(srcPollen%iRipeningType == prCDLinear .or. srcPollen%iRipeningType == prCDNormal &
                   & .or. srcPollen%iRipeningType == prCDNormDrn)then
                fTmp = fTmp * fu_sec(timestep) / 86400.0 / &
                & (disp_buf%p2d(indStartCDThr)%present%ptr(iDisp)-disp_buf%p2d(indStartDay)%present%ptr(iDisp))
              else 
                call set_error('No idea how to normalize plant growth', 'compute_emission_for_pollen')
              endif
              if(fTmp < 0.0)then
                call set_error('Negative growth','compute_emission_for_pollen')
                return
              endif
              disp_buf%p2d(indPlantGrowth)%future%ptr(iDisp) = disp_buf%p2d(indPlantGrowth)%present%ptr(iDisp) + fTmp
            endif
            disp_buf%p2d(indPollenTotal)%future%ptr(iDisp) = disp_buf%p2d(indPollenTotal)%future%ptr(iDisp) * &
                                                           & disp_buf%p2d(indPlantGrowth)%future%ptr(iDisp)
          endif  ! if pre-season growth of the plant
          if(srcPollen%ifHSGrowth)then  ! Linear to accumulated heatsum
              disp_buf%p2d(indPollenTotal)%future%ptr(iDisp) = disp_buf%p2d(indPollenTotal)%future%ptr(iDisp) * &
                         & disp_buf%p2d(indHS)%present%ptr(iDisp) / disp_buf%p2d(indStartHSThr)%present%ptr(iDisp)
          endif
        endif

        ! Buffered model. Two steps are to be made:
        ! 1. Produce more ready-to-fly-pollen
        ! 2. In case of favorable other parameters, release it.
        !
        if(disp_buf%p2d(indPollenTotal)%future%ptr(iDisp) < 1.0e-6)then
          disp_buf%p2d(indPollenLeft)%future%ptr(iDisp) = 0.0
          cycle
        endif  ! No pollen in this gridcell

        if(ifEms .and. disp_buf%p2d(indPollenLeft)%present%ptr(iDisp) > 1.0e-6)then ! Produce new pollen (pollen ready to fly will fly anyway)
          
          fTmp = 0.0
          select case(srcPollen%iRipeningType)
          case(prCDLinear)
             fTmp = fu_mk_new_pollen_linear_dual(real(fu_julian_date_real(now) - &
                                                            & fu_days(srcPollen%timeMapShift), &
                                                    & DEFAULT_REAL_KIND), &
                                          & disp_buf%p2d(indStartCDThr)%present%ptr(iDisp), &
                                          & disp_buf%p2d(indEndCDThr)%present%ptr(iDisp), &
                                          & fu_sec(timestep) / 86400., &
                                          & srcPollen%fUncertainty_CD_days_start / &
                                                 & disp_buf%p2d(indStartCDThr)%present%ptr(iDisp), &
!                                                    & disp_buf%p2d(indCDDiff)%present%ptr(iDisp), &
                                           & srcPollen%fUncertainty_tot_pollen_relat, &
                                           & disp_buf%p2d(indPollenTotal)%future%ptr(iDisp), &
                                           & srcPollen%UncertaintyParams)
          case(prHSLinear)
            if(dHSdt > 0.0) &
               & fTmp = fu_mk_new_pollen_linear(disp_buf%p2d(indHS)%present%ptr(iDisp), &
                                        & disp_buf%p2d(indStartHSThr)%present%ptr(iDisp), &
                                        & disp_buf%p2d(indEndHSThr)%present%ptr(iDisp), &
                                        & dHSdt,  &
                                        & srcPollen%fUncertainty_HS_relative_start, &
                                        & srcPollen%fUncertainty_tot_pollen_relat, &
                                        & disp_buf%p2d(indPollenTotal)%future%ptr(iDisp), &
                                        & srcPollen%UncertaintyParams)

          case(prCDNormal) ! start & end by 2.5% criterion would give you 2 sigma from middle
             fTmp = fu_mk_new_pollen_normal(real(fu_julian_date_real(now) - &
                                                          & fu_days(srcPollen%timeMapShift), &
                                               & DEFAULT_REAL_KIND), &
                                    & 0.5 * (disp_buf%p2d(indStartCDThr)%present%ptr(iDisp) + &
                                           & disp_buf%p2d(indEndCDThr)%present%ptr(iDisp)), &
                                    & 0.25 * disp_buf%p2d(indCDDiff)%present%ptr(iDisp), & !srcPollen%fUncertainty_CD_days, &
                                    & fu_sec(timestep)/86400., &
                                    & disp_buf%p2d(indPollenTotal)%future%ptr(iDisp))
                                    
          case(prHSNormal) ! start & end by 2.5% criterion would give you 2 sigma from middle
            if(dHSdt > 0.0) &
               & fTmp = fu_mk_new_pollen_normal(disp_buf%p2d(indHS)%present%ptr(iDisp), &
                                    & 0.5 * (disp_buf%p2d(indStartHSThr)%present%ptr(iDisp) + &
                                          & disp_buf%p2d(indEndHSThr)%present%ptr(iDisp)), &
                                    & 0.25 * disp_buf%p2d(indHSDiff)%present%ptr(iDisp), &   !srcPollen%fUncertainty_HS_relative, &
                                    & dHSdt, &
                                    & disp_buf%p2d(indPollenTotal)%future%ptr(iDisp))                          

          case(prCDNormDrn) ! start & end by 2.5% criterion would give you 2 sigma from middle
            fTmp = fu_mk_new_pollen_normal_diurnal(now, fu_sec(timestep), &
                                          & fu_lon_geographical_from_grid(real(ix),  real(iy), dispersion_grid), &
                                          & fu_lat_geographical_from_grid(real(ix),  real(iy), dispersion_grid), &
                                          & 0.5 * (disp_buf%p2d(indStartCDThr)%present%ptr(iDisp) + &
                                                               & disp_buf%p2d(indEndCDThr)%present%ptr(iDisp)), &
                                          & 0.25 * disp_buf%p2d(indCDDiff)%present%ptr(iDisp), & !srcPollen%fUncertainty_CD_days, &
                                          & srcPollen%dayMid1, srcPollen%dayMid2, &
                                          & srcPollen%daySgm1, srcPollen%daySgm2, &
                                          & srcPollen%dayFrac1, &
                                          & disp_buf%p2d(indPollenTotal)%future%ptr(iDisp))
            
          case(prCDGammaWTails) ! season described by Gamma-type distribution with elevated tails
                              ! It requires relative terms for time axis: starting_%=0, ending_%=1
!            fTmp = fu_mk_new_pollen_gamma_w_tails( &
            fTmp = fu_mk_new_pollen_gamma_tails_2( &
                            & (real(fu_julian_date_real(now) - &    ! distance from start_% ...
                                         & fu_days(srcPollen%timeMapShift), &
                                  & DEFAULT_REAL_KIND) - &
                                            & disp_buf%p2d(indStartCDThr)%present%ptr(iDisp)) / &
                            & (disp_buf%p2d(indEndCDThr)%present%ptr(iDisp) - &  ! ...relative to season length
                                            & disp_buf%p2d(indStartCDThr)%present%ptr(iDisp)), &
                            & disp_buf%p2d(indPollenTotal)%future%ptr(iDisp), &             !fTotalPollen, 
                            & fu_sec(timestep) / 86400.0 / (disp_buf%p2d(indEndCDThr)%present%ptr(iDisp) - &  ! ...relative to season length
                                            & disp_buf%p2d(indStartCDThr)%present%ptr(iDisp)), &   !delta_t, &
                            & srcPollen%gf_beta, &       ! coef at exponent in gamma distr.function
                            & srcPollen%gf_start_times, &  ! start moments for power terms
                            & srcPollen%gf_scales, &     ! scales for power terms
                            & srcPollen%gf_powers, &     ! powers of the power terms
                            & srcPollen%gf_nTerms)       ! number of power terms
          end select 

          if(fTmp < 0.0)then
            call set_error('Negative emission','compute_emission_for_pollen')
            call msg('fTmp', fTmp)
            call msg('iy, ix', iy, ix)
            return
          endif
          
          fTmp = min(fTmp, disp_buf%p2d(indPollenLeft)%future%ptr(iDisp) *  &
                                                            & disp_buf%p2d(indPollenTotal)%future%ptr(iDisp))
          disp_buf%p2d(indPollenLeft)%future%ptr(iDisp) =  disp_buf%p2d(indPollenLeft)%future%ptr(iDisp) - &
                                                      & fTmp / disp_buf%p2d(indPollenTotal)%future%ptr(iDisp)
          disp_buf%p2d(indPollenRdyToFly)%future%ptr(iDisp) = &
                                         & disp_buf%p2d(indPollenRdyToFly)%future%ptr(iDisp) + fTmp
          
          ! If necessary, allergen in new pollen needs to be added to allergen production
          if((.not. (srcPollen%indPolAlrg == int_missing)) .and. (fTmp > 0.0))then
            fPolAlrg = fPolAlrg + fu_make_allergen_in_new_pollen(fTmp,  &
                                                      & disp_buf%p2d(indPotency)%future%ptr(iDisp))

          endif
!if(ix == 20 .and. iy == 25)then
!  call msg('Season check 50-40 end of pollen production. Produced:', fTmp)
!endif
        endif ! if make new pollen
        leftPollen = leftPollen + disp_buf%p2d(indPollenLeft)%future%ptr(iDisp) * &
                              & disp_buf%p2d(indPollenTotal)%future%ptr(iDisp) * &
                              & xSize(iDisp) * ySize(iDisp) * pSrcMask(iDisp) ! for reporting
        
        ! If necessary, new allergen needs to be added to buffer
        if((.not. (srcPollen%indPolAlrg == int_missing)) .and. (fPolAlrg > 0.0))then
            disp_buf%p2d(indAlrgRdyToFly)%future%ptr(iDisp) = &
                                         & disp_buf%p2d(indAlrgRdyToFly)%future%ptr(iDisp) + fPolAlrg
        endif
        
        !
        ! Release whatever needs to be released. 
        if(disp_buf%p2d(indPollenRdyToFly)%future%ptr(iDisp) < 1.0e-6)cycle ! speedup

        if(srcPollen%iReleaseType == relInst)then 
          fTmp = disp_buf%p2d(indPollenRdyToFly)%future%ptr(iDisp)
          disp_buf%p2d(indPollenRdyToFly)%future%ptr(iDisp) = 0.0
        elseif(met_buf%p2d(indPrecip)%present%ptr(iMeteo) < srcPollen%precipThreshold) then ! Precip
                  
          ! Now check humidity threshold (expensive, do it only if rest is OK)
         !          call msg("met_buf%buffer_quantities", met_buf%buffer_quantities )
          !fTmp = fu_get_value(met_buf, indHumid, ixMeteo, iyMeteo, 1, nx_meteo, now)
          !Humidity from the first meteo layer
          fTmp = met_buf%p2d(indHumid)%present%ptr(iMeteo)

          if(error)return
          if(fTmp < srcPollen%HighHumidThresh)then
            !
            ! All thresholds are OK, can make the release of the ready-made pollen
            ! An uncertainty here is that we have little idea on "reference" rate of
            ! release in case of some "reference" conditions.
            ! So, we take the linear discharge model: dm/dt=-alpha(meteo)*m
            ! and roughly estimate alpha from a statement: in case of very good conditions,
            ! all ready-to-fly pollen must be released within one time step.
            fTmp = fu_pollen_emis_from_buffer(fTmp, &                                             ! humidity
                                          & met_buf%p2d(indWind10m)%present%ptr(iMeteo), &      ! 10m wind speed
                                          & met_buf%p2d(indConvVelocity)%present%ptr(iMeteo), & ! conv velocity scale
                                          & met_buf%p2d(indPrecip)%present%ptr(iMeteo), &       ! precipitation
                                          & srcPollen%LowHumidThresh, srcPollen%HighHumidThresh, & ! humidity thresholds
                                          & srcPollen%precipThreshold, &                           ! precipitation threshold
                                          & srcPollen%windSaturation, srcPollen%windMaxScale, &    ! windspeed saturation
                                          & srcPollen%fBufReleaseLimit, &                          ! max rate of release
                                          & disp_buf%p2d(indPollenRdyToFly)%future%ptr(iDisp), &  ! just-updated amt of ready pollen
                                          & fu_sec(timestep), &
                                          & srcPollen%iReleaseType)
                    
            ! If we emit some pollen, less is ready to fly for the future
            if(fTmp < 1.0e-6)then                                                    ! skip very low release
              fTmp = 0.
            elseif(fTmp > disp_buf%p2d(indPollenRdyToFly)%future%ptr(iDisp))then  ! cut too high release
              fTmp = disp_buf%p2d(indPollenRdyToFly)%future%ptr(iDisp)
              disp_buf%p2d(indPollenRdyToFly)%future%ptr(iDisp) = 0.0
            else                                                                  ! general case
              disp_buf%p2d(indPollenRdyToFly)%future%ptr(iDisp) = &
                                         & disp_buf%p2d(indPollenRdyToFly)%future%ptr(iDisp) - fTmp
            endif  ! values of fTmp
          else
            fTmp = 0.
          endif  ! Humidity is low enough for the release
        else
            fTmp = 0.
        endif  ! no precipitation
        rdyPollen = rdyPollen + disp_buf%p2d(indPollenRdyToFly)%future%ptr(iDisp) * &
                                    & xSize(iDisp) * ySize(iDisp) * pSrcMask(iDisp) ! for reporting                    
        if(fTmp < 1.0e-6)cycle 

!if(ix == 20 .and. iy == 25)then
!  call msg('Season check 50-40 final release prior to scaling with cell size:',fTmp)
!endif

         ! Pollen allergen emission 
         ! has to be called before scaling the pollen emission to grid size (same fraction emitted)!!
        if(.not. srcPollen%indPolAlrg == int_missing)then
          fPolAlrg = fu_allergen_emis_from_buffer(disp_buf%p2d(indAlrgRdyToFly)%future%ptr(iDisp), &
                                                & fTmp, disp_buf%p2d(indPollenRdyToFly)%future%ptr(iDisp))
          disp_buf%p2d(indAlrgRdyToFly)%future%ptr(iDisp) = &
                                         & disp_buf%p2d(indAlrgRdyToFly)%future%ptr(iDisp) - fPolAlrg
          rdyAllergen = rdyAllergen + disp_buf%p2d(indAlrgRdyToFly)%future%ptr(iDisp) *  &
                                    & xSize(iDisp) * ySize(iDisp) * pSrcMask(iDisp) ! for reporting !!!incorrect because cycle if no ems !!!
          fPolAlrg = fPolAlrg * xSize(iDisp) * ySize(iDisp) * pSrcMask(iDisp)
        endif
        
#ifdef CAMS_DUMP
if(ix == ixCAMSdisp .and. iy == iyCAMSdisp)then
 fMassInjected_out_ = fTmp
endif
#endif        
        
        ! Finally, scale from m2 to grid cell size and bring the emitted amount of pollen into 
        ! the emission mass map
        fTmp = fTmp * xSize(iDisp) * ySize(iDisp) * pSrcMask(iDisp)

        ! Fill-in the emission map, not forgetting the number to mass convertion where needed
        check = 0.0
        do iLev = 1, srcPollen%nLevsDispVert

          if(srcPollen%ems2wholeABL)then
            if(srcPollen%fEmissionBottom > disp_layer_top_m(iLev) .or. &
              & met_buf%p2d(indBLH)%present%ptr(iMeteo) < disp_layer_top_m(iLev-1))cycle
         
            overlapBottom = max(srcPollen%fEmissionBottom, disp_layer_top_m(iLev-1))
            overlapTop = min(met_buf%p2d(indBLH)%present%ptr(iMeteo), disp_layer_top_m(iLev))
            fzDisp = (overlapBottom + overLapTop) / 2.0

            fzDisp =  (fzDisp - (disp_layer_top_m(iLev) + disp_layer_top_m(iLev-1))/2) / &
                   &  (disp_layer_top_m(iLev) - disp_layer_top_m(iLev-1))
            
            
            levFractDispVert = (overlapTop - overlapBottom) / &
                             & (met_buf%p2d(indBLH)%present%ptr(iMeteo) - srcPollen%fEmissionBottom)
            if(abs(fzDisp) > 0.5)then
              call set_error('emission to wrong layer', 'compute_emission_for_pollen')
              call msg('iLev, fzDisp', iLev, fzDisp)
              call msg('dispersion_layer_boundaries, bottom of emissoin and ABL top: ', &
                     & (/disp_layer_top_m(iLev-1), disp_layer_top_m(iLev), &
                     & srcPollen%fEmissionBottom, met_buf%p2d(indBLH)%present%ptr(iMeteo)/))
            endif
          else
            fzDisp = srcPollen%fzDisp(iLev)
            levFractDispVert = srcPollen%levFractDispVert(iLev)
          endif
          !First do the check for the overlap: speed-up
          if(levFractDispVert < 1.0e-6)cycle  ! nothing for this dispersion layer
          check = check + levFractDispVert
          fCellTotal = 0.0
          do iSp = 1, srcPollen%nSpecies
            ! Allergen emission
            if(iSp == srcPollen%indPolAlrg)then
              emsAmntLev = fPolAlrg * levFractDispVert
            elseif (iSp == srcPollen%indFreeAlrg)then
              emsAmntLev = fu_free_allergen_ems(fTmp, srcPollen%freeAlrgFract) * &
                         & levFractDispVert
            else ! Everything else emitted as pollen, to allow for instance passive tracer 
              emsAmntLev = fTmp * levFractDispVert
            endif
            !
            ! Here, if the source index mapping given, it should be used instead of the source ID number
            !
            if(ifUseSourceIdMapping)then
              iSrcNbr = arSourceIdMapping(ix,iy)
            else
              iSrcNbr = srcPollen%id_nbr
            endif
            emisMap%arM(srcPollen%adaptor%iSp(iSp),iSrcNbr,iLev,ix,iy) = &
                        & emisMap%arM(srcPollen%adaptor%iSp(iSp),iSrcNbr,iLev,ix,iy) + emsAmntLev
            emisMap%ifColumnValid(iSrcNbr,ix,iy) = .true.
            emisMap%ifGridValid(iLev,iSrcNbr) = .true.
 
            fCellTotal = fCellTotal + emsAmntLev
            emsAmnt(iSp) = emsAmnt(iSp) + emsAmntLev
            if (ifSpeciesMoment) then
              mapCoordZ%arm(srcPollen%adaptor%iSp(iSp),iSrcNbr, ilev, ix,iy) = &
                             & mapCoordZ%arm(srcPollen%adaptor%iSp(iSp),iSrcNbr, ilev, ix,iy) + &
                             & emsAmntLev * fzDisp
            end if
          end do  ! species

          if (.not. ifSpeciesMoment) then
            mapCoordZ%arM(1,iSrcNbr, iLev, ix, iy) = fzDisp * fCellTotal + &
                                                     & mapCoordZ%arM(1,iSrcNbr, iLev, ix, iy)
          end if
        end do  ! iLev
        if(abs(check-1.0)>0.001)then
            call set_error('level fractions dont add up', 'compute_emission_for_pollen')
            call msg('blh, check', met_buf%p2d(indBLH)%present%ptr(iMeteo), check)
            return
        endif
      end do  ! ix
    end do  ! iy
    !
    ! Have to update the valid time of the amount of pollen left for the future emission
    ! Do not forget that it will be valid for the next time step.
    !
    call set_valid_time(disp_buf%p2d(indPollenLeft)%future%idPtr, now+timestep)
    call set_valid_time(disp_buf%p2d(indPollenTotal)%future%idPtr, now+timestep)
    call set_valid_time(disp_buf%p2d(indPollenRdyToFly)%future%idPtr, now+timestep)
    if(.not. srcPollen%indPolAlrg == int_missing) &
      & call set_valid_time(disp_buf%p2d(indAlrgRdyToFly)%future%idPtr, now+timestep)
    if(srcPollen%ifSWGrowth) &
      & call set_valid_time(disp_buf%p2d(indPlantGrowth)%future%idPtr, now+timestep)

    ! Report the emitted and ready to fly amounts
!    call msg('Emission of pollen source:'+srcPollen%src_nm+';sector:'+srcPollen%sector_nm)
    do iSp = 1, srcPollen%nSpecies
!      call msg(fu_name(fu_material(srcPollen%species(iSp)))+':', emsAmnt(iSp))
      fMassInjected(srcPollen%adaptor%iSp(iSp)) = fMassInjected(srcPollen%adaptor%iSp(iSp)) + &
                                                & emsAmnt(iSp)
    enddo
    call msg('Pollen ready to fly:', rdyPollen)
    if(iSp == srcPollen%indPolAlrg)then
      call msg('Allergen ready to fly:', rdyAllergen)
    endif
    call msg('Pollen left: ', leftPollen)

#ifdef CAMS_DUMP
PollenLeft_out_ = disp_buf%p2d(indPollenLeft)%present%ptr(iCAMSdisp)
PollenRdyToFly_out_ = disp_buf%p2d(indPollenRdyToFly)%present%ptr(iCAMSdisp)

write(55,'(6I4,12(F,1x),D)') &
             & fu_year(now), fu_mon(now), fu_day(now), fu_hour(now), fu_min(now), iDayInYear_, &   ! Julian day of the year. IN
             & now_sec_since_sunrise_, &      ! time moment within the day, sec since sunrise, IN
             & mdl_timestep_sec_, &           ! model time step, seconds, IN
             & dayLength_hours_, &            ! day length, hours,      IN
             & T2m_, &            ! tempr at 2m and its low threshold, IN
             & DailyTempr_, &  ! daily-mean temperature and its low threshold. IN
             & StartCDThr_, &    ! calendar day of year thresholds for season start and end, IN
             & HS_in_, HS_out_, &                 ! heat sum         IN-OUT
             & PollenLeft_in_, PollenLeft_out_, &         ! fraction of pollen remained for future, IN-OUT
             & PollenRdyToFly_in_, PollenRdyToFly_out_, &     ! amnt of pollen ready to fly now, IN-OUT
             & fMassInjected_out_
#endif


  CONTAINS

    !================================================================================

    integer function fu_get_buffer_index(src, buf, quantity, ifSpeciesMandatory)
      !
      ! Just to shorten the above calls
      !
      implicit none
      
      ! Imported parameters
      type(silam_pollen_source), intent(in) :: src
      type(Tfield_buffer), pointer :: buf
      integer, intent(in) :: quantity
      logical, intent(in) :: ifSpeciesMandatory
      
      if(error)then
         call msg_warning('ENTERED with error','fu_get_buffer_index')
         return
      endif

      fu_get_buffer_index = fu_index(buf, quantity, fu_species_src(src, quantity, &
                                                                 & ifSpeciesMandatory), .true.) 
      if(error .or. fu_get_buffer_index < 1)then
        call set_error('Failed to find buffer index for:>>' + fu_quantity_string(quantity) + &
                     & ',' + fu_str(fu_species_src(src, quantity, ifSpeciesMandatory)) + '<<', & 
                     & 'fu_get_buffer_index')
        fu_get_buffer_index = int_missing
        return
      endif
      
    end function fu_get_buffer_index

    !================================================================================

    real function fu_plant_growth_SW(SW)
      !
      ! plant growth depending on whether the meteorological conditions 
      !(temperature, soil moisture) are favorable for plant growth
      !
      implicit none

      ! Imported parameters
      real, intent(in) ::SW

     fu_plant_growth_SW = 1.0 / (1.0 + exp(srcPollen%SWGrowthSigma *(SW-srcPollen%SWGrowthMid)))

    end function fu_plant_growth_SW

    !================================================================================

    real function fu_mk_new_pollen_normal(fNow, fMid, stDev, delta, fTotalPollen)
      !
      ! Returns the emission mass assuming normally distributed release
      !
      implicit none

      ! Imported parameters
      real, intent(in) :: fNow, fMid, stDev, delta, fTotalPollen
      
      ! Local variables
      real(r8k) :: f1, f2, s
      
      fu_mk_new_pollen_normal = 0.0
      if(stDev <= 0.0)return
      
      f1 = fNow - fMid
      f2 = f1 + delta
      s = stDev * sqrt(2.)
      fu_mk_new_pollen_normal = fTotalPollen * 0.5 * (ERF(f2/s)-ERF(f1/s))

    end function fu_mk_new_pollen_normal
    
    !================================================================================

    real function fu_mk_new_pollen_normal_diurnal(now, tStep, lon, lat, &
                                        & seasonMid, seasonSgm, &
                                        & dayMid1, dayMid2, daySgm1, daySgm2, dayFrac1, &
                                        & fTotalPollen)
      !
      ! Returns the emission mass assuming normally distributed season and diurnal pattern
      ! described by 2 normal modes parameterised for time after sunrise
      !
      
      implicit none

      ! Imported parameters
      type(silja_time), intent(in) :: now
      real, intent(in) :: lon, lat, tStep, &
                        & seasonMid, seasonSgm, &
                        & dayMid1, dayMid2, daySgm1, daySgm2, dayFrac1, &
                        & fTotalPollen
      real :: dayT
      real*8 :: t1, t2, s
    
      fu_mk_new_pollen_normal_diurnal = 0.0
      if(seasonSgm <= 0.0)return
      dayT = fu_sec(now - fu_sunrise_utc(lon, lat, now)) ! seconds after sunrise
 
      ! Speedup - no nighttime emission (might need to be reconsidered ..)
      if(dayT < 0. .or. dayT > 36000.0)return  
      
      ! Bimodal daily release
      t1 = dayT - dayMid1
      t2 = t1 + tStep
      s = daySgm1 * sqrt(2.)
      fu_mk_new_pollen_normal_diurnal = dayFrac1 * 0.5 * (ERF((t2)/s) - ERF(t1/s))       
    
      if(fu_mk_new_pollen_normal_diurnal < 0.0)then
        call msg('t1, t2', real(t1), real(t2))
        call msg('dayFrac1, s', dayFrac1, real(s))
        call msg('fu_mk_new_pollen_normal_diurnal', fu_mk_new_pollen_normal_diurnal)
      endif
    
      t1 = dayT - dayMid2
      t2 = t1 + tStep
      s = daySgm2 * sqrt(2.)
      fu_mk_new_pollen_normal_diurnal = fu_mk_new_pollen_normal_diurnal + &
                                      & (1.0 - dayFrac1) * 0.5 * (ERF((t2)/s) - ERF(t1/s))      
      if(fu_mk_new_pollen_normal_diurnal < 0.0)then
        call msg('t1, t2', real(t1), real(t2))
        call msg('dayFrac1, s', dayFrac1, real(s))
        call msg('fu_mk_new_pollen_normal_diurnal', fu_mk_new_pollen_normal_diurnal)
      endif
 
      ! Normal season
      t1 = real(fu_julian_date(now)) - seasonMid ! Whole day
      t2 = t1 + 1.0
      s = seasonSgm * sqrt(2.)
      fu_mk_new_pollen_normal_diurnal = fu_mk_new_pollen_normal_diurnal * fTotalPollen * &  
                              & 0.5 * (ERF(t2/s) - ERF(t1/s))       ! normal distribution
      
      if(fu_mk_new_pollen_normal_diurnal < 0.0)then
        call msg('t1, t2', real(t1), real(t2))
        call msg('s, fu_mk_new_pollen_normal_diurnal', real(s), fu_mk_new_pollen_normal_diurnal)
        call msg('fTotalPollen', fTotalPollen)
      endif
                                    

    end function fu_mk_new_pollen_normal_diurnal
    
    !================================================================================

    real function fu_mk_new_pollen_linear(fNow, fStartThres, fEndThres, delta, fURelStart, fURelEnd, &
                                        & fTotalPollen, fade_inout_params)
      !
      ! Returns the emission mass assuming linear release. Both start and end are counted from the same point,
      ! such as the heatsum start moment.
      !
      implicit none

      ! Imported parameters
      real, intent(in) :: fNow, fStartThres, fEndThres, delta, fURelStart, fURelEnd, fTotalPollen
      type(Tfade_inout), intent(in) :: fade_inout_params
     
      fu_mk_new_pollen_linear = fTotalPollen *  delta / (fEndThres - fStartThres) * & 
                    & fu_fade_in(fNow / fStartThres, fURelStart, fade_inout_params) * &   ! start treshold fade-in
                    & fu_fade_out(fNow / fEndThres, fURelEnd, fade_inout_params)      ! end treshold fade-out

if(fu_mk_new_pollen_linear < 0.)then
  call msg('Negative fu_mk_new_pollen_linear:',fu_mk_new_pollen_linear )
  call msg('fTotalPollen', fTotalPollen)
  call msg('delta, fNow:', delta, fNow)
  call msg('fEndThres fStartThres',fEndThres, fStartThres)
  call msg('fURelStart',fURelStart, fURelEnd)
  call msg('fade-in, fade_out:', fu_fade_in(fNow / fStartThres, fURelStart, fade_inout_params), &
                               & fu_fade_out(fNow / fEndThres, fURelEnd, fade_inout_params) )
endif

    end function fu_mk_new_pollen_linear
    
    !================================================================================

    real function fu_mk_new_pollen_linear_dual(fNow, fStartThres, fEndThres, delta, fURelStart, fURelEnd, &
                                             & fTotalPollen, fade_inout_params)
      !
      ! Returns the emission mass assuming linear release.
      ! Difference from the above one is that the uncertainties of the season start and end refer
      ! to different integration points: start proceeds from some early point while end proceeds 
      ! from the start point.
      !
      implicit none

      ! Imported parameters
      real, intent(in) :: fNow, fStartThres, fEndThres, delta, fURelStart, fURelEnd, &
                        & fTotalPollen
      type(Tfade_inout), intent(in) :: fade_inout_params
     
      fu_mk_new_pollen_linear_dual = fTotalPollen *  delta / (fEndThres - fStartThres) * & 
                                   & fu_fade_in(fNow / fStartThres, fURelStart, fade_inout_params) * &   ! start treshold fade-in
                                   & fu_fade_out((fNow-fStartThres) / (fEndThres-fStartThres), &
                                                      & fURelEnd, fade_inout_params) ! end fade-out

if(fu_mk_new_pollen_linear_dual < 0.)then
  call msg('Negative fu_mk_new_pollen_linear:',fu_mk_new_pollen_linear_dual )
  call msg('fTotalPollen', fTotalPollen)
  call msg('delta, fNow:', delta, fNow)
  call msg('fEndThres fStartThres',fEndThres, fStartThres)
  call msg('fURelStart',fURelStart, fURelEnd)
  call msg('fade-in, fade-out:',fu_fade_in(fNow / fStartThres, fURelStart, fade_inout_params), &
              & fu_fade_out((fNow-fStartThres) / (fEndThres-fStartThres), fURelEnd, fade_inout_params))
endif

    end function fu_mk_new_pollen_linear_dual

    !================================================================================

    real function fu_mk_new_pollen_gamma_w_tails(fNowRel, fTotalPollen, delta_t, &
                                               & beta, timesRel, scales, powers, nTerms)
      !
      ! Returns the pollen prepared for release assuming the modified "taily" Gamma distribution
      ! of the season. Tails are the reason for many parameters: have to describe the main peak
      ! via gamma-type distribution, and both elevated tails via add-on corrections.
      ! formula: rate(x)=exp(-a_1/beta)* sum(scale_i * a_i^power_i), i=1:3
      ! where a_i = max(x-timesRel_i,0)
      ! Fitting for grass showed that:
      ! beta = 0.155, 
      ! timesRel_1=0.164, scale1=13.1, power1=1.16, 
      ! timesRel_2=-0.26, scale2=12.6, power2=2.8, 
      ! timesRel_3=-0.7, scale3=0.25, power3=1.3
      !
      implicit none
      
      ! Imported parameters
      real, intent(in) :: fNowRel, fTotalPollen, delta_t, beta ! normalised time, total scaling, timestep, beta
      real, dimension(:), pointer:: timesRel, scales, powers
      integer, intent(in) :: nTerms
      
      ! Local variables
      real :: a1, a, sum, time
      integer :: iTmp, jTmp
      
      a1 = max(0.0,fNowRel-timesRel(1))
      if(a1 > beta*10.)then
        fu_mk_new_pollen_gamma_w_tails = 0.0  ! too far from the season peak (as decided by timesRel 1)
        return
      else
        ! Be careful: the rise of the function can be quite steep. First-order integration
        !
        sum = 0.0
        time = fNowRel
        do jTmp = 1,2
          do iTmp = nTerms, 1, -1
            if(time > timesRel(iTmp)) sum = sum + scales(iTmp) * (time - timesRel(iTmp))**powers(iTmp)
          end do
          sum = sum * exp(-a1 / beta)
          time = time + delta_t
        end do ! trapezoid integration cycle

        fu_mk_new_pollen_gamma_w_tails = 0.5 * fTotalPollen * delta_t * sum

      endif  ! if after the season
      
    end function fu_mk_new_pollen_gamma_w_tails
                                             
    
    !================================================================================

    real function fu_mk_new_pollen_gamma_tails_2(fNowRel, fTotalPollen, delta_t, &
                                               & beta, timesRel, scales, powers, nTerms)
      !
      ! The more rigorous implementation of the formal sum of three pieces with overall scaling.
      !
      ! Returns the pollen prepared for release assuming the modified "taily" Gamma distribution
      ! of the season. Tails are the reason for many parameters: have to describe the main peak
      ! via gamma-type distribution, and both elevated tails via add-on corrections.
      ! formula: rate(x)=exp(-a_1/beta)* sum(scale_i * a_i^power_i), i=1:3
      ! where a_i = max(x-timesRel_i,0)
      ! Fitting for grass showed that:
      ! beta = 0.155, 
      ! timesRel_1=0.164, scale1=13.1, power1=1.16, 
      ! timesRel_2=-0.26, scale2=12.6, power2=2.8, 
      ! timesRel_3=-0.7, scale3=0.25, power3=1.3
      !
      implicit none
      
      ! Imported parameters
      real, intent(in) :: fNowRel, fTotalPollen, delta_t, beta ! normalised time, total scaling, timestep, beta
      real, dimension(:), pointer:: timesRel, scales, powers
      integer, intent(in) :: nTerms
      
      ! Local variables
      real :: sum, time
      integer :: iTmp, jTmp
      
      fu_mk_new_pollen_gamma_tails_2 = 0.0  

      if(fNowRel-timesRel(1) > beta*10.) return ! too far from the season peak (as decided by timesRel 1)

      ! Be careful: the rise of the function can be quite steep. First-order integration should do however
      !
      time = fNowRel
      do jTmp = 1,2    ! trapezoid edges
        sum = 0.0
        do iTmp = 1, nTerms
          if(time > timesRel(iTmp)) sum = sum + scales(iTmp) * (time - timesRel(iTmp))**powers(iTmp)
        end do
        if(time > timesRel(1)) sum = sum * exp((timesRel(1)-time) / beta)   ! after the season peak, suppressing
        fu_mk_new_pollen_gamma_tails_2 = fu_mk_new_pollen_gamma_tails_2 + sum
        time = time + delta_t
      end do ! trapezoid integration cycle

      fu_mk_new_pollen_gamma_tails_2 = 0.5 * fTotalPollen * delta_t * sum

    end function fu_mk_new_pollen_gamma_tails_2
                                             
    
    !================================================================================
    
    real function fu_make_allergen_in_new_pollen(newPollen, potencyBasic) ! , &
                                     ! & fHumid, fWind10m, fConvVelocity, fPrecip)
      !
      ! Returns the allergen in newly made pollen
      ! to be called before new pollen is put to buffer
      ! This far just fixed potency map  
      ! Meteo dependence TO BE PARAMETERISED
      !
      implicit none

      ! Imported parameters
      real, intent(in) :: newPollen, potencyBasic ! , fHumid, fWind10m, fConvVelocity, fPrecip
     
      fu_make_allergen_in_new_pollen = newPollen * potencyBasic           ! fraction of pollen emission  
                          
    end function fu_make_allergen_in_new_pollen

    !================================================================================

    real function fu_make_allergen_in_rdy_pollen(rdyPollen, rateBasic) ! , &
                                     ! & fHumid, fWind10m, fConvVelocity, fPrecip)
      !
      ! Returns the allergen growth in waiting pollen 
      ! to be called before new pollen is put to buffer
      ! This far just linear growth in time
      ! Meteo dependence TO BE PARAMETERISED
      !
      implicit none

      ! Imported parameters
      real, intent(in) :: rdyPollen, rateBasic ! , fHumid, fWind10m, fConvVelocity, fPrecip
     
      fu_make_allergen_in_rdy_pollen = rdyPollen * rateBasic       ! fraction of pollen emission  

    end function fu_make_allergen_in_rdy_pollen
    
     !================================================================================

    real function fu_pollen_emis_from_buffer(fHumid, fWind10m, fConvVelocity, fPrecip, &
                                         & fLowHumidThresh, fHighHumidThresh, fPrecipThresh, &
                                         & fWindSatur, fWindMaxScale, &
                                         & fMaxReleaseRate_FromBuffer, &
                                         & fPollenRdyToFly, timestep_sec, rType) result(fMass)
      !
      ! Puts the ripe pollen from buffer to emission map
      !
      implicit none

      ! Imported parameters
      real, intent(in) :: fHumid, fWind10m, fConvVelocity, fPrecip, &
                        & fLowHumidThresh, fHighHumidThresh, fPrecipThresh, &
                        & fWindSatur, fWindMaxScale, &
                        & fMaxReleaseRate_FromBuffer, &
                        & fPollenRdyToFly, timestep_sec
      integer, intent(in) :: rType
      !
      ! Main components taken into account one-by-one:
      ! - wind forcing with saturation
      ! - precipitation dumping
      ! - humidity correction
      !
      ! For intermediate humidity bwteen two thresholds, decrease the emission rate
      ! Note: for above-the-high-threshold values this function is not called at all
      !
      if(rType == relInst)then
          fMass = fPollenRdyToFly
          return
      endif
      if(fHumid <= fLowHumidThresh)then  
        fMass = (fWindMaxScale - exp(-(sqrt(fWind10m**2 + (1.2*fConvVelocity)**2)) / fWindSatur)) * &   ! wind
           & (1. - fPrecip / fPrecipThresh)                                           ! precip
      else
        fMass = (fWindMaxScale - exp(-(sqrt(fWind10m**2 + (1.2*fConvVelocity)**2)) / fWindSatur)) * &   ! wind
           & (1. - fPrecip / fPrecipThresh) * &                                       ! precip
           & (fHighHumidThresh - fHumid) / (fHighHumidThresh - fLowHumidThresh)       ! humidity
      endif
      
      if(rType == relExp)then  ! 2-nd order implicit scheme
          fMass = fPollenRdyToFly *(1.0-exp(- fMaxReleaseRate_FromBuffer * timestep_sec * fMass))
      elseif(rType == relLim)then  ! Limited maximum release rate 
          fMass = fPollenRdyToFly * min(1.0, fMass)
          if(fMass > fMaxReleaseRate_FromBuffer * timestep_sec) fMass = fMaxReleaseRate_FromBuffer * timestep_sec
      else
          call set_error('Unknown emission type', 'fu_pollen_emis_from_buffer')
      endif
      
    end function fu_pollen_emis_from_buffer


    !================================================================================

    real function fu_allergen_emis_from_buffer(bufAl, emsPol, bufPol) result(fMass)
      !
      ! Puts the allergen from buffer to emission map proportionally to released pollen fraction
      !
      implicit none

      ! Imported parameters
      real, intent(in) :: bufAl, emsPol, bufPol
 
      fMass = bufAl * emsPol / (emsPol + bufPol)
      if(fMass > bufAl) fMass = bufAl
      
    end function fu_allergen_emis_from_buffer

    !================================================================================

    real function fu_free_allergen_ems(emsPollen, fractBasic) ! , &
                                     ! & fHumid, fWind10m, fConvVelocity, fPrecip)
      !
      ! Returns the free allergen emission mass 
      ! NO separate buffer for free allergen this far!
      ! This far just fixed fraction of released pollen
      ! Meteo dependence TO BE PARAMETERISED
      !
      implicit none

      ! Imported parameters
      real, intent(in) :: emsPollen, fractBasic ! , fHumid, fWind10m, fConvVelocity, fPrecip
     
      fu_free_allergen_ems = emsPollen * fractBasic           ! fraction of pollen emission  
                         
    end function fu_free_allergen_ems

    !================================================================================

    subroutine update_chill_sum(met_buf, disp_buf, &
                              & indCS, indTempr, &
                              & aver_interval, mdl_now, mdl_timestep, heatsum_params)
      !
      ! Updates the heat sum of all types
      !
      implicit none

      ! Imported parameters
      integer, intent(in) :: indCS, indTempr
      type(Tfield_buffer), pointer :: met_buf, disp_buf
      type(silja_interval), intent(in) :: aver_interval, mdl_timestep
      type(silja_time), intent(in) :: mdl_now
      type(Theat_sum), intent(in) :: heatsum_params

      ! Local variables
      integer :: ixDsp, iyDsp, ixMet, iyMet, iDisp, iMet
      real :: fTempr, fCU

!open (51,file="d:\tmp\heatsum_before.txt")
!iDisp = 0
!do iyDsp = 1, ny_dispersion
!  write(51,'(1000F7.3)') (disp_buf%p2d(indHS)%present%ptr(iDisp + ixDsp), ixDsp = 1, nx_dispersion)
!  iDisp = iDisp + nx_dispersion
!end do
!close(51)
      select case(heatsum_params%iChillSumType)

        case(csSarvasGeneralised)
          do iyDsp = 1, ny_dispersion
            do ixDsp = 1, nx_dispersion
              iDisp = ixDsp + (iyDsp-1) * nx_dispersion
              !
              ! If we are late enough in the course of the year, sum up the
              ! temperature into the chilling units.
              !
              iMet =  fu_grid_index(nx_meteo, ixDsp, iyDsp, pHorizInterpMet2DispStruct)

              fTempr = met_buf%p2d(indTempr)%present%ptr(iMet)

              if(fTempr < heatsum_params%chill_Tmin) cycle       ! too cold

              if(fTempr <= heatsum_params%chill_Topt)then        ! below the optimal temperature
                fCU = (fTempr - heatsum_params%chill_Tmin) / &
                    & (heatsum_params%chill_Topt - heatsum_params%chill_Tmin)
              else                                               ! above the optimal temperature. 
                if(fTempr <= heatsum_params%chill_Tmax .or. heatsum_params%chill_ifNegativeAboveMax)then
                  fCU = (heatsum_params%chill_Tmax - fTempr) / &                ! can be negative
                      & (heatsum_params%chill_Tmax - heatsum_params%chill_Topt)
                else
                  fCU = 0.                      ! too hot
                endif
              endif    ! fTempr > tOpt

              ! Update the chill sum in the grid cell
              ! past, present and future all the same !
              !
              disp_buf%p2d(indCS)%present%ptr(iDisp) = disp_buf%p2d(indCS)%present%ptr(iDisp) + fCU
            end do   ! ix
          end do   ! iy

        case default
          call set_error('Unknown chill sum type:'+fu_str(heatsum_params%iHeatSumType),'update_chill_sum')
          return
      end select
        
      call set_valid_time(disp_buf%p2d(indCS)%present%idPtr, mdl_now)
      call set_field_kind(disp_buf%p2d(indCS)%present%idPtr, accumulated_flag)
      call set_accumulation_length(disp_buf%p2d(indCS)%present%idPtr, &
                                 & fu_valid_time(disp_buf%p2d(indCS)%present%idPtr) - &
                                 & fu_analysis_time(disp_buf%p2d(indCS)%present%idPtr))
      call set_validity_length(disp_buf%p2d(indCS)%present%idPtr, aver_interval)  ! until the next update

call msg('Updated chill sum. Present and future sum:', sum(disp_buf%p2d(indCS)%present%ptr)/fs_dispersion, &
                                                    & sum(disp_buf%p2d(indCS)%future%ptr)/fs_dispersion)
!open (51,file="d:\tmp\heatsum_after.txt")
!iDisp = 0
!do iyDsp = 1, ny_dispersion
!  write(51,'(1000F7.3)') (disp_buf%p2d(indHS)%present%ptr(iDisp + ixDsp), ixDsp = 1, nx_dispersion)
!  iDisp = iDisp + nx_dispersion
!end do
!close(51)
!pause
    end subroutine update_chill_sum

    !================================================================================

    subroutine update_heat_sum(met_buf, disp_buf, &
                             & indHS, indHS_StartDay, indTempr, &
                             & aver_interval, mdl_now, mdl_timestep, heatsum_params)
      !
      ! Updates the heat sum of all types
      !
      implicit none

      ! Imported parameters
      integer, intent(in) :: indHS, indHS_StartDay, indTempr
      type(Tfield_buffer), pointer :: met_buf, disp_buf
      type(silja_interval), intent(in) :: aver_interval, mdl_timestep
      type(silja_time), intent(in) :: mdl_now
      type(Theat_sum), intent(in) :: heatsum_params

      ! Local variables
      integer :: ixDsp, iyDsp, ixMet, iyMet, iDayInYear, iDisp, iMet
      type(Tfield_buffer), pointer :: mb, db
      real :: fTempr, day_fraction, BD, r_T, h, l_h

      mb => met_buf
      db => disp_buf
      iDayInYear = fu_julian_date(now)
      day_fraction = fu_days(aver_interval)

!open (51,file="d:\tmp\heatsum_before.txt")
!iDisp = 0
!do iyDsp = 1, ny_dispersion
!  write(51,'(1000F7.3)') (disp_buf%p2d(indHS)%present%ptr(iDisp + ixDsp), ixDsp = 1, nx_dispersion)
!  iDisp = iDisp + nx_dispersion
!end do
!close(51)
      select case(heatsum_params%iHeatSumType)

        case(hsDegreeDay, hsDegreeHour)
          do iyDsp = 1, ny_dispersion
            do ixDsp = 1, nx_dispersion
              iDisp = ixDsp + (iyDsp-1) * nx_dispersion
              ! Too early to do anything?
              if(iDayInYear < disp_buf%p2d(indHS_StartDay)%present%ptr(iDisp))cycle
              !
              ! If we are late enough in the course of the year, sum up the
              ! temperature into the degree-day or degree-hour. 
              iMet =  fu_grid_index(nx_meteo, ixDsp, iyDsp, pHorizInterpMet2DispStruct)

              if(heatsum_params%iHeatSumType == hsDegreeDay)then
                fTempr = disp_buf%p2d(indTempr)%present%ptr(iDisp)
              else
                fTempr = met_buf%p2d(indTempr)%present%ptr(iMet)
              endif

              if(fTempr >= disp_buf%p2d(indTemprCutOff)%present%ptr(iDisp)) &
                  & disp_buf%p2d(indHS)%present%ptr(iDisp) = &                       ! past, present and future all the same !
                                   & disp_buf%p2d(indHS)%present%ptr(iDisp) + &
                                   & fTempr - disp_buf%p2d(indTemprCutOff)%present%ptr(iDisp)
            end do   ! ix
          end do   ! iy
      
        case(hsBioDay)
          do iyDsp = 1, ny_dispersion
            do ixDsp = 1, nx_dispersion
              iDisp = ixDsp + (iyDsp-1) * nx_dispersion
              ! Too early to do anything?
              if(iDayInYear < disp_buf%p2d(indHS_StartDay)%present%ptr(iDisp))cycle
              iMet =  fu_grid_index(nx_meteo, ixDsp, iyDsp, pHorizInterpMet2DispStruct)
              ! T2m
              fTempr = met_buf%p2d(indTempr)%present%ptr(iMet)
              ! Biotime
              BD = disp_buf%p2d(indHS)%present%ptr(iDisp)
              ! temperature response
              r_T = fu_btime_tempr_response(fTempr)
              h = fu_hours(fu_daylength(fu_lat_geographical_from_grid(real(ixDsp), real(iyDsp), &
                                                                    & dispersion_grid), &
                                     & mdl_now))
              l_h = fu_btime_daylen_response(h, BD)    
              disp_buf%p2d(indHS)%present%ptr(iDisp) = BD + l_h * r_T * day_fraction
            end do   ! ix
          end do   ! iy

        case(hsSigmoidPeriodUnits)
          do iyDsp = 1, ny_dispersion
            do ixDsp = 1, nx_dispersion
              iDisp = ixDsp + (iyDsp-1) * nx_dispersion
              !
              ! Too early to do anything?
              if(iDayInYear < disp_buf%p2d(indHS_StartDay)%present%ptr(iDisp))cycle
              !
              ! If we are late enough in the course of the year, sum up the heatsum
              !
              iMet =  fu_grid_index(nx_meteo, ixDsp, iyDsp, pHorizInterpMet2DispStruct)
              !
              ! Sigmoid is multiplied with a coefficient that pushes it to zero when cut-off is reached
              !
              if(met_buf%p2d(indTempr)%present%ptr(iMet) >= heatsum_params%TemprCutOff)then
                disp_buf%p2d(indHS)%present%ptr(iDisp) = disp_buf%p2d(indHS)%present%ptr(iDisp) + &
                         & day_fraction * &
                         & fu_sigmoid_hsum_response_tempr(met_buf%p2d(indTempr)%present%ptr(iMet), &
                                                        & heatsum_params)
              endif
            end do   ! ix
          end do   ! iy

        case default
          call set_error('Unknown heatsum type:'+fu_str(heatsum_params%iHeatSumType),'update_heat_sum')
          return
      end select
        
      call set_valid_time(disp_buf%p2d(indHS)%present%idPtr, mdl_now)
      call set_field_kind(disp_buf%p2d(indHS)%present%idPtr, accumulated_flag)
      call set_accumulation_length(disp_buf%p2d(indHS)%present%idPtr, &
                                 & fu_valid_time(disp_buf%p2d(indHS)%present%idPtr) - &
                                 & fu_analysis_time(disp_buf%p2d(indHS)%present%idPtr))
      call set_validity_length(disp_buf%p2d(indHS)%present%idPtr, aver_interval)  ! until the next update

call msg('Updated heat sum. Present and future sum:', sum(disp_buf%p2d(indHS)%present%ptr)/fs_dispersion, &
                                                    & sum(disp_buf%p2d(indHS)%future%ptr)/fs_dispersion)
!open (51,file="d:\tmp\heatsum_after.txt")
!iDisp = 0
!do iyDsp = 1, ny_dispersion
!  write(51,'(1000F7.3)') (disp_buf%p2d(indHS)%present%ptr(iDisp + ixDsp), ixDsp = 1, nx_dispersion)
!  iDisp = iDisp + nx_dispersion
!end do
!close(51)
!pause
    end subroutine update_heat_sum

    !================================================================================

    real function fu_sigmoid_hsum_response_tempr(fTempr_K, sig_params) result(dHSdt)
      !
      ! Gives the dHSdt function that is to be integrated over time to get the increase of the 
      ! sigoid-based heatsum
      !
      implicit none
      
      ! Improted parameters
      real, intent(in) :: fTempr_K
      type(Theat_sum), intent(in) :: sig_params
      
      dHSdt = sig_params%maxVal * &
            & ((fTempr_K - sig_params%TemprCutOff) / (fTempr_K - sig_params%TemprCutOff + &
                                                    & sig_params%smoother)) ** sig_params%sm_power / &
            & (1.0 + EXP (sig_params%exp_power * (sig_params%TemprMidPoint - fTempr_K)))
      
    end function fu_sigmoid_hsum_response_tempr
                             
    !================================================================================

    real function fu_btime_tempr_response(T)
      ! Temperature response for biotime accumulation
      implicit none
      real, intent(in) :: T
      
      fu_btime_tempr_response = 0.
      
      select case(srcPollen%iTRfunc)
      case(trfLinear)      ! Linear 
        if (T >= srcPollen%loTemp .and. T <= srcPollen%optTemp)then
          fu_btime_tempr_response = (T - srcPollen%loTemp) / (srcPollen%optTemp - srcPollen%loTemp)
        else if (T > srcPollen%optTemp .and. T <= srcPollen%hiTemp)then
          fu_btime_tempr_response =(srcPollen%hiTemp - T) / (srcPollen%hiTemp - srcPollen%optTemp)
        endif
      
      case(trfBeta)!Beta (Yan & Hunt, 1999)
        if (T > srcPollen%loTemp .and. T < srcPollen%hiTemp)then
          fu_btime_tempr_response = (srcPollen%hiTemp - T) / &
                                        &(srcPollen%hiTemp - srcPollen%optTemp) * &
                                        &((T - srcPollen%loTemp) / &
                                         &(srcPollen%optTemp - srcPollen%loTemp)) ** &
                                        &((srcPollen%optTemp - srcPollen%loTemp) / &
                                         &(srcPollen%hiTemp - srcPollen%optTemp))
        endif
      case default
        call set_error('Unknown temperature response function','fu_btime_tempr_response')
      end select
      
    end function fu_btime_tempr_response   
 
  !========================================================================================

    real function fu_btime_daylen_response(daylen, bioday)    
      ! Photoperiod response for biotime accumulation
      implicit none
      real, intent(in) :: daylen, bioday
      fu_btime_daylen_response = 1.
      if(bioday > 11.5 .and. bioday < 16.)then
        if(daylen > srcPollen%photoperiod)then
          fu_btime_daylen_response = exp((daylen-srcPollen%photoperiod) * log(1.-0.5))
        endif 
      else if(bioday > 16. .and. bioday < 20.5)then
        if(daylen > srcPollen%photoperiod)then
          fu_btime_daylen_response = exp((daylen-srcPollen%photoperiod) * log(1.-0.6))
        endif 
      endif
    end function fu_btime_daylen_response


  end subroutine compute_emission_for_pollen   !***********************************************

  
  !********************************************************************************************
  !********************************************************************************************
  !********************************************************************************************
  
  subroutine compute_HS_sigmoid_parameters(fHSTemprSatur, fHSMaxRate, fExpTempResp, pSigmoid_params)
    !
    ! Computes parameters for sigmoid shape resembling the linear heat-sum function where
    ! possible.
    ! Ordinary sigmoid has 3 params: max value, mid-point, and slope
    ! Heatsum additionally needs cut-off, which cannot be below zero (no ripening in frost)
    ! To simplify the life of user and further development, these parameters are derived here
    ! from the more physically clear quantities.
    !
    implicit none
  
    real, intent(in) :: fHSTemprSatur, fHSMaxRate, fExpTempResp
    type(Theat_sum), intent(out) :: pSigmoid_params
    
call set_error('Obsolete sigmoid function setup','compute_HS_sigmoid_parameters')
call unset_error('compute_HS_sigmoid_parameters')

    ! max value is the max rate of HS accumulation, degree-day per day or alike
    pSigmoid_params%maxVal = fHSMaxRate*1.1  ! +1-2 DD/day to handle T > fHSTemprSatur
    
    ! cut-off is above freesing temperature but close to it. Significantly positive only for
    ! very warm-loving species
    pSigmoid_params%TemprCutOff = 0.05 * (fHSTemprSatur-273.15) + 273.15
    pSigmoid_params%smoother = 5.0
    pSigmoid_params%sm_power = 1.0
    
    ! Mid-point is in the middle between the cut-off and saturation, a bit to the right
    ! to handle the asymmetric fade-in curve
    pSigmoid_params%TemprMidPoint = (pSigmoid_params%TemprCutOff + fHSTemprSatur-273.15) / 2. + 1.0 +273.15
    
!    ! Slope should be 1 around the mid-point to meet the linear heatsum curve
!    !
!    pSigmoid_params%exp_power = 5.0 / pSigmoid_params%maxVal
    
    pSigmoid_params%exp_power = fExpTempResp  ! from ini file
    
  end subroutine compute_HS_sigmoid_parameters
  
  !********************************************************************************************
  
  subroutine compute_HS_sigmoid_parameters2(fHS_midpoint_tempr, fHS_cutoff_tempr, &
                                          & fSmoother, fSM_power, &
                                          & fHS_max_rate, fHS_exp_temperature_resp, &
                                          & pSigmoid_params)
    !
    ! Computes parameters for sigmoid shape resembling the linear heat-sum function where
    ! possible.
    ! Ordinary sigmoid has 3 params: max value, mid-point, and slope
    ! Heatsum additionally needs cut-off, which cannot be below zero (no ripening in frost)
    ! To simplify the life of user and further development, these parameters are derived here
    ! from the more physically clear quantities.
    !
    implicit none
  
    real, intent(in) :: fHS_midpoint_tempr, fHS_cutoff_tempr, fSmoother, fSM_power, &
                      & fHS_max_rate, fHS_exp_temperature_resp
    type(Theat_sum), intent(out) :: pSigmoid_params
    
call msg('Sigmoid parameters new setup')

    ! max value is the max rate of HS accumulation, degree-day per day or alike - for very high temperature
    pSigmoid_params%maxVal = fHS_max_rate
    
    ! cut-off better be above freesing temperature but close to it. Significantly positive only for
    ! very warm-loving species
    pSigmoid_params%TemprCutOff = fHS_cutoff_tempr
    pSigmoid_params%smoother = fSmoother
    pSigmoid_params%sm_power = fSM_power
    
    ! Mid-point is in the middle of the logistic function rise range 
    !
    pSigmoid_params%TemprMidPoint = fHS_midpoint_tempr
    
    pSigmoid_params%exp_power = fHS_exp_temperature_resp  ! from ini file
    
  end subroutine compute_HS_sigmoid_parameters2
  

  !********************************************************************************************
  !********************************************************************************************

  logical function fu_pollen_emis_owned_quantity(pollen_src, quantity)
    !
    ! Checsk whether the specific quantity is handled by the emission sources exclusively
    ! Such quantities may be related to emission species somewhere else
    !
    implicit none

    ! Imported parameters
    type(silam_pollen_source), intent(in) :: pollen_src
    integer, intent(in) :: quantity
    !
    ! The pollen source has plenty of own quantities
    !
    select case(quantity)
      case(heatsum_flag, chillsum_flag, & 
         & start_calday_threshold_flag, end_calday_threshold_flag, &
         & start_heatsum_threshold_flag, end_heatsum_threshold_flag, &
         & growth_season_start_day_flag, heatsum_cutoff_tempr_flag, &
         & temperature_threshold_flag, daily_temp_threshold_flag, soil_moisture_threshold_flag, &
         & pollen_left_relative_flag, pollen_total_per_m2_flag, pollen_correction_flag, pollen_rdy_to_fly_flag, &
         & heatsum_start_end_diff_flag, calday_start_end_diff_flag, plant_growth_flag)
        fu_pollen_emis_owned_quantity = .true.
      case default
        fu_pollen_emis_owned_quantity = .false.
    end select

  end function fu_pollen_emis_owned_quantity


  !*****************************************************************

  integer function fu_source_id_nbr_of_pollen_src(pollen_src)
    !
    ! Returns the source number. Reason: all sources are enumerated
    ! sequencially, so that the source can, in fact, be refered by its
    ! number without other information. But so far this number is 
    ! copied to the particles in the pollution cloud. 
    !
    ! NOTE. One and only index may be reasonable. The other one MUST be
    ! negative or zero
    !
    implicit none

    ! Imported parameters with intent IN
    type(silam_pollen_source), intent(in) :: pollen_src

    ! Stupidity check
    if(.not.(pollen_src%defined == silja_false))then
      call set_error('Undefined source given','fu_source_id_nbr_of_pollen_src')
      return
    endif
    fu_source_id_nbr_of_pollen_src = pollen_src%id_nbr

  end function fu_source_id_nbr_of_pollen_src


  !*************************************************************************

  integer function fu_source_nbr_of_pollen_src(pollen_src)
    !
    ! Returns the source number. Reason: all sources are enumerated
    ! sequencially, so that the source can, in fact, be refered by its
    ! number without other information. But so far this number is 
    ! copied to the particles in the pollution cloud. 
    !
    ! NOTE. One and only index may be reasonable. The other one MUST be
    ! negative or zero
    !
    implicit none

    ! Imported parameters with intent IN
    type(silam_pollen_source), intent(in) :: pollen_src

    ! Stupidity check
    if(.not.(pollen_src%defined == silja_false))then
      fu_source_nbr_of_pollen_src = pollen_src%src_nbr
    else
      fu_source_nbr_of_pollen_src = int_missing
      call set_error('Undefined source given','fu_source_nbr_of_pollen_src')
      return
    endif

  end function fu_source_nbr_of_pollen_src


  !*************************************************************************

  subroutine typical_species_cnc_pollen(pollen_src, &
                                      & species, nSpecies, arConc)  ! output
    !
    ! Guesses a typical level of concentration and divides it with the given accuracy factor
    !
    implicit none

    ! Imported parameters
    type(silam_pollen_source), intent(in) :: pollen_src
    type(silam_species), dimension(:), pointer :: species
    integer, intent(out) :: nSpecies
    real, dimension(:), pointer :: arConc
    integer :: iSpecies

    species => pollen_src%species
    nSpecies = pollen_src%nSpecies
    
    do iSpecies = 1, nSpecies
      if(index(fu_name(fu_material(species(ispecies))), 'ALLERGEN') > 0)then
        arConc(iSpecies) = 1e-12     ! no idea, actually
      else ! everything else emitted as pollen
        arConc(iSpecies) = 20.0   ! Not too high, neither too low, suitable for nearly all taxa
      endif
    enddo
      
  end subroutine typical_species_cnc_pollen


  !*************************************************************************

  function fu_pollen_source_name(pollen_src)result(chNm)
    implicit none
    type(silam_pollen_source), intent(in) :: pollen_src
    character(len=clen) :: chNm
    chNm = pollen_src%src_nm
  end function fu_pollen_source_name


  !*************************************************************************

  subroutine report_pollen_src(pSrc)
    implicit none
    type(silam_pollen_source), intent(in) :: pSrc
    integer :: iSpecies, iTmp
    
    call msg('')
    call msg('=========================================================') 
    call msg('Reporting pollen source:' + pSrc%src_nm +';'+ pSrc%sector_nm)
    
    call msg('Emitted species:')
    do iSpecies = 1, pSrc%nSpecies
      call report(pSrc%species(iSpecies))
    end do
    
    call msg('Emitted standard total and uncertainty:', &
            & pSrc%standardPollenTotal, pSrc%fUncertainty_tot_pollen_relat)
    if(pSrc%ifHSGrowth) call msg('Heatsum dependent plant size')
    if(pSrc%ifSWGrowth) call msg('Soil moisture dependent plant size')
    
    if(.not. pSrc%indPolAlrg==int_missing)call msg('Pollen allergen growth rate:',  pSrc%alrgGrowthRate)
    if(.not. pSrc%indFreeAlrg==int_missing)call msg('Free allergen fraction:',  pSrc%freeAlrgFract)
    
    call msg('Flowering thresholds:')
    if(pSrc%ifStartHSThr)call msg('Start heatsum')
    if(pSrc%ifStartCDThr)call msg('Start calendar day')
    if(pSrc%ifEndHSThr)call msg('End heatsum')
    if(pSrc%ifEndCDThr)call msg('End calendar day')
    if(pSrc%ifStartEndGammaThr) call msg('Start-end gamma-distribution percentiles')
    if(pSrc%ifTempThr)call msg('temperature')
    if(pSrc%ifDayTempThr)call msg('daily temperature')
    if(pSrc%ifSWthr)call msg('soil moisture')
    
    call report(pSrc%UncertaintyParams)
 
    if(pSrc%ifStartCDThr .or. pSrc%ifEndCDThr) &
        & call msg('Calendar day shift and uncertainty of the start:', fu_days(pSrc%timeMapShift), &
                                                                  & pSrc%fUncertainty_CD_days_start)
    if(pSrc%ifStartEndGammaThr)then
      call msg('Gamma shape parameters')
      call msg('Exponent scale and number of power terms:', pSrc%gf_beta, pSrc%gf_nTerms)
      call msg('Power terms scales:',pSrc%gf_scales)
      call msg('Power terms start moments, relative:', pSrc%gf_start_times)
      call msg('Powre terms powers:', pSrc%gf_powers)
    endif
    
    if(pSrc%ifStartHSThr .or. pSrc%ifEndHSThr)then
       select case(pSrc%heatsum_params%iHeatSumType)
       case(hsDegreeDay)
         call msg('Heatsum type: degree day')
       case(hsDegreeHour)
         call msg('Heatsum type: degree hour')
       case(hsBioDay)
         call msg('Heatsum type: bioday')
         call msg('lo temp, hi temp:', pSrc%LoTemp, pSrc%HiTemp)
         call msg('opt temp, photoperiod:', pSrc%OptTemp, pSrc%photoperiod)
         call msg('Temperature response function',  pSrc%iTRfunc)
       case(hsSigmoidPeriodUnits)
         call msg('Heatsum type: sigmoid')
         call msg('Sigmoid max daily rate and exponent power', pSrc%heatsum_params%maxVal, &
                                                             & pSrc%heatsum_params%exp_power)
         call msg('Sigmoid cut-off amd mid-point temperature', pSrc%heatsum_params%TemprCutOff, &
                                                             & pSrc%heatsum_params%TemprMidPoint)
         call msg('Sigmoid smoother and smoothing power', pSrc%heatsum_params%smoother, &
                                                        & pSrc%heatsum_params%sm_power)
       case default
         call msg('Unknown heatsum type:', pSrc%heatsum_params%iHeatSumType)
       end select
       call msg('Uncertainty of the start:', pSrc%fUncertainty_HS_relative_start)
    endif
    
    select case(pSrc%iRipeningType)
    case(prCDLinear)
      call msg('Pollen ripening type: calendar day; linear')
    case(prHSLinear)
      call msg('Pollen ripening type: heatsum; linear')
    case(prCDNormal)
      call msg('Pollen ripening type: calendar day; normal')
    case(prHSNormal)
      call msg('Pollen ripening type: heatsum; normal')
    case(prCDNormDrn)
      call msg('Pollen ripening type: calendar day; normal; diurnal')
      call msg('First diurnal mode time, sigma: ', pSrc%dayMid1, pSrc%daySgm1)
      call msg('Second diurnal mode time, sigma: ', pSrc%dayMid2, pSrc%daySgm2)
      call msg('First diurnal mode fraction: ', pSrc%dayFrac1)
    case(prCDGammaWTails)
      call msg('Pollen ripening type: gamma-distribution with tails')
    case default
      call msg('Unknown ripening type:', pSrc%iRipeningType)
    end select
      
    call msg('Pollen release parameters:')
      select case(pSrc%iReleaseType)
    case(relExp)
      call msg('Pollen release type: exponential')
      call msg('MaxReleaseRate_FromBuffer', pSrc%fBufReleaseLimit)
    case(relLim)
      call msg('Pollen release type: limited')
      call msg('MaxReleaseRate_FromBuffer', pSrc%fBufReleaseLimit)
    case(relInst)
      call msg('Pollen release type: instant')
      case default
      call msg('Unknown release type:', pSrc%iReleaseType)
    end select
    
    if(.not. pSrc%iReleaseType == relInst)then
      call msg('low_humidity_threshold', pSrc%LowHumidThresh)
      call msg('high_humidity_threshold', pSrc%HighHumidThresh)
      call msg('precipitation_threshold', pSrc%precipThreshold)
      call msg('wind_speed_saturation_level', pSrc%windSaturation)
      call msg('wind_speed_max_impact', pSrc%windMaxScale)
    endif
    
    if(pSrc%ems2wholeABL)then
      call msg('Emission to whole abl')
    else
      call msg('Emission to layer between', pSrc%fEmissionBottom, pSrc%fEmissionTop)
    endif
    
    call msg('End of pollen source report')
    call msg('=========================================================') 
    call msg('')
    
  end subroutine report_pollen_src


END MODULE source_terms_pollen



