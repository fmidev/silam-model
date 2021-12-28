MODULE source_terms_time_params

  ! This module contains the general description of the time parameter of the SILAM
  ! source-term. It describes the temporal distribution of the release to atmosphere, 
  ! and the amount of the species released.
  !
  ! In SILAM the source term consists of an unlimited number of
  ! release sources, each describing one release
  !
  ! In this module we define: temporal and chemical characteristics of the source terms
  !
  ! The basic term is the time_slot, which is the exact composition, rate, and physical
  ! characteritics of the release - all taken at specific time.
  !
  ! Units: not necessarily SI
  ! 
  ! All degrees: real numbers in degrees (degrees and parts per
  ! hundred/parts per thousand , NOT in degrees and minutes)
  ! 
  ! Language: ANSI standard Fortran 90
  !
  ! Code owner: Mikhail Sofiev, FMI
  ! 

  USE dispersion_server !chemistry_server
!  use silam_namelist
  
  IMPLICIT NONE

  ! The public functions and subroutines available in this module:

  public set_missing
  public fu_time
  public enlarge_param_vector
  public decode_time_parameter_v4
  public decode_time_parameter_v5
  public decode_time_params_volc_v1
  public integrate_dynamic_slot_descr
  public store_time_param_as_namelist
  public count_descriptor_names
  public add_descriptor_names
  public complete_descriptor_names
  public get_source_aer_species
  public getTimeSlots_from_params
  public get_overlapping_slots
  public add_input_needs_meteo_dependence
  public fu_meteodep_model
  public get_temp_dep_from_namelist

  private set_release_parameters
  private reset_release_parameters
  private fu_time_from_time_param
  private set_missing_time_param
  private fu_descr_name
  private store_time_param_as_nl_params
  private store_time_param_as_nl_time
  private fu_meteodep_model_time_params


  interface fu_time
    module procedure fu_time_from_time_param
  end interface

  interface set_missing
    module procedure set_missing_time_param
  end interface

  interface store_time_param_as_namelist
    module procedure store_time_param_as_nl_params
    module procedure store_time_param_as_nl_time
  end interface

  interface fu_meteodep_model
    module procedure fu_meteodep_model_time_params
  end interface

  ! Release cell types:  
  INTEGER, PARAMETER, PUBLIC :: point_source = 2101
  INTEGER, PARAMETER, PUBLIC :: line_source = 2102   ! void so far
  INTEGER, PARAMETER, PUBLIC :: area_source = 2103
  INTEGER, PARAMETER, PUBLIC :: fire_source = 2104
  INTEGER, PARAMETER, PUBLIC :: bomb_source = 2105
  INTEGER, PARAMETER, PUBLIC :: sea_salt_source = 2106
  INTEGER, PARAMETER, PUBLIC :: pollen_source = 2107
  INTEGER, PARAMETER, PUBLIC :: bio_voc_source = 2108
  INTEGER, PARAMETER, PUBLIC :: wind_blown_dust_source = 2109
  INTEGER, PARAMETER, PUBLIC :: dms_source = 2110
  INTEGER, PARAMETER, PUBLIC :: volcano_source = 2111
 

  ! Vertical distribution:
  INTEGER, PARAMETER, PUBLIC :: vertically_even_distribution = 1
  INTEGER, PARAMETER, PUBLIC :: vertically_normal_distribution = 2
  INTEGER, PARAMETER, PUBLIC :: vertically_single_level = 3

  !
  ! First, introduce the release parameter. This class describes the releaase rate 
  ! and composition at some time moment. Note that the source may be in some non-
  ! standard unit - depending only on user definitions
  !
  type silam_release_parameters
!    private 
    type(silja_time) :: time       ! Time moment when the values are defined
    real, dimension(:), allocatable :: rate_descr_unit     ! (nDescriptors): release rate but SI unit,
    real :: xy_size                ! horizontal linear size - for plume rise, stack diameter
    real :: z_size                 ! vertical size in m, thickness of the layerDynamic
    real :: z_velocity              ! gas vertical velocity - for plume rise
    real :: tempr                  ! gas temperature - for plume rise
    TYPE(silja_level) :: layerDynamic   !must be a layer with bottom and top
    INTEGER :: vertical_distribution
    real, dimension(:), allocatable :: levFractDispVert, fzDisp
    integer :: nDescriptors
  end type silam_release_parameters

  !
  ! Some sources can depend on meteorological conditions. 
  !
  type TmeteoDepRules
    integer :: indMeteoDepModelSwitch  ! Various dependenciess can be considered
    !
    ! Heating as a function of temperature deficit: linear regression with saturation.
    ! 
    ! EmissionScaling = fMinEmissionScaling + fTemprDeficit_2_emis * max(0, T2mDaily - fT0)
    !
    real :: fT0    !  Saturation: no heating when T2m >fT0 fT0=288 should be fine for all:) 
    real :: fdEdT   ! slope, default. Dynamic part of dependence.
    real :: fMinEmission         ! intercept: non-HDD part
    !
    ! If the source emission should be injected into ABL homogeneously from surface to the top
    !
    logical :: ems2wholeABL=.false.
  end type TmeteoDepRules

  type (TmeteoDepRules), public, parameter :: meteoDepRulesMissing = &
      & TmeteoDepRules(int_missing, real_missing, real_missing, real_missing)
  
  integer, public, parameter :: max_descriptors = 30

  integer, private, parameter :: emis_heating_vs_T_linregr_mdl = 3101 ! linear regres. with saturation

  real, dimension(:), private, save, pointer :: fldDailyTempr2m => Null()
  
CONTAINS

  !====================================================================
  !====================================================================
  !
  !  Release parameter class functions
  !
  !====================================================================
  !====================================================================

  subroutine set_release_parameters(the_params, & 
                                  & ParTime, & 
                                  & fRate_descr_unit, &  ! emission rate of descriptors
                                  & xy_size, &
                                  & bottom_lev, top_lev, vert_distrib, &
                                  & z_velocity,& 
                                  & fTempr, &
                                  & cocktail_descr_lst_of_param, &  ! order corresponds to fRate_descr_unit
                                  & cocktail_descr_lst_of_source, & ! all descriptors of the source
                                  & nDescrMax) !, &
!                                  & indDescriptorOfParLine, &
!                                  & nDescriptorInParLine)
    !
    ! This function creates one alias of the release parameters class.
    !
    ! All units: SI
    !
    implicit none

    ! Return value of the function
    type(silam_release_parameters), intent(inout), target :: the_params
    
    ! Imported values with the intent IN
    type(silja_time), intent(in) :: ParTime
    real, intent(in) :: xy_size, z_velocity, fTempr
    real, dimension(:), intent(in) :: fRate_descr_unit   ! (up to nDescrMax)
    type(silja_level), intent(in) :: bottom_lev, top_lev
    integer, intent(in) :: vert_distrib, nDescrMax
    type(Tcocktail_descr), dimension(:), intent(in) :: cocktail_descr_lst_of_param, & ! same as fRate_descr_unit
                                                     & cocktail_descr_lst_of_source
    ! Local variables
    integer :: iTmp, jTmp
    type(silam_release_parameters), pointer :: paramsPtr

    paramsPtr => the_params  ! for debugging purposes

    if(fu_fails(defined(ParTime), 'Undefined time given','set_release_parameters'))return
    the_params%time = ParTime

    if(the_params%nDescriptors > nDescrMax)then  ! used to be >=, do not know why. MAS 30.12.2010
      call msg('nDescriptors > nDescrMax:', the_params%nDescriptors, nDescrMax)
      call set_error('Wrong descriptor vector size','set_release_parameters')
    endif

    if(the_params%nDescriptors <= 0)then
      the_params%nDescriptors = 0.
      allocate(the_params%rate_descr_unit(nDescrMax), stat=iTmp)
      if(fu_fails(iTmp == 0, 'Failed memory allocation for rate of parameter line','set_release_parameters'))return
      the_params%rate_descr_unit(1:nDescrMax) = 0.0
    endif
    the_params%xy_size = xy_size
    the_params%z_velocity = z_velocity
    the_params%tempr = fTempr
    the_params%nDescriptors = nDescrMax

    if(fu_fails(defined(bottom_lev), 'Undefined bottom level given','set_release_parameters'))return
    if(fu_fails(defined(top_lev), 'Undefined top level given','set_release_parameters'))return

    if(fu_cmp_levs_eq(bottom_lev, top_lev))then
      the_params%vertical_distribution = vertically_single_level
    else
      the_params%vertical_distribution = vert_distrib
    endif
    the_params%layerDynamic = fu_set_layer_between_two(top_lev, bottom_lev)
    the_params%z_size = fu_layer_thickness_m(the_params%layerDynamic)
    !
    ! Now store the descriptor rates. Note that the descriptors are the same for all params
    ! and stored elsewhere but rates are the features of the params
    !
    do iTmp = 1, min(nDescrMax, size(cocktail_descr_lst_of_param))
      if(.not. defined(cocktail_descr_lst_of_param(iTmp)))exit
      !
      ! Find the right place for each descriptor in the list of source descriptors.
      ! The rate value should be stored into that place
      !
      do jTmp = 1, nDescrMax
        if(cocktail_descr_lst_of_param(iTmp) == cocktail_descr_lst_of_source(jTmp))then
          the_params%rate_descr_unit(jTmp) = fRate_descr_unit(iTmp)
          exit
        endif
      end do
    end do

    if(allocated(the_params%levFractDispVert)) deallocate(the_params%levFractDispVert)
    if(allocated(the_params%fzDisp)) deallocate(the_params%fzDisp)

  end subroutine set_release_parameters


  !*************************************************************************

  subroutine set_missing_time_param(params, nParams)
    !
    ! Nullifies the time parameter making it ready to accept the values viia set_release_parameters
    ! Note that the same parameter can be set several times, so several of the below values are 
    ! for repeated use afterwards
    !
    implicit none

    ! Imported parameter
    type(silam_release_parameters), dimension(:), intent(inout) :: params
    integer, intent(in) :: nParams

    ! Local vairable
    integer :: iParam

    do iParam = 1, nParams
      params(iParam)%time = time_missing
      if(allocated(params(iParam)%rate_descr_unit)) deallocate(params(iParam)%rate_descr_unit)
      params(iParam)%xy_size = real_missing
      params(iParam)%z_velocity = real_missing
      params(iParam)%tempr = real_missing
      params(iParam)%layerDynamic = level_missing
      params(iParam)%vertical_distribution = int_missing
      if(allocated(params(iParam)%levFractDispVert)) deallocate(params(iParam)%levFractDispVert)
      if(allocated(params(iParam)%fzDisp)) deallocate(params(iParam)%fzDisp)
      params(iParam)%nDescriptors = 0
    end do
  end subroutine set_missing_time_param


  !**************************************************************************

  subroutine reset_release_parameters(params, nParams)
    !
    ! Similar to set_missing but does not destroy the allocated pointers,
    ! only sets them to zero values
    !
    implicit none
    ! Imported parameter
    type(silam_release_parameters), dimension(:), pointer :: params
    integer, intent(in) :: nParams

    ! Local vairable
    integer :: iParam, iDescr

    do iParam = 1, nParams
      params(iParam)%time = time_missing
      params(iParam)%rate_descr_unit = 0.0
      params(iParam)%xy_size = real_missing
      params(iParam)%z_velocity = real_missing
      params(iParam)%tempr = real_missing
      params(iParam)%layerDynamic = level_missing
      params(iParam)%vertical_distribution = int_missing
      params(iParam)%levFractDispVert = 0.0
      params(iParam)%fzDisp = 0.0
!      do iDescr = 1, params(iParam)%nDescriptors
!        call set_missing(params(iParam)%cocktail_descr(iDescr))
!      end do
!      params(iParam)%nDescriptors = 0
    end do

  end subroutine reset_release_parameters


  !*************************************************************************

  function fu_time_from_time_param(param) result(time)
    implicit none
    type(silja_time) :: time
    type(silam_release_parameters), pointer :: param
    time = param%time
  end function fu_time_from_time_param


  !*************************************************************************

  subroutine decode_time_parameter_v4(chContent, &  ! parameter line
                                    & iParam, &     ! index of the time slot to decode
                                    & chUnitVert, & ! unit of vertical axis
                                    & chSrcRateUnit, & ! unit of the release as in the source file
                                    & params, &     ! vector of time slots
                                    & cocktail_descr_lst_of_src, & ! list of source descriptors
                                    & iSrcType, &   ! point/area source
                                    & nDescrMax, &  ! absolute max of what can be
                                    & indDescriptorOfParLine, & ! OUT: index of par_str descriptors in param descr array
                                    & nDescriptorInParLine, &  ! OUT: number of such descriptors
                                    & chSubstNm, iModeNbr, & ! emitted substance and mode number
                                    & ifSpecies)   ! OUT : if species or cocktail is emitted
    !
    ! Based on the name of the content decodes it and writes to the
    ! given p_arc/a_src.
    ! Note: grid must be defined prior to the first value is stored.
    ! For its actual creation we shall use grib_io module.
    ! Do not overlook that geographical parameters are in milliseconds,
    ! in order to keep them integer. 
    ! In SILAM v4, the time parameter string cannot have more than one release and one
    ! cocktail descriptor. To keep the backward compatibility, we have this sub here.
    ! For multi-descriptor-multi-release, see decode_time_parameter_v5
    ! 
    ! Note that the source values are usually given in some derived
    ! units. In order to avoid unnecessary conversion later on, we
    ! convert these units here to the descriptor-default values. Then
    ! the source-term release rate unit is not stored. Instead, we
    ! deal with descriptor default units, which are then explored into
    ! species default units.
    !
    implicit none

    ! Imported parameters
    character(len=*), intent(in) :: chContent, chUnitVert, chSrcRateUnit, chSubstNm
    integer, intent(in) :: iParam, iSrcType, nDescrMax, iModeNbr
    type(silam_release_parameters), dimension(:), pointer :: params
    TYPE(Tcocktail_descr), dimension(:), pointer :: cocktail_descr_lst_of_src
    integer, dimension(:), pointer :: indDescriptorOfParLine ! index of par_str descriptors in param descr array
    integer, intent(out) :: nDescriptorInParLine             ! OUT: number of such descriptors
    logical, dimension(:), allocatable, intent(out) :: ifSpecies

    ! Local variables
    integer io_status, iYr4,iMon,iDay,iHr,iMin,iDur, iDescr
    real :: fRate, height_a, height_b, fSec, xySize, zVelocity, Tempr, fSrcUnit2DescrUnit, fTmp
    logical :: ifNew
    type(silam_sp) :: chCocktail
    character(len=clen) :: chTmp
    character(len=unitNmLen) :: chDescrUnit
    TYPE(Tcocktail_descr) :: cocktail_descr
    TYPE(silja_level) :: top_level, bottom_level

    do while(iParam > size(params)) 
      call enlarge_param_vector(params)
      if(error)return
    end do
    
    chCocktail%sp => fu_work_string()

    !point source:  = 2003 4 15 1 10 0.    400.  1.  925. 900. 10. 273.   PASSIVE_COCKTAIL  
    !area source:   = 2003 4 15 1 10 0.    400.  925. 900.    PASSIVE_COCKTAIL  
    !----------------------------------------
    !
    ! Time variation of the main release parameters:
    !   <rate> - emission rate [mass]/[time], units must correspond to chRateUnit
    !   <bottom> <top>  - vertical boundaries,[hpa] or [m] correspond to chUnit_vertical
    !   <cocktail_name> - name of the nuclear cocktail or name of a cocktail file
    ! Idea is:
    ! 1. If there is only one line starting from word NOW then wall-clock time as 
    !    start time, then duration in minutes follows and then the set of parameters 
    !    They are all considered to be fixed in time
    ! 2. If there are several lines - then the first line determines UTC time start of 
    !    release and sets its initial parameters. 
    !    Each next line determines the parameters of the release at some other time.
    !    All characteristics are interpolated in time between the sequencial moments
    !    The last line determines the end of the release and final parameters
    !
    if(index(chContent,'NOW') > 0 .or. index(chContent,'LAST_METEO_TIME') > 0)then
      if(iSrcType == point_source)then
        read(unit=chContent,fmt=*,iostat=io_status)chTmp,iDur,fRate, xySize, &
                                                 & height_a,height_b,zVelocity, Tempr, chCocktail%sp
      elseif(iSrcType == area_source)then
        read(unit=chContent,fmt=*,iostat=io_status)chTmp,iDur,fRate,height_a,height_b,chCocktail%sp
        xySize = 0.
        zVelocity = 0.
        Tempr = 273.15
      else
        call msg('Unknown source type:',iSrcType)
        call set_error('Unknown source type','decode_time_parameter_v4')
        return
      endif
    else
      if(iSrcType == point_source)then
        read(unit=chContent,fmt=*,iostat=io_status) iYr4,iMon,iDay,iHr,iMin,fSec,fRate, &
                                                  & xySize,height_a,height_b, zVelocity, Tempr, &
                                                  & chCocktail%sp
      elseif(iSrcType == area_source)then
        read(unit=chContent,fmt=*,iostat=io_status) iYr4,iMon,iDay,iHr,iMin,fSec, fRate, &
                                                  & height_a,height_b, chCocktail%sp
        xySize = 0.
        zVelocity = 0.
        Tempr = 273.15
      else
        call msg('Unknown source type:',iSrcType)
        call set_error('Unknown source type','decode_time_parameter_v4')
        return
      endif
    endif
    if(io_status /= 0)then
      call set_error(fu_connect_strings('Strange parameter string:',chContent), &
                   & 'decode_time_parameter_v4_5')
      return
    endif

    chCocktail%sp = fu_descr_name(chCocktail%sp, chSubstNm, iModeNbr)
    !
    ! Negative rates are not allowed
    !
    if(fRate < 0.0)then
      call msg_warning(fu_connect_strings('Negative rate found in:',chContent))
      call set_error('Negative rate is not allowed','decode_time_parameter_v4')
      return
    endif
    !
    ! Vertical boundaries for the just-read line
    !
    if (fu_unit_type(chUnitVert) == size_unit) then
      fTmp = fu_conversion_factor(chUnitVert, "m")
      IF (height_a > height_b) THEN
        top_level = fu_set_constant_height_level(height_a*fTmp)
        bottom_level = fu_set_constant_height_level(height_b*fTmp)
      ELSE
        top_level = fu_set_constant_height_level(height_b*fTmp)
        bottom_level = fu_set_constant_height_level(height_a*fTmp)
      END IF

    elseif (fu_unit_type(chUnitVert) == pressure_unit) then
      fTmp = fu_conversion_factor(chUnitVert, "Pa")
      IF (height_a > height_b) THEN
        top_level = fu_set_pressure_level(height_b*fTmp)
        bottom_level = fu_set_pressure_level(height_a*fTmp)
      ELSE
        top_level = fu_set_pressure_level(height_a*fTmp)
        bottom_level = fu_set_pressure_level(height_b*fTmp)
      END IF
    ELSE
      call msg("Could not parse unit '"+chUnitVert+"'")
      CALL set_error('unit must be height or pressure', 'decode_time_parameter_v4')
      RETURN
    END IF

    !
    ! Set the single descriptor looking for it in the standard descriptor list
    !
    nDescriptorInParLine = 1  ! always for this type of parameter line
    allocate(ifSpecies(nDescriptorInParLine), stat=io_status)
    if(fu_fails(io_status==0,'Failes ifSpecies allocation','decode_time_parameter_v4'))return
    call set_cocktail_description(chCocktail%sp, cocktail_descr, ifSpecies(1))
    if(error)return
    !
    ! Find out whether this descriptor already exists in the source descriptor list. Add if not
    !
    ifNew = .true.
    do iDescr = 1, nDescrMax
      if(.not. defined(cocktail_descr_lst_of_src(iDescr)))exit
      if(cocktail_descr == cocktail_descr_lst_of_src(iDescr))then
        ifNew = .false.
        indDescriptorOfParLine(1) = iDescr
        exit
      endif
    end do
    !
    ! New descriptor: find the empty place at the end of list and add it there
    !
    if(ifNew)then
      do iDescr = 1, nDescrMax
        if(.not. defined(cocktail_descr_lst_of_src(iDescr)))then
          call copy_cocktail_descriptor(cocktail_descr, cocktail_descr_lst_of_src(iDescr))
          indDescriptorOfParLine(1) = iDescr
          exit
        endif
      end do
    endif
    if(error)return
    !
    ! Now deal with release rates
    !
    if(index(fu_str_u_case(chSrcRateUnit),'DESCRIPTOR_DEFAULT') == 1)then
      fSrcUnit2DescrUnit = 1.0
    else
      chDescrUnit = fu_basic_mass_unit(cocktail_descr) + '/sec'
      if(error)then
        call report(cocktail_descr)
        call set_error('Failed to get the default unit of descriptor','decode_time_parameter_v4')
        return
      endif
      fSrcUnit2DescrUnit = fu_conversion_factor(cocktail_descr, chSrcRateUnit, chDescrUnit)
      if(error .or. (fSrcUnit2DescrUnit .eps. real_missing))then
        call set_error('Cannot convert from:'+chSrcRateUnit+', to:'+chDescrUnit,'decode_time_parameter_v4')
        return
      endif
    endif

    !
    ! The order of the parameters depends on the line format
    !
    if(index(chContent,'NOW') > 0 .or. index(chContent,'LAST_METEO_TIME') > 0)then
      call set_release_parameters(params(1), & 
                                & fu_io_string_to_time(chTmp),&
                                & (/fRate * fSrcUnit2DescrUnit/), &
                                & xySize,bottom_level,top_level, &
                                & vertically_even_distribution, &
                                & zVelocity,&
                                & Tempr, &
                                & (/cocktail_descr/), &
                                & cocktail_descr_lst_of_src, &
                                & nDescrMax) !, &
!                                & indDescriptorOfParLine, & ! OUT: index of par_str descriptors in param descr array
!                                & nDescriptorInParLine)  ! OUT: number of such descriptors
      call enlarge_param_vector(params)  ! To create one more parameter
      call set_release_parameters(params(2), &
                                & fu_io_string_to_time(chTmp)+fu_set_interval_min(iDur),&
                                & (/fRate * fSrcUnit2DescrUnit/), &
                                & xySize,bottom_level,top_level, &
                                & vertically_even_distribution, &
                                & zVelocity,&
                                & Tempr, &
                                & (/cocktail_descr/), &
                                & cocktail_descr_lst_of_src, &
                                & nDescrMax) !, &
!                                & indDescriptorOfParLine, & ! OUT: index of par_str descriptors in param descr array
!                                & nDescriptorInParLine)  ! OUT: number of such descriptors
    else
      call set_release_parameters(params(iParam), &
                                & fu_set_time_UTC(iYr4,iMon,iDay,iHr,iMin,fSec),&
                                !local to UTC
                                & (/fRate * fSrcUnit2DescrUnit/), &
                                & xySize,bottom_level,top_level, &
                                & vertically_even_distribution, &
                                & zVelocity,&
                                & Tempr, &
                                & (/cocktail_descr/), &
                                & cocktail_descr_lst_of_src, &
                                & nDescrMax) !, &
!                                & indDescriptorOfParLine, & ! OUT: index of par_str descriptors in param descr array
!                                & nDescriptorInParLine)  ! OUT: number of such descriptors

    endif ! Type of the line - NOW or real time

    call free_work_array(chCocktail%sp)

  end subroutine decode_time_parameter_v4


  !*************************************************************************

  subroutine decode_time_parameter_v5(chContent, &  ! parameter line
                                    & iParam, &     ! index of the time slot to decode
                                    & chUnitVert, & ! unit of vertical axis
                                    & chSrcRateUnit, & ! the source rate unit
                                    & params, &     ! vector of time slots
                                    & cocktail_descr_lst_of_src, & ! list of source descriptors
                                    & iSrcType, &   ! point/area source
                                    & nDescrMax, &
                                    & indDescriptorOfParLine, & ! OUT: index of par_str descriptors in param descr array
                                    & nDescriptorInParLine, &  ! OUT: number of such descriptors
                                    & chDescrTotals, pArDescrTotals, & ! total emission: input string and output array
                                    & ifSpecies)               ! OUT: species or cocktail is emitted
    !
    ! Based on the name of the content decodes it and writes to the
    ! given p_arc/a_src.
    ! Note: grid must be defined prior to the first value is stored.
    ! For its actual creation we shall use grib_io module.
    ! Do not overlook that geographical parameters are in milliseconds,
    ! in order to keep them integer. 
    !
    ! This version can handle the multi-descriptor-multi-release strings.
    !
    implicit none

    ! Imported parameters
    character(len=*), intent(in) :: chContent, chUnitVert, chSrcRateUnit, chDescrTotals
    integer, intent(in) :: iParam, iSrcType, nDescrMax
    type(silam_release_parameters), dimension(:), pointer :: params
    TYPE(Tcocktail_descr), dimension(:), pointer :: cocktail_descr_lst_of_src
    integer, dimension(:), pointer :: indDescriptorOfParLine ! OUT: index of par_str descriptors in param descr array
    integer, intent(out) :: nDescriptorInParLine  ! OUT: number of such descriptors
    real, dimension(:), pointer :: pArDescrTotals
    logical, dimension(:), allocatable, intent(out) :: ifSpecies    ! this-far, not actually used...

    ! Local variables
    integer io_status, iYr4,iMon,iDay,iHr,iMin, iDescr, iTmp, nPars, iStart, iEnd
    real :: height_a, height_b, fSec, xySize, zVelocity, Tempr, fSrcUnit2DescrUnit, fTmp
    logical :: ifNew
    real, dimension(:), pointer :: fRates
    type(silja_interval) :: timeShiftUTC
    type(silja_time) :: timeTmp
    character(len=clen), dimension(:), allocatable :: chCocktail
    character(len=clen) :: chTmp, chTime, chDurationUnit, chDuration
    TYPE(Tcocktail_descr), dimension(:), allocatable :: cocktail_descr
    TYPE(silja_level) :: top_level, bottom_level
    character(len=unitNmLen) :: chDescrUnit
    CHARACTER(len=*), PARAMETER :: sub_name="decode_time_parameter_v5"

    do while(iParam > size(params)) 
      call enlarge_param_vector(params)
      if(error)return
    end do
    fRates => fu_work_array(max_species)

    !            yr mon day hr min sec diam bottom top vel tempr  descr1  rate1 descr2 rate2 ...
    !point src: 2003  4 15  1  10  0.    1.  925.  900. 10. 273.  PASSIVE  400   so2   100   ...

    !           yr mon day hr min sec bottom top  descr1  rate1 descr2 rate2  ...
    !area src: 2003  4 15  1  10  0.  925.  900. PASSIVE  400    so2    100   ...
    !----------------------------------------
    !
    ! Time variation of the main release parameters:
    !   <rate> - emission rate [mass]/[time], units must correspond to chRateUnit
    !   <bottom> <top>  - vertical boundaries,[hpa] or [m] correspond to chUnit_vertical
    !   <cocktail_name> - name of the nuclear cocktail or name of a cocktail file
    ! Idea is:
    ! 1. If there is only one line starting from word NOW then wall-clock time as 
    !    start time, then duration in minutes follows and then the set of parameters 
    !    They are all considered to be fixed in time
    ! 2. If there are several lines - then the first line determines UTC time start of 
    !    release and sets its initial parameters. 
    !    Each next line determines the parameters of the release at some other time.
    !    All characteristics are interpolated in time between the sequencial moments
    !    The last line determines the end of the release and final parameters
    !
    !
    ! First of all, we have to find out how many descriptors are in the string
    !
    do iTmp = 1, nDescrMax+1
      if(index(chContent,'NOW') > 0 .or. index(chContent,'LAST_METEO_TIME') > 0)then
        if(iSrcType == point_source)then
          read(unit=chContent,fmt=*,iostat=io_status)chTime,chDuration,chDurationUnit,xySize, &
                                                   & height_a,height_b,zVelocity, Tempr, &
                                                   & (chTmp, fRates(1), iDescr = 1, iTmp)
        elseif(iSrcType == area_source)then
          read(unit=chContent,fmt=*,iostat=io_status)chTime,chDuration,chDurationUnit,height_a,height_b, &
                                                   & (chTmp, fRates(1), iDescr = 1, iTmp)
        else
          call set_error('Unknown source type:' + fu_str(iSrcType),sub_name)
          return
        endif

      elseif(index(chContent,'VOID_TIME') > 0)then
        if(iSrcType == point_source)then
          read(unit=chContent,fmt=*,iostat=io_status)chTmp, xySize, height_a,height_b,zVelocity, Tempr, &
                                                   & (chTmp, fRates(1), iDescr = 1, iTmp)
        elseif(iSrcType == area_source)then
          read(unit=chContent,fmt=*,iostat=io_status)chTmp, height_a,height_b, &
                                                   & (chTmp, fRates(1), iDescr = 1, iTmp)
        else
          call set_error('Unknown source type:' + fu_str(iSrcType),sub_name)
          return
        endif

      else
        if(iSrcType == point_source)then
          read(unit=chContent,fmt=*,iostat=io_status) iYr4,iMon,iDay,iHr,iMin,fSec,&
                                                    & xySize,height_a,height_b, zVelocity, Tempr, &
                                                    & (chTmp, fRates(1), iDescr = 1, iTmp)
        elseif(iSrcType == area_source)then
          read(unit=chContent,fmt=*,iostat=io_status) iYr4,iMon,iDay,iHr,iMin,fSec, height_a,height_b, &
                                                    & (chTmp, fRates(1), iDescr = 1, iTmp)

        else
          call set_error('Unknown source type:' + fu_str(iSrcType),sub_name)
          return
        endif
      endif
      if(io_status /= 0)then
        nDescriptorInParLine = iTmp-1
        exit
      endif
    end do  ! up to nDescrMax
    if(fu_fails(nDescriptorInParLine >= 1, 'Strange number of descriptors in par line',sub_name))return
    
    !
    ! Prepare the space
    !
    allocate(chCocktail(nDescriptorInParLine), cocktail_descr(nDescriptorInParLine), &
           & ifSpecies(nDescriptorInParLine), stat=io_status)
    if(fu_fails(io_status == 0, 'Failed allocation for temporary',sub_name))return

    !
    ! And now can read the string
    !
    if(index(chContent,'NOW') == 0 .and. index(chContent,'LAST_METEO_TIME') == 0)then
      if(iSrcType == point_source)then
        read(unit=chContent,fmt=*,iostat=io_status) iYr4, iMon, iDay, iHr, iMin, fSec, &
                              & xySize,height_a, height_b, zVelocity, Tempr, &
                              & (chCocktail(iDescr), fRates(iDescr), iDescr=1, nDescriptorInParLine)
      elseif(iSrcType == area_source)then
        read(unit=chContent,fmt=*,iostat=io_status) iYr4, iMon, iDay, iHr, iMin, fSec,  &
                              & height_a, height_b, &
                              & (chCocktail(iDescr), fRates(iDescr), iDescr=1, nDescriptorInParLine)
        xySize = 0.
        zVelocity = 0.
        Tempr = 273.15
      else
        call msg('Unknown source type:',iSrcType)
        call set_error('Unknown source type', sub_name)
        return
      endif
    else
      if(iSrcType == point_source)then
        read(unit=chContent,fmt=*,iostat=io_status)chTmp, chDuration, chDurationUnit, xySize, &
                           & height_a, height_b, zVelocity, Tempr, &
                           & (chCocktail(iDescr), fRates(iDescr), iDescr=1, nDescriptorInParLine)
      elseif(iSrcType == area_source)then
        read(unit=chContent,fmt=*,iostat=io_status) chTmp, chDuration, chDurationUnit, &
                           & height_a, height_b, &
                           & (chCocktail(iDescr), fRates(iDescr), iDescr=1, nDescriptorInParLine)
        xySize = 0.
        zVelocity = 0.
        Tempr = 273.15
      else
        call msg('Unknown source type:',iSrcType)
        call set_error('Unknown source type', sub_name)
        return
      endif
    endif
    if(io_status /= 0)then
      call set_error(fu_connect_strings('Strange parameter string:',chContent), sub_name)
      return
    endif
    !
    ! Vertical boundaries for the just-read line
    ! 
    if (fu_unit_type(chUnitVert) == size_unit) then
      fTmp = fu_conversion_factor(chUnitVert, "m")
      IF (height_a > height_b) THEN
        top_level = fu_set_constant_height_level(height_a*fTmp)
        bottom_level = fu_set_constant_height_level(height_b*fTmp)
      ELSE
        top_level = fu_set_constant_height_level(height_b*fTmp)
        bottom_level = fu_set_constant_height_level(height_a*fTmp)
      END IF

    elseif (fu_unit_type(chUnitVert) == pressure_unit) then
      fTmp = fu_conversion_factor(chUnitVert, "Pa")
      IF (height_a > height_b) THEN
        top_level = fu_set_pressure_level(height_b*fTmp)
        bottom_level = fu_set_pressure_level(height_a*fTmp)
      ELSE
        top_level = fu_set_pressure_level(height_a*fTmp)
        bottom_level = fu_set_pressure_level(height_b*fTmp)
      END IF
    ELSE
      call msg("Could not parse unit '"+chUnitVert+"'")
      CALL set_error('unit must be height or pressure', sub_name)
      RETURN
    END IF
    !
    ! We have no cocktails so far because we may have different sources
    ! with different speciacion. Here we store the cocktail description for 
    ! all time slots.
    ! The chSubstNm is either the name of one of cocktail species or WHOLE_COCKTAIL.
    ! If it is a specific substance, its fraction is set to 1 and others are zeroed,
    ! otherwise the default fractions from the cocktail file are used
    !
    do iDescr = 1, nDescriptorInParLine
      !
      ! Create the descriptor (actually, select it from the database)
      !
      call set_cocktail_description(chCocktail(iDescr), cocktail_descr(iDescr), ifSpecies(iDescr))
      if(error)return
      !
      ! Find out whether this descriptor already exists in the source descriptor list. Add if not
      !
      ifNew = .true.
      do iTmp = 1, nDescrMax
        if(.not. defined(cocktail_descr_lst_of_src(iTmp)))exit
        if(cocktail_descr(iDescr) == cocktail_descr_lst_of_src(iTmp))then
          ifNew = .false.
          indDescriptorOfParLine(iDescr) = iTmp
          exit
        endif
      end do
      !
      ! New descriptor found: find the empty place at the end of the list and put it there.
      !
      if(ifNew)then
        do iTmp = 1, nDescrMax
          if(.not. defined(cocktail_descr_lst_of_src(iTmp)))then
            call copy_cocktail_descriptor(cocktail_descr(iDescr), cocktail_descr_lst_of_src(iTmp))
            indDescriptorOfParLine(iDescr) = iTmp
            exit
          endif
        end do
      endif
if(indDescriptorOfParLine(iDescr) == int_missing)then
  call set_error('Missing descriptor',sub_name)
endif
      if(error)return
      !
      ! Make up the conversion factor from the source global rate unit to the descriptor unit
      !
      if(index(fu_str_u_case(chSrcRateUnit),'DESCRIPTOR_DEFAULT') == 1)then
        fSrcUnit2DescrUnit = 1.0
      else
        chDescrUnit = fu_basic_mass_unit(cocktail_descr(iDescr)) + '/sec'
        if(error)then
          call report(cocktail_descr(iDescr))
          call set_error('Failed to get the default unit of descriptor',sub_name)
          return
        endif

        fSrcUnit2DescrUnit = fu_conversion_factor(cocktail_descr(iDescr), chSrcRateUnit, chDescrUnit)
        if(error .or. (fSrcUnit2DescrUnit .eps. real_missing))then
          call set_error('Cannot convert from:'+chSrcRateUnit+', to:'+chDescrUnit, sub_name)
          return
        endif
      endif
      fRates(iDescr) = fRates(iDescr) * fSrcUnit2DescrUnit
      !
      ! Find this descriptor in the line of emission totals and fill it in
      !
      iTmp = index(chDescrTotals,trim(chCocktail(iDescr)))
      if(iTmp == 0)then  ! may be, there is no such line?
        pArDescrTotals(indDescriptorOfParLine(iDescr)) = real_missing
      else
        read(unit=chDescrTotals(iTmp:),fmt=*)chTmp, pArDescrTotals(indDescriptorOfParLine(iDescr))
        pArDescrTotals(indDescriptorOfParLine(iDescr)) = fSrcUnit2DescrUnit * &
                                                   & pArDescrTotals(indDescriptorOfParLine(iDescr))
      endif

    end do  ! iDescr 


    !
    ! The order of the parameters depends on the line format
    !
    if(index(chContent,'NOW') > 0 .or. index(chContent,'LAST_METEO_TIME') > 0)then
      call enlarge_param_vector(params)  ! To create one more parameter
      timeTmp = fu_io_string_to_time(chTmp)
      iStart = 1
      iEnd = 2
    else
      timeTmp = fu_set_time_UTC(iYr4,iMon,iDay,iHr,iMin,fSec)
      iStart = iParam
      iEnd = iParam
    endif ! Type of the line - NOW or real time
    do iTmp = iStart, iEnd
      call set_release_parameters(params(iTmp), & 
                                & timeTmp,&
                                & fRates, &
                                & xySize,bottom_level,top_level, &
                                & vertically_even_distribution, &
                                & zVelocity,&
                                & Tempr, &
                                & cocktail_descr, &
                                & cocktail_descr_lst_of_src, &
                                & nDescrMax)
      if ( iStart == iEnd) exit          
      timeTmp = timeTmp + fu_set_named_interval(chDuration+chDurationUnit)
    end do

    deallocate(chCocktail,cocktail_descr)
    call free_work_array(fRates)

  end subroutine decode_time_parameter_v5

                                    
  !*************************************************************************

  subroutine decode_time_params_volc_v1(chContent, chVertUnit, params, &   ! item to decode and structure to put it into
                                      & cocktail_descr_lst, nDescriptors, &    ! list of source cocktail descriptors
                                      & ifSpecies)
    !
    ! Decondes one line of volcano time parameters
    !
    implicit none
    
    ! Importe parameters
    character(len=*), intent(in) :: chContent, chVertUnit
    type(silam_release_parameters), dimension(:), intent(inout) :: params
    TYPE(Tcocktail_descr), dimension(:), intent(inout) :: cocktail_descr_lst
    logical, dimension(:), allocatable, intent(out) :: ifSpecies
    integer, intent(out) :: nDescriptors
    
    ! Local variables
    integer :: io_status, nItems, iParam, iDescr, iTmp
    character(len=clen), dimension(100) :: arStr

    ! Split the content, get the items
    ! format: 1 2010 4 14 0 0 0.  COCKTAIL_SO2 8000.   COCKTAIL_ASH 5000.
    ! Note that the vertical unit is given separately
    !
    call split_string(chContent, ' ', arStr, nItems)
    if(error)return
    !
    ! Reference sizes: params index and number of descriptors
    !
    nDescriptors = (nItems-7) / 2
    read(unit=arStr(1), fmt=*)iParam
    params(iParam)%time = fu_io_string_to_time(chContent(index(chContent,' ') : ))
    if(error)return
    !
    ! If rate_descr_inut is allocated. we already treated this time step: error. 
    !
    if(allocated(params(iParam)%rate_descr_unit))then
      call set_error('Re-allocation of already allocated params element. Dupilicated time index?','decode_time_params_volc_v1')
      return
    else
      allocate(params(iParam)%rate_descr_unit(nDescriptors), ifSpecies(nDescriptors), stat=io_status)
      if(fu_fails(io_status==0,'Failed params rate allocation','decode_time_params_volc_v1'))return
    endif
    !
    ! Knowing the number of items, we can get the number of cocktail descriptors: (nItems-7/)2
    !
    do iDescr = 1, nDescriptors
      iTmp = (iDescr-1)*2 + 8   ! start of the new <cocktail> <height> tupla
      call set_cocktail_description(arStr(iTmp), cocktail_descr_lst(iDescr), ifSpecies(iDescr))
      if(error)return
      ! rate_descr_unit is a placeholder - I need an array of the proper dimension to store the injection height.
      params(iParam)%rate_descr_unit(iDescr) = fu_set_named_value(arStr(iTmp+1) // ' ' // chVertUnit)
    end do

  end subroutine decode_time_params_volc_v1


  !*************************************************************************

  subroutine enlarge_param_vector(params, nDescriptors)
    !
    ! Just increases the size of the param vector for 1. The reason for so
    ! small increase - in a several places the size(params) is used without 
    ! checking of reasonable parameter value
    !
    implicit none

    ! Imported parameters
    type(silam_release_parameters), dimension(:), pointer :: params
    integer, intent(in), optional :: nDescriptors

    ! Local variables
    integer :: iTmp, jTmp
    type(silam_release_parameters), dimension(:), pointer :: parTmp

    if(associated(params))then
      allocate(parTmp(size(params)),stat=jTmp)
      if(fu_fails(jTmp == 0,'Failed to allocate memory for params.1','enlarge_param_vector'))return
      do iTmp = 1, size(params)
        allocate(parTmp(iTmp)%rate_descr_unit(params(1)%nDescriptors), stat=jTmp)
        if(fu_fails(jTmp == 0, 'Failed to allocate memory for params.2','enlarge_param_vector'))return
      end do
      do iTmp=1,size(params)
        parTmp(iTmp)%time = params(iTmp)%time
        parTmp(iTmp)%xy_size = params(iTmp)%xy_size
        parTmp(iTmp)%z_velocity = params(iTmp)%z_velocity
        parTmp(iTmp)%tempr= params(iTmp)%tempr
        parTmp(iTmp)%layerDynamic = params(iTmp)%layerDynamic
        parTmp(iTmp)%vertical_distribution = params(iTmp)%vertical_distribution
        parTmp(iTmp)%nDescriptors = params(1)%nDescriptors
        do jTmp = 1, params(1)%nDescriptors
          parTmp(iTmp)%rate_descr_unit(jTmp) = params(iTmp)%rate_descr_unit(jTmp)
        end do
        deallocate(params(iTmp)%rate_descr_unit)
      end do
      deallocate(params); 
      allocate(params(size(parTmp)+1),stat=jTmp)
      if(fu_fails(jTmp == 0, 'Failed to allocate memory for values.3','enlarge'))return
      do iTmp = 1, size(params)
        allocate(params(iTmp)%rate_descr_unit(params(1)%nDescriptors), stat=jTmp)
        if(fu_fails(jTmp == 0, 'Failed to allocate memory for params.2','enlarge_param_vector'))return
      end do
      do iTmp=1,size(parTmp)
        params(iTmp)%time = parTmp(iTmp)%time
        params(iTmp)%xy_size = parTmp(iTmp)%xy_size
        params(iTmp)%z_velocity = parTmp(iTmp)%z_velocity
        params(iTmp)%tempr= parTmp(iTmp)%tempr
        params(iTmp)%layerDynamic = parTmp(iTmp)%layerDynamic
        params(iTmp)%vertical_distribution = parTmp(iTmp)%vertical_distribution
        params(iTmp)%nDescriptors = parTmp(1)%nDescriptors
        do jTmp = 1, params(1)%nDescriptors
          params(iTmp)%rate_descr_unit(jTmp) = parTmp(iTmp)%rate_descr_unit(jTmp)
        end do
      end do
      params(size(parTmp)+1)%nDescriptors = params(1)%nDescriptors

      deallocate(parTmp)

    else
      !
      ! Entirely new parameter vector - allocate it. But for that we need nDescriptors
      !
      if(.not. present(nDescriptors))then
        call set_error('Need nDescriptors for new param vector creation','enlarge_param_vector')
        nullify(params)
        return
      endif
      allocate(params(2), stat=jTmp)
      if(fu_fails(jTmp == 0, 'Failed to allocate memory for values.3','enlarge'))return
      do iTmp=1, 2
        allocate(params(iTmp)%rate_descr_unit(nDescriptors), stat=jTmp)
        if(fu_fails(jTmp == 0, 'Failed to allocate memory for params.2','enlarge_param_vector'))return
        params(iTmp)%nDescriptors = nDescriptors
      end do
    endif

  end subroutine enlarge_param_vector


  !********************************************************************************************
  
  subroutine integrate_dynamic_slot_descr(params, timeStart, timeEnd, iSlot, &
                                        & cocktail_descr_lst, nDescr, &
                                        & fDescrScaling, &
                                        & fVertFraction1, fVertFraction2, &
                                        & ifRateVary, &
                                        & amounts)
    !
    ! Integrates dynamic slot. Depending on the number of varying dimensions, have to use 
    ! either descriptor-based mapping or can integrate over the descriptors and use species-wise
    ! Descriptor amount start/end
    !
    implicit none

    ! Imported parameters
    type(silam_release_parameters), dimension(:), pointer :: params
    type(silja_time), intent(in) :: timeStart, timeEnd
    integer, intent(in) :: iSlot, nDescr
    TYPE(Tcocktail_descr), dimension(:), pointer :: cocktail_descr_lst
    real, intent(in) :: fVertFraction1, fVertFraction2
    logical, intent(in) :: ifRateVary
    real, dimension(:), intent(in) :: fDescrScaling    ! 
    real, dimension(:), intent(out) :: amounts          ! output

    ! Local variables
    integer :: iDescr
    real :: fSlot_sec, fTime_sec_start, fTime_sec_end, fA, fB, fC

    fSlot_sec = fu_sec(params(iSlot+1)%time - params(iSlot)%time)
    fTime_sec_start = fu_sec(timeStart - params(iSlot)%time)
    fTime_sec_end = fu_sec(timeEnd - params(iSlot)%time)
    if(error)return

    !
    ! Use every possible chance to reduce the amount of computations but start from brute-force
    !
    if(ifRateVary)then
      !
      ! Either rates and, may be, vertical. Brute force
      !
      do iDescr = 1, nDescr
        !
        ! rates, may be, vertical
        !
        fA = params(iSlot)%rate_descr_unit(iDescr) * &  ! descr.rate*composition, spiecies unit
           & fVertFraction1                ! vertical overlap

        fB = (params(iSlot)%rate_descr_unit(iDescr) * &
                & (fVertFraction2 - fVertFraction1) + &
            & (params(iSlot+1)%rate_descr_unit(iDescr) - params(iSlot)%rate_descr_unit(iDescr)) * &
                & fVertFraction1) / fSlot_sec

        fC = ((params(iSlot+1)%rate_descr_unit(iDescr) - params(iSlot)%rate_descr_unit(iDescr)) * &
                & (fVertFraction2 - fVertFraction1)) / (fSlot_sec*fSlot_sec)

        amounts(iDescr) = amounts(iDescr) + fDescrScaling(iDescr) * &
                         & (fA * (fTime_sec_end - fTime_sec_start) + &
                          & fB * (fTime_sec_end*fTime_sec_end - fTime_sec_start*fTime_sec_start) / 2. + &
                          & fC * (fTime_sec_end*fTime_sec_end*fTime_sec_end - &
                                      & fTime_sec_start*fTime_sec_start*fTime_sec_start) / 3.)
       end do

    else                                                     
      !
      ! Only vertical varies
      !
      do iDescr = 1, nDescr
        fA = params(iSlot)%rate_descr_unit(iDescr) * fVertFraction1  ! descr.rate*composition * vertical overlap
        fB = (params(iSlot)%rate_descr_unit(iDescr) * (fVertFraction2 - fVertFraction1)) / fSlot_sec
        fC = 0.
        amounts(iDescr) = amounts(iDescr) + fDescrScaling(iDescr) * &
                         & (fA * (fTime_sec_end - fTime_sec_start) + &
                          & fB * (fTime_sec_end*fTime_sec_end - fTime_sec_start*fTime_sec_start) / 2.)
      end do

    endif   ! what varies while only rate is requested

  end subroutine integrate_dynamic_slot_descr


  !***********************************************************************************

  subroutine store_time_param_as_nl_params(param, cocktail_descr_of_src, uOut, iSourceType)
    !
    ! Stores the parameter into the namelist line, so that it could be read later on
    ! by the above subroutines. Note that we always write the v3 for area source and v5 for point
    !
    implicit none

    ! Imported parameters
    type(silam_release_parameters), intent(in) :: param
    TYPE(Tcocktail_descr), dimension(:), pointer :: cocktail_descr_of_src
    integer, intent(in) :: uOut, iSourceType

    ! local variables
    integer :: iDescr

    if(iSourceType == point_source)then
      write(uOut, fmt='(A,1x,i4,4I3,1x,F4.1,1x,5(F15.7,1x),50(A,1x,F15.7,1x))') 'par_str_point = ', &
!      write(uOut, fmt='(A,i4,1x,4I3,1x,6(F15.7,1x),50(A,1x,F15.7))') 'par_str_point = ', &
                     & fu_year(param%time), &
                     & fu_mon(param%time), &
                     & fu_day(param%time), &
                     & fu_hour(param%time), &
                     & fu_min(param%time), &
                     & fu_sec(param%time), &
                     & param%xy_size, &
                     & fu_top_of_layer_value(param%layerDynamic), &
                     & fu_bottom_of_layer_value(param%layerDynamic), &
                     & param%z_velocity, &
                     & param%tempr, &
                     & (trim(fu_name(cocktail_descr_of_src(iDescr))), &
                      & param%rate_descr_unit(iDescr), &
                      & iDescr = 1, param%nDescriptors)
    elseif(iSourceType == area_source)then
      write(uOut, fmt='(A,1x,i4,4I3,1x,F4.1,1x,2(F15.7,1x),50(A,1x,F15.7,1x))') 'par_str_area = ', &
!      write(uOut, fmt='(A,i4,1x,4I3,1x,3(F15.7,1x),50(A,1x,F15.7,1x))') 'par_str_area = ', &
                     & fu_year(param%time), &
                     & fu_mon(param%time), &
                     & fu_day(param%time), &
                     & fu_hour(param%time), &
                     & fu_min(param%time), &
                     & fu_sec(param%time), &
                     & fu_top_of_layer_value(param%layerDynamic), &
                     & fu_bottom_of_layer_value(param%layerDynamic), &
                     & (trim(fu_name(cocktail_descr_of_src(iDescr))), param%rate_descr_unit(iDescr), &
                                                                    & iDescr = 1, param%nDescriptors)
    else
      call msg('Unknown type of the source:',iSourceType)
      call set_error('Unknown type of the source','store_time_param_as_nl_params')
    endif

  end subroutine store_time_param_as_nl_params


  !***********************************************************************************

  subroutine store_time_param_as_nl_time(timePar, fTop, fBottom, species, nSpecies, &
                                       & uOut, iSourceType)
    !
    ! Stores the parameter into the namelist line, so that it could be read later on
    ! by the above subroutines. Note that we always write the v3 for area source and v5 for point
    !
    implicit none

    ! Imported parameters
    type(silja_time), intent(in) :: timePar
    TYPE(silam_species), dimension(:), intent(in) :: species
    real, intent(in) :: fTop, fBottom
    integer, intent(in) :: nSpecies, uOut, iSourceType

    ! local variables
    integer :: iDescr
    type(silam_sp) :: sp

    sp%sp => fu_work_string()
    if(error)return

    if(iSourceType == point_source)then
      call set_error('Does not work yet for point sources','store_time_param_as_nl_time')
      call free_work_array(sp%sp)
      return
    elseif(iSourceType == area_source)then
!      write(unit=sp%sp, fmt='(A,1x,i4,4I3,1x,F4.1,1x,2(F15.7,1x),50(A,1x,2A,1x,E9.3,1x))') &
      write(unit=sp%sp, fmt='(A,1x,i4,4I3,1x,F4.1,1x,2(E9.3,1x),50(A,1x,E9.3,1x))') &
                                   & 'par_str_area = ', &
                                   & fu_year(timePar), &
                                   & fu_mon(timePar), &
                                   & fu_day(timePar), &
                                   & fu_hour(timePar), &
                                   & fu_min(timePar), &
                                   & fu_sec(timePar), &
                                   & fTop, &  !top 
                                   & fBottom , & !bottom
                                   & (trim(fu_name(fu_material(species(iDescr)))), &
!                                    & trim(fu_basic_mass_unit(fu_material(species(iDescr)))),'/sec', &
                                    & 1.0, &
                                    & iDescr = 1, nSpecies)
    endif
    write(uOut,fmt='(A)')trim(sp%sp)

    call free_work_array(sp%sp)

  end subroutine store_time_param_as_nl_time


  !*****************************************************************************************

  subroutine count_descriptor_names(line, nDescriptors, chDescriptorNames)
    !
    ! Gets the par_str* line and counts the descriptors in it. The newly found descriptors are
    ! embedded into the existing list and the grand total is returned
    !
    implicit none

    ! Imported parameters
    character(len=*), intent(in) :: line
    integer, intent(inout) :: nDescriptors
    character(len=*), dimension(:), intent(inout) :: chDescriptorNames

    ! Local variables
    character(len=clen) :: chTmp
    character(len=clen), dimension(max_descriptors) :: chArTmp
    integer :: iDescr, iTmp, i, nDescrInArTmp, io_status
    real :: f

    !
    ! There are four types of the lines: par_str for point or area source, par_str_point, and par_str_area
    ! Their treatment is different
    !
    if(index(adjustl(line),'par_str_point') == 1)then
      !
      ! Try to read the line with as many descriptor names as possible
      !
      do iTmp = 1, 1000
        if(index(line,'NOW') > 0 .or. index(line,'LAST_METEO_TIME') > 0)then
          read(unit=line,fmt=*,iostat=io_status) chTmp, chTmp, chTmp,f,chTmp,f, f,f,f,f, &
                                                   & (chArTmp(iDescr), f, iDescr = 1, iTmp)
        elseif(index(line,'VOID_TIME') > 0)then
          read(unit=line,fmt=*,iostat=io_status) chTmp, chTmp, chTmp, f, f,f,f, f, &
                                                   & (chArTmp(iDescr), f, iDescr = 1, iTmp)
        else
          read(unit=line,fmt=*,iostat=io_status) chTmp, chTmp, i,i,i,i,i,f,f,f,f, f, f, &
                                                    & (chArTmp(iDescr), f, iDescr = 1, iTmp)
        endif
        if(io_status /= 0)then
          nDescrInArTmp = iTmp-1
          exit
        endif
      end do  ! up to nDescrMax
      if(nDescrInArTmp< 1)then
        call set_error('Strange number of descriptors in par_str_point line','count_descriptor_names')
        return
      endif

    elseif(index(adjustl(line),'par_str_area') == 1)then
      !
      ! Try to read the line with as many descriptor names as possible
      !
      do iTmp = 1, 1000
        if(index(line,'NOW') > 0 .or. index(line,'LAST_METEO_TIME') > 0)then
          read(unit=line,fmt=*,iostat=io_status) chTmp, chTmp, chTmp,f,chTmp,f,f, &
                                               & (chArTmp(iDescr),f, iDescr=1,iTmp)

        elseif(index(line,'VOID_TIME') > 0)then
          read(unit=line,fmt=*,iostat=io_status) chTmp, chTmp, chTmp, f,f, &
                                               & (chArTmp(iDescr), f, iDescr=1,iTmp)

        else
          read(unit=line,fmt=*,iostat=io_status) chTmp, chTmp, i,i,i,i,i,f,f,f, &
                                               & (chArTmp(iDescr), f, iDescr=1,iTmp)
        endif
        if(io_status /= 0)then
          nDescrInArTmp = iTmp-1
          exit
        endif
      end do  ! up to nDescrMax
      if(nDescrInArTmp< 1)then
        call set_error('Strange number of descriptors in par_str_area line','count_descriptor_names')
        return
      endif
      
    elseif(index(adjustl(line),'par_str_volcano') == 1)then
      !
      ! Split the line and count the number of the corresponding fields
      ! Line is of this kind 
      ! par_str_volcano = 1 2017 7 22 0 0 0.  COCKTAIL_SO2 8. COCKTAIL_ASH 5.
      !
      call split_string(line, ' ', chArTmp, nDescrInArTmp)
      if(error)return
      nDescrInArTmp = (nDescrInArTmp - 9) / 2
      if(nDescrInArTmp < 1)then
        call set_error('Strange number of descriptors in par_str_volcano line:' + line,'count_descriptor_names')
        return
      endif
      do iDescr = 1, nDescrInArTmp
        chArTmp(iDescr) = chArTmp((iDescr-1)*2 + 9 + 1)
      enddo
      
    elseif(index(adjustl(line),'par_str ') == 1)then
      !
      ! Only one descriptor at the end of the line. Check if it is new and add to the list if it is.
      !
      chArTmp(1) = line(index(trim(line),' ',.true.)+1:)
!call msg('Descriptor in the line:' + chArTmp(1) + '<=end')
      nDescrInArTmp = 1
    else
      call set_error('No par_str in the line:'+line,'count_descriptor_names')
      return
    endif

    !
    ! Now copy the temporary array to the main one checking for duplicates
    !
    if(nDescriptors == 0)then
      !
      ! New source, first par_str line, just copy
      !
      if(size(chDescriptorNames) == nDescrInArTmp)then
        call set_error('Too small descriptor names array','count_descriptor_names')
        return
      endif
      do iTmp = 1, nDescrInArTmp
        chDescriptorNames(iTmp) = chArTmp(iTmp)
      end do
      nDescriptors = nDescrInArTmp
    else
      !
      ! Check that new and existing stuff is identical 
      !
      do iTmp = 1, nDescrInArTmp
        if(chArTmp(iTmp) /= chDescriptorNames(iTmp))then
          call msg('Descriptor in line:'+chArTmp(iTmp)+',existing one:'+chDescriptorNames(iTmp)+ &
                 & '- at position:',iTmp)
          call set_error('New line has different descriptors:'+line,'count_descriptor_names')
          return
        endif
      end do  ! over the temporary descriptor names array
    endif

  end subroutine count_descriptor_names


  !******************************************************************************************

  subroutine add_descriptor_names(nExistingDescriptors, chExistingDescriptorNames, &
                                & nNewDescriptors, chNewDescriptorNames, &
                                & chNewSubstance, iNewMode)
    !
    ! Adds the descriptor names collected from diffeent files for the single source.
    ! Checks that the new descriptors are not the same as the existing ones
    !
    implicit none

    ! Imported parameters
    integer, intent(inout) :: nExistingDescriptors
    character(len=*), dimension(:), intent(inout) :: chExistingDescriptorNames
    integer, intent(in) :: nNewDescriptors
    character(len=*), dimension(:), intent(in) :: chNewDescriptorNames
    character(len=*), intent(in) :: chNewSubstance
    integer, intent(in) :: iNewMode

    ! Local stuff
    integer :: iTmp, iTmp2

    
    if(nExistingDescriptors < 1)then
      call msg_warning('No existing descriptors in existing source','add_descriptor_names')
      if(size(chNewDescriptorNames) < nNewDescriptors)then
        call set_error('Too small descriptor names array','count_descriptor_names')
        return
      endif
      do iTmp = 1, nNewDescriptors
        chExistingDescriptorNames(iTmp) = fu_descr_name(chNewDescriptorNames(iTmp), &
                                                      & chNewSubstance, &
                                                      & iNewMode)
      end do
      nExistingDescriptors = nNewDescriptors
    else

      do iTmp = 1, nNewDescriptors
        do iTmp2 = 1, nExistingDescriptors   ! Check that new is not already there 
          if(chExistingDescriptorNames(iTmp2) == fu_descr_name(chNewDescriptorNames(iTmp), &
                                                             & chNewSubstance, &
                                                             & iNewMode))then
            call msg("New:"+fu_descr_name(chNewDescriptorNames(iTmp), chNewSubstance, iNewMode))
            call msg("Existing:"+chExistingDescriptorNames(iTmp2))

            call set_error('source emits same descriptor twice or duplicated source name-sector pair', &
                         & 'add_descriptor_names')
            return
          endif
        enddo
        nExistingDescriptors = nExistingDescriptors +1
        if(size(chNewDescriptorNames) < nExistingDescriptors)then
          call set_error('Too small descriptor names array','count_descriptor_names')
          return
        endif
        chExistingDescriptorNames(nExistingDescriptors) = fu_descr_name(chNewDescriptorNames(iTmp), &
                                                                      & chNewSubstance, &
                                                                      & iNewMode)
      end do
    endif

  end subroutine add_descriptor_names


  !**************************************************************************************

  subroutine complete_descriptor_names(chDescrNames, chSubstNm, iModeNbr, nDescriptors)
    !
    ! Backward compatibility only !
    ! This procedure ensures generation of something resembling the descriptor name
    ! using the v.4-style area and point sources with separate "big" cocktail descriptors
    ! in the par_str and additional substance name and mode number as separate lines
    ! For v.5, the descriptor will be compiled as a combination of these three.
    !
    implicit none

    ! Imported parameters
    integer, intent(inout) :: nDescriptors
    character(len=*), dimension(:), intent(inout) :: chDescrNames
    character(len=*), intent(in) :: chSubstNm
    integer, intent(in) :: iModeNbr

    ! Local variables
    integer :: iDescr

    do iDescr = 1, nDescriptors
      chDescrNames(iDescr) = fu_descr_name(chDescrNames(iDescr), chSubstNm, iModeNbr)
    end do

  end subroutine complete_descriptor_names


  !**************************************************************************************

  function fu_descr_name(chDescrIni, chSubstIni, iModeIni) result(name)
    !
    ! Generates the v.5 descriptor name out of given descriptor, substance and mode
    !
    implicit none

    ! result
    character(len=clen) :: name

    ! Imported parameters
    character(len=*), intent(in) :: chDescrIni, chSubstIni
    integer, intent(in) :: iModeIni

    if(len_trim(chSubstIni) == 0)then
      name = chDescrIni
    elseif(fu_str_u_case(chSubstIni) == 'WHOLE_COCKTAIL')then
      name = chDescrIni
    else
      if(iModeIni < 1)then
        name = chDescrIni + '_' + chSubstIni
      else
        name = chDescrIni + '_' + chSubstIni + '_' + fu_str(iModeIni)
      endif
    endif

  end function fu_descr_name

  
  !***************************************************************************************
  
  subroutine get_source_aer_species(speciesSrc, nSpeciesSrc, expected_species, chSubstNm, nlSetup)
    !
    ! The list of aerosol modes that will be emitted. To set them up, we will
    ! use the standard aerosol procedure - for the sake of unification.
    ! Note that there are two potentially concurring definitions: one coming from aerosol dynamics
    ! the other - written in the ini file. The first one prevails, if it exists
    !
    implicit none
    
    ! Imported parameters
    type(silam_species), dimension(:), pointer :: speciesSrc
    type(silam_species), dimension(:), intent(in), allocatable :: expected_species
    integer, intent(out) :: nSpeciesSrc
    character(len=*), intent(in) :: chSubstNm
    type(Tsilam_namelist),  intent(in) :: nlSetup
    
    ! Local variables
    integer, dimension(max_species) :: indices
    integer :: iTmp, jTmp, iSpecies
    logical :: ifFound
    type(silam_species) :: speciesTmp
    type(Taerosol) :: aerosolSrc
    
    nullify(speciesSrc)
    nSpeciesSrc = 0

    if(allocated(expected_species))then
      !
      ! Species have been ordered by the aerosol dynamics. Use that array
      !
      !
      ! Aerosol dynamics can (and usually does) request e.g. sea salt in specific modes
      !
      call select_species(expected_species, size(expected_species), chSubstNm, aerosol_mode_missing, &
                        & real_missing, indices, iTmp)
      if(error)return

      do iSpecies = 1, iTmp
        call addSpecies(speciesSrc, nSpeciesSrc, (/expected_species(indices(iSpecies))/), 1)
      end do
      !
      ! If someone is interested in number emission directly from this source, it has to be said
      ! via adding the special material 'nbr_aer', which must, of course, have the proper modes
      !
      call select_species(expected_species, size(expected_species), 'nbr_aer', aerosol_mode_missing, &
                        & real_missing, indices, iTmp)
      if(error)return
      do iSpecies = 1, iTmp
        !
        ! The nbr_aer mode must be one of the modes requested for the mass emission
        !
        ifFound = .false.
        do jTmp = 1, nSpeciesSrc
          ifFound = (fu_mode(speciesSrc(jTmp)) == fu_mode(expected_species(indices(iSpecies))))
          if(ifFound)exit
        end do
        if(ifFound)then
          !
          ! Mode found, can add this number-emission species to the list
          !
          call addSpecies(speciesSrc, nSpeciesSrc, (/expected_species(indices(iSpecies))/), 1)
        else
          call msg_warning('nbr_aer species have strange mode:','get_source_aer_species')
          call report(expected_species(indices(iSpecies)))
          call msg('Available modes are:')
          do jTmp = 1, nSpeciesSrc
            call report(speciesSrc(jTmp))
          end do
          call set_error('nbr_aer species have strange mode:','get_source_aer_species')
          return
        endif
      end do

    else
      !
      ! Species are up to the source itself
      !
      call set_aerosol(nlSetup, aerosolSrc)
      if(error)return
      !
      ! With the modes stored, we can create the list of species for the source
      !
      do iTmp = 1, fu_nModes(aerosolSrc)
        call set_species(speciesTmp, fu_get_material_ptr(chSubstNm), fu_mode(aerosolSrc, iTmp))
        if(error)return
        call addSpecies(speciesSrc, nSpeciesSrc, (/speciesTmp/), 1)
        if(error)return
      end do

    endif  ! whether the species list is up to source or requested by aerosol dyanmics
    
  end subroutine get_source_aer_species

  ! ***************************************************************

  subroutine getTimeSlots_from_params(params, time, iSlot1, iSlot2, fWeight_past)
    !
    ! Finds the time slots surrounding the given time
    ! 
    IMPLICIT NONE

    ! Imported parameters
    type(silam_release_parameters), dimension(:), intent(in) :: params
    type(silja_time), intent(in) :: time
    integer, intent(out) :: iSlot1, iSlot2
    real, optional, intent(out) :: fWeight_past
    !
    ! Stupidity check
    !
    if(present(fWeight_past)) fWeight_past = real_missing
    if(time < params(1)%time .or. time > params(size(params))%time)then
      call set_error('Given time is outside the source activity period', 'getTimeSlots_from_params')
      iSlot1=-1
      iSlot2=-1
      return
    endif
    !
    ! Search
    !
    do iSlot2 = 1, size(params)
      if(time < params(iSlot2)%time) then
        iSlot1 = iSlot2 -1
        if(present(fWeight_past)) &
            & fWeight_past = (params(iSlot2)%time - time) / (params(iSlot2)%time - params(iSlot1)%time)
        return
      endif
    end do
    !
    ! Last time moment ?
    !
    if (time == params(size(params))%time) then
      iSlot1 = size(params)
      iSlot2 = iSlot1
      if(present(fWeight_past)) fWeight_past = 1.0
      return 
    end if
    !
    ! One should never come to here
    !
    call set_error('Failed to find the proper time slots','getTimeSlots_from_params')
    iSlot1=-1
    iSlot2=-1
  END subroutine  getTimeSlots_from_params

  
  !**********************************************************************
  
  subroutine get_overlapping_slots(params, start, duration, iStartSlot, iEndSlot)
    implicit none
    type(silam_release_parameters), dimension(:), intent(in) :: params
    type(silja_time), intent(in) :: start
    type(silja_interval), intent(in) :: duration
    integer, intent(out) :: iStartSlot, iEndSlot
    
    integer :: iSlt1
    
    iStartSlot = int_missing
    iEndSlot = int_missing
    call getTimeSlots_from_params(params, start, iStartSlot, iSlt1)
    if(error)return
    call getTimeSlots_from_params(params, start+duration, iSlt1, iEndSlot)
    if(error)return
    
  end subroutine get_overlapping_slots

  !**************************************************************************************
  !
  ! Weather dependence of emission
  !
  !**************************************************************************************

  !*********************************************************************
  
  
  subroutine get_temp_dep_from_namelist(rulesMeteoDep, nlSrc)
    !
    ! Reads the parameters of the temperature dependence of emission for the given source term
    !
    implicit none
  
    type(TmeteoDepRules), intent(out) :: rulesMeteoDep
    type(Tsilam_namelist), intent(in) :: nlSrc
    CHARACTER(len=*), PARAMETER :: sub_name='get_temp_dep_from_namelist'

    if(fu_str_u_case(fu_content(nlSrc,'tempr_dependence_model')) == 'TEMPR_DEFICIT_REGRESSION')then
      !
      ! Need standard temperature indoor, fraction of non-heating combustion and insulation quality
      !
      rulesMeteoDep%indMeteoDepModelSwitch = emis_heating_vs_T_linregr_mdl
      rulesMeteoDep%fT0 = fu_content_real(nlSrc,'T0_K') 
      ! slope, default. Dynamic part of dependence.
      rulesMeteoDep%fdEdT = fu_content_real(nlSrc,'dEdT') 
      ! intercept, default. Temperature-independent fraction
      rulesMeteoDep%fMinEmission = fu_content_real(nlSrc,'frac_const')
      !
      ! All three should present
      !
      if(any ((/ rulesMeteoDep%fT0, rulesMeteoDep%fdEdT,  rulesMeteoDep%fMinEmission/) & 
                &  == real_missing)) then
        call msg_warning('Failed setup of temperature dependence',sub_name)
        call msg('From namelist: T0_K, dEdT, frac_const' , &
               & (/rulesMeteoDep%fT0, rulesMeteoDep%fdEdT, rulesMeteoDep%fMinEmission/))
        call set_error('Failed setup of temperature dependence',sub_name)
      endif
    else
      ! Unknown model
      call set_error('Unknown Tempr dependence model:' + fu_content(nlSrc,'tempr_dependence_model'), sub_name)
      call msg('Only TEMPR_DEFICIT_REGRESSION is allowed')
    endif  !  Model of dependence
      
  end subroutine get_temp_dep_from_namelist
  
  
  ! *******************************************************************
  
  subroutine add_input_needs_meteo_dependence(rulesMeteoDep, q_met_dynamic, q_met_st, &
                                                           & q_disp_dynamic, q_disp_st)
    !
    ! Adds the needed input data, dynamic and static meteo fields
    !
    implicit none
    
    type(TmeteoDepRules), intent(in) :: rulesMeteoDep
    integer, dimension(:), intent(inout) :: q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st
    
    ! Local variables
    integer :: iTmp

    select case(rulesMeteoDep%indMeteoDepModelSwitch)
      case(emis_heating_vs_T_linregr_mdl)
        iTmp = fu_merge_integer_to_array(day_mean_temperature_2m_flag, q_disp_st)
      case(int_missing)
      case default
        call set_error('Unknown meteo-depence model:' + fu_str(rulesMeteoDep%indMeteoDepModelSwitch), &
                     & 'add_input_needs_meteo_dependence')
    end select
  
  end subroutine add_input_needs_meteo_dependence
  

  !*********************************************************************************************

  integer function fu_meteodep_model_time_params(rules)
    implicit none
    type(TmeteoDepRules), intent(in) :: rules
    fu_meteodep_model_time_params = rules%indMeteoDepModelSwitch
  end function fu_meteodep_model_time_params

  
  !*********************************************************************************************
  
  subroutine prepare_injection_meteodep(dispBuf, arMeteoDepModels)
    !
    ! Stores the indices of the meteodata field(s). Called for each source type that has such
    ! dependence. So, multiple calls are probable
    !
    implicit none
    
    ! Input parameters
    type(Tfield_buffer), intent(in) :: dispBuf
    integer, dimension(:), intent(in) :: arMeteoDepModels
    
    ! Local variables
    integer, dimension(:), pointer :: met_q
    integer :: iQ, iMdl, quant

    ! nullify the pointer
    nullify(fldDailyTempr2m)

    ! Scan the meteo buffer and the dependence models
    ! Each model may need own meteo data
    !
    met_q => dispBuf%buffer_quantities    
    do iQ = 1, size(met_q)
      quant = met_q(iQ)
      if(quant == int_missing)exit  ! no more quantities to check
      do iMdl = 1, size(arMeteoDepModels)   
        select case(arMeteoDepModels(iMdl))
          case(emis_heating_vs_T_linregr_mdl)
            if(quant == day_mean_temperature_2m_flag) &
               & fldDailyTempr2m => dispBuf%p2d(iQ)%present%ptr
          case(int_missing)  ! all models scanned
            exit
          case default
            call set_error('Unknown meteo-depence model:' + fu_str(arMeteoDepModels(iMdl)), &
                         & 'prepare_injection_meteodep')
            return
        end select   ! dependence model type
      enddo   ! dependence models
    enddo  ! quantities
    
  end subroutine prepare_injection_meteodep
  
  
  !********************************************************************************************
  
  real function fu_meteodep_emis_scaling(rulesMD, iDisp)
    !
    ! Determines the scaling of emission in what concerns the daily-mean temperature
    !
    implicit none
    
    ! Imported parameters
    type(TmeteoDepRules), intent(in) :: rulesMD
    integer, intent(in) :: iDisp    ! = ixDisp + nx_disp * (iyDisp-1)

    ! Local variables
    integer :: iMdl
    real :: fTmp
    
    ! Action depends on the model
    !
    select case(rulesMD%indMeteoDepModelSwitch)
      case(emis_heating_vs_T_linregr_mdl)
        !
        ! regression from temperature deficit to emission. A small hole is that the emission
        ! is scaled with regard to yesterday's temperature but it might even make sense: houses
        ! have inertia.
        !
        fTmp = rulesMD%fMinEmission
        fu_meteodep_emis_scaling = fTmp  + (1. - fTmp) * &                ! min emission level
                                 & rulesMD%fdEdT * max(0., rulesMD%fT0 - fldDailyTempr2m(iDisp)) 
      case(int_missing)  ! all models scanned
        return
      case default
        call set_error('Unknown meteo-depence model:' + fu_str(rulesMD%indMeteoDepModelSwitch), &
                     & 'fu_meteodep_emis_scaling')
        return
    end select   ! dependence model type
      
  end function fu_meteodep_emis_scaling
  
  
END MODULE source_terms_time_params

