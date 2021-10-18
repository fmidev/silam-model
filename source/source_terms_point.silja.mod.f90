MODULE source_terms_point

  ! This module contains the general description source-term of SILAM
  ! Source-term describes the spatial- and time- distribution of the 
  ! release to atmosphere, and the amount of chemical/radioactive 
  ! materials released.
  !
  ! In SILAM the source term consists of an unlimited number of
  ! release sources, each describing one release
  !
  ! In this module we define:
  !
  ! Words SOURCE and SOURCE TERM REFER TO FULL DESCRIPTION OF THE
  ! RELEASE SITUATION AND CIRCUMSTANCES (location, amounts of
  ! materials released, heat contents, spatial/vertical/time
  ! distribution etc...)
  !
  ! WORD "RELEASE" REFERS ONLY TO THE MATERIALS RELEASED TO THE
  ! ATMOSPHERE (material and their amounts)
  ! 
  ! NOTE. Co-ordinates of the point source (position) are always stored
  !       in the geographical grid as in ini file. System grid is not yet
  !       defined. BUT the cloud from this source will already be in
  !       system grid. This very point is a border between general
  !       geographical grids and the system grid, usually derived by
  !       the meteorological ones.
  !
  ! Currently the module contains descriptions for POINT SOURCE, BOMB SOURCE
  ! and AREA_SOURCE
  ! There can be two point source description formats - simple for the SILAM v.2
  ! and more sophisticated for the SILAM version 3. They called point_source_2 and
  ! point_source_3 correspondingly.
  !
  ! All units: SI
  ! 
  ! All degrees: real numbers in degrees (degrees and parts per
  ! hundred/parts per thousand , NOT in degrees and minutes)
  ! 
  ! Language: ANSI standard Fortran 90
  !
  ! Code owner: Mikhail Sofiev, FMI
  ! 
  ! Modules used:

  USE source_terms_time_params
!  USE source_terms_bomb  !time_params

  IMPLICIT NONE
  public                ! the only source term that has this open

  ! The public functions and subroutines available in this module:
  
  ! Functions universal for all types of the sources
  !
  public reserve_point_source
  public add_input_needs_point_source
  PUBLIC report
  public defined
  PUBLIC fu_name
  public fu_sector
  PUBLIC fu_start_time
  PUBLIC fu_end_time
  public fu_SlotTime
  PUBLIC fu_duration
  public total_amt_species_unit
  public fu_NbrOfTimeSlots
  PUBLIC source_vertical
  public fu_source_id_nbr
  public fu_source_nbr
  public fu_useTimeVarCoef  ! logical - if use the coefs
  public getTimeSlots
  public fu_cocktail_descr
  public fill_p_src_from_namelist
  public store_source_as_namelist
  public link_source_to_species
  public add_source_species
  public source_2_map
  public project_p_src_2_grids
  public prepare_src_vert_params_p_src
  public create_src_cont_grd_p_src
  public prepare_inject_p_src
  public inject_emission_euler_p_src
  public inject_emission_lagr_p_src
  public add_inventory_p_src
  public store_as_namelist
  public fu_meteodep_model


  PRIVATE report_point_source
  private fu_point_source_defined

  private store_p_src_as_namelist
  PRIVATE fu_point_source_start_time
  PRIVATE fu_point_source_end_time
  PRIVATE fu_point_source_duration
  PRIVATE point_source_vertical
  private fu_name_p_src
  private fu_sector_p_src
  private fu_NbrOfTimeSlots_of_p_src
  private total_from_p_src_descr_unit
  private total_from_p_src_species_unit
  private getTimeSlots_of_point_source
  private fu_SlotTime_of_point_source
  private fu_cocktail_descr_of_p_src
  private fu_source_id_nbr_of_p_src
  private point_source_rate_descr_unit
  
  private source_2_map_point_source
  private fu_useTimeVarCoef_p_src
  private link_p_src_to_species
  private add_source_species_p_src
  private fu_source_nbr_of_p_src
  private calc_plume_rise
  private determine_release_params
  private fu_meteodep_model_p_src

  ! Generic names and operator-interfaces of some functions:

  INTERFACE report
    MODULE PROCEDURE report_point_source
  END INTERFACE

  INTERFACE store_as_namelist
    MODULE PROCEDURE store_p_src_as_namelist
  END INTERFACE
 
  interface defined
    module procedure fu_point_source_defined
  end interface

  INTERFACE fu_name
    module procedure fu_name_p_src
  END INTERFACE

  INTERFACE fu_sector
    module procedure fu_sector_p_src
  END INTERFACE

  INTERFACE fu_start_time
    MODULE PROCEDURE fu_point_source_start_time
  END INTERFACE

  INTERFACE fu_end_time
    MODULE PROCEDURE fu_point_source_end_time
  END INTERFACE

  INTERFACE fu_duration
    MODULE PROCEDURE fu_point_source_duration
  END INTERFACE

  interface total_amt_species_unit
    module procedure total_from_p_src_species_unit
  end interface

  interface fu_NbrOfTimeSlots
    module procedure fu_NbrOfTimeSlots_of_p_src
  end interface 

  interface getTimeSlots
    module procedure getTimeSlots_of_point_source
  end interface

  interface fu_SlotTime
    module procedure fu_SlotTime_of_point_source
  end interface

  interface fu_cocktail_descr
    module procedure fu_cocktail_descr_of_p_src
  end interface

  INTERFACE source_vertical
    MODULE PROCEDURE point_source_vertical
  END INTERFACE

  interface fu_source_id_nbr
    module procedure fu_source_id_nbr_of_p_src
  end interface

  interface fu_source_nbr
    module procedure fu_source_nbr_of_p_src
  end interface

  interface link_source_to_species
    module procedure link_p_src_to_species
  end interface

  interface add_source_species
    module procedure add_source_species_p_src
  end interface

  interface source_2_map
    module procedure source_2_map_point_source
  end interface

  interface fu_useTimeVarCoef
    module procedure fu_useTimeVarCoef_p_src
  end interface

  interface store_source_as_namelist
    module procedure store_p_src_as_namelist
  end interface

  interface fu_meteodep_model
    module procedure fu_meteodep_model_p_src
  end interface

  integer, private, parameter :: singleLevelDynamic = 7771
  integer, private, parameter :: multiLevelFixed = 7772  
  integer, private, parameter :: plumeRise = 7773

  integer, private, parameter :: LocalTimeIndex = -55555 ! Not actual index
  integer, private, parameter :: SolarTimeIndex = -66666
  
  
  ! Public types with private components defined in this module:
  ! The main SILAM source term. Contains pointers to all possible
  ! source types. There may be plenty of the sources of each type
  !
  ! ATTENTION!! Each source (except for bomb one) is uniquely identified
  !             by means of the following parameters: src_nm + sector_nm
  !             or src_nbr. src_nbr is absolutely unique, while src_nm and
  !             sector_nm each may show up several times (but their
  !             combination is unique). This trick allows the following:
  !          1. Complete separation of all source-sector pairs
  !          2. Split of sources, with sectors summed-up
  !          3. Split of sectors, with sources summed-up
  !             This informaion is stored in iSrcIdType field, one for all sources.
  !
  ! Point source capable to handle varying strength and composition,
  ! as well as parameters for the plume rise routines
  !
  TYPE silam_point_source
    PRIVATE
    CHARACTER(len=clen) :: src_nm, sector_nm ! Name of the point source and sector
    integer :: src_nbr, id_nbr               ! Just a source and id numbers in a WHOLE source list
    integer :: nDescriptors                  ! Total over all time slots
    integer :: vertDispType
    logical :: ifUseTimeVarCoef, ifRateVary, ifEmissionLagrangian, if_inside_domain
    real :: lon, lat, fXDispGrd, fYDispGrd, stackHeight
    type(silam_release_parameters),dimension(:), pointer :: params
    type(Tcocktail_descr), dimension(:), pointer :: cocktail_descr_lst   ! for the whole source
    type(silam_vertical) :: vertLevs, vertLevsDispVert
    integer, dimension(:), pointer :: nSpeciesInDescr
    real, dimension(:,:), pointer :: fDescr2SpeciesUnit
    integer, dimension(:,:), pointer :: pEmisSpeciesMapping
    real, dimension(:), pointer :: levFraction, levFractDispVert, fzDisp, dz_m
    real, dimension(:,:), pointer :: indHour     ! Hourly rate indices
    real, dimension(:,:), pointer :: indDay      ! Daily rate indices
    real, dimension(:,:), pointer :: indMonth    ! Monthly rate indices
    integer :: ixDispGrd, iyDispGrd, izStackDispVert, nzDispVert
    character(len=tz_name_length) :: tz_name  ! Timezone name 
    integer :: tz_index ! Index of timezone to use with hourly/daily indices
                    ! can be LocalTimeIndex  or SolarTimeIndex
    ! Environment for meteo-dependent variation
    logical, public :: ifMeteoDependent            ! monthly variation coefficients are ignored
    type(TmeteoDepRules), pointer :: rulesMeteoDep => null()
    type(silja_logical) :: defined
  END TYPE silam_point_source

  !             
  type p_src_ptr
    TYPE(silam_point_source) :: p_src
  end type p_src_ptr

  type(field_4d_data_ptr), private, pointer, save :: fldTempr, fldWind, fldN2, fldPotTemp, &
                                                   & fldPressure, fldRdown, fldZ
  type(field_2d_data_ptr), private, save, pointer :: fldAblHeight, fldTZindex
  integer, private, save :: ind_height

CONTAINS

  !=========================================================================  
  !=========================================================================  
  !
  !           POINT source functions 
  !
  !=========================================================================  
  !=========================================================================  


  !**************************************************************************

  subroutine reserve_point_source(p_src, &        ! Src to initialise
                                & iSrcNbr, &      ! Src number in the grand list
                                & iSrcIdNbr, &    ! SrcID number
                                & nDescriptors, & ! number of chemical descriptors to be reserved
                                & iDynamicEnvironment)
    !
    ! Initialises the source:
    ! - stores the reference information: source number and source ID number
    ! - stores the total number of chemical descriptors that will be stored in the source
    ! - nullifies the source dynamic arrays 
    ! - set a few basic internal variables of the source
    !
    implicit none

    ! Imported parameters
    type(silam_point_source), intent(inout) :: p_src
    integer, intent(in) :: iSrcNbr, iSrcIdNbr, nDescriptors, iDynamicEnvironment

    !
    ! Nullify the basic variables
    !
    p_src%src_nm = ''
    p_src%sector_nm = ''
    p_src%stackHeight = real_missing
    p_src%vertDispType = int_missing
    nullify(p_src%params)
    nullify(p_src%levFraction)
    nullify(p_src%dz_m)
    nullify(p_src%fzDisp)
    nullify(p_src%levFractDispVert)
    nullify(p_src%cocktail_descr_lst); allocate(p_src%cocktail_descr_lst(nDescriptors))
    allocate(p_src%indHour(nDescriptors,24), p_src%indDay(nDescriptors,7), &
           & p_src%indMonth(nDescriptors,12))
    p_src%indHour = int_missing
    p_src%indDay = int_missing
    p_src%indMonth = int_missing
    p_src%tz_index = 0
    p_src%tz_name = ''
    nullify(p_src%nSpeciesInDescr)
    nullify(p_src%fDescr2SpeciesUnit)
    nullify(p_src%pEmisSpeciesMapping)
    !
    ! Main source parameters - enough to identify it in the global information list
    !
    p_src%src_nbr = iSrcNbr
    p_src%id_nbr = iSrcIdNbr
    p_src%nDescriptors = nDescriptors
    !
    ! Emission goes into Eulerian or Lagrangian environments
    !
    p_src%ifEmissionLagrangian = fu_if_lagrangian_present(iDynamicEnvironment)
    !
    ! Finally, mark the source as incomplete
    !
    p_src%defined = silja_false

  end subroutine reserve_point_source


  !*************************************************************************

  subroutine add_input_needs_point_source(p_src, q_met_dynamic, q_met_static, &
                                            & q_disp_dynamic, q_disp_static)
    !
    ! Returns input needs for the emission/transformation routines. 
    !
    implicit none

    ! Imported parameters
    type(silam_point_source), intent(in) :: p_src
    integer, dimension(:), intent(inout) :: q_met_dynamic, q_met_static, &
                                            & q_disp_dynamic, q_disp_static

    ! Local variables
    integer :: iTmp

    !
    ! Meteo data may be needed if the plume rise computations show up
    !
    if (p_src%vertDispType == plumeRise) then
      iTmp = fu_merge_integer_to_array(abl_height_m_flag, q_met_dynamic)
      iTmp = fu_merge_integer_to_array(temperature_flag, q_met_dynamic)
      iTmp = fu_merge_integer_to_array(windspeed_flag, q_met_dynamic)
      iTmp = fu_merge_integer_to_array(brunt_vaisala_freq_flag, q_met_dynamic)
      iTmp = fu_merge_integer_to_array(potential_temperature_flag, q_met_dynamic)
    endif
    !
    ! ... or for Lagrangian pre-distribution of particles
    !
    if(p_src%ifEmissionLagrangian)then
      iTmp = fu_merge_integer_to_array(abl_height_m_flag, q_met_dynamic)
      iTmp = fu_merge_integer_to_array(height_flag, q_met_dynamic)
      iTmp = fu_merge_integer_to_array(R_down_meteo_flag, q_met_dynamic)
      iTmp = fu_merge_integer_to_array(pressure_flag, q_met_dynamic)
    endif
    !
    ! ... or for Local timezone for source location
    !
    if(p_src%tz_index == LocalTimeIndex)then !Use local time
      iTmp = fu_merge_integer_to_array(timezone_index_flag, q_met_static)
    endif
    return

  end subroutine add_input_needs_point_source


  ! ***************************************************************

  subroutine fill_p_src_from_namelist(nlSrc, p_src, chPointSrcFileVersion)
    !
    ! Reads and sets one area source term v.2 from an external file.
    ! A trick:
    ! Since the area source is potentially huge, it can be moved to a 
    ! separate file. In this case in the main file we will just have 
    ! a name of the file to read. 
    ! Note also that we rely on the nDescriptors value to be set, as well the 
    ! number of the source in the source info list (available for the general
    ! source only).
    !
    ! All units: SI, unless otherwise is stated
    !
    IMPLICIT NONE
    !
    ! Imported parameters
    TYPE(silam_point_source), intent(inout) :: p_src
    type(Tsilam_namelist), pointer :: nlSrc
    character(len=*), intent(in) :: chPointSrcFileVersion   ! at present, v2 or v3

    ! Local variables
    type(silam_sp) :: spContent
    integer :: iLocal, iTmp, nx, ny, nItems, iStat, iCell, iDescr, &
             & nDescriptorInParLine, iModeNbr, iItem, iPointSrcFileVersion
    integer, dimension(:), pointer :: indDescriptorOfParLine    ! places of descriptors in descriptor list
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: ptrItems
    real :: fLon, fLat
    logical :: ifGeoCoords ! Whether values are in geographical or grid co-oprdinates
    logical, dimension(:), allocatable :: ifSpecies
    type(silam_grid_position) :: posTmp
    type(silja_level) :: level
    character(len=substNmLen) :: chSubstNm
    real, dimension(:), pointer :: arTmp

    ! A bit of preparations
    !
    iPointSrcFileVersion = int_missing
    spContent%sp => fu_work_string()
    nullify(ptrItems)
    indDescriptorOfParLine => fu_work_int_array()
    arTmp => fu_work_array()
    if(error)return
    !
    ! If the source header has not yet been defined, define it. Otherwise check that we 
    ! got the right one
    !
    if(.not. (p_src%defined == silja_true))then
      call set_header()
      if(error)return
    else

      if(.not. fu_same_header())then
        call msg_warning('The header and the namelist are different','fill_p_src_from_namelist')
        call report(p_src)
        call report(nlSrc)
        call set_error('The header and the namelist are different','fill_p_src_from_namelist')
        return
      endif
    endif

    !
    ! Timeslot parameters
    !
    select case(chPointSrcFileVersion)
      case('4')
        iPointSrcFileVersion = 4
        call get_items(nlSrc, 'par_str', ptrItems, nItems)
      case('5')
        iPointSrcFileVersion = 5
        call get_items(nlSrc, 'par_str_point', ptrItems, nItems)
      case default
        call msg('Wrong point source file version (must be 4 or 5): "'// &
                    & trim(chPointSrcFileVersion)// '"')
        call set_error('Wrong point source file version','fill_p_src_from_namelist')
        return
    end select

    if(nItems < 1)then
      call report(nlSrc)
      call set_error('No parameter items in namelist','fill_p_src_from_namelist')
      return
    endif

    if(.not. p_src%defined == silja_true)then
      !
      ! If nothing has been read yet, param vector is null. Have to allocate it
      !
      allocate(p_src%params(nItems), stat=iTmp)
      if(iTmp /= 0)then
        call set_error('Failed to allocate parameters for area source', &
                     & 'fill_p_src_from_namelist')
        return
      endif
      call set_missing(p_src%params, nItems)
    endif

    spContent%sp=fu_content(nlSrc,'vertical_unit')
    
    if(error)return
!call msg("p_src%nDescriptors",p_src%nDescriptors)
    indDescriptorOfParLine(1:p_src%nDescriptors) = int_missing
    chSubstNm = fu_content(nlSrc,'emitted_substance')
    iModeNbr = fu_content_int(nlSrc,'emitted_size_mode_nbr')

    do iTmp = 1, nItems
      select case(iPointSrcFileVersion)
        case(4)
          call decode_time_parameter_v4 (fu_content(ptrItems(iTmp)), &  ! line to decode
                                       & iTmp, &                        ! param element index to write to
                                       & spContent%sp, &                ! unit of the vertical
                                       & fu_content(nlSrc,'release_rate_unit'), & ! source release rate unit
                                       & p_src%params, &                ! param array
                                       & p_src%cocktail_descr_lst, &    ! list of source cocktail descriptors
                                       & point_source, &                 ! type of the source
                                       & p_src%nDescriptors, &          ! max nbr of descriptors in the source
                                       & indDescriptorOfParLine, nDescriptorInParLine, & 
                                       & chSubstNm, & 
                                       & iModeNbr, ifSpecies)
        case(5)
          call decode_time_parameter_v5 (fu_content(ptrItems(iTmp)), &
                                       & iTmp, &
                                       & spContent%sp, &
                                       & fu_content(nlSrc,'release_rate_unit'), & ! source release rate unit
                                       & p_src%params, &
                                       & p_src%cocktail_descr_lst, &    ! list of source cocktail descriptors
                                       & point_source, &   ! type of the source
                                       & p_src%nDescriptors, &
                                       & indDescriptorOfParLine, nDescriptorInParLine, &
                                       & '', arTmp, ifSpecies)
        case default
          call msg('Wrong area source file version (must be 4 or 5):',iPointSrcFileVersion)
          call set_error('Wrong area source file version','fill_p_src_from_namelist')
          return
      end select
      if(error)then
        call set_error('Failed source:' + p_src%src_nm + '_' + p_src%sector_nm,'fill_p_src_from_namelist')
        return
      endif
    end do   ! time slots

    !
    ! Start location & time. 
    ! Trick. Position deemed to be in relative units, but if input_grid 
    ! is omitted then it just stores what is given as it is. So, the
    ! POINT_SOURCE IS IN GEOGRAPHICAL GRID - as in input file
    ! Reason - we do not have system_grid yet => nowhere to project
    !
!    p_src%position = fu_set_pos_in_geo_global_grid(fu_content_real(nlSrc,'source_longitude'), &
!                                                 & fu_content_real(nlSrc,'source_latitude'), &
!                                                 & 92500.,p_src%params(1)%time)
    p_src%lon = fu_content_real(nlSrc,'source_longitude')
    p_src%lat = fu_content_real(nlSrc,'source_latitude')


    
    p_src%tz_name = fu_content(nlSrc,'source_timezone')
    if( p_src%tz_name/= '')then 
            if(fu_str_u_case(p_src%tz_name) == 'LOCAL')then
                     ! Use local time at the source location
                     p_src%tz_index = LocalTimeIndex
            elseif (fu_str_u_case(p_src%tz_name) == 'SOLAR')then
                     p_src%tz_index = SolarTimeIndex
            else 
                     p_src%tz_index =  timezone_index_by_name(p_src%tz_name)
            endif
    else 
            p_src%tz_name = 'Solar'
            p_src%tz_index = SolarTimeIndex
            !            call msg('No timezone for source. Assuming solar time.')
!            call msg('No timezone for source:'+p_src%src_nm+','+ p_src%sector_nm +'. Assuming solar time.')
    endif 


    !
    ! It is useful to know if the composition varies with time - just a speedup
    !
    SLOTS: do iTmp = 1, size(p_src%params)-1
      do iDescr = 1, p_src%nDescriptors
        p_src%ifRateVary = .not. (p_src%params(iTmp)%rate_descr_unit(iDescr) .eps. &
                                & p_src%params(iTmp+1)%rate_descr_unit(iDescr))
        if(p_src%ifRateVary)exit SLOTS
      end do
    end do SLOTS

    !
    ! If time variation indices are present - read and set them, otherwise - set all to 1
    !
    select case(iPointSrcFileVersion)
      case(4)
        !
        ! Same species can already be met in other files. Check the the same variation
        !
        if(p_src%indHour(indDescriptorOfParLine(1),1) >= 0.0)then ! something reasonable found
          call set_error('The species:' + fu_name(p_src%cocktail_descr_lst(indDescriptorOfParLine(1))) + &
                 & '- has already been assigned time variation for the source:' + p_src%src_nm + '_' + &
                 & p_src%sector_nm,'fill_p_src_from_namelist')
        endif
        !
        ! Put the variation to the corresponding indDescriptorOfParLine
        !
        spContent%sp = fu_content(nlSrc,'hour_in_day_index')
        if(spContent%sp /= '')then
          read(unit=spContent%sp, iostat=iLocal, fmt=*) &
                                 & (p_src%indHour(indDescriptorOfParLine(1),iTmp),iTmp=1,24)
          if(iLocal /= 0)then
            call set_error('Failed to read  hourly variation:' + spContent%sp,'set_header')
            return
          endif
!          call msg('Hourly variation indices are read. sum=',real_value=sum(a_src%indHour(1:24)))
        else
          p_src%indHour(indDescriptorOfParLine(1),:) = 1
        endif
        spContent%sp = fu_content(nlSrc,'day_in_week_index')
        if(spContent%sp /= '')then
          read(unit=spContent%sp, iostat=iLocal, fmt=*) &
                                    & (p_src%indDay(indDescriptorOfParLine(1),iTmp),iTmp=1,7)
          if(iLocal /= 0)then
            call set_error('Failed to read  day-in-week variation:' + spContent%sp,'set_header')
            return
          endif
!          call msg('Daily variation indices are read. sum=',real_value=sum(a_src%indDay(1:24)))
        else
          p_src%indDay(indDescriptorOfParLine(1),:) = 1
        endif
        spContent%sp = fu_content(nlSrc,'month_in_year_index')
        if(spContent%sp /= '')then
          read(unit=spContent%sp, iostat=iLocal, fmt=*) &
                                         & (p_src%indMonth(indDescriptorOfParLine(1),iTmp),iTmp=1,12)
          if(iLocal /= 0)then
            call set_error('Failed to read  monthly variation:' + spContent%sp,'set_header')
            return
          endif
!          call msg('Monthly variation indices are read. sum=',real_value=sum(a_src%indMonth(1:24)))
        else
          p_src%indMonth(indDescriptorOfParLine(1),:) = 1
        endif
        p_src%ifUseTimeVarCoef = p_src%ifUseTimeVarCoef .or. &
                               & any(p_src%indHour(indDescriptorOfParLine(1),1:24) < 0.999) .or. &
                               & any(p_src%indHour(indDescriptorOfParLine(1),1:24) > 1.0001) .or. &
                               & any(p_src%indDay(indDescriptorOfParLine(1),1:7) < 0.999) .or. &
                               & any(p_src%indDay(indDescriptorOfParLine(1),1:7) > 1.0001) .or. &
                               & any(p_src%indMonth(indDescriptorOfParLine(1),1:12) < 0.999) .or. &
                               & any(p_src%indMonth(indDescriptorOfParLine(1),1:12) > 1.0001)
      case(5)
        !
        ! POINT_SOURCE_5: several descriptors in one file.
        ! Same species can already be met in other files. Check the the same variation
        !
        do iDescr = 1, nDescriptorInParLine
          if(p_src%indHour(indDescriptorOfParLine(iDescr),1) >= 0.0)then ! something reasonable found
            call set_error('The species:' + &
                   & fu_name(p_src%cocktail_descr_lst(indDescriptorOfParLine(iDescr))) + &
                   & '- has already been assigned time variation for the source:' + p_src%src_nm + '_' + &
                   & p_src%sector_nm,'fill_p_src_from_namelist')
          endif
        enddo
        !
        ! Put the variation to the corresponding indDescriptorOfParLine
        !
        call get_items(nlSrc,'hour_in_day_index',ptrItems,nItems)
        do iDescr = 1, nDescriptorInParLine
          p_src%indHour(indDescriptorOfParLine(iDescr),:) = 1.
          do iItem = 1, nItems
            spContent%sp = fu_content(ptrItems(iItem))
            read(unit=spContent%sp,fmt=*)chSubstNm
            if(chSubstNm == trim(fu_name(p_src%cocktail_descr_lst(indDescriptorOfParLine(iDescr))))) then
              read(unit=spContent%sp, iostat=iLocal, fmt=*) chSubstNm, &
                                     & (p_src%indHour(indDescriptorOfParLine(iDescr),iTmp),iTmp=1,24)
              if(iLocal /= 0)then
                call set_error('Failed to read  hourly variation:' + spContent%sp,'set_header')
                return
              endif
              exit
            endif
          end do
        end do

        call get_items(nlSrc,'day_in_week_index',ptrItems,nItems)
        do iDescr = 1, nDescriptorInParLine
          p_src%indDay(indDescriptorOfParLine(iDescr),:) = 1.
          do iItem = 1, nItems
            spContent%sp = fu_content(ptrItems(iItem))
            read(unit=spContent%sp,fmt=*)chSubstNm
            if(chSubstNm == trim(fu_name(p_src%cocktail_descr_lst(indDescriptorOfParLine(iDescr))))) then
              read(unit=spContent%sp, iostat=iLocal, fmt=*) chSubstNm, &
                                     & (p_src%indDay(indDescriptorOfParLine(iDescr),iTmp),iTmp=1,7)
              if(iLocal /= 0)then
                call set_error('Failed to read  daily-in-week variation:' + spContent%sp,'set_header')
                return
              endif
              exit
            endif
          end do
        end do

        call get_items(nlSrc,'month_in_year_index',ptrItems,nItems)
        do iDescr = 1, nDescriptorInParLine
          p_src%indMonth(indDescriptorOfParLine(iDescr),:) = 1.
          do iItem = 1, nItems
            spContent%sp = fu_content(ptrItems(iItem))
            read(unit=spContent%sp,fmt=*)chSubstNm
            if(chSubstNm == trim(fu_name(p_src%cocktail_descr_lst(indDescriptorOfParLine(iDescr))))) then
              read(unit=spContent%sp, iostat=iLocal, fmt=*) chSubstNm, &
                                     & (p_src%indMonth(indDescriptorOfParLine(iDescr),iTmp),iTmp=1,12)
              if(iLocal /= 0)then
                call set_error('Failed to read month-in-year variation:' + spContent%sp,'set_header')
                return
              endif
              exit
            endif
          end do
        end do

        do iDescr = 1, nDescriptorInParLine
          p_src%ifUseTimeVarCoef = p_src%ifUseTimeVarCoef .or. &
                           & any(p_src%indHour(indDescriptorOfParLine(iDescr),1:24) < 0.999) .or. &
                           & any(p_src%indHour(indDescriptorOfParLine(iDescr),1:24) > 1.0001) .or. &
                           & any(p_src%indDay(indDescriptorOfParLine(iDescr),1:7) < 0.999) .or. &
                           & any(p_src%indDay(indDescriptorOfParLine(iDescr),1:7) > 1.0001) .or. &
                           & any(p_src%indMonth(indDescriptorOfParLine(iDescr),1:12) < 0.999) .or. &
                           & any(p_src%indMonth(indDescriptorOfParLine(iDescr),1:12) > 1.0001)
        enddo

      case default
        call set_error('Only POINT_SOURCE_4 and POINT_SOURCE_5 are supported','fill_p_src_from_namelist')
        return
    end select

    call destroy_items(ptrItems)
    if(error)return

    !
    ! Things have been done but do NOT destroy the namelist - it was NOT you who made it, stupid!
    !
    call free_work_array(spContent%sp)
    call free_work_array(indDescriptorOfParLine)
    call free_work_array(arTmp)

    if(p_src%sector_nm == '')then
      call msg('Point source name:' + p_src%src_nm)
    else
      call msg('Point source name & sector:' + p_src%src_nm + ',' + p_src%sector_nm)
    endif

    p_src%defined = silja_true


    CONTAINS

    !====================================================================================

    subroutine set_header()
      !
      ! Sets the basic parameters of the source
      !
      implicit none

      !
      ! Source name
      !
      p_src%src_nm = fu_Content(nlSrc, 'source_name')
      p_src%sector_nm = fu_Content(nlSrc, 'source_sector_name')

!      call msg(fu_connect_strings('setting the source:',p_src%src_nm))

!      if(p_src%src_nm == 'FI')then
!        call msg('')
!      endif
       p_src%ifMeteoDependent = fu_str_u_case(fu_content(nlSrc,'if_temperature_dependent_emission')) == 'YES'
       if (p_src%ifMeteoDependent) then
         call msg("Source: "//trim(p_src%src_nm))
         call msg("Sector: "//trim(p_src%sector_nm))
         call set_error("Meteo_dependent point source not implemented yet", 'set_header')
         return
       endif


  
      select case(fu_str_u_case(fu_content(nlSrc,'vertical_distribution')))
        case ('SINGLE_LEVEL_DYNAMIC')
          p_src%vertDispType = singleLevelDynamic
          call set_missing(p_src%vertLevs, .true.)
          nullify(p_src%levFraction)
          nullify(p_src%dz_m)
      
        case ('MULTI_LEVEL_FIXED')
          p_src%vertDispType = multiLevelFixed
          !
          ! Read levels and set the source vertical
          !
          call get_items(nlSrc, 'vert_level', ptrItems, nItems)
          spContent%sp=fu_content(nlSrc,'vertical_unit')
          if(nItems < 1)then
            call set_error('No levels for vertical_distribution = MULTI_LEVEL_FIXED', &
                         & 'set_header')
            return
          else
            allocate(p_src%levFraction(nItems), p_src%dz_m(nItems))
            call set_named_level_with_fract(ptrItems(1), level, p_src%levFraction(1), spContent%sp)
            if(error)return
            call set_vertical(level, p_src%vertLevs)
            if(error)return
            do iTmp = 2, nItems
              call set_named_level_with_fract(ptrItems(iTmp), level, p_src%levFraction(iTmp), spContent%sp)
              call add_level(p_src%vertLevs, level)
              if(error)return
            end do
          endif
        
        case ('PLUME_RISE')
          p_src%vertDispType = plumeRise          
          p_src%stackHeight = fu_set_named_value(fu_content(nlSrc,'stack_height'))
          if (p_src%stackHeight .eps. real_missing) then
          call set_error('stack_height needed in plume rise mode.', &
                          & 'fill_p_src_from_namelist')
          call set_missing(p_src%vertLevs, .true.)
          nullify(p_src%levFraction)
          nullify(p_src%dz_m)
        endif
        case default
          p_src%vertDispType = singleLevelDynamic ! Default value
          ! call set_error('Unknown vertical distripution type:' + fu_content(nlSrc,'vertical_distribution'), &
          !        & 'fill_p_src_from_namelist')
      end select    
      

    
    end subroutine set_header


    !==============================================================================

    logical function fu_same_header()
      !
      ! Compares the current source header and the one defined in the namelist
      !
      implicit none

      ! Local variables
      type(silam_release_parameters), dimension(1), target :: paramTmp
      type(silam_release_parameters), dimension(:), pointer :: paramPtr
      type(silja_level) :: levTmp
      real :: fTmp
      fu_same_header = .false.

      !
      ! Name
      !
      if(.not. (fu_str_u_case(p_src%src_nm) == fu_str_u_case(fu_Content(nlSrc, 'source_name'))))then
        call msg_warning(fu_connect_strings('Source names are different:', &
                                         & p_src%src_nm,',',fu_Content(nlSrc, 'source_name')), &
                       & 'fu_same_header')
        call set_error('Cannot collect species between different sources','fu_same_header')
        return
      endif

      ! Sector
      if(.not. (fu_str_u_case(p_src%sector_nm) == &
                                      & fu_str_u_case(fu_Content(nlSrc, 'source_sector_name'))))then
        call msg_warning(fu_connect_strings('Source sector names are different:', &
                                         & p_src%sector_nm,',',fu_Content(nlSrc, 'source_sector_name')), &
                       & 'fu_same_header')
        call set_error('Cannot collect species between different sources','fu_same_header')
        return
      endif

     ! Vertical distribution
      if(.not. ((p_src%vertDispType == singleLevelDynamic .and. fu_str_u_case(fu_Content(nlSrc, 'vertical_distribution')) == 'SINGLE_LEVEL_DYNAMIC') .or. &
              & (p_src%vertDispType == multiLevelFixed .and. fu_str_u_case(fu_Content(nlSrc, 'vertical_distribution')) == 'MULTI_LEVEL_FIXED') .or. &
              & (p_src%vertDispType == plumeRise .and. fu_str_u_case(fu_Content(nlSrc, 'vertical_distribution')) == 'PLUME_RISE')))then
        call msg_warning('Vertical distribution types are different:' + &
                       & fu_str(p_src%vertDispType) + ',' + fu_Content(nlSrc, 'vertical_distribution'), &
                       & 'fu_same_header')
        call set_error('Cannot collect species between different sources','fu_same_header')
        return
      endif
      select case(p_src%vertDispType)
        case(singleLevelDynamic)
          continue
        case(multiLevelFixed)
          call get_items(nlSrc, 'vert_level', ptrItems, nItems)
          spContent%sp=fu_content(nlSrc,'vertical_unit')
          if(.not. (size(p_src%levFraction) == nItems))then
            call msg_warning('Different number of levels in vertical distributions:' + &
                           & fu_str(size(p_src%levFraction)) + ', and,' + fu_str(nItems), &
                           & 'fu_same_header')
            call set_error('Cannot collect species between different sources','fu_same_header')
            return
          endif
          do iTmp = 1, nItems
            call set_named_level_with_fract(ptrItems(iTmp), levTmp, fTmp, spContent%sp)
            if(.not. (fu_cmp_levs_eq(fu_level(p_src%vertLevs, iTmp), levTmp) .and. p_src%levFraction(iTmp) == fTmp))then
               call msg('Vertical distributions are different', p_src%levFraction(iTmp), fTmp )
               call report(fu_level(p_src%vertLevs, iTmp))
               call report(levTmp)
               call set_error('Cannot collect species between different sources','fu_same_header')
            endif
          enddo
        case(plumeRise)
          if(.not. (p_src%stackHeight .eps.  fu_content_real(nlSrc, 'stack_height')))then
            call msg('Stack heights are different:', &
                                         & p_src%stackHeight,fu_content_real(nlSrc, 'stack_height'))
            call set_error('Cannot collect species between different sources','fu_same_header')
            return
          endif
        case default
          call set_error('Strange vertical distribution types','fu_same_header')
          return
        end select     

      ! Timeslot parameters - just check the number of slots. The list of descriptors will be
      ! enlarged later on
      !
      select case(iPointSrcFileVersion)
        case(4)
          call get_items(nlSrc, 'par_str', ptrItems, nItems)
        case(5)
          call get_items(nlSrc, 'par_str_point', ptrItems, nItems)
        case default
          call msg('Wrong area source file version (must be 2 or 3):',iPointSrcFileVersion)
          call set_error('Wrong area source file version','fu_same_header')
          return
      end select
      if(nItems /= size(p_src%params))then
        call msg('Different number of time slots in p_src and namelist:', nItems, real(size(p_src%params)))
        call set_error('Cannot collect species between different sources','fu_same_header')
        return
      endif

!      posTmp = fu_set_pos_in_geo_global_grid(fu_content_real(nlSrc,'source_longitude'), &
!                                           & fu_content_real(nlSrc,'source_latitude'), &
!                                           & fu_pressure(p_src%position), &
!                                           & fu_time(p_src%position))
!      if(.not. (fu_x(p_src%position) .eps. fu_x(posTmp)) .or. &
!       & .not. (fu_y(p_src%position) .eps. fu_y(posTmp)))then
      if(.not. (p_src%lon .eps. fu_content_real(nlSrc,'source_longitude')) .or. &
       & .not. (p_src%lat .eps. fu_content_real(nlSrc,'source_latitude')))then
        call msg('Longitude in the source and namelist:', p_src%lon, &
                                                        & fu_content_real(nlSrc,'source_longitude'))
        call msg('Latitude in the source and namelist:', p_src%lat, &
                                                       & fu_content_real(nlSrc,'source_latitude'))
        call set_error('Position in namelist and point source are different','fu_same_header')
        return
      endif

      !
      ! And now - the parameters themselves
      !
      paramPtr => paramTmp
      call set_missing(paramPtr, 1)

      do iTmp = 1, nItems
        select case(iPointSrcFileVersion)
          case(4)
            call decode_time_parameter_v4(fu_content(ptrItems(iTmp)), &  ! parameter itself
                                        & 1, &                           ! always 1 to get to paramTmp
                                        & 'm', &                         ! arbitrary, vertical is NOT CHECKED
                                        & fu_content(nlSrc,'release_rate_unit'), & ! source release rate unit
                                        & paramPtr, &                    ! receiving temporary
                                        & p_src%cocktail_descr_lst, &    ! list of source cocktail descriptors
                                        & point_source, &                 ! type of the source
                                        & p_src%nDescriptors, &
                                        & indDescriptorOfParLine, &
                                        & nDescriptorInParLine, &
                                        & fu_content(nlSrc,'emitted_substance'), &
                                        & fu_content_int(nlSrc,'emitted_size_mode_nbr'), ifSpecies)
          case(5)
            call decode_time_parameter_v5(fu_content(ptrItems(iTmp)), &  ! parameter itself
                                        & 1, &                           ! always 1 to get to paramTmp
                                        & 'm', &                         ! arbitrary, vertical is NOT CHECKED
                                        & fu_content(nlSrc,'release_rate_unit'), & ! source release rate unit
                                        & paramPtr, &                    ! receiving temporary
                                        & p_src%cocktail_descr_lst, &    ! list of source cocktail descriptors
                                        & point_source, &                 ! type of the source
                                        & p_src%nDescriptors, &
                                        & indDescriptorOfParLine, &
                                        & nDescriptorInParLine, &
                                        & '', arTmp, ifSpecies)
          case default
            call msg('Wrong area source file version (must be 4 or 5):',iPointSrcFileVersion)
            call set_error('Wrong point source file version','fill_p_src_from_namelist')
            return
        end select
        if(.not.(paramTmp(1)%time == p_src%params(iTmp)%time))then
          call msg_warning(fu_connect_strings('Failed time comparison for par_str:',fu_content(ptrItems(iTmp))), &
                         & 'fu_same_header')
          call report(paramTmp(1)%time)
          call report(p_src%params(iTmp)%time)
          call set_error('Cannot collect chemicals between sources with different parameters', &
                       & 'fu_same_header')
          return
        endif
        
        select case(p_src%vertDispType)
          case(singleLevelDynamic)
            if(.not. fu_cmp_levs_eq(paramTmp(1)%layerDynamic, p_src%params(iTmp)%layerDynamic))then
              call msg_warning(fu_connect_strings('Failed level comparison for par_str:',fu_content(ptrItems(iTmp))), &
                         & 'fu_same_header')
              call report(paramTmp(1)%layerDynamic)
              call report(p_src%params(iTmp)%layerDynamic)
              call set_error('Cannot collect chemicals between sources with different parameters', &
                       & 'fu_same_header')      
            endif
          case(PlumeRise)
            if(.not.(paramTmp(1)%xy_size .eps. p_src%params(iTmp)%xy_size) .or. &
             & .not.(paramTmp(1)%z_velocity .eps. p_src%params(iTmp)%z_velocity) .or. &
             & .not.(paramTmp(1)%tempr .eps. p_src%params(iTmp)%tempr))then
              call msg_warning(fu_connect_strings('Failed plumerise parameters comparison for par_str:',fu_content(ptrItems(iTmp))), &
                         & 'fu_same_header')
              call report(p_src)
              call set_error('Cannot collect chemicals between sources with different parameters', &
                       & 'fu_same_header')
              return
            endif
          case default
            continue
        end select
        
      end do  ! params items

      fu_same_header = .true.

    end function fu_same_header

  end subroutine fill_p_src_from_namelist


  ! ***************************************************************

  subroutine getTimeSlots_of_point_source(p_src, time, iSlot1, iSlot2)
    !
    ! Finds the time slots surrounding the given time
    ! 
    IMPLICIT NONE

    ! Imported parameters
    TYPE(silam_point_source), INTENT(in) :: p_src
    type(silja_time), intent(in) :: time
    integer, intent(out) :: iSlot1, iSlot2
    
    call getTimeSlots_from_params(p_src%params, time, iSlot1, iSlot2)
    
  END subroutine  getTimeSlots_of_point_source
  

  !=========================================================================  

  function fu_SlotTime_of_point_source(p_src, iSlot)
    implicit none
    ! Return value of the function
    type(silja_time) :: fu_SlotTime_of_point_source
    !
    ! Imported parameters with intent IN
    TYPE(silam_point_source), INTENT(in) :: p_src
    integer, intent (in) :: iSlot

    if(iSlot > size(p_src%params))then
      call set_error('Too big slot number','fu_SlotTime_of_point_source')
      return
    endif
    fu_SlotTime_of_point_source = p_src%params(iSlot)%time

  end function fu_SlotTime_of_point_source


  !=========================================================================  
  function fu_cocktail_descr_of_p_src (p_src) result(descrLst)
    implicit none
    type(Tcocktail_descr), dimension(:), pointer :: descrLst
    type(silam_point_source), intent(in), target :: p_src
    descrLst => p_src%cocktail_descr_lst
  end function fu_cocktail_descr_of_p_src


  !*****************************************************************

  subroutine point_source_rate_descr_unit(p_src, iTimeSlot, arRate)
    !
    ! Returns the release rate of the specific point source 
    ! pointed by p_src, at the time moment pointed by iTimeSlot
    !
    ! Note that now the rate is returned for each descriptor, NOT total.
    ! Neither it is the species-wise rate, i.e. these are the source-defined mixtures,
    ! not the emission mass map species (and certainly not the transport species).
    ! To get the emission species-wise rate use arep_src_rate_species.
    !
    ! Time variation coefs are NOT HERE. Rate_local_unit is just for 
    ! the overall slot rate, not for time fluctuations. Slot can be long
    !
    implicit none

    type(silam_point_source), intent(in) :: p_src
    integer, intent(in) :: iTimeSlot
    real, dimension(:), pointer :: arRate

    ! Local variables
    integer :: iDescr, iTmp

    if(iTimeSlot < 1 .or. iTimeSlot > size(p_src%params))then
      call msg('Time slot number: ', iTimeSlot)
      call set_error('Strange time slot number','fu_point_source_rate')
      return
    endif
    do iDescr = 1, p_src%nDescriptors
      arRate(iDescr) = p_src%params(iTimeSlot)%rate_descr_unit(iDescr)
    end do  ! iDescr

  end subroutine point_source_rate_descr_unit


  !*****************************************************************

  subroutine link_p_src_to_species(species_list, p_src)
    !
    ! Having the cocktails created, we should establish the shortcut links between the 
    ! source descriptors and cocktail. The link goes via descr%iEmisCocktSpeciesMapping 
    ! and  descr%factor_to_basic_unit.
    ! Note that cocktail is "anything" and the only requirement is that it has all the species
    ! existing in the source inventory.
    ! That has to happen in two steps. Firstly, we establish these links using the 
    ! single descriptor per source. Then, these connections are distributed to each time slot
    !
    implicit none

    ! Imported parameters
    type(silam_point_source), intent(inout) :: p_src
    type(silam_species), dimension(:), pointer :: species_list

    call link_src_species_to_given_list(p_src%nDescriptors, &
                                      & p_src%nSpeciesInDescr, &
                                      & p_src%cocktail_descr_lst, &
                                      & p_src%fDescr2SpeciesUnit, &
                                      & p_src%pEmisSpeciesMapping, &
                                      & species_list)

  end subroutine link_p_src_to_species


  !*******************************************************************

  subroutine add_source_species_p_src(p_src, speciesLst, nSpecies)
    !
    ! Fills-in the given list with the own species. Checks for the duplicates
    !
    implicit none

    ! Improted parameters
    type(silam_point_source), intent(inout) :: p_src
    type(silam_species), dimension(:), pointer :: speciesLst
    integer, intent(inout) :: nSpecies

    ! Local variables
    integer :: iDescr, nSpeciesDescr
    type(silam_species), dimension(:), pointer :: pSpecies

    do iDescr = 1, p_src%nDescriptors
      call get_inventory(p_src%cocktail_descr_lst(iDescr), pSpecies, nSpeciesDescr)
      if(error)return
      call addSpecies(speciesLst, nSpecies, pSpecies, nSpeciesDescr)
      if(error)return
    end do

  end subroutine add_source_species_p_src


  !*****************************************************************

  subroutine total_from_p_src_descr_unit(p_src, &                   ! mandatory, input
                                       & amounts, &                 ! mandatory, output
                                       & start_, duration_, layer_) ! optional, input
    !
    ! Returns the amount of the released material IN DESCRIPTOR UNIT starting from 
    ! start during the duration time interval.
    !
    ! We have to integrate the release rate, which is given in descriptor basic units - for each
    ! descriptor. 
    !
    ! Four cases are considered (notebook 11 pp.21-23)
    ! 1. slot rates are constant: descriptor rate and vertical fraction
    ! 2. place-connected time variation coefficients are constant (fHour, fDay, fMonth)
    ! 3. both sets are varying
    ! 4. None of them is time dependent
    !
    ! This subroutine produces a total descriptor-wise emission and returns a vector.
    !
    ! Units: SI basic and derived ones
    !
    implicit none

    ! Imported parameters with intent IN
    type(silam_point_source), intent(in) :: p_src
    real, dimension(:), intent(out) :: amounts
    type(silja_time), intent(in), optional :: start_
    type(silja_interval), intent(in), optional :: duration_
    type(silja_level), intent(in), optional :: layer_  ! if defined, include only fraction emitted in it

    ! Local variables
    integer :: iSlot, iStartSlot, iEndSlot, iTmp, iDescr, iCell
    real :: fTimeConversion, fTmp, fWeightPast, fVertFraction1, fVertFraction2, fTime_sec, &
          & fTime_sec_start, fTime_sec_end, fSlot_sec
    type(silja_time) :: start, start_integr, end_integr, startTmp, endTmp, timeTmp, timeTmpLocal
    type(silja_interval) :: duration, timeStep
    type(silja_interval) :: myTZoffset ! Offset from UTC
    logical :: ifSlotRateVary
    real, dimension(:), pointer :: fCells
    type(silja_level) :: layer

    !
    ! Time period for integration, then start and end slots
    !
    if(present(start_))then
      start = start_
      if(start > p_src%params(size(p_src%params))%time)then
        call report(p_src)
        call msg_warning('Source is earlier then the computation starts','total_from_p_src_descr_unit')
        return
      endif
      do iStartSlot = 1, size(p_src%params)-1
        if(p_src%params(iStartSlot+1)%time > start)exit ! Skip whole slots before the start time
      end do
    else
      start = p_src%params(1)%time
      iStartSlot = 1
    endif
    if(present(duration_))then
      duration = duration_
      if(start+duration > p_src%params(size(p_src%params))%time) &
                                      & duration = p_src%params(size(p_src%params))%time - start
      if(start+duration < p_src%params(1)%time)then
        call report(p_src)
        call msg_warning('Source is later then the computation ends','total_from_p_src_descr_unit')
        return
      endif
      do iEndSlot = iStartSlot+1, size(p_src%params)
        if(p_src%params(iEndSlot-1)%time < start+duration)exit ! Skip whole slots after the end time
      end do
      if(iEndSlot > size(p_src%params)) iEndSlot = size(p_src%params)
    else
      iEndSlot = size(p_src%params)
      duration = p_src%params(iEndSlot)%time - start
    endif
    if (present(layer_)) then
      layer = layer_
    else
      layer = level_missing
    end if
 
    amounts(1:p_src%nDescriptors) = 0.0
    !
    ! To have non-zero released amount, there must be at least two slots in the given range
    !
    if(iStartSlot >= iEndSlot)then
      return
    endif

    if(error)return
    fCells => fu_work_array()
    if(error)return
    fCells(1:p_src%nDescriptors) = 1.0

    !--------------------------------------------------------------------------------
    !
    ! Go with time-wise integration
    !
    do iSlot = iStartSlot, iEndSlot-1
      !
      ! Check the vertical layer - may be, we do not have to do anything - if there is no emission
      ! into this layer ?
      ! Note that the layer can be dynamic or static.
      !
      if(defined(layer))then
        if(p_src%vertDispType == multiLevelFixed)then
          fVertFraction1 = 0.
          do iTmp = 1, fu_NbrOfLevels(p_src%vertLevs)
            fVertFraction1 = fVertFraction1 + p_src%levFraction(iTmp) * &
                                    & fu_vert_overlap_fraction(fu_level(p_src%vertLevs,iTmp), layer)
          end do
          if(fVertFraction1 .eps. 0.)then
            call free_work_array(fCells)
            return  ! No emission into this output layer
          endif
          fVertFraction2 = fVertFraction1
        else if (p_src%vertDispType == singleLevelDynamic) then
          fVertFraction1 = fu_vert_overlap_fraction(p_src%params(iSlot)%layerDynamic, layer)
          fVertFraction2 = fu_vert_overlap_fraction(p_src%params(iSlot+1)%layerDynamic, layer)
          if((fVertFraction2 .eps. 0.) .and. (fVertFraction1 .eps. 0.)) cycle
        else if (p_src%vertDispType == plumeRise) then
          call set_error('Plume rise not supported.','total_from_p_src_descr_unit')
        else
          call set_error('Undefined vertical distripution type','total_from_p_src_descr_unit')
        endif
      else
        fVertFraction1 = 1.0
        fVertFraction2 = 1.0
      endif

      !
      ! Compute the time overlap
      !
      if(p_src%params(iSlot)%time < start)then ! Skip part of slot if needed
        start_integr = start
      else
        start_integr = p_src%params(iSlot)%time
      endif
      if(p_src%params(iSlot+1)%time > start+duration)then ! Skip part of the slot if needed
        end_integr = start+duration
      else
        end_integr = p_src%params(iSlot+1)%time
      endif

      !
      ! Check the rate / composition variation
      !
      ifSlotRateVary = p_src%ifRateVary .or. .not. (fVertFraction1 .eps. fVertFraction2)
        
       ! What is the timezone offset
      if ((p_src%tz_index == SolarTimeIndex) .or. &
                  &p_src%tz_index == LocalTimeIndex) then 
                  ! Use solar time here: it is not an actual emission
            myTZoffset =  fu_set_interval_sec( 3600. * (modulo(p_Src%lon + 180.,360.)-180.) / 15. )
      else
            myTZoffset = fu_set_interval_sec( 1.0*tz_offset(p_src%tz_index))
      endif

      !--------------------------------------------------------------------------------
      !
      ! Here we can finally decide on the integration procedure - along one of the four
      ! options outlined above in the function description. See notebook 11, p.24.
      !
      if(ifSlotRateVary .and. p_src%ifUseTimeVarCoef)then  !-------------------- all dynamic
        !
        ! The most-difficult case. Everything varies, have to go cell-by-cell and hour-by-hour
        ! (or day-by-day, if hourly variation is void). Case 3 of the above function description
        !
        call integrate_dynamic_slot_descr(p_src%params, &    
                                        & start_integr, end_integr, iSlot, &
                                        & p_src%cocktail_descr_lst, p_src%nDescriptors, &
                                        & fCells, &
                                        & fVertFraction1, fVertFraction2, &
                                        & p_src%ifRateVary, &
                                        & amounts)

        fSlot_sec = fu_sec(p_src%params(iSlot+1)%time - p_src%params(iSlot)%time)
        fTime_sec_start = fu_sec(start_integr - p_src%params(iSlot)%time)
        fWeightPast = 1. - fTime_sec_start/fSlot_sec
        !
        ! First, cover the start period until the first round hour, then go ahead with hourly step
        !
        startTmp = start_integr
        endTmp = fu_round_closest_hour(startTmp, forwards)
        if(endTmp > end_integr) endTmp = end_integr

        fTime_sec = fu_sec(endTmp - startTmp)

        ! In some cases, start is the round hour, so can skip this step
        !
        do while(startTmp < end_integr)

          fWeightPast = fWeightPast + 0.5*fTime_sec/fSlot_sec
        
          timeTmpLocal = startTmp + myTZoffset

          do iDescr = 1, p_src%params(iSlot)%nDescriptors
            amounts(iDescr) = amounts(iDescr) + &
                                 & fTime_sec * &                ! time period
                                 & (fWeightPast*fVertFraction1 + (1-fWeightPast)*fVertFraction2) * & !vert
                                 & (fWeightPast * p_src%params(iSlot)%rate_descr_unit(iDescr) + &    ! descr rate
                                     & (1.-fWeightPast) * p_src%params(iSlot)%rate_descr_unit(iDescr)) * &
                                 & p_src%indHour(iDescr, fu_hour(timeTmpLocal)+1) * &
                                 & p_src%indDay(iDescr, fu_weekday(timeTmpLocal)) * &
                                 & p_src%indMonth(iDescr, fu_mon(timeTmpLocal))
          end do  ! iDescr
          !
          ! The time step is completed. The next one is one hour or shorter if this is the end of integration
          !
          fWeightPast = fWeightPast + 0.5*fTime_sec/fSlot_sec
          startTmp = endTmp

          fTime_sec = 3600.
          endTmp = startTmp + one_hour
          if(endTmp > end_integr)then
            endTmp = end_integr
            fTime_sec = fu_sec(endTmp - startTmp)
          endif

        end do  ! startTmp < end_integr


      elseif(ifSlotRateVary)then  !-------------------------------------------- time variation coef void
        !
        ! Hourly/daily/monthly time variation coefficients =1 but slot rate varies.
        ! Can integrate analytically over the whole slot, then multiply with cell-specific
        ! constant rates (case 2 in the above function description). See notebook 11, p.24.
        !
        call integrate_dynamic_slot_descr(p_src%params, &    
                                        & start_integr, end_integr, iSlot, &
                                        & p_src%cocktail_descr_lst, p_src%nDescriptors, &
                                        & fCells, &
                                        & fVertFraction1, fVertFraction2, &
                                        & p_src%ifRateVary, &
                                        & amounts)
        if(error) then 
          call free_work_array(fCells)
          return
        endif

      elseif(p_src%ifUseTimeVarCoef)then !----------------------------------------- slot rates const
        !
        ! Slot-wise rates are constant but time variation coefficients are non-unity. Case 1 in
        ! the function description. Go hour by hour or day by day
        !
        startTmp = start_integr
        endTmp = fu_round_closest_hour(startTmp, forwards)
        if(endTmp > end_integr) endTmp = end_integr

        ! First, sum up over the whole period with time coefficients taken into account
        ! For that, detection of the whether hourly coefficient is useful may bring extra speedup
        !
        timeStep = one_day
        do iDescr = 1, p_src%nDescriptors
          do iTmp = 1, 24
            if(.not. (p_src%indHour(iDescr,iTmp) .eps. 1.0))then
              timeStep = one_hour
              exit
            endif
          end do
        end do

        fTime_sec = fu_sec(endTmp - startTmp)
        fWeightPast = 0.0
        fSlot_sec = fu_sec(p_src%params(iSlot+1)%time - p_src%params(iSlot)%time)
        !
        ! In some cases, start is the round hour, so can skip this step
        !
        do while(startTmp < end_integr)

          fWeightPast = fWeightPast + 0.5*fTime_sec/fSlot_sec

          timeTmpLocal = startTmp + myTZoffset

          do iDescr = 1, p_src%nDescriptors
            amounts(iDescr) = amounts(iDescr) + &
                                 & fTime_sec * &                ! time period
                                 & fVertFraction1 * &           ! vertical overlap
                                 & p_src%params(iSlot)%rate_descr_unit(iDescr) * & ! descr.rate
                                 & p_src%indHour(iDescr,fu_hour(timeTmpLocal)+1) * &
                                 & p_src%indDay(iDescr,fu_weekday(timeTmpLocal)) * &
                                 & p_src%indMonth(iDescr,fu_mon(timeTmpLocal))
          end do  ! iDescr

          !
          ! The time step is completed. The next one is one hour or
          ! shorter if this is the end of integration
          !
          fWeightPast = fWeightPast + 0.5*fTime_sec/fSlot_sec
          startTmp = endTmp

          fTime_sec = fu_sec(timeStep)
          endTmp = startTmp + timeStep
          if(endTmp > end_integr)then
            endTmp = end_integr
            fTime_sec = fu_sec(endTmp - startTmp)
          endif

        end do  ! startTmp < end_integr

      else  !----------------------------------------------------------------------- all static
        !
        ! All is static, i.e. time-independent. The fastest option (case 4 above)
        !
        fTime_sec = fu_sec(end_integr - start_integr)

        do iDescr = 1, p_src%nDescriptors
          amounts(iDescr) = amounts(iDescr) + &
                             & fTime_sec * &            ! time period
                             & fVertFraction1 * &       ! vertical overlap
                             & p_src%params(iSlot)%rate_descr_unit(iDescr) ! descr.rate * sp.fraction, species unit
        end do
      endif  ! four types of temporal dependence

    end do  ! time slots

    call free_work_array(fCells)

  end subroutine total_from_p_src_descr_unit


  !*****************************************************************

  subroutine total_from_p_src_species_unit(p_src, &                      ! mandatory, input
                                         & species, nSpecies, amounts, & ! mandatory, output
                                         & start_, duration_, layer_)    ! optional, input
    !
    ! Returns the amount of the released material IN SPECIES UNIT starting from 
    ! start during the duration time interval. Species must be initialised by that moment
    !
    ! We have to integrate the release rate, which is given in descriptor basic units - for each
    ! descriptor. Species may have different units, so the return value is an array of species
    !
    ! Four cases are considered (notebook 11 pp.21-23)
    ! 1. slot rates are constant: descriptor rate, species fractions, and vertical fraction
    ! 2. cell-connected time variation coefficients are constant (fHour, fDay, fMonth)
    ! 3. both sets are varying
    ! 4. None of them is time dependent
    !
    ! This subroutine produces a total species-wise emission and returns a vector.
    ! If requested, it can skip the cell-wise spatial integration leaving the space for
    ! further map- or total- aiming integration over cells.
    !
    ! Units: SI basic and derived ones
    !
    implicit none

    ! Imported parameters with intent IN
    type(silam_point_source), intent(in) :: p_src
    type(silam_species), dimension(:), pointer :: species
    integer, intent(out) :: nSpecies
    real, dimension(:), intent(out) :: amounts
    type(silja_time), intent(in), optional :: start_
    type(silja_interval), intent(in), optional :: duration_
    type(silja_level), intent(in), optional :: layer_  ! if defined, include only fraction emitted in it

    ! Local variables
    integer :: iDescr, iSpecies, nSpeciesDescr
    type(silam_species), dimension(:), pointer :: speciesDescr
    type(silja_time) :: start
    type(silja_interval) :: duration
    real, dimension(max_descriptors_in_source) :: amountsDescr
    type(silja_level) :: layer
    type(chemical_adaptor) :: adaptor

    nSpecies = 0

    !
    ! Time period for integration, then start and end slots
    !
    if(present(start_))then
      start = start_
      if(start > p_src%params(size(p_src%params))%time)then
        call report(p_src)
        call msg_warning('Source is earlier then the computation starts','total_from_p_src_species_unit')
        return
      endif
    else
      start = p_src%params(1)%time
    endif
    if(present(duration_))then
      duration = duration_
      if(start+duration > p_src%params(size(p_src%params))%time) &
                                      & duration = p_src%params(size(p_src%params))%time - start
      if(start+duration < p_src%params(1)%time)then
        call report(p_src)
        call msg_warning('Source is later then the computation ends','total_from_p_src_species_unit')
        return
      endif
    else
      duration = p_src%params(size(p_src%params))%time - start
    endif
    if (present(layer_)) then
      layer = layer_
    else
      layer = level_missing
    end if

    !
    ! Get the total amounts in descriptor unit
    !
    call total_from_p_src_descr_unit(p_src, &
                                   & amountsDescr, &
                                   & start, duration, layer)
    if(error)return

!do iDescr = 1, p_src%nDescriptors
!call msg('Emission for descriptor: ' + fu_name(p_src%cocktail_descr_lst(iDescr)) + ':',amountsDescr(iDescr))
!end do

    if(sum(amountsDescr(1:p_src%nDescriptors)) == 0.) then 
            return
    endif


    !
    ! Conversion into the species units is comparatively straightforward.
    !
    nullify(species)
    nSpecies = 0
    call add_inventory_p_src(p_src, species, nSpecies)
    amounts(1:nSpecies) = 0.
    !
    ! Now explore the descriptors
    !
    do iDescr = 1, p_src%nDescriptors
      call get_inventory(p_src%cocktail_descr_lst(iDescr), speciesDescr, nSpeciesDescr)
      call create_adaptor(speciesDescr, species, adaptor)
      if (error) return
      do iSpecies = 1, nSpeciesDescr
!call msg('Factor to species unit for:'+fu_name(p_src%cocktail_descr_lst(iDescr)), p_src%fDescr2SpeciesUnit(iSpecies,iDescr))
        amounts(adaptor%iSp(iSpecies)) = amounts(adaptor%iSp(iSpecies)) + &
                                & amountsDescr(iDescr) * p_src%fDescr2SpeciesUnit(iSpecies,iDescr)
!call msg('Amounts for species:' + fu_name(fu_material(species(adaptor%iSp(iSpecies)))),amounts(adaptor%iSp(iSpecies)))
      end do
    end do

  end subroutine total_from_p_src_species_unit


  !*****************************************************************

  integer function fu_source_id_nbr_of_p_src(p_src)
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
    type(silam_point_source), intent(in) :: p_src

    ! Stupidity check
    if(.not.(p_src%defined == silja_false))then
      call set_error('Undefined source given','fu_source_id_nbr_of_p_src')
      return
    endif
    fu_source_id_nbr_of_p_src = p_src%id_nbr

  end function fu_source_id_nbr_of_p_src


  !*****************************************************************

  integer function fu_source_nbr_of_p_src(p_src)
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
    type(silam_point_source), intent(in) :: p_src

    ! Stupidity check
    if(.not.(p_src%defined == silja_false))then
      fu_source_nbr_of_p_src = p_src%src_nbr
    else
      fu_source_nbr_of_p_src = int_missing
      call set_error('Undefined source given','fu_source_nbr_of_point_source')
      return
    endif

  end function fu_source_nbr_of_p_src


  !*****************************************************************

  function  fu_name_p_src(p_src) result(nm)
    implicit none

    character(len=clen) :: nm

    type(silam_point_source), intent(in) :: p_src

    nm = p_src%src_nm

  end function fu_name_p_src


  !*****************************************************************

  function  fu_sector_p_src(p_src) result(nm)
    implicit none

    character(len=clen) :: nm

    type(silam_point_source), intent(in) :: p_src

    nm = p_src%src_nm

  end function fu_sector_p_src


  ! ***************************************************************


  subroutine store_p_src_as_namelist(p_src, uOut, ifOriginalGrd)
    !
    ! Writes one area source term v.3 to an external file. Can do it either in 
    ! original or in dispersion grid and vertical. Will expand all val values
    ! using the emission cocktail descriptors, etc.
    !
    ! All units: SI, unless otherwise is stated
    !
    IMPLICIT NONE
    !
    ! Imported parameters
    TYPE(silam_point_source), intent(in) :: p_src
    integer, intent(in) :: uOut
    logical, intent(in) :: ifOriginalGrd

    ! Local variables
    integer :: iTmp, iDescr

    write(uOut,fmt='(A)')'POINT_SOURCE_5'
    ! Basic parameters: Source name, sector and grid
        !
    write(uOut,fmt='(2A)')'source_name = ', trim(p_src%src_nm)
    write(uOut,fmt='(2A)')'source_sector_name = ', trim(p_src%sector_nm)

    write(uOut,fmt='(2A)')'source_timezone = ', trim(p_src%tz_name)
    write(uOut,fmt='(A,F15.7)')'source_longitude = ', p_src%lon
    write(uOut,fmt='(A,F15.7)')'source_latitude= ', p_src%lat
    write(uOut,fmt='(A,F15.7)')'stack_height= ', p_src%stackHeight

    !
    ! Release rate parameters
    ! We will write the whole source in one file, the whole cocktail at once
    !
    write(uOut,fmt='(2A)')'release_rate_unit = DESCRIPTOR_DEFAULT'

    !
    ! Time slot and variation parameters
    !
    do iTmp = 1, size(p_src%params)
      call store_time_param_as_namelist(p_src%params(iTmp), p_src%cocktail_descr_lst, uOut, point_source)
    end do
    do iDescr = 1, p_src%nDescriptors
      write(uOut,fmt='(2A,1x,24(F10.6,1x))')'hour_in_day_index = ', &
                               & trim(fu_name(p_src%cocktail_descr_lst(iDescr))),(p_src%indHour(iDescr,iTmp),iTmp=1,24)
      write(uOut,fmt='(2A,1x,7(F10.6,1x))')'day_in_week_index = ', &
                               & trim(fu_name(p_src%cocktail_descr_lst(iDescr))),(p_src%indDay(iDescr,iTmp),iTmp=1,7)
      write(uOut,fmt='(2A,1x,12(F10.6,1x))')'month_in_year_index = ', &
                               & trim(fu_name(p_src%cocktail_descr_lst(iDescr))),(p_src%indMonth(iDescr,iTmp),iTmp=1,12)
    end do
    !
    ! If  multi-level structure is constant in time, write it
    !
    if(p_src%vertDispType == multiLevelFixed)then
      write(uOut,fmt='(A)')'vertical_distribution = MULTI_LEVEL_FIXED'
      if(ifOriginalGrd)then
        call report_as_namelist(uOut, p_src%vertLevs, p_src%levFraction, .true.) ! if skip 0
      else
        call report_as_namelist(uOut, p_src%vertLevsDispVert, p_src%levFractDispVert, .true.)
      end if  ! original grid
    else if (p_src%vertDispType == singleLevelDynamic) then
      write(uOut,fmt='(A)')'vertical_distribution = SINGLE_LEVEL_DYNAMIC'
    else if (p_src%vertDispType == plumeRise) then
       write(uOut,fmt='(A)')'vertical_distribution = PLUME_RISE'
    else
      ! not possible
    endif  ! if multi level time-fixed emission
    
    write(uOut,fmt='(A)')'END_POINT_SOURCE_5'
    write(uOut,fmt=*)

  end subroutine store_p_src_as_namelist


  !****************************************************************

  subroutine source_2_map_point_source(p_src, dataPtr, id, ifWholeVertical)
    !
    ! Projects the point source to the map, having the given id as a 
    ! template: the id deterines the grid, level and substance name.
    ! Since the map does not have any aerosol stuff inside - the modes are either
    ! sumed-up or a given one is picked.
    !
    implicit none

    ! Imported parameters
    type(silam_point_source), intent(inout) :: p_src
    real, dimension(:), pointer :: dataPtr
    type(silja_field_id), intent(in) :: id
    logical, intent(in) :: ifWholeVertical
    ! There was iMode, however, it was never called in that form.
    !integer, intent(in), optional :: iMode

    ! Local variables
    real :: xOut, yOut
    integer :: nSpecies, nx, ny, i, iDescr
    integer, dimension(max_species) :: species_index
    real, dimension(max_species) :: amounts
    type(silam_species), dimension(:), pointer :: species
    logical :: ifFound
    !
    ! Stupidity check
    !
    if(.not.(p_src%defined == silja_true))then
      call set_error('Undefined point source given','source_2_map_point_source')
      return
    endif
    if(.not.defined(id))then
      call set_error('Undefined id given','source_2_map_point_source')
      return
    endif
    if(.not.defined(fu_grid(id)))then
      call set_error('Undefined grid given','source_2_map_point_source')
      return
    endif
    if(.not.associated(dataPtr))then
      call set_error('Data array is not associated','source_2_map_point_source')
      return
    endif
    if(fu_number_of_gridpoints(fu_grid(id)) > size(dataPtr))then
      call msg_warning('Too small data array given for the grid of the id')
      call report(id)
      call msg('Size of the data array: ', size(dataPtr))
      call set_error('Too small data array given','source_2_map_point_source')
      return
    endif

    !
    ! Find the grid location of the source. A trick: the point source is always
    ! in geographical co-ordinates, while the output grid can be whatever
    !
    call project_point_to_grid(p_src%lon, p_src%lat, fu_grid(id), xOut, yOut)
    call grid_dimensions(fu_grid(id),nx, ny)

    if(xOut <= 0.5 .or. xOut >= nx+0.5)then
      call msg('Point source is outside the x-limits.x=', xOut)
      call msg_warning('Point source is outside the x-limits.Skipping','source_2_map_point_source')
      return
    endif
    if(yOut <= 0.5 .or. yOut >= ny+0.5)then
      call msg('Point source is outside the y-limits;y=', yOut)
      call msg_warning('Point is outside the y-limits.Skipping','source_2_map_point_source')
      return
    endif
    !
    ! Rate has to be integrated over time. We have three independently varying
    ! parameters - substance mass fraction in the cocktail, aerosol size distribution
    ! (but not the number of modes) and vertical overlap of the layers, all linear in time. 
    !
    ! Linear variation goes between time slots. Let's read the start time
    !
    !
    ! Rate has to be integrated over time. We have such function:
    ! total_from_p_src_species_unit(a_src, start, duration, species, amounts, ifRateOnly)
    !
    nullify(species)
    nSpecies = 0
    call total_from_p_src_species_unit(p_src, &
                                     & species, nSpecies, amounts, &
                                     & fu_accumulation_start_time(id), &
                                     & fu_accumulation_length(id), &
                                     & fu_level(id))
    if(error)return
    if(sum(amounts(1:nSpecies)) < 1e-10)then  ! quite arbitrary number but 1e-10 is indeed small
      return
    endif
    nullify(species)
    nSpecies = 0
    !
    ! Once the full inventory is returned, have to select the specific species for each descriptor.
    ! In fact, this is almost the worst possible solution but we cannot afford keeping all the
    ! maps for all the sources
    !
    ifFound = .false.
    do iDescr = 1, p_src%nDescriptors
      call get_inventory(p_src%cocktail_descr_lst(iDescr), species, nSpecies)
      species_index(iDescr) = fu_index(fu_species(id), species, nSpecies)
      if(species_index(iDescr) >= 1 .and. species_index(iDescr) <= nSpecies) ifFound = .true.
    end do

    if(.not. ifFound) return
    !
    ! Preparatory work is over
    !
    i = int(int(xOut+0.5) + (int(yOut+0.5)-1)*nx + 0.5)
    do iDescr = 1, p_src%nDescriptors
      if(species_index(iDescr) < 1)cycle    ! Some descriptors may not have the needed species
      dataPtr(i) = dataPtr(i) + amounts(species_index(iDescr)) * &
                              & p_src%fDescr2SpeciesUnit(species_index(iDescr),iDescr)
    end do

  end subroutine source_2_map_point_source


  !*****************************************************************************************

  subroutine project_p_src_2_grids(ps, gridMeteo, gridDisp)
    !
    ! Just finds the grid-position in the given grid and stores it to the
    ! position_disp_grd
    !
    implicit none

    ! Imported parameters
    type(silam_point_source), intent(inout) :: ps
    type(silja_grid), intent(in) :: gridMeteo, gridDisp

    ! Local variables
    real :: fX, fY
    integer :: nx, ny

    !
    ! Just create new position from the original one with modified horizontal coordinates
    ! ATTENTION. If the source is emitting to Lagrangian environment, it takes meteo grid
    ! If it emits in Eulerian environment - take dispersion grid.
    !
    if(ps%ifEmissionLagrangian)then
      call project_point_to_grid(ps%lon, ps%lat, gridMeteo, ps%fXDispGrd, ps%fYDispGrd)
      call grid_dimensions(gridMeteo,nx,ny)
    else
      call project_point_to_grid(ps%lon, ps%lat, gridDisp, ps%fXDispGrd, ps%fYDispGrd)
      call grid_dimensions(gridDisp,nx,ny)
    endif

    if(error)then
      ps%fXDispGrd = real_missing
      ps%fYDispGrd = real_missing
      ps%ixDispGrd = int_missing
      ps%iyDispGrd = int_missing
      call set_error('Failed to project the source','project_p_src_2_grids')
    else
      ps%ixDispGrd = nint(ps%fXDispGrd)
      ps%iyDispGrd = nint(ps%fYDispGrd)
    endif
    !
    ! Check that it is inside the grid
    !
    if(ps%ixDispGrd < 0.5 .or. ps%ixDispGrd > nx + 0.5 .or. &
     & ps%iyDispGrd < 0.5 .or. ps%iyDispGrd > ny + 0.5)then
      ps%if_inside_domain = .false.
!      if (smpi_global_tasks > 1) then
!              !Dirty hack for MPI runs
!              !Just turn the source off
!              ps%params(1)%time = really_far_in_future
!
!      else
!              call msg_warning('Source is outside dispersion grid','project_p_src_2_grids')
!              call report(ps)
!              call report(gridDisp)
!              call msg('Grid coordinates are: ixDisp, iyDisp',ps%ixDispGrd, ps%iyDispGrd)
!              call set_error('Source is outside dispersion grid','project_p_src_2_grids')
!              ps%fXDispGrd = real_missing
!              ps%fYDispGrd = real_missing
!              ps%ixDispGrd = int_missing
!              ps%iyDispGrd = int_missing
!              return
!      endif
    else
      ps%if_inside_domain = .true.
    endif

  end subroutine project_p_src_2_grids


  !*******************************************************************

  subroutine prepare_src_vert_params_p_src(p_src, vertMeteo, vertDisp, vertProj)
    !
    ! Projects the vertical parameters of the source to the given vertical
    !
    implicit none

    ! Imported parameter
    type(silam_point_source), intent(inout) :: p_src
    type(silam_vertical), intent(in) :: vertMeteo, vertDisp, vertProj

    ! Local variable
    integer :: indTime, iz, nz, i
    real :: fTopInd, fBottomInd, fTmp, met_data_surf, fTopInd2
    type(silam_vertical) :: vertTmp
    real, dimension(:), pointer :: met_data_column
    type (silja_level) :: levTmp, layer, top, bottom
    character (len=*), parameter :: sub_name="prepare_src_vert_params_p_src"

    !
    ! The major cplit is here: Lagrangian environment works in meteo grid, Eulerian - in dispersion
    !
    if(p_src%ifEmissionLagrangian)then
      p_src%vertLevsDispVert = vertMeteo
    else
      p_src%vertLevsDispVert = vertDisp
    endif
    !
    ! The rest depends on the source vertical arrangements
    !
    if(p_src%vertDispType == multiLevelFixed)then
      !
      ! Nothing to project in parameter strings - all is fixed
      !
      do indTime = 1, size(p_src%params)
        if(allocated(p_src%params(indTime)%levFractDispVert)) deallocate(p_src%params(indTime)%levFractDispVert)
        if(allocated(p_src%params(indTime)%fzDisp)) deallocate(p_src%params(indTime)%fzDisp)
      end do

      if(p_src%ifEmissionLagrangian)then
        !
        ! Lagrangian environment. Layers are kept but stored in terms of the given vertical
        !
        p_src%nzDispVert = fu_NbrOfLevels(p_src%vertLevs)
        allocate(p_src%levFractDispVert(p_src%nzDispVert), &  ! fractions
               & p_src%fzDisp(p_src%nzDispVert + 1), &        ! bottoms and tops => n+1
               & p_src%dz_m(p_src%nzDispVert), stat=i)        ! thickness, m, of the layers
        if(fu_fails(i==0,'Failed to allocate dispersion-vertical lagrangian level fractions', &
                       &sub_name))return
        p_src%levFractDispVert(1:p_src%nzDispVert) = 0.0
        p_src%fzDisp(1:p_src%nzDispVert + 1) = 0.0
        do iz = 1, p_src%nzDispVert
          levTmp = fu_level_to_vertical_crude(fu_level(p_src%vertLevs,iz), p_src%vertLevsDispVert)
          if(error)return
          p_src%levFractDispVert(iz) = p_src%levFraction(iz)
          if(iz==1)p_src%fzDisp(iz) = fu_level_index(fu_lower_boundary_of_layer(levTmp), &
                                                   & p_src%vertLevsDispVert)
          p_src%fzDisp(iz+1) = fu_level_index(fu_upper_boundary_of_layer(levTmp), &
                                            & p_src%vertLevsDispVert)
          p_src%dz_m(iz) = fu_layer_thickness_m(fu_level(p_src%vertLevs,iz))
        end do
      else
        !
        ! Eulerian environment. Re-project the source vertical grid to the given vertical
        ! Since vertical distributions are very poorly known and crude, we do not need
        ! any precise meteo-dependent projection, crude will do the job
        !
        allocate(p_src%levFractDispVert(fu_NbrOfLevels(p_src%vertLevsDispVert)), &
               & p_src%fzDisp(fu_NbrOfLevels(p_src%vertLevsDispVert)), stat=i)
        if(fu_fails(i==0,'Failed to allocate dispersion-vertical level fractions', &
                       &sub_name))return
        p_src%levFractDispVert(:) = 0.0
        p_src%fzDisp(:) = 0.0

        call reproject_verticals(p_src%vertLevs, p_src%levFraction, &   ! vertical from, fractions from
                               & vertProj, p_src%levFractDispVert, &  ! vertical to, fractions to
                               & p_src%fzDisp, p_src%nzDispVert, &  ! mass centres, number of non-zero levels
                               & ifMassCentreInRelUnit=.true.)
        if (error) then
            call msg_warning("error after reproject_verticals", sub_name)
            call report(p_src)
            call set_error("error after reproject_verticals", sub_name)
            return
        endif          
      endif  ! if emission lagrangian

    elseif(p_src%vertDispType == singleLevelDynamic)then
      !
      ! Have to do the job for parameter strings
      !
      nullify(p_src%levFractDispVert,p_src%fzDisp)

      if(p_src%ifEmissionLagrangian)then
        !
        ! In Lagrangian environment the layer is kept at each slot, just recorded in terms
        ! of the given vertical.
        !
        nz = 1
        do indTime = 1, size(p_src%params)
          allocate(p_src%params(indTime)%levFractDispVert(1), p_src%params(indTime)%fzDisp(2))
          p_src%params(indTime)%fzDisp(1) = fu_project_level_crude( &
                           & fu_lower_boundary_of_layer(p_src%params(indTime)%layerDynamic), &
                           & p_src%vertLevsDispVert)
          p_src%params(indTime)%fzDisp(2) = fu_project_level_crude( &
                           & fu_upper_boundary_of_layer(p_src%params(indTime)%layerDynamic), &
                           & p_src%vertLevsDispVert)
          p_src%params(indTime)%levFractDispVert(1) = 1.0
          p_src%params(indTime)%z_size = fu_layer_thickness_m(p_src%params(indTime)%layerDynamic)
        end do  ! indTime

      else
        !
        ! In Eulerian advection, quite detailed reprojection of the dispersion vertical onto 
        ! that single layer is neeed
        !
        nz = fu_NbrOfLevels(p_src%vertLevsDispVert)
        met_data_column => fu_work_array()

        do indTime = 1, size(p_src%params)
          !
          ! Create the temporary vertical from a single layer of the slot and project the 
          ! dispersion_vertical layers onto it
          !
          allocate(p_src%params(indTime)%levFractDispVert(nz), p_src%params(indTime)%fzDisp(nz))
          call set_vertical(p_src%params(indTime)%layerDynamic,vertTmp)
          if(error)return
          call vert_interp_data_crude(p_src%vertLevsDispVert, vertTmp, met_data_column, met_data_surf)
          if (error) return
          do iz = 1, nz
            !
            ! Project the layer of the dispersion vertical to the source vertical 
            !
            layer = fu_level(p_src%vertLevsDispVert,iz)
            fTopInd = fu_project_level(fu_upper_boundary_of_layer(layer), &
                                     & vertTmp, met_data_column, met_data_surf)
            fBottomInd = fu_project_level(fu_lower_boundary_of_layer(layer), vertTmp, &
                                        & met_data_column, met_data_surf)
            if(error)return
            !
            ! Having upper and lower indices of the current dispersion vertical layer with regard to 
            ! original source layer, get the fractions that are within this layer
            ! 
            ! Find out the overlap between the levels (temporary use of xNew variable)
            !
            fTmp = min(fTopInd,1.5) - max(fBottomInd,0.5)

            if(fTmp > 0.0)then
              !
              ! Mass fraction is the overlap multiplied with mass fraction into the original layer
              !
              p_src%params(indTime)%levFractDispVert(iz) = fTmp
              !call msg('ilev, vertical fraction:', iz, ftmp)
              !
              ! Centre of the emitted mass in dispersion vertical requires opposite projection
              ! of source layer to dispersion vertical:
              !
              top = fu_upper_boundary_of_layer(p_src%params(indTime)%layerDynamic)
              bottom = fu_lower_boundary_of_layer(p_src%params(indTime)%layerDynamic)
              fTopInd = fu_project_level_crude(top, vertProj)
              fBottomInd = fu_project_level_crude(bottom, vertProj)

              ! Centre of mass for the emission, in dispersion vertical unit, either m or unitless.
              !
              p_src%params(indTime)%fzDisp(iz) = (0.5 * (max(fBottomInd,real(iz)-0.5) + &
                                                       & min(fTopInd,real(iz)+0.5)) - iz) !* &
!                                                       & fu_layer_thickness_local_unit(layer)
            else
              p_src%params(indTime)%fzDisp(iz) = 0.
              p_src%params(indTime)%levFractDispVert(iz) =0.
            endif  ! fTmp > 0 - overlap

          end do   ! iz in dispersion_vertical
        end do   ! time slots
        call free_work_array(met_data_column)
      endif ! ifEmissionLagrangian

    elseif(p_src%vertDispType == plumeRise)then
      !
      ! For plume-rise, we have to just project the stack height
      !
      if(p_src%ifEmissionLagrangian)then
        p_src%izStackDispVert = fu_project_level_crude(fu_set_constant_height_level(p_src%stackHeight), &
                                                     & p_src%vertLevsDispVert)
      else
        p_src%izStackDispVert = fu_project_level_crude(fu_set_constant_height_level(p_src%stackHeight), &
                                                     & p_src%vertLevsDispVert)
      endif
    endif  ! source vertical arrangements

  end subroutine prepare_src_vert_params_p_src


  !****************************************************************

  logical function fu_useTimeVarCoef_p_src(p_src)
    implicit none
    type(silam_point_source), intent(in) :: p_Src
    fu_useTimeVarCoef_p_src = p_src%ifUseTimeVarCoef
  end function fu_useTimeVarCoef_p_src


  !****************************************************************

  subroutine create_src_cont_grd_p_src(p_src, grid_template, ifVerbose, ifMinimal, ifExtended)
    !
    ! This routine creates the grid similar to the given grid_template
    ! covering all sources of the em_source
    !
    implicit none

    ! Imported parameters
    type(silam_point_source), intent(in) :: p_src
    type(silja_grid), intent(inout) :: grid_template
    logical, intent(in) :: ifVerbose, ifMinimal
    logical, intent(out) :: ifExtended

    ! Local variables
    integer :: iSrc, nx, ny, iCell
    real :: x, y

    if(.not. (p_src%defined == silja_true))then
      call set_error('Undefined source','create_src_cont_grd_p_src')
      return
    endif

    if(defined(grid_template))then
      !
      ! For each source we scan its central points and extending the template_grid 
      ! if these points are outside. For point and bomb sources it is only one point,
      ! while area source has plenty of them.
      !
      if(ifMinimal)then
        !
        ! The requested grid is to be minimal-size just covering the current source
        ! We ignore the current size of the grid and its resolution but still keep the 
        ! projection and grid type
        !
        call project_point_to_grid(p_src%lon, p_src%lat, grid_template, x, y)
        if(error)return
        
        call make_minimal_grid(grid_template, x, y)
        
      else
        !
        ! Just check the source to be inside the grid
        !
        call grid_dimensions(grid_template, nx,ny)

        call project_point_to_grid(p_src%lon, p_src%lat, grid_template, x, y)
        if(x<1 .or. x>nx .or. y<1 .or. y>ny) then
          ifExtended = .true.
          if(ifVerbose)then
            call msg('Grid nx and source x:',nx,x)
            call msg('Grid ny and source y:',ny,y)
            call msg('Extending the grid for the point source:',iSrc)
            call report(p_src)
          endif
          call extend_grid_to_coordinates(grid_template, x, y)
        else
          ifExtended = .false.
        end if
      endif  ! if minimal grid is requested
    else
      !
      ! Grid_template is undefined. Create it using this source as the starting point
      !
      grid_template = fu_set_grid('', lonlat, pole_geographical, &
                                & p_src%lon-0.01, p_src%lat-0.01, &
                                & 3, 3, 0.01, 0.01)   ! nx, ny, dx, dy
      if(error)return
    endif  ! if defined grid_template
  end subroutine create_src_cont_grd_p_src
  
  
  ! *************************************************************************
  
  subroutine prepare_inject_p_src(met_buf)
    !
    ! The subroutine prepares the private module pointers to the fields possibly requested 
    ! for injecting the point sources. As this is called only once before all point sources
    ! their actual needs are unknown. All we can do here is to check the meteo buffer for any 
    ! possibly needed quantity and set the pointer to all existing ones.
    !   
    implicit none

    ! Imported parameters
    type(Tfield_buffer), pointer :: met_buf

    ! Local variables
    integer, dimension(:), pointer :: met_q
    integer :: iQ, iTmp

    ! nullify the pointers 
    nullify(fldTempr)
    nullify(fldWind)
    nullify(fldN2)
    nullify(fldPotTemp)
    nullify(fldAblHeight)
    nullify(fldZ)
    nullify(fldRdown)
    nullify(fldTZindex)
    nullify(fldPressure)
    ind_height = int_missing

    ! Scan the meteo buffer   
    met_q => met_buf%buffer_quantities    
    do iQ = 1, size(met_q)
      if(met_q(iQ) == int_missing)exit
      if(fu_dimension(met_buf, iQ) == 4)then !4D

        select case(met_q(iQ))

          case(temperature_flag)
            fldTempr => met_buf%p4d(iQ)

          case(windspeed_flag)
            fldWind => met_buf%p4d(iQ)
          
          case(brunt_vaisala_freq_flag)
            fldN2 => met_buf%p4d(iQ)          
          
          case(potential_temperature_flag)
            fldPotTemp => met_buf%p4d(iQ)
            
          case(R_down_meteo_flag)
            fldRdown => met_buf%p4d(iQ)

          case(pressure_flag)
            fldPressure => met_buf%p4d(iQ)

          case(height_flag)
            ind_height = iQ
            fldZ => met_buf%p4d(iQ)

          case default
            cycle
            
        end select
      else !2D

        select case(met_q(iQ))
          
          case(abl_height_m_flag)
            fldAblHeight => met_buf%p2d(iQ)
          
          case(timezone_index_flag)
            fldTZindex => met_buf%p2d(iQ)
          
          case default
            cycle
            
        end select
      endif
    enddo
  end subroutine prepare_inject_p_src
  
  
  !**************************************************************************

  subroutine calc_plume_rise(p_src, met_buf, exhTemperature, exhSpeed, stackSize, &
                             & interpCoefMeteo2DispHoriz, &
                             & ifMetHorizInterp, interpCoefMeteo2DispVert, ifMetVertInterp, &
                             & plumeBottom, plumeTop)
    implicit none
    TYPE(silam_point_source), INTENT(in) :: p_src
    TYPE(Tfield_buffer), POINTER :: met_buf
    real, intent(in) :: exhTemperature, exhSpeed, stackSize
    type(TVertInterpStruct), intent(in) ::  interpCoefMeteo2DispVert
    type(THorizInterpStruct), intent(in) ::  interpCoefMeteo2DispHoriz
    logical, intent(in) :: ifMetVertInterp
    logical, intent(in) :: ifMetHorizInterp
    real, intent(inout) :: plumeBottom, plumeTop
    
    ! Local variables
    integer :: stackLevel, iTmp
    real :: ablHeight, temperature, windSpeed
    real :: bruntVaisalaFreq, potentialTemp
    real :: stackHeight, potTempDif, buoyancyFlux
    real :: deltaH, tmp, fSplit, rSemiT, p, effectiveHeight

    ! For debug output
    logical :: debug_output
    integer, save :: plume_rise_funit = int_missing
    debug_output = .true.
    
    ! If necessary, open the debug file 
    if (plume_rise_funit == int_missing .and. debug_output) then  
        plume_rise_funit = fu_next_free_unit()
        open(plume_rise_funit, file = 'plume_rise.csv', iostat = iTmp)
        if(fu_fails(iTmp == 0,'Failed to open plume_rise.csv file','calc_plume_rise'))return

        write(plume_rise_funit,"(25A)")"Name;","Year;", "Month;","Day;","Hour;","Minute;","Sec;", &
                                     & "Temperature;","Exhaust_temp;","Stack_diameter;","Grav_ac;", &
                                     & "Exhaust_speed;","Wind_speed;","Potential_temp;",&
                                     & "Brunt_Vaisala_freq;","ABL_height;","Height_lower;",&
                                     & "Height_upper;","Buoyancy_flux;","Pot_temp_dif;","-;", &
                                     & "f;","r;","Delta_h:;","Stack_height:;"
    endif

    stackHeight = p_src%stackHeight
    stackLevel = p_src%izStackDispVert
    
    ! Get meteo data
    if (.not.(associated(fldTempr) .and. associated(fldWind) .and. associated(fldN2) .and. &
            & associated(fldPotTemp) .and. associated(fldAblHeight)))then
      call set_error('Not all meteo found','calc_plume_rise')
      return
    endif
    ablHeight = fu_get_value(fldAblHeight, nx_meteo, p_src%ixDispGrd, p_src%iyDispGrd, &
                           & met_buf%weight_past, &
                           & interpCoefMeteo2DispHoriz, ifMetHorizInterp)
    
    temperature = fu_get_value(fldTempr, nx_meteo, p_src%ixDispGrd, p_src%iyDispGrd, &
                             & stackLevel, met_buf%weight_past, interpCoefMeteo2DispHoriz, &
                             & interpCoefMeteo2DispVert, ifMetHorizInterp, ifMetVertInterp)

    windspeed = fu_get_value(fldWind, nx_meteo, p_src%ixDispGrd, p_src%iyDispGrd, &
                             & stackLevel, met_buf%weight_past, interpCoefMeteo2DispHoriz, &
                             & interpCoefMeteo2DispVert, ifMetHorizInterp, ifMetVertInterp)

    bruntVaisalaFreq = fu_get_value(fldN2, nx_meteo, p_src%ixDispGrd, p_src%iyDispGrd, &
                             & stackLevel, met_buf%weight_past, interpCoefMeteo2DispHoriz, &
                             & interpCoefMeteo2DispVert, ifMetHorizInterp, ifMetVertInterp)

    potentialTemp = fu_get_value(fldPotTemp, nx_meteo, p_src%ixDispGrd, p_src%iyDispGrd, &
                             & stackLevel, met_buf%weight_past, interpCoefMeteo2DispHoriz, &
                             & interpCoefMeteo2DispVert, ifMetHorizInterp, ifMetVertInterp)

    
    ! Calculate
    potTempDif = bruntVaisalaFreq * potentialTemp / g   ! bruntVaisalaFreq in silam is already squared
    
    
    buoyancyFlux = g * exhSpeed * stackSize * exhSpeed * &
      (exhTemperature - temperature) / (4.0 * exhTemperature)

    if (potTempDif > 0.01) then
      deltaH = 2.6 * (buoyancyFlux * potentialTemp / &
            & (windSpeed * g * potTempDif))** (1.0/3.0)
    else
      deltaH = 40.0 * (log(1 + buoyancyFlux)**2.0) / &
            & (1.0 + 160.0/buoyancyFlux)**0.5 / windSpeed
    endif

    tmp = ablHeight - stackHeight
    fSplit = 0.0

    if (tmp <= 0.5 * deltaH) then
      fSplit = 0.0
    else if (tmp >= 1.5 * deltaH) then
      fSplit = 1.0
    else
      ! Division by zero should be impossible
      fSplit = tmp / deltaH - 0.5
    endif

    p = 1.0 - fSplit
    rSemiT = 0.6 * deltaH / (pi)**0.5
    
    if (abs(p) < 1e-3 .or. abs(p - 1.0) < 1e-3) then
       effectiveHeight = stackHeight + deltaH 
    else 
       effectiveHeight = ablHeight + rSemiT * (p - 2.0)
    endif

    plumeTop = effectiveHeight + rSemiT

    if (abs(p - 1.0) < 1e-3) then
        plumeBottom = max(ablHeight, stackHeight)
    else
        plumeBottom = effectiveHeight - rSemiT
    endif

    ! Error checking
    if (plumeBottom < 0) then
        plumeBottom = 0
        !call set_error('Lower plume height less than 0', 'calc_heights')
         !return
    endif

    if (plumeBottom > plumeTop) then
        plumeTop = plumeBottom
        !call set_error('Lower plume height is greather than upper plume height', 'calc_heights')
         !return
    endif
    
   if (debug_output) then
      write (plume_rise_funit, "(2A, 6(I12,A), 17(F15.7,A))") p_src%src_nm, ";", 0, ";", 0, ";", 0, ";", &
                                 & 0, ";", 0, ";", 0, ";", &
                                 & temperature,";",exhTemperature,";",stackSize,";", &
                                 & g,";",exhSpeed,";",windSpeed,";",potentialTemp,";",bruntVaisalaFreq,";", &
                                 & ablHeight,";", plumeBottom,";", plumeTop,";", buoyancyFlux,";", &
                                 & potTempDif,";-1.0;", fSplit,";", rSemiT,";", deltaH, ";", stackHeight, ";"
    endif
       
  end subroutine calc_plume_rise


  !**************************************************************************

  subroutine inject_emission_euler_p_src(p_src, &
                                       & mapEmis, mapCoordX, mapCoordY, mapCoordZ, & ! Eulerian
                                       & met_buf, disp_buf, &
                                       & now, timestep, &
                                       & ifSpeciesMoment, &
                                       & fMassInjected, &
                                       & interpCoefMeteo2DispHoriz, ifMetHorizInterp, &
                                       & interpCoefMeteo2DispVert, ifMetVertInterp)
    !
    ! Adds the emission flux to the emission mass map and stores the position
    ! of the injected mass to mapCoord. All maps may contain some data before,
    ! so no overwriting - just adding. Note that for successfull trick one has to 
    ! sum-up momentum, not coordinates. Prior conversion to momentum has been done in
    ! source_general for the whole mass map.
    ! Heavily utilises the precomputed variables: 
    ! - reprojected cells to the dispersion grid
    ! - reprojected vertical to the dispersion vertical
    ! Interpolation goes along these precomputed variables ONLY, no redundancy here
    !
    implicit none

    TYPE(silam_point_source), INTENT(in) :: p_src
    type(Tmass_map), intent(inout) :: mapEmis, mapCoordX, mapCoordY, mapCoordZ
    TYPE(Tfield_buffer), POINTER :: met_buf, disp_buf
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep
    logical, intent(in) :: ifSpeciesMoment
    real(r8k), dimension(:), intent(inout) :: fMassInjected
    type(TVertInterpStruct), intent(in) ::  interpCoefMeteo2DispVert
    type(THorizInterpStruct), intent(in) ::  interpCoefMeteo2DispHoriz
    logical, intent(in) :: ifMetVertInterp
    logical, intent(in) :: ifMetHorizInterp

    ! Local variables
    integer :: iSlot, iLev, ix, iy, iSrc, nLevs, iDescr, iSpecies, &
             & iSlotStart, iSlotEnd, ispeciesEmis, iSpeciesDescr
    real :: fWeightPastSrc, fLevFraction, fTimeScale, fCellTotal, factor, timestep_sec
    real, dimension(max_descriptors_in_source) :: fMassTimeCommon
    real, dimension(3) :: ptrCoord
    logical :: ifFound
    type(silam_species), dimension(:), pointer :: species
    real :: plumeBottom, plumeTop, overlapTop, overlapBottom

    !
    ! Have to be very careful: starting and ending of the source slots can be pretty close
    ! especially in case of adjoint run. The same is true for time variation coefficients,
    ! which can vary sharply. Therefore, the total mass injected to the grid will be
    ! integrated along time.
    !
    nullify(species)
    if(.not. p_src%if_inside_domain)return  ! or outside the whole domain if not MPI
    !
    ! General parameters of the release
    !
    call determine_release_params(p_src, now, timestep, mapEmis%nSpecies, &
                                & met_buf, interpCoefMeteo2DispHoriz, ifMetHorizInterp, &
                                & interpCoefMeteo2DispVert, ifMetVertInterp, &
                                & fMassTimeCommon, plumeTop, plumeBottom, fWeightPastSrc, ifFound, &
                                & iSlotStart, iSlotEnd)
    if(error)return

    if(.not. ifFound) return                   ! if any emission found

    if (iSlotEnd - iSlotStart > 2 .and. p_src%vertDispType == singleLevelDynamic ) then
            call msg("singleLevelDynamic can't emit more than two timeslots at once!")
            call set_error("Pending rewriting...","inject_emission_euler_p_src")
            return
    endif
    !
    ! We do the injection layer by layer
    !
    do iLev = 1, nz_dispersion
      !
      ! First do the check for the overlap: speed-up
      !

!call msg('Lev:',iLev)


      if(p_src%vertDispType == multiLevelFixed)then

        ptrCoord(3) = p_src%fzDisp(iLev)
        fLevFraction = p_src%levFractDispVert(iLev)
        
      elseif (p_src%vertDispType == singleLevelDynamic) then

        ifFound = .false.
        do iSlot = iSlotStart, iSlotEnd-1
          if( .not. ((p_src%params(iSlot)%levFractDispVert(iLev) .eps. 0.0) .and. &
                   & (p_src%params(iSlot+1)%levFractDispVert(iLev) .eps. 0.0)))then
            ifFound = .true.
            exit
          endif
        end do
        if(.not. ifFound) cycle
        ptrCoord(3) = p_src%params(iSlotStart)%fzDisp(iLev) * fWeightPastSrc + &
                    & p_src%params(iSlotEnd)%fzDisp(iLev) * (1. - fWeightPastSrc) 
        fLevFraction = p_src%params(iSlotStart)%levFractDispVert(iLev) * fWeightPastSrc + &
                     & p_src%params(iSlotEnd)%levFractDispVert(iLev) * (1. - fWeightPastSrc)
   
      elseif (p_src%vertDispType == plumeRise) then           !======================= PLUME RISE

        overlapBottom = max(plumeBottom, disp_layer_top_m(iLev-1))
        overlapTop = min(plumeTop, disp_layer_top_m(iLev))

        if (overlapBottom <= overlapTop) then
          ! CM relative to the cell center 
          ptrCoord(3) = (overlapBottom + overLapTop) / 2
          ptrCoord(3) =  (ptrCoord(3) - (disp_layer_top_m(iLev) + disp_layer_top_m(iLev-1))/2) / &
                       &  (disp_layer_top_m(iLev) - disp_layer_top_m(iLev-1))
          if(abs(ptrCoord(3)) > 0.5)then
            call msg('Relative centre of mass position is strange in layer:',iLev, ptrCoord(3))
            call msg('Plume top, bottom:', plumeTop, plumeBottom)
            call msg('Dispersion layer top [m]:', disp_layer_top_m(iLev), disp_layer_top_m(iLev-1))
            call set_error('Wrong calculated vertical plume positon','inject_emission_point_source')
            if(abs(ptrCoord(3)) - 0.5 < 0.0001)then
              call unset_error('inject_emission_point_source')
              ptrCoord(3)  = sign(0.498, ptrCoord(3))
            else
              return
            endif
          endif
          fLevFraction = (overlapTop - overlapBottom) / (plumeTop- plumeBottom)
          !                        &(disp_layer_top_m(iLev) - disp_layer_top_m(iLev-1))
        else
          ptrCoord(3) = 0 !(disp_layer_top_m(iLev) + disp_layer_top_m(iLev-1)) / 2
          fLevFraction = 0.
        endif
          
      else
        call set_error('Vertical distripution type not supported.' ,'inject_emission_point_source')
      endif

      if(fLevFraction .eps. 0.0)cycle  ! nothing for this dispersion layer
      !   call msg('ilev, vertical fraction:', ilev, fLevFraction)
      ptrCoord(1) = p_src%fXDispGrd - p_src%ixDispGrd
      ptrCoord(2) = p_src%fYDispGrd - p_src%iyDispGrd
      fCellTotal = 0.0
      if(error)return

      !
      ! Emit the mass going species by species
      !
      do iDescr = 1, p_src%nDescriptors
        !
        ! Time variation coefficients can be handled this way:
        !

!call msg('iDescr:',iDescr)

        do iSpeciesDescr = 1, p_src%nSpeciesInDescr(iDescr)

!call msg('iSpeciesDescr:',iSpeciesDescr)

          iSpeciesEmis = p_src%pEmisSpeciesMapping(iSpeciesDescr,iDescr)
          factor = p_src%fDescr2SpeciesUnit(iSpeciesDescr,iDescr) * fLevFraction 
          mapEmis%arM(iSpeciesEmis, p_src%id_Nbr, iLev,p_src%ixDispGrd,p_src%iyDispGrd) = & 
                 & mapEmis%arM(iSpeciesEmis, p_src%id_Nbr, iLev, p_src%ixDispGrd, p_src%iyDispGrd) + &
                 & fMassTimeCommon(iDescr) * factor
          fCellTotal = fCellTotal + fMassTimeCommon(iDescr) * factor
          fMassInjected(iSpeciesEmis) = fMassInjected(iSpeciesEmis) + fMassTimeCommon(iDescr) * factor

          if (ifSpeciesMoment) then
            mapCoordX%arm(iSpeciesEmis,p_src%id_nbr,ilev,p_src%ixDispGrd,p_src%iyDispGrd) = &
               & mapCoordX%arm(iSpeciesEmis, p_src%id_nbr, ilev, p_src%ixDispGrd, p_src%iyDispGrd) + &
               & fMassTimeCommon(iDescr) * factor * ptrCoord(1)
            mapCoordY%arm(iSpeciesEmis, p_src%id_nbr, ilev, p_src%ixDispGrd, p_src%iyDispGrd) = &
               & mapCoordY%arm(iSpeciesEmis, p_src%id_nbr, ilev, p_src%ixDispGrd, p_src%iyDispGrd) + &
               & fMassTimeCommon(iDescr) * factor * ptrCoord(2)
            mapCoordZ%arm(iSpeciesEmis, p_src%id_nbr, ilev, p_src%ixDispGrd, p_src%iyDispGrd) = &
               & mapCoordZ%arm(iSpeciesEmis, p_src%id_nbr, ilev, p_src%ixDispGrd, p_src%iyDispGrd) + &
               & fMassTimeCommon(iDescr) * factor * ptrCoord(3)
          end if

        end do  ! species inside the descriptor
      end do  ! iDescr

      ! If we use bulk moment, add it of the injected masses and the total injected mass
      !
      if (.not. ifSpeciesMoment) then
        mapCoordx%arM(1,p_src%id_nbr, iLev, p_src%ixDispGrd,p_src%iyDispGrd) = &
                     & mapCoordx%arM(1,p_src%id_nbr, iLev, p_src%ixDispGrd,p_src%iyDispGrd) + &
                     & ptrCoord(1) * fCellTotal
        mapCoordy%arM(1,p_src%id_nbr, iLev, p_src%ixDispGrd,p_src%iyDispGrd) = &
                     & mapCoordy%arM(1,p_src%id_nbr, iLev, p_src%ixDispGrd,p_src%iyDispGrd) + &
                     & ptrCoord(2) * fCellTotal
        mapCoordz%arM(1,p_src%id_nbr, iLev, p_src%ixDispGrd,p_src%iyDispGrd) = &
                     & mapCoordz%arM(1,p_src%id_nbr, iLev, p_src%ixDispGrd,p_src%iyDispGrd) + &
                     & ptrCoord(3) * fCellTotal
      end if

      mapEmis%ifColumnValid(p_src%id_nbr,p_src%ixDispGrd,p_src%iyDispGrd) = .true.
      mapEmis%ifGridValid(ilev, p_src%id_nbr) = .true.

    end do  ! iLev dispersion


  end subroutine inject_emission_euler_p_src


  !**************************************************************************

  subroutine inject_emission_lagr_p_src(p_src, &
                                      & lpset, arParticleMass, & ! Lagrangian
                                      & ChemRunSetup, &  ! translate emission species to transport
                                      & met_buf, disp_buf, &
                                      & now, timestep, &
                                      & fMassInjected)
    !
    ! Adds the emission flux to Lagrangian structure by starting new particles.
    ! Note that particles fly in the meteorological grid to utilise max of available 
    ! dynamic information and also to have a simpler connection to pressure and omega-wind.
    !
    ! Therefore, all source variables with "dispersion-grid" meaning here mean "meteo grid"
    ! Refer to the project_p_src_2_grids, where the selection is done.
    !
    implicit none

    TYPE(silam_point_source), INTENT(in) :: p_src
    type(Tlagrange_particles_set), INTENT(inout), target :: lpSet
    real, dimension(:), intent (in) :: arParticleMass
    type(TchemicalRunSetup), pointer :: ChemRunSetup
    TYPE(Tfield_buffer), POINTER :: met_buf, disp_buf
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep
    real(r8k), dimension(:), intent(inout) :: fMassInjected

    ! Local variables
    integer :: iSlot, iLev, ix, iy, iSrc, nLevs, iDescr, nSpEmis, iSpDescr, nP, iParticle, iP, &
             & nLevsToStart, iSpEmis, iSpTo, iSpTransp
    real :: fWeightPastSrc, fLevFraction, fTimeScale, fCellTotal, factor, timestep_sec, fDx, fDy
    integer :: nSpeciesTrn
    real, dimension(:), pointer :: xCellSize, yCellSize
    real, dimension(max_levels) :: xSize, ySize, zSize
    real, dimension(max_species) :: fMassTmp
    real, dimension(max_descriptors_in_source) :: fMassTimeCommon
    integer, dimension(max_levels) :: nPartInLev
    logical :: ifFound
    type(TspeciesReference), dimension(:), pointer :: references
    integer :: iSlotStart, iSlotEnd
    real, dimension(max_levels) :: fLevBottom, fLevTop, meteo_heights
    type(TVertInterpStruct), pointer ::  interpCoefVert_void
    type(THorizInterpStruct), pointer ::  interpCoefHoriz_void
    real, dimension(:,:), pointer :: arDyn, arMass
    integer, dimension(:), pointer :: arStatus

    arDyn    => lpSet%lpDyn
    arMass   => lpSet%lpMassTrn
    arStatus => lpSet%lpStatus
    !
    ! Have to be very careful: starting and ending of the source slots can be pretty close
    ! especially in case of adjoint run. The same is true for time variation coefficients,
    ! which can vary sharply. Therefore, the total mass injected to the grid will be
    ! integrated along time.
    !

!call msg('1')
    nSpeciesTrn = lpset%nSpeciesTrn
    
    references => chemRunSetup%refEmis2Transp_mass
    nSpEmis = size(references)
    !
    ! General parameters of the release
    !
    
!call msg('2')

    call determine_release_params(p_src, now, timestep, nSpEmis, &
                                & met_buf, interpCoefHoriz_void, .false., &
                                & interpCoefVert_void, .false., &
                                & fMassTimeCommon, fLevTop(1), fLevBottom(1), fWeightPastSrc, ifFound, &
                                & iSlotStart, iSlotEnd)
    if(error)return

!call msg('3')

    if(.not. ifFound) return                   ! if any emission found
    !
    ! Get the amount of each species to be released
    !
    fMassTmp(1:nSpeciesTrn) = 0.0

!call msg('4')

    do iDescr = 1, p_src%nDescriptors  ! all descriptors
      do iSpDescr = 1, p_src%nSpeciesInDescr(iDescr)    ! all species over descriptor
        iSpEmis = p_src%pEmisSpeciesMapping(iSpDescr,iDescr)
        !
        ! This emission species contributes to nRefSpecies transport species
        !
        do iSpTo = 1, ChemRunSetup%refEmis2Transp_mass(ispEmis)%nRefSpecies  
          factor = ChemRunSetup%refEmis2Transp_mass(ispEmis)%fract(iSpTo)   ! fractionation
          iSpTransp = ChemRunSetup%refEmis2Transp_mass(iSpEmis)%indSpeciesTo(iSpTo)  ! to whom
          fMassTmp(iSpTransp) = fMassTmp(iSpTransp) + &
                 & fMassTimeCommon(iDescr) * p_src%fDescr2SpeciesUnit(iSpDescr,iDescr) * factor
        end do
      end do
    end do  ! descr

!call msg('5')


    !
    ! Lagrangian injection goes particle by particle. At least one particle is
    ! always released, so the particle masses are NOT identical - but an effort is 
    ! made to have them close. Within one release step they are identical.
    !
    ! Mass of a single particle and number of released particles: take the 
    ! low-mass threshold for each emission species, compare it to the rate of
    ! emission for this species, and take the largest ratio as the number of particles.
    ! Other species are distributed among these particles
    !
    nP = 1
    do iSpTransp = 1, nSpeciesTrn
      nP = max(nP, int(fMassTmp(iSpTransp) / arParticleMass(iSpTransp)))
!call msg('fMassTmp(iSpTransp), arParticleMass(iSpTransp)', fMassTmp(iSpTransp), arParticleMass(iSpTransp))
!call msg('ratio:',fMassTmp(iSpTransp) / arParticleMass(iSpTransp), nP)
    end do

    !
    ! Let's get the dispersion cell size
    !
    xCellSize => fu_grid_data(meteo_cell_x_size_fld)  ! temporary use of xSize variable
    yCellSize => fu_grid_data(meteo_cell_y_size_fld)  ! temporary use of ySize variable
    iP = p_src%ixDispGrd + (p_src%iyDispGrd-1) * nx_meteo
    fDx = 1./xCellSize(iP)
    fDy = 1./yCellSize(iP)     ! end of temporary use of variables
    !
    ! Particles are injected layer-by-layer of the p_src (no connection to dispersion layer)
    ! First determine the parameters and layers
    !
    timestep_sec = fu_sec(timestep)
    if(error)return

    if(p_src%vertDispType == multiLevelFixed)then                 !================== MULTI-LAYER FIXED
      !
      ! many levels. Set the variables to be used further
      !
      nLevsToStart = 0
      do iLev = 1, nz_meteo
        if(p_src%levFractDispVert(iLev) .eps. 0.0)cycle ! Disp is meteo here
        nLevsToStart = nLevsToStart + 1

        xSize(nLevsToStart) = p_src%params(iSlotStart)%xy_size * fWeightPastSrc + &            ! size of the term
                            & p_src%params(iSlotEnd)%xy_size * (1. - fWeightPastSrc) + &
               & sqrt( get_kz(fldZ, fldRdown, iLev,p_src%ixDispGrd, p_src%iyDispGrd, fWeightPastSrc)  * &
                                 & abs(timestep_sec) / 2.)
        ySize(nLevsToStart) = xSize(nLevsToStart)
        zSize(nLevsToStart) = p_src%dz_m(iLev)
        fLevBottom(nLevsToStart) = p_src%fzDisp(iLev)   ! relative, in meteo grid
        fLevTop(nLevsToStart) = p_src%fzDisp(iLev+1)    ! relative, in meteo grid
        nPartInLev(nLevsToStart) = nint(nP * p_src%levFractDispVert(iLev))  ! in meteo grid
      end do
      if(nLevsToStart == 0)then
        call set_error('No layers to emit anything','inject_emission_p_src')
        return
      endif
      !
      ! Ensure that the number of particles to start is exact
      !
      do while(sum(nPartInLev(1:nLevsToStart)) /= nP)
        iP = 1
        do iLev = 1, nLevsToStart   ! find the level with the largest number of particles injected
          if(nPartInLev(iLev) > nPartInLev(iP)) iP = iLev
        end do
        nPartInLev(iP) = max(nPartInLev(iP) + (nP - sum(nPartInLev(1:nLevsToStart))), 0)
      enddo
        
    else if (p_src%vertDispType == singleLevelDynamic) then         !================= SINGLE LAYER DYNAMIC

!call msg('7')
      nLevsToStart = 1
      nPartInLev(1) = nP
      fLevTop(1) = p_src%params(iSlotStart)%fzDisp(2) * fWeightPastSrc + &      ! relative coord
                 & p_src%params(iSlotEnd)%fzDisp(2) * (1. - fWeightPastSrc)
      fLevBottom(1) = p_src%params(iSlotStart)%fzDisp(1) * fWeightPastSrc + &   ! relative coord
                    & p_src%params(iSlotEnd)%fzDisp(1) * (1. - fWeightPastSrc) 
      
!call msg('8',iSlotStart,iSlotEnd)
!call msg('fWeightPast, p_src%params(iSlotStart)%xy_size',fWeightPast, p_src%params(iSlotStart)%xy_size)
!call msg('p_src%ixDispGrd, p_src%iyDispGrd',p_src%ixDispGrd, p_src%iyDispGrd)
!call msg('fLevTop(1) + fLevBottom(1)',fLevTop(1), fLevBottom(1))
!print *, 'p_src%params(iSlotStart)%fzDisp', p_src%params(iSlotStart)%fzDisp
!print *, 'p_src%params(iSlotEnd)%fzDisp', p_src%params(iSlotEnd)%fzDisp

      iLev = nint((fLevTop(1) + fLevBottom(1)) * 0.5 + 0.5) ! Level just above middle

      xSize(1) = p_src%params(iSlotStart)%xy_size * fWeightPastSrc + &            ! size of the term
               & p_src%params(iSlotEnd)%xy_size * (1. - fWeightPastSrc) + &
               & sqrt( get_kz(fldZ, fldRdown, iLev,p_src%ixDispGrd, p_src%iyDispGrd, &
                            & met_buf%weight_past) * abs(timestep_sec) / 2.)

      ySize(1) = xSize(1)
      zSize(1) = p_src%params(iSlotStart)%z_size * fWeightPastSrc + &
               & p_src%params(iSlotEnd)%z_size * (1. - fWeightPastSrc)

!call msg('9')

    elseif(p_src%vertDispType == plumeRise) then                  !=================== PLUME RISE
      !
      ! Note that the plume top and bottom are to be converted to relative indices
      !
      nLevsToStart = 1
      zSize(1) = fLevTop(1) - fLevBottom(1)  ! absolute, m
      call column_from_buffer(met_buf, ind_height, p_src%ixDispGrd, p_src%iyDispGrd, nx_meteo, &
                            & meteo_heights, &
                            & interpCoefHoriz_void, interpCoefVert_void, &  ! whatever
                            & .false., .false., met_buf%weight_past)
      if(error)return
      fLevTop(1) = max(0.5, fu_value_index_in_array(fLevTop(1), &   ! m, turn to relative in meteo grid
                                                  & meteo_heights, nz_meteo) - 0.5)
      fLevBottom(1) = max(0.5, fu_value_index_in_array(fLevBottom(1), & ! m, -> relative in meteo grid
                                                  & meteo_heights, nz_meteo) - 0.5)
       xSize(1) = sqrt( &
          & get_kz(fldZ, fldRdown,  p_src%izStackDispVert,p_src%ixDispGrd, p_src%iyDispGrd, &
                 & met_buf%weight_past) * &
!       fu_get_value(fldKz, nx_meteo, p_src%ixDispGrd, p_src%iyDispGrd, &
!                                 & p_src%izStackDispVert, fWeightPast, &
!                                 & interpCoefHoriz_void, interpCoefVert_void, & ! whatever
!                                 & .false., .false.) * &
                     & abs(timestep_sec) / 2.)
      ySize(1) = xSize(1)
      nPartInLev(1) = nP
    else
      call set_error('Vertical distribution type not supported:' + fu_str(p_src%vertDispType), &
                   & 'inject_emission_p_src')
      return
    endif
!------------------------------------------------------------------------
!
! UNcomment this in order to get the particles in the absolute pressure coordinates
!    !
!    ! Having computed fLevBottom(1:nLevsToStart) and fLevBottom(1:nLevsToStart) in relative indices
!    ! of meteo vertical, now we have to turn them into absolute presure
!    !
!    do iLev = 1, nLevsToStart
!      fLevBottom(iLev) = fu_4d_interpolation(fldPressure, &
!                                           & p_src%fxDispGrd, p_src%fyDispGrd, fLevBottom(iLev), &
!                                           & nx_meteo, ny_meteo, nz_meteo, &
!                                           & fWeightPast, &
!                                          & linear, linear, notallowed)
!      fLevTop(iLev) =  fu_4d_interpolation(fldPressure, &
!                                         & p_src%fxDispGrd, p_src%fyDispGrd, fLevTop(iLev), &
!                                         & nx_meteo, ny_meteo, nz_meteo, &
!                                         & fWeightPast, &
!                                         & linear, linear, notallowed)
!    end do
    !
    ! Make-up the particles
    !
!call msg('10')

    iParticle = lpset%iFirstEmptyParticle
    do iLev = 1, nLevsToStart
      do iP = 1, nPartInLev(iLev)
        !
        ! Find free space
        !
!call msg('11, iLev,iP',iLev,iP)

        do while(arStatus(iParticle) /= int_missing)
          iParticle = iParticle + 1
          if(iParticle >= lpset%nop)then
      if(iParticle >= lpset%nop)then
call msg('Enlarging the number of particles, 1:',  lpset%nop + max(lpset%nop*5/4,nPartInLev(iLev)-iP))
          call enlarge_lagrange_particles_set(lpset, lpset%nop + max(lpset%nop*5/4,nPartInLev(iLev)-iP))
      endif
            exit
          endif
        end do
        !
        ! Put the masses
        ! Emit the mass going species by species. Note multi-step re-indexing: from species index in 
        ! descriptor to species index in emission list, then to species index in transport list.
        ! Simialrly, MassTimeCommon is mass of a descriptor cocktail, which needs to be explored
        !
        arMass(1:nSpeciesTrn, iParticle) = fMassTmp(1:nSpeciesTrn) / real(nP)
!call msg('12')


        !
        ! Position, turbulent movement, starting size.
        ! Note that vertical has to be converted to relative
        !
        arDyn(lp_x, iParticle) = fu_random_number_center(p_src%fxDispGrd, max(0.1,xSize(iLev)) * fDx)
        arDyn(lp_y, iParticle) = fu_random_number_center(p_src%fyDispGrd, max(0.1,ySize(iLev)) * fDy)
        arDyn(lp_z, iParticle) = fu_random_number_boundaries(fLevBottom(iLev), fLevTop(iLev))
        
!call msg('13')

        arDyn(lp_uT:lp_wT, iParticle) = 0.0    ! turbulent-wind motion
        arDyn(lp_dx, iParticle) = xSize(iLev)    ! metres
        arDyn(lp_dy, iParticle) = ySize(iLev)    ! metres

        !!!!Zero-sized particels cause troubles
        arDyn(lp_dz, iParticle) = max(zSize(iLev), 100.)   ! 100 meters -- minimum size 
        arStatus(iParticle) = p_src%id_nbr

        lpset%nNewPart = lpset%nNewPart + 1
        lpset%newPartList(lpset%nNewPart) = iParticle

        iParticle = iParticle + 1

!call msg('14')

        if(iParticle >= lpset%nop)then
                    call msg('Enlarging lpset, 2:',  lpset%nop + max(lpset%nop*5/4,nPartInLev(iLev)-iP))
          call enlarge_lagrange_particles_set(lpset, lpset%nop + max(lpset%nop*5/4,nPartInLev(iLev)-iP))
          exit
        endif
      end do  ! iP
    end do  ! iLev of the p_src

!call msg('15')

    !
    ! Fix the starting particle. Not exactly (place can be occupied) but better than starting from 1
    !
    lpset%iFirstEmptyParticle = iParticle

  contains
                                      
    real function get_kz(Z, Rdown, iLev, ix, iy, fweightpast)
      integer, intent(in) :: iLev, ix, iy
      real, intent(in) :: fweightpast
      type(field_4d_data_ptr), pointer :: Z, Rdown

      integer :: iCellIndex
      real :: dz

      iCellIndex = ix+(iy-1)*nx_meteo

      if (iLev==1) then 
        dz = (Z%past%p2d(iLev)%ptr(iCellIndex) * fweightpast + &
               & Z%future%p2d(iLev)%ptr(iCellIndex) * (1.-fweightpast))
      else 
        dz = (  (Z%past%p2d(iLev)%ptr(iCellIndex)- Z%past%p2d(iLev-1)%ptr(iCellIndex)) * fweightpast + &
              & Z%future%p2d(iLev)%ptr(iCellIndex) - Z%future%p2d(iLev-1)%ptr(iCellIndex) * (1.-fweightpast))
      endif

      get_kz =  dz / (Rdown%past%p2d(iLev)%ptr(iCellIndex) * fweightpast + &
                                  & Rdown%future%p2d(iLev)%ptr(iCellIndex) * (1.-fweightpast)) 
    end function get_kz

  end subroutine inject_emission_lagr_p_src


  !************************************************************************
  
  subroutine determine_release_params(p_src, now, timestep, nSpecies, &
                                    & met_buf, interpCoefMeteo2DispHoriz, ifMetHorizInterp, &
                                    & interpCoefMeteo2DispVert, ifMetVertInterp, &
                                    & fMassTimeCommon, fLevTop, fLevBottom, fWeightPastSrc, ifFound, &
                                    & iSlotStart, iSlotEnd)
    !
    ! Determines the main release parameters: spatial profile and absolute release of 
    ! each species
    !
    implicit none
    
    ! Imported parameters
    TYPE(silam_point_source), INTENT(in) :: p_src
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep
    type(silja_time) :: local_now
    integer, intent(in) :: nSpecies
    TYPE(Tfield_buffer), POINTER :: met_buf
    type(TVertInterpStruct), intent(in) ::  interpCoefMeteo2DispVert
    type(THorizInterpStruct), intent(in) ::  interpCoefMeteo2DispHoriz
    logical, intent(in) :: ifMetVertInterp
    logical, intent(in) :: ifMetHorizInterp
    real, dimension(:), intent(out) :: fMassTimeCommon
    real, intent(out) :: fLevTop, fLevBottom, fWeightPastSrc
    logical, intent(out) :: ifFound
    integer, intent(out) :: iSlotStart, iSlotEnd

    ! Local variables
    integer :: iSlot, iLev, nLevs, iDescr, iTmp, iMeteo
    real :: fLevFraction, fTimeScale, fCellTotal, factor, timestep_sec, exhTempr, exhSpeed, stackSize
    type(silja_time) :: timeStart, timeEnd

!call msg('in')

    timeStart = now
    timeEnd = timeStart + fu_abs(timestep)
    ifFound = .false.
    
    ! Check that this source is active during this time period and cut the start-end time if needed
    !
    if(timeEnd <= p_src%params(1)%time .or. timeStart >= p_src%params(size(p_src%params))%time) return

    if(timeStart < p_src%params(1)%time) timeStart = p_src%params(1)%time
    if(timeEnd > p_src%params(size(p_src%params))%time) timeEnd = p_src%params(size(p_src%params))%time
    if(error)return
    !
    ! Find the range of slots to be looked at
    !
    call get_overlapping_slots(p_src%params, timeStart, timeEnd-timeStart, iSlotStart, iSlotEnd)
    if(error .or. iSlotStart >= iSlotEnd)return  ! nothing, for whatever reasons

    fWeightPastSrc = (p_src%params(iSlotEnd)%time - now) / &
                   & (p_src%params(iSlotEnd)%time - p_src%params(iSlotStart)%time)
    !
    ! Rate has to be integrated over time. We have such function:
    ! total_from_p_src_species_unit(a_src, start, duration, layer, ifRateOnly, &
    !                                 species, nSpecies, amounts, z_moment)
    !
    call total_from_p_src_descr_unit(p_src, &
                                   & fMassTimeCommon, &
                                   & timeStart, timeEnd - timeStart) !, &
!                                     & fu_level(dispersion_vertical,iLev))
    !
    ! If zero emission, nothing to emit...
    !
    if(all(fMassTimeCommon(1:nSpecies) <= 0.0))then
      ifFound = .false.
      return  ! nothing, for whatever reasons
    endif
    ifFound = .true.
    !
    ! The above function ignores the temporal variation. Have to include it here
    !
    
    ! What time is now?
    if (p_src%tz_index == SolarTimeIndex) then
        local_now = now + fu_set_interval_sec( 3600. * (modulo(p_Src%lon + 180.,360.)-180.) / 15. )
    elseif  (p_src%tz_index == LocalTimeIndex) then
        !if local timezone -- find proper timezone index
        ! ideally, should be done once per source...
        if (.not.(associated(fldTZindex)))then
                call set_error('No tz index field found','determine_release_params')
                return
        endif

        iMeteo = fu_grid_index(nx_meteo, p_src%ixDispGrd, p_src%iyDispGrd, interpCoefMeteo2DispHoriz)
        iTmp = nint(fldTZindex%present%ptr(iMeteo))

!        = nint(fu_get_value(fldTZindex, nx_meteo, p_src%ixDispGrd, p_src%iyDispGrd, &
!                               & met_buf%weight_past, &
!                               & interpCoefMeteo2DispHoriz, ifMetHorizInterp))
        if (iTmp > number_of_time_zones .or. iTmp < 0) then
                call msg("Got wrong timezone index for the source", iTmp)
                call set_error("Can't get local timezone index for the source",'determine_release_params')
                return
        endif 

        local_now = now + fu_set_interval_sec(1.0*tz_offset(iTmp))
        call msg("Timezone index for the source", iTmp)
    else
        local_now = now + fu_set_interval_sec(1.0*tz_offset(p_src%tz_index))
    endif

    if(p_src%ifUseTimeVarCoef)then
      do iDescr = 1, p_src%nDescriptors
        fMassTimeCommon(iDescr) = fMassTimeCommon(iDescr) * &
                                            & p_src%indHour(iDescr, fu_hour(local_now)+1) * &
                                            & p_src%indDay(iDescr, fu_weekday(local_now)) * &
                                            & p_src%indMonth(iDescr, fu_mon(local_now))
      enddo
    endif
    timestep_sec = fu_sec(timestep)

    ! Three types of vertical structure: fixed multi-;ayer, dynamic single-layer, 
    ! and plume-rise single layer. Note that static structure does not have to be projected
    !
    if (p_src%vertDispType == plumeRise) then
      exhTempr = p_src%params(iSlotStart)%tempr * fWeightPastSrc + &
                          & p_src%params(iSlotEnd)%tempr * (1. - fWeightPastSrc) 
      
      exhSpeed = p_src%params(iSlotStart)%z_velocity * fWeightPastSrc + &
                          & p_src%params(iSlotEnd)%z_velocity * (1. - fWeightPastSrc) 
      
      stackSize = p_src%params(iSlotStart)%xy_size * fWeightPastSrc + &
                          & p_src%params(iSlotEnd)%xy_size * (1. - fWeightPastSrc)

      call calc_plume_rise(p_src, met_buf, exhTempr, exhSpeed, stackSize, &
                         & interpCoefMeteo2DispHoriz, ifMetHorizInterp, &
                         & interpCoefMeteo2DispVert, ifMetVertInterp, &
                         & fLevBottom, fLevTop)
    else
      fLevBottom = real_missing
      fLevTop = real_missing
    endif

!call msg('out')
  end subroutine determine_release_params


  ! ***************************************************************

  SUBROUTINE report_point_source(ps)

    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silam_point_source), INTENT(in) :: ps
    
    integer :: j
    real, dimension(:), pointer :: fTmp

    ! Local declarations:
    call msg('================= POINT SOURCE =======================')
    call msg(fu_connect_strings('Point source:',ps%src_nm))
    call msg('start lon & lat: ', ps%lon, ps%lat)
    call msg(fu_connect_strings('start time of release: ', &
                              & fu_str(ps%params(1)%time)))
    call msg(fu_connect_strings('duration of release: ', &
                              & fu_str(ps%params(size(ps%params))%time-ps%params(1)%time)))
    call msg(fu_connect_strings('end time of release: ',&
                              & fu_str(ps%params(size(ps%params))%time)))
    call msg(' Parameters at the beginning of the release: ')
    if(ps%vertDispType == multiLevelFixed)then
      call msg('Source vertical:')
      call report(ps%vertLevs, .true.)
    else if (ps%vertDispType == singleLevelDynamic) then
      IF (ps%params(1)%vertical_distribution == vertically_single_level) THEN
        call msg('single start level: ')
        CALL report(ps%params(1)%layerDynamic)
      ELSE
        call msg('bottom & top: ')
        CALL report(ps%params(1)%layerDynamic)
      END IF
    else if (ps%vertDispType == plumeRise) then
      call msg('Vertical Distripution: Plume rise')
      call msg('stack height: ', ps%stackHeight)
    else
      call msg('Unknown vertical distripution type.')
    endif
    
    fTmp => fu_work_array()
    if(error)return
    fTmp(1:ps%nDescriptors)=0.
    call msg('Value dimensions (nDescriptors):',ps%nDescriptors)
   
    do j = 1, ps%nDescriptors
      call msg('Starting release rate,' + fu_name(ps%cocktail_descr_lst(j)) + ',' + &
                                 & fu_basic_mass_unit(ps%cocktail_descr_lst(j)) + '/sec:', &
           & ps%params(1)%rate_descr_unit(j) * fTmp(j))
    end do
    call msg('********************* end of POINT source **********************')

    call free_work_array(fTmp)

  END SUBROUTINE report_point_source


  !********************************************************************

  subroutine add_inventory_p_src(p_src, species_list, nSpecies)
    ! Get the species emitted by this source and make a list of
    ! them. The source must be initialized, of course.
    implicit none
    type(silam_point_source), intent(in) :: p_src
    type(silam_species), dimension(:), pointer :: species_list
    integer, intent(inout) :: nspecies

    integer :: iPar, iDesc, n
    ! There is a semantic difference between add_inventory for sources
    ! and get_inventory for descriptors: in the latter case, only a pointer is
    ! assigned.
    type(silam_species), dimension(:), pointer :: speciesPtr

    do iDesc = 1, p_src%nDescriptors
      call get_inventory(p_src%cocktail_descr_lst(iDesc), speciesPtr, n)  ! Getting the pointers!!
      if(error)return
      call addSpecies(species_list, nSpecies, speciesPtr, n)
      if(error)return
    end do

  end subroutine add_inventory_p_src


  !==================================================================
  !==================================================================
  !==================================================================
  !
  !  Some encapsulation
  !
  !==================================================================
  !==================================================================
  !==================================================================

  !==================================================================
  SUBROUTINE point_source_vertical(p_src,iTimeSlot, lyr,vertical_distr)
    IMPLICIT NONE
    TYPE(silam_point_source), INTENT(in) :: p_src
    TYPE(silja_level), INTENT(out) :: lyr
    integer, intent(in) :: iTimeSlot
    INTEGER, INTENT(out) :: vertical_distr
    if(iTimeSlot > size(p_src%params))then
      call set_error('Time slot number is too large','point_source_vertical')
      return
    endif
    lyr = p_src%params(iTimeSlot)%layerDynamic
    vertical_distr = p_src%params(iTimeSlot)%vertical_distribution
  END SUBROUTINE point_source_vertical

  !==================================================================
  FUNCTION fu_point_source_start_time(p_src)
    IMPLICIT NONE
    TYPE(silja_time) :: fu_point_source_start_time
    TYPE(silam_point_source), INTENT(in) :: p_src
    fu_point_source_start_time = p_src%params(1)%time
  END FUNCTION fu_point_source_start_time

  !==================================================================
  FUNCTION fu_point_source_end_time(p_src)
    IMPLICIT NONE
    TYPE(silja_time) :: fu_point_source_end_time
    TYPE(silam_point_source), INTENT(in) :: p_src
    fu_point_source_end_time = p_src%params(size(p_src%params))%time
  END FUNCTION fu_point_source_end_time

  !==================================================================
  FUNCTION fu_point_source_duration(p_src)
    IMPLICIT NONE
    TYPE(silja_interval) :: fu_point_source_duration
    TYPE(silam_point_source), INTENT(in) :: p_src
    fu_point_source_duration = p_src%params(size(p_src%params))%time - p_src%params(1)%time
  END FUNCTION fu_point_source_duration

  !==================================================================
  FUNCTION fu_NbrOfTimeSlots_of_p_src(p_src)
    IMPLICIT NONE
    integer :: fu_NbrOfTimeSlots_of_p_src
    TYPE(silam_point_source), INTENT(in) :: p_src
    fu_NbrOfTimeSlots_of_p_src = size(p_src%params)
  END FUNCTION fu_NbrOfTimeSlots_of_p_src

  !==================================================================
  logical function fu_point_source_defined(p_src)
    implicit none
    type(silam_point_source), intent(in) :: p_src
    fu_point_source_defined = p_src%defined == silja_true
  end function fu_point_source_defined

  !=========================================================================
  integer function fu_meteodep_model_p_src(p_src)
    implicit none
    type(silam_point_source), intent(in) :: p_src
    if(p_src%ifMeteoDependent)then
      fu_meteodep_model_p_src = fu_meteodep_model(p_src%rulesMeteoDep)
    else
      fu_meteodep_model_p_src = int_missing
    endif
  end function fu_meteodep_model_p_src

END MODULE source_terms_point

