MODULE input_data_rules

  ! Description:
  ! This module defined the weather data rules to be used by the
  ! models when getting data from the dataserver.
  !
  ! Original code : Mika Salonoja
  ! Author: Mikhail Sofiev
  !
  ! All units: SI
  ! 
  ! All degrees: real numbers in degrees (degrees and parts per
  ! hundred/parts per thousand , NOT in degrees and minutes
  ! 
  ! Language: ANSI standard Fortran 90
  !
  ! Original code: Mika Salonoja
  ! Author: Mikhail Sofiev, FMI email mikhail.sofiev@fmi.fi
   
  ! Modules used:

  USE grads_templates
  USE silam_levels
  use grids_geo

  IMPLICIT NONE

  ! The public functions and subroutines available in this module:
  PUBLIC fu_set_wdr
  public set_storage_region
  public set_top_level
  public set_bottom_level

  public defined
  PUBLIC fu_horizontal_interp_method
  PUBLIC fu_vertical_interp_method
  PUBLIC fu_time_interp_method
  PUBLIC fu_met_src  ! Former fu_source. Returns derived type meteo_data_source
  public fu_mds      ! Returns the pointer to the met_src from the particular wdr
  PUBLIC fu_fname_template
  PUBLIC fu_fname_oro_template
  PUBLIC fu_file_format
  PUBLIC fu_file_oro_format
  PUBLIC fu_fname_templates
  PUBLIC fu_fname_oro_templates
  PUBLIC fu_file_formats
  PUBLIC fu_file_oro_formats
  public fu_land_use_descriptors
  public fu_if_input_file
  PUBLIC fu_top_level
  PUBLIC fu_bottom_level
  PUBLIC fu_start_time
  PUBLIC fu_obstime_interval
  public fu_period_to_compute
  public fu_index     ! Index of mds in wdr/stackvector - see interface
  public fu_NbrOfMetSrcs ! Number of the meteo data sources
  public fu_storage_area
  public fu_storage_grid
  public fu_ablh_method 
  public fu_abl_param  
  public fu_kz_param 
  public fu_abl_min_m
  public fu_precipitation_low_limit
  public fu_if_zero_fc_len_allowed
  public fu_ifWait
  public fu_ifDisregardMDS
  public fu_obstimes
  public fu_closest_obstime
  public fu_if_randomise

  public fu_set_met_src
  public add_mds
  public fu_name
  public fu_centre
  public fu_subCentre
  public fu_model
  public fu_meteo_time_shift
  public set_meteo_time_shift

  public report

  public FNm_from_single_template
  public FNm_from_template_array
  public fu_sample_file_name

  public fu_3d_leveltype
  public fu_number_of_precip_flds

  public fu_residence_interval
  public fu_next_available_meteo_time

  ! The private functions and subroutines not to be used elsewhere:
  private fu_set_wdr_from_namelist
  private fu_set_wdr_from_params
  PRIVATE fu_wdr_defined
  PRIVATE fu_compare_wdrs_eq
  PRIVATE fu_met_src_of_wdr ! Former fu_source_of_wdr
  private fu_mds_from_mds
  private fu_mds_from_index
  PRIVATE fu_fname_template_of_wdr
  PRIVATE fu_fname_templates_of_wdr
  PRIVATE fu_fname_oro_template_of_wdr
  PRIVATE fu_fname_oro_templates_of_wdr
  PRIVATE fu_file_format_of_wdr
  PRIVATE fu_file_formats_of_wdr
  PRIVATE fu_file_oro_format_of_wdr
  PRIVATE fu_file_oro_formats_of_wdr
  PRIVATE fu_top_of_wdr
  PRIVATE fu_obstime_int_wdr
  PRIVATE fu_start_time_wdr

  private fu_set_met_src_from_params
  private fu_set_met_src_from_namelist
  private fu_compare_met_srcs_eq ! For individual meteo data sources
  private fu_name_of_mds
  private fu_centre_of_mds
  private fu_subcentre_of_mds
  private fu_model_of_mds
  private fu_mds_defined
  private mds_report
  private fu_index_met_src_of_mds
  private fu_bottom_of_wdr

  ! Generic names and operator-interfaces of some functions:\

  interface fu_set_wdr
    module procedure fu_set_wdr_from_namelist
    module procedure fu_set_wdr_from_params
  end interface

  INTERFACE defined
    MODULE PROCEDURE fu_wdr_defined
    MODULE PROCEDURE fu_mds_defined
  END INTERFACE

  INTERFACE operator(==)
    MODULE PROCEDURE fu_compare_wdrs_eq
    module procedure fu_compare_met_srcs_eq
  END INTERFACE

  INTERFACE fu_met_src  ! Former source
    MODULE PROCEDURE fu_met_src_of_wdr
  END INTERFACE

  INTERFACE fu_mds    ! Former source
    MODULE PROCEDURE fu_mds_from_mds
    MODULE PROCEDURE fu_mds_from_index
  END INTERFACE

  interface fu_index ! Index of the meteo data source
    module procedure fu_index_met_src_of_mds ! from meteo_data_source
  end interface

  INTERFACE fu_fname_template
    MODULE PROCEDURE fu_fname_template_of_wdr
  END INTERFACE

  INTERFACE fu_fname_oro_template
    MODULE PROCEDURE fu_fname_oro_template_of_wdr
  END INTERFACE

  INTERFACE fu_fname_templates
    MODULE PROCEDURE fu_fname_templates_of_wdr
  END INTERFACE

  INTERFACE fu_fname_oro_templates
    MODULE PROCEDURE fu_fname_oro_templates_of_wdr
  END INTERFACE

  INTERFACE fu_file_format
    MODULE PROCEDURE fu_file_format_of_wdr
  END INTERFACE

  INTERFACE fu_file_formats
    MODULE PROCEDURE fu_file_formats_of_wdr
  END INTERFACE

  INTERFACE fu_file_oro_format
    MODULE PROCEDURE fu_file_oro_format_of_wdr
  END INTERFACE

  INTERFACE fu_file_oro_formats
    MODULE PROCEDURE fu_file_oro_formats_of_wdr
  END INTERFACE

  INTERFACE fu_top_level  
    MODULE PROCEDURE fu_top_of_wdr
  END INTERFACE

  INTERFACE fu_bottom_level  
    MODULE PROCEDURE fu_bottom_of_wdr
  END INTERFACE

  INTERFACE fu_obstime_interval
    MODULE PROCEDURE fu_obstime_int_wdr
  END INTERFACE

  INTERFACE fu_start_time
    MODULE PROCEDURE fu_start_time_wdr
  END INTERFACE


  !
  ! Meteo data source
  !
  interface fu_set_met_src
    module procedure  fu_set_met_src_from_params
    module procedure  fu_set_met_src_from_namelist
  end interface

  interface fu_name
    module procedure fu_name_of_mds
  end interface

  interface fu_centre
    module procedure fu_centre_of_mds
  end interface

  interface fu_subCentre
    module procedure fu_subCentre_of_mds
  end interface

  interface fu_model
    module procedure fu_model_of_mds
  end interface

  interface report
    module procedure mds_report
  end interface

  !
  ! Public types with private components defined in this module:
  !
  type meteo_data_source
    private
    integer :: SrcIndex ! From 1 to max_sources - just index in stack array
                        ! as well as in the wdr 
    integer :: CentreCode ! e.g. 86 for FMI - like in GRIB
    integer :: subCentreCode ! someting, which exists in GRIB 2 specification
    integer :: ModelCode  ! e.g. 96(?) for HIRLAM - like in GRIB
    character(len=clen) :: MetSrcNm ! e.g. "Atlantic_HIRLAM_FMI"
  end type meteo_data_source
  
  type silam_fformat
    integer :: iFormat=int_missing
    character(len=clen) :: title=''
  end type silam_fformat
  
  TYPE silja_wdr
    PRIVATE
    INTEGER :: NbrOfMetSrcs
    type(meteo_data_source), dimension(max_met_srcs) :: met_srcs
    type(grads_template), dimension(max_met_files) :: fname_template, fname_oro_template
    type(silam_fformat), dimension(max_met_files) :: FileFormat, OrogrFileFormat
    character(len=fnlen), dimension(max_met_files) :: fname_land_use_descriptor
    TYPE(silja_level) :: top, bottom
    TYPE(silja_time) :: start_time
    TYPE(silja_interval) :: obstime_interval, period_to_compute, max_hole_length
    type(silja_interval) :: meteo_time_shift   ! if EnKF perturbs meteo by shifting it
    INTEGER :: horizontal_interp
    INTEGER :: vertical_interp
    INTEGER :: time_interp
    type(silam_area) :: storage_area ! Define the area/grid to be covered by the
    type(silja_grid) :: storage_grid ! meteorological data
    integer :: ablh_method   ! Way of computation of ABL height
    integer :: abl_param    ! Method of ABL parameerization
    real :: fABLlowLimit    ! Minimal ABL height for Kz
    integer :: kz_param     ! Method to treat Kz profile
    integer :: number_of_precip_flds
    integer :: LAIsrc
    real :: fPrecipitationLowLimit
    logical :: ifAllowZeroFcLen  ! If 0-long forecast may be used for simulations
    logical :: ifWait ! in case of absent data
    logical :: ifDisregardMDS  ! if several met_src are allowed for a single run
    logical :: ifRandomise = .true.
    TYPE(silja_logical) :: defined
  END TYPE silja_wdr

  type wdr_ptr
    type(silja_wdr), pointer :: ptr
  end type wdr_ptr


  ! silja_wdr%LAIsrc
  INTEGER, PUBLIC, PARAMETER :: LAI_dynamic_1 = 2001  !  single LAI
  INTEGER, PUBLIC, PARAMETER :: LAI_dynamic_2 = 2002  !  LAI_hv + LAI_lv
  INTEGER, PUBLIC, PARAMETER :: LAI_static_1 = 2003   !
  INTEGER, PUBLIC, PARAMETER :: LAI_static_2 = 2004   !

  


  !------------------------------------------------------------------------------------------

  ! Public parameters:

  INTEGER, PUBLIC, PARAMETER :: start_at_top = 1001  ! e.g. HIRLAM layer order
  INTEGER, PUBLIC, PARAMETER :: start_at_surface = 1002 ! opposite to HIRLAM order

  type(silam_fformat),parameter, public :: format_missing = silam_fformat(int_missing,'')
  type(silam_fformat),parameter, public :: grads_file_format = silam_fformat(grads_file_flag,'')



  type(meteo_data_source), parameter, public :: met_src_missing = &
    meteo_data_source(int_missing, int_missing, int_missing, int_missing, 'Undefined source')

  TYPE(silja_wdr),PARAMETER, PUBLIC :: wdr_missing=silja_wdr(0, &
                                                           & met_src_missing, &
                                                           & template_missing, template_missing, &
                                                           & format_missing, format_missing, &
                                                           & '', &
                                                           & level_missing, level_missing, &
                                                           & time_missing, &
                                                           & interval_missing, interval_missing, &
                                                           & interval_missing, interval_missing, &
                                                           & int_missing, &
                                                           & int_missing, &
                                                           & int_missing, &
                                                           & area_missing, &
                                                           & grid_missing, &
                                                           & int_missing, &
                                                           & int_missing, &
                                                           & real_missing, &
                                                           & int_missing, &
                                                           & int_missing, &
                                                           & int_missing,&  ! LAIsrc
                                                           & -1., &
                                                           & .false., &
                                                           & .false., &
                                                           & .false., &
                                                           & .true., &     ! ifRandomise
                                                           & silja_false)

  INTEGER, PRIVATE :: i

  TYPE(meteo_data_source), PARAMETER, PUBLIC :: fmi_hirlam_src = &
                & meteo_data_source( int_missing, &
                                   & centre_fmi, &
                                   & int_missing, &
                                   & model_hirlam, & 
                                   & 'FMI_HIRLAM_data_source')
  
  TYPE(meteo_data_source), PARAMETER, PUBLIC :: fmi_silam_src = &
                & meteo_data_source( int_missing, &
                                   & centre_fmi, &
                                   & int_missing, &
                                   & model_silam, & 
                                   & 'FMI_SILAM_data_source')

  TYPE(meteo_data_source), PARAMETER, PUBLIC :: silam_internal_src = &
                & meteo_data_source( int_missing, &
                                   & centre_fmi, &
                                   & int_missing, &
                                   & model_silam_internal, & 
                                   & 'SILAM_internal_data_source')
  
  TYPE(meteo_data_source), PARAMETER, PUBLIC :: centre_model_ecmwf_src = &
                & meteo_data_source( int_missing, &
                                   & centre_ecmwf, &
                                   & int_missing, &
                                   & model_ecmwf, & 
                                   & 'ECMWF_data_source')

  TYPE(meteo_data_source), PARAMETER, PUBLIC :: hacked_cosmo_src = &
                & meteo_data_source( int_missing, &
                                   & centre_moscow, &
                                   & 255, &
                                   & 200, &  !Originally 132
                                   & 'HACKED_COSMO7_data_source')

  TYPE(meteo_data_source), PARAMETER, PUBLIC :: smhi_echam_src = &
                & meteo_data_source( int_missing, &
                                   & centre_smhi, &
                                   & 98, &
                                   & 1, & 
                                   & 'SMHI_ECHAM_data_source')

  TYPE(meteo_data_source), PARAMETER, PUBLIC :: MetCoop_MEPS_src = &
                & meteo_data_source( int_missing, &
                                   & centre_metcoop, &
                                   & 255, &
                                   & 0, & 
                                   & 'MetCoop_MEPS_data_source')

  integer, parameter, private :: max_forecast_steps = 999

!  TYPE(silja_interval),PARAMETER,PRIVATE :: smallest_obstime_interval= one_hour


CONTAINS
    

  !*********************************************************************

  FUNCTION fu_set_wdr_from_namelist(nlSetup, &
                                  & nlStandardSetup, &
                                  & start_time,&
                                  & period_to_calculate, &
                                  & storage_area, &
                                  & storage_grid, &
                                  & meteo_time_shift) result(wdr)
    !
    ! Sets value for a package of weather data rules: what source of
    ! data, lengths of forecasts, interpolation methods etc. etc.
    ! Specific data sources are not set at this time - they have to be 
    ! explicitly set one-by-one using fu_set_met_src.
    ! Information is taken from the setup namelist
    !
    ! Note: Either storage_area or storage_grid MUST be present (however,
    ! may be undefined - this is a small sequirity hole left to the 
    ! responsibility of the programmer). But they definitely can not be 
    ! defined simultaneously.
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_wdr) :: wdr
    !
    ! Imported parameters with intent(in):
    type(Tsilam_namelist), pointer :: nlSetup, nlStandardSetup
    TYPE(silja_time), INTENT(in) :: start_time
    TYPE(silja_interval), INTENT(in) :: period_to_calculate
    type(silam_area), intent(in), optional :: storage_area
    type(silja_grid), intent(in), optional :: storage_grid
    TYPE(silja_interval), INTENT(in), optional :: meteo_time_shift

    ! Local declarations:
    INTEGER :: i, j, nItems
    character(len=worksize_string) :: chTmp
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: pItems ! namelist items
    character(len=*), parameter :: sub_name ="fu_set_wdr_from_namelist"
    
    ! Simple things:
    !
    wdr%start_time = start_time
    wdr%period_to_compute = period_to_calculate
    if(present(meteo_time_shift)) then
      wdr%meteo_time_shift = meteo_time_shift
    else
      wdr%meteo_time_shift = zero_interval
    endif
    wdr%NbrOfMetSrcs = 0  ! No data sources so far
    wdr%met_srcs = met_src_missing

    ! Grid and area
    !
    if(present(storage_area))then
      if(present(storage_grid))then
        if(defined(storage_area).and.defined(storage_grid))then
          call set_error('Both storage_area and storage_grid defined ',sub_name)
          return
        end if
        wdr%storage_area = storage_area  ! At least one is undefined
      wdr%storage_grid = storage_grid
        
      else
        wdr%storage_area = storage_area
        wdr%storage_grid = grid_missing
      end if

    elseif (present(storage_grid))then
      wdr%storage_area = area_missing
      wdr%storage_grid = storage_grid

    else
      call set_error('Neither storage_area nor storage_grid are present',sub_name)
      return
    end if

    ! interval btw meteotimes
    !
    wdr%obstime_interval = fu_set_named_interval(fu_content(nlSetup,'meteo_time_step'))
    ! Number of precipitation fields
    !
    i = fu_content_int(nlSetup, 'number_of_precipitation_fields')
    if (.not. (i == 1 .or. i == 2)) then
      call set_error('number_of_precipitation_fields must be either 1 or 2', &
                   & sub_name)
      return
    end if
    wdr%number_of_precip_flds = i

    ! Number of precipitation fields
    !
    i = fu_content_int(nlSetup, 'number_of_precipitation_fields')
    if (.not. (i == 1 .or. i == 2)) then
      call set_error('number_of_precipitation_fields must be either 1 or 2', &
                   & sub_name)
      return
    end if
    wdr%number_of_precip_flds = i

    ! ABL height method
    !
    SELECT CASE (fu_content(nlStandardSetup,'abl_height_method'))
      CASE ('CONSTANT')
        wdr%ablh_method= constant_abl_height
      CASE ('RICHARDSON')
        wdr%ablh_method= richardson_method
      CASE ('COMBINATION')
        wdr%ablh_method= combination_method
      CASE ('PARCEL')
        wdr%ablh_method= parcel_method
      CASE ('KZ_CORIOLIS')
        wdr%ablh_method= coriolis_method
      CASE ('NWP_ABL')
        wdr%ablh_method= nwp_abl
    CASE default
      CALL msg_warning('abl-height must be CONSTANT,RICHARDSON,PARCEL,COMBINATION or NWP_ABL',&
                     & sub_name)
      wdr%ablh_method = combination_method
    END SELECT

    ! LAI type
    !
    SELECT CASE (fu_content(nlSetup,'use_lai'))
      CASE ('STATIC1')
         wdr%LAIsrc=LAI_static_1
      CASE ('STATIC2')
         wdr%LAIsrc=LAI_static_2
      CASE ('DYNAMIC1')
         wdr%LAIsrc=LAI_dynamic_1
      CASE ('DYNAMIC2')
         wdr%LAIsrc=LAI_dynamic_2
      case ('NONE')
         wdr%LAIsrc=int_missing
      CASE default
         CALL msg_warning('use_lai should be STATIC1 STATIC2 DYNAMIC1 DYNAMIC2 or NONE;  Will use: STATIC1',&
                        & sub_name)
         wdr%LAIsrc=LAI_static_1
    END SELECT

    ! ABL parameterization method
    !
    select case (fu_content(nlSetup,'abl_parameterization_method'))
      case('FULL_PARAM')
        wdr%abl_param = abl_full_param
      case('DRY_ABL')
        wdr%abl_param = abl_dry_param
      case default
        call set_error(fu_connect_strings('Must be FULL_PARAM or DRY_ABL, not:', &
                                        & fu_content(nlSetup,'abl_parameterization_method')), &
                     & sub_name)
        return
    end select

    ! Lower limit for ABL heigt
    !
    wdr%fABLlowLimit = fu_content_real(nlStandardSetup,'abl_minimal_height')
    if(error .or. wdr%fABLlowLimit < 0 .or. wdr%fABLlowLimit > 10000. )then
        call msg('Strange or no abl_minimal_height in StandardSetup:'+ &
                & fu_content(nlStandardSetup,'abl_minimal_height'))
        call msg('Setting it to 30 m')
        wdr%fABLlowLimit = 30
        call msg_warning('abl_minimal_height must be real [0..10000]',sub_name)
        if (error) call unset_error ("fu_set_wdr_from_namelist")
    endif

!    ! Kz parameterization method
!    !
    select case (fu_content(nlStandardSetup,'kz_profile_method'))
      case('SILAM_KZ_EMULATOR')
        wdr%kz_param = silam_kz_emulator
      case('SIMPLE_KZ')
        wdr%kz_param = simple_kz
      case('SILAM_ABL_EC_FT_KZ')
        wdr%kz_param = silam_abl_ec_ft_kz
      case('SIMPLE_ABL_EC_FT_KZ')
        wdr%kz_param = simple_abl_ec_ft_kz
      case('EC_KZ')
        wdr%kz_param = ec_kz
      case('HUNTEN_KZ')
        wdr%kz_param = hunten_kz
      case('ZERO_KZ')
        wdr%kz_param = zero_kz
      case default
        CALL msg_warning("requested kz_profile_method:",fu_content(nlSetup,'kz_profile_method'))
        CALL msg_warning('kz_profile_method must be SILAM, SILAM_RESISTANCE or SILAM_EC_KZ',&
                     & sub_name)
        CALL msg_warning('Using kz_profile_method = SILAM_ABL_EC_FT_KZ')
        wdr%kz_param = silam_abl_ec_ft_kz
    end select

    ! So far, we do not take into account upper and lower level limits
    !
    wdr%top = level_missing
    wdr%bottom = level_missing

    ! Interpolation methods: horizontal
    !
    select case(fu_content(nlStandardSetup,'horizontal_interpolation'))
      case('LINEAR')
        wdr%horizontal_interp = linear

      case('NEAREST_POINT')
        wdr%horizontal_interp = nearest_point

      case('CUBIC')
        wdr%horizontal_interp = cubic

      case default
        call msg_warning(fu_connect_strings('Strange horizontal interpolation type:', &
                  & fu_content(nlStandardSetup, 'horizontal_interpolation')), &
                  & sub_name)
        wdr%horizontal_interp = linear
    end select

    ! Interpolation methods: vertical
    !
    select case(fu_content(nlStandardSetup,'vertical_interpolation'))
      case('LINEAR')
        wdr%vertical_interp = linear

      case('NEAREST_POINT')
        wdr%vertical_interp = nearest_point

      case('CUBIC')
        wdr%vertical_interp = cubic

      case default
        call msg_warning(fu_connect_strings('Strange vertical interpolation type:', &
                  & fu_content(nlStandardSetup, 'horizontal_interpolation')), &
                  & sub_name)
        wdr%vertical_interp = linear
    end select

    ! Interpolation methods: time
    !
    select case(fu_content(nlStandardSetup,'time_interpolation'))
      case('LINEAR')
        wdr%time_interp = linear

      case('NEAREST_POINT')
        wdr%time_interp = nearest_point

      case('CUBIC')
        wdr%time_interp = cubic

      case default
        call msg_warning(fu_connect_strings('Strange time interpolation type:', &
                  & fu_content(nlStandardSetup, 'horizontal_interpolation')), &
                  & sub_name)
        wdr%time_interp = linear
    end select

    !
    ! Template strings for dynamic data. Each string now has own file format specifier
    !
    nullify(pItems)
    call get_items(nlSetup, 'dynamic_meteo_file', pItems, nItems)
    if(nItems < 1)then
      call set_error('No dynamic_meteo_file string in control file',sub_name)
      return
    endif
    if(fu_content(pItems(1)) == '')then
      call set_error('Empty dynamic_meteo_file string in control file',sub_name)
      return
    endif

    do i=1,nItems
      chTmp = fu_content(pItems(i))
      if(chTmp == '')exit
      if(i > size(wdr%fname_template))then
        call set_error('Too many dynamic input templates',sub_name)
        return
      endif
      !
      ! Set the file format and decode the template
      !
      wdr%FileFormat(i) = fu_input_file_format(chTmp)
      if(error)return
      if(wdr%FileFormat(i)%iFormat == test_field_value_flag)then
        call set_collection(wdr%fname_template(i), adjustl(chTmp(index(adjustl(chTmp),' ')+1:))) ! just to store the stuff
      else
        call decode_template_string(adjustl(chTmp(index(adjustl(chTmp),' ')+1:)), wdr%fname_template(i))
      endif
      call clean(pItems(i))
      if(error)return
    end do

    !
    ! Static file template (orography). Each string now has own file format specifier
    !
    call get_items(nlSetup, 'static_meteo_file', pItems, nItems)

    if(nItems < 1)then
      call set_error('No static_meteo_file string in control file',sub_name)
      return
    endif

    if(fu_content(pItems(1)) == '-')then
      wdr%fname_oro_template = wdr%fname_template
      wdr%OrogrFileFormat = wdr%FileFormat
    else
      do i=1, nItems
        chTmp = fu_content(pItems(i))
        if(chTmp == '')exit
        if(i > size(wdr%fname_oro_template))then
          call set_error('Too many static input templates',sub_name)
          return
        endif
        !
        ! Set the file format and decode the template
        !
        wdr%OrogrFileFormat(i) = fu_input_file_format(chTmp)
        if(error)return
        if(wdr%OrogrFileFormat(i)%iFormat == test_field_value_flag)then
          call set_collection(wdr%fname_oro_template(i), adjustl(chTmp(index(adjustl(chTmp),' ')+1:))) ! just to store the stuff
          wdr%fname_oro_template(i) = template_missing
        else
          call decode_template_string(adjustl(chTmp(index(adjustl(chTmp),' ')+1:)), wdr%fname_oro_template(i))
        endif
        call clean(pItems(i))
        if(error)return
      end do
    end if

    !
    ! Land use data descriptor file. Inside each file, there is a land use file reference and a set of rules
    ! how to treat it. 
    !
    call get_items(nlSetup, 'land_use_descriptor_file', pItems, nItems)

    if(nItems < 1) call msg('No land_use_descriptor_file string in control file')

    wdr%fname_land_use_descriptor(1) = ''
    do i=1, nItems
      chTmp = fu_content(pItems(i))
      if(chTmp == '')exit
      if(i > size(wdr%fname_land_use_descriptor))then
        call set_error('Too many land_use_descriptor_file lines:' + fu_str(nItems),sub_name)
        return
      endif
      wdr%fname_land_use_descriptor(i) = fu_process_filepath(chTmp)
      call clean(pItems(i))
      if(error)return
    end do
    
    !
    ! Set the cut-off limit for very small precipitation intensity
    !
    chTmp = fu_content(nlStandardSetup,'precipitation_low_limit')
    if(chTmp == '')then
      wdr%fPrecipitationLowLimit = -1.
    else
      call setNamedValue(chTmp,'mm/sec',wdr%fPrecipitationLowLimit)
    endif

    if(error)then
      call unset_error(sub_name)
      call msg_warning('Failed to set precipitation low limit',sub_name)
      wdr%fPrecipitationLowLimit = -1.
    endif
    !
    ! Set the permit for zero-long forecast horizon. If the parameter is absent or not YES/yes,
    ! the value is set to false.
    !
    wdr%ifAllowZeroFcLen = fu_str_u_case(fu_content(nlStandardSetup, &
                                                  & 'allow_zero_forecast_length')) == 'YES'
    !
    ! Set the action in case of absent weather data. Default is break the run. However,
    ! long-lasting runs might be made so that the model will wait for the data AS LONG AS NEEDED
    !
    wdr%ifWait = fu_str_u_case(fu_content(nlSetup,'if_wait_for_data')) == 'YES'
    !
    ! Some files can be missing in the meteodata. We can allow some hole of limited length
    !
    ! Moved max_hole_in_meteo_data  to meteo parameters, left backward compatibility
    chTmp = fu_content(nlSetup,'max_hole_in_meteo_data')
    if ( chTmp == '') then
      ! Try to get it from the old place 
      chTmp = fu_content(nlStandardSetup,'max_hole_in_meteo_data')
      if (chTmp == '') then
         wdr%max_hole_length = wdr%obstime_interval
         call msg("Setting max_hole_in_meteo_data to meteo_time_step: "//fu_str(wdr%max_hole_length))
      else
        wdr%max_hole_length = fu_set_named_interval(chTmp)
        call msg("Setting max_hole_in_meteo_data from standard_setup: "//fu_str(wdr%max_hole_length))
        call msg_warning("Please consider moving max_hole_in_meteo_data to meteorological_parameters", &
                                    & sub_name )
      endif
    else 
      wdr%max_hole_length = fu_set_named_interval(chTmp)
      ! Check if both namelists have it. This is an excuse to crash
      if ( fu_content(nlStandardSetup,'max_hole_in_meteo_data') /= '') then
       call set_error("max_hole_in_meteo_data in both standard_setup and meteorological_parameters", &
           &sub_name)
      endif
    endif

    !
    ! If several meteo data sources are allowed for the run. A safe answer is NO but
    ! sometimes it is safe to allow this.
    !
    wdr%ifDisregardMDS = fu_str_u_case(fu_content(nlStandardSetup,'disregard_meteo_data_sources')) == 'YES'
    !
    ! The reprojection aliasing may be smoothed by randomization
    !
    wdr%ifRandomise = .not. fu_str_u_case(fu_content(nlStandardSetup,'randomise_reprojection')) == 'NO'
    if(wdr%ifRandomise)then
      call msg('Reprojection with randomization')
    else
      call msg('Reprojection WITHOUT randomization')
    endif
    
    wdr%defined = silja_true

  END FUNCTION fu_set_wdr_from_namelist


  !****************************************************************

  function fu_set_wdr_from_params(NbrOfMetSrcs, met_srcs, fnm_multitime_tmpl, fnm_singletime_tmpl, &
                                & fFmt_multitime, fFmt_singletime, top, bottom, start_time, &
                                & obstime_interval, period_to_compute, max_hole_length,  &
                                & horizontal_interp, vertical_interp, time_interp, ablh_method, &
                                & abl_param, fABLlowLimit,  kz_param, fPrecipitationLowLimit, &
                                & ifAllowZeroFcLen, ifWait, ifDisregardMDS, storage_area, storage_grid) &
          & result(wdr)
  
    IMPLICIT NONE
    !
    ! result:
    TYPE(silja_wdr) :: wdr
    !
    ! Imported parameters with intent(in):
    !
    INTEGER, intent(in) :: NbrOfMetSrcs
    type(meteo_data_source), dimension(:), intent(in) :: met_srcs
    type(grads_template), dimension(:), intent(in)  :: fnm_multitime_tmpl, fnm_singletime_tmpl
    type(silam_fformat), dimension(:), intent(in)  :: fFmt_multitime, fFmt_singletime
    TYPE(silja_level), intent(in) :: top, bottom
    TYPE(silja_time), intent(in) :: start_time
    TYPE(silja_interval), intent(in) :: obstime_interval, period_to_compute, max_hole_length
    INTEGER, intent(in) :: horizontal_interp, vertical_interp, time_interp
    integer, intent(in) :: ablh_method, abl_param, kz_param
    real, intent(in) :: fPrecipitationLowLimit, fABLlowLimit
    logical, intent(in) :: ifAllowZeroFcLen, ifWait, ifDisregardMDS
    type(silam_area), intent(in), optional :: storage_area 
    type(silja_grid), intent(in), optional :: storage_grid


    ! Simple things:
    !
    wdr%NbrOfMetSrcs = NbrOfMetSrcs
    wdr%start_time = start_time
    wdr%period_to_compute = period_to_compute
    wdr%max_hole_length = max_hole_length
    wdr%top = top
    wdr%bottom = bottom
    wdr%obstime_interval = obstime_interval
    wdr%horizontal_interp = horizontal_interp
    wdr%vertical_interp = vertical_interp
    wdr%time_interp = time_interp
    wdr%ablh_method = ablh_method
    wdr%abl_param  = abl_param
    wdr%fABLlowLimit = fABLlowLimit
    wdr%kz_param  = kz_param
    wdr%fPrecipitationLowLimit = fPrecipitationLowLimit
    wdr%ifAllowZeroFcLen = ifAllowZeroFcLen
    wdr%ifWait = ifWait
    wdr%ifDisregardMDS = ifDisregardMDS

    ! Grid and area
    !
    if(present(storage_area))then
      if(present(storage_grid))then
        if(defined(storage_area).and.defined(storage_grid))then
          call set_error('Both storage_area and storage_grid defined ','fu_set_wdr_from_params')
          return
        end if
        wdr%storage_area = storage_area  ! At least one is undefined
        wdr%storage_grid = storage_grid

    else
        wdr%storage_area = storage_area
        wdr%storage_grid = grid_missing
      end if

    elseif (present(storage_grid))then
      wdr%storage_area = area_missing
      wdr%storage_grid = storage_grid

    else
      call set_error('Neither storage_area nor storage_grid are present','fu_set_wdr_from_params')
      return
    end if

    ! Meteo sources array
    !
    wdr%met_srcs = met_src_missing
    do i = 1, wdr%NbrOfMetSrcs
      if (i > size(met_srcs))exit
      wdr%met_srcs(i) = met_srcs(i)
    enddo


    ! Filename and format arrays
    do i = 1, size(fnm_multitime_tmpl)
      if(i>size(wdr%fname_template))then
        call set_error('Too many filenames','fu_set_wdr_from_params')
        return
      endif
      wdr%fname_template(i) = fnm_multitime_tmpl(i)
      wdr%FileFormat(i) = fFmt_multitime(i)
    enddo

    do i = 1, size(fnm_singletime_tmpl)
      if(i>size(wdr%fname_oro_template))then
        call set_error('Too many filenames','fu_set_wdr_from_params')
        return
      endif
      wdr%fname_oro_template(i) = fnm_singletime_tmpl(i)
      wdr%OrogrFileFormat(i) = fFmt_singletime(i)
    enddo

    wdr%defined = fu_set_true()

  end function fu_set_wdr_from_params


  !****************************************************************

  subroutine set_storage_region (wdr, area, grid)
    !
    ! Sets storage_area and/or storage_grid.
    ! The only condition is that at least one of them MUST be either
    ! absent or undefined. Previous values will be overwritten
    !
    ! IMPORTANT. 
    ! This is a dangerous procedure because setting / unsetting
    ! of the storage grid / area have wide implications for treatment
    ! of the GRIB files during reading / writing. All sorts of errors
    ! are waiting for you to use this function !! But if you want to be safe
    ! just set them once by fu_set_wdr and then do not touch.
    !
    implicit none
    !
    ! Imported parameters with intent IN
    type(silam_area), intent(in), optional :: area
    type(silja_grid), intent(in), optional :: grid

    ! Imported parameters with intent INOUT
    type(silja_wdr), intent(inout) :: wdr

    if(present(area))then
      if(present(grid))then ! Both present => one must be undefined
        if(defined(area).and.defined(grid))then
          call set_error('Both storage_area and storage_grid defined ','set_storage_region')
          return
        end if
        wdr%storage_area = area  ! At least one is undefined
        wdr%storage_grid = grid

      else ! grid is absent => nullify it and set the area
        wdr%storage_area = area
        wdr%storage_grid = grid_missing

      end if
    else
      if(present(grid))then ! area is absent, grid present
        wdr%storage_area = area_missing
        wdr%storage_grid = grid

      else ! both absent => nullify all, but with warning
        call msg_warning('Both storage_area and grid are nullified')

        wdr%storage_area = area_missing
        wdr%storage_grid = grid_missing
      end if
    end if

  end subroutine set_storage_region


   !*****************************************************************

  subroutine set_top_level(wdr, level)
    !
    ! Sets the top level for defined wdr.  Level may be undefined
    !
    implicit none
    !
    ! Imported parameters with intent IN
    type(silja_level), intent (in) :: level

    ! Imported parameters with intent INOUT
    type(silja_wdr), intent(inout) :: wdr

    if(defined(wdr))then
      wdr%top = level
    else
      call set_error('Undefined wdr given','set_top_level')
    end if

  end subroutine set_top_level


   !*****************************************************************

  subroutine set_bottom_level(wdr, level)
    !
    ! Sets the bottom level for defined wdr. Level may be undefined
    !
    implicit none
    !
    ! Imported parameters with intent IN
    type(silja_level), intent (in) :: level

    ! Imported parameters with intent INOUT
    type(silja_wdr), intent(inout), target :: wdr

    ! Local variables
    type(silja_wdr), pointer :: wdr_ptr

    wdr_ptr => wdr
    if(.not.defined(wdr_ptr))then
      call set_error('Undefined wdr given','set_bottom_level')
      return
    end if

    wdr_ptr%bottom = level

  end subroutine set_bottom_level


   ! ***************************************************************

  subroutine add_mds(wdr, mds)
    !
    ! Adds one more source of data into the wdr object.
    ! 
    implicit none

    ! Imported variables
    type(silja_wdr), intent(inout), target :: wdr
    type(meteo_data_source), intent(in) :: mds


    ! Local variables
    type(silja_wdr), pointer :: wdr_ptr

    wdr_ptr => wdr
    if(.not.defined(wdr_ptr))then
      call set_error('Undefined wdr','add_mds')
      return
    end if

    if(wdr_ptr%NbrOfMetSrcs > max_met_srcs - 1)then
      call set_error('Too many meteo data sources','add_mds')
      return
    end if

    wdr_ptr%NbrOfMetSrcs = wdr_ptr%NbrOfMetSrcs + 1
    wdr_ptr%met_srcs(wdr%NbrOfMetSrcs) = mds

  end subroutine add_mds

  
  ! ***************************************************************

  LOGICAL FUNCTION fu_compare_wdrs_eq(wdr1, wdr2) result(eq)

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    TYPE(silja_wdr), INTENT(in), target :: wdr1, wdr2

    ! Local variables
    type(silja_wdr), pointer :: wdr_ptr1, wdr_ptr2
    integer :: i

    wdr_ptr1 => wdr1
    wdr_ptr2 => wdr2

    IF (defined(wdr_ptr1)) THEN
      IF (defined(wdr_ptr2)) THEN
        eq = .true. 
        if((wdr_ptr1%NbrOfMetSrcs == wdr2%NbrOfMetSrcs).and.&
            & (fu_cmp_levs_eq(wdr_ptr1%top, wdr_ptr2%top)).and.&
            & (wdr_ptr1%horizontal_interp == wdr_ptr2%horizontal_interp).and.&
            & (wdr_ptr1%vertical_interp == wdr_ptr2%vertical_interp).and.&
            & (wdr_ptr1%time_interp == wdr_ptr2%time_interp).and. &
            & (wdr_ptr1%storage_area == wdr_ptr2%storage_area).and. &
            & (wdr_ptr1%storage_grid == wdr_ptr2%storage_grid)) then
          do i=1,wdr_ptr1%NbrOfMetSrcs
            eq = eq .and. wdr_ptr1%met_srcs(i) == wdr_ptr2%met_srcs(i)
          end do
        else
          eq = .false. 
        end if

      ELSE
        eq = .false.
      END IF
    ELSE
      eq = .NOT.defined(wdr_ptr2)
    END IF

  END FUNCTION fu_compare_wdrs_eq
 

  ! ***************************************************************

  FUNCTION fu_met_src_of_wdr(wdr, indexMDS) result(met_src)
    
    IMPLICIT NONE
  
    type(meteo_data_source) :: met_src

    ! Imported parameters with intent(in):
    TYPE(silja_wdr), INTENT(in) :: wdr
    integer, intent(in) :: indexMDS
    
    if(indexMDS < 0 .or. indexMDS > wdr%NbrOfMetSrcs)then
      met_src = met_src_missing
      call Msg_warning ('Invalid index of the met_src','fu_met_src_of_wdr')
    else
      met_src = wdr%met_srcs(indexMDS)
    end if

  END FUNCTION fu_met_src_of_wdr


  ! ***************************************************************

  FUNCTION fu_mds_from_mds (wdr, mds) result(mds_out)
    
    IMPLICIT NONE
  
    type(meteo_data_source) :: mds_out

    ! Imported parameters with intent(in):
    TYPE(silja_wdr), INTENT(in), target :: wdr
    type(meteo_data_source), intent(in) :: mds
    
    ! Local variables
    integer :: i
    
    do i=1, wdr%NbrOfMetSrcs
      if(mds == wdr%met_srcs(i))then
        mds_out = wdr%met_srcs(i)
        return
      end if
    end do
    
!    print *, trim(fu_name(mds))
!    call msg_warning('Unknown data source','fu_mds_ptr')
    mds_out = met_src_missing

  END FUNCTION fu_mds_from_mds


  ! ***************************************************************

  FUNCTION fu_mds_from_index (wdr, indexMDS) result(mds_out)
    
    IMPLICIT NONE
  
    type(meteo_data_source) :: mds_out

    ! Imported parameters with intent(in):
    TYPE(silja_wdr), INTENT(in), target :: wdr
    integer, intent(in) :: indexMDS
    
    if(indexMDS < 0 .or. indexMDS > wdr%NbrOfMetSrcs)then
      mds_out = met_src_missing
      call Msg_warning ('Invalid index of the met_src','fu_mds_ptr_from_index')
    else
      mds_out = wdr%met_srcs(indexMDS)
    end if

  END FUNCTION fu_mds_from_index


  ! ***************************************************************

  integer FUNCTION fu_index_met_src_of_mds(mds, wdr) 
    !
    ! Returns the index of the meteo data source. If it is undefined - 
    ! it makes an attempt to find it and only if fails - return int_missing
    
    IMPLICIT NONE
  
    ! Imported parameters with intent(in):
    TYPE(meteo_data_source), INTENT(in) :: mds
    type(silja_wdr), intent(in), optional :: wdr

    ! Local variables
    integer :: i
    
    fu_index_met_src_of_mds = mds%SrcIndex

    if(mds%SrcIndex == int_missing)then
      if(mds == met_src_missing) return ! Really undefined mds
      if(present(wdr))then
        do i = 1, wdr%NbrOfMetSrcs  ! Scan existing mds-s in wdr
          if(mds == wdr%met_srcs(i))then
            fu_index_met_src_of_mds = i  ! Found. Return the index
            return
          end if
        end do
      end if
    end if

  END FUNCTION fu_index_met_src_of_mds


  ! ***************************************************************

  FUNCTION fu_fname_template_of_wdr(wdr, ind) result(fnm_template)
    
    IMPLICIT NONE
  
    ! Results value
    type(grads_template) :: fnm_template

    ! Imported parameters with intent(in):
    TYPE(silja_wdr), INTENT(in) :: wdr
    integer, intent(in) :: ind
    
    if(ind <= 0 .or. ind <= size(wdr%fname_template))then
      fnm_template = wdr%fname_template(ind)
    else
      call msg('Strange index:',ind)
      call set_error('Strange index','fu_fname_template_of_wdr')
      return
    endif

  END FUNCTION fu_fname_template_of_wdr


  ! ***************************************************************

  FUNCTION fu_fname_oro_template_of_wdr(wdr,ind) result(fnm_oro_template)
    
    IMPLICIT NONE
  
    ! Results value
    type(grads_template) :: fnm_oro_template
    integer, intent(in) :: ind

    ! Imported parameters with intent(in):
    TYPE(silja_wdr), INTENT(in) :: wdr
    
    if(ind <= 0 .or. ind <= size(wdr%fname_oro_template))then
      fnm_oro_template = wdr%fname_oro_template(ind)
    else
      call msg('Strange index:',ind)
      call set_error('Strange index','fu_fname_oro_template_of_wdr')
      return
    endif

  END FUNCTION fu_fname_oro_template_of_wdr


  ! ***************************************************************

  FUNCTION fu_fname_templates_of_wdr(wdr) result(fnm_templates)
    IMPLICIT NONE
    type(grads_template), dimension(:), pointer :: fnm_templates
    TYPE(silja_wdr), INTENT(in), target :: wdr
    fnm_templates => wdr%fname_template
  END FUNCTION fu_fname_templates_of_wdr


  ! ***************************************************************

  FUNCTION fu_fname_oro_templates_of_wdr(wdr) result(fnm_oro_templates)
    IMPLICIT NONE
    type(grads_template), dimension(:), pointer :: fnm_oro_templates
    TYPE(silja_wdr), INTENT(in), target :: wdr
    fnm_oro_templates => wdr%fname_oro_template
  END FUNCTION fu_fname_oro_templates_of_wdr


  ! ***************************************************************

   function fu_file_format_of_wdr(wdr, ind)
    
    IMPLICIT NONE
  
    ! Imported parameter 
    TYPE(silja_wdr), INTENT(in) :: wdr
    integer, intent(in) :: ind
    type(silam_fformat) :: fu_file_format_of_wdr
    
    if(ind <= 0 .or. ind <= size(wdr%FileFormat))then
      fu_file_format_of_wdr = wdr%FileFormat(ind)
    else
      call msg('Strange index:',ind)
      call set_error('Strange index','fu_file_format_of_wdr')
      return
    endif

  end function fu_file_format_of_wdr


  !*****************************************************************

  logical function fu_if_input_file(fFormat)
    !
    ! Returns true if the submitted format implies input file
    ! as a source of information
    !
    implicit none
    type(silam_fformat), intent(in) :: fFormat

    select case(fFormat%iformat)
      case(grib_file_flag, grads_file_flag, ascii_file_flag, netcdf_file_flag)
        fu_if_input_file = .true.
      case(test_field_value_flag)
        fu_if_input_file = .false.
      case default
        call msg('Unknown input format:',fFormat%iformat)
        call set_error('Unknown input format','fu_if_input_file')
        fu_if_input_file = .false.
    end select

  end function fu_if_input_file


  ! ***************************************************************

  FUNCTION fu_closest_obstime(time, direction, obstime_interval) !, meteo_time_shift)
    !
    ! The returned time is a time, that there should be data available:
    ! an integer number of observation intervals from the midnight (this is how
    ! round_closest works).
    !
    IMPLICIT NONE

    ! Return value of this function:
    TYPE(silja_time) :: fu_closest_obstime

    ! Imported parameters with intent IN:
    TYPE(silja_time), INTENT(in) :: time
    INTEGER, INTENT(in) :: direction
    TYPE(silja_interval), INTENT(in) :: obstime_interval

    fu_closest_obstime = fu_round_closest(time, obstime_interval, direction) 

  END FUNCTION fu_closest_obstime


  !****************************************************************

  logical function fu_if_randomise(wdr)
    IMPLICIT NONE
    TYPE(silja_wdr), INTENT(in), target :: wdr
    fu_if_randomise = wdr%ifRandomise
  end function fu_if_randomise


  ! ***************************************************************

  function fu_file_formats_of_wdr(wdr) result(FileFormats)
    IMPLICIT NONE
    type(silam_fformat), dimension(:), pointer :: FileFormats
    TYPE(silja_wdr), INTENT(in), target :: wdr
    FileFormats => wdr%FileFormat
  end function fu_file_formats_of_wdr


  ! ***************************************************************

  function fu_file_oro_format_of_wdr(wdr, ind)
    
    IMPLICIT NONE
  
    ! Imported parameter 
    type(silam_fformat) :: fu_file_oro_format_of_wdr
    TYPE(silja_wdr), INTENT(in) :: wdr
    integer, intent(in) :: ind
    
    if(ind <= 0 .or. ind <= size(wdr%OrogrFileFormat))then
      fu_file_oro_format_of_wdr = wdr%OrogrFileFormat(ind)
    else
      call msg('Strange index:',ind)
      call set_error('Strange index','fu_file_oro_format_of_wdr')
      return
    endif

  end function fu_file_oro_format_of_wdr


  ! ***************************************************************

  function fu_file_oro_formats_of_wdr(wdr) result(FileFormats)
    IMPLICIT NONE
    type(silam_fformat), dimension(:), pointer :: FileFormats
    TYPE(silja_wdr), INTENT(in), target :: wdr
    FileFormats => wdr%OrogrFileFormat
  end function fu_file_oro_formats_of_wdr

  !****************************************************************
  
  function fu_land_use_descriptors(wdr) result(descrs)
    IMPLICIT NONE
    character(len=fnlen), dimension(:), pointer :: descrs
    TYPE(silja_wdr), INTENT(in), target :: wdr
    descrs => wdr%fname_land_use_descriptor
  end function fu_land_use_descriptors
  

  ! ***************************************************************

  FUNCTION fu_top_of_wdr(wdr) 
    IMPLICIT NONE
    type(silja_level) :: fu_top_of_wdr
    TYPE(silja_wdr), INTENT(in) :: wdr
    fu_top_of_wdr = wdr%top
  END FUNCTION fu_top_of_wdr


  ! ***************************************************************

  FUNCTION fu_bottom_of_wdr(wdr) 
    IMPLICIT NONE
    type(silja_level) :: fu_bottom_of_wdr
    TYPE(silja_wdr), INTENT(in) :: wdr
    fu_bottom_of_wdr = wdr%bottom
  END FUNCTION fu_bottom_of_wdr


   ! ***************************************************************

  FUNCTION fu_obstime_int_wdr(wdr) 
    IMPLICIT NONE
    TYPE(silja_interval) :: fu_obstime_int_wdr
    TYPE(silja_wdr), INTENT(in) :: wdr
    fu_obstime_int_wdr = wdr%obstime_interval
  END FUNCTION fu_obstime_int_wdr


   ! ***************************************************************

  FUNCTION fu_start_time_wdr(wdr) 
    IMPLICIT NONE
    TYPE(silja_time) :: fu_start_time_wdr
    TYPE(silja_wdr), INTENT(in) :: wdr
    fu_start_time_wdr = wdr%start_time
  END FUNCTION fu_start_time_wdr


  !*******************************************************************

  function fu_period_to_compute(wdr) 
    IMPLICIT NONE
    TYPE(silja_interval) :: fu_period_to_compute
    TYPE(silja_wdr), INTENT(in) :: wdr
    fu_period_to_compute = wdr%period_to_compute
  END FUNCTION fu_period_to_compute

  !*******************************************************************

  function fu_meteo_time_shift(wdr) 
    IMPLICIT NONE
    TYPE(silja_interval) :: fu_meteo_time_shift
    TYPE(silja_wdr), INTENT(in) :: wdr
    fu_meteo_time_shift = wdr%meteo_time_shift
  END FUNCTION fu_meteo_time_shift

  !*******************************************************************

  subroutine set_meteo_time_shift(wdr, meteo_time_shift) 
    IMPLICIT NONE
    TYPE(silja_interval), intent(in) :: meteo_time_shift
    TYPE(silja_wdr), INTENT(inout) :: wdr
    wdr%meteo_time_shift = meteo_time_shift
  END subroutine set_meteo_time_shift


  ! ***************************************************************

  INTEGER FUNCTION fu_horizontal_interp_method(rules) 
    IMPLICIT NONE
    TYPE(silja_wdr), INTENT(in) :: rules
    fu_horizontal_interp_method = rules%horizontal_interp
  END FUNCTION fu_horizontal_interp_method

  ! ***************************************************************

  INTEGER FUNCTION fu_vertical_interp_method(rules) 
    IMPLICIT NONE
    TYPE(silja_wdr), INTENT(in) :: rules
    fu_vertical_interp_method = rules%vertical_interp
  END FUNCTION fu_vertical_interp_method

  ! ***************************************************************

  INTEGER FUNCTION fu_time_interp_method(rules) 
    IMPLICIT NONE
    TYPE(silja_wdr), INTENT(in) :: rules
    fu_time_interp_method = rules%time_interp
  END FUNCTION fu_time_interp_method


  ! ***************************************************************

  INTEGER FUNCTION fu_ablh_method(rules) 
    IMPLICIT NONE
    TYPE(silja_wdr), INTENT(in) :: rules
    fu_ablh_method = rules%ablh_method
  END FUNCTION fu_ablh_method


  ! ***************************************************************

  INTEGER FUNCTION fu_abl_param(rules) 
    IMPLICIT NONE
    TYPE(silja_wdr), INTENT(in) :: rules
    fu_abl_param= rules%abl_param
  END FUNCTION fu_abl_param

  ! ***************************************************************

  real FUNCTION fu_abl_min_m(rules) 
    IMPLICIT NONE
    TYPE(silja_wdr), INTENT(in) :: rules
    fu_abl_min_m= rules%fABLlowLimit
  END FUNCTION fu_abl_min_m


  ! ***************************************************************
  INTEGER FUNCTION fu_kz_param(rules) 
    IMPLICIT NONE
    TYPE(silja_wdr), INTENT(in) :: rules
    fu_kz_param = rules%kz_param
  END FUNCTION fu_kz_param


  ! ***************************************************************

  integer function fu_NbrOfMetSrcs(wdr) ! Number of the meteo data sources
    implicit none
    type(silja_wdr), intent(in) :: wdr
    fu_NbrOfMetSrcs = wdr%NbrOfMetSrcs
  end function fu_NbrOfMetSrcs


  !****************************************************************
  
  function fu_storage_area(wdr) result (area)
    !
    ! Returns the storage area (possibly, area_missing)
    !
    implicit none

    ! Return value:
    type(silam_area) :: area

    ! Imported parameters with intent IN
    type(silja_wdr), intent(in), target :: wdr
     
    ! Local variables
    type(silja_wdr), pointer :: wdr_ptr

    wdr_ptr => wdr
    if(.not.defined(wdr_ptr))then
      call set_error('Undefined wdr given','fu_storage_area')
      return
    end if

    area = wdr_ptr%storage_area

  end function fu_storage_area


  !****************************************************************
  
  function fu_storage_grid(wdr) result (grid)
    !
    ! Returns the storage area (possibly, area_missing)
    !
    implicit none

    ! Return value:
    type(silja_grid), pointer :: grid

    ! Imported parameters with intent IN
    type(silja_wdr), intent(in), target :: wdr

    ! Local variables
    type(silja_wdr), pointer :: wdr_ptr

    wdr_ptr => wdr
     
    if(.not.defined(wdr_ptr))then
      call set_error('Undefined wdr given','fu_storage_grid')
      nullify(grid)
      return
    end if

    grid => wdr_ptr%storage_grid

  end function fu_storage_grid


  !******************************************************************

  real function fu_precipitation_low_limit(wdr)
    implicit none
    type(silja_wdr), intent(in), target :: wdr
    fu_precipitation_low_limit = wdr%fPrecipitationLowLimit
  end function fu_precipitation_low_limit


  !******************************************************************

  logical function fu_if_zero_fc_len_allowed(wdr)
    !
    ! Returns the permission to use anaysis fields
    !
    implicit none

    ! Imported parameters 
    type(silja_wdr), intent(in), target :: wdr

    fu_if_zero_fc_len_allowed = wdr%ifAllowZeroFcLen

  end function fu_if_zero_fc_len_allowed


  !******************************************************************

  logical function fu_ifWait(wdr)
    !
    ! Returns the requirement to wait for the meteo data if they are not found
    !
    implicit none

    ! Imported parameters 
    type(silja_wdr), intent(in), target :: wdr

    fu_ifWait = wdr%ifWait

  end function fu_ifWait

  
  !******************************************************************

  function fu_max_hole_length(wdr)
    !
    ! Returns the permission to have up to XX hours holes in forecast, 0=> no holes
    !
    implicit none

    ! return value
    type(silja_interval) :: fu_max_hole_length
    
    ! Imported parameters 
    type(silja_wdr), intent(in), target :: wdr

    fu_max_hole_length = wdr%max_hole_length

  end function fu_max_hole_length

  
  !******************************************************************

  logical function fu_ifDisregardMDS(wdr)
    !
    ! Returns the requirement to wait for the meteo data if they are not found
    !
    implicit none

    ! Imported parameters 
    type(silja_wdr), intent(in), target :: wdr

    fu_ifDisregardMDS = wdr%ifDisregardMDS

  end function fu_ifDisregardMDS


  !******************************************************************

  logical function fu_wdr_defined (wdr)
    ! 
    ! Returns a true value, if the wdr is defined
    !
    IMPLICIT NONE
    ! Imported parameters with intent IN:
    TYPE(silja_wdr), intent(in) :: wdr
    
    fu_wdr_defined = fu_true(wdr%defined)

  end function fu_wdr_defined

  !*****************************************************
  
  integer function fu_number_of_precip_flds(wdr)
    implicit none
    type(silja_wdr), intent(in) :: wdr
    fu_number_of_precip_flds = wdr%number_of_precip_flds
  end function fu_number_of_precip_flds

  !*****************************************************
  
  integer function fu_LAIsrc(wdr)
    implicit none
    type(silja_wdr), intent(in) :: wdr
    fu_LAIsrc = wdr%LAIsrc
  end function fu_LAIsrc



  ! ***************************************************************
  ! ***************************************************************
  !
  !   Meteo data source class
  !
  ! ***************************************************************



   ! ***************************************************************

  function fu_set_met_src_from_params (SrcIdx, centre, subcentre, model, name) result(mds)
    !
    ! Sets a meteo data source
    !
    implicit none

    ! Result value of the function
    type(meteo_data_source) :: mds
    
    ! Imported parameters
    integer, intent(in) :: srcIdx, centre, subcentre, model
    character(len=*), intent(in) :: name

    mds%SrcIndex = SrcIdx
    mds%CentreCode = centre
    mds%subCentreCode = subcentre
    mds%ModelCode = model
    mds%MetSrcNm = name

  end function fu_set_met_src_from_params


   ! ***************************************************************

  function fu_set_met_src_from_namelist (nlMetSrc, srcIdx) result(mds)
    !
    ! Sets a meteo data source
    !
    implicit none

    ! Result value of the function
    type(meteo_data_source) :: mds
    
    ! Imported parameters
    type(Tsilam_namelist), pointer :: nlMetSrc
    integer, intent(in), optional :: srcIdx ! Index of the data source in wdr and data stack

    if(present(SrcIdx)) then
      mds%SrcIndex = SrcIdx
    else
      mds%SrcIndex = int_missing  ! Quite bad, of course, but still better than nothing
    endif
    mds%CentreCode = fu_content_int(nlMetSrc, 'data_producing_centre_code')
    mds%ModelCode = fu_content_int(nlMetSrc, 'data_producing_model_code')
    mds%MetSrcNm = fu_content(nlMetSrc, 'data_source_name')

  end function fu_set_met_src_from_namelist


  !*******************************************************************

  logical function fu_compare_met_srcs_eq(src1, src2)
    !
    ! Compares two different sources - if they are equivalent.
    ! Criterium for the comparison is just centre, model and name.
    ! Index in the stack array does not play any role.
    !
    implicit none

    type(meteo_data_source), intent(in) :: src1, src2

    if(defined(src1))then
      if(defined(src2))then
        fu_compare_met_srcs_eq = ((src1%CentreCode == src2%CentreCode).and. &
                                & (src1%subCentreCode == src2%subCentreCode).and. &
                                & (src1%ModelCode == src2%ModelCode)) ! .and. &
!                                & (src1%MetSrcNm == src2%MetSrcNm))
      else
        fu_compare_met_srcs_eq = .false.
      end if
    else
      fu_compare_met_srcs_eq = .not. defined(src2)
    end if

  end function fu_compare_met_srcs_eq

  !****************************************************
  function fu_name_of_mds(mds) result(nm)
    implicit none
    character(len=clen) :: nm
    type(meteo_data_source), intent(in) :: mds
    nm = mds%MetSrcNm
  end function fu_name_of_mds

  !*****************************************************
  integer function fu_centre_of_mds(mds)
    implicit none
    type(meteo_data_source), intent(in) :: mds
    fu_centre_of_mds = mds%CentreCode
  end function fu_centre_of_mds

  !*****************************************************
  integer function fu_subCentre_of_mds(mds)
    implicit none
    type(meteo_data_source), intent(in) :: mds
    fu_subCentre_of_mds = mds%subCentreCode
  end function fu_subCentre_of_mds

  !*****************************************************
  integer function fu_model_of_mds(mds)
    implicit none
    type(meteo_data_source), intent(in) :: mds
    fu_model_of_mds = mds%ModelCode
  end function fu_model_of_mds


  !*****************************************************

  logical function fu_mds_defined (mds)
    ! 
    ! Returns a true value, if the meteodata source pointer is
    ! initialized and equals to non-missing value
    !
    IMPLICIT NONE
    ! Imported parameters with intent IN:
    TYPE(meteo_data_source) :: mds
    
    fu_mds_defined = mds%CentreCode /= int_missing .or. &
                   & mds%ModelCode /= int_missing
  end function fu_mds_defined


  !******************************************************

  subroutine mds_report(mds)
    !
    ! Prints the report of the meteo data source
    !
    implicit none

    ! Imported parameters with intent IN
    type(meteo_data_source), intent(in) :: mds

    if (defined(mds))then
      call msg(fu_connect_strings('Meteo data source:',mds%MetSrcNm))
      call msg('Centre:',mds%centreCode)
      call msg('subCentre:',mds%subCentreCode)
      call msg(' Model:',mds%ModelCode)
      if(mds%SrcIndex /= int_missing) call msg(' Source index: ',mds%SrcIndex)
    else
      call msg('Undefined meteo data source')
    end if

  end subroutine mds_report


  !*********************************************************************************
  !
  ! Supplementary routines
  !
  !*********************************************************************************


  ! ***************************************************************


 character(len=fnlen) function fu_sample_file_name(template, tstart, tend,  chSource, chSector)

    !
    !  Finds a file that 
    !      1. matches the template and
    !      2. exists  and
    !      3. with timestep (more or less) within time range
    !
    implicit none
    type(grads_template), intent(in) :: template
    type (silja_time),    intent(in) :: tstart, tend
    character(len=*), intent(in) :: chSource, chSector


    ! local variables
    type (silja_interval) :: timestep
    type (silja_time) :: tTry
    character(len=fnlen) :: fnTmp
    logical :: exists
    integer :: iStep, nSteps


    !!!First try to avoid limitations of fu_template_timestep
    fu_sample_file_name = fu_FNm(template, tstart, tstart, zero_interval, &
        & chSource=chSource, chSector= chSector)
!!    call msg("Trying: "//trim(fu_sample_file_name))
    INQUIRE (file=fu_sample_file_name, EXIST=exists)
    if (exists) return

    timestep = fu_template_timestep(template)

    if (timestep > one_day*366) then  !i.e.  very_long_interval
      nsteps = -1  !!No iterations for the loop below
    else
      nsteps = (tend - tstart) / timestep
    endif

    do iStep = 0, nsteps
      if (iStep == 0) then !!File timestep might be beyond range, File for tstart missing
         tTry = tstart - timestep
      else
         tTry = tstart + timestep*iStep
      endif
      fu_sample_file_name = fu_FNm(template, tTry, tTry, zero_interval, &
         & chSource=chSource, chSector= chSector)
 !     call msg("Trying: "//trim(fu_sample_file_name))
      INQUIRE (file=fu_sample_file_name, EXIST=exists)
      if (exists) return
    enddo

    call msg_warning("Failed to find sample filename ", "fu_sample_file_name")
    call msg("last tried file: "//trim(fu_sample_file_name))
    fu_sample_file_name = fu_FNm(template, tstart, tstart, zero_interval, chSource=chSource, chSector= chSector)
    call msg("first tried file: "//trim(fu_sample_file_name))
    fu_sample_file_name = ""  !return empty string 

  end function fu_sample_file_name


  ! ***************************************************************

  subroutine FNm_from_single_template(templ,valid_time,fnames, nbrOfFiles, anal_time, &
                                    & ifStrict,ifAdd,ifAllowZeroFcLen,ifWait, &
                                    & fc_step, max_hole_length)
    ! 
    ! Creates the name of the grib-coded forecast file following the given 
    ! template, valid time (mandatory) and analysis time (optional). If analysis
    ! time is present - it is considered to be fixed, which means that 
    ! everything is fixed (valid time is of course fixed). Otherwise there
    ! will be a search of existing files starting from zero-length forecast.
    !
    ! The way to soften a bit the analysis_time situation is - ifStrict variable.
    ! If it is absent or present and set to true, the analysis_time (if present)
    ! is considered to be totally fixed. Should no files happend to be with this
    ! analysis time - the error is set. But if ifStrict is present and set to false
    ! the analysis_time is considered as a recomendation only - appropriate 
    ! combination is checked first, but if nothing is found - the analysis_time
    ! is disregarded and a full search is performed.
    ! Another problem in this area is: a zero-length forecast may be dangerous because
    ! many meteo models report e.g. zero precipitation amount for the first output
    ! Therefore, one more switch ifAllowZeroFcLen allows or forbids it. Should
    ! the zero length of the forecast is forbidden, the search will be started from
    ! 1-hr forecast horizon.
    ! 
    ! ifAdd is a possibility not to destroy the existing information in fnames
    ! but rather add new names to the end of the list
    !
    ! A trick - ini_time is a totally stupid GrADS feature, so here
    ! it is substituted with something more or less reasonable - analysis_time
    !
    ! Template may contain wildcards, so this subroutine returns a list 
    ! of possible names
    !
    ! ifWait allows for absence of the input files. The system will wait until they come
    ! AS LONG AS NEEDED. This actually means that the simulations may stuck in this cycle.
    ! Dangerous thing, yes, but this seems to be the only way to handle network breaks or
    ! necessity to delete used files and copy (or unzip) the new ones. Such cycle may stay
    ! over the weekend, so in theory the waiting can be as much as 50-60 hours.
    ! 
    implicit none

    ! Imported parameters
    type(grads_template), intent(in) :: templ
    type(silja_time), intent(in) :: valid_time
    type(silam_sp), dimension(:), pointer :: fnames ! inout
    integer, intent(out), optional :: nbrOfFiles
    type(silja_time), intent(in), optional :: anal_time
    logical, intent(in) :: ifStrict, ifAdd,ifAllowZeroFcLen, ifWait
    type(silja_interval), intent(in), optional :: fc_step
    type(silja_interval), intent(in), optional :: max_hole_length  ! if some time slot is missing

    ! Local declarations:
    integer :: i, j, iName, iNbrOfFiles, iStart, max_forecast_steps, iProceed
    logical :: file_exists, lFileOK
    real :: fTmp
    type(silja_interval) :: f_length,  fc_step_my
    character(len=20) :: chReadOK
    type(silja_time)  :: anal_time_tmp
    character(len=fnlen) :: strFNmOld, strFNm
    type(fnlen_str), dimension(:), allocatable :: ptrFNames
    integer, parameter :: max_wait_count = 100
    character(len=*), parameter :: sub_name="FNm_from_single_template"


    ! Some stupidity check...
    !
!call msg('Looking for the input data files...')

    if(fu_fails(defined(valid_time),'undefined valid time',sub_name))return

    if(present(nbrOfFiles)) nbrOfFiles = 0

    ! Determine the starting point to add the new names
    !
    iName = 1
    if(associated(fnames))then
      if(ifAdd)then
        do i=1,size(fnames)
          if(len_trim(fnames(i)%sp) > 0) iName = i+1 ! Place for a new name
        enddo
      else
        do i=1,size(fnames)
          fnames(i)%sp = '' ! Destroy the name
        enddo
      endif  ! ifAdd
    endif ! associated fnames

    
    if (present(fc_step)) fc_step_my = fc_step
    
    ! try to guess reasonable forecast step
    if (.not. defined(fc_step_my)) then
      !!Trick to handle strange start times
      fTmp = modulo(silja_time_to_real8(valid_time)/3.6D3,1D0)  !!!fraction of hour
      if (fTmp < 0.5/3.6e3) then
        fc_step_my = one_hour  
      else
        !!Find denominator of rational fraction of hour
        do i= 2,1000
             if (abs(nint(fTmp*i) - fTmp*i) < 0.5/3.6e3 ) exit
        enddo
        if (i <= 1000) then !denominator found
              fc_step_my = one_hour / real(i)
        else
          call msg("Expanding template for valid time, failed to guess reasonable fc_step")
          call report(valid_time)
          call set_error("Trouble ahead", sub_name)
          return
        endif
      endif
    endif

    !
    ! Start the possibly endless cycle of trying-and-waiting if no files found
    !
    iProceed = 0
    max_forecast_steps = int(fu_set_named_interval('999 hr') / fc_step_my)

    do while(iProceed <= max_wait_count)
      !
      ! Is analysis time present and defined ? If yes - just check the only
      ! possible time combination (still may be a few files) and return.
      !
      if(present(anal_time))then
        if(defined(anal_time))then

          ! This function re-allocates the ptrFNames to a correct size
          ! and puts found file names in it
          !
!call msg('String to fnames')
          call string_to_fnames(fu_FNm(templ, anal_time, anal_time, valid_time - anal_time), &
                              & ptrFNames, &
                              & iNbrOfFiles)
!call msg('String to fnames done')
          !
          ! Checking the availability of the files for reading. Should we found at least
          ! something - fine
          !
          lFileOK = .false.
          if(iNbrOfFiles > 0) then
            do i=1, iNbrOfFiles

              inquire(file=ptrFNames(i)%s, exist=file_exists, read=chReadOK)

              if((.not.file_exists) .or. chReadOK == 'NO') then
                ptrFNames(i)%s=''
              else
                lFileOK = .true.  ! At least one file is OK for reading
              endif

            end do ! Cycle over fnames
          endif

          if(lFileOK) then  ! At least something is found

            call enlarge_array(fnames, iName + iNbrOfFiles - 1) ! More than may be in reality

            do i=1, size(ptrFNames)
              if(len_trim(ptrFNames(i)%s) > 0) then
                fnames(iName)%sp = ptrFNames(i)%s
                iName = iName + 1
              endif
            end do
!            print *, 'Done'
            if(present(nbrOfFiles)) nbrOfFiles = iName -1
            return
          end if
        
          !
          ! No data available for the given analysis time. What to do ?
          !
          if(ifStrict)then
            call set_error('Data do not exists for given analysis and valid times', &
                         & sub_name)
            return
          end if

        end if ! defined anal_times
      end if ! present anal_time

      !---------------------------------------------------------------
      !
      ! Analysis time either does not exist or is not defined, or data are absent and
      ! ifStrict is set to false, so we should search through analysis time + forecast
      ! pairs keeping just valid time.
      !
      ! A trick: it may happen that the file names are not changing every hour. Then
      ! it is profitable to check that the next name to try differs from the previous 
      ! one and only then call string_to_fnames.
      !
!      print *, 'No analysis time is given - trying reasonable combinations'

      strFNm = ''
      strFNmOld = ''
      lFileOK = .false.
      if(ifAllowZeroFcLen)then
        iStart=0
      else
        iStart=1
      endif
      
      do i = iStart, max_forecast_steps

         f_length = fc_step_my * i 
         anal_time_tmp = valid_time-f_length
         

         if (fu_min(anal_time_tmp) /= 0) cycle  ! No analysis time with minutes
              !! Ditrty hack to make it working when starting run at not whole hour
              !! 

        strFNm = fu_FNm(templ, anal_time_tmp, anal_time_tmp, f_length)
        if(error)then
          call msg('Failed to decode template for analysis time=' + fu_str(anal_time_tmp) + &
                 & ', and forecast length=' + fu_str(f_length))
          call unset_error(sub_name)
          cycle
        endif
        !
        ! Did new forecast length resulted in a new file name template?
        !
        if(strFNm /= strFNmOld)then
          call string_to_fnames (strFNm, ptrFNames, iNbrOfFiles)
          if(error)then
            call set_error('Failed to handle the string:' + strFNm, sub_name)
            call unset_error(sub_name)
            cycle
          endif
          strFNmOld = strFNm
        else
          cycle
        endif
        !
        ! Got some new names. Check for the file existence
        !
        do j=1,iNbrOfFiles

          if(ptrFNames(j)%s == '')exit ! short-cut

          inquire(file=ptrFNames(j)%s, exist = file_exists, read = chReadOK)

          if((.not.file_exists) .or. chReadOK == 'NO') then
            ptrFNames(j)%s=''
          else
            lFileOK = .true.  ! At least one file is OK for reading
          endif
        end do ! cycle over fnames

        if(lFileOK) then  ! At least something is found

          call enlarge_array(fnames, iName + iNbrOfFiles - 1)

          do j=1, iNbrOfFiles
            if(len_trim(ptrFNames(j)%s) > 0) then
              fnames(iName)%sp = ptrFNames(j)%s
              iName = iName + 1
            endif
          end do
!            print *, 'Done'
          if(allocated(ptrFNames)) deallocate(ptrFNames)
          if(present(nbrOfFiles)) nbrOfFiles = iName - 1
          return
        else
        end if  ! ifFileOK

      end do ! cycle through the forecast length
      
      !
      ! If we are here, the search for the time has failed: no files at all. Three options: 
      ! - waiting
      ! - allow hole and do nothing
      ! - start breaking dishes
      !
      if(ifWait)then
        call msg_warning('Local time now: ' + fu_computer_time_string() + ', No files for:' + &
                       & fu_str(valid_time) + ', retry in 5 min', sub_name)
        !
        ! wait for 10 minutes checking if someone sent apkill signal meanwhile. That will set error
        !
        if(iProceed == max_wait_count)then
          call set_error('Cannot find any readable input file for:' + fu_str(valid_time) + &
                         & ', too much waiting', sub_name)
          call msg("Last tried file: "//trim(strFNm))
        else
          do j = 1, 60
            call wait(5.)                         ! waiting
            if(error)return
          enddo
          iProceed = iProceed + 1
        end if
      else
        if(present(max_hole_length))then
          if(.not. (max_hole_length > zero_interval))then
            call set_error('Cannot find any readable input file for:' + fu_str(valid_time) + &
                         & ', zero hole tolerance', sub_name)
            call msg("Last tried file: "//trim(strFNm))
          endif
        else
          call set_error('Cannot find any readable input file 2 for:' + fu_str(valid_time) + &
                       & ', no hole tolerance stated', sub_name)
          call msg("Last tried file: "//trim(strFNm))
        endif
        if(allocated(ptrFNames)) deallocate(ptrFNames)
        return
      endif   ! wait?

    enddo ! while ifProceed
    
!call msg('exit FNm_from_single_template')

  end subroutine FNm_from_single_template

  
  !*******************************************************************

  subroutine FNm_from_template_array(arTemplate, arTemplFormats, &
                                   & valid_time, fnames, nFileNames, fFormats, anal_time, &
                                   & ifStrict, ifAllowZeroFcLen, ifWait, &
                                   & fc_step, max_hole_length)
    !
    ! Same as above but there is an array of templates, not a single template
    ! For real decoding, it just calls the above function for each non-missing 
    ! template and then merges the file names.
    !
    implicit none

    ! Imported parameters
    type(grads_template), dimension(:), intent(in) :: arTemplate
    type(silam_fformat), dimension(:), intent(in) :: arTemplFormats
    type(silja_time), intent(in) :: valid_time
    type(silja_time), intent(in) :: anal_time                     ! can be time_missing
    logical, intent(in) :: ifStrict, ifAllowZeroFcLen, ifWait
    type(silam_sp), dimension(:), pointer :: fnames         ! inout
    type(silam_fformat), dimension(:), pointer :: fFormats  ! inout, save size as fnames
    integer, intent(out) :: nFileNames
    type(silja_interval), intent(in) :: fc_step, max_hole_length

    ! Local variables
    integer :: i, iFilesReceived
    integer, dimension(:), pointer :: iFormatsTmp
    !
    ! Go through the whole template array, get the file names for each template
    ! and record them all to fnames
    !
    if(associated(fnames))then
      do i=1,size(fnames)
        fnames(i)%sp=''
        fFormats(i)%iFormat = int_missing
      end do
    else
      nullify(fFormats)
    endif
    nFileNames = 0
    
    iFormatsTmp => fu_work_int_array()
    
    do i=1, size(arTemplate)
      if(arTemplate(i) == template_missing)exit    ! all done
      !
      ! Some input may be not from file but e.g. be a test field
      !
      select case(arTemplFormats(i)%iFormat)
        case(grib_file_flag, ascii_file_flag, netcdf_file_flag)
          !
          ! A meaningful data file
          !
          iFilesReceived = nFileNames !!Save it
          call FNm_from_single_template(arTemplate(i), &
                                      & valid_time, &
                                      & fnames, nFileNames, &
                                      & anal_time, &
                                      & ifStrict = ifStrict, &
                                      & ifAdd = .true., & 
                                      & ifAllowZeroFcLen = ifAllowZeroFcLen, &
                                      & ifWait = ifWait, &
                                      & fc_step = fc_step, &
                                      & max_hole_length = max_hole_length)
          iFormatsTmp(iFilesReceived+1 : nFileNames) = i   ! store what format these files are 
          
        case(test_field_value_flag)  ! No files at all
          !
          ! Template stores a request of some non-file input in its collection
          !
          nFileNames = nFileNames + 1
          call enlarge_array(fnames,nFileNames)
          fnames(nFileNames)%sp = fu_collection(arTemplate(i))
          !if(size(fnames) > 1) fnames(2)%sp = ''
          iFormatsTmp(nFileNames) = i

        case default
          call set_error('Unknown input format','FNm_from_template_array')
          return
      end select

    end do  ! Cycle over template array
    !
    ! Having all file names received and their format indices stored, fill-in the formats arracy
    !
    if(nFileNames > 0)then
      if(associated(fFormats))then
        if(size(fnames) /= size(fFormats))then
          deallocate(fFormats)
          allocate(fFormats(size(fnames)))
        endif
      else
        allocate(fFormats(size(fnames)))
      endif  ! if fFormats are associated
      do i = 1, nFileNames
        fFormats(i) = arTemplFormats(iFormatsTmp(i))
      end do
    endif  ! found something

    call free_work_array(iFormatsTmp)
    
  end subroutine FNm_from_template_array


  !**********************************************************************************

  function fu_input_file_format(chFileString)
    !
    ! Returns the format of the input file
    !
    implicit none

    ! Imported parameters
    character(len=*), intent(in) :: chFileString
    type(silam_fformat) :: fu_input_file_format
    integer :: iTmp

    fu_input_file_format%title = ''
    if(index(chFileString,'GRIB') == 1)then
      fu_input_file_format%iFormat = grib_file_flag
    elseif(index(chFileString,'GRADS') == 1)then
      fu_input_file_format%iFormat = grads_file_flag
    elseif(index(chFileString,'NETCDF') == 1)then
      fu_input_file_format%iFormat = netcdf_file_flag
      if(index(chFileString,'NETCDF:') == 1)then
         iTmp = index(chFileString,' ')
          if(iTmp  == 0)then
              fu_input_file_format%title = chFileString(8:)
          else
              fu_input_file_format%title = chFileString(8:iTmp-1)
          endif
      endif
    elseif(index(chFileString,'ASCII_V1') == 1)then
      fu_input_file_format%iFormat = ascii_file_flag
    elseif(index(chFileString,'TEST_FIELD') == 1)then ! No template, only quantity name and value
      fu_input_file_format%iFormat = test_field_value_flag
    else
      call set_error(fu_connect_strings('Unknown or missing file format in:',chFileString), &
                   & 'fu_input_file_format')
      return
    endif

  end function fu_input_file_format


  !*****************************************************************************

  LOGICAL FUNCTION fu_3d_leveltype(leveltype)
    !
    ! Returns true if quantities may be available on vertical levels of given type. 
    !
    IMPLICIT NONE 
    !
    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: leveltype

    fu_3d_leveltype = .not. (leveltype == surface .or. &
                           & leveltype == top_of_the_atmosphere .or. &
                           & leveltype == mean_sea .or. &
                           & leveltype == entire_atmosphere_single_layer .or. &
                           & leveltype == no_level .or. &
                           & leveltype == any_level)

  END FUNCTION fu_3d_leveltype


  ! ***************************************************************

  FUNCTION fu_obstimes(time_lim1,time_lim2, between_times_only, obstime_interval) result(times)

    ! Description:
    ! Finds all observation times (valid times) needed to cover the
    ! whole period between time-boundaries. 
    !
    ! If between_times_only is set true, then only obstimes between
    ! boundaries are accepted. 
    !
    ! If between_times_only is set false, then also
    ! the closest obstime backwards from time_start
    ! and closest obstime forwards from time_end are
    ! included in the list. This is the normal recommended behaviour for
    ! dispersion models, when time interpolation is needed to cover the
    ! whole period from time_start to time_end.
    !
    ! If the given boundaries
    ! are obstimes themselves, then they are also the
    ! first and last obstime on the list, no matter how
    ! between_times_only is set.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_time), DIMENSION(max_times) :: times

    ! Imported parameters with intent(in):
    TYPE(silja_time), INTENT(in) :: time_lim1, time_lim2
    LOGICAL, INTENT(in) :: between_times_only ! see above
    type(silja_interval), intent(in) :: obstime_interval

    ! Local declarations:
    TYPE(silja_time) :: now, obs_start, obs_end
    TYPE(silja_time) :: time_start, time_end
    INTEGER :: i, ord
    !
    ! 1. Check
    !
    IF ((.NOT.defined(time_lim1)).or.(.NOT.defined(time_lim2))) THEN
      CALL set_error('undefined time boundaries','fu_obstimes')
      RETURN
    END IF

    times = time_missing

    IF (time_lim1 < time_lim2) THEN
      time_start = time_lim1
      time_end = time_lim2
    ELSE
      time_start = time_lim2
      time_end = time_lim1
    END IF

    !
    ! 2. Set limits of timeperiod of potential obstimes.
    !
    IF (between_times_only) THEN
      obs_start = fu_closest_obstime(time_start, forwards, obstime_interval)
      obs_end = fu_closest_obstime(time_end, backwards, obstime_interval)
    ELSE
      obs_start = fu_closest_obstime(time_start, backwards, obstime_interval)
      obs_end = fu_closest_obstime(time_end, forwards, obstime_interval)
    END IF

    IF (error) RETURN

    !
    ! 3. Find all obstimes within these limits.
    !
    now = obs_start
    
    DO i = 1, max_times
      times(i) = fu_closest_obstime(now, backwards, obstime_interval)
      IF (error) RETURN
      IF (now >= obs_end) EXIT
      now = now + obstime_interval  !smallest_obstime_interval
    END DO
    !
    ! 4. Times in order from past to future, get rid of same times.
    !
    CALL times_to_ascending_order(times, ord)
    IF (error) RETURN

  END FUNCTION fu_obstimes


  !*******************************************************************************

  function fu_residence_interval(grid) result(residenceInterval)
    !
    ! Estimates the mean time period of a mass (or particle) staying inside the domain.
    ! The interval is simply the size of the domain divided by some 20km/hr
    ! wind velocity ~6m/s. Experience of the year 2000 showed that it is 
    ! generally OK - most of the year 35-50% of largrangian particles were moving, while
    ! at some time in autumn a few had to be added.
    !
    implicit none

    ! Imported parameter and output value
    type(silja_grid), intent(in) :: grid
    type(silja_interval) :: residenceInterval

    ! Local variables
    real :: dx, dy
    integer :: nx, ny
     
    if (fu_ifLonGlobal(grid)) then
       residenceInterval = one_day * 365.25    !very_long_interval
    else
       call grid_dimensions(grid, nx, ny)
       dx = fu_dx_cell_m(grid, int(0.5*nx), int(0.5*ny))*nx
       dy = fu_dy_cell_m(grid, int(0.5*nx), int(0.5*ny))*ny
       if(dx < dy) dx = dy  ! dx is now the biggest grid dimension, [m]

       residenceInterval = fu_set_interval_sec(0.15*dx) !0.15 sec m-1 = 1/6 m/s
    endif

  end function fu_residence_interval

  
  !******************************************************************************
  
  function fu_next_available_meteo_time(meteo_time, direction, fnames, &
                                      & ifStrict, ifAllowZeroFcLen, ifWait, wdr, &
                                      & max_hole_length)
    !
    ! Checks that meteodata are available for the given meteo_time and, if not, 
    ! tries to find the next meteo time up to 6 hours forwards, which is available.
    !
    implicit none
    
    ! Return value
    type(silja_time) :: fu_next_available_meteo_time
    
    ! Imported parameters
    TYPE(silja_wdr), intent(in) :: wdr
    TYPE(silja_time), intent(inout) :: meteo_time
    integer, intent(in) :: direction
    logical, intent(in) :: ifStrict, ifAllowZeroFcLen, ifWait
    type(silja_interval), intent(in) :: max_hole_length
    type(silam_sp), dimension(:), pointer :: fnames ! inout
    
    ! local variables
    real  :: fSign     ! silja_interval can be multiplied with real only
    !
    ! Do brute force: resolve templates and check existence of at least one file
    !
    if(direction == forwards)then
      fSign = 1
    else
      fSign = -1
    endif
    fu_next_available_meteo_time = meteo_time   ! start from this one
    do while(fu_abs(fu_next_available_meteo_time - meteo_time) <= max_hole_length)
      nullify(fnames)
!!      call FNm_from_template_array(fu_fname_templates_of_wdr(wdr), fu_file_formats_of_wdr(wdr), &
 !                                & fu_next_available_meteo_time, &
 !                                & fnames, &
 !                                & ifStrict = ifStrict, ifAllowZeroFcLen = ifAllowZeroFcLen, &
 !                                & ifWait = ifWait, max_hole_length = max_hole_length)

    
      
      
      if(associated(fnames))then
        if(size(fnames) > 0)then
          deallocate(fnames)     ! OK, time exists. Return it, whether same or different one.
          return
        endif
      endif
      fu_next_available_meteo_time = fu_next_available_meteo_time + wdr%obstime_interval * fSign
    end do    ! cycle till the end of the day
    
  end function fu_next_available_meteo_time
  
END MODULE input_data_rules
