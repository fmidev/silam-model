module photolysis
  ! 
  ! Lookup table based photolysis rates
  ! 
  ! This module implements the calculation of photolysis rates following the approach and
  ! data of FinROSE. The lookup tables were regenerated August 2018 by R. Hanninen using
  ! the PHODIS (pseudo-spherical version of disort, psphodis) to re-generate the tables.
  ! The standard atmospheres were same than used in FinROSE. 
  ! 
  ! The tables are given in text files, and have the following features:
  ! 
  ! - the tables are given for 2 seasons and 5 latitude bands
  ! - each table includes its standard atmosphere (67 levels, height, pressure and ozone density are used)
  ! - the LUT dimensions are zenith angle, total O3 column relative to standard and albedo.
  ! - for albedo one may use either single albedo datafile (albedo = 0.3), or multialbedo file (currently 7 values)
  ! - the seasons and regions are hard-coded. The standard atmospheres and dimensions are defined 
  !   from the input.
  ! - the rates in the table are base-10 logarithms.
  ! 
  ! UNITS NOT ALWAYS SI. Observe variable names and/or comments.

  use cocktail_basic
  
  implicit none
  
  private

  public test_lut_init
 !  public test_lut_rates
  public test_lut_rates_full

  public init_photolysis_lut
  public get_photorates_column

  public test_effective_albedo_v2
 
  public photolysis_input_needs

  public compute_water_cloud_properties


  ! Indices for each reaction in the photodissociation rate array:
  ! 
  integer, parameter, public :: &
       & pd_o2 = 1, &
       & pd_o3 = 2, & 
       & pd_h2o = 3, &
       & pd_n2o = 4, &
       & pd_ch4 = 5, &
       & pd_no2 = 6, &
       & pd_hno3 = 7, &
       & pd_hocl = 8, &
       & pd_ho2no2 = 9, &
       & pd_clono2 = 10, &
       & pd_n2o5 = 11, &
       & pd_o3_o1d = 12, &
       & pd_h2o2 = 13, &
       & pd_oclo = 14, &
       & pd_cl2o2 = 15, &   !RH2018: NOTE that FinROSE states that this results Cl+ClOO (which should further break to Cl+Cl+O2), and not ClO+ClO
       & pd_hcl = 16, &
       & pd_cl2 = 17, &
       & pd_co2 = 18, &
       & pd_clno2 = 19, &
       & pd_brono2 = 20, &
       & pd_brcl = 21, &
       & pd_hobr = 22, &
       & pd_ch3br = 23, &
       & pd_ch3cl = 24, &
       & pd_cfc11 = 25, &
       & pd_cfc12 = 26, &
       & pd_ccl4 = 27, &
       & pd_ch3ccl3 = 28, &
       & pd_hono = 29, &
       & pd_hcho_2h = 30, & ! CO + 2H
       & pd_hcho_h2 = 31, & ! CO + H2
       & pd_no = 32, &
       & pd_ho2no2_oh_no3 = 33, &
       & pd_no3_no_o2 = 34, &
       & pd_no3_no2_o = 35, &
       & pd_bro = 36, &
       & pd_clono2_cl_no3 = 37, &
!Adding the remaining Phodis photolysis rates.
       & pd_ch3ooh = 38, &
       !& pd_ch3cho = 39, &
       & pd_ald2 = 39, &          !Use acetaldehyde CH3CHO as presentative of aldehydes (CB05: ALD2 = CH3CHO)
       & pd_ch3coc2h5 = 40, &
       & pd_nacl = 41, &
       & pd_ccl2o = 42, &
       & pd_cclfo = 43, &
       & pd_cf2br2 = 44, &        !Halon-1202
       & pd_cf2brcf2br = 45, &    !Halon-2402
       & pd_cf2clbr = 46, &       !Halon-1211
       & pd_cf2clcf2chfcl = 47, & !HCFC-225cb
       & pd_cf2clcf2cl = 48, &    !CFC-114
       & pd_cf2clcfcl2 = 49, &    !CFC-113
       & pd_cf2o = 50, &
       & pd_cf3br = 51, &         !Halon-1301
       & pd_cf3cf2chcl2 = 52, &   !HCFC-225ca
       & pd_cf3cf2cl = 53, &      !CFC-115
       & pd_cf3chcl2 = 54, &      !HCFC-123
       & pd_cf3chfcl = 55, &      !HCFC-124
       & pd_ch3cf2cl = 56, &      !HCFC-142b
       & pd_ch3cfcl2 = 57, &      !HCFC-141b
       & pd_ch3cococh3 = 58, &
       !& pd_ch3cohco = 59, &
       & pd_mgly = 59, &          !NOTE: Phodis states this as ch3cohco, but likely ch3cocho = MGLY (correct order for J-value)
       & pd_chocho = 60, &
       & pd_chclf2 = 61, &        !HCFC-22
       & pd_chbr3 = 62, &
       & pd_cl2o = 63, &
       & pd_cl2o3 = 64, &
       & pd_cl2o4 = 65, &
       & pd_cl2o6 = 66, &
       & pd_clno = 67, &
       & pd_clono = 68, &
       & pd_cloo = 69, &
       & pd_ocs = 70, &
       & pd_cf3i = 71, &
       & pd_pan = 72, &
       & pd_fno = 73, &
       & pd_ch3ocl = 74, &
       & pd_o2_o_o1d = 75, &
       & pd_br2 = 76, &
       & pd_brno2 = 77, &
       & pd_c2h5cho = 78, &        !CB05: use for ALDX
       & pd_c3h7ono2 = 79, &       !isopropyl nitrate. Needed for CB05.
       & pd_ch2chcho = 80, &       !acrolein = propenal. Needed for CB05.
       & pd_ch3coooh = 81, &       !PACD: peroxyacetic acid/peracetic acid. Needed for CB05
       & pd_clo = 82, &
       & pd_fmcl = 83, &
       & pd_panx = 84, &
       & pd_open = 85
  integer, parameter, public :: maxPhotoIndex = 85 !Remember to increase this if number of possible photolysis reactions increases      

  ! Cloud models for photolysis
  integer, parameter, public :: simple_cloud = 1001
  integer, parameter, public :: detailed_cloud = 1002
  integer, parameter, public :: fake_cloud = 1003

  integer, private, save :: num_reactions=int_missing, num_lut_levs=int_missing, lut_size_sza=int_missing, &
       & lut_size_alb=int_missing, lut_size_o3=int_missing

  ! lut indices: (ind_react, ind_lat, ind_season, ind_alb, ind_o3, ind_sza, level)
  real, dimension(:,:,:,:,:,:,:), save, allocatable, private :: lut_data 
  real, dimension(:), pointer, save, private :: lut_sza_rad, lut_o3, lut_alb ! sza is in radians
  real, dimension(:,:,:), allocatable, save, private :: atm_o3_cuml, atm_press

  integer, parameter, private :: num_seasons = 2, num_regions = 5
  character(len=3), dimension(num_seasons), parameter, private :: seasons = (/'sum', 'win'/)
  !character(len=2), dimension(num_regions), parameter, private :: regions = (/'np', 'nm', 'tr', 'sm', 'sp'/) !ORIGINAL ver5.5
  character(len=2), dimension(num_regions), parameter, private :: regions = (/'sp', 'sm', 'tr', 'nm', 'np'/)  !Corrected July2017 by R.H.

  integer, private, save :: ind_summer, ind_winter

  integer, public, pointer, save :: imet_albedo, imet_cwc3d, imet_press, imet_cwcol, imet_lat, &
       & imet_tcc, imet_cc3d, imet_temp, imet_lcwc3d, imet_lcwcol, imet_dx_size, imet_dy_size, imet_dz_size, &
       & imet_airmass, imet_airdens, imet_prec, imet_cwc, imet_cic

  logical, private, save :: initialized = .false.

contains

  !************************************************************************************

  subroutine init_photolysis_lut(filename_lut)
      ! 
      ! Read the lookup data from a given file, and set the module variables. After this, the rates
      ! can be used.
      implicit none
      character(len=*), intent(in) :: filename_lut

      integer :: stat, file_unit, ind_season, ind_region, j
      type(Tsilam_namelist), pointer :: nl_lut
      character(len=*), parameter :: sub_name = 'init_photolysis_lut'!, lut_id_req = 'hammo-std-atm'
      !integer, parameter :: num_reactions_req = 37
      character(len=worksize_string) :: content
      
      call msg('')
      call msg('Initializing photolysis lookup tables from '//trim(filename_lut))
      call msg('')

      ! read in the LUT data
      file_unit = fu_next_free_unit()
      !filename_lut = fu_process_filepath(fu_content(nl_setup, 'photolysis_data_file'), must_exist=.true.)
      if (error) return
      open(file=filename_lut, unit=file_unit, form='formatted', iostat=stat)
      if (fu_fails(stat == 0, 'Failed to open LUT data file', 'init_photolysis_lut')) return

      nullify(lut_o3, lut_alb, lut_sza_rad)
      !do
      do j=1,num_seasons*num_regions !just to avoid extra warning when reaching the end of the data.
        !nl_lut => fu_read_namelist(file_unit, .false., 'BEGIN_LUT_DATA')
        nl_lut => fu_read_namelist(file_unit, .false., 'BEGIN_LUT_DATA',chNotNamedLine='BEGIN_LUT') !to avoid warnings
        if (.not. associated(nl_lut) .or. empty(nl_lut)) exit
        
        ! check & allocate lut axes o3, albedo and zenith:
        if (check_lut_val(nl_lut, 'lut_o3', lut_size_o3, lut_o3)) return
        if (check_lut_val(nl_lut, 'lut_alb', lut_size_alb, lut_alb)) return
        if (check_lut_val(nl_lut, 'lut_sza', lut_size_sza, lut_sza_rad)) return


        ! check the size of reference atmosphere and number of reactions. The order of
        ! reactions not checked.
        if (check_param(nl_lut, 'num_atm_levs', num_lut_levs)) return
        if (check_param(nl_lut, 'num_reactions', num_reactions)) return
        !if (fu_fails(num_reactions == num_reactions_req, 'num_reactions not matching', sub_name)) return

        if (.not. allocated(lut_data)) then
          allocate(lut_data(num_reactions, num_regions, num_seasons, lut_size_alb, lut_size_o3, &
                          & lut_size_sza, num_lut_levs), stat=stat)
          if (fu_fails(stat == 0, 'Allocate failed (lut_data)', sub_name)) return
        end if
        
        content = fu_content(nl_lut, 'season_id')
        ind_season = fu_index(content, seasons)
        if (fu_fails(ind_season > 0, 'Cannot read season from namelist', sub_name)) return
        
        content = fu_content(nl_lut, 'region_id')
        ind_region= fu_index(content, regions)
        !!!call msg('RHtest: ind_region = ',ind_region) !RISTO TEST
        if (fu_fails(ind_region > 0, 'Cannot read region from namelist', sub_name)) return
          
        if (error) return
        call set_atmosphere(nl_lut, ind_region, ind_season)
        if (error) return
        call set_lut(file_unit, nl_lut, lut_data(:, ind_region, ind_season, :, :, :, :))            
        if (error) return
        if (error) return
        call destroy_namelist(nl_lut)
      end do
      close(file_unit)
      
      ! convert sza entries to radians in zenith angle, not inverse cosine
      lut_sza_rad = acos(1.0 / lut_sza_rad)
      !print *, 'zenith angles finally:', lut_sza_rad
      
      ind_summer = fu_index('sum', seasons)
      if (fu_fails(ind_summer > 0, 'Failed to set ind_summer', sub_name)) return
      ind_winter = fu_index('win', seasons)
      if (fu_fails(ind_summer > 0, 'Failed to set ind_winter', sub_name)) return

      initialized = .true.

  contains

    subroutine set_lut(file_unit_lut, nl_lut, lut_values)
      implicit none
      type(Tsilam_namelist), pointer :: nl_lut
      integer, intent(in) :: file_unit_lut
      real, dimension(:,:,:,:,:), intent(out) :: lut_values

      integer :: ind_item, ind_react, ind_alb, ind_sza, iostat, ind_o3, num_items, num_items_req, &
           & stat, line_count
      character(len=1024), pointer :: lut_line
      real, dimension(:), pointer :: lut_val_tmp
      logical :: eof

      lut_values = real_missing
      
      ! the application of LUT assumes at least two points for each dimension, even if only
      ! one was ever used. For albedo there's special treatment, but at the only either/or
      ! is supported.
      !if (fu_fails(lut_size_alb == 1, 'Only one albedo supported', 'set_lut')) return
      !if (fu_fails(lut_size_alb > 1, 'Multiple values for albedo needed', 'set_lut')) return
      if (fu_fails(lut_size_alb > 0, 'At least one albedo needed', 'set_lut')) return
      if (fu_fails(lut_size_o3 > 1, 'Multiple values for O3 variation needed', 'set_lut')) return
      if (fu_fails(lut_size_sza > 1, 'Multiple sza values are required', 'set_lut')) return
            
      num_items_req = lut_size_o3*lut_size_alb*lut_size_sza*num_reactions

      allocate(lut_line, stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed (lut_line)', 'set_lut')) return

      lut_val_tmp => fu_work_array()
      line_count = 0
      do
        call next_line_from_input_file(file_unit_lut, lut_line, eof)
        if (eof .or. lut_line == 'END_LUT_DATA' .or. error) exit
        read(unit=lut_line, fmt=*, iostat=iostat) ind_react, ind_alb, ind_o3, &
             & ind_sza, lut_val_tmp(1:num_lut_levs)
        if (fu_fails(iostat == 0, 'Failed to parse lut_val', 'set_lut')) return
        lut_values(ind_react, ind_alb, ind_o3, ind_sza, 1:num_lut_levs) = lut_val_tmp(1:num_lut_levs)
        line_count = line_count + 1
      end do
      if (fu_fails(num_lut_levs > 1, 'More than one level required for LUT', 'set_lut')) return
      if (fu_fails(line_count == num_items_req, 'Wrong number of lut data lines', 'set_lut')) return
      if (any(lut_values .eps. real_missing)) then
        call set_error('Not all LUT values set', 'set_lut')
      end if
      
      call free_work_array(lut_val_tmp)
      deallocate(lut_line)

    end subroutine set_lut

    logical function check_lut_val(nl_lut, key, val_size, values) result(bad_value)
      implicit none
      type(Tsilam_namelist), pointer :: nl_lut
      character(len=*), intent(in) :: key
      integer, intent(out) :: val_size
      real, dimension(:), pointer :: values
      
      character(len=worksize_string) :: line
      real, dimension(:), pointer :: values_new
      integer :: num_values_new, stat

      bad_value = .true.
      values_new => fu_work_array()
      line = fu_content(nl_lut, key)
      call split_string(line, ' ', values_new, num_values_new)
      if (error) return

      if (fu_fails(num_values_new > 0, 'No values for ' // trim(key), 'check_lut_val')) return
      
      if (.not. associated(values)) then
        ! first entry processed
        allocate(values(num_values_new), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed', 'check_lut_val')) return
        val_size = num_values_new
        values(1:val_size) = values_new(1:val_size)
      else
        if (num_values_new /= size(values)) then
          call set_error('Sizes do not match for ' // trim(key), 'check_lut_val')
          return
        end if
        if (.not. all(values .eps. values_new(1:num_values_new))) then
          call set_error('Values do not match for ' // trim(key), 'check_lut_val')
          return
        end if
      end if
      bad_value = .false.
      call free_work_array(values_new)
      
    end function check_lut_val

    subroutine set_atmosphere(nl_lut, ind_region, ind_season)
      implicit none
      type(Tsilam_namelist), pointer :: nl_lut
      integer, intent(in) :: ind_region, ind_season

      type(Tsilam_nl_item_ptr), dimension(:), pointer :: items
      integer :: num_items, ind_lev, ind_item, stat
      real, dimension(:), pointer :: o3_profile, hgt
      real :: dz_cm, hgt_km, press_hpa, o3_molec_cm3, scale_hgt_cm, tempr
      real, parameter :: molar_mass_o3 = 48.0
      character(len=worksize_string) :: content
      
      nullify(items)
      call get_items(nl_lut, 'std_atm', items, num_items)
      if (fu_fails(num_items == num_lut_levs, 'Wrong number of levels in lut', 'set_atmosphere')) return
      
      o3_profile => fu_work_array()
      o3_profile(1:num_lut_levs) = real_missing
      if (.not. allocated(atm_o3_cuml)) then
        allocate(atm_o3_cuml(num_regions, num_seasons, num_lut_levs), &
               & atm_press(num_regions, num_seasons, num_lut_levs), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed', 'set_atmosphere')) return
      end if

      hgt => fu_work_array()
      do ind_item = 1, num_items
        content = fu_content(items(ind_item))
        read(unit=content, fmt=*, iostat=stat) ind_lev, hgt_km, tempr, press_hpa, o3_molec_cm3
        if (fu_fails(stat == 0, 'Failed to parse std_atm entry', 'set_atmosphere')) return
        atm_press(ind_region, ind_season, ind_lev) = press_hpa * 100
        o3_profile(ind_lev) = o3_molec_cm3
        hgt(ind_lev) = hgt_km * 1000
      end do
      if (any(o3_profile(1:num_lut_levs) .eps. real_missing)) then
        call set_error('Profiles not fully set', 'set_atmosphere') 
        return
      end if
      
      ! Make the cumulative O3 profile (in molec/cm2)
      !
      !scale_hgt_cm = 7.0e5  !Where does this come from (must depend on the uppermost layer)?
      !The following is taken from the phodis code. Height of the uppermost leyer must be in km, therefore 1e-3.
      scale_hgt_cm = hgt(num_lut_levs)*1e-3 &
             & *alog(atm_press(ind_region,ind_season,num_lut_levs)/atm_press(ind_region,ind_season,1))
      ! the following is from finrose...
      atm_o3_cuml(ind_region, ind_season, num_lut_levs) &
           & = o3_profile(num_lut_levs) * scale_hgt_cm * molecular_weight_air/molar_mass_o3
      do ind_lev = num_lut_levs - 1, 1, -1
        dz_cm = (hgt(ind_lev+1) - hgt(ind_lev)) * 1e2
        atm_o3_cuml(ind_region, ind_season, ind_lev) &
             & = atm_o3_cuml(ind_region, ind_season, ind_lev+1) &
             & + 0.5 * (o3_profile(ind_lev)+o3_profile(ind_lev+1)) * dz_cm
      end do
      
      !RISTO TEST 15Nov2018
      call msg('RISTO TEST O3: region, season, total O3 column [DU]', &
           & (/real(ind_region), real(ind_season), atm_o3_cuml(ind_region, ind_season, 1)/2.68683701e+16  /) )
      !END: RISTO TEST 15Nov2018
      call free_work_array(o3_profile)
      call free_work_array(hgt)
      deallocate(items)

    end subroutine set_atmosphere

    logical function check_param(nl, param_name, param) result(param_bad)
      implicit none
      type(Tsilam_namelist), pointer :: nl
      character(len=*), intent(in) :: param_name
      integer, intent(inout) :: param
      
      integer :: new_value
      
      param_bad = .true.
      
      new_value = fu_content_int(nl, param_name)
      if (new_value == int_missing) then
        call set_error('Missing ' // trim(param_name), 'check_param')
        return
      end if

      if (param == int_missing) then
        param = new_value
      else if (param /= new_value) then
        call msg('lut id:' // fu_content(nl, 'lut_id'))
        call msg('season: ' // fu_content(nl, 'season_id'))
        call msg('region: ' // fu_content(nl, 'region_id'))
        call set_error('Inconsistent value for ' // trim(param_name), 'check_param')
        return
      end if
        
      param_bad = .false.

    end function check_param
  end subroutine init_photolysis_lut

  !************************************************************************************


  !************************************************************************************

!  subroutine get_photorates_column(lon, lat, now, cos_zen_angle, col_press, col_cwcabove, &
!        &cwc_totcol, aodcolumn, alb_sfc, fCloudCover, rates)
  subroutine  get_photorates_column(metdat_col, zenith_cos, alb_sfc_fixed, now, &
         & aod_ext, aod_scat, o3_col, rates, ifPhotoAOD, ifPhotoO3col, tau_above_bott, ssa, cloud_model)
    ! 
    ! Interpolate the LUT to obtain photolysis rates for each reaction in each vertical
    ! level defined by pressures in col_press. 
    ! Standard O3 profile is taken.
    implicit none
    real, dimension(:,:), intent(in) ::  metdat_col !(nMet, nLev)
    real, intent(in) :: zenith_cos, &
                     & alb_sfc_fixed !! real_missing or  constant prescribed albedo
    type(silja_time), intent(in) :: now
    real, dimension(:), intent(in) :: aod_ext, aod_scat !(nZ) extinction  of each layer (m2/m2)
    real, dimension(:), intent(in) :: o3_col   ! ozone column (DU) above each level
    real, dimension(:,:), intent(out) :: rates ! (react, level)
    logical, intent(in) :: ifPhotoAOD   !! Account for aod_ext, and aod_scat for photolysis
    logical, intent(in) :: ifPhotoO3col !! Account for total O3 column above for photolysis
    real, dimension(:), intent(in) :: tau_above_bott
    real, intent(in) :: ssa !, tcc
    integer, intent(in) :: cloud_model
    
!    real, intent(in) :: lon, lat
!    real, intent(in) :: cos_zen_angle, alb_sfc, cwc_totcol, fCloudCover
!    real, dimension(:), intent(in) :: col_press, col_cwcabove, aodcolumn! (level)

    integer :: ind_sza, ind_o3, ind_alb, ind_lat, ind_lat2, ind_lev !!!Indices in LUT
    integer :: ind_press, ind_season, ind_region, ind_react, num_levs, istepalb, icld
    real, dimension(4) :: p, q
    real, dimension(2,2,2,2) :: w4
    real :: sza, yfr, earth_dist_corr, fTmp, fTmp2
    real, parameter :: sharpness = 0.35 !sharpness of the jumps between the regions, default 0.35
    real :: wlat1, wlat2 !Weights for smoothing the jumps between the different regions. 
    real, dimension(max_levels) :: cld_att, & ! Cloud attenuation
                                 & aer_att    ! aerosol_att in in-cloud area 
    real :: lat, alb_sfc, fCloudCover, cwc_totcol   !!sfc meteo input
    real, dimension(max_levels) :: col_press, col_cwcabove !! column meteo input
    real, dimension(max_levels) :: cwc_cloud  ! Cloudy-area cwc_above [kg/m2]
    real :: albedo_aer, albedo_cld    ! Effective albedo of surface without/with clouds
    real :: o3_col_std ! ozone column above for standard atmosphere (DU)
    real :: o3factor ! ozone column above / the ozone column above for the standard atmosphere

    real, parameter :: min_cc = 0.02 !! MInimum cloud cover to calculate the attenuation
    character (len=*), parameter :: sub_name = "get_photorates_column"
    
    sza = acos(zenith_cos)

    ! Do not extrapolate to zenith angles greater than last LUT value
    !
    if (sza > lut_sza_rad(size(lut_sza_rad))) then
      rates = 0.0
      return
    end if

    num_levs = size(metdat_col,2)

    lat = metdat_col(imet_lat,1)
    if (alb_sfc_fixed == real_missing) then
      alb_sfc = metdat_col(imet_albedo,1)
    else 
      alb_sfc = alb_sfc_fixed
    endif

    col_press(1:num_levs) = metdat_col(imet_press, 1:num_levs)
    fCloudCover = metdat_col(imet_tcc, 1)

    !call test()
    !call msg('*** get_photorates_column ***', lon, lat)

    !call msg('Risto test press4column = ',col_press)
    !call msg('Risto test albedo = ',albedo)

    if (fCloudCover > min_cc) then
      if (cloud_model == fake_cloud) then
        !No science here. Clouds assumed in 900-700hPa    
        albedo_cld = 0.9
        fTmp = (alb_sfc*alb_sfc + 1.) * 0.5  !! attenuation under the cloud                 
        cld_att(1:num_levs) = fTmp + (1.-fTmp)*max(0.,min(1., (90000.-col_press(1:num_levs))/(20000.)))

      else if (cloud_model == simple_cloud .or. cloud_model == detailed_cloud) then
        cwc_cloud(1:num_levs) = metdat_col(imet_cwc3d, 1:num_levs)/fCloudCover
        cwc_totcol = metdat_col(imet_cwc3d, 1) !col_cwcabove(1)
        fTmp = cwc_totcol/fCloudCover  
        call effective_albedo_cld(cwc_cloud, fTmp, alb_sfc, zenith_cos, num_levs, cld_att, &
             & albedo_cld, cloud_model, tau_above_bott, ssa)
      else
        albedo_cld = alb_sfc  !!! Actually, should not be used 
      endif
    end if

    if (ifPhotoAOD) then
      call effective_albedo_aer(aod_ext, aod_scat, alb_sfc, zenith_cos,  num_levs, aer_att, albedo_aer)
    else
      aer_att(1:num_levs) = 1.
      albedo_aer = alb_sfc
    endif

    !RISTO TEST 8JAN2019: Investigate the effect of aerosols and clouds on the photolysis rates
    !call msg("ZZtop", (/alb_sfc, albedo_aer, albedo_cld, fCloudCover, aer_att(1:num_levs), cld_att(1:num_levs)/))
    !call msg("ZZtop1", (/alb_sfc, albedo_aer, albedo_cld, fCloudCover/))
    !call msg("ZZtop2", (/aer_att(1:num_levs)/))
    !call msg("ZZtop3", (/cld_att(1:num_levs)/))
    !END: RISTO TEST 8JAN2019: Investigate the effect of aerosols and clouds on the photolysis rates
    
    ! Find appropriate LUT from lon, lat, now
    ind_season = month_to_season(fu_mon(now))
!!$    if (lat < -60) then
!!$      ind_lat = 1
!!$    else if (lat < -30) then
!!$      ind_lat = 2
!!$    else if (lat < 30) then
!!$      ind_lat = 3
!!$    else if (lat < 60) then
!!$      ind_lat = 4
!!$    else
!!$      ind_lat = 5
!!$    end if
    !In order to avoid sharp steps in the reaction rates due to different rates at different regions
    !we weight the neighboring regions with tanh profile. NOTE: This is just a TEMPORARY FIX. 
    !Index ind_lat2 denotes the nearest neighbor regions with whom we make the smoothing.
    if (lat < -60) then !South-pole region (sp symbol in LUT-data file)
       ind_lat = 1
       ind_lat2 = 2
       wlat1=0.5*(tanh(sharpness*(abs(lat)-60.0))+1.0)
    else if (lat > 60) then !North-pole region (np symbol in LUT-data file)
       ind_lat = 5
       ind_lat2 = 4
       wlat1=0.5*(tanh(sharpness*(lat-60.0))+1.0)
    else if (lat < -30) then !south medium latitudes (sm symbol in the LUT-data file)
       ind_lat = 2
       if (lat < -45) then
          ind_lat2 = 1
          wlat1=0.5*(tanh(sharpness*(60.0+lat))+1.0)
       else
          ind_lat2 = 3
          wlat1=0.5*(tanh(sharpness*(abs(lat)-30.0))+1.0)
       end if
    else if(lat > 30) then !north medium latitudes (nm symbol in the LUT-data file)
       ind_lat = 4
       if (lat > 45) then
          ind_lat2 = 5
          wlat1=0.5*(tanh(sharpness*(60.0-lat))+1.0)
       else
          ind_lat2 = 3
          wlat1=0.5*(tanh(sharpness*(lat-30))+1.0)
       end if
    else !tropical region between -30..30 degrees (tr symbol in the LUT-data file).
       ind_lat = 3
       if (lat < 0.0) then
          ind_lat2 = 2
          wlat1=0.5*(tanh(sharpness*(30-abs(lat)))+1.0)
       else
          ind_lat2 = 4
          wlat1=0.5*(tanh(sharpness*(30.0-lat))+1.0)
       end if
    end if
    wlat2=1.0-wlat1

    !call msg('ind_lat, ind_season:', ind_lat, ind_season)
    !call msg('sza', acos(cos_zen_angle)/pi*180)
    !call msg('pressure request', col_press(1))
    ! Find indices in LUT

    !RISTO TEST 15Nov2018
    !call msg('RISTO TEST O3: ind_season', ind_season )
    !call msg('RISTO TEST O3: ind_lat', ind_lat )
    !call msg('RISTO TEST O3: num_levs', num_levs )
    !call msg('RISTO TEST O3: size(o3_col)', size(o3_col) )
    !call msg('RISTO TEST O3: size(atm_o3_cuml(ind_region, ind_season, :))', size(atm_o3_cuml(ind_lat, ind_season, :)) )
    !END: RISTO TEST 15Nov2018

    do ind_lev = 1, num_levs
      call get_interp_weights_rev(col_press(ind_lev), atm_press(ind_lat, ind_season, :), ind_press, q(4))
      !print *, 'press:', atm_press(ind_lat, ind_season, :), 'req:', col_press(ind_lev), ind_press, q(4)
      call get_interp_weights(sza, lut_sza_rad, ind_sza, q(3))
      !print *, 'sza:', sza, ind_sza, q(3)
      if (ifPhotoO3col) then
         !RISTO TEST 15Nov2018
         !Amount of ozone column above the cell, compared with standard atmosphere. Note the different units (DU vs molec/cm3)!
         !Interpolate the ozone column from standard atmosphere. Note the different altitude levels used in atm_o3_cuml 
         o3_col_std = q(4)*atm_o3_cuml(ind_lat, ind_season, ind_press) &
            & + (1.0-q(4))*atm_o3_cuml(ind_lat, ind_season, ind_press+1)
         o3_col_std = o3_col_std/2.68683701e+16 !Convert to Dobson Units
         o3factor = o3_col(ind_lev)/o3_col_std
         !call msg('RISTO TEST O3: iz, o3_col(ind_lev), o3_col_std, o3factor', (/real(ind_lev), o3_col(ind_lev), o3_col_std, o3factor /) )
         call get_interp_weights(o3factor, lut_o3, ind_o3, q(2))
         !call msg('RISTO TEST O3: iz, ind_press, ind_o3', (/ind_lev, ind_press, ind_o3 /) )
         !Since the smoothing along latitude direction becomes less important and since we should,
         !in principle, use different ozone factor for different region we remove the smoothing.
         wlat1 = 1.0
         wlat2 = 0.0
         !END: RISTO TEST 15Nov2018
      else
         ! ozone variation == 1.0 for now...
         call get_interp_weights(1.0, lut_o3, ind_o3, q(2)) !NOTE NOTE NOTE RISTO RISTO RISTO RISTO
         !print *, 'ozone:', ind_o3, q(2)
      end if
      ! albedo == 0.3 in finrose


      do icld = 0,1 !Cloudy/non-cloudy fraction
        !Allow using original single albedo lut_data files for albedo:
        if (size(lut_alb) == 1) then
           q(1) = 1.0
           istepalb = 0
           ind_alb = 1
        else
           !Interpolation for effective surface albedo
           if (icld == 0) then
             call get_interp_weights(albedo_aer, lut_alb, ind_alb, q(1))
           else if (fCloudCover > min_cc) then
             call get_interp_weights(albedo_cld, lut_alb, ind_alb, q(1))
           endif
           istepalb = 1
        end if

        !print *, 'albedo:', ind_alb, q(1)
        p = 1 - q
        w4(1,1,1,1) = product(q)
        w4(1,1,1,2) = q(1)*q(2)*q(3)*p(4)
        w4(1,1,2,1) = q(1)*q(2)*p(3)*q(4)
        w4(1,1,2,2) = q(1)*q(2)*p(3)*p(4)
        w4(1,2,1,1) = q(1)*p(2)*q(3)*q(4)
        w4(1,2,1,2) = q(1)*p(2)*q(3)*p(4)
        w4(1,2,2,1) = q(1)*p(2)*p(3)*q(4)
        w4(1,2,2,2) = q(1)*p(2)*p(3)*p(4)

!!$      w4(2,:,:,:) = 0.0 !If no albedo dependence.
        w4(2,1,1,1) = p(1)*q(2)*q(3)*q(4)
        w4(2,1,1,2) = p(1)*q(2)*q(3)*p(4)
        w4(2,1,2,1) = p(1)*q(2)*p(3)*q(4)
        w4(2,1,2,2) = p(1)*q(2)*p(3)*p(4)
        w4(2,2,1,1) = p(1)*p(2)*q(3)*q(4)
        w4(2,2,1,2) = p(1)*p(2)*q(3)*p(4)
        w4(2,2,2,1) = p(1)*p(2)*p(3)*q(4)
        w4(2,2,2,2) = product(p)

        do ind_react = 1, num_reactions
!!$        !Rates without any smoothing between different (5) regions:
!!$        rates(ind_react, ind_lev) = 10**(sum(w4(1,1:2,1:2,1:2)*lut_data(ind_react, &
!!$                                                                 & ind_lat, ind_season, &
!!$                                                                 & 1, & ! ind_alb
!!$                                                                 & ind_o3:ind_o3+1, &
!!$                                                                 & ind_sza:ind_sza+1, &
!!$                                                                 & ind_press:ind_press+1)))
          !Smooth the reaction rates near the region boundaries
          fTmp= 10**(wlat1*sum(w4(1:1+istepalb,1:2,1:2,1:2)*lut_data(ind_react, &
                                                             & ind_lat, ind_season, &
                                                             & ind_alb:ind_alb+istepalb, & ! ind_alb
                                                             & ind_o3:ind_o3+1, &
                                                             & ind_sza:ind_sza+1, &
                                                             & ind_press:ind_press+1)) &
                                    +wlat2*sum(w4(1:1+istepalb,1:2,1:2,1:2)*lut_data(ind_react, &
                                                             & ind_lat2, ind_season, &
                                                             & ind_alb:ind_alb+istepalb, & ! ind_alb
                                                             & ind_o3:ind_o3+1, &
                                                             & ind_sza:ind_sza+1, &
                                                             & ind_press:ind_press+1)))
          if (icld == 0) then
            rates(ind_react, ind_lev) = fTmp * aer_att(ind_lev) * (1-fCloudCover)  !!Cloud-free subcell
          else
            !!Cloudy subcell
            rates(ind_react, ind_lev) =  rates(ind_react, ind_lev) + fTmp * cld_att(ind_lev)  * fCloudCover 
          endif


!!!#ifdef DEBUG      
          if (.not. rates(ind_react, ind_lev) >= 0.) then
            !$OMP CRITICAL (bark_rates)
            call msg("(ind_react, ind_lev)", ind_react, ind_lev)
            call msg("rates(ind_react, ind_lev)", rates(ind_react, ind_lev) )
            call msg("(/icld,ind_lat,ind_season,ind_alb,istepalb,ind_o3,ind_sza,ind_press/)",(/icld,ind_lat,ind_season,ind_alb,istepalb,ind_o3,ind_sza,ind_press/))
            call msg("cld_att(ind_lev), alb_sfc, albedo_cld, albedo_aer, fTmp", (/cld_att(ind_lev), alb_sfc, albedo_cld, albedo_aer, fTmp/))
            call msg("lat,cos_zen_angle, cwc_totcol,fCloudCover", (/lat,zenith_cos, cwc_totcol,fCloudCover/))
            call msg("tau_above_bott", tau_above_bott)
            call set_error("Trouble with LUT-photorates", sub_name)
            !$OMP END CRITICAL (bark_rates)
            return
          endif
!!!!#endif
        enddo  !!!ind_react

        !!Reactions missing from FinRose LUT -> No albedo dependence
        if (icld == 0) then
            !Additional reaction rates not covered by Phodis:
            !NOTE: Cosine of the zenith angle can be negative (below horizon)
            !      The LUT tables cover some negative values but here these must be cutted.
            if (fCloudCover > min_cc) then
              fTmp = (aer_att(ind_lev)*(1-fCloudCover) + cld_att(ind_lev)  * fCloudCover) !!Total att factor
            else
              !cld_att  uninitialized
              fTmp = aer_att(ind_lev)
            endif
            
            !Rates for MGLY (methylglyoxal) Br2, BrNO2, ClO etc. are not in standard (ver5.5 with 37 reactions) FinROSE LUT tables
            !Take the old value in case when the table does not cover the species:
            if (pd_br2 > num_reactions) then
               rates(pd_br2, ind_lev)   = photo(4.37e-2,0.0369,0.2606,zenith_cos) * fTmp !function photo correctly cuts the negative cos_zen_angle.
                                       !rate=0.0437*cos_theta**0.0369*exp(-0.2606/cos_theta), Kanaya et al., Atm. Env. 37, 2463 (2003).
            end if
            
            fTmp = fTmp * MAX(zenith_cos,0.0) !!Total att factor * function of sza            
            if (pd_mgly > num_reactions) then
               rates(pd_mgly, ind_lev)  = 1.654E-04*fTmp !NOTE: Phodis rate pd_ch3cohco = 59, is likely this!
            end if
            if (pd_ald2 > num_reactions) then
               rates(pd_ald2, ind_lev)  = 4.000E-06*fTmp 
            end if
            if (pd_brno2 > num_reactions) then
               rates(pd_brno2, ind_lev) = 2.000E-2*fTmp  !tau is order of minute: Roberts at al., Atmos. Chem. Phys. 14, 11201 (2014)
            end if
            if (pd_clo > num_reactions) then
               rates(pd_clo, ind_lev)   = 1.300E-4*fTmp  !25km value from Table B.4 on page 600 on Jacobson's book. !Currently not needed in strato.
            end if
            if (pd_fmcl > num_reactions) then
               rates(pd_clo, ind_lev)   = 1.0E-4*fTmp !TEST stetsoni-arvo!
            end if
            if (pd_panx > num_reactions) then
               rates(pd_panx, ind_lev)   = 1.0E-4*fTmp !TEST stetsoni-arvo!
            end if
            !rates(pd_open, ind_lev)  = 5.334E-05*fTmp
            
!#ifdef DEBUG
            !! Quite a few rates are uninitialised doe far. The check triggers
             !! randomly on them
!            In my case (forced ununutialized to NaN) the output of below was:
!num_reactions+1 ,maxPhotoIndex-1, ind_lev        38        84         1
!rates(num_reactions+1:maxPhotoIndex-1, ind_lev)          NaN  0.39307E-08          NaN          NaN          NaN          NaN          NaN          NaN          NaN          NaN          NaN          NaN          NaN          NaN          NaN          NaN          NaN          NaN          NaN          NaN          NaN  0.16254E-06          NaN          NaN          NaN          NaN          NaN          NaN          NaN          NaN          NaN          NaN          NaN          NaN          NaN          NaN          NaN          NaN  0.72568E-06  0.19654E-04          NaN          NaN          NaN          NaN  0.98268E-07          NaN  0.98268E-07
!

            !Checked only up to maxPhotoIndex-1 since the last one (pd_open) is currently set by using the value of pd_no2.
            if (.not. all(rates(num_reactions+1:maxPhotoIndex-1, ind_lev) >= 0.)) then
              call msg("num_reactions+1 ,maxPhotoIndex-1, ind_lev", (/num_reactions+1 ,maxPhotoIndex-1, ind_lev/))
              call msg("rates(num_reactions+1:maxPhotoIndex-1, ind_lev)",rates(num_reactions+1:maxPhotoIndex-1, ind_lev))
              call set_error("Trouble with non-LUT-photorates", sub_name)
              return
            endif
!#endif
        endif


        if (.not. fCloudCover > min_cc) exit

      enddo !!!icld
      !Extra photolysis rates that are not in LUT but that can be obtained from other photolysis rates:
      !For OPEN (aromatic opening ring product) we use the CB06 recommendation J_OPEN = 0.028*J_NO2
      rates(pd_open,ind_lev) = 0.028*rates(pd_no2,ind_lev)
      
    end do   !!!!ind_lev

    ! Correction for distance between earth and sun, from finrose. About +-3.5% max.
    ! yfr = fu_julian_date_real(now) / 365.0
    ! eart_dist_corr = 1.000110+0.034221*cos(yfr)+0.001280*sin(yfr)+0.000719*cos(2.*yfr)+0.000077*sin(2.*yfr)
    ! rates = earth_dist_corr * 10**rates

    
  contains

    !Simple photolysis rate function used by some photolysis reactions
    !This is copied from $KPP_root_interface module.
    !TODO: Eventually one may remove this function from those modules!
    real function photo(l, m, n, cos_theta)
      implicit none
      real, intent(in) :: l, n, m, cos_theta
      if (cos_theta < 1e-5) then
        photo = 0.0
        return
      end if
      photo = l * cos_theta**(m) * exp(-n / cos_theta)
    end function photo


    subroutine test()
      implicit none
      real :: a
      integer :: k

      call get_interp_weights(0.5, (/1.,2.,4./), k, a)
      print *, 'Should be 1.0, 1:', a, k
      call get_interp_weights(3.4, (/1.,2.,4./), k, a)
      print *, 'Should be 0.3, 2:', a, k
      call get_interp_weights(4.0, (/1.,2.,4./), k, a)
      print *, 'Should be 2, 0.0:', a, k

      call get_interp_weights_rev(4.0, (/4.,2.,1./), k, a)
      print *, 'Should be 1, 1.0:', a, k
      call get_interp_weights_rev(3.0, (/4.,2.,1./), k, a)
      print *, 'Should be 1, 0.5:', a, k
      call get_interp_weights_rev(.05, (/4.,2.,1./), k, a)
      print *, 'Should be 2, 0.0:', a, k

      
    end subroutine test

    integer function month_to_season(month) result(ind)
      implicit none
      integer, intent(in) :: month

      select case(month)
      case (3,4,5,6,7,8)
        ind = ind_summer
      case (1,2,9,10,11,12)
        ind = ind_winter
      end select
      !!!! With line below uncommented scores for ozone improve
      !ind = ind_winter
    end function month_to_season
  end subroutine get_photorates_column
  
  !************************************************************************************

  subroutine effective_albedo_cld(col_cwcabove, cwc_totcol, alb_sfc, &
       & cosza, num_levs, cld_att, alb_eff, cloud_model, tau_above_bott, ssa)

    implicit none
    real, dimension(:), intent(in) :: col_cwcabove
    real, intent(in) :: cwc_totcol, alb_sfc, cosza
    integer, intent(in) :: num_levs
    integer, intent(in) :: cloud_model
    real, dimension(max_levels) :: tau_above_bott
    real :: ssa
    real, dimension(:), intent(out) :: cld_att
    real, intent(out) :: alb_eff

    integer :: iLev

    !! From "Optical properties of terrestrial clouds" by Kokhanovsky                                                                        
    !real, parameter :: dp = 10e-6 !! Effective droplet diameter for water and ice clouds combined                                           
    real, parameter :: dp = 6e-6 !! Droplet diameter for liquid water clouds                                                                 
                                 ! !(ok value based on e.g. Kokhanovsky)                                                                     
    
    real, parameter :: cwc_to_tau = 3/(2*1000*dp)  !geometrical extinction per unit of liquid cwc     
            !Kokhanovsky eq. 3.9                                                                      
    real, parameter :: g = 0.8 !asymmetry , single-scattering albedo                                  
          !!A. Kokhanovsky / Earth-Science Reviews 64 (2004) 189-241                                  
          !! doi:10.1016/S0012-8252(03)00042-4                                                        
          !!  0.75 ice, 0.85 water

    if (cloud_model == simple_cloud) then
      tau_above_bott(1) = cwc_totcol*cwc_to_tau                                                         
      do iLev=2,num_levs                                                                               
        tau_above_bott(ilev) = 0.5*cwc_to_tau*(col_cwcabove(iLev-1)+col_cwcabove(iLev))                 
      enddo
      ssa = 0.9999
    end if

    call effective_albedo_v2(ssa, g, tau_above_bott, alb_sfc, cosza, num_levs, cld_att, alb_eff)

  end subroutine effective_albedo_cld
  
  !************************************************************************************

  subroutine effective_albedo_aer(col_ext, col_scat, alb_sfc, cosza, num_levs, aer_att, alb_eff)
    !
    ! Simple aerosol attenuation scheme
    ! aerosol properties are considered uniform over vertical
    !

    implicit none
    real, dimension(:), intent(in) :: col_ext, col_scat
    real, intent(in) :: alb_sfc, cosza
    integer, intent(in) :: num_levs
    
    real, dimension(:), intent(out) :: aer_att
    real, intent(out) :: alb_eff
    real, dimension(max_levels) :: tau_above_bott
    integer :: iLev

    real :: ssa
    real, parameter :: g = 0.6  !!Asymmetry
    ! Lose median from 
    ! Andrews, E., et al. (2006), Comparison of methods for deriving aerosol asymmetry parameter, J. Geophys. Res., 111,
    !D05S04, doi:10.1029/2004JD005734.
    
    !!Model breaks at SSA=1.
    ssa  = 0.999 * sum(col_scat(1:num_levs))/sum(col_ext(1:num_levs))  !! ssa mean over column

    !make cumulative aod above bottom
    tau_above_bott(num_levs+1) = 0.
    do iLev=num_levs,1,-1
      tau_above_bott(iLev) = tau_above_bott(iLev+1)+col_ext(iLev)
    enddo
    call effective_albedo_v2(ssa, g, tau_above_bott, alb_sfc, cosza, num_levs, aer_att, alb_eff)


  end subroutine effective_albedo_aer

            
  !************************************************************************************

  subroutine effective_albedo_v2(ssa, g, tau_above_bott, alb_sfc, cosza_in, num_levs, att_prof, alb_eff)
    !
    ! Full 2-stream implementation
    !
    !
    ! Simple cloud attanuation scheme
    !

    implicit none
    real, intent(in) :: ssa, g, alb_sfc, cosza_in
    real, dimension(:), intent(in) :: tau_above_bott
    integer, intent(in) :: num_levs
    
    real, dimension(:), intent(out) :: att_prof
    real, intent(out) :: alb_eff

    integer :: iLev
    real(8) :: f, ssap, taufactor, gprime ! parameters for delta-function adjustment
    real(8) :: H, K, u, v, eps, gama, epsgamafact  !! Coeffs from (Liu 6.5.29ab)  
    real(8) gama1, gama2, gama3, fk !!coeffs 6.5.30ab
    real(8) :: att0, attdir, d, dp, taup0, taup, vpuK, vpuH, cosza
    real(8) , parameter :: maxtau = 1000.

    cosza = max(0.01,cosza_in) !! small negatives can come here, but cause
      !    trouble for two-stream treatment of a layer
    !
    ! delta-func correcton
    f = g*g !Fraction of energy in forward peak
    taufactor = 1-f*ssa ! 6.5.31a
    ssap = (1-f)*ssa/(1-f*ssa)
    gprime = (g-f)/(1-f)


    !! Eeddington: !!! Table 6.2 in Liou 2002
    gama1 =   0.25 * ( 7. - (4.+3.*gprime)*ssap )
    gama2 = - 0.25 * ( 1. - (4.-3.*gprime)*ssap )
    gama3 =   0.25 * ( 2. - 3.*ssap*cosza )
    fk = sqrt(gama1*gama1 - gama2*gama2) !!Liou (6.5.30ab)
    v = 0.5*(1+(gama1-gama2)/fk)
    u = 0.5*(1-(gama1-gama2)/fk)

    epsgamafact = cosza*ssap / (1. - cosza*cosza*fk*fk)  ! Common part of eps and gama
    eps  = ( gama3*(1.-gama1*cosza) - gama2*cosza*(1.-gama3) ) * epsgamafact
    gama = -1.*( (1 - gama3)*(1+gama1*cosza) +  gama2*gama3*cosza)* epsgamafact
      
    taup0 = min(tau_above_bott(1)*taufactor, maxtau)

    if (fk*taup0 > 1e-3) then
      d = exp(fk*taup0)
      att0 = exp(-taup0/cosza)

      alb_eff =  -(((d**2-att0*d)*gama*u**2+(1-d**2)*eps*u*v+(att0*d-1)*gama*v**2+(att0*d*cosza*v**2-att0*d*cosza*u**2))*alb_sfc+(d**2-att0*d)*eps*v**2+(1-d**2)*gama*u*v+(att0*d-1)*eps*u**2)/((cosza*u**2-d**2*cosza*v**2)+(d**2-1)*cosza*u*v*alb_sfc)

      H =  (((att0*d-d**2)*gama*u+att0*d*cosza*u)*alb_sfc+d**2*gama*v-att0*d*eps*u)/(u**2-d**2*v**2+(d**2-1)*u*v*alb_sfc)
      K =  -(((att0*d-1)*gama*v+att0*d*cosza*v)*alb_sfc-att0*d*eps*v+gama*u)/(u**2-d**2*v**2+(d**2-1)*u*v*alb_sfc)

      vpuK = (v+u)*K
      vpuH = (v+u)*H
      
      do iLev=1,num_levs
        taup = min(tau_above_bott(iLev)*taufactor, maxtau)
        dp = exp(fk*taup)
        attdir = exp(-taup/cosza)  !!Direct beam attenuation
        att_prof(iLev) = (vpuK*dp + vpuH/dp + (eps+gama+1)*attdir) / (1+cosza*alb_eff)
        !!Factor to the profile with given albedo_eff at surface
      enddo
!!!FIXME DBG
!      if (.not. all(att_prof(1:num_levs ) .ge. 0.)) then
!      !$OMP CRITICAL (bark_effective_albedo_v2)
!      call msg("taufactor, fk, cosza, vpuK, vpuH, eps, gama, alb_eff, alb_sfc", (/taufactor, fk, cosza, vpuK, vpuH, eps, gama, real(alb_eff,kind=8), real(alb_sfc,kind=8)/))
!      call msg("u,v,k,h,taup0,d,att0",(/u,v,k,h,taup0,d,att0/))
!      do iLev=1,num_levs
!        taup = tau_above_bott(iLev)*taufactor
!        dp = exp(fk*taup)
!        attdir = exp(-taup/cosza)  !!Direct beam attenuation
!        att_prof(iLev) = (vpuK*dp + vpuH/dp + (eps+gama+1)*attdir) / (1+cosza*alb_eff)
!        call msg("lev:"//trim(fu_str(iLev))//" tau_above_bott att_prof  taup dp ", (/tau_above_bott(ilev), att_prof(ilev), real(taup), real(dp)/))
!        !!Factor to the profile with given albedo_eff at surface
!      enddo
!      !$OMP END CRITICAL (bark_effective_albedo_v2)
        
!      endif

    else
       alb_eff = alb_sfc
       att_prof(1:num_levs) = 1.
    endif

  end subroutine effective_albedo_v2
  
  !***********************************************************************************


  subroutine compute_water_cloud_properties(masses, metdat_col, col_ext, col_scat, density, diam, &
       & soot_col, ix, iy, pwc_col, tau_above_bott, ssa)

    ! Used for both the detailed cloud model for photolysis and the precipitable water
    ! for 2020_scavenging
    ! TRANSFORMATION_PARAMETERS:
    ! For photolysis:
    ! cloud_model_for_photolysis = DETAILED_CLOUD
    ! For precipitable water:
    ! wet_deposition_scheme = 2020_SCAVENGING

    real, dimension(:,:,:), intent(in) :: masses !nSp, nSrc, nLev  slice of a massmap                                      
    real, dimension(:,:), intent(in) :: metdat_col
    real, dimension(:), intent(in) :: col_ext, col_scat, density, diam, soot_col
    integer, intent(in) :: ix, iy

    real, intent(out) :: ssa !, tcc_photo, tcc_scav
    real, dimension(:), intent(out):: pwc_col, tau_above_bott

    !real, dimension(max_levels) :: col_cwcabove, col_lcwcabove, col_cwcabove2
    real :: lwc_d, cdnc, B, T, iwc_normalized, dp_dyn, tau, ssa0, aerosol_nbr_cnc
    real :: norm, tcc, area
    integer :: iLev, iSpecies, nSpecies, num_levs
    integer, parameter :: mu = 3 ! the asymmetry parameter of the gamma distribution
    real, parameter :: m_liq = mu*14., m_ice = mu*1800. ! Integration limits for the formation of
                                                       ! precipitable water in 2019_scav
                                                       ! The numerical value is in micrometers

    real, dimension(:), allocatable :: clw_dens_col, ciw_dens_col, lwc_col, iwc_col, ccover, cwc, cic

    tau_above_bott(:) = 0.0
    pwc_col(:) = 0.0

    if (maxval(metdat_col(imet_cwc,:) + metdat_col(imet_cic,:)) < 1e-25) then
      ssa = 0.9999
      return
    end if

    nSpecies = size(masses,1)
    num_levs = size(masses,3)

    area = metdat_col(imet_dx_size,1) * metdat_col(imet_dy_size,1)

    allocate(clw_dens_col(num_levs), ciw_dens_col(num_levs), lwc_col(num_levs), &
         & iwc_col(num_levs), ccover(num_levs), cwc(num_levs), cic(num_levs))

    ccover(:) = metdat_col(imet_cc3d, 1:num_levs)

    where (ccover < 0.01) ccover = 0.01
    where (ccover > 0.99) ccover = 0.99

    cwc(1:num_levs) = metdat_col(imet_cwc,1:num_levs)
    cic(1:num_levs) = metdat_col(imet_cic,1:num_levs)

    where (cwc(:) < 1e-10) cwc = 1e-10
    where (cic(:) < 1e-10) cic = 1e-10

    clw_dens_col(:) = cwc(1:num_levs)*metdat_col(imet_airdens,1:num_levs)/ccover(1:num_levs)
    ciw_dens_col(:) = cic(1:num_levs)*metdat_col(imet_airdens,1:num_levs)/ccover(1:num_levs)

    lwc_col(:) = cwc(1:num_levs)*metdat_col(imet_airmass,1:num_levs)/area/ccover(1:num_levs)
    iwc_col(:) = cic(1:num_levs)*metdat_col(imet_airmass,1:num_levs)/area/ccover(1:num_levs)

    ssa = 0.0
    norm = 0.0

    do iLev=num_levs,1,-1      

      lwc_d = 1e-3*clw_dens_col(iLev)

      if (lwc_d > 0) then
        ! Sum over aerosols only: diam>0; real_missing**3 is -Inf in gfortran, 
        ! leads to formal zero, but triggers floating-point exception
        aerosol_nbr_cnc = sum(&
            &  masses(:,1,iLev)/(density(1:nSpecies)*1.33*3.14*(0.5*diam(1:nSpecies))**3),&
            & MASK = (diam(1:nSpecies) > 0)&
            &) *  metdat_col(imet_airdens,iLev)/metdat_col(imet_airmass,iLev)

        if (aerosol_nbr_cnc > 9e7) then
          cdnc = max(-765.5 + 395.71 * log10(aerosol_nbr_cnc*1e-6), 15.0)
        else
          cdnc = 15.0
        end if

        dp_dyn = 7*(1e9*lwc_d/cdnc)**0.33
        dp_dyn = max(dp_dyn, 4.0)

        pwc_col(iLev) = lwc_col(iLev) * exp(-m_liq/dp_dyn) * (1 + m_liq/dp_dyn + &
             & (1/2.)*(m_liq/dp_dyn)**2 + (1/6.)*(m_liq/dp_dyn)**3)

        dp_dyn = min(dp_dyn, 15.)
        
        tau = 3/(2*1000*dp_dyn*1e-6) * lwc_col(iLev)
      else
        tau = 0.0
      end if
      
      ! Ice water contribution                                                                                             
      ! "The effective radius in Ice Clouds", Klaus Wyser                                                                  
      ! https://doi.org/10.1175/1520-0442(1998)011%3C1793:TERIIC%3E2.0.CO;2                                                
      iwc_normalized = ciw_dens_col(iLev)/0.05
      if (iwc_normalized > 0.0) then
        T = min(metdat_col(imet_temp, iLev), 272.99)
        B = -2 + 1e-3 * (273.-T)**1.5 * log10(iwc_normalized)
        dp_dyn = 377.4 + 203.3*B + 37.91*B**2 + 2.3696*B**3
        dp_dyn = max(dp_dyn, 5.0)

        pwc_col(iLev) = pwc_col(iLev) + iwc_col(iLev)*exp(-m_ice/dp_dyn) * (1 + m_ice/dp_dyn + &
             & (1/2.)*(m_ice/dp_dyn)**2 + (1/6.)*(m_ice/dp_dyn)**3)
        tau = tau + 3/(2*910*dp_dyn*1e-6) * iwc_col(iLev)     
      end if

      if (tau < col_ext(iLev))  then
        tau = col_ext(iLev)
        ssa0 = 0.999 * col_scat(iLev)/col_ext(iLev)
      else
        ssa0 = max(0.8, 0.9999 - (0.005*soot_col(iLev)/(iwc_col(iLev)+lwc_col(iLev)))/8e-6)
      endif

      !Computation of the average ssa weighted by tau                                                                      
      ssa = ssa + ssa0*tau !*ccover
      norm = norm + tau !*ccover

      if (iLev < num_levs) then 
        tau_above_bott(iLev) = tau + tau_above_bott(iLev+1)
      else
        tau_above_bott(iLev) = tau
      end if

      !$OMP CRITICAL
      if (pwc_col(iLev) < 0.0) then
        call msg_warning('Negative precipitable water!', 'compute_water_cloud_properties')
        call msg('pwc_col(iLev), lwc_col(iLev), iwc_col(iLev)', &
             & (/pwc_col(iLev), lwc_col(iLev), iwc_col(iLev)/))
        call msg('metdat_col(imet_airdens,iLev), metdat_col(imet_airmass,iLev)', &
             & metdat_col(imet_airdens,iLev), metdat_col(imet_airmass,iLev))
        call msg('metdat_col(imet_cwc,iLev), metdat_col(imet_cic,iLev)', &
             & metdat_col(imet_cwc,iLev), metdat_col(imet_cic,iLev))
        call msg('aerosol_nbr_cnc', aerosol_nbr_cnc)
        call msg('lwc_d', lwc_d)
        call msg('iwc_normalized', iwc_normalized)
        call msg('iLev', iLev)
        call set_error('Negative precipitable water!', 'compute_water_cloud_properties')
      end if
      !$OMP END CRITICAL
    end do

    ssa = ssa/norm

    pwc_col(:) = pwc_col(:)*ccover(:)

  end subroutine compute_water_cloud_properties


  subroutine get_interp_weights(value_req, values, ind_low, weight_low)
    implicit none
    real, intent(in) :: value_req
    real, dimension(:), intent(in) :: values
    integer, intent(out) :: ind_low
    real, intent(out) :: weight_low

    real :: value_avail
  
    if (value_req < values(1)) then
      ind_low = 1
      weight_low = 1.0
    else if (value_req > values(size(values))) then
      ind_low = size(values) - 1
      weight_low = 0.0
    else
      do ind_low = 1, size(values)-2
        if (value_req <= values(ind_low+1)) exit
      end do
      weight_low = (values(ind_low+1) - value_req) / (values(ind_low+1) - values(ind_low))
    end if

  end subroutine get_interp_weights

  !************************************************************************************

  subroutine photolysis_input_needs(ifDynamicAlbedo, cloud_model,  meteo_input_local)
    !
    ! Fills meteo input for in-transformation 
    !
    implicit none

    ! Imported parameters
    logical, intent(in) :: ifDynamicAlbedo
    integer, intent(in) :: cloud_model
    type(Tmeteo_input),  intent(out), target :: meteo_input_local
    character(len=*), parameter :: subname="photolysis_input_needs"

    ! Local variables
    integer :: iQ, iTmp, nq
    
    meteo_input_local = meteo_input_empty
    meteo_input_local%quantity(1:2) = (/cell_size_x_flag, cell_size_y_flag/)
    meteo_input_local%q_type(1:2) = dispersion_single_time_flag
    imet_dx_size =>  meteo_input_local%idx(1)
    imet_dy_size =>  meteo_input_local%idx(2)

    nq = 2

    nq = nq +1
    meteo_input_local%quantity(nq) = latitude_flag
    meteo_input_local%q_type(nq) = meteo_single_time_flag
    imet_lat => meteo_input_local%idx(nq) 

    nq = nq +1
    meteo_input_local%quantity(nq) = pressure_flag
    meteo_input_local%q_type(nq) = meteo_dynamic_flag
    imet_press => meteo_input_local%idx(nq) 

    imet_cc3d => null()
    imet_cwc3d => null()
    imet_cwcol => null()
    imet_dz_size => null()
    imet_lcwc3d => null()
    imet_lcwcol => null()
    imet_temp => null()

    if (cloud_model /= fake_cloud) then
      ! Common for non-fake cloud

      nq = nq +1
      meteo_input_local%quantity(nq) = total_cloud_cover_flag
      meteo_input_local%q_type(nq) = meteo_dynamic_flag
      imet_tcc => meteo_input_local%idx(nq)

      nq = nq + 1
      meteo_input_local%quantity(nq) = cwcabove_3d_flag
      meteo_input_local%q_type(nq) = meteo_dynamic_flag
      imet_cwc3d => meteo_input_local%idx(nq)

      nq = nq +1
      meteo_input_local%quantity(nq) = cloud_ice_flag
      meteo_input_local%q_type(nq) = meteo_dynamic_flag
      imet_cic => meteo_input_local%idx(nq)

      if (cloud_model == simple_cloud) then
        nq = nq +1
        meteo_input_local%quantity(nq) = cwcolumn_flag
        meteo_input_local%q_type(nq) = meteo_dynamic_flag
        imet_cwcol => meteo_input_local%idx(nq) 

      
        nq = nq + 1
        meteo_input_local%quantity(nq) = lcwcabove_3d_flag
        meteo_input_local%q_type(nq) = meteo_dynamic_flag
        imet_lcwc3d => meteo_input_local%idx(nq)

        nq = nq + 1
        meteo_input_local%quantity(nq) = lcwcolumn_flag
        meteo_input_local%q_type(nq) = meteo_dynamic_flag

        imet_lcwcol => meteo_input_local%idx(nq)
      elseif (cloud_model == detailed_cloud) then

        nq = nq +1
        meteo_input_local%quantity(nq) = cloud_cover_flag
        meteo_input_local%q_type(nq) = meteo_dynamic_flag
        imet_cc3d => meteo_input_local%idx(nq)

        nq = nq +1
        meteo_input_local%quantity(nq) = temperature_flag
        meteo_input_local%q_type(nq) = meteo_dynamic_flag
        imet_temp => meteo_input_local%idx(nq)

        nq = nq +1
        meteo_input_local%quantity(nq) = cell_size_z_flag
        meteo_input_local%q_type(nq) = dispersion_dynamic_flag
        imet_dz_size => meteo_input_local%idx(nq)
        
        nq = nq +1
        meteo_input_local%quantity(nq) = disp_cell_airmass_flag
        meteo_input_local%q_type(nq) = dispersion_dynamic_flag
        imet_airmass => meteo_input_local%idx(nq)

        nq = nq +1
        meteo_input_local%quantity(nq) = air_density_flag
        meteo_input_local%q_type(nq) = dispersion_dynamic_flag
        imet_airdens => meteo_input_local%idx(nq)

        nq = nq +1
        meteo_input_local%quantity(nq) = cloud_water_flag
        meteo_input_local%q_type(nq) = meteo_dynamic_flag
        imet_cwc => meteo_input_local%idx(nq)


      end if
    endif

    imet_albedo => null()
    if (ifDynamicAlbedo) then 
      nq = nq + 1
      meteo_input_local%quantity(nq) = albedo_flag
      meteo_input_local%q_type(nq) = meteo_dynamic_flag
      imet_albedo => meteo_input_local%idx(nq) 
    else
      imet_albedo => null()
    endif

    meteo_input_local%nQuantities = nq

  end subroutine photolysis_input_needs

  !***********************************************************************

  subroutine get_interp_weights_rev(value_req, values, ind_low, weight_low)
    implicit none
    real, intent(in) :: value_req
    real, dimension(:), intent(in) :: values
    integer, intent(out) :: ind_low
    real, intent(out) :: weight_low

    if (value_req < values(size(values))) then
      ind_low = size(values)-1
      weight_low = 0.0
    else if (value_req > values(1)) then
      ind_low = 1
      weight_low = 1.0
    else
      do ind_low = size(values)-1, 2, -1
        if (value_req <= values(ind_low)) exit
      end do
      weight_low = (values(ind_low+1) - value_req) / (values(ind_low+1) - values(ind_low))
    end if

  end subroutine get_interp_weights_rev

!!!!  subroutine test_lut_rates()
!!!!    implicit none
!!!!    
!!!!    real, dimension(1) :: col_press = (/1e5/) ! Pa!!
!!!!    real, dimension(1) :: col_cwc = (/1e-3/) ! kg/m3!!
!!!!    real, dimension(1) :: col_aod = (/1e-3/) ! m2/m2
!!!!    real, dimension(2, 1) :: rates
!!!!    real :: cwc_totcol = 1e-3, alb_sfc = 0.3, cld_cover = 0.5
!!!!    type(silja_time) :: now
!!!!
!!!!            
!!!!    now = fu_set_time_utc(2001, 7, 1, 0, 0, 0.0)
!!!!    call report_lut(1, 1)
!!!!    call get_photorates_column(0.0, 0.0, now, 1.0, col_press, col_cwc, cwc_totcol, col_aod,  alb_sfc, cld_cover, rates)
!!!!    print *, 'rates 1000 hPa'
!!!!    print *, rates
!!!!
!!!!    ! the second sza bin in the test lut, cos(theta) = 0.76
!!!!    call get_photorates_column(0.0, 0.0, now, 0.76, col_press, col_cwc, cwc_totcol, col_aod,  alb_sfc, cld_cover, rates)
!!!!    print *, 'rates 1000 hPa, 0.76'
!!!!    print *, rates
!!!!
!!!!    call get_photorates_column(0.0, 0.0, now, 0.9397, col_press, col_cwc, cwc_totcol, col_aod,  alb_sfc, cld_cover, rates)
!!!!    print *, 'rates 1000 hPa, 0.94'
!!!!    print *, rates
!!!!
!!!!
!!!!    col_press(1) = 0.75e5
!!!!    call get_photorates_column(0.0, 0.0, now, 1.0, col_press, col_cwc, cwc_totcol, col_aod,  alb_sfc, cld_cover, rates)
!!!!    print *, 'rates 750 hPa'
!!!!    print *, rates
!!!!
!!!!    col_press(1) = 0.75e5
!!!!    call get_photorates_column(0.0, 0.0, now, 0.76, col_press, col_cwc, cwc_totcol, col_aod,  alb_sfc, cld_cover, rates)
!!!!    print *, 'rates 750 hPa, 0.76'
!!!!    print *, rates ! 316.22775      3.16227768E+10
!!!!
!!!!    col_press(1) = 0.75e5
!!!!    call get_photorates_column(0.0, 0.0, now, 0.9397, col_press, col_cwc, cwc_totcol, col_aod,  alb_sfc, cld_cover, rates)
!!!!    print *, 'rates 750 hPa, 0.94' 
!!!!    print *, rates ! 31.606142      3.16061158E+09
!!!!
!!!!  contains
!!!!    subroutine report_lut(ind_lat, ind_season)
!!!!      implicit none
!!!!      integer, intent(in) :: ind_lat, ind_season
!!!!      integer :: ii, jj,kk
!!!!
!!!!      print *, 'Reporting lut data:'
!!!!           
!!!!      do ii = 1, size(lut_data, 1)
!!!!        print *, 'Reaction:', ii
!!!!        do jj = 1, size(lut_data, 5)
!!!!          print *, 'O3:', jj
!!!!          do kk = 1, size(lut_data, 6)
!!!!            print *, 'sza:', kk
!!!!            print *, lut_data(ii, ind_lat, ind_season, 1, jj, kk, :)
!!!!          end do
!!!!        end do
!!!!      end do
!!!!
!!!!      print *, ''
!!!!      
!!!!    end subroutine report_lut
!!!!
!!!!  end subroutine test_lut_rates

  subroutine test_lut_rates_full()
    implicit none
    

    integer, parameter :: nLev = 3, nMet = 5

    integer, dimension(9), target :: metindex
    real, dimension(9,3) :: met_dat
    real, dimension(num_reactions, 3) :: rates
    real, dimension(3), parameter :: aod_ext =  (/1. , 1., 1./), &
                                  & aod_scat = (/ .9,  .9, .9/), &
                                  & thickness = (/100., 100., 100./)
    real, dimension(3), parameter :: o3_col =  (/346. , 300., 2.8e-2/) ! in DU (values approximate for nm-sum) CHECK!

    type(silja_time) :: now 
    integer :: ind_lev, ind_react
    
    metindex(:) = (/1,2,3,4,5,6,7,8,9/)
    imet_albedo => metindex(1)
    imet_cwc3d  => metindex(2) 
    imet_press  => metindex(3) 
    imet_cwcol  => metindex(4) 
    imet_lat    => metindex(5) 
    imet_tcc    => metindex(6)
    imet_temp    => metindex(7)
    imet_lcwc3d  => metindex(8)
    imet_lcwcol  => metindex(9)

     met_dat(imet_albedo,:) = 0.3
     met_dat(imet_cwc3d,:) =  (/3e-3, 1e-3, 0.0/) ! kg/m3!!
     met_dat(imet_press,:) = (/1e5, 2500.0, 280.0/) ! Pa!!
     met_dat(imet_cwcol,:) =  3e-3
     met_dat(imet_lat,:) =    0. 
     met_dat(imet_tcc,:) =    0.3    
     met_dat(imet_temp,:) =    (/273, 263, 253/)
     met_dat(imet_lcwc3d,:) =  (/3e-3, 1e-3, 0.0/) ! kg/m3!!                                                 
     met_dat(imet_lcwcol,:) =  3e-3
    
    now = fu_set_time_utc(2001, 7, 1, 0, 0, 0.0)
    
    ! equator, time as above, pressure as in col_press, sun at zenith.
    !call get_photorates_column(met_dat, 1., real_missing, now, aod_ext, &
    !     & aod_scat, o3_col, rates, .true., .true., thickness, version)
    do ind_lev = 1, nMet
      call msg('Level, pressure:', ind_lev, met_dat(imet_press,ind_lev))
      do ind_react = 1, num_reactions
        call msg('Reaction, rate', ind_react, rates(ind_react, ind_lev))
      end do
    end do

    ! equator, time as above, pressure as in col_press, 30 degree zenith angle.
    !call get_photorates_column(met_dat, cos(30./180*pi), real_missing, now, aod_ext, &
    !     & aod_scat, o3_col, rates, .true., .true., thickness)
    do ind_lev = 1, nMet
      call msg('Level, pressure:', ind_lev, met_dat(imet_press,ind_lev))
      do ind_react = 1, num_reactions
        call msg('Reaction, rate', ind_react, rates(ind_react, ind_lev))
      end do
    end do

    ! equator, time as above, pressure as in col_press, 91 degree zenith angle.
    !call get_photorates_column(met_dat,  cos(91./180*pi), real_missing, now, aod_ext, &
    !     & aod_scat, o3_col, rates, .true., .true., thickness)
    call msg('Twilight:')
    do ind_lev = 1, nMet
      call msg('Level, pressure:', ind_lev, met_dat(imet_press,ind_lev))
      do ind_react = 1, num_reactions
        call msg('Reaction, rate', ind_react, rates(ind_react, ind_lev))
      end do
    end do


  end subroutine test_lut_rates_full


  subroutine test_lut_init(lut_file_name)
    implicit none
    character(len=*), intent(in) :: lut_file_name
    type(Tsilam_namelist), pointer :: nl_fake
    integer :: num_items
    real, parameter :: molec_to_du = 1/2.69e16 ! du = 2.69e16 molecules/cm2
    
    call init_photolysis_lut(lut_file_name)
    if (error) return
    print *, 'lut_sza_rad:', lut_sza_rad
    print *, 'lut_o3', lut_o3
    print *, 'lut_alb', lut_alb
    !print *, 'atm_o3_cuml, nm, win', atm_o3_cuml(2, 2, :)*molec_to_du
    print *, 'atm_o3_cuml, sm, win', atm_o3_cuml(2, 2, :)*molec_to_du !RISTO: after re-reordering the indeces
    
  end subroutine test_lut_init

  subroutine test_effective_albedo_v2()
      integer, parameter :: nlevs=10
      real, parameter :: ssa=0.999, g=0.78  !!!aerosol/drop parameters

      real, dimension(nlevs) :: tau_above_bott, photoatt
      real :: alb_eff

      real, parameter :: tautot=30, alb_sfc=0.1, mu0=1.0


      integer :: i

      do i=1,nlevs
        tau_above_bott(i)=tautot*(nlevs-i+1)*1.0/nlevs
      enddo
      call effective_albedo_v2(ssa, g, tau_above_bott, alb_sfc, mu0, nlevs, photoatt, alb_eff)

      print *, "tautot=", tautot, "alb_sfc=",alb_sfc, "mu0=", mu0, "albeff=", alb_eff
      print *, "tau_above_bott(i), photoatt(i)"
      do i=1,nlevs
        print *, tau_above_bott(i), photoatt(i)
      enddo

  endsubroutine test_effective_albedo_v2





end module photolysis
