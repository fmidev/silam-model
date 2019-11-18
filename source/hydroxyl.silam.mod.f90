!
! Lookup-table based parameterization for OH concentration.
! 
! This module provides a function for evaluating OH concentration (mol/m3) from surface up
! to ~20 km. For troposphere, a zonal-mean climatology is used. Above the climatology top,
! the the zenith angle fit of Hanisco et al. (2001, J Phys Chem) is used. The fit is valid
! for ~17-19 km; the code allows extrapolation up to 25 km. Above, the OH concentration
! begins to increase rapidly (see eg. Pickett and Peterson (1996, JGR)).
! 
! The climatology is from Spivakovsky et al. (2000, JGR). The monthly means are
! disaggregated temporally using the rate of O3-O1D reaction as proxy. The rate is
! evaluated with the parameterization of MCM (Sanders 2003).
!

module hydroxyl
  use silam_levels
  use silam_times
  implicit none
  private

  public setup_oh_lut
  public fu_oh_cnc_clim
  public interp1d

  integer :: num_levs_in, num_lats_in, num_months_in
  real, dimension(:), allocatable :: lut_press_in
  real, dimension(:), allocatable :: lut_lats_in
  integer, dimension(:), allocatable :: lut_jdays, lut_months
  real, dimension(:,:,:), allocatable :: lut_data_in

  real, dimension(:), allocatable :: lut_alt, lut_lat
  real, dimension(:,:,:), allocatable :: lut_data, lut_norm ! lat, lev, month
  ! lut lats = -90,-89...89,90
  integer, private, save :: num_lats
  real, private, parameter :: lat_incr = 1.0
  integer, private, save :: num_levs ! equal to nz_dispersion
  logical, private, save :: initialized = .false., allow_strato
  integer, private, parameter :: num_lut_months = 12

contains
  
  subroutine destroy_lut()
    implicit none
    call msg('Deallocating OH LUT data...')
    deallocate(lut_press_in, lut_lats_in, lut_jdays, lut_data_in, lut_alt, lut_data, lut_lat, lut_months)
    num_levs_in = int_missing 
    num_lats_in = int_missing
    num_levs = int_missing
    num_lats = int_missing
    initialized = .false.
  end subroutine destroy_lut

  subroutine read_lut(filename)
    ! Read the lookup table data (OH concentration, mol/m3) from a namelist file and
    ! perform a few checks. The data can be either 4-monthly or monthly.
    implicit none
    character(len=*), intent(in) :: filename

    integer :: file_unit
    character(len=*), parameter :: sub_name = 'read_lut'
    type(Tsilam_namelist), pointer :: nlptr
    integer :: num_items, ind_item, ind_lat, ind_month
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: p_items

    logical, dimension(:,:), allocatable :: items_ok
    integer :: stat
    character(len=fnlen) :: content

    call msg('Reading OH LUT from '//trim(filename))

    file_unit = fu_next_free_unit()
    open(file_unit, file=fu_process_filepath(filename), iostat=stat, status='old', action='read')
    if (fu_fails(stat == 0, 'Failed to open: ' // trim(filename), sub_name)) return
    
    nlptr => fu_read_namelist(file_unit, .false.)
    close(file_unit)
    if (error) return
    if (fu_fails(associated(nlptr), 'nlptr not associated', sub_name)) return

    num_levs_in = fu_content_int(nlptr, 'number_of_levels')
    if (fu_fails(num_levs_in /= int_missing .and. num_levs_in > 0, 'Bad number_of_levels', sub_name)) return
    num_lats_in = fu_content_int(nlptr, 'number_of_latitudes')
    if (fu_fails(num_lats_in /= int_missing .and. num_lats_in > 0, 'Bad number_of_latitudes', sub_name)) return
    num_months_in = fu_content_int(nlptr, 'number_of_months')
    if (fu_fails(num_months_in /= int_missing .and. num_months_in > 0, 'Bad number_of_months', sub_name)) return

    !if (fu_fails(num_months_in == 4, 'Number of months must be 4', sub_name)) return
    allocate(lut_data_in(num_lats_in, num_levs_in, num_months_in), lut_press_in(num_levs_in), &
           & lut_lats_in(num_lats_in), lut_jdays(num_months_in), items_ok(num_lats_in, num_months_in), &
           & lut_months(num_months_in), stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return

    select case(num_months_in)
    case (4)
      lut_months = (/1, 4, 7, 10/)
    case (12)
      lut_months = (/1,2,3,4,5,6,7,8,9,10,11,12/)
    case default
      call msg('Number of months:', num_months_in)
      call set_error('Number of months not supported', 'read_lut')
      return
    end select
    do ind_month = 1, num_months_in
      lut_jdays(ind_month) = fu_julian_date(fu_set_time_utc(2001,lut_months(ind_month),15, 0, 0, 0.0))
    end do
    
    content = fu_content(nlptr, 'pressure')
    if (fu_fails(content /= '', 'Missing pressure entry', sub_name)) return
    read(content, fmt=*) lut_press_in
    content = fu_content(nlptr, 'latitude')
    if (fu_fails(content /= '', 'Missing latitude entry', sub_name)) return
    read(content, fmt=*) lut_lats_in
    
    nullify(p_items)
    call get_items(nlptr, 'val', p_items, num_items)
    if (fu_fails(associated(p_items) .and. num_items > 0, 'Failed to get items', sub_name)) return

    items_ok = .false.
    do ind_item = 1, num_items
      content = fu_content(p_items(ind_item))
      read(content, fmt=*, iostat=stat) ind_lat, ind_month, lut_data_in(ind_lat, :, ind_month)
      if (fu_fails(stat == 0, 'Failed to parse: ' // trim(content), sub_name)) return
      if (fu_fails(.not. items_ok(ind_lat,ind_month), 'Something wrong reading', sub_name)) return
      items_ok(ind_lat,ind_month) = .true.
    end do
    if (fu_fails(all(items_ok), 'Not all values read', sub_name)) return
    if (fu_fails(all(lut_data_in >= 0.0), 'Negative lut values read', sub_name)) return
    deallocate(items_ok)
    
  end subroutine read_lut
  
  subroutine setup_oh_lut(lut_filename, vertical)
    ! Setup the lookup table: 
    ! - interpolate to monthly level (if needed)
    ! - interpolate to dispersion vertical
    ! - interpolate to a finer latitude grid (determined by lat_incr).
    ! 
    implicit none
    character(len=*), intent(in) :: lut_filename
    type(silam_vertical), intent(in) :: vertical

    integer :: ind_lev, ind_lat, ind_n, ind_month
    real, dimension(:,:,:), allocatable :: lut_tmp, lut_monthly
    type(silam_vertical) :: vert_lut, vert_alt, vert_tmp
    real :: flev, lat, weight_north, weight_up, lat_lut_n, lat_lut_s, lev_hgt, w1, w2, div
    integer :: stat, ind_item, jday, i1, i2
    character(len=fnlen) :: content
    character(len=*), parameter :: sub_name = 'setup_oh_lut'
    real, dimension(:), pointer :: disp_hgt, lut_hgt_in, lut_hgt_full, data_tmp

    if (initialized) call destroy_lut()
    
    call read_lut(lut_filename)
    if (error) return
        
    ! Interpolate the LUT to 12 months
    ! 
    allocate(lut_monthly(num_lats_in, num_levs_in, 12), stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed1', sub_name)) return
    do ind_month = 1, 12
      if (num_months_in == 4) then
        jday = fu_julian_date(fu_set_time_utc(2000, ind_month, 15, 12, 0, 0.0))
        if (jday > lut_jdays(4)) then
          i1 = 4
          i2 = 1
          div = real(365-lut_jdays(4) + lut_jdays(1))
        else if (jday > lut_jdays(3)) then
          i1 = 3
          i2 = 4
          div = lut_jdays(i2) - lut_jdays(i1)
        else if (jday > lut_jdays(2)) then
          i1 = 2
          i2 = 3
          div = lut_jdays(i2) - lut_jdays(i1)
        else
          i1 = 1
          i2 = 2
          div = lut_jdays(i2) - lut_jdays(i1)
        end if
        w2 = (jday-lut_jdays(i1)) / div
        w1 = 1 - w2
        lut_monthly(:,:,ind_month) = lut_data_in(:,:,i1)*w1 + lut_data_in(:,:,i2)*w2
      else
        lut_monthly(:,:,ind_month) = lut_data_in(:,:,ind_month)
      end if
    end do

    
    ! Pre-interpolate the LUT to dispersion levels that are covered. Higher are handled by
    ! another method.
    !
    call set_vertical((/(fu_set_level(constant_pressure, lut_press_in(ind_lev)), ind_lev=1, num_levs_in)/), &
                    & vert_lut)
    if (error) return
    call vert_to_metric(vert_lut, vert_tmp)
    if (error) return
    call vert_to_metric(vertical, vert_alt)
    if (error) return

    call report(vert_tmp, .true.)
    call report(vert_alt, .true.)

    disp_hgt => fu_work_array()
    lut_hgt_in => fu_work_array()
    lut_hgt_full => fu_work_array()
    lut_hgt_in => lut_hgt_full(2:)
    lut_hgt_in(1:num_levs_in) = (/(fu_level_height(fu_level(vert_tmp, ind_lev)), ind_lev=1, num_levs_in)/)
    do ind_lev = 1, fu_NbrOfLevels(vertical)
      lev_hgt = fu_level_height(fu_level(vert_alt, ind_lev))
      if (lev_hgt > lut_hgt_in(num_levs_in)) exit
      num_levs = ind_lev
      disp_hgt(ind_lev) = lev_hgt
    end do
    call set_missing(vert_tmp, ifNew=.false.)
    ! Store altitudes from the (metric or not) dispersion vertical:
    call set_missing(vert_alt, ifNew=.false.)
    allocate(lut_tmp(num_lats_in, num_levs, num_lut_months), &
           & lut_alt(num_levs), stat=stat)
    ! store the possibly non-metric dispersion level heights
    lut_alt = disp_hgt(1:num_levs)
    call msg('LUT covers up to', lut_alt(num_levs))

    if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
    lut_hgt_full(1) = 0.0
    data_tmp => fu_work_array()
    do ind_month = 1, 12
      do ind_lat =1, num_lats_in
        data_tmp(2:num_levs_in+1) = lut_monthly(ind_lat,:,ind_month)
        data_tmp(1) = data_tmp(2)
        call interp1d(lut_hgt_full(1:num_levs_in+1), data_tmp(1:num_levs_in+1), &
                    & disp_hgt(1:num_levs), lut_tmp(ind_lat, :, ind_month))
        if (error) return
      end do
    end do
    call free_work_array(disp_hgt)
    call free_work_array(data_tmp)
    call free_work_array(lut_hgt_full)
    deallocate(lut_monthly)
    
    ! The stratospheric parametrization is used only if the LUT ends below 19 km.
    ! 
    allow_strato = lut_alt(num_levs) < 19e3
    
    ! Interpolate to the latitude grid:
    ! 
    num_lats = 180 / lat_incr + 1
    allocate(lut_lat(num_lats), lut_data(num_lats, num_levs, num_lut_months), stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
    lut_lat = (/(-90 + lat_incr*ind_lat, ind_lat = 0, num_lats-1)/)

    do ind_month = 1, num_lut_months
      do ind_lev = 1, num_levs
        call interp1d(lut_lats_in(1:num_lats_in), lut_tmp(:,ind_lev,ind_month), &
                    & lut_lat, lut_data(:,ind_lev,ind_month))
        if (error) return
      end do
    end do
    deallocate(lut_tmp)
    ! Cleanup numerics
    if (any(lut_data < -1e8)) then
      call set_error('Negative values in the final LUT', sub_name)
      return
    end if
    lut_data = max(lut_data, 0.0)

    call normalize_lut()
    
    call msg('OH LUT init done')
    initialized = .true.

  end subroutine setup_oh_lut

  subroutine normalize_lut()
    implicit none
    integer :: ind_month, ind_lat, counter
    real :: rate_o3_o1d, rate_mean, sza_cos
    type(silja_time) :: now

    do ind_lat = 1, num_lats
      do ind_month = 1, 12
        now = fu_set_time_utc(2001, ind_month, 1, 0, 0, 0.0)
        rate_mean = 0.0
        counter = 0
        do while (fu_mon(now) == ind_month)
          call solarsetup(now)
          sza_cos = fu_solar_zenith_angle_cos(0.0, lut_lat(ind_lat), now)
          rate_o3_o1d = photo(6.073e-5, 1.743, 0.474, sza_cos)
          rate_mean = rate_mean + rate_o3_o1d
          counter = counter + 1
          now = now + one_hour
        end do
        rate_mean = rate_mean / counter
        if (rate_mean < 1e-12) then
          lut_data(ind_lat,:,ind_month) = 0.0
        else
          lut_data(ind_lat,:,ind_month) = lut_data(ind_lat,:,ind_month) / rate_mean
        end if
      end do
    end do
    

  end subroutine normalize_lut

  subroutine interp1d(x_in, val_in, x_out, val_out)
    ! Interpolate from arbitrary x_in, val_in to x_out, val_out. The axes may be arbitrary
    ! but must be increasing. x_in must cover x_out.
    implicit none
    real, dimension(:), intent(in) :: x_in, val_in, x_out
    real, dimension(:), intent(out) :: val_out
    character(len=*), parameter :: sub_name = 'interp1d'

    integer :: nx_in, nx_out, ind_out, ind_upper
    real :: weight_up
    
    nx_in = size(x_in)
    if (fu_fails(size(val_in) >= nx_in, 'val_in too small', sub_name)) return
    nx_out = size(x_out)
    if (fu_fails(size(val_out) >= nx_out, 'val_out too small', sub_name)) return
    if (fu_fails(x_out(nx_out) <= x_in(nx_in) .and. x_out(1) >= x_in(1), 'x_out not covered', sub_name)) return

    ind_upper = 2
    do ind_out = 1, nx_out
      do while (x_out(ind_out) > x_in(ind_upper) .and. ind_upper < nx_in)
        ind_upper = ind_upper + 1
      end do
      weight_up = (x_out(ind_out) - x_in(ind_upper-1)) / (x_in(ind_upper)-x_in(ind_upper-1))
      val_out(ind_out) = val_in(ind_upper-1)*(1-weight_up) + val_in(ind_upper)*weight_up
    end do
    
  end subroutine interp1d

  real function fu_oh_cnc_clim(lat, lon, alt, when) result(cnc)
    implicit none
    real, intent(in) :: lat, lon, alt
    type(silja_time), intent(in) :: when

    real :: cnc_lut_top, cnc_19km, slope

    if (alt <= lut_alt(num_levs)) then
      ! Table lookup
      cnc = fu_cnc_oh_lut(lat, lon, alt, when)
    else if (alt <= 19e3 .and. allow_strato) then
      ! Between the LUT top and the nominal level of the stratospheric fit.
      cnc_lut_top = fu_cnc_oh_lut(lat, lon, lut_alt(num_levs), when)
      cnc_19km = fu_cnc_oh_strato(lat, lon, when)
      slope = (cnc_19km - cnc_lut_top) / (19e3 - lut_alt(num_levs))
      cnc = cnc_lut_top + (alt - lut_alt(num_levs)) * slope
    else if (alt <= 25e3 .and. allow_strato) then
      ! Above stratospheric fit, but not too much.
      cnc = fu_cnc_oh_strato(lat, lon, when)
    else
      call msg('Request:', alt)
      call set_error('Cannot produce OH this high', 'fu_cnc_oh_clim')
    end if

  contains

    real function fu_cnc_oh_strato(lat, lon, when) result(cnc)
      ! Apply the formula of Hanisco et al., 2001. 
      ! 
      implicit none
      real, intent(in) :: lat, lon
      type(silja_time), intent(in) :: when
      
      real, parameter :: molec_to_mol = 1e6 / avogadro
      real, parameter :: j01 = 8e-5, j02 = 6.4e-5, &
           & z1 = 2.2, z2 = 0.67, &
           & c1 = 1.5e10*molec_to_mol, c2 = 1.75e10*molec_to_mol
      real :: sza, cos_chi1, cos_chi2, j1, j2
      
      sza = acos(fu_solar_zenith_angle_cos(lon, lat, when))
      cos_chi1 = cos(0.85*sza)
      if (cos_chi1 > 1e-12) then
        j1 = j01 * exp(-z1 * (1.0/cos_chi1 - 1))
      else 
        j1 = 0.0
      end if

      cos_chi2 = cos(0.86*sza)
      if (cos_chi2 > 1e-12) then
        j2 = j02 * exp(-z2 * (1.0/cos_chi2 - 1))
      else 
        j2 = 0.0
      end if

      cnc = c1*j1 + c2*j2
      
    end function fu_cnc_oh_strato

    real function fu_cnc_oh_lut(lat, lon, alt, when) result(cnc)
      implicit none
      real, intent(in) :: lat, lon, alt
      type(silja_time), intent(in) :: when

      integer :: jday, i1, i2, num_levs, ind_lat, ind_lev, ind_month
      real :: w1, w2, div, lut_val, rate_o3_o1d, sza_cos, dist, mindist
      real, parameter :: molec2mol = 1e6 / avogadro
      character(len=*), parameter :: sub_name = 'fu_cnc_oh_clim'

      sza_cos = fu_solar_zenith_angle_cos(lon, lat, when)

      ind_month = fu_mon(when)

      ind_lat = int((lat+90)/lat_incr + 1 + 0.5)

      mindist = 1e12
      num_levs = size(lut_data, 2)
      do ind_lev = 1, num_levs
        dist = abs(alt - lut_alt(ind_lev))
        if (dist < mindist) mindist = dist
        if (dist > mindist) exit
      end do
      ind_lev = ind_lev - 1
      if (fu_fails(ind_lev > 0 .and. ind_lev <= num_levs, 'Failed to find level', sub_name)) return

      lut_val = lut_data(ind_lat, ind_lev, ind_month)

      ! Rate from MCM (Sanders 2003), slightly different from the EMEP values in CB4.
      rate_o3_o1d = photo(6.073e-5, 1.743, 0.474, sza_cos)

      cnc = rate_o3_o1d * lut_val

    end function fu_cnc_oh_lut

  end function fu_oh_cnc_clim

  real function photo(l, m, n, cos_theta)
    implicit none
    real, intent(in) :: l, n, m, cos_theta
    if (cos_theta < 1e-5) then
      photo = 0.0
      return
    end if
    photo = l * cos_theta**m * exp(-n / cos_theta)
  end function photo

end module hydroxyl
