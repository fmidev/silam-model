program test_enkf
  use enkf
  use toolbox, only : random_normal
  use correlations
  use silam_levels
  use grids_geo
  use optimisation, only : report_time 
  implicit none
  
  run_log_funit = 17
  
  !call test_1d()
  call test_3d(400)

  call report_time()
  
contains
  
  subroutine test_1d()
    implicit none
    integer, parameter :: ens_size = 10000

    !real, dimension(mdl_size, ens_size) :: ens
    real, dimension(:,:), allocatable :: ens, ens_old
    type(silam_vertical) :: vert
    real :: obs_height, flev, weight_top, weight_bot
    integer :: ilev, ind_top, ind_bot, file_unit
    character(len=*), parameter :: filename_conf='test_enkf.nl', sub_name='test_1d'
    type(Tsilam_namelist), pointer :: nlptr
    type(t_spatial_correlation), dimension(:), pointer :: sp_corr
    type(t_correlation) :: corr
    real, dimension(:), allocatable :: corr_wrk
    integer :: ind_ens, iostat, mdl_size_uncorr
    integer :: mdl_size

    real, dimension(1, ens_size) :: ens_obs
    real, dimension(2,1) :: obs_loc = reshape((/0.0, 0.0/), (/2,1/))
    real, dimension(1) :: obs_data = 10
    real, dimension(1) :: obs_var = 0.5**2
    real, dimension(:,:), allocatable :: mdl_loc

    if (error) return

    file_unit = fu_next_free_unit()
    open(file_unit, file=filename_conf, iostat=iostat, action='read', status='old')
    if (fu_fails(iostat == 0, 'Failed to open '//trim(filename_conf), sub_name)) return
    nlptr => fu_read_namelist(file_unit, ifEnv=.true.)
    close(file_unit)
    
    call set_vertical(nlptr, vert)
    mdl_size = fu_NbrOfLevels(vert)
    allocate(ens(mdl_size, ens_size), ens_old(mdl_size, ens_size), corr_wrk(mdl_size), mdl_loc(2,mdl_size))
    mdl_loc = 0.0

    allocate(sp_corr(1))
    call set_spatial_corr_from_nl(nlptr, grid_missing, sp_corr(1))
    if (error) return
    call set_total_correlation(sp_corr, grid_missing, vert, 1, corr)
    if (error) return
    
    obs_height = 25.0
    flev = fu_project_level_crude(fu_set_level(constant_height, fval1=obs_height), vert)
    if (error) return
    
    ind_top = ceiling(flev)
    ind_bot = floor(flev)
    weight_top = (flev - ind_bot)
    weight_bot = 1 - weight_top

    call msg('indices', ind_top, ind_bot)
    call msg('weights', weight_top, weight_bot)
    
    do ind_ens = 1, ens_size
      call random_normal(ens(:, ind_ens))
      mdl_size_uncorr = fu_dimension_control(corr)
      call transf_fwd(corr, ens(1:mdl_size_uncorr, ind_ens), corr_wrk)
      if (error) return
      ens(:,ind_ens) = corr_wrk
      ens_obs(1,ind_ens) = ens(ind_bot,ind_ens)*weight_bot + ens(ind_top,ind_ens)*weight_top
    end do
    
    !call msg('First member:', ens(:,1))
    !call msg('First ens_obs', ens_obs(1,1))
    !call msg('ens_obs cross', ens_obs(:,1))
    
    ens_old = ens
    !call msg('stdev before', fu_ens_stdev(ens))
    !call msg('mean before', fu_ens_mean(ens))
    write(25, fmt=*) ens
    !call enkf_update(ens, ens_obs, obs_data, obs_var, obs_loc, mdl_loc, loc_dist_m = 100.0, &
    !               & loc_type=loc_none, enkf_flavor=enkf_enkf, diagn_out_dir='.', rfactor=1.0)
    call enkf_update_openmp(ens, ens_obs, obs_data, obs_var, obs_loc, mdl_loc, loc_dist_m = 100.0, &
                          & loc_type=loc_none, enkf_flavor=enkf_enkf, diagn_out_dir='.', rfactor=1.0)

    write(26, fmt=*) ens
    !call msg('stdev after', fu_ens_stdev(ens))
    !call msg('mean after', fu_ens_mean(ens))
    !call msg('First member after:', ens(:,1))
    !call msg('Increment on the mean:',  fu_ens_mean(ens_old) - fu_ens_mean(ens))
        
  end subroutine test_1d

  subroutine test_3d(num_obs)
    implicit none
    integer, parameter :: ens_size = 200
    !integer, parameter :: num_obs = 1
    integer, parameter :: num_species = 1

    integer, intent(in) :: num_obs

    !real, dimension(mdl_size, ens_size) :: ens
    real, dimension(:,:), allocatable :: ens, ens_old
    type(silam_vertical) :: vert
    real :: obs_height, flev, weight_top, weight_bot, lon, lat
    integer :: ilev, ind_top, ind_bot, file_unit
    character(len=*), parameter :: filename_conf='test_enkf.nl', sub_name='test_1d'
    type(Tsilam_namelist), pointer :: nlptr
    type(t_spatial_correlation), dimension(:), pointer :: sp_corr
    type(t_correlation) :: corr
    real, dimension(:), allocatable :: corr_wrk
    integer :: ind_ens, iostat, mdl_size_uncorr
    integer :: mdl_size, nx, ny, ind_obs, iz, isp, i1d, ix, iy, nz

    real, dimension(num_obs,ens_size) :: ens_obs
    real, dimension(2,num_obs) :: obs_loc
    real, dimension(num_obs) :: obs_data
    real, dimension(num_obs) :: obs_var
    real, dimension(:,:), allocatable :: mdl_loc
    real, dimension(2) :: rand

    real, dimension(:), pointer :: mean, stdev, row

    type(silja_grid) :: grd

    obs_data = 10.0
    obs_var = 2.0**2

    file_unit = fu_next_free_unit()
    open(file_unit, file=filename_conf, iostat=iostat, action='read', status='old')
    if (fu_fails(iostat == 0, 'Failed to open '//trim(filename_conf), sub_name)) return
    nlptr => fu_read_namelist(file_unit, ifEnv=.true.)
    close(file_unit)

    call set_vertical(nlptr, vert)
    if (error) return
    grd = fu_set_grid(nlptr)
    if (error) return
    call grid_dimensions(grd, nx, ny)
    mdl_size = fu_NbrOfLevels(vert)*nx*ny
    print *, 'mdl_size:', mdl_size
    allocate(ens(mdl_size, ens_size), corr_wrk(mdl_size), mdl_loc(2,mdl_size))
    
    call get_grid_loc(grd, vert, 1, mdl_loc)
    
    allocate(sp_corr(1))
    call set_spatial_corr_from_nl(nlptr, grd, sp_corr(1))
    if (error) return
    call set_total_correlation(sp_corr, grd, vert, 1, corr)
    if (error) return

    do ind_ens = 1, ens_size
      call random_normal(ens(:, ind_ens))
      mdl_size_uncorr = fu_dimension_control(corr)
      call transf_fwd(corr, ens(1:mdl_size_uncorr, ind_ens), corr_wrk)
      if (error) return
      ens(:,ind_ens) = corr_wrk
    end do

    iz = 1
    isp = 1
    nz = fu_NbrOfLevels(vert)

    do ind_obs = 1, num_obs
      call random_number(rand)
      ix = max(nint(rand(1)*nx), 1)
      iy = max(nint(rand(2)*ny), 1)
      lon = fu_lon_geographical_from_grid(real(ix), real(iy), grd)
      lat = fu_lat_geographical_from_grid(real(ix), real(iy), grd)
      obs_loc(1, ind_obs) = lon
      obs_loc(2, ind_obs) = lat
      print *, 'obs:', lon, lat, ix, iy
      i1d = (iy-1)*(nx*nz*num_species) + (ix-1)*(nz*num_species) + (iz-1)*num_species + isp
      ens_obs(ind_obs,:) = ens(i1d,:)
    end do
    
!!$    ind_obs = 0
!!$    do iy = 1, ny, ny/int(sqrt(real(num_obs)))
!!$      do ix = 1, nx, nx/int(sqrt(real(num_obs)))
!!$        ind_obs = ind_obs + 1
!!$        if (ind_obs > num_obs) exit
!!$        lon = fu_lon_geographical_from_grid(real(ix), real(iy), grd)
!!$        lat = fu_lat_geographical_from_grid(real(ix), real(iy), grd)
!!$        obs_loc(1, ind_obs) = lon
!!$        obs_loc(2, ind_obs) = lat
!!$        print *, 'obs:', lon, lat, ix, iy
!!$        i1d = (iy-1)*(nx*nz*num_species) + (ix-1)*(nz*num_species) + (iz-1)*num_species + isp
!!$        ens_obs(:,ind_obs) = ens(i1d,:)
!!$      end do
!!$    end do

    call dump_arr(obs_loc(1,:), 'obs_lon.dat')
    call dump_arr(obs_loc(1,:), 'obs_lat.dat')

    call dump_arr(fu_ens_stdev(ens), 'fc_stdev.dat')
    call dump_arr(fu_ens_mean(ens), 'fc_mean.dat')

    mean => fu_work_array(mdl_size)
    stdev => fu_work_array(mdl_size)
    row => fu_work_array(mdl_size)
    ix = 2
    iy = 1
    i1d = (iy-1)*(nx*nz*num_species) + (ix-1)*(nz*num_species) + (iz-1)*num_species + isp
    call get_row(ens, i1d, row(1:mdl_size), stdev(1:mdl_size), mean(1:mdl_size), 'corr')
    call dump_arr(row(1:mdl_size), 'corr_fc.dat')
    
    !call enkf_update(ens, ens_obs, obs_data, obs_var, obs_loc, mdl_loc, loc_dist_m = 500.0e3, &
    !               & loc_type=loc_step, enkf_flavor=enkf_enkf, diagn_out_dir='.', rfactor=1.0)

    call enkf_update_openmp(ens, ens_obs, obs_data, obs_var, obs_loc, mdl_loc, loc_dist_m = 500.0e3, &
                   & loc_type=loc_step, enkf_flavor=enkf_enkf, diagn_out_dir='.', rfactor=1.0)


    call dump_arr(fu_ens_mean(ens), 'an_mean.dat')
    call dump_arr(fu_ens_stdev(ens), 'an_stdev.dat')
    
    call get_row(ens, i1d, row(1:mdl_size), stdev(1:mdl_size), mean(1:mdl_size), 'corr')
    call dump_arr(row(1:mdl_size), 'corr_an.dat')
    
  end subroutine test_3d


  subroutine get_grid_loc(grid, vertical, num_species, loc)
    ! Set the lon/lat localisation for a control vector in mass map order.
    implicit none
    type(silja_grid), intent(in) :: grid
    type(silam_vertical), intent(in) :: vertical
    integer, intent(in) :: num_species
    real, dimension(:,:), intent(out) :: loc

    integer :: nx, ny, nz
    integer :: ix, iy, iz, isp, i1d

    call grid_dimensions(grid, nx, ny)
    nz = fu_NbrOfLevels(vertical)
    
    if (fu_fails(nz*nx*ny*num_species <= size(loc,2), 'loc too small', 'get_grid_loc')) return
    
    do iy = 1, ny
      do ix = 1, nx
        do iz = 1, nz
          do isp = 1, num_species
            i1d = (iy-1)*(nx*nz*num_species) + (ix-1)*(nz*num_species) + (iz-1)*num_species + isp
            loc(1,i1d) = fu_lon_geographical_from_grid(real(ix), real(iy), grid)
            loc(2,i1d) = fu_lat_geographical_from_grid(real(ix), real(iy), grid)
          end do
        end do
      end do
    end do
    
    
  end subroutine get_grid_loc

  function fu_ens_mean(ens) result(mean)
    implicit none
    real, dimension(:,:) :: ens
    real, dimension(size(ens, 1)) :: mean

    mean = sum(ens, dim=2) / size(ens, 2)    
    
  end function fu_ens_mean

  function fu_ens_stdev(ens) result(stdev)
    implicit none
    real, dimension(:,:) :: ens
    real, dimension(size(ens, 1)) :: stdev

    real, dimension(size(ens, 1)) :: mean
    integer :: ii

    mean = sum(ens, dim=2) / size(ens, 2)
    stdev = 0.0
    do ii = 1, size(ens, 2)
      stdev(:) = stdev + (mean - ens(:,ii))**2
    end do
    stdev = sqrt(stdev / size(ens, 2))
    

  end function fu_ens_stdev
  
  
  subroutine get_vert(nlevs, dz, vert)
    implicit none
    integer, intent(in) :: nlevs
    real, intent(in) :: dz
    type(silam_vertical), intent(out) :: vert

    type(silja_level), dimension(nlevs) :: levels
    integer :: ilev

    levels = (/(fu_set_level(constant_height, fval1=dz*ilev), ilev=1, nlevs)/)
    
    call set_vertical(levels, vert)
    
  end subroutine get_vert

  subroutine dump_arr(arr, filename)
    implicit none
    real, dimension(:), intent(in) :: arr
    character(len=*), intent(in) :: filename
    
    integer :: file_unit, iostat

    file_unit = fu_next_free_unit()
    open(file_unit, file=filename, access='stream', form='unformatted', iostat=iostat)
    if (fu_fails(iostat == 0, 'Failed to open: ' // trim(filename), 'dump_arr')) return
    
    write(file_unit) arr
    close(file_unit)

  end subroutine dump_arr


  
  
  
end program test_enkf
