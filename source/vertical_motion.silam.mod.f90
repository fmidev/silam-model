module vertical_motion
  use physiographies
  use omp_lib

  private
  
  public diagnose_eta_dot
  public eta_from_omega

  private solve_poisson
  ! The vertical_motion_precision = vm_p parameter. This is the
  ! working (single/double) precision for the module. It needs to be
  ! synchronized with the fishpack routines.
  ! 
  ! While single precision is surely ok for integrating the continuity
  ! equation, the poisson solver seems to work better in double
  ! precision. On the other hand, the difference in the flux
  ! correction (= the gradient of the solution) does not seem to
  ! change much.
  !
  ! The difference in performance seems to be insignificant.
#ifdef DOUBLE_PRECISION
  integer, parameter, public :: vm_p = r8k
#else
  integer, parameter, public :: vm_p = r4k
#endif

contains
  
  subroutine diagnose_eta_dot(u3d, v3d, ps, spt, layers, u_grid, v_grid, p_grid, &
                            & use_pressure_fixer, go_top_down, eta_dot)
    !
    ! Diagnose the vertical motion in eta-coordinates. Has separate
    ! schemes for regular and staggered grids. Pressure fixer can be
    ! requested if a staggered grid is used. A non-staggered version
    ! has been around, but did not perform too well.
    implicit none
    real(vm_p), dimension(:,:), intent(inout) :: u3d, v3d ! horizontal winds
    real(vm_p), dimension(:), intent(in) :: ps, spt ! surface pressure + tendency
    real(vm_p), dimension(:,:), intent(out) :: eta_dot
    type(silam_vertical), intent(in) :: layers ! The vertical to compute in
    type(silja_grid), intent(in) :: p_grid, u_grid, v_grid ! The grids - either the same or Ara C.
    logical, intent(in) :: use_pressure_fixer, go_top_down

    real, dimension(:), pointer :: hybrid_a, hybrid_b
    real(vm_p), dimension(:), allocatable :: massflux_u, massflux_v, div_x, div_y, residual, fx, fy, &
         & dps_dx, dps_dy, div_sum, ps_sum
    integer :: nx, ny, fs, nlevs, ilev, i2d, sgn_x_shift, sgn_y_shift
    type(silja_level) :: lev_bot, lev_top
    real(vm_p) :: da, db, dp
    logical :: staggered

    call init()
    call start_count('pressure_fixer')
    if (use_pressure_fixer .and. .not. staggered) then 
      call msg_warning('Pressure fixer only available for staggered grids')
    end if
    if (use_pressure_fixer .and. staggered) call pressure_fixer_stg()
    call stop_count('pressure_fixer')
    call start_count('integrate_to_eta')

    if (go_top_down) then
      if (staggered) then
        call get_eta_dot_down_stg()
      else
        call get_eta_dot_down()
      end if
    else
      if (staggered) then
        call get_eta_dot_up_stg()
      else
        call get_eta_dot_up()
      end if
    end if

    !call get_eta_dot_down()
    call stop_count('integrate_to_eta')
    
    call cleanup()

  contains

    subroutine init()
      implicit none
      real :: shift_lon, shift_lat

      call grid_dimensions(p_grid, nx, ny)
      fs = nx*ny
      nlevs = fu_nbrOfLevels(layers)

      sgn_y_shift = 0
      sgn_x_shift = 0

      call grid_shift_indices(u_grid, p_grid, shift_lon, shift_lat)
      if (.not. (abs(shift_lon) .eps. 0.0)) then
        if (.not. (abs(shift_lat) .eps. 0.0) .or. .not. (abs(shift_lon) .eps. 0.5)) then
          call set_error('Bad staggering', 'diganose_eta_dot')
          return
        else
          staggered = .true.
          sgn_x_shift = sign(1.0, shift_lon)
        end if
      else
        staggered = .false.
      end if

      call grid_shift_indices(v_grid, p_grid, shift_lon, shift_lat)
      if (.not. (abs(shift_lat) .eps. 0.0)) then
        if (.not. (abs(shift_lon) .eps. 0.0) .or. .not. (abs(shift_lat) .eps. 0.5) .or. .not. staggered) then
          call set_error('Bad staggering', 'diganose_eta_dot')
          return
        end if
        sgn_y_shift = sign(1.0, shift_lat)
      end if
                
      allocate(massflux_u(fs), massflux_v(fs), &
             & hybrid_a(nlevs+1), hybrid_b(nlevs+1), &
             & div_x(fs), div_y(fs), &
             & residual(fs), &
             & fx(fs), fy(fs))
      allocate(dps_dx(fs), dps_dy(fs), div_sum(fs), ps_sum(fs))
      
      do ilev = 1, nlevs
        lev_bot = fu_lower_boundary_of_layer(fu_level(layers, ilev))
        hybrid_a(ilev) = fu_hybrid_level_coeff_a(lev_bot)
        hybrid_b(ilev) = fu_hybrid_level_coeff_b(lev_bot)
      end do
      lev_top = fu_upper_boundary_of_layer(fu_level(layers, nlevs))
      hybrid_a(nlevs+1) = fu_hybrid_level_coeff_a(lev_top)
      hybrid_b(nlevs+1) = fu_hybrid_level_coeff_b(lev_top)
      
    end subroutine init

    subroutine cleanup()
      implicit none
      integer :: stat
      deallocate(massflux_u, massflux_v, hybrid_a, hybrid_b, &
               & div_x, div_y, residual, fx, fy, &
               & dps_dx, dps_dy, div_sum, ps_sum, stat=stat)
      if (stat /= 0) then
        call set_error('Deallocate failed', 'cleanup')
        return
      end if
    end subroutine cleanup

    subroutine pressure_fixer_stg()
      ! Computes a 2d correction field to the horizontal winds using
      ! the "Poisson-equation method". This version assumes staggerd
      ! grids. Poisson equation is solved using the fishpack
      ! subroutine hwssp. 
      implicit none
      
      integer :: iter, niter = 10, ix, iy, ixu, iyv, ind_u, ind_v
      real(vm_p) :: dp_u, dp_v, ps_u, ps_v
      logical :: iterate
      character(len=fnlen) :: name
      real(vm_p), dimension(:), allocatable :: div_fx, div_fy
      
      allocate(div_fx(fs), div_fy(fs))

      fx = 0.0
      fy = 0.0
      call msg('Pressure fixer (staggered)')
      call msg('Mean spt from met:', sum(spt(1:fs))/fs)

      ! Compute the column-integrated mass flux, store temporarily into div_fx and div_fy.
      !
      do ilev = 1, nlevs
        da = hybrid_a(ilev) - hybrid_a(ilev+1)
        db = hybrid_b(ilev) - hybrid_b(ilev+1)
        do iy = 1, ny
          do ix = 1, nx
            i2d = (iy-1)*nx + ix
            ixu = ix + sgn_x_shift
            if (ixu > nx) ixu = nx
            if (ixu < 1) ixu = 1
            iyv = iy + sgn_y_shift
            if (iyv > ny) iyv = ny
            if (iyv < 1) iyv = 1

            ind_u = (iy-1)*nx + ixu
            ind_v = (iyv-1)*nx + ix
            ps_u = 0.5 * (ps(i2d) + ps(ind_u))
            ps_v = 0.5 * (ps(i2d) + ps(ind_v))
            dp_u = da + db * ps_u
            dp_v = da + db * ps_v
            div_fx(i2d) = div_fx(i2d) + dp_u*u3d(i2d,ilev)
            div_fy(i2d) = div_fy(i2d) + dp_v*v3d(i2d,ilev)
          end do
        end do
      end do
      
      ! Compute the divergence of the original mass flux:
      call ddx_of_field_dp(div_fx, u_grid, p_grid, div_x)
      call div_y_of_field_dp(div_fy, v_grid, p_grid, div_y)

      iter = 1
      iterate = .true.
      
      ! Next: fx, fy are the components of the correction. The idea is to iterate updating
      ! fx and fy until div_x + div_y + div(fx,fy) = spt (spt = surface pressure
      ! tendency). Note that this part only deals with the column-integrated fluxes, the
      ! attribution to levels is done below.
      !
      do while(iterate)
        iterate = iter <= niter
        call msg('Iteration', iter)

        call start_count('pres_fix_div')
        call ddx_of_field_dp(fx, u_grid, p_grid, div_fx)
        call div_y_of_field_dp(fy, v_grid, p_grid, div_fy)
        call stop_count('pres_fix_div')

        do i2d = 1, fs
          residual(i2d) = -(div_fx(i2d) + div_fy(i2d) + div_x(i2d) + div_y(i2d))
        end do

        call msg('Mean spt from div:', sum(residual) / fs)
        residual = residual - spt(1:fs)
        call msg('||r||^2 = ', sum(residual**2))
        if (.not. iterate) exit
        call start_count('poisson')
        call solve_poisson(residual, p_grid)
        ! use temporarily div_x, div_y, the solution is in residual.
        call ddx_of_field_dp(residual, p_grid, u_grid, div_fx)
        call ddy_of_field_dp(residual, p_grid, v_grid, div_fy)
        fx(1:fs) = fx + div_fx
        fy(1:fs) = fy + div_fy
        call stop_count('poisson')
        iter = iter + 1
      end do
      
      call flux_to_wind(fx, fy, u3d, v3d, ps, 'even')

    end subroutine pressure_fixer_stg
    
    subroutine flux_to_wind(flux_x, flux_y, u3d, v3d, ps, method)
      ! Translate the correction for the column-integrated mass flux
      ! into horizontal wind increments. There is a method for doing
      ! this to minimize the fractional difference (method =
      ! 'fractional'), but this doesn't seem to work too
      ! well. Otherwise, just weights the increments with the
      ! thickness of layer (in Pa).
      !
      implicit none
      real(vm_p), dimension(:), intent(in) :: flux_x, flux_y, ps
      real(vm_p), dimension(:,:), intent(inout) :: u3d, v3d
      !real(vm_p), dimension(:) :: work
      character(len=*), intent(in) :: method

      integer :: ix, iy, ind_u, ind_v, ixu, iyv
      real(vm_p), dimension(:), allocatable :: sum_u, sum_v
      real(vm_p) :: ps_mid, dp_u, dp_v
      

      if (method == 'even') then
        do ilev = 1, nlevs
          do iy = 1, ny
            do ix = 1, nx
              i2d = (iy-1)*nx + ix
              ixu = ix + sgn_x_shift
              if (ixu > nx) ixu = nx
              if (ixu < 1) ixu = 1
              iyv = iy + sgn_y_shift
              if (iyv > ny) iyv = ny
              if (iyv < 1) iyv = 1
              
              ind_u = (iy-1)*nx + ixu
              ind_v = (iyv-1)*nx + ix
              ps_mid = 0.5 * (ps(i2d) + ps(ind_u))
              dp_u = da + db * ps_mid
              ps_mid = 0.5 * (ps(i2d) + ps(ind_v))
              dp_v = da + db * ps_mid
              u3d(i2d,ilev) = u3d(i2d,ilev) + fx(i2d) / ps_mid
              v3d(i2d,ilev) = v3d(i2d,ilev) + fy(i2d) / ps_mid
            end do
          end do
        end do ! level

      else if (method == 'fractional') then
        ! First, the normalization.
        allocate(sum_u(fs), sum_v(fs))
        sum_u = 0.0
        sum_v = 0.0
        do ilev = 1, nlevs
          do iy = 1, ny
            do ix = 1, nx
              i2d = (iy-1)*nx + ix
              ixu = ix + sgn_x_shift
              if (ixu > nx) ixu = nx
              if (ixu < 1) ixu = 1
              iyv = iy + sgn_y_shift
              if (iyv > ny) iyv = ny
              if (iyv < 1) iyv = 1
              
              ind_u = (iy-1)*nx + ixu
              ind_v = (iyv-1)*nx + ix
              ps_mid = 0.5 * (ps(i2d) + ps(ind_u))
              dp_u = da + db * ps_mid
              ps_mid = 0.5 * (ps(i2d) + ps(ind_v))
              dp_v = da + db * ps_mid
              sum_u(i2d) = sum_u(i2d) + dp_u*dp_u*u3d(i2d,ilev)*u3d(i2d,ilev)
              sum_v(i2d) = sum_v(i2d) + dp_v*dp_v*v3d(i2d,ilev)*v3d(i2d,ilev)

            end do
          end do
        end do ! level

        do ilev = 1, nlevs
          do iy = 1, ny
            do ix = 1, nx
              i2d = (iy-1)*nx + ix
              ixu = ix + sgn_x_shift
              if (ixu > nx) ixu = nx
              if (ixu < 1) ixu = 1
              iyv = iy + sgn_y_shift
              if (iyv > ny) iyv = ny
              if (iyv < 1) iyv = 1
              
              ind_u = (iy-1)*nx + ixu
              ind_v = (iyv-1)*nx + ix
              ps_mid = 0.5 * (ps(i2d) + ps(ind_u))
              dp_u = da + db * ps_mid
              ps_mid = 0.5 * (ps(i2d) + ps(ind_v))
              dp_v = da + db * ps_mid
              u3d(i2d,ilev) = u3d(i2d,ilev) + fx(i2d) * u3d(i2d,ilev)*u3d(i2d,ilev)*dp_u / sum_u(i2d)
              v3d(i2d,ilev) = v3d(i2d,ilev) + fy(i2d) * v3d(i2d,ilev)*v3d(i2d,ilev)*dp_v / sum_v(i2d)
            end do
          end do
        end do ! level
        
        deallocate(sum_u, sum_v)
      
      else
        call set_error('Bad flux_to_wind method:' // trim(method), 'flux_to_wind')
        return
      end if
    
    end subroutine flux_to_wind

    subroutine get_eta_dot_up_stg()
      ! Cmopute eta-dot by integrating the continuity equation from
      ! ground to up. Assumes staggerd grids.
      implicit none
      real(vm_p) :: dp_u, dp_v
      integer :: ix, iy, ii, ind_u, ind_v, ixu, iyv
      
      call msg('...staggered, up')

      massflux_u = 0.0
      massflux_v = 0.0
      
      eta_dot(:,1) = 0.0

      do ilev = 1, nlevs
        da = hybrid_a(ilev) - hybrid_a(ilev+1)
        db = hybrid_b(ilev) - hybrid_b(ilev+1)
          
        do iy = 1, ny
          do ix = 1, nx
            ii = (iy-1)*nx + ix
            
            ixu = ix + sgn_x_shift
            if (ixu > nx) ixu = nx
            if (ixu < 1) ixu = 1
            iyv = iy + sgn_y_shift
            if (iyv > ny) iyv = ny
            if (iyv < 1) iyv = 1

            ind_u = (iy-1)*nx + ixu
            ind_v = (iyv-1)*nx + ix
            
            dp_u = da + db * 0.5 * (ps(ii) + ps(ind_u)) 
            dp_v = da + db * 0.5 * (ps(ii) + ps(ind_v))

            massflux_u(ii) = massflux_u(ii) + u3d(ii,ilev)*dp_u
            massflux_v(ii) = massflux_v(ii) + v3d(ii,ilev)*dp_v
          end do
        end do

        call ddx_of_field_dp(massflux_u, u_grid, p_grid, div_x)
        call div_y_of_field_dp(massflux_v, v_grid, p_grid, div_y)
        eta_dot(:,ilev+1) = (div_x + div_y) +  (1.0 - hybrid_b(ilev))*spt

      end do
            
    end subroutine get_eta_dot_up_stg

    subroutine get_eta_dot_down_stg()
      ! Compute eta-dot by integrating the continuiry equation
      ! starting from the top of atmosphere (p = 0). Assumes staggered
      ! grids.
      implicit none
      real(vm_p) :: dp_u, dp_v
      integer :: ix, iy, ii, ind_u, ind_v, ixu, iyv
      
      call msg('...staggered, down')

      massflux_u = 0.0
      massflux_v = 0.0
      
      eta_dot(:,nlevs+1) = 0.0

      do ilev = nlevs, 1, -1
        da = hybrid_a(ilev) - hybrid_a(ilev+1)
        db = hybrid_b(ilev) - hybrid_b(ilev+1)
          
        do iy = 1, ny
          do ix = 1, nx
            ii = (iy-1)*nx + ix
            
            ixu = ix + sgn_x_shift
            if (ixu > nx) ixu = nx
            if (ixu < 1) ixu = 1
            iyv = iy + sgn_y_shift
            if (iyv > ny) iyv = ny
            if (iyv < 1) iyv = 1

            ind_u = (iy-1)*nx + ixu
            ind_v = (iyv-1)*nx + ix
            dp_u = da + db * 0.5 * (ps(ii) + ps(ind_u)) 
            dp_v = da + db * 0.5 * (ps(ii) + ps(ind_v))

            massflux_u(ii) = massflux_u(ii) + u3d(ii,ilev)*dp_u
            massflux_v(ii) = massflux_v(ii) + v3d(ii,ilev)*dp_v
          end do
        end do

        call ddx_of_field_dp(massflux_u, u_grid, p_grid, div_x)
        call div_y_of_field_dp(massflux_v, v_grid, p_grid, div_y)
        eta_dot(:,ilev) = -(div_x + div_y) - (hybrid_b(ilev))*spt(1:fs)
        
      end do
            
    end subroutine get_eta_dot_down_stg

    subroutine get_eta_dot_up()
      ! Cmopute eta-dot by integrating the continuity equation from
      ! ground to up.
      implicit none
      
      integer :: i1d
      real :: dp

      massflux_u = 0.0
      massflux_v = 0.0

      eta_dot(1:fs,1) = 0.0

      do ilev = 1, nlevs
        da = hybrid_a(ilev) - hybrid_a(ilev+1)
        db = hybrid_b(ilev) - hybrid_b(ilev+1)
        do i1d = 1, fs
          dp = da + db*ps(i1d)
          massflux_u(i1d) = massflux_u(i1d) + dp*u3d(i1d, ilev)
          massflux_v(i1d) = massflux_v(i1d) + dp*v3d(i1d, ilev)
        end do
        call ddx_of_field_dp(massflux_u, u_grid, p_grid, div_x)
        call div_y_of_field_dp(massflux_v, v_grid, p_grid, div_y)
        do i1d = 1, fs
          dp = da + db*ps(i1d)
          eta_dot(i1d,ilev+1) = ((div_x(i1d) + div_y(i1d)) + (1.0-hybrid_b(ilev))*spt(i1d))!/dp
        end do
        
      end do

    end subroutine get_eta_dot_up

    subroutine get_eta_dot_down()
      ! Compute eta-dot by integrating the continuiry equation
      ! starting from the top of atmosphere (p = 0).
      implicit none
      massflux_u = 0.0
      massflux_v = 0.0
      eta_dot(:,nlevs+1) = 0.0
      do ilev = nlevs, 1, -1
        da = hybrid_a(ilev) - hybrid_a(ilev+1)
        db = hybrid_b(ilev) - hybrid_b(ilev+1)
        do i2d = 1, fs
          dp = da + db*ps(i2d)
          massflux_u(i2d) = massflux_u(i2d) + dp*u3d(i2d, ilev)
          massflux_v(i2d) = massflux_v(i2d) + dp*v3d(i2d, ilev)
        end do
        call ddx_of_field_dp(massflux_u, u_grid, p_grid, div_x)
        call div_y_of_field_dp(massflux_v, v_grid, p_grid, div_y)
        eta_dot(:,ilev) = -(div_x + div_y) - hybrid_b(ilev)*spt(1:fs)
        
      end do
      
    end subroutine get_eta_dot_down

  end subroutine diagnose_eta_dot

  !************************************************************************************

  subroutine solve_poisson(phi, grid)
    ! Solve the poisson equation /\y = phi, where phi gives initially
    ! the forcing and on return contains the solution. 
    implicit none
    real(vm_p), dimension(:), intent(inout) :: phi
    type(silja_grid), intent(in) :: grid

    real(vm_p), dimension(:,:), allocatable:: phi2d
    real(vm_p), dimension(:), allocatable :: d_at_ts, d_at_tf, d_at_ps, d_at_pf
    real(vm_p), dimension(:), pointer :: work
    integer :: nx, ny, leading_dim_f, ierror, n, m, ix, iy, boundary_flag_x, boundary_flag_y, &
         & ext, im, in
    real :: x0, y0, dx, dy, southpole_lon, southpole_lat
    real(vm_p) :: x1, y1, dummy, lambda, perturbation, &
         & ts, tf, ps, pf, R2

    logical :: corner_in_geo_lonlat

    ext = 0
    
    call set_error('Poisson not working', 'solve_poisson')
    return

    if (fu_gridtype(grid) /= lonlat) then
      call set_error('The pressure fixer is only available for lon/lat grids.', &
           &  'solve_poisson')
      return
    end if
    if (fu_ifLonGlobal(grid) .or. fu_ifLatGlobal(grid)) then
      call set_error('Global grids not yet working', 'solve_poisson')
      return
    end if

    call lonlat_grid_parameters(grid, x0, y0, corner_in_geo_lonlat, nx, ny, &
         & southpole_lon, southpole_lat, dx, dy)

    ts = (y0 - dy*(ext+1)) * degrees_to_radians + pi/2
    tf = (y0 + (ny+ext+1)*dy) * degrees_to_radians + pi/2
    ps = (x0 - dx*(ext+1)) * degrees_to_radians + pi
    pf = (x0 + (nx+ext+1)*dx) * degrees_to_radians + pi

    boundary_flag_x = 3 ! neumann
    boundary_flag_y = 3 ! neumann
    m = ny - 1
    n = nx - 1

    boundary_flag_x = 1 ! neumann
    boundary_flag_y = 1 ! neumann
    m = ny + 2*ext + 1
    n = nx + 2*ext + 1

    allocate(phi2d(m+1, n+1), d_at_ts(n+1), d_at_tf(n+1), d_at_ps(m+1), d_at_pf(m+1), work(n*m))
    phi2d = 0.0
    d_at_ts = 0.0
    d_at_tf = 0.0
    d_at_ps = 0.0
    d_at_pf = 0.0
    R2 = earth_radius**2
    do ix = 1, nx 
      do iy = 1, ny
        phi2d(iy+ext+1,ix+ext+1) = phi((iy-1)*nx + ix) * R2
      end do
    end do

    leading_dim_f = m+1
    lambda = 0.0
    !call msg('...Solving the poisson equation')
    !call hwsssp(ts, tf, m, boundary_flag_y, d_at_ts, d_at_tf, &
    !     & ps, pf, n, boundary_flag_x, d_at_ps, d_at_pf, &
    !     & 0.0, &
    !     & phi2d, leading_dim_f, &
    !     & perturbation, &
    !     & ierror, work)
    if (perturbation > 0.0) call msg('...perturbation = ', perturbation)

    if (ierror /= 0) then
      call msg('ierror = ', ierror)
      call set_error('hwsssp had an error', 'solve_poisson')
      return
    end if

    do ix = 1, nx
      do iy = 1, ny
        phi((iy-1)*nx + ix) = phi2d(iy+ext+1,ix+ext+1)
      end do
    end do

  end subroutine solve_poisson

  !************************************************************************************

  subroutine eta_from_omega(u3d, v3d, ps, spt, omega, layers, u_grid, v_grid, p_grid, eta_dot)
    !
    ! Compute eta-dot from omega. Eta dot is required on half-levels while the omega and
    ! all other fields are on full levels. Consequently, eta dot is now first computed on
    ! full levels, then interpolated into half levels.
    ! 
    implicit none
    real(vm_p), dimension(:,:), intent(inout) :: u3d, v3d ! horizontal winds
    real(vm_p), dimension(:), intent(in) :: ps, spt ! surface pressure + tendency
    real(vm_p), dimension(:,:), intent(in) :: omega
    real(vm_p), dimension(:,:), intent(out) :: eta_dot
    type(silam_vertical), intent(in) :: layers ! The vertical to compute in
    type(silja_grid), intent(in) :: p_grid, u_grid, v_grid ! The grids - either the same or Ara C.
    
    real, dimension(:), pointer :: hybrid_b_full
    real :: press_term
    integer :: nlevs, ilev, fs, i1d
    type(silja_level) :: level_full
    real(vm_p), dimension(:), allocatable :: dps_dy, dps_dx

    if (.not. (p_grid == u_grid .and. u_grid == v_grid)) then
      call set_error('Staggering not supported', 'eta_from_omega')
      return
    end if

    call init()
    if (error) return

    ! eta dot * dp/d_eta = omega - dp/dt - grad(p) dot V_horiz 
    ! dp/dt = B*dps/dt
    ! grad(p) = B*grad(ps)

    ! note eta_dot is on half levels (total nlevs + 1), everything else on full levels. 
    eta_dot(:,1) = 0.0
    do ilev = 1, nlevs
      do i1d = 1, fs
        press_term = spt(i1d) + dps_dx(i1d)*u3d(i1d,ilev) + dps_dy(i1d)*v3d(i1d,ilev)
        eta_dot(i1d,ilev+1) = omega(i1d,ilev) - hybrid_b_full(ilev) * press_term
      end do
      if (ilev > 1) then
        ! interpolate to half level
        do i1d = 1, fs
          eta_dot(i1d,ilev) = 0.5 * (eta_dot(i1d,ilev) + eta_dot(i1d,ilev+1))
        end do
      end if
    end do
    eta_dot(:,nlevs+1) = 0.0
    call cleanup()

  contains
    
    subroutine init()
      implicit none
      integer :: nx, ny, stat

      nlevs = fu_NbrOfLevels(layers)
      hybrid_b_full => fu_work_array()
      call grid_dimensions(p_grid, nx, ny)
      fs = nx*ny

      allocate(dps_dx(fs), dps_dy(fs), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', 'eta_from_ometa')) return

      !dps_dx => fu_work_array()
      !dps_dy => fu_work_array()

      do ilev = 1, nlevs
        level_full = fu_level(layers, ilev, .true.)
        hybrid_b_full(ilev) = fu_hybrid_level_coeff_b(level_full)
      end do

      call ddx_of_field_dp(ps(1:fs), p_grid, p_grid, dps_dx)
      if (error) return
      call ddy_of_field_dp(ps(1:fs), p_grid, p_grid, dps_dy)
    end subroutine init

    subroutine cleanup()
      call free_work_array(hybrid_b_full)
      deallocate(dps_dx, dps_dy)
    end subroutine cleanup
    
  end subroutine eta_from_omega

  !************************************************************************************

  ! The following contains copies of the ddx/ddy methods in the grids
  ! module. If single precision seems sufficient, these might be just
  ! removed.

  SUBROUTINE ddx_of_field_dp(grid_data, grd, grd_out, ddx)
    !
    ! Calculates d/dx of field (which has grid grd)
    ! in every gridpoint of grid grd_out.
    !
    ! The two grids must either be exactly the same or match one
    ! of the Arakawa grid systems (half a gridpoint difference in either
    ! direction).
    !
    ! Positive derivative to Silja's positive X-direction eastwards.
    ! 
    ! Unit is [unit-of-grid-data]/m.
    !
    ! Handles any +-1/2 cell shift.
    !
    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL(vm_p), DIMENSION(:) :: grid_data ! the data itself
    TYPE(silja_grid), INTENT(in) :: grd ! the grid the data is in
    TYPE(silja_grid), INTENT(in) :: grd_out ! the grid in which' points
    ! we want to have the values of derivatives

    ! Imported parameters with intent(out):
    REAL(vm_p), DIMENSION(:), INTENT(out) :: ddx

    ! Local declarations:
    REAL(vm_p) :: lat, dist_x
    REAL(vm_p) :: lat1, dist_x1
    REAL(vm_p) :: lat2, dist_x2
    real :: shift_lon, shift_lat
    real(vm_p) :: fTmp
    INTEGER :: i, j, ind1, ind2, nx, ny, nxout, nyout, d, ixtr
    logical :: lon_global

    ! Check and initialize.
    
    IF (.NOT.fu_grids_match_closely(grd, grd_out)) THEN
      CALL report(grd)
      CALL report(grd_out)
      CALL set_error('cannot differentiate, no match','ddx_of_field')
      RETURN
    END IF
    
    lon_global = fu_ifLonGlobal(grd)

!    IF (grd%gridtype /= lonlat ) THEN
!      CALL set_error('sorry, can only handle lonlat grids','ddx_of_field')
!      RETURN
!    END IF

    IF (SIZE(grid_data) > SIZE(ddx)) THEN
      CALL set_error('ddx vector too small','ddx_of_field')
      RETURN
    END IF

    CALL grid_dimensions(grd, nx, ny)
    call grid_dimensions(grd_out, nxout, nyout)
    IF (error) RETURN
  
    if(nx /= nxout .or. ny /= nyout) then
      call set_error('Grids are of different size.', 'ddx_of_field')
      return
    end if
    
    call grid_shift_indices(grd_out, grd, shift_lon, shift_lat)
    

    ! First check the possible shift in x direction. This affects the
    ! finite differencing: if grids are shifted, the derivative is
    ! defined naturally in between two points in the "from" grid.
    
    !$OMP PARALLEL DEFAULT(NONE) SHARED(shift_lon, shift_lat, nx, ny, grd, grid_data, ddx, lon_global) &
    !$OMP & PRIVATE(i, j, lat, dist_x, d)

    if(shift_lon .eps. 0.)then

      ! The grids are identical in x direction. Use central
      ! differences. Derivative is defined in the middle of the three
      ! points.
      
      !$OMP DO
      do j = 1, ny
        !lat = real(grd%lonlat%sw_corner_modlat + (real(j-1) * grd%lonlat%dy_deg), vm_p)
        dist_x = 2. * real(fu_dx_cell_m(grd,i,j), vm_p)
        ! The interior of grid.
        do i = 2, nx-1
          dist_x = 2*real(fu_dx_cell_m(grd, i, j), vm_p)
          ddx((j-1)*nx + i) = (grid_data(nx*(j-1) + i + 1) - grid_data(nx*(j-1) + i - 1)) &       
                              / dist_x
        end do
        
        if (lon_global) then
          ! ix = nx
          dist_x = 2 * real(fu_dx_cell_m(grd,nx,j), vm_p)
          ddx(j*nx) = (grid_data((j-1)*nx + 1) - grid_data(j*nx - 1)) / dist_x
          ! ix = 1
          dist_x = 2 * real(fu_dx_cell_m(grd,1,j), vm_p)
          ddx((j-1)*nx + 1) = (grid_data((j-1)*nx + 2) - grid_data((j)*nx)) / dist_x
        else
          ! Boundaries, use forward or backward differences (hence half of dist_x).
          ddx((j-1)*nx + 1) = 2. * (grid_data((j-1)*nx + 2) - grid_data((j-1)*nx + 1)) &
                                   / dist_x
          ddx(j*nx) = 2. * (grid_data(j*nx) - grid_data(j*nx - 1)) / dist_x
        end if

      end do
      !$OMP END DO

    elseif(abs(shift_lon) .eps. 0.5) then
      ! Use central differences, the derivative is defined between two
      ! adjacent points in grd. Depending of the sign of shift, the
      ! first or last cell has to be extrapolated.


      ! Global shifted grids are not handled accurately (as above). Not difficult, but
      ! since there is no data to test it, I'll drop it for the time being.
      if (lon_global) call msg_warning('*** Global shifted grids are handled inaccurately!')

      if (shift_lon .eps. 0.5) then
        d = 1
      else
        d = 0
      end if

      !$OMP DO
      do j = 1, ny
        !lat = grd%lonlat%sw_corner_modlat + (real(j-1+d) * grd%lonlat%dy_deg)
        !dist_x = !fu_dx_deg_to_m(grd%lonlat%dx_deg, lat)
        dist_x = 2. * real(fu_dx_cell_m(grd,i,j+d), vm_p)
        do i = 1, nx - 1
          dist_x = fu_dx_cell_m(grd, i, j)
          ddx((j-1)*nx + i + d) = (grid_data((j-1)*nx + i + 1) - grid_data((j-1)*nx + i)) &
                                  / dist_x
        end do

        ! Check which cell is missing:
        if (d == 1) then
          ddx((j-1)*nx + 1) = ddx((j-1)*nx + 2)
        else
          ddx(j*nx) = ddx(j*nx - 1)
        end if

      end do
      !$OMP END DO

    else
      call set_error('Grids do not match', 'ddx_of_field')
      !return

    endif

    ! Check the shift in y direction. If there is one, the derivatives
    ! are now defined between the gridpoints, and we need to
    ! interpolate. The last value is extrapolated from the two
    ! previous.
    !
    ! To avoid temporary variables, the grid is scanned from north to
    ! south in the case of positive shift, and south to north for the
    ! negative.

    if(shift_lat .eps. 0.)then
      continue
        
    elseif (shift_lat .eps. 0.5)then
      !$OMP DO
      do j = ny, 2, -1
        do i = 1, nx
          ddx((j-1) * nx + i) = 0.5 * (ddx((j-1) * nx + i) + ddx((j-2)*nx + i))
        end do
      end do
      !$OMP END DO
      !$OMP DO
      do i = 1, nx
        ddx(i) = 2. * ddx(nx + i) - ddx(2*nx + i)
      end do
      !$OMP END DO
    elseif (shift_lat .eps. -0.5)then
      !$OMP DO
      do j = 1, ny-1
        do i = 1, nx
          ddx((j-1) * nx + i) = 0.5 * (ddx((j-1) * nx + i) + ddx(j*nx + i))
        end do
      end do
      !$OMP END DO
      !$OMP DO
      do i = 1, nx
        ddx((ny-1)*nx + i) = 2. * ddx((ny-2)*nx + i) - ddx((ny-2)*nx + i)
      end do
      !$OMP END DO
    else
      call set_error('Grids do not match', 'ddx_of_field')
      !return

    endif
    
    !$OMP END PARALLEL
    
  END SUBROUTINE ddx_of_field_dp

  SUBROUTINE div_y_of_field_dp(grid_data, grd, grd_out, ddy)
    !
    ! Calculates d/dx of field (which has grid grd)
    ! in every gridpoint of grid grd_out.
    !
    ! The two grids must either be exactly the same or match one
    ! of the Arakawa grid systems (half a gridpoint difference in either
    ! direction).
    !
    ! Positive derivative to Silja's positive X-direction eastwards.
    ! 
    ! Unit is [unit-of-grid-data]/m.
    !
    ! Handles a +- half-cell shift in either direction.
    !
    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL(vm_p), DIMENSION(:) :: grid_data ! the data itself
    TYPE(silja_grid), INTENT(in) :: grd ! the grid the data is in
    TYPE(silja_grid), INTENT(in) :: grd_out ! the grid in which' points
    ! we want to have the values of derivatives

    ! Imported parameters with intent(out):
    REAL(vm_p), DIMENSION(:), INTENT(out) :: ddy

    ! Local declarations:
    REAL(vm_p) :: lat, dist_y
    REAL(vm_p) :: lat1, dist_x1
    REAL(vm_p) :: lat2, dist_x2
    real :: fTmp, shift_lon, shift_lat
    INTEGER :: i, j, ind1, ind2, nx, ny, nxout, nyout, d, ixtr
    real(vm_p) :: metric_m, metric_s, metric_n

    ! Check and initialize.
    
    IF (.NOT.fu_grids_match_closely(grd, grd_out)) THEN
      CALL report(grd)
      CALL report(grd_out)
      CALL set_error('cannot differentiate, no match','ddy_of_field')
      RETURN
    END IF

!    IF (grd%gridtype /= lonlat ) THEN
!      CALL set_error('sorry, can only handle lonlat grids','ddy_of_field')
!      RETURN
!    END IF

    IF (SIZE(grid_data) > SIZE(ddy)) THEN
      CALL set_error('ddy vector too small','ddy_of_field')
      RETURN
    END IF

    CALL grid_dimensions(grd, nx, ny)
    call grid_dimensions(grd_out, nxout, nyout)
    IF (error) RETURN
  
    if(nx /= nxout .or. ny /= nyout) then
      call set_error('Grids are of different size.', 'ddy_of_field')
      return
    end if
    
    call grid_shift_indices(grd_out, grd, shift_lon, shift_lat)

    !$OMP PARALLEL DEFAULT(NONE) SHARED(shift_lon, shift_lat, nx, ny, grd, grid_data, ddy) &
    !$OMP & PRIVATE(i, j, dist_y, d, metric_n, metric_s, metric_m)

    if(shift_lat .eps. 0.)then
      ! The grids are identical in x direction. Use central
      ! differences. Derivative is defined in the middle of the three
      ! points.

!      dist_y = 2. * fu_dy_deg_to_m(grd%lonlat%dy_deg)

      !$OMP DO
      do j = 2, ny-1
        ! The interior of grid.
        do i = 1, nx
          metric_n = real(fu_dx_cell_m(grd, i, j+1), vm_p)
          metric_s = real(fu_dx_cell_m(grd, i, j-1), vm_p)
          metric_m = real(fu_dx_cell_m(grd, i, j), vm_p)
          
          !dist_y = 2*metric_m!fu_dy_cell_m(grd, i, j)
          dist_y = 2*fu_dy_cell_m(grd, i, j)
          
          ddy((j-1)*nx + i) = (  grid_data(nx*(j) + i)*metric_n &
                             & - grid_data(nx*(j-2)+i)*metric_s) &       
                              / (dist_y*metric_m)
        end do
      end do
      !$OMP END DO

      ! Boundaries, use forward or backward differences (hence half of dist_y).
      !$OMP DO
      do i = 1, nx
        dist_y = real(fu_dy_cell_m(grd, i, 1), vm_p)
        metric_n = real(fu_dx_cell_m(grd, i, 2), vm_p)
        metric_s = real(fu_dx_cell_m(grd, i, 1), vm_p)
        ddy(i) = (grid_data(nx+i)*metric_n - grid_data(i)*metric_s) / (dist_y*metric_s)
        
        dist_y = real(fu_dy_cell_m(grd, i, ny), vm_p)
        metric_n = real(fu_dx_cell_m(grd, i, ny), vm_p)
        metric_s = real(fu_dx_cell_m(grd, i, ny-1), vm_p)
        ddy((ny-1)*nx + i) = (  grid_data((ny-1)*nx + i)*metric_n &
                            & - grid_data((ny-2)*nx + i)*metric_s) &
                            / (dist_y*metric_n)
      end do
      !$OMP END DO

    elseif (shift_lat .eps. 0.5)then

      ! Use central differences, the derivative is defined between two
      ! adjacent points in grd. Depending of the sign of shift, the
      ! first or last cell has to be extrapolated.
      ! 
      ! The column j=1 has to be extrapolated (actually, copied).

      !dist_y = fu_dy_deg_to_m(grd%lonlat%dy_deg)
      !$OMP DO
      do j = 1, ny - 1
        do i = 1, nx
          ! Here we take the meric coefs from grid of the differentiated data.
          !
          metric_n = real(fu_dx_cell_m(grd, i, j+1), vm_p)
          metric_m = real(fu_dx_cell_m(grd, i, j), vm_p)
          metric_s = real(fu_dx_cell_m(grd, i, j-1), vm_p)
          dist_y = real(fu_dy_cell_m(grd, i, j), vm_p)
          ddy(j*nx + i) = (  grid_data(j*nx + i)*metric_n &
                         & - grid_data((j-1)*nx + i)*metric_s) &
                          / (dist_y*metric_m)
        end do
      end do
      !$OMP END DO
      !$OMP DO
      do i = 1, nx
        ddy(i) = ddy(nx + i)
      end do
      !$OMP END DO
    elseif (shift_lat .eps. -0.5)then
      ! As above, but copy the column j=ny.
      
!      dist_y = fu_dy_deg_to_m(grd%lonlat%dy_deg)
      !$OMP DO
      do j = 1, ny - 1
        do i = 1, nx
          dist_y = real(fu_dy_cell_m(grd, i, j), vm_p)
          metric_n = real(fu_dx_cell_m(grd, i, j+1), vm_p)
          metric_m = real(fu_dx_cell_m(grd, i, j), vm_p)
          metric_s = real(fu_dx_cell_m(grd, i, j-1)          , vm_p)

          ddy((j-1)*nx + i) = (  grid_data(j*nx + i)*metric_n &
                             & - grid_data((j-1)*nx + i)*metric_s) &
                              / (dist_y*metric_m)
        end do
      end do
      !$OMP END DO
      !$OMP DO
      do i = 1, nx
        ddy((ny-1)*nx + i) = ddy((ny-2)*nx + i)
      end do
      !$OMP END DO

    else
      call set_error('Grids do not match', 'ddy_of_field')
      !return

    endif ! shift_lat

    ! Check the shift in x direction. If there is one, the derivatives
    ! are now defined between the gridpoints, and we need to
    ! interpolate. The last value is extrapolated from the two
    ! previous.
    !
    ! To avoid temporary variables, the grid is scanned from west to
    ! east in the case of negative shift, and east to west for the
    ! positive.

    if(shift_lon .eps. 0.)then
        continue
        
    elseif (shift_lon .eps. 0.5)then
      !$OMP DO
      do j = 1, ny
        do i = nx, 2, -1
          ddy((j-1) * nx + i) = 0.5 * (ddy((j-1)*nx + i) + ddy((j-1)*nx + i - 1))
        end do
      end do
      !$OMP END DO
      !$OMP DO
      do j = 1, ny
        ddy((j-1)*nx + 1) = 2. * ddy((j-1)*nx + 2) - ddy((j-1)*nx + 3)
      end do
      !$OMP END DO

    elseif (shift_lon .eps. -0.5)then
      !$OMP DO
      do j = 1, ny
        do i = 1, nx - 1
          ddy((j-1) * nx + i) = 0.5 * (ddy((j-1) * nx + i) + ddy((j-1)*nx + i + 1))
        end do
      end do
      !$OMP END DO
      !$OMP DO
      do j = 1, ny
        ddy(j*nx) = 2. * ddy(j*nx - 1) - ddy(j*nx - 2)
      end do
      !$OMP END DO

    else
      call set_error('Grids do not match', 'ddy_of_field')
      !return

    endif
    !$OMP END PARALLEL
    
  END SUBROUTINE div_y_of_field_dp

  SUBROUTINE ddy_of_field_dp(grid_data, grd, grd_out, ddy)
    !
    ! Calculates d/dx of field (which has grid grd)
    ! in every gridpoint of grid grd_out.
    !
    ! The two grids must either be exactly the same or match one
    ! of the Arakawa grid systems (half a gridpoint difference in either
    ! direction).
    !
    ! Positive derivative to Silja's positive X-direction eastwards.
    ! 
    ! Unit is [unit-of-grid-data]/m.
    !
    ! Handles a +- half-cell shift in either direction.
    !
    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL(vm_p), DIMENSION(:) :: grid_data ! the data itself
    TYPE(silja_grid), INTENT(in) :: grd ! the grid the data is in
    TYPE(silja_grid), INTENT(in) :: grd_out ! the grid in which' points
    ! we want to have the values of derivatives

    ! Imported parameters with intent(out):
    REAL(vm_p), DIMENSION(:), INTENT(out) :: ddy

    ! Local declarations:
    REAL(vm_p) :: lat, dist_y
    REAL(vm_p) :: lat1, dist_x1
    REAL(vm_p) :: lat2, dist_x2
    real :: fTmp, shift_lon, shift_lat
    INTEGER :: i, j, ind1, ind2, nx, ny, nxout, nyout, d, ixtr
    real(vm_p) :: metric_m, metric_s, metric_n

    ! Check and initialize.
    
    IF (.NOT.fu_grids_match_closely(grd, grd_out)) THEN
      CALL report(grd)
      CALL report(grd_out)
      CALL set_error('cannot differentiate, no match','ddy_of_field')
      RETURN
    END IF

    IF (SIZE(grid_data) > SIZE(ddy)) THEN
      CALL set_error('ddy vector too small','ddy_of_field')
      RETURN
    END IF

    CALL grid_dimensions(grd, nx, ny)
    call grid_dimensions(grd_out, nxout, nyout)
    IF (error) RETURN
  
    if(nx /= nxout .or. ny /= nyout) then
      call set_error('Grids are of different size.', 'ddy_of_field')
      return
    end if
    
    call grid_shift_indices(grd_out, grd, shift_lon, shift_lat)

    !$OMP PARALLEL DEFAULT(NONE) SHARED(shift_lon, shift_lat, nx, ny, grd, grid_data, ddy) &
    !$OMP & PRIVATE(i, j, dist_y, d, metric_n, metric_s, metric_m)

    if(shift_lat .eps. 0.)then
      ! The grids are identical in x direction. Use central
      ! differences. Derivative is defined in the middle of the three
      ! points.
      !$OMP DO
      do j = 2, ny-1
        ! The interior of grid.
        do i = 1, nx
          
          !dist_y = 2*metric_m!fu_dy_cell_m(grd, i, j)
          dist_y = real(fu_dy_cell_m(grd,i,j-1),vm_p)
          ddy((j-1)*nx + i) = (  grid_data(nx*(j) + i) &
                             & - grid_data(nx*(j-2)+i)) &       
                              / (dist_y)
        end do
      end do
      !$OMP END DO

      ! Boundaries, use forward or backward differences (hence half of dist_y).
      !$OMP DO
      do i = 1, nx
        dist_y = real(fu_dy_cell_m(grd, i, 1), vm_p)
!!$        metric_n = real(fu_dx_cell_m(grd, i, 2), vm_p)
!!$        metric_s = real(fu_dx_cell_m(grd, i, 1), vm_p)
        ddy(i) = (grid_data(nx+i) - grid_data(i)) / (dist_y)
        
        dist_y = real(fu_dy_cell_m(grd, i, ny), vm_p)
!!$        metric_n = real(fu_dx_cell_m(grd, i, ny), vm_p)
!!$        metric_s = real(fu_dx_cell_m(grd, i, ny-1), vm_p)
        ddy((ny-1)*nx + i) = (  grid_data((ny-1)*nx + i) &
                            & - grid_data((ny-2)*nx + i)) &
                            / (dist_y)
      end do
      !$OMP END DO

    elseif (shift_lat .eps. 0.5)then

      ! Use central differences, the derivative is defined between two
      ! adjacent points in grd. Depending of the sign of shift, the
      ! first or last cell has to be extrapolated.
      ! 
      ! The column j=1 has to be extrapolated (actually, copied).

      !dist_y = fu_dy_deg_to_m(grd%lonlat%dy_deg)
      !$OMP DO
      do j = 1, ny - 1
        do i = 1, nx
          ! Here we take the meric coefs from grid of the differentiated data.
          !
          dist_y = real(fu_dy_cell_m(grd, i, j), vm_p)
          ddy(j*nx + i) = (  grid_data(j*nx + i) &
                         & - grid_data((j-1)*nx + i)) &
                          / (dist_y)
        end do
      end do
      !$OMP END DO
      !$OMP DO
      do i = 1, nx
        ddy(i) = ddy(nx + i)
      end do
      !$OMP END DO
    elseif (shift_lat .eps. -0.5)then
      ! As above, but copy the column j=ny.
      
!      dist_y = fu_dy_deg_to_m(grd%lonlat%dy_deg)
      !$OMP DO
      do j = 1, ny - 1
        do i = 1, nx
          dist_y = real(fu_dy_cell_m(grd, i, j), vm_p)
          ddy((j-1)*nx + i) = (  grid_data(j*nx + i) &
                             & - grid_data((j-1)*nx + i)) &
                              / (dist_y)
        end do
      end do
      !$OMP END DO
      !$OMP DO
      do i = 1, nx
        ddy((ny-1)*nx + i) = ddy((ny-2)*nx + i)
      end do
      !$OMP END DO

    else
      call set_error('Grids do not match', 'ddy_of_field')
      !return

    endif ! shift_lat

    ! Check the shift in x direction. If there is one, the derivatives
    ! are now defined between the gridpoints, and we need to
    ! interpolate. The last value is extrapolated from the two
    ! previous.
    !
    ! To avoid temporary variables, the grid is scanned from west to
    ! east in the case of negative shift, and east to west for the
    ! positive.

    if(shift_lon .eps. 0.)then
        continue
        
    elseif (shift_lon .eps. 0.5)then
      !$OMP DO
      do j = 1, ny
        do i = nx, 2, -1
          ddy((j-1) * nx + i) = 0.5 * (ddy((j-1)*nx + i) + ddy((j-1)*nx + i - 1))
        end do
      end do
      !$OMP END DO
      !$OMP DO
      do j = 1, ny
        ddy((j-1)*nx + 1) = 2. * ddy((j-1)*nx + 2) - ddy((j-1)*nx + 3)
      end do
      !$OMP END DO

    elseif (shift_lon .eps. -0.5)then
      !$OMP DO
      do j = 1, ny
        do i = 1, nx - 1
          ddy((j-1) * nx + i) = 0.5 * (ddy((j-1) * nx + i) + ddy((j-1)*nx + i + 1))
        end do
      end do
      !$OMP END DO
      !$OMP DO
      do j = 1, ny
        ddy(j*nx) = 2. * ddy(j*nx - 1) - ddy(j*nx - 2)
      end do
      !$OMP END DO

    else
      call set_error('Grids do not match', 'ddy_of_field')
      !return

    endif
    !$OMP END PARALLEL
    
  END SUBROUTINE ddy_of_field_dp


end module vertical_motion
