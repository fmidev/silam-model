module correlations
  !
  ! correlations : subroutines for correlation-length based 3D spatial correlation models.
  ! 
  ! In the SILAM variational analysis system, the background correlation matrix is defined
  ! implicitly via the transformation x -> Lx, where the correlation matrix equals R =
  ! LL^T. This module defines the forward transformation x -> Lx from control to physical
  ! space, and the adjoint transformation y -> L^Ty, which projects the physical space to
  ! the control space, as required for computing the gradient.
  ! 
  ! The full transformation includes also the diagonal matrix of variances, which is not
  ! included here.
  
use grids_geo
use silam_levels
!$use omp_lib

implicit none

private

public t_correlation
public t_spatial_correlation
public set_spatial_corr_from_nl
public transf_fwd
public transf_adj
public test_correlations
public test_correlations_grid_vert_given
public set_total_correlation
public kill_total_correlation
public fu_dimension_physical
public fu_dimension_control
public defined
public set_missing
public fu_dim_control_3d
public fu_is_diagonal
public fu_spatial_corr_none
public set_corr_times
public fu_corr_times
public get_correlation_matrix

interface defined
   module procedure t_correlation_defined
end interface

interface set_missing
   module procedure set_missing_t_correlation
end interface

type t_correlation
   ! A correlation of several species. Since chemical correlations are not considered,
   ! this is just a combination of spatial correlations for each variable. 
   integer :: num_species ! chemical species
   integer :: field_size  ! nx*ny*nz, possibly nz=1
   type(t_spatial_correlation), dimension(:), pointer :: spatial_correlations => null() ! for each species
   real, dimension(:), pointer :: corr_times => null() ! for each species
   logical :: defined = .false.
end type t_correlation

type t_spatial_correlation
   ! Correlation model for a single variable in three space dimensions. 2D case given by
   ! nz=1.
   integer :: method = int_missing
   
   ! grid, vertical - not defined for no-correlation
   type(silja_grid) :: grid
   type(silam_vertical) :: vertical

   ! Data for the decomposition method
   real, dimension(:,:), pointer :: eigvec_x, eigvec_y, eigvec_z
   real, dimension(:), pointer :: sqrt_eigval_x, sqrt_eigval_y, sqrt_eigval_z
   integer :: nev_x, nev_y, nev_z

   ! Data for the diffusion method
   !type(t_diffusion_data) :: diff_data

   logical :: defined = .false.
end type t_spatial_correlation

real, dimension(:), allocatable, target, private, save :: transf_work, transf_work_2, values_corr_ind

! When the decomposition method is used, the correlation matrix may eigenvalues near
! zero. These correspond to the high wavenumbers "forbidden" by the correlation radius, so
! the spectrum is truncated starting from some threshold. If the threshold is too low, the
! iteration will converge slowly, if too high, the solution will start to suffer. 
! 
! In the following, the threshold is chosen so that eigval_max / eigval_min < max_cond_nbr.
! In practice, the condition number can be quite high before causing problems.
real, parameter, private :: max_cond_nbr = 1e12

integer, parameter, public :: corr_decomp = 12001, corr_none = 12000

contains

  subroutine set_total_correlation(spatial_correlations, grid, vertical, num_species, total_correlation)
    implicit none
    type(t_spatial_correlation), dimension(:), pointer :: spatial_correlations
    type(silja_grid), intent(in) :: grid
    type(silam_vertical), intent(in) :: vertical
    integer, intent(in) :: num_species
    type(t_correlation), intent(out) :: total_correlation

    integer :: ind_species, nx, ny, nz

    if (size(spatial_correlations) /= num_species) then
      call set_error('Inconsistent num_species', 'set_total_correlation')
      return
    end if

    do ind_species = 1, num_species
      if (defined(spatial_correlations(ind_species)%grid)) then
        if (.not.(spatial_correlations(ind_species)%grid == grid)) then
          call set_error('Inconsistent grid', 'set_total_correlation')
          return
        end if
      end if
      if (defined(spatial_correlations(ind_species)%vertical)) then
        if (.not.(fu_cmp_verts_eq(spatial_correlations(ind_species)%vertical, vertical))) then
          call msg_warning('Inconsistent vertical', 'set_total_correlation')
          call msg('Spatial correlation for species:',ind_species)
          call report(spatial_correlations(ind_species)%vertical)
          call msg('External vertical:')
          call report(vertical)
          call set_error('Inconsistent vertical', 'set_total_correlation')
          return
        end if
      end if
    end do

    if (defined(grid)) then
      call grid_dimensions(grid, nx, ny)
    else
      ! vertical-only correlation
      nx = 1
      ny = 1
    end if
    if (fu_fails(defined(vertical), 'vertical not defined', 'set_total_correlation')) return
    nz = fu_nbrOfLevels(vertical)

    total_correlation%spatial_correlations => spatial_correlations
    total_correlation%num_species = num_species
    total_correlation%field_size = nx*ny*nz
    nullify(total_correlation%corr_times)
    total_correlation%defined = .true.

  end subroutine set_total_correlation

  subroutine set_corr_times(corr, corr_times)
    implicit none
    type(t_correlation), intent(inout) :: corr
    real, dimension(:), intent(in), target :: corr_times
    
    if (size(corr_times) /= corr%num_species*corr%field_size) then
      call set_error('corr_times has wrong size', 'set_corr_times')
      return
    end if
    corr%corr_times => corr_times

  end subroutine set_corr_times

  function fu_corr_times(corr) result(times)
    type(t_correlation), intent(in) :: corr
    real, dimension(:), pointer :: times

    times => corr%corr_times
    if (fu_fails(associated(corr%corr_times), 'corr_times not associated', 'fu_corr_times')) return

  end function fu_corr_times
  
  subroutine kill_total_correlation(corr)
    implicit none
    type(t_correlation), intent(inout) :: corr

    integer :: ind_species

    if (.not. corr%defined) return
    do ind_species = 1, corr%num_species
      call kill_spatial_correlation(corr%spatial_correlations(ind_species))
    end do

    corr%defined = .false.

  end subroutine kill_total_correlation
  
  logical function t_correlation_defined(corr) result(is_defined)
    implicit none
    type(t_correlation), intent(in) :: corr
    is_defined = corr%defined
  end function t_correlation_defined

  logical function fu_is_diagonal(corr) result(is_diagonal)
    implicit none
    type(t_correlation), intent(in) :: corr

    integer :: ind_species
    
    if (.not. corr%defined) then
      is_diagonal = .true.
      return
    end if
    is_diagonal = .true.
    do ind_species = 1, corr%num_species
      if (corr%spatial_correlations(ind_species)%method /= corr_none) is_diagonal = .false.
    end do

  end function fu_is_diagonal

  subroutine set_missing_t_correlation(corr)
    implicit none
    type(t_correlation), intent(out) :: corr
    ! default initialization
    corr%defined = .false.
  end subroutine set_missing_t_correlation

  subroutine kill_spatial_correlation(corr)
    implicit none
    type(t_spatial_correlation), intent(inout) :: corr

    if (.not. corr%defined) return
    select case (corr%method)
    case (corr_decomp)
      call set_missing(corr%vertical, ifNew=.false.)
      deallocate(corr%eigvec_x, corr%eigvec_y, corr%eigvec_z, &
               & corr%sqrt_eigval_x, corr%sqrt_eigval_y, corr%sqrt_eigval_z)
    case (corr_none)
      continue ! nothing to deallocate
    end select
    corr%defined = .false.

  end subroutine kill_spatial_correlation

  function fu_spatial_corr_none() result(corr)
    ! return a spatial correlation corresponding to corr_none method
    implicit none
    type(t_spatial_correlation) :: corr
 
    corr%method = corr_none
    corr%grid = grid_missing
    call set_missing(corr%vertical, ifNew=.true.)
    corr%defined = .true.
   
  end function fu_spatial_corr_none

  subroutine set_spatial_corr_from_nl(nl, grid, correlation)
    implicit none
    type(Tsilam_namelist), intent(in), target :: nl
    type(silja_grid), intent(in) :: grid
    type(t_spatial_correlation), intent(out) :: correlation
    
    type(Tsilam_namelist), pointer :: nlptr
    
    nlptr => nl

    select case(fu_content(nlptr, 'correlation_method'))
    case ('decomposition')
      if (defined(grid)) then
        call set_corr_from_nl_decomp(nlptr, grid, correlation)
      else
        ! vertical only
        call set_corr_from_nl_1d(nlptr, correlation)
      end if
      correlation%method = corr_decomp
    case ('none')
      correlation%method = corr_none
      correlation%grid = grid_missing
      call set_missing(correlation%vertical, ifNew=.true.)
    case ('')
      call set_error('Missing correlation_method', 'read_cov_config')
      return
    case default 
      call set_error('Bad correlation_method: ' // trim(fu_content(nlptr, 'correlation_method')), &
                   & 'set_corr_from_nl')
      return
    end select

    correlation%defined = .true.
  end subroutine set_spatial_corr_from_nl

  
  !******************************************************************************
  
  subroutine set_corr_from_nl_decomp(nl, grid, correlation)
    ! 
    ! Initialize the correlation data for use with the decomposition method. Assumes given
    ! x/y correlation distances and an explicit vertical correlation matrix. The resulting
    ! correlation matrices for each dimensions are diagonalised separately, and the
    ! decompositions are stored.
    ! 
    ! The correlation matrix may be rank-deficient, in this case only the eigenvectors
    ! corresponding to eigenvalues greater that eigval_lim are stored, and the dimension
    ! of the control space will be less than that of the physical space.
    ! 
    implicit none
    type(Tsilam_namelist), intent(in) :: nl
    type(silja_grid), intent(in) :: grid
    type(t_spatial_correlation), intent(out) :: correlation
    
    real :: distance, dx, dy, unit_conv
    character(len=*), parameter :: sub_name = 'set_corr_from_nl_decomp'
    integer :: status, i, nx, ny
    real, dimension(:), pointer :: gridval
    character(len=fnlen) :: content
    character(len=clen) :: unit

    correlation%grid = grid
    if (error) return
    call grid_dimensions(correlation%grid, nx, ny)

    gridval => fu_work_array()

    content = fu_content(nl, 'correlation_distance_x')
    if (fu_fails(content /= '', 'Missing correlation_distance_x', sub_name)) return
    read(unit=content, fmt=*, iostat=status) distance, unit
    if (fu_fails(status == 0, 'Failed to parse correlation distance: ' // trim(content), sub_name)) return
    if (.not. (unit == 'm' .or. unit == 'deg')) then
      unit_conv = fu_conversion_factor(unit, 'm')
      if (fu_fails(.not. error, 'Cannot convert correlation_distance_x unit', sub_name)) return
      unit = 'm'
    else
      unit_conv = 1.0
    end if
    distance = distance * unit_conv
    if (unit == 'm') then
      dx = fu_dx_cell_m(correlation%grid, nx/2, ny/2)
      call msg('Correlation distance x: (gridcells, m) ', distance, int(distance/dx))
    else if (unit == 'deg') then
      dx = fu_dx_cell_deg(correlation%grid, nx/2, ny/2)
      call msg('Correlation distance x: (gridcells, deg) ', distance, int(distance/dx))
    end if
    gridval(1:nx) = (/(i*dx, i=1, nx)/)
    call get_corr_from_dist(distance, nx, gridval, &
                          & correlation%eigvec_x, correlation%sqrt_eigval_x, correlation%nev_x)
    if (error) return
    call msg('x-dim condition number: ', &
           & (correlation%sqrt_eigval_x(correlation%nev_x) / correlation%sqrt_eigval_x(1))**2)
    call msg('nx, nev(x):', nx, real(correlation%nev_x))

    content = fu_content(nl, 'correlation_distance_y')
    if (fu_fails(content /= '', 'Missing correlation_distance_y', sub_name)) return
    read(unit=content, fmt=*, iostat=status) distance, unit
    if (fu_fails(status == 0, 'Failed to parse correlation distance: ' // trim(content), sub_name)) return
    if (.not. (unit == 'm' .or. unit == 'deg')) then
      unit_conv = fu_conversion_factor(unit, 'm')
      if (fu_fails(.not. error, 'Cannot convert correlation_distance_y unit', sub_name)) return
      unit = 'm'
    else
      unit_conv = 1.0
    end if
    distance = distance * unit_conv
    if (unit == 'm') then
      dy = fu_dy_cell_m(correlation%grid, nx/2, ny/2)
      call msg('Correlation distance y: (gridcells, m) ', distance, int(distance/dy))
    else if (unit == 'deg') then
      dy = fu_dy_cell_deg(correlation%grid, nx/2, ny/2)
      call msg('Correlation distance y: (gridcells, deg) ', distance, int(distance/dy))
    end if
    gridval(1:ny) = (/(i*dy, i=1, ny)/)
    call get_corr_from_dist(distance, ny, gridval, &
                          & correlation%eigvec_y, correlation%sqrt_eigval_y, correlation%nev_y)
    if (error) return
    call msg('y-dim condition number: ', &
           & (correlation%sqrt_eigval_y(correlation%nev_y) / correlation%sqrt_eigval_y(1))**2)
    call msg('ny, nev(y):', ny, real(correlation%nev_y))

    call get_vertical_correlation(nl, correlation%vertical, &
                                & correlation%eigvec_z , correlation%sqrt_eigval_z, correlation%nev_z)
    call free_work_array(gridval)

    correlation%defined = .true.
    correlation%method = corr_decomp

  contains

    !==========================================================
  
    subroutine get_corr_from_dist(distance, dimlen, dimval, eigvec, sqrt_eigval, nev)
      ! Compute the gaussian correlation function and diagonalize.
      implicit none
      real, intent(in) :: distance
      integer, intent(in) :: dimlen
      real, dimension(:), intent(in) :: dimval
      real, dimension(:,:), pointer :: eigvec
      real, dimension(:), pointer :: sqrt_eigval
      integer, intent(out) :: nev

      ! precision. It seems preferable to do the eigenvalue computation in double
      ! precision, even if the rest goes in single.
      integer, parameter :: eigval_kind = r8k
      integer :: i, j, stat, info
      real :: dist_native
      real(eigval_kind), dimension(:,:), allocatable :: corr
      real(eigval_kind), dimension(:), allocatable :: eigval_tmp
      real(eigval_kind), dimension(:), allocatable :: work
      !real, dimension(:), pointer :: work

      allocate(corr(dimlen, dimlen), eigval_tmp(dimlen), work(worksize), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', 'get_corr_from_dist')) return
      call msg('no 2*distance used')
      do i = 1, dimlen
        do j = 1, dimlen
          corr(i,j) = exp(-(dimval(i)-dimval(j))**2 / (distance**2))
        end do
      end do
      
      !work => fu_work_array()
      
      if (eigval_kind == r8k) then
        call dsyev('V', &    ! Compute eigenvalues and eigenvectors
                 & 'U', &    ! Upper part of A is stored (actually both...)
                 & dimlen, &      ! order of the matrix
                 & corr, &   ! the matrix, overwritten with eigenvectors
                 & dimlen, &      ! leading dimension of svec
                 & eigval_tmp, &  ! eigenvalues in ascending order
                 & work, &   ! work array
                 & size(work),&   ! length of work
                 & info)     ! info
      else
        call ssyev('V', &    ! Compute eigenvalues and eigenvectors
                 & 'U', &    ! Upper part of A is stored (actually both...)
                 & dimlen, &      ! order of the matrix
                 & corr, &   ! the matrix, overwritten with eigenvectors
                 & dimlen, &      ! leading dimension of svec
                 & eigval_tmp, &  ! eigenvalues in ascending order
                 & work, &   ! work array
                 & size(work),&   ! length of work
                 & info)     ! info   
      end if

      !call free_work_array(work)

      if (info /= 0) then
        call msg('dsyev info = ', info)
        call set_error('Eigenvalue computation failed', 'compute_svds')
      end if
#ifdef DOUBLE_PRECISION
      call process_eigval(corr, eigval_tmp, &
                        & real(eigval_tmp(dimlen)) / max_cond_nbr, dimlen, eigvec, sqrt_eigval, nev)
#else
      call process_eigval(real(corr,4), real(eigval_tmp,4), &
                        & real(eigval_tmp(dimlen)) / max_cond_nbr, dimlen, eigvec, sqrt_eigval, nev)
#endif
      if (error) return
      deallocate(corr, eigval_tmp, work)
    end subroutine get_corr_from_dist
    
  end subroutine set_corr_from_nl_decomp

  
  !*****************************************************************************
  
  subroutine set_corr_from_nl_1d(nlptr, correlation)
    ! set a vertical correlation operator. Horizontal components are set to 1.
    implicit none
    type(Tsilam_namelist), pointer :: nlptr
    type(t_spatial_correlation), intent(out) :: correlation
    
    integer :: stat
    
    call get_vertical_correlation(nlptr, correlation%vertical, &
                                & correlation%eigvec_z, correlation%sqrt_eigval_z, &
                                & correlation%nev_z)
    if (error) return
    allocate(correlation%eigvec_x(1,1), correlation%sqrt_eigval_x(1), &
           & correlation%eigvec_y(1,1), correlation%sqrt_eigval_y(1), stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', 'set_corr_from_nl_1d')) return
    correlation%eigvec_x = 1.0
    correlation%sqrt_eigval_x = 1.0
    correlation%eigvec_y = 1.0
    correlation%sqrt_eigval_y = 1.0
    correlation%nev_x = 1
    correlation%nev_y = 1
    correlation%grid = grid_missing
    correlation%defined = .true.

  end subroutine set_corr_from_nl_1d
  
  
  !***************************************************************************************
  
  subroutine get_vertical_correlation(nl, vertical, eigvec, sqrt_eigval, nev)
    ! Parse the vertical correlation matrix from the namelist and diagnoalize.
    implicit none
    type(Tsilam_namelist), intent(in), target :: nl
    type(silam_vertical), intent(out) :: vertical
    real, dimension(:,:), pointer :: eigvec
    real, dimension(:), pointer :: sqrt_eigval
    integer, intent(out) :: nev

    type(Tsilam_nl_item_ptr), dimension(:), pointer :: itemptr
    !real, dimension(:,:), allocatable :: work_2d
    !real, dimension(:), pointer :: work, eigval_tmp
    integer, parameter :: eigval_kind = r8k
    real(eigval_kind), dimension(:), allocatable :: eigval_tmp, work
    real(eigval_kind), dimension(:,:), allocatable :: work_2d
    integer :: nz, iz, ind, info, stat, num_items
    character(len=fnlen) :: content
    type(Tsilam_namelist), pointer :: nlptr
    integer :: i

    nlptr => nl
    !work => fu_work_array()

    call set_vertical(nlptr, vertical)
    if (error) return
    nz = fu_NbrOfLevels(vertical)

    nullify(itemptr)
    allocate(work_2d(nz,nz), eigval_tmp(nz), work(worksize), stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', 'init_vertical_correlation')) return
    call get_items(nlptr, 'vert_correlation', itemptr, num_items)
    if (num_items /= fu_nbrOfLevels(vertical)) then
      call msg('Number of vert_correlation items for following vertical', num_items)
      call report(vertical)
      call set_error('Inconsistent vertical', 'init_vertical_correlation')
      return
    end if

    do iz = 1, fu_NbrOfLevels(vertical)
      content = fu_content(itemptr(iz)) 
      read(content, *, iostat=stat) ind, work(1:nz)
      if (stat /= 0) then
        call msg('vert_correlation: ' // trim(fu_content(itemptr(iz))))
        call set_error('Problem reading vert_correlation', 'init_vertical_correlation')
        return
      end if
      ! Read the matrix temporarily. It will be replaced with the eigenvectors.
      !
      work_2d(ind, 1:nz) = work(1:nz)
    end do
    deallocate(itemptr)

    ! Done reading. Compute eigendecomposition.
    !
    if (eigval_kind == r8k) then
      call dsyev('V', 'U', nz, work_2d, nz, eigval_tmp, work, size(work), info)
    else
      call ssyev('V', 'U', nz, work_2d, nz, eigval_tmp, work, size(work), info)
    end if
    if (info /= 0) then
      call msg('dsyev info = ', info)
      call set_error('Eigenvalue computation failed', 'init_vertical_correlation')
      return
    end if

    call process_eigval(real(work_2d, kind(sqrt_eigval)), real(eigval_tmp, kind(sqrt_eigval)), &
                      & real(eigval_tmp(nz) / max_cond_nbr, kind(sqrt_eigval)), nz, eigvec, sqrt_eigval, nev)

    call msg('Z-dim condition number: ', (sqrt_eigval(nev) / sqrt_eigval(1))**2)
    call msg('nz, nev(z):', nz, real(nev))

    deallocate(work_2d, eigval_tmp, work)
    !call free_work_array(work)

  end subroutine get_vertical_correlation

  
  !********************************************************************************************
  
  subroutine process_eigval(eigvec, eigval, eigval_lim, dim_full, &
       & eigvec_trunc, sqrt_eigval_trunc, nev_trunc)
    !
    ! Given eigenvectors and -values, do following:
    ! - truncate the spectrum so lambda_min > eigval_lim
    ! - compute their square roots and align them in the arrays as used by rest of this module.
    implicit none
    real, dimension(:,:), intent(in) :: eigvec
    real, dimension(:), intent(in) :: eigval
    integer, intent(in) :: dim_full
    real, intent(in) :: eigval_lim
    real, dimension(:,:), pointer :: eigvec_trunc
    real, dimension(:), pointer :: sqrt_eigval_trunc
    integer, intent(out) :: nev_trunc

    integer :: ind_ev, stat

    ! eigval = eigenvalues in ascending order. Find the smallest above threshold and
    ! truncate.
!    call msg('dim_full, eigval_lim:',dim_full, eigval_lim)
!    call msg('Eigen values:', eigval(1:dim_full))
    if (fu_fails(eigval(dim_full) > 1e-15, 'The largest eigenvalue vanishes', 'process_eigval')) return
    do ind_ev = dim_full, 1, -1
      !print *, ind_ev, eigval(ind_ev), eigval_lim
      if (eigval(ind_ev) < eigval_lim .or. eigval(ind_ev) < 0) exit
    end do
    !print *, ind_ev
    nev_trunc = dim_full - ind_ev

    allocate(sqrt_eigval_trunc(nev_trunc), eigvec_trunc(nev_trunc, dim_full), stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', 'init_vertical_correlation')) return
    sqrt_eigval_trunc = sqrt(eigval(ind_ev+1:dim_full))
    eigvec_trunc = transpose(eigvec(1:dim_full, ind_ev+1:dim_full))

  end subroutine process_eigval

       
       
  !************************************************************************************


  ! The transformation routines. The summation order is adjusted to allow n^4/3 complexity
  ! instead of n^2.

  subroutine transf_fwd_3d(corr, values_in, values_out)
    ! fwd == control to physical
    implicit none
    real, dimension(:), intent(inout) :: values_in ! will be lost on return
    type(t_spatial_correlation), intent(in) :: corr
    real, dimension(:), intent(out) :: values_out
    
    integer :: nx, ny, nz, stat
    logical :: allocate_work
    
    select case(corr%method)
    case(corr_none)
      ! no correlation, do nothing
      values_out(1:size(values_in)) = values_in
    case (corr_decomp)
      if (defined(corr%grid)) then
        call grid_dimensions(corr%grid, nx, ny)
      else
        ! vertical only
        nx = 1
        ny = 1
      end if

      if (error) return
      nz = fu_NbrOfLevels(corr%vertical)

      call apply_decomp(values_in, values_out, &
                      & corr%eigvec_x, corr%eigvec_y, corr%eigvec_z, &
                      & corr%sqrt_eigval_x, corr%sqrt_eigval_y, corr%sqrt_eigval_z, &
                      & corr%nev_x, corr%nev_y, corr%nev_z, &
                      & nx, ny, nz)

    case default
      call set_error('Strange correlation method', 'transf_fwd')
    end select

    
  contains 

    subroutine apply_decomp(values, values_out, svecx, svecy, svecz, sqrt_svx, sqrt_svy, sqrt_svz, &
                          & nsv_x, nsv_y, nsv_z, nx, ny, nz)
      implicit none
      real, dimension(:), intent(inout) :: values, values_out

      real, dimension(:,:), intent(in) :: svecx, svecy, svecz
      real, dimension(:), intent(in) :: sqrt_svx, sqrt_svy, sqrt_svz
      integer, intent(in) :: nsv_x, nsv_y, nsv_z, nx, ny, nz
      
      integer :: isx, isy, isz, ix, iy, iz, ind_out, ind_in, nn
      real, dimension(:), pointer :: work

      nn = nx*ny*nz

      values_out = 0.0
      
      if (.not. allocated(transf_work)) then
        allocate(transf_work(nn), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed', 'transf_adj_3d')) return
      else if (size(transf_work) < size(values_in)) then
        deallocate(transf_work)
        allocate(transf_work(nn), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed', 'transf_adj_3d')) return
      end if
      
      work => transf_work

      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(isx, isy, isz, ix, iy, iz, ind_out, ind_in)

      !$OMP DO
      do isx = 1, nsv_x
        do isy = 1, nsv_y
          do iz = 1, nz
            ! isz -> iz
            ind_out = (iz-1)*nsv_x*nsv_y + (isy-1)*nsv_x + isx
            values_out(ind_out) = 0.0
            do isz = 1, nsv_z
              ind_in = (isz-1)*nsv_x*nsv_y + (isy-1)*nsv_x + isx
              values_out(ind_out) = values_out(ind_out) + values(ind_in)*svecz(isz, iz)*sqrt_svz(isz)
            end do
          end do
        end do
      end do
      !$OMP END DO

      !$OMP DO
      do isx = 1, nsv_x
        do iz = 1, nz
          do iy = 1, ny
            ! isy -> iy
            ind_out = (iz-1)*nsv_x*ny + (iy-1)*nsv_x + isx
            work(ind_out) = 0.0
            do isy = 1, nsv_y
              ind_in = (iz-1)*nsv_x*nsv_y + (isy-1)*nsv_x + isx
              work(ind_out) = work(ind_out) + values_out(ind_in)*svecy(isy, iy) * sqrt_svy(isy)
            end do
          end do
        end do
      end do
      !$OMP END DO

      !$OMP DO
      do ix = 1, nx
        do iy = 1, ny
          do iz = 1, nz
            ! silam indexing for single-species slice
            ind_out = iz + (ix-1)*nz + (iy-1)*nx*nz
            !ind_out = (iz-1)*nx*ny + (iy-1)*nx + ix
            values_out(ind_out) = 0.0
            do isx = 1, nsv_x
              ind_in = (iz-1)*nsv_x*ny + (iy-1)*nsv_x + isx
              values_out(ind_out) = values_out(ind_out) + work(ind_in)*svecx(isx, ix)*sqrt_svx(isx)
            end do
          end do
        end do
      end do
      !$OMP END DO

      !$OMP END PARALLEL
     end subroutine apply_decomp
  end subroutine transf_fwd_3d

  subroutine transf_adj_3d(corr, values_in, values_out)
    ! physical -> control, dim(values_in) >= dim(values_out)
    implicit none
    type(t_spatial_correlation), intent(in) :: corr
    real, dimension(:), intent(inout) :: values_in, values_out

    integer :: nx, ny, nz, stat
    !real, dimension(:), pointer :: mywork

    select case(corr%method)
    case(corr_none)
      ! no correlation, do nothing
      values_out(1:size(values_in)) = values_in
    case (corr_decomp)
      if (defined(corr%grid)) then
        call grid_dimensions(corr%grid, nx, ny)
      else
        nx = 1
        ny = 1
      end if
      nz = fu_NbrOfLevels(corr%vertical)
        
      !call msg('sum before:', sum(values_in(1:nx*ny*nz)))
      call apply_decomp(values_in, values_out, &
                      & corr%eigvec_x, corr%eigvec_y, corr%eigvec_z, &
                      & corr%sqrt_eigval_x, corr%sqrt_eigval_y, corr%sqrt_eigval_z, &
                      & corr%nev_x, corr%nev_y, corr%nev_z, &
                      & nx, ny, nz)
      !call msg('sum after:', sum(values_out(1:corr%nev_x*corr%nev_y*corr%nev_z)))
    case default
      call set_error('Strange correlation method', 'transf_fwd')
    end select

  contains

    subroutine apply_decomp(values, values_out, svecx, svecy, svecz, sqrt_svx, sqrt_svy, sqrt_svz, &
                          & nsv_x, nsv_y, nsv_z, nx, ny, nz)
      real, dimension(:), intent(in) :: values   !nx*ny*nz
      real, dimension(:), intent(inout) :: values_out !nsx*nsy*nsz
      real, dimension(:,:), intent(in) :: svecx, svecy, svecz
      real, dimension(:), intent(in) :: sqrt_svx, sqrt_svy, sqrt_svz
      integer, intent(in) :: nsv_x, nsv_y, nsv_z, nx, ny, nz
      
      integer :: isx, isy, isz, ix, iy, iz, ind_out, ind_in, nn, nn_phys
      real, dimension(:), pointer :: work1, work2

      nn = nsv_x*nsv_y*nsv_z
      nn_phys = nx*ny*nz

      if (.not. allocated(transf_work)) then
        call msg('Allocating transformation work array')
        allocate(transf_work(nn_phys), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed', 'transf_adj_3d')) return
      else if (size(transf_work) < size(values)) then
        deallocate(transf_work)
        call msg('Reallocating transformation work array')
        allocate(transf_work(nn_phys), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed', 'transf_adj_3d')) return
      end if

      if (.not. allocated(transf_work_2)) then
        call msg('Allocating transformation work array 2')
        allocate(transf_work_2(nn_phys), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed', 'transf_adj_3d')) return
      else if (size(transf_work_2) < size(values)) then
        deallocate(transf_work_2)
        call msg('Reallocating transformation work array 2')
        allocate(transf_work_2(nn_phys), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed', 'transf_adj_3d')) return
      end if

      work1 => transf_work
      work2 => transf_work_2
      
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(isx, isy, isz, ix, iy, iz, ind_in, ind_out)


      !$OMP DO
      do isz = 1, nsv_z
        do ix = 1, nx
          do iy = 1, ny
            ind_out = (isz-1)*nx*ny + (iy-1)*nx + ix
            work1(ind_out) = 0.0
            do iz = 1, nz
              ! original TDV indexing
              !ind_in = (iz-1)*nx*ny + (iy-1)*nx + ix
              ! silam indexing for single-species slice
              ind_in = iz + (ix-1)*nz + (iy-1)*nx*nz
              work1(ind_out) = work1(ind_out) + values(ind_in)*svecz(isz,iz)*sqrt_svz(isz)
            end do
          end do
        end do
      end do
      !$OMP END DO

      ! indexing now (ix,iy,isz)
      
      !$OMP DO
      do isy = 1, nsv_y
        do ix = 1, nx
          do isz = 1, nsv_z
            ind_out = (isz-1)*nx*nsv_y + (isy-1)*nx + ix
            work2(ind_out) = 0.0
            do iy = 1, ny
              ind_in = (isz-1)*nx*ny + (iy-1)*nx + ix
              work2(ind_out) = work2(ind_out) + work1(ind_in)*svecy(isy,iy)*sqrt_svy(isy)
            end do
          end do
        end do
      end do
      !$OMP END DO
      ! now (ix,isy,isz)
      
      !$OMP DO
      do isx = 1, nsv_x
        do isz = 1, nsv_z
          do isy = 1, nsv_y
            ind_out = (isz-1)*nsv_x*nsv_y + (isy-1)*nsv_x + isx
            values_out(ind_out) = 0.0
            do ix = 1, nx
              ind_in = (isz-1)*nx*nsv_y + (isy-1)*nx + ix
              values_out(ind_out) = values_out(ind_out) + work2(ind_in)*svecx(isx,ix)*sqrt_svx(isx)
            end do
          end do
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

    end subroutine apply_decomp
  end subroutine transf_adj_3d

  subroutine transf_fwd(corr, values_contr, values_phys)
    implicit none
    type(t_correlation), intent(in) :: corr
    real, dimension(:), intent(in), target :: values_contr
    real, dimension(:), intent(out), target  :: values_phys

    integer :: ind_species, ind_start_contr, ind_start_phys
    real, dimension(:), pointer :: species_vals_contr, species_vals_phys
     
    if (fu_fails(size(values_phys) == corr%num_species*corr%field_size, 'Size mismatch p', 'transf_fwd')) return
    if (fu_fails(size(values_contr) == fu_dimension_control(corr), 'Size mismatch c', 'transf_fwd')) return

    ind_start_contr = 1
    ind_start_phys = 1
    do ind_species = 1, corr%num_species
      ! There is a different number of values for each species in physical and control
      ! space. Also, the count in control space varies between species. Hence:
      species_vals_contr => values_contr(ind_start_contr:)

      ! stride = corr%num_species
      !species_vals_phys => values_phys(ind_start_phys:)
      species_vals_phys => values_phys(ind_species::corr%num_species)

      ind_start_contr = ind_start_contr + fu_dim_control_3d(corr, ind_species)
      ind_start_phys = ind_start_phys + corr%field_size
      call transf_fwd_3d(corr%spatial_correlations(ind_species), &
                       & species_vals_contr(1:fu_dim_control_3d(corr, ind_species)), &
                       & species_vals_phys(1:corr%field_size))
      if (error) return
    end do

  end subroutine transf_fwd
  
  !*******************************************************************************
  
  subroutine transf_adj(corr, values_in, values_out)
    ! from physical to control
    implicit none
    type(t_correlation), intent(in) :: corr
    real, dimension(:), intent(in), target :: values_in
    real, dimension(:), intent(out), target :: values_out

    integer :: ind_species, ind_start_contr, ind_start_phys
    real, dimension(:), pointer :: species_vals_contr, species_vals_phys

    ind_start_phys = 1
    ind_start_contr = 1
    do ind_species = 1, corr%num_species
      species_vals_contr => values_out(ind_start_contr:)
 
      !species_vals_phys => values_in(ind_start_phys:)
      species_vals_phys => values_in(ind_species::corr%num_species)

      !call msg('size(values_in)', size(values_in))
      !call msg('ind_species', ind_species)
      !call msg('corr%num_species', corr%num_species)
      !call msg('size(values_out)', size(values_out))
      !call msg('ind_start_contr', ind_start_contr)

      ind_start_contr = ind_start_contr + fu_dim_control_3d(corr, ind_species)
      ind_start_phys = ind_start_phys + corr%field_size
      call transf_adj_3d(corr%spatial_correlations(ind_species), &
                       & species_vals_phys(1:corr%field_size), &
                       & species_vals_contr)
      if (error) return
    end do
    
  end subroutine transf_adj

  !************************************************************************************

  subroutine get_correlation_matrix(spat_corr, matrix)
    ! Evaluate the elements of correlation matrix as R = LL^T . I.
    ! Not a good idea except for 1D or very small 2D correlations.
    implicit none
    type(t_spatial_correlation), intent(in) :: spat_corr
    real, dimension(:,:) :: matrix

    integer :: dim_phys, ii, stat
    real, dimension(:), allocatable :: workarr
    character(len=*), parameter :: sub_name = 'get_correlation_matrix'

    if (fu_fails(spat_corr%defined, 'Correlation not defined', sub_name)) return
    
    if (spat_corr%method == corr_none) then
      ! this is a bit dirty, but I don't know the proper dimension!
      matrix = 0.0
      do ii = 1, min(size(matrix, 1), size(matrix,2))
        matrix(ii,ii) = 1.0
      end do
      return
    end if

    if (fu_fails(defined(spat_corr%vertical), 'Correlation has no vertical', sub_name)) return
    if (defined(spat_corr%grid)) then
      dim_phys = fu_number_of_gridpoints(spat_corr%grid) * fu_NbrOfLevels(spat_corr%vertical)
    else
      dim_phys = fu_NbrOfLevels(spat_corr%vertical)
    end if
    if (dim_phys > 1000) then
      call set_error('Cannot produce correlation matrix', sub_name)
      return
    end if
    
    allocate(workarr(dim_phys), stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
    
    if (fu_fails(size(matrix,1) >= dim_phys, 'Matrix too small (i)', sub_name)) return
    if (fu_fails(size(matrix,2) >= dim_phys, 'Matrix too small (j)', sub_name)) return

    do ii = 1, dim_phys
      matrix(1:dim_phys,ii) = 0.0
      matrix(ii,ii) = 1.0
      call transf_adj_3d(spat_corr, matrix(1:dim_phys,ii), workarr)
      call transf_fwd_3d(spat_corr, workarr, matrix(1:dim_phys,ii))
    end do

    deallocate(workarr)
  
  end subroutine get_correlation_matrix

  !*******************************************************************************

  integer function fu_dimension_control(corr_4d) result(num_dims)
    ! Gives the dimension of the control variable with a given 4d-correlation. This may be
    ! smaller than the physical dimension due to truncation of the spectrum.
    implicit none
    type(t_correlation), intent(in) :: corr_4d
    
    type(t_spatial_correlation) :: corr
    integer :: ind_species

    if (fu_fails(corr_4d%defined, 'Correlation not defined', 'fu_dimension_control')) return

    num_dims = 0
    do ind_species = 1, corr_4d%num_species
      num_dims = num_dims + fu_dim_control_3d(corr_4d, ind_species)
      if (error) then
        num_dims = int_missing
        return
      end if
    end do
      
  end function fu_dimension_control

  integer function fu_dim_control_3d(corr_4d, ind_species) result(num_dims)
    type(t_correlation), intent(in) :: corr_4d
    integer, intent(in) :: ind_species
    
    type(t_spatial_correlation) :: corr_3d

    corr_3d = corr_4d%spatial_correlations(ind_species)
    select case(corr_3d%method)
    case (corr_decomp)
      num_dims = size(corr_3d%sqrt_eigval_x)*size(corr_3d%sqrt_eigval_y)*size(corr_3d%sqrt_eigval_z)
    case (corr_none)
      num_dims = corr_4d%field_size
    case default
      call set_error('strange spatial corr method', 'fu_dimension_control')
    end select
    
  end function fu_dim_control_3d

  integer function fu_dimension_physical(corr) result(num_dims)
    implicit none
    type(t_correlation), intent(in) :: corr

    integer :: nx, ny

    if (fu_fails(corr%defined, 'Correlation not defined', 'fu_dimension_physical')) return
    
    num_dims = corr%field_size * corr%num_species

  end function fu_dimension_physical

  !************************************************************************************

  subroutine test_correlations(distance, nlptr, correlation, surface_only)
    implicit none
    real, intent(in) :: distance
    type(Tsilam_namelist), pointer :: nlptr
    type(t_spatial_correlation), intent(out) :: correlation
    logical, intent(in) :: surface_only

    type(silja_grid) :: grid
    integer :: ilev, nn, nc
    character(len=80) :: distance_str
    type(t_spatial_correlation), dimension(:), pointer :: corr_list_ptr
    type(t_correlation) :: total_correlation

    real, dimension(:), allocatable :: val1, val2, val3, val1p, val2c

    write(distance_str, '(F8.1, A)') distance, ' deg'

    nlptr => fu_create_namelist()
    call add_namelist_item(nlptr, 'correlation_method', 'decomposition')
    call add_namelist_item(nlptr, 'correlation_distance_x', distance_str)
    call add_namelist_item(nlptr, 'correlation_distance_y', distance_str)
    
    if (surface_only) then
      call add_namelist_item(nlptr, 'vertical_method', 'SURFACE_LEVEL')
      call add_namelist_item(nlptr, 'vert_correlation', '1 1.0')
    else
      call add_namelist_item(nlptr, 'number_of_levels', '3')
      call add_namelist_item(nlptr, 'vertical_method', 'CUSTOM_LAYERS')
      call add_namelist_item(nlptr, 'layer_thickness', '50. 100. 150.')
      call add_namelist_item(nlptr, 'level_type', 'HEIGHT_FROM_SURFACE')

      call add_namelist_item(nlptr, 'vert_correlation', '1 1.0 0.5 0.1')
      call add_namelist_item(nlptr, 'vert_correlation', '2 0.5 1.0 0.2')
      call add_namelist_item(nlptr, 'vert_correlation', '3 0.1 0.2 1.0')
    end if
    !call add_namelist_item(nlptr, 'vert_correlation', '1 1.0 0.0 0.0')
    !call add_namelist_item(nlptr, 'vert_correlation', '2 0.0 1.0 0.0')
    !call add_namelist_item(nlptr, 'vert_correlation', '3 0.0 0.0 1.0')


    grid = fu_set_lonlat_grid('testgrid', 20.0, 40.0, .true., 20, 20, pole_geographical, 0.5, 0.5)

    if (error) return
    
    call set_corr_from_nl_decomp(nlptr, grid, correlation)
    if (error) return
    
    allocate(corr_list_ptr(1))
    corr_list_ptr(1) = correlation
    call set_total_correlation(corr_list_ptr, grid, correlation%vertical, 1, total_correlation)

    nn = fu_dimension_physical(total_correlation)
    nc = fu_dimension_control(total_correlation)
    allocate(val1(nn), val2(nn), val3(nn), val1p(nn), val2c(nn))
    call random_number(val1)
    call random_number(val2)
    val3 = val1
    call transf_fwd(total_correlation, val1, val1p)
    call msg('left:', sum(val1p*val2))
    call transf_adj(total_correlation, val2, val2c)
    call msg('right:', sum(val3(1:nc)*val2c(1:nc)))
    ! note values of val1, val2 are now invalid!
  end subroutine test_correlations

  subroutine test_correlations_grid_vert_given(distance, nlptr, correlation, grid, vertical)
    implicit none
    real, intent(in) :: distance

    type(t_spatial_correlation), intent(out) :: correlation
    type(silja_grid), intent(in) :: grid
    type(silam_vertical), intent(in) :: vertical

    integer :: ilev, nn, nc, tempfile_unit, iostat, nlevs
    character(len=255) :: distance_str, corr_line
    
    real, dimension(:), allocatable :: val1, val2, val3, val1p, val2c
    real, dimension(:), pointer :: corr
    type(Tsilam_namelist), pointer :: nlptr, nl_vertical
    type(t_correlation) :: total_correlation
    type(t_spatial_correlation), dimension(:), pointer :: corr_list_ptr

    write(distance_str, '(F8.1, A)') distance, ' deg'

    !nlptr => fu_create_namelist()


    tempfile_unit = fu_next_free_unit()
    open(unit=tempfile_unit, status='SCRATCH', iostat=iostat)
    if (fu_fails(iostat == 0, 'Tempfile creation failed', 'test_correlations_grid_vert_given')) return
    call report_as_namelist(vertical, tempfile_unit)
    rewind(tempfile_unit)
    nl_vertical => fu_read_namelist(tempfile_unit, .false.)
    if (fu_fails(associated(nl_vertical), 'Failed to read vert namelist', 'test_correlations_etc')) return
    close(tempfile_unit)

    nlptr => nl_vertical
    call add_namelist_item(nlptr, 'correlation_method', 'decomposition')
    call add_namelist_item(nlptr, 'correlation_distance_x', distance_str)
    call add_namelist_item(nlptr, 'correlation_distance_y', distance_str)
    
    
    !call add_namelist_item(nlptr, 'number_of_levels', '3')
    !call add_namelist_item(nlptr, 'vertical_method', 'CUSTOM_LAYERS')
    !call add_namelist_item(nlptr, 'layer_thickness', '50. 100. 150.')
    !call add_namelist_item(nlptr, 'level_type', 'HEIGHT_FROM_SURFACE')
    nlevs = fu_NbrOfLevels(vertical)
    corr => fu_work_array()
    do ilev = 1, nlevs
      corr(1:nlevs) = 0.0
      corr(ilev) = 1.0
      if (ilev > 1) corr(ilev-1) = 0.5
      if (ilev < nlevs) corr(ilev+1) = 0.5
      write(corr_line, fmt='(I4, 100F4.1)') ilev, corr(1:nlevs)
      call add_namelist_item(nlptr, 'vert_correlation', corr_line)
    end do
    call free_work_array(corr)
    !grid = fu_set_lonlat_grid('testgrid', 20.0, 40.0, .true., 20, 20, pole_geographical, 0.5, 0.5)

    if (error) return
    
    call set_corr_from_nl_decomp(nlptr, grid, correlation)
    if (error) return
    if (fu_fails(fu_cmp_verts_eq(correlation%vertical, vertical), 'Vertical changed', 'test_correlations')) return

    allocate(corr_list_ptr(1))
    corr_list_ptr(1) = correlation
    call set_total_correlation(corr_list_ptr, grid, vertical, 1, total_correlation)
    nn = fu_dimension_physical(total_correlation)
    nc = fu_dimension_control(total_correlation)
    allocate(val1(nn), val2(nn), val3(nn), val1p(nn), val2c(nn))
    call random_number(val1)
    call random_number(val2)
    val3 = val1
    call transf_fwd(total_correlation, val1, val1p)
    ! val1 is lost now
    call msg('left:', sum(val1p*val2))
    call transf_adj(total_correlation, val2, val2c)
    ! val2 now lost
    call msg('right:', sum(val3(1:nc)*val2c(1:nc)))
  end subroutine test_correlations_grid_vert_given


end module correlations
