module observations_vertical
!  use dispersion_server
!  use optical_density
  use observations_in_situ  !, only : observation_scaling

  implicit none

!  private

  public test_vert_obs
  public test_vert_kern
  public set_vert_obs_from_nc

!  public observe_lidar

  interface inject
     module procedure inject_vertical
  end interface
  public inject

  !interface observe
  !   module procedure observe_vertical
  !end interface
  !public observe

  public observe_vertical


  interface destroy
     module procedure destroy_vertical
  end interface
  public destroy

  interface fu_size
     module procedure fu_size_vertical
  end interface
  public fu_size

  interface obs_to_file
     module procedure obs_to_file_vertical
  end interface
  public obs_to_file

  interface set_data
     module procedure set_data_vertical
  end interface
  public set_data

  interface get_data
     module procedure get_data_vertical
  end interface
  public get_data

  public restart_vertical
  
  interface set_missing
     module procedure set_missing_vertical
  end interface
  public set_missing
  
  interface defined
     module procedure vert_obs_defined
  end interface
  public defined

  interface get_localisation
     module procedure get_localisation_vertical
  end interface
  public get_localisation

  type t_vertical_observation
     ! Vertical averaging kernel. The first 2 dimensions are iz_disp and iz_obs. The last
     ! dimension depends on type: for regular observations it is the pixel, for AOD it is
     ! (transport) species.  The kernel depends on pixel also for AOD, but it is never
     ! read from file but calculated as it gets used.
     real, dimension(:,:,:), pointer :: vert_kernel ! (ilev_disp, ilev_obs, ind_col)
     ! Chemical kernel: for aod, just mask of the contributing species, for others, also
     ! unit conversion. The unit conversion of AOD is handle by the optical module. 
     real, dimension(:), pointer :: chem_kernel
     real, dimension(:), pointer :: lon, lat ! (ind_col)
     real, dimension(:,:), pointer :: values, values_mdl, variance
     real, dimension(:,:), pointer :: weight, disp_area! ind_coef, ind_col
     integer, dimension(:,:), pointer :: ind_x, ind_y ! ind_coef, ind_col
     ! indices in dispersion grid for taking met data
     integer, dimension(:), pointer :: ind_x_mean, ind_y_mean
     ! optical species for each observed transport species
     integer, dimension(:), pointer :: isp_opt 

     ! if averaging kernel is given in time-dependent vertical, it needs to be stored in
     ! the original form and projected as it gets used.
     real, dimension(:,:,:), pointer :: vert_kern_orig ! (ind_kern_lev, ind_obs_lev, ind_col)
     real, dimension(:,:), pointer :: vert_kern_levs_orig ! (ind_kern_lev, ind_col)
     
     ! if the averaging kernel is OK or it has to be refined using the current meteo
     logical :: vert_kernel_valid

     integer :: num_columns = int_missing
     integer :: num_levels = int_missing ! in observation space
     integer :: num_coefs = int_missing

     logical :: defined = .false.

     character(len=clen) :: label = ''
     character(len=clen) :: vert_kern_unit = ''
     type(silja_time) :: valid_time = time_missing
     
     logical :: is_aod = .false.
     logical :: is_lidar = .false.
     ! if kernel vertical is given in Pa, whether scale the pressures to match model surface pressure:
     logical :: scale_to_surface
  end type t_vertical_observation
  public t_vertical_observation

  character(len=*), parameter, public :: var_name_aod = 'aot'
  character(len=*), parameter, public :: var_name_lidar = 'lidar'

  type(tVertInterpStruct), private, save, pointer :: pVertInterp4LevBnds => null()

contains
  
  subroutine create_vert_obs(val_3d, std_3d, lon_2d, lat_2d, nx, ny, &
                           & chem_kernel, vert_kernel, &
                           & have_projected_vert_kern, &
                           & kern_lev_bnds_orig, &
                           & valid_time, label, vert_kern_unit, is_aod, &
                           & is_lidar, disp_grid, interp_method, obs)
    implicit none
    real, dimension(:,:,:), intent(in) :: val_3d, std_3d ! (ix,iy,ilev)
    real, dimension(:) :: lon_2d, lat_2d
    integer, intent(in) :: nx, ny, interp_method
    ! not used for AOD, chem dimension handled through vert_kernel.
    real, dimension(:), intent(in) :: chem_kernel 
    ! in model or orginal vertical, (ilev_kern, ilev_obs, ix, iy)
    real, dimension(:,:,:,:), intent(in) :: vert_kernel 
    logical, intent(in) :: have_projected_vert_kern
    real, dimension(:,:,:), intent(in) :: kern_lev_bnds_orig ! not used if have_projected_vert_kern
    character(len=*) :: label
    type(silja_grid), intent(in) :: disp_grid
    type(t_vertical_observation), intent(out) :: obs
    type(silja_time), intent(in) :: valid_time
    character(len=*), intent(in) :: vert_kern_unit ! the original unit, irrelevant for preprojected kernels
    logical, intent(in) :: is_aod, is_lidar

    type(silja_grid) :: frame_grid
    type(THorizInterpStruct), pointer :: interp
    type(ThorizInterpCells) :: interp_coefs
    integer :: num_valid, num_contrib, num_levs_obs, num_levs_disp, ilev, ind_coef, ind_disp_1d, ind_in, &
         & ix, iy,  ind_out, ix_disp, iy_disp
    integer :: nx_disp, ny_disp, stat, nz_kern_orig
    logical, dimension(:,:), allocatable :: valid
    real :: fx, fy

    nullify(interp)

    frame_grid = fu_set_any_grid(label, nx, ny, lon_2d, lat_2d)
    if (error) return
    call grid_dimensions(disp_grid, nx_disp, ny_disp)

    select case (interp_method)
    case (average)                  ! beware: turned off the randomization
      interp => fu_horiz_interp_struct(frame_grid, disp_grid, interp_method, .false., ifMakeRotation = .false.)
    case (linear)
      interp => fu_horiz_interp_struct(disp_grid, frame_grid, interp_method, .false., ifMakeRotation = .false.)
    case default
      call set_error('Only linear or average interpolation methods allowed', 'create_vert_obs')
      return
    end select
    if (error) return

    num_levs_obs = size(vert_kernel, 2)
    num_levs_disp = nz_dispersion !size(vert_kernel, 1) bad if kernel in orig vertical

    num_valid = 0
    allocate(valid(nx,ny), stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', 'create_vert_obs')) return
    valid = .false.
    do iy = 1, ny
      do ix = 1, nx
        ind_in = (iy-1) * nx + ix
        call project_point_to_grid(lon_2d(ind_in), lat_2d(ind_in), disp_grid, fx, fy)
        if (fx < 1 .or. fx > nx_disp .or. fy < 1 .or. fy > ny_disp) cycle
        if (any(val_3d(ix,iy,:) .eps. real_missing)) cycle
        num_valid = num_valid + 1
        valid(ix,iy) = .true.
      end do
    end do

    if (num_valid == 0) then
      call msg_warning('No valid data for observation: ' // trim(label) &
                     & // trim(fu_str(valid_time,.false.)))
      call set_missing(obs)
      call cleanup()
      return
    end if
    
    num_contrib = fu_nCoefs(interp)

    allocate(obs%values(num_levs_obs, num_valid), &
           & obs%values_mdl(num_levs_obs, num_valid), obs%variance(num_levs_obs, num_valid), &
           & obs%weight(num_contrib, num_valid), &
           & obs%ind_x(num_contrib, num_valid), obs%ind_y(num_contrib, num_valid), &
           & obs%lon(num_valid), obs%lat(num_valid), &
           & obs%disp_area(num_contrib, num_valid), obs%ind_x_mean(num_valid), obs%ind_y_mean(num_valid), &
           & stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', 'create_vert_obs')) return

    ! chem kernel size == num_transp_species
    allocate(obs%chem_kernel(size(chem_kernel)), stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', 'create_vert_obs')) return
    if (is_aod) then
      ! for aod, vertical kernel has last dimension equal to number of species. The kernel
      ! is recomputed for each point on the fly - for each species.
      allocate(obs%vert_kernel(num_levs_disp, num_levs_obs, size(chem_kernel)), &
             & obs%isp_opt(size(chem_kernel)), stat=stat)
    else
      ! for non-aod, vertical kernel holds the data for all observation points at a time,
      ! but does not vary for species.
      allocate(obs%vert_kernel(num_levs_disp, num_levs_obs, num_valid), stat=stat)
      nullify(obs%isp_opt)
    end if
    if (fu_fails(stat == 0, 'Allocate failed', 'create_vert_obs')) return
    if (.not. have_projected_vert_kern) then
      ! averaging kernel has not been pre-projected and the original kernel must be
      ! stored.
      obs%vert_kernel_valid = .false.
      nz_kern_orig = size(vert_kernel, 1)
      allocate(obs%vert_kern_orig(nz_kern_orig, num_levs_obs, num_valid), &
             & obs%vert_kern_levs_orig(nz_kern_orig+1, num_valid), stat=stat)
      obs%vert_kern_unit = vert_kern_unit
      obs%scale_to_surface = vert_kern_unit == 'Pa'
    else
      nullify(obs%vert_kern_orig, obs%vert_kern_levs_orig)
      obs%vert_kernel_valid = .true.
      obs%vert_kern_unit = char_missing
    end if
    if (fu_fails(stat == 0, 'Allocate failed', 'create_vert_obs')) return

    call get_coefs(interp, interp_coefs)
    ind_out = 0
    do iy = 1, ny
      do ix = 1, nx
        if (.not. valid(ix,iy)) cycle
        ind_in = (iy-1) * nx + ix
        ind_out = ind_out + 1
          
        obs%values(:,ind_out) = val_3d(ix,iy,:)
        obs%variance(:,ind_out) = std_3d(ix,iy,:)**2
        
        if (.not. is_aod) then
          if (have_projected_vert_kern) then
            obs%vert_kernel(:,:,ind_out) = vert_kernel(:,:,ix,iy)
          else
            obs%vert_kern_orig(:,:,ind_out) = vert_kernel(:,:,ix,iy)
            obs%vert_kern_levs_orig(:,ind_out) = kern_lev_bnds_orig(:,ix,iy)
            obs%vert_kernel(:,:,ind_out) = real_missing
          end if
        end if
        
        obs%ind_x(:,ind_out) = interp_coefs%indx(:,ix,iy)
        obs%ind_y(:,ind_out) = interp_coefs%indy(:,ix,iy)
        if (any(obs%ind_y(:,ind_out) < 1)) then
          print *, obs%ind_y(:,ind_out)
          print *, ix, iy
          print *, interp_coefs%indy(:,ix,iy)
          call set_error('Bad ind_y', 'here')
          return
        end if

        obs%weight(:,ind_out) = interp_coefs%weight(:,ix,iy)
        do ind_coef = 1, num_contrib
          obs%disp_area(ind_coef, ind_out) = fu_cell_size(disp_grid, &
                                                        & obs%ind_x(ind_coef,ind_out), &
                                                        & obs%ind_y(ind_coef,ind_out))
        end do
        obs%lon(ind_out) = lon_2d(ind_in)
        obs%lat(ind_out) = lat_2d(ind_in)
          
        ind_disp_1d = fu_grid_index(nx_disp, ix, iy, interp)
        iy_disp = ind_disp_1d / nx_disp
        if (mod(ind_disp_1d, nx_disp) > 0) iy_disp = iy_disp + 1
        ix_disp = ind_disp_1d - (iy_disp-1) * nx_disp
        obs%ind_x_mean(ind_out) = ix_disp
        obs%ind_y_mean(ind_out) = iy_disp

        !call msg('coefs: ix_disp, iy_disp', ix_disp, iy_disp)
        !call msg('ind_1d', ind_disp_1d)
        

      end do
    end do
    if (fu_fails(ind_out == num_valid, 'Counts don''t match', 'create_vert_obs')) return
    obs%valid_time = valid_time
    obs%chem_kernel = chem_kernel
    obs%label = label
    obs%is_aod = is_aod
    obs%is_lidar = is_lidar
    obs%num_levels = num_levs_obs
    obs%num_columns = num_valid
    obs%num_coefs = num_contrib
    if (.not. is_aod) obs%defined = .true.
    
    call msg('Observation defined: ' // trim(label) // ' - valid/total:', num_valid, nx*ny)
    call msg('Nonzero: ', count(obs%values > 0))
    call cleanup()

    
!!$    num_levs_obs = size(val_3d, 1)
!!$    do iy = 1, ny
!!$      do ix = 1, nx
!!$        do ilev = 1, num_levs_obs
!!$          if (.not. (val_3d(ilev,ix,iy) .eps. real_missing)) then
!!$            column_valid(ix,iy) = .true.
!!$        end do
!!$      end do
!!$    end do
  contains
    
    subroutine cleanup()
      implicit none
      if (allocated(valid)) deallocate(valid)
      if (associated(interp)) call release_horiz_interp_struct(interp)
      if (defined(frame_grid)) call release_grid(frame_grid)
    end subroutine cleanup

  end subroutine create_vert_obs
  
  
  !******************************************************************************************
  
  subroutine set_aod_obs(obs, wavelength, species_transp, species_opt)
    implicit none
    type(t_vertical_observation), intent(inout) :: obs
    type(silam_species), dimension(:), intent(in) :: species_transp, species_opt
    !type(TspeciesReference), dimension(:), intent(in) :: transp_to_opt_refs
    real, intent(in) :: wavelength

    integer :: isp_transp, isp_opt, n_selected
    integer, dimension(:), pointer :: indices

    if (fu_fails(obs%is_aod, 'Cannot set_aod_obs for non-aod obs', 'set_aod_obs')) return
    
    ! Need to find the optical species for this wavelength:
    indices => fu_work_int_array()
    do isp_transp = 1, size(species_transp)
      call report(species_transp(isp_transp))
      call msg('wave:', wavelength)
      call select_species(species_opt, size(species_opt), &
                        & fu_substance_name(species_transp(isp_transp)), &
                        & fu_mode(species_transp(isp_transp)), wavelength, indices, n_selected)
      if (n_selected == 0) then
        obs%isp_opt(isp_transp) = int_missing
      else if (n_selected > 1) then
        call set_error('Several optical species match', 'set_aod_obs')
        return
      else
        obs%isp_opt(isp_transp) = indices(1)
        call msg('set_aod_obs, transport species:')
        call report(species_transp(isp_transp))
        call msg('set_aod_obs, optical species:')
        call report(species_opt(indices(1)))
      end if
    end do

    if (all(obs%isp_opt(:) == int_missing)) then
      call msg('Creating aod observations with wavelength:', wavelength)
      call set_error('No optical species found for observation: ' // trim(obs%label), &
                   & 'set_aod_obs')
      return
    end if
    
    obs%defined = .true.
    call free_work_array(indices)

  end subroutine set_aod_obs
  
  
  !******************************************************************************************
  
  subroutine get_chem_kern(observed_species, transport_species, obs_unit, kernel, if_profile)
    ! observed species < transport_species
    implicit none
    type(silam_species), dimension(:), intent(in) :: observed_species
    type(silam_species), dimension(:), intent(in) :: transport_species
    character(len=*), intent(in) :: obs_unit
    real, dimension(:), pointer :: kernel
    logical , intent(in) :: if_profile

    integer :: num_obs_species, num_transp_species, ind_transp_species, ind_obs_species, &
         & stat, ind_sqr, ind_slash
    real :: conversion, conversion_area
    type(chemical_adaptor) :: adaptor
    type(silam_material), pointer :: material
    character(len=clen) :: area_unit, mass_unit
    character(1) :: pow

    num_obs_species = size(observed_species)
    num_transp_species = size(transport_species)

    allocate(kernel(num_transp_species), stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', 'get_chem_kernel')) return
    
    call create_adaptor(observed_species, transport_species, adaptor)
    if (error) return

    ind_slash = index(obs_unit, '/')
    if (fu_fails(ind_slash > 0, 'Observation unit not per area: ' // trim(obs_unit), 'get_chem_kern')) return
    mass_unit = obs_unit(1:ind_slash-1)
    area_unit = obs_unit(ind_slash+1:)
    ind_sqr = len_trim(area_unit)
    
    pow = area_unit(ind_sqr:ind_sqr)
    if (if_profile) then
      if (fu_fails(pow == '3', 'Volume unit must end with 3', 'get_chem_kern')) return
    else  
      if (fu_fails(pow == '2', 'Area unit must end with 2', 'get_chem_kern')) return
    end if
    area_unit = area_unit(1:ind_sqr-1)

    ! All conversions are from model unit to observed unit!
    !
    conversion_area = fu_conversion_factor('m', area_unit)**2
    if (error) return

    kernel = 0.0
    do ind_obs_species = 1, num_obs_species
      ind_transp_species = adaptor%isp(ind_obs_species)
      material => fu_material(transport_species(ind_transp_species))
      conversion = fu_conversion_factor(fu_basic_mass_unit(material), mass_unit, material) * conversion_area
      if (error) return
      call msg('Observed substance, conversion ' // trim(fu_name(material)), conversion)
      kernel(ind_transp_species) = conversion
    end do

  end subroutine get_chem_kern

  
  !******************************************************************************************
  
  subroutine test_vert_kern()
    implicit none
    integer, parameter :: nzdisp=5, nzkern=4, nzobs=2
    real, dimension(nzdisp+1) :: disp_vert_bnds_pa = (/1000.0, 950.0, 850.0, 700.0, 500.0, 350.0/), &
         & disp_vert_bnds_m = (/0.0, 100.0, 200.0, 400.0, 700.0, 1000.0/)
    real, dimension(nzkern+1) :: kern_vert_bnds_pa = (/1013.0, 800.0, 550.0, 300.0, 150.0/), &
         & kern_vert_bnds_m = (/20.0, 400.0, 800.0, 1400.0, 2000.0/)
    real, dimension(nzobs, nzkern) :: kern_values
    real, dimension(nzdisp, nzobs) :: kernel

    integer :: ilev_kern, ilev_disp, ilev_obs


    print *, 'Testing av kernel: Pa'
    kern_values = 0.0
    kern_values(1,1:2) = 1.0
    kern_values(2,3:4) = 1.0
    call get_vert_kern(disp_vert_bnds_pa, kern_vert_bnds_pa, kern_values, 'Pa', kernel)
    
    do ilev_disp = 1, nzdisp
      do ilev_obs = 1, nzobs
        print *, 'ilev_disp, ilev_obs, value', ilev_disp, ilev_obs, kernel(ilev_disp, ilev_obs)
      end do
    end do

    print *, 'Testing av kernel: m'
    kern_values = 0.0
    kern_values(1,1:2) = 1.0
    kern_values(2,3:4) = 1.0
    call get_vert_kern(disp_vert_bnds_m, kern_vert_bnds_m, kern_values, 'm', kernel)
    
    do ilev_disp = 1, nzdisp
      do ilev_obs = 1, nzobs
        print *, 'ilev_disp, ilev_obs, value', ilev_disp, ilev_obs, kernel(ilev_disp, ilev_obs)
      end do
    end do

    print *, 'Testing av kernel for column: m'

    kern_values(1,:) = (/0.2, 0.3, 0.4, 0.5/)
    call get_vert_kern(disp_vert_bnds_m, kern_vert_bnds_m, kern_values(1:1,:), 'm', kernel(:,1:1))
    
    print *, 'Kernel in original vert'
    do ilev_kern = 1, nzkern
      print *, 'meters, value:', 0.5 * (kern_vert_bnds_m(ilev_kern)+kern_vert_bnds_m(ilev_kern+1)), &
           & kern_values(1,ilev_kern)
    end do

    print *, 'Kernel in disp vert'
    do ilev_disp = 1, nzdisp
      print *, 'meters, value:', 0.5 * (disp_vert_bnds_m(ilev_disp)+disp_vert_bnds_m(ilev_disp+1)), &
           & kernel(ilev_disp, 1)
    end do

  end subroutine test_vert_kern

  
  !******************************************************************************************
  
  subroutine get_vert_kern(disp_vert_bnds_kern_unit, kern_vert_bnds, kern_values, kern_unit, kernel)
    ! Get the vertical averaging kernel from dispersion grid to observed values (1 for
    ! column, many for profile).
    ! 
    ! The subroutine processes an individual datapoint with given vertical bounds and
    ! kernel values.
    implicit none
    real, dimension(:), intent(in) :: disp_vert_bnds_kern_unit, kern_vert_bnds
    real, dimension(:,:), intent(in) :: kern_values ! (ind_obs_lev, ind_kern_lev)
    real, dimension(:,:), intent(out) :: kernel ! (ind_disp_lev, ind_obs_lev)
    character(len=*), intent(in) :: kern_unit

    real, dimension(:,:), pointer :: disp_to_kern
    integer :: ilev_kern, ilev_disp
    integer :: nz_kern, nz_disp, start_disp, end_disp, step, nz_obs
    real :: bottom_disp, top_disp, bottom_kern, top_kern
    real :: thickness_disp, match_top, match_bottom, fract_to_lev
    integer :: ind_bottom, ind_mid, ind_top, start_kern, end_kern

    nz_disp = size(disp_vert_bnds_kern_unit) - 1
    nz_kern = size(kern_vert_bnds) - 1

    !print *, '*** get_vert_kern ***'
    !print *, 'kern_vert_bnds:', kern_vert_bnds
    !print *, 'kern_values', kern_values

    ! Dispersion to kernel: vertical fractions
    !
    disp_to_kern => fu_work_array_2d()
    disp_to_kern(1:nz_kern, 1:nz_disp) = 0.0

    if (kern_unit == 'm') then
      do ilev_disp = 1, nz_disp
        bottom_disp = disp_vert_bnds_kern_unit(ilev_disp)
        top_disp = disp_vert_bnds_kern_unit(ilev_disp + 1)
        thickness_disp = top_disp - bottom_disp
        
        if (bottom_disp > kern_vert_bnds(nz_kern+1)) then
          call set_error('Dispersion vertical is outside the observed vertical', 'get_vert_kern')
          return
        end if
        if (top_disp < kern_vert_bnds(1)) cycle

        do ilev_kern = 1, nz_kern
          bottom_kern = kern_vert_bnds(ilev_kern)
          top_kern = kern_vert_bnds(ilev_kern + 1)
          match_bottom = max(bottom_disp, bottom_kern)
          match_top = min(top_disp, top_kern)
          if (match_top > match_bottom) then
            fract_to_lev = (match_top-match_bottom) / thickness_disp
            disp_to_kern(ilev_kern, ilev_disp) = fract_to_lev
          !else
          !  disp_to_kern(ilev_kern, ilev_disp) = 0.0
          end if
        end do
      end do

    else if (kern_unit == 'Pa') then
      
      do ilev_disp = 1, nz_disp
        bottom_disp = disp_vert_bnds_kern_unit(ilev_disp)
        top_disp = disp_vert_bnds_kern_unit(ilev_disp + 1)
        thickness_disp = top_disp - bottom_disp

        if (bottom_disp < kern_vert_bnds(nz_kern+1)) then
          call set_error('Dispersion vertical is outside the observed vertical', 'get_vert_kern')
          return
        end if
        if (top_disp > kern_vert_bnds(1)) cycle

        do ilev_kern = 1, nz_kern
          bottom_kern = kern_vert_bnds(ilev_kern)
          top_kern = kern_vert_bnds(ilev_kern + 1)
          match_bottom = min(bottom_disp, bottom_kern)
          match_top = max(top_disp, top_kern)
          if (match_top < match_bottom) then
            fract_to_lev = (match_top-match_bottom) / thickness_disp ! negative/negative
            disp_to_kern(ilev_kern, ilev_disp) = fract_to_lev
          !else
          !  disp_to_kern(ilev_kern, ilev_disp) = 0.0
          end if
        end do
      end do

    else
      call set_error('Bad kernel vertical unit', 'get_vert_kern')
      return
    end if

    do ilev_disp = 1, nz_disp
      do ilev_kern = 1, nz_kern
        if (disp_to_kern(ilev_kern, ilev_disp) == 0.0) cycle
        !print *, 'ilev_disp, ilev_kern', ilev_disp, ilev_kern
        !print *, 'Input:', disp_vert_bnds_kern_unit(ilev_disp), disp_vert_bnds_kern_unit(ilev_disp + 1)
        !print *, 'Output:', kern_vert_bnds(ilev_kern), kern_vert_bnds(ilev_kern + 1)
        !print *, 'Fraction:', disp_to_kern(ilev_kern, ilev_disp)
      end do
    end do

    ! Kernel to observation: matrix product
    !
    nz_obs = size(kernel, 2)
    !kernel(1:nz_obs, 1:nz_disp) = matmul(kern_values, disp_to_kern(1:nz_kern,1:nz_disp))
    kernel(1:nz_disp, 1:nz_obs) = transpose(matmul(kern_values, disp_to_kern(1:nz_kern,1:nz_disp)))
    call free_work_array(disp_to_kern)

  end subroutine get_vert_kern

  
  !******************************************************************************************
  
  subroutine refine_vert_kern(obs, ptr_press, ptr_hgt, ptr_surf_press, &
                            & weight_past, pHorizInterp, ifHorizInterp)
    implicit none
    type(t_vertical_observation), intent(inout) :: obs
    type(field_4d_data_ptr), pointer :: ptr_press, ptr_hgt
    type(field_2d_data_ptr), pointer :: ptr_surf_press
    real, intent(in) :: weight_past
    type(THorizInterpStruct), pointer :: pHorizInterp
    logical, intent(in) :: ifHorizInterp

    type(field_4d_data_ptr), pointer :: ptr_lev_data
    logical :: have_meters
    real, dimension(:), pointer :: column
    integer :: ix, iy, i1d, ind_col
    real :: model_sp, kern_sp

    if (obs%vert_kernel_valid) return
    
    call msg('Refining averaging kernel')
    
    ! check for the interp struct and create if needed
    if (.not. associated(pVertInterp4LevBnds)) then
      call get_interp()
      if (error) return
    end if
    
    ! check the kernel vertical type
    select case(obs%vert_kern_unit)
    case ('m')
      have_meters = .true.
      ptr_lev_data => ptr_hgt
    case ('Pa')
      ptr_lev_data => ptr_press
      have_meters = .false.
    case default
      call set_error('Bad vert_kern_unit for obs', 'refine_vert_kern')
      return
    end select
      
    column => fu_work_array()
    do ind_col = 1, obs%num_columns
      if (error) cycle
      ix = obs%ind_x_mean(ind_col)
      iy = obs%ind_y_mean(ind_col)
      i1d = (iy-1)*nx_meteo + ix
      call get_column(ptr_lev_data, ix, iy, column(2:))
      if (.not. have_meters) then
        model_sp = fu_get_value(ptr_surf_press, nx_meteo, ix, iy, weight_past, &
                              & pHorizInterp, ifHorizInterp)
        column(1) = model_sp
        if (obs%scale_to_surface) then
          ! if the first kernel layer is directly above surface, scale the pressures to
          ! match the model surface pressure.
          kern_sp = obs%vert_kern_levs_orig(1,ind_col)
          obs%vert_kern_levs_orig(:,ind_col) = obs%vert_kern_levs_orig(:,ind_col) * model_sp / kern_sp
        end if
      else
        column(1) = 0.0
      end if
      call get_vert_kern(column(1:nz_dispersion+1), obs%vert_kern_levs_orig(:,ind_col), &
                       & transpose(obs%vert_kern_orig(:,:,ind_col)), &
                       & obs%vert_kern_unit, &
                       & obs%vert_kernel(:,:,ind_col))
    end do
    call free_work_array(column)

    obs%vert_kernel_valid = .true.

  contains
    
    subroutine get_column(ptr_data, ix, iy, column)
      implicit none
      type(field_4d_data_ptr), intent(in) :: ptr_data
      integer, intent(in) :: ix, iy
      real, dimension(:), intent(out) :: column

      integer :: ilev
      type(THorizInterpStruct), pointer :: null_struct
      
      nullify(null_struct)
      do ilev = 1, nz_dispersion
        column(ilev) = fu_get_value(ptr_data, nx_meteo, ix, iy, ilev, weight_past, &
                                  & pHorizInterp, pVertInterp4LevBnds, &
                                  & ifHorizInterp, ifVertInterp=.true.)
      end do
    end subroutine get_column
    
    subroutine get_interp()
      ! create/request interpolation structure from meteo to dispersion level boundaries
      implicit none
      type(silam_vertical) :: vert_disp_bnds
      
      type(silja_level), dimension(nz_dispersion+1) :: levels
      integer :: ilev

      call msg('Requesting vertical interpolation for refining averaging kernels')
      ! obs the vertical interpolation is for vertical starting from the top boundary!
      do ilev = 1, nz_dispersion
        levels(ilev) = fu_upper_boundary_of_layer(fu_level(dispersion_vertical, ilev))
      end do
      call set_vertical(levels, vert_disp_bnds)
      if (error) return
      call arrange_levels_in_vertical(vert_disp_bnds)
      if (error) return
      pVertInterp4LevBnds => fu_vertical_interp_struct(meteo_vertical, &
                                                     & vert_disp_bnds, &
                                                     & dispersion_grid, &
                                                     & linear, &
                                                     & very_long_interval, &
                                                     & 'vertInterp4LevBnds')
      if (error) return
      call set_missing(vert_disp_bnds, ifNew=.false.)
    end subroutine get_interp

  end subroutine refine_vert_kern

  
  !******************************************************************************************
  
  subroutine get_vert_kern_trivial(vertical, nx, ny, kernel)
    implicit none
    type(silam_vertical), intent(in) :: vertical
    integer, intent(in) :: nx, ny ! of observation
    real, dimension(:,:,:,:), pointer :: kernel

    integer :: nlevs, ilev, stat

    nlevs = fu_NbrOfLevels(vertical)
    allocate(kernel(nlevs, 1, nx, ny), stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', 'get_vert_kern_trivial')) return
    
    kernel = 1.0
    
  end subroutine get_vert_kern_trivial

  
  !******************************************************************************************
  
  subroutine get_vert_kern_dummy(vertical, nx, ny, kernel)
    ! A convenience subroutine for creating AOD observations - in this case the vertical
    ! kernel argument is not used.
    implicit none
    type(silam_vertical), intent(in) :: vertical
    integer, intent(in) :: nx, ny ! of observation
    real, dimension(:,:,:,:), pointer :: kernel

    integer :: nlevs, ilev, stat

    allocate(kernel(0,0,0,0), stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', 'get_vert_kern_trivial')) return
    
  end subroutine get_vert_kern_dummy

  
  !******************************************************************************************
  
  subroutine destroy_vertical(obs)
    implicit none
    type(t_vertical_observation), intent(inout) :: obs
    deallocate(obs%vert_kernel, obs%chem_kernel, obs%values, obs%values_mdl, obs%variance, obs%weight, &
             & obs%ind_x, obs%ind_y, &
             & obs%lon, obs%lat, obs%disp_area, obs%ind_x_mean, obs%ind_y_mean)
    if (associated(obs%vert_kern_orig)) deallocate(obs%vert_kern_orig, obs%vert_kern_levs_orig)
    if (associated(obs%isp_opt)) deallocate(obs%isp_opt)
    call set_missing(obs)
  end subroutine destroy_vertical
  
  
  !******************************************************************************************
  
  subroutine obs_to_file_vertical(obs, species_trn, file_unit)
    implicit none
    type(t_vertical_observation), intent(in) :: obs
    type(silam_species), dimension(:), pointer :: species_trn
    integer, intent(in) :: file_unit
    
    integer :: ind_val, ind_lev
    character(len=*), parameter :: fmt = '(A, " ", A, " ", F7.2, F7.2, I5, G12.3, G12.3, G12.3)'
    
    do ind_val = 1, obs%num_columns
      do ind_lev = 1, obs%num_levels
        write(file_unit, fmt=fmt) trim(obs%label), trim(fu_str(obs%valid_time)), &
             & obs%lon(ind_val), obs%lat(ind_val), ind_lev, &
             & obs%values(ind_lev, ind_val), obs%values_mdl(ind_lev, ind_val), &
             & obs%variance(ind_lev, ind_val)
      end do
    end do
  end subroutine obs_to_file_vertical

  
  !******************************************************************************************
  
  subroutine observe_vertical(obs, map_c, now, timestep, &
                            & p_rel_hum, p_tempr, p_press, p_hgt, p_surf_press, &
                            & met_disp_interp_horiz, if_met_disp_interp_horiz, &
                            & met_disp_interp_vert, if_met_disp_interp_vert, &
                            & weight_past, rules_opt_dens, n_opt_species)
    implicit none
    type(t_vertical_observation), intent(inout) :: obs
    type(Tmass_map), intent(in) :: map_c
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep
    type(field_4d_data_ptr), pointer :: p_rel_hum, p_tempr, p_press, p_hgt
    type(field_2d_data_ptr), pointer :: p_surf_press
    type(THorizInterpStruct), pointer :: met_disp_interp_horiz
    logical, intent(in) :: if_met_disp_interp_horiz
    type(TVertInterpStruct), pointer :: met_disp_interp_vert
    logical, intent(in) :: if_met_disp_interp_vert
    real, intent(in) :: weight_past
    type(Toptical_density_rules), intent(in) :: rules_opt_dens
    integer, intent(in) :: n_opt_species
    !type(TspeciesReference), dimension(:), intent(in) :: transp_to_opt_refs

    real :: mass_transp
    integer, target :: isp_tr, ind_col
    integer :: ilev_obs, ilev_disp
    integer :: num_species_transp
    integer, pointer :: ind_vert_kern
    integer, parameter :: ind_src = 1

    if (fu_fails(obs%defined, 'Undefined vertical observation', 'observe_vertical')) return

    if (.not. fu_between_times(obs%valid_time, now, now+timestep, accept_boundaries_too=.true.)) return
    ! above both boundaries are accepted, but we need to choose one to avoid double
    ! counting. Since observe is called at the end of step, it seems reasonable to include that.

    if (obs%valid_time == now) return

    if (obs%is_aod) then
      ! vertical kernel updated for each point, second dimension is transport
      ! species
      ind_vert_kern => isp_tr
    else
      ! vertical kernel predefined, second dimension is the observed column
      ind_vert_kern => ind_col
    end if
    
    num_species_transp = map_c%nspecies

    call msg('Observe vertical: ' // obs%label)
    
    if (.not. obs%vert_kernel_valid) then 
      call refine_vert_kern(obs, p_press, p_hgt, p_surf_press, weight_past, &
                          & met_disp_interp_horiz, if_met_disp_interp_horiz)
      if (error) return
    end if
    
    !if (.not. obs%is_aod) then
    do ind_col = 1, obs%num_columns
      if (obs%is_aod) call update_vert_kern_aod(obs, ind_col, p_rel_hum, p_tempr, &
                                                & met_disp_interp_horiz, if_met_disp_interp_horiz, &
                                                & met_disp_interp_vert, if_met_disp_interp_vert, &
                                                & weight_past, rules_opt_dens, n_opt_species)
      do ilev_obs = 1, obs%num_levels
        obs%values_mdl(ilev_obs,ind_col) = 0.0
        do isp_tr = 1, num_species_transp
          do ilev_disp = 1, map_c%n3d
            mass_transp = fu_get_mass(ind_col, ilev_disp, isp_tr)
            ! vert_kernel includes conversion to concentration
            ! chem_kernel includes conversion to observed unit
!!$            call msg('mass_transp, chem_kern', mass_transp, obs%chem_kernel(isp_tr))
!!$            call msg('isp_tr, vert_kernel', isp_tr, obs%vert_kernel(ilev_disp, ilev_obs, ind_vert_kern))
            obs%values_mdl(ilev_obs,ind_col) = obs%values_mdl(ilev_obs,ind_col) &
                 & + obs%vert_kernel(ilev_disp, ilev_obs, ind_vert_kern)*obs%chem_kernel(isp_tr) &
                 &   * mass_transp
            if (ind_vert_kern == 1) then
            end if
          end do ! dispersion level
        end do ! transport species
        !call msg('finalval', obs%values_mdl(ilev_obs,ind_col))
      end do ! observation level
    end do ! column
    !end if ! not aod
    
  contains
    
    real function fu_get_mass(ind_col, ilev_disp, isp_tr) result(mass)
      implicit none
      integer, intent(in) :: ind_col, ilev_disp, isp_tr

      integer :: ind_coef, ix, iy
      
      mass = 0.0
      do ind_coef = 1, obs%num_coefs
        ix = obs%ind_x(ind_coef,ind_col)
        iy = obs%ind_y(ind_coef,ind_col)
        !call msg('ix,iy', ix, iy)
        !call msg('weight,area', obs%weight(ind_coef,ind_col), obs%disp_area(ind_coef, ind_col))
        !call msg('value, ilev', map_c%arm(isp_tr,ind_src,ilev_disp,ix,iy), ilev_disp)
        if (iy < 1) then
          print *, 'iy:', iy, ind_coef
          print *, obs%ind_y(ind_coef,ind_col)
          call set_error('Bad iy', 'observe_vertical')
          return
        end if
        mass = mass + map_c%arm(isp_tr,ind_src,ilev_disp,ix,iy) &
             & * obs%weight(ind_coef,ind_col) / obs%disp_area(ind_coef, ind_col)
      end do

    end function fu_get_mass

  end subroutine observe_vertical

  
  !******************************************************************************************
  
  ! subroutine observe_lidar(obs, map_c, now, timestep, &
  !      & p_rel_hum, p_tempr, p_press, p_hgt, p_surf_press, &
  !      & met_disp_interp_horiz, if_met_disp_interp_horiz, &
  !      & met_disp_interp_vert, if_met_disp_interp_vert, &
  !      & weight_past, rules_opt_dens, n_opt_species, &
  !      & ChemRunSetup)


  !   implicit none
  !   type(t_vertical_observation), intent(inout) :: obs
  !   type(Tmass_map), intent(in) :: map_c
  !   type(silja_time), intent(in) :: now
  !   type(silja_interval), intent(in) :: timestep
  !   type(field_4d_data_ptr), pointer :: p_rel_hum, p_tempr, p_press, p_hgt
  !   type(field_2d_data_ptr), pointer :: p_surf_press
  !   type(THorizInterpStruct), pointer :: met_disp_interp_horiz
  !   logical, intent(in) :: if_met_disp_interp_horiz
  !   type(TVertInterpStruct), pointer :: met_disp_interp_vert
  !   logical, intent(in) :: if_met_disp_interp_vert
  !   real, intent(in) :: weight_past
  !   type(Toptical_density_rules), intent(in) :: rules_opt_dens
  !   integer, intent(in) :: n_opt_species

  !   type(TChemicalRunSetup), pointer, optional, intent(in) :: ChemRunSetup
  !   character(len=*), parameter :: sub_name = 'observe_vertical'

  !   real :: mass_transp

  !   real, dimension(:), pointer :: alpha, beta, beta_lidrat!, alpha
  !   real, dimension(:), pointer :: ExtCoef, BackScatterCoef
  !   real, dimension(:), pointer :: thickness, relHumid, tmpr
  !   real, dimension(:), pointer :: alpha_Rayleigh, beta_Rayleigh

  !   integer, target :: isp_tr, ind_col
  !   integer :: ilev_obs, ilev_disp
  !   integer :: num_species_transp
  !   integer, pointer :: ind_vert_kern
  !   integer, parameter :: ind_src = 1
  !   integer :: ix, iy, ilev
  !   integer :: isp_opt, iSpTo
  !   logical :: ifRelHumidDep,ifTDep
  !   integer :: stat
  !   real :: area

    
  !   ExtCoef => fu_work_array()
  !   BackScatterCoef => fu_work_array()
  !   alpha => fu_work_array()
  !   beta => fu_work_array()
  !   beta_lidrat => fu_work_array()

  !   alpha_Rayleigh => fu_work_array()    
  !   beta_Rayleigh => fu_work_array()
    
  !   thickness => fu_work_array()
  !   thickness(1:map_c%n3d) = (/(fu_layer_thickness_m(fu_level(map_c%vertTemplate, ilev_disp)), &
  !        & ilev_disp=1, map_c%n3d)/)
    
  !   do ind_col = 1, obs%num_columns
  !     ix = obs%ind_x_mean(ind_col)
  !     iy = obs%ind_y_mean(ind_col)
      
  !     area = fu_cell_size(dispersion_grid, ix, iy)
      
  !     alpha(1:map_c%n3d) = 0.0
  !     beta(1:map_c%n3d) = 0.0
  !     beta_lidrat(1:map_c%n3d) = 0.0
  !     alpha_Rayleigh(1:map_c%n3d) = 0.0
  !     beta_Rayleigh(1:map_c%n3d) = 0.0

  !     do ilev_disp = map_c%n3d, 1, -1

  !       ExtCoef(1:n_opt_species) = 0.0
  !       BackScatterCoef(1:n_opt_species) = 0.0
        
  !       call compute_optical_dens_coefs(p_rel_hum, p_tempr, n_opt_species, &
  !            & met_disp_interp_horiz, met_disp_interp_vert, &
  !            & if_met_disp_interp_horiz, if_met_disp_interp_vert, &
  !            & weight_past, rules_opt_dens, &
  !            & ix, iy, ilev_disp, ExtCoef, BackScatterCoef)
        
  !       do isp_tr = 1, num_species_transp  
  !         isp_opt = obs%isp_opt(isp_tr)
  !         if (isp_opt /= int_missing) then
              
  !           alpha_0(ilev_disp) = alpha_0(ilev_disp) + ExtCoef(isp_opt) * map_c%arm(isp_tr,ind_src,ilev_disp,ix,iy) / &
  !                &  area
            
  !           !alpha(ilev_disp) = alpha(ilev_disp) + ExtCoef(isp_opt) * map_c%arm(isp_tr,ind_src,ilev_disp,ix,iy) / &
  !           !     &  (area * thickness(ilev_disp))
            
  !           beta(ilev_disp) = beta_Mie(ilev_disp) + BackScatterCoef(isp_opt) * &
  !                  & map_c%arm(isp_tr,ind_src,ilev_disp,ix,iy) / (4*pi*area*thickness(ilev_disp))
            
  !         end if
  !       end do

  !       !Backscatter coefficient from the lidar ratio:
  !       beta_lidrat(ilev_disp) = (1.0/65.0) * alpha_0(ilev_disp) / thickness(ilev_disp)
 
  !       !Extinction and backscatter from Rayleigh scattering:
  !       alpha_Rayleigh = 1.16e-5 * (550e-9/ wavelength)**4.09 * p_press(ind_col, ilev) /p0 * t0/p_tempr(ind_col, ilev)
        
  

  !       beta = beta + alpha_Rayleigh  * 3/(8*3.14159)
  !       !beta = beta_lidrat + alpha_Rayleigh  * 3/(8*3.14159)

  !      if (ilev_disp .lt. map_c%n3d) then
  !         alpha_0(ilev_disp) = alpha_0(ilev_disp) + alpha_0(ilev_disp + 1)
  !       end if
  !     end do
        
  !     obs%values_mdl(1:map_c%n3d,ind_col) = beta(1:map_c%n3d) * exp(-2.0*alpha_0(1:map_c%n3d))

  !     !print *, 'cumulative AOD', alpha_0(1:map_c%n3d)
  !     !print *, 'beta_ lidrat', beta_lidrat(1:map_c%n3d)
  !     !print *, 'beta_Mie', beta_Mie(1:map_c%n3d)
  !     !print *, 'Mie_lidrat', alpha(1:map_c%n3d)/beta_Mie(1:map_c%n3d) 
  !     !print *, 'transmittance', exp(-alpha_0(1))
  !   end do
  !   call free_work_array(thickness)
  !   call free_work_array(ExtCoef)
  !   call free_work_array(BackScatterCoef)
  !   call free_work_array(alpha_0)
  !   !call free_work_array(alpha)
  !   call free_work_array(alpha_Rayleigh)
  !   call free_work_array(beta_Rayleigh)
  !   call free_work_array(beta_Mie)
  !   call free_work_array(beta_lidrat)
    
  ! end subroutine observe_lidar

  
  !******************************************************************************************
  
  subroutine update_vert_kern_aod(obs, ind_col, p_rel_hum, p_tempr, &
                                & met_disp_interp_horiz, if_met_disp_interp_horiz, &
                                & met_disp_interp_vert, if_met_disp_interp_vert, &
                                & weight_past, rules_opt_dens, n_opt_species)
    type(t_vertical_observation), intent(inout) :: obs
    integer, intent(in) :: ind_col
    type(field_4d_data_ptr), pointer :: p_rel_hum, p_tempr
    type(THorizInterpStruct), pointer :: met_disp_interp_horiz
    logical, intent(in) :: if_met_disp_interp_horiz
    type(TVertInterpStruct), pointer :: met_disp_interp_vert
    logical, intent(in) :: if_met_disp_interp_vert
    real, intent(in) :: weight_past
    type(Toptical_density_rules), intent(in) :: rules_opt_dens
    integer, intent(in) :: n_opt_species
    !type(TspeciesReference), dimension(:), intent(in) :: transp_to_opt_refs

    integer :: ilev, isp_transp, isp_opt
    integer :: nlevs_disp, nsp_transp
    real, dimension(:), pointer :: coefs

    nsp_transp = size(obs%chem_kernel)
    nlevs_disp = size(obs%vert_kernel, 1)

    coefs => fu_work_array()
    do ilev = 1, nlevs_disp
      call compute_optical_dens_coef(p_rel_hum, p_tempr, n_opt_species, &
                                   & met_disp_interp_horiz, met_disp_interp_vert, &
                                   & if_met_disp_interp_horiz, if_met_disp_interp_vert, &
                                   & weight_past, rules_opt_dens, &
                                   & obs%ind_x_mean(ind_col), obs%ind_y_mean(ind_col), &
                                   & ilev, coefs)
      do isp_transp = 1, nsp_transp
        isp_opt = obs%isp_opt(isp_transp)
        if (isp_opt /= int_missing) then
          ! 1 == number of observed layers
          obs%vert_kernel(ilev,1,isp_transp) = coefs(isp_opt)
        else
          obs%vert_kernel(ilev,1,isp_transp) = 0.0 ! not optically contributing species
        end if
      end do
    end do
    call free_work_array(coefs)
    
  end subroutine update_vert_kern_aod

  
  !******************************************************************************************
  
  subroutine inject_vertical(obs, map_c, map_px, map_py, map_pz, &
                           & now, timestep, p_rel_hum, p_tempr, &
                           & met_disp_interp_horiz, if_met_disp_interp_horiz, &
                           & met_disp_interp_vert, if_met_disp_interp_vert, &
                           & weight_past, rules_opt_dens, n_opt_species)
    implicit none
    type(t_vertical_observation), intent(inout) :: obs
    type(Tmass_map), intent(inout) :: map_c, map_px, map_py, map_pz
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep
    type(field_4d_data_ptr), pointer :: p_rel_hum, p_tempr
    type(THorizInterpStruct), pointer :: met_disp_interp_horiz
    logical, intent(in) :: if_met_disp_interp_horiz
    type(TVertInterpStruct), pointer :: met_disp_interp_vert
    logical, intent(in) :: if_met_disp_interp_vert
    real, intent(in) :: weight_past
    type(Toptical_density_rules), intent(in) :: rules_opt_dens
    integer, intent(in) :: n_opt_species
    !type(TspeciesReference), dimension(:), intent(in) :: transp_to_opt_refs

    integer :: ilev_obs, ilev_disp, num_species_transp, ind_src
    real :: inject_val, inject_val_chem_vert
    integer, pointer :: ind_vert_kern
    integer, target :: isp_tr, ind_col
    real, dimension(:), pointer :: thickness

    if (fu_fails(obs%defined, 'Undefined vertical observation', 'inject_vertical')) return

    if (.not. fu_between_times(obs%valid_time, now, now+timestep, accept_boundaries_too=.true.)) return
    ! above both boundaries are accepted, but we need to choose one to avoid double
    ! counting. Since inject is called at the start of (adjoint) step, it seems reasonable
    ! to include that.(remember negative timestep).
    if (obs%valid_time == now+timestep) return

    call msg('Inject vertical: ' // obs%label)

    if (obs%is_aod) then
      ! vertical kernel updated for each point, second dimension is transport
      ! species
      ind_vert_kern => isp_tr
    else
      ! vertical kernel predefined, second dimension is the observation datapoint
      ind_vert_kern => ind_col
    end if

    thickness => fu_work_array()
    thickness(1:map_c%n3d) = (/(fu_layer_thickness_m(fu_level(map_c%vertTemplate, ilev_disp)), &
                              & ilev_disp=1, map_c%n3d)/)
    
    num_species_transp = map_c%nspecies
        
    do ind_col = 1, obs%num_columns
      if (obs%is_aod) call update_vert_kern_aod(obs, ind_col, p_rel_hum, p_tempr, &
                                              & met_disp_interp_horiz, if_met_disp_interp_horiz, &
                                              & met_disp_interp_vert, if_met_disp_interp_vert, &
                                              & weight_past, rules_opt_dens, n_opt_species)
      do ilev_obs = 1, obs%num_levels
        inject_val = (obs%values_mdl(ilev_obs, ind_col) - obs%values(ilev_obs, ind_col)) &
             & * observation_scaling / obs%variance(ilev_obs, ind_col)
        if (inject_val < 0) then
          ind_src = DA_NEGATIVE
          inject_val  = inject_val * DA_NEGT_COEF_OBS
        else
          ind_src = DA_POSITIVE
        end if

        do isp_tr = 1, num_species_transp
          if (obs%chem_kernel(isp_tr) == 0.0) cycle
          do ilev_disp = 1, map_c%n3d
            inject_val_chem_vert = inject_val &
                 & * obs%vert_kernel(ilev_disp,ilev_obs,ind_vert_kern) &
                 & * obs%chem_kernel(isp_tr) * thickness(ilev_disp)
            call inject_horiz(inject_val_chem_vert, ind_src, ind_col, ilev_disp, isp_tr)
          end do ! dispersion level
        end do ! transport species
      end do ! observation level
    end do ! column
    call free_work_array(thickness)

  contains
    
    subroutine inject_horiz(inject_val, ind_src, ind_col, ilev_disp, isp_tr)
      real, intent(in) :: inject_val
      integer, intent(in) :: ind_col, ilev_disp, isp_tr, ind_src

      integer :: ind_coef, ix, iy

      do ind_coef = 1, obs%num_coefs
        ix = obs%ind_x(ind_coef,ind_col)
        iy = obs%ind_y(ind_coef,ind_col)
        map_c%arm(isp_tr, ind_src, ilev_disp, ix, iy) = map_c%arm(isp_tr,ind_src, ilev_disp, ix, iy) &
             & + inject_val*obs%weight(ind_coef,ind_col)
      end do
    end subroutine inject_horiz
  
  end subroutine inject_vertical

  
  !******************************************************************************************
  
  integer function fu_size_vertical(obs) result(obs_size)
    implicit none
    type(t_vertical_observation), intent(in) :: obs
    
    obs_size = obs%num_levels * obs%num_columns
  end function fu_size_vertical
  
  
  !******************************************************************************************
  
  subroutine get_data_vertical(obs, values_obs, values_mdl, variance)
    implicit none
    type(t_vertical_observation), intent(in) :: obs
    real, dimension(:), intent(out), optional :: values_obs, values_mdl, variance

    integer :: ind_col, ind_lev, ind
    
    do ind_col = 1, obs%num_columns
      do ind_lev = 1, obs%num_levels
        ind = (ind_col-1) * obs%num_levels + ind_lev
        if (present(values_obs)) then
          !if (obs%values(ind_lev,ind_col) > 0) call msg('val', obs%values(ind_lev,ind_col))
          values_obs(ind) = obs%values(ind_lev,ind_col)
        end if
        if (present(values_mdl)) then
          !if (obs%values_mdl(ind_lev,ind_col) > 0) call msg('mdlval', obs%values_mdl(ind_lev,ind_col))
          values_mdl(ind) = obs%values_mdl(ind_lev,ind_col)
        end if
        if (present(variance)) then
          variance(ind) = obs%variance(ind_lev,ind_col)
        end if
      end do
    end do

  end subroutine get_data_vertical

  
  !******************************************************************************************
  
  subroutine set_data_vertical(obs, values_obs, values_mdl, variance)
    implicit none
    type(t_vertical_observation), intent(inout) :: obs
    real, dimension(:), intent(in), optional :: values_obs, values_mdl, variance

    integer :: ind_col, ind_lev, ind
    
    do ind_col = 1, obs%num_columns
      do ind_lev = 1, obs%num_levels
        ind = (ind_col-1) * obs%num_levels + ind_lev
        if (present(values_obs)) then
          obs%values(ind_lev,ind_col) = values_obs(ind)
        end if
        if (present(values_mdl)) then
          obs%values_mdl(ind_lev,ind_col) = values_mdl(ind)
        end if
        if (present(variance)) then
          obs%variance(ind_lev,ind_col) = variance(ind)
        end if
      end do
    end do

  end subroutine set_data_vertical

  
  !******************************************************************************************
  
  subroutine get_localisation_vertical(obs, lonlats)
    implicit none
    type(t_vertical_observation), intent(in) :: obs
    real, dimension(:,:), intent(out) :: lonlats
    
    integer :: ind_col, ind_lev, ind
    
    do ind_col = 1, obs%num_columns
      do ind_lev = 1, obs%num_levels
        ind = (ind_col-1) * obs%num_levels + ind_lev
        lonlats(1, ind) = obs%lon(ind_col)
        lonlats(2, ind) = obs%lat(ind_col)
      end do
    end do
    
  end subroutine get_localisation_vertical
  
  
  !******************************************************************************************
  
  subroutine restart_vertical(obs, ifNullifyModel)
    implicit none
    type(t_vertical_observation), intent(inout) :: obs
    logical, intent(in) :: ifNullifyModel

    if(ifNullifyModel) obs%values_mdl = 0.0
  end subroutine restart_vertical

  
  !******************************************************************************************
  
  subroutine set_missing_vertical(obs)
    implicit none
    type(t_vertical_observation), intent(out) :: obs
    nullify(obs%vert_kernel, obs%chem_kernel, obs%lon, obs%lat, obs%values, obs%values_mdl, obs%variance)
    nullify(obs%weight, obs%disp_area, obs%ind_x, obs%ind_y, obs%ind_x_mean, obs%ind_y_mean, obs%isp_opt)
    obs%defined = .false.
  end subroutine set_missing_vertical
  
  !************************************************************************************

  subroutine set_vert_obs_from_nc(filename, var_name, obs_species, transp_species, opt_species, &
                                & time_start, time_end, obs_list, num_obs)
    ! Read the observation data from Netcdf. The assumption is that the observation
    ! cocktail, eg. obs_species has been named elsewhere (this avoids changing the
    ! observation if the model configuration does, or editing the cocktails). 

    use netcdf
    implicit none
    character(len=*), intent(in) :: filename, var_name
    ! contributing observed species + transport species
    type(silam_species), dimension(:), intent(in) :: transp_species
    type(silam_species), dimension(:), pointer :: obs_species, opt_species ! can be null()
    type(silja_time), intent(in) :: time_start, time_end
    type(t_vertical_observation), dimension(:), intent(out) :: obs_list
    integer, intent(out) :: num_obs

    real, dimension(:,:,:), allocatable :: val_3d, std_3d
    real, dimension(:,:,:), pointer :: kern_lev_bnds
    real, dimension(:), pointer :: lon_2d, lat_2d, chem_kernel, disp_lev_bnds
    real, dimension(:,:,:,:), pointer :: vert_kernel
    integer :: ncid, stat, dim_id_time, dim_id_x, dim_id_y, nx, ny, dim_id_z, nz, var_id_kern_levs
    integer :: var_id_val, var_id_var, var_id_time, unit_len, ind_time, var_id_lat, var_id_lon, &
         & label_len, var_id_vert_kern, dim_id_kern_layer
    character(len=*), parameter :: sub_name = 'set_vert_obs_from_nc'
    character(len=clen) :: label,  unit, timestr, kern_unit
    integer, dimension(:), pointer :: seconds_arr
    type(silja_time) :: ref_time, valid_time
    real :: fillValue, wavelength
    logical :: is_aod, have_kern_layers, can_project_kernel, need_vert_kern, is_lidar

    nullify(vert_kernel, kern_lev_bnds, seconds_arr, chem_kernel)

    stat = nf90_open(filename, 0, ncid)
    if (fu_fails(stat == NF90_NOERR, 'Failed to open: ' // trim(filename), sub_name)) return
    
    stat = nf90_inq_dimid(ncid, 'time', dim_id_time)
    if (fu_fails(stat == NF90_NOERR, 'Failed to inquire', sub_name)) return
    stat = nf90_inquire_dimension(ncid, dim_id_time, len=num_obs)
    if (fu_fails(stat == NF90_NOERR, 'Failed to inquire', sub_name)) return
    
    if (num_obs == 0) then
      call msg_warning('File contains no data')
      return
    end if

    if (fu_fails(num_obs <= size(obs_list), 'Too many column observations', sub_name)) return
 
    stat = nf90_inq_dimid(ncid, 'x', dim_id_x)
    if (fu_fails(stat == NF90_NOERR, 'Failed to inquire', sub_name)) return
    stat = nf90_inquire_dimension(ncid, dim_id_x, len=nx)
    if (fu_fails(stat == NF90_NOERR, 'Failed to inquire', sub_name)) return

    stat = nf90_inq_dimid(ncid, 'y', dim_id_y)
    if (fu_fails(stat == NF90_NOERR, 'Failed to inquire', sub_name)) return
    stat = nf90_inquire_dimension(ncid, dim_id_y, len=ny)
    if (fu_fails(stat == NF90_NOERR, 'Failed to inquire', sub_name)) return

    stat = nf90_inq_dimid(ncid, 'z', dim_id_z)
    if (fu_fails(stat == NF90_NOERR, 'Failed to inquire', sub_name)) return
    stat = nf90_inquire_dimension(ncid, dim_id_z, len=nz)
    if (fu_fails(stat == NF90_NOERR, 'Failed to inquire', sub_name)) return
    stat = nf90_inq_dimid(ncid, 'time', dim_id_z)

    allocate(val_3d(nx,ny,nz), std_3d(nx,ny,nz), stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return

    stat = nf90_inq_varid(ncid, var_name, var_id_val)
    if (fu_fails(stat == NF90_NOERR, 'Failed to get variable id: ' // trim(var_name), sub_name)) return
    stat = nf90_inq_varid(ncid, 'stdev', var_id_var)
    if (fu_fails(stat == NF90_NOERR, 'Failed to get variable id: stdev', sub_name)) return

    stat = nf90_inq_varid(ncid, 'lon', var_id_lon)
    if (fu_fails(stat == NF90_NOERR, 'Failed to get variable id: lon', sub_name)) return
    stat = nf90_inq_varid(ncid, 'lat', var_id_lat)
    if (fu_fails(stat == NF90_NOERR, 'Failed to get variable id: lat', sub_name)) return

    stat = nf90_inquire_attribute(ncid, nf90_global, 'label', len=label_len)
    if (fu_fails(stat == NF90_NOERR, 'Failed to inquire label', sub_name)) return
    if (fu_fails(label_len <= clen, 'Label attribute too long', sub_name)) return
    stat = nf90_get_att(ncid, nf90_global, 'label', label)
    if (fu_fails(stat == NF90_NOERR, 'Failed to get label', sub_name)) return
    
    stat = nf90_get_att(ncid, nf90_global, 'silamtime', timestr)
    if (fu_fails(stat == NF90_NOERR, 'Failed to get silamtime', sub_name)) return
    ref_time = fu_io_string_to_time(timestr)
    if (error) return
    
    stat = nf90_inq_varid(ncid, 'time', var_id_time)
    if (fu_fails(stat == NF90_NOERR, 'Failed to get variable id', sub_name)) return
    seconds_arr => fu_work_int_array()
    stat = nf90_get_var(ncid, var_id_time, seconds_arr, count=(/num_obs/))
    if (fu_fails(stat == NF90_NOERR, 'Failed to get seconds', sub_name)) return

    ! Is there an averaging kernel available?
    !
    !need_vert_kern = nz > 1
 
    stat = nf90_inq_varid(ncid, 'kernel_layer_interfaces', var_id_kern_levs)
    have_kern_layers = stat == NF90_NOERR
    if (.not. have_kern_layers) then
      if (nz > 1) then
        need_vert_kern = .true.
      else 
        need_vert_kern = .false.
      end if
      if (fu_fails(.not. need_vert_kern, 'nz > 1 but averaging kernel not given', sub_name)) return
      call get_vert_kern_trivial(dispersion_vertical, nx, ny, vert_kernel)
      var_id_vert_kern = int_missing
      nullify(disp_lev_bnds)
    else
      need_vert_kern = .true.
      stat = nf90_get_att(ncid, var_id_kern_levs, 'units', kern_unit)
      if (fu_fails(stat == NF90_NOERR, 'Failed to get kern_unit', sub_name)) return
      can_project_kernel = kern_unit == 'm' &
           & .and. fu_leveltype(dispersion_vertical) == layer_btw_2_height
      if (can_project_kernel) then
        call msg('Vertical averaging kernel defined in dispersion levels')
        disp_lev_bnds => fu_work_array()
        call get_disp_lev_bnds(kern_unit, disp_lev_bnds)
        if (error) return
        ! allocate zero-size kern_lev_bnds so we can still give it to create_vert_obs
        allocate(vert_kernel(fu_NbrOfLevels(dispersion_vertical), nz, nx, ny), &
               & kern_lev_bnds(0,0,0), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
      else
        call msg('Vertical averaging kernel in foreign units, reprojection is deferred')
        ! projecting the kernel will be delayed first use. The kernel is read in original
        ! vertical, and the corresponding buffers are allocated while reading.
        nullify(vert_kernel, kern_lev_bnds, disp_lev_bnds)
      end if
    end if

    stat = nf90_inquire_attribute(ncid, var_id_val, 'units', len=unit_len)
    if (fu_fails(stat == NF90_NOERR, 'Failed to inquire unit', sub_name)) return
    if (fu_fails(unit_len <= clen, 'Unit attribute too long', sub_name)) return
    stat = nf90_get_att(ncid, var_id_val, 'units', unit)
    if (fu_fails(stat == NF90_NOERR, 'Failed to unit', sub_name)) return
    stat = nf90_get_att(ncid, var_id_val, '_FillValue', fillValue)
    if (fu_fails(stat == NF90_NOERR, 'Failed to get _FillValue', sub_name)) return

    call msg('var_name:' // var_name)
    is_aod = var_name == var_name_aod
    is_lidar = var_name == var_name_lidar

    if (is_aod) then
      stat = nf90_get_att(ncid, nf90_global, 'wavelength', wavelength)
      if (fu_fails(stat == NF90_NOERR, 'Failed to get wavelength', sub_name)) return
      allocate(chem_kernel(size(transp_species)), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
      ! Set chem_kernel to 1.0: all species contribute. Those with no optical counterpart
      ! are ignored, but this is handled later.
      chem_kernel = 1.0
    else if (is_lidar) then
      stat = nf90_get_att(ncid, nf90_global, 'wavelength', wavelength)
      if (fu_fails(stat == NF90_NOERR, 'Failed to get wavelength', sub_name)) return
      chem_kernel = 1.0
    else
      if (fu_fails(unit /= '', 'Unit attribute is empty', sub_name)) return
      call get_chem_kern(obs_species, transp_species, unit, chem_kernel, if_profile=nz > 1)
      if (error) return
    end if

    lon_2d => fu_work_array()
    lat_2d => fu_work_array()

    ! *** loop over observation times ***
    ! 
    do ind_time = 1, num_obs
      valid_time = ref_time + fu_set_interval_sec(real(seconds_arr(ind_time)))
      call msg('Valid time: ' // fu_str(valid_time, .false.))
      if (valid_time < time_start .or. valid_time > time_end) then
        call set_missing(obs_list(ind_time)) 
        cycle
      end if
      stat = nf90_get_var(ncid, var_id_val, val_3d, &
                        & count = (/nx, ny, nz, 1/), start=(/1, 1, 1, ind_time/))
      !print *, nf90_strerror(stat)
      !call inq_var(var_id_val)
      if (fu_fails(stat == NF90_NOERR, 'Failed to read values', sub_name)) return
      stat = nf90_get_var(ncid, var_id_var, std_3d, &
                        & count = (/nx, ny, nz, 1/), start=(/1, 1, 1, ind_time/))
      if (fu_fails(stat == NF90_NOERR, 'Failed to read stdev', sub_name)) return
      
      stat = nf90_get_var(ncid, var_id_lon, lon_2d, &
                        & count = (/nx, ny/), start=(/1, 1, ind_time/))
      if (fu_fails(stat == NF90_NOERR, 'Failed to read lon', sub_name)) return
      stat = nf90_get_var(ncid, var_id_lat, lat_2d, &
                        & count = (/nx, ny/), start=(/1, 1, ind_time/))
      if (fu_fails(stat == NF90_NOERR, 'Failed to read lat', sub_name)) return
      
      if (need_vert_kern) then      
        if (can_project_kernel) then
          call msg('setting kernel')
          call set_aver_kern_disp_vert(disp_lev_bnds, vert_kernel)
        else
          call set_aver_kern_orig_vert(vert_kernel, kern_lev_bnds, &
               & allocate_arrays=.not. associated(vert_kernel))
        end if
        if (error) return
      end if
      
      where (val_3d .eps. fillValue) val_3d = real_missing

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !std_3d = 0.1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call create_vert_obs(val_3d, std_3d, lon_2d, lat_2d, nx, ny, chem_kernel, vert_kernel, &
                         & can_project_kernel .or. .not. need_vert_kern, kern_lev_bnds, &
                         & valid_time, label, kern_unit, is_aod, is_lidar, &
                         & dispersion_grid, linear, obs_list(ind_time))
      if (error) return
      if (obs_list(ind_time)%is_aod) then
        call msg('Creating aod observation, ind_time', ind_time)
        if (fu_fails(associated(opt_species), 'Optical species not associated', sub_name)) return
        call set_aod_obs(obs_list(ind_time), wavelength, transp_species, opt_species)
      else if (.not. defined(obs_list(ind_time))) then
        continue
      else
        call msg('Regular observation...')
      end if
    end do ! ind_time

    if (fu_fails(nf90_close(ncid) == NF90_NOERR, 'Failed to close netcdf', sub_name)) continue

    call free_work_array(lon_2d)
    call free_work_array(lat_2d)
    call free_work_array(seconds_arr)
    if (associated(disp_lev_bnds)) call free_work_array(disp_lev_bnds)
    if (associated(kern_lev_bnds)) deallocate(kern_lev_bnds)
    if (associated(vert_kernel)) deallocate(vert_kernel)
    if (associated(chem_kernel)) deallocate(chem_kernel)

    deallocate(val_3d, std_3d)
    
  contains
    
    subroutine get_disp_lev_bnds(kern_unit, bnds)
      implicit none
      character(len=*), intent(in) :: kern_unit
      real, dimension(:), intent(out) :: bnds
      
      type(silam_vertical) :: vert_in_kern_unit
      integer :: ilev

      select case(kern_unit)
      case ('m')
        call vert_to_metric(dispersion_vertical, vert_in_kern_unit)
        if (error) return
        bnds(1) = fu_level_height(fu_lower_boundary_of_layer(fu_level(vert_in_kern_unit, 1)))
        do ilev = 1, fu_NbrOfLevels(dispersion_vertical)
          bnds(ilev+1) = fu_level_height(fu_upper_boundary_of_layer(fu_level(vert_in_kern_unit, ilev)))
        end do
      case ('Pa')
        call vert_to_pressure(dispersion_vertical, vert_in_kern_unit)
        bnds(1) = fu_pr_level_pressure(fu_lower_boundary_of_layer(fu_level(vert_in_kern_unit, 1)))
        do ilev = 1, fu_NbrOfLevels(dispersion_vertical)
          bnds(ilev+1) = fu_pr_level_pressure(fu_upper_boundary_of_layer(fu_level(vert_in_kern_unit, ilev)))
        end do
      case default
        call set_error('Unsupported kernel vertical unit: ' // trim(kern_unit), sub_name)
        return
      end select

      call set_missing(vert_in_kern_unit, ifNew=.false.)

    end subroutine get_disp_lev_bnds

    !=====================================================
    
    subroutine set_aver_kern_orig_vert(kernels, kern_lev_bnds, allocate_arrays)
      implicit none
      real, dimension(:,:,:,:), pointer :: kernels
      real, dimension(:,:,:), pointer :: kern_lev_bnds
      logical, intent(in) :: allocate_arrays

      integer :: stat, var_id_lev_bnds, var_id_kernel, nz_kern
      character(len=clen) :: kern_type
      character(len=*), parameter :: sub_name = 'set_aver_kern_orig_vert'
      real, dimension(:), pointer :: lev_bnds_1d
      real, dimension(:,:), pointer :: kern_values_2d
      real, dimension(:,:,:,:), allocatable :: kern_values_4d
      real, dimension(:,:,:), allocatable :: lev_bnds_3d
      integer :: ix, iy, iz

      ! inquire kernel type
      stat = nf90_get_att(ncid, nf90_global, 'aver_kern_type', kern_type)
      if (fu_fails(stat == NF90_NOERR, 'Failed to inquire aver_kern_type', sub_name)) return
      select case(kern_type)
      case ('single_vertical_const', 'single_vertical_var', &
           & 'multi_vertical_var')
        stat = nf90_inq_varid(ncid, 'kernel_layer_interfaces', var_id_lev_bnds)
        if (fu_fails(stat == NF90_NOERR, 'Failed to inquire kernel_layer_interfaces', sub_name)) return
        stat = nf90_inq_varid(ncid, 'kernel', var_id_kernel)
        if (fu_fails(stat == NF90_NOERR, 'Failed to inquire kernel', sub_name)) return
        
        stat = nf90_inq_dimid(ncid, 'kernel_layer', dim_id_kern_layer)
        if (fu_fails(stat == NF90_NOERR, 'Failed to inquire kernel_layer dimension', sub_name)) return
        stat = nf90_inquire_dimension(ncid, dim_id_kern_layer, len=nz_kern)
        if (fu_fails(stat == NF90_NOERR, 'Failed to inquire kernel_layer dimension size', sub_name)) return

      case ('single_vertical_trivial')
        stat = nf90_inq_varid(ncid, 'kernel_layer_interfaces', var_id_lev_bnds)
        if (fu_fails(stat == NF90_NOERR, 'Failed to inquire kernel_layer_interfaces', sub_name)) return
        stat = nf90_inq_dimid(ncid, 'kernel_layer', dim_id_kern_layer)
        if (fu_fails(stat == NF90_NOERR, 'Failed to inquire kernel_layer dimension', sub_name)) return
        stat = nf90_inquire_dimension(ncid, dim_id_kern_layer, len=nz_kern)
        if (fu_fails(stat == NF90_NOERR, 'Failed to inquire kernel_layer dimension size', sub_name)) return

      case default
        call set_error('Bad aver_kern_type: ' // trim(kern_type), sub_name)
        return
      end select
      
      if (allocate_arrays) then
        allocate(kernels(nz_kern,nz,nx,ny), kern_lev_bnds(nz_kern+1,nx,ny), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
      end if
      if (fu_fails(associated(kernels), 'kernels not associated', sub_name)) return
      if (fu_fails(associated(kern_lev_bnds), 'kern_lev_bnds not associated', sub_name)) return
      
      select case(kern_type)
      case ('single_vertical_trivial')
        kernels = 0.0
        lev_bnds_1d => fu_work_array()
        stat = nf90_get_var(ncid, var_id_lev_bnds, lev_bnds_1d, count=(/nz_kern+1/), start=(/1,ind_time/))
        if (fu_fails(stat == NF90_NOERR, 'Failed to read level boundaries', sub_name)) return
        call msg('nz, nz_kern', nz, nz_kern)
        if (fu_fails(nz_kern == nz, 'Bad kernel shape for single_vertical_trivial', sub_name)) return
        do iy = 1, ny
          do ix = 1, nx
            do iz = 1, nz
              kernels(iz,iz,ix,iy) = 1.0 / (lev_bnds_1d(iz+1)-lev_bnds_1d(iz))
            end do
            kern_lev_bnds(:,ix,iy) = lev_bnds_1d(1:nz_kern+1)
          end do
        end do
        call free_work_array(lev_bnds_1d)

      case ('single_vertical_const')
        lev_bnds_1d => fu_work_array()
        stat = nf90_get_var(ncid, var_id_lev_bnds, lev_bnds_1d, count=(/nz_kern+1/), start=(/1,ind_time/))
        if (fu_fails(stat == NF90_NOERR, 'Failed to read level boundaries', sub_name)) return
        
        kern_values_2d => fu_work_array_2d()
        stat = nf90_get_var(ncid, var_id_kernel, kern_values_2d, &
                          & count=(/nz_kern, nz/), start=(/1,1,ind_time/))
        if (fu_fails(stat == NF90_NOERR, 'Failed to read kernel', sub_name)) return
              
        do iy = 1, ny
          do ix = 1, nx
            kernels(:,:,ix,iy) = kern_values_2d(1:nz_kern, 1:nz)
            kern_lev_bnds(:,ix,iy) = lev_bnds_1d(1:nz_kern+1)
          end do
        end do
 
        call free_work_array(kern_values_2d)
        call free_work_array(lev_bnds_1d)

      case ('single_vertical_var')
        lev_bnds_1d => fu_work_array()
        stat = nf90_get_var(ncid, var_id_lev_bnds, lev_bnds_1d, count=(/nz_kern+1/), start=(/1,ind_time/))
        if (fu_fails(stat == NF90_NOERR, 'Failed to read level boundaries', sub_name)) return
        stat = nf90_get_var(ncid, var_id_kernel, kernels, &
                          & count=(/nz_kern, nz, nx, ny/), start=(/1,1,1,1,ind_time/))
        if (fu_fails(stat == NF90_NOERR, 'Failed to read kernel', sub_name)) return
        do iy = 1, ny
          do ix = 1, nx
            kern_lev_bnds(:,ix,iy) = lev_bnds_1d(1:nz_kern+1)       
          end do
        end do
        call free_work_array(lev_bnds_1d)

      case ('multi_vertical_var')
        allocate(lev_bnds_3d(nz_kern+1, nx, ny), kern_values_4d(nz_kern, nz, nx, ny), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed', 'set_aver_kern')) return
        stat = nf90_get_var(ncid, var_id_lev_bnds, kern_lev_bnds, count=(/nz_kern+1, nx, ny/), &
                          & start=(/1, 1, 1, ind_time/))
        if (fu_fails(stat == NF90_NOERR, 'Failed to read level boundaries', sub_name)) return
        stat = nf90_get_var(ncid, var_id_kernel, kernels, &
                          & count=(/nz_kern, nz, nx, ny/), start=(/1,1,1,1,ind_time/))
        if (fu_fails(stat == NF90_NOERR, 'Failed to read kernel', sub_name)) return
        
      case default
        call set_error('Unsupported kernel type', 'set_aver_kern')
      end select
      
    end subroutine set_aver_kern_orig_vert

    !=====================================================
    
    subroutine set_aver_kern_disp_vert(disp_lev_bnds, kernels)
      implicit none
      real, dimension(:), intent(in) :: disp_lev_bnds ! kern_unit
      real, dimension(:,:,:,:), intent(out) :: kernels !(ind_disp_lev, ind_obs_lev, ix, iy)
      
      integer :: stat, var_id_lev_bnds, var_id_kernel, nz_kern
      character(len=clen) :: kern_type
      real, dimension(:), pointer :: lev_bnds_1d
      real, dimension(:,:), pointer :: kern_values_2d
      real, dimension(:,:,:,:), allocatable :: kern_values_4d
      real, dimension(:,:,:), allocatable :: lev_bnds_3d
      integer :: ix, iy, iz

      ! inquire kernel type
      stat = nf90_get_att(ncid, nf90_global, 'aver_kern_type', kern_type)
      if (fu_fails(stat == NF90_NOERR, 'Failed to inquire aver_kern_type', sub_name)) return
      select case(kern_type)
      case ('single_vertical_const', 'single_vertical_var', 'multi_vertical_var')
        stat = nf90_inq_varid(ncid, 'kernel_layer_interfaces', var_id_lev_bnds)
        if (fu_fails(stat == NF90_NOERR, 'Failed to inquire kernel_layer_interfaces', sub_name)) return
        stat = nf90_inq_varid(ncid, 'kernel', var_id_kernel)
        if (fu_fails(stat == NF90_NOERR, 'Failed to inquire kernel', sub_name)) return
        
        stat = nf90_inq_dimid(ncid, 'kernel_layer', dim_id_kern_layer)
        if (fu_fails(stat == NF90_NOERR, 'Failed to inquire kernel_layer dimension', sub_name)) return
        stat = nf90_inquire_dimension(ncid, dim_id_kern_layer, len=nz_kern)
        if (fu_fails(stat == NF90_NOERR, 'Failed to inquire kernel_layer dimension size', sub_name)) return
      
      case ('single_vertical_trivial')
        stat = nf90_inq_varid(ncid, 'kernel_layer_interfaces', var_id_lev_bnds)
        if (fu_fails(stat == NF90_NOERR, 'Failed to inquire kernel_layer_interfaces', sub_name)) return
        stat = nf90_inq_dimid(ncid, 'kernel_layer', dim_id_kern_layer)
        if (fu_fails(stat == NF90_NOERR, 'Failed to inquire kernel_layer dimension', sub_name)) return
        stat = nf90_inquire_dimension(ncid, dim_id_kern_layer, len=nz_kern)
        if (fu_fails(stat == NF90_NOERR, 'Failed to inquire kernel_layer dimension size', sub_name)) return

      case default
        call set_error('Bad aver_kern_type: ' // trim(kern_type), sub_name)
      end select
      
      select case(kern_type)
      case ('single_vertical_const')
        lev_bnds_1d => fu_work_array()
        stat = nf90_get_var(ncid, var_id_lev_bnds, lev_bnds_1d, count=(/nz_kern+1/), start=(/1,ind_time/))
        if (fu_fails(stat == NF90_NOERR, 'Failed to read level boundaries', sub_name)) return
        
        kern_values_2d => fu_work_array_2d()
        stat = nf90_get_var(ncid, var_id_kernel, kern_values_2d, &
                          & count=(/nz_kern, nz/), start=(/1,1,ind_time/))
        if (fu_fails(stat == NF90_NOERR, 'Failed to read kernel', sub_name)) return

        do iy = 1, ny
          do ix = 1, nx
            call get_vert_kern(disp_lev_bnds(1:nz_dispersion+1), lev_bnds_1d(1:nz_kern+1), &
                             & transpose(kern_values_2d(1:nz_kern,1:nz)), kern_unit, kernels(:,:,ix,iy))
            if (error) return
          end do
        end do
 
        call free_work_array(kern_values_2d)
        call free_work_array(lev_bnds_1d)


      case ('single_vertical_trivial')
        if (fu_fails(nz_kern == nz, 'Bad kernel shape for single_vertical_trivial', sub_name)) return
        
        lev_bnds_1d => fu_work_array()
        stat = nf90_get_var(ncid, var_id_lev_bnds, lev_bnds_1d, count=(/nz_kern+1/), start=(/1,ind_time/))
        if (fu_fails(stat == NF90_NOERR, 'Failed to read level boundaries', sub_name)) return
        
        kern_values_2d => fu_work_array_2d()
        kern_values_2d(1:nz, 1:nz) = 0

        do iz = 1, nz
          kern_values_2d(iz,iz) = 1.0 / (lev_bnds_1d(iz+1)-lev_bnds_1d(iz))
        end do
        
        do iy = 1, ny
          do ix = 1, nx
            call get_vert_kern(disp_lev_bnds(1:nz_dispersion+1), lev_bnds_1d(1:nz_kern+1), &
                             & kern_values_2d(1:nz_kern,1:nz), kern_unit, kernels(:,:,ix,iy))
            if (error) return
          end do
        end do
        call free_work_array(kern_values_2d)
        call free_work_array(lev_bnds_1d)


      case ('single_vertical_var')
        lev_bnds_1d => fu_work_array()
        stat = nf90_get_var(ncid, var_id_lev_bnds, lev_bnds_1d, count=(/nz_kern+1/), start=(/1,ind_time/))
        if (fu_fails(stat == NF90_NOERR, 'Failed to read level boundaries', sub_name)) return
        allocate(kern_values_4d(nz_kern, nz, nx, ny), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed', 'set_aver_kern')) return
        stat = nf90_get_var(ncid, var_id_kernel, kern_values_4d, &
                          & count=(/nz_kern, nz, nx, ny/), start=(/1,1,1,1,ind_time/))
        if (fu_fails(stat == NF90_NOERR, 'Failed to read kernel', sub_name)) return
        do iy = 1, ny
          do ix = 1, nx
            call get_vert_kern(disp_lev_bnds(1:nz_dispersion+1), lev_bnds_1d(1:nz_kern+1), &
                             & transpose(kern_values_4d(1:nz_kern,1:nz,ix,iy)), kern_unit, &
                             & kernels(:,:,ix,iy))
            if (error) return
          end do
        end do
        deallocate(kern_values_4d)
        call free_work_array(lev_bnds_1d)

      case ('multi_vertical_var')
        allocate(lev_bnds_3d(nz_kern+1, nx, ny), kern_values_4d(nz_kern, nz, nx, ny), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed', 'set_aver_kern')) return
        stat = nf90_get_var(ncid, var_id_lev_bnds, lev_bnds_3d, count=(/nz_kern+1, nx, ny/), &
                          & start=(/1, 1, 1, ind_time/))
        if (fu_fails(stat == NF90_NOERR, 'Failed to read level boundaries', sub_name)) return
        
        stat = nf90_get_var(ncid, var_id_kernel, kern_values_4d, &
                          & count=(/nz_kern, nz, nx, ny/), start=(/1,1,1,1,ind_time/))
        if (fu_fails(stat == NF90_NOERR, 'Failed to read kernel', sub_name)) return
        do iy = 1, ny
          do ix = 1, nx
            call get_vert_kern(disp_lev_bnds(1:nz_dispersion+1), lev_bnds_3d(1:nz_kern+1,ix,iy), &
                             & transpose(kern_values_4d(1:nz_kern,1:nz,ix,iy)), kern_unit, &
                             & kernels(:,:,ix,iy))
            if (error) return
          end do
        end do
        deallocate(lev_bnds_3d, kern_values_4d)

      case default
        call set_error('Unsupported kernel type', 'set_aver_kern')
      end select

    end subroutine set_aver_kern_disp_vert

    !=====================================================
    
    subroutine inq_var(var_id)
      implicit none
      integer :: var_id
      
      integer, dimension(:), pointer :: int_wrk
      integer :: ind_dim, ndims, dimlen, stat

      int_wrk => fu_work_int_array()
      
      stat = nf90_inquire_variable(ncid, var_id, dimids=int_wrk, ndims=ndims)
      do ind_dim = 1, ndims
        stat = nf90_inquire_dimension(ncid, int_wrk(ind_dim), len=dimlen)
        print *, ind_dim, dimlen
      end do
      
      call free_work_array(int_wrk)
    end subroutine inq_var

  end subroutine set_vert_obs_from_nc

  
  !******************************************************************************************
  
  logical elemental function vert_obs_defined(obs) result(is_defined)
    implicit none
    type(t_vertical_observation), intent(in) :: obs
    is_defined = obs%defined
  end function vert_obs_defined


  
!!$  subroutine get_profile_kern_trivial(dispersion_vertical, nx, ny, nz, kernel)
!!$    implicit none
!!$    type(silam_vertical), intent(in) :: vertical
!!$    integer, intent(in) :: nx, ny, nz ! of observation
!!$    real, dimension(:,:,:,:), pointer :: kernel
!!$
!!$    integer :: nz_disp, ilev, stat, dim_id_kern_layer_intrfc, nz_kern_intrfc
!!$    character(len=*), parameter :: subname = 'get_profile_kernel_trivial'
!!$    character(len=clen) :: kern_unit
!!$    logical :: can_project_kernel
!!$
!!$    stat = nf90_inq_varid(ncid, 'kernel_layer_interfaces', var_id_lev_bnds)
!!$    if (fu_fails(stat == NF90_NOERR, 'Failed to inquire kernel_layer_interfaces', subname)) return
!!$    stat = nf90_inq_dimid(ncid, 'kernel_layer_interfaces', dim_id_kern_layer_intrfc)
!!$    if (fu_fails(stat == NF90_NOERR, 'Failed to inq kernel_layer_interfaces dim', subname)) return
!!$    stat = nf90_inquire_dimension(ncid, dim_id_kern_layer_intrfc, len=nz_kern_intrfc)
!!$    if (fu_fails(stat == NF90_NOERR, 'Failed to inquire kernel_layer dimension size', subname)) return
!!$    if (fu_fails(nz_kern_intrfc == nz+1, 'Kernel layers sizes don''t match', subname)) return
!!$    
!!$    stat = nf90_get_var(ncid, var_id_lev_bnds, lev_bnds_1d, count=(/nz_kern_intrfc/), start=(/1/))
!!$    if (fu_fails(stat == NF90_NOERR, 'Failed to read level boundaries', subname)) return
!!$
!!$    nz_disp = fu_NbrOfLevels(dispersion_vertical)
!!$    allocate(kernel(nz_disp, nz, nx, ny), stat=stat)
!!$    if (fu_fails(stat == 0, 'Allocate failed', subname)) return
!!$
!!$    stat = nf90_get_att(ncid, var_id_lev_bnds, 'units', kern_unit)
!!$    if (fu_fails(stat == NF90_NOERR, 'Failed to get kern_unit', subname)) return
!!$    
!!$    can_project_kernel = kern_unit == 'm' .and. fu_leveltype(dispersion_vertical) == constant_height
!!$
!!$    if (can_project_kernel) then
!!$      ! Projection done immediately
!!$      do iy = 1, ny
!!$        do ix = 1, nx
!!$          kernels(:,:,ix,iy) = kern_values_2d(1:nz_kern, 1:nz)
!!$          kern_lev_bnds(:,ix,iy) = lev_bnds_1d(1:nz_kern+1)
!!$        end do
!!$      end do
!!$    else
!!$      ! Kernel stored in native vertical, projection done later
!!$
!!$    end if
!!$    kernel = 1.0
!!$    
!!$  end subroutine get_profile_kern_trivial


  !************************************************************************************
  !
  ! Test codes
  !
  !************************************************************************************

  subroutine make_test_obs(which, valid_time, transport_species, optical_species, obs)
    implicit none
    character(len=*), intent(in) :: which
    type(silam_species), dimension(:), intent(in) :: transport_species
    type(silam_species), dimension(:), pointer :: optical_species
    type(silja_time), intent(in) :: valid_time
    type(t_vertical_observation), intent(out) :: obs
    
    real, dimension(:,:,:), allocatable :: val_3d, std_3d
    real, dimension(:), allocatable:: lon_2d, lat_2d
    real, dimension(:,:,:,:), pointer :: vert_kern
    real, dimension(:), pointer :: chem_kern, disp_vert_bnds_m, kern_vert_bnds_m, kern1d
    real, dimension(:,:), pointer :: kern_values, kern, work2d
    real, dimension(:,:,:,:), allocatable :: kern2d
    real, dimension(0,0,0) :: kern_lev_bnds ! won't use
    type(silam_species), dimension(:), pointer :: observed_species, dummy_species_ptr
    integer nx, ny, num_obs, ilev, nzdisp
    type(t_vertical_observation), dimension(:), pointer :: obs_list
    
    select case(which)
    case ('column_so2')
      nx = 2
      ny = 2
      call get_vert_kern_trivial(dispersion_vertical, nx, ny, vert_kern)
      allocate(observed_species(1))
      call set_species(observed_species(1), fu_get_material_ptr('SO2'), in_gas_phase) 
      call get_chem_kern(observed_species, transport_species, 'mole/m2', chem_kern, if_profile=.false.)
      !allocate(val_3d(1,nx,ny), std_3d(1,nx,ny), lon_2d(nx*ny), lat_2d(nx*ny))
      allocate(val_3d(nx,ny,1), std_3d(nx,ny,1), lon_2d(nx*ny), lat_2d(nx*ny))

      val_3d(1,1,1) = 1.0
      val_3d(2,1,1) = 2.0
      val_3d(1,2,1) = 0.0
      val_3d(2,2,1) = real_missing
      
      val_3d = 0.0

      std_3d = 1.0
      
      lon_2d(1) = 10.0 ! 1, 1
      lon_2d(2) = 10.5 ! 2, 1
      lon_2d(3) = 10.0
      lon_2d(4) = 10.5
      lat_2d = (/40.0, 40.0, 40.5, 40.5/)
      
      call create_vert_obs(val_3d, std_3d, lon_2d, lat_2d, nx, ny, chem_kern, vert_kern, .false., &
                         & kern_lev_bnds, &
                         & valid_time, which, '', .false., .false., dispersion_grid, linear, obs)

    case ('column_so2_kern')
      nx = 2
      ny = 2
      
      kern => fu_work_array_2d()
      kern1d => fu_work_array()

      disp_vert_bnds_m => fu_work_array()
      nzdisp = fu_NbrOfLevels(dispersion_vertical)
      do ilev = 1, nzdisp
        disp_vert_bnds_m(ilev+1) = &
             & fu_level_height(fu_upper_boundary_of_layer(fu_level(dispersion_vertical, ilev)))
        disp_vert_bnds_m(ilev) = &
             & fu_level_height(fu_lower_boundary_of_layer(fu_level(dispersion_vertical, ilev)))
      end do
      kern_vert_bnds_m => fu_work_array()
      kern_vert_bnds_m(1:5) = (/0.0, 500.0, 1000.0, 5000.0, 100000.0/)
      kern_values => fu_work_array_2d()
      !kern_values(1,1:4) = (/0.1, 0.2, 0.3, 0.4/)
      kern_values(1,1:4) = (/1.0, 1.0, 1.0, 1.0/)
      work2d => fu_work_array_2d()
      kern => work2d(1:1, 1:nzdisp)
      call get_vert_kern(disp_vert_bnds_m(1:nzdisp+1), kern_vert_bnds_m(1:5), kern_values(1:1,1:4), 'm', kern)
      call free_work_array(kern_vert_bnds_m)
      call free_work_array(kern_values)
      call free_work_array(disp_vert_bnds_m)

      allocate(observed_species(1))
      call set_species(observed_species(1), fu_get_material_ptr('SO2'), in_gas_phase)
      call get_chem_kern(observed_species, transport_species, 'mg/m2', chem_kern, .false.)
      !allocate(val_3d(1,nx,ny), std_3d(1,nx,ny), lon_2d(nx*ny), lat_2d(nx*ny))
      allocate(val_3d(nx,ny,1), std_3d(nx,ny,1), lon_2d(nx*ny), lat_2d(nx*ny))

      val_3d(1,1,1) = 1.0
      val_3d(2,1,1) = 2.0
      val_3d(1,2,1) = 0.0
      val_3d(2,2,1) = real_missing
      
      std_3d = 0.5
      
      lon_2d(1) = 10.0 ! 1, 1
      lon_2d(2) = 10.5 ! 2, 1
      lon_2d(3) = 10.0
      lon_2d(4) = 10.5
      lat_2d = (/40.0, 40.0, 40.5, 40.5/)
      
      allocate(kern2d(nzdisp, 1, nx,ny))
      kern2d(:,:,1,1) = transpose(kern)
      kern2d(:,:,2,1) = transpose(kern)
      kern2d(:,:,1,2) = transpose(kern)
      kern2d(:,:,2,2) = transpose(kern)

      call create_vert_obs(val_3d, std_3d, lon_2d, lat_2d, nx, ny, chem_kern, kern2d, .false., &
                         & kern_lev_bnds, valid_time, which, '', .false., .false., dispersion_grid, linear, obs)
      

    case ('file')
      allocate(observed_species(1), obs_list(100))
      call set_species(observed_species(1), fu_get_material_ptr('SO2'), in_gas_phase)
      nullify(dummy_species_ptr)
      call set_vert_obs_from_nc('test.nc', 'col_dens', observed_species, &
                              & transport_species, dummy_species_ptr, &
                              & really_far_in_past, really_far_in_future, obs_list, num_obs)
      if (error) return
      if (fu_fails(num_obs > 0, 'Failed to read observations', 'make_test_obs')) return
      obs = obs_list(1)

    case ('file_kern')
      allocate(observed_species(1), obs_list(100))
      call set_species(observed_species(1), fu_get_material_ptr('SO2'), in_gas_phase)
      nullify(dummy_species_ptr)
      call set_vert_obs_from_nc('testdata_kern.nc', 'SO2', &
                              & observed_species, transport_species, dummy_species_ptr, &
                              & really_far_in_past, really_far_in_future, obs_list, num_obs)
      if (error) return
      if (fu_fails(num_obs > 0, 'Failed to read observations', 'make_test_obs')) return
      obs = obs_list(1)

    case ('file_omi_kern')
      allocate(observed_species(1), obs_list(100))
      call set_species(observed_species(1), fu_get_material_ptr('NO2'), in_gas_phase)
      nullify(dummy_species_ptr)
      call set_vert_obs_from_nc('testdata_omi_kern.nc', 'NO2', &
                              & observed_species, transport_species, dummy_species_ptr, &
                              & really_far_in_past, really_far_in_future, obs_list, num_obs)
      if (error) return
      if (fu_fails(num_obs > 0, 'Failed to read observations', 'make_test_obs')) return
      obs = obs_list(1)

    case ('file_aod')
      allocate(observed_species(1))
      call set_species(observed_species(1), fu_get_material_ptr('SO2'), in_gas_phase) 
      call set_vert_obs_from_nc('test_aod.nc', 'aot', observed_species, &
                              & transport_species, optical_species, &
                              & really_far_in_past, really_far_in_future, obs_list, num_obs)
      if (error) return
      if (fu_fails(num_obs > 0, 'Failed to read observations', 'make_test_obs')) return
      obs = obs_list(1)

      
    case default 
      call set_error('No such option: ' // trim(which), 'make_test_obs')
    end select

  end subroutine make_test_obs

  
  !******************************************************************************************
  
  subroutine test_vert_obs(map_c, map_px, map_py, map_pz, now, timestep, &
                         & p_rel_hum, p_tempr, p_press, p_hgt, p_surf_press, &
                          & met_disp_interp_horiz, if_met_disp_interp_horiz, &
                          & met_disp_interp_vert, if_met_disp_interp_vert, &
                          & weight_past, rules_opt_dens, optical_species)
    implicit none
    !type(t_vertical_observation), intent(inout) :: obs
    type(Tmass_map) :: map_c, map_px, map_py, map_pz
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep
    type(field_4d_data_ptr), pointer :: p_rel_hum, p_tempr, p_press, p_hgt
    type(field_2d_data_ptr), pointer :: p_surf_press
    type(THorizInterpStruct), pointer :: met_disp_interp_horiz
    logical, intent(in) :: if_met_disp_interp_horiz
    type(TVertInterpStruct), pointer :: met_disp_interp_vert
    logical, intent(in) :: if_met_disp_interp_vert
    real, intent(in) :: weight_past
    type(Toptical_density_rules), intent(in) :: rules_opt_dens
    type(silam_species), dimension(:), pointer :: optical_species
    type(t_vertical_observation) :: obs
    real, dimension(:), allocatable :: val_obs, val_mdl, var
    integer :: n_obs, ilev
    real :: right, left, area, dz
    real, dimension(:,:,:,:,:), allocatable :: x_rhs
    real, dimension(:,:,:,:,:), pointer :: map_c_in 
    integer :: ix, iy, n_opt_species
    
    if (associated(optical_species)) then
      n_opt_species = size(optical_species) 
    else
      n_opt_species = 0
    end if
    
    !call make_test_obs('column_so2', now+timestep/2.0, map_c%species, obs)
    !call make_test_obs('column_so2_kern', now+timestep/2.0, map_c%species, optical_species, obs)
    !call make_test_obs('file_kern', now+timestep/2.0, map_c%species, optical_species, obs)
    call make_test_obs('file_omi_kern', now+timestep/2.0, map_c%species, optical_species, obs)
    !call make_test_obs('file', now+timestep/2.0, map_c%species, optical_species, obs)
    !call make_test_obs('file_aod', now+timestep/2.0, map_c%species, optical_species, obs)
    if (error) return
    
    map_c%arm = 0.0
    
    n_obs = fu_size(obs)
    allocate(val_obs(n_obs), val_mdl(n_obs), var(n_obs), &
           & x_rhs(map_c%nspecies, map_c%nsrc, map_c%n3d, map_c%nx, map_c%ny))
    
    call random_number(x_rhs)
    map_c%arm(1:map_c%nspecies, 1:map_c%nsrc, 1:map_c%n3d, 1:map_c%nx, 1:map_c%ny) = x_rhs    
    call random_number(val_obs)

    ! check the identity < Rinv.(y - Hx), Hx > = < H*.Rinv.(y - Hx), x > 
    
    ! note: the control variable x in concentration, but when given to observe_vertical,
    ! map_c is assumed to have cell mass. In the final computation, the values in x are
    ! hence divided by volume to evaluate the expression on RHS.
    
    !obs%disp_area = 1.0
    call observe_vertical(obs, map_c, &
                        & now, timestep, p_rel_hum, p_tempr, p_press, p_hgt, p_surf_press, &
                        & met_disp_interp_horiz, if_met_disp_interp_horiz, &
                        & met_disp_interp_vert, if_met_disp_interp_vert, &
                        & weight_past, rules_opt_dens, n_opt_species)
    if (error) return
    call obs_to_file(obs, optical_species, 6)
    call get_data(obs, values_obs = val_obs)
    call get_data(obs, values_mdl = val_mdl)
    call get_data(obs, variance = var)

    ! now: 
    ! map_c%arm = x
    ! val_mdl = Hx
    ! val_obs = y
    ! 1/var = Rinv
    
    left = sum( (val_mdl - val_obs)/var * val_mdl  ) ! < Rinv . (y - Hx), Hx >

    map_c%arm = 0.0
    call inject_vertical(obs, map_c, map_px, map_py, map_pz, &
                       & now, timestep, p_rel_hum, p_tempr, &
                       & met_disp_interp_horiz, if_met_disp_interp_horiz, &
                       & met_disp_interp_vert, if_met_disp_interp_vert, &
                       & weight_past, rules_opt_dens, n_opt_species)
    ! now map_c%arm = H* Rinv.(y - Hx)

    ! after injection we need to involve volume:
    right = 0.0
    do ilev = 1, map_c%n3d
      dz = fu_layer_thickness_m(fu_level(map_c%vertTemplate, ilev))
      do iy = 1, map_c%ny
        do ix = 1, map_c%nx
          area = fu_cell_size(dispersion_grid, ix, iy)
          right = right + sum(map_c%arm(:,1,ilev,ix,iy) * x_rhs(:,1,ilev,ix,iy) / (area*dz))
        end do
      end do
    end do
    !right = sum( map_c%arm(1:map_c%nspecies, 1:map_c%nsrc, 1:map_c%n3d, 1:map_c%nx, 1:map_c%ny) *  x_rhs )
    call msg('LEFT', left)
    call msg('RIGHT', right)

    ! hack the cell area and chem kernel to 1 so we should get exactly the kernel values.
    obs%disp_area = 1.0
    where (obs%chem_kernel > 0) 
      obs%chem_kernel = 1.0
    end where
    do ilev = 1, map_c%n3d
      map_c%arm = 0.0
      map_c%arm(:,1,ilev,:,:) = 1.0
      call msg('level:', ilev)
      call msg('range', fu_level_height(fu_lower_boundary_of_layer(fu_level(dispersion_vertical, ilev))), &
             & fu_level_height(fu_upper_boundary_of_layer(fu_level(dispersion_vertical, ilev))))
      call observe_vertical(obs, map_c, &
                          & now, timestep, p_rel_hum, p_tempr, p_press, p_hgt, p_surf_press, &
                          & met_disp_interp_horiz, if_met_disp_interp_horiz, &
                          & met_disp_interp_vert, if_met_disp_interp_vert, &
                          & weight_past, rules_opt_dens, n_opt_species)
      call get_data(obs, values_mdl=val_mdl)
      print *, 'val_mdl:', val_mdl
    end do

  end subroutine test_vert_obs

end module observations_vertical
