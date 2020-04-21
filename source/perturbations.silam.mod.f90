!
! Module for generating perturbations. 
! 
! The use sequence is as follows:
! 1. set DA rules to define the perturbed variables
! 2. set_perturbations: reads the statistical configuration and initializes the correlation operators
! 3. apply_perturbations: call get_perturbation and sends it to the appropriate place in model.
! 
! Note that some perturbations might also come from giving different input data for each
! ensemble member.

module perturbations
  use io_server
  use da_common
  use source_apportionment

  implicit none
  
!  private

  interface defined
    module procedure pert_defined
  end interface

  public set_perturbations
  public get_perturbation
  public defined
  public apply_perturbations

  public perturb_mass_map
  public perturb_emis_proc_volc
  public perturb_emission_xy
  public store_perturbation_val
  public retrieve_perturbation_val
  public fu_eruption_height

  integer, parameter, private :: distr_normal = 451700, &
                               & distr_normal_squared = 451701, &
                               & distr_normal_posit = 451702, &
                               & distr_lognormal = 451703
  
  type t_perturbation
    private
     integer :: perturb_type = int_missing ! concentration, etc.
     integer :: distribution = distr_normal
     real, dimension(:), pointer :: scale  ! "stdev"
     real, dimension(:), pointer :: offset ! "mean"
     real, dimension(:), pointer :: storage ! auxiliary data for some perturbs
     real, dimension(:,:), pointer :: L ! lower triangular matrix, cov = LL^T
     type(t_correlation) :: correlation
     ! correlation time - with AR(1) method, later
     real :: correlation_time = real_missing
     integer, dimension(:), pointer :: address
     integer :: pert_size ! of scale and offset
     integer :: n_gridpoints
     logical :: have_cov = .false.
     logical :: first_time = .true.
     logical :: defined = .false.
  end type t_perturbation
  public t_perturbation
  !
  ! Types of perturbations known for this module
  !
  integer, parameter, public :: perturb_concentration = 441700, &
!                              & perturb_emis = 441701, &
                              & perturb_emis_volc = 441702, &
                              & perturb_emis_zt = 441703, &
                              & perturb_emis_xy = 441704

CONTAINS

  logical function pert_defined(pert)
    implicit none
    type(t_perturbation), intent(in) :: pert
    pert_defined = pert%defined
  end function pert_defined

  !************************************************************************************

  subroutine set_perturbations(perturbations, species_emission, species_transport, &
                             & rules, analysis_grid, analysis_vertical, analysis_time)
    implicit none
    ! note intent inout: perturbations(:)%perturb_type sets the request.
    type(t_perturbation), dimension(:), intent(out) :: perturbations
    type(silam_species), dimension(:), intent(in), target :: species_emission, species_transport
    type(da_rules), intent(in) :: rules
    type(silja_grid), intent(in) :: analysis_grid
    type(silam_vertical), intent(in) :: analysis_vertical
    type(silja_time), intent(in) :: analysis_time
    
    type(t_background_covariance) :: bgr_cov
    type(da_control) :: backgr_invalid
    type(silam_species), dimension(:), pointer :: p_species_emission, p_species_transport
    character(len=*), parameter :: sub_name = 'set_perturbations'
    integer :: ind_pert, num_contr_species, pert_size, stat, num_pert, n
    logical :: need_cov

    call msg('Setting perturbations')

    ! Get the perturbation types
    num_pert = 0
    need_cov = .false.
    if (fu_have_initial(rules%perturbVariable)) then
      num_pert = num_pert + 1
      perturbations(num_pert)%perturb_type= perturb_concentration
      need_cov = .true.
    end if
    if (fu_have_emission_xy(rules%perturbVariable)) then
      num_pert = num_pert + 1
      perturbations(num_pert)%perturb_type= perturb_emis_xy
    end if
    if (fu_have_emission_zt(rules%perturbVariable)) then
      num_pert = num_pert + 1
      perturbations(num_pert)%perturb_type= perturb_emis_zt
      need_cov = .true.
    end if
    if (fu_have_emission_volc(rules%perturbVariable)) then
      num_pert = num_pert + 1
      perturbations(num_pert)%perturb_type = perturb_emis_volc
    end if
    
    if (num_pert == 0) then
      call msg('No forced perturbations requested')
    end if

    p_species_emission => species_emission
    p_species_transport => species_transport
    
    if (need_cov) then
      ! Make a dummy "background". The stdev_methods which use background value should not
      ! be used for perturbations. Todo: way to check this.
      call init_control_from_par(backgr_invalid, analysis_grid, analysis_vertical, &
                      & p_species_transport, p_species_emission, &
                      & rules, physical_space, initial_value=0.0, controlvar=rules%perturbVariable)
      if (error) return
      
      ! Use a background covariance. This can be conveniently read from a config file. Later
      ! the content will be wrapped into a t_perturbation though.
      call set_cov_mdl(bgr_cov, backgr_invalid, species_emission, species_transport,  rules, &
                     & analysis_grid, analysis_vertical, analysis_time, rules%perturbVariable)
      if (error) return
    end if

    do ind_pert = 1, size(perturbations)
      
      select case(perturbations(ind_pert)%perturb_type)
        case (perturb_concentration)
          call msg('Setting perturbations: concentration')
          perturbations(ind_pert)%scale => bgr_cov%stdev_initial
          pert_size = size(bgr_cov%stdev_initial)
          perturbations(ind_pert)%pert_size = pert_size
          allocate(perturbations(ind_pert)%offset(pert_size), stat=stat)
          if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
          perturbations(ind_pert)%offset = 0.0
          perturbations(ind_pert)%correlation = bgr_cov%correlation_initial

          num_contr_species = fu_nsp_initial(backgr_invalid)
          if (fu_fails( num_contr_species > 0, 'No init species', sub_name)) return

          allocate(perturbations(ind_pert)%address(num_contr_species), stat=stat)
          if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return

          perturbations(ind_pert)%address = fu_isp_initial(backgr_invalid)
          perturbations(ind_pert)%defined = .true.

       case (perturb_emis_xy)
          call msg('Setting perturbations: emission_xy')
          nullify(perturbations(ind_pert)%address)
          num_contr_species = size(species_emission)
          n = fu_number_of_gridpoints(analysis_grid)
          pert_size = n*num_contr_species
          perturbations(ind_pert)%pert_size = pert_size
          perturbations(ind_pert)%n_gridpoints = n
          allocate(perturbations(ind_pert)%scale(pert_size), perturbations(ind_pert)%offset(pert_size), &
               & perturbations(ind_pert)%storage(pert_size), stat=stat)
          if (fu_fails(stat == 0, 'Allocate failed 0', sub_name)) return
          perturbations(ind_pert)%offset(1:pert_size) = 0.0
          perturbations(ind_pert)%scale(1:pert_size) = rules%emission_stdev
          perturbations(ind_pert)%correlation_time = rules%emis_corr_time

          if (fu_fails( num_contr_species > 0, 'No emis species', sub_name)) return

          allocate(perturbations(ind_pert)%address(num_contr_species), stat=stat)
          if (fu_fails(stat == 0, 'Allocate failed 1', sub_name)) return
          perturbations(ind_pert)%distribution = distr_lognormal

          if (rules%emis_corr_dist > 0) then
            perturbations(ind_pert)%have_cov = .true.
            allocate(perturbations(ind_pert)%L(n, n), stat=stat)
            if (fu_fails(stat == 0, 'Allocate failed 2', sub_name)) return
            call sqrt_cov(perturbations(ind_pert), rules%emission_stdev, rules%emis_corr_dist, analysis_grid, rules%emis_max_corr)
          end if
          perturbations(ind_pert)%defined = .true.


!!$      case (perturb_emis)
!!$        call msg('Setting perturbations: emission')
!!$        perturbations(ind_pert)%scale => bgr_cov%stdev_emission
!!$        pert_size = size(bgr_cov%stdev_emission)
!!$        perturbations(ind_pert)%pert_size = pert_size
!!$        allocate(perturbations(ind_pert)%offset(pert_size), stat=stat)
!!$        if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
!!$        perturbations(ind_pert)%offset = 0.0
!!$        perturbations(ind_pert)%correlation = bgr_cov%correlation_emission
!!$        if (fu_fails(associated(backgr_invalid%ind_species_emis), 'No emis species', sub_name)) return
!!$        num_contr_species = size(backgr_invalid%ind_species_emis)
!!$        allocate(perturbations(ind_pert)%address(num_contr_species), stat=stat)
!!$        if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
!!$        perturbations(ind_pert)%address = backgr_invalid%ind_species_emis
!!$        perturbations(ind_pert)%defined = .true.        
!!$        perturbations(ind_pert)%distribution = distr_normal_squared

        case (perturb_emis_zt)
          call msg('Setting perturbations: emission_zt')
          perturbations(ind_pert)%scale => bgr_cov%stdev_emission
          call msg('scale:', perturbations(ind_pert)%scale)
          pert_size = size(bgr_cov%stdev_emission)
          perturbations(ind_pert)%pert_size = pert_size
          allocate(perturbations(ind_pert)%offset(pert_size), stat=stat)
          if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
          perturbations(ind_pert)%offset = 1.0
          perturbations(ind_pert)%correlation = bgr_cov%correlation_emission
      
          num_contr_species = fu_nsp_emis_ctrl(backgr_invalid)
          if (fu_fails( num_contr_species > 0, 'No emis species', sub_name)) return
          allocate(perturbations(ind_pert)%address(num_contr_species), perturbations(ind_pert)%storage(1), &
                 & stat=stat)
          if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
          perturbations(ind_pert)%address = fu_isp_emis_ctrl(backgr_invalid)
          perturbations(ind_pert)%defined = .true.        
          perturbations(ind_pert)%distribution = distr_lognormal
        

!!$      case (perturb_emis_volc)
!!$        call msg('Setting perturbations: volcanic emission')
!!$        pert_size = size(species_emission)
!!$        allocate(perturbations(ind_pert)%scale(pert_size), perturbations(ind_pert)%offset(pert_size), &
!!$               & perturbations(ind_pert)%storage(2), stat=stat)
!!$        if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
!!$        nullify(perturbations(ind_pert)%address)
!!$        perturbations(ind_pert)%pert_size = pert_size
!!$        perturbations(ind_pert)%scale = rules%volc_massflux_scaling
!!$        perturbations(ind_pert)%offset = rules%volc_massflux_sigma
!!$        perturbations(ind_pert)%storage(1) = real_missing ! will be the previous height
!!$        perturbations(ind_pert)%storage(2) = rules%volc_hgt_corr_time
!!$        perturbations(ind_pert)%defined = .true.        

      case (perturb_emis_volc)
        call msg('Setting perturbations: volcanic emission')
        pert_size = size(species_emission)
        allocate(perturbations(ind_pert)%scale(pert_size), perturbations(ind_pert)%offset(pert_size), &
               & perturbations(ind_pert)%storage(7), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
        nullify(perturbations(ind_pert)%address)
        perturbations(ind_pert)%pert_size = pert_size
        perturbations(ind_pert)%scale = rules%volc_massflux_scaling
        perturbations(ind_pert)%offset(1) = rules%volc_volflux_mode
        perturbations(ind_pert)%storage(1) = real_missing ! will be logarithm of previous flux
        perturbations(ind_pert)%storage(2) = rules%volc_volflux_corr_time
        perturbations(ind_pert)%storage(3) = rules%volc_volflux_sigma
        perturbations(ind_pert)%storage(4) = rules%volc_height_exp
        perturbations(ind_pert)%storage(5) = rules%volc_height_scaling
        perturbations(ind_pert)%storage(6) = rules%volc_fract_mass_top
        perturbations(ind_pert)%storage(7) = rules%volc_fract_height_top
        perturbations(ind_pert)%defined = .true.        
        
        case (int_missing)
          exit

        case default
          call set_error('Strange perturbation', sub_name)
      end select
    end do

    ! take care. The "background" is not needed so it can be deallocated. But the data in
    ! bgr_cov is now used via the perturbations, so it must not be deallocated.
    if (defined(backgr_invalid)) call destroy(backgr_invalid)
    
    ! get perturbation items

    ! allocate perturbation rules list

    ! get the parameters from the separate file(s)
    
  contains

    subroutine init_perturb(perturbation, total_size)
      implicit none
      integer :: total_size
      type(t_perturbation), intent(inout) :: perturbation
      
    end subroutine init_perturb
    
  end subroutine set_perturbations

  
  !************************************************************************************

  subroutine apply_perturbations(cloud, disp_buffer, perturbations, timestep, wdr, OutDef)
    !
    ! Applies the required perturbations
    !
    implicit none
    ! Imported parameters
    type(silam_pollution_cloud), intent(inout) :: cloud
    type(TField_buffer), intent(inout) :: disp_buffer
    type(t_perturbation), dimension(:), intent(inout) :: perturbations
    type(silja_interval), intent(in) :: timestep
    type(silja_wdr), intent(inout) :: wdr
    type(silam_output_definition), intent(inout) :: OutDef

    ! Local variables
    integer :: ind_pert
    type(Tmass_map), pointer :: p_cnc_mm
    
    call msg('Apply perturbations')

    do ind_pert = 1, size(perturbations)
      select case(perturbations(ind_pert)%perturb_type)

        case (perturb_concentration)
          call msg('...concentration')
          p_cnc_mm => fu_concMM_ptr(cloud)
          call perturb_mass_map(p_cnc_mm, perturbations(ind_pert))
          
        case (perturb_emis_xy)
          call msg('...emission map')
          if(.not. associated(fu_emission_processor_ptr(cloud))) then
            call set_error('emission processor not associated', 'apply_perturbations')
            return
          end if
          call perturb_emission_xy(fu_emission_processor_ptr(cloud), perturbations(ind_pert), timestep)

        case (perturb_emis_zt)
          call msg('...z-t emission')
          if (.not. associated(fu_emission_processor_ptr(cloud))) then
            call set_error('emission processor not associated', 'apply_perturbations')
            return
          end if
          call perturb_emis_proc_zt(fu_emission_processor_ptr(cloud), perturbations(ind_pert), timestep)

        case (perturb_emis_volc)
          call msg('...parametric emission for volcanoes')
          if (fu_fails(associated(fu_emission_processor_ptr(cloud)), &
                       & 'emission processor not associated', 'apply_perturbations'))return
          !call perturb_emis_proc_zt_param(fu_emission_processor_ptr(cloud), perturbations(ind_pert), timestep)
          call perturb_emis_proc_volc(fu_emission_processor_ptr(cloud), perturbations(ind_pert), timestep)
          
        case (int_missing)
          exit
        case default
          call set_error('Bad perturbation', 'apply_perturbations')
      end select
      if(error)return  
    end do  ! perturbations

  end subroutine apply_perturbations
  
  
  !************************************************************************************
  
  subroutine perturb_mass_map(mass_map, perturbation)
    !
    ! Perturbs the concentration: add noise to mass in the mass map
    !
    implicit none
    type(TMass_map), intent(inout) :: mass_map
    type(t_perturbation), intent(inout) :: perturbation
    
    integer :: size_total
    real :: area, volume
    real, dimension(:), pointer :: thickness, values
    integer :: iy, ix, iz, isp, i1d, isp_transp, num_species
    integer :: nx, ny, nz

    num_species = size(perturbation%address) 
    if (fu_fails(num_species <= mass_map%nSpecies, 'Size mismatch', 'perturb_mass_map')) return

    nx = mass_map%nx
    ny = mass_map%ny
    nz = mass_map%n3d
    size_total = nx * ny * nz * num_species

    values => fu_work_array(size_total)
    thickness => fu_work_array()
    thickness(1:mass_map%n3d) = (/(fu_layer_thickness_m(fu_level(mass_map%vertTemplate, iz)),&
                                  & iz=1, mass_map%n3d)/)

    call get_perturbation(perturbation, values)

    do iy = 1, mass_map%ny
      do ix = 1, mass_map%nx
        area = fu_cell_size(mass_map%gridTemplate, ix, iy)
        do iz = 1, mass_map%n3d
          volume = thickness(iz) * area
        
          do isp = 1, num_species
            isp_transp = perturbation%address(isp)
            i1d = (iy-1)*nx*nz*num_species + (ix-1)*nz*num_species + (iz-1)*num_species + isp
            mass_map%arm(isp_transp,1,iz,ix,iy) = mass_map%arm(isp_transp,1,iz,ix,iy) + values(i1d)*volume
            if (mass_map%arm(isp_transp,1,iz,ix,iy) < 0) mass_map%arm(isp_transp,1,iz,ix,iy) = 0.0
          end do

        end do
      end do
    end do
    
    call free_work_array(values)
    call free_work_array(thickness)

  end subroutine perturb_mass_map

  
  !************************************************************************************

  subroutine perturb_emis_proc_zt(proc, pert, timestep)
    implicit none
    type(Temission_processor) :: proc
    type(t_perturbation), intent(inout) :: pert
    type(silja_interval), intent(in) :: timestep

    real, dimension(:), pointer :: values
    character(len=*), parameter :: subname = 'perturb_emis_proc_zt'

    real, dimension(:,:), allocatable :: sample
    integer :: ii, ind_sample
    real, parameter :: std_logn = sqrt((exp(1.0) - 1.0)*exp(1.0)), exp_05 = exp(0.5)
    real, parameter :: std_req = 1.0
    integer, parameter :: sample_size = 10000

    if (.not. pert%first_time) return

    values => fu_data(proc)
    if (fu_fails(associated(values), 'values not associated', subname)) return
    
    call msg('...first time. Perturbing emission processor z-t')
    call get_perturbation(pert, values)
    call msg('...min, max of perturbed values', minval(values), maxval(values))
    call msg('values', values)
!!$
!!$    call msg('will generate 10000 samples')
!!$    allocate(sample(size(values), sample_size))
!!$    do ii = 1, sample_size
!!$      call get_perturbation(pert, values)
!!$      sample(:,ii) = values
!!$    end do
!!$
!!$    call msg('Mean:', sum(sample, 2) / sample_size)
!!$
!!$    call msg('Now the same in-line')
!!$
!!$     do ind_sample = 1, sample_size
!!$      call random_normal(sample(:,ind_sample))
!!$      call transf_fwd(pert%correlation, sample(1:fu_dimension_control(pert%correlation), ind_sample), values)
!!$      sample(:,ind_sample) = values
!!$      values = exp(sample(:, ind_sample))
!!$      sample(:,ind_sample) = values * std_req / std_logn
!!$      values = exp_05 * std_req/std_logn
!!$      sample(:,ind_sample) = sample(:,ind_sample) - values + 1
!!$    end do
!!$    call msg('Mean:', sum(sample, 2) / sample_size)

    pert%first_time = .false.

  end subroutine perturb_emis_proc_zt
 
  
  !************************************************************************************                                                                                                                                                                                                                                    

  subroutine perturb_emission_xy(proc, pert, timestep)
    implicit none
    type(Temission_processor) :: proc
    type(t_perturbation), intent(inout) :: pert
    type(silja_interval), intent(in) :: timestep

    real, dimension(:), pointer :: values_pert, proc_data, values_pert_sp, rand_vect
    type(silam_species), dimension(:), pointer :: species_list

    character(len=*), parameter :: sub_name = 'perturb_emis_xy'

    real :: alpha, corr_time
    integer :: nn, stat, n_gridpoints
    integer :: ix, iy, nx, ny, ind_species, nspecies, i
    integer :: ind_PAR = -999, ind_OLE = -999, ind_P, ind_O
    integer :: ind_SO2 = -999, ind_SO4 = -999, ind_S2, ind_S4
    character(len=substNmLen) :: material_name
    real, parameter :: fact = 0.5*log(10.)**2

    ! The solution for perturbing PAR, OLE, SO2, and SO4 is intended to be temporary only                                                                                                                                                                                                                                  

    call msg('perturbing emis xy', pert%pert_size)
    call msg('storage size', size(pert%storage))

    nn = pert%pert_size
    n_gridpoints = pert%n_gridpoints
    call grid_dimensions(fu_grid(proc), nx, ny)

    species_list => fu_species_list(proc)
    nspecies = size(species_list)

    call msg('average storage before perturbation', sum(pert%storage(1:nn))/nn)
    call msg('max storage before perturbation', maxval(pert%storage(1:nn)))

    proc_data => fu_data(proc)
    if (fu_fails(associated(proc_data), 'proc_data not associated', sub_name)) return

    values_pert => fu_work_array(nn)
    rand_vect => fu_work_array(nn)

    call random_normal(rand_vect)

    if (pert%have_cov) then
      !$OMP PARALLEL default(none) private(ind_species, values_pert_sp, ix, iy) &                                                                                                                                                                                                                                          
      !$OMP & shared(values_pert, pert, nspecies, n_gridpoints, nx, ny, rand_vect)                                                                                                                                                                                                                                         
      values_pert_sp => fu_work_array(n_gridpoints)
      !OMP DO                                                                                                                                                                                                                                                                                                              
      do ind_species = 1, nspecies
        values_pert_sp(1:n_gridpoints) =  rand_vect((ind_species-1)*n_gridpoints+1:ind_species*n_gridpoints)
        call strmv('L','N', 'N', n_gridpoints, pert%L, n_gridpoints, values_pert_sp(1:n_gridpoints), 1)
        do ix = 1, nx
          do iy = 1, ny
            values_pert((iy-1)*nspecies*nx + (ix-1)*nspecies + ind_species) = values_pert_sp((iy-1)*nx+ix)
          end do
        end do
      end do
      !OMP END DO                                                                                                                                                                                                                                                                                                          
      call free_work_array(values_pert_sp)
      !$OMP END PARALLEL                                                                                                                                                                                                                                                                                                   
    end if

    call free_work_array(rand_vect)

    if (pert%first_time) then
      if (.not. pert%have_cov) then
        call random_normal(values_pert(1:nn))
        values_pert(1:nn) = values_pert(1:nn) * pert%scale(1:nn)
      end if
      pert%storage(1:nn) = values_pert(1:nn)
      pert%first_time = .false.
    else
      alpha = 1 - fu_sec(timestep) / pert%correlation_time

      if (.not. pert%have_cov) then
        call random_normal(values_pert(1:nn))
        values_pert(1:nn) = values_pert(1:nn) * pert%scale(1:nn)
      end if
      pert%storage(1:nn) = alpha*pert%storage(1:nn) + sqrt(1-alpha**2) * values_pert(1:nn)
    end if

    if (pert%distribution == distr_lognormal) then
      proc_data(1:nn) = 10**(pert%storage(1:nn)) / exp(fact*pert%scale(1:nn)**2)
    else if (pert%distribution == distr_normal) then
      !$OMP PARALLEL DO default(none) private(i) shared(pert, nn)                                                                                                                                                                                                                                                          
      do i = 1, nn
        if (pert%storage(i) < -1.) then
          pert%storage(i) = -1.
        end if
      end do
      !$OMP END PARALLEL DO                                                                                                                                                                                                                                                                                                
      proc_data(1:nn) = 1. + pert%storage(1:nn)
    end if

    nspecies = size(species_list)

    do ind_species = 1, nspecies
      material_name = fu_substance_name(species_list(ind_species))
      if (material_name == 'PAR') then
        ind_PAR = ind_species
      else if (material_name == 'OLE') then
        ind_OLE = ind_species
      else if (material_name == 'SO2') then
        ind_SO2 = ind_species
      else if (material_name == 'SO4') then
        ind_SO4 = ind_species
      end if
    end do

    !$OMP PARALLEL default(none) private(ix, iy, ind_P, ind_O, ind_S2, ind_S4) &                                                                                                                                                                                                                                           
    !$OMP & shared(proc_data, nx, ny, nspecies, ind_PAR, ind_OLE, ind_SO2, ind_SO4)                                                                                                                                                                                                                                        
    if (ind_PAR > 0 .and. ind_OLE > 0) then
      !$OMP DO                                                                                                                                                                                                                                                                                                             
      do ix = 1, nx
        do iy = 1, ny
          ind_P = (iy-1)*nspecies*nx + (ix-1)*nspecies + ind_PAR
          ind_O = (iy-1)*nspecies*nx + (ix-1)*nspecies + ind_OLE
          proc_data(ind_P) = proc_data(ind_O)
        end do
      end do
      !$OMP END DO                                                                                                                                                                                                                                                                                                         
    end if
    if (ind_SO2 > 0 .and. ind_SO4 > 0) then
      !$OMP DO                                                                                                                                                                                                                                                                                                              
      do ix = 1, nx
        do iy = 1, ny
          ind_S2 = (iy-1)*nspecies*nx + (ix-1)*nspecies + ind_SO2
          ind_S4 = (iy-1)*nspecies*nx + (ix-1)*nspecies + ind_SO4
          proc_data(ind_S4) = proc_data(ind_S2)
        end do
      end do
      !$OMP END DO                                                                                                                                                                                                                                                                                                         
    end if
    !$OMP END PARALLEL                                                                                                                                                                                                                                                                                                      
    call free_work_array(values_pert)

    call msg('...min, max of perturbed values', minval(proc_data), maxval(proc_data))
    call msg('...average emission scaling factor', sum(proc_data(1:nn)) / nn)

  end subroutine perturb_emission_xy


  !************************************************************************************

  subroutine perturb_emis_proc_volc(proc, pert, timestep)
    implicit none
    type(Temission_processor) :: proc
    type(t_perturbation), intent(inout) :: pert
    type(silja_interval), intent(in) :: timestep

    real :: height, bottom, top, vol_sigma, height_exp, height_scaling, height_sigma, vol_exp
    character(len=*), parameter :: sub_name = 'perturb_emis_proc_volc'
    real, dimension(pert%pert_size) :: species_flux                                                                                 
    real, dimension(1) :: sample                                                                                                                                                                             
    integer :: num_levs, ind_lev, num_species, ind_species, ind_proc
    !real, parameter :: fract_height_top = 0.25, fract_height_stem = 1-fract_height_top, &                                                                                                                         
    !     & fract_mass_top = 0.75, fract_mass_stem = 1-fract_mass_top                                                                                                                                              
    real :: mass, fract_stem, fract_top, dz_stem, dz_top, vol_flux, top_stem 
    real, dimension(:), pointer :: proc_data
    type(silam_vertical) :: vert
    real :: alpha, corr_time
    real :: fract_height_top, fract_mass_top, fract_height_stem, fract_mass_stem, fract_height_sigma

    vert = fu_vertical(proc)
    if (fu_fails(defined(vert), 'processor has no vertical', sub_name)) return
    bottom = fu_level_height(fu_lower_boundary_of_layer(fu_level(vert, 1)))
    if (error) return
    num_levs = fu_NbrOfLevels(vert)
    top = fu_level_height(fu_lower_boundary_of_layer(fu_level(vert, num_levs)))

    vol_sigma = pert%storage(3)
    height_exp = pert%storage(4) ! 0.241 for Mastin 2009                                                                                                                                                           
    fract_mass_top = pert%storage(6)
    fract_height_top = pert%storage(7)

    fract_mass_stem = 1 - fract_mass_top

    call random_normal(sample)

    ! The model for rate perturbations in in form 10**X, where X ~ N(0,s^2). Convert this                                                                                                                         
    ! into the form exp(X') so we can use usual log-normal statistics with X' ~ N(0, ln10^2*s^2).                                                                                                                  
    sample(1) = sample(1) * vol_sigma

    if (pert%first_time) then
      pert%storage(1) = sample(1)
      pert%first_time = .false.
    else
      corr_time = pert%storage(2)
      alpha = 1 - fu_sec(timestep) / corr_time
      call msg('zt_param: alpha, sample', alpha, sample(1))
      pert%storage(1) = alpha*pert%storage(1) + sqrt(1-alpha**2) * sample(1)
    end if

    vol_exp = pert%storage(1) + pert%offset(1)
    height_scaling = pert%storage(5)

    fract_height_stem = 1 - fract_height_top

    call msg('volcano exponent:', vol_exp)
    call msg('height exponent:', height_exp)
    call msg('height scaling:', height_scaling)
    call msg('fract height top:', fract_height_top)
    call msg('fract mass top:', fract_mass_top)

    vol_flux = 10**(vol_exp)/(3600*24)
    species_flux = vol_flux*pert%scale
    height = height_scaling * vol_flux**height_exp

    call msg('Eruption height:', height)
 
    call msg('Volume flux m3/sec', vol_flux)
    call msg('Mass flux kgDRE/sec', vol_flux*2500)
    call msg('Species flux kg|mol / sec ', species_flux*64e-3)
 
    top_stem = fract_height_stem * height

    ! vertical attribution                                                                                                                                                                                         
    num_species = size(fu_species_list(proc))
    proc_data => fu_data(proc)

    if (fu_fails(size(proc_data) == num_species*num_levs, 'proc is wrong size', sub_name)) return
    do ind_species = 1, num_species
      !call msg('species:', ind_species)                                                                                                                                                                           
      do ind_lev = 1, num_levs
        top = fu_level_height(fu_upper_boundary_of_layer(fu_level(vert, ind_lev)))
        bottom = fu_level_height(fu_lower_boundary_of_layer(fu_level(vert, ind_lev)))

        if (bottom > height) then
          dz_top = 0.0
          dz_stem = 0.0
        else if (bottom > top_stem) then
          dz_stem = 0.0
          dz_top = min(height-top_stem, top-bottom)
        else if (top > top_stem) then
          dz_stem = top_stem - bottom
          dz_top = top - top_stem
        else
          dz_stem = top-bottom
          dz_top = 0.0
        end if

        ! fractions of stem and top emission that go to this layer                                                                                                                                                 
        dz_top = min(fract_height_top*height, dz_top)
        dz_stem = min(fract_height_stem*height, dz_stem)

        if (height < 1e-5) then
          if (ind_lev == 1) then
            fract_top = 1.0
            fract_stem = 0.0
          else
            fract_top = 0.0
            fract_stem = 0.0
          end if
        else
          if (fract_height_stem > 0) then
            fract_stem = dz_stem / (fract_height_stem*height)
          else
            fract_stem = 0.0
          end if
          if (fract_height_top > 0) then
            fract_top = dz_top / (fract_height_top*height)
          else
            fract_top = 0.0
          end if
        end if
                                                                                                                                                                          
        mass = (fract_stem*fract_mass_stem + fract_top*fract_mass_top) * species_flux(ind_species)
        !call msg('mass to layer', ind_lev, mass)                                                                                                                                                                  
        ind_proc = ind_species + (ind_lev-1)*num_species
        proc_data(ind_proc) = mass / tm_hgt_force_scaling
      end do
    end do

  end subroutine perturb_emis_proc_volc    
  
  !************************************************************************************

  subroutine get_perturbation(rule, values)
    ! Get a gaussian 2d or 3d random field.
    implicit none
    type(t_perturbation), intent(in) :: rule
    real, dimension(:), intent(out) :: values
    
    real, dimension(:), pointer :: work
    integer :: ii, size_uncorr, nn
    real, parameter :: std_logn = sqrt((exp(1.0) - 1.0)*exp(1.0)), exp_05 = exp(0.5)
    real :: new_mean
    
    if (fu_fails(defined(rule), 'Perturbation rule not defined', 'get_perturbation')) return
    if (fu_fails(rule%pert_size <= size(values), 'values too small', 'get_perturbation')) return

    nn = rule%pert_size

    work => fu_work_array(size(values))
    call random_normal(values(1:rule%pert_size))

    if (defined(rule%correlation)) then
      size_uncorr = fu_dimension_control(rule%correlation) ! just to avoid error message
      call transf_fwd(rule%correlation, values(1:size_uncorr), work(1:rule%pert_size))
      if (error) return
      values = work(1:size(values))
    end if

    if (rule%distribution == distr_normal) then
      do ii = 1, rule%pert_size
        values(ii) = values(ii) * rule%scale(ii) + rule%offset(ii)
      end do

    else if (rule%distribution == distr_normal_squared) then
      do ii = 1, rule%pert_size
        values(ii) = rule%scale(ii) * values(ii)**2 + rule%offset(ii)
      end do

    else if (rule%distribution == distr_normal_posit) then
      ! doesn't work don't use
      call random_normal(values(1:rule%pert_size))
      do ii = 1, rule%pert_size
        values(ii) = rule%scale(ii)*values(ii) + rule%offset(ii)
      end do
      ! reject negative values
      do while (any(values(1:rule%pert_size) < 0)) 
        call random_normal(work(1:rule%pert_size))
        do ii = 1, rule%pert_size
          if (values(ii) < 0) values(ii) = work(ii)*rule%scale(ii) + rule%offset(ii)
        end do
      end do
      if (defined(rule%correlation)) then
        call msg('Applying correlation operator')
        size_uncorr = fu_dimension_control(rule%correlation) ! just to avoid error message
        call transf_fwd(rule%correlation, values(1:size_uncorr), work(1:rule%pert_size))
        if (error) return
        values = work(1:size(values))
      end if

    else if (rule%distribution == distr_lognormal) then
      call msg('distr_lognormal')

      values(1:nn) = rule%offset(1:nn) + values(1:nn) * rule%scale(1:nn)
      values(1:nn) = exp(values(1:nn))

    else
      call set_error('Bad distribution', 'get_perturbation')
      return
    end if

    call free_work_array(work)

  end subroutine get_perturbation

  
  !*****************************************************************************************
  
  subroutine store_perturbation_val(perturbations, pert_type, inData, nData, nSpecies)
    !
    ! Stores the perturbation values as a result of control-2-cloud conversion
    !
    implicit none
    
    ! imported parameters
    type(t_perturbation), dimension(:), intent(inout) :: perturbations
    integer, intent(in) :: pert_type, nData, nSpecies
    real, dimension(:), intent(in) :: inData
    
    ! local variables
    integer :: indPert

    indPert = fu_index(pert_type, perturbations(:)%perturb_type)
    if(fu_fails(indPert /= int_missing, 'Failed to find perturbation index:' + fu_str(pert_type), 'store_perturbation_val'))return

    select case(pert_type) 
      case(perturb_emis_xy)
        perturbations(indPert)%storage = inData(1:nData)
        
      case(perturb_emis_volc)  ! some unpacking is needed
        if (fu_fails(nData == nSpecies*2 + fu_num_params_volc(), 'bad nData', 'store_perturbation_val')) return
        perturbations(indPert)%scale(1:nSpecies)  = inData(1:nSpecies) 
        perturbations(indPert)%offset(1:nSpecies)  = inData(nSpecies+1:2*nSpecies)
        perturbations(indPert)%storage(1:fu_num_params_volc()) = inData(2*nSpecies+1:nData)
      case default
        call set_error('Unknown perturbation type:' + fu_str(pert_type),'store_perturbation_val')
        return
    end select
    
  end subroutine store_perturbation_val

  
  !*****************************************************************************************
  
  subroutine retrieve_perturbation_val(perturbations, pert_type, outData, nData, nSpecies)
    !
    ! Retrieves the perturbation values as a result of control-2-cloud conversion
    !
    implicit none
    
    ! imported parameters
    type(t_perturbation), dimension(:), intent(in) :: perturbations
    integer, intent(in) :: pert_type, nData, nSpecies
    real, dimension(:), intent(out) :: outData
    
    ! local variables
    integer :: indPert

    indPert = fu_index(pert_type, perturbations(:)%perturb_type)
    if(fu_fails(indPert /= int_missing, 'Failed to find perturbation index:' + fu_str(pert_type), 'retrieve_perturbation_val'))return

    select case(pert_type) 
      case(perturb_emis_xy)
        outData(1:nData) = perturbations(indPert)%storage(1:nData)
        
      case(perturb_emis_volc)  ! some unpacking is needed
        if (fu_fails(nData == nSpecies*2 + fu_num_params_volc(), 'bad nData', 'retrieve_perturbation_val')) return
        outData(1:nSpecies)  = perturbations(indPert)%scale(1:nSpecies)
        outData(nSpecies+1:2*nSpecies) = perturbations(indPert)%offset(1:nSpecies)
        outData(2*nSpecies+1:nData) = perturbations(indPert)%storage(1:fu_num_params_volc())
      case default
        call set_error('Unknown perturbation type:' + fu_str(pert_type),'retrieve_perturbation_val')
        return
    end select
    
  end subroutine retrieve_perturbation_val

  !*******************************************************************************************

  subroutine sqrt_cov(pert, stdev, dist_0, grid, max_corr)
    implicit none
    type(t_perturbation), intent(inout) :: pert
    type(silja_grid), intent(in) :: grid
    real, intent(in) :: stdev, dist_0, max_corr                                                                                                                                                                                                                                                                          

    real, dimension(:,:), allocatable :: cov
    real, dimension(:,:), pointer :: L
    integer :: i, j, k, n, stat
    character(len=*), parameter :: sub_name = 'sqrt_cov'
    real :: d
    real, dimension(:), pointer :: lats, lons

    ! Creates a distance-based covariance matrix cov and                                                                                                                                                                                                                                
    ! computes the lower triangular matrix L that satisfies cov = LL^T                                                                                                     

    n = fu_number_of_gridpoints(grid)
    lats => fu_geolats_fld(grid)
    lons => fu_geolons_fld(grid)

    allocate(cov(n, n), stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
    !$OMP PARALLEL DO default(none) private(i, j) &                                                                                                                                                                                                                                                                        
    !$OMP & shared(cov, n, stdev, dist_0, lats, lons, max_corr, L)                                                                                                                                                                                                                                                         
    do i=1, n
      do j=i, n
        if (i /= j) then
          cov(i, j) = stdev**2 * max_corr * exp(-0.5*(fu_gc_distance(lons(i), lons(j), lats(i), lats(j))/dist_0)**2)
        else
          cov(i, j) = 0.5*stdev**2
        end if
      end do
    end do
    !$OMP END PARALLEL DO                                                                                                                                                                                                                                                                                                 
    cov = cov + transpose(cov)
    call msg('Covariance matrix ready')

    L => pert%L
    L(:,:)=0.0
    do k = 1, n
      d = cov(k,k) - dot_product(L(k,1:k-1),L(k,1:k-1))
      if (d <= 0) then
        call msg('Error with Cholesky')
        call set_error('Unable to perform Cholesky decomposition', sub_name)
        return
      end if
      L(k,k) = sqrt(d)

      !$OMP PARALLEL DO default(none) private(i) &                                                                                                                                                                                                                                                                         
      !$OMP & shared(cov, n, L, k)                                                                                                                                                                                                                                                                                         
      do i = k+1, n
        L(i,k)  = (cov(i,k) - dot_product(L(i,1:k-1),L(k,1:k-1)) ) / L(k,k)
      end do
      !$OMP END PARALLEL DO                                                                                                                                                                                                                                                                                                
    end do
    call msg('Square root of the covariance matrix computed')                                                                                                                                                                                           

    deallocate(cov, stat=stat)
    if (fu_fails(stat == 0, 'Deallocate failed', sub_name)) return
  end subroutine sqrt_cov
  
  !*******************************************************************************************
  
  real function fu_eruption_height(perturbations)
    !
    ! Returns the height of eruption if the perturbations invclude volcanic source
    !
    implicit none
    
    ! imported parameters
    type(t_perturbation), dimension(:), pointer :: perturbations
    
    ! Local variables
    integer :: ind_pert_volc
    
    fu_eruption_height = real_missing
    if(.not. associated(perturbations))return
    
    ind_pert_volc = fu_index(perturb_emis_volc, perturbations(:)%perturb_type)
    if (ind_pert_volc /= int_missing) then
      fu_eruption_height = perturbations(ind_pert_volc)%storage(5) * &
                         & (10**(perturbations(ind_pert_volc)%storage(1) + &
                         & perturbations(ind_pert_volc)%offset(1)) / (3600*24))** &
                         & perturbations(ind_pert_volc)%storage(4)
    end if
  end function fu_eruption_height

end module perturbations
