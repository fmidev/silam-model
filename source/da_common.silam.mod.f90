module da_common
  use correlations
  use pollution_cloud
  use cocktail_basic
  use globals
  use observation_server
  use source_apportionment, only : slot_start, slot_end, &
       & processor_tm_hgt_force, processor_tm_hgt_force_wgt, processor_tm_hgt_scale
  use dispersion_server, only : update_mass_map_from_file
  use enkf
  use io_server, only : fu_species_output_name

  implicit none

  interface init_control
     module procedure init_control_from_cloud
  end interface

  interface init_control_array
     module procedure init_control_array_from_cloud
     module procedure init_control_array_from_par
  end interface

  interface destroy
     module procedure destroy_control
  end interface

  interface fu_size
     module procedure fu_size_control
  end interface

  interface fu_mode
     module procedure fu_mode_control
  end interface

  interface defined
     module procedure fu_defined_control
     module procedure fu_defined_bgr_cov
  end interface

  interface report
     module procedure report_control
  end interface

  interface get_localisation
     module procedure get_localisation_control
  end interface
  public get_localisation

  public init_control_array
  public fu_mode
  public defined
  public report
  public fu_size
!  public init_control
  public init_control_from_par
  public copy_control
  public fu_n_emission
  public fu_n_initial
  public fu_in_physical_space
  !public set_physical_space
  !public set_control_space

  public pick_species

  public fu_values_ptr
  public fu_initial_ptr
  public fu_emission_ptr
  public fu_species_initial
  public fu_species_emission_ctrl
  public fu_nsp_initial     !number
  public fu_nsp_emis_ctrl
  public fu_isp_initial    !Indices
  public fu_isp_emis_ctrl
  public report_norm
  public set_bgr_cov_ptr
  public get_emis_time_slots
  public fu_num_emis_time_slots
  public get_posit_constr
  public to_control
  public to_physical
  public recover_control
  
  public fu_have_initial
  public fu_have_emission_xy
  public fu_have_emission_zt
  public fu_have_emission_volc
!  public fu_control_includes

  public get_fieldvars_for_control
  public get_non_fieldvars_for_control
  public control_to_fieldset_get_sizes
  public control_to_fieldset
  public control_non_fields_to_files
  
  public fu_num_params_volc

  ! options for the DA output content
  !
  integer, public, parameter :: da_out_none_flag = 190001         ! silent DA, not even bacground and analysis, just forecast after DA
  integer, public, parameter :: da_out_first_last_flag = 190002   ! background and analysis in ^/iteration directory
  integer, public, parameter :: da_out_trajectory_flag = 190003   ! all steps of DA in ^/iteration directory
  integer, public, parameter :: da_out_full_dump_flag = 190004    ! all DA steps and, for 4D-Var, all time evolution of all runs
  !
  ! Meteorology can be perturbed for EnKF and EnKS. This type defines, how
  !
  type Tperturb_meteo
    type(silja_interval) :: time_shift_stdev, max_time_shift  ! standard deviation and max (can be undefined) allowed time shift
    integer :: time_shift_distribution = int_missing         ! e.g. gaussian, ...
    type(silja_logical) :: defined = silja_false
  end type Tperturb_meteo
  !
  ! Some basic numerical parameters of DA
  !
  type DA_numericalParameters
     ! line search method, steepest descent only
     integer :: lineSearchMethod = LINE_SEARCH_MS
     integer :: searchMethod = steepest_descent_flag
     ! steepest descent only:
     real :: minimumStepLength = 1e-9
     ! stop if cost function reduced by this factor
     real :: cost_rel_tol = 1e-6
     ! stop if relative change in cost function below this
     real :: cost_incr_rel_tol = 1e-6
     ! stop if gradient reduced by this factor
     real :: grad_rel_tol = 1e-2
     integer :: maxIterations = 40
     integer :: maxLineSearchSteps = 25
     real :: stepIncreaseCoefMS = 3.0
     real :: stepDecreaseCoefMS = 5.1
     real :: quasi_newton_df1 = 1.0
  end type DA_numericalParameters
  !
  ! The main set of rules for performing the data assimilation
  !
  type DA_rules
     integer :: method = int_missing ! 3d|4d
     type(silja_time) :: DA_begin, first_analysis_time, last_analysis_time
     type(silja_interval) :: da_window ! assimilation window, 4dvar: periodToCompute,  (zero for 3dvar)
     type(silja_interval) :: assim_interval ! 3dvar (or seq. 4dvar): time between assimilations
     logical :: ifRandomise
     logical :: defined = .false.
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
     logical :: use_log_cloud = .false.
     logical :: use_log_obs = .false.
     real :: observation_stdev = -1.0
     real :: emission_stdev = 0.6
     real :: emis_corr_dist = -1.0
     real :: emis_max_corr = -1.0
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!111
     ! Introducing support for arbitrary mixes of control variables:
     ! define the control variable as a bitmask - see above.
     integer :: controlVariable = int_missing
     integer:: perturbVariable = int_missing 
     real :: tolerance, weightNonzero = 1, weightZero = 1 
     
     type(grads_template) :: cov_setup_templ_init, cov_setup_templ_emis
     character(len=fnlen) :: station_path
     ! specifications for observations
     character(len=fnlen), dimension(:), pointer :: obs_items
     ! a pointer to file with the specifications above
     character(len=fnlen) :: obs_list_path = ''
     integer :: num_obs_items
     type(grads_template), dimension(:), pointer :: obs_templ_list
     ! 4dvar: analysis species always full transport/emission, 3dvar - can be a subset of transport.
     character(len=substNmLen), dimension(:), pointer :: analysis_subst_list_3d
     !integer, dimension(:), pointer :: ind_species_in_transp 
     type(silja_logical) :: have_analysis_species = silja_undefined
     ! 4dvar for time/height profile:
     !type(silja_time), dimension(:,:), pointer :: time_slots
     type(silja_interval) :: emis_time_slot
     ! forcing/scaling for time height emission control variable.
     ! either processor_tm_hgt_scale or processor_tm_hgt_force
     integer :: time_height_mode = int_missing
     ! unit release when computing the H-matrix
     real :: h_matrix_unit = 1.0

     character(len=fnlen) :: outputDir = char_missing
     character(len=fnlen) :: outdir_templ_str

     type(DA_numericalParameters) :: numerics
     integer :: output_level = da_out_first_last_flag, ind_obs

     ! Allow negatives in mass map. Only for 3dvar, which doesn't run the model.
     logical :: allow_negative = .false.
     
     ! background fields
     character(len=fnlen) :: initialStateBgrFile = ''
     character(len=fnlen) :: emissionCorrectionBgrFile = ''

     logical :: have_init_background = .false., have_emis_background = .false.
     integer :: adjointMethod = int_missing

     logical :: log_transform = .false.

     ! The background term can be disabled altogether (in a kludgy way). This is done by
     ! removing the contribution from gradient and cost function, and by setting the
     ! background value to 0 and standard deviation to 1.0 so that the control-to-physical
     ! transformation becomes identity.
     logical :: no_background_term = .false.
     logical :: use_zero_emis_backgr = .false.
     ! Option to restart iteration: search the last control_??.grads and initialize
     ! iteration from that.
     logical :: restart_iteration = .false.
     ! 4Dvar: automatically run forecast after assimilation
     logical :: run_forecast = .true.

     ! enkf: ensemble size
     integer :: ens_size = int_missing

     ! enkf: correlation time for general emissions                   
     real :: emis_corr_time = 30*60.0

     ! enkf: parametric volcano emission perturbation

     ! new version: scaling from total volume to species mass, [mass unit] / [m3 tephra]
     real :: volc_massflux_scaling = 1.0
     ! correlation time for logarithm of volume flux
     real :: volc_volflux_corr_time = 30*60.0
     
     ! standard deviation of the logarithm of total volume flux
     real :: volc_volflux_sigma = 2.0
     ! the mode of the lognormal distribution for total volume flux, m3/s
     real :: volc_volflux_mode = real_missing 
     ! scale and exponent in the mastin-type power law for computing eruption height
     real :: volc_height_exp = 0.241, volc_height_scaling = 2.0 * 1e3 ! to meters, mastin has km
     ! sigma (meters) for eruption height deviation of the Mastin fit
     real :: volc_height_sigma = 0.0
     ! fraction of mass emitted to the top of eruption embrella
     real :: volc_fract_mass_top = 0.75
     ! fraction of height covered by the top
     real :: volc_fract_height_top = 0.25
     ! deviation of height covered by top
     real :: volc_fract_height_sigma = 0.0
     
     ! enkf filter options
     ! localisation function
     integer :: loc_type = int_missing
     ! enkf or denkf ("deterministic enkf")
     integer :: enkf_flavor = int_missing
     ! localisation distance
     real :: loc_dist_m = real_missing
     ! inflation factor
     real :: rfactor = real_missing
     ! localisation of the z-t emission 
     real :: zt_loc_lon = real_missing, zt_loc_lat = real_missing
     
     ! If set, only master task does the analysis computations
     logical :: enkf_master_only = .false.

     ! EnKS setup. If the EnKS method is used, also past times can be included in the
     ! state vector by setting a smoother_step and one or more smoother_output_steps.
     ! smoother_output_step = 0 means analysis time - assimilation interval.
     ! Later steps are by adding smoother_steps to step 0.
     type(silja_interval) :: smoother_step = interval_missing
     integer, dimension(:), pointer :: smoother_output_steps
     integer :: num_smoother_output_steps = 0
     
     ! A group of parameters controlling meteorological perturbations
     type(Tperturb_meteo), allocatable :: perturb_meteo

  end type DA_rules

  !-------------------------------------------------------------------------------------

  type da_control
     ! 
     ! Container for the control variable. Consider the following:
     ! - the values are all in one real array. The content is split (currently) into
     ! emission and initial state components. There are functions that return pointers to
     ! these, as real arrays.
     ! - A control can be either in physical or control space - the mapping between these
     ! with an affine transformation (see eg. the 3d-var modules). The transf_state flag
     ! allows guarding that we are in the correct representation before using the values.
     !
     private
     real, dimension(:), pointer :: vals
     integer :: nvals = int_missing
     integer :: mode = int_missing
     integer :: transf_state = int_missing
     integer :: nx, ny, nz_init, nz_emis, nsp_init, nsp_emis
     integer :: shift_initial = int_missing, shift_emission = int_missing
     integer, dimension(max_species) :: ind_species_emis, ind_species_init
     type(silam_species), dimension(max_species) :: species_initial, species_emission
     integer :: size_initial = 0, size_emission = 0
     type(t_background_covariance), pointer :: bgr_cov_ptr => null()! only if in control space
     type(silja_grid) :: grid
     ! there's no vertical. Currently not needed...should probably be pointer if
     ! introduced, but that has its own complications.
     logical :: defined = .false.
  end type da_control
  public da_control
  
  type da_control_ptr
     type(da_control), pointer :: ptr
  end type da_control_ptr
  public da_control_ptr

  type(da_control), public, save :: control_missing ! with default initialization

  !---------------------------------------------------------------------------------------
  
  type t_background_covariance
     !private
     real, dimension(:), pointer :: stdev_initial, stdev_emission

     type(t_correlation) :: correlation_initial, correlation_emission

     logical :: defined = .false.
     
  end type t_background_covariance

  public t_background_covariance
  type (t_background_covariance), target :: background_covariance_missing !Should get all undefined 

  private init_control_from_cloud
  private destroy_control
  
  ! Values for the transf_state flag above.
  integer, parameter, public :: physical_space = 1001, control_space = 1002

  ! Options for type-flag in t_background_covariance
  !
  integer, parameter, public :: constant = 20001 ! constant variance for each species
  integer, parameter, public :: variable_1d = 20002 ! variance depending on level and species
  ! variance depending on location and species but not level
  integer, parameter, public :: variable_2d = 20003
  integer, parameter, public :: fraction_of_background = 20004
  ! in the future:
  integer, parameter, public :: homogeneous_gaussian = 20005 ! gaussian spatial correlation
  integer, parameter, public :: diagonal = 20006, & ! no spatial correlation
                               & diagonal_logarithmic = 20007

  real, dimension(:), private, allocatable, save, target ::  transf_work

  ! The control variable mode. Originally, we had separate flags for each combination of
  ! control (emission, initial state, etc). Now this is done more flexibly and each
  ! possible component can be enabled or disabled (although some combinations are still
  ! forbidden). This is handled by the bitmasks below.
  integer, parameter, public :: control_none = 0, control_init = ibset(0,0), &
                                                & control_emis_xy = ibset(0,1), &
                                                & control_emis_zt = ibset(0,2), &
                                                & control_emis_volc = ibset(0,3), &
                                                & control_meteo_time = ibset(0,4)
  ! "IBSET(I,POS) returns the value of I with the bit at position POS set to one."

  ! for backwards compatibility, still used by some subroutines.
  integer, parameter, public :: da_initial_state = control_init, &
                              & da_emission_correction = control_emis_xy, &
                              & da_emission_and_initial = ior(control_init, control_emis_xy), &
                              & da_emission_time_height = control_emis_zt

  CONTAINS

  logical function fu_defined_control(control) result(is_defined)
    implicit none
    type(da_control), intent(in) :: control
    is_defined = control%defined
  end function fu_defined_control

  !************************************************************************************

  subroutine report_control(control)
    implicit none
    type(da_control), intent(in) :: control

    integer :: n_data
    real, dimension(:), pointer :: p_data

    call msg('Reporting control..-')
    if (.not. defined(control)) then
      call msg('Undefined control')
      return
    end if
    call msg('Transformation state:', control%transf_state)

    p_data => fu_values_ptr(control)
    
    select case (control%transf_state)
    case (physical_space)
      n_data = fu_size(control)
    case (control_space)
    end select

    n_data = fu_size(control)
    call msg('Whole control:')
    call msg('  min:', minval(p_data(1:n_data)))
    call msg('  max:', maxval(p_data(1:n_data)))
    call msg('  mean:', sum(p_data(1:n_data))/n_data)
    call msg('  norm:', sqrt(sum(p_data(1:n_data)**2)))

    if (control%shift_initial /= int_missing) then
      p_data => fu_initial_ptr(control)
      n_data = fu_n_initial(control)
      call msg('Initial state component:')
      call msg('  min:', minval(p_data(1:n_data)))
      call msg('  max:', maxval(p_data(1:n_data)))
      call msg('  mean:', sum(p_data(1:n_data))/n_data)
      call msg('  norm:', sqrt(sum(p_data(1:n_data)**2)))
    end if
    if (control%shift_emission /= int_missing) then
      p_data => fu_emission_ptr(control)
      n_data = fu_n_emission(control)
      call msg('Emission component:')
      call msg('  min:', minval(p_data(1:n_data)))
      call msg('  max:', maxval(p_data(1:n_data)))
      call msg('  mean:', sum(p_data(1:n_data))/n_data)
      call msg('  norm:', sqrt(sum(p_data(1:n_data)**2)))
    end if
    
  end subroutine report_control

  !************************************************************************************

  function fu_values_ptr(control) result(p_vals)
    implicit none
    type(da_control), intent(in) :: control
    real, dimension(:), pointer :: p_vals
    p_vals => control%vals
  end function fu_values_ptr

  ! The following pair of routines returns pointers to the array sections corresponding to
  ! the emission and initial state components. The dimensions of these components are
  ! returned by the fu_n_* functions.
  ! 
  function fu_initial_ptr(control, ind_species) result(p_init)
    implicit none
    type(da_control), intent(in) :: control
    real, dimension(:), pointer :: p_init
    integer, intent(in), optional :: ind_species

    integer :: ind_start, ind_end,  field_size
    
    if (fu_fails(control%shift_initial /= int_missing, 'initial_ptr not available', 'fu_intial_ptr')) return

    if (control%nvals == 0) then
      p_init => null()
    else
      if (present(ind_species)) then
        field_size = control%nx * control%ny * control%nz_init
        ind_start = control%shift_initial + (ind_species-1)*field_size
        ind_end = ind_start + field_size - 1
      else
        ind_start = control%shift_initial
        ind_end = ind_start + fu_n_initial(control) - 1
      end if

      p_init => control%vals(ind_start:ind_end)
    endif
  end function fu_initial_ptr
    
  function fu_emission_ptr(control, ind_species) result(p_emis)
    implicit none
    type(da_control), intent(in) :: control
    integer, intent(in), optional :: ind_species
    real, dimension(:), pointer :: p_emis

    integer :: ind_start, ind_end, field_size, nx, ny
    character(len=*), parameter :: subname = 'fu_emission_ptr'
    if (fu_fails(control%shift_emission /= int_missing, 'emission_ptr not available', subname)) return

    if (fu_have_emission_zt(fu_mode(control))) then
      nx = 1
      ny = 1
    else
      nx = control%nx
      ny = control%ny
    end if

    if (present(ind_species)) then
      if (fu_have_emission_volc(fu_mode(control))) then
        call set_error('ind_species not allowed with emission_volc', subname)
        return
      end if
      field_size = nx*ny*control%nz_emis 
      ind_start = control%shift_emission + (ind_species-1)*field_size
      ind_end = ind_start + field_size - 1
    else
      ind_start = control%shift_emission
      ind_end = ind_start + fu_n_emission(control) - 1
    end if
    
    p_emis => control%vals(ind_start:ind_end)
  end function fu_emission_ptr

  !************************************************************************************

  subroutine get_localisation_control(control, lonlat, darules)
    ! The control vector will be localised horizontally and used for EnKF local analysis.
    implicit none
    type(da_control), intent(in) :: control
    real, dimension(:,:), intent(out) :: lonlat
    type(da_rules) :: darules

    integer :: shift, count
    character(len=*), parameter :: sub_name = 'get_localisation_control'

    if (fu_fails(size(lonlat, 1) > 1, 'lonlat 1st dim too small', sub_name)) return 
    if (fu_fails(defined(control), 'control not defined', sub_name)) return

    if (fu_have_initial(control%mode)) then
      shift = control%shift_initial
      count = fu_n_initial(control)
      if (error) return
      if (fu_fails(size(lonlat,2) >= shift+count-1, 'lonlat too small', sub_name)) return
      call get_grid_loc(control%grid, control%nz_init, control%nsp_init, &
                      & lonlat(:,shift:shift+count-1))
      if (error) return
    end if
    
    if (fu_have_emission_xy(control%mode)) then
      shift = control%shift_emission
      count = fu_n_emission(control)
      if (fu_fails(size(lonlat,2) >= shift+count-1, 'lonlat too small', sub_name)) return
      call get_grid_loc(control%grid, control%nz_emis, control%nsp_emis, &
                      & lonlat(:,shift:shift+count-1))
      if (error) return
    end if
    
    if (fu_have_emission_zt(fu_mode(control)) .or. fu_have_emission_volc(fu_mode(control))) then
      call msg('localisation for control: emission zt/volc', darules%zt_loc_lon, darules%zt_loc_lat)
      shift = control%shift_emission
      count = fu_n_emission(control)
      if (fu_fails(size(lonlat,2) >= shift+count-1, 'lonlat too small', sub_name)) return
      lonlat(1, shift:shift+count-1) = darules%zt_loc_lon
      lonlat(2, shift:shift+count-1) = darules%zt_loc_lat
      if (error) return
    end if

  contains
    
    !=============================================================
  
    subroutine get_grid_loc(grid, nz, num_species, loc)
      ! Set the lon/lat localisation for a control vector in mass map order.
      implicit none
      type(silja_grid), intent(in) :: grid
      integer, intent(in) :: nz
      integer, intent(in) :: num_species
      real, dimension(:,:), intent(out) :: loc

      integer :: nx, ny
      integer :: ix, iy, iz, isp, i1d

      call grid_dimensions(grid, nx, ny)

      !if (fu_fails(nz*nx*ny*num_species <= size(loc,2), 'loc too small', 'get_grid_loc')) return

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

  end subroutine get_localisation_control


  !************************************************************************************

  subroutine init_control_from_cloud(control, cloud, rules, space_flag, initial_value, bgr_cov)
    implicit none
    type(da_control), intent(out) :: control
    type(silam_pollution_cloud), target, intent(in) :: cloud
    type(da_rules), intent(in) :: rules
    ! Specify whether the control is initially in physical or control space:
    integer, intent(in) :: space_flag
    real, intent(in), optional :: initial_value
    type(t_background_covariance), intent(in), optional, target :: bgr_cov

    type(Tmass_map), pointer :: p_map_transp
    type(silam_pollution_cloud), pointer :: cloud_ptr
    type(silam_species), dimension(:), pointer :: transp_species, emis_species
    type(silja_grid) :: analysis_grid
    type(silam_vertical) :: analysis_vertical
    type(t_background_covariance), pointer :: bgr_cov_ptr
    real :: initial_value_

    cloud_ptr => cloud
    p_map_transp => fu_concMM_ptr(cloud)

    if (smpi_adv_rank == 0) then
      analysis_grid = wholeMPIdispersion_grid !p_map_transp%gridTemplate
    else
      analysis_grid = grid_missing !Make da_control with no data
    endif

    analysis_vertical = p_map_transp%vertTemplate
    transp_species => fu_species_transport(cloud_ptr)
    emis_species => fu_species_emission(cloud_ptr)

    bgr_cov_ptr => background_covariance_missing
    if (present(bgr_cov)) then
      bgr_cov_ptr => bgr_cov
    endif

    initial_value_ = real_missing
    if (present(initial_value)) then
      initial_value_ = initial_value
    endif

    call init_control_from_par(control, analysis_grid, analysis_vertical, transp_species, emis_species, &
                               & rules, space_flag, initial_value_, bgr_cov_ptr)
    
    call set_missing(analysis_vertical, .False. )

  end subroutine init_control_from_cloud

  
  !************************************************************************************

  subroutine init_control_array_from_cloud(control_arr, cloud, rules, space_flag, num_controls, storage, &
                                         & initial_value, bgr_cov)
    implicit none
    type(da_control), dimension(:), intent(out) :: control_arr
    type(silam_pollution_cloud), target, intent(in) :: cloud
    type(da_rules), intent(in) :: rules
    ! Specify whether the control is initially in physical or control space:
    integer, intent(in) :: space_flag, num_controls
    real, dimension(:), pointer :: storage
    
    real, intent(in), optional :: initial_value
    type(t_background_covariance), intent(in), optional, target :: bgr_cov

    type(Tmass_map), pointer :: p_map_transp
    type(silam_pollution_cloud), pointer :: cloud_ptr
    type(silam_species), dimension(:), pointer :: transp_species, emis_species
    type(silja_grid) :: analysis_grid
    type(silam_vertical) :: analysis_vertical
    type(t_background_covariance), pointer :: p_bgr_cov
    real :: initial_value_

    cloud_ptr => cloud
    p_map_transp => fu_concMM_ptr(cloud)
    
    if (smpi_adv_rank==0) then
       analysis_grid = wholeMPIdispersion_grid
    else
      analysis_grid = grid_missing
    endif
    analysis_vertical = p_map_transp%vertTemplate
    transp_species => fu_species_transport(cloud_ptr)
    emis_species => fu_species_emission(cloud_ptr)

    if (present(bgr_cov)) then
      p_bgr_cov => bgr_cov
    else
      if (fu_fails(space_flag == physical_space, 'control space requires bgr_cov', 'init_control')) return
      ! avoid instant crash calling below
      p_bgr_cov => background_covariance_missing
    end if

    initial_value_ = real_missing
    if (present(initial_value)) then
          initial_value_ = initial_value
    endif

    call init_control_array_from_par(control_arr, analysis_grid, analysis_vertical, &
                                     & transp_species, emis_species, &
                                     & rules, space_flag, num_controls, storage, initial_value, p_bgr_cov)

  end subroutine init_control_array_from_cloud


  !************************************************************************************

  subroutine init_control_from_par(control, analysis_grid, analysis_vertical, &
                                 & p_species_transport, p_species_emission, rules, &
                                 & space_flag, initial_value, bgr_cov, controlvar)
    implicit none
    type(da_control), intent(out) :: control
    type(silja_grid), intent(in) :: analysis_grid
    type(silam_vertical), intent(in) :: analysis_vertical
    type(silam_species), dimension(:), pointer :: p_species_transport, p_species_emission
    type(da_rules), intent(in) :: rules
    ! Specify whether the control is initially in physical or control space:
    integer, intent(in) :: space_flag

    real, intent(in) :: initial_value
    type(t_background_covariance), intent(in), optional, target :: bgr_cov
    integer, intent(in), optional :: controlvar ! to override what is in the rules
    
    integer :: controlvar_, size_req, stat
    real, dimension(:), pointer :: storage
    character(len=*), parameter :: subname = 'init_control_from_par'
    
    if (present(controlvar)) then
      controlvar_ = controlvar
    else
      controlvar_ = rules%controlVariable
    end if
    
    
    ! Initialization in two parts: first do everything but allocate, then allocate here,
    ! then finalize with set_control_storage. The reason is to enable allocation of
    ! several control variables sharing a contiguous storage vector (see init_control_array).
    ! 
    call init_control_no_alloc(control, analysis_grid, analysis_vertical, &
                             & p_species_transport, p_species_emission, rules, &
                             & space_flag, bgr_cov, controlvar_, size_req)
    if (error) return
    if (size_req > 0) then
      allocate(storage(size_req), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', subname)) return
      if (initial_value /= real_missing) storage = initial_value
      call set_control_storage(control, storage)
    else
      call set_control_storage(control, null())
    endif
    
  end subroutine init_control_from_par

  !************************************************************************************

  subroutine init_control_array_from_par(control_arr, analysis_grid, analysis_vertical, &
                                       & p_species_transport, p_species_emission, rules, &
                                       & space_flag, num_controls, storage, &
                                       & initial_value, bgr_cov, controlvar)
    !
    ! Initialize num_controls identical control variables with a joint, contiguous storage
    ! vector.
    ! 
    implicit none
    type(da_control), dimension(:), intent(out) :: control_arr
    type(silja_grid), intent(in) :: analysis_grid
    type(silam_vertical), intent(in) :: analysis_vertical
    type(silam_species), dimension(:), pointer :: p_species_transport, p_species_emission
    type(da_rules), intent(in) :: rules
    ! Specify whether the control is initially in physical or control space:
    integer, intent(in) :: space_flag
    integer, intent(in) :: num_controls
    real, dimension(:), pointer :: storage
    
    real, intent(in), optional :: initial_value
    type(t_background_covariance), intent(in), target :: bgr_cov
    integer, intent(in), optional :: controlvar ! to override what is in the rules
    
    integer :: controlvar_, size_req, stat, ind_control, ind_start
    type(t_background_covariance), pointer :: bgr_cov_
    type(t_background_covariance), target :: bgr_cov_missing
    real, dimension(:), pointer :: storage_subset
    character(len=*), parameter :: subname = 'init_control_arr'
    
    if (present(controlvar)) then
      controlvar_ = controlvar
    else
      controlvar_ = rules%controlVariable
    end if
    
    
    nullify(storage)
    if (fu_fails(num_controls > 0, 'Bad num_controls', subname)) return
    if (fu_fails(size(control_arr) >= num_controls, 'control_arr too small', subname)) return

    ! Use the two-step initialization (see init_control_from_par).
    ! 
    do ind_control = 1, num_controls
      call msg('init control arr:', ind_control)
      call init_control_no_alloc(control_arr(ind_control), analysis_grid, analysis_vertical, &
                               & p_species_transport, p_species_emission, rules, &
                               & space_flag, bgr_cov, controlvar_, size_req)
      if (error) return
      if (.not. associated(storage)) then
        if (fu_fails(size_req > 0, 'Bad size_req', subname)) return
        allocate(storage(num_controls*size_req), stat=stat)
        if (initial_value /= real_missing) storage = initial_value
        if (fu_fails(stat == 0, 'Allocate failed', subname)) return
      end if
      ind_start = (ind_control-1)*size_req + 1
      call msg('::start:', ind_start)
      call msg('::end:', ind_start+size_req-1)

      storage_subset => storage(ind_start:ind_start+size_req-1)
      call set_control_storage(control_arr(ind_control), storage_subset)
      if (error) return
    end do
    
  end subroutine init_control_array_from_par


  !************************************************************************************

  subroutine test_init_control_array(darules, cloud)
    implicit none
    type(silam_pollution_cloud) :: cloud
    type(da_rules) :: darules

    type(da_control) :: controls(3)
    real, dimension(:), pointer :: storage, rp
    character(len=*), parameter :: subname = 'test_init_control_array'
    
    call init_control_array_from_par(controls, wholeMPIdispersion_grid, dispersion_vertical, &
                                   & fu_species_transport(cloud), fu_species_emission(cloud), &
                                   & darules, physical_space, 3, storage, initial_value=0.0, &
                                   & bgr_cov=background_covariance_missing)

    rp => fu_values_ptr(controls(1))
    rp = 1.0
    rp => fu_values_ptr(controls(2))
    rp = 2.0
    if (fu_fails(all(fu_values_ptr(controls(1)) .eps. 1.0), 'bad1', subname)) call unset_error(subname)
    if (fu_fails(all(fu_values_ptr(controls(2)) .eps. 2.0), 'bad2', subname)) call unset_error(subname)
    if (fu_fails(all(fu_values_ptr(controls(3)) .eps. 0.0), 'bad3', subname)) call unset_error(subname)
    call msg(subname // ' finished')
    
  end subroutine test_init_control_array


  !****************************************************************************

  subroutine init_control_no_alloc(control, analysis_grid, analysis_vertical, &
                                 & p_species_transport, p_species_emission, rules, &
                                 & space_flag, bgr_cov, controlvar, size_req)
    ! Internal subroutine. Initialize control without allocation. Vertical and grid
    ! required for dimension, species for chemicals. Control space control needs always a
    ! bgr_cov because the dimensions are not determined without one.
    ! 
    ! The calling subroutine will then obtain the storage needed (size_req) and call
    ! set_control_storage, which finishes the initialization.
    ! 
    implicit none
    type(da_control), intent(out) :: control
    type(silja_grid), intent(in) :: analysis_grid
    type(silam_vertical), intent(in) :: analysis_vertical
    type(silam_species), dimension(:), pointer :: p_species_transport, p_species_emission
    type(da_rules), intent(in) :: rules
    ! Specify whether the control is initially in physical or control space:
    integer, intent(in) :: space_flag
    type(t_background_covariance), intent(in), target :: bgr_cov
    integer, intent(in) :: controlvar
    integer, intent(out) :: size_req

    integer :: nvals, status, nvals_emis, nvals_init, ii
    integer :: nz, stat, nvals_contr, nvals_phys, num_times, shift
    integer, dimension(:), pointer :: indices

    character(len=*), parameter :: subname = 'init_control_no_alloc'

    nz = fu_NbrOfLevels(analysis_vertical)
    if (defined(analysis_grid)) then
      call grid_dimensions(analysis_grid, control%nx, control%ny)
    else
      control%nx = 0
      control%ny = 0
    endif
    control%grid = analysis_grid

    if (defined(bgr_cov)) then
      nvals_contr = fu_control_space_dim(bgr_cov)
    else
      nvals_contr = int_missing
    end if

    control%nsp_init = 0  
    control%nsp_emis = 0

    ! These will be set above for physical space, in control space the values will be
    ! overriden below.
    control%size_initial = 0
    control%size_emission = 0
        
    if (fu_have_emission_xy(controlvar) .and. fu_have_emission_zt(controlvar)) then
      call set_error('Unsupported combination of control variables', subname)
      return
    end if

    if (fu_have_emission_volc(controlvar) .and. ( fu_have_emission_zt(controlvar) &
         & .or. fu_have_emission_xy(controlvar) ) ) then
      call set_error('Unsupported combination of control variables', subname)
      return
    end if

    nvals_phys = 0
    nvals_init = 0
    nvals_emis = 0
    shift = 0

    if (fu_have_initial(controlvar)) then
      call msg('Init control: initial')
      if (fu_true(rules%have_analysis_species)) then
        call pick_species(rules%analysis_subst_list_3d, p_species_transport, &
                        & control%ind_species_init, control%nsp_init)
        if (error) return

      else if (fu_false(rules%have_analysis_species)) then
        control%nsp_init = size(p_species_transport)
        control%ind_species_init(1:control%nsp_init) = (/(ii, ii=1, control%nsp_init)/)
      else
        call set_error('Substance set switch not defined for initial', subname)
        return
      end if ! have analysis species                                                                                                                                                                                                                                                                                        
      control%species_initial = (/(p_species_transport(control%ind_species_init(ii)), &
                                 & ii=1, control%nsp_init)/)
      if (debug_level > 1) then
        call msg('Number of control species for initial:', control%nsp_init)
        call msg('The species are:')
        do ii = 1, control%nsp_init
          call report(control%species_initial(ii))
        end do
      end if
      if (defined(analysis_grid)) then
        call grid_dimensions(analysis_grid, control%nx, control%ny)
        nvals_init = control%nx*control%ny*nz*control%nsp_init
        nvals_phys = nvals_phys + nvals_init
      else
        nvals_init = 0
        nvals_phys = 0
      endif
      control%shift_initial = shift + 1
      shift = shift + nvals_init
      control%size_initial = nvals_init
      control%nz_init = fu_NbrOfLevels(analysis_vertical)
    end if

    if (fu_have_emission_xy(controlvar)) then
      call msg('Init control: emission xy')

      if (fu_fails(.not. fu_have_emission_zt(controlvar), 'Unsupported control variable', subname)) return
      control%nsp_emis = size(p_species_emission)
      if (defined(analysis_grid)) then
        call grid_dimensions(analysis_grid, control%nx, control%ny)
      end if
      nvals_emis = control%nx*control%ny*control%nsp_emis
      nvals_phys = nvals_phys + nvals_emis
      control%shift_emission = shift + 1
      shift = shift + nvals_emis
      control%ind_species_emis(1:control%nsp_emis) = (/(ii, ii=1, control%nsp_emis)/)

      control%species_emission = (/(p_species_transport(control%ind_species_emis(ii)), &
                                  & ii=1, control%nsp_emis)/)
      control%size_emission = nvals_emis
      control%nz_emis = 1
    end if
      
    if (fu_have_emission_zt(controlvar)) then
      call msg('Init control: emission zt')
      if (fu_fails(.not. fu_have_emission_xy(controlvar), 'Unsupported control variable', subname)) return
      control%nsp_emis = size(p_species_emission)
      control%ind_species_emis = (/(ii, ii=1, control%nsp_emis)/)
      control%species_emission = (/(p_species_emission(control%ind_species_emis(ii)), &
                                  & ii=1, control%nsp_emis)/)

      num_times = fu_num_emis_time_slots(rules)
      nvals_emis = num_times*nz*control%nsp_emis
      nvals_phys = nvals_phys + nvals_emis
      control%size_emission = nvals_emis
      control%nz_emis = fu_NbrOfLevels(analysis_vertical)
      control%shift_emission = shift + 1
      shift = shift + nvals_emis
      !control%nx = 1
      !control%ny = 1
    end if

    if (fu_have_emission_volc(controlvar)) then
      call msg('Init control: emission volc')
      if (fu_fails(.not. fu_have_emission_xy(controlvar), 'Unsupported control variable', subname)) return
      control%nsp_emis = size(p_species_emission)
      control%ind_species_emis = (/(ii, ii=1, control%nsp_emis)/)
      control%species_emission = (/(p_species_emission(control%ind_species_emis(ii)), &
                                  & ii=1,  control%nsp_emis)/)
      nvals_emis = fu_num_params_volc()+2*control%nsp_emis
      nvals_phys = nvals_phys + nvals_emis
      control%size_emission = nvals_emis
      control%nz_emis = int_missing
      control%shift_emission = shift + 1
      shift = shift + nvals_emis
    end if

    select case(space_flag)
    case (physical_space)
      control%nvals = nvals_phys
      nullify(control%bgr_cov_ptr)
    case (control_space)
      if (smpi_adv_rank ==0 ) then
        if (fu_fails(defined(bgr_cov), 'Control space requested but no bgr_cov given', subname)) return
      endif
      control%nvals = nvals_contr
      control%bgr_cov_ptr => bgr_cov
      call control_space_dim(bgr_cov, control%size_initial, control%size_emission)
      if (control%size_initial > 0) control%shift_initial = 1
      if (control%size_emission > 0) control%shift_emission = control%size_initial + 1
    case default
      call set_error('Bad space_flag', 'init_control')
      return
    end select
    control%transf_state = space_flag
    
    control%mode = controlvar

    size_req = control%nvals

  end subroutine init_control_no_alloc

                                 
  !************************************************************************************
                               
  subroutine set_control_storage(control, storage)
    ! 
    ! Internal subroutine: after allocating the storage, point the control values into
    ! it. The control is now fully defined.
    !
    implicit none
    type(da_control), intent(inout) :: control
    real, dimension(:), pointer :: storage
    character(len=*), parameter :: subname = 'set_control_storage'

    if (fu_fails(.not. defined(control), 'Control already defined', subname)) return
    if (control%nvals > 0) then
      if (fu_fails(size(storage) == control%nvals, 'Wrong storage size', subname)) return
    endif
    
    control%vals => storage
    control%defined = .true.
    
  end subroutine set_control_storage

  !************************************************************************************

  subroutine pick_species(subst_name_list, full_species_list, indices, num_contr_species)
    implicit none
    character(len=*), dimension(:), intent(in) :: subst_name_list
    type(silam_species), dimension(:), intent(in) :: full_species_list
    integer, dimension(:), intent(out):: indices
    integer, intent(out) :: num_contr_species
    
    integer, dimension(max_species) :: indices_wrk
    integer :: num_for_subst, ind_subst, stat, ind_sp
    
    num_contr_species = 0
    do ind_subst = 1, size(subst_name_list)
      call select_species(full_species_list, size(full_species_list), subst_name_list(ind_subst), &
                        & aerosol_mode_missing, real_missing, indices_wrk(num_contr_species+1:), num_for_subst)
      if (num_for_subst == 0) then
        call msg('Species available:')
        do ind_sp = 1, size(full_species_list)
          call report(full_species_list(ind_sp))
        end do
        call msg('Substance requested: ' // trim(subst_name_list(ind_subst)))
        call set_error('Analysis substance not available', 'init_control')
        return
      end if
      num_contr_species = num_contr_species + num_for_subst
    end do

    indices(1:num_contr_species) = indices_wrk(1:num_contr_species)
    
    
  end subroutine pick_species

  !************************************************************************************

  subroutine copy_control(control_from, control_to, if_allocate)
    implicit none
    type(da_control), intent(in) :: control_from
    type(da_control), intent(inout) :: control_to
    logical, intent(in), optional :: if_allocate
   
    real, dimension(:), pointer :: p_to_tmp
    logical :: if_allocate_loc
    integer :: stat

    if (present(if_allocate)) then
      if_allocate_loc = if_allocate
    else
      if_allocate_loc = .false.
    end if

    if (fu_fails(control_from%defined, 'control_from not defined', 'copy_control')) return
    
    if (if_allocate_loc) then
      allocate(control_to%vals(control_from%nvals), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', 'copy_control')) return
    else
      if (fu_fails(control_to%defined, 'control_to not defined', 'copy_control')) return
      if (fu_fails(control_to%nvals == control_from%nvals, 'controls not same size', 'copy_control')) return
    end if
    
    ! Use the default assignment to copy non-pointer data. Now control_to%vals =>
    ! control_from%vals, but we stored the original pointer - copy the data there and set
    ! control_to%vals to point at it.
    !
    p_to_tmp => control_to%vals
    
    control_to = control_from
!    do stat = 1, control_to%nvals     ! if stack overflow occurs at vector operation, use this cycle
!      p_to_tmp(stat) = control_to%vals(stat)
!    end do
    p_to_tmp(1:control_to%nvals) = control_to%vals(1:control_to%nvals) ! rhs points to the control_from data
    control_to%vals => p_to_tmp
    
  end subroutine copy_control

  !************************************************************************************

  subroutine destroy_control(control)
    implicit none
    type(da_control), intent(inout) :: control

    if (control%nvals > 0)  deallocate(control%vals)
    
    control = control_missing

  end subroutine destroy_control
  
  !************************************************************************************

  ! The total number of values in the control vector:
  !
  integer function fu_size_control(control) 
    type(da_control), intent(in) :: control
    fu_size_control = control%nvals
  end function fu_size_control

  !*********************************************************
  
  integer function fu_n_initial(control) result(n)
  ! The dimensions of the sub-components:
  !
    type(da_control), intent(in) :: control
  
    if (.not. fu_have_initial(fu_mode(control))) then
      call set_error('Bad DA mode', 'fu_n_initial')
      return
    end if
    n = control%size_initial
    
  end function fu_n_initial

  !**********************************************************
  
  integer function fu_n_emission(control) result(n)
    type(da_control), intent(in) :: control
    
    if (.not. (fu_have_emission_xy(fu_mode(control)) .or. &
      & fu_have_emission_volc(fu_mode(control)) .or. &
      & fu_have_emission_zt(fu_mode(control)))) then
      call set_error('Bad DA mode', 'fu_n_emission')
      return
    end if
    n = control%size_emission

  end function fu_n_emission

  !************************************************************************************
  
  function fu_species_initial(control) result(p_species)
    implicit none
    type(da_control), target, intent(in) :: control
    type(silam_species), dimension(:), pointer :: p_species
    p_species => null()
    if (control%nsp_init > 0) p_species => control%species_initial(1:control%nsp_init)
  end function fu_species_initial

  function fu_species_emission_ctrl(control) result(p_species)
    implicit none
    type(da_control), target, intent(in) :: control
    type(silam_species), dimension(:), pointer :: p_species
    p_species => null()
    if (control%nsp_emis >0) p_species => control%species_emission(1:control%nsp_emis)
  end function fu_species_emission_ctrl

  function fu_isp_initial(control) result(iptr)
    implicit none
    type(da_control), target, intent(in) :: control
    integer, dimension(:), pointer :: iptr
    iptr => null()
    if (control%nsp_init > 0) iptr => control%ind_species_init(1:control%nsp_init)
  end function fu_isp_initial

  function fu_isp_emis_ctrl(control) result(iptr)
    implicit none
    type(da_control), target, intent(in) :: control
    integer, dimension(:), pointer :: iptr
    iptr => null()
    if (control%nsp_emis > 0) iptr =>  control%ind_species_emis(1:control%nsp_emis)
  end function fu_isp_emis_ctrl


  integer function fu_nsp_initial(control)
    implicit none
    type(da_control), target, intent(in) :: control
    fu_nsp_initial = control%nsp_init
  end function fu_nsp_initial

  integer function fu_nsp_emis_ctrl(control)
    implicit none
    type(da_control), target, intent(in) :: control
    fu_nsp_emis_ctrl = control%nsp_emis
  end function fu_nsp_emis_ctrl

  subroutine set_bgr_cov_ptr(control, bgr_cov_ptr)
    implicit none
    type(da_control), intent(inout) :: control
    type(t_background_covariance), target :: bgr_cov_ptr

    control%bgr_cov_ptr => bgr_cov_ptr
  end subroutine set_bgr_cov_ptr

  !************************************************************************************

  integer function fu_mode_control(control) result(mode)
    implicit none
    type(da_control), intent(in) :: control
    mode = control%mode
  end function fu_mode_control

  logical function fu_in_physical_space(control) result(in_physical)
    implicit none
    type(da_control), intent(in) :: control
    in_physical = control%transf_state == physical_space
  end function fu_in_physical_space

  !************************************************************************************
  
  subroutine control_from_files(control, cloud, file_initial, update_initial, default_initial, &
                              & file_emission, update_emission, default_emission, valid_time, &
                              & ifRandomise)
    implicit none
    type(da_control), intent(inout) :: control
    type(silam_pollution_cloud), intent(in) :: cloud
    character(len=*), intent(in) :: file_initial, file_emission
    logical, intent(in) :: update_initial, update_emission
    real, intent(in) :: default_initial, default_emission
    type(silja_time), intent(in) ::valid_time
    logical, intent(in) :: ifRandomise

    type(Tmass_map_ptr), dimension(1) :: mass_map_ptr_array
    type(Tmass_map), pointer :: map_c
    type(Tmass_map), target :: map_cnc
    integer :: ix, iy, iz, isp, nx, ny, nz, nsp, idim, ind, isp_transp, nsp_contr, nFields_updated
    real, dimension(:), pointer :: init_p, emis_p
    character(len=*), parameter :: sub_name = 'control_from_file'

    map_c => fu_concMM_ptr(cloud)
    
    nx = map_c%nx
    ny = map_c%ny
    nz = map_c%n3d
    
    if (iand(fu_mode(control),control_init) /= 0) then 
      init_p => fu_initial_ptr(control)
      !if (rules%have_init_background) then
      if (update_initial) then
        !
        ! Initial state must come from a file. We borrow the
        ! concentration map for reading it.
        nFields_updated = 0
        call set_mass_map(map_cnc, concentration_flag, 1, 0, map_c%gridTemplate, map_c%vertTemplate, &
                             & map_c%species, val=0.0)
        mass_map_ptr_array(1)%ptrMassMap => map_cnc
        call update_mass_map_from_file(file_initial, &
                                     & grads_file_format, &
                                     & mass_map_ptr_array,&
                                     & 1, &
                                     & valid_time, ifRandomise, &
                                     & nFields_updated=nFields_updated)
        

        if(nFields_updated == 0)call msg_warning('DA_INITIAL_STATE: No fields have been red','control_from_files')
        if (error) return

        nsp_contr = size(control%ind_species_init)
        do isp = 1, nsp_contr
          call msg('isp_contr, real(isp_transp)', isp, real(control%ind_species_init(isp)))
        end do
        do iy = 1, ny
          do ix = 1, nx
            do iz = 1, nz
              do isp = 1, nsp_contr
                isp_transp = control%ind_species_init(isp)
                ind = (iy-1)*nx*nz*nsp_contr + (ix-1)*nz*nsp_contr + (iz-1)*nsp_contr + isp
                init_p(ind) = map_cnc%arm(isp_transp, 1, iz, ix, iy)
              end do
            end do
          end do
        end do
        call dealloc_mass_map(map_cnc)
      else
        init_p = default_initial
      end if
    endif

      !call msg('sum of background init_p', sum(init_p))
      !call msg('sum of arm', sum(map_cnc%arm))
    if (iand(fu_mode(control),control_emis_xy) /= 0) then 
      emis_p => fu_emission_ptr(control)
      if (update_emission) then
        call emission_from_grads(file_emission, fu_species_emission(cloud), emis_p)
      else
        emis_p = default_emission
      end if
    endif  
      

    if (iand(fu_mode(control),control_emis_zt) /= 0) then
      emis_p => fu_emission_ptr(control)
      if (update_emission) then
        call set_error('Cannot update emission_time_height', sub_name)
      else
        emis_p = default_emission
      end if
     endif

    
  contains

    subroutine emission_from_grads(file_name, emission_species, p_emis)
      implicit none
      character(len=*), intent(in) :: file_name
      type(silam_species), dimension(:), intent(in) :: emission_species
      real, dimension(:), intent(out) :: p_emis

      real, dimension(:), pointer :: input_grid_data
      integer :: isp, nsp_emis, ind_gf, ind_start, ind_end, step
      integer :: nx, ny, offx, offy, gnx, gny
      type(silja_field_id) :: fid
      
      ind_gf = fu_open_gradsfile_i(file_name)
      if (error) return

      if (.not. (fu_silamGrid_of_grads(ind_gf) == wholeMPIdispersion_grid)) then
        call msg('Grid for emission correction background field:')
        call report(fu_silamGrid_of_grads(ind_gf))
        call msg('Dispersion grid:')
        call report(dispersion_grid)
        call msg('wholeMPIdispersion grid:')
        call report(wholeMPIdispersion_grid)
        call set_error('Bad emission correction background grid', 'emission_from_grads')
        return
      end if
      !if (.not. fu_cmp_verts_eq(fu_silamVert_of_grads(ind_gf), dispersion_vertical)) then
      !  call msg('Vertical for emission correction background field:')
      !  call report(fu_silamVert_of_grads(ind_gf))
      !  call msg('Dispersion vertical:')
      !  call report(dispersion_vertical)
      !  call set_error('Bad emission correction background vertical', 'emission_from_grads')
      !  return
      !end if
      call smpi_get_decomposition(nx, ny, offx, offy, gnx, gny)
      if (nx /= gnx .or. ny /= gny) then
        call set_error("MPI emission_from_grads not implemented yet", "emission_from_grads")
        return
      endif

      input_grid_data => fu_work_array(gnx*gny)

      nsp_emis = size(emission_species)
      do isp = 1, nsp_emis
        fid = fu_set_field_id(met_src_missing, emission_scaling_flag, valid_time, zero_interval, &
                            & wholeMPIdispersion_grid, surface_level, &
                            & species=emission_species(isp))
        if (error) return
        call read_field_from_grads_id(ind_gf, fid, input_grid_data)
        if (error) return
        ind_start = isp
        step = nsp_emis
        ind_end = ind_start + (fs_dispersion-1)*step
          p_emis(ind_start:ind_end:step) = max(input_grid_data(1:fs_dispersion), 0.0)
      end do
      
      call close_gradsfile_i(ind_gf)
      call free_work_array(input_grid_data)
      
    end subroutine emission_from_grads

    !=================================================================
    
    subroutine test_time_hgt_from_file()
      implicit none
      integer, parameter :: num_slots = 4, num_levs = 2
      type(silja_time), dimension(2, num_slots) :: time_slots
      real, dimension(:), pointer :: p_emis
      type(silam_species), dimension(1) :: emission_species
      type(silam_vertical) :: vertical
      character(len=*), parameter :: file_name = 'tst', sub_name = 'test_time_hgt_from_file'
      integer :: stat, unit, ind_slot
      real, dimension(num_levs, num_slots) :: values

      time_slots = reshape((/fu_set_time_utc(2010,1,1,0,0, 0.0), fu_set_time_utc(2010,1,1,12,0,0.0), &
                          &  fu_set_time_utc(2010,1,1,12,0, 0.0), fu_set_time_utc(2010,1,2,0,0,0.0), &
                          &  fu_set_time_utc(2010,1,2,0,0, 0.0), fu_set_time_utc(2010,1,2,12,0,0.0), &
                          &  fu_set_time_utc(2010,1,2,12,0, 0.0), fu_set_time_utc(2010,1,3,0,0,0.0)/), &
                          & (/2, num_slots/))
      p_emis => fu_work_array()
      call set_species(emission_species(1), fu_get_material_ptr('SO2'), in_gas_phase)
      call set_vertical((/fu_set_level(constant_height, fval1=20.0), &
                       &  fu_set_level(constant_height, fval1=50.0)/), vertical)
      if (error) return
      
      unit = fu_next_free_unit()
      open(unit, file=file_name, action='write', form='formatted', iostat=stat)
      if (fu_fails(stat == 0, 'Cannot open test file', sub_name)) return
      call random_number(values)
      do ind_slot = 1, num_slots
        write(unit, fmt='(A, 1x, A, 2G16.6)') fu_time_to_io_string(time_slots(1, ind_slot)), &
             & 'SO2_gas', values(1, ind_slot), values(2, ind_slot)
      end do
      close(unit)
      
      call time_hgt_from_file(file_name, emission_species, time_slots, vertical, p_emis)

      if (.not. all(p_emis(1:num_levs*num_slots) .eps. reshape(values, (/num_levs*num_slots/)))) then
        call set_error('Values are not same', sub_name)
      end if

      call free_work_array(p_emis)
      
    end subroutine test_time_hgt_from_file
    
    !==============================================================
    
    subroutine time_hgt_from_file(file_name, emission_species, time_slots, vertical, p_emis)
      implicit none
      character(len=*), intent(in) :: file_name
      type(silam_species), dimension(:), intent(in) :: emission_species
      type(silja_time), intent(in), dimension(:,:) :: time_slots
      type(silam_vertical), intent(in) :: vertical
      real, dimension(:), intent(out) :: p_emis

      integer :: file_unit, stat
      character(len=fnlen) :: file_name_pr
      character(len=*), parameter :: sub_name = 'time_hgt_from_file'
      logical, dimension(size(time_slots, 2), size(emission_species)) :: slot_ok
      integer :: year, month, day, hour, minute, num_words, quantity, ind_slot, ind_species, ilev, &
           & ind_time_height, num_levs 
      real :: sec
      real, dimension(:), pointer :: values
      character(len=16), dimension(:), allocatable :: words
      type(silam_species) :: species
      logical :: eof
      character(len=worksize_string) :: species_name, line
      type(silja_time) :: when
      character(len=10) :: utc
      integer :: num_words_expect 

      file_name_pr = fu_process_filepath(file_name, must_exist=.true.)
      if (error) return
      file_unit = fu_next_free_unit()

      num_levs = fu_NbrOfLevels(vertical)

      num_words_expect = num_levs + 8
      allocate(words(num_words_expect), stat=stat) ! 7 words from time + species name
      if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
      call msg('Reading time-height background from ' // trim(file_name_pr))
      open(file_unit, file=trim(file_name_pr), form='formatted', action='read', iostat=stat)
      if (fu_fails(stat == 0, 'Error opening time height background file', sub_name)) return

      slot_ok = .false.
      values => fu_work_array()
      do
        call next_line_from_input_file(file_unit, line, eof)
        if (error) return
        if (eof) exit
        read(line, fmt=*, iostat=stat) year, month, day, hour, minute, sec, utc
        if (fu_fails(stat == 0, 'Failed to read times', sub_name)) return
        if (fu_fails(utc == 'UTC', 'utc is not UTC', sub_name)) return
        when = fu_set_time_utc(year, month, day, hour, minute, sec)
        if (error) return
        ind_slot = fu_index(when, time_slots(1,:))
        if (fu_fails(ind_slot /= int_missing, 'Failed to find time slot', sub_name)) return
        call split_string(line, ' ', words, num_words)
        
        if (fu_fails(num_words == num_words_expect, &
                   & 'Bad number of words in time-height background field', sub_name)) return
        
        if (error) return
        species_name = words(8)
        ! add cnc_ and use the grads-io subroutine to find species.
        call decode_id_params_from_io_str('cnc_' // species_name, num_levs > 1, quantity, species, .true.)
        if (error) return
        ind_species = fu_index(species, emission_species)
        if (fu_fails(ind_species /= int_missing, 'Species not found in emission', sub_name)) return
        do ilev = 1, num_levs
          ind_time_height = fu_time_height_index(ilev, ind_slot, ind_species, num_levs, &
                                               & size(emission_species))
          read(unit=words(8+ilev), fmt=*, iostat=stat) p_emis(ind_time_height)
          if (fu_fails(stat == 0, 'Failed to parse value', sub_name)) return
          !p_emis(ind_time_height) = 
        end do
        slot_ok(ind_slot, ind_species) = .true.
      end do
      
      call free_work_array(values)
      
      if (fu_fails(all(slot_ok), 'Not all timeslots OK', sub_name)) then
        do ind_slot = 1, size(slot_ok, 1)
          do ind_species = 1, size(slot_ok, 2)
            call report(emission_species(ind_species))
            call msg(fu_str(time_slots(1,ind_slot)))
            if (slot_ok(ind_slot, ind_species)) then
              call msg('--> Slot OK')
            else
              call msg('--> Slot not OK')
            end if
          end do
        end do
      end if
    end subroutine time_hgt_from_file
  end subroutine control_from_files

  !*****************************************************************************
  
  integer function fu_time_height_index(ind_lev, ind_slot, ind_species, num_levs, num_species) result(ind)
    implicit none
    integer, intent(in) :: ind_lev, ind_slot, ind_species, num_levs, num_species
    ind = (ind_slot-1)*num_species*num_levs + (ind_lev-1)*num_species + ind_species
  end function fu_time_height_index

  !************************************************************************************
  
  subroutine set_da_rules(nlptr, nlStandardSetup, darules, startTime, periodToCompute)
    ! 
    ! Set DA rules from the namelist. Exception: 3dvar rules needs chemicals to be set.
    !
    implicit none
    type(Tsilam_namelist), pointer :: nlPtr, nlStandardSetup
    type(da_rules), intent(inout), target :: darules
    type(silja_time), intent(in) :: startTime
    type(silja_interval), intent(in) :: periodToCompute

    !type(DA_rules), intent(out) :: rules
    type(da_rules), pointer :: rules
    character(len=worksize_string) :: nl_content
    character(len=*), parameter :: sub_name = 'set_da_rules'
    integer :: ind_obs, num_obs_items, stat
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: items
    logical :: need_emission, need_initial
    integer :: ind_subst, num_subst, num_items, ind_step, int_content

    rules => darules
    nullify(items)
    
    ! Method
    nl_content = fu_str_u_case(fu_content(nlptr, 'method'))
    select case(nl_content)
    case ('3D')
      rules%method = flag_3dvar
    case ('4D')
      rules%method = flag_4dvar
    case ('4D_SEQ')
      rules%method = flag_4dvar_seq
    case ('MATRIX')
      rules%method = flag_h_matrix
    case ('ENKF')
      rules%method = flag_enkf
    case ('ENKS')
      rules%method = flag_enks
    case default
      call set_error('Invalid method', sub_name)
      return
    case ('')
      call set_error('Missing method', sub_name)
      return
    end select

    rules%allow_negative = rules%method == flag_3dvar


    ! Parse parameters without any analysis
    nl_content = fu_content(nlptr, 'assimilation_window')
    rules%da_window = interval_missing
    if (nl_content /= '') rules%da_window = fu_set_named_interval(nl_content)
    if (error) then
      call set_error("assimilation_window",sub_name)
      return
    endif
    nl_content = fu_content(nlptr, 'assimilation_interval')
    rules%assim_interval = interval_missing
    if (nl_content /= '') rules%assim_interval = fu_set_named_interval(nl_content)
    if (error) then
      call set_error("assimilation_interval",sub_name)
      return
    endif
    nl_content = fu_content(nlptr, 'smoother_step')
    rules%smoother_step = interval_missing
    if (nl_content /= '') rules%smoother_step = fu_set_named_interval(nl_content)
    if (error) then
      call set_error("smoother_step",sub_name)
      return
    endif

    ! These two are used only for ensemble  DAs so far
    nl_content = fu_content(nlptr, 'first_analysis_time')
    rules%first_analysis_time = time_missing
    if (nl_content /= '') rules%first_analysis_time = fu_io_string_to_time(nl_content)
    if (error) then
      call set_error("first_analysis_time",sub_name)
      return
    endif
    nl_content = fu_content(nlptr, 'last_analysis_time')
    rules%last_analysis_time = time_missing
    if (nl_content /= '') rules%last_analysis_time = fu_io_string_to_time(nl_content)
    if (error) then
      call set_error("last_analysis_time",sub_name)
      return
    endif



    ! ini
    ! dump                   <---+
    ! assimilate                 |  
    ! run                    ----+


    
    ! Define the assimilation window. For 4D, the
    ! assimilation window is just taken from the time window of the run.
    select case(rules%method)
      case(flag_4dvar,flag_h_matrix) 
        rules%da_begin = startTime
        if (defined(rules%da_window)) then
          call set_error("da_window should not appear in the DA namelist for 4dvar",sub_name)
          return
        else
          rules%da_window = periodToCompute
        endif
        if (.not. fu_sec(rules%da_window) > 0) then
            call set_error('Bad assimilation window: '//fu_str(rules%da_window), sub_name)
            return
        endif

      case (flag_3dvar)
          
       if (defined(rules%da_window)) then
          call msg_warning("Defined da_window. Please use assimilation_interval instead", sub_name)
          if (.not. defined(rules%assim_interval)) then !! workaround for obsolete syntax
             rules%assim_interval = rules%da_window
          else
             call set_error("da_window should not appear for 3dvar",sub_name)
             return
          endif
       endif
       rules%da_window = zero_interval
       if (.not. fu_sec(rules%assim_interval) > 0) then
           call set_error('Bad assimilation interval: '//fu_str(rules%assim_interval), sub_name)
           return
       endif
       if (.not. defined(rules%first_analysis_time)) rules%first_analysis_time = startTime
       if (.not. defined(rules%last_analysis_time) ) &
                                  &rules%last_analysis_time =  startTime + periodToCompute



      case (flag_4dvar_seq)
        call set_error("flag_4dvar_seq not implemented (YET?)", sub_name)
        return
      case (flag_enkf,flag_enks)

        rules%da_begin = startTime ! will be later set according to the current step
        

        nl_content = fu_content(nlptr, 'smoother_step')
        if (rules%method == flag_enks .and. .not. defined(rules%smoother_step) ) then
            call msg('No smoother_step -- only single-time analysis')
        end if
        
        call get_items(nlptr, 'smoother_output_step', items, num_items)
        if (num_items > 0) then
          if (.not. defined(rules%smoother_step)) then
            call set_error('smoother_output_step requested, no smoother_step', sub_name)
            return
          end if
          allocate(rules%smoother_output_steps(num_items), stat=stat)
          if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
          do ind_step = 1, num_items
            int_content = fu_content_int(items(ind_step))
            if (fu_fails(int_content >= 0, 'Bad smoother_output_step', sub_name)) return
            !if (rules%smoother_step*int_content > rules%assim_interval) then
            !  call set_error('Requested smoother_output_step exceeds available', sub_name)
            !  return
            !end if
            rules%smoother_output_steps(ind_step) = int_content
          end do
          rules%num_smoother_output_steps = num_items
        else
          nullify(rules%smoother_output_steps)
          rules%num_smoother_output_steps = 0
        end if

        if (.not. defined(rules%last_analysis_time)) then
          ! No assimilation at simulation end step; this would require extending beyond
          ! simulation end if 4d-var is used. 
          rules%last_analysis_time = startTime + periodToCompute - rules%assim_interval
        end if

    case default
      call set_error('Strange assimilation method', sub_name)
      return
    end select

    nl_content = fu_process_filepath(fu_content(nlPtr,'output_directory'))
    if (fu_fails(nl_content /= '', 'Missing output_directory', sub_name)) return
    !call get_main_output_dir(simrules, main_output_dir)
    !call replace_string(nl_content, '%main_output_dir', main_output_dir, case_sensitive=.false.)

    call msg('D/A output to ' // trim(nl_content))
    rules%outdir_templ_str = trim(nl_content)
    !
    ! Get rid of obsolete statement and enforce extended list format
    !
    if(len_trim(fu_content(nlptr, 'output_level')) > 0)then
      call set_error('output_level is obsolete. Use da_output = none / first_last / da_trajectory / full_dump',sub_name)
      return
    endif
    nl_content = fu_content(nlptr,'da_output')
    select case(nl_content)
      case ('full_dump')
        rules%output_level = da_out_full_dump_flag
      case ('da_trajectory')
        rules%output_level = da_out_trajectory_flag
      case ('first_last')
        rules%output_level = da_out_first_last_flag
      case ('none')
        rules%output_level = da_out_none_flag
      case default
        call set_error('Bad output_level:' // trim(nl_content), sub_name)
        return
    end select
    !
    ! The reprojection aliasing may be smoothed by randomization
    !
    rules%ifRandomise = .not. fu_str_u_case(fu_content(nlStandardSetup,'randomise_reprojection')) == 'NO'

    ! observation templates
    !
    if (error) return
    call get_items(nlptr, 'obs_data', items, num_obs_items)
    rules%num_obs_items = num_obs_items
    if (num_obs_items > 0) then
      allocate(rules%obs_items(num_obs_items), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
      do ind_obs = 1, num_obs_items
        rules%obs_items(ind_obs) = fu_content(items(ind_obs))
      end do
    end if
    rules%obs_list_path = fu_content(nlptr, 'obs_list') ! may be empty
    rules%station_path = fu_content(nlptr, 'station_list') ! may be empty if no surface obs


    ! assimilated species
    !
    nullify(items)
    call get_items(nlptr, 'analysis_substance', items, num_subst)
    if (error) return
    !if (num_subst == 0 .or. rules%method == flag_4dvar) then
    if (num_subst == 0) then
      nullify(rules%analysis_subst_list_3d)
      rules%have_analysis_species = silja_false
    else 
      allocate(rules%analysis_subst_list_3d(num_subst), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
      do ind_subst = 1, num_subst
        rules%analysis_subst_list_3d(ind_subst) = fu_content(items(ind_subst))
      end do
      rules%have_analysis_species = silja_true
    end if


    !************************************************************************************
    ! call the specific setup for enkf / variational

    if (rules%method == flag_enkf .or. rules%method == flag_enks) then
      call set_enkf_rules()
    else
      call set_xdvar_rules()
    end if

    rules%defined = .true.

  contains

    subroutine set_enkf_rules()
      implicit none
      integer :: nl_content_int
      character(len=*), parameter :: sub_name = 'set_enkf_rules'
      type(Tsilam_nl_item_ptr), dimension(:), pointer :: p_items
      integer :: num_pert_vars, ind_var, request, num_state_vars
      real :: nl_content_real
      character(len=clen) :: unit

      !rules%controlVariable = da_initial_state
      nl_content_int = fu_content_int(nlptr, 'ensemble_size')
      if (nl_content_int == int_missing) then
        call msg('Ensemble size from mpi size:', smpi_ens_size)
        rules%ens_size = smpi_ens_size
      else
        if (fu_fails(nl_content_int > 0, 'Bad ensemble_size', sub_name)) return
        if (fu_fails(nl_content_int == smpi_ens_size, 'Ensemble size not matching mpi', sub_name)) return
        rules%ens_size = nl_content_int
      end if
      !
      ! Perturbed variables
      !
      nullify(p_items)
      call get_items(nlptr, 'perturb_variable', p_items, num_pert_vars)
      rules%perturbVariable = control_none
      do ind_var = 1, num_pert_vars
        request = fu_get_control_bits(fu_content(p_items(ind_var)))
        rules%perturbVariable = ior(rules%perturbVariable, request)
      end do
      if (rules%perturbVariable == 0) then
        call msg('No forced perturbations requested')
      else
        call msg('perturbvar:', rules%perturbVariable)
      end if
      if (num_pert_vars > 0) deallocate(p_items)
      !
      ! State variables
      !
      call get_items(nlptr, 'state_variable', p_items, num_state_vars)
      if (num_state_vars == 0) then
        ! default
        rules%controlVariable = da_initial_state
      else
        rules%controlVariable = control_none
        do ind_var = 1, num_state_vars
          request = fu_get_control_bits(fu_content(p_items(ind_var)))
          rules%controlVariable = ior(rules%controlVariable, request)
        end do

        if (fu_have_emission_zt(rules%controlVariable)) then 
          rules%emis_time_slot = fu_set_named_interval(fu_content(nlptr, 'emis_time_slot'))
          if (fu_fails(.not. error, 'Failed to parse emis_time_slot', sub_name)) return
        end if
        if (fu_have_emission_volc(rules%controlVariable) .or. fu_have_emission_zt(rules%controlVariable)) then 
          nl_content_real = fu_content_real(nlptr, 'zt_loc_lon')
          if (fu_fails(.not. (nl_content_real .eps. real_missing), 'Missing zt_loc_lon', sub_name)) return
          rules%zt_loc_lon = nl_content_real
          nl_content_real = fu_content_real(nlptr, 'zt_loc_lat')
          if (fu_fails(.not. (nl_content_real .eps. real_missing), 'Missing zt_loc_lat', sub_name)) return
          rules%zt_loc_lat = nl_content_real
        end if

        if (rules%controlVariable == 0) then
          call set_error('No control variable', sub_name)
          return
        end if
        deallocate(p_items)
      end if  !num_state_vars > 0

      if (fu_have_emission_volc(rules%controlVariable) .and. &
                         & .not. fu_have_emission_volc(rules%perturbVariable)) then
        call set_error('Volc emission not included in perturbVariable', sub_name)
        return
      end if

      if (fu_have_emission_volc(rules%perturbVariable)) then
        ! set some ze-emission specific rules
        rules%time_height_mode = processor_tm_hgt_force
        rules%emis_time_slot = rules%da_window
      else if (fu_have_emission_zt(rules%perturbVariable)) then
        rules%time_height_mode = processor_tm_hgt_scale
      else if (fu_have_emission_xy(rules%perturbVariable)) then
        rules%emis_time_slot = rules%da_window
      end if

      ! Background/perturbation covariances
      !
      need_initial = fu_have_initial(rules%perturbVariable)
      if (need_initial) then
        nl_content = fu_content(nlptr, 'cov_setup_initial')
        if (fu_fails(nl_content /= '', 'Missing cov_setup_initial', sub_name)) return
        call decode_template_string(nl_content, rules%cov_setup_templ_init)
        if (error) return
      end if
      !******************************************************************
      !need_emission =  fu_have_emission_xy(rules%perturbVariable) &
      !     & .or. fu_have_emission_zt(rules%perturbVariable)
      need_emission = .false.
      !******************************************************************
      if (need_emission) then
        nl_content = fu_content(nlptr, 'cov_setup_emission')
        if (fu_fails(nl_content /= '', 'Missing cov_setup_emission', sub_name)) return
        call decode_template_string(nl_content, rules%cov_setup_templ_emis)
        if (error) return        
      end if
      
      call get_real_param(nlptr, 'observation_stdev', rules%observation_stdev, .false., .false., sub_name)
      call get_real_param(nlptr, 'emission_stdev', rules%emission_stdev, .true., .false., sub_name)
      call get_real_param(nlptr, 'emis_corr_dist', rules%emis_corr_dist, .true., .false., sub_name)
      call get_real_param(nlptr, 'emis_max_corr', rules%emis_max_corr, .true., .false., sub_name)

      ! volcanic emission parametric perturbations
      !
      if (iand(control_emis_volc, rules%perturbVariable) /= 0) then
        call get_real_param(nlptr, 'volc_massflux_scaling', rules%volc_massflux_scaling, &
                          & if_positive=.true., if_required=.true., caller=sub_name)
        call get_real_param(nlptr, 'volc_volflux_sigma', rules%volc_volflux_sigma, .true., .true., sub_name)
        call get_real_param(nlptr, 'volc_volflux_mean_base_10_exponent', rules%volc_volflux_mode, .true., .true., sub_name)
        call get_real_param(nlptr, 'volc_height_sigma', rules%volc_height_sigma, .true., .false., sub_name)
        call get_real_param(nlptr, 'volc_height_exp', rules%volc_height_exp, .false., .false., sub_name)
        call get_real_param(nlptr, 'volc_height_scaling', rules%volc_height_scaling, &
                          & .true., .false., sub_name)
        call get_real_param(nlptr, 'volc_fract_mass_top', rules%volc_fract_mass_top, &
                          & .true., .false., sub_name)
        call get_real_param(nlptr, 'volc_fract_height_top', rules%volc_fract_height_top, &
                          & .true., .false., sub_name)
        call get_real_param(nlptr, 'volc_fract_height_sigma', rules%volc_fract_height_sigma, &
                          & .true., .false., sub_name)
        call get_real_param(nlptr, 'volc_', rules%volc_height_exp, .false., .false., sub_name)
        if (error) return
        
        nl_content = fu_str_l_case(fu_content(nlptr, 'use_log_cloud'))
        rules%use_log_cloud = nl_content == 'yes'

        nl_content = fu_str_l_case(fu_content(nlptr, 'use_log_obs'))
        rules%use_log_obs = nl_content == 'yes'

        nl_content = fu_content(nlptr, 'volc_volflux_corr_time')
        if (fu_fails(nl_content /= '', 'missing volc_volflux_corr_time', sub_name)) return
        rules%volc_volflux_corr_time = fu_sec(fu_set_named_interval(nl_content))
        if (error) return
        if (fu_fails(rules%volc_volflux_corr_time > 0, 'Correlation time < 0', sub_name)) return

      else if (fu_have_emission_xy(rules%perturbVariable)) then ! correlation time for xy emission
        nl_content = fu_content(nlptr, 'emis_corr_time')
        if (fu_fails(nl_content /= '', 'missing emis_corr_time', sub_name)) return
        rules%emis_corr_time = fu_sec(fu_set_named_interval(nl_content))
        if (error) return
        if (fu_fails(rules%emis_corr_time > 0, 'Emission correlation time < 0', sub_name)) return

      end if ! have volcanic emission perturbations

      ! Filter parameters
      nl_content = fu_content(nlptr, 'localisation_distance')
      if (fu_fails(nl_content /= '', 'Missing localisation_distance', sub_name)) return
      call set_named_value_and_unit(nl_content, nl_content_real, unit)
      rules%loc_dist_m = nl_content_real * fu_conversion_factor(unit, 'm')
      if (error) return
      nl_content = fu_content(nlptr, 'enkf_flavor')
      select case(nl_content)
      case ('enkf')
        rules%enkf_flavor = enkf_enkf
      case ('denkf')
        rules%enkf_flavor = enkf_denkf
      case default
        call set_error('enkf_flavor must enkf or denkf', sub_name)
        return
      end select
      nl_content = fu_content(nlptr, 'loc_type') 
      select case(nl_content)
      case ('step')
        rules%loc_type = loc_step
      case ('none')
        rules%loc_type = loc_none
      case ('gaspari_cohn')
        rules%loc_type = loc_gaspari_cohn
      case default
        call msg('Allowed values for loc_type: none step gaspari_cohn')
        call set_error('Bad loc_type', sub_name)
        return
      end select
      nl_content_real = fu_content_real(nlptr, 'rfactor')
      if (nl_content_real .eps. real_missing) then
        call set_error('Missing rfactor', sub_name)
        return
      end if
      rules%rfactor = nl_content_real

    end subroutine set_enkf_rules

    !===========================================================
    
    subroutine get_real_param(nlptr, param_name, param, if_positive, if_required, caller)
      implicit none
      type(Tsilam_namelist), pointer :: nlptr
      character(len=*), intent(in) :: param_name, caller
      real, intent(inout) :: param
      logical, intent(in) :: if_positive, if_required
      
      real :: nl_content_real
      
      nl_content_real = fu_content_real(nlptr, param_name)
      if (nl_content_real .eps. real_missing) then
        if (fu_fails(.not. if_required, 'Missing '//trim(param_name), caller)) return
      else
        if (if_positive .and. nl_content_real < 0) then
          call set_error('Bad '//trim(param_name), caller)
          return
        end if
        param = nl_content_real
      end if

    end subroutine get_real_param

    !=============================================================
    
    subroutine set_xdvar_rules()
      implicit none

      character(len=*), parameter :: sub_name = 'set_xdvar_rules'
      real :: nl_content_real
      integer :: nl_content_int
    
      call set_numerics(nlptr, rules%numerics)
      if (error) return
      
      ! Control variable and the background variances
      ! 
      nl_content = fu_content(nlptr, 'control_variable')
      select case(nl_content)
      case ('initial_state')
        rules%controlVariable = da_initial_state
      case ('emission_correction')
        rules%controlVariable = da_emission_correction
      case ('emission_and_initial')
        rules%controlVariable = da_emission_and_initial
      case ('emission_time_height')
        rules%controlVariable = da_emission_time_height
        if (fu_fails(fu_content(nlptr, 'emis_time_slot') /= '', 'Missing emis_time_slot', sub_name)) return
        rules%emis_time_slot = fu_set_named_interval(fu_content(nlptr, 'emis_time_slot'))
        if (error) then
          call set_error('Failed to parse emis_time_slot', sub_name)
          return
        end if
        select case(fu_str_l_case(fu_content(nlptr, 'time_height_mode')))
        case ('force')
          rules%time_height_mode = processor_tm_hgt_force
        case ('scale', '')
          rules%time_height_mode = processor_tm_hgt_scale
        case ('force_weighted')
          rules%time_height_mode = processor_tm_hgt_force_wgt
        case default
          call set_error('Bad time_height_mode', sub_name)
          return
        end select

      case ('')
        call set_error('Missing control_variable:', sub_name)
      case default
        call set_error('Strange control_variable:' // trim(nl_content), sub_name)
      end select

      if (error) return

      if (rules%method == flag_3dvar .and. rules%controlVariable /= da_initial_state) then
        call set_error('Bad control variable for 3dvar', 'set_xdvar_rules')
        return
      end if

      ! Background fields
      !
      nl_content = fu_content(nlPtr,'initial_state_background_file')
      if (nl_content /= '') then
        nl_content = fu_process_filepath(nl_content, must_exist=.true.)
        rules%have_init_background = .true.
        if (error) return
      else
        call msg_warning('No background given for initial state')
        rules%have_init_background = .false.
      end if
      rules%initialStateBgrFile = nl_content

      nl_content = fu_content(nlptr, 'emission_background_file')
      if (nl_content /= '') then
        nl_content = fu_process_filepath(nl_content, must_exist=.true.)
        rules%have_emis_background = .true.
        if (error) return
      else
        call msg_warning('No background given for emission correction')
        rules%have_emis_background = .false.
      end if
      rules%emissionCorrectionBgrFile = nl_content

      nl_content = fu_str_l_case(fu_content(nlptr, 'disable_background'))
      rules%no_background_term = nl_content == 'yes'
      if (rules%no_background_term) call msg_warning('Background term disabled in cost function!')

      nl_content = fu_str_l_case(fu_content(nlptr, 'zero_emis_backgr'))
      rules%use_zero_emis_backgr = nl_content == 'yes'

      ! Background/perturbation covariances
      !
      need_initial = fu_have_initial(rules%controlVariable)
      if (need_initial) then
        nl_content = fu_content(nlptr, 'cov_setup_initial')
        if (fu_fails(nl_content /= '', 'Missing cov_setup_initial', sub_name)) return
        call decode_template_string(nl_content, rules%cov_setup_templ_init)
        if (error) return
      end if
      need_emission = fu_have_emission_xy(rules%controlVariable) &
           & .or. fu_have_emission_zt(rules%controlVariable)
      if (need_emission) then
        nl_content = fu_content(nlptr, 'cov_setup_emission')
        if (fu_fails(nl_content /= '', 'Missing cov_setup_emission', sub_name)) return
        call decode_template_string(nl_content, rules%cov_setup_templ_emis)
        if (error) return        
      end if
    
      nl_content = fu_str_l_case(fu_content(nlptr, 'restart_iteration'))
      rules%restart_iteration = nl_content == 'yes'

    end subroutine set_xdvar_rules

    !=======================================================
    
    subroutine set_numerics(nlptr, numerics)
      implicit none
      type(Tsilam_namelist), pointer :: nlptr
      type(DA_numericalParameters), intent(out) :: numerics

      character(len=fnlen) :: nl_content
      real :: nl_content_real
      character(len=*), parameter :: sub_name = 'set_nemerics'
      ! Minimization parameters. Defaults are in the type definition.
      !
      nl_content_real = fu_content_int(nlPtr,'max_iterations')
      if (nl_content_real /= int_missing) numerics%maxIterations = nl_content_real

      nl_content_real = fu_content_real(nlPtr, 'cost_rel_tol')
      if (.not. (nl_content_real .eps. real_missing)) numerics%cost_rel_tol = nl_content_real

      nl_content_real = fu_content_real(nlPtr, 'minimum_step')
      if (.not. (nl_content_real .eps. real_missing)) numerics%minimumStepLength = nl_content_real
      
      !nl_content_real = fu_content_real(nlPtr, 'min_move')
      !if (.not. (nl_content_real .eps. real_missing)) numerics%minmax = nl_content_real

      nl_content_real = fu_content_real(nlPtr, 'grad_rel_tol')
      if (.not. (nl_content_real .eps. real_missing)) numerics%grad_rel_tol = nl_content_real
      
      nl_content_real = fu_content_real(nlPtr, 'quasi_newton_df1')
      if (.not. (nl_content_real .eps. real_missing)) numerics%quasi_newton_df1 = nl_content_real

      nl_content_real = fu_content_real(nlPtr, 'cost_incr_rel_tol')
      if (.not. (nl_content_real .eps. real_missing)) numerics%cost_incr_rel_tol = nl_content_real

      nl_content_real = fu_content_real(nlPtr, 'cost_rel_tol')
      if (.not. (nl_content_real .eps. real_missing)) numerics%cost_rel_tol = nl_content_real

      nl_content = fu_content(nlptr, 'search_method')
      select case(nl_content)
      case ('steepest_descent')
        numerics%searchMethod = steepest_descent_flag
      case ('m1qn3')
        numerics%searchMethod = m1qn3_flag
      case ('l_bfgs_b')
        numerics%searchMethod = l_bfgs_b_flag
      case ('')
        call msg_warning('No search_method, will use default')
      case default
        call set_error('Strange search_method: ' // trim(nl_content), sub_name)
        return
      end select

    end subroutine set_numerics

  end subroutine set_da_rules

  
  !*****************************************************************************
  
  integer function fu_num_emis_time_slots(darules) result(num_slots)
    implicit none
    type(da_rules), intent(in) :: darules

    num_slots = ceiling(darules%da_window / darules%emis_time_slot)
    
  end function fu_num_emis_time_slots

  
  !******************************************************************************
  
  subroutine get_emis_time_slots(darules, ptr_slots, lst_slots)
    implicit none
    type(da_rules), intent(in) :: darules
    ! ptr_slots present => allocate list
    ! lst_slots present => use the given list
    type(silja_time), dimension(:,:), pointer, optional :: ptr_slots
    type(silja_time), dimension(:,:), target, optional :: lst_slots

    integer :: num_slots, ind_slot, stat
    type(silja_interval) :: slot_interval
    type(silja_time), dimension(:,:), pointer :: slots
    logical :: allocate_slots

    allocate_slots = present(ptr_slots)
    if (present(lst_slots)) then
      if (allocate_slots) then
        call set_error('Only one of ptr_slots or lst_slots can be given', 'get_emis_time_slots')
        return
      end if
      slots => lst_slots
    end if
    
    slot_interval = darules%emis_time_slot
    if (fu_fails(.not. (slot_interval == interval_missing), &
               & 'slot_interval not defined', 'get_emis_time_slots')) return
    num_slots = fu_num_emis_time_slots(darules) !ceiling(darules%da_window / slot_interval)
    if (allocate_slots) then
      allocate(ptr_slots(2,num_slots), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', 'get_emis_time_slots')) return
      slots => ptr_slots
    end if
    slots(slot_start, 1) = darules%da_begin !simrules%startTime
    slots(slot_end, 1) = darules%da_begin + slot_interval
    do ind_slot = 2, num_slots
      slots(slot_start, ind_slot) = slots(slot_end, ind_slot-1)
      slots(slot_end, ind_slot) = slots(slot_start, ind_slot) + slot_interval
      if (slots(slot_end, ind_slot) > darules%da_begin + darules%da_window) then
        slots(slot_end, ind_slot) = darules%da_begin + darules%da_window
      end if
      if (fu_fails(.not. slots(slot_end, ind_slot) == slots(slot_start, ind_slot), &
                 & 'Bad time slot', 'get_emis_time_slots')) return
    end do
      
  end subroutine get_emis_time_slots

  !**********************************************************************************

  subroutine to_physical(control, background, control_phys)
    implicit none
    type(da_control), intent(in) :: control
    type(da_control), intent(inout) :: control_phys
    type(da_control), intent(in) :: background
    
    integer :: ii
    !real, dimension(:), pointer :: bgrvals, values
    type(t_background_covariance), pointer :: p_bgr_cov
        
    if (fu_fails(.not.fu_in_physical_space(control), 'Control already in physical space', 'to_physical')) return
    if (fu_size(control_phys) == 0 .and.  fu_size(background) == 0) return

    p_bgr_cov => control%bgr_cov_ptr
    if (fu_fails(associated(p_bgr_cov), 'bgr_cov not associated', 'to_physical')) return

    if (fu_fails(fu_size(control) <= fu_size(control_phys), 'physical space too small', 'to_physical')) return
    if (fu_fails(fu_size(background) == fu_size(control_phys), 'wrong background size', 'to_physical')) return

    select case(fu_mode(control))
    case(DA_INITIAL_STATE)
      call apply_cov(p_bgr_cov%correlation_initial, p_bgr_cov%stdev_initial, &
                   & fu_initial_ptr(control), fu_initial_ptr(background), fu_initial_ptr(control_phys))
      if (error) then
        call set_error('Error with initial state background covariance', 'to_physical')
        return
      end if

    case (DA_EMISSION_CORRECTION)
      call apply_cov(p_bgr_cov%correlation_emission, p_bgr_cov%stdev_emission, &
                   & fu_emission_ptr(control), fu_emission_ptr(background), fu_emission_ptr(control_phys))
      if (error) then
        call set_error('Error with emission correction background covariance', 'to_physical')
        return
      end if
      
    case (DA_EMISSION_AND_INITIAL)
      call apply_cov(p_bgr_cov%correlation_initial, p_bgr_cov%stdev_initial, &
                   & fu_initial_ptr(control), fu_initial_ptr(background), fu_initial_ptr(control_phys))
      if (error) then
        call set_error('Error with initial state background covariance', 'to_physical')
        return
      end if
      call apply_cov(p_bgr_cov%correlation_emission, p_bgr_cov%stdev_emission, &
                   & fu_emission_ptr(control), fu_emission_ptr(background), fu_emission_ptr(control_phys))
      if (error) then
        call set_error('Error with emission correction background covariance', 'to_physical')
        return
      end if

    case (DA_EMISSION_TIME_HEIGHT)
      ! like emission correction in space, but so far no support for correlations:
      if (defined(p_bgr_cov%correlation_emission)) then
        call set_error('Correlation for emission must be void for emission_time_height', &
                     & 'to_physical')
        return
      end if
      call apply_cov(p_bgr_cov%correlation_emission, p_bgr_cov%stdev_emission, &
                   & fu_emission_ptr(control), fu_emission_ptr(background), fu_emission_ptr(control_phys))
      if (error) then
        call set_error('Error with emission correction background covariance', 'to_physical')
        return
      end if

    case default
      call set_error('Unsupported control mode', 'to_physical')
      return
    end select
    
  contains

    subroutine apply_cov(correlation, stdev, controlvals, bgrvals, physvals)
      ! Application of 
      ! x = SL x' + xb 
      ! where S is standard deviations, LL^T is the correlation matrix and xb is background.
      implicit none
      type(t_correlation), intent(in) :: correlation
      real, dimension(:), intent(in) :: stdev
      real, dimension(:) :: physvals, controlvals, bgrvals 
      !real, dimension(:), pointer :: stdev, bgrvals, modelvals
      
      integer :: ii, stat
      real, dimension(:), pointer :: work
      
!      call msg('size(physvals), size(controlvals)', size(physvals), size(controlvals))
!      call msg('size(bgrvals)', size(bgrvals))
      
      if (.not. allocated(transf_work)) then
        allocate(transf_work(size(physvals)), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed', 'to_physical')) return
      else if (size(transf_work) < size(physvals)) then
        deallocate(transf_work)
        allocate(transf_work(size(physvals)), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed', 'to_physical')) return
      end if

      work => transf_work(1:size(controlvals))
      if (defined(correlation)) then
        work(:) = controlvals(:) ! to avoid changing controlvals
        call transf_fwd(correlation, work, physvals)
      else
        physvals = controlvals
      end if
      ! standard deviation + background
      do ii = 1, size(physvals)
        physvals(ii) = stdev(ii)*physvals(ii) + bgrvals(ii)
      end do

    end subroutine apply_cov

  end subroutine to_physical

  !************************************************************************************
  
  subroutine to_control(control_phys, control)
    implicit none
    type(da_control), intent(in) :: control_phys
    type(da_control), intent(inout) :: control
    
    integer :: ii, stat
    real, dimension(:), pointer :: modelvals, values
    type(t_background_covariance), pointer :: p_bgr_cov
        
    p_bgr_cov => control%bgr_cov_ptr
    if (fu_fails(associated(p_bgr_cov), 'bgr_cov not associated', 'to_control')) return

    if (fu_fails(fu_in_physical_space(control_phys), 'Control not in physical space', 'to_control')) return
    if (fu_fails(.not.fu_in_physical_space(control), 'Control not in control space', 'to_control')) return

    if (fu_fails(fu_size(control) <= fu_size(control_phys), 'physical space too small', 'to_control')) return

    select case(fu_mode(control))
    case (DA_INITIAL_STATE)
      call apply_cov(p_bgr_cov%correlation_initial, p_bgr_cov%stdev_initial, &
                   & fu_initial_ptr(control), fu_initial_ptr(control_phys))
      if (error) then
        call set_error('Problem with initial state background covariance', 'to_control')
        return
      end if
      
    case (DA_EMISSION_CORRECTION)
      call apply_cov(p_bgr_cov%correlation_emission, p_bgr_cov%stdev_emission, &
                   & fu_emission_ptr(control), fu_emission_ptr(control_phys))
      if (error) then
        call set_error('Problem with emission background covariance', 'to_control')
        return
      end if

    case (DA_EMISSION_AND_INITIAL)
      call apply_cov(p_bgr_cov%correlation_initial, p_bgr_cov%stdev_initial, &
                   & fu_initial_ptr(control), fu_initial_ptr(control_phys))
      if (error) then
        call set_error('Problem with initial state background covariance', 'to_control')
        return
      end if
      call apply_cov(p_bgr_cov%correlation_emission, p_bgr_cov%stdev_emission, &
                   & fu_emission_ptr(control), fu_emission_ptr(control_phys))
      if (error) then
        call set_error('Problem with emission background covariance', 'to_control')
        return
      end if

    case (DA_EMISSION_TIME_HEIGHT)
      call apply_cov(p_bgr_cov%correlation_emission, p_bgr_cov%stdev_emission, &
                   & fu_emission_ptr(control), fu_emission_ptr(control_phys))
      if (error) then
        call set_error('Problem with emission background covariance', 'to_control')
        return
      end if

    case default
      call set_error('Unsupported control mode', 'to_control')
      return

    end select
    
  contains

    subroutine apply_cov(correlation, stdev, controlvals, physvals)
      ! Application of 
      ! x -> L^T Sx
      ! where S is standard deviations, LL^T is the correlation matrix.
      implicit none
      type(t_correlation) :: correlation
      real, dimension(:), intent(in) :: stdev
      real, dimension(:) :: controlvals, physvals ! no intent to allow pointers (above)

      real, dimension(:), pointer :: work
      integer :: ii

      if (.not. allocated(transf_work)) then
        allocate(transf_work(size(physvals)), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed', 'to_control')) return
      else if (size(transf_work) < size(physvals)) then
        deallocate(transf_work)
        allocate(transf_work(size(physvals)), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed', 'to_control')) return
      end if

      work => transf_work(1:size(physvals))
      !call msg('sum of modelvals before stdev:', sum(physvals))
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)
      !$OMP DO
      do ii = 1, size(physvals)
        work(ii) = stdev(ii)*physvals(ii)
      end do
      !$OMP END DO
      !$OMP END PARALLEL
      !call msg('sum of modelvals before:', sum(work(1:size(physvals))))
      if (defined(correlation)) then
        call transf_adj(correlation, work(1:size(physvals)), controlvals)
      else
        controlvals = work
      end if
      !call msg('sum of modelvals after:', sum(controlvals))
    end subroutine apply_cov

  end subroutine to_control

  !************************************************************************************

  subroutine recover_control(control_phys, control, background)
    ! This subroutine attempts to perform inverse of to_physical. Works only if diagonal B
    ! is used.
    implicit none
    type(da_control), intent(in) :: control_phys, background
    type(da_control), intent(inout) :: control
    type(t_background_covariance), pointer :: p_bgr_cov

    character(len=*), parameter :: sub_name = 'recover_control'
    real, dimension(:), pointer :: p_values_phys, p_values_ctrl, p_values_bgr, p_stdev
    type(da_control) :: control_tmp
    type(t_background_covariance) :: bgr_cov_inv

    integer :: ii

    p_bgr_cov => control%bgr_cov_ptr
    if (fu_fails(associated(p_bgr_cov), 'bgr_cov not associated', 'to_control')) return

    if (fu_fails(defined(control_phys), 'control_phys not defined', sub_name)) return
    if (fu_fails(defined(control), 'control not defined', sub_name)) return
    if (fu_fails(fu_in_physical_space(control_phys), 'Control not in physical space', sub_name)) return
    if (fu_fails(.not. fu_in_physical_space(control), 'Control not in control space', sub_name)) return

    call copy_control(control_phys, control_tmp, if_allocate=.true.)

    ! subtract background
    p_values_phys => fu_values_ptr(control_tmp)
    p_values_phys = p_values_phys - fu_values_ptr(background)
    
    ! call to_control. The indexing will be changed. Now control is OK except that it has
    ! been multiplied with stdev instead of dividing.
    call to_control(control_tmp, control)
    if (error) return

    if (fu_mode(control) == DA_INITIAL_STATE .or. fu_mode(control) == DA_EMISSION_AND_INITIAL) then
      if (.not. fu_is_diagonal(p_bgr_cov%correlation_initial)) then
        call set_error('Initial state correlation not diagonal', sub_name)
        return
      end if
      p_values_ctrl => fu_initial_ptr(control)
      p_stdev => p_bgr_cov%stdev_initial
      do ii = 1, fu_n_initial(control)
        p_values_ctrl(ii) = p_values_ctrl(ii) / (p_stdev(ii)*p_stdev(ii))
      end do
    end if
    if (fu_mode(control) == DA_EMISSION_CORRECTION .or. fu_mode(control) == DA_EMISSION_AND_INITIAL &
      & .or. fu_mode(control) == DA_EMISSION_TIME_HEIGHT) then
      if (.not. fu_is_diagonal(p_bgr_cov%correlation_emission)) then
        call set_error('Initial state correlation not diagonal', sub_name)
        return
      end if
      p_values_ctrl => fu_emission_ptr(control)
      p_stdev => p_bgr_cov%stdev_emission
      do ii = 1, fu_n_emission(control)
        p_values_ctrl(ii) = p_values_ctrl(ii) / (p_stdev(ii)*p_stdev(ii))
      end do
    end if

    call destroy(control_tmp)
        
  end subroutine recover_control

  
  !************************************************************************************
  
  subroutine report_norm(control)
    implicit none
    type(da_control), intent(in) :: control
    
    logical :: no_correlation

    call msg('4DVAR: reporting control norm by species')

    if (control%shift_initial /= int_missing) then
      no_correlation = fu_in_physical_space(control)
      if (.not. no_correlation) then
        ! check if the correlation operator is void
        if (fu_fails(associated(control%bgr_cov_ptr), 'bgr_cov_ptr not associated', 'report_norm')) return
        no_correlation = .not. defined(control%bgr_cov_ptr%correlation_initial)
      end if
      if (no_correlation) then
        call msg('4DVAR: Initial state phys, nocorr:')
        call report_component_phys(fu_initial_ptr(control), fu_species_initial(control))
      else
        call msg('4DVAR: Initial state contr, (with corr):')
        call report_component_contr(fu_initial_ptr(control), fu_species_initial(control), &
                                  & control%bgr_cov_ptr%correlation_initial)
      end if
    end if

    if (control%shift_emission /= int_missing) then
      no_correlation = fu_in_physical_space(control)
      if (.not. no_correlation) then
        ! check if the correlation operator is void
        if (fu_fails(associated(control%bgr_cov_ptr), 'bgr_cov_ptr not associated', 'report_norm')) return
        no_correlation = .not. defined(control%bgr_cov_ptr%correlation_emission)
      end if
      if (no_correlation) then
        call msg('4DVAR: Emission  phys, nocorr:')
        call report_component_phys(fu_emission_ptr(control), fu_species_emission_ctrl(control) )
      else 
        call msg('4DVAR: Emission contr, (with corr):')
        call report_component_contr(fu_emission_ptr(control), fu_species_emission_ctrl(control), &
                                  & control%bgr_cov_ptr%correlation_emission)
      end if
    end if

  contains

    
    subroutine report_component_phys(values, species_list)
      implicit none
      real, dimension(:), intent(in) :: values
      type(silam_species), dimension(:), intent(in) :: species_list

      integer :: ind_species, ind_start, ind_end, field_size, num_species

      field_size = size(values) / size(species_list)
      ind_start = 1
      num_species = size(species_list)
      do ind_species = 1, num_species
        ind_end = ind_start + field_size - 1
        call msg('4DVAR: Species '//trim(fu_str(ind_species))//'--> ', sum(values(ind_start:ind_end)**2))
        ind_start = ind_end + 1
      end do      

    end subroutine report_component_phys

    !================================================================
    
    subroutine report_component_contr(values, species_list, correlation)
      implicit none
      real, dimension(:), intent(in) :: values
      type(silam_species), dimension(:), intent(in) :: species_list
      type(t_correlation), intent(in) :: correlation

      integer :: ind_species, ind_start, ind_end, field_size, num_species

      ind_start = 1
      num_species = size(species_list)
      do ind_species = 1, num_species
        field_size = fu_dim_control_3d(correlation, ind_species)
        ind_end = ind_start + field_size - 1
        call msg('4DVAR: Species '//trim(fu_str(ind_species))//'--> ', sum(values(ind_start:ind_end)**2))
        !call msg('4DVAR: ind_species, ind_start, field_size:', (/ind_species,ind_start, field_size/))
        !call msg('4DVAR: --> ', sum(values(ind_start:ind_end)**2))
        ind_start = ind_end + 1
      end do
      

    end subroutine report_component_contr
  end subroutine report_norm

  
  !************************************************************************************
    
  subroutine set_cov_mdl(bgr_cov, background, species_emission, species_transport, &
                       & rules, analysis_grid, analysis_vertical, analysis_time, control_var)
    ! 
    ! Set the background covariance structures following da_rules. This includes both
    ! correlation operators (per species) and standard deviations. The correlation
    ! operators are required to conform to analysis_grid and
    ! analysis_vertical. Analysis_time is used to find the covariance setup file from
    ! template. Species_emission and species_transport become the analysis_species unless
    ! a subset is required (only 3d-var).
    ! 
    implicit none
    type(t_background_covariance), intent(out) :: bgr_cov
    type(da_control), intent(in) :: background
    type(silam_species), dimension(:), intent(in), target :: species_transport, species_emission
    
    type(DA_Rules), intent(in) :: rules
    type(silja_grid), intent(in) :: analysis_grid
    type(silam_vertical), intent(in) :: analysis_vertical
    type(silja_time), intent(in) :: analysis_time
    ! rules%controlVariable not used - may use either controlVariable or perturbVariable
    integer, intent(in) :: control_var

    type(silam_species), dimension(:), pointer :: analysis_species
    type(Tsilam_namelist_group), pointer :: nlgrp_cov_setup
    type(Tsilam_namelist), pointer :: nlptr
    integer, dimension(max_species) :: indices
    integer :: file_unit, iostat, num_substances, num_species, stat, nx, ny, nz, num_analysis_species, ii
    integer :: stdev_size
    logical :: need_emission_xy, need_initial, need_emission_zt, need_volcano
    character(len=fnlen) :: filename
    type(silam_vertical) :: surface_vertical
    type(silja_time), dimension(:,:), pointer :: ptr_slots

    !
    ! Three controllers are supported (by at least one of the techniques): 
    ! emissino in plain, emission vertical-time, emission volcanoes, initial conditions
    !
    need_emission_xy = fu_have_emission_xy(control_var)
    need_emission_zt = fu_have_emission_zt(control_var)
    need_initial = fu_have_initial(control_var)
    need_volcano = fu_have_emission_volc(control_var)

    if (need_emission_zt .and. (need_emission_xy .or. need_initial)) then
      call set_error('Unsupported control variable combination', 'set_cov_mdl')
      return
    end if
    if (.not. (need_emission_xy .or. need_initial .or. need_emission_zt)) then
      call msg('controlvar:', control_var)
      call msg('Supported only: emission_zt, emission_xy, initial')
      call set_error('Unsupported control variable configuration', 'set_cov_mdl')
      return
    end if

    if(fu_fails(.not. need_volcano,'Setting background error not supported for volcano emission', 'set_cov_mdl'))return

    call grid_dimensions(analysis_grid, nx, ny)
    nz = fu_NbrOfLevels(analysis_vertical)
    
    nullify(bgr_cov%stdev_emission, bgr_cov%stdev_initial)
    
    if (rules%no_background_term) then
      call set_no_background(need_emission_xy, need_initial, need_emission_zt, bgr_cov)
      return
    end if
    !
    ! Emission 2D xy
    !
    if (need_emission_xy) then
      call expand_template(rules%cov_setup_templ_emis, analysis_time, filename)
      if (error) return
      file_unit = fu_next_free_unit()
      open(file_unit, file=filename, action='read', iostat=iostat)
      if (fu_fails(iostat == 0, 'Failed to open: ' // trim(filename), 'set_cov_mdl')) return
      nlgrp_cov_setup => fu_read_namelist_group(file_unit, .false.)
      if (error) return
      close(file_unit)

      ! note emission -> 2d, vertical = surface level
      num_analysis_species = size(species_emission)
      stdev_size = num_analysis_species * nx*ny
      allocate(bgr_cov%stdev_emission(stdev_size), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', 'set_cov_mdl')) return

      call set_vertical(surface_level, surface_vertical)
      call msg('Setting background covariances for emission')
      call set_cov(emission_scaling_flag, nlgrp_cov_setup, species_emission, &
                 & analysis_grid, surface_vertical, bgr_cov%correlation_emission, &
                 & fu_emission_ptr(background), bgr_cov%stdev_emission)
      call set_missing(surface_vertical, .false.)
      call destroy_namelist_group(nlgrp_cov_setup)
    end if
    !
    ! Concentratiobns (initial conditions)
    !
    if (need_initial) then 
      call expand_template(rules%cov_setup_templ_init, analysis_time, filename)
      if (error) return
      file_unit = fu_next_free_unit()
      open(file_unit, file=filename, action='read', iostat=iostat)
      if (fu_fails(iostat == 0, 'Failed to open: ' // trim(filename), 'set_cov_mdl')) return
      nlgrp_cov_setup => fu_read_namelist_group(file_unit, .false.)
      if (error) return
      close(file_unit)

      if (fu_true(rules%have_analysis_species)) then
        call pick_species(rules%analysis_subst_list_3d, species_transport, indices, num_analysis_species)
        if (fu_fails(num_analysis_species > 0, 'No analysis species', 'set_cov_mdl')) return
        allocate(analysis_species(num_analysis_species), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed', 'set_cov_mdl')) return
        analysis_species = (/(species_transport(indices(ii)), ii=1, num_analysis_species)/)
      else
        analysis_species => species_transport
      end if
      num_analysis_species = size(analysis_species)
      stdev_size = num_analysis_species * nx*ny*nz
      allocate(bgr_cov%stdev_initial(stdev_size), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', 'set_cov_mdl')) return
      call msg('Setting background covariances for concentration')
      call set_cov(concentration_flag, nlgrp_cov_setup, analysis_species, &
                 & analysis_grid, analysis_vertical, &
                 & bgr_cov%correlation_initial, fu_initial_ptr(background), bgr_cov%stdev_initial)
      if (error) return
      if (.not. associated(analysis_species, species_transport)) deallocate(analysis_species)
      call destroy_namelist_group(nlgrp_cov_setup)
    end if
    !
    ! Time-vertical emission profile
    !
    if (need_emission_zt) then
      stdev_size &
           & = fu_NbrOfLevels(analysis_vertical) * size(species_emission) * fu_num_emis_time_slots(rules)
      allocate(bgr_cov%stdev_emission(stdev_size), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', 'set_cov_mdl')) return
      call get_emis_time_slots(rules, ptr_slots=ptr_slots)
      if (error) return
      call set_cov_time_height(rules%cov_setup_templ_emis, analysis_vertical, species_emission, &
                             & ptr_slots, bgr_cov%stdev_emission, bgr_cov%correlation_emission)
      if (error) return
      deallocate(ptr_slots)
    end if

    bgr_cov%defined = .true.

  contains

    !========================================================

    subroutine set_no_background(need_emission, need_initial, need_time_height_emis, bgr_cov)
      implicit none
      logical, intent(in) :: need_emission, need_initial, need_time_height_emis
      type(t_background_covariance), intent(inout) :: bgr_cov
      
      integer :: stdev_size
      
      if (need_time_height_emis) then
        stdev_size &
             & = fu_NbrOfLevels(analysis_vertical) * size(species_emission) * fu_num_emis_time_slots(rules)
        allocate(bgr_cov%stdev_emission(stdev_size), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed', 'set_cov_mdl')) return
        bgr_cov%stdev_emission = 1.0
      else if (need_emission) then
        num_analysis_species = size(species_emission)
        stdev_size = num_analysis_species * nx*ny
        allocate(bgr_cov%stdev_emission(stdev_size), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed', 'set_cov_mdl')) return
        bgr_cov%stdev_emission = 1.0
      end if
      
      if (need_initial) then
        num_analysis_species = size(analysis_species)
        stdev_size = num_analysis_species * nx*ny*nz
        allocate(bgr_cov%stdev_initial(stdev_size), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed', 'set_cov_mdl')) return
        bgr_cov%stdev_initial = 1.0
      end if
      
      call set_missing(bgr_cov%correlation_emission)
      call set_missing(bgr_cov%correlation_initial)
      bgr_cov%defined = .true.

    end subroutine set_no_background

    !===================================================================
    
    subroutine expand_template(templ, now, filename)
      implicit none
      type(grads_template), intent(in) :: templ
      type(silja_time), intent(in) :: now
      character(len=*), intent(out) :: filename
      
      type(silam_sp), dimension(:), pointer :: filenames
      
      nullify(filenames)
      
      call fnm_from_single_template(templ, now, filenames, &
                                  & ifStrict = .true., &
                                  & ifadd = .false., &
                                  & ifWait = .false., &
                                  & ifAllowZeroFcLen = .true.)
      if (size(filenames) < 1 .or. error) then
        call set_error('No configuration found for time:' // fu_str(now), &
                     & 'set_cov_mdl')
      end if
      if (size(filenames) > 1) then
        call set_error('Multiple configurations found for time:' // fu_str(now), &
                     & 'set_cov_mdl')
      end if
      if (error) return

      filename = fu_process_filepath(filenames(1)%sp, must_exist=.true.)
      deallocate(filenames(1)%sp)
      deallocate(filenames)
    end subroutine expand_template

    !=================================================================

    subroutine set_cov(quantity, nlgrp, analysis_species, analysis_grid, analysis_vertical, &
                     & correlation, bgr_values, stdev)
      implicit none
      integer, intent(in) :: quantity
      type(Tsilam_namelist_group) :: nlgrp
      type(silja_grid), intent(in) :: analysis_grid
      type(silam_vertical), intent(in) :: analysis_vertical
      type(silam_species), dimension(:), intent(in) :: analysis_species
      type(t_correlation), intent(inout) :: correlation
      real, dimension(:), target, intent(in) :: bgr_values
      real, dimension(:), target, intent(out) :: stdev
      
      logical, dimension(max_species) :: species_ok
      type(Tsilam_nl_item_ptr), dimension(:), pointer :: p_items
      integer :: num_substances, num_species, nx, ny, ind_lev, num_vals_3d, ind_nl, &
           & ind_selected, ind_species
      integer :: ind_subst, quantity_flag, num_selected, ind_start, ind_start_lev, stride, stride_lev, stat
      integer, dimension(max_species) :: subst_indices
      real, dimension(:), pointer :: stdev_species, bgr_species
      character(len=fnlen) :: file_name, quantity_name, stdev_input
      character(len=substNmLen) :: subst_name
      real :: stdev_const, stdev_a, stdev_b
      type(t_spatial_correlation), dimension(:), pointer :: correlations

      p_items => null()
      num_species = size(analysis_species)

      allocate(correlations(num_species), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', 'set_cov_mdl')) return
      species_ok(:) = .false.

      call grid_dimensions(analysis_grid, nx, ny)
      num_vals_3d = nx*ny*fu_NbrOfLevels(analysis_vertical)
      
      ! Find the namelists valid for the quantity
      do ind_nl = 1, fu_nbr_of_namelists(nlgrp_cov_setup)
        nlptr => fu_namelist(nlgrp_cov_setup, ind_nl)
        quantity_name = fu_content(nlptr, 'quantity')
        if (fu_fails(quantity_name /= '', 'Missing quantity', 'set_cov_mdl')) return
        quantity_flag = fu_get_silam_quantity(quantity_name)
        if (error) return
        if (quantity_flag /= quantity) cycle
        call get_items(nlptr, 'substance', p_items, num_substances)
        if (num_substances < 1) then
          call msg('ind_nl', ind_nl)
          call set_error('Cov namelist has no substance items',  'set_cov_mdl')
          return
        end if

        ! Fill in the correlation & stdev for the species found.
        do ind_subst = 1, num_substances
          subst_name = fu_content(p_items(ind_subst))
          if (fu_fails(subst_name /= '', 'Empty substance', 'set_cov_mdl')) return
          if (subst_name == '*') subst_name = char_missing ! accept all
          call select_species(analysis_species, num_species, subst_name, aerosol_mode_missing, &
                            & real_missing, subst_indices, num_selected)

          do ind_selected = 1, num_selected
            ind_species = subst_indices(ind_selected)
            ! do not override with '*' as substance
            if (species_ok(ind_species) .and. subst_name == char_missing) cycle 
            call set_spatial_corr_from_nl(nlptr, analysis_grid, correlations(ind_species))
            if (error) return

            !ind_start_species = (ind_species-1) * num_vals_3d + 1
            !stdev_species => stdev(ind_start_species:)

            stdev_species => stdev(ind_species::num_species)

            select case(fu_content(nlptr, 'stdev_method'))
            case ('constant')
              stdev_const = fu_content_real(nlptr, 'stdev_const')
              if (stdev_const .eps. real_missing) then
                call set_error('Missing stdev_const', 'set_cov_mdl')
                return
              else if (stdev_const <= 0) then
                call set_error('Strange stdev_const', 'set_cov_mdl')
                return
              end if
              stdev_species(1:num_vals_3d) = stdev_const
              if (debug_level > -1) then
                call report(analysis_species(ind_species))
                call msg('-- constant stdev:', stdev_const)
              end if

            case ('from_file')
              file_name = fu_content(nlptr, 'stdev_file')
              if (file_name == '') then
                call set_error('Missing stdev_file', 'set_cov_mdl')
                return
              end if
              do ind_lev = 1, fu_NbrOfLevels(analysis_vertical)
                ind_start_lev = (ind_lev-1)*nx*ny + 1
                stride_lev = fu_NbrOfLevels(analysis_vertical)
                call stdev_from_grads(file_name, quantity, analysis_species(ind_species), &
                                    & fu_level(analysis_vertical, ind_lev), analysis_grid, &
                                    & stdev_species(ind_start_lev::stride_lev))
                if (error) return
              end do
              if (debug_level > -1) then
                call report(analysis_species(ind_species))
                call msg('-- variable stdev, mean:', sum(stdev_species(1:num_vals_3d)/(num_vals_3d)))
              end if

            case ('linear')
              ! linear stdev = a*x_backgr + b
              stdev_input = fu_content(nlptr, 'stdev_linear')
              if (fu_fails(stdev_input /= '', 'Missing stdev_linear', 'set_cov_mdl')) return
              read(unit=stdev_input, fmt=*, iostat=stat) stdev_a, stdev_b
              if (fu_fails(stat == 0, 'Failed to parse: '//trim(stdev_input), 'set_cov_mdl')) return
              bgr_species => bgr_values(ind_species::num_species)
              stdev_species = stdev_a*bgr_species + stdev_b
              if (debug_level > -1) then
                call report(analysis_species(ind_species))
                call msg('-- linear stdev, mean:', sum(stdev_species(1:num_vals_3d)/(num_vals_3d)))
                call msg('-- mean background:', sum(bgr_species(1:num_vals_3d)/num_vals_3d))
              end if

            case ('')
              call set_error('Missing stdev_method', 'set_cov_mdl')
              return

            case default
              call set_error('Strange stdev_method: ' // trim(fu_content(nlptr, 'stdev_method')), &
                           & 'set_cov_mdl')
              return

            end select
            
            species_ok(ind_species) = .true.
          end do ! selected species
          
        end do ! substance
        deallocate(p_items)
        
      end do ! namelist
      do ind_species = 1, num_species
        if (.not. species_ok(ind_species)) then
          call report(analysis_species(ind_species))
          call set_error('No covariance setup found for species', 'set_cov_mdl')
        end if
      end do

      ! Finally combine the separate correlations into one:
      call set_total_correlation(correlations, analysis_grid, analysis_vertical, num_species, correlation)
      !!FIXME Should be deallocated, but used elsewhere via pointers
      !!deallocate (correlations) 
    end subroutine set_cov

    !========================================================
    
    subroutine stdev_from_grads(file_name, quantity, species, level, analysis_grid, stdev)
      implicit none
      character(len=*), intent(in) :: file_name
      integer, intent(in) :: quantity
      type(silam_species), intent(in) :: species
      type(silja_level), intent(in) :: level
      type(silja_grid), intent(in) :: analysis_grid
      real, dimension(:), intent(out) :: stdev

      integer :: ind_gf, ind_species, ind_start, fs, ii, ind_lev, stat, num_species, nx, ny, nz
      real, dimension(:), pointer :: stdev_species, data_in
      type(silja_time) :: when
      type(silja_field_id) :: fid

      ind_gf = fu_open_gradsfile_i(file_name)
      if (error) return

      if (.not. (fu_silamGrid_of_grads(ind_gf) == analysis_grid)) then
        call msg('Grid for emission correction background field:')
        call report(fu_silamGrid_of_grads(ind_gf))
        call msg('Analysis grid:')
        call report(analysis_grid)
        call set_error('Bad emission correction background grid', 'stdev_from_grads')
        return
      end if
      when = fu_time_of_grads(ind_gf, 1)
      fid = fu_set_field_id(met_src_missing, quantity, when, zero_interval, &
                          & analysis_grid, level, &
                          & species=species)
      call read_field_from_grads_id(ind_gf, fid, stdev)
      
      call close_gradsfile_i(ind_gf)

    end subroutine stdev_from_grads

    !========================================================

    subroutine set_cov_time_height(cov_template, analysis_vertical, analysis_species, time_slots, &
                                 & stdev, correlation_emission)
      ! The time_height case works differently. The template is expanded with regard to
      ! each time slot, and standard deviations are then read with a scheme similar to the
      ! gridded case.
      implicit none
      type(grads_template), intent(in) :: cov_template
      type(silam_vertical), intent(in) :: analysis_vertical
      type(silam_species), dimension(:), intent(in) :: analysis_species
      type(silja_time), dimension(:,:), intent(in) :: time_slots
      real, dimension(:), intent(out), target :: stdev
      type(t_correlation), intent(out) :: correlation_emission ! corr times only

      integer :: ind_slot, file_unit, num_substances, quantity_flag, num_species, ind_lev, num_vals_1d
      integer :: ind_nl, ind_species, ind_start, ind_stdev_item, num_levs, num_stdev_items, stride
      integer :: iostat, ind_selected, ind_subst, num_selected
      character(len=fnlen) :: filename, content
      character(len=clen) :: quantity_name
      character(len=substNmLen) :: subst_name
      !type(Tsilam_namelist_group), pointer :: nlgrp
      logical, dimension(size(analysis_species)) :: species_ok
      real, dimension(:), pointer :: stdev_1d, corr_times
      type(Tsilam_nl_item_ptr), dimension(:), pointer :: p_items, p_stdev_items
      real :: stdev_default, stdev_lev, corr_time, content_real
      type(silja_time) :: time_start
      integer, dimension(:), pointer :: subst_indices
      type(t_spatial_correlation), dimension(:), pointer :: spatial_correlations
      character(len=*), parameter :: sub_name = 'set_cov_time_height'
      real :: thickness_total, trg_a

      nullify(corr_times)
      species_ok = .false.
      num_species = size(analysis_species)
      num_levs = fu_NbrOfLevels(analysis_vertical)
      call msg('Setting background covariance for time_height_emissions')
  
      allocate(spatial_correlations(num_species), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', 'set_cov_time_height')) return
  
      do ind_slot = 1, size(time_slots, 2)
        call msg('Time slot:', ind_slot)
        time_start = time_slots(slot_start, ind_slot)
        call expand_template(cov_template, time_start, filename)
        if (error) return
        file_unit = fu_next_free_unit()
        open(file_unit, file=filename, action='read', iostat=iostat)
        if (fu_fails(iostat == 0, 'Failed to open: ' // trim(filename), sub_name)) return
        nlgrp_cov_setup => fu_read_namelist_group(file_unit, .false.)
        close(file_unit)
        if (error) return

        species_ok = .false.
        
        do ind_nl = 1, fu_nbr_of_namelists(nlgrp_cov_setup)
          nlptr => fu_namelist(nlgrp_cov_setup, ind_nl)
          quantity_name = fu_content(nlptr, 'quantity')
          if (fu_fails(quantity_name /= '', 'Missing quantity', sub_name)) return
          quantity_flag = fu_get_silam_quantity(quantity_name)
          if (error) return
          if (.not. (quantity_flag == emission_scaling_flag &
                   & .or. quantity_flag == emission_intensity_flag)) cycle
          call get_items(nlptr, 'substance', p_items, num_substances)
          if (num_substances < 1) then
            call msg('ind_nl', ind_nl)
            call set_error('Cov namelist has no substance items',  sub_name)
            return
          end if

          ! Fill in the correlation & stdev for the species found.
          subst_indices => fu_work_int_array()
          
          do ind_subst = 1, num_substances
            subst_name = fu_content(p_items(ind_subst))
            if (fu_fails(subst_name /= '', 'Empty substance', sub_name)) return
            if (subst_name == '*') subst_name = char_missing ! accept all
            call select_species(analysis_species, num_species, subst_name, aerosol_mode_missing, &
                              & real_missing, subst_indices, num_selected)

            call msg('Substance: ' // trim(subst_name) // ', n. analysis species:', num_selected)
            do ind_selected = 1, num_selected
              ind_species = subst_indices(ind_selected)
              ! do not override with '*' as substance
              if (species_ok(ind_species) .and. subst_name == char_missing) cycle 

              call set_spatial_corr_from_nl(nlptr, grid_missing, spatial_correlations(ind_species))
              if (error) return

              ! the stdev has indexing (ispecies,iz,ind_time)
              ind_start = (ind_slot-1) * num_species*num_levs + ind_species 
              stride = num_species
              stdev_1d => stdev(ind_start::stride)
              stdev_default = fu_content_real(nlptr, 'stdev_default')
              if (.not. (stdev_default .eps. real_missing)) then 
                stdev_1d(1:num_levs) = stdev_default
              else
                stdev_1d(1:num_levs) = -1
              end if
              call msg('Default stdev:', stdev_default, ind_start)
              call get_items(nlptr, 'stdev', p_stdev_items, num_stdev_items)
              
              if (fu_content(nlptr, 'stdev_method') == 'gauss_squared') then
                ! perturbation is on emission density (per vertical meter). Expectation is
                ! density of one for each layer. Scale stdev_1d so that the expected
                ! column emission is equatl to column_prior_mean.
                content_real = fu_content_real(nlptr, 'column_prior_mean')
                if (content_real .eps. real_missing) then
                  call set_error('trunc_gauss requires column_prior_mean', sub_name)
                  return
                end if
                stdev_1d = content_real &
                     & / sum((/(fu_layer_thickness_m(fu_level(analysis_vertical, ind_lev)), &
                              & ind_lev=1, num_levs)/)) 
              end if
              
              do ind_stdev_item = 1, num_stdev_items
                content = fu_content(p_stdev_items(ind_stdev_item))
                read(unit=content, fmt=*, iostat=iostat) ind_lev, stdev_lev
                if (fu_fails(iostat == 0, 'Failed to parse: ' // trim(content), sub_name)) return
                stdev_1d(ind_lev) = stdev_lev
              end do
              if (any(stdev_1d(1:num_levs) < 0)) then
                call msg('Missing stdev for ' // fu_substance_name(analysis_species(ind_species)))
                call set_error('stdev not found', 'set_cov_time_height')
                return
              end if

              content = fu_content(nlptr, 'correlation_time')
              if (content /= '') then
                corr_time = fu_sec(fu_set_named_interval(content))
                if (error) return
                ! this only works in ensemble mode, where we have not slots but correlation time.
                if (size(time_slots, 2) /= 1) then
                  call set_error('Cannot handle correlation times if > 1 time slot', sub_name)
                  return
                end if
                if (.not. associated(corr_times)) then
                  ! corr times indexes as ind_species, ind_level 
                  allocate(corr_times(num_levs*num_species), stat=stat)
                  corr_times = real_missing
                  if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
                end if
                ind_start = ind_species 
                stride = num_species
                corr_times(ind_start::stride) = corr_time
              end if

              species_ok(ind_species) = .true.
            end do ! selected species
            
          end do ! substance
          
          call free_work_array(subst_indices)
          
        end do ! nl
        
        do ind_species = 1, num_species
          if (.not. species_ok(ind_species)) then
            call report(analysis_species(ind_species))
            call set_error('No covariance setup found for species', sub_name)
          end if
        end do
        
      end do ! slot

      call set_total_correlation(spatial_correlations, grid_missing, analysis_vertical, num_species, &
                               & correlation_emission)
      if (error) return

      if (associated(corr_times)) then
        call msg('Correlation times defined')
        if (any(corr_times .eps. real_missing)) then
          call set_error('Not all correlation times set', 'set_cov_time_height')
          return
        end if
        call set_corr_times(correlation_emission, corr_times)
      end if

    end subroutine set_cov_time_height

    !============================================================                             
                                 
    real function fu_trg_stdev(trg_a, corr, corr_size) result(stdev)
      ! Computes the standard deviation for the column emission with the
      ! truncated-gaussian model. This is
      ! Var(sum(Y_i)) = a**2 * sx**2 sum(r_ij h_i h_j), where
      ! sx**2 = 1 - 2/pi
      ! r_ij is correlation and h_i, h_j are layer thicknesses.
      implicit none
      real, intent(in) :: trg_a
      type(t_spatial_correlation) :: corr
      integer, intent(in) :: corr_size

      real, dimension(:), pointer :: h_i
      real, dimension(:,:), pointer :: r_ij
      
      integer :: ii, jj
      real, parameter :: sx2 = 1 - 2/pi

      h_i => fu_work_array()
      r_ij => fu_work_array_2d()
      
      call get_correlation_matrix(corr, r_ij)
      if (error) return
      
      h_i(1:fu_NbrOfLevels(corr%vertical)) &
           & = (/(fu_layer_thickness_m(fu_level(corr%vertical, ii)), ii=1, fu_NbrOfLevels(corr%vertical))/)
      stdev = 0.0
      do jj = 1, corr_size
        do ii = 1, corr_size
          stdev = stdev + r_ij(ii,jj) * h_i(ii) * h_i(jj)
        end do
      end do
      stdev = sqrt(stdev * trg_a**2 * sx2)
      
      call free_work_array(h_i)
      call free_work_array(r_ij)

    end function fu_trg_stdev

    !===============================================================

    subroutine next_selected_species(subst_items, num_subst_items, species_list, ind_species, species_done)
      implicit none
      type(Tsilam_nl_item_ptr), dimension(:), intent(in) :: subst_items
      integer, intent(in) :: num_subst_items
      type(silam_species), dimension(:), intent(in) :: species_list
      integer, intent(inout) :: ind_species
      logical, dimension(:), intent(inout) :: species_done

      logical, save :: initialized = .false., list_done
      logical :: begin, skip_species
      integer, dimension(:), pointer, save :: subst_indices
      character(len=substNmLen), save :: subst_name
      integer :: ind_subst, ind_selected, num_selected

      if (fu_fails(ind_species > -1, 'Bad initial ind_species', 'iterate_selected_species')) return
      
      begin = ind_species == 0
      if (begin) then ! begin iteration
        if (initialized .and. associated(subst_indices)) call free_work_array(subst_indices)
        subst_indices => fu_work_int_array()
        ind_subst = 0
        ind_selected = 0
        initialized = .true.
        list_done = .false.
      end if

      do ! only once unless skip species
        do while ((begin .or. num_selected < 1 .or. list_done))
          ind_subst = ind_subst + 1

          if (ind_subst > num_subst_items) then
            ! Terminate:
            ind_species = -1
            call free_work_array(subst_indices)
            nullify(subst_indices)
            return
          end if

          subst_name = fu_content(subst_items(ind_subst))
          if (fu_fails(subst_name /= '', 'Empty substance', 'set_cov_mdl')) return
          if (subst_name == '*') subst_name = char_missing ! accept all
          call select_species(species_list, size(species_list), subst_name, aerosol_mode_missing, &
                            & real_missing, subst_indices, num_selected)
          if (error) return
          ind_selected = 0
        end do

        ind_selected = ind_selected + 1
        
        ! pick index for species counter and return
        ind_species = subst_indices(ind_selected)
        skip_species = (species_done(ind_species) .and. subst_name == char_missing)
        species_done(ind_species) = .true.
        
        list_done = ind_selected == num_selected
        
        if (.not. skip_species) exit ! skip loop

      end do ! skip loop

    end subroutine next_selected_species

  end subroutine set_cov_mdl
  

  !*************************************************************************
  
  subroutine kill_bgr_cov(bgr_cov)
    implicit none
    type(t_background_covariance), intent(inout) :: bgr_cov

    if (.not. bgr_cov%defined) return

    if (associated(bgr_cov%stdev_initial)) deallocate(bgr_cov%stdev_initial)
    if (associated(bgr_cov%stdev_emission)) deallocate(bgr_cov%stdev_emission)
    call kill_total_correlation(bgr_cov%correlation_initial)
    call kill_total_correlation(bgr_cov%correlation_emission)
  end subroutine kill_bgr_cov


  !*************************************************************************
  
  logical function fu_defined_bgr_cov(bgr_cov) result(is_defined)
    implicit none
    type(t_background_covariance), intent(in) :: bgr_cov
    
    is_defined = bgr_cov%defined

  end function fu_defined_bgr_cov

  !************************************************************************************  

  integer function fu_control_space_dim(bgr_cov) result(nn)
    implicit none
    type(t_background_covariance), intent(in) :: bgr_cov
    
    integer :: dim_init, dim_emis
    
    nn = 0

    call control_space_dim(bgr_cov, dim_init, dim_emis)
    nn = dim_init + dim_emis

  end function fu_control_space_dim


  !*************************************************************************
  
  subroutine control_space_dim(bgr_cov, dim_init, dim_emis)
    implicit none
    type(t_background_covariance), intent(in) :: bgr_cov
    integer, intent(out) :: dim_init, dim_emis

    dim_init = 0
    dim_emis = 0

    if (defined(bgr_cov%correlation_initial)) then 
      dim_init = fu_dimension_control(bgr_cov%correlation_initial)
    else if (associated(bgr_cov%stdev_initial)) then
      dim_init = size(bgr_cov%stdev_initial)
    end if

    if (defined(bgr_cov%correlation_emission)) then 
      dim_emis = fu_dimension_control(bgr_cov%correlation_emission)
    else if (associated(bgr_cov%stdev_emission)) then
      dim_emis = size(bgr_cov%stdev_emission)
    end if
    
  end subroutine control_space_dim

  !************************************************************************************

  subroutine get_posit_constr(control_template, background, bgr_cov, lower_bounds, bounds_defined)
    ! Get the positivity constraint (lower boundary) in control space. Depends on the
    ! background standard deviation, and only possible if the correlation matrix is
    ! diagonal.
    ! 
    ! For diagonal B, 0 < x_physical = stdev * x_control + x_background, so the positivity
    ! constraint is x_control > -x_background / stdev.
    ! 
    ! The diagonality is required blockwise; this is now by control variable component
    ! (emission/initial), however, in absence of chemical correlations it would be
    ! possible to extend this to species level.
    implicit none
    type(da_control), intent(in) :: background, control_template
    type(t_background_covariance), intent(in) :: bgr_cov
    
    ! this is used by l_bfgs_b, so real(8)
    real(r8k), dimension(:), intent(out) :: lower_bounds
    integer, dimension(:), intent(out) :: bounds_defined
    
    real, dimension(:), pointer :: stdev, mask, set_mask, values
    
    type(da_control) :: control_bnd, control_stdev

    call msg('Checking positivity constraint...')

    ! The trick is to use the operations on control variables: knowing their mathematical
    ! definition we don't care of the physical meaning.
    
    if (fu_fails(.not. fu_in_physical_space(control_template), 'Need template in control space', &
               & 'get_posit_constr')) return
    call copy_control(control_template, control_bnd, if_allocate=.true.) ! in control space
    call copy_control(background, control_stdev, if_allocate=.true.) ! in physical space
    values => fu_values_ptr(control_bnd)
    values = 1.0
    call to_physical(control_bnd, background, control_stdev)
    values => fu_values_ptr(control_stdev)
    ! now values = stdev + background
    values = values - fu_values_ptr(background)
    stdev => values
    where (stdev > 0)
      lower_bounds = -fu_values_ptr(background) / stdev
    elsewhere
      ! may happen if B is non-diagonal. Lower boundary will not be set (see below).
      lower_bounds = real_missing
    end where
    ! check that constraint can be used for each component: if correlation is defined, use
    ! the control_bnd to mask out elements corresponding to that component.
    bounds_defined = 1
    mask => fu_values_ptr(control_bnd)
    mask = 0

    if (.not. fu_is_diagonal(bgr_cov%correlation_initial)) then
      call msg('...constraint not available for initial state')
      set_mask => fu_initial_ptr(control_bnd)
      set_mask = 1
    else
      call msg('...constraint set for intitial state')
    end if

    if (.not. fu_is_diagonal(bgr_cov%correlation_emission)) then
      call msg('...constraint not available for emission')
      set_mask => fu_emission_ptr(control_bnd)
      set_mask = 1
    else
      call msg('...constraint set for emission')
    end if
        
    where (mask > 0) bounds_defined = 0

    call destroy(control_bnd)
    call destroy(control_stdev)
    
  end subroutine get_posit_constr


  !*************************************************************************
  

  ! Inquiry functions to check whether the control variable includes a particular
  ! component.
  logical function fu_have_emission_xy(controlvar)
    integer, intent(in) :: controlvar
    fu_have_emission_xy = iand(controlvar, control_emis_xy) /= 0
  end function fu_have_emission_xy
  
  logical function fu_have_emission_zt(controlvar)
    integer, intent(in) :: controlvar
    fu_have_emission_zt = iand(controlvar, control_emis_zt) /= 0
  end function fu_have_emission_zt

  logical function fu_have_initial(controlvar)
    integer, intent(in) :: controlvar
    fu_have_initial = iand(controlvar, control_init) /= 0
  end function fu_have_initial

  logical function fu_have_emission_volc(controlvar)
    integer, intent(in) :: controlvar
    fu_have_emission_volc = iand(controlvar, control_emis_volc) /= 0
  end function fu_have_emission_volc

  logical function fu_have_meteo_time(controlvar)
    integer, intent(in) :: controlvar
    fu_have_meteo_time = iand(controlvar, control_meteo_time) /= 0
  end function fu_have_meteo_time
  
  function fu_get_control_bits(name) result(bits)
    character(len=*), intent(in) :: name
    integer :: bits
    select case(name)
    case ('initial_state')
      bits = control_init
    case ('emission_correction_xy')
      bits = control_emis_xy
    case ('emission_zt')
      bits = control_emis_zt
    case ('emission_volc')
      bits = control_emis_volc
    case default
      call set_error('Bad control variable: ' // trim(name), 'fu_get_control_bits')
    end select
  end function fu_get_control_bits


  !*************************************************************************
  
!  logical function fu_control_includes(controlvar, name) result(includes)
!    implicit none
!    integer, intent(in) :: controlvar
!    character(len=*), intent(in) :: name
!    includes = iand(fu_get_control_bits(name), controlvar) /= 0
!  end function fu_control_includes

  !************************************************************************************

  subroutine control_to_file(control, background, darules, filename, valid_time)
    ! 
    ! Subroutine for writing control variables to files. Supports only grads. A newer
    ! version could be made based on the subroutines below this subroutine (currently only
    ! used for the EnK smoother output).
    implicit none
    type(da_control), intent(in) :: control, background
    type(da_rules), intent(in) :: darules
    character(len=*), intent(in) :: filename
    type(silja_time), intent(in), optional :: valid_time

    integer :: ind_file, iostat
    type(da_control) :: control_output
    type(silja_grid) :: grid


    logical :: destr_control_output
    character(len=*), parameter :: sub_name = 'control_to_file'
    type(silja_time) :: valid_time_

    grid=control%grid

!    if (fu_control_includes(fu_mode(control), 'emission_volc')) then
    if(fu_have_emission_volc(fu_mode(control)))then
      call set_error('Control_to_file not supported for volc emis assimilation', sub_name)
      return
    end if

    if (present(valid_time)) then
      valid_time_ = valid_time
    else
      valid_time_ = darules%da_begin
    end if

    if (fu_in_physical_space(control)) then
      control_output = control
      destr_control_output = .false.
    else
      if (fu_fails(defined(background), 'Undefined background but control not physical', sub_name)) return
      call copy_control(background, control_output, if_allocate=.true.)
      if (error) return
      destr_control_output = .true.
      call to_physical(control, background, control_output)
      if (error) return
    end if

    select case (fu_mode(control))
    case (DA_INITIAL_STATE)
      ind_file = open_gradsfile_o('', filename, grid)
      if (error) return
      call msg('::initial to grads')
      call initial_to_grads(fu_initial_ptr(control_output), grid, &
                          & dispersion_vertical, fu_species_initial(control), &
                          & valid_time_, ind_file)
      call close_gradsfile_o(ind_file,"")
    
    case (DA_EMISSION_CORRECTION)
      ind_file = open_gradsfile_o('', filename, grid)
      if (error) return
      call emission_to_grads(fu_emission_ptr(control_output), grid, &
                           & fu_species_emission_ctrl(control), &
                           & valid_time_, ind_file)
      call close_gradsfile_o(ind_file,"")

    case (DA_EMISSION_AND_INITIAL)
      ind_file = open_gradsfile_o('', filename, grid)
      if (error) return
      call initial_to_grads(fu_initial_ptr(control_output), grid, &
                          & dispersion_vertical, fu_species_initial(control), &
                          & valid_time_, ind_file)
      call emission_to_grads(fu_emission_ptr(control_output), grid, &
                           & fu_species_emission_ctrl(control), &
                           & valid_time_, ind_file)
      call close_gradsfile_o(ind_file,"")
      
    case (DA_EMISSION_TIME_HEIGHT)
      ind_file = fu_next_free_unit()
      open(file=filename, unit=ind_file, form='formatted', action='write', iostat=iostat)
      if (fu_fails(iostat == 0, 'Failed to open: ' // trim(filename), 'control_to_file')) return
      call emis_time_hgt_to_file(fu_emission_ptr(control_output), dispersion_vertical, &
                               & fu_species_emission_ctrl(control), darules, ind_file)
      close(ind_file)

    case default
      call set_error('Unsupported control_mode', 'control_to_file')
      return

    end select
    
    if (destr_control_output) call destroy(control_output)

  contains
    !======================================================
    
    subroutine initial_to_grads(data_4d, grid, vertical, species_list, valid_time, ind_grads_file)
      implicit none
      real, dimension(:), intent(in), target :: data_4d
      type(silja_grid), intent(in) :: grid
      type(silam_vertical), intent(in) :: vertical
      type(silam_species), dimension(:), intent(in) :: species_list
      type(silja_time) :: valid_time
      integer, intent(in) :: ind_grads_file
      
      integer :: nspecies, fs, nx, ny, ind_start, ind_end, step, ilev, isp, nlevs
      real, dimension(:), pointer :: data_2d
      type(silja_field_id) :: fid

      nspecies = size(species_list)
      nlevs = fu_NbrOfLevels(vertical)
      call grid_dimensions(grid, nx, ny)
      fs = nx*ny
      
      do isp = 1, nspecies
        do ilev = 1, nlevs
          ind_start = isp + (ilev-1)*nspecies
          step = nlevs*nspecies
          ind_end = ind_start + (fs-1)*step

          data_2d => data_4d(ind_start:ind_end:step)
          fid = fu_set_field_id(met_src_missing, concentration_flag, &
                              & valid_time, zero_interval, &
                              & grid, fu_level(vertical, ilev), &
                              & species=species_list(isp))
          call write_next_field_to_gradsfile(ind_grads_file, fid, data_2d)
          if (error) return
        end do
      end do

    end subroutine initial_to_grads

    !===============================================================
    
    subroutine emission_to_grads(data_3d, grid, species_list, valid_time, ind_grads_file)
      implicit none
      real, dimension(:), intent(in), target :: data_3d
      type(silja_grid), intent(in) :: grid
      type(silam_species), dimension(:), intent(in) :: species_list
      type(silja_time), intent(in) :: valid_time
      integer, intent(in) :: ind_grads_file
      
      integer :: nspecies, fs, nx, ny, ind_start, ind_end, step, isp
      real, dimension(:), pointer :: data_2d
      type(silja_field_id) :: fid

      nspecies = size(species_list)
      call grid_dimensions(grid, nx, ny)
      fs = nx*ny
      !call msg('::nspecies', nspecies)
      do isp = 1, nspecies
        ind_start = isp
        step = nspecies
        ind_end = ind_start + (fs-1)*step
        data_2d => data_3d(ind_start:ind_end:step)
        fid = fu_set_field_id(met_src_missing, emission_scaling_flag, &
                            & valid_time, zero_interval, &
                            & grid, surface_level, &
                            & species=species_list(isp))
        call write_next_field_to_gradsfile(ind_grads_file, fid, data_2d)
        if (error) return
      end do

    end subroutine emission_to_grads

  end subroutine control_to_file
  
  
  !*********************************************************************************
  
  subroutine emis_time_hgt_to_file(values, vertical, species_list, darules, ind_file)
    implicit none
    real, dimension(:), intent(in) :: values
    type(silam_vertical), intent(in) :: vertical
    type(silam_species), dimension(:), intent(in) :: species_list
    type(da_rules), intent(in) :: darules
    integer, intent(in) :: ind_file

    type(silja_time), dimension(:,:), pointer :: time_slots
    integer :: ind_slot, ind, ind_lev, ind_species, num_levs, num_species
    character(len=clen) :: time_str, species_str
    real, dimension(:), pointer :: line

    call get_emis_time_slots(darules, time_slots)
    if (error) return

    line => fu_work_array()

    num_levs = fu_NbrOfLevels(vertical)
    num_species = size(species_list)
    do ind_slot = 1, size(time_slots, 2)
      time_str = fu_str(time_slots(slot_start,ind_slot))
      do ind_species = 1, num_species
        species_str = fu_species_output_name(species_list(ind_species))
        do ind_lev = 1, num_levs
          ind = (ind_slot-1)*num_species*num_levs + (ind_lev-1)*num_species + ind_species
          line(ind_lev) = values(ind)
          !call msg('contr_to_file:', ind, values(ind))
        end do
        write(ind_file, fmt='(A, " ", A, 100G10.2)') trim(time_str), trim(species_str), line(1:num_levs)
      end do
    end do

    deallocate(time_slots)
    call free_work_array(line)

  end subroutine emis_time_hgt_to_file

  !************************************************************************************
  !
  ! Below is a newer implementation of control_to_file, although not yet entirely in use.
  ! 
  ! The control may contain both field-like components (like initial state) and others
  ! (like z-t emission). The field-like components can be output using the standard SILAM
  ! IO using the control_to_fieldset subroutine. The presence of non-field components can
  ! be inquired using the get_non_fieldvars_for_control subroutine, and if any exist, they
  ! can be written to a text file(s) using control_non_fields_to_files. 

  subroutine get_fieldvars_for_control(control, num_vars, quantities, species_list, if3d)
    ! Inquire the type and species for field-type output variables from the control.
    implicit none
    type(da_control), intent(in) :: control
    integer, intent(out), dimension(:), optional :: quantities
    type(silam_species), dimension(:), intent(out), optional :: species_list
    logical, intent(out), dimension(:), optional :: if3d
    integer, intent(out) :: num_vars

    character(len=*), parameter :: subname = 'get_fieldvars_for_control'
    type(silam_species), dimension(:), pointer :: species_initial, species_emission
    integer :: ind_out_list, ind_species, max_out_list


    max_out_list = huge(kind(num_vars))
    if (present(quantities)) max_out_list = min(size(quantities), max_out_list)
    if (present(species_list)) max_out_list = min(size(species_list), max_out_list)
    if (present(if3d)) max_out_list = min(size(if3d), max_out_list)
    if (fu_fails(max_out_list > 0, 'Something wrong with quantites, species_list or if3d', subname)) return

    ind_out_list = 0
    if (fu_have_initial(fu_mode(control))) then
      species_initial => fu_species_initial(control)
      do ind_species = 1, size(species_initial)
        ind_out_list = ind_out_list + 1
        if (fu_fails(ind_out_list <= max_out_list, 'One of output arguments too small', subname)) return
        if (present(species_list)) species_list(ind_out_list) = species_initial(ind_species)
        if (present(quantities)) quantities(ind_out_list) = concentration_flag
        if (present(if3d)) if3d(ind_out_list) = .true.
      end do
    end if
      
    if (fu_have_emission_xy(fu_mode(control))) then
      species_emission => fu_species_emission_ctrl(control)
      do ind_species = 1, size(species_emission)
        ind_out_list = ind_out_list + 1
        if (fu_fails(ind_out_list <= max_out_list, 'One of output arguments too small_xy', subname)) return
        if (present(species_list)) species_list(ind_out_list) = species_emission(ind_species)
        if (present(quantities)) quantities(ind_out_list) = emission_scaling_flag
        if (present(if3d)) if3d(ind_out_list) = .false.
      end do
    end if

    num_vars = ind_out_list

  end subroutine get_fieldvars_for_control

  !*************************************************************************
  
  subroutine get_non_fieldvars_for_control(control, num_vars, names)
    ! As above, for non-field variables.
    implicit none
    type(da_control), intent(in) :: control
    integer, intent(out) :: num_vars
    character(len=*), dimension(:), optional, intent(out) :: names
    
    character(len=*), parameter :: subname = 'get_non_fieldvars_for_control'

    num_vars = 0
    
    if (fu_have_emission_zt(fu_mode(control))) then
      num_vars = num_vars + 1
      if (present(names)) then
        if (fu_fails(size(names) >= num_vars, 'Names not big enough', subname)) return
        names(num_vars) = 'emission_zt'
      end if
    end if
    if (fu_have_emission_volc(fu_mode(control))) then
      num_vars = num_vars + 1
      if (present(names)) then
        if (fu_fails(size(names) >= num_vars, 'Names not big enough', subname)) return
        names(num_vars) = 'emission_volc'
      end if
    end if

  end subroutine get_non_fieldvars_for_control


  !*************************************************************************
  
  subroutine control_non_fields_to_files(control, file_units, darules, step_header)
    implicit none
    type(da_control), intent(in) :: control
    integer, dimension(:), intent(in) :: file_units
    type(da_rules), intent(in) :: darules
    character(len=*), intent(in), optional :: step_header

    character(len=*), parameter :: subname = 'control_non_fields_to_files'
    integer :: num_vars, ind_file
    integer :: iostat
    real, dimension(:), pointer :: values

    if (fu_fails(fu_in_physical_space(control), 'Control not in physical space', subname)) return
    
    call get_non_fieldvars_for_control(control, num_vars)
    if (error) return
    if (fu_fails(num_vars > 0, 'No non-field variables in control', subname)) return
    if (fu_fails(size(file_units) >= num_vars, 'file_units not big enough', subname)) return
    
    ! Note same order as in get_nonfieldvars_for_control
    
    ind_file = 0
    if (fu_have_emission_zt(fu_mode(control))) then
      ind_file = ind_file + 1
      if (present(step_header)) write(file_units(ind_file), fmt=*) trim(step_header)
      values => fu_emission_ptr(control)
      call emis_time_hgt_to_file(values, dispersion_vertical, fu_species_emission_ctrl(control), darules, &
                               & file_units(ind_file))
      if (error) return
    end if

    if (fu_have_emission_volc(fu_mode(control))) then
      ind_file = ind_file + 1
      if (present(step_header)) write(file_units(ind_file), fmt=*) trim(step_header)
      values => fu_emission_ptr(control)
      write(file_units(ind_file), fmt=*, iostat=iostat) values 
      call msg('volc_emission smoother output:', values)
      call msg('was written to file', file_units)
      if (fu_fails(iostat == 0, 'failed to write control to file', subname)) return
      ! call set_error('Not implemented', subname)
      return
    end if
    
  end subroutine control_non_fields_to_files


  !*************************************************************************
  
  subroutine control_to_fieldset_get_sizes(control, num_flds, fld_size)
    ! 
    ! Get the number and size of field-type variables from this control. Use for
    ! allocating the arguments to control_to_fieldset.
    ! 
    implicit none
    type(da_control), intent(in) :: control
    integer, intent(out) :: num_flds, fld_size
    
    character(len=*), parameter :: subname = 'control_to_fieldset_get_sizes'
    !type(silam_species), dimension(:), pointer :: species_list
    integer :: num_out_vars, num_levs, nx, ny, ind_var, stat
    integer, parameter :: max_variables = 10000 ! hopefully enough
    logical, dimension(:), allocatable :: if3d

    allocate(if3d(max_variables), stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', subname)) return
    
    if (fu_fails(fu_in_physical_space(control), 'Control not in physical space', subname)) return

    call get_fieldvars_for_control(control, num_out_vars, if3d=if3d)
    if (error) return
    if (num_out_vars == 0) then
      ! Only non-field variables.
      num_flds = 0
      fld_size = int_missing
      return
    end if

    num_levs = fu_NbrOfLevels(dispersion_vertical)
    num_flds = 0
    do ind_var = 1, num_out_vars
      if (if3d(ind_var)) then 
        num_flds = num_flds + num_levs
      else
        num_flds = num_flds + 1
      end if
    end do
    deallocate(if3d)

    call grid_dimensions(control%grid, nx, ny)
    fld_size = nx*ny
      
  end subroutine control_to_fieldset_get_sizes


  !*************************************************************************
  
  subroutine control_to_fieldset(control, fids, data2d, valid_time)
    ! 
    ! Creates a set of fields representing the field data in this control. Valid time is
    ! given, since the control has no timestamp. The output arguments need to be
    ! pre-allocated, the subroutine above can be used to inquire the sizes.
    ! 
    implicit none
    type(da_control), intent(in) :: control
    type(silja_field_id), dimension(:), intent(out) :: fids
    real, dimension(:,:), intent(out) :: data2d
    type(silja_time), intent(in) :: valid_time

    integer :: num_flds, fld_size, ilev, isp, ind_fld, num_species, num_levs, ind_end, ind_start, step
    character(len=*), parameter :: subname = 'control_to_fieldset'
    type(silam_species), dimension(:), pointer :: species_list
    real, dimension(:), pointer :: dataptr, data1d

    call control_to_fieldset_get_sizes(control, num_flds, fld_size)
    if (error) return

    if (num_flds == 0) then
      call set_error('This control has no field output variables', subname)
      return
    end if

    if (fu_fails(num_flds <= size(fids), 'fids too small', subname)) return
    if (fu_fails(num_flds <= size(data2d, 2), 'data2d 2nd dimension too small', subname)) return
    if (fu_fails(fld_size <= size(data2d, 1), 'data2d 1st dimension too small', subname)) return

    ind_fld = 0

    if (fu_have_initial(fu_mode(control))) then
      species_list => fu_species_initial(control)
      num_species = size(species_list)
      num_levs = fu_NbrOfLevels(dispersion_vertical)
      dataptr => fu_initial_ptr(control)
      if (error) return
      do isp = 1, num_species
        do ilev = 1, num_levs
          ind_fld = ind_fld + 1
          ind_start = isp + (ilev-1)*num_species
          step = num_levs * num_species
          ind_end = ind_start + (fld_size-1)*step
          data1d => dataptr(ind_start:ind_end:step)
          if (fu_fails(size(data2d, 1) >= size(data1d), 'data2d wrong size', subname)) return
          data2d(1:size(data1d), ind_fld) = data1d
          fids(ind_fld) = fu_set_field_id(met_src_missing, concentration_flag, &
                                        & valid_time, zero_interval, &
                                        & control%grid, fu_level(dispersion_vertical, ilev), &
                                        & species=species_list(isp))

        end do
      end do
    end if

    if (fu_have_emission_xy(fu_mode(control))) then
      species_list => fu_species_emission_ctrl(control)
      num_species = size(species_list)
      dataptr => fu_emission_ptr(control)
      if (error) return
      do isp = 1, num_species
        ind_fld = ind_fld + 1
        ind_start = isp
        step = num_species
        ind_end = ind_start + (fld_size-1)*step
        data1d => dataptr(ind_start:ind_end:step)
        if (fu_fails(size(data2d, 1) >= size(data1d), 'data2d wrong size', subname)) return
        data2d(1:size(data1d),ind_fld) = data1d
        fids(ind_fld) = fu_set_field_id(met_src_missing, emission_scaling_flag, &
                                      & valid_time, zero_interval, &
                                      & control%grid, surface_level, &
                                      & species=species_list(isp))
        
      end do
    end if
    
  end subroutine control_to_fieldset

  
  !***************************************************************
  
  integer function fu_num_params_volc()
     fu_num_params_volc=7
  end function  fu_num_params_volc

end module da_common
