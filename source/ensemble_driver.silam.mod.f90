!
! Driver module for ensemble data assimilation.
! 
! Currently runs the EnKF with a scheme similar to the 3DVar cycling in da_driver.
! 
! Each ensemble member is an indenpendently running MPI task (MPI must be enabled to use
! EnKF). At analysis times, the state vector is collected into the t_ensemble type, which
! collectively form the ensemble.  Two parallel options for analysis: 1. collect the full
! ensemble into task 0 2. transpose the ensemble: each task analyses a chunk of
! gridpoints.

! In practice, option 2 seems to work fine. It does imply that each member needs to
! allocate space worth one more state vector, and for the full ensemble of observations. 

! There's some confusing naming. When the ensemble is distributed, t_ensemble, ens, etc
! tend to refer to the data of a single member (the ensemble is formed cross the
! processed). When the ensemble is in collected form, ens means the full ensemble. 
! 
! This module sets the rules for perturbations. The perturbations themselves will be
! generated and applied by the perturbations module.

module ensemble_driver
  
  !use da_interface
  use da_driver!, only : update_da_outdir, control_to_file
!  use silam_mpi
!  use silam_partitioning
!  use perturbations
  use enkf
!  use toolbox

  implicit none
  
  private

  public run_enkf
  public run_enks
  public get_chunk
  public transpose_ensemble_in
  public transpose_ensemble_out

  interface defined
     module procedure fu_ens_defined
  end interface

  type t_ensemble
     integer :: ens_size = int_missing
     integer :: state_size = int_missing
     integer :: num_steps = int_missing

     type(da_control) :: state
     type(da_control), dimension(:), pointer :: states => null()
     
     real, dimension(:), pointer :: full_state => null()
     
     logical :: defined = .false.
  end type t_ensemble

  
  CONTAINS

  !************************************************************************************
  
  logical function fu_ens_defined(ens)
    type(t_ensemble), intent(in) :: ens
    fu_ens_defined = ens%defined
  end function fu_ens_defined

  
  !**********************************************************************************
  
  subroutine get_ens_mean(ens, simrules)
    implicit none
    type(t_ensemble), intent(in) :: ens
    type(general_dispersion_rules), intent(in) :: simrules

    integer :: state_size, stat
    type(da_control) :: ct_mean
    real, dimension(:), pointer :: p_mean, p_state, p_mean_chunk
    logical :: master, reduce_ok
    integer :: ind_start, ind_end, remains
    integer, parameter :: chunk_size = 50000000

    if (fu_fails(defined(ens), 'ens not defined', 'get_ens_mean')) return
    state_size = fu_size(ens%state)

    master = smpi_ens_rank == 0
    if (master) then
      call copy_control(ens%state, ct_mean, if_allocate=.true.)
      if (error) return
      p_mean => fu_values_ptr(ct_mean)
      p_mean = 0.0
    else
      p_mean => fu_work_array() ! dummy only
    end if
    
    p_state => fu_values_ptr(ens%state)
    
    ! some memory problems arise with doing the whole reduction at once
    ! ...probably that was just a bug. But now this is here anyway. 
    call smpi_reduce_add(p_state, p_mean, &
                       & 0, smpi_enkf_comm, reduce_ok)
    if (fu_fails(reduce_ok, 'Error with smpi_reduce_add', 'get_ens_mean')) return
    if (master) then
      p_mean = p_mean / ens%ens_size
      call control_to_file(ct_mean, control_missing, simrules%darules, 'mean.grads')
    end if
    call smpi_global_barrier()
    if (master) then
      call destroy(ct_mean)
    else
      call free_work_array(p_mean)
    end if
    
  end subroutine get_ens_mean

  
  !**********************************************************************************

  subroutine get_ens_mean_and_stdev(ens, simrules, filename_base)
    implicit none
    type(t_ensemble), intent(in) :: ens
    type(general_dispersion_rules), intent(in) :: simrules
    character(len=*), intent(in) :: filename_base

    integer :: state_size, stat
    type(da_control) :: ct_mean, ct_stdev
    real, dimension(:), pointer :: p_mean, p_stdev, p_state
    logical :: master, reduce_ok
    integer :: ind_start, ind_end, remains
    integer, parameter :: chunk_size = 50000000
    character(len=*), parameter :: sub_name = 'get_ens_mean_and_stdev'

    if (fu_fails(defined(ens), 'ens not defined', 'get_ens_mean')) return
    call start_count('ens_mean_and_stdev')
    state_size = fu_size(ens%state)

    master = smpi_ens_rank == 0
    call copy_control(ens%state, ct_mean, if_allocate=.true.)
    if (error) return
    p_mean => fu_values_ptr(ct_mean)
    p_mean = 0.0
        
    p_state => fu_values_ptr(ens%state)
    
    call smpi_reduce_add(p_state, p_mean, &
                       & 0, smpi_enkf_comm, reduce_ok)
    if (fu_fails(reduce_ok, 'Error with smpi_reduce_add', sub_name)) return
    if (master) then
      p_mean = p_mean / ens%ens_size
      call control_to_file(ct_mean, control_missing, simrules%darules, trim(filename_base) // '_mean')
    end if
    call smpi_global_barrier()
#ifdef SILAM_MPI
    call mpi_bcast(p_mean(1), size(p_mean), smpi_real_type, 0, smpi_enkf_comm, stat)
    if (fu_fails(stat == 0, 'mpi_bcast failed', sub_name)) return
#endif    
    ! p_mean -> (state - mean)**2
    p_mean = (p_mean - p_state)**2
    if (master) then
      call copy_control(ct_mean, ct_stdev, if_allocate=.true.)
      if (error) return
      p_stdev => fu_values_ptr(ct_stdev)
      p_stdev = 0.0
      call smpi_reduce_add(p_mean, p_stdev, &
                         & 0, smpi_enkf_comm, reduce_ok)
      if (fu_fails(reduce_ok, 'Error with smpi_reduce_add', sub_name)) return
      p_stdev = p_stdev / (ens%ens_size-1)
      p_stdev = sqrt(p_stdev)
      call control_to_file(ct_stdev, control_missing, simrules%darules, trim(filename_base) // '_stdev')
    else
      call smpi_reduce_add(p_mean, p_state, &
                         & 0, smpi_enkf_comm, reduce_ok)
      if (fu_fails(reduce_ok, 'Error with smpi_reduce_add', sub_name)) return
    end if
    
    if (master) then
      call destroy(ct_stdev)
      call destroy(ct_mean)
    else
      call destroy(ct_mean)
    end if
    
    call smpi_global_barrier()
    call stop_count('ens_mean_and_stdev')

  end subroutine get_ens_mean_and_stdev
  
  
  !**********************************************************************************

  subroutine get_ens_mean_v2(ens, simrules)
    implicit none
    type(t_ensemble), intent(in) :: ens
    type(general_dispersion_rules), intent(in) :: simrules

    integer :: state_size, stat
    type(da_control) :: ct_mean
    real, dimension(:), pointer :: p_mean, p_state, p_mean_chunk
    logical :: master, reduce_ok
    integer :: ind_start, ind_end, remains
    integer, parameter :: chunk_size = 50000000

    if (fu_fails(defined(ens), 'ens not defined', 'get_ens_mean')) return
    state_size = fu_size(ens%state)

    master = smpi_ens_rank == 0
    if (master) then
      call copy_control(ens%state, ct_mean, if_allocate=.true.)
      if (error) return
      p_mean => fu_values_ptr(ct_mean)
      p_mean = 0.0
    else
      p_mean_chunk => fu_work_array() ! dummy only
    end if
    
    p_state => fu_values_ptr(ens%state)
    
    ! some memory problems arise with doing the whole reduction at once
    
    remains = size(p_state)
    ind_start = 1
    do while(remains > 0)
      ind_end = min(ind_start + chunk_size - 1, size(p_state))
      call msg('remains, ind_end', remains, ind_end)
      if (master) then 
        p_mean_chunk => p_mean(ind_start:ind_end)
      end if
      call smpi_reduce_add(p_state(ind_start:ind_end), p_mean_chunk, &
                         & 0, smpi_enkf_comm, reduce_ok)
      if (fu_fails(reduce_ok, 'Error with smpi_reduce_add', 'get_ens_mean')) return
      ind_start = ind_start + chunk_size
      remains = remains - chunk_size
    end do
    if (master) then
      p_mean = p_mean / ens%ens_size
      call control_to_file(ct_mean, control_missing, simrules%darules, 'mean.grads')
    end if
    call smpi_global_barrier()
    if (master) then
      call destroy(ct_mean)
    else
      call free_work_array(p_mean_chunk)
    end if
    
  end subroutine get_ens_mean_v2

  
  !************************************************************************************

  subroutine propagate_ensemble(ens, model)
    implicit none
    type(t_ensemble), intent(inout) :: ens
    type(model_container), intent(inout) :: model

    integer :: ind_ens
    type(da_control) :: control
    real, dimension(:), pointer :: p_control
    character(len=*), parameter :: sub_name = 'collect_ensemble'
 
    if (fu_fails(ens%ens_size > 0, 'ens%size < 1', sub_name)) return

    ! call control2perturb(control, model%perturbation)
    
    call msg('ENKF: run model')
    call run_forecast_from_control(ens%state, model)
    
    call vector_from_model(ens%state, model, handle_negatives=.false.)
    call msg('ENKF: run model done')
      
    ! call control_from_perturb(model%perturbation, control)

  end subroutine propagate_ensemble

  
  !************************************************************************************

  subroutine run_enkf(model)
    implicit none
    type(model_container), intent(inout) :: model
    !type(da_control) :: analysis, background
    !real :: cost_obs, cost
    !real, dimension(:), pointer :: obs_values, obs_variance
    !real, dimension(2) :: cost_bgr
    !character(len=128) :: log_str
    type(da_rules) :: darules, darules_missing
    type(silam_pollution_cloud), pointer :: cloud
    !character(len=fnlen) :: analysis_file_name
    integer :: n_obs_values
    
    type(model_container) :: model_integr, & ! the model for forward integration, no DA
         & model_3dvar, &                    ! the dummy model for 3dvar minimization 
         & model_3dvar_adj                   ! as above but "adjoint"
    integer :: stat, chunk_size, chunk_offset
!    type(model_container) :: model_fwd, model_adj
!    type(silam_pollution_cloud), pointer :: cloud_adj
    type(Tsilam_namelist) :: nl_dummy
    type(observationPointers), target :: obs_ptr_missing
    type(general_dispersion_rules), target :: simrules_integr, simrules_3d
    type(da_rules) :: da_rules_step
    type(silja_time) :: now, end_time, time_forecast_end
    logical :: last_forecast
    character(len=*), parameter :: sub_name = 'run_enkf'
    type(t_ensemble) :: ens
    real, dimension(:,:), allocatable, save :: ens_full
    integer, parameter :: count_max_perturbations = 10
    logical :: fatal

    if (.not. smpi_is_mpi_version()) then
      call set_error('MPI needed for EnKF', sub_name)
      return
    end if

    call set_emis_proc_by_ctrl(model%cloud, model%rules, model%rules%darules%perturbVariable)
    if (error) return

    cloud => model%cloud
    darules = model%rules%darules
    !model%rules%if_collect_linearization = .true.
    !model%rules%if_collect_linearization = .false.
    if (fu_fails(darules%method == flag_enkf, 'Method not enkf', sub_name)) return

    model_integr = model
    simrules_integr = model%rules
    model_integr%rules => simrules_integr ! no longer pointing to the original rules
    model_integr%obs_ptr => obs_ptr_missing

    ! Take period to compute from assimilation interval, start time updated after each
    ! analysis.
    model_integr%rules%periodToCompute = model%rules%darules%interval_btwn_assim
    model_3dvar = model
    simrules_3d = model%rules
    model_3dvar%rules => simrules_3d
    
    if (error) return
    
    now = model%rules%startTime
    if (now < model%rules%darules%first_analysis_time) then
      call set_error('doesn''t work yet', sub_name)
      return
    end if

    end_time = model%rules%startTime + model%rules%periodToCompute

    call init_ensemble(ens, model%cloud, darules)
    
    ! Allocate work storage depending on whether the filter is run master-only or with
    ! mpi.
    call msg('ENKF: state and ensemble sizes:', ens%state_size, ens%ens_size)
    if (((darules%enkf_master_only .and. smpi_ens_rank == 0)  .and. .not. allocated(ens_full))) then
      call msg('ENKF: allocating work arrays for master-only mode')
      allocate(ens_full(ens%state_size, ens%ens_size), stat=stat)
    else if (.not. darules%enkf_master_only .and. .not. allocated(ens_full)) then
      call msg('ENKF: allocating work arrays for the MPI mode')
      call get_enkf_chunk_size_offset(ens%state_size, ens%ens_size, chunk_size, chunk_offset)
      call msg('ENKF: my chunk size and offset:', chunk_size, chunk_offset)
      allocate(ens_full(chunk_size, ens%ens_size), stat=stat)
    else if (.not. allocated(ens_full)) then
      ! dummy array, non-master process with master-only analysis
      allocate(ens_full(1,1), stat=stat)
    end if
    if (fu_fails(stat == 0, 'Allocating full ensemble failed', 'run_enkf')) return

    allocate(model_integr%perturbations(count_max_perturbations), stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', 'run_enkf')) return
    call set_perturbations(model_integr%perturbations, &
                         & fu_species_emission(model%cloud), &
                         & fu_species_transport(model%cloud), &
                         & darules, dispersion_grid, dispersion_vertical, now)
    if (error) return
    call init_random_seed(smpi_global_rank)
    simrules_integr%if_finalize = .false.

    do while (now < end_time)
      if (model%rules%darules%output_level > da_out_none_flag) call update_enkf_outdir(model, simrules_3d%darules, now)
      call smpi_global_barrier()
      if (error) return

      call start_count('EnKF analysis total')
      if (now >= model%rules%darules%first_analysis_time) then
        call msg('ENKF: analysis time:' // trim(fu_str(now)))
        if (darules%enkf_master_only) then
          call assimilate_master_only(now, ens, fatal)
        else
          call assimilate_mpi(now, ens, fatal)
        end if
        if (fatal) return
        if (error) then
          call msg('ENKF: assimilation failed, ensemble not changed')
          call unset_error('run_enkf')
        else
          call msg('ENKF: assimilation done')
        end if
        call msg('ENKF: ===================')
      else
        call msg('ENKF: spinup, skip analysis')
      end if
      call smpi_global_barrier()
      call stop_count('EnKF analysis total')
      
      model_integr%rules%startTime = now
      time_forecast_end =  now + model_integr%rules%periodToCompute 
      if (time_forecast_end > darules%last_analysis_time) then
        call msg('ENKF: Last analysis done')
        time_forecast_end = end_time - model_integr%rules%timestep
      end if
      last_forecast = (time_forecast_end + model_integr%rules%timestep >= end_time)
      if (last_forecast) then 
        call msg('ENKF: Will do final forecast')
        ! Make the last integration to go until end_time to get output there.
        model_integr%rules%periodToCompute = end_time - now 
        simrules_integr%if_finalize = .true.
        !model_integr%rules%periodToCompute + model_integr%rules%timestep
      end if

      call propagate_ensemble(ens, model_integr)
      now = now + model_integr%rules%periodToCompute + model_integr%rules%timestep
    end do

    ! the analysis step

    da_rules_step = model_3dvar%rules%darules
  
  contains
    
    !======================================================================
  
    subroutine assimilate_master_only(now, ens, fatal_error)
      implicit none
      type(silja_time), intent(in) :: now
      type(t_ensemble), intent(inout) :: ens
      logical, intent(out) :: fatal_error

      real, dimension(:), pointer :: mdl_obs_val, obs_data, obs_var
      real, dimension(:,:), allocatable :: obs_loc
      integer :: iPurpose
      integer, dimension(2) :: obs_size
      character(len=*), parameter :: sub_name = 'assimilate_master_only'
      real, dimension(:,:), allocatable, save :: model_localisation
      real, dimension(:,:), allocatable :: obs_localisation, mdl_obs_ens
      character(len=12) :: time_str
      type(observationPointers), save :: obs_ptr
      logical :: scatter_ok
      type(DA_rules), pointer :: pDA_rules

      ! Normal errors can be unset and the integration proceeds with the ensemble
      ! unchanged. However, if scatter fails, the run must crash.
      fatal_error = .false.

      model_3dvar%rules%darules%da_begin = now
      
      pDA_rules => model_3dvar%rules%darules
      call set_observations(pDA_rules%station_path, &
                          & pDA_rules%DA_begin, pDA_rules%assim_window, pDA_rules%obs_files, &
                          & pDA_rules%num_obs_files_assim, pDA_rules%num_obs_files_eval, &
                          & fu_species_transport(cloud), &
                          & fu_species_optical(cloud), model_3dvar%obs_ptr)
      fatal_error = .not. sync_errors()
      if (fatal_error .or. error) then
        call msg('ENKF: failed with observations, returning...')
        return
      end if
      !call set_error('crash!', sub_name)
      !call destroy_observations(model_3dvar%obs_ptr)
      !return
            
      !if (sum(model_3dvar%obs_ptr%obs_size) == 0) then
      !  call set_error('ENKF: no observations', sub_name)
      !end if
      fatal_error = .not. sync_errors()
      if (fatal_error .or. error) return

      obs_size = model_3dvar%obs_ptr%obs_size
      mdl_obs_val => fu_work_array(sum(obs_size))
      obs_var => fu_work_array(sum(obs_size))
      
      if (smpi_ens_rank == 0) then
        allocate(mdl_obs_ens(sum(obs_size), ens%ens_size), obs_localisation(2, sum(obs_size)), stat=stat)
        if (fu_fails(stat == 0, 'Allocate master obs arrays failed', sub_name)) continue
      else
        ! dummy allocation for mpi_gather
        allocate(mdl_obs_ens(1,1), stat=stat)
        if (fu_fails(stat == 0, 'Allocate master obs arrays failed', sub_name)) continue
      end if
      fatal_error = .not. sync_errors()
      if (error .or. fatal_error) return

      call msg('ENKF: evaluate observations')
      call forward_3d(model_3dvar)
      if (error) return

      call collect_model_data(model_3dvar%obs_ptr, mdl_obs_val) !, obs_size)
      call msg('ENKF: gather ensemble')
      call start_count('ensemble_scatter_gather')
      call gather_ensemble(ens, ens_full, mdl_obs_val, mdl_obs_ens, obs_size)
      call stop_count('ensemble_scatter_gather')
      if (error) return
      
      call get_obs_data(model_3dvar%obs_ptr, obs_data, obs_size)
      call collect_variance(model_3dvar%obs_ptr, obs_var) !, obs_size)
      if (simrules_3d%darules%output_level > da_out_none_flag) then
        call msg('Dumping FC observations...')
        write(time_str, fmt='(I4,I2.2,I2.2,I2.2,I2.2)') fu_year(now), fu_mon(now), fu_day(now), &
             & fu_hour(now), fu_min(now)
        do iPurpose = 1, 2
          call dump_observations(model_3dvar%obs_ptr, fu_species_transport(model%cloud), &
                               & simrules_3d%darules%outputdir + dir_slash + 'obs_fc_' &
                               & + time_str + chObsPurpose(iPurpose), iPurpose)
        end do
        call msg('Done')
      end if
      fatal_error = .not. sync_errors()
      if (fatal_error .or. error) return

      if (smpi_ens_rank == 0) call get_localisation(model_3dvar%obs_ptr, obs_localisation)
      if (.not. allocated(model_localisation)) then
        allocate(model_localisation(2, ens%state_size), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed (model loc)', sub_name)) return
        call get_localisation(ens%state, model_localisation, darules)
      end if
      
      !if (simrules_3d%darules%output_level > da_out_first_last_flag) then
      !  call msg('ENKF: call get_ens_mean')
      !  call get_ens_mean_and_stdev(ens, model_3dvar%rules, trim&
      !                            & (model_3dvar%rules%darules%outputdir) // dir_slash // 'fc')
      !end if
      if (smpi_ens_rank == 0 .and. .not. error) then
        call msg('ENKF: update')
        call enkf_update(ens_full, mdl_obs_ens, obs_data(1:obs_size(1)), obs_var(1:obs_size(1)), &
                       & obs_localisation, model_localisation, &
                       & model%rules%darules%loc_dist_m, &
                       & model%rules%darules%loc_type, &
                       & model%rules%darules%enkf_flavor, &
                       & model%rules%darules%rfactor, &
                       & diagn_out_dir=char_missing)
      end if
      fatal_error = .not. sync_errors()
      if (fatal_error) return

      ! error is now synced with all tasks
      if (.not. error) then
        if (simrules_3d%darules%output_level > da_out_first_last_flag) then
          call forward_3d(model_3dvar)
          do iPurpose = 1, 2
            call dump_observations(model_3dvar%obs_ptr, fu_species_transport(model%cloud), &
                                 & simrules_3d%darules%outputdir + dir_slash + 'obs_an_' &
                                     & + fu_str(smpi_ens_rank) + chObsPurpose(iPurpose), iPurpose)
          end do
          !call get_ens_mean_and_stdev(ens, model_3dvar%rules, &
          !                          & trim(model_3dvar%rules%darules%outputdir) // dir_slash // 'an')
        end if
        call msg('ENKF: scatter ensemble')
        call start_count('ensemble_scatter_gather')
        call scatter_ensemble(ens, ens_full, scatter_ok)
        call stop_count('ensemble_scatter_gather')
        fatal_error = .not.scatter_ok
      end if

      ! cleanup 
      call free_work_array(mdl_obs_val)
      call free_work_array(obs_var)
      if (allocated(model_localisation)) deallocate(model_localisation)
      if (allocated(obs_localisation)) deallocate(obs_localisation)
      if (allocated(mdl_obs_ens)) deallocate(mdl_obs_ens)
      call destroy_observations(model_3dvar%obs_ptr)
        
    end subroutine assimilate_master_only

    !======================================================================

    subroutine assimilate_mpi(now, ens, fatal_error)
      implicit none
      type(silja_time), intent(in) :: now
      type(t_ensemble), intent(inout) :: ens
      logical, intent(out) :: fatal_error

      real, dimension(:), pointer :: mdl_obs_val_single, obs_data, obs_var, ptr_val
      real, dimension(:,:), allocatable :: obs_loc, mdl_obs_val_all
      integer :: my_size, my_offset, iPurpose
      character(len=*), parameter :: sub_name = 'assimilate_master_only'
      real, dimension(:,:), allocatable, save :: model_localisation_full, model_localisation_chunk
      real, dimension(:,:), allocatable :: obs_localisation, mdl_obs_ens
      character(len=12) :: time_str
      type(observationPointers), save :: obs_ptr
      logical :: scatter_ok
      integer, dimension(:), pointer :: offsets, sizes
      type(DA_rules), pointer :: pDA_rules
      integer, dimension(2) :: obs_size

      call msg('ENKF: mpi analysis')

      ! Normal errors can be unset and the integration proceeds with the ensemble
      ! unchanged. However, if scatter fails, the run must crash.
      fatal_error = .false.

      model_3dvar%rules%darules%da_begin = now
      
      call start_count('EnKF load observations')
      pDA_rules => model_3dvar%rules%darules
      call set_observations(pDA_rules%station_path, &
                          & pDA_rules%DA_begin, pDA_rules%assim_window, pDA_rules%obs_files, &
                          & pDA_rules%num_obs_files_assim, pDA_rules%num_obs_files_eval, &
                          & fu_species_transport(cloud), &
                          & fu_species_optical(cloud), model_3dvar%obs_ptr)
      fatal_error = .not. sync_errors()
      if (fatal_error .or. error) then
        call msg('ENKF: failed with observations, returning...')
        return
      end if
      call stop_count('EnKF load observations')
            
      !if (sum(model_3dvar%obs_ptr%obs_size) == 0) then
      !  call set_error('ENKF: no observations', sub_name)
      !end if
      fatal_error = .not. sync_errors()
      if (fatal_error .or. error) return

      obs_size = model_3dvar%obs_ptr%obs_size
      mdl_obs_val_single => fu_work_array(sum(obs_size))
      obs_var => fu_work_array(sum(obs_size))
      
      allocate(mdl_obs_ens(sum(obs_size), ens%ens_size), obs_localisation(2, sum(obs_size)), stat=stat)
      if (fu_fails(stat == 0, 'Allocate master obs arrays failed', sub_name)) continue
      
      fatal_error = .not. sync_errors()
      if (error .or. fatal_error) return

      call msg('ENKF: evaluate observations')
      call forward_3d(model_3dvar)
      if (error) return

      call start_count('EnKF evaluate observations')
      call collect_model_data(model_3dvar%obs_ptr, mdl_obs_val_single) !, obs_size)

      call get_localisation(model_3dvar%obs_ptr, obs_localisation)
      call stop_count('EnKF evaluate observations')

      call msg('ENKF: gather ensemble')
      call start_count('ensemble_transpose')
      call transpose_ensemble_in(fu_values_ptr(ens%state), ens_full, ens%ens_size, ens%state_size, &
                               & mdl_obs_val_single(1:sum(obs_size)), mdl_obs_ens, sum(obs_size))
      call stop_count('ensemble_transpose')
      if (error) return
      
      call get_obs_data(model_3dvar%obs_ptr, obs_data, obs_size)
      call collect_variance(model_3dvar%obs_ptr, obs_var) !, obs_size)
      if (simrules_3d%darules%output_level > da_out_none_flag) then
        call msg('Dumping FC observations...')
        write(time_str, fmt='(I4,I2.2,I2.2,I2.2,I2.2)') fu_year(now), fu_mon(now), fu_day(now), &
             & fu_hour(now), fu_min(now)
        do iPurpose = 1, 2
          call dump_observations(model_3dvar%obs_ptr, fu_species_transport(model%cloud), &
                               & simrules_3d%darules%outputdir + dir_slash + 'obs_fc_' &
                               & + trim(time_str) + chObsPurpose(iPurpose), iPurpose)
        end do
        call msg('Done')
      end if
      fatal_error = .not. sync_errors()
      if (fatal_error .or. error) return
      
      if (.not. allocated(model_localisation_chunk)) then
        ! do this only once, model localisation and chunking won't change.
        ! one-based offsets
        offsets => fu_work_int_array()
        sizes => fu_work_int_array()
        call get_chunk(ens%state_size, ens%ens_size, offsets, sizes, if_zero_based=.false.)
        if (error) return
        my_offset = offsets(smpi_ens_rank+1)
        my_size = sizes(smpi_ens_rank+1)
        call free_work_array(offsets)
        call free_work_array(sizes)
        allocate(model_localisation_full(2, ens%state_size), &
               & model_localisation_chunk(2, my_size), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed (model loc)', sub_name)) return
        call get_localisation(ens%state, model_localisation_full, darules)
        model_localisation_chunk  = model_localisation_full(:, my_offset: my_offset+my_size-1)
        deallocate(model_localisation_full)
      end if
      
      !if (simrules_3d%darules%output_level > da_out_first_last_flag) then
      !  call msg('ENKF: call get_ens_mean')
      !  call get_ens_mean_and_stdev(ens, model_3dvar%rules, trim&
      !                            & (model_3dvar%rules%darules%outputdir) // dir_slash // 'fc')
      !end if
      if (.not. error) then
        call msg('ENKF: update')
        !
        ! Guess, for the ensemble update one needs only assimilated observations
        !
        call enkf_update(ens_full, mdl_obs_ens, obs_data(1:obs_size(1)), obs_var(1:obs_size(1)), &
                       & obs_localisation, model_localisation_chunk, &
                       & model%rules%darules%loc_dist_m, &
                       & model%rules%darules%loc_type, &
                       & model%rules%darules%enkf_flavor, &
                       & model%rules%darules%rfactor, &
                       & diagn_out_dir=char_missing)

      end if
      fatal_error = .not. sync_errors()
      if (fatal_error) return

      ! error is now synced with all tasks
      if (.not. error) then
        call msg('ENKF: scatter ensemble')
        call start_count('ensemble_transpose')

        ptr_val => fu_values_ptr(ens%state)

        call transpose_ensemble_out(ptr_val, ens_full, ens%ens_size, &
                                  & ens%state_size, scatter_ok)
        call stop_count('ensemble_transpose')
        fatal_error = .not.scatter_ok
        if (simrules_3d%darules%output_level > da_out_first_last_flag) then
          call forward_3d(model_3dvar)
          do iPurpose = 1, 2
            call dump_observations(model_3dvar%obs_ptr, fu_species_transport(model%cloud), &
                                 & simrules_3d%darules%outputdir + dir_slash + 'obs_an_' &
                                 & + fu_str(smpi_ens_rank) + chObsPurpose(iPurpose), iPurpose)
          enddo
          call get_ens_mean_and_stdev(ens, model_3dvar%rules, &
                                    & trim(model_3dvar%rules%darules%outputdir) // dir_slash // 'an')
        end if
      end if

      ! cleanup 
      call free_work_array(mdl_obs_val_single)
      call free_work_array(obs_var)
      !if (allocated(model_localisation)) deallocate(model_localisation_full)
      if (allocated(obs_localisation)) deallocate(obs_localisation)
      if (allocated(mdl_obs_ens)) deallocate(mdl_obs_ens)
      call destroy_observations(model_3dvar%obs_ptr)
        
    end subroutine assimilate_mpi
    
  end subroutine run_enkf
  
  
  !************************************************************************************

  logical function sync_errors() result(ok)
    implicit none
    integer :: stat
    logical :: any_error

    if (.not. smpi_is_mpi_version()) then
      ok = .true.
      return
    end if
      
#ifdef SILAM_MPI
    call mpi_allreduce(error, any_error, 1, mpi_logical, mpi_lor, smpi_enkf_comm, stat)
    ok = stat == MPI_SUCCESS
    error = any_error
#endif
  end function sync_errors

  
  !************************************************************************************

  subroutine transpose_ensemble_in(mbr_val, ens_chunk, ens_size, mdl_size, &
                                 & mdl_obs_single, mdl_obs_all, obs_size)
    implicit none
    real, dimension(:), intent(in) :: mbr_val ! local, all tasks
    real, dimension(:,:), intent(out) :: ens_chunk ! own chunk sent for each task
    real, dimension(:), intent(in) :: mdl_obs_single ! local, all tasks
    real, dimension(:,:), intent(out) :: mdl_obs_all ! broadcasted to all tasks
    integer, intent(in) :: mdl_size, ens_size, obs_size
    
    integer, dimension(:), pointer :: displs_send, sizes_send, displs_recv, sizes_recv
    integer :: stat, my_size, ii
    character(len=*), parameter :: sub_name = 'transpose_ensemble_in'

    if (.not. smpi_is_mpi_version()) then
      call set_error('Running without MPI', sub_name)
      call unset_error(sub_name)
      return
    end if

    displs_send => fu_work_int_array()
    sizes_send => fu_work_int_array()
    displs_recv => fu_work_int_array()
    sizes_recv => fu_work_int_array()
    if (error) return
    
    call get_chunk(mdl_size, ens_size, displs_send, sizes_send, .true.)
    if (error) return

    if (fu_fails(size(ens_chunk, 1) <= mdl_size, 'ens_chunk too small (mdl)', sub_name)) return
    if (fu_fails(size(ens_chunk, 2) <= ens_size, 'ens_chunk too small (ens)', sub_name)) return
    
    ! Model state
    !
    my_size = sizes_send(smpi_ens_rank+1)
    displs_recv(1:ens_size) = (/(ii*my_size, ii=0, ens_size-1)/)
    sizes_recv(1:ens_size) = my_size
    ! MPI_ALLTOALLV(SENDBUF, SENDCOUNTS, SDISPLS, SENDTYPE,
    !               RECVBUF, RECVCOUNTS, RDISPLS, RECVTYPE, COMM, IERROR)

#ifdef SILAM_MPI
    call mpi_alltoallv(mbr_val(1:), sizes_send(1:), displs_send(1:), smpi_real_type, &
                     & ens_chunk(1:,1), sizes_recv(1:), displs_recv(1:), smpi_real_type, &
                     & smpi_enkf_comm, stat)
    if (fu_fails(stat == MPI_SUCCESS, 'Failed mpi_alltoallv', sub_name)) return
    ! Observations
    !
    ! MPI_ALLGATHER(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT, 
    !               RECVTYPE, COMM, IERROR)
    call mpi_allgather(mdl_obs_single(1), obs_size, smpi_real_type, mdl_obs_all(1,1), obs_size, smpi_real_type, &
                     & smpi_enkf_comm, stat)
    if (fu_fails(stat == MPI_SUCCESS, 'Failed mpi_allgather', sub_name)) return
#endif

    call free_work_array(displs_send)
    call free_work_array(displs_recv)
    call free_work_array(sizes_send)
    call free_work_array(sizes_recv)
        
  end subroutine transpose_ensemble_in

                                 
  !**********************************************************************************

  subroutine transpose_ensemble_out(mbr_val, ens_chunk, ens_size, mdl_size, transp_ok)
    implicit none
    real, dimension(:), intent(out) :: mbr_val ! local, all tasks
    real, dimension(:,:), intent(in) :: ens_chunk ! own chunk sent by each task
    integer, intent(in) :: mdl_size, ens_size
    logical, intent(out) :: transp_ok

    integer, dimension(:), pointer :: displs_send, sizes_send, displs_recv, sizes_recv
    integer :: stat, my_size, ii
    character(len=*), parameter :: sub_name = 'transpose_ensemble_out'

    if (.not. smpi_is_mpi_version()) then
      call set_error('Running without MPI', sub_name)
      call unset_error(sub_name)
      transp_ok = .true.
      return
    end if

    transp_ok = .false.

    displs_send => fu_work_int_array()
    sizes_send => fu_work_int_array()
    displs_recv => fu_work_int_array()
    sizes_recv => fu_work_int_array()
    if (error) return
    
    call get_chunk(mdl_size, ens_size, displs_recv, sizes_recv, .true.)
    if (error) return

    if (fu_fails(size(ens_chunk, 1) <= mdl_size, 'ens_chunk too small (mdl)', sub_name)) return
    if (fu_fails(size(ens_chunk, 2) <= ens_size, 'ens_chunk too small (ens)', sub_name)) return

    my_size = sizes_recv(smpi_ens_rank+1)
    displs_send(1:ens_size) = (/(ii*my_size, ii=0, ens_size-1)/)
    sizes_send(1:ens_size) = my_size

#ifdef SILAM_MPI
    ! send ens_chunk(:,ind_ens) to each receiver
    ! take appropriate chunk size on the receiving side
    call mpi_alltoallv(ens_chunk(1:,1), sizes_send(1:), displs_send(1:), smpi_real_type, &
                     & mbr_val(1:), sizes_recv(1:), displs_recv(1:), smpi_real_type, &
                     & smpi_enkf_comm, stat)
    if (fu_fails(stat == MPI_SUCCESS, 'Failed mpi_alltoallv', sub_name)) return
#endif
    
    call free_work_array(displs_send)
    call free_work_array(displs_recv)
    call free_work_array(sizes_send)
    call free_work_array(sizes_recv)

    transp_ok = .true.
    
  end subroutine transpose_ensemble_out


  !************************************************************************************

  subroutine gather_ensemble(member, ens_full, mdl_obs_single, mdl_obs_all, obs_size)
    ! Collect ensemble state to a 2D array on task 0. Other tasks allocate the array as
    ! zero-sized.
    implicit none
    type(t_ensemble), intent(in) :: member ! local, all tasks
    real, dimension(:,:), intent(out) :: ens_full ! task 0
    real, dimension(:), intent(in) :: mdl_obs_single ! all tasks
    real, dimension(:,:), intent(out) :: mdl_obs_all ! task 0
    integer, dimension(:), intent(in) :: obs_size

    character(len=*), parameter :: sub_name = 'gather_ensemble'
    logical :: is_master
    integer :: ind_task, stat
    real, dimension(:), pointer :: p_values

    is_master = smpi_ens_rank == 0

    if (is_master) then
      if (fu_fails(size(ens_full, 2) == member%ens_size, 'Ensemble sizes don''t match', sub_name)) return
      if (fu_fails(size(ens_full, 1) == member%state_size, 'State sizes don''t match', sub_name)) return
      if (fu_fails(size(mdl_obs_all,1) >= sum(obs_size), 'mdl_obs_all too small', sub_name)) return
      if (fu_fails(size(mdl_obs_all,2) == member%ens_size, 'Ens-obs sizes don''t match', sub_name)) return
    end if

    p_values => fu_values_ptr(member%state)
    if (fu_fails(associated(p_values), 'p_values not associated', sub_name)) return
#ifdef SILAM_MPI
    call mpi_gather(p_values(1), member%state_size, &
                  & smpi_real_type, ens_full(1,1), member%state_size, smpi_real_type, &
                  & 0, smpi_enkf_comm, stat)
    if (fu_fails(stat == MPI_SUCCESS, 'Failed smpi_gather state', sub_name)) return

    call mpi_gather(mdl_obs_single(1), obs_size(1), &
                  & smpi_real_type, mdl_obs_all(1,1), obs_size(1), smpi_real_type, &
                  & 0, smpi_enkf_comm, stat)
    if (fu_fails(stat == MPI_SUCCESS, 'Failed smpi_gather state', sub_name)) return
#endif    
    
  end subroutine gather_ensemble

  
  !**********************************************************************************

  subroutine scatter_ensemble(member, ens_full, success)
    implicit none
    type(t_ensemble), intent(inout) :: member
    real, dimension(:,:), intent(in) :: ens_full
    ! additional success flag so assimilate() can return a "fatal error" flag. If scatter
    ! fails, cannot rely on the members having a valid state, and the run must crash.
    logical :: success

    character(len=*), parameter :: sub_name = 'scatter_ensemble'
    integer :: ind_task, stat
    real, dimension(:), pointer :: p_values

    success = .false.
    
    if (smpi_ens_rank == 0) then
      if (fu_fails(size(ens_full, 2) == member%ens_size, 'Ensemble sizes don''t match', sub_name)) return
      if (fu_fails(size(ens_full, 1) == member%state_size, 'State sizes don''t match', sub_name)) return
    end if

    p_values => fu_values_ptr(member%state)
    if (fu_fails(associated(p_values), 'p_values not associated', sub_name)) return
#ifdef SILAM_MPI
    call mpi_scatter(ens_full(1,1), member%state_size, smpi_real_type, &
                   & p_values(1), member%state_size, smpi_real_type, &
                   & 0, smpi_enkf_comm, stat)
    if (fu_fails(stat == MPI_SUCCESS, 'Failed smpi_gather', sub_name)) return
#endif
    success = .true.

  end subroutine scatter_ensemble

  !************************************************************************************

  subroutine init_ensemble(ens, cloud, rules)
    implicit none
    type(t_ensemble), intent(out) :: ens
    type(silam_pollution_cloud), intent(in) :: cloud
    type(da_rules), intent(in) :: rules
    
    integer :: ind_ens
    integer :: stat

    call init_control_from_cloud(ens%state, cloud, rules, physical_space, initial_value=0.0)
    if (error) return

!    call report(ens%state)
    
    ens%state_size = fu_size(ens%state)
    ens%ens_size = rules%ens_size
    ens%defined = .true.
    
  end subroutine init_ensemble

  
  !**********************************************************************************

  subroutine init_ensemble_multitime(ens, cloud, rules, num_steps)
    implicit none
    type(t_ensemble), intent(out) :: ens
    type(silam_pollution_cloud), intent(in) :: cloud
    type(da_rules), intent(in) :: rules
    integer, intent(in) :: num_steps
    
    integer :: ind_ens
    integer :: stat
    character(len=*), parameter :: subname = 'ini_ensemble_multitime'

    if (fu_fails(num_steps > 0, 'Bad num_steps', subname)) return
    allocate(ens%states(num_steps), stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', subname)) return

    call init_control_array_from_cloud(ens%states, cloud, rules, physical_space, num_steps, ens%full_state, &
                          & initial_value=0.0)
    if (error) return

    call msg('report ens states 1')
    call report(ens%states(1))
    ens%num_steps = num_steps
    ens%state_size = fu_size(ens%states(1)) * num_steps
    ens%ens_size = rules%ens_size
    ens%defined = .true.
    
  end subroutine init_ensemble_multitime

  !************************************************************************************

  subroutine update_enkf_outdir(model, darules, now)
    implicit none
    type(model_container), intent(inout) :: model 
    type(da_rules), intent(inout) :: darules
    type(silja_time), intent(in) :: now

    character(len=fnlen) :: template_string, main_out_dir
    type(grads_template) :: gr_template
    type(silja_interval) :: forecast_length
    character(len=fnlen), save :: prev_out_dir = ''

    ! The 4dvar "iteration" model has no outDef, here it is needed though!
    if (fu_fails(associated(model%out_def), 'model%out_def not associated', 'update_da_outdir')) return
    !if (fu_fails(defined(model%out_def), 'model%out_def not defined', 'update_da_outdir')) return

    forecast_length = now - model%rules%startTime
    
    template_string = darules%outdir_templ_str
    if (fu_fails(template_string /= '', 'Empty output template string', 'update_da_outdir')) return
    ! Get main output directory according to current time
    main_out_dir = fu_FNm(fu_output_template(model%out_def), &
                        & model%rules%startTime, & ! ini_time
                        & model%rules%startTime, & ! anal_time
                        & forecast_length, &
                        & model%rules%chCaseNm, &
                        & fu_get_id_str_from_id_nbr(model%source,1)) ! same as in io_server
    if (error) return
    main_out_dir = fu_dirname(main_out_dir)
    call replace_string(template_string, '%main_output_dir', trim(main_out_dir))
    if (error) return
    call msg('Updating output directory: ' // trim(template_string) // '-', len_trim(template_string))
    call decode_template_string(template_string, gr_template)
    darules%outputdir = fu_FNm(gr_template, &
                           & model%rules%startTime, &
                           & now, &
                           & zero_interval, &
                           & model%rules%chCaseNm, &
                           & fu_get_id_str_from_id_nbr(model%source,1)) ! same as in io_server
    !call da_msg('Current output directory:' // trim(darules%outputdir))
    if (prev_out_dir /= darules%outputdir) then
      call create_directory_tree(darules%outputdir)
      if (error) return
      prev_out_dir = darules%outputdir
    end if

  end subroutine update_enkf_outdir

  
!**********************************************************************************
  
  subroutine run_enks(model)
    !
    ! Ensemble Kalman Smoother, as referred by Evensen (2003). Uses the same EnKF code as
    ! run_enkf, but implements the following extensions:
    ! 
    ! - the observations are evaluated during a forward integration (similar to
    ! 4D-Var). They may depend on multiple timesteps (averaging), and need not be forced
    ! to the analysis time.  
    ! 
    ! - multiple time steps can be included in the state vector. In addition to the
    ! current timestep, the analysis can update previous times based on the newly
    ! evaluated observations. The smoother steps will be written into separate output
    ! files. 

    implicit none
    type(model_container), intent(inout) :: model
    !type(da_control) :: analysis, background
    !real :: cost_obs, cost
    !real, dimension(:), pointer :: obs_values, obs_variance
    !real, dimension(2) :: cost_bgr
    !character(len=128) :: log_str
    type(da_rules) :: darules, darules_missing
    type(silam_pollution_cloud), pointer :: cloud
    !character(len=fnlen) :: analysis_file_name
    integer :: n_obs_values
    
    type(model_container) :: model_integr ! the model for forward integration, no DA
    integer :: stat, chunk_size, chunk_offset, num_smtr_steps
!    type(model_container) :: model_fwd, model_adj
!    type(silam_pollution_cloud), pointer :: cloud_adj
    type(Tsilam_namelist) :: nl_dummy
    type(observationPointers), target :: obs_ptr_missing
    type(general_dispersion_rules), target :: simrules_integr, simrules_3d
    type(da_rules) :: da_rules_step
    type(silja_time) :: now, end_time, time_forecast_end, time_analysis_next, time_get_control_next
    type(silja_interval) :: step_integr, smoother_window
    logical :: last_forecast
    character(len=*), parameter :: sub_name = 'run_enks'
    type(t_ensemble) :: ensemble_multitime
    real, dimension(:,:), allocatable, save :: ens_full
    integer, parameter :: count_max_perturbations = 10
    logical :: fatal
    type(da_control_ptr), dimension(:), allocatable :: smtr_step_ptrs
    type(silja_time), dimension(:), allocatable :: smtr_times
    integer, dimension(:), pointer :: smtr_file_units, smtr_file_units_nfld
    real, dimension(:,:), pointer :: smtr_fieldstorage
    type(Tmass_map), pointer :: cnc_ptr
    real, dimension(:), pointer :: ptr1d
    integer :: num_fld_out_files, num_nfld_out_files
    integer, dimension(8) :: values 

    if (fu_fails(smpi_is_mpi_version(), 'MPI needed for EnKF', sub_name)) return
    if (fu_fails(smpi_enkf_comm /= int_missing, 'Bad smpi_enkf_comm', sub_name)) return

    nullify(smtr_fieldstorage)

    call set_emis_proc_by_ctrl(model%cloud, model%rules, model%rules%darules%perturbVariable)
    if (error) return

    cloud => model%cloud
    darules = model%rules%darules
    if (fu_fails(darules%method == flag_enks, 'Method not enks', sub_name)) return

    model_integr = model
    simrules_integr = model%rules
    model_integr%rules => simrules_integr ! no longer pointing to the original rules
    model_integr%obs_ptr => obs_ptr_missing

    if (defined(darules%smoother_step) .and. darules%num_smoother_output_steps > 0) then
      ! The number of steps required is determined by the earliest step requested for
      ! output. After the fields need to be stored with the interval of smoother_step.
      smoother_window = darules%smoother_step*maxval(darules%smoother_output_steps)
      call msg('EnKS: smoother window required:' // fu_str(smoother_window))
      ! outputsteps + 1 for the analysis time
      num_smtr_steps = smoother_window / darules%smoother_step + 1
    else if (defined(model%rules%darules%interval_btwn_assim)) then
      smoother_window = interval_missing
      num_smtr_steps = 1
    else
      call set_error('Neither smoother or analysis step is defined', sub_name)
      return
    end if
    step_integr = model%rules%darules%interval_btwn_assim

    allocate(smtr_step_ptrs(num_smtr_steps), smtr_times(num_smtr_steps), stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
    smtr_times = time_missing
    call init_ensemble_multitime(ensemble_multitime, model%cloud, darules, num_smtr_steps)
    if (error) return
    call set_smtr_ptrs()
    if (error) return

    ! Allocate work storage, this version only includes MPI-parallel analysis.
    !
    call msg('ENKF: state and ensemble sizes:', ensemble_multitime%state_size, ensemble_multitime%ens_size)
    if (.not. allocated(ens_full)) then
      call msg('ENKF: allocating work arrays for the MPI mode')
      call get_enkf_chunk_size_offset(ensemble_multitime%state_size, ensemble_multitime%ens_size, &
                                    & chunk_size, chunk_offset)
      call msg('ENKF: my chunk size and offset:', chunk_size, chunk_offset)
      allocate(ens_full(chunk_size, ensemble_multitime%ens_size), stat=stat)
    end if
    if (fu_fails(stat == 0, 'Allocating full ensemble failed', sub_name)) return

    model_integr%rules%periodToCompute = step_integr

    simrules_3d = model%rules
    if (error) return
    now = model%rules%startTime
    !if (now < model%rules%darules%first_analysis_time) then
    !  call set_error('doesn''t work yet', sub_name)
    !  return
    !end if

    end_time = model%rules%startTime + model%rules%periodToCompute

    ! init enkf etc
    allocate(model_integr%perturbations(count_max_perturbations), stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
    call set_perturbations(model_integr%perturbations, &
                         & fu_species_emission(model%cloud), &
                         & fu_species_transport(model%cloud), &
                         & darules, dispersion_grid, dispersion_vertical, now)
    if (error) return

    !call init_random_seed(smpi_global_rank)
    call date_and_time(VALUES=values)
    call init_random_seed(smpi_global_rank*values(8)*values(7))
    simrules_integr%if_finalize = .false.
    
    call msg('EnKS: da_window:' // fu_str(model%rules%darules%assim_window))
    call msg('EnKS: assimilation_interval:' // fu_str(model%rules%darules%interval_btwn_assim))
    call msg('EnKS: smoother_step:' // fu_str(model%rules%darules%smoother_step))
    call msg('EnKS: last analysis time:' // fu_str(darules%last_analysis_time))
    if (num_smtr_steps > 0) then
      smtr_file_units => fu_work_int_array()
      smtr_file_units_nfld => fu_work_int_array()
      smtr_file_units_nfld = int_missing
      call init_smtr_output(ensemble_multitime%states(1), model%out_def, model%rules, &
                          & smtr_file_units, smtr_file_units_nfld, smtr_fieldstorage)
      if (error) return
    else
      nullify(smtr_fieldstorage, smtr_file_units)
    end if
    
    time_analysis_next = now + model%rules%darules%interval_btwn_assim
    time_get_control_next = now
    call start_cycle(now, time_analysis_next)

    !call test_init_control_array(darules, model%cloud)
    !stop
    cnc_ptr => fu_concMM_ptr(model%cloud)

!!$    ptr1d => fu_work_array()
!!$    call random_normal(ptr1d)
!!$    cnc_ptr%arm(1,1,:,:,:) = reshape(ptr1d(1:nz_dispersion*ny_dispersion*nx_dispersion), &
!!$                                   & (/nz_dispersion,nx_dispersion,ny_dispersion/))

    !call run_forecast_dummy(ensemble_multitime%states(ensemble_multitime%num_steps), model_integr)
    do while (now < end_time)
      !if (model%rules%darules%output_level > da_out_none_flag) call update_enkf_outdir(model, simrules_3d%darules, now)
      call msg('ENKS: now = ' // fu_str(now))
      !call msg('sum cnc_ptr', sum(cnc_ptr%arm))
      if (now == time_get_control_next) then
        call get_control(now)
        if (error) return
        time_get_control_next = now + step_integr
      end if
      if (now == time_analysis_next) then
        !call make_analysis()
        call msg('Reporting control before assimilate')
        call report(smtr_step_ptrs(ensemble_multitime%num_steps)%ptr)

        call assimilate(now, ensemble_multitime, fatal)
        
        call msg('Reporting control after assimilate')
        call report(smtr_step_ptrs(ensemble_multitime%num_steps)%ptr)
        
        if (fu_fails(.not. fatal, 'Fatal error in assimilate', sub_name)) return
        if (error) call unset_error(sub_name)
        call finish_cycle(now)
        time_analysis_next = now + model%rules%darules%interval_btwn_assim
        call start_cycle(now, time_analysis_next)
      end if
      
      call msg('Run forecast')
      model_integr%rules%startTime = now
      if (time_analysis_next > darules%last_analysis_time) then
        call msg('Assimilation done, run forecast until end')
        model_integr%rules%periodToCompute = end_time - now
        call run_forecast_from_control(smtr_step_ptrs(ensemble_multitime%num_steps)%ptr, model_integr)
        exit
      else
        call msg('Run forecast until next assimilation')
        !call run_forecast_dummy(smtr_step_ptrs(ensemble_multitime%num_steps)%ptr, model_integr)
        call run_forecast_from_control(smtr_step_ptrs(ensemble_multitime%num_steps)%ptr, model_integr)
        now = now + step_integr
      end if
    end do

    if (num_smtr_steps > 0) then 
      call finish_smtr_output(smtr_file_units(1:num_fld_out_files), smtr_file_units_nfld(1:num_nfld_out_files))
      call free_work_array(smtr_file_units)
      call free_work_array(smtr_file_units_nfld)
      if (associated(smtr_fieldstorage)) deallocate(smtr_fieldstorage)
    end if

  contains

    !======================================================================
 
    subroutine start_cycle(time_start, time_end)
      implicit none
      type(silja_time), intent(in) :: time_start, time_end
      type(DA_rules), pointer :: pDA_rules

      call msg('Start cycle')
      model_integr%rules%darules%da_begin = time_start
      model_integr%rules%darules%assim_window = time_end - time_start
      
      pDA_rules => model_integr%rules%darules
      call set_observations(pDA_rules%station_path, &
                          & pDA_rules%DA_begin, pDA_rules%assim_window, pDA_rules%obs_files, &
                          & pDA_rules%num_obs_files_assim, pDA_rules%num_obs_files_eval, &
                          & fu_species_transport(cloud), &
                          & fu_species_optical(cloud), model_integr%obs_ptr)
      if (error) then
        call msg('Failed to set observations, will clean up')
        call unset_Error('start_cycle')
        call destroy_observations(model_integr%obs_ptr)
      end if

    end subroutine start_cycle

    !======================================================================
    
    subroutine finish_cycle(now)
      implicit none
      type(silja_time), intent(in) :: now
      character(len=fnlen) :: filename_out

      character(len=*), parameter :: subname = 'finish_cycle'
      type(silja_field_id), dimension(:), allocatable :: fid_list
      integer :: ind_smtr_step, ind_fld, ind_time, ind_first_file, ind_last_file, num_nfld_vars, stat
      character(len=worksize_string) :: header
      real, dimension(:), pointer :: work

      call msg('Finish cycle')
      call destroy_observations(model_integr%obs_ptr)
      if (num_fld_out_files > 0) then
         work => fu_work_array(size(smtr_fieldstorage, 1))
      endif

      call start_count('enks_output')
      
      call msg('output files field', smtr_file_units(1:4))
      call msg('output files non-field', smtr_file_units_nfld(1:4))

      allocate(fid_list(size(smtr_fieldstorage, 2)), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', subname)) return

      do ind_smtr_step = 1, darules%num_smoother_output_steps
        ind_time = num_smtr_steps - darules%smoother_output_steps(ind_smtr_step)
        if (fu_fails(ind_time > 0, 'bad ind_time', subname)) return
        if (.not. defined(smtr_times(ind_time))) cycle ! maybe too early for this step
        call msg('Writing smoother step ' // fu_str(smtr_times(ind_time)), ind_smtr_step)
        
        if (num_fld_out_files > 0) then

          call control_to_fieldset(smtr_step_ptrs(ind_time)%ptr, fid_list, &
                                 & smtr_fieldstorage, smtr_times(ind_time))
          if (error) return

          do ind_fld = 1, size(smtr_fieldstorage, 2)
            work(1:size(smtr_fieldstorage, 1)) = max(smtr_fieldstorage(:,ind_fld), 0.0)
            call msg('fieldmean', ind_fld, sum(smtr_fieldstorage(:,ind_fld)) / size(smtr_fieldstorage,1))
            call write_next_field_to_netcdf_file(smtr_file_units(ind_smtr_step), fid_list(ind_fld), work)
            if (error) return
          end do
        end if
        if (num_nfld_out_files > 0) then
          call msg('will call non_fields_to_files')
          call get_non_fieldvars_for_control(smtr_step_ptrs(ind_time)%ptr, num_nfld_vars)
          ind_first_file = (ind_smtr_step-1)*num_nfld_vars + 1
          ind_last_file = ind_first_file + num_nfld_vars
          header = 'BEGIN_STEP ' // trim(fu_str(smtr_times(ind_time)))
          call control_non_fields_to_files(smtr_step_ptrs(ind_time)%ptr, &
                                         & smtr_file_units_nfld(ind_first_file:ind_last_file), &
                                         & darules, header)
          if (error) return
        end if
        
      end do
      call stop_count('enks_output')
      if (num_fld_out_files > 0) then
         call free_work_array(work)
      endif
      if (allocated(fid_list)) deallocate(fid_list)

    end subroutine finish_cycle

    !======================================================================
    
    subroutine get_control(now)
      implicit none
      type(silja_time), intent(in) :: now

      type(da_control), pointer :: ctrlptr_tmp
      integer :: ind_step
      real, dimension(:), pointer :: vp

      call msg('Get control')

!!$      vp => fu_values_ptr(ensemble_multitime%states(1))
!!$      vp = 1
!!$      do ind_step = 1, num_smtr_steps
!!$        call msg('step sum before:', sum(fu_values_ptr(smtr_step_ptrs(ind_step)%ptr)))
!!$      end do

      ! Shift the pointers in time. The new time replaces the oldest.
      ! 
      ctrlptr_tmp => smtr_step_ptrs(1)%ptr ! the oldest
      do ind_step = 1, num_smtr_steps - 1
        !if (.not. defined(smtr_times(ind_step+1))) cycle
        smtr_times(ind_step) = smtr_times(ind_step+1)
        smtr_step_ptrs(ind_step)%ptr => smtr_step_ptrs(ind_step+1)%ptr
      end do
      smtr_step_ptrs(num_smtr_steps)%ptr => ctrlptr_tmp
      smtr_times(num_smtr_steps) = now

      do ind_step = 1, num_smtr_steps
        call msg('smoother timestep:' // fu_str(smtr_times(ind_step)), ind_step)
      end do
      
      ! Run control from cloud on the newest
      call start_count('enkf_control_from_cloud')
      call vector_from_model(smtr_step_ptrs(num_smtr_steps)%ptr, model_integr, handle_negatives=.false.)
      call stop_count('enkf_control_from_cloud')
      call msg('Control from cloud')
      !call report(smtr_step_ptrs(num_smtr_steps)%ptr)
      if (error) return
      
!!$      vp => fu_values_ptr(ensemble_multitime%states(1))
!!$      vp = 1
!!$      do ind_step = 1, num_smtr_steps
!!$        call msg('step mean:', &
!!$               & sum(fu_values_ptr(smtr_step_ptrs(ind_step)%ptr))/fu_size(smtr_step_ptrs(ind_step)%ptr))
!!$      end do
      
    end subroutine get_control

    !======================================================================

    subroutine set_smtr_ptrs()
      ! Initialize step pointers to ensemble. Later they will be cycled on each smoother
      ! timestep.
      implicit none
      integer :: ind_step

      do ind_step = 1, num_smtr_steps
        smtr_step_ptrs(ind_step)%ptr => ensemble_multitime%states(ind_step)
      end do

    end subroutine set_smtr_ptrs

    !======================================================================

    subroutine make_analysis()
      implicit none
      call msg('Analysis')
    end subroutine make_analysis

    !======================================================================
    
    subroutine assimilate(now, ens, fatal_error)
      implicit none
      type(silja_time), intent(in) :: now
      type(t_ensemble), intent(inout) :: ens
      logical, intent(out) :: fatal_error

      real, dimension(:), pointer :: mdl_obs_val_single, obs_data, obs_var, ptr_val
      real, dimension(:,:), allocatable :: obs_loc, mdl_obs_val_all
      integer :: my_size, my_offset
      character(len=*), parameter :: sub_name = 'assimilate'
      real, dimension(:,:), allocatable, save :: model_localisation_full, model_localisation_chunk
      real, dimension(:,:), allocatable :: obs_localisation, mdl_obs_ens
      character(len=12) :: time_str
      type(observationPointers), save :: obs_ptr
      logical :: scatter_ok
      integer, dimension(:), pointer :: offsets, sizes
      integer, dimension(2) :: obs_size

      call msg('EnKS: mpi analysis')
      ! Normal errors can be unset and the integration proceeds with the ensemble
      ! unchanged. However, if scatter fails, the run must crash.
      fatal_error = .false.

      model_integr%rules%darules%da_begin = now

      !if (sum(model_integr%obs_ptr%obs_size) == 0) then
      !  call set_error('ENKF: no observations', sub_name)
      !end if
      fatal_error = .not. sync_errors()
      if (fatal_error .or. error) return
      
      obs_size = model_integr%obs_ptr%obs_size
      mdl_obs_val_single => fu_work_array(sum(obs_size))
      obs_var => fu_work_array(sum(obs_size))

      allocate(mdl_obs_ens(sum(obs_size), ens%ens_size), obs_localisation(2, sum(obs_size)), stat=stat)
      if (fu_fails(stat == 0, 'Allocate master obs arrays failed', sub_name)) continue
      fatal_error = .not. sync_errors()
      if (error .or. fatal_error) return

      call collect_model_data(model_integr%obs_ptr, mdl_obs_val_single) !, obs_size)
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (darules%use_log_obs) then
        mdl_obs_val_single = log(mdl_obs_val_single+0.02)
      end if

      call get_localisation(model_integr%obs_ptr, obs_localisation)
      fatal_error = .not. sync_errors()
      if (error .or. fatal_error) return

      call msg('ENKS: gather ensemble')
      call start_count('ensemble_transpose')
      call transpose_ensemble_in(ens%full_state, ens_full, ens%ens_size, ens%state_size, &
                               & mdl_obs_val_single, mdl_obs_ens, sum(obs_size))
      call stop_count('ensemble_transpose')
      if (error) return

      call get_obs_data(model_integr%obs_ptr, obs_data, obs_size)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (darules%use_log_obs) then
        obs_data = log(obs_data+0.02)
      end if

      call collect_variance(model_integr%obs_ptr, obs_var) !, obs_size)
      fatal_error = .not. sync_errors()
      if (fatal_error .or. error) return

      if (.not. allocated(model_localisation_chunk)) then
        ! do this only once, model localisation and chunking won't change.
        ! one-based offsets
        offsets => fu_work_int_array()
        sizes => fu_work_int_array()
        ! The chunk is from the multistep state.
        call get_chunk(ens%state_size, ens%ens_size, offsets, sizes, if_zero_based=.false.)
        if (error) return
        my_offset = offsets(smpi_ens_rank+1)
        my_size = sizes(smpi_ens_rank+1)
        call free_work_array(offsets)
        call free_work_array(sizes)
        allocate(model_localisation_full(2, ens%state_size), &
               & model_localisation_chunk(2, my_size*ens%num_steps), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed (model loc)', sub_name)) return
        call get_localisation(ens%states(1), model_localisation_full, darules)
        !call test_repl_loc()
        call replicate_localisation(model_localisation_full, model_localisation_chunk, &
                                  & ens%num_steps, fu_size(ens%states(1)), &
                                  & my_offset, my_size)
        deallocate(model_localisation_full)
        if (error) then
          ! This shouldn't fail...
          call set_error('Something wrong with handling localisation', sub_name)
          fatal_error = .true.
          return
        end if
      end if

!      call msg('ave mdl_obs_ens', sum(mdl_obs_ens)/size(mdl_obs_ens))
!      call msg('ave ens_full', sum(ens_full)/size(ens_full))
!      call msg('ave obs_data', sum(obs_data(1:obs_size))/size(obs_data(1:obs_size)))
!      call msg('ave obs_var', sum(obs_var(1:obs_size))/size(obs_var(1:obs_size)))

      !call msg('mean of my chunk before:', sum(ens_full) / size(ens_full))
      if (.not. error) then
        call msg('ENKF: update openmp')
        !
        ! Guess, only assimilated observations are needed for the ensemble update
        !
        call enkf_update_openmp(ens_full, mdl_obs_ens, obs_data(1:obs_size(1)), obs_var(1:obs_size(1)), &
                              & obs_localisation, model_localisation_chunk, &
                              & model%rules%darules%loc_dist_m, &
                              & model%rules%darules%loc_type, &
                              & model%rules%darules%enkf_flavor, &
                              & model%rules%darules%rfactor, &
                              & diagn_out_dir=char_missing)
      end if

!      call msg('ave mdl_obs_ens', sum(mdl_obs_ens)/size(mdl_obs_ens))
!      call msg('ave ens_full', sum(ens_full)/size(ens_full))
!      call msg('ave obs_data', sum(obs_data(1:obs_size))/size(obs_data(1:obs_size)))
!      call msg('ave obs_var', sum(obs_var(1:obs_size))/size(obs_var(1:obs_size)))



      fatal_error = .not. sync_errors()
      if (fatal_error) return
      !call msg('mean of my chunk after:', sum(ens_full) / size(ens_full)) 

      call free_work_array(mdl_obs_val_single)
      call free_work_array(obs_var)

      if (.not. error) then
        call msg('ENKF: scatter ensemble')
        call start_count('ensemble_transpose')
        ptr_val => ens%full_state
        call transpose_ensemble_out(ptr_val, ens_full, ens%ens_size, &
                                  & ens%state_size, scatter_ok)
        call stop_count('ensemble_transpose')
        fatal_error = .not.scatter_ok
!!$        if (simrules_3d%darules%output_level > da_out_none_flag) then
!!$          call forward_3d(model_3dvar)
!!$          call dump_observations(model_3dvar%obs_ptr, fu_species_transport(model%cloud), &
!!$                               & trim(simrules_3d%darules%outputdir) // dir_slash // 'obs_an_' &
!!$                               & // trim(fu_str(smpi_ens_rank)))
!!$          call get_ens_mean_and_stdev(ens, model_3dvar%rules, &
!!$                                    & trim(model_3dvar%rules%darules%outputdir) // dir_slash // 'an')
!!$        end if
      end if


      
    end subroutine assimilate

    !======================================================================

    subroutine replicate_localisation(loc_full_step, loc_chunk, num_steps, size_full_step, my_offset, my_size)
      ! Copy the localisation from the single-step, full domain array to multistep, single-chunk array.
      implicit none
      real, dimension(:,:), intent(in) :: loc_full_step
      real, dimension(:,:), intent(out) :: loc_chunk
      integer, intent(in) :: num_steps, size_full_step, my_offset, my_size

      integer :: ind_chunk, ind_full, ind_step
      
      do ind_chunk = 1, my_size
        ind_full = ind_chunk + my_offset - 1
        ind_step = mod(ind_full-1, size_full_step) + 1
        loc_chunk(:,ind_chunk) = loc_full_step(:,ind_step)
      end do

    end subroutine replicate_localisation

    !======================================================================

    subroutine init_smtr_output(control_template, outdef, simrules_fullrun, file_units, file_units_nfld, &
                              & fieldstorage)
      ! Initialize the output for the smoother steps. The output will be to a single
      ! NetCDF4 file - add other formats or more complex arrangements if needed.
      implicit none
      type(da_control), intent(in) :: control_template ! to obtain variable list, etc
      type(silam_output_definition), intent(in) :: outdef
      type(general_dispersion_rules), intent(in) :: simrules_fullrun ! not just one step btw assimilations
      integer, dimension(:), intent(out) :: file_units, file_units_nfld ! non-field output units
      real, dimension(:,:), pointer :: fieldstorage

      character(len=fnlen) :: main_out_dir, filename_step
      integer :: ind_smtr_step, ind_var, fld_size, num_flds, num_out_vars, ind_time, ind_nfld_var, num_nfld_vars,&
           & ind_1d
      character(len=*), parameter ::  subname = 'init_smtr_output'
      type(silja_time) :: time_start, time_first_out
      type(silam_species), dimension(:), allocatable :: list_sp
      integer, dimension(:), allocatable :: list_q
      logical, dimension(:), allocatable :: if3d
      type(ToutputList) :: output_list
      type(ToutputList), dimension(3) :: output_lists
      character(len=fnlen), dimension(100) :: names
      
      time_start = simrules_fullrun%startTime

      nullify(fieldstorage)

      main_out_dir = fu_FNm(fu_output_template(model%out_def), &
                          & time_start, & ! ini_time
                          & time_start, & ! anal_time
                          & zero_interval, &
                          & model%rules%chCaseNm, &
                          & fu_get_id_str_from_id_nbr(model%source,1)) ! same as in io_server
      
      call control_to_fieldset_get_sizes(control_template, num_flds, fld_size)
      if (error) return
      if (num_flds > 0) then
        if (fu_fails(fld_size > 0, 'Strange num_flds or fld_size', subname)) return
        ! Simplify the allocation, num_flds > num_vars surely.
        allocate(fieldstorage(fld_size, num_flds), list_q(num_flds), list_sp(num_flds), if3d(num_flds), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed', subname)) return
        call get_fieldvars_for_control(control_template, num_out_vars, list_q, list_sp, if3d)
        if (error) return
        allocate(output_list%ptrItem(num_flds))
        call msg('Field varibles for output from control variable')
        do ind_var = 1, num_out_vars
          call msg('Quantity:', list_q(ind_var))
          call msg('Species:' // fu_str(list_sp(ind_var)))
          output_list%ptritem(ind_var) = outputLstItem_missing
          output_list%ptritem(ind_var)%quantity = list_q(ind_var)
          output_list%ptritem(ind_var)%if3d = if3d(ind_var)
          output_list%ptritem(ind_var)%species = list_sp(ind_var)
        end do
        output_list%ptritem(num_out_vars+1:) = outputLstItem_missing
        deallocate(list_q, list_sp, if3d)

        file_units = int_missing

        ! open_netcdf_o needs 3 output lists ("dispersion", "meteo", "mass map"). Missing
        ! lists don't work, they must be valid but zero-sized.
        output_lists = (/outputList_missing, outputList_missing, output_list/)
        allocate(output_lists(1)%ptritem(0), output_lists(2)%ptritem(0), stat=stat)
        if (fu_fails(stat == 0, 'Allocate failed wtf', subname)) return
        do ind_smtr_step = 1, simrules_fullrun%darules%num_smoother_output_steps
          ind_time = simrules_fullrun%darules%smoother_output_steps(ind_smtr_step)
          ! the first time this smoother step is available:
          time_first_out = simrules_fullrun%startTime + simrules_fullrun%darules%smoother_step*ind_time
          filename_step = trim(main_out_dir) // dir_slash // 'smoother_' // trim(fu_str(ind_time)) // '.nc4'
          call msg('EnKS output file:' // filename_step)
          call msg('Main out dir:' // main_out_dir)
          file_units(ind_smtr_step) = open_netcdf_file_o(filename_step, &
               & dispersion_grid, dispersion_vertical, & 
               & time_first_out, output_lists, '', &
               & ifAllInOne=.true., ncver=4, &
               & ifMPIIO=.false., fMissingVal=real_missing)
          if (error) return
        end do
        num_fld_out_files = simrules_fullrun%darules%num_smoother_output_steps
      else
        num_fld_out_files = 0
      end if

      call get_non_fieldvars_for_control(control_template, num_nfld_vars, names)
      if (error) return
      num_nfld_out_files = num_nfld_vars * simrules_fullrun%darules%num_smoother_output_steps
      if (num_nfld_vars > 0) then
        ! One file per smtr step and variable
        do ind_smtr_step = 1, simrules_fullrun%darules%num_smoother_output_steps
          do ind_nfld_var = 1, num_nfld_vars
            ind_1d = ind_nfld_var + (ind_smtr_step-1)*num_nfld_vars
            file_units_nfld(ind_1d) = fu_next_free_unit()
            write(unit=filename_step, fmt='(3A,I0,A)') trim(main_out_dir), dir_slash, 'smoother_', &
                 & ind_smtr_step, trim(names(ind_nfld_var))
            call create_directory_tree(fu_dirname(filename_step))
            open(file_units_nfld(ind_1d), file=filename_step, &
               & form='formatted', action='write', status='replace', iostat=stat)
            if (fu_fails(stat == 0, 'Failed to open:'//trim(filename_step), subname)) return
            call  msg('ENKS text output file:' // trim(filename_step), file_units_nfld(ind_1d))
          end do
        end do
      else
        call msg('ENKS: No non-field output variables')
      end if
      
    end subroutine init_smtr_output

    !======================================================================

    subroutine finish_smtr_output(file_units_fld, file_units_nfld)
      implicit none
      integer, dimension(:), intent(in) :: file_units_fld, file_units_nfld
      integer :: ind_file
      
      do ind_file = 1, size(file_units_fld)
        if (file_units_fld(ind_file) == int_missing) cycle
        call close_netcdf_file(file_units_fld(ind_file))
        if (error) return
      end do
      do ind_file = 1, size(file_units_nfld)
        if (file_units_nfld(ind_file) == int_missing) cycle
        close(file_units_nfld(ind_file))
      end do
      
    end subroutine finish_smtr_output

    !======================================================================

    subroutine test_repl_loc()
      implicit none
      
      real, dimension(1,8) :: loc_full_step
      integer, parameter :: num_steps = 2
      real, target :: loc_chunk1(1,3*num_steps), loc_chunk2(1,3*num_steps), loc_chunk3(1,2*num_steps)
      type(silja_rp_2d), dimension(3) :: chunkptr

      integer, dimension(3) :: sizes = (/6,6,4/), offsets=(/1, 7, 13/)
      integer ::ii, ind_chunk

      chunkptr(1)%pp => loc_chunk1
      chunkptr(2)%pp => loc_chunk2
      chunkptr(3)%pp => loc_chunk3
      
      loc_full_step(1,:) = (/(ii, ii=1, 8)/)
      
      do ind_chunk = 1, 3
        call replicate_localisation(loc_full_step, chunkptr(ind_chunk)%pp, num_steps, 8, &
                                  & offsets(ind_chunk), sizes(ind_chunk))
      end do
      
      print *, 'loc_full_step', loc_full_step
      print *, ''

      do ind_chunk = 1, 3
        print *, 'loc', ind_chunk, chunkptr(ind_chunk)%pp
      end do

      !stop
      
    end subroutine test_repl_loc

  end subroutine run_enks
  
  !**********************************************************************************
  
  subroutine get_enkf_chunk_size_offset(state_size, num_chunks, chunk_size, chunk_offset)
    implicit none
    integer, intent(in) :: state_size, num_chunks
    integer, intent(out) :: chunk_size, chunk_offset
    
    integer, dimension(:), pointer :: offsets, sizes
    
    ! just some convenience
    offsets => fu_work_int_array()
    sizes => fu_work_int_array()
    call get_chunk(state_size, num_chunks, offsets, sizes, .false.)
    if (error) return
    chunk_offset = offsets(smpi_ens_rank+1)
    chunk_size = sizes(smpi_ens_rank+1)
    call free_work_array(offsets)
    call free_work_array(sizes)
    
  end subroutine get_enkf_chunk_size_offset


end module ensemble_driver
