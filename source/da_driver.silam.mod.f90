module da_driver
  ! 
  ! 4DVAR data assimilation, driver routines.
  !
  ! This module implements the "top-level" part of 4D variational assimilation. This
  ! includes the optimization algorithm and background covariance model, and the
  ! functional form of the cost function. The operators M and H are effectively defined in
  ! lower-level modules.
  !
  ! Below is an implementation of 4D-VAR algorithm:
  ! - the quadratic cost function is mimimized using gradient descent
  ! - The usual square-root transformation is applied: the cost function is mimized 
  !   with respect to the transformed variable x' = B^-1/2 (x - x_b) (see the 3D-Var modules
  !   for more explanations. In this code, x' is referred as the 'control space' representation, 
  !   and x = B^1/2 x' + x_b is the 'physical space' representation. The transformation itself
  !   is in the da_interface module together with other routines on the control variables.
  !
  ! JV, Jan 2012
  !
  use da_interface
  use optimiser
  use chemical_setup
  use grids_geo
  use silam_levels
  !use pollution_cloud
  use observation_server  !, only : da_rules

  implicit none

!  public set_da_rules
  public fu_DA_time_window

  public gradient_test
  public gradient_test_emis
  public gradient_test_emis_2
  
  private run_xdvar_seq
  private run_4dvar
  private da_msg
  private set_background
  private to_physical
  private to_control
  private set_cov_mdl
  private evaluate
  private get_gradient
  private reload_control
  private fu_next_iter_file
  !  private fu_next_grad_file
  private fu_next_obs_file

  !real, dimension(:,:,:), allocatable, save, private :: emission_correction

  ! Debug/testing
  !
  logical, private, parameter :: compute_transport = .true., compute_transformation = .false., &
       & compute_emission = .false.
  
  type(da_control), private, save :: control_physical
  integer, parameter :: optimizer_precision = 8

contains

  !************************************************************************************

  subroutine da_msg(what, fval, ival)
    implicit none
    character(len=*), intent(in) :: what
    real, intent(in), optional :: fval
    integer, intent(in), optional :: ival

    if (present(fval) .and. present(ival)) then
      call msg('4DVAR: ' // what, fval, ival)
    else if (present(fval)) then
      call msg('4DVAR: ' // what, fval)
    else if (present(ival)) then
      call msg('4DVAR: ' // what, ival)
    else
      call msg('4DVAR: ' // what)
    end if

  end subroutine da_msg

  !************************************************************************************

  subroutine run_xdvar(model)
    implicit none
    type(model_container), intent(inout) :: model
    
    call set_emis_proc_by_ctrl(model%cloud, model%rules, model%rules%darules%controlVariable)
    if (error) return

    call report(fu_emission_processor_ptr(model%cloud), 'emission processor at the start of run_xdvar')

    select case(model%rules%darules%method)
      case (flag_3dvar)  !, flag_4dvar_seq)
        call msg('Starting 3DVAR')
        call run_xdvar_seq(model)
      case (flag_4dvar)
        call msg('Starting 4DVAR')
        call run_4dvar(model)
      case (flag_4dvar_seq)
        call msg('Starting 4DVAR_SEQ')
        call run_xdvar_seq(model)   ! it will control the sequence
      case (flag_h_matrix)
        call msg('Starting h-matrix')
        call run_h_matrix(model)
      case default
        call set_error('Strange DA method', 'run_xdvar')
    end select

  end subroutine run_xdvar

  
  !*******************************************************************************
  
  subroutine run_xdvar_seq(model)
    implicit none
    type(model_container), intent(inout) :: model

    type(da_control) :: analysis, background
    real :: cost
    !real, dimension(:), pointer :: obs_values, obs_variance
    !real, dimension(2) :: cost_bgr
    !character(len=128) :: log_str
    type(da_rules) :: darules, darules_missing
    type(t_background_covariance), target :: bgr_cov
    type(silam_pollution_cloud), pointer :: cloud
    character(len=fnlen) :: analysis_file_name
    integer :: n_obs_values
    
    type(model_container) :: model_integr, & ! the model for forward integration, no DA
                           & model_3dvar, &  ! the dummy model for xdvar minimization 
                           & model_3dvar_adj ! as above but "adjoint"

    type(observationPointers), target :: obs_ptr_missing
    type(general_dispersion_rules), target :: simrules_integr, & ! rules for running the forecasts
         & simrules_3d ! rules for the assimilation model
!    type(da_rules) :: da_rules_step
    type(silja_time) :: now, end_time, time_forecast_end
    character(len=*), parameter :: sub_name = 'run_xdvar_seq'
    type(t_tla_trajectory), target :: tla_traj_missing
    real, dimension(1) :: dummyarr
    logical :: ifProgressFile
    integer :: assimilations_total,  assimilations_failed
    
    cloud => model%cloud
    darules = model%rules%darules
    if (fu_fails(any(darules%method == (/flag_3dvar,flag_4dvar_seq/)), &
                          & 'Bad da method', sub_name)) return

    assimilations_total = 0
    assimilations_failed = 0

    !call test_var_da(model)
    !stop
    model_integr = model
    ! Get rid of all da definitions in model_integr, it is for runs between assimilations
    simrules_integr = model%rules
    model_integr%rules => simrules_integr ! no longer pointing to the original rules
    !simrules_integr%darules = darules_missing
    simrules_integr%darules%allow_negative = .false.
    model_integr%obs_ptr => obs_ptr_missing
    model_integr%tla_traj => tla_traj_missing
    
    ifProgressFile = model%rules%ifWriteProgressFile
    model_integr%rules%ifWriteProgressFile = .false.

    model_3dvar = model
    simrules_3d = model%rules
    if (simrules_3d%darules%method == flag_4dvar_seq) then
      simrules_3d%darules%method = flag_4dvar  !! Not "sequential" so 4dvar does not get confused"
      simrules_3d%periodtocompute = simrules_3d%darules%assim_window
      simrules_3d%darules%run_forecast = .false. !! No need to run "final forecast"
    endif
    model_3dvar%rules => simrules_3d

    !!Here  used to be a code that reads input to control overriding initial conditions, emission processor
    !! of the cloud. If something has to be set in the control it should be done via cloud.
    ! 

    !
    ! Spin-up. 
    ! Run_forecast applies control_physical and integrate until first assimilation if needed
    call da_msg('Spin-up, from:'+fu_str(model%rules%startTime) + ', to:' + fu_str(model%rules%darules%da_begin))
    model_integr%rules%startTime = model%rules%startTime
    model_integr%rules%periodToCompute = model%rules%darules%da_begin - model%rules%startTime
    call model_forward(model_integr, model_integr%obs_ptr, model_integr%rules, make_output=.true.)

    call da_msg('Spin-up complete, Now: '//trim(fu_str(model%rules%darules%da_begin)))

    if(fu_fails(.not.error, 'assimilation stopped due to errors in spin-up','run_xdvar_seq'))return

    ! Now, we are at the DA_begin time
    if(ifProgressFile)then
      call write_progress_file(model_integr%rules%chProgressFNmTemplate, &
                             & fu_compute_progress_time(model%rules%startTime, &
                                            & model%rules%startTime + &
                                                  & model%rules%periodToCompute, &
                                            & model%rules%darules%da_begin))
    endif
    !
    ! Set the two controls, for the assimilation and for between-assimilation runs.
    ! In MPI, only master handles vectors, others get dummy vectors with grid_missing
    !

    call msg ("run_xdvar_seq memusage 1 (kB)", fu_system_mem_usage())
    ! Create control
    call init_control_from_cloud(background,      cloud, darules, physical_space, 0.0) !!! Background/fallback
    call init_control_from_cloud(control_physical, cloud, darules, physical_space, 0.0) !! workspace for DA
    if(error)return

!    ! Set control's values, either from files or defaults
!    call msg("Setting background")
!    call set_background(control_physical,  cloud, darules)
!    call msg("Setting control_physical")
!    call copy_control(background, control_physical, if_allocate=.false.)
!    if(error)return
    !
    !----------------------------------------------------------------------
    !
    ! Main time cycle over the computational period with a sequence of assimilations
    !
    now = model%rules%darules%da_begin
    end_time = model%rules%startTime + model%rules%periodToCompute
    simrules_integr%if_finalize = .false.
    do while (now < end_time)
      call msg ("run_xdvar_seq memusage loop_beg (kB)", fu_system_mem_usage())
      call vector_from_model(control_physical, model, handle_negatives=.true.)  !!! Save the state for the case DA fails 
      if (.not. darules%no_background_term) then
        call copy_control(control_physical, background, if_allocate=.false.)
      end if
      call msg ("run_xdvar_seq memusage after setting controls(kB)", fu_system_mem_usage())

      call update_da_outdir(model_integr, simrules_3d%darules, now)
      if (error) return
      !
      ! The assimilation itself: minimization of the deviation from obs
      !      
      if (darules%method == flag_3dvar) then
        call msg ("run_xdvar_seq memusage before assimilate (kB)", fu_system_mem_usage())
        call start_count('Assimilate')
        call assimilate(background)
        call stop_count('Assimilate')
        call msg ("run_xdvar_seq memusage after assimilate (kB)", fu_system_mem_usage())
      else if (darules%method == flag_4dvar_seq) then
        call assimilate_4d(background)
      else
        call set_error('Bad darules%method', sub_name)
        return
      end if
      !
      ! Assimilation finished. Next steps
      !
     
      assimilations_total = assimilations_total + 1
      if(error) then
        assimilations_failed = assimilations_failed + 1
        call msg('assimilate returned with error, fallback to background'// &
               & " failed "//trim(fu_str(assimilations_failed))// &
               & " out of " // trim(fu_str(assimilations_total)))
        call copy_control(background, control_physical, if_allocate=.false.)
        if (defined(analysis)) call destroy(analysis)
        if (defined(bgr_cov)) call kill_bgr_cov(bgr_cov)
        call unset_error('run_xdvar_seq')
      else
        !
        ! 4Dvar returns analysis in control space, while 3Dvar seems to set it in physical
        !
        if(.not. fu_in_physical_space(analysis)) &
                      & call to_physical(analysis, background, control_physical)
        call msg ("run_xdvar_seq memusage after to_physical (kB)", fu_system_mem_usage())
        call destroy(analysis)
        call msg ("run_xdvar_seq memusage after destroy(analysis) (kB)", fu_system_mem_usage())
        call kill_bgr_cov(bgr_cov)
        call msg ("run_xdvar_seq memusage after kill_bgr_cov (kB)", fu_system_mem_usage())
      end if
      !
      ! Run from the assimilation time till next assimilation period - or till the end
      ! of the whole run if all assimilations are done
      !
      

      !! Period for coming forecast
      time_forecast_end =  now + model%rules%darules%interval_btwn_assim
      if (time_forecast_end  > darules%last_analysis_time) then
        call da_msg('Last analysis done')
        time_forecast_end = end_time
        simrules_integr%if_finalize = .true.
        call da_msg('Will do final forecast')
      end if
      ! Times for the coming forecast
      model_integr%rules%startTime = now
      model_integr%rules%periodToCompute = time_forecast_end - now
      call da_msg('time_forecast_end:' // trim(fu_str(time_forecast_end)))
      call da_msg('last analysis time:' // trim(fu_str(darules%last_analysis_time)))

      !!Dirty hack: feed forward model with the trajectory, so non-reentrant sources still can emit
      if (if_store_forward(model_3dvar%tla_traj)) model_integr%tla_traj => model_3dvar%tla_traj
      !
      ! The main-model run, eventually calls run_dispersion
      ! For the case progress file is needed, calculate the fraction of the period done
      !
      call run_forecast_from_control(control_physical, model_integr)
      call msg ("run_xdvar_seq memusage after run_forecast (kB)", fu_system_mem_usage())
      if (error) return
      now = time_forecast_end

      if(ifProgressFile)then
        call write_progress_file(model_integr%rules%chProgressFNmTemplate, &
                               & fu_compute_progress_time(model%rules%startTime, &
                                              & model%rules%startTime + &
                                                    & model%rules%periodToCompute, &
                                              & now))
      endif

    end do  ! main cycle over times

    if (assimilations_failed == assimilations_total) then
      call set_error("All assimilations failed!", 'run_xdvar_seq')
      return
    endif
    !
    ! All done, close the progress file
    !
    if(ifProgressFile) call write_progress_file(model_integr%rules%chProgressFNmTemplate, 100)

  contains
    
    subroutine assimilate(background)
      implicit none
      type(da_control), intent(in) :: background
      !type(da_control), intent(inout) :: control

      real, dimension(:), pointer :: p_values
      logical :: ok
      integer :: iPurpose
      character(len=fnlen) :: analysis_file_name, background_file_name
      type(Tmass_map), save, target :: mass_map_adjoint
      type(da_rules), pointer :: pDA_rules

      model_3dvar%rules%darules%da_begin = now
      model_3dvar_adj = model_3dvar

      if (.not. defined(mass_map_adjoint)) then
        call msg('Creating adjoint mass map')
        call set_aux_disp_massmap(mass_map_adjoint, 1,  fu_species_transport(cloud), &
            &fu_nbr_of_species_transport(cloud), concentration_flag)
        if (error) return
      end if
      model_3dvar_adj%mass_map_adj => mass_map_adjoint

      call msg ("run_xdvar_seq memusage before assimilate (kB)", fu_system_mem_usage())

      pDA_rules => model_3dvar%rules%darules
      !
      ! Read observations
      !
      call set_observations(pDA_rules%station_path, &
                          & pDA_rules%DA_begin, pDA_rules%assim_window, pDA_rules%obs_files, &
                          & pDA_rules%num_obs_files_assim, pDA_rules%num_obs_files_eval, &
                          & fu_species_transport(cloud), &
                          & fu_species_optical(cloud), model%obs_ptr)
      if (error) return
      call msg ("run_xdvar_seq memusage after set_observations (kB)", fu_system_mem_usage())

      if (smpi_adv_rank == 0) then
        call set_cov_mdl(bgr_cov, background, fu_species_emission(cloud), fu_species_transport(cloud), &
                       & model_3dvar%rules%darules, wholeMPIdispersion_grid, &
                       & dispersion_vertical, now, model_3dvar%rules%darules%controlVariable)
        if (error) return
        call msg ("run_xdvar_seq memusage after set_cov_mdl (kB)", fu_system_mem_usage())
      else
        bgr_cov = background_covariance_missing
      endif 
      call init_control_from_cloud(analysis, cloud, darules, control_space, 0.0, bgr_cov)
      if (error) return
        call msg ("run_xdvar_seq memusage after init_control (kB)", fu_system_mem_usage())

      p_values => fu_values_ptr(analysis)
      
      !call gradient_test(model_3dvar)
      !stop
      if (smpi_adv_rank == 0) then
        if (simrules_3d%darules%output_level >= da_out_first_last_flag) then
          background_file_name = simrules_3d%darules%outputdir + dir_slash + 'background_da_' + &
                               & fu_time_fname_string_utc(simrules_3d%darules%da_begin) + '.grads'
          call control_to_file(background, background, model_3dvar%rules%darules, background_file_name)
          call msg ("run_xdvar_seq memusage after control_to_file (kB)", fu_system_mem_usage())
        end if
      endif

      select case(model%rules%darules%numerics%searchMethod)
      case (steepest_descent_flag)
        call steepest(model_3dvar, model_3dvar_adj, analysis, background, bgr_cov, report_iterates=.false.)
      case (m1qn3_flag)
        call quasi_newton(model_3dvar, model_3dvar_adj, analysis, background, bgr_cov, .false.)
      case (l_bfgs_b_flag)
        call l_bfgs_b(model_3dvar, model_3dvar_adj, analysis, background, bgr_cov, .false.)
      case default
        call set_error('Strange search method', 'run_xdvar_seq__assimilate')
        return
      end select
      call msg ("run_xdvar_seq memusage after optimize (kB)", fu_system_mem_usage())
    
      if (error) then
        call msg('Error on analysis step, reverting to background')
        call unset_error('assimilate')
        p_values = 0.0
      elseif (simrules_3d%darules%output_level > da_out_none_flag .and. smpi_adv_rank == 0) then
        !
        ! Store the new analysis to file
        !
        analysis_file_name = simrules_3d%darules%outputdir + dir_slash + 'analysis_da_' + &
                           & fu_time_fname_string_utc(simrules_3d%darules%da_begin) + '.grads' 
        !        call ooops("writing analysis.grads, valid time: "//fu_str(model%rules%darules%da_begin)) 
        call control_to_file(analysis, background, model%rules%darules, analysis_file_name)
        call msg ("run_xdvar_seq memusage after control_to_file (kB)", fu_system_mem_usage())
        !
        ! Dump observations
        !
        do iPurpose = 1, 2
          call dump_observations(model%obs_ptr, fu_species_transport(cloud), &
                & simrules_3d%darules%outputdir + dir_slash + 'obs_final' + chObsPurpose(iPurpose) + &
                & 'da_' + fu_time_fname_string_utc(model%rules%darules%da_begin), &
                & iPurpose)
        end do
        call msg ("run_xdvar_seq memusage after dump_observations (kB)", fu_system_mem_usage())
      endif

      call destroy_observations(model%obs_ptr)
      call msg ("run_xdvar_seq memusage after destroy_observations (kB)", fu_system_mem_usage())
      
    end subroutine assimilate

    !=================================================================

    subroutine assimilate_4d(background)
      implicit none
      type(da_control), intent(in) :: background

      !call msg('BEFORE_ASSIMILATE_4D')
      !call report_real_work_arrays()
      
      ! should force recomputation of low mass thresholds. This implies that the arrays become
      ! distinct from those used for model_integr.
      !call init_low_mass_threshold(simrules_3d%chemicalrules, &
      !                          & fu_nbr_of_species_transport(model_3dvar%cloud), &
      !                          & if_def_lagr=.false.)
      !if (error) return
      if (fu_true(simrules_3d%darules%have_analysis_species)) then
        call set_default_cnc(fu_concMM_ptr(model_integr%cloud))
      end if
      model_3dvar%rules%darules%da_begin = now
      model_3dvar%rules%startTime = now

      call run_4dvar(model_3dvar, background, analysis, bgr_cov)

      call set_bgr_cov_ptr(analysis, bgr_cov)

    end subroutine assimilate_4d

  end subroutine run_xdvar_seq

  
  !***********************************************************************************
  
  subroutine run_4dvar(model, background_in, analysis_out, bgr_cov_out)
    implicit none
    ! imported: 
    type(model_container), intent(inout), target :: model
    ! optional arguments for the 4d-sequential use
    type(da_control), intent(in), optional :: background_in
    type(da_control), intent(out), optional :: analysis_out
    type(t_background_covariance), intent(out), optional :: bgr_cov_out

    type(da_control) :: control, gradient, background
    real :: cost
    real, dimension(:), pointer :: obs_values => null(), obs_variance =>null(), background_values =>null()
    real, dimension(2) :: cost_bgr
!    character(len=128) :: log_str
    type(da_rules), pointer :: darules
    type(t_background_covariance) :: bgr_cov
    type(silam_pollution_cloud), pointer :: cloud => null()
    character(len=fnlen) :: analysis_file_name, background_file_name
    integer :: iPurpose
    type(model_container) :: model_adj
    type(silam_pollution_cloud), pointer :: cloud_adj =>null()
    type(Tsilam_namelist) :: nl_dummy
    integer, save :: num_timesteps_prev = -1
    type(DA_rules), pointer :: pDA_rules
    type(TMass_map), pointer :: mapConc

    cloud => model%cloud
    darules => model%rules%darules
    model%rules%if_finalize = .false.
    !model%rules%if_collect_linearization = .true.
    !model%rules%if_collect_linearization = .false.
    
    call msg('Start of 4DVAR, reporting rules')
    call report(darules)

    call da_msg('Read observations')
    pDA_rules => model%rules%darules
    call set_observations(pDA_rules%station_path, &
                        & pDA_rules%DA_begin, pDA_rules%assim_window, pDA_rules%obs_files, &
                        & pDA_rules%num_obs_files_assim, pDA_rules%num_obs_files_eval, &
                        & fu_species_transport(model%cloud), &
                        & fu_species_optical(model%cloud), model%obs_ptr)
    if (error) return

    call da_msg('===========================================================')
    call da_msg('4dvar starting')
    call da_msg('Analysis time: ' // trim(fu_str(darules%da_begin)))
    !
    ! Get observation data and their variances
    !
    obs_variance => fu_work_array(sum(model%obs_ptr%obs_size))
    obs_values => fu_work_array(sum(model%obs_ptr%obs_size))
    call collect_obs_data(model%obs_ptr, obs_values) !, n_obs_values)
    if (error) return

    call collect_variance(model%obs_ptr, obs_variance) !, n_obs_values)
    if (error) return
    !
    ! Inititalise controls in physical space
    ! Background can be given from above (sequential 4D-var) or made here
    ! Control_physical is a storage place, will be used later for reloading the
    ! best-iteration control
    !
    call init_control_from_cloud(control_physical, cloud,  darules, physical_space)
    call init_control_from_cloud(background, cloud, darules, physical_space)  ! must be initialized to be copied to
    if (error) return
    if (present(background_in)) then
      ! Copy input and do not touch it further on
      call msg('Reporting control: background_in')
      call copy_control(background_in, background, if_allocate=.false.)
    else
      ! Gets initial and emission corrections from files or sets their defaults from rules
      call msg('Reporting control: from set_background')
      call set_background(background, cloud, darules)
    end if
    if (error) return
    call report(background)
    call msg('Reporting control: control_physical')
    call report(control_physical)
    !
    ! Allow for no-background assimilation. Note that the default background is 0 for model state
    ! and 1 for emission assimilation
    !
    if (darules%no_background_term) then
      !
      ! no background term in cost function -> set background to zero after setting the
      ! initial condition. See observation_server.
      !
      if (fu_have_initial(fu_mode(background))) then
        background_values => fu_initial_ptr(background)
        background_values(1:fu_n_initial(background)) = 0.0
      endif
      if (fu_have_emission_xy(fu_mode(background))) then
        background_values => fu_emission_ptr(background)
        background_values(1:fu_n_emission(background)) = 1.0
      endif
    end if   ! no background term

    call update_da_outdir(model, darules, darules%da_begin)
    if (error) return

    if (darules%output_level > da_out_none_flag) then
      background_file_name = model%rules%darules%outputdir + dir_slash + 'background_da_' + &
                           & fu_time_fname_string_utc(model%rules%darules%da_begin) + '.grads' 
      call control_to_file(background, background, model%rules%darules, background_file_name)
      call da_msg('Background written to ' // trim(background_file_name))
    end if
    call msg('Reporting background')
    call report(background)

    call set_cov_mdl(bgr_cov, background, fu_species_emission(cloud), &
                   & fu_species_transport(cloud), darules, &
                   & wholeMPIdispersion_grid, dispersion_vertical, &
                   & darules%da_begin, darules%controlVariable)
    if (error) return
    call init_control_from_cloud(control, cloud, darules, control_space, 0.0, bgr_cov)
    call msg('Repotring control: control')
    call report(control)

    if(darules%restart_iteration) &
              & call reload_control(cloud, darules, control, background, int_missing, control_space)
    if(error)return
    !
    ! Init TLA at the first entry (num_timesteps_prev = -1 at initialization)
    ! or if the number of assimilation timesteps has changed - seemingly, never the case now
    !
    if (num_timesteps_prev /= fu_nbr_of_timesteps(model%rules)) then
      !
      ! TLA can be requested for chemistry and for emission
      !
      call get_tla_chemistry(model%rules%chemicalRules, model%tla_traj) !! Depositions hndled in there

      mapConc => fu_emisMM_ptr(model%cloud)  !! Emission massmap, just reuse pointer here
      call add_tla_emission(model%source, model%tla_traj, size(mapConc%arM, 1))
      if (error .or. .not. defined(model%tla_traj)) return
      !
      ! get the storage space
      !
      mapConc => fu_concMM_ptr(model%cloud)
      call alloc_tla_trajs(model%tla_traj, mapConc%nx, mapConc%ny, mapConc%n3d, &
                         & fu_nbr_of_timesteps(model%rules), 10000000000_8 / smpi_adv_tasks, & 
                         & model%rules%darules%outputdir ) 
                       !!! Keep 10G memory storage for all memebers together
      if (error) return
      num_timesteps_prev = fu_nbr_of_timesteps(model%rules)
    end if

    model_adj = model
    !allocate(cloud_adj)
    !call init_pollution_cloud(cloud_adj, fu_advection_method(model%cloud), model%rules%timestep, &
    !                        & nl_dummy, model%source, model%rules%chemicalRules)
    !if (error) return
    !call source_to_initial_cloud(model%source, cloud_adj, model%rules%chemicalRules, &
    !                           & model%rules%startTime, model%rules%periodToCompute, &
    !                           & model%meteo_market_ptr, dispersion_gridptr, & ! lagrangian output grid
    !                           & 10)!, 2) ! accuracy, number of sources
    !if (error) return
    !call set_meteo2disp_interp(cloud_adj, model%wdr)
    !if (error) return
    !model_adj%cloud => cloud_adj
    !call msg('Adjoint cloud:')
    !call report(cloud_adj)
    !call msg('Forward cloud:')
    !call report(model%cloud)


    !call gradient_test(model)
    !call gradient_test_emis(model)
    !call gradient_test_emis_2(model)
    !stop

    call start_count('Minimization cycle')
    
#ifdef DEBUG    
    call report(fu_emission_processor_ptr(cloud), 'Minimization cycle starts')
    call msg('Reporting control: control')
    call report(control)
    call msg('Reporting control: background')
    call report(background)
#endif

    select case(model%rules%darules%numerics%searchMethod)
      case (steepest_descent_flag)
        call msg('4DVar steepest descent')
        call steepest(model, model_adj, control, background, bgr_cov, report_iterates=.true.)
      case (m1qn3_flag)
        call msg('4DVAR m1qn3')
        call quasi_newton(model, model_adj, control, background, bgr_cov, .true.)
      case (l_bfgs_b_flag)
        call msg('4DVAR l_bfgs_b')
        call l_bfgs_b(model, model_adj, control, background, bgr_cov, .true.)
      case default
        call set_error('Strange search method', 'run_4dvar')
        return
    end select
    if (error) return
    call da_msg('===========================================================')
    call stop_count('Minimization cycle')

    call to_physical(control, background, control_physical)
    !
    ! Store the analysis
    !
    if (model%rules%darules%output_level > da_out_none_flag) then
      analysis_file_name = model%rules%darules%outputdir + dir_slash + 'analysis_da_' + &
                         & fu_time_fname_string_utc(model%rules%darules%da_begin) + '.grads' 
      call control_to_file(control_physical, background, model%rules%darules, analysis_file_name)
      if (error) return
      call da_msg('Analysis written to ' // trim(analysis_file_name))
    end if
    !
    ! Run the forecast (!!FIXME ) The forecast will be run anyway, should we do it here?
    !
    !
    if (darules%run_forecast) then
      call da_msg('run_4dvar Minimization finished, will run forecast')
      call run_forecast_from_control(control_physical, model)
    else
      call da_msg('run_4dvar Minimization finished no run_forecast')
    endif
    if (present(analysis_out)) call copy_control(control, analysis_out, if_allocate=.true.)
    if (present(bgr_cov_out)) bgr_cov_out = bgr_cov
    if (error) return

    if (model%rules%darules%output_level > da_out_none_flag) then    
      do iPurpose = 1, 2
        call dump_observations(model%obs_ptr, fu_species_transport(cloud), &
                & model%rules%darules%outputdir + dir_slash + 'obs_final'  + chObsPurpose(iPurpose) + &
                & 'da_' + fu_time_fname_string_utc(model%rules%darules%da_begin), &
                & iPurpose)
      end do
      call da_msg('Observations from forecast written to ' &
                & // trim(model%rules%darules%outputdir + dir_slash + 'obs_final' + '_da_' + &
                & fu_time_fname_string_utc(model%rules%darules%da_begin)))
    end if

    call destroy_observations(model%obs_ptr)
    call destroy(control)
    call destroy(background)
    call free_work_array(obs_variance)
    call free_work_array(obs_values)

  end subroutine run_4dvar

  
  !************************************************************************************

  subroutine update_da_outdir(model, darules, now)
    !
    ! Must be called by all MPI members
    !
    implicit none
    type(model_container), intent(in) :: model 
    type(da_rules), intent(inout) :: darules
    type(silja_time), intent(in) :: now

    character(len=fnlen) :: template_string, main_out_dir
    type(grads_template) :: gr_template
    type(silja_interval) :: forecast_length

    ! The 4dvar "iteration" model has no outDef, here it is needed though!
    if (fu_fails(associated(model%out_def), 'model%out_def not associated', 'update_da_outdir')) return

    if (darules%output_level == da_out_none_flag) then   !Nothing will be output. No directory needed
      darules%outputdir=char_missing
      return
    endif 

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
    call da_msg('Updating output directory: ' // trim(template_string) // '-', ival=len_trim(template_string))
    call decode_template_string(template_string, gr_template)
    darules%outputdir = fu_FNm(gr_template, &
                           & model%rules%startTime, &
                           & now, &
                           & zero_interval, &
                           & model%rules%chCaseNm, &
                           & fu_get_id_str_from_id_nbr(model%source,1)) ! same as in io_server

    !Only one creates directory, otherwie the call might get stuck (did at Mahti)
    if (smpi_adv_rank == 0) call create_directory_tree(trim(darules%outputdir))
    call smpi_advection_barrier()

  end subroutine update_da_outdir

  
  !***********************************************************************************
  
  subroutine reload_control(cloud, darules, control, background, iter_to_recover, space_switch)
    !
    ! Resets the 4D-VAR environment to the last-processed interation, or to the given one
    !
    implicit none

    ! Imported parameters
    type(silam_pollution_cloud), pointer :: cloud
    type(da_rules), intent(in) :: darules
    type(da_control), target, intent(inout) :: control
    type(da_control), intent(in) :: background
    integer, intent(in) :: iter_to_recover, space_switch
      
    integer :: iter, last_iter
    integer, parameter :: max_iter_to_chk = 500
    logical :: exist
    type(da_control), target :: control_phys
    type(da_control), pointer :: control_ptr
    character(len=fnlen) :: iterfilename, output_dir

    output_dir = darules%outputdir
    call msg('Trying to reload control from directory:' + output_dir)
    if(iter_to_recover == int_missing)then
      !
      ! get the last successful one
      !
      do iter = 1, max_iter_to_chk
        iterfilename = fu_next_iter_file(output_dir, 'control_da_', darules%da_begin, iter)
        inquire(file=iterfilename, exist=exist)
        if (.not. exist) exit
      end do
      if (iter == 1 .or. iter > max_iter_to_chk) then
        call set_error('reload_control requested but no files found', 'reload_control')
        return
      endif
      last_iter = iter-1
    else
      last_iter = iter_to_recover
    endif
    !
    ! Exists? Add GRADS as format and proceed
    !
    iterfilename = fu_next_iter_file(output_dir, 'control_da_', darules%da_begin, last_iter) + '.super_ctl'
    call msg('Checking:' + iterfilename)
    inquire(file=iterfilename, exist=exist)
    if (exist) then
      call msg('Trying to reload_control, iteration from file: ' // trim(iterfilename))
      iterfilename = 'GRADS ' // trim(iterfilename)
    else
      call set_error('reload_control is requested but file not found: '//trim(iterfilename), 'reload_control')
      return
    end if
    !
    ! use set_background to read the control. If the control is needed in the control space
    ! use thermodynamic_tools intermediate to read the stores control - that one is in physical space
    !
    if(space_switch == control_space)then
      control_ptr => control_phys
    elseif(space_switch == physical_space)then
      control_ptr => control
    endif
    call init_control_from_cloud(control_ptr, cloud, darules, physical_space)
    if (error) return
    call background_from_files(control_ptr, cloud, &
                             & iterfilename, .true., real_missing, &    ! initial-condition file
                             & iterfilename, .true., real_missing, &    ! emission file
                             & darules%da_begin, darules%ifRandomise)
    if (any(fu_values_ptr(control_ptr) .eps. real_missing)) then
      call set_error('Failed setting the control from reload', 'reload_control')
      return
    end if
    if(space_switch == control_space)then
      ! try to recover it using recover_control
      call recover_control(control_phys, control, background)
      call destroy(control_phys)
    elseif(space_switch == physical_space)then
      ! just return
    else
      call set_error('Unknown space_switch:' + fu_str(space_switch),'reload_control')
    endif
      
  end subroutine reload_control

  !************************************************************************************

  subroutine run_h_matrix(model)
    use netcdf
    implicit none
    type(model_container), intent(inout), target :: model

    type(da_control) :: control, background
    type(silam_pollution_cloud), pointer :: cloud
    type(da_rules), pointer :: darules

    integer :: contr_space_dim, ind_contr, n_obs_values
    real, dimension(:), pointer :: contr_val, mdl_obs_values, obs_values, obs_variance
    real :: cost, unit
    real, dimension(2) :: cost_obs
    real, dimension(5) :: cost_bgr
    integer, dimension(2) :: ncid, obs_space_dim
    integer :: stat, nc_var_id, iPurpose
    type(t_background_covariance) :: bgr_cov
    character(len=*), parameter :: sub_name = 'h_matrix'
    character(len=fnlen) :: filepath
    type(DA_rules), pointer :: pDA_rules
    
    cloud => model%cloud
    darules => model%rules%darules
    
    call da_msg('Read observations')
    pDA_rules => model%rules%darules
    call set_observations(pDA_rules%station_path, &
                        & pDA_rules%DA_begin, pDA_rules%assim_window, pDA_rules%obs_files, &
                        & pDA_rules%num_obs_files_assim, pDA_rules%num_obs_files_eval, &
                        & fu_species_transport(model%cloud), &
                        & fu_species_optical(model%cloud), model%obs_ptr)
    if (error) return

    call init_control_from_cloud(background, cloud, darules, physical_space, 0.0)
    if (error) return
    call set_cov_mdl(bgr_cov, background, fu_species_emission(cloud), &
                   & fu_species_transport(cloud), darules, &
                   & wholeMPIdispersion_grid, dispersion_vertical, &
                   & darules%da_begin, darules%controlVariable)
    if (error) return
    call init_control_from_cloud(control, cloud, darules, control_space, bgr_cov=bgr_cov)
    if (error) return

    contr_space_dim = fu_size(control)
    contr_val => fu_values_ptr(control)
    obs_space_dim = model%obs_ptr%obs_size   ! (assim,eval)

    call da_msg('Creating H matrix')
    call da_msg('Control space dimension:', ival=contr_space_dim)
    call da_msg('Observation space assimilation dimension:', ival=obs_space_dim(1))
    call da_msg('Observation space evaluation dimension:', ival=obs_space_dim(2))

    mdl_obs_values => fu_work_array(sum(model%obs_ptr%obs_size))
    !obs_values => fu_work_array()
    !obs_variance => fu_work_array()
    
    obs_values => fu_obs_data(model%obs_ptr)
    obs_variance => fu_obs_variance(model%obs_ptr)
        
    call update_da_outdir(model, darules, model%rules%startTime)
    if (error) return
    filepath = trim(darules%outputdir) // dir_slash // 'hmatrix.nc'
    do iPurpose = 1, 2
      call init_nc(filepath + chObsPurpose(iPurpose), obs_values, obs_variance, &
                 & obs_space_dim(iPurpose), contr_space_dim, ncid(iPurpose), nc_var_id)
      if (error) return
    end do

    if (darules%time_height_mode == processor_tm_hgt_force_flag) then
      unit = darules%h_matrix_unit / tm_hgt_force_scaling
    else
      unit = 1.0
    end if

    do ind_contr = 1, contr_space_dim
      contr_val = 0.0
      contr_val(ind_contr) = unit
      call da_msg('Processing H matrix column', ival=ind_contr)
      call evaluate(control, background, bgr_cov, model, cost_obs, cost_bgr, cost, mdl_obs_values)
      if (error) return
      call column_to_nc(ncid, nc_var_id, ind_contr, mdl_obs_values, obs_space_dim)
      if (error) return
    end do

    ! Close the nc file
    do iPurpose = 1, 2
      if (fu_fails(nf90_close(ncid(iPurpose)) == NF90_NOERR, 'Failed to close netcdf',sub_name))return
    end do

  contains
    
    subroutine init_nc(filename, obs_val, obs_var, n_obs, n_contr, ncid, mdl_var_id)
      implicit none
      character(len=*), intent(in) :: filename
      real, dimension(:), intent(in) :: obs_val, obs_var
      integer, intent(in) :: n_obs, n_contr
      
      integer, intent(out) :: ncid, mdl_var_id
      
      integer :: obs_val_id, obs_var_id, stat, obs_dim_id, contr_dim_id

      stat = nf90_create(filename, 0, ncid)
!      stat = nf90_create(filename, NF90_CLASSIC_MODEL, ncid)
      !stat = nf90_open(filename, NF90_WRITE, ncid)
      call msg('stat', stat)
      if (fu_fails(stat == NF90_NOERR, 'Failed to open: ' // trim(filename), sub_name)) return

      stat = nf90_def_dim(ncid, 'observation', n_obs, obs_dim_id)
      if (fu_fails(stat == NF90_NOERR, 'Failed to create obs dimension', sub_name)) return

      stat = nf90_def_dim(ncid, 'control', n_contr, contr_dim_id)
      if (fu_fails(stat == NF90_NOERR, 'Failed to create contr dimension', sub_name)) return

      stat = nf90_def_var(ncid, 'mdl_value', nf90_float, (/obs_dim_id, contr_dim_id/), mdl_var_id)
      if (fu_fails(stat == NF90_NOERR, 'Failed to create mdl val', sub_name)) return
      
      stat = nf90_def_var(ncid, 'obs_value', nf90_float, (/obs_dim_id/), obs_val_id)
      if (fu_fails(stat == NF90_NOERR, 'Failed to create obs val', sub_name)) return
      stat = nf90_def_var(ncid, 'obs_variance', nf90_float, (/obs_dim_id/), obs_var_id)
      if (fu_fails(stat == NF90_NOERR, 'Failed to create obs var', sub_name)) return

      stat = nf90_enddef(ncid)
      if (fu_fails(stat == NF90_NOERR, 'Failed to enddef netcdf', sub_name)) return

      stat = nf90_put_var(ncid, obs_val_id, obs_val(1:n_obs), start=(/1/), count=(/n_obs/))
      if (fu_fails(stat == NF90_NOERR, 'Failed to write obs val', sub_name)) return
      stat = nf90_put_var(ncid, obs_var_id, obs_var(1:n_obs), start=(/1/), count=(/n_obs/))
      if (fu_fails(stat == NF90_NOERR, 'Failed to write obs var', sub_name)) return

    end subroutine init_nc

    !========================================================================

    subroutine column_to_nc(ncid, nc_var_id, ind_contr, values, num_values)
      implicit none
      integer, intent(in) :: nc_var_id, ind_contr
      integer, dimension(:), intent(in) :: ncid
      real, dimension(:), intent(in) :: values
      integer, dimension(:), intent(in) :: num_values

      integer :: stat, iPurpose, iStart
      
      iStart = 0
      do iPurpose = 1, 2
        call da_msg('Column mean for:' + chObsPurpose(iPurpose) + ':', &
                  & sum(values(iStart+1 : iStart+num_values(iPurpose))) / num_values(iPurpose))
        stat = nf90_put_var(ncid(iPurpose), nc_var_id, values(iStart+1:iStart+num_values(iPurpose)), &
                          & start=(/1, ind_contr/))
        if (fu_fails(stat == NF90_NOERR, 'Failed to write values', sub_name)) return
        iStart = num_values(iPurpose)
      enddo

    end subroutine column_to_nc
    
  end subroutine run_h_matrix


  !*****************************************************************************

  subroutine steepest(model, model_adj, control, background, bgr_cov, report_iterates)
    implicit none
    type(model_container), intent(inout) :: model, model_adj
    type(da_control), intent(inout) :: control
    type(da_control), intent(in) :: background
    type(t_background_covariance), intent(in) :: bgr_cov
    logical, intent(in) :: report_iterates

    type(da_control) :: control_linesearch, gradient
    integer :: ii, iter, iline, n_obs_values, control_size, outputlev, iPurpose, itercount
    real, dimension(2) :: cost_bgr, cost_obs
    real :: norm, cost_try, delta, delta_old, cost_old, grad_norm_first, cost, largest_change, cost_first
    logical :: have_moved
    type(DA_numericalParameters) :: numerics
    character(len=128) :: strCount   !log_str, 
    character(len=fnlen) :: output_dir, chOutFNm
    real, dimension(:), pointer :: grad_vals, control_linsrch_vals, control_vals, mdl_obs_values, &
                                 & mdl_obs_values_prev

    outputlev = model%rules%darules%output_level

    numerics = model%rules%darules%numerics
    output_dir = model%rules%darules%outputdir
    if (outputlev > da_out_none_flag) then 
      if (fu_fails(output_dir /= char_missing, 'outputdir not set', 'steepest')) return
    endif

    call init_control_from_cloud(control_linesearch, model%cloud, model%rules%darules, control_space, 0.0, bgr_cov)
    call init_control_from_cloud(gradient, model%cloud, model%rules%darules, control_space, 0.0, bgr_cov)

    mdl_obs_values => fu_work_array(sum(model%obs_ptr%obs_size))
    mdl_obs_values_prev => fu_work_array(sum(model%obs_ptr%obs_size))

!    call collect_variance(model%obs_ptr, mdl_obs_values) !, n_obs_values) ! just to get n_obs_values
    n_obs_values = sum(model%obs_ptr%obs_size)
    
    !**********************************************************************************
    !
    ! The first forward run
    !
    !**********************************************************************************

    if (outputlev > da_out_none_flag) then
      chOutFnm = fu_next_iter_file(output_dir, 'control_da_', model%rules%darules%da_begin, 1) 
      call control_to_file(control, background, model%rules%darules, chOutFnm)
    end if
    strCount = 'initial forward run'
    call start_count(chCounterNm = strCount)

    call evaluate(control, background, bgr_cov, model, cost_obs, cost_bgr, cost, mdl_obs_values)
    if (error) return
    if (outputlev > da_out_none_flag) then
      chOutFNm = fu_next_obs_file(output_dir, model%rules%darules%da_begin, 1)
      do iPurpose = 1, 2
        call dump_observations(model%obs_ptr, fu_species_transport(model%cloud), &
                             & chOutFNm + chObsPurpose(iPurpose), iPurpose)
      end do
    endif
    if (error) return

    call stop_count(chCounterNm = strCount)

    call da_msg('===========================================================')
    call da_msg('Initial forward run complete')
    call report_cost(0, cost_obs, cost_bgr, cost, fu_mode(control)) !, log_str)
!    call da_msg(log_str)
    if (da_test) return
    delta = -1.0
    cost_first = cost

    grad_vals => fu_values_ptr(gradient)
    control_linsrch_vals => fu_values_ptr(control_linesearch)
    control_vals => fu_values_ptr(control)

    control_size = fu_control_space_dim(bgr_cov)

    strCount = 'gradient search'
    call start_count(chCounterNm = strCount)

    itercount = 0
    GRADIENT_SEARCH: do iter = 1, numerics%maxIterations
      itercount = itercount + 1
      have_moved = .false.

      if (defined(model%tla_traj) .and. delta > 0) then
        call da_msg('Re-evaluate forward for tla trajectory')
        call evaluate(control, background, bgr_cov, model, cost_obs, cost_bgr, cost, mdl_obs_values)
        if (error) return
      end if

      !call evaluate(control, background, bgr_cov, model, cost_obs, cost_bgr, cost, mdl_obs_values)
      call get_gradient(control, gradient, background, bgr_cov, model_adj, mdl_obs_values)
      if (error) return
      if (report_iterates) then
        call msg('Reporting gradient')
        call report(gradient)
      end if
      if (outputlev > da_out_first_last_flag) then
        chOutFnm = fu_next_iter_file(output_dir, 'gradient_da_', model%rules%darules%da_begin, iter)
        call control_to_file(gradient, background, model%rules%darules, chOutFnm)
      end if

      if (error) return

      norm = sqrt(sum(grad_vals(1:control_size)**2))
      if (iter == 1) grad_norm_first = norm

      call da_msg('===========================================================')
      call da_msg('Adjoint run complete, ||grad||: ', norm)

      if (norm .eps. 0.0) then
        call da_msg('Iteration finished, vanishing gradient')
        exit 
      else if (norm < grad_norm_first*numerics%grad_rel_tol) then
        call da_msg('Iteration finished, ||grad|| < a*||grad_0||, a = ', numerics%grad_rel_tol)
        exit
      end if

      grad_vals = grad_vals/norm

      if (delta < 0) then
        delta = cost / norm * 20
      else
        delta = delta / numerics%stepDecreaseCoefMS
      end if

      strCount = 'line_search'
      call start_count(chCounterNm = strCount)

      LINE_SEARCH: do iline = 1, numerics%maxLineSearchSteps
        itercount = itercount + 1
        call da_msg('Trying step: ', delta)
        do ii = 1, fu_size(control)
          control_linsrch_vals(ii) = control_vals(ii) - delta * grad_vals(ii)
        end do
        if (report_iterates) then
          call msg('Report control_linesearch')
          call report(control_linesearch)
        end if
        if (outputlev > da_out_first_last_flag) then
          chOutFnm = fu_next_iter_file(output_dir,  'control_da_', model%rules%darules%da_begin, itercount)
          call control_to_file(control_linesearch, background, model%rules%darules, chOutFnm)
          if (error) return
        end if

        call evaluate(control_linesearch, background, bgr_cov, model, &
                    & cost_obs, cost_bgr, cost_try, mdl_obs_values)
        if (error) return

        call report_cost(itercount, cost_obs, cost_bgr, cost_try, fu_mode(control)) !, log_str)
!        call da_msg(log_str)
        if (outputlev > da_out_first_last_flag)then
          chOutFNm = fu_next_obs_file(output_dir, model%rules%darules%da_begin, itercount)
          do iPurpose = 1, 2
            call dump_observations(model%obs_ptr, fu_species_transport(model%cloud), &
                                 & chOutFNm + chObsPurpose(iPurpose), iPurpose)
          end do
        endif
        if (error) return

        if (have_moved .and. cost_try < cost_old .or. .not. have_moved .and. cost_try < cost) then
          if (2 * abs(cost - cost_try) / (cost + cost_try) < numerics%cost_incr_rel_tol) then
            call da_msg('===========================================================')
            call da_msg('Iteration finished: 2 * (cost - cost_try) / (cost + cost_try) < a, a = ', &
                      & numerics%cost_incr_rel_tol)
            exit GRADIENT_SEARCH
          else if (cost_try < numerics%cost_rel_tol * cost_first) then
            call da_msg('===========================================================')
            call da_msg('Iteration finished: J(u) < a*J_0, a = ', numerics%cost_rel_tol)
            exit GRADIENT_SEARCH
          end if

          delta_old = delta
          delta = delta * numerics%stepIncreaseCoefMS
          cost_old = cost_try
          mdl_obs_values_prev(1:n_obs_values) = mdl_obs_values(1:n_obs_values)
          !call getModelDataContainerArray(model%obs_ptr, arrayMDC)
          have_moved = .true.
          cycle LINE_SEARCH

        else
          ! Not improving. If first linesearch step, try shorter
          ! delta. Otherwise revert to the last improving value and
          ! return to gradient search.
          if (.not. have_moved) then
            delta = delta / numerics%stepDecreaseCoefMS
            if (delta < numerics%minimumStepLength) then
              call da_msg('===========================================================')
              call da_msg('Iteration finished: step length below', numerics%minimumStepLength)
              exit GRADIENT_SEARCH
            end if
            cycle LINE_SEARCH
          else
            delta = delta_old
            call set_mdl_data(model%obs_ptr, mdl_obs_values_prev)
            control_vals = control_vals - delta*grad_vals
            cost = cost_old
!!$            if (largest_change < numerics%minmax) then
!!$              call da_msg('===========================================================')
!!$              call da_msg('Iteration finished: max(abs(d*grad/(x + d*grad))) < ',numerics%minmax)
!!$              exit GRADIENT_SEARCH
!!$            end if
            exit LINE_SEARCH
          end if
        end if

      end do LINE_SEARCH

      call stop_count(chCounterNm = strCount)

    end do GRADIENT_SEARCH

    call free_work_array(mdl_obs_values)
    call free_work_array(mdl_obs_values_prev)

    strCount = 'gradient search'
    call stop_count(chCounterNm = strCount)

  end subroutine steepest

  !**********************************************************************************

  subroutine quasi_newton(model, model_adj, control, background, bgr_cov, report_iterates)
    implicit none
    type(model_container), intent(inout) :: model, model_adj
    type(da_control), intent(inout) :: control
    type(da_control), intent(in) :: background
    type(t_background_covariance), intent(in) :: bgr_cov
    logical, intent(in) :: report_iterates

    real :: refcost
    integer :: stat, iter, iz(5), nsubiter = 3000, iterstat, &
         & nsim = 2500, icom, niter = 2500, imode(2)
    !integer :: optimization_problem_size, worksize_optimizer
    real, dimension(1) :: rzs, rcomarr
    real(optimizer_precision), dimension(1) :: dzs
    integer, dimension(1) :: izs
    real(optimizer_precision) :: epsg, df1

    real, dimension(2) :: evalstat
    real, dimension(:,:),allocatable ::  cost_bgr,  cost_obs, control_vals_store
    real, dimension(:), allocatable :: costs, gradnorms

    integer :: contr_size, worksize_optimizer
    real(optimizer_precision), dimension(:), allocatable :: control_dp, gradient_dp, optimizer_work
    real, dimension(:), pointer :: mdl_obs_values, control_values, gradient_values
    real(optimizer_precision) :: cost_dp
    type(DA_numericalParameters) :: numerics
    type(da_control) :: gradient
    character(len=fnlen) :: output_dir, chFnm
    integer :: outputlev, impres, iPurpose, best_iter, iTmp

    logical :: ifMaster
    integer (kind=8):: count_start,count1, count2, eval_cnt,getgrad_cnt,copy2dp_cnt,copy2sp_cnt,opt_cnt,gradnorm_cnt
    character(len=*), parameter :: sub_name = 'quasi_newton'

     eval_cnt=0; getgrad_cnt=0; copy2dp_cnt=0; copy2sp_cnt=0; opt_cnt=0;
     gradnorm_cnt=0;
     CALL SYSTEM_CLOCK(count_start)


    best_iter = int_missing

    ifMaster = (smpi_adv_rank==0)

    outputlev = model%rules%darules%output_level

    output_dir = model%rules%darules%outputdir
    !
    if (outputlev > da_out_none_flag .and. ifMaster) then
      if (fu_fails(output_dir /= char_missing, 'outputdir not set', 'quasi_newton')) return
    endif

    contr_size = fu_size(control)
    worksize_optimizer = contr_size*20
    
    numerics = model%rules%darules%numerics


    mdl_obs_values => fu_work_array(sum(model%obs_ptr%obs_size))
    control_values => fu_values_ptr(control)
    nsim = numerics%maxIterations
    niter = numerics%maxIterations
    nsubiter = numerics%maxIterations
    epsg = numerics%grad_rel_tol

    
    call msg("Allocating optimizer stuff for contr_size,nsim", contr_size,nsim)
    if(ifMaster) then 
      call da_msg("allocating control store, size: " &
                &//trim(fu_str(contr_size))//"x"//trim(fu_str(nsim+1)))
      allocate(control_dp(contr_size), gradient_dp(contr_size), optimizer_work(worksize_optimizer), &
          & control_vals_store(contr_size,0:nsim),  cost_bgr(2,0:nsim),  cost_obs(2,0:nsim), & 
          & costs(0:nsim), gradnorms(0:nsim), stat=iTmp)
    else 
      allocate(cost_bgr(2,0:nsim),  cost_obs(2,0:nsim),  costs(0:nsim), stat=iTmp)
    endif
    if (iTmp /= 0) then
        call set_error("Allocation failed..", sub_name)
        return
    endif

    call init_control_from_cloud(gradient, model%cloud, model%rules%darules, control_space, 0.0, bgr_cov)
    gradient_values => fu_values_ptr(gradient)
    call da_msg('Initial cost function evaluation')

    do iter = 0, niter
      
      CALL SYSTEM_CLOCK(count1)
      if (iter > 0 .and. ifMaster) control_values(:) = control_dp(:) !! Update control vector from optimizer
      if (ifMaster) control_vals_store(:,iter) = control_values(:) ! Store control vector
      CALL SYSTEM_CLOCK(count2)
      copy2sp_cnt = copy2sp_cnt + count2 - count1

      call evaluate(control, background, bgr_cov, model, cost_obs(:,iter), &
                     & cost_bgr(:,iter), costs(iter), mdl_obs_values)
      CALL SYSTEM_CLOCK(count1)
       eval_cnt =  eval_cnt + count1 - count2
      if (error) return
      call report_cost(iter, cost_obs(:,iter), cost_bgr(:,iter), costs(iter), fu_mode(control)) !, log_str)
!      call da_msg(log_str)
      if (outputlev > da_out_none_flag .and. ifMaster) then
        chFnm = fu_next_iter_file(output_dir, 'background_da_', model%rules%darules%da_begin, 0) !iter
        call control_to_file(control, background, model%rules%darules, chFnm)
      endif

      if (outputlev > da_out_none_flag) then !! dump_observations is a collective call
        do iPurpose =1, 2
          chFnm = fu_next_obs_file(output_dir, model%rules%darules%da_begin, iter) !obs
          call dump_observations(model%obs_ptr, fu_species_transport(model%cloud), chFnm + chObsPurpose(iPurpose), iPurpose)
        end do
      end if

      CALL SYSTEM_CLOCK(count1)
      call get_gradient(control, gradient, background, bgr_cov, model_adj, mdl_obs_values)
      CALL SYSTEM_CLOCK(count2)
      getgrad_cnt =  getgrad_cnt + count2 - count1
      if (error) return

      if (ifMaster) then

        if (outputlev > da_out_none_flag) then
          chFnm = fu_next_iter_file(output_dir, 'control_da_', model%rules%darules%da_begin, iter)
          call control_to_file(control, background, model%rules%darules, chFnm)
        end if


        CALL SYSTEM_CLOCK(count1)
        gradnorms(iter) = sqrt(dot_product(gradient_values, gradient_values)) 
        CALL SYSTEM_CLOCK(count2)
        gradnorm_cnt =  gradnorm_cnt + count2 - count1
        call da_msg('||grad|| = ', gradnorms(iter))
        if (outputlev > da_out_first_last_flag) then
          chfnm = fu_next_iter_file(output_dir, 'gradient_da_', model%rules%darules%da_begin, iter)
          call control_to_file(gradient, background, model%rules%darules, chFnm)
        end if

        if (iter==0) then
          df1 = costs(iter) * numerics%quasi_newton_df1
          call da_msg("Cost reduction guess", real(df1, kind=4))
          call da_msg("Cost reduction factor to initial cost",  numerics%quasi_newton_df1)

          CALL SYSTEM_CLOCK(count1)
          control_dp(:) = control_values(:)
          CALL SYSTEM_CLOCK(count2)
          copy2dp_cnt =  copy2dp_cnt + count2 - count1

          icom = 1
          imode = (/0,0/)
          impres = 5


          !  !!FIXME Have no clue what is this comment about
          ! Careful with FP constants - the interface is implicit!
          !
        else  
          icom = 2
          imode = (/0,1/)
          impres=2
        endif

        best_iter = minloc(costs(0:iter),1) - 1 !!! Minimal-cost iteration
        if (iter > 3) then 
          if (best_iter + 2 < iter) then !! Two iterations worse than old one
              icom = -1  !! Do not call m1qn3
              iterstat = 10 !! Termination status (not one from m1ql)
          endif
          if ((icom > 0))then
              refcost = 1.05*costs(best_iter) !! Last two within 5% of the best one
                                         !! otherwise no substantial improvement
              if (all(refcost > costs(iter-2:iter-1))) then !! 
                icom = -1  !! Do not call m1qn3
                iterstat = 11 !! Termination status (not one from m1ql)
                best_iter = iter-2  
              endif
          endif
        endif

        if (icom > 0) then  !! New iteration neded
          CALL SYSTEM_CLOCK(count1)
          gradient_dp = gradient_values
          CALL SYSTEM_CLOCK(count2)
          copy2dp_cnt =  copy2dp_cnt + count2 - count1
          cost_dp = costs(iter)

          CALL SYSTEM_CLOCK(count1)
          call m1qn3(simul_rc, euclid, ctonbe, ctcabe, & ! routines for optimization, use default
                 & contr_size, & ! problem dimension
                 & control_dp, & ! initial value
                 & cost_dp, &    ! cost function 
                 & gradient_dp, & ! gradient
                 & 1d-6, &      ! dxmin
                 & df1, &       ! df1 !!Expected cost function reduction in _absolute_ values
                 & epsg, &      ! epsg
                 & 'two', &     ! normtype
                 & impres, &    ! impres (printing)
                 & 6, &         ! io channel for printing
                 & imode, &     ! imode: scalar updating, warm start
                 & iterstat, &  ! omode: state on return
                 & nsubiter, &  ! niter
                 & nsim, &      ! number of function evaluations
                 & iz, &        ! address array
                 & optimizer_work, & ! work array
                 & worksize_optimizer, & ! size of work (ndz)
                 & .true., &    ! use reverse communication
                 & icom, &      ! (indic) communication indicator
                 & izs,rzs,dzs) ! working area addresses, not used in reverse communication mode
          CALL SYSTEM_CLOCK(count2)
          opt_cnt =  opt_cnt + count2 - count1
        endif
      endif !! ifMaster
    
      ! Bcast iteration status from master
      if (smpi_adv_tasks > 1) then
        evalstat(1) = real(icom)
        evalstat(2) = real(iterstat)
        call smpi_bcast_aray(evalstat,2,smpi_adv_comm,0)
        if (smpi_adv_rank > 0) then
          icom = nint(evalstat(1))
          iterstat = nint(evalstat(2))
        endif
      endif

       !Check after initial only 
      if (fu_fails(iter>0 .or. icom == 4, 'Minimizer had an error', 'optimize')) return

      if (icom <= 0) then
        select case(iterstat)
          case(1)
            call da_msg('Normal termination')
          case(2)
            call set_error('Error in m1qn3 input', 'main')
          case(3,4,5,6,7)
            call da_msg('Iterstat = ', ival=iterstat)
            call da_msg('Abnormal termination, but continuing')
          case(10)
            call da_msg('No more improvement. Breaking now.')
          case(11)
            call da_msg('Too slow improvement. Breaking now.')
          case(0) 
            call da_msg('Termination with omode = 0')
          case default
            call set_error('Strange m1qn3 omode', 'quasi_newton')
        end select
        exit
      end if
    end do

    call da_msg('Best iteration found   = ', ival=best_iter)
    control_values(:) = control_vals_store(:,best_iter) !!restore best_iter control
   

    CALL SYSTEM_CLOCK(count1)
    call da_msg('Iterations over, gradient evaluations: '+ fu_str(iter+1))
    call msg("Counters: eval_cnt,getgrad_cnt,copy2dp_cnt,copy2sp_cnt,opt_cnt,gradnorm_cnt, total")
    call msg ("ms",int(1e-6*(/eval_cnt,getgrad_cnt,copy2dp_cnt,copy2sp_cnt,opt_cnt,gradnorm_cnt, count1-count_start/)))
    
    if (error) return 
    if (iter > niter) then
      call da_msg('Iteration terminated: maximum number of iterations')
    end if
    call free_work_array(mdl_obs_values)
    call destroy(gradient)

    if(ifMaster) then 
      deallocate(control_dp, gradient_dp, optimizer_work,  control_vals_store,  cost_bgr,  cost_obs, & 
          & costs, gradnorms)
    else 
      deallocate(cost_bgr,  cost_obs,  costs)
    endif

    call da_msg('quasi_newton finished, iter',ival=iter)
  end subroutine quasi_newton


  !**********************************************************************************************
  
  subroutine l_bfgs_b(model, model_adj, control, background, bgr_cov, report_iterates)
    !
    ! Runs BFGS assimilation algorithm for 4D-VAR
    !
    use bfgs
    implicit none
    
    ! Imported parameters
    type(model_container), intent(inout) :: model, model_adj
    type(da_control), intent(inout) :: control
    type(da_control), intent(in) :: background
    type(t_background_covariance), intent(in) :: bgr_cov
    logical, intent(in) :: report_iterates

    ! Local variables
    integer :: best_iter, eval_count
    integer, parameter :: max_iter=100
    integer :: iterCost, outputlev, iPurpose
    !integer :: optimization_problem_size, worksize_optimizer
    real, dimension(1) :: rzs
    real(optimizer_precision), dimension(1) :: dzs
    integer, dimension(1) :: izs

    integer :: contr_size, worksize_optimizer, unitPrint
    integer, parameter :: num_corrects = 15, print_opt = 40
    real(optimizer_precision), dimension(:), allocatable :: control_dp, gradient_dp, optimizer_work, &
         & lower_bounds, upper_bounds
    real, dimension(:), pointer :: mdl_obs_values, control_values, gradient_values
    real(optimizer_precision) :: cost_dp
    real :: gradient_norm, gradient_norm_first, xscale
    character(len=60) :: task, l_bfgs_b_csave
    real(optimizer_precision), parameter :: l_bfgs_b_factr = 0, l_bfgs_b_pgtol = 0 !! Sippres internal BFGS criteria
    real(optimizer_precision), dimension(29) :: l_bfgs_b_dsave
    integer, dimension(:), allocatable :: l_bfgs_b_iwa, bound_def
    integer, dimension(44) :: l_bfgs_b_isave
    logical, dimension(4) :: l_bfgs_b_lsave 
    type(DA_numericalParameters) :: numerics
    type(da_control) :: gradient
    real, dimension(2,0:max_iter) :: cost_bgr, cost_obs
    real, dimension(0:max_iter) :: cost
    character(len=fnlen) :: output_dir, chOutFNm   !log_str, 
    real, dimension(10) :: fit_params
    real, dimension(:,:), allocatable :: contol_vals_store
    character(len=*), parameter :: sub_name = 'l_bfgs_b'
    !
    ! Organise the DA output
    !
    outputlev = model%rules%darules%output_level

    if (outputlev >= da_out_first_last_flag) then 
      output_dir = model%rules%darules%outputdir
      if (fu_fails(output_dir /= char_missing, 'outputdir not set', sub_name)) return
      unitPrint = fu_next_free_unit()
      open(unitPrint, status = 'unknown', &
         & file = output_dir + dir_slash + 'optimization_summary_' + model%rules%chCaseNm + '_' + &
         & fu_time_fname_string_utc(model%rules%darules%da_begin) + '.txt', action='write')
    else
      unitPrint = run_log_funit
    endif
    !
    ! Prepare working structures
    !
    contr_size = fu_size(control)
    worksize_optimizer = (2*num_corrects + 5) * contr_size + 11*num_corrects**2 + 8 * num_corrects
    numerics = model%rules%darules%numerics
    best_iter = int_missing
    fit_params(:) = real_missing

    allocate(control_dp(contr_size), gradient_dp(contr_size), optimizer_work(worksize_optimizer), &
           & l_bfgs_b_iwa(3*contr_size), &
           & lower_bounds(contr_size), upper_bounds(contr_size), bound_def(contr_size))

    call da_msg("allocating control store, size: " &
                    &//trim(fu_str(contr_size))//"x"//trim(fu_str(max_iter+1)))

    allocate(contol_vals_store(contr_size,0:max_iter))

    mdl_obs_values => fu_work_array(sum(model%obs_ptr%obs_size))
    control_values => fu_values_ptr(control)

    call get_posit_constr(control, background, bgr_cov, lower_bounds, bound_def)
    if (error) return

    gradient_norm_first = -9999.0
    cost = 1e30  ! something big enough
    xscale = 100./sqrt(1.*contr_size) !! 
    control_dp(:) = control_values(:)*xscale !!!Copy to double precision array
    !lower_bounds = lower_bounds * xscale 

    call init_control_from_cloud(gradient, model%cloud, model%rules%darules, control_space, 0.0, bgr_cov)
    gradient_values => fu_values_ptr(gradient)
    
    task = 'START'

    call da_msg('To iteration loop...')

    eval_count = 0 
    do iterCost = 0, numerics%maxIterations
      !
      ! Run BFGS iteration
      !
      call setulb(contr_size, num_corrects, &
                & control_dp, lower_bounds, upper_bounds, bound_def, &
                & cost_dp, gradient_dp, l_bfgs_b_factr, l_bfgs_b_pgtol, &
                & optimizer_work, l_bfgs_b_iwa, task, print_opt, &
                & l_bfgs_b_csave, l_bfgs_b_lsave, l_bfgs_b_isave, l_bfgs_b_dsave, &
                & unitPrint)
      call da_msg('Task on exit from setulb:' // trim(task(1:5)))
      if (task(1:5) == 'NEW_X') then
        !! New_x is usually very close to the last-evaluated. 
        !! Not much sense to evaluate again,
        !! Re-evaluating the same vector triggers termination on too small change of cost or gradient
        !! So just re-launch it again
        call setulb(contr_size, num_corrects, &
                  & control_dp, lower_bounds, upper_bounds, bound_def, &
                  & cost_dp, gradient_dp, l_bfgs_b_factr, l_bfgs_b_pgtol, &
                  & optimizer_work, l_bfgs_b_iwa, task, print_opt, &
                  & l_bfgs_b_csave, l_bfgs_b_lsave, l_bfgs_b_isave, l_bfgs_b_dsave, &
                  & unitPrint)
        call da_msg('Task on exit from setulb on newX retry:' // trim(task(1:5)))
      endif

      !
      ! There are several possible regimes BFGS finishes with
      !
      if (task(1:2) == 'FG') then   ! <<<<<<-----------------------------------------------
        !
        ! Next state obtained, evaluate the costs to see if job done
        !
        if  (iterCost /= 0) then 
          control_values = control_dp / xscale !! extract to single-precision
        endif
        contol_vals_store(1:contr_size, iterCost) =  control_values(1:contr_size) !! Save control
        
        !
        ! If requested, dump control and observations.
        ! Note that minimization is done in control space x'=B^-1/2(x-xb), whereas storage is needed
        ! in physical space x = B^1/2 * x' + xb. Sizes are the same but meaning differs.
        !
        if (outputlev > da_out_first_last_flag .or. &
          & outputlev == da_out_first_last_flag .and. iterCost == 1) then
          ! physical space - the solution
          chOutFNm = fu_next_iter_file(output_dir, 'control_da_', model%rules%darules%da_begin, iterCost)
          call control_to_file(control, background, model%rules%darules, chOutFNm)
          ! Observations
          chOutFNm = fu_next_obs_file(output_dir, model%rules%darules%da_begin, iterCost)
          do iPurpose = 1, 2
            call dump_observations(model%obs_ptr, fu_species_transport(model%cloud), &
                                 & chOutFNm  + chObsPurpose(iPurpose), iPurpose)
          end do
        end if
        if(error)return

        
!        if (task(1:2) == 'FG') then   ! <<<<<<-----------------------------------------------
          !! New evaluation and gradient requested -- calculate it!
          call evaluate(control, background, bgr_cov, model, &
                      & cost_obs(:,iterCost), cost_bgr(:,iterCost), cost(iterCost), mdl_obs_values)
          if (error) return
          !
          ! The forward run has finished, report costs for assimilated sites and, if any, evaluation ones
          !
          call report_cost(iterCost, cost_obs(:,iterCost), cost_bgr(:,iterCost), cost(iterCost), fu_mode(control))
          !
          ! Next round: get the gradient
          !
          call get_gradient(control, gradient, background, bgr_cov, model_adj, mdl_obs_values)
          if (error) return
          
          !!! Prepare for next iteration
          cost_dp = cost(iterCost)
          gradient_dp = gradient_values / xscale 
          eval_count = eval_count + 1
        
!        else !!! New X
!          cost(iterCost) = cost_dp !! Projcted cost
!          call da_msg("Iter   "//trim(fu_str(iterCost))//" Projected cost = ", cost(iterCost) )
!          gradient_values = gradient_dp(:)  !! projected gradient
!          cost_obs(:,iterCost) = cost_obs(:,iterCost-1) ! Placeholder costs not to break l-curve fit
!          cost_bgr(:,iterCost) = cost_bgr(:,iterCost-1) !
!        endif

        gradient_norm = sqrt(sum(gradient_values**2))
        call da_msg('||grad|| = ', gradient_norm)
        !
        ! Should we break the iterations?
        !
        !
        ! Check if L-curve computations suggest to end. For now just evaluates
        ! Should be done before other breaks, so best_iter is valid
        !
        if(iterCost > 5)then
          !! OBS! We do not evaluate 0th iteration here, then best_iter index is the right one
          call best_iter_L_curve(iterCost, cost_obs(1,1:iterCost), cost_bgr(1,1:iterCost)+cost_bgr(2,1:iterCost), &
                               & best_iter, fit_params)
          call da_msg('Best iteration this-far found:' + fu_str(best_iter))
        endif
        if(error)return

        if (outputlev > da_out_first_last_flag) then
          chOutFNm = fu_next_iter_file(output_dir, 'gradient_da_', model%rules%darules%da_begin, iterCost)
          call control_to_file(gradient, background, model%rules%darules, chOutFNm)
        end if

        if (gradient_norm_first < 0) then 
          gradient_norm_first = gradient_norm
        else if (gradient_norm < numerics%grad_rel_tol*gradient_norm_first) then
          call da_msg('Iteration finished FG, ||grad|| < a*||grad_0||, a = ', numerics%grad_rel_tol)
          exit         
        end if

        !! Relative tolerance of cost
        if (cost(iterCost) < numerics%cost_rel_tol * cost(0)) then
          call da_msg('Iteration finished NEW_X: J(u) < a*J_0, a = ', numerics%cost_rel_tol)
          exit
        endif

        ! Cost increment 
        if (iterCost > 5)then
          if(abs(cost(iterCost) - cost(iterCost-1)) &
           & < numerics%cost_incr_rel_tol * 0.5*(cost(iterCost) + cost(iterCost))) then
            call da_msg('Brekaing iterations: 2 * (cost - cost_prev) / (cost + cost_prev) < a, a = ', &
                      & numerics%cost_incr_rel_tol)
            call msg('cost, cost_prev', cost(iterCost), cost(iterCost-1))
            exit
          endif
        end if

      elseif (task(1:5) == 'NEW_X') then
        call set_error("NEW_X returned on retry", sub_name)
        return
      else if (task(1:4) == 'CONV') then    ! <<<<<<-----------------------------------------------
        call da_msg('Normal termination CONV')
        exit
      else if (task(1:4) == 'ABNO') then    ! <<<<<<-----------------------------------------------
        call set_error('Abnormal termination ABNO', sub_name)
        return
      else if (task(1:5) == 'ERROR') then   ! <<<<<<-----------------------------------------------
        call da_msg('L-BFGS-B had an error')
        call set_error('L-BFGS-B had an error', sub_name)
        return
      else                                  ! <<<<<<-----------------------------------------------
        call set_error('Strange task:' + task, sub_name)
        return
      end if
    end do  ! iterations

    if (task(1:2) /= 'FG')  iterCost = iterCost - 1  !!Last iteration not evaluated
    
    call da_msg('Iteration:' + fu_str(iterCost) + ', out of:' + fu_str(numerics%maxIterations))
    call da_msg('Iterations over, gradient evaluations: '+ fu_str(eval_count))

    if (best_iter > 0) then
      call da_msg('Applying iteration (L-curve):' + fu_str(best_iter))
    else
      if (iterCost > 0) then
        best_iter  = minloc(cost(0:iterCost), 1)-1 !! Not enough iters for L-curve: use the smallest
        call da_msg('Applying iteration (cost):' + fu_str(best_iter))
      else
        best_iter = 0  !! Use background 
        call da_msg("Fallback to  background")
      endif
    endif
      !! Restore the control
    control_values(1:contr_size) = contol_vals_store(1:contr_size, best_iter)

    !Cleanup
    deallocate(contol_vals_store, control_dp, gradient_dp, optimizer_work, l_bfgs_b_iwa, &
            & lower_bounds, upper_bounds, bound_def)
    call free_work_array(mdl_obs_values)
    if ( .not. any(unitPrint == (/run_log_funit, 6/)))    close(unitPrint) !! Close if it was a file

  end subroutine l_bfgs_b

  
  !****************************************************************************
  
  subroutine report_cost(itNo, cost_obs, cost_bgr, cost_total, mode) !, log_str)
    implicit none
    integer, intent(in) :: itNo
    real, intent(in) :: cost_total
    real, dimension(:), intent(in) :: cost_obs, cost_bgr
    integer, intent(in) :: mode
!    character(len=*), intent(out) :: log_str

    character(len=clen) :: fmt, log_str
    !
    ! Basic report from assimilation
    !
    select case(mode)
    case (DA_INITIAL_STATE)
      fmt='("Iter ",I3,": [obs] + [b_init] = ", G12.6, " + ", G12.6, " = ", G12.6)'
      write(log_str, fmt=fmt)  itNo, cost_obs(1), cost_bgr(1), cost_total

    case (DA_EMISSION_CORRECTION, DA_EMISSION_TIME_HEIGHT)
      fmt = '("Iter ",I3,": [obs] + [b_emis] = ", G12.6, " + ", G12.6, " = ", G12.6)'
      write(log_str, fmt=fmt) itNo, cost_obs(1), cost_bgr(1), cost_total

    case (DA_EMISSION_AND_INITIAL)
      fmt = '("Iter ",I3,"[obs] + [b_emis] + [b_init] = ", 2(G12.6, " + "), G12.6, " = ", G12.6)'
      write(log_str, fmt=fmt) itNo,  cost_obs(1), cost_bgr(1), cost_bgr(2), cost_total

    case default
      call set_error('Unsupported control mode', 'report_cost')
      return
    end select
    
    call da_msg(log_str)
    !
    ! Auxiliary report, from evaluation
    !
    write(log_str, fmt='(A, G12.6)') '[obs_eval] = ', cost_obs(2)
    call da_msg(log_str)

  end subroutine report_cost

  !**********************************************************************************

  function fu_next_iter_file(dir, pref, da_begin, iterNo) result(filename)
    implicit none
    character(len=*), intent(in) :: dir
    character(len=*), intent(in) :: pref
    type(silja_time), intent(in) :: da_begin
    integer, intent(in) :: iterNo

    character(len=fnlen) :: filename
    ! replacement for thre functions:
    !fu_next_iter_file 'control_da_'
    !fu_next_grad_file 'gradient_da_'
    !fu_next_obs_file  'obs_da_'


    filename = dir + dir_slash + pref + &
             & fu_time_fname_string_utc(da_begin) + '_' + fu_str(iterNo,2) + '.grads'

    call msg('Next iter file: ' // trim(filename))
  end function fu_next_iter_file

  !**********************************************************************************

  function fu_next_obs_file(dir, da_begin, iterNo) result(filename)
    implicit none
    character(len=*), intent(in) :: dir
    type(silja_time), intent(in) :: da_begin
    integer, intent(in) :: iterNo

    character(len=fnlen) :: filename

    if (smpi_adv_tasks > 1) then
      filename = dir + dir_slash + 'obs_da_' + &
             & fu_time_fname_string_utc(da_begin) + '_' + fu_str(iterNo, 2)+ '_MPI'+ fu_str(smpi_adv_rank, 2)
    else
      filename = dir + dir_slash + 'obs_da_' + &
             & fu_time_fname_string_utc(da_begin) + '_' + fu_str(iterNo, 2)
    endif
    call msg('Next obs file: ' // trim(filename))
  end function fu_next_obs_file
  
  !**********************************************************************************
  
  real function fu_largest_relative_change(control_vals, gradient_vals, delta, n) result(largest)
    implicit none
    real, dimension(:), intent(in) :: control_vals, gradient_vals
    real, intent(in) :: delta
    integer, intent(in) :: n

    integer :: ii
    real :: num

    largest = 1e35

    do ii = 1, n
      num = control_vals(ii) + delta*gradient_vals(ii)
      if (num .eps. 0.0) cycle
      num = abs(delta*gradient_vals(ii) / num)
      if (num > largest) largest = num
    end do

  end function fu_largest_relative_change

  !************************************************************************************

  subroutine evaluate(control, background, bgr_cov, model, cost_obs, cost_bgr, cost_total, &
                    & mdl_obs_values)
    !
    ! Runs the forward model and gets its predictions at the observation points
    !
    implicit none
    type(da_control), intent(inout) :: control
    type(da_control), intent(in) :: background
    type(t_background_covariance), intent(in) :: bgr_cov
    type(model_container), intent(inout) :: model
    
    real, intent(out) :: cost_total
    real, dimension(:), intent(out) :: cost_obs, cost_bgr, mdl_obs_values
    
    real, dimension(:), pointer :: control_values, obs_values, variance, bgr_values
    integer :: ii,jj
    integer, dimension(2) :: n_obs

    if (.not. defined(control_physical)) then
      call init_control_from_cloud(control_physical, model%cloud, model%rules%darules, physical_space)
      if (error) return
    end if

    call to_physical(control, background, control_physical)
    !call report(control_physical)
    if (error) return
    n_obs(:) = model%obs_ptr%obs_size(:)    ! (assim, eval)
    if (fu_fails(size(mdl_obs_values) > sum(n_obs), 'mdl_obs_values too small', 'evaluate')) return
    !
    !------------------------------------------------------------------------------
    !
    ! Run the model and get observations from it
    !
    call obs_from_mdl(control_physical, model, mdl_obs_values)
    if (error) return
    call msg("mdl_obs_values ASSIM from model, nobs="//trim(fu_str(n_obs(1)))//": min, max, mean", &
           & (/minval(mdl_obs_values(1:n_obs(1))), maxval(mdl_obs_values(1:n_obs(1))), &
           & sum(mdl_obs_values(1:n_obs(1)))/))
    if (n_obs(2) > 0) then
      call msg("mdl_obs_values EVAL from model, nobs="//trim(fu_str(n_obs(2)))//" : min, max, mean", &
             & (/minval(mdl_obs_values(n_obs(1)+1:n_obs(1)+n_obs(2))), &
             & maxval(mdl_obs_values(n_obs(1)+1:n_obs(1)+n_obs(2))), &
             & sum(mdl_obs_values(n_obs(1)+1:n_obs(1)+n_obs(2)))/))
    else
      call msg("mdl_obs_values EVAL from model, nobs="//trim(fu_str(n_obs(2))) )
    endif

    !------------------------------------------------------------------------------
    !
    obs_values => fu_obs_data(model%obs_ptr)
!    call msg("obs_values min, max, mean", (/minval(obs_values(1:n_obs)), maxval(obs_values(1:n_obs)), sum(obs_values(1:n_obs))/))
    variance => fu_obs_variance(model%obs_ptr)
    if (fu_fails(associated(obs_values), 'obs_values not associated', 'evaluate')) return
    if (fu_fails(associated(variance), 'obs_variance not associated', 'evaluate')) return
    
    !!sum only stations within domain                                       
    ii = 1; jj = n_obs(1) !! Assimilation
    cost_obs(ii) = 0.5 * sum((obs_values(ii:jj) - mdl_obs_values(ii:jj))**2 / variance(ii:jj), &
                               & mdl_obs_values(ii:jj) /= real_missing)
    call mod_obs_stats("Assim", mdl_obs_values(ii:jj), obs_values(ii:jj), variance(ii:jj), .TRUE.)

    ii = n_obs(1) + 1; jj = n_obs(1) + n_obs(2) !! Validation
    cost_obs(2) = 0.5 * sum((obs_values(ii:jj) - mdl_obs_values(ii:jj))**2 / variance(ii:jj), &
                               & mdl_obs_values(ii:jj) /= real_missing)
    call mod_obs_stats("Valid", mdl_obs_values(ii:jj), obs_values(ii:jj), variance(ii:jj), .FALSE.)
    !
    ! total cost includes background cost and assimilated-observations deviation
    !
    cost_bgr = 0.0
    select case(fu_mode(control))
    case (DA_INITIAL_STATE, DA_EMISSION_CORRECTION, DA_EMISSION_TIME_HEIGHT)
      control_values => fu_values_ptr(control)
      cost_bgr(1) = 0.5 * sum(control_values(1:fu_size(control))**2)
      bgr_values => fu_values_ptr(background)
      cost_total = cost_obs(1) + cost_bgr(1)

    case (DA_EMISSION_AND_INITIAL)
      control_values => fu_initial_ptr(control)
      cost_bgr(2) = 0.5 * sum(control_values(1:fu_n_initial(control))**2)
      control_values => fu_emission_ptr(control)
      cost_bgr(1) = 0.5 * sum(control_values(1:fu_n_emission(control))**2)
      cost_total = sum(cost_bgr(1:2)) + cost_obs(1)

    case default
      call set_error('Unsupported control mode', 'evaluate')
      return
    end select
    !
    ! If no background, its costs are not included in the evaluation 
    !
    if (model%rules%darules%no_background_term) cost_total = cost_obs(1)

    if (debug_level >= 1) then
      call msg('Reporting control norm by species')
      call report_norm(control)
    end if

  end subroutine evaluate

  
  !************************************************************************************

  subroutine get_gradient(control, gradient, background, bgr_cov, model, mdl_obs_values)
    !
    ! Computes the gradient of the cost function via adjoint run through the assimilation window
    !
    implicit none
    type(da_control), intent(in) :: control, background
    type(da_control), intent(inout) :: gradient
    type(t_background_covariance), intent(in) :: bgr_cov
    type(model_container), intent(inout) :: model
    real, dimension(:), intent(in) :: mdl_obs_values

    integer :: idim, ii
    
    real, dimension(:), pointer :: grad_ptr, control_ptr
    integer(kind=8) :: count1, count2, count3

    if (.not. defined(control_physical)) then
      call init_control_from_cloud(control_physical, model%cloud, model%rules%darules, physical_space)
      if (error) return
    end if
    !call to_physical(control, background, control_physical)
    !if (error) return

    CALL SYSTEM_CLOCK(count1)
    call get_obs_grad(model, control_physical)
    CALL SYSTEM_CLOCK(count2)

    if (smpi_adv_rank > 0) return !! Only master gets proper control vector

    ! Gradient is now in physical space: M*H*(y-HMx)
    !
    if (error) return
    ! Apply the inverse square root of B:
    call to_control(control_physical, gradient)
    CALL SYSTEM_CLOCK(count3)

    call msg("get_obs_grad, to_control, ms", int(1e-6*(count2-count1)), int(1e-6*(count3 - count2)))

    if (debug_level > 0) then
      call msg('Reporting gradient norm by species before background')
      call report_norm(gradient)
    end if
    ! Add the background term, in control space
    !
    grad_ptr => fu_values_ptr(gradient)
    control_ptr => fu_values_ptr(control)
    if (model%rules%darules%log_transform) then
      grad_ptr = grad_ptr * exp(control_ptr)
    end if
   
    if (.not. model%rules%darules%no_background_term) then
      grad_ptr = grad_ptr + control_ptr
    end if
    if (debug_level > 0) then
      call msg('Reporting gradient norm by species after background')
      call report_norm(gradient)
    end if
  end subroutine get_gradient


  !**********************************************************************************
  !
  ! TEST CODES
  !
  !**********************************************************************************
  
  subroutine test_var_da(model)
    ! test the da subsystem, model must defined with some emission, etc.
    implicit none
    type(model_container), intent(inout) :: model
    ! test_init creates:
    type(DA_rules) :: rules_fake
    type(t_background_covariance) :: bgr_cov_fake 
    type(t_spatial_correlation) :: correlation

    call test_init_v2('test', fu_species_transport(model%cloud), fu_species_emission(model%cloud), &
                    & wholeMPIdispersion_grid, dispersion_vertical, correlation, bgr_cov_fake, rules_fake)

  end subroutine test_var_da

  subroutine gradient_test(model)
    implicit none
    type(model_container), intent(inout) :: model
    
    type(da_control) :: control, gradient, background, random
    type(t_background_covariance) :: bgr_cov
    real, dimension(:), pointer ::mdl_obs_values, initial_values, grad_values, rand_values
    real :: cost_obs(2), cost_bgr(2), cost, norm

    integer :: ix, iy, iz, isp, nx, ny, nz, nsp, ind, n_control
    type(Temission_processor), pointer :: proc_ptr

!    real :: delta = 1e1 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5
    real :: delta = 0.3e-1
    real :: fdiff, cost2, cost1
    real :: background_val = 0.0, background_init = 5e-6

    mdl_obs_values => fu_work_array(sum(model%obs_ptr%obs_size))
    
    ! Kill emissions
    !
    proc_ptr => fu_emission_processor_ptr(model%cloud)
    call set_emission_processor(fu_emisMM_ptr(model%cloud), fu_species_transport(model%cloud), &
                              & processor_scaling_flag, proc_ptr)
    call zero_processor(proc_ptr)

    call init_control_from_cloud(background, model%cloud, model%rules%darules, physical_space, background_init)
    
    call set_cov_mdl(bgr_cov, background, &
                   & fu_species_emission(model%cloud), fu_species_transport(model%cloud), &
                   & model%rules%darules, wholeMPIdispersion_grid, dispersion_vertical, &
                   & model%rules%darules%da_begin, model%rules%darules%controlVariable)
    call init_control_from_cloud(control, model%cloud, model%rules%darules, control_space, 0.0, bgr_cov)
    call init_control_from_cloud(random, model%cloud, model%rules%darules, control_space, 0.0, bgr_cov)
    call init_control_from_cloud(gradient, model%cloud, model%rules%darules, control_space, 0.0, bgr_cov)

    !model%rules%if_collect_linearization = .true.
    rand_values => fu_initial_ptr(random)
    call random_number(rand_values)
    n_control = fu_size(control) !fu_dimension_control(bgr_cov%correlation_initial)
    norm = sum(rand_values(1:n_control)**2)
    rand_values = rand_values / norm


    ix = 28
    iy = 23
    iz = 1
    isp = 1
    
    nx = nx_dispersion
    ny = ny_dispersion
    nz = nz_dispersion
    nsp = fu_size(control) / (nx*ny*nz) !fu_nbr_of_species_transport(model%cloud)

    !call msg('delta*1e-6*volume', delta*1e-6*fu_cell_size(wholeMPIdispersion_grid, ix, iy)*50)
    !call msg('bgr*1e-6*volume', background_val*1e-6*fu_cell_size(wholeMPIdispersion_grid, ix, iy)*50)
    !call msg('bgr+delta...', (background_val+delta)*1e-6*fu_cell_size(dispersion_grid, ix, iy)*50)
    ind = (iy-1)*nx*nz*nsp + (ix-1)*nz*nsp + (iz-1)*nsp + isp
    initial_values => fu_initial_ptr(control)
    initial_values(1:fu_n_initial(control)) = background_val

    ! Evaluate at -delta
    !initial_values(ind) = background_val - delta
    initial_values = background_val - delta*rand_values
    call evaluate(control, background, bgr_cov, model, cost_obs, cost_bgr, cost1, mdl_obs_values)

    ! Evaluate at delta
    !initial_values(ind) = background_val + delta
    initial_values = background_val + delta*rand_values
    call evaluate(control, background, bgr_cov, model, cost_obs, cost_bgr, cost2, mdl_obs_values)

    ! Gradient at zero
    initial_values(ind) = background_val
    call evaluate(control, background, bgr_cov, model, cost_obs, cost_bgr, cost, mdl_obs_values)
    call finalise_output(model%cloud, model%out_def, model%wdr, model%source)
    call get_gradient(control, gradient, background, bgr_cov, model, mdl_obs_values)
    !call to_physical(gradient, background, bgr_cov)
    call control_to_file(gradient, background, model%rules%darules, 'gradient')


    call msg('Cost1', cost1)
    call msg('Cost2', cost2)

    grad_values => fu_initial_ptr(gradient)
    call msg('physical:', size(grad_values))
    call msg('control:', n_control)
    call msg('max(abs(grad))', maxval(abs(grad_values(1:fu_n_initial(gradient)))))
    !call msg('grad', grad_values(ind))
    call msg('grad_proj', sum(grad_values(1:n_control)*rand_values(1:n_control)))
    fdiff = (cost2 - cost1) / (2*delta)
    
    call msg('fdiff', fdiff)
    !call msg('err', abs(fdiff-grad_values(ind)))
    stop
  end subroutine gradient_test

  subroutine gradient_test_emis(model)
    implicit none
    type(model_container), intent(inout) :: model
    
    type(da_control) :: control, gradient, background
    type(t_background_covariance) :: bgr_cov
    real, dimension(:), pointer ::mdl_obs_values, grad_values, proc_values
    real :: cost_obs(2), cost_bgr(2), cost

    integer :: ix, iy, iz, isp, nx, ny, nz, nsp, ind, iPurpose
    type(Temission_processor), pointer :: proc_ptr

!    real :: delta = 1e1 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5
    real :: delta = 0.1e2
    real :: fdiff, cost2, cost1
    real :: background_val = 0.05 ! standard deviations
    


    mdl_obs_values => fu_work_array(sum(model%obs_ptr%obs_size))
    
    proc_ptr => fu_emission_processor_ptr(model%cloud)
    call set_emission_processor(fu_emisMM_ptr(model%cloud), fu_species_transport(model%cloud), &
                              & processor_scaling_flag, proc_ptr)
    proc_values => fu_data(proc_ptr)
    proc_values = 1.0

    call init_control_from_cloud(control, model%cloud, model%rules%darules, control_space, 0.0)
    call init_control_from_cloud(gradient, model%cloud, model%rules%darules, physical_space)
    call init_control_from_cloud(background, model%cloud, model%rules%darules, physical_space, 1.0)
    
    call set_cov_mdl(bgr_cov, background, &
                   & fu_species_emission(model%cloud), fu_species_transport(model%cloud), &
                   & model%rules%darules, dispersion_grid, dispersion_vertical, &
                   & model%rules%darules%da_begin, model%rules%darules%controlVariable)
    !model%rules%if_collect_linearization = .false.

    ix = 23
    iy = 29
    isp = 1
    
    nx = nx_dispersion
    ny = ny_dispersion
    nz = nz_dispersion
    nsp = fu_nbr_of_species_transport(model%cloud)

    call msg('delta*1e-6*volume', delta*1e-6*fu_cell_size(dispersion_grid, ix, iy)*50)
    call msg('bgr*1e-6*volume', background_val*1e-6*fu_cell_size(dispersion_grid, ix, iy)*50)
    call msg('bgr+delta...', (background_val+delta)*1e-6*fu_cell_size(dispersion_grid, ix, iy)*50)
    ind = (iy-1)*nx*nsp + (ix-1)*nsp + isp

    proc_values => fu_emission_ptr(control)
    proc_values(1:fu_n_emission(control)) = background_val

    ! Evaluate at -delta
    proc_values(ind) = background_val - delta
    call evaluate(control, background, bgr_cov, model, cost_obs, cost_bgr, cost1, mdl_obs_values)
    do iPurpose = 1, 2
      call dump_observations(model%obs_ptr, fu_species_transport(model%cloud), &
                           & 'obs-delta' + chObsPurpose(iPurpose), iPurpose)
    end do
    ! Evaluate at delta
    proc_values(ind) = background_val + delta
    call evaluate(control, background, bgr_cov, model, cost_obs, cost_bgr, cost2, mdl_obs_values)
    do iPurpose = 1, 2
     call dump_observations(model%obs_ptr, fu_species_transport(model%cloud), &
                          & 'obs+delta' + chObsPurpose(iPurpose), iPurpose)
    end do
    ! Gradient at zero
    proc_values(ind) = background_val
    call evaluate(control, background, bgr_cov, model, cost_obs, cost_bgr, cost, mdl_obs_values)
    call finalise_output(model%cloud, model%out_def, model%wdr, model%source)
    call get_gradient(control, gradient, background, bgr_cov, model, mdl_obs_values)
    !call to_physical(gradient, background, bgr_cov)
    call control_to_file(gradient, background, model%rules%darules, 'gradient')


    call msg('Cost1', cost1)
    call msg('Cost2', cost2)

    grad_values => fu_emission_ptr(gradient)
    call msg('max(abs(grad))', maxval(abs(grad_values(1:fu_n_emission(gradient)))))
    call msg('grad', grad_values(ind))
    
    fdiff = (cost2 - cost1) / (2*delta)
    
    call msg('fdiff', fdiff)
    call msg('err', abs(fdiff-grad_values(ind)))
    stop
  end subroutine gradient_test_emis


  subroutine gradient_test_emis_2(model)
    ! 
    ! Yet another test routine.
    ! Here we'll check the following identity:
    !
    ! < H^T R^-1 (Hx - y), x > = < R^-1 (Hx - y), Hx >.
    !
    ! Here H is taken to mean the compound action of model and observation operator. The
    ! terms can be evaluated when Hx, y and the cost function gradient are available. The
    ! requirement is that model and observations are strictly linear (otherwise the
    ! tangent-linear model would be needed). 
    !
    ! This also assumes that there is no additional forcing besides the control variable:
    ! for emission_correction, initial condition must be void, for initial_state, there
    ! must be no emission (this is handled below). Also when the control variable does not
    ! include all species, the others must be zero. 
    ! 
    ! For the joint initial-emission mode, emission should exist, but initial condition of
    ! the non-control-species must be zero.
    ! 

    implicit none
    type(model_container), intent(inout) :: model
    
    type(da_control) :: control, gradient, background, control_phys
    type(t_background_covariance) :: bgr_cov
    real, dimension(:), pointer ::mdl_obs_values, grad_values, proc_values, bgr_values
    real :: cost_obs(2), cost_bgr(2), cost

    integer :: ix, iy, iz, isp, nx, ny, nz, nsp, ind, ii
    integer, dimension(2) :: n_obs_val
    type(Temission_processor), pointer :: proc_ptr
    real, dimension(:), pointer :: mdl_data, obs_data, obs_var, phys_values

!    real :: delta = 1e1 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5
    real :: delta = 0.1e5
    real :: fdiff, cost2, cost1
    real :: background_val = 0.05 ! standard deviations

    call random_seed(put=(/(12, ii=1,120)/))

    nx = nx_dispersion
    ny = ny_dispersion
    nz = nz_dispersion
    nsp = fu_nbr_of_species_transport(model%cloud)
    

    ix = 55
    iy = 3
    ind = (iy-1)*nx + ix

    mdl_data => fu_work_array()
    obs_data => fu_work_array()
    obs_var => fu_work_array()
    mdl_obs_values => fu_work_array()

    call init_control_from_cloud(background, model%cloud, model%rules%darules, physical_space, 0.0)
    call init_control_from_cloud(control_phys, model%cloud, model%rules%darules, physical_space) 
    
    call set_cov_mdl(bgr_cov, background, &
                   & fu_species_emission(model%cloud), fu_species_transport(model%cloud), &
                   & model%rules%darules, dispersion_grid, dispersion_vertical, &
                   & model%rules%darules%da_begin, model%rules%darules%controlVariable)
    call init_control_from_cloud(control, model%cloud, model%rules%darules, control_space, 0.0, bgr_cov)
    call init_control_from_cloud(gradient, model%cloud, model%rules%darules, control_space, 0.0, bgr_cov)

    proc_values => fu_values_ptr(control)
    call random_number(proc_values)

    if (fu_mode(control) == DA_INITIAL_STATE) then
      ! Kill emissions
      !
      proc_ptr => fu_emission_processor_ptr(model%cloud)
      call set_emission_processor(fu_emisMM_ptr(model%cloud), fu_species_transport(model%cloud), &
                                & processor_scaling_flag, proc_ptr)
      call zero_processor(proc_ptr)
    end if

    call evaluate(control, background, bgr_cov, model, cost_obs, cost_bgr, cost1, mdl_obs_values)
    call collect_model_data(model%obs_ptr, mdl_data) !, n_obs_val)
    call get_obs_data(model%obs_ptr,  obs_data, n_obs_val)
    call collect_variance(model%obs_ptr, obs_var) !, n_obs_val)
    ! can now compute R^{-1}(Hx-y)

    call get_gradient(control, gradient, background, bgr_cov, model, mdl_obs_values)
    grad_values => fu_values_ptr(gradient)
    !call copy_control(control, control_phys, .false.)
    !call to_physical(control_phys, background)
    call to_physical(control, background, control_phys)
    phys_values => fu_values_ptr(control_phys)

    !call msg('from proc:', proc_values(ind))
    !call msg('in grad:', grad_values(ind)-proc_values(ind))
    call msg('LEFT', sum((proc_values) * (grad_values-proc_values)))
    call msg('RIGHT', sum(mdl_data(1:sum(n_obs_val)) &
                        & * ((mdl_data(1:sum(n_obs_val))-obs_data(1:sum(n_obs_val))) / obs_var(1:sum(n_obs_val)))))

    print *, 'MDLDATA', mdl_data(1:sum(n_obs_val))
    print *, 'OBSDAT', obs_data(1:sum(n_obs_val))
    print *, 'OBSVAR', obs_var(1:sum(n_obs_val))
    print *, 'COUNT', count(abs(grad_values-proc_values) > 0)
    !call msg('LEFT', sum(grad_values))
    !call msg('RIGHT', sum(mdl_data(1:n_obs_val) &
                        !& * ((mdl_data(1:n_obs_val)-obs_data(1:n_obs_val)) / obs_var(1:n_obs_val))))

    stop
  end subroutine gradient_test_emis_2

  subroutine test_init_v2(test_file_name, species_transport, species_emission, grid, vertical, &
                        & correlation, bgr_cov, darules)
    use correlations
    implicit  none
    character(len=*), intent(in) :: test_file_name
    type(Tsilam_namelist), pointer :: nlptr_corr
    type(Tsilam_namelist_group), pointer :: p_nlgrp
    type(silam_species), dimension(:), intent(in), target :: species_transport, species_emission
    type(silja_grid), intent(in) :: grid
    type(silam_vertical), intent(in) :: vertical

    real, dimension(:), pointer :: wrk
    integer :: nx, ny, fs, nlevs, igf, ilev
    real :: stdev
    type(t_spatial_correlation) :: correlation
    type(silja_field_id) :: fid
    type(silam_species) :: species_so2, species_no2
    type(silja_time) :: sometime = ref_time_01012000
    character(len=fnlen) :: test_cov_ini_file
    type(silam_species), dimension(:), pointer :: p_species_transport, p_species_emission
    
    type(DA_rules) :: darules
    type(t_background_covariance) :: bgr_cov
    type(da_control) :: ctrl_background
    logical :: is_open
    integer :: file_unit

    p_species_transport => species_transport
    p_species_emission => species_emission

    !call set_species(species_so2, fu_get_material_ptr('SO2'), in_gas_phase)
    !if (error) return
    !call set_species(species_no2, fu_get_material_ptr('NO2'), in_gas_phase)
    !if (error) return

    !species_transport = (/species_so2, species_no2/)
    !species_emission = species_transport

    test_cov_ini_file = trim(test_file_name) + '.nl'
    file_unit = fu_next_free_unit()
    open(file_unit, file=test_cov_ini_file, action='write')
    call make_correlation(file_unit, emission_scaling_flag)
    call make_correlation(file_unit, concentration_flag)
    close(file_unit)
   
    call grid_dimensions(correlation%grid, nx, ny)
    fs = nx*ny
    nlevs = fu_nbrOfLevels(correlation%vertical)

    ! let's manufacture some variance data
    !
    call msg('Generating stdev data: emission')
    igf = open_gradsfile_o('', test_file_name, correlation%grid, &
        & time_label_position = instant_fields) !! Not exactly true, but shoudl hit timestamp match
    if (error) return
    wrk => fu_work_array()
    stdev = 0.5
    wrk(1:fs) = stdev
    call msg('SO2, stdev', stdev)

    fid = fu_set_field_id(met_src_missing, emission_scaling_flag, sometime, zero_interval, &
                        & correlation%grid, surface_level, species=species_transport(1))
    call write_next_field_to_gradsfile(igf, fid, wrk(1:fs))
    if (error) return
    stdev = 1.5
    wrk(1:fs) = stdev
    fid = fu_set_field_id(met_src_missing, emission_scaling_flag, sometime, zero_interval, &
                        & correlation%grid, surface_level, species=species_transport(2))
    call write_next_field_to_gradsfile(igf, fid, wrk(1:fs))
    if (error) return
    call msg('NO2, stdev', stdev)
    
    call msg('Generating stdev data: concentration')
    do ilev = 1, nlevs
      stdev = 1e-6 + ilev*1e-6
      wrk(1:fs) = stdev
      fid = fu_set_field_id(met_src_missing, concentration_flag, sometime, zero_interval, &
                          & correlation%grid, fu_level(correlation%vertical, ilev), &
                          & species=species_emission(1))
      call write_next_field_to_gradsfile(igf, fid, wrk(1:fs))
      if (error) return
      call msg('SO2, lev, stdev', ilev, stdev)
    end do
    do ilev = 1, nlevs
      stdev = 1e-5 + ilev*1e-5
      call msg('NO2, lev, stdev', ilev, stdev)
      fid = fu_set_field_id(met_src_missing, concentration_flag, sometime, zero_interval, &
                          & correlation%grid, fu_level(correlation%vertical, ilev), &
                          & species=species_emission(2))
      call write_next_field_to_gradsfile(igf, fid, wrk(1:fs))
      if (error) return
    end do
    
    call close_gradsfile_o(igf,"")

    darules%controlVariable = DA_EMISSION_AND_INITIAL
    !darules%background_var_mode_initial = constant_variance
    !darules%background_var_mode_emission = constant_variance
    !darules%subst_set_init = silja_false
    !darules%subst_set_emis = silja_false
    !call decode_template_string(test_file_name // '.grads.super_ctl', darules%template_var)
    call decode_template_string(test_cov_ini_file, darules%cov_setup_templ_emis)
    call decode_template_string(test_cov_ini_file, darules%cov_setup_templ_init)
    if (error) return
    
    call init_control_from_par(ctrl_background, correlation%grid, correlation%vertical, &
                    & p_species_emission, p_species_transport, darules, physical_space, 0.0)
    inquire(10, opened=is_open)
    print *, 'is_open', is_open

    call set_cov_mdl(bgr_cov, ctrl_background, species_emission, species_transport, &
                   & darules, correlation%grid, correlation%vertical, sometime, darules%controlVariable)
    
  contains
    
    subroutine make_correlation(file_unit, quantity_flag)
      implicit none
      integer, intent(in) :: file_unit
      integer, intent(in) :: quantity_flag
      
      type(Tsilam_namelist), pointer :: nlptr_corr
      logical :: surface_only

      call msg('Making correlation tests...')
      surface_only = quantity_flag /= concentration_flag

      call test_correlations(0.9, nlptr_corr, correlation, surface_only)
      !file_unit = fu_next_free_unit()
      !open(file_unit, file=filename_nl, action='write', form='formatted')
      call add_namelist_item(nlptr_corr, 'stdev_method', 'from_file')
      call add_namelist_item(nlptr_corr, 'stdev_file', test_file_name // '.super_ctl')
      call add_namelist_item(nlptr_corr, 'substance', 'SO2')
      call add_namelist_item(nlptr_corr, 'quantity', fu_quantity_short_string(quantity_flag))
      write(file_unit, *)'LIST = cov_so2'
      call write_namelist(file_unit, nlptr_corr, ifName=.false., ifPublic=.false.)
      write(file_unit, *)'END_LIST = cov_so2'

      call test_correlations(1.5, nlptr_corr, correlation, surface_only)
      call add_namelist_item(nlptr_corr, 'stdev_method', 'constant')
      call add_namelist_item(nlptr_corr, 'stdev_const', '5.67e-6')
      call add_namelist_item(nlptr_corr, 'substance', '*')
      call add_namelist_item(nlptr_corr, 'quantity', fu_quantity_short_string(quantity_flag))
      write(file_unit, *)'LIST = cov_def'
      call write_namelist(file_unit, nlptr_corr, ifName=.false., ifPublic=.false.)
      write(file_unit, *)'END_LIST = cov_def'
      
      !close(file_unit)
      
    end subroutine make_correlation

!!$    subroutine make_correlation(filename_nl)
!!$      implicit none
!!$      character(len=*), intent(in) :: filename_nl
!!$
!!$      type(Tsilam_namelist), pointer :: nlptr_corr
!!$      integer :: file_unit
!!$      
!!$      call msg('Making correlation tests...')
!!$
!!$      call test_correlations_grid_vert_given(0.9, nlptr_corr, correlation, grid, vertical)
!!$      file_unit = fu_next_free_unit()
!!$      open(file_unit, file=filename_nl, action='write', form='formatted')
!!$      call add_namelist_item(nlptr_corr, 'standard_deviation', '1.23e-6')
!!$      write(file_unit, *)'LIST = default'
!!$      call write_namelist(file_unit, nlptr_corr, ifName=.false., ifPublic=.false.)
!!$      write(file_unit, *)'END_LIST = default'
!!$
!!$      call test_correlations_grid_vert_given(1.5, nlptr_corr, correlation, grid, vertical)
!!$      call add_namelist_item(nlptr_corr, 'standard_deviation', '5.67e-6')
!!$      write(file_unit, *)'LIST = SO2'
!!$      call write_namelist(file_unit, nlptr_corr, ifName=.false., ifPublic=.false.)
!!$      write(file_unit, *)'END_LIST = SO2'
!!$      
!!$      close(file_unit)
!!$      
!!$    end subroutine make_correlation

  end subroutine test_init_v2

end module da_driver
