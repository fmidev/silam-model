module dispersion_supplementary

  ! This module was originally isolated from the dispersion_models module in order
  ! to allow the data_assimilation module to run dispersion simulations independently.
  ! Therefore, the subroutine run_dispersion basically contains the time loop of the simulation.
  ! Eventually, we could remove it entirely from dispersion_models and have it call this module.

!  use supermarket_of_fields
  use io_server
  use da_common
  use depositions
  use ini_boundary_conditions
!  use source_apportionment
  use perturbations

  !$ use OMP_LIB
  
  implicit none

  public run_dispersion

  ! This structure conveniently stores the different sets of rules needed in different
  ! parts of the model, as well as individual options and parameters. The implementation
  ! is not complete yet, much of the code in dispersion_models still does not use 
  ! this.

  ! Note that some of these have been moved from the wdr structure - namely, startingTime
  ! and periodToCompute.

  type general_dispersion_rules
     character(len=clen) :: chCaseNm
     type(DA_rules) :: daRules
     type(Tchem_rules) :: chemicalRules
     type(Tini_boundary_rules) :: IniBoundaryRules
     type(Tdiagnostic_rules) :: diagnosticRules
     type(Tdynamics_rules) :: dynamicsRules
     type(silja_time) :: startTime
     type(silja_interval) :: timestep, residenceInterval, periodToCompute
     integer :: iComputationAccuracy
     logical :: ifWriteProgressFile
     character(len=fnlen) :: chProgressFNmTemplate
     ! debug only, will remove, maybe
     logical :: ifRunDispersion
     logical :: if_invert_substeps = .false.
     ! Finalize output etc -- normally true but false for 3d-var integrations, which
     ! continue from where they left.
     logical :: if_finalize = .true.
     ! Whether to write and collect output, false for 4dvar iterations
     logical :: if_make_output = .true.
  end type general_dispersion_rules

  
contains 

  !*************************************************************************

  subroutine run_dispersion(cloud, em_source, wdr, simrules, &
                          & meteo_input_dyn_shopping_list, meteo_full_dyn_shopping_list, &
                          & disp_dyn_shopping_list, disp_stat_shopping_list, &
                          & outDef, &
                          & meteo_ptr, disp_buf_ptr, obs_ptrs, pMeteo_input, output_buf_ptr, &
                          & meteoMarketPtr, dispersionMarketPtr, BCMarketPtr, outputMarketPtr, &
                          & tla_traj, perturbations)

    ! This is the main time-integration cycle, which can be forward or adjoint. It pushes
    ! the run from the start to end time using the given structures and rules, which all
    ! must be initialised. 
    ! The only exception is the low-mass threshold, which is decided here.
    ! It is called either from dispersion_models after the initilization or from data assimilation
    ! modules. In the later case there may be several calls, both forward and adjoint in various
    ! combinations.
    !
    implicit none
    
    ! Imported parameters
    !
    type(silam_pollution_cloud), pointer :: cloud
    type(general_dispersion_rules), intent(inout) :: simrules
    type(silja_wdr), pointer :: wdr
    type(silja_shopping_list), intent(inout) :: meteo_input_dyn_shopping_list, &
                                              & meteo_full_dyn_shopping_list, &
                                              & disp_dyn_shopping_list, disp_stat_shopping_list
    type(silam_output_definition), pointer :: outDef
    type(silam_source), pointer :: em_source
    type(Tfield_buffer), pointer :: meteo_ptr, disp_buf_ptr, output_buf_ptr
    type(Tmeteo_input), pointer :: pMeteo_input
    type(observationPointers), intent(inout) :: obs_ptrs
    type(mini_market_of_stacks), pointer :: meteoMarketPtr, dispersionMarketPtr, &
                                          & BCMarketPtr, outputMarketPtr
    type(t_tla_trajectory), intent(inout), optional :: tla_traj
    type(t_perturbation), dimension(:), pointer, optional :: perturbations

    ! Local variables
    !
    real ::  fTmp
    type(silja_time) :: now
    real, dimension(max_species) ::  fInjectedMassTotal 
    real(r8k), dimension(max_species) :: fInjectedMass 
    CHARACTER (LEN=fnlen) :: command_string = ' '
    integer :: nParticlesToReset, iTmp, jTmp, num_steps, ind_step, step_count, number_of_times
    logical :: first_step,  ifGotNewMeteoData, have_tla_traj, have_perturbs, ifCalculate,&
         & first_output_step
    type(silja_shopping_list), pointer :: pShpLst, pOutput_dyn_shopping_list, pOutput_stat_shopping_list
    type(Tmass_map), pointer :: pOptDepth, pTranspMass, pTranspXm,  pTranspYm,  pTranspZm
!    type(silja_time), save :: nwp_lim1, nwp_lim2, nwp_lim1_prev, nwp_lim2_prev, &
!         & bc_lim1, bc_lim2, bc_lim1_prev, bc_lim2_prev
    type(silja_time) :: reftime
    type(silja_time), dimension(max_times) :: valid_times
    type(silja_interval) :: periodTmp
    type(silam_species), dimension(:), pointer :: pSpecies
    integer, save :: call_counter = 0
    logical, save :: have_done_output = .false.
    type(Tmoment_mapping) :: cm_to_moment
    type(t_tla_step) :: tla_step
    type(Tmass_map), pointer :: null_pointer
    type(silam_trajectory_set), pointer :: traj_set
    character (len=*), parameter :: sub_name = "run_dispersion"

    nullify(null_pointer)
    call_counter = call_counter + 1

    if (present(tla_traj)) then
      have_tla_traj = defined(tla_traj)
    else
      have_tla_traj = .false.
    end if
    if (present(perturbations)) then
      have_perturbs = associated(perturbations)
      if(have_perturbs) have_perturbs = size(perturbations) > 0
    else
      have_perturbs = .false.
    end if

    ! Is this the first simulation step, ever?
    first_step = call_counter == 1
    ! Is this the first time output is to be written/collected?
    first_output_step = (.not. have_done_output) .and. simrules%if_make_output

    pTranspMass => fu_concMM_ptr(cloud)
    pTranspXm   => fu_advection_moment_X_MM_ptr(cloud)
    pTranspYm   => fu_advection_moment_Y_MM_ptr(cloud)
    pTranspZm   => fu_advection_moment_Z_MM_ptr(cloud)

  !  call msg("Run_dispersion got dyn shopping list")
  !  call report(disp_dyn_shopping_list)
  !  call msg("Run_dispersion got ST shopping list")
  !  call report(disp_stat_shopping_list)
  !  call msg ("Oops..")

    !-----------------------------------------
    !
    ! Start grand loop over time.
    !

    command_string = 'Dispersion loop'
    call start_count(chCounterNm = command_string)

    fInjectedMassTotal(:) = 0.0

    num_steps = nint(fu_sec8(simrules%periodToCompute) / fu_sec8(simrules%timestep))

    !
    ! Write new progress file
    !
    if(simrules%ifwriteprogressfile)then
      call write_progress_file(simrules%chprogressfnmtemplate, 0)
      if(error)call unset_error(sub_name)
    endif

    loop_over_time: DO step_count = 1,num_steps
      now = simrules%startTime + simRules%timestep * (step_count - 1)



      call msg(' ')
      call msg('Time in loop now:' + fu_str(now) + ', run to: ' + &
             & fu_str(simrules%startTime + simrules%periodToCompute))
#ifdef DEBUG
    if(simRules%ifRunDispersion)then
      call msg("In run_dispersion1")
      call collect_total_masses(cloud)
      call report_total_masses(cloud, step_count, .false.)

      if (simrules%ifRunDispersion .and. fu_interval_positive(simrules%timestep) &
          & .and. fu_if_eulerian_present(simRules%dynamicsRules%simulation_type) &
          & .and. associated(fu_low_mass_threshold(simrules%chemicalRules))) then
        call check_masses(cloud,'beginning of step', fu_low_mass_threshold(simrules%chemicalRules))
        call check_mass_centres(cloud, 'beginning of step',0.9999)
      end if
    endif
#endif
       
      
      if(error)then
        call msg_warning('Problem at the beginning of time step',sub_name)
        exit loop_over_time
      endif


!call msg('before_output')
!if (fu_interval_positive(simrules%timestep)) then
!  call check_mass_centres(cloud, 'before output',0.9999)
!  call check_masses(cloud,'before output')
!endif
      
!call msg('Checking range before output')
!call check_supermarket_fields_range(meteoMarketPtr)
!call check_supermarket_fields_range(dispersionMarketPtr)
!call check_supermarket_fields_range(outputMarketPtr)

      !-----------------------------------------------------------------------
      !
      ! Collection of the data for output should be done every time step
      !
      command_string = 'Output'
      call start_count(chCounterNm = command_string)
      if (simrules%if_make_output) then 
        call collect_output( meteo_ptr,  &
                          & disp_buf_ptr, output_buf_ptr, cloud, now, &
                          & OutDef, wdr, simRules%timestep, simrules%dynamicsRules%simulation_type, &
                          & first_output_step, .false.)
        first_output_step = .false.
        have_done_output = .true.
      end if

      call stop_count(chCounterNm = command_string)

      !call da_debug_output()

      if(error) exit loop_over_time 

!call msg('Checking range after output')
!call check_supermarket_fields_range(meteoMarketPtr)
!call check_supermarket_fields_range(dispersionMarketPtr)
!call check_supermarket_fields_range(outputMarketPtr)

      ! Handle tangent linearization if needed
      ! 
      if (have_tla_traj) then
        if (fu_interval_positive(simrules%timestep)) then
          ind_step = step_count
        else
          ind_step = num_steps - step_count + 1
        end if
        call get_tla_step(tla_traj, tla_step, ind_step)
        call msg('Set TLA step, ind_step:', ind_step)
        if (error) return
      end if
      
      ! Handle perturbations, if any:
      if (have_perturbs) then
        call apply_perturbations(cloud, disp_buf_ptr, perturbations, simrules%darules%assim_interval, wdr, OutDef)
        if (error) return
        have_perturbs = .false.
      end if

      !
      ! Decide if we need new meteo or boundary data: check the available times in the supermarket
      ! and request more if needed
      !
      call supermarket_times(meteoMarketPtr, met_src_missing, valid_times, number_of_times)  !, analysis_time)
      if(error)return
      !
      ! Reftime refers to the middle of the timestep in the normal run but the whole run can be shifted
      ! if EnKF perturbs meteodata by shifting them
      !
      reftime = now + simrules%timestep / 2. + fu_meteo_time_shift(wdr)

      ! check if should read meteo
      ifGotNewMeteoData = (.not. fu_between_times(reftime , &
                               & fu_earliest_time(valid_times), &
                               & fu_latest_time(valid_times), accept_boundaries_too = .true.))

        call  prepare_meteo( wdr, simrules%iComputationAccuracy, simRules%diagnosticRules, &
                          & meteo_input_dyn_shopping_list, meteo_full_dyn_shopping_list, &
                          & disp_dyn_shopping_list, disp_stat_shopping_list, &
                          & fu_output_dyn_shopping_list(OutDef, now),&
                          & fu_output_st_shopping_list(OutDef, now), & 
                          & meteo_ptr, disp_buf_ptr, output_buf_ptr, &
                          & now, simrules%timestep, &
                          & meteoMarketPtr, dispersionMarketPtr, outputMarketPtr, ifGotNewMeteoData)

        
         if(error)EXIT loop_over_time
!      endif

!!
!! DEBUG only: a massive dump of the meteofield features
!!
!call msg('')
!call msg('Meteo buffer after setting pointers')
!call report(meteo_ptr,2)
!call msg('')

!call msg('Checking range after diagnostic fields')
!call check_supermarket_fields_range(meteoMarketPtr)
!call check_supermarket_fields_range(dispersionMarketPtr)
!call check_supermarket_fields_range(outputMarketPtr)

      if(simRules%ifRunDispersion)then
        !
        ! This is the right moment to (re)set low-mass threshold: current procedure may need meteodata
        ! Of course, this is all only at the first time step
        ! Essentially, we can have 4 cases: thresholds can be based on emission or model-measurement
        ! discrepancy, and they can be for forward and adjoint runs.
        ! To get the Finally, we might need to 
        !
        ! 4DVAR adjoint run emits discrepancies btw modelled and observed concentrations via injectAll
        ! which is not an emission injection. Mass threshold has to be different too.
        ! We shall do it via injecting all discrepancies into a just-nullified mass map, then collect 
        ! totals and call low_mass threshold setting sub.
        !
        ifCalculate = .not. fu_low_mass_thres_defined(simrules%chemicalRules, &
                                                    & ifforward=.not. simrules%if_invert_substeps)
        if(ifCalculate)then
          if(simrules%if_invert_substeps)then
            if(obs_ptrs%hasObservations)then
              !
              ! Adjoint from mdl-meas discrepancy
              !
              !              call reset_map_to_val(fu_concMM_ptr(cloud), 0.0)
              call reset_map_to_val(pTranspMass, 0.0)

              call injectAll(obs_ptrs, cloud, meteo_ptr, &            ! Inject all observations at once
                           & simRules%chemicalRules, simRules%dynamicsRules, &
                           ! the whole assimilation window, negative period length
                           & simrules%periodToCompute, now)  
              call reset_all(obs_ptrs, .false.)  ! reset the counters but not nullify model data
              call collect_total_masses(cloud, fInjectedMassTotal)
              call msg('===> Making the low-mass threshold for adjoint simulations. Initial RMSE:')
              call report_total_masses(cloud, step_count, .true.)

              iTmp= fu_nbr_of_species_transport(cloud)
              fInjectedMassTotal(1:iTmp) = abs(fInjectedMassTotal(1:iTmp))

              call reset_map_to_val(pTranspMass, 0.0)
              call reset_map_to_val(pTranspXm, 0.0)
              call reset_map_to_val(pTranspYm, 0.0)
              call reset_map_to_val(pTranspZm, 0.0)

!              if (debug_level > 0) then
                call msg('Total discrepances over the whole assimilation window.')
                call report_total_masses(cloud, step_count, .false.)
!              endif
            endif  ! adjoint inject observation discrepancy
            reftime = now
            periodTmp = simrules%periodToCompute
          else
            reftime = fu_earliest_start(cloud)
            periodTmp = fu_period_to_compute(wdr)
          endif   ! if invert timesteps

          call set_low_mass_threshold(simRules%chemicalRules, &
                                    & em_source, &
                                    & fu_species_transport(cloud), fu_species_emission(cloud), &
                                    & fu_nbr_of_species_transport(cloud), &
                                    & simRules%ResidenceInterval, &
                                    & periodTmp, &
                                    & simRules%timestep, &
                                    & simRules%iComputationAccuracy, &
                                    & reftime, &                            ! earliest_start
                                    & .not. simrules%if_invert_substeps, &  ! ifForwardRun
                                    & .not. (simrules%if_invert_substeps .and. &
                                           & obs_ptrs%hasObservations), &   ! ifFromEmission
                                    & fInjectedMassTotal, &                      ! discrepancy total
                                    ! 4 adjoint: instead of n emission 
                                    & fu_number_of_observed_cells(obs_ptrs), &  
                                    & fu_if_lagrangian_present(simRules%dynamicsRules%simulation_type))

          fInjectedMassTotal(:) = 0.0  !Reset it back Will be needed for accounting
        endif  ! ifCalculate threshold
        
        ! Point the thresholds either to forward or adjoint ones, calculated now or not.
        ! 
        call point_low_mass_thresholds(simrules%chemicalRules, ifforward=.not.simrules%if_invert_substeps)
        if(error)return

 !call msg('===> dispersion_supplementary reporting cnc thresholds:', simrules%chemicalRules%low_cnc_trsh)
 !call msg('===> dispersion_supplementary reporting mass thresholds:', simrules%chemicalRules%low_mass_trsh)
        !
        ! Read boundary conditions if needed
        !
        ! They either not needed by definition or have to be excluded in adjoint runs of 4DVAR
        ! Note that just-adjoint might still have boundaries
        !
        if(fu_ifReadBoundaries(simRules%IniBoundaryRules))then
          
          call supermarket_times(BCMarketPtr, met_src_missing, valid_times, number_of_times)  !, analysis_time)
          if(error)return
          !
          ! Reftime refers to the middle of the timestep in the normal run but the whole run can be shifted
          ! if EnKF perturbs meteodata by shifting them
          !
          reftime = now + simrules%timestep / 2. + fu_meteo_time_shift(wdr)
          if(now == simrules%startTime .or. &    ! start of the run, read no matter what
           & (.not. fu_between_times(reftime , &
                                   & fu_earliest_time(valid_times), &
                                   & fu_latest_time(valid_times), accept_boundaries_too = .true.)))then
            command_string = 'Boundary_condition_acquisition'
            call msg(command_string)
            call start_count(chCounterNm = command_string)
            do iTmp = 1, fu_nr_boundary_inFiles(simRules%IniBoundaryRules)
              pShpLst => fu_shplst(simRules%IniBoundaryRules, iTmp)
              CALL fix_shopping_time_boundaries(pShpLst, now, (now + simRules%timestep))
!call msg('')
!call msg('BC Shopping list:')
!call report(fu_shplst(simRules%IniBoundaryRules, iTmp))
!call msg('')
            enddo
            CALL fill_boundary_market(BCMarketPtr, simRules%IniBoundaryRules, first_step)
            IF (error) EXIT loop_over_time
            call stop_count(chCounterNm = command_string)
          
          
!          bc_lim1 = fu_closest_obstime(reftime, backwards, fu_obstime_interval(simRules%IniBoundaryRules))
!          bc_lim2 = fu_closest_obstime(reftime, forwards, fu_obstime_interval(simRules%IniBoundaryRules))
!          IF(now == simrules%startTime)THEN 
!            !--- For the first time data must always be acquired
!            bc_lim1_prev = bc_lim1 + simRules%timestep
!            bc_lim2_prev = bc_lim2 + simRules%timestep
!          endif
!          IF (error) EXIT loop_over_time
!          IF(.not.(fu_between_times(bc_lim1, bc_lim1_prev, bc_lim2_prev, .true.) &
!           & .and. fu_between_times(bc_lim2, bc_lim1_prev, bc_lim2_prev, .true.)))THEN
!
!            command_string = 'Boundary_condition_acquisition'
!            call msg(command_string)
!            call start_count(chCounterNm = command_string)
!            !
!            ! ATTENTION. There used to be now, now+timestep*2. in both time boundaries
!            ! However, I have not found any reason to keep it.
!            !
!            do iTmp = 1, fu_nr_boundary_inFiles(simRules%IniBoundaryRules)
!              CALL fix_shopping_time_boundaries(fu_shplst(simRules%IniBoundaryRules, iTmp), &
!                                              & now, (now + simRules%timestep))
!!call msg('')
!!call msg('BC Shopping list:')
!!call report(fu_shplst(simRules%IniBoundaryRules, iTmp))
!!call msg('')
!            enddo
!            CALL fill_boundary_market(BCMarketPtr, simRules%IniBoundaryRules, first_step)
!            IF (error) EXIT loop_over_time
!
!            !command_string = 'Boundary_condition_acquisition'
!            call stop_count(chCounterNm = command_string)

            command_string = 'Boundary_condition_processing'
            call msg(command_string)
            call start_count(chCounterNm = command_string)

            CALL fillBoundaryStruct(simRules%IniBoundaryRules, BCMarketPtr, &
                                  & now, fu_boundaryStructures(cloud), meteo_ptr, first_step)
            IF (error) EXIT loop_over_time

            !command_string = 'Boundary_condition_processing'
            call stop_count(chCounterNm = command_string)

          endif  ! if time to read boundary conditions

!          bc_lim2_prev = bc_lim2
!          bc_lim1_prev = bc_lim1

        endif  ! if boundaries are needed at all
        ! kludge pt. 2
        call set_time_direction_sm(meteoMarketPtr, fu_interval_positive(simRules%timestep))
        call set_time_direction_sm(dispersionMarketPtr, fu_interval_positive(simRules%timestep))
        call set_time_direction_sm(BCMarketPtr, fu_interval_positive(simRules%timestep))

        !call set_disp_buffer_pointers(now, disp_buf_ptr)
        command_string = 'Transformation'
        call start_count(chCounterNm = command_string)

        !call msg('after_boundaries')
        if (simrules%ifRunDispersion .and. fu_interval_positive(simrules%timestep) &
          & .and. fu_if_eulerian_present(simRules%dynamicsRules%simulation_type)) then
!PAR-OLE tweaks are now applied: no need for the following
!          if(any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm42_strato) .or. &
!           & any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm42_strato_SOA) .or. &
!           & any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm4_SOA) .or. &
!           & any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm5_strato_SOA) .or. &
!           & any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm5_SOA) .or. &
!           & any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm4)) &
!           &   call check_species(cloud,'Before transformation', fu_low_mass_threshold(simrules%chemicalRules), &
!                                & 'OLE_gas', 'PAR_gas')
        end if

        !        call msg_test('Prepare advection-tranformation')
        !
        ! Prepare to the advection-diffusion-transformation-deposition computations
        !
        call prepare_deposition(meteo_ptr, disp_buf_ptr, simRules%chemicalRules%rulesDeposition)
        if(error)return

        call pick_meteo_pointers(pMeteo_input, meteo_ptr, disp_buf_ptr)
        if(error) exit loop_over_time

        call stop_count(chCounterNm = command_string)
        if(error) EXIT loop_over_time

      else
        ! No dispersion run
        ! kludge pt. 2
        call set_time_direction_sm(meteoMarketPtr, fu_interval_positive(simRules%timestep))
        call set_time_direction_sm(dispersionMarketPtr, fu_interval_positive(simRules%timestep))
        call set_time_direction_sm(BCMarketPtr, fu_interval_positive(simRules%timestep))
        
      end if ! ifRunDispersion


      if(simRules%ifRunDispersion)then
        !-----------------------------------------------------------------------
        !
        ! Before making the dump output and starting dispersion computations,
        ! compute the now-time emission. For many sources and species this is void
        ! but for some others it is not.
        !
        
        
        if (debug_level > 0 .and. fu_interval_positive(simrules%timestep)) then
          call check_masses(cloud,'before dynamic emission', &
                          & fu_low_mass_threshold(simRules%chemicalRules))
          call smpi_advection_barrier()
        end if

        IF (error) EXIT loop_over_time
        
        !-----------------------------------------------------------------------------
        !
        ! Before running into the main time step, we have to convert a few variables
        ! to their moments: coordinates of centres of masses and aerosol particle sizes.
        ! So far, there is only aerosol stuff
        !
        if(associated(simRules%chemicalRules%ChemRunSetup%mapVolume2NumberCnc))then
          if(fu_if_eulerian_present(simRules%dynamicsRules%simulation_type))then
            call flip_moment_and_basic_var_Euler( &
                                 & simRules%chemicalRules%ChemRunSetup%mapVolume2NumberCnc, &
                                 & fu_aerosolMM_ptr(cloud), &
                                 & fu_concMM_ptr(cloud), &
                                 & to_moment)
          elseif(fu_if_lagrangian_present(simRules%dynamicsRules%simulation_type))then
            call flip_moment_and_basic_var_Lagr( &
                                 & simRules%chemicalRules%ChemRunSetup%mapVolume2NumberCnc, &
                                 & fu_aerosolLP_ptr(cloud), &
                                 & fu_concLP_ptr(cloud), &
                                 & fu_nbr_of_lagr_particles(cloud, .false.), &
                                 & fu_nbr_of_species_transport(cloud), &
                                 & to_moment)
          endif
        endif
        if(error)return 

        fInjectedMass(:) = 0.0

        if (.not. simrules%if_invert_substeps) then

#ifdef DEBUG           
          if (simrules%ifRunDispersion .and. fu_interval_positive(simrules%timestep) &
            & .and. fu_if_eulerian_present(simRules%dynamicsRules%simulation_type)) then
            call check_masses(cloud,'before emission', fu_low_mass_threshold(simrules%chemicalRules))
            call check_mass_centres(cloud, 'before emission',0.9999)
!PAR-OLE tweaks are now applied: no need for the following
!            if(any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm42_strato) .or. &
!             & any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm42_strato_SOA) .or. &
!             & any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm4_SOA) .or. &
!             & any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm5_strato_SOA) .or. &
!             & any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm5_SOA) .or. &
!             & any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm4)) &
!             &   call check_species(cloud,'before emission', fu_low_mass_threshold(simrules%chemicalRules), &
!                                  & 'OLE_gas', 'PAR_gas')
          end if
#endif
          if (smpi_global_tasks > 1) then
            command_string = 'Before-emission barrier'
            call start_count(chCounterNm = command_string)
            call smpi_advection_barrier()
            call stop_count(chCounterNm = command_string)
          endif

          command_string = 'now_time_emission'
          call start_count(chCounterNm = command_string)

          call msg_test('Computing now-time emission (dispersion supplementary)')
          call set_dynamic_emission(cloud, &
                                  & em_source, &
                                  & meteo_ptr, &
                                  & disp_buf_ptr, &
                                  & simRules%chemicalRules, simRules%dynamicsRules, &
                                  & now, &
                                  & simRules%timestep, simRules%residenceInterval, &
                                  & fInjectedMass, &
                                  & dispersionMarketPtr)
                                  
        call stop_count(chCounterNm = command_string)
      !Another estimate of load imbalance
          if (smpi_global_tasks > 1) then
            command_string = 'After-emission barrier'
            call start_count(chCounterNm = command_string)
            call smpi_advection_barrier()
            call stop_count(chCounterNm = command_string)
          endif
!        call msg('Reporting emission mass after set_dynamic_emission')
!        call report_emission_mass()
#ifdef DEBUG           
          if (simrules%ifRunDispersion .and. fu_interval_positive(simrules%timestep) &
            & .and. fu_if_eulerian_present(simRules%dynamicsRules%simulation_type)) then
            call check_masses(cloud,'after emission', fu_low_mass_threshold(simrules%chemicalRules))
            call check_mass_centres(cloud, 'after emission',0.9999)
!PAR-OLE tweaks are now applied: no need for the following
!            if(any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm42_strato) .or. &
!             & any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm42_strato_SOA) .or. &
!             & any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm4_SOA) .or. &
!             & any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm5_strato_SOA) .or. &
!             & any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm5_SOA) .or. &
!             & any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm4)) then
!              call check_species(cloud,'after emission', fu_low_mass_threshold(simrules%chemicalRules), &
!                               & 'OLE_gas', 'PAR_gas')
!            end if
          end if
#endif
!         call msg('Reporting emission mass after checks')
!        call report_emission_mass()

          !   call msg_warning('Skip emission')
          if(fu_if_dump_emission_flux(em_source) .and. simrules%if_make_output) then
            !
            ! Each time the emission mass map is filled-in, we check whether it should be dumped
            !
            call dump_emission_mass_map(em_source, now, &
                                      & simRules%chCaseNm, fu_output_template(OutDef), &
                                      & fu_emisMM_ptr(cloud), &
                                      ! & null_pointer, null_pointer, null_pointer
                                      & fu_emission_moment_X_MM_ptr(cloud), &
                                      & fu_emission_moment_Y_MM_ptr(cloud), &
                                      & fu_emission_moment_Z_MM_ptr(cloud), &
                                      & simrules%timestep)
          end if

!          call msg('Reporting emission mass after dump')
          call report_emission_mass()
          
        endif


        IF (error) EXIT loop_over_time

        !
        !  Trajectories positions must be saved after emission to get first point at the source
        !
        traj_set => fu_traj_set(outDef)
        if (defined(traj_set)) &
           & call store_traj_time(traj_set, meteo_ptr, fu_lpset_ptr(cloud), now)

        
        ! In adjoint runs for DA, we inject the discrepancy between measurement and model.
        
        if (obs_ptrs%hasObservations .and. simrules%if_invert_substeps) then
          command_string = 'da_observations'
          call start_count(chCounterNm = command_string)
          call injectAll(obs_ptrs, cloud, meteo_ptr,&
                       & simRules%chemicalRules, simRules%dynamicsRules, simRules%timestep, now)
          if (debug_level > 0) then
            call report_total_masses(cloud, step_count, .false.)
          end if
          !            call msg('low mass thr', fu_low_mass_threshold(simRules%chemicalRules))
          call stop_count(chCounterNm = command_string)
        end if  ! DA obs & invert_step
        
        IF (error) EXIT loop_over_time
  
        !-----------------------------------------------------------------------
        !
        ! Calculate advection and random walk for one timestep
        ! Now it is in advections and random_walk modules
        !
        !call msg_test('SKIP advection****************')
        if (simrules%if_invert_substeps) then
          call transform_pollution_cloud_v5(cloud, now, simRules%timestep, &
                                          & tla_step, &
                                          & meteo_ptr, disp_buf_ptr, &
                                          & pMeteo_input, & 
                                          & wdr, &
                                          & simRules%chemicalRules, simRules%dynamicsRules, &
                                          & fu_ifDryDep_cumulative(OutDef), fu_ifWetDep_cumulative(OutDef))
          !call flip_moment_and_basic_var(cm_to_moment, fu_advection_moment_X_MM_ptr(cloud), &
          !                             & fu_concMM_ptr(cloud), to_moment)
          !call flip_moment_and_basic_var(cm_to_moment, fu_advection_moment_Y_MM_ptr(cloud), &
          !                             & fu_concMM_ptr(cloud), to_moment)
          !call flip_moment_and_basic_var(cm_to_moment, fu_advection_moment_Z_MM_ptr(cloud), &
          !                             & fu_concMM_ptr(cloud), to_moment)
          if (error) return
          if (debug_level > 0) then
            call msg('')
            call msg('After adjoint chemistry:')
            call collect_total_masses(cloud)
            call report_total_masses(cloud, step_count, .true.) 
            call msg('')
          end if
        end if    ! invert steps

#ifdef DEBUG
        if (simrules%ifRunDispersion .and. fu_interval_positive(simrules%timestep) &
          & .and. fu_if_eulerian_present(simRules%dynamicsRules%simulation_type)) then
          call check_masses(cloud,'before advection', fu_low_mass_threshold(simrules%chemicalRules))
          call check_mass_centres(cloud, 'before advection',0.9999)
!PAR-OLE tweaks are now applied: no need for the following
!          if(any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm42_strato) .or. &
!           & any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm42_strato_SOA) .or. &
!           & any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm4_SOA) .or. &
!           & any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm5_strato_SOA) .or. &
!           & any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm5_SOA) .or. &
!           & any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm4)) &
!           &   call check_species(cloud,'before advection', fu_low_mass_threshold(simrules%chemicalRules), &
!                                & 'OLE_gas', 'PAR_gas')
        end if
#endif

        command_string = 'Advection'
        call start_count(chCounterNm = command_string)

        CALL advect_pollution_cloud_v4(cloud, &
                                     & now, &
                                     & simRules%timestep,&
                                     & simrules%dynamicsRules%diffusionMethod, &
                                     & simrules%if_invert_substeps, & ! have_centre_of_mass
                                     & simrules%IniBoundaryRules, &
                                     & meteo_ptr, & 
                                     & disp_buf_ptr, &
                                     & wdr, &
                                     & simRules%chemicalRules, simRules%dynamicsRules)

        if (debug_level > 0) then
          call msg('')
          call msg('After advection:')
          call collect_total_masses(cloud)
          call report_total_masses(cloud, step_count, .true.) 
          call msg('')
        end if

#ifdef DEBUG
        if (simrules%ifRunDispersion .and. fu_interval_positive(simrules%timestep) &
          & .and. fu_if_eulerian_present(simRules%dynamicsRules%simulation_type)) then
          call check_masses(cloud,'after advection', fu_low_mass_threshold(simrules%chemicalRules))
          call check_mass_centres(cloud, 'after advection',0.9999)
!PAR-OLE tweaks are now applied: no need for the following
!          if(any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm42_strato) .or. &
!           & any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm42_strato_SOA) .or. &
!           & any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm4_SOA) .or. &
!           & any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm5_strato_SOA) .or. &
!           & any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm5_SOA) .or. &
!           & any(simrules%chemicalRules%iTransformTypes(:) == transformation_cbm4)) &
!           &   call check_species(cloud,'after advection', fu_low_mass_threshold(simrules%chemicalRules), &
!                                & 'OLE_gas', 'PAR_gas')
        end if
#endif

        IF (error) EXIT loop_over_time

        call stop_count(chCounterNm = command_string)

        !
        ! Now we can make a radioactive decay or chemical transformations - 
        ! depending on the type of the cocktails.
        !
          if (smpi_global_tasks > 1) then
            command_string = 'Before-transformation barrier'
            call start_count(chCounterNm = command_string)
            call smpi_advection_barrier()
            call stop_count(chCounterNm = command_string)
          endif
        call msg_test('Transformation')

        command_string = 'Transformation'
        call start_count(chCounterNm = command_string)

        if (.not. simrules%if_invert_substeps) then
          call transform_pollution_cloud_v5(cloud, now, simRules%timestep, tla_step, &
                                          & meteo_ptr, disp_buf_ptr, &
                                          & pMeteo_input, & 
                                          & wdr, simRules%chemicalRules, simRules%dynamicsRules, &
                                          & fu_ifDryDep_cumulative(OutDef), fu_ifWetDep_cumulative(OutDef))
        end if


        call stop_count(chCounterNm = command_string)

          if (smpi_global_tasks > 1) then
            command_string = 'After-transformation barrier'
            call start_count(chCounterNm = command_string)
            call smpi_advection_barrier()
            call stop_count(chCounterNm = command_string)
          endif

        if(error) exit loop_over_time
        if (obs_ptrs%hasObservations .and. .not. simrules%if_invert_substeps) then
          command_string = 'da_observations'
          call start_count(chCounterNm = command_string)
          !
          ! Here is terrible violation of all rules.
          ! eruption height is simultaneously observable quantity and non-cloud emission
          ! parameter, which can be assimilated directly. A formal solution is here, to be
          ! removed asap.
          !
          call observeAll(obs_ptrs, cloud, meteo_ptr, disp_buf_ptr, simRules%chemicalRules, &
                        & simRules%timestep, now, fu_eruption_height(perturbations))

          !call observeAll(obs_ptrs, cloud, meteo_ptr, &
          !              & simRules%chemicalRules, simRules%timestep, now )
          if (debug_level > 0) then
            call collect_total_masses(cloud)
            call report_total_masses(cloud, step_count, .false.)
            !            call msg('low mass thr', fu_low_mass_threshold(simRules%chemicalRules))
          end if
          call stop_count(chCounterNm = command_string)
        end if
        
        if (simrules%if_invert_substeps) then
          command_string = 'now_time_emission'
          call start_count(chCounterNm = command_string)
          call msg_test('Computing now-time emission, inverted substeps (dispersion supplementary)')
          call set_dynamic_emission(cloud, &
                                  & em_source, &
                                  & meteo_ptr, &
                                  & disp_buf_ptr, &
                                  & simRules%chemicalRules, simRules%dynamicsRules, &
                                  & now, &
                                  & simRules%timestep, simRules%residenceInterval, &
                                  & fInjectedMass, &
                                  & dispersionMarketPtr)
          !   call msg_warning('Skip emission')
          call stop_count(chCounterNm = command_string)
          call report_emission_mass()

        end if

      endif ! if run dispersion

      !
      ! Store passed time and set new time of loop. If the file is not needed,
      ! the unit will be int_missing
      !
      if(info_funit /= int_missing)then
        write(info_funit,*) step_count,'  ', trim(fu_str(now))
        call flush(info_funit)
      endif

      first_step = .false.
    
!      !Some estimate of load imbalance
!     if (smpi_global_tasks > 1) then
!        command_string = 'End-loop barrier'
!        call start_count(chCounterNm = command_string)
!        call smpi_advection_barrier()
!        call stop_count(chCounterNm = command_string)
!      endif


      if(simrules%ifwriteprogressfile)then
        call write_progress_file(simrules%chprogressfnmtemplate, step_count*100/num_steps)
      endif
      
    END DO loop_over_time !************* TIME LOOP
    command_string = 'Dispersion loop'
    call stop_count(chCounterNm = command_string)

  contains

    !===============================================================================

    subroutine da_debug_output()
      implicit none
      character(len=fnlen) :: filename
      integer :: isrc

      if (fu_interval_positive(simrules%timestep)) then
        write(filename, '(A,A,A,I0,A,A,A)') trim(simrules%darules%outputDir), dir_slash, &
             & 'fwd', call_counter, '_', trim(fu_str(now,.false.)), '.grads'
        call mass_map_to_grads_file(fu_concMM_ptr(cloud), 1, filename, now, 1.0)
      else
        do isrc = 1, fu_nbr_of_sources(cloud)
          write(filename, '(A,A,A,I0,A,I0,A,A,A)') trim(simrules%darules%outputDir), dir_slash, &
               & 'adj_', isrc, '_', call_counter, '_', &
               & trim(fu_str(now,.false.)), '.grads'
          call mass_map_to_grads_file(fu_concMM_ptr(cloud), isrc, filename, now, 1.0)
        end do
      end if
    end subroutine da_debug_output

    !===============================================================
    
    subroutine report_emission_mass()
      implicit none
      integer :: iTmp, jTmp
      real, dimension(:), pointer :: work, mass_src, mass_cloud, mass_cum
      type(silam_species), dimension(:), pointer :: pSpecies
      logical :: ifOk

      iTmp = fu_nbr_of_species_emission(cloud)
      work => fu_work_array(iTmp*6)
      mass_src => work(1:iTmp) 
      mass_cloud     =>  work(  iTmp+1:2*iTmp)
      mass_cum       =>  work(2*iTmp+1:3*iTmp)

      fInjectedMassTotal(1:iTmp) = fInjectedMassTotal(1:iTmp) + fInjectedMass(1:iTmp)
      mass_src(1:iTmp) = fInjectedMass(1:iTmp) !Type conversion
      mass_cum(1:iTmp) = fInjectedMassTotal(1:iTmp) !Type conversion



      if (smpi_global_tasks > 1) then 
        call emis_mass_report('Total subdomain injected emission mass as counted from sources', &
                    & fu_species_emission(cloud), mass_src, '__src_my '//fu_str(now))

        if (debug_level > 0) then
          call get_cloud_inventory(cloud, .true., species_emission, pSpecies, iTmp, fInjectedMass)
          mass_cloud(1:iTmp) = fInjectedMass(1:iTmp)
          call emis_mass_report('Total subdomain injected emission mass as counted from map', &
                   &  pSpecies, mass_cloud, '__mp_my '//fu_str(now) )

          call emis_mass_report('Cumulative subdomain injected emission mass', &
               & fu_species_emission(cloud), mass_cum, '__cum_my '//fu_str(now))
        end if

        !Exchange and reset pointers to WHOLE_MPI
        call  smpi_reduce_add(work(1:3*iTmp), work(3*iTmp+1:6*iTmp), 0, smpi_adv_comm, ifOk)
        if (fu_fails(ifOk, "smpi_reduce_add failed", "report_emission_mass")) return
        mass_src   =>  work(3*iTmp+1:4*iTmp)
        mass_cloud =>  work(4*iTmp+1:5*iTmp)
        mass_cum   =>  work(5*iTmp+1:6*iTmp)

      endif

      ! Master reports for all
      if (smpi_global_rank == 0) then
        call msg("")
        call emis_mass_report('Total injected emission mass as counted from sources', &
                      & fu_species_emission(cloud), mass_src, '__src '//fu_str(now))

        if (debug_level > 0) then
          call get_cloud_inventory(cloud, .true., species_emission, pSpecies, iTmp, fInjectedMass)
          mass_cloud(1:iTmp) = fInjectedMass(1:iTmp)
          call emis_mass_report('Total emission mass as counted from map', &
                     &  pSpecies, mass_cloud, '__mp '//fu_str(now) )

          call emis_mass_report('Cumulative emission mass', &
                 & fu_species_emission(cloud), mass_cum, '__cum '//fu_str(now))
        end if
      endif
      call free_work_array(work)

    end subroutine report_emission_mass

    !=========================================================================
    
    subroutine emis_mass_report(head, species, amounts, timestr)
        character(len=*), intent(in) :: head,timestr
        type(silam_species), dimension(:), intent(in) :: species
        real, dimension(:), intent(in) :: amounts

        integer :: jTmp
        if(.not. all(amounts ==  0.0)) then
          call msg(head)
          do jTmp = 1, size(species)
            call msg(trim(fu_str(species(jTmp)))//timestr, amounts(jTmp))
          enddo
        endif
    end subroutine emis_mass_report
    
  end subroutine run_dispersion
  
  integer function fu_nbr_of_timesteps(simrules) result(nbr)

    !!! The function actually returns one too much steps
    implicit none
    type(general_dispersion_rules), intent(in) :: simrules
    nbr = abs(simrules%periodToCompute / simrules%timestep) + 1
  end function fu_nbr_of_timesteps
  
end module dispersion_supplementary
