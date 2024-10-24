MODULE dispersion_models

  !
  ! This is the starting module of SILAM.
  ! The actual work is done by in dispersion_supplementary,
  ! here the control file is consumed and first steps towards the run 
  ! initialization are performed.
  !
  ! Author: Mikhail Sofiev, mikhail.sofiev@fmi.fi
  !
  ! All units: SI
  ! 
  ! All degrees: real numbers in degrees (degrees and parts per
  ! hundred/parts per thousand , NOT in degrees and minutes
  ! 
  ! Language: ANSI standard Fortran 90
  ! 
  use da_driver
  use ensemble_driver
!  use silam_partitioning

  
  
!use source_terms_time_params
  

  IMPLICIT NONE

  private

  ! The public functions and subroutines available in this module:
  PUBLIC start_silam_v5

  ! The private functions and subroutines not to be used elsewhere:
  PRIVATE silam_v5
  PRIVATE set_standard_setup

!  private create_cocktail_map_set
!  private distribute_cocktail_maps

CONTAINS

  ! ***************************************************************

  SUBROUTINE start_silam_v5 (control_param_fname, had_error)

    ! Description:
    ! Read the setup, control parameters, and source term files,
    ! creates the correct input parameters and then
    ! runs the SILAM dispersion model.
    !
    ! A split:
    ! Setup file determines internal parameters of the model: advection type,
    ! rw-type, ABL method, etc. Ordinary user should not know about this file.
    ! Control file describes the specific run parameters - grid, period, etc.
    ! Source term file describes the source of particles.
    !
    ! Input:
    ! path to the control file.
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    CHARACTER(LEN=*),INTENT(in):: control_param_fname
    logical, intent(out) :: had_error

    ! Local declarations:
    TYPE(silam_pollution_cloud),target :: cloud
    TYPE(silam_pollution_cloud),pointer :: cloudPtr => null()
    type(meteo_data_source) :: met_src
    TYPE(silja_wdr),target :: wdr
    TYPE(silja_wdr), pointer :: wdrPtr => null()
    type(silam_output_definition), target :: outDef
    type(silam_output_definition),pointer :: outDefPtr => null()
    TYPE(silam_source),target :: em_source ! emission source
    TYPE(silam_source),pointer :: em_sourcePtr => null() ! emission source
    type(observationPointers) :: daPointers
    type(general_dispersion_rules) :: simRules
    INTEGER :: i, control_file, status, n_levs, time_sign
    CHARACTER (LEN=clen) :: line, setup_FNm, chOutTimeFileSeries, chOutType
    CHARACTER(len=fnlen) :: fname_string, fname_oro_string, chOutVarLst, chOutTemplate, &
                          & chGribCodeTableFNm, source_term_fname
    LOGICAL :: eof
    type(Tsilam_namelist), pointer :: nlStandardSetup =>null()
    type(Tsilam_namelist_group), pointer :: nlSetupGrp => null()

    outDef = output_definition_missing

    had_error = .true. ! if all goes ok, will be set false in the end of time loop
    !
    ! Time counter...
    !
    line = 'Overall_initialisation'
    call start_count(chCounterNm = line)

    ! The main buffer of temporary fields
    !
    CALL initialize_work_arrays()

    !call exec_tests()

    !----------------------------------------
    !
    ! Open the control parameter file
    ! 
    control_file = fu_next_free_unit()

    OPEN(file = control_param_fname,&
        & unit = control_file,&
        & action = 'read',&
        & status = 'old',&
        & iostat = status)

    IF (status /= 0) THEN
      CALL set_error('cannot open input file:' + control_param_fname,'start_silam_v5')
      RETURN
    END IF

    !-------------------------------------------------------------
    !
    ! It must start from the CONTROL_V... line
    !
    CALL next_line_from_input_file(control_file, line, eof)
    call str_u_case(line)
    IF (error.or.eof) RETURN
    if(index(line,'CONTROL_V5_3') == 0)then
      call set_error('File must begin from  CONTROL_V5_3, not:' + line, 'start_silam_v5')
      return
    endif
    
    em_sourcePtr => em_source
    cloudPtr => cloud
    OutDefPtr => OutDef

    !
    ! There can be several formats of the control file. Each
    ! has to be treated appropriately
    !
    if(line(1:len_trim(line)) == 'CONTROL_V5_3')then
      !
      ! Setup v5_3
      !
      call setup_v5_3()   !------------------------------------------ setup v5.2
      if(error)return

    else
      call msg("Failed to interpret the line")
      call msg(">"+line+"<", len_trim(line))
      call set_error('Only CONTROL_V5_3 is allowed','start_silam_v5')
      return
    endif

    !
    ! When the basic model parts are intialised, start the run.
    ! Here the select the type of the model: Random-walk-Lagrangian or SCD-Eulerian.
    ! Later more options may come.
    !
    CALL silam_v5(simRules, &        ! A set of all simulation rules controlling the run
                & em_sourcePtr,&     ! emission source
!                & srcMarketPtr, &    ! emission mini-market of fields
                & cloudPtr, &        ! pollution cloud - whether Lagrangian or Eulerian
                & wdrPtr, &          ! weather data rules
                & daPointers, &
                & outDefPtr, &       ! output definitions
                & nlSetupGrp, &        ! A group of setup namelists
                & had_error) ! if an error occurred


  CONTAINS


    !=====================================================================================
    !=====================================================================================

    subroutine setup_v5_3()
      !
      ! The main advantage of this type of initialisation is that it does not do
      ! anything by itself but rather relies on appropriate procedures in each 
      ! class and silam_namelist class.
      ! The method is: every part of the model should be able to initialise 
      ! itself from the namelist. 
      !
      implicit none

      ! Local variables
      type(Tsilam_namelist), pointer :: nlPtr, nlMultiRun
      integer :: iSrcIdType
      logical :: ifEmissionGiven
      integer, dimension(:), pointer :: iTypes
      type(Toptical_density_rules), pointer :: pOpticalRules

      had_error = .true.
      !
      ! Create the namelist by reading the control file till the end
      ! 
      nlSetupGrp => fu_read_namelist_group(control_file, .true., 'END_CONTROL_V5_3')
      close(control_file)
      IF (error) RETURN

      call msg('')
      call msg('The following control namelist group v5_3 has been read:')
      call msg("")
      call msg("CONTROL_V5_3")
      call write_namelist_group(run_log_funit, nlSetupGrp, .true.)
      call msg("END_CONTROL_V5_3")
      call msg('')

      call set_missing(simRules%chemicalRules)

      !-------------------------------------------------------
      !
      ! Check the MPI related parameters for MPI runs, sets missing values for non-MPI.
      call smpi_setup_parameters(nlSetupGrp)
      if (error) return

      ! 
      ! Determine the general run parameters
      !
      nlPtr => fu_namelist(nlSetupGrp,'general_parameters')
      if(.not. associated(nlPtr))then
        call set_error('Cannot find the general_parameters namelist','setup_v5_3')
        return
      endif
      !
      ! Case name
      !
      simRules%chCaseNm = fu_expand_environment(fu_content(nlPtr,'case_name'))
      !
      ! Type of the run: Eulerian / Lagrangian / hybrid. Settings are strongly different, so no
      ! default is allowed: operator must understand its request.
      !
      select case(fu_str_u_case(fu_content(nlPtr,'simulation_type')))
        case('EULERIAN')
          simRules%dynamicsRules%simulation_type = eulerian_flag
        CASE ('LAGRANGIAN')
          simRules%dynamicsRules%simulation_type = lagrangian_flag
        CASE ('HYBRID')
          simRules%dynamicsRules%simulation_type = hybrid_flag
        case default
          call set_error('simulation_type must be LAGRANGIAN / EULERIAN / HYBRID, not:' + &
                         & fu_content(nlPtr,'simulation_type'),'setup_v5_3')
          return
      end select
      call msg("Simulation type", simRules%dynamicsRules%simulation_type)

      if(fu_if_lagrangian_present(simRules%dynamicsRules%simulation_type))then
        select case(fu_str_u_case(fu_content(nlPtr,'lagrangian_particles_removal')))
          case('ABSOLUTE_AGE')
            simRules%dynamicsRules%projectionRules%criterionType = maximum_age_flag
            simRules%dynamicsRules%projectionRules%LP_age_max = &
                                       & fu_sec(fu_set_named_interval(fu_content(nlPtr,'max_age')))
          case('RESIDENCE_INTERVAL_FRACTION')
            simRules%dynamicsRules%projectionRules%criterionType = residence_interval_fraction
            simRules%dynamicsRules%projectionRules%LP_age_max = &
                                       & fu_content_real(nlPtr,'residence_interval_fraction')
          case('RELATIVE_SIZE_INCREMENT')
            simRules%dynamicsRules%projectionRules%criterionType = maximum_relative_size_flag
            simRules%dynamicsRules%projectionRules%LP_relative_size_max = &
                                       & fu_content_real(nlPtr,'max_relative_size')
          case('NONE')
            simRules%dynamicsRules%projectionRules%criterionType = never_project_flag
          case default
            call set_error('Unknown lagrangian_particles_removal' + &
                         & fu_content(nlPtr,'lagrangian_particles_removal'),'setup_v5_3')
            return
        end select
      endif
      !
      ! Computational accuracy
      !
      simRules%iComputationAccuracy = fu_content_int(nlPtr,'computation_accuracy')
      if(error .or. simRules%iComputationAccuracy == int_missing)then
        call msg_warning('Wrong computaiton_accuracy. Must be from 0 to 10. Default 5 is taken')
        if(error)call unset_error('setup_v5_3')
        simRules%iComputationAccuracy = 5
      endif
      !
      !------------------------------------------------------------------------ STANDARD SETUP
      !
      ! Read the standard setup file: dynamic rules are known
      !
      nlStandardSetup => fu_namelist(nlSetupGrp,'STANDARD_SETUP')
      if(fu_fails(associated(nlStandardSetup),'Absent STANDARD_SETUP namelist','setup_v5_3'))return
      !
      ! Expanding the path is somewhat tricky: there may be a setup directory. Otherwise ^ means
      ! the control file directory
      !
      if(len_trim(fu_content(nlStandardSetup,'standard_setup_directory')) > 0)then
        call expand_paths(nlStandardSetup, fu_content(nlStandardSetup,'standard_setup_directory') + &
                                                    & dir_slash + 'tmp.tmp')
      else
        call expand_paths(nlStandardSetup, control_param_fname)
      endif
      if (error) return
      call set_standard_setup(nlStandardSetup, simRules%dynamicsRules, .false.)
      IF (error) RETURN
      !
      ! Initalise input-output structures
      !
      nlPtr => fu_namelist(nlSetupGrp,'output_parameters')
      call init_io_files(nlStandardSetup, nlPtr, max_nbr_of_grads_files)
!      !
!      ! Chemicals are needed for source creation, since many
!      ! substance features are defined through them. 
!      !
!      call init_chemical_materials(nlStandardSetup)
!      if(error)return
      !
      ! Upon initialisation of chemical materials, store chemical transformation rules
      !
      nlPtr => fu_namelist(nlSetupGrp, 'transformation_parameters')
      if(fu_fails(associated(nlPtr),'Absent transformation_parameters namelist','setup_v5_3'))return

      call set_chemistry_rules(nlPtr, nlStandardSetup, simRules%chemicalRules, .false., &
                      & fu_if_lagrangian_present(simRules%dynamicsRules%simulation_type))
      if(error)return
      !
      ! Initalise timezones
      !
      call init_timezones(nlStandardSetup) 
      if(error)return
      !
      ! Decide whether & how the sources are to be split.
      ! The idea is: the split can be absent, by-name, by-sector or by-name-and-sector
      ! I made it here, which looks clumsy. Indeed, this is an output parameter but
      ! it is stored in the source structure and, consequently, I have to read it at 
      ! wrong time moment. But this gives me a complete encapsulation of the source id
      ! principle. All other modules simply call for the source name and even do not know
      ! that this source id is, while it can consolidate tens of individual sources.
      !

      nlPtr => fu_namelist(nlSetupGrp,'output_parameters')
      !call report(nlPtr)
      if(.not.associated(nlPtr))then
        call set_error('output_parameters namelist is absent','setup_v5_3')
        return
      endif
      if(fu_str_u_case(fu_content(nlPtr,'source_id')) == 'NO_SOURCE_SPLIT')then
        iSrcIdType = iNoId  ! Mix sources
      elseif(fu_str_u_case(fu_content(nlPtr,'source_id')) =='SOURCE_NAME')then
        iSrcIdType = iSrcNmId
      elseif(fu_str_u_case(fu_content(nlPtr,'source_id')) == 'SOURCE_SECTOR')then
        iSrcIdType = iSrcSectorId
      elseif(fu_str_u_case(fu_content(nlPtr,'source_id')) == 'SOURCE_NAME_AND_SECTOR')then
        iSrcIdType = iSrcNmSectorId
      elseif(index(fu_str_u_case(fu_content(nlPtr,'source_id')), 'TIME_ZONE') > 0)then
        iSrcIdType = iSrcTimeZoneId
      else
        call set_error('Strange source_identification string:' + fu_content(nlPtr,'source_id'), &
                     & 'setup_v5_3')
        return
      endif


      !-------------------------------------------------------------------
      !
      !   The SOURCE TERMS
      !
      ! Now ready to read the sources, which is controlled in the soruce_general
      !
      nlPtr => fu_namelist(nlSetupGrp,'emission_parameters')
      if(.not. associated(nlPtr))then
        call set_error('Cannot find the emission_parameters namelist','setup_v5_3')
        return
      endif

      call msg('Setting source terms')

      call set_source_terms(em_source, nlPtr, iSrcIdType, ifEmissionGiven, &
                          & simRules%chemicalRules%expected_species, &
                          & simRules%dynamicsRules%simulation_type, .false.)
      if(error)return


      !--------------------------------------------------------------------------
      !
      ! Time variables
      ! 1. Having the forward/inverse run type - determine the time variables
      ! 2. Default values are set for the operational radioactivity setup.
      ! 3. User can define start_time and end_time rather than *_time and duration.
      !
      nlPtr => fu_namelist(nlSetupGrp,'general_parameters')

      !
      ! Forward/inverse type of run
      ! 
      select case(fu_str_u_case(fu_content(nlPtr, 'direction_in_time')))
        case('FORWARD')
          time_sign = 1      ! Forward direction
        case('BACKWARD', 'INVERSE')
          time_sign = -1     ! Inverse task - adjoint model version
        case default
          CALL set_error('Unknown time direction:' + fu_content(nlPtr, 'direction_in_time'),'setup_v5_3')
          RETURN
      end select
      !
      ! Model time step
      !
      line = fu_content(nlPtr,'time_step')   !  Model time step
      if(fu_str_u_case(line) == 'DEFAULT')then
        simRules%timestep = fu_set_named_interval('15 MIN') * real(time_sign)
      else
        simRules%timestep = fu_set_named_interval(line) * real(time_sign)
      endif
      if(error)then
        call set_error(fu_connect_strings('Strange model time_step:',line),'setup_v5_3')
        return
      endif

      ! Simualtion period
      !
      line = fu_content(nlPtr, 'start_time')  ! mandatory
      IF(TRIM(line)=="-")THEN
        !
        ! Start time is determined from the source
        !
        if(ifEmissionGiven)then
          if(time_sign == 1)then
            simRules%startTime = fu_round_closest(fu_start_time(em_source), &
                                                & simRules%timestep, backwards)  ! Forward run
          else
            simRules%startTime = fu_round_closest(fu_end_time(em_source), &
                                                & simRules%timestep * (-1.), forwards)  ! Inverse run
          endif
        else
          call set_error('Emission source is void, simulation time must be defined explicitly', &
                       & 'setup_v5_3')
          return
        endif
        call msg("start_time taken from sources:  "//trim(fu_str(fu_start_time(em_source))))
        call msg("start_time used for the   run:  "//trim(fu_str(simRules%startTime)))
      ELSE
        simRules%startTime = fu_io_string_to_time(line)
      END IF
      IF (error) RETURN
      !
      ! Duration of the run defined either explicitly or via end_time
      !
      line = fu_content(nlPtr,'computed_period')   ! Simulation duration 
      if(fu_str_u_case(line) == 'DEFAULT')then
        simRules%periodToCompute = fu_set_named_interval('48 hr') * real(time_sign)
      elseif(len_trim(line) > 1)then
        simRules%periodToCompute = fu_set_named_interval(line, simRules%startTime) * real(time_sign)
      else
        !
        ! Void duration. Take end_time.
        !
        line = fu_content(nlPtr,'end_time')
        if(len_trim(line) > 1)then
          simRules%periodToCompute = fu_io_string_to_time(line) - simRules%startTime
          if(simRules%periodToCompute * real(time_sign) < zero_interval)then
            call set_error('start_time and end_time are mixed-up or computed_period has wrong sign', &
                         & 'setup_v5_3')
            call msg('Start_time is a start of the simulations, which may go back in time')
            call msg('in case of inverse run. Computed_period is always a positive interval')
          endif
        else
          call set_error('Neither computed_period nor end_time are given','setup_v5_3')
        endif
      endif
      if(error)return

      !
      ! Check whether the special file is to be created with progress of the task in %
      !
      simrules%chProgressFNmTemplate = fu_content(nlPtr, 'progress_file_name')
      simrules%ifWriteProgressFile = len_trim(simrules%chProgressFNmTemplate) > 1


      !-----------------------------------------------------------------
      !
      ! Initial and boundary conditions. So far, only rules are to be set
      !
      nlPtr => fu_namelist(nlSetupGrp,'initial_and_boundary_conditions')

      call set_ini_boundary_rules(nlPtr, nlStandardSetup, simRules%IniBoundaryRules, simRules%startTime)


      !--------------------------------------------------------------
      !
      ! Initialisation of the output rules and checking their consistency with ifEmissionGiven
      !
      nlPtr => fu_namelist(nlSetupGrp,'output_parameters')
      if(.not. associated(nlPtr))then
        call set_error('Cannot find the output_parameters namelist','setup_v5_3')
        return
      endif

      call msg('Intialising output rules')
      call set_outDef_basic(nlPtr, nlStandardSetup, OutDef, simRules%timestep, &
                          & simrules%dynamicsRules%simulation_type, simRules%chCaseNm)
      if(error)return

      call msg('Basic output parameters intialized from namelist')
      call report(outDef)

      !
      ! Set parameters of pollution cloud - directly from the namelist
      ! 
      if(ifEmissionGiven)then
        call msg('Initialising pollution cloud')
        CALL init_pollution_cloud(cloud, &                    ! cloud to initialiaze
                                & em_source, &                ! emission source
                                & simrules%chemicalRules, &   ! chemical rules
                                & simRules%startTime, &       ! start of the run
                                & simRules%timestep, &        ! main time step
                                & fu_timestep(OutDef))        ! output time step
        IF (error) RETURN
      else
        if(fu_ifRunDispersion(OutDef))then
          call set_error('Output definitions require dispersion output but source is void', &
                       & 'setup_v5_3')
        endif
      endif

!      OutDefPtr => OutDef
      if(error)return

      !
      ! Set weather data rules wdr:
      !
      nlPtr => fu_namelist(nlSetupGrp,'meteo_parameters')
      if(fu_fails(associated(nlPtr),'Cannot find the meteo_parameters namelist','setup_v5_3'))return
      wdr = fu_set_wdr(nlPtr, nlStandardSetup, &
                     & simRules%startTime, simRules%periodToCompute, &
                     & fu_area_from_grid(output_grid))
      IF (error) RETURN
      wdrPtr => wdr
      if(len_trim(fu_content(nlPtr,'meteo_data_time_shift')) > 0)then
        call set_meteo_time_shift(wdr, fu_set_named_interval(fu_content(nlPtr,'meteo_data_time_shift')))
        call msg_warning('Shifted meto time in the run:' + fu_str(fu_meteo_time_shift(wdr)),'setup_v5_3')
      endif
    
      !
      ! Set diagnostic rules: vertical wind, etc.
      !
      !w is always on half-levels
      call set_diagnostic_rules(nlPtr, nlStandardSetup, simrules%diagnosticRules)
      if (error) return

      !
      ! Set optical rules. The part may be not mandatory and not even needed.
      ! Therefore, no alarm is set now if the namelist does not exist.
      ! Later we will see if it is needed or not
      !
      pOpticalRules => fu_optical_rules(simRules%chemicalRules)
      if(error)return

      nlPtr => fu_namelist(nlSetupGrp, 'optical_density_parameters')
      if(error)then
        call unset_error('setup_v5_3')
        nullify(nlPtr)
      endif

      if(associated(nlPtr))then
        call set_optical_density_rules(pOpticalRules, nlPtr)
        if(error)return
      else
        call set_missing(pOpticalRules)
      endif

      !
      ! Data assimilation list
      !
      nlPtr => fu_namelist(nlSetupGrp, 'data_assimilation')
      if (associated(nlPtr)) then
        if(fu_fails(simrules%dynamicsRules%simulation_type == eulerian_flag, &
                          & 'Data assimilation run requires Eulerian dynamics','setup_v5_3'))return
        call msg('Setting DA rules')
        call set_da_rules(nlptr, nlStandardSetup, simrules%darules, &
                        & simrules%startTime, simrules%periodToCompute)
        if (simrules%darules%defined) call msg('DA rules are defined')
        if (error) return
        if (len_trim(fu_content(nlPtr,'time_shift_stdev')) > 0) then
          if (.not. (fu_meteo_time_shift(wdr) == zero_interval))  then
            call set_error('Both meteo data time shift and time shift stdev given', &
                 & 'setup_v5_3')
          else
            call create_perturbed_meteo_time(simRules%timestep, fu_set_named_interval(fu_content(nlPtr, 'time_shift_stdev')), wdr)
            call msg_warning('Shifted meteo time in the run:' + fu_str(fu_meteo_time_shift(wdr)),'setup_v5_3')
          endif
        endif
      else
        simRules%daRules%defined = .false.
      end if
      
    end subroutine setup_v5_3



    !==============================================================
    subroutine expand_paths(nlptr, filename_setup)
      !
      ! Process the (possible) filepaths in the setup namelist:
      ! replace environment variables
      ! replace ^ with the directory setupfile is in
      ! Careful with possibly many items with same names
      !
      implicit none
      type(Tsilam_namelist), pointer :: nlptr
      character(len=*), intent(in) :: filename_setup
      
      integer :: ind_item, num_items
      type(Tsilam_nl_item_ptr), pointer :: item
      character(len=worksize_string) :: content

      if (fu_fails(associated(nlptr), 'nlptr not associated', 'expand_paths')) return
      call msg('Expanding paths: ' // trim(filename_setup))

      do ind_item = 1, fu_nbr_of_items(nlptr)
        item => fu_get_item(nlptr, ind_item)
        content = fu_content(item)
        content = fu_expand_environment(content)
        if (error) return
        content = fu_extend_grads_hat(content, filename_setup)
        if (error) return
        call replace_namelist_item(item, fu_name(item), fu_name(item), content)
        if (error) return
      end do
      
    end subroutine expand_paths

    subroutine create_perturbed_meteo_time(timestep, time_shift_stdev, wdr)
      implicit none
      type(silja_wdr), intent(inout) :: wdr
      type(silja_interval), intent(in) :: timestep, time_shift_stdev

      type(silja_interval) ::time_delta
      real, dimension(1) :: time_shift
      integer, dimension(8) :: values

      !call date_and_time(VALUES=values)
      !call init_random_seed(smpi_global_rank*values(8)*values(7))
      call init_random_seed(smpi_global_rank)                                                                                                                                                                                                                                                                              
      call random_normal(time_shift)

      time_shift(1) = time_shift(1) * fu_mins(time_shift_stdev)
      time_delta = fu_set_interval_min(int(fu_mins(timestep) * nint(time_shift(1) / fu_mins(timestep))))

      call set_meteo_time_shift(wdr, time_delta)
      if(error)return

    end subroutine create_perturbed_meteo_time

  END SUBROUTINE start_silam_v5



  ! ***************************************************************
  ! ***************************************************************
  !
  !
  !     Private subroutines.
  !
  !
  ! ***************************************************************
  ! ***************************************************************


  SUBROUTINE silam_v5(simRules, &
                    & em_source, & ! initial state
                    & cloud, &
                    & wdr, & 
                    & daPointers, &
                    & outDef, &
                    & nlSetupGrp, &
                    & had_error) 

    ! The main dispersion routine. Does the follwing:
    ! - Initializes all the stuff
    ! - Executes the main time cycling
    ! - Performs the output
    ! IMPORTANT: This routine uses direct methods of access to data, which 
    ! speeds up the simulations. The payment is - complicated initialization 
    ! scheme and data acquisition applied here.
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(general_dispersion_rules), intent(inout) :: simRules
    TYPE(silam_source), pointer :: em_source
    TYPE(silam_pollution_cloud), pointer :: cloud
    TYPE(silja_wdr), pointer :: wdr
    TYPE(silam_output_definition), pointer :: outDef
    type(Tsilam_namelist_group), pointer :: nlSetupGrp
    type(observationPointers), intent(inout) :: daPointers
    logical, intent(out) :: had_error

    ! Local declarations:
    !
    TYPE(silja_shopping_list) :: met_dyn_shopping_list, met_dyn_full_shopping_list, &
                               & met_st_shopping_list, met_st_full_shopping_list, &
                               & disp_dyn_shopping_list, disp_st_shopping_list
    TYPE(silja_time) :: now, out_time, nwp_lim1, nwp_lim2, nwp_lim1_prev, nwp_lim2_prev, &
                      & ParticleResetTime
    INTEGER :: i, t, grib_unit, iSource, list, met_src, &
             & u_ind, v_ind, nParticlesToReset
    CHARACTER (LEN=fnlen) :: command_string = ' '
    LOGICAL :: first_step = .true., ifOK, lTmp
    ! quantities can have species attached, so size can be large
    INTEGER, DIMENSION(:), pointer :: q_met_dyn, q_met_st, q_disp_dyn, q_disp_st, q_out_met, q_out_disp
    type(Tsilam_namelist), pointer :: nlPtr
    type(Tmass_map), pointer :: pOptDepth
    type(wdr_ptr), dimension(:), pointer :: wdrPtr
    real, dimension(:), pointer :: fTmp
    type(t_tla_trajectory) :: tla_traj
    !
    ! Data structures: meteorological and dispersion stacks for field data
    ! and cocktail_map for substance-specific data (masses, emission fluxes,
    ! concentrations, depositions, etc.)
    !
    type(mini_market_of_stacks), target :: meteoMarket, dispersionMarket, &
                                         & boundaryMarket, outputMarket
    type(mini_market_of_stacks), pointer :: meteoMarketPtr, dispersionMarketPtr, &
                                          & boundaryMarketPtr, outputMarketPtr

    TYPE(Tfield_buffer), POINTER :: meteo_ptr, disp_buf_ptr, out_buf_ptr

    TYPE(Tmass_map_ptr), dimension(:), pointer :: ptrCocktailMap ! Substance-split data

    type(Tmeteo_input), target :: Meteo_input
    type(Tmeteo_input), pointer :: pMeteo_input
    type(model_container) :: model_pack
    type(silam_trajectory_set) :: traj_set  
    
!    !
!    ! Temporary for Eulerian advection
!    !
!    type(Tmass_map), pointer :: pDispFlds, pDispFldsXC, pDispFldsYC, pDispFldsZC


    call msg('')
    call msg('')
    call msg('Starting SILAM v.5 dispersion run')

    !
    ! Clean out the progress files from previous runs, if any
    ! make zero progress
    !
    if(simrules%ifWriteProgressFile) call init_progress_file(simrules%chProgressFNmTemplate)

    !
    ! Stupidity check: does the source emit withing the computation period ?
    ! Do not forget to distinguish between forward and adjoint runs
    !
    if(fu_ifRunDispersion(OutDef))then
      if(fu_period_to_compute(wdr) > zero_interval)then
        if((fu_start_time(wdr) > fu_end_time(em_source)) .or. &
         & (fu_start_time(wdr)+fu_period_to_compute(wdr) < fu_start_time(em_source)))then
           call set_error('Source (' + &
                        & fu_str(fu_start_time(em_source))+'-'+fu_str(fu_end_time(em_source)) + &
                        & ') is fully outside the computation time(' + fu_str(fu_start_time(wdr)) + &
                        & '-' + fu_str(fu_start_time(wdr)+fu_period_to_compute(wdr))+')','silam_v5')
           return
         endif
      else
        if((fu_start_time(wdr) < fu_start_time(em_source)) .or. &
         & (fu_start_time(wdr)+fu_period_to_compute(wdr) > fu_end_time(em_source)))then
           call set_error('Source (' + &
                        & fu_str(fu_start_time(em_source))+'-'+fu_str(fu_end_time(em_source)) + &
                        & ') is fully outside the computation time:' + fu_str(fu_start_time(wdr) + &
                        & fu_period_to_compute(wdr)) + '-' + fu_str(fu_start_time(wdr)),'silam_v5')
           return
         endif
      endif
    endif  ! ifRunDispersion

    meteoMarketPtr => meteoMarket
    dispersionMarketPtr => dispersionMarket
    boundaryMarketPtr => boundaryMarket
    pMeteo_input => meteo_input
    outputMarketPtr => outputMarket


    ! Force the number of sources to 2 or 3 if data assimilation is used
    !
    !----------------------------------------
    nlPtr => fu_namelist(nlSetupGrp, 'data_assimilation')
    !call msg('Not setting 3 sources for D/A')
    call set_3_sources(nlptr)

    !----------------------------------------
    !
    ! Collect data needs from all the model units
    ! Important: parts of the model ask only those variables
    ! which are mandatory for them. But io_server in its
    ! output_input_needs can ask also for desirable but not
    ! mandatory variables.
    ! 
    ! There are two types of quantities: dynamic and fixed-in-time. The first ones
    ! are to be in dynamic stack, while the fixed ones should be read once
    ! and stored in permanent stack - done in set_physiography routine
    !
    ! Another trick for the dynamic fields
    ! Temporarily store the requested quantities into the input_shopping_list.
    ! Based on them, the GRIB analyser will set the true input_shopping_list,
    ! which will be used for the data acquisition. The full_shopping_list
    ! will contain all fields needed for the model run and used for making 
    ! the derived fields.
    !
    ! In addition, we should know the number of model fields and cocktail maps 
    ! to be placed in dispersion_stack of dispersion_buffer. They will be 
    ! initialized a bit later, by each of the model units. Do NOT mix them with the
    ! meteo quantities ! They have nothing common.
    !
    ! Some of these quantities can have attached species, then the size grows quickly and random
    ! max_smth is not good. Also, this sub is never left, so its stack is costly
    !
    q_met_dyn => fu_work_int_array()
    q_met_st => fu_work_int_array()
    q_disp_dyn => fu_work_int_array()
    q_disp_st => fu_work_int_array()
    q_out_met => fu_work_int_array()
    q_out_disp => fu_work_int_array()
    
    q_met_dyn(1) = int_missing
    q_met_st(1) = int_missing
    q_disp_dyn(1) = int_missing
    q_disp_st(1) = int_missing

    met_dyn_shopping_list = fu_set_shopping_list (met_src_missing, &
                                                & q_met_dyn,&
                                                & fu_start_time(wdr),&
                                                & (fu_start_time(wdr)+simRules%timestep),&
                                                & level_missing, & ! Accept all levels of meteodata
                                                & fu_top_level(wdr)) ! Not used if all levels accepted
    IF (error) RETURN
    met_st_shopping_list = met_dyn_shopping_list  ! Empty too
    disp_dyn_shopping_list = met_dyn_shopping_list  ! Empty too
    disp_st_shopping_list = met_dyn_shopping_list  ! Empty too
    
    !-------------------------------------------------------------------------
    !
    ! Collect the requests for input variables from all model units
    ! and also record all internal model fields, which will be used during the run.
    ! There are two types of input needs - needs of input data and needs in internal fields
    ! Input data are requested from meteorology part, while internal are stored into 
    ! dispersion buffer. Consequently, below the input needs are stored into the 
    ! shopping lists, while internal field requests are sent to global_io_init
    !
    ! The baseline - OutDef
    call get_output_quantities(OutDef, q_out_met, q_out_disp)
    !
    ! Model dispersion part
    !
    if(fu_ifRunDispersion(OutDef)) then

      ! ADVECTION
      !
      if(fu_if_eulerian_present(simRules%dynamicsRules%simulation_type))then
        CALL euler_adv_input_needs(simRules%dynamicsRules%advMethod_Eulerian, &
                                 & ifUseMassfluxes(simRules%diagnosticRules), &
                                 & q_met_dyn, q_met_st, q_disp_dyn, q_disp_st)
      endif
      if(fu_if_lagrangian_present(simRules%dynamicsRules%simulation_type))then
        CALL add_lagr_advection_input_needs(simRules%dynamicsRules%advMethod_Lagrangian, &
                                          & q_met_dyn, q_met_st, q_disp_dyn, q_disp_st)
        ! ... and Langrangian diffusion
        CALL add_random_walk_input_needs(simRules%dynamicsRules%diffusionMethod, &
                                       & q_met_dyn, q_met_st, q_disp_dyn, q_disp_st)
      endif
      if(error)return

      ! SOURCE TERMS
      !
      call add_source_term_input_needs(em_source, q_met_dyn, q_met_st, q_disp_dyn, q_disp_st, wdr)
      if(error)return

      ! TRANSFORMATION & CHEMISTRY
      !
      call add_transformation_input_needs(simRules%chemicalRules, &
                                        & q_met_dyn, q_met_st, q_disp_dyn, q_disp_st, &
                                        & pMeteo_input)
      if(error)return

      call add_deposition_input_needs(q_met_dyn, q_met_st, q_disp_dyn, q_disp_st, &
                                    & simRules%chemicalRules%rulesDeposition, wdr)
      if(error)return

      ! pollution_cloud may need something as well because it has some
      ! short-cuts, which, violating the encapsulation, still make the program faster
      !
      call add_cloud_driver_input_needs(em_source,simRules%chemicalRules, q_met_dyn, q_met_st, &
                                                                        & q_disp_dyn, q_disp_st)
      if(error)return

      !
      ! Diagnostic variables.
      !
      ! ATTENTION. The full outputlist is not yet ready, it is prepared in global_io_uinit
      ! However, here the call is safe because we need just a list of quantitites, not
      ! the full list of variables. And the list of quantities is complete already here. 
      !
      call add_optical_dens_input_needs(q_out_disp, &
                                      & q_met_dyn, q_met_st, &
                                      & q_disp_dyn, q_disp_st, &
                                      & fu_optical_rules(simRules%chemicalRules))
      if(error)return

      ! Boundaries might have input needs to convert quantities in the input files to concentrations
      !
      if(fu_if_eulerian_present(simRules%dynamicsRules%simulation_type))then
         call add_ini_boundary_input_needs(simRules%IniBoundaryRules, q_met_dyn, q_met_st, &
                                                                    & q_disp_dyn, q_disp_st)
         if(error)return
      endif
    endif ! ifRunDispersion
    if(error)return
    !
    ! All input needs are collected. Can add them to the shopping lists
    !
    call add_shopping_quantities(met_dyn_shopping_list, q_met_dyn)
    call add_shopping_quantities(met_st_shopping_list, q_met_st)
    call add_shopping_quantities(disp_dyn_shopping_list, q_disp_dyn)
    call add_shopping_quantities(disp_st_shopping_list, q_disp_st)
    if(error)return
    !
    ! Info collected, get rid of temporaries
    call free_work_array(q_met_dyn)
    call free_work_array(q_met_st)
    call free_work_array(q_disp_dyn)
    call free_work_array(q_disp_st)
    call free_work_array(q_out_met)
    call free_work_array(q_out_disp)
    !
    ! Having model input needs, we can initialize all IO structures,
    ! pre- and post-processors
    !
    call global_io_init(met_dyn_shopping_list, met_dyn_full_shopping_list, &
                      & met_st_shopping_list, met_st_full_shopping_list, &
                      & disp_dyn_shopping_list, disp_st_shopping_list, &
                      & meteoMarketPtr, dispersionMarketPtr, outputMarketPtr, &
                      & meteo_ptr, disp_buf_ptr, out_buf_ptr, &
                      & OutDef, &
                      & traj_set, &  !Container for trajectories
                      & wdr, simRules%diagnosticRules, &
                      & simRules%chemicalRules, &
                      & simRules%dynamicsRules, &
                      & em_source, &   ! initial state
                      & cloud, &
                      & simRules%timestep, &          ! of advection
                      & simRules%ResidenceInterval, & ! Can be negative in adjoint runs
                      & simRules%PeriodToCompute, &
                      & simRules%iComputationAccuracy, &
                      & fu_DA_time_window(simRules%DArules), & ! time window covered by data assimilation
                      & nlSetupGrp)                 ! Group of setup namelists
    if(error)return

    now = fu_start_time(wdr)
    out_time = now
    simrules%ifRunDispersion = fu_ifRunDispersion(OutDef)  ! May be, model quantities disappeared

    !! Should be run before meteo
    call update_btypes(wholeMPIdispersion_grid,simRules%IniBoundaryRules)
    !
    ! Get initial meteo
    !
    call msg("Acquiring meteo after init (in silam_v5)")
    lTmp=.True.  !!Force reading of new meteo
    call prepare_meteo(wdr, simrules%iComputationAccuracy, simRules%diagnosticRules, &
                     & met_dyn_shopping_list, met_dyn_full_shopping_list,  &
                     & disp_dyn_shopping_list, disp_st_shopping_list, &
                     & fu_output_dyn_shopping_list(OutDef, now),&
                     & fu_output_st_shopping_list(OutDef, now), & 
                     & meteo_ptr, disp_buf_ptr, out_buf_ptr, &
                     & now, simrules%timestep, &
                     & meteoMarketPtr, dispersionMarketPtr, outputMarketPtr, lTmp, &
                     & simrules%ifRunDispersion)
    !
    ! For Eulerian advection we have to define the minimum mass threshold, below which
    ! computations are not performed
    ! For the Lagrangian advection we have to reset particles every while. 
    !
    if(simrules%ifRunDispersion)then
!      !
!      ! Determine the minimum considered mass in a grid cell. Note the ambiguity of the value:
!      ! near the surface we have thin layers thus higher concentrations are reached at lower
!      ! amount per grid cell. Still, working with actual concentrations per m3 costs a lot
!      ! because this checking will be done millions of times. Therefore we assume the typical
!      ! vertical size of 100m, which is about what we have near the surface.
!      !
!      call set_low_mass_threshold(simRules%chemicalRules, &
!                                & em_source, &
!                                & fu_species_transport(cloud), fu_species_emission(cloud), &
!                                & fu_nbr_of_species_transport(cloud), &
!!                                & nlSetupGrp, &
!                                & simRules%ResidenceInterval, &
!                                & fu_period_to_compute(wdr), &
!                                & simRules%timestep, &
!                                & simRules%iComputationAccuracy, &
!                                & fu_earliest_start(cloud), &
!                                & .false., &
!                                & (/real_missing/), 0, &  ! DA parameters, void here
!                                & fu_if_lagrangian_present(simRules%dynamicsRules%simulation_type))
!      if (error) return
!
!      fTmp => fu_low_mass_threshold(simRules%chemicalRules)
!      do i = 1, fu_nbr_of_species_transport(cloud)
!        if(fTmp(i) <= 0.0)then
!          call msg('BAD Low-mass per grid cell threshold=',i,fTmp(i))
!          call set_error('Non-positive low-mass per grid cell threshold','silam_v5')
!          return
!        endif
!      end do

      !
      ! Create the necessary fields
      !
      if(fu_if_eulerian_present(simRules%dynamicsRules%simulation_type))then
        call InitEulerAdvectionFields(simRules%dynamicsRules%advMethod_Eulerian, &
                                    & simRules%dynamicsRules%advection_variant, &
                                    & simRules%dynamicsRules%smoother_factor,& 
                                    & fu_nbr_of_sources(cloud), &
                                    & fu_nbr_of_species_transport(cloud), &
                                    & fu_nbr_of_species_aerosol(cloud), &
                                    &  simRules%dynamicsRules%ifMolecDiff, &
                                    &  simRules%dynamicsRules%ifSubgridDiff )
      endif
      if(fu_if_lagrangian_present(simRules%dynamicsRules%simulation_type))then
        call InitLagrAdvectionFields(simRules%dynamicsRules%advMethod_Lagrangian, &
                                    & fu_nbr_of_species_transport(cloud), &
                                    & fu_nbr_of_species_aerosol(cloud))
      endif
      if(error)return

!      call InitChemistryFields(dispersionMarketPtr, simRules%chemicalRules, wdr, &
!                             & fu_cocktail_transport(cloud))
!      if(error)return

      call init_cloud_internalFields(dispersionMarketPtr, simRules%chemicalRules, wdr)
      if(error)return

    endif   ! ifRunDispersion

    !
    ! Set the boundary conditions - now here, maybe eventually elsewhere. Note that this is
    ! called always, since part of the boundaries do not depend only on the rules, but whether the 
    ! run is global or not.
    !
    if(simrules%ifRunDispersion)then
      if(fu_if_eulerian_present(simRules%dynamicsRules%simulation_type))then
        call init_boundary_structures(simRules%IniBoundaryRules, &
                                    & wdr, &
                                    & fu_species_transport(cloud), fu_nbr_of_species_transport(cloud), &
                                    & fu_concMM_ptr(cloud), &
                                    & fu_DA_time_window(simRules%DArules), &
                                    & boundaryMarketPtr, &
                                    & fu_boundaryStructures(cloud), &  ! boundary structures themsemves
                                    & fu_boundaryBuffer(cloud))       ! buffer for data flow to advection
        if (fu_fails(.not.error,"error after init_boundary_structures", 'silam_v5')) return
      endif
      !
      ! Having everything created, we can apply initial conditions.
      !
      if(defined(simRules%IniBoundaryRules) .and. &
       & fu_if_eulerian_present(simRules%dynamicsRules%simulation_type))then

        call msg('Setting initial conditions')

        call set_cloud_initial_conditions(cloud, &            ! mass maps, boundary structures
                                        & disp_buf_ptr, &
                                        & simRules%IniBoundaryRules, &
                                        & now, simRules%timestep, &
                                        & meteoMarketPtr, dispersionMarketPtr)
        if(error)return

        call msg("Initial conditions set")

        !
        ! If emission scaling data was loaded, it will be used to set the emission processor.
        !
        call set_emission_processor(dispersionMarketPtr, &
                                  & now, &
                                  & fu_emisMM_ptr(cloud), &
                                  & fu_species_transport(cloud), &
                                  & processor_scaling_flag, &
                                  & fu_emission_processor_ptr(cloud))
        if(error)return

#ifdef DEBUG
        call report(fu_emission_processor_ptr(cloud), 'emission processor at the start')
#endif

      endif
    endif
    !
    ! Check the species in the cloud: there must be no garbage
    !
    call check_cloud_species(cloud)
    
    command_string = 'Overall_initialisation'
    call stop_count(chCounterNm = command_string)

!    !-----------------------------------------------------------------
!    !
!    ! Set data assimilation rules. I was planning to include this earlier, but DA needs to know
!    ! the dispersion grid as well as cocktails, so it has to be around here.
!
!    nlPtr => fu_namelist(nlSetupGrp, 'data_assimilation')
    
    ! So we are initialized. If defined, data assimilation will happen here.
    !
    if (simRules%daRules%defined) then
      call msg('Starting data assimilation')
      call pack_model(cloud, em_source, simrules, outDef, &
                    & meteo_ptr, disp_buf_ptr, out_buf_ptr, &
                    & meteoMarketPtr, dispersionMarketPtr, boundaryMarketPtr, outputMarketPtr, &
                    & wdr, met_dyn_shopping_list, met_dyn_full_shopping_list, &
                    & disp_dyn_shopping_list, disp_st_shopping_list, &
                    & pMeteo_input, daPointers, tla_traj, &
                    & model_pack)
      !call test_var_da(model_pack)
      !stop
      if (simrules%darules%method == flag_enkf) then
        call msg('Starting ENKF')
        call run_enkf(model_pack)
      else if (simrules%darules%method == flag_enks) then
        call msg('Starting ENKS')
        call run_enks(model_pack)
      else
        call msg('Starting xDVAR')
        call run_xdvar(model_pack, ifOk)
      end if
      !call run_4dvar(model_pack)

      if (error) then
        had_error = .TRUE.
        call unset_error('silam_v5 after assimilation run')
        !!! Only unsets in non-MPI version
      elseif (ifOk) then
        had_error = .FALSE.
        call msg('Finalising files after data assimilation run')
      else
        had_error = .TRUE.
        call msg_warning('assimilation run failed', 'silam_v5')
      end if
    else
      !
      ! Ordinary run, no data assimilation
      !
      call run_dispersion(cloud, em_source, &
                        & wdr, simRules, &
                        & met_dyn_shopping_list, met_dyn_full_shopping_list, &
                        & disp_dyn_shopping_list, &
                        & disp_st_shopping_list, &
                        & outDef, &
                        & meteo_ptr, disp_buf_ptr, daPointers, &
                        & pMeteo_input, out_buf_ptr, &
                        & meteoMarketPtr, dispersionMarketPtr, boundaryMarketPtr, outputMarketPtr, tla_traj)
      !
      ! If we left without errors - the last time step has to be stored into the file
      !
      had_error = error
      IF(error) then
        CALL unset_error('silam_v5 after time loop')
      else
        call msg('Finalising files after time loop')
      endif

    endif ! DA run or normal dispersion 
    !
    ! Last output
    call msg("Final output after-step")

    command_string = 'Output'
    call start_count(chCounterNm = command_string)
     if (simrules%if_make_output) then 
      call collect_output(meteo_ptr,  disp_buf_ptr, out_buf_ptr, cloud, &
                        & simrules%startTime + simrules%periodToCompute, &
                        & OutDef, wdr, simRules%timestep, simrules%dynamicsRules%simulation_type, &
                        & .false., & ! Not the first output
                        & .true.)   ! Very last output  
    end if

      call stop_count(chCounterNm = command_string)
    !
    ! All done, flush buffers, close files, store multi-run information, 
    ! write needed ctls and logs 
    !
    call finalise_output(cloud, OutDef, wdr, em_source)
    !
    ! If we use "infinity", we go until the end of data and the progress never reaches 100. 
    ! So we have to force it at the end of the run.
    !
    if(simrules%ifWriteProgressFile)then
      if(simrules%periodToCompute == fu_set_named_interval('INFINITY'))then
        call write_progress_file(simrules%chProgressFNmTemplate, 100)
      endif
    endif
    
    ! We must not exit clean if there were errors, so others would know that
    ! we are in trouble
    if (had_error) call msg_warning("had_error at the end", "silam_v5")

    had_error = had_error .or. error

  END SUBROUTINE silam_v5


  ! ***************************************************************

  SUBROUTINE set_standard_setup(nlSetup, rulesDynamics, ifOldFileFormat)
    !
    ! Reads the setup file containing a set of internal model parameters
    !
    ! Idea of the file is not to force a user to set the parameters, which
    ! he does not understand. This file is mainly hidden from users
    !
    IMPLICIT NONE

    ! IN-parameters
    type(Tsilam_namelist), pointer :: nlSetup
    logical, intent(in) :: ifOldFileFormat

    ! INOUT-parameters
    type(Tdynamics_rules), intent(inout) :: rulesDynamics

    ! Local variables
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: pItems
    integer :: nItems
    character(len=clen) :: line

    character(len=*), parameter :: sub_name = 'set_standard_setup'
   
    pItems=>null() 
    test_messages = .false. ! For debug information

    if(.not. associated(nlSetup))then
      call set_error('Undefined pointer to standard setup namelist',sub_name)
      return
    endif
    !
    ! Advection: eulerian, lagrangian, and default. 
    ! The treatment of species coming from specific sources can be by either of them. 
    ! If nothing stated, the default will be taken.
    !
    if(ifOldFileFormat)then
      !
      ! Old file format: Lagrangian environment does not exist, Eulerian defined as the only one
      !
      rulesDynamics%advMethod_Lagrangian = int_missing
      rulesDynamics%diffusionMethod = int_missing
      !
      ! Advection method
      !
      line = fu_content(nlSetup,'advection_method')
      if (line=='') then
         line = fu_content(nlSetup,'horizontal_advection_method')
         if (line=='') then
          CALL set_error('Unknown advection_method',  'read_setup_file')
          call msg("advection_method must be one of: EULERIAN_V4/EULERIAN_3D_BULK/EULERIAN_V5")
          call report(nlSetup)
         else
          call msg_warning("horizontal_advection_method is depricated, use advection_method instead", &
               &  'read_setup_file')
         endif
      endif

      select case(line)
        case('NO_ADVECTION')
          rulesDynamics%advMethod_Eulerian = no_advection
        CASE ('EULERIAN_HORIZ_V4', 'EULERIAN_V4')
          rulesDynamics%advMethod_Eulerian = adv_euler_Galperin_v5
          call msg_warning("V4 is depricated, using EULERIAN_V5 instead")
        CASE ('EULERIAN_3D_BULK')
          rulesDynamics%advMethod_Eulerian = adv_euler_Galperin_3d_bulk
        CASE ('EULERIAN_V5')
          rulesDynamics%advMethod_Eulerian = adv_euler_Galperin_v5
        CASE default
          CALL set_error('Unknown old-format advection_method:' +  line, &
                         & 'read_setup_file')
          call msg("horizontal_advection_method must be: EULERIAN_3D_BULK EULERIAN_V5")
          call report(nlSetup)
      end select
      if (fu_content(nlSetup,'vertical_advection_method') /= '') then
          call msg("vertical_advection_method is depricated, use advection_method instead")
          call msg_warning("vertical_advection_method was ignored",  'read_setup_file')
      endif

    else
      !
      ! New file format: Lagrangian and Eulerian environments are different
      !
      select case(fu_content(nlSetup,'advection_method_lagrangian'))
        case('LAGRANGIAN_WIND_MIDPOINT_2D')
          rulesDynamics%advMethod_Lagrangian = adv_lagrange_2D_wind_midpoint
        CASE ('LAGRANGIAN_WIND_MIDPOINT_3D')
          rulesDynamics%advMethod_Lagrangian = adv_lagrange_3D_wind_midpoint
        CASE ('LAGRANGIAN_WIND_ENDPOINT_3D')
          rulesDynamics%advMethod_Lagrangian = adv_lagrange_3D_wind_endpoint
        CASE ('LAGRANGIAN_IMPLICIT_3D')
          rulesDynamics%advMethod_Lagrangian = adv_lagrange_3D_implicit
        case default
          CALL msg_warning('Unknown advection_method_lagrangian:'+&
                     & fu_content(nlSetup,'advection_method_lagrangian'), &
                         & 'read_setup_file')
          call msg("advection_method_lagrangian must be one of:")
          call msg("LAGRANGIAN_WIND_MIDPOINT_2D")
          call msg("LAGRANGIAN_WIND_MIDPOINT_3D")
          call msg("LAGRANGIAN_WIND_ENDPOINT_3D")
          call msg("LAGRANGIAN_IMPLICIT_3D")
          call msg("Setting to missing and hoping for the best...")
          rulesDynamics%advMethod_Lagrangian = int_missing
      end select

      line = fu_content(nlSetup,'advection_method_eulerian')
      if (line=='') then
         line = fu_content(nlSetup,'horizontal_advection_method_eulerian')
         if (line=='') then
          CALL set_error('Unknown advection_method_eulerian',  'read_setup_file')
          call msg("advection_method_eulerian must be one of: EULERIAN_V4/EULERIAN_3D_BULK/EULERIAN_V5")
          call report(nlSetup)
         else
          call msg_warning("horizontal_advection_method_eulerian is depricated, use advection_method_eulerian instead", &
               &  'read_setup_file')
         endif
      endif

      select case(line)
        case('NO_ADVECTION')
          rulesDynamics%advMethod_Eulerian = no_advection
        CASE ('EULERIAN_HORIZ_V4', "EULERIAN_V4")
          rulesDynamics%advMethod_Eulerian = adv_euler_Galperin_v5
          call msg_warning("V4 is depricated, using EULERIAN_V5 instead")
        CASE ('EULERIAN_3D_BULK')
          rulesDynamics%advMethod_Eulerian = adv_euler_Galperin_3d_bulk
        CASE ('EULERIAN_V5')
          rulesDynamics%advMethod_Eulerian = adv_euler_Galperin_v5
        CASE default
          CALL set_error('Unknown old-format advection_method_eulerian' +  line, &
                         & 'read_setup_file')
          call msg("horizontal_advection_method_eulerian must be: EULERIAN_V4 EULERIAN_3D_BULK EULERIAN_V5")
          call report(nlSetup)
      end select
      if (fu_content(nlSetup,'vertical_advection_method_eulerian') /= '') then
          call msg("vertical_advection_method_eulerian is depricated, use advection_method_eulerian instead")
          call msg_warning("vertical_advection_method_eulerian was ignored",  'read_setup_file')
      endif

      !
      ! Random walk
      !
      select case(fu_content(nlSetup,'random_walk_method'))
        CASE ('NONE')
          rulesDynamics%diffusionMethod = void_rw
        CASE ('FIXED')
          rulesDynamics%diffusionMethod = fixed_rw
        CASE ('FULLY_MIXED')
          rulesDynamics%diffusionMethod = fully_mixed_abl_rw
        CASE ('BULK_GAUSSIAN')
          rulesDynamics%diffusionMethod = bulk_gaussian_rw
        CASE default
          call msg("Strange random_walk_method ="+fu_content(nlSetup,'random_walk_method'))
          call report(nlSetup)
          call set_error('Unknown random_walk_method. Must be NONE/FIXED/FULLY_MIXED/BULK_GAUSSIAN',sub_name)
          return
      end select
    endif ! if old file format

    select case(fu_content(nlSetup,'advection_method_default'))
        CASE ('EULERIAN')
          rulesDynamics%advectionType_default = eulerian_flag
        CASE ('LAGRANGIAN')
          rulesDynamics%advectionType_default = lagrangian_flag
        CASE default
          if (rulesDynamics%advMethod_Lagrangian == int_missing) then
               call msg_warning("No lagrangian method was specified, &
                  & using eulerian advection_method_default", "set_standard_setup")
               rulesDynamics%advectionType_default = eulerian_flag
            else
              CALL set_error('Default advection must be EULERIAN or LAGRANGIAN, not:' + &
                         & fu_content(nlSetup,'advection_method_default'), sub_name)
                  return
            endif
    end select

    ! Galperin v5-specific parameters
    select case(fu_content(nlSetup,'mass_distributor'))
      case ("TRIANGLE_SLAB")
            rulesDynamics%advection_variant = advect_tri
      case ("RECTANGLE_SLAB")
            rulesDynamics%advection_variant = advect_rect
      case ("STEP_SLAB")
            rulesDynamics%advection_variant = advect_step
      case default
            call msg("No or strange mass_distributor:" + fu_content(nlSetup,'mass_distributor'))
            call msg("Can be TRIANGLE_SLAB, RECTANGLE_SLAB, or STEP_SLAB. Using rectangle")
            rulesDynamics%advection_variant = advect_rect
    end select

    rulesDynamics%ifMolecDiff = (fu_content(nlSetup,'grav_separation') == 'YES')
    rulesDynamics%ifSubgridDiff = (fu_content(nlSetup,'diffuse_vert_cm') == 'YES')

    rulesDynamics%smoother_factor=fu_content_real(nlSetup,'smoother_factor')
    if (rulesDynamics%smoother_factor < 0.) then
      call msg("No or strange smoother_factor:"+fu_content(nlSetup,'smoother_factor'))
      call msg("Setting smoother_factor = 1. (no smoothing)")
      rulesDynamics%smoother_factor = 1.
    endif
    
    !! Cloud geometry (cell size or air_mass or ones tracer)
    line = fu_str_l_case(fu_content(nlSetup, 'cloud_metric'))
    select case (line)
      case ('cell_size')
        rulesDynamics%cloud_metric = cloud_metric_geometry_flag
      case ('cell_mass')
        rulesDynamics%cloud_metric =  cloud_metric_cellmass_flag
      case ('ones_mass')
        rulesDynamics%cloud_metric = cloud_metric_ones_flag
      case ('')
        call msg_warning("cloud_metric not found in the namelist, using SIZE", sub_name)
        rulesDynamics%cloud_metric = cloud_metric_geometry_flag
      case default
        call set_error("Strange cloud_metric in the namelist: '"//trim(line)//"'", sub_name)
        return
    end select 
   
    !
    ! Read the chemical database
    !
    call init_chemical_materials(fu_expand_environment(fu_content(nlSetup,'chemical_database_fnm')), &
                               & fu_expand_environment(fu_content(nlSetup,'nuclide_database_fnm')))
    if(error)return

    call init_optical_metadata(nlSetup)
    if (error) return

    ! Set the standard cocktail descriptions
    !
    call get_items(nlSetup, 'standard_cocktail_fnm', pItems, nItems)
    if(error)return
    call read_standard_cocktail_descrs(pItems, nItems)
    if(error)return


!call msg('Creating the GRIB-2 code table')
!call create_grib2_code_table()

    !
    ! Debug messages
    !
    SELECT CASE (fu_content(nlSetup,'print_debug_info'))
      CASE ('DEBUG_INFO_YES')
        test_messages = .true.
      CASE ('DEBUG_INFO_NO')
        test_messages = .false.
      case default 
        call msg_warning('Debug switch may be DEBUG_INFO_YES or DEBUG_INFO_NO', sub_name)
        test_messages = .false.
    end select

  END SUBROUTINE set_standard_setup



!subroutine create_grib2_code_table()
!  !
!  ! Just makes a unified table out of a set of files given in the definitions directory of GRIB-2 package
!  !
!  implicit none
!
!  ! Local variables
!  integer :: iStat, uFNmParamId, uFNmName, uFNmShortName, uFNmUnit, iParamId, levType, qSILAM, iTmp, &
!           & iDiscipline, iParamCategory, iParamNbr, indexIn, uOut, iItem, iGlobalIndex, tabVersion, centre, model
!  character(len=clen) :: line, chName, chUnit, chShortName
!  integer, dimension(:), pointer :: arParamId, arDiscipline, arParamCategory, arParamNbr
!  character(len=unitNmLen), dimension(:), pointer :: arUnit, arShortName
!  character(len=clen), dimension(:), pointer :: arName
!  type(Tsilam_namelist_ptr), dimension(:), pointer :: nlArPtr
!  character(len=fnlen) :: chFNmParamId, chFNmName, chFNmShortName, chFNmUnit
!  logical :: eof, outOfTable
!  real :: levVal, fModeSILAM
!  type(silja_level) :: levSILAM
!  character(len=12) :: chSubstSILAM
!  type(silam_sp) :: sp, sp2
!  type(Tsilam_nl_item_ptr), pointer :: nlItemPtr
!  type(Tsilam_namelist_ptr) :: nlTmp
!  type(Tsilam_namelist_group), pointer :: nlGrp
!  integer, dimension(:), pointer :: arQuantity
!  !
!  ! prepare the place
!  !
!  allocate(arParamId(1000), arUnit(1000), arShortName(1000), &
!         & arName(1000), nlArPtr(1000), arDiscipline(1000), arParamCategory(1000), arParamNbr(1000), stat=iStat)
!  if(iStat /= 0)then
!    call set_error('Allocaton failed','create_grib2_code_table')
!    stop
!  endif
!  !
!  ! Swallow the stupid files of GRIB-2
!  !
!  chFNmParamId = 'd:\!tools\lib\grib_api\definitions\grib2\paramId.def'
!  chFNmName = 'd:\!tools\lib\grib_api\definitions\grib2\name.def'
!  chFNmShortName = 'd:\!tools\lib\grib_api\definitions\grib2\shortName.def'
!  chFNmUnit = 'd:\!tools\lib\grib_api\definitions\grib2\units.def'
!  uFNmParamId = fu_next_free_unit()
!  open(unit=uFNmParamId, file=chFNmParamId, status='old')
!  uFNmName = fu_next_free_unit()
!  open(unit=uFNmName, file=chFNmName, status='old')
!  uFNmShortName = fu_next_free_unit()
!  open(unit=uFNmShortName, file=chFNmShortName, status='old')
!  uFNmUnit = fu_next_free_unit()
!  open(unit=uFNmUnit, file=chFNmUnit, status='old')
!  
!  eof = .false.
!  sp%sp => fu_work_string()
!  sp2%sp => fu_work_string()
!  iGlobalIndex = 1
!  
!  do while(.false.) !  (.not.eof)
!    call next_line_from_input_file(uFNmParamId, line, eof)
!    if(error .or. eof)exit
!    if(index(line,'{') > 0)then
!      !
!      ! Opened new parameter
!      !
!      sp%sp = line(2:index(line,'=')-3)
!      read(unit=sp%sp, fmt=*)iParamId
!      call next_line_from_input_file(uFNmName, line, eof)
!      chName = line(2:index(line,'=')-3)
!      call next_line_from_input_file(uFNmShortName, line, eof)
!      chShortName = line(2:index(line,'=')-3)
!      call next_line_from_input_file(uFNmUnit, line, eof)
!      chUnit = line(2:index(line,'=')-3)
!      do iTmp = 1, len_trim(chUnit)
!        if(chUnit(iTmp:iTmp) == ' ') chUnit(iTmp:iTmp) = '_'   ! for the time being
!      end do
!    else
!      call next_line_from_input_file(uFNmName, line, eof)
!      call next_line_from_input_file(uFNmShortName, line, eof)
!      call next_line_from_input_file(uFNmUnit, line, eof)
!      
!      if(index(line,'discipline') > 0)then
!        sp%sp = line(index(line,'=')+1: index(line,';')-1)
!        read(unit=sp%sp, fmt=*)iDiscipline
!        
!      elseif(index(line,'parameterNumber') > 0)then
!        sp%sp = line(index(line,'=')+1: index(line,';')-1)
!        read(unit=sp%sp, fmt=*)iParamNbr
!        
!      elseif(index(line,'parameterCategory') > 0)then
!        sp%sp = line(index(line,'=')+1: index(line,';')-1)
!        read(unit=sp%sp, fmt=*)iParamCategory
!        
!      elseif(index(line,'}') > 0)then
!        !
!        ! List over: write the stuff down
!        !
!        if(iDiscipline > 10 .or. iDiscipline < 0)then
!          call msg('Wrong discipline:',iDiscipline)
!          call set_error('Wrong discipline','create_grib2_code_table')
!          stop
!        endif
!        if(iParamNbr > 255 .or. iParamNbr < 0)then
!          call msg('Wrong parmeter number:',iParamNbr)
!          call set_error('Wrong parameter number','create_grib2_code_table')
!          stop
!        endif
!        if(iParamCategory > 255 .or. iParamCategory < 0)then
!          call msg('Wrong parameter category:',iParamCategory)
!          call set_error('Wrong parameter category','create_grib2_code_table')
!          stop
!        endif
!        arDiscipline(iGlobalIndex) = iDiscipline
!        arParamCategory(iGlobalIndex) = iParamCategory
!        arParamNbr(iGlobalIndex) = iParamNbr
!        arParamId(iGlobalIndex) = iParamId
!        arUnit(iGlobalIndex) = chUnit
!        arShortName(iGlobalIndex) = chShortname
!        arName(iGlobalIndex) = chName
!        nlArPtr(iGlobalIndex)%nl => nlTmp%nl
!        nullify(nlTmp%nl)  ! do not destroy
!        iGlobalIndex = iGlobalIndex + 1
!
!      elseif(index(line,'=') > 0)then
!        !
!        ! Something extra. Cut out the namelist item and story to the namelist
!        !
!        sp%sp = adjustl(line)
!        do iTmp = 1, len_trim(sp%sp)
!          if(sp%sp(iTmp:iTmp) == ';') sp%sp(iTmp:iTmp) = ' '
!        end do
!        if(.not. associated(nlTmp%nl)) nlTmp%nl => fu_create_namelist()
!        read(unit=sp%sp,fmt=*)sp2%sp   ! name
!        call add_namelist_item(nlTmp%nl, sp2%sp, sp%sp(index(sp%sp,'=')+1:))
!      endif
!    endif
!
!  end do  ! cycle over the stupid GRIB-2 files
!  
!  !
!  ! Once the information is consumed from these files, co-locate the parameters
!  ! with the SILAM quantities using the 128 code table and create the namelist group 
!  !
!  nlGrp => fu_create_namelist_group('GRIB2_CODE_TABLE_V5')
!  arQuantity => fu_work_int_array()
!  arQuantity(:) = int_missing
!  
!  outOfTable = .false.
!  indexIn = 1
!  uOut = fu_next_free_unit()
!  
!!  tabVersion = 128    ! ECMWF
!!  centre = 98
!!  model = int_missing
!!  open(unit=uOut, file='d:\!model\2011\silam_v5_0\ini\grib_code_table.txt')  
!  
!!  tabVersion = 1      ! HIRLAM
!!  centre = 96           ! HIRLAM, also 94(DMI) and 86(FMI)
!!  model = int_missing
!!  open(unit=uOut, file='d:\!model\2011\silam_v5_0\ini\grib_code_table_HIRLAM_general.txt')  
!  
!!  tabVersion = 1      ! HIRLAM  SMHI - ISBA is specific
!!  centre = 82
!!  model = int_missing
!!  open(unit=uOut, file='d:\!model\2011\silam_v5_0\ini\grib_code_table_HIRLAM_SMHI.txt')  
!  
!  tabVersion = 130   ! SILAM
!  centre = 86        ! FMI
!  model = 3
!  open(unit=uOut, file='d:\!model\2011\silam_v5_0\ini\grib_code_table_SILAM.txt')  
!  
!  do while(.not. outOfTable)
!    !
!    ! Get one entry of SILAM grib code table
!    !
!    call get_all_params(indexIn, tabVersion, centre, int_missing, outOfTable, &
!                      & iParamId, levType, levVal, qSILAM, levSILAM, chSubstSILAM, fModeSILAM)
!    if(outOfTable)exit
!    call msg(fu_print_silam_quantity(qSILAM))
!    do iTmp = 1, indexIn+1
!      if(arQuantity(iTmp) == qSILAM)exit
!    end do
!    if(iTmp > indexIn)then
!      !
!      ! This SILAM quantity has not yet been met. Create new namelist in the group
!      !
!      arQuantity(indexIn) = qSILAM
!      nlTmp%nl => fu_add_namelist(nlGrp, fu_print_silam_quantity(qSILAM))  
!      call add_namelist_item(nlTmp%nl,'centre', fu_str(centre))
!      call add_namelist_item(nlTmp%nl,'grib1_table_version', fu_str(tabVersion))
!      call add_namelist_item(nlTmp%nl,'grib1_model_identification',fu_str(model))
!      !
!      ! Find the correspondence with the consumed 3D array
!      !
!      eof = .true.
!      do iGlobalIndex = 1, size(arParamId)
!        if(arParamId(iGlobalIndex) == iParamId)then
!          eof = .false.
!          exit 
!        endif
!      end do
!      if(eof)then
!!        iDiscipline = int_missing
!!        iParamCategory = int_missing
!!        iParamNbr = int_missing
!        call add_namelist_item(nlTmp%nl,'discipline',fu_str(int_missing))
!        call add_namelist_item(nlTmp%nl,'parameterCategory',fu_str(int_missing))
!        call add_namelist_item(nlTmp%nl,'parameterNumber',fu_str(int_missing))
!        call add_namelist_item(nlTmp%nl,'name','')
!        call add_namelist_item(nlTmp%nl,'short_name','')
!        call add_namelist_item(nlTmp%nl,'unit','')
!      else
!        call add_namelist_item(nlTmp%nl,'discipline',fu_str(arDiscipline(iGlobalIndex)))
!        call add_namelist_item(nlTmp%nl,'parameterCategory',fu_str(arParamCategory(iGlobalIndex)))
!        call add_namelist_item(nlTmp%nl,'parameterNumber',fu_str(arParamNbr(iGlobalIndex)))
!        call add_namelist_item(nlTmp%nl,'name',arName(iGlobalIndex))
!        call add_namelist_item(nlTmp%nl,'short_name',arShortName(iGlobalIndex))
!        call add_namelist_item(nlTmp%nl,'unit',arUnit(iGlobalIndex))
!        do iItem = 1, fu_nbr_of_items(nlArPtr(iGlobalIndex)%nl)
!          nlItemPtr => fu_get_item(nlArPtr(iGlobalIndex)%nl, iItem)
!          call add_namelist_item(nlTmp%nl, nlItemPtr)
!        end do
!      endif
!      call add_namelist_item(nlTmp%nl,'parameterId',fu_str(iParamId))
!      call add_namelist_item(nlTmp%nl,'level_type',fu_str(levType))
!      call add_namelist_item(nlTmp%nl,'level_value',fu_str(levVal))
!      call add_namelist_item(nlTmp%nl,'silam_quantity_name', fu_print_silam_quantity(qSILAM))
!      call level_to_short_string(levSILAM,sp%sp)
!      call add_namelist_item(nlTmp%nl,'silam_level',sp%sp)
!      call add_namelist_item(nlTmp%nl,'silam_species_io_string = ','')
!      call msg('... added: ',iParamCategory-1,iParamNbr-1)
!    else
!      call msg('met before...')
!      nlTmp%nl => fu_namelist(nlGrp, fu_print_silam_quantity(qSILAM))
!      if(.not. associated(nlTmp%nl))then
!        call set_error('Failed to find namelist','create_grib2_code_table')
!        stop
!      endif
!      sp%sp = fu_content(nlTmp%nl,'level_type')
!      if(levType == accept_all)then
!        write(unit=sp2%sp,fmt=*) trim(sp%sp),'  *'
!      else
!        write(unit=sp2%sp,fmt=*) trim(sp%sp),'  ',levType
!      endif
!      call msg('level_type:' + sp2%sp)
!      call replace_namelist_item(nlTmp%nl,'level_type','level_type', sp2%sp)
!
!      sp%sp = fu_content(nlTmp%nl,'level_value')
!      if(levType == accept_all)then
!        write(unit=sp2%sp,fmt=*) trim(sp%sp),'  *'
!      else
!        write(unit=sp2%sp,fmt=*) trim(sp%sp),'  ',levVal
!      endif
!      call msg('level_value:' + sp2%sp)
!      call replace_namelist_item(nlTmp%nl,'level_value','level_value', trim(sp2%sp))
!
!    endif  ! if the quantity namelist has already been set
!
!    indexIn = indexIn + 1
!  end do
!
!  write(uOut,*)'GRIB2_CODE_TABLE_V5'
!  call write_namelist_group(uOut, nlGrp)
!  write(uOut,*)'END_GRIB2_CODE_TABLE_V5'
!  stop
!  
!end subroutine create_grib2_code_table
!
  
  subroutine exec_tests()
    ! 
    ! Add here whatever tests to execute before SILAM starts.
    ! 
    use photolysis
    use silam_levels
    implicit none
    
    continue

  end subroutine exec_tests


END MODULE dispersion_models

