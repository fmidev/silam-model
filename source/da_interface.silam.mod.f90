module da_interface
  ! 
  ! The lower level subroutines of the DA system. In terms of variational DA, this
  ! module implements the model operator M and applies the observator operator H. The
  ! operator H itself is defined by the observation modules.
  !
  !use dispersion_supplementary
  use observations_in_situ!, only :  observationStationPtr, inSituObservation, columnObservation
!  use correlations
!  use silam_partitioning, only : smpi_ens_size
!  use da_common
  use dispersion_supplementary

  !private
  
  public vector_from_cloud
  public init_control_from_par
  public run_forecast
  public obs_from_mdl
  public get_obs_grad
  public set_background
  public set_da_rules
!  public init_da
  public fu_da_time_window

  public fu_values_ptr
  public fu_initial_ptr
  public fu_emission_ptr
  public fu_species_initial
  public fu_species_emission
  public set_default_cnc
  public pack_model
  public unpack_model

  public adj_test
  public adj_test_obs
  public run_forecast_dummy


  ! Public derived types
  !
  public model_container


  private init_control_from_cloud
  private fu_species_emission_ctrl
  private destroy_control
  private control2cloud
  
  
  type model_container
     ! A convenience type for passing around the model structures. This facilitates the
     ! implementation of top-level data assimilation subroutines.
     type(silam_source), pointer :: source
     type(silja_shopping_list), pointer :: disp_dyn_shopping_list, disp_stat_shopping_list, &
          & meteo_input_dyn_shopping_list, meteo_full_dyn_shopping_list
     type(silja_shopping_list), pointer :: input_shopping_list, full_shopping_list
     type(silam_output_definition), pointer :: out_def
     type(Tfield_buffer), pointer :: meteo_ptr, disp_buf_ptr, output_buf_ptr
     type(mini_market_of_stacks), pointer :: meteo_market_ptr, disp_market_ptr, bc_market_ptr, &
                                           & output_market_ptr
     type(silja_wdr), pointer :: wdr
     type(Tmeteo_input), pointer :: meteo_input
     type(silam_pollution_cloud), pointer :: cloud
     type(observationPointers), pointer :: obs_ptr
     type(general_dispersion_rules), pointer :: rules
     ! Storage for 4DVar tangent linear
     type(t_tla_trajectory), pointer :: tla_traj
     ! a Helper mass map for 3dvar. Allows using a subset of observed species as analysis species.
     type(TMass_map), pointer :: mass_map_adj
     ! perturbations for ensemble methods
     type(t_perturbation), dimension(:), pointer :: perturbations => null()
  end type model_container

  integer, parameter, private :: fill_none = 31001, fill_all = 31002, fill_non_analysis = 31003, &
       & reset_only = 31004

  ! Storage for the initial condition when it is not included in the control variable.
  !
  real, dimension(:,:,:,:,:), private, allocatable, save :: emis_only_init_state


  logical, parameter, public :: da_test = .false.

contains

  subroutine pack_model(cloud, source, simrules, outdef, &
                      & meteo_ptr, dispBufPtr, outputBufPtr, &
                      & meteoMarketPtr, dispersionMarketPtr, BCMarketPtr, outputMarketPtr, &
                      & wdr, meteo_input_dyn_shopping_list, meteo_full_dyn_shopping_list, &
                      & disp_dyn_shopping_list, disp_stat_shopping_list, meteo_input, obs_ptr, &
                      & tla_traj, &
                      & model)
    implicit none
    type(silam_pollution_cloud), pointer :: cloud
    type(general_dispersion_rules), target :: simrules
    type(silam_source), pointer :: source
    type(Tfield_buffer), pointer :: dispBufPtr, meteo_ptr, outputBufPtr
    type(silja_shopping_list), intent(in), target :: &
         & disp_dyn_shopping_list, &
         & disp_stat_shopping_list, &
         & meteo_input_dyn_shopping_list, &
         & meteo_full_dyn_shopping_list
    type(silam_output_definition), pointer :: outDef ! might be removed eventually
    type(silja_wdr), pointer :: wdr
    type(mini_market_of_stacks), pointer :: meteoMarketPtr, dispersionMarketPtr, BCMarketPtr, &
                                          & outputMarketPtr
    type(Tmeteo_input), intent(in), target :: meteo_input
    type(observationPointers), target :: obs_ptr
    type(t_tla_trajectory), intent(in), target :: tla_traj

    type(model_container), intent(out) :: model

    model%cloud => cloud
    model%source => source
    model%disp_dyn_shopping_list => disp_dyn_shopping_list
    model%disp_stat_shopping_list => disp_stat_shopping_list
    model%meteo_input_dyn_shopping_list => meteo_input_dyn_shopping_list
    model%meteo_full_dyn_shopping_list => meteo_full_dyn_shopping_list
    model%out_def => outDef
    model%meteo_ptr => meteo_ptr
    model%disp_buf_ptr => dispBufPtr
    model%output_buf_ptr => outputBufPtr
    model%meteo_market_ptr => meteoMarketPtr
    model%disp_market_ptr => dispersionMarketPtr
    model%output_market_ptr => outputMarketPtr
    model%bc_market_ptr => BCMarketPtr
    model%wdr => wdr
    model%meteo_input => meteo_input
    model%rules => simrules
    model%obs_ptr => obs_ptr
    model%tla_traj => tla_traj
    
    ! Set by only 3d-var if needed:
    nullify(model%mass_map_adj) 
    ! Set by EnKF if needed
    nullify(model%perturbations)

  end subroutine pack_model

  subroutine unpack_model(model, cloud, source, wdr, dispBufPtr, meteo_ptr, outDef, & 
                        & meteo_input_dyn_shopping_list, meteo_full_dyn_shopping_list, &
                        & disp_dyn_shopping_list, disp_stat_shopping_list, & 
                        & meteoMarketPtr, dispersionMarketPtr, BCMarketPtr, meteo_input, tla_traj, perturbations)
    implicit none
    type(silam_pollution_cloud), pointer :: cloud
    type(general_dispersion_rules) :: simrules
    type(silam_source), pointer :: source
    type(Tfield_buffer), pointer :: dispBufPtr, meteo_ptr
    type(silja_shopping_list), pointer :: &
         & disp_dyn_shopping_list, &
         & disp_stat_shopping_list, &
         & meteo_input_dyn_shopping_list, &
         & meteo_full_dyn_shopping_list
    type(silam_output_definition), pointer :: outDef ! might be removed eventually
    type(silja_wdr), pointer :: wdr
    type(mini_market_of_stacks), pointer :: meteoMarketPtr, dispersionMarketPtr, BCMarketPtr
    type(Tmeteo_input), pointer :: meteo_input
    type(t_tla_trajectory), pointer :: tla_traj
    type(t_perturbation), dimension(:), pointer :: perturbations

    type(model_container), intent(in) :: model
    cloud => model%cloud
    source => model%source
    disp_dyn_shopping_list => model%disp_dyn_shopping_list
    disp_stat_shopping_list => model%disp_stat_shopping_list
    meteo_input_dyn_shopping_list => model%meteo_input_dyn_shopping_list
    meteo_full_dyn_shopping_list => model%meteo_full_dyn_shopping_list
    outDef => model%out_def
    meteo_ptr => model%meteo_ptr
    dispBufPtr => model%disp_buf_ptr
    meteoMarketPtr => model%meteo_market_ptr
    dispersionMarketPtr => model%disp_market_ptr
    BCMarketPtr => model%bc_market_ptr
    wdr => model%wdr
    meteo_input => model%meteo_input
    tla_traj => model%tla_traj
    perturbations => model%perturbations
  end subroutine unpack_model


  !************************************************************************************

  subroutine control2cloud(control, model, handle_negatives)
    ! 
    ! Apply the control to the model: set the initial state and/or set the emission
    ! processor.
    implicit none
    type(da_control), intent(in) :: control
    type(model_container), intent(inout) :: model
    ! A trick for PM assimilation in 3dvar mode. Each substance is chekced for existence a
    ! __NEG__ counterpart. If it exists, any negative mass is transferred to that species
    ! multplied by -1.
    logical, intent(in) :: handle_negatives
    !type(silam_pollution_cloud), intent(inout) :: cloud
    
    integer :: ix, iy, iz, isp, nz, nsp, idim, ind, stat, n_emis, isp_transp, nsp_contr
    integer :: nx, ny, offx, offy, gnx, gny, wrk_size 
    integer :: iMPI, jMPI, niMPI, njMPI
    integer, dimension(:), pointer :: iPtr
    integer :: ind_pert_volc, ind_pert_xy
    type(Tmass_map), pointer :: map_c
    real, dimension(:), pointer :: p_init, dx_ptr, dy_ptr, p_emis, emis_data, fWork
    real, dimension(max_levels) :: thickness
    type(Tmass_map), pointer :: map_px, map_py, map_pz, map_sl
    real :: volume, mass_neg, mass_pos
    type(Temission_processor), pointer :: proc_ptr
    type(silja_time), dimension(:,:), pointer :: time_slots
    logical :: need_initial_state, loc_handle_negatives, ifKalmanEnsemble, ifOk, if_use_work
    integer, dimension(max_divisions) :: sendcounts, offsets
    real, dimension (10) :: reportMy, reportAll ! To synchronize reporting

    if (fu_fails(fu_in_physical_space(control), 'Control not in physical space', 'control2cloud')) return

    map_c => fu_concMM_ptr(model%cloud)

    call smpi_get_decomposition(nx, ny, offx, offy, gnx, gny)
    call smpi_get_process_xy_topology(iMPI, jMPI, niMPI, njMPI)

    nz = map_c%n3d
    nsp = map_c%nspecies
    
    thickness(1:nz_dispersion) = (/(fu_layer_thickness_m(fu_level(dispersion_vertical, iz)),&
                                  & iz=1, nz_dispersion)/)
    
    loc_handle_negatives = handle_negatives .and. fu_mode(control) == DA_INITIAL_STATE
    

    ! If the control doesn't completely determine the initial state, apply the default first:
    ! 
    if (debug_level > 0) call msg('fill mode:', fu_initial_fill_mode(model%rules%darules))
    select case(fu_initial_fill_mode(model%rules%darules))
      case (fill_all, fill_non_analysis)
        call reset_cloud_data_structures(model%cloud, ifMomentsToo=.true.)
        ! Apply initial condition where available. In the fill_non_anlaysis case, some
        ! species will be overwritten later.
        if (.not. allocated(emis_only_init_state)) then
          call msg('Loading initial condition...')
          call set_cloud_initial_conditions(model%cloud, &            ! mass maps, boundary structures
                                          & model%disp_buf_ptr, &
                                          & model%rules%IniBoundaryRules, &
                                          & model%rules%startTime, &
                                          & model%rules%timestep, &
                                          & model%meteo_market_ptr, model%disp_market_ptr)

          call set_default_cnc(map_c)
          if (error) return
          call msg('Initial condition: max:', maxval(emis_only_init_state(:, 1, 1:nz, 1:nx, 1:ny)))
          call msg('Initial condition: min:', minval(emis_only_init_state(:, 1, 1:nz, 1:nx, 1:ny)))
          call msg('Initial condition: total:', sum(emis_only_init_state(:, 1, 1:nz, 1:nx, 1:ny)))
        else
          call msg('Setting initial condition from cached...')
          map_c%arm = 0.0
          map_c%arm(:,1,1:nz,1:nx,1:ny) = emis_only_init_state(:,1,1:nz,1:nx,1:ny)
          call msg('Grand total:', sum(map_c%arm(:,1,1:nz,1:nx,1:ny)))
        end if
      case (fill_none)
        continue

      case (reset_only)
        call msg('Reset cloud...')
        call reset_cloud_data_structures(model%cloud, ifMomentsToo=.true.)
        
      case default
        call set_error('Strange fill_mode', 'control2cloud')
    end select

    if (fu_have_initial(fu_mode(control))) then
      
      nsp_contr = fu_nsp_initial(control)
      p_init => null() !! Stupidity check
      
      if (iMPI == 0 .and. jMPI == 0) then !Do scatter mpi rows from control if needed
         p_init => fu_initial_ptr(control)
         sendcounts(1:njMPI) = smpi_y_sizes(1:njMPI)*gnx*nz*nsp_contr
         offsets(1:njMPI)    = smpi_y_offsets(1:njMPI)*gnx*nz*nsp_contr
      else  !Need work for exchane
         sendcounts(1:njMPI)=0 ! Does not matter
         offsets(1:njMPI)=0!
       endif

       wrk_size = ny*gnx*nz*nsp_contr !MPI-row-size 


       !Master needs work only for y exchange
       !Others -- always
       if_use_work = (njMPI > 1 .or. iMPI /= 0 .or.  jMPI /= 0)

       if (if_use_work)  fWork => fu_work_array(wrk_size)

       if (iMPI == 0 .and. njMPI > 1) then !!Do stripe scatter in column 0
         call smpi_scatterv_real(p_init, sendcounts, offsets, fWork, wrk_size,  smpi_adv_y_comm)

       endif
       
       if (if_use_work) p_init => fWork !redirect to our stripe

       !!Can not get X stripe due to wrong order
       !! broadcast MPU row if needed
       if ( niMPI > 1) call smpi_bcast_aray(p_init, wrk_size, smpi_adv_x_comm,0) 

       !Now p_init has our MPI row, 
     
       reportMy(:) = 0.

      if (debug_level > 0) reportMy(1) = sum(map_c%arm)

      iPtr => fu_isp_initial(control)
      mass_neg = 0
      mass_pos = 0
      do iy = 1, ny !Local massmap indices
        do ix = 1, nx
          do iz = 1, nz
            volume = thickness(iz) * fu_cell_size(dispersion_grid, ix, iy)
            do isp = 1, nsp_contr
              isp_transp = iPtr(iSp)   !control%ind_species_init(isp)
              !indexing inside MPI row
              ind = (((iy-1)*gnx + ix-1+offx)*nz + iz-1)*nsp_contr + isp
              !if (p_init(ind) > 0) then
                !call msg('ix, iy', ix, iy)
                !call msg('iz, isp_transp', iz, isp_transp)
                !call msg('pinit', p_init(ind))
              !end if
              if (model%rules%darules%use_log_cloud) then
                map_c%arm(isp_transp, 1, iz, ix, iy) = (exp(p_init(ind)) - 0.02) * volume
              else
                map_c%arm(isp_transp, 1, iz, ix, iy) = p_init(ind) * volume
                if (p_init(ind) < 0) then
                  mass_neg = mass_neg + p_init(ind)*volume
                else
                  mass_pos = mass_pos + p_init(ind)*volume
                end if
               
              end if
         
            end do
          end do
        end do
      end do


      if (loc_handle_negatives) call move_neg_mass(map_c, 'neg_pos')

      if (debug_level > 0) reportMy(2) = sum(map_c%arm) !before truncating

      if (.not. model%rules%darules%allow_negative) map_c%arm = max(map_c%arm, 0.0)
      reportMy(3) = mass_neg
      reportMy(4) = mass_pos
      if (debug_level > 0) reportMy(5) = sum(map_c%arm)

      if (smpi_adv_tasks>1 .and. debug_level > 0) then !report local domain and reduce
        call msg('LOCAL Sum of map_c before control', reportMy(1))
        call msg('LOCAL Sum of map_c before truncating', reportMy(2))
        call msg('LOCAL Total negative mass:', reportMy(3))
        call msg('LOCAL Total positive mass:', reportMy(4))
        call msg('LOCAL Sum of map_c after control', reportMy(5))
      endif

      if (smpi_adv_tasks>1) then
        call smpi_allreduce_add(reportMy, reportAll, smpi_adv_comm, ifOk)
        if (fu_fails(ifOk, "failed MPI COMM report", "control2cloud")) return
      else
        reportAll(:) = reportMy(:)
      endif

      if (debug_level > 0) then
        call msg('Sum of map_c before control', reportAll(1))
        call msg('Sum of map_c before truncating', reportAll(2))
        if (smpi_adv_rank == 0) then
          call msg('Sum of p_init', sum(p_init))
          call msg('Number of negative values:', count(p_init(1:wrk_size) < 0))
        else
          call msg('Sum of p_init (not master) XXXX')
          call msg('Number of negative values (not master) XXXX')
        endif
      endif

      call msg('Total negative mass:', reportAll(3))
      call msg('Total positive mass:', reportAll(4))
      if (debug_level > 0) call msg('Sum of map_c after control', reportAll(5))


      if (if_use_work) call free_work_array(fWork)
    end if ! have initial
    proc_ptr => fu_emission_processor_ptr(model%cloud)

    !if (fu_have_emission_xy(fu_mode(control))) then
    !  if (fu_fails(associated(proc_ptr), 'Emission processor not associated', 'control2cloud')) return
    !  if (fu_fails(defined(proc_ptr), 'Emission processor not defined', 'control2cloud')) return
    !  if (fu_fails(fu_type(proc_ptr) == processor_scaling, &
    !                     & 'Emission processor is wrong type', 'control2cloud')) return
    !  p_emis => fu_emission_ptr(control)
    !  n_emis = fu_n_emission(control)
    !  emis_data => fu_data(proc_ptr)
    !  emis_data(1:n_emis) = max(p_emis(1:n_emis), 0.0)
    !  call msg('Emission to cloud, max = ', maxval(p_emis(1:n_emis)))
    !  call msg('Number of negative values ignored:', count(p_emis(1:n_emis) < 0))
    !end if ! have emission xy
    
    ifKalmanEnsemble = any((/flag_enkf, flag_enks/) == model%rules%daRules%method)
    
    if (fu_have_emission_zt(fu_mode(control))) then
      if (fu_fails(associated(proc_ptr), 'Emission processor not associated', 'control2cloud')) return
      if (fu_fails(defined(proc_ptr), 'Emission processor not defined', 'control2cloud')) return
      p_emis => fu_emission_ptr(control)
      n_emis = fu_n_emission(control)
      emis_data => fu_data(proc_ptr)
      !call msg('About to copy data...')
      if (fu_fails(associated(p_emis), 'p_emis not associated', 'control2cloud')) return
      if (fu_fails(associated(emis_data), 'emis_data not associated', 'control2cloud')) return
      emis_data(1:n_emis) = max(p_emis(1:n_emis), 0.0)
      call msg('Emission to cloud, max = ', maxval(p_emis(1:n_emis)))
      call msg('Number of negative values ignored:', count(p_emis(1:n_emis) < 0))
    end if ! have emission zt


    if (fu_have_emission_volc(fu_mode(control))) then
      p_emis => fu_emission_ptr(control)
      n_emis = fu_n_emission(control)
      if(ifKalmanEnsemble)then
        call store_perturbation_val(model%perturbations, perturb_emis_volc, p_emis, n_emis, &
                                & size(fu_species_emission_ctrl(control)))
      else
        call set_error('Volcano source is for EnKF/EnKS only but perturbations not associated', 'control2cloud')
      endif
      if(error)return
    end if     ! volcano source

    if (fu_have_emission_xy(fu_mode(control))) then
      if (fu_fails(associated(proc_ptr), 'Emission processor not associated', 'control2cloud')) return
      if (fu_fails(defined(proc_ptr), 'Emission processor not defined', 'control2cloud')) return
      if (fu_fails(fu_type(proc_ptr) == processor_scaling, &
                         & 'Emission processor is wrong type', 'control2cloud')) return

      p_emis => fu_emission_ptr(control)
      n_emis = fu_n_emission(control)
      emis_data => fu_data(proc_ptr)
      emis_data(1:n_emis) = max(p_emis(1:n_emis), 0.0)
      !
      ! This can be done via 4D-Var or EnKF but here it matters only because EnKF needs extra
      ! structure: emission perturbation. 
      !
!      if(model%rules%daRules%method == flag_enkf .or. model%rules%daRules%method == flag_enks)then
      if(ifKalmanEnsemble)then
        call store_perturbation_val(model%perturbations, perturb_emis_xy, p_emis, n_emis, &
                                  & size(fu_species_emission_ctrl(control)))
        if(error)return
      endif
    end if


!!$    case default
!!$      call set_error('Unsupported control mode', 'control2cloud')
!!$      return
!!$    end select
    

  end subroutine control2cloud

  !************************************************************************************

  subroutine set_emis_proc_by_ctrl(cloud, disprules, control_mode)
    implicit none
    type(silam_pollution_cloud), intent(inout) :: cloud
    integer, intent(in) :: control_mode
    type(general_dispersion_rules), intent(in) :: disprules

    type(Temission_processor), pointer :: proc_ptr
    character(len=*), parameter :: sub_name = 'set_emis_proc_by_ctrl' 
    type(silja_time), dimension(:,:), pointer :: time_slots
    integer :: stat
    

    allocate(proc_ptr, stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
    
    if (fu_have_emission_zt(control_mode)) then
      if (incompatible(control_mode, (/control_emis_xy, control_emis_volc/))) return
      call get_emis_time_slots(disprules%darules, time_slots)
      if (error) return
      call set_emission_processor(fu_emisMM_ptr(cloud), fu_species_transport(cloud), &
                                & disprules%darules%time_height_mode,  proc_ptr, time_slots)
      if (error) return
      deallocate(time_slots)
    else if (fu_have_emission_xy(control_mode)) then
      if (incompatible(control_mode, (/control_emis_zt, control_emis_volc/))) return
      call set_emission_processor(fu_emisMM_ptr(cloud), fu_species_transport(cloud), &
                                & processor_scaling,  proc_ptr)
      if (error) return
    else if (fu_have_emission_volc(control_mode)) then
      if (incompatible(control_mode, (/control_emis_zt, control_emis_xy/))) return
      allocate(time_slots(2,1), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
      time_slots(slot_start, 1) = disprules%startTime
      time_slots(slot_end, 1) = really_far_in_future
      call set_emission_processor(fu_emisMM_ptr(cloud), fu_species_transport(cloud), &
                                & processor_tm_hgt_force,  proc_ptr, time_slots)
      if (error) return
      deallocate(time_slots)
    else
      call set_emission_processor(fu_emisMM_ptr(cloud), fu_species_transport(cloud), &
                                & processor_void,  proc_ptr)
      if (error) return
    end if
    
    call set_emission_processor_ptr(cloud, proc_ptr)

  contains
    
    logical function incompatible(controlvar, incompatibles)  
      implicit none
      integer, intent(in) :: controlvar
      integer, dimension(:), intent(in) :: incompatibles
      
      integer :: ii
      
      do ii = 1, size(incompatibles)
        if (iand(controlvar, incompatibles(ii)) /= 0) then
          call set_error('Unsupported control configuration', sub_name)
          incompatible = .true.
          return
        end if
      end do
      incompatible = .false.
    end function incompatible
    
  end subroutine set_emis_proc_by_ctrl
  
  !************************************************************************************

  subroutine set_default_cnc(map_cnc)
    ! If the control does not include some (or any) initial state species, the defaults
    ! are read from the array set by this sub. In a single-window 4d-var this will be the
    ! initial condition as set in control file, for sequential 4d-var this will be the
    ! previous forecast field.
    implicit none
    type(TMass_map), intent(in) :: map_cnc

    integer :: nsp, nz, nx, ny
    integer :: stat

    character(len=*), parameter :: sub_name = 'set_default_cnc'
    
    if (fu_fails(map_cnc%nsrc == 1, 'Cannot handle multiple sources', sub_name)) return
    if (fu_fails(defined(map_cnc), 'Mass map not defined', sub_name)) return

    nsp = map_cnc%nSpecies
    nx = map_cnc%nx
    ny = map_cnc%ny
    nz = map_cnc%n3d
    if (.not. allocated(emis_only_init_state)) then
      allocate(emis_only_init_state(nsp,1,nz,nx,ny), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
    end if
    emis_only_init_state(:, 1, 1:nz, 1:nx, 1:ny) = map_cnc%arm(:, 1, 1:nz, 1:nx, 1:ny)

  end subroutine set_default_cnc

  
  !************************************************************************************
  
  subroutine move_neg_mass(mass_map, direction)
    ! a helper subroutine for control2cloud, control_from_cloud to handle negative values
    ! in PM assimilation.
    implicit none
    type(Tmass_map), intent(inout) :: mass_map
    character(len=*), intent(in) :: direction

    integer :: ind_species, ind_neg, ind_species_pos
    character(len=substNmLen) :: subst_name_neg, subst_name_pos
    type(silam_material), pointer :: material_pos
    type(silam_species) :: species_pos

    do ind_species = 1, size(mass_map%species)
      subst_name_neg = fu_substance_name(mass_map%species(ind_species))
      ind_neg = index(subst_name_neg, '__NEG__')
      if (ind_neg < 1) then
        cycle
      else if (ind_neg == 1) then
        call msg('Substance name: ' // trim(subst_name_neg))
        call set_error('Cannot handle substance name', 'move_neg_mass')
        return
      else if (ind_neg > 1) then
        subst_name_pos = subst_name_neg(1:ind_neg-1)
        material_pos => fu_get_material_ptr(subst_name_pos)
        if (error) return
        call set_species(species_pos, material_pos, fu_mode(mass_map%species(ind_species)))
        ind_species_pos = fu_index(species_pos, mass_map%species)
        if (ind_species_pos <= 0) then
          call set_error('Failed match found for' // trim(subst_name_neg), 'move_neg_mass')
          return
        end if
      end if

      if (debug_level > 0) then
        call msg('Total negative mass:', sum(mass_map%arm(ind_species_pos,:,:,:,:), &
             & mass_map%arm(ind_species_pos,:,:,:,:) < 0))
      end if

      select case (direction)
      case ('neg_pos')
        ! move from main species to __NEG__ if < 0
        call msg('Transfer negative mass, to, from', ind_species, ind_species_pos)
        where (mass_map%arm(ind_species_pos,:,:,:,:) < 0.0)
          ! no addition since the negative is zero or has been zeroed below.
          mass_map%arm(ind_species,:,:,:,:)  = -1*mass_map%arm(ind_species_pos,:,:,:,:)
        end where
      case ('pos_neg')
        ! add from __NEG__ species to main and multiply by -1
        call msg('Transfer negative mass, to, from', ind_species_pos, ind_species)
        mass_map%arm(ind_species_pos,:,:,:,:) &
             & = mass_map%arm(ind_species_pos,:,:,:,:) - mass_map%arm(ind_species,:,:,:,:)
        ! must zero negative species: after assimilation negative masses are put back but
        ! some cells might not be negative anymore.
        mass_map%arm(ind_species,:,:,:,:) = 0.0
      case default
        call set_error('Bad direction', 'move_neg_mass')
        return
      end select
    end do

  end subroutine move_neg_mass


  !************************************************************************************

  subroutine vector_from_model(control, model, handle_negatives, massmap)
    ! Intended to be replacement for gradient_from_cloud( and 
    !                                control_from_cloud
    implicit none
    ! Note: the gradient is not normalized!
    type(da_control), intent(inout) :: control
    type(model_container), intent(in) :: model
    logical, intent(in) :: handle_negatives ! control_from_cloud 
    type(Tmass_map), intent(in), target, optional :: massmap ! if Present -- do gradient
                                                  ! else -- take massmap from model%cloud and do control

    integer :: ix, iy, iz, isp, iTask, lenChunk,  nz,  ind, n_emis, isp_transp, nsp_contr
    integer :: nx, ny, offx, offy, gnx, gny, offglob, offloc, wrk_size
    integer :: nx1, ny1, offx1, offy1 
    integer, dimension(:), pointer :: iPtr
    type(Tmass_map), pointer :: map_c
    real, dimension(:), pointer :: p_init, dx_ptr, dy_ptr, p_emis, emis_data, wrk_local
    real, dimension(max_levels) :: thickness
    real :: volume, inv_scaling
    type(Temission_processor), pointer :: processor_ptr
    logical :: ifGradient, ifKalmanEnsemble
    character(len=*), parameter :: sub_name = 'vector_from_model'
    integer :: mode

    wrk_local => null() !! to be used only in MPI configuration


    mode = fu_mode(control)
    if (fu_fails( mode == model%rules%darules%controlVariable,"mode /= rules%controlVariable",  sub_name)) return
    if (present(massmap)) then !Do gradient
       map_c => massmap
       ifGradient = .True. 
       inv_scaling = 1./observation_scaling
    else
       ! do control vector
       if (fu_fails(fu_in_physical_space(control), 'Control not in physical space', sub_name)) return    
       map_c => fu_concMM_ptr(model%cloud)
       ifGradient = .False.
       inv_scaling = 1.
    endif

    !call msg('Sum map_c', sum(map_c%arm))
    call smpi_get_decomposition(nx, ny, offx, offy, gnx, gny)
    nz = map_c%n3d
    
    thickness(1:nz_dispersion) = (/(fu_layer_thickness_m(fu_level(dispersion_vertical, iz)),&
                                  & iz=1, nz_dispersion)/)

    if (iand(mode,control_init) /= 0) then ! Init present
      if (handle_negatives) call move_neg_mass(map_c, 'pos_neg')
      !call msg('Sum map_c 2', sum(map_c%arm(control%ind_species_init(1),:,:,:,:)))
      nsp_contr = fu_nsp_initial(control)
      iPtr => fu_isp_initial(control)
      lenChunk = nz*nsp_contr
      call smpi_get_maxsubdomain_size(nx1, ny1) !!  reuse of vars
      wrk_size = nx1*ny1*nz*nsp_contr + 4
      !call msg('nsp_contr', nsp_contr)
      !call msg('isp_contr', control%ind_species_init(1))

      if (smpi_adv_tasks == 1) then
        p_init => fu_initial_ptr(control)
      else
        wrk_local => fu_work_array(wrk_size) !
        p_init => wrk_local(5:wrk_size)  !keep first indices for geometry
        wrk_local(1:4) = (/nx, ny, offx, offy/)
      endif
      do iy = 1, ny
        do ix = 1, nx
          do iz = 1, nz
            volume = thickness(iz) * fu_cell_size(dispersion_grid, ix, iy)
            do isp = 1, nsp_contr
              isp_transp = iPtr(isp)
              ind = (((iy-1)*nx + ix-1)*nz + iz-1) *nsp_contr + isp
              if (ifGradient) then
                p_init(ind)  = (map_c%arm(isp_transp, DA_POSITIVE, iz, ix, iy) &
                            & + DA_ZERO_COEF_GRAD * map_c%arm(isp_transp, DA_ZERO, iz, ix, iy) &
                            & + DA_NEGT_COEF_GRAD * map_c%arm(isp_transp, DA_NEGATIVE, iz, ix, iy)) &
                          & * inv_scaling !/ observation_scaling
              else
                p_init(ind) = map_c%arm(isp_transp, 1, iz, ix, iy) / volume

              end if
            end do
          end do
        end do
      end do
      !call msg("p_init(1:1000)", p_init(1:1000))
      
      if ((.not. ifGradient) .and. model%rules%darules%use_log_cloud) then
         p_init(:) = log(p_init(:) + 0.02)
      endif

      if (smpi_adv_tasks > 1) then !! Do exchange
        if(smpi_adv_rank == 0) then   !! Gather the vector
          p_init => fu_initial_ptr(control)   
          do iTask = 0,smpi_adv_tasks-1 !Receive vector -- only from other tasks
            if (error) cycle
            if (itask /=0)  call smpi_recv(wrk_local, wrk_size, iTask) 
            if (error) cycle
            nx1 = wrk_local(1); ny1=wrk_local(2);  offx1= wrk_local(3);  offy1 = wrk_local(4);
            do iy = 1, ny1
               do ix = 1, nx1
                  offglob = ((offy1 + iy-1)*gnx + (offx1+ix-1)) *nz*nsp_contr
                  offloc  = ((iy-1)*nx1 + ix-1) * nz*nsp_contr 
                  p_init(offglob+1 : offglob+lenChunk) = wrk_local(offloc+5:offloc+4+lenChunk)
               enddo
            enddo
          enddo
        else
          !not Master -- just send our stuff
          call smpi_send(wrk_local, wrk_size, 0)
        endif
        call free_work_array(wrk_local) 
      endif !do exchange
    end if  !Inital state present
    
    ifKalmanEnsemble = any((/flag_enkf, flag_enks/) == model%rules%daRules%method)
   
    if (fu_have_emission_zt(fu_mode(control))) then
      p_emis => fu_emission_ptr(control)
      if (fu_fails(associated(fu_emission_processor_ptr(model%cloud)), 'proc not associated', sub_name)) return
      emis_data => fu_data(fu_emission_processor_ptr(model%cloud))
      if (fu_fails(size(p_emis) == size(emis_data), 'sizes don''t match', sub_name)) return
      p_emis(:) = emis_data(:) * inv_scaling !
    end if

    if (fu_have_emission_xy(fu_mode(control))) then
      p_emis => fu_emission_ptr(control)
      n_emis = fu_n_emission(control)
      if(ifKalmanEnsemble)then
        call retrieve_perturbation_val(model%perturbations, perturb_emis_xy, p_emis, n_emis, int_missing)
        if(error)return
      endif
      
!      if (fu_fails(associated(fu_emission_processor_ptr(model%cloud)), 'proc not associated', sub_name)) return
!      ind_pert_xy = fu_index(perturb_emis_xy, model%perturbations(:)%target)
!      p_emis(1:n_emis) = model%perturbations(ind_pert_xy)%storage
    end if


    if (iand(mode,control_emis_volc) /= 0) then 
      n_emis = fu_n_emission(control)
      p_emis => fu_emission_ptr(control)
      if(ifKalmanEnsemble)then
        call retrieve_perturbation_val(model%perturbations, perturb_emis_volc, p_emis, n_emis, &
                                     & size(fu_species_emission_ctrl(control)))
      else
        call set_error('Volcano source is for EnKF/EnKS only but perturbations not associated', 'vector_from_model')
      endif
      if(error)return

!      
!      if (fu_fails(associated(model%perturbations), 'perturbations not associated', sub_name)) return
 !     if (fu_fails(size(model%perturbations) > 0, 'perturbations is empty', sub_name)) return
 !     ind_pert_volc = fu_index(perturb_emis_volc, model%perturbations(:)%target)
 !     if (fu_fails(ind_pert_volc > 0, 'no perturb_emis_volc', sub_name)) then !return
 !       call msg('perturbations:', model%perturbations(1)%target)
 !       call msg('...')
 !       return
 !     end if
 !     nsp_contr = size(fu_species_emission_ctrl(control)) 
 !     if (fu_fails(n_emis == nsp_contr*2 + num_params_volc, 'bad n_emis', sub_name)) return
 !     p_emis(1:nsp_contr) = model%perturbations(ind_pert_volc)%scale
 !     p_emis(nsp_contr+1:2*nsp_contr) = model%perturbations(ind_pert_volc)%offset
 !     p_emis(2*nsp_contr+1:n_emis) = model%perturbations(ind_pert_volc)%storage
    end if
      
  end subroutine vector_from_model



  !************************************************************************************

  subroutine model_forward(model, obs_ptr, simrules, make_output)
    !
    ! After applying the control, run the model.
    !
    implicit none
    type(model_container) :: model
    type(ObservationPointers), intent(inout) :: obs_ptr
    type(general_dispersion_rules), intent(inout) :: simrules
    logical, intent(in), optional :: make_output

    type(silam_source), pointer :: em_source
    type(Tfield_buffer), pointer :: disp_buf_ptr, meteo_ptr, output_buf_ptr
    type(silja_shopping_list), pointer :: &
         & pMeteo_input_dyn_shopping_list, &
         & pMeteo_full_dyn_shopping_list, &
         & pDisp_dyn_shopping_list, &
         & pDisp_stat_shopping_list
    type(silam_output_definition), pointer :: outDef ! might be removed eventually
    type(silja_wdr), pointer :: wdr
    type(mini_market_of_stacks), pointer :: meteoMarketPtr, dispersionMarketPtr, BCMarketPtr, &
                                          & srcMarketPtr, outputMarketPtr
    real :: weight_past
    character(len=clen) :: command_string
    type(Tmeteo_input), pointer :: meteo_input_ptr
!    type(general_dispersion_rules) :: rules_fwd
    !type(Tmass_map), pointer :: map_px, map_py, map_pz, map_sl
    type(silja_time) :: now
    type(silam_pollution_cloud), pointer :: cloud
    logical :: ifWriteProgressFileTmp, make_output_tmp
    logical, save :: if_first = .true.
    type(t_perturbation), dimension(:), pointer :: perturbations


    ! Unpack the model structures:
    !
    em_source => model%source

    pMeteo_input_dyn_shopping_list => model%meteo_input_dyn_shopping_list
    pMeteo_full_dyn_shopping_list => model%meteo_full_dyn_shopping_list
    pDisp_dyn_shopping_list => model%disp_dyn_shopping_list
    pDisp_stat_shopping_list => model%disp_stat_shopping_list

    outDef => model%out_def
    meteo_ptr => model%meteo_ptr
    disp_buf_ptr => model%disp_buf_ptr
    output_buf_ptr => model%output_buf_ptr
    meteoMarketPtr => model%meteo_market_ptr
    dispersionMarketPtr => model%disp_market_ptr
    outputMarketPtr => model%output_market_ptr
    BCMarketPtr => model%bc_market_ptr
    wdr => model%wdr
    meteo_input_ptr => model%meteo_input
    cloud => model%cloud
    perturbations => model%perturbations

!
! simrules cannot be copied now! run_dispersion will modify them, in particular,
! by setting the low_cnc and low_mass thresholds.
!
!    rules_fwd = simrules
!    !rules_fwd%startTime = simrules%darules%da_begin
!    !rules_fwd%periodToCompute = simrules%darules%da_window
!    rules_fwd%ifWriteProgressFile = .false.
!    now = rules_fwd%startTime

    now = simrules%startTime
    ifWriteProgressFileTmp = simrules%ifWriteProgressFile  ! store temporary
    make_output_tmp = simrules%if_make_output
    if (da_test) then
      simrules%if_make_output = .true.
    else
      if (present(make_output)) then
        simrules%if_make_output = make_output
      else
        simrules%if_make_output = .false.
      end if
    end if
    
    call start_count('forward')
!    call run_dispersion(cloud, em_source, wdr, rules_fwd, &

    ! below seems no longer needed?
!!$    if (.not. if_first) then
!!$      ! If several forward runs are made after each other, the buffer times are wrong
!!$      ! unless set properly.
!!$      call arrange_buffer(dispersionMarketPtr, met_src_missing, now, simrules%timestep, & 
!!$                        & disp_buf_ptr,  dispersion_verticalPtr)
!!$      if (error) return
!!$    end if
      
    call run_dispersion(cloud, em_source, wdr, simrules, &
                      & pMeteo_input_dyn_shopping_list, pMeteo_full_dyn_shopping_list, &
                      & pDisp_dyn_shopping_list, pDisp_stat_shopping_list, &
                      & outDef, meteo_ptr, disp_buf_ptr, obs_ptr, meteo_input_ptr, output_buf_ptr, &
                      & meteoMarketPtr, dispersionMarketPtr, BCMarketPtr, outputMarketPtr, &
                      & model%tla_traj, perturbations)
    call stop_count('forward')

    simrules%ifWriteProgressFile = ifWriteProgressFileTmp
    simrules%if_make_output = make_output_tmp

    if_first = .false.
  end subroutine model_forward

  !************************************************************************************

  subroutine model_adjoint(model, obs_ptr, simrules)
    ! Run the model in adjoint mode - the H* operator will be active.
    !
    implicit none
    type(model_container), intent(inout) :: model
    type(ObservationPointers), intent(inout) :: obs_ptr
    type(general_dispersion_rules), intent(inout) :: simrules

    type(silam_source), pointer :: em_source
    type(Tfield_buffer), pointer :: disp_buf_ptr, meteo_ptr, output_buf_ptr
    type(silja_shopping_list), pointer :: &
         & pMeteo_input_dyn_shopping_list, &
         & pMeteo_full_dyn_shopping_list, &
         & pDisp_dyn_shopping_list, &
         & pDisp_stat_shopping_list
    type(silam_output_definition), pointer :: outDef ! might be removed eventually
    type(silja_wdr), pointer :: wdr
    type(mini_market_of_stacks), pointer :: meteoMarketPtr, dispersionMarketPtr, BCMarketPtr, &
                                          & outputMarketPtr
    type(silam_pollution_cloud), pointer :: cloud
    type(Tmass_map), pointer :: map_px, map_py, map_pz
    character(len=clen) :: command_string
    type(Tmeteo_input), pointer :: meteo_input_ptr
    !type(general_dispersion_rules) :: rules_adj
    type(silja_time) :: now, startTimeTmp
    type(silja_interval) :: timestepTmp, periodToComputeTmp
    logical :: ifWriteProgressFileTmp, if_invert_substepsTmp, if_make_output_tmp
    type(Tini_boundary_rules) :: iniBoundaryRulesTmp
    
    ! Unpack the model structures:
    !
    em_source => model%source

    pMeteo_input_dyn_shopping_list => model%meteo_input_dyn_shopping_list
    pMeteo_full_dyn_shopping_list => model%meteo_full_dyn_shopping_list
    pDisp_dyn_shopping_list => model%disp_dyn_shopping_list
    pDisp_stat_shopping_list => model%disp_stat_shopping_list

    outDef => model%out_def
    meteo_ptr => model%meteo_ptr
    disp_buf_ptr => model%disp_buf_ptr
    output_buf_ptr => model%output_buf_ptr
    meteoMarketPtr => model%meteo_market_ptr
    dispersionMarketPtr => model%disp_market_ptr
    outputMarketPtr => model%output_market_ptr
    BCMarketPtr => model%bc_market_ptr
    wdr => model%wdr
    meteo_input_ptr => model%meteo_input
    cloud => model%cloud

    !rules_adj = simrules
    !rules_adj%startTime = simrules%darules%da_begin + simrules%darules%da_window
    !rules_adj%periodToCompute = fu_opposite(simrules%darules%da_window)
    !rules_adj%timestep = fu_opposite(simrules%timestep)
    !rules_adj%ifWriteProgressFile = .false.
    !rules_adj%if_invert_substeps = .true.
    !rules_adj%if_make_output = .false.

!    call set_missing(rules_adj%iniBoundaryRules)

    !
    ! Store temporaries
    !
    startTimeTmp = simrules%startTime
    periodToComputeTmp = simrules%periodToCompute
    timestepTmp = simrules%timestep
    ifWriteProgressFileTmp = simrules%ifWriteProgressFile
    iniBoundaryRulesTmp = simrules%iniBoundaryRules
    if_invert_substepsTmp = simrules%if_invert_substeps
    if_make_output_tmp = simrules%if_make_output
    !
    ! Reset values for adjoint run
    !
    simrules%startTime = simrules%darules%da_begin + simrules%darules%da_window
    simrules%periodToCompute = fu_opposite(simrules%darules%da_window)
    simrules%timestep = fu_opposite(simrules%timestep)
    simrules%ifWriteProgressFile = .false.
    call set_missing(simrules%iniBoundaryRules)
    simrules%if_invert_substeps = .true.
    simrules%if_make_output = .false.
    
    call start_count('adjoint')
    call msg('')
    call msg('Running adjoint')
    
    call run_dispersion(cloud, em_source, wdr, simrules, &
                      & pMeteo_input_dyn_shopping_list, pMeteo_full_dyn_shopping_list, &
                      & pDisp_dyn_shopping_list, pDisp_stat_shopping_list, &
                      & outDef, meteo_ptr, disp_buf_ptr, obs_ptr, meteo_input_ptr, output_buf_ptr, &
                      & meteoMarketPtr, dispersionMarketPtr, BCMarketPtr, outputMarketPtr, &
                      & model%tla_traj)
    call stop_count('adjoint')

    !
    ! Return temporaries
    !
    simrules%startTime = startTimeTmp
    simrules%periodToCompute = periodToComputeTmp
    simrules%timestep = timestepTmp
    simrules%ifWriteProgressFile = ifWriteProgressFileTmp
    simrules%iniBoundaryRules = iniBoundaryRulesTmp
    simrules%if_invert_substeps = if_invert_substepsTmp
    simrules%if_make_output = if_make_output_tmp
    
  end subroutine model_adjoint

  !************************************************************************************

  subroutine run_forecast(control, model)
    ! Run the forecast (after finishing assimilation): output will be produced. Emission
    ! processors will be active, if defined.
    !
    implicit none
    type(model_container), intent(inout) :: model
    type(da_control), intent(in) :: control

    if (.not. fu_in_physical_space(control)) then
      call set_error('Control variable not in physical space', 'get_obs_cost')
      return
    end if

    if (model%obs_ptr%hasObservations) call reset_all(model%obs_ptr, .true.)
    call control2cloud(control, model, handle_negatives=fu_mode(control) == da_initial_state)
    call model_forward(model, model%obs_ptr, model%rules, make_output=.true.) 
    
  end subroutine run_forecast

  !************************************************************************************

  subroutine forward_3d(model)
    ! Run the 3D-Var forward model: only observations. The assimilation window is passed
    ! to observeAll - observations inside this interval will contribute.
    implicit none
    type(model_container) :: model

    type(silja_time) :: assim_begin
    type(silja_interval) :: assim_window

    assim_window = model%rules%darules%da_window
    assim_begin = model%rules%darules%da_begin - assim_window*0.5

    if (debug_level > 0) then
      call collect_total_masses(model%cloud)
      call report_total_masses(model%cloud, 1, .false.)
    end if
    call msg('FWD_3D: timestep' // trim(fu_str(assim_window)))
    call observeAll(model%obs_ptr, model%cloud,&
                  & model%meteo_ptr, model%disp_buf_ptr, model%rules%chemicalRules, &
                  & assim_window, assim_begin)

  end subroutine forward_3d

  !************************************************************************************

  subroutine adjoint_3d(model)
    ! Run the 3D-Var adjoint model: only adjoint observation operators. The assimilation
    ! window is passed to observeAll - observations inside this interval will contribute.
    implicit none
    type(model_container) :: model

    type(silja_time) :: assim_begin
    type(silja_interval) :: assim_window

    assim_window = model%rules%darules%da_window
    assim_begin = model%rules%darules%da_begin + assim_window*0.5

    ! so injection interval is now assim_begin-assim_window...assim_begin and assim_window
    ! < 0 as required.
    call injectAll(model%obs_ptr, model%cloud,&
                  & model%meteo_ptr, model%rules%chemicalRules, model%rules%dynamicsRules, &
                  & fu_opposite(assim_window), assim_begin, model%mass_map_adj) 
  end subroutine adjoint_3d

  !************************************************************************************

  subroutine obs_from_mdl(control, model, mdl_obs_values)
    !
    ! The composite action of the model and observation operators: apply control, run the
    ! model and collect the observations into a vector.
    !
    implicit none
    type(model_container), intent(inout) :: model
    type(da_control), intent(in) :: control
    real, dimension(:), intent(out) :: mdl_obs_values

    integer ::n_obs_values

    if (.not. fu_in_physical_space(control)) then
      call set_error('Control variable not in physical space', 'get_obs_cost')
      return
    end if
    
    call reset_all(model%obs_ptr, .true.)
    call control2cloud(control, model, handle_negatives=.false.)
    if (model%rules%darules%method == flag_3dvar) then
        call forward_3d(model)  ! Observe every model datum
    else if (model%rules%darules%method == flag_4dvar .or. model%rules%darules%method == flag_h_matrix) then
      call model_forward(model, model%obs_ptr, model%rules)
    else
      call set_error('Bad assimilation method', 'obs_from_mdl')
      return
    end if
    if (error) return
    call collect_model_data(model%obs_ptr, mdl_obs_values, n_obs_values)
    
  end subroutine obs_from_mdl

  !************************************************************************************

  subroutine get_obs_grad(model, gradient)
    ! The combined action of H* and M* on the y - Hx residual. The gradient is collected
    ! into a control variable object.
    implicit none
    type(model_container), intent(inout) :: model
    !type(da_control), intent(in) :: control
    type(da_control), intent(inout) :: gradient
    ! use this mass map for adjoint injection in adjoint_3d (needed for 3dvar)

    type(Tmass_map), pointer :: map_c
    type(Tmass_map), pointer :: map_px, map_py, map_pz
    logical, save :: processor_is_set = .false. 
    type(Temission_processor), pointer :: processor_ptr
    integer :: ind_spec_anl, ind_spec_transp

    !if (fu_fails(associated(fu_emission_processor_ptr(model%cloud)), &
    !           & 'Emis proc not associated', 'get_obs_grad')) return
    !if (fu_fails(defined(fu_emission_processor_ptr(model%cloud)), &
    !           & 'Emis proc not defined', 'get_obs_grad')) return
!!$    if (fu_mode(control) == DA_EMISSION_AND_INITIAL .or. fu_mode(control) == DA_EMISSION_CORRECTION &
!!$      & .and. .not. processor_is_set) then
!!$      call set_error('No emis proc', 'get_obs_grd') 
!!$      return
!!$      call set_emis_proc_in_cloud(model)
!!$      processor_is_set = .true.
!!$    end if

    processor_ptr => fu_emission_processor_ptr(model%cloud)
    if (defined(processor_ptr)) call zero_processor(processor_ptr)
    if (debug_level > 0) then
      call msg("In the beginning of get_obs_grad")
      call collect_total_masses(model%cloud)
      call report_total_masses(model%cloud, 1, .false.)
    end if
    if (model%rules%darules%method == flag_4dvar) then
      call reset_cloud_data_structures(model%cloud, ifMomentsToo=.true.)
      call model_adjoint(model, model%obs_ptr, model%rules)
     !      call gradient_from_cloud(gradient, model%cloud, model%rules%darules)
      call vector_from_model(gradient, model, .false., fu_concMM_ptr(model%cloud))

    else
      if (fu_fails(associated(model%mass_map_adj), 'mass_map_adj not associated', 'get_obs_grad')) return
      ! the observations might change values also for non-control species. But these are
      ! not used ever, so no need to bother.
      call reset_control_species(model%mass_map_adj, gradient)
      call adjoint_3d(model)
      !call gradient_from_cloud(gradient, model%cloud, model%rules%darules, model%mass_map_adj)
      call vector_from_model(gradient, model, .false., model%mass_map_adj)
    end if
    if (error) return
    if (debug_level > 1) call report_total_masses(model%cloud, 1, .false.)
    if (error) return

  contains

    subroutine set_emis_proc_in_cloud(model)
      implicit none
      type(model_container), intent(inout) :: model

      integer :: stat

      allocate(processor_ptr, stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', 'set_emission_proc_in_cloud')) return
      call set_emission_processor(fu_emisMM_ptr(model%cloud), &
                                & fu_species_transport(model%cloud), &
                                & processor_scaling,  processor_ptr)
      call set_emission_processor_ptr(model%cloud, processor_ptr)

    end subroutine set_emis_proc_in_cloud

    subroutine reset_control_species(mass_map, control)
      ! set control species to 0, leave all others untouched.
      implicit none
      type(Tmass_map), intent(inout) :: mass_map
      type(da_control), intent(in) :: control

      integer :: isp_contr
      integer, dimension (:),pointer :: iPtr
      
      iPtr => fu_isp_initial(control)
      do isp_contr = 1, fu_nsp_initial(control)
        mass_map%arm(iPtr(isp_contr), :,:,:,:) = 0.0
      end do
    end subroutine reset_control_species

  end subroutine get_obs_grad


  subroutine set_background(background, cloud, rules)
    ! Set the background value for a control variable according to the DA rules.
    ! 
    implicit none
    type(da_control), intent(inout) :: background
    type(silam_pollution_cloud), intent(in) :: cloud
    type(da_rules), intent(in) :: rules

    type(Tmass_map_ptr), dimension(1) :: mass_map_ptr_array
    type(Tmass_map), pointer :: map_c, map_cnc
    integer :: ix, iy, iz, isp, nx, ny, nz, nsp, idim, ind, isp_transp, nsp_contr
    real, dimension(:), pointer :: init_p
    real :: init_bgr_def, emis_bgr_def

    init_bgr_def = 0.0
    if (fu_mode(background) /= DA_INITIAL_STATE) then
      if (rules%use_zero_emis_backgr) then
        emis_bgr_def = 0.0
        call msg_warning('Will use 0.0 as background emission correction')
      else
        emis_bgr_def = 1.0
        call msg_warning('Will use 1.0 as background emission correction')
      end if
    else
      emis_bgr_def = real_missing
    end if
    call control_from_files(background, cloud, &
                          & rules%initialStateBgrFile, rules%have_init_background, init_bgr_def, &
                          & rules%emissionCorrectionBgrFile, rules%have_emis_background, emis_bgr_def, &
                          & rules%da_begin, rules%ifRandomise)
    if (error) return

  end subroutine set_background

  
  
  !!************************************************************************************
  !
  !subroutine init_da(nlptr, simrules, rules, obs_pointers, cloud)
  !  ! 
  !  ! Initialize the DA modules. At this point everything else needs to be defined,
  !  ! including grids, verticals and chemical setup.
  !  !
  !  implicit none
  !  type(Tsilam_namelist), pointer :: nlPtr
  !  type(general_dispersion_rules), intent(in) :: simrules
  !  type(DA_rules), intent(in) :: rules
  !  type(observationPointers), intent(out) :: obs_pointers
  !  type(silam_pollution_cloud), target :: cloud
  !  
  !  select case(fu_content(nlptr, 'method'))
  !  case ('4D')
  !    call set_observations(rules, fu_species_transport(cloud), fu_species_optical(cloud), obs_pointers)
  !  case ('3D')
  !    !call init_3dvar(dispersion_grid, dispersion_vertical)
  !  case ('')
  !    call set_error('Missing da_method', 'init_da')
  !  case default
  !    call set_error('Strange da_method', 'init_da')
  !  end select
  !  
  !end subroutine init_da

  !*********************************************************************************************

  function fu_DA_time_window(rules) result(DA_time_window)
    !
    ! Returns the time window processed by the data assimilation. Should no assimilation needed
    ! the window is intrval_missing
    !
    implicit none
    
    ! Return value
    type(silja_interval) :: DA_time_window
    
    ! Imported parameters
    type(DA_rules), intent(in) :: rules
    
    if(rules%defined)then
      DA_time_window = rules%DA_window
      if(error)then
        !call unset_error('fu_DA_time_window')
        ! what, this can't be allowed...
        DA_time_window = interval_missing
      endif
    else
      DA_time_window = interval_missing
    endif
  end function fu_DA_time_window

  !************************************************************************************

  integer function fu_initial_fill_mode(rules) result(fill)
    ! possibilities:
    ! fill_none : only touch analysis species.
    ! fill_non_analysis : set non-analysis species from initial in forward, zero all in adjoint
    ! fill_all : fill everything from initial, zero all in adjoint
    implicit none
    type(da_rules), intent(in) :: rules
    logical :: have_initial

    if (fu_fails(rules%defined, 'DA rules not defined', 'fu_initial_fill_mode')) return

    have_initial = fu_have_initial(rules%controlVariable)

    if (rules%method == flag_4dvar .or. rules%method == flag_4dvar_seq) then
      if (have_initial) then
        if (fu_true(rules%have_analysis_species)) then
          ! only some species included in control
          fill = fill_non_analysis
        else
          ! control completely determines initial condition.
          fill = reset_only
        end if
      else
        ! emission-only control
        fill = fill_all
      end if
      
    !else if (rules%method == flag_enkf) then
    !  fill = fill_non_analysis
    else ! 3dvar or undefined
      fill = fill_none
    end if
    
  end function fu_initial_fill_mode



  !**********************************************************************************
  !
  ! Test codes
  !
  !**********************************************************************************


  subroutine adj_test_obs(model)
    implicit none
    type(model_container) :: model
    type(Tmass_map), pointer :: map_c
    integer :: nx, ny, nz, nsp
    real, dimension(:,:,:,:), allocatable :: arr1, arr2, m_arr

    integer :: n_obs_val
    real, dimension(:), pointer :: mdl_data, obs_data, obs_var
    type(silja_time) :: now, in_past, in_future
    type(silja_interval) :: timestep, long_interval
    type(Tmass_map), pointer :: MMptr

    now = model%rules%startTime
    timestep = model%rules%timestep
    model%rules%periodToCompute = timestep
    model%rules%darules%da_window = timestep
    !model%rules%darules%da_end = model%rules%startTime + (timestep)
    model%obs_ptr%hasObservations = .false.

    call msg('Will run model for one timestep')
    call model_forward(model, model%obs_ptr, model%rules)
    call msg('')
    call msg('done')
    call msg('')

    call msg('Now testing adjoint...')

    mdl_data => fu_work_array()
    obs_data => fu_work_array()
    obs_var => fu_work_array()

    map_c => fu_concMM_ptr(model%cloud)
    nx = map_c%nx
    ny = map_c%ny
    nz = map_c%n3d
    nsp = map_c%nspecies

    call reset_map_to_val(map_c, 0.0)
    MMptr => fu_advection_moment_X_MM_ptr(model%cloud)
    call reset_map_to_val(MMptr, 0.0)
    MMptr => fu_advection_moment_Y_MM_ptr(model%cloud)
    call reset_map_to_val(MMptr, 0.0)
    MMptr => fu_advection_moment_Z_MM_ptr(model%cloud)
    call reset_map_to_val(MMptr, 0.0)

    allocate(arr1(nsp,nz,nx,ny), arr2(nsp,nz,nx,ny), m_arr(nsp,nz,nx,ny))
    in_past = fu_set_time_utc(2000, 1, 1, 0, 0, 0.0)
    in_future = fu_set_time_utc(2020, 1, 1, 0, 0, 0.0)
    long_interval = in_future - in_past

    call random_number(arr1)
    call random_number(arr2)
    arr1 = (arr1+10)!*1000
    arr2 = (arr2+10)!*1000
    call reset_all(model%obs_ptr, .true.)
    map_c%arm = 0.0
    map_c%arm(1:nsp,1,1:nz,1:nx,1:ny) = arr1
    !call by_volume(map_c, 'mult')
    
    ! include timing
    !call observeAll(model%obs_ptr, model%cloud, model%meteo_ptr, model%disp_buf_ptr, model%rules%chemicalRules, &
    !              & model%rules%timestep, model%rules%startTime)
    
    ! ignore timing
    call observeAll(model%obs_ptr, model%cloud, model%meteo_ptr, model%disp_buf_ptr, model%rules%chemicalRules, &
                  & long_interval, in_past)

    call collect_model_data(model%obs_ptr, mdl_data, n_obs_val)
    call collect_obs_data(model%obs_ptr,  obs_data, n_obs_val)
    call collect_variance(model%obs_ptr, obs_var, n_obs_val)
!    call reset_map_to_val(fu_concMM_ptr(model%cloud), 0.0)
    MMptr => fu_concMM_ptr(model%cloud)
    call reset_map_to_val(MMptr, 0.0)
    !call reset_all(model%obs_ptr)
    !call injectAll(model%obs_ptr, model%cloud, model%meteo_ptr, &
    !             & model%rules%chemicalRules, model%rules%dynamicsRules, &
    !             & fu_opposite(model%rules%timestep), model%rules%startTime)

    ! ignore timing
    call injectAll(model%obs_ptr, model%cloud, model%meteo_ptr, &
                 & model%rules%chemicalRules, model%rules%dynamicsrules, &
                 & fu_opposite(long_interval), in_future)

    call by_volume(map_c, 'div')
    m_arr = map_c%arm(1:nsp,1,1:nz,1:nx,1:ny) / observation_scaling
    

    call msg('Mean model', sum(mdl_data(1:n_obs_val)/n_obs_val))
    call msg('Mean obs', sum(obs_data(1:n_obs_val)/n_obs_val))
    
    call msg('LEFT', sum(arr1*m_arr))
    call msg('RIGHT', sum(mdl_data(1:n_obs_val) &
                        & * ((mdl_data(1:n_obs_val)-obs_data(1:n_obs_val)) / obs_var(1:n_obs_val))))
    stop

  contains

    subroutine by_volume(mass_map, what)
      implicit none
      type(Tmass_map), intent(inout) :: mass_map
      character(len=*), intent(in) :: what

      integer :: ilev, ix, iy, nlevs
      real :: volume
      real, dimension(:), pointer :: thickness
      
      
      nlevs = mass_map%n3d
      
      thickness => fu_work_array()
      thickness(1:nlevs) = (/(fu_layer_thickness_m(fu_level(mass_map%vertTemplate, ilev)), ilev=1, nlevs)/)
      
      if (what == 'div') then
        do iy = 1, mass_map%ny
          do ix = 1, mass_map%nx
            do ilev = 1, mass_map%n3d
              volume = thickness(ilev) * fu_cell_size(mass_map%gridTemplate, ix, iy)
              mass_map%arm(:,:,ilev,ix,iy) = mass_map%arm(:,:,ilev,ix,iy) / volume 
            end do
          end do
        end do
      else
        do iy = 1, mass_map%ny
          do ix = 1, mass_map%nx
            do ilev = 1, mass_map%n3d
              volume = thickness(ilev) * fu_cell_size(mass_map%gridTemplate, ix, iy)
              mass_map%arm(:,:,ilev,ix,iy) = mass_map%arm(:,:,ilev,ix,iy) * volume 
            end do
          end do
        end do
      end if
      call free_work_array(thickness)
      
    end subroutine by_volume
  end subroutine adj_test_obs

  subroutine adj_test(model)
    implicit none
    type(model_container) :: model

    integer :: nx, ny, nz, nsp, istep
    type(tmass_map), pointer :: map_c
    real, dimension(:,:,:,:), allocatable :: arr1, arr2, m_arr

    type(silja_time) :: now
    type(silja_interval) :: timestep
    integer, parameter :: nsteps = 48
    type(Temission_processor), pointer :: proc_ptr
    integer :: ii
    type(Tmoment_mapping) :: cm_to_moment
    type(Tmass_map), pointer :: MMptr

    call create_moment_2_cm_mapping(fu_advection_moment_X_MM_ptr(model%cloud),&
                                  & fu_concMM_ptr(model%cloud), &
                                  & fu_low_mass_threshold(model%rules%chemicalRules), &
                                  & centre_mass_vs_moment_species, cm_to_moment)

    call random_seed(put=(/(12, ii=1,120)/))
    
    timestep = model%rules%timestep
    model%rules%periodToCompute = timestep
    model%rules%darules%da_window = timestep
    !model%rules%darules%da_end = model%rules%startTime + (timestep)
    model%obs_ptr%hasObservations = .false.

    ! Kill emissions
    !
    proc_ptr => fu_emission_processor_ptr(model%cloud)
    call set_emission_processor(fu_emisMM_ptr(model%cloud), fu_species_transport(model%cloud), &
                              & processor_scaling, proc_ptr)
    call zero_processor(proc_ptr)

    call msg('Will run model for one timestep')
    call model_forward(model, model%obs_ptr, model%rules)
    call msg('')
    call msg('done')
    call msg('')

    call msg('Now testing adjoint...')

    map_c => fu_concMM_ptr(model%cloud)
    nx = map_c%nx
    ny = map_c%ny
    nz = map_c%n3d
    nsp = map_c%nspecies

    MMptr => fu_concMM_ptr(model%cloud)
    call reset_map_to_val(MMptr, 0.0)
    MMptr => fu_advection_moment_X_MM_ptr(model%cloud)
    call reset_map_to_val(MMptr, 0.0)
    MMptr => fu_advection_moment_Y_MM_ptr(model%cloud)
    call reset_map_to_val(MMptr, 0.0)
    MMptr => fu_advection_moment_Z_MM_ptr(model%cloud)
    call reset_map_to_val(MMptr, 0.0)

    allocate(arr1(nsp,nz,nx,ny), arr2(nsp,nz,nx,ny), m_arr(nsp,nz,nx,ny))
    
    arr1 = 0.0
    arr2 = 0.0
    m_arr = 0.0
    
    !call random_number(arr1(:,:,5:nx-5,5:ny-5))
    !call random_number(arr1(:,:,:,:))
    !call random_number(arr2)
    !call random_number(arr2(:,:,5:nx-5,5:ny-5))
    !arr1 = arr1 + 1
    !arr2 = arr2 + 1

    !call random_number(arr1(:,:,20,20))
    !call random_number(arr2(:,:,20,20))

!!$    arr1(:,:,20,20) = 1.0
!!$    arr2(:,:,20,20) = 1.0
!!$
!!$    arr1(:,:,10:20,10:20) = 1.0
!!$    arr2(:,:,10:20,10:20) = 1.0

    call random_number(arr1(:,:,10:20,10:20))
    call random_number(arr2(:,:,10:20,10:20))
    arr1 = (arr1+10)*1000
    arr2 = (arr2+10)*1000



    !call random_number(arr1(:,:,:,:))
    !call random_number(arr2(:,:,:,:))


    !arr1 = (arr1 + 10)*1000
    !arr2 = (arr2 + 10)*1000

    !arr1(:,3,:,:) = 1.0
    !arr2(:,3,:,:) = 1.0

    model%rules%periodToCompute = timestep
    model%rules%darules%da_window = timestep*nsteps
    !model%rules%darules%da_end = model%rules%startTime + timestep*nsteps
    !model%rules%if_collect_linearization = .true.
    !allocate(linearizations(nsteps))
    !do istep = 1, nsteps
    !  call alloc_transp_t_lin(linearizations(istep), dispersion_grid, dispersion_vertical, &
    !                        & fu_nbr_of_species_transport(model%cloud), 1)
    !end do
    map_c%arm = 0.0
    map_c%arm(1:nsp,1,1:nz,1:nx,1:ny) = arr1
    timestep = timestep
    ! Now run the model process on arr1, forward

!!$    do istep = 1, nsteps
!!$      call msg('')
!!$      call msg('', istep)
!!$      call msg('')

!!$      CALL advect_pollution_cloud_v4(model%cloud, &
!!$                                   & now, &
!!$                                   & timestep,&
!!$                                   & model%rules%diffusionMethod, &
!!$                                   & model%rules%ifThermodiffusion, &
!!$                                   & 1.0, & ! weight_past
!!$                                   & model%meteo_ptr, & 
!!$                                   & model%disp_buf_ptr, &
!!$                                   & model%wdr, &
!!$                                   & model%rules%chemicalRules, &
!!$                                   & linearizations(istep))

!!$      call transform_pollution_cloud_v5(model%cloud, now, timestep, &
!!$                                      & model%meteo_ptr, model%disp_buf_ptr, &
!!$                                      & model%Meteo_input, & 
!!$                                      & model%wdr, model%rules%chemicalRules)

!!$      call msg('sum',sum(abs(linearizations(istep)%cm_x_before_adv)) / (nx*ny*nz*nsp))
!!$    
!!$      call flip_moment_and_basic_var(cm_to_moment, fu_advection_moment_X_MM_ptr(model%cloud), &
!!$                                   & fu_concMM_ptr(model%cloud), to_moment)
!!$      call flip_moment_and_basic_var(cm_to_moment, fu_advection_moment_Y_MM_ptr(model%cloud), &
!!$                                   & fu_concMM_ptr(model%cloud), to_moment)
!!$      call flip_moment_and_basic_var(cm_to_moment, fu_advection_moment_Z_MM_ptr(model%cloud), &
!!$                                   & fu_concMM_ptr(model%cloud), to_moment)
!!$      call collect_total_masses(model%cloud)
!!$      call report_total_masses(model%cloud, 99, .true.)
!!$      
!!$    end do

    call model_forward(model, model%obs_ptr, model%rules)

    m_arr = map_c%arm(1:nsp,1,1:nz,1:nx,1:ny)
    call msg('RIGHT', fu_inner(m_arr, arr2))
    
    ! Adjoint
    
    MMptr => fu_advection_moment_X_MM_ptr(model%cloud)
    call reset_map_to_val(MMptr, 0.0)
    MMptr => fu_advection_moment_Y_MM_ptr(model%cloud)
    call reset_map_to_val(MMptr, 0.0)
    MMptr => fu_advection_moment_Z_MM_ptr(model%cloud)
    call reset_map_to_val(MMptr, 0.0)

    map_c%arm(1:nsp,1,1:nz,1:nx,1:ny) = arr2

!!$    do istep = 1, nsteps
!!$      CALL advect_pollution_cloud_v4(model%cloud, &
!!$                                   & now, &
!!$                                   & fu_opposite(timestep),&
!!$                                   & model%rules%diffusionMethod, &
!!$                                   & model%rules%ifThermodiffusion, &
!!$                                   & 1.0, & ! weight_past
!!$                                   & model%meteo_ptr, & 
!!$                                   & model%disp_buf_ptr, &
!!$                                   & model%wdr, &
!!$                                   & model%rules%chemicalRules, &
!!$                                   & linearizations(nsteps-istep+1))

!!$      call transform_pollution_cloud_v5(model%cloud, now, fu_opposite(timestep), &
!!$                                      & model%meteo_ptr, model%disp_buf_ptr, &
!!$                                      & model%Meteo_input, & 
!!$                                      & model%wdr, model%rules%chemicalRules)
!!$
!!$      
!!$    end do
    
    call model_adjoint(model, model%obs_ptr, model%rules)
    m_arr = map_c%arm(1:nsp,1,1:nz,1:nx,1:ny)    
    call msg('LEFT', fu_inner(m_arr, arr1))
    
    stop
    
  contains
    
    real function fu_inner(xx, y) result(inner)
      implicit none
      real, dimension(:,:,:,:), intent(in) :: xx, y
      
      integer :: ind_1, ind_2, ind_3, ind_sp
      real :: dz, dxdy
      real(r8k) :: inner_8

      inner = 0.0
      inner_8 = 0

      do ind_1 = 1, size(xx, 3)
        do ind_2 = 1, size(xx, 4)
          do ind_3 = 1, size(xx, 2)
            do ind_sp = 1, size(xx,1)
              dz = 1.0!dispersion_layer_z_size(ind_3)
              !dz = dispersion_layer_z_size(ind_3)
              dxdy = 1.0!fu_cell_size(dispersion_grid, ind_1, ind_2)
              inner_8 = inner_8 + xx(ind_sp,ind_3,ind_1,ind_2)*y(ind_sp,ind_3,ind_1,ind_2)/(dz*dxdy)
            end do
          end do
        end do
      end do
      inner = inner_8
    end function fu_inner
  end subroutine adj_test


  subroutine run_forecast_dummy(control, model)
    ! Doesn't run the model but perturbs the fields with gaussian noise and processes the
    ! observations.
    implicit none
    type(model_container), intent(inout) :: model
    type(da_control), intent(in) :: control
    
    type(silja_time) :: now
    type(Tmass_map), pointer ::cnc_ptr
    real, parameter :: autocorr = 0.9
    real, dimension(1) :: value, noise
    real, dimension(:), pointer :: thickness
    integer :: ilev
    logical, save :: first = .true.

    if (.not. fu_in_physical_space(control)) then
      call set_error('Control variable not in physical space', 'get_obs_cost')
      return
    end if
    
    model%rules%darules%allow_negative = .true.
    if (model%obs_ptr%hasObservations) call reset_all(model%obs_ptr, .true.)
    cnc_ptr => fu_concMM_ptr(model%cloud)
    call msg('sum cloud before:', sum(cnc_ptr%arm(8,:,:,:,:)))
    call msg('sum control:', sum(fu_values_ptr(control)))
    call control2cloud(control, model, handle_negatives=fu_mode(control) == da_initial_state)
    call msg('sum cloud after:', sum(cnc_ptr%arm(8,:,:,:,:)))

    now = model%rules%startTime

    thickness => fu_work_array()
    thickness(1:nz_dispersion) = (/(fu_layer_thickness_m(fu_level(dispersion_vertical, ilev)), &
                                  & ilev = 1, nz_dispersion)/)
    !call msg('thick:', thickness(1:nz_dispersion))
    if (first) then
      call random_normal(noise)
      call putcnc(noise(1), 0.0, value(1))
      first = .false.
      return
    end if
    do while(fu_between_times(now, &
                            & model%rules%startTime, &
                            & model%rules%startTime + model%rules%periodToCompute, &
                            & .true.))
      call random_normal(noise)
      !call msg('New sample:', sample_new(1))
            
      call putcnc(noise(1), autocorr, value(1))
      call msg('Concentration:', value(1))
      call msg('sum cloud:', sum(cnc_ptr%arm(8,:,:,:,:)))
      call msg('Run forecast dummy: ' // fu_str(now))
      call observeAll(model%obs_ptr, model%cloud,&
                    & model%meteo_ptr, model%disp_buf_ptr, model%rules%chemicalRules, &
                    & model%rules%timestep, now)
      now = now + model%rules%timestep

    end do

    call free_work_array(thickness)

  contains

    subroutine putcnc(noise, weight_past, new_value)
      implicit none
      real :: noise, new_value, weight_past

      integer :: ix, iy, iz
      real :: volume

      do iy = 1, ny_dispersion
        do ix = 1, nx_dispersion
          do iz = 1, nz_dispersion
            volume = thickness(iz) * fu_cell_size(dispersion_grid, ix, iy)
            cnc_ptr%arm(:,:,iz,ix,iy) = &
                 & (weight_past*cnc_ptr%arm(:,:,iz,ix,iy)/volume + (1-weight_past)*noise)*volume
          end do
        end do
      end do

      volume = thickness(1) * fu_cell_size(dispersion_grid, 1, 1)
      new_value = cnc_ptr%arm(1,1,1,1,1)/volume

    end subroutine putcnc

  end subroutine run_forecast_dummy


end module da_interface
