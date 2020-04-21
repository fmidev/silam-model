module observation_server
  
  ! This module contains routines for dealing with the set of observations - mainly, the
  ! structures for storing them, and subroutines for processing them during forward and 
  ! adjoint runs during data assimilation. Currently, also DA_rules is defined here, even
  ! though it is mainly used by the module data_assimilation. 
  !
  ! Them main idea in the SILAM 4D-VAR implementation is that every observation interacts 
  ! with the model only through its routines 'observe' and 'inject'. These are called for 
  ! every observation at every timestep during forward and adjoint runs. Everything behind 
  ! these routines, including data storage etc., should be invisible for the model. Observation 
  ! server acts as a proxy between model (dispersion_supplementary) and the observations.

!  use observations_in_situ
  use chemistry_manager
  use pollution_cloud
  use source_apportionment
!  use optical_density
  use observations_vertical 
  use observations_dose_rate

  implicit none
  
  private

  public observeAll
  public injectAll

  public destroy_observations
  public dump_observations
  public reset_all

  public fu_obs_data
  public fu_obs_variance
  public collect_model_data
  public collect_obs_data
  public get_obs_data
  public collect_variance
  public set_mdl_data
  public set_observations
  public fu_number_of_observed_cells

  public test_read_timeseries
  public get_localisation
  
  private read_observations
  private searchStationWithID
  private stations_from_namelist_ptr
  private stations_from_namelist_no_ptr

  interface get_localisation
     module procedure get_localisation_all
  end interface

  !integer, parameter, public :: DA_INITIAL_STATE = 13001, DA_EMISSION = 13002, DA_MODEL_ONLY = 13000, &
  !     & DA_INITIAL_STATE_MOMENT = 13003, DA_EMISSION_CORRECTION = 13004, DA_EMISSION_AND_INITIAL = 13005, &
  !     & DA_EMISSION_CORRECTION_LOG = 13006, DA_EMISSION_TIME_HEIGHT = 13007
  integer, parameter, public :: LINE_SEARCH_ARMIJO = 1, LINE_SEARCH_MS = 2
  integer, parameter, public :: steepest_descent_flag = 91, m1qn3_flag = 92, l_bfgs_b_flag = 93
  integer, parameter, public :: max_column_observations = 100000
  integer, parameter, public :: max_eruption_observations = 1000

  !**********************************************************************************
  !
  ! observationPointers
  !
  !**********************************************************************************
  ! These are the main datatypes defined in this file. 
  ! This structure will contain pointers to all observations loaded. It should be mainly
  ! handled by the observation_server.

  type observationPointers
     type(inSituObservation), dimension(:), pointer :: observationsInSitu
     type(t_vertical_observation), dimension(:), pointer :: observationsVertical
     type(t_eruptionObservation), dimension(:), pointer :: observationsEruption 
     type(inSituObservation), dimension(:), pointer :: observationsDoseRate
     type(t_dose_rate_obs_addition), dimension(:), pointer :: DoseRateAddition

     real, dimension(:), pointer :: obs_values
     real, dimension(:), pointer :: obs_variance

     integer :: nInSituObservations = 0, nVerticalObservations = 0, nDoseRateObservations = 0
     integer :: nEruptionObservations = 0
     integer :: obs_size = 0
     logical :: hasObservations = .false. 
  end type ObservationPointers
  public observationPointers

    
contains
 
  !************************************************************************************

  subroutine stations_from_namelist_no_ptr(nlPtr, station_list, nstations, grid)
    implicit none
    type(Tsilam_namelist), pointer :: nlPtr
    type(observationStation), dimension(:), pointer :: station_list
    integer, intent(out) :: nstations
    type(silja_grid), intent(in) :: grid

    integer :: nitems, stat, i, iStat
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: items
    type(silam_sp) :: workstring
    real :: lat, lon, hgt
    character(len=STATION_LABEL_LENGTH) :: label
    character(len=clen) :: name

    nullify(items)
    workstring%sp => fu_work_string()
    
    call get_items(nlPtr, 'station', items, nstations)
    
    if (nstations < 1) then
      call set_error('No stations found', 'stations_from_namelist')
      nullify(station_list)
      return
    end if

    allocate(station_list(nstations), stat=stat)
    if (stat /= 0) then
      call set_error('Allocate failed', 'stations_from_namelist')
      return
    end if
    
    !
    ! Have to take care of the stations outside the modellign domain
    ! Shrink the list kicking them out
    !
    iStat = 1
    do i = 1, nstations
      workstring%sp = fu_content(items(i))
      read(workstring%sp, fmt=*, iostat=stat) label, lat, lon, hgt, name
      if (stat /= 0) then
        call msg('Problem with record: ' // workstring%sp)
        call set_error('Problem parsing station list', 'stations_from_namelist')
        return
      end if
      station_list(iStat) = fu_initObservationStation(label, name, lon, lat, hgt, grid)
      if(defined(station_list(iStat))) iStat = iStat + 1  ! if initialization was successful
    end do
    nstations = iStat - 1 ! the number of actually initialised stations
    
    call free_work_array(workstring%sp)

  end subroutine stations_from_namelist_no_ptr

  !************************************************************************************

  subroutine stations_from_namelist_ptr(nlPtr, station_list, nstations, grid)
    implicit none
    type(Tsilam_namelist), pointer :: nlPtr
    type(observationStationPtr), dimension(:), pointer :: station_list
    integer, intent(out) :: nstations
    type(silja_grid), intent(in) :: grid

    integer :: nitems, stat, ii, iStat
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: items
    type(silam_sp) :: workstring
    real :: lat, lon, hgt
    character(len=STATION_LABEL_LENGTH) :: label
    character(len=clen) :: name
    logical :: ifOK

    nullify(items)
    workstring%sp => fu_work_string()
    
    call get_items(nlPtr, 'station', items, nstations)
    
    if (nstations < 1) then
      call set_error('No stations found', 'stations_from_namelist')
      nullify(station_list)
      return
    end if

    allocate(station_list(nstations), stat=stat)
    if (stat /= 0) then
      call set_error('Allocate failed', 'stations_from_namelist')
      return
    end if
    
    iStat = 1
    do ii = 1, nstations
      workstring%sp = fu_content(items(ii))
      read(workstring%sp, fmt=*, iostat=stat) label, lat, lon, hgt, name
      if (stat /= 0) then
        call msg('Problem with record: ' // workstring%sp)
        call set_error('Problem parsing station list', 'stations_from_namelist')
        return
      end if
      allocate(station_list(iStat)%ptr, stat=stat)
      if (stat /= 0) then
        call set_error('Allocate failed', 'stations_from_namelist')
        return
      end if
      station_list(iStat)%ptr = fu_initObservationStation(label, name, lon, lat, hgt, grid)
      if(defined(station_list(iStat)%ptr))then
        iStat = iStat + 1
      else
        deallocate(station_list(iStat)%ptr)
      endif
    end do
    nstations = iStat - 1

    call free_work_array(workstring%sp)
    
  end subroutine stations_from_namelist_ptr

  !************************************************************************************
  
 subroutine read_observations(obs_start, obs_period, obs_list_path, obs_items, nObsItems, &
                            & station_list, nstations, &
                            & in_situ_array, n_in_situ, &
                            & vert_obs_array, n_vert_obs, &
                            & eruption_obs_array, n_eruption, &
                            & dose_rate_addition_array, n_dose_rate, &
                            & transport_species, optical_species, n_transp_species, n_opt_species)
    implicit none
    type(silja_time), intent(in) :: obs_start
    type(silja_interval), intent(in) :: obs_period
    character(len=*), intent(in) :: obs_list_path
    character(len=*), dimension(:), pointer :: obs_items
    type(observationStationPtr), dimension(:), intent(in) :: station_list
    integer, intent(in) :: nstations
    integer, intent(inout) :: nObsItems
    type(silam_species), dimension(:), intent(in) :: transport_species
    type(silam_species), dimension(:), pointer :: optical_species
    integer, intent(in) :: n_transp_species, n_opt_species
    ! out
    integer, intent(out) :: n_in_situ, n_dose_rate, n_eruption
    type(inSituObservation), dimension(:), pointer :: in_situ_array
    type(t_dose_rate_obs_addition), dimension(:), pointer :: dose_rate_addition_array
    type(t_vertical_observation), dimension(:), pointer :: vert_obs_array

    type(t_eruptionObservation), dimension(:), pointer :: eruption_obs_array

    integer :: status, file_unit, nread, i, icounter=-1, n_vert_obs, ind_obs, n_vert_obs_def
    character(len=64) :: obs_type, obs_unit, var_name
    character(len=fnlen) :: file_name, obs_item, filename_templ, file_name_old
    type(inSituObservation), dimension(:), allocatable :: observations_tmp, dose_rate_obs_tmp
    type(t_dose_rate_obs_addition), dimension(:), allocatable :: dose_rate_addition_tmp
    type(t_vertical_observation), dimension(:), allocatable :: vert_obs_tmp
    type(t_eruptionObservation), dimension(:), allocatable :: eruption_obs_tmp
    type(silam_species), dimension(:), pointer :: p_opt_species, p_obs_species
    type(grads_template) :: obs_templ
    character(len=fnlen), dimension(:), pointer :: p_obs_items
    integer :: obs_index
    logical :: dealloc_obs_items, ifGroundToo, force_instant
    type(silja_time) :: time_in_filename
    !!FIXME some smarter allocation needed here
    integer, parameter :: max_n_obs_per_time_window = 1
    character(len=*), parameter :: sub_name = 'read_observations'

    if (n_opt_species > 0 .and. .not. associated(optical_species)) then
      call set_error('Inconsistent arguments', sub_name)
      return
    end if
    
    call msg('Read observations, time limits: ' // fu_str(obs_start) // ' ' &
           & // fu_str(obs_start+obs_period))
    call start_count(sub_name)

    force_instant = (nint(fu_sec(obs_period)) == 0)

    n_in_situ = 0
    n_vert_obs = 0
    n_dose_rate = 0
    n_eruption = 0

!     call msg("Before allocating observations memory usage (kB)",  fu_system_mem_usage())
    allocate(observations_tmp(max_n_obs_per_time_window*nstations*n_transp_species), &
           & dose_rate_obs_tmp(nstations*n_transp_species), &
           & dose_rate_addition_tmp(nstations*n_transp_species), &
           & vert_obs_tmp(max_column_observations), &
           & eruption_obs_tmp(max_eruption_observations), &
           & stat=status)
!     call msg("After allocating observations memory usage (kB)",  fu_system_mem_usage())
    
    if (fu_fails(status == 0, 'Allocate failed', sub_name)) return

    file_unit = fu_next_free_unit()

    if (obs_list_path == '') then
      p_obs_items => obs_items
      dealloc_obs_items = .false.
    else
      call decode_template_string(obs_list_path, obs_templ)
      if (error) return
      call expand_template(obs_templ, obs_start, file_name)
      if (error) return
      call add_obs_items_from_list(obs_items, nObsItems, file_name, p_obs_items, nObsItems)
      if (error) return
      dealloc_obs_items = .true.
    end if
    
    do i = 1, nObsItems
      obs_item = p_obs_items(i)
      read(obs_item, fmt=*, iostat=status) obs_type
      if (fu_fails(status == 0, 'Failed to parse:' // trim(obs_item), sub_name)) return

      if (obs_type == 'cnc') then
        time_in_filename = obs_start
        if (fu_fails(nstations > 0, 'Surface observation given but no stations', sub_name)) return
        read(obs_item, fmt=*, iostat=status) obs_type, obs_unit, file_name
        if (fu_fails(status == 0, 'Failed to parse: ' // trim(obs_item), sub_name)) return
        call decode_template_string(file_name, obs_templ)
        if (error) return

        do while (time_in_filename <= obs_start + obs_period)
          call expand_template(obs_templ, time_in_filename, file_name)
          if (file_name /= file_name_old) then
            if (.not. error) then
              open(file_unit, file=file_name, status='old', action='read', iostat=status)
              if (fu_fails(status == 0, 'Failed to open: ' // trim(file_name), sub_name)) return
              call msg('Reading cnc observations from ' // file_name)
              call timeseries_from_file(file_unit, obs_unit, observations_tmp(n_in_situ+1:), nread, force_instant)
              close(file_unit)
              n_in_situ = n_in_situ + nread
              if (error) return
            end if
            file_name_old = file_name
          else
            call unset_error(sub_name)
          end if
          time_in_filename = time_in_filename + fu_set_interval_h(1)
        end do
        
      elseif (obs_type == 'dose_rate_no_ground') then
        ifGroundToo = .false.
        read(obs_item, fmt=*, iostat=status) obs_type, file_name
        if (fu_fails(status == 0, 'Failed to parse: ' // trim(obs_item), sub_name)) return
        call decode_template_string(file_name, obs_templ)
        if (error) return
        call expand_template(obs_templ, obs_start, file_name)
        if (.not. error) then
          open(file_unit, file=file_name, status='old', action='read', iostat=status)
          if (fu_fails(status == 0, 'Failed to open: ' // trim(file_name), sub_name)) return
          call msg('Reading dose rate no-ground observations from ' // file_name)
          call timeseries_from_file_dose_rate(ifGroundToo, file_unit, obs_unit, dose_rate_obs_tmp(n_dose_rate+1:), &
                                            & dose_rate_addition_tmp(n_dose_rate+1:), nread)
          close(file_unit)
          n_dose_rate = n_dose_rate + nread
        else
          call unset_error(sub_name)
        end if
         
      elseif (obs_type == 'dose_rate_with_ground') then
        ifGroundToo = .true.
        read(obs_item, fmt=*, iostat=status) obs_type, file_name
        if (fu_fails(status == 0, 'Failed to parse: ' // trim(obs_item), sub_name)) return
        call decode_template_string(file_name, obs_templ)
        if (error) return
        call expand_template(obs_templ, obs_start, file_name)
        if (.not. error) then
          open(file_unit, file=file_name, status='old', action='read', iostat=status)
          if (fu_fails(status == 0, 'Failed to open: ' // trim(file_name), sub_name)) return
          call msg('Reading dose rate with-ground observations from ' // file_name)
          call timeseries_from_file_dose_rate(ifGroundToo, file_unit, obs_unit, dose_rate_obs_tmp(n_dose_rate+1:), &
                                            & dose_rate_addition_tmp(n_dose_rate+1:), nread)
          close(file_unit)
          n_dose_rate = n_dose_rate + nread
        else
          call unset_error(sub_name)
        end if


      elseif (obs_type == var_name_aod) then
        time_in_filename = obs_start
        read(obs_item, fmt=*, iostat=status) obs_type, file_name
        if (fu_fails(status == 0, 'Failed to parse: ' // trim(obs_item), sub_name)) return
        call decode_template_string(file_name, obs_templ)
        if (error) return

        do while (time_in_filename <= obs_start + obs_period)
          call expand_template(obs_templ, time_in_filename, file_name)
          if (.not. error) then
            call msg('Reading AOT from '  // trim(file_name))
            ! obs_species = null() - not used                                                                                          
            nullify(p_obs_species)
            call set_vert_obs_from_nc(file_name, var_name_aod, p_obs_species, &
                 & transport_species, optical_species, &
                 & obs_start, obs_start+obs_period, &
                 & vert_obs_tmp(n_vert_obs+1:), nread)
            !if ((rules%observation_stdev > 0) .and. (nread > 0)) then                                                                 
            !  call msg('nread', nread)                                                                                                
            !  do obs_index = n_vert_obs+1, nread                                                                                      
            !    vert_obs_tmp(obs_index)%variance = (rules%observation_stdev)**2                                                       
            !  end do                                                                                                                  
            !end if                                                                                                                    
            if (error) return
            n_vert_obs = n_vert_obs + nread
          else
            call unset_error(sub_name)
          end if
          time_in_filename = time_in_filename + fu_set_interval_h(1)
        end do
         
      else if (obs_type == 'eruption') then
        read(obs_item, fmt=*, iostat=status) obs_type, obs_unit, file_name
        if (fu_fails(status == 0, 'Failed to parse: ' // trim(obs_item), sub_name)) return
        call decode_template_string(file_name, obs_templ)
        if (error) return
        call expand_template(obs_templ, obs_start, file_name)
        if (.not. error) then
          open(file_unit, file=file_name, status='old', action='read', iostat=status)
          if (fu_fails(status == 0, 'Failed to open: ' // trim(file_name), sub_name)) return
          call msg('Reading eruption observations from ' // file_name)
          call eruption_timeseries_from_file(file_unit, obs_unit, eruption_obs_tmp(n_eruption+1:), nread)
          close(file_unit)
          n_eruption = n_eruption + nread
          if (error) return
        else
          call unset_error(sub_name)
        end if

      else if (obs_type == 'vertical') then
        time_in_filename = obs_start
        read(obs_item, fmt=*, iostat=status) obs_type, var_name, file_name
        if (fu_fails(status == 0, 'Failed to parse: ' // trim(obs_item), sub_name)) return
        call decode_template_string(file_name, obs_templ)
        if (error) return
        do while (time_in_filename <= obs_start + obs_period)
          call expand_template(obs_templ, time_in_filename, file_name)
          if (.not. error) then
            ! var_name is also the observation cocktail:                                                         
            call msg('Reading vertical observation from '  // trim(file_name))
            call get_observed_species(var_name, p_obs_species)
            if (error) return
            call set_vert_obs_from_nc(file_name, var_name, p_obs_species, &
                 & transport_species, optical_species, &
                 & obs_start, obs_start+obs_period, &
                 & vert_obs_tmp(n_vert_obs+1:), nread)
            if (error) return
            deallocate(p_obs_species)
            n_vert_obs = n_vert_obs + nread
         else
            call unset_error(sub_name)
         end if
         time_in_filename = time_in_filename + fu_set_interval_h(1)
       end do
     else
       call set_error('Unknown observation type:' // trim(obs_type), sub_name)
       return
     end if
        
 end do

   if (fu_fails(n_in_situ + n_vert_obs + n_dose_rate + n_eruption > 0, 'No observations read succesfully', sub_name)) return
    
    if (n_in_situ > 0) then
      allocate(in_situ_array(n_in_situ), stat=status)
      if (fu_fails(status == 0, 'Allocate failed', sub_name)) return
      in_situ_array = observations_tmp(1:n_in_situ)
    end if
    
    if (n_dose_rate > 0) then
      allocate(dose_rate_addition_array(n_dose_rate), stat=status)
      if (fu_fails(status == 0, 'Allocate failed', sub_name)) return
      dose_rate_addition_array = dose_rate_addition_tmp(1:n_dose_rate)
    end if

    if (n_eruption > 0) then
      allocate(eruption_obs_array(n_eruption), stat=status)
      if (fu_fails(status == 0, 'Allocate failed', sub_name)) return
      eruption_obs_array = eruption_obs_tmp(1:n_eruption)
    end if

    if (n_vert_obs > 0) then
      ! some observations could be missing because they are outside, etc. Ignore them.
      n_vert_obs_def = count(defined(vert_obs_tmp(1:n_vert_obs)))
      allocate(vert_obs_array(n_vert_obs_def), stat=status)
      if (fu_fails(status == 0, 'Allocate failed', sub_name)) return
      n_vert_obs_def = 0
      do ind_obs = 1, n_vert_obs
        if (defined(vert_obs_tmp(ind_obs))) then
          n_vert_obs_def = n_vert_obs_def + 1
        else
          cycle
        end if
        vert_obs_array(n_vert_obs_def) = vert_obs_tmp(ind_obs)
      end do
      call msg('Vertical observation loaded/defined:', n_vert_obs, n_vert_obs_def)
      n_vert_obs = n_vert_obs_def ! return value
    end if
    
    deallocate(observations_tmp, vert_obs_tmp, dose_rate_addition_tmp, dose_rate_obs_tmp, eruption_obs_tmp)
    if (dealloc_obs_items) deallocate(p_obs_items)
    
    call stop_count(sub_name)
    !call report_time(icounter, chCounterNm=sub_name)

    CONTAINS

    !==================================================================================
                            
    subroutine get_observed_species(cockt_name, species_list_ptr)
      implicit none
      character(len=*), intent(in) :: cockt_name
      type(silam_species), dimension(:), pointer :: species_list_ptr

      type(Tcocktail_descr) :: descr
      type(silam_species), dimension(:), pointer :: cockt_species_ptr
      integer :: n_obs_species, stat
      logical :: ifSpecies

      nullify(species_list_ptr)
      call set_cocktail_description(cockt_name, descr, ifSpecies)
      if (error) return
      call get_inventory(descr, cockt_species_ptr, n_obs_species)
      if (fu_fails(n_obs_species > 0, 'Descriptor has no species', 'get_observed_species')) return

      ! Copy the species to allow deallocating the cocktail
      !
      allocate(species_list_ptr(n_obs_species), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', 'get_observed_species')) return
      species_list_ptr(1:n_obs_species) = cockt_species_ptr(1:n_obs_species)
      call destroy(descr)
      ! cockt_species_ptr not deallocated
      
    end subroutine get_observed_species

    !================================================================================

    subroutine timeseries_from_file(file_unit, obs_unit, observations, nread, force_instant)
      implicit none
      integer, intent(in) :: file_unit
      character(len=*), intent(in) :: obs_unit
      type(inSituObservation), dimension(:), intent(inout) :: observations
      logical, intent(in) :: force_instant !Disregard duration and force it to zero

      integer, intent(out) :: nread
      
      logical :: eof, same_series, found, same_cocktail, store_value, flush_values
      integer :: iostat

      character(len=station_id_length) :: id, prev_id
      character(substNmLen) :: cockt_name, prev_cockt_name
      character(len=255) :: line
      integer :: year, month, day, hour, itime, quantity, n_cockt_species
      real :: duration_hours, value, modeval, wavelength, stdev, variance
      type(silja_interval) :: duration
      type(silja_time) :: time, prev_time
      type(observationStation), pointer :: station
      real, dimension(:), pointer :: values_tmp, variances_tmp
      type(silja_time), dimension(:), allocatable :: times_tmp
      type(silja_interval), dimension(:), allocatable :: durations_tmp
      type(silam_species), dimension(:), pointer :: cockt_species
      type(Tcocktail_descr) :: cockt_descr

      !
      integer, dimension(max_species) :: ind_obs2transp
      real, dimension(max_species) :: scale_transp2obs
      logical, dimension(max_species) :: is_aerosol
      integer :: status, nx, ny, i, isp_transp, isp_obs
      type(silam_material), pointer :: material
      type(chemical_adaptor) :: adaptor
      logical :: ifSpecies

      prev_time = time_missing
      prev_cockt_name = ''
      prev_id = ''
      nread = 0
      itime = 0
      values_tmp => fu_work_array()
      variances_tmp => fu_work_array()
      allocate(times_tmp(worksize), durations_tmp(worksize), stat=iostat)
      if (fu_fails(iostat == 0,'Allocate failed', 'timeseries_from_file'))return

      !read(file_unit,fmt=*) line

      do
        call next_line_from_input_file(file_unit, line, eof)
        if (eof) iostat = -1
        !read(unit, fmt='(A)', iostat=iostat) line
        if (iostat > 0) then
          call set_error('Failed to read record', 'timeseries_from_file')
          return
        end if
        store_value = .false.
        flush_values = .false.
        if (iostat == 0) then
          read(line, fmt=*, iostat=iostat) &
               & id, cockt_name, year, month, day, hour, duration_hours, value, stdev
          variance = stdev**2
          if (iostat /= 0) then
            call set_error('Failed to parse record: ' // trim(line), 'timeseries_from_file')
            return
          end if

          if (force_instant) then
            duration = zero_interval
          else
            duration = fu_set_interval_sec(duration_hours*3600.0)
          endif
          time = fu_set_time_utc(year, month, day, hour, 0, 0.0)

          ! No duct tape here, please!!!
!
!          if (duration == zero_interval) then
!            ! "instant" observation. To avoid needing special treatment later, we'll set
!            ! the duration to 20 sec, but no shifting due to averaging, since there wasn't any.
!            duration = fu_set_interval_sec(20.0)
!          end if
          store_value = .true.

          if (time < obs_start) then 
            store_value = .false.
          else
            store_value = .true.
          end if
          if (time > obs_start + obs_period) then
            flush_values = .true.
            store_value = .false.
          end if

          if (defined(prev_time)) then
            if (prev_id == id .and. prev_time >= time) then
              call set_error('Cannot read: timeseries not sorted', 'timeseries_from_file')
              return
            end if
          end if
          same_series = (cockt_name == prev_cockt_name .and. id == prev_id)
          ! disregard same_series if this is the first line read:
          if (prev_cockt_name /= '' .and. .not. same_series) flush_values = .true.

          if (prev_cockt_name /= cockt_name) then
            !!if (debug_level > 0) 
            call msg("New cocktail:"+cockt_name)
            call set_cocktail_description(cockt_name, cockt_descr, ifSpecies)
            call get_inventory(cockt_descr, cockt_species, n_cockt_species)
            call create_adaptor(cockt_species, transport_species, adaptor)

            ind_obs2transp(1:n_cockt_species) = adaptor%isp(1:n_cockt_species)
            is_aerosol(1:n_cockt_species) = .false.

            do isp_obs = 1, n_cockt_species
              material => fu_material(transport_species(ind_obs2transp(isp_obs)))
              if (error) return
              ! From basic unit to the unit of observation - preparation to the obs operator.
              scale_transp2obs(isp_obs) = fu_conversion_factor(fu_basic_mass_unit(material), obs_unit, material)
              if (error) return
              if (fu_nominal_d(fu_mode(transport_species(ind_obs2transp(isp_obs)))) /= int_missing) then   
                is_aerosol(isp_obs) = .true.    
              end if

              if (fu_fails(scale_transp2obs(isp_obs) > 0.0, 'Bad conversion factor to obs unit', 'timeseries_from_file')) return
!!                call msg('Observed species, transport species', isp_obs, ind_obs2transp(isp_obs))
              call msg('Factor transp2obs' //trim(fu_str(cockt_species(isp_obs))), scale_transp2obs(isp_obs))
            end do
          end if
        else ! eof
          flush_values = itime > 0
        end if ! read ok
        
        flush_values = flush_values .and. itime > 0

        if (flush_values) then
          !call msg(id)
          call searchStationWithID(prev_id, station_list, nstations, station, found)
          if (.not. found) then
            call msg_warning('Station ' // trim(prev_id) // ' not found')
          else
            nread = nread + 1
            if (nread > size(observations)) then
              call set_error('Too many observations', 'timeseries_from_file')
              return
            end if
            observations(nread) = fu_initInSituObservation(times_tmp, &
                                                         & durations_tmp, &
                                                         & values_tmp, &
                                                         & variances_tmp, &
                                                         & itime, &
                                                         & variable_variance, &
                                                         & ind_obs2transp, &
                                                         & scale_transp2obs, &
                                                         & is_aerosol, &
                                                         & n_cockt_species, &
                                                         & station, &
                                                         & level_missing, &
                                                         & dispersion_vertical, &
                                                         & cockt_name)
          end if ! not found
          itime = 0
        end if ! have new series

        if (iostat < 0) exit
        
        if (store_value) then
          itime = itime + 1
          values_tmp(itime) = value
          durations_tmp(itime) = duration
          times_tmp(itime) = time
          variances_tmp(itime) = variance
        end if
        prev_time = time
        prev_id = id
        prev_cockt_name = cockt_name
      end do ! loop over lines
      
      deallocate(times_tmp, durations_tmp)
      call free_work_array(values_tmp)
      call free_work_array(variances_tmp)
    end subroutine timeseries_from_file

    !========================================================================================

    subroutine eruption_timeseries_from_file(file_unit, obs_unit, observations, nread)
      implicit none
      integer, intent(in) :: file_unit
      character(len=*), intent(in) :: obs_unit
      type(t_eruptionObservation), dimension(:), intent(inout) :: observations

      integer, intent(out) :: nread

      logical :: eof
      integer :: iostat

      character(len=255) :: line
      integer :: year, month, day, hour
      real :: duration_hours, value, stdev, variance, lat, lon
      real, dimension(:), pointer :: value_tmp, variance_tmp, lat_tmp, lon_tmp
      type(silja_interval) :: duration
      type(silja_interval), dimension(:), allocatable :: duration_tmp
      !type(silja_interval), pointer :: duration_tmp
      type(silja_time) :: time
      type(silja_time), dimension(:), allocatable :: time_tmp
      !type(silja_time), pointer :: time_tmp
      !                                                                                                                                                                                                                                                                                                                    
      integer :: status, nx, ny, i
      integer, parameter :: tmp_arr_size = 20

      value_tmp => fu_work_array()
      variance_tmp => fu_work_array()
      lat_tmp => fu_work_array()
      lon_tmp => fu_work_array()

      allocate(time_tmp(tmp_arr_size), duration_tmp(tmp_arr_size), stat=iostat)
      if (iostat /= 0) then
        call set_error('Allocate failed', 'timeseries_from_file')
        return
      end if

      nread = 0
     
      do
        call next_line_from_input_file(file_unit, line, eof)
        if (eof) iostat = -1
        !read(unit, fmt='(A)', iostat=iostat) line                                                                                                                                                                                                                                                                          
        if (iostat > 0) then
          call set_error('Failed to read record', 'timeseries_from_file')
          return
        end if

        if (iostat == 0) then
          read(line, fmt=*, iostat=iostat) &
               & year, month, day, hour, duration_hours, value, stdev, lat, lon
          variance = stdev**2
          if (iostat /= 0) then
            call set_error('Failed to parse record: ' // trim(line), 'eruption_timeseries_from_file')
            return
          end if

          duration = fu_set_interval_sec(duration_hours*3600.0)
          time = fu_set_time_utc(year, month, day, hour, 0, 0.0)
          
          if (duration == zero_interval) then
            ! "instant" observation. To avoid needing special treatment later, we'll set                                                                                                                                                                                                                                   
            ! the duration to 20 sec, but no shifting due to averaging, since there wasn't any.                                                                                                                                                                                                                            
            duration = fu_set_interval_sec(20.0)
          end if

          value_tmp(1) = value
          duration_tmp(1) = duration
          time_tmp(1) = time
          variance_tmp(1) = variance
          lat_tmp(1) = lat
          lon_tmp(1) = lon

          !if ((time > rules%da_begin) .and. (time < rules%da_begin + rules%da_window)) then
 
            nread = nread + 1
            if (nread > size(observations)) then
              call set_error('Too many observations', 'timeseries_from_file')
              return
            end if
            observations(nread) = fu_initEruptionObservation(time_tmp, &
                 & duration_tmp, value_tmp, variance_tmp, lat_tmp, lon_tmp, 1, 'eruption')                                                                                                                                                                                                                                                       
          !end if                                                                                                                                                                                                
        end if
        
        if (iostat < 0) exit
      end do ! loop over lines
      !deallocate(value_tmp, variance_tmp, lat_tmp, lon_tmp, duration_tmp, time_tmp)
      deallocate(time_tmp, duration_tmp)
      call free_work_array(value_tmp)
      call free_work_array(variance_tmp)
      call free_work_array(lon_tmp)
      call free_work_array(lat_tmp)

    end subroutine eruption_timeseries_from_file

    !==========================================================================

    subroutine timeseries_from_file_dose_rate(ifGroundToo, file_unit, obs_unit, observations_dose_rate, dose_rate_addition, nread)
      !edit: no cocktail here compared to normal timeseries from file!
      !dose rate measurement doesn't know what is emitting
      implicit none
      logical, intent(in) :: ifGroundToo
      integer, intent(in) :: file_unit
      character(len=*), intent(in) :: obs_unit

      type(inSituObservation), dimension(:), intent(inout) :: observations_dose_rate
      
      type(t_dose_rate_obs_addition), dimension(:), intent(out) :: dose_rate_addition
      integer, intent(out) :: nread
      
      logical :: eof, same_series, found, store_value, flush_values
      integer :: iostat
      character(len=station_id_length) :: id, prev_id
      character(len=255) :: line
      integer :: year, month, day, hour, itime, quantity
      real :: duration_hours, value, modeval, wavelength, stdev, variance
      type(silja_interval) :: duration
      type(silja_time) :: time, prev_time
      type(observationStation), pointer :: station
      real, dimension(:), pointer :: values_tmp, variances_tmp
      type(silja_time), dimension(:), allocatable :: times_tmp
      type(silja_interval), dimension(:), allocatable :: durations_tmp
      integer, parameter :: tmp_arr_size = 20
      
      prev_time = time_missing
      prev_id = ''
      nread = 0
      itime = 0
      values_tmp => fu_work_array()
      variances_tmp => fu_work_array()
      allocate(times_tmp(tmp_arr_size), durations_tmp(tmp_arr_size), stat=iostat)
      if (iostat /= 0) then
        call set_error('Allocate failed', 'timeseries_from_file')
        return
      end if

      !read(file_unit,fmt=*) line

      do
        call next_line_from_input_file(file_unit, line, eof)
        if (eof) iostat = -1
        !read(unit, fmt='(A)', iostat=iostat) line
        if (iostat > 0) then
          call set_error('Failed to read record', 'timeseries_from_file')
          return
        end if
        store_value = .false.
        flush_values = .false.
        if (iostat == 0) then
          read(line, fmt=*, iostat=iostat) &
               & id, year, month, day, hour, duration_hours, value, stdev
          variance = stdev**2
          if (iostat /= 0) then
            call set_error('Failed to parse record: ' // trim(line), 'timeseries_from_file')
            return
          end if

          duration = fu_set_interval_sec(duration_hours*3600.0)
          time = fu_set_time_utc(year, month, day, hour, 0, 0.0)

          if (duration == zero_interval) then
            ! "instant" observation. To avoid needing special treatment later, we'll set
            ! the duration to 20 sec, but no shifting due to averaging, since there wasn't any.
            duration = fu_set_interval_sec(20.0)
          end if
          store_value = .true.

          if (time < obs_start) then 
            store_value = .false.
          else
            store_value = .true.
          end if
          if (time > obs_start + obs_period) then
            flush_values = .true.
            store_value = .false.
          end if

          if (defined(prev_time)) then
            if (prev_id == id .and. prev_time >= time) then
              call set_error('Cannot read: timeseries not sorted', 'timeseries_from_file')
              return
            end if
          end if
        else ! eof
          flush_values = itime > 0
        end if ! read ok
        
        flush_values = flush_values .and. itime > 0

        if (flush_values) then
          !call msg(id)
          call searchStationWithID(prev_id, station_list, nstations, station, found)
          if (.not. found) then
            call msg_warning('Station ' // trim(prev_id) // ' not found')
          else
            nread = nread + 1
            if (nread > size(dose_rate_addition)) then
              call set_error('Too many observations', 'timeseries_from_file')
              return
            end if
            observations_dose_rate(nread) = fu_init_dose_rate_obs(times_tmp, &
                                                         & durations_tmp, &
                                                         & values_tmp, &
                                                         & variances_tmp, &
                                                         & itime, &
                                                         & variable_variance, &
                                                         & transport_species, &
                                                         & obs_unit, &
                                                         & station, &
                                                         & level_missing, &
                                                         & dispersion_vertical, &
                                                         & 'dose_rate')
            dose_rate_addition(nread) = fu_init_dose_rate_addition(ifGroundToo, &
                                                                 & observations_dose_rate(nread), &
                                                                 & dispersion_vertical, &
                                                                 & transport_species)
            
          end if ! not found
          itime = 0
        end if ! have new series

        if (iostat < 0) exit
        
        if (store_value) then
          itime = itime + 1
          values_tmp(itime) = value
          durations_tmp(itime) = duration
          times_tmp(itime) = time
          variances_tmp(itime) = variance
        end if
        prev_time = time
        prev_id = id
      end do ! loop over lines
      
      deallocate(times_tmp, durations_tmp)
      call free_work_array(values_tmp)
      call free_work_array(variances_tmp)
    end subroutine timeseries_from_file_dose_rate

    !============================================================

    subroutine add_obs_items_from_list(obs_items_fixed, num_obs_items_fixed, &
                                 & list_file_path, p_obs_items, num_obs_items)
      implicit none
      character(len=fnlen), dimension(:), intent(in) :: obs_items_fixed
      integer, intent(in) :: num_obs_items_fixed
      character(len=*), intent(in) :: list_file_path
      character(len=fnlen), dimension(:), pointer :: p_obs_items
      integer, intent(out) :: num_obs_items
      
      type(Tsilam_namelist), pointer :: nlptr
      integer :: iostat, file_unit, ind_obs_item
      type(Tsilam_nl_item_ptr), dimension(:), pointer :: p_items
      character, parameter :: sub_name = 'obs_items_from_list'

      call msg('Reading observation defs from ' // trim(list_file_path))
      file_unit = fu_next_free_unit()
      open(file=list_file_path, unit=file_unit, iostat=iostat)
      if (fu_fails(iostat == 0, 'Failed to open:'//trim(list_file_path), sub_name)) return
      nlptr => fu_read_namelist(file_unit, .false.)
      close(file_unit)
      if (fu_fails(associated(nlptr), 'No namelist read from obs_list', sub_name)) return
      call get_items(nlptr, 'obs_data', p_items, num_obs_items)
      if (fu_fails(num_obs_items > 0, 'No obs_items in obs_list', sub_name)) return
      
      allocate(p_obs_items(num_obs_items_fixed + num_obs_items), stat=iostat)
      if (fu_fails(iostat == 0, 'Allocate failed', sub_name)) return
      do ind_obs_item = 1, num_obs_items_fixed
        p_obs_items(ind_obs_item) = obs_items_fixed(ind_obs_item)
      end do
      do ind_obs_item = 1, num_obs_items
        p_obs_items(ind_obs_item+num_obs_items_fixed) = fu_content(p_items(ind_obs_item))
      end do
      
      call destroy_namelist(nlptr)

    end subroutine add_obs_items_from_list

  end subroutine read_observations
  
  !****************************************************************************************
  !****************************************************************************************
  !****************************************************************************************
    
  subroutine verify_in(value, allowed_values, label, owner)
    implicit none
    character(len=*), intent(in) :: value
    character(len=*), dimension(:), intent(in) :: allowed_values
    character(len=*), intent(in) :: label, owner

    integer :: ii

    if (.not. any(allowed_values == value)) then
      call msg('Allowed values for ' // trim(label))
      do ii = 1, size(allowed_values)
        call msg(allowed_values(ii))
      end do
      call msg('Value given: ' // trim(value))
      call set_error('Invalid value for ' // trim(label), owner)
      return
    end if
  end subroutine verify_in

  !*****************************************************************************************
  
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
      call set_error('No files found for time:' // fu_str(now), &
                   & 'expand_template')
      return
    end if
    if (size(filenames) > 1) then
      call set_error('Multiple files found for time:' // fu_str(now), &
                   & 'expand_template')
    end if
    if (error) return

    filename = fu_process_filepath(filenames(1)%sp, convert_slashes=.true., must_exist=.true.)
    deallocate(filenames)
  end subroutine expand_template


  !**********************************************************************************
  !
  ! InjectAll and observeAll
  !
  !**********************************************************************************
  ! These subroutines are basically similar. In forward run, we loop through all observations
  ! (stored in the pointers structure) and ask them to observe the cloud. The observed value is
  ! saved, and in adjoint run, inject all will inject the difference between the model and 
  ! the observation.

  subroutine observeAll(pointers, cloud, metBuf, disp_buf_ptr, chemRules, timestep, now, eruption_height)
    implicit none
    type(observationPointers), intent(inout) :: pointers
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep
    type(silam_pollution_cloud), pointer :: cloud
    type(TChem_rules), intent(in) :: chemRules
    type(TField_buffer), intent(in) :: metBuf
    type(TField_buffer), intent(in) :: disp_buf_ptr
    real, optional :: eruption_height
    ! local 

    type(THorizInterpStruct), pointer :: p_met_disp_interp_horiz
    logical :: if_met_disp_interp_horiz
    type(TVertInterpStruct), pointer :: p_met_disp_interp_vert
    logical :: if_met_disp_interp_vert
    
    type(TchemicalRunSetup), pointer :: chemRunSetup
    type(Toptical_density_rules), pointer :: opticalRules
    type(Tmass_map), pointer :: mapConc, mapWetdep, mapDrydep, map_px, map_py
    type(field_4d_data_ptr), pointer :: p_rel_hum, p_tempr, p_press, p_hgt
    type(field_2d_data_ptr), pointer :: p_surf_press
    integer :: i, iQ, ind_tempr, ind_rh, ind_surf_press, ind_hgt, ind_press
    
    type(TchemicalRunSetup), pointer :: chem_run_setup
    
    ! To observe mean-grid-cell concentrations, use nearest_point
    ! To observe from trapezoid with mass centre as a proxy for within-cell interpolation, use toMassCentreLinear
    ! 
    integer, parameter :: observe_interpolation = nearest_point ! toMassCentreLinear  ! nearest_point
    character(len=*), parameter :: sub_name = 'observeAll'

!    select case(observe_interpolation)
!      case(nearest_point)
!        call msg('Observing cell-mean concentrations')
!      case(toMassCentreLinear)
!        call msg('Observing profile built by mass-centre and interpolate to site location')
!      case default
!        call set_error('Unknown observation interpolation method','observAll')
!        return
!    end select
    
    
    p_met_disp_interp_horiz => fu_interpCoefMeteo2DispHoriz(cloud)
    p_met_disp_interp_vert => fu_interpCoefMeteo2DispVert(cloud)
    if_met_disp_interp_horiz = fu_ifMeteo2DispHorizInterp(cloud)
    if_met_disp_interp_vert = fu_ifMeteo2DispVertInterp(cloud)    

    mapConc => fu_concMM_ptr(cloud)
    map_PX => fu_advection_moment_X_MM_ptr(cloud)
    map_PY => fu_advection_moment_Y_MM_ptr(cloud)
    mapWetdep => fu_wetdepMM_ptr(cloud)
    mapDrydep => fu_drydepMM_ptr(cloud)
    chemRunSetup => fu_ChemRunSetup(chemRules)
    opticalRules => fu_optical_rules(chemRules)

    if (pointers%nVerticalObservations > 0) then
      ind_rh = fu_index(metbuf, relative_humidity_flag)
      if (fu_fails(ind_rh > 0, 'No relative humidity', sub_name)) return
      p_rel_hum => metbuf%p4d(ind_rh)
      ind_tempr = fu_index(metbuf, temperature_flag)
      if (fu_fails(ind_tempr > 0, 'No temperature', sub_name)) return
      p_tempr => metbuf%p4d(ind_tempr)
      ind_press = fu_index(metbuf, pressure_flag)
      if (fu_fails(ind_press > 0, 'No pressure', sub_name)) return
      p_press => metbuf%p4d(ind_press)
      ind_hgt = fu_index(metbuf, height_flag)
      if (fu_fails(ind_hgt > 0, 'No height', sub_name)) return
      p_hgt => metbuf%p4d(ind_hgt)
      ind_surf_press = fu_index(metbuf, surface_pressure_flag)
      if (fu_fails(ind_surf_press > 0, 'No surface pressure', sub_name)) return
      p_surf_press => metbuf%p2d(ind_surf_press)

      ! A quick hack to avoid hard crash if assimilation is attempted in 3D-Var before first
      ! model integration - meteo data could not be ready yet.
      if (.not. associated(p_rel_hum%past%p2d(1)%ptr)) then
        call set_error('RH pointer not associated', sub_name)
        return
      end if
    end if

    chem_run_setup => fu_chemRunSetup(chemrules)

!!$    call test_vert_obs(mapConc, & 
!!$                     & fu_advection_moment_X_MM_ptr(cloud), &
!!$                     & fu_advection_moment_Y_MM_ptr(cloud), &
!!$                     & fu_advection_moment_Z_MM_ptr(cloud), &
!!$                     & now, timestep, &
!!$                     & p_rel_hum, p_tempr, &
!!$                     & p_met_disp_interp_horiz, if_met_disp_interp_horiz, &
!!$                     & p_met_disp_interp_vert, if_met_disp_interp_vert, &
!!$                     & metbuf%weight_past, fu_optical_rules(chemrules), fu_species_optical(cloud))
!!$    stop

!    call msg("observe_in_situ before parallel called nthreads", omp_get_num_threads(), omp_get_thread_num())
    
    
    !$OMP PARALLEL &
#ifdef DEBUG_OBS    
    !$OMP   & if (.FALSE.)  &
#endif    
    !$OMP  & default (none)   private(i) & 
    !$OMP &  shared(pointers,  mapConc, map_px, map_py, mapWetdep, mapDrydep, now, timestep, p_rel_hum, p_tempr, &
    !$OMP &  p_press, p_hgt, p_surf_press, p_met_disp_interp_horiz, &
    !$OMP if_met_disp_interp_horiz, p_met_disp_interp_vert, chem_run_setup, &
    !$OMP &  if_met_disp_interp_vert, metbuf, disp_buf_ptr, chemrules, cloud, eruption_height) 
!    call msg("observe_in_situ called nthreads", omp_get_num_threads(), omp_get_thread_num())
    !$OMP DO
    do i = 1, pointers%nInSituObservations
      call observe_in_situ(pointers%observationsInSitu(i), mapConc, map_px, map_py, disp_buf_ptr, &
                          & observe_interpolation, now, timestep)
    end do
    !$OMP END DO
    
    !$OMP DO
    do i = 1, pointers%nDoseRateObservations
       call observe_dose_rate(pointers%DoseRateAddition(i), pointers%observationsDoseRate(i), &
                            & mapConc, mapWetdep, mapDrydep, now, timestep)
    end do
    !$OMP END DO

    !$OMP DO
    do i = 1, pointers%nVerticalObservations
      if (pointers%observationsVertical(i)%is_lidar) then
        cycle
        !call observe_lidar(pointers%observationsVertical(i), mapConc, now, timestep, &
        !     & p_rel_hum, p_tempr, p_press, p_hgt, p_surf_press, &
        !     & p_met_disp_interp_horiz, if_met_disp_interp_horiz, &
        !     & p_met_disp_interp_vert, if_met_disp_interp_vert, &
        !     & metbuf%weight_past, fu_optical_rules(chemrules), &
        !     & fu_nbr_of_species_optical(cloud), chem_run_setup)
      else
        call observe_vertical(pointers%observationsVertical(i), mapConc, now, timestep, &
             & p_rel_hum, p_tempr, p_press, p_hgt, p_surf_press, &
             & p_met_disp_interp_horiz, if_met_disp_interp_horiz, &
             & p_met_disp_interp_vert, if_met_disp_interp_vert, &
             & metbuf%weight_past, fu_optical_rules(chemrules), fu_nbr_of_species_optical(cloud))
      end if
    end do
    !$OMP END DO
                                                                                                                                                                                                                                                                                                               
    if (present(eruption_height) .and. (.not. eruption_height == real_missing)) then
      !$OMP DO
      do i = 1, pointers%nEruptionObservations
        call observe_eruption(pointers%observationsEruption(i), eruption_height)
      end do
      !$omp End Do
    end if
    !$omp End Parallel
    
  End subroutine observeAll


  !**************************************************************************************
  
  subroutine injectAll(pointers, cloud, metBuf, chemRules, dynRules, &
                     & timestep, now, injectMap)
    use advection_eulerian, only: fu_if_bulk_eulerian_advection
    implicit none
    type(observationPointers), intent(inout) :: pointers
    type(silam_pollution_cloud), pointer :: cloud
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep
    type(TChem_rules), intent(in) :: chemRules
    type(TDynamics_rules), intent(in) :: dynRules
    type(TField_buffer), intent(in) :: metBuf
    ! A kludgy solution to support 3D-var: instead of cloud's mapConc, use this. The
    ! reason is to allow observed species to be a superset of control species (PM
    ! assimilation). Moments are not dealt with: currently the observations don't touch
    ! them. If this is changed, need to add the moment maps!
    type(Tmass_map), pointer, optional :: injectMap

    !local
    type(THorizInterpStruct), pointer :: p_met_disp_interp_horiz
    logical :: if_met_disp_interp_horiz
    type(TVertInterpStruct), pointer :: p_met_disp_interp_vert
    logical :: if_met_disp_interp_vert
    type(TChemicalRunSetup), pointer :: chemRunSetup
    type(Toptical_density_rules), pointer :: opticalRules
    type(Tmass_map), pointer :: mapConc, mapMomX, mapMomY, mapMomZ, mapWetdep, mapDrydep
        integer :: i, iQ , ind_tempr, ind_rh
    type(TchemicalRunSetup), pointer :: chem_run_setup
    type(field_4d_data_ptr), pointer :: p_rel_hum, p_tempr
    integer :: ithread, nthreads
    
    ! To inject into centre of the cell, use nearest_point
    ! To inject into actual site location, use toMassCentreLinear
    ! 
    integer, parameter :: inject_where = nearest_point ! toMassCentreLinear  ! nearest_point
    character(len=*), parameter :: sub_name = 'injectAll'

!    select case(inject_where)
!      case(nearest_point)
!        call msg('Injecting discrepancy into cell centre')
!      case(toMassCentreLinear)
!        call msg('Injecting discrepancy at the site location')
!      case default
!        call set_error('Unknown injection location','observAll')
!        return
!    end select

    p_met_disp_interp_horiz => fu_interpCoefMeteo2DispHoriz(cloud)
    p_met_disp_interp_vert => fu_interpCoefMeteo2DispVert(cloud)
    if_met_disp_interp_horiz = fu_ifMeteo2DispHorizInterp(cloud)
    if_met_disp_interp_vert = fu_ifMeteo2DispVertInterp(cloud)    

    chemRunSetup => fu_ChemRunSetup(chemRules)
    opticalRules => fu_optical_rules(chemRules)
    if (present(injectMap)) then
      mapConc => injectMap
    else
      mapConc => fu_concMM_ptr(cloud)
    end if
    mapWetdep => fu_wetdepMM_ptr(cloud)
    mapDrydep => fu_drydepMM_ptr(cloud)
    mapMomX => fu_advection_moment_X_MM_ptr(cloud)
    mapMomY => fu_advection_moment_Y_MM_ptr(cloud)
    mapMomZ => fu_advection_moment_Z_MM_ptr(cloud)
    if (pointers%nVerticalObservations > 0) then
      ind_rh = fu_index(metbuf, relative_humidity_flag)
      if (fu_fails(ind_rh /= int_missing, 'No relative humidity', sub_name)) return
      ind_tempr = fu_index(metbuf, temperature_flag)
      if (fu_fails(ind_tempr /= int_missing, 'No temperature', sub_name)) return
      p_rel_hum => metbuf%p4d(ind_rh)
      p_tempr => metbuf%p4d(ind_tempr)
    end if
    chem_run_setup => fu_chemRunSetup(chemrules)


    if (debug_level > 0) then
      call msg('sum before obs. injection...', &
           &   sum(mapConc%arm(:, :, 1:nz_dispersion,1:nx_dispersion,1:ny_dispersion)))
      call msg('Max px', maxval(mapMomX%arm))
      call msg('Max py', maxval(mapMomY%arm))

    end if
    
    nthreads = 1  ! for the non-OMP case
    ithread = 0
    
    !$OMP PARALLEL default (none) private(ithread,nthreads,i) &
    !$OMP shared(pointers, mapConc, mapMomX, mapMomY, mapMomZ,now, timestep)
    !$ ithread  = omp_get_thread_num()
    !$ nthreads = omp_get_num_threads()
    
    !! THIS IS NOT OMP DO!!!!!
    do i = 1, pointers%nInSituObservations
      call inject_in_situ(pointers%observationsInSitu(i), mapConc, mapMomX, mapMomY, mapMomZ, &
                & now, timestep, ithread, nthreads, inject_where)
    end do
    !$OMP END PARALLEL 
    
    do i = 1, pointers%nDoseRateObservations
       call inject_dose_rate(pointers%DoseRateAddition(i), pointers%observationsDoseRate(i), &
                           & mapConc, mapMomX, mapMomY, mapMomZ, mapWetdep, mapDrydep, now, timestep)
    end do

    do i = 1, pointers%nVerticalObservations
       call inject(pointers%observationsVertical(i), mapConc, mapMomX, mapMomY, mapMomZ, &
                 & now, timestep, &
                 & p_rel_hum, p_tempr, &
                 & p_met_disp_interp_horiz, if_met_disp_interp_horiz, &
                 & p_met_disp_interp_vert, if_met_disp_interp_vert, &
                 & metbuf%weight_past, fu_optical_rules(chemrules), fu_nbr_of_species_optical(cloud))
    end do
    if (debug_level > 0) then
      do iQ = 1, mapConc%nSpecies
        call msg('Species: ' + fu_str(mapConc%species(iQ)))
        do i = 1, mapConc%nSrc
          call msg('sum after obs. injection, iSrc', i, sum(mapConc%arm(iQ,i,1:nz_dispersion,1:nx_dispersion,1:ny_dispersion)))
        end do
      end do
      call msg('Max px', maxval(abs(mapMomX%arm)))
      call msg('Max py', maxval(abs(mapMomY%arm)))
    end if
  end subroutine injectAll

  !***********************************************************************************

  subroutine set_observations(station_path, obs_list_path, obs_start, obs_period, obs_items, nObsItems, &
                            & transport_species, optical_species, observations)
    !
    !Gets needed observations from files ans stores them to "observations"
    implicit none
    character(len=*), intent(in) :: station_path, obs_list_path
    character(len=*), dimension(:), pointer :: obs_items
    type(silja_time), intent(in) :: obs_start
    type(silja_interval), intent(in) :: obs_period
    integer, intent(inout) :: nObsItems
    type(silam_species), dimension(:), intent(in) :: transport_species
    type(silam_species), dimension(:), pointer :: optical_species
    type(observationPointers), intent(out) :: observations
    
    !character(len=worksize) :: nl_content
    integer :: file_unit, status, ii, n_obs_items, n_observations, n_stations
    type(Tsilam_namelist), pointer :: nl_stations
    type(observationStationPtr), dimension(:), pointer :: station_list
    integer :: n_transp_species, n_opt_species
    type(grads_template) :: station_template
    character(len=fnlen) :: station_list_file
    
    if (associated(optical_species)) then
      n_opt_species = size(optical_species)
    else
      n_opt_species = 0
    end if
    
    call msg('Reading station descriptions')

    if (station_path == '') then
      n_stations = 0
      nullify(station_list)
    else
      file_unit = fu_next_free_unit()
      call decode_template_string(station_path, station_template)
      if (error) return
      call expand_template(station_template, obs_start, station_list_file)
      if (error) return
      open(file_unit, file=station_list_file, status='old', action='read', iostat=status)
      if (status /= 0) then
        call set_error('Failed to open ' // trim(station_list_file), 'set_observations')
        return
      end if
      nl_stations => fu_read_namelist(file_unit, .false.)
      close(file_unit)
      call stations_from_namelist_ptr(nl_stations, station_list, n_stations, dispersion_grid) ! Sic!!
                    ! Set sstation indices for current ssubdomain
      call destroy_namelist(nl_stations)
      if (error) return
    end if
        
    n_transp_species = size(transport_species)
    
    if (nObsItems > 0 .or. obs_list_path /= '') then
      call read_observations(obs_start, obs_period, obs_list_path, obs_items, nObsItems, &
                           & station_list, n_stations, &
                           & observations%observationsInSitu, observations%nInSituObservations, &
                           & observations%observationsVertical, observations%nVerticalObservations, &
                           & observations%observationsEruption, observations%nEruptionObservations, &
                           & observations%DoseRateAddition, observations%nDoseRateObservations, &
                           & transport_species, optical_species, n_transp_species, n_opt_species)
      if (error) return
      call msg('N. of in situ observations loaded:', observations%nInSituObservations)
      call msg('N. of column observations loaded:', observations%nVerticalObservations)
      call msg('N. of eruption observations loaded:', observations%nEruptionObservations)
      call msg('N. of dose rate observations loaded:', observations%nDoseRateObservations)
    else
      call msg('No observations given')
    end if
    
    observations%obs_size = fu_count_obs_size(observations)
    allocate(observations%obs_values(observations%obs_size), &
           & observations%obs_variance(observations%obs_size), stat=status)
    if (fu_fails(status == 0, 'Allocate failed', 'set_observations')) return
    call collect_obs_data(observations, observations%obs_values, observations%obs_size)
    call collect_variance(observations, observations%obs_variance, observations%obs_size)
    observations%hasObservations = .true.

    !!!!!!!!!!!!!!!!!!!!!!!!!!
   ! print *, 'Max obs value: ', maxval(observations%obs_values)
   ! print *, 'Number of obs: ', observations%obs_size
    call msg('Max obs value: ', maxval(observations%obs_values))
    call msg('Ave obs value: ', sum(observations%obs_values)/size(observations%obs_values))
    call msg('Number of obs: ', observations%obs_size)
   ! call ooops("End obs read")
    !!!!!!!!!!!!!!!!!!!!!!!!

  contains

    integer function fu_count_obs_size(obs) result(how_many)
      implicit none
      type(observationPointers), intent(in) :: obs
      
      integer :: ii, subtotal

      how_many = 0
      do ii = 1, obs%nInSituObservations
        how_many = how_many + fu_size(obs%observationsInSitu(ii))
      end do
      do ii = 1, obs%nEruptionObservations
         how_many = how_many + fu_size(obs%observationsEruption(ii))
      end do
      do ii = 1, obs%nVerticalObservations
        how_many = how_many + fu_size(obs%observationsVertical(ii))
      end do
      do ii = 1, obs%nDoseRateObservations
         how_many = how_many + fu_size(obs%observationsDoseRate(ii))
      end do

    end function fu_count_obs_size

  end subroutine set_observations


  !************************************************************************************

  subroutine collect_model_data(observations, values, n_values_total)
    !
    ! Picks model-subdomain-observed data from observations structure into values array
    ! makes exchange to form whole-MPI values array
    !
    implicit none
    type(observationPointers), intent(in) :: observations
    real, dimension(:), intent(out) :: values
    integer, intent(out) :: n_values_total

    integer :: ind_obs, ind_start, n_values
    real, dimension(:), pointer :: wrk
    logical :: ok
    
    ind_start = 1
    n_values_total = 0

    call msg('Collecting model data')

    do ind_obs = 1, observations%nInSituObservations
      n_values = fu_size(observations%observationsInSitu(ind_obs))
      call get_data(observations%observationsInSitu(ind_obs), values_mdl=values(ind_start:))
      ind_start = ind_start + n_values
      n_values_total = n_values_total + n_values
      if (fu_fails(ind_start <= size(values), 'values too small', 'get_model_data')) return
    end do

    do ind_obs = 1, observations%nEruptionObservations 
      n_values = fu_size(observations%observationsEruption(ind_obs))  
      call get_data(observations%observationsEruption(ind_obs), values_mdl=values(ind_start:))
      ind_start = ind_start + n_values
      n_values_total = n_values_total + n_values
      if (fu_fails(ind_start <= size(values), 'values too small', 'get_model_data')) return
    end do
    
    do ind_obs = 1, observations%nDoseRateObservations
      n_values = fu_size(observations%observationsDoseRate(ind_obs))
      call get_data(observations%observationsDoseRate(ind_obs), values_mdl=values(ind_start:))
      ind_start = ind_start + n_values
      n_values_total = n_values_total + n_values
      if (fu_fails(ind_start <= size(values), 'values too small', 'get_model_data')) return
    end do

    do ind_obs = 1, observations%nVerticalObservations
      n_values = fu_size(observations%observationsVertical(ind_obs))
!      call msg('Calling get data for vertical observation', ind_obs)
      call get_data(observations%observationsVertical(ind_obs), values_mdl=values(ind_start:))
!      call msg('ave values', sum(values(ind_start:ind_start + n_values))/size(values(ind_start:ind_start + n_values)))
!      call msg('max values', maxval(values(ind_start:ind_start + n_values)))
!      call msg('ind_start, n_values', ind_start, n_values)
      ind_start = ind_start + n_values
      n_values_total = n_values_total + n_values
      if (fu_fails(ind_start <= size(values), 'values too small', 'get_model_data')) return
    end do

!    call msg('n values total', n_values_total)
!    call msg('size values', size(values))


    if (smpi_global_tasks > 1) then
#ifdef DEBUG
      call msg("LOCAL Sum of observarions from model", sum(values(1:n_values_total)))
      call msg("LOCAL values(1:10)", values(1:10))
#endif
      ! Synchronize model observationss
      wrk=>fu_work_array(n_values_total)
      wrk(1:n_values_total)=values(1:n_values_total)
!      print *, "n_values_total", n_values_total
!      print *, wrk(1:10)
!      print *, values(1:10)
      call smpi_reduce_add(wrk(1:n_values_total), values(1:n_values_total), 0, smpi_adv_comm, ok)
      call free_work_array(wrk)
      if (.not. ok) call set_error("failed MPI COMM", "collect_model_data")
      
      ! call ooops("Reducing obs")

    endif
#ifdef DEBUG
    call msg("Sum of observarions from model , n_values_total", sum(values(1:n_values_total)), n_values_total)
    call msg("values(1:10)", values(1:10))
#endif
  end subroutine collect_model_data


  !************************************************************************************

  subroutine collect_obs_data(observations, values, n_values_total)
    implicit none
    type(observationPointers), intent(in) :: observations
    real, dimension(:), intent(out) :: values
    integer, intent(out) :: n_values_total

    integer :: ind_obs, ind_start, n_values
    
    ind_start = 1
    n_values_total = 0

    do ind_obs = 1, observations%nInSituObservations
      n_values = fu_size(observations%observationsInSitu(ind_obs))
      call get_data(observations%observationsInSitu(ind_obs), values_obs=values(ind_start:))
      ind_start = ind_start + n_values
      n_values_total = n_values_total + n_values
    end do

    do ind_obs = 1, observations%nEruptionObservations
      n_values = fu_size(observations%observationsEruption(ind_obs))
      call get_data(observations%observationsEruption(ind_obs), values_obs=values(ind_start:))
      ind_start = ind_start + n_values
      n_values_total = n_values_total + n_values
    end do
    
    do ind_obs = 1, observations%nDoseRateObservations
      n_values = fu_size(observations%observationsDoseRate(ind_obs))
      call get_data(observations%observationsDoseRate(ind_obs), values_obs=values(ind_start:))
      ind_start = ind_start + n_values
      n_values_total = n_values_total + n_values
    end do
    
    do ind_obs = 1, observations%nVerticalObservations
      n_values = fu_size(observations%observationsVertical(ind_obs))
      call get_data(observations%observationsVertical(ind_obs), values_obs=values(ind_start:))
      ind_start = ind_start + n_values
      n_values_total = n_values_total + n_values
    end do
    
  end subroutine collect_obs_data
  

  !************************************************************************************

  function fu_obs_data(observations) result(obs_data)
    implicit none
    type(observationPointers), intent(in) :: observations
    real, dimension(:), pointer :: obs_data
    obs_data => observations%obs_values
  end function fu_obs_data


  !************************************************************************************

  function fu_obs_variance(observations) result(obs_variance)
    implicit none
    type(observationPointers), intent(in) :: observations
    real, dimension(:), pointer :: obs_variance
    obs_variance => observations%obs_variance
  end function fu_obs_variance


  !************************************************************************************

  subroutine get_obs_data(observations, values, n_values_total)
    implicit none
    type(observationPointers), intent(in) :: observations
    real, dimension(:), pointer :: values
    integer, intent(out) :: n_values_total
    
    n_values_total = 0

    values => observations%obs_values
    n_values_total = observations%obs_size
#ifdef DEBUG
    call msg("get_obs_data, sum, ndata", sum(values(1:n_values_total)), real(n_values_total))
#endif
    
  end subroutine get_obs_data


  !************************************************************************************

  subroutine set_mdl_data(observations, values)
    implicit none
    type(observationPointers), intent(inout) :: observations
    real, dimension(:), intent(in) :: values

    integer :: ind_obs, ind_start, n_values
    
    ind_start = 1

    do ind_obs = 1, observations%nInSituObservations
      n_values = fu_size(observations%observationsInSitu(ind_obs))
      call set_data(observations%observationsInSitu(ind_obs), values_mdl=values(ind_start:))
      ind_start = ind_start + n_values
    end do

    do ind_obs = 1, observations%nEruptionObservations
      n_values = fu_size(observations%observationsEruption(ind_obs))
      call set_data(observations%observationsEruption(ind_obs), values_mdl=values(ind_start:))
      ind_start = ind_start + n_values
    end do

    
    do ind_obs = 1, observations%nDoseRateObservations
      n_values = fu_size(observations%observationsDoseRate(ind_obs))
      call set_data(observations%observationsDoseRate(ind_obs), values_mdl=values(ind_start:))
      ind_start = ind_start + n_values
    end do
    
    do ind_obs = 1, observations%nVerticalObservations
      n_values = fu_size(observations%observationsVertical(ind_obs))
      call set_data(observations%observationsVertical(ind_obs), values_mdl=values(ind_start:))
      ind_start = ind_start + n_values
    end do
    
  end subroutine set_mdl_data


  !************************************************************************************

  subroutine collect_variance(observations, values, n_values_total)
    implicit none
    type(observationPointers), intent(in) :: observations
    real, dimension(:), intent(out) :: values
    integer, intent(out) :: n_values_total

    integer :: ind_obs, ind_start, n_values
    
    ind_start = 1
    n_values_total = 0

    do ind_obs = 1, observations%nInSituObservations
      n_values = fu_size(observations%observationsInSitu(ind_obs))
      call get_data(observations%observationsInSitu(ind_obs), variance=values(ind_start:))
      ind_start = ind_start + n_values
      n_values_total = n_values_total + n_values
    end do
    
    do ind_obs = 1, observations%nEruptionObservations
      n_values = fu_size(observations%observationsEruption(ind_obs))
      call get_data(observations%observationsEruption(ind_obs), variance=values(ind_start:))
      ind_start = ind_start + n_values
      n_values_total = n_values_total + n_values
    end do

    do ind_obs = 1, observations%nDoseRateObservations
      n_values = fu_size(observations%observationsDoseRate(ind_obs))
      call get_data(observations%observationsDoseRate(ind_obs), variance=values(ind_start:))
      ind_start = ind_start + n_values
      n_values_total = n_values_total + n_values
    end do
    
    do ind_obs = 1, observations%nVerticalObservations
      n_values = fu_size(observations%observationsVertical(ind_obs))
      call get_data(observations%observationsVertical(ind_obs), variance=values(ind_start:))
      ind_start = ind_start + n_values
      n_values_total = n_values_total + n_values
    end do

    
  end subroutine collect_variance

  !************************************************************************************

  subroutine get_localisation_all(obs_ptr, locations)
    implicit none
    type(observationPointers), intent(in) :: obs_ptr
    real, dimension(:,:) :: locations
    
    integer :: ind_obs, ind_start, num_values

    if (fu_fails(size(locations, 1) > 1, 'locations too small 1st dim', 'get_localisation')) return
    if (fu_fails(size(locations, 2) >= obs_ptr%obs_size, 'locations too small 2st dim', &
               & 'get_localisation')) return

    ind_start = 1
    num_values = 0
    do ind_obs = 1, obs_ptr%nInSituObservations
      num_values = fu_size(obs_ptr%observationsInSitu(ind_obs))
      call get_localisation(obs_ptr%observationsInSitu(ind_obs), locations(:,ind_start:))
      ind_start = ind_start + num_values
    end do
    
    do ind_obs = 1, obs_ptr%nEruptionObservations
      num_values = fu_size(obs_ptr%observationsEruption(ind_obs))
      call get_localisation(obs_ptr%observationsEruption(ind_obs), locations(:,ind_start:))
      ind_start = ind_start + num_values
    end do

    do ind_obs = 1, obs_ptr%nDoseRateObservations
      num_values = fu_size(obs_ptr%observationsDoseRate(ind_obs))
      call get_localisation(obs_ptr%observationsDoseRate(ind_obs), locations(:,ind_start:))
      ind_start = ind_start + num_values
    end do

    do ind_obs = 1, obs_ptr%nVerticalObservations
      num_values = fu_size(obs_ptr%observationsVertical(ind_obs))
      call get_localisation(obs_ptr%observationsVertical(ind_obs), locations(:,ind_start:))
      ind_start = ind_start + num_values
    end do
    
  end subroutine get_localisation_all

  !************************************************************************************

  subroutine reset_all(observations, ifNullifyModel)
    !
    ! Resets the observation counters to the beginning, nullifies the model values on request
    !
    implicit none
    type(observationPointers), intent(inout) :: observations
    logical, intent(in) :: ifNullifyModel

    integer :: ii

    do ii = 1, observations%nInSituObservations
      call restart_in_situ(observations%observationsInSitu(ii), ifNullifyModel)
    end do
    
    do ii = 1, observations%nDoseRateObservations
      call restart_dose_rate(observations%observationsDoseRate(ii), observations%DoseRateAddition(ii), ifNullifyModel)
    end do

    do ii = 1, observations%nVerticalObservations
      call restart_vertical(observations%observationsVertical(ii), ifNullifyModel)
    end do

  end subroutine reset_all


  !************************************************************************************

  subroutine destroy_observations(pointers)
    implicit none
    type(observationPointers) :: pointers
    
    integer :: i

    do i = 1, pointers%nInSituObservations
       call destroy(pointers%observationsInSitu(i))
    end do
    if (pointers%nInSituObservations > 0) deallocate(pointers%observationsInSitu)
    pointers%nInSituObservations = 0

    do i = 1, pointers%nEruptionObservations
       call destroy(pointers%observationsEruption(i))
    end do
    if (pointers%nEruptionObservations > 0) deallocate(pointers%observationsEruption)
    pointers%nInSituObservations = 0
    
    do i = 1, pointers%nDoseRateObservations
       call destroy(pointers%observationsDoseRate(i))
    end do
    if (pointers%nDoseRateObservations > 0) deallocate(pointers%observationsDoseRate)
    pointers%nDoseRateObservations = 0
    
    do i = 1, pointers%nVerticalObservations
       call destroy(pointers%observationsVertical(i))
    end do
    if (pointers%nVerticalObservations > 0) deallocate(pointers%observationsVertical)
    pointers%nVerticalObservations = 0

    if (pointers%obs_size > 0) deallocate(pointers%obs_values, pointers%obs_variance)
    
    pointers%obs_size = 0
    pointers%hasObservations = .false.
         
  end subroutine destroy_observations

  !**********************************************************************************
  !
  ! Auxiliary routines
  !
  !**********************************************************************************

  subroutine searchStationWithID(id, station_list, nstations, station, wasFound)
    implicit none

    ! Walk through the station list and return 
    ! the station with a given ID (which is assumed to be unique).

    !integer, intent(in) :: id
    character(len=*), intent(in) :: id
    type(observationStationPtr), dimension(:), intent(in) :: station_list
    integer, intent(in) :: nstations
    type(observationStation), pointer :: station
    logical, intent(out) :: wasFound

    integer :: i
    wasfound = .false.

    do i = 1, nstations
!call msg('Checking:' + station_list(i)%ptr%id, i)
      if (station_list(i)%ptr%id == id) then
        station => station_list(i)%ptr
        wasFound = .true.
        return
      end if
    end do

    wasFound = .false.
    nullify(station)

  end subroutine searchStationWithID

  !************************************************************************************

  subroutine dump_observations(pointers, species_trn, file_name)
    implicit none
    type(observationPointers), intent(in) :: pointers
    character(len=*), intent(in) :: file_name
    type(silam_species), dimension(:), pointer :: species_trn

    integer :: file_unit, ind_obs, iostat
    
    call msg('Observations to file: ' // trim(file_name))
    
    file_unit = fu_next_free_unit()
    
    if(pointers%nInSituObservations > 0)then
      open(file_unit, file=trim(file_name) // '.in_situ', iostat=iostat, action='write')
!      call msg('iostat', iostat)
      if (fu_fails(iostat == 0, 'Failed to open file: ' // trim(file_name), 'dump_observations')) return
      write(file_unit, fmt='(A)') '# In-situ observations # station_id, time, obs_data, model_data, obs_variance'
      do ind_obs = 1, pointers%nInSituObservations
        call obs_to_file(pointers%observationsInSitu(ind_obs), species_trn, file_unit)
        if (error) exit
      end do
      close(file_unit)
    endif
    
    if(pointers%nDoseRateObservations > 0)then
      open(file_unit, file=trim(file_name) // '.dose_rate', iostat=iostat, action='write')
      if (fu_fails(iostat == 0, 'Failed to open file: ' // trim(file_name), 'dump_observations')) return
      do ind_obs = 1, pointers%nDoseRateObservations
        call obs_to_file(pointers%observationsDoseRate(ind_obs), species_trn, file_unit)
        if (error) exit
      end do
      close(file_unit)
    endif
    
    if(pointers%nVerticalObservations > 0)then
      open(file_unit, file=trim(file_name) // '.vertical', iostat=iostat, action='write')
      if (fu_fails(iostat == 0, 'Failed to open file: ' // trim(file_name), 'dump_observations')) return
      do ind_obs = 1, pointers%nVerticalObservations
        call obs_to_file(pointers%observationsVertical(ind_obs), species_trn, file_unit)
        if (error) exit
      end do
      close(file_unit)
    endif
    
  end subroutine dump_observations

  
  !***************************************************************************************
  
  integer function fu_number_of_observed_cells(obs_pointers)
    !
    ! A funny function that says how many grid points are filled with observations.
    ! Used in low_mass_threshold setting routine
    !
    implicit none

    type(observationPointers), intent(in) :: obs_pointers
  
    fu_number_of_observed_cells = obs_pointers%nInSituObservations &
                               & + obs_pointers%nVerticalObservations &
                               & + obs_pointers%nDoseRateObservations
    
  end function fu_number_of_observed_cells
    
  
!!$  subroutine refineStationList(stationListPtr)
!!$    implicit none
!!$    type(DA_listNode), pointer :: stationListPtr
!!$    
!!$    type(DA_listNode), pointer :: auxPtr, auxPtr2, auxPtr3
!!$    integer, dimension(:,:), allocatable :: stationsInCells
!!$    integer :: allocStat
!!$
!!$    allocate(stationsInCells(nx_dispersion, ny_dispersion), stat=allocStat)
!!$    if (allocStat /= 0) then
!!$      call set_error('Allocate failed','refineStationList')
!!$      return
!!$    end if
!!$    
!!$    stationsInCells = 0
!!$    
!!$    
!!$    auxPtr => stationListPtr
!!$    auxPtr3 => null()
!!$    do while(associated(auxPtr))
!!$      
!!$      auxPtr2 => auxPtr%next ! next
!!$      
!!$      if (stationsInCells(auxPtr%station%iX_dispersionGrid, auxPtr%station%iY_dispersionGrid) == 1) then
!!$        ! kill this node
!!$        call deallocatePointInterpStruct(auxPtr%station%interpStruct)
!!$        call msg('Removing station ', auxPtr%station%id)
!!$        deallocate(auxPtr%station)
!!$        ! if it was not the first, move the previous to point to the next one
!!$        if (associated(auxPtr3)) auxPtr3%next => auxPtr2
!!$        if (.not. associated(auxPtr2)) auxPtr3%hasNext = .false.
!!$        deallocate(auxPtr)
!!$      else 
!!$        stationsInCells(auxPtr%station%iX_dispersionGrid, auxPtr%station%iY_dispersionGrid) = 1
!!$        auxPtr3 => auxPtr ! previous
!!$      end if
!!$
!!$      auxPtr => auxPtr2 ! current
!!$    end do
!!$    deallocate(stationsInCells)
!!$    
!!$    
!!$  end subroutine refineStationList

  subroutine test_read_timeseries()
    implicit none
    
    integer, parameter :: num_transp_sp = 2, num_stations = 2
    type(silam_species), dimension(2) :: transp_species

    type(observationStationPtr), dimension(2) :: station_list
    type(inSituObservation), dimension(:), pointer :: in_situ_list, dose_rate_list
    type(t_dose_rate_obs_addition), dimension(:), pointer :: addition_list
    integer :: num_in_situ, num_vert_obs, num_dose_rate, num_eruption_obs
    character(len=fnlen), dimension(:), pointer :: obs_items
    character(len=fnlen) :: obs_list_path
    character(len=*), parameter :: tmpfilename = 'tmp.dat'
    type(silja_grid) :: grid
    integer :: obs_size, ind_obs, num_val, nObsItems
    character(len=*), parameter :: sub_name = 'test_read_timeseries'
    logical :: failure
    real, dimension(:), pointer :: dataptr
    real, parameter :: dummy_obs_val = 1.0
    real :: conv_to_mol
    type(t_vertical_observation), pointer, dimension(:) :: null_vert_obs_lst
    type(silam_species), dimension(:), pointer :: null_species_lst
    type(t_eruptionObservation), pointer, dimension(:) :: null_eruption_obs_lst
    type(silja_time) :: obs_start
    type(silja_interval) :: obs_period

    grid = fu_set_lonlat_grid('testgrid', 20.0, 40.0, .true., &
                            & 100, 100, pole_geographical, 1.0, 1.0)
    if (error) return
    call set_species(transp_species(1), fu_get_material_ptr('SO2'), in_gas_phase)
    call set_species(transp_species(2), fu_get_material_ptr('NO2'), in_gas_phase)
    if (error) return
    
    call set_vertical(fu_set_level(layer_btw_2_height, 0.0, 20.0), dispersion_vertical)
    if (error) return

    nObsItems = 1
    allocate(obs_items(1))
    obs_items(1) = 'cnc g ' // tmpfilename
    conv_to_mol = 1/64.0

    call make_stations(station_list, 2)
    !station_list(1) = fu_initObservationStation('A1', 'station_a1', 30.0, 55.0, 5.0, wholeMPIdispersion_grid)
    !station_list(2) = fu_initObservationStation('A2', 'station_a2', 31.0, 56.0, 5.0, wholeMPIdispersion_grid)
    if(error) return
    
    obs_start = fu_set_time_utc(2012, 6, 1, 0, 0, 0.0)
    obs_period = fu_set_interval_h(36)

    nullify(null_vert_obs_lst, null_species_lst)

    ! what if obs not foud?


    call msg('Expect failure:')
    call read_observations(obs_start, obs_period, obs_list_path, obs_items, nObsItems, &
                         & station_list, size(station_list), in_situ_list, num_in_situ, &
                         & null_vert_obs_lst, num_vert_obs, &
                         & null_eruption_obs_lst, num_eruption_obs, &
                         & addition_list, num_dose_rate,  &
                         & transp_species, null_species_lst, num_transp_sp, 0)
    if (fu_fails(error, 'Shouldn''t work!', 'test_read_timeseries')) continue
    call unset_error('test_read_timeseries')


    call msg('Expect success:')

    call msg('Hourly obs')
    call make_observations(tmpfilename, station_list, obs_start, obs_start+obs_period, &
                         & 'SO2', one_hour, obs_size)
    call read_observations(obs_start, obs_period, obs_list_path, obs_items, nObsItems, &
                         & station_list, size(station_list), in_situ_list, num_in_situ, &
                         & null_vert_obs_lst, num_vert_obs, &
                         & null_eruption_obs_lst, num_eruption_obs, &
                         & addition_list, num_dose_rate, &
                         & transp_species, null_species_lst, num_transp_sp, 0)

    call msg('series read:', num_in_situ)
    if (fu_fails(num_in_situ == 2, 'Failed to read 2 obs', sub_name)) return 

    call check_observations(in_situ_list, obs_size, station_list)

    call msg('Daily obs')
    call make_observations(tmpfilename, station_list, obs_start, obs_start+obs_period, &
                         & 'SO2', one_day, obs_size)
    call read_observations(obs_start, obs_period, obs_list_path, obs_items, nObsItems, &
                         & station_list, size(station_list), in_situ_list, num_in_situ, &
                         & null_vert_obs_lst, num_vert_obs, & 
                         & null_eruption_obs_lst, num_eruption_obs, &
                         & addition_list, num_dose_rate, &
                         & transp_species, null_species_lst, num_transp_sp, 0)
    call msg('series read:', num_in_situ)
    if (fu_fails(num_in_situ == 2, 'Failed to read 2 obs', sub_name)) return 

    call check_observations(in_situ_list, obs_size, station_list)
    
    
    !call cleanup(tmpfilename)

  contains
    
    subroutine check_observations(obs_list, obs_size, station_ptr_list)
      implicit none
      type(insituObservation), dimension(:), intent(in) :: obs_list
      type(observationStationPtr), dimension(:), intent(in) :: station_ptr_list
      integer, intent(in) :: obs_size

      real, dimension(:), pointer :: dataptr
      real, dimension(obs_size) :: expect_values
      integer :: ii

      dataptr => fu_work_array()
      do ind_obs = 1, size(obs_list)
        expect_values = (/(ind_obs + ii - 1, ii = 1, obs_size)/)
        num_val = fu_size(obs_list(ind_obs))
        call msg('num_val, obs_size', num_val, obs_size)
        failure = fu_fails(num_val == obs_size, 'Obs are wrong size', sub_name)
        failure = fu_fails(obs_list(ind_obs)%station%id == station_ptr_list(ind_obs)%ptr%id, &
             & 'wrong id', sub_name)
        call get_data(in_situ_list(ind_obs), values_obs=dataptr)
        if (fu_fails(all(dataptr(1:obs_size) .eps. expect_values), 'Wrong values', sub_name)) then
          call msg('ind_obs:', ind_obs)
          print *, 'Expected:', expect_values(1:obs_size)
          print *, 'Obtained:', dataptr(1:obs_size)
        end if
        
      end do
      
      call free_work_array(dataptr)
      
    end subroutine check_observations

    subroutine make_stations(station_ptr_list, how_many)
      implicit none
      type(observationStationPtr), dimension(:), intent(out) :: station_ptr_list
      integer, intent(in) :: how_many
      
      integer :: ind_station
      type(observationStation), pointer :: stationptr
      character(clen) :: code, name
      
      do ind_station = 1, how_many
        allocate(stationptr)
        write(code, fmt='("A",I0)') ind_station
        write(name, fmt='("station_", I0)') ind_station
        stationptr = fu_initObservationStation(trim(code), trim(name), &
                                             & 30.0 + ind_station, 55.0 + ind_station, 2.0, &
                                             & grid)
        station_ptr_list(ind_station)%ptr => stationptr
      end do
      
    end subroutine make_stations

    subroutine make_observations(filename, station_ptr_list, time_start, time_end, cockt_name, duration, &
                               & num_val)
      implicit none
      character(len=*), intent(in) :: filename, cockt_name
      type(observationStationPtr), dimension(:), intent(in) :: station_ptr_list
      type(silja_time), intent(in) :: time_start, time_end
      !real, intent(in) :: value
      type(silja_interval), intent(in) :: duration
      integer, intent(out) :: num_val

      integer :: unit
      real :: hours
      integer :: ind_station
      type(silja_time) :: now, obs_end
      real :: value

      unit = fu_next_free_unit()
      
      open(unit, action='write', form='formatted', file=filename)

      hours = fu_sec(duration) / 3600.0

      do ind_station = 1, size(station_list)
        now = time_start
        num_val = 0
        do while (now < time_end)
          value = ind_station + num_val
          obs_end = now + duration
          write(unit, fmt='(A, A, I5, I3, I3, I3, 3G12.3)') station_ptr_list(ind_station)%ptr%id, cockt_name, &
               & fu_year(obs_end), fu_mon(obs_end), fu_day(obs_end), fu_hour(obs_end), &
               & hours, value, 1.0
          now = now + duration
          ! The values that end outside assimilation window are not read, hence don't count them in:
          if (obs_end <= time_end) num_val = num_val + 1 
        end do
      end do

      close(unit)
      
    end subroutine make_observations
    
    subroutine cleanup(tmpfilename)
      implicit none
      character(len=*), intent(in) :: tmpfilename
      
      integer :: unit

      unit = fu_next_free_unit()

      open(unit=unit, file=tmpfilename)
      close(unit, status='delete')
    end subroutine cleanup

  end subroutine test_read_timeseries



end module observation_server

