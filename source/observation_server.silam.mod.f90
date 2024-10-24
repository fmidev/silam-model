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
  public dump_observation_stations
  public reset_all

  public collect_model_data
  public get_obs_pointers
  public set_observations
  public fu_number_of_observed_cells

  public test_read_timeseries
  public get_localisation
  
  private collect_obs_data
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
  integer, parameter, public :: m1qn3_flag = 92, l_bfgs_b_flag = 93
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
     
     !! Aggregate values in a single array
     real, dimension(:), allocatable :: obs_values, obs_variance, mdl_values
     ! observations are split to asssimilation and evalution subsets put one after another, in
     ! this order. Hence, dimension = 2
     integer, dimension(2) :: nInSituObsID = 0, nVerticalObsID = 0, nDoseRateObsID = 0, &
                            & nEruptionObsID = 0
     integer, dimension(2) :: obs_size = (/0,0/)
     logical :: hasObservations = .false. 
     logical :: mdlCollected = .false.  !! safeguard against double collecting causes trouble in MPI
                           
  end type ObservationPointers
  public observationPointers

    
contains
 
  !************************************************************************************

  subroutine stations_from_namelist_no_ptr(nlPtr, station_list, nstations)
    implicit none
    type(Tsilam_namelist), pointer :: nlPtr
    type(observationStation), dimension(:), pointer :: station_list
    integer, intent(out) :: nstations

    integer :: nitems, stat, i, iStat
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: items
    character(len=worksize_string) :: workstring
    real :: lat, lon, hgt
    character(len=STATION_LABEL_LENGTH) :: label
    character(len=clen) :: name

    nullify(items)
    
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
      workstring = fu_content(items(i))
      read(workstring, fmt=*, iostat=stat) label, lat, lon, hgt, name
      if (stat /= 0) then
        call msg('Problem with record 1: ' + workstring + ', retrying without station name')
        read(workstring, fmt=*, iostat=stat) label, lat, lon, hgt
        if (stat /= 0) then
          call set_error('No way. Problem parsing station list', 'stations_from_namelist')
          return
        end if
        name = label
      end if
      station_list(iStat) = fu_initObservationStation(label, name, lon, lat, hgt)
      if(defined(station_list(iStat))) iStat = iStat + 1  ! if initialization was successful
    end do
    nstations = iStat - 1 ! the number of actually initialised stations
    
    deallocate(items)

  end subroutine stations_from_namelist_no_ptr

  !************************************************************************************

  subroutine stations_from_namelist_ptr(nlPtr, station_list, nstations)
    !
    ! Reads the stations metadata from namelist
    ! ATTENTION!
    ! Fixed fields order is assumed:
    ! read(workstring, fmt=*) label, lat, lon, hgt, name
    !
    implicit none
    type(Tsilam_namelist), intent(in) :: nlPtr
    type(observationStation), dimension(:), allocatable :: station_list
    integer, intent(out) :: nstations

    integer :: nitems, stat, ii, iStat
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: items
    character(len=worksize_string) :: workstring
    real :: lat, lon, hgt
    character(len=STATION_LABEL_LENGTH) :: label
    character(len=clen) :: name
    logical :: ifOK

    nullify(items)
    
    call get_items(nlPtr, 'station', items, nstations)
    
    if (nstations < 1) then
      call set_error('No stations found', 'stations_from_namelist')
      return
    end if

    allocate(station_list(nstations), stat=stat)
    if (stat /= 0) then
      call set_error('Allocate failed', 'stations_from_namelist')
      return
    end if
    
    iStat = 1
    do ii = 1, nstations
      workstring = fu_content(items(ii))
      read(workstring, fmt=*, iostat=stat) label, lat, lon, hgt, name
      if (stat /= 0) then
        call msg('Problem with record 2: ' + workstring + ', retrying without station name')
        read(workstring, fmt=*, iostat=stat) label, lat, lon, hgt
        if (stat /= 0) then
          call set_error('No way. Problem parsing station list', 'stations_from_namelist')
          return
        end if
        name = label
      end if
      station_list(iStat) = fu_initObservationStation(label, name, lon, lat, hgt)
      if(defined(station_list(iStat)))then
        iStat = iStat + 1
      endif
    end do
    nstations = iStat - 1

    deallocate(items)
    
  end subroutine stations_from_namelist_ptr

  !************************************************************************************
  
  subroutine read_observations(obs_start, obs_period, obs_files, nObsFiles_assim, nObsFiles_eval, &
                             & station_list, nstations, &
                             & in_situ_array, n_in_situ, &
                             & vert_obs_array, n_vert_obs, &
                             & eruption_obs_array, n_eruption, &
                             & dose_rate_addition_array, n_dose_rate, &
                             & transport_species, optical_species, n_transp_species, n_opt_species)
    !  
    ! Reads the given obs_files, each is an individual data file from the control file
    ! Note that here we do not distinguish between assimilation data and validation data
    ! Read them all and just pass the info on who is hwo above
    !
    implicit none
    type(silja_time), intent(in) :: obs_start
    type(silja_interval), intent(in) :: obs_period
    character(len=*), dimension(:), intent(in) :: obs_files
    type(observationStation), dimension(:), intent(in) :: station_list
    integer, intent(in) :: nstations
    integer, intent(in) :: nObsFiles_assim, nObsFiles_eval
    type(silam_species), dimension(:), intent(in) :: transport_species
    type(silam_species), dimension(:), intent(in) :: optical_species
    integer, intent(in) :: n_transp_species, n_opt_species
    ! out
    integer, dimension(2), intent(out) :: n_in_situ, n_vert_obs, n_dose_rate, n_eruption
    type(inSituObservation), dimension(:), pointer :: in_situ_array
    type(t_dose_rate_obs_addition), dimension(:), pointer :: dose_rate_addition_array
    type(t_vertical_observation), dimension(:), pointer :: vert_obs_array
    type(t_eruptionObservation), dimension(:), pointer :: eruption_obs_array

    ! Local variables
    integer :: status, uFile, nread, i,j, icounter=-1, ind_obs, n_vert_obs_def
    character(len=64) :: obs_type, obs_unit, var_name
    character(len=fnlen) :: file_name, obs_file, filename_templ, file_name_old
    type(inSituObservation), dimension(:), allocatable :: observations_tmp, dose_rate_obs_tmp
    type(t_dose_rate_obs_addition), dimension(:), allocatable :: dose_rate_addition_tmp
    type(t_vertical_observation), dimension(:), allocatable :: vert_obs_tmp
    type(t_eruptionObservation), dimension(:), allocatable :: eruption_obs_tmp
    type(silam_species), dimension(:), pointer :: p_opt_species, p_obs_species
    type(grads_template) :: obs_templ
    integer :: obs_index, iObsPurpose, iFileStart, iFileEnd
    logical :: ifGroundToo, force_instant, ifAssim
    type(silja_time) :: time_in_filename
    type(silja_interval) :: template_step
    !!FIXME some smarter allocation needed here
    integer, parameter :: max_n_obs_per_time_window = 24
    character(len=*), parameter :: sub_name = 'read_observations'


    call msg('Read observations, time limits: ' // fu_str(obs_start) // ' ' &
           & // fu_str(obs_start + obs_period))
    call start_count(sub_name)

    force_instant = (nint(fu_sec(obs_period)) == 0)

    n_in_situ = 0
    n_vert_obs = 0
    n_dose_rate = 0
    n_eruption = 0

    call msg("Before allocating observations memory usage (kB)",  fu_system_mem_usage())
    !! Trick: observations_tmp allocated for nstations+1, to allow for tsmatr, that needs no external stations
    allocate(observations_tmp(max_n_obs_per_time_window*(nstations+1)*n_transp_species), &
           & dose_rate_obs_tmp(nstations*n_transp_species), &
           & dose_rate_addition_tmp(nstations*n_transp_species), &
           & vert_obs_tmp(max_column_observations), &
           & eruption_obs_tmp(max_eruption_observations), &
           & stat=status)
    call msg("After allocating observations memory usage (kB)",  fu_system_mem_usage())

    if (fu_fails(status == 0, 'Allocate failed', sub_name)) return

    uFile = fu_next_free_unit()
    ifAssim = .true.
    iFileStart = 1
    iFileEnd = nObsFiles_assim
    !
    ! The observations files may be for two purposes: assimilation and evaluation
    ! The first set is always assimilation, the second - evaluation
    ! Process them one by one
    !
    do iObsPurpose = 1, 2
      !
      ! Scan over assimilation and evaluation lists of files
      !
      do i = iFileStart, iFileEnd
        !
        ! Process the obs_files one by one
        ! Each obs_file element has the following format:
        ! <obs_type> [<obs_unit>/<var_name>] <obs_file_name>
        !
        obs_file = obs_files(i)
        j = index(obs_file, ' ')
        obs_type = obs_file(1:j)
        file_name = trim(obs_file(j+1:))  ! together with obs_unit/var_name or just file_name
        !
        ! Processing depends on the type of observations
        !
        if (any(obs_type==(/'cnc    ','cncEU  ','cncWHO ','perkilo','permole', 'TSmatr '/))) then
          !
          ! In-situ concentrations
          ! file_name is <obs_unit> <obs_file_name>
          !
          time_in_filename = obs_start
          obs_unit = trim(file_name(1:index(file_name,' '))) 
          file_name = trim(file_name(index(file_name,' ')+1:))

          call decode_template_string(file_name, obs_templ)
          template_step = fu_template_timestep(obs_templ)
          if (error) return

          file_name_old = ""
          do while (time_in_filename <= obs_start + obs_period)
            call expand_template(obs_templ, time_in_filename, file_name)
            if (file_name /= file_name_old) then
              if (len_trim(file_name) > 0) then
                call timeseries_from_file(file_name, obs_type, obs_unit, observations_tmp(sum(n_in_situ)+1:), &
                                        & nread, force_instant, obs_start, obs_period, &
                                        & station_list, nstations, &
                                        & transport_species, n_transp_species)
                n_in_situ(iObsPurpose) = n_in_situ(iObsPurpose) + nread
                if (error) return
              end if
              file_name_old = file_name
            end if
            time_in_filename = time_in_filename + template_step
          end do
        
        elseif (obs_type == 'dose_rate_no_ground') then
          ifGroundToo = .false.
          call decode_template_string(file_name, obs_templ)
          if (error) return
          call expand_template(obs_templ, obs_start, file_name)
          if (len_trim(file_name) > 0) then
            open(uFile, file=file_name, status='old', action='read', iostat=status)
            if (fu_fails(status == 0, 'Failed to open: ' // trim(file_name), sub_name)) return
            call msg('Reading dose rate no-ground observations from ' // file_name)
            call timeseries_from_file_dose_rate(ifGroundToo, uFile, obs_unit, &
                                              & dose_rate_obs_tmp(sum(n_dose_rate)+1:), &
                                              & dose_rate_addition_tmp(sum(n_dose_rate)+1:), nread, &
                                              & station_list, nstations, obs_start, obs_period, &
                                              & transport_species, n_transp_species)
            close(uFile)
            n_dose_rate(iObsPurpose) = n_dose_rate(iObsPurpose) + nread
          end if
         
        elseif (obs_type == 'dose_rate_with_ground') then
          ifGroundToo = .true.
          call decode_template_string(file_name, obs_templ)
          if (error) return
          call expand_template(obs_templ, obs_start, file_name)
          if (len_trim(file_name) > 0) then
            open(uFile, file=file_name, status='old', action='read', iostat=status)
            if (fu_fails(status == 0, 'Failed to open: ' // trim(file_name), sub_name)) return
            call msg('Reading dose rate with-ground observations from ' // file_name)
            call timeseries_from_file_dose_rate(ifGroundToo, uFile, obs_unit, &
                                              & dose_rate_obs_tmp(sum(n_dose_rate)+1:), &
                                              & dose_rate_addition_tmp(sum(n_dose_rate)+1:), nread, &
                                              & station_list, nstations, obs_start, obs_period, &
                                              & transport_species, n_transp_species)
            close(uFile)
            n_dose_rate(iObsPurpose) = n_dose_rate(iObsPurpose) + nread
          end if

        elseif (obs_type == var_name_aod) then
          time_in_filename = obs_start
          call decode_template_string(file_name, obs_templ)
          if (error) return

          do while (time_in_filename <= obs_start + obs_period)
            call expand_template(obs_templ, time_in_filename, file_name)
            if (len_trim(file_name) > 0) then
              call msg('Reading AOT from '  // trim(file_name))
              ! obs_species = null() - not used                                                                                          
              nullify(p_obs_species)
              call set_vert_obs_from_nc(file_name, var_name_aod, p_obs_species, &
                                      & transport_species, optical_species, &
                                      & obs_start, obs_start+obs_period, &
                                      & vert_obs_tmp(sum(n_vert_obs)+1:), nread)
              !if ((rules%observation_stdev > 0) .and. (nread > 0)) then                                                                 
              !  call msg('nread', nread)                                                                                                
              !  do obs_index = n_vert_obs+1, nread                                                                                      
              !    vert_obs_tmp(obs_index)%variance = (rules%observation_stdev)**2                                                       
              !  end do                                                                                                                  
              !end if                                                                                                                    
              if (error) return
              n_vert_obs(iObsPurpose) = n_vert_obs(iObsPurpose) + nread
            end if
            time_in_filename = time_in_filename + fu_set_interval_h(1)
          end do

        else if (obs_type == 'eruption') then
          ! file_name is together with obs_unit 
          obs_unit = trim(file_name(1:index(file_name,' ')))
          file_name = trim(file_name(index(file_name,' ')+1:))
          call decode_template_string(file_name, obs_templ)
          if (error) return
          call expand_template(obs_templ, obs_start, file_name)
          if (len_trim(file_name) > 0) then
            open(uFile, file=file_name, status='old', action='read', iostat=status)
            if (fu_fails(status == 0, 'Failed to open: ' // trim(file_name), sub_name)) return
            call msg('Reading eruption observations from ' // file_name)
            call eruption_timeseries_from_file(uFile, obs_unit, eruption_obs_tmp(sum(n_eruption)+1:), &
                                             & nread)
            close(uFile)
            n_eruption(iObsPurpose) = n_eruption(iObsPurpose) + nread
            if (error) return
          end if

        else if (obs_type == 'vertical') then
          ! file_name is together with var_name 
          time_in_filename = obs_start

          var_name = trim(file_name(1:index(file_name,' ')))
          file_name = trim(file_name(index(file_name,' ')+1:))
          call decode_template_string(file_name, obs_templ)
          if (error) return
          do while (time_in_filename <= obs_start + obs_period)
            call expand_template(obs_templ, time_in_filename, file_name)
            if (len_trim(file_name) > 0) then
              ! var_name is also the observation cocktail:                                                         
              call msg('Reading vertical observation from '  // trim(file_name))
              call get_observed_species(var_name, p_obs_species)
              if (error) return
              call set_vert_obs_from_nc(file_name, var_name, p_obs_species, &
                   & transport_species, optical_species, &
                   & obs_start, obs_start+obs_period, &
                   & vert_obs_tmp(sum(n_vert_obs)+1:), nread)
              if (error) return
              deallocate(p_obs_species)
              n_vert_obs(iObsPurpose) = n_vert_obs(iObsPurpose) + nread
            end if
            time_in_filename = time_in_filename + fu_set_interval_h(1)
          end do
        else
          call set_error('Unknown observation type:' // trim(obs_type), sub_name)
          return
        end if

      end do ! types of observations
      !
      ! switch the purpose
      !
      ifAssim = .false. 
      iFileStart = nObsFiles_assim + 1
      iFileEnd = nObsFiles_assim + nObsFiles_eval

    end do ! purpose of observations: assimilation or evaluation

    if ( sum(n_in_situ + n_vert_obs + n_dose_rate + n_eruption) < 1) then
       call msg_warning( 'No observations read succesfully', sub_name) 
       return
    endif
    !
    ! Allocate memory
    !
    if (sum(n_in_situ) > 0) then
      allocate(in_situ_array(sum(n_in_situ)), stat=status)
      if (fu_fails(status == 0, 'Allocate failed', sub_name)) return
      in_situ_array = observations_tmp(1:sum(n_in_situ))
    end if
    
    if (sum(n_dose_rate) > 0) then
      allocate(dose_rate_addition_array(sum(n_dose_rate)), stat=status)
      if (fu_fails(status == 0, 'Allocate failed', sub_name)) return
      dose_rate_addition_array = dose_rate_addition_tmp(1:sum(n_dose_rate))
    end if

    if (sum(n_eruption) > 0) then
      allocate(eruption_obs_array(sum(n_eruption)), stat=status)
      if (fu_fails(status == 0, 'Allocate failed', sub_name)) return
      eruption_obs_array = eruption_obs_tmp(1:sum(n_eruption))
   end if
   ! FIXME: the following does not work, and no defined vertical observations seem to be
   ! ever found for assimilation. Also, n_vert_obs is a vector.
    if (sum(n_vert_obs) > 0) then
      ! some observations could be missing because they are outside, etc. Ignore them.
      i = 0
      do iObsPurpose = 1,2
        n_vert_obs_def = count(defined(vert_obs_tmp(i+1:i+n_vert_obs(iObsPurpose))))
        allocate(vert_obs_array(n_vert_obs_def), stat=status)
        if (fu_fails(status == 0, 'Allocate failed', sub_name)) return
        n_vert_obs_def = 0
        do ind_obs = i + 1, i + n_vert_obs(iObsPurpose)
          if (defined(vert_obs_tmp(ind_obs))) then
            n_vert_obs_def = n_vert_obs_def + 1
          else
            cycle
          end if
          vert_obs_array(n_vert_obs_def) = vert_obs_tmp(ind_obs)
        end do
        call msg('Vertical observation loaded/defined:', n_vert_obs, n_vert_obs_def)
        i = n_vert_obs(iObsPurpose) ! evaluation obs start from this index
        n_vert_obs = n_vert_obs_def ! now reset the number of assimilation
      end do ! iObsPurpose
    end if

    deallocate(observations_tmp, vert_obs_tmp, dose_rate_addition_tmp, &
             & dose_rate_obs_tmp, eruption_obs_tmp)

    call stop_count(sub_name)
    !call report_time(icounter, chCounterNm=sub_name)

  end subroutine read_observations

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

    subroutine timeseries_from_TSmatr(ncfilename, cockt_name, newObservation, force_instant,&
          & obs_start, obs_period, transport_species, n_transp_species)
      !
      ! Get timeseries from TSmatrix netcdf file
      ! The file must have two variables: 
      !  (vmr and vmrErr) or  (cnc and cncErr)
      !  
      use netcdf
      implicit none
      
      ! Imported parameters
      character(len=*), intent(in) :: ncfilename, cockt_name !! Cocktail name not stored in TSmatrix
      type(inSituObservation), intent(out) :: newObservation
      logical, intent(in) :: force_instant !Disregard duration and force it to zero
      type(silja_time), intent(in) :: obs_start
      type(silja_interval), intent(in) :: obs_period
      type(silam_species), dimension(:), intent(in) :: transport_species
      integer, intent(in) :: n_transp_species

      type(Tcocktail_descr) :: cockt_descr
      type(silam_species), dimension(:), pointer :: cockt_species

      integer :: ncid, var_id, dim_id !! ids
      integer :: iStat, nSt, nStr, nVar, nT, nSpecies, nStout, nTout !! sizes
      integer :: iT, iTout, iSt, iStout, isp_obs, ii, iTstart, iTend !! Counters

      !!Stuff needed fo NetCDF parsing
      character (len=fnlen) :: strTmp1, strTmp2
      character (len=clen) :: obs_amt_unit
      TYPE(silja_time) :: origin
      TYPE(silja_interval) :: deltat
      real :: fMissVal, scale_to_silam, fScale
      real, dimension(:), allocatable :: lat, lon, alt, timesd
      logical, dimension(:), allocatable :: timesselect, stationselect
      logical :: ifmmr, ifSpecies
      !      real (kind=8), dimension(:), allocatable :: timesd
      real, dimension(:,:,:), allocatable :: valvar
      character, dimension(:,:), allocatable :: station_name, station_code, varnames, units
      TYPE(silja_time), dimension(:), allocatable :: times
      TYPE(observationStation), dimension(:), allocatable :: station
      type(chemical_adaptor) :: adaptor
      type(silam_material), pointer :: material

      character(len=*), parameter :: sub_name = 'timeseries_from_TSmatr'

      !! How many species?
      call msg("Reading cocktail '"//trim(cockt_name)//"' from '"//trim(ncfilename)//"'")

      call set_cocktail_description(cockt_name, cockt_descr, ifSpecies)
      call get_inventory(cockt_descr, cockt_species, nSpecies)


      
      !! Check tsmatrix dimensions
      iStat = nf90_open(ncfilename, 0, ncid)
      if (fu_fails(iStat == NF90_NOERR, 'Failed to open nc file: ' // trim(ncfilename), sub_name)) return

      iStat = nf90_inq_dimid(ncid, 'station', dim_id)
      if (fu_fails(iStat == NF90_NOERR, 'Failed to inquire station dimID', sub_name)) return
      iStat = NF90_INQUIRE_DIMENSION(ncid, dim_id, len=nSt)
      if (fu_fails(iStat == NF90_NOERR, 'Failed to inquire station dim', sub_name)) return

      iStat = nf90_inq_dimid(ncid, 'time', dim_id)
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get time dimID', sub_name)) return
      iStat = NF90_INQUIRE_DIMENSION(ncid, dim_id, len=nT)
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get time dim size', sub_name)) return

      iStat = nf90_inq_dimid(ncid, 'name_strlen', dim_id)
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get name_strlen dimID', sub_name)) return
      iStat = NF90_INQUIRE_DIMENSION(ncid, dim_id, len=nStr)
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get name_strlen dim size', sub_name)) return

      iStat = nf90_inq_dimid(ncid, 'variable', dim_id)
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get variable dimID', sub_name)) return
      iStat = NF90_INQUIRE_DIMENSION(ncid, dim_id, len=nVar)
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get variable dim size', sub_name)) return
    
      if (fu_fails(nVar == 2, "file must contain 2 vars", sub_name) )return

      !! Temporary structures fot TSmatrix, essentially just arrays to read NetCDF
      allocate(times(nT),timesd(nT),timesselect(nT), &
         & stationselect(nSt), lat(nSt), lon(nSt), alt(nSt), station_code(nStr,nSt), &
         & station_name(nStr,nSt), units(nStr,nVar), station(nSt), &
         & varnames(nStr,nVar), valvar(nSt,nT,nVar), stat=iStat)
      if (fu_fails(iStat == NF90_NOERR, 'Failed to allocate temporaries', sub_name)) return


      !! Parse and select times
      iStat = nf90_inq_varid(ncid, 'time', var_id)
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get time var_id', sub_name)) return
      iStat = nf90_get_att(ncid, var_id, 'units', strTmp1)
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get time units', sub_name)) return
      iStat = nf90_get_att(ncid, var_id, 'calendar', strTmp2)
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get time units', sub_name)) return
      iStat = NF90_get_var(ncid, var_id, timesd, start=(/1/), count=(/nT/))
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get time data', sub_name)) return

      ! Time list
      call parse_time_units_and_origin(strTmp1, strTmp2, origin, deltat)
      do iT  =  1, nT
          times(iT) = origin + deltat * timesd(iT)
      enddo

      ! Times must be equally-spaced unless the values are instant
      if (force_instant) then 
          deltat = zero_interval
      else
        deltat = (times(nT) - times(0)) * (1. / (nT - 1))
        do iT  =  2, nT !! Check for equal spacing
            if ( times(iT) - times(iT-1) == deltat) cycle
            call report(times(iT-1))
            call report(times(iT))
            call report(deltat)
            call set_error("Unequally spaced times in "//trim(ncfilename), sub_name)
        enddo
      endif
      ! Select times within th eobservation range
      do iT=1,nT
        if (times(iT) - deltat >= obs_start) exit
        timesselect(iT) = .FALSE.
      enddo
      iTstart = iT
      iTend = 0 !! (iTstart > iTend) on no valid timestamps
      do iT=iT,nT !! All iT have to be iterated, observation must be fully covered!
        timesselect(iT) = (times(iT) <= obs_start + obs_period)
        if (timesselect(iT)) iTend = iT
      enddo

      !! Get station features (area_type, source_type ignored here)
      iStat = nf90_inq_varid(ncid, 'lon', var_id)
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get lon var_id', sub_name)) return
      iStat = NF90_get_var(ncid, var_id, lon, start=(/1/), count=(/nSt/))
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get lon data', sub_name)) return

      iStat = nf90_inq_varid(ncid, 'lat', var_id)
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get lat var_id', sub_name)) return
      iStat = NF90_get_var(ncid, var_id, lat, start=(/1/), count=(/nSt/))
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get lat data', sub_name)) return

      iStat = nf90_inq_varid(ncid, 'alt', var_id)
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get alt var_id', sub_name)) return
      iStat = NF90_get_var(ncid, var_id, alt, start=(/1/), count=(/nSt/))
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get alt data', sub_name)) return

      iStat = nf90_inq_varid(ncid, 'station_code', var_id)
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get station_code var_id', sub_name)) return
      iStat = NF90_get_var(ncid, var_id, station_code, start=(/1,1/), count=(/nStr,nSt/))
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get station_code data', sub_name)) return

      iStat = nf90_inq_varid(ncid, 'station_name', var_id)
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get station_name var_id', sub_name)) return
      iStat = NF90_get_var(ncid, var_id, station_name, start=(/1,1/), count=(/nStr,nSt/))
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get station_name data', sub_name)) return

      !! Make stations and select them
      do iSt = 1, nSt
          strTmp1 = nc_char_arr_to_string(station_code(:,iSt))
          strTmp2 = nc_char_arr_to_string(station_name(:,iSt))
          station(iSt) = fu_initObservationStation(strTmp1, strTmp2, lon(iSt), lat(iSt), alt(iSt)) 
          stationselect(iSt) = defined(station(iSt))
      enddo

      !! Figure out variable and unit
      iStat = nf90_inq_varid(ncid, 'variable', var_id)
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get variable var_id', sub_name)) return
      iStat = NF90_get_var(ncid, var_id, varnames, start=(/1,1/), count=(/nStr,2/))
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get variable data', sub_name)) return
      strTmp1 = nc_char_arr_to_string(varnames(:,1)) !! Variable
      strTmp2 = nc_char_arr_to_string(varnames(:,2)) !! its std (err)
      ii = len_trim(strTmp1)
      if (strTmp2(1:ii) /= strTmp1 .or. strTmp2(ii+1:) /= 'err' ) then
        call set_error("Inconsistent 'variable' var: "//trim(strTmp1)//" and "//trim(strTmp2), sub_name)
        return
      endif

      ! Check unit
      iStat = nf90_inq_varid(ncid, 'unit', var_id) !! reuse varnames array
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get variable var_id', sub_name)) return
      iStat = NF90_get_var(ncid, var_id, varnames, start=(/1,1/), count=(/nStr,2/))
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get variable data', sub_name)) return
      if (.not. all(varnames(1,:) == varnames(1,:) )) then
        call set_error("Variable and error units mismatch", sub_name)
        return
      endif
      strTmp2 = nc_char_arr_to_string(varnames(:,1))  !! Unit

      call  parse_obs_unit(strTmp1, strTmp2, obs_amt_unit, scale_to_silam, ifmmr)
      call msg("Variabe from file: "//trim(strTmp1)//", unit: "//trim(strTmp2) )
      if (ifmmr) then
         call msg("Conversion factor to "//trim(obs_amt_unit)//"/kg", scale_to_silam)
      else
        call msg("Conversion factor to "//trim(obs_amt_unit)//"/m3", scale_to_silam)
      endif

      !! read the data
      iStat = nf90_inq_varid(ncid, 'val', var_id)
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get sval var_id', sub_name)) return
      !! Read only needed times

      if (iTend >= iTstart) then !! Some data to read
        iStat = NF90_get_var(ncid, var_id, valvar(:,iTstart,1), start=(/1,iTstart,1/), count=(/nSt,iTend-iTstart+1,1/))
        if (fu_fails(iStat == NF90_NOERR, 'Failed to get val data1', sub_name)) return
        iStat = NF90_get_var(ncid, var_id, valvar(:,iTstart,2), start=(/1,iTstart,2/), count=(/nSt,iTend-iTstart+1,1/))
        if (fu_fails(iStat == NF90_NOERR, 'Failed to get val data2', sub_name)) return
      endif
      iStat = nf90_get_att(ncid, var_id, '_FillValue', fMissVal)
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get _FillValue', sub_name)) return
      iStat =  nf90_close(ncid)
      if (fu_fails(iStat == NF90_NOERR, 'Failed to get _FillValue', sub_name)) return


      !! Now finally init the newobservation structure
      nTout = count(timesselect(1:nT))
      nStout = count(stationselect(1:nSt))


      !! Finally allocate the observations
      allocate(newobservation%endTimes(nTout), newObservation%durations(nTout), &
               & newObservation%station(nStout), &
               & newObservation%cell_volume(nStout), &
               & newObservation%obsData(nStout,nTout), newObservation%variance(nStout,nTout), &
               & newObservation%modelData(nStout,nTout), newObservation%inject_scaling(nStout,nTout), &
               & newObservation%ind_obs2transp(nSpecies), &
               & newObservation%scale_transp2obs(nSpecies), &
               & stat=iStat)

      if(fu_fails(iStat == 0, 'Allocate failed', sub_name))return

      
      newObservation%tag =  trim(cockt_name)//':'//trim(strTmp1)//':'//trim(strTmp2) 
      newObservation%ifmmr = ifmmr

      newObservation%nStations = nStout
      newObservation%nTimes = nTout

      iStout = 0
      newObservation%inThisSubdomain = .FALSE.
      do iSt = 1, nSt 
        if (.not. stationselect(iSt)) cycle
        iStout = iStout + 1
        newObservation%station(iStout) = station(iSt)
        newObservation%inThisSubdomain = newObservation%inThisSubdomain .OR. station(iSt)%inThisSubdomain 
        newObservation%cell_volume(iStOut) = newObservation%station(iStOut)%cell_area * &
                                            & fu_layer_thickness_m(fu_level(dispersion_vertical, 1))
      enddo 
      
      newObservation%dataLength = 0
      iTout = 0
      do iT = iTstart, iTend
        if (.not. timesselect(iT)) cycle
        iTout = iTout + 1
        newObservation%endTimes(iTout) = times(iT)
        newObservation%durations(iTout) = deltat
        iStout = 0
        do iSt = 1, nSt 
          if (.not. stationselect(iSt)) cycle
          iStout = iStout + 1
          if ( valvar(iSt,iT,1) == fMissVal) then
              newObservation%obsData(iStout,iTout)  = real_missing
              newObservation%variance(iStout,iTout) = F_NAN !!! Should never be used
          else
              newObservation%obsData(iStout,iTout)  = valvar(iSt,iT,1) 
              newObservation%variance(iStout,iTout) = &  !! STD stored there 
                  &   valvar(iSt,iT,2) * valvar(iSt,iT,2) 
              newObservation%dataLength = newObservation%dataLength  + 1 ! Only valid points count
          endif
        enddo
      enddo

    newObservation%num_obs_species = nSpecies

    !! Map observation to massmap
    call create_adaptor(cockt_species, transport_species, adaptor)
    newObservation%ind_obs2transp(1:nSpecies) = adaptor%isp(1:nSpecies)

    do isp_obs = 1, nSpecies
      material => fu_material(transport_species(adaptor%isp(isp_obs)))
      if (error) return
      fScale = fu_conversion_factor(fu_basic_mass_unit(material), obs_amt_unit, material) / scale_to_silam
      ! From basic unit to the unit of observation - preparation to the obs operator.
      newObservation%scale_transp2obs(isp_obs) = fScale
      if (error) return
      if (fu_fails( fScale > 0.0, 'Bad conversion factor to obs unit', sub_name)) return
      call msg('Factor transp2obs' //trim(fu_str(cockt_species(isp_obs))), fScale )
    end do

    !! Clean up the mess
    deallocate(times,timesd,timesselect, &
         & stationselect, lat, lon, alt, station_code, &
         & station_name, units, station, &
         & varnames, valvar)

    end subroutine timeseries_from_TSmatr

    !================================================================================

    

    subroutine timeseries_from_file(file_name, obs_param, obs_amt_unit, observations, nread,&
         & force_instant,  obs_start, obs_period, station_list, nstations, &
         & transport_species, n_transp_species)
      !
      ! Get in_situ timeseries from a file
      ! Dispatcher for text and netcdf files
      !
      implicit none
      
      ! Imported parameters
      character(len=*), intent(in) :: file_name, obs_param 
         !! can be cnc (per m3), cncEU (per 1.204kg of air), cncWHO(per 1.184 kg of air ), 
         !! vmr (per mole of air), mmr (per kilo of air), or TSmatr (quantity comes from NetCDF)
      character(len=*), intent(in) :: obs_amt_unit !! kg or mole or whatever
                                 !! cocktail name for TSmatr
      type(inSituObservation), dimension(:), intent(inout) :: observations
      logical, intent(in) :: force_instant !Disregard duration and force it to zero
      integer, intent(out) :: nread
      type(observationStation), dimension(:), intent(in) :: station_list
      integer, intent(in) :: nstations
      type(silja_time), intent(in) :: obs_start
      type(silja_interval), intent(in) :: obs_period
      type(silam_species), dimension(:), intent(in) :: transport_species
      integer, intent(in) :: n_transp_species

      integer :: uFile, iStat

      character(len=*), parameter :: sub_name = 'timeseries_from_file'

      if (obs_param == "TSmatr") then
        !!  TSmatrix -- multistation observation with regular times
        if (1 > size(observations)) then
          call set_error('Too many observations', sub_name)
          return
        end if
        !!!!cockt_name = obs_amt_unit !! The second word in the line is used for cocktail
        call timeseries_from_TSmatr(file_name, obs_amt_unit, observations(1), &
           & force_instant, obs_start, obs_period, &
           & transport_species, n_transp_species)
        nread = 1
      else
        !!  MMAS file with a bunch of single-station timeseries
         if (fu_fails(nstations > 0, 'Surface observation given but no stations', sub_name)) return
           
        uFile = fu_next_free_unit()
        open(uFile, file=file_name, status='old', action='read', iostat=iStat)
        if (fu_fails(iStat == 0, 'Failed to open: ' // trim(file_name), sub_name)) return
        call msg('Reading point observations from ' // file_name)

        call  timeseries_from_text_file(uFile, obs_param, obs_amt_unit, observations, nread,&
         & force_instant,  obs_start, obs_period, station_list, nstations, &
         & transport_species, n_transp_species)
         close(uFile)
      endif
    end subroutine timeseries_from_file
    !================================================================================

    

    subroutine timeseries_from_text_file(uFile, obs_param, obs_amt_unit, observations, nread,&
         & force_instant,  obs_start, obs_period, station_list, nstations, &
         & transport_species, n_transp_species)
      !
      ! Get timeseries from text file. The MMAS format
      ! site_id, cocktail/species_name, year, month, day, hour, duration_hours, value, stdev
      !
      implicit none
      
      ! Imported parameters
      integer, intent(in) :: uFile
      character(len=*), intent(in) :: obs_param 
         !! can be cnc (per m3), cncEU (per 1.204kg of air), cncWHO(per 1.184 kg of air ), 
         !! vmr (per mole of air), mmr (per kilo of air), or TSmatr (quantity comes from NetCDF)
      character(len=*), intent(in) :: obs_amt_unit !! kg or mole or whatever
                                 !! cocktail name for TSmatr
      type(inSituObservation), dimension(:), intent(inout) :: observations
      logical, intent(in) :: force_instant !Disregard duration and force it to zero
      integer, intent(out) :: nread
      type(silja_time), intent(in) :: obs_start
      type(silja_interval), intent(in) :: obs_period
      type(observationStation), dimension(:), intent(in) :: station_list
      integer, intent(in) :: nstations
      type(silam_species), dimension(:), intent(in) :: transport_species
      integer, intent(in) :: n_transp_species

      ! Local variables
      logical :: eof, same_series, found, same_cocktail, store_value, flush_values
      integer :: iostat

      character(len=station_id_length) :: id, prev_id
      character(substNmLen) :: cockt_name, prev_cockt_name
      character(len=255) :: line
      integer :: year, month, day, hour, itime, quantity, n_cockt_species, istation
      real :: duration_hours, value, modeval, wavelength, stdev, variance
      type(silja_interval) :: duration
      type(silja_time) :: time, prev_time
      real, dimension(:), pointer :: values_tmp, variances_tmp
      type(silja_time), dimension(:), allocatable :: times_tmp
      type(silja_interval), dimension(:), allocatable :: durations_tmp
      type(silam_species), dimension(:), pointer :: cockt_species
      type(Tcocktail_descr) :: cockt_descr

      integer, dimension(max_species) :: ind_obs2transp
      real, dimension(max_species) :: scale_transp2obs
      integer :: status, nx, ny, i, isp_transp, isp_obs
      type(silam_material), pointer :: material
      type(chemical_adaptor) :: adaptor
      logical :: ifSpecies
      real :: scale_to_silam
      logical :: ifmmr

      character(len=*), parameter :: sub_name = 'timeseries_from_text_file'

      prev_time = time_missing
      prev_cockt_name = ''
      prev_id = ''
      nread = 0
      itime = 0
      values_tmp => fu_work_array()
      variances_tmp => fu_work_array()
      allocate(times_tmp(worksize), durations_tmp(worksize), stat=iostat)
      if (fu_fails(iostat == 0,'Allocate failed', sub_name))return

      !! Figure out what we observe
      call parse_obs_param(obs_param, scale_to_silam, ifmmr)

      !read(uFile,fmt=*) line

      do
        call next_line_from_input_file(uFile, line, eof)
        if (eof) iostat = -1
        !read(unit, fmt='(A)', iostat=iostat) line
        if (iostat > 0) then
          call set_error('Failed to read record', sub_name)
          return
        end if
        store_value = .false.
        flush_values = .false.
        if (iostat == 0) then
          read(line, fmt=*, iostat=iostat) &
               & id, cockt_name, year, month, day, hour, duration_hours, value, stdev
          variance = stdev**2
          if (iostat /= 0) then
            call set_error('Failed to parse record: ' // trim(line), sub_name)
            return
          end if

          if (force_instant) then
            duration = zero_interval
          else
            duration = fu_set_interval_sec(duration_hours*3600.0)
          endif
          time = fu_set_time_utc(year, month, day, hour, 0, 0.0)

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
              call msg_warning('Cannot read: timeseries not sorted', sub_name)
              call msg(id + ', prev_time > time:' + fu_time_fname_string_utc(prev_time) + &
                                            & '>' + fu_time_fname_string_utc(time))
              call set_error('Cannot read: timeseries not sorted', sub_name)
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

            do isp_obs = 1, n_cockt_species
              material => fu_material(transport_species(ind_obs2transp(isp_obs)))
              if (error) return
              ! From basic unit to the unit of observation - preparation to the obs operator.
              scale_transp2obs(isp_obs) = fu_conversion_factor(fu_basic_mass_unit(material), obs_amt_unit, material) / scale_to_silam
              if (error) return
              if (fu_fails(scale_transp2obs(isp_obs) > 0.0, 'Bad conversion factor to obs unit', sub_name)) return
!!                call msg('Observed species, transport species', isp_obs, ind_obs2transp(isp_obs))
              call msg('Factor transp2obs' //trim(fu_str(cockt_species(isp_obs))), scale_transp2obs(isp_obs))
            end do
          end if
        else ! eof
          flush_values = itime > 0
        end if ! read ok
        
        flush_values = flush_values .and. itime > 0

        if (flush_values) then
          
          call searchStationWithID(prev_id, station_list, nstations, istation)
          if (istation < 1) then
            call msg_warning('Station ' // trim(prev_id) // ' not found')
          else
            nread = nread + 1
            if (nread > size(observations)) then
              call set_error('Too many observations', sub_name)
              return
            end if
            observations(nread) = fu_initInSituObservation1(times_tmp, &
                                                         & durations_tmp, &
                                                         & values_tmp, &
                                                         & variances_tmp, &
                                                         & itime, &
                                                         & variable_variance, &
                                                         & ifmmr, &
                                                         & ind_obs2transp, &
                                                         & scale_transp2obs, &
                                                         & n_cockt_species, &
                                                         & station_list(istation), &
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
    end subroutine timeseries_from_text_file

    !========================================================================================

    subroutine eruption_timeseries_from_file(uFile, obs_unit, observations, nread)
      implicit none
      integer, intent(in) :: uFile
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
      character(len=*), parameter :: sub_name = 'eruption_timeseries_from_file'

      value_tmp => fu_work_array()
      variance_tmp => fu_work_array()
      lat_tmp => fu_work_array()
      lon_tmp => fu_work_array()

      allocate(time_tmp(tmp_arr_size), duration_tmp(tmp_arr_size), stat=iostat)
      if (iostat /= 0) then
        call set_error('Allocate failed', sub_name)
        return
      end if

      nread = 0
     
      do
        call next_line_from_input_file(uFile, line, eof)
        if (eof) iostat = -1
        !read(unit, fmt='(A)', iostat=iostat) line                                                                                                                                                                                                                                                                          
        if (iostat > 0) then
          call set_error('Failed to read record', sub_name)
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
              call set_error('Too many observations', sub_name)
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

    subroutine timeseries_from_file_dose_rate(ifGroundToo, uFile, obs_unit, observations_dose_rate, &
          & dose_rate_addition, nread, station_list, nstations, obs_start, obs_period, &
         & transport_species, n_transp_species)
      !edit: no cocktail here compared to normal timeseries from file!
      !dose rate measurement doesn't know what is emitting
      implicit none
      logical, intent(in) :: ifGroundToo
      integer, intent(in) :: uFile
      character(len=*), intent(in) :: obs_unit

      type(inSituObservation), dimension(:), intent(inout) :: observations_dose_rate
      
      type(t_dose_rate_obs_addition), dimension(:), intent(out) :: dose_rate_addition
      integer, intent(out) :: nread
      type(observationStation), dimension(:), intent(in) :: station_list
      integer, intent(in) :: nstations
      type(silja_time), intent(in) :: obs_start
      type(silja_interval), intent(in) :: obs_period
      type(silam_species), dimension(:), intent(in) :: transport_species
      integer, intent(in) :: n_transp_species
      
      logical :: eof, same_series, found, store_value, flush_values
      integer :: iostat
      character(len=station_id_length) :: id, prev_id
      character(len=255) :: line
      integer :: year, month, day, hour, itime, quantity, istation
      real :: duration_hours, value, modeval, wavelength, stdev, variance
      type(silja_interval) :: duration
      type(silja_time) :: time, prev_time
      real, dimension(:), pointer :: values_tmp, variances_tmp
      type(silja_time), dimension(:), allocatable :: times_tmp
      type(silja_interval), dimension(:), allocatable :: durations_tmp
      integer, parameter :: tmp_arr_size = 20
      character(len=*), parameter :: sub_name = 'timeseries_from_file_dose_rate'
      
      prev_time = time_missing
      prev_id = ''
      nread = 0
      itime = 0
      values_tmp => fu_work_array()
      variances_tmp => fu_work_array()
      allocate(times_tmp(tmp_arr_size), durations_tmp(tmp_arr_size), stat=iostat)
      if (iostat /= 0) then
        call set_error('Allocate failed', sub_name)
        return
      end if

      !read(uFile,fmt=*) line

      do
        call next_line_from_input_file(uFile, line, eof)
        if (eof) iostat = -1
        !read(unit, fmt='(A)', iostat=iostat) line
        if (iostat > 0) then
          call set_error('Failed to read record', sub_name)
          return
        end if
        store_value = .false.
        flush_values = .false.
        if (iostat == 0) then
          read(line, fmt=*, iostat=iostat) &
               & id, year, month, day, hour, duration_hours, value, stdev
          variance = stdev**2
          if (iostat /= 0) then
            call set_error('Failed to parse record: ' // trim(line), sub_name)
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
              call msg_warning('Cannot read: timeseries not sorted', sub_name)
              call msg(id + ', prev_time > time:' + fu_time_fname_string_utc(prev_time) + &
                                            & '>' + fu_time_fname_string_utc(time))
              call set_error('Cannot read: timeseries not sorted', sub_name)
              return
            end if
          end if
        else ! eof
          flush_values = itime > 0
        end if ! read ok
        
        flush_values = flush_values .and. itime > 0

        if (flush_values) then
          !call msg(id)
          call searchStationWithID(prev_id, station_list, nstations, istation)
          if (istation < 1) then
            call msg_warning('Station ' // trim(prev_id) // ' not found')
          else
            nread = nread + 1
            if (nread > size(dose_rate_addition)) then
              call set_error('Too many observations', sub_name)
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
                                                         & station_list(istation), &
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
    character(len=*), intent(out) :: filename  !!Empty string for missing file

    type(silam_sp), dimension(:), pointer :: filenames
    integer :: i

    filename = ""
    nullify(filenames)

    if(fu_fails(defined(templ),'Undefined template','expand_template'))return
    call fnm_from_single_template(templ, now, filenames, &
                                & ifStrict = .true., &
                                & ifadd = .false., &
                                & ifWait = .false., &
                                & max_hole_length = one_hour, & !!Actually whatever to suppress set_error on missing file 
                                & ifAllowZeroFcLen = .true.)
    if (.not. associated(filenames)) then
      !call set_error('No files found for time:' // fu_str(now), &
      !             & 'expand_template')
      call msg('No files found for time:' // fu_str(now)) 
      return
    end if
    if (size(filenames) /= 1) then
      do i=1,size(filenames)
        call msg("fname "//trim(fu_str(i))//": '"//trim(filenames(i)%sp)//"'")
      enddo
      call set_error('Only one fname allowed for time:' // fu_str(now), &
                   & 'expand_template')
    end if
    if (error) return

    filename = fu_process_filepath(filenames(1)%sp, must_exist=.true.)
    deallocate(filenames(1)%sp)
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

  subroutine observeAll(pointers, cloud, metBuf, disp_buf_ptr, chemRules, dynRules, timestep, now, eruption_height)
    implicit none
    type(observationPointers), intent(inout) :: pointers
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep
    type(silam_pollution_cloud), pointer :: cloud
    type(TChem_rules), intent(in) :: chemRules
    type(TDynamics_rules), intent(in) :: dynRules
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
    integer :: i, iQ, ind_tempr, ind_rh, ind_surf_press, ind_hgt, ind_press, iSpOnes
    
    type(TchemicalRunSetup), pointer :: chem_run_setup
    
    ! To observe mean-grid-cell concentrations, use nearest_point
    ! To observe from trapezoid with mass centre as a proxy for within-cell interpolation, use toMassCentreLinear
    ! 
    character(len=*), parameter :: sub_name = 'observeAll'
    
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

    if (sum(pointers%nVerticalObsID) > 0) then
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

    iSpOnes = int_missing
    if (dynRules%cloud_metric == cloud_metric_ones_flag) then
      call msg(sub_name // " will use ones tracer as cell size")
      iSpOnes = select_single_species(mapConc%species, mapConc%nSpecies, &
                        & 'ones', in_gas_phase, real_missing)
      if (error .or. (iSpOnes < 1)) then
        call set_error("Couldn't find 'ones' species", sub_name)
        return
      endif
    endif

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
    !$OMP &  p_press, p_hgt, p_surf_press, p_met_disp_interp_horiz, dynRules, iSpOnes, &
    !$OMP if_met_disp_interp_horiz, p_met_disp_interp_vert, chem_run_setup, &
    !$OMP &  if_met_disp_interp_vert, metbuf, disp_buf_ptr, chemrules, cloud, eruption_height)
!    call msg("observe_in_situ called nthreads", omp_get_num_threads(), omp_get_thread_num())
    !$OMP DO
    do i = 1, pointers%nInSituObsID(1) + pointers%nInSituObsID(2)
      call observe_in_situ(pointers%observationsInSitu(i), mapConc, map_px, map_py, disp_buf_ptr, &
                         & dynRules%cloud_metric, iSpOnes,  now, timestep)
    end do
    !$OMP END DO
    
    !$OMP DO
    do i = 1, pointers%nDoseRateObsID(1) + pointers%nDoseRateObsID(2)
       call observe_dose_rate(pointers%DoseRateAddition(i), pointers%observationsDoseRate(i), &
                            & mapConc, mapWetdep, mapDrydep, now, timestep)
    end do
    !$OMP END DO

    !$OMP DO
    do i = 1, pointers%nVerticalObsID(1) + pointers%nVerticalObsID(2)
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
                                                                                                                                                                                                                                                                                                               
    if (present(eruption_height) )then
      if  (.not. eruption_height == real_missing) then
        !$OMP DO
        do i = 1, pointers%nEruptionObsID(1) + pointers%nEruptionObsID(2)
          call observe_eruption(pointers%observationsEruption(i), eruption_height)
        end do
        !$omp End Do
      endif
    end if
    !$omp End Parallel
    
  end subroutine observeAll


  !**************************************************************************************
  
  subroutine injectAll(pointers, cloud, metBuf, dispBuf, chemRules, dynRules, &
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
    type(TField_buffer), intent(in) :: dispBuf
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
    integer :: ithread, nthreads, iSpOnes
    
    ! To inject into centre of the cell, use nearest_point
    ! To inject into actual site location, use toMassCentreLinear
    ! 
    character(len=*), parameter :: sub_name = 'injectAll'

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
    if (sum(pointers%nVerticalObsID) > 0) then
      ind_rh = fu_index(metbuf, relative_humidity_flag)
      if (fu_fails(ind_rh /= int_missing, 'No relative humidity', sub_name)) return
      ind_tempr = fu_index(metbuf, temperature_flag)
      if (fu_fails(ind_tempr /= int_missing, 'No temperature', sub_name)) return
      p_rel_hum => metbuf%p4d(ind_rh)
      p_tempr => metbuf%p4d(ind_tempr)
    end if

    iSpOnes = int_missing
    if (dynRules%cloud_metric == cloud_metric_ones_flag) then
      call msg(sub_name // " will use ones tracer as cell size")
      iSpOnes = select_single_species(mapConc%species, mapConc%nSpecies, &
                        & 'ones', in_gas_phase, real_missing)
      if (error .or. (iSpOnes < 1)) then
        call set_error("Couldn't find 'ones' species", sub_name)
        return
      endif
    endif

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
    do i = 1, pointers%nInSituObsID(1)
      call inject_in_situ(pointers%observationsInSitu(i), mapConc, mapMomX, mapMomY, mapMomZ, &
                & now, timestep, ithread, nthreads)
    end do
    !$OMP END PARALLEL 
    
    do i = 1, pointers%nDoseRateObsID(1)
       call inject_dose_rate(pointers%DoseRateAddition(i), pointers%observationsDoseRate(i), &
                           & mapConc, mapMomX, mapMomY, mapMomZ, mapWetdep, mapDrydep, now, timestep)
    end do

    do i = 1, pointers%nVerticalObsID(1)   ! only assimilated observations
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

  subroutine set_observations(station_path, obs_start, obs_period, obs_files, &
                            & nObsFiles_assim, nObsFiles_eval, &
                            & transport_species, optical_species, observations)
    !
    !Gets needed observations from files ans stores them to "observations"
    implicit none
    character(len=*), intent(in) :: station_path
    character(len=*), dimension(:), intent(in) :: obs_files
    type(silja_time), intent(in) :: obs_start
    type(silja_interval), intent(in) :: obs_period
    integer, intent(in) :: nObsFiles_assim, nObsFiles_eval
    type(silam_species), dimension(:), intent(in) :: transport_species
    type(silam_species), dimension(:), pointer :: optical_species
    type(observationPointers), intent(out) :: observations
    
    !character(len=worksize) :: nl_content
    integer :: uFile, status, ii, n_obs_items, n_observations, n_stations
    type(Tsilam_namelist), pointer :: nl_stations
    type(observationStation), dimension(:), allocatable :: station_list
    integer :: n_transp_species, n_opt_species
    type(grads_template) :: station_template
    character(len=fnlen) :: station_list_file
    
    if (associated(optical_species)) then
      n_opt_species = size(optical_species)
    else
      n_opt_species = 0
    end if
    
    call msg('Reading station descriptions')
    !
    ! Get the in-situ stations from the station_list
    !
    if (station_path == '') then
      n_stations = 0
    else
      uFile = fu_next_free_unit()
      call decode_template_string(station_path, station_template)
      if (error) return
      call expand_template(station_template, obs_start, station_list_file)
      if (error) return
      open(uFile, file=station_list_file, status='old', action='read', iostat=status)
      if (status /= 0) then
        call set_error('Failed to open ' // trim(station_list_file), 'set_observations')
        return
      end if
      nl_stations => fu_read_namelist(uFile, .false.)
      close(uFile)
      ! Set station indices for current subdomain
      call stations_from_namelist_ptr(nl_stations, station_list, n_stations)
      call msg('Number of found stations: ', n_stations)
      call destroy_namelist(nl_stations)
      if (error) return
    end if
        
    n_transp_species = size(transport_species)
    !
    ! Get observations. Note that for each station we require the data to be ordered in time
    ! Some stations can be for assimilation (first ones in the list of items), rest for evaluation
    ! The observation ID is the observation type (implicitly) + station code
    !
    if (nObsFiles_assim + nObsFiles_eval > 0) then
      call read_observations(obs_start, obs_period, obs_files, nObsFiles_assim, nObsFiles_eval, &
                           & station_list, n_stations, &
                           & observations%observationsInSitu, observations%nInSituObsID, &
                           & observations%observationsVertical, observations%nVerticalObsID, &
                           & observations%observationsEruption, observations%nEruptionObsID, &
                           & observations%DoseRateAddition, observations%nDoseRateObsID, &
                           & transport_species, optical_species, n_transp_species, n_opt_species)
      if (error) return
      call msg('N. of in situ observation IDs loaded   (assi,vali):', observations%nInSituObsID)
      call msg('N. of column observation IDs loaded    (assi,vali):', observations%nVerticalObsID)
      call msg('N. of eruption observation IDs loaded  (assi,vali):', observations%nEruptionObsID)
      call msg('N. of dose rate observation IDs loaded (assi,vali):', observations%nDoseRateObsID)
    else
      call msg('No observations given')
    end if
    !
    ! Stations are consumed in their entirety, both for assimilation and evaluation
    ! Set the common metadata
    !
!    observations%obs_size = fu_count_obs_size(observations)
    call count_obs_size(observations)
    ii = sum(observations%obs_size)
    allocate(observations%obs_values(ii), observations%obs_variance(ii), &
                & observations%mdl_values(ii), stat=status)
    if (fu_fails(status == 0, 'Allocate failed', 'set_observations')) return
    call collect_obs_data(observations, observations%obs_values, observations%obs_variance) 
#ifdef DEBUG
    observations%mdl_values(1:ii) = F_NAN !! Reset modlevalues
#endif

    observations%mdlCollected = .false.

    if (size(observations%obs_values) > 0) then
      call msg('Max obs value: ', maxval(observations%obs_values))
      call msg('Ave obs value: ', sum(observations%obs_values)/size(observations%obs_values))
      call msg('Number of obs (assi,vali):', observations%obs_size)
      observations%hasObservations = .true.
    else
      call msg_warning("No observations read", 'set_observations')
      observations%hasObservations = .false.
    endif
    if (allocated(station_list)) deallocate(station_list)

  end subroutine set_observations

  !************************************************************************************

    subroutine count_obs_size(obs)
      implicit none
      type(observationPointers), intent(inout) :: obs
      
      integer :: ii, iPurpose, iShift

      iShift = 0
      obs%obs_size = 0
      do iPurpose = 1, 2
        do ii = 1, obs%nInSituObsID(iPurpose)
          obs%obs_size(iPurpose) = obs%obs_size(iPurpose) &
                               & + fu_size(obs%observationsInSitu(ii + iShift))
        end do
        iShift = obs%nInSituObsID(iPurpose)
      end do

      iShift = 0
      do iPurpose = 1, 2
        do ii = 1, obs%nEruptionObsID(iPurpose)
          obs%obs_size(iPurpose) = obs%obs_size(iPurpose) &
                               & + fu_size(obs%observationsEruption(ii + iShift))
        end do
        iShift = obs%nEruptionObsID(iPurpose)
      end do

      iShift = 0
      do iPurpose = 1, 2
        do ii = 1, obs%nVerticalObsID(iPurpose)
          obs%obs_size(iPurpose) = obs%obs_size(iPurpose) &
                               & + fu_size(obs%observationsVertical(ii + iShift))
        end do
        iShift = obs%nVerticalObsID(iPurpose)
      end do

      iShift = 0
      do iPurpose = 1, 2
        do ii = 1, obs%nDoseRateObsID(iPurpose)
          obs%obs_size(iPurpose) = obs%obs_size(iPurpose) &
                               & + fu_size(obs%observationsDoseRate(ii + iShift))
        end do
        iShift = obs%nDoseRateObsID(iPurpose)
      end do

    end subroutine count_obs_size



  !************************************************************************************

  subroutine collect_model_data(obs)
    !
    ! Picks model-subdomain-observed data from observations structure into obs%mdl_values array
    ! makes exchange to form whole-MPI obs%mdl_values array
    !
    implicit none
    type(observationPointers), intent(inout) :: obs

    integer :: ind_obs, ind_start, n_values, n_values_total, max_values
    real, dimension(:), pointer :: wrk, obsvals, obsvar
    character(len=clen) :: chTmp
    logical :: ok
    character(len = *), parameter :: sub_name = 'collect_model_data'
     
    
    if ( .not. obs%hasObservations) then
        call msg("No observations to collect in "//sub_name)
        return
    endif

    if (obs%mdlCollected) then
       call set_error("Doble-collect model data", sub_name)
       return
     endif

    ind_start = 1
    n_values_total = 0
    max_values = 0

    call msg('Collecting model data')

    do ind_obs = 1, sum(obs%nInSituObsID)
      n_values = fu_size(obs%observationsInSitu(ind_obs))
      call get_data_in_situ(obs%observationsInSitu(ind_obs), values_mdl=obs%mdl_values(ind_start:))
      ind_start = ind_start + n_values
      n_values_total = n_values_total + n_values
      max_values = max(max_values, n_values)
    end do

    do ind_obs = 1, sum(obs%nEruptionObsID)
      n_values = fu_size(obs%observationsEruption(ind_obs))  
      call get_data_eruption(obs%observationsEruption(ind_obs), values_mdl=obs%mdl_values(ind_start:))
      ind_start = ind_start + n_values
      n_values_total = n_values_total + n_values
      max_values = max(max_values, n_values)
    end do
    
    do ind_obs = 1, sum(obs%nDoseRateObsID)
      n_values = fu_size(obs%observationsDoseRate(ind_obs))
      call get_data_in_situ(obs%observationsDoseRate(ind_obs), values_mdl=obs%mdl_values(ind_start:))
      ind_start = ind_start + n_values
      n_values_total = n_values_total + n_values
      max_values = max(max_values, n_values)
    end do

    do ind_obs = 1, sum(obs%nVerticalObsID)
      n_values = fu_size(obs%observationsVertical(ind_obs))
      call get_data_vertical(obs%observationsVertical(ind_obs), values_mdl=obs%mdl_values(ind_start:))
      ind_start = ind_start + n_values
      n_values_total = n_values_total + n_values
      max_values = max(max_values, n_values)
    end do

    if (fu_fails(ind_start == size(obs%mdl_values)+1, 'wrong size obs%mdl_values', sub_name)) return

    !! Do exchange and refit observations with modle data
    !! So they can be reported/evaluated by any MPI memeber 
    if (smpi_adv_tasks > 1) then
#ifdef DEBUG
      call msg("LOCAL Sum of observarions from model", sum(obs%mdl_values(1:n_values_total)))
      call msg("LOCAL obs%mdl_values(1:10)", obs%mdl_values(1:10))
#endif
      ! Synchronize model observationss
      wrk=>fu_work_array(n_values_total)
      wrk(1:n_values_total)=obs%mdl_values(1:n_values_total)
      call smpi_reduce_add(wrk(1:n_values_total), obs%mdl_values(1:n_values_total), 0, smpi_adv_comm, ok)
      call free_work_array(wrk)
      if (.not. ok) call set_error("failed MPI COMM", sub_name)

      ind_start = 1
      do ind_obs = 1, sum(obs%nInSituObsID)
        n_values = fu_size(obs%observationsInSitu(ind_obs))
        call set_data_in_situ(obs%observationsInSitu(ind_obs), obs%mdl_values(ind_start:))
        ind_start = ind_start + n_values
      end do

      do ind_obs = 1, sum(obs%nEruptionObsID)
        n_values = fu_size(obs%observationsEruption(ind_obs))  
        call set_data_eruption(obs%observationsEruption(ind_obs), obs%mdl_values(ind_start:))
        ind_start = ind_start + n_values
      end do

      do ind_obs = 1, sum(obs%nDoseRateObsID)
        n_values = fu_size(obs%observationsDoseRate(ind_obs))
        call set_data_in_situ(obs%observationsDoseRate(ind_obs), obs%mdl_values(ind_start:))
        ind_start = ind_start + n_values
      end do

      do ind_obs = 1, sum(obs%nVerticalObsID)
        n_values = fu_size(obs%observationsVertical(ind_obs))
        call set_data_vertical(obs%observationsVertical(ind_obs), obs%mdl_values(ind_start:))
        ind_start = ind_start + n_values
      end do
    endif

    !!! Separate report over in_situ TSmatrices
    if (max_values > 100) then 
      call msg("In situ scores:")
      ind_start= 1
      do ind_obs = 1, sum(obs%nInSituObsID)
        n_values = fu_size(obs%observationsInSitu(ind_obs))
        if (n_values >  20) then !! No wasting space for small observations
  !        call ooops("Here")
          if (ind_obs <= obs%nInSituObsID(1)) then 
              chTmp = 'A:'//trim(obs%observationsInSitu(ind_obs)%tag)
          else
              chTmp = 'V:'//trim(obs%observationsInSitu(ind_obs)%tag)
          endif
  !        call msg (chTmp)
  !        call msg("ind_obs Nvalues, max_values", (/ ind_obs, n_values, max_values /))
          call mod_obs_stats( chTmp, obs%mdl_values(ind_start: ind_start + n_values - 1),  &
                               &     obs%obs_values(ind_start: ind_start + n_values - 1), &
                               &   obs%obs_variance(ind_start: ind_start + n_values - 1),&
                               ind_obs == 1)
        endif
        ind_start = ind_start + n_values
      end do
    endif
#ifdef DEBUG
    call msg("Sum of observarions from model , n_values_total", sum(obs%mdl_values(1:n_values_total)), n_values_total)
    call msg("obs%mdl_values(1:10)", obs%mdl_values(1:10))
#endif
    obs%mdlCollected = .TRUE.
  end subroutine collect_model_data


  !************************************************************************************

  subroutine collect_obs_data(observations, values, variances)
    !
    !  Collects observations and variances into a single array. Should be called only once
    !  immediately after observations red. Just for convenient access to the values
    !  together with model data 

    implicit none
    type(observationPointers), intent(in) :: observations
    real, dimension(:), intent(out) :: values, variances

    integer :: ind_obs, ind_start, n_values
    
    ind_start = 1

    do ind_obs = 1, sum(observations%nInSituObsID)
      n_values = fu_size(observations%observationsInSitu(ind_obs))
      call get_data_in_situ(observations%observationsInSitu(ind_obs), & 
                     & values_obs=values(ind_start:), variance = variances(ind_start:))
      ind_start = ind_start + n_values
    end do

    do ind_obs = 1, sum(observations%nEruptionObsID)
      n_values = fu_size(observations%observationsEruption(ind_obs))
      call get_data_eruption(observations%observationsEruption(ind_obs), &
                     & values_obs=values(ind_start:), variance = variances(ind_start:))
      ind_start = ind_start + n_values
    end do
    
    do ind_obs = 1, sum(observations%nDoseRateObsID)
      n_values = fu_size(observations%observationsDoseRate(ind_obs))
      call get_data_in_situ(observations%observationsDoseRate(ind_obs), &
                     & values_obs=values(ind_start:), variance = variances(ind_start:))
      ind_start = ind_start + n_values
    end do
    
    do ind_obs = 1, sum(observations%nVerticalObsID)
      n_values = fu_size(observations%observationsVertical(ind_obs))
      call get_data_vertical(observations%observationsVertical(ind_obs), &
                     & values_obs=values(ind_start:), variance = variances(ind_start:))
      ind_start = ind_start + n_values
    end do
    
  end subroutine collect_obs_data
  

  !************************************************************************************

  subroutine get_obs_pointers(observations, obsvals, mdlvals, obsvar, n_values_total)
    !! Unified getter for values arrays
    !! If you modify them -- you get what you deserve
    implicit none
    type(observationPointers), intent(in), target :: observations
    real, dimension(:), pointer, intent(out)  :: obsvals, mdlvals, obsvar
    integer, dimension(2), intent(out) :: n_values_total
    character(len = *), parameter :: sub_name = 'get_obs_pointers'

    if (.not. observations%mdlCollected) then 
        call set_error("Observations not collected by collect_obs_data", sub_name)
        return
    endif
    
    obsvals => observations%obs_values
    mdlvals => observations%mdl_values
    obsvar => observations%obs_variance
    n_values_total = observations%obs_size
#ifdef DEBUG
    call msg(sub_name//", sum, ndata", sum(obsvals(:)), real(sum(observations%obs_size)))
#endif
    
  end subroutine get_obs_pointers
  
  !************************************************************************************

  subroutine get_localisation_all(obs_ptr, locations)
    implicit none
    type(observationPointers), intent(in) :: obs_ptr
    real, dimension(:,:), intent(out) :: locations
    
    integer :: ind_obs, ind_start, num_values

    if (fu_fails(size(locations, 1) > 1, 'locations too small 1st dim', 'get_localisation')) return
    if (fu_fails(size(locations, 2) >= sum(obs_ptr%obs_size), 'locations too small 2st dim', &
               & 'get_localisation')) return

    ind_start = 1
    num_values = 0
    do ind_obs = 1, sum(obs_ptr%nInSituObsID)
      num_values = fu_size(obs_ptr%observationsInSitu(ind_obs))
      call get_localisation_in_situ(obs_ptr%observationsInSitu(ind_obs), locations(:,ind_start:))
      ind_start = ind_start + num_values
    end do
    
    do ind_obs = 1, sum(obs_ptr%nEruptionObsID)
      num_values = fu_size(obs_ptr%observationsEruption(ind_obs))
      call get_localisation_eruption(obs_ptr%observationsEruption(ind_obs), locations(:,ind_start:))
      ind_start = ind_start + num_values
    end do

    do ind_obs = 1, sum(obs_ptr%nDoseRateObsID)
      num_values = fu_size(obs_ptr%observationsDoseRate(ind_obs))
      call get_localisation_in_situ(obs_ptr%observationsDoseRate(ind_obs), locations(:,ind_start:))
      ind_start = ind_start + num_values
    end do

    do ind_obs = 1, sum(obs_ptr%nVerticalObsID)
      num_values = fu_size(obs_ptr%observationsVertical(ind_obs))
      call get_localisation_vertical(obs_ptr%observationsVertical(ind_obs), locations(:,ind_start:))
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

    do ii = 1, sum(observations%nInSituObsID)
      call restart_in_situ(observations%observationsInSitu(ii), ifNullifyModel)
    end do
    
    do ii = 1, sum(observations%nDoseRateObsID)
      call restart_dose_rate(observations%observationsDoseRate(ii), observations%DoseRateAddition(ii), ifNullifyModel)
    end do

    do ii = 1, sum(observations%nVerticalObsID)
      call restart_vertical(observations%observationsVertical(ii), ifNullifyModel)
    end do
    observations%mdl_values(:) = F_NAN
    observations%mdlCollected = .FALSE.

  end subroutine reset_all


  !************************************************************************************

  subroutine destroy_observations(pointers)
    implicit none
    type(observationPointers) :: pointers
    
    integer :: i

    do i = 1, sum(pointers%nInSituObsID)
       call destroy(pointers%observationsInSitu(i))
    end do
    if (sum(pointers%nInSituObsID) > 0) deallocate(pointers%observationsInSitu)
    pointers%nInSituObsID = 0

    do i = 1, sum(pointers%nEruptionObsID)
       call destroy(pointers%observationsEruption(i))
    end do
    if (sum(pointers%nEruptionObsID) > 0) deallocate(pointers%observationsEruption)
    pointers%nInSituObsID = 0
    
    do i = 1, sum(pointers%nDoseRateObsID)
       call destroy(pointers%observationsDoseRate(i))
    end do
    if (sum(pointers%nDoseRateObsID) > 0) deallocate(pointers%observationsDoseRate)
    pointers%nDoseRateObsID = 0
    
    do i = 1, sum(pointers%nVerticalObsID)
       call destroy(pointers%observationsVertical(i))
    end do
    if (sum(pointers%nVerticalObsID) > 0) deallocate(pointers%observationsVertical)
    pointers%nVerticalObsID = 0

    if (sum(pointers%obs_size) > 0) deallocate(pointers%obs_values, pointers%obs_variance, pointers%mdl_values)
    
    pointers%obs_size = 0
    pointers%hasObservations = .false.
         
  end subroutine destroy_observations

  !**********************************************************************************
  !
  ! Auxiliary routines
  !
  !**********************************************************************************

  subroutine searchStationWithID(id, station_list, nstations, istation)
    implicit none

    ! Walk through the station list and return 
    ! the station with a given ID (which is assumed to be unique).

    !integer, intent(in) :: id
    character(len=*), intent(in) :: id
    type(observationStation), dimension(:), intent(in) :: station_list
    integer, intent(in) :: nstations
    integer, intent(out) :: istation

    do istation = 1, nstations
!call msg('Checking:' + station_list(istation)%id, istation)
!          if(trim(station_list(istation)%id) == 'RUMOSC')then
!            call msg(id)
!          endif
      if (station_list(istation)%id == id) return
    end do
    istation = int_missing

  end subroutine searchStationWithID

  !************************************************************************************

  subroutine dump_observation_stations(pointers, file_name)
    !
    ! Dumpr the observation data to a file, separating the assimilated / evaluation data
    !
    implicit none
    type(observationPointers), intent(in) :: pointers
    character(len=*), intent(in) :: file_name
    type(silam_species), dimension(:), pointer :: species_trn

    ! Local variables
    integer :: uFile, ind_obs, iostat, iStart, nobs
    character(len = *), parameter :: sub_name = 'dump_observation_stations'
    
    if (smpi_adv_rank /= 0) return !! Only rank  0 writes 
    call msg('Observation stations to file: ' // trim(file_name))
    
    uFile = fu_next_free_unit()
    iStart = 0
   
    nobs = sum(pointers%nInSituObsID(:))
    if (nobs > 0) then
      open(uFile, file=file_name + '.in_situ', iostat=iostat, action='write')
      if (fu_fails(iostat == 0, 'Failed to open file: ' // trim(file_name), sub_name)) return

      do ind_obs = 1, nobs
        call obs_station_to_file_in_situ(pointers%observationsInSitu(ind_obs), uFile, ind_obs == 1)
        if (error) exit
      end do
      close(uFile)
    endif
    
    nobs = sum(pointers%nDoseRateObsID(:))
    if (nobs > 0) then
      call set_error("Dose station dump not yet implemented", sub_name)
    endif
    
    nobs = sum(pointers%nVerticalObsID(:))
    if (nobs > 0) then
      call set_error("Dose station dump not yet implemented", sub_name)
    endif
    
  end subroutine dump_observation_stations

  !************************************************************************************

  subroutine dump_observations(pointers, species_trn, file_name, iPurpose)
    !
    ! Dumpr the observation data to a file, separating the assimilated / evaluation data
    ! Only master should do it. Exchange doen bu collect_observations
    !
    implicit none
    type(observationPointers), intent(in) :: pointers
    character(len=*), intent(in) :: file_name
    type(silam_species), dimension(:), pointer :: species_trn
    integer, intent(in) :: iPurpose

    ! Local variables
    integer :: uFile, ind_obs, iostat, iStart
    character(len = *), parameter :: sub_name = 'dump_observations'
    
    call msg('Observations to file: ' // trim(file_name))

    if (smpi_adv_rank /= 0) then
      call set_error("Only advection master can dump obserations", sub_name)
      return
    endif
    
    uFile = fu_next_free_unit()
    iStart = 0
    
    if(pointers%nInSituObsID(iPurpose) > 0)then
      if(iPurpose == 2) iStart = pointers%nInSituObsID(1)  ! evaluation data are after asssimilated 
      open(uFile, file=file_name + '.in_situ', iostat=iostat, action='write')
!      call msg('iostat', iostat)
      if (fu_fails(iostat == 0, 'Failed to open file: ' // trim(file_name), sub_name)) return
      write(uFile, fmt='(A)') '# In-situ observations # station_id, time, obs_data, model_data, obs_variance, obs_inj_scaling'
      do ind_obs = 1, pointers%nInSituObsID(iPurpose)
        call obs_to_file_in_situ(pointers%observationsInSitu(ind_obs + iStart), species_trn, uFile)
        if (error) exit
      end do
      close(uFile)
    endif
    
    if(pointers%nDoseRateObsID(iPurpose) > 0)then
      if(iPurpose == 2) iStart = pointers%nDoseRateObsID(1)  ! evaluation data are after asssimilated 
      open(uFile, file=file_name + '.dose_rate', iostat=iostat, action='write')
      if (fu_fails(iostat == 0, 'Failed to open file: ' // trim(file_name), sub_name)) return
      do ind_obs = 1, pointers%nDoseRateObsID(iPurpose)
        call obs_to_file_in_situ(pointers%observationsDoseRate(ind_obs + iStart), species_trn, uFile)
        if (error) exit
      end do
      close(uFile)
    endif
    
    if(pointers%nVerticalObsID(iPurpose) > 0)then
      if(iPurpose == 2) iStart = pointers%nVerticalObsID(1)  ! evaluation data are after asssimilated 
      open(uFile, file=file_name + '.vertical', iostat=iostat, action='write')
      if (fu_fails(iostat == 0, 'Failed to open file: ' // trim(file_name), sub_name)) return
      do ind_obs = 1, pointers%nVerticalObsID(iPurpose)
        call obs_to_file(pointers%observationsVertical(ind_obs + iStart), species_trn, uFile)
        if (error) exit
      end do
      close(uFile)
    endif
    
  end subroutine dump_observations

  !*******************************************************

  subroutine parse_obs_param(param, scale_to_silam, ifmmr)
    !! Unit conversion of text-style observations
    !! Assumes observed mass unit is basic one, tells if it is cnc or mmr
    !! Treats concentration and mixing-ratio observations
    !! and converts them to "per m3" or "per kg"
     implicit none
     character(len = *), intent(in) ::  param
     real, intent(out) :: scale_to_silam  !! Scale to silam unit (per kg or per m3)
     logical, intent(out) :: ifmmr !! "true" == "per kg of air"
    !the air density as universally defined by the EC, i.e. 1.2041 kg/m**3 at 20 deg C and 1013.25 mbar                                          
    real, parameter :: air_dens_EC =  1.2041

    !the air density as universally defined by the US EPA and the WHO, i.e. 1.2041 kg/m**3 at 25 deg C and 1013.25 mbar                          
    real, parameter :: air_dens_EPA =  1.1839
     character(len = *), parameter :: sub_name = 'parse_obs_param'
     !!     (/'cnc    ','cncEU  ','cncWHO ','perkilo','permole'/))
   
     scale_to_silam = real_missing
     if  (param == 'cnc') then !!! Concentration
        ifmmr = .FALSE.
        scale_to_silam =  1. !! per m3
     else
        ifmmr = .TRUE.
      select case (param)
        case ('cncEU')
          scale_to_silam = 1./air_dens_EC  !! 1/air_dens
        case ('cncWHO')
          scale_to_silam = 1./air_dens_EPA  
        case ('perkilo')
          scale_to_silam = 1.
        case ('permole')
          scale_to_silam = 1./molecular_weight_air  !! per air molar mass
        case default
          call set_error("Strange param '"//trim(param)// &
              &"', can be one of: cnc, cncEU, cncWHO, perkilo, permole", sub_name)
          
      end select 
     endif
 end subroutine parse_obs_param

  !***************************************************************************************
 
  subroutine parse_obs_unit(obsvar, obsunit, obs_amt_unit, scale_to_silam, ifmmr)
    !! Parses obs unit (comes from NetCDF time series),
    !! and gets the basic amount unit and scaling to it
    !! moreflexible wrapper for  parse_obs_param 
     implicit none
     character(len = *), intent(in) ::  obsvar, obsunit
     character(len = *), intent(out) :: obs_amt_unit !! Needed for further massmap tp cocktail conversion
     real, intent(out) :: scale_to_silam  !! Scale of obs to silam unit (per kg or per m3)
     logical, intent(out) :: ifmmr !! "true" == "per kg of air"

     integer :: iTmp
     real :: fFactor
     character(len = clen) :: obsuntitstuff, obsunitair, param
     character(len = 1) :: sepchar
    character(len=*), parameter :: sub_name = 'parse_obs_unit'

     !! slash or space separated unit
     iTmp = index(obsunit,"/")
     if (iTmp  == 0) iTmp = index(obsunit," ")
     
     if (iTmp > 1) then 
       obsuntitstuff = obsunit(1:iTmp-1)
       obsunitair = obsunit(iTmp+1:)
       sepchar = obsunit(iTmp:iTmp)
     else
       call set_error("Strange obsunit '"//trim(obsunit)//"'", sub_name)
       return
     endif

     !! What is the quantity?
     !! Can be cnc with e.g. units of EUug/m3
     if (any(obsvar== (/'cnc   ','cncEU ','cncWHO'/))) then 
          param = obsvar
          if (obsuntitstuff(1:2) == "EU")  then
              param = "cncEU"
              obsuntitstuff = obsuntitstuff(3:)
          elseif (obsuntitstuff(1:3) == "WHO")  then
              param = "cncWHO"
              obsuntitstuff = obsuntitstuff(4:)
          endif
          if ( ( sepchar=='/' .and. all(obsunitair /= (/'m3  ','m^3 ','m**3'/) ) ) .or. &
             & ( sepchar==' ' .and. all(obsunitair /= (/'m-3  ','m^-3 ','m**-3'/) ) ) ) then
             call set_error("Unparsable unit of concentration //'"//trim(obsunit)//"'", sub_name)
             return
          endif
     elseif (obsvar == "vmr") then
       if (.not. any(obsunitair == (/'mole', 'mol '/))) then
          call set_error("Unparsable unit of vmr //'"//trim(obsunit)//"'", sub_name)
       endif
       param = "permole"
     elseif (obsvar == 'mmr') then 
       if (obsunitair /= 'kg') then
          call set_error("Unparsable unit of mmr //'"//trim(obsunit)//"'", sub_name)
        endif
        param = "perkilo"
     else
        call set_error("Unknown variable name //'"//trim(obsunit)//"'", sub_name)
     endif

     !! Conversion from e.g ug
     obs_amt_unit = fu_SI_unit(obsuntitstuff) !! Shold be mole kg or something else
     fFactor =  fu_conversion_factor(obsuntitstuff, obs_amt_unit)

     call  parse_obs_param(param, scale_to_silam, ifmmr)
     if (error) return

     scale_to_silam = scale_to_silam * fFactor

 end subroutine parse_obs_unit

  !***************************************************************************************
  
  integer function fu_number_of_observed_cells(obs_pointers)
    !
    ! A funny function that says how many grid points are filled with observations.
    ! Used in low_mass_threshold setting routine
    ! A silent assumption is that each observationID is in own grid cell
    !
    implicit none

    type(observationPointers), intent(in) :: obs_pointers
  
    fu_number_of_observed_cells = sum(obs_pointers%nInSituObsID &
                                  & + obs_pointers%nVerticalObsID &
                                  & + obs_pointers%nDoseRateObsID)
    
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

    type(observationStation), dimension(2) :: station_list
    type(inSituObservation), dimension(:), pointer :: in_situ_list, dose_rate_list
    type(t_dose_rate_obs_addition), dimension(:), pointer :: addition_list
    integer, dimension(2) :: num_in_situ, num_vert_obs, num_dose_rate, num_eruption_obs
    character(len=fnlen), dimension(:), pointer :: obs_items
    character(len=fnlen) :: obs_list_path
    character(len=*), parameter :: tmpfilename = 'tmp.dat'
    type(silja_grid) :: grid
    integer :: obs_size, ind_obs, num_val, nObsItems
    logical :: failure
    real, dimension(:), pointer :: dataptr
    real, parameter :: dummy_obs_val = 1.0
    real :: conv_to_mol
    type(t_vertical_observation), pointer, dimension(:) :: null_vert_obs_lst
    type(silam_species), dimension(:), pointer :: null_species_lst
    type(t_eruptionObservation), pointer, dimension(:) :: null_eruption_obs_lst
    type(silja_time) :: obs_start
    type(silja_interval) :: obs_period
    character(len=*), parameter :: sub_name = 'test_read_timeseries'

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
    call read_observations(obs_start, obs_period,  obs_items, nObsItems, 0, &
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
    call read_observations(obs_start, obs_period, obs_items, nObsItems, 0, &
                         & station_list, size(station_list), in_situ_list, num_in_situ, &
                         & null_vert_obs_lst, num_vert_obs, &
                         & null_eruption_obs_lst, num_eruption_obs, &
                         & addition_list, num_dose_rate, &
                         & transp_species, null_species_lst, num_transp_sp, 0)

    call msg('series read:', num_in_situ)
    if (fu_fails(sum(num_in_situ) == 2, 'Failed to read 2 obs', sub_name)) return 

    call check_observations(in_situ_list, obs_size, station_list)

    call msg('Daily obs')
    call make_observations(tmpfilename, station_list, obs_start, obs_start+obs_period, &
                         & 'SO2', one_day, obs_size)
    call read_observations(obs_start, obs_period, obs_items, nObsItems, 0, &
                         & station_list, size(station_list), in_situ_list, num_in_situ, &
                         & null_vert_obs_lst, num_vert_obs, & 
                         & null_eruption_obs_lst, num_eruption_obs, &
                         & addition_list, num_dose_rate, &
                         & transp_species, null_species_lst, num_transp_sp, 0)
    call msg('series read:', num_in_situ)
    if (fu_fails(sum(num_in_situ) == 2, 'Failed to read 2 obs', sub_name)) return 

    call check_observations(in_situ_list, obs_size, station_list)
    
    
    !call cleanup(tmpfilename)

  contains
    
    subroutine check_observations(obs_list, obs_size, station_list)
      implicit none
      type(insituObservation), dimension(:), intent(in) :: obs_list
      type(observationStation), dimension(:), intent(in) :: station_list
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
        failure = fu_fails(obs_list(ind_obs)%station(1)%id == station_list(ind_obs)%id, &
             & 'wrong id', sub_name)
        call get_data_in_situ(in_situ_list(ind_obs), values_obs=dataptr)
        if (fu_fails(all(dataptr(1:obs_size) .eps. expect_values), 'Wrong values', sub_name)) then
          call msg('ind_obs:', ind_obs)
          print *, 'Expected:', expect_values(1:obs_size)
          print *, 'Obtained:', dataptr(1:obs_size)
        end if
        
      end do
      
      call free_work_array(dataptr)
      
    end subroutine check_observations

    subroutine make_stations(station_list, how_many)
      implicit none
      type(observationStation), dimension(:), intent(out) :: station_list
      integer, intent(in) :: how_many
      
      integer :: ind_station
      type(observationStation), pointer :: stationptr
      character(clen) :: code, name
      
      do ind_station = 1, how_many
        write(code, fmt='("A",I0)') ind_station
        write(name, fmt='("station_", I0)') ind_station
        station_list(ind_station) = fu_initObservationStation(trim(code), trim(name), &
                                             & 30.0 + ind_station, 55.0 + ind_station, 2.0)
      end do
      
    end subroutine make_stations

    subroutine make_observations(filename, station_list, time_start, time_end, cockt_name, duration, &
                               & num_val)
      implicit none
      character(len=*), intent(in) :: filename, cockt_name
      type(observationStation), dimension(:), intent(in) :: station_list
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
          write(unit, fmt='(A, A, I5, I3, I3, I3, 3G12.3)') station_list(ind_station)%id, cockt_name, &
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

