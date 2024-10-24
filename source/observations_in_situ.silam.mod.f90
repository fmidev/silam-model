module observations_in_situ
  !
  ! This module defines subroutines and datatypes for time-series type observations of
  ! either in-situ or column-integral type. The observations act on the model fields
  ! either in forward (observe) or adjoint (inject) mode. In the forward mode, the
  ! modelled observations are stored, in the adjoint mode, the discrepancy y - Hx between
  ! modelled and real observations is injected as a forcing.
  !
  ! The cost-function related stuff is no longer here. Instead, the observed values are
  ! packed into a vector and sent to the upper level subroutines.
  !
  ! A numerical trick: the injected mass is multiplied with a number >>1 and divided when
  ! the gradient is collected. This is because the adjoint runs have no source terms, and
  ! hence the much lower concentrations would end up negligible compared to the low-mass
  ! threshold. This is needed even though the observations and variances given in units
  ! /m3 or /m2 are converted into cell or cell-area integrated values.

  use dispersion_server !chemistry_manager !source_apportionment
  use optical_density
!  use perturbations

  implicit none

!  private

  public inject_in_situ
  public observe_in_situ
  public observe_eruption

  public fu_initInSituObservation1
  public fu_initObservationStation
  public fu_initEruptionObservation
  public fu_station_label
  public fu_lat
  public fu_lon

  !public destroyObservation
  public destroy
  public defined

  public observationStation
  public inSituObservation
  public t_eruptionObservation
  public test_in_situ
  public obs_to_file_in_situ
  public obs_station_to_file_in_situ
  public fu_size
  public get_data_in_situ
  public get_data_eruption
  public report
  public get_localisation
  public restart_in_situ

  !  public fu_init_aot_observation

  ! PRIVATE subroutines start here
  !

  private fu_getValueAtStation
  private fu_size_in_situ
  private fu_station_defined


  interface destroy
     module procedure destroyInSitu
     module procedure destroyEruption
  end interface

  interface fu_size
     module procedure fu_size_in_situ
     module procedure fu_size_eruption
  end interface

  interface defined
     module procedure fu_station_defined
  end interface

  interface get_localisation
     module procedure get_localisation_in_situ
     module procedure get_localisation_eruption
  end interface

  interface report
     module procedure report_obs
  end interface
  integer, parameter, public :: STATION_ID_LENGTH = 35
  integer, parameter, public :: STATION_LABEL_LENGTH = 64

  integer, parameter, public :: constant_variance = 1, variable_variance = 2

  real, public, parameter :: observation_scaling = 1e5


  !**********************************************************************************
  !
  ! Type definition for an in situ observation
  !
  !**********************************************************************************
  type observationStation
     !
     ! Observation station: defines a site where AOD or concentration
     ! observations are made.
     !
     character(len=STATION_ID_LENGTH) :: ID
     character(len=STATION_LABEL_LENGTH) :: label
     real :: lat, lon !
!     integer, dimension (2) :: ix_dispersion, iy_dispersion !Indices 4 interp in this_domain - dispersion_grid !!!!
!     integer, dimension (2) :: ix_whole_dispersion, iy_whole_dispersion !Indices 4 interp in this_domain - dispersion_grid !!!!
     integer :: ix_dispersion, iy_dispersion !Indices of the cell in this_domain - dispersion_grid !!!!
     integer :: ix_whole_dispersion, iy_whole_dispersion !Indices 4 interp in this_domain - dispersion_grid !!!!
     real :: xLoc_dispGrid, yLoc_dispGrid              ! local shift of station from the cell centre: (-0.5, 0.5)
!     real :: xpos, ypos !! 0..1 relative position between two coordinates in x and y
     logical :: inThisSubdomain !! Our massmap contributes to this station
     real :: cell_area
     integer :: ind_lev = 1
     logical :: defined = .false.
  end type observationStation

  type(observationStation), parameter, public ::  station_missing = &
        & observationStation('', '', real_missing, real_missing, &
                            & int_missing, int_missing, int_missing, int_missing, &
                            & real_missing, real_missing, .false.,  &
                            & real_missing,  int_missing, &
                            & .false.)
  ! This structure defines a time series type observation from a single station. Each
  ! station has a fixed postion and measures a single substance only.

  type inSituObservation
     ! something to identify this observation in the dump file.
     character(len=clen) :: tag = '???'
     integer :: iInterpolation =  nearest_point
     ! for each measurement:
     type(silja_time), dimension(:), allocatable :: endTimes   !startingTimes ! => null()
     type(silja_interval), dimension(:), allocatable :: durations
     ! the real observed data, and that "observed" from the model
     real, dimension(:, :), allocatable :: obsData, modelData !!! silamMassUnit/kgAir or silamMassUnit/m3
                                                                  !! depending of ifmmr
     real, dimension(:,:), allocatable :: variance  !! (nt,nst) Square of Data unit
     real, dimension(:,:), allocatable ::inject_scaling  !! (nt,nst)  Scaling for inject depending on 
                                    ! unit and massmap metrics 
     logical :: ifmmr = .false. !! true: Treat observation as mixing ratio (silamMassUnit/kgAir)
                                !! false: Treat observation as concentrations (silamMassUnit/m3)

     type(observationStation), dimension(:), allocatable :: station
     logical :: inThisSubdomain !! True if any station belongs to this subdomain
     integer :: ilevel = 1

     ! The mapping from (possibly several) transport to the observed species. Also
     ! contains factors for unit conversion.
     integer, dimension(:), allocatable :: ind_obs2transp
     real, dimension(:),allocatable :: scale_transp2obs
     integer :: num_obs_species

     ! the length of data and modelData
     integer :: dataLength  !! valid observed data points in array. Those will be exposed to assimilator
     ! The number of obsData that are not real_missing
     ! For single-station observation it should be 
     integer :: nStations, nTimes !! Array size 

     ! the counter. As the model proceeds, and observer/inject routines are
     ! called, the observation keeps a count of the index of current measurement.
     ! In forward run, it runs from 1 to nTimes, and the opposite in the adjoint
     ! case.
     integer :: observationTimestep = 1
     integer :: time_direction
     ! Volume of the grid cell, for conversion from mass to concentration (only if geometrical metric used).
     real, dimension(:), allocatable :: cell_volume
  end type inSituObservation

  type t_eruptionObservation
    character(len=substNmLen) :: tag = '???'
    type(silja_time), dimension(:), allocatable :: endTime   !startingTimes ! => null()
    type(silja_interval), dimension(:), allocatable :: duration
    real, dimension(:), allocatable :: model_eruption_height, obs_eruption_height, variance
    real, dimension(:), allocatable :: lat, lon
  end type t_eruptionObservation


  !***********************************************************************************************

contains

  ! Initialization. Returns a new observation object with given data,
  ! locations, and time windows. The returned object is ready for use
  ! at any time.  Array pointers endTimes and durations should be
  ! allocated and dipersionGrid should be defined in advance.
  !! WARNING: Initializes a single-station observation


  !***********************************************************************************************


  function fu_initInSituObservation1(endTimes, durations, obsData, variance, &
                                  & n_values, variance_flag, ifmmr,&
                                  & ind_obs2transp, scale_transp2obs, nSpecies, &
                                  & station, &
                                  & vertical, tag) &
                                  & result(newObservation)
    implicit none
    ! input
    integer, intent(in) :: variance_flag, n_values, nSpecies
    type(silja_time), dimension(n_values), intent(in)  :: endTimes
    type(silja_interval), dimension(n_values), intent(in) :: durations
    real, dimension(n_values), intent(in) :: obsData, variance
    integer, dimension(:), intent(in) :: ind_obs2transp
    real, dimension(:), intent(in) :: scale_transp2obs
    logical, intent(in) :: ifmmr
    ! the durations are currently assumed to be given in hours. True,
    ! silja_interval would make sense in this context.
    type(observationStation), intent(in) :: station
    type(silam_vertical), intent(in) :: vertical
    character(len=*), intent(in) :: tag

    type(inSituObservation) :: newObservation

    !
    ! This sub only sets a single-station observation
    ! for multistation observations NetCDF input of TSMartix to be used

    ! local
    integer :: status, isp_obs

    allocate(newobservation%endTimes(n_values), newObservation%durations(n_values), &
           & newObservation%obsData(1,n_values), newObservation%variance(1,n_values), &
           & newObservation%modelData(1,n_values), newObservation%inject_scaling(1,n_values), &
           & newObservation%ind_obs2transp(nSpecies), &
           & newObservation%station(1), &
           & newObservation%scale_transp2obs(nSpecies),stat=status)
    if(fu_fails(status == 0, 'Allocate failed', 'fu_initInSituObservation1'))return

    newObservation%tag = tag
    newObservation%endTimes = endTimes(1:n_values)
    newObservation%durations = durations(1:n_values)
    newObservation%obsData(1,:) = obsData(1:n_values)
    newObservation%dataLength = n_values
    newObservation%nTimes = n_values 
    newObservation%nStations = 1
    newObservation%station(1) = station
    newObservation%inThisSubdomain = station%inThisSubdomain !! Here we have the only station
    newObservation%ifmmr = ifmmr

    newObservation%modelData(:,:) = 0.0
    newObservation%inject_scaling(:,:) = 0.0 !! They will be filled

    select case(variance_flag)
    case (variable_variance)
      newObservation%variance(1,:) = variance(1:n_values)
    case (constant_variance)
      newObservation%variance(1,:) = variance(1)
    case default
      call set_error('Strange variance_flag', 'fu_initInSituObservation1')
      return
    end select


    newObservation%num_obs_species = nSpecies
    newObservation%ind_obs2transp(1:nSpecies) = ind_obs2transp(1:nSpecies)
    newObservation%scale_transp2obs(1:nSpecies) = scale_transp2obs(1:nSpecies)


    newObservation%cell_volume = newObservation%station%cell_area * &
                                            & fu_layer_thickness_m(fu_level(vertical, 1))

    newObservation%obsData(1,:) = obsData(1:n_values)

  end function fu_initInSituObservation1

  function fu_initEruptionObservation(endTime, duration, obsData, variance, &
                                  & lat, lon, n_values, tag) &
                                  & result(newObservation)
    implicit none
    ! input
    type(silja_time), dimension(:), intent(in)  :: endTime
    type(silja_interval), dimension(:), intent(in) :: duration
    real, dimension(:), intent(in) :: obsData, variance, lat, lon
    ! the durations are currently assumed to be given in hours. True,
    ! silja_interval would make sense in this context.
    integer, intent(in) :: n_values
    character(len=*), intent(in) :: tag

    type(t_eruptionObservation) :: newObservation

    integer :: status

    allocate(newobservation%endTime(n_values), newObservation%duration(n_values), &
           & newObservation%obs_eruption_height(n_values), newObservation%variance(n_values), &
           & newObservation%lon(n_values), newObservation%lat(n_values), &
           & newObservation%model_eruption_height(n_values), stat=status)
    if (status /= 0) then
      call set_error('Allocate failed', 'fu_initInSituObservation1')
      return
    end if

    newObservation%lat = lat(1:n_values)
    newObservation%lon = lon(1:n_values)

    newObservation%endTime = endTime(1:n_values)
    newObservation%duration = duration(1:n_values)
    newObservation%obs_eruption_height = obsData(1:n_values)
    newObservation%variance = variance(1:n_values)

    newObservation%model_eruption_height = 0.0

  end function fu_initEruptionObservation

  !***********************************************************************************************

  function fu_initObservationStation(ID, label, lon, lat, z) result(new_station)

    ! Create an observation station. Dispersion grid is needed for computing the
    ! projection from the geo grid. No interpolation is implemented due to
    ! trouble with domain split and other consustency issues.
    ! Since grids in SILAM are inaccurate we project always to wholeMPIdispersion_grid
    ! and then subtract offsets
    implicit none
    character(len=*), intent(in) :: id, label
    real, intent(in) :: lat, lon, z ! z doesn't do anything yet

    type(observationStation) :: new_station

    real :: x,y, xpos, ypos
    integer :: status, ix, iy, ixloc, iyloc, i,j
    INTEGER :: nx, ny, offx, offy, gnx, gny

    logical :: ifLglob

    character(len=fnlen) :: str

    call smpi_get_decomposition(nx, ny, offx, offy, gnx, gny)
    call project_point_to_grid(lon, lat, wholeMPIdispersion_grid, x, y)
    ifLglob = fu_ifLonGlobal(wholeMPIdispersion_grid)
    ix = floor(x)
    iy = floor(y)
    xpos = x - real(ix)    ! position with regard to the grid cell centre
    ypos = y - real(iy)

    !!! For nearest neighbour instead of linear: Some x2 faster assimilation
    !!! due to much better convergence :: inject also uses nearest neighbour!!!!
    !--------------------
   !        call msg_warning("Trying to emulate nearest neighbour","fu_initObservationStation")
!      xpos=real(nint(xpos))
!      ypos=real(nint(ypos))
    !------------------


    if ( ( (ix < 1 .or. ix >=  gnx) .and. (.not. ifLglob)) .or. &
        & iy < 1 .or. iy >= gny  ) then
      ! Only stations within the grid are considered
      ! Stations at poles will be killed:
      !          To handle them we need at least information on the domain closure
      ! Information that station is out of grid is important in normal run too. It is not an excuse to break the table format
#ifdef DEBUG
      write(str, fmt='("OUT ", A10, 1x, F8.2, F8.2," ", A20, I6, I6)') id, lat, lon, trim(label), ix, iy
#else
      write(str, fmt='("OUT of grid: ", A10, 1x, F8.2, F8.2," ", A20, I6, I6)') id, lat, lon, trim(label), ix, iy
#endif
      call msg(str)
      new_station = station_missing
    else
      ixloc = ix - offx
      iyloc = iy - offy
      new_station= observationStation(id, label, lat, lon, &
                                            & ixLoc, iyLoc, ix, iy, &
!                                            & (/ixloc, ixloc+1/), (/iyloc, iyloc+1/), & !ix_dispersion, iy_dispersion
!                                            & (/ix,ix+1/), (/iy,iy+1/), &
                                            & xpos, ypos, &
                                            & .True., &
                                            & real_missing, 1, &! cell_volume, ind_lev
                                            & .false.) !not yet defined

      ! iy_whole_dispersion are okay, ix_whole_dispersion has to be adjusted
      !  only for lon_global "in the cut" situation
      !  For ix_dispersion and iy_dispersion  need to check if it is still our domain
      !
!      do i=1,2
!       if (ifLglob) then
!          new_station%ix_dispersion(i) = modulo(new_station%ix_dispersion(i)-1, gnx) + 1
!          new_station%ix_whole_dispersion(i) = new_station%ix_dispersion(i) !! No offsets for global
!        else
!          ixloc = new_station%ix_dispersion(i)
!          if (ixloc < 1 .or. ixloc > nx) new_station%ix_dispersion(i) = int_missing
!        endif
!        iyloc = new_station%iy_dispersion(i)
!        if (iyloc < 1 .or. iyloc > ny) new_station%iy_dispersion(i) = int_missing
!      enddo
      if (ifLglob) then
        new_station%ix_dispersion = modulo(new_station%ix_dispersion-1, gnx) + 1
        new_station%ix_whole_dispersion = new_station%ix_dispersion !! No offsets for global
      else
        ixloc = new_station%ix_dispersion
        if (ixloc < 1 .or. ixloc > nx) new_station%ix_dispersion = int_missing
      endif
      iyloc = new_station%iy_dispersion
      if (iyloc < 1 .or. iyloc > ny) new_station%iy_dispersion = int_missing

      new_station%cell_area = fu_dx_cell_m(wholeMPIdispersion_grid, nint(x), nint(y)) * &
                           &  fu_dy_cell_m(wholeMPIdispersion_grid, nint(x), nint(y))

      ! Check if the station is completely out
      if(new_station%ix_dispersion == int_missing .or. new_station%iy_dispersion == int_missing)then
        new_station%inThisSubdomain = .False.
      endif
#ifdef DEBUG
        write(str, fmt='("IN  ", A10, 1x, F8.2, F8.2," ", A20, I6, I6, " My ix,iy:", I8, I8, " Whole ix,iy:", I8, I8)') &
              & id, lat, lon, trim(label), ix, iy, new_station%ix_dispersion, new_station%iy_dispersion, &
              & new_station%ix_whole_dispersion, new_station%iy_whole_dispersion
        call msg(str)
#endif
      new_station%defined = .true.
    end if

  end function fu_initObservationStation

  !************************************************************************************
  ! Writing the observations into text files:
  !
  subroutine obs_to_file_in_situ(obs, species_trn, file_unit)
    implicit none
    type(inSituOBservation), intent(in) :: obs
    type(silam_species), dimension(:), intent(in) :: species_trn
    integer, intent(in) :: file_unit

    integer :: ii, iT, iSt
    character (len=worksize_string) :: sp

    sp = fu_str(species_trn(obs%ind_obs2transp(1)))
    do ii = 2, size(obs%ind_obs2transp)
      sp = sp + "__" + fu_str(species_trn(obs%ind_obs2transp(ii)))
    enddo

    do iSt = 1, obs%nStations
      do iT = 1, obs%nTimes
        if (obs%obsData(iSt,iT) == real_missing) cycle
        write(file_unit, fmt='(2(A,1x), I4, 2I3, 1x, 2(F9.2,1x), 4G12.3)') &
                    & trim(obs%station(iSt)%id), trim(sp), &
                    & fu_year(obs%endTimes(iT)), fu_mon(obs%endTimes(iT)), fu_day(obs%endTimes(iT)), &
                    & fu_hour(obs%endTimes(iT)) + fu_min(obs%endTimes(iT)) / 60. + &
                                                & fu_sec(obs%endTimes(iT)) / 3600., &
                    & fu_hours(obs%durations(iT)), &
                    & obs%obsData(iSt,iT), obs%modelData(iSt,iT), obs%variance(iSt,iT), obs%inject_scaling(iSt,iT)
      enddo !iT
    enddo !iSt
  end subroutine obs_to_file_in_situ

  !************************************************************************************
  subroutine obs_station_to_file_in_situ(obs, file_unit,ifHead)
    implicit none
    type(inSituOBservation), intent(in) :: obs
    logical, intent(in) :: ifHead
    integer, intent(in) :: file_unit
    integer :: iSt

    if (ifHead) write(file_unit, fmt='(A)') '# ID lon lat label'

    do iSt = 1, obs%nStations
      write(file_unit, fmt='(A35,X,F10.5, X,F10.5,X,A)') &
                 & obs%station(iSt)%id, obs%station(iSt)%lon, obs%station(iSt)%lat, obs%station(iSt)%label
    enddo

  end subroutine obs_station_to_file_in_situ


  !**********************************************************************************
  !
  ! Inject and observe
  !
  !**********************************************************************************

  subroutine inject_in_situ(obs, map_c, map_px, map_py, map_pz, now, timestep, iThread, nThreads)

    ! Inject all obs that happened between now and now+timestep
    ! timestep should be non-positive                              
    implicit none
    type(inSituObservation), intent(inout) :: obs
    type(Tmass_map), intent(inout) :: map_c, map_px, map_py, map_pz !!!!Centers of mass here
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep
    integer, intent(in) :: iThread, nThreads

    integer :: iTobs, iSeconds, ind_src, ix_obs, iy_obs, ilev_obs, isp_obs, iSt
    real ::  obs_val, mdl_val, inject_mass, mass_fract
    real ::  mOld, mAdd ! mass in the map prior to inhection and it increment
    integer :: isp_transp
    type(silja_interval) :: overlap
    real ::  fObsEnd,fObsDur,fOverStart,fOverEnd, fWeight
    type(silja_time) :: obs_start, step_start !Earliest time of the step (timestep+now)
    character(len=*), parameter :: sub_name = 'inject_in_situ'

    if (.not. obs%inThisSubdomain) return ! All statons from foreign subdomain

    ix_obs = obs%station(1)%ix_dispersion
    if (obs%nstations==1 .and.  mod(ix_obs, nThreads) /= iThread) return !!!speedup single-station injection

    if (obs%time_direction == forwards) then
      obs%observationTimestep = obs%nTimes
      obs%time_direction = backwards
    end if

    iSeconds = nint(fu_sec(timestep))
    if (iSeconds > 0 ) then 
      call set_error("Injecting on non-negative timesteps!", sub_name)
      return
    endif

    !!
    step_start = now + timestep !negative timesep! earliest time to inject
    ! fast-backward through later observations
    do iTobs = obs%observationTimestep,1,-1
      obs_start =  obs%endTimes(iTobs) - obs%durations(iTobs)
      if (obs_start <= now) exit  !observation might affect current time step
    enddo

    if (iTobs < 1) return ! no more obs

    obs%observationTimestep = iTobs !save for future use

    
    do iTobs = obs%observationTimestep,1,-1

      fWeight = 0.
      if ( iSeconds == 0) then  !! instant massmap (e.g. 3Dvar) -> only instant obs. possible
        if (obs%endTimes(iTobs) < step_start) exit  !observation is fully before 
        if (obs%endTimes(iTobs) ==  now) then
            if (.not. obs%durations(iTobs) == zero_interval) then
              call set_error("Non-instant observation from instant step", sub_name)
            endif
            fWeight = 1.

        endif
      else  
        !finite  time steps
        if (obs%endTimes(iTobs) <= step_start) exit  !observation is fully before or at 

       if ( fu_sec(obs%durations(iTobs)) < 1.) then
         !instant
         fWeight = 1.
       else

       ! normalize all times so that (now+timestep:now] maps to (-1 : 0] 
       ! note that now+timestep is excluded



         fObsEnd = (fu_sec(obs%endTimes(iTobs) - now)) / (-iSeconds)
         if (fObsEnd < -1.) call set_error("Must never happen",sub_name)
         fObsDur = fu_sec(obs%durations(iTobs)) / (-iSeconds)
             !                                                                !
             !                          |                                     !
             !                          |                                     !
             !                  ___ __1_|                                     !
             !                 |  |  |  |                                     !
             !                 |  |  |  |                                     !
             !                 |  |  |  |                                     !
             !                 |  |  | 0|                                     !
             !                ______________________                          !
             !                -1  L  R  0        1     -> Normalized time     !
             ! overlap of [-1 : 1] and normalized observation range [L:R]          !
         fOverStart = max(-1.,fObsEnd - fObsDur) 
         fOverEnd   = max(fObsEnd, 0.)
         fWeight = fOverEnd - fOverStart
         fWeight = fWeight / fObsDur ! Fraction of the integral
        endif !instant obs
      endif  !!!iSeconds == 0

      !
      ! Nothing to add
#ifdef DEBUG      
      if (.not. fWeight >= 0) then
        call set_error("Negative weight", sub_name)
      endif
#endif 

      if ( fWeight < 1e-5) cycle

      
      do iSt = 1, obs%nStations
        if (.not. obs%station(iSt)%inThisSubdomain) cycle 
        !Actual value to the massmap


        ix_obs = obs%station(iSt)%ix_dispersion
        if (mod(ix_obs, nThreads) /= iThread) cycle !!!Make it possible for parallel injection
        iy_obs = obs%station(iSt)%iy_dispersion
        ilev_obs = obs%station(iSt)%ind_lev
        !!! separate positive, zeros and negative
        obs_val = obs%obsData(iSt,iTobs)
        if (obs_val == real_missing) cycle
        mdl_val = obs%modelData(iSt,iTobs)
        if (abs(obs_val) <  5e-15) then  ! nasty: dimentional variable
          ind_src = DA_ZERO
        else if (obs_val > mdl_val) then
          ind_src = DA_NEGATIVE
          fWeight = fWeight*DA_NEGT_COEF_OBS
        else
          ind_src = DA_POSITIVE
        end if

        inject_mass = fWeight * (mdl_val - obs_val)  * observation_scaling /obs%variance(iSt, iTobs)
        !! Handles inconsistency between massmap metrics and observed quantity
        inject_mass = inject_mass * obs%inject_scaling(iSt, iTobs)
       !!call msg("INJ: "//trim(obs%station(iSt)%ID), inject_mass)
#ifdef DEBUG_OBS
obs_start =  obs%endTimes(iTobs) - obs%durations(iTobs)
call msg("")
call msg(trim(fu_str(iTobs))//' obs_start: ' // fu_str(obs_start) // ', obs_end: ' // fu_str(obs%endTimes(iTobs) ) // ', now: '// fu_str(now))
call msg('Station ' // trim(obs%station(iSt)%id)//': obs , model',  obs%obsData(iSt,iTobs), obs%modeldata(iTobs,iTobs))
call msg('Fraction to inject:', fWeight )

#endif 

          !call msg('weight, var', weight, obs%variance(iTobs))
          do isp_obs = 1, obs%num_obs_species
            isp_transp = obs%ind_obs2transp(isp_obs)
            mass_fract = obs%scale_transp2obs(isp_obs)
            if (mass_fract > 0.0) then

              !call msg('Injecting mass in:' + fu_str(ix_obs) + ',' &
              !       & + fu_str(iy_obs), mass_fract*weight*inject_mass)

              map_c%ifColumnValid(ind_src, ix_obs, iy_obs) = .true.
              map_c%ifGridValid(1, ind_src) = .true.
              mOld = map_c%arm(isp_transp,ind_src,ilev_obs,ix_obs, iy_obs) !What is there
              mAdd = mass_fract*inject_mass !increment
              map_c%arm(isp_transp,ind_src,ilev_obs,ix_obs, iy_obs)  = mOld + mAdd
            end if   ! mass_fract > 0
         end do   ! obs_species
      enddo !! iSt

    end do ! loop over measurements

  end subroutine inject_in_situ


  !************************************************************************

  subroutine observe_in_situ(obs, map_c, map_cmX, map_cmY, disp_buf_ptr, cloud_metric_flag, iSpOnes, now, timestep)
    !
    ! Scans through the list of time intervals adding the model predictions from mass map.
    ! Assumes that model valuse are picewize linear between timesteps
    ! Adds contribution from "now" state to relevant slots
    ! as a result all slots from now-timestep to now+timestep are affected
    ! In this formulation "now" can comtribute to motre than one observation slot
    ! Observations assuned to be non-overlapping in time (gaps possible)
    ! 
    ! Gotcha: every state of the massmap should be observed once and only once
    ! No account for incomplete averages made
    implicit none
    type(inSituObservation), intent(inout) :: obs
    type(Tmass_map), intent(in) :: map_c, map_cmX, map_cmY
    integer, intent(in) :: cloud_metric_flag, iSpOnes
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep
    type(TField_buffer), intent(in) :: disp_buf_ptr

    integer :: iTobs, iSeconds
    real, dimension(obs%nStations) :: val, mass2obs, delta2injmass
    real :: fObsEnd,fObsDur,fOverStart,fOverEnd,fL,fR, fWeight
    type(silja_time) :: time_start, time_end !Debug output
    character(len=*), parameter :: sub_name = 'observe_in_situ'

    if (.not. obs%inThisSubdomain) return ! No stations from obs to observe here

    if (obs%time_direction == backwards) then
      obs%observationTimestep = 1
      obs%time_direction = forwards
    end if

    iSeconds = nint(fu_sec(timestep))
    if (iSeconds < 0 ) then 
      call set_error("Observing on negative timesteps!", sub_name)
      return
    endif


    time_start = now-timestep !Sic! earliest time that can be affected by "now"
    do iTobs = obs%observationTimestep, obs%nTimes
        if (obs%endTimes(iTobs) >= time_start) exit
    enddo

    if (iTobs > obs%nTimes) return ! no more obs

    obs%observationTimestep = iTobs !save for future use


    ! Distribute it over slots
    val(1) = real_missing
    do iTobs = obs%observationTimestep, obs%nTimes
      fWeight = 0.
      if ( iSeconds == 0) then  !! instant massmap (e.g. 3Dvar) -> only instant obs. possible
        if (obs%endTimes(iTobs) ==  now) then
            if (.not. obs%durations(iTobs) == zero_interval) then
              call set_error("Non-instant observation from instant step", sub_name)
            endif
            fWeight = 1.

        endif
        if (obs%endTimes(iTobs) > now) exit
      else

       !normalize all times so that [now-timestep:now+timestep] maps to [-1 : 1] 
         fObsEnd = (fu_sec(obs%endTimes(iTobs) - now)) / iSeconds
         if (fObsEnd < -1.) call set_error("Must never happen",sub_name)

         fObsDur = fu_sec(obs%durations(iTobs)) / iSeconds
             !                        1_|                                     !
             !                          ^                                     !
             !                         /|\                                    !
             !                        / | \                                   !
             !                       /  | |\                                  !
             !                      /   | | \                                 !
             !                     /|   | |  \                                !
             !                    / |   | |   \                               !
             !                   /  |   | |    \                              !
             !                  /   |  0| |     \                             !
             !                ______________________                          !
             !                -1    L   0 R      1     -> Normalized time     !
             ! overlap of [-1 : 1] and normalized observation range           !
         fOverStart = max(-1.,fObsEnd - fObsDur) 
         if (fOverStart > 1.) exit
         fOverEnd   = min(fObsEnd, 1.)

         if (fObsDur < 1e-4) then ! obs shorter than .01% of time step
            fWeight = 1. - 0.5*abs(fOverStart+fOverEnd)
         else

             !Left trapezoid bounds
             fL = min(fOverStart, 0.)
             fR = min(fOverEnd, 0.)
             fWeight = ((fL+fR)*0.5+1.)*(fR-fL) !area
             
             ! Right trapezoid bounds
             fL = max(fOverStart, 0.)
             fR = max(fOverEnd, 0.)
             fWeight = fWeight + ((fL+fR)*0.5+1.)*(fR-fL) !area

             fWeight = fWeight / fObsDur ! Fraction of the integral
        endif !instant obs
      endif  !!!iSeconds == 0

      !Actual value from the massmap
      if (fWeight > 0) then
        if (val(1) == real_missing) then
             call   mass_to_obs_scalingAtStation(mass2obs, delta2injmass, obs, map_c, disp_buf_ptr, &
                                   & cloud_metric_flag, iSpOnes, now)
           val = fu_getValueAtStation(obs, map_c, map_cmX, map_cmY)
           val = val * mass2obs
        endif
!        if (obs%modeldata(iTobs) /= 0)then
!          call msg('obs%modelData( iTobs, obs%modeldata(iTobs,iTobs) /= 0:'))
!          call msg_warning('obs%modelData( 'observe_in_situ',iTobs) /= 0:')
!        endif
        obs%modelData(:,iTobs) = obs%modelData(:,iTobs) + fWeight*val(:)
        !! Here we assume that scaling also averages
        obs%inject_scaling(:,iTobs)   = obs%inject_scaling(:,iTobs) + fWeight * delta2injmass


#ifdef DEBUG_OBS
time_end = obs%endTimes(iTobs)
time_start = time_end - obs%durations(iTobs)
call msg('obs_start: ' // fu_str(time_start) // ', obs_end: ' // fu_str(time_end) // ', now: '// fu_str(now))
call msg('obs(1): ' // trim(obs%station(1)%id))
call msg('step, obs data tStep='//fu_str(iTobs), obs%obsData(iTobs,:))
call msg('modeldata', obs%modelData(iTobs,:))
call msg('val',  val)
call msg('Mean over concentration', sum(map_c%arM))
#endif
      endif

    enddo !iTobs

  end subroutine observe_in_situ

!******************************************************************************************

  subroutine observe_eruption(obs, eruption_height)
    implicit none
    type(t_eruptionObservation), intent(inout) :: obs
    real :: eruption_height
    obs%model_eruption_height(1) = eruption_height
  end subroutine observe_eruption

  !************************************************************************************
  ! get_data/set_data: collect the observed model values into arrays, or set them from an
  ! array (the latter is used by the line search algorithm).

  subroutine get_data_in_situ(obs, values_obs, values_mdl, variance)
    implicit none
    type(inSituObservation), intent(in) :: obs
    real, dimension(:), intent(out), optional :: values_obs, values_mdl, variance

    integer :: ii, iT, iSt

   !!   call msg("obs%dataLength "//trim(obs%tag), (/obs%nTimes, obs%nStations, obs%dataLength/))
    ii = 0
    do iT = 1, obs%nTimes
      do iSt = 1, obs%nStations
        if (obs%obsData(iSt,iT) == real_missing) cycle
        ii = ii + 1
        if (present(values_obs)) values_obs(ii) = obs%obsData(iSt,iT)
        if (present(values_mdl)) values_mdl(ii) = obs%modelData(iSt,iT)
        if (present(variance)) variance(ii) = obs%variance(iSt,iT)
      end do
    end do

    if (ii /= obs%dataLength) &
      & call set_error("Mismatch in valid data size", "get_data_in_situ")

  end subroutine get_data_in_situ

  !************************************************************************************

  subroutine set_data_in_situ(obs, values_mdl)
    !! Set obs data: to be done after sub_domain exchange
    implicit none
    type(inSituObservation), intent(inout) :: obs
    real, dimension(:), intent(in), optional :: values_mdl

    integer :: ii, iT, iSt

    ii = 0
    do iT = 1, obs%nTimes
      do iSt = 1, obs%nStations
        if (obs%obsData(iSt,iT) == real_missing) cycle
        ii = ii + 1
        obs%modelData(iSt,iT) = values_mdl(ii)
      end do
    end do

    if (ii /= obs%dataLength) &
      & call set_error("Mismatch in valid data size", "set_data_in_situ")

  end subroutine set_data_in_situ

  !************************************************************************************

  subroutine get_data_eruption(obs, values_obs, values_mdl, variance)
    implicit none
    type(t_eruptionObservation), intent(in) :: obs
    real, dimension(:), intent(out), optional :: values_obs, values_mdl, variance

    integer :: ind_param

    if (present(values_obs)) then
      values_obs(1) = obs%obs_eruption_height(1)
    end if
    if (present(values_mdl)) then
      values_mdl(1) = obs%model_eruption_height(1)
    end if
    if (present(variance)) then
      variance(1) = obs%variance(1)
    end if

  end subroutine get_data_eruption

  !************************************************************************************

  subroutine set_data_eruption(obs, values_mdl)
    implicit none
    type(t_eruptionObservation), intent(inout) :: obs
    real, dimension(:), intent(in) :: values_mdl

    integer :: ind_param

      obs%model_eruption_height(1) = values_mdl(1)

  end subroutine set_data_eruption

  !************************************************************************************

  subroutine get_localisation_in_situ(obs, lonlats)
    implicit none
    type(inSituObservation), intent(in) :: obs
    real, dimension(:,:), intent(out) :: lonlats

    integer :: ii, iSt, iT

    ii = 0
    do iT = 1, obs%nTimes
      do iSt = 1, obs%nStations
        if (obs%obsData(iSt,it) == real_missing) cycle
        ii = ii + 1
        lonlats(1,ii) = obs%station(iSt)%lon
        lonlats(2,ii) = obs%station(iSt)%lat
      end do
    end do
    if (ii /= obs%dataLength) &
      & call set_error("Mismatch in valid data size", "get_localisation_in_situ")

  end subroutine get_localisation_in_situ

  !************************************************************************************

  subroutine get_localisation_eruption(obs, lonlats)
    implicit none
    type(t_eruptionObservation), intent(in) :: obs
    real, dimension(:,:), intent(out) :: lonlats
    integer :: ii=1

    lonlats(1,ii) = obs%lon(1)
    lonlats(2,ii) = obs%lat(1)

  end subroutine get_localisation_eruption


  !************************************************************************************
  ! Restart: call before starting a forward run.
  ! or any other time if counters are to be reset - but then do not nullify model data
  !
  subroutine restart_in_situ(obs, ifNullifyModel)
    implicit none
    type(inSituObservation), intent(inout) :: obs
    logical, intent(in) :: ifNullifyModel

    obs%observationTimestep = 1
    obs%time_direction = forwards
    if(ifNullifyModel) then
      obs%modeldata = 0.0
      obs%inject_scaling = 0.0
    endif
  end subroutine restart_in_situ


  !************************************************************************************

  subroutine test_in_situ(grid, vert, species, species_transp, time_start, shift, timestep, obs_len, obs)
    implicit none
    type(inSituObservation), intent(out) :: obs
    type(silja_time), intent(in) :: time_start!, time_end
    type(silja_interval), intent(in) :: timestep, shift
    type(silja_grid), intent(in) :: grid
    type(silam_vertical), intent(in) :: vert
    type(silam_species), intent(in), target :: species ! observed
    type(silam_species), dimension(:), intent(in) :: species_transp
    integer, intent(in) :: obs_len

    !integer, parameter :: obs_len = 18

    integer :: ii, file_unit
    real, dimension(obs_len) :: obs_data, variance
    type(silja_time), dimension(obs_len) :: obs_times_start
    type(silja_interval), dimension(obs_len) :: obs_durations
    !type(silam_species) :: species
    type(observationStation) :: station
    type(silja_time) :: now, time_end
    type(tmass_map), pointer :: map_c, map_px, map_py, map_pz
    integer :: nx, ny
    real :: lon, lat
    integer, dimension(max_species) :: ind_obs2transp
    real, dimension(max_species) :: scale_transp2obs
    type(silam_material), pointer :: material
    type(TField_buffer), pointer :: disp_buf_ptr


    time_end = time_start + shift + one_hour*obs_len + shift

    ind_obs2transp(1) = fu_index(species,species_transp)
    material => fu_material(species)
    scale_transp2obs(1) = fu_conversion_factor(fu_basic_mass_unit(material), 'mole', material)

    !time_start = fu_set_time_utc(2011, 1, 1, 0, 0, 0.)
    !time_end = fu_set_time_utc(2011, 1, 2, 0, 0, 0.)
    !timestep = fu_set_interval_min(15)
    !shift = fu_set_interval_min(20)

    !grid = fu_set_grid('test_grid', lonlat, pole_geographical, 5.0, 35.0, 10, 10, 1.0, 1.0)
    !call set_vertical(level, vert)

    obs_times_start = (/(time_start + shift + one_hour*ii, ii=1, obs_len)/)
    obs_durations = one_hour
    obs_data = 1.0
    variance = 100.0

    !call set_species(species, fu_get_material_ptr('NO2'), in_gas_phase)

    call grid_dimensions(grid, nx, ny)
    lon = fu_lon_native_from_grid(real(nx/2), real(ny/2), grid)
    lat = fu_lat_native_from_grid(real(nx/2), real(ny/2), grid)
    station = fu_initObservationStation('S1', 'S1', lon, lat, 10.0)
    obs = fu_initInSituObservation1(obs_times_start, obs_durations, obs_data, variance, &
      & obs_len, constant_variance, .false., & !!ifmmr = .false.
                                 & ind_obs2transp,scale_transp2obs, 1, station, vert, 'test')
    if (error) return
    obs%obsData(1,:) = obs_data
    call set_mass_map(map_c,  concentration_flag, 3, 0, grid, vert, (/species/), val=2.0)
    call set_mass_map(map_px, advection_moment_x_flag, 3, 0, grid, vert, (/species/), val=0.0)
    call set_mass_map(map_py, advection_moment_y_flag, 3, 0, grid, vert, (/species/), val=0.0)
    call set_mass_map(map_pz, advection_moment_z_flag, 3, 0, grid, vert, (/species/), val=0.0)
    if (error) return
    now = time_start
    do while (now <= time_end)
      disp_buf_ptr => null()
      call observe_in_situ(obs, map_c, map_px, map_py, disp_buf_ptr, &
               & cloud_metric_geometry_flag, int_missing, now, timestep)
      if (error) return
      now = now + timestep
    end do

    file_unit = fu_next_free_unit()
    open(file_unit, file='test.observations.1')
    call obs_to_file_in_situ(obs, (/species/), file_unit)

    map_c%arm = 0.0
    do while (now >= time_start)
      print *, 'Sum of mass map:', sum(map_c%arm) / observation_scaling
      call inject_in_situ(obs, map_c, map_px, map_py, map_pz, now, fu_opposite(timestep), 0, 1)
      if (error) return
      now = now - timestep
    end do

  end subroutine test_in_situ



  !**********************************************************************************
  !
  ! Auxiliary routines
  !
  !**********************************************************************************

  subroutine destroyInSitu(obs)
    implicit none
    type(inSituObservation) :: obs
    integer :: iStat
    
    deallocate(obs%durations, obs%obsData, obs%modelData, obs%endTimes, &
         obs%variance, obs%ind_obs2transp, obs%scale_transp2obs, stat=iStat)

    if (iStat /= 0) call set_error('Deallocation error','destroyObservationInSitu')

  end subroutine destroyInSitu

  !**********************************************************************************

  subroutine destroyEruption(obs)
    implicit none
    type(t_eruptionObservation) :: obs

    integer :: iStat

    deallocate(obs%duration, obs%obs_eruption_height, obs%model_eruption_height, obs%endTime, &
         & obs%lat, obs%lon, obs%variance, stat=iStat)
    if (iStat /= 0) then
      call set_error('Deallocation error','destroyEruption')
    end if
  end subroutine destroyEruption

  !**********************************************************************************

  function fu_getValueAtStation(obs, mapConc, map_cmX, map_cmY) result(v)
    !gets a value corresponding to the  observation from the massmap

    implicit none
    type(InSituObservation), intent(in) :: obs
    type(Tmass_map), intent(in)  :: mapConc, map_cmX, map_cmY ! concentrations and mass centres
    real, dimension(obs%nStations) :: v  !! Vector of values

    integer :: ispecies_transp, ispecies_obs, iTmp
    real :: fract, weight, weight_past, cmX, cmY, fStdX, fStdY
    real, parameter :: fZcTrimax = 1.0/6.0  ! Maximum CM for trapezoid slab

    real :: air_mass, air_dens, ones_mass
    logical :: ifNeedDens, ifNeedCellmass 
    integer :: iSt !! Station index of the observation


    integer :: ind_air_mass, ind_air_dens
    character(len=*), parameter :: sub_name = 'fu_getValueAtStation'

    v(:) = 0.  !! Should stay zero if the station fully out of our domain
    ! Whole domain is obtained by smpi_allreduce_add
    ! Warning! Stations not filled by anyone will be silently reported as zeros
    !
    ! There can be several ways to deal with data extraction.
    ! For the time being, we allow either grid-cell average (a.k.a. nearest_neighbour)
    ! or interpolated to the mass centre point. In the latter case, we limit the variation
    ! of the mass-centre to the extent when trapezoid turns into triangle: empty space
    ! in the grid cell is not allowed. Hence, -1/6 < cm < 1/6
    !
    ! For the stations from outside the domain coordinates set 
    ! Check moved one level up
    !if (any((/obs%station%ix_dispersion, obs%station%iy_dispersion/) == int_missing)) return

    select case(obs%iInterpolation)
    case(nearest_point)
      do iSt = 1, obs%nStations
        if (.not. obs%station(iSt)%inThisSubdomain) cycle 
        do ispecies_obs = 1, obs%num_obs_species
          ispecies_transp = obs%ind_obs2transp(ispecies_obs)
          fract = obs%scale_transp2obs(ispecies_obs)
          
          !!!! 
          v(iSt) = v(iSt) + mapconc%arm(ispecies_transp, 1, obs%iLevel, obs%station(iSt)%ix_dispersion, &
                 & obs%station(iSt)%iy_dispersion)  * fract 

        enddo !iSp
      enddo ! iSt
      
    case(linear)
      call set_error('Linear interpolation between the grid cells is considered wrong',sub_name)
      return
      
    case(toMassCentreLinear)
      do iSt = 1, obs%nStations
        if (.not. obs%station(iSt)%inThisSubdomain) cycle 
        do ispecies_obs = 1, obs%num_obs_species
          ispecies_transp = obs%ind_obs2transp(ispecies_obs)
          fract = obs%scale_transp2obs(ispecies_obs)
          cmX = max(-fZcTrimax, min(fZcTrimax, map_cmX%arm(ispecies_transp, 1, obs%ilevel, &
               & obs%station(iSt)%ix_dispersion, &
               & obs%station(iSt)%iy_dispersion)))
          cmY = max(-fZcTrimax, min(fZcTrimax, map_cmY%arm(ispecies_transp, 1, obs%ilevel, &
               & obs%station(iSt)%ix_dispersion, &
               & obs%station(iSt)%iy_dispersion)))
          fStdX = obs%station(iSt)%xLoc_dispGrid
          fStdY = obs%station(iSt)%yLoc_dispGrid
          ! Slope on each dimension is: 12 * cm                                                                                                  
          weight = (1. + 12. * cmX * fStdX) * (1. + 12. * cmY * fStdY)
          v(iSt) = v(iSt) + mapconc%arm(ispecies_transp, 1, obs%iLevel, obs%station(iSt)%ix_dispersion, &
               & obs%station(iSt)%iy_dispersion) * fract * weight
        enddo !iSp
      enddo ! iSt

    case default
      call set_error('Unknown interpolation type:' + fu_str(toMassCentreLinear),sub_name)
      return
    end select

  end function fu_getValueAtStation

  !**********************************************************************************

  subroutine mass_to_obs_scalingAtStation(mass2obs, delta2injmass, obs, mapConc, buf, &
                  & cloud_metric_flag, iSpOnes, now) 
    !gets a value corresponding to the  observation from the massmap

    implicit none
    type(InSituObservation), intent(in) :: obs
    real, dimension(obs%nStations), intent(out) :: mass2obs, & !! Scale massmap quantity to observed units
                     & delta2injmass !! scale injected discrepancy  (obs - mod)/variance according 
                                     !! to the used massmap metric
    type(Tmass_map), intent(in)  :: mapConc ! concentrations, inly ones might be needed
    integer, intent(in) :: cloud_metric_flag, iSpOnes
    type(TField_buffer), intent(in) :: buf !Dispersion buffer
    type(silja_time), intent(in) :: now

    integer :: ispecies_transp, ispecies_obs, iDisp, iSt
    real :: weight_past

    real :: air_mass, air_dens, ones_mass
    logical :: ifNeedDens, ifNeedCellmass 

    integer :: ind_air_mass, ind_air_dens
    character(len=*), parameter :: sub_name = 'mass_to_obs_scalingAtStation'

    ifNeedDens = (obs%ifmmr .eqv. (cloud_metric_flag == cloud_metric_geometry_flag))
    ifNeedCellmass = (cloud_metric_flag == cloud_metric_cellmass_flag)

    if (ifNeedCellmass .or. ifNeedDens) then
      !! weight_past is normally for the mid-timestep here we are at the 'moment' 
      weight_past = (buf%time_future - now) / (buf%time_future - buf%time_past)
    endif

    if (ifNeedDens) then
      !! air_dens needed below
      ind_air_dens = fu_index(buf, air_density_flag)
      if (ind_air_dens < 1) then
        call set_error('Failed to find the air density field', sub_name)
        return
      end if
    endif

    if (ifNeedCellmass) then
       ind_air_mass = fu_index(buf, disp_cell_airmass_flag)
       if (ind_air_mass < 1) then                              
         call set_error('Failed to find the air mass field', sub_name)
         return
       end if
    endif
    
    do iSt = 1, obs%nStations
      if (.not. obs%station(iSt)%inthissubdomain) then
          mass2obs(iSt) = real_missing !! Just not to leave uninitialized
          delta2injmass(iSt) = real_missing
          cycle !! cannot do anything beyond our domain
      endif

      if (ifNeedCellmass .or. ifNeedDens) then 
        iDisp = obs%station(iSt)%ix_dispersion+(obs%station(iSt)%iy_dispersion-1)*mapConc%nx !!Dispersion index
      endif
      if (ifNeedDens) then
        air_dens = buf%p4d(ind_air_dens)%past%p2d(obs%iLevel)%ptr(iDisp) * weight_past + &
                 & buf%p4d(ind_air_dens)%future%p2d(obs%iLevel)%ptr(iDisp) * (1. - weight_past)
      endif


      select case (cloud_metric_flag)
        case (cloud_metric_geometry_flag)
           if (obs%ifmmr) then
             mass2obs(iSt) = 1. / (obs%cell_volume(iSt) * air_dens)
             delta2injmass(iSt) =  1. / air_dens !! vmr in geometrical metric
           else
             mass2obs(iSt) = 1./ obs%cell_volume(iSt)
             delta2injmass(iSt) = 1. !! concentration in gemetrical metric
           endif

         case (cloud_metric_cellmass_flag)
           air_mass = buf%p4d(ind_air_mass)%past%p2d(obs%iLevel)%ptr(iDisp) * weight_past + &
                    & buf%p4d(ind_air_mass)%future%p2d(obs%iLevel)%ptr(iDisp) * (1. - weight_past)

           if (obs%ifmmr) then
             mass2obs(iSt) = 1. / air_mass
             delta2injmass(iSt) = 1.
           else
             mass2obs(iSt) = air_dens / air_mass !! Same inverse cell volume, but without geometry.
             delta2injmass(iSt) =  air_dens
           endif

         case(cloud_metric_ones_flag)
           ones_mass =  mapconc%arm(iSpOnes, 1, obs%iLevel, obs%station(iSt)%ix_dispersion, obs%station(iSt)%iy_dispersion)
           if (obs%ifmmr) then
             mass2obs(iSt) = 1. / ones_mass  !! As above, but ones_mass instead of air_mass
             delta2injmass(iSt) = 1.
           else
             mass2obs(iSt) = air_dens / ones_mass   
             delta2injmass(iSt) =  air_dens
           endif

        case default
          call set_error("Strange cloud_metric_flag: "//trim(fu_str(cloud_metric_flag)), sub_name)
      end select 
    enddo   !! Stations
  end subroutine mass_to_obs_scalingAtStation
  

  !**********************************************************************************
  !
  ! Accessor functions
  !
  !**********************************************************************************
  function fu_station_label(stat)
    implicit none
    type(observationStation), intent(in) :: stat
    character(len=STATION_LABEL_LENGTH) :: fu_station_label
    fu_station_label = stat%label
  end function fu_station_label

  integer function fu_size_in_situ(obs) result(n)
    implicit none
    type(inSituObservation), intent(in) :: obs
    n = obs%dataLength
  end function fu_size_in_situ

  integer function fu_size_eruption(obs) result(n)
    implicit none
    type(t_eruptionObservation), intent(in) :: obs
    n = 1
  end function fu_size_eruption

  real function fu_lat(stat)
    implicit none
    type(observationStation), intent(in) :: stat
    fu_lat = stat%lat
  end function fu_lat

  real function fu_lon(stat)
    implicit none
    type(observationStation), intent(in) :: stat
    fu_lon = stat%lon
  end function fu_lon

  elemental logical function fu_station_defined(station)
    implicit none
    type(observationStation), intent(in) :: station
    fu_station_defined = station%defined
  end function fu_station_defined

  subroutine report_obs(obs)
     type(inSituObservation) :: obs
     integer :: i
     call msg("")
     call msg("Observation: "+obs%tag)

     call msg("obs%dataLength", obs%dataLength)
     call msg("obs%endtimes:"+fu_str(obs%endTimes(1)) + "_to_" + fu_str(obs%endTimes(obs%dataLength) ))
     call msg("obs%durations:"+fu_str(obs%durations(1)) + "_to_" + fu_str(obs%durations(obs%dataLength) ))

  end subroutine report_obs

end module Observations_in_situ
