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

  private

  public inject_in_situ
  public observe_in_situ
  public observe_eruption

  public fu_initInSituObservation
  public fu_initObservationStation
  public fu_initEruptionObservation
  public fu_label
  public fu_lat
  public fu_lon

  !public destroyObservation
  public destroy
  public defined

  public observationStation
  public inSituObservation
  public t_eruptionObservation
  public observationStationPtr
  public t_column_observation
  public test_in_situ
  public obs_to_file
  public fu_size
  public get_data
  public set_data
  public report
  public get_localisation
  public restart_in_situ
  public restart_column

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

  interface obs_to_file
     module procedure obs_to_file_in_situ
  end interface

  interface set_data
     module procedure set_data_in_situ
     module procedure set_data_column
     module procedure set_data_eruption
  end interface

  interface get_data
     module procedure get_data_in_situ
     module procedure get_data_column
     module procedure get_data_eruption
  end interface


  interface get_localisation
     module procedure get_localisation_in_situ
     module procedure get_localisation_eruption
  end interface

  interface report
     module procedure report_obs
  end interface
  integer, parameter, public :: STATION_ID_LENGTH = 16
  integer, parameter, public :: STATION_LABEL_LENGTH = 64

  integer, parameter, public :: constant_variance = 1, variable_variance = 2

  real, public, parameter :: observation_scaling = 1e5


  !**********************************************************************************
  !
  ! Type definition for an in situ observation
  !
  !**********************************************************************************
  ! This structure defines a time series type observation from a single station. Each
  ! station has a fixed postion and measures a single substance only.

  type inSituObservation
     ! something to identify this observation in the dump file.
     character(len=substNmLen) :: tag = '???'
     integer :: iInterpolation = linear    ! nearest_point
     ! for each measurement:
     type(silja_time), dimension(:), pointer :: endTimes   !startingTimes ! => null()
     type(silja_interval), dimension(:), pointer :: durations
     ! the real observed data, and that "observed" from the model
     real, dimension(:), pointer :: obsData, modelData, variance
     type(observationStation), pointer :: station => null()
     !
     ! Level. Currently not used.
     type(silja_level) :: level = ground_level
     integer :: ilevel = 1

     ! The mapping from (possibly several) transport to the observed species. Also
     ! contains factors for unit conversion.
     integer, dimension(:), pointer :: ind_obs2transp
     real, dimension(:), pointer :: scale_transp2obs
     integer :: num_obs_species

     ! the length of data and modelData
     integer :: dataLength

     ! the counter. As the model proceeds, and observer/inject routines are
     ! called, the observation keeps a count of the index of current measurement.
     ! In forward run, it runs from 1 to dataLength, and the opposite in the adjoint
     ! case.
     integer :: observationTimestep = 1
     integer :: time_direction

     ! The statistical weight for measured zero concentrations
     ! and others. This will most likely be the same for every observation,
     ! but in principle, each observation should know its own uncertainty.
     ! (the weights are equal to the 1/variance of the observation)
     ! real :: weightZero = 1 , weightNonzero = 1
     ! the measured concentration is considered zero if under
     real :: measurementEpsilon = 1.0e-19

     ! Volume of the grid cell, for conversion from mass to concentration.
     !
     real :: cell_volume
  end type inSituObservation

  type t_eruptionObservation
    character(len=substNmLen) :: tag = '???'
    type(silja_time), dimension(:), pointer :: endTime   !startingTimes ! => null()
    type(silja_interval), dimension(:), pointer :: duration
    real, dimension(:), pointer :: model_eruption_height, obs_eruption_height, variance
    real, dimension(:), pointer :: lat, lon
  end type t_eruptionObservation

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

  type observationStationPtr
     type(observationStation), pointer :: ptr
  end type observationStationPtr

  type(observationStation), parameter, public ::  station_missing = &
        & observationStation('', '', real_missing, real_missing, &
                            & int_missing, int_missing, int_missing, int_missing, &
!                            & (/int_missing, int_missing/), (/int_missing, int_missing/), &
!                            & (/int_missing, int_missing/), (/int_missing, int_missing/), &
                            & real_missing, real_missing, .false.,real_missing, int_missing, &
                            & .false.)

  type t_column_observation
     ! A column based observation. In contrast to in-situ, the column observations are
     ! taken to be instant, and no time integration/averaging is done.
     real, dimension(:), pointer :: fx, fy
     real, dimension(:), pointer :: value, modelvalue, variance, cell_area
     real, dimension(:,:,:), pointer :: kernel_val
     integer, dimension(:), pointer :: kernel_vert_ind, kernel_chem_ind, kernel_opt_ind
     type(silja_time), dimension(:), pointer :: timelist
     logical :: has_variable_kernel
     logical :: defined = .false.
     logical :: is_aod = .false.
     integer :: n_values, n_obs_species, n_obs_levels
     integer :: ind_step
     integer :: time_direction
     character(len=STATION_LABEL_LENGTH) :: label
  end type t_column_observation

  !***********************************************************************************************

contains

  ! Initialization. Returns a new observation object with given data,
  ! locations, and time windows. The returned object is ready for use
  ! at any time.  Array pointers endTimes and durations should be
  ! allocated and dipersionGrid should be defined in advance.


  !***********************************************************************************************


  function fu_initInSituObservation(endTimes, durations, obsData, variance, &
                                  & n_values, variance_flag, &
                                  & ind_obs2transp, scale_transp2obs, nSpecies, &
                                  & station, &
                                  & level, vertical, tag) &
                                  & result(newObservation)
    implicit none
    ! input
    type(silja_time), dimension(:), intent(in)  :: endTimes
    type(silja_interval), dimension(:), intent(in) :: durations
    real, dimension(:), intent(in) :: obsData, variance
    integer, dimension(:), intent(in) :: ind_obs2transp
    real, dimension(:), intent(in) :: scale_transp2obs
    integer, intent(in) :: variance_flag, n_values, nSpecies
    ! the durations are currently assumed to be given in hours. True,
    ! silja_interval would make sense in this context.
    type(observationStation), intent(in), target:: station
    type(silam_vertical), intent(in) :: vertical
    type(silja_level), intent(in) :: level
    character(len=*), intent(in) :: tag

    type(inSituObservation) :: newObservation

    ! local
    integer :: status, isp_obs

    allocate(newobservation%endTimes(n_values), newObservation%durations(n_values), &
           & newObservation%obsData(n_values), newObservation%variance(n_values), &
           & newObservation%modelData(n_values), &
           & newObservation%ind_obs2transp(nSpecies), &
           & newObservation%scale_transp2obs(nSpecies), stat=status)
    if(fu_fails(status == 0, 'Allocate failed', 'fu_initInSituObservation1'))return

    newObservation%tag = tag
    newObservation%endTimes = endTimes(1:n_values)
    newObservation%durations = durations(1:n_values)
    newObservation%obsData = obsData(1:n_values)
    newObservation%dataLength = n_values
    newObservation%station => station

    newObservation%modelData = 0.0

    select case(variance_flag)
    case (variable_variance)
      newObservation%variance = variance(1:n_values)
    case (constant_variance)
      newObservation%variance(:) = variance(1)
    case default
      call set_error('Strange variance_flag', 'fu_initInSituObservation1')
      return
    end select


    newObservation%num_obs_species = nSpecies
    newObservation%ind_obs2transp(1:nSpecies) = ind_obs2transp(1:nSpecies)
    newObservation%scale_transp2obs(1:nSpecies) = scale_transp2obs(1:nSpecies)

    ! We assume the measurement is on the first model layer. These go
    ! with the defaults.
    !
    if (defined(level)) then
      call set_error('No support for defined levels yet', 'fu_initInSituObservation1')
      return
    end if

    newObservation%cell_volume = newObservation%station%cell_area * &
                                            & fu_layer_thickness_m(fu_level(vertical, 1))

    newObservation%obsData = obsData(1:n_values)

  end function fu_initInSituObservation

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

  function fu_initObservationStation(ID, label, lon, lat, z, grid) result(new_station)

    ! Create an observation station. Dispersion grid is needed for computing the
    ! projection from the geo grid. No interpolation is implemented due to
    ! trouble with domain split.
    ! Since grids in SILAM are inaccurate we project always to wholeMPIdispersion_grid
    ! and then subtract offsets
    implicit none
    character(len=*), intent(in) :: id, label
    real, intent(in) :: lat, lon, z ! z doesn't do anything yet
    type(silja_grid), intent(in) :: grid

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

#ifdef DEBUG
      write(str, fmt='("OUT ", A10, 1x, F8.2, F8.2," ", A20, I6, I6)') id, lat, lon, trim(label), ix, iy
      call msg(str)
#endif
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
                                            & real_missing, 1, &! cell_area, ind_lev
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

    integer :: ii
    character (len=worksize_string) :: sp

    sp = fu_str(species_trn(obs%ind_obs2transp(1)))
    do ii = 2, size(obs%ind_obs2transp)
      sp = sp + "__" + fu_str(species_trn(obs%ind_obs2transp(ii)))
    end do

    do ii = 1, obs%datalength
      write(file_unit, fmt='(2(A,1x), I4, 2I3, 1x, 2(F5.2,1x), 3G12.3)') &
                  & trim(obs%station%id), trim(sp), &
                  & fu_year(obs%endTimes(ii)), fu_mon(obs%endTimes(ii)), fu_day(obs%endTimes(ii)), &
                  & fu_hour(obs%endTimes(ii)) + fu_min(obs%endTimes(ii)) / 60. + &
                                              & fu_sec(obs%endTimes(ii)) / 3600., &
                  & fu_hours(obs%durations(ii)), &
                  & obs%obsData(ii), obs%modeldata(ii), obs%variance(ii)
    end do
  end subroutine obs_to_file_in_situ



  !**********************************************************************************
  !
  ! Inject and observe
  !
  !**********************************************************************************

  subroutine inject_in_situ(obs, map_c, map_px, map_py, map_pz, now_, timestep, iThread, nThreads, iInjectionPlace)
    implicit none
    type(inSituObservation), intent(inout) :: obs
    type(Tmass_map), intent(inout) :: map_c, map_px, map_py, map_pz !!!!Centers of mass here
    type(silja_time), intent(in) :: now_
    type(silja_interval), intent(in) :: timestep
    integer, intent(in) :: iThread, nThreads, iInjectionPlace

    integer :: ind_obs, ind_src, ix_obs, iy_obs, ilev_obs, isp_obs
    real :: overlap_sec, obs_val, mdl_val, weight, inject_mass, mass_fract
    real :: cmX, cmY !! Station location with respect to the cell-center (+- 0.5)
    real ::  mOld, mAdd ! mass in the map prior to inhection and it increment
    integer :: isp_transp
    type(silja_interval) :: overlap
    type(silja_time) :: time_start, time_end, step_start, now

    now = now_  !- timestep

    if (.not. obs%station%inThisSubdomain) return ! Station from foreign subdomain

    ix_obs = obs%station%ix_dispersion
    if (mod(ix_obs, nThreads) /= iThread) return !!!Make it possible for parallel injection
    iy_obs = obs%station%iy_dispersion
    cmX = obs%station%xLoc_dispGrid     ! injection goes into the centre of mass == location of station
    cmY = obs%station%yLoc_dispGrid
    ilev_obs = obs%station%ind_lev

    if (obs%time_direction == forwards) then
      obs%observationTimestep = obs%datalength
      obs%time_direction = backwards
    end if

    step_start = now + timestep
    ind_obs = obs%observationTimestep

    time_end = obs%endTimes(ind_obs)
    time_start = time_end - obs%durations(ind_obs)
    do while (time_start > now .and. ind_obs > 1)
      call advance()
    end do

    if (time_start >= now) return ! used all observations

    do while (time_end > step_start)
      overlap = fu_time_overlap(time_start, time_end, step_start, now)
      obs_val = obs%obsData(ind_obs)
      mdl_val = obs%modelData(ind_obs)
      if (obs_val .eps. 0.0) then
        ind_src = DA_ZERO
        weight = 1.0
      else if (obs_val > mdl_val) then
        ind_src = DA_NEGATIVE
        weight = DA_NEGT_COEF_OBS
      else
        ind_src = DA_POSITIVE
        weight = 1.0  ! obs%weightNonZero
      end if
      !call msg('obs_val and mdl_val', obs_val, mdl_val)
      if (fu_fails(.not.(overlap == zero_interval), 'Impossible error', 'inject_in_situ')) return

      inject_mass = fu_sec(overlap) / fu_sec(obs%durations(ind_obs)) * (mdl_val - obs_val) &
           & * observation_scaling !* obs%cell_volume
      inject_mass = inject_mass / obs%variance(ind_obs)
      weight = weight
      !call msg('weight, var', weight, obs%variance(ind_obs))
      do isp_obs = 1, obs%num_obs_species
        isp_transp = obs%ind_obs2transp(isp_obs)
        mass_fract = obs%scale_transp2obs(isp_obs)
        if (mass_fract > 0.0) then

          !call msg('Injecting mass in:' + fu_str(ix_obs) + ',' &
          !       & + fu_str(iy_obs), mass_fract*weight*inject_mass)

          map_c%ifColumnValid(ind_src, ix_obs, iy_obs) = .true.
          map_c%ifGridValid(1, ind_src) = .true.
          mOld = map_c%arm(isp_transp,ind_src,ilev_obs,ix_obs, iy_obs) !What is there
          mAdd = mass_fract*weight*inject_mass !increment
          map_c%arm(isp_transp,ind_src,ilev_obs,ix_obs, iy_obs)  = mOld + mAdd

          select case(iInjectionPlace)
            case(nearest_point)   ! means cell-centre in the current context => no change to moment
            case(toMassCentreLinear)    ! use the station location to modify the moment
              map_px%arm(isp_transp,ind_src,ilev_obs,ix_obs,iy_obs) = &
                          & (map_px%arm(isp_transp,ind_src,ilev_obs,ix_obs,iy_obs) * mOld + &
                                                                  & mAdd * cmX) / (mOld + mAdd)
              map_py%arm(isp_transp,ind_src,ilev_obs,ix_obs,iy_obs) = &
                          & (map_py%arm(isp_transp,ind_src,ilev_obs,ix_obs,iy_obs) * mOld + &
                                                                  & mAdd * cmY) / (mOld + mAdd)
!!$          map_pz%arm(isp_transp,ind_src,ilev_obs,ix_obs,iy_obs) &
!!$               & = map_pz%arm(isp_transp,ind_src,ilev_obs,ix_obs,iy_obs) + mass_fract*weight*inject_mass*crdz
!              if (abs(map_px%arm(isp_transp,ind_src,ilev_obs,ix_obs,iy_obs)) > 0.5 .or.  &
!                & abs(map_py%arm(isp_transp,ind_src,ilev_obs,ix_obs,iy_obs)) > 0.5) then
!                call msg('Old mass and added mass:', mOld, mAdd)
!                call msg('Site coordinates:', ix_obs, iy_obs)
!                call msg('New X_cm and Y_cm:', map_px%arm(isp_transp,ind_src,ilev_obs,ix_obs,iy_obs), &
!                                             & map_py%arm(isp_transp,ind_src,ilev_obs,ix_obs,iy_obs))
!                call set_error ("Trouble injecting observation", "inject_in_situ")
!              end if
            case default
              call set_error('Unknown interpolatino type:' + fu_str(iInjectionPlace),'inject_in_situ')
              return
          end select  ! iInjectionPlace
        end if   ! mass_fract > 0
      end do   ! obs_species
      if (time_start < step_start) then
        exit
      else if (ind_obs == 1) then
        exit
      else
        call advance()
      end if

    end do ! while over time
!call msg("Done injection")

    obs%observationTimestep = ind_obs

  contains

    subroutine advance()
      implicit none
#ifdef DEBUG_MORE
      call msg('Advancing from:' + fu_str(time_end) + ', towards:' + fu_str(time_start) + ', ind_obs_times=', ind_obs)
#endif
      ind_obs = ind_obs - 1
      time_end = obs%endTimes(ind_obs)
      time_start = time_end - obs%durations(ind_obs)
    end subroutine advance

  end subroutine inject_in_situ

  !************************************************************************

  subroutine observe_in_situ(obs, map_c, map_cmX, map_cmY, iInterpolation, now, timestep)
    !
    ! Scans through the lsit of stations ordering the model predictions from mass map
    !
    implicit none
    type(inSituObservation), intent(inout) :: obs
    type(Tmass_map), intent(in) :: map_c, map_cmX, map_cmY
    integer, intent(in) :: iInterpolation
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep

    integer :: ind_obs
    real :: seconds, val, time_fract
    type(silja_time) :: time_start, time_end, step_end
    type(silja_interval) :: overlap

    seconds = fu_sec(timestep)

    if (obs%time_direction == backwards) then
      obs%observationTimestep = 1
      obs%time_direction = forwards
    end if

    ind_obs = obs%observationTimestep

    time_end = obs%endTimes(ind_obs)
    time_start = time_end - obs%durations(ind_obs)

    !call msg('ind_obs', ind_obs, obs%datalength)
    !call msg(fu_str(time_start))
    !call msg(fu_str(time_end))

    do while (time_end <= now .and. ind_obs < obs%datalength)
      call advance()
    end do

    if (time_end <= now) then
      !call msg('obs return')
      ! used all observations
      return
    end if

    if (obs%cell_volume > 0.) then !Station is in our domain and has a valid cell volume
      step_end = now + timestep
      do while (time_start < step_end)
#ifdef DEBUG_MORE
call msg('obs_start:' + fu_str(time_start) + ', obs_end: ' + fu_str(time_end) + ', now:' + fu_str(now))
#endif
        ! observe
        overlap = fu_time_overlap(time_start, time_end, now, step_end)
        if (fu_fails(.not.(overlap == zero_interval), 'Impossible error', 'observe_in_situ')) return

        time_fract = fu_sec(overlap) / fu_sec(obs%durations(ind_obs))
        val = fu_getValueAtStation(obs, map_c, map_cmX, map_cmY, iInterpolation) / obs%cell_volume
        obs%modeldata(ind_obs) = obs%modeldata(ind_obs) + time_fract * val
#ifdef DEBUG_MORE
call msg('obs:' + obs%station%id)
call msg('step, obs data', ind_obs, obs%obsData(ind_obs))
call msg('modeldata, new val', obs%modeldata(ind_obs), val)
#endif
        if (time_end > step_end) then
          ! timestep ends before observation -> leave the remainder for the next time step.
          exit
        else if (ind_obs < obs%datalength) then
          call advance()
        else
          ! all obs times processed
          exit
        end if
      end do
    else
      ! No observations from our domain
     !!!!obs%modeldata(ind_obs) = real_missing ! incompatible with reduce_add
      obs%modeldata(ind_obs) = 0.
    endif

    obs%observationTimestep = ind_obs
    !call msg('timestep set to', ind_obs)
  contains

    subroutine advance()
      implicit none
      ind_obs = ind_obs + 1
      time_end = obs%endTimes(ind_obs)
      time_start = time_end - obs%durations(ind_obs)
      obs%modeldata(ind_obs) = 0.0
    end subroutine advance

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

    integer :: ii

    if (present(values_obs)) then
      do ii = 1, obs%datalength
        values_obs(ii) = obs%obsData(ii)
      end do
    end if

    if (present(values_mdl)) then
      do ii = 1, obs%datalength
        values_mdl(ii) = obs%modeldata(ii)
      end do
    end if

    if (present(variance)) then
      do ii = 1, obs%datalength
        variance(ii) = obs%variance(ii)
      end do
    end if

  end subroutine get_data_in_situ

  subroutine get_data_column(obs, values_obs, values_mdl, variance)
    implicit none
    type(t_column_observation), intent(in) :: obs
    real, dimension(:), intent(out), optional :: values_obs, values_mdl, variance

    integer :: ii

    if (present(values_obs)) then
      do ii = 1, obs%n_values
        values_obs(ii) = obs%value(ii)
      end do
    end if

    if (present(values_mdl)) then
      do ii = 1, obs%n_values
        values_mdl(ii) = obs%modelvalue(ii)
      end do
    end if

    if (present(variance)) then
      do ii = 1, obs%n_values
        variance(ii) = obs%variance(ii)
      end do
    end if

    if (present(values_mdl) .and. present(values_obs)) then
      do ii = 1, obs%n_values
         !call msg('obs', values_obs(ii))
         !call msg('mdl', values_mdl(ii))
      end do
    end if
  end subroutine get_data_column

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

  subroutine set_data_eruption(obs, values_obs, values_mdl, variance)
    implicit none
    type(t_eruptionObservation), intent(inout) :: obs
    real, dimension(:), intent(in), optional :: values_obs, values_mdl, variance

    integer :: ind_param

    if (present(values_obs)) then
      obs%obs_eruption_height(1) = values_obs(1)
    end if
    if (present(values_mdl)) then
      obs%model_eruption_height(1) = values_mdl(1)
    end if
    if (present(variance)) then
      obs%variance(1) = variance(1)
    end if

  end subroutine set_data_eruption

  subroutine set_data_in_situ(obs, values_mdl)
    implicit none
    type(inSituObservation), intent(inout) :: obs
    real, dimension(:), intent(in) :: values_mdl

    integer :: ii

    do ii = 1, obs%datalength
      obs%modeldata(ii) = values_mdl(ii)
    end do

  end subroutine set_data_in_situ

  subroutine set_data_column(obs, values_mdl)
    implicit none
    type(t_column_observation), intent(inout) :: obs
    real, dimension(:), intent(in) :: values_mdl

    integer :: ii

    do ii = 1, obs%n_values
      obs%modelvalue(ii) = values_mdl(ii)
    end do

  end subroutine set_data_column

  subroutine get_localisation_in_situ(obs, lonlats)
    implicit none
    type(inSituObservation), intent(in) :: obs
    real, dimension(:,:), intent(out) :: lonlats

    integer :: ii

    do ii = 1, obs%datalength
      lonlats(1,ii) = obs%station%lon
      lonlats(2,ii) = obs%station%lat
    end do

  end subroutine get_localisation_in_situ

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
    if(ifNullifyModel) obs%modeldata = 0.0
  end subroutine restart_in_situ

  subroutine restart_column(obs, ifNullifyModel)
    implicit none
    type(t_column_observation), intent(inout) :: obs
    logical, intent(in) :: ifNullifyModel

    obs%ind_step = 1
    obs%time_direction = forwards
    if(ifNullifyModel) obs%modelvalue = 0.0
  end subroutine restart_column

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
    type(silja_level) :: level
    !type(silam_species) :: species
    type(observationStation) :: station
    type(silja_time) :: now, time_end
    type(tmass_map), pointer :: map_c, map_px, map_py, map_pz
    integer :: nx, ny
    real :: lon, lat
    integer, dimension(max_species) :: ind_obs2transp
    real, dimension(max_species) :: scale_transp2obs
    type(silam_material), pointer :: material



    time_end = time_start + shift + one_hour*obs_len + shift

    ind_obs2transp(1) = fu_index(species,species_transp)
    material => fu_material(species)
    scale_transp2obs(1) = fu_conversion_factor(fu_basic_mass_unit(material), 'mole', material)

    !time_start = fu_set_time_utc(2011, 1, 1, 0, 0, 0.)
    !time_end = fu_set_time_utc(2011, 1, 2, 0, 0, 0.)
    !timestep = fu_set_interval_min(15)
    !shift = fu_set_interval_min(20)

    !grid = fu_set_grid('test_grid', lonlat, pole_geographical, 5.0, 35.0, 10, 10, 1.0, 1.0)
    !level = fu_set_level(layer_btw_2_height, fval1=0.0, fval2=100.0)
    !call set_vertical(level, vert)

    obs_times_start = (/(time_start + shift + one_hour*ii, ii=1, obs_len)/)
    obs_durations = one_hour
    obs_data = 1.0
    variance = 100.0

    !call set_species(species, fu_get_material_ptr('NO2'), in_gas_phase)

    call grid_dimensions(grid, nx, ny)
    lon = fu_lon_native_from_grid(real(nx/2), real(ny/2), grid)
    lat = fu_lat_native_from_grid(real(nx/2), real(ny/2), grid)
    station = fu_initObservationStation('S1', 'S1', lon, lat, 10.0, grid)
    obs = fu_initInSituObservation(obs_times_start, obs_durations, obs_data, variance, &
                                 & obs_len, constant_variance, &
                                 & ind_obs2transp,scale_transp2obs,1, station, level_missing, vert, 'test')
    if (error) return
    obs%obsData = obs_data
    call set_mass_map(map_c,  concentration_flag, 3, 0, grid, vert, (/species/), val=2.0)
    call set_mass_map(map_px, advection_moment_x_flag, 3, 0, grid, vert, (/species/), val=0.0)
    call set_mass_map(map_py, advection_moment_y_flag, 3, 0, grid, vert, (/species/), val=0.0)
    call set_mass_map(map_pz, advection_moment_z_flag, 3, 0, grid, vert, (/species/), val=0.0)
    if (error) return
    now = time_start
    do while (now <= time_end)
      call observe_in_situ(obs, map_c, map_px, map_py, nearest_point, now, timestep)
      if (error) return
      now = now + timestep
    end do

    file_unit = fu_next_free_unit()
    open(file_unit, file='test.observations.1')
    call obs_to_file_in_situ(obs, (/species/), file_unit)

    map_c%arm = 0.0
    do while (now >= time_start)
      print *, 'Sum of mass map:', sum(map_c%arm) / observation_scaling
      call inject_in_situ(obs, map_c, map_px, map_py, map_pz, now, fu_opposite(timestep), 0, 1, nearest_point)
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
    nullify(obs%durations, obs%obsData, obs%modelData, obs%endTimes, &
         obs%variance, obs%ind_obs2transp, obs%scale_transp2obs, obs%station)

    if (iStat /= 0) then
      call set_error('Deallocation error','destroyObservationInSitu')
    end if
  end subroutine destroyInSitu

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

   function fu_getValueAtStation(observation, mapConc, map_cmX, map_cmY, iInterpolation) result(v)
    !gets a value corresponding to the  observation from the massmap

    implicit none
    real :: v
    type(InSituObservation), intent(in) :: observation
    type(Tmass_map), intent(in)  :: mapConc, map_cmX, map_cmY ! concentrations and mass centres
    integer, intent(in) :: iInterpolation

    !
    integer :: ispecies_transp, ispecies_obs
    real :: fract, weight, cmX, cmY, fStdX, fStdY
    real, parameter :: fZcTrimax = 1.0/6.0  ! Maximum CM for trapezoid slab

    v = 0.  !! Should stay zero if the station fully out of our domain
    ! Whole domain is obtained by smpi_allreduce_add
    ! Warning! Stations not filled by anyone will be silently reported as zeros
    !
    ! There can be several ways to deal with data extraction.
    ! For the time being, we allow either grid-cell average (a.c.a. nearest_neighbour)
    ! or interpolated to the mass centre point. In the latter case, we limit the variation
    ! of the mass-centre to the extent when trapecoid turns into triangle: empty space
    ! in the grid cell is not allowed. Hence, -1/6 < cm < 1/6
    !
    if (any((/observation%station%ix_dispersion, observation%station%iy_dispersion/) == int_missing)) return
    select case(iInterpolation)
      case(nearest_point)
        weight = 1.
        do ispecies_obs = 1, observation%num_obs_species
          ispecies_transp = observation%ind_obs2transp(ispecies_obs)
          fract = observation%scale_transp2obs(ispecies_obs)
          v = v + mapconc%arm(ispecies_transp, 1, observation%iLevel, observation%station%ix_dispersion, &
                                                  & observation%station%iy_dispersion) * fract * weight
        enddo !iSp

      case(linear)
        call set_error('Linear interpolation between the grid cells is considered wrong','fu_getValueAtStation')
        return

      case(toMassCentreLinear)
        do ispecies_obs = 1, observation%num_obs_species
          ispecies_transp = observation%ind_obs2transp(ispecies_obs)
          fract = observation%scale_transp2obs(ispecies_obs)
          cmX = max(-fZcTrimax, min(fZcTrimax, map_cmX%arm(ispecies_transp, 1, observation%ilevel, &
                                                            & observation%station%ix_dispersion, &
                                                            & observation%station%iy_dispersion)))
          cmY = max(-fZcTrimax, min(fZcTrimax, map_cmY%arm(ispecies_transp, 1, observation%ilevel, &
                                                            & observation%station%ix_dispersion, &
                                                            & observation%station%iy_dispersion)))
          fStdX = observation%station%xLoc_dispGrid
          fStdY = observation%station%yLoc_dispGrid
          ! Slope on each dimension is: 12 * cm
          weight = (1. + 12. * cmX * fStdX) * (1. + 12. * cmY * fStdY)
          v = v + mapconc%arm(ispecies_transp, 1, observation%iLevel, observation%station%ix_dispersion, &
                                                  & observation%station%iy_dispersion) * fract * weight
        enddo !iSp
      case default
        call set_error('Unknown interpolation type:' + fu_str(toMassCentreLinear),'fu_getValueAtStation')
        return
    end select

  end function fu_getValueAtStation


  !**********************************************************************************
  !
  ! Accessor functions
  !
  !**********************************************************************************

  function fu_stationPtr(obs) result(station)
    implicit none
    type(inSituObservation) :: obs
    type(observationStation), pointer :: station
    station => obs%station
  end function fu_stationPtr

  function fu_label(stat)
    implicit none
    type(observationStation), intent(in) :: stat
    character(len=4) :: fu_label
    fu_label = stat%label
  end function fu_label

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

  integer function fu_size_column(obs) result(n)
    implicit none
    type(t_column_observation), intent(in) :: obs
    n = obs%n_values
  end function fu_size_column

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
