MODULE silam_times ! times and time-intervals + toolbox


  ! Description:
  ! Contains the silja_time and silja_interval.
  !
  ! Author: Mika Salonoja, FMI email Mika.Salonoja@fmi.fi
  ! 
  ! Language: ANSI Fortran 90
  !
  ! All units: SI
  ! 
  ! Modules used:
!  USE globals
!  USE toolbox
  use globals !, only : r4k, r8k
  use silam_namelist

  IMPLICIT NONE

  ! The public functions and subroutines available in this module:
  PUBLIC defined
  PUBLIC fu_set_time_utc
  public fu_start_of_day_utc
  PUBLIC get_time_components_utc
  public fu_str
  public fu_computer_time_string ! The only function deals with local time
                                 !  Use for debug/report proposes only
  PUBLIC fu_time_fname_string_utc ! general-purpose function
  PUBLIC fu_time_to_date_string ! silja_time to yyyymmdd
  PUBLIC fu_time_to_grads_string ! silja_time to  00:00Z10mar2001
  public fu_time_to_netcdf_string
  public fu_time_to_thredds_string
  PUBLIC fu_io_string_to_time ! yyyy mm dd hh min sec to silja_time
  public fu_time_to_io_string ! reverse to the above
  PUBLIC fu_hirlam_string_to_time ! yymmddhh to silja_time
  PUBLIC fu_interval_string_filename
  public fu_interval_to_named_string
  PUBLIC fu_interval_to_grads_string
  PUBLIC fu_wallclock
  PUBLIC msg_time
  PUBLIC msg_time_test
  PUBLIC time_tests
  PUBLIC fu_set_interval_sec
  PUBLIC fu_set_interval_min
  PUBLIC fu_set_interval_h
  PUBLIC fu_set_interval_d
  public fu_set_named_interval
  PUBLIC fu_opposite
  PUBLIC fu_sec
  PUBLIC fu_sec8
  PUBLIC fu_min
  PUBLIC fu_hour
  PUBLIC fu_day
  public fu_mon
  public fu_year
  PUBLIC fu_mins
  PUBLIC fu_hours
  PUBLIC fu_days
  public fu_weekday
  PUBLIC fu_abs
  PUBLIC fu_interval_positive
  PUBLIC fu_round_closest_hour
  PUBLIC fu_julian_date
  public fu_julian_date_real
  public fu_julian_date_to_time
  PUBLIC fu_real_hours_utc
  PUBLIC fu_total_cpu_time ! since the start of program
  PUBLIC times_to_ascending_order
  PUBLIC fu_closest_time
  PUBLIC fu_latest_time
  PUBLIC fu_earliest_time
  PUBLIC fu_between_times
  public fu_time_overlap
  public fu_next_special_time
  public fu_days_in_month
  PUBLIC SolarSetup
  public fu_solar_zenith_angle_cos
  public fu_daylength
  public fu_sunrise_utc
  PUBLIC report
  PUBLIC fu_round_closest
  public update_timezone_offsets
  public fu_index

  public real8_to_silja_time
  public silja_time_to_real8

  ! The private functions and subroutines not to be used elsewhere:
  private fu_time_to_string
  private fu_interval_string
  PRIVATE fu_time_defined
  private fu_set_interval_sec_from_real
  private fu_set_interval_sec_from_double
  PRIVATE fu_interval_defined
  PRIVATE fu_interval_fraction
  PRIVATE fu_divide_interval
  PRIVATE fu_compare_intervals_gt
  PRIVATE fu_compare_intervals_ge
  PRIVATE fu_compare_intervals_lt
  PRIVATE fu_compare_intervals_le
  PRIVATE fu_compare_intervals_eq ! interface ==
  PRIVATE fu_sub_time_and_interval ! interface -
  PRIVATE fu_add_time_and_interval ! interface +
  private fu_add_intervals
  private fu_multiply_interval_int
  private fu_multiply_interval_real
  PRIVATE fu_time_difference ! interface -
  PRIVATE fu_compare_times_eq ! interface ==
  PRIVATE fu_compare_times_ne ! interface /=
  PRIVATE fu_compare_times_gt ! interface >
  PRIVATE fu_compare_times_ge ! interface >=
  PRIVATE fu_compare_times_lt ! interface <
  PRIVATE fu_compare_times_le ! interface <=
  PRIVATE fu_compare_times ! Used only for the comparision above.
  PRIVATE fu_time_ok
  PRIVATE fu_closest_time_from_times
  private fu_sec_time
  private fu_sec_interval
  PRIVATE report_time
  PRIVATE report_times
  PRIVATE report_interval
  private fu_solar_zenith_angle_cos_time


  ! Generic names and operator-interfaces of some functions:

  INTERFACE defined
    MODULE PROCEDURE  fu_time_defined, fu_interval_defined
  END INTERFACE

  interface fu_str
    module procedure fu_time_to_string
    module procedure fu_interval_string
  end interface
  
  interface fu_set_interval_sec
    module procedure fu_set_interval_sec_from_real
    module procedure fu_set_interval_sec_from_double
  end interface

  INTERFACE report
    MODULE PROCEDURE  report_time, report_interval, report_times
  END INTERFACE

  interface fu_index
     module procedure fu_index_of_time
  end interface
  
  INTERFACE fu_closest_time
    MODULE PROCEDURE  fu_closest_time_from_times
  END INTERFACE

  INTERFACE operator(+)
    MODULE PROCEDURE  fu_add_time_and_interval, fu_add_intervals
  END INTERFACE

  INTERFACE operator(-)
    MODULE PROCEDURE  fu_time_difference,&
                    & fu_sub_time_and_interval,&
                    & fu_subtract_intervals
  END INTERFACE

  INTERFACE operator(/)
    MODULE PROCEDURE fu_interval_fraction, fu_divide_interval
  END INTERFACE

  INTERFACE operator(*)
    MODULE PROCEDURE fu_multiply_interval_int, fu_multiply_interval_real
  END INTERFACE

  INTERFACE operator(/=)
    MODULE PROCEDURE fu_compare_times_ne
  END INTERFACE

  INTERFACE operator(==)
    MODULE PROCEDURE fu_compare_intervals_eq,&
                   & fu_compare_times_eq
  END INTERFACE

  INTERFACE operator(>)
    MODULE PROCEDURE fu_compare_times_gt, &
                   & fu_compare_intervals_gt
  END INTERFACE

  INTERFACE operator(<)
    MODULE PROCEDURE fu_compare_intervals_lt,&
                   & fu_compare_times_lt
  END INTERFACE

  INTERFACE operator(>=)
    MODULE PROCEDURE fu_compare_times_ge,&
                   & fu_compare_intervals_ge
  END INTERFACE

  INTERFACE operator(<=)
    MODULE PROCEDURE fu_compare_times_le,&
                   & fu_compare_intervals_le
  END INTERFACE

  interface fu_sec
    module procedure fu_sec_time
    module procedure fu_sec_interval
  end interface
  
  interface fu_solar_zenith_angle_cos
    module procedure fu_solar_zenith_angle_cos_time
  end interface

 ! Public types with private components defined in this module:

   TYPE silja_interval
    PRIVATE
    REAL(r8k) :: ii = real_missing ! in seconds, of course
    TYPE(silja_logical) :: defined = silja_false
  END TYPE silja_interval

  TYPE silja_time ! date and time
    PRIVATE
    INTEGER :: year, month, day, hour, minute
    REAL(r8k) :: sec
    TYPE(silja_logical) :: defined = silja_false
  END TYPE silja_time


  ! Public parameters:
  INTEGER*8, PARAMETER, PUBLIC :: i_seconds_in_min = 60
  INTEGER*8, PARAMETER, PUBLIC :: i_seconds_in_hour = 3600
  INTEGER*8, PARAMETER, PUBLIC :: i_seconds_in_day = 86400
  INTEGER*8, PARAMETER, PUBLIC :: i_seconds_in_week = 604800

  REAL(4), PARAMETER, PUBLIC :: seconds_in_day = 86400.
  REAL(4), PARAMETER, PUBLIC :: seconds_in_week = 604800.

  integer, dimension(13), parameter :: &
       & days_before_month = (/0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365/), &
       & days_before_month_leapyear = (/0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366/)

  real(r8k), parameter, public :: d_zero=0.
  real(r8k), PARAMETER, PUBLIC :: d_seconds_in_min = 60.
  real(r8k), PARAMETER, PUBLIC :: d_seconds_in_hour = 3600.
  real(r8k), PARAMETER, PUBLIC :: d_seconds_in_day = 86400.
  real(r8k), PARAMETER, PUBLIC :: d_seconds_in_week = 604800.
  real(r8k), PARAMETER, PUBLIC :: d_seconds_in_month = 2592000.  ! Roughly, 30 days
  real(r8k), PARAMETER, PUBLIC :: d_seconds_in_year = 31557600.   ! Roughly, 365.25 days

  TYPE(silja_interval), PARAMETER, PUBLIC :: zero_interval = silja_interval(0., silja_true)
  TYPE(silja_interval), PARAMETER, PUBLIC :: one_second = silja_interval(1., silja_true)
  TYPE(silja_interval), PARAMETER, PUBLIC :: ten_seconds = silja_interval(10., silja_true)
  TYPE(silja_interval), PARAMETER, PUBLIC :: one_minute = silja_interval(60., silja_true)
  TYPE(silja_interval), PARAMETER, PUBLIC :: one_hour = silja_interval(3600., silja_true)
  TYPE(silja_interval), PARAMETER, PUBLIC :: two_hour = silja_interval(7200., silja_true)
  TYPE(silja_interval), PARAMETER, PUBLIC :: three_hour = silja_interval(10800., silja_true)
  TYPE(silja_interval), PARAMETER, PUBLIC :: six_hours = silja_interval(21600., silja_true)
  TYPE(silja_interval), PARAMETER, PUBLIC :: twelve_hour = silja_interval(43200., silja_true)
  TYPE(silja_interval), PARAMETER, PUBLIC :: one_day = silja_interval(seconds_in_day, silja_true)
  TYPE(silja_interval), PARAMETER, PUBLIC :: one_week = silja_interval(seconds_in_week, silja_true)
  TYPE(silja_interval), PARAMETER, PUBLIC :: very_long_interval = silja_interval(12622780416., silja_true) !!Some 400 years
 
  !
  ! Some flags of the time type: there can be normal or regular time and can be some special
  ! time moments, which happen regular but cannot be described in terms of ordinary time
  !
  INTEGER, PARAMETER, PUBLIC :: iRegular = 8000
  INTEGER, PARAMETER, PUBLIC :: iStartOfHour = 8001
  INTEGER, PARAMETER, PUBLIC :: iStartOfDay = 8002
  INTEGER, PARAMETER, PUBLIC :: iStartOfWeek = 8003
  INTEGER, PARAMETER, PUBLIC :: iStartOfMonth = 8004
  INTEGER, PARAMETER, PUBLIC :: iStartOfYear = 8005
  INTEGER, PARAMETER, PUBLIC :: iEndOfDay = 8006
  INTEGER, PARAMETER, PUBLIC :: iEndOfWeek = 8007
  INTEGER, PARAMETER, PUBLIC :: iEndOfMonth = 8008
  INTEGER, PARAMETER, PUBLIC :: iEndOfYear = 8009
  INTEGER, PARAMETER, PUBLIC :: iMeteoTime = 8010


 
!!!!!!!!!!!!!#define XXX ABCDR   Uncoment to enable buggy time difference

#ifdef XXX
 integer(8), parameter, private :: TM_YEAR_BASE = 1900  !Dull legacy of unix struct_tm
#endif
  integer, parameter, private :: EPOCH_YEAR=1970   !UNIX epoch
  integer, parameter, private :: MIN_ALLOWED_YEAR=0      !There was nothing  before Christ!
  integer, parameter, private :: MAX_ALLOWED_YEAR=4000   !We all will be dead and forgotten by then...

  ! The missing codes of times:
  TYPE(silja_interval), PARAMETER, PUBLIC :: interval_missing = silja_interval(real_missing, &
                                                                             & silja_false)

  TYPE(silja_time), PARAMETER :: time_missing = silja_time(int_missing,&
                                                         & int_missing,&
                                                         & int_missing,&
                                                         & int_missing,&
                                                         & int_missing,&
                                                         & real_missing,&
                                                         & silja_false)

  TYPE(silja_time), PARAMETER, public :: epoch = silja_time(1970,&
                                                          & 1,&
                                                          & 1,&
                                                          & 0,&
                                                          & 0,&
                                                          & 0.,&
                                                          & silja_true)

  TYPE(silja_time), PARAMETER, public :: really_far_in_past = silja_time(1800,&
                                                                       & 1,&
                                                                       & 1,&
                                                                       & 0,&
                                                                       & 0,&
                                                                       & 0.,&
                                                                       & silja_true)
  TYPE(silja_time), PARAMETER, public :: really_far_in_future = silja_time(2200,&
                                                                         & 1,&
                                                                         & 1,&
                                                                         & 0,&
                                                                         & 0,&
                                                                         & 0.,&
                                                                         & silja_true)
  type(silja_time), parameter, public :: ref_time_01012000 = silja_time(2000,&
                                                                       & 1,&
                                                                       & 1,&
                                                                       & 12,&  ! mid-day
                                                                       & 0,&
                                                                       & 0.,&
                                                                       & silja_true)

  type(silja_time), parameter, public :: ref_time_30121974 = silja_time(1974,&
                                                                       & 12,&
                                                                       & 30,&
                                                                       & 0,&  
                                                                       & 0,&
                                                                       & 0.,&
                                                                       & silja_true)
  type(silja_time), public, parameter :: tDebugTime = silja_time(2000,&
                                                               & 1,&
                                                               & 21,&
                                                               & 5,&
                                                               & 0,&
                                                               & 0.,&
                                                               & silja_true)
  type(silja_interval), public, parameter :: max_forecast_length = silja_interval(3.6D3 * 9.99D2, &
                                                                                & silja_true)
  integer, parameter, private :: weekday_01012000 = 6 ! Saturday
  !
  ! For averaged and period-valid variables, the period can be labeled with just one time
  ! but we would need to know whether it marks the start or end of the period.
  !
  integer, public, parameter :: start_of_period = 1101
  integer, public, parameter :: mid_of_period = 1102
  integer, public, parameter :: end_of_period = 1103
  integer, public, parameter :: instant_fields = 1104


  ! arrays corresponding to tz indices
  ! To be filled from namelist
  integer, public :: number_of_time_zones
  integer, private, dimension(:), allocatable :: TZ_offset_summer
  integer, private, dimension(:), allocatable :: TZ_offset_winter
  integer, public, parameter :: tz_name_length = 256
  character(len= tz_name_length), dimension(:), public, allocatable :: TZ_name          ! TZdata names
                                                              ! e.g.
                                                              ! Europe/Helsinki
  character(len= tz_name_length), dimension(:), public, allocatable :: TZ_windows_name
  character(len= 10), dimension(:), public, allocatable :: TZ_code          ! TZdata codes, e.g. FIN

  ! current offset, to be updated every day or something
  integer, public, dimension(:), allocatable :: TZ_offset
  integer, public, dimension(24) :: solar_TZ_index !! Indices in above atrray
                                                   !! For Solar TZones 
                                                   !! UTC, 3600 ...

  integer, private :: timezone_method
!  integer, private, parameter :: mean_timezone_method = 66601
  integer, private, parameter :: summer_timezone_method = 66602
  integer, private, parameter :: winter_timezone_method = 66603
  integer, private, parameter :: lastsunday_timezone_method = 66604
  integer, private, parameter :: system_timezone_method = 66605


  !================== Ashrae/Iqbal stuff  to simulate Clear-sky radiation
  ! Some variables which are dependent only on day of year and GMT time
  ! - hence they do not vary from grid to grid and can safely be stored here

  real, private, save :: rdecl,sinrdecl,cosrdecl
  real, private, save :: eqtime, eqt_h, tan_decl

  type, public :: ashrae_tab
        integer :: nday
        real    :: a
        real    :: b
        real    :: c
  end type ashrae_tab

  type(ashrae_tab), save, public ::  Ashrae    ! Current values

  type(ashrae_tab), parameter, dimension(14), private ::  ASHRAE_REV = (/ &
    !              nday    a      b      c
         ashrae_tab(  1, 1203.0, 0.141, 0.103 ) &
        ,ashrae_tab( 21, 1202.0, 0.141, 0.103 ) &
        ,ashrae_tab( 52, 1187.0, 0.142, 0.104 ) &
        ,ashrae_tab( 81, 1164.0, 0.149, 0.109 ) &
        ,ashrae_tab(112, 1130.0, 0.164, 0.120 ) &
        ,ashrae_tab(142, 1106.0, 0.177, 0.130 ) &
        ,ashrae_tab(173, 1092.0, 0.185, 0.137 ) &
        ,ashrae_tab(203, 1093.0, 0.186, 0.138 ) &
        ,ashrae_tab(234, 1107.0, 0.182, 0.134 ) &
        ,ashrae_tab(265, 1136.0, 0.165, 0.121 ) &
        ,ashrae_tab(295, 1136.0, 0.152, 0.111 ) &
        ,ashrae_tab(326, 1190.0, 0.144, 0.106 ) &
        ,ashrae_tab(356, 1204.0, 0.141, 0.103 ) &
        ,ashrae_tab(366, 1203.0, 0.141, 0.103 ) &
      /)
  !=============================



CONTAINS 

  ! ****************************************************************

  ! ****************************************************************
  ! 
  !
  !      HERE STARTS THE SETTING, RETURNING AND CHECKING OF TIMES
  !
  !
  ! ***************************************************************
  ! ***************************************************************



  ! ****************************************************************


  FUNCTION fu_set_time_utc(year, month, day, hour, minute, sec) result(time)

    ! Sets a value for a timevariable, and checks if the values of
    ! timecomponents are legal. The time is always converted to utc,
    ! which is the internal format of silja.
    !
    ! Mika Salonoja, FMI, 05-1995
    !
    IMPLICIT NONE

    ! Return value of this function:
    TYPE(silja_time) :: time

    ! Imported parameters with intent(in)::
    INTEGER, INTENT(in) :: year, month, day, hour, minute
    REAL, INTENT(in) :: sec

    time%year = year
    time%month = month
    time%day = day
    time%hour = hour
    time%minute = minute
    time%sec = sec

    IF (.NOT.fu_time_ok(time)) THEN
      CALL set_error('Illegal time', 'fu_set_time')
      RETURN
    END IF

    time%defined = fu_set_true()

  END FUNCTION fu_set_time_utc



  !************************************************************


  LOGICAL FUNCTION fu_time_ok(time)
    !
    ! Checks, if the time given as a parameter is reasonable
    !
    ! Mika Salonoja, FMI, 05-1995
    !
    IMPLICIT NONE

    ! Imported parameters with intent(in):
    TYPE(silja_time), INTENT(in) :: time

    !
    ! 1. Check if the time components make sense.

    fu_time_ok = (time%month >= 1) .and. (time%month <= 12)
    if(fu_time_ok) fu_time_ok = (time%day >= 1) .and. &
                               &(time%day <= (fu_days_in_month(time%month, time%year)))
    if(fu_time_ok) fu_time_ok = (time%hour >= 0) .and. (time%hour <= 23)
    if(fu_time_ok) fu_time_ok = (time%minute >= 0) .and. (time%minute <= 59)
    if(fu_time_ok) fu_time_ok = (time%sec >= 0.) .and. (time%sec < 60.)

    if(.not. fu_time_ok) CALL msg_time ('Not a legal time', time)

    !----------------------------------------------------------
    ! 2. Check that year is given in four digits (2000 is here soon!)
    ! NO CAN DO In netcdf files time count can start from year 0. Year 99 literally means year 99.

!    IF (time%year < 100) THEN
!      CALL msg('Give year in four digits!')
!      fu_time_ok = .false.
!      RETURN
!    END IF

  END FUNCTION fu_time_ok


  ! ***************************************************************


  LOGICAL FUNCTION fu_time_defined(time)

    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_time), INTENT(in) :: time

    fu_time_defined = fu_true(time%defined)

  END FUNCTION fu_time_defined


  ! ***************************************************************


  SUBROUTINE get_time_components_utc (time, year, month, day, hour, minute, sec)

    ! Description:
    ! Return the components of a time-variable in utc. Since this is
    ! the internal format of time in silja, no conversions are needed.
    !
    ! Method:
    ! 
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    TYPE(silja_time), INTENT(in) :: time

    ! Imported parameters with intent OUT:
    INTEGER, INTENT(out) :: year, month, day, hour, minute
    REAL, INTENT(out) :: sec

    IF (.NOT.defined(time)) THEN
      CALL set_error('time not defined', 'get_time_components_utc')
      RETURN
    END IF

    year = time%year
    month = time%month
    day = time%day
    hour = time%hour
    minute = time%minute
    sec = time%sec

  END SUBROUTINE get_time_components_utc



  !****************************************************************

  FUNCTION fu_time_to_string (time, ifBlanks_)
    ! Converts a given time to a human readable string without spaces. 
    ! -- if time%sec < 10, last character will be blank!
    IMPLICIT NONE

    ! Return value of this function:  
    CHARACTER (LEN=21) :: fu_time_to_string 

    ! Imported parameters
    TYPE(silja_time), INTENT(in) :: time
    logical, intent(in), optional :: ifBlanks_    ! if
    
    ! Local variables
    logical :: ifBlanks
    character(len=4) :: seconds_ljust
    
    if(present(ifBlanks_))then
      ifBlanks = ifBlanks_
    else
      ifBlanks = .true.
    endif

    if ( .not. fu_true(time%defined)) then
      WRITE(fu_time_to_string,*) "undefined_time"
      return
    endif

    write(seconds_ljust, fmt='(F4.1)') time%sec
    seconds_ljust = adjustl(seconds_ljust)

    if(ifBlanks)then
      WRITE(fu_time_to_string, fmt = '(I4.4, A, I2.2, A, I2.2, A, I2.2, A, I2.2, A, A)') &
           & time%year, '-', &
           & time%month, '-', &
           & time%day, ' ', &
           & time%hour, ':', &
           & time%minute, ':', &
           & seconds_ljust
    else

      WRITE(fu_time_to_string, fmt = '(I4.4, A, I2.2, A, I2.2, A, I2.2, A, I2.2, A, A)') &
          & time%year, '_', &
          & time%month, '_', &
          & time%day, '_', &
          & time%hour, '.', &
          & time%minute, '.', &
          & seconds_ljust
    endif

  END FUNCTION fu_time_to_string



  !****************************************************************


  FUNCTION fu_time_fname_string_utc(time)

    ! Converts a given time to a string suitable for filename
    ! creation. The value is returned in utc.
    ! 
    IMPLICIT NONE

    ! Return value of this function:  
    CHARACTER (LEN=30) :: fu_time_fname_string_utc

    ! Imported parameters
    TYPE(silja_time), INTENT(in) :: time

    ! Local  declarations:
    CHARACTER(LEN=30) :: long
    integer :: i

    IF (.NOT.defined(time)) THEN
      CALL set_error('undefined time','fu_time_fname_string_utc')
      RETURN
    END IF
    long = fu_str(time)

    ! Replace empty spots with zeroes:
    DO i=1,len_trim(long)
      if(long(i:i) == ' ' .or. long(i:i) == ':' .or. long(i:i) == '.') long(i:i) = '_'
    END DO

    fu_time_fname_string_utc = long + '_UTC'
    
  END FUNCTION fu_time_fname_string_utc


  !****************************************************************


  FUNCTION fu_time_to_date_string (time) result(string)
    !
    ! Converts a given time to a 8-character date string YYYYMMDD
    !
    IMPLICIT NONE

    ! Return value of this function:  
    CHARACTER (LEN=8) :: string

    ! Imported parameters
    TYPE(silja_time), INTENT(in) :: time
    
    ! Local variables
    integer :: i

    IF (.NOT.defined(time)) THEN
      CALL set_error('undefined time given','fu_time_to_date_string')
      RETURN
    END IF

    WRITE(string, fmt = '(I4.4,2I2.2)') time%year, time%month, time%day

  END FUNCTION fu_time_to_date_string


  !****************************************************************


  FUNCTION fu_time_to_daily_time_string (time) result(string)
    !
    ! Converts a given time to a 6-character date string HHMMSS
    !
    IMPLICIT NONE

    ! Return value of this function:  
    CHARACTER (LEN=6) :: string

    ! Imported parameters
    TYPE(silja_time), INTENT(in) :: time

    ! Local variables
    integer :: i

    IF (.NOT.defined(time)) THEN
      CALL set_error('undefined time given','fu_time_to_daily_time_string')
      RETURN
    END IF

    WRITE(string, fmt = '(3I2)') time%hour, time%minute, int(time%sec+0.5)

    ! Replace empty spots with zeroes:
    DO i=1,6
      IF (string(i:i) == ' ') string(i:i) = '0'
    END DO

  END FUNCTION fu_time_to_daily_time_string



  !****************************************************************


  FUNCTION fu_time_to_grads_string (time) result(string)
    
    ! Converts a given time to a GrADS string. The type is: hh:mmZddmmmyyyy
    ! For example 06:00Z25apr1996 for 25.4.1996 at 06.00 UTC
    !
    IMPLICIT NONE

    ! Return value of this function:  
    CHARACTER (LEN=15) :: string

    ! Imported parameters
    TYPE(silja_time), INTENT(in) :: time

    ! Local variables:
    INTEGER :: year_of_century, i
    
    IF (.NOT.defined(time)) THEN
      CALL set_error('undefined time given' ,'fu_time_to_grads_string')
      RETURN
    END IF

    WRITE(string, fmt='(I2,A1,I2,A1,I2,A3,I4)') &
            & time%hour,':', &
            & time%minute,'Z', &
            & time%day, &
            & chMonthNames_3chr(time%month), &
            & time%year

      ! Replace empty spots with zeroes:
    DO i=1,15
      IF (string(i:i) == ' ') string(i:i) = '0'
    END DO

  !  PRINT*,string,' GrADS'

  END FUNCTION fu_time_to_grads_string


  !*****************************************************************

    FUNCTION fu_time_to_netcdf_string(time)
    ! Converts a given time to udunits-legal string
    ! Roux
    IMPLICIT NONE

    ! Return value of this function:  
    CHARACTER (LEN=23) :: fu_time_to_netcdf_string

    ! Imported parameters
    TYPE(silja_time), INTENT(in) :: time

    WRITE(fu_time_to_netcdf_string, fmt = '(I4.4, A, I2.2, A, I2.2, A, I2.2, A, I2.2, A, I2.2, 4A)') &
        & time%year, '-', &
        & time%month, '-', &
        & time%day, ' ', &
        & time%hour, ':', &
        & time%minute, ':', &
        & nint(time%sec), ' UTC'

  END FUNCTION fu_time_to_netcdf_string
  !*****************************************************************

    FUNCTION fu_time_to_thredds_string(time)
    ! Converts a given time to thredds string
    ! Roux
    IMPLICIT NONE

    ! Return value of this function:  
    CHARACTER (LEN=23) :: fu_time_to_thredds_string

    ! Imported parameters
    TYPE(silja_time), INTENT(in) :: time

    WRITE(fu_time_to_thredds_string, fmt = '(I4.4, A, I2.2, A, I2.2, A, I2.2, A, I2.2, A, I2.2, 4A)') &
        & time%year, '-', &
        & time%month, '-', &
        & time%day, 'T', &
        & time%hour, ':', &
        & time%minute, ':', &
        & nint(time%sec), 'Z'

  END FUNCTION fu_time_to_thredds_string


  !*****************************************************************

  FUNCTION fu_hirlam_string_to_time (string) result(time)

    ! Converts a given 8-character string used in hirlam
    ! and other nwp-filenames to time. The type of input-string is: 
    ! yymmddhh in which year is given in two numbers (year of
    ! century). For example 96042506 for 25.4.1996 at 06.00 UTC
    ! 
    ! Mika Salonoja, FMI
    !
    IMPLICIT NONE

    ! Return value of this function:  
    TYPE(silja_time) :: time

    ! Imported parameters
    CHARACTER (LEN=*), INTENT(in) :: string

    ! Local variables:
    INTEGER :: i, io_status
    INTEGER, DIMENSION(4) :: components

    READ (unit=string, fmt='(4I2)', iostat=io_status) (components(i), i=1,4)

    IF (io_status /= 0) THEN
      CALL set_error(fu_connect_strings('something strange in given time:',&
                                       & string), 'fu_hirlam_string_to_time')
    END IF


    ! Check year:
    IF (components(1) > 50) THEN
      components(1) = components(1) + 1900
    ELSE
      components(1) = components(1) + 2000
    END IF

    time = fu_set_time_utc(components(1),&
                         & components(2),&
                         & components(3),&
                         & components(4),&
                         & 0,&
                         & 0.)

    IF (error) THEN
      CALL set_error(fu_connect_strings('something strange in given time:',&
                                       & string), 'fu_hirlam_string_to_time')
    END IF

  END FUNCTION fu_hirlam_string_to_time


  !****************************************************************


  FUNCTION fu_io_string_to_time (string) result(time)

    ! Converts a given string to utc-time. The type of input-string can be: 
    ! 1.  [yyyy mm dd hh min sec] where seconds are real and the rest integer
    ! 2.  single word "NOW", which means that the release has started just at 
    !     the last meteoperiod (one of: 00, 06, 12, 18)
    ! 
    ! Original code: Mika Salonoja
    ! Author: Mikhail Sofiev, FMI
    !
    IMPLICIT NONE

    ! Return value of this function:  
    TYPE(silja_time) :: time

    ! Imported parameters
    CHARACTER (LEN=*), INTENT(in) :: string

    ! Local variables:
    INTEGER :: i, io_status
    INTEGER, DIMENSION(5) :: components
    REAL :: sec

    IF(fu_str_u_case(string) == 'NOW')THEN
      time = fu_wallclock()  ! time will be just NOW, but in UTC, of course

    ELSEIF(fu_str_u_case(string) == 'LAST_METEO_TIME')THEN
      time = fu_wallclock()  ! time will be the last-most meteoperiod, in UTC, of course
      DO i = 0,18,6
        IF(time%hour >= i .and. time%hour < i+6)THEN
          time%hour = i
          EXIT
        END IF
      END DO
      time%minute = 0
      time%sec = 0.

    elseif(fu_str_u_case(string) == 'FAR_IN_PAST')then
      time = really_far_in_past

    elseif(fu_str_u_case(string) == 'FAR_IN_FUTURE')then
      time = really_far_in_future

    elseif(fu_str_u_case(string) == 'VOID_TIME')then
      time = time_missing

    ELSE
      READ(unit=string, fmt=*, iostat=io_status) (components(i),i=1,5), sec

      IF (io_status /= 0) THEN
        CALL set_error(fu_connect_strings('something strange in given string:',&
                                         & string), 'fu_io_string_to_time')
        RETURN
      END IF

      time = fu_set_time_utc(components(1),&
                           & components(2),&
                           & components(3),&
                           & components(4),&
                           & components(5),&
                           & sec)
      IF (error) THEN
        CALL set_error(fu_connect_strings('cannot make time out of:', string), &
                    & 'fu_io_string_to_time')
      END IF
    END IF
  END FUNCTION fu_io_string_to_time


  !******************************************************************

  function fu_time_to_io_string(time)result(str)
    !
    ! Converts SILAM time to yyyy mm dd hh mm sec. UTC
    !
    implicit none

    ! Return value of this function:  
    CHARACTER (LEN=25) :: str

    ! Imported parameters
    TYPE(silja_time), INTENT(in) :: time

    if(time == really_far_in_past)then
      str = 'FAR_IN_PAST'
    elseif(time == really_far_in_future)then
      str = 'FAR_IN_FUTURE'
    else
      WRITE(str, fmt = '(I4.4, 1x, I2.2, 1x, I2.2, 1x, I2.2, 1x, I2.2, 1x, F4.1, A4)') &
                  & time%year, time%month, time%day, time%hour, time%minute, time%sec, ' UTC'
    endif
  end function fu_time_to_io_string

  ! ***************************************************************


  REAL FUNCTION fu_real_hours_utc(time)
    !
    ! Description:
    ! Returns the clokctime of time as a real value of hours. For
    ! example 15.45:00 becomes 15,75000.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_time), INTENT(in) :: time

    fu_real_hours_utc = REAL(time%hour) + (REAL(time%minute)/60.) + (time%sec/3600.)

  END FUNCTION fu_real_hours_utc


  ! ***************************************************************


  INTEGER FUNCTION fu_julian_date(time)

    ! Description:
    ! Returns the Julian date of time, which is the order number of
    ! the day inside the current year.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_time), INTENT(in) :: time
    !
    ! Local declarations:
    INTEGER :: i

    january: IF (time%month == 1) THEN

      fu_julian_date = time%day

    ELSE

      fu_julian_date = 0

      DO i = 1, time%month - 1
          fu_julian_date = fu_julian_date + fu_days_in_month(i, time%year)
      END DO

      fu_julian_date = fu_julian_date + time%day

    END IF january


  END FUNCTION fu_julian_date


  ! ***************************************************************


  function fu_julian_date_real(time)

    ! Description: Returns the Julian date of time as real number, or more properly, the
    ! number of days elapsed from 1 January 00:00. Eg. 1.5 is 2 January, 12:00. The integer
    ! Julian date is ceiling(real_julian_date).
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_time), INTENT(in) :: time

    real(r8k) :: fu_julian_date_real

    ! Local declarations:
    INTEGER :: i

    fu_julian_date_real = time%day-1 + real(time%hour, r8k)/24 &
         & + real(time%minute, r8k)/(24*60) + time%sec/(24*60*60)
    do i = 1, time%month - 1
      fu_julian_date_real = fu_julian_date_real + fu_days_in_month(i, time%year)
    enddo

  end function fu_julian_date_real


  ! ***************************************************************


  function fu_julian_date_to_time (julian_date, year)
    ! Returns SILAM time from the 'real' julian date defined above. The reverse
    ! calculation is sensitive to roundoff errors: take care.
    
    real(r8k), intent(in) :: julian_date
    integer, intent(in) :: year
    TYPE(silja_time) :: fu_julian_date_to_time

    integer :: month, day, hour, min, days_in_month
    real(r8k) :: sec, jDate
    real(r8k), parameter :: epsilon = 1e-5

    if (julian_date < 0 .or. 365 <= julian_date) then
      call msg('Julian date given:', julian_date)
      call set_error('Bad Julian date', 'fu_julian_date_to_time')
      return
    end if

    jDate = julian_date
    do month = 1, 12
      if (fu_days_in_month(month, year) <= jDate)then
        jDate = jDate - fu_days_in_month(month, year)
      else
        exit
      endif
    enddo

    day = ceiling(jDate + epsilon)
    jDate = jDate - floor(jDate) ! == fractional part of jDate

    jDate = jDate * 24
    hour = floor(jDate + epsilon)
    jDate = jDate - hour
    jDate = jDate * 60
    min = floor(jDate + epsilon)
    jDate = jDate - min
    jDate = jDate * 60

    !hour = jDate * 24
    !jDate = jDate - (real(hour, 8)+epsilon)/24

    !min = jDate * 1440
    !jDate = jDate - (real(min, 8)+epsilon)/1440

    !sec = jDate * 86400


    !print *, 'sec', jDate
    if (fu_fails(jdate > -epsilon, 'jDate too inaccurate', 'fu_julian_date_to_time')) return
    sec = max(0.0_r8k, jDate)
    !sec = jDate
    
    !hour = jDate * 24
    !jDate = jDate - real(hour, 8)/24

    !min = jDate * 1440
    !jDate = jDate - real(min, 8)/1440

    !sec = jDate * 86400

    fu_julian_date_to_time = fu_set_time_utc(year, month, day, hour, min, real(sec))

  end function fu_julian_date_to_time

  ! ****************************************************************





  integer function fu_weekday(time)
    !
    ! Just computes the local time adding the UTC difference
    !
    implicit none

    type(silja_time), intent(in) :: time

    real :: fTmp
    integer :: iTmp

    iTmp = int(abs(fu_sec(time - ref_time_01012000))/86400. + 0.5)+weekday_01012000
    fu_weekday = mod(iTmp,7)
    if(fu_weekday == 0) fu_weekday = 7


!    fTmp = abs(fu_hour(time - ref_time_01012000))/24  ! in days
!    iTmp = int(fTmp + 0.5)+weekday_01012000
!    fu_weekday = mod(iTmp,7)
!!    fu_weekday = mod(abs(fu_day(time - ref_time_01012000)) + weekday_01012000,7)
!    if(fu_weekday == 0) fu_weekday = 7

  end function fu_weekday


  ! *****************************************************************
  ! ****************************************************************
  ! 
  !
  !        TIME-INTERVAL STUFF: SET AND RETURN VALUES
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  FUNCTION fu_set_interval_sec_from_real(real_interval)

    ! Description:
    ! Sets a value for a time interval given in seconds. The value
    ! must be given either as a real-type variable.
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_interval) :: fu_set_interval_sec_from_real
    !
    ! Imported parameters with intent(in):
    REAL(4), INTENT(in):: real_interval

    fu_set_interval_sec_from_real%ii = real_interval
    fu_set_interval_sec_from_real%defined = fu_set_true()

  END FUNCTION fu_set_interval_sec_from_real



  ! ***************************************************************

  FUNCTION fu_set_interval_sec_from_double(double_interval)

    ! Description:
    ! Sets a value for a time interval given in seconds. The value
    ! must be given either as a real-type variable.
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_interval) :: fu_set_interval_sec_from_double
    !
    ! Imported parameters with intent(in):
    real(r8k), INTENT(in):: double_interval

    fu_set_interval_sec_from_double%ii = double_interval
    fu_set_interval_sec_from_double%defined = fu_set_true()

  END FUNCTION fu_set_interval_sec_from_double



  ! *****************************************************************


  FUNCTION fu_set_interval_min(interval)

    ! Description:
    ! Sets a value for a time interval given in minutes. The value
    ! must be given as an integer-type variable (full minutes).
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_interval) :: fu_set_interval_min
    !
    ! Imported parameters with intent(in):
    INTEGER, INTENT(in):: interval

    fu_set_interval_min%ii = REAL(interval)*60.
    fu_set_interval_min%defined = fu_set_true()


  END FUNCTION fu_set_interval_min



  ! *****************************************************************


  FUNCTION fu_set_interval_h(interval)

    ! Description:
    ! Sets a value for a time interval given in hours. The value
    ! must be given as an integer-type variable (full hours).
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_interval) :: fu_set_interval_h
    !
    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: interval

    fu_set_interval_h%ii = REAL(interval)*3600.
    fu_set_interval_h%defined = fu_set_true()


  END FUNCTION fu_set_interval_h



  ! *****************************************************************


  FUNCTION fu_set_interval_d(interval)

    ! Description:
    ! Sets a value for a time interval given in days. The value
    ! must be given as an integer-type variable (full days).
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_interval) :: fu_set_interval_d
    !
    ! Imported parameters with intent(in):
    real, INTENT(in) :: interval

    fu_set_interval_d%ii = interval * seconds_in_day
    fu_set_interval_d%defined = fu_set_true()

  END FUNCTION fu_set_interval_d


  ! *****************************************************************


  FUNCTION fu_set_named_interval(chInterval, start_time) result(interval)
    !
    ! Sets a value for a time interval from the character string: "named 
    ! interval". The named interval includes a number and the unit, separated
    ! by space. Unit can be omitted, then it is assumed to be seconds, 
    ! according to SI system
    !
    ! Such a named interval can be e.g. "1 mon", which means that it will
    ! depend on start moment due to varying number of days. Therefore, one can
    ! supply this moment so that more exact consideration can be done.
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_interval) :: interval
    !
    ! Imported parameters with intent(in):
    character(len=*), INTENT(in) :: chInterval
    type(silja_time), intent(in), optional :: start_time
    
    !
    ! Local variables
    integer :: io_status
    real(r8k) :: fInt
    character(len=10) :: chUnit
    type(silja_time) :: end_time

    !
    ! First, check for specific key word "INFINITY". Then set something very long and exit
    !
    if(fu_str_u_case(chInterval) == 'INFINITY')then
      interval = very_long_interval
      return
    endif
    
    !
    ! Try to get the number and the unit
    !
    fInt = -1.
    chUnit = ''
    read(unit=chInterval,fmt=*,iostat=io_status) fInt, chUnit

    if(io_status /= 0)then
      call set_error(fu_connect_strings('Strange named interval:',chInterval), &
                   & 'fu_set_named_interval')
    elseif(chUnit == '')then
      interval%ii = fInt
      interval%defined = fu_set_true()
    elseif(fu_str_u_case(chUnit) == 'YR')then
      interval%ii = fInt * seconds_in_day * 365.25    ! Approx days in year
      interval%defined = fu_set_true()
      if(present(start_time))then      ! Can we compute it exactly ?
        if(defined(start_time) .and. (real(int(fInt+0.5),8).eps.fInt))then
          end_time = start_time
          end_time%year = end_time%year + int(fInt+0.5)
          interval = end_time - start_time
        endif
      endif
    elseif(fu_str_u_case(chUnit) == 'MON')then
      interval%ii = fInt * seconds_in_day * 30.4366667  ! Approx. days in month
      interval%defined = fu_set_true()
      if(present(start_time))then      ! Can we compute it exactly ?
        if(defined(start_time) .and. (real(int(fInt+0.5),8).eps.fInt))then
          end_time = start_time
          end_time%month = end_time%month + int(fInt+0.5)
          io_status = int(end_time%month / 12)
          end_time%year = end_time%year + io_status
          end_time%month = end_time%month - io_status * 12
          interval = end_time - start_time
        endif
      endif
    elseif(fu_str_u_case(chUnit) == 'WEEK')then
      interval%ii = fInt * 7. * seconds_in_day        ! Exact days in week
      interval%defined = fu_set_true()
    elseif(fu_str_u_case(chUnit) == 'DAY')then
      interval%ii = fInt * seconds_in_day
      interval%defined = fu_set_true()
    elseif(fu_str_u_case(chUnit) == 'HR')then
      interval%ii = fInt * 3600.
      interval%defined = fu_set_true()
    elseif(fu_str_u_case(chUnit) == 'MIN')then
      interval%ii = fInt * 60.
      interval%defined = fu_set_true()
    elseif(fu_str_u_case(chUnit) == 'SEC')then
      interval%ii = fInt
      interval%defined = fu_set_true()
    else
      call set_error(fu_connect_strings("Can't parse named interval:",chInterval), &
                   & 'fu_set_named_interval')

    endif
    
  END FUNCTION fu_set_named_interval


  ! ***************************************************************


  FUNCTION fu_opposite(interval)

    ! Description:
    ! Changes the sign of a time-interval.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_interval) :: fu_opposite
    !
    ! Imported parameters with intent(in):
    TYPE(silja_interval), INTENT(in) :: interval

    fu_opposite%ii = -1. * interval%ii
    fu_opposite%defined = fu_set_true()

  END FUNCTION fu_opposite


  ! ***************************************************************


  FUNCTION fu_abs(interval)

    ! Description:
    ! Returns the absolut non-negative value of interval.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_interval) :: fu_abs
    !
    ! Imported parameters with intent(in):
    TYPE(silja_interval), INTENT(in) :: interval

    IF (defined(interval)) THEN
      fu_abs%ii = ABS(interval%ii)
      fu_abs%defined = fu_set_true()
    ELSE
      CALL set_error('undefined interval given','fu_abs')
    END IF

  END FUNCTION fu_abs


  ! ***************************************************************


  real(r8k) FUNCTION fu_sec8(interval)

    ! Description:
    ! Returns the value of a silja_interval in seconds.
    !
    ! Imported parameters with intent(in):
    TYPE(silja_interval), INTENT(in) :: interval

    IF (defined(interval)) THEN
      fu_sec8 = interval%ii ! it is the internal
      ! format already
    ELSE
      CALL set_error('undefined interval given','fu_sec8')
      RETURN
    END IF

  END FUNCTION fu_sec8


  ! ***************************************************************


  REAL FUNCTION fu_sec_interval(interval)
    !
    ! Returns the value of a silja_interval in seconds.
    !
    ! Imported parameters with intent(in):
    TYPE(silja_interval), INTENT(in) :: interval

    IF (defined(interval)) THEN
      fu_sec_interval = real(interval%ii) ! it is the internal
      ! format already
    ELSE
      fu_sec_interval = real_missing
      CALL set_error('undefined interval given','fu_sec_interval')
      RETURN
    END IF
    
  END FUNCTION fu_sec_interval


  ! ***************************************************************

  real FUNCTION fu_mins(interval)

    ! Description:
    ! Returns the value of a silja_interval in full minutes.
    !
    ! Imported parameters with intent(in):
    TYPE(silja_interval), INTENT(in) :: interval

    IF (defined(interval)) THEN
      fu_mins = interval%ii/60.
    ELSE
      CALL set_error('undefined interval given','fu_interval_min')
      fu_mins = real_missing
      RETURN
    END IF
    
  END FUNCTION fu_mins
  
  
  ! ***************************************************************

  real FUNCTION fu_hours(interval)

    ! Description:
    ! Returns the value of a silja_interval in full hours.
    !
    ! Imported parameters with intent(in):
    TYPE(silja_interval), INTENT(in) :: interval

    IF (defined(interval)) THEN
      fu_hours = interval%ii/3600.
    ELSE
      CALL set_error('undefined interval given','fu_interval_hour')
      fu_hours = real_missing
      RETURN
    END IF
    
  END FUNCTION fu_hours
  
  
  ! ***************************************************************

  real FUNCTION fu_days(interval)
    !
    ! Returns the value of a silja_interval in full days.
    !
    ! Imported parameters with intent(in):
    TYPE(silja_interval), INTENT(in) :: interval
    
    IF (defined(interval)) THEN
      fu_days = interval%ii/seconds_in_day
    ELSE
      fu_days = real_missing
      CALL set_error('undefined interval given','fu_interval_day')
      RETURN
    END IF
    
  END FUNCTION fu_days
  
  
  ! ***************************************************************

  REAL FUNCTION fu_sec_time(time)
    !
    ! Returns the value of a silja_time in seconds.
    !
    ! Imported parameters with intent(in):
    TYPE(silja_time), INTENT(in) :: time

    IF (defined(time)) THEN
      fu_sec_time = time%sec
    ELSE
      CALL set_error('undefined time given','fu_sec')
      RETURN
    END IF
    
  END FUNCTION fu_sec_time

  
  ! ***************************************************************

  INTEGER FUNCTION fu_min(time)

    ! Description:
    ! Returns the value of a silja_time in full minutes.
    !
    ! Imported parameters with intent(in):
    TYPE(silja_time), INTENT(in) :: time

    IF (defined(time)) THEN
      fu_min = time%minute
    ELSE
      CALL set_error('undefined time given','fu_time_min')
      RETURN
    END IF
    
  END FUNCTION fu_min
  
  
  ! ***************************************************************

  INTEGER FUNCTION fu_hour(time)

    ! Description:
    ! Returns the value of a silja_time in full hours.
    !
    ! Imported parameters with intent(in):
    TYPE(silja_time), INTENT(in) :: time

    IF (defined(time)) THEN
      fu_hour = time%hour
    ELSE
      CALL set_error('undefined time given','fu_time_hour')
      RETURN
    END IF
    
  END FUNCTION fu_hour
  
  
  ! ***************************************************************

  
  INTEGER FUNCTION fu_day(time)
    
    ! Description:
    ! Returns the value of a silja_time in full days.
    !
    ! Imported parameters with intent(in):
    TYPE(silja_time), INTENT(in) :: time
    
    IF (defined(time)) THEN
      fu_day = time%day
    ELSE
      CALL set_error('undefined time given','fu_time_day')
      RETURN
    END IF
    
  END FUNCTION fu_day
  
  ! ***************************************************************

  
  INTEGER FUNCTION fu_mon(time)
    
    ! Description:
    ! Returns the value of a silja_time in full days.
    !
    ! Imported parameters with intent(in):
    TYPE(silja_time), INTENT(in) :: time
    
    IF (defined(time)) THEN
      fu_mon = time%month
    ELSE
      CALL set_error('undefined time given','fu_time_mon')
      RETURN
    END IF
    
  END FUNCTION fu_mon
  
  ! ***************************************************************

  
  INTEGER FUNCTION fu_year(time)
    
    ! Description:
    ! Returns the value of a silja_time in full days.
    !
    ! Imported parameters with intent(in):
    TYPE(silja_time), INTENT(in) :: time
    
    IF (defined(time)) THEN
      fu_year = time%year
    ELSE
      CALL set_error('undefined time given','fu_time_year')
      RETURN
    END IF
    
  END FUNCTION fu_year
  

  !****************************************************************


  FUNCTION fu_interval_string_filename (interval) result(string)

    ! Converts a given interval to a 2-character string in full hours
    ! used in hirlam
    ! and other nwp-filenames. 
    ! 
    ! Mika Salonoja, FMI, 04-1996
    !
    IMPLICIT NONE

    ! Return value of this function:  
    CHARACTER (LEN=2) :: string

    ! Imported parameters
    TYPE(silja_interval), INTENT(in) :: interval

    ! Local declarations:    
    INTEGER :: i

    IF (.NOT.defined(interval)) THEN
      CALL set_error('undefined interval given','fu_interval_string_filename')
      RETURN
    END IF

    WRITE(string, fmt = '(I2)') fu_hours(interval)

    ! Replace empty spots with zeroes:
    DO i=1,2
      IF (string(i:i) == ' ') string(i:i) = '0'
    END DO

  END FUNCTION fu_interval_string_filename


  !****************************************************************


  FUNCTION fu_interval_string (interval) result(string)

    ! Returns a pretty string of an interval. 
    ! 
    ! Mika Salonoja, FMI
    !
    IMPLICIT NONE

    ! Return value of this function:  
    CHARACTER (LEN=30) :: string

    ! Imported parameters
    TYPE(silja_interval), INTENT(in) :: interval

    ! Local declarations:    
    INTEGER :: i
    INTEGER :: hours, minutes
    real(r8k) :: seconds

    IF (.NOT.defined(interval)) THEN
      string = 'undefined interval'
    elseif(interval == very_long_interval)then
      string = 'INFINITY'
    else

      IF (mod(interval%ii, d_seconds_in_hour) .deps. d_zero) THEN
        hours = NINT(interval%ii/3600.)
      ELSE
        hours = INT(interval%ii/3600.)
      END IF

      IF (mod(interval%ii, d_seconds_in_min) .deps. d_zero) THEN
        minutes = NINT(interval%ii/60.) - (hours*60)
      ELSE
        minutes = INT(interval%ii/60.) - (hours*60)
      END IF

      seconds = interval%ii - REAL(minutes)*60. - REAL(hours)*3600.

      IF (hours == 0) THEN

        IF (minutes == 0) THEN
          WRITE(string, fmt = '(F7.4, A)') seconds, 'sec.'
        ELSE
          WRITE(string, fmt = '(I3, A, F7.4, A)') minutes, 'min ', seconds, 'sec.'
        END IF

      ELSE

        IF ((hours < 100).and.(hours>-10)) THEN
          WRITE(string, fmt = '(I3, A, I2, A, F7.4, A)') hours, 'h ', minutes, 'min ', seconds, 'sec.'
        ELSE
          WRITE(string, fmt = '(I8, A, I2, A)') hours, 'h ', minutes, 'min '
        END IF
      END IF
    END IF

  END FUNCTION fu_interval_string

  !****************************************************************


  FUNCTION fu_interval_to_named_string (interval) result(string)

    ! Returns a string of an interval, which later can be consumed by set_named_interval function
    ! 
    IMPLICIT NONE

    ! Return value of this function:  
    CHARACTER (LEN=30) :: string

    ! Imported parameters
    TYPE(silja_interval), INTENT(in) :: interval

    IF (.NOT.defined(interval)) THEN
      string = 'undefined interval'
    elseif(interval == very_long_interval)then
      string = 'INFINITY'
    elseif(interval%ii < 60.)then
      write(unit=string, fmt='(F10.6,2x,A3)')interval%ii,'SEC'
    elseif(interval%ii < 3600.)then
      write(unit=string, fmt='(F10.6,2x,A3)')interval%ii/60.,'MIN'
    elseif(interval%ii < 86400.)then
      write(unit=string, fmt='(F10.6,2x,A2)')interval%ii/3600.,'HR'
    elseif(interval%ii < 3153600000.0)then
      write(unit=string, fmt='(F12.6,2x,A3)')interval%ii/86400.,'DAY'
    else
      write(unit=string, fmt='(D15.8,2x,A3)')interval%ii/86400.,'DAY'
    endif
    if(index(string,'*') > 0)then
      call report(interval)
    endif
  END FUNCTION fu_interval_to_named_string
  
  
  !****************************************************************


  recursive FUNCTION fu_interval_to_grads_string (interval) result(string)

    ! Returns a GrADS short string of an interval up to days.
    ! NOTE. There is no tool to handle monthly increment as it is 
    ! irregular and can not be handled without explicit information
    ! about the dates surrounding the given interval
    ! 
    ! Mikhail Sofiev, FMI
    !
    IMPLICIT NONE

    ! Return value of this function:  
    CHARACTER (LEN=10) :: string

    ! Imported parameters
    TYPE(silja_interval), INTENT(in) :: interval

    IF (.NOT.defined(interval)) THEN
      string = '?????'
    ELSEIF(interval%ii < -1.)THEN ! Negative interval, but it's forbidden in GrADS
      string=fu_interval_to_grads_string (interval * (-1.))
    ELSEIF(interval%ii < 1.)THEN ! Zero interval, but it's forbidden in GrADS
      string='  1hr'
    ELSEIF (mod(interval%ii, d_seconds_in_day) .deps. d_zero) THEN
      WRITE(string,fmt='(I3,A2)')ABS(NINT(interval%ii/86400.)), 'dy'
    ELSEIF (mod(interval%ii, d_seconds_in_hour) .deps. d_zero) THEN
      WRITE(string,fmt='(I3,A2)')ABS(NINT(interval%ii/3600.)), 'hr'
    ELSEIF (mod(interval%ii, d_seconds_in_min) .deps. d_zero) THEN
      WRITE(string,fmt='(I4,A2)')ABS(NINT(interval%ii/60.)), 'mn'
    ELSE
      CALL set_error('Strange time step','fu_interval_to_grads_string')
      RETURN
    END IF

  END FUNCTION fu_interval_to_grads_string


  ! ****************************************************************


  LOGICAL FUNCTION fu_compare_intervals_eq(int1, int2)

    ! Compares two time-intervals. True if equal.
    !
    ! Mika Salonoja, FMI, 05-1995
    !
    IMPLICIT NONE

    TYPE(silja_interval), INTENT(in) :: int1, int2

    fu_compare_intervals_eq = (int1%ii .deps. int2%ii)

  END FUNCTION fu_compare_intervals_eq


  ! ****************************************************************


  LOGICAL FUNCTION fu_compare_intervals_gt(int1, int2)

    ! Compares two time-intervals. True if first greater..
    !
    IMPLICIT NONE

    TYPE(silja_interval), INTENT(in) :: int1, int2

    fu_compare_intervals_gt = (int1%ii > int2%ii)

  END FUNCTION fu_compare_intervals_gt


  ! ****************************************************************


  LOGICAL FUNCTION fu_compare_intervals_ge(int1, int2)

    ! Compares two time-intervals. True if first greater or equal.
    !
    ! Mika Salonoja, FMI, 05-1995
    !
    IMPLICIT NONE

    TYPE(silja_interval), INTENT(in) :: int1, int2

    fu_compare_intervals_ge = (int1%ii >= int2%ii)

  END FUNCTION fu_compare_intervals_ge


  ! ****************************************************************


  LOGICAL FUNCTION fu_compare_intervals_lt(int1, int2)

    ! Compares two time-intervals. True if first smaller.
    !
    ! Mika Salonoja, FMI, 05-1995
    !
    IMPLICIT NONE

    TYPE(silja_interval), INTENT(in) :: int1, int2

    IF (int1%ii < int2%ii) THEN
      fu_compare_intervals_lt = .true.
    ELSE
      fu_compare_intervals_lt = .false.
    END IF

  END FUNCTION fu_compare_intervals_lt


  ! ****************************************************************


  LOGICAL FUNCTION fu_compare_intervals_le(int1, int2)

    ! Compares two time-intervals. True if first smaller or equal.
    !
    IMPLICIT NONE

    TYPE(silja_interval), INTENT(in) :: int1, int2

    fu_compare_intervals_le = (int1%ii <= int2%ii)

  END FUNCTION fu_compare_intervals_le


  ! ****************************************************************


  LOGICAL FUNCTION fu_interval_positive(interval)

    ! Checks the sign of an interval.
    !
    ! Mika Salonoja, FMI, 05-1995
    !
    IMPLICIT NONE

    TYPE(silja_interval), INTENT(in) :: interval

    fu_interval_positive = interval%ii >= 0.

  END FUNCTION fu_interval_positive



  ! ***************************************************************


  REAL FUNCTION fu_interval_fraction(int1, int2) result(frac)

    ! Description:
    ! Calculates int1/int2.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_interval), INTENT(in) :: int1, int2

    IF (int2%ii /= 0.) THEN
      frac = int1%ii/int2%ii
    ELSE
      CALL set_error('tried to divide by a zero time interval',&
                   & 'fu_interval_fraction')
    END IF

  END FUNCTION fu_interval_fraction


  ! ***************************************************************


  FUNCTION fu_divide_interval(interval, divider)

    ! Description:
    ! Calculates int1/int2.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_interval) :: fu_divide_interval

    ! Imported parameters with intent(in):
    TYPE(silja_interval), INTENT(in) :: interval
    REAL, INTENT(in) :: divider

    IF (divider /= 0.) THEN
      fu_divide_interval = fu_set_interval_sec(interval%ii/divider)
    ELSE
      CALL set_error('tried to divide by a zero time interval',&
          & 'fu_divide_interval')
    END IF

  END FUNCTION fu_divide_interval


  ! ***************************************************************


  FUNCTION fu_add_intervals(int1, int2) 

    ! Description:
    ! Calculates int1 + int2.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_interval) :: fu_add_intervals

    ! Imported parameters with intent(in):
    TYPE(silja_interval), INTENT(in) :: int1, int2
    !

    fu_add_intervals%ii = int1%ii + int2%ii
    fu_add_intervals%defined = fu_set_true()

  END FUNCTION fu_add_intervals


  ! ***************************************************************


  FUNCTION fu_subtract_intervals(int1, int2) 

    ! Description:
    ! Calculates int1 - int2.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_interval) :: fu_subtract_intervals

    ! Imported parameters with intent(in):
    TYPE(silja_interval), INTENT(in) :: int1, int2
    !

    fu_subtract_intervals%ii = int1%ii - int2%ii
    fu_subtract_intervals%defined = fu_set_true()

  END FUNCTION fu_subtract_intervals


  ! ***************************************************************


  FUNCTION fu_multiply_interval_real(interval, x) result(int_out)

    ! Description:
    ! Multiplies interval with x.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Return value of this function:
    TYPE(silja_interval) :: int_out

    ! Imported parameters with intent(in):
    TYPE(silja_interval), INTENT(in) :: interval
    REAL, INTENT(in) :: x

    int_out%ii = interval%ii * x
    int_out%defined = interval%defined

  END FUNCTION fu_multiply_interval_real


  ! ***************************************************************


  FUNCTION fu_multiply_interval_int(interval, i) result(int_out)

    ! Multiplies interval with i.
    !
    IMPLICIT NONE
    !
    ! Return value of this function:
    TYPE(silja_interval) :: int_out

    ! Imported parameters with intent(in):
    TYPE(silja_interval), INTENT(in) :: interval
    integer, INTENT(in) :: i

    int_out%ii = interval%ii * i
    int_out%defined = interval%defined

  END FUNCTION fu_multiply_interval_int


  ! ***************************************************************


  LOGICAL FUNCTION fu_interval_defined(interval)

    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_interval), INTENT(in) :: interval

    if(fu_true(interval%defined))then
      if(real(interval%ii) .eps. real_missing)then
        fu_interval_defined = .false.
        call msg_warning('Defined interval == real_missing')
        call msg_warning('Defined interval == real_missing')
        call msg_warning('Defined interval == real_missing')
        call msg_warning('Defined interval == real_missing')
        call msg_warning('Defined interval == real_missing')
        call msg_warning('Defined interval == real_missing')
      else
        fu_interval_defined = .true.
      endif
    else
      fu_interval_defined = .false.
    endif

  END FUNCTION fu_interval_defined






  ! ****************************************************************

  ! ****************************************************************
  ! 
  !
  !           ADDING AND SUBTRACTING TIMES AND TIME-INTERVALS
  !    
  !
  ! ***************************************************************
  ! ***************************************************************

  ! (The following functions are used through interfaces, check the
  ! beginning of this module.)


  FUNCTION fu_sub_time_and_interval(time, interval)
    !
    ! Adds to a time an interval of SI-unit (seconds, that is!)
    !
    ! Mika Salonoja, FMI, 05-1995
    !
    IMPLICIT NONE

    ! Return value of this function:
    TYPE(silja_time) :: fu_sub_time_and_interval

    ! Imported parameters with intent(in): 
    TYPE(silja_time), INTENT(in) :: time
    TYPE(silja_interval), INTENT(in) :: interval

    ! Local declarations    
    fu_sub_time_and_interval = time + (fu_opposite(interval))

  END FUNCTION fu_sub_time_and_interval


  !****************************************************************


  FUNCTION fu_add_time_and_interval(time, interval) result(res)
    !
    ! Adds to a time an interval of SI-unit (seconds, that is!)
    !
    ! Mika Salonoja, FMI
    !
    IMPLICIT NONE

    ! Return value of this function:
    TYPE(silja_time) :: res

    ! Imported parameters with intent(in): 
    TYPE(silja_time), INTENT(in) :: time
    TYPE(silja_interval), INTENT(in) :: interval

    ! Local declarations:
    REAL(r8k) :: rel
    INTEGER :: minutes

    REAL(r8k), PARAMETER :: time_epsilon = 1.0e-08

    IF (.NOT.defined(time)) THEN
      CALL set_error('undefined time given','fu_add_time_and_interval')
      call msg("time%hour:", time%hour)
      RETURN
    END IF

    IF (.NOT.defined(interval)) THEN
      CALL set_error('undefined interval given','fu_add_time_and_interval')
      call msg("Interval:", interval%ii)
      RETURN
    END IF


    !-------------------------------------------------------------
    ! Calculate seconds

    res = time
    res%sec = res%sec + interval%ii

    ! Check if seconds are legal:

    IF ((res%sec < 60.) .and. (res%sec >= 0.)) RETURN

    ! Normalize seconds between 0...59 by adding/subtracting minutes
    !
    ! The amount of minutes is rounded downwards, but
    ! extremely small differences caused by dividing real numbers are
    ! checked and in that case rounded to closest full minute:

    rel = res%sec/60.

    IF (ABS(REAL(NINT(rel)) - rel) < time_epsilon) THEN
      minutes = NINT(rel)
      res%sec = 0.0
    ELSE
      minutes = INT(rel)
      res%sec = MOD(res%sec, 60.)
   END IF

    res%minute = res%minute + minutes

    IF (res%sec < 0.) THEN
      res%sec = res%sec + 60.
      res%minute = res%minute - 1
    END IF
    if(res%sec .eps. 60.)then
      res%sec = 0.
      res%minute = res%minute + 1
    endif

    ! Check if minutes are legal:
    IF ((res%minute >= 0) .and. (res%minute < 60)) RETURN


    ! Normalize minutes between 0...59  by adding/subtracting hours
    res%hour = res%hour + res%minute/60
    res%minute = MOD(res%minute, 60)

    IF (res%minute < 0) THEN
      res%minute = res%minute + 60
      res%hour = res%hour - 1
    END IF

    ! Check if hours are legal:
    IF ((res%hour >= 0) .and. (res%hour < 23)) RETURN


    ! Normalize hours between 0...23 by adding/subtracting days
    res%day = res%day + res%hour/24
    res%hour = MOD(res%hour, 24)

    IF (res%hour < 0) THEN
      res%hour = res%hour + 24
      res%day = res%day - 1
    END IF

    IF ((res%day >= 1) .and. & 
      & (res%day <= fu_days_in_month(res%month, res%year))) RETURN

    ! The month and maybe also the year has to be changed! -------

    !   Forwards in month & year:
    IF (res%day > fu_days_in_month(res%month, res%year)) THEN
      month_fw:      DO ! Move one month at a time forwards until day
        ! is legal:

        res%day = res%day - fu_days_in_month(res%month, res%year)

        res%month = res%month + 1

        ! Next year:
        IF (res%month > 12) THEN
          res%month = 1 ! It has to be January!
          res%year = res%year + 1
        END IF ! Next year

        IF (res%day <= fu_days_in_month(res%month, res%year)) EXIT

      END DO month_fw

      RETURN      
    END IF

    IF (res%day < 1) THEN    !   Backwards in month & year:
      month_bw:      DO ! Move one month at a time backwards until day
        ! is legal:
        res%month = res%month - 1

        ! Previous year:
        IF (res%month < 1) THEN
          res%month = 12 ! It has to be December!
          res%year = res%year - 1
        END IF

        res%day = res%day + fu_days_in_month(res%month, res%year)

        IF (res%day >= 1) EXIT

      END DO month_bw

    END IF

  END FUNCTION fu_add_time_and_interval

  !***********************************************************************
  
  FUNCTION fu_round_closest(time, interval, direction) result(time1)
    !
    ! Description:
    ! Returns the latest time that is representable
    ! as  integer number of intrvals since dayStart
    ! Works within a day
    ! Should work porperly when time is turned to integer
    ! 
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Roux, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function: 
    TYPE(silja_time) :: time1
    ! Imported parameters with intent(in):
    TYPE(silja_time), INTENT(in) :: time
    INTEGER, INTENT(in) :: direction
    TYPE(silja_interval),  INTENT(in) :: interval
    
    TYPE(silja_time) :: dayStart ! Dirty workaround
    real ::  fSteps

    dayStart = fu_start_of_day_utc(time)
    fSteps = fu_sec(time - dayStart) / fu_sec(interval)
!    call msg("*************************truncating time")
!    call msg("Request:  "+fu_str(time))
!    call msg("Interval: "+fu_interval_string(interval))
!    call msg("Daystart: "+fu_str(dayStart))
   !OBS fSteps can be zero
    if(abs(fSteps - nint(fSteps,4)) <=  fSteps * 1.e-5)then  ! Hit exactly the obstime
      time1 = time
!        call msg("Match!")
    else
      if(direction == forwards)then   ! left or right edge of the range?
        time1 = dayStart + (interval * ceiling(fSteps))
!        call msg("Forwards")
      elseif(direction == backwards)then
!        call msg("Backwards")
        time1 = dayStart + (interval * floor(fSteps))
      elseif(direction == back_and_forwards)then
!        call msg("Closest")
        time1 = dayStart + (interval * nint(fSteps,4))
      else
        call set_error('Strange direction in time:' + fu_str(direction),'fu_round_closest')
        return
      endif
    endif

! call msg(fu_str(time1))
  
  END FUNCTION fu_round_closest

  !***********************************************************************
  
  function real8_to_silja_time(val) result (time)
    ! Code inspired by minix gmtime, extended for crazy values
    ! http://www.cise.ufl.edu/~cop4600/cgi-bin/lxr/http/source.cgi/lib/ansi/gmtime.c
    ! Roux
    TYPE(silja_time) :: time
    real(r8k), INTENT(in) :: val  
    real(r8k) :: dayclock
    integer :: idayno, idaymin, imon
    integer :: year
    integer :: iLeapYear 

    integer, dimension(12,2), parameter :: ytab = &
         & reshape((/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,&
         &            31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /), (/12,2/))
    integer, dimension(2), parameter :: YEARSIZE = (/365,366/) \

    time%defined = silja_true
        
    idayno = floor(val / d_seconds_in_day, kind=8)
    dayclock = val - idayno*i_seconds_in_day
    idaymin = floor(dayclock / 60.) 
 
    time%sec    =  dayclock - idaymin*60
    time%minute =  mod(idaymin, 60)
    time%hour   =  idaymin / 60;

    !timep->tm_wday = (dayno + 4) % 7;       /* day 0 was a thursday */
    if (idayno >= 0) then
       do year = EPOCH_YEAR, MAX_ALLOWED_YEAR
               if (mod(year,4) /= 0) then
                 iLeapYear = 1 !non-leap
               elseif (mod(year,100) == 0)then
                  iLeapYear = 1 !Non-leap
                  if (mod(year, 400) == 0) iLeapYear = 2 !Leap
               else
                  iLeapYear = 2 !Leap
               endif
               if (idayno < YEARSIZE(iLeapYear)) exit
               idayno = idayno - YEARSIZE(iLeapYear)
       enddo
    else !Before epoch
       do year = EPOCH_YEAR-1,MIN_ALLOWED_YEAR,-1 
               if (mod(year,4) /= 0) then
                 iLeapYear = 1 !non-leap
               elseif (mod(year,100) == 0)then
                  iLeapYear = 1 !non-leap
                  if (mod(year, 400) == 0) iLeapYear = 2 !Leap
               else
                  iLeapYear = 2 !Leap
               endif
               idayno = idayno + YEARSIZE(iLeapYear)
               if (idayno >= 0) exit
       enddo
    endif

    if (year < MIN_ALLOWED_YEAR .or. year >  MAX_ALLOWED_YEAR) then
       time%defined = silja_false
       call msg("real8_to_silja_time val=", val )
       call set_error("Got crazy value for year:"+fu_str(year),"real8_to_silja_time")
    endif
    time%year = year

    do imon = 1,12
            if (idayno < ytab(imon,ileapyear)) exit
            idayno = idayno - ytab(imon,ileapyear)
    enddo
    time%month = imon
    time%day   = idayno + 1;

  end function real8_to_silja_time

  
  ! ***************************************************************
  
   function silja_time_to_real8(time) result (val)

! Code adapted from 
! http://www.opensource.apple.com/source/lukemftp/lukemftp-3/lukemftp/libukem/timegm.c
! with dropped leap seconds
! Roux
    real(r8k) :: val  
    TYPE(silja_time), INTENT(in) :: time
    integer :: y, nleapdays
    integer, dimension(12), parameter :: moff = &
         & (/0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334/)

   !y = tm->tm_year + TM_YEAR_BASE - (tm->tm_mon < 2);
   y = time%year
   if (time%month <= 2) y = y - 1

   nleapdays = y / 4 - y / 100 + y / 400 - &
       & ((EPOCH_YEAR-1) / 4 - (EPOCH_YEAR-1) / 100 + (EPOCH_YEAR-1) / 400);

   !t = ((((time_t) (tm->tm_year - (EPOCH_YEAR - TM_YEAR_BASE)) * 365 +
   !      moff[tm->tm_mon] + tm->tm_mday - 1 + nleapdays) * 24 +
   !   tm->tm_hour) * 60 + tm->tm_min) * 60 + tm->tm_sec;
   val = (((int(time%year - EPOCH_YEAR, kind=8) * 365 + &
        & moff(time%month) + time%day - 1 + nleapdays) * 24 +  &
        & time%hour) * 60 + time%minute) * 60 + time%sec;

   end  function silja_time_to_real8


  !***********************************************************************

  FUNCTION fu_time_difference(time1, time2) result(interval)
    !
    ! Description:
    ! Returns the difference of two times. This function
    ! is used through the interface -
    !
    ! the difference
    ! is positive, if time1 is greater, that is, in future when
    ! compared to time2.
    !
    ! The additional time-difference caused by local/utc-time is
    ! also taken into account. 
    !
    ! Method:
    ! If the month and year are the same for both times, the
    ! calculation is trivial. If not, AIKAERO is not used any more.
    !
    ! The time+interval -tools provided elsewhere in this module are
    ! used for changing the months
    ! 
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function: 
    TYPE(silja_interval) :: interval

    ! Imported parameters with intent(in):
    TYPE(silja_time), INTENT(in) :: time1, time2

    ! Local declarations:
    REAL :: diff = 0. ! The difference in seconds
    REAL :: the_sign
    INTEGER :: days, hours, minutes, ierr, tm_year_1, tm_year_2
    integer(8) :: itime1, itime2
    real(r8k) :: dtime1,dtime2


  

    IF (.NOT.defined(time1)) THEN
      CALL set_error('first time not defined','fu_time_difference')
      RETURN
    END IF

    IF (.NOT.defined(time2)) THEN
      CALL set_error('second time not defined','fu_time_difference')
      RETURN
    END IF

#ifndef XXX
    dtime1 = silja_time_to_real8(time1)
    dtime2 = silja_time_to_real8(time2)
    interval = fu_set_interval_sec_from_double(dtime1 - dtime2)



#else
    ! The code below has a bug 
    ! I do not want to debug it. Use system functions instead
    IF ((time1%year == time2%year).and.&
      & (time1%month == time2%month)) THEN

      diff = (time1%sec - time2%sec) + &
          &  REAL((time1%minute - time2%minute)*60 + & 
          &  (time1%hour - time2%hour)*3600 + &
          &  (time1%day - time2%day)*24*3600)

    ELSE

      tm_year_1 = time1%year - TM_YEAR_BASE
      tm_year_2 = time2%year - TM_YEAR_BASE

      if(tm_year_1 < 0)then
        call set_error('Time 1:' + fu_str(time1) + '- is less than the base time:' + &
                     & fu_str(int(TM_YEAR_BASE,4)), 'fu_time_difference')
        return
      endif
      if(tm_year_2 < 0)then
        call set_error('Time 2:' + fu_str(time2) + '- is less than the base time:' + &
                     & fu_str(int(TM_YEAR_BASE,4)), 'fu_time_difference')
        return
      endif
      
      IF (time1%year > time2%year .or. &
        & (time1%year == time2%year .and. time1%month > time2%month)) THEN
        !the_sign = 1.
        diff = timediff_sec(tm_year_1, &
                          & julian_day(tm_year_1, time1%month, time1%day), &
                          & time1%hour, time1%minute, 0, &
                          & tm_year_2, &
                          & julian_day(tm_year_2, time2%month, time2%day), &
                          & time2%hour, time2%minute, 0)
      ELSE
        diff = -timediff_sec(tm_year_2, &
                           & julian_day(tm_year_2, time2%month, time2%day), &
                           & time2%hour, time2%minute, 0, &
                           & tm_year_1, &
                           & julian_day(tm_year_1, time1%month, time1%day), &
                           & time1%hour, time1%minute, 0)

        !the_sign = -1.
      END IF

      diff = diff + time1%sec - time2%sec
      !diff = REAL(minutes) * 60. * the_sign + time1%sec - time2%sec

    END IF


    ! -------------------------------------------------------
    !
    !  2. Set the value for the interval.
    !     ------------------------------

    interval = fu_set_interval_sec(diff)

  contains
    
    integer function julian_day(year, month, mday)
      implicit none
      integer :: year, month, mday
      
      logical :: leap
    
      leap = (iand(year, 3) == 0) &
           & .and. (mod(year, 100) /= 0 &
                  & .or. (iand(year/100, 3) == iand(-TM_YEAR_BASE / 100, 3)))
      
      if (leap) then
        julian_day = days_before_month_leapyear(month) + mday
      else
        julian_day = days_before_month(month) + mday
      end if
    end function julian_day

    !===========================================================================

    integer(8) function timediff_sec(year1, yday1, hour1, min1, sec1, &
       & year0, yday0, hour0, min0, sec0)
      implicit none
      integer(4) :: year1, yday1, hour1, min1, sec1, &
           & year0, yday0, hour0, min0, sec0
      
      integer(4) :: a4, b4, a100, b100, a400, b400, intervening_leap_days, &
           & tyear1, years, days, hours, subtract1, subtract0
      
      integer(8) :: minutes
      
      a4 = ishft(year1, -2) + ishft(TM_YEAR_BASE, -2)
      b4 = ishft(year0, -2) + ishft(TM_YEAR_BASE, -2)
      
      if (iand(year1, 3) == 0) a4 = a4 - 1
      if (iand(year0, 3) == 0) b4 = b4 - 1
      !print *, '1', a4, b4, ishft(year1, -2), ishft(TM_YEAR_BASE, -2), year1
      !print *, '2', a4, b4
      a100 = a4 / 25
      b100 = b4 / 25
      !print *, '3', a100, b100
      if (mod(a4, 25) < 0) a100 = a100 - 1
      if (mod(b4, 25) < 0) b100 = b100 - 1
      !print *, '4', a100, b100
      a400 = ishft(a100, -2)
      b400 = ishft(b100, -2)
      !print *, '5', a100, b100
      intervening_leap_days = (a4 - b4) - (a100 - b100) + (a400 - b400)
      
      tyear1 = year1
      years = tyear1 - year0
      days = 365 * years + yday1 - yday0 + intervening_leap_days
      hours = 24 * days + hour1 - hour0
      minutes = 60 * int(hours,8) + int(min1,8) - int(min0,8)
      timediff_sec = 60 * minutes + int(sec1,8) - int(sec0,8)
      
    end function timediff_sec
#endif

  END FUNCTION fu_time_difference


  ! *****************************************************************

  ! ****************************************************************
  ! 
  !
  !          HERE STARTS THE SECTION FOR COMPARING TWO TIMES
  !
  !
  ! ***************************************************************
  ! ***************************************************************
  !  
  ! (The following functions are used through interfaces, check the
  ! beginning of this module.)


  LOGICAL FUNCTION fu_compare_times_eq(time1, time2)
    !
    ! Compares two times of silja_time. True if equal.
    !
    ! Mikhail Sofiev, FMI, 2003
    !
    IMPLICIT NONE

    TYPE(silja_time), INTENT(in) :: time1, time2

    IF (fu_true(time1%defined)) THEN
      IF (fu_true(time2%defined)) THEN

        ! both defined, compare. To speed-up, first try to go without subtraction
        !
        if(time1%year == time2%year) then
          if(time1%month == time2%month) then
            if(time1%day == time2%day) then
              if(time1%hour == time2%hour)then
                if(time1%minute == time2%minute)then
                  fu_compare_times_eq = (time1%sec .eps. time2%sec)
                  return
                endif
              endif
            endif
          endif
        endif

        fu_compare_times_eq = abs(fu_sec(time1-time2)) < 0.01

      ELSE
        ! 1 defined, 2 not:
        fu_compare_times_eq = .false.
      END IF

    ELSE ! 1 not defined

      IF (fu_true(time2%defined)) THEN
        ! 2 defined, 1 not:
        fu_compare_times_eq = .false.
      ELSE
        ! both undefined so they're equal:
        fu_compare_times_eq = .true.
      END IF

    END IF

  END FUNCTION fu_compare_times_eq


  !*****************************************************************


  LOGICAL FUNCTION fu_compare_times_ne(time1, time2)

    ! Compares two times of silja_time. True if not equal.
    !
    ! Mika Salonoja, FMI, 05-1995
    !
    IMPLICIT NONE

    TYPE(silja_time), INTENT(in) :: time1, time2

    fu_compare_times_ne = .NOT.(fu_compare_times_eq(time1, time2))

  END FUNCTION fu_compare_times_ne


  !*****************************************************************


  LOGICAL FUNCTION fu_compare_times_gt(time1, time2)
    !
    ! Compares two times of silja_time, true if time1 is greater (time1
    ! is in the future compared to time2)
    !
    ! Mika Salonoja, FMI, 05-1995
    !
    IMPLICIT NONE

    TYPE(silja_time), INTENT(in) :: time1, time2

    IF (fu_compare_times(time1, time2) == 1) THEN
      fu_compare_times_gt = .true.
    ELSE
      fu_compare_times_gt = .false.
    END IF

  END FUNCTION fu_compare_times_gt


  !****************************************************************


  LOGICAL FUNCTION fu_compare_times_ge(time1, time2)
    !
    ! Compares two times of silja_time, true if time1 is greater (time1
    ! is in the future compared to time2) or if time1 and time2 are
    ! equal times.
    !
    ! Mika Salonoja, FMI, 05-1995
    !
    IMPLICIT NONE

    TYPE(silja_time), INTENT(in) :: time1, time2

    INTEGER :: comp

    comp = fu_compare_times(time1, time2)

    IF ((comp == 1) .or. (comp == 0)) THEN
      fu_compare_times_ge = .true.
    ELSE
      fu_compare_times_ge = .false.
    END IF

  END FUNCTION fu_compare_times_ge


  !****************************************************************


  LOGICAL FUNCTION fu_compare_times_lt(time1, time2)
    !
    ! Compares two times of silja_time, true if time2 is greater (time1
    ! is in the past compared to time2)
    !
    ! Mika Salonoja, FMI, 05-1995
    !
    IMPLICIT NONE

    TYPE(silja_time), INTENT(in) :: time1, time2

    IF (fu_compare_times(time1, time2) == 2) THEN
      fu_compare_times_lt = .true.
    ELSE
      fu_compare_times_lt = .false.
    END IF

  END FUNCTION fu_compare_times_lt


  !***************************************************************


  LOGICAL FUNCTION fu_compare_times_le(time1, time2)
    !
    ! Compares two times of silja_time, true if time2 is greater (time1
    ! is in the past compared to time2) or if time1 and time2 are
    ! equal times.
    !
    ! Mika Salonoja, FMI, 05-1995
    !
    IMPLICIT NONE

    TYPE(silja_time), INTENT(in) :: time1, time2

    INTEGER :: comp

    comp = fu_compare_times(time1, time2)

    IF ((comp == 2) .or. (comp == 0)) THEN
      fu_compare_times_le = .true.
    ELSE
      fu_compare_times_le = .false.
    END IF

  END FUNCTION fu_compare_times_le


  !***************************************************************


  INTEGER FUNCTION fu_compare_times(time1, time2)
    !
    ! Compares two times of silja_time and return a value indicating
    ! the values of times:
    ! 0: the two times are equal
    ! 1: time1 is greater (time1 is in the future compared to time2)
    ! 2: time2 is greater (time2 is in the future compared to time1)
    !
    ! Later correction:
    ! There is a 0.001 sec insensitive range, inside which times are
    ! considered to be equal
    !
    IMPLICIT NONE
    !
    ! Importes parameters:
    TYPE(silja_time), INTENT(in) :: time1, time2

    ! Local variables
    real :: fDiff

    ! Check definition of times:
    IF ((.NOT.fu_true(time1%defined)).or.(.NOT.fu_true(time2%defined))) THEN
      CALL set_error('><== comparisions possible for defined times','fu_compare_times')
      RETURN
    END IF

    fDiff = fu_sec(time1-time2)

    if(fDiff < 0.001 .and. fDiff > -0.001)then
      fu_compare_times = 0  ! Inside uncertainty range - times are equal
    elseif(fDiff > 0)then
      fu_compare_times = 1
    else
      fu_compare_times = 2
    endif

  END FUNCTION fu_compare_times


  ! *****************************************************************

  ! ****************************************************************
  ! 
  !
  !        TOOLS FOR ROUNDING TIMES
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  FUNCTION fu_round_closest_hour(time, direction) result(newtime)

    ! Description:
    ! Rounds the given time to closest full hour in given direction
    ! (forwards, backwards or absolutely closest in both directions).
    ! The direction values are given in module globals.
    !
    ! If given time happens to be a full hour, then the same given time
    ! is returned unaltered.
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_time) :: newtime

    ! Imported parameters with intent(in):
    TYPE(silja_time), INTENT(in) :: time
    INTEGER, intent(in) :: direction

    IF (.NOT.defined(time)) THEN
      CALL set_error('undefined time given','fu_round_closest_hour')
      RETURN
    END IF

    ! Check full hour
    IF ((time%minute == 0).and.(time%sec.eps.0.)) THEN
      newtime = time
      RETURN
    END IF

    SELECT CASE (direction)

      CASE (backwards)
      newtime = time
      CALL round_backwards(newtime)


      CASE (forwards)
      newtime = time + one_hour
      CALL round_backwards(newtime)


      CASE (back_and_forwards)
      IF (time%minute >= 30) THEN
        ! Round forwards:
        newtime = time + one_hour
        CALL round_backwards(newtime)
      ELSE
        ! Round backwards:
        newtime = time
        CALL round_backwards(newtime)
      END IF

    CASE default
      CALL set_error('unknown direction value','fu_round_closest_hour')

    END SELECT


  CONTAINS

    SUBROUTINE round_backwards(time)

      IMPLICIT NONE

      ! Imported parameters with intent INOUT:
      TYPE(silja_time), INTENT(inout) :: time

      time%minute = 0
      time%sec = 0.

    END SUBROUTINE round_backwards

  END FUNCTION fu_round_closest_hour


  ! *****************************************************************

  ! ****************************************************************
  ! 
  !
  !        THE WALLCLOCK STUFF (THE REAL TIME NOW, THAT IS)
  !
  !
  ! ***************************************************************
  ! ***************************************************************
  FUNCTION fu_computer_time_string ()

    ! Returns machine's time as string.  For reporting only.
    IMPLICIT NONE
    
    ! Local declarations:
    CHARACTER (LEN=29) :: fu_computer_time_string 
    CHARACTER (LEN=8) :: c_date
    CHARACTER (LEN=10) :: c_time
    CHARACTER (LEN=5) :: c_zone
    INTEGER, DIMENSION(8) :: values

    CALL DATE_AND_TIME(c_date, c_time, c_zone, values)

    IF ((c_date == ' ').or.(values(1) == 0)) THEN
      CALL set_error('No real-time clock available!', &
                   & 'fu_wallclock/silja_times')
      RETURN
    END IF

    WRITE(fu_computer_time_string, &
        !         YYYY -  mm    -   dd  ' '  hh   :   mm   :   ss   zone
        &fmt = '(I4.4, A, I2.2, A, I2.2, A, I2.2, A, I2.2, A, I2.2, A5)') &
        & values(1), '-', &
        & values(2), '-', &
        & values(3), '_', &
        & values(5), ':', &
        & values(6), ':', &
        & values(7),c_zone

  END FUNCTION fu_computer_time_string


  FUNCTION fu_wallclock ()

    ! Returns machine's UTC date and time into a timetype-variable. 
    !
    ! Uses the f90 intrinsic function date_and_time.
    ! The returning vector values
    ! contains the data in integer numbers. The values(4) should,
    ! according to the standard, contain the time difference of local
    ! time compared to the UTC-time, but this is good to check when
    ! going to a new machine, or when the time changes between summer
    ! - and wintertime.
    ! 
    !
    ! Obsolete
    ! !!Also wallclock-time is internally converted to utc, of course.
    ! !!Use functions fu_time_string_original and
    ! !!fu_str to return its value in local and
    ! !!utc, respectively.
    !
    IMPLICIT NONE

    ! Return value of this function:  
    TYPE(silja_time) :: fu_wallclock

    ! Local declarations:
    CHARACTER (LEN=8) :: c_date
    CHARACTER (LEN=10) :: c_time
    CHARACTER (LEN=5) :: c_zone
    INTEGER, DIMENSION(8) :: values

    CALL DATE_AND_TIME(c_date, c_time, c_zone, values)

    IF ((c_date == ' ').or.(values(1) == 0)) THEN
      CALL set_error('No real-time clock available!', &
                   & 'fu_wallclock/silja_times')
      RETURN
    END IF

    fu_wallclock = fu_set_time_utc(values(1), &
                             & values(2), &
                             & values(3), &
                             & values(5), &  
                             & values(6), &   
                             & (REAL(values(7)) + (REAL(values(8)) / 1000.))) &
                             - fu_set_interval_min(values(4))

    ! Values (7) contains seconds of minute
    ! Values (8) contains milliseconds of second from 0...999

  END FUNCTION fu_wallclock


  ! *****************************************************************

  ! ****************************************************************
  ! 
  !
  !        CPU-time calculations
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  FUNCTION fu_total_cpu_time()

    ! Description:
    ! Returns the total cpu time used by all the routines since the
    ! start of the program run. Can be also used to calculate cpu
    ! time used by a single routine, loop or section by calling
    ! before and after the execution, and then subtracting the
    ! cputimes:
    ! before = fu_total_cpu_time()
    ! CALL big_job()
    ! after = fu_total_cpu_time()
    ! usage = after - before
    !
    ! Method:
    ! SECOND on Cray
    ! 
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_interval) :: fu_total_cpu_time

    fu_total_cpu_time = fu_set_interval_sec(fu_total_cpu_time_sec())

  END FUNCTION fu_total_cpu_time




  ! *****************************************************************

  ! ****************************************************************
  ! 
  !
  !        SOME OTHER USEFUL SMALL STUFF
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  INTEGER FUNCTION fu_days_in_month (month, year)

    ! Returns the number of days in given month on a given year.
    ! proleptic georgian assumed
    IMPLICIT NONE

    INTEGER, intent(in) :: month, year

    integer, dimension(12), parameter :: nDaysPerMonCrude = (/31,28,31,30,31,30,31,31,30,31,30,31/)

    fu_days_in_month = nDaysPerMonCrude(month)

    ! Check for the extra day in February:
    IF (month == 2 .and. MOD(year, 4) == 0) then
      if ( all(MOD(year, 400) /= (/100,200,300/))) fu_days_in_month = 29 !Leap year
    endif


  END FUNCTION fu_days_in_month


  ! ****************************************************************


  SUBROUTINE  msg_time(message, time)

    ! Tells a time message to the user. The first dummy version of
    ! doing it.
    ! Example of a call:
    ! CALL msg_time(' Now it is: ', time_now)
    ! 
    !
    ! Mika Salonoja, FMI 06-1995

    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(in) :: message 
    TYPE(silja_time), INTENT(in) :: time
    ! Local variables:
!    CHARACTER (LEN=24) :: fu_time_to_string

    PRINT *, message, fu_str(time)

  END SUBROUTINE msg_time


  ! ****************************************************************

  SUBROUTINE  msg_time_test(message, time)

    ! Tells a time message to the user. The first dummy version of
    ! doing it.
    ! Example of a call:
    ! CALL msg_time_test(' Now it is: ', time_now)
    ! 
    !
    ! Mika Salonoja, FMI 06-1995

    IMPLICIT NONE

    CHARACTER (LEN=*) :: message 
    TYPE(silja_time) :: time

    IF (test_messages) CALL msg_time(message, time)


  END SUBROUTINE msg_time_test





  ! ***************************************************************

  ! ***************************************************************
  !
  !
  !         Arrange time-vectors
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  SUBROUTINE times_to_ascending_order(times, ordered_count)

    ! Description
    ! Arranges a set of times to ascending order, so that the
    ! earliest time is placed to vector's location 1 and latest to
    ! the last. 
    !
    ! If there are undefined times mixed in times-vector, they are
    ! skipped and the defined ones are put in order at the
    ! beginning of the resulting vector.
    ! 
    ! If the same time is found several times in the list, then in
    ! the returning list this time occurs only once, in its correct
    ! place.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent INOUT:
    TYPE(silja_time), DIMENSION(:), INTENT(inout) :: times

    ! Optional parameters with intent OUT:
    INTEGER, INTENT(out), OPTIONAL :: ordered_count ! the number of
    ! defined, different times in the list after ordering

    ! Local declarations:
    TYPE(silja_time), DIMENSION(size(times)):: newtime
    INTEGER, DIMENSION(size(times)) :: order
    INTEGER :: i, j, n
    INTEGER :: status, earlier_count
    LOGICAL :: found_earlier

    ! -----------------------------------------------------------
    !
    ! 2. Get rid of undefined times and several occurrences of the
    !    same time.
    !    -----------------------------------------------

    n = 0 ! number of defined different times
    newtime = time_missing

    outer: DO i = 1, SIZE(times)

      IF (.NOT.defined(times(i))) CYCLE outer

      ! Check that it hasn't been found earlier in the list:
      found_earlier = .false.
      IF (i > 1) THEN
        inner: DO j = 1, (i-1)
          IF (times(i) == times(j)) THEN
            found_earlier = .true.
            EXIT inner
          END IF
        END DO inner
      END IF

      ! If defined and not earlier onhte list, add it to new:
      IF (.NOT.found_earlier) THEN
        n = n+1
        newtime(n) = times(i)
      END IF

    END DO outer


    ! -----------------------------------------------------------
    !
    ! 3. Count the earlier times for each defined time.
    !    ---------------------------------------------

    order = 0
    DO i = 1, n
      DO j = 1, n
        IF (i==j) CYCLE 
        IF (newtime(i) > newtime(j)) order(i) = order(i) + 1
      END DO
    END DO


    ! -----------------------------------------------------------
    !
    ! 4. Create the output list by order
    !    -------------------------------

    DO i = 1, n
      times(order(i)+1) = newtime(i)
    END DO

    ! Fill the end with missing time:
    IF (n < SIZE(times)) times((n+1):) = time_missing

    IF (PRESENT(ordered_count)) ordered_count = n

  END SUBROUTINE times_to_ascending_order


  ! ***************************************************************


  FUNCTION fu_earliest_time(times)
    
    ! Description:
    ! From a list of times find the earliest defined time.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_time) :: fu_earliest_time
    !
    ! Imported parameters with intent(in):
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: times
  
    ! Local declarations:
    INTEGER :: i

    fu_earliest_time = really_far_in_future

    DO i = 1, SIZE(times)
      IF (.NOT.defined(times(i))) CYCLE
      IF (times(i) < fu_earliest_time) fu_earliest_time = times(i)
    END DO

  END FUNCTION fu_earliest_time



  ! ***************************************************************


  FUNCTION fu_latest_time(times)
    
    ! Description:
    ! From a list of times find the earliest defined time.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_time) :: fu_latest_time
    !
    ! Imported parameters with intent(in):
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: times
  
    ! Local declarations:
    INTEGER :: i

    fu_latest_time = really_far_in_past

    DO i = 1, SIZE(times)
      IF (.NOT.defined(times(i))) CYCLE
      IF (times(i) > fu_latest_time) fu_latest_time = times(i)
    END DO

  END FUNCTION fu_latest_time

 !************************************************************************************
  
  integer function fu_index_of_time(time, time_list) result(ind)
    implicit none
    type(silja_time), intent(in) :: time
    type(silja_time), dimension(:), intent(in) :: time_list

    integer :: ii

    ind = int_missing

    do ii = 1, size(time_list)
      if (time_list(ii) == time) then
        ind = ii
        return
      end if
    end do
    
  end function fu_index_of_time
 

  ! ***************************************************************


  LOGICAL FUNCTION fu_between_times(time, lim1, lim2, accept_boundaries_too)
    
    ! Description:
    ! Returns true value, if given time is between given boundaries,
    ! or equal to either of them.
    ! The limits can be in any order. All must be defined.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_time), INTENT(in), target :: time, lim1, lim2
    LOGICAL, INTENT(in) :: accept_boundaries_too ! if false, then this
    ! function returns true value only if given time is between limits and
    ! not equal to either of them.
    !
    ! Local declarations:
    TYPE(silja_time), pointer :: back, fwd

    IF (lim1 < lim2) THEN
      back => lim1
      fwd => lim2
    ELSE
      back => lim2
      fwd => lim1
    END IF

    IF (accept_boundaries_too) THEN
      fu_between_times = ((time>=back).and.(time<=fwd))
    ELSE
      fu_between_times = ((time>back).and.(time<fwd))
    END IF

  END FUNCTION fu_between_times

  
  !*********************************************************************************
  
  function fu_time_overlap(earlier_1, later_1, earlier_2, later_2) result(overlap)
    !
    ! Computes overlap of two time intervals 1 and 2
    !
    implicit none
    type(silja_time), intent(in) :: earlier_1, later_1, earlier_2, later_2
    type(silja_interval) :: overlap
    
    ! INternal variables
    type(silja_time) :: start, end

    if (earlier_2 >= later_1 .or. later_2 <= earlier_1) then
      overlap = zero_interval
      return
    end if
    
    if (later_1 > later_2) then
      end = later_2
    else
      end = later_1
    endif
    
    if (earlier_2 > earlier_1) then
      start = earlier_2
    else
      start = earlier_1
    end if

    overlap = end - start

  end function fu_time_overlap


  !*********************************************************************************

  function fu_next_special_time(timeNow, iOccasionType, iDirection, timeQuantum_) result(timeOccasion)
    !
    ! Having the given time, one can easily compute the next "special occasion", such
    ! as the end of this very week or beginning of the next month. This very task is 
    ! done here. Note: if the current time is not such an occasion, the next one may belong
    ! to the same period. Thus, search for end of the day will give you this very day if done
    ! in the morning. We will get the next day only if time-now is closer to midnight that
    ! one time quantum (e.g. one model time step)
    !
    implicit none

    ! Return value of the function
    type(silja_time) :: timeOccasion

    ! Imported parameters
    type(silja_time), intent(in) :: timeNow
    integer, intent(in) :: iOccasionType, iDirection
    type(silja_interval), intent(in) :: timeQuantum_

    ! Local variables
    integer :: nQuanta
    logical :: ifClose
    type(silja_interval) :: step
    type(silja_interval) :: timeQuantum

    !
    ! Stupidity check
    !
    timeOccasion = time_missing
    if(.not. defined(timeNow))then
      call set_error('Undefined now-time','fu_next_special_time')
      return
    endif
    if(iDirection /= forwards .and. iDirection /= backwards)then
      call msg('Strange time direction:',iDirection)
      call set_error('Strange time direction','fu_next_special_time')
      return
    endif
    if(timeQuantum_%ii > 0)then
      timeQuantum = timeQuantum_
    else
      timeQuantum = timeQuantum_ * (-1)
    endif
    !
    ! All start- special occasions are related to beginning of the current day.
    ! All end- special occasions are trickier because they have to be made exact.
    ! Indeed, we cannot allow the model to overpass this next occasion and show up
    ! at e.g. the next day - the output file name, if time-dependent, will be wrong
    !
    timeOccasion = timeNow
    !
    ! Now let's add / subtract necessary number of days
    !
    select case(iOccasionType)
      case (iStartOfDay)

        timeOccasion%hour = 0    ! the start of the current day
        timeOccasion%minute = 0
        timeOccasion%sec = 0.
        step = one_day
        if(iDirection == forwards)then
          timeOccasion = timeOccasion + one_day   ! start of the next day
        else
          if(timeOccasion == timeNow) timeOccasion = timeOccasion - one_day  ! if already there, jump back
        endif

      case (iStartOfWeek)   ! Week starts on Monday at 00:00
        timeOccasion%hour = 0    ! the very start of the current day
        timeOccasion%minute = 0
        timeOccasion%sec = 0.
        step = one_day * 7
        if(iDirection == forwards)then
          !
          ! If forward in time - get to the start of the next week
          !
          timeOccasion = timeOccasion + (one_day * (8 - fu_weekday(timeNow)))
        else
          !
          ! If backward in time - depending on the current time we should search for start 
          ! of this or previous week
          !
          timeOccasion = timeOccasion - (one_day * (fu_weekday(timeNow) - 1))  ! start of this week
          if(timeOccasion == timeNow) timeOccasion = timeOccasion - step  ! if already there, jump back
        endif

      case (iStartOfMonth)

        timeOccasion%day = 1
        timeOccasion%hour = 0    ! the very start of the current month
        timeOccasion%minute = 0
        timeOccasion%sec = 0.
        if(iDirection == forwards)then
          !
          ! start of the next month
          !
          timeOccasion = timeOccasion + one_day * fu_days_in_month(timeOccasion%month, &
                                                                        & timeOccasion%year)
          step = one_day * fu_days_in_month(timeOccasion%month, timeOccasion%year)
        else
          !
          ! If backward in time - depending on the current time we should search for start 
          ! of this or previous month
          !
          if(timeOccasion == timeNow)then
            if(timeNow%month == 1)then
              timeOccasion = timeOccasion - one_day * 31
              step = one_day * 30
            else
              timeOccasion = timeOccasion - one_day * fu_days_in_month(timeOccasion%month-1, &
                                                                            & timeOccasion%year)
              step = one_day * fu_days_in_month(timeOccasion%month-1, timeOccasion%year)
            endif
          endif
        endif

      case (iStartOfYear)
        timeOccasion%month = 1
        timeOccasion%day = 1
        timeOccasion%hour = 0    ! the very start of the current month
        timeOccasion%minute = 0
        timeOccasion%sec = 0.

        if(iDirection == forwards)then
          timeOccasion%year = timeOccasion%year + 1            ! start of the next year
          step = one_day * (365 + fu_days_in_month(2, timeOccasion%year) - 28)
        else
          !
          ! Depending on the current time we should search for start of this or previous year
          !
          if(timeOccasion == timeNow) timeOccasion%year = timeOccasion%year - 1
          step = one_day * (365 + fu_days_in_month(2, timeOccasion%year-1) - 28)
        endif

      case (iEndOfDay)

        timeOccasion%hour = 23    ! the very end of the current day
        timeOccasion%minute = 59
        timeOccasion%sec = 59.9
        step = one_day
        if(iDirection == forwards)then
          !
          ! If forward in time - get to the end of this or next day
          !
          if(timeOccasion == timeNow) timeOccasion = timeOccasion + one_day
        else
          !
          ! If backward in time - search for end of previous day
          !
          timeOccasion = timeOccasion - one_day
        endif

      case (iEndOfWeek)

        timeOccasion%hour = 23    ! the very end of the current day
        timeOccasion%minute = 59
        timeOccasion%sec = 59.9
        step = one_day * 7
        if(iDirection == forwards)then
          !
          ! If forward in time - get to the end of this week or the next one if now == end of this week
          !
          timeOccasion = timeOccasion + one_day * (7-fu_weekday(timeNow))
          if(timeOccasion == timeNow) timeOccasion = timeOccasion + step
        else
          !
          ! If backward in time - always end of previous week
          !
          timeOccasion = timeOccasion - one_day * fu_weekday(timeNow)  ! end of the previous week
        endif

      case (iEndOfMonth)
        timeOccasion%day = fu_days_in_month(timeNow%month, timeNow%year)
        timeOccasion%hour = 23          ! the very end of the current month
        timeOccasion%minute = 59
        timeOccasion%sec = 59.9

        if(iDirection == forwards)then
          !
          ! If forward in time - get to the end of this or next month
          !
          if(timeOccasion == timeNow) then
            if(timeNow%month == 12)then
              timeOccasion%day = 31
              timeOccasion%month = 1
              timeOccasion%year = timeOccasion%year + 1
              step = one_day * fu_days_in_month(2, timeOccasion%year)
            else
              timeOccasion = timeOccasion + (one_day * fu_days_in_month(timeNow%month+1, &
                                                                             & timeNow%year))
              if(timeOccasion%month == 12)then
                step = one_day * 31
              else
                step = one_day * fu_days_in_month(timeOccasion%month+1, timeOccasion%year)
              endif
            endif
          endif
        else
          !
          ! If backward in time - search for end of previous month
          !
          if(timeNow%month == 1)then
            timeOccasion%day = 31
            timeOccasion%month = 12
            timeOccasion%year = timeOccasion%year - 1
            step = one_day * 31
          else
            timeOccasion = timeOccasion - (one_day * fu_days_in_month(timeNow%month - 1, &
                                                                           & timeNow%year))
            if(timeOccasion%month == 1)then
              step = one_day * 31
            else
              step = one_day * fu_days_in_month(timeOccasion%month - 1, timeOccasion%year)
            endif
          endif
        endif

      case (iEndOfYear)

        timeOccasion%month = 12
        timeOccasion%day = 31
        timeOccasion%hour = 23    ! the very end of the current year
        timeOccasion%minute = 59
        timeOccasion%sec = 59.9
        if(iDirection == forwards)then
          !
          ! If forward in time - get to start of the next year
          !
          if(timeOccasion == timeNow) timeOccasion%year = timeOccasion%year + 1
          step = one_day * (365 + fu_days_in_month(2, timeOccasion%year + 1) - 28)
        else
          !
          ! If backward in time - always search for end of previous year
          !
          timeOccasion%year = timeOccasion%year - 1
          step = one_day * (365 + fu_days_in_month(2, timeOccasion%year) - 28)
        endif

      case default
        call msg('Strange special occasion type:',iOccasionType)
        call set_error('Strange special occasion type','fu_next_special_time')
        return
    end select
    !
    ! Having intermediate variables decided, make the final answer
    !
    if(timeQuantum%ii < 0.1)return ! zero timeQuantum requires exact timeOccasion
    !
    ! If timeQuantum is not zero interval, the distance from timeNow to this new occasion
    ! must be an integer number of quanta. Correspondence to the occasion is then up to
    ! this quantum of accuracy
    !
    if(iDirection == forwards)then
      !
      ! Forwards in time
      !
      nQuanta = nint((timeOccasion - timeNow) / timeQuantum)
      do while(nQuanta == 0)
        timeOccasion = timeOccasion + step
        nQuanta = nint((timeOccasion - timeNow) / timeQuantum)
      end do
      timeOccasion = timeNow + timeQuantum * nQuanta
    else
      !
      ! Backwards in time
      !
      nQuanta = nint((timeNow - timeOccasion) / timeQuantum)
      do while(nQuanta == 0)
        timeOccasion = timeOccasion - step
        nQuanta = nint((timeNow - timeOccasion) / timeQuantum)
      end do
      timeOccasion = timeNow - timeQuantum * nQuanta
    endif

  end function fu_next_special_time


  !********************************************************************************


  INTEGER FUNCTION fu_closest_time_from_times(time, times, direction, ifSilent) result(closest)

    ! Description:
    ! Finds the time that is the closest from the given time, in
    ! the given direction (backwards or forwards). Result is
    ! returned as index to the original time list
    !
    ! If not defined times in given direction are found, then error
    ! is set.
    !
    ! If direction has value single_time, then only the exact time
    ! is searched and its index returned. If it is not found, an error
    ! is set.
    !
    ! Search is stopped when first undefined time is found in the list
    !
    ! If some of the daily_times happens to have same time as given time,
    ! then it is always set to be the closest in all directions.
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in):
    TYPE(silja_time), INTENT(in) :: time
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: times
    INTEGER, INTENT(in) :: direction ! in time
    logical, intent(in), optional :: ifSilent  ! suppresses the error if nothing found. BE CAREFUL !!

    ! Local declarations:
    INTEGER :: i
    REAL :: smallest_so_far, diff 
    LOGICAL :: nextday
   ! nextday = .false.
    !----------------------------------------
    !
    ! 1. Check
    !    -----
    closest = int_missing

    IF (.NOT.defined(time)) THEN
      CALL set_error('undefined time given','fu_closest_time_from_times')
      RETURN
    END IF

    IF (.NOT.defined(times(1))) THEN
      CALL set_error('undefined times given','fu_closest_time_from_times')
      RETURN
    END IF

    closest = 0
    smallest_so_far = 1.0e16

    SELECT CASE (direction)

      !----------------------------------------
      !
      ! 1. Search forwards
      !
      CASE (forwards)

        DO i = 1, SIZE(times)
          IF (.NOT.defined(times(i))) EXIT
          IF (times(i) >= time) THEN
            diff = ABS(fu_sec(time - times(i)))
            IF (diff < smallest_so_far) THEN
              smallest_so_far = diff
              closest = i
!              PRINT*,'closest', closest
            END IF
          END IF
        END DO

        IF (closest == 0) THEN
          if(present(ifSilent))then
            if(ifSilent)return
          endif
          CALL report(time)
          DO i = 1, SIZE(times)
            call msg('i=',i)
            CALL report(times(i))
          END DO
          CALL set_error('all times in past from given time','fu_closest_time_from_times')
        END IF

      !----------------------------------------
      !
      ! 2. Search backwards
      !
      CASE (backwards)

        DO i = 1, SIZE(times)
          IF (.NOT.defined(times(i))) EXIT
          IF (times(i) <= time) THEN
            diff = ABS(fu_sec(time - times(i)))
            IF (diff < smallest_so_far) THEN
              smallest_so_far = diff
              closest = i
            END IF
          END IF
        END DO

        IF (closest == 0) THEN
          if(present(ifSilent))then
            if(ifSilent)return
          endif
          CALL report(time)
          DO i = 1, SIZE(times)
            call msg('i=',i)
            CALL report(times(i))
          END DO
          CALL set_error('all times in future from given time','fu_closest_time_from_times')
        END IF

      !----------------------------------------
      !
      ! 3. Search from both directions
      !
      CASE (back_and_forwards)
        DO i = 1, SIZE(times)
          IF (.NOT.defined(times(i))) EXIT
          diff = ABS(fu_sec(time - times(i)))
          IF (diff < smallest_so_far) THEN
            smallest_so_far = diff
            closest = i
          END IF
        END DO

        IF (closest == 0) THEN
          if(present(ifSilent))then
            if(ifSilent)return
          endif
          CALL report(time)
          DO i = 1, SIZE(times)
            call msg('i=',i)
            CALL report(times(i))
          END DO
          CALL set_error('no times in future or past from given time','fu_closest_time_from_times')
        END IF

      !----------------------------------------
      !
      ! 4. Search only for given time.
      !
      CASE (single_time)
        DO i = 1, SIZE(times)
          IF (.NOT.defined(times(i))) EXIT
          IF (times(i) == time) THEN
            closest = i
            EXIT
          END IF
        END DO

        IF (closest == 0) THEN
          if(present(ifSilent))then
            if(ifSilent)return
          endif
          call msg("looking for time: " +fu_str(time))
          DO i = 1, SIZE(times)
            IF (.NOT.defined(times(i))) EXIT
            call msg("Vector of times: "+fu_str(times(i))+", position:", i)
          END DO
          CALL set_error('time not found in time-vector','fu_closest_time_from_times')
        END IF

      CASE default
        CALL set_error('direction in time wrong','fu_closest_time_from_times')
        RETURN

    END SELECT

  END FUNCTION fu_closest_time_from_times


  ! ***************************************************************



  function fu_start_of_day_utc(time) result(newtime)
   !
   !Replacement for 
   !  fu_set_time_as_daily_time(ordered_times(iTime), midnight_utc)
   ! to get rid of daily_times
   !
   ! Language: ANSI Fortran 90
   !
   ! Author: Roux, FMI

    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_time) :: newtime
    !
    ! Imported parameters with intent(in):
    TYPE(silja_time), INTENT(in) :: time
    !
    ! Local declarations:
    INTEGER :: year, month, day, hour, minute
    REAL :: sec

    IF(.NOT.defined(time)) THEN
      CALL set_error('undefined time given' ,'fu_start_of_day_utc')
      RETURN
    END IF


    ! Components of original time:
    CALL get_time_components_utc(time, year, month, day, hour, minute, sec)

    newtime = fu_set_time_utc(year, month, day, 0, 0, 0.)


  end function fu_start_of_day_utc

  ! ***************************************************************


  real function fu_solar_zenith_angle_cos_time(lon_deg, lat_deg, now)
     ! Uses solar setup. Call SolarSetup before use to get decl terms and eqt_H
!    real, intent(in) :: now
    real, intent(in) :: lat_deg
    real, intent(in) :: lon_deg
    type(silja_time), intent(in) :: now
!    real, intent(out) :: Z                  ! Zenith Angle (degrees)
!    real, intent(out) :: CosZ               ! cos(Z)

    real :: rlt,lzgmt,zpt,lbgmt,thour

        thour = fu_hour(now) + fu_min(now)/60.
        rlt = lat_deg * degrees_to_radians
        lbgmt = 12.0 - eqt_h - lon_deg *24.0/360.0
        lzgmt = 15.0*(thour - lbgmt)
        zpt = lzgmt * degrees_to_radians
        fu_solar_zenith_angle_cos_time = sin(rlt)*sinrdecl + cos(rlt)*cosrdecl*cos(zpt)
        !Check for numerics
        fu_solar_zenith_angle_cos_time = min( 1.0, fu_solar_zenith_angle_cos_time)
        fu_solar_zenith_angle_cos_time = max(-1.0, fu_solar_zenith_angle_cos_time)

 !        Z = acos(CosZ)*degrees_to_radians

  end function fu_solar_zenith_angle_cos_time

!  real function fu_solar_zenith_angle_cos_time(lon_deg, lat_deg, now)
!    !
!    ! Computes the solar zenith angle from basic parameters and LOCAL time
!    !
!    implicit none
!    
!    ! Imported parameters
!    real, intent(in) :: lon_deg, lat_deg
!    type(silja_time), intent(in) :: now
!
!    ! Local variables
!    real :: lat_rad, d, Eqt, Tsolar, w
!    integer :: julian_day
!
!    julian_day = fu_julian_date(now)
!    lat_rad = lat_deg * degrees_to_radians
!
!    !
!    ! Declination of the sun
!    !
!    d = 23.45 * pi / 180. * sin(2. * pi * (284. + julian_day) / 365.)
!
!    !
!    ! Equation of time in minutes
!    ! http://solardat.uoregon.edu/SolarRadiationBasics.html
!    !
!    if(julian_day < 106)then
!      Eqt = -14.2 * sin(pi * (julian_day + 7.) / 111.)
!    elseif(julian_day < 166)then
!      Eqt = 4.0 * sin(pi * (julian_day - 106.) / 59.)
!    elseif(julian_day < 246)then
!      Eqt = -6.5 * sin(pi * (julian_day - 166.) / 80.)
!    else
!      Eqt = 16.4 * sin(pi * (julian_day - 247.) / 113.)
!    endif
!
!    !
!    ! Solar time in hours. Longsm -s the longitude of the "standard meridian" for the given time zone,
!    ! while longLocal is the actual longitude needed. The differrence is then less than 1 hour.
!    ! If we count from Greenwich, Longsm=0, of course
!    !
!!    Tsolar = Tlocal + Eqt / 60. + (Longsm - Longlocal) / 15
!    Tsolar = now%hour + (now%minute+Eqt) / 60. + now%sec / 3600. + lon_deg/15.
!    !
!    ! Hour angle is:
!    !
!    w = pi * (12. - Tsolar) / 12.  ! in radians
!    !
!    ! Cosine of zenith angle
!    !
!    fu_solar_zenith_angle_cos_time = sin(lat_rad) * sin(d) + cos(lat_rad) * cos(d) * cos(w)
!
!  end function fu_solar_zenith_angle_cos_time
!
!  !***********************************************************************

!<===========================================================================
  subroutine SolarSetup(now)
    !
    ! Calculates solar declination and solar-constant related quantities
    ! and sets corresponding private variables
    ! Original code from EMEP model Radiation_ml module
    ! 

    ! Sets up decelention and related terms, as well as Ashrae coefficients
    ! Should be called before other routines.

!    integer, intent(in) :: year,month,day
!    real, intent(in) :: hour
    type(silja_time), intent(in) :: now

    real :: d,ml,rml,w,wr,ec,epsi,yt,pepsi,cww
    real :: sw,ssw, eyt, feqt1, feqt2, feqt3, feqt4, feqt5, &
            feqt6, feqt7, feqt,ra,reqt
    real :: dayinc
    integer :: i

    logical, parameter :: MY_DEBUG = .false.


!* count days from dec.31,1973 
!    d = julian_date(year,month,day) - julian_date(1973,12,31) + 1
!    d = d + hour/24.0

    d = fu_days(now - ref_time_30121974)


!* calc geom mean longitude

    ml = 279.2801988 + 0.9856473354*d + 2.267e-13*d*d
    rml = ml*degrees_to_radians

!* calc equation of time in sec
!*  w = mean long of perigee
!*  e = eccentricity
!*  epsi = mean obliquity of ecliptic

    w = 282.4932328 + 4.70684e-5*d + 3.39e-13*d*d
    wr = w*degrees_to_radians
    ec = 1.6720041e-2 - 1.1444e-9*d - 9.4e-17*d*d
    epsi = 23.44266511 - 3.5626e-7*d - 1.23e-15*d*d
    pepsi = epsi*degrees_to_radians
    yt = tan(pepsi/2.0)
    yt = yt*yt
    cww = cos(wr)
    sw = sin(wr)
    ssw = sin(2.0*wr)
    eyt = 2.0*ec*yt
    feqt1 = sin(rml)*(-eyt*cww - 2.0*ec*cww)
    feqt2 = cos(rml)*(2.0*ec*sw - eyt*sw)
    feqt3 = sin(2.0*rml)*(yt - (5.0*ec*ec/4.0)*(cww*cww-sw*sw))
    feqt4 = cos(2.0*rml)*(5.0*ec*ec*ssw/4.0)
    feqt4 = cos(2.0*rml)*(5.0*ec*ec*ssw/4.0)
    feqt5 = sin(3.0*rml)*(eyt*cww)
    feqt6 = cos(3.0*rml)*(-eyt*sw)
    feqt7 = -sin(4.0*rml)*(0.5*yt*yt)
    feqt = feqt1 + feqt2 + feqt3 + feqt4 + feqt5 + feqt6 + feqt7

    eqtime = feqt*13751.0

!* equation of time in hrs:

    eqt_h = eqtime/3600.0

!* convert eq of time from sec to deg

    reqt = eqtime/240.0

!* calc right ascension in rads

    ra = ml - reqt

!* calc declination in rads, deg

    tan_decl = 0.43360*sin( ra * degrees_to_radians )

    rdecl    = atan(tan_decl)
    sinrdecl = sin(rdecl)
    cosrdecl = cos(rdecl)

 !-----------------------------------------------------------------

     ! Find coefficients for Iqbal/Ashrae algorith, used in ClearSkyRadn routine
     ! first, perform the table look up

      d = fu_julian_date(now)!day_of_year(year,month,day)
      do i = 1, 14
        if (d <=   ASHRAE_REV(i)%nday ) exit
      end  do

      if ( ASHRAE_REV(i)%nday == 1) then
        Ashrae%a = ASHRAE_REV(1)%a
        Ashrae%b = ASHRAE_REV(1)%b
        Ashrae%c = ASHRAE_REV(1)%c
      else
        dayinc = real( d-ASHRAE_REV(i-1)%nday ) / &
                 real( ASHRAE_REV(i)%nday-ASHRAE_REV(i-1)%nday )
        Ashrae%a = ASHRAE_REV(i-1)%a + &
             ( ASHRAE_REV(i)%a - ASHRAE_REV(i-1)%a )*dayinc
        Ashrae%b = ASHRAE_REV(i-1)%b + &
             ( ASHRAE_REV(i)%b - ASHRAE_REV(i-1)%b )*dayinc
        Ashrae%c = ASHRAE_REV(i-1)%c + &
             ( ASHRAE_REV(i)%c - ASHRAE_REV(i-1)%c )*dayinc
      end if

 end subroutine SolarSetup








  function fu_daylength(lat, now) !returns hours
    implicit none
    real, intent(in) :: lat  ! degrees
    type(silja_time), intent(in) :: now
    type(silja_interval) :: fu_daylength
    ! real, intent(in) :: julday ! should be integer?
    integer :: julday
    real :: ftmp
    
    julday = fu_julian_date(now)
    ! http://mathforum.org/library/drmath/view/56478.html
    !   D = daylength
    !   L = latitude
    !   J = day of the year
    !
    !   P = asin[.39795*cos(.2163108 + 2*atan{.9671396*tan[.00860(J-186)]})]
    !
    !                          _                                         _
    !                         / sin(0.8333*pi/180) + sin(L*pi/180)*sin(P) \
    !   D = 24 - (24/pi)*acos{  -----------------------------------------  }
    !                         \_          cos(L*pi/180)*cos(P)           _/
    !
    !Roux:  added check for saturation
    !
    ftmp = .39795*cos(.2163108 + 2.*(.00860*(julday-186.)))
    ftmp = (sin(0.8333/57.6) + sin(lat/57.6)*sin(ftmp))/(cos(lat/57.6)*cos(ftmp))
    ftmp = min(max(ftmp, -1.), 1.)
        
    fu_daylength = fu_set_interval_sec(24.*( 1. -  acos(ftmp) / 3.1415927)*3600.)
      
  end function fu_daylength
  
  
  !***********************************************************************

  function fu_sunrise_utc(lon, lat, now) 
    !
    ! Computes the sunrise time in UTC
    ! http://williams.best.vwh.net/sunrise_sunset_algorithm.htm
    !
    !zenith: offical 90deg50'; civil 96deg; nautical 102deg; astronomical 108deg
    !
    implicit none
    
    real, intent(in) :: lon, lat
    type(silja_time), intent(in) :: now
    type(silja_time) :: fu_sunrise_utc
    real :: T, M, L, RA, dec, cosH
      
   T = real(fu_julian_date(now)) + ((6. - lon / 15.) / 24.)  !(Jd+((18-lon/15)/24) for sunset)
   
   M = (0.9856 * T) - 3.289  !Sun's mean anomaly
   
   L = M + (1.916 * sin(M * degrees_to_radians)) + &
         & (0.020 * sin(2 * M * degrees_to_radians)) + 282.634 !Sun's true longitude
   if(L > 360.) L = L - 360.
   if(L < 0.) L = L + 360. 

   RA = atan(0.91764 * tan(L * degrees_to_radians)) * radians_to_degrees  ! Sun's right ascension
    if(RA > 360.) RA = RA - 360.
   if(RA < 0.) RA = RA + 360.
   RA = RA + ((floor( L/90.)) * 90. - (floor(RA/90.)) * 90.) !same quadrant as L

   dec = asin(0.39782 * sin(L*degrees_to_radians)) * radians_to_degrees !Sun's declination

   cosH = (cos(90.833 * degrees_to_radians) - &
          & (sin(dec * degrees_to_radians) * sin(lat * degrees_to_radians))) / &
          & (cos(dec * degrees_to_radians) * cos(lat * degrees_to_radians))  !Sun's local hour angle
      
   if (cosH >  1)then ! the sun never rises 
       fu_sunrise_utc = fu_set_time_utc(fu_year(now), fu_mon(now), fu_day(now), 0, 0, 0.0) + one_day
   elseif (cosH < -1)then !the sun never sets 
        fu_sunrise_utc = fu_set_time_utc(fu_year(now), fu_mon(now), fu_day(now), 0, 0, 0.0)
    else
        T = (360. - acos(cosH) * radians_to_degrees + RA - lon) / 15 - (0.06571 * T) - 6.622  !UTC time of rising/setting
        if(T > 24.) T = T - 24.
        if(T < 0.) T = T + 24.   
        
        fu_sunrise_utc = fu_set_time_utc(fu_year(now), fu_mon(now), fu_day(now), 0, 0, 0.0) + &
                   & fu_set_interval_sec(T * 3600.0)
    endif
    
  end function fu_sunrise_utc


  ! ***************************************************************




  SUBROUTINE report_time(time)
  
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silja_time), INTENT(in) :: time

    character(len=clen) :: timestr

    timestr =  fu_str(time)

    call msg(trim(timestr))
  
  END SUBROUTINE report_time


  ! ***************************************************************


  SUBROUTINE report_times(times)
  
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: times
    
    INTEGER :: i

    DO i = 1, SIZE(times)
      IF (.NOT.defined(times(i))) CYCLE
      call report(times(i))
    END DO

  END SUBROUTINE report_times


  ! ***************************************************************


  SUBROUTINE report_interval(interval)
  
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silja_interval), INTENT(in) :: interval

    character(len=clen) :: timestr

    timestr = fu_interval_string(interval)

    call msg(timestr)
  
  END SUBROUTINE report_interval



 ! ***************************************************************


 function timezone_index_by_name(myTZName)
         CHARACTER (LEN=*), INTENT(in) :: myTZName
         integer :: timezone_index_by_name
         integer :: i
         logical :: Found
         Found = .False.
        ! if ( trim(adjustl(fu_str_u_case(myTZname))) == "UTC") then
        !        timezone_index_by_name = 1
        !        Found = .True.
        ! else
                do i = 1, number_of_time_zones
                     if (trim(adjustl(TZ_name(i))) == trim(adjustl(myTZname)))then
                        timezone_index_by_name = i
                        Found = .True.
                        exit
                      end if
                end do
        !end if
        if (.Not. Found) then
               call set_error("Can't find index for timezone "+ myTZName, &
                        &'timezone_index_by_name')
        endif

 end function timezone_index_by_name

 !*****************************************************************************

 subroutine init_timezones(nlStdSetup) ! initializes timezones form a file

    type(Tsilam_namelist), pointer, intent(in) :: nlStdSetup

    integer :: iTmp,iZone,iStat
    integer :: nVals
    type(Tsilam_namelist_group), pointer :: nmTblNlGrp
    type(Tsilam_namelist), pointer :: nlPtr
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: ptrItems
    character(len=1024) :: chContent
    integer :: iInd, iSummerOff, iWinterOff
    character(len=1024) :: chTZName, chTZWName
    character(len=10) :: chTmp

    logical :: if_have_country_codes

    if_have_country_codes = .False.

    select case (fu_str_u_case(fu_content(nlStdSetup,'timezone_method')) )
   !        case('MEAN')  
   !             call msg("init_timezones: No timezone update methos: using mean")
   !             timezone_method = mean_timezone_method
        
        case('JANUARY')
                call msg("init_timezones: Using  January time")
                timezone_method = winter_timezone_method
        case('JULY')
                call msg("init_timezones: Using  July time")
                timezone_method = summer_timezone_method
        case('','LASTSUNDAY')  !DEFAULT
                call msg("init_timezones: Using  last sunday of March/October to swith timezone")
                timezone_method = lastsunday_timezone_method
        case('SYSTEM')
                call msg("init_timezones: System timezones not implemented yet")
                call msg("init_timezones: Falling back to lastsunday method")
                timezone_method = lastsunday_timezone_method
        case default
                call set_error("Strange timezone method" + &
                        &fu_content(nlStdSetup,'timezone_method'), 'init_timezones' )
    end select


    ! No timezone file -- only UTC zone 
    if (fu_content(nlStdSetup,'timezone_list_fnm') /= '') then


            ! Read the timezone table file. 
            iTmp = fu_next_free_unit()
            open(file=fu_content(nlStdSetup,'timezone_list_fnm'), &
                   & unit=iTmp, action='read', status='old', iostat=iStat)
            if(iStat /= 0)then
              call set_error('Failed to open timezone list:' + &
                           & fu_content(nlStdSetup,'timezone_list_fnm'), 'timezones_init')
              return
            endif

            nmTblNlGrp => fu_read_namelist_group(iTmp, .false., 'END_TIMEZONE_INDEX_0')
            close(iTmp)
            IF (error) RETURN


            ! Process timezone list
            if(.not. associated(nmTblNlGrp))then
              call msg_warning('TIMEZONE INDEX namelist group is not associated','timezones_init')
            endif

            if (fu_nbr_of_namelists(nmTblNlGrp) /= 1) then
                call set_error("Only one timezone index is allowed",'timezones_init')
                return
            endif

            nlPtr => fu_namelist(nmTblNlGrp, 1)
            nullify(ptrItems)
            call get_items(nlPtr, 'zone', ptrItems, nVals)
            if(nVals < 1)then
                 call set_error('No zones in namelist','timezones_init')
                 return
            endif
            number_of_time_zones = nVals

            if_have_country_codes = &
             & (fu_str_u_case(fu_content(nlPtr,'have_country_code')) == 'YES')
           


      else
            call msg("No timezone list: using only 24 UTC-related zones")
            nVals = 0
            number_of_time_zones = 24
      endif  !Timezone file


      allocate(TZ_offset_summer(number_of_time_zones), &
              & TZ_offset_winter(number_of_time_zones), &
              & TZ_name(number_of_time_zones), &
              & TZ_code(number_of_time_zones), &
              & TZ_windows_name(number_of_time_zones), &
              & TZ_offset(number_of_time_zones),stat = iStat)
      if(iStat /= 0)then
        call set_error('Failed to allocate TZ arrays :', 'timezones_init')
        return
      endif

     if (nVals == 0) then
            ! No timezone file -- fill in the info to enable at least solar times
            do iZone = 1,24
              solar_TZ_index(iZone) = iZone
              iSummerOff = modulo(iZone -13,24) - 12  ! UTC first!
              TZ_offset_summer(iZone) = iSummerOff*3600
              TZ_offset_winter(iZone) = TZ_offset_summer(iZone)
              TZ_offset(iZone) =  TZ_offset_summer(iZone)
              write (TZ_name(iZone),'(A,I2)') "Etc/GMT", -iSummerOff ! Sic!!!
              write (TZ_windows_name(iZone),'(A,I2)') "UTC", iSummerOff 
            enddo
            TZ_name(1) = "UTC"
            TZ_windows_name(1) = "UTC"

    else
       do iZone = 1, nVals
            chContent=fu_content(ptrItems(iZone))

            if (if_have_country_codes) then
                read(unit=chContent,fmt=*,iostat=iStat) iInd, chTmp, iWinterOff, iSummerOff, chTZName
            else
               read(unit=chContent,fmt=*,iostat=iStat) iInd,  iWinterOff, iSummerOff, chTZName
               chTmp = "ALL" !Placeholder for timezone country code
            endif

            if (iZone /= iInd) then
              call msg( trim(adjustl(fu_content(ptrItems(iZone)))) // ' Expected index: ', iZone)
              call set_error('Wrong timezone index :', 'timezones_init')
              return
            endif
            TZ_offset_summer(iZone) = iSummerOff
            TZ_offset_winter(iZone) = iWinterOff

            TZ_offset(iZone) =   (iSummerOff+iWinterOff) / 2 ! fill with default value: failsafe  
            
            TZ_code (iZone) = chTmp
            ! Fortran can't read beyond salsh somehow...
            iTmp = index(chContent,trim(chTZName))
            chTZName=chContent(iTmp:)
            iTmp = index(chTZname,' ')
            TZ_name(iZone) = chTZname(1:(iTmp-1)) !before space
            TZ_windows_name(iZone) = trim(adjustl(chTZname(iTmp:))) ! everything after space
            
!            write (*,*) "TZ Name: ", trim(TZ_name(iZone))
!            write (*,*) "TZ W Name: ",trim(TZ_windows_name(iZone))    
!            call msg( "Index", iInd)
!            call msg( "Summer", iSummerOff)
!            call msg( "Summer", TZ_offset_summer(iZone))
!            call msg( "Winter", iWinterOff)
!            call msg( "Winter", TZ_offset_summer(iZone))
!            call msg( "Default", (iZone))

!            call msg( "")

       enddo

       !Find proper zones for Solar times
       do iTmp = 1,24
            iSummerOff = 3600 * (modulo(iTmp -13,24) - 12)  ! UTC first!
            do iZone = 1, nVals
              if (TZ_offset_summer(iZone) == iSummerOff .and. &
                 & TZ_offset_winter(iZone) == iSummerOff) then
                 solar_TZ_index(iTmp) = iZone
                 exit
              endif
            enddo
            if (iZone > nVals) then
              call msg("Failed Solar time offset", iSummerOff)
              call set_error("Could not find proper timezone for Solar time", &
                        & "init_timezones")
            endif
       enddo
       deallocate(ptrItems)
     endif

 end subroutine init_timezones

 !*****************************************************************************************

 subroutine update_timezone_offsets(time)
         !
         ! updates global timezones 
         ! They get updated at 00 UTC over whole globe
         !
    implicit none
    TYPE(silja_time),  INTENT(in) :: time
    integer :: y,m,d,hh,mm, dow
    real :: s
    integer :: iZone
    logical :: Winter


    select case (timezone_method)
!        case(mean_timezone_method)
!                call msg("No timezone update method: winter/summer time")
        case(winter_timezone_method)
                TZ_offset(1:number_of_time_zones) =  TZ_offset_Winter(1:number_of_time_zones)
                call msg("Using fixed January local times")
        case(summer_timezone_method)
                TZ_offset(1:number_of_time_zones) =  TZ_offset_Winter(1:number_of_time_zones)
                call msg("Using  fixed July local times")
        case(lastsunday_timezone_method)
                ! Switch on last Sunday of March and October
                CALL get_time_components_utc(time,y,m,d,hh,mm,s)
                Winter = .False.
                if (m <= 3 .or. m > 10) then! March and October
                        Winter = .True.
                endif
                if ((m == 3 .or. m == 10) .and. (d >= 25)) then ! need to check weekday 
                        dow=fu_weekday(time)
                        !Sunday today or later then end of month
                        ! new time should be in force
                        if (dow == 7 .or. (7-dow) > (31-d) ) then
                                Winter = .not. Winter 
                        endif
                endif
                if (Winter) then
                  TZ_offset(1:number_of_time_zones) =  TZ_offset_Winter(1:number_of_time_zones)
                  call msg("Applying wintertime for zones for " // trim(fu_str(time)))
                else
                  TZ_offset(1:number_of_time_zones) =  TZ_offset_summer(1:number_of_time_zones)
                  call msg("Applying summertime for zones for " // trim(fu_str(time)))
                endif
        case default
                call set_error("Strange timezone method " + &
                        & fu_str(timezone_method), &
                        & "update_timezone_offsets")
    end select
    

 end subroutine update_timezone_offsets



  ! ***************************************************************

  ! ***************************************************************
  ! 
  !
  !                 TESTS
  !
  !
  ! ***************************************************************
  ! ***************************************************************


  SUBROUTINE time_tests()
    !
    ! Description:
    ! Test stuff for times.
    !
    ! Author: Mika Salonoja, FMI
    !
    ! Modules used:
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters:
    !
    ! Local declarations:

    TYPE(silja_time) :: time1, time2, time
    TYPE(silja_interval) :: diff, interval, interval1, interval2

    TYPE(silja_time) :: summertime_change
    integer :: it, itmp
    real(kind=8) :: fTmp8

    !----------------------------------------
    !
    ! 1. 
    !
    ! This  is the demonstartion of a bug
    print *, "Conversion to seconds and back"
    time = fu_set_time_utc(1969, 1, 1, 0, 0, 0.0)
    fTmp8  = silja_time_to_real8(time)
    time1 = real8_to_silja_time(fTmp8)
    ! call ooops("")
    print *, "Time   "//fu_str(time), silja_time_to_real8(time)
    print *, "Time1  "//fu_str(time1), silja_time_to_real8(time1)


!!!$
!!!$
!!!$    PRINT*, ' Time tests:'
!!!$    diff = fu_set_interval_min(20)    
!!!$
!!!$    time1 = fu_set_time(1995, 6, 30, 23, 0, 0., fu_set_interval_h(3))
!!!$    time2 = fu_set_time(1995, 8, 1, 0, 0, 0., fu_set_interval_h(3))  
!!!$    
!!!$    PRINT*, 'time1 :        ', fu_time_string_original(time1)
!!!$    PRINT*, 'time1 in utc:  ', fu_str(time1)
!!!$    
!!!$    PRINT*, 'time1 + diff : ', fu_time_string_original(time1 +&
!!!$  & diff)
!!!$    PRINT*, 'time2 :        ', fu_time_string_original(time2)
!!!$    PRINT*, 'time-diff:     ', time1 - time2
!!!$    
!!!$    PRINT*, 'time now :     ',&
!!!$  & fu_time_string_original(fu_wallclock())
!!!$    PRINT*, 'utc-time now : ',&
!!!$  & fu_str(fu_wallclock())
!!!$    
!!!$
!!!$    PRINT*, ' End of time tests:'
!!!$    PRINT *, '****************************************************'
!!!$
!!!$    interval = fu_set_interval_sec(120.)
!!!$
!!!$    time = fu_set_time_utc(1995, 12, 1, 0, 0, 0.)
!!!$    PRINT*, 'time in utc:  ', fu_str(time)
!!!$    !PRINT*, 'time in utc:  ', time
!!!$    PRINT*,  interval
!!!$    time = time + interval
!!!$    PRINT*, 'time in utc:  ', fu_str(time)
!!!$    !PRINT*, 'time in utc:  ', time
!!!$    
!!!$    PRINT *, '****************************************************'
!!!$    
!!!$    time1 = fu_set_time_utc(1995, 12, 1, 0, 0, 0.)
!!!$    time2 = fu_set_time_utc(1995, 12, 1, 0, 1, 0.)
!!!$    
!!!$    interval = time2 - time1
!!!$!    PRINT*,  interval
!!!$
!!!$    PRINT *, fu_time_string_original(fu_wallclock())
!!!$    PRINT *, fu_str(fu_wallclock())
!!!$    PRINT *, fu_str(fu_wallclock() + one_hour)

    PRINT *, fu_str(fu_wallclock())
    PRINT *, fu_computer_time_string()
    PRINT*, fu_real_hours_utc(fu_wallclock())

    PRINT *, fu_julian_date(fu_wallclock())

    time1 = fu_set_time_utc(1997, 2, 28, 0, 0, 0.) 
    time2 = fu_set_time_utc(1997, 3, 1, 0, 0, 0.) 

    PRINT *, fu_str(time1)
    PRINT *, fu_str(time2)

    PRINT *, fu_interval_string(time1 - time2)
    PRINT *, fu_sec(time1 - time2)

    print *, 'Testing Julian dates'
    call test_julian_dates()



    ! This  is the demonstartion of a bug
    print *, "Trying one second accuracy for time + interval"
    time = fu_set_time_utc(2009, 5, 31, 12, 51, 0.0)


    iTmp = 1

    do while (.not. error)
       time1 = time + (one_week*iTmp)
       time2 = time1 + one_second
       interval1 = time1 - time
       interval2 = time2 - time
       if ( abs(fu_sec(((time + interval1) - (time + interval2)) + one_second)) > 0.001 ) then
            print *, "1s accuracy broke at interval:",  fu_str(interval1)
            exit
       endif
       if (fu_sec(interval1)< 0.) then
          print *, "Broke interval counter at: ", fu_str(one_week * (iTmp/2) )
          print *, "Test passed"
          exit
       endif
       print *, interval1%ii,  interval2%ii, iTmp
       iTmp = iTmp * 2
       if (iTmp <0) then
          print *, "loop counter overflow, accuracy not broken"
          print *, "Test passed"
          exit
       endif

    enddo

    

    print *, "Time   "//fu_str(time), silja_time_to_real8(time)
    print *, "Time1  "//fu_str(time1), silja_time_to_real8(time1)
    print *, "Time2  "//fu_str(time2), silja_time_to_real8(time2)
    
    print *, "Time2 - Time1 = ", fu_sec8(time2 - time1)
    print *, "Time1 - Time  = ", fu_sec8(interval1)
    print *, "Time2 - Time  = ",fu_sec8(interval2)


!    print *, fu_str(time)
!    time = real8_to_silja_time(silja_time_to_real8(time))
!    print *, fu_str(time)
!    print *, ""
!    print *, fu_str(time1)
!    time = real8_to_silja_time(silja_time_to_real8(time1))
!    print *, fu_str(time)
!    print *, ""
!    print *, fu_str(time2)
!    time = real8_to_silja_time(silja_time_to_real8(time2))
!    print *, fu_str(time)



  contains
    
    subroutine test_julian_dates()
      implicit none
      integer :: jul_day_int
      real(r8k) :: jul_day_real
      type(silja_time) :: time, time2
      character(len=*), parameter :: sub_name = 'test_julian_dates'
      
      time = fu_set_time_utc(2009, 5, 31, 12, 51, 0.0)
      if (fu_fails(fu_julian_date(time) == 151, 'T1', sub_name)) continue
      if (fu_fails(ceiling(fu_julian_date_real(time)) == 151, 'T2', sub_name)) then 
        print *, fu_julian_date(time), fu_julian_date_real(time)
      end if
      
      time2 = fu_julian_date_to_time(fu_julian_date_real(time), fu_year(time))
      print *, fu_str(time2)
      if (fu_fails(fu_julian_date_to_time(fu_julian_date_real(time), fu_year(time)) == time, 'T3', sub_name)) &
           & continue

    end subroutine test_julian_dates
    

  END SUBROUTINE time_tests



END MODULE silam_times

