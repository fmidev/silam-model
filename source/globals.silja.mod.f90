MODULE globals
  !
  ! Description:
  ! All kinds of global stuff. This
  ! module has to be used by all the modules of silja-program.!
  !
  ! List of contents:
  !
  ! 1. Definition and tools for the type silja_logical
  !
  ! 2. Global parameters, like missing codes and default length of
  ! character strings etc.
  !
  ! 3. Global error-flag and tools for error situations.
  !
  ! 4. Global test-messages -flag and msg_test -subroutine.
  !
  ! 5. Names of physical quantities and the corresponding flag-values
  ! (integer-type).
  ! 
  ! 6. Names of sources of data and the corresponding flag-values
  ! (integer-type). Also tools for checking for the type of a
  ! quantity.
  ! 
  ! 7. Names of interpolation methods and the corresponding flag
  ! -values (integer-type).
  !
  ! 8. Name of host computer (integer-type) and its int2string conversion
  !
  ! Author: Mika Salonoja, FMI
  ! Addition 8 - Mikhail Sofiev, FMI
  !
  ! Language: ANSI standard Fortran 90

  ! Modules used:
  USE natural_constants
  USE max_sizes_and_limits
  USE iso_fortran_env, only : real32, real64, int64, int32
  !$ use OMP_LIB
#ifdef SILAM_MPI
  USE mpi
#endif
#ifdef VS2012
  use ifport
  use ifcore
#endif

  IMPLICIT NONE

  INTEGER, PARAMETER, PUBLIC :: r4k = real32
  INTEGER, PARAMETER, PUBLIC :: r8k = real64

  ! The public functions and subroutines available in this module:
  PUBLIC fu_set_true
  PUBLIC fu_set_false
  PUBLIC fu_true
  PUBLIC fu_false
  PUBLIC fu_undefined
  public fu_fails
  PUBLIC set_error
  PUBLIC set_error_nomsg
  PUBLIC unset_error
  PUBLIC msg_warning
  PUBLIC msg
  PUBLIC msg_test
  PUBLIC fu_connect_strings
  public fu_index_non_delim
  public fu_str
  public ooops
  PUBLIC fu_name
  PUBLIC smpi_is_mpi_version
  public report_vectorstat


  PRIVATE fu_compare_logicals_eq ! interface ==
  private fu_not_true_silja_logical
  private msg_txt
  private msg_int
  private msg_real4
  private msg_real8
  private msg_int_real4
  private msg_int_real8
  private msg_real4_int
  private msg_real8_int
  private msg_int_int
  private msg_real4_real4
  private msg_real8_real8
  private msg_int_array
  private msg_real_array
  private msg_logical_array
  private msg_real8_array
  private msg_int8_array
  private msg_test_txt
  private msg_test_int
  private msg_test_real
  private msg_test_int_real
  private msg_test_real_int
  private msg_test_int_int
  private msg_test_real_real
  private fu_connect_two_strings
  private fu_int2str
  private fu_real2str
  private fu_name_of_integer
  
  interface msg
    module procedure msg_txt
    module procedure msg_int
    module procedure msg_real4
    module procedure msg_real8
    module procedure msg_real4_int
    module procedure msg_real8_int
    module procedure msg_int_real4
    module procedure msg_int_real8
    module procedure msg_int_int
    module procedure msg_real4_real4
    module procedure msg_real8_real8
    module procedure msg_int_array
    module procedure msg_int8_array
    module procedure msg_logical_array
    module procedure msg_real_array
    module procedure msg_real8_array
  end interface

  interface msg_test
    module procedure msg_test_txt
    module procedure msg_test_int
    module procedure msg_test_real
    module procedure msg_test_real_int
    module procedure msg_test_int_real
    module procedure msg_test_int_int
    module procedure msg_test_real_real
  end interface

  interface fu_str
    module procedure fu_int2str
    module procedure fu_real2str
    module procedure fu_real82str
  end interface

  interface fu_name
    module procedure fu_name_of_integer
  end interface
  
  INTERFACE operator(+)
    MODULE PROCEDURE fu_connect_two_strings
  END INTERFACE

  !
  !        SILJA-LOGICAL STUFF
  !
  TYPE silja_logical
    PRIVATE
    INTEGER :: flag = 0
  END TYPE silja_logical

  INTEGER, PARAMETER, PRIVATE :: false_flag = -5
  INTEGER, PARAMETER, PRIVATE :: true_flag = 5

  TYPE(silja_logical), PARAMETER, PUBLIC :: silja_false = silja_logical(false_flag)

  TYPE(silja_logical), PARAMETER, PUBLIC :: silja_true = silja_logical(true_flag)

  TYPE(silja_logical), PARAMETER, PUBLIC :: silja_undefined = silja_logical(0)

  INTERFACE operator(==)
    MODULE PROCEDURE fu_compare_logicals_eq
  END INTERFACE

  interface operator (.not.)
    module procedure fu_not_true_silja_logical
  end interface


  ! The SILJA_LOGICAL -type is defined here. It differs from the
  ! normal logical variable so that it has three possible values:
  ! ture, false or not any value. It is created so that we can check
  ! the value of a logical variable even when it has been declared in
  ! a subroutine, but not set any value yet. It is especially used in
  ! higher-level types to check, if the variable has already been
  ! given a value or initialized correctly. The type is private and
  ! it is handled with functions found in this module. 
  !
  ! Here is an example of usage of silja_logical:
  !
  !  TYPE(silja_logical) :: defined
  !  ...
  !  ...
  !  
  !  IF (fu_true(defined)) THEN
  !    ...
  !
  !  ELSE IF (fu_false(defined)) THEN
  !    ...
  !
  !  ELSE ! not any value yet!
  !    ...
  !
  !  END IF
  !
  !        Default sizes
  !
  INTEGER, PARAMETER, PUBLIC :: clen = 80 ! The typical length of string
  INTEGER, PARAMETER, PUBLIC :: fnlen = 700 ! The typical length of
                                            ! file name with full path
  integer, parameter, public :: substNmLen = 30  ! Length of the substance name
  integer, parameter, public :: unitNmLen = 20  ! Length of the unit name


  !
  !        MISSING CODES
  !
  ! The missing codes of intrinsic types:
  INTEGER, PARAMETER, PUBLIC :: int_missing = -999999
  REAL, PARAMETER, PUBLIC :: real_missing = -999999.e9   !1.0E35
  ! Quiet NAN, double precision.
  REAL(r4k), PARAMETER, PUBLIC :: F_NAN = TRANSFER(2143289344,1.0_r4k)
  REAL(r8k), PARAMETER, PUBLIC :: D_NAN =  TRANSFER(-2251799813685248_int64,1.0_r8k)
  real, PARAMETER, PUBLIC :: F_EPS = EPSILON(real_missing)

  integer(kind=8), parameter, PUBLIC :: MAX_INT32 = 2**30+(2*30-1)
#ifdef DOUBLE_PRECISION
  REAL(r8k), PARAMETER :: CONST_NAN = D_NAN
  INTEGER, PARAMETER :: DEFAULT_REAL_KIND = r8k
#else
  REAL(r4k), PARAMETER ::  MIN_FLOAT = TRANSFER(INT(Z'00800000'),1.0_r4k) !!Minimum float at full precision
  REAL(r4k), PARAMETER :: CONST_NAN = F_NAN
  INTEGER, PARAMETER :: DEFAULT_REAL_KIND = r4k
#endif
  CHARACTER (LEN=4), PARAMETER, PUBLIC :: char_missing = '-999999'

  integer, parameter, public :: accept_all = -888888

  !
  !        TEST MESSAGES
  !
  LOGICAL, SAVE, PUBLIC :: test_messages ! if true, test messages are
                                         ! printed  during the run

!!!!!!!!!#defin DEBUG
  
#ifdef DEBUG
#ifdef DEBUG_MORE
  integer, save, public :: debug_level = 2
#else
  integer, save, public :: debug_level = 1
#endif
#else
  integer, save, public :: debug_level = 0
#endif

  !
  !        ERROR VARIABLES
  !
  INTEGER, PARAMETER, PUBLIC :: error_msglength = 200 
  ! moved to silam_mpi! LOGICAL, SAVE, PUBLIC :: error

  !
  !        LOG FILE UNITS
  !
  integer, public, save :: info_funit = int_missing, run_log_funit = int_missing
  character(len=fnlen), public, save :: run_log_name='',  run_log_tmp_name=''
  character(len=fnlen), public, save :: proc_ID_string = '' !! Set in silam_main, unique for the MPI process

  logical, public, save :: ifPrintDump = .false.

  !
  !   GRIB CODE DEFINITIONS.
  !
  !  IMPORTANT. In HIRLAM grib files the centre is 96, although must be 86
  !  Table version is 1, although must be 2 (list of parameters is from there)
  !  As a result, there is some mess in the code_table definitions,
  !  requiring funny cente_hirlam and table_fmi.
  !  They must NOT exist, but appear for whatever reasons
  !
  ! centre
  INTEGER, PARAMETER, PUBLIC :: centre_fmi = 86
  INTEGER, PARAMETER, PUBLIC :: centre_ecmwf = 98  
  INTEGER, PARAMETER, PUBLIC :: centre_hirlam = 96
  INTEGER, PARAMETER, PUBLIC :: centre_smhi = 82
  INTEGER, PARAMETER, PUBLIC :: centre_moscow = 4
  INTEGER, PARAMETER, PUBLIC :: centre_UERRA = 233
  INTEGER, PARAMETER, PUBLIC :: centre_metcoop =  251  !Harmonie 

  INTEGER, PARAMETER, PUBLIC :: table_harmonie = 253 ! Version
  INTEGER, PARAMETER, PUBLIC :: table_wmo_internat = 2 ! Version of code_table 2
  INTEGER, PARAMETER, PUBLIC :: table_fmi = 1          ! Version of code_table 2
!  INTEGER, PARAMETER, PUBLIC :: table_silam = 5        ! Version of code_table 2

! generatingProcessIdentifier
  INTEGER, PARAMETER, PUBLIC :: model_hirlam  = 1 ! I do not know but...
  INTEGER, PARAMETER, PUBLIC :: model_hirlam_eno = 2 ! I do not know but...
  INTEGER, PARAMETER, PUBLIC :: model_harmonie = 3
  INTEGER, PARAMETER, PUBLIC :: model_meps = 40 !!Harmonie MEPS
  INTEGER, PARAMETER, PUBLIC :: model_meps1 = 0 !!Harmonie MEPS
  INTEGER, PARAMETER, PUBLIC :: model_silam = 1001  ! 
  INTEGER, PARAMETER, PUBLIC :: model_silam_internal = 1002  ! 
  integer, parameter, public :: model_ecmwf = 5  ! 
  INTEGER, PARAMETER, PUBLIC :: model_hirlam_ata = 6 ! I do not know but...
  INTEGER, PARAMETER, PUBLIC :: model_cosmo = 132

  integer, dimension(7), parameter, public :: silam_grib_tables = &
                                   & (/130,131,132,133,134,135,136/) ! Version of code_table 2

  !
  !   EXTERNAL SOURCES OF METEOROLOGICAL DATA
  !
  !INTEGER, PARAMETER, PUBLIC :: atlantic_hirlam = 230106
  !INTEGER, PARAMETER, PUBLIC :: european_hirlam = 230107
  !INTEGER, PARAMETER, PUBLIC :: nordic_hirlam = 230108
  !INTEGER, PARAMETER, PUBLIC :: ecmwf_model = 230109
  !
  !INTEGER, PARAMETER, PUBLIC :: synoptic_observations = 230150
  !INTEGER, PARAMETER, PUBLIC :: powerplant_mast = 230151

  !
  ! Allowed formats of input data
  !
  integer, parameter, public :: grib_file_flag = 20001
  integer, parameter, public :: ascii_file_flag = 20002
  integer, parameter, public :: grads_file_flag = 20003
  integer, parameter, public :: test_field_value_flag = 20004
  integer, parameter, public :: netcdf_file_flag = 20005

  !
  !   INTERNAL SOURCES OF MODEL OUTPUT DATA
  !
  INTEGER, PARAMETER, PUBLIC :: silja_pasi = 230110

  !
  !        NAMES OF THE FORECAST-LENTGH-RULES
  !
  ! These are the rules, that user can give to the dataserver for
  ! defining what kind forecasts or analyses are used. The default
  ! value is shortest_forecasts. Long forecasts can be used for
  ! example test purposes.
  !
  INTEGER, PARAMETER, PUBLIC :: use_analyses = 230200
  INTEGER, PARAMETER, PUBLIC :: use_short_forecasts = 230201 
  INTEGER, PARAMETER, PUBLIC :: use_any_data = 230202
  INTEGER, PARAMETER, PUBLIC :: use_6h_forecasts = 230206 
  INTEGER, PARAMETER, PUBLIC :: use_12h_forecasts = 230212
  INTEGER, PARAMETER, PUBLIC :: use_24h_forecasts = 230224 
  INTEGER, PARAMETER, PUBLIC :: use_48h_forecasts = 230248 
  INTEGER, PARAMETER, PUBLIC :: default_forecasts = use_short_forecasts

  !
  !        INTERPOLATION METHODS and some constants
  !
  INTEGER, PARAMETER, PUBLIC :: nearest_point = 230300
  INTEGER, PARAMETER, PUBLIC :: linear = 230301
!  INTEGER, PARAMETER, PUBLIC :: second_order = 230302 !!!! Makes no sense! Was killed
  INTEGER, PARAMETER, PUBLIC :: cubic = 230303
!!  INTEGER, PARAMETER, PUBLIC :: log_linear = 230304 !! Never implemented
  integer, parameter, public :: average = 230305
  integer, parameter, public :: summation = 230306
  integer, parameter, public :: sigmoid = 230307
  integer, parameter, public :: double_sigmoid = 230308
  integer, parameter, public :: toMassCentreLinear = 230309  ! inside the grid cell to its mass centre

  !
  !        OUT OF GRID HANDLING 
  !
  integer, public, parameter :: notAllowed = 230320
  integer, public, parameter :: nearestPoint = 230321
  integer, public, parameter :: setZero = 230322
  integer, public, parameter :: setMissVal = 230323
  integer, public, parameter :: handleGlobalGrid = 230324

  !
  !        DIRECTIONS IN TIME AND HEIGHT
  !
  INTEGER, PARAMETER, PUBLIC :: backwards = 230400
  INTEGER, PARAMETER, PUBLIC :: forwards = 230401
  INTEGER, PARAMETER, PUBLIC :: back_and_forwards = 230402
  INTEGER, PARAMETER, PUBLIC :: single_time = 230403

  INTEGER, PARAMETER, PUBLIC :: upwards = 230404
  INTEGER, PARAMETER, PUBLIC :: downwards = 230405
  INTEGER, PARAMETER, PUBLIC :: up_and_downwards = 230406

  !  
  !       NAMES FOR DATA-ACCESS METHODS
  !
  INTEGER, PARAMETER, PUBLIC :: gribfile_access = 230500
  INTEGER, PARAMETER, PUBLIC :: fmi_database_access = 230502

  !
  !        NAMES FOR DIRECTIONS and borders
  !
  INTEGER, PARAMETER, PUBLIC :: north_flag = 230601
  INTEGER, PARAMETER, PUBLIC :: south_flag = 230602
  INTEGER, PARAMETER, PUBLIC :: east_flag = 230603
  INTEGER, PARAMETER, PUBLIC :: west_flag = 230604
  INTEGER, PARAMETER, PUBLIC :: northern_boundary = 1   ! these are indices in arrays 
  INTEGER, PARAMETER, PUBLIC :: southern_boundary = 2   ! Warning! This order is used in few places!
  INTEGER, PARAMETER, PUBLIC :: eastern_boundary  = 3
  INTEGER, PARAMETER, PUBLIC :: western_boundary  = 4
  INTEGER, PARAMETER, PUBLIC :: top_boundary      = 5
  INTEGER, PARAMETER, PUBLIC :: bottom_boundary   = 6
  character (len=6), dimension(6), parameter, public :: boundary_name= &
      &(/'North ', 'South ', 'East  ', 'West  ', 'Top   ', 'Bottom' /)

  !
  ! Mass budget includes incoming and outgoing masses: these are indices in arrays
  ! Funny order is because incoming are non-existent in many cases
  !
  integer, parameter, public :: outgoing = 1
  integer, parameter, public :: incoming = 2

  !
  ! Indices for boundary exchange MUST be 1 for our  and 2 for their 
  ! otherwise SILAM turns into a pumpkin !!! 
  integer, parameter, public :: left = 1
  integer, parameter, public :: right = 2
  character (len=5), dimension(2), parameter, public :: side_name=(/"left ","right"/)

  ! Indices for boundary exchange MUST be 1 for our  and 2 for their 
  ! otherwise SILAM turns into a pumpkin !!! 
  integer, parameter, public :: our = 1
  integer, parameter, public :: their = 2

  integer, parameter, public :: iPast = 1
  integer, parameter, public :: iFuture = 2
  integer, parameter, public :: iRealTime = 3

  ! Advection mechanisms:
  integer, parameter, public :: no_advection = 100001
!  integer, parameter, public :: adv_euler_Galperin_v4 = 100009
  integer, parameter, public :: adv_euler_Galperin_v5 = 100011
  integer, parameter, public :: adv_euler_Galperin_3d_bulk = 100010


  !  
  !       NAMES FOR GROUND SURFACE TYPES
  !
  INTEGER, PARAMETER, PUBLIC :: open_water_surface = 230701
  INTEGER, PARAMETER, PUBLIC :: open_ground_surface = 230702


  ! ***************************************************************

  !
  !
  !    SOME FLAGS FOR RETRIEVING MET. DATA 
  !    (used in shopping lists etc.)
  !
  ! ***************************************************************
  !
  !    Typical thickness of the layer considered to be subject for 
  !    dry deposition. Actual only for particle models, because in
  !    Eulerain models it is naturally the lowest model level
  !
  real, parameter, public :: Lagr_dry_dep_lyr = 50.  ! Meters

  !
  ! Types of the ABL height assessment methods
  !
  INTEGER, PARAMETER, PUBLIC :: constant_abl_height = 41
  INTEGER, PARAMETER, PUBLIC :: parcel_method = 42
  INTEGER, PARAMETER, PUBLIC :: richardson_method = 43
  INTEGER, PARAMETER, PUBLIC :: combination_method = 44
  INTEGER, PARAMETER, PUBLIC :: coriolis_method = 45
  INTEGER, PARAMETER, PUBLIC :: nwp_abl = 46

  
  ! Types of the ABL parameterizations
  !
  INTEGER, PARAMETER, PUBLIC :: abl_full_param = 51  ! Full-blown parameterization
  INTEGER, PARAMETER, PUBLIC :: abl_dry_param = 52   ! Dry ABL assumption

  ! Types of Kz parameterizations (essentially R_to_1m)
  !
  !INTEGER, PARAMETER, PUBLIC :: silam_old_kz = 61     ! New scheme with old Kz 
  !!!!No one can do  silam_old_kz anymore...

  INTEGER, PARAMETER, PUBLIC :: silam_kz_emulator = 62 ! Plain old scheme 
  INTEGER, PARAMETER, PUBLIC :: silam_abl_ec_ft_kz = 63  ! Old Kz within PBL, ECMWF stability above it 
  INTEGER, PARAMETER, PUBLIC :: ec_kz = 64 ! ECMWF stability above it 
  INTEGER, PARAMETER, PUBLIC :: zero_kz = 65 ! No vertical exchange
  INTEGER, PARAMETER, PUBLIC :: hunten_kz = 66 ! Hunten 1975 fixed Kz profile
  INTEGER, PARAMETER, PUBLIC :: simple_kz = 67 ! Weakly taken from chimere
  INTEGER, PARAMETER, PUBLIC :: simple_abl_ec_ft_kz = 68 !Simple Kz within PBL, ECMWF stability above it 


  !
  !  Types of the fields with regard to time resolution
  !
  ! For single-time fields
  integer, public, parameter :: static_value = 610        ! just a number -- forever valid
  integer, public, parameter :: static_climatology = 611  ! fixed in time map -- forever valid
  integer, public, parameter :: instant_map = 613         ! single map of instant values (almost never valid)
  integer, public, parameter :: period_valid_map = 614    ! updatable map
  integer, public, parameter :: monthly_climatology = 615 ! monthly updatable, same month is ok
  ! For multitime fields
  integer, public, parameter :: dynamic_map = 612         ! set of maps of instant values -- interpolatable

  !
  ! Short names of the months
  !
  CHARACTER (len=3), DIMENSION(12), PARAMETER, PUBLIC :: chMonthNames_3chr = &
    & (/'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'/)


  ! The following parameters can be used to control the compromise between
  ! accuracy and execution time of dose calculation:

  INTEGER, PARAMETER :: decay_method = 2  ! 1 = single-nuclide, 2 = chains

  ! Actions with the output files:
  !
  integer, public, parameter :: DoNothing = 220 ! Continue to write to current files
  integer, public, parameter :: SwitchBinary = 221 ! Keep ctl but change binary
  integer, public, parameter :: OpenFiles = 222 ! Totally new output files

  !
  ! Arrangement of the output files (GRIB and GrADS - the others are not affected)
  !
  integer, public, parameter :: all_in_one = 231 ! Just one GRIB and one GrADS file 
  integer, public, parameter :: hourly_new_file = 232 ! New file is started every hour
  integer, public, parameter :: daily_new_file = 233 ! New file is started every day
  integer, public, parameter :: monthly_new_file = 234 ! New file is started every month
  integer, public, parameter :: yearly_new_file = 235 ! New file is started every month
  !
  ! Type of averaging of the output variables
  !
  integer, public, parameter :: iAsIs = 501 ! No averaging, instant fields at output time
  integer, public, parameter :: iInstant = 502 ! No averaging, instant fields at output time
  integer, public, parameter :: iAverage = 503 ! Average betweent the output times
  integer, public, parameter :: iCumulative =504 ! Cumulative since the start of simulations
  integer, public, parameter :: iMeanLastHrs=505 ! Mean over last X hours up to output time
  integer, public, parameter :: iTotalWholePeriod=506 ! Total over whole computation peiord
  integer, public, parameter :: iTechnicalOriginalGrid = 507    ! as namelist, native grid
  integer, public, parameter :: iTechnicalDispersionGrid = 508  ! as namelist, dispersion grid


  !*****************************************************************
  !
  ! Common types and constants for defining the parameter assimilation 
  ! The interface for assimilated parameter
  !
  type assimParameter
    integer :: dim
    integer, dimension(:), allocatable :: conditions
    real, dimension(:), allocatable :: val
    real, dimension(:,:), allocatable :: cov
  end type assimParameter
  public assimParameter
  !
  ! Standard conditions for the above interface
  !
  integer, public, parameter :: non_negative = 801
  integer, public, parameter :: positive = 802
  integer, public, parameter :: range_zero_one = 803
  
  integer, parameter, public :: flag_3dvar = 3003
  integer, parameter, public :: flag_4dvar = 3004
  integer, parameter, public :: flag_h_matrix = 3005
  integer, parameter, public :: flag_4dvar_seq = 3006
  integer, parameter, public :: flag_enkf = 3007
  integer, parameter, public :: flag_enks = 3008

  
  ! ***************************************************************

  !
  ! Some module variables private to this module:

  ! When run in a batch queue, the program could easily generate millions of
  ! lines of error messages (in case that something goes wrong) without the
  ! user knowing anything about the situation. To prevent this at least
  ! partially, as much printing as possible is directed through the small
  ! subroutines of this module. The total number of lines printed by them will
  ! never exceed max_print_count.

  INTEGER, SAVE :: print_count = 0
  INTEGER, PARAMETER :: max_print_count = 50

  !*******************************************************************
  !
  ! Sometimes one needs more than 1- or 2-D arrays for temporary use.
  ! Here is the way to make it
  !
  type TrealPtr
    real, dimension(:), pointer :: ptr
  end type TrealPtr
  public TrealPtr



!*****************************Leftover of foremer silam_mpi module  

#ifdef SILAM_MPI   
  INTEGER, PARAMETER, PUBLIC :: smpi_proc_null = MPI_PROC_NULL
  LOGICAL, PARAMETER :: is_mpi_silam_version = .TRUE.

  ! These mpi constants are exported to next level
  PUBLIC MPI_LOGICAL, MPI_LOR, MPI_STATUS_IGNORE, MPI_SUCCESS, MPI_COMM_WORLD, MPI_PROC_NULL, &
       & MPI_UNDEFINED, MPI_OFFSET_KIND, MPI_INFO_NULL, MPI_ORDER_FORTRAN, MPI_REAL, MPI_INTEGER, &
       & MPI_DOUBLE_PRECISION, MPI_CHARACTER, MPI_STATUS_SIZE, MPI_ANY_TAG, MPI_ANY_SOURCE, &
       & MPI_INTEGER8, MPI_MODE_RDONLY, MPI_ERRORS_RETURN, MPI_MAX_ERROR_STRING, MPI_THREAD_SINGLE, &
       & MPI_MAX, MPI_SUM, MPI_MODE_WRONLY, MPI_MODE_CREATE
#else


  INTEGER, PUBLIC, PARAMETER :: MPI_UNDEFINED = -1, MPI_REAL = -1, &
       & MPI_SUCCESS = -1, MPI_LOGICAL = -1, MPI_LOR = -1, &
       & MPI_COMM_WORLD = -1,  MPI_INFO_NULL = -1
  INTEGER, PUBLIC, PARAMETER :: smpi_proc_null = -1
  LOGICAL, PRIVATE, PARAMETER :: is_mpi_silam_version = .FALSE.
#endif

#ifdef DOUBLE_PRECISION
  INTEGER, PARAMETER, PUBLIC :: smpi_real_type = MPI_DOUBLE_PRECISION
#else
  INTEGER, PARAMETER, PUBLIC :: smpi_real_type = MPI_REAL
#endif

  INTEGER, PUBLIC :: smpi_global_tasks = 1, smpi_global_rank = 0
  INTEGER, PUBLIC :: smpi_adv_tasks = 1, smpi_adv_rank =0
  INTEGER, PUBLIC :: smpi_adv_cart_rank = 0, smpi_io_rank = 0, smpi_ens_rank = 0
  ! MPI communicators:
  INTEGER, PUBLIC :: smpi_io_comm, smpi_ensmember_comm, smpi_enkf_comm 
  INTEGER, PUBLIC :: smpi_adv_comm, smpi_adv_cart_comm, smpi_adv_x_comm, smpi_adv_y_comm 

  !****************************END of MPI stuff

  !**********************Cinversion between concentrations mixing ratios and masses
  integer, public, parameter :: cloud_metric_geometry_flag = 190011       !! Use geometric cell size to convert mass to cnc/mixing_ratio
  integer, public, parameter :: cloud_metric_cellmass_flag = 190012       !! Use air mass in the cell -- acounts for compressibility
  integer, public, parameter :: cloud_metric_ones_flag = 190013        !! Use ones tracer -- accounts for compressibility and
                                                                       !!  advection artifacts

  ! Global error status variable is set here
  LOGICAL, SAVE, PUBLIC :: error

CONTAINS

  ! ****************************************************************

  ! 
  !
  !          TOOLS FOR SILJA_LOGICAL -TYPE
  !
  !
  ! ***************************************************************

  FUNCTION fu_set_true()

    ! Description:
    ! Sets a true value for a variable of silja_logical -type.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_logical) :: fu_set_true
    integer :: i

    fu_set_true = silja_true

  END FUNCTION fu_set_true



  ! ****************************************************************


  FUNCTION fu_set_false()

    ! Description:
    ! Sets a false value for a variable of silja_logical -type.
    ! 
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_logical) :: fu_set_false

    fu_set_false = silja_false

  END FUNCTION fu_set_false



  ! ****************************************************************


  LOGICAL FUNCTION fu_true (switch)

    ! Description:
    ! Return a true value, if the parameter switch is true.
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in): 
    TYPE(silja_logical) :: switch

    select case (switch%flag)
       case (true_flag)
         fu_true = .true.
       case  (false_flag, 0)
         fu_true = .false.
       case default
         fu_true = .false.
         call set_error("strange silja_logical flag", "fu_true")
     end select

  END FUNCTION fu_true




  ! ****************************************************************


  LOGICAL FUNCTION fu_false (switch)

    ! Description:
    ! Return a true value, if the parameter switch is false.
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in): 
    TYPE(silja_logical) :: switch

    select case (switch%flag)
       case (true_flag, 0)
         fu_false = .false.
       case  (false_flag)
         fu_false = .true.
       case default
         fu_false = .false.
         call set_error("strange silja_logical flag", "fu_false")
     end select

  END FUNCTION fu_false




  ! ****************************************************************


  LOGICAL FUNCTION fu_undefined (switch)

    ! Description:
    ! Return a true value, if the parameter switch is undefined, (that
    ! means it hasn't been given any value yet)
    !.
    ! Returns a false value, if switch is either true or false.
    ! 
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in): 
    TYPE(silja_logical) :: switch

    SELECT CASE (switch%flag)

      CASE (false_flag, true_flag)
      fu_undefined = .false.

    CASE default
      fu_undefined = .true.

    END SELECT

  END FUNCTION fu_undefined




  ! ****************************************************************


  LOGICAL FUNCTION fu_compare_logicals_eq (log1, log2)

    ! Description:
    ! Returns a true value, if same logicals.
    ! 
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in): 
    TYPE(silja_logical), INTENT(in) :: log1, log2

    IF ((log1%flag == true_flag).and.(log2%flag == true_flag)) THEN
      fu_compare_logicals_eq = .true.
    ELSE IF ((log1%flag == false_flag).and.(log2%flag == false_flag)) THEN
      fu_compare_logicals_eq = .true.
    ELSE IF ((fu_undefined(log1)).and.(fu_undefined(log2))) THEN
      fu_compare_logicals_eq = .true.
    ELSE
      fu_compare_logicals_eq = .false.
    END IF


  END FUNCTION fu_compare_logicals_eq


  !******************************************************************

  logical function fu_not_true_silja_logical(log1)
    implicit none

    type(silja_logical), intent(in) :: log1

    fu_not_true_silja_logical = log1%flag /= true_flag

  end function fu_not_true_silja_logical


  ! ****************************************************************

  ! 
  !
  !          TOOLS FOR HANDLING ERROR SITUATIONS
  !
  !
  ! ***************************************************************

  !**********************************************************************************************

  function fu_fails(condition, message, place) result(fails)
    implicit none
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message, place
    logical :: fails
    
    if (.not. condition) then
      call set_error(message, place)
      fails = .true.
    else
      fails = .false.
    end if
    
  end function fu_fails


  !**********************************************************************************************

  SUBROUTINE set_error (message, place)

    ! Sets an error situation.
    ! This subroutine is called whenever there occurs
    ! an error. 
    !
    ! Sets the values of global variables error, error_message and
    ! place_of_error.
    !
    ! Sends the error message to the user.
    ! Mika Salonoja, FMI 06-1995

    IMPLICIT NONE

    ! Imported parameters with intent(in): 
    CHARACTER (LEN=*), INTENT(in) :: message, place
    
    ! Local variables
    integer, save :: error_counter = 0
    logical :: fresh_error
    

    if (error) then 
#ifdef VS2012
      call sleepqq(1)
#else
      call sleep(1)
#endif
      fresh_error = .FALSE.
    else
      fresh_error = .TRUE.
    endif

    !$OMP CRITICAL(globals_error)
    error = .true.
    !$OMP END CRITICAL(globals_error)
    
    !$OMP CRITICAL(globals_io)
    IF(message /= '' .or. place /= '')THEN
      IF(smpi_global_rank == 0)THEN
      PRINT *, '#####################################################'
      if (fresh_error) PRINT *, 'Fresh'
      PRINT *, 'Error: ', TRIM(message)
      PRINT *,' Place of error: ', TRIM(place)
      PRINT *, '#####################################################'
      END IF
      write(run_log_funit,*)'#####################################################'
      if (fresh_error) write(run_log_funit,*)'Fresh'
      write(run_log_funit,*)'Error: ',trim(message)
      write(run_log_funit,*)'Place of error: ',trim(place)
      write(run_log_funit,*)'#####################################################'
    END IF
    !$OMP END CRITICAL (globals_io)
    call flush(run_log_funit)
#ifdef VS2012
    !$OMP CRITICAL(globals_error)
    call tracebackqq(user_exit_code=-1)
    !$OMP END CRITICAL(globals_error)
#elif defined(__GFORTRAN__) &&  (__GNUC__ > 4 || __GNUC__ == 4 &&  __GNUC_MINOR__ > 7)
    ! gfortran must be at least 4.8 for this
    !$OMP CRITICAL(globals_error)
    if(error_counter < 20 .and. fresh_error)then      ! due to Cray GNU bug, backtrace leaks ~40 MB of memory per call
      if (smpi_global_rank == 0 ) call backtrace()
      error_counter = error_counter + 1
    endif
    !$OMP END CRITICAL(globals_error)
#else
    call msg('(backtrace() not available in this compiler)')
#endif
  END SUBROUTINE set_error


  ! ****************************************************************


  SUBROUTINE set_error_nomsg (message, place)

    ! Sets an error situation.
    ! This subroutine is called whenever there occurs
    ! an error. 
    !
    ! Sets the values of global variables error, error_message and
    ! place_of_error.
    !
    ! Mika Salonoja, FMI 06-1995

    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(in) :: message, place

    !$OMP CRITICAL(globals_error)
    error = .true.
    !$OMP END CRITICAL(globals_error)

  END SUBROUTINE set_error_nomsg


  ! ****************************************************************

  SUBROUTINE unset_error(place)

    ! Unsets an error in case if the calling routine decides that the
    ! error is not fatal, and the running can continue.
    ! Sets new values for global variables error, error_message and
    ! place_of_error.
    !
    ! Mika Salonoja, FMI

    IMPLICIT NONE

    ! Imported parameters with intent(in): 
    CHARACTER (LEN=*), OPTIONAL, INTENT(in) :: place
    character (len=fnlen) :: msg

    !$OMP CRITICAL(globals_io)
    if (smpi_adv_tasks == 1) then
       msg = 'Unsetting error' 
    else
       msg = 'MPI: NOT Unsetting error' 
    endif

    IF (PRESENT(place)) THEN
      msg = trim(msg) // ' by: '//TRIM(place)
    else
      msg = trim(msg) // '.'
    endif

    IF (smpi_global_rank == 0) PRINT *,  TRIM(msg)
    write(run_log_funit,*) TRIM(msg)

    !$OMP END CRITICAL(globals_io)
   
    if (smpi_adv_tasks == 1) then
    !$OMP CRITICAL(globals_error)
    error = .false.
    !$OMP END CRITICAL (globals_error)
     endif

  END SUBROUTINE unset_error


  ! ****************************************************************


  SUBROUTINE msg_warning(message, place)

    ! Tells a warning to the user. 
    ! One integer or one real may be supplied with the
    ! message. 
    !
    ! Mika Salonoja, FMI 06-1995

    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(in) :: message 
    CHARACTER (LEN=*), INTENT(in), OPTIONAL :: place
    
    !$OMP CRITICAL (globals_io)
    IF (smpi_global_rank==0) PRINT*, '*** WARNING: ', TRIM(message)
    write(run_log_funit,*)'*** WARNING: ', TRIM(message)

    IF (PRESENT(place))then
       IF (smpi_global_rank == 0) PRINT*, 'Place of warning: ', TRIM(place)
       write(run_log_funit,*)'Place of warning: ', TRIM(place)
    endif
    !$OMP END CRITICAL (globals_io)

  END SUBROUTINE msg_warning



  ! ****************************************************************

  ! 
  !
  !          TOOLS FOR MESSAGES AND TEST MESSAGES 
  !
  !
  ! ***************************************************************

  subroutine ooops(message)
     IMPLICIT NONE
     CHARACTER (LEN=*), intent(in) :: message
     
#ifdef VS2012
    !$OMP CRITICAL(globals_error)
    call tracebackqq(user_exit_code=-1)
    !$OMP END CRITICAL(globals_error)
#elif defined(__GFORTRAN__) &&  (__GNUC__ > 4 || __GNUC__ == 4 &&  __GNUC_MINOR__ > 7)
    ! gfortran must be at least 4.8 for this
    !$OMP CRITICAL(globals_error)
    call backtrace()
    !$OMP END CRITICAL(globals_error)
#else
    call msg('(backtrace() not available in this compiler)')
#endif
     call msg(message)
  end subroutine ooops



  SUBROUTINE msg_txt(message)

    ! Tells a test message to the user

    IMPLICIT NONE

    CHARACTER (LEN=*), intent(in) :: message 

    !$OMP CRITICAL (globals_io)
    IF (smpi_global_rank==0) PRINT*,  trim(message)
    write(run_log_funit,fmt='(A)')trim(message)
    !$OMP END CRITICAL (globals_io)
  END SUBROUTINE msg_txt


  !************************************************************************************

  SUBROUTINE msg_int(message, int_value)

    ! Tells a test message to the user

    IMPLICIT NONE

    CHARACTER (LEN=*), intent(in) :: message 
    INTEGER, intent(in) :: int_value

    !$OMP CRITICAL (globals_io)
    IF (smpi_global_rank==0) PRINT*, trim(message), int_value
    if(run_log_funit > 0) write(run_log_funit,fmt='(A,1x,I12)')trim(message), int_value
    !$OMP END CRITICAL (globals_io)
  END SUBROUTINE msg_int


  !************************************************************************************

  SUBROUTINE msg_real4(message, real4_value)

    ! Tells a test message to the user

    IMPLICIT NONE

    CHARACTER (LEN=*), intent(in) :: message 
    REAL(4), intent(in) :: real4_value

    !$OMP CRITICAL (globals_io)
    if (smpi_global_rank==0) PRINT*, trim(message), real4_value
    if (real4_value /= real4_value) then
      write(run_log_funit,fmt='(A,1x,A15)')trim(message), "NaN"
    elseif(abs(real4_value) < 1.0e6 .and. abs(real4_value) > 1.e-6)then
      write(run_log_funit,fmt='(A,1x,F15.7)')trim(message), real4_value
    else
      write(run_log_funit,fmt='(A,1x,E15.7)')trim(message), real4_value
    endif
    !$OMP END CRITICAL (globals_io)
  END SUBROUTINE msg_real4


  !************************************************************************************

  SUBROUTINE msg_real8(message, real8_value)

    ! Tells a test message to the user

    IMPLICIT NONE

    CHARACTER (LEN=*), intent(in) :: message 
    REAL(r8k), intent(in) :: real8_value

    !$OMP CRITICAL (globals_io)
    if (smpi_global_rank==0) PRINT*, trim(message), real8_value
    if (real8_value /= real8_value) then
      write(run_log_funit,fmt='(A,1x,A25)')trim(message), "NaN"
    elseif(abs(real8_value) < 1.0e6 .and. abs(real8_value) > 1.e-6)then
      write(run_log_funit,fmt='(A,1x,F25.16)')trim(message), real8_value
    else
      write(run_log_funit,fmt='(A,1x,D25.16)')trim(message), real8_value
    endif
    !$OMP END CRITICAL (globals_io)
  END SUBROUTINE msg_real8


  !***********************************************************************************

  SUBROUTINE msg_int_real4(message, int_value, real4_value)

    ! Tells a test message to the user

    IMPLICIT NONE

    CHARACTER (LEN=*), intent(in) :: message 
    REAL(4), intent(in) :: real4_value
    INTEGER, intent(in) :: int_value

    !$OMP CRITICAL (globals_io)
    if (smpi_global_rank==0) PRINT*, trim(message), int_value, real4_value
    if(abs(real4_value) < 1.0e6 .and. abs(real4_value) > 1.e-6)then
      write(run_log_funit,fmt='(A,1x,I12,1x,F15.7)')trim(message), int_value, real4_value
    else
      write(run_log_funit,fmt='(A,1x,I12,1x,E15.7)')trim(message), int_value, real4_value
    end if
    !$OMP END CRITICAL (globals_io)
  END SUBROUTINE msg_int_real4


  !***********************************************************************************

  SUBROUTINE msg_int_real8(message, int_value, real8_value)

    ! Tells a test message to the user

    IMPLICIT NONE

    CHARACTER (LEN=*), intent(in) :: message 
    real(r8k), intent(in) :: real8_value
    INTEGER, intent(in) :: int_value

    !$OMP CRITICAL (globals_io)
    if (smpi_global_rank==0) PRINT*, trim(message), int_value, real8_value
    if(abs(real8_value) < 1.0e6 .and. abs(real8_value) > 1.e-6)then
      write(run_log_funit,fmt='(A,1x,I12,1x,F25.16)')trim(message), int_value, real8_value
    else
      write(run_log_funit,fmt='(A,1x,I12,1x,D25.16)')trim(message), int_value, real8_value
    end if
    !$OMP END CRITICAL (globals_io)
  END SUBROUTINE msg_int_real8


  !***********************************************************************************

  SUBROUTINE msg_real4_int(message, real_value, int_value)

    ! Tells a test message to the user

    IMPLICIT NONE

    CHARACTER (LEN=*), intent(in) :: message 
    REAL(4), intent(in) :: real_value
    INTEGER, intent(in) :: int_value

    !$OMP CRITICAL (globals_io)
    if (smpi_global_rank==0) PRINT*, trim(message), real_value, int_value
    if(abs(real_value) < 1.0e6 .and. abs(real_value) > 1.e-6)then
      write(run_log_funit,fmt='(A,1x,F15.7,1x,I12)')trim(message), real_value, int_value
    else
      write(run_log_funit,fmt='(A,1x,E15.7,1x,I12)')trim(message), real_value, int_value
    end if
    !$OMP END CRITICAL (globals_io)
  END SUBROUTINE msg_real4_int


  !***********************************************************************************

  SUBROUTINE msg_real8_int(message, real_value, int_value)

    ! Tells a test message to the user

    IMPLICIT NONE

    CHARACTER (LEN=*), intent(in) :: message 
    REAL(r8k), intent(in) :: real_value
    INTEGER, intent(in) :: int_value

    !$OMP CRITICAL (globals_io)
    if (smpi_global_rank==0) PRINT*, trim(message), real_value, int_value
    if(abs(real_value) < 1.0e6 .and. abs(real_value) > 1.e-6)then
      write(run_log_funit,fmt='(A,1x,F15.7,1x,I12)')trim(message), real_value, int_value
    else
      write(run_log_funit,fmt='(A,1x,D15.7,1x,I12)')trim(message), real_value, int_value
    end if
    !$OMP END CRITICAL (globals_io)
  END SUBROUTINE msg_real8_int


  !***********************************************************************************

  SUBROUTINE msg_int_int(message, int1, int2)

    ! Tells a test message to the user

    IMPLICIT NONE

    CHARACTER (LEN=*), intent(in) :: message 
    INTEGER, intent(in):: int1, int2

    !$OMP CRITICAL (globals_io)
    if (smpi_global_rank==0) PRINT*, trim(message), int1, int2
    write(run_log_funit,fmt='(A,1x,I12,1x,I12)')trim(message), int1, int2
    !$OMP END CRITICAL (globals_io)
  END SUBROUTINE msg_int_int


  !*********************************************************************

  SUBROUTINE msg_real4_real4(message, real1, real2)

    ! Tells a test message to the user

    IMPLICIT NONE

    CHARACTER (LEN=*), intent(in) :: message 
    REAL(r4k), intent(in) :: real1, real2

    !$OMP CRITICAL (globals_io)
    if (smpi_global_rank==0) PRINT*, trim(message), real1, real2
    if((abs(real1) < 1.0e6 .and. abs(real1) > 1.e-6) .and. &
     & (abs(real2) < 1.0e6 .and. abs(real2) > 1.e-6))then
      write(run_log_funit,fmt='(A,1x,F15.7,1x,F15.7)')trim(message), real1, real2
    else
      write(run_log_funit,fmt='(A,1x,E15.7,1x,E15.7)')trim(message), real1, real2
    endif
    !$OMP END CRITICAL (globals_io)
  END SUBROUTINE msg_real4_real4


  !*********************************************************************

  SUBROUTINE msg_real8_real8(message, real1, real2)

    ! Tells a test message to the user

    IMPLICIT NONE

    CHARACTER (LEN=*), intent(in) :: message 
    REAL(r8k), intent(in) :: real1, real2

    !$OMP CRITICAL (globals_io)
    if (smpi_global_rank==0) PRINT*, trim(message), real1, real2
    if((abs(real1) < 1.0e6 .and. abs(real1) > 1.e-6) .and. &
     & (abs(real2) < 1.0e6 .and. abs(real2) > 1.e-6))then
      write(run_log_funit,fmt='(A,1x,F25.16,1x,F25.16)')trim(message), real1, real2
    else
      write(run_log_funit,fmt='(A,1x,D25.16,1x,D25.16)')trim(message), real1, real2
    endif
    !$OMP END CRITICAL (globals_io)
  END SUBROUTINE msg_real8_real8


  !*****************************************************************************
  
  subroutine msg_logical_array(message, log_array, counter)

    ! Tells a test message to the user

    IMPLICIT NONE

    CHARACTER (LEN=*), intent(in) :: message 
    logical, dimension(:), intent(in) :: log_array
    integer, intent(in), optional :: counter
    
    ! Local parameters
    integer :: iTmp, jTmp, iSizeTmp, iSizeTmp1
    character(len=99) :: strtmp

    !$OMP CRITICAL (globals_io)
    if(present(counter))then
      iSizeTmp = min(size(log_array),counter)
    else
      iSizeTmp = size(log_array)
    endif

    if (smpi_global_rank==0)  write(6,fmt='(A)',advance="no") trim(message)
    write(run_log_funit,fmt='(A)',advance="no") trim(message)
    do iTmp = 1, iSizeTmp, 99
      iSizeTmp1=min(99, iSizeTmp - iTmp + 1) ! Size of a string
      !Replace "F" with "." top make it mode visible
      write(unit=strtmp,fmt='(99L1)')(log_array(iTmp:min(itmp+98,iSizeTmp)))
      do jTmp = 1, iSizeTmp1 
          if (strtmp(jTmp:jTmp)=='F') strtmp(jTmp:jTmp)='o'
      enddo

      if (smpi_global_rank==0) write(6,fmt='(A)',advance='no') strtmp(1:iSizeTmp1)
      write(run_log_funit,fmt='(A)',advance='no') strtmp(1:iSizeTmp1)
!      if (smpi_global_rank==0)  &
!                 &write(6,fmt='(100L1)',advance='no')(log_array(iTmp:min(itmp+99,iSizeTmp)))
!      write(run_log_funit,fmt='(100L1)',advance='no')(log_array(iTmp:min(itmp+99,iSizeTmp)))
    enddo
    ! put EOL
    write(run_log_funit,*) ""
    if (smpi_global_rank==0) write(6,*) ""

    !$OMP END CRITICAL (globals_io)
  end subroutine msg_logical_array

  !*****************************************************************************
  subroutine msg_int_array(message, int_array, counter)

    ! Tells a test message to the user

    IMPLICIT NONE

    CHARACTER (LEN=*), intent(in) :: message 
    integer, dimension(:), intent(in) :: int_array
    integer, intent(in), optional :: counter
    
    ! Local parameters
    integer :: iTmp, iSizeTmp

    !$OMP CRITICAL (globals_io)
    if(present(counter))then
      iSizeTmp = min(size(int_array),counter)
    else
      iSizeTmp = size(int_array)
    endif

    if (smpi_global_rank==0) write(6,fmt='(A)',advance="no") trim(message)
    write(run_log_funit,fmt='(A)',advance="no") trim(message)
      do iTmp = 1, iSizeTmp, 100
      if (smpi_global_rank==0) write(6,fmt='(100(1x,I9))')(int_array(iTmp:min(itmp+100,iSizeTmp)))
      write(run_log_funit,fmt='(100(1x,I9))')(int_array(iTmp:min(itmp+100,iSizeTmp)))
    enddo
    !$OMP END CRITICAL (globals_io)
  end subroutine msg_int_array

  !*****************************************************************************
  
  subroutine msg_int8_array(message, int_array, counter)

    ! Tells a test message to the user

    IMPLICIT NONE

    CHARACTER (LEN=*), intent(in) :: message 
    integer (kind=8), dimension(:), intent(in) :: int_array
    integer, intent(in), optional :: counter
    
    ! Local parameters
    integer :: iTmp, iSizeTmp

    !$OMP CRITICAL (globals_io)
    if(present(counter))then
      iSizeTmp = min(size(int_array),counter)
    else
      iSizeTmp = size(int_array)
    endif

    if (smpi_global_rank==0) write(6,fmt='(A)',advance="no") trim(message)
    write(run_log_funit,fmt='(A)',advance="no") trim(message)
      do iTmp = 1, iSizeTmp, 100
      if (smpi_global_rank==0) write(6,fmt='(100(1x,I10))')(int_array(iTmp:min(itmp+100,iSizeTmp)))
      write(run_log_funit,fmt='(100(1x,I10))')(int_array(iTmp:min(itmp+100,iSizeTmp)))
    enddo
    !$OMP END CRITICAL (globals_io)
  end subroutine msg_int8_array

  !*****************************************************************************
  
  subroutine msg_real_array(message, real_array, counter)

    ! Tells a test message to the user

    IMPLICIT NONE

    CHARACTER (LEN=*), intent(in) :: message 
    real(r4k), dimension(:), intent(in) :: real_array
    integer, intent(in), optional :: counter
    
    ! Local parameters
    integer :: iTmp, iSizeTmp

    !$OMP CRITICAL (globals_io)
    if(present(counter))then
      iSizeTmp = min(size(real_array),counter)
    else
      iSizeTmp = size(real_array)
    endif
!    if (smpi_global_rank==0) PRINT '(A)', trim(message)
!    if (smpi_global_rank==0) PRINT '(E9.3,1x)', real_array(1:iSizeTmp)

    if (smpi_global_rank==0) then
      write(6,fmt='(A)',advance="no") trim(message)
      do iTmp = 1, iSizeTmp, 100
        write(6,fmt='(100(2x,E11.5))')(real_array(iTmp:min(itmp+100,iSizeTmp)))
      enddo
    end if

    write(run_log_funit,fmt='(A)',advance="no") trim(message)
    do iTmp = 1, iSizeTmp, 100
      write(run_log_funit,fmt='(100(2x,E11.5))')(real_array(iTmp:min(itmp+100,iSizeTmp)))
    enddo
    !$OMP END CRITICAL (globals_io)
  end subroutine msg_real_array


  ! *****************************************************************
  subroutine msg_real8_array(message, real_array, counter)

    ! Tells a test message to the user

    IMPLICIT NONE

    CHARACTER (LEN=*), intent(in) :: message 
    real(kind=8), dimension(:), intent(in) :: real_array
    integer, intent(in), optional :: counter
    
    ! Local parameters
    integer :: iTmp, iSizeTmp

    !$OMP CRITICAL (globals_io)
    if(present(counter))then
      iSizeTmp = min(size(real_array),counter)
    else
      iSizeTmp = size(real_array)
    endif

    if (smpi_global_rank==0) then
      write(6,fmt='(A)',advance="no") trim(message)
      do iTmp = 1, iSizeTmp, 100
        write(6,fmt='(100(2x,E11.5))')(real_array(iTmp:min(itmp+100,iSizeTmp)))
      enddo
    end if

    write(run_log_funit,fmt='(A)',advance="no") trim(message)
    do iTmp = 1, iSizeTmp, 100
      write(run_log_funit,fmt='(100(2x,E11.5))')(real_array(iTmp:min(itmp+100,iSizeTmp)))
    enddo
    !$OMP END CRITICAL (globals_io)
  end subroutine msg_real8_array


  ! *****************************************************************
  !
  ! Same as above but conditional to test_messages switch
  !
  ! *****************************************************************

  SUBROUTINE msg_test_txt(message)

    ! Tells a test message to the user

    IMPLICIT NONE

    CHARACTER (LEN=*), intent(in) :: message 

    if(.not.test_messages)return

    !$OMP CRITICAL (globals_io)
    if (smpi_global_rank==0) PRINT*,  trim(message)
    write(run_log_funit,fmt='(A)')trim(message)
    !$OMP END CRITICAL (globals_io)
  END SUBROUTINE msg_test_txt


  !************************************************************************************

  SUBROUTINE msg_test_int(message, int_value)

    ! Tells a test message to the user

    IMPLICIT NONE

    CHARACTER (LEN=*), intent(in) :: message 
    INTEGER, intent(in) :: int_value

    if(.not.test_messages)return

    !$OMP CRITICAL (globals_io)
    if (smpi_global_rank==0) PRINT*, trim(message), int_value
    write(run_log_funit,fmt='(A,1x,I12)')trim(message), int_value
    !$OMP END CRITICAL (globals_io)
  END SUBROUTINE msg_test_int


  !************************************************************************************

  SUBROUTINE msg_test_real(message, real_value)

    ! Tells a test message to the user

    IMPLICIT NONE

    CHARACTER (LEN=*), intent(in) :: message 
    REAL, intent(in) :: real_value

    if(.not.test_messages)return

    !$OMP CRITICAL (globals_io)
    if (smpi_global_rank==0) PRINT*, trim(message), real_value
    if(abs(real_value) < 1.0e6 .and. abs(real_value) > 1.e-6)then
      write(run_log_funit,fmt='(A,1x,F15.7)')trim(message), real_value
    else
      write(run_log_funit,fmt='(A,1x,E15.7)')trim(message), real_value
    endif
    !$OMP END CRITICAL (globals_io)
  END SUBROUTINE msg_test_real


  !***********************************************************************************

  SUBROUTINE msg_test_int_real(message, int_value, real_val)

    ! Tells a test message to the user

    IMPLICIT NONE

    CHARACTER (LEN=*), intent(in) :: message 
    REAL, intent(in) :: real_val
    INTEGER, intent(in) :: int_value

    if(.not.test_messages)return

    !$OMP CRITICAL (globals_io)
    if (smpi_global_rank==0) PRINT*, trim(message), int_value, real_val
    if(abs(real_val) < 1.0e6 .and. abs(real_val) > 1.e-6)then
      write(run_log_funit,fmt='(A,1x,I12,1x,F15.7)')trim(message), int_value, real_val
    else
      write(run_log_funit,fmt='(A,1x,I12,1x,E15.7)')trim(message), int_value, real_val
    end if
    !$OMP END CRITICAL (globals_io)
  END SUBROUTINE msg_test_int_real


  !***********************************************************************************

  SUBROUTINE msg_test_real_int(message, real_value, int_value)

    ! Tells a test message to the user

    IMPLICIT NONE

    CHARACTER (LEN=*), intent(in) :: message 
    REAL, intent(in) :: real_value
    INTEGER, intent(in) :: int_value

    if(.not.test_messages)return

    !$OMP CRITICAL (globals_io)
    if (smpi_global_rank==0) PRINT*, trim(message), real_value, int_value
    if(abs(real_value) < 1.0e6 .and. abs(real_value) > 1.e-6)then
      write(run_log_funit,fmt='(A,1x,F15.7,1x,I12)')trim(message), real_value, int_value
    else
      write(run_log_funit,fmt='(A,1x,E15.7,1x,I12)')trim(message), real_value, int_value
    end if
    !$OMP END CRITICAL (globals_io)
  END SUBROUTINE msg_test_real_int


  !***********************************************************************************

  SUBROUTINE msg_test_int_int(message, int1, int2)

    ! Tells a test message to the user

    IMPLICIT NONE

    CHARACTER (LEN=*), intent(in) :: message 
    INTEGER, intent(in):: int1, int2

    if(.not.test_messages)return

    !$OMP CRITICAL (globals_io)
    if (smpi_global_rank==0) PRINT*, trim(message), int1, int2
    write(run_log_funit,fmt='(A,1x,I12,1x,I12)')trim(message), int1, int2
    !$OMP END CRITICAL (globals_io)
  END SUBROUTINE msg_test_int_int


  !*********************************************************************

  SUBROUTINE msg_test_real_real(message, real1, real2)

    ! Tells a test message to the user

    IMPLICIT NONE

    CHARACTER (LEN=*), intent(in) :: message 
    REAL, intent(in) :: real1, real2

    if(.not.test_messages)return

    !$OMP CRITICAL (globals_io)
    if (smpi_global_rank==0) PRINT*, trim(message), real1, real2
    if((abs(real1) < 1.0e6 .and. abs(real1) > 1.e-6) .and. &
     & (abs(real2) < 1.0e6 .and. abs(real2) > 1.e-6))then
      write(run_log_funit,fmt='(A,1x,F15.7,1x,F15.7)')trim(message), real1, real2
    else
      write(run_log_funit,fmt='(A,1x,E15.7,1x,E15.7)')trim(message), real1, real2
    endif
    !$OMP END CRITICAL (globals_io)
  END SUBROUTINE msg_test_real_real

  
  ! ***************************************************************

  function fu_connect_strings(str1, str2, str3, str4, str5, str6) result(connected)

    ! Description:
    ! Connects character strings and returns them in one. The strings
    ! 3, 4 and 5 are optional.
    ! Trick. Digital debugger does not like the following sequence:
    ! write(...) TRIM(ADJUSTL(str)),TRIM(ADJUSTL(str)),TRIM(ADJUSTL(str))
    ! It leads to an immediate crush of debugger. So, a compromise
    ! solution is applied below.
    !
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    CHARACTER (LEN=worksize_string) :: connected

    ! Imported parameters with intent(in):
    CHARACTER (LEN=*), INTENT(in) :: str1, str2
    CHARACTER (LEN=*), INTENT(in), OPTIONAL :: str3, str4, str5, str6
    integer :: iTmp 

    connected = str1(fu_index_non_delim(str1):len_trim(str1)) // &
              & str2(fu_index_non_delim(str2):len_trim(str2))

    if(present(str3))then
      iTmp=len_trim(connected)
      connected = connected(1:iTmp) // str3(fu_index_non_delim(str3):len_trim(str3))
      if(present(str4))then
        iTmp=len_trim(connected)
        connected = connected(1:iTmp) // str4(fu_index_non_delim(str4):len_trim(str4))
        if(present(str5))then
          iTmp=len_trim(connected)
          connected = connected(1:iTmp) // str5(fu_index_non_delim(str5):len_trim(str5))
          if(present(str6)) then
            iTmp=len_trim(connected)
            connected = connected(1:iTmp) // str6(fu_index_non_delim(str6):len_trim(str6))
          endif
        endif
      endif
    endif
    
  end function fu_connect_strings

  
  !*****************************************************************************
  
  integer function fu_index_non_delim(strIn)
    implicit none
    ! imported parameter
    character(len=*), intent(in) :: strIn
    ! Local declarations
                                                ! space-tab-linefeed-carriage_return
    character(len=4), parameter :: chDelimiters = ' ' // char(11) // char(10) // char(13) 

    do fu_index_non_delim = 1, len(strIn)
      if(index(chDelimiters, strIn(fu_index_non_delim:fu_index_non_delim)) == 0) return
    end do
  end function fu_index_non_delim



   !*******************************************************
  
  subroutine report_vectorstat(str,arr, n)
      implicit none
      character(len = *) :: str
      real, dimension(:), intent(in) :: arr
      integer, intent(in) :: n
      character(len = fnlen) :: chTmp
      character(len = *), parameter :: sub_name = 'report_vectorstat'

      if (n > 0) then
        write (chTmp, *) trim(str), ': len=', n, ' min=', minval(arr(1:n)),  ' max=', maxval(arr(1:n)), &
            &' mean=', sum(arr(1:n))/n, ' norm=',  sum(arr(1:n)**2)/n
       elseif (n == 0) then
         write (chTmp, *) trim(str), "len= 0"
       else
         call set_error("Strange vector length n="//trim(fu_str(n)), sub_name)
         return
       endif
       call msg(chTmp)

  end subroutine report_vectorstat

  !******************************************************************

  function fu_connect_two_strings(str1,str2)result(connected)
    !
    ! Connects two character strings and returns them in one. 
    ! In fact, just calls more universal fu_connect_strings
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    CHARACTER (LEN=fnlen) :: connected

    ! Imported parameters with intent(in):
    CHARACTER (LEN=*), INTENT(in) :: str1, str2

    connected = fu_connect_strings(str1,str2)

  END FUNCTION fu_connect_two_strings

  
  !*****************************************************************************


  FUNCTION fu_int2str(int, iLengthRequested, ifForce) result(str)
    !
    ! Prints integer to a string and adjusts it to the left
    !
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    CHARACTER (LEN=20) :: str
    !
    ! Imported parameters
    integer, intent(in) :: int
    integer, intent(in), optional :: iLengthRequested  ! desired length of the field
    logical, intent(in), optional :: ifForce  ! whether the iLength can be ignored

    ! Local variables
    integer :: iTmp, iLength

    write(str,fmt='(I12)')int
    str = adjustl(str)
    iLength = len_trim(str)
    !
    ! If the length is requested, we can add leading zeroes and set error if the number is too long
    !
    if(present(iLengthRequested))then
      if(iLength > iLengthRequested)then  ! too long value
        if(present(ifForce))then
          if(ifForce)then           ! the length is enforced: error
            call msg('The integer is longer than the given space:',int, iLengthRequested)
            call set_error('The integer is longer than the given space','fu_int2str')
            str = ''
          endif
        endif
      elseif(iLength < iLengthRequested)then  ! too short string. Shift the string and add leading zeroes
        do iTmp = 0, iLength-1
          str(iLengthRequested-iTmp:iLengthRequested-iTmp) = str(iLength-iTmp:iLength-iTmp)
        end do
        do iTmp = 1, iLengthRequested - iLength
          str(iTmp:iTmp) = '0'
        end do
      endif
    endif

  END FUNCTION fu_int2str


  ! ****************************************************************
  FUNCTION fu_real82str(r) result(str)
    !
    ! Prints integer to a string and adjusts it to the left
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    CHARACTER (LEN=20) :: str
    !
    ! Imported parameter
    real (kind=8), intent(in) :: r
    
    if (r /= r) then
      write(str,fmt=*)r
    elseif(abs(r)>1.e5 .or. abs(r)<1.e-3)then
      write(str,fmt='(E15.7)')r
    else
      write(str,fmt='(F15.7)')r
    endif
    str = adjustl(str)

  END FUNCTION fu_real82str


  FUNCTION fu_real2str(r) result(str)
    !
    ! Prints integer to a string and adjusts it to the left
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    CHARACTER (LEN=20) :: str
    !
    ! Imported parameter
    real, intent(in) :: r

    if (r /= r) then
      write(str,fmt=*)r
    elseif(abs(r)>1.e5 .or. abs(r)<1.e-3)then
      write(str,fmt='(E15.7)')r
    else
      write(str,fmt='(F15.7)')r
    endif
    str = adjustl(str)

  END FUNCTION fu_real2str


  ! ****************************************************************

  ! 
  !
  !  TOOLS FOR RETURNING THE NAMES OF SOME PARAMETERS FOR IO-USE
  !
  !
  ! ***************************************************************


  FUNCTION fu_name_of_integer(intVal) result(string)
    ! 
    ! We have plenty of various switches in the model. They are deemed to be unique. So, this function
    ! gives names to these switches, also serving as a single point of references for the empty slots, 
    ! for instance.
    ! ATTENTION. 
    ! Not all switches are included below. It is up to everyone to keep that list populated
    !
    !  OBS: constants are hard-coded here. They might or might not match the parameter definitions...
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    CHARACTER(LEN=clen) :: string

    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: intVal
    !if(fu_known_qnaitity(intVal))then
    !  string = fu_quantity_short_string(intVal)
    !else
      select case(intVal)
        case(5001)               !   transformation_passive = 5001
          string = 'trnsf_passive'
        case(5002)               !   transformation_radioactive
          string = 'trnsf_radioact'
        case(5003)               !   transformation_polen
          string = 'trnsf_pollen'
        case(5004)               !   transformation_PM
          string = 'trnsf_pm'
        case(5005)               !   transformation_S_DMAT
          string = 'trnsf_sulphur_dmat'
        case(5006)               !   transformation_acid_basic
          string = 'trnsf_acid_basic'

        case(5009)               !   transformation_CB4
          string = 'trnsf_cb4'
        case(5010)               !   transformation_CB42_strato
          string = 'trnsf_cb42_strato'
        case(5011)               !   transformation_CB42_SOA
          string = 'trnsf_cb4_SOA'
          
        case(5020)
          string = 'aer_dyn_basic'
        case(5021)                 !aerosol_dynamics_simple
          string = 'aer_dyn_simple'
        case(5022)                 !aerosol_dynamics_mid_atm)
          string = 'aer_dyn_mid_atm'
        case(5023)               !aerosol_dynamics_mid_atm)
          string = 'aer_dyn_VBS'          

  !
  ! Types of interpolation and other shapes of smoothers
  !
        case(nearest_point)
          string = 'nearest_point'
        case(linear)
          string = 'linear'
        case(cubic)
          string = 'cubic'
!        case(log_linear)
!          string = 'log_linear'
        case(average)
          string = 'average'
        case(summation)
          string = 'summation'
        case(sigmoid)
          string = 'sigmoid'
        case(double_sigmoid)
          string = 'double_sigmoid'
  
          ! out of grid interpolation
         case(10010) 
                string = 'notAllowed'
         case(10011) 
                string = 'notAllowed_lon_global'
         case(10012) 
                string = 'nearestPoint'
         case(10013) 
                string = 'setZero'
         case(10014) 
                string = 'setMissVal'
         case(10015) 
                string = 'handleGlobalGrid'

        case(7000) 
          string = 'gas_phase'       ! chemical_setup
        case(7001)
          string = 'fixed_diameter'  ! chemical_setup
        case(7002)
          string = 'gamma_function'  ! chemical_setup
        case(7003) 
          string = 'moving_diameter' ! chemical_setup
        case(7004)
          string = 'lognormal'       ! chemical_setup
        case(7005)
          string = 'sea_salt_mode'  ! chemical_setup
        !
        !case(100004)
        !  string = 'adv_euler_Galperin_horiz_v2'
        !case(100005)
        !  string = 'adv_euler_Galperin_horiz_v3'
        !case(100006)
        !  string = 'adv_euler_Galperin_horiz_v4'
        !case(100007)
        !  string = 'adv_euler_Galperin_vert_v2'
        !case(100008)
        !  string = 'adv_euler_Galperin_vert_v3'
        !case(100009)
        !  string = 'adv_euler_Galperin_vert_v4'
        !case(100010)
        !  string = 'adv_euler_Galperin_3d_bulk'
        !case(100011)
        !  string = 'adv_euler_Galperin_v5'

  !
  ! Types of the ABL height assessment methods
  !
          case(constant_abl_height) 
              string = 'constant_abl_height'
          case(parcel_method) 
              string = 'parcel_method'
          case(richardson_method) 
              string = 'richardson_method'
          case(combination_method) 
              string = ' combination_method'
          case(coriolis_method) 
              string = 'coriolis_method'
          case(nwp_abl) 
              string = 'nwp_abl'

  
  ! Types of the ABL parameterizations
  !
          case(abl_full_param) 
              string = ' Full-blown parameterization'
          case(abl_dry_param) 
              string = ' Dry ABL assumption'

  ! Types of Kz parameterizations (essentially R_to_1m)
  !
!          case(silam_old_kz) 
!              string = ' New scheme with old Kz'
          case(silam_kz_emulator) 
              string = ' Plain old scheme '
          case(silam_abl_ec_ft_kz) 
              string = ' Old Kz within PBL, ECMWF stability above it '
          case(simple_abl_ec_ft_kz) 
              string = ' Simple Kz within PBL, ECMWF stability above it '
          case(ec_kz) 
              string = ' ECMWF stability above'
          case(zero_kz) 
              string = ' No vertical exchange'

  !
  !  Types of the fields with regard to time resolution
  !
          case(static_value) 
              string = ' just a number'
          case(static_climatology) 
              string = ' fixed in time map'
          case(monthly_climatology) 
              string = ' monthly updated map'
          case(dynamic_map) 
              string = ' map updated each <whatever> time step'

  !control the compromise between accuracy and execution time of dose calculation:

          case(decay_method) 
              string = ' 1 = single-nuclide, 2 = chains'

  ! Actions with the output files:
  !
          case(DoNothing) 
              string = ' Continue to write to current files'
          case(SwitchBinary) 
              string = ' Keep ctl but change binary'
          case(OpenFiles) 
              string = ' Totally new output files'

  !
  ! Arrangement of the output files (GRIB and GrADS - the others are not affected)
  !
          case(all_in_one) 
              string = ' Just one GRIB and one GrADS file '
          case(hourly_new_file) 
              string = ' New file is started every hour'
          case(daily_new_file) 
              string = ' New file is started every day'
          case(monthly_new_file) 
              string = ' New file is started every month'
          case(yearly_new_file) 
              string = ' New file is started every month'
  !
  ! Type of averaging of the output variables
  !
          case(iAsIs) 
              string = ' No averaging, instant fields at output time'
          case(iInstant) 
              string = ' No averaging, instant fields at output time'
          case(iAverage) 
              string = ' Average betweent the output times'
          case(iCumulative) 
              string = ' Cumulative since the start of simulations'
          case(iMeanLastHrs) 
              string = ' Mean over last X hours up to output time'
          case(iTotalWholePeriod) 
              string = ' Total over whole computation peiord'
          case(iTechnicalOriginalGrid) 
              string = ' as namelist, native grid'
          case(iTechnicalDispersionGrid) 
              string = ' as namelist, dispersion grid'

          case(no_advection)       ! = 100001
              string = ' no advection'
!          case(adv_euler_Galperin_v4)       ! = 100009
!              string = ' as namelist, dispersion grid'
          case(adv_euler_Galperin_v5)       ! = 100011
              string = ' advection Galperin v.5'
          case(adv_euler_Galperin_3d_bulk)       ! = 100010
              string = ' advection Galperin v.3 bulk'

        case default
          string = fu_str(intVal)
        end select
    !endif

  END FUNCTION fu_name_of_integer

  ! **********************************************************************************

  FUNCTION smpi_is_mpi_version() RESULT(val)
    !
    ! This routine can be used to check if the binary was compiled
    ! with MPI support
    !
    LOGICAL :: val
    val = is_mpi_silam_version
  END FUNCTION smpi_is_mpi_version

! ***********************************************************************************
END MODULE globals
